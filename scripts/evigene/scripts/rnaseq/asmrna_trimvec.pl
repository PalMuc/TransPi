#!/usr/bin/env perl
# evigene/scripts/rnaseq/asmrna_trimvec.pl

=item about asmrna_trimvec 

  merge of veccutmrna2.pl and parts from evigene/scripts/evgmrna2tsa.pl
  - revised and complex script for mRNA with vector screen, removal and NNN gap trimming
  - detects vectors in mRNA assemblies (using vecscreen/blastn -d UniVec_core data)
  - trims end gaps, per NCBI TSA requirements
  - uses protein/CDS info with ref-protein to prevent spurious/poor screen/trim that damages good mRNA
  - pulled simpler version in evgmrna2tsa, tested in veccutmrna2 to work right

  - recomputes protein,cds of trimmed mRNA and merges for evgmrna2tsa uses;
  - logs ambigious/problem cases for expert inspection while trying to minimize this,
    tries to keep valid protein while trimming vectors/gaps.

=item update for publicset uses

  replace some of evgmrna2tsa2.pl work, -nodeferupdate does merge of pubset seq + trimset seq
  $evigene/scripts/rnaseq/asmrna_trimvec.pl -nodeferupdate  -mrna publicset/kfish2evgxx11.mrna -log -debug

=item FIXME 201807 tsancbigapck.pl

 damn ever changing ncbi gap policies;
 need to add changable gap policy scanner, for mrna/cds to fsa/tbl
 ncbi tsa submit now has different endgap policy 
   not same as their current tbl2asn checker
   and it takes hours and hours to get thru the ncbi web tsa checker to find this out..
   and days of rechecking to figure out computation of that text-defined gap policy  

 /bio/bio-grid/verts/sra2genes/pig18evgsub
 tsancbigapck.pl

=item update 2017

  minor updates for evgpipe_sra2genes.pl
  See also asmrna_trimvecsiam.pl : should merge these, vecsiam has newer bug fixes
  
=item 2014/15 bugs
  
  * bug: leaves terminal nnnnnnnnnnn on several trs, both ends, no reason to not trim (in utrs)
    cant see why ; hasBadGaps() or endTrimNNN() are failing to catch 
      
=cut

use constant VERSION => '2018.07.04'; # updates for evgpipe_sra2genes.pl, -degaponly, other small changes
  # '2017.12.10'; # updates for evgpipe_sra2genes.pl
  # See also asmrna_trimvecsiam.pl : should merge these, vecsiam has newer bug fixes
  # '2015.01.05'; # more bug fixes, improvements;
  # '2014.12.26'; # upd for vecscreen v2/ncbic++ variant; '2013.06.01'; # 05.28

use FindBin;
use lib ($FindBin::Bin,"$FindBin::Bin/.."); # assume evigene/scripts/rnaseq/thisscript.pl

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname);
use cdna_proteins;
use cdna_evigenesub;

# cdna_evigenesub globals:
our $EVIGENES="$FindBin::Bin/..";  
our $EGAPP='egtrimvec';  
our $dryrun=0; ## $DRYRUN ?
our $DEBUG= $ENV{debug}|| 0;
#drop# my $UVGAP=  $ENV{uvgap}||0;   ## NO: leave nnn gaps for cut vec, preserve codon%3 size?

## Evigene tr2aacds subdirs: see tidyup
## add for trimvec, mrna2tsa:  trimset?  publicset? tsasubmit/submitset ?
## change trimvec,mrna2tsa output subdir: publicset? adding pubids, main2alt, ...
## separate tsasubmit subdir ..
our @evgdirs = qw(okayset dropset inputset tmpfiles erasefiles publicset);
our (@okayset,@dropset,@inputset,@tmpfiles,@erasefiles,@publicset); # tidyup file sets

our $vecoutdir='publicset'; # trimset now ? or leave in same folder as input mrna ?

use constant{ HasGapNone => 0, HasGapVector => 1, 
    HasGapEnd5 => 2, HasGapEnd3 => 4, HasGapMaxSpan => 8, HasGapTooManyXs => 16, HasTooShortTr => 32, 
    HasGapEndBig => 64, # upd18 flag for ENDGAP2
    NochangeVecOrNNN => -99,
    }; # sub hasBadGaps

my $NCBICXX=$ENV{ncbicxx}||0; # upd for vecscreen v2/ncbic++ variant;
my $UniVecDB= $ENV{UniVec} || "UniVec_Core"; # Not UniVec
my $MAXGAP=$ENV{maxgap}|| 15; # NCBI decides, changes..

my $ENDGAP20=$ENV{endgap}|| 20; ## was 10; # trim ends if gaps w/in this of ends; NCBI changed again.
# dang ncbi, changed their forbidden gap set again, 2018 it is 
# "more than 5 Ns in the last 10 bases or more than 15 Ns in the last 50 bases"
my $ENDGAP=$ENV{endgap}|| 10; ## 5/10 OR ; was 10; # trim ends if gaps w/in this of ends; NCBI changed again.
my $ENDGAP2=$ENV{endgap2}|| 50; ## 15/50 ; was 10; # trim ends if gaps w/in this of ends; NCBI changed again.
my $NCBIpolicy= 0; #ie DEFAULT ON? strict ncbi tsa policy on gaps, uvectors, etc.

my $MINSIZE_NCBI=200; # NCBI
my $MINSIZE= $ENV{min} || 150; # lower to check cuts w/ valid homology?
my $GAPSOK=1; # 2012-dec .. for now, NCBI may change again (2013-may?) BUT see endgap policy, different
my $NCPU= 1; 
my $tidyup=0;
my $DEFER_UPDATE_FILESET= 1; #... 2014.12upd: make this option!
my $TRANSER_ANNOTS= 1;  # use option? or assume always transfer?
my $ONLYDEGAP=0; # 1807 new opt: asmrna_trimvec.pl -degaponly -mrna my.mrna ; == -novecscreen 

my(%vecscreen);
my($genenames,$cdnaseq,$trclass,$vecscreenf, $logfile, $VECSUF,$outsuf,$namesuf,$parsevecscreen);  ## ,$sufaa,$sufcds,$sufuncut
my(@mergeadd); # upd18

## ARGS
## mixup: locust1all5asm.locust1all5asm.mrna.vector.tab2
$VECSUF= "vector.tab"; #fixed in  vecscreen: makename($cdnaseq,".vector.tab");
$outsuf= $ENV{outsuf}|| "uvcut"; # was uvcut.mrna; 
# upd18: outsuf also now: offcds updates (from merge)
$namesuf= $ENV{namesuf}|| "names";

my $optok= GetOptions(
  "mrna|cdna=s", \$cdnaseq,
  "class|trclass=s", \$trclass,
  "vectors|vecscreen=s", \$vecscreenf,  
  "parsevecscreen=s", \$parsevecscreen,  # parsevec=vecscreen.tmp eg from ncbic++ out
  "names|genenames=s", \$genenames, ## ? allow for 2 files: myspecies.namerefids + allrefprot.names
  "mergeadd=s", \@mergeadd, # -merge xxx.byhand.{mrna,cds,aa,ids} for merge_updates
  "deferupdate!", \$DEFER_UPDATE_FILESET, 
  "logfile:s", \$logfile,
  "MINSIZE=i", \$MINSIZE,  
  "MAXGAP=i", \$MAXGAP,  
  "ENDGAP=i", \$ENDGAP,  
  "NCPU=i", \$NCPU,## "MAXMEM=i", \$MAXMEM,  
  "dryrun|n!", \$dryrun, 
  "degaponly|novecscreen!", \$ONLYDEGAP, 
  "ncbipolicy!", \$NCBIpolicy, 
  "tidyup!", \$tidyup, 
  "DEBUG!", \$DEBUG, 
  );

die "EvidentialGene asmrna_trimvec VERSION ",VERSION,"
mRNA vector screen/removal and gap trimming
  - detects vectors in mRNA transcripts using vecscreen/blastn -d $UniVecDB 
  - trim end gaps per NCBI TSA requirements
  - uses protein/CDS info with ref-protein to prevent spurious/poor screen/trim that damages good mRNA
  - recomputes protein,cds of trimmed mRNA and merges for evgmrna2tsa uses;
  - logs ambigious/problem cases for expert inspection while trying to minimize this,
    tries to keep valid protein while trimming vectors/gaps.
  - ncbipolicy = NCBI current 2018 end-gap policy (should be default?)
  
Usage: asmrna_trimvec.pl -mrna mrna.fasta OR -class name.trclass ...
opts: -genenames mrna.names -degaponly -ncbipolicy -log  -debug
    -NCPU=$NCPU -MINSIZE=$MINSIZE  -MAXGAP=$MAXGAP\n"
  unless($optok and ($cdnaseq or $trclass or $parsevecscreen));     

$tidyup= 1 unless($dryrun||$DEBUG); # default on unless debug|dryrun ?
my $GAPSMAX = ('N') x $MAXGAP;
$TRANSER_ANNOTS= not $DEFER_UPDATE_FILESET; # is this okay? or use option? or assume always transfer?
if($NCBIpolicy) {
  $MINSIZE= $MINSIZE_NCBI if($MINSIZE < $MINSIZE_NCBI); # lower to check cuts w/ valid homology?
}

openloggit($logfile,$cdnaseq||$trclass);
loggit(1, "EvidentialGene asmrna_trimvec.pl VERSION",VERSION);
loggit(1, "ERR: unused arguments:",@ARGV) if(@ARGV>0);

our $APPvecscreen= findapp("vecscreen", 1); # unless($ONLYDEGAP or $parsevecscreen);  
#unused now: our $APPtraa2cds= findevigeneapp("prot/traa2cds.pl"); # move to cdna_evigenesub for get_mRNA
#unused now: our $APPcdnabest= findevigeneapp("cdna_bestorf.pl"); # allow ENV/path substitutions?
#-------------------------------------

##.... REWRITE HERE for single input .mrna + vector.tab + genenames ..........
#cdna_evigenesub: my( %genenames,%genenamepct,%genedbxref,%namedgenes,%cddnames); 

=item steps for asmrna_trimvec 

  0. input.mrna of evigene tr2aacds 
    - requires okayset/name.mrna file creation (revcomp -strand); not done in tr2aacds

  q0? should outputs be relocated to new subdir? eg. tsasubmit/ or leave in okayset/
    - probably leave in okayset/ replacing old versions
    
  1. vecscreen() : run ncbi c-- vecscreen, alternately c++ blastn (later)
        using UniVec_Core (or other) db
  
  2. mrna_trimvec() : remove vector spans in mRNA, and end gaps, recalc orfs of trim set
  
  3. update_mrna_fileset : should this be located in folder or okayset/ ?
      ^^ move to evgmrna2tsa, so can make final pubset merging trimset + okayset + pubids + annot ..
      -- fixme : okboth.aa, .cds may not exist where okboth.mrna is input
      
  Next: evigene/scripts/evgmrna2tsa.pl  now looking for these outputs as input.
  
=cut

=item trimvec preserve hdr annots

  asmrna_trimvec trimset/ 
    kfish2evgxx11.uvcut.PROBLEMCUT	kfish2evgxx11.uvcut.cds		kfish2evgxx11.uvcut.mrna
    kfish2evgxx11.uvcut.aa		kfish2evgxx11.uvcut.ids
  update_mrna_fileset(trimset):  tmpfiles/
    kfish2evgxx11.aa		kfish2evgxx11.cds		kfish2evgxx11.mrna == pubset+trimset merge
    kfish2evgxx11.aa.untrim		kfish2evgxx11.cds.untrim	kfish2evgxx11.mrna.untrim  == orig pubset

trimset/kfish2evgxx11.uvcut.ids
  Funhe2Exx11m012599t5	OKCUT,	uvcut=Strong,27,1-27,end5;	aalen=353,72%,complete; clen=1462; strand=+; offs=343-1404;
  Funhe2Exx11m002577t2	OKCUT,	uvcut=end3trim,28,0-0,end5;	aalen=1108,94%,complete; clen=3519; strand=+; offs=175-3501;
  Funhe2Exx11m050168t1	OKCUT,	uvcut=Moderate,20,1-20,cds5end5;	aalen=119,96%,partial5; clen=372; strand=+; offs=3-362;
tmpfiles/kfish2evgxx11.aa: ** lost annots, need to preserve
  >Funhe2Exx11m012599t5 aalen=353,72%,complete; clen=1462; strand=+; offs=343-1404;; uvcut=Strong,27,1-27,end5;
  >Funhe2Exx11m002577t2 aalen=1108,94%,complete; clen=3519; strand=+; offs=175-3501;; uvcut=end3trim,28,0-0,end5;
  >Funhe2Exx11m050168t1 aalen=119,96%,partial5; clen=372; strand=+; offs=3-362;; uvcut=Moderate,20,1-20,cds5end5;
tmpfiles/kfish2evgxx11.aa.untrim:
  >Funhe2Exx11m012599t5 type=protein; aalen=353,71%,complete; clen=1489; offs=370-1431; oid=Fungr1EG3m010897t1; organism=Fundulus_heteroclitus; evgclass=altmid;
  >Funhe2Exx11m002577t2 type=protein; Name=ACF7 protein; Dbxref=CDD:215849,TrEMBL:UniRef50_Q13696,TrEMBL:Q13696_HUMAN; aalen=1108,93%,complete; clen=3547; offs=175-3501; oid=Fungr1EG3m001095t1; organism=Fundulus_heteroclitus; evgclass=althi;
  >Funhe2Exx11m050168t1 type=protein; aalen=126,97%,partial5; clen=392; offs=2-382; oid=Fungr1EG3m028904t1; organism=Fundulus_heteroclitus; evgclass=main;
     
=cut


sub MAIN_start {}
MAIN: {
  loggit(0, "BEGIN with input=",$cdnaseq||$trclass,"date=",`date`);

  if($parsevecscreen) {
    ## -vectors=$vecscreenf  -parsevec=$parsevecscreen; need both as file names?
    my $vectab= makename($parsevecscreen,".$VECSUF");
    vecscreen_parse($parsevecscreen,$vectab); # special fuction, only do this
    loggit(0, "DONE vecscreen_parse at date=",`date`);
    exit;
  }
  
  my($mrnaseq,$trpath,$trname)= get_evgtrset($trclass,$cdnaseq,$vecoutdir);
    loggit(0, "get_evgtrset=",$mrnaseq,$trpath,$trname); ## facount($cdsseqnr)
    loggit(LOG_DIE, "Missing -mrna",$mrnaseq) unless($mrnaseq and -s $mrnaseq);

  # FIXME here? 2014.12.27: pull Name= other pubset annots from mrna if genenames missing..
  
	unless($genenames) { #? put in get_evgtrset
		my $gnt="$trpath/$trname.names";  $genenames=$gnt if(-s $gnt); 
		loggit(LOG_WARN, "Missing product -names",$gnt) unless($genenames);
	}

  # unless($ONLYDEGAP)  
  ($vecscreenf)= vecscreen( $mrnaseq, $vecscreenf||"", $dryrun||$ONLYDEGAP);


  #here instead of in mrna_trimvec() before/after/with putVecOrNNN():
  # UPD18: add *.offcds.{aa,cds,mrna} update files to @trimfiles for update_mrna_fileset
  #? my($nmcerr,$nupoff,$nrecalc,$mcerrtab)= check_mrna_cds_agree($mrnaseq);
  
  ## mrnaseq = input .mrna; outf, outuncut = updated .mrna << DROP outuncut for update_fileset
  # UPD18: @trimfiles, %$trimids now may be from merge_updates() of several update mrna sets
  # passed to update_mrna_fileset() where %$trimids take precedence over prior input mrna set,
  #  including DROPped seqs not in @trimfiles seq set, but in *.id tables with DROP action
  ## @mergeadd included now here
  
  my($ntrim, $trimids, @trimfiles)  ## $outf, $outaa, $outcds
      = mrna_trimvec($mrnaseq,$vecscreenf,$genenames);

  ## FIXME: add tidyup: remove or move to tmpfiles 
  ##   vecscreen ncpu "_vecsplit/" if ok: -s mrna.vecscreen.tmp or -s mrna.vector.tab
  ## merge new uvcut files + old uvuncut:
  # $nam  .uvcut.mrna + .uvuncut.mrna >> replace input.mrna and flag action 
  # ditto: .uvcut.aa .uvcut.cds ; need pull uncut .aa,.cds from okayset/*
  # need file? of ids for uvcut, uvuncut ?

  #... defer update_mrna_fileset() to mrna2tsa, push @trimfiles to trimset/ dir 
  #... 2014.12upd: make this option!
  #use constant DEFER_UPDATE_FILESET => 1;
  
  my $upstatus=0; 
  my($upstatus,$upfiles,$uptemp,$upokids,$tmpfolder); # my @upfiles=(); $upokids == \%okids
if($DEFER_UPDATE_FILESET) {
	push @tmpfiles, @trimfiles;  ## $uptemp = \@trimfiles; 
	$tmpfolder="trimset";
	$upstatus=1;
	$tidyup= 1; #? always
	
} else {
  ## FIXME in mrna2tsa?  publicset/mrna,aa,cds need >pubid not >oid 
  ## .. rewrite fasta hdr again? or make pubids before this? 
  ## need some adjustments for DROPped mrna to mainalt,pubid 
  
  # NOW: mrnaseq in $vecoutdir, get path for others from it .. but need okayset param
	my $trimflag='trimvec_done';# ($SKIPTRIMSET) ? 'trimvec_SKIP' : 
  ($upstatus, $upfiles, $uptemp, $upokids)  
    = update_mrna_fileset($trpath, $mrnaseq, $trimflag, $trimids, @trimfiles) if(@trimfiles>0); 
	
	## bug fix: dont move upfiles, should be same path as input cdnaseq
	#OFF# push @publicset, @$upfiles if(ref $upfiles); #?? THIS IS PROBLEM
	
	push @tmpfiles, @$uptemp if(ref $uptemp);
	push @tmpfiles, @trimfiles;  ## $uptemp =? \@trimfiles; NO, uvcut == @trimfiles separate
	$tmpfolder="tmpfiles";
}
  
  if( $tidyup and $upstatus > 0 ) { ## not: and -s $upfiles->[2] 
    # FIXME: tidyup also : publicset/name.mrna_vecsplit/  publicset/name.vecscreen.tmp
    my $vtmp= makename($mrnaseq,".vecscreen.tmp");
    push @tmpfiles, $vtmp if(-f $vtmp);
    my $tmpsplit= $mrnaseq."_vecsplit";
    if( -d $tmpsplit) { push @tmpfiles, $tmpsplit; } # remove if empty dir
    
    tidyupFileset($tmpfolder,@tmpfiles);  
	  # $vecoutdir= $trpath; #? solve PROBLEM of next step missing new mrna inputs from this.
    tidyupFileset($vecoutdir,@publicset) if(@publicset);  #? "publicset" not used
    my @rmlist;
    foreach my $fn (@erasefiles) { if(-f $fn) { unlink($fn); push @rmlist,$fn; } }  
    if(@rmlist) { my $nrm=@rmlist; my $rml=join" ",@rmlist[0..4]; loggit(0,"tidyup erase: n=$nrm, $rml .."); } 
  }
  loggit(0, "DONE at date=",`date`);
  loggit(0, "======================================"); # log may append runs
}

#...................................
#... test loop
# my @inmrna= @ARGV;
# foreach my $mf (@inmrna) {
#   (my $nam = $mf) =~ s/\.mrna.*//;
#   my $vtab="$nam.mrna.$VECSUF";
#   my $genenamef="$nam.$namesuf"; # what? evg path is above okayset/mrna, evgmrna2tsa.pl gets it
#   mrna_trimvec($mf,$vtab,$genenamef);
# } # end inmrna

#---------------------------------------------------------------------------


# ($cdnaseq,$trpath,$trname)= get_evgtrset($trclass,$cdnaseq,$outdir);
# from  evigene/scripts/evgmrna2tsa.pl:getmRNA() : should this be here?

sub get_evgtrset {
  my($trclass,$cdnaseq,$outdir)= @_;
  my($trpath,$trname,$nsra,$sradatah)=("","",0,undef);
  my $notokay=0;
  
  if($cdnaseq) { 
    $notokay=1; # dont look in okayset/? look in $outdir now?
    $trclass= makename($cdnaseq,".trclass") unless($trclass); 
  }
  
  if($trclass) { # dont require this exists. just trpath/okayset
    my $trpname= makename($trclass,"","trclass"); 
    if($trpname =~ m,/,) { ($trpath,$trname)= $trpname =~ m,(.*)/([^/]+)$,; } # BADDDDD
    else { $trname=$trpname; }
    $trpath ||= '.';  
     
    ## ?? fixme for update_mrna_fileset : merge okay.aa, okalt.aa,.cds also
    #my $okpath= "$trpath/$outdir" if($outdir); # need to check both, 3 paths?
    my $okpath= ($notokay) ? $trpath :"$trpath/okayset"; # should I check if curdir has okayset files?
    ($cdnaseq)= getmRNA($okpath,$trname,$outdir) if(!$cdnaseq and -d $okpath);
  }
  
  return($cdnaseq,$trpath,$trname,); ## $sradatah);
}

# sub getOkFileset ## moved to cdna_evigenesub.pm
# sub openRead # moved to cdna_evigenesub.pm
# sub getmRNA_OLD   ## moved to cdna_evigenesub.pm  
# sub update_mrna_fileset_OLD  
## moved to cdna_evigenesub.pm, expecting getmRNA:ALSOMAKE_AACDS
## see tr2aacds.pl:asmdupclass_fileset

=item check_mrna_cds_agree

  upd1807: check cds offset + seq agrees in mrna & cds seq files
  wrong for some utrorf, others (e.g. mrna utrorf cut after cds annot)

  ncbi tbl2asn picked up problem cases, mostly utrorf, where cds offset is
  wrong:      406 publicset_we/pig4321ew.tsa.err.ids

  CDSer1 + tsa.err: Susscr4EVm000079t31 Susscr4EVm000168t120 Susscr4EVm000331t135 Susscr4EVm000473t65
  CDSer3 + tsa.err: Susscr4EVm003853t82 Susscr4EVm004297t6 Susscr4EVm008588t23 ..

  
  publicset_we/pigevg4we.mrna_cdserr.tab   
                        m.mlen   m.cdsoff ertype  c.mlen  c.cdsoff oid
  Susscr4EVm000168t120    1026    3085-4107       CDSer1  4107    3085-4107       Susscrtrvelo4b_sRn3l1SRR6236876velvk73Loc10691t10utrorf
  Susscr4EVm000079t31     6559    4338-8048       CDSer1  10893   4338-8048       Susscrtrvelo4b_sRn3l1SRR6236888velvk63Loc8265t1utrorf
  Susscr4EVm003853t82     687     4-531   CDSer3  3880    684-157         Susscrtrvelo1a_sBn1l1SRR6236889velvk45Loc272t2utrorf
  Susscr4EVm004297t6      3672    4-1626  CDSer3  6608    3669-2047       Susscrtrvelo1a_sBn1l1SRR6236889velvk75Loc34t2268utrorf

  some of these are offset-errors in utrorf-cut-mrna, e.g. off-by-1 here
    mrna 4-453 offs maybe is 3-452 of cds, but utrcut at 455, then revcomp makes off-bys common
  Susscr4EVm008595t60     455     4-453   CDSer3  2441    452-3   Susscrtrvelo1a_sBn3l2ERR972387velvk65Loc4531t1utrorf

.. buggy utrorf seqs here?
evg17pig4w/tmpfiles/missparts/                      evg17pig4w/tmpfiles/pig4321ew.cds.missutrorf.ids
evg17pig4w/tmpfiles/pig4321ew.aa.missutrorf         evg17pig4w/tmpfiles/pig4321ew.mrna.missutrorf
evg17pig4w/tmpfiles/pig4321ew.cds.missutrorf        evg17pig4w/tmpfiles/pig4321ew.mrna.missutrorf.ids

method of buggy utrorf.mrna ?
  set pt=pig4321ew
  comm -32 pig4321ew.{cds,mrna}.ids > pig4321ew.cds.missutrorf.ids
  $evigene/scripts/prot/traa2cds.pl -utrorf -trout -aa $pt.aa.missutrorf -cdna $pt.tr -out $pt.mrna.missutrorf
  $evigene/scripts/prot/traa2cds.pl -utrorf -trout -aa $pt.aa.missutrorf -cdna allutrorfmiss.tr -out $pt.mrna.missutrorf2b

orig: .. looks like bug from utrorf-cut of mrna .. wrong offset cut? see below  offs=1990-2439;
>Susscr4EVm008595t60 type=cdna; aalen=150,18%,partial3-utrbad; clen=455; offs=4-453; evgclass=alt; oid=Susscrtrvelo1a_sBn3l2ERR972387velvk65Loc4531t1utrorf;
bestorf: ugh, is revcomp, shorter .. what happened to orig cds?
>Susscr4EVm008595t60 type=protein; aalen=141,92%,partial3; clen=455; strand=-; offs=424-2;  evgclass=alt; oid=Susscrtrvelo1a_sBn3l2ERR972387velvk65Loc4531t1utrorf;

orig full rna:clen=2441; has utrorf Susscr4EVm008595t60
>Susscr4EVm008595t2 type=cdna; Name=ICOS ligand; Dbxref=cow18ncrf:XP_005202100.1,human18nc:XP_006723962.1/45,pig18ncrf:XM_005657167.3/99; 
  aalen=471,58%,complete-utrpoor; clen=2441; offs=32-1447; evgclass=althim; oid=Susscrtrvelo1a_sBn3l2ERR972387velvk65Loc4531t1;

$evigene/scripts/cdna_bestorf.pl  -cdna Susscr4EVm008595t2.mrna
>Susscr4EVm008595t2 type=protein; aalen=471,58%,complete-utrpoor; clen=2441; strand=+; offs=32-1447;  
  Name=ICOS ligand; Dbxref=cow18ncrf:XP_005202100.1,human18nc:XP_006723962.1/45,pig18ncrf:XM_005657167.3/99; evgclass=althim; 
  oid=Susscrtrvelo1a_sBn3l2ERR972387velvk65Loc4531t1;
MGSRGPRSLVLVLPRPGGPGGRGLWGSEGSRGPSSPAARGLRPTMRLRSPRLLLLLFCGL
RAVVAVSQEQEVRAMVGSDVRLACVSLEEGSFDLNDVFVYWQISEPGKANAKSVVTYYLP
ENSSAGHSDNHYRDRAWLSLDSMRQGDFSLHLHNVTPQDEQKFNCLVFRKSLELRKILDV
VVSLHVAANYSMPVVSGPSQDEELTFTCTSTNGYPRPNVYWINRTDGSLLDGALQSSTVS
LNARGLYDVVSVLRIGRAPSVNVGCCIENVLLHQNLTSNQPEMFTGDKKSVPDGPAHDTL
EAXXXXXXXXXXXXXXXTGVPTEATQMPGLRGRSWNLRSSCDESTRPGRGRAHRRDRQQG
HRGGRRRPTGAGGAGAPCPPELQGHGASQPHGLPPAGPADPGEWVLGAAPIPVALMKGQQ
LETCVPKRVGSAPSLSDSEQGGLAGPQGPGRCSASQPCSEAEEGRPVQKEA*

>Susscr4EVm008595t2utrorf type=protein; aalen=150,18%,partial3-utrbad; clen=2441; strand=+; offs=1990-2439; 
MLGSSLAVECFQKVLHRRLLGTHCAHPQENGCVRPAPPRPRLCANQRDGSARSPSQGRSK
KLCVCFLYQKSHTALNACDAGPTRAFNTGALEDSAPERGSGSESRPPARPPGGPPRAPHA
GTLAVWDAMCGRFWKPGFGEYRPSALCQKK
----

  #.. not tsa.err as m.cdsoff fits in m.mlen?                        
  Susscr4EVm096999t1      934     566-931 CDSer1  2037    566-931 Susscrtrvelo1a_sNn9l1ERR789444velvk55Loc23652t1utrorf
  Susscr4EVm042961t2      482     4-480   CDSer3  1061    479-3   Susscrtrvelo1a_sNn9l1ERR789444velvk59Loc51944t1utrorf
        
  # not err, this 691/691 is only problem, from this uvcut annot
  Susscr4EVm109914t1      691/691 302-625 CDSer1  691     302-625 Susscrtrvelo1a_sBn1l1SRR6236889velvk35Loc1331t1

  # not err, NOTE: 'rev' cases by-hand made rev-orfs, are ok-seq, but cds annot is mismatched
  Susscr4EVm005853t4      1281    7-810   CDSer3  2629    1373-93  Susscrtridba2b_sBn2l1SRR1519321idbaidbtk25Loc5932rev

method

/bio/bio-grid/verts/sra2genes/evg17pig4wc
publicset_we/cdsfix.info

# part1: tab of cds<>mrna offs
perl -ne 'if(/^>(\S+)/){ $id=$1; $h=$_; ($tp,$aw,$tw,$ofs)=map{ ($v)=
$h=~m/$_=([^;\s]+)/?$1:0; $v; } qw(type aalen clen offs);
($oid)=m/oid=(\w+)/?$1:"noid";  if(uc($tp) eq "CDS"){
$cann{$id}=[$tw,$ofs,$aw]; } elsif($tp =~ m/^(mRNA|cdna)/i){
if($can=$cann{$id}) { $er=0; $er |=1 if($tw ne $$can[0]); $er |= 2 if($ofs
ne $$can[1]); ($oid)=m/oid=(\w+)/?$1:"noid"; print
join("\t",$id,$tw,$ofs,"CDSer$er",$$can[0],$$can[1],$oid)."\n" if($er); }
else { print "# miss CDS\t$id\t$oid\n";  } } } ' \
  pigevg4we.{cds,mrna}_pub.fa > pigevg4we.mrna_cdserr.tab

# part2: index(mrna,cds) 
perl -ne 'if(/^>(\S+)/){ $d=$1; $infa=1; putfa($id,$fa,$tp) if($fa and
$id);  $id=$fa=""; next unless($tab{$d}); $id=$d;  $h=$_;
($tp,$aw,$tw,$ofs)=map{ ($v)= $h=~m/$_=([^;\s]+)/?$1:0; $v; } qw(type aalen
clen offs); $tp=lc($tp); ($oid)=m/oid=(\w+)/?$1:"noid";  if(uc($tp) eq
"CDS"){ $cann{$id}=[$tw,$ofs,$aw]; } elsif($tp =~ m/^(mRNA|cdna)/i){
if($can=$cann{$id}) { $er=0; $er |=1 if($tw ne $$can[0]); $er |= 2 if($ofs
ne $$can[1]); ($oid)=m/oid=(\w+)/?$1:"noid"; print
join("\t",$id,$tw,$ofs,"CDSer$er",$$can[0],$$can[1],$oid)."\n" if($er and
$ENV{dotab}); } else { print "# miss CDS\t$id\t$oid\n";  } } }
elsif($infa){ if(/^\w/){ chomp; $fa.=$_; } } elsif(/^Sus\w+\t/){
($pd,$tw,$ofs,$cer,$ctw,$cofs,$oid)=@v=split; $tab{$pd}=[@v]; } sub putfa{
my($id,$fa,$tp)=@_; if($tp eq "cds") { $cds{$id}=$fa; } else { my
$cds=$cds{$id}||"nocds"; my $tab=$tab{$id}; my $cb=index($fa,$cds); my
$cw=length($cds); my $tw=length($fa); my $coff= ($cb>=0) ? ($cb+1) ."-".
($cb+$cw) : "misscds";  print join("\t",$id, $tw, $coff, @$tab)."\n";  } } '\
  pigevg4we.mrna_cdserr.tab pigevg4we.{cds,mrna}_pub.fa \
   > pigevg4we.mrna_cdsfix.tab


=cut

# my($nmcerr,$nmcok,$mcerrids)= check_mrna_cds_agree($mrnaf);
# .. use other sub to fix errors, recalc bestorf(mrna) > cds,aa ? handle utrorfs, main problem
sub check_mrna_cds_agree {
  my($inmrna)=@_;
  return 0 unless($inmrna and -f $inmrna);
	my($inaa,$incds,$nmcerr,$nupoff,$nrecalc,)=(0) x 9;
  my $UPDATE_MCSEQ=0; # $nupoff > 0 .. debug, on if nupdateable 
  
  # BUG from uvcut=nocut repeated calls..
  # mrna>Susscr4EVm000316t88 clen=3272/3272/3272/3272/3272; offs=243-2927;...uvcut=nocut,0,0-0,end5; uvcut=nocut,0,0-0,end5; uvcut=nocut,0,0-0,end5; uvcut=nocut,0,0-0,end5; uvcut=nocut,0,0-0,end5;

  #(my $mrnapath= $inmrna) =~ s,[^/]+$,,; 
  # unless($mrnapath) { $mrnapath="./"; } else { $mrnapath =~ s,/$,,; }
  my $mrnapath= dirname( $inmrna);  
  unless($mrnapath) { $mrnapath="./"; } else { $mrnapath =~ s,/$,,; }

  my $mname= basename( $inmrna); 
  # $mname=~s/\.\w+$//; # no good for mrna_pub.fa 
  unless($mname=~s/\.mrna.*//) { $mname=~s/\.\w+$//; }
  my($pubdir)= getFileset($mrnapath);
  
	if($inmrna =~ m/mrna_\w+\./) { # publicset/mrna_(pub|cull|xxx).fa
	  ($incds = $inmrna) =~ s/mrna_/cds_/;
	  ($inaa  = $inmrna) =~ s/mrna_/aa_/;
  } 
  unless($incds and -f $incds) {
    ($incds) = grep /$mname.*cds/, @$pubdir; 
    ($inaa)  = grep /$mname.*aa/, @$pubdir; 
  }
  return 0 unless($incds and -f $incds);
  loggit(0, "check_mrna_cds_agree($inmrna contains $incds)\n"); 
  
  our(%tab,%cann,%cds);
  use constant { tyCDS => 1, tyRNA => 2, tySeqAlign=>3, tySeqUpdate=> 4};
  
  sub putfa {
    my($id,$fa,$tp)=@_; our(%tab,%cds); 
    if($tp == tyCDS) { $cds{$id}= uc($fa); } 
    elsif($tp == tyRNA) { 
      my $cds= $cds{$id}||""; my $tab=$tab{$id}; 
      my $cw= length($cds); my $tw=length($fa); 
      $fa= uc($fa); # BUG input seq case
      my $cb= ($cds) ? index($fa,$cds) : -1; # test revcomp? yes, at least 'rev' set are here
      my $coff= ($cb>=0) ? ($cb+1) ."-". ($cb+$cw) : "misscds";  
      if($cb < 0 and $cds) {
        my $rfa= revcomp($fa);
        my $rcb= index($rfa,$cds);
        if($rcb>=0) {
          $cb=$rcb; # FIXME need rev offs
          $coff= ($cb+$cw)  ."-". ($cb+1); # revc=1
        }
      }
      my($tabtw,$tabcoff)= @{$tab}; # tab = [$tw,$ofs,"CDSer$er",$ctw,$cofs,$oid];
      my $seqer=0; 
      $seqer |= 1 if($tw ne $tabtw); # should be eq
      $seqer |= 2 if($coff ne $tabcoff and $cb>=0);
      $seqer |= 4 if($coff eq "misscds"); # cb == 0
      $seqer |= 8 if($cds eq ""); #? same as misscds? not?
      my @newtab= @$tab; 
      unshift(@newtab, $tw, $coff, "SEQer$seqer"); 
      $tab{$id} = \@newtab;
      } 
  }

  sub readseqhd {
    my($inh,$intype,$seqtest,$outseqh)=@_; 
    our(%tab,%cann); $seqtest||=0;
    my($id,$fa,$infa,$putfa)=(0) x 9;
    while(<$inh>) {
      if(/^>(\S+)/){ my $d=$1; 
        if($seqtest == tySeqAlign) {
          if($fa and $id) { putfa($id,$fa,$intype);  }
          $id=$fa=""; $infa=0; 
          next unless($tab{$d});
          $infa=1;  
        } elsif($seqtest == tySeqUpdate) {
          $putfa=0;
          next unless($tab{$d});
          $putfa=1;  
        }
      
      $id=$d; my $h=$_; 
      my($tp,$aw,$tw,$ofs)=map{ my($v)= $h=~m/$_=([^;\s]+)/?$1:0; $v; } qw(type aalen clen offs);
      (my $twn=$tw)=~s,/\d+$,,; # clen=nnnn/mmmm from last uvcut, drop '/mmmm'
      my($oid)= $h=~m/oid=(\w+)/?$1:"noid";  
      if($seqtest == tySeqUpdate) { # write updated header from tab[ clen offs ]
        my $tval= $tab{$id};
        my ($twb,$sofs,$serr)= @$tval; # my ($twb,$sofs,$serr,$tw,$mofs,$cerr,$ctw,$cofs,$oid)=  
        my $hnew=$h;
        # error unless $soff =~ /\d+-\d+/
        unless($hnew=~s/clen=$tw/clen=$twb/) { $hnew=~s/$/ clen=$twb;/; }
        unless($hnew=~s/offs=$ofs/offs=$sofs/) { $hnew=~s/$/ offs=$sofs;/; }
        print $outseqh $hnew;
      }
      elsif($seqtest == tySeqAlign) { } # have %cann,%tab already, need only id, fa
      elsif($intype == tyCDS){ $cann{$id}=[$twn,$ofs,$aw]; } # uc($tp) eq "CDS"
      elsif($intype == tyRNA) { # $tp =~ m/^(mRNA|cdna)/i
        my($ctw,$cofs,$er)=(0,0,0);
        if(my $can= $cann{$id}) { 
          ($ctw,$cofs)= @$can; # ($$can[0],$$can[1]);
          $er |= 1 if($twn ne $ctw); ### twn NOT tw here
          $er |= 2 if($ofs ne $cofs); 
        } else { $er="miss"; } # or $er |= 4;
        if($er){ 
          #x $tab{$id}= join("\t",$id,$tw,$ofs,"CDSer$er",$ctw,$cofs,$oid);  
          $tab{$id}= [$twn,$ofs,"CDSer$er",$ctw,$cofs,$oid]; # twn NOT tw here
          }
        } 
      } elsif($infa) {
        chomp; $fa.=$_;
      } elsif($putfa) {
        print $outseqh $_;
      }
    }
    if($seqtest == tySeqAlign) { if($fa and $id) { putfa($id,$fa,$intype);  } }
  }
  
  # read headers .. compare cds offset, mrna len
  # also check cds seq contained in mrna seq
  my($ok,$hseq);
  ($ok,$hseq)= openRead($incds);
  if($ok){ readseqhd($hseq,tyCDS); close($hseq); }
  
  ($ok,$hseq) = openRead($inmrna);
  if($ok){ readseqhd($hseq, tyRNA); close($hseq); }
  
  ($ok,$hseq)= openRead($incds);
  if($ok){ readseqhd($hseq,tyCDS,tySeqAlign); close($hseq); }
  
  # now readseqall($hcds, $hrna) and index( cds.seq in mrna.seq)
  # OR just recalc bestorf for each taberr ?
  ($ok,$hseq) = openRead($inmrna);
  if($ok){ readseqhd($hseq, tyRNA, tySeqAlign); close($hseq); }

  my @tabids= sort keys %tab; # only error set, no ok set
  my $oname= makename($inmrna,""); # (my $oname = $mrnaf) =~ s/\.mrna.*//;   
  my $outsuf= "offcds";
  my $outmc= "$oname.$outsuf.ids"; # same as below outd == ids
  my $houtmc= undef;
  my $okmc= (@tabids) ? open($houtmc,'>',$outmc) : 0;
  print $houtmc "#".join("\t",qw(ID Act Oclass sLen sOff sErr mLen mOff cErr cLen cOff OID))."\n" if($okmc);
  
  for my $id (@tabids) {
    my @tval= @{$tab{$id}}; 
    my ($twb,$sofs,$serr,$tw,$mofs,$cerr,$ctw,$cofs,$oid)= (0) x 19;
    if(@tval > 7) { ($twb,$sofs,$serr,$tw,$mofs,$cerr,$ctw,$cofs,$oid)= @tval; }
    else{ ($tw,$mofs,$cerr,$ctw,$cofs,$oid)=@tval; } 
    my $cla="other";
    $nmcerr++;
    # SEQer0 == update_off, only CDS/AA
    # SEQer2 == update_off, both CDS/AA and mRNA
    if($serr eq "SEQer0"){ $cla="update1offs"; $nupoff++; } 
    elsif($serr eq "SEQer2"){ $cla="update2offs"; $nupoff++; } 
    elsif($serr =~ m/SEQer\d/){ $cla="recalc_orfs"; $nrecalc++; }
    # else what? recalc?

    # change to match uvcut.ids : id, action = OKCUT/DROPCUT/PROBLEMCUT
    # cla = OK_update1offs, upd2offs, recalc_orfs ??
    my $act= "OKoffcds";
    print $houtmc join("\t",$id,$act,$cla,@tval)."\n" if($okmc);
    #x print $houtmc join("\t",$id,@tval,$cla)."\n" if($okmc);
  }
  close($houtmc) if($okmc);
  loggit(0, "check_mrna_cds_agree nerr=$nmcerr, nupdate=$nupoff, nrecalc=$nrecalc, table $outmc\n"); 

    ## add here? write changed seq files: header changed, seq changed..
    ## as per uvcut.mra,cds,aa,ids 
  my @outfiles=();
  $UPDATE_MCSEQ=1 if($nmcerr and $nupoff > 0); # recalcs? do we need class recalc/update?
  if($UPDATE_MCSEQ) {
    my($outm,$outa,$outc,$outd)= @outfiles = map{ "$oname.$outsuf.$_" } ("mrna", "aa", "cds", "ids");
    my($outseq,$houtseq);
    my @iof= ( [$inmrna,$outm], [$incds,$outc], [$inaa,$outa] );
    foreach my $iof (@iof) {
      my($inf,$outf)= @$iof;
      my($hout);
      my($ok,$hseq) = openRead($inf);
      $ok = open($hout,'>',$outf) if($ok);
      if($ok){ 
        readseqhd($hseq, tyRNA, tySeqUpdate, $hout);         
        close($hseq); close($hout); 
        }
    }
  }
  # return outnames if have them
  # $nmcerr,$nupoff,$nrecalc,$mcerrtab,$mcoutfiles
  return($nmcerr,$nupoff,$nrecalc,\%tab, \@outfiles);
}

# my($nout,@outf)= merge_updates($oname, \@uvcut, \@offcds);
sub merge_updates
{
  my($oname, @outsets)= @_;
  my @ups= ("mrna", "aa", "cds", "ids");
  my $outsuf= "updates"; # uvcut, offcds, ..
  my($nok)=(0); my @outf; my %mergeids;
  ## update here %trimids? make new as %did
  ## require *.ids? or pick out of seq >id lists? * 2nd way now
  
  for my $ups (@ups) {
    my @ms=(); 
    my $isidtab=($ups eq "ids");
    
    ## check for existance of xxx.ups if not in $oset ; ie -merge xxx.mrna , add xxx.aa,cds,ids here
    ## check for more than $$oset[0] ? may have @merge=(subseta.mrna subb.mrna subc.mrna)
    ## then want to add all the other parts for suba,b,c
    for my $oset (@outsets) { 
      my @m= grep/\.$ups$/, @$oset; 
      unless(@m) {
        for my $osf (@$oset) {
          (my $osu= $osf) =~ s/\.\w+$/.$ups/;  
          push @m, $osu if(-f $osu);
        }
        ## my $osf= $$oset[0] || "nonesuch"; # BUG need all of @$oset -merge list
      }
      push @ms, @m if(@m); 
    }
    next unless(@ms);
    
    my $mout= "$oname.$outsuf.$ups"; 
    my $hmout; my $oko= open($hmout,'>', $mout);
    next unless($oko);
    push @outf, $mout;
    my %did=(); $nok=0 unless($isidtab);
    for my $ms (@ms) {
      my($ok,$hms)= openRead($ms);
      while(<$hms>) {
        if($isidtab){ 
          my($id,$act)= split; 
          $mergeids{$id}=2; # yes? may not exist in seq file but as DROP in id file
          print $hmout $_; 
        } else {
          if(/^>(\S+)/) { my $id=$1; $ok= ($did{$id}++)?0:1; if($ok) { $mergeids{$id}=2; $nok++; } }
          print $hmout $_ if ($ok);
        } 
      } close($hms); 
    } close($hmout);  
    pop(@outf) unless($nok);
    loggit(0,"merge n=$nok of @ms to $mout");
  }
  
  ## $nok= @outf;
  $nok= scalar keys %mergeids;
  return($nok,\%mergeids,@outf); # @ups
}


sub mrna_trimvec 
{  
  my($mrnaf,$vtab,$genenamef)= @_;
  unless ( -f $mrnaf ) { loggit(1, "ERR: uvcut missing mRNA $mrnaf\n"); return -1; }  
  unless ( $ONLYDEGAP or -f $vtab  ) { loggit(1, "ERR: uvcut missing $VECSUF $vtab\n"); return -1; } 
  #FIXME: keep & return id list of cut,uncut for update_mrna_fileset

  #was IN OUT AAOUT CDSOUT OUTUNCUT, $houtuncut
  my $oname= makename($mrnaf,""); # (my $oname = $mrnaf) =~ s/\.mrna.*//; 

  our %trimids=();
  our($hin, $houtf, $houtaa, $houtcds, $hidlist) = (undef) x 9; # file handles replace hard names
  our($outf,$outaa,$outcds,$outids)= map{ "$oname.$outsuf.$_" } ("mrna", "aa", "cds", "ids");
  ## change these to hash of handle/filenames by suffix keys= ("mrna", "aa", "cds", "ids")
  my @outf=($outf, $outaa, $outcds, $outids);
  my @outh=($houtf, $houtaa, $houtcds, $hidlist);
  loggit(0, "uvcut $mrnaf with $vtab to $outf\n"); 
  
  my($namgot,$namin)= parse_genenames($genenamef);  # do in caller?
  loggit(0, "uvcut names $genenamef n=$namgot\n"); 
  
  my $nvecid= readVectab($vtab);  
  loggit(0, "uvcut tr with vector n=$nvecid\n"); 

  # NO/FIXME: add OUTUNCUT for non-uvector set, later combined
  # add gapclean() / trimNNNends() ? from evgmrna2tsa.pl:putseq()
  # .. needs to adjust CDSoffset, maybe CDS if partial with NNN in ENDGAP; easier to recalc prot as w/ uvcut
  
  #add before/after/with putVecOrNNN():  @$mcoutfiles are added via merge_updates()
  my($nmcerr,$nupoff,$nrecalc,$mcerrtab,$mcoutfiles)= check_mrna_cds_agree($mrnaf);

  my $ok=0; 
  ($ok,$hin)= openRead($mrnaf);
  if($ok) { 
    for(my $i=0; $i<@outf; $i++) { $ok= open($outh[$i],'>',$outf[$i]); last unless($ok); } 
    ($houtf, $houtaa, $houtcds, $hidlist)= @outh; # for below... fixme
    }   # DROP# $ok= open($houtuncut,'>',$outuncut); 
  unless($ok) { loggit(1,"ERR: bad files in:$mrnaf out:$outf,.aa,.cds,.."); return -1; }

  ## FIXME TOO MANY 'uvcut=nocut' in .uvcut, from hasBadGaps() == 8 == has $GAPSMAX, now allowed...
  ## FIXME1: 1 case of single 'N' at end5 now is skipped. this came from utrorf.
  ## SEQ_INST.TerminalNs  N at beginning of sequence: RhipulEGm009616t2 
  
  ## FIXME15: add new check just for $dropit=($clen<$MINSIZE) for otherwise clean tr
  ## constant HasTooShortTr => nnn

  our($hasVecOrNNN,$didput,$nmin,$nput,$nskip,$nbadcut,$nnochange)=(0) x 19;
  my($id,$hd,$fa)=(0);
  
  sub cut_vec_or_gap { 
    my($id,$hd,$fa)=@_; 
    our($hasVecOrNNN,$didput,$nmin,$nput,$nskip,$nbadcut,$nnochange);
    our($houtf, $houtaa, $houtcds, $hidlist);
    if($id) { 
      $hasVecOrNNN |= hasBadGaps($fa);
      ## PROBLEM skipping HasGapMaxSpan : CDShasTooManyXs or not GAPSOK
      unless($hasVecOrNNN == HasGapNone) {  ##  or $hasVecOrNNN == HasGapMaxSpan
        $didput= putVecOrNNN($houtf, $houtaa, $houtcds, $hidlist, 
              $id,$hd,$fa,$hasVecOrNNN); 
        if($didput == NochangeVecOrNNN) { # NochangeVecOrNNN == -99 means nochange
          $nnochange++;
        } else {
          $trimids{$id}=1; # regardless of $didput return
          if($didput>0) { $nput++; } elsif($didput<0) { $nbadcut++; } else { $nskip++; }
        }
      } else { 
        $nnochange++; #DROPuncut# print $houtuncut $hd,$fa; 
      }
    } 
  }    

  while(<$hin>) {
    if(/^>(\S+)/) { 
      my $d=$1; $nmin++;
      cut_vec_or_gap($id,$hd,$fa);
      $id=$d; $hd=$_; $fa=""; 
      $hasVecOrNNN= ($vecscreen{$id})? HasGapVector: HasGapNone; 
      } 
    elsif(/\w/) { $fa .= $_; } # dont chomp now, for OUTUNCUT/ASIS
  } 
  cut_vec_or_gap($id,$hd,$fa);

  close($hin); for(my $i=0; $i<@outh; $i++) { $ok= close($outh[$i]); }
  #DROPuncut# close($houtuncut);
  #?? write %trimids to file? should be in .log .. maybe not.

	##add this: 
	if($nbadcut>0) { (my $po=$outids)=~s/ids/PROBLEMCUT/; runcmd("grep PROBLEMCUT $outids > $po"); push @outf, $po; }
  my $nalltrim= $nput + $nskip + $nbadcut; # badcut also counts
  my $trimids= \%trimids;
  
  #* merge ORDER important, assume @mergeadd is best, goes in first.. the uvec @outf, then mcout
  if(@mergeadd or ($nmcerr and @$mcoutfiles)) {
    # push @outf, @$mcoutfiles; # no
    # want to *merge* uvcut.mrna + mcout.mrna, aa, cds, ids for seq update
    # FIXME: %trimids need added mcout ids
    my($nout,$mergeids,@mergeout)= merge_updates($oname, \@mergeadd, \@outf, $mcoutfiles); # opt? -mergeadd xxx*.{ids,aa,cds,mrna} 
    if($nout) { 
      @outf= @mergeout; $trimids= $mergeids; 
      push @tmpfiles, @$mcoutfiles; #UPD1905
      }

    # if($nalltrim == 0) { # but add @mergeadd
    #   @outf= @$mcoutfiles; 
    #   my %mcids= map{ $_ => 2 } keys %$mcerrtab;
    #   $trimids= \%mcids;  # no change in uvcut files
    # } else {
    #   my($nout,$mergeids,@mergeout)= merge_updates($oname, \@outf, $mcoutfiles, \@mergeadd); # opt? -mergeadd xxx*.{ids,aa,cds,mrna} 
    #   if($nout) { @outf= @mergeout; $trimids= $mergeids; }
    # }
  }
  
  my $trimlog="nin=$nmin, nochange=$nnochange, ncut=$nput, ndrop=$nskip, nbadcut=$nbadcut to $outf";
  loggit(0, "uvcut: $trimlog"); 
  return($nalltrim, $trimids, @outf); ## @outf=$outf, $outaa, $outcds, $outids 
} # end sub mrna_trimvec(mrna)


sub hasBadGaps {
  my($fa)= @_;
  $fa =~ s/\s+//g; # do this once not each use?
  my $hasNNN=HasGapNone; # 0
  my $clen= length($fa);
  my $nnn= $fa =~ tr/N/N/; 

  ## FIXME15: add new check just for $dropit=($clen<$MINSIZE) for otherwise clean tr
  $hasNNN |= HasTooShortTr if($clen - $nnn < $MINSIZE); #FIXME15, clen - nnngaps < $MINSIZE
  return $hasNNN unless($nnn);

  # these tests need some NNN
  my $nbig= index($fa,$GAPSMAX); 
  my $n1= index($fa,'N'); 
  my $ne= rindex($fa,'N'); 
  
  # $hasNNN |= HasGapEnd5 if($n1 >= 0 and $n1 < $ENDGAP);
  # $hasNNN |= HasGapEnd3 if($ne >= $clen - $ENDGAP);
  
  # 2018 endgap change .. 5nn/10-end or 15nn/50-end; check endgap2 even if has endgap1?
  if($n1 >= 0 and $n1 < $ENDGAP2) {
    if($n1 >= 0 and $n1 < $ENDGAP) {$hasNNN |= HasGapEnd5; }
    if($nbig >=0 and $nbig < $ENDGAP2) { $hasNNN |= HasGapEndBig + HasGapEnd5; }
    else {
      my $endfa= substr($fa,0, $ENDGAP2); 
      my $nend= $endfa =~ tr/N/N/; 
      if($nend >= $MAXGAP and index($endfa,'NNNNN') >= 0) { $hasNNN |= HasGapEndBig + HasGapEnd5; }
    }
  }
    
  if($ne >= $clen - $ENDGAP2) {
    if($ne >= $clen - $ENDGAP) { $hasNNN |= HasGapEnd3; }
    # got a case of NNN crossing ENDGAP2 bound, ncbi counts NNs past bound
    my $bigend=0;
    if($nbig >=0 ) { my $enbig= rindex($fa,$GAPSMAX); 
      $bigend=1 if ($enbig+$MAXGAP > $clen - $ENDGAP2);
    }
    if($bigend) { $hasNNN |= HasGapEndBig + HasGapEnd5; }
    else { 
      my $endfa= substr($fa, -$ENDGAP2); 
      my $nend= $endfa =~ tr/N/N/; 
      if($nend >= $MAXGAP and index($endfa,'NNNNN') >= 0) { $hasNNN |= HasGapEndBig + HasGapEnd3; }
    }
  }
  
  $hasNNN |= HasGapMaxSpan if($nbig >= 0);
  
  return $hasNNN;
}

sub readVectab
{
  my($vtab)= @_;
  %vecscreen=(); # global now
  my $nvecid= 0;
  my($ok,$hin)= openRead($vtab); #open(F,$vtab) 
  if($ok) {
    while(<$hin>){ next if(/^\W/); my($id,$b,$e,$vt)=split"\t"; $vecscreen{$id} .= "$b\t$e\t$vt\n"; } close($hin);
    ## compress overlaps here ; see also below new vecscreen_parse()
    foreach my $oid (keys %vecscreen) {
      my @vec=split"\n", $vecscreen{$oid};  
      if(@vec>1) { 
        @vec= sort { $a <=> $b } @vec; 
        my @vec2=(); my $ncut=0;
        for(my $i=@vec - 1; $i>0; $i--) { 
           my($ub,$ue,$vty)=split"\t",$vec[$i-1];
           my($xb,$xe,$xty)=split"\t",$vec[$i]; 
           next if($xb<1 or $xe < $xb or $xty =~ /Weak/i); # bad data?
           if($xb <= $ue) { $ub=$xb if($xb<$ub); $vec[$i-1]=join("\t",$ub,$xe,$vty); $ncut++; }  
           else { unshift @vec2, $vec[$i]; }
           }
        unshift @vec2, $vec[0];  $vecscreen{$oid}= join("\n",@vec2);
      }  
    }
  }
  # return %vecscreen;
  $nvecid= scalar(keys %vecscreen);
  return $nvecid;
}

 
sub putVecOrNNN { 
  my($houtf, $houtaa, $houtcds, $houtidlist,
     $oid,$hdr,$fain,$hasVecOrNNN)= @_;
  ## param of outhandles: $houtf, $houtaa, $houtcds  for OUT AAOUT CDSOUT
  #now: hasVecOrNNN & 1 == uvec; & 2 == end5gap; & 4 == end3gap; & 8 == biggap
  ##c{ HasGapNone => 0, HasGapVector => 1, HasGapEnd5 => 2, HasGapEnd3 => 4, HasGapMaxSpan => 8, }; # hasBadGaps

  my $retval=0;
  $fain =~ s/\s+//g;
  my $olen=length($fain);
  
  my $tblinfo= parse_evgheader($oid,$hdr,$olen);
  my $pubid= $tblinfo->{'pubid'};
  my $cdsoff= $tblinfo->{'cdsoff'}; 
  my ($cdsb,$cdse)= split/[-]/,$cdsoff; 
  my $oldaaq= $tblinfo->{'aaqual'};    
  my($aafull)= $oldaaq =~ m/(complete|partial\w*)/; 
  my $aastop=  ($aafull =~ /partial3|partial$/)? 0 : 1; # use as !$aastop == $partial3
  my $aastart= ($aafull =~ /partial5|partial$/)? 0 : 1;
  my $oldorflen = 1 + (($cdsb>$cdse)? $cdsb - $cdse : $cdse - $cdsb);
	my $oldcdsnn=0;
  my $gnamed= $tblinfo->{'name'} || $genenames{$oid} || 0;
  my $namepct= $genenamepct{$oid} || 0; # is structured: 99%,123/345,678
  my $nameref= $genedbxref{$oid} || 0; #  
  my ($okname,$uniqname)= (0,0);
  if($gnamed) { #?? use nameref instead of gnamed to decide uniqness ?? 
    my($npct,$naln)= $namepct=~m/^(\d+)\D+(\d+)/;
    $okname= (not $gnamed or $gnamed =~ /^CDD:|^hypothetical|^uncharacterized|^na$/i or $npct<50) ? 0 : 1;
    my $allids= $namedgenes{$gnamed}||"na"; 
    $uniqname=1 if($okname and $allids eq "$oid,");
  }
        
  
  my($ucut,$rue,$rub)= (0) x 9;
  my($fac,$facgap,$fav,$vectype)=("") x 10;
  
  my $checkNNN= ( $hasVecOrNNN > HasGapVector )?1:0;
  
  if($vecscreen{$oid}) {
    # FIXME: need special case for vector cutting start-codon, insert new ATG (or part missing) in facgap
    ($fac,$facgap,$fav,$ucut,$rub,$rue,$vectype)= cutVector($oid,$fain,$cdsb,$cdse); # if( $hasVecOrNNN & 1 > 0);    
    if($ucut > 0.5*$olen and not $fac) { $checkNNN= 0; $hasVecOrNNN = HasGapVector; } #cut ALL ; dont test NNN ; dont change to fain
    else { $checkNNN= 1; }  # ??? not here FIXME18: recalc hasBadGaps .. cutVector may insert 'nnnnn'
  }
  
  # FIXME: bug DROPS are not being dropped, or kept no change by mistake:
  ## ^^ Problem here, fac becomes blank, but that test below changes to fain : 
  #   sowhiteflyv8k61loc90865t1  UVector! really big uvcut,
  #   uvcut4b #BAD neworf uvtype=Strong match, sowhiteflyv8k61loc90865t1: loss orflen=0-1176; loss named=92%,293/320,391,CDD
  #   uvcut5f sowhiteflyv8k61loc90865t1 uvcut=Strong,1736,1-1736,cdsinend3; clen=1736/1736; << NO CUT

  # FIXME2: lots of ENDGAP in (partial) CDS, but have valid cds/prot before/after those endgaps
  # .. instead of trimming to end and chopping valid cds-bases, should squeeze these gaps down to 3+frame minimum
  ## vectype="endtrim" for non-UVector but endtrim; ..
  ## dang, do we need to trim both fac, facgap ?
  ## FIXME: which uvec vals need endtrim update? ucut? rub,rue?
  my $trimtype=""; my $ntrim= 0; my $nNcut=0; 
  
  ## DEBUG 1807 end gap

  if( $checkNNN ) { #was ($hasVecOrNNN > HasGapVector) 
    my $fadegap= $fac || $fain; # NOT for fac == cutVector entirely
    my $gappy= hasBadGaps($fadegap); # check again vector may have cut
    unless($gappy == HasGapNone) # NO!! or ($GAPSOK && $gappy == HasGapMaxSpan)
    {
      my ($trimgap); 
      ## change this to not endTrim in CDS ! is chopping too many valid bases of partials, but need new cds after vectrim
      ## eg: litovavel3k55Loc3448t5,cdsoffs=2-400: 1-=TGAAATCGTCCGTNNNNNNNNNNNGCTTTGCTA trim>GCTTTGCTACA
      my @cdsvalid= ($fac)? () : ($cdsb,$cdse);
      
      #upd1807: patch partial3 endgap left by endTrimNNN .. aa.end == N, mrna.end extends a bit over that
      #o:($fac,$trimtype,$ntrim,$nNcut,$oldcdsnn)= endTrimNNN($oid,$fadegap,@cdsvalid);
      ($fac,$trimtype,$ntrim,$nNcut,$oldcdsnn)= endTrimNNN($oid,$fadegap, !$aastop,@cdsvalid);
      ($facgap,$trimgap)= endTrimNNN($oid,$facgap, !$aastop) if($facgap);
      $vectype.=$trimtype if($trimtype);
      ## what when ntrim==0 && trimtype == "" ?? nocut
      # $ucut += $ntrim; # below
    }
  }
  ## FIXME: badcut should not count cds-nnn-squeezes, inframe reduction of nnn shouldnt affect prot aligns
  ## .. subtract $nNcut  from tests for bad changes = num NNN cut
  
  my $vecNotStrong = ($vectype =~ /Strong/i)?0:1;
  my $vecNotSqueeze= ($vectype eq "cdsns")?0:1; # special trimtype for cds gaps
  
  ## fixme .. old: rub,rue
  my $clen=length($fac); my $fl=""; 
  if($clen == $olen and not $vectype) { $vectype = "nocut"; } # vectype="nocut" # nochange? uncut?
  elsif(not $vectype) { $vectype = "errcut"; } # what? err?
  
  #BAD# $fac= $fac || $fain; # # NOT for fac == cutVector entirely
  unless($fac =~ /\w/) {  # skip to dropit?? $fac ||= $fain; is this ever right?
    $fl.="allcut"; $clen=0; # clen should == 0
    if( $ucut+$ntrim < 0.5*$olen) { } #problem
  }

  if($rub < $cdse and $rue > $cdsb) { if($rub <= $cdsb+2 and $rue<$cdse) { $fl.="cds5"; }
  elsif($rue >= $cdse-2 and $rub > $cdsb) { $fl.="cds3"; } else { $fl.="cdsin"; } }
  if($rue > $olen-9) { $fl.="end3"; } elsif($rub <= 9) { $fl.="end5"; }
  
  #FIXME2: need bestorf -partial option (always?) so that uvcut doesnt trigger further chop for complete-aa
  #FIXME: here? add protein qual check: cdna_bestorf($fac) : 
  ## if( aasize << cutsize/3? AND vectype = Moderate? and aahomol > minalign) 
  ##   then cancel cutvec; keep orig w/ note; check homol-align also?
  
  my $note=""; my $badcut=0;
  my $expect_neworflen= $oldorflen - $ucut; # not quite right for ucut in UTR
  my ($newaahdr,$newaa,$newcdshdr,$newcds,$newmrnahdr) = ("") x 10;
  my ($newaalen,$newpcds,$newcompl,$neworflen)= (0) x 10; 
  
  ## FIXME here, 2014.12.27: transfer hdr annots into newaahdr,newcdshdr : Name=, other pubset annots..
  ## .. mrna,cds,aa have ~same pub annots; old hdr should have it.
  ## sub getbestorf() here should transfer annots, hdr > newhdr
  
  if($clen > 0) {
  ($newaahdr,$newaa,$newcdshdr,$newcds,$newmrnahdr, $newaalen,$newpcds,$newcompl,$neworflen)
      = getbestorf($oid,$hdr,$fac,$clen, $expect_neworflen);

  if($neworflen < $oldorflen) {
    my $cglen= length($facgap);
    # warn "# fagap.$oid=$cglen,$facgap\n" if $DEBUG;

    if($cglen > $clen) {
    my ($newaahdr1,$newaa1,$newcdshdr1,$newcds1,$newmrnahdr1,
       $newaalen1,$newpcds1,$newcompl1,$neworflen1)= getbestorf($oid,$hdr,$facgap,$cglen, $expect_neworflen);
    my $nga= $newaa1 =~ tr/X/X/; # my $ngc= $newcds1 =~ tr/N/N/;   
    my $newleng= $newaalen1 - $nga;
    # warn "# aacut=$newaahdr; lens=$newaalen1-$newaalen; aacgap.$newaahdr1\n" if $DEBUG;

    if( $newleng > $newaalen) { 
      $note .= "framefixlen=$newleng-$newaalen; ";
      ($newaahdr,$newaa,$newcdshdr,$newcds,$newmrnahdr, $newaalen,$newpcds,$newcompl,$neworflen) = 
        ($newaahdr1,$newaa1,$newcdshdr1,$newcds1, $newmrnahdr1, $newaalen1,$newpcds1,$newcompl1,$neworflen1);
      ($fac,$clen)= ($facgap,$cglen);   
      } 
    }
  }
  ## FIXME-maybe: old aahdr has evgclass tags:  evgclass=main,okay,match:locust1sop4p4k23loc30t48,pct:99/96;
  ## .. copy to new??  in update_mrna_fileset ?
  } # fac/clen
  
  ## ?? here, do after recompute cds  
  my $CDShasTooManyXs=0; my $newcdsnn=0;
  if(1) {
    my($newcdsb,$newcdse)= $newaahdr =~ m/offs=(\d+).(\d+)/; 
    my $cdsfa= substr($fac,$newcdsb-1,1+$newcdse-$newcdsb); ## isnt this $newcds == formatted \n 
    my $cdsw= length($cdsfa);
 		$newcdsnn= $cdsfa =~ tr/N/N/; 
    if($cdsw > 1 and $newcdsnn > 0.48*$cdsw) { # ERR
      $CDShasTooManyXs= $newcdsnn; # flag it for below .. uniqname ?? check below
      $fl .= ",CDShasTooManyXs:$newcdsnn/$cdsw";
      $hasVecOrNNN |= HasGapTooManyXs;
    }
  }
  
  my $cdsnncut= ($oldcdsnn > 0) ? $oldcdsnn - $newcdsnn : 0; # use to adjust badcut tests
	my $neworflentest=$neworflen;  $neworflentest += $cdsnncut if($cdsnncut > 2);
  
  #FIXME here? for uvcut=nocut .. dont keep updating.. skip w/ nochange
  my $noCutNochange= 0;
  if( $vectype eq "nocut" and $fain eq $fac) {
    #?? and newcdsb,e == $cdsb,$cdse
    return(NochangeVecOrNNN); # $noCutNochange=1; # 
  }
  
  #FIXME here for no cut-able gaps?
  #above: use constant NochangeVecOrNNN => -99;  
  if($hasVecOrNNN == HasGapNone or ($GAPSOK && $hasVecOrNNN == HasGapMaxSpan)) { 
    $retval= NochangeVecOrNNN; # == -99  # return and print asis
    return($retval);
  } 
  
  loggit(0, "getbestorf: $oid oldaa=$oldaaq; new.$newaahdr\n"); 
        
  if($neworflen > $oldorflen) { # DEBUG check these
    # my $namepct= $genenamepct{$oid} || 0; # is structured: 99%,123/345,678
    # my $nameref= $genedbxref{$oid} || 0; #  
    # my $named= $genenames{$oid} || 0; #  same as $gnamed
    my($npct,$naln)= $namepct=~m/^(\d+)\D+(\d+)/;
    my $odiff= $neworflen-$oldorflen;
    loggit(0, "named:  $oid LONGER dorf=+$odiff, named=$namepct,$nameref,$gnamed;\n") if($npct>=50 or $nameref =~ /:/);  
  }
  
  if($neworflentest < $oldorflen - 9) { # check named align; also debug check names for neworflen >> oldorflen ?
    my $odiff= $neworflentest-$oldorflen;
    
    ## ?? ignore 'CDD:' as main name, it lacks related species homolog
    ## .. not sure ignore CDD: is right yet.
    my($npct,$naln)= $namepct=~m/^(\d+)\D+(\d+)/;
    loggit(0, "named:  $oid dorf=$odiff, named=$namepct,$nameref,$gnamed;\n") if($npct>=50 or $nameref =~ /:/); 
    
    my($minpct,$minaln)= ($vecNotStrong) ? (66,50) : (90,90);
    ## for vecIsStrong, want to measure size of uvcut vs size of align? really need to know if uvcut is align-supported
    ## but dont yet have ref-align spans as inputs
    
    #FIXME: adjust badcut to accept Strong when also partial-cds at cut end?)
    # .. eg: shrimpt nbadcut=26 but 20 are uv=Strong, cut ~30 cds; probably accept all Strong uvcut
    # .. also have some perfect huge vector matches w/ complete prot: NUBP1 = Ecoli vec
    if($okname and not $vecNotStrong) {
      if($oldaaq =~ /partial/ and $neworflentest >= 0.95 * $expect_neworflen) { $okname=0; }
      elsif($oldaaq !~ /partial/  and $ucut > 99) { $okname=0; } # this is a long-Strong uvec, named likely vector protein
    }    
    $okname=0 if(not $vecNotSqueeze); # if($vecIsSqueeze); # not badcut for these, I hope..
     
    ## ?? need levels of npct test here, dont know unless npct ~100 if uvcut affects alignment
    
    if($okname and $npct >= $minpct and $naln >= $minaln and $neworflentest/3 < 0.99*$naln) { ##  ??
      $fl .= ",refalignlosscut:$odiff";
      $badcut++; $note.="loss named=$namepct,$nameref,$gnamed; "; # d=$odiff, 
    }
  }

  ## Disallow  $vectype =~ /Strong/ here?   ; skip this if refalignlosscut ?
  ## FIXME: add note for any largish change in orflen, smaller/bigger, whether uvec or trimNNN
  if(not $badcut and ($newaalen<1 or $neworflentest < 0.90 * $expect_neworflen)) {
    my $odiff= $neworflentest-$oldorflen; # neworflentest or neworflen ??
    $fl .= ",shortorfcut:$odiff"; # note even for Strong? yes. 
    if($vecNotStrong and $vecNotSqueeze) { $badcut++; $note.="loss orflen=$neworflentest-$oldorflen; "; }
  }
  
  $ucut += $ntrim; #?
  my $idflag="";  
  my $uvhdr="uvcut=$vectype,$ucut,$rub-$rue,$fl;";
  if($badcut) {  
    $idflag.="PROBLEMCUT,";
    # complain, maybe cancel uvcut; not badcut if $vectype=~/Strong/ always?
    loggit(1, "BAD neworf $uvhdr $oid: $note\n");
    $retval= -1;
  }
  
  ## FIXME1412: bad mRNA cds off= ** transfer from newaahdr: offs|aalen ..
  ## update mRNA hdr from $newaahdr : 
  # $hdr=$newmrnahdr; # chomp($hdr); 
  #bad# $newmrnahdr =~ s,clen=,clen=$clen/,; #?**BUG clen=nnn/mmm/ppp/qqq ; transfer to $newmrnahdr
  $newmrnahdr =~ s,\bclen=[^;\s]+,clen=$clen,; # fix, add clenold=$oldclen?

  my $uvh2=$uvhdr; 
  if($fav and $DEBUG) {
    if(length($fav)>45) { $fav=substr($fav,0,20)."..".substr($fav,-20); } ## dont stick LONG fav in header ?
    $uvh2.=" uvfa=$fav;" ; ##?? add uvfa= only for DEBUG ? drop uvcut= for pubset ?
  }
  
  ## BUG subtle, if hdr has ONLY >ID, append "; xxx" changes ID !!! was: $newmrnahdr=~s/$/; $uvh2/;
  $newmrnahdr=~s/$/ $uvh2/; # fixme, not 1st, at end.. was: s/ / $uvh2 /;
  my $nnn= $fac =~ tr/N/N/; ## which is MINSIZE limit ? w/ or w/o nnn?
  $fac =~ s/(.{60})/$1\n/g; $fac.="\n" unless($fac=~/\n$/);
  
  ## MINSIZE should change w/ named value, as per evgmrna2tsa.pl:putseq() 
  ## unique(name)/strong homol-name should keep shorter clen; but w/ annotation to support short len
  ## FIXME: dropit conflict with uniqnamed genes. from mrna2tsa

  my $dropit= ($clen - $nnn < $MINSIZE or $CDShasTooManyXs)?1:0;  
  
  if($dropit and $uniqname) { 
    if($CDShasTooManyXs) { } # drop anyway, uniqname is suspect?
    elsif(($newaalen > 15 and $newcompl =~ /complete/) or $newaalen > 20) { $dropit=0; }
    my $dval= ($dropit)? "DROP.":"KEEP.";
    $dval .= ($CDShasTooManyXs) ? "gaps:$CDShasTooManyXs" : "short:$clen-$nnn";
    $idflag.="PROBLEMCUT.uniquename:$gnamed,"; $idflag.="error:$dval," unless($dropit);
    loggit(1,"PROBLEM: unique name '$gnamed' but error:$dval"); 
  }
  
  if($dropit) { 
    my $dval= ($CDShasTooManyXs) ? "gaps:$CDShasTooManyXs" : "short:$clen-$nnn";
    $idflag .= "DROPCUT.$dval,";
    loggit(1, "DROPCUT.$dval: $hdr\n"); 
    $retval= ($badcut)? -1 : 0;  # change from 0 to ?? -2
  } else { 
    $idflag.="OKCUT," unless($idflag =~ /CUT/);
    # $houtf, $houtaa, $houtcds, .. change to hash of handles?  $houts{mrna}, $houts{aa} ..
    print $houtf ">$oid " unless($newmrnahdr =~ /^>/);
    print $houtf   $newmrnahdr,"\n",$fac;  #?? has old/ newhdr >ID ??
    print $houtaa  ">$oid $newaahdr; $uvhdr\n$newaa\n";    
    print $houtcds ">$oid $newcdshdr; $uvhdr\n$newcds\n";   
    $retval= ($badcut)? -1 : 1;
  }
  print $houtidlist  join("\t",$oid,$idflag,$uvhdr,$newaahdr)."\n"; # if($houtidlist);
  return $retval;
}    


=item UVCUT > ORFloss too big: bad bestorf
    
# ** Problem w/ getbestorf() only gets new starts at ATG/M, but cut over start can damage this
# .. need getbestorf() to look for partial5 after each stop codon.. patch here when cut thru ATG startcodon
# .. should add back that plus cds-inframe gap:  s/$uvectorwithstart/UtrnnnATGnnCds/; check 'cds5' flags for cdsloss > uvcut

grep BAD log.uvcut5f | perl -ne'($uvc)=m/uvcut=\w+,(\d+)/; ($oc,$ob)=m/orflen=(\d+).(\d+)/; $oloss
1=$ob-$oc; ($oloss)=m/cut:-(\d+);/; $oloss||=$oloss1; print if($oloss > 1.1*$uvc); ' | head
#BAD neworf uvcut=Moderate,39,144-182,cds5,refalignlosscut:-810; socatfishv1k95loc34193t1: loss named=100%,450/450,446,CDD:191181,DRERI:ENSDARG00000052408,UniProt:A4IG58,,Vertebrate mannosyl (Alpha-1,6-)-glycoprotein beta-1,2-N-acetylglucosaminyltransferase (MGAT2, zgc:162268); 
#BAD neworf uvcut=Moderate,50,134-1779,cds3end3,shortorfcut:-219; socatfishv1k25loc6796t5: loss orflen=1356-1575; 

.. not too many cases ..
/bio/bio-grid/aabugs4/tsaevgc/tsaoutz

cacao3all7f/log.uvcut5f
#BAD neworf uvcut=Moderate,26,429-454,cdsin,shortorfcut:-120; cacao3vel14sc9k35Loc1726t4: loss orflen=489-609; 
  >> has cdsgaps, problem likely is those gaps, not cut but mangle new orf call, same prot start.
#BAD neworf uvcut=Moderate,26,695-720,cdsin,shortorfcut:-105; cacao3sopcsc2k29loc6700t1: loss orflen=522-627; 
  >> ditto, NNNN cdsgaps trail uvcut, mangle new orf call
  
catfish1all4cf/log.uvcut5f
#BAD neworf uvcut=Moderate,39,144-182,cds5,refalignlosscut:-810; socatfishv1k95loc34193t1: loss named=100%,450/450,446,CDD:191181,DRERI:ENSDARG00000052408,UniProt:A4IG58,,Vertebrate mannosyl (Alpha-1,6-)-glycoprotein beta-1,2-N-acetylglucosaminyltransferase (MGAT2, zgc:162268); 
  below>> DEFINITELY cancel uvcut socatfishv1k95loc34193t1 (or replace ATG start in cutspan)
#BAD neworf uvcut=Moderate,50,134-1779,cds3end3,shortorfcut:-219; socatfishv1k25loc6796t5: loss orflen=1356-1575; 
  below>> SHOULD ok this uvcut; >> 2 uvcuts, damage cds, which also has NNN spans;

whitefly1evgcf/log.uvcut5f
#BAD neworf uvcut=Moderate,25,308-332,cds5,shortorfcut:-207; whitefly1vel6k39Loc14435t1: framefixlen=40-24; loss orflen=129-336; 
  >> uvcut is in startcodon-base3, replace 'G' corrects this overcut
  orig>  aalen=111,37%,complete-utrpoor; clen=906; strand=+; offs=306-641;
  G+cut> aalen=104,35%,complete-utrpoor; clen=885; strand=+; offs=306-620;  same prot
  G-cut> aalen=42,14%,complete-utrbad; clen=885; strand=+; offs=259-387;    diff prot
  
banana1all3cf/log.uvcut5f
litova1all3f/log.uvcut5f
locust1evgcf/log.uvcut5f
pogonus1all3cf/log.uvcut5f
shrimpt1evgf/log.uvcut5f
zticktr2acf/log.uvcut5f

  end5.squeezeNNN
        ($cdsfac,$ccut)= squeezeNNN($cdsfac,0,'NNNN',0,2*$ENDGAP);
  2018 ncbi gap policy:
  more than 5 Ns in the first 10 bases or more than 15 Ns in the first 50 bases
  
=cut

sub squeezeNNN {
  my($fa,$nlower,$gapsmax,$keepnnn,$inmax)= @_;
  my $ncut=0;
  my $gapw= length( $gapsmax);
  # keepnnn = 0 or 3 for squeeze keep
  my $NNNs= 'NNNNNNNNNNN'; my $N1='N';
  if($nlower) { map{ $_= lc($_) } ($NNNs,$N1,$gapsmax); } # NOT USED here
  my $dorev=0; if($inmax < 0) { $fa=reverse($fa); $dorev=1; $inmax= -$inmax; }

  for (my $in= index($fa,$gapsmax); $in >= 0; $in=index($fa,$gapsmax)) {
    last if($inmax>0 and $in>=$inmax);
    my $w=length($fa); my $en=$in+$gapw; 
    $en++ while($en<$w and substr($fa,$en,1) eq $N1); 
    my $wn= $en-$in; 
    my $keep= $keepnnn + ($wn % 3); # PROBLEM: got +10 nnn in cases; keep shold not be that big
    $keep=3 if($keep==0);
    my $cut= $wn-$keep; $ncut+=$cut; 
    my $facut= substr($fa,0,$in).substr($NNNs,0,$keep).substr($fa,$en); 
    $fa=$facut; 
  } 
  if($dorev) { $fa=reverse($fa); }
  return($fa,$ncut);
}

#    ($fain2,$trimtype)= endTrimNNN($oid,$fain2);
sub endTrimNNN {
  my($oid,$fa,$partial3,$cdsb,$cdse) = @_;
  my $trimtype="";
  my $ncut=0; 
  my $olen=length($fa);
  $fa =~ s/n/N/g;
  my $nNold= $fa =~ tr/N/N/; my $cdsNold= 0;
  
  my $EndGapThis= ($NCBIpolicy)? $ENDGAP2 : $ENDGAP; #?? upd18 this way or not
  
  $partial3||=0;
  ## drop cds stuff? or not? dont want endtrim over valid cds bases
  my $validcds= (defined $cdse and $cdse>0)?1:0;
  my ($lcdsb,$lcdse)= ($validcds)?($cdsb,$cdse):(0,0);

=item upd1807 partial3 cdsend-N gap fix

 upd1807: patch partial3 endgap left by endTrimNNN .. aa.end == N, mrna.end extends a bit over that
  uvcut=cdsnsend3trim,31,0-0,end5;  uvcut=cdsnsend3trim,13,0-0,end5; uvcut=nocut,0,0-0,end5; uvcut=cdsnsend3trim,28,0-0,end5;
  uvcut=Strong,63,1881-1943,cds3end3;  uvcut=Strong,63,1392-1454,cdsinend3; uvcut=Strong,63,620-682,cdsinend3;
  now fixed:
  uvcut=cdsnsend3trim,31,0-0,end5; uvcut=cdsnsend3trim,2,0-0,end5;<< fix2 here
  uvcut=Strong,63,1881-1943,cds3end3; uvcut=cdsnsend3trim,2,0-0,end5;

  < >Susscr4EVm006998t15 type=cdna; Name=Protein SET, partial; Dbxref=human18nc:NP_003002.2,pig18ncrf:NM_001244090
  .1/99; aalen=257,66%,partial3; clen=1156/1156; offs=385-1155; evgclass=alt; oid=Susscrtrvelo1a_trvelo2b_sBn1l1SR
  R6236876velvk45Loc4618t6; uvcut=cdsnsend3trim,28,0-0,end5; uvcut=nocut,0,0-0,end5;
  < GATGAGGATGAAGGNC
  
  upd trims 'NC' end gap, CDS-N ending
  > >Susscr4EVm006998t15 type=cdna; Name=Protein SET, partial; Dbxref=human18nc:NP_003002.2,pig18ncrf:NM_001244090
  .1/99; aalen=256,66%,partial3; clen=1154/1156; offs=385-1152; evgclass=alt; oid=Susscrtrvelo1a_trvelo2b_sBn1l1SR
  R6236876velvk45Loc4618t6; uvcut=cdsnsend3trim,28,0-0,end5; uvcut=cdsnsend3trim,2,0-0,end5;
  > GATGAGGATGAAGG

=cut

  ## ?? FIXME for in-validcds, ncRNA .. or is it ok here
  if($validcds) {
    my $fixit=0;
    my $cdsfa= substr($fa,$cdsb-1,1+$cdse-$cdsb); 
    my $cdsw=length($cdsfa);
		$cdsNold= $cdsfa =~ tr/N/N/;
		my $fixone=0;
    my $partial5= (substr($cdsfa,0,3) eq 'ATG')?0:1;
		
    my $NNs= ($partial3)?'N':'NN'; #upd1807: patch partial3 endgap
    my $ne= rindex($fa,$NNs); 
    ## NOT SEEN for problem cases..
    ## if($partial3 and $ne>=0){ loggit(0, "partial3NN $oid $cdsb-$cdse, Ne=$ne\n");  }

    #upd18: need ENDGAP2 checks
    
    if($ne >= $olen - $EndGapThis and $olen - $cdse < 3*$EndGapThis) {
      my $cne= rindex($cdsfa,$NNs); 
      # if($partial3 and $cne>=0){ loggit(0, "partial3NN $oid $cdsb-$cdse, Ne=$ne, cdsNe=$cne\n");  }
      if($partial3 and $cne >= 0 and $cdsw - $cne < 3) { $fixone |=2; $fixit |=2; } # 1807 test 
      if($cne >=0 and $cdsw - $cne < 2*$EndGapThis) { 
        #o# my $endc= substr($cdsfa,0,-$EndGapThis); # this is 1..450 for cdsfa = 1..500
        #? bug^^ 
        my $endc= substr($cdsfa,-$EndGapThis);  # this is 451..500 for cdsfa = 1..500
        $endc =~ s/[^N][^N]/Z/g;  
        my $zz= $endc =~ tr/Z/Z/;
        $fixone |=2; $fixit |=2 if($zz >= 3);#?? is this zz test bad? bad if fixit 1 or 2 and other zz fails, need fix both
      }
    }

    ## maybe partial5 treat same way as partial3 end gaps?
    my $n1= index($fa,'NN'); 
    if( $n1 >= 0 and $n1 <= $EndGapThis and $cdsb < 3*$EndGapThis) {
      $n1= index($cdsfa,'NN'); 
      if($n1 >=0 and $n1 < 2*$EndGapThis) { 
        my $endc= substr($cdsfa,0,$EndGapThis);
        $endc =~ s/[^N][^N]/Z/g;  my $zz= $endc =~ tr/Z/Z/;
        $fixone |=1; $fixit |=1 if($zz >= 3); #?? is this zz test bad?
        if($partial5 and $n1 < $EndGapThis) { $fixone |=1; $fixit |=1; } # 1809 test 
      }
    }
    ## zz test is to chop end NNN unless have enough non-N for squeeze.
    if($fixit and $fixone) { $fixit= $fixone; } # ignore zz test if one end passes
    elsif($fixone and $partial3) { $fixit |= $fixone; } #??
    # or $fixit= $fixone if($fixone); ??
    
    if($fixit) { 
      my $cdsfac= $cdsfa;
      my $ccut=0;
      ## fixme: endgap only for 1..end5gap, need reverse(cdsfa) for end3gap
      ## maxgap 4 here, must be bigger than squeezemax of 3
      
      ## BAD here ?? got no SqueezeNNN .. long 'nnnnnnnn' spans in cdsns same as orig mRNA.untrim
      ## loggit(ERR) if ccut == 0
      
      if($fixit & 2 ) {## == 2
        ($cdsfac,$ccut)= squeezeNNN($cdsfac,0,'NNNN',0,-2*$EndGapThis);
        $ncut += $ccut;
      }
      if($fixit & 1) { ##  == 1
        ($cdsfac,$ccut)= squeezeNNN($cdsfac,0,'NNNN',0,2*$EndGapThis);
        $ncut += $ccut;
      }
      $cdsfac =~ s/^N+//; $cdsfac =~ s/N+$//;
      $trimtype.="cdsns"; # $partial3 end gaps .. see this now
      $cdsfac =~ s/N/n/g; # lower to prevent following trim
      # if($partial3) dont append  substr($fa,$cdse); ??
      my $faend= ($partial3) ? "" : substr($fa,$cdse);
      #  my $fabeg= ($partial5) ? "" : substr($fa,0,$cdsb-1); # maybe
      $fa= substr($fa,0,$cdsb-1). $cdsfac . $faend;
    }
    
  }
        
    ## not right; NN may start in UTR, end in CDS
    # my $incds= ( $n1 >= 0 and $n1 <= $EndGapThis and $n1 >= $cdsb)
    #  or ($ne >= $olen - $EndGapThis and $ne <= $cdse);

    # //no//FIXME4: No NN End trim in CDS for aastart,aastop: keep stop/start codon 
    # FIXME: cdsb,cdse adjust for inner gaps
    # fixme2: must adjust cdsb,e when cut BEFORE cdsb
    # YES: fixme3: this is a mess; better to a. cut NNN, b. rerun cdna_bestorf for new cds offset?
    # FIXME: single end N happens.. need 
 
  #upd18: need ENDGAP2 checks
   
  $fa=~s/(N+)$//;
  my $curlen= length($fa); 
  my $NNs='NN';
  #x: my $NNs= ($partial3)?'N':'NN'; #upd1807: patch partial3 endgap
  my $ne= rindex($fa,$NNs);
  for(my $iter=0; $ne >= $curlen - $EndGapThis and $curlen > $EndGapThis and $iter<5; $iter++) {
    $fa= substr($fa,0,$ne); 
    if($fa=~s/(N+)$//) {  my $ncut=length($1); $ne-=$ncut; }
    $curlen= length($fa); 
    $ne= rindex($fa,$NNs);  
  }
  $trimtype.="end3trim" if($curlen < $olen);
  
  ## FIXME: cds-phase/codon_start changes w/ mod 3 of n1   
  ## FIX2: see above  nnnAnnn < need to recurse chop endgaps ?
  my $ol2= length($fa); 
  $fa=~s/^(N+)//; 
  $curlen= length($fa);  
  my $n1= index($fa,'NN'); 
  for( my $iter=0; $n1 >= 0 and $n1 <= $EndGapThis and $iter<5; $iter++) {
    $n1++; $fa= substr($fa,$n1);  
    if($fa=~s/^(N+)//) { my $ncut=length($1); $n1+=$ncut; }
    $curlen= length($fa); 
    $n1= index($fa,'NN');
  }
  $trimtype.="end5trim" if($curlen < $ol2);

  unless($GAPSOK) {
    my($fac,$cut)= squeezeNNN($fa,0,$GAPSMAX,3,0);
    $fa= $fac; $ncut+=$cut; 
  }
  
  my $nNnew= $fa =~ tr/N/N/;
	my $nNcut = $nNold - $nNnew;
  $ncut = $olen - length($fa);
  return($fa,$trimtype,$ncut,$nNcut,$cdsNold);
}


=item endTrim tests
    
#   if($validcds and ($cdsb < 3*$ENDGAP or $olen - $cdse < 3*$ENDGAP)) {
#     my $cdsfa= substr($fa,$cdsb-1,1+$cdse-$cdsb); my $cdsw=length($cdsfa);
#     my $n1= index($cdsfa,'NN'); 
#     my $ne= rindex($cdsfa,'NN'); 
#     
#     ## this isnt right yet... now other mistake: chopping too many valid cds bases
#     ## unfixit if %N > 50%? or >80%? at ends
#     if($n1 >=0 and $n1 < 2*$ENDGAP) { 
# ## 5d:    
#       my $endc= substr($cdsfa,0,$ENDGAP);
#       $endc =~ s/[^N][^N]/Z/g;  my $zz= $endc =~ tr/Z/Z/;
#       $fixit |=1 if($zz >= 3);
# ## 5c:      
# #       my $endc= substr($cdsfa,0,2*$ENDGAP);
# #       my $nn= $endc =~ tr/N/N/; 
# #       $fixit |=1 unless($nn >= $ENDGAP);
#      }
#     if($ne >=0 and $cdsw - $ne < 2*$ENDGAP) {
# ## 5d:    
#       my $endc= substr($cdsfa,0,-$ENDGAP);
#       $endc =~ s/[^N][^N]/Z/g;  my $zz= $endc =~ tr/Z/Z/;
#       $fixit |=2 if($zz >= 3);
# #       my $endc= substr($cdsfa,-2*$ENDGAP);
# #       my $nn= $endc =~ tr/N/N/; 
# #       $fixit |=2 unless($nn >= $ENDGAP);
#     }
#       
#     if($fixit) { 
#       my $cdsfac= $cdsfa;
#       my $ccut=0;
#       ## fixme: endgap only for 1..end5gap, need reverse(cdsfa) for end3gap
#       ## maxgap 4 here, must be bigger than squeezemax of 3
#       if($fixit & 2 == 2) {
#         ($cdsfac,$ccut)= squeezeNNN($cdsfac,0,'NNNN',0,-2*$ENDGAP);
#         $ncut += $ccut;
#       }
#       if($fixit & 1 == 1) {
#         ($cdsfac,$ccut)= squeezeNNN($cdsfac,0,'NNNN',0,2*$ENDGAP);
#         $ncut += $ccut;
#       }
#       $cdsfac =~ s/^N+//; $cdsfac =~ s/N+$//;
#       $trimtype.="cdsns";
#       $cdsfac =~ s/N/n/g; # lower to prevent following trim
#       $fa= substr($fa,0,$cdsb-1). $cdsfac . substr($fa,$cdse);
#     }
#   }

## problem 5c: chopping too many valid cds bases; should squeeze this instead.
#BAD neworf uvcut=end5trim,34,0-0,end5,refalignlosscut:-33; litovavel2rk35Loc20449t1: loss named=72%,144/200,148,CDD:201192,UniRef50_B1PT29,UniProt:B1PT29_ARTSF,,Cuticle protein; 
# uncut>litovavel2rk35Loc20449t1 type=cdna; aalen=148,99%,partial; clen=447;  strand=+; offs=2-445;
# CCCTGCTTNNNNNNNNNNNNNNNNNNNNNNNNNN GCCCCTGCTCCTGCTTACAAAGCCCC
# cut5c>litovavel2rk35Loc20449t1 uvcut=end5trim,34,0-0,end5,refalignlosscut:-33; type=cdna; aalen=148,99%,partial; clen=413/447;  strand=+; offs=2-445;
#  GCCCCTGCTCCTGCTTACAAAGCCCCTGAGCCTACCTACTCTGCCCCTTCCCCTAGCTAC

## problem 5b case: mostly NNN at cds end
##BAD neworf uvcut=cdsnsend5trim,26,0-0,end5,refalignlosscut:-27; litovavel1k25Loc19454t5: loss named=75%,308/411,311,CDD:201393,UniRef50_Q9GZS9,UniProt:CHST5_HUMAN,,Carbohydrate sulfotransferase 5; 
# uncut>litovavel1k25Loc19454t5 type=cdna; aalen=311,99%,partial; clen=937;  strand=+; offs=3-935;
# ANNNNNNNNNNNNNNNNNNNNNNNANNNNNN AGCAACGGCCGACGAAGGAGAGAAAGCCA
# cut5b>litovavel1k25Loc19454t5 uvcut=cdsnsend5trim,26,0-0,end5,refalignlosscut:-27; type=cdna; aalen=311,99%,partial; clen=911/937;  strand=+; offs=3-935;
# nAnnn AGCAACGGCCGACGAAGGAGAGAAAGCCAACAGCACCGAGAGCATCATCGCCTCC
#...
#BAD neworf uvcut=cdsnsend5trim,29,0-0,end5,refalignlosscut:-30; litovavel2rk21Loc21431t4: loss named=95%,346/365,342,CDD:173624,UniRef50_E9Q3W1,UniProt:E9Q3W1_MOUSE,,Casein kinase I; 
# uncut>litovavel2rk21Loc21431t4 type=cdna; aalen=342,99%,partial; clen=1029;  strand=+; offs=3-1028;
# CNNNNNNNNNNNNNNNNNNNTNNNNNNNNNNNNNNN CTTCCTCCTCTTCCCTCCTCTCCG
# cut5b>litovavel2rk21Loc21431t4 uvcut=cdsnsend5trim,29,0-0,end5,refalignlosscut:-30; type=cdna; aalen=342,99%,partial; clen=1000/1029;  strand=+; offs=3-1028;
# nnnTnnn CTTCCTCCTCTTCCCTCCTCTCCGCCAAGATGTCTTCGGGAATCATGGGGTGC

=cut
  
# sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
# sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }

sub cutVector {
  my($oid,$fain,$cdsb,$cdse) = @_;
  my $olen=length($fain);
  my($fac,$facgap,$fav,$vectype,$lue,$lub,$ucut,$rue,$rub);
  $fac=$facgap=$fav=$vectype=""; $lue=$lub=$ucut=$rue=$rub=0; 
  return($fac,$facgap,$fav,$ucut,$rub,$rue,$vectype) unless($vecscreen{$oid});
    
  my @vec=split"\n", $vecscreen{$oid}; my $nv=@vec;
  
## moved up....
#  ## fixmed: not wrong order, but overlapped;
#  #MISORDER catfishvel3ik45Loc262t5 uvec 1/2, 528-561 .. 549-567
#  #MISORDER catfishvel3ik45Loc513t4 uvec 1/2, 1-27 .. 15-34
  
  # FIXMEd: special case for vector cuts start-codon, keep ATG (or part missing) in facgap
  my $hasstart= ($cdsb>0 and uc(substr($fain,$cdsb-1,3)) eq 'ATG')?1:0;
  
  for(my $i=0; $i<$nv; $i++) { 
    my($ub,$ue,$vty)=split"\t",$vec[$i]; 
    
      # special case: retain ATG; need also handle uvcut in middle of codon..
      # BUT there are uvectors w/ proteins, drop-all is right answer, this way preserves only 3-base cds and is dropped
      #DROP too short clen=3: >sowhiteflyv8k61loc90865t1 uvcut=Strongc1,1733,1-1736,cdsinend3,shortorfcut:-1176; uvfa=TTCTTATCTCCTTTTGTAGT..TGGCGAGCTGGATGATGAGC; type=cdna; aalen=391,67%,complete; clen=3/1736;  strand=+; offs=132-1307;
      # n=12 cases in whitefly.
      
    my $overstart=($hasstart and $ue >= $cdsb and $ub <= $cdsb+2)?1:0; # special case: retain ATG in facgap
    $overstart=0 if($overstart and $ub <= $cdsb and $ue >= $cdse); # skip prot contained in vector
    if($overstart) {
      my $ue1= $cdsb-1; my $ub2= $cdsb+3; 
      my @add=();
      if($ub < $ue1) { push(@add, "$ub\t$ue1\t$vty"."c1"); }
      if($ub2 < $ue) { push(@add, "$ub2\t$ue\t$vty"."c1"); }
      if(@add) { 
        splice(@vec,$i,1,@add); # replace [$i] + add 1 
        $nv = @vec; 
        ($ub,$ue,$vty)=split"\t",$vec[$i]; # continue on w/ new vec[i]
      } else { next; }
     # now push these into @vec to replace overstart?
    }  
    
    my $uw=1+$ue-$ub; 
    $ucut += $uw;  ## add ucutInCDS += nb for ue - ub in ce - cb
    $vectype= $vty unless($vectype); # take 1st only?
    # $vectype .= $vty unless($vectype =~ /$vty/); # or this?
    if($i == 0) {  $rub=$ub; if($ub>1) { my $s=substr($fain,0,$ub - 1); $fac.=$s; $facgap.=$s; } }
    else { 
      my $wnu= $ub - 1 - $lue; if($wnu>0) { my $s= substr($fain,$lue,$wnu); $fac.=$s; $facgap.=$s; }
      loggit(1, "MISORDER $oid uvec $i/$nv, $lub-$lue .. $ub-$ue\n") if($ub < $lue);
      }
    
    ##?? add UVGAP only if in cds-span? if($ub < $cdse and $ue > $cdsb) 
    ## For ~3 cases, this is WRONG, leave out UVgap and get better new aa
    ## In 1 case, moderate uvec cuts cds5 makes much worse prot; 
    ## socatfishv1k95loc34193t1 (origaa=446,full and 100% match zfish gene; cutaa=176,utrbad;)
    ## vecscreen finds ~39 bp align w/ 2+ mismatch in cds5span
    ## w/o uvcut, full align to alpha-1,6-mannosyl-glycoprotein 2-beta-N-acetylglucosaminyltransferase [Danio rerio]
    ## .. use mrna.ann.txt homol-align vs vecsreen-medium to decide not to cut?
    ## FIXmaybe: some vector spans surrounded by NNN ; cut those?

    ## try both, use neworf to pick best : seems to work  # if($UVGAP)
    if(1) { 
      my $nug=3 + ($uw % 3); my $uvg= substr("nnnnnnnnn",0,$nug); $facgap .= $uvg; 
    }
    
    if($i == $nv-1) { $rue=$ue; if($ue < $olen) { my $s= substr($fain,$ue); $fac.=$s; $facgap.=$s; } }
    $fav.="n" if($i>0 and $ub>$lue+1); $fav .= substr($fain,$ub-1,$uw); 
    $lue=$ue; $lub=$ub; 
  }

  ##? do this here or wait for endtrimNNN
  $fac =~ s/[Nn]+$//; $fac =~ s/^[Nn]+//; #? yes
  $facgap =~ s/[Nn]+$//; $facgap =~ s/^[Nn]+//; #? yes
  $vectype =~ s/ match//g;  

  return($fac,$facgap,$fav,$ucut,$rub,$rue,$vectype);
}
    

sub getbestorf {
  my($id, $hdr,$cdnain,$cdnasize, $oldorfcutlen)= @_;

  ## FIXME here, 2014.12.27: transfer hdr annots into newaahdr,newcdshdr : Name=, other pubset annots..
  ## .. mrna,cds,aa have ~same pub annots; old hdr should have it.
  ## sub getbestorf() here should transfer annots, hdr > newhdr
  # old:>Funhe2Exx11m002577t2 type=protein; Name=ACF7 protein; Dbxref=CDD:215849,TrEMBL:UniRef50_Q13696,TrEMBL:Q13696_HUMAN; aalen=1108,93%,complete; clen=3547; offs=175-3501; oid=Fungr1EG3m001095t1; organism=Fundulus_heteroclitus; evgclass=althi;
  # new:>Funhe2Exx11m002577t2 aalen=1108,94%,complete; clen=3519; strand=+; offs=175-3501;; uvcut=end3trim,28,0-0,end5;
  # keepan: Name|Dbxref|oid|organism|evgclass
  # replace: aalen|clen|strand|offs ;  uvcut added by caller.
  
  my @keepan=();
  if($TRANSER_ANNOTS) { 
    @keepan= grep{ not m/(type|aalen|clen|strand|offs)=/ } $hdr =~ m/\b(\w+=[^;\n]+)/g; 
    }
  ## update mRNA hdr from $newaahdr : return here as newmrnahdr
  my $newmrnahdr= $hdr; chomp($newmrnahdr);
  
  # my $cmd="$APPcdnabest -nostop -cdna $cdnaseq -aaseq $aaseq -cdsseq $cdsseq"; # -minaa=$MINAA 
  my $MINSIZE4CDS= 60; # what?
  
  if($cdnasize > $MINSIZE4CDS) {
  ## from cdna_bestorf.pl:cdna_bestorf()
  my $fullpart= "best.fwdstrand"; # bothstrand or only fwdstrand ? this is mrna, oriented
  ## ^^ need longorf vs fullorf option
  ## -- use if uvcut-long-partial > uvcut-best but < origaa (dont want > origaa)
  
  # FIXME: set bestorf:
  my $savenostop= $NoStopCodon; $NoStopCodon=1; # for proteindoc()
  
  my( $bestorf)= getBestProt2( $fullpart, $cdnain);  
  # my( $orfprot, $orient)= orfParts($bestorf, [qw(protein orient)]);
  if(ref($bestorf)) {
    my $orient= undef;
    my $asfasta=1; # ($action =~ /fasta/)
    my($aalen,$pcds,$compl,$orflen,$orfprothdr,$orfprotfa)
          = proteindoc($bestorf,$cdnasize,$orient,$asfasta);
    # warn "best.fwdstrand: $id $aalen aa, cdnasize=$cdnasize\n" if($DEBUG);
          
    if($orflen < 0.98*$oldorfcutlen and $compl =~ /complete/) {
      my($longorf)= getBestProt2( "long.fwdstrand", $cdnain);  # ask for longest partial orf
      my($aalen1,$pcds1,$compl1,$orflen1,$orfprothdr1,$orfprotfa1)
           = proteindoc($longorf,$cdnasize,$orient,$asfasta);
    # warn "long.fwdstrand: $id $aalen1 aa, cdnasize=$cdnasize\n" if($DEBUG);
           
     if($orflen1 > $orflen) {
      $bestorf= $longorf;
      ($aalen,$pcds,$compl,$orflen,$orfprothdr,$orfprotfa)=
        ($aalen1,$pcds1,$compl1,$orflen1,$orfprothdr1,$orfprotfa1);
     }      
    }
    
    if(@keepan) {
      my %ohdr= map{ $_ => 1 } $orfprothdr =~ m/\b(\w+)=[^;\n]+/g;
      @keepan= grep{ my($k)=split"="; not $ohdr{$k}; } @keepan;
      $orfprothdr.= " ".join "; ",@keepan; 
    } 

    ## FIXME1412: bad mRNA cds off= ** transfer from newaahdr: offs|aalen ..
    for my $nk (qw(aalen strand offs)) {  # clen or not?
      my($nkv)= $orfprothdr =~ m/\b($nk=[^;\s]+)/;  $nkv||=""; #?
      $newmrnahdr =~ s/\b$nk=[^;\s]+/$nkv/;
    }

    $NoStopCodon= $savenostop; # for proteindoc()
       
    my $cdsfa= $bestorf->{sequence};  $cdsfa  =~ s/(.{60})/$1\n/g;
    (my $cdshdr=$orfprothdr) =~ s/ / type=cds; /;
    ## NOTE: >id is not in these hdr; BUT for input hdr/newmrnahdr
    return($orfprothdr,$orfprotfa,$cdshdr,$cdsfa,$newmrnahdr,
       $aalen,$pcds,$compl,$orflen);          
    }
  }
  return("","","","",$newmrnahdr,0,0,0,0);          
}


## move from  evigene/scripts/evgmrna2tsa.pl
sub vecscreen
{
  my($cdnaseq,$vectab,$skiprun)=@_;
  $vectab= makename($cdnaseq,".$VECSUF") unless($vectab);
  ## look various folder for data?? trimset? publicset?  vecoutdir??
  unless(-f $vectab) { my($dt,$ft)= getFileset($vecoutdir,$VECSUF); $vectab=$ft if($ft); 
    unless(-f $vectab) { my($dt,$ft)= getFileset("trimset",$VECSUF); $vectab=$ft if($ft); }
    unless(-f $vectab) { my($dt,$ft)= getFileset('.',$VECSUF); $vectab=$ft if($ft); }
    }
  return($vectab) if( -s $vectab or $skiprun); # or dryrun .. or ONLYDEGAP
  
  # my($id,$vb,$ve,$ty,$vd,$outh,$inh);  
  ## FIXME: can we use ncbic++ instead? output not same as c-vecscreen... need to check curr ncbi source
  ## ncbic++ doesnt yet support vecscreen .. need own blastn -db UniVec parser to match.. but ncbic- seems obsolete
  
  ##  $ncbi/bin/vecscreen
  unless($NCBICXX) {
    (my $ncbicxxapp=$APPvecscreen) =~ s,/vecscreen,/makeblastdb,; 
    $NCBICXX=1 if( -x $ncbicxxapp);
  }
  
  my $univecdb="";
  (my $ncbid=$APPvecscreen) =~ s,/vecscreen,/..,; ## want option for db UniVec path
  if( -f "$UniVecDB.nsq") {
    $univecdb= "$UniVecDB"; 
  } elsif( -f "$ncbid/data/$UniVecDB.nsq") {
    $univecdb= "$ncbid/data/$UniVecDB"; 
  } else {
    loggit(1,"ERR: $APPvecscreen missing ../data/$UniVecDB.nsq"); return; 
  }
  
  ## lots of this warn: [vecscreen] WARNING:  [000.000]  Blast: No valid letters to be indexed on context 0
  #old# my $ok= open($inh,"$APPvecscreen -i $cdnaseq -d $univecdb -f3 |");
  
  my ($cmddone,$err)=(0,0);
  my $vectmp= makename($cdnaseq,".vecscreen.tmp");
  my $veclog= makename($cdnaseq,".vecscreen.log");
  
  my $cmd0="$APPvecscreen -f3 -d $univecdb"; # add -i cdna.split1.fa -o cdna.split1.vec 2> log
  my $cmd1= $cmd0 . " -i $cdnaseq -o $vectmp 2> $veclog";
  if($NCBICXX) { #NCBI Cxx version has diff opts
    $cmd0="$APPvecscreen -text_output -outfmt 0 -db $univecdb"; # add -query input.fa -out name.vec
    $cmd1= $cmd0 . " -query $cdnaseq -out $vectmp 2> $veclog";
  }
  # if($VECuseBLASTN) { # and $NCBICXX ; later..
  #  $cmd0="$APPblastn -task blastn -outfmt 7 -db $univecdb -penalty -5 -gapopen 4 -gapextend 4 -soft_masking true -dust yes -evalue 700 -searchsp 1750000000000 ";
  #  $cmd1= $cmd0 . " -query $cdnaseq -out $vectmp 2> $veclog";
  #}
  
  if(-s $vectmp) { $cmddone=1; } # skip runcmd if have old data
  elsif( $NCPU > 1 ) {
    my $ccount= facount($cdnaseq); # use this, not fasize
    if($ccount >= 50*$NCPU) {   
      ($err)= vecscreen_ncpu($NCPU,$cmd0,$ccount,$cdnaseq,$vectmp); # cat outparts.
      $cmddone=1; 
    }
  } 
  unless($cmddone) {
    $err= runcmd($cmd1); # "$cmd0 -i $cdnaseq -o $vectmp 2> $veclog"
  }

  my($nvid)= vecscreen_parse($vectmp,$vectab);
  if($nvid<0) { loggit(1,"ERR: $APPvecscreen -i $cdnaseq -d $univecdb TO $vectab"); return; }
  return($vectab);
}


=item NCBICXX vecscreen format

   #>Vector Funhe2Exx11m011522t4 Screen: VecScreen Database: UniVec_Core201404  # v1 starts this way
   Database: UniVec_Core (build 6.0) # new v2
   Strong match             # .. same here + filler lines to ignore
   3       33
   Suspect origin
   1       2
   34      73
   Query= Funhe2Exx11m127091t3 type=mRNA; aalen=103,44%,complete-utrpoor; # new v2, get q.ID here
   > gnl|uv|NGB00362.1:1-61 Illumina Paired End PCR Primer 2.0            # new v2, get hit.ID here
   
=cut


sub vecscreen_parse
{
  my($vectmp,$vectab)= @_;
  my($id,$vb,$ve,$ty,$vd,$outh);  
  my($nvid, $so, $indb, @vbe, %vtbe); $nvid=0; 
  my $nvs1=not $NCBICXX; my $nvs2=$NCBICXX; ## V2: Revise for NCBIXX outfmpt ..

use constant VPARSE3 => 1; # 201501 BUG fix for multiple hits..
    %vtbe=(); #VPARSE v3
    @vbe=();  #VPARSE v2;

  sub putv2 { # VPARSE2
    my($outh,$id,$ty,$vd,@vbe)= @_; 
    return 0 unless($id and $ty and not $ty=~/Weak/i);
    my $no=0; foreach my $vbe (@vbe) { my($b,$e)= @$vbe;
      print $outh join("\t",$id,$b,$e,$ty,$vd)."\n"; $no++;
    } return $no;
  }

=item bugs here adding strong + susp + weak

  >Vector Funhe2Exx11m015740t1 Screen: VecScreen Database: UniVec_Core201404
  Strong match
  1731    1762
  Weak match
  1766    1782
  Suspect origin
  1763    1765

=item compress overlap spans? NO, let cutVector() handle

  Funhe2Exx11m134113t1
  Strong match
  18	51
  Moderate match
  1	20
  
  Funhe2Exx11m134113t1	1	20	Moderate match	UniVec_Core201404
  Funhe2Exx11m134113t1	18	51	Strong match	UniVec_Core201404    # 18b < 20e
  
  >Vector Funhe2Exx11m045311t1
  Moderate match
  8	30
  59	83
  140	165
  Weak match
  168	191  
  Suspect origin
  1	7
  31	58
  166	167
  
  Funhe2Exx11m045311t1	1	58	Moderate match	UniVec_Core201404
  Funhe2Exx11m045311t1	59	83	Moderate match	UniVec_Core201404
  #Funhe2Exx11m045311t1	168	191	Weak match	UniVec_Core201404,Added1 to Moderate match
  Funhe2Exx11m045311t1	140	191	Moderate match	UniVec_Core201404

=cut
  
  sub putv3 { # VPARSE3
    my($outh,$id,$vd,$vtbe)= @_; 
    my $no=0; 
    my @vt = sort keys %$vtbe; # have Suspect origin in ty? NO, already applied below
    my @vts= grep(/Strong|Moderate/,@vt);  my($wt)= grep(/Weak/,@vt);
    return 0 unless($id and @vts);
     ## add Weak if it hits edges of Strong/Mod/Suspect : treat like SuspectOrigin **
    if( @vts and $wt ) {
      my $wadd=0; my $st="";
      my($wb,$we,@wvbemore)= split",",$vtbe->{$wt};      
      for my $t (@vts) {
        my @vbe= split",",$vtbe->{$t}; 
        for(my $i=0; $i<@vbe; $i+=2) { 
          if($we == $vbe[$i]-1) { $vbe[$i]= $wb; $wadd++; $st=$t; last; }
          elsif($wb == $vbe[$i+1]+1) { $vbe[$i+1]= $we; $wadd++; $st=$t; last; }  
          }
        $vtbe->{$t}= join",",@vbe;  
        }
      print $outh "#".join("\t",$id,$wb,$we,$wt,$vd.",Added$wadd to $st")."\n" if($wadd>0);  # log it??
    }
    
    for my $ty (@vts) { 
      my @vbe= split",",$vtbe->{$ty}; #  can be list of (b1,e1),(b2,e2) ??
      while( my($b,$e)=splice(@vbe,0,2)){ print $outh join("\t",$id,$b,$e,$ty,$vd)."\n"; $no++; }
    }
    return $no;
  }

  my($ok,$hin)= openRead($vectmp);  
  $ok= open($outh,'>',$vectab) if($ok);
  return(-1) unless($ok); # unless($ok) { loggit(1,"ERR: $APPvecscreen -i $cdnaseq -d $univecdb TO $vectab"); return; }

  while(<$hin>) {
    chomp; 
    if($nvs1) { $indb=1; }
    elsif($nvs2) { 
      if(/^Database: (\S+)/) { 
if(VPARSE3) { $nvid += putv3($outh,$id,$vd,\%vtbe) if($id and %vtbe);  }
       else { $nvid += putv2($outh,$id,$ty,$vd,@vbe) if(@vbe and $id); }
        $indb=$vd=$1; $id=0;  %vtbe=();
      } elsif(/^Query= (\S+)/) { $id=$1; $indb=0; } 
    }

    if($nvs2 and not $indb and /^>/) { # V2: vec hit, last needed info
      my($vid)=m/>\s*(\S+)/; $vid=~s/gnl\|//; $vid=~s/\|/:/; $vd.=":$vid" if($vid);
if(VPARSE3) { $nvid += putv3($outh,$id,$vd,\%vtbe) if($id and %vtbe);  }
       else { $nvid += putv2($outh,$id,$ty,$vd,@vbe) if(@vbe and $id); }
      $vb=$ve=$so=$ty=0; @vbe=(); %vtbe=();
      
    } elsif($nvs1 and /^>/) {  # V1: this is next hit, >Vector mytrid..
if(VPARSE3) { $nvid += putv3($outh,$id,$vd,\%vtbe) if($id and %vtbe);  }
       else { $nvid += putv2($outh,$id,$ty,$vd,@vbe) if(@vbe and $id); }
      ($id)=m/>Vector (\S+)/; ($vd)=m/Database: (\S+)/; 
      $vb=$ve=$so=$ty=0; @vbe=(); %vtbe=();
      
    } elsif(/^No hits/ or /No hits found/) {  ## V1 & V2
      $id=0; 
      
    } elsif($indb and /^\w/) {  
      if(/ match/) {  # V2: indb b/n ^Database and ^Query
        $ty=$_; $so=0; # BUG here, need multiple match types, and associated spans.
        if(VPARSE3) { $vtbe{$ty}="" unless($vtbe{$ty}); }
      } elsif(/^Suspect origin/) { 
        $so=1;  # So follows all Match, need to append following b,e to appropriate @vbe offby1 of match range.
        
      } elsif(/^(\d+)\s+(\d+)/) { 
        my($b,$e)=($1,$2); 
        
        if($e < $b) { } # bug? skip or warn
        elsif($so) {  
if(VPARSE3){ ## v3: use all types in %vtbe
          for my $t (sort keys %vtbe) {
            my @vbe= split",",$vtbe{$t}; #  list of (b1,e1),(b2,e2) ??
            for(my $i=0; $i<@vbe; $i+=2) { 
              if($e == $vbe[$i]-1) { $vbe[$i]= $b; $b=$e=-6;  last ; }
              elsif($b == $vbe[$i+1]+1) { $vbe[$i+1]= $e; $b=$e=-9; last ; } 
              }
            $vtbe{$t}= join",",@vbe;  
            }
} else { ## v2:
          for(my $i=0; $i<@vbe; $i++) {
            my($mb,$me)= @{$vbe[$i]};
            if($e == $mb-1) { $vbe[$i]= [$b,$me]; last; }
            elsif($b == $me+1) { $vbe[$i]= [$mb,$e]; last; } # ($b >= $me and $b <= $me+2) is range possible ?
          }
}
        } elsif($ty) { 
if(VPARSE3){ $vtbe{$ty}.="," if($vtbe{$ty}); $vtbe{$ty}.= "$b,$e"; } # or push hash, [$b,$e] ??
      else { push @vbe, [$b,$e]; }
        } 
         
      } 
    }
  } 
  
if(VPARSE3) { $nvid += putv3($outh,$id,$vd,\%vtbe) if($id and %vtbe);  }
       else { $nvid += putv2($outh,$id,$ty,$vd,@vbe) if(@vbe and $id); }
  close($outh); close($hin);

  loggit(0, "vectors found in ntr=",$nvid,$vectab); 
  return($nvid); # return($vectab);
}


sub vecscreen_ncpu
{
  my($npart,$cmd0,$ccount,$cdnaseq,$vecout)=@_;
  
  # my $ccount= facount($cdnaseq); # use this, not fasize
  my $splcount= int(0.99 + $ccount/$npart);
  my $spldir= makename($cdnaseq,"_vecsplit/",'cdna|fasta|fsa|fa');  # use _tsasubmit/ instead?
  mkdir($spldir); # dryrun?
  
  # push @tmpfiles, $spldir;  ## tmpfiles or erasefiles ?
  
  ## my $err= runcmd("$APPvecscreen -f3 -d $univecdb -i $cdnaseq -o $vectmp 2> $veclog");
  my @splset= fasplitcount( $cdnaseq, $spldir, $npart, $splcount,"fa"); 
  my @vecset;
  my $npartgot= @splset;
  my $icpu= 0;   my $err=0;
  
  #NO: chdir($spldir); ## BAD!! chdir("../"); ?
  ## forkCMD= /home/ux455375/bio/ncbic11/bin/vecscreen -f3 -d /home/ux455375/bio/ncbic11/bin/../data/UniVec_Core\
  #  -i ./okayset/litova1all3.mrna_vecsplit/litova1all3.mrna.split28.fa \
  #  -o ./okayset/litova1all3.mrna_vecsplit/litova1all3.mrna.split28.vecout \
  #  2> ./okayset/litova1all3.mrna_vecsplit/litova1all3.mrna.split28.veclog
  ## errlog: ./okayset/litova1all3.mrna_vecsplit/litova1all3.mrna.split22.veclog : No such file or directory
  
  for(my $ip=0; $ip< $npartgot; $ip++) {
    my $cdna1= $splset[$ip];
    push @erasefiles, $cdna1;
    (my $veco1=$cdna1) =~ s/\.\w+$/.vecout/;
    (my $dlog1=$cdna1) =~ s/\.\w+$/.veclog/;
    push @vecset, $veco1;
    my $cmd1= $cmd0 . " -i $cdna1 -o $veco1 2> $dlog1";
    if($NCBICXX) { #NCBI Cxx version has diff opts
       $cmd1= $cmd0 . " -query $cdna1 -out $veco1 2> $dlog1";
    }

    my $pid= forkcmd($cmd1);    
    if(++$icpu > $npartgot) { while (wait() != -1) { }; $icpu= 0; }
  }
  while (wait() != -1) { };
  
  # my $cmd= "cat ".join(' ',@vecset)." > $vecout"; runcmd($cmd);
  # preserve dlogs also?
  open(my $outh,'>',$vecout);
  foreach my $vf (@vecset) { 
    if(open(FT,$vf)) { while(<FT>) { print $outh $_; } close(FT);} 
    
    ## if($err||$DEBUG||$dryrun) { push @erasefiles, $vf; } else { unlink $vf if(-f $vf); }
    push @erasefiles, $vf;
    $vf=~s/vecout/veclog/; push @erasefiles, $vf;
    }
  close($outh);
  # rmdir($spldir) unless($err||$DEBUG||$dryrun); ## if($nerase == $npartgot);
  return($err); # $vecout
}

__END__


=item OPTION 201501 add ncbi c++ blastn parsing, as vecscreen app hard to find..
    
    blastn1 vecs opts:  nucmis -q -5, gapopen -G 3, gapextend -E 3, filter -F "m D", expect -e 700, searchsp -Y 1.75e12
    blastn2 vecs opts: (gap open/ex 4 because 3 not allowed..)
    $nbinp/blastn -penalty -5 -gapopen 4 -gapextend 4 -soft_masking true -dust yes -evalue 700 -searchsp 1750000000000 \
        -task blastn -db $nbin/../data/UniVec_Core -query Funhe2Exx11m044139t3.mrna.untrim

    # vecscreen1 ...............
    head -30 tsa5kfish2rae5x11.uvechitmiss.mrna.vecscreen2
    >Vector Funhe2Exx11m127720t5 Screen: VecScreen Database: UniVec_Core201404
    Moderate match  # eval=1.4
    1202    1221
    Weak match      # eval=81
    1       17 
    >Vector Funhe2Exx11m017337t2 Screen: VecScreen Database: UniVec_Core201404
    Moderate match  # eval=2e-06
    40      73
    Weak match      # eval=21, but same vector as mod, other half
    4       21
    Suspect origin  # should cut full span 1..73 of uv|NGB00362
    1       3
    22      39
    >Vector Funhe2Exx11m011522t4 Screen: VecScreen Database: UniVec_Core201404
    Weak match      # eval=81
    1704    1728
    >Vector Funhe2Exx11m038259t1 Screen: VecScreen Database: UniVec_Core201404
    Moderate match  # eval=1.4
    446     465
    Weak match
    1       17
    >Vector Funhe2Exx11m015740t1 Screen: VecScreen Database: UniVec_Core201404
    Strong match    # eval=3e-08
    1731    1762
    Weak match
    1766    1782

  # blastn v2 equiv ...............
  # .. Strong/Moderate cutoff in e-value? b/n e-08 and e-06  ?;   Mod/Weak cut b/n eval 1.4 and >10?
  % ggrep -A56 '^# Query' tsa5kfish2rae5x11.uvechitmiss.mrna.blastn2t | egrep -v 'Datab|Field' | head -30
  # Query: Funhe2Exx11m127720t5 aalen=335,82%,complete; clen=1221; offs=81-1088; pubid=Funhe2EKm001196t1; oid=Fungr1EG3m008237t5; Name=Glyceraldehyde-3-phosphate dehydrogenase, testis-specific (100%M); genegroup=FISH11G_G896.s1; Dbxref=kfish2gene:Funhe2EKm001196t1,CDD:223135,TrEMBL:UniRef50_G3HPT5,TrEMBL:G3HPT5_CRIGR,family:FISH11G_G896;
  # 2 hits found
  Funhe2Exx11m127720t5    gnl|uv|NGB00362.1:1-61  100.00  20      0       0       1202    1221    61      42      1.4     40.2
  Funhe2Exx11m127720t5    gnl|uv|NGB00362.1:1-61  100.00  17      0       0       1       17      45      61         81   34.3
  # BLASTN 2.2.29+
  # Query: Funhe2Exx11m017337t2 aalen=358,63%,partial5; clen=1701; offs=3-1079; pubid=Funhe2EKm003428t1; oid=Fungr1EG3m008299t4; Name=DnaJ subfamily C member (97%P); genegroup=FISH11G_G4189; Dbxref=kfish2gene:Funhe2EKm003428t1,CDD:223560,TrEMBL:UniRef50_Q9H1X3,TrEMBL:DJC25_HUMAN,family:FISH11G_G4189;
  # 3 hits found
  Funhe2Exx11m017337t2    gnl|uv|NGB00362.1:1-61  97.06   34      1       0       40      73      28      61      2e-06   59.9
  Funhe2Exx11m017337t2    gnl|uv|NGB00362.1:1-61  100.00  18      0       0       4       21      40      57         21   36.3
  Funhe2Exx11m017337t2    gnl|uv|M28829.1:1-8684-49       100.00  16      0       0       1591    1606    1307    1292      315   32.4
  # BLASTN 2.2.29+
  # Query: Funhe2Exx11m011522t4 aalen=418,72%,complete; clen=1728; offs=7-1263; pubid=Funhe2EKm003823t1; oid=Fungr1EG3m006328t3; Name=Arrestin domain-containing protein 3 (100%P); genegroup=FISH11G_G1836.s2; Dbxref=kfish2gene:Funhe2EKm003823t1,CDD:215866,TrEMBL:UniRef50_Q96B67,TrEMBL:ARRD3_HUMAN,family:FISH11G_G1836;
  # 2 hits found
  Funhe2Exx11m011522t4    gnl|uv|NGB00362.1:1-61  100.00  17      0       0       1712    1728    61      45         81   34.3
  Funhe2Exx11m011522t4    gnl|uv|J02459.1:1-48502-49      100.00  16      0       0       1704    1719    36763   36778     315   32.4
  # BLASTN 2.2.29+
  # Query: Funhe2Exx11m038259t1 aalen=155,100%,partial; clen=465; offs=1-465; pubid=Funhe2EKm004243t1; oid=Fungr1EG3m023720t1; Name=Microtubule-associated proteins 1A/1B light chain 3A (98%P); genegroup=FISH11G_G6021; Dbxref=kfish2gene:Funhe2EKm004243t1,CDD:176355,TrEMBL:UniRef50_A6NCE7,TrEMBL:MP3B2_HUMAN,family:FISH11G_G6021;
  # 20 hits found
  Funhe2Exx11m038259t1    gnl|uv|NGB00360.1:1-58  100.00  20      0       0       446     465     58      39      1.4     40.2
  Funhe2Exx11m038259t1    gnl|uv|NGB00360.1:1-58  100.00  17      0       0       1       17      42      58         81   34.3
  Funhe2Exx11m038259t1    gnl|uv|NGB00846.1:1-65  100.00  20      0       0       446     465     65      46      1.4     40.2
  Funhe2Exx11m038259t1    gnl|uv|NGB00846.1:1-65  100.00  17      0       0       1       17      49      65         81   34.3
  # Query: Funhe2Exx11m015740t1 aalen=399,67%,complete; clen=1782; offs=248-1447; pubid=Funhe2EKm004566t1; oid=Fungr1EG3m009288t2; organism=Fundulus_heteroclitus; type=mRNA; isoform=1; Name=DNA-directed RNA polymerases I and III subunit RPAC1 (100%P); genegroup=FISH11G_G3442; Dbxref=kfish2gene:Funhe2EKm004566t1,CDD:132910,TrEMBL:UniRef50_UPI0002C64BBF,UPI0002C64BBF,family:FISH11G_G3442;
  Funhe2Exx11m015740t1    gnl|uv|NGB00362.1:1-61  97.30   37      1       0       1731    1767    61      25      3e-08   65.8
  Funhe2Exx11m015740t1    gnl|uv|NGB00362.1:1-61  100.00  17      0       0       1766    1782    57      41         81   34.3
  Funhe2Exx11m015740t1    gnl|uv|J02459.1:1-48502-49      100.00  16      0       0       1214    1229    36520   36505     315   32.4
      
=item BUG parsing vecscreen v1 201501

  201501, BUG now parsing vecscreen v1 output: NOT No hits cases lost including Strong + Medium, 
  .. mixed Strong + Weak hits? no, should be in output vectab
  .. ** bug is multiple match types not handled yet ** last only used
  
#>Funhe2Exx11m015740t1 aalen=399,67%,complete; clen=1782; offs=248-1447; pubid=Funhe2EKm004566t1; oid=Fungr1EG3m009288t2; 
#  Name=DNA-directed RNA polymerases I and III subunit RPAC1 (100%P); 
#  genegroup=FISH11G_G3442; Dbxref=kfish2gene:Funhe2EKm004566t1,CDD:132910,TrEMBL:UniRef50_UPI0002C64BBF,UPI0002C64BBF,family:FISH11G_G3442;
>Vector Funhe2Exx11m015740t1 Screen: VecScreen Database: UniVec_Core201404
Strong match
1731	1762
Weak match
1766	1782     # Weak but end bases probably should also be trimmed
Suspect origin
1763	1765
--
#>Funhe2Exx11m011238t11 aalen=444,64%,complete; clen=2057; offs=87-1421; pubid=Funhe2EKm017474t1; oid=Fungr1EG3m007309t4; 
#  Name=Solute carrier family 39 (zinc transporter), member (89%P); 
#  genegroup=FISH11G_G5206; Dbxref=kfish2gene:Funhe2EKm017474t1,CDD:217089,TrEMBL:UniRef50_Q92504,TrEMBL:S39A7_HUMAN,family:FISH11G_G5206;
>Vector Funhe2Exx11m011238t11 Screen: VecScreen Database: UniVec_Core201404
Strong match
1	24
Weak match
2017	2033
Suspect origin
2034	2057

=cut



=item vecscreen ncbi v2/c++

  vecscreen [-h] [-help] [-query input_file] [-db dbname] [-out output_file]
    [-outfmt format] [-text_output] [-version]
   Vector screening tool, version 2.0.0
   -outfmt      0 = Show alignments pairwise, default
     1 = Do not show alignments, just contaminated range offsets

  outfmt 0 needs more parsing..
..
Database: UniVec_Core (build 6.0)
           1,788 sequences; 412,552 total letters

Strong match
1       26

Query= Funhe2Exx11m049371t1 type=mRNA; aalen=128,59%,partial5-utrpoor;
clen=653; offs=1-387; oid=Fungr1EG3m028593t1;
organism=Fundulus_heteroclitus;

Length=653

> gnl|uv|NGB00362.1:1-61 Illumina Paired End PCR Primer 2.0
Length=61

 Score = 52.6 bits (26),  Expect = 3e-04
 Identities = 26/26 (100%), Gaps = 0/26 (0%)
 Strand=Plus/Plus

Query  1   TTCCTGCTGAACCGCTCTTCCGATCT  26
           ||||||||||||||||||||||||||
Sbjct  36  TTCCTGCTGAACCGCTCTTCCGATCT  61
..

Database: UniVec_Core (build 6.0)
           1,788 sequences; 412,552 total letters

Strong match
3       33
Moderate match
74      101
163     188
Suspect origin
1       2
34      73

Query= Funhe2Exx11m127091t3 type=mRNA; aalen=103,44%,complete-utrpoor;
clen=706; offs=224-535; oid=Fungr1EG3m035176t1;
organism=Fundulus_heteroclitus;

Length=706

> gnl|uv|NGB00362.1:1-61 Illumina Paired End PCR Primer 2.0
Length=61

 Score = 62.6 bits (31),  Expect = 3e-07
 Identities = 31/31 (100%), Gaps = 0/31 (0%)
 Strand=Plus/Plus

Query  3   CGGCATTCCTGCTGAACCGCTCTTCCGATCT  33
           |||||||||||||||||||||||||||||||
Sbjct  31  CGGCATTCCTGCTGAACCGCTCTTCCGATCT  61

=cut

=item FIXME: bad vecscreen parsing

    ## FIXME: bad parse here for 2+ hits : getting max range, not vec spans
    ## bug in NNN trim: cut to 0 len : n=17 cut to zero for pogonus
    ## -- these are all BIG vectrim = all of cds; need to check that.
    # pogonusvel1pk45Loc1488t3	1	356	Strong match	UniVec_Core   << ** BAD parse of vecscreen.output
    # >Vector pogonusvel1pk45Loc1488t3 Screen: VecScreen Database: UniVec_Core
    # Strong match
    # 3	36
    # 320	356
    # Suspect origin
    # 1	2
    #............
    # pogonusvel1pk45Loc18227t1	1	631	Strong match	UniVec_Core
    ## eg name=Cytochrome c oxidase subunit IV ; PogonEG0004606t2	oid=pogonusvel1pk45Loc1488t3
    ##    oldcds=1-507,vectrim=509,
    ## CDShasTooManyXs:621 #PogonEG0011277t7	oid=pogonusvel1pk45Loc1897t  name=NEL-like	cutcds=339--1,oldcds=339-959,vectrim=961,
    # CDShasTooManyXs:497 #PogonEG0011080t2	oid=pogonusvel1pk45Loc18227t1	len=0; olen=631; nnn=0/631;	name=Glutathione S-transferase, putative	cutcds=134--1,oldcds=134-631,vectrim=631
    # #er2g: PROBLEM: keep unique name but CDShasTooManyXs:354 #PogonEG0004606t2      oid=pogonusvel1pk45Loc1488t3
    # len=0; olen=356; nnn=0/356;     name=Cytochrome c oxidase subunit IV    cutcds=2--1,oldcds=2-355,vectrim=356,
    #.........................
    
=item eg vecscreen 

  >Vector sobeetlepogo1ak39loc7758t1 Screen: VecScreen Database: UniVec (build 6.0)
  Strong match
  1387    1428
  Suspect origin
  1429    1431
  
  >Vector sobeetlepogo1ak31loc100502t1 Screen: VecScreen Database: UniVec (build 6.0)
  Weak match
  9       24    : r1
  563     578   : r2
  Suspect origin
  1       8     << add to 1st range
  579     588   << add to 2nd range

  >Vector sobeetlepogo1ak39loc88670t1 Screen: VecScreen Database: UniVec (build 6.0)
  Strong match
  1       32
  461     490
  Suspect origin
  491     491
  
  updated output vectable:
  sobeetlepogo1ak39loc88670t1     1       32      Strong match    UniVec
  sobeetlepogo1ak39loc88670t1     461     490     Strong match    UniVec << MISSED origin due to b == e

=cut
  

