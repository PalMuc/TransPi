#!/usr/bin/env perl
# evgmrna2tsa.pl; was evgrna2genbanktsa.pl  

=item notes

  EvidentialGene evgmrna2tsa.pl
  process tr2aacds.pl outputs for ncbi tsa  submit
  
  -- main/alternate id table
  -- pubids version of main/alt ids
  -- vecscreen mrna.tr
  -- asmrna2ncbitsa.pl process vecscreen data, annot table?
  -- tbl2asn project.trclass ...

  parts from evigene/scripts/
   evigene2genbanktbl.pl
   asmrna2ncbitsa.pl
   bestgenes_update.pl 
   bestgenes_puban_kfish.pl ??

=item MERGED 201405 asmrna_altreclass.pl 

  $evigene/scripts/evgmrna2tsa2.pl -species Fundulus_heteroclitus -idprefix Funhe2Exx9m -novectrim \
  -class $pt.trclass -names $pt.names -log -tidy >& log.evpub9

  this step now is merged into evgmrna2tsa2.pl
  swaps poor main, good alt cases using scores of more qualities than longest orf
  e.g. complete-aa > partial-aa, utr-good > utr-bad, ..
  
  $evigene/scripts/rnaseq/asmrna_altreclass.pl -trclass $pt.trclass -noclass=70 -maxalt=20 -altrenum \
  -out -debug >& publicset/$pt.realt.log 

=item old steps

  1. get_evgtrset() : collect tr2aacds file set, getmRNA, genenames, sra_result_cvs metadata
  2. trclass2maintab() : tables of main-alt tr, pubids
  3. vecscreen(mRNA)  : tabulates vec-locs for later putseq that NNN's out vectors
  4. trprocess()/putseq() : link tr,names,pubids; make tbl2asn, annotation file set; trim seq NNN, 
  5. call tbl2asn (as desired)
  
=item fixme vecscreen

  vecscreen can be/is cutting major CDS parts now, with reason..
  should preceed most of this processing, with re-call of .aa, .cds for veccut.mrna

=item new pre-vecscreen steps

  0. get_evgtrset() : collect tr2aacds file set, getmRNA, genenames, sra_result_cvs metadata
    -- need genenames/stats for vecscreen counterbalance, ie some vecs are wrong!
    ## socatfishv1k95loc34193t1 (origaa=446,full and 100% align zfish gene; cutaa=176,utrbad;)
    ## vecscreen finds Moderate ~39 bp align w/ 2+ mismatch in cds5span
    ## w/o uvcut, full align to alpha-1,6-mannosyl-glycoprotein 2-beta-N-acetylglucosaminyltransferase [Danio rerio]
  1a. getmRNA() only   : for vecscreen()  
  1b. vecscreen(mRNA)  : cut out vectors, recall .aa,.cds for those, revise evgtrset 
     .. reprocess evgtrset, new mRNA fileset, ..
     -- add trim NNN ends here? so new mRNA set is clean of that
  ^^ Now in evigene/scripts/rnaseq/asmrna_trimvec.pl   
  
  2. get_evgtrset() : collect tr2aacds file set, getmRNA, genenames, sra_result_cvs metadata
  3. trclass2maintab() : tables of main-alt tr, pubids
  4. trprocess()/putAnnoSeq() : link tr,names,pubids; make tbl2asn, annotation file set; trim seq NNN, 
      -- ? reannotate .mrna,.cds,.aa fileset with pubids, select .ann.txt (gene names, dbxref,..)?
      -- for public use outside ncbi tsa submit
  5. call tbl2asn (as desired)

=item new step script

	#!/bin/bash
	# evgmrna2tsa2.sh
	
	#..... test input set from /bio/bio-grid/aabugs4/tsaevgc/pogonus1all3cf/ .......
	#. /bio/bio-grid/aabugs4/tsaevgc/tsaoutz/pogonus1test3
	#. okayset/                      pogonus1all3.sra_result.csv@  pogonus1all3.vector.tab@
	#. pogonus1all3.names@           pogonus1all3.trclass@
	#. ./okayset:
	#. pogonus1all3.okalt.aa.gz@               pogonus1all3.utrorf.aa@
	#. pogonus1all3.okalt.cds.gz@              pogonus1all3.utrorf.aain@
	#. pogonus1all3.okalt.tr.gz@               pogonus1all3.utrorf.cds@
	#. pogonus1all3.okay.aa.gz@                pogonus1all3.utrorf.mrna@
	#. pogonus1all3.okay.cds.gz@               pogonus1all3.utrorf.mrna.traa2cds.log@
	#. pogonus1all3.okay.tr.gz@
	#.........
	
	evigene=/bio/bio-grid/mb/evigene
	nbin=/bio/bio-grid/mb/ncbic/bin
	export vecscreen=$nbin/vecscreen
	
	if [ "X" = "X$trclass" ]; then "echo env trclass=path/to/name.trclass"; exit -1; fi
	trpath=`dirname $trclass`
	trname=`basename $trclass .trclass`
	## need chdir?# 
	cd $trpath
	
	## Step1. trim
	##.. this will reuse trpath/trimset/$trname.vecscreen.tmp or vector.tab
	$evigene/scripts/rnaseq/asmrna_trimvec.pl -class $trname.trclass -log
	
	## Step2. pubset + submitset; add -runtbl2asn for submitset/*.sqn 
	$evigene/scripts/evgmrna2tsa2.pl -class $trname.trclass -log
#...................
  
=item old script

  $evigene/scripts/rnaseq/asmrna2ncbitsa.pl -GAPSOK -idpre Thecc1ER_ \
  -cdna ../tr5parts/pub3ig.$pt.tab4g.tr.gz -vec ../tr5parts/pub3ig.trasm.tab4g.vector.tab \
  -geneinfo ../tr5parts/pub3ig.trasm.tab4g.geneinfo1.tab  -log tr4g.$pt.log \
  -out $pt/TCM01.tsa_rasm.$pt.fsa -tbl $pt/TCM01.tsa_rasm.$pt.tbl

=item tbl2asn test

  need these input templates; generate some from other configs?
    evgr_tsamethods.cmt evgr_tsadesc.cmt evgr_tsasubmit.sbt
  
  pt=tsasub1 
  cp -p ../tsasubmit/evgr_*.{cmt,sbt} $pt/   
  # .. edit cmt per needs .. 
  # Should add PRJ of data source: PRJNA73443 .. where? in evgr_tsamethods.cmt ?
  
  org='[organism=Litopenaeus vannamei] [bioproject=PRJNA12345]'
  sra='[SRA=SRR346404]' 
  
  $ncbin/tbl2asn -p $pt/ -Z $pt/$pt.discrep.log \
     -w $pt/evgr_tsamethods.cmt -Y  $pt/evgr_tsadesc.cmt -t  $pt/evgr_tsasubmit.sbt \
     -a r10k  -l paired-ends -Vtb -Mt -XE \
     -j "[moltype=mRNA] [tech=TSA] $org $sra"

=item

  update 2017.12
  minor updates for evgpipe_sra2genes.pl

  1510: NCBI requires mrna-cds partial to abut ends of mrna,unless gaps, 
    -- problem with this, of course, when mrna-*end* is +2 from partial3-cds end,
       cannot extend cds to mrna end w/o calling 2-base codon, adding to protein :((
    -- answer: trim mrna-end to cds-end (only for +2 case or +1 also?)
       
=cut

use constant VERSION => '2018.06.18';  # bugfixes during rewrite to trclass2pubset (replaces MAIN_pubsetonly)
  # '2017.12.10'; # updates for evgpipe_sra2genes.pl
  # '2015.10.08'; #  NCBI mrna-cds partial offset by 1,2 fix in putCDSloc
  # '2014.12.31'; # various updates for separate use of submitset parts, tbl2asn, trimvec,..
  # '2014.05.21'; #adding mrna_utrorf, rnaseq/asmrna_altreclass.pl
  # '2013.09.13'; # 08.10'; # 06.21'; # 05.31; # 28; 20; 07; # 05; # '04.20'; # .16'; # 03.20
  # 13.05.07: fix big vecscreen parsing bug .. redo all data from that

#** 2018.06 OBSOLETE getmRNA_utrorf(), now utrorfs in cdna_evigenesub:getmRNA()
use constant OBSOLETE_getmRNA_utrorf => 1;

use FindBin;
use lib ("$FindBin::Bin","$FindBin::Bin/prot"); # assume evigene/scripts/ layout: main, prot/, rnaseq/, ...
# ? should this script move to rnaseq/ subdir?

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname);
use cdna_evigenesub; #added # 
#maybe# use cdna_proteins;
#maybe# use protein_names;

   
## evigene path = self path
##my $EVIGENES="$FindBin::Bin"; #??
# cdna_evigenesub globals:
our $EVIGENES="$FindBin::Bin";  
our $EGAPP='mrna2tsa';  
our $EGLOG='m2t';
our $dryrun=0; ## $DRYRUN ?
our $DEBUG= $ENV{debug}|| 0;
our (%genenames, %genedbxref, %genenamepct, %namedgenes, %cddnames, %pubids, %pubidinfo); # cdna_evigenesub.pm globals

## Evigene tr2aacds subdirs: see tidyup
## add for trimvec, mrna2tsa:  trimset?  publicset? tsasubmit/submitset ?
## change trimvec,mrna2tsa output subdir: publicset? adding pubids, main2alt, ...
## separate tsasubmit subdir ..  ?? add trimset?
our @evgdirs = qw(okayset dropset inputset trimset tmpfiles erasefiles publicset submitset);
our (@okayset,@dropset,@inputset,@tmpfiles,@erasefiles,@publicset,@submitset); # tidyup file sets
#? our $outdir='publicset'; # pubset ?
my $tidyup=0;

our $okayclass  ='main|noclass';
our $okaltclass ='alt|part';
our $dropclass  ='drop|cull';
# problem class: maybeok.althi1 (trclass); some have best refaa quals, but this is large subset
#  s/maybeok/okay/;

## these move to asmrna_trimvec
my $MAXGAP=15; # NCBI 
my $ENDGAP=20; ## was 10; # trim ends if gaps w/in this of ends; NCBI changed again.
  ## FIXME, NCBI tbl2asn changed this to 20 !
  ## [SEQ_INST.HighNContentStretch] Sequence has a stretch of at least 10 Ns within the first 20 bases BIOSEQ: lcl|MusacuEGm017348t1: delta, rna len= 509
my $MINSIZE=200; # NCBI
my $GAPSOK=1; # default on? new policy 2012Dec for TSA tbl2asn: -a r10u -l paired-ends
my $UniVecDB= $ENV{UniVec} || "UniVec_Core"; # Not UniVec
#..............

my $AAMIN_NOCLASS=$ENV{aaminnoclass}||60; # asmrna_altreclass -noclasscut=$AAMIN_NOCLASS; drop noclass (needs user opt), but rescued by aaref

my $MINGENEIDENT=85; # for asm == gene identity
my $DEFAULTidpre= 'EVGmRNA';  ## this seems to fix -idprefix ignored bug
my $IDPREFIX= $ENV{idprefix} || $DEFAULTidpre; ## "evgr"; #  opt
my $GDB_PREFIX='gnl|Evigene|';  #see below; 'gnl|Evigene|'; # use IDPREFIX ? or not, since ID has this also
my $pubidnum_start=0;
my $NCPU= 1; 

my $DOtbl2asn=0; 
my($organism,$sraids,$BioProject)=("Noname","SRR000000","");  ## SRR 346404
my $TSADESC="-w evgr_tsamethods.cmt -Y evgr_tsadesc.cmt -t evgr_tsasubmit.sbt"; # FIXME: where from?
    ## ^^ need config for these; generate some of this?
my $DATE=`date '+%Y%m%d'`; chomp($DATE); # default is today; use perl func?

## namegenes messed up, didn't screen out loqualnames:
##  1%,309/22971,716        Dumpy
##  0%,126/34350,248        Titin
our $MIN_NAMEIDENT = 35;  # min similar% for protein evid align/equivalence, for naming; JCVI uses 35%  
our $MIN_IDLIKE   = 15;  # low for now; 15..20 seems right;add Note with loqualname
use constant NAME_NOEDITS => 1;

## use protein_names.pm
use constant FOR_NCBI => 1;
our $NAME_NONE = "Unknown|Uncharacterized conserved|Uncharacterized|Hypothetical";  
our $NAME_UNK  = 'Uncharacterized protein'; # uniprot  
our $NAME_UNKNCBI= 'hypothetical protein';  # NCBI prefers this to Uncharacterized .. special handling

my $SKIPTRIMSET=0;
my $SHOWDROPS=0; 
my $skipdropseqs=1; # make_pubseq fixup flag, default on/off? special cases of changing trclass 
my $FAHDRisValid= 0; # for annotab2tblinfo annot.tab vs mRNA header conflicts
my $preserveOldIds=0;

my ($vecscreenf,$trclass,$genenames,$cdnaseq,$output,$logfile,$tblfile,$runsteps); ## ,$dryrun above
my $DBRECODE= $ENV{dbrecode};

my %DEFAULT_SETTINGS= ( #? add LOCUSTAG
  IDPREFIX=>$DEFAULTidpre, LOCUSTAG => $DEFAULTidpre,
  TSADESC=>$TSADESC, DATE=>$DATE, 
  MINSIZE=>$MINSIZE, MAXGAP=>$MAXGAP,
  sraids => $sraids, organism => $organism, BioProject => $BioProject, #fixme bioproj lost
  trclass => '', mrna => '', genenames=>'', 
  ); # vecscreen => '',
  

# FIXME below: save Options to trname.info.stash : TSADESC, organism, ... 

my $MAIN_pubsetonly=0; #:  $SKIPTRIMSET=1; $skipTSAparts=1;  $DOtbl2asn= 0;
my $MAIN_submitonly=0; #:  $SKIPTRIMSET=0; $skipTSAparts=0;  $DOtbl2asn= 1;

my @saveopt= grep /^\-/, @ARGV;
my $optok= GetOptions(
  # "config=s", \$config, "cadd=s", \@configadd,
  "mrna|cdna=s", \$cdnaseq,
  "class|trclass=s", \$trclass,
  "names|genenames=s", \$genenames, ## ? allow for 2 files: myspecies.namerefids + allrefprot.names
  "organism|species=s", \$organism,   
  "sraids=s", \$sraids,   
##settings#  "TSADESC=s", \$TSADESC,   
  "output:s",  \$output,
  "tblfile:s", \$tblfile,  
  "logfile:s", \$logfile,
  "idprefix=s", \$IDPREFIX,  # FIXME: idpre option  overwritten by spppref
  "gdbprefix=s", \$GDB_PREFIX, #?? IDPREFIX default 
  "keepoldids|preserveOldIds=s", \$preserveOldIds,  # upd1806; allow -idpre=XXX -keepold w/o 2nd param
  "DBRECODE=s", \$DBRECODE,  
  "DATE=s", \$DATE,  
  "MINSIZE=i", \$MINSIZE,  
  "MAXGAP=i", \$MAXGAP,  
  "NCPU=i", \$NCPU,## "MAXMEM=i", \$MAXMEM,  
  "runsteps=s", \$runsteps,   
  "runtbl2asn!", \$DOtbl2asn, 
  "onlypubset!", \$MAIN_pubsetonly, 
  "onlysubmit!", \$MAIN_submitonly, 
  "notrimvec|novectrim", \$SKIPTRIMSET,
  # "MINGENEIDENT=i", \$MINGENEIDENT,  ## not used?
  "GAPSOK!", \$GAPSOK, # default on always?
  "dropshow!", \$SHOWDROPS,  
  "skipdropseqs|onlyokaypublicseq!", \$skipdropseqs,  # default on?
  "dryrun|n!", \$dryrun, 
  "tidyup!", \$tidyup, 
  "debug!", \$DEBUG, 
  );


# OPTION: add -outdir opt ; keep list of created files, move there, or work there..
# OPTIONed: -class evg_tr2aacds.trclass only requirement, gives paths to mrna, names.tab ?
# OPTIONed: here or caller? make mrna/cdnaseq  from evg traa2cds.pl okayset/name.{okay,okalt}.{tr,aa}  
# $evigene/scripts/prot/traa2cds.pl -trout -cdna $pt.tr.gz -aa $pt.aa.gz -out -log
##  -trout : cdna output, not cds, as name.mrna.tr
#?No# $cdnaseq= shift @ARGV unless($cdnaseq);
## caller script bug: #er2g: ERR: unused arguments: vannamei'
## from  opts="$opts -species='$species'"

die "EvidentialGene mrna2tsa VERSION ",VERSION,"
  make mRNA fasta and annotation table for TSA submit via tbl2asn
  using mRNA classes from Evigene tr2aacds, gene product names, vecscreen
Usage: evgmrna2tsa.pl -trclass mrna.trclass -mrna mrna.fasta ...
opts: -idprefix Thecc1EG -names mrna.names -log -dryrun -debug
    -runtbl2asn -organism=$organism -SRAids=$sraids 
    -NCPU=$NCPU -MINSIZE=$MINSIZE  -MAXGAP=$MAXGAP\n"
  unless($optok and ($cdnaseq or $trclass) and (-s $trclass or -s $cdnaseq));  
#  -out=out.fsa  -conf=evigene.conf  -vecscreen=infile -TSADESC='$TSADESC'

## openloggit(), loggit() now in subs.pm

$tidyup= 1 unless($dryrun); # ||$DEBUG default on unless debug|dryrun ?
# evigene_config($config, \@configadd); # always even if $config null
unless($DATE) { $DATE=`date '+%Y%m%d'`; chomp($DATE); } # fixme

 
# #*** defer IDPREFIX use till read sra_result.cvs ..
# #  get species/organism from that, make IDPREFIX from Spp.org. abbev. 
# my $ChangeDefaultIdPrefix= ($IDPREFIX eq "EVGmRNA")?1:0; 
# my $nd= ( $IDPREFIX =~ s/(0\d+)$// ) ? length($1) : 6;
# my $pubid_format = $IDPREFIX.'%0'.$nd.'d'; # $public_options{'publicid'} || "evgr000000";
# my $altid_format = 't%d'; # t%d or t%02d # $public_options{'altid'} || "t00";
# my $GBPROID= $IDPREFIX."_".$DATE; # "cacao11evigene_20120827";
# #? No?# unless($GDB_PREFIX) { $GDB_PREFIX= "gnl|$IDPREFIX|"; } #?

openloggit($logfile,$trclass||$cdnaseq); # which first? publicset/logfile or ./logfile ?
loggit(1, "EvidentialGene mrna2tsa.pl VERSION",VERSION);
loggit(1, "ERR: unused arguments:",@ARGV) if(@ARGV>0);

## also:  $MAIN_pubsetonly=1 if($SKIPTRIMSET and not $DOtbl2asn);
if($MAIN_pubsetonly){ $SKIPTRIMSET=1;  $DOtbl2asn= 0;} #$skipTSAparts=1; 
elsif($MAIN_submitonly){ $SKIPTRIMSET=0;  $DOtbl2asn= 1;} #$skipTSAparts=0; $FAHDRisValid=1 ??

my $APPtbl2asn  =  findapp("tbl2asn") if($DOtbl2asn); # add this, needs configs + data doc files (.cmt, .sbt,..)
## .. tbl2asn can be real slow .. use fasplit(ncpu); forkcmd(tbl2asn,ncpu,..)
my $APPtraa2cds= findevigeneapp("prot/traa2cds.pl");

my $APPtrimvec= findevigeneapp("rnaseq/asmrna_trimvec.pl"); # NOT run here, see get_trimvecset()
## part of asmrna_trimvec: findapp("vecscreen");  

my $APPaltreclass= findevigeneapp("rnaseq/asmrna_altreclass.pl",1); # upd 201405 ;  NODIE
if($APPaltreclass =~ /MISSING/) { $APPaltreclass=""; }

## FAIL at this point if any apps missing?
#-------------------------------------

		
#above# our (%genenames, %genedbxref, %genenamepct, %namedgenes, %cddnames, %pubids); # cdna_evigenesub.pm globals

#gone: %vecscreen, %geneinfo, 
my($trpath,$trname,$sradatah, %settings);
my($pubid_format,$altid_format,$GBPROID);
my %didpubid=(); #only for make_tblfsaset ??

## add OPTION $steps : what to do/not here ??
$runsteps||="";
# use e.g. to make only trclass2maintab
#unused# my $skipmaintabrun=($runsteps and $runsteps !~ /main|ids|pub/) ? 1 : 0;
#gone# my $skipvecrun=($runsteps and $runsteps !~ /vec/) ? 1 : 0; ## see SKIPTRIMSET

my $skiptrrun= ($runsteps and $runsteps !~ /tblout|process/) ? 1 : 0; # annot and ncbi input/output files
$DOtbl2asn= ($runsteps and $runsteps !~ /asn|ncbi|tbl2/) ? 0 : $DOtbl2asn;
$DOtbl2asn=0 if($skiptrrun);
my $skipTSAparts=($SKIPTRIMSET and not $DOtbl2asn)?1:0;
# ^^upd1712: maybe $skipTSAparts=1 unless requested via onlysubmit ?

our %DBXREF_OK=(); # ncbi/db_xref; should get from evigene_config();
{
	my $dbxref_ok="AFTOL AntWeb APHIDBASE ApiDB ASAP ATCC Axeldb BEETLEBASE BGD Bold CABRI CCAP CDD dbEST
dbProbe dbSNP dbSTS dictyBase EcoGene ENSEMBL EPD ESTLIB FANTOM_DB FBOL FLYBASE GABI GDB
GenBank GeneDB GeneID GI GO GOA Greengenes GRIN H-InvDB HGNC HMP HOMD HSSP IKMC IMGT
InterPro IntrepidBio IRD ISFinder JCM JGIDB LocusID MaizeGDB MGI MIM
miRBase MycoBank NBRC NextDB niaEST NMPDR NRESTdb OrthoMCL Osa1 Pathema PBmice PFAM PGN
Phytozome PIR PomBase PSEUDO PseudoCap RAP-DB RATMAP RFAM RGD RiceGenes RZPD SEED SGD SGN
SK-FST SoyBase SubtiList TAIR taxon TIGRFAM UNILIB UniProtKB/Swiss-Prot
Swiss-Prot SwissProt UniProtKB/TrEMBL TrEMBL UNITE VBASE2 ViPR WorfDB WormBase ZFIN";
	my @pg= split/[,\s]+/, $dbxref_ok; 
	map{ my($k,$v)=split/[=:]/,$_; $v=1 unless(defined $v); $DBXREF_OK{$k}=$v; } @pg; 
}       

our %DBXREF_RECODE=();
if($DBRECODE) {
  my @pg= split/[,\s]+/, $DBRECODE; 
  map{ my($v,$k)=split/[=:]/,$_; my @k=split/[|]/,$k; foreach my $k1 (@k) { $DBXREF_RECODE{$k1}=$v; } } @pg; 
}


=item db_xref

http://www.insdc.org/db_xref.html
http://www.insdc.org/documents/dbxref-qualifier-vocabulary

HGNC Human Gene Nomenclature Database /db_xref="HGNC:2041"
ENSEMBL /db_xref="ENSEMBL:HUMAN-Gene-ENSG00000007102"
CDD   Conserved Domain Database  /db_xref="CDD:02194
OrthoMCL Ortholog Groups of Protein Sequences	/db_xref="OrthoMCL:OG5_130679"
UniProtKB/Swiss-Prot  /db_xref="UniProtKB/Swiss-Prot:P12345"
UniProtKB/TrEMBL  /db_xref=" UniProtKB/TrEMBL:Q00177"
GeneID  Entrez Gene Database (replaces NCBI Locus Link) /db_xref="GeneID:3054987"
	
=cut

#..............................................

=item revise MAIN for separate submitset use

 2014.12.27: revise MAIN for 2 separate uses: 
   a. current all in one option (w/ or w/o submitset parts, often w/o submitset part)
   b. submitset only, from prior publicset (which is more usual need)
      -- use fewer assumptions of file sets, check more
      -- BUT want update of publicset after trimvec..

  * trimvecset merge w/ pubset seq seems only problem to separate these.
     $name.trimvec_done is flag file this has been done on publicset seqs.
     asmrna_trimvec.pl has code for update_mrna_fileset(pubset,trimfiles) .. reinstate that
     add a./b. option: pubset only, submit only, or both.

  * this works now, creates new publicset/name.{mrna,aa,cds} set
    $evigene/scripts/rnaseq/asmrna_trimvec.pl -nodeferupdate  -mrna publicset/kfish2evgxx11.mrna -log
     
  submitset parts:
    get_srainfo()
    get_trimvecset() : assumes prior run of asmrna_trimvec .. or do that here?
    update_mrna_fileset(pubfiles,@trimfiles) ?? need this to merge trim/uvcut mrna/.. with prior publicset
     #.. update pubset per below? maintab,altreclass,make annot,pubseq ??
    make_tblfsaset()  : NOT skipTSAparts
    tsa_tbl2asn()  : DOtbl2asn
      
  publicset, not submitset parts:
    make_IDPREFIX()  
    update_mrna_fileset(..) 
    trclass2maintab() .. convert trclass to mainalt.tab
    APPaltreclass .. shuffles alts/locus for better choices
    make_annotab(): ann.txt
    make_pubseq(): mrna, aa, cds

=cut

#..............................................

if($MAIN_pubsetonly) { MAIN_pubsetonly(); }
elsif($MAIN_submitonly) { MAIN_submitonly(); }
else { MAIN_sub(); }

# FIXME: replace this MAIN_sub w/ above two: pubset + submit, depending on flags
sub MAIN_sub {  
 
  loggit(0, "BEGIN mrna2tsa with input=",$trclass||$cdnaseq,"date=",`date`);

	## FIXMEd: write/restore collected info to $trname.info.stash or such, so don't have to refind: SRR ids, IDPREFIX, .. @saveopt
	do_settings("restore",$trclass||$cdnaseq,); # ("log|restore|save");

	($cdnaseq,$trpath,$trname)= get_evgtrset($trclass,$cdnaseq,"publicset"); # subs.pm; cdnaseq => mrnaseq here
		loggit(0, "get_evgtrset=",$cdnaseq,$trpath,$trname); ## facount($cdnaseq)
		loggit(LOG_DIE, "Missing -mrna",$cdnaseq) unless($cdnaseq and -s $cdnaseq);
	($sradatah)= get_srainfo($trpath,$trname);
	# ^^ makes files for submitset/
	
	## FIXMEd below: trclass2maintab $pubids may exist ; use that for IDPREFIX if so ..
	($pubid_format,$altid_format,$GBPROID)= make_IDPREFIX(); # default abbrev of $organism now
	
## BUGGERS: bad pubid for utrorfs; some/most? are not in okayset/*.tr but in okayset/*.aa,cds only; not in .mrna 
## how many can we throw out? = have althi that are not utrorf.
# oid=sobeetlepogo1ak21loc83t3utrorf ; sobeetlepogo1ak21loc20766t1utrorf ; sobeetlepogo1ak21loc34353t1utrorf
#  pogonusvel1pk45Loc2448t2utrorf pogonusvel1pk45Loc124t7utrorf
## only .aa,.cds these are in utr-parts of .mrna
# > type=protein; Name=Ras-related protein Rab-20; Dbxref=CDD:31297,UniRef50_Q9NX57,UniProt:RAB20_HUMAN; aalen=; clen=0; offs=0; oid=sobeetlepogo1ak21loc83t3utrorf; evgclass=maina2,okay,match:sobeetlepogo1lk21loc18244t2,pct:99/100;
# grep -c ^> type= publicset/pogonus1all3.aa_pub.fa   n=1795 == num utrorf in aa.untrim,okayset/aa,.trclass
# grep -c utrorf.okay pogonus1all3.trclass  n=1795
#  grep utrorf.okay pogonus1all3.trclass | grep okay.altmid | wc = 1365 << may be tandem dups
#  grep utrorf.okay pogonus1all3.trclass | egrep -v 'okay.alt|okay.main' | wc = 11 << no alts, few

## FIXME: utrorf : made okayset/*.utrorf.{mrna,aa,cds} ; merge into update_mrna_fileset() or getmRNA/okayset ??
	
	#................. make_publicset steps ...................
  		## update_mrna_fileset should be merged w/ trclass2maintab, make_annotab()
	my $trimids= undef;
	my @trimfiles= get_trimvecset($cdnaseq, $SKIPTRIMSET);
  my($upstatus,$upfiles,$uptemp,$upokids); # my @upfiles=(); $upokids == \%okids
	$upstatus=0; 
	
	my $trimflag=($SKIPTRIMSET) ? 'trimvec_SKIP' : 'trimvec_done';# novectrim
  ($upstatus, $upfiles, $uptemp, $upokids)  
    = update_mrna_fileset($trpath, $cdnaseq, $trimflag, $trimids, @trimfiles);#? unless($SKIPTRIMSET);   
	push @publicset, @$upfiles if(ref $upfiles); #??
	push @tmpfiles, @$uptemp if(ref $uptemp); # these are .untrim.{mrna,aa,cds} NOT @trimfiles; need for redo/errcheck
	# exit() if($runsteps =~ /trimvec/); #DEBUG
	
			## trclass2maintab should be done AFTER merge/drops of trimvecset
			## but nicer to have pubids first, then update fileset including add >pubid annot headers
## FIXME: utrorf not in upokids ; only mrna ids ..  
			
	my($maintab,$pubids,$nmaintr,$nalltr)= trclass2maintab($trclass,"publicset",$upokids);
	loggit(0, "trclass2maintab primary n=",$nmaintr,"allntr=",$nalltr,$pubids); 
	$upstatus++ if($nmaintr>0);
	push @publicset, $maintab, $pubids;

  ## FIX 201405: insert here asmrna_altreclass.pl, .realt_pubids to become new .pubids ..

  altreclass_block($trclass,$pubids);
  
  # 	if(not $ENV{norealt} and $APPaltreclass and -x $APPaltreclass) {
  # 	  # ADD OPTIONS: -noclass=$MINAA == drop noclass short things if no other good qualities
  # 	  #  -MAXALTSAME=n == drop excessive number of althi class that appear same by aa size/aa cluster, 
  # 	  #  asmrna_altreclass to be run or not?   other opts?
  # 	  
  # 	  my $realtids="$pubids.realt";
  # 	  ## needs also $path/inputset/$trname.aa.qual
  # 	  my $cmd="$APPaltreclass -trclass $trclass -pubids $pubids -altrenum -debug -noclasscut=$AAMIN_NOCLASS -out $realtids > $realtids.log 2>&1";
  #     my $runerr= runcmd($cmd);
  #     if(-s $realtids) {
  #       rename($pubids,"$pubids.old"); rename($realtids,$pubids); 
  #       push @publicset, "$pubids.old", "$realtids.log"; # or tmpfiles ?
  #     } else {
  #       loggit(LOG_WARN,"ERR: failed update $realtids by $APPaltreclass");  # 2014 add notice...
  #       push @publicset, $realtids, "$realtids.log"; # or tmpfiles ?
  #     }
  #     my $realtstat=`tail -n1 $realtids.log`; chomp($realtstat); $realtstat ||= "$APPaltreclass updated $pubids";
  #   	loggit(0, $realtstat, "log=$realtids.log"); # log tail has stats 
  # 
  # 	} else {
  #     loggit(LOG_WARN,"PLEASE ALSO RUN publicset fixup:\n",  # 2014 add notice...
  #       "evigene/scripts/rnaseq/asmrna_altreclass.pl -trclass $trclass -altrenum -out -debug");
  #   }	
	
	
	read_pubids($pubids, $cdnaseq); #  if($pubids);

  #FOXME: annotab had WRONG aalen/qual from mrna,  not updated for uvcut.aa or utrorf.aa
  # BUT okayset/okay,okalt.aa has wrong (rev) cdsoffs versus .mrna from makemRNA()/APPtraa2cds
  # fixed in make_annotab();  parse pubset/*.aa  as well as .mrna
  # .. maybe should start annotab earlier, with aaqual/cdsoffs updated as needed?
  
	#FIXME: also re-write publicset/mrna,aa,cds w/ >pubid  old=oldid; name=xxx, other annots in header ...
	my($annotab, $ntra1, $ncdsa1) 
		= make_annotab($cdnaseq,$trclass,$skiptrrun); # add main/alt pub ids, other geneinfo 
	push @publicset, $annotab;
  
	## FIXME: read annot.tab into %annot for following subs ......
  $FAHDRisValid=0; #??
  my %annothash= read_annotab($annotab);
  
	my($pubmrna,$npm)	= make_pubseq($cdnaseq,'mRNA',\%annothash);
	my($pubaa,$npa) 	= make_pubseq(makename($cdnaseq,'.aa'),'protein',\%annothash);
	my($pubcds,$npc)	= make_pubseq(makename($cdnaseq,'.cds'),'CDS',\%annothash);
	#? push @publicset, $pubmrna,$pubaa,$pubcds; ## do in make_pubseq
	#? push @tmpfiles, $cdnaseq,$aaseq,$cdsseq; ## do ditto
	loggit(0,"publicset: ",$pubmrna,$pubaa,$pubcds,$annotab); 
    
	#?? make_publicset: merge make_annotab, update_mrna_fileset(), trclass2maintab ??
	#................. make_publicset end ...................
	
# 	## moved to asmrna_trimvec.pl
# 	($vecscreenf)= vecscreen($cdnaseq,$vecscreenf||"",$skipvecrun);

	## FIXMEd: write/restore collected info to $trname.info.stash or such, so don't have to refind: SRR ids, IDPREFIX, .. @saveopt
	do_settings("log|save",$trclass||$cdnaseq,); # or after last call?  ("log|restore|save");

# FIXME! : rework trprocess()/tsa_tbl2asn and public output files
#   .. a. primary public output = .mrna,.cds,.aa,.anno.txt with pubids, useful annot >pubid name,dbxref,size,qual
#   .. b. tsa submit files are simple transform of primary publicset/  
#       ... submit.fsa == public.mrna stripped of header annots
#       ... submit.tbl == modest transform of .anno.txt to tbl format
#   .. then as needed recreate .tbl,.fsa and run tbl2asn w/ small perl script 
#      so that one-off edits, revisions by hand can be done w/ less hassle.
#       (drop script to submitset/runtbl2asn.pl? or submitset_remake.pl)
	
# REVISE trprocess to NOT do vecscreen or gaptrim .. now in asmrna_trimvec
# REVISE trprocess/putseq into two steps: make_publicset() and make_tblfsaset()
	# my($outfa, $tblout, $annotab, $ntrout, $ncdsout) 
	#	   = trprocess($cdnaseq,$trclass,$skiptrrun); # add main/alt pub ids, other geneinfo 

  my($outfa, $tblout, $ntrout, $ncdsout, $pepout)=("") x 9;
  unless($skipTSAparts) {
	($outfa, $tblout, $ntrout, $ncdsout, $pepout) 
			= make_tblfsaset($cdnaseq,\%annothash,$skiptrrun);  
	push @submitset, $outfa, $tblout, $pepout;  # <<added pep, may be missing

	$upstatus++ if($ntrout>0); # 0 if files already made
	loggit(0,"submitset: ntr=$ntrout, ncds=$ncdsout in $outfa, $tblout"); 
	}
	
	if($DOtbl2asn) {
		my($npartgot,$spldir,$sqnoutlist)= tsa_tbl2asn($outfa,$tblout,$organism,$sraids);
		loggit(0,"tsa_tbl2asn nparts=$npartgot, submitset=$sqnoutlist"); #?? $spldir  
		## ?? push @submitset, $sqnoutlist; # leave in $spldir !!
	}

  if( $tidyup and $upstatus > 0 ) {  
    tidyupFileset("publicset",@publicset);  
  	tidyupFileset("submitset",@submitset);  #? after or before tsa_tbl2asn? need path in fileset
    tidyupFileset("tmpfiles",@tmpfiles);  
    my @rmlist;
    foreach my $fn (@erasefiles) { if(-f $fn) { unlink($fn); push @rmlist,$fn; } }  
    if(@rmlist) { my $nrm=@rmlist; my $rml=join" ",@rmlist[0..4]; loggit(0,"tidyup erase: n=$nrm, $rml .."); } 
  }

  loggit(0, "DONE at date=",`date`);
  loggit(0, "======================================"); # log may append runs
} # end MAIN


## BUGGERS: bad pubid for utrorfs; some/most? are not in okayset/*.tr but in okayset/*.aa,cds only; not in .mrna 
## how many can we throw out? = have althi that are not utrorf.
# oid=sobeetlepogo1ak21loc83t3utrorf ; sobeetlepogo1ak21loc20766t1utrorf ; sobeetlepogo1ak21loc34353t1utrorf
#  pogonusvel1pk45Loc2448t2utrorf pogonusvel1pk45Loc124t7utrorf
## only .aa,.cds these are in utr-parts of .mrna
# > type=protein; Name=Ras-related protein Rab-20; Dbxref=CDD:31297,UniRef50_Q9NX57,UniProt:RAB20_HUMAN; aalen=; clen=0; offs=0; oid=sobeetlepogo1ak21loc83t3utrorf; evgclass=maina2,okay,match:sobeetlepogo1lk21loc18244t2,pct:99/100;
# grep -c ^> type= publicset/pogonus1all3.aa_pub.fa   n=1795 == num utrorf in aa.untrim,okayset/aa,.trclass
# grep -c utrorf.okay pogonus1all3.trclass  n=1795
#  grep utrorf.okay pogonus1all3.trclass | grep okay.altmid | wc = 1365 << may be tandem dups
#  grep utrorf.okay pogonus1all3.trclass | egrep -v 'okay.alt|okay.main' | wc = 11 << no alts, few

## FIXME: utrorf : made okayset/*.utrorf.{mrna,aa,cds} ; merge into update_mrna_fileset() or getmRNA/okayset ??
# FIXME! : rework trprocess()/tsa_tbl2asn and public output files
#   .. a. primary public output = .mrna,.cds,.aa,.anno.txt with pubids, useful annot >pubid name,dbxref,size,qual
#   .. b. tsa submit files are simple transform of primary publicset/  
#       ... submit.fsa == public.mrna stripped of header annots
#       ... submit.tbl == modest transform of .anno.txt to tbl format
#   .. then as needed recreate .tbl,.fsa and run tbl2asn w/ small perl script 
#      so that one-off edits, revisions by hand can be done w/ less hassle.
#       (drop script to submitset/runtbl2asn.pl? or submitset_remake.pl)
# REVISE trprocess to NOT do vecscreen or gaptrim .. now in asmrna_trimvec
# REVISE trprocess/putseq into two steps: make_publicset() and make_tblfsaset()
# my($outfa, $tblout, $annotab, $ntrout, $ncdsout) 
#	   = trprocess($cdnaseq,$trclass,$skiptrrun); # add main/alt pub ids, other geneinfo 
	
sub MAIN_submitonly 
{
  # MAIN_submitonly:  $SKIPTRIMSET=0; $skipTSAparts=0;  $DOtbl2asn= 1;
  loggit(0, "BEGIN submitonly with input=",$trclass||$cdnaseq,"date=",`date`);

	## FIXMEd: write/restore collected info to $trname.info.stash or such, so don't have to refind: SRR ids, IDPREFIX, .. @saveopt
	do_settings("restore",$trclass||$cdnaseq,); # ("log|restore|save");

	($cdnaseq,$trpath,$trname)= get_evgtrset($trclass,$cdnaseq,"publicset"); # subs.pm; cdnaseq => mrnaseq here
		loggit(0, "get_evgtrset=",$cdnaseq,$trpath,$trname); ## facount($cdnaseq)
		loggit(LOG_DIE, "Missing -mrna",$cdnaseq) unless($cdnaseq and -s $cdnaseq);

  ($sradatah)= get_srainfo($trpath,$trname); #  makes files for submitset/
  # .. how is srainfo used now ??

  my($upstatus,$upfiles,$uptemp,$upokids); # my @upfiles=(); $upokids == \%okids
  $upstatus=0; 
  
  #..... trimvecset parts : assume this is already run via asmrna_trimvec.pl ?? or check??
  my $flagfile= makename($cdnaseq,'.trimvec_done'); 
  if( -f $flagfile) {
  
  } else {
    my $trimids= undef;
    my @trimfiles= get_trimvecset($cdnaseq, 0); # $SKIPTRIMSET
    my $trimflag='trimvec_done';# ($SKIPTRIMSET) ? 'trimvec_SKIP' : 
    ($upstatus, $upfiles, $uptemp, $upokids)  
      = update_mrna_fileset($trpath, $cdnaseq, $trimflag, $trimids, @trimfiles);#? unless($SKIPTRIMSET);   
    push @publicset, @$upfiles if(ref $upfiles); #??
    push @tmpfiles, @$uptemp if(ref $uptemp); # these are .untrim.{mrna,aa,cds} NOT @trimfiles; need for redo/errcheck
  }
  
  # FIXME: Many of trimvec/uvcut updates are wrong CDS annots (cdsoff bad) for tbl2asn ......
  # .. update uvcut.mrna has wrong cdsoff vs uvcut.aa ; problem in asmrna_trimvec.pl **
  
  #....... end trimvecset .........

  # FIXME: need %pubids for tblfsaset() ; see pubset sub..
  $FAHDRisValid=1; #??
  use constant ANNOT_SAVEPUBIDS => 0; # use read_pubids instead

  ## PROBLEM here where ann.txt out-of-date wrong annots, add checks w/ cdnaseq hdr info, at least oid check
  ## reuse %upokids as %pubids  ?? or read thru cdnaseq for pubids, oids ??
  ## DONT assume annotab is accurate for all, esp trimvec set. assume Name/Dbxref ok, but not cdsoff, trlen, ..

  read_pubids( makename($cdnaseq,".pubids"), $cdnaseq); # CHECK for file, or read from annotab ??

  my($nnamed,$namin)= parse_genenames($genenames,NAME_NOEDITS); # here or in read_annotab? other call in make_annotab
  loggit(0, "names $genenames n=$nnamed\n"); 

  my %annothash= read_annotab( makename($cdnaseq,".ann.txt"), ANNOT_SAVEPUBIDS);
  
  my($outfa, $tblout, $ntrout, $ncdsout, $pepout)=("") x 9;
	($outfa, $tblout, $ntrout, $ncdsout, $pepout) 
			= make_tblfsaset($cdnaseq,\%annothash,$skiptrrun);  

	push @submitset, $outfa, $tblout, $pepout;  # <<added pep, may be missing
	$upstatus++ if($ntrout>0); # 0 if files already made
	loggit(0,"submitset: ntr=$ntrout, ncds=$ncdsout in $outfa, $tblout"); 
	
  my($npartgot,$spldir,$sqnoutlist)= tsa_tbl2asn($outfa,$tblout,$organism,$sraids);
  loggit(0,"tsa_tbl2asn nparts=$npartgot, submitset=$sqnoutlist"); #?? $spldir  
	## FIXME: submitset tbl2asn: name.val,.sqn,.discrep,.fixedproducts errorsummary.val

	do_settings("log|save",$trclass||$cdnaseq,); # or after last call?  ("log|restore|save");
  if( $tidyup and $upstatus > 0 ) {  
    tidyupFileset("publicset",@publicset);  
  	tidyupFileset("submitset",@submitset);  #? after or before tsa_tbl2asn? need path in fileset
    tidyupFileset("tmpfiles",@tmpfiles);  
    my @rmlist;
    foreach my $fn (@erasefiles) { if(-f $fn) { unlink($fn); push @rmlist,$fn; } }  
    if(@rmlist) { my $nrm=@rmlist; my $rml=join" ",@rmlist[0..4]; loggit(0,"tidyup erase: n=$nrm, $rml .."); } 
  }
  loggit(0, "DONE at date=",`date`);
  loggit(0, "======================================"); # log may append runs
} # MAIN_submitonly


sub altreclass_block {
  my($trclass,$pubids)= @_;
  
  ## FIX 201405: insert here asmrna_altreclass.pl, .realt_pubids to become new .pubids ..
  
	if(not $ENV{norealt} and $APPaltreclass and -x $APPaltreclass) {
	  # ADD OPTIONS: -noclass=$MINAA == drop noclass short things if no other good qualities
	  #  -MAXALTSAME=n == drop excessive number of althi class that appear same by aa size/aa cluster, 
	  #  asmrna_altreclass to be run or not?   other opts?
	  
	  my $realtids="$pubids.realt";
	  my $aopts="-debug -noclasscut=$AAMIN_NOCLASS";
	  
	  (my $mapqual=$trclass)=~ s/\.\w+$/.align.tab/;
	  if(! -f $mapqual and -d "geneval") { $mapqual="geneval/$mapqual"; }
	  if(-s $mapqual) { $aopts.=" -mapqual $mapqual"; }
	  
	  my $cmd="$APPaltreclass $aopts -altrenum -trclass $trclass -pubids $pubids -out $realtids > $realtids.log 2>&1";
    my $runerr= runcmd($cmd);
    
    if(-s $realtids) {
      rename($pubids,"$pubids.old"); rename($realtids,$pubids); 
      push @publicset, "$pubids.old", "$realtids.log"; # or tmpfiles ?
    } else {
      loggit(LOG_WARN,"ERR: failed update $realtids by $APPaltreclass");  # 2014 add notice...
      push @publicset, $realtids, "$realtids.log"; # or tmpfiles ?
    }
    my $realtstat=`tail -n1 $realtids.log`; chomp($realtstat); $realtstat ||= "$APPaltreclass updated $pubids";
  	loggit(0, $realtstat, "log=$realtids.log"); # log tail has stats 

	} else {
    loggit(LOG_WARN,"PLEASE ALSO RUN publicset fixup:\n",  # 2014 add notice...
      "evigene/scripts/rnaseq/asmrna_altreclass.pl -trclass $trclass -altrenum -out -debug");
  }	

}


=item read_pubids

  %pubids is global in cdna_evigenesub
  -- move read_pubids to cdna_evigenesub ?
  -- extend to read from mrna/aa publicset seq files instead, if no .pubids?
  
  realtids has extended fields, new drop/keep classing .. dmagset56i.realt_pubids
  now publicset/name.pubids looks this way, preserve cols 1..4 as prior version; 5=Class is also useful

  #Public_mRNA_ID           originalID                      PublicGeneID            AltNum  Class   AAqual  pIdAln  Notes
  # Dapma6iEVm000001t1      dmag4vel4ifik65Loc753t140       Dapma6iEVm000001        1       main    18640,93%,complete,44X  99/99/. .
  # Dapma6iEVm000001t2      hsICO2931no9gvelvk53Loc431t59   Dapma6iEVm000001        2       althi   17246,98%,complete      99/67/./dmag4vel4ifik65Loc753t141       oldid:Dapma6iEVm000001t10,
  # Dapma6iEVm000001t28     dmag4vel4ipak35Loc6068t3        Dapma6iEVm000001        28      althi   51,73%,complete 99/83/./dmag4vel4ifik65Loc753t141       oldid:Dapma6iEVm000001t8,
  # Dapma6iEVm000001t29     dmag4vel4ipak65Loc1003t146      Dapma6iEVm000001        29      dropalthi1      13944,65%,partial3      99/99/./dmag4vel4ifik65Loc753t141       oldid:Dapma6iEVm000001t18,

=cut

sub read_pubids 
{
  my($pubidtab,$pubseq)= @_;
  # -- extend to read from mrna/aa publicset seq files instead, if no .pubids?
  #  read_pubids( makename($cdnaseq,".pubids"), $cdnaseq); # CHECK for file, or read from annotab ??
  $pubseq||="";
  
  my($ok,$hin,$isseq,$nred,$isrealt,$realtdrops,$nalltr)=(0) x 9;
  if($pubseq and not -s $pubidtab and -s $pubseq) {
    ($ok,$hin)= openRead($pubseq); $isseq=$ok;
  } else {
    ($ok,$hin)= openRead($pubidtab); 
  }
	loggit(1,"ERR read_pubids fail from $pubidtab,$pubseq") unless($ok); 
  
  %pubids=(); # clear global .. or read 1st to local then decide..
  # my %newpubids=();
  
  while(<$hin>) { 
    unless($nred or $isseq or /^\w/) { $isseq=1 if(/^>/); }
    if($isseq) {
      if(/^>(\S+)/) { 
        my $pubid=$1; my($oid)= (m/\boid=([^;,\s]+)/)?$1:"noid";
        my($gid,$alti)= $pubid=~m/(\w+t)(\d+)$/; unless($alti) { $gid||=$pubid; $alti||=0; }
        my($aqual)=(m/aalen=([^;,\s]+)/)?$1:0;
        my($reclass)=($alti>1)?"alt":"main"; # look for evgclass=??
        
        #?? split this global %pubids into 2 so no oid/pubid confusion
        $pubids{$oid}= $pubid; # PROBLEMS, need only pubids for user.
        $pubids{$pubid}= $pubid; # PROBLEMS, need only pubids for user.
        $pubidinfo{$pubid}=join("\t", $oid,$gid,$alti,$reclass,$aqual); #  
        
        if(++$nred == 1 and $nalltr == 0) { # reset to existing IDPREFIX
          ($pubid_format,$altid_format,$GBPROID)= make_IDPREFIX($pubid); 
        }
      }
      next;
    }
    unless(/^\w/) { 
      if(not $nred and /^#Public_mRNA_ID/) { $isrealt=1 if(/AltNum\tClass/); }
      next;
    } 
    chomp;  my @v=split"\t"; next unless(@v > 1);
    my($pubid,$oid,$gid,$alti,$reclass,$aqual,$pialn,$notes)=@v;  # extended for realt_pubids
    ## ^^ save extended data: reclass, aqual?
    
    if($reclass =~ /^(drop|cull)/) { # $isrealt and .. dont need isrealt set?
      $realtdrops++; # should these go into new publicset/name.culled.files ?? keep from primary publicset files
      next;
    } 
    $pubids{$oid}= $pubid; 
    $pubids{$pubid}= $pubid; 
    ## buggers, cant put both into pubid key. only oids and pubid=pubid ..
    $pubidinfo{$pubid}=join("\t", splice(@v,1,5)); # $oid,$gid,$alti,$reclass,$aqual 
    # "$oid\t$gid\t$alti"; ## does this need gid,alti extra?
    
    if(++$nred == 1 and $nalltr == 0) { # reset to existing IDPREFIX
      ($pubid_format,$altid_format,$GBPROID)= make_IDPREFIX($pubid); 
    }
  } close($hin);

  # %pubids=(); # clear global .. or read 1st to local then decide..
  # map{ $pubids{$_}= $newpubids{$_} } keys %newpubids;
  
	loggit(0,"read_pubids n=$nred form $pubid_format from $pubidtab,$pubseq"); 
  return($nred);
}

sub MAIN_pubsetonly 
{
  # MAIN_pubsetonly:  $SKIPTRIMSET=1; $skipTSAparts=1;  $DOtbl2asn= 0;
  
  loggit(0, "BEGIN pubsetonly with input=",$trclass||$cdnaseq,"date=",`date`);

	## FIXMEd: write/restore collected info to $trname.info.stash or such, so don't have to refind: SRR ids, IDPREFIX, .. @saveopt
	do_settings("restore",$trclass||$cdnaseq,); # ("log|restore|save");

	($cdnaseq,$trpath,$trname)= get_evgtrset($trclass,$cdnaseq,"publicset"); # subs.pm; cdnaseq => mrnaseq here
		loggit(0, "get_evgtrset=",$cdnaseq,$trpath,$trname); ## facount($cdnaseq)
		loggit(LOG_DIE, "Missing -mrna",$cdnaseq) unless($cdnaseq and -s $cdnaseq);

	## FIXMEd below: trclass2maintab $pubids may exist ; use that for IDPREFIX if so ..
	($pubid_format,$altid_format,$GBPROID)= make_IDPREFIX(); # default abbrev of $organism now
	
	#................. make_publicset steps ...................
  ## update_mrna_fileset should be merged w/ trclass2maintab, make_annotab()
  
  ## ? check for trimset/ and/or name.uvcut.mrna .. use if found
	my $trimids= undef;
	my @trimfiles= (); # get_trimvecset($cdnaseq, $SKIPTRIMSET);
	my $trimflag= ($SKIPTRIMSET) ? 'trimvec_SKIP' : 'trimvec_done';# novectrim

  my($upstatus,$upfiles,$uptemp,$upokids);  
	$upstatus=0; 	
  ($upstatus, $upfiles, $uptemp, $upokids)  
    = update_mrna_fileset($trpath, $cdnaseq, $trimflag, $trimids, @trimfiles);
    
	push @publicset, @$upfiles if(ref $upfiles);  
	push @tmpfiles, @$uptemp if(ref $uptemp); # these are .untrim.{mrna,aa,cds} NOT @trimfiles; need for redo/errcheck
	
			## trclass2maintab should be done AFTER merge/drops of trimvecset
			## but nicer to have pubids first, then update fileset including add >pubid annot headers
      ## FIXME: utrorf not in upokids ; only mrna ids ..  
			
	my($maintab,$pubids,$nmaintr,$nalltr)= trclass2maintab($trclass,"publicset",$upokids);
	loggit(0, "trclass2maintab primary n=",$nmaintr,"allntr=",$nalltr,$pubids); 
	$upstatus++ if($nmaintr>0);
	push @publicset, $maintab, $pubids;

  ## FIX 201405: insert here asmrna_altreclass.pl, .realt_pubids to become new .pubids ..

  altreclass_block($trclass,$pubids);
  
  # 	if($APPaltreclass and -x $APPaltreclass and not $ENV{norealt} ) { 
  # 	  # ADD OPTIONS: -noclass=$MINAA == drop noclass short things if no other good qualities
  #     # skip option? $ENV{norealt} 
  # 	  #  -MAXALTSAME=n == drop excessive number of althi class that appear same by aa size/aa cluster, 
  # 	  #  asmrna_altreclass to be run or not?   other opts?
  # 	  
  # 	  my $realtids="$pubids.realt";
  # 	  ## needs also $path/inputset/$trname.aa.qual
  # 	  my $cmd="$APPaltreclass -trclass $trclass -pubids $pubids -altrenum -debug -noclass=60 -out $realtids > $realtids.log 2>&1";
  #     my $runerr= runcmd($cmd);
  #     if(-s $realtids) {
  #       rename($pubids,"$pubids.old"); rename($realtids,$pubids); 
  #       push @publicset, "$pubids.old", "$realtids.log"; # or tmpfiles ?
  #     } else {
  #       loggit(LOG_WARN,"ERR: failed update $realtids by $APPaltreclass");  # 2014 add notice...
  #       push @publicset, $realtids, "$realtids.log"; # or tmpfiles ?
  #     }
  #     my $realtstat=`tail -n1 $realtids.log`; chomp($realtstat); $realtstat ||= "$APPaltreclass updated $pubids";
  #   	loggit(0, $realtstat, "log=$realtids.log"); # log tail has stats 
  # 	} else {
  #     loggit(LOG_WARN,"PLEASE ALSO RUN publicset fixup:\n",  # 2014 add notice...
  #       "evigene/scripts/rnaseq/asmrna_altreclass.pl -trclass $trclass -altrenum -out -debug");
  #   }	

	read_pubids($pubids) if($pubids); # NOT cdnaseq here it isnt publicset/pubid file

  #FOXME: annotab had WRONG aalen/qual from mrna,  not updated for uvcut.aa or utrorf.aa
  # BUT okayset/okay,okalt.aa has wrong (rev) cdsoffs versus .mrna from makemRNA()/APPtraa2cds
  # fixed in make_annotab();  parse pubset/*.aa  as well as .mrna
  # .. maybe should start annotab earlier, with aaqual/cdsoffs updated as needed?
	#FIXME: also re-write publicset/mrna,aa,cds w/ >pubid  old=oldid; name=xxx, other annots in header ...

	my($annotab, $ntra1, $ncdsa1) 
		= make_annotab($cdnaseq,$trclass,$skiptrrun); # add main/alt pub ids, other geneinfo 
	push @publicset, $annotab;
  $FAHDRisValid=0; #??

	## FIXME: read annot.tab into %annot for following subs ......
  #  my %annothash= read_annotab(makename($cdnaseq,".ann.txt"));
  my %annothash= read_annotab($annotab);
  
	my($pubmrna,$npm)	= make_pubseq($cdnaseq,'mRNA',\%annothash);
	my($pubaa,$npa) 	= make_pubseq(makename($cdnaseq,'.aa'),'protein',\%annothash);
	my($pubcds,$npc)	= make_pubseq(makename($cdnaseq,'.cds'),'CDS',\%annothash);
	loggit(0,"publicset: ",$pubmrna,$pubaa,$pubcds,$annotab); 
    
	#  make_publicset: merge make_annotab, update_mrna_fileset(), trclass2maintab ??
	#................. make_publicset end ...................
	
	## FIXMEd: write/restore collected info to $trname.info.stash or such, so don't have to refind: SRR ids, IDPREFIX, .. @saveopt
	do_settings("log|save",$trclass||$cdnaseq,); # or after last call?  ("log|restore|save");

  if( $tidyup and $upstatus > 0 ) {  
    tidyupFileset("publicset",@publicset);  
  	tidyupFileset("submitset",@submitset);  #? after or before tsa_tbl2asn? need path in fileset
    tidyupFileset("tmpfiles",@tmpfiles);  
    my @rmlist;
    foreach my $fn (@erasefiles) { if(-f $fn) { unlink($fn); push @rmlist,$fn; } }  
    if(@rmlist) { my $nrm=@rmlist; my $rml=join" ",@rmlist[0..4]; loggit(0,"tidyup erase: n=$nrm, $rml .."); } 
  }

  loggit(0, "DONE at date=",`date`);
  loggit(0, "======================================"); # log may append runs
} # end MAIN_pubsetonly

#---------------------------------------------------  

sub do_settings {
  my($action,$pathname)= @_;
  ## write these to work dir; reread next go
  ## action == 'log|save|restore' ; restore called AFTER read new options, shouldnt replace
  my $PRES='m2t';
  my $trpname= makename($pathname,".mrna2tsa.info"); 
  $sraids=~s/ +;/;/g; $sraids=~s/ +/;/g; # ??

	## merge this and get_srainfo() from sra.cvs .. allow fields in both?
	## $sradatah->{'assemblers'} == "Velvet/Oases v; SOAPDenovoTrans v; Trinity v;"
  
  ## add for .sbt : m2t.BioProject=PRJNA201485
  ## do we need global vars for each of these?? 
  ## current; global is %settings, and DEFAULT_SETTINGS before GetOptions
## add locustagprefix and make new lotag for each mrna  
#     SubmissionID:           SUB763256
#     BioProject ID:          PRJNA269174
#     Organism name/label:    Fundulus heteroclitus/MDIBL+WHOI
#     Title:                  Fundulus heteroclitus Transcriptome or Gene expression
#m2t.LOCUSTAG =     Locus tag prefix:       RC70 (SAMN02116571)    <<dgg is this ok? conflict w/ other kfish EST project?
  
  my %mysettings=();
	map{ $mysettings{$_}=$DEFAULT_SETTINGS{$_} } (keys %DEFAULT_SETTINGS); # ensure these are in??
	# but replace default w/ local vals
  %mysettings=( ## use %settings more instead of individual vars ..
		IDPREFIX=>$IDPREFIX, TSADESC=>$TSADESC, DATE=>$DATE, 
		MINSIZE=>$MINSIZE, MAXGAP=>$MAXGAP,
		sraids => $sraids, organism => $organism, BioProject => $BioProject,  
		trclass => $trclass, mrna => $cdnaseq, genenames=>$genenames,
		);
	
  ## dangit, preserve global sets; bad now missing .mrna2tsa.info sets for m2t.TSADESC, other ??
  #** BUG, need to reset single vars also.??
#   map{ my $ov= $settings{$_}; 
#     $mysettings{$_}= $ov if($ov and $ov ne $DEFAULT_SETTINGS{$_}); 
#     } sort keys %settings;
#   map{ my $ov= $ivset{$_}; 
#     $mysettings{$_}= $ov if($ov and $ov ne $DEFAULT_SETTINGS{$_}); 
#     } sort keys %ivset;
  
  if($action =~ /restore|read/ and -f $trpname) {
    open(my $inh, $trpname); # or loggit(warn ..);
    while(<$inh>) { chomp; if(s/^$PRES.//) { 
      my($k,$v)=split /\s*=\s*/,$_,2; 
      $v=~s/;*$//; $v=~s/\s*\#.*$//; # should I chomp trailing '#' comments? 
      my $ov= $mysettings{$k}; 
      unless($ov and $ov ne $DEFAULT_SETTINGS{$k}) { 
        ## fixme: need to reset global defaults
        $organism=  $v if($k eq 'organism');
        do { $v=~s/ +;/;/g; $v=~s/ +/;/g; $sraids= $v; } if($k eq 'sraids');
        $IDPREFIX=  $v if($k eq 'IDPREFIX');
        $TSADESC=   $v if($k eq 'TSADESC');
        $DATE=      $v if($k eq 'DATE');
        $MINSIZE=   $v if($k eq 'MINSIZE'); # require int
        $MAXGAP=    $v if($k eq 'MAXGAP'); # require int
        $genenames=  $v if($k eq 'genenames');
        $BioProject= $v if($k eq 'BioProject'); # FIXME lost
        # $trclass=  $v if($k eq 'trclass');
        # $cdnaseq=  $v if($k eq 'mrna'); # fixme: key != varname
        $mysettings{$k}=$v; # after possible changes
        } 
      }
    } close($inh);
    
    %settings = %mysettings; # make global now; ONLY after restore|read ??

  } else { # copy new single vars to settings..
    for my $ks (sort keys %mysettings) { 
      my $ov= $mysettings{$ks};  
      $settings{$ks}= $ov if($ov and $ov ne $DEFAULT_SETTINGS{$ks});    
      } ;
  }

  my $settings= join "\n", map{ "$PRES.$_=".$settings{$_} } sort keys %settings;
  if($action =~ /log/) { loggit(0, "mrna2tsa.info:\n$settings");  }
  if($action =~ /save/) { open(my $outh, '>', $trpname); print $outh $settings,"\n"; close($outh); }
}


=item FIXME: add evgrutrorf.sh method to get_evgtrset()

  #!/bin/bash
  ## evgrutrorf.sh  : fix make okayset/*.utrorf.mrna to go w/ utrorf.aa,cds
  
  evigene=/bio/bio-grid/mb/evigene/
  ptnames=`/bin/ls -1 {banana,catfish,litova,locust,pogonus,shrimp,whitefly,ztick}*/*.names`
  ptdirs=`echo $ptnames | sed 's,/.*,,g;'`
  thisdir=`pwd`
  
  for ptn in $ptnames; do {
   ptd=`echo $ptn | sed 's,/.*,,g;'`
   pt=`basename $ptn .names`
   ptar=outz/$ptd.tar
   if [ ! -f $ptar ]; then 
     if [ -f outz/$ptd.tar.gz ]; then ptar=outz/$ptd.tar.gz; fi
   fi
  
   if [ -f $ptd/okayset/$pt.utrorf.mrna ]; then continue; fi
   echo "# make $ptd/okayset/$pt.utrorf.mrna";
  
   if [ ! -f $ptd/dropset/$pt.drop.tr.gz ]; then
    gtar -xvf $ptar $ptd/dropset/$pt.drop.tr.gz
   fi
  
   cd $ptd/okayset/
   gunzip -c $pt*{okay,okalt}.aa.gz | perl -ne'if(/^>/) { $ok=(m/utrorf /)?1:0; } print if($ok);' > $pt.utrorf.aain
  
   $evigene/scripts/prot/traa2cds.pl -utrorf -trout -aa $pt.utrorf.aain -out $pt.utrorf.mrna  -log \
       -cdna $pt.{okay,okalt}.tr.gz ../dropset/$pt.drop.tr.gz
  
     ## and remake .aa,.cds for hopefully small changes..
   $evigene/scripts/cdna_bestorf.pl -nostop -noutrorf -act fwdfasta -cdna $pt.utrorf.mrna \
        -outaa $pt.utrorf.aa -outcds $pt.utrorf.cds
  
    cd $thisdir
  } 
  done

=cut

sub getmRNA_utrorf
{
  my($okpath,$trname)= @_;
  # input: okayset/{okay,okalt}.aa.gz; okayset/{okay,okalt}.tr.gz
  # what of dropset/drop.tr ??
  # output  okayset/xxx.utrorf.mrna .. maybe also regen xxx.utrorf.{aa,cds}, droping old okay.aa utrorf

  #** 2018.06 OBSOLETE getmRNA_utrorf(), now done in getmRNA() via
  #   traa2cds -utrorf  .. updated traa2cds does this along with -mrnaout for both cds/utr orfs
  my $mrnaset="$okpath/$trname.mrna"; #** 2018.06 okayset/mrna should contain all utrorf.mrna
  return($mrnaset) if(-s $mrnaset);
  
  ##.. looks like getmRNA()
  my $cdnatmp="$okpath/$trname.utrorf.mrna";
  return($cdnatmp) if(-s $cdnatmp);
  
  my ($ok,$err,$ninutr,$nin,$hin);
  my $TRSUFFIX='tr|cdna|fasta|fna'; # is this enough? NOT .mrna|cds|aa,  
  #was: getOkFileset($okpath,$TRSUFFIX);
  my($oktr,$alttr,$okd)= getOkFileset($okpath,$TRSUFFIX,undef,$trname);
  
  ## ** FIXME: need dropset/trname.drop.tr input also, many utrorf tr are here..
  ## FIXME2: drop.tr can be huge file, any way to cut down parsing? 
  ##  -- pull okay + okalt utrorf ids 1 time, scan drop.tr 1 time
  
  my $droptr=$oktr; $droptr =~ s/okay/drop/g; # okayset/name.okay.tr >> dropset/name.drop.tr
  $droptr="" unless(-s $droptr);

  ## * SHOULD replace okay.aa/utrorfs with newly gen cdna_bestorf utrorf.mrna ***
  system("touch $cdnatmp"); #? 
  my($okaa) = grep /.okay\.aa$|.okay\.aa.gz$/, grep/$trname/, @$okd;  
  if($okaa) {
  $nin=0; ($ok,$hin)= openRead($okaa); open(FO,">$okaa.utrorf");
  while(<$hin>){ if(/^>/){ $ok=(m/^>\S+utrorf/)?1:0; $nin++; } print FO $_ if($ok); } close(FO); close($hin);
  if($nin>0) { $err= runcmd("$APPtraa2cds -utrorf -trout -aa $okaa.utrorf -out stdout -cdna $oktr $droptr >> $cdnatmp"); }
  $ninutr+=$nin; unlink("$okaa.utrorf"); # $okall++ unless($err);
  }
  
  my($altaa)= grep /.okalt\.aa$|.okalt\.aa.gz$/, grep/$trname/, @$okd; # upd: grep/$trname/, 
  if($altaa) {
  $nin=0; ($ok,$hin)= openRead($altaa); open(FO,">$altaa.utrorf");
  while(<$hin>){ if(/^>/){ $ok=(m/^>\S+utrorf/)?1:0; $nin++; } print FO $_ if($ok); } close(FO); close($hin);
  if($nin>0) { $err= runcmd("$APPtraa2cds -utrorf -trout -aa $altaa.utrorf -out stdout -cdna $alttr $droptr >> $cdnatmp"); }
  $ninutr+=$nin; unlink("$altaa.utrorf"); # $okall++ unless($err);
  # push(@tmpfiles,"$okaa.utrorf","$altaa.utrorf");
  }

## ??? need this: and remake .aa,.cds for hopefully small changes.. but need to remove same IDs from okayset/*.{aa,cds} 
## $evigene/scripts/cdna_bestorf.pl -nostop -noutrorf -act fwdfasta -cdna $cdnatmp -outaa -outcds 
#  my $APPcdnabest= findevigeneapp("$EVIGENES/cdna_bestorf.pl"); # allow ENV/path substitutions?
#  if($ninutr>0) { $err= runcmd("$APPcdnabest -nostop -noutrorf -act fwdfasta -cdna $cdnatmp -outaa -outcds "); }
  
	loggit(0,"getmRNA_utrorf n=$ninutr in",$cdnatmp);
  return($cdnatmp,$ninutr);
}


=item get_evgtrset

  my($cdnaseq,$trpath,$trname)= get_evgtrset($trclass,$cdnaseq,);

  add: parse sra_result.csv if exists, for species, sraids .. other metadata
  see new in asmrna_trimvec.pl
    
=cut

# move get_evgtrset ?? to subs.pm also getmRNA()
sub get_evgtrset { 
  my($trclass,$cdnaseq,$pubdir)= @_;
  my($trpath,$trname,$nsra,$sradatah)=("","",0,undef);
	my $notokay=0;
  
  if($cdnaseq) { 
    $notokay=1; # dont look in okayset/? look in $pubdir now?
    $trclass= makename($cdnaseq,".trclass") unless($trclass); 
  }

  if($trclass) {
    my $trpname= makename($trclass,"","trclass"); 
    if($trpname =~ m,/,) { ($trpath,$trname)= $trpname =~ m,(.*)/([^/]+)$,; } # was BAD?
    else { $trname=$trpname; }
    $trpath ||= '.';  
    
   	# my $okpath="$trpath/okayset"; # fixme change to publicset for output, okayset for input
    my $okpath= ($notokay) ? $trpath :"$trpath/okayset"; # should I check if curdir has okayset files?

    # FIXME: 201405 here? insert getmRNA_utrorf() before getmRNA() which looks for .utrorf.mrna
    #** 2018.06 OBSOLETE getmRNA_utrorf(), now done in getmRNA() via
    #a.use constant OBSOLETE_getmRNA_utrorf => 1;
    unless(OBSOLETE_getmRNA_utrorf) {
    my($mrnaOfUtrorf,$nutrorf)= getmRNA_utrorf($okpath,$trname);
    }
    
    #? check if input cdnaseq is actually mRNA seq ?    
    ($cdnaseq)= getmRNA($okpath,$trname,$pubdir) if(!$cdnaseq and -d $okpath);
  }
  
  return($cdnaseq,$trpath,$trname,$sradatah);
}


sub get_srainfo {
	my($trpath,$trname)= @_;
	
	my $sracvs="$trpath/$trname.sra_result.csv";    
	$sracvs="$trpath/sra_result.csv" unless(-f $sracvs);  
	my($nsra,$sradatah)= parse_sra_result_cvs($sracvs) if(-f $sracvs);
	loggit(0,"sra_result from",$sracvs,"nresult=",$nsra);
	
	## add from where? mrna2tsa.info global %settings
	## $sradatah->{'assemblers'} == "Velvet/Oases v; SOAPDenovoTrans v; Trinity v;"
  # my @SRAK=("Experiment Accession","Organism Name","Instrument", "Submitter", "Total Size, Mb","Total Spots");

  # for my $ks ("Assemblers", "Instrument", "Total Size, Mb","Total Spots",'Total Assemblies') {
  #   if($settings{$ks} and not $sradatah->{$ks}) { $sradatah->{$ks} = $settings{$ks}; }
  #   elsif($sradatah->{$ks}) { $settings{$ks} = $sradatah->{$ks}; }
  # }
	my @SRAK=();
	if($sradatah->{cvsformat} == 1) {
  	@SRAK= ("Assemblers", "Instrument", "Total Size, Mb","Total Spots",'Total Assemblies');
	} else {
  	@SRAK= ("Assemblers", "Platform", "size_MB","spots","Total Assemblies", "BioProject");
	}
	for my $ks (@SRAK) {
  	if($settings{$ks} and not $sradatah->{$ks}) { $sradatah->{$ks} = $settings{$ks}; }
  	elsif($sradatah->{$ks}) { $settings{$ks} = $sradatah->{$ks}; }
	}
	
	 ## also use sradatah to edit template evgr_tsamethods.cmt, evgr_tsadesc.cmt, evgr_tsasubmit.sbt ??
	my($nupinfo,$tsamethf,$tsadescf,$tsasubf)= tsa_infotemplates($trpath, $trname, $sradatah);
	loggit(0,"info updated $nupinfo TSADESC=",$TSADESC);
	#NOT here# push @submitset, $tsamethf,$tsadescf,$tsasubf; #?? dont move unless tbl2asn gets path
 
	unless($genenames) { #? put in get_evgtrset
		my $gnt="$trpath/$trname.names";  $genenames=$gnt if(-s $gnt); 
		loggit(LOG_WARN, "Missing product -names",$gnt) unless($genenames);
	}
	
	return($sradatah);
}


sub get_trimvecset {
  my($mrnaf,$skipit)= @_;
    
  my $trimtag= "uvcut";  # $ENV{trimset}|| xxx  global?? in subs.pm ?
	my $trimdir="trimset";  ## look here? and tmpfiles? and mrna dir ?
	my @trimsuf= ("mrna", "aa", "cds", "ids");
  my $tname= makename($mrnaf,"");  
  my($outf,$outaa,$outcds,$outids)= map{ "$tname.$trimtag.$_" } @trimsuf;

  my $flagfile= makename($mrnaf,'.trimvec_done'); 
	return() if( -f $flagfile); # dont need trim files
  
  unless(-f $outf and -f $outids) {
  	my($okd,@files)= getFileset($trimdir,$trimtag); ## OK
  	#getFileset(trimset,uvcut,)= dir:6, suf:trimset/locust1all5asm.uvcut.aa trimset/locust1all5asm.uvcut.cds trimset/locust1all5asm.uvcut.ids trimset/locust1all5asm.uvcut.mrna 
    my @okd=@$okd;
  	($outf,$outaa,$outcds,$outids)= map{ my $s=$_; my($f)= grep /\.$s/, @files; $f||""; } @trimsuf;  
	}
	  
  if(-f $outf and -f $outids) {
		return($outf, $outaa, $outcds, $outids);
  } elsif($skipit) {
		return(); # return($outf, $outaa, $outcds, $outids); # or (); what?
	} else {
		loggit(1,"ERR: missing trimvec set: $outf,$outids; run asmrna_trimvec= $APPtrimvec");
		return;
	}
}

# see new in asmrna_trimvec.pl # move to cdnasubs.pm ?
# sub getmRNA_OLD  ## moved to cdna_evigenesub.pm

# my($pubid_format,$altid_format,$GBPROID)= make_IDPREFIX();
sub make_IDPREFIX
{
  my($existingID)= @_;
  my $digits=6; # or 7?
  #?? defer IDPREFIX use till read sra_result.cvs ..
  #  get species/organism from that, make IDPREFIX from Spp.org. abbev. 
  
  ## FIXME here? IDPREF option ignored for sppname parsing.. DEFAULT_SETTINGS wrong? has new?
  # my $ChangeDefaultIdPrefix= ($IDPREFIX eq "EVGmRNA")?1:0; 
  my $ChangeDefaultIdPrefix= ($IDPREFIX eq $DEFAULT_SETTINGS{'IDPREFIX'}) ? 1:0;
  
  if($existingID and $existingID !~ m/^$IDPREFIX/) {
    $existingID=~ s/t\d+$//;
    my($prefix,$nums)= $existingID =~ m/^(\w+)(\d\d\d+)$/; ## prefix may have numbers, how many? 
    if($prefix) { $IDPREFIX= $prefix; $digits= length($nums) || 6; }
    $ChangeDefaultIdPrefix=0;
  }
  
  if($ChangeDefaultIdPrefix and not($organism eq "Noname" or $organism !~ m/\w\w/)) {
    my $prefix="";
    my($gen,$spp)=split /[\s\_]/, $organism; # my($gen,$spp)=split" ",$organism; 
    if($spp and length($gen)>1) { $prefix= ucfirst(substr($gen,0,3)) . lc(substr($spp,0,3)) . "EGm"; }
    else { $prefix= ucfirst(substr($organism,0,6)) . "EGm"; }
    ##else { $prefix=$IDPREFIX; } # make default?
    $IDPREFIX= $prefix;
  }
  my $nd= ( $IDPREFIX =~ s/(0\d+)$// ) ? length($1) : $digits;
  my $pubid_format = $IDPREFIX.'%0'.$nd.'d'; # $public_options{'publicid'} || "evgr000000";
  my $altid_format = 't%d'; # t%d or t%02d # $public_options{'altid'} || "t00";
  my $GBPROID= $IDPREFIX."_".$DATE; # "cacao11evigene_20120827";
  #? No?# unless($GDB_PREFIX) { $GDB_PREFIX= "gnl|$IDPREFIX|"; } #?
  return($pubid_format,$altid_format,$GBPROID);
}





=item sra_ftp2srr from FTPpath

  curl  'ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX124/SRX124618/'
  dr-xr-xr-x 1073741824 ftp      anonymous        0 Aug 12  2012 SRR424344
  
  >> best:
  curl -s -l 'ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX160/SRX160070/'
  SRR521835
  SRR521944
  
  wget -A 'sra' -r -l 2 -nv -nd -nH -np --spider \
   'ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX124/SRX124618/'
  14:25:43 URL: ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX124/SRX124618/ [74] -> ".listing" [1]
  14:25:44 URL: ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX124/SRX124618/SRR424344/ [74] -> ".listing" [1]
  14:25:44 URL: ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX124/SRX124618/SRR424344/SRR424344.sra [0] -> "SRR424344.sra" [1]

=cut

sub sra_ftp2srr {
  my(@ftps)= @_;  return () unless(@ftps);
  my $APPcurl= findapp('curl'); return () if($APPcurl =~ /MISSING/);
  my @srrs=(); 
  foreach my $ftppath (grep /^ftp:/, @ftps) {
  	my $cmd="$APPcurl -s -l $ftppath/"; 
  	loggit( ($dryrun) ? 1 : 0,"CMD=",$cmd);  
    my $srrs= `$cmd`; # or http: ?? ## runcmd() instead ?
    push @srrs, grep /SRR/, map{ m/(SRR\w+)/; $1; } split " ",$srrs;
    }
  return @srrs;
}


=item parse_sra_result_cvs

  add: parse sra_result.csv if exists, for species, sraids .. other metadata
  evgr2tsa/litova1all3f/
    whiteshrimp_sra_result.csv # rename litova1all3.sra_result.csv

  "Experiment Accession","Experiment Title","Organism Name","Instrument",
    "Submitter","Study Accession","Study Title","Sample Accession","Sample Title",
    "Total Size, Mb","Total RUNs","Total Spots","Total Bases","FTP Path to Experiment",
    "Library Name","Library Strategy","Library Source","Library Selection"
    
  "SRX098246","Litopenaeus vannamei  transcriptome","Litopenaeus vannamei","Illumina HiSeq 2000","BGI",
    "SRP008317","BGI Litopenaeus vannamei transcriptome sequencing","SRS265043","Pacific white shrimp",
    "1515.37","1","13697473","2465545140","ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX098/SRX098246",
    "Ex-zai-2_l1","RNA-Seq","TRANSCRIPTOMIC","cDNA"
  
=item sra_cvs 2017 format

Run,ReleaseDate,LoadDate,spots,bases,spots_with_mates,avgLength,size_MB,AssemblyName,download_path,
  Experiment,LibraryName,LibraryStrategy,LibrarySelection,LibrarySource,LibraryLayout,InsertSize,InsertDev,
  Platform,Model,SRAStudy,BioProject,Study_Pubmed_id,ProjectID,Sample,BioSample,SampleType,
  TaxID,ScientificName,SampleName,g1k_pop_code,source,g1k_analysis_group,Subject_ID,
  Sex,Disease,Tumor,Affection_Status,Analyte_Type,Histological_Type,Body_Site,CenterName,Submission,
  dbgap_study_accession,Consent,RunHash,ReadHash

SRR3247180,2016-10-04,2016-06-13,26018777,6504694250,26018777,250,3026,,
  https://sra-download.ncbi.nlm.nih.gov/srapub/SRR3247180,SRX1645097,,RNA-Seq,
  PCR,TRANSCRIPTOMIC,PAIRED,0,0,ILLUMINA,Illumina HiSeq 2500,
  SRP072048,PRJNA315720,2,315720,SRS1350203,SAMN04569208,simple,743457,
  Daphnia similoides,Parthenogenetic female,,,,,female,,no,,,,,HUAIBEI NORMAL UNIVERSITY,
  SRA388003,,public,49CAA8290AC1FD00A6A5F0C61EDD392E,1EFA3DFE6E3818211435373BFC191AE2
...

=cut

sub parse_sra_result_cvs
{
  my($sracvs)= @_;
  my($ngot,$nin,$nerr,$cvsformat)=(0) x 9;
  my (%sradata, @hd);
  
  # revise for this format
  #? are these fields fixed now? need to check if no hdr in sracvs
  my $SRA2CVSHDR='Run,ReleaseDate,LoadDate,spots,bases,spots_with_mates,avgLength,size_MB,AssemblyName,download_path,Experiment,LibraryName,LibraryStrategy,LibrarySelection,LibrarySource,LibraryLayout,InsertSize,InsertDev,Platform,Model,SRAStudy,BioProject,Study_Pubmed_id,ProjectID,Sample,BioSample,SampleType,TaxID,ScientificName,SampleName,g1k_pop_code,source,g1k_analysis_group,Subject_ID,Sex,Disease,Tumor,Affection_Status,Analyte_Type,Histological_Type,Body_Site,CenterName,Submission,dbgap_study_accession,Consent,RunHash,ReadHash';
  
  open( my $inh, $sracvs) or loggit(1,"ERR: parse_sra_result_cvs reading $sracvs");
  while(<$inh>) { 
    next unless(/,/ or /\t/); chomp; 
    if($cvsformat == 0) { # always/no header?
      if(/^Run,ReleaseDate/){  $cvsformat=2; }
      elsif(/^(SRR|ERR)\d+,/){ $cvsformat=2; 
        @hd=split",",$SRA2CVSHDR; # guess this col set
      }
      elsif(/Experiment Accession/){ $cvsformat=1; } # "Experiment..","xxx","yyy"
      elsif(/^\w+\t/){ $cvsformat=3; } # tabbed ?
      $sradata{cvsformat}= $cvsformat;
    }
    
    my @cols;
    if($cvsformat == 3) {
      @cols= split"\t",$_;
      if(/^Run/){ @hd=@cols; next; }
    } elsif($cvsformat == 2) {
      s/","/\t/g;  s/^"//; s/"$//; #quotes or not?
      s/,/\t/g;  @cols= split"\t",$_;
      if(/^Run/){ @hd=@cols; next; }
    } elsif($cvsformat == 1) {
      @cols= map{ s/^"//; s/"$//; $_; } split /\",\s*\"/, $_;
      if(/^Experiment Accession/){ @hd=@cols; next; }
    } else {
      $nerr++; # warn/die ??
      next;
    }
    
    $ngot++;
    for(my $i=0; $i<@cols; $i++) { 
      my $hd=$hd[$i]||$i; my $v=$cols[$i]; 
      $sradata{$hd}.="$v;" unless($sradata{$hd} and $sradata{$hd} eq "$v;"); 
      }
  } close($inh);
  
  foreach my $k (sort keys %sradata) { $sradata{$k}=~s/;$//; }
  
  #reset defaults:  
  # my($deforg,$defsra)=("Noname","SRR000000"); 
  my($deforg,$defsra)= ($DEFAULT_SETTINGS{'organism'},$DEFAULT_SETTINGS{'sraids'});
  my $sorg= $sradata{"ScientificName"} || $sradata{"Organism Name"};
  my $sids= $sradata{"Run"} || $sradata{"Experiment Accession"};
  # globals: set here?
  $organism= $sorg if($sorg and ($organism eq $deforg or $organism !~ m/\w\w/));
  $sraids=   $sids if($sids and ($sraids eq $defsra or $sraids !~ m/SRR/)); #? or always use sids?
  if($IDPREFIX eq $DEFAULTidpre and $sorg) {
    my($gn,$sp)= split/[_\s]/,$organism,2;
    if($sp and $gn){ $IDPREFIX= ucfirst(substr($gn,0,3)) . lc(substr($sp,0,3)) . 'EVm'; }
  }
  ## also use sradatah to edit template evgr_tsamethods.cmt, evgr_tsadesc.cmt, evgr_tsasubmit.sbt ??
  ## rewrite template .cmt, .sbt unless updated already.

  ## ALSO use 'FTP Path to Experiment' to turn SRX into SRR for picky ncbi tsa : wget or curl calls?
  my $HAVEsrrids= ($sraids =~ /SRR|ERR/ and $sraids ne $defsra)?1:0; ## Need to save to info/config file
  unless($HAVEsrrids) { # FIXME or 
    my @urls;
    @urls= grep /http|ftp/, split/;/, $sradata{"download_path"} || $sradata{"FTP Path to Experiment"};  
    #o @ftps= grep /ftp/, split/;/, $sradata{"FTP Path to Experiment"};  
    my @srr= sra_ftp2srr(@urls);
    if(@srr>0) {
      $sraids= join(";",@srr);
      $sradata{'SRAids'}= $sraids;
      loggit(1,"sra_id=",$sraids);
      }
  }
   
  return($ngot, \%sradata, $sraids); # sids ?
}

# sub parse_sra_result_cvs_OLD
# {
#   my($sracvs)= @_;
#   my($ngot,$nin)=(0,0);
#   my (%sradata, @hd);
#   
#   open( my $inh, $sracvs) or loggit(1,"ERR: parse_sra_result_cvs reading $sracvs");
#   while(<$inh>) { 
#     next unless(/,/); chomp; 
#     my @cols= map{ s/^"//; s/"$//; $_; } split /\",\s*\"/, $_;
#     if($cols[0] =~ /Experiment Accession/) { @hd= @cols; }
#     else { 
#       $ngot++;
#       for(my $i=0; $i<@cols; $i++) { 
#         my $hd=$hd[$i]||$i; my $v=$cols[$i]; 
#         $sradata{$hd}.="$v;" unless($sradata{$hd} and $sradata{$hd} eq "$v;"); 
#         }
#     }
#   } close($inh);
#   
#   foreach my $k (keys %sradata) { $sradata{$k}=~s/;$//; }
#   
#   #reset defaults:  
#   # my($deforg,$defsra)=("Noname","SRR000000"); 
#   my($deforg,$defsra)= ($DEFAULT_SETTINGS{'organism'},$DEFAULT_SETTINGS{'sraids'});
#   my($sorg,$sids)    = ($sradata{"Organism Name"},$sradata{"Experiment Accession"});
#   $organism= $sorg if($sorg and ($organism eq $deforg or $organism !~ m/\w\w/));
#   $sraids=   $sids if($sids and ($sraids eq $defsra or $sraids !~ m/SRR/));
# 
# 	## log all these?
#   my @SRAK=("Experiment Accession","Organism Name","Instrument",
#             "Submitter", "Total Size, Mb","Total Spots");
#   if($ngot>0) { my @v= map{ $_ .'='. $sradata{$_}.';'; } @SRAK; loggit(0,"sra_result:",@v); }
#   
#   ## also use sradatah to edit template evgr_tsamethods.cmt, evgr_tsadesc.cmt, evgr_tsasubmit.sbt ??
#   ## rewrite template .cmt, .sbt unless updated already.
# 
#   ## ALSO use 'FTP Path to Experiment' to turn SRX into SRR for picky ncbi tsa : wget or curl calls?
#   my $NEEDsrrids=0; 
#   $NEEDsrrids= ($sraids =~ /SRR/ and $sraids ne $defsra)?0:1; ## Need to save to info/config file
#   if($NEEDsrrids) {
#     my @ftps=  grep /ftp/, split/;/, $sradata{"FTP Path to Experiment"};  
#     my @srr= sra_ftp2srr(@ftps);
#     if(@srr>0) {
#       $sraids= join(";",@srr);
#       $sradata{'SRAids'}= $sraids;
#       loggit(1,"sra_id=",$sraids);
#       }
#   }
#    
#   return($ngot, \%sradata);
# } # parse_sra_result_cvs_OLD

sub read_annotab
{
  my($annotab,$savePubids)=@_;
	# my $annotab =  makename($mrnaseq,".ann.txt"); ## FIXME: change suffix: was .annotab
  my %annothash=(); ## make this global? use same for 2+ subs
	if($annotab) {
		my ($ok,$annoth)= openRead($annotab);
		while(<$annoth>) { 
		  next unless(/^\w/); chomp; 
			my ($pubid,$oid,@cols)= split"\t";  
			  ##? maybe change annotab format to kvtab, ncbiannot=value, allow special fields like note=xxxx
			# OPT to fill in %pubids from this table
			#? is annot tab bad? should save pubids, use input mrna instead??
			if($savePubids) {
			  my($gid,$alti)= $pubid=~m/(\w+t)(\d+)$/;
        $pubids{$oid}= $pubid; 
        $pubids{$pubid}= $pubid; 
        $pubidinfo{$pubid}=join("\t", $oid,$gid,$alti); #? ,$reclass,$aqual); #  
        }
			$annothash{$oid}= $annothash{$pubid}= $_;  #? both? mrnaseq in should be >oid; might be >pubid
		}	close($annoth);
	}
  return %annothash; # or \%hash ?
}

sub make_annotab # or make_publicset
{
  my($mrnaseq,$trclass,$skiprun)=@_;
	#FIXME: also re-write publicset/mrna,aa,cds w/ >pubid  old=oldid; name=xxx, other annots in header ...
  
  my($notr,$nocds)=(0,0);
	my $annotab =  makename($mrnaseq,".ann.txt"); ## FIXME: change suffix: was .annotab
	my $pubdir="publicset";
  if(not -f $annotab and $pubdir and -d $pubdir) {
  	 my($pubd,$ft)= getFileset($pubdir,'.ann.txt');  $annotab=$ft if($ft);
  }

  return($annotab,$notr,$nocds) if( -s $annotab or $skiprun); # or dryrun ..

  my($nnamed,$namin)= parse_genenames($genenames,NAME_NOEDITS); #?? NAME_NOEDITS option?
  loggit(0, "names $genenames n=$nnamed\n"); 

  my ($inh,$outh,$hd,$oid,$fa,$ok);
	($ok,$inh)= openRead($mrnaseq);
  $ok= open($outh,'>',$annotab) if($ok);
  unless($ok) { loggit(1,"ERR: make_annotab $mrnaseq TO $annotab"); return; }

  #FOXME: annotab had WRONG aalen/qual, from mrna but not updated uvcut.aa or utrorf.aa
  # .. need to parse pubset/*.aa not or as well as .mrna
  
  ## BUG 201405  here missed utrorf aaqual, from missing input aaseq next to mrnaseq .. makename() suffix bug?
  our %aahdr=(); 
  my $aaseq  =  makename($mrnaseq,".aa"); # aa.gz ?
  $aaseq="$aaseq.gz" if(! -f $aaseq and -f "$aaseq.gz");
  if(-f $aaseq) { # warn if missing ??
	  my($aok,$aah)= openRead($aaseq);
    while(<$aah>) { if(/^>(\S+)\s+(.+)$/) { $aahdr{$1}= $2; } } close($aah);
  } else {
    loggit(1,"Warn: make_annotab missing $aaseq");
  }

  sub maketblinfo {
    my($oid,$hd,$fa)= @_; our(%aahdr);
    my $lenfa= length($fa);
    my $tblinfo= parse_evgheader($oid,$hd,$lenfa); # uses %genenames ...
    if(($oid =~ /utrorf/ or $hd=~/uvcut=/) and (my $aahdr= $aahdr{$oid})) {
      if($hd=~/uvcut=([^;\s]+)/){ $tblinfo->{'uvcut'}=$1; }
      my $aainfo= parse_evgheader($oid,$aahdr,$lenfa);
      ## DANG, this is bad for orig ok.aa set (rev.tr), only do this for utrorf.aa  AND? trim/uvcut.aa  
      ## .. bad case for strand=- // cdsor=; check aaqual=\d+ is same on both?
      if($aainfo->{'cdsor'} eq '+') { 
        map{ $tblinfo->{$_}= $aainfo->{$_}; } qw(aaqual cdsoff cdsor);
      }
    }
    my $trgaps= $fa =~ tr/Nn/Nn/; $tblinfo->{trgaps}= $trgaps; #? 
    return $tblinfo;
  }  
  
  my %didoid=();
  my($itr,$otr,$ocds,$oerr)= (0) x 10;
  while(<$inh>) { 
    if(/^>(\S+)/) { my $d=$1; 
    	if($oid and ($pubids{$oid} or not $skipdropseqs))
      { # check for $didoid{$oid} ?
        my $tblinfo= maketblinfo($oid,$hd,$fa);
       	my($no)= putAnnot($outh,$oid,$tblinfo,$itr); 
	      $notr+= $no; $nocds+= $no; 
	      $didoid{$oid}++;
      }
      $oid=$d; $hd=$_; chomp($hd); 
      $fa=""; $itr++;
      }
 		elsif(/^\w/) { chomp; $fa.=$_; }
  } 
  
	if($oid and ($pubids{$oid} or not $skipdropseqs))
	{
    my $tblinfo= maketblinfo($oid,$hd,$fa);
  	my($no)= putAnnot($outh,$oid,$tblinfo,$itr); 
	  $notr+= $no; $nocds+= $no; $didoid{$oid}++;
	}
  close($inh); 
  
	## FIXME for utrorf.aa also in pubids{}, but not in mrnaseq ...
  ## check  $didoid{$oid} vs pubids{} 
  my @undone= grep { not( m/^$IDPREFIX/ or $didoid{$_}) } keys %pubids;
  if(@undone) {
		# my $pubaa =  makename($mrnaseq,".aa"); 
  	# if( ! -f $pubaa and -f "$pubaa.gz" ) { $pubaa.=".gz"; }
 		# ($ok,$inh)= openRead($pubaa); .. parse @undone aa and putAnnot ??
 	 	foreach $oid (sort @undone) {  
 	 	  my $tblinfo= annotab2tblinfo($oid,"",0,0);
 	 	  my($no)= putAnnot($outh,$oid,$tblinfo,++$itr); 
	 	  $nocds+= $no; $didoid{$oid}++; ## NOT: $notr+= $otr; 
		}
  }
  
  close($outh);
  return($annotab,$notr,$nocds); 
}


sub make_pubseq
{
  my($inseq,$intype,$annothash)=@_;  # call for each input mrna,aa,cds
  my($notr,$nocds)=(0,0);
	
	my $pubdir="publicset"; my $pubd=undef; my $ft;
#   if(not -f $annotab and $pubdir and -d $pubdir) {
#   	($pubd,$ft)= getFileset($pubdir,'.ann.txt',$pubd);  $annotab=$ft if($ft);
#   }

  # fixme: look in publicset/ for annot, submitset/ for outfa,tblout
	my $outfa= $inseq; $outfa =~ s/.gz$//; 
	$outfa.="_pub.fa"; my($pubsuf)= $outfa=~m/(\.\w+.fa)$/;
  unless(-s $outfa or not -d $pubdir) {	($pubd,$ft)= getFileset($pubdir,$pubsuf,$pubd);  $outfa=$ft if($ft); }
  return($outfa,$notr) if( -s $outfa); # or dryrun ..

  my ($inh,$outh,$tblh,$annoth,$hd,$oid,$fa,$ok);
	($ok,$inh)= openRead($inseq);
  $ok= open($outh,'>',$outfa) if($ok);
  ## unless($ok) { loggit(1,"ERR: make_pubseq $inseq TO $outfa"); return; }
  return(undef,0) unless($ok); # or dryrun ..
  
  my($itr,$otr,$ocds,$oerr)= (0) x 10;
  while(<$inh>) { 
    if(/^>(\S+)/) { my $d=$1; # is this oid or pubid ???
			if($fa and $oid and ($pubids{$oid} or not $skipdropseqs))
      {
     		my $annorec= annotab2tblinfo( $oid, $annothash->{$oid}, $hd, $FAHDRisValid); ## = $annothash{$oid}||{}; 
     	  ($otr,$oerr)= putPubSeq($outh,$intype,$oid,$hd,$fa,$annorec); 
	      $notr+= $otr; 
     	}
      $oid=$d; $hd=$_; chomp($hd); 
      $fa=""; $itr++;
      }
    elsif(/^\w/) { chomp; $fa.=$_; }
  } 
  
	if($fa and $oid and ($pubids{$oid} or not $skipdropseqs))
	{
 		my $annorec= annotab2tblinfo( $oid, $annothash->{$oid},$hd,$FAHDRisValid); ## = $annothash{$oid}||{}; 
  	($otr,$ocds,$oerr)= putPubSeq($outh,$intype,$oid,$hd,$fa,$annorec); 
	  $notr+= $otr; $nocds+= $ocds;
	}
  close($inh); close($outh);
  push @publicset, $outfa; push @tmpfiles, $inseq; #??
  return($outfa,$notr); 
}


sub make_tblfsaset
{
  my($mrnaseq,$annothash,$skiprun)=@_;
  my($notr,$nocds)=(0,0);
  my $outfa = ($output)  ? $output  : makename($mrnaseq,".fsa"); # was .fna ; tbl2asn wants fsa
  my $tblout= ($tblfile) ? $tblfile : makename($mrnaseq,".tbl");
  
	my $tsadir="submitset";
  if(not -f $tblout and $tsadir and -d $tsadir) {
  	 my($tsad,@files)= getFileset($tsadir,'tbl|fsa');   
   	 my($ft)= grep /\.tbl/, @$tsad; $tblout=$ft if($ft);
     ($ft)= grep /\.fsa/, @$tsad; $outfa=$ft if($ft);
  }

  # fixme: look in publicset/ for annot, submitset/ for outfa,tblout
  return($outfa,$tblout,$notr,$nocds) if( -s $outfa or $skiprun); # or dryrun ..

  my ($inh,$outh,$tblh,$annoth,$hd,$oid,$fa,$ok, $inaah, $peph);
  #FIXME: 1412: add .pep output, from mrnaseq>aaseq, using same protids in tbl
  my $aaseq= makename($mrnaseq,".aa"); # or .aa.gz ?
  my $pepout= (-s $aaseq) ? makename($tblout,".pep"):""; ## was BAD pepout Noname27640.pep

	($ok,$inh)= openRead($mrnaseq);
  $ok= open($outh,'>',$outfa) if($ok);
  $ok= open($tblh,'>',$tblout) if($ok);
  unless($ok) { loggit(1,"ERR: make_tblfsaset $mrnaseq TO $outfa,$tblout"); return; }
 
  ## do this at end? collect oids, protids from putTblFsa()
  my %protids=();
  # if($pepout) { ($ok,$inaah)= openRead($aaseq); $ok= open($peph,'>',$pepout) if($ok); }
  
  my($itr,$otr,$ocds,$oerr,$protid)= (0) x 10;
  while(<$inh>) { 
    if(/^>(\S+)/) { my $d=$1; # is this oid or pubid : EITHER **
			if($fa and $oid and ($pubids{$oid} or not $skipdropseqs))
      {
				# FIXME: annorec becomes tblinfo, not from parse_evgheader(fa.hdr)
				# FIXME2: ** check annotrec is uptodate, match oid at least w/ evgheader
				# .. add this check to annotab2tblinfo() .. option to return faheader version?
				# my $tblinfo= parse_evgheader($oid,$hd,length($fa)); 

     		my $annorec= annotab2tblinfo( $oid, $annothash->{$oid}, $hd, $FAHDRisValid); ## = $annothash{$oid}||{}; 
     	  ($otr,$ocds,$oerr,$protid)= putTblFsa($outh,$tblh,$oid,$hd,$fa,$annorec,$itr); 

     	  if($oerr ne "OK") { }
     	  $protids{$oid}= $protid if($protid);
	      $notr+= $otr; $nocds+= $ocds;
      }
      $oid=$d; $hd=$_; chomp($hd); 
      $fa=""; $itr++;
      }
    elsif(/^\w/) { chomp; $fa.=$_; }
  } 
  
  ##($otr,$ocds,$oerr)= putAnnoSeq($outh,$tblh,$annoth,$oid,$hd,$fa,$itr) if($fa); # last
	if($fa and $oid and ($pubids{$oid} or not $skipdropseqs))
	{
    # my $tblinfo= parse_evgheader($oid,$hd,length($fa)); 
 		my $annorec= annotab2tblinfo( $oid, $annothash->{$oid},$hd,$FAHDRisValid); ## = $annothash{$oid}||{}; 
		($otr,$ocds,$oerr,$protid)= putTblFsa($outh,$tblh,$oid,$hd,$fa,$annorec,$itr); 
    $protids{$oid}= $protid if($protid);
	  $notr+= $otr; $nocds+= $ocds;
	}
  close($inh); close($outh);
  
  ## do this at end? collect oids, protids from putTblFsa()
  if($pepout and scalar(%protids)) { 
    ($ok,$inh)= openRead($aaseq);  my %didp; # block dupids
    $ok= open($outh,'>',$pepout) if($ok); 
    while(<$inh>) {
      if(/^>(\S+)/){ my $d=$1; $ok=0; if(my $p= $protids{$d}){ s/>$d.*$/>$p/; $ok=($didp{$p}++)?0:1; } }
      print $outh $_ if $ok;    
    }
    close($inh); close($outh);
  }
  
  return($outfa,$tblout,$notr,$nocds,$pepout); # pepout ??
}



=item tr2main

 aabugs4/aabugs4qual/tsaevg/daphmag5_evgt2c_2013mar7.txt
 ## .. this main{md} isnt right either, as alt-2ndid can be other alt.

  grep okay $pt.trclass | perl -ne'($td,$ok,$cl,$md,$piad,$aq,$fl)=split; $cl=~s/a2$//; \
  ($pi,$pa,$pd)=split"/",$piad; $md=$pd if($pd); if($cl =~ /main|noclass/) { $main{$td}=$cl; } \
  else { $alt{$md}{$td}= $cl; $balt{$td}=$md; } $n++; \
  END { @amain= grep { not $main{$_} } sort keys %alt; \
  for $am (@amain) { $md= $balt{$am}; \
  if($main{$md}) { @at=keys %{$alt{$am}}; map{ $alt{$md}{$_}=$alt{$am}{$_} } @at; } \
  elsif($md) {  $main{$am}="NOMAIN";  } } \
  foreach $md (sort keys %main) { @ad= sort{$alt{$md}{$a} cmp $alt{$md}{$b}} keys %{$alt{$md}};\
  $ad=join",",map{ "$_/".$alt{$md}{$_} } @ad; $mc=$main{$md}; print join("\t",$md,$mc,$ad)."\n"; } }'\
    > okayset/$pt.mainalt.tab

=cut

sub trclass2maintab
{
  my($trclass,$pubdir, $okids)=@_;
  my $ntr=0;  my $nerr=0;
  my $mainindex= $pubidnum_start;
  ## okids = \%validids after merge filesets
  my $hasokids= ($okids and ref($okids) and scalar(%$okids))?1:0;
## FIXME4: utrorf buggered, ids not in okids or .mrna, but are okay in .aa,.cds
  
  my $maintab = makename($trclass,".mainalt.tab","trclass");  # > $pt.mainalt.tab
  my $pubidtab= makename($trclass,".pubids","trclass");   
  if(not -f $pubidtab and $pubdir and -d $pubdir) {
		# my(undef,undef,$pubd)= getOkFileset($pubdir);  
  	my($pubd,$ft);
  	($pubd,$ft)= getFileset($pubdir,'pubids',$pubd);  $pubidtab=$ft if($ft);  
  	($pubd,$ft)= getFileset($pubdir,'mainalt.tab',$pubd);  $maintab=$ft if($ft);  
  }
  return($maintab,$pubidtab,$mainindex,$ntr) if( -s $maintab and -s $pubidtab);# or dryrun ..

  my($ok,%main,%mainsize,%alt,%altsize,%maindrops,%altdrops,%didaltdrops,%balt,%drop,$outh,$outpubidh,$inh);
  ## $ok= open($inh,$trclass);
  ($ok,$inh)= openRead($trclass);
  $ok= open($outh,'>',$maintab) if($ok);
  $ok= open($outpubidh,'>',$pubidtab) if($ok);
  unless($ok) { loggit(1,"ERR: parse $trclass TO $maintab"); return; }

  ## FIXME: only althi are reliably locus alternates; altmid .. are more likely paralogs
  while(<$inh>) {
    next unless(/^\w/); chomp;
    my($td,$ok,$cl,$md,$piad,$aq,$fl)=split "\t";
    unless($cl and $md and $aq) { $nerr++; loggit(1,"ERR: trclass inline:$_"); next; } # what?? report?

		## fix in asmrna_dupfilter2.pl, these are bogus entries (dupl map locus)
		#m2t: ERR: trclass inline:Funhe2Exy3m004019t1_G2        drop    altmap  Fungr1EG3m004332t1      99/48           0,0,pflag:3
		## these are from new asmrna_dupfilter using xxx.eqgene table, should not be in trclass .. kfish2evg367mixx.trclass9:821

		my $dropit=0;
		## maybeok needs handling, drop? or keep/mark/cull ??
		## should revise asm dupfliter to avoid these maybes, includes 'refbest', others some good/bad
		## now all are from exoneq flag, all 'althi1|part|frag' classes; 5702 refbest(dups?) of 13059
		## ** maybe keep all maybes ; found some refbest/good being dropped here.. but maybes are large subset
		## should filter by hoscore if avail
		if($ok eq 'maybeok') { 
		   if($fl=~/refbest|refgood/) { $ok="okay"; $cl.="maybe"; } 
		   else { $ok="okay"; $cl.="maybe"; } # 
		}
		
    if($ok ne 'okay') { $cl=$ok.$cl; $drop{$td}=$cl; $dropit=1;  } # next unless($SHOWDROPS);  OPTION: include drops?
    elsif($hasokids and not $okids->{$td}) { $dropit=1; $cl='dropid'.$cl;  $drop{$td}=$cl; } # next unless($SHOWDROPS);  unless $td =~ /utrorf/ ??

    #old# my($pi,$pa,$pd)=split"/",$piad; $md=$pd if($pd =~ /^\w/); # update BUG here -sense not pd
    #bad#my($pi,$pa,$pd,$asense)=split"/",$piad; # THIS IS WRONG, asense before pd
    #bad#if($pd =~ /.sense/ and not $asense) { $asense=$pd; $pd=""; }
    #bad#else { $md=$pd if($pd =~ /^\w/); } # update BUG here -sense not pd

    ## 98/100/-sense/PitaEaR000975t24 << buggers.
    ## all piad entries should have pd and sense/asense fields, placeholder where needed '.' ?
    my($pi,$pa,$asense,$pd)=split"/",$piad; # asense before pd ** NEED To revise asm dupfilter to clarify
    map{ $_ ||=""; } ($pi,$pa,$asense,$pd);
    # piad == 100/100/self1 << self1 misused here.
    ## OR: my($pi,$pa,@pmore)= split"/",$piad;
    ## my $pd= ($pmore[-1] =~ /^\w/) ? $pmore[-1] : "";
    ## my ($asense)= grep /sense/, @pmore; $asense||="";
    
    if($asense =~ /sense/) { $pd="" unless($pd =~ /^\w/); }
    elsif($asense =~ /^\w/) { $pd=$asense; $asense=""; }
    $md=$pd if($pd =~ /^\w/ and not $pd=~/self/);
   	$cl=~s/a2$//;  #? dont need 

		## FIXME here?  new altmap class from eqgene hides main class now; all mains are reclassed altmap, 1 should be left as main.
 		if($dropit and not $SHOWDROPS) {
 			next unless($cl =~ /main|noclass/); # this is enough, skip dropalts
#  			if($cl =~ /main/) {  # keep in dropmain{} for NOCLASS resolves?  swap td,md as new main, if not drop?
#  				$maindrops{$td}= $md; #? main cause of NOCLASS ?? can this use SHOWDROPS instead to collect dropmain in alt{}?
#  				$alt{$md}{$td}= $cl; $balt{$td}=$md; # treat as if SHOWDROPS was on
#  			}
 		}
 		
 		my $aasize= ($aq =~ m/(\d+).(\d*)/) ? "$1.$2" : 0; # size.pCDS for best sort
    if($cl =~ /^main|^noclass/) { 
    	$main{$td}=$cl; $balt{$td}=$td;
    	$mainsize{$td}= $aasize; # ($aq =~ m/(\d+).(\d*)/) ? "$1.$2" : 0;
    } else { 
    	$alt{$md}{$td}= $cl; $balt{$td}=$md; 
    	$altsize{$td}= $aasize; 
    	# NOMAIN fix: always add dropmain here 
    }  
  }

	if($SHOWDROPS) {
		# drops from dropset/*.aa headers for perfect_dups, perfect_frags info not in trclass ..
		my $ndr=0;
  	my($dset,$droptr)= getFileset("$trpath/dropset",'drop.tr|drop.aa|drop.cdna');  # drop.cdna/.fa/.xxx ??
  	my($ok,$hin)= ($droptr) ? openRead($droptr) : (0,0);
  	if($ok) { 

  	## FIXME.160911: MISSING evgclass=,drop, listings in mainalt.tab .. should include all dropset/*.aa hdr ids

  		# >sobananav1k21loc12t1 ... evgclass=althia2,drop,match:bananavel2k55Loc28462t1,pct:100/59/bananavel2k55Loc28462t2; 
  		# >fungrvel2k35Loc50t1 ...  evgclass=perfectdup,drop,match:rfung1savelk1ptr061k35Loc45t1
  		# FIXME2: spurious drops getting to pubid table via this $md == main of drop td, but may also be a drop!
  		# below via @amain from alt{md} != this main?
  		## skip this for drops?? $balt{$td}=$md
  		## add new hash:  altdrops{md}{td}= drop.cl ?
  		## new problem: md here may well not be in mainlist or linked there .. will then leave
  		## these dropped main/alt out of mainalt.tab .. maybe right, these are items of no value?

  		while(<$hin>) { if(/^>(\S+)/) {  
  		  my $td=$1;
  			my($cl,$ok1,$md)= m/evgclass=(\w+),(\w+),match:([^\s;,]+)/;
  		 	if($cl and $md and not $balt{$td}) { 
    			$ok1="drop"; # ensure no bad cases
  			  my $mmd= $balt{$md}||$md; # locate new main
  		 	  $drop{$td}="$ok1.$cl,$md"; 
  		 	  $altdrops{$mmd}{$td}= $ok1.$cl; $ndr++; 
  		 	 }
  		 	   # BUG160911: altdrops{md} miss when md reclass to not-main;  not: $balt{$td}=$md; 
  		} 
  	} close($hin); }
		loggit(0,"trclass2maintab: dropset adds $ndr"); 
	}
	
  ## Fix MISSING main links, from alt to other alts ..  
  ## FIXME2: some of these NOMAIN are drops *** dont retain;
  ##   .. NOMAIN can be avoided for some, due to drop of main2 == main1, need to switch link to main1
  ## FIXME3: adding drop{xxx} has screwed up alt-links somewhere; ** STILL MESSED UP
  # ... use drop{xx} only on output?
	## noDROPS  : grep -c NOMAIN publicset/v9b/kfish2evg367mixx9b.mainalt.tab = 21350
	## SHOWDROPS: grep -c NOMAIN kfish2evg367mixx.mainalt.tab = 13653 << get rid of more 
	## .. most remaining NOMAIN from eqgene altmap, a few from missing utrorf
  
  my %hasmain;
  my @amain= grep { not $main{$_} } sort keys %alt; # dropmain here now only for SHOWDROPS !
  foreach my $am (@amain) { 
    my $md= $balt{$am} || $am; ## $md=$am if($drop{$md});
    	## check for drop-main2 > ok-main1 link before NOMAIN : is this safe?
    	## dropmain here now only for SHOWDROPS ! fix above
  	if(!$main{$md} and $drop{$md}) { my $md1= $balt{$md}||""; if($md1 and $main{$md1}) { $md=$md1; } }
  	
    if($main{$md}) { my @at=keys %{$alt{$am}}; map{ $alt{$md}{$_}=$alt{$am}{$_}; $balt{$_}=$md; } @at; } 
    elsif($md) { 
     $main{$am}="NOMAIN";  # FIXME: get rid of these by finding alt main
     } 
  }
  
  foreach my $td (keys %balt) {
    my $md= $balt{$td} || $td; 
    # $md=$td if($drop{$md}); # what ??
    	## check2 for drop-main2 > ok-main1 link before NOMAIN : is this safe?
    	## need $alt{$md1}= $td;
  	##? if(!$main{$md} and $drop{$md}) { my $md1= $balt{$md}||""; if($md1 and $main{$md1}) { $balt{$td}= $md= $md1; } }
    $main{$md}="NOMAIN" unless($main{$md});# FIXME: get rid of these by finding alt main
  }
     
  ## add headers to these:
  #originalID     MainClass  Alternates
  #Public_mRNA_ID         originalID      PublicGeneID    AltNum
  print $outh '#'.join("\t",qw(originalID MainClass Alternates))."\n";
  print $outpubidh '#'.join("\t",qw(Public_mRNA_ID originalID PublicGeneID AltNum))."\n"
    if($outpubidh);

## BUGGERs got some dups
# grep locust1Svel1K23L10145t2 publicset/*pubids
# LocmigEGm000098t1       locust1Svel1K23L10145t2 LocmigEGm000098 1
# LocmigEGm000098t4       locust1Svel1K23L10145t2 LocmigEGm000098 4
# grep locust1Svel1K23L10145t2 publicset/locust1all5asm.mainalt.tab << dup main + alt
# locust1Svel1K23L10145t2 NOMAIN  locust1tri1loc42162c0t1/althi,locust1sop4p2k39loc5402t1/althi,locust1Svel1K23L10145t2/althi1
# grep ^locust1Svel1K23L10145t2 *.trclass
# locust1Svel1K23L10145t2 okay    althi1  locust1tri1loc42363c0t1 100/99/locust1tri1loc42162c0t1  819,67%,complete        0,0,pflag:0
  

	my %doneid=();
  #above# my $mainindex= $pubidnum_start;
  #FIXME? this is where pubids are assigned, by sort on oid name, maybe change that to sort on aasize?
  # .. on otherhand sometimes oid name sort is informative
  # my @mainlist= sort keys %main;
  my @mainlist= sort{ $mainsize{$b} <=> $mainsize{$a} or $a cmp $b } keys %main;

  #upd1806 from evigene/scripts/genes/altbest2pubset.pl
  # 1708: preserveOldIds option ..  keep input gene ids as best feasible, including altnum?
  #  need all alts per main, check evg gene ids for 1 x many, also need check across loci for 2+ new loci w/ old geneids
  #  where cant preserve geneid for all of locus, do what? new geneid, or modify old?
  #  UPD: want/need main = t1, cant preserve old alt nums if alt classing changes ..
  
  my ($nnewids,$newidh)= ($preserveOldIds) ? preserveOldIds(\@mainlist, \%alt, \%drop, \%altsize) : (0, undef);
  
  foreach my $md (@mainlist) { 

    ## FIXME 201405 sort alts for pubid by aasize, here sort is by class (althi,drop,..), keep that, add aasize  
    #o# my @ad= sort{$alt{$md}{$a} cmp $alt{$md}{$b}} keys %{$alt{$md}}; # not here?  grep { ! $drop{$_} } 
    my @ad= sort{$alt{$md}{$a} cmp $alt{$md}{$b}
      or $altsize{$b} <=> $altsize{$a} or $a cmp $b } keys %{$alt{$md}}; 

    my $ad=join",",map{ "$_/".$alt{$md}{$_} } @ad; 
    my $mc=$main{$md}; 
    
    # FIXME 1608: change $mc eq "noclass" to "main" if alts exist, and vversa
    if($mc =~ /^NOMAIN/) { 
      # do below
    } elsif(@ad>0 and $mc !~ /^main/) {
      $mc =~ s/^\w+/main/;
      $main{$md}= $mc; 
    } elsif( @ad==0 and $mc !~ /^noclass/) {
      $mc =~ s/^\w+/noclass/;
      $main{$md}= $mc; 
    }
    
    # SHOWDROPS: add altdrops{} as per alt{}, only for mainalt.tab, keep out of pubids
    ## BUG HERE for altdrops md reclassified not main (alt or drop..)
    if($SHOWDROPS) {
    	my @add= sort{$altdrops{$md}{$a} cmp $altdrops{$md}{$b}} keys %{$altdrops{$md}};  
    	map{ $didaltdrops{$_}=1 } @add; # later dump not didaltdrops
    	my $add=join",",map{ "$_/".$altdrops{$md}{$_} } @add; 
    	$ad .= ",$add" if ($add);
    }
    
    print $outh join("\t",$md,$mc,$ad)."\n";  # mainalt.tab
    
    if($outpubidh) { # should be required ??
      # $mainindex++; # BUG: SHOWDROPS : move below drop{}, but alt needs this?
      my $ialt= 0; my $needmain=0;
      if($drop{$md} or $doneid{$md}) { $needmain=1; }
      else {
      	$mainindex++; $needmain=0; # BUG: move below drop{}
      	my ($pubmrnaid,$pubgeneid,$pubti)= make_pubid($md, $mainindex, ++$ialt, $newidh);
      	print $outpubidh join("\t",$pubmrnaid,$md,$pubgeneid,$pubti)."\n"; 
      	$ntr++; $doneid{$md}++;
      	}
      
      ## FIXME 201405 sort alts for pubid by aasize, here sort is by class (althi,drop,..), keep that, add aasize  
      my @sad= sort{ $altsize{$b} <=> $altsize{$a} or $a cmp $b } @ad;      
      foreach my $ad (@sad) {
        unless($drop{$ad} or $doneid{$ad}) {
      	if($needmain) { $mainindex++; $needmain=0; }  
        my ($altmrnaid,$altgeneid,$alti)= make_pubid($ad, $mainindex, ++$ialt, $newidh);
        print $outpubidh join("\t",$altmrnaid,$ad,$altgeneid,$alti)."\n"; 
        $ntr++; $doneid{$ad}++;
        }
      }
    }
  }
  
  if($SHOWDROPS) { # UPD.160911
    my @mains= sort keys %altdrops;
    for my $md (@mains) {
      my @misdrop= grep{ not $didaltdrops{$_} } keys %{$altdrops{$md}};
      if(@misdrop) {
        @misdrop= sort{$altdrops{$md}{$a} cmp $altdrops{$md}{$b}} @misdrop;  
        map{ $didaltdrops{$_}=1 } @misdrop; # later dump not didaltdrops
        my $misdrop= join",",map{ "$_/".$altdrops{$md}{$_} } @misdrop; 
        my $mc= $main{$md}||"NOMAINd"; 
        print $outh join("\t",$md,$mc,$misdrop)."\n";  # mainalt.tab
      }
    }
  }
    

  close($inh); close($outh);
  close($outpubidh) if($outpubidh);
  
  # id order bug w/ preserveOldIds .. need pubids sorted by geneid, tinum
  if($preserveOldIds and $nnewids) {
    ## print $outpubidh '#'.join("\t",qw(Public_mRNA_ID originalID PublicGeneID AltNum))."\n"
    ## dang,  need '#Public..' hdr at top
    my $err= runcmd("sort -k3,3 -k4,4n -k1,1 $pubidtab | grep -v '#' > $pubidtab.srt");    
    if(-s "$pubidtab.srt" and not $err) { 
      system("head -1 $pubidtab >  $pubidtab.hdr");
      unlink($pubidtab);
      system("cat $pubidtab.hdr $pubidtab.srt > $pubidtab");
      unlink("$pubidtab.hdr"); unlink("$pubidtab.srt");
    }
  }
  
#	# fill global %pubids, here or above ??
# 	my($ok,$hin)= openRead($pubids); 
# 	my $first=0;
# 	while(<$hin>) { next unless(/^\w/); chomp; 
# 		my($id,$oid,$gid,$alti)=split"\t";     
# 		$pubids{$oid}= $id; $pubids{$id}="$oid\t$gid\t$alti";
# 		if(++$first == 1 and $nalltr == 0) { # reset to existing IDPREFIX
# 			($pubid_format,$altid_format,$GBPROID)= make_IDPREFIX($id); 
# 		}
# 	} close($hin);

	push @publicset, $maintab, $pubidtab; #?
  return($maintab,$pubidtab,$mainindex,$ntr);  # return main,alt id hashes ....
}

#upd1806 from evigene/scripts/genes/altbest2pubset.pl
# change for evgmrna2tsa, read table of pubid,oid from old gene set to preserve
# @$mainlist is of oids, not old pubids;
#upd1708: preserveOldIds option ..  keep input gene ids as best feasible, including altnum?
#  need all alts per main, check evg gene ids for 1 x many, also need check across loci for 2+ new loci w/ old geneids
#  where cant preserve geneid for all of locus, do what? new geneid, or modify old?

sub preserveOldIds {
  my($mainlist, $altsOfMain, $drop, $altsize)= @_;
  my(%gids, %gnums, %gdone, %newids, %newpubids, $gprefix);
  my $nids=0;
  #above# 
  my $NEWIDpre=''; # 'n';
  
  ## prefix, gnum problem: have mixed gprefices .. diff gnum set for each
  ## FIXME  gnum only for $gprefat =~ /$IDPREOK/
  # change for evgmrna2tsa, read table of pubid,oid from old gene set to preserve
  # @$mainlist is of oids, not old pubids;

  my $GNEXTNUM=0;
  my (%pod,%idparts);
  if( -f $preserveOldIds and open(F, $preserveOldIds)) {
    while(<F>){ if(/^\w/){ 
      my($pd,$oid)=split; next unless($oid =~ /\w/);
      $pod{$oid}=$pd; 
      my($gd,$gpre,$gnum,$ti)=(0,0,0,0);
      if($pd =~ m/^(\w+[A-Za-su-z])(\d\d+)t(\d+)$/) { # basic evg id form
        ($gpre,$gnum,$ti)=($1,$2,$3);
        $gd= $gpre.$gnum;
        $gnums{$gnum}++;
        $idparts{$pd}= [$gpre,$gnum,$ti]; # unless($idparts{$pd}); # dups? shouldnt be
        $gprefix= $gpre unless($gprefix);
        # $gids{$gnum}{$ti}= $pd;
        # $newids{$oid}= $pd;
      }      
      } 
    } close(F);
  }
  return unless(%pod);
  
  # foreach my $md (@$mainlist) {
  #   my @okd = grep{ not($drop->{$_}) } ($md, keys %{$altsOfMain->{$md}} );
  #   for my $id (@okd) {
  #     my($gpre,$gnum,$gd,$ti,$isplit)= evigene_idparts($id);
  #     next unless($gnum and $gpre =~  m/$IDPREOK/);
  #     $gnums{$gnum}++; # $gids{$gd}{$ti}= $id;
  #     $gprefix= $gpre unless($gprefix);
  #     }
  # }
  
  my($glast)= sort{ $b <=> $a } keys %gnums;
  $GNEXTNUM= 9 + $glast;
  
  my $idformat= $NEWIDpre . $gprefix . '%06d'; #?  use pubid_format of make_IDPREFIX  ?
  my(%havepubg);

  foreach my $md (@$mainlist) {
  
    my @ad= sort{ $a cmp $b } keys %{$altsOfMain->{$md}}; 
    my @okd = grep{ not($drop->{$_}) } ($md,@ad);
  
    my ($mnum,$nd,$timax)=(0,0,0);
    for my $oid (@okd) {
      my $pubid= $pod{$oid} or next;
      # check all alts for diff old pubgene id?
      my($gpre,$gnum,$ti)= @{$idparts{$pubid}};
      $mnum= $gnum unless($mnum);
      if($ti > $timax) { $timax = $ti; } else { $ti= ++$timax; }
      if($gnum ne $mnum or $newpubids{$pubid}) {
        do {
          $pubid= $gprefix . $mnum . sprintf( $altid_format, $ti); # FIXME ti clash
          if( $newpubids{$pubid} ) { $ti= ++$timax; }
        } while( $newpubids{$pubid} );
      }
      $newids{$oid}= $pubid; $nd++; $nids++;
      $newpubids{$pubid}= $oid;
    }
    
    if($nd == 0 and @okd > 0) {
      $mnum= ++$GNEXTNUM;
    }
    if($nd < @okd) { # finish, some or all new per locus / mnum 
      for my $oid (@okd) {
        next if($newids{$oid});
        my($pubid,$ti);
        do {
          $ti= ++$timax;
          $pubid= $gprefix . $mnum . sprintf( $altid_format, $ti);
        } while( $newpubids{$pubid} );
        $newids{$oid}= $pubid; $nd++; $nids++;
        $newpubids{$pubid}= $oid;
      }
    }
  } # mainlist
  
  return($nids,\%newids,\%newpubids);
} # sub preserveOldIds


sub make_pubid # add preseveOld %newids here? use altnum w/ it or preserve old altnum?
{
  my($oid, $pubidnum, $altnum, $preservedIds)= @_;
  $pubidnum_start= $pubidnum; #? pubidnum == mainindex
  my $alti=$altnum;
  my($pubid,$pubgene);
  if($preservedIds and ref($preservedIds)) {
    if(my $pid= $preservedIds->{$oid}) {
      ($pubid,$pubgene)=($pid,$pid);
      if(my($pgene,$palt)= $pid=~m/(\w+)t(\d+)$/) {
        $pubgene= $pgene;
        $alti= $palt; # or caller's altnum ?
        # if($palt ne $altnum) ..
        #NOT here: $pubid   = $pubgene . sprintf( $altid_format, $alti); # or palt ?
      }
      return($pubid,$pubgene,$alti);
    }
  }
  
  $pubgene = sprintf( $pubid_format, $pubidnum); 
  $pubid   = $pubgene . sprintf( $altid_format, $alti);
  return($pubid,$pubgene,$alti);
}

# sub make_pubid_OLD
# {
#   my($oid, $mainindex, $altnum)= @_;
# 
#   ## use/check oid? keep global hash for pubid <=> oid ?
#   ## my $pubid_format = $IDPREFIX.'%06d'; # $public_options{'publicid'} || "evgr000000";
#   ## my $altid_format = 't%d'; # t%d or t%02d # $public_options{'altid'} || "t00";
# 
#   # $pubidnum_start++ if($altnum == 1); ## or $mainindex == last index?
#   # my $pubidnum= $pubidnum_start;  # ONLY if altnum == 1 ? or if not seen this case..
#   my $pubidnum= $mainindex;
#   $pubidnum_start= $pubidnum; #?
#   
#   my $pubgene = sprintf( $pubid_format, $pubidnum); 
#   my $pubid   = $pubgene . sprintf( $altid_format, $altnum);
#   return($pubid,$pubgene);
# }



sub putPubSeq # for mrna,aa,cds.public w/ rewritten headers????
{
  my ($outh, $stype, $oid, $hdr, $faseq, $tblinfo)=@_; 
	## stype = mRNA|cDNA, CDS, amino|protein|polypeptide
  my($ntrout,$ncdsout)=(0) x 10;
  my $pubid= $tblinfo->{'pubid'} || $oid; ## or what? err?
  # OPTION : return/drop any seq lacking pubid .. ie source okayset has extras..
  
  my($cdsoff,$gname,$dbxref,$aaqual,$trlen,$protid,$lotag,$namepct, $cddname)= 
      @{$tblinfo}{qw(cdsoff name dbxref aaqual trlen protid locustag namepct cdd)};
  map{ $_='' if($_ eq "na"); } ($protid,$lotag,$gname,$cddname,$dbxref,$aaqual);  
	$gname='' if($gname =~ /^($NAME_NONE)/i);

  my $DROPxtrakeys= 'i|val|flag|path|strand|organism|species';
	my $DROPevgkeys = 'type|Name|Dbxref|db_xref|aalen|aaqual|aaSize|clen|offs|cdsoff|oid';  # evigene hdr keys, will replace
  # evigene set: do below# aalen|aaqual|aaSize|clen|offs|strand';
	
	my $hdr2="";
 	$hdr =~ s/>\S+//; 
 	$hdr =~ s/\s*\b($DROPxtrakeys)=[^;\n]+;?//g;
 	$hdr =~ s/\s*\b($DROPevgkeys)=[^;\n]+;?//g;
 	$hdr2 .="type=$stype; ";
 	$hdr2 .="Name=$gname; " if($gname);
 	$hdr2 .="Dbxref=$dbxref; " if($dbxref);
 	$hdr2 .="aalen=$aaqual; clen=$trlen; offs=$cdsoff; oid=$oid; ";
 	$hdr2 .="organism=$organism; " if($organism and $organism ne $DEFAULT_SETTINGS{'organism'});
#  	#............
# 	$hdr =~ s/\s*\btype=\S+;?//i; $hdr2 .="type=$stype; ";
# 	$hdr =~ s/\s*\bName=[^;\t\n]*;?//i; $hdr2 .="Name=$gname; " if($gname); #
# 	$hdr =~ s/\s*\b(Dbxref|db_xref)=[^;]+;?//i; $hdr2 .="Dbxref=$dbxref; " if($dbxref);
# 	$hdr =~ s/\s*\b(aaqual|aalen|aaSize)=[^;]+;?//i;  $hdr2 .="aalen=$aaqual; ";
# 	$hdr =~ s/\s*\b(clen)=[^;]+;?//i;  $hdr2 .="clen=$trlen; ";
# 	$hdr =~ s/\s*\b(offs|cdsoff)=[^;]+;?//i;  $hdr2 .="offs=$cdsoff; ";
# 	$hdr =~ s/\s*\boid=[^;]+//i; $hdr2 .="oid=$oid; ";
#   ## FIXME: option add to hdr: organism=Xxx yyy;
# 	$hdr =~ s/\s*(organism|species)=[^;]+;?//i;  $hdr2 .="organism=$organism; " if($organism and $organism ne $DEFAULT_SETTINGS{'organism'});

	$hdr ="" unless($hdr =~ /\w/); # or drop all?
	$hdr2 =~ s/;\s*$//; $hdr=~s/^\W+/ /;
  $faseq =~ s/(.{60})/$1\n/g; 
  print $outh ">$pubid $hdr2;$hdr\n$faseq\n";  $ntrout++;
}


sub putCDSloc 
{
  my($cdsoff,$partial,$cdsphase,$mrnalen,$mrna)= @_;
  
  ## is cdsoff == "<123->456" allowed here?
  $cdsoff =~ s/:[+.-]$//; # drop strand if there.
  my($start,$stop)= split/[-]/,$cdsoff; # or $cdsoff =~ m/(\d+)-(\d+)/; # 
  my($pstart,$pstop)= ($start,$stop);
  my($p5,$p3,$codonstart,$mrnatrim,$mrnadif)= ("<",">",0,0,0);
  if($partial =~ /complete/) { }
  else {
    ## FIXME.1510: NCBI requires mrna-cds partial to abut ends of mrna,unless gaps, 
    ## ie. ofs=mrnastart+2,3 wrong and ofs=mrnaend-1,2 wrong; for start shift, also shift cdsphase/codonstart
    ## -- problem with this, of course, when mrna-*end* is +2 from partial3-cds end,
    ##   cannot extend cds to mrna end w/o calling 2-base codon, adding to protein :((
    ## -- answer: trim mrna-end to cds-end (only for +2 case or +1 also?)
    
    unless($partial =~ /partial[53]/) { $partial.="53"; } # both
    
    if($partial =~ /3/) {  # no problem if called here already?
      #bad# if($stop < $mrnalen and $stop > $mrnalen - 3) { $stop= $mrnalen; } ## FIXME.1510: cant do this way due to ncbi desire for partial codons at cds/mrna-end
      if($stop < $mrnalen and $stop > $mrnalen - 3) { 
        $mrnatrim= $mrnalen - $stop;  ## FIXME.1510:
      if(0) { 
        # Causes more problems : Warn: SEQ_FEAT.TerminalXDiscrepancy mismatch prot,cds..
        # without get  Warn: SEQ_INST.TerminalNs for same prot/cds set; 
        # problem derives from prior fix, gap trims of AGCTnnnnnGnnnnnA left  AGCTnGnA un-cut single n's
        if(substr($mrna,$stop-1,1) eq 'N') { $stop--; $mrnatrim++; $cdsoff="$start-$stop"; } ## was below: if($faseq=~m/[ACGT]N$/)
      }
        }
      elsif($stop <= $mrnalen - 3 and $mrna) { 
        ##?? problem w/ this case, trim inner=mrna-NNN to cds-end may be large ??
        # if(substr($mrna,$stop,1) eq 'N') { $eshift=1; } ## == stop-1 + 1, zero-origin
        # elsif(substr($mrna,$stop+1,1) eq 'N') { $eshift=2;  }  
        # $stop += $eshift; $cdsoff="$start-$stop"; 
        }
      # $mrna= substr($mrna, 0, $mrnalen - $mrnatrim);  # leave to caller.
      $pstop="$p3$stop"; 
      }
      
    if($partial =~ /5/) { # no problem if called here already?
      if($start > 1 and $start < 4) { $cdsphase += $start-1; $start=1; } ## FIXME.1510:
      elsif($start >= 4 and $mrna and substr($mrna,$start-2,1) ne 'N') {
        ## buggggg this is too messy, still ncbi checha errs; dont shift mrna but change leading 1,2 to NN?
        ## bugg here, offby1, NNNN preceding, only shift for what?
        if(substr($mrna,$start-3,1) eq 'N') {  ## == start-1 -1, zero-origin
          #x# $start -= 1; $cdsphase= ($cdsphase-1) % 3; 
          $mrnadif=1; substr($mrna,$start-2,1)= 'N'; # try this way
          } 
        elsif(substr($mrna,$start-4,1) eq 'N') { # shift-2
          #x# $start -= 2; $cdsphase= ($cdsphase-2) % 3; 
          $mrnadif=2; substr($mrna,$start-2,1)= 'N'; # try this way
          substr($mrna,$start-3,1)= 'N'; # try this way
          } 
      }
      $cdsoff="$start-$stop";
      
#////// bad vvv /////////////
#       my $bshift=0; 
#       if($start > 1 and $start < 4) { $bshift=$start-1; } ## FIXME.1510:
#       elsif($start >= 4 and $mrna) {
#         if(substr($mrna,$start-2,1) eq 'N') { $bshift=1; } ## == start-1 -1, zero-origin
#         elsif(substr($mrna,$start-3,1) eq 'N') { $bshift=2;  } # shift-2
#       }
#       if($bshift>0) {
#         # $cdsphase += $bshift; # ugh wrong way? no; ugh ugh .. bshift wrong way for start>4 vs start 2,3
#         $cdsphase -= $bshift; $cdsphase = $cdsphase % 3; # convert -2 > 1, -1 > 2 << only for start>3
#         $start -= $bshift; $cdsoff="$start-$stop"; 
#       }
#--------------      
      $codonstart=$cdsphase+1; # is this right?
      $pstart="$p5$start";  
    }
  }
  my $tbl= join("\t",$pstart,$pstop,"CDS\n");
  $tbl .= "\t\t\tcodon_start\t$codonstart\n" if($codonstart>0);
  return ($tbl,$mrnatrim,$cdsoff,$cdsphase,$mrnadif,$mrna);## NOT: ,$mrna
}


use constant FOR_NOTE => 2;
sub putTblFsa
{
  my ($outh, $tblh, $oid, $hdr, $faseq, $tblinfo)=@_; 

	# FIXME: annorec becomes tblinfo, not from faseq.hdr parsing here
  # my $tblinfo= parse_evgheader($oid,$hdr,length($faseq));
	# ADD global   $protids{$oid}= $protid; for .pep output
  # OPTION: isoform attrs, use only pubid t[1..n] syntax? 
  #   need also protein == and not= info (isoform A,B,C) vs same isoform/diff mrna
  #  pubgenes/kfish2rae5h.puballnr.isoforms .. has needed attr
  ## fixme: ann.txt cdsoff has dang :+/- screws up tbl and mrna<>tbl .. fix in ann.txt
  
  my($ntrout,$ncdsout,$nfatrim)=(0) x 10;
  my $pubid= $tblinfo->{'pubid'} || $oid; ## or what? err?
  my $falen= length($faseq);  
  
  my($cdsoff,$gname,$dbxref,$aaqual,$trlen,$protid,$lotag,$namepct, $cddname)= 
      @{$tblinfo}{qw(cdsoff name dbxref aaqual trlen protid locustag namepct cdd)};
  my $Selcstop= $tblinfo->{'Selcstop'} ||"";
  $namepct =~ s/,.*//; $namepct =~ s/%//;
  my $gnameref= $tblinfo->{nameref} || $genedbxref{$oid}; #?? $pubid  # also require dbx DBXREF_OK limited set
  
  $cdsoff =~ s/:[+.-]$//; # drop strand if there.
  my($cdsb,$cdse)= split/[-]/,$cdsoff; # dang, :[+-.] strand
  my $cdsphase=0; # unused here?
  my($aafull)= $aaqual =~ m/(complete|partial\w*)/; 
  my $aastop=  ($aafull =~ /partial3|partial$/)? 0 : 1;
  my $aastart= ($aafull =~ /partial5|partial$/)? 0 : 1;

  map{ $_='' if($_ eq "na"); } ($protid,$lotag,$gname,$cddname,$dbxref,$aaqual);  

  $aafull ||= $aaqual;
  my($cdnap,$cdsym);# prodname: $cdsym also ? add to comments?
  my $hascddname= ($cddname and $cddname ne $gname)?1:0;
  ($gname,$namepct)= productname($gname,$aaqual,$namepct, FOR_NCBI); #? add FOR_NCBI flag?
  ($cddname,$cdnap,$cdsym)= productname($cddname,"100,99%,complete",0, FOR_NCBI|FOR_NOTE) if($hascddname); #? add FOR_NCBI flag?
  
  if($didpubid{$pubid}++) { # catch dup seq ids !
    loggit(1,"ERR: dup seqid $pubid/$oid");
    return($ntrout,$ncdsout,"ERR",0);
  }
  unless($cdse>0) { # or $cdsoff =~ m/\d+\-/  # err, missing cds span for tbl; print fsa seq anyway?
    loggit(1,"ERR: missing cds-offset $pubid/$oid:",$cdsoff,$aaqual,$gname);
  }

=item fix shift partial cds at +1/+2 of mrna-ends, or mrna-gaps

    FIXME.1510: NCBI requires mrna-cds partial to abut ends of mrna,unless gaps, 
    ie. ofs=mrnastart+2,3 wrong and ofs=mrnaend-1,2 wrong; for start shift, also shift cdsphase/codonstart
    -- problem with this, of course, when mrna-*end* is +2 from partial3-cds end,
      cannot extend cds to mrna end w/o calling 2-base codon, adding to protein :((
    -- answer: trim mrna-end to cds-end (only for +2 case or +1 also?)
    
   double darn ncbi again, trim mrna end when +2 or +1 of cds-partial3 end
   .. need putCDSloc() logic before print faseq. 
   .. mrna length doesnt appear in outputs here to ncbi asn except in print faseq.

=item also fix inner cds +1,2 of NNN gaps
      
      inner gap cds-off by 1,2            ofs     trlen obeg         v'5  oend     v'3
  Dapma6tiEVm002532t2	  814,67%,partial5;808-3252	3646	808	NNNNNNNCA.GGA	3252	TGA.ACACTTTTG
  Dapma6tiEVm003622t3	  520,71%,partial5;344-1906	2183	344	NNNNNNNAT.CAT	1906	TGA.AATAAAAAT
  Dapma6tiEVm015968t20	130,98%,partial;7-396	     397	  7    AANNNA.AAA  396	AAA.
      -- extend inner partial cds'5 -1,-2 if Not N, but precedes NNN by 1,2
      
      3'end-N-bug                         ofs     trlen obeg         v'5  oend     v'3
  Dapma6tiEVm009340t53	264,71%,partial3;316-1107	1108	316	GTATTTTTA.ATG	1107	NAN.A
  Dapma6tiEVm001867t281	235,99%,partial3;1-705     706    1	         .ATG  705	AAN.C
  Dapma6tiEVm009340t60	203,92%,partial3;49-657    658   49	AAGAAGATA.ATG  657	NGN.G
      -- trim -1 cds'3 if == 1-N only
=cut
  
  $faseq= uc($faseq); # ugh, some fix gaps are 'nnn'.
  #x# my($cdsline,$fatrim,$cdsoff2,$cdsphase2)= putCDSloc($cdsoff,$aafull,$cdsphase,$falen,$faseq); 
  my($cdsline,$fatrim,$cdsoff2,$cdsphase2,$mrna2dif,$mrna2)= putCDSloc($cdsoff,$aafull,$cdsphase,$falen,$faseq); 
  if($cdsoff2 and $cdsoff2 ne $cdsoff) { $cdsoff= $cdsoff2; ($cdsb,$cdse)= split/[-]/,$cdsoff; } ## always reset?
  $cdsphase= $cdsphase2;
    ## mrna2dif problem: faseq change inner5'part, leading 1,2 AGCT before NNN, change to NN also,
  if($mrna2dif) { $faseq=$mrna2; $falen= length($faseq); $nfatrim++; } 
  if($fatrim>0) {
    $nfatrim++; # *should* log any changed mrna/cds for later check, as aa.qual table? id/alen/gap/aqual/tlen/cdsoff/oid
    $faseq= substr($faseq, 0, $falen - $fatrim);
    $falen= length($faseq);  
    ## now putCDSloc() set fatrim++: if($faseq=~m/[ACGT]N$/)   
  }
  ## if($trlen != $falen) { } #?? warn? use falen
  
  $faseq =~ s/(.{60})/$1\n/g; 
  print $outh ">$pubid\n$faseq\n";  $ntrout++;

  if($cdsoff =~ m/\d+\-/) { # allow for <123->456
    $ncdsout++; #? always same count as ntrout ?
    print $tblh ">Features\t$pubid\t$GBPROID\n\n"; #?? is this used now
    # >Features       Thecc1ER_sopcsc10rk23loc140t1     cacao11evigene_20120827

		unless($protid) { ## dont need protid in tblinfo if it is only this transform
			$protid= $pubid;
			$protid=~s/t(\d+)$/p$1/; # genbank requires diff protid from mrnaid
			$protid= $GDB_PREFIX.$protid if($GDB_PREFIX); #  protein_id  gnl|CacaoGD|Thecc1EG016762p1
			}
    #no# $protids{$oid}= $protid;

# test w/ locus_tag:    26 WARNING: SEQ_FEAT.GeneXrefWithoutGene
# WARNING: valid [SEQ_FEAT.GeneXrefWithoutGene] Feature has gene locus_tag cross-reference 
#  but no equivalent gene feature exists FEATURE: CDS: Dimethylaniline monooxygenase (N-oxide-forming) [lcl|Funhe2Exx11m008798t1:109-1893] [lcl|Funhe2Exx11m008798t1: raw, rna len= 2716] -> [gnl|Evigene|Funhe2Exx11m008798p1]
    
    my $ctab="\t\t\t";
    #above.now# my($cdsline)= putCDSloc($cdsoff,$aafull,$cdsphase,$falen); # add fakeb
    print $tblh $cdsline;
    print $tblh $ctab, "protein_id\t$protid\n" if($protid);

    # print $tblh $ctab,"locus_tag\t$lotag\n" if($lotag); # problems, needs gene entry for this..
    if($lotag) {
      my($alt)= ($pubid=~m/t(\d+)$/)?$1:0; # check for >1 alt?
      my $v="locus_tag:$lotag"; $v.=", isoform $alt" if($alt);
      print $tblh $ctab,"note\t$v\n";
    }
    print $tblh $ctab,"note\toriginal_id:$oid\n" if($oid);
    ## productname does: my $gnn=($gname eq $NAME_UNK or not $gname)?$NAME_UNKNCBI:$gname;
    print $tblh $ctab,"product\t$gname\n";
    print $tblh $ctab,"note\tconserved domain $cddname\n" if($hascddname); # use ORIG cdd not prodname
      ## cdd note: XXX domain-containing protein >> XXX 
      
    if($gname and $namepct >= $MIN_IDLIKE) {  # MIN_IDLIKE MIN_NAMEIDENT
    
      #UPD108: change to "product alignment is 98 to human gene [GSYM?] NP_nnnnn" or similar, merge db_xrefs
      my $pnote="product alignment is $namepct";
      if($gnameref) { $pnote.=" to $gnameref"; }
      print $tblh $ctab,"note\t$pnote\n"; 
      #old1805: print $tblh $ctab,"note\tproduct alignment:blastp is $namepct\n"; # NO: \%
      ## SUSPECT_PRODUCT_NAMES:  8 cds comments or protein descriptions contain '%'
      #see below# inference  similar to AA sequence:Phytozome:PGSC0003DMP400040846
    }
    
    if($Selcstop) { # add /transl_except=(pos:1002..1004,aa:Sec);
      my @ss=split",",$Selcstop;
      for my $sb (@ss) { my $se=$sb+2; print $tblh $ctab,"transl_except\t(pos:$sb..$se,aa:Sec)\n" if($sb>0); }
    }

=item DBXREF config

    see evigene2genbanktbl.pl
		evigene.conf:
			dbxref_recode = TAIR=arath Phytozome=poptr|vitvi|soybn|soltu|sorbi
			dbxref_other DBXMISSING  

   		dbxref_ok  CDD taxon DDBJ EMBL NCBI GeneID GI dbEST dbSNP TrEMBL Swiss-Prot  UNILIB InterPro ENSEMBL GO 
  	   + JGIDB ESTLIB APHIDBASE AntWeb BEETLEBASE dictyBase FLYBASE GDB GeneDB GRIN MGI PDB PFAM PGN SGD SGN TAIR VectorBase WormBase ZFIN

		dbxref_ok from http://www.ncbi.nlm.nih.gov/genbank/collab/db_xref
			-- this html not easily parsible.
    see also ncbicsrc/api/asn2gnb6.c superset of http://www.ncbi.nlm.nih.gov/collab/db_xref.html
			
		init db_xref:
			my $DBXOTHER= $configh->{general}->{dbxref_other} || "DBXMISSING"; # special code or blank?
			my %DBXREF_RECODE= (); # dbxref_recode = TAIR=arath Phytozome=poptr|vitvi|soybn|soltu|sorbi
			my %DBXREF_OK= ( taxon=>1, TAIR=>1, TrEMBL=>1, SwissProt=>1, ENSEMBL=>1, ); 
			if( my $dbxrefok= $configh->{general}->{dbxref_ok}) {
				my @pg= split/[,\s]+/, $dbxrefok; 
				map{ my($k,$v)=split/[=:]/,$_; $v=1 unless(defined $v); $DBXREF_OK{$k}=$v; } @pg; 
			}
			
			if( my $dbxrefre= $configh->{general}->{dbxref_recode}) {
				my @pg= split/[,\s]+/, $dbxrefre; 
				map{ my($v,$k)=split/[=:]/,$_; my @k=split/[|]/,$k; foreach my $k1 (@k) { $DBXREF_RECODE{$k1}=$v; } } @pg; 
			}

		process db_xref:
      my($d)= m/^(\w+):/; 
      if( $DBXREF_RECODE{$d} ) { $d= $DBXREF_RECODE{$d}; }
      unless($d) { $_= "$DBXOTHER:$_"; } elsif(not $DBXREF_OK{$d}) { s/$d:/$DBXOTHER:$d./; } 

=cut
	
    if($dbxref =~ /\w/) {
    	my @dbother=(); # dont pass NCBI valid db_xref:
      foreach my $dx (split",",$dbxref) { 
        # DBXREF_RECODE fix for ncbi's  dbx: restrictions
        # FIXME999: more problems w/ gene.names table having odd/local DBprefix:ID
        #   .. fix where? should have here list of valid NCBI/ISxxx db prefixes. from where?
        # .. isnameref not seen , namepct '%' problem? or gnameref wrong/miss?
        ## must be bad:  $dx eq $gnameref; $gnameref == "RepID:xxx,RefID:yyy"
        my $isnameref= ($namepct >= $MIN_NAMEIDENT and $gnameref =~ m/^$dx/)?1:0; # $dx eq $gnameref
        
        # ? is HUGO:symbol valid at ncbi? hugo HGNC:idnum is valid, not symbol?
        
        $dx =~ s/^UniRef/SwissProt:UniRef/; ## UniProtKB:UniRef/; # or /TrEMBL:UniRef/;  SwissProt
        $dx =~ s/UniProt:/SwissProt:/; ## TrEMBL: .. UniProtKB: # need list to check/replace as per 
        # is CDD:id ok? accepted by tbl2asn .. YES
        
        ## ==> catfish1all4cf/okayset/catfish1all4.mrna.tbl2asn.sumval <==
        ## 60800 WARNING: SEQ_FEAT.IllegalDbXref
        ##  db_xref DRERI:ENSDARG00000019365 ## >> ENSEMBL:  FIXME in where? catfish/zfish.names
        ## anything:ENS\w+\d\d+ should become ENSEMBL:
        
        if($dx =~ /:ENS/ and $dx !~ /^ENSEMBL:/) { $dx =~ s/^\w+:(ENS\w+\d\d+)$/ENSEMBL:$1/; }
        
        ## fixme also: somedb:XP_ NP_ are ncbi prot ids;   db_xref_other:mayzebr:XP_004557409.1;
        # my $NCBIdb="GeneID"; # "LocusID";#? PIR? LocusID>GeneID, ok
        if($dx =~ /:[XN]P_/ and $dx !~ /^GenBank:/) { $dx =~ s/^\w+:([XN]P_\w+[\.\d]+)$/GenBank:$1/; }
        #old:if($dx =~ /:[XN]P_/ and $dx !~ /^GeneID:/) { $dx =~ s/^\w+:([XN]P_\w+[\.\d]+)$/GeneID:$1/; }
	      	  # ^^ UPD201705: new tbl2asn whines GeneID for AA inference, but not for db_xref ??

        ## FIX2: arath:AT5G60040.1 << should be TAIR: dammit ; fix .names instead of here..
        # evigene2genbanktbl.pl sub reformatval()
 
        ## print $tblh $ctab,"db_xref\t$dx\n" if($dx=~/\w/ and $dx ne "na"); 
	      if($dx=~/\w/ and $dx ne "na") { 
		      my($db,$did)= $dx =~ m/^(\w[^:]+):(.+)/; 
	        if( $DBXREF_RECODE{$db} ) { $db= $DBXREF_RECODE{$db}; }
	        $did=~s,/\d+,,; # ugh, dromel:FBgn0036451/46, dappu1:E9GZG3_DAPPU/100 extra junk
	        my $dxr="$db:$did";
	        ##upd1805: GenBank: is bad for db_xref, okay for inference ??
      	  if( $DBXREF_OK{$db} ) { 
	      	  print $tblh $ctab,"db_xref\t$dxr\n"
	      	    unless($db eq "GenBank");  # ugh...
	      	  # UPD201705: new tbl2asn whines about db_xref GenBank|NCBI|GeneID:NP_nnn prot dbxrefs
	      	  # .. change db_xref to note 
            print $tblh $ctab,"inference\tsimilar to AA sequence:$dxr\n" if($isnameref);
	      	} else { 
	      	  push @dbother, $dxr; 
	      	}
	      }
      }
      if(@dbother) { # print as notes
    		map{ print $tblh $ctab,"note\tdb_xref_other:$_\n"; } @dbother;
      }
    }
    print $tblh "\n";
  }
  
  return ($ntrout,$ncdsout,"OK",$protid);  
}

# annotab2tblinfo parse putAnnot table row; make same rec as parse_evgheader()
# FIXME: use updated annotab2tblinfo() in evigene/scripts/evigene_pubsets.pm
sub annotab2tblinfo 
{
  my($oidin,$tabrow,$fahdr,$FAHDRisValid)= @_;
	my($pubid,$oid,$trlen,$cdsoff,$aaqual,$annogaps,$dbxref,$namepct,$gname,$cddname)
			= (0) x 19; 
  $fahdr||="";

## FIXME: allow variations in cols, use header, ID=PublicID
# 	#print $annoth join("\t",qw(PublicID OrigID TrLen CDSoff AAqual TrGaps Dbxref Namealign Product_Name CDD_Name))."\n" if($itr==1);
# 	#print $annoth join("\t",$pubid,$oid,$trlen,$cdsoff,$aaqual,$annogaps,$dbxref,$namepct,$gname,$cddname)."\n";
	## ?? ERR if no annorec? == utrorf in .aa,cds not in mrna

  my $noann=0;
  $oid= $oidin; $pubid= $pubids{$oid} || $oid; # ERR if missing? FIXME, oidin may be pubid, global %pubids has both pub/oid
  $cddname= $aaqual= $dbxref= 'na'; 
  my $genenameinfo=0;
  if( $genenames{$oid} ) { # give these precedence ??
    $gname  =  $genenames{$oid};
    $namepct=  $genenamepct{$oid} || 0;
    $dbxref =  $genedbxref{$oid}||"na";
    $cddname=  $cddnames{$oid}||"na";
    my %tblinfo= (pubid => $pubid, oid => $oid, name => $gname, 
      namepct => $namepct, nameref => $dbxref, cdd =>  $cddname);
    $genenameinfo= \%tblinfo;
  }

  $tabrow ||=""; # undef for some?
	if(ref($tabrow) =~ /ARRAY/) {
		($pubid,$oid,$trlen,$cdsoff,$aaqual,$annogaps,$dbxref,$namepct,$gname,$cddname)
 			= @$tabrow;
	} elsif($tabrow =~ /\t/) {
		chomp($tabrow);
		($pubid,$oid,$trlen,$cdsoff,$aaqual,$annogaps,$dbxref,$namepct,$gname,$cddname)
 			= split"\t",$tabrow;
	} else {
	  $noann=1; $oid= $oidin; 
	  $pubid= $pubids{$oid} || $oid; # ERR if missing? FIXME, oidin may be pubid, global %pubids has both pub/oid
# 		$cddname= $aaqual= $dbxref= 'na'; 
# 		if( $genenames{$oid} ) {  # done above now..
# 			$gname  =  $genenames{$oid};
# 			$namepct=  $genenamepct{$oid} || 0;
# 			$dbxref =  $genedbxref{$oid}||"na";
# 			$cddname=  $cddnames{$oid}||"na";
# 		}
	}		

  my $cdsor=1; if($cdsoff =~ s/:([+.-])$//) { $cdsor=$1; } # drop strand if there. set cdsor  ?
  $dbxref =~ s/,$//;	 		
	my %tblinfo= (pubid => $pubid, oid => $oid,  protid => 0, locustag => 0,
      aaqual => $aaqual, trlen => $trlen, trgaps => $annogaps, cdsoff => $cdsoff, cdsor => $cdsor, 
      name => $gname, namepct => $namepct, dbxref => $dbxref, cdd =>  $cddname,
      );  # added trgaps =>  $annogaps
  my $anninfo= \%tblinfo;
  
  # FIXME2: ** check annotrec is uptodate, match oid at least w/ evgheader
  # .. add this check to annotab2tblinfo() .. option to return faheader version?
    ## FIXME: Selcstop=  only in fahdr now ?
    # >Funhe2Exx11m129535t7 aalen=321,92%,complete,selcstop; Selcstop=655,676,694,742,748,754,802,832,838,883,889,910,934,940,961,967,; clen=1039; strand=+; offs=16-981;  pubid=Funhe2EKm029571t1; Name=Selenoprotein P,  
  if($fahdr) {
    my $fatinfo= parse_evgheader($oidin,$fahdr); 
    #?? special case for kfish2 x11 == pubid + kfish2 gene pubid=Funhe2EKm... want pubid= as other dbxref
    #.. should revise mRNA Dbxref= entry w/ it..
    if($fatinfo->{Selcstop}) { 
      $anninfo->{Selcstop}=$fatinfo->{Selcstop}; 
      $anninfo->{aaqual}=  $fatinfo->{aaqual}; # qual,selcstop
      }
    my $dif="";
    my @KS=qw(oid pubid cdsoff aaqual trlen); # fckn.mess cdsoff/off has :+/- sometimes or not
    @KS= grep{ not m/cdsoff/ }@KS if($fahdr=~/strand=/); # if($fatinfo->{strand} eq "-");
    ## mrna: aalen=1884,83%,complete; clen=6736; strand=+; offs=99-5753; 
    ## ann.txt: 6736    99-5753:+       1884,83%,complete 
    ## ** cdsoff ignore :+ end or fix fahdr reader
      # which are required to match? cdsoff, aaqual?, trlen? 
      # * add aaqual, need proper complete/partial flag
    for my $ks (@KS) {
      my $fvs=$fatinfo->{$ks}; my $avs=$anninfo->{$ks};
      if($fvs and $fvs ne $avs) {
       my $dok=($ks =~ /^(oid|cdsoff)/ and $avs and ($fvs=~m/$avs/ or $avs=~m/$fvs/))?1:0;
       $dif .="$ks=$fvs/$avs," unless($dok);
       }
    }
  #m2t: ERR: annot mismatch Funhe2Exx11m002607t1/Funhe2Exx11m002607t1: 
  #  oid=Funhe2Exx11m002607t1,Fungr1EG3m001115t1 / Funhe2EKm000004t1,Fungr1EG3m001115t1,
  ##  ^^ these are dif due to update .ann.txt vs .mrna diff policy ; fixme pubid=Funhe2EKm000004t1 vs xx11 id

  ## annot mismatch problems: non-mrna strand- cdsoff diff mrna cdsoff
  ## utrorf trlen differs from input.. should be fixed somewhere, but not in seqfile.fa hdr
  ## ?? okayset/xxxx.utrorf.mrna
      
    if($dif) {
      my $fatok=($FAHDRisValid or ($noann and $dif=~/cdsoff=/))?1:0;
      if($fatok) { 
        my $annold=$anninfo; # or use %tblinfo
        $anninfo=$fatinfo; #?? merge best of both ?? where dif is missing/0
        ## fix2: copy bits from wrong anninfo: dbxref,  name???? missing genenameinfo ?
        if( length($annold->{dbxref}) > length($anninfo->{dbxref}) and $annold->{dbxref}=~m/:/) {
          $anninfo->{dbxref}= $annold->{dbxref}; # ok?
          }
      }
      ## both info are missing cdsoff for some .. recalc?
      #** cancel this err for mrna uvcut= annots
      my $nodif=($noann or $fahdr=~/uvcut=/)?1:0; #?? $fatinfo->{'uvcut'}
      #Off: loggit(1,"ERR: annot mismatch $pubid/$oidin:",$dif) unless($nodif); 
        #**TOO MUCH annot mismatch n=686532 : due to tr<>mrna revcomp ? cdsoff changed 
    }
  }
  
  ## fix3: names precedence; copy cdd, name? namepct? nameref=dbxref
  if(ref $genenameinfo) {
    for my $k (qw(name namepct nameref cdd)) {
      if($genenameinfo->{$k}=~/\w/) { $anninfo->{$k}= $genenameinfo->{$k}; }
    }
    
    ## this needs to be in parse_evgheader, also
    my $adx= $anninfo->{dbxref} || ""; my @adx=split",",$adx;
    my $ndx= $genenameinfo->{nameref} || ""; my @ndx=split",",$ndx;
    ## fixme.. drop adx when gnameref is similar:  gnameref=SwissProt:TDRD7_HUMAN; adxref=TrEMBL:TDRD3_HUMAN,TrEMBL:UniRef50_Q9H7E2
    my @addx=(); my @dropx=();
    for my $nd (reverse @ndx) { 
      unless($adx =~ /$nd/) { unshift(@addx,$nd); 
        if($nd=~/SwissProt|TrEMBL|UniProt/) { my($sp)= $nd=~m/(_\w+)$/; 
          if($sp and $adx =~ /$sp/){ push(@dropx, grep( /$sp/, @adx)); @adx= grep{ not m/$sp/ } @adx; }
        }
      } 
      # unless($adx =~ /$nd/) { unshift(@adx,$nd); } 
    }
    unshift(@adx, @addx);
    $anninfo->{dbxref}= join",",@adx;
  }

  ## add pubid to locustagid for ncbi tsa
  #m2t.LOCUSTAG =     Locus tag prefix:       RC70 (SAMN02116571)    <<dgg is this ok? conflict w/ other kfish EST project?
  unless($anninfo->{locustag}) {
    if(my $lotagpre= $settings{'LOCUSTAG'}) {
      my($pubidnum)= $pubid =~ m/$IDPREFIX(\d+)/;
      $anninfo->{locustag}= $lotagpre.'_'.$pubidnum if($pubidnum);
    }
  }
  
	return $anninfo; # add $genenameinfo? adx dbxref vs ndx nameref? nameref should be 1st of dbxref
}

=item mixed dbxref confusion
      note    original_id:Funhe2Exx11m011441t6
      product Porcupine fragment
      note    conserved domain CDD: Membrane bound O-acyl transferase (MBOAT) family protein
      note    product alignment:blastp is 68
      db_xref ENSEMBL:ENSXMAP00000012887
      db_xref SwissProt:PORCN_HUMAN     << gnameref best
      db_xref CDD:215190
      db_xref TrEMBL:UniRef50_Q6ZNC8
      db_xref TrEMBL:MBOA1_HUMAN        << older dbxref
=cut
  
=item CDD2name

  convert CD name to ncbi/uniprot acceptable protein name
  NCBI eukgenosub_annotation "<domain|repeat>-containing protein". e.g. "PAS domain-containing protein 5".

   PH_IRS, Insulin receptor substrate (IRS) pleckstrin (PH) domain protein, putative       lcl|Funhe2Exx11m084853t1:24-236 
   ^^ need new CDD to NCBI name routine: cut out leading symbol XXX_yyy, for note?
      .. remove fragment,partial,-like,putative and such terms
      .. add ' domain-containing protein' unless has domain
   FOG: bug, whatisit? CDD: FOG: Zn-finger   CDD:227381 
   bad: "ZP, Zona pellucida  domain-like fragment domain-containing protein"

=cut

sub CDD2name
{
  local $_= shift;  my $fornote= shift or 0; 
  my $csym="";
  s/CDD:\s*//; s/FOG:\s*//;
  if(s/^(\S+),\s*(\w)/$2/) { $csym=$1; }
  s/(fragment|putative|\-like)\s*/ /g; s/_/ /g;
  ## Changed 'Domain of unknown function DUF4218' to 'protein of unknown function DUF4218' ..
  ## better productname: Unknown function DUFnnn domain-containing protein
  if(/^Domain of unknown function/i and not $fornote) {
    s/Domain of unknown function (\w.*)/Unknown function $1 domain-containing protein/i ; # use key $UNK_CDD ??
    s/Domain of unknown function/Protein of unknown function/i ; # use key $UNK_CDD ??
  }
  
  if(m/\bdomain/i) { s/\bdomain protein// if($fornote); }
  else { s/\bprotein//; $_.=" domain-containing protein" unless($fornote); }
  return (wantarray) ? ($_,$csym) : $_;
}

sub productname 
{
  my($gname,$aaqual,$napct, $FOR_NCBI)=@_;  #? add FOR_NCBI flag?
  local $_= $gname || ""; 
  $napct ||= 0; $FOR_NCBI ||=0; 
  my $cdsym="";
  ## see evigene/scripts/prot/protein_names.pm : nameclean() .. BUT these are supposed to have been run thru that.
  ## add option to call nameclean() here .. again?
  
  use constant NAMEFIX1_FRAG => 1; # FIXME: NCBI hates fragment now, wants ', partial' ??
  ## fix3? -like should never end name, or follow punctuation .. always "blah blah-like protein" ??

  #noname? should this be done already?# 
  $_='' if($napct>0 and $napct < $MIN_IDLIKE);
  
  if($_ eq 'na' or /^($NAME_NONE)$/i or /^($NAME_NONE) protein\W*$/i) { 
    $_='';
  } elsif(/^($NAME_NONE)\s*\w+/) { 
    if(/^$NAME_UNK\s*\w+/) { s/^$NAME_UNK\s*//; }
    else { s/^($NAME_NONE)\s*//; s/^protein\s*//i; }
  }
  if(/\w\w/) {
    if(s/\s*\((\d+)\%\w*\)\W*$//) { my $p=$1; $napct=$p if($p>$napct); } # dang remove namepct: (99%P)
   
    ## TE: CDD: RP_RTVL_H_like, Retropepsin .. names bug, should leave off TE: 
    s/^TE:\s*//; # names bug, leave out TE: here
    my $fornote = ($FOR_NCBI & FOR_NOTE)?1:0;
    ($_,$cdsym)= CDD2name($_,$fornote) if(/^CDD:/);
    # if(/^CDD:/) { s/CDD:\s*//; s/\-like\s*/ /; unless(m/\bdomain/i) { s/\bprotein//; $_.=" domain-containing protein"; } }  
      
    # 45 WARNING: SEQ_FEAT.BadTrailingCharacter  Funhe2Exx11m001464t1
      # WARNING: valid [SEQ_FEAT.BadTrailingCharacter] Protein name ends with undesired character FEATURE:
      # "Roundabout, axon guidance receptor," n=32 are Roundabout
      #  "LATS, large tumor suppressor," << ','  n=13 rest are LATS,
      
    # 8 WARNING: SEQ_FEAT.ProteinNameEndsInBracket .. Funhe2Exx11m008798t7/Funhe2EKm026612t + other alts of same
      # Protein name ends with bracket and may contain organism name FEATURE: 
      # Dimethylaniline monooxygenase [N-oxide-forming]  << this appears to be Human HUGO approved name FMO3_HUMAN/FMO5_HUMAN
      # Funhe2EKm026612t6/Funhe2Exx11m008798t7 == TrEMBL:UniRef50_G7MFB6,TrEMBL:G7MFB6_MACMU,mayzebr:XP_004557286.1,Omcl:FISH11G_G24186

if(NAMEFIX1_FRAG) {  
    s/[,]?\s*(partial|fragment)\b//i; # always drop these? add , partial below     
} else {
    # 4 WARNING: SEQ_DESCR.InconsistentProteinTitle : "partial" name on complete aa
    if($aaqual =~ /complete/ and m/\b(partial|fragment)\b/i) { 
      s/[,]?\s*(partial|fragment)\b//i; 
    } elsif($aaqual =~ /partial/ and not m/(partial|fragment)/i) {
      $_ .= " fragment"; # FIXME: NCBI wants ', partial' ??  fragment is Uniprot's preferred syntax?
    }
}    
    # see protein_names:nameclean()
    s/^\W+//; # no leading crap
    s/\s*[\/\.,;_-]+\s*$//; # trailing punc; BadTrailingCharacter
if(NAMEFIX1_FRAG) { # add after trim; BUT NOT FOR_NCBI .. they add it ?? 
   unless($FOR_NCBI) { $_ .= ", partial" if(/\w\w/ and $aaqual =~ /partial/); } # FIXME: NCBI hates fragment now, wants ', partial' ??
}    
	}
	## change b/n NAME_UNKUNIP and NAME_UNKNCBI depending on outputs?  prefer NAME_UNKUNIP but for TSA submit
	unless(m/\w\w/) { $_= ($FOR_NCBI)? $NAME_UNKNCBI: $NAME_UNK; }  # SEQ_FEAT.MissingCDSproduct; need some name
	return($_,$napct,$cdsym);
}


sub putAnnot 
{
  my ($annoth, $oid, $tblinfo, $itr)=@_;  ## $hdr, 

	# FIXME: annorec becomes tblinfo, not from faseq.hdr parsing here
  # my $tblinfo= parse_evgheader($oid,$hdr,length($faseq));

  my $pubid= $tblinfo->{'pubid'} || $oid; ## or what? err?
  
  my($cdsoff,$gname,$dbxref,$aaqual,$trlen,$trgaps,$protid,$lotag,$namepct, $cddname)= 
      @{$tblinfo}{qw(cdsoff name dbxref aaqual trlen trgaps protid locustag namepct cdd)};

  # my($cdsb,$cdse)= split/[-]/,$cdsoff; 
  # my $cdsphase=0; # unused here?
	# my($aafull)= $aaqual =~ m/(complete|partial\w*)/; 
  # my $aastop=  ($aafull =~ /partial3|partial$/)? 0 : 1;
  # my $aastart= ($aafull =~ /partial5|partial$/)? 0 : 1;

# trgaps ==	my $annogaps="0"; # FIXME: check mrna/cds for NNN per old vers, for annotable
#   my $annogaps= ($nl==$ol and $nn==0) ? "0" : "gaps=$nn1/$nn";
#   $annogaps.= ",oldcds=".$tblinfo->{'cdsold'} if($tblinfo->{'cdsold'});
#   $annogaps.= ",cut=$ncut" if($ncut>0); 
#   $annogaps.= ",vectrim=$vectrimw" if($vectrimw); 
#   $annogaps = "TOOSHORT=$nl,$annogaps,cdsnn$CDShasTooManyXs" if($nl<$MINSIZE or $CDShasTooManyXs);

  ($gname,$namepct)= productname($gname,$aaqual,$namepct);

  ## FIXME: always annotate CDD name if exists.
  #always now: if($annoth) ..
    # usable output table ; FIXMEd: add header at top
    # NOTE: add vectrim info, other?; will need to regen .aa, .cds from .fsa output to be accurate..
    # .. user choice: ncbi submit restrictions may not be desired.
    # FIXmaybe: recode dbxref tags: as for .tbl ?? UniProt => SwissProt, xxxx:ENS => ENSEMBL:ENS ..
    
	print $annoth join("\t",qw(PublicID OrigID TrLen CDSoff AAqual TrGaps Dbxref Namealign Product_Name CDD_Name))."\n" if($itr==1);
	print $annoth join("\t",$pubid,$oid,$trlen,$cdsoff,$aaqual,$trgaps,$dbxref,$namepct,$gname,$cddname)."\n";
 
 	# FIXME here? add rewrite mrna, cds, aa seq files w/ pubid, some annots in headers
 	return(1);
}




#.......... tsa/tbl2asn subs : move out to separate app ?? ...................

=item tsa_infotemplates

  - Make these if needed; update for trasm name as needed

evgr_tsamethods.cmt:
  StructuredCommentPrefix ##Assembly-Data-START##
  Assembly Method EvidentialGene v2013.03.02; << from VERSION
     Velvet/Oases v1.2.03/o0.2.06 (2012.02); SOAPdenovo-Trans v2011.12.22; Trinity v2012.03.17  << USER inputs
  Assembly Name   evigeneR13     << from local myspecies.trclass name
  Sequencing Technology   Illumina   << from sra_result

evgr_tsadesc.cmt : insert parts of sra_result
  RNA-Seq data of [[this species]] are
  assembled de-novo with [[various RNA assemblers]], using multiple options
  for kmer fragmenting, insert sizes, read coverage, quality and
  abundance filtering. 
  EvidentialGene tr2aacds pipeline software is used
  to process the [[several million]] resulting assemblies by coding
  sequences, translate to proteins, score gene evidence including
  CDS/UTR quality, homology, and classify/reduce into a biologically
  informative transcriptome of primary and alternate transcripts.

evgr_tsasubmit.sbt
Submit-block ::= {
  contact { contact { name name { ... } }
  ..
  Seqdesc ::= user { type str "DBLink",
  data { { label str "BioProject", num 1, data strs { "PRJNA99999" } } } << need PRJNA here..
  }
  

=cut

sub tsa_infotemplates  # make if not found as files; edit w/ sra_result info, other
{
  my($trpath, $trname, $sradata)= @_; ## add sraresult hash info
  my($methtxt,$desctxt,$subtxt)=("") x 3;
  my($tsamethf,$tsadescf,$tsasubf)=("") x 3;
  my $update=0; my $nupdate=0;
  my $NEEDUPDATE = ! $skipTSAparts;
  
  # FIXME: must have \tabs in StructuredComment
  # FIXME paths: trpath may be submitset/ for 2nd runs..
  my $tpath=$trpath;
  my($okd,@files)= getFileset($tpath,'cmt|sbt');  # check for trname in file set
  unless(@files) {
    $tpath =~ s,\w+$,submitset,; 
    $tpath= $trpath."/submitset" unless(-d $tpath);
    ($okd,@files)= getFileset($tpath,'cmt|sbt');  
  }
  
  # ($tsamethf)= grep/$trname.tsamethods.cmt/,@files; $update=0;
  # unless($tsamethf) { ($tsamethf)= grep/evgr_tsamethods.cmt/,@files; $update=$NEEDUPDATE; }
  
  $tsamethf= "$tpath/$trname.tsamethods.cmt";  $update=0;
  #? unless( -s $tsamethf) { ($tsamethf)= grep /tsamethods.cmt/, @files; $update=$NEEDUPDATE; }
  unless( -s $tsamethf) {  $tsamethf= "$tpath/evgr_tsamethods.cmt"; $update=$NEEDUPDATE; }
  if( -s $tsamethf) {
    open(F, $tsamethf); $methtxt= join "", <F>; close(F);
  } else {
    $update=$NEEDUPDATE; $methtxt=""; # TEMPlATE
    my $Illumina= $sradata->{"Instrument"} || "Illumina"; # from sradata
    my $asmsoft= $sradata->{'Assemblers'}||""; $asmsoft="; $asmsoft" if($asmsoft); 
    my $asmname="$trname.evigene"; # from $trname
    my $egvers= 'v'.VERSION;
    $methtxt=<<"EOT";
StructuredCommentPrefix\t##Assembly-Data-START##
Assembly Method\tEvidentialGene $egvers$asmsoft
Assembly Name\t$asmname
Sequencing Technology\t$Illumina
EOT
  }
  if($update) {
    $tsamethf= "$tpath/$trname.tsamethods.cmt"; $nupdate++;
    open(O,'>',$tsamethf); print O $methtxt,"\n"; close(O);
	  push @submitset, $tsamethf;  
  }
  
  $tsadescf= "$tpath/$trname.tsadesc.cmt"; $update=0;
  #? unless( -s $tsadescf) { ($tsadescf)= grep /tsadesc.cmt/, @files; $update=$NEEDUPDATE; }
  unless( -s $tsadescf) {  $tsadescf= "$tpath/evgr_tsadesc.cmt"; $update=$NEEDUPDATE; }
  if( -s $tsadescf) {
    open(F, $tsadescf); $desctxt= join "", <F>; close(F);
  } else {
    $update=$NEEDUPDATE; $desctxt=""; # TEMPlATE
    # m2t.Assemblers=Velvet/Oases v1.2.07/o0.2.08 (2012.10); SOAPdenovo-Trans v2011.12.22; Trinity v2012.03.17 ; Cufflinks v1.0.3 (2012.07)
    my $asmsoft= $sradata->{'Assemblers'} || "various RNA assemblers"; #?
    if($asmsoft) { $asmsoft =~ s/\s+v\d[^;,]+//g; $asmsoft =~ s/;/,/g; } ## drop version info, kept above..
    my $datasize= $sradata->{"Total Size, Mb"} || ""; ## "Total Size, Mb"
    my $nspots= $sradata->{"Total Spots"} || "";
    my @ds=split";",$datasize;  my @ns=split";",$nspots; 
    if(@ds>1) { my $ds=0; my $ns=0; 
      for(my $i=0; $i<@ds; $i++) { my $n;
        ($n)= $ds[$i] =~ m/(\d+)/; $ds+=$n; 
        ($n)= $ns[$i] =~ m/(\d+)/; $ns+=$n; 
        } 
      $datasize=$ds; $nspots= $ns; }
    # calc mb from nspots x 100 x 2 ?
    if($nspots>1000 and not $datasize) { $datasize= int($nspots * 200/1000000); } # ok or not?
    if($datasize > 0) {  $nspots||=""; $datasize=", $datasize Mb in $nspots read-pairs,"; }
    else { $datasize=""; }
    my $asmcount= $sradata->{'Total Assemblies'} || "many"; #?
    # m2t.Total Assemblies=5583471

    $desctxt=<<"EOT";
RNA-Seq data$datasize of $organism are assembled with
$asmsoft, using multiple options. 
EvidentialGene tr2aacds pipeline is used to process the $asmcount
resulting assemblies by coding sequences, translated proteins, and gene
evidence, then classify/reduce to a biologically informative
transcriptome of primary and alternate transcripts.
EOT
  }
  if($update) {
    $tsadescf= "$tpath/$trname.tsadesc.cmt"; $nupdate++;
    open(O,'>',$tsadescf); print O $desctxt,"\n"; close(O);
	  push @submitset, $tsadescf;  
  }

  # check for $ENV{tsasubmit_sbt} ?
  $tsasubf= "$tpath/$trname.tsasubmit.sbt"; $update=0;
  unless( -s $tsasubf) { $tsasubf= $ENV{tsasubmit_sbt} || "$tpath/evgr_tsasubmit.sbt";  $update=$NEEDUPDATE; }
  if( -s $tsasubf) {
    open(F, $tsasubf); $subtxt= join "", <F>; close(F);
  } else {
    $update=$NEEDUPDATE; $subtxt=""; # TEMPlATE cant really guess this. but empty sbt isnt useful
  }
	## need "BioProject",  "PRJNA99999" from somewhere ..
  ##   data { { label str "BioProject", num 1, data strs { "PRJNA99999" } } }
  if(my $bioprj= $sradata->{"BioProject"}) { $update=$NEEDUPDATE if($subtxt=~s/PRJNA99999/$bioprj/); }
  
  if($update and $subtxt =~ /\w/) {
    $tsasubf= "$tpath/$trname.tsasubmit.sbt"; $nupdate++;
    open(O,'>',$tsasubf); print O $subtxt,"\n"; close(O);
	  push @submitset, $tsasubf;  
  }

  ## edit TSADESC if($nupdate or paths not in $TSADESC)
  ## my $TSADESC="-w evgr_tsamethods.cmt -Y evgr_tsadesc.cmt -t evgr_tsasubmit.sbt"; # FIXME: where from?
  unless($TSADESC =~ m/$tsamethf/ and $TSADESC =~ m/$tsadescf/ and $TSADESC =~ m/$tsasubf/ ) {
    #BAD# $TSADESC="-w $tsamethf -Y $tsadescf -t $tsasubf"; 
    $TSADESC.=" -w $tsamethf" unless($TSADESC =~ s/\-w \S+/\-w $tsamethf/);
    $TSADESC.=" -Y $tsadescf" unless($TSADESC =~ s/\-Y \S+/\-Y $tsadescf/);
    $TSADESC.=" -t $tsasubf" unless($TSADESC =~ s/\-t \S+/\-t $tsasubf/);
    #see.above# loggit(0,"info updated $nupdate TSADESC=",$TSADESC);
  }
  
  return($nupdate,$tsamethf,$tsadescf,$tsasubf);
}


sub tsa_tbl2asn
{
  my($cdnaseq,$cdnatbl,$organism,$sraids)=@_;
  return unless(-s $cdnaseq);
 
  ## FIXME 2015.01, missing .pep/.pep.report for NCPU case .. need to split those also
   
  ##.. this file part input works..
  # pt=litova1all3.mrnatop2000
  # $ncbin/tbl2asn -i tsa.$pt.fsa -f tsa.$pt.tbl -Z discrep.$pt.log \
  # -a r10k -l paired-ends -Vtb -Mt -XE \
  # -w evgr_tsamethods.cmt -Y evgr_tsadesc.cmt -t evgr_tsasubmit.sbt \
  # -j "[moltype=mRNA] [tech=TSA] $org $sra"
  ##  -r  Path for Results [String]  Optional

## fix2: test cdnaseq.gz cdnatbl.gz - gunzip, dont think tbl2asn can read gz
## FIXME FIXME FIXME damn paths again
## .info has basedir ./ for info files; but above moves them to ./submitset/
## should resave info after tidyup w/ new paths..
## m2t.TSADESC=-w ./ztick4cvel_pt3.tsamethods.cmt -Y ./ztick4cvel_pt3.tsadesc.cmt -t ./ztick4cvel_pt3.tsasubmit.sbt
  
  # my $wdir= dirname($cdnaseq);
  # my $cdnabase= basename($cdnaseq);
  
  my $spldir= makename($cdnaseq,"_tsasubmit/",'cdna|fasta|fsa|fa');  # ok for _split dir
  # ^^ UNUSED NOW? 
  
  ## -Vb = gen genbank gbf, not needed, option; -Vtb == genbank + tsa
  my $tsaopts="-a r10k -l paired-ends -Vt -Mt -XE"; # read from config.; use TSADESC ??
  my $tsadesc= $TSADESC; # "-w evgr_tsamethods.cmt -Y evgr_tsadesc.cmt -t evgr_tsasubmit.sbt"; # FIXME: where from?
    ## ^^ need config for these; generate some of this?
    ## FIXME: tbl2asn dies silently if missing files -w .. -Y .. or -t ..
    ## UGH: qsub gave this via env sra='xxx':
    ##   -j '[moltype=mRNA] [tech=TSA] [organism=Pogonus chalceus] [SRA=\'SRR424344;SRR424342;SRR424340\']'
    ## egmrna2tsa.19433.err sh: -c: line 0: unexpected EOF while looking for matching `''
#er2g: forkCMD= /home/ux455375/bio/ncbi2227/bin/tbl2asn -a r10k -l paired-ends -Vt -Mt -XE \
#  -w ./pogonus1all3.tsamethods.cmt -Y ./pogonus1all3.tsadesc.cmt -t ./pogonus1all3.tsasubmit.sbt \
# -j '[moltype=mRNA] [tech=TSA] [organism=Pogonus chalceus] [SRA=\'SRR424344;SRR424342;SRR424340\']' \
# -i ./okayset/pogonus1all3.mrna_tsasubmit/pogonus1all3.mrna.split1.fsa \
# -f ./okayset/pogonus1all3.mrna_tsasubmit/pogonus1all3.mrna.split1.tbl \
# -Z ./okayset/pogonus1all3.mrna_tsasubmit/pogonus1all3.mrna.split1.discrep

  # map{ s/^['"]+//; s/['"]+$//; } ($organism, $sraids); ## $sraids=~s/ +;/;/g; $sraids=~s/ +/;/g; ??
  map{ s/^\W+//; s/\W+$//; } ($organism, $sraids); ## $sraids=~s/ +;/;/g; $sraids=~s/ +/;/g; ??
  my $tsaqual="-j \'[moltype=mRNA] [tech=TSA] [organism=$organism] [SRA=$sraids]\'";
  my $cmd0="$APPtbl2asn $tsaopts $tsadesc $tsaqual "; # add in loop: -i $pt.fsa -f $pt.tbl -Z discrep.$pt.log
  
  my ($cmddone, $nout, $sqnout)=(0) x 10;
  my $runok= ($APPtbl2asn =~ /MISSING/) ? 0 : 1;

  #?? check for tbl2asn file set?  
  ## FIXME: may be in "submitset" subdir : use instead makename($cdnaseq,"_tsasubmit/")
   
  $sqnout= makename($cdnaseq,".sqn",'cdna|fasta|fsa|fa');   
  return(1,"",$sqnout) if( -s $sqnout ); ## all sqnout have full path ?
  if( -d $spldir ) { ## ! -s $sqnout and 
    opendir(D,$spldir); 
    my @sqnout= map{ "$spldir/$_" } grep /\.sqn/, readdir(D); 
    closedir(D);
    $sqnout= join ", ", @sqnout; $nout= @sqnout;
    return($nout,$spldir,$sqnout) if($nout>0); 
  }
  
  $sqnout="";
  #only>1c# mkdir($spldir); # dryrun? only ncpu or use for all cases?
  if($runok and $NCPU > 1) {
    my $ccount= facount($cdnaseq); # use this, not fasize; ERR if ccount<??
    if($ccount >= 50*$NCPU) {
      mkdir($spldir); # dryrun? only ncpu or use for all cases?
      ($nout,$sqnout)= tbl2asn_ncpu($NCPU,$cmd0,$ccount,$spldir,$cdnaseq,$cdnatbl);
      $cmddone=1;
    	#?? push @submitset, (split /,\s*/, $sqnout); # leave in name_tsasubmit ; 
    	# move all of name_tsasubmit to submitset/ ??
    }
  } 
  
    # put into $spldir ?
  if($runok and ! $cmddone) {
    $spldir="./"; # cwd?
    my $dlog1=$cdnatbl; $dlog1 =~ s/\.\w+$//; $dlog1.=".discrep";
    my $cmd1= $cmd0 . " -i $cdnaseq -f $cdnatbl -Z $dlog1"; # want -r $spldir or not ?
    my $err= runcmd($cmd1);    
    my $sqnt= makename($cdnaseq,".sqn",'cdna|fasta|fsa|fa'); # add $spldir ??
    if( -s $sqnt ) { 
      my $f; $nout=1; $sqnout=$sqnt; 
      ## rename errorsummary so not overwrittne by others
      my $esum0= dirname($sqnout)."/errorsummary.val";
      my $esumf= makename($sqnout,".sumval"); # or cdnaseq.errorsummary.val ?
      rename($esum0, $esumf); 
      push @submitset, $sqnout, $dlog1, $esumf; 
      for my $suf (qw(val gbf fixedproducts pep.report)) { 
        ($f=$sqnt) =~ s/\.sqn/.$suf/; push @submitset,$f if(-f $f);
      }
      # t2a outputs include discrep, xxx.sqn, xxx.val, errorsum.val ..
      # FIXME: submitset tbl2asn: name.val,.sqn,.discrep,.fixedproducts errorsummary.val
    }
  }
  
  return($nout,$spldir,$sqnout); # list all?
}


sub tbl2asn_ncpu
{
  my($npart,$cmd0,$ccount,$spldir,$cdnaseq,$cdnatbl)=@_;
  # $NCPU,$cmd0,$ccount,$spldir,$cdnaseq,$cdnatbl
  
  # my $ccount= facount($cdnaseq); # use this, not fasize; ERR if ccount<??
  my $splcount= int(0.99 + $ccount/$npart);
  # my $spldir= makename($cdnaseq,"_split/");  # use _tsasubmit/ instead?
  # mkdir($spldir); # dryrun?
  
  ## NOTE: fasplit adds spldir to set path
  my @splset= fasplitcount( $cdnaseq, $spldir, $npart, $splcount,"fsa"); 
  my @tblset= fasplitcount( $cdnatbl, $spldir, $npart, $splcount,"tbl");  # need to split into same parts.
  
  ## FIXME 2015.01, missing .pep/.pep.report for NCPU case .. need to split those also
  # * PROBLEM for pep split, out of order ids ? got 2 in wrong split set..
  # .. fsa,tbl above are same id-order as tbl is made from cdna/fsa ..
  # .. need idlist from fsa, or each part, and new fasplitordered(...) to manage those.
  
  my @pepset=(); 
  # my $cdnapeps=$cdnaseq; $cdnapeps =~ s/\.\w+$/.pep/;
  my $cdnapeps= makename($cdnaseq,".pep"); 
  if(-s $cdnapeps) { 
    if(1) { # split by idlist..
      my $ns= @splset;
      for(my $ip=0; $ip < $ns; $ip++) {
        my $cdna1= $splset[$ip];
        my $pep1= makename($cdna1,".pep"); 
        my $idh = faidlist($cdna1,{},"ashash");

        # Ah, damn diff ids for .fsa > .pep Id000t1 > gnl|Evigene|Id000p1  
        # FIXME in %$idh ??
			  # $protid= $GDB_PREFIX.$protid if($GDB_PREFIX); #  protein_id  gnl|CacaoGD|Thecc1EG016762p1
        for my $m (keys %$idh) { my $n=$m; $n=~s/t(\d+)$/p$1/;
          $n=$GDB_PREFIX.$n if($GDB_PREFIX); $idh->{$n}=1; 
        }

        $pep1= faextract( $cdnapeps, $pep1, $idh, ); # FAILed empty peps
        push @pepset, $pep1;
      }
    } else {
      @pepset= fasplitcount( $cdnapeps, $spldir, $npart, $splcount,"pep"); 
    }
  }
  
  my $npartgot= @splset;
  my $icpu= 0;  
  # NO: chdir($spldir); # bad!, need files in starting path ... use -i spldir/in -r spldir/
  for(my $ip=0; $ip< $npartgot; $ip++) {
    my $cdna1= $splset[$ip];
    my $tbl1 = $tblset[$ip];
    ## dont need to add, tbl2asn looks for it# my $pep1 = $pepset[$ip];
    (my $dlog1=$tbl1) =~ s/\.\w+$/.discrep/; # $dlog1= makename($cdna1,".discrep");
    my $cmd1= $cmd0 . " -i $cdna1 -f $tbl1 -Z $dlog1"; 
      ## DONT need -r spldir ? w/o put parts in same dir as -i, -f
    my $pid= forkcmd($cmd1);    
    if(++$icpu > $npartgot) { while (wait() != -1) { }; $icpu= 0; }
  }
  while (wait() != -1) { };
  
  ## One part problem, errorsummary.val parts overwrite .. rebuild from $name.val ? or
  # my @valout= # my @gbfout= 
  # my @sqnout= map{ my $f=$_; $f=~s/\.fsa/.sqn/; $f; } @splset;
  
  ## deal with $name.split*.fixedproducts also?
  ## deal with $name.split*.discrep also?
  
  ## remake errorsummary.val from parts  .. maybe for 1cpu also? count ERROR: and loggit ?
  my %et;
  my $esumf= makename($cdnaseq,".sumval"); # or cdnaseq.errorsummary.val ?
  foreach my $sf (@splset) { 
    (my $valf=$sf) =~ s/\.fsa/.val/;
    if(open(VF,$valf)) { 
      while(<VF>) { my($err,$typ)= m/^(\S+)[\w\s]+.([\w\.]+)/; $et{$err}{$typ}++ if($typ); } close(VF); 
    }
  }
  open(my $esumh,'>',$esumf);
  foreach my $e (sort keys %et) { my @t=sort keys %{$et{$e}}; 
    foreach my $t (@t) { my $c=$et{$e}{$t}; printf $esumh "%6d %s %s\n",$c,$e,$t; } 
  } close($esumh); 
	#?? push @submitset, $esumf; 

  ## check existance -s sqn     # or readdir(D) as above  
  my @sqnout= grep { -s $_ } map{ my $f=$_; $f=~s/\.fsa/.sqn/; $f; } @splset; 
  my $sqnout= join ", ", @sqnout;
  my $nsqn= @sqnout;
  return($nsqn,$sqnout); # list all?
}

sub END_here {}


__END__

=item tbl2asn log

#er2g: app=tbl2asn, path=/home/gilbertd/bin/tbl2asn
#er2g: evigeneapp=prot/traa2cds.pl, path=/bio/bio-grid/mb/evigene/scripts/prot/traa2cds.pl
#er2g: get_evgtrset= ./okayset/allstrimpt1.mrna . allstrimpt1
#er2g: trclass2maintab primary n=  allntr=  allstrimpt1.pubids
#er2g: vectors found in ntr= 962 ./okayset/allstrimpt1.mrna.vector.tab
#er2g: DONE output ntr=0, ncds=0 in files allstrimpt1.mainalt.tab, allstrimpt1.pubids, ./okayset/allstrimpt1.mrna.fsa, ./okayset/allstrimpt1.mrna.tbl, ./okayset/allstrimpt1.mrna.ann.txt
#er2g: forkCMD= /home/gilbertd/bin/tbl2asn -a r10k -l paired-ends -Vt -Mt -XE \
 -w evgr_tsamethods.cmt -Y evgr_tsadesc.cmt -t evgr_tsasubmit.sbt \
 -j '[moltype=mRNA] [tech=TSA] [organism=Penaeus monodon] [SRA=SRX110652; SRX110651; SRX110649]' \
 -i ./okayset/allstrimpt1.mrna_tsasubmit/allstrimpt1.mrna.split1.fsa \
 -f ./okayset/allstrimpt1.mrna_tsasubmit/allstrimpt1.mrna.split1.tbl \
 -Z ./okayset/allstrimpt1.mrna_tsasubmit/allstrimpt1.mrna.split1.discrep
#er2g: forkCMD= /home/gilbertd/bin/tbl2asn -a r10k -l paired-ends -Vt -Mt -XE -w evgr_tsamethods.cmt -Y evgr_tsadesc.cmt -t evgr_tsasubmit.sbt -j '[moltype=mRNA] [tech=TSA] [organism=Penaeus monodon] [SRA=SRX110652; SRX110651; SRX110649]'  -i ./okayset/allstrimpt1.mrna_tsasubmit/allstrimpt1.mrna.split2.fsa -f ./okayset/allstrimpt1.mrna_tsasubmit/allstrimpt1.mrna.split2.tbl -Z ./okayset/allstrimpt1.mrna_tsasubmit/allstrimpt1.mrna.split2.discrep
#er2g: forkCMD= /home/gilbertd/bin/tbl2asn -a r10k -l paired-ends -Vt -Mt -XE -w evgr_tsamethods.cmt -Y evgr_tsadesc.cmt -t evgr_tsasubmit.sbt -j '[moltype=mRNA] [tech=TSA] [organism=Penaeus monodon] [SRA=SRX110652; SRX110651; SRX110649]'  -i ./okayset/allstrimpt1.mrna_tsasubmit/allstrimpt1.mrna.split3.fsa -f ./okayset/allstrimpt1.mrna_tsasubmit/allstrimpt1.mrna.split3.tbl -Z ./okayset/allstrimpt1.mrna_tsasubmit/allstrimpt1.mrna.split3.discrep
#er2g: DONE tsa_tbl2asn nparts=3, submitset=./okayset/allstrimpt1.mrna_tsasubmit/  ./okayset/allstrimpt1.mrna_tsasubmit/allstrimpt1.mrna.split1.sqn, ./okayset/allstrimpt1.mrna_tsasubmit/allstrimpt1.mrna.split2.sqn, ./okayset/allstrimpt1.mrna_tsasubmit/allstrimpt1.mrna.split3.sqn

remake errorsummary.val

cat okayset/allstrimpt1.mrna_tsasubmit/allstrimpt1.mrna.split?.val | perl -ne\
'($err,$typ)= m/^(\S+)[\w\s]+.([\w\.]+)/; $et{$err}{$typ}++; END{ foreach $e (sort keys %et) { @t=sort keys %{$et{$e}}; foreach $t (@t) { $c=$et{$e}{$t}; printf "%6d %s %s\n",$c,$e,$t; } } }' | head

     1 ERROR: SEQ_FEAT.NoStop
     2 WARNING: SEQ_FEAT.CDShasTooManyXs
     2 WARNING: SEQ_FEAT.EcNumberProblem
 33777 WARNING: SEQ_FEAT.PartialProblem
     8 WARNING: SEQ_INST.HighNContentStretch


=item tbl2asn work

# ## tbl2asn outputs:  .sqn == asn1, can it be catted ??? NO.. need submit each part file?
# ##  .val == errors list, can be catted
# ##  .gbf if made can be catted
# # my $cmd;
# # $cmd= "cat ".join(' ',@aaset)." > $aaseq";   runcmd($cmd);
# # $cmd= "cat ".join(' ',@cdsset)." > $cdsseq"; runcmd($cmd);
## keep as part files for tsa submit for now
#   if($DEBUG||$dryrun) {
#     push @erasefiles, @splset, @aaset, @cdsset; # which?
#   } else {
#     foreach my $fn (@splset, @aaset, @cdsset) {  unlink($fn) if(-f $fn); } 
#     rmdir($spldir); 
#   }
  
=cut


## in cdna_evigenesub.pm
# sub parse_genenames {}
# sub parse_evgheader {}
# sub trprocess_OLD {}
# sub putAnnoSeq_OLD { } # rename: putseq > putAnnoSeq ; merge w/ update_mrna_fileset

## moved to rnaseq/asmrna_trimvec.pl
# sub vecscreen { }
# sub vecscreen_ncpu {}

## move to cdna_evigenesub.pm
# sub revcomp {}
# sub facount {}
# sub fasize {}
# sub fasplit {}
# sub fasplitcount {}
# sub findapp {}
# sub findevigeneapp {}
# sub runcmd {}
# sub forkcmd {}
# sub makename {}


 
 
=item tsasub11 cacao

  evigene=/Users/gilbertd/Desktop/dspp-work/genes2/evigene
  ncbin=/Users/gilbertd/Desktop/dspp-work/genomesoft//ncbi/bin
  nwbsra='[SRA=SRR531454; SRR531455; SRA058778; SRA058779]'
  vel3sra='[SRA=SRA058777; SRA058780; SRA058781; SRA058782; SRR531454; SRR531455; SRA058778; SRA058779]'
  rnasra='[SRA=SRA058777; SRA058780; SRA058781; SRA058782]'
    ## redo asmrna2ncbitsa.pl -GAPSOK
    
  cd tsasub11
  
  for pt in cacao[345]{cuf[28],nwb,vel,sop,tri}; do {
  
  if [ -f $pt/TCM01.tsa_rasm.$pt.tbl ]; then continue; fi
  echo asmrna2ncbitsa $pt/TCM01.tsa_rasm.$pt
  
  $evigene/scripts/rnaseq/asmrna2ncbitsa.pl -GAPSOK -idpre Thecc1ER_ \
  -cdna ../tr5parts/pub3ig.$pt.tab4g.tr.gz -vec ../tr5parts/pub3ig.trasm.tab4g.vector.tab \
  -geneinfo ../tr5parts/pub3ig.trasm.tab4g.geneinfo1.tab  -log tr4g.$pt.log \
  -out $pt/TCM01.tsa_rasm.$pt.fsa -tbl $pt/TCM01.tsa_rasm.$pt.tbl
  
  } done
  
  
    # was -a r10u; unknown gap len; YES, r10k got rid of warning
  
  for pt in cacao[345]{cuf[28],nwb,vel,sop,tri}; do {
   if [ $pt == cacao3nwb ]; then sra=$nwbsra; 
   elif [ $pt == cacao3vel ]; then sra=$vel3sra; 
   else sra=$rnasra; fi
   
   if [ -f $pt/TCM01.tsa_rasm.$pt.sqn ]; then continue; fi
  
   echo $pt : $sra
   
   $ncbin/tbl2asn -p $pt/ -Z $pt/$pt-discrep.log  -w tsamethods.$pt.cmt \
   -Y tsadesc.$pt.cmt -t cacao3i_tsasubmit.sbt \
   -a r10k  -l paired-ends -Vtb -Mt -XE \
   -j "$sra [bioproject=PRJNA51633] [moltype=mRNA] [tech=TSA] [organism=Theobroma cacao]"
  
  } done
  
=cut

=item run_evgmrna2tsa.sh

  #! /bin/bash
  ### env idprefix=MysppEGm trclass=myspecies.trclass datad=`pwd`  qsub -q shared run_evgmrna2tsa.sh
  #PBS -N evgmrna2tsa 
  #PBS -l nodes=1:ppn=8,walltime=5:55:00
  #PBS -o egmrna2tsa.$$.out
  #PBS -e egmrna2tsa.$$.err
  #PBS -V
  
  ncpu=8; # most 100K trsets run on 8cpu in 10m-20min; 300K set can take 1hr.

  evigene=$HOME/bio/evigene/scripts
  export vecscreen=$HOME/bio/ncbic11/bin/vecscreen
  export tbl2asn=$HOME/bio/ncbi2227/bin/tbl2asn
  
  ## FIXME: -TSADESC=tbl2asn flags for path/to/tsa.cmt files
  ## default $TSADESC="-w evgr_tsamethods.cmt -Y evgr_tsadesc.cmt -t evgr_tsasubmit.sbt";
  
  opts="-debug -runtbl2asn -NCPU $ncpu"
  if [ "X" = "X$datad" ]; then echo "missing env datad=path/to/data"; exit -1; fi
  if [ "X" = "X$trclass" ]; then "echo env trclass=path/to/name.trclass"; exit -1; fi
  ## .. these are now read via sra_result.csv, species => idprefix
  if [ "X" != "X$species" ]; then opts="$opts -species=\'$species\'"; fi
  if [ "X" != "X$sra" ]; then opts="$opts -sraids=\'$sra\'"; fi
  if [ "X" != "X$idprefix" ]; then opts="$opts -idprefix $idprefix"; fi
  
  cd $datad/
  ## add -outdir opt or : mkdir tsasubmit; cd tsasubmit
  ##FIXME: template files ; need these  or generate defaults?
  if [ ! -f evgr_tsasubmit.sbt ]; then
    if [ -d ../tsasubmit ]; then cp ../tsasubmit/*.{cmt,sbt} ./; fi
  fi
    
  echo $evigene/evgmrna2tsa.pl  $opts -log -class $trclass
  $evigene/evgmrna2tsa.pl  $opts -log -class $trclass

=cut

