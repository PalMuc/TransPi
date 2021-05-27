#!/usr/bin/env perl
# evgs/genes/trclass2pubset.pl 

=item trclass2pubset.pl

 revision of parts of  evgmrna2tsa2.pl
 following genes/altbest2pubset.pl (pubset from annot GFF)
 roughly same as MAIN_pubsetonly() of evgmrna2tsa2
 
 updated parts:
   trclass2maintab() =  public class/id table from trclass, okayset of tr2aacds
   
=item usage (maybe)

 set pt=pig4321ew; 
 $evigene/scripts/genes/trclass2pubset.pl -debug -outbase pig4ewc \
  -keepdrop $pt.keepdrop.tab -keepoldids $pt.keepids -names $pt.names \
  -mrna tmpfiles/$pt.mrna -trclass $pt.trclass

 .. this is fairly messy set of opts, need all?
 .. expect/require? inputs of both -mrna (with okayset oids) and -trclass table
    name.mrna should also have name.aa, cds seq files associated
    -mrna could be from okayset/, or from vec trimset result
 .. want output of all publicset/, pubids, ann.txt
 .. -outbase should be option, default to publicset/$trname
 .. fix -keepdrop <> -keepoldids ?
     -keepoldids is unusual option, but want for now
     -keepdrop table may be common input, for culled and dropped okayset
   
=cut

use FindBin;
use lib ("$FindBin::Bin","$FindBin::Bin/.."); # assume evigene/scripts/ layout: main, prot/, rnaseq/, ...
use strict;
use warnings;
use Getopt::Long;
use cdna_evigenesub; 
use evigene_pubsets; # now has some of below subs

use constant VERSION => '2018.06.18'; #  from evgmrna2tsa2.pl & altbest2pubset

our $EVIGENES= "$FindBin::Bin/.."; # ="$FindBin::Bin";  # WRONG place ; this is in evgs/genes/ subdir ..
our $DEBUG= $ENV{debug}|| 0;
our $EGAPP='class2pub';  # or okay2pub ?
our $IDPREFIX= $ENV{idprefix} || 'EVGm'; 
our $ORGANISM= $ENV{ORGANISM} || ""; # no default
my  $NEWIDpre= $ENV{idtag}||'n'; #??
my  $SRC= $ENV{source}||"evm"; # fixme: for gff output, source.class col2 src_(main/alt/other);
our $SHOWDROPS=1; #default on? only for mainalt.tab outputs
our $SIZESORT=1; 
our $FAHDRisValid=1; # export from pubset.pm for annots from mrna hdr
our ($pubid_format,$altid_format,$GBPROID);
our $pubidnum_start=0;
our $preserveOldIds=undef; # change to ok pattern, IDPREOK
our $SORT_PUBSET= 0;
our $DATE;  
our $RNATYPES; # our $RNATYPES='mRNA|ncRNA'; # in evigene_pubsets.pm; transcribed_RNA mrnatypes exontypes 
our ($KEEPDROPX,%keepdropxtab); # evigene_pubsets.pm
my $ONLY_PUBIDS= 0; my $CLASSbySOURCE= 0;
my $CULLOTHERS= 0; # how to treat unclassed input genes

# main opt:
my $AAMIN_NOCLASS=$ENV{aaminnoclass}||60; # asmrna_altreclass -noclasscut=$AAMIN_NOCLASS; drop noclass (needs user opt), but rescued by aaref
my $NOALTDROPS=0; # turn on if have keepdrop table?

my ($trclass,$output,$logfile,$keepdropin);  
my ($trpath,$trname, $sradatah, %settings);
my ($oname, @genes, @seqset, $genenames);

our @evgdirs = qw(okayset dropset inputset trimset tmpfiles erasefiles publicset submitset);
our (@okayset,@dropset,@inputset,@tmpfiles,@erasefiles,@publicset,@submitset); # tidyup file sets
my $tidyup=1; my $NCPU=1;

my $optok= GetOptions(
  "trclass|class|input=s", \$trclass, # many gff inputs : use perl input pipe while(<>) ?
  "mrna|cdna|sequences=s", \@seqset, # many inputs? 
  "oname|outbase|output=s", \$oname,  # FIXME: idpre option  overwritten by spppref
  "idprefix=s", \$IDPREFIX,  
  "tagidnew=s", \$NEWIDpre,   
  "source=s", \$SRC,   
  "names=s", \$genenames,   
  "keepoldids|preserveOldIds=s", \$preserveOldIds,  # allow -idpre=XXX -keepold w/o 2nd param
  "keepdropin=s", \$keepdropin,   # other opt key? reclassin?
  "logfile:s", \$logfile,
  "dropshow!", \$SHOWDROPS,  
  "sizesort!", \$SIZESORT,  # default:on ; may be required to get proper ids w/ preserveOld
  "pubsortseq!", \$SORT_PUBSET,  # default:off?
  "CULLOTHERS!", \$CULLOTHERS,  # or otherclass=cull|skip|remainder name?
  "onlypubids|ONLY_PUBIDS!", \$ONLY_PUBIDS,  #CLASSbySOURCE also : same opt?
  #gff: "classbysource!", \$CLASSbySOURCE,  # default:off?
  # "dryrun|n!", \$dryrun, 
  "NCPU=s", \$NCPU, # not used here, leave in opts   
  "tidyup!", \$tidyup, # default on ?
  "debug!", \$DEBUG, 
  );

#? push @genes, grep(/gff/, @ARGV); # remainder gff ? trclass here?

die "EvidentialGene trclass2pubset -trclass myspp.trclass [-idprefix $IDPREFIX ] 
  makes tables of public ids and main-alt linkage, from results of tr2aacds
  opts: -idprefix Thecc1EG  -mrna myspp.mrna  -names mrna.names 
     -keepdrop keep_drop_ids.table -preserveOldIds=old.pubids
     -nosizesort -[no]pubsortseq -debug
  version ", VERSION, "\n"
  unless($optok and $trclass); # (@genes or @seqset) 

$NEWIDpre="" unless($NEWIDpre=~/^[A-Za-z]/); # allow same old id?

my $IDPREOK=$IDPREFIX;
# if(defined $preserveOldIds) {
#   if($preserveOldIds and $preserveOldIds=~/^[a-zA-Z]/) { $IDPREOK= $preserveOldIds; }
#   else { $preserveOldIds=$IDPREOK=$IDPREFIX; }
# } else {
#   $preserveOldIds=0;
# }  

# data globals
my($mainindex,$ntr)= (0) x 9;
my(%annotes, %main, %mainsize,%alt,%oclass,%altsize,%altdrops,%didaltdrops,%balt,%drop);  # globals of trclass2maintab
#unused here: my(%aaqual, %piad, %newmain, %altdefer,%maindrops);

my($nkeepdrop, $keepdroph, )= (0);
my($nmapqual, $mapqualh, $alntabh)= (0); #? globals from map align.tab, if exists
# %dropids, %dropid not used yet; change to culldropkeep table: id, actclass, actinfo, oid?, other..
# orig dropids,dropid is  short list of gff.source tags, e.g. 'cull' in source 'zf17evgcull'

openloggit($logfile,$trclass);  
loggit(1, "EvidentialGene trclass2pubset VERSION",VERSION);

MAIN();

sub MAIN
{
  ## default output to publicset/oname
  ## maybe default input seqs from okayset/ ; see mrna2tsa:get_evgtrset()
  unless($oname){ $oname=$trclass || $seqset[0]; $oname=~s/\.\w+$//; }

  $DATE=`date '+%Y%m%d'`; chomp($DATE);
  loggit(0, "BEGIN $0  input=",$oname,"date=",$DATE);
	#? do_settings("restore",$trclass);  
	
  my($upokgenes, $upokids, $upstatus)= (0) x 9; 
  my($pubids,$maltab)= map{ "$oname.$_" } qw(pubids mainalt.tab);
  # my $generef= (@genes > 0)? \@genes : undef;
  my $mrnaseq= (@seqset>0)? $seqset[0] : ""; # many?? all input .mrna .cds and/or .aa

  ($mrnaseq,$trpath,$trname)= get_evgtrset($trclass,$mrnaseq,"publicset"); # subs.pm; cdnaseq => mrnaseq here
  loggit(0, "get_evgtrset=",$mrnaseq,$trpath,$trname);  
  loggit(LOG_DIE, "Missing -mrna",$mrnaseq) unless($mrnaseq and -s $mrnaseq);
  
  ($nkeepdrop,$keepdroph)= ($keepdropin) ? readKeepDrop($keepdropin) : (0); # fill  global %keepdrop
  $NOALTDROPS=1 if($nkeepdrop);

  ## FIXME19: pubset wants trname.align.tab , gmap now makes trname-chrname.align.tab ..

  # (my $mapqual=$trclass)=~ s/\.\w+$/.align.tab/;
  # if(! -f $mapqual and -d "geneval") { $mapqual="geneval/$mapqual"; }
  my($mok,$mapqual)= getMapqual($trname);
  if($mok) {
    ($nmapqual, $mapqualh, $alntabh)= readAlignTab($mapqual);
  }
  
  # Separate steps here, call w/ input pubid tab to make pubgff, pubseq sets
  unless( -s $pubids ) {
    ($upokids)= makePubIdtabFromTrclass($pubids,$maltab,$trclass); # read attr from mrnaseq?
	  $upstatus++ if($upokids>0);
	  #see makePubIdtab:  altreclass_block($trclass,$pubids); # reclassifies alts in pubids table, drops some
  }
  
  if($ONLY_PUBIDS) { # want this opt?
    loggit(0, "done -onlypubids makePubIdtab n=",$upokids);
    #x exit(0); # go to tidyup below
  } else {
     
    ## change this mrnaseq to array= grep 'mrna|cdna|xxx' @seqset
    ($upokgenes)= makePublicset($pubids, $oname, $mrnaseq); # , $generef);
    $upstatus++ if($upokgenes>0);
  }
  
	#? do_settings("log|save",$trclass||$mrnaseq,); # or after last call?  ("log|restore|save");

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
}
  
# ---------------------------

## readKeepDrop moved to evigene_pubsets.pm
# my($nkeepdrop,$keepdrop_idhash)= readKeepDrop($keepdropin) if($keepdropin); # fill  global %keepdrop
use constant { kdDROP => -999, kdCULL => -1, kdOK => 1, kdMUSTKEEP => 2, kdOther => 0 };  

sub getMapqual {
  my($trname)= @_;
  ## FIXME19: pubset wants trname.align.tab , gmap now makes trname-chrname.align.tab ..
  # (my $mapqual=$trclass)=~ s/\.\w+$/.align.tab/;
  my $mapqual = $settings{'mapqual'} || finddata("$trname*.align.tab") || finddata("geneval/$trname*.align.tab") ; 
  my $ok= ($mapqual and -s $mapqual)?1:0;
  return ($ok,$mapqual);
}

# sub makePubIdtabFromGeneGFF
# add sub makePubIdtabFromTrclass

sub makePubIdtabFromTrclass
{
  my($pubids,$maltab,$trclass)= @_; # ,$genes
  
  # FIXME: 18.03/04 got 300 dup pubids in 300k now, mixup recent..
  
	# ($pubid_format,$altid_format)= make_IDPREFIX(); # default abbrev of $organism now
    #^^ wait for input gff IDs as existing id?
  loggit(0, "makePubIdtabFromTrclass($pubids,$maltab)");

  my($ok, $ndone, $nmiss)= (0) x 9; 
  my($upstatus,$upfiles,$uptemp,$upokids) = (0) x 9;  
  
  # add trclass reading parts of evgmrna2tsa2.pl MAIN_pubsetonly() 
  #	$upstatus=0; 	
  #  ($upstatus, $upfiles, $uptemp, $upokids)  
  #    = update_mrna_fileset($trpath, $cdnaseq, $trimflag, $trimids, @trimfiles); #? leave out here?
  
  ($upokids)= readTrclass($trclass,"publicset");
	# my($maintab,$pubids,$nmaintr,$nalltr)= trclass2maintab($trclass,"publicset",$upokids);
  
  reclassGenes(); # link alts to mains via %main,%alt,%balt id hashes
  
  ($pubid_format,$altid_format,$GBPROID)= make_IDPREFIX() unless($pubid_format); # default abbrev of $organism now
  
  $ok= open(OPD,'>',$pubids); my $outpubidh=*OPD;
  $ok= open(OMA,'>',$maltab); my $outh=*OMA;
  $nmiss=$upokids;
  if($ok) {
    ($ndone,$nmiss)= putPubidTab($outh,$outpubidh); # $mainindex,$ntr
    close($outh); close($outpubidh);
    push @publicset, $pubids, $maltab;
  } else {
    loggit(1, "#ERR: cant putPubidTab($pubids)");
  }
  loggit(0, "done makePubIdtab($pubids,$maltab) nin=$upokids, npub=$ndone, nmiss=$nmiss");
  
	#?add here?  
	## OPTION, dont want as is when have keepdrop cull table
	## BUT do want its added mapqual
	# $ENV{norealt}=1 stops it..
	altreclass_block($trclass,$pubids); # reclassifies alts in pubids table, drops some

  return($ndone,$upokids);
}

=item altreclass or not

  no altreclass: $ENV{norealt} = 1 
  altreclass w/o dropalt:
    -noclasscut= 1
    -maxaltsame= 99999
    -nodrops : added
    >> not yet ENV/opts
    DROPTINYALT = 0
    DROPSHORT_AAPART = 0
    DROPSHORT_AAFULL = 0
    DROPSHORT_ANTISENSE = 0
    
cat publicset/pigevg4wc.pubids.old | cut -f5 | sed 's/a2//; s/hi1/hi/; s/midfrag/frag/; s/cull[12]/cull/; s/cullalt.*
/cullalt/; s/maybe//; ' | sort | uniq -c | head -40

before altreclass
376671 althi
  4264 altmid
  3840 altfrag
 32299 cullalt
 33160 cullmain
 61469 cullnoclass
     3 cullparthi
 30668 main
 11970 noclass
  3515 parthi

after altreclass : is cutting cull tag
357756 althi
  2984 altfrag
 13687 althim       | reclass mains
  4371 altmid
   971 dropaltfrag  | want
 36619 dropalthi    | these 
   373 dropaltmid   | altdrops ?
 64259 main    : but has culls
 73321 noclass : but ann.txt has cull tags
  3518 parthi

=cut

sub altreclass_block {
  my($trclass,$pubids)= @_;

  ## FIX 201405: insert here asmrna_altreclass.pl, .realt_pubids to become new .pubids ..
  my $APPaltreclass= findevigeneapp("rnaseq/asmrna_altreclass.pl",1); # upd 201405 ;  NODIE
  if($APPaltreclass =~ /MISSING/) { $APPaltreclass=""; }
  
	if(not $ENV{norealt} and $APPaltreclass and -x $APPaltreclass) {
	  # ADD OPTIONS: -noclass=$MINAA == drop noclass short things if no other good qualities
	  #  -MAXALTSAME=n == drop excessive number of althi class that appear same by aa size/aa cluster, 
	  #  asmrna_altreclass to be run or not?   other opts?
	  
	  my $realtids="$pubids.realt";
	  my $aopts="-debug -noclasscut=$AAMIN_NOCLASS";
	  $aopts.=" -nodrops" if($NOALTDROPS);
	  
	  # class calls: pubids now override trclass class, from culls, etc. 
	  # need flag to altreclass
	  $ENV{trustpubids}=1; # $aopts.= " -trustpubids";
	  
	  (my $trname=$trclass)=~ s/\.trclass.*//;
	  my($mok,$mapqual)= getMapqual($trname);
	  if($mok) { $aopts.=" -mapqual $mapqual"; }
    # (my $mapqual=$trclass)=~ s/\.\w+$/.align.tab/;
	  # if(! -f $mapqual and -d "geneval") { $mapqual="geneval/$mapqual"; }
	  # if(-s $mapqual) { $aopts.=" -mapqual $mapqual"; }
	  
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



sub makePublicset {
  my($pubids, $outname, $cdnaseq, $genegffset)= @_;
  # from evgmrna2tsa sub MAIN_pubsetonly()
  # $cdnaseq => \@mrnaseq ?
  my $skiptrrun=0;
  #global opt now# my $genenames=""; #?? look for
  our @publicset; #global in pubsets.pm
	my $pubdir="publicset";

  loggit(0,"makePublicset($pubids,$outname)"); # LOG_NOTE,LOG_DEBUG // if($DEBUG);

  ## change this cdnaseq/mrnaseq to array= grep 'mrna|cdna|xxx' @seqset
  my @mrna= grep /\.(mrna|ncrna|rna|cdna|tr)$/, @seqset;
  my @aa  = grep /\.(aa|pep)$/, @seqset; # if empty but @mrna try s/.mrna/.aa/ ?
  my @cds = grep /\.(cds)$/, @seqset;
  #if(not @mrna and $cdnaseq) { @mrna=($cdnaseq); } elsif(@mrna) { $cdnaseq= $mrna[0]; }
  if(@mrna) { $cdnaseq= $mrna[0]; } elsif($cdnaseq) {  @mrna=($cdnaseq); } 
  
  if(@mrna and not @aa) { @aa= map{ (my $f=$_) =~ s/\.\w+/.aa/; $f; } @mrna; }
  if(@mrna and not @cds) { @cds= map{ (my $f=$_) =~ s/\.\w+/.cds/; $f; } @mrna; }
  
  if($genenames and -f $genenames) { }
  elsif( -f "$outname.names") { $genenames= "$outname.names"; } #? need -option
  else { my $pn= makename($pubids,'.names'); $genenames=$pn if(-f $pn); }
  
	my($npubid, $pubidh)= read_pubids($pubids, $cdnaseq); # if($pubids); 
	  # NOT cdnaseq here it isnt publicset/pubid file?? but make_annotab uses cdnaseq annots/IDs
    # return($nred, \%pubids);

	my($annotab, $ntra1, $ncdsa1)=(0) x 9;
	
	## ? add $trclass to make_annot(), or not needed: pubids has that info
	
  use constant OLDANN => 0;
  if(OLDANN) {
    # ** FIXME: make_annotab input  @mrna  AND/OR @genegff annot input
    # UPD: make_annotab() uses pubset.pm global pubids hash info : ADD gff mRNA annots for names, other
    # uses pack global %pubids from read_pubids
	  ($annotab, $ntra1, $ncdsa1) 
		  = make_annotab($cdnaseq, $genenames, $skiptrrun, $outname); # add main/alt pub ids, other geneinfo 
  
  } else {
    $annotab =  makename($outname,".ann.txt"); 
    if(not -f $annotab and -d $pubdir) {
       my($pubd,$ft)= getFileset($pubdir,'.ann.txt');  
       $annotab=$ft if($ft);
    }
    if(ref($genegffset)) { ## not in trclass2pubset
      # DUPID fix: may need keepdrop{ oid } handling here to skip dropped dups by oid..
      ($ntra1, $ncdsa1)= make_annot_from("gff", $annotab, $genegffset,  $genenames);
    } elsif(@mrna) {
      ($ntra1, $ncdsa1)= make_annot_from("mrna", $annotab, \@mrna,  $genenames);
    }
  }
  #NO, make_annot does: push @publicset, $annotab; # if($ntra1);
   
  $FAHDRisValid=1; # export from pubset.pm for annots from mrna hdr
  # my %annothash= read_annotab($annotab); #** fill in annothash from pubidh if no table
  my($nann,$annothash)= read_annotab2($annotab);
  $FAHDRisValid=2 if($npubid and $nann > 0); ## upd1809; FLAG 2 == sensible merge only missing vals

  use evgpubsetsum; # debug; move up
  my(@sumtables)= geneset_summary($pubidh, $annothash);
  if(@sumtables) { 
    my $gs = makename($outname,".genesum.txt");  
    open(TO,'>',$gs); for my $tab (@sumtables) { print TO $tab,"\n"; } close(TO); 
    push @publicset, $gs; 
    }

  # may already have .aa, .cds, .mrna from prior steps. is this ok?
  # add input -seqset params ??
  #     our $RNATYPES='mRNA|ncRNA';

	my($pubmrna,$npm,$minfo) = make_pubseq(\@mrna,'mRNA',$annothash, "$outname.mrna");
	my($pubaa,$npa,$ainfo) 	= make_pubseq(\@aa,'protein',$annothash, "$outname.aa"); # makename($cdnaseq,'.aa')
	my($pubcds,$npc,$cinfo)	= make_pubseq(\@cds,'CDS',$annothash, "$outname.cds"); # makename($cdnaseq,'.cds')
  
  #make_pub does: push @publicset, $pubmrna,$pubaa,$pubcds; # if($npm);
  #UPD: now make_pub puts *.checktab to tmpfiles if no err
  
      # DUPID fix: may need keepdrop{ oid } handling here to skip dropped dups by oid..
  my($pubgff,$ngenes,$ginfo)= make_pubgff($genegffset,$SRC,$annothash, $outname);
  
  my $nout= $npm || $npa || $npc || 1;
  loggit(0,"publicset: ",$pubmrna,$minfo,$pubaa,$ainfo,$pubcds,$cinfo,$pubgff,$ginfo,$annotab); 
  return($nout);
}



sub reclassGenes {
  #?? keepdrop/cull handle here?
  my @amain= grep { not $main{$_} } sort keys %alt; # dropmain here now only for SHOWDROPS !
  foreach my $am (@amain) { 
    my $md= $balt{$am} || $am;  
  	if(!$main{$md} and $drop{$md}) { my $md1= $balt{$md}||""; if($md1 and $main{$md1}) { $md=$md1; } }
  	
    if($main{$md}) { my @at=keys %{$alt{$am}}; map{ $alt{$md}{$_}=$alt{$am}{$_}; $balt{$_}=$md; } @at; } 
    elsif($md) { 
     $main{$am}="NOMAIN";  # FIXME: get rid of these by finding alt main
     } 
  }
  
  foreach my $td (keys %balt) {
    my $md= $balt{$td} || $td; 
    #?need: unless($md) { (undef,undef,$md)= evigene_idparts($td); }
    $main{$md}="NOMAIN" unless($main{$md});# FIXME: get rid of these by finding alt main
  }
}


# evgmrna2tsa2.pl get_evgtrset(); move to package?
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
    
    my $okpath= ($notokay) ? $trpath :"$trpath/okayset"; # should I check if curdir has okayset files?

    #OBSOLETE# my($mrnaOfUtrorf,$nutrorf)= getmRNA_utrorf($okpath,$trname);
        
    ($cdnaseq)= getmRNA($okpath,$trname,$pubdir) if(!$cdnaseq and -d $okpath);
  }
  
  return($cdnaseq,$trpath,$trname,$sradatah);
}


sub readTrclass # from sub trclass2maintab
{
  my($trclass,$pubdir, $okids)=@_;
  my $ntr=0;  my $nerr=0;
  my $mainindex= $pubidnum_start;
  my $hasokids= ($okids and ref($okids) and scalar(%$okids))?1:0;
  
  #... not here?? no output written in this sub
  my $maintab = makename($trclass,".mainalt.tab","trclass");  # > $pt.mainalt.tab
  my $pubidtab= makename($trclass,".pubids","trclass");   
  # if(not -f $pubidtab and $pubdir and -d $pubdir) {
  #   my($pubd,$ft);
  #   ($pubd,$ft)= getFileset($pubdir,'pubids',$pubd);  $pubidtab=$ft if($ft);  
  #   ($pubd,$ft)= getFileset($pubdir,'mainalt.tab',$pubd);  $maintab=$ft if($ft);  
  # }
  # return($maintab,$pubidtab,$mainindex,$ntr) if( -s $maintab and -s $pubidtab);# or dryrun ..

  #globals: (%annotes, %main,%mainsize,%alt,%oclass,%altsize,%altdrops,%didaltdrops,%balt,%drop);  # globals of trclass2maintab
  #unused here: %newmain, %altdefer, %maindrops [unused],  %aaqual, %piad [in annots]    from mrna2tsa
  
  my($ok,$inh,$outh,$outpubidh);
  ($ok,$inh)= openRead($trclass);
  #o.$ok= open($outh,'>',$maintab) if($ok);
  #o.$ok= open($outpubidh,'>',$pubidtab) if($ok);
  unless($ok) { loggit(1,"ERR: parse $trclass TO $maintab"); return; }

  while(<$inh>) {
    next unless(/^\w/); chomp;
    my($td,$ok,$cl,$md,$piad,$aq,$notes)=split "\t";
    unless($cl and $md and $aq) { $nerr++; loggit(1,"ERR: trclass inline:$_"); next; } # what?? report?

		## maybeok needs handling, drop? or keep/mark/cull ??
		## should revise asm dupfliter to avoid these maybes, includes 'refbest', others some good/bad
		## now all are from exoneq flag, all 'althi1|part|frag' classes; 5702 refbest(dups?) of 13059
		## ** maybe keep all maybes ; found some refbest/good being dropped here.. but maybes are large subset
		## should filter by hoscore if avail
		
		## UPD1807 keepdropin should override trclass drop/okay : not working below ***
		
		my $dropit=0;
		my $clorig=$cl;
		if($ok eq 'maybeok') { 
		   if($notes=~/refbest|refgood/) { $ok="okay"; $cl.="maybe"; } 
		   else { $ok="okay"; $cl.="maybe"; } 
		}
		
    if($ok ne 'okay') { $cl=$ok.$cl; $drop{$td}=$cl; $dropit=1;  } #  OPTION: include drops?
    elsif($hasokids and not $okids->{$td}) { $dropit=1; $cl='dropid'.$cl;  $drop{$td}=$cl; } 
      #^^ replace okids with keepdroph 
 
    my($kdxact, $kdxid)=(0,0);
    ## if($KEEPDROPX) { ($kdxact, $kdxid)= check_keepdrop('seqid',"$td,$alloids"); }
    ## if($kdxact =~ /cull|drop/) { $evgclass=$kdxact.$evgclass unless($evgclass =~ /^(cull|drop)/); }
     
    if($keepdropin and (my $kdv= $keepdroph->{$td})) {      
      #? my $oid= $tblinfo->{oid}; # for gff
      # for my $od (split",",$oid){ if(my $kdvo= $keepdroph->{$od}){ $kdv=$kdvo; last; } }
      
      my($kdi,$kact)= split" ",$kdv; #no need kdi?
      if($kdi < 0) {   # drop or cull
        if($kdi <= kdDROP){ 
          $dropit=1; $cl='dropid'.$cl;  $drop{$td}=$cl;
          $kdxact="drop"; $kdxid=$td;
          # $mclass=cDROP; $ok=0; $ndrop++; #  $kact=~/drop/ else drop anything else?
          }
        elsif($kdi == kdCULL){ 
          # $cl= ($isalt or $mainid)?"cullalt":"cull"; # variants?
          # $cl .= "cull"; # suffix? prefix
          $cl = "cull$cl"; # suffix? prefix
          $kdxact="cull"; $kdxid=$td;
          # $mclass= cCULL; 
        } 
      } elsif($kdi>0) { 
        if($dropit) { $dropit=0; delete $drop{$td}; $cl=$clorig; $ok="okay"; }
        # $mclass= cMAIN unless($mclass == cALT || $mclass == cNEWLOC);  
      }  
      
    }
      
    ## all piad entries should have pd and sense/asense fields, placeholder where needed '.' ?
    #o.my($pi,$pa,$asense,$pd)=split"/",$piad; # asense before pd ** NEED To revise asm dupfilter to clarify
    my($pi,$pa,$asense,$pd)= split"/",$piad; # asense before pd ** NEED To revise asm dupfilter to clarify
    map{ $_ ||=""; } ($pi,$pa,$asense,$pd);
    # piad == 100/100/self1 << self1 misused here. means mainid == td ?
    
    if($asense =~ /sense/) { $pd="" unless($pd =~ /^\w/); }
    elsif($asense =~ /^\w/) { $pd=$asense; $asense=""; }
    # if($pd eq "self1") { $md=$td; } # is this right? only for noclass, where td == md already above
    $md=$pd if($pd =~ /^\w/ and not $pd=~/self/);
   	$cl=~s/a2$//;  #? dont need 

		## FIXME here?  new altmap class from eqgene hides main class now; all mains are reclassed altmap, 1 should be left as main.
 		if($dropit and not $SHOWDROPS) {
 			next unless($cl =~ /main|noclass/); # this is enough, skip dropalts
 		}
 		
 		$ntr++;
 		my $aasize= ($aq =~ m/(\d+).(\d*)/) ? "$1.$2" : 0; # size.pCDS for best sort
 		
 		my $mapq= ($nmapqual) ? $mapqualh->{$td} : "";
 		# my $maptab= ($nmapqual) ? $alntabh->{$td} : "";
    $notes.= ",chrmap:$mapq" if($mapq);
 		
 		my($gannot, $tblinfo)= gene_annot_brief($td,"mRNA",$notes);
 		
		$annotes{$td}{aaqual}= $aq;
		$annotes{$td}{piad}= $piad;
	  $annotes{$td}{notes}= $gannot;
		$annotes{$td}{oid}= $tblinfo->{oid};
  	$mainsize{$td}= $altsize{$td}= $aasize; #? dont need both?  aasize{$id}= $aasize;
 		
    #x if($cl =~ /^main|^noclass/) #?? cullmain|cullnoclass 
    if($cl =~ /main|noclass/) 
    { 
    	$main{$td}=$cl; $balt{$td}=$td;
    } else { 
    	$alt{$md}{$td}= $cl; $balt{$td}=$md; 
    	# NOMAIN fix: always add dropmain here 
    }  
  }

	if($SHOWDROPS) {
		# drops from dropset/*.aa headers for perfect_dups, perfect_frags info not in trclass ..
		my $ndr=0;
  	my($dset,$droptr)= getFileset("$trpath/dropset",'drop.tr|drop.aa|drop.cdna');  # dangit need fixed tr/cdna/fa suffix here
  	my($ok,$hin)= ($droptr) ? openRead($droptr) : (0,0);
  	if($ok) { 

  	  ## FIXME.160911: MISSING evgclass=,drop, listings in mainalt.tab .. should include all dropset/*.aa hdr ids
  		# FIXME2: spurious drops getting to pubid table via this $md == main of drop td, but may also be a drop!
  		# below via @amain from alt{md} != this main?
  		## new problem: md here may well not be in mainlist or linked there .. will then leave
  		## these dropped main/alt out of mainalt.tab .. maybe right, these are items of no value?
  		
  		while(<$hin>) { if(/^>(\S+)/) {  
  		  my $td=$1;
  			my($cl,$ok1,$md)= m/evgclass=(\w+),(\w+),match:([^\s;,]+)/;
  			# others:  evgclass=noclass,drop;  and no evgclass=
  			# unless($md) { }
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
	
  return($ntr); 
	  
  # reclassGenes(); # called after this reader     
}



=item gannot pubid annots

  * moved to evigene_pubsets gene_annot_brief
  reproduce this annotation (trclass2mainalt.pl) in notes hash
    aaref:5767,dapsim:Dapsim1EVm000004t1,chrmap:100a,98i,25555l,33x,29sc:324762-356329:+,pflag:0
    
  my @ANKEY_MRNA= qw(aalen  cov pid nexon clen namealn Dbxref scoresum);
  leave out special annots of altbest, alttr, altid/mainid/newlocusid

  gannot moved to evigene_pubsets gene_annot_brief($id,@mrnarow)

=cut

=item readGeneIDs

  read gene class info + annots from evigene processed gff gene sets,
  i.e., main + alternate genes processed for map overlaps, best locus representatives
  
=cut

use constant KEEP_C => 0; # split ID _C[12..n] tag ; this is problematic, have keep_C1, cull_C2

# sub readGeneIDs {
#   my($inh,$infile)= @_; #? or while(<>)
# 
#     ## add cCULL, cDROP? cKEEP => main or alt
#   use constant { cNone=> 0, cMAIN => 1, cALT => 2, cNEWLOC => 3, cCULL => 4, cDROP => 5, cOther => 6 };
#   
#   sub mainofaltid { my($aid,$ismain)=@_; 
#     unless($ismain) { 
#       my($t)= $aid=~s/(_[CG]\d+)//?$1:""; 
#       $aid=~s/t\d+$//; 
#       $aid.="t1" unless($aid=~/t\d+$/);  
#       if(KEEP_C or $t=~/_G/){ $aid.=$t; } } 
#     return($aid); 
#   }
#   
#   my($src,$ok,$inid,$nok,$ndrop,$nskip,$firstID)=(0) x 9;
#   my(%did); # FIXME need global over all in.gff >> use %balt{id} global
#   my $OTHERCLASS= cOther; # == cCULL for CULLOTHERS
#   
#   while(<$inh>) {
#     if(/^\W/){
#       # ? collect scoresum of '#x.mRNA' of bestmain.gff ?
#       # %mainids has to be global, and read main.gff before alts.gff
#       if(/^#x./ and /\t($RNATYPES)/) { 
#         my($id)= m/\bID=([^;\s]+)/; 
#         ## KEEP_C here ?? 
#         my($ovd)=m/\boverids=([^;\s]+)/g; 
#         if($ovd){ $mainids{$id}= (split",",$ovd)[0]; }
#         # my($gscore)= m/\bscoresum=([^;\s]+)/g; # gscore or scoresum or col5
#         # my($skip)= m/\bskip=([^;\s]+)/g; # gscore or scoresum or col5
#         #my @an= m/((?:scoresum|overids|skip)=[^;\s]+)/g; 
#         } 
#       next;
#     }
# 
#     if(/\t($RNATYPES)/) { # mRNA
#     $ok=1;
#     my @gff=split"\t";
#     my $gsrc= $gff[1];
#     # gannot(): my $gtype= $gff[2]; # mRNA | ncRNA .. preserve this in pubid tab?
#     my $gann= $gff[8];
#     
#     # FIXME: add cull/skip class, ($class)= $gsrc =~ m/(alt|cull)$/ ;
#     # culls may have been main or alt, ie have newaltid
#     # use separate input keep_drop_cull.table ?
#     
#     $src= $SRC || $gsrc; #? use input gff[1] col unless defined? not used here now?
#     my($id)= m/\bID=([^;\s]+)/;  $inid= $id;
#     #o.my($gpre,$gnum,$gd,$ti,$isplit)= evigene_idparts($id);
#     my($gpre,$gnum,$ti,$isplit,$gdup)= evgIdParts($id);
#     
#     # should add mainid=XXX to problem cases, ie byhand.gff alts
#     my($mainid)=m/\bmainid=(\w+)/?$1:0; # trust this annot? require isalt also?
#     unless($mainid) { $mainid= $mainids{$id} || 0; }
#     my($oad)=m/\boldaltid=(\w+)/?$1:0; my($nad)=m/\bnewaltid=(\w+)/?$1:0;
#     my($old)=m/\boldlocus=(\w+)/?$1:0; my($nld)=m/\bnewlocus=(\w+)/?$1:0;
#     my($scoresum)= (m/scoresum=([\d-]+)/)?$1:0;
#     my($aq)= (m/aalen=([^;\s]+)/)?$1:0;
#     my $aasize= ($aq =~ m/(\d+).(\d*)/) ? "$1.$2" : 0; # size.pCDS for best sort
#   
#     unless(KEEP_C){ $id=~s/_C\d// if($isplit); $mainid=~s/_C\d$//; } # want this?
# 
#     # dang _C1,2 parts separate? ** FIXME, need to see both _C1,_C2
#     if($did{$inid}) { $ok=0; next; } 
#     elsif($balt{$id}) {   #  or $did{$id}
#       # BUG below? getting dup main entries when shouldn't, from extra in.gff?
#       if($inid eq $id) { $ok= 0; }
#       else { if(my $ocl= $oclass{$id}) { $ok=($ocl=~/^cull/)?1:0; } }
#       next if($ok == 0); 
#     } 
#     
#     #* may need annot in pubid table for _C/_G edits, exists : oid retains?    
# ## ugh altbest id t1_C2 bug : need to cut out 'tNNtNN' 2nd alt tag?
# # Danrer6pEVm022697t4_C2	noclass	
# # Danrer6pEVm022697t4t2_C2t1	NOMAIN	Danrer6pEVm022697t4t2_C2/alt
# # Danrer6pEVm022697t4t3_C2t1	NOMAIN	Danrer6pEVm022697t4t3_C2/alt
# # grep ID=Danrer6pEVm022697t chr21gn.malt.altbest.gff | cut -f9 | sed 's/aamap=.*//; '   
# # ID=Danrer6pEVm022697t4t2_C2;oldaltid=Danrer6pEVm022697t2_C2;trg=Danrer6pEVm022697t2 67 1018;
# # ID=Danrer6pEVm022697t4t3_C2;oldaltid=Danrer6pEVm022697t5_C2;trg=Danrer6pEVm022697t5 78 650;
# # ID=Danrer6pEVm022697t4t4_C2;oldaltid=Danrer6pEVm022697t7;trg=Danrer6pEVm022697t7 91 495;
# # ID=Danrer6pEVm022697t4t5_C2;oldaltid=Danrer6pEVm022697t9_C2;trg=Danrer6pEVm022697t9 98 426;
# 
# =item split bug 180330
# 
#  grep '^vDanrer6pEVm004866t1     ' chreset-map6na.pubids | head
#  vDanrer6pEVm004866t1	Danrer6pEVm005074t4	vDanrer6pEVm004866	1	main	788,86p,complete	0	aaref:726,zfish16nc:NP_571403.1,chrmap:87a,99i,2733l,11x,chr20:54192943-54206663:-	Danrer1a_sRn1l1SRR4994225velvk79Loc1037t4,Danrer6pEVm005074t4
#  vDanrer6pEVm004866t1	Danrer6pEVm005074t152	vDanrer6pEVm004866	1	main	725,74p,complete	0	0,0,chrmap:99a,99i,2930l,12x,chr20:54192871-54222720:-	tridba1a_sNn12l1SRR1524240ridk81Loc231231,Danrer6pEVm005074t152
# 
#  grep '^Danrer6pEVm005074t4_' $pt.map6k.keepdrop.gscore     
#  Danrer6pEVm005074t4_C2	keep	gscore=4267	vin=10	vho=726	pho=100	altpar=1	aalen=788	cov=87	acov=95	lspan=1.14	locus=chr20:54192943-54206663:-	gclass=alt4	flags=split,evuncull:h0/i1	gsold=4095
#  Danrer6pEVm005074t4_C1	cull	gscore=368	vin=0	vho=726	pho=100	altpar=1	aalen=788	cov=10	acov=4	lspan=0	locus=chr20:54190464-54192154:-	gclass=main1	flags=cull3.ovEVGm001440t9/99.96,split	gsold=368
# 
#  chreset-map6na_pub.gff
#  chr20	zf17evgm6v	mRNA	54192871	54222720	4409	-	.	ID=vDanrer6pEVm004866t1;trg=Danrer6pEVm005074t152 1 2925;cov=99.8;nexon=12;pid=99.4;clen=2930;offs=544-2721;oid=Danrer6pEVm005074t152,tridba1a_sNn12l1SRR1524240ridk81Loc231231;aalen=725,74p,complete;cdsoff=544-2721;Name=Uncharacterized protein;gscore=4335;vin=9;vho=725;altpar=1;acov=100;cxlen=2178/2925,74%;inexon=12/12/10;scoresum=4409
#  chr20	zf17evgm6v	mRNA	54190464	54192154	413	-	.	ID=vDanrer6pEVm004866t1_C1;Split=1;trg=Danrer6pEVm005074t4 1 278;aamap=38;cov=10.2;nexon=2;pid=98.9;path=1/6;chim2=chr20:54192943-54206663:-;Name=Heat shock protein HSP 90-alpha 1;clen=2733;offs=167-2533;oid=Danrer6pEVm005074t4_C1,Danrer6pEVm005074t152,Danrer1a_sRn1l1SRR4994225velvk79Loc1037t4;aalen=788,86p,complete;cdsoff=167-2533;namealn=100p,726/726,788;Dbxref=zfish16nc:NP_571403.1,;gscore=368;vin=0;vho=726;altpar=1;acov=4;flags=split,;cxlen=112/277,40%;inexon=2/2/1;scoresum=413
#       ^^ should be cull
#  chr20	zf17evgm6v	mRNA	54192943	54206663	4095	-	.	ID=vDanrer6pEVm004866t1_C2;Split=2;trg=Danrer6pEVm005074t4 279 2665;aamap=750;cov=87.3;nexon=11;pid=99.2;path=2/6;chim1=chr20:54190464-54192154:-;Name=Heat shock protein HSP 90-alpha 1;clen=2733;offs=167-2533;oid=Danrer6pEVm005074t4_C2,Danrer6pEVm005074t152,Danrer1a_sRn1l1SRR4994225velvk79Loc1037t4;aalen=788,86p,complete;cdsoff=167-2533;namealn=100p,726/726,788;Dbxref=zfish16nc:NP_571403.1,;gscore=4095;vin=10;vho=726;altpar=1;acov=95;flags=split,
#       ^^ should be vDanrer6pEVm004866t2 or other alt
# =cut
# 
#     my $cl="none";
#     my $isalt=(/;alttr=(\d+)/)?$1:0; # ;isalt=1 ?
#     my $isother=0; # NOT isalt or ismain, ie iscull + other
#     my $mclass = cNone; 
#     my $md= $mainid || $id; ## $gd; # evigene_idOfparts($gpre,$gnum,$gd,1,0); #? use $gd? drop ti?
#     
#     my($gannot,$tblinfo)= gene_annot_brief($id, @gff); # was gannot()
# 
#     ## opt to reparse of last round of pub,cull.gff : no oad/scoresum/..
#     if($CLASSbySOURCE) {
#       my $gcla= ($gsrc=~m/(alt|cull|nc)$/)?$1:"";
#       my $iscull= ($gcla eq "cull")?1:0;
#       $md=  mainofaltid($id); # $mainid || : outof date mainid= tag skip or update
#       $isalt= ($md ne $id or $mainid)? 1 : 0; # for cullalt
#       $mclass= ($iscull)? cCULL # need cullalt also
#           : ($gcla eq "alt")? cALT 
#           : ($gcla eq "nc") ? cMAIN 
#           : ($md eq $id) ? cMAIN 
#           : cOther;
#     }      
#     elsif($oad or $nad) { $mclass = cALT; $md= $mainid || mainofaltid($nad || $id); } # EVm000t2t3 hack format, or EVm00t3
#     elsif($old or $nld) { $mclass= cNEWLOC; $md= $mainid || mainofaltid($nld || $id, not $isalt); }
#     elsif($scoresum) { $mclass= cMAIN; $md= $mainid || $id; } # or $ti < 2 ?? cant assume?
#     else {
#       # my ($INS)= $infile =~ m/(\w+).gff/; #? file name not good for this
#       # $isalt=1 if($INS =~ /altbest/); $isnewloc=1 if($INS =~ /newloci/);
#       ## mainofaltid() guess doesnt work well enough; add mainid=XXX annot; use overids=?
#       if($mainid) {
#         if(my $mcl= $main{$mainid}) {
#           $isalt=2; # guess this .. for keep.other group w/o altbest info
#         } else { # check if mainid now is alt?
#           if(my $mmd= $balt{$mainid}) {
#             if(my $acl= $alt{$mmd}{$mainid}) {
#               $isalt=2; $mainid= $mmd; # maybe.
#             }
#           }
#         }
#       }
#       if($isalt) { $mclass = cALT; $md= $mainid || mainofaltid($id); } 
#       else { 
#         $mclass= cOther; #or $OTHERCLASS; #was $mclass= cMAIN;  # FIXME: option for -otherclass=cull|xxx|remainder..
#         $md= $mainid || $id; 
#         }
#     }
#     # if($mclass == cALT and $id eq $mainid) { $mclass= cMAIN; } #?? check here, abive?
#     
#     # if($iscull) { $mclass= cCULL; } #? like drop? but output all data w/ new ids, separate?
#     ## KEEP_C=0 bug : keepdrop id can have _C tag
#     ## FIXME: are loosing some _C keep due to _C1 keep, _C2 cull ;  ie KEEP_C=0 bug
#     if($keepdropin and (my $kdv= $keepdroph->{$inid} || $keepdroph->{$id})) {
#       ## add dup pubid check ; add dupid flag in kdv?
#       ## may need to keep row in pubid tab for dup id cases, even drops, so dup data are handled right
#       ## but then need uniqof(dupid) ie. {pubid}d2,d3..
#       my $oid= $tblinfo->{oid}; 
#       for my $od (split",",$oid){ if(my $kdvo= $keepdroph->{$od}){ $kdv=$kdvo; last; } }
#       
#       my($kdi,$kact)= split" ",$kdv; #no need kdi?
#       if($kdi < 0) {   # drop or cull
#         if($kdi <= kdDROP){ $mclass=cDROP; $ok=0; $ndrop++; } #  $kact=~/drop/ else drop anything else?
#         elsif($kdi == kdCULL){ $mclass= cCULL; } #  $kact=~/cull/
#       } elsif($kdi>0) { 
#         $mclass= cMAIN unless($mclass == cALT || $mclass == cNEWLOC); ## cKEEP ? treat like cMAIN if no other
#       } # keep/ok
#     }
#     # if($dropids) { # fixme for map5best too liberal gscore, remove no(ho+in)evidence loci
#     #   if(my $ds=$dropid{$id}) { if($gsrc =~ m/$ds/){ $ok=0; $ndrop++; } }
#     # }
#     
#     if($mclass == cALT) { 
#       if($oad and $did{$oad}) { $ok= 0; }
#       elsif($nad and $did{$id}) { $ok= 0; }
#       elsif($oad) { $id= $oad; } #?? need this to stop tNNtNN hack bugs?
#       
#       ## this id change is problem: tNtN; when main is t2, alt is t1, altbest assumes main new id t1..
#       ## ID=Danrer6pEVm004586t2t7; oldaltid=Danrer6pEVm004586t1;
#       $isalt=1;
#       $cl="alt"; # variants?
#       $src .= "alt";
#     	# $alt{$md}{$td}= $cl; $balt{$td}=$md; $altsize{$td}= $aasize; 
#     
#     } elsif($mclass == cMAIN) { 
#       $isalt=0;
#       $cl="main"; # variants?
#       # $main{$md}=$cl; $balt{$td}=$md;   $mainsize{$td}= $aasize;  
#  
#     } elsif($mclass == cNEWLOC) { 
#       $ok= ($id=~/^newloc/ or $nld =~ /^newloc/)?1:0; #? others? what are ids
#       if($old and $did{$old}) { $ok= 0; }
#       elsif($nld and $did{$id}) { $ok= 0; }
#       #X elsif($old) { $id= $old; } #??
#       # $isalt=(/;alttr=1/)?1:0; # see above
#       $src .= "alt" if($isalt);
#       $cl=($isalt)?"altnew":"mainnew"; # variants?
# 
#     } elsif($mclass == cDROP) {  
#       $ok=0; $nskip++;
#     
#     } elsif($mclass == cCULL or $CULLOTHERS) {  # CULLOTHERS = opt to treat cNone/Cother as cCULL ?
#       #? $isalt=0; #? maybe yes/no
#       # $isalt=2; #? $md= $id; # isolate is noalt/nomain locus?
#       #** cull fails unless isalt ? need some check of main{id} == cull class
#       $cl= ($isalt or $mainid)?"cullalt":"cull"; # variants?
#       # $cl="cull";  
#       $src .= "cull";
#       $mclass= cCULL;
#       $isother=1; # separate from noncull alt/main
#       $ok=1;
#       
#     } else {  
#       # change old ID t2..N to t1 ?
#       $ok=0; $nskip++;
#     }
#     
#     # 0330: still bugs with dang splits, ie. cull _C1, keep _C2 as alt
#     #  .. then pubids dup new ids, n=47, one is _C2 alt, other is main of _C2 alt
#     if($ok and $did{$inid}) {  
#       $ok= 0; # unless $dids eq $INS for splits ? no splits same ID ?
#     }
#      
#     if($ok) {
#       $nok++;
#       # ugh ID changes from oldaltid are problem..
#       # $aaqual{$id}= $aq; #  for pubtable
#       # $notes{$id}= $gannot; #  for pubtable, add annots
#       # $piad{id}= $piad; #  for pubtable, align to main mrna: ident%,aln% == 99,80
#       $annotes{$id}{aaqual}= $aq;
#       $annotes{$id}{piad}= 0;
#       $annotes{$id}{notes}= $gannot;
#       $annotes{$id}{oid}= $tblinfo->{oid};
#       
#       $mainsize{$id}= $altsize{$id}= $aasize; #? dont need both?  aasize{$id}= $aasize;
#       $oclass{$id}= $cl; # use for all classes?
#       if($isother) {
#     	## oclass adds problems, treat like %main? add to 
#     	$oclass{$id}= $cl; $balt{$id}=$md; $main{$md}="main";
#     	#?? maybe need:  $alt{$md}{$id}= $cl if($id ne $md);
#       } elsif($isalt) {
#     	$alt{$md}{$id}= $cl; $balt{$id}=$md; 
#       } else {
#       # FIXME cullalts here cause problems (=> cull main)
#       $firstID= $id unless($firstID); # but not _C/G suffix
#       $main{$md}=$cl; $balt{$id}=$md; # dang, should be md == id 
#       }
#        
#       $did{$id} = $did{$inid} = $src;
#     }
#     
#   }
#   
#   # elsif(/\texon\t/){
#   #
#   # } 
#   
#   }  
#   
#   return($nok,$ndrop,$nskip,$firstID); # $nmain,$nalt,...
# }


sub putPubidTab {
  my($outh, $outpubidh)= @_; # globals?
 
  #globals: (%annotes, %main,%mainsize,%alt,%oclass,%altsize,%altdrops,%didaltdrops,%balt,%drop);  # globals of trclass2maintab
 
  ## headers 
  #mainalt: originalID     MainClass  Alternates
  #pubid  : Public_mRNA_ID originalID      PublicGeneID    AltNum  Class   AAqual  pIdAln  Notes

  print $outh '#'.join("\t",qw(originalID MainClass Alternates))."\n";
  print $outpubidh '#'.join("\t",qw(Public_mRNA_ID originalID PublicGeneID AltNum Class AAqual pIdAln Notes Oids))."\n"
    if($outpubidh);
  #upd1802: pubid add Oids col
  
  my %doneid=();
  my @mainlist;
  # %oclass adds problems == main or alt
  
  #* BUG bad %main id set?? %balt has all; this doesnt fix missing doneid
  # problem in %alt{md} **
  foreach my $ad (keys %balt) {
    unless( $drop{$ad} ) {
    my $md= $balt{$ad} || $ad;  
    unless( $main{$md} ) { $main{$md}="NOMAIN"; } # NOMAIN ?
    unless( $md eq $ad or $alt{$md}{$ad} ) { 
      my $omc= $oclass{$ad} || "altmiss"; # $cull= $omc if($omc); #?
      $alt{$md}{$ad}= $omc; } # culls here ? also mains only for ad ne md?
      }
  }

  if($SIZESORT) { # or NOT IDSORT ?
  @mainlist= sort{ ($mainsize{$b}||0) <=> ($mainsize{$a}||0) or $a cmp $b } keys %main;
  } else { # IDsort, only special cases, like merge old/update, keeping old order for most
  @mainlist= sort{  $a cmp $b } keys %main;
  }

  #upd1708: preserveOldIds option ..  keep input gene ids as best feasible, including altnum?
  #  need all alts per main, check evg gene ids for 1 x many, also need check across loci for 2+ new loci w/ old geneids
  #  where cant preserve geneid for all of locus, do what? new geneid, or modify old?
  #  UPD: want/need main = t1, cant preserve old alt nums if alt classing changes ..

  my ($nnewids,$newidh)= ($preserveOldIds) ? preserveOldIds(\@mainlist, \%alt, \%drop, \%altsize) : (0, undef);
      #,$newpubidh ret
      
  foreach my $md (@mainlist) { 

    my @ad= sort{ $alt{$md}{$a} cmp $alt{$md}{$b}
      or $altsize{$b} <=> $altsize{$a} or $a cmp $b } keys %{$alt{$md}}; 

    my $culls={}; #? ($CULLXEQ) ? cullExonEq($md,\@ad,\%alt,\%notes) : {}; #?? here
    my $ad= join",",map{ "$_/".$alt{$md}{$_} } @ad; 
    my $mc= $main{$md}; 
    ## oclass cullalt mix up main class
    # if($mc=~/other/){ $mc=$oclass{$md}; }
    
    # FIXME: change $mc eq "noclass" to "main" if alts exist, and vversa
    # FIXME2: preserve $mc eq "cull" .. use %culls{id} ?
    my($cull)= ($mc =~ m/^cull/) ? "cull":""; # this is problem? replace w/ omc=oclass{md}
    my $omc= $oclass{$md} || ""; # $cull= $omc if($omc); #?
    
    # FIXME cullalt in @ad dont count to main/noclass
    my $altcount= @ad;
    $altcount= scalar( grep{ my $c= $oclass{$_}||""; $c !~ m/cull/ } @ad );
    
    if($mc =~ /^NOMAIN/) { 
      # do below
    } elsif($altcount>0 and $mc !~ /^main/) {
      $mc =~ s/^\w+/main/;
      $main{$md}= $mc= $cull.$mc; 
    } elsif( $altcount==0 and $mc !~ /^noclass/) {
      $mc =~ s/^\w+/noclass/;
      $main{$md}= $mc= $cull.$mc; 
    }
    
    if($SHOWDROPS) {
    	my @add= sort{$altdrops{$md}{$a} cmp $altdrops{$md}{$b}} keys %{$altdrops{$md}};  
    	map{ $didaltdrops{$_}=1 } @add; # later dump not didaltdrops
    	my $add=join",",map{ "$_/".$altdrops{$md}{$_} } @add; 
    	$ad .= ",$add" if ($add);
    }
    
    my $mcout=$mc;
    if($omc) { $mcout=$omc.$mc if($omc=~/^(cull|drop)/); } #? both?
    print $outh join("\t",$md,$mcout,$ad)."\n";  # mainalt.tab
    
    if($outpubidh) { # should be required ??
      my $ialt= 0; my $needmain=0;
      my($cla,$aaq,$pida,$nots,$oids);
      $cla= $mcout;  # cla=$main{$td}=$cl;  ; changed above?
      # ** $cla=  $main{$md}||"main"; # need main == cull now
      # $aaq= $aaqual{$md}||"noaa";
      # $pida=$piad{$md}||0;  
      # $nots=$notes{$md}||"nonote";  
      $aaq= $annotes{$md}{aaqual}||"noaa";
      $pida=$annotes{$md}{piad}||0;  
      $nots=$annotes{$md}{notes}||"nonote";  
      $oids=$annotes{$md}{oid}||"noid";  

      if($mc eq "NOMAIN") { $cla=$cull . (($altcount>0)?"main":"noclass"); } ## needs to change, to main? to noclass?
      elsif($mc =~ /^alt/) { } # is this were nomain show up? or @ad?

      my @sad = @ad;      
      if(1) {  # $SIZESORT .. want to use altbest sort order, from tNN altnum
        @sad= sort{ $altsize{$b} <=> $altsize{$a} or $a cmp $b } @ad;      
      }
      
      if($drop{$md} or $doneid{$md}) { $needmain=1; }
      else {
      	$mainindex++; $needmain=0; # BUG: move below drop{}
      	my($pubmrnaid,$pubgeneid,$pubti)= make_pubid($md, $mainindex, ++$ialt, $newidh);
      	print $outpubidh join("\t",$pubmrnaid,$md,$pubgeneid,$pubti,$cla,$aaq,$pida,$nots,$oids)."\n";  #n
      	$ntr++; $doneid{$md}++;
      	}
      
      foreach my $ad (@sad) {
        unless($drop{$ad} or $doneid{$ad}) {
        $cla=  $alt{$md}{$ad}||"nocl"; 
        $aaq= $annotes{$ad}{aaqual}||"noaa";
        $pida=$annotes{$ad}{piad}||0;  
        $nots=$annotes{$ad}{notes}||"nonote";  
        $oids=$annotes{$ad}{oid}||"noid";  
        # $cull= $culls->{$ad}||""; # $cla == "cullalt..";
      	if($needmain) { 
      	  $mainindex++; $needmain=0; 
      	  if($cla=~/^alt/){ $cla=(@sad>1)?"main":"noclass"; } 
      	}  
        my($altmrnaid,$altgeneid,$alti)= make_pubid($ad, $mainindex, ++$ialt, $newidh);
        print $outpubidh join("\t",$altmrnaid,$ad,$altgeneid,$alti,$cla,$aaq,$pida,$nots,$oids)."\n"; 
        $ntr++; $doneid{$ad}++;
        }
      }
    }
  }
  
  # double check not-done .. missing what? LOTS .. missing from %main ??
  my($notdone)=(0);
  foreach my $ad (keys %balt) {
    unless( $drop{$ad} or $doneid{$ad} ) {
    my $md= $balt{$ad} || $ad; 
    $altdrops{$md}{$ad}++;  $notdone++;
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
        my $mc= $main{$md}||"NOMAINd";  # fixme: oclass{$md}
        print $outh join("\t",$md,$mc,$misdrop)."\n";  # mainalt.tab
      }
    }
  }
  return($ntr,$notdone);  
  # close($inh); close($outh);
  #? return($maintab,$pubidtab,$mainindex,$ntr);  # return main,alt id hashes ....
}


sub preserveOldIds {  # for trclass/okayseq set with input -preserve oldpubid.tab 
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

sub preserveOldIds_FOR_GFF {
  my($mainlist, $altsOfMain, $drop, $altsize)= @_;
  my(%gids, %gnums, %gdone, %newids, %newpubids, $gprefix);
  my $nids=0;
  #above# my $NEWIDpre='n';
  
  ## prefix, gnum problem: have mixed gprefices .. diff gnum set for each
  ## FIXME  gnum only for $gprefat =~ /$IDPREOK/

  my $GNEXTNUM=0;
  foreach my $md (@$mainlist) {
    my @okd = grep{ not($drop->{$_}) } ($md, keys %{$altsOfMain->{$md}} );
    for my $id (@okd) {
      #o.my($gpre,$gnum,$gd,$ti,$isplit)= evigene_idparts($id);
      my($gpre,$gnum,$ti,$isplit,$gdup)= evgIdParts($id);
      next unless($gnum and $gpre =~  m/$IDPREOK/);
      $gnums{$gnum}++; # $gids{$gd}{$ti}= $id;
      $gprefix= $gpre unless($gprefix);
      }
  }
  my($glast)= sort{ $b <=> $a } keys %gnums;
  $GNEXTNUM= 9 + $glast;
  
  my $idformat= $NEWIDpre . $gprefix . '%06d'; #?  use pubid_format of make_IDPREFIX  ?
  my(%havepubg);

  foreach my $md (@$mainlist) {
  
    # my @ad= sort{ $altsOfMain->{$md}{$a} cmp $altsOfMain->{$md}{$b} # class sort
    #  or $altsize->{$b} <=> $altsize->{$a} or $a cmp $b } keys %{$altsOfMain->{$md}}; 
    my @ad= sort{ $a cmp $b } keys %{$altsOfMain->{$md}}; 
    my @okd = grep{ not($drop->{$_}) } ($md,@ad);
  
    %gnums=(); %gids=();
    my $gnumfirst=0; my $gprefat=0;
    for my $id (@okd) {
      #o.my($gpre,$gnum,$gd,$ti,$isplit)= evigene_idparts($id);
      my($gpre,$gnum,$ti,$isplit,$gdup)= evgIdParts($id);
      next unless($gnum); # for _G2,n.. new loci
      $gnums{$gnum}++; $gids{$gnum}{$ti}= $id;
      $gnumfirst=$gnum unless($gnumfirst or $gdone{$gnum}); 
      $gprefat= $gpre unless($gprefat); #?? problems? yes
    }
    
    my @gnums= sort{ $gnums{$b}<=>$gnums{$a}  or $a <=> $b } keys %gnums;
    unless($gnumfirst) { # pick most or first in (main) <<
      for(my $i=0; $i<=$#gnums; $i++) {
        unless($gdone{$gnums[$i]}) { $gnumfirst=$gnums[$i]; last; }
      }
    }
    # usage assumption: -keepids=IdPreA -idpre IdPreA, so inval gprefat will go to ++GNEXTNUM 
    unless($gprefat and $gprefat =~ /$IDPREOK/){ $gprefat=$IDPREFIX;  $gnumfirst=0; }
    $idformat= $NEWIDpre . $gprefat . '%06d'; # maybe new for each main id??
    my($pubgn,$haveit)=(0,0);
    do { 
      unless($gnumfirst) { $gnumfirst= ++$GNEXTNUM; } #??
      $pubgn = sprintf( $idformat, $gnumfirst);
      $haveit= $havepubg{$pubgn} || 0;
      ## fixme dup gene ids now, check with t1
      unless($haveit) { $haveit=1 if($newpubids{$pubgn . 't1'}); }
      $gnumfirst=0 if($haveit);
    } while ($haveit);
    $havepubg{$pubgn}=1; 

    # try2
    my %tidone=(); my $timax= @okd; #? not same using orig ti > @okd
    for my $id (@okd) {
      #o.my($gpre,$gnum,$gd,$ti,$isplit)= evigene_idparts($id);
      my($gpre,$gnum,$ti,$isplit,$gdup)= evgIdParts($id);
      my $tinew=0;
      if($gnum eq $gnumfirst and not $tidone{$ti}) { $tinew=$ti; }
      if($tinew==0) { do { $tinew++; } while( $tidone{$tinew} ); } # can put lowqual alts at ti top
      $tidone{$tinew}++; 
      $gdone{$gnumfirst}++;
      #above# my $pubgn = sprintf( $idformat, $gnumfirst); 
      my $pubti = sprintf( $altid_format, $tinew); # same ti or new? cant use same ti w/o checks
      my $pubid=  $pubgn . $pubti;
      $newids{$id}= $pubid;
      $newpubids{$pubid}= $id;
      $nids++;
    }
  } # mainlist
  
  return($nids,\%newids,\%newpubids);
} # sub preserveOldIds




1;

__END__
