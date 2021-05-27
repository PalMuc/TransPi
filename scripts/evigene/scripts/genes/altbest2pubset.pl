#!/usr/bin/env perl
# evgs/genes/altbest2pubset.pl 

=item altbest2pubset.pl

 evigene/scripts/genes/altbest2pubset.pl 
  from evigene/scripts/prot/trclass2mainalt.pl + map2bestmerge.pl + other?
 see also evigene/scripts/genes/evigenegff.pm parts?

=item need output pub.gff, pub.seq outputs w/ annots
  here? or better separate perl: read pubid table, gff or seq inputs
  output updated publicset/ new name
  see evgmrna2tsa pipe for parts, use those

=item add input table keep/drop/cull reclass table

=item main - cullalt = noclass (no alts)

=item FIXME mainid= needed in altbest.pl gff out

 id tags from altbest are a mess now; use alt.tab?
 perl -p -e 'if(/\tmRNA/){ if(($ad)=m/;newaltid=(\w+)/) { 
  $md=$ad; $md=~s/(EVm\d+)t(\d+)$/${1}t1/;  $md=~s/(EVm\d+)t(\d+)_C/${1}t1_C/; $md=~s/(\d[dt]\d+)t\d+/$1/; 
    $md=~s/(_[CG]\d+)t\d+$/$1/;  s/;newaltid=$ad/;mainid=$md/; }  } ' \
  chr21t.map5best.gff > chr21t.map5bfix.gff
   
=cut

use FindBin;
use lib ("$FindBin::Bin","$FindBin::Bin/.."); # assume evigene/scripts/ layout: main, prot/, rnaseq/, ...
use strict;
use warnings;
use Getopt::Long;
use cdna_evigenesub; 
use evigene_pubsets; # now has some of below subs

use constant VERSION => '2018.03.31'; # '2018.02.18'; # from trclass2mainalt + main,altbest gff merge

our $EVIGENES= "$FindBin::Bin/.."; # ="$FindBin::Bin";  # WRONG place ; this is in evgs/genes/ subdir ..
our $DEBUG= $ENV{debug}|| 0;
our $IDPREFIX= $ENV{idprefix} || 'EVGm'; 
my  $NEWIDpre= $ENV{idtag}||'n';
my  $SRC= $ENV{source}||"evm"; # fixme: for gff output, source.class col2 src_(main/alt/other);
our $SHOWDROPS=1; #default on? only for mainalt.tab outputs
our $SIZESORT=1; 
our $FAHDRisValid=1; # export from pubset.pm for annots from mrna hdr
our ($pubid_format,$altid_format,$GBPROID);
our $pubidnum_start=0;
our $preserveOldIds=undef; # change to ok pattern, IDPREOK
my $CULLOTHERS= 0; # how to treat unclassed input genes
our $SORT_PUBSET= 0;
#? our $CULLXEQ=0;
our $DATE;  $DATE=`date '+%Y%m%d'`; chomp($DATE);
our $RNATYPES; # our $RNATYPES='mRNA|ncRNA'; # in evigene_pubsets.pm; transcribed_RNA mrnatypes exontypes 
our ($KEEPDROPX,%keepdropxtab); # evigene_pubsets.pm
my $ONLY_PUBIDS= 0; my $CLASSbySOURCE= 0;

my ($trclass,$output,$logfile,$keepdropin);  
my ($trpath,$trname, $sradatah, %settings);
my ($oname, @genes, @seqset, $genenames);

my $optok= GetOptions(
  "genes|input=s", \@genes, # many gff inputs : use perl input pipe while(<>) ?
  "mrna|cdna|sequences=s", \@seqset, # many inputs? 
  "oname|outbase|output=s", \$oname,  # FIXME: idpre option  overwritten by spppref
  "idprefix=s", \$IDPREFIX,  
  "tagidnew=s", \$NEWIDpre,   
  "source=s", \$SRC,   
  "names=s", \$genenames,   
  "keepoldids|preserveOldIds:s", \$preserveOldIds,  # allow -idpre=XXX -keepold w/o 2nd param
  "keepdropin=s", \$keepdropin,   # other opt key? reclassin?
  "dropshow!", \$SHOWDROPS,  
  "sizesort!", \$SIZESORT,  # default:on ; may be required to get proper ids w/ preserveOld
  "pubsort!", \$SORT_PUBSET,  # default:off?
  "CULLOTHERS!", \$CULLOTHERS,  # or otherclass=cull|skip|remainder name?
  "onlypubids|ONLY_PUBIDS!", \$ONLY_PUBIDS,  #CLASSbySOURCE also : same opt?
  "classbysource!", \$CLASSbySOURCE,  # default:off?
  # "CULLXEQ!", \$CULLXEQ,  
  # "dryrun|n!", \$dryrun, 
  "debug!", \$DEBUG, 
  );

push @genes, grep(/gff/, @ARGV); # remainder gff ? grep gff ?

die "EvidentialGene altbest2pubset -genes main.gff alts.gff [-idprefix $IDPREFIX ] 
  makes tables of public ids and main-alt linkage, from main,alt.gff of alt2best.
  version ", VERSION, "\n"
  unless($optok and (@genes or @seqset));  

# my $IDPREOK= ($preserveOldIds and $preserveOldIds=~/^[a-zA-Z]/) ? $preserveOldIds : $IDPREFIX;
$NEWIDpre="" unless($NEWIDpre=~/^[A-Za-z]/); # allow same old id?

my $IDPREOK=$IDPREFIX;
if(defined $preserveOldIds) {
  if($preserveOldIds and $preserveOldIds=~/^[a-zA-Z]/) { $IDPREOK= $preserveOldIds; }
  else { $preserveOldIds=$IDPREOK=$IDPREFIX; }
} else {
  $preserveOldIds=0;
}  

# data globals
my($mainindex,$ntr)=(0) x 0;
my(%main,%mainsize,%alt,%oclass,%altsize,%maindrops,%altdrops,%didaltdrops,%balt,%drop);  # globals of trclass2maintab
my(%aaqual, %piad, %annotes, %newmain, %altdefer);

my(%mainids, $keepdroph, $nkeepdrop);  
# %dropids, %dropid not used yet; change to culldropkeep table: id, actclass, actinfo, oid?, other..
# orig dropids,dropid is  short list of gff.source tags, e.g. 'cull' in source 'zf17evgcull'

MAIN();

sub MAIN
{
  unless($oname){ $oname=$genes[0] || $seqset[0]; $oname=~s/\.\w+$//; }
  #? openloggit(undef,$oname);
  loggit(0, "BEGIN $0  input=",$oname,"date=",$DATE);
	# do_settings("restore",$trclass);  
	
  my($upokgenes, $upokids)= (0) x 9; 
  my($pubids,$maltab)= map{ "$oname.$_" } qw(pubids mainalt.tab);
  my $generef=(@genes > 0)? \@genes : undef;
  
  ($nkeepdrop,$keepdroph)= ($keepdropin) ? readKeepDrop($keepdropin) : (0); # fill  global %keepdrop
  
  # Separate steps here, call w/ input pubid tab to make pubgff, pubseq sets
  unless( -s $pubids ) {
    ($upokids)= makePubIdtab($pubids,$maltab,$generef);
  }
  if($ONLY_PUBIDS) {
    loggit(0, "done -onlypubids makePubIdtab n=",$upokids);
    exit(0); # loggit?
  }
   
  my $mrnaseq= (@seqset)? $seqset[0] : ""; # many?? all input .mrna .cds and/or .aa
  ## change this mrnaseq to array= grep 'mrna|cdna|xxx' @seqset
  ($upokgenes)= makePublicset($pubids, $oname, $mrnaseq, $generef);
  
  # add tidyup(@publicset,@tmpfiles) ?
}
  
# ---------------------------

## readKeepDrop moved to evigene_pubsets.pm
# my($nkeepdrop,$keepdrop_idhash)= readKeepDrop($keepdropin) if($keepdropin); # fill  global %keepdrop
# use constant { kdDROP => -999, kdCULL => -1, kdOK => 1, kdMUSTKEEP => 2, kdOther => 0 };  
# 
# sub readKeepDrop { # find other subs for this
#   my($intable)=@_;
#   my($nin,$ndup)= (0) x 9; 
#   my(%lkeepdrop); # return for global %keepdrop
#   my($nkdx,%lkeepdropx); # special case dup pubids x oid sets; for global package %keepdropxtab
#   my($ok,$hin)= openRead($intable); 
#   unless($ok) { loggit(1, "#ERR: missing readKeepDrop($intable)"); return 0; }
#   while(<$hin>) {
#     next if(/^\W/); chomp;
#     my($id,@v)=split; 
#     my $act=$v[0];  ## ID here can have _Cn split tag
#     my ($oid,$moid)= ($id,$id); # allow dup pub ids, diff oid in col2, act in col3, need dup id fix for keepdrop{$id}
#     
#     ## upd to for dup pubid table, now use 3 id cols: pubid, seqid, mapid (_C1/2), action, ..
#     ## pass this special keepdrop table to evigene_pubsets.pm, make_pubgff, make_pubseq, ..
#     unless($act =~ /^(ok|keep|drop|cull)/) {
#       my $no=0;
#       if($v[2] =~ /^(ok|keep|drop|cull)/) { # and $v[0]=~/^\w+/ no: $v[0]=~/$IDPREFIX/
#         $no=2; ($oid)= shift @v; ($moid)= shift @v; $act= $v[0];  
#       } elsif($v[1] =~ /^(ok|keep|drop|cull)/) { # and $v[0]=~/^\w+/ no: $v[0]=~/$IDPREFIX/
#         $no=1; ($oid)= shift @v; $moid= $oid; $act= $v[0];  
#       }
#       if($no>0) {
#         $nkdx++;  
#         $lkeepdropx{act}{$moid}=$act;
#         $lkeepdropx{mapid}{$moid}=$id;
#         $lkeepdropx{act}{$oid}=$act unless($lkeepdropx{act}{$oid});
#         $lkeepdropx{seqid}{$oid}=$id unless($lkeepdropx{seqid}{$oid});
#       }
#     }
#     
#     my $kdok= ($act=~/^drop/i) ? kdDROP : ($act=~/^cull/i) ? kdCULL :
#              ($act=~/^(ok|keep)/i) ? kdOK : kdOther; # other action verbs?
#           
#     # NOTE dup id entries allowed, last is active
#     # FIXME: need dup pubid / diff oid keep/drop handling, need oid column (3?) defined in keepdrop table
#     my $kdval= join"\t",$kdok,$oid,@v; # @v has oid now if id,oid,act,.. table
#     my $dupskip=0;
#     if(my $oldv= $lkeepdrop{$id}){ # warn?
#       $dupskip=1 if($kdok ne kdOK and $oldv =~ m/\b(ok|keep)/);
#       $ndup++; } 
#     $nin++;
#     $lkeepdrop{$id} = $kdval unless($dupskip); 
#     $lkeepdrop{$oid}= $kdval if($oid ne $id and not $lkeepdrop{$oid}); #? skip for keepdropxtab
#   } close($hin);
#   
#   if($nkdx>0) {
#     $KEEPDROPX= $nkdx;  %keepdropxtab= %lkeepdropx; #?? bugs?
#     #   for my $s (qw(act mapid seqid)) {
#     #     for my $k (keys %{$lkeepdropx{$s}}) { $keepdropxtab{$s}{$k}= $lkeepdropx{$s}{$k}; }
#     #   }
#   } else {
#     $KEEPDROPX=0; %keepdropxtab= (); 
#   }
#   loggit(0, "readKeepDrop($intable), n=$nin, dupentries=$ndup, keepdropx=$nkdx");
#   return($nin,\%lkeepdrop);
# }

# sub makePubIdtabFromGeneGFF
# add sub makePubIdtabFromTrclass
sub makePubIdtab 
{
  my($pubids,$maltab,$genes)= @_;
  my($ok, $upokids, $ndone, $nmiss)= (0) x 9; 
  
  # FIXME: 18.03/04 got 300 dup pubids in 300k now, mixup recent..
  
	# ($pubid_format,$altid_format)= make_IDPREFIX(); # default abbrev of $organism now
    #^^ wait for input gff IDs as existing id?
  loggit(0, "makePubIdtab($pubids,$maltab)");
  unless(ref($genes)) { loggit(1, "#ERR: missing genes.gff, cant makePubIdtab($pubids)"); return 0; }

  for my $ggf (@$genes) {
    loggit(0, "readGeneIDs($ggf)");
    my $inh;
    if($ggf =~ /^stdin|^-/) { $inh=*STDIN; $ok='stdin'; } 
    else { $ok= open(GIN, $ggf); $inh=*GIN; }
    if($ok) {
      my($nok,$ndrop,$nskip,$firstID)= readGeneIDs($inh,$ggf); # if $ok
      close($inh) unless($ok eq 'stdin'); # unless stdin
      $upokids += $nok;

      if($nok and $firstID and $preserveOldIds and not $pubid_format) {
	      ($pubid_format,$altid_format)= make_IDPREFIX($firstID); # default abbrev of $organism now
	      ## can have new IDPRE here   
        $IDPREOK= $IDPREFIX;
      }
      
    } else {
      loggit(1, "#ERR: cant readGeneIDs($ggf)");
    } 
  }

  reclassGenes(); # link alts to mains via %main,%alt,%balt id hashes
  
  ($pubid_format,$altid_format)= make_IDPREFIX() unless($pubid_format); # default abbrev of $organism now
  $ok= open(OPD,'>',$pubids); my $outpubidh=*OPD;
  $ok= open(OMA,'>',$maltab); my $outh=*OMA;
  $nmiss=$upokids;
  if($ok) {
    ($ndone,$nmiss)= putPubidTab($outh,$outpubidh); # $mainindex,$ntr
    close($outh); close($outpubidh);
  } else {
    loggit(1, "#ERR: cant putPubidTab($pubids)");
  }
  loggit(0, "done makePubIdtab($pubids,$maltab) nin=$upokids, npub=$ndone, nmiss=$nmiss");
  return($upokids);
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
    if(ref($genegffset)) {
      # DUPID fix: may need keepdrop{ oid } handling here to skip dropped dups by oid..
      ($ntra1, $ncdsa1)= make_annot_from("gff", $annotab, $genegffset,  $genenames);
    } elsif(@mrna) {
      ($ntra1, $ncdsa1)= make_annot_from("mrna", $annotab, \@mrna,  $genenames);
    }
  }
   
  
  $FAHDRisValid=1; # export from pubset.pm for annots from mrna hdr
  my %annothash= read_annotab($annotab); #** fill in annothash from pubidh if no table
  
  # may already have .aa, .cds, .mrna from prior steps. is this ok?
  # add input -seqset params ??
  #     our $RNATYPES='mRNA|ncRNA';

	my($pubmrna,$npm,$minfo) = make_pubseq(\@mrna,'mRNA',\%annothash, "$outname.mrna");
	my($pubaa,$npa,$ainfo) 	= make_pubseq(\@aa,'protein',\%annothash, "$outname.aa"); # makename($cdnaseq,'.aa')
	my($pubcds,$npc,$cinfo)	= make_pubseq(\@cds,'CDS',\%annothash, "$outname.cds"); # makename($cdnaseq,'.cds')

      # DUPID fix: may need keepdrop{ oid } handling here to skip dropped dups by oid..
  my($pubgff,$ngenes,$ginfo)= make_pubgff($genegffset,$SRC,\%annothash, $outname);
  
  loggit(0,"publicset: ",$pubmrna,$minfo,$pubaa,$ainfo,$pubcds,$cinfo,$pubgff,$ginfo,$annotab); 
}

sub reclassGenes {
  # my %hasmain;
  my @amain= grep { not $main{$_} } sort keys %alt; # dropmain here now only for SHOWDROPS !
  #?? keepdrop/cull handle here?
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



=item fixme mainofaltid/readGeneIDs problems

 -- split genes, keep as 1 record, not 2
 nDanrer6pEVm000603t1  Danrer6pEVm000603t1_C1  nDanrer6pEVm000603      1       main 
 nDanrer6pEVm081432t1  Danrer6pEVm000603t1_C2  nDanrer6pEVm081432      1       main 
  -- should add mainid=XXX to problem cases, ie byhand.gff alts
  e.g. this locus now has 3 mains:
 nDanrer6pEVm000603t1    Danrer6pEVm000603t1_C1  nDanrer6pEVm000603      1       main   #hobest add
 nDanrer6pEVm081432t1    Danrer6pEVm000603t1_C2  nDanrer6pEVm081432      1       main   #hobest add
 nDanrer6pEVm081438t1    Danrer6pEVm000603t4     nDanrer6pEVm081438      1       main   #map5best

  -- split hassles:
  ** dont drop _C1,2 from id, but get two alts now for 0603t1_C1,2 from hobest, want one.

>> one locus not two
nDanrer6pEVm000126t1    Danrer6pEVm000126t3_C1  nDanrer6pEVm000126      1       main    3346,88p,complete       0       nonote
nDanrer6pEVm081335t1    Danrer6pEVm000126t3_C2  nDanrer6pEVm081335      1       main    3346,88p,complete       0       nonote
  
grep  Danrer6pEVm000603t zfish17m6c_map5beste.pubids
nDanrer6pEVm000603t1    Danrer6pEVm000603t4     nDanrer6pEVm000603      1       main    1948,82p,complete       0       nonote
nDanrer6pEVm000603t2    Danrer6pEVm000603t1_C1  nDanrer6pEVm000603      2       alt     1984,89p,complete       0       nonote
nDanrer6pEVm000603t3    Danrer6pEVm000603t1_C2  nDanrer6pEVm000603      3       alt     1984,89p,complete       0       nonote
nDanrer6pEVm000603t4    Danrer6pEVm000603t2     nDanrer6pEVm000603      4       alt     1962,79p,complete       0       nonote
  
  >> orig 2 loci of 055t in map5best.anntab
grep Danrer6pEVm000055t zfish17m6c_map5best.anntab     
Danrer6pEVm000055t1/same_C2     tridba1a_sNn12l1SRR1524240ridk71Loc7237 main    chr12:47178653:47340608:+       97/92.ci        79/83.ix        93p,4611/4968,4582      catfish16nc:XP_017339435.1,
Danrer6pEVm000055t1d1/Danrer1a_sBn2l2SRR4026141velvk45Loc3455t1 Danrer1a_sBn2l2SRR4026141velvk45Loc3455t1       alt     chr12:47250263:47367232:+       99/99.ci        50.ix   na      na
Danrer6pEVm000055t1d3/Danrertridba1a_sBn2l2SRR1524238idbaidbtk41Loc55430        Danrertridba1a_sBn2l2SRR1524238idbaidbtk41Loc55430      alt     chr12:47247639:47334005:+       100/99.ci       42.ix   na      na
Danrer6pEVm000055t3/Danrertrvelo1a_sBn2l2SRR1524238velvk55Loc1702t1     Danrertrvelo1a_sBn2l2SRR1524238velvk55Loc1702t1 alt     chr12:47250275:47367232:+       99/99.ci        49.ix   62p,3091/4968,3040      catfish16nc:XP_017339435.1,
Danrer6pEVm000055t6/same_C1     Danrer1a_sNn7l1SRR4026141velvk85Loc2302t1       main    chr12:47168344:47174291:-       33/99.ci        2/2.ix  2p,109/4762,196 sacavefish16nc:XP_016323833.1,
  >> update after $id=~ s/_C// adds another 055t1 locus for 055t1_C2
nDanrer6pEVm000055t1    Danrer6pEVm000055t1     nDanrer6pEVm000055      1       noclass 4582,99p,partial3       0       nonote
nDanrer6pEVm092495t1    Danrer6pEVm000055t6     nDanrer6pEVm092495      1       noclass 196,92p,complete        0       nonote
nDanrer6pEVm102106t1    Danrer6pEVm000055t1_C2  nDanrer6pEVm102106      1       main    noaa    0       nonote
nDanrer6pEVm102106t2    Danrer6pEVm000055t1d1   nDanrer6pEVm102106      2       alt     3076,87p,complete       0       nonote
nDanrer6pEVm102106t3    Danrer6pEVm000055t3     nDanrer6pEVm102106      3       alt     3040,87p,complete       0       nonote
nDanrer6pEVm102106t4    Danrer6pEVm000055t1d3   nDanrer6pEVm102106      4       alt     2435,98p,partial3       0       nonote

=item mainofaltid probs

  Danrer6pEVm079420t3_C2	noclass	
  Danrer6pEVm079420t3	NOMAIN	Danrer6pEVm079420t1/alt,Danrer6pEVm079420t4/alt
  ..
  ID=Danrer6pEVm005545t31_G3t2;oldaltid=Danrer6pEVm005545t2_C1
  Danrer6pEVm005545t31_G3	noclass	
  Danrer6pEVm005545t31_G3t1	NOMAIN	Danrer6pEVm005545t2_C1/alt

=item gannot pubid annots

  * moved to evigene_pubsets gene_annot_brief
  reproduce this annotation (trclass2mainalt.pl) in notes hash
    aaref:5767,dapsim:Dapsim1EVm000004t1,chrmap:100a,98i,25555l,33x,29sc:324762-356329:+,pflag:0
    
  my @ANKEY_MRNA= qw(aalen  cov pid nexon clen namealn Dbxref scoresum);
  leave out special annots of altbest, alttr, altid/mainid/newlocusid

=cut

## gannot moved to evigene_pubsets gene_annot_brief($id,@mrnarow)
# sub gannot {
#   my($id, @gff)=@_;
#   my $at= $gff[8];
#   my $stype= $gff[2]; # mRNA|ncRNA.. add to tblinfo? to an
#   
#   #? also make tblinfo for ann.txt table, now reads from mrna.seq hdr ..
#   # tblinfo: $pubid,$oid,$trlen,$cdsoff,$aaqual,$annogaps,$dbxref,$namepct,$gname,$cddname
#   my $tblinfo= parse_evgheader($id,$at,0); # this works right; parses evg.gff or mrna attr
#   
#   my @mapq= map{ my($v)= $at=~m/\b$_=([^;\s]+)/?$1:0; int($v); } qw(cov pid clen nexon);
#   
#   my($ndx,$nal)=(0,0);
#   if($tblinfo) {
#     $ndx= $$tblinfo{'dbxref'}||0; $ndx=~s/,$//;
#     $nal= $$tblinfo{'namealn'}||0;
#   } else {
#     ($ndx)= $at =~ m/Dbxref=([^;\s]+)/?$1:0; $ndx=~s/,$//;
#     if( my($nap)= $at =~ m/namealn=([^;\s]+)/) {  # namealn=58p,17578/30278,17622;
#       ($nal)= $nap=~m/^\d+..(\d+)/; unless($nal) { ($nal)= $nap=~m/(\d+)/; }
#     }
#   }
#   
#   my $an= ($nal>0)?"aaref:$nal,$ndx,":"0,0,";
#   my $cm= sprintf "chrmap:%da,%di,%dl,%dx,%s:%d-%d:%s", @mapq, @gff[0,3,4,6];
#   unless($stype=~/mRNA/){ 
#     $cm .= ",mol:$stype"; #? stype= seqt= rnat= rtype= mol= ?
#     $tblinfo->{'seqtype'}= $stype;
#     } 
#   $an .= $cm;
#   
#   return ($an, $tblinfo);
# }

=item readGeneIDs

  read gene class info + annots from GFF of evigene processed gene sets,
  i.e., main + alternate genes processed for map overlaps, best locus representatives
  
=cut


use constant KEEP_C => 0; # split ID _C[12..n] tag ; this is problematic, have keep_C1, cull_C2

sub readGeneIDs {
  my($inh,$infile)= @_; #? or while(<>)

    ## add cCULL, cDROP? cKEEP => main or alt
  use constant { cNone=> 0, cMAIN => 1, cALT => 2, cNEWLOC => 3, cCULL => 4, cDROP => 5, cOther => 6 };
  
  sub mainofaltid { my($aid,$ismain)=@_; 
    unless($ismain) { 
      my($t)= $aid=~s/(_[CG]\d+)//?$1:""; 
      $aid=~s/t\d+$//; 
      $aid.="t1" unless($aid=~/t\d+$/);  
      if(KEEP_C or $t=~/_G/){ $aid.=$t; } } 
    return($aid); 
  }
  
  my($src,$ok,$inid,$nok,$ndrop,$nskip,$firstID)=(0) x 9;
  my(%did); # FIXME need global over all in.gff >> use %balt{id} global
  my $OTHERCLASS= cOther; # == cCULL for CULLOTHERS
  
  while(<$inh>) {
    if(/^\W/){
      # ? collect scoresum of '#x.mRNA' of bestmain.gff ?
      # %mainids has to be global, and read main.gff before alts.gff
      if(/^#x./ and /\t($RNATYPES)/) { 
        my($id)= m/\bID=([^;\s]+)/; 
        ## KEEP_C here ?? 
        my($ovd)=m/\boverids=([^;\s]+)/g; 
        if($ovd){ $mainids{$id}= (split",",$ovd)[0]; }
        # my($gscore)= m/\bscoresum=([^;\s]+)/g; # gscore or scoresum or col5
        # my($skip)= m/\bskip=([^;\s]+)/g; # gscore or scoresum or col5
        #my @an= m/((?:scoresum|overids|skip)=[^;\s]+)/g; 
        } 
      next;
    }

    if(/\t($RNATYPES)/) { # mRNA
    $ok=1;
    my @gff=split"\t";
    my $gsrc= $gff[1];
    # gannot(): my $gtype= $gff[2]; # mRNA | ncRNA .. preserve this in pubid tab?
    my $gann= $gff[8];
    
    # FIXME: add cull/skip class, ($class)= $gsrc =~ m/(alt|cull)$/ ;
    # culls may have been main or alt, ie have newaltid
    # use separate input keep_drop_cull.table ?
    
    $src= $SRC || $gsrc; #? use input gff[1] col unless defined? not used here now?
    my($id)= m/\bID=([^;\s]+)/;  $inid= $id;
    #o.my($gpre,$gnum,$gd,$ti,$isplit)= evigene_idparts($id);
    my($gpre,$gnum,$ti,$isplit,$gdup)= evgIdParts($id);
    
    # should add mainid=XXX to problem cases, ie byhand.gff alts
    my($mainid)=m/\bmainid=(\w+)/?$1:0; # trust this annot? require isalt also?
    unless($mainid) { $mainid= $mainids{$id} || 0; }
    my($oad)=m/\boldaltid=(\w+)/?$1:0; my($nad)=m/\bnewaltid=(\w+)/?$1:0;
    my($old)=m/\boldlocus=(\w+)/?$1:0; my($nld)=m/\bnewlocus=(\w+)/?$1:0;
    my($scoresum)= (m/scoresum=([\d-]+)/)?$1:0;
    my($aq)= (m/aalen=([^;\s]+)/)?$1:0;
    my $aasize= ($aq =~ m/(\d+).(\d*)/) ? "$1.$2" : 0; # size.pCDS for best sort
  
    unless(KEEP_C){ $id=~s/_C\d// if($isplit); $mainid=~s/_C\d$//; } # want this?

    # dang _C1,2 parts separate? ** FIXME, need to see both _C1,_C2
    if($did{$inid}) { $ok=0; next; } 
    elsif($balt{$id}) {   #  or $did{$id}
      # BUG below? getting dup main entries when shouldn't, from extra in.gff?
      if($inid eq $id) { $ok= 0; }
      else { if(my $ocl= $oclass{$id}) { $ok=($ocl=~/^cull/)?1:0; } }
      next if($ok == 0); 
    } 
    
    #* may need annot in pubid table for _C/_G edits, exists : oid retains?    
## ugh altbest id t1_C2 bug : need to cut out 'tNNtNN' 2nd alt tag?
# Danrer6pEVm022697t4_C2	noclass	
# Danrer6pEVm022697t4t2_C2t1	NOMAIN	Danrer6pEVm022697t4t2_C2/alt
# Danrer6pEVm022697t4t3_C2t1	NOMAIN	Danrer6pEVm022697t4t3_C2/alt
# grep ID=Danrer6pEVm022697t chr21gn.malt.altbest.gff | cut -f9 | sed 's/aamap=.*//; '   
# ID=Danrer6pEVm022697t4t2_C2;oldaltid=Danrer6pEVm022697t2_C2;trg=Danrer6pEVm022697t2 67 1018;
# ID=Danrer6pEVm022697t4t3_C2;oldaltid=Danrer6pEVm022697t5_C2;trg=Danrer6pEVm022697t5 78 650;
# ID=Danrer6pEVm022697t4t4_C2;oldaltid=Danrer6pEVm022697t7;trg=Danrer6pEVm022697t7 91 495;
# ID=Danrer6pEVm022697t4t5_C2;oldaltid=Danrer6pEVm022697t9_C2;trg=Danrer6pEVm022697t9 98 426;

=item split bug 180330

 grep '^vDanrer6pEVm004866t1     ' chreset-map6na.pubids | head
 vDanrer6pEVm004866t1	Danrer6pEVm005074t4	vDanrer6pEVm004866	1	main	788,86p,complete	0	aaref:726,zfish16nc:NP_571403.1,chrmap:87a,99i,2733l,11x,chr20:54192943-54206663:-	Danrer1a_sRn1l1SRR4994225velvk79Loc1037t4,Danrer6pEVm005074t4
 vDanrer6pEVm004866t1	Danrer6pEVm005074t152	vDanrer6pEVm004866	1	main	725,74p,complete	0	0,0,chrmap:99a,99i,2930l,12x,chr20:54192871-54222720:-	tridba1a_sNn12l1SRR1524240ridk81Loc231231,Danrer6pEVm005074t152

 grep '^Danrer6pEVm005074t4_' $pt.map6k.keepdrop.gscore     
 Danrer6pEVm005074t4_C2	keep	gscore=4267	vin=10	vho=726	pho=100	altpar=1	aalen=788	cov=87	acov=95	lspan=1.14	locus=chr20:54192943-54206663:-	gclass=alt4	flags=split,evuncull:h0/i1	gsold=4095
 Danrer6pEVm005074t4_C1	cull	gscore=368	vin=0	vho=726	pho=100	altpar=1	aalen=788	cov=10	acov=4	lspan=0	locus=chr20:54190464-54192154:-	gclass=main1	flags=cull3.ovEVGm001440t9/99.96,split	gsold=368

 chreset-map6na_pub.gff
 chr20	zf17evgm6v	mRNA	54192871	54222720	4409	-	.	ID=vDanrer6pEVm004866t1;trg=Danrer6pEVm005074t152 1 2925;cov=99.8;nexon=12;pid=99.4;clen=2930;offs=544-2721;oid=Danrer6pEVm005074t152,tridba1a_sNn12l1SRR1524240ridk81Loc231231;aalen=725,74p,complete;cdsoff=544-2721;Name=Uncharacterized protein;gscore=4335;vin=9;vho=725;altpar=1;acov=100;cxlen=2178/2925,74%;inexon=12/12/10;scoresum=4409
 chr20	zf17evgm6v	mRNA	54190464	54192154	413	-	.	ID=vDanrer6pEVm004866t1_C1;Split=1;trg=Danrer6pEVm005074t4 1 278;aamap=38;cov=10.2;nexon=2;pid=98.9;path=1/6;chim2=chr20:54192943-54206663:-;Name=Heat shock protein HSP 90-alpha 1;clen=2733;offs=167-2533;oid=Danrer6pEVm005074t4_C1,Danrer6pEVm005074t152,Danrer1a_sRn1l1SRR4994225velvk79Loc1037t4;aalen=788,86p,complete;cdsoff=167-2533;namealn=100p,726/726,788;Dbxref=zfish16nc:NP_571403.1,;gscore=368;vin=0;vho=726;altpar=1;acov=4;flags=split,;cxlen=112/277,40%;inexon=2/2/1;scoresum=413
      ^^ should be cull
 chr20	zf17evgm6v	mRNA	54192943	54206663	4095	-	.	ID=vDanrer6pEVm004866t1_C2;Split=2;trg=Danrer6pEVm005074t4 279 2665;aamap=750;cov=87.3;nexon=11;pid=99.2;path=2/6;chim1=chr20:54190464-54192154:-;Name=Heat shock protein HSP 90-alpha 1;clen=2733;offs=167-2533;oid=Danrer6pEVm005074t4_C2,Danrer6pEVm005074t152,Danrer1a_sRn1l1SRR4994225velvk79Loc1037t4;aalen=788,86p,complete;cdsoff=167-2533;namealn=100p,726/726,788;Dbxref=zfish16nc:NP_571403.1,;gscore=4095;vin=10;vho=726;altpar=1;acov=95;flags=split,
      ^^ should be vDanrer6pEVm004866t2 or other alt
=cut

    my $cl="none";
    my $isalt=(/;alttr=(\d+)/)?$1:0; # ;isalt=1 ?
    my $isother=0; # NOT isalt or ismain, ie iscull + other
    my $mclass = cNone; 
    my $md= $mainid || $id; ## $gd; # evigene_idOfparts($gpre,$gnum,$gd,1,0); #? use $gd? drop ti?
    
    my($gannot,$tblinfo)= gene_annot_brief($id, @gff); # was gannot()

    ## opt to reparse of last round of pub,cull.gff : no oad/scoresum/..
    if($CLASSbySOURCE) {
      my $gcla= ($gsrc=~m/(alt|cull|nc)$/)?$1:"";
      my $iscull= ($gcla eq "cull")?1:0;
      $md=  mainofaltid($id); # $mainid || : outof date mainid= tag skip or update
      $isalt= ($md ne $id or $mainid)? 1 : 0; # for cullalt
      $mclass= ($iscull)? cCULL # need cullalt also
          : ($gcla eq "alt")? cALT 
          : ($gcla eq "nc") ? cMAIN 
          : ($md eq $id) ? cMAIN 
          : cOther;
    }      
    elsif($oad or $nad) { $mclass = cALT; $md= $mainid || mainofaltid($nad || $id); } # EVm000t2t3 hack format, or EVm00t3
    elsif($old or $nld) { $mclass= cNEWLOC; $md= $mainid || mainofaltid($nld || $id, not $isalt); }
    elsif($scoresum) { $mclass= cMAIN; $md= $mainid || $id; } # or $ti < 2 ?? cant assume?
    else {
      # my ($INS)= $infile =~ m/(\w+).gff/; #? file name not good for this
      # $isalt=1 if($INS =~ /altbest/); $isnewloc=1 if($INS =~ /newloci/);
      ## mainofaltid() guess doesnt work well enough; add mainid=XXX annot; use overids=?
      if($mainid) {
        if(my $mcl= $main{$mainid}) {
          $isalt=2; # guess this .. for keep.other group w/o altbest info
        } else { # check if mainid now is alt?
          if(my $mmd= $balt{$mainid}) {
            if(my $acl= $alt{$mmd}{$mainid}) {
              $isalt=2; $mainid= $mmd; # maybe.
            }
          }
        }
      }
      if($isalt) { $mclass = cALT; $md= $mainid || mainofaltid($id); } 
      else { 
        $mclass= cOther; #or $OTHERCLASS; #was $mclass= cMAIN;  # FIXME: option for -otherclass=cull|xxx|remainder..
        $md= $mainid || $id; 
        }
    }
    # if($mclass == cALT and $id eq $mainid) { $mclass= cMAIN; } #?? check here, abive?
    
    # if($iscull) { $mclass= cCULL; } #? like drop? but output all data w/ new ids, separate?
    ## KEEP_C=0 bug : keepdrop id can have _C tag
    ## FIXME: are loosing some _C keep due to _C1 keep, _C2 cull ;  ie KEEP_C=0 bug
    if($keepdropin and (my $kdv= $keepdroph->{$inid} || $keepdroph->{$id})) {
      ## add dup pubid check ; add dupid flag in kdv?
      ## may need to keep row in pubid tab for dup id cases, even drops, so dup data are handled right
      ## but then need uniqof(dupid) ie. {pubid}d2,d3..
      my $oid= $tblinfo->{oid}; 
      for my $od (split",",$oid){ if(my $kdvo= $keepdroph->{$od}){ $kdv=$kdvo; last; } }
      
      my($kdi,$kact)= split" ",$kdv; #no need kdi?
      if($kdi < 0) {   # drop or cull
        if($kdi <= kdDROP){ $mclass=cDROP; $ok=0; $ndrop++; } #  $kact=~/drop/ else drop anything else?
        elsif($kdi == kdCULL){ $mclass= cCULL; } #  $kact=~/cull/
      } elsif($kdi>0) { 
        $mclass= cMAIN unless($mclass == cALT || $mclass == cNEWLOC); ## cKEEP ? treat like cMAIN if no other
      } # keep/ok
    }
    # if($dropids) { # fixme for map5best too liberal gscore, remove no(ho+in)evidence loci
    #   if(my $ds=$dropid{$id}) { if($gsrc =~ m/$ds/){ $ok=0; $ndrop++; } }
    # }
    
    if($mclass == cALT) { 
      if($oad and $did{$oad}) { $ok= 0; }
      elsif($nad and $did{$id}) { $ok= 0; }
      elsif($oad) { $id= $oad; } #?? need this to stop tNNtNN hack bugs?
      
      ## this id change is problem: tNtN; when main is t2, alt is t1, altbest assumes main new id t1..
      ## ID=Danrer6pEVm004586t2t7; oldaltid=Danrer6pEVm004586t1;
      $isalt=1;
      $cl="alt"; # variants?
      $src .= "alt";
    	# $alt{$md}{$td}= $cl; $balt{$td}=$md; $altsize{$td}= $aasize; 
    
    } elsif($mclass == cMAIN) { 
      $isalt=0;
      $cl="main"; # variants?
      # $main{$md}=$cl; $balt{$td}=$md;   $mainsize{$td}= $aasize;  
 
    } elsif($mclass == cNEWLOC) { 
      $ok= ($id=~/^newloc/ or $nld =~ /^newloc/)?1:0; #? others? what are ids
      if($old and $did{$old}) { $ok= 0; }
      elsif($nld and $did{$id}) { $ok= 0; }
      #X elsif($old) { $id= $old; } #??
      # $isalt=(/;alttr=1/)?1:0; # see above
      $src .= "alt" if($isalt);
      $cl=($isalt)?"altnew":"mainnew"; # variants?

    } elsif($mclass == cDROP) {  
      $ok=0; $nskip++;
    
    } elsif($mclass == cCULL or $CULLOTHERS) {  # CULLOTHERS = opt to treat cNone/Cother as cCULL ?
      #? $isalt=0; #? maybe yes/no
      # $isalt=2; #? $md= $id; # isolate is noalt/nomain locus?
      #** cull fails unless isalt ? need some check of main{id} == cull class
      $cl= ($isalt or $mainid)?"cullalt":"cull"; # variants?
      # $cl="cull";  
      $src .= "cull";
      $mclass= cCULL;
      $isother=1; # separate from noncull alt/main
      $ok=1;
      
    } else {  
      # change old ID t2..N to t1 ?
      $ok=0; $nskip++;
    }
    
    # 0330: still bugs with dang splits, ie. cull _C1, keep _C2 as alt
    #  .. then pubids dup new ids, n=47, one is _C2 alt, other is main of _C2 alt
    if($ok and $did{$inid}) {  
      $ok= 0; # unless $dids eq $INS for splits ? no splits same ID ?
    }
     
    if($ok) {
      $nok++;
      # ugh ID changes from oldaltid are problem..
      # $aaqual{$id}= $aq; #  for pubtable
      # $notes{$id}= $gannot; #  for pubtable, add annots
      # $piad{id}= $piad; #  for pubtable, align to main mrna: ident%,aln% == 99,80
      $annotes{$id}{aaqual}= $aq;
      $annotes{$id}{piad}= 0;
      $annotes{$id}{notes}= $gannot;
      $annotes{$id}{oid}= $tblinfo->{oid};
      
      $mainsize{$id}= $altsize{$id}= $aasize; #? dont need both?  aasize{$id}= $aasize;
      $oclass{$id}= $cl; # use for all classes?
      if($isother) {
    	## oclass adds problems, treat like %main? add to 
    	$oclass{$id}= $cl; $balt{$id}=$md; $main{$md}="main";
    	#?? maybe need:  $alt{$md}{$id}= $cl if($id ne $md);
      } elsif($isalt) {
    	$alt{$md}{$id}= $cl; $balt{$id}=$md; 
      } else {
      # FIXME cullalts here cause problems (=> cull main)
      $firstID= $id unless($firstID); # but not _C/G suffix
      $main{$md}=$cl; $balt{$id}=$md; # dang, should be md == id 
      }
       
      $did{$id} = $did{$inid} = $src;
    }
    
  }
  
  # elsif(/\texon\t/){
  #
  # } 
  
  }  
  
  return($nok,$ndrop,$nskip,$firstID); # $nmain,$nalt,...
}


sub putPubidTab {
  my($outh, $outpubidh)= @_; # globals?
 
  # uses set globals of trclass2maintab
  # my(%main,%mainsize,%alt,%oclass,%altsize,%maindrops,%altdrops,%didaltdrops,%balt,%drop);  # globals of trclass2maintab
 
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


#upd1708: preserveOldIds option ..  keep input gene ids as best feasible, including altnum?
#  need all alts per main, check evg gene ids for 1 x many, also need check across loci for 2+ new loci w/ old geneids
#  where cant preserve geneid for all of locus, do what? new geneid, or modify old?
## moved to evigene_pubsets.pm:evgIdParts()
# sub evigene_idparts {  # evgIdParts($id)
#   my($id)=@_;
#   my($gpre,$gnum,$gd,$ti,$gdup,$isplit)=(0) x 9;
#   $gd=$id; 
#   ## FIXME for _G2 _C2? tags, G2,n is diff locus id
#   ($isplit)= ($gd=~s/_(C\d+)$//)?$1:0;
#   ($gdup)  = ($gd=~s/_(G\d+)$//)?$1:0;
#   
#   # bug: Danrer6pEVm003029t41/main Danrer6pEVm003029t41t2, Danrer6pEVm003029t41t3 altids
#   if($gd =~ m/^(\w+[A-Za-su-z])(\d\d+)t(\d+)$/) { 
#     ($gpre,$gnum,$ti)=($1,$2,$3);
#     $gd= $gpre.$gnum;
#   } else {
#     $gd =~ s/[a-su-z]\d+$//; # non-t suffix ??? ie. t2d33
#     ($ti) = ($gd=~s/t(\d+)$//)?$1:0;  # BUG: m00001t12t3 hack format
#     if( $gd =~m/^(\w+[A-Za-su-z])(\d\d+)/) { ($gpre,$gnum)=($1,$2); $gd=$gpre.$gnum; }
#     else { # shouldnt be here?
#       $gd =~ s/t\d+$//; # extra tNtN..
#       ($gnum) = ($gd=~m/(\d+)$/) ? $1:0;  # BUG.. ? set gnum=0 meaning need new id?
#       ($gpre=$gd) =~ s/$gnum$//;
#       $gnum=0; # BUG fix??
#     }
#   }
#   
#   if($gdup) { $gd.="_$gdup"; $gnum=0; } # gnum invalid?
#   return($gpre,$gnum,$gd,$ti,$isplit);
#   #OR:   return($gpre,$gnum,$ti,$isplit,$gdup); #?
# }

# sub evigene_idOfparts {  # evgIdOfParts(@parts)
#   my($gpre,$gnum,$gd,$ti,$isplit)= @_; # from evigene_idparts
#   #n my($gpre,$gnum,$ti,$isplit,$gdup)= @_; # new evigene_idparts
#   $gpre ||= $IDPREFIX;
#   $ti ||=1;
#   $gnum= sprintf( '%06d', int($gnum)); # leading 000
#   my $id= $gpre . $gnum . "t$ti"; 
#   $id.="_C$isplit" if($isplit); 
#   # $id.="_G$gdup" if($gdup);
#   return($id);
# }

sub preserveOldIds {
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


# sub OLD_make_pubid # moved to evigene_pubsets.pm
# {
#   my($oid, $pubidnum, $altnum, $preservedIds)= @_;
#   $pubidnum_start= $pubidnum; #?
#   my $alti=$altnum;
# 
#   my($pubid,$pubgene);
#   if($preservedIds and ref($preservedIds)) {
#     if(my $pid= $preservedIds->{$oid}) {
#       ($pubid,$pubgene)=($pid,$pid); #? yes? no?
#       if(my($pgene,$palt)= $pid=~m/(\w+)t(\d+)$/) {
#         $pubgene= $pgene;
#         $alti= $altnum; # $palt; #? preserve altnum? or caller's altnum 
#         # if($palt ne $altnum) ..
#         $pubid   = $pubgene . sprintf( $altid_format, $alti);  
#         return($pubid,$pubgene,$altnum);
#       }
#     }
#   }
#   $pubgene = sprintf( $pubid_format, $pubidnum); 
#   $pubid   = $pubgene . sprintf( $altid_format, $alti);
#   return($pubid,$pubgene,$alti);
# }

# sub OLD_make_IDPREFIX
# {
#   my($existingID)= @_;
#   my $digits=6; # or 7?
#   my $ChangeDefaultIdPrefix= 0; # ($IDPREFIX eq $DEFAULT_SETTINGS{'IDPREFIX'}) ? 1:0;
#   
#   if($existingID and $existingID !~ m/^$IDPREFIX/) {
#     $existingID=~ s/t\d+$//;
#     my($prefix,$nums)= $existingID =~ m/^(\w+)(\d\d\d+)$/; ## prefix may have numbers, how many? 
#     if($prefix) { $IDPREFIX= $prefix; $digits= length($nums) || 6; }
#     $ChangeDefaultIdPrefix=0;
#   }
# 
#   my $nd= ( $IDPREFIX =~ s/(0\d+)$// ) ? length($1) : $digits;
#   my $pubid_format = $IDPREFIX.'%0'.$nd.'d';  
#   my $altid_format = 't%d'; 
#   return($pubid_format,$altid_format);  
# }

1;

__END__
