#!/usr/bin/env perl
# overgenedup.pl

=item about

  overgenedup.pl -act [drop|keep|mark] -over refgenes.gff -in geneset{1,2,3}.gff  > output.gff

  simple filter to drop|keep|mark genes that are identical to other input genes
  identity = same location/orientation for mRNA+exon+CDS

  ** ASSUMES gene model is mRNA > exon > CDS
  ** ASSUMES input gff is ordered by gene records (mRNA/exon,CDS all together per ID)

=item updates
  
  revised to work like overlapfilter -in tomark.gff -over overlap.gff -act markid -mark XXX

  now has fairly accurate _similargene method, for identifying equivalent
  gene transcript models at same locus.  Equivalence is from exon and CDS_exon
  overlaps (same strand), with numeric result for percent equivalence of both:
      I100    = perfect identity (+/- slopexon=N)
      C95.75  = SAME_CDS at >= 90% identity, 75% exon equivalence.
      50.30   = 50% CDS equivalence, 30% all exon equivalence.

  Primary output is the -input geneset.gff with annotations added to mRNA rows
  of equivalent genes ID/score
  
=item mem overloading from large genes.gff

  now reads all genes.gff, same data for overlaps if self overgenedup,
  easily overloaded w/ large gene sets,
  add -sortedgff opt to say input.gff is sorted by location/chr
  .. process per chr/scaffold, also reading overlaps.gff as needed
  
=item equalgene or similargene usage

  commonest use is to tabulate equivalences in two gene.gff files
  add below tabulation filter, and symlink this with equalgene.pl:
    equalgene.pl -in $ingff -over $ovgff > $inam.over$onam.tab
  
=item similargene usage

  #1
  set ingff=aphid2pub8f.gff ; set inam=evigene8f
  set ovgff=../acyr1-ACYPImRNA.anmap1.gff ; set onam=acypi1
  set ovflags=""
  
  #2
  set ovgff=aphid2pub8f.gff ; set onam=evigene8f
  set ingff=../acyr1-ACYPImRNA.anmap1.gff ; set inam=acypi1
  set ovflags=""
  
  #3
  set ingff=aphid2pub8f.gff ; set inam=evigene8f
  set ovgff=stdin ; set onam=ncbigene2
  set ovflags="-wrongmrna"
  cat ../acyr2_gnomonxref.gff ../acyr2_ncbirefgene.gff | sed 's/^#a.//;' | grep -v '^#' |\
  
  $evigene/scripts/overgenedup.pl -in $ingff -over $ovgff $ovflags \
  -slopexon=8 -type similarCDS -mincds=10 -minutr=33 -act markid -mark overg | \
  grep mRNA | perl -ne '($g,$ov)= m/(?:ID|overg)=([^;\s]+)/g;  $ov||="na"; \
  $og= (m/oid=([^;\s]+)/) ?$1:"noid"; \
  ($r,$b,$e)=(split"\t")[0,3,4]; $r=~s/^Scaffold(\S+)/${1}sc/i; \
  print join("\t",$g,$og,$ov,"$r:$b-$e"),"\n";' \
  > equal/$inam.over$onam.tab5

=item similargene equal.ids table

  Reciprocal X.overY.tables are created above. 
  Gene1ID                 Gene1AltId              Gene2ID/equalscore     Location
  acyp2eg0000001t1        AUGepi4p1s1g2t1         na      1sc:35326-36818
  acyp2eg0000002t1        aphid_cuf8r27Gsc1.248.1 XM_003239978.1/C100.83  1sc:72526-88761
  acyp2eg0000002t2        aphid_cuf8r27Gsc1.248.1t2       XM_003239978.1/4.83     1sc:71599-88761
  acyp2eg0000006t1        aphid2Trinity34495p1    XM_003239980.1/C99.81   1sc:153906-158411
  acyp2eg0000007t1        AUGepir10p1s1g14t1      XM_001945736.2/C99.70,XM_003239981.1/C99.65     1sc:159336-162347

  NM_001126134.2  noid    acyp2eg0000191t1/C99.85,acyp2eg0000191t2/C98.76,acyp2eg0000191t4/80.77,acyp2eg0000191t3/74.83   1sc:2365517-2373872
  NM_001161963.1  noid    acyp2eg0000010t1/C99.49 1sc:174502-177800
  NM_001162137.1  noid    acyp2eg0000173t1/C98.83 1sc:1912915-1917649

  
  This perl combines such to one final table of XY.equal.ids
  
  cat ncbigene2.overevigene8f.tab5 evigene8f.overncbigene2.tab5 | cut -f1,3,4 | cat ../ncbiref2.altids - | \
  env LEF=acyp2eg perl -ne 'chomp; if(/^alt\t(\S+)/) { $alt{$1}=1; next;} \
  ($d,$e,$loc)=split"\t"; @e=split",",$e; @v1=(); \
  foreach (@e) { ($_,$v)=split"/"; push @v1, $v||"."; }  \
  $eg=$e1=shift @e; $v1=shift @v1; $eg =~ s,t\d+|\-R\w,,;  @v=@f=(); \
  for $i (0..$#e) { $e=$e[$i]; unless($e=~/$eg/ or $alt{$e}) { push @f,$e; push @v,$v1[$i]; } } \
  $d1=$d; @e=($e1,@f); @v=($v1,@v); for $i (0..$#e) { $e=$e[$i]; $v=$v[$i]; $d=$d1; \
  unless($d =~ /$left/) { ($d,$e)=($e,$d); } $de="$d\t$e"; \
  $dlo{$de}=$loc; $dv{$de}=$v unless($dv{$de} and $v eq "."); } \
  BEGIN{ $left=$ENV{LEF}; } END { @de= sort keys %dv; foreach $d (@de) { $lo=$dlo{$d}; $v=$dv{$d} || "."; \
  print "$d\t$v\t$lo\n"; }  }'\
  > ncbigene2-evigene8f.equal.ids


  
=item author
  
  don gilbert, gilbertd near indiana edu, ca 2010-2012
  part of EvidentialGene, evigene/scripts/

=cut

use strict;
use warnings;
use Getopt::Long;

use constant VERSION => '2017.04.01'; #upd17

use constant { ACT_DROP=>1, ACT_KEEP=>2, ACT_MARK=>3, ACT_MARK_WITH_ID=>4,
      ACT_MARK_OVERBASE=> 5, ACT_MARK_OVERBASEID=> 6, ACT_NULL => 9};
use constant { kOVERLAP=>1, kSAMELOC=>2, kNEARLOC=>3, kCUTOVER=>4 };

# revise to work like overlapfilter -in tomark.gff -over overlap.gff
use constant COLLECT_OVER => 1;

my $UTRSLOP = 49; # <=  when samecds, UTR difference < this is same gene
  my $UTRSLOPN = 0; # UTR exon count diff
my $TINYCDS = 40; # ignore tiny CDS tests for typeover == CDS
my $EXONSLOP = 0;
my $SAME_CDS= 90; # ** option
my $SHOW_CDSUTR= 0; # similar only ?
my $MINID_CDS= 0;
my $MINID_UTR= 0;
my $EQUAL_CUT= 66; # ** option for summary

my $SELFKEEP= 0; # - option
my $MRNA_SPAN_WRONG = 0; # FIXME;
my $BINSIZE   = 5000 ; 
my $debug=0;
my $pctover= 0;

my $SKIPREF=$ENV{skipref} || "NOPATH"; # new 2013 gmap2gff fake mRNA locations for unmapped
my ($overlaps,$dosum,$passtypes,@input,$itype,$action,$actid,$typeover,$ok,$mark,$nin);
my $mrnatypes='mRNA'; # FIXME: allow ncRNA|... others
my $exontypes='exon|CDS';

# collect_overlaps globals; see drop_overlaps(), addgene()
my $overlaplist= {}; my @overgenes= (); my %overlapids=();  my $overgeneid= 0; my $lastgeneid= 0;
  
my (%sums,%sumid);
my $IgnoreOutOfOrder=0;
my $orientstrict=0;
my $nostrand=0;
my $EqualGene=0;
my $NostrandOneExon=0; # 17.03 make default maybe: -oneexonstrandless vs -nooneexonstrandless
my $THISSIZE=0; # 16.02, revise over/pCDS.pEXON relative to this-gene, not max(this,overgene), sub _similargene
my $addAltsNotOver= 0; # 17.03 self equalgene, add altsNotOver column
my $SORTEDGFF= 0; # upd 17.03
my $SELFSAME= 0; # upd 17.03

$overlaps="";
# $action="drop";

my $optok= GetOptions(
  "overlaps=s", \$overlaps, 
  "input=s", \@input,  
  "typeover=s", \$typeover, 
  "exontypes=s", \$exontypes, 
  "mrnatype=s", \$mrnatypes, 
  ##"pctover=i", \$pctover, #? use this or not
  "action=s", \$action, 
  "mark=s", \$mark, 
  "passtypes=s", \$passtypes,  #??
  "slopexon=i", \$EXONSLOP,  
  "sloputr=i", \$UTRSLOP,  #??
  "samecds=i", \$SAME_CDS,  
  "mincds=i", \$MINID_CDS,  
  "minutr=i", \$MINID_UTR,  
  "wrongmrnaspan!", \$MRNA_SPAN_WRONG,
  "equalgene!", \$EqualGene,
  "nostrand", \$nostrand,
  "oneexonstrandless!", \$NostrandOneExon,
  "strictstrand|orientstrict!", \$orientstrict,  # change to general strict/loose flag
  "IgnoreOutOfOrder!", \$IgnoreOutOfOrder,
  "SELFSAME!", \$SELFSAME, # overlaps == input
  "keepself!", \$SELFKEEP,
  #old#"SELFKEEP!", \$SELFKEEP,
  "THISSIZE!", \$THISSIZE,
  "SORTEDGFF!", \$SORTEDGFF,
  "SKIPREF=s", \$SKIPREF,
  "summary:s", \$dosum,
  "debug!", \$debug, 
  );

# upd17: option to say overlaps == input for  efficiency : reuse -SELF? -THIS?; 
#  -sortedgff opt to process per chr/scaff, not reading all data before

if($EqualGene or $0 =~ /geneequal|equalgene/) {
  $EqualGene=1;
  $action="markid";  $mark="overg";
  $typeover= "similarCDS" unless($typeover);
  $addAltsNotOver= 1 if($THISSIZE); #? need own opt, reuse -self opt ???
  $EXONSLOP= 8 unless($EXONSLOP);
  $MINID_CDS=10 unless($MINID_CDS); # higher?
  $MINID_UTR=33 unless($MINID_UTR);
}

push @input, @ARGV;
push @input, "stdin" unless(@input);
#  $SELFSAME upd 17.03
if($SELFSAME) { if($overlaps) {} else { $overlaps=$input[0]; } } #?? 

die "usage:
  overgenedup.pl -act [drop|keep|mark|markid] -over refgene.gff -in geneset{1,2,3}.gff > ingenes.over.gff
   simple filter to drop|keep|mark -input genes.gff that are identical or similar to -overlap genes
OR
  equalgene.pl -in \$iname.gff -over \$oname.gff > \$iname.over\$oname.tab
  equalgene.pl -selfsame -in \$iname.gff > \$iname.overself.tab
  
  options
    -equalgene  : output table of equivalences of -in to -over
    -typeover=$exontypes  : feature equivalence type, also -type similarCDS|CDSonly 
    -mark=mark-tag
    -slopexon=$EXONSLOP : mismatch bases per exon
    -sloptutr=$UTRSLOP  : mismatch bases for all of UTR, given typeover=CDS,exon
    -samecds=$SAME_CDS  : minimal CDS overlap to call equal
    -mincds=$MINID_CDS -minutr=$MINID_UTR : minimal overlap to call any equivalence
    -self : use when -over and -in contain same gene IDs
    -nostrand  -strictstrand -summary  -debug 
    
  identity = same location/orientation for gene = (mRNA/exon,CDS)
  -type similar marks equivalence with 'overid/%cdsident.%exonident', 
        with overid/I100 = perfectly equal, and overid/C95.85 = at samecds level
  ** ASSUMES gene model is mRNA > exon > CDS
  ** ASSUMES input gff is ordered by gene records (mRNA/exon/CDS all together per ID)
" unless($optok and $action and $overlaps and @input);


if(defined $dosum) {
  $EQUAL_CUT= $dosum if($dosum =~ /\d/ and $dosum >= 5);
  $dosum=1;
}

$mark  ||= "overg";
my $cdsmark="cds".$mark; #?

## $pctover= $pctover/100.0 if($pctover);  # not used now

$passtypes ||="";
$passtypes =~ s/[,\s]+/\|/g; ## drop for mrnatypes,exontypes
$exontypes =~ s/[,\s]+/\|/g; $mrnatypes =~ s/[,\s]+/\|/g;
$typeover=$exontypes unless($typeover);

if($action =~ /^sum/) { $dosum=1; $action.="null"; }
$actid= ($action =~ /keep/) ? ACT_KEEP : ($action =~ /mark/) ? ACT_MARK 
      : ($action =~ /drop/) ? ACT_DROP : ($action =~ /null/) ? ACT_NULL : ACT_KEEP;
$actid= ACT_MARK_WITH_ID if($actid == ACT_MARK && $action =~ /id/i);

my $onegeneNostrand= $nostrand; # setting per gene, 1-exon things are ambiguous

sub MAIN{}

my $overgffh= undef; 
#upd17: opt to defer collect_over till needed; -SELF? over == input; -sorted per chr

sub open_overlaps {
  my($overlaps)= @_;
  my $overgffh= undef; 
  my $ok= 0;
  if($overlaps =~ /.gz$/){ $ok= open(OVR,"gunzip -c $overlaps |");  $overgffh= *OVR; }
  elsif($overlaps =~ /^(stdin|-)/) { $overgffh= *STDIN; $ok=1; }
  else { $ok= open(OVR,$overlaps); $overgffh= *OVR; }
  die "bad -overlaps=$overlaps" unless($ok);
  return($overgffh,$ok);
}
# .. need to reopen overgffh if SORTEDGFF and @input > 1 ..
#    if($overlaps =~ /.gz$/) { $ok= open(OVR,"gunzip -c $overlaps |");  $overgffh= *OVR; }
# elsif($overlaps =~ /^(stdin|-)/) { $overgffh= *STDIN; $ok=1; }
# else { $ok= open(OVR,$overlaps); $overgffh= *OVR; }
# die "bad -overlaps=$overlaps" unless($ok);

my($ngov, $nxov)=(0,0); 
# unless($SORTEDGFF) { # defer now
#   $overgffh= open_overlaps($overlaps);
#   ($ngov, $nxov)= collect_overlaps($overgffh); 
#   close($overgffh); $overgffh=undef;
# }

foreach my $input (@input) {
  my $inh= *STDIN;
  $ok = ($input =~ /.gz$/) ? open($inh,"gunzip -c $input |") 
        : ($input =~ /^(stdin|-)/) ? $inh= *STDIN
        : open($inh,$input);
  die "bad -input=$input" unless($ok);

  if($SELFSAME) {  drop_overlaps(); $overlaps= $input; }
  unless($SORTEDGFF or @overgenes) {
    ($overgffh,$ok)= open_overlaps($overlaps);
    ($ngov, $nxov)= collect_overlaps($overgffh); 
    close($overgffh); $overgffh=undef;
  }
  
  my $selfover= ($input eq $overlaps)?1:0;
  my($ngin, $nxin, $ngsame)=(0) x 0;
  %sums=(); %sumid=(); 
  
  if($SORTEDGFF) { ($overgffh,$ok)= open_overlaps($overlaps); } # for each input, now read via filter_gff
  
  if($selfover) { 
    ($ngin, $nxin, $ngsame) = self_filter_gff($inh, $overgffh);
  } else {
    ($ngin, $nxin, $ngsame) = filter_gff($inh, $overgffh);
  }
  
  if($SORTEDGFF) { close($overgffh); $overgffh=undef; drop_overlaps(); }
  
  warn"#overlaps over=$overlaps in=$input genes=$ngin same=$ngsame\n" if $debug;
  
  summary($input, $ngov, $nxov, $ngin, $nxin, $ngsame) if ($dosum);
}

#..................


# add summary stats option per
# ngene=10194; genehit=7020; gperfect=912; nexon=43016; exonhit=31932
# sensitivity/specificity ?
sub summary
{
  my($input, $ngov, $nxov, $ngin, $nxin, $ngsame)= @_;
  
  my $ncut= $sums{"C".$EQUAL_CUT}||0; 
  my $xt  = $sums{"X".$EQUAL_CUT}||0; $ncut= $xt if($xt>$ncut);
  my $spec= $ncut/$ngin; # $ngsame/$ngin;  
  my $sens= $ncut/$ngov; # $ngsame/$ngov;  # ?? < change ngsame to opt C/X level
  
  my %sumovc=();
  foreach my $cid (sort keys %sumid) {
    my $noid= scalar(keys %{$sumid{$cid}});
    $sumovc{$cid}= $noid;
    $sens= $noid / $ngov if($cid eq "C".$EQUAL_CUT);
    #? should spec= $noid/$ngin here? NO, counts ids matched not how many ngenein have accurate evidence
  }
  
  # FIXME: need uniq IDs hit of overgenes for sens stat
  # for my $p (10,33,50,66,75,90) { $sumid{"C".$p}{$sd}++ if($c>=$p); }
  
  map{ $_= sprintf "%0.3f",$_; } ($sens,$spec);
  # ngsame == issimilar, any equiv, for type=similar,
     # $sums{Csum}+=$c; $sums{Xsum}+=$x; $sums{NHIT}++; $sums{$ci}++ if($ci);
     # for my $p (33,50,66,80,90) { $sums{"X".$p}++ if($x>=$p); $sums{"C".$p}++ if($c>=$p); }

  my $nperf= ($sums{I}||0) + ($sums{C}||0);
  print "#overgenedup over=$overlaps in=$input equalcut=$EQUAL_CUT\n";
  print "# novergene=$ngov noverexon=$nxov \n";
  print "# ningene=$ngin ninexon=$nxin genehit=$ngsame gperfect=$nperf equalcut=$ncut sens=$sens spec=$spec \n";#exonhit= 
  my $nhit= delete $sums{NHIT} || 1; # same as nghit ?
  my $csum= delete $sums{Csum}||0;  my $cave= int(10*$csum/$nhit)/10;
  my $xsum= delete $sums{Xsum}||0;  my $xave= int(10*$xsum/$nhit)/10;
  my @sk= sort keys %sums;
  print "# identity ave. CDS=$cave Exon=$xave \n";
  print "# identity levels: ", join " ", map{ "$_=$sums{$_}" } sort keys %sums;
  print "\n";
  print "# identity ovgenes: ", join " ", map{ "$_=$sumovc{$_}" } sort keys %sumovc;
  print "\n";
}

sub testgene
{
  my($generecIN, $geneother)= @_;
  
  my @generec= sort _sortgene @$generecIN;
  my $generec= \@generec;

  my $mrna= $generec->[0]; #?? or check type?
  # my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid)= @$mrna;
 
  my ($issame,$samecds, $hasidtag)= samegene($generec) ;
  #? add opt to mark for mostly, partly same : using exon bases overlapped
  # e.g., ${mark}same=ID1,.. or ${mark}66=ID2,66%same) or ${mark}33=ID3, 33%same

  if($samecds) {
    if($actid == ACT_MARK) { $mrna->[8] =~ s/$/;$cdsmark=1/ ;  } 
    elsif($actid == ACT_MARK_WITH_ID) { $mrna->[8] =~ s/$/;$cdsmark=$samecds/ ;  } 
  }
  
  if($issame) {
    if($dosum) {
      my @cuts= (10,33,50,66,75,90);
      push @cuts, $EQUAL_CUT unless( grep{ $EQUAL_CUT == $_ } @cuts );
      
      $sums{NHIT}++;
      unless($hasidtag) { $sums{I}++; }
      else {
      # FIXME: add overID counter for sensitiv calc? must collect all IDs in issame, per CX level
      my @s= split",",$issame; my $i=0;      
      #my $s1=$s[0];
      foreach my $s1 (@s) { $i++;
      my $sd= ($s1 =~ s,^([^/]+)/,,) ? $1 : 0;
      my ($c,$x)= split /\./,$s1; 
      my($ci)= ($c=~s/^([CI])//) ? $1 : "";
      if($c =~ /^\d/) { $x ||= 0;
        if($i==1) {
          $sums{Csum}+=$c; $sums{Xsum}+=$x;  $sums{$ci}++ if($ci);
          for my $p (@cuts) { $sums{"X".$p}++ if($x>=$p); $sums{"C".$p}++ if($c>=$p); }
          }
        $c=$x if($x>$c); # max only here
        for my $p (@cuts) { $sumid{"C".$p}{$sd}++ if($c>=$p); }
        }
      }
      }
    }
    
    return 1 if($actid == ACT_DROP);
    if($actid == ACT_MARK) { $mrna->[8] =~ s/$/;$mark=1/; } 
    elsif($actid == ACT_MARK_WITH_ID) { $mrna->[8] =~ s/$/;$mark=$issame/ ;  } 

    putgene($generec, $geneother) unless($actid == ACT_NULL); #?? gff lines?

    
  } else { 

    return 0 if($actid == ACT_KEEP);     
    putgene($generec, $geneother) unless($actid == ACT_NULL); #?? gff lines?
  }

  return ($issame)?1:0;
}

=item EqualGene table

  EqualGene : run as this tabulator
   $evigene/scripts/overgenedup.pl -in $ingff -over $ovgff $ovflags \
   -slopexon=8 -type similarCDS -mincds=10 -minutr=33 -act markid -mark overg | \
   grep mRNA | perl -ne '($g,$ov)= m/(?:ID|overg)=([^;\s]+)/g;  $ov||="na"; \
   $og= (m/oid=([^;\s]+)/) ?$1:"noid"; \
   ($r,$b,$e)=(split"\t")[0,3,4]; $r=~s/^Scaffold(\S+)/${1}sc/i; \
   print join("\t",$g,$og,$ov,"$r:$b-$e"),"\n";' \
   > equal/$inam.over$onam.tab5

=item gmap2gff mRNA quals

chr4	evg5atcds	mRNA	9805652	9814910	100	-	.	ID=Arath5EVm000112t1;Target=Arath5EVm000112t1 1 5466;
aalen=1821;cov=100.0;match=5466;nexon=35;pid=100.0;qlen=5466;
aalen=1821,94%,complete;offs=127-5592;oid=Arath3EVm000114t1

mapqual:  85a,100i,2469l,4x,Spl:29%,KB668900 == %align,%ident,bases,exons,Spl:split

=cut

use constant EXTEND17EQTAB =>  1;

sub evGeneAltID {
  my($id)=@_;
  $id=~s/_[CG]\d+$//; # split,dup tags
  my($gd,$ti)= ($id=~m/(.+)(\d+)$/) ? ($1,$2) : ($id,1);
  return($gd,$ti);
}

## this is slow...
# my @overlapids=();
# sub altsNotOverOld { # only for equalgene -in self.gff -over self.gff
#   my($gid,$ovg)= @_;
#   my($gd,$ti)=evGeneAltID($gid);
#   my @ovg= map{ s,/.*,,; s/_[CG]\d+$//; $_ } split",",$ovg;
#   my %haveid= map { $_,1 } ($gid,@ovg);
#   unless(@overlapids) { @overlapids= sort keys %overlapids; }
#   my @otheralts= grep{ /^$gd/ and not $haveid{$_} } @overlapids;
#   return @otheralts;
# }

#? faster than altsNotOver[A] ? yes
my %overlapgenealt=();
sub altsNotOverB { # only for equalgene -in self.gff -over self.gff
  my($gid,$ovg)= @_;
  my($gd,$ti)=evGeneAltID($gid); $gid=$gd.$ti; # remove split,dup tags
  my @ovg= map{ s,/.*,,; s/_[CG]\d+$//; $_ } split",",$ovg;
  my %haveid= map { $_,1 } ($gid,@ovg);
  unless(%overlapgenealt) { 
    my @overlapids= sort keys %overlapids; 
    for my $otd (@overlapids) {
      $otd=~s/_[CG]\d+$//; # yuck .. remove these or not ? split,dup tags
      my($og,$ot)= evGeneAltID($otd);
      push @{$overlapgenealt{$og}}, $otd;
    }
  }
  my @otheralts=();
  if(my $allalts= $overlapgenealt{$gd}) { # should always have self gid
    @otheralts= grep{ not $haveid{$_} } @$allalts;
  }
  return @otheralts;
}

#  mrnaMapQual($tattr, $nx, $tp); # tp == genescore default, maybe invalid
sub mrnaMapQual {
  my($attr,$nxin,$gscorein)=@_;
  # missing vals?  dup map? _G2;path=2/2
  $nxin||=1;
  my $gs= $gscorein||0; # add gff/annot genescore if available
  my($cov)  = ($attr=~m/cov=([\d\.]+)/)?$1:0; $cov=int(0.5+$cov);
  my($pid)  = ($attr=~m/pid=([\d\.]+)/)?$1:0; $pid=int(0.5+$pid);
  my($qlen) = ($attr=~m/qlen=([\d\.]+)/)?$1:0;
  my($nexon)= ($attr=~m/nexon=([\d\.]+)/)?$1:$nxin;
  if($attr=~m/scoresum=([\d\.-]+)/) { $gs=$1; }
  # ID=Arath5EVm007498t1_C1; path=1/2;chimera=breakpoint at 796;chim2=chr1:1959033-1959652:.
  # ID=Arath5EVm007498t1_C2; path=2/2;chimera=breakpoint at 796;chim1=chr4:10609405-10610200:.
  my $sp="";
  # my($npa)= ($attr=~m/path=([^;\s]+)/) ? $1 : 0;
  my($chimat)= ($attr=~m/chim[12]=([^;\s]+)/) ? $1 : 0;
  if($chimat) {
    my $scov= 100 - $cov; $scov=1 if($scov<1); # not right?
    $sp="$scov%,$chimat";
  }
  my $mapq= join",", $cov.'a', $pid.'i', $qlen.'l', $nexon.'x';
  $mapq.=",Spl:$sp" if($sp); # Split tag
  $mapq.=",".$gs."gs" if($gs);
  return($mapq);
}


sub putEqualTab
{
  my ($generec, $geneother)= @_;
  my $mrna= $generec->[0]; # or check type? 
  # my ($mrna)= grep { $_->[2] eq "mRNA" } @$generec;
  my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid)= @$mrna;  
  
  my $ov= ($tattr =~ m/;$mark=([^;\n]+)/) ? $1 : "na";
  my $og= ($tattr =~ m/;oid=([^;\s]+)/) ? $1 :"noid";
  (my $rk=$ref) =~ s/^scaffold(\S+)/${1}sc/i; $rk =~ s/^contig(\S+)/${1}ct/i;
  if(EXTEND17EQTAB) {
    #u 1703: extend equaltab with mapqual (basic map.attr) and altnotover cols
    #u eqgene mapqual: 85a,100i,2469l,4x,Spl:29%,KB668900 == %align,%ident,bases,exons,Spl:split

    my $altnotover="na";
    if($addAltsNotOver) {
      my @altnotover= altsNotOverB($gid,$ov);
      $altnotover= join",", @altnotover if(@altnotover);
    }
    my @xn= grep { $_->[2] eq "exon" } @$generec;
    unless(@xn) { @xn= grep { $_->[2] eq "CDS" } @$generec; }
    my $nx= @xn;
    my $mapqual= mrnaMapQual($tattr, $nx, $tp); # tp == genescore default, maybe invalid
    print join("\t",$gid,$og,$ov,"$rk:$tb-$te:$to",$mapqual,$altnotover),"\n";
  } else {
    print join("\t",$gid,$og,$ov,"$rk:$tb-$te:$to"),"\n";
  }
}

sub putgene
{
  return putEqualTab(@_) if($EqualGene);
  my ($generec, $geneother)= @_;
  ##  my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid]; 
  foreach my $ft (@$generec, @$geneother) { 
    if(ref $ft) { my @v= @$ft; print join("\t",@v[0..8]),"\n" if(@v>4); }
    }
}


sub _sortgene  
{
  # ($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid)
  # my($ta,$tb)= map{  (/gene/)?1:(/mRNA/)?2:(/CDS|exon/)?3:4; } ($a->[2],$b->[2]); #? is CDS==exon sort bad?
  my($ta,$tb)= map{ (/gene/)?1:(/mRNA/)?2:(/exon/)?3:(/CDS/)?4:5; } ($a->[2],$b->[2]); 
  return ($ta <=> $tb)
      || ($a->[0] cmp $b->[0]) # ref : should all be same for gene record
      || ($a->[3] <=> $b->[3]) # begin
      || ($a->[4] <=> $b->[4]) # end
      ;
}

# fixme for "." strand compare

sub _eqstrand { 
   if($orientstrict){ return($_[0] eq $_[1])?1:0; }
   elsif($onegeneNostrand) { return 1; }
   else { return($_[0] eq $_[1] || $_[0] eq '.' || $_[1] eq '.' )?1:0; }
}

sub _samefeat 
{
  my($a,$b)= @_;
  return 0 unless( ref($a) and ref($b));
# ($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid)

  my($ab,$ae)=($a->[3],$a->[4]);
  my($bb,$be)=($b->[3],$b->[4]);

if($MRNA_SPAN_WRONG) { # ncbi gff mRNA == gene span, not mRNA, for alt-tr; fix source
  ($bb,$be)= ($ab,$ae)
    if($a->[2] eq "mRNA" and $b->[2] eq "mRNA" and
      ((abs($ab - $bb) <= $EXONSLOP) or (abs($ae - $be) <= $EXONSLOP)) );
}
  if($EXONSLOP) {
  return ($a->[0] eq $b->[0]) # ref
      && ($a->[2] eq $b->[2]) # typ
      && (abs($ab - $bb) <= $EXONSLOP) # begin
      && (abs($ae - $be) <= $EXONSLOP) # end
      && _eqstrand($a->[6], $b->[6]) # ($a->[6] eq $b->[6]) # orient; allow "." to match + or - ?? ** FIXME yes
      ;
  } else {
  return ($a->[0] eq $b->[0]) # ref
      && ($a->[2] eq $b->[2]) # typ
      && ($ab == $bb) # begin
      && ($ae == $be) # end
      &&  _eqstrand($a->[6], $b->[6]) # ($a->[6] eq $b->[6]) # orient; allow "." to match + or - ??
      ;
  }
}

sub _samebases 
{
  my($a,$b)= @_;
  return (0,0) unless( ref($a) and ref($b));
# ($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid)
 
  my $wmax= 1 + _max( $b->[4] - $b->[3], $a->[4] - $a->[3]);
  if ( ($a->[0] eq $b->[0]) # ref
      && ($a->[2] eq $b->[2]) # typ
      &&  _eqstrand($a->[6], $b->[6]) # ($a->[6] eq $b->[6]) # orient; allow "." to match + or - ??
      && ($a->[3] < $b->[4]) # begin - end
      && ($a->[4] > $b->[3]) # end - beg
      )
    {
      my $omin= _min( $b->[4], $a->[4]) - _max( $b->[3], $a->[3]);
      # my $wmax= 1 + _max( $b->[4] - $b->[3], $a->[4] - $a->[3]);
      return( $omin, $wmax); #?
    }
  return (0,$wmax);
}


sub _samegene  
{
  my($agene,$bgene)= @_;
  my $na= @$agene;
  my $nb= @$bgene;
  return 0 unless($na == $nb); ## fixme: allow some end-exon slop
  for(my $i=0; $i<$na; $i++) {
    my $issame= _samefeat( $agene->[$i], $bgene->[$i]);
    return 0 unless $issame;
  }
  return 1;
}

sub _samegenecds  # see instead _similargene now better
{
  my($agene,$bgene)= @_;  # agene=input; bgene=over
  my $na= @$agene;
  my $nb= @$bgene;
# ($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid)
  
  # option: classify _sameCDS, do both?
  my $gdiff= 0;
  my $nn= _max($na,$nb);
  my(@acds,@bcds);
  for(my $i=0; $i<$nn; $i++) {
    my $issame= ($i>=$na or $i>=$nb)? 0 : _samefeat( $agene->[$i], $bgene->[$i]);
    $gdiff++ unless($issame); # return 0 unless $issame;
    
    if($i < $na and $agene->[$i]->[2] eq "CDS") { push(@acds,$i); }
    if($i < $nb and $bgene->[$i]->[2] eq "CDS") { push(@bcds,$i); }
  }

  return(1,1,1,"I100") if($gdiff==0);
  
  my($nca,$ncb)= (scalar(@acds),scalar(@bcds));
  #??# return(0,0) if($nca != $ncb); #? allow CDS num mismatch ?
  $onegeneNostrand= ($NostrandOneExon and ($nca==1 or $ncb==1)) ? 1 : $nostrand;

  $nn= _max($nca,$ncb);
  my($cdiff, $csame, $cba,$cbb, $cea,$ceb, $acdsw, $bcdsw)= (0) x 10;
  my($lasti,$lastj)= (-1,-1);

if(1) {
# FIXME: bug here, similar* : need to look inside acds,bcds for matching subset, not start from 0
  for ( my ($i,$j)=(0,0); $i < $nca or $j < $ncb; ) {
    my $ia= $acds[$i];  
    my $ib= $bcds[$j];
    my ($aok,$bok)= (defined $ia, defined $ib); # ? ($i<$nca, $b<$ncb)
    # my $issame= ($aok and $bok) ? _samefeat( $agene->[$ia], $bgene->[$ib]) : 0;
    my ($sameb)= ($aok and $bok) ? _samebases( $agene->[$ia], $bgene->[$ib]) : 0;
    $csame += $sameb;

    if($aok) { 
      $acdsw += 1 + $agene->[$ia]->[4] - $agene->[$ia]->[3] if($i > $lasti);
      $cba= $agene->[$ia]->[3] if($i==0) ;
      $cea= $agene->[$ia]->[4] if($i==$nca-1);
      }
    if($bok) { 
      $bcdsw += 1 + $bgene->[$ib]->[4] - $bgene->[$ib]->[3] if($j > $lastj); 
      $cbb= $bgene->[$ib]->[3] if($j==0) ;
      $ceb= $bgene->[$ib]->[4] if($j==$ncb-1);
      }
    ($lasti,$lastj)= ($i,$j);
    
    unless($sameb) { 
      $cdiff++;
      if(($aok and $bok) and $agene->[$ia]->[4] < $bgene->[$ib]->[3]) { $i++; }
      elsif(($aok and $bok) and $bgene->[$ib]->[4] < $agene->[$ia]->[3]) { $j++; }
      else { $i++; $j++; } #??
    } else { 
      $i++; $j++; 
    }
    
  }
} 

  my $cdsmax= _max( $acdsw, $bcdsw);
  my $cdsmin= _min( $acdsw, $bcdsw);
  my $cdsident= ($cdsmax>0) ? int( 0.5 + 100 * $csame / $cdsmax) : 0;
  # add: issimilar= cds_aencloses if $csame == $bcdsw;  cds_ainside if $csame == $acdsw
  
  my $utrdiff="";
  my($autrw, $butrw, $utrwdiff)= (0,0,0);
  if( 1 or $cdiff==0) { 
    # option: measure utr diff: utr5, utr3 n and size
    # this allows utr to be at any location .. given cds are same, should be ok
  my ($ab,%utr);  
  $utr{w5}{a}= $utr{w5}{b}= $utr{w3}{a}= $utr{w3}{b}= 0;
  
  for(my $i=0, $ab='a'; $i<$na; $i++) {
    next if($agene->[$i]->[2] ne "exon");
    my $xb= $agene->[$i]->[3];
    my $xe= $agene->[$i]->[4];
    if($xb < $cba) { $utr{w5}{$ab} += 1 + _min($xe,$cba) - $xb; $utr{n5}{$ab}++; }
    elsif($xe > $cea) { $utr{w3}{$ab} += 1 + $xe - _max($xb,$cea) ; $utr{n3}{$ab}++; }
    }  
  for(my $i=0, $ab='b'; $i<$nb; $i++) {
    next if($bgene->[$i]->[2] ne "exon");
    my $xb= $bgene->[$i]->[3];
    my $xe= $bgene->[$i]->[4];
    if($xb < $cbb) { $utr{w5}{$ab} += 1 + _min($xe,$cbb) - $xb; $utr{n5}{$ab}++; }
    elsif($xe > $ceb) { $utr{w3}{$ab} += 1 + $xe - _max($xb,$ceb) ; $utr{n3}{$ab}++; }
    }
  
  $autrw= $utr{w5}{a} + $utr{w3}{a};
  $butrw= $utr{w5}{b} + $utr{w3}{b};
  
    # ** many of these differ by only a few bases of utr; call same if wdiff < UTRSLOP
  my ($maxw, $maxn)= (0,0);
  foreach my $k (sort keys %utr) { 
    my($ua, $ub)= ($utr{$k}{a}||0, $utr{$k}{b}||0);
    my $d= $ua - $ub;  # diff is input - overlap, +=excess, -=missing
    if($d != 0) {
      $utrdiff.="$k=$d,"; #? need this
      if($k =~ /w/) { $maxw += abs($d); } elsif($k =~ /n/) { $maxn += abs($d); }
      }
    }
  $utrwdiff= $maxw;
  $utrdiff = "$maxw,$utrdiff"; # maxw == 0 ok
  $utrdiff =~ s/,$//;
  
    ## not this: and $maxn <= $UTRSLOPN ; use UTRSLOP only
  if($cdiff==0 and $maxw <= $UTRSLOP) { $gdiff= 0; } ## ( and $utrloc ~ sameloc); ??
  #? elsif($acdsw < $TINYCDS or $bcdsw < $TINYCDS) { $cdiff=1; } # dont call same
  }
  
  my $umax  = _max( $autrw, $butrw);
  my $usame = ($umax - $utrwdiff); $usame=0 if($usame<0);
  my $trident = ($umax>0 or $cdsmax>0) ? int( 0.5 + 100 * ($csame + $usame) / ($cdsmax + $umax) ) : 0;

  my $genesize= int( $cdsmax + $umax);
  my $genesame= ($cdsident>0 or $trident>0) ? "$cdsident.$trident" : 0;
  my $issimilar= ($cdsident > $MINID_CDS && $cdsident > 0.01)?1:0; # _max(minid,0.01)
  
  if($MINID_CDS>0 and $cdsident < $MINID_CDS) { $genesame=0; }

  my $flag="";
  if($cdsident >= $SAME_CDS) { $flag = "C"; }  
  return ($gdiff==0, $cdiff==0, $issimilar, $flag.$genesame, $genesize);  
  # return ($gdiff==0, $cdiff==0, $utrdiff);
}


sub checkgene { # debug sub
  my($agene,$info)= @_;
  my $na= @$agene;
# ($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid)
  return 1 unless($debug);
  
  my(@acds, @aexon, @amrna, @aother);
  my $err=0;
  
  for(my $i=0; $i<$na; $i++) {
    ##my $issame= ($i>=$na or $i>=$nb)? 0 : _samefeat( $agene->[$i], $bgene->[$i]);
    $err++ if($agene->[$i]->[3] < 1 or $agene->[$i]->[4] < 1);
    if($i < $na and $agene->[$i]->[2] eq "CDS") { push(@acds,$i); }
    elsif($i < $na and $agene->[$i]->[2] eq "exon") { push(@aexon,$i); }
    elsif($i < $na and $agene->[$i]->[2] eq "mRNA") { push(@amrna,$i); }
    elsif($i < $na) { push(@aother,$i); }
  }

    ## have some ncbi no-cds genes
  if($err==0 and @amrna == 1 and @aexon > 0 ) {
    return 1;
  } else { #error
    warn "\n#ERR.checkgene $info nerr=$err, nmrna=",scalar(@amrna)," nexon=",scalar(@aexon),"gn=",$agene->[0]->[8],"\n";    
    print "\n";
    print "\n#ERR.checkgene $info nerr=$err, nmrna=",scalar(@amrna)," nexon=",scalar(@aexon),"gn=",$agene->[0]->[8],"\n";    
    print "#ERR.checkgene=",join(", ",@{$agene->[0]}),"\n";
    print "\n";
    return 0;
  }
}



sub _similargene  
{
  my($agene,$bgene)= @_;  # agene=input; bgene=over
  my $na= @$agene;
  my $nb= @$bgene;
# ($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid)
  
  # option: classify _sameCDS, do both?
  my $gdiff= 0;
  my $nn= _max($na,$nb);
  my(@acds,@bcds, @aexon, @bexon);
  
  for(my $i=0; $i<$nn; $i++) {
    my $issame= ($i>=$na or $i>=$nb)? 0 : _samefeat( $agene->[$i], $bgene->[$i]);
    $gdiff++ unless($issame);  
    
    if($i < $na and $agene->[$i]->[2] eq "CDS") { push(@acds,$i); }
    if($i < $nb and $bgene->[$i]->[2] eq "CDS") { push(@bcds,$i); }
    if($i < $na and $agene->[$i]->[2] eq "exon") { push(@aexon,$i); }
    if($i < $nb and $bgene->[$i]->[2] eq "exon") { push(@bexon,$i); }
  }
  
  return(1,1,1,"I100") if($gdiff==0);
  
  my($nca,$ncb)= (scalar(@acds),scalar(@bcds));
  my($nxa,$nxb)= (scalar(@aexon),scalar(@bexon));
  # OPTION: ignore strand when nxa==1 or nxb==1 
  $onegeneNostrand= ($NostrandOneExon and ($nxa==1 or $nxb==1)) ? 1 : $nostrand;
  
  $nn= _max($nca,$ncb);
  my($cdiff, $csame, $cba,$cbb, $cea,$ceb, $acdsw, $bcdsw,$lasti,$lastj)= (0) x 10;

# FIXME2: Probably bug, when overlap exon sizes differ greatly,
#         Aexon1tiny - Bexon1long = 1% over ; Aexon2long - Bexon2tiny = 1%over, final score is 1-2%over
#         but total base overlap is 99%, from exon split disagreement (e.g. in TE genes, but others too)
#          A:    --  ----------------------------------
#          B:   ---------------------------   ---     -      xsame/xmax is tiny here from exon shift
#
# >> redo exon overlap as max align transcript base overlap regardless of exon breaks
#    but keep cds overlap method same? diff introns mean diff protein

# FIXME: bug here, similar* : need to look inside acds,bcds for matching subset, not start from 0

  ($lasti,$lastj)= (-1,-1);
  for ( my ($i,$j)=(0,0); $i < $nca or $j < $ncb; ) {  # USE i or j; need full acdsw,bcdsw
    my $ia= $acds[$i];  
    my $ib= $bcds[$j];
    #  my ($aok,$bok)= (defined $ia, defined $ib); # ? ($i<$nca, $b<$ncb)
    my $aok= ($i < $nca)?1:0; my $bok= ($j < $ncb)?1:0;
    my $abok= ($aok and $bok)?1:0; ## this may have been bug: need 1:0 ??
    my ($sameb, $maxb)= ($abok) ? _samebases( $agene->[$ia], $bgene->[$ib]) : (0,0);
    $csame += $sameb;


=item upd17: cds-overlap bugfix
  cure for obscure CDS align bug, e.g. CDS overlap but offset starts: case ($sameb>0) below
  new: Daplx6cgEVm001230t1	  Daplx5cEVm003968t3/54.63,
  old: Daplx6cgEVm001230t1    Daplx5cEVm003968t3/3.63,
for these exons
  scaffold_352	dpxevg5a	  CDS	45291	46040	1	-	.	Parent=Daplx6cgEVm001230t1    c1
  scaffold_352	dpxevg5a	  CDS	46420	47433	1	-	.	Parent=Daplx6cgEVm001230t1  | c2
  scaffold_352	evg3daplx5	CDS	47411	49374	1	-	.	Parent=Daplx5cEVm003968t3   | d1 .. overlap offset
  scaffold_352	dpxevg5a	  CDS	47485	49374	1	-	.	Parent=Daplx6cgEVm001230t1  | c3 .. here
  scaffold_352	dpxevg5a	  CDS	49459	49567	1	-	.	Parent=Daplx6cgEVm001230t1
  scaffold_352	evg3daplx5	CDS	49459	49565	1	-	.	Parent=Daplx5cEVm003968t3

=cut
    
if(1) { #upd17: revise as per exons 
    my($ab,$ae,$bb,$be)= (0) x 10;
    if($aok) { 
      ($ab,$ae)= ($agene->[$ia]->[3], $agene->[$ia]->[4]);
      $acdsw += 1 + $ae - $ab if($i > $lasti); 
      $cba= $ab if($i==0) ;
      $cea= $ae if($i==$nca-1);
      }
    if($bok) { 
      ($bb,$be)= ($bgene->[$ib]->[3], $bgene->[$ib]->[4]);
      $bcdsw += 1 + $be - $bb if($j > $lastj); 
      $cbb= $bb if($j==0) ;
      $ceb= $be if($j==$ncb-1);
      }
    ($lasti,$lastj)= ($i,$j);

    if($maxb > 0 and $sameb >= $maxb) {
      $i++; $j++;
      
    } elsif($sameb>0) { # and $sameb < $maxb .. look for more align to maxb
      if($abok and $ae < $be) { $i++; } # see if next aexon hits this bexon
      elsif($abok and $be < $ae) { $j++; } # other way      
      else { $i++; $j++; }

    } else { # if($sameb < 1)
      $cdiff++;
      if($abok and $ae < $bb) { $i++; }
      elsif($abok and $be < $ab) { $j++; }
      else { $i++; $j++; }
    }
    
} else {    #old
    if($aok) { 
      $acdsw += 1 + $agene->[$ia]->[4] - $agene->[$ia]->[3] if($i > $lasti);
      $cba= $agene->[$ia]->[3] if($i==0) ;
      $cea= $agene->[$ia]->[4] if($i==$nca-1);
      }
    if($bok) { 
      $bcdsw += 1 + $bgene->[$ib]->[4] - $bgene->[$ib]->[3] if($j > $lastj); 
      $cbb= $bgene->[$ib]->[3] if($j==0) ;
      $ceb= $bgene->[$ib]->[4] if($j==$ncb-1);
      }
    ($lasti,$lastj)= ($i,$j);
    
    unless($sameb) { 
      $cdiff++;
      if($abok and $agene->[$ia]->[4] < $bgene->[$ib]->[3]) { $i++; }
      elsif($abok and $bgene->[$ib]->[4] < $agene->[$ia]->[3]) { $j++; }
      else { $i++; $j++; } #??
    } else { 
      $i++; $j++; 
    }
} # old upd17
    
  }

  my $cdsmax= _max( $acdsw, $bcdsw);
  my $cdsmin= _min( $acdsw, $bcdsw);
  $cdsmax= $acdsw if($THISSIZE); # upd16.02
  my $cdsident= ($cdsmax>0) ? int( 0.5 + 100 * $csame / $cdsmax) : 0;
  # add: issimilar= cds_aencloses if $csame == $bcdsw;  cds_ainside if $csame == $acdsw

  $nn= _max($nxa,$nxb);
  my($xdiff, $xsame, $axw, $bxw)= (0) x 10;
  
# FIXME: find align subset, same as cds
# FIXME2: >> redo exon overlap as max align transcript base overlap regardless of exon breaks; not same as cds

# FIXME ODD: 30jul wacky new results from same data, esp exon = 0 when CDS = 99;
  # die unless(checkgene($bgene, "_similargene.ov"));
  # ok here, still get wacky  C100.00  = no exon over when all CDS exon over
  
if(1) {

  ($lasti,$lastj)= (-1,-1);
  for ( my ($i,$j)=(0,0); $i < $nxa or $j < $nxb; ) { # USE i< OR j <
    my $ia= $aexon[$i];  
    my $ib= $bexon[$j];
    #?? my ($aok,$bok)= (defined $ia, defined $ib); # ? ($i<$nca, $b<$ncb)
    my $aok= ($i < $nxa)?1:0;
    my $bok= ($j < $nxb)?1:0;
    my $abok= ($aok and $bok)?1:0; ## this may have been bug: need 1:0 ??
    
    my ($sameb, $maxb)= ($abok) ? _samebases( $agene->[$ia], $bgene->[$ib]) : (0,0);
    $xsame += $sameb;
    
    my($ab,$ae,$bb,$be)= (0) x 10;
    if($aok) { 
      ($ab,$ae)= ($agene->[$ia]->[3], $agene->[$ia]->[4]);
      $axw += 1 + $ae - $ab if($i > $lasti); # WHOO, only if $i is NEW
      }
    if($bok) { 
      ($bb,$be)= ($bgene->[$ib]->[3], $bgene->[$ib]->[4]);
      $bxw += 1 + $be - $bb if($j > $lastj); 
      }
    ($lasti,$lastj)= ($i,$j);
            
    # new: 
    if($maxb > 0 and $sameb >= $maxb) {
      $i++; $j++;
      
    } elsif($sameb>0 ) { # and $sameb < $maxb .. look for more align to maxb
      if($abok and $ae < $be) { $i++; } # see if next aexon hits this bexon
      elsif($abok and $be < $ae) { $j++; } # other way      
      else { $i++; $j++; }

    } else { # if($sameb < 1)
      $xdiff++;
      if($abok and $ae < $bb) { $i++; }
      elsif($abok and $be < $ab) { $j++; }
      else { $i++; $j++; }
    }
      
  }
} 


  my $xmax  = _max( $axw, $bxw);
  $xmax= $axw if($THISSIZE); # upd16.02
  my $trident = ($xmax>0)   ? int( 0.5 + 100 * $xsame / $xmax) : 0;

  my $usame = _max( 0, $xsame - $csame);
  my $umax  = _max( $usame, _max( $axw - $acdsw, $bxw - $bcdsw));
  
  my $genesize= int( $cdsmax + $umax); # not same as xmax
  my $genesame= ($genesize>0) ? int( 0.5 + 100 * ($usame + $csame) / $genesize) : 0;
  my $issimilar= $genesame > 0 ? 1 : 0;

## this gets wacky results: C0 = cds same > 95, but xsame == 0
#  my $genesize= int( 0.5 * ($axw + $bxw));
#  my $genesame= ($xsame  < 10) ? 0 : int( 0.5 +  100 * ( $xsame ) / $genesize ); #? max or min or ave ?
  
  # * Option: split genesame to cdsident, utrident, percent of overlaps
  # ? Option: min-ident (cds,utr) to call similar: here or in caller?
  # want to call same-locus using only cdsident ?
  # return genesame = "cdsident.utrident" as pct max

if($SHOW_CDSUTR) {
  # FIXME: need to format pad .trident w/ 0 for numeric compare
  if($trident < 0) {  $trident = "00"; } # buggy data
  elsif($trident < 10) { $trident = "0$trident"; } elsif($trident >= 100) { $trident="99"; } # do we ever get 100.100 ? should be all I100
  $genesame= ($cdsident>0 or $trident>0) ? "$cdsident.$trident" : 0;
}


  my $flag="";
  if($cdsident >= $SAME_CDS) { $flag = "C"; } ## $csame/$cdsmax ;   $cdiff == 0 and 
  elsif($MINID_CDS>0 and $cdsident < $MINID_CDS and ($MINID_UTR == 0 or $trident < $MINID_UTR))
  {
  # * use trident here not utrident, where exon overlap can be noncds x cds
  # check for large-cds over most of small-cds gene, but pct ID to cdsmax is small:
  $issimilar=0 if($cdsmin < 1 or ($MINID_CDS > 100 * $csame / $cdsmin) );
  }

  return ($gdiff==0, $cdiff==0, $issimilar, $flag.$genesame, $genesize); ## $utrdiff
}




sub _sortbestid {
  # best == smallest base difference >= 0
  my($aa,$av)= split"/",$a; $av=~s/,.*// if($av); $av ||= 0; 
  my($bb,$bv)= split"/",$b; $bv=~s/,.*// if($bv); $bv ||= 0;
  return ($av <=> $bv) || $aa cmp $bb; 
}

## Cnum _sortbestbases bug:
# equal/ncbigene2.overacypi1.tab3:XM_001951885.2          ACYPI000496-RB/78.55,ACYPI000496-RA/C99.69
# equal/ncbigene2.overacypi1.tab3:XM_001951896.2          ACYPI000496-RA/78.65,ACYPI000496-RB/C99.64
#    C99.69 should sort > 78.55 : s/\D+// is bad for 12.34 
# other sort bug: GNOMON:1740193/C100.75,NM_001126134.2/I100   I100 > C100.75
## and  50.11 > 50.2; need to format .digits with 0 pad

sub _sortbestbases {
  # best == largest bases same >= 0
  my($ac,$bc)=(0,0);
  my($aa,$av)= split"/",$a; if($av){ $av=~s/,.*//; $ac=$1 if($av=~s/^([^\d\.]+)//);  $av+=10 if($ac eq "I"); } $av ||= 0;  # was s/\D+//
  my($bb,$bv)= split"/",$b; if($bv){ $bv=~s/,.*//; $bc=$1 if($bv=~s/^([^\d\.]+)//);  $bv+=10 if($bc eq "I"); } $bv ||= 0;
  return ($bv <=> $av) || $aa cmp $bb; 
}

sub samegene
{
  my($generec)= @_;   
  my (@lid, @cid, %didr);

  # die unless(checkgene($generec, "samegene.in")); # ok here
  my $mrna= $generec->[0]; #?? or check type?
  my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid)= @$mrna;
  
  return 0 unless($overlaplist->{$ref});
  my $cdstest= ($typeover =~ /CDS/i) ? 1 : 0;
  my $cdsonly= ($typeover =~ /CDSonly/i) ? 1 : 0;
  my $similartest= ($typeover =~ /similar/) ? 1 : 0;
  $SHOW_CDSUTR=1 if($similartest and $cdstest); #?? use other option? always for similar??
  
  my $hasidenttag=  $cdstest||$similartest;
  
  my @bins= (int($tb/$BINSIZE) .. int($te/$BINSIZE));
  foreach my $ib (@bins) {
    $overlaplist->{$ref}{$ib} or next;
    my @rgenes= @{$overlaplist->{$ref}{$ib}}; # 30jul11: change to IDlist, @overgenes index
    foreach my $rgene (@rgenes) {
      next if($didr{$rgene}++);
      my $ovgene= $overgenes[$rgene];
  # die unless(checkgene($ovgene, "samegene.ov=$rgene")); # all ok here
      my($rref,$rsrc,$rtyp,$rb,$re,$rp,$ro,$rph,$rattr,$lid)= @{$ovgene->[0]}; #was @{$rgene->[0]};
      next unless($rb < $te and $re > $tb); # no gene overlap
      
      # FIXME: WRONG SOMETIMES SELF-ID match desired for bestgenes including over test set
      next if(!$SELFKEEP and $lid eq $gid); # truly same gene.. 
      
      my($issame, $issimilar, $cdssame, $cxident, $genesize)= (0) x 10;
      $onegeneNostrand= $nostrand; # setting per gene, 1-exon things are ambiguous
      
      if($similartest) {
        ## FIXME: this misses real similargenes
        ($issame, $cdssame, $issimilar, $cxident, $genesize)= _similargene($generec, $ovgene); ## $rgene
        
      } elsif($cdstest) { 
          # now returns same as _similargene() ; cxident was utrdiff == [IC]cdsident.trident
        ($issame, $cdssame, $issimilar, $cxident, $genesize) = _samegenecds($generec, $ovgene); 
        $issame= $cdssame if($cdsonly and not $issame);
        
      } else { 
        $issame = _samegene($generec, $ovgene); $cdssame=$issame; 
      }
      
      # my $lid= $rgene->[0]->[-1];
      $lid .="/$cxident" if($cxident);
      push @lid, $lid if($issame or $issimilar);
      push @cid, $lid if($cdstest && $cdssame); # drop this?  
      }
    }
    
  my($lid,$cid)=(0,0);
  if(@lid) { 
    my %lid= map{$_,1}@lid; 
    if($similartest or $cdstest) {    
    @lid= sort _sortbestbases keys %lid; # ? only 1st by most same bases, or need many matches, i.e. 1 input == 2 over genes
    
    $lid= join",", @lid;
    } else {
    $lid= join",", sort _sortbestid keys %lid;     #? drop off 2ndary matches if score not same??
    }
    
  } elsif($cdsonly and @cid) { # only if not @lid ?
    my %cid= map{$_,1}@cid; 
    # $cid= join",", sort _sortbestid keys %cid; 
    @cid= sort _sortbestbases keys %cid; # ? only 1st by most same bases, or need many matches, i.e. 1 input == 2 over genes
    $cid= join",", @cid;  
  }

    #? return both samegene and samecds id lists ?
  return ($lid, $cid, $hasidenttag);
}


sub filter_gff {
  my( $ingffh, $overgffh )= @_;
  
  my ($atchr, $lastchr, $ng,$nx,$nr,$nsame,$nhit,$nwarn)= (0) x 10;
  my $printpass=1;
  $printpass=0 if($EqualGene or $actid == ACT_DROP or $actid == ACT_KEEP or $actid == ACT_NULL);
  my $nocomm= ($EqualGene or $actid == ACT_NULL) ? 1 : 0;
  my @generec=();
  my @geneother=();
  my $gid=undef;
  
  while(<$ingffh>) {
    unless(/^\w/){ print unless($nocomm or /^(#n |$)/); next; }
    my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr)= split"\t";
    $nr++; chomp($tattr);

    ## should collect and put after gene/exons; put all non-exons by default?
    # if($passtypes and "$typ.$src" !~ m/$passtypes/) { print if $printpass; next; } 
    next if($SKIPREF and $ref =~ m/^($SKIPREF)/);

    if($SORTEDGFF and $ref ne $lastchr) {
      drop_overlaps(); # last ref
      ($atchr, $ng, $nx, $nwarn)= collect_overlaps_bychr($overgffh, $ref); 
    }
    $lastchr= $ref;
    
    my($pid,$id)=(0,0); 
    if($tattr =~ m/\bID=([^;]+)/) { $id=$1; $gid=$id if($typ =~ /^($mrnatypes)$/); }
    if($tattr =~ m/\bParent=([^;]+)/) {  $pid=$1; 
      $gid=$pid unless($gid and ($typ =~ /^($mrnatypes)$/));  }
    unless(defined $gid) { $gid = "N".$ng; }

    ## need to collec FULL gene record to test gene overlap...
    #?? for -SELF? overgffh == ingffh could reuse this gff record
    my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid]; 

    if($typ =~ /^($mrnatypes)$/) {  
      # allow gene and mRNA types .. and/or keep other types in generec
      $nsame += testgene(\@generec, \@geneother) if(@generec);
      
      $ng++;
      @generec = ($rloc); # maybe best store as array of [gffcols], sorted by type-loc
      @geneother= ();
      
    } elsif($typ =~ /^($exontypes)$/) {
      my $rid= $generec[0]->[9] || 0;
      if(not $IgnoreOutOfOrder and $rid ne $pid) { warn "#ERR: out of order gene record: mRNA=$rid, $typ=$pid\n" if($nwarn++ < 9); } # limit warns
      else { push @generec, $rloc; $nx++; }
      
    } elsif( ! $passtypes or "$typ.$src" =~ m/$passtypes/ ) {
      push @geneother, $rloc;
    }
  }
  
  $nsame += testgene(\@generec, \@geneother) if(@generec);
  return ($ng,$nx,$nsame,$nhit);
}

sub self_filter_gff { # $ingffh == $overgffh
  my( $ingffh, $overgffh )= @_;
  
  my ($atchr, $lastchr, $ingenes, $ng,$nx,$nr,$nsame,$nhit,$nwarn)= (0) x 10;
  my $more=1;
  while ($more) {
  
    if($SORTEDGFF) { # get first chr
      drop_overlaps(); # last ref
      ($atchr, $ng, $nx, $nwarn)= collect_overlaps_bychr($overgffh, 'NEXTCHR', $lastchr); 
      $ingenes= \@overgenes;
    } else {
      ($atchr, $ng, $ingenes)= copy_overlaps_bychr('NEXTCHR', $lastchr); # assumes all in @overgenes
    }
    unless(ref($ingenes) and scalar(@$ingenes) > 0) { $more=0; last; } #  warn "ERR: filtergff no ingenes\n";
    
    foreach my $generec (@$ingenes) {
      if(ref($generec)) {
        my $mrna= $generec->[0]; #?? or check type?
        # my($mrna)= grep { $$_[2] eq 'mRNA' } @$generec;
        my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid)= @$mrna;
        my $geneother= undef; #later
        
        $nsame += testgene( $generec, $geneother) ;
        $nx += scalar( grep { $$_[2] eq 'exon' } @$generec);
        $ng++;
        $lastchr= $ref;
      }
    }
  
  } # done batch
  
  return ($ng,$nx,$nsame);
}




# convert to array handling? _max(@list) 
sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }

sub drop_overlaps {
  # clear global memory for overlaps
  # my $overlaplist= {}; my @overgenes= (); my $overgeneid= 0;  my %overlapids=(); 
  @overgenes=(); $overgeneid=0;  $lastgeneid= 0;
  %overlapids=();
  for my $ref (keys %{$overlaplist}) {
    for my $ib (keys %{$overlaplist->{$ref}}) { delete $overlaplist->{$ref}{$ib}; }
    delete $overlaplist->{$ref};
  }
}

sub addgene {
  my($generecIN, $geneother)= @_;
  my @generec= sort _sortgene @$generecIN;
  my $generec= \@generec;

  # 30jul11: change to IDlist, @overgenes index : dont need, bug was elsewhere
  my $geneid= $overgeneid= scalar(@overgenes);
  push @overgenes, $generec;
  #?? push @overgeneother, $geneother; # needs to exist? same order as @overgenes ?
  #or $geneother{$gid}= $geneother || [];
  
    # add_overlap($generec);
  my $mrna= $generec->[0]; #?? or check type?
  my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid)= @$mrna;
  $overlapids{$gid}++; # not same as above..
  my @bins= (int($tb/$BINSIZE) .. int($te/$BINSIZE));
  foreach my $ib (@bins) { push( @{$overlaplist->{$ref}{$ib}}, $geneid); } # $generec
}

sub copy_overlaps_bychr { # upd17 
  my( $thischr, $lastchr)= @_;
  my( $ng,$nx,$nwarn)=(0) x 10;
  $thischr||="";  $lastchr||="";
  if($thischr eq 'NEXTCHR') {
    my @chr= sort keys %$overlaplist;
    if($lastchr) { 
      for my $i (1..$#chr) { if($chr[$i-1] eq $lastchr) { $thischr= $chr[$i]; last; } }
    } else {
      $thischr= $chr[0];
    }
  }
  my @generecs=(); my %gids;
  my @bins= sort{$a<=>$b} keys %{$overlaplist->{$thischr}};
  foreach my $ib (@bins) { map{ $gids{$_}++ } @{$overlaplist->{$thischr}{$ib}};  }
  # for my $id (sort{$a<=>$b} keys %gids) { push @generecs, $overgenes[$id]; $ng++; }
  my @ids= (sort{$a<=>$b} keys %gids); @generecs= @overgenes[ @ids ]; $ng= @generecs;
  $lastchr= $thischr;
  return ($lastchr,$ng,\@generecs);
}

sub collect_overlaps_bychr {
  my( $gffinh, $thischr, $lastchr)= @_;
  my( $ng,$nx,$nwarn)=(0) x 10;
  my @generec=(); my @geneother=();
  $thischr||="";  $lastchr||="";
  while(<$gffinh>) {
    next unless(/^\w/); chomp;
    my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr)= split"\t";
    
    if($thischr and $ref ne $thischr) {
      if($thischr eq 'NEXTCHR' and $ref ne $lastchr) { $thischr= $ref; }
      elsif($lastchr eq $thischr) { last; } 
      else { next; }
    }
    # last if($thischr and $ref ne $thischr and $lastchr eq $thischr); # return, finished thischr
    # next if($thischr and $ref ne $thischr);
    
    next if($SKIPREF and $ref =~ m/^($SKIPREF)/);
    # next if($passtypes and "$typ.$src" !~ m/$passtypes/); # pass other types
    #^ drop passtypes for mrnatypes,exontypes
    $lastchr= $ref;
    
    my $gid= undef; my($id,$pid)=(0) x 9; 
    $tattr ||="";
    if($tattr =~ m/\bID=([^;]+)/) {  $id=$1; $gid=$id if($typ =~ /^($mrnatypes)$/); }
    if($tattr =~ m/\bParent=([^;]+)/) {  $pid=$1; 
      $gid=$pid unless($gid and ($typ =~ /^($mrnatypes)$/));  }
    unless(defined $gid) { $gid = "N".$ng; }

    my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid]; 
    if($typ =~ /^($mrnatypes)$/) {  
      addgene(\@generec, \@geneother) if(@generec);
      $ng++; @generec = ($rloc); # maybe best store as array of [gffcols], sorted by type-loc
      @geneother=();
      
    } elsif($typ =~ /^($exontypes)$/) {
      my $rid= $generec[0]->[9] or 0;
      if(not $IgnoreOutOfOrder and $rid ne $pid) { 
        warn "#ERR: out of order gene record: mRNA=$rid, $typ=$pid\n" if($nwarn++ < 9); } # limit warns
      else { push @generec, $rloc;  $nx++; }
      
    } elsif( ! $passtypes or "$typ.$src" =~ m/$passtypes/ ) {
      push @geneother, $rloc;
    }
  }
  addgene(\@generec, \@geneother) if(@generec);
  
  $ngov += $ng; $nxov += $nx; # upd17: need to upd global now
  return ($lastchr, $ng, $nx, $nwarn);
}

sub collect_overlaps
{
  my($gffinh)= @_;
  my($lastchr, $ng, $nx, $nwarn)= collect_overlaps_bychr($gffinh,""); #all
  warn"#collect_overlaps ngene=$ng, nexon=$nx, lastchr=$lastchr, nerrOutOfOrder=$nwarn\n" if $debug;
  return ($ng, $nx, $nwarn, $lastchr);
}

sub OLDcollect_overlaps {
  my($gff)= @_;
  my ($ng,$nx,$nwarn)=(0) x 10;
  my @generec=();  my $gid;
  while(<$gff>) {
    next unless(/^\w/); chomp;
    my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr)= split"\t";
    if($passtypes and "$typ.$src" !~ m/$passtypes/) { next; } # pass other types
    #^ drop passtypes for mrnatypes,exontypes
    if($SKIPREF and $ref =~ m/^($SKIPREF)/) { next; }
    $tattr ||="";
    
    my($id,$pid)=(0,0); 
    if($tattr =~ m/\bID=([^;]+)/) {  $id=$1; $gid=$id if($typ =~ /^($mrnatypes)$/); }
    if($tattr =~ m/\bParent=([^;]+)/) {  $pid=$1; 
      $gid=$pid unless($gid and ($typ =~ /^($mrnatypes)$/));  }
    unless(defined $gid) { $gid = "N".$ng; }

    my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid]; 

    if($typ =~ /^($mrnatypes)$/) {  
      addgene(\@generec) if(@generec);
      $ng++; @generec = ($rloc); # maybe best store as array of [gffcols], sorted by type-loc
    
    } elsif($typ =~ /^($exontypes)$/) {
      my $rid= $generec[0]->[9] or 0;
      if(not $IgnoreOutOfOrder and $rid ne $pid) { warn "#ERR: out of order gene record: mRNA=$rid, $typ=$pid\n" if($nwarn++ < 9); } # limit warns
      else { push @generec, $rloc;  $nx++; }
    }
  }
  
  addgene(\@generec) if(@generec);
  
  warn"#collect_overlaps ngene=$ng, nexon=$nx, nerrOutOfOrder=$nwarn\n" if $debug;
  return ($ng, $nx);
}


__END__

=item example

 $td/overgenedup.pl -act markid -mark gdup dpx_epit27-augmap.gff dpx_ri27-augmap.gff | grep -c gdup= 
   7566

 $td/overgenedup.pl -type CDS  -act markid -mark gdup dpx_epit27-augmap.gff dpx_ri27-augmap.gff | grep -c gdup=
  16386  **

 $td/overgenedup.pl -act markid -mark gdup dpx_epit27-augmap.gff dpx_ti_esex27-augmap.gff | grep -c gdup= 
    139

 $td/overgenedup.pl -type CDS -act markid -mark gdup dpx_epit27-augmap.gff dpx_ti_esex27-augmap.gff | grep -c gdup=
   6321 **


grep -c 27.mRNA dpx_*27-augmap.gff
  dpx_epit27-augmap.gff:29183
  dpx_ri27-augmap.gff:28416
  dpx_ti_ecad27-augmap.gff:15730
  dpx_ti_esex27-augmap.gff:18970

=cut
