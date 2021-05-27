#!/usr/bin/env perl
# overlapeval.pl from overlapfilter.perl

=item usage

  perl overlapeval.pl -overlap predict.gff -input stdin|reference.gff 

  evaluate overlap match to input (gene models, other exon-like groups)
  
=cut

use strict;
use warnings;
use Getopt::Long;

use constant SAMEBASE => 0; # for _sameloc, slop allowed in loca == locb
use constant { ACT_DROP=>1, ACT_KEEP=>2, ACT_MARK=>3, ACT_MARK_WITH_ID=>4,
      ACT_MARK_OVERBASE=> 5};
use constant { kOVERLAP=>1, kSAMELOC=>2, kNEARLOC=>3, kCUTOVER=>4 };
use constant DOT => '.';

our $BINSIZE  = 1000 ; #was# 5000;
our $NEARDIST =  500; # was 15k; needs to be < BINSIZE

our $debug=1;
my ($overlaps,$overlaplist,$markidtype,$input,$itype,$action,$actid,$ok,$mark,$baseover);
my %overids;
my ($overtype,$typeover,$sametypes,$pctover,$passtypes)= (kOVERLAP,"",1,0,"");
my ($save_pctover,$save_baseover)=(0,0);
my ($save_miss,$n_overlaps)=(0,0);
my $UseStrand= 0; # or not?
$input= $overlaps="";

my $optok= GetOptions(
  "overlaps=s", \$overlaps, 
  #"typeover=s", \$typeover, 
  "pctover=i", \$pctover, 
  "strand!", \$UseStrand, 
  # "NEARDIST=i", \$NEARDIST, "BINSIZE=i", \$BINSIZE, 
  # "mark=s", \$mark, 
  "input=s", \$input,  
  #"midtype=s", \$markidtype, # return not ID= but other attribute or score
  "action=s", \$action, 
  "passtypes=s", \$passtypes,  
  "debug!", \$debug, 
  "baseover!", \$baseover, # show base,pct overlap counts per item and total
#  "itype=s", \$itype,  
  );

$action ||= "markbase";

die "
usage:  perl overlapeval -overlaps predict.gff  -input stdin|ref.gff 
" unless($optok and $action and (-f $overlaps or -f $input));

#         -act keep|drop|mark|markid|markbase -mark=terepeat 
#         -typeover overlap|sameloc|samefeat|near
#         -baseover : show overlap base count 
#         -neardist=$NEARDIST (for near type, base distance)
#         -passtypes='CDS' | 'mRNA,CDS' | 'gene.Gnomon,gene.Chainer' : act on these types only
#         -midtype=ID|Name|score|source|... (for markid, attribute or score to mark, ID default)
#         -itype=gff|blast -input 
#        -pctover=0 (for overlap type, min % overlap)

# if($action =~ /^cut/) { $typeover= "cut"; $action= "mark"; }
## for matching IDs b/n updates; need also match type: gene/mRNA/exon/CDS/...

# force pctover
$pctover= 1 unless($pctover);

$mark  ||= "terepeat";
$itype ||= "gff";
$actid= ($action =~ /keep/) ? ACT_KEEP : ($action =~ /mark/) ? ACT_MARK : ACT_DROP;
$actid= ACT_MARK_WITH_ID if($actid == ACT_MARK && $action =~ /id/i);
$actid= ACT_MARK_OVERBASE if($actid == ACT_MARK && $action =~ /base/i);

$overtype= ($typeover =~ /same/) ? kSAMELOC : ($typeover =~ /near/) ? kNEARLOC :
           ($typeover =~ /cut/) ? kCUTOVER : kOVERLAP;
$sametypes= ($typeover =~ /feat/) ? 1 : 0; #??
$pctover= $pctover/100.0 if($pctover);
$passtypes =~ s/[,]/\|/g;

$BINSIZE= int($NEARDIST*2) if($overtype == kNEARLOC and $NEARDIST > $BINSIZE);


my $ovh; 
   if($overlaps =~ /.gz$/) { $ok= open(OVR,"gunzip -c $overlaps |");  $ovh= *OVR; }
elsif($overlaps =~ /^(stdin|-)/) { $ovh= *STDIN; $ok=1; }
else { $ok= open(OVR,$overlaps); $ovh= *OVR; }
die "bad -overlaps=$overlaps" unless($ok);

$overlaplist= collect_overlaps($ovh); close($ovh);
die "bad -overlaps=$overlaps" unless(scalar(%$overlaplist));

my $inh= *STDIN;
$ok = ($input =~ /.gz$/) ? open($inh,"gunzip -c $input |") 
      : ($input =~ /^(stdin|-)/) ? $inh= *STDIN
      : open($inh,$input);
die "bad -input=$input" unless($ok);

my($sum_baseover,$sum_pctover,$sum_basetotal,$sum_pcttotal,$sum_miss)=(0) x 9;

my $nin=0;
my $nr=0;
if($itype =~ /blast/i) { die "not supported"; }
else { $nr= filter_gff($inh); }
warn"#overlaps found=$nr\n" if $debug;

# $sum_basetotal||=1;
# $sum_pcttotal||=1;
# warn "# base statistics: overlaps n=$nr , input n=$nin ,  overset n=$n_overlaps\n"
#   if ($baseover || $debug);
# if ($baseover) {  $nr||= 1;
# my $ave=sprintf("ave_baseover=%.3f, ave_pctover=%.3f, ave_miss=%.3f",
#         $sum_baseover/$sum_basetotal,$sum_pctover/$sum_pcttotal, $sum_miss/$nr);
# warn "# $ave 
# # sum_baseover=$sum_baseover, sum_pctover=$sum_pctover, sum_miss=$sum_miss 
# # sum_basetotal=$sum_basetotal, sum_pcttotal=$sum_pcttotal\n" 
# } ;

#..................

sub filter_gff
{
  my($inh)= @_;
  my $nover=0;
  my %overref;
  my($exonbase,$exonmax, $exonmiss)=(0,0,0);
  my (%genes, $geneid, $hasgene); 
  $geneid="none";
  while(<$inh>){
    unless(/^\w/){ next if(/^(#n |$)/); next; } # print and 
    # my $line=$_;
    my($ref,$src,$typ,$tb,$te,$tp,$to,$tx,$tattr)= split"\t";
  
    my($cid)= m/\bID=([^;]+)/;
    my($pid)= m/\bParent=([^;]+)/;
    if($typ eq "mRNA") { $geneid= $pid || "none";  $hasgene=($geneid ne "none")?1:0; } # special case; hope it has pid and data sorted
    if($passtypes and "$typ.$src" !~ m/$passtypes/) { next; } # print; pass other types

    $nin++; # exons w/ dup gene/tr ?? dont count twice ..
    $pid ||= $cid || $nin;
    unless($hasgene) { $geneid= (m/Name=([^;:\s]+)/) ? $1 : "none"; }
    
    #? pick longest tr among alternates?
    my $tw= 1+$te-$tb;
    my %pid= map{ $genes{$geneid}{$_} += $tw; $_,1 }split",",$pid; # many?
    my @pid= sort keys %pid;
    
    my %overstat = overlaps( $ref, $tb, $te, $to);
    # overids per id: baseover, missover, pctover?
    
      # transfer to input model pid
    my $xbase= 1 + $te-$tb;
    my ($xmax,$xmiss)= (0,0);
    foreach my $pid (@pid) {
      $overref{$pid}{self}{hit} += $xbase;
      $overref{$pid}{self}{miss} = 0;
      my($maxhit,$maxmiss)=(0,0);
      foreach my $oid (sort keys %overstat) {
        my $hit = $overstat{$oid}{hit};
        my $miss= $overstat{$oid}{miss};
        $overref{$pid}{$oid}{hit}  += $hit;
        $overref{$pid}{$oid}{miss} += $miss;
        if( $hit > $maxhit or ($hit == $maxhit and $miss < $maxmiss)) {
          $maxhit= $hit; $maxmiss= $miss;
        }  
      }
      $overref{$pid}{max}{hit}  += $maxhit;
      $overref{$pid}{max}{miss} += $maxmiss;
      $xmax= $maxhit; $xmiss=$maxmiss;
    }
    $nover++ if ($xmax/$xbase >= $pctover); #?
    $exonmax+=$xmax; $exonmiss+= $xmiss; $exonbase+= $xbase;
  }
  
  my %gene1id=(); # "best"/longest tr per gene
  foreach my $gid (keys %genes) {
    next if($gid eq "none");
    my @trid= keys %{$genes{$gid}};
    my($t1id, $glen)=(0,0);
    foreach my $trid (@trid) {
      my $tlen= $genes{$gid}{$trid};
      if($tlen > $glen) { $glen= $tlen; $t1id= $trid; }
    }
    $gene1id{$t1id}= $gid;
  }
  
  my $alloverbase= 0;
  foreach my $oid (keys %overids) { $alloverbase += overbases($oid); }
  
  my($nperf, $genehit, $genemiss, $genebase,
     $tophit,$topmiss,$allhit,$allmiss,$allbase,$ng,$nghit)= (0) x 30;
  foreach my $pid (sort keys %overref) {  
    my @oid= sort keys %{$overref{$pid}};
    my $gbase= $overref{$pid}{self}{hit};
    my $gbest= ($gene1id{$pid}) ? 1: 0;
    my($besthit,$bestmiss)=(0,0);
    foreach my $oid (@oid) {
      next if ($oid eq "self" or $oid eq "max"); # also max ?
      my $hit = $overref{$pid}{$oid}{hit};  # true pos
      my $miss= $overref{$pid}{$oid}{miss}; # false pos per exon, also get FP per oid
      my $obase= overbases($oid);
      my $FP= _max(0, $obase - $hit);  ## or $miss ? false pos ??
      # ^^ this is much larger than miss at exon matching
      if($hit > $besthit or ($hit == $besthit and $FP < $bestmiss)) {
        $besthit= $hit; $bestmiss= $FP;
      }
    }
    # besides all, want to know best gene model match
    $tophit += $besthit;
    $topmiss+= $bestmiss;
    $nperf++ if( $besthit/$gbase >= 0.95 and $bestmiss/$gbase <= 0.5);
    $allbase+= $gbase;
    $allhit += $overref{$pid}{max}{hit};  # adds to more than $alloverbase, dupl locations/counts ?
    $allmiss+= $overref{$pid}{max}{miss};

    $genebase+= $gbase if $gbest;
    $genehit += $besthit if $gbest;
    $genemiss+= $bestmiss if $gbest;
    $ng++;  $nghit++ if($besthit/$gbase >= $pctover);
  }
  
  print "ngene=$ng; genehit=$nghit; gperfect=$nperf; nexon=$nin; exonhit=$nover\n";
  return $nover unless($allbase);
  
  my($sens,$spec,$spec2,$avss);
  # sens = TP/(TP+FN);  spec = TP/(TP+FP);
  $sens=int(1000 * $allhit/$allbase)/10;
  $spec=int(1000 * $allhit/(1+$allhit+$allmiss))/10; # FP not right here
  $spec2=int(1000 * $allhit/$alloverbase)/10; # ?? right; not quite: allhit > alloverbase dupl locations?
      # spec2 should be sum(hits-per-overid) / alloverbase
      # other counts are based on all transcripts (w/ dupl. exons)
      # genehit, genemiss, genebase, alloverbase should be counting locations once only
      # ditto for exonmax,exonbase
  $avss= ($sens+$spec)/2;
## this one not useful:
##  print "All: Sn=$sens; Sp=$spec; Sp2=$spec2; Ave=$avss; hit=$allhit; miss=$allmiss; allbase=$allbase; allover=$alloverbase\n";

  # exon-only Sn/Sp
  if($exonbase>0) {
  $sens=int(1000 * $exonmax/$exonbase)/10;
  $spec=int(1000 * $exonmax/(1+$exonmax+$exonmiss))/10; # FP not right here
  $spec2=int(1000 * $exonmax/$alloverbase)/10; # ?? right
  $avss= ($sens+$spec)/2;
  print "Exon: Sn=$sens; Sp=$spec; Sp2=$spec2; Ave=$avss; hit=$exonmax; miss=$exonmiss; exonbase=$exonbase; allover=$alloverbase\n";
  }
  
  $sens=int(1000 * $tophit/$allbase)/10;
  $spec=int(1000 * $tophit/(1+$tophit+$topmiss))/10; # FP not right here
  $spec2=int(1000 * $tophit/$alloverbase)/10; # ?? right
  $avss= ($sens+$spec)/2;
  print "BestTr: Sn=$sens; Sp=$spec; Sp2=$spec2; Ave=$avss; hit=$tophit; miss=$topmiss; allbase=$allbase;\n";

  if($genebase>0) {
  $sens=int(1000 * $genehit/$genebase)/10;
  $spec=int(1000 * $genehit/(1+$genehit+$genemiss))/10; # FP not right here
  $spec2=int(1000 * $genehit/$alloverbase)/10; # ?? right
      # genehit, genemiss, genebase, alloverbase should be counting locations once only
  $avss= ($sens+$spec)/2;
  print "Gene: Sn=$sens; Sp=$spec; Sp2=$spec2; Ave=$avss; hit=$genehit; miss=$genemiss; genebase=$genebase;\n";
  }
  
  return $nover;
}


# sub filter_blast # ncbi format=8,9 blast table

my $warns=0;

# convert to array handling? _max(@list) 
sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }
sub _samestrand { my($a,$b)=@_; return ($b eq $a || $b eq DOT || $a eq DOT || $b eq "" || $a eq "")?1:0; }
sub _sort_over { return ($a->[0] <=> $b->[0]) || ($b->[1] <=> $a->[1]);} # @[b,e] min-b, max-e


sub overlaps
{
  my($ref,$tb,$te,$to)= @_;    

  $save_baseover= $save_pctover= 0;
  $sum_basetotal += $te - $tb + 1;
  $sum_pcttotal  += 1;
  $save_miss = 0;
  my($missb, $misse)=(undef,undef);
  my %didid=();
  my %overstat=();
  
  return %overstat unless($overlaplist->{$ref});
  my @bins= (int($tb/$BINSIZE) .. int($te/$BINSIZE));
  foreach my $ib (@bins) {
    $overlaplist->{$ref}{$ib} or next;
    my @locs= @{$overlaplist->{$ref}{$ib}};
    foreach my $rloc (@locs) {
      my ($lb,$le,$lo,$lid,$oid)= @{$rloc}[0,1,4,3,6];
      next if($didid{$oid.$lb.$le}++);
      
      my $over= ($tb <= $le && $te >= $lb) ? 1 : 0;
      $over= 0 if($over && $UseStrand && ! _samestrand($to,$lo) );
      
      ## ?set overstat{miss} if not over?

      if($over) {   #  and $pctover
        my ($bb,$be)= ( _max($tb,$lb), _min($te,$le) );
        
    # change for eval: save per lid: baseover==maxo, missb+misse
        my $maxo= 1 + abs($be - $bb);
        # my $misso=  abs($tb - $lb) + abs($te - $le);
        my $misso= _max(0, $tb - $lb) + _max(0, $le - $te); # amount OUTSIDE ?

        my %lid= map{$_,1}split",",$lid; # many?
        foreach my $d (sort keys %lid) {
          $overstat{$d}{hit}  += $maxo;        
          $overstat{$d}{miss} += $misso; #?? should be false-pos only, outside of tb,te
          }
        }
      }
    }

  return %overstat;    
}


# sub nearloc

# sub sameloc
# {
#   my($ref,$tb,$te,$to,$typ)= @_;   
#   my @lid;
#   return 0 unless($overlaplist->{$ref});
#   $to ||= '+';
#   ($tb,$te,$to)= ($te,$tb,'-') if($tb>$te);
#   my @bins= (int($tb/$BINSIZE) .. int($te/$BINSIZE));
#   foreach my $ib (@bins) {
#     $overlaplist->{$ref}{$ib} or next;
#     my @locs= @{$overlaplist->{$ref}{$ib}};
#     foreach my $rloc (@locs) {
#       my ($lb,$le,$lid,$lo,$ltyp)= @{$rloc}[0,1,3,4,5];
#       next if($sametypes and not($typ eq $ltyp and $to eq $lo)); #?
#       push @lid, $lid if(abs($tb-$lb) <= SAMEBASE and abs($te-$le) <= SAMEBASE);
#       }
#     }
#   if(@lid) { my %lid= map{$_,1}@lid; return join",", sort keys %lid; }
#   return 0;
# }

sub overbases
{
  my($oid)= @_;
  my $nb=0;
  my @locs= @{$overids{$oid}};
  foreach my $x (@locs) {
    my($tb,$te,$ref,$gid)= @$x;
    $nb += 1 + $te-$tb;
  }
  return $nb;
}

sub collect_overlaps
{
  my($gff)= @_;
  my %overlaps=(); my $nr=0;
  while(<$gff>){
    next unless(/^\w/); chomp;
    my($ref,$src,$typ,$tb,$te,$tp,$to,@gffmore)= split"\t";
    
    if($passtypes and "$typ.$src" !~ m/$passtypes/) { next; } # pass other types

    ## save these always?
    my($cid)= m/\bID=([^;]+)/;
    my($pid)= m/\bParent=([^;]+)/;
    $pid ||= $cid;
    
    $nr++;
    my $oid= "N".$nr; # save always for uniq feature test
    my($gid); # fixme: separate $gid from markid : want two fields, one for id tests
    if($markidtype) {
      if($markidtype =~ /^score/i) { $gid=$tp; }
      elsif($markidtype =~ /^source/i) { $gid=$src; }
      elsif($markidtype =~ /^type/i) { $gid=$typ; }
      elsif($gffmore[-1] =~ m/\b$markidtype=([^;]+)/) { $gid=$1; }
    } else {
      $gid= $pid;
      #if($gffmore[-1] =~ m/\bID=([^;]+)/) {  $gid=$1; }
      #elsif($gffmore[-1] =~ m/\bParent=([^;]+)/) {  $gid=$1; }
    }
    unless(defined $gid) { $gid = $oid; }
    $pid ||= $oid;
    
    my $rloc= [$tb,$te,$ref,$gid,$to,$typ,$oid,$pid,$cid]; # change to string; save mem

    my %pid= map{$_,1}split",",$pid; # many?
    foreach my $p (sort keys %pid) { push( @{$overids{$p}}, $rloc); }
    my @bins= (int($tb/$BINSIZE) .. int($te/$BINSIZE));
    foreach my $ib (@bins) {  push( @{$overlaps{$ref}{$ib}}, $rloc); }
    }
  $n_overlaps=$nr; # global
  warn"#collect_overlaps=$nr\n" if $debug;
  return \%overlaps;
}

__END__

=item results for rgasp

melon.% cat  dmel-3L.exons.gff | grep exon | sort -k4,4n -k5,5nr -k2,2 | $td/overlapeval.pl -pass exon -pct 50 -in
dmel-3L.pasavalid.exons.gff -over stdin
#collect_overlaps=12758
ngene=2730; genehit=2676; gperfect=1621; nexon=17385; exonhit=16778
all: Sn=91.9; Sp=96.1; hit=8312432; miss=329334; allbase=9037063;
top: Sn=90.3; Sp=96.8; hit=8164844; miss=266200; allbase=9037063;

melon.% cat pasa_dmel.3L.exons.gff | grep exon | sort -k4,4n -k5,5nr -k2,2 | $td/overlapeval.pl -pass exon -pct 50
-in dmel-3L.pasavalid.exons.gff -over stdin
#collect_overlaps=19941
ngene=2730; genehit=1913; gperfect=1309; nexon=17385; exonhit=15013
all: Sn=79.3; Sp=94.5; hit=7172666; miss=415678; allbase=9037063;
top: Sn=64; Sp=98.3; hit=5790555; miss=97122; allbase=9037063;
#overlaps found=15013
 
melon.% cat $rgasp/gmaps/dm*.asm10s8.gff | grep exon | sort -k4,4n -k5,5nr -k2,2 | $td/overlapeval.pl -pass exon -p
ct 50 -in dmel-3L.pasavalid.exons.gff -over stdin
#collect_overlaps=156362
ngene=2730; genehit=1922; gperfect=463; nexon=17385; exonhit=14315
all: Sn=98.3; Sp=69.2; hit=8884780; miss=3953708; allbase=9037063;
top: Sn=83.4; Sp=62; hit=7544596; miss=4608582; allbase=9037063;
#overlaps found=14315


melon.% cat $rgasp/gmaps/dmKc167.asm10s8.gff | grep exon | sort -k4,4n -k5,5nr -k2,2 | $td/overlapeval.pl -pass exo
n -pct 50 -in dmel-3L.pasavalid.exons.gff -over stdin
#collect_overlaps=43126
ngene=2730; genehit=1684; gperfect=287; nexon=17385; exonhit=12777
all: Sn=76.9; Sp=70.2; hit=6954652; miss=2949868; allbase=9037063;
top: Sn=68.2; Sp=63; hit=6167395; miss=3610219; allbase=9037063;
#overlaps found=12777

melon.% cat $rgasp/gmaps/dmML.asm10s8.gff | grep exon | sort -k4,4n -k5,5nr -k2,2 | $td/overlapeval.pl -pass exon -
pct 50 -in dmel-3L.pasavalid.exons.gff -over stdin
#collect_overlaps=32878
ngene=2730; genehit=1529; gperfect=305; nexon=17385; exonhit=11281
all: Sn=70.3; Sp=77.4; hit=6360750; miss=1855963; allbase=9037063;
top: Sn=64; Sp=74.5; hit=5790840; miss=1972489; allbase=9037063;
#overlaps found=11281

melon.% cat $rgasp/gmaps/dmS2all.asm10s8.gff | grep exon | sort -k4,4n -k5,5nr -k2,2 | $td/overlapeval.pl -pass exo
n -pct 50 -in dmel-3L.pasavalid.exons.gff -over stdin
#collect_overlaps=50127
ngene=2730; genehit=1724; gperfect=314; nexon=17385; exonhit=12987
all: Sn=82.4; Sp=69.2; hit=7446582; miss=3305520; allbase=9037063;
top: Sn=72.6; Sp=57.5; hit=6561015; miss=4846092; allbase=9037063;
#overlaps found=12987

#.. more versions of samasm : 10.8, 11.0 (gsignals)

# pasa EST
melon.% cat pasa_dmel.3L.exons.gff | $td/overlapeval.pl -pass exon -pct 50 -in dmel-3L.exons.gff -over stdin -strand
#collect_overlaps=19941
ngene=4099; genehit=2513; gperfect=1421; nexon=12758; exonhit=9623
all: Sn=69.8; Sp=84.2; hit=6348097; miss=1188626; allbase=9087529;
top: Sn=56.1; Sp=83.9; hit=5103774; miss=973006; allbase=9087529;
#overlaps found=9623

# samasm 10.5?
melon.% cat $rgasp/gmaps/dm*.asm10s8.gff | $td/overlapeval.pl -pass exon -pct 50 -in dmel-3L.exons.gff -over stdin -strand
#collect_overlaps=156362
ngene=4099; genehit=2850; gperfect=962; nexon=12758; exonhit=9593
all: Sn=94.3; Sp=68; hit=8577304; miss=4035545; allbase=9087529;
top: Sn=80.2; Sp=47.5; hit=7290636; miss=8032411; allbase=9087529;
#overlaps found=9593

# samasm 10.8 (~ final w/o genefinder signals)
melon.% cat $rgasp/gmaps/dm*.asm10s9.gff | $td/overlapeval.pl -pass exon -pct 50 -in dmel-3L.exons.gff -over stdin -strand
#collect_overlaps=133018
ngene=4099; genehit=2884; gperfect=1155; nexon=12758; exonhit=9661
all: Sn=96.9; Sp=61.6; hit=8811965; miss=5473502; allbase=9087529;
top: Sn=83.2; Sp=59.7; hit=7562100; miss=5090747; allbase=9087529;
#overlaps found=9661

## CME only

# samasm 10.5?
melon.% cat $rgasp/gmaps/dmCME.asm10s8.gff | $td/overlapeval.pl -pass exon -pct 50 -in dmel-3L.exons.gff -over stdin -strand
#collect_overlaps=30231
ngene=4099; genehit=2152; gperfect=630; nexon=12758; exonhit=6769
all: Sn=64.4; Sp=77.8; hit=5858790; miss=1666114; allbase=9087529;
top: Sn=58.6; Sp=64.3; hit=5331363; miss=2956609; allbase=9087529;
#overlaps found=6769

# samasm 10.8 (~ final w/o genefinder signals)
melon.% cat $rgasp/gmaps/dmCME.asm10s9.gff | $td/overlapeval.pl -pass exon -pct 50 -in dmel-3L.exons.gff -over stdin -strand
#collect_overlaps=27526
ngene=4099; genehit=2187; gperfect=699; nexon=12758; exonhit=6868
all: Sn=66.4; Sp=71.9; hit=6037660; miss=2354464; allbase=9087529;
top: Sn=60.6; Sp=71.7; hit=5511048; miss=2168288; allbase=9087529;
#overlaps found=6868

# samasm 11.0 (genefinder signals, limited) << Score drops from samasm 10.8
melon.% cat $rgasp/gmaps/dmCME.asm11.gff | $td/overlapeval.pl -pass exon -pct 50 -in dmel-3L.exons.gff -over stdin -strand
#collect_overlaps=31054
ngene=4099; genehit=2196; gperfect=645; nexon=12758; exonhit=6903
all: Sn=64.5; Sp=78.2; hit=5867950; miss=1635606; allbase=9087529;
top: Sn=59.2; Sp=63.9; hit=5381169; miss=3034219; allbase=9087529;
#overlaps found=6903


## Per Cell-line, samasm 10.8

# CME (pe)  : Highest Spec, low Sens
melon.% cat $rgasp/gmaps/dmCME.asm10s9.gff | $td/overlapeval.pl -pass exon -pct 50 -in dmel-3L.exons.gff -over stdin -strand
#collect_overlaps=27526
ngene=4099; genehit=2187; gperfect=699; nexon=12758; exonhit=6868
all: Sn=66.4; Sp=71.9; hit=6037660; miss=2354464; allbase=9087529;
top: Sn=60.6; Sp=71.7; hit=5511048; miss=2168288; allbase=9087529;
#overlaps found=6868

# Kc167 (all reads, pe, sr)
melon.% cat $rgasp/gmaps/dmKc167.asm10s9.gff | $td/overlapeval.pl -pass exon -pct 50 -in dmel-3L.exons.gff -over stdin -strand
#collect_overlaps=37738
ngene=4099; genehit=2459; gperfect=732; nexon=12758; exonhit=8202
all: Sn=74.3; Sp=64.8; hit=6759719; miss=3662516; allbase=9087529;
top: Sn=66; Sp=57.8; hit=5999543; miss=4365174; allbase=9087529;
#overlaps found=8202

# ML (pe)
melon.% cat $rgasp/gmaps/dmML.asm10s9.gff | $td/overlapeval.pl -pass exon -pct 50 -in dmel-3L.exons.gff -over stdin -strand
#collect_overlaps=27856
ngene=4099; genehit=2211; gperfect=726; nexon=12758; exonhit=6981
all: Sn=66.7; Sp=72.5; hit=6069493; miss=2293453; allbase=9087529;
top: Sn=60.5; Sp=68.7; hit=5505514; miss=2506907; allbase=9087529;
#overlaps found=6981

# S2 (all reads, pe, sr, lr)  : Highest Sens, lowest Spec
melon.% cat $rgasp/gmaps/dmS2all.asm10s9.gff | $td/overlapeval.pl -pass exon -pct 50 -in dmel-3L.exons.gff -over stdin -strand
#collect_overlaps=39898
ngene=4099; genehit=2541; gperfect=775; nexon=12758; exonhit=8353
all: Sn=78.1; Sp=62; hit=7103720; miss=4345303; allbase=9087529;
top: Sn=70.8; Sp=56.7; hit=6435566; miss=4912267; allbase=9087529;
#overlaps found=8353


## Test combined group supermodels versus separately:

# S2 only, first 2MB
melon.% cat $rgasp/drosasm10/dmS2all/dm*_3Lb1.asm10s.sugff | $td/overlapeval.pl -strand -pass exon -pct 50 -in dmel-3L.exons.gff -over stdin
#collect_overlaps=5369
ngene=4099; genehit=326; gperfect=81; nexon=12758; exonhit=1025
all: Sn=10.3; Sp=53.3; hit=941609; miss=824192; allbase=9087529;
top: Sn=9.3; Sp=51.4; hit=851729; miss=804637; allbase=9087529;
#overlaps found=1025

# Mixed group Supermodels
# note: this now has ~same no. exons as separate groups; not merging much ..

melon.% cat $rgasp/drosasm10/dmAll/dm*_3Lb1.asm10s.sugff | $td/overlapeval.pl -strand -pass exon -pct 50 -in dmel-3L.exons.gff -over stdin
#collect_overlaps=18113
ngene=4099; genehit=355; gperfect=129; nexon=12758; exonhit=1165
all: Sn=12.4; Sp=55.5; hit=1131429; miss=906954; allbase=9087529;
top: Sn=10.5; Sp=43.4; hit=962692; miss=1251143; allbase=9087529;
#overlaps found=1165

# Separate, all group supermodels  << slightly better than mixed group sumodels

melon.% cat $rgasp/drosasm10/{dmML,dmKc167,dmS2all,dmCMEj}/dm*_3Lb1.asm10s.sugff | $td/overlapeval.pl -strand -pass exon -pct 50 -in dmel-3L.exons.gff -over stdin
#collect_overlaps=18190
ngene=4099; genehit=354; gperfect=130; nexon=12758; exonhit=1164
all: Sn=12.9; Sp=56; hit=1174987; miss=920171; allbase=9087529;
top: Sn=11; Sp=55.2; hit=1003162; miss=812453; allbase=9087529;
#overlaps found=1164

# 21aug09: matecov2gene rough models:

melon.% cat matecov3genes.3L.gff | $td/overlapeval.pl -strand -pass exon -pct 50 -in $dpa/dmel-3L.exons.gff -over stdin
#collect_overlaps=10872
ngene=4099; genehit=2194; gperfect=461; nexon=12758; exonhit=6754
All: Sn=59; Sp=87.5; Ave=73.25; hit=5369073; miss=766959; allbase=9087529;
Exon: Sn=54.4; Sp=85.8; Ave=70.1; hit=3390007; miss=557951; allbase=6229845;
BestTr: Sn=56; Sp=76.6; Ave=66.3; hit=5091436; miss=1549019; allbase=9087529;
#overlaps found=6754

# versus pasa_dmel as above
#   #collect_overlaps=19941
#   ngene=4099; genehit=2513; gperfect=1421; nexon=12758; exonhit=9623
#   All: Sn=69.8; Sp=84.2; Ave=77; hit=6348097; miss=1188626; allbase=9087529;
#   Exon: Sn=66.8; Sp=83; Ave=74.9; hit=4167383; miss=849154; allbase=6229845;
#   BestTr: Sn=56.1; Sp=83.9; Ave=70; hit=5103774; miss=973006; allbase=9087529;



=cut 

