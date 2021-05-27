#!/usr/bin/env perl
# countblatalign.pl

# updates:
# parallel runs? on teragrid; 
# mods for query=cirad1asm non-align blocks?
#
# change to megablast? blat crawls on some multi-matches
# or use blatopts to speed up
# $nb/megablast -a2 -e1e-50 -FF -m8 -d blast/cacao41asm -i mars41s69k83.fa 


use strict;
use Getopt::Long;

my $bbin="/bio/bio-grid/mb/bin";

my $qtag="mars41";
my $qgenome="cacao41asm.fa";
my $tgenome="tcacao_cirad1asm.fa";

my $didhead=0;
my $debug=0;
my $blatopt=" -tileSize=12 -minMatch=4 -minScore=100  "; # for speed
my $fastaonly=0;
my $expand=0;

my $optok= GetOptions(
  "qname=s", \$qtag, 
  "qgenome=s", \$qgenome, 
  "tgenome=s", \$tgenome, 
  "biobin=s", \$bbin, 
  "blatopts=s", \$blatopt, 
  "fastaonly!", \$fastaonly, 
  "expandspan=i",\$expand,
  "debug!", \$debug, 
  );

die "usage: 
countblatalign.pl -qname=$qtag -qgenome $qgenome -tgenome $tgenome 
  < qregions.gff > qalign.table
  using -biobin=$bbin blat,faFrag -blatopt=$blatopt
  -fastaonly -expand=$expand
" unless($optok and -f $qgenome and -f $tgenome);

print "#blatalign $qtag regions to genomes same=$qgenome, other=$tgenome\n";

while(<>) { #nongap.gaps
  next unless(/^\w/); # and /\tgap/; option \tgap : region, any type?

  my($ref,$src,$typ,$b,$e,$w)=split;
  #? use $qtag=$src:$typ ?
  my $qlen= 1 + $e-$b;
  my ($nnn)= (m/nnn=(\d+)/) ? $1 : 0;
  
  # fixme faFrag uses b-1 start?
  my ($b0,$e0)= ($b-1, $e);
  if($expand>0) { ($b0,$e0)= ($b0 - $expand, $e0 + $expand); $b0=0 if($b0 < 0); }
  # ^ add to qloc?

  (my $na= $ref) =~ s/(caffold|uper_)//i;
  my $query= $qtag.$na."k".int($b/1000); # qname
  if($expand>0) { $query.="x".int($expand/100); }
  
  my $qloc="$qtag:$ref:$b-$e";
  if($expand>0) { $qloc.=":+-".$expand; }
  
  my $queryfa="$query.fa";
  
  # print "#blatalign $query, len=$qlen, $typ:$src, loc=$qloc to same=$qgenome, other=$tgenome\n" if $debug;
  print "#fasta $queryfa, len=$qlen, nnn=$nnn, $typ:$src, loc=$qloc\n" 
    if ($fastaonly or $debug);

  system("$bbin/faOneRecord $qgenome $ref | $bbin/faFrag stdin $b0 $e0 $queryfa > /dev/null 2>&1");
  
  if( $fastaonly ) {
  
  } elsif( -f $queryfa ) {
    my ($cmd,$res);
    
    $cmd="$bbin/blat $blatopt $qgenome $queryfa stdout | ";
    $res= countAlign( "same", $qloc, $qlen, $ref,$b,$e, $cmd);
    
    $cmd= "$bbin/blat $blatopt $tgenome $queryfa stdout | ";
    $res= countAlign( "other", $qloc, $qlen, $ref,$b,$e, $cmd);
    
    #?? rewrite queryfa header >ID $res ...
  }
  
  # unlink($queryfa); # keep around for other uses; annotate w/ align level
}



sub countAlign { # parse blat
  my($type,$qname,$qlen, $qref, $qb, $qe, $blatcmd)= @_;
  my $minmat= ($qlen > 7000) ? 800 : 400; # int($qlen / 10);
  my $fullmat= int( 0.95 * $qlen);
 
  my($nfull,$fmat,$fspan,$nmat,$smat,$sspan,$maxmat,$maxloc,$nb)= (0) x 10;
  # warn "#blatalign $type,$qname,$qlen,$blatcmd\n" if $debug;

  open(B, $blatcmd) or die "FAIL: $blatcmd\n";  # return; die/warn ?
  while( <B> ) {
    next unless(/^\d/);     
    chomp; my @v= split"\t";
    my( $mat, $mis, $rep_matches, $orient,
      $qid, $qsize, $qstart, $qend,
      $tid, $tsize, $tstart, $tend,
      $blocksizes, $qstarts, $tstarts,
      )= @v[0..2, 8..16, 18..20];
    $qstart++; $tstart++; # move to 1-origin
      
    $nb++; 
    last if($nb > 70 and ($nfull>0 or $nmat > 1)); #? need or not
    next if($mat < $minmat);
    
    my $span=($mat + $mis);
    if($mat >= $fullmat) { 
      $nfull++;
      $fmat += $mat;
      $fspan+= $span;
      if($type ne "same" and $mat>$maxmat) { $maxmat= $mat; $maxloc="$tid:$tstart-$tend"; }      
    } else {
      # next 
      if($type eq "same" and $qref eq $tid and $tstart < $qe and $tend > $qb) {
        # ^^ but have some part matches from spans w/ NNN.
        $nfull++;
        $fmat += $mat;
        $fspan+= $span;
      } else {
        $sspan += $span;
        $smat += $mat;
        if($mat>$maxmat) { $maxmat= $mat; $maxloc="$tid:$tstart-$tend"; }      
        $nmat++;
      }
    }
    
  } close(B);
  
  
  my $pquery= ($sspan>0) ? int(1000*$maxmat/ $qlen)/10 : 0; 
  my $pmat  = ($sspan>0) ? int(1000*$smat/$sspan)/10 : 0;
  my $pfull = ($fspan>0) ? int(1000*$fmat/$fspan)/10 : 0;
  
  my $hd= join("\t",qw(genome qlen maxaln nmatch bmatch pmatch  nfull bfull pfull qloc tloc));
  my $res= join("\t",$type,$qlen,$pquery,$nmat,$smat,$pmat,$nfull,$fmat,$pfull,$qname,$maxloc);
  print $hd,"\n" unless($didhead++);
  print $res,"\n";
  
  $res= "$type:$qname:$qlen,$pquery%:$maxloc";
  return $res;
}


__END__

    # blat parse @v for locations; check outside query if same, but not full align
    # from /bio/bio-grid/dpulex/prots/blat2gff.pl
#     my @blocksizes = split( /,/ , $blocksizes );
#     my @qstarts    = split( /,/ , $qstarts );
#     my @tstarts    = split( /,/ , $tstarts );
#     my $npart      = @qstarts;


==> ../work/cacao41asm_cirad1.align.nongap.gaps <==
super_35        marcir  gap     3301000 3339489 38490   .       .
super_44        marcir  gap     1499513 1522801 23289   .       .
super_39        marcir  gap     2088004 2104952 16949   .       .
super_4 marcir  gap     9942804 9959147 16344   .       .
super_28        marcir  gap     3775758 3791712 15955   .       .


# sample blat out: keep only high-match aligns

melon2.% $bg/mb/bin/blat tcacao_cirad1asm.fa cacao9s7k415.fa stdout | more
psLayout version 3

match   mis-    rep.    N's     Q gap   Q gap   T gap   T gap   strand  Q               Q       Q       Q       T               
T       T       T       block   blockSizes      qStarts  tStarts
        match   match           count   bases   count   bases           name            size    start   end     name            
size    start   end     count
---------------------------------------------------------------------------------------------------------------------------
------------------------------------
42      1       0       0       1       3       1       2       +       super_7:4149340-4159426 10086   9430    9476    sca
ffold04560      882910  417951  417996  2       36,7,   9430,9469,      417951,417989,
52      5       0       0       0       0       1       2       +       super_7:4149340-4159426 10086   9430    9487    sca
ffold02410      87164   86831   86890   2       31,26,  9430,9461,      86831,86864,
42      4       0       0       0       0       0       0       +       super_7:4149340-4159426 10086   8651    8697    sca
ffold03282      2083    203     249     1       46,     8651,   203,
61      5       0       0       0       0       1       1       +       super_7:4149340-4159426 10086   9430    9496    sca
ffold03074      458712  281799  281866  2       60,6,   9430,9490,      281799,281860,
50      3       0       0       1       148     1       150     +       super_7:4149340-4159426 10086   9436    9637    sca
ffold02637      1826201 343627  343830  2       36,17,  9436,9620,      343627,343813,
65      5       0       0       0       0       0       0       +       super_7:4149340-4159426 10086   9427    9497    sca
ffold02002      682909  354157  354227  1       70,     9427,   354157,
75      2       0       0       2       42      3       28      +       super_7:4149340-4159426 10086   9430    9549    sca
ffold01386      945499  493743  493848  4       9,32,31,5,      9430,9444,9513,9544,    493743,493756,493811,493843,
65      3       0       0       1       3       2       3       +       super_7:4149340-4159426 10086   9427    9498    sca
ffold01001      2267    552     623     3       19,4,45,        9427,9449,9453, 552,573,578,
32      0       0       0       0       0       0       0       +       super_7:4149340-4159426 10086   8694    8726    sca
ffold00794      337618  332441  332473  1       32,     8694,   332441,
222     22      0       0       1       1       0       0       +       super_7:4149340-4159426 10086   8524    8769    sca
ffold00709      695430  520221  520465  2       169,75, 8524,8694,      520221,520390,
60      2       0       0       1       204     1       204     +       super_7:4149340-4159426 10086   9445    9711    sca
ffold00598      362676  281376  281642  2       31,31,  9445,9680,      281376,281611,
37      3       0       0       0       0       0       0       +       super_7:4149340-4159426 10086   9584    9624    sca
ffold00447      620610  581188  581228  1       40,     9584,   581188,

keep>>
4720    287     0       0       32      4545    39      132922  +       super_7:4149340-4159426 10086   507     10059   sca
ffold01453      1266781 460845  598774  48      118,64,16,43,75,69,48,56,198,103,137,29,150,83,44,4,357,78,147,72,91,34,106
,97,47,58,7,42,68,8,4,627,33,17,64,99,199,216,227,72,271,27,380,95,99,66,54,8,  507,631,700,718,787,862,931,983,1079,1281,1
384,1526,1555,1705,1790,1834,1840,2197,2282,2434,2843,2934,2968,3087,3198,3249,3314,3321,6201,6863,6871,6876,7504,7538,7555
,7626,7725,7925,8414,8649,8727,9019,9046,9714,9829,9929,9995,10051,     460845,460976,461041,461060,461126,461203,461282,46
1330,464672,464873,464977,465127,468567,468718,468801,468846,468851,469212,469324,469483,469555,469647,469682,469796,469910
,469962,470026,470034,470077,470145,470162,470166,470793,470826,470844,470908,471010,471211,471434,471670,471748,475175,475
320,478370,478489,478590,598698,598766,
<< 
71      4       0       0       1       75      2       90      +       super_7:4149340-4159426 10086   5739    5889    sca
ffold01027      2987    108     273     3       27,28,20,       5739,5841,5869, 108,224,253,
66      7       0       0       0       0       0       0       +       super_7:4149340-4159426 10086   8694    8767    sca
