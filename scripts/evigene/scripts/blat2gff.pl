#!/usr/bin/env perl
# blat2gff.pl 

use strict;
use Getopt::Long;

my ($swap,$addmatch);
my $matchtype=undef;
my $matchpart="match_part"; # no other SO choice ?
my $source="BLAT";
my %matchids;
my $prefixmatch="_mid";
my $gffver=3;
my $percentscore= 1; # default yes?

my $optok= GetOptions(
  "swapquerytarget!",\$swap, 
  "source=s", \$source, 
  "version:i", \$gffver, 
  "matchtype:s", \$matchtype, 
  "prefixmatch=s", \$prefixmatch, 
  "percentscore!", \$percentscore,
);

die "
usage: blat2gff < blat.psl > blat.gff
options: -source $source -matchtype EST_match -swapQueryTarget -nopercentscore -version [3|2]
" unless($optok);
  
if(defined $matchtype) { $addmatch=1; $matchtype="match" unless($matchtype); }
print "##gff-version $gffver\n";
my $ispsl=0;

while(<>) {
  $ispsl=1 ; ## if m/^psLayout/; # should do but want some slack here?
  next unless(/^\d/);
  chomp; 
  my @v= split "\t";
  if(@v==22) { shift(@v);} # ucsc psl starts with a 'bin' field?
  unless($ispsl and @v==21){ die "# error: doesnt look like psl format I know" ; }

  my( $matchscore, $mismatches, $rep_matches, $orient,
    $qid, $qsize, $qstart, $qend,
    $tid, $tsize, $tstart, $tend,
    $blocksizes, $qstarts, $tstarts,
    )= @v[0..2, 8..16, 18..20];
  
  $qstart++; $tstart++; # move to 1-origin
  if($swap) {
    ($tid, $tsize, $tstart, $tend,  $qid, $qsize, $qstart, $qend,)
      = ($qid, $qsize, $qstart, $qend, $tid, $tsize, $tstart, $tend,);
    ($qstarts,$tstarts)= ($tstarts,$qstarts);
  }
		
  my $matchid    = $qid .$prefixmatch. ++$matchids{$qid};
  my @blocksizes = split( /,/ , $blocksizes );
  my @qstarts    = split( /,/ , $qstarts );
  my @tstarts    = split( /,/ , $tstarts );
  my $npart      = @qstarts;

  $matchscore    = sprintf "%.2f",	(100 * ( $matchscore ) / $qsize) if($percentscore); 
  ## BioPerl Tools/Blat.pm does this, which seems wrong:  $matchscore + $mismatches + $rep_matches 
  
  if($addmatch) {
    ## with - orient need to back-calculate Query location (gff target) 
    if($orient eq '-') { ($qstart,$qend)= ($qend,$qstart); } # print reversed for Target
    gffput( 0, $tid, $source, $matchtype, $tstart, $tend, $matchscore, $orient, ".",
      $matchid, "$qid $qstart $qend");
    $matchscore="."; # dont duplicate in match_parts; score is for entire feature
    }
  for(my $p=0; $p<$npart; $p++) {
    my($qstart,$tstart,$len)= ($qstarts[$p], $tstarts[$p], $blocksizes[$p]);
    my $qend= $qstart+$len;
    my $tend= $tstart+$len;
    $qstart++; $tstart++; # move to 1-origin
    if($orient eq '-') { # reverse all exons for Target
      $len   = $blocksizes[$npart-1-$p]; 
      $qend  = $qstarts[$npart-1-$p]; # lower
      $qstart= $qend + $len; # higher
      $qend++;
      } 
    gffput( 1, $tid, $source, $matchpart, $tstart, $tend, $matchscore, $orient, ".",
      $matchid, "$qid $qstart $qend");
  }
  
}

sub gffput {
  my($part,$r,$s,$t,$b,$e,$p,$o,$f,$id,$tg)=@_;
  if($gffver<3){
    my $IDtag="gene" ;
    print join("\t", $r, $s, $t, $b, $e, $p, $o, $f, "$IDtag \"$id\"; Target \"$tg\""),"\n";
  }else{
    my $IDtag=($part) ? "Parent" : "ID";
    print join("\t", $r, $s, $t, $b, $e, $p, $o, $f, "$IDtag=$id;Target=$tg"),"\n";
  }
}

=head1 ABOUT blat2gff

Convert blat psl to gff v3, preserving exon/match_part info, and
separate distinct matches for same query id onto target.  This is
the usual case for EST/mRNA matches onto genome backbone, with multiple
duplicate but distinct genes, each with exon structure.  Blat provides
this info and it is important for genome analyses.

Note that bioperl SearchIO/psl.pm doesn't do this, it just uses the psl
start,end fields, and it lumps any same query-id into one common match group,
even for distinct locations (separate rows from psl output).

I missed the other BioPerl Bio/Tools/Blat.pm which does the qstarts,tstarts exon
calculations. It however gives a common ID and superfeature to all same-query-id
matches, which is a problem (same as below combining multiple gene hits into one).
Someone should rewrite BioPerl SearchIO/psl.pm to use the Tools/Blat.pm package.


=head1 USAGE

 $td/blat2gff.pl < dgri-gnoest.blat > dgri-gnoest.blat2gff

=head1 AUTHOR
  
  don gilbert, 2007, gilbertd@indiana.edu

=item Blat/PSL header

  0..8
  match   mis-    rep.    N's     Q gap   Q gap   T gap   T gap   strand 
          match   match           count   bases   count   bases        
  
  9..16
  Q               Q       Q       Q       T               T       T       T       
  name            size    start   end     name            size    start   end     
  
  17..20
  block   blockSizes      qStarts  tStarts
  count

Parse qStarts, tStarts for match_part structure

=item examples

  blat -fastMap  $em/dgri1/dgri_caf060210.fa dgri-est.fa dgri-gnoest.blat
  cat dgri-gnoest.blat | sort -k14,14 -k10,10 -k16,16n | more

  grep EB637050.1 dgri-gnoest.blat
  563     0       0       0       0       0       3       1912    -       
  gi|93002893|gb|EB637050.1  830 72  635    scaffold_15110  24565398  10715419  10717894        
  4       162,231,97,73,  195,357,588,685,        10715419,10716978,10717307,10717821,

  grep EB637050.1 dgri-gnoest.blat | perl -pe's/gi.\d+.(gb.\w+).\d+.//;' | $td/blat2gff.pl -match EST_match

  scaffold_15110  BLAT    EST_match   10715419   10717894  563  -  .   ID=EB637050_mid1;Target=EB637050 72 635
  scaffold_15110  BLAT    match_part  10715419   10715581  563  -  .   Parent=EB637050_mid1;Target=EB637050 195 357
  scaffold_15110  BLAT    match_part  10716978   10717209  563  -  .   Parent=EB637050_mid1;Target=EB637050 357 588
  scaffold_15110  BLAT    match_part  10717307   10717404  563  -  .   Parent=EB637050_mid1;Target=EB637050 588 685
  scaffold_15110  BLAT    match_part  10717821   10717894  563  -  .   Parent=EB637050_mid1;Target=EB637050 685 758

  # example tandem match pair
  grep EB634440.1 dgri-gnoest.blat  
  
  670   ..    +  gi|93000283|gb|EB634440.1|EB634440  753  0  683     scaffold_14830  6267026 2489511 2490194 2   398,272,        0,411,  2489511,2489922,
  683   ..    -  gi|93000283|gb|EB634440.1|EB634440  753  0  683     scaffold_14830  6267026 2484805 2485488 1   683,    70,     2484805,
  646   ..    -  gi|93000283|gb|EB634440.1|EB634440  753  0  683     scaffold_14830  6267026 2480511 2481194 3   272,71,305,     70,355,448,     2480511,2480796,2480889,

  grep EB634440.1 dgri-gnoest.blat | perl -pe's/gi.\d+.(gb.\w+).\d+.//;' | \
   $td/blat2gff.pl -match EST_match

  ##gff-version 3
  scaffold_14830  BLAT    EST_match    2489512 2490194 670  +  .  ID=EB634440_mid1;Target=EB634440 1 683
  scaffold_14830  BLAT    match_part   2489512 2489909 670  +  .  Parent=EB634440_mid1;Target=EB634440 1 398
  scaffold_14830  BLAT    match_part   2489923 2490194 670  +  .  Parent=EB634440_mid1;Target=EB634440 412 683
  scaffold_14830  BLAT    EST_match    2484806 2485488 683  -  .  ID=EB634440_mid2;Target=EB634440 1 683
  scaffold_14830  BLAT    match_part   2484806 2485488 683  -  .  Parent=EB634440_mid2;Target=EB634440 71 753
  scaffold_14830  BLAT    EST_match    2480512 2481194 646  -  .  ID=EB634440_mid3;Target=EB634440 1 683
  scaffold_14830  BLAT    match_part   2480512 2480783 646  -  .  Parent=EB634440_mid3;Target=EB634440 71 342
  scaffold_14830  BLAT    match_part   2480797 2480867 646  -  .  Parent=EB634440_mid3;Target=EB634440 356 426
  scaffold_14830  BLAT    match_part   2480890 2481194 646  -  .  Parent=EB634440_mid3;Target=EB634440 449 753

=item equivalent Bioperl  bp_search2gff3.pl

  grep EB634440.1 dgri-gnoest.blat | perl -pe's/gi.\d+.(gb.\w+).\d+.//;' | \
  $in/lib/Bio/script/bp_search2gff3.pl  -f psl -m -ver 3 -t hit -i -

  ##gff-version 3
  scaffold_14830  BLAT    match_part   2489512 2490194 .  +   0   Parent=EB634440;Target=Sequence:EB634440 1 683
  scaffold_14830  BLAT    match_part   2484806 2485488 .  -   0   Parent=EB634440;Target=Sequence:EB634440 1 683
  scaffold_14830  BLAT    match_part   2480512 2481194 .  -   0   Parent=EB634440;Target=Sequence:EB634440 1 683
  scaffold_14830  BLAT    match        2480512 2490194 .  -   .   ID=EB634440



  cat dgri-gnoest.blat | sort -k14,14 -k10,10 -k16,16n | more
  
=cut
