#!/usr/bin/env perl
# mpile2tiletab.pl

# usage:
# samtools mpileup bams/Manak_20101102_lib[1-8].bam | perl mpile2tiletab.pl > readtiles.mk_metal.tab


use strict;

my $debug=$ENV{debug} || 0;

my $tloc="/bio/bio-grid/dpulex/tilex//dpxtiles0907/locs/daphnia_tilelocs.tab.gz";
open(LOC,"gunzip -c $tloc |") or die "$tloc";

my %tileloc; my $nloc=0;
while(<LOC>){ 
  my($r,$b,$e)=split; ## $r="scaffold_".$r; 
  push(@{$tileloc{$r}}, [$b,$e]); 
  $nloc++; }  
close(LOC);  print STDERR "# read nloc=$nloc\n";
# gzcat $bg/dpulex/tilex//dpxtiles0907/locs/daphnia_tilelocs.tab.gz | head
# scaf    start   end
# 10000   1       50
# 10000   31      80


my %nmap = (lib1 => 38916668, lib2 => 45681714, lib3 => 50565265,
  lib4 => 47340781, lib5 => 42540573, lib6 => 36261629, lib7 => 44106626, 
  lib8 => 45946087);

# make this relative to max val, so fac is always >=1 ?
my ($maxn)= sort{ $b <=> $a } values %nmap;
# my @libfac = map{ 1000000/$nmap{"lib$_"} } (1..8);
my @libfac = map{ $maxn/$nmap{"lib$_"} } (1..8);

my @libgrp= qw( cCon cCha  mCon mCd mAs mCu mNi mZn );
my @icount= map{ $_ * 3 }(0..7);

# input is loc-sorted as is tileloc
my ($btile, $etile, $otile, $rtile, $nin, $nout,$nskip)= (0) x 10;
my (@tiles,@ntiles,@itile);

putheader();
while(<>) {   
  my($r,$b,$cref,@v)=split"\t";  
  my @rd= @v[@icount];  $nin++;
  $r =~ s/scaffold_//;

  if($debug) {
  print STDERR "." if($nin % 1000 == 1);
  print STDERR "# $nin / $nskip / $nout\n" if($nin % 100000 == 0);
  }
  
  ## fixme overlapping tiles need dup counts
  ## > @itile = 1 or 2 values to add into
  
  if($r ne $rtile) {
    # puttiles( \@matiles) if(@matiles);
    sumtiles( $rtile, $btile, @itile);
    @tiles=();
    ($rtile, $btile, $etile, $otile, @itile)=  tileat($r,$b,0);
  } elsif( $b > $etile) {
    sumtiles( $rtile, $btile, $itile[0]);
    ($rtile, $btile, $etile, $otile, @itile)=  tileat($r,$b,$itile[0]);
  }
  
  if( $b < $btile) {  $nskip++; next; } # skip, no tile ??
  
  for my $it (@itile) {
    $ntiles[$it] ++; 
    for my $j (0..7) { $tiles[$it][$j] +=  $rd[$j]; } # should this be max not sum?
  }
}

sumtiles( $rtile, $btile, @itile);
print STDERR "#done: $nin / $nskip / $nout\n";

#............................

sub tileat {
  my($r, $b, $i)=@_;
  my $rloc= $tileloc{$r} 
    or return; ## do{ warn "tileat($r,$b,$i) tileloc $r bad \n"; return; }

  $i ||= 0;
  for( ; $i>=0 ; $i++) {
    # die "tileat($r,$b,$i) rloc $i bad \n" 
    return undef unless($rloc->[$i]);
    my($tb,$te)= @{$rloc->[$i]};
    if($b < $tb) { 
      return ($r,$tb,$te,0,$i); # will skip read
      } # die "tileat($r,$b,$i) $b < $tb\n"; return undef; } # error
    if($b >= $tb and $b <= $te) { 
      # print STDERR "# tileat($r,$b,$i) = $tb,$te\n" if($nin < 60);
      if($rloc->[$i+1]) {
      my($ob,$oe)= @{$rloc->[$i+1]}; 
      if($b >= $ob and $b <= $oe) { return ($r,$tb,$te,$oe,$i,$i+1); } # overlap 2 tiles
      }
      return ($r,$tb,$te,0,$i); 
      }
  }
}

sub putheader {
  my @rexp= map{ ($_."A",$_."M") } @libgrp[1,3..7];  
  print join("\t","loc",@rexp)."\n";
}

sub sumtiles {
  my($r,$b,@itile)= @_;
  # print STDERR "# sumtiles($r,$b) \n" if($nin < 600);
  return unless($r and @itile);
  
  foreach my $it (@itile) {
  # my @rd= @{$tiles[$it][0..7]};
  # my @rd = $tiles[$it][0..7];

  my $nc = $ntiles[$it] || 1;
  my @rd;
  for my $j (0..7) {
    my $rd= $tiles[$it][$j];
    my $rf= ( $rd/$nc ) * $libfac[$j];
    $rd[$j]= ($rf < 1) ? 0 : log($rf); 
    }
   
  my @rcon= @rd[0,2];
  my @rexp= @rd[1,3..7];

  ## do this after sum to tiles, add @libfac before log()
  # my @MA = ("$r:$b"); # tile loc
  # $r =~ s/scaffold_//;
  print "$r:$b";
  for my $j (0..$#rexp) { 
    my $jc=($j>0)? 1: 0; 
    # my $jma= 1 + $j*2;
    my $aval= 0.5*($rexp[$j] + $rcon[$jc]); 
    my $mval= $rexp[$j] - $rcon[$jc];
    # push(@MA, $aval, $mval);
    printf "\t%.4f\t%.4f",$aval,$mval;
    }
    
  print "\n";  $nout++;
  # $matiles[$it]= \@MA;  
  }
  
}

