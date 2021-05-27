#!/usr/bin/env perl
# genomewigdens.pl

=item about

genome feature densities from bigwig.bw data,
using kent/bin/bigWigSummary
density count per sliding base window

  env bw=bigwig.bw bin=50000 slide=1 feat=rnanc1 stat=mean chrs=genome/cacao11chrs.fa.count genomewigdens.pl
  kent/bin/bigWigSummary -type=mean data.bw  chr start end
  
=cut

use strict;

my $usage=
"\nenv bw=bigwig.bw bwpath=/path/to/bigWigSummary bin=50000 slide=1 feat=rnanc1 stat=mean chrs=genome/cacao11chrs.fa.count genomewigdens.pl\n";

my $bigWigSummary=  $ENV{bwpath} ||"bigWigSummary"; # in path?
my $ok=`which $bigWigSummary`; die "Need bwpath= ; $usage" unless($ok =~ m,/,);

my $bwfile= $ENV{bw} or die "Need bw=bigwigfile; $usage";
my $chrs  = $ENV{chrs} or die "Need chrs=chrsizesfile: name\tsize; $usage";
my $bins  = $ENV{bin} || 50000;
my $slide = $ENV{slide} || 0; 
my $stat  = $ENV{stat} || "mean";
my $fname = $ENV{feat} || "reads";
my $debug = $ENV{debug} || 0;
#my $typesource= defined $ENV{src} ? $ENV{src} : 0;

#my $bin2 = 2*$bins;
my $binstep= ($slide>0) ? int($bins/(1+$slide)) : $bins;

my (@chrs,%chrs);
open(C,$chrs) or die "$chrs";
while(<C>){ next unless(/^\w/); my($c,$clen)=split; push @chrs, $c; $chrs{$c}= $clen; } close(C);

print join("\t","scaf","loc",$fname),"\n";

foreach my $c (@chrs) {
  my $clen= $chrs{$c};
  my @stepbounds=(); # my($cstart, $cend)=(0,0);
  for(my $k=1; $k<$clen; ) {
    my($cb,$ce)= ($k, $k+$bins); $ce=$clen if($ce>=$clen);
    push(@stepbounds, $cb-1,$ce-1); # use 0-origin; need only start?
    $k= $k + $binstep;
  }
  
  my $parts= int(scalar(@stepbounds) / 2);
  my @vals= bigval( $c, 0, $clen-1, $parts );
  my $cnum=$c; $cnum =~ s/(contig|chr|scaffold|super)\D*//i; # numeric part only?

  for(my $i=0; $i<@vals; $i++) {
    my $val= $vals[$i];  $val=0 unless($val =~ /\d/); # n/a
    my $j= 2*$i;
    my($cb,$ce)= @stepbounds[$j,$j+1];
    print join("\t",$cnum,$cb,$val),"\n";  # use 0-origin for .bed
  }
}

sub bigval {
  my($chr,$cb,$ce,$parts)= @_;
  my $cmd="$bigWigSummary -type=$stat $bwfile $chr $cb $ce $parts"; # 1 or many steps?
  my $res=`$cmd`; # errors?
  print STDERR "# $bigWigSummary $bwfile $chr $cb $ce $parts = [$res]\n" if($debug);
  my @vals= split " ",$res;
  return @vals;
}


