#!/usr/bin/env perl
# bam2bw.pl

use strict;

# change to getOptions() ..
my $debug=$ENV{debug}||0;
my $dolog= (defined $ENV{log}) ? $ENV{log} : 1;
my $doave= (defined $ENV{ave}) ? $ENV{ave} : 1;  # of n-bams
my $roundfac= $ENV{round} || 0;
my $normfac= (defined $ENV{norm}) ? $ENV{norm} : 0;  # of n-bams
my $samtools=$ENV{samtools} || "samtools";
my $sopt= $ENV{opt} || "-d 49999";
my $allstats= $ENV{all} || 0;

my $dgenomesize= shift @ARGV;
#? add genome.faidx to get genome base vs read bases?
my @inbam= grep /\.bam/,  @ARGV;
my $inbam= join " ", @inbam;
my $nbam=@inbam;

### this is true, but the seq col has <<>> for introns (N) and other info to count..
# warn "FIXME: new samtools mpileup DOES NOT count reads/base; it is counting read span? i.e. NO INTRONS\n
# need older samtools pileup one.bam, or other tool to turn reads.bam to readcover.bed counts";

## mpileup opt -d maxreads-per-base-bam may be needed; default autosets (depend on bam count?) 
##    [mpileup] 7 samples in 7 input files : Set max per-file depth to 1142 
if($dgenomesize =~ m/\.fai$/ and $sopt !~ m/\-f/) { (my $dg=$dgenomesize) =~ s/.fai$//;  $sopt .= " -f $dg"; }

warn "bam2bw opts: log=$dolog norm=$normfac ave=$doave opt=$sopt all=$allstats\n
    samtools=$samtools genosize=$dgenomesize nbam=$nbam $inbam[0]..\n" if($debug);
die "usage: bam2bw.pl /path/to/genome.chr_size.tab in1.bam [in2.bam ..]"
  unless(-f $dgenomesize and -f $ARGV[0]);

open(S,$dgenomesize) or die; 
my %rsize=(); while(<S>){ my($r,$n)=split; $rsize{$r}=$n; } close(S);

if($normfac>0) { $doave=0; } elsif($doave and $nbam>1) { $normfac= 1/$nbam; }
## my $makeave= ($doave and $nbam>1)?1:0;
my $hasref= ($sopt =~ m/-f \S/) ? 1 : 0;

open(IN, "$samtools mpileup $sopt $inbam |") or die "samtools mpileup $sopt $inbam";

  # separate loop by options so it is speedy
if($allstats) { allstats(); }  # readcov, introncov, bothcov
else { readcover(); }

close(IN);

sub readcover
{
  my($lr,$lb,$le,$lc,$li,$lm);
  my $halfr= ($roundfac > 0) ? int($roundfac/2) : 0;
  while(<IN>) {
    my($r,$b,$xn,@v)=split"\t";
    my $c=0;
    if($hasref) {
      for(my $i=1; $i<@v; $i+=3) { $c +=  $v[$i] =~ tr/.,/.,/;  }
    } else {
      for(my $i=1; $i<@v; $i+=3) { $c +=  $v[$i] =~ tr/ACGTacgt/ACGTacgt/;  }
    }
    if($dolog) { 
      $c *= $normfac if($normfac>0); #gone  elsif($makeave) { $c /= $nbam; }
      $c= int(10 * log(1+$c))/10; 
      } 
    elsif($normfac>0) { $c = int($c * $normfac); } #gone elsif($makeave) { $c = int($c/$nbam); }
## add round factor: c = $fac * int( $c/$fac) ; fac = 10 or 5 or whatever
    elsif($roundfac>0) { $c = $roundfac * int( ($halfr + $c) / $roundfac); }
    if($r eq $lr and $c == $lc and $le == $b-1 ) { $le=$b; }
    else { putb1($lr,$lb-1,$le,$lc);  $le=$lb=$b; }
    ($lr,$lc)= ($r,$c);
  }
  putb1($lr,$lb-1,$le,$lc); # last
}

sub allstats 
{
  my($lr,$lb,$le,$lc,$li,$lm);
  while(<IN>) {
  my($r,$b,$xn,@v)=split"\t"; 
  my ($c,$in,$mc)=(0,0,0); 

  for(my $i=1; $i<@v; $i+=3) { # i=1 offset to seq
    $mc += $v[$i-1];  # mpileup: triplets per in.bam file: count,seq,qual
    my $s= $v[$i];
    my($nc,$ni);
    if($hasref) { $nc = $s=~tr/.,/.,/; }
    else { $nc= $s=~tr/ACGTacgt/ACGTacgt/; }  # OR count "." if using genome.faidx
    $ni= $s=~tr/<>/<>/;
    $c += $nc;
    $in += $ni;
    }

  if($dolog) {
    map{
     $_ *= $normfac if($normfac>0); # $_ /= $nbam if($makeave); 
     $_= int(10 * log(1+$_))/10; } ($c,$in,$mc);
  } elsif($normfac>0) {  
    map{ $_ = int($_ * $normfac); } ($c,$in,$mc); 
  }

  if($r eq $lr and $c == $lc and $in == $li and $le == $b-1 ) { $le=$b; } 
  else { putb2($lr,$lb-1,$le,$lc,$li,$lm);  $le=$lb=$b; }  
  ($lr,$lc,$li,$lm)= ($r,$c,$in,$mc); 
  } 

  putb2($lr,$lb-1,$le,$lc,$li,$lm); # last
}

sub putb1 { 
  if($_[0] && $_[2] > $_[1] && $_[3] >= 0.5 ) { 
  my $s= $rsize{$_[0]}||0; return if($_[2] > $s or $_[1] >= $s); 
  print join("\t",@_),"\n"; }  
}

sub putb2 {
  if($_[0] && $_[2] > $_[1] && ($_[3] >= 0.5 || $_[4] >= 0.5)) {
  my $s= $rsize{$_[0]}||0; return if($_[2] > $s or $_[1] >= $s);
  print join("\t",@_),"\n"; }
}

__END__

# merge subset.bam to one bed/bw set
for i in 3 4 16; do   { 
 for j in 1 8 C; do {
  brna="merge-$i-$j" 
  inbam=`ls $i-${j}R*-dmag2.bam`
  echo bam2bw $brna FROM $inbam
  
  bam2bw.pl $dgenomesize $inbam > $brna.bed
  ## nasty fails when read beyond end if dgenomesize.. filter
  bedGraphToBigWig $brna.bed $dgenomesize  $brna.bw
  
 } done
} done  
