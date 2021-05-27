#!/usr/bin/env perl
# fasample.pl : sample subset of _1.fa, _2.fa paired rnaseq
# fasample.pl -parts 5 -nnnskip input_1.fa[.gz] (assume _2.fa or ignore)
# .. want variants for .fastq (1 file, 2 pair files) and .fasta2 pair file

use strict;
use Getopt::Long;

my $parts=5;
my $skipnnn=1;
my $debug=1;

my $optok= GetOptions(
  "parts=i", \$parts, 
  "skipnnn!", \$skipnnn, 
  "debug!", \$debug, 
  );
unless($optok and @ARGV) {
  die "usage: fasample.pl reads_1.fa[.gz]\n : create subset samples of input reads.fasta, paired if _1.fa and _2.fa \n";
}

foreach my $fain ( @ARGV ) {
  unless($fain and -f $fain) { warn "# missing $fain\n"; next; }
  
  my $paired=0; my $fa2in=""; 
  if($fain =~ /_1\.|\W1\.fa/) { 
    $fa2in= $fain;  $fa2in=~ s/_1\./_2./; $fa2in=~ s/1\.fa/2.fa/;
    $paired= (-f $fa2in)?1:0;
  }
  
  (my $fname=$fain) =~ s/.gz//; $fname=~s/\.(fa|fasta|fq|fastq)//;
  my($ok,$nin,$nskip,$nerr,$npair)= (0) x 9;
  
  if($fain=~/.gz/) { $ok=open(F,"gunzip -c $fain |"); } else { $ok=open(F,$fain); }
  if($paired) { if($fa2in=~/.gz/) { $ok=open(F2,"gunzip -c $fa2in |"); } else { $ok=open(F2,$fa2in); } }

  my %outh;
  for(my $i=0; $i<$parts; $i++) { 
    my $outna= "$fname.fa.pt$i";
    my $out2= $outna; $out2=~ s/_1\./_2./; $out2=~ s/1\.fa/2.fa/;
    my $fh=undef; open( $fh, ">$outna") or die "write $outna"; $outh{$i}{L}= $fh;
    if($paired) { my $fhr=undef; open( $fhr, ">$out2") or die "write $out2";  $outh{$i}{R}= $fhr; }
   }
   
  while(<F>) { 
    my($fh,$fs,$lh,$fhr,$fsr,$rh)=("") x 10;
    $fh=$_; $fs=<F>;  ($lh)=$fh=~/^.(\S+)/;  # $qh=<F>; $qs=<F>;
    if($paired) {
      $fhr=<F2>; $fsr=<F2>; ($rh)=$fhr=~/^.(\S+)/;  #  $qhr=<F>; $qsr=<F>;  
      if($lh ne $rh) { $lh=~s,/[12],,; $rh=~s,/[12],,; }
      if($lh ne $rh) { $nerr++; 
        warn "# mismatch pair: $lh ne $rh\n";  
        die "ERR: too many mismatch:$nerr" if($nerr>99); 
        next;
      } 
    }
    
    if($skipnnn and ($fs =~ /N/ or $fsr =~ /N/)) { $nskip++; next; }
    $nin++;
    my $ipart= $nin % $parts;
    my $outl = $outh{$ipart}{L};
    print $outl $fh,$fs;  # ,$qh,$qs
    if($paired) { my $outr= $outh{$ipart}{R}; print $outr $fhr,$fsr; } # ,$qhr,$qsr
    $npair++;
    } 

  close(F); close(F2) if($paired);
  for(my $i=0; $i<$parts; $i++) { close($outh{$i}{L}); close($outh{$i}{R}) if $paired; }
  warn "# done: nin=$nin; npair=$npair, nskip=$nskip, nerr=$nerr for $fain to $parts parts\n";  

}

__END__

./fasample.pl fastq/peshorti/SRR058495_1.fa
# done: nin=9023643; npair=9023643, nskip=23348, nerr=0 for fastq/peshorti/SRR058495_1.fa to 5 parts

921M Jun  9 19:44 SRR058495_1.fa
184M Jun 10 13:13 SRR058495_1.fa.pt0
184M Jun 10 13:13 SRR058495_1.fa.pt1
184M Jun 10 13:13 SRR058495_1.fa.pt2
184M Jun 10 13:13 SRR058495_1.fa.pt3
184M Jun 10 13:13 SRR058495_1.fa.pt4
526M Jun  7 16:56 SRR058495_1.fq.gz
921M Jun  9 19:45 SRR058495_2.fa
184M Jun 10 13:13 SRR058495_2.fa.pt0
184M Jun 10 13:13 SRR058495_2.fa.pt1
184M Jun 10 13:13 SRR058495_2.fa.pt2
184M Jun 10 13:13 SRR058495_2.fa.pt3
184M Jun 10 13:13 SRR058495_2.fa.pt4
