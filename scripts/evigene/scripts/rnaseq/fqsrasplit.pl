#!/usr/bin/env perl
# fqsrasplit.pl fastq , split read pairs /1 /2 in same file into _1.fq, _2.fq

foreach $fq ( @ARGV ) {

($fn=$fq) =~ s/.gz//; $fn=~s/\.(fq|fastq)//;
$fq1=$fn."_1.fq";
$fq2=$fn."_2.fq";
my($ok,$nerr,$npair)= (0) x 9;

if($fq and -f $fq) {
  next if( -f $fq1 and -f $fq2);
  warn "# split $fq to $fq1,$fq2\n";
  if($fq=~/.gz/) { $ok=open(F,"gunzip -c $fq |"); } else { $ok=open(F,$fq); }
  open(L,">$fq1") or die "write $fq1"; 
  open(R,">$fq2") or die "write $fq2";
  while(<F>) { 
    $fh=$_; $fs=<F>; $qh=<F>; $qs=<F>; ($lh)=$fh=~/^.(\S+)/;  
    $fhr=<F>; $fsr=<F>; $qhr=<F>; $qsr=<F>; ($rh)=$fhr=~/^.(\S+)/;  
    if($lh ne $rh) { $lh=~s,/[12],,; $rh=~s,/[12],,; }
    if($lh ne $rh) { $nerr++; warn "# mismatch pair: $lh ne $rh\n";  
        die "ERR: too many mismatch:$nerr" if($nerr>99); 
    } else {
      print L $fh,$fs,$qh,$qs;
      print R $fhr,$fsr,$qhr,$qsr;  $npair++;
    }
    } 
  close(L); close(R); close(F);  
  warn "# done: npair=$npair, nerr=$nerr for $fq to $fq1,$fq2\n";  
} else {
  warn "usage: fqsrasplit.pl readpairs.fastq[.gz]\n : split read pairs /1 /2 in same file into _1.fq, _2.fq\n";
  exit(1);
}

}
