#!/usr/bin/env perl
# fastq , merge /1 /2 into velvet.fa2 format
# input .fq1, fq2 already revcomp/1, sorted.

($fq1,$fq2)= @ARGV;
# ($nam=$fq1) =~ s/\.\w+$//;

if($fq2 and -f $fq2) {
open(F,$fq1); open(R,$fq2);
while(<F>) { 
  $f=$_;  $fs=<F>; ($fh)=$f=~/^.(\S+)/; print ">$fh\n$fs"; $fq=<F>; $fq=<F>; 
  $f=<R>; $fs=<R>; ($fh)=$f=~/^.(\S+)/; print ">$fh\n$fs"; $fq=<R>; $fq=<R>; 
  } 
} elsif($fq1 and -f $fq1) {
open(F,$fq1); 
while(<F>) {
  $f=$_; $fs=<F>; ($fh)=$f=~/^.(\S+)/; print ">$fh\n$fs"; $fq=<F>; $fq=<F>; 
  } 
} else {
  warn "usage: velfq2fa.pl seqs.fq1 seqs.fq2; paste fastq pairs to sequentially fasta for velvet\n";
  exit(1);
}


