#!/bin/bash
# env fastq=SSRnnnnn fastq2faq.sh

if [ ! -f $fastq ]; then  echo "usage: env fastq=fq.gz fastq2faq.sh"; die; fi
  
pt=`echo $fastq | sed 's/.fastq.gz//'`

gunzip -c $fastq | perl -ne\
'$k++; if($k % 4 == 3 and /^\+\w/) { $p=1; s/^./>/; print; } elsif($p) { $p=0; chomp; @q=split""; 
foreach $i (0..$#q) { $n=ord($q[$i]) - 32; print " $n"; print "\n" if($i % 20 == 19); } print "\n";}'\
 > $pt.fa.qual

gunzip -c $fastq | perl -ne\
'$k++; if($k % 4 == 1 and /^\@\w/) { $p=1; s/^./>/; print; } elsif($p) { $p=0; chomp; s/(.{60})/$1\n/g;  print $_,"\n";}'\
 > $pt.fa
