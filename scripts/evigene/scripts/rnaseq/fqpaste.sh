#!/bin/bash
# fqpaste.sh : paste 2 cut parts of rnase fastq to one full read set

fq1=004_R1_1_40.fastq.gz
fq2=004_R1_41_80.fastq.gz
fout=004_R1.fastq
TMPDIR=/var/tmp

gunzip -c $fq1 $fq2 | perl -ne 's/\n/\t/; 
if(s/^\@//){ s/_([12])(\-\S+)/$2\t$1/; print "\n",$_; }
elsif(/^\+/){ } else { print; }  END{ print"\n"; } ' \
| sort | perl -ne \
'($g,$p,$s,$q)=split; END{ puts(); }
if($p == 1) { puts(); ($lg,$ls,$lq)=($g,$s,$q); } 
elsif($p == 2) { if($g eq $lg) { $ls.=$s;  $lq.=$q; } 
else { puts(); ($lg,$ls,$lq)=($g,$s,$q); }  puts(); }  
sub puts { print join("\n", "\@".$lg, $ls, "\+".$lg, $lq),"\n" if($ls); ($lg,$ls,$lq)=("") x 4; }' \
> $fout


