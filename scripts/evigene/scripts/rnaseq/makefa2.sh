#!/bin/tcsh
# fastq , merge /1 /2 into gsnap.fa2 format

foreach fq (*_1.txt.gz )
  set nam=`basename $fq _1.txt.gz`
  if( -f $nam.fa2.gz ) continue
  if( -f $nam.fa2 ) continue
  echo -n "# ${nam}_[12].fastq to gsnap.fa2 "

  gzcat ${nam}_[12].txt.gz | perl -ne \
  'if(/^\@(\w\S+)/){ $d=$1; $d=~s/[:]/x/g; print ">$d "; $s=1; } elsif($s) { $s=0; print;}' | \
   sort -T /var/tmp | perl -pe 'if(s,/1,,){ s/ /\n/; $n1++; }elsif(s,/2,,){s/^>\S+\s*//; $n2++;} \
   END{ warn"# n1=$n1 n2=$n2\n";}' > $nam.fa2

  echo
end

#   Paired-end RNAseq
# ------> 3-1_R2_1.txt.gz <------
# @HWI-EAS222_0025:3:1:1331:1963#0/1
# CGGCAAGCTTCAAGCAACTCANTATCACNNGACGTATTGATAGTGTCACGC
# +HWI-EAS222_0025:3:1:1331:1963#0/1
# ------> 3-1_R2_2.txt.gz <------
# @HWI-EAS222_0025:3:1:1331:1963#0/2
# TGGAAATANTCTGAGNNGCGTGACACTATCAATACGTCTTGTGTTATTGCG
# +HWI-EAS222_0025:3:1:1331:1963#0/2
