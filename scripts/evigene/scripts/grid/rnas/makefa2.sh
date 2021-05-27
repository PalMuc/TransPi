#!/bin/tcsh
# fastq , merge /1 /2 into gsnap.fa2 format

# setenv TMPDIR /var/tmp
setenv TMPDIR /scratch/dgilbert

## fixme: no /1 /2 on SRA ids

foreach fq (fastq/*_1.fastq.gz )
  set nam=`basename $fq _1.fastq.gz`
  set fa2="aphidpe_$nam.fa2"  
  if( -f $fa2 ) continue
  echo -n "# ${nam}_[12].fastq to $fa2 "

## out a mem on did{$d} ... split tasks
gunzip -c fastq/${nam}_1.fastq.gz | perl -ne \
'BEGIN{$t=1;} if(/^\@(SR\S+)/){ $d=$1; print ">$d/$t "; $s=1; } elsif($s) { $s=0; print;}' \
> $TMPDIR/$fa2.1.tmp
 
gunzip -c fastq/${nam}_2.fastq.gz | perl -ne \
'BEGIN{$t=2;} if(/^\@(SR\S+)/){ $d=$1; print ">$d/$t "; $s=1; } elsif($s) { $s=0; print;}' \
> $TMPDIR/$fa2.2.tmp

   sort $TMPDIR/$fa2.[12].tmp | \
   perl -pe 'if(s,/1,,){ s/ /\n/; $n1++; }elsif(s,/2,,){s/^>\S+\s*//; $n2++;} \
   END{ warn"# n1=$n1 n2=$n2\n";}' > $fa2

   /bin/rm $TMPDIR/$fa2.[12].tmp

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
