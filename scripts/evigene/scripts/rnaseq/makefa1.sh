#!/bin/tcsh
# for SRR data

##foreach faqz (*.fastq.gz)
foreach faqz (SRR063707.fastq.gz )
  set nam=`basename $faqz .fastq.gz`
  set fa="aphidrs_$nam.fa"  
  echo $fa
  # or? pick by line count: print 0,1 skip 2,3; ..
  gzcat $faqz | perl -ne 'if(/^\@(SRR\S+)/){print ">$1\n";$s=1;} elsif($s==1){$s=0; print;}' \
   > $fa
end

### fastq
## @SRR063706.4506 GA-B_0010:1:100:10342:15278 length=36
## GCGGTCGTCGTCCGGGTGGAACCCCGCACCGGAACC
## +SRR063706.4506 GA-B_0010:1:100:10342:15278 length=36
## GFGFGGGGGGGGGGEFAFFAFGGGG=GGDDDAGGGD

