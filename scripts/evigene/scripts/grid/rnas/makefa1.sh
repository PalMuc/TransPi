#!/bin/tcsh
# for SRR data

foreach faqz (fastq/*.fastq.gz)
  set nam=`basename $faqz .fastq.gz`
  if( $nam =~ *_[12] ) continue
  set fa="aphidrs_$nam.fa"  
  echo $fa
  gunzip -c $faqz | perl -ne 'if(/^\@(SRR\S+)/){print ">$1\n";$s=1;} elsif($s==1){$s=0; print;}' > $fa
end

### fastq
## @SRR063706.4506 GA-B_0010:1:100:10342:15278 length=36
## GCGGTCGTCGTCCGGGTGGAACCCCGCACCGGAACC
## +SRR063706.4506 GA-B_0010:1:100:10342:15278 length=36
## GFGFGGGGGGGGGGEFAFFAFGGGG=GGDDDAGGGD

