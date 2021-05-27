#!/bin/bash
# makeparts.sh

workd=/export/udisk3/work/aphid/
scripts=$workd/scripts/
bamset="bams/*.bam aphid_est.sam"
partlist=sparts.list

for bam in $bamset; {
  case "$bam" in 
  *aphidpe_*) 
  $scripts/sam2seqparts.pl -nodupl  -format "fasta,sam" -type bam_pair -in $bam -part $partlist -debug
  $scripts/sam2seqparts.pl -nodupl  -type splitpairs -in $bam -part $partlist -debug  
  ;;  
  *est.sam)
  $scripts/sam2seqparts.pl -nodupl -format "fasta" -type sam_est -in $bam -part $partlist -debug 
  ;;
  *)
  $scripts/sam2seqparts.pl -nodupl -format "fasta,sam" -type bam_single -in $bam -part $partlist -debug 
  ;; 
  esac
}


