#! /bin/bash -l
### qsub -q batch makepjob.sh
#PBS -N makepj
#PBS -A ind114
#   #debug# -l nodes=1:ppn=8,walltime=23:55:00
#PBS -l nodes=1:ppn=2
#PBS -o makepj.$$.out
#PBS -e makepj.$$.err
#PBS -V

# makeparts.sh

#workd=/export/udisk3/work/aphid/
workd=$HOME/scratch/chrs/aphid2
scripts=$workd/scripts/
export PATH=$HOME/bio/bin:$PATH

cd $workd/rnas/
## split bamset to 4 parts for 8-core 1 node job
# bamset="bams/aphidrs_*_1.bam bams/bams/aphidrs_*W.bam"
# bamset="bams/aphidrs_SRR*.bam  aphid_est.sam"
# bamset="bams/aphidpe_SRR07580?*.bam"
## test case
bamset="bams/aphidrs_SRR073576.bam"

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


