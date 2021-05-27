#! /bin/bash -l
### qsub -q normal rnavelv.sh
#PBS -N rnavel4
#PBS -A TG-MCB100147
## for -q normal ember.ncsa ncpus=6 6cpu/node mem=NNNmb
#PBS -l ncpus=6,mem=63gb,walltime=10:55:00
#PBS -o velv4.out
#PBS -e velv4.err
#PBS -V

# acyr_est201009.Scaffold2.fa  merged_pe.Scaffold2.fa2  merged_rs.Scaffold2.fa  rnavelv.sh*

#workd=$HOME/scratch/chrs/daphmag
workd=/export/udisk3/work/aphid
datad=$workd/rnas/velv
#bin_dir=$HOME/bio/velvet/bin
bin_dir=/bio/bio-grid/mb/rnaseq/velvet/bin

echo "start rnavelv : `date`"
cd $datad

$bin_dir/velveth vel4r 21 \
 -fasta -shortPaired merged_pe2.Scaffold2.fa2 \
 -short merged_rs.Scaffold2.fa merged_pe2.Scaffold2.fa1 \
 -long acyr_est201009.Scaffold2.fa \
> log.velr1 2>&1

$bin_dir/velvetg vel4r -read_trkg yes > log.velr2 2>&1

$bin_dir/oases vel4r -ins_length 200 > log.velr3 2>&1

echo "end rnavelv : `date`"

