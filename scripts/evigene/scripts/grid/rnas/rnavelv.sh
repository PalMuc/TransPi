#! /bin/bash -l
### qsub -q large rnavelv.sh
#PBS -N rnavelv2
#PBS -A ind114
#... -l nodes=1:ppn=2,walltime=43:55:00
#PBS -l nodes=1:ppn=1,walltime=43:55:00
#PBS -o velv2.out
#PBS -e velv2.err
#PBS -V

workd=$HOME/scratch/chrs/daphmag
datad=$workd/rnas/velv
bin_dir=$HOME/bio/velvet/bin

echo "start rnavelv : `date`"
cd $datad
## this should be /scratch/$PBS_JOBID  ***
# cd /ramfs/$PBS_JOBID
# cd /scratch/$PBS_JOBID

$bin_dir/velveth vel2 21 \
 -fasta.gz -shortPaired $datad/mags.pairs.s01361.fa.gz \
 -fasta.gz -long $datad/daphnia_plxmag_est.sread.s01361.fa.gz \
 -fasta.gz -longPaired $datad/daphnia_plxmag_est.pairs.s01361.fa.gz \
> log.vel1 2>&1

$bin_dir/velvetg vel2 -read_trkg yes > log.vel2 2>&1

$bin_dir/oases vel2 -ins_length 200 > log.vel3 2>&1

#cp -r /ramfs/$PBS_JOBID $workd/rnas/
#cp -r /scratch/$PBS_JOBID $workd/rnas/
echo "end rnavelv : `date`"

