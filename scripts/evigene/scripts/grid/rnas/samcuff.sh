#! /bin/bash -l
### qsub -q batch samcuff.sh
#PBS -N samcuff
#PBS -A ind114
#PBS -l nodes=1:ppn=2,walltime=33:55:00
#PBS -o samcuf1.out
#PBS -e samcuf1.err
#PBS -V

workd=$HOME/scratch/chrs/daphmag
bin_dir=$HOME/bio/tophat/bin

samall=rnaseqall.sam
opts="--inner-dist-mean 100"

echo "start samcuff: `date`"
## data too big for /scratch; do in /work
# cd /scratch/$PBS_JOBID
cd $workd/rnas

gunzip sams/*.sam.gz
sort -T ./ -m -k 3,3 -k 4,4n  sams/*.sam > $samall

mkdir cuff.$drna
cd cuff.$drna
$bin_dir/cufflinks $opts ../$samall  > log.cuf$samall 2>&1

# cp -rp /scratch/$PBS_JOBID/cuff.$drna $workd/rnas/
echo "end samcuff: `date`"

