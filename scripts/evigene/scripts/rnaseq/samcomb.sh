#! /bin/bash -l
### qsub -q batch samcomb.sh
#PBS -N samcomb 
#PBS -A ind114
#PBS -l nodes=1:ppn=2,walltime=23:55:00
#PBS -o samc1.out
#PBS -e samc1.err
#PBS -V

workd=$HOME/scratch/chrs/daphmag
outna=rnaseqall

echo "start samcom : `date`"

cd /scratch/$PBS_JOBID

## 11 Gb total compressed, 16 files : sort merge to 1
cp -rp $workd/rnas/sams .

gunzip sams/*.sam.gz

#sort -T ./ -m -k 3,3 -k 4,4n  sams/*.sam > $outna.sam
sort -T /var/tmp -m -k 3,3 -k 4,4n  sams/*.sam > $outna.sam

gzip $outna.sam

cp /scratch/$PBS_JOBID/$outna.sam.gz $workd/rnas/

echo "end samcom: `date`"

