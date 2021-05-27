#!/bin/bash
##  qsub -q batch velsubt.sh
#PBS -N velsubt
#PBS -l mem=127gb,nodes=1:ppn=8,walltime=23:55:00
#PBS -m abe
#PBS -j oe
#PBS -k o
#PBS -M gilbertd@indiana.edu

## note velv OPENMP threads trades w/ memory, more mem w/ more cpu; also uses +1 threads
ncpu=7
export OMP_NUM_THREADS=$ncpu
export OMP_THREAD_LIMIT=$ncpu

workd=$HOME/scratch/chrs/nasv1

cd $workd/rnas/velv/

echo "#.. start subsetvel : `date`"

export subset=nvit1_trsub1 ; $workd/rnas/velv-gg-subset.sh > log.$subset 2>&1 &

wait
echo "#.. end subsetvel : `date`"

