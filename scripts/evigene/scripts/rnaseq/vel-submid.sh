#!/bin/bash
##  qsub -q batch velsubt.sh
#PBS -N velmid
#PBS -l mem=256gb,nodes=1:ppn=32,walltime=23:55:00
#PBS -m abe
#PBS -j oe
#PBS -k o
#PBS -M gilbertd@indiana.edu

## note velv OPENMP threads trades w/ memory, more mem w/ more cpu; also uses +1 threads
## note2: oases doesnt use threads; not much gain, use 1 core/subset or 2? : 10 subs/part
## factor in memory total = 10x subset; how many cores/node??
ncpu=3
export OMP_NUM_THREADS=$ncpu
export OMP_THREAD_LIMIT=$ncpu

workd=$HOME/scratch/chrs/nasv1

cd $workd/rnas/velmid/

echo "#.. start subsetvel : `date`"

export subset=subset.mid0 ; $workd/rnas/velv-gg-subset.sh > log.$subset 2>&1 &
export subset=subset.mid1 ; $workd/rnas/velv-gg-subset.sh > log.$subset 2>&1 &
export subset=subset.mid2 ; $workd/rnas/velv-gg-subset.sh > log.$subset 2>&1 &
export subset=subset.mid3 ; $workd/rnas/velv-gg-subset.sh > log.$subset 2>&1 &
export subset=subset.mid4 ; $workd/rnas/velv-gg-subset.sh > log.$subset 2>&1 &
export subset=subset.mid5 ; $workd/rnas/velv-gg-subset.sh > log.$subset 2>&1 &
export subset=subset.mid6 ; $workd/rnas/velv-gg-subset.sh > log.$subset 2>&1 &
export subset=subset.mid7 ; $workd/rnas/velv-gg-subset.sh > log.$subset 2>&1 &
export subset=subset.mid8 ; $workd/rnas/velv-gg-subset.sh > log.$subset 2>&1 &
export subset=subset.mid9 ; $workd/rnas/velv-gg-subset.sh > log.$subset 2>&1 &

wait
echo "#.. end subsetvel : `date`"

