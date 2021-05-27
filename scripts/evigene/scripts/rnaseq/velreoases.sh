#!/bin/bash
##  qsub -q batch velreoases.sh
#PBS -N velreoas
#PBS -l mem=127gb,nodes=1:ppn=3,walltime=35:55:00
#PBS -m abe
#PBS -j oe
#PBS -k o
#PBS -M gilbertd@indiana.edu

ncpu=2
export OMP_NUM_THREADS=$ncpu
export OMP_THREAD_LIMIT=$ncpu

workd=$HOME/scratch/chrs/nasv1
cd $workd/rnas/velv/

## rerun oases out-of-time at 24hr
export subset=nvit1_trsub1 ;
velbin=$HOME/bio/velvet/bin

echo "#.. start reoases $subset : `date`"
if [ 1 == 1 ]; then
  $velbin/oases vel$subset -ins_length 200 -ins_length_sd 50 -min_trans_lgth 40 
  /bin/rm vel$subset/{Graph2,LastGraph,PreGraph,Roadmaps,Sequences}
fi

if [ -f vel$subset/transcripts.fa ]; then
  echo "#.. OK "
else
  echo "#.. ERROR: no vel/transcripts.fa"
fi
echo "#.. end reoases : `date`"

