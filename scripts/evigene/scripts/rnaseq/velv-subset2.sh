#! /bin/bash -l
##  qsub -q normal subset2.sh
#PBS -N velsub2
#PBS -A TG-MCB100147
#PBS -l ncpus=6,mem=127gb,walltime=23:55:00
#PBS -o velsub.$$.out
#PBS -e velsub.$$.err
#PBS -V

workd=$HOME/scratch/chrs/aphid2
velbin=$HOME/bio/velvet/bin

cd $workd/rnas/subset2/

echo "#.. start subsetvel : `date`"

export subset=subset.mid0 ; $workd/rnas/gg-subset.sh > log.$subset 2>&1
export subset=subset.mid1 ; $workd/rnas/gg-subset.sh > log.$subset 2>&1
export subset=subset.mid2 ; $workd/rnas/gg-subset.sh > log.$subset 2>&1
export subset=subset.mid3 ; $workd/rnas/gg-subset.sh > log.$subset 2>&1
export subset=subset.mid4 ; $workd/rnas/gg-subset.sh > log.$subset 2>&1
export subset=subset.mid5 ; $workd/rnas/gg-subset.sh > log.$subset 2>&1
export subset=subset.mid6 ; $workd/rnas/gg-subset.sh > log.$subset 2>&1
export subset=subset.mid7 ; $workd/rnas/gg-subset.sh > log.$subset 2>&1
export subset=subset.mid8 ; $workd/rnas/gg-subset.sh > log.$subset 2>&1
export subset=subset.mid9 ; $workd/rnas/gg-subset.sh > log.$subset 2>&1

echo "#.. end subsetvel : `date`"

