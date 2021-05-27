#! /bin/bash
## env trset=mygenes.tr datad=`pwd` prog=./runtr2cds.sh sbatch srun_prog.sh
## comet.sdsc sbatch update; show_accounts ?? not ind114
## .. redo as 2 script set as aprun/bigred hack, call w/ 2nd script param
#SBATCH -A ind114
#SBATCH --job-name="tr2cds"
#SBATCH --output="tr2cds.%j.%N.out"
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH -t 39:55:00
#SBATCH --export=ALL

## 2/20cpu blastn failed evg2, try ncpu=16
export ncpuall=24
#x#export ncpu=20
## BIG DATA atrsetevg2all.tr blastn fail 20,12 cpu cut to 8 cpu?
if [ "X" = "X$ncpu" ]; then export ncpu=8; fi
export maxmem=120000

## FIXME: cd-hit at least is OMP threaded app, add the ENV for it, use mrun_comet.sh??
export OMP_NUM_THREADS=$ncpu

if [ "X" = "X$datad" ]; then echo "ERR:datad=what?"; exit -1; fi
if [ "X" = "X$prog" ]; then echo "ERR:prog=what?"; exit -1; fi
cd $datad
if [ "X" != "X$name" ]; then echo "#runname=$name datad=$datad"; fi

## * this is wacky like aprun/ bigred in that it calls $prog ncpu times !! need opt to not do that
echo ibrun --npernode 1 -v $prog
ibrun --npernode 1 -v $prog

#---------
