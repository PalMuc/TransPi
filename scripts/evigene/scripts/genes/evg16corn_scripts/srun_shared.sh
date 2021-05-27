#! /bin/bash
## srun_redo.sh : shared que rerun a few jobs
#SBATCH -A ind114
#SBATCH --job-name="tr2cds"
#SBATCH --output="tr2cds.%j.%N.out"
#... SBATCH --partition=compute
#... SBATCH --nodes=1
#... SBATCH --ntasks-per-node=24
#... SBATCH -t 39:55:00
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=30G
#SBATCH -t 04:00:00
#SBATCH --export=ALL

export ncpu=8
export maxmem=20000

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
