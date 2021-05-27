#! /bin/bash
## env prog=velrun9.sh inpe=xxx name=xxx datad=`pwd` sbatch mrun_comet.sh
#SBATCH -A ind114
#SBATCH --job-name="runlsc"
#SBATCH --output="runlsc.%j.%N.out"
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH -t 39:55:00
#SBATCH --export=ALL

export ncpuall=24
# export ncpu=12
export ncpu=24
export mem=128G
export OMP_NUM_THREADS=$ncpu

if [ "X" = "X$datad" ]; then echo "ERR:datad=what?"; exit -1; fi
if [ "X" = "X$prog" ]; then echo "ERR:prog=what?"; exit -1; fi
cd $datad

if [ "X" != "X$name" ]; then echo "#runname=$name datad=$datad"; fi
echo ibrun --npernode 1 -v $prog
ibrun --npernode 1 -v $prog

