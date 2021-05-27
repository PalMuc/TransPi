#! /bin/bash -l
### qsub -q batch makepjob.sh
#PBS -N makepj
#PBS -A ind114
#PBS -l nodes=1:ppn=8,walltime=11:55:00
#PBS -o makepj.$$.out
#PBS -e makepj.$$.err
#PBS -V

workd=$HOME/scratch/chrs/aphid2
bindir=$HOME/bio/bin

cd $workd/rnas/

./makepartr.sh > log.part1r 2>&1 &
# ./makepart2.sh > log.part2 2>&1 &
./makepart3.sh > log.part3 2>&1 &
./makepart4.sh > log.part4 2>&1 &

wait

