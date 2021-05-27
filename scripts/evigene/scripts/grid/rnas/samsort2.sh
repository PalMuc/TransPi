#! /bin/bash -l
### env rnain=xxxxx qsub -q batch samsort.sh
#PBS -N samsort
#PBS -A ind114
#PBS -l nodes=1:ppn=8,walltime=03:55:00
#PBS -o sams1.$$.out
#PBS -e sams1.$$.err
#PBS -V

workd=$HOME/scratch/chrs/aphid2
bindir=$HOME/bio/bin
dgenome=aphid2asm

cd $workd/rnas/parts

# merge merged parts
samset=SRR071347
# samset=$rnain
# samset=`basename $samset .gsnap1.samu`

$bindir/samtools merge $samset.bam ${samset}a?.bam

## redo this w/ samstats.pl, this not useful
$bindir/samtools flagstat $samset.bam > $samset.stats

