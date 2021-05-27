#! /bin/bash -l
### qsub -q batch bamcuff.sh
#PBS -N bamcuff
#PBS -A ind114
#PBS -l nodes=1:ppn=8,walltime=33:55:00
#PBS -o bamcuf1.$$.out
#PBS -e bamcuf1.$$.err
#PBS -V

workd=$HOME/scratch/chrs/aphid2
bin_dir=$HOME/bio/bin

drnabam=$workd/rnas/bams/merge_all6.bam
if [ "${rnain}x" != "x" ]; then drnabam=$rnain; fi
drna=`basename $drnabam .bam`

## try higher premrna to reduce spurious joins?
opts="-p 8 -v --pre-mrna-fraction 0.25 "

echo "start bamcuff: `date`"
cd $workd/rnas

mkdir cuff.$drna
cd cuff.$drna
echo $bin_dir/cufflinks09 $opts $drnabam 
$bin_dir/cufflinks09 $opts $drnabam  > log.cuf  2>&1

echo "end bamcuff: `date`"

