#!/bin/bash

bin_dir=$HOME/bio/tophat/bin
opts="-p 2 --inner-dist-mean 100"
samall=all.sam
partdir=$PWD

echo "#.. start samcuff: `date`"

if [ -d $partdir/cuff8 ]; then
  echo "#.. Already DONE $partdir/cuff8"
  exit
fi

## parts run under same PBS_JOBID
mypart=pt$$
cd /scratch/$PBS_JOBID
mkdir $mypart ; cd $mypart
## dont use sort -m cuz cufflinks whines if scaf order isnt to its like
export LC_ALL="C"
sort -T ./ -k 3,3 -k 4,4n  $partdir/*.sam > $samall

mkdir cuff8 ; cd cuff8
$bin_dir/cufflinks $opts ../$samall  > log.cuf$samall 2>&1

cd ../ ; cp -rp cuff8 $partdir/
# clear out, job is still filling up this scratch disk
/bin/rm -rf $samall cuff8 
echo "# .. end samcuff: `date`"

