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

samset=$rnain
samset=`basename $samset .gsnap1.samu`

for i in  0 1 2 3 4 5 6 7 ;
{
  isamu="$samset.gsnap$i.samu"
  inam=`basename $isamu .samu`
  if ! test -f $isamu ; then  echo "missing $isamu"; continue;  fi

  ( $bindir/samtools view -u -t $workd/genome/$dgenome.fa.count  $isamu \
    | $bindir/samtools sort - $inam.sort ) &
}

wait

for i in  8 9 10 11 12 13 14 15 ;
{
  isamu="$samset.gsnap$i.samu"
  inam=`basename $isamu .samu`
  if ! test -f $isamu ; then  echo "missing $isamu"; continue;  fi

  ( $bindir/samtools view -u -t $workd/genome/$dgenome.fa.count $isamu \
    | $bindir/samtools sort - $inam.sort ) &
}

wait

$bindir/samtools merge $samset.bam $samset*.sort.bam

## redo this w/ samstats.pl, this not useful
$bindir/samtools flagstat $samset.bam > $samset.stats
# echo -n 'stranded reads: ' >> $samset.stats
# $bindir/samtools view $samset.bam | grep -c 'XS:A:' >> $samset.stats

