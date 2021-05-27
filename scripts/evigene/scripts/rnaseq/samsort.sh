#! /bin/bash -l
### qsub -q batch samsort.sh
#PBS -N samsort
#PBS -A ind114
#PBS -l nodes=1:ppn=8,walltime=01:55:00
#PBS -o sams1.out
#PBS -e sams1.err
#PBS -V

workd=$HOME/scratch/chrs/daphplx
bindir=$HOME/bio/bin

dgenome=daphplx06

cd $workd/rnas/parts

samset=`ls *.gsnap1.samu`
samset=`basename $samset .gsnap1.samu`

for i in  0 1 2 3 4 5 6 7 ;
{
  isamu="$samset.gsnap$i.samu"
  inam=`basename $isamu .samu`
  if ! test -f $isamu ; then  echo "missing $isamu"; continue;  fi

  ( $bindir/samtools view -u -t $workd/genome/$dgenome-sizes.txt $isamu \
    | $bindir/samtools sort - $inam.sort ) &
}

wait

for i in  8 9 10 11 12 13 14 15 ;
{
  isamu="$samset.gsnap$i.samu"
  inam=`basename $isamu .samu`
  if ! test -f $isamu ; then  echo "missing $isamu"; continue;  fi

  ( $bindir/samtools view -u -t $workd/genome/$dgenome-sizes.txt $isamu \
    | $bindir/samtools sort - $inam.sort ) &
}

wait

$bindir/samtools merge $samset.merge.bam $samset*.sort.bam

$bindir/samtools flagstat $samset.merge.bam > $samset.merge.stats
echo -n 'stranded reads: ' >> $samset.merge.stats
$bindir/samtools view $samset.merge.bam | grep -c 'XS:A:' >> $samset.merge.stats

