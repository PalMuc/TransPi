#! /bin/bash
### env trdb=xxx reads=fqreaddir qsub -q normal bowtrmap2.sh
#PBS -N bowtrmap
#PBS -A ind114
#PBS -l nodes=4:ppn=32,walltime=21:55:00
#PBS -o bowtrmap.$$.out
#PBS -e bowtrmap.$$.err
#PBS -V

# ncpu=32
ncpu=128
ncpuPerMap=8

bowd=$HOME/bio/bowtie2
samd=$HOME/bio/bin
export PATH=$PATH:$bowd:$samd
# trestles
datad=$HOME/scratcht
workd=$datad/chrs/daphmag/rnas/

tnam=`basename $trdb .tr`
rnam=`basename $reads`
odir=bowtr$rnam
outna=$odir/$tnam-$rnam

## FIXME tr.count for sam
trsize=`echo $trdb | sed 's/$/.count/'`

opts="-q --threads $ncpuPerMap --fast "
## --end-to-end : default;  -M 5 : default
## --un-conc <path>      write pairs that didn't align concordantly

cd $workd/
mkdir $odir

readset=`ls $reads/*_1.fastq.gz`
## .. handle diff namings _1.fastq, _R1_001.fastq, ..no need, renamed input fastq to _[12].fastq.gz
# readset=`ls $reads/*_R1_*.fastq.gz`
# if [ "X" = "X$readset" ]; then  readset=`ls $reads/*_1.fastq.gz`;  fi

i=0;
for lreads in $readset; do {
  rreads=`echo $lreads | sed 's/_1\./_2./; '`
  rdnam=`basename $lreads .fastq.gz | sed 's/_1\..*//; '`
  if [ ! -f $rreads ]; then 
    echo "ERR: missing right _2.fq of $lreads"; continue;
  fi 

  echo bowtie2 $opts -x $trdb -1 $lreads -2 $rreads -S $odir/$tnam-$rdnam.sam 2>$odir/$tnam-$rdnam.blog
  bowtie2 $opts -x $trdb -1 $lreads -2 $rreads \
     --un-conc $odir/$rdnam.unc.fq -S $odir/$tnam-$rdnam.sam 2>$odir/$tnam-$rdnam.blog  &

  i=$(( $i + $ncpuPerMap ))
  ## this will hang up till all in batch done; really want: wait1; i=i-1
  if [ $i -ge $ncpu ]; then wait; i=0; fi

} done
wait
# reads loop

## exit
#..............................................
##...... sam sort/merge........
samset=`ls $odir/$tnam-*.sam`

i=0;
for isamu in $samset; do {
  inam=`echo $isamu | sed 's/.sam//'`
  if ! test -f $isamu ; then  echo "missing $isamu"; continue;  fi

  ( samtools view -u -t $trsize $isamu | samtools sort - $inam.sort; samtools index $inam.sort.bam; ) &

  i=$(( $i + 1 ))
  if [ $i -ge $ncpu ]; then wait; i=0; fi
  
} done
wait

#.. want per readset .bam more than merge ..
samtools merge $outna.bam $odir/$tnam*.sort.bam
samtools index $outna.bam
samtools idxstats $outna.bam > $outna.bam.count 

# echo /bin/rm $odir/$tnam-*.sam
# echo /bin/rm $odir/$tnam-*.sort.bam

