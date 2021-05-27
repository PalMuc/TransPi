#! /bin/bash
### env trdb=trseqs/xxxx reads=fq06/ qsub -q batch bow2trsets.sh
#PBS -N bowtrmap
#PBS -l nodes=1:ppn=28,walltime=3:55:00
#PBS -o bowtrmap.$$.out
#PBS -e bowtrmap.$$.err
#PBS -V
## bow2trmapc.sh for cacao/genes

#.. runs fast; 5..10min using 4cpu per readset (2 Mill reads)
#.. add reads loop .. add samtools sort/merge

ncpu=28
ncpuPerMap=4

module add bowtie/2.0.0_b6
module add samtools

## input opt
# trdb=trseqs/cacao3tri4asm_cd.goodset

# .. 7 read sets x 4cpu/set
datad=$HOME/scratch
workd=$datad/chrs/cacao/genes
odir=bowtr2

## .. 7 read pair sets per fqxx/
## fq06/tc06l1_GCCAAT_1.fastq.gz tc06l1_GCCAAT_2.fastq.gz tc06l2_GCCAAT_1.fastq.gz tc06l2_GCCAAT_2.fastq.gz

#bowtie2: .. only sam output, -S outfile.sam?; keep unal in .sam, skip -un x.fq
opts="-q --threads $ncpuPerMap --fast -k 20 --no-head "

cd $workd/
mkdir $odir

trsize=`echo $trdb | sed 's/.tr$// s/$/.count/'`
trdb=`echo $trdb | sed 's/.tr$//;'`
tnam=`basename $trdb`
outna=$odir/$tnam-tc67a

# add reads dir loop here .. or put all in 1 dir? need cpu counter..
for reads in fq06 fq07 fq10 ; do {

readset=`ls $reads/*_1.fastq.gz`

i=0;
for lreads in $readset; do {
  rreads=`echo $lreads | sed 's/_1/_2/;'`
  rdnam=`basename $lreads _1.fastq.gz`
  
  echo bowtie2 $opts -x $trdb -1 $lreads -2 $rreads -S $odir/$tnam-$rdnam.sam 2>$odir/$tnam-$rdnam.blog
  bowtie2 $opts -x $trdb -1 $lreads -2 $rreads -S $odir/$tnam-$rdnam.sam 2>$odir/$tnam-$rdnam.blog  &

  i=$(( $i + $ncpuPerMap ))
  if [ $i -ge $ncpu ]; then wait; i=0; fi

} done
wait
# reads loop

} done
# ..redir loop

##...... sam sort/merge........
samset=`ls $odir/$tnam-*.sam`

i=0;
for isamu in $samset; do {
  inam=`echo $isamu | sed 's/.sam//'`
  if ! test -f $isamu ; then  echo "missing $isamu"; continue;  fi

  ( samtools view -u -t $trsize $isamu | samtools sort - $inam.sort; ) &
  # rm $isamu .. later

  i=$(( $i + 1 ))
  if [ $i -ge $ncpu ]; then wait; i=0; fi
  
} done
wait

samtools merge $outna.bam $odir/$tnam*.sort.bam
samtools index $outna.bam
samtools idxstats $outna.bam > $outna.bam.count 

echo /bin/rm $odir/$tnam-*.sam
echo /bin/rm $odir/$tnam-*.sort.bam
