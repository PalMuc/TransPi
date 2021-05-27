#! /bin/bash
### env  subd=bamq5treat qsub -q batch bam2bed.sh
#PBS -N bam2bed
#... -A ind114
#PBS -l nodes=1:ppn=24,walltime=22:55:00
#PBS -o bam2bed.$$.out
#PBS -e bam2bed.$$.err
#PBS -V

ncpu=24

## run1: subdirs by treatment/clone group ; combine replicates
subd=bamq5treat

datad=$HOME/scratch
workd=$datad/chrs/daphmag/
rund=$workd/rnas/$subd
dgenomesize=$workd/genome/dmagna20100422assembly.fa.fai

bindir=$HOME/bio/bin
sdir=$HOME/bio/evigene/scripts/rnaseq

# export samtools or put on path
module add samtools

cd $rund
treats=`ls -d HS*{X,I} ND*X`

echo "start bam2bed.$subd : `date`"  
i=0; j=1; 

for tdir in $treats; do {
  nam=`basename $tdir`
  echo "bam2bw $tdir/*.bam TO $nam.$j.bed " 

  # ** add normfactor : calc relative min/max of groups = maxgroupreads/thisgroupreads
  opts="debug=1 log=1 ave=1"
  if test -f $tdir.norm; then
    tnorm=`cat $tdir.norm`
    opts="debug=1 log=1 ave=0 norm=$tnorm"
  fi
  (env $opts $sdir/bam2bw.pl $dgenomesize  $tdir/*.bam > $nam.$j.bed ; gzip --fast $nam.$j.bed; ) &

  j=$(( $j + 1 ))
  i=$(( $i + 1 ))
  if [ $i -ge $ncpu ]; then wait; i=0; fi

} done
wait

echo "end bam2bed: `date` "
