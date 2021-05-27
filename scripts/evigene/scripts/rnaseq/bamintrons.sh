#! /bin/bash
### env rnain=bams/xxx.bam qsub -q normal bamintrons.sh
#PBS -N bamintr
#... -A ind114
#PBS -l nodes=1:ppn=32,walltime=2:55:00
#PBS -o bamintr.$$.out
#PBS -e bamintr.$$.err
#PBS -V

# takes only few minutes/25Mill bam
ncpu=32

# seq align identity to collect intron
inident=90

datad=$HOME/scratch
workd=$datad/chrs/daphmag/rnas/bam4
bindir=$HOME/bio/bin
sdir=$HOME/bio/evigene/scripts/rnaseq
logf=$workd/log.bamintr$$

cd $workd/

# drnabam=`ls tpout*/accepted_hits.bam`
# drna=tophatb
drnabam=`ls *.bam`
drna=gsnap4

outdir=$workd/introni$inident
mkdir $outdir
touch $logf
echo "start bamintr.$drna : `date`"  >> $logf
i=0; j=1; 

for bam in $drnabam ;  do { 

  nam=`basename $bam .bam`
  #tophat# nam=`echo $bam | sed 's,/accepted_hits.bam,,;'`
  echo "$bam TO $nam.$j.intron.gff " >> $logf
  ( $bindir/samtools view -F 0x104 $bam | grep 'XS:A:' | \
    $sdir/samintrons.pl -identity $inident -source rs$nam -sort > $outdir/$nam.$j.intron.gff ; ) &
  
  j=$(( $j + 1 ))
  i=$(( $i + 1 ))
  if [ $i -ge $ncpu ]; then wait; i=0; fi

} 
done
wait

echo "end bamintr: `date` " >> $logf

