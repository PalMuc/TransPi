#! /bin/bash
### env bamset=bams  datad=`pwd` qsub -q shared bamresort.sh
#... -A ind114
#PBS -N bamresort
#PBS -l vmem=32gb,nodes=1:ppn=16,walltime=7:55:00
#PBS -V

ncpu=16
MEM=5000000000
sortop="-m $MEM"
## convert namesort bam back to location sort, count reads/loc

if [ "X" = "X$bamset" ]; then echo "ERR: miss bamset=gsnapdataprefix"; exit -1; fi
if [ "X" = "X$datad" ]; then echo "ERR: miss datad=path/to/gsnapset "; exit -1; fi
bindir=$HOME/bio/bin
export PATH=$bindir:$PATH

cd $datad/
i=0;
for bamin in $bamset; do {
  outs=`echo $bamin | sed 's/.bam//; s/$/lc/;'`
  ( if [ ! -f $outs.bam ]; then samtools sort $sortop  $bamin $outs; fi; \  
  samtools index $outs.bam; samtools idxstat $outs.bam > $outs.count; ) &

  i=$(( $i + 1 )); if [ $i -ge $ncpu ]; then wait; i=0; fi

} done
wait

## sam2bam
#  samtools view -u -t $genosize -S $insam | samtools sort $sortop - $outs ;\
#   samtools flagstat $outs.bam > $outs.fstat; \
#   samtools index $outs.bam; samtools idxstat $outs.bam > $outs.count; 
