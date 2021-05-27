#! /bin/bash
### env fastdir=fastq/ qsub -q normal gsnp2genome.sh
#PBS -A ind114
#PBS -N gsnp4
#PBS -l nodes=2:ppn=16,walltime=19:55:00
#PBS -o gsnp1.$$.out
#PBS -e gsnp1.$$.err
#PBS -V

ncpu=32
npart=32
jstart=0
# jend=$(( $jstart + $ncpu ))

# run on gordo, use persist datad n
datad=$HOME/scratchn
workd=$datad/chrs/kfish
bindir=$HOME/bio/bin
gmapd=$HOME/bio/gmap1206
rund=$workd/rnas
gdb=gmap12
gdbpath=$workd/genome/

# snpd=snp10
# dgenome=dmagna20100422assembly; gtag=dmag2
dgenome=killifish20121129asm; gtag=kfish2
dgenosize=$dgenome.fa.count

#5 for genome map, .sam out, with introns, no snps print
# .. should add known intron database
# .. test for chasity: both, either, none
optsnp=""
# "--use-snps=$snpd "
## dmag fq4: optreads="--pairexpect=475 --pairdev=50"
## fundgr: dont know guess 200
optreads="--pairexpect=200 --pairdev=50"
optgenome="-N 1  " 

snapopt1="$optgenome $optsnp $optreads" 

cd $workd/rnas/
if [ "X$fastdir" = "X" ]; then
  fastdir="fastq"
fi

## tie outdir to fastdir or not?
outs=`basename $fastdir`
outdir=$workd/rnas/gsnout$outs
mkdir $outdir

#.. drop/add overloop for fastdir sets; edit per submit ..
# for fastdir in fastq4; do

#... START LOOP rnain : limit to suffix .fq, .fastq, sequence.txt ..

for rnain in $fastdir/*_1.fq.gz ;  do 
{ 
  cd $rund/
  if ! test -f "$rnain" ; then continue; fi
  
  drna=`basename $rnain .gz | sed 's/\.fastq//; s/\.fq//; '`
  outna=$drna-$gtag

  snapopts="--gunzip $snapopt1"
  
  rnain2=`echo $rnain | sed 's/1.fq/2.fq/; s/1.fastq/2.fastq/; '`
  if test -f "$rnain2" ; then
    inset="$rnain $rnain2"
  else 
    inset=$rnain
  fi

  # i=1,ncpu, j=jstart,jend  subset of npart
  i=0;  
  j=$jstart;
  while [ $i -lt $ncpu ]; do { 
  
   $gmapd/bin/gsnap $snapopts  -A "sam"  --part=$j/$npart \
    -D $gdbpath/$gdb -d $dgenome $inset > $outdir/$drna.gsnap$i.samu &
  
   i=$(( $i + 1 ))
   j=$(( $j + 1 ))
  }
  done
  wait

  ##...... sam sort/merge........

  samset=$drna
  cd $outdir
  i=0; while [ $i -lt $ncpu ]; do {
    isamu="$samset.gsnap$i.samu"
    inam=`basename $isamu .samu`
    if ! test -f $isamu ; then  echo "missing $isamu"; continue;  fi
  
    ( $bindir/samtools view -u -t $gdbpath/$dgenosize $isamu | $bindir/samtools sort - $inam.sort; )&
    # .sort may fail#  rm $isamu 
  
    i=$(( $i + 1 ))
  }
  done
  wait

  $bindir/samtools merge $outna.bam $samset*.sort.bam
  $bindir/samtools flagstat $outna.bam > $outna.fstat
  $bindir/samtools index $outna.bam 
} 
done

#... END LOOP rnain

