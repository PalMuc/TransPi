#! /bin/bash
### qsub -q normal gsnapset.sh
#PBS -N gsnapset
#PBS -A ddp138
#PBS -l nodes=1:ppn=32,walltime=23:55:00
#PBS -o gsnapset.$$.out
#PBS -e gsnapset.$$.err
#PBS -V

ncpu=32

datad=/oasis/$USER
workd=$datad/chrs/cacao
bindir=$datad/bio/bin
gmapd=$datad/bio/gmap118
rund=/scratch/$USER/$PBS_JOBID
gdb=gmap118

dgenome=cacao11allasm ; gtag=mars11
dgenosize=$dgenome.chr_size.txt

# dgenome=xxx ; gtag=cirad1c
# dgenosize=$dgenome.chr_size.txt

# drna=`basename $rnain .fastq | sed 's/.gz//; s/.fq//; s/.fastq//; s/_sequence.txt//; ' `
notef=$workd/rnas/gsnapset.$$.RUNNING
donef=`echo $notef | sed 's/RUNNING/DONE/'`

#  --gmap-mode=none maybe will cure seg faults?
snapopt1="-N 1 --gmap-mode=none --quality-protocol=illumina "

if [ "X$fastdir" = "X" ]; then
  fastdir="fastq"
fi

touch $notef
echo "START " >> $notef
echo `date`   >> $notef

mkdir -p $rund
mkdir $workd/rnas/bamout

cd $rund/
cp -p $workd/genome/$dgenosize $rund/
cp -rp $workd/rnas/$fastdir $rund/fastq
cp -rp $workd/genome/$gdb $rund/

du -h >> $notef
ls -l >> $notef
ls -l fastq/ >> $notef

#... START LOOP rnain : limit to suffix .fq, .fastq, sequence.txt ..

for rnain in fastq/*{.1.fastq,_1_sequence.txt}.gz ;  do 
{ 

  cd $rund/
  /bin/rm -rf $rund/sams
  mkdir $rund/sams
  
  if ! test -f "$rnain" ; then continue;  fi
  
  drna=`basename $rnain .fastq.gz | sed 's/.gz//; s/.fastq//; s/.fq//; s/_sequence.txt//;'`
  outna=$drna-$gtag
  echo "START gsnap : $rnain to $outna.bam" >> $notef
  echo `date`  >> $notef

  snapopts="--gunzip $snapopt1"
  
  rnain2=`echo $rnain | sed 's/.1.fastq/.2.fastq/; s/_1_sequence/_2_sequence/;'`
  if test -f "$rnain2" ; then
    inset="$rnain $rnain2"
  else 
    inset=$rnain
  fi

  i=0; while [ $i -lt $ncpu ]; do { 
  
   $gmapd/bin/gsnap $snapopts -A "sam" --part=$i/$ncpu \
    -D ./$gdb -d $dgenome $inset > sams/$drna.gsnap$i.samu &
  
   i=$(( $i + 1 ))
  }
  done
  echo "wait gsnap : $drna" >> $notef
  
  wait

  echo "DONE gsnap : $drna" >> $notef
  du -h >> $notef
  ls -l sams >> $notef
#..................
  
  samset=$drna
  cd $rund/sams/
  
  i=0; while [ $i -lt $ncpu ]; do {
    isamu="$samset.gsnap$i.samu"
    inam=`basename $isamu .samu`
    if ! test -f $isamu ; then  echo "missing $isamu"; continue;  fi
  
    ( $bindir/samtools view -u -t $rund/$dgenosize  $isamu | $bindir/samtools sort - $inam.sort ) &
  
    i=$(( $i + 1 ))
  }
  done
  
  wait

  echo "DONE bam sort" >> $notef
  du -h >> $notef
  #..................
  
  $bindir/samtools merge $outna.bam $samset*.sort.bam
  echo "DONE merge to $outna.bam" >> $notef
  cp -p $outna.bam $workd/rnas/bamout/
  ls -l *.bam  >> $notef
  
  $bindir/samtools flagstat $outna.bam > $outna.bam.flagstat
  cp -p $outna.bam.flagstat $workd/rnas/bamout/
  cd $rund/

}
done
#... END LOOP rnain

echo "DONE " >> $notef
echo `date`  >> $notef
mv $notef $donef
