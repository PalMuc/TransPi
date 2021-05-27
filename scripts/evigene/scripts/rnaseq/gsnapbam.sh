#! /bin/bash
### env rnain=aphidrs_SRR071347  qsub -q normal gsnapbam.sh
#PBS -N gsnap1
#PBS -A ddp138
#PBS -l nodes=1:ppn=32,walltime=12:55:00
#PBS -o gsnap1.$$.out
#PBS -e gsnap1.$$.err
#PBS -V

ncpu=32

datad=/oasis/$USER
workd=$datad/chrs/nasv1
bindir=$datad/bio/bin
gmapd=$datad/bio/gmap11
rund=/scratch/$USER/$PBS_JOBID

dgenome=nasvit1asm
dgenosize=$dgenome.chr_size.txt

drna=`basename $rnain .fastq | sed 's/.gz//; s/.fq//; s/.fastq//;' `
notef=$workd/rnas/$drna.$$.RUNNING

# snapopts="--gunzip -N 1 -k 14 --quality-protocol=sanger "
# snapopts="-N 1 -k 14 "
snapopts="-N 1 --quality-protocol=illumina "
if [ `basename $rnain .gz` != $rnain ]; then 
  snapopts="--gunzip $snapopts"
fi

# cd $workd/rnas/
mkdir -p $rund
mkdir $workd/rnas/bamout/

cd $rund/
touch $notef
mkdir sams

cp -p $workd/rnas/fastq/$rnain $rund/
cp -p $workd/genome/$dgenosize $rund/
cp -rp $workd/genome/gmap $rund/

echo 'START gsnap' >> $notef
echo `date`  >> $notef
du -h >> $notef
ls -l >> $notef

i=0; while [ $i != $ncpu ]; do { 

 $gmapd/bin/gsnap $snapopts -A "sam" --part=$i/$ncpu \
  -D ./gmap -d $dgenome $rnain > sams/$drna.gsnap$i.samu &

 i=$(( $i + 1 ))
}
done

wait

echo 'DONE gsnap' >> $notef
du -h >> $notef
ls -l sams >> $notef
#..................

samset=$drna
cd $rund/sams/

i=0; while [ $i != $ncpu ]; do {
  isamu="$samset.gsnap$i.samu"
  inam=`basename $isamu .samu`
  if ! test -f $isamu ; then  echo "missing $isamu"; continue;  fi

  ( $bindir/samtools view -u -t $rund/$dgenosize  $isamu | $bindir/samtools sort - $inam.sort ) &

  i=$(( $i + 1 ))
}
done

wait

echo 'DONE bam sort' >> $notef
du -h >> $notef
#..................

$bindir/samtools merge $samset.bam $samset*.sort.bam
echo 'DONE merge' >> $notef
cp -p $samset.bam $workd/rnas/bamout/
ls -l *.bam  >> $notef

$bindir/samtools flagstat $samset.bam > $samset.bam.flagstat
cp -p $samset.bam.flagstat $workd/rnas/bamout/

echo 'DONE ' >> $notef
echo `date`  >> $notef

