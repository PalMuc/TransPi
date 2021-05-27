#! /bin/bash
### env icut=0..3 rnain=xxx.1.fastq.gz rnain2=xxx.2.fastq.gz qsub -q normal gsnapbig.sh
#PBS -N gsnapbig
#PBS -A ddp138
#PBS -l nodes=1:ppn=32,walltime=23:55:00
#PBS -o gsnapbig.$$.out
#PBS -e gsnapbig.$$.err
#PBS -V

# .. took 8hr to run on 24 GB read set (2 x 12 GB pair files); then failed at end
# .. fails again after 8hr run; gsnap completed on all parts; 
# .. samsort filled up 128 GB ss disk
# >> cut to 4 parts; do those as sep. complete runs to .bam output

# env icut=0,1,2,3
# leave open 1 cpu for overhead
ncut=4
ncpu=31
npart=$(( $ncpu * $ncut))

if [ "X$icut" = "X" ]; then 
  echo "usage: env icut=0,1,2,3 rnain=xxx qsub gsnapbig.sh"; exit -1; 
fi
if [ "X$rnain" = "X" ]; then 
  echo "usage: env icut=0,1,2,3 rnain=xxx qsub gsnapbig.sh"; exit -1; 
fi

ecut=$(( $icut + 1 ))
jstart=$(( $ncpu * $icut ))
# note jend= +1 of final j, using -lt $jend below
jend=$(( $ncpu * $ecut ))
# echo i=$icut jb=$jstart je=$jend np=$npart


datad=/oasis/$USER
workd=$datad/chrs/cacao
bindir=$datad/bio/bin
gmapd=$datad/bio/gmap118
rund=/scratch/$USER/$PBS_JOBID
gdb=gmap118

dgenome=cacao11allasm
dgenosize=$dgenome.chr_size.txt

# drop this for inputs=path/to/rnain,rnain2 
fastqdir=$workd/rnas/fastqbig

drna=`basename $rnain .fastq | sed 's/.gz//; s/.fastq//; s/.fq//; s/_sequence.txt//; ' `
outf=$drna-i$icut
notef=$workd/rnas/$drna-i$icut.$$.RUNNING
donef=`echo $notef | sed 's/RUNNING/DONE/'`

# --gmap-mode=none  # maybe fix memory faults
snapopts="-N 1 --gmap-mode=none --quality-protocol=illumina "
if [ `basename $rnain .gz` != $rnain ]; then 
  snapopts="--gunzip $snapopts"
fi

touch $notef
echo "START to $outf" >> $notef
echo `date`  >> $notef
echo icut=$icut jb=$jstart je=$jend np=$npart >> $notef

mkdir -p $rund
mkdir $rund/sams
mkdir $workd/rnas/bamout

cd $workd/rnas/

cd $fastqdir/
if test -f "$rnain2" ; then
  cp -p  $fastqdir/$rnain2 $rund/
  inset="$rnain $rnain2"
else 
  inset=$rnain
fi

cp -p   $fastqdir/$rnain  $rund/
cp -p   $workd/genome/$dgenosize $rund/
cp -rp  $workd/genome/$gdb $rund/

cd $rund/
echo "START gsnap $inset" >> $notef
echo `date`  >> $notef
du -h >> $notef
ls -l >> $notef

j=$jstart; 
##while [ $j -lt $jend ]; do { # now jend-jstart == ncpu ; dont need this loop
  i=0; while [ $i -lt $ncpu ]; do { 
  
   $gmapd/bin/gsnap $snapopts -A "sam" --part=$j/$npart \
    -D ./$gdb -d $dgenome $inset > sams/$drna.gsnap$j.samu &
  
   i=$(( $i + 1 ))
   j=$(( $j + 1 ))
  } 
  done
  echo "wait gsnap to $j/$npart" >> $notef
  
  wait
  echo "DONE gsnap to $j/$npart" >> $notef
  echo `date`  >> $notef

##} done

echo "DONE gsnap" >> $notef
du -h >> $notef
ls -l sams >> $notef
#..................

# ** failing here now .. dont know why ..
# ** was TOO MANY CPU - i..npart samtools calls, not ncpu **
# ** fail2 .. filled up 128 GB ss disk; gsnap to .samu completed ok
# remove $isamu when done ..

samset=$drna
cd $rund/sams/

j=$jstart; 
##while [ $j -lt $jend ]; do { # now jend-jstart == ncpu ; dont need this loop
  i=0; while [ $i -lt $ncpu ]; do { 
  
    isamu="$samset.gsnap$j.samu"
    inam=`basename $isamu .samu`
    if ! test -f $isamu ; then  echo "missing $isamu"; continue;  fi
  
    ( $bindir/samtools view -u -t $rund/$dgenosize $isamu | $bindir/samtools sort - $inam.sort ; rm $isamu ) &
  
   i=$(( $i + 1 ))
   j=$(( $j + 1 ))
  } 
  done
  
  wait
  echo "done part samsort $j/$npart" >> $notef
  echo `date`  >> $notef

## } done

echo "DONE samsort" >> $notef
du -h >> $notef
#..................

$bindir/samtools merge $outf.bam $samset*.sort.bam
echo "DONE merge" >> $notef
cp -p $outf.bam $workd/rnas/bamout/
ls -l *.bam  >> $notef
rm $samset*.sort.bam

$bindir/samtools flagstat $outf.bam > $outf.bam.flagstat
cp -p $outf.bam.flagstat $workd/rnas/bamout/

cd $rund/

echo "DONE " >> $notef
echo `date`  >> $notef
mv $notef $donef
