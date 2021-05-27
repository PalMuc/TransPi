#! /bin/bash -l
### env rnain=aphidrs_SRR071347  qsub -q normal rnagsnap11.sh
#PBS -N gsnap1
#PBS -A ddp138
#PBS -l nodes=1:ppn=32,walltime=2:55:00
#PBS -o gsnap1.$$.out
#PBS -e gsnap1.$$.err
#PBS -V

ncpu=32
workd=$HOME/scratch/chrs/aphid2
gmapd=$HOME/bio/gmap11
# dgenome=aphid2asm ; gtag=asm2
dgenome=aphid2DGILpub8tr ; gtag=a2pub8tr

drna=`basename $rnain .fastq.gz`
drna=`basename $drna .fa2.gz`

# snapopts="-N 1 -k 14 --quality-protocol=sanger "
snapopts="--gunzip --quality-protocol=sanger "

cd $workd/rnas/

# if ! test -f "$rnain" ; then  echo "missing $drna.$suf"; exit;  fi

if test -f "$rnain2" ; then
  inset="$rnain $rnain2"
else 
  inset=$rnain
fi

i=0; while [ $i != $ncpu ]; do { 

 $gmapd/bin/gsnap $snapopts -A "sam" --part=$i/$ncpu \
  -D $workd/genome/gmap -d $dgenome $inset > sams/$drna-$gtag.gsnap$i.samu &

 i=$(( $i + 1 ))
}
done

wait


