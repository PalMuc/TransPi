#! /bin/bash -l
### env rnain=aphidrs_SRR071347  qsub -q normal rnagsnap11.sh
#PBS -N gsnap1
#PBS -A ddp138
#PBS -l nodes=1:ppn=32,walltime=2:55:00
#PBS -o gsnap1.$$.out
#PBS -e gsnap1.$$.err
#PBS -V

ncpu=32
## ^ 17min split 32 ways for SRR063707
# workd=$HOME/scratch/chrs/nasv1
# gmapd=$HOME/bio/gmap11

datad=/oasis/$USER
workd=$datad/chrs/nasv1
gmapd=$datad/bio/gmap11
rund=/scratch/$USER/$PBS_JOBID

dgenome=nasvit1asm

# suf=txt  fa2 fastq fq
drna=`basename $rnain .fastq | sed 's/.gz//; s/.fq//; s/.fastq//;' `

# snapopts="-N 1 -k 14 --max-mismatches=0.07 --pairmax=15000"
## gmap 2011: new opts, using fastq sr input;  --gunzip
## damn crap dont know which qualprot; how am i supposed to know?
#  snapopts="--gunzip -N 1 -k 14 --quality-protocol=sanger "
snapopts="--gunzip -N 1 -k 14 "

# cd $workd/rnas/
mkdir -p $rund
cd $rund/
mkdir sams
mkdir gmap
cp -p $workd/rnas/fastq/$rnain $rund/
cp -rp $workd/genome/gmap $rund/gmap

i=0; while [ $i != $ncpu ]; do { 

 $gmapd/bin/gsnap $snapopts -A "sam" --part=$i/$ncpu \
  -D gmap -d $dgenome $rnain > sams/$drna.gsnap$i.samu &

 i=$(( $i + 1 ))
}
done

wait

cp -rp sams  $workd/rnas/

