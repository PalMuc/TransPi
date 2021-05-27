#! /bin/bash -l
### env rnain=bams/xxx.bam qsub -q normal bamcuff.sh
#PBS -N bamcuff
#PBS -A ddp138
#PBS -l nodes=1:ppn=32,walltime=19:55:00
#PBS -o bamcuf1.$$.out
#PBS -e bamcuf1.$$.err
#PBS -V

# less than max for slop
ncpu=30

datad=/oasis/$USER
workd=$datad/chrs/nasv1
# bindir=$datad/bio/bin
# rund=/scratch/$USER/$PBS_JOBID
bindir=$datad/bio/cufflinks
logf="log.bamcuf$$"

dgenome=nasvit1asm

drnabam=$workd/rnas/bams/*.bam
if [ "${rnain}x" != "x" ]; then drnabam=$rnain; fi
drna=`echo $drnabam | sed 's,^.*/,,; s/.bam//;' `
if [ "${rname}x" != "x" ]; then drna=$rname; fi

# opts="-p $ncpu -v --no-update-check "
# ^ defaults, too many retained introns? try higher -j/--pre-mrna-fraction 0.25 <0.0-1.0> (0.15 default)
#  and --min-intron-length 35 (50 default)
#  maybe: --trim-3-dropoff-frac 0.20  (0.1 default)
opts="--pre-mrna-fraction 0.25 --min-intron-length 40 -p $ncpu -v --no-update-check "

cd $workd/rnas
mkdir cuff.$drna
cd cuff.$drna

touch $logf
echo "start bamcuff.$drna : `date`"  >> $logf

echo $bindir/cufflinks $opts $drnabam  >> $logf
$bindir/cufflinks $opts $drnabam  >> $logf 

echo "end bamcuff: `date`" >> $logf

