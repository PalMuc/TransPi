#! /bin/bash
### env rnain=bams/one.bam qsub -q batch bamcuff.sh
#PBS -N bamcuff
#PBS -A ind114
#PBS -l nodes=2:ppn=16,walltime=21:55:00
#PBS -o bamcuf1.$$.out
#PBS -e bamcuf1.$$.err
#PBS -V

# less than max for slop?
ncpu=32
# gmap set # bamdir=bama
# tophat set # bamdir=bamt

datad=$HOME/scratchn
workd=$datad/chrs/kfish/rnas
# module add cufflinks
# bindir=/home/diag/opt/cufflinks/1.3.0/bin
bindir=$HOME/bio/cufflinks202
logf="log.bamcuf$$"
cuft=cuf2gs

# --max-intron-length 300000 def;  --min-intron-length 50  def;
# have some intr>300k but maybe not real; minint 30 seems real
opts="-p $ncpu -v --no-update-check --min-intron-length 30"

#** cufflinks aint reading all .bam files: need  merged allbam3.bam
# samtools merge -1 bama/allbam3.bam bam3/*.bam >&log.bama
# .. and cut by chr subsets .. all.bam too big.

## needs to be full path to one.bam
drnabam=$workd/$rnain

drna=`echo $drnabam | sed 's,^.*/,,; s/.bam.*//;' `
if [ "${rname}x" != "x" ]; then drna=$rname; fi
odir=$cuft$drna

cd $workd/
mkdir $odir
cd $odir

touch $logf
echo "start bamcuff.$drna : `date`"  >> $logf

echo $bindir/cufflinks $opts $drnabam  >> $logf
$bindir/cufflinks $opts $drnabam  >> $logf 

echo "end bamcuff: `date`" >> $logf

