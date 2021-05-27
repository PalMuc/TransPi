#! /bin/bash
### env allsam=xxx/*.sam genoasm=path/to/gasm.fa qsub -q normal sambambed.sh
#PBS -N sambambed
#PBS -l nodes=1:ppn=24,walltime=9:55:00
#PBS -o sambambed.$$.out
#PBS -e sambambed.$$.err
#PBS -V

module add samtools
evigene=$HOME/bio/evigene
# bam2bw=$evigene/scripts/rnaseq/bam2bw.pl
bedopts="debug=1 log=1 ave=1 opt=-A"
samsortop="-l 9 --threads 4 -m 4G"

## allsam == insam
if [ "X" = "X$insam" ]; then echo "ERR: miss insam=gsodaplx/xxx*_.{concordant_uniq,concordant_mult,..}"; exit -1; fi
if [ "X" = "X$genoasm" ]; then echo "ERR: miss genoasm=genoasm.fa"; exit -1; fi
if [ "X" = "X$datad" ]; then echo "ERR: datad=what?"; exit -1; fi

cd $datad

if [ "X" = "X$samo" ]; then
  samo=($insam)
  samo=$samo.srt
fi

# ?? need sort -T tmpfilename
samsortop="$samsortop -T $samo.ss"

cat $insam | samtools view -T $genoasm -S -u - | samtools sort $samsortop -o $samo.bam - ;
env $bedopts $evigene/scripts/rnaseq/bam2bw.pl $genoasm.fai  $samo.bam > $samo.bed ;
gzip --fast  $samo.bed; 

## usage
# env insam=gsodaplx16*rs/*halfmapping*{mult,uniq} samo=daplx16halfmaprs genoasm=$genoasm \
#    datad=`pwd` prog=./sam1bambed.sh sbatch ../../srun_share2.sh
# env insam=gsodaplx16*rs/*mult samo=daplx16multrs ..
# env insam=gsodaplx16*rs/*concordant_uniq samo=daplx16uniqrs ..

## new samtoo sort opts:
#  --threads 4 ; -m 4G = mem/thread
