#! /bin/bash
### env fain=allpe.fa qsub -q batch dignorm.sh
#PBS -N dignorm
#PBS -A ind114
#PBS -l mem=36gb,nodes=1:ppn=2,walltime=28:55:00
#PBS -o dignorm.$$.out
#PBS -e dignorm.$$.err
#PBS -V

# use env param? need '-p' for pe, not for sr
# fain=allpe.fa; kopt="-p";
fain=SRR392744.pe.fa; kopt="-p";
# fain=allsr1.fa; kopt="";

kmer=30
# 8e9 x 4 = 32Gb
xhash=8e9
rkeep=20
ncpu=2

#** requires python2.6+ not on gordon.sdsc, on trestles module add python
module() { eval `/opt/modules/Modules/3.2.5/bin/modulecmd bash $*`; }
module add python

export PYTHONPATH=$HOME/bio/khmer/python
normapp=$HOME/bio/khmer/scripts/normalize-by-median.py
datad=$HOME/scratchn
workd=$datad/chrs/aabugs/tsa/ztick
fastad=$workd/sraf

cd $fastad

# ~/bio/sratoolkit/fastq-dump --fasta 0 SRR346404.sra  : paired spots
#  --split-spot, have 60bp  pair reads now; use -p

if [ -f $fain ]; then
 $normapp $kopt -C $rkeep -x $xhash  -k $kmer  $fain
 fao=`echo $fain | sed 's/\.fasta//; s/\.fa//;'`
 mv $fain.keep $fao.dnorm$kmer.fa
fi

