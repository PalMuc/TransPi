#! /bin/bash
### env fain=allpe.fa datad=`pwd` qsub -q batch dignorm.sh
#PBS -N dignorm
#PBS -A ind114
#PBS -l mem=40gb,nodes=1:ppn=2,walltime=28:55:00
#... ask for mem=64gb or not? is que set for that?
#PBS -o dignorm.$$.out
#PBS -e dignorm.$$.err
#PBS -V

# use env param? need '-p' for pe, not for sr
if [ "X" = "X$fain" ]; then
  echo "env fain=pair.fasta missing "; exit -1;
fi

if [ "X" = "X$datad" ]; then
  echo "env datad=path/to/data missing"; exit -1;
  # datad=$HOME/scratchn/chrs/aabugs/tsa/ztick/sraf;
fi

ncpu=2
kopt="-p";
# mem:  8e9 x 4 = 32Gb ; 10e9 x 4 = 40Gb; 64gb = 16e9
#62gb# xhash=15.6e9
xhash=9.8e9
#v1# xhash=8e9
#v2: reduce size
kmer=23
rkeep=10
#v1# kmer=30
#v1# rkeep=20

#** requires python2.6+ not on gordon.sdsc, on trestles module add python
module() { eval `/opt/modules/Modules/3.2.5/bin/modulecmd bash $*`; }
module add python

export PYTHONPATH=$HOME/bio/khmer/python
normapp=$HOME/bio/khmer/scripts/normalize-by-median.py

if [ "X" = "X$datad" ]; then
  echo "env datad=path/to/data missing"; exit -1;
  # datad=$HOME/scratchn/chrs/aabugs/tsa/ztick/sraf;
fi

cd $datad

if [ -f $fain ]; then
 $normapp $kopt -C $rkeep -x $xhash  -k $kmer  $fain
 fao=`echo $fain | sed 's/\.fa2//; s/\.fasta//; s/\.fa//;'`
 mv $fain.keep $fao.dnorm$kmer.fa
fi

