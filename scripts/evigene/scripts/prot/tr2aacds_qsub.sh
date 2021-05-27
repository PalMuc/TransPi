#! /bin/bash
### env trset=allshrimpt1.tr datad=`pwd` qsub -q normal tr2aacds_qsub.sh
#PBS -N tr2cds
#PBS -l nodes=1:ppn=32,walltime=18:55:00
#PBS -o tr2cds.$$.out
#PBS -e tr2cds.$$.err
#PBS -V

ncpu=30
maxmem=50000
evigene=$HOME/bio/evigene/scripts

#t2ac: app=cd-hit-est, path= echo MISSING_cd-hit-est
export PATH=$HOME/bio/cdhit/bin:$PATH
#t2ac: app=fastanrdb, path= echo MISSING_fastanrdb
export fastanrdb=$HOME/bio/exonerate/bin/fastanrdb
#t2ac: app=blastn, path= echo MISSING_blastn
export PATH=$HOME/bio/ncbi2227/bin:$PATH

if [ "X" = "X$trset" ]; then
  echo "missing env trset=xxxx.tr"; exit -1
fi
if [ "X" = "X$datad" ]; then
  echo "missing env datad=/path/to/data"; exit -1
fi

cd $datad/

echo $evigene/prot/tr2aacds.pl -debug -NCPU $ncpu -MAXMEM $maxmem -log -cdna $trset
$evigene/prot/tr2aacds.pl -debug -NCPU $ncpu -MAXMEM $maxmem -log -cdna $trset
