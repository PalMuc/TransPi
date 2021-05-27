#!/bin/bash
##  env dv=2a PI=90 aaset=aadir/xxx.aa qsub -q normal cdhitrun.sh
#PBS -N cdhitr
#PBS -A ind114
#PBS -l nodes=1:ppn=16,walltime=19:55:00
#PBS -o cdhitr.$$.out
#PBS -e cdhitr.$$.err
#PBS -V

#.. see evigene/scripts/prot/aabest.sh
## PI=percent ident; PA=percent align
##  aa-cdhit:   cd-hit -c 0.$PI -l 30 -d 0 -i name.aa -o name.cd$PI.aa 
## cds-cdhit:   cd-hit-est -c 0.99 -G 0 -aS 0.$PA -l 90 -d 0 -i name.cds -o name.cde$PA.cds

# runversion
if [ "X" = "X$dv" ]; then dv=1; fi
if [ "X" = "X$PI" ]; then PI=90; fi

## doesnt use all this mem; 800 MB default, 8000 MB likely enough
mem=50000 # MB units
ncpu=15
export OMP_NUM_THREADS=$ncpu
export OMP_THREAD_LIMIT=$ncpu

cdhitapp=$HOME/bio/cdhit/bin/cd-hit
evigene=$HOME/bio/evigene

datad=$HOME/scratchg/
workd=$datad/chrs/aabugs/tsa/evigene
rund=$workd/ocdhit$dv
suf="_cd$PI"

cd $workd
mkdir $rund
#NO# cd $rund/

## add loop for all aadir/*.aa
if [ -d $aaset ]; then
  aaset=$aaset/*.aa
fi

for aaone in $aaset; do {

 onam=`basename $aaone .aa`
 onam=$rund/${onam}$suf

 echo $cdhitapp -c 0.$PI -T $ncpu -M $mem -d 0 -i $aaone -o $onam.aa 
 $cdhitapp -c 0.$PI -T $ncpu -M $mem -d 0 -i $aaone -o $onam.aa >& $onam.log
 rm $onam*.bak.clstr; 
 ## sort: write failed: standard output: Broken pipe; from head -1000
 env stat=1 $evigene/scripts/prot/aaqual.sh  $onam.aa
 ## gzip --fast $onam.aa

} done


