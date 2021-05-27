#!/bin/bash
##  env subset=sc6 qsub -q normal velrun4g.sh
#PBS -N velrun 
#PBS -A ind114
#PBS -l nodes=1:ppn=16,walltime=29:55:00
#PBS -o velrun.$$.out
#PBS -e velrun.$$.err
#PBS -V

# runversion
dv=3
inpe=all*.dnorm23.fa.gz
inlong=ztickmf1velsop_ok.tr.gz

datad=$HOME/scratchg/
velbin=$HOME/bio/velvet127s/bin2

iopts="-ins_length 300"
vopts="-conserveLong yes $iopts"
oopts="-merge no -scaffolding yes -min_pair_count 3 -min_trans_lgth 180 $iopts"
kset="95 85 75 65 55 45 35 25 21"

ncpu=7
export OMP_NUM_THREADS=$ncpu
export OMP_THREAD_LIMIT=$ncpu

workd=$datad/chrs/aabugs/tsa/ztick
cd $workd

velfad=$workd/sraf
rund=$workd/velv$dv

echo "START `date` " 
mkdir $rund
cd $rund/

#..... run loop for vel kmer steps
shopt -s nullglob

kseqdir=vel${dv}_seq
if [ ! -f $kseqdir/Sequences ]; then
 $velbin/velveth $kseqdir 27 -noHash \
    -shortPaired -fasta.gz -interleaved $velfad/$inpe \
    -long -fasta.gz $velfad/$inlong 

fi

for k in $kset;  do { 
  ksubdir=vel${dv}_$k
  echo "#.. start velrun $ksubdir : `date`"
  mkdir $ksubdir
  ln -s ../$kseqdir/Sequences $ksubdir/
  $velbin/velveth $ksubdir $k  -reuse_Sequences
  $velbin/velvetg $ksubdir $vopts -read_trkg yes 
  $velbin/oases   $ksubdir $oopts 
  
  /bin/rm $ksubdir/{Graph2,LastGraph,PreGraph,Roadmaps,Sequences}
  echo "#.. end velrun $ksubdir : `date`"
}
done

echo "DONE `date` " 

