#!/bin/bash
##  env inpe="sraf/xxx.fa.gz" name="xxx"  qsub -q normal velrun3x.sh
#PBS -N velrun 
#PBS -A ind114
#PBS -l nodes=1:ppn=16,walltime=29:55:00
#PBS -o velrun.$$.out
#PBS -e velrun.$$.err
#PBS -V

# runversion
dv=3
#below# inpe=all*.dnorm23.fa.gz
#below# inlong=ztickmf1velsop_ok.tr.gz

datad=$HOME/scratchg/
workd=$datad/chrs/aabugs/tsa/ztick
velbin=$HOME/bio/velvet127s/bin2

if [ "X" = "X$name" ]; then name=ztick; fi
rund=$workd/velv$name$dv

# all PE data for ztick; now .fa nogz
if [ "X" = "X$inpe" ]; then echo "env inpe=MISSING"; exit -1; fi
inpe=`echo $inpe | sed "s,sraf,$workd/sraf,g;"`
## FIXME:
inlong="$workd/sraf/ztickmf1velsop_ok.tr.gz"

iopts="-ins_length 300"
vopts="-conserveLong yes $iopts"
oopts="-merge no -scaffolding yes -min_pair_count 3 -min_trans_lgth 180 $iopts"
#3x.kset="85 75 65 55 45 35 31"
kset="71 61 51 41 33"

ncpu=7
export OMP_NUM_THREADS=$ncpu
export OMP_THREAD_LIMIT=$ncpu

cd $workd

echo "START `date` " 
mkdir $rund
cd $rund/

#..... run loop for vel kmer steps
shopt -s nullglob

kseqdir=vel${dv}_seq
if [ ! -f $kseqdir/Sequences ]; then
 $velbin/velveth $kseqdir 27 -noHash \
    -shortPaired -fasta -interleaved $inpe \
    -long -fasta.gz $inlong 

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

