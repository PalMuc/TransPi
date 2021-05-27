#!/bin/bash
##  env subset=sc6 qsub -q normal velrun4g.sh
#PBS -N velrun 
#PBS -A ind114
#PBS -l nodes=1:ppn=16,walltime=29:55:00
#PBS -o velrun.$$.out
#PBS -e velrun.$$.err
#PBS -V

# runversion
# v1: all raw reads, fq
# v2: diginorm, 30%, fa
dv=2

ncpu=7
export OMP_NUM_THREADS=$ncpu
export OMP_THREAD_LIMIT=$ncpu

datad=$HOME/scratchg/
## bin1 for simple, kmer<60 .. bin2/ for kmer<106
velbin=$HOME/bio/velvet127s/bin2
workd=$datad/chrs/aabugs/tsa/ztick

velfad=$workd/sraf
rund=$workd/velv$dv
cd $workd

#v1 : fail out of mem but for k95 = 25k tr
# inpe=$velfad/pe/SRR*.fastq.gz
#v2
inpe=$velfad/SRR*.dnorm30.fa

## ave ins 200 for this 
vopts="-ins_length 300 "
oopts="-min_pair_count 2 -scaffolding yes -min_trans_lgth 180 -ins_length 300"
kset="95 85 75 65 55 45 35 29 25 21"

echo "START `date` " 
mkdir $rund
cd $rund/

#..... run loop for vel kmer steps
shopt -s nullglob

kseqdir=vel${dv}_seq
if [ ! -f $kseqdir/Sequences ]; then
 $velbin/velveth $kseqdir 27 -noHash \
    -fasta -shortPaired $inpe 
fi
#v1:    -fastq.gz -shortPaired $inpe  -short $insr
#   -shortPaired -fasta.gz -interleaved $velfad/$inreads
#   -shortPaired -fastq.gz  -separate $velfad/SRR346404_[12].fastq.gz 

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

