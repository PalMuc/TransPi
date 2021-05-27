#!/bin/bash
##  env qsub -q normal velrun1.sh
#PBS -N velrun 
#PBS -A ind114
#PBS -l nodes=1:ppn=16,walltime=29:55:00
#PBS -o velrun.$$.out
#PBS -e velrun.$$.err
#PBS -V

# runversion
# tiger shrimp, 72bp pe rnaseq illum., approx 33M x 3 pairs, outer insertsize=300
# zebra tick, 101bp pe illum, 100b inner insert

dv=1
ncpu=7
export OMP_NUM_THREADS=$ncpu
export OMP_THREAD_LIMIT=$ncpu

datad=$HOME/scratchg/
## bin1 for simple, kmer<60 .. bin2/ for kmer<106
velbin=$HOME/bio/velvet127s/bin2
workd=$datad/chrs/aabugs/tsa/ztick

velfad=$workd/sraf
rund=$workd/velv$dv
# all PE data for ztick
inpe=$velfad/pe/SRR*.fastq

## ave ins 300 (100rd x 2 + 100 inside)
vopts="-ins_length 300 "
oopts="-min_pair_count 2 -scaffolding yes -min_trans_lgth 180 -ins_length 300"
#kset="65 61 55 51 45 41 35 31 25 21"
kset="95 85 75 65 55 45 35 29"

echo "START `date` " 
cd $workd
mkdir $rund
cd $rund/

#..... run loop for vel kmer steps
shopt -s nullglob

kseqdir=vel${dv}_seq
if [ ! -f $kseqdir/Sequences ]; then
 $velbin/velveth $kseqdir 27 -noHash \
    -fastq -shortPaired $inpe  
fi
#  -fastq.gz -shortPaired $inpe  -short $insr
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

