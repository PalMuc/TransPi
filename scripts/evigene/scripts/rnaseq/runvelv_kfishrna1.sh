#!/bin/bash
##  env subset=sc6 qsub -q normal velrun4g.sh
#PBS -N velrun 
#...... PBS -A ind114
#PBS -l vmem=256gb,nodes=1:ppn=16,walltime=29:55:00
#PBS -o velrun.$$.out
#PBS -e velrun.$$.err
#PBS -V

# runversion
dv=1h

inreads=kf_*_rnaseq_reads_[12].fq.gz
inlong=
# kf_mdibl_rnaseq_reads_[12].fq.gz kf_whoi_rnaseq_reads_[12].fq.gz
# inreads=allkfishgr.dnorm30.fa.gz

ncpu=9
export OMP_NUM_THREADS=$ncpu
export OMP_THREAD_LIMIT=$ncpu

# velbin1: cat=3; kmax=61
# velbin2: cat=9; kmax=110
velbin=$HOME/bio/velvet127s/bin2
workd=$HOME/projectd/chrs/kfish/rnas
velfad=$workd/fastq
rund=$workd/vel$dv

vopts="-ins_length 200 "
oopts="-scaffolding yes -min_pair_count 3 -min_trans_lgth 180 -ins_length 200"

kset="95 85 75 65 53 43 35 29 25"

echo "START `date` " 
cd $workd
mkdir $rund
cd $rund/

#..... run loop for vel kmer steps
shopt -s nullglob

kseqdir=vel${dv}_seq
if [ ! -f $kseqdir/Sequences ]; then
 $velbin/velveth $kseqdir 27 -noHash \
    -shortPaired -fastq.gz -separate $velfad/$inreads 
fi
##    -long -fasta.gz $velfad/$inlong

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

