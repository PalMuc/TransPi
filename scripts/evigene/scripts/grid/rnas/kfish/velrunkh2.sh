#!/bin/bash
##  env subset=sc6 qsub -q normal velrun4g.sh
#PBS -N velrun 
#PBS -A ind114
#PBS -l nodes=1:ppn=16,walltime=29:55:00
#PBS -o velrun.$$.out
#PBS -e velrun.$$.err
#PBS -V

# runversion
dv=4mi; inreads=kf_mdibl_rnaseq_reads_[12].fq.gz
# dv=2w; inreads=kf_whoi_rnaseq_reads_[12].fq.gz
# kf_mdibl_rnaseq_reads_1.fq.gz  kf_whoi_rnaseq_reads_1.fq.gz

ncpu=7
export OMP_NUM_THREADS=$ncpu
export OMP_THREAD_LIMIT=$ncpu

datad=$HOME/scratchg/
# datad=$HOME/scratchn/
# kfishg=/oasis/scratch/ux455375/temp_project/chrs/kfish
velbin=$HOME/bio/velvet127s/bin2
workd=$datad/chrs/kfish/rnas
velfad=$workd/fastq
rund=$workd/vel$dv

## FIXME: inner insert = 180, but add 200 bp reads for outer insert span = 380 !!!
#v2#iopts="-ins_length 180 "
iopts="-ins_length 380 "
vopts="$iopts"
oopts=" -scaffolding yes -min_pair_count 3 -min_trans_lgth 180 $iopts"
#2m# kset="91 81 71 63 53 43 35 29 25"
#2w# kset="93 83 73 65 55 45 35 29 25"
#4mi# ## outamem at k53
kset="95 85 75 69 65 59 55 "

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
#   -shortPaired -fastq.gz  -separate $velfad/SRR346404_[12].fastq.gz 
#   -shortPaired -fasta.gz -interleaved $velfad/$inreads

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

