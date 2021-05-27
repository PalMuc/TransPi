#!/bin/bash
##  env subset=sc6 qsub -q normal velrun4g.sh
## vel4g/trs/velrun4g.sh
#PBS -N velrun 
#sdscg .. PBS -A ind114
#sdscg .. PBS -l nodes=1:ppn=16,walltime=21:55:00
#PBS -l vmem=128gb,nodes=1:ppn=10,walltime=21:55:00
#PBS -o velrun.$$.out
#PBS -e velrun.$$.err
#PBS -V

#.. vel40: rerun with old velvet10 (1.0..) similar to vel4 run last fall .. is it soft vers bug?
#.. vel43: velvet10 is not answer, revert newest vel/o, use only .fa2 paired data
#.. vel43 best so far for new runs; nearing cuff08, trin qual.
#.. vel44: higher kmers, maybe add .longfa but not short .fa1 (error prone)
#.. vel45: up kmermax=101 velvet/bin2, cats=3 for short/long/longer inserts
## ins av/sd cgb:204/57; ncgr 090511=115/71; 090922_8=292/240; 091005_2=306/196; 090609_3=382/199
#v4# -conserveLong no  -ins_length 200 -ins_length_sd 40 -ins_length2 350 -ins_length2_sd 90
#-------

# runversion
dv=12
# data: cacao/rnas/velw/vel4cfa/ ; as per trsoap7, trin1
datasub=trim.$subset
if [ "X$subset" = "X" ]; then
  echo "ERR: missing env subset=XXX"; exit -1
fi

ncpu=8
export OMP_NUM_THREADS=$ncpu
export OMP_THREAD_LIMIT=$ncpu

datad=$HOME/scratch/
# datad=/oasis/scratch/$USER/temp_project
# rund=/oasis/scratch/$USER/$PBS_JOBID

## bin2: kmermax=105 insert cat=3
velbin=$HOME/bio/velvet127s/bin2
workd=$datad/chrs/cacao/rnas/velw
velfad=$datad/chrs/cacao/rnas/velw/vel4cfa
rund=$workd/vel$dv

## ** FIXME ?? : velvet wants read len added to inslen, for cgb 106b pairs, full inslen = 412

vopts=""; oopts="";
#v43,4 oopts="-min_pair_count 4 -min_trans_lgth 100 -ins_length 410 -ins_length2 550"
# vopts="-ins_length 300 -ins_length2 410 -ins_length3 550"
# oopts="-min_pair_count 4 -min_trans_lgth 100 -ins_length 300 -ins_length2 410 -ins_length3 550"
#v12, per trsoap 
vopts="-ins_length 200 -ins_length2 200 -ins_length3 450"
oopts="-min_pair_count 4 -min_trans_lgth 180 -ins_length 200 -ins_length2 200 -ins_length3 450"

#v46 kset="89 49 29 25 23 21"
kset="91 51 45 35 25 21"
#soap kset="89 49 35 29 25 23"

echo "START `date` " 

cd $workd
mkdir $rund
cd $rund/

#..... run loop for vel kmer steps
shopt -s nullglob

# check existence; premake seqdir on 1cpu
kseqdir=vel$dv${subset}_seq
if [ ! -f $kseqdir/Sequences ]; then
$velbin/velveth $kseqdir 27 -fasta.gz \
    -shortPaired  $velfad/$datasub.*.sifa2.gz  \
    -shortPaired2 $velfad/$datasub.*.fa2.gz  \
    -shortPaired3 $velfad/$datasub.*.lifa2.gz  \
    -noHash
fi

# v12 data: only fa2 paired rnaseq
# v12: no:  -long $datasub.longfa.gz

for k in $kset;  do { 
  ksubdir=vel$dv${subset}_$k
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

