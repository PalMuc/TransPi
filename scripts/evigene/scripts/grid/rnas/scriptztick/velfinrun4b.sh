#!/bin/bash
## env ksetfin="55 49 43" workd=`pwd` qsub -q batch velfinrun4a.sh
#PBS -N velrun 
#PBS -A ind114
#PBS -l nodes=1:ppn=16,walltime=39:55:00
#...... vmem=124gb,nodes=1:ppn=2,walltime=39:55:00
#PBS -o velrun.$$.out
#PBS -e velrun.$$.err
#PBS -V

## ../velfinrun4b.sh : gordo.sdsc, try kset=vel4_55 vel4_49 vel4_43 with fork/wait for 3 at one run
ksetfin="55 49 43"
# if [ "X" = "X$ksetfin" ]; then echo "err missing ksetfin="55 49"; exit -1; fi

## ../velfinrun4a.sh = finish oases that died outa time.. test opts to speed up; 1 cpu, <127GB shared?
## at mason.iu this is taking forever... not finishing; oases is logjam;
## at gordon.sdsc, try 2-4 ksubdir for 16 cpu, 64GB ? Oases is using only 10G-20G at iu..
 
## usage: env ksubdir==vel4_49 workd=path/above/ksubdir qsub ..

# runversion
# zebra tick, 101bp pe illum, 100b inner insert, approx 2-sex x 100M reads
## velrun4yu.sh : overlay run dv=4 w/ alternate kset, lower kmers for bigmem

dv=4
ncpu=1
# export OMP_NUM_THREADS=$ncpu
# export OMP_THREAD_LIMIT=$ncpu
## bin1 for simple, kmer<60 .. bin2/ for kmer<106
# velbin=$HOME/bio/velvet127s/bin2
velbin=$HOME/bio/velvet127s/bin

if [ "X" = "X$workd" ]; then echo "err missing workd=path/to/data"; exit -1; fi
## if [ "X" = "X$ksubdir" ]; then echo "err missing ksubdir=xxx_nn"; exit -1; fi
## rund=$workd/velv$name$dv

## ave ins 300 (100rd x 2 + 100 inside)
vopts="-ins_length 300 "
#oopt0="-min_pair_count 2 -scaffolding yes -min_trans_lgth 180 -ins_length 300"
#oopt1="-paired_cutoff 0.2 -edgeFractionCutoff 0.1 -scaffolding yes -min_trans_lgth 180 -ins_length 300"
oopts="-paired_cutoff 0.2 -edgeFractionCutoff 0.1 -scaffolding yes -min_trans_lgth 200 -ins_length 300"
## boosting: edgeFractionCutoff paired_cutoff min_pair_count to remove data/reduce time?
#x1: kset="95 85 75 65 61 55 51 45 41 35 29 25"
#y2: kset="49 43 33 27 37 47 57 67 77"

echo "START `date` " 
cd $workd
# mkdir $rund; cd $rund/

#..... run loop for vel kmer steps
shopt -s nullglob

#. kseqdir=vel${dv}_seq
#. if [ ! -f $kseqdir/Sequences ]; then
#.  $velbin/velveth $kseqdir 27 -noHash -fasta.gz -shortPaired -interleaved $inpe  
#. fi

## add kset loop for gordo 16cpu ?
for k in $ksetfin; do {
  ksubdir=vel${dv}_$k
  echo "$velbin/oases   $ksubdir $oopts"
  $velbin/oases   $ksubdir $oopts &
} done

wait

#.............
# for k in $kset;  do { 
  #param# ksubdir=vel${dv}_$k
  # echo "#.. start velrun $ksubdir : `date`"
  # mkdir $ksubdir
  # ln -s ../$kseqdir/Sequences $ksubdir/
  # $velbin/velveth $ksubdir $k  -reuse_Sequences
  # $velbin/velvetg $ksubdir $vopts -read_trkg yes 
  # $velbin/oases   $ksubdir $oopts 
  # echo "# /bin/rm $ksubdir/{Graph2,LastGraph,PreGraph,Roadmaps,Sequences} "
  # echo "#.. end velrun $ksubdir : `date`"
# } done

echo "DONE `date` " 

