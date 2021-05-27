#!/bin/bash
## env name=itick inpe=sraf/allpe.fa2 datad=`pwd` qsub -q normal velrun2.sh
#PBS -N velrun
#PBS -A ind114
#PBS -l nodes=1:ppn=32,walltime=39:55:00
#PBS -V

# velmerge.sh == retry velv/oases -merge of all velv*.tr 
# ** velbin5 for  LONGSEQ merge only (kmax 31)
# 9sm= mouse trin data 75b reads, strspp, insert=300?
# 9am= mouse, OFF strspp * BUG, not off, 9an == OFF
# 9sg= 	corn jgi ss data, 151b reads, 
# velvh: -strand_specific	: for strand specific transcriptome sequencing data (default: off)

INSIZE=290
STRANDS=0
dv=merge

if [ "X" = "X$ncpu" ]; then ncpu=8; fi
if [ "X" = "X$datad" ]; then echo "err missing datad=path/to/data"; exit -1; fi
if [ "X" = "X$name" ]; then name=ztick; fi
if [ "X" != "X$stranded" ]; then STRANDS=1; fi

if [ "X" = "X$intr" ]; then echo "ERR miss intr=invelvall.tr"; exit -1; fi
# if [ "X" = "X$inpe" ]; then echo "ERR miss inpe=sraf/allpe*.fa2"; exit -1; fi
# inpe=`echo $inpe | sed "s,sraf,$datad/sraf,g;"`

rund=$datad/velv$name$dv
noasesfork=1
##......................

export OMP_NUM_THREADS=$ncpu
# try ncpu: export OMP_NUM_THREADS=2
# export OMP_THREAD_LIMIT=2

#longseq, k<=31# 
velbin=$HOME/bio/velvet1210/bin5

# kset="73 63 57 53 47 43 37 33 27"

## for norm, ? -cov_cutoff 3 or -cov_cutoff 4, -min_pair_count 3, lower edgeFrac 0.05 ?
vopth=""
## DANG BUG if $STRANDS always true for 0/1 in bash **
if [ $STRANDS -eq 1 ]; then vopth="-strand_specific"; fi
# iopts="-ins_length $INSIZE -ins_length_sd 50"
# vopts="$iopts -max_gap_count 5 " 
# oopts="-scaffolding yes -min_pair_count 3 -edgeFractionCutoff 0.03 -min_trans_lgth 200 $iopts"

# merge
kset="27"
iopts=""
vopts="-conserveLong yes -read_trkg yes $iopts "
oopts="-scaffolding yes -min_trans_lgth 200 $iopts"
#.......

echo "START `date` " 
cd $datad
#no# mkdir $rund
#no# cd $rund/

shopt -s nullglob

# kseqdir=vel${dv}_seq
# ksubdir=vel${dv}
ksubdir=$rund
kseqdir=$ksubdir
mkdir $ksubdir

# for k in $kset;  do { ksubdir=vel${dv}_$k;  mkdir $ksubdir; ln -s ../$kseqdir/Sequences $ksubdir/; }

$velbin/velveth $kseqdir 27 $vopth -fmtAuto -long $intr

$velbin/velvetg $ksubdir -conserveLong yes  -read_trkg yes 

$velbin/oases $ksubdir -scaffolding yes -merge yes

echo "DONE `date` " 

