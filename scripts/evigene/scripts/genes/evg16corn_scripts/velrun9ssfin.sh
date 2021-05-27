#!/bin/bash
## env name=itick inpe=sraf/allpe.fa2 datad=`pwd` qsub -q normal velrun2.sh
#PBS -N velrun
#PBS -A ind114
#PBS -l nodes=1:ppn=32,walltime=39:55:00
#... -l vmem=96gb,nodes=1:ppn=16,walltime=19:55:00
#PBS -V

# velrun9ssfin.sh = oases finish only, velvetg fail outamem, but for 4 ks
# 9sg= 	corn jgi ss data, 151b reads, 
# velvh: -strand_specific	: for strand specific transcriptome sequencing data (default: off)

## 9ag = 9af with kmer shuffle, high only to k43? 
## 9tb: redo 9t, dapx adult9 set best asm, but modify kmers, fail below k47, add above k73
##  .. mod INSIZE, old/bad 150, new 270, split diff at 200? .. same as v9b
## 9t: redo successful 8t, remove/reduce velg covcut (=1 or none)
## .. also update velvet127s to velvet1210 

INSIZE=250
STRANDS=1
dv=9sg
if [ "X" = "X$ncpu" ]; then ncpu=8; fi
if [ "X" = "X$datad" ]; then echo "err missing datad=path/to/data"; exit -1; fi
if [ "X" = "X$name" ]; then name=ztick; fi
rund=$datad/velv$name$dv
if [ "X" = "X$inpe" ]; then echo "ERR miss inpe=sraf/allpe*.fa2"; exit -1; fi
inpe=`echo $inpe | sed "s,sraf,$datad/sraf,g;"`

noasesfork=1
##......................

export OMP_NUM_THREADS=$ncpu

#k>61# 
velbin2=$HOME/bio/velvet1210/bin2
#k<=61# 
velbin1=$HOME/bio/velvet1210/bin
velbin=$velbin2; 

kset="97 93 83 73 63 57 47 37 33"
#kset1="93 83 73 63 57 47 37"
#kset2="87 77 67 59 53 43 35"

## for norm, ? -cov_cutoff 3 or -cov_cutoff 4, -min_pair_count 3, lower edgeFrac 0.05 ?
vopth=""
if [ $STRANDS ]; then vopth="-strand_specific"; fi
iopts="-ins_length $INSIZE -ins_length_sd 50"
vopts="$iopts -max_gap_count 5 " 
oopts="-scaffolding yes -min_pair_count 3 -edgeFractionCutoff 0.03 -min_trans_lgth 200 $iopts"

echo "START `date` " 
cd $datad
mkdir $rund
cd $rund/

#..... run loop for vel kmer steps
shopt -s nullglob

kseqdir=vel${dv}_seq
if [ ! -f $kseqdir/Sequences ]; then
 $velbin/velveth $kseqdir 27 $vopth -noHash -shortPaired -fmtAuto -separate $inpe
fi

## fin: DONE HERE..
#d for k in $kset;  do { 
#d   ksubdir=vel${dv}_$k
#d   echo "#.. start velrun $ksubdir : `date`"
#d   velbin=$velbin2; if [ $k -le 61 ]; then velbin=$velbin1; fi
#d   mkdir $ksubdir
#d   ln -s ../$kseqdir/Sequences $ksubdir/
#d   $velbin/velveth $ksubdir $k $vopth -reuse_Sequences
#d   $velbin/velvetg $ksubdir $vopts -read_trkg yes 
#d 
#d   # $velbin/oases   $ksubdir $oopts 
#d   # /bin/rm $ksubdir/{Graph2,LastGraph,PreGraph,Roadmaps,Sequences}
#d   echo "#.. end velvetg $ksubdir : `date`"
#d } done
#d ## end loop 1
 
## loop 2 oases, 1cpu; can we fork some at same time? enough mem?
export OMP_NUM_THREADS=2
export OMP_THREAD_LIMIT=2
i=0; 
for k in $kset;  do {
  ksubdir=vel${dv}_$k
  if [ -f $ksubdir/transcripts.fa ]; then continue; fi
  if [ -f $ksubdir/contigs.fa ]; then
    ## DONT fork here save mem ..
    velbin=$velbin2; if [ $k -le 61 ]; then velbin=$velbin1; fi
    $velbin/oases   $ksubdir $oopts 
    i=$(( $i + 1 ))
    #D# if [ $i -ge $noasesfork ]; then wait; i=0; fi
  fi
} done

wait

echo "DONE `date` " 

