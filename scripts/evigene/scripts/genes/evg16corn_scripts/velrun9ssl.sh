#!/bin/bash
## env name=itick inpe=sraf/allpe.fa2 datad=`pwd` qsub -q normal velrun2.sh
#PBS -N velrun
#PBS -A ind114
#PBS -l nodes=1:ppn=32,walltime=39:55:00
#... -l vmem=96gb,nodes=1:ppn=16,walltime=19:55:00
#PBS -V

## test kmer>99 for rdsize=151, * really odd result w/ -strand_specific, 1000 longest aa at kmer 80..95
# 9sl= kmer>99 to rdsize=151, and -strand_specific
# 9sg= 	corn jgi ss data, 151b reads, 
# velvh: -strand_specific	: for strand specific transcriptome sequencing data (default: off)

## 9ag = 9af with kmer shuffle, high only to k43? 
## 9tb: redo 9t, dapx adult9 set best asm, but modify kmers, fail below k47, add above k73
##  .. mod INSIZE, old/bad 150, new 270, split diff at 200? .. same as v9b
## 9t: redo successful 8t, remove/reduce velg covcut (=1 or none)
## .. also update velvet127s to velvet1210 

INSIZE=250
STRANDS=1
dv=9sl
if [ "X" = "X$ncpu" ]; then ncpu=8; fi
if [ "X" = "X$datad" ]; then echo "err missing datad=path/to/data"; exit -1; fi
if [ "X" = "X$name" ]; then name=ztick; fi
rund=$datad/velv$name$dv
if [ "X" = "X$inpe" ]; then echo "ERR miss inpe=sraf/allpe*.fa2"; exit -1; fi
inpe=`echo $inpe | sed "s,sraf,$datad/sraf,g;"`

noasesfork=1
##......................

export OMP_NUM_THREADS=$ncpu

#k>99 bin4
velbin4=$HOME/bio/velvet1210/bin4
#k>61# 
velbin2=$HOME/bio/velvet1210/bin2
#k<=61# 
velbin1=$HOME/bio/velvet1210/bin

velbin=$velbin4; 

kset="145 135 125 115 105 95 85 75 65 55 45 35"
## k151 = read size is no good
#kset1="93 83 73 63 57 47 37"
#kset2="87 77 67 59 53 43 35"

## for norm, ? -cov_cutoff 3 or -cov_cutoff 4, -min_pair_count 3, lower edgeFrac 0.05 ?
vopth=""
if [ $STRANDS ]; then vopth="-strand_specific"; fi
## sd was 50, let velv calc w/o -ins ?
iopts="-ins_length $INSIZE -ins_length_sd 20"
vopts="$iopts -max_gap_count 5 " 
oopts="-scaffolding yes -min_pair_count 3 -edgeFractionCutoff 0.03 -min_trans_lgth 200 $iopts"
#n#oopts="-scaffolding yes -min_pair_count 3 -min_trans_lgth 200 $iopts"

echo "START `date` " 
cd $datad
mkdir $rund
cd $rund/

#..... run loop for vel kmer steps
shopt -s nullglob

kseqdir=vel${dv}_seq
if [ ! -f $kseqdir/Sequences ]; then
 $velbin/velveth $kseqdir 37 $vopth -noHash -shortPaired -fmtAuto -separate $inpe
fi

for k in $kset;  do { 
  ksubdir=vel${dv}_$k
  echo "#.. start velrun $ksubdir : `date`"
  velbin=$velbin4; 
  if [ $k -le 99 ]; then velbin=$velbin2; fi
  if [ $k -le 61 ]; then velbin=$velbin1; fi
  mkdir $ksubdir
  ln -s ../$kseqdir/Sequences $ksubdir/
  $velbin/velveth $ksubdir $k $vopth -reuse_Sequences
  $velbin/velvetg $ksubdir $vopts -read_trkg yes 

  # $velbin/oases   $ksubdir $oopts 
  # /bin/rm $ksubdir/{Graph2,LastGraph,PreGraph,Roadmaps,Sequences}
  echo "#.. end velvetg $ksubdir : `date`"
} done
## end loop 1

## loop 2 oases, 1cpu; can we fork some at same time? enough mem?
export OMP_NUM_THREADS=2
export OMP_THREAD_LIMIT=2
i=0; 
for k in $kset;  do {
  ksubdir=vel${dv}_$k
  if [ -f $ksubdir/transcripts.fa ]; then continue; fi
  if [ -f $ksubdir/contigs.fa ]; then
    ## DONT fork here save mem ..
    velbin=$velbin4;
    if [ $k -le 99 ]; then velbin=$velbin2; fi
    if [ $k -le 61 ]; then velbin=$velbin1; fi
    $velbin/oases   $ksubdir $oopts 
    i=$(( $i + 1 ))
    #D# if [ $i -ge $noasesfork ]; then wait; i=0; fi
  fi
} done

wait

echo "DONE `date` " 

