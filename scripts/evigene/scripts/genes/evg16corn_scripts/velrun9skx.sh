#!/bin/bash
## env name=itick inpe=sraf/allpe.fa2 datad=`pwd` qsub -q normal velrun2.sh
#PBS -N velrun
#PBS -A ind114
#PBS -l nodes=1:ppn=32,walltime=39:55:00
#... -l vmem=96gb,nodes=1:ppn=16,walltime=19:55:00
#PBS -V

# 9skx = try extend stranded using 9ssk long,partialaa seqs + new shortpaired 
# need velbin5 (but kmax = 31!) .. add velbin6, kmax = 151 + longseq **
## test kmer>99 for rdsize=151, * really odd result w/ -strand_specific, 1000 longest aa at kmer 80..95
# 9sk= sl with lower rd count opts, try to get -ss to finish transcripts (many partial3 in -ss)
# 9sl= kmer>99 to rdsize=151, and -strand_specific
# 9sg= 	corn jgi ss data, 151b reads, 
# velvh: -strand_specific	: for strand specific transcriptome sequencing data (default: off)
## 9ag = 9af with kmer shuffle, high only to k43? 

INSIZE=250
STRANDS=1
dv=9skx

if [ "X" = "X$ncpu" ]; then ncpu=8; fi
if [ "X" = "X$datad" ]; then echo "err missing datad=path/to/data"; exit -1; fi
if [ "X" = "X$name" ]; then name=ztick; fi
# if [ "X" != "X$stranded" ]; then STRANDS=1; fi

if [ "X" = "X$inlong" ]; then echo "ERR miss inlong=sraf/longseq.tr"; exit -1; fi
if [ "X" = "X$inpe" ]; then echo "ERR miss inpe=sraf/allpe_[12].fa"; exit -1; fi
inpe=`echo $inpe | sed "s,sraf,$datad/sraf,g;"`
inlong=`echo $inlong | sed "s,sraf,$datad/sraf,g;"`

rund=$datad/velv$name$dv
noasesfork=1
##......................

export OMP_NUM_THREADS=$ncpu

#bin6: longseq + k<=151 or 99?
velbin6=$HOME/bio/velvet1210/bin6
#bin5: longseq>32k + k<=31 ** kmax too low here, boost to 61 **
velbin5=$HOME/bio/velvet1210/bin5
#bin4: k<=151
velbin4=$HOME/bio/velvet1210/bin4
#bin2: k<=99
velbin2=$HOME/bio/velvet1210/bin2
#bin1: k<=61
velbin1=$HOME/bio/velvet1210/bin

velbin=$velbin6; 

#kset4="145 135 125 115 105 95 85 75 65 55 45 35"
kset="93 83 73 63 57 47 37 27"
#kset2="87 77 67 59 53 43 35"

## for norm, ? -cov_cutoff 3 or -cov_cutoff 4, -min_pair_count 3, lower edgeFrac 0.05 ?
vopth=""
if [ $STRANDS -gt 0 ]; then vopth="-strand_specific"; fi

## sd was 50, let velv calc w/o -ins ?
iopts="-ins_length $INSIZE -ins_length_sd 30"
vopts="-conserveLong yes -max_gap_count 5 -min_pair_count 2 $iopts" 
#n#oopts="-scaffolding yes -cov_cutoff 2 -min_pair_count 2 -edgeFractionCutoff 0.03 -min_trans_lgth 200 $iopts"
#skx: merge long+short extend
oopts="-scaffolding yes -merge yes -min_pair_count 2 -edgeFractionCutoff 0.03 -min_trans_lgth 200 $iopts"

echo "START `date` " 
cd $datad
mkdir $rund
cd $rund/

#.... merge long only, kmer ineffective
# $velbin/velveth $kseqdir 27 $vopth -fmtAuto -long $intr
# $velbin/velvetg $ksubdir -conserveLong yes  -read_trkg yes 
# $velbin/oases $ksubdir -scaffolding yes -merge yes

#..... run loop for vel kmer steps
shopt -s nullglob

kseqdir=vel${dv}_seq
if [ ! -f $kseqdir/Sequences ]; then
 $velbin/velveth $kseqdir 31 $vopth -noHash -fmtAuto -long $inlong -shortPaired -separate $inpe
fi

for k in $kset;  do { 
  ksubdir=vel${dv}_$k
  echo "#.. start velrun $ksubdir : `date`"
## longseq velbin fix
  velbin=$velbin6; 
  # if [ $k -le 99 ]; then velbin=$velbin2; fi
  # if [ $k -le 61 ]; then velbin=$velbin5; fi
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
## longseq velbin fix
    velbin=$velbin6;
    # if [ $k -le 99 ]; then velbin=$velbin2; fi
    # if [ $k -le 61 ]; then velbin=$velbin5; fi
    $velbin/oases   $ksubdir $oopts 
    i=$(( $i + 1 ))
    #D# if [ $i -ge $noasesfork ]; then wait; i=0; fi
  fi
} done

wait

echo "DONE `date` " 

