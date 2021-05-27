#! /bin/bash
##  export subset=subset.Scaffold16; gg-subset.sh > log.$subset 2>&1
## call from qsub list.sh

velbin=$HOME/bio/velvet/bin
shopt -s nullglob
# need * in filenames and nullglob to work for null files
# .fa2, .longfa2 pairs as >id/1 aaa; >id/2 ccc;

# Set these in caller:  OPENMP variant: set number threads/cpus
# http://listserver.ebi.ac.uk/pipermail/velvet-users/2011-April/001325.html
#     $OMP_NUM_THREADS $OMP_THREAD_LIMIT
# export OMP_NUM_THREADS=number 
# export OMP_THREAD_LIMIT=number 


# kmer=21
# kmer="19,29,2" # does it work for oases/rna?
kmer=27   # faster, better for rnas? for high-cover tr; save UnusedReads to run w/ low kmer=17

echo "#.. start rnavelv $subset : `date`"
##bad.fn# if [ -f $subset.fa.gz ]; then
if [ 1 == 1 ]; then

  # velh opts: kmer= min,max,step instead of one value .. does it work for rna?
  $velbin/velveth vel$subset $kmer  -fasta.gz \
     -shortPaired ${subset}*.fa2.gz -short ${subset}*.fa1.gz ${subset}*.fa.gz \
     -longPaired ${subset}*.longfa2.gz -long ${subset}*.longfa1.gz ${subset}*.longfa.gz 

  $velbin/velvetg vel$subset -read_trkg yes 
  
  # oases opt defaults: -min_pair_count 4  -min_trans_lgth $kmer  -conserveLong no
  # use -conserveLong yes for long ESTs??
  $velbin/oases vel$subset -min_trans_lgth 40 \
    -ins_length 200 -ins_length_long 400 
    
  /bin/rm vel$subset/{Graph2,LastGraph,PreGraph,Roadmaps,Sequences}
fi

if [ -f vel$subset/transcripts.fa ]; then
  echo "#.. OK "
else
  echo "#.. ERROR: no vel/transcripts.fa"
fi
echo "#.. end rnavelv : `date`"

## sample caller
# #----------------------------------------------------------------
# #! /bin/bash -l
# ##  qsub -q normal subset2.sh
# #PBS -N velsub2
# #PBS -A TG-MCB100147
# #PBS -l ncpus=8,mem=127gb,walltime=23:55:00
# #PBS -o velsub.$$.out
# #PBS -e velsub.$$.err
# #PBS -V
# 
# ncpu=7  # at least cpu-1 for thread max; also >>mem for >cpu
# export OMP_NUM_THREADS=$ncpu 
# export OMP_THREAD_LIMIT=$ncpu 
# workd=$HOME/scratch/chrs/aphid2
# velbin=$HOME/bio/velvet/bin
# 
# cd $workd/rnas/subset2/
# echo "#.. start subsetvel : `date`"
# export subset=subset.mid0 ; $workd/rnas/gg-subset.sh > log.$subset 2>&1
# export subset=subset.mid1 ; $workd/rnas/gg-subset.sh > log.$subset 2>&1
# ...
# echo "#.. end subsetvel : `date`"
