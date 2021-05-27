#! /bin/bash
### env trfasta=genes.fasta gsnapdir="gsodm13xpa gsodm13xtr" datad=`pwd` qsub -q normal gsnpbam2xp3n.sh
## daphmag/rnas/scriptg/gsnpbam2xp3n.sh
#PBS -N gsnap2rnax
#PBS -l nodes=1:ppn=16,walltime=19:55:00
#PBS -V

ncpu=16
ncpu1=4
# gsnpbam2xp3n.sh : rnaexpress runs #  gsnapdir="gsodm13xpa gsodm13xtr" ..
## NO good: rnaexpress ... in1.bam in2.bam : only handles 1 dammit.. << NO, uses "in1.bam,in2.bam,.."
## .. still no good: fails w/ other error for 2+ bams ;; 
# problem may be paired_ = ERROR: Input BAM file contains no valid alignments.

xprbig="concordant_mult"
xpradd="concordant_uniq"
## cc mult,uniq include ~95% mapped reads
# xprset="concordant_mult concordant_uniq paired_mult paired_uniq_long paired_uniq_inv"
#NO xpradd="concordant_uniq paired_mult paired_uniq_long paired_uniq_inv"

if [ "X" = "X$gsnapdir" ]; then echo "ERR: miss gsnapdir=gsnap dir list "; exit -1; fi
if [ "X" = "X$datad" ]; then echo "ERR: miss datad=path/to/gsnapdirs"; exit -1; fi
if [ "X" = "X$trfasta" ]; then echo "ERR: missing env trfasta=mydata.tr"; exit -1; fi

biobin=$HOME/bio/bin
cd $datad
tnam=`basename $trfasta .mrna | sed 's/\.tr$//; s/\.fasta//; s/\.fa$//; s/_ok.*//;'`
# datad=xxx/gsodm13sd  gdir=gsodm13cxpa_dn1/Dman_80_GCCAAT_L006_1-gsnap-concordant_mult.bam

echo "START rnax: `date`"
icpu=0;
for gdir in $gsnapdir ; do { 
  cd $datad
  glist=`ls $gdir/*$xprbig.bam`; 
  for bam1 in $glist; do { 
    bnam=`echo $bam1| sed "s/.$xprbig.bam//;"`
    bams=$bam1; for bsuf in $xpradd; do { bams="$bams,$bnam-$bsuf.bam"; } done
    bnamo=`basename $bnam | sed 's/.gsnap//; s/_1$//;'`
    outrx=rx$tnam-$bnamo
    ## no.need.app.does# echo mkdir $outrx

    echo $biobin/rnaexpress -p $ncpu1 -o $outrx --no-update-check $trfasta $bams 
    $biobin/rnaexpress -p $ncpu1 -o $outrx --no-update-check $trfasta $bams  &
   
    ## innerloop per readfile
    icpu=$(( $icpu + $ncpu1 )); 
    if [ $icpu -ge $ncpu ]; then wait; icpu=0; fi

   } done
   # glist

  cd $datad
} done
# gsnapdir
wait

echo "DONE rnax: `date`"
