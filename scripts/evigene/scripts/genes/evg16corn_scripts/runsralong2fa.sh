#!/bin/bash
## env sra=sraf/SRX*/SRR*/SRR*.sra datad=`pwd` qsub -q normal runsra2fa.sh
#PBS -N sra2fa
#PBS -A ind114
#PBS -l nodes=1:ppn=16,walltime=9:55:00
#PBS -V

## PacBio fix: $fqdump --minReadLen 400  -O spotfa  --fasta 0 SRX1472849/SRR2983686/SRR2983686.sra

fqdump=$HOME/bio/sratoolkit/fastq-dump

if [ "X" = "X$ncpu" ]; then ncpu=16; fi
if [ "X" = "X$datad" ]; then echo "ERR: datad=/path/to/data"; exit -1; fi
if [ "X" = "X$sra" ]; then echo "ERR: sra=what?"; exit -1; fi
if [ "X" = "X$minlen" ]; then minlen=400; fi

cd $datad
mkdir spotfa
i=0; for sr in $sra; do { 
  # skip if done
  sfa=`basename $sr .sra`
  if [ -f spotfa/$sfa.fasta ]; then echo "already have spotfa/$sfa.fasta"; continue; fi
  echo $fqdump --minReadLen $minlen -O spotfa --fasta 0 $sr
  $fqdump --minReadLen $minlen -O spotfa --fasta 0 $sr  &
  i=$(( $i + 1 )); if [ $i -ge $ncpu ]; then wait; i=0; fi 
} done

wait

## NO pairfa here..
