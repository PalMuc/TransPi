#!/bin/bash
## env fasta=xxx datad=`pwd` qsub -q shared bowbuild.sh
#PBS -N bowb
#PBS -A ind114
#PBS -l nodes=1:ppn=2,walltime=3:55:00
#PBS -o bowb.$$.out
#PBS -e bowb.$$.err
#PBS -V

bowd=$HOME/bio/bowtie2

if [ "X" = "X$datad" ]; then
  echo "ERR: missing env datad=/path/to/fasta"; exit -1
fi
if [ "X" = "X$fasta" ]; then
  echo "ERR: missing env fasta=mydata.fasta"; exit -1
fi

cd $datad
mkdir bow

$bowd/bowtie2-build $fasta bow/$fasta 
# $bowd/bowtie2-build daphmag5xall13f_okall.tr bow/daphmag5xall13f_okall.tr >& bow/log.bbx

