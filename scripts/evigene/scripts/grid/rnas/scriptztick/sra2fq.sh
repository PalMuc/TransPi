#!/bin/bash

# use --split-3, --split-spot for 2 reads/spot ; soaptr needs 1.fq, 2.fq; velv  either
# -E|--qual-filter ??  drop NNN-tail reads; do for --fasta 0 output

sraset=$*
# sraset=SRX*/SRR*/SRR*.sra

for sra in $sraset; do {
  echo $HOME/bio/sratoolkit/fastq-dump --split-spot $sra
  $HOME/bio/sratoolkit/fastq-dump --split-spot --split-3 $sra
  # $HOME/bio/sratoolkit/fastq-dump --split-spot  $sra
} done

## --gzip here too slow; use gzip --fast after dump.
