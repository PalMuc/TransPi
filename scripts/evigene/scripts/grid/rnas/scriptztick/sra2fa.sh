#!/bin/bash

# for diginorm: fasta pairs in same .fa, qual-filter
# process pe and sr separately
# use --split-spot for 2 reads/spot ?
# -E|--qual-filter ??  drop NNN-tail reads << DONT for diginorm pe that wants pairs always; no drop1/2
## filter only if -nosplitspot, then split pairs for dignorm

for sra in SRX*/SRR*/SRR*.sra; do {
  echo $HOME/bio/sratoolkit/fastq-dump --split-spot --fasta 0 $sra
  $HOME/bio/sratoolkit/fastq-dump --qual-filter --split-spot --split-3 --fasta 0 $sra
  # $HOME/bio/sratoolkit/fastq-dump --split-spot --fasta 0 $sra
} done

