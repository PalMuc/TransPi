#!/bin/tcsh
# runcuffdiff.sh

setenv LD_LIBRARY_PATH /usr/lib:/usr/local/gcc4/lib/amd64:/usr/local/gcc4/lib:/usr/local/lib:/usr/sfw/lib
set bindir=/bio/bio-grid/mb/rnaseq/bin

# cd bams/

$bindir/cuffdiff09 -p 2 -o cuffd9rep/ \
--labels 16-1,16-8,16-C,3-1,3-8,3-C,4-1,4-8,4-C \
dmaggenedum.gtf \
16-1-R2.bam,16-1_R3.bam \
16-8-R2.bam,16-8-R3.bam \
16-C-R2.bam,16-C_R3.bam \
3-1_R2.bam \
3-8_R2.bam \
3-C_R1.bam,3-C_R2.bam \
4-1_R2.bam \
4-8_R1.bam,4-8_R2.bam,4-8_R3.bam \
4-C-R2.bam,4-C-R3.bam \

