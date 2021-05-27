#!/bin/tcsh

set rnas=/bio/bio-grid/mb/rnaseq
## sort big file.tmpfa
setenv TMPDIR /var/tmp

foreach bam ( ../bam2/*.bam )
set rsa=`basename $bam .bam`
samtools view -f 4 $bam | $rnas/rgsoft/sam2velv.pl -debug -name $rsa.nomap
end
# $rnas/rgsoft/sam2velv.pl -debug -nodup=1 -name $rsa.nomap

