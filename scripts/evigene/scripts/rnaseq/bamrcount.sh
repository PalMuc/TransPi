#!/bin/tcsh

foreach bam ( bams/*.bam )
 samtools idxstats $bam > $bam.count
end
