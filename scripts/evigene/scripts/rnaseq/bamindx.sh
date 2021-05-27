#!/bin/tcsh

# foreach bam (004_R1.bam VV*.bam)
foreach bam (*.bam)
 if( -f $bam.bai ) continue
 echo samtools index $bam
 samtools index $bam
 if( -f $bam.count) continue
 samtools idxstats $bam > $bam.count
end

