#!/bin/tcsh

foreach bam ( $* )
set samset=`basename $bam .bam`
samtools flagstat $bam > $samset.stats
echo -n 'stranded reads: ' >> $samset.stats
samtools view $bam | grep -c 'XS:A:' >> $samset.stats
end
