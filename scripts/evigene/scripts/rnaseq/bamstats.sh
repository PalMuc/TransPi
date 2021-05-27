#!/bin/tcsh

set rnas=/bio/bio-grid/mb/rnaseq/
set bstat=aphidrs_bam.stats
touch $bstat

foreach bam ( bams/aphidrs_*.bam )
 if ( $bam =~ *_mapt.bam ) continue
 echo "==== $bam =======================" >> $bstat
 samtools view $bam | $rnas/rgsoft/samstats.pl >> $bstat
end
