#!/bin/tcsh

set rnas=/bio/bio-grid/mb/rnaseq/
set scaf=Scaffold2
set sloc="Scaffold2:1-2384549"

#set bamset=aphidrs_*.bam
set bamset="$*"

## for rs_ needs index bam.bai
foreach bam ( $bamset )
 set nam=`basename $bam .bam`
 if ( -f $nam.$scaf.sam ) continue
 echo samtools view $bam $sloc 
 samtools view $bam $sloc > $nam.$scaf.sam
end

exit

## for pe_
foreach bam ( aphidpe_*.bam )
 set nam=`basename $bam .bam`
 if ( -f $nam.$scaf.sam ) continue
 echo "==== $bam =======================" 
 samtools view $bam | grep "$scaf	" > $nam.$scaf.sam
end

