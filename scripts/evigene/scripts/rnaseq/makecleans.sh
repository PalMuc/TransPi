#!/bin/tcsh
# make clean.sam, rrna.sam with hardclip for cufflinks v0.6
# samtools -F 0x104 = drop unmapped(4), duplicates(x100)

set rnas=/bio/bio-grid/mb/rnaseq
setenv TMPDIR /var/tmp

foreach bam ( ../bam2/*.bam )
set reads=`basename $bam .bam`
samtools view -F 0x104 $bam | \
$rnas/rgsoft/samfiltergff.pl -hardclip -expandover=50 -mate \
    -in stdin -over ../../misc/repeat_rrna.gff \
    -filter $reads.rrna.sam -notfilter $reads.clean.sam \
    >& log.samfilt.$reads

end

