#!/bin/tcsh

set rnas=/bio/bio-grid/mb/rnaseq
set workd=/export/udisk2/work/aphid/
set filtfile=$workd/misc/repeat_rrna.gff
## sort big file.tmpfa
setenv TMPDIR /var/tmp

##foreach bam ( $workd/rnas/bams/aphidpe_*.bam )
foreach bam ( $workd/rnas/bams/aphidrs_*.bam)
 set rsa=`basename $bam .bam`
 echo "#======= $rsa.mapclean ============="
 samtools view -F 4 $bam | \
 $rnas/rgsoft/samfiltergff.pl -expandover=50 -in stdin -over $filtfile  | \
 $rnas/rgsoft/sam2velv.pl  -debug -nodup=4 -name $rsa.mapclean
end

