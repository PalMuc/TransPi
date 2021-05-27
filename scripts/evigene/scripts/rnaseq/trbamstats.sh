#!/bin/tcsh
# trbamstats.sh

set i=0
set maxi=99999
set evigene=/bio/bio-grid/mb/evigene/
set outf="trstats.tab"
touch $outf

## skip 0,1 read entries to speed this up, also drop alttr (only t1)? # arg but need to do this for all .bam 
# samtools idxstats SRR064409-ap2pub8tr.bam | cut -f1,3 | grep ' 0' | wc = 9320
# samtools idxstats SRR064409-ap2pub8tr.bam | cut -f1,3 | grep ' 1$' | wc = 2145

set eglist=`samtools idxstats SRR064409-ap2pub8tr.bam | cut -f1 | grep 't1$'`
## | grep '^EG' `

foreach eg ( $eglist )
  @ i = $i + 1
 if ($i > $maxi ) break
 touch $eg.sam
 foreach bm (*tr.bam)
   samtools view $bm  $eg >> $eg.sam
 end
 cat $eg.sam | env name=$eg chrlen=1 keepids=0 table=$i $evigene/scripts/rnaseq/samstats.pl >> $outf
 /bin/rm $eg.sam
end

