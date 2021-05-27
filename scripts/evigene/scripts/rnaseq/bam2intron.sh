#!/bin/tcsh
# samtools -F 0x104 = drop unmapped(4), duplicates(x100)

#set sbin=/bio/bio-grid/mb/rnaseq/rgsoft
set evigene=/bio/bio-grid/mb/evigene/
set sbin=$evigene/scripts/rnaseq/
setenv TMPDIR /var/tmp

# mkdir introns; cd introns/

foreach bam ( ../bams/*.bam )
  set nam=`basename $bam .bam`
  samtools view -F 0x104 $bam | grep 'XS:A:' | $sbin/samintrons.pl -source rs$nam -sort > $nam.intron.gff
end

