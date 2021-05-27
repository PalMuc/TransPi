#!/bin/tcsh

set sbin=/bio/bio-grid/mb/rnaseq/rgsoft
# set sbin=../scripts

# foreach samz (sams/cacao01.gsnap.sam.gz)
# foreach samz (sams/cacao0[2,4-9].gsnap.sam.gz sams/cacao10.gsnap.sam.gz sams/cacao03.gsnap.sam.gz)
foreach samz (cleans/aphid*.clean.sam.gz)
  set nam=`basename $samz .sam.gz`
  echo "# $nam.intron1.gff"
  gzgrep 'XS:A:' $samz | $sbin/samintrons.pl -sort > cleans/$nam.intron1.gff
end
