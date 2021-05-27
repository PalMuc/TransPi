#!/bin/tcsh

set cacfad=/export/disk2/work/cacao2d/est2
set gffdir=estgff_mars
set fagff2fa="./fagff2sam.pl -out=fa -skip='exon,intron,CDS'"

foreach gf ($gffdir/reads.*-mars11.gmap.gff.gz)
  set nam=`basename $gf  -mars11.gmap.gff.gz`
  echo "$fagff2fa -fa=$cacfad/$nam.fa.gz -gff=$gf TO  estfa/$nam.longfa"
  $fagff2fa -fa=$cacfad/$nam.fa.gz -gff=$gf > estfa/$nam.longfa
end

