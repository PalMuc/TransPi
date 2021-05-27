#!/bin/bash
# makealt.sh

MINTR=200
altopts="-mincds=50 -minexon=39 -mincoding=40 -minpcds=60 -debug "
#v1.altopts="-minexon=50 -mincod=59 -minpcds=89 "

## 2013sep updates, using gmap loci info, equalgene, ..
#......
## fixme: merge locussametab.pl and diffloci.eqalt2main tabulators to one script
# .. equalgene: need both of alt2main (for diffloci) and all.eqgene (for sameloci)
# $evigene/scripts/equalgene.pl -in $nam.analt.gff -ov $nam.anmain.gff > $nam.alt2main.eqgene &
# $evigene/scripts/equalgene.pl -in $nam.an2.gff -ov $nam.an2.gff > $nam.all.eqgene &
#
# locussametab.pl -eqgene $nam.all.eqgene -idtab $nam.newids -diffloc $nam.diffloci.eqalt2main2 -alntab $nam.alnsense.tab > $nam.sameloci9.tab
# cat $nam.diffloci.eqalt2main2 | egrep ' (Df1|Df2)       ' | sed 's/^/difflocus  /;' | cat - $nam.sameloci9.tab $nam.missid9 $nam.newids | perl locussamediff2pubid.pl > $nam.renewid9
# cat $nam.an2.gff missmain2a.gff | env idtab=$nam.renewid9 nam=${nam}2e ../scripts/splitmainalt.pl
# .. altbest
# set altopts="-noCHANGEALTID -minintron=30 -mincds=33 -minexon=33 -mincoding=40 -minpcds=20 -debug "
# $evigene/scripts/altbest.pl $altopts -main $nam.main.gff -alt $nam.alt.gff -eqtab $nam.equalalt.tab 
#.........

# step1. select subset of all 6 asmrnaest.gff with long enough transcript; use mRNA cxlen= all have
# step2.  annotate introns intr=
# step3.  altbest.pl -main $main.gff -alt allalt.gff -eqtab $main.equalalt.tab
# input exons must be marked with intr=..N1,N2, intron splice ids ; intron chain annotation:
#   $workd/scripts/overlapfilter.perl -intron2splice=error -pass 'exon,intron' -act markid -midtype scoresum \
#   -mark intr -over $workd/intron/intron_good.gff.gz -in genes.gff \

# ** FIXED altbest.pl : reduce for partial prot vs shorter full prot of ~same rna models
# .. problem from pasa_asm > longest-orf-partials not as good as full prot alt models.


## rewrite params ..
# evigene=/bio/bio-grid/mb/evigene/
# workd=/bio/bio-grid/nasv4
# introns=$workd/intron/intron_good.gff.gz
# maingff=../nvit2_evigenes.pub11d.gff.gz
if [ "X" = "X$evigene" ]; then echo "ERR: evigene=what evigene path?"; exit -1; fi
if [ "X" = "X$workd" ]; then echo "ERR: workd=what datapath/?"; exit -1; fi
if [ "X" = "X$genes" ]; then echo "ERR: genes=what genes.gff?"; exit -1; fi
if [ "X" = "X$introns" ]; then echo "WARN: introns=$workd/intron/intron_good.gff.gz ?"; introns=$workd/intron/intron_good.gff.gz; fi
if [ "X" = "X$asmset" ]; then echo "WARN: asmset=*.gff.gz ?"; asmset=*.gff.gz; fi

namain=`basename $genes .gff.gz`

for gfz in $asmset ; do {
  nam=`basename $gfz .gff.gz`
  echo $gfz to $nam.analt.gff

  gzcat $gfz | env mint=$MINTR perl -ne \
  'BEGIN{$MT=$ENV{mint};} if(/\tmRNA/){ unless( ($c,$x)=m/cxlen=(\d+).(\d+)/) { $x=$MT; } $p=($x>=$MT)?1:0; } print if $p;' |\
  $evigene/scripts/overlapfilter.perl -intron2splice=error -pass 'exon,intron' -act markid -midtype scoresum \
    -mark intr -over $introns -in stdin > $nam.analt.gff
} done

#? exit

cat *.analt.gff > allalt.gff

echo $evigene/scripts/altbest.pl $altopts -main $genes -alt allalt.gff -eqtab $namain.equalalt.tab
$evigene/scripts/altbest.pl $altopts -main $genes -alt allalt.gff -eqtab $namain.equalalt.tab


