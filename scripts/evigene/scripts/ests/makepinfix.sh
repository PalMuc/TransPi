#!/bin/bash

evigene=/bio/bio-grid/mb/evigene
caca=/bio/bio-grid/cacao3
dna=$caca/genome/cacao11allasm.fa
intron=$caca/intron/intron_good.gff.gz

for gff in cgbAssembly*.gmap.gff.gz; do 
{
 gset=`echo $gff | sed 's/.gmap.gff.gz//'`
 echo $evigene/scripts/genefindcds.pl -dna $dna -intron $intron -genes $gff TO $gset.pinfix.gff
 $evigene/scripts/genefindcds.pl -dna $dna -intron $intron -genes $gff > $gset.pinfix.gff

} done

