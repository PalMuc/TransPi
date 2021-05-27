#!/bin/bash

## change to this:
#  pt=kfish2rae5h.locov90t1.mblastn.split.15; 
# $evigene/scripts/rnaseq/splign2gff2.pl -in $pt.splign -log $pt.splog -trinfo kfish2rae5h.locov90t1.mrna -debug -nosorted -out $pt.splgff2b

upd=15h
oname=kfish2nsplign$upd
tsrc=splkf3n$upd
sgopt="-MINIDENT=0.80  -src $tsrc"

evigene=/bio/bio-grid/mb/evigene
trhead=kfish2rae5h.locov90t1.mrna

for spt in spldt1p*/kf*.splign ; do {

 spl=`echo $spt | sed 's/.splign/.splog/;'`
 $evigene/scripts/rnaseq/splign2gff2.pl $sgopt -in $spt -log $spl -trinfo $trhead -out

} done

cat spldt1p*/kf*.gff | perl -pe 's/Target=/trg=/;' > $oname.gff

# ./makesplgff.sh
# ./aligngmap.sh $oname.gff

#fixed# perl -pi -e 's/$/\tgspl15/;' $oname.map.attr

## rename gff gff.old spldt1p*/kf*.gff 
# perl -pi -e's/\tsplign/\tsplkf2p67vs/; s/Target=/trg=/;' $oname.gff

#.................
# info for kfish2p67vs.mrna gene set
# gmap nopath:14950  splign nomap:30106 of ntr=137660 ids
# .. check 'good' gmap chimera on 2 scaffolds versus splign

