#!/bin/bash

noview=1
splign=1
STAG=gspl15

## FIX This, option..
## if [ "X" = "X$noview" ]; then noview=0; fi
## if [ "X" = "X$splign" ]; then splign=0; fi
## fix for kfish2n_splign15.gff ..

evigene=/bio/bio-grid/mb/evigene
keepat="aalen,oid,offs"
inmap=$*

for gmapz in $inmap; do {

## may be gff not gz
pt=`basename $gmapz .gz | perl -pe 's/\.gff.*//; s/\.gmap.*//; s/\.(mrna_pub|mrna|tr|fa)-/-/; s/\-/pub-/ unless(/pub/);'`

if [ -f $pt.align.tab ]; then continue ; fi

## FIXME: need som best/nobest; loosing best parts of broken map genes.. eg Funhe5EG022567t1
## off: nobest=1 strand=1; add  best=0 == dups of equal map qual, not lower
## FIXME: skiplow=0 added so no skipped records

if  [ "X" = "X$splign" ]; then
  ptgff=$pt.gmap.gff

  gunzip -c $gmapz | env keepat=$keepat src=$pt skiplow=0 nopath=1 best=0 strand=0 noerrspan=1 intron=-1  \
   $evigene/scripts/gmap_to_gff.pl  > $ptgff

else
  ptgff=$gmapz
fi

cat $ptgff | $evigene/scripts/ests/gff2aligntab.pl | sort -k1,1 -k2,2n -k4,4 \
 > $pt.align.tab

## fixme add tag: perl -pi -e 's/$/\tgspl15/;' $oname.map.attr
## add this brief align.tab; BUT that gd,ti split only good for evg000t123 ids.. need option
cut -f1-4,8,9,10,12,13 $pt.align.tab | env src=$STAG perl -ne \
'BEGIN{$STAG=$ENV{src}||"gspl0";} 
@v=split; $m=join":",@v[0,1,2]; ($td,$cov,$pi,$np,$nx,$sp)=@v[3..8]; if(/AQueryID/) { $spx=$sp; $gd=$td; $ti=0; } 
else { $cov=~s/,.*//; $spx=int(0.5 + $sp/2); ($gd,$ti)= $td=~/(\w+t)(\d+)$/; } 
print join("\t",$td,$cov,$pi,"$spx/$nx",$m,$np,$STAG)."\n";' | sort \
 > $pt.map.attr

# print join("\t",$gd,$ti,$cov,$pi,"$spx/$nx",$m,$np)."\n";' | \
# sort -k1,1 -k2,2n | sed 's/t	/t/;' > $pt.map.attr

if [ "X" = "X$noview" ]; then
cat $ptgff | perl -ne \
'if(/^\W/){ print unless(/^#i/); next; }  ($r,$src,$tp,$b,$e,$v,$o)=split; 
@al= m/;(aalen=)/g; if(@al>1) { s/aalen=\d+;//; } $isp=0; 
if(/\t(ID|Parent)=(\w+)_C([12])/){ $at=$1; $d=$2; $isp=$3; s/${d}_C$isp/$d;Split=$isp/; 
if(/mRNA/) { $sat=$sat{$d}; if($sat) { s/aalen=\d+;//; s/$/;$sat/; } 
else { $sat=join";", m/((?:$aspl)=[^;\s]+)/g; $sat{$d}=$sat if($sat); } } } 
($ap,$d)=/\t(ID|Parent)=([^;\s]+)/;  $mid=($d=~/t1$/)?1:0; $sn=$src; $sn.="pub" unless($sn=~/pub/); 
$sn=~s/pub/alt/ unless($mid);  s/\t$src/\t$sn/; s/;($acut)=[^;\n]+//g; print; 
BEGIN{ $acut="trg|Target|match|chimera|cdnabest|nover|cdsover|ocds|svec"; $aspl="aalen|offs|oid"; $IDP=""; }' \
 > $pt.view.gff
fi

} done

