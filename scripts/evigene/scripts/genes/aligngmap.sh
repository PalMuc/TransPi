#!/bin/bash
# aligngmap.sh my.gmap.out  -S format
# rewrite in perl ..

if [ "X" = "X$evigene" ]; then
  if [ -d /bio/bio-grid/mb/evigene/ ]; then
    evigene=/bio/bio-grid/mb/evigene
  else
    echo "ERR: need evigene=/path/to/evigene"; exit -1; 
  fi
fi

#off: noview=1
#off: doequal=0

#201703: keepat="aalen,oid,offs,clen"
#201705upd gmap2evgff: -cututrchim and -cututrx=2 can be useful but also make mistakes ..
eggopt="-nocututrchim -cututrx=0 -keepat 'aalen,oid,offs,clen' -nopath -noerrspan -best=0 -intron=-1"

for gmapz in $*; do {

if [ $gmapz = "-h" ]; then 
  echo "usage: aligngmap.sh mygenes.gmap.out ; gmap -S format"; exit 0;
fi

gn=`basename $gmapz`
pt=`basename $gmapz .gz`
CAT='gunzip -c'; if [ $gn = $pt ]; then CAT=cat; fi
pt=`echo $pt | perl -pe 's/\.gmap.*//; s/\.(mrna_pub|mrna|cdna|tr|fa)-/-/g; '`

##? keep dirpath on pt?
src=$pt
gd=`dirname $gmapz`
pt=$gd/$pt

if [ -f $pt.align.tab ]; then continue ; fi

## FIXME: need som best/nobest; loosing best parts of broken map genes.. eg Funhe5EG022567t1
## off: nobest=1 strand=1; add  best=0 == dups of equal map qual, not lower
## FIXME: skiplow=0 added so no skipped records
## FIXME align.tab or map.attr, trim out C1:xxx,.,0.8 tiny overhangs .. option?
## 2017.03: switch to genes/gmap2evgff.pl from gmap_to_gff.pl
## gmap2evgff: cututr=2 add? 
## skiplow=0 nopath=1 best=0 noerrspan=1 intron=-1 << make more of these defaul
## NOTE cututr is problematic, not all cuts are gene joins (half not?) some are CDS exons that have map indel
#     and/or true long utr.

#above# eggopt="-nocututrchim -cututrx=0 -keepat 'aalen,oid,offs,clen' -nopath -noerrspan -best=0 -intron=-1"
$CAT $gmapz | $evigene/scripts/genes/gmap2evgff.pl $eggopt -src $src  > $pt.gmap.gff

#o $CAT $gmapz | env keepat=$keepat src=$src skiplow=0 nopath=1 best=0 strand=0 noerrspan=1 intron=-1  \
#o $evigene/scripts/gmap_to_gff.pl  > $pt.gmap.gff

cat $pt.gmap.gff | $evigene/scripts/ests/gff2aligntab.pl | sort -k1,1 -k2,2n -k4,4 \
 > $pt.align.tab

## add -f7 = clen
## add this brief align.tab; BUT that gd,ti split only good for evg000t123 ids.. need option
## drop "spx/" part of "$spx/$nx" exon col, too messy here
## add ,anti[sense] flag to npath
#o: if(/AQueryID/){ $spx=$sp; $gd=$td; $ti=0; } else { $spx=int(0.5 + $sp/2); ($gd,$ti)= $td=~/(\w+t)(\d+)$/; } 
#o: cut -f1-4,8,9,10,12,13,18,19 $pt.align.tab | ..

perl -ne 'BEGIN{ @IA=map{$_ -1}(4,7,8,9,10,12,13,17,18,19); }
@v=split; $mloc=join":",@v[0,1,2]; ($td,$clen,$cov,$pi,$np,$nx,$sp,$asens,$oid,$tag)=@v[@IA]; 
unless(/AQueryID/){ ($cov,$pi)=map{ int(0.5+$_) }($cov,$pi); $np.=",anti" if($asens); }
print join("\t",$td,$cov,$pi,$nx,$mloc,$np,$oid,$tag,$clen)."\n";' \
  $pt.align.tab | sort  > $pt.map.attr

if [ $doequal ]; then
 #old# $evigene/scripts/equalgene.pl -this -in $pt.gmap.gff -over $pt.gmap.gff > $pt.ovself.eqgene
 $evigene/scripts/equalgene.pl -selfsame -mrnatype mRNA,ncRNA  -this -oneexonstrandless -in $pt.gmap.gff > $pt.ovself.eqgene

 perl -ne '@v=split; ($td,$ovd)=@v[0,2]; $td=~s/_C\d//;  ($gd=$td)=~s/t\d+//; $gda=$gd; 
 @ov=grep/\w/, map{ ($d,$v)=split"/"; $d=~s/_C\d//;  ($g=$d)=~s/t\d+//; $v=~s/[CI]//; $ok=($d=~/$gda/)?0:1; 
 $gda.="|$g"; ($ok)?$_:""; } split",",$ovd; @ov=("na") unless(@ov); $v[2]=join",",@ov; 
 print join("\t",@v)."\n"; ' \
    $pt.ovself.eqgene > $pt.ovselflocs.eqgene

fi

} done

#.. for gbrowse maps only..
#m if [ "X" = "X$noview" ]; then
#m cat $pt.gmap.gff | perl -ne \
#m 'if(/^\W/){ print unless(/^#i/); next; }  ($r,$src,$tp,$b,$e,$v,$o)=split; 
#m @al= m/;(aalen=)/g; if(@al>1) { s/aalen=\d+;//; } $isp=0; 
#m if(/\t(ID|Parent)=(\w+)_C([12])/){ $at=$1; $d=$2; $isp=$3; s/${d}_C$isp/$d;Split=$isp/; 
#m if(/mRNA/) { $sat=$sat{$d}; if($sat) { s/aalen=\d+;//; s/$/;$sat/; } 
#m else { $sat=join";", m/((?:$aspl)=[^;\s]+)/g; $sat{$d}=$sat if($sat); } } } 
#m ($ap,$d)=/\t(ID|Parent)=([^;\s]+)/;  $mid=($d=~/t1$/)?1:0; $sn=$src; $sn.="pub" unless($sn=~/pub/); 
#m $sn=~s/pub/alt/ unless($mid);  s/\t$src/\t$sn/; s/;($acut)=[^;\n]+//g; print; 
#m BEGIN{ $acut="trg|Target|match|chimera|cdnabest|nover|cdsover|ocds|svec"; $aspl="aalen|offs|oid"; $IDP=""; }' \
#m  > $pt.view.gff
#m fi
#m 
