#!/bin/bash
# gmapgfftab.sh / aligngmap.sh : convert GMAP --summary to GFF3 + table per gene of alignment

if [ "X" = "X$evigene" ]; echo "ERR evigene=where?"; exit -1; fi

ingmap=$*
keepat="aalen,oid,offs"

for gmapz in $ingmap; do {
  TCAT=cat;
  nogz=`echo $gmapz | sed 's/.gz//;'`; if [ $gmapz != $nogz ]; then TCAT="gunzip -c"; fi
  pt=`basename $gmapz .gz | sed 's/\.gmap.*//; s/\.mrna_pub-/-/; s/\.mrna-/-/; s/\.fa-/-/; s/\.tr-/-/;'`
	if [ -f $pt.align.tab ] ; then continue; fi

	if [ ! -f $pt.gmap.gff] ; then 
	 $TCAT $gmapz | env keepat=$keepat src=$pt skiplow=0 nopath=1 best=0 strand=0 noerrspan=1 intron=-1  \
    $evigene/scripts/gmap_to_gff.pl  > $pt.gmap.gff
	fi
	
	grep mRNA $pt.gmap.gff | perl -ne \
'BEGIN{ @k=qw(match qlen cov pid path indels cdsindel);
print join("\t","GeneID","gespan","geor","QueryID","quspan",@k)."\n"; }
@v=split"\t"; ($gid,$gb,$ge,$go,$at)=@v[0,3,4,6,-1];
 @at=map{ $v=($at=~m/\b$_=([^;\s]+)/)?$1:0; $v; } ("ID",@k);
($tid,$tb,$te)=m/Target=(\S+) (\d+) (\d+)/; $id=shift @at;
print join("\t",$gid,"$gb-$ge",$go,$tid,"$tb-$te",@at)."\n";' \
 | sort -k1,1 -k2,2n -k4,4 > $pt.align.tab

} done


