#!/bin/bash
# aaqual.sh any*.aa.gz : count protein sizes with faCount (dna), with tr for XXX gap count
# output 3-col table: ID, aatotal, XXgaps
# add on >aalen qual columns; ** Only good for evigene cdna_bestorf.aa files
# >whitefly:vel2k25Loc100074t1 aalen=41,63%,complete; clen=199; strand=-; offs=144-19; 
# ** FIXME2: w/ new aaqual, set col2 == size - gap, not sizetotal, and let user add gaps in
## FIX for ensembl pep.fa with gene: instead of oid= gene:ENSDARG00000030494
## FIX3: add strand to span/offs !!
## FIX4: Selcstop=1 transfer to aaqual string
## FIX4: cdsoff= alias offs=; cxlen=cdslen/trlen alias to clen=; dang also cdsoffs=nnn

## option: offs= column, strand= ??
dostat=0; if [ "X$stat" != "X" ]; then dostat=$stat; fi
outd=0; if [ "X$outdir" != "X" ]; then outd=$outdir; fi
#x if [ "X$off" != "X" ]; then off=$off; elif ..
if [ "X$span" != "X" ]; then export off=$span; fi; 
#dont need# if [ "X$oid" = "X" ]; then export oid=0; fi
#dont need# if [ "X$ismrna" = "X" ]; then export ismrna=0; fi

inaa=$*
for az in $inaa; do {
  TCAT=cat; 
  nogz=`echo $az | sed 's/.gz//;'`; if [ $az != $nogz ]; then TCAT="gunzip -c"; fi
  nam=`echo $az | sed 's/.gz//; s/\.qual//;'`; 
  if [ $outd != 0 ]; then nam=`basename $nam | sed "s,^,$outd/,;"`; fi

  if [ ! -f $nam.qual ]; then 
  $TCAT $az | perl -ne \
'if(/^>(\S+)/) { puta() if($d); $d=$1; ($al)=m/aalen=([^;\s]+)/; $al||="na"; 
if(/Selcstop=\w/) { $al.=",selcstop" unless($al=~/selc/); }
($cl)=m/clen=(\d+)/; unless($cl) { ($cdl,$cl)=m/cxlen=(\d+).(\d+)/; } $cl||=0; 
$aas=$aam=$aat=$aag=0; unless(($oid)=m/oid=([^;\s]+)/) { ($oid)=m/gene[=:]([^;\s]+)/; } $oid||="noid";
if($doff){($or)=m/strand=(.)/; $or||="."; ($ofs)=m/\b(?:offs|cdsoff\w*)=([\d-]+)/; $cl.= ($ofs)?"\t$ofs:$or":"\t0";}}
else { $aas=(s/\*$//)?1:0; if($ismrna) { $aat += tr/ACGTacgt/ACGTacgt/; $aag+= tr/Nn/Nn/; } 
else { $aat += tr/A-WYZa-wyz/A-WYZa-wyz/; $aag += tr/Xx\*/Xx\*/; if($aam==0){ $aam=(/^M/)?1:-1; } } }
END{ puta(); } BEGIN{ $doff=$ENV{off}||0; $doid=$ENV{oid}||0; $ismrna=$ENV{ismrna}||$ENV{mrna}||0; }
sub puta { if($al eq "na"){ $part=($aas==1 && $aam==1)?"complete":($aam==1)?"partial3":"partial"; 
 $al=($aat+$aag).",na,$part"; } 
$cl.="\t$oid" if($doid); print join("\t",$d,$aat,$aag,$al,$cl)."\n"; }' \
  > $nam.qual

  fi

  # add stat here if desired:
  if [ $dostat != 0  ]; then
  pnam=`basename $nam | sed 's/\.aa//; s/.allcd//;'`
echo "# aa-quality for $pnam : top 1k summary"
cat $nam.qual | egrep -v '^#|^total' | sort -k2,2nr | head -1000 | env nam=$pnam perl -ne \
'($aw,$nn)=(split)[1,2]; $n++; $sw+=$aw; $sn+=$nn; push @aw,$aw;
END{ $aw=int($sw/$n); $an=int(10*$sn/$n)/10; @aw=sort{$b <=> $a}@aw ; ($mx,$md,$mi)=@aw[0,int($n/2),-1];
print "#$ENV{nam}\t n=$n; average=$aw; median=$md; min,max=$mi,$mx; sum=$sw; gaps=$sn,$an\n"; }'
  fi

} done

