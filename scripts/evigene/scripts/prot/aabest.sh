#!/bin/bash
### bestaa.sh : find best proteins (ORFs) in transcript set, and cluster to non-redundant longest aa set
### revise to subset ok,poor.aa by aaqual:%CDS per rnaseq/maketraa3.sh

# cd-hit cluster identity cutoff: 90, 80, 60 good, but may need adjust other param for <75
# * evigene/rnaseq/maketraa.sh now uses PI=95 for 1st cluster; need that to compare single methds.
if [ "X$PI" = "X" ]; then
PI=95
fi

idir=1
if [ "X$name" = "X" ]; then
name="bestof";
idir=0
fi

function usage() {
  echo bestaa.sh : find best proteins, ORFs, in transcript set, and cluster to non-redundant longest aa set
  echo usage: env name=bestof PI=95 bestaa.sh  proteins1.fasta.gz proteins2.fasta
  echo part of http://eugenes.org/EvidentialGene/
  exit -1
}

## pbs qsub: use env aaset=xxx instead of 
if [ "X$aaset" = "X" ]; then
aaset=$*
fi

if [ "X$evigene" = "X" ]; then
evigene=/bio/bio-grid/mb/evigene
fi

cdhitapp=`which cd-hit`
cdhit2dapp=`which cd-hit-2d`
if [ "X$cdhitapp" = "X" ]; then
cdhitapp=/bio/bio-grid/mb/bin/cd-hit
cdhit2dapp=/bio/bio-grid/mb/bin/cd-hit-2d
fi

if [ ! -d $evigene/scripts/ -o ! -x $cdhitapp ]; then
  echo "ERR: need path to evigene/scripts/ or cd-hit "; usage;
fi

suf="_cd$PI"
onam=${name}$suf

# idir=0
for faa in $aaset; do { 
  if [ ! -f $faa ]; then
    echo "ERR: missing input proteins: $faa"; usage;
  fi

  idir=$(( $idir + 1 ))
  if [ $idir == 1 ]; then
    tdir=`dirname $faa`
    cd $tdir
  fi
  touch $name.aa

  TCAT=cat;
  fan=`echo $faa | sed 's/\.gz//'`; if [ $faa != $fan ]; then TCAT="gunzip -c"; fi
  $TCAT $faa  >> $name.aa

} done

$cdhitapp -c 0.$PI -d 0 -i $name.aa -o $onam.aa >& $onam.log
/bin/rm $onam*.bak.clstr; 
  
grep '^>'  $onam.aa | sed 's/>//; s/ .*//;' > $onam.ids
env stat=1 $evigene/scripts/prot/aaqual.sh  $onam.aa

