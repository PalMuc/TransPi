#!/bin/bash
### bestaa3.sh : find best proteins (ORFs) in transcript set, and cluster to non-redundant longest aa set
### revise to subset ok,poor.aa by aaqual:%CDS per rnaseq/maketraa3.sh

## use aa flag utrbad/utrpoor here, 
#  1. segregate: ok.aa, poor.aa , 
#  2. cd-hit -i ok.aa -o ok_cd.aa
#  3. cd-hit-2d -i ok_cd.aa -i2 poor.aa -o poor_cd.aa
# aminbad = 200 may be too high for some real genes, protein-coding + long ncrna-utr; check this.

function usage() {
  echo bestaa.sh : find best proteins, ORFs, in transcript set, and cluster to non-redundant longest aa set
  echo usage: env name=bestof PI=90  bestaa.sh  proteins1.fasta.gz proteins2.fasta
  echo part of http://eugenes.org/EvidentialGene/
  exit -1
}

# AAMIN=30; AMINBAD=200; AMINPOO=100
if [ "X$AAMIN" = "X" ]; then AAMIN=30; fi
if [ "X$AMINBAD" = "X" ]; then AMINBAD=200; fi
if [ "X$AMINPOO" = "X" ]; then AMINPOO=100; fi

# cd-hit cluster identity cutoff: 90, 80, 60 good, but may need adjust other param for <75
# * evigene/rnaseq/maketraa.sh now uses PI=95 for 1st cluster; need that to compare single methds.
if [ "X$PI" = "X" ]; then PI=90; fi
if [ "X$PI2" = "X" ]; then PI2=$PI; fi

idir=1
if [ "X$name" = "X" ]; then
name="bestof";
idir=0
fi

## pbs qsub: use env aaset=xxx instead of 
if [ "X$aaset" = "X" ]; then aaset=$*; fi

if [ "X$evigene" = "X" ]; then evigene=/bio/bio-grid/mb/evigene; fi

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
  touch $name.ok.aa $name.poor.aa

  ## FIXME: allow for aa w/o evigene header info, if no aalen, count them, ok if no aaqual
  TCAT=cat;
  fan=`echo $faa | sed 's/\.gz//'`; if [ $faa != $fan ]; then TCAT="gunzip -c"; fi
  $TCAT $faa | env aamin=$AAMIN mpoo=$AMINPOO mbad=$AMINBAD nam=$name perl -ne\
'BEGIN { $aamin=$ENV{aamin}||20; $mpoo=$ENV{mpoo}; $mbad=$ENV{mbad}; $nam=$ENV{nam}; 
open(OK,">>$nam.ok.aa"); open(BAD,">>$nam.poor.aa"); }
if(/^>/) { $bad=(/utr(poor|bad)/)?1:0; ($al)=m/aalen=(\d+)/; $skip=0; 
if($al>0) { ($ap)= m/aalen=$al,(\d+)/; $skip=1 if($al < $aamin);
if($bad) { $skip=1 if((/utrbad/ and $al<$mbad) or (/utrpoor/ and $al<$mpoo)); }
elsif($ap) { $bad=1 if($ap<60); $skip=1 if(($ap<=33 and $al<$mbad) or ($ap<60 and $al<$mpoo)); } } }
if($skip){ } elsif($bad) { print BAD $_; } else { print OK $_; }'

} done

#old# $cdhitapp -c 0.$PI -d 0 -i $name.aa -o $onam.aa >& $onam.log
$cdhitapp -c 0.$PI -d 0 -i $name.ok.aa -o $onam.ok_cd.aa >& $onam.cd1.log
$cdhitapp -c 0.$PI2 -d 0 -i $name.poor.aa -o $onam.poorcd1.aa  >& $onam.cd2a.log 
$cdhit2dapp -d 0 -c 0.$PI2 -i $name.ok.aa -i2 $onam.poorcd1.aa -o $onam.poor_cd.aa >& $onam.cd2.log

/bin/rm $onam*.bak.clstr; /bin/rm $onam.poorcd1.aa*

cat $onam.{ok,poor}_cd.aa > ${onam}_allcd.aa
  
grep '^>' ${onam}_allcd.aa | sed 's/>//; s/ .*//;' > ${onam}_allcd.ids
env stat=1 $evigene/scripts/prot/aaqual.sh   ${onam}_allcd.aa

mkdir cdhits
mv $onam.cd*.log cdhits/
mv $onam.{ok,poor}* cdhits/
mv $name.{ok,poor}.aa cdhits/
gzip --fast ${onam}_allcd.aa cdhits/*.aa

