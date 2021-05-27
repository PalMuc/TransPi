#!/bin/bash
# ./veltrmake.sh  veltrdir

trdir=$1
spt=$2
if [ ! -d $trdir ];   then echo "usage: $0  trdir  speciestag"; exit -1; fi
if [ "X" = "X$spt" ]; then echo "usage: $0  trdir  speciestag"; exit -1; fi

## rename here or before?
trdd=${trdir}trs
mkdir $trdd
trdd=`basename $trdd`
cd $trdir
for pt in vel*_??; do { mv $pt/transcripts.fa ../$trdd/$pt-transcripts.fa; } done; 
cd ../
gzip --fast $trdd/*transcripts.fa

scs=""

trs=`ls $trdd/vel*_??-transcripts.fa.gz`
if [ "X$trs" = "X" ]; then exit -1; fi

trset=($trs)[1]
if [ ! -f $trset ]; then exit -1; fi
rnam=`basename $trset -transcripts.fa.gz | sed "s/_..//; s/^/$spt/;"`
echo "$trset to $rnam.tr"

if test -f $rnam.tr; then { mv $rnam.tr $rnam.tr.old; } fi
touch $rnam.tr

for tr in $trs ; do {
  nam=`basename $tr -transcripts.fa.gz | sed -e "s/_/k/; s/^/$spt/;" `
  echo $tr TO $nam

  gzcat $tr | env rs="$nam" perl -pe 'BEGIN{ $rs=$ENV{rs};} if(/^>/) {
s,/(\d+)_Confidence_, nt=$1; cf=,; s/_Length_/; len=/; s/Locus_/Loc/; s/_Transcript_/t/; s/>/>$rs/;}' \
  >> $rnam.tr
  
} done
