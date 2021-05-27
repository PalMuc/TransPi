#!/bin/bash

# trs=`ls trout*/Trinity.fasta`
trs=$*

if [ "X$prefix" = "X" ]; then
prefix=trin
fi

for tr in $trs ; do {
  trnogz=`echo $tr | sed 's/\.gz//'`; 
  TCAT=cat; if [ $tr != $trnogz ]; then TCAT="gunzip -c"; fi
  dnam=`echo $trnogz | sed -e"s,/Trinity.fasta,,; s/\.fa//; s/\./_/g; s,trout,,; s,/,,; s,/,-,g; " `
  onam=$prefix$dnam
  #or this# onam=$dnam$prefix
  echo $tr TO $onam

  $TCAT $tr | env s=$onam  perl -pe 'BEGIN{ $rs=$ENV{s}||"s"; } 
  if(m/^>/){ s/[\[\]]//g; s/(path=\S+) /$1,/g; s/comp/loc/; s/_c/c/; s/_seq/t/; s/>/>$rs/; }' \
  > $onam.tr
  
}
done

