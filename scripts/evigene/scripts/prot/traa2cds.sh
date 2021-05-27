#!/bin/bash
# traa2cds.sh

if [ "X" = "X$evigene" ]; then
evigene=/bio/bio-grid/mb/evigene/
fi

trset=$*

for trf in $trset; do {
  pt=`echo $trf | sed 's/.gz//; s/\.tr//; s/\.fa.*//;'`

  if [ -f $pt.cds -o -f $pt.cds.gz ]; then 
    echo "# done already $pt.cds"
  else
    aaf=$pt.aa
    if [ -f $pt.aa.gz ]; then aaf=$pt.aa.gz; fi
    if [ -f $aaf ]; then
      echo "# traa2cds $pt.cds"
      $evigene/scripts/prot/traa2cds.pl -cdna $trf -aa $aaf -out $pt.cds -log 
    else
      echo "ERR: missing $pt.aa";
    fi
  fi

} done


