#!/bin/bash
# env refname=../aaset/uniref50arptop16.aa.names ./makename.sh outz/sd-uniref50arptop16*deblastp.gz

if [ "X" = "X$evigene" ]; then evigene=/bio/bio-grid/mb/evigene; fi
if [ ! -d $evigene ]; then "ERR: env evigene=path/to/evigene"; exit -1; fi

# refname=../aaset/uniref13arpod.aa.names
if [ "X" = "X$refname" ]; then echo "ERR: env refname=../aaset/uniref50.aa.names"; exit -1; fi

# cddname=../aaset/info.cdd.txt
if [ "X" = "X$cddname" ]; then cddname=`dirname $refname | sed 's,$,info.cdd.txt,;'`; fi
if [ -f $cddname ]; then 
 cddopt="-cddname $cddname"
else
 cddopt=""; echo "NO env cddname=..., skipping CDD naming"; 
fi

if [ "X" = "X$suf" ]; then suf="names.tab"; fi

blist=$*

for bz in $blist; do {
  pt=`basename $bz .deblastp.gz | sed 's/sd-//; s/\.aa.*//; s/top16/top/; s/\.okay/okay/; s/\.okboth/okboth/;'`
  if [ -f $pt.$suf ]; then echo "DONE already $pt.$suf"; continue; fi

  echo $evigene/scripts/prot/namegenes.pl -form consen $cddopt -refname $refname -blast $bz -out $pt.$suf
  $evigene/scripts/prot/namegenes.pl -form consen $cddopt -refname $refname -blast $bz -out $pt.$suf

} done

