#!/bin/bash 
# NOT gg_job script; run as per est/gmap.sh

cd $HOME/scratch/chrs/cacao/rnas
set dgenome=cacao1asm
# this is PE data .fa2
#paired-end  --pairlength=200 default; --pairmax=1000 default << increase to 5000 or 10000 ?

foreach rseq (M*.fa2.gz)
  set nam=`basename $rseq .fa2.gz`
  gzcat $rseq |\
  gsnap -N 1 -k 14 --local-splice-penalty=1 --indel-penalty=1 --pairmax=5000 \
    -A "sam" -D ../genome/gmap -d $dgenome > $nam.gsnap.samu

  cat $nam.gsnap.samu | grep -v '^@SQ' | sort -k 3,3 -k 4,4n > $nam.gsnap.sam
  echo /bin/rm $nam.gsnap.samu
end
