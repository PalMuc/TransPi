#!/bin/bash
# aacount.sh any*.aa.gz : count protein sizes with faCount (dna), with tr for XXX gap count
# output 3-col table: ID, aatotal, XXgaps

inaa=$*
for az in $inaa; do { 
  nam=`basename $az .gz`; 
  # if .gz then gunzip .. else cat ..
  if [ -f $nam.count ]; then continue; fi
  echo aaCount $az
  gunzip -c $az | perl -pe 'unless(/^>/){ s/X/n/g; s/[A-Z]/a/g;}' |\
  faCount stdin | cut -f1,2,7 > $nam.count; 
} done

## added calcs: cd-hit -c 0.90 .. aacount .. stats of top 1000 prots:
# cd-hit -c 0.9 -d 0 -i $pt.aa -o ${pt}_cd90.aa > & log.cd9$pt 
# $evigene/scripts/prot/aacount.sh $pt.aa.gz
# cat $pt.aa.count | egrep -v '^#|^total' | sort -k2,2nr | head -1000 | env nam=$pt perl -ne \
#   '($aw,$nn)=(split)[1,2]; $aw= $aw-$nn;  $n++; $sw+=$aw; $sn+=$nn; push @aw,$aw; END{ $aw=int($sw/$n); $an=int(10*$sn/$n)/10; @aw=sort{$b <=> $a}@aw ; $md=$aw[int($n/2)]; print "$ENV{nam}\t  n=$n; aw=$aw; med=$md; sw=$sw; sn=$sn,$an\n"; }'
#
