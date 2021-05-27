#!/bin/bash
## 2014.04: add ncomplete count

if [ "X" = "X$top" ]; then top=1000; fi
if [ "X" = "X$main" ]; then main=0; fi
inaa=$*

for atab in $inaa; do {
  tmpin=0
  if [ $atab = "stdin" ]; then atab="-"; fi
  if [ $atab = "-" ]; then 
    if [ "X" = "X$nam" ]; then nam="stdin"; fi
    tmpfile=/tmp/aastat$$.in
    cat - > $tmpfile
    atab=$tmpfile ; tmpin=1
  else
    nam=`basename $atab .aa.qual | sed 's/\.count//; s/\.aa//;'`
  fi

# naa=`cat $atab | egrep -v '^#|^total' | wc -l`
if [ $main != 0 ]; then

sort -k2,2nr $atab | env top=$top main=$main nam=$nam perl -ne \
'BEGIN{ $top=$ENV{top}||1000; $main=$ENV{main}||"t1"; $GTOP=($main eq "t")?1:0; $nt=0; }
next unless(/^\w/); ($id,$aw,$nn)=split; if($GTOP) { $gd=$id; $gd=~s/utrorf//; $gd=~s/t\d+$//; next if($did{$gd}++); } 
else { next unless($id=~/$main$/); } $nt++; next if($n >= $top);  
$n++; $sw+=$aw; $sn+=$nn; push @aw,$aw; $nc++ if(/complete/);
END{ $aw=int($sw/$n); $an=int(10*$sn/$n)/10; @aw=sort{$b <=> $a}@aw ; ($mx,$md,$mi)=@aw[0,int($n/2),-1];
print "#$ENV{nam}\t  nt=$nt; average=$aw; median=$md; min,max=$mi,$mx; nfull=$nc; sum=$sw; gaps=$sn,$an\n" ; }'

else

naa=`egrep -cv '^#|^total' $atab`
#skip# echo "# aa-quality for $nam : longest $top of $naa"
cat $atab | egrep -v '^#|^total' | sort -k2,2nr | head -$top | env nt=$naa top=$top nam=$nam perl -ne \
'next unless(/^\w/); ($aw,$nn)=(split)[1,2]; $n++; $sw+=$aw; $sn+=$nn; push @aw,$aw; $nc++ if(/complete/);
END{ $aw=int($sw/$n); $an=int(10*$sn/$n)/10; @aw=sort{$b <=> $a}@aw ; ($mx,$md,$mi)=@aw[0,int($n/2),-1];
$top=$ENV{top}; $nt=$ENV{nt}||$n; 
print "#$ENV{nam}\t  nt=$nt; average=$aw; median=$md; min,max=$mi,$mx; nfull=$nc; sum=$sw; gaps=$sn,$an\n"; }'

fi

if [ $tmpin = 1 ]; then rm $tmpfile; fi

} done

