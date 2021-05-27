#!/bin/bash
# env spp=cacao3sop maketrsoap.sh  *scafSeq.gz
# sotri3sub15k45/sotri3sub15.k45.scafSeq.gz
# >scaffold1 Locus_0_0 48.4 FORK 
# >C21245 29.0


MINTR=160
if [ "X$spp" == "X" ]; then spp="species1sop"; fi
#spp=cacao3sop

# trs=`ls $trset*scafSeq.gz`
trs=$*
tr1=$1

rnam=`basename $tr1 .scafSeq.gz | sed "s/\.//; s/sotri3/${spp}/;"`

mv $rnam.tr $rnam.tr.old
touch $rnam.tr

for tr in $trs ; do {
  nam=`basename $tr .scafSeq.gz | sed -e "s/\.//;  s/sotri3/${spp}/; " `
  echo $tr TO $nam

  gzcat $tr | env MINTR=$MINTR rs="$nam" perl -ne 'BEGIN{ $MIN=$ENV{MINTR}||0; $rs=$ENV{rs}; $maxl=0; } 
if(/^>/) { putr() if($hd); $hd=$sq="";
if(m,>scaffold(\d+)\s+Locus_(\d+)_(\d+)\s+(\S+)\s*(\S*),) { ($i,$l,$t,$val,$fl)=($1,$2,$3,$4,$5); 
$maxl=$l if($l>$maxl); ++$t; $fl="; flag=$fl" if($fl); } elsif(m/^>C(\d+)\s*(\S*)/) { ($i,$val)=($1,$2);
$l=++$maxl; $t=1; $fl=""; } else { warn "#Fixme: $_"; }
s/>.*$/>${rs}loc${l}t${t} i=$i; val=$val$fl;/; $hd=$_; }
else { chomp; $sq .= $_; } END{ putr() if($hd); }
sub putr{ my $sl= length($sq); if($sl>=$MIN) { $hd=~s/$/ len=$sl;/; $sq=~s/(.{60})/$1\n/g; 
print "$hd$sq\n"; } $hd=$sq=""; } ' \
  >> $rnam.tr
  
}
done

