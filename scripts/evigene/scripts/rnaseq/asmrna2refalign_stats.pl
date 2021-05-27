#!/usr/bin/env perl
# asmrna2refalign_stats.pl
# see: evigene_docwork/evigene_bugs10_aaqual.txt
# cat $pt.aligny.tab | sort -k1,1 -k14,14nr -k6,6nr -k4,4r | env bestlist=0 g1=0 nt=30643 pc=80 na=$pt perl asmrna2refalign_stats.pl
# ref table in: "evg3itr.cd90.aa.qual"

# aligny.tab cols
#      1                      4             6          8                         12
# GeneID  gespan  geor  QueryID quspan  match  qlen  cov  pid  path indels  cdsindel
#     13       14       15     16       17
#   join    aalen   aaqual  aalign  aaholn # aalign : use max of aalign ref, alholn from plantho

BEGIN{ $NRT=$ENV{nt}||30643; $PCOV=$ENV{pc}||80; $QTY=$ENV{qt}||0; $GENET1=$ENV{g1}||0; $BESTID=$ENV{bestlist};
open(F,"evg3itr.cd90.aa.qual"); while(<F>){ ($d,$aw,$ag,$aq,$tw)=split; $refw{$d}=$tw; $refaaw{$d}=$aw; $refaaq{$d}=$aq; } 
close(F); }

while(<>) {
next if(/^GeneID|^\W/);
($gdt,$gbe,$go,$td,$qbe,$ma,$ql,$cov,$pid,$np,$ind,$cind,$join,$aaw,$aaq,$aalign,$aaho)=split;
$refw=$refw{$gdt}; $refaaw=$refaaw{$gdt}; next unless($refw>0 and $refaaw>0); next unless($aalign>0);
$aalign=$aaho if($aaho>$aalign); next unless($aalign>0); $aalign=$aaw if($aalign>$aaw); 
$gd=$gdt; $gd=~s/t\d+$// if($GENET1); next if($did{$gd});
if($QTY) { $qty=substr($td,0,10); $qty=~s/cacao\d+//; $qty="esta" if($qty=~/^(B|L|P1|P2)_/); $qty=~s/\d.*$//;  $qty{$qty}++; }
($qb,$qe)=split/[-]/,$qbe;  $qw=1+$qe-$qb;
$cov=100*$qw/$refw; $scov+= $cov; $nc80++ if($cov>=$PCOV);
$sp+=$pid; $snp+=$np; $sw+=$qw; $sm+=$ma; $strw+= $ql;
$saalen+=$aaw; $saalign+=$aalign; $naac++ if($aaq =~ /complete/);
$cind= abs($cind); $scind+=$cind; $ncind++ unless($cind % 3 == 0); $njoin++ if($join =~ /join/);
$aasizep=100*$aaw/$refaaw; $aacov=100*$aalign/$refaaw; $sacov+= $aacov; $saasizep+= $aasizep; $ac80++ if($aacov>=$PCOV);
$n++; $did{$gd}++; 
$bestid{"$gdt\t$td"}= join("\t",int(0.5+$cov),int(0.5+$aacov));
}

END {
$pcov= int(10*$scov/$n)/10; $pCov= int(1000*$nc80/$NRT)/10;
$pacov= int(10*$sacov/$n)/10; $paCov= int(1000*$ac80/$NRT)/10; 
$pid= int(10*$sp/$n)/10; $aln=int($sw/$n); $trlen=int($strw/$n);  
$pac= int(1000*$sm/$sw)/10; $pasize= int(10*$saasizep/$n)/10;
$aalen=int($saalen/$n); $pCDS=int(3000*$saalen/$strw)/10; $paac=int(1000*$naac/$n)/10;
$cind=int(10*$scind/$n)/10; $pindel=int(1000*$ncind/$n)/10; $pjn=int(1000*$njoin/$n)/10;
$qty=""; if($QTY) { $qty= " qtyp=".join",", map{ "$_:$qty{$_}" } (sort keys %qty); }

$na=substr($ENV{na},0,11); print "# " if($BESTID);
print "$na : naln=$n; Tcov$PCOV=$pCov%; pTCov=$pcov%; pIden=$pid%; ";
print " ACov$PCOV=$paCov%; pACov=$pacov%; pCDS=$pCDS%; pAFull=$paac%; pJoin=$pjn%; pIndel=$pindel%,$cind; ";
print " aln=$aln; aalen=$aalen; trlen=$trlen; $qty\n";
if($BESTID) { print "#BESTID_list\n#GeneID\tQueryID\tpTRcov\tpAAcov\n"; 
foreach $rt (sort keys %bestid) { print join("\t",$rt,$bestid{$rt}),"\n"; } }
} 

__END__

## out tables:
## single-method, sort +trlen +tralign -trid;  effects: -5 ACov/pAcov; -pCDS; +1 pJoin; +pIndel; -10..20 aalen; 
# cacao3cuf13 : naln=15695; Tcov80=20.8%; pTCov=65.9%; pIden=99.3%;  ACov80=34.9%; pACov=81.8%; pCDS=50.9%; pAFull=84.1%; pJoin=21.4%; pIndel=16.1%,10.4;  aln=1327; aalen=390; trlen=2299; 
## single out: sort +aalen +tralign -trid; added pAFull, pIndel here is CDS-frameshift count
# cacao3vel12 : naln=20988; Tcov80=24.9%; pTCov=62.2%; pIden=99%;    ACov80=49.5%; pACov=83.5%; pCDS=68.6%; pAFull=74.4%; pJoin=12.5%; pIndel=30%,9.4;     aln=1243; aalen=400; trlen=1751;  
## bestof
# best6aan    : naln=25984; Tcov80=42.3%; pTCov=68.8%; pIden=98.6%;  ACov80=67.9%; pACov=88%;   pCDS=59.5%; pAFull=82.6%; pJoin=20.6%; pIndel=28.7%,12.8;  aln=1331; aalen=419; trlen=2117;  qtyp=cuf:5820,esta:1452,sopc:2204,tri:5861,velv:3765,velw:6867
#
# common opts above: bestlist=0 nt=30643 pc=80
pt=best6laa # select by longest aa ( .. and trmatch sort ?)
cat *.aligny.tab | sort -k1,1 -k14,14nr -k6,6nr -k4,4r | env qt=1 bestlist=0 nt=30643 pc=80 na=$pt perl asmrna2refalign_stats.pl

pt=best6aan # best aa align, for any orthology selection
cat *.aligny.tab | sort -k1,1 -k16,16nr -k14,14nr -k6,6nr -k4,4r | env qt=1 bestlist=0 nt=30643 pc=80 na=$pt perl asmrna2refalign_stats.pl

pt=best6qtr  # select by longest tr, -k7, and match sort
cat *.aligny.tab | sort -k1,1 -k7,7nr -k6,6nr -k4,4r | env qt=1 bestlist=0 nt=30643 pc=80 na=$pt perl asmrna2refalign_stats.pl

for atab in *.aligny.tab; do { 
  pt=`basename $atab .aligny.tab`; 
  cat $pt.aligny.tab | sort -k1,1 -k7,7nr -k6,6nr -k4,4r | \
  env bestlist=0 nt=30643 pc=80 na=$pt perl asmrna2refalign_stats.pl
} done
