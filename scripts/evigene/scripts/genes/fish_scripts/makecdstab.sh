#!/bin/bash

fcvers="_fcds7"

kfish2=/bio/bio-grid/kfish2
gasm2a=$kfish2/genome/killifish20130322asm.fa
gasm3n=$kfish2/submitf/pubgenome/ncbifunhe302scaf.fa

pubaaqual=kfish2rae5h.*.pub.aa.qual
# pubaaqual=$kfish2/submitf/pubgenes/kfish2rae5h.*.pub.aa.qual
xidtab=../kfish2x11tsa.pubidtab

if [ "X" = "X$evigene" ]; then evigene=/bio/bio-grid/mb/evigene; fi

if [ "X" = "X$suf" ]; then suf=cds2aatab7; fi
if [ "X" = "X$findcds" ]; then findcds=0; fi
if [ "X" = "X$maketab" ]; then maketab=0; fi
if [ "X" = "X$makesum" ]; then makesum=0; fi
if [ "X" = "X$makeendgaps" ]; then makeendgaps=0; fi
if [ "X" = "X$doalt" ]; then doalt=0; fi

## add source tag col: gmap2pub5h ..
#i ==> kfish2nsplign15h_fcds4.cds2aatab <== #i Funhe2EKm000014t1       gspl3n5h
#i ==> kfish2nsplign15n_fcds4.cds2aatab <== #i Funhe2EKm000004t1       gspl3n5n
#i ==> kfish2rae5h_asm3n_fcds4.cds2aatab <== #i Funhe2EKm000003t1       gmap3npub5h
#i ==> kfish2rae5h_kfish2a_t1fcds4.cds2aatab <== #i Funhe2EKm027243l2t1     gmap2pub5h
#i ==> kfish2x11tsa_pubid_asm3n_fcds4.cds2aatab <== #i Funhe2EKm003823t1       gmap3nx11tsa
## add cov,pid map stats cols
## handle errs: ocds=cancel:NEn,aatrans-failed, ie no cds2aa; 
## now get cds2alen = 0, but only 1 has protein= byhand, so most are missing from output
#.. aatrans-failed:
# kfish2rae5h_asm3n_fcds4.gff.gz:9018
# kfish2rae5h_kfish2a_t1fcds4.gff.gz:2268
# kfish2rae5h_kfish2a_tafcds4.gff.gz:6712
# kfish2x11tsa_asm3n_fcds4.gff.gz:19725
##
## drop field= tags for header col?
## Funhe2EKm008403l2t1	0.da	aalen=78,10%,complete	aaold=78,10%,complete	cdsoff=546-782	offp=0	clen=0	cxlen=234/2332	0	gmap2pub5h
## AQueryID ADiff cds2aalen  pubaalen  cdsoff  puboff  pubclen  cxlen cov pid path maptag

if [ $maketab = 1 ]; then

# for gfz in *$fcvers.gff.gz; do
for gfz in $*; do 
{
  pt=`basename $gfz .gff.gz`
  echo "$gfz TO $pt.$suf"
  gunzip -c $gfz | grep mRNA | grep protein= | \
  cat $pubaaqual - | env pt=$pt xidtab=$xidtab perl -ne \
'BEGIN{ $pt=$ENV{pt}; $STAG="none"; 
if($pt=~/splign15h/) { $STAG="gspl3n5h"; } 
elsif($pt=~/splign15n/) { $STAG="gspl3n5n"; } 
elsif($pt=~/5h_asm3n/) { $STAG="gmap3npub5h"; } 
elsif($pt=~/5h_kfish2a/) { $STAG="gmap2pub5h"; }  # fcds5 update to kfish2rae5ht1asm2_fcds5.*
elsif($pt=~/5ht.asm2/) { $STAG="gmap2pub5u"; }  # fcds5 update to kfish2rae5ht1asm2_fcds5.*
elsif($pt=~/x11tsa/) { $STAG="gmap3nx11tsa"; 
 if( open(F,$ENV{xidtab}) ){ while(<F>){ ($xd,$pd)=split; $xpd{$xd}=$pd; $nxpd++; } close(F);}
 } 
} 
if(/\tmRNA/) { ($id)=m/ID=([^;\s]+)/; ($ids=$id)=~s/_C\d+$//; 
if($nxpd) { $pd=$xpd{$ids} or next; $id=~s/$ids/$pd/; $ids=$pd; }
@at=(); for $k (qw(aalen aaold cdsoff offs clen cxlen cov pid)) { 
($v)=m/\b$k=([^;\s]+)/; $v||=0; $v=~s/%?,.*// if($k eq "cov"); push @at,"$k=$v"; } 
unless(($aw)=m/\baalen=(\d+),\w/){ ($aw)=m/\baalen=(\d+)/; }
($aa)=m/protein=([^;\s]+)/; ($xa)= $aa=~tr/Xx/Xx/; if($xa) { $at[0].=",gaps:$xa"; $aw -= $xa; }
if($pubv=$pubaq{$ids}) { my($d,$w,$g,$aq,$cl,$off)=@$pubv; 
$ao=$w; $at[1]="aaold=$aq"; $at[4]="clen=$cl"; $at[3]="offp=$off"; }
else { ($ao)=(m/\baaold=(\d+)/)?$1:$aw; }
unless(($sp)=m/(Split=[^;\s]+)/){ $sp=($id=~/_C(\d+)/)?"Split=$1":0; }  
map{ s/^\w+=//; } @at; 
print join("\t", qw(AQueryID diffa cds2alen pubalen cdsoff puboff pubclen cxlen cov pid path maptag))."\n" 
 unless($didhdr++);
print join("\t",$id,($aw - $ao).".da",@at,$sp,$STAG)."\n"; 
} elsif(/^\w+\t\d/) { ($td,$aw,$gp,$aq,$clen,$offs)= @v= split"\t"; 
if($gp>0 and not($aq =~ /gap/)) { $aq.=",gaps:$gp"; $v[3]=$aq; } 
$offs=~s/:.//; $v[5]=$offs;
$pubaw{$td}=$aw; $pubaq{$td}=[@v]; }' \
 > $pt.$suf

} done
fi 

#---------------------------------------

if [ $makesum = 1 ]; then

  nb=`/bin/ls -1 *.$suf | wc -l | sed 's/ *//;'`
  alts=t1; if [ $doalt = 1 ]; then alts=ta; fi

  sumf=kfish2$suf$alts.bestof$nb.tab
  echo "# makesum to $sumf"

  sort -k1,1 -k2,2nr -k10,10 *.$suf | env doalt=$doalt perl -ne \
'BEGIN{ $SLOP=$ENV{slop}||3; $DOALT=$ENV{doalt}||0; }
if(/^AQueryID/) { s/\t\w+$/\tgoodmap\tpoormap/; print unless($hdr++); next; }
chomp; @v=split"\t"; ($td,$alen)=@v[0,2]; if($td=~s/_C\d$//) { $v[0]=$td; }
if($DOALT) { next if($td =~ /t1$/); } else { next unless($td =~ /t1$/); }
if($td ne $ltd) { putb() if($ltd); @av=(); }
($aw)= $alen=~m/(\d+)/; if(($ag)= $alen=~m/gaps.(\d+)/){ $aw -= $ag; }
unshift(@v,$aw); push @av, [@v]; $ltd=$td; END { putb(); } 
sub putb { my($bv,$bvw,$av,$avw,$tg,@topt,@poort); 
@av=sort{ $$b[0] <=> $$a[0] } @av; $bv=$av[0]; $bvw= $$bv[0];
@ov= @{$bv}; shift(@ov); pop(@ov);
for $av (@av) { $tg= $$av[-1]; $avw= $$av[0];
if($avw + $SLOP < $bvw) { push @poort,$tg; } else { push @topt, $tg; } } 
@poort=qw(none) unless(@poort);
print join("\t" ,@ov, join(",",@topt), join(",",@poort))."\n";  }' \
  > $sumf

fi

#----------- findcds4 inputs -------------

if [ $findcds = 1 ]; then

  fcopt=" -debug -nostopcodon -full=0.75"
  fixopt=" -fixcds -debug -nostopcodon -full=0.75"
  fixvers="${fcvers}x"

  for ingff in in.*.gff*; do {
    ogff=`basename $ingff .gff.gz | sed 's/^in.//; s/\.gff.*//;'`; 
    ogff=$ogff$fcvers.gff
    gasm=$gasm3n 
    fasm2=`echo $ingff | sed 's/asm2.gff//;'`
    if [ $fasm2 != $ingff ]; then gasm=$gasm2a; fi
    echo $evigene/scripts/genefindcds2.pl $fcopt -dna $gasm -genes $ingff -outgff $ogff -outcdna -outaa -outcds
    $evigene/scripts/genefindcds2.pl $fcopt -dna $gasm -genes $ingff -outgff $ogff -outcdna -outaa -outcds

  } done

  for ingff in fixset/in.*.gff*; do {
    ogff=`basename $ingff .gff.gz | sed 's/^in.//; s/\.gff.*//;'`; 
    ogff=$ogff$fixvers.gff
    gasm=$gasm3n
    echo $evigene/scripts/genefindcds2.pl $fixopt -dna $gasm -genes $ingff -outgff fixset/$ogff -outcdna -outaa -outcds
    $evigene/scripts/genefindcds2.pl $fixopt -dna $gasm -genes $ingff -outgff fixset/$ogff -outcdna -outaa -outcds
  } done
  
fi

if [ $makeendgaps = 1 ]; then
  fixvers="${fcvers}x"
  #and fixset/*$fixvers.cdna*
  
  # for cdnaz in fixset/*$fixvers.cdna*; do 
  for cdnaz in *$fcvers.cdna* fixset/*$fixvers.cdna*; do 
  {
  pt=`echo $cdnaz | sed 's/.gz//;'`
  CAT="gunzip -c" ;  if [ $pt = $cdnaz ]; then CAT=cat; fi
  echo "#Make $pt.endgaps"
  
  $CAT $cdnaz | perl -ne'if(/^>(\S+)/){$d=$1; putn() if($id and $fa); $hd=$_; $id=$d; $fa=""; } 
elsif(/\w/) { chomp; $fa.=$_; } END{ putn(); } sub putn{ $nb=index($fa,"N"); $ne=rindex($fa,"N"); $w=length($fa); 
$gb=$ge=-1; $gb1=$ge0=0; if($nb>=0 and $nb<15) { if($nb==0 or substr($fa,$nb+1,1) eq "N") { $gb=$gb1=$nb; 
$gb1++ while(substr($fa,$gb1+1,1) eq "N"); } } if($ne > $w-15) { if($ne == $w-1 or substr($fa,$ne-1,1) eq "N") 
{ $ge=$ge0=$ne; $ge0-- while(substr($fa,$ge0-1,1) eq "N"); } } 
if($gb> -1 or $ge > -1) { $gb=($gb<0)?0:"$gb-$gb1"; $ge=($ge<0)?0:"$ge0-$ge"; 
($offs)= $hd=~m/cdsoff=([^;\s]+)/; $offs||=0; print join("\t",$id,$w,$offs,$gb,$ge)."\n";}} 
BEGIN{printf "%-16s","AGeneID"; print "\t".join("\t",qw(TrSize Cdsoff Gapbeg Gapend))."\n"; } ' \
  > $pt.endgaps
  } done

fi

exit;
#--------
## input.gff links: renamed for outputs; ** add fixset/ ** opt -fixcds
# in.kfish2nsplign15h.gff.gz@	in.kfish2rae5h_asm3n.gff.gz@	in.kfish2rae5htaasm2.gff@
# in.kfish2nsplign15n.gff.gz@	in.kfish2rae5ht1asm2.gff@
# 
# kfish2nsplign15h_fcds6.gff.gz@	kfish2rae5h_asm3n_fcds6.gff.gz	kfish2rae5htaasm2_fcds6.gff.gz
# kfish2nsplign15n_fcds6.gff.gz@	kfish2rae5ht1asm2_fcds6.gff.gz
## fixset/
# pt=kfish2nsplign15n;  $evigene/scripts/genefindcds2.pl -fixcds -nostopcodon -full=0.75 -dna $gasm3n -genes $pt.gff.gz -outgff $pt.fcds6x.gff -outcdna -outaa -outcds
# pt=kfish2nsplign15h;  $evigene/scripts/genefindcds2.pl -fixcds -nostopcodon -full=0.75 -dna $gasm3n -genes $pt.gff.gz -outgff $pt.fcds6x.gff -outcdna -outaa -outcds 
# pt=kfish2rae5h_gmap3n;  $evigene/scripts/genefindcds2.pl -fixcds -nostopcodon -full=0.75 -dna $gasm3n -genes $pt.gff.gz -outgff $pt.fcds6x.gff -outcdna -outaa -outcds

#......
# DROP -nogoodlen ..
# x11-gmap3n : same as kfish2rae5hpub-kfish2n.gmap/kfish2rae5h_asm3n_fcds4
# $evigene/scripts/genefindcds2.pl -debug -nostopcodon -nogoodlen -full=0.75 -dna $kfish2/submitf/pubgenome/ncbifunhe302scaf.fa -genes kfish2rae5x11t6_tsasubmitpub-kfish2n.gmap.gff.gz -outgff kfish2x11tsa_asm3n_fcds4.gff -outcdna -outaa -outcds &

# gmap/kfish2rae5h_asm3n
$evigene/scripts/genefindcds2.pl -debug -nostopcodon -nogoodlen -full=0.75 -dna $kfish2/submitf/pubgenome/ncbifunhe302scaf.fa -genes kfish2rae5hpub-kfish2n.gmap.gff.gz -outgff kfish2rae5h_asm3n_fcds4.gff -outcdna -outaa -outcds &

# gmap/kfish2rae5h_kfish2a main
$evigene/scripts/genefindcds2.pl -debug -nostopcodon -nogoodlen -full=0.75 -dna $kfish2/genome/killifish20130322asm.fa -genes kfish2rae5h.main.gff -outgff kfish2rae5h_kfish2a_t1fcds4.gff -outcdna -outaa -outcds &

# gmap/kfish2rae5h_kfish2a alts
$evigene/scripts/genefindcds2.pl -debug -nostopcodon -nogoodlen -full=0.75 -dna $kfish2/genome/killifish20130322asm.fa -genes kfish2rae5h.alt.gff -outgff kfish2rae5h_kfish2a_tafcds4.gff -outcdna -outaa -outcds &

# gsplign15b/
$evigene/scripts/genefindcds2.pl -debug -nostopcodon -full=0.75 -dna $kfish2/submitf/pubgenome/ncbifunhe302scaf.fa -genes kfish2nsplign15n.gff.gz -outgff kfish2nsplign15n_fcds4g.gff -outcdna -outaa -outcds &

# gsplign15/
$evigene/scripts/genefindcds2.pl -debug -nostopcodon -nogoodlen -full=0.75 -dna $kfish2/submitf/pubgenome/ncbifunhe302scaf.fa -genes kfish2nsplign15h.gff.gz -outgff kfish2nsplign15h_fcds4.gff -outcdna -outaa -outcds &

