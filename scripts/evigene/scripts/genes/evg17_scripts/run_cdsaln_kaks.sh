#! /bin/bash
### env dpt=subject.cds qpt=query.cds datad=path/to/data qsub -q normal runtr2cds.sh
#PBS -N tr2cds
#PBS -A ind114
#PBS -l nodes=1:ppn=16,walltime=47:55:00
#PBS -V

## add to path?
kkbin=/bio/bio-grid/tmpd/gnohist/kaks_calculator/bin
# self cdsaln

if [ "X" = "X$dpt" ]; then
  if [ "X" != "X$trset" ]; then 
    dpt=`echo $trset | sed 's/.cds//;'`; qpt=$dpt; 
  else 
    echo missing dpt=subject.cds qpt=query.cds; exit -1; 
  fi
else
  dpt=`echo $dpt | sed 's/.cds//;'`; 
  qpt=`echo $qpt | sed 's/.cds//;'`; 
fi

# blastcds2axt opt: env oneonly=1 .. limit to 1/best pair align per query
axtopt=""; psuf=""
if [ "X" != "X$one" ]; then axtopt="oneonly=1"; psuf=".one"; fi

if [ "X" = "X$ncpu" ]; then ncpu=1; fi
if [ "X" = "X$maxmem" ]; then maxmem=5000; fi

if [ "X" = "X$evigene" ]; then
   if [ -d $HOME/bio/evigene ]; then evigene=$HOME/bio/evigene; fi
   if [ -d /bio/mb/evigene ]; then evigene=/bio/mb/evigene; fi
   if [ "X" = "X$evigene" ]; then echo missing evigene=/path/to/evigene; exit -1; fi
fi
# if [ "X" = "X$trset" ]; then echo "missing env trset=xxxx.tr"; exit -1; fi
if [ "X" = "X$datad" ]; then echo "missing env datad=/path/to/data "; exit -1; fi

# pn=fastanrdb;  pp=`which $pn`; pp=`basename $pp`; if [ "$pp" != "$pn" ]; then echo missing path to $pn; exit -1; fi
# pn=cd-hit-est; pp=`which $pn`; pp=`basename $pp`; if [ "$pp" != "$pn" ]; then echo missing path to $pn; exit -1; fi
pn=blastn;  pp=`which $pn`; pp=`basename $pp`; if [ "$pp" != "$pn" ]; then echo missing path to $pn; exit -1; fi
# pn=KaKs_Calculator;  pp=`which $pn`; pp=`basename $pp`; if [ "$pp" != "$pn" ]; then echo missing path to $pn; exit -1; fi

cd $datad/

echo "#START $trset : `date`"

if [ ! -f $dpt.cds.nsq ]; then
  makeblastdb -dbtype nucl -in $dpt.cds -logfile /dev/null
fi

if [ ! -f $qpt-$dpt.cds.dcblast1n ]; then
blastn -evalue 1e-9 -task dc-megablast -template_type coding -template_length 18 \
 -db $dpt.cds -query $qpt.cds -out $qpt-$dpt.cds.dcblast1n 
fi

ptin=$qpt-$dpt
pt=$qpt-$dpt$psuf
cat $ptin.cds.dcblast1n | env nop=0 $axtopt $evigene/scripts/prot/blastcds2axt.pl > $pt.axtc

$kkbin/KaKs_Calculator -m MYN -i $pt.axtc -o $pt.kaks > $pt.klog 2>&1

cat $pt.kaks | perl -pe \
'($sd)=@v=split"\t"; if($sd=~/\slen=\d/){ ($td)= $sd=~/^([\w\.-]+)/; ($lw)= $sd=~m/len=(\d+)/; 
($qd)= $sd=~m,/([\w\.-]+),; ($aw)= $sd=~m/aln=(\d+)/; $tda="$td,$aw/$lw,$qd"; }
elsif($v[1]=~/^len=\d/){ ($td,$awl,$qd)=@v[0,1,2]; ($lw,$aw)= $awl=~m/len=(\d+).aln=(\d+)/;
$qd=~s,^[\d:]+/,,; $qd=~s/:\d.*$//; $tda="$td,$aw/$lw,$qd"; splice(@v,1,2); }
if($lw and $qd) { $v[0]=$tda; } else { $v[0]= "0$sd"; } $_=join"\t",@v; ' \
 | cut -f1,3,4,5,6 | grep -v 'NA' > $pt.sigkaks.tab

echo "#DONE $trset : `date`"

#... freq distr for kaks.tab
# set pt=tsaGDQRpinealb-tsaGBYRok.one
# env kamax=9990.2 all=1 perl -ne 'BEGIN{ $KAMAX=$ENV{kamax}||999; $TOP=$ENV{top}||0; $ALL=$ENV{all}||0; } ($d)=@v=split; if(@v==1) { $nok++; $ok{$d}=1; } else { ($davdb,$ka,$ks,$kas,$pv)=@v; next if($ka > $KAMAX or $ks=~/NaN/); ($da,$al,$db)=split",",$davdb; if(($ALL or not $nok or ($ok{$da} and $ok{$db}))){  $rks=int($ks*10)/10; $fks{$rks}++; $sks+=$ks; push @ks,$ks;  $nt++; } } END{ @fks=sort{$a <=> $b} keys %fks; @sks=sort{$a <=> $b}@ks; $mdks=$sks[int($nt/2)]; $aks=sprintf "%.3f",$sks/$nt; print "nt=$nt; aveks=$aks, medianks=$mdks \n",join("\t",qw(Ks pKs nKs))."\n"; for $ks (@fks) { $c=$fks{$ks}; $p=int(0.5 + 1000*$c/$nt)/10; print "$ks\t$p\t$c\n"; } } ' \
#  sp2pairs.tabb.commidsNOT $pt.sigkaks.tab | head -25 ; echo $pt

