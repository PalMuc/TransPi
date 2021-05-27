#!/bin/bash
# processtr.sh  for soap,velvet,trinity .. with variants

bioapp=/bio/bio-grid/mb
if [ "X$evigene" = "X" ]; then
evigene=$bioapp/evigene
fi

fastanrdb=`which fastanrdb`
if [ "X$fastanrdb" = "X" ]; then
fastanrdb=$bioapp/exonerate/bin/fastanrdb
fi

function usage() {
  echo 'processtr.sh : process mRNA assembly transcripts to proteins'
  echo 'usage: processtr.sh trfile  speciestag '
  echo '  trfile = transcripts.fa[.gz] OR trdir/xxx [for velvet/soap/trinity dirs]'
  echo '  speciestag add to transcript ids for trdir/ '
  echo 'part of http://eugenes.org/EvidentialGene/'
  exit -1
}

if [ ! -x $evigene/scripts/rnaseq/maketraa.sh ]; then
  echo "ERR: need path to evigene/scripts/rnaseq/maketraa.sh "; usage;
fi

trdir=$1
spt=$2

# FIXME: collect these asmbly reformatters into one perl, detect input format, write common output form
#.. assembly.tr to trs/commonform.tr
## replace w/ $evigene//scripts/rnaseq/trformat.pl
if [ -d $trdir ]; then 
  if [ "X" = "X$spt" ]; then usage; fi
  $evigene/scripts/rnaseq/veltrmake.sh $trdir  $spt
  # $evigene/scripts/rnaseq/maketrsoap.sh $trdir  $spt
  # $evigene/scripts/rnaseq/maketrtrin.sh $trdir  $spt
  
  trfile=`ls $spt*.tr`
  if [ ! -s $trfile ]; then  echo "missing veltrmake.sh output $spt *.tr"; exit -1; fi
  trset=`echo $trfile | sed 's/\.tr//'`
  
elif [ -f $trdir ]; then
  trfile=$trdir
  trset=`echo $trfile | sed 's/\.tr//; s/\.fasta//; s/\.fa//;'` 

else
  usage;
fi


trnogz=`echo $trfile | sed 's/\.gz//'`; 
# TCAT=cat; if [ $trfile != $trnogz ]; then TCAT="gunzip -c"; fi

#.. remove identicals; typically < 5% even for multi-kmer tr sets.
if [ -x $fastanrdb -a $trfile = $trnogz ]; then 
  $fastanrdb $trfile > ${trset}nr.tr
  if [ -s ${trset}nr.tr ]; then trfile=${trset}nr.tr; trset=${trset}nr; fi
fi

#.. protein work
$evigene/scripts/rnaseq/maketraa.sh  $trfile 

## outputs: {trset}_allcd.aa, _allcd.aa.qual, _allcd.ids,  _allcd.tr
## collate best aa,tr,ids : now in maketraa.sh 

## for gff input only (cufflinks,..) add aa-annot to mrna
if [ -f $trset.gff -a -f ${trset}_allcd.aa ]; then 
  grep '^>' ${trset}_allcd.aa | cat - $trset.gff | perl -ne \
'if(/^>(\S+) +(\w.*)/) { $d=$1; $an=$2; $an=~s/ +//g; $an{$d}=$an; } else { \
if(/^\W/) { print; next; } if(/\tmRNA/) { ($d)=m/ID=([^;\s]+)/; $ok=0; if($an=$an{$d}) {\
$ok=1; s/;?$/;$an/; } } print if $ok; }'\
   > ${trset}_allcd.gff
fi

##.. tabulate  prots; do for full.aa and {ok,poor}_cd.aa
##.. maketraa.sh does this now via aaqual; drop this stat sum
# if [ -f ${trset}_allcd.aa.qual ]; then 
#   echo "# SUMMARY $trset : Top1k AA"
#   cat ${trset}_allcd.aa.qual | egrep -v '^#|^total' | sort -k2,2nr | head -1000 | env nam=$trset perl -ne \
#   '($aw,$nn)=(split)[1,2]; $n++; $sw+=$aw; $sn+=$nn; push @aw,$aw;
#   END{ $aw=int($sw/$n); $an=int(10*$sn/$n)/10; @aw=sort{$b <=> $a}@aw ; ($mx,$md,$mi)=@aw[0,int($n/2),-1];
#   print "# $ENV{nam}\t  n=$n; aw=$aw; med=$md; min,max=$mi,$mx; sw=$sw; sn=$sn,$an\n"; }'
# fi
##.. drop this kmer stat also?
# if [ -f ${trset}_allcd.aa.gz ]; then 
#   echo "# SUMMARY $trset : kmer qualities"
#   gunzip -c ${trset}_allcd.aa.gz | grep  '^>' | perl -ne \
# '($d)=m/>(\w+)/; ($al,$pc,$ac,$ut)=m/aalen=(\d+).(\d+)..(\w+)[-]?(\w*)/; ($cl)=m/clen=(\d+)/; 
# unless( ($k)=$d=~m/k(\d+)[Ll]oc/ ) { ($k)= $d =~ m/sub(\d+)[Ll]oc/; $k||=0; } $ut||="okutr"; 
# for $k ($k,"all") { $ku{$k}{Tot}++;  $ku{$k}{$ut}++; $kul{$k}{$ut}+= $al; } $ut{$ut}++; 
# END{ @ut=sort keys %ut; print "# ",join("\t","kmer","total",@ut)."\n"; foreach $k (sort keys %ku) { 
# $tc=$ku{$k}{Tot}; print "# $k\t$tc"; foreach $u (@ut) { $sw=$kul{$k}{$u}; $c=$ku{$k}{$u}||1; 
# $aw=int($sw/$c); $pc=int(100*$c/$tc);  print "\t$c,$pc%,${aw}aa"; } print "\n"; } }' 
# 
# fi

