#!/bin/bash
# evgmrna2asn.sh/evgmrnakfish2.sh

if [ "X" = "X$datad" ]; then echo "missing env datad=path/to/data"; exit -1; fi
if [ "X" = "X$mrna" ]; then "echo env mrna=path/to/name.mrna"; exit -1; fi
cd $datad/
namepath=`echo $mrna | sed 's/\.mrna//; s/\.fasta//; s/\.fa//; s/\.tr//;'`
namebase=`basename $namepath`

evigene=/bio/bio-grid/mb/evigene
export evigenes=$evigene/scripts/
nbin=/bio/bio-grid/mb/ncbic/bin
export vecscreen=$nbin/vecscreen
export tbl2asn=$nbin/tbl2asn

# fixme: should be in name.info, option or look for, not here..
# fixme2: name.mrna2tsa.info has right val in m2t.TSADESC= but ignored !!!
export tsasubmit_sbt=submitset/evgr_tsasubmit.sbt

ncpu=3
trimopts="-nodeferupdate -NCPU $ncpu -log"
tblopts="-onlysubmit -runtbl2asn -NCPU 1 -log -debug"
# tidier submit package w/ ncpu=1, run a few same time on dgmac

echo "# Step1: " $evigene/scripts/rnaseq/asmrna_trimvec.pl $trimopts -mrna $mrna
if [ ! -f $namepath.trimvec_done ]; then
  $evigene/scripts/rnaseq/asmrna_trimvec.pl $trimopts -mrna $mrna
  if [ -f publicset/$namebase.mrna ]; then mrna=publicset/$namebase.mrna; fi
else
  echo "# Step1: $namepath.trimvec_done"
fi

# IF ok...
echo "# Step2: " $evigene/scripts/evgmrna2tsa2.pl $tblopts -mrna $mrna 
if [ -f $namepath.trimvec_done ]; then
  $evigene/scripts/evgmrna2tsa2.pl $tblopts -mrna $mrna 
else
  echo "# Step1: not done, $namepath.trimvec_done"
fi

