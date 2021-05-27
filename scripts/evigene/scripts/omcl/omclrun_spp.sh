#! /bin/bash
### env sppall='APISU,TCAST,DMELA,daphplx,HSAPI' sppone='AMELL CFLOR HSALT wasp bombusimp' \
###    fasta=myspp.aa blastp=myspp.blastp.gz datad=`pwd` qsub -q normal omclrun_spp.sh
#PBS -N omclrun_spp
#PBS -A ind114
#PBS -l nodes=1:ppn=16,walltime=47:55:00
#PBS -V

ncpu=16

# datad=$HOME/scratchg/chrs/nasv1/prot
# blin=sd-hym11set12b-hym11set12b.aa.blastp.gz
# fasta=aa12b/hym11set12b.aa
if [ "X" = "X$datad" ]; then echo "ERR: env datad=path/to/data"; exit -1; fi
if [ "X" = "X$blastp" ]; then echo "ERR: env blastp=myspecies.blastp.gz"; exit -1; fi
if [ "X" = "X$fasta" ]; then echo "ERR: env fasta=myspecies.prot.fasta"; exit -1; fi

## should be options
# sppall='APISU,TCAST,DMELA,daphplx,HSAPI,DRERI'
# sppone='AECHI AMELL CFLOR HSALT wasp bombusimp bombusterr ACEPH LHUMI PBARB SINVI'
if [ "X" = "X$sppall" ]; then echo "ERR: env sppall='APISU,TCAST,DMELA,daphplx,HSAPI,DRERI'"; exit -1; fi
if [ "X" = "X$sppone" ]; then echo "ERR: env sppone='AMELL CFLOR HSALT wasp bombusimp bombusterr ACEPH LHUMI'"; exit -1; fi

evigene=$HOME/bio/evigene
orlib=$HOME/bio/orthomcl3
mcl=$HOME/bio/mcl9/bin/mcl

domake=1
#orig#bopts=""
bopts="-identmin 33"
onam=`basename $fasta .fasta | sed 's/\.aa.*//; s/\.pep.*//; s/\.fa.*//; s/$/_omcl/;'`
odirbase="${onam}_f$$"

echo "START omcl_spp `date` "
i=0; 
for sp1 in $sppone; do {

  orgs="$sp1,$sppall"
  odir="$odirbase$sp1"
  bopt1="$bopts -organism=$orgs"
  
  export MCL=$mcl 
  export ORTHOMCL=$datad/$odir/
  
  cd $datad
  mkdir $odir

  echo "# $odir: blast92orthomcl10.pl $bopt1 -fasta=$fasta -in=$blastp -out=$odir/$onam ; perl -I$orlib $orlib/orthomcl.pl --mode 4 --bpo=$onam.bpo --gg=$onam.gg"
 
  ( $evigene/scripts/blast92orthomcl10.pl $bopt1  -fasta=$fasta -in=$blastp -out=$odir/$onam; \
    cd $datad/$odir/; perl -I$orlib $orlib/orthomcl.pl --mode 4 --bpo=$onam.bpo --gg=$onam.gg ) &
    
  i=$(( $i + 1 ))
  if [ $i -ge $ncpu ]; then wait; i=0; fi

} done
wait

echo "DONE omcl_spp : `date`" 

