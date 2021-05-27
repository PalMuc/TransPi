#! /bin/bash
## env prog=blastallc.sh prodb=xxx protin=prot.aa datad=`pwd` sbatch srun_prog.sh
## comet.sdsc sbatch update; show_accounts ?? not ind114
## .. redo as 2 script set as aprun/bigred hack, call w/ 2nd script param
## this is 2nd script blastallc.sh
#... old pbs ..
## env prodb=xxx protin=prot.aa datad=`pwd` qsub -q normal blastallg.sh
#..PBS -N blasta1 
#..PBS -A ind114
#..PBS -l nodes=4:ppn=16,walltime=47:55:00
#..PBS -V

#env caller sets  ncpu=24

nbin=$HOME/bio/ncbi2230/bin
evigene=$HOME/bio/evigene

if [ "X" = "X$datad" ]; then echo "ERR:datad=what?"; exit -1; fi
if [ "X" = "X$ncpu" ]; then echo "ERR:ncpu=what?"; exit -1; fi

## selfblast uses ident=98 
## lower for cross-spp cds blastn
CDSBLAST_IDENT=65
CDSBLAST_EVALUE=1e-19

pronam=`basename $prodb .aa`
bltag=sdcx; 
blopt="-evalue 1e-5"
selfblopt="-task megablast -ungapped -xdrop_ungap 4 -dust no -perc_identity $CDSBLAST_IDENT -evalue $CDSBLAST_EVALUE "; 
blnopt="-task megablast -perc_identity $CDSBLAST_IDENT -evalue $CDSBLAST_EVALUE "; 

if [ "X" = "X$name" ]; then name=`basename $protin .tr | sed "s/\.[a-z]*//;"`; fi;
odir=bloutc1$name

cd $datad
if [ "X" = "X$protin" ]; then echo "ERR: missing protin=what?"; exit -1; fi
if [ "X" = "X$prodb" ]; then echo "ERR: missing prodb=$prodb"; exit -1; fi
#x if [ ! -f $prodb.psq ]; then echo "ERR: missing prodb=$prodb"; exit -1; fi
mkdir $odir
echo "#START: blastn $protin `date`"

if [ ! -f $protin.split.1.fa ]; then
 pindir=`dirname $protin`
 splitsize=`grep -v '^>' $protin | wc -c | sed 's/ .*//' `
 splitbp=$(( $splitsize / $ncpu ))
 $evigene/scripts/splitMfasta.pl --outputpath=$pindir --minsize=$splitbp $protin
fi

qset=`/bin/ls $protin.split.*.fa`

for qfile in $qset
{
  qnam=`basename $qfile .fa`
  onam=$odir/$bltag-$pronam-$qnam
  echo $nbin/blastn $blnopt -outfmt 7 -db $prodb -query $qfile -out $onam.blastn
  $nbin/blastn $blnopt -outfmt 7 -db $prodb -query $qfile -out $onam.blastn &
}

wait

opack=`echo $bltag-$pronam-$qnam | sed 's/.split.*//'`
cat $odir/$opack.*.blastn > $opack.blastn
gzip --fast $opack.blastn
echo "#DONE: blastn $opack  `date`"

