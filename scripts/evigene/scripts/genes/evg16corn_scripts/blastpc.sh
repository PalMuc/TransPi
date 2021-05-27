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

pronam=`basename $prodb .aa`
bltag=sdc; blopt="-evalue 1e-5"
if [ "X" = "X$name" ]; then name=`basename $protin .aa`; fi;
odir=bloutc1$name

cd $datad
if [ "X" = "X$protin" ]; then echo "ERR: missing protin=what?"; exit -1; fi
if [ "X" = "X$prodb" ]; then echo "ERR: missing prodb=$prodb"; exit -1; fi
#x if [ ! -f $prodb.psq ]; then echo "ERR: missing prodb=$prodb"; exit -1; fi
mkdir $odir
echo "#START: blastp $protin `date`"

if [ ! -f $protin.split.1.fa ]; then
 pindir=`dirname $protin`
 # dgg: protin.aa.gz allowed now ? need splitsize calc : updated splitMfasta.pl
 $evigene/scripts/splitMfasta.pl --outputpath=$pindir --nparts=$ncpu $protin
 ## old way
 # splitsize=`grep -v '^>' $protin | wc -c | sed 's/ .*//' `
 # splitbp=$(( $splitsize / $ncpu ))
 # $evigene/scripts/splitMfasta.pl --outputpath=$pindir --minsize=$splitbp $protin
fi

qset=`/bin/ls $protin.split.*.fa`

for qfile in $qset
{
  qnam=`basename $qfile .fa`
  onam=$odir/$bltag-$pronam-$qnam
  echo $nbin/blastp $blopt -outfmt 7 -db $prodb -query $qfile -out $onam.blastp
  $nbin/blastp $blopt -outfmt 7 -db $prodb -query $qfile -out $onam.blastp  &
}

wait

opack=`echo $bltag-$pronam-$qnam | sed 's/.split.*//'`
cat $odir/$opack.*.blastp > $opack.blastp
gzip --fast $opack.blastp
echo "#DONE: blastp $opack  `date`"

