#! /bin/bash
### env datad=path/to/data trset=myspecies_all.tr.gz qsub -q normal runtr2cds.sh
#PBS -N tr2cds
#PBS -A ind114
#PBS -l nodes=1:ppn=16,walltime=47:55:00
#PBS -V

evigene=$HOME/bio/evigene
# blastp/n:
export PATH=$HOME/bio/ncbi2230/bin:$PATH

if [ "X" = "X$ncpu" ]; then ncpu=15; fi
if [ "X" = "X$maxmem" ]; then maxmem=50000; fi
if [ "X" = "X$aaset" ]; then echo "missing env aaset=xxxx.aa"; exit -1; fi
if [ "X" = "X$refaa" ]; then echo "missing env refaa=refaa_blastdb"; exit -1; fi
if [ "X" = "X$datad" ]; then echo "missing env datad=/path/to/data "; exit -1; fi

cd $datad/

echo "#START $aaset : `date`"

# pt=`basename $trset .tr | sed 's/\.fa.*//; s/\.tr.*//; s/\.cdna/;'`
# cat okayset/$pt.okay.aa okayset/$pt.okalt.aa > $qname.aa

# env protin=$pt.aa prodb=$refaa datad=`pwd` prog=./blastpc.sh sbatch ./srun_comet.sh
#.......
qname=`basename $aaset .aa`
refnam=`basename $refaa .aa`

blopt="-evalue 1e-5"
odir=blout1$qname
mkdir $odir

if [ ! -f $aaset.split.1.fa ]; then
 pindir=`dirname $aaset`
 splitsize=`grep -v '^>' $aaset | wc -c | sed 's/ .*//' `
 splitbp=$(( $splitsize / $ncpu ))
 $evigene/scripts/splitMfasta.pl --outputpath=$pindir --nparts $ncpu --minsize=$splitbp $aaset
fi

qset=`/bin/ls $aaset.split.*.fa`

for qfile in $qset
{
  qnamspl=`basename $qfile .fa`
  onam=$odir/$refnam-$qnamspl
  echo blastp $blopt -outfmt 7 -db $refaa -query $qfile -out $onam.blastp
  blastp $blopt -outfmt 7 -db $refaa -query $qfile -out $onam.blastp  &
}

wait

opack=$refnam-$qname
aablast=$refnam-$qname.blastp
cat $odir/$opack.*.blastp > $aablast
/bin/rm $qset

if [ ! -s $aaset.qual ]; then
env oid=1 off=1 $evigene/scripts/prot/aaqual.sh $aaset
fi

mbaopts="-tall -aasize $aaset.qual,$refaa.qual"
aabltab=refaa-$qname.aa.btall

$evigene/scripts/makeblastscore3.pl $mbaopts $aablast > $aabltab 

if [ -f $refaa.names ]; then
$evigene/scripts/prot/namegenes.pl -blast $aabltab -refnames $refaa.names -out $qname.names
fi

gzip --fast $aablast

echo "#DONE $trset : `date`"

