#! /bin/bash
### env prodb=xxx protin=prot.aa datad=`pwd` qsub -q normal deblastref.sh
#PBS -N deblasta1 
#PBS -A ind114
#PBS -l nodes=1:ppn=32,walltime=28:55:00
#PBS -o deblasta1.$$.out
#PBS -e deblasta1.$$.err
#PBS -V

## * deltablast to replace blastp, and get cdd domains.
## * -q sppgenes.aa -db ref.aa -rspdb cdb/cdd_delta
## takes about 1.5 hr for cdd_delta + ref3all.aa x aphid2 aa (35kaa) x 16 parts,cpus
## but large 300k aa set, whitefly8.aa, takes 5-6hr
## .. too long for large prodb (250k UniRef50 vert..) + large 100k mRNA.asm okay set : < 1/4 done in 11hr
## .. reduce protdb **  best if cant reduce trasm set
## .. 100k of 250k have only 1 species ref, drop those as outliers?

ncpu=32
nbin=$HOME/bio/ncbi2228/bin
evigene=$HOME/bio/evigene

if [ "X" = "X$datad" ]; then
  datad=$HOME/scratchn/chrs/kfish/prot
fi

# NCBI blastdb/cdd_delta.tar.gz
cddb=cdb/cdd_delta
prodbz=$prodb
prodba=`echo $prodb | sed 's/.gz//;'`
pronam=`basename $prodba .aa`

odir=bldref
bltag=sd; 
blopt="-rpsdb $cddb -show_domain_hits -evalue 1e-5"

cd $datad/

if [ ! -f $prodba.psq ]; then 
   if [ $prodb = "$prodba.gz" ]; then
    gunzip -c $prodb | $nbin/makeblastdb -dbtype prot -title $prodba -out $prodba
   else
    $nbin/makeblastdb -dbtype prot -in $prodb; 
   fi
fi
prodb=$prodba

mkdir $odir

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
  echo $nbin/deltablast $blopt -outfmt 7 -db $prodb -query $qfile -out $onam.deblastp
  $nbin/deltablast $blopt -outfmt 7 -db $prodb -query $qfile -out $onam.deblastp  &
}

wait

opack=`echo $bltag-$pronam-$qnam | sed 's/.split.*//'`
cat $odir/$opack.*.deblastp > $opack.deblastp
gzip --fast $opack.deblastp
