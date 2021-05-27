#! /bin/bash
## blast2genome_msk.sh
### env protin=prot.aa qsub -q normal blast2genome_msk.sh
#PBS -N tblastn1 
#PBS -A ind114
#PBS -l nodes=1:ppn=32,walltime=39:55:00
#PBS -o tblastn1.$$.out
#PBS -e tblastn1.$$.err
#PBS -V

ncpu=32
# runtime ~6h x 32 cpu for fish prot x kfish2asm

nbin=$HOME/bio/ncbi2228/bin
evigene=$HOME/bio/evigene

# shoud be params:
workd=$HOME/scratchn/chrs/kfish
dgenome=killifish2asm_repmask; dnam=kfish2brm

# # tblastn opt -lcase_masking
blopt="-lcase_masking -outfmt 7 -evalue 1e-5"
# blopt="-outfmt 7 -evalue 1e-5"

cd $workd/prot/
outdir=run$dnam
mkdir $outdir

echo "START blast2genome $prots `date`" 

prots=$protin
if [ ! -f $prots.split.1.fa ]; then
 pindir=`dirname $prots`
 splitsize=`grep -v '^>' $prots | wc -c | sed 's/ .*//' `
 splitbp=$(( $splitsize / $ncpu ))
 $evigene/scripts/splitMfasta.pl --outputpath=$pindir --minsize=$splitbp $prots
fi

qset=`/bin/ls $prots.split.*.fa`
for qfile in $qset
{
  qnam=`basename $qfile .fa`
  onam=bl-$dnam-$qnam
  echo $nbin/tblastn $blopt -db $dgenome -query $qfile -out $outdir/$onam.tblastn
  $nbin/tblastn $blopt -db $workd/genome/$dgenome -query $qfile -out $outdir/$onam.tblastn &
}

wait

protf=`basename $prots .aa`
outf="sd-$dnam-$protf.tblastn"
cat $outdir/bl-$dnam*$protf*.tblastn > $outf
gzip --fast $outf

echo "DONE `date`" 
