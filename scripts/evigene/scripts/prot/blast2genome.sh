#! /bin/bash
### env protin=prot.aa qsub -q normal blastgenome.sh
#PBS -N tblastn1 
#PBS -A ddp138
#PBS -l nodes=1:ppn=32,walltime=23:55:00
#PBS -o blast1.$$.out
#PBS -e blast1.$$.err
#PBS -V

ncpu=32

nbin=/home/diag/opt/blast/2.2.24/bin/
workd=$HOME/scratch/chrs/nasv1
dgenome=nasvit1asm
dgsize=$workd/genome/$dgenome.chr_size.txt

db=$workd/genome/$dgenome
dnam=$dgenome

# want tblastn opt -lcase_masking
blopt="-lcase_masking -outfmt 7 -evalue 1e-5"

cd $workd/prot/

if [ ! -f $protin.split.1.fa ]; then
 splitsize=`wc -c $protin | sed 's/ .*//' `
 splitsize=`echo $splitsize / $ncpu | bc`
 $workd/scripts/splitMfasta.pl --minsize=$splitsize $protin
fi

qset=`/bin/ls $protin.split.*.fa`

for qfile in $qset
{
  qnam=`basename $qfile .fa`
  onam=bl-$dnam-$qnam
  echo $nbin/tblastn $blopt -db $db -query $qfile -out $onam.tblastn
  $nbin/tblastn $blopt -db $db  -query $qfile -out $onam.tblastn &
}

wait

