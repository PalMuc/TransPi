#! /bin/bash
### env protin=prot.aa qsub -q normal blastself.sh
#PBS -N blasts1 
#PBS -l nodes=1:ppn=32,walltime=4:55:00
#PBS -o blasts1.$$.out
#PBS -e blasts1.$$.err
#PBS -V

ncpu=32

nbin=/home/diag/opt/blast/2.2.24/bin/
workd=$HOME/scratch/chrs/aphid2

prodb=$protin
pronam=`basename $prodb .aa`

cd $workd/prot/

if [ ! -f $protin.psq ]; then
  $nbin/makeblastdb -dbtype prot -in $protin
fi

if [ ! -f $protin.split.1.fa ]; then
 splitsize=`wc -c $protin | sed 's/ .*//' `
 splitsize=`echo $splitsize / $ncpu | bc`
 $workd/scripts/splitMfasta.pl --minsize=$splitsize $protin
fi

qset=`/bin/ls $protin.split.*.fa`

for qfile in $qset
{
  qnam=`basename $qfile .fa`
  onam=bl-$pronam-$qnam
  echo $nbin/blastp -outfmt 7 -evalue 1e-5 -db $prodb -query $qfile -out $onam.blastp
  $nbin/blastp -outfmt 7 -evalue 1e-5 -db $prodb -query $qfile -out $onam.blastp  &
}

wait

