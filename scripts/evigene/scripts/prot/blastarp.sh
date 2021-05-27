#! /bin/bash
### env protin=prot.aa qsub -q normal blastarp.sh
#PBS -N blasta1 
#PBS -l nodes=1:ppn=32,walltime=4:55:00
#PBS -o blast1.$$.out
#PBS -e blast1.$$.err
#PBS -V

ncpu=32

nbin=/home/diag/opt/blast/2.2.24/bin/
workd=$HOME/scratch/chrs/aphid2

prodb=arp5hum.aa
pronam=`basename $prodb .aa`

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
  onam=bl-$pronam-$qnam
  echo $nbin/blastp -outfmt 7 -evalue 1e-5 -db $prodb -query $qfile -out $onam.blastp
  $nbin/blastp -outfmt 7 -evalue 1e-5 -db $prodb -query $qfile -out $onam.blastp  &
}

wait

