#! /bin/bash
### env protin=prot.aa qsub -q normal blast2genome.sh
#PBS -N tblastn1 
#PBS -A ddp138
#PBS -l nodes=2:ppn=32,walltime=39:55:00
#PBS -o tblastn1.$$.out
#PBS -e tblastn1.$$.err
#PBS -V

#* FIXME: general problem w/ sd-ssd disk use: fails return no info
#  .. including over time limit when partial results can be useful
#  .. add forked shell sub for periodic (sleep 1hr? 30min?) update on outputs: ls, copy home?
#  .. do copy per fork: ( blast ... -out xx1 ; cp -p xx1 $home/ ) &
#  plant8.aa x cacao11 takes > 32cpu x 22hr wall ; retry 64cpu x 40hr split

ncpu=64

nbin=/home/diag/opt/blast/2.2.24/bin/
workd=/oasis/$USER/chrs/cacao
rund=/scratch/$USER/$PBS_JOBID
outdir=$workd/prot

dgenome=cacao11allasm_repmask
dnam=cacao11all
notef=$workd/prot/$dnam.$$.RUNNING


# want tblastn opt -lcase_masking
blopt="-lcase_masking -outfmt 7 -evalue 1e-5"

touch $notef
mkdir -p $rund
cp -p $workd/prot/$protin $rund/
cp -p $workd/genome/$dgenome.n?? $rund/
prots=`basename $protin`

# $nbin/makeblastdb -in $dgenome.fa -dbtype nucl  -out $dgenome

# cd $workd/prot/
cd $rund/
echo "START blast2genome $prots" >> $notef
echo `date`  >> $notef
du -h >> $notef
ls -l >> $notef

# if [ ! -f $prots.split.1.fa ]; then
 splitsize=`wc -c $prots | sed 's/ .*//' `
 splitsize=`echo $splitsize / $ncpu | bc`
 $workd/scripts/splitMfasta.pl --minsize=$splitsize $prots
# fi

qset=`/bin/ls $prots.split.*.fa`
for qfile in $qset
{
  qnam=`basename $qfile .fa`
  onam=bl-$dnam-$qnam
  echo $nbin/tblastn $blopt -db $dgenome -query $qfile -out $onam.tblastn >> $notef
  ( $nbin/tblastn $blopt -db $dgenome  -query $qfile -out $onam.tblastn ; cp -p $onam.tblastn $outdir/ ) &
}

wait

#done# cp -p *.tblastn $outdir/
du -h >> $notef
ls -l >> $notef

echo 'DONE ' >> $notef
echo `date`  >> $notef
