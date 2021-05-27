#! /bin/bash -l
### env estin=xxxxx qsub -q batch estgmap1.sh
#PBS -N estgmap1b
#PBS -A ind114
#PBS -l nodes=1:ppn=8,walltime=11:55:00
#PBS -o estgmap1b.out
#PBS -e estgmap1b.err
#PBS -V

workd=$HOME/scratch/chrs/aphid2
gmapd=$HOME/bio/gmap

#ests
opts="-n 4 -S"
# assemblies
opts="-n 50 -f 2"

dest=$estin
dgenome=aphid2asm
#dest=acyr_est201009
#dest=agossypi_est
#dest=mpersica_est
#dest=soyaphid_est

cd $workd/est/
if ! test -f $dest.fa ; then  echo "missing $dest.fa"; exit;  fi

for i in  0 1 2 3 4 5 6 7 ; 
{
  $gmapd/bin/gmap $opts -D $workd/genome/gmap -d $dgenome --jobdiv=$i/16 $dest.fa > $dest.gmap$i.out &
}

wait


