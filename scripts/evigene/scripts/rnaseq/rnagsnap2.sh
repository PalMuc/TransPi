#! /bin/bash -l
### qsub -q batch rnagsnap2.sh
#PBS -N gsnap1
#PBS -A ind114
#PBS -l nodes=1:ppn=8,walltime=11:55:00
#PBS -o gsnap2.out
#PBS -e gsnap2.err
#PBS -V

workd=$HOME/scratch/chrs/daphmag
gmapd=$HOME/bio/gmap
dgenome=dmagna20100422assembly

drna=`basename $rnain .fa2`
#drna=16-1-R2

cd $workd/rnas/

if ! test -f "$drna.fa2" ; then  echo "missing $drna.fa2"; exit;  fi

for i in  8 9 10 11 12 13 14 15 ; {
 $gmapd/bin/gsnap -N 1 -k 14 --local-splice-penalty=1 --indel-penalty=1 --pairmax=5000 -A "sam" \
    --part=$i/16 -D $workd/genome/gmap -d $dgenome $drna.fa2  > $drna.gsnap$i.samu &

}

wait

