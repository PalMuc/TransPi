#! /bin/bash -l
### qsub -q batch rnacuff.sh
#PBS -N rcuff1
#PBS -A ind114
#PBS -l nodes=1:ppn=2,walltime=22:55:00
#PBS -o cuff1.out
#PBS -e cuff1.err
#PBS -V

workd=$HOME/scratch/chrs/daphmag
bin_dir=$HOME/bio/tophat/bin
dgenome=dmagna20100422assembly
drna=`basename $rnain .gsnap.sam`

## prepare:
# cat $drna.gsnap*.samu | grep -v '^@SQ' | sort -T /scratch/dgilbert -k 3,3 -k 4,4n > $drna.gsnap.sam
## merged set:
#  sort -m -T /scratch/dgilbert -k 3,3 -k 4,4n Matina_1.gsnap.sam Matina_2.gsnap.sam > Matina_m.gsnap.sam

opts="--inner-dist-mean 100"
sam=$drna.gsnap.sam

cd $workd/rnas/

if ! test -f "$sam" ; then  echo "missing $sam"; exit;  fi

mkdir cuff.$drna
cd cuff.$drna

$bin_dir/cufflinks $opts ../$sam  > log.cuf$drna 2>&1

