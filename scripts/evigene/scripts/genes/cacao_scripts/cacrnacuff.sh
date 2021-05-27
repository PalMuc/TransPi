#! /bin/bash -l
### qsub -q batch cacrnagsnap1.sh
#PBS -N cuff2
#PBS -A ind114
#PBS -l nodes=1:ppn=2,walltime=22:55:00
#PBS -o cuff2.out
#PBS -e cuff2.err
#PBS -V

workd=$HOME/scratch/chrs/cacao
bin_dir=$HOME/bio/tophat/bin

dgenome=cacao9asm
#drna=Matina_1
drna=Matina_2

## prepare:
# cat $drna.gsnap*.samu | grep -v '^@SQ' | sort -T /scratch/dgilbert -k 3,3 -k 4,4n > $drna.gsnap.sam
## try merged set:
#  sort -m -T /scratch/dgilbert -k 3,3 -k 4,4n Matina_1.gsnap.sam Matina_2.gsnap.sam > Matina_m.gsnap.sam
drna=Matina_m

opts="--inner-dist-mean 100"
sam=$drna.gsnap.sam

cd $workd/rnas/

## dang cant run 2 at once, writes same file names ... subdir
mkdir cuff.$drna
cd cuff.$drna

$bin_dir/cufflinks $opts ../$sam > log.cuf$drna

# wait

