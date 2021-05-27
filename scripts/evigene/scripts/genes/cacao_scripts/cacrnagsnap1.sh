#! /bin/bash -l
### qsub -q batch cacrnagsnap1.sh
#PBS -N gsnap1
#PBS -A ind114
#PBS -l nodes=1:ppn=8,walltime=11:55:00
#PBS -o gsnap1.out
#PBS -e gsnap1.err
#PBS -V

workd=$HOME/scratch/chrs/cacao
gmapd=$HOME/bio/gmap

dgenome=cacao9asm
drna=Matina_1
#drna=Matina_2 #done

cd $workd/rnas/

for i in  0 1 2 3 4 5 6 7 ; {
 $gmapd/bin/gsnap -N 1 -k 14 --local-splice-penalty=1 --indel-penalty=1 --pairmax=5000 -A "sam" \
    --part=$i/16 -D $workd/genome/gmap -d $dgenome $drna.fa2 > $drna.gsnap$i.samu &

}

wait

# for i in  8 9 10 11 12 13 14 15 ;

