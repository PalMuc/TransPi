#! /bin/bash -l
### qsub -q batch cacestgmap1.sh
#PBS -N estgmap1b
#PBS -A ind114
#PBS -l nodes=1:ppn=8,walltime=11:55:00
#PBS -o estgmap1b.out
#PBS -e estgmap1b.err
#PBS -V

workd=$HOME/scratch/chrs/cacao
gmapd=$HOME/bio/gmap

dgenome=cacao9asm
# cacao1asm
# dest=cacao1est # done
# dest=cacao2est # done
# dest=cacao1leafreads_s1r1clp # done
# dest=cacao1leafreads_s1r2clp # done
# dest=cacao1leafreads_s2r1clp
# dest=cacao1floralreads_s3r1clp  # done
# dest=cacao1floralreads_s3r2clp  # done
dest=dbest091211_wsuclean

cd $workd/est/

#pt2: for i in  8 9 10 11 12 13 14 15 ;
for i in  0 1 2 3 4 5 6 7 ; {

  $gmapd/bin/gmap -D $workd/genome/gmap -d $dgenome --jobdiv=$i/16 -n 4 -S  $dest.fa > $dest.gmap$i.out &

}

wait


