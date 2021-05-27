#! /bin/bash -l
### env rnain=aphidrs_SRR071347 qsub -q batch rnagsnap2.sh
#PBS -N gsnap2
#PBS -A ind114
#PBS -l nodes=1:ppn=8,walltime=11:55:00
#PBS -o gsnap2.$$.out
#PBS -e gsnap2.$$.err
#PBS -V

workd=$HOME/scratch/chrs/aphid2
gmapd=$HOME/bio/gmap
dgenome=aphid2asm

# input fn update for # fastq/SRA026584.1:  Dec 22 08:44 24hr_hd_cr_1.txt.gz
suf=txt
# suf=fa2
drna=`basename $rnain .$suf`

# maybe these gsnap opts=-m 3?? -N 1 -k 14
# drop these:  -local-splice-penalty=1 --indel-penalty=1 
snapopts="-N 1 -k 14 --max-mismatches=0.07 --pairmax=15000"
#  ^^ max-mis 0.0n slighly > default;  works for all read len: 37bp = 2; 50bp = 3; 72bp = 5

cd $workd/rnas/

if ! test -f "$drna.$suf" ; then  echo "missing $drna.$suf"; exit;  fi

for i in  8 9 10 11 12 13 14 15 ; {
 $gmapd/bin/gsnap $snapopts -A "sam" --part=$i/16 \
   -D $workd/genome/gmap -d $dgenome $drna.$suf > $drna.gsnap$i.samu &

}

wait

