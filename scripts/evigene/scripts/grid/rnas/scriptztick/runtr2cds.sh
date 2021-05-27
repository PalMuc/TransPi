#! /bin/bash
### env datad=path/to/data trset=myspecies_all.tr.gz qsub -q normal runtr2cds.sh
#PBS -N tr2cds
#PBS -A ind114
#PBS -l nodes=1:ppn=16,walltime=8:55:00
#PBS -o tr2cds.$$.out
#PBS -e tr2cds.$$.err
#PBS -V

ncpu=16
maxmem=50000
evigene=$HOME/bio/evigene/scripts

# cd-hit
export PATH=$HOME/bio/cdhit/bin:$PATH
# fastanrdb
export fastanrdb=$HOME/bio/exonerate/bin/fastanrdb
# blastn:
export PATH=$HOME/bio/ncbi2227/bin:$PATH

if [ "X" = "X$trset" ]; then echo "missing env trset=xxxx.tr"; exit -1; fi
if [ "X" = "X$datad" ]; then echo "missing env datad=/path/to/data "; exit -1; fi

cd $datad/

echo $evigene/prot/tr2aacds.pl -debug -NCPU $ncpu -MAXMEM $maxmem -tidy -log -cdna $trset
$evigene/prot/tr2aacds.pl -debug -NCPU $ncpu -MAXMEM $maxmem -tidy -log -cdna $trset


