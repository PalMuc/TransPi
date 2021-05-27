#! /bin/bash
### env datad=path/to/data trset=myspecies_all.tr.gz qsub -q normal runtr2cds.sh
#PBS -N fareduce
#PBS -A ind114
#PBS -l nodes=1:ppn=16,walltime=47:55:00
#PBS -V

maxmem=120000
# comet.sdsc: caller srun_prog.sh sets ncpu, 24 cores/node, 128G mem
if [ "X" = "X$ncpu" ]; then ncpu=15; fi
if [ "X" = "X$maxmem" ]; then maxmem=50000; fi
evigene=$HOME/bio/evigene/scripts
evapp=$evigene/prot/tr2aacds2.pl

## 2014 new opts for asmrna_dupfilter2:
## tcas4/tribol beetle: 2014 ncbi has 10 prot under 40aa, 5 are NP curated/refseq; 
export aamin=30
export aapart=120
export aaminbad=90
export aaminpoo=60

traopts="-tidy -log -debug"

# lastz, etc in bio/bin
export PATH=$HOME/bio/bin:$PATH
# cd-hit
export PATH=$HOME/bio/cdhit/bin:$PATH
# fastanrdb
export fastanrdb=$HOME/bio/exonerate/bin/fastanrdb
# blastn:
## upd: ncbi2230, was ncbi2227 ?
export PATH=$HOME/bio/ncbi2230/bin:$PATH

if [ "X" = "X$trset" ]; then echo "missing env trset=xxxx.tr"; exit -1; fi
if [ "X" = "X$datad" ]; then echo "missing env datad=/path/to/data "; exit -1; fi
cd $datad/

trname=`echo $trset | sed 's/\.[a-z]*$//;'`
echo "#START $trset : `date`"

#t2ac: CMD= /home/ux455375/bio/exonerate/bin/fastanrdb evg1corn.cds > evg1cornnr.cds
#t2ac: CMD= /home/ux455375/bio/cdhit/bin/cd-hit-est  -c 1.00 -T 20 -M 120000 -l 89 -d 0 -i evg1cornnr.cds -o evg1cornnrcd1.cds 1> evg1cornnrcd1.log 2>&1

$fastanrdb $trset > $trname.nr.fa

cd-hit-est -c 1.00 -T $ncpu -M $maxmem -l 89 -d 0 -i $trname.nr.fa -o $trname.cd1.fa  

# echo $evapp -NCPU $ncpu -MAXMEM $maxmem $traopts -cdna $trset
# $evapp -NCPU $ncpu -MAXMEM $maxmem $traopts -cdna $trset

echo "#DONE $trset : `date`"

