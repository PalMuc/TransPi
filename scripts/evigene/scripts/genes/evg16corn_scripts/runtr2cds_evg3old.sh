#! /bin/bash
### env datad=path/to/data trset=myspecies_all.tr.gz qsub -q normal runtr2cds.sh
#PBS -N tr2cds
#PBS -A ind114
#PBS -l nodes=1:ppn=16,walltime=47:55:00
#PBS -V

## FIXMEd: upd to tr2aacds2.pl **
# comet.sdsc: caller srun_prog.sh sets ncpu, 24 cores/node, 128G mem
if [ "X" = "X$ncpu" ]; then ncpu=15; fi
if [ "X" = "X$maxmem" ]; then maxmem=50000; fi
evigene=$HOME/bio/evigene/scripts
evapp=$evigene/prot/tr2aacds2.pl

## 2016 new opts for asmrna_dupfilter: aminbad/poo should be same as aamin, as tiny aa usually have "normal" size utr
## new opt only in asmrna_dupfilter3.pl update: aagapmax = BAD_GAPS = 5% better try, was 25% too high
## this aaminbad change adds lots of new, excess fragment loci, need to be filtered out by other means after homol test

export asmrna_dupfilter2=$evigene/rnaseq/asmrna_dupfilter3.pl

export aagapmax=5
export aamin=30
export aapart=120
export aaminbad=30
export aaminpoo=30

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

echo "#START $trset : `date`"
echo $evapp -NCPU $ncpu -MAXMEM $maxmem $traopts -cdna $trset
$evapp -NCPU $ncpu -MAXMEM $maxmem $traopts -cdna $trset
echo "#DONE $trset : `date`"

