#! /bin/bash
### env datad=path/to/data trset=myspecies_all.tr.gz qsub -q normal runtr2cds.sh
#PBS -N tr2cds
#PBS -A ind114
#PBS -l nodes=1:ppn=16,walltime=47:55:00
#PBS -V

if [ "X" = "X$ncpu" ]; then ncpu=2; fi
if [ "X" = "X$maxmem" ]; then maxmem=5000; fi
if [ "X" = "X$evigene" ]; then
   if [ -d $HOME/bio/evigene ]; then evigene=$HOME/bio/evigene; fi
   if [ -d /bio/mb/evigene ]; then evigene=/bio/mb/evigene; fi
   if [ "X" = "X$evigene" ] & [ -d $evigene ]; then 
     evgok=1; 
   else 
    echo missing evigene=/path/to/evigene; exit -1; 
   fi
fi
if [ "X" = "X$trset" ]; then echo "missing env trset=xxxx.tr"; exit -1; fi
if [ "X" = "X$datad" ]; then echo "missing env datad=/path/to/data "; exit -1; fi

evapp=$evigene/scripts/prot/tr2aacds2.pl
export asmrna_dupfilter2=$evigene/rnaseq/asmrna_dupfilter3.pl

pn=fastanrdb;  pp=`which $pn`; pp=`basename $pp`; if [ "$pp" != "$pn" ]; then echo missing path to $pn; exit -1; fi
pn=cd-hit-est; pp=`which $pn`; pp=`basename $pp`; if [ "$pp" != "$pn" ]; then echo missing path to $pn; exit -1; fi
pn=blastn;     pp=`which $pn`; pp=`basename $pp`; if [ "$pp" != "$pn" ]; then echo missing path to $pn; exit -1; fi

## 2016 new opts for asmrna_dupfilter: aminbad/poo should be same as aamin, as tiny aa usually have "normal" size utr
export aagapmax=10
export aamin=30
export aapart=120
export aaminbad=30
export aaminpoo=30

traopts="-tidy -log -debug"
if [ "X" != "X$aablast" ]; then traopts="$traopts -ablastab $aablast"; fi

# lastz, etc in bio/bin: export PATH=$HOME/bio/bin:$PATH
# cd-hit: export PATH=$HOME/bio/cdhit/bin:$PATH
# fastanrdb: export fastanrdb=$HOME/bio/exonerate/bin/fastanrdb
# blastn: export PATH=$HOME/bio/ncbi2230/bin:$PATH

cd $datad/

echo "#START $trset : `date`"
echo $evapp -NCPU $ncpu -MAXMEM $maxmem $traopts -cdna $trset
$evapp -NCPU $ncpu -MAXMEM $maxmem $traopts -cdna $trset
echo "#DONE $trset : `date`"

