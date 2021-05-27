#! /bin/bash
### env fain=allpe.fa qsub -q batch dignorm.sh
#PBS -N dignorm
#PBS -A ind114
#PBS -l mem=36gb,nodes=1:ppn=2,walltime=28:55:00
#PBS -o dignorm.$$.out
#PBS -e dignorm.$$.err
#PBS -V

# use env param? need '-p' for pe, not for sr
# redo, 2 runs
# fainpe=allpe.fa; fainsr=allsr.fa

#run# fain=allpe.fa; kopt="-p";
fain=allsr1.fa; kopt="";

kmer=30
# 8e9 x 4 = 32Gb
xhash=8e9
rkeep=20
ncpu=2

#** requires python2.6+ not on gordon.sdsc, on trestles module add python
module() { eval `/opt/modules/Modules/3.2.5/bin/modulecmd bash $*`; }
module add python

export PYTHONPATH=$HOME/bio/khmer/python
# kmerapp=$HOME/bio/khmer/scripts/normalize-by-median.py
datad=$HOME/scratchn
workd=$datad/chrs/aabugs/tsa/plants
fastad=$workd/sraf

cd $fastad

# ~/bio/sratoolkit/fastq-dump --fasta 0 SRR346404.sra  : paired spots
#  --split-spot, have 60bp  pair reads now; use -p

# fain=$fainpe
# $HOME/bio/khmer/scripts/normalize-by-median.py -p -C $rkeep -x $xhash  -k $kmer  $fain
# mv $fain.keep $fain.dnorm$kmer.fa

# fain=$fainsr
if [ -f $fain ]; then
 $HOME/bio/khmer/scripts/normalize-by-median.py $kopt -C $rkeep -x $xhash  -k $kmer  $fain
 mv $fain.keep $fain.dnorm$kmer.fa
fi

## NO .. split to fa /1,/2 ; measure rdsize if not given?
# cat $fain.keep | perl -ne \
# 'BEGIN{ $RL=0; $MAXN=9; } if(/^>(\S+)/) { $h=$1; } else { chomp; 
# $nn= tr/N/N/; next if($nn>$MAXN);
# $RL=int(length($_)/2); $a=substr($_,0,$RL); $b=substr($_,$RL); 
# print ">$h/1\n$a\n",">$h/2\n$b\n"; }' > $fain.dnorm$kmer.fa12


