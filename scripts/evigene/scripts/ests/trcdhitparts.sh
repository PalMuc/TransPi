#! /bin/bash
### env  trset=myspecies_all.tr datad=path/to/data prog=./trcdhitparts.sh sbatch srun_comet.sh
#PBS -N trcdhitparts
#PBS -A ind114
#PBS -l nodes=1:ppn=24,walltime=39:55:00
#PBS -V

## test effect of cd-hit-est trasm filter on reducing cds quality (select for broken cds)
## cd-hit-est -in trasm.fasta > trasm_cd90.tr | bestorf > trasm_cd90.cds > blastn -q ref.cds -db trcd90.cds + trfull.cds
## see also /bio-grid/plants/arabidopsis/geneval/tr2cdsreduce.sh

## ** cdhit -c 0.95 runs fast (<2hr, same data), -c 0.90 outatime with 40hr/24cpu/124gb (weed genes)
## ARGHH.. way too long to run: 24 cpu * 120 GB mem * 40hr wont finish basic velvet or trin trasm set
## .. smaller trsets finish quicker; do data split: 6 parts x 4cpu/part x 20 Gb/part, then cat part_cdtr and finish

RUN=1

cdhitbin=$HOME/bio/cdhit466/bin
# old: $HOME/bio/cdhit/bin
# option: cdhit -c 0.90 default, test and other levels? 1.00, 0.99, 0.95, ..

cdhitcut_default=90
maxmem=120000

if [ "X" = "X$ncpu" ]; then ncpu=24; fi
if [ "X" = "X$mem" ]; then mem=$maxmem; fi
if [ "X" = "X$trset" ]; then echo "missing env trset=xxxx.tr"; exit -1; fi
if [ "X" = "X$datad" ]; then echo "missing env datad=/path/to/data "; exit -1; fi
if [ "X" = "X$cdcut" ]; then cdcut=$cdhitcut_default; fi
CDCUT="0.$cdcut";  echo "opt cdhitcut=$CDCUT";

ncpurun=4
# ncpupart=6
ncpupart=$(( $ncpu / $ncpurun ))
memrun=20000
# memrun=$(( $mem / $ncpupart ))

evigene=$HOME/bio/evigene;
biohome=$HOME/bio;
export PATH=$cdhitbin:$biohome/bin:$biohome/exonerate/bin:$biohome/ncbi/bin:$PATH

#x qname=`echo $trset | sed 's/\.gz//; s/\.[a-z]*$//;'`
qname=`basename $trset .tr | sed 's/\.[a-z]*$//;'`
odir=trsplit$qname

cd $datad/
mkdir $odir

echo "#START $trset : `date`"

if [ ! -f $trset.split.1.fa ]; then
 pindir=`dirname $trset`
 splitsize=`grep -v '^>' $trset | wc -c | sed 's/ .*//' `
 splitbp=$(( $splitsize / $ncpupart ))
 $evigene/scripts/splitMfasta.pl --outputpath=$pindir --nparts $ncpupart --minsize=$splitbp $trset
fi

i=0;
qset=`/bin/ls $trset.split.*.fa`
for qfile in $qset
{
  qnamspl=`basename $qfile .fa`
  cdcutname=$odir/$qnamspl-cd$cdcut.cdna
  echo cd-hit-est -c $CDCUT -T $ncpurun -M $memrun -l 89 -d 0 -i $qfile  -o $cdcutname
  if [ $RUN = 1 ]; then
    cd-hit-est -c $CDCUT -T $ncpurun -M $memrun -l 89 -d 0 -i $qfile -o $cdcutname &
    #x if [ ! -s $cdcutname ]; then echo "fail cd-hit-est $cdcutname"; exit -1; fi
    # grep '^>' $cdcutname | sed 's/>//; s/ .*//;' | sort -u > $cdcutname.ids
  fi
  i=$(( $i + 1 ));  if [ $i -ge $ncpupart ]; then wait; i=0; fi
}

wait

## catenate all part cd$cdcut.cdna and rerun
trinfile=$qname-cd$cdcut-parts.tr
cdcutname=$qname-cd$cdcut-all.cdna
cat $odir/$qname*-cd$cdcut.cdna > $trinfile

# all ncpu, mem here
echo cd-hit-est -c $CDCUT -T $ncpu -M $mem -l 89 -d 0 -i $trinfile -o $cdcutname
if [ $RUN = 1 ]; then
  cd-hit-est -c $CDCUT -T $ncpu -M $mem -l 89 -d 0 -i $trinfile -o $cdcutname
  grep '^>' $cdcutname | sed 's/>//; s/ .*//;' | sort -u > $cdcutname.ids
  echo /bin/rm -r $odir
fi

echo "#DONE $trset : `date`"
#------------------------------

