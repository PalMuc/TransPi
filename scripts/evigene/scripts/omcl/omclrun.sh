#! /bin/bash
### env fasta=myspp.aa blastp=myspp.blastp.gz datad=`pwd` qsub -q shared omclrun.sh
#PBS -N omclrun
#PBS -A ind114
#PBS -l nodes=1:ppn=2,walltime=47:55:00
#PBS -V

ncpu=2

# see also evigene/scripts/omcl/omclrun_spp.sh
#  to calc omcl with loop over species subsets

# datad=$HOME/scratchn/chrs/kfish/prot
# blin=outz/sd-vert8kfish9fnoalt.blastp.gz
# fasta=aaset2/vert8kfish9fnoalt.aa
if [ "X" = "X$datad" ]; then echo "ERR: env datad=path/to/data"; exit -1; fi
if [ "X" = "X$blastp" ]; then echo "ERR: env blastp=myspecies.blastp.gz"; exit -1; fi
if [ "X" = "X$fasta" ]; then echo "ERR: env fasta=myspecies.prot.fasta"; exit -1; fi

evigene=$HOME/bio/evigene
orlib=$HOME/bio/orthomcl3
mcl=$HOME/bio/mcl9/bin/mcl

domake=1
#orig#bopts=""
bopts="-identmin 33"
# -identmin will filter out trivial matches that orthomcl drops anyway

## make onam, odir from fasta or blastp name
onam=`basename $fasta .fasta | sed 's/\.aa.*//; s/\.pep.*//; s/\.fa.*//; s/$/_omcl/;'`
odir="${onam}_f$$"
#onam=vert8kfna_omcl
#odir=omclkf2

echo "START omcl `date` "
cd $datad
mkdir $odir

if [ $domake ]; then
 $evigene/scripts/blast92orthomcl10.pl $bopts  -fasta=$fasta -in=$blastp -out=$odir/$onam
fi

cd $datad/$odir/
export MCL=$mcl 
export ORTHOMCL=$datad/$odir/
perl -I$orlib $orlib/orthomcl.pl --mode 4 --bpo=$onam.bpo --gg=$onam.gg 

echo "DONE omcl : `date`" 

