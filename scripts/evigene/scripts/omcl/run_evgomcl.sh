#! /bin/bash
### env ntaxa=8 fasta=prots.fa blastp=prots.blastp.gz onam=xxx odir=yyy datad=`pwd` qsub -q shared omclrun.sh
#PBS -N omclevgrun
#PBS -A ind114
#PBS -l nodes=1:ppn=24,walltime=37:55:00
#PBS -V

# ntaxa=11
# omclevgrun13d.sh : full 1st run

if [ "X" = "X$ncpu" ]; then ncpu=16; fi 
if [ "X" = "X$ntaxa" ]; then echo "ERR: env ntaxa=how many?";  exit -1; fi
if [ "X" = "X$datad" ]; then echo "ERR: env datad=path/to/data"; exit -1; fi
if [ "X" = "X$blastp" ]; then echo "ERR: env blastp=myspecies.blastp.gz"; exit -1; fi
if [ "X" = "X$fasta" ]; then echo "ERR: env fasta=myspecies.prot.fasta"; exit -1; fi

evigene=$HOME/bio/evigene
orlib=$evigene/scripts/omcl
mcl=$HOME/bio/mcl9/bin/mcl

domake=1
bopts="-identmin 33"

## make onam, odir from fasta or blastp name
if [ "X" = "X$onam" ]; then 
  # echo "ERR: need env onam=xxx.bpo "; exit -1;
  onam=`basename $fasta .fasta | sed 's/\.aa.*//; s/\.pep.*//; s/\.fa.*//; s/$/_omcl/;'`
fi
if [ "X" = "X$odir" ]; then
  # echo "ERR: need env odir=path with all spp.bbh "; exit -1;
  odir="${onam}$$"
fi

echo "START omcl : `date` "
cd $datad

if [ $domake = 1 ]; then
 $evigene/scripts/blast92orthomcl10.pl $bopts  -fasta=$fasta -in=$blastp -out=$onam
fi

export MCL=$mcl 
export ORTHOMCL=$datad/

perl -I$orlib $orlib/orthomcl_evg.pl --mode par_start --former_run_dir $odir --bpo=$onam.bpo --gg=$onam.gg 

t=0; i=0; 
while [ $t -lt $ntaxa ]; do { 
 perl -I$orlib $orlib/orthomcl_evg.pl --mode par_part$t --former_run_dir $odir --bpo=$onam.bpo --gg=$onam.gg & 
 
 t=$(( $t + 1 )); i=$(( $i + 1 )); 
 if [ $i -ge $ncpu ]; then wait; i=0; fi;
} done 
wait

## now do connect$t loop
t=0; i=0;
while [ $t -lt $ntaxa ]; do {
 perl -I$orlib $orlib/orthomcl_evg.pl --mode par_connect$t --former_run_dir $odir --bpo=$onam.bpo --gg=$onam.gg &

  t=$(( $t + 1 )); i=$(( $i + 1 ));
  if [ $i -ge $ncpu ]; then wait; i=0; fi;
} done
wait

perl -I$orlib $orlib/orthomcl_evg.pl --mode par_end --former_run_dir $odir --bpo=$onam.bpo --gg=$onam.gg 

##.. add this tabulation, if run ok; what configs?
# env sppmap=FISH11 myspecies=Fundulus orthoutisokay=1 xml=0 date=20131225 clade=Fish \
#  goodname='mayzebr|human|kfish2|platyfish' poorname='stickleback|medaka|tetraodon' \
#  sppindex="0,2,3,4,5,7,8,9,10" tmin=7 $evigene/scripts/omcl/orthomcl_tabulate.pl \
#   -debug -idprefix="FISH${pt}_G" -gtag="FISH$pt" -omclpath ./ -namepath ../names/

echo "DONE omcl : `date`" 


