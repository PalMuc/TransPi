#! /bin/bash
### env fasta=xxx.fa blastp=xxx.blastp.gz onam=xxx odir=yyy datad=`pwd` qsub -q normal omcl_evgrun.sh
#PBS -N omcl_evgrun
#PBS -A ind114
#PBS -l nodes=1:ppn=32,walltime=27:55:00
#PBS -V

ncpu=32
ntaxa=11
domake=1
bopts="-identmin 33"

## omcl_evgrun.sh : full run, ~280,000 proteins, 11 species; runtime 10hr x 11cpu used == ntaxa
## versus 1cpu run of original omcl: pushing 48hr.
## blastp,fasta needed for domake step only
## for reruns after domake, input onam=xxx and odir=bbhdir to reuse data  

if [ "X" = "X$datad" ]; then echo "ERR: env datad=path/to/data"; exit -1; fi
if [ "X" = "X$blastp" ]; then echo "ERR: env blastp=myspecies.blastp.gz"; exit -1; fi
if [ "X" = "X$fasta" ]; then echo "ERR: env fasta=myspecies.prot.fasta"; exit -1; fi

export evigene=$HOME/bio/evigene
orlib=$evigene/scripts/omcl
#^ perl: export PERL5LIB=$evigene/scripts/omcl or perl -I$orlib as below..
export MCL=$HOME/bio/mcl9/bin/mcl
export ORTHOMCL=$datad/

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

echo "DONE omcl : `date`" 