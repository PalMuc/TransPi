#! /bin/bash -l
### qsub -q batch cablastp.sh
#PBS -N cablastp
#PBS -A ind114
#PBS -l nodes=1:ppn=8,walltime=11:55:00
#PBS -o blastp1.out
#PBS -e blastp1.err
#PBS -V

workd=$HOME/scratch/chrs/cacao
ncbi_bin=/opt/bio/ncbi/bin

runs=(cacao9_cacao0 cacao9_epi3 cacao9_epi5 cacao9_epit2 cacao9_cacao4 cacao9_epi4 cacao9_epir6 cacao9_mix6)
ri=7 # iterate this 0..7
runset=${runs[$ri]}

##drop.gene2# runset=cacao9_arb # done1 done2 dn3
## redo this as loop : NO? prevents multi que jobs
### runset=cacao9_epit1  # skip
### runset=cacao9_fgenesh_fix

# _self set:
# query=$runset

query=uniprot_sprot_plants
## plant set
# query=plant_arab # done
# query=plant_vitis # done
# query=plant_populus #done
# query=plant_prunus

cd $workd/genes/
# mkdir splitin; # mkdir splitout # $aug/scripts/splitMfasta.pl --mins 1000000 --out splitin $query.fa

if ! test -f $runset.aa.psq ; then
  ## if ! test -f $runset.aa ; then exit 1 fi
  $ncbi_bin/formatdb -pT -i $runset.aa
fi

for i in  1 2 3 4 5 6 7 8 ;  {
  queryiaa=$query.split.$i.fa
  queryinam=$query.$i
  if test -f splitin/$queryiaa ; then
  $ncbi_bin/blastall -p blastp -e1e-5 -m9  \
    -i splitin/$queryiaa -d $runset.aa -o splitout/$runset-$queryinam.blastp &
  fi
}

wait

