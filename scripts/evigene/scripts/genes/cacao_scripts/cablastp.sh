#! /bin/bash -l
### qsub -q batch cablastp.sh
#PBS -N gsnap1
#PBS -A ind114
#PBS -l nodes=1:ppn=8,walltime=11:55:00
#PBS -o blastp1.out
#PBS -e blastp1.err
#PBS -V

workd=$HOME/scratch/chrs/cacao
ncbi_bin=/opt/bio/ncbi/bin

runset=cacao9_epi4
# runset=cacao9_epi3
# runset=cacao9_epit2
# runset=cacao9_cacao0
# runset=cacao9_arb
# ### runset=cacao9_epit1  # dropped

# _self set:
query=$runset

# plant set:
# query=plant_arab
# query=plant_vitis
# query=plant_populus # not done
# query=plant_prunus  # not done


cd $workd/genes/

# mkdir splitin
# mkdir splitout
# $aug/scripts/splitMfasta.pl --mins 1000000  $query.fa
# tidy: put split*fa in splitin/ ; put $i.blastp in splitout/


for i in  1 2 3 4 5 6 7 8 ;  {
# for i in  9 10 11 12 13 14 15 16 17 18; {

  queryiaa=$query.split.$i.fa
  queryinam=$query.$i
  if test -f splitin/$queryiaa ; then
  $ncbi_bin/blastall -p blastp -e1e-5 -m9  \
    -i splitin/$queryiaa -d $runset.aa -o splitout/$runset-$queryinam.blastp &
  fi

}

wait

