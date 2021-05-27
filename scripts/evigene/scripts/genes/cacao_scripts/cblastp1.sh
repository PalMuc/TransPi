#! /bin/bash -l
### qsub -q batch blastpcu1.sh
#PBS -N blastpcu1.n
#PBS -A ind114
#PBS -l nodes=1:ppn=8,walltime=11:55:00
#PBS -o blastpcu1.out
#PBS -e blastpcu1.err
#PBS -V

workd=$HOME/scratch/chrs/cacao
ncbi_bin=/opt/bio/ncbi/bin

## allxall blast of 300+k prots; 109 parts; 14 jobs/nodes
runset=plant9
query=$runset

cd $workd/genes/

for i in  1 2 3 4 5 6 7 8 ;  
{
  queryiaa=$query.split.$i.fa
  queryinam=$query.$i
  if test -f splitin/$queryiaa ; then
  $ncbi_bin/blastall -p blastp -e1e-5 -m9  \
    -i splitin/$queryiaa -d $runset.aa -o splitout/$runset-$queryinam.blastp &
  fi
}

wait

#EOF ------------------------
# cat plant8uniprot.fa cacao9_mix6.aa > plant9.fa
# $aug/scripts/splitMfasta.pl --mins 1000000 --out splitin $query.fa
# = 119 parts / 8 = 15 nodes
#cblastp[1..n].sh jobs:
# ~2000 seq/splitin; ~100 blastp/5 min; approx 1.5 hr to finish each/all 15 jobs?
#
# cat cblastp1.sh | perl -ne\
# 'last if(/^#EOF/); $s.=$_; END{ $na="cublastp"; @a=(1..119); $j=0; 
# while( @b=splice(@a,0,8) ) { $t=$s; $j++; $t=~s/ blastpcu1.(n|out|err)/ blastpcu$j.$1/sg;
# $it=join(" ",@b); $t=~s/(for i in) [\d ]+;/$1 $it ;/; 
# open(F,">$na$j.sh"); print F $t; close(F); } }'
#
