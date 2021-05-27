#!/bin/bash
##  env subset=sc6 qsub -q normal vel4merge.sh
#PBS -A ind114
#PBS -N velmerge
#PBS -l nodes=1:ppn=16,walltime=41:55:00
#PBS -o velmerge.$$.out
#PBS -e velmerge.$$.err
#PBS -V

#==> vel4/vel4merge.sh <==
# FIXME: bin3/velvet = LONGSEQ for merge of tr > 32k
#... NOTE: merge can take long time, oases not threaded, velvetg uses cpus
# ** Wrong ^^ this was really quick. 
# [0.000000] Reading FastA file alltranscripts.fa;
# [1.649595] 82171 sequences found
# [38.313805] Counted 3595 mRNA loci # BUT makes 63201 mRNA total..
# [40.149259] Found 713 missing transcripts
# [50.524146] Finished extracting transcripts, used 63224/82171 reads

## REDO: preselect besttranscripts.fa to merge using CDS/exon ratios, 
#        so merge doesnt return mostly biggest but frameshifted error transcripts

dv=4
ncpu=14
export OMP_NUM_THREADS=$ncpu
export OMP_THREAD_LIMIT=$ncpu

# gordon.sdsc
datad=/oasis/scratch/$USER/temp_project
rund=/oasis/scratch/$USER/$PBS_JOBID

velbin=$HOME/bio/velvet/bin3
workd=$datad/chrs/cacao/rnas/vel4
homed=$HOME/work/

kmerge=27
kmergedir=velm$dv${subset}
trset=$workd/vel*${subset}_*/transcripts.fa

notef=$workd/$kmergedir.$$.RUNNING
donef=`echo $notef | sed 's/RUNNING/DONE/'`
touch $notef
echo "START `date` " >> $notef
echo "mergeset $kmergedir: $trset " >> $notef

mkdir -p $rund
cd $rund/

# need subdirs or renames .. use input param?
cat $workd/vel*${subset}_*/transcripts.fa > alltranscripts.fa

# for trs in $workd/vel*${subset}_*/transcripts.fa ; do {
#   nam=`echo $trs | sed "s,$workd/,,; s,/,-,g;"` ; cp -p $trs  $nam  
# } done

# shopt -s nullglob

$velbin/velveth $kmergedir $kmerge -long alltranscripts.fa 
$velbin/velvetg $kmergedir -read_trkg yes -conserveLong yes 
$velbin/oases $kmergedir -merge yes

cp -p $kmergedir/transcripts.fa $homed/$kmergedir-transcripts.fa
/bin/rm $kmergedir/{Graph2,LastGraph,PreGraph,Roadmaps,Sequences}
cp -rp $kmergedir $workd/

/bin/rm -rf $kmergedir

echo "DONE `date` " >> $notef
mv $notef $donef

#.............................................................
#... multi-kmer runs
# kset="41 33 25"
# oopts="-min_trans_lgth 100 "
# for k in $kset; do {
#   ksubdir=vel$dv${subset}_$k
#   ( $velbin/velveth $ksubdir $k -fastq.gz -shortPaired fastq/*.fastq.gz ; \
#     $velbin/velvetg $ksubdir -read_trkg yes ; $velbin/oases $ksubdir $oopts ) &
# } done
# wait

#........... test set for merge, pick all kmer variants, data variants (v9=bams, v4,8=sub2.fa)
# vel4/cacao3vel4[34]sc6.tr vel6/cacao3vel6sc6.tr vel8/cacao3vel8sc6.tr vel9/cacao3vel92sc6.tr
# == 120k trs
