#!/bin/bash
##  env subset=sc8 kmer=39 qsub -q normal velrun6.sh
#PBS -A ind114
#PBS -N velrun 
#PBS -l nodes=1:ppn=32,walltime=31:55:00
#PBS -o velrun.$$.out
#PBS -e velrun.$$.err
#PBS -V

# using velvet velvet_1.2.03/oases_0.2.06 2012-april
# .. rnaseq-pe is partitioned to long inserts (400..500bp) and shortin (100..300bp), + long ESTs (unpaired)
# .. this is similar to successful locust rna-pe short/long ins
# trestles.sdsc, limit of 64GB mem per node is cramped for velvet 
# .. try normal not shared for mem/ssd : 32 cores but all 64GB mem + SSD disk (120GB? should be enough)
# .. run 1 kmer at a time, limit cores to 8 so memory isnt overloaded (more cpu = more memory)
# .. small data set, 3 kmers, cacao3sc8:  5hr;  k39=1:20 ; k27=1:45 ; k23=2hr ; memuse=30GB
# ..   cacaosc10r: 3:30 hr, 22089681/22669501 reads for k23; 

# runversion
dv=6c
ncpu=8
export OMP_NUM_THREADS=$ncpu
export OMP_THREAD_LIMIT=$ncpu

rund=/scratch/$USER/$PBS_JOBID
datad=/phase1/$USER
velbin=$datad/bio/velvet/bin
workd=$datad/chrs/cacao/rnas/vel6
homed=$HOME/work/

kset="39 27 23"
#k=$kmer
k=39
ksubdir=vel$dv${subset}_$k

notef=$workd/$ksubdir.$$.RUNNING
donef=`echo $notef | sed 's/RUNNING/DONE/'`

touch $notef
echo "START " >> $notef
echo `date`   >> $notef

mkdir -p $rund
cd $rund/
cp $workd/velfa/sub3.$subset.* $rund/

#..... run loop for vel kmer steps
shopt -s nullglob

vopts="-ins_length 450 -ins_length_sd 80 -ins_length2 200 -ins_length2_sd 40"
oopts="-min_trans_lgth 100 -ins_length 450 -ins_length_sd 80 -ins_length2 200 -ins_length2_sd 40"

for k in $kset;  do { 

ksubdir=vel$dv${subset}_$k
echo "#.. start velrun $ksubdir : `date`"
du -h >> $notef

echo "velveth $ksubdir" >> $notef
$velbin/velveth $ksubdir $k -fasta.gz  -shortPaired sub3.$subset.*.lifa2.gz \
    -shortPaired2 sub3.$subset.*.sifa2.gz sub3.$subset.*.fa2.gz  \
    -short sub3.$subset.*.fa1.gz -long sub3.$subset.longfa.gz 

echo "velvetg : `date`" >> $notef
ls -l $ksubdir >> $notef
$velbin/velvetg $ksubdir $vopts -read_trkg yes ; 

echo "oases : `date`" >> $notef
ls -l $ksubdir >> $notef
$velbin/oases $ksubdir $oopts 

# 1st save trans.fa ; flaky scratch disks..
cp -p $ksubdir/transcripts.fa $homed/$ksubdir-transcripts.fa
/bin/rm $ksubdir/{Graph2,LastGraph,PreGraph,Roadmaps,Sequences}
cp -rp $ksubdir $workd/

ls -l $ksubdir >> $notef
/bin/rm -rf $ksubdir

echo "#.. end velrun $ksubdir : `date`"

}
done

# not wait

echo "START `date` " >> $notef
mv $notef $donef

#... flaky scratch disk.. maybe havent lost any runs?
## for trs in $caca/rnas/vel6/vel6csc[x256]_*/transcripts.fa ; do 
# for trs in $caca/rnas/vel6/vel6csc[135]_*/transcripts.fa ; do 
# {
#   nam=`echo $trs | sed 's,/phase1/ux455375/chrs/cacao/rnas/vel6/,,; s,/,-,g;'`
#   cp -p $trs $HOME/work/vel6trs/$nam  
# } done
