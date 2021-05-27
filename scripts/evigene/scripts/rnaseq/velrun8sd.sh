#!/bin/bash
#==> vel8/velrun8.sh <==
##  env subset=sc8 qsub -q normal velrun8.sh
#PBS -A ind114
#PBS -N velrun 
#PBS -l nodes=1:ppn=32,walltime=21:55:00
#.. lower time=22:  all vel6 done 3kmers in < 12hr
#PBS -o velrun.$$.out
#PBS -e velrun.$$.err
#PBS -V

# no longPaired in cacao data; shortPaired + long
## vel8: switch to sub2.* data from sub3. .. has stricter error filter, likely aid for vel1-4 results
#...
## note: run6 is poorer result than prior vel1-5, maybe due to 
#  a. higher kmer? (prior used 21,27, here k23 is longest), b. opts  -conserveLong yes/no ? 
#  .. vel1 has most prot omcl best matches by 2x ; vel3 next best
#
#v0# kmer=21 #v1# kmer=27 #v2# kmer=29 #v3# kmer=27 #v4# kmer=31 #v5# kmer=31
## -min_trans_lgth 40 on all, 100 should be ok
## -min_pair_count 4  default used on all
#
#.. Data note: vel1-4 used sub1,sub2*.fa w/ diff data selection; 
#.. vel5,6,7 use sub3 with maxdup=8?,other err cuts
#..
#v1# -conserveLong no -ins_length 200 
#v2# -conserveLong no -ins_length 250 -ins_length_sd 180
## ins av/sd cgb:204/57; ncgr 090511=115/71; 090922_8=292/240; 091005_2=306/196; 090609_3=382/199
#v3: -kmer 27; -ins_len 200 ; -ins_sd 90; -conserveLong yes ; -minpair 3
#v3: shortPaired = ins_length 200, ins_length_sd 60; shortPaired2 = ins_length2 340, ins_length2_sd 200
#v3# -conserveLong yes -ins_length 200 -ins_length_sd 60 -ins_length2 350 -ins_length2_sd 150
#v4# -conserveLong no  -ins_length 200 -ins_length_sd 40 -ins_length2 350 -ins_length2_sd 90
#v5# -ins_length 200 -ins_length2 350 -ins_length2_sd 60
#..........
#v6# best orthogroup kmers: 3162 k23, 2419 k27, 4881 k39  <<? k39 best for many; maybe extreme kmers?
#-------

# runversion
dv=8
# vel5,6,7: datasub=sub3.$subset
datasub=sub2.$subset

ncpu=8
export OMP_NUM_THREADS=$ncpu
export OMP_THREAD_LIMIT=$ncpu

rund=/scratch/$USER/$PBS_JOBID
datad=/phase1/$USER
velbin=$datad/bio/velvet/bin
workd=$datad/chrs/cacao/rnas/vel8
homed=$HOME/work/

#v6
#vopts="-ins_length 450 -ins_length_sd 80 -ins_length2 200 -ins_length2_sd 40"
#oopts="-min_trans_lgth 100 -ins_length 450 -ins_length_sd 80 -ins_length2 200 -ins_length2_sd 40"
#v7,v8
vopts="-ins_length 200 -ins_length2 400"
oopts="-min_trans_lgth 100 -ins_length 200 -ins_length2 400"

#v6#kset="39 27 23"
#v7#kset="33 21"
kset="35 23"
ksubdir=vel$dv${subset}

notef=$workd/$ksubdir.$$.RUNNING
donef=`echo $notef | sed 's/RUNNING/DONE/'`

touch $notef
echo "START " >> $notef
echo `date`   >> $notef

mkdir -p $rund
cd $rund/
cp $workd/velfa/$datasub.* $rund/

#..... run loop for vel kmer steps
shopt -s nullglob

for k in $kset;  do { 

ksubdir=vel$dv${subset}_$k
echo "#.. start velrun $ksubdir : `date`"

#. du -h >> $notef
#. echo "velveth $ksubdir" >> $notef
$velbin/velveth $ksubdir $k -fasta.gz \
    -shortPaired $datasub.*.sifa2.gz $datasub.*.fa2.gz  \
    -shortPaired2 $datasub.*.lifa2.gz \
    -short $datasub.*.fa1.gz   -long $datasub.longfa.gz 

#. echo "velvetg : `date`" >> $notef
#. ls -l $ksubdir >> $notef
$velbin/velvetg $ksubdir $vopts -read_trkg yes ; 

#. echo "oases : `date`" >> $notef
#. ls -l $ksubdir >> $notef
$velbin/oases $ksubdir $oopts 

cp -p $ksubdir/transcripts.fa $homed/$ksubdir-transcripts.fa
/bin/rm $ksubdir/{Graph2,LastGraph,PreGraph,Roadmaps,Sequences}
cp -rp $ksubdir $workd/

#. ls -l $ksubdir >> $notef
/bin/rm -rf $ksubdir

echo "#.. end velrun $ksubdir : `date`"

}
done

# not wait

echo "DONE" >> $notef
echo `date`   >> $notef
mv $notef $donef


