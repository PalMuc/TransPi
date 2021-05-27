#!/bin/bash
#==> vel4/velrun4g.sh <==
##  env subset=sc8 qsub -q normal velrun9.sh
#PBS -A ind114
#PBS -N velrun 
#PBS -l nodes=1:ppn=16,walltime=21:55:00
#.. lower time=22:  all vel6 done 3kmers in < 12hr
#PBS -o velrun.$$.out
#PBS -e velrun.$$.err
#PBS -V

#.. vel40: rerun with old velvet10 (1.0..) similar to vel4 run last fall .. is it soft vers bug?
#.. note no OPENMP in vel10; kmer=31 max ; use that only
#.. vel43: velvet10 is not answer, revert newest vel/o, use only .fa2 paired data
# vel43 ESTbest kmers: 12251 sc6k25    32847 sc6k35
#
# *** add Merge step to condense all kmer sets; see prior, or oases_0.2.06/scripts/oases_pipeline.py
# kmerge=27
# vmergedir=velm$dv${subset}
# $velbin/velveth $vmergedir $kmerge -long vel$dv${subset}*/transcripts.fa 
# $velbin/velvetg $vmergedir -read_trkg yes -conserveLong yes 
# $velbin/oases $vmergedir -merge yes
#
#.. vel9: problem appears to be lowqual data, test w/ higher qual subset, .bam format
#.. v9a: k35 ok, k23 failed outofment (64G)  at velvetg; using more mem than fasta variant. more reads?
#.. v9b: remove .ncgr bams, keep .cgb only: tr errors may be in part due to strain snps ..
#
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
# vel5,6,7: datasub=sub3.$subset
dv=47
datasub=sub2.$subset

# reads cleaned/trimmed of lowqual (sickle)
dv=4t1
datasub=trim.$subset

ncpu=8
export OMP_NUM_THREADS=$ncpu
export OMP_THREAD_LIMIT=$ncpu

#tres#datad=/phase1/$USER
datad=/oasis/scratch/$USER/temp_project

#tres#rund=/scratch/$USER/$PBS_JOBID
rund=/oasis/scratch/$USER/$PBS_JOBID

velbin=$HOME/bio/velvet/bin
#OLD#velbin=$HOME/bio/velvet10/bin
workd=$datad/chrs/cacao/rnas/vel4
homed=$HOME/work/

## ** FIXME FIXME : velvet wants read len added to inslen, for cgb 106b pairs, full inslen = 412
## ** FIXME222222 : .bams need to be ID-sorted, so mates follow each other, like .fa2: samtools sort -n
#v6
#vopts="-ins_length 450 -ins_length_sd 80 -ins_length2 200 -ins_length2_sd 40"
#oopts="-min_trans_lgth 100 -ins_length 450 -ins_length_sd 80 -ins_length2 200 -ins_length2_sd 40"
#v7,v8,v9
#vopts="-ins_length 200 -ins_length2 400"
#oopts="-min_trans_lgth 100 -ins_length 200 -ins_length2 400"
#
#v9c : cgb ins is 400 (212b-rd + 200 ins); nc short 108-rd+ 200ins = 300; nc long=108r+400i=550
#v9 # vopts="-min_pair_count 4 -ins_length 410 -ins_length2 550"
#v4 : drop velvetg opts .. maybe helps
# vopts=""
#v43# oopts="-min_pair_count 4 -min_trans_lgth 100 -ins_length 410 -ins_length2 550"
# oopts="-min_pair_count 4 -min_trans_lgth 100 -ins_length 350 -ins_length2 550"

#v45.46  bin2/velvet has cat=3 for this
vopts="-ins_length 300 -ins_length2 410 -ins_length3 550"
oopts="-min_pair_count 4 -min_trans_lgth 100 -ins_length 300 -ins_length2 410 -ins_length3 550"

#v6#kset="39 27 23"
#v7#kset="33 21"
#v8; k35,23 > v9 k35,25
#v9# kset="35 25"
#v42# kset="31 25"
#v43# kset="35 25"
#.. k55 is >readsize of ncgr (54) .. should have used 53 max
#v44 kset="55 45 39 29"
#v45 kset="101 91 81 71 61 51 41 31"
#.. v45 k91 is single best; v45 combined nearly matches ESTqual of trinity
#v46 kset="89 49 29 25 23 21"

#v47 .. all data, from k-bestest set ** most long aa from k49, but k23,29,35 are similarly high
kset="93 89 79 69 55 49 35 29 25 23"
#v4t1: 
kset="69 51 47 41 35 29 25"

ksubdir=vel$dv${subset}

notef=$workd/$ksubdir.$$.RUNNING
donef=`echo $notef | sed 's/RUNNING/DONE/'`

touch $notef
echo "START `date` " >> $notef

mkdir -p $rund
cd $rund/
cp $workd/velfa/$datasub.*gz $rund/
# cp $workd/bamsc/$datasub.*bam $rund/

#..... run loop for vel kmer steps
shopt -s nullglob

#v44: add back longfa
#v43NOT.#    -short $datasub.*.fa1.gz   -long $datasub.longfa.gz 
## revise above to use velveth -noHash >> Seqs only
##  then per k, velveth -reuse_Sequences in ksubdir/  # this works
# sp2 = cgb = 410 ; sp3 = nclong = 550 ; sp1 = ncshort = 300
kseqdir=vel$dv${subset}_seq
$velbin/velveth $kseqdir 27 -fasta.gz \
    -shortPaired $datasub*.sifa2.gz \
    -shortPaired2 $datasub*.fa2.gz  \
    -shortPaired3 $datasub*.lifa2.gz \
    -short $datasub*.fa1.gz  -long $datasub.longfa.gz \
    -noHash

for k in $kset;  do { 

ksubdir=vel$dv${subset}_$k
echo "#.. start velrun $ksubdir : `date`"

mkdir $ksubdir
ln -s ../$kseqdir/Sequences $ksubdir/
$velbin/velveth $ksubdir $k  -reuse_Sequences

#. echo "velvetg : `date`" >> $notef
#. ls -l $ksubdir >> $notef
#NOT.vopts/not used before# $velbin/velvetg $ksubdir $vopts -read_trkg yes ; 
$velbin/velvetg $ksubdir -read_trkg yes ; 

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

echo "DONE `date` " >> $notef
mv $notef $donef
