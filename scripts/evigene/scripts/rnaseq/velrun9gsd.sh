#!/bin/bash
#==> velrun9g.sh <==
##  env subset=sc8 qsub -q normal velrun9.sh
#PBS -A ind114
#PBS -N velrun 
#PBS -l nodes=1:ppn=16,walltime=21:55:00
#.. lower time=22:  all vel6 done 3kmers in < 12hr
#PBS -o velrun.$$.out
#PBS -e velrun.$$.err
#PBS -V

#.. vel9: bam selection: MAPQ 20, properpairs, notSecondary
# pt=sub9.sc6.nc1; 
# $bindir/samtools view -b -q 20 -f 0x2 -F 0x100 -o $pt-hiq.bam $pt.bam
# $bindir/samtools sort -n $pt-hiq.bam $pt-hiqsn >& log.s$pt &
#.. note this results in many 2x? MORE reads than prior sam2velv.pl using align qual and end clipping..
#..
#.. vel9: problem appears to be lowqual data, test w/ higher qual subset, .bam format
#.. v91: k35 ok, k23 failed outofment (64G)  at velvetg; using more mem than fasta variant. more reads?
#.. v92: remove .ncgr bams, keep .cgb only: tr errors may be in part due to strain snps ..
#.. v93: remove velvetg --opts; not used before in vel1..5;

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
dv=94
# vel5,6,7: datasub=sub3.$subset
datasub=sub9.$subset

ncpu=8
export OMP_NUM_THREADS=$ncpu
export OMP_THREAD_LIMIT=$ncpu

#tres#datad=/phase1/$USER
datad=/oasis/scratch/$USER/temp_project

#tres#rund=/scratch/$USER/$PBS_JOBID
rund=/oasis/scratch/$USER/$PBS_JOBID

#velbin=$datad/bio/velvet/bin
velbin=$HOME/bio/velvet/bin
workd=$datad/chrs/cacao/rnas/vel9
homed=$HOME/work/

## ** FIXME FIXME : velvet wants read len added to inslen, for cgb 106b pairs, full inslen = 412
## ** FIXME222222 : .bams need to be ID-sorted, so mates follow each other, like .fa2: samtools sort -n
#v6
#vopts="-ins_length 450 -ins_length_sd 80 -ins_length2 200 -ins_length2_sd 40"
#oopts="-min_trans_lgth 100 -ins_length 450 -ins_length_sd 80 -ins_length2 200 -ins_length2_sd 40"
#v7,v8,v9
#vopts="-ins_length 200 -ins_length2 400"
#oopts="-min_trans_lgth 100 -ins_length 200 -ins_length2 400"
#v9c : cgb ins is 400 (212b-rd + 200 ins); nc short 108-rd+ 200ins = 300; nc long=108r+400i=550

#NO.off.v93#vopts="-ins_length 410 -ins_length2 550"
vopts=""
#v93: oopts="-min_trans_lgth 100 -min_pair_count 6 -cov_cutoff 6  -ins_length 410"
# ^ v93 new opts -minpair,cutoff not helpful vs v92
#v94
oopts="-min_trans_lgth 100 -ins_length 380"

#v6#kset="39 27 23"
#v7#kset="33 21"
#v8; k35,23 > v9 k35,25
#v93 kset="35 25"
#v94
kset="63 53 45 39 29"
## FAIL: vel94sc6_53   ok: vel94sc6_63 vel94sc6_45

ksubdir=vel$dv${subset}

notef=$workd/$ksubdir.$$.RUNNING
donef=`echo $notef | sed 's/RUNNING/DONE/'`

touch $notef
echo "START " >> $notef
echo `date`   >> $notef
#..... run loop for vel kmer steps
shopt -s nullglob

for k in $kset;  do { 

ksubdir=vel$dv${subset}_$k
echo "#.. start velrun $ksubdir : `date`"

#. du -h >> $notef
#. echo "velveth $ksubdir" >> $notef

#v9: only -shortPaired .bam and -long .longfa.gz
$velbin/velveth $ksubdir $k -bam -shortPaired $datasub.*.bam  -fasta.gz -long $datasub.longfa.gz

# $velbin/velveth $ksubdir $k -fasta.gz \
#    -shortPaired $datasub.*.sifa2.gz $datasub.*.fa2.gz  \
#    -shortPaired2 $datasub.*.lifa2.gz \
#    -short $datasub.*.fa1.gz   -long $datasub.longfa.gz 

#. echo "velvetg : `date`" >> $notef
#. ls -l $ksubdir >> $notef
$velbin/velvetg $ksubdir $vopts -read_trkg yes 

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

