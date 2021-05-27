#! /bin/bash
### qsub -q batch gsnap2mag4trmap.sh
#PBS -N gsnaptrmap
#... -A ind114
#PBS -l nodes=1:ppn=8,walltime=21:55:00
#PBS -o gsnatrmap.$$.out
#PBS -e gsnatrmap.$$.err
#PBS -V

trdb=gmapgoodalt3m8veltribest2
trdir=trsets
odir=mag4trmapgs

# per trdb x 6 rdsets
ncpu=8
datad=$HOME/scratch
workd=$datad/chrs/daphmag/rnas/
bindir=$HOME/bio/gmap1204/bin

lreads=$workd/fastq4/Dman21_TTAGGC_L006_R1_001.fastq.gz
#below#rdnam=rdDman21

## for tr/cds mapping: no introns, no gmap, no failed reads
# optsnp="--use-snps=$snpd --show-refdiff --print-snps"
# optreads=" --antistranded-penalty=1 --filter-chastity=both "

optout="-A sam --no-sam-headers"
#  --split-output=$rdnam.gsn"
## maybe: -N 1 for trmap? trs have sometimes retained introns
optcds="-N 0 --gmap-mode=none "
snapopt1="$optcds $optout"

cd $workd/
mkdir $odir
tnam=`basename $trdb .fa | sed 's/gmap//;'`

# for lreads in fastq3/*_1.sanfastq; do {

rreads=`echo $lreads | sed 's/_R1/_R2/;'`
rdnam=`basename $lreads _R1_001.fastq.gz `
snapopts="$snapopt1 --gunzip"

$bindir/gsnap $snapopts --nthreads=$ncpu --split-output=$odir/$tnam-$rdnam.gs -D $trdir -d $trdb \
  $lreads $rreads  > $odir/$tnam-$rdnam.gsnap.out & 

wait

# i=0;   # use --nthreads= instead of this loop?
# while [ $i -lt $ncpu ]; do { 
# 
#   $bindir/gsnap $snapopts --split-output=$odir/$tnam-$rdnam.gs$i  --part=$i/$ncpu -D $trdir -d $trdb \
#     $lreads $rreads  > $odir/$tnam-$rdnam.gsnap$i.samu &
#    i=$(( $i + 1 ))
# 
# } done
# wait

# } done  # infile
# wait

exit


# -----------------

#! /bin/bash
### qsub -q normal gsnapset.sh
#PBS -N gsnp4
#PBS -l mem=64gb,nodes=1:ppn=32,walltime=9:55:00
#PBS -o gsnp1.$$.out
#PBS -e gsnp1.$$.err
#PBS -V

# 280 rna files .. split into 10 ? parts/nodes
ncpu=32
npart=32
#ncpu=64 npart=64
jstart=0
jend=$(( $jstart + $ncpu ))

datad=/N/dc/scratch/$USER
workd=$datad/chrs/cacao
bindir=$HOME/bio/bin
gmapd=$HOME/bio/gmap118
rund=$workd/rnas
gdb=gmap118
gdbpath=$workd/genome/

# run5 on genome
snpd=snp5pub3hcds
dgenome=cacao11allasm ; gtag=mars11
dgenosize=$dgenome.chr_size.txt

# .. run thru gene.cds.fa only
# snpd=snp10pub3hcds1
# dgenome=genes3h_goodcds ; gtag=cds3h

# run3, cg3 data, 4 snp groups (all in 1 db, uniq ids)
# snpd=snp4pub3hcds
# dgenome=genes3h_goodcds ; gtag=cds3h

#1 snapopt1="--use-snps=$snpd -N 1 --quality-protocol=illumina "
#2 snapopt1="--use-snps=$snpd --show-refdiff --print-snps -N 1 --quality-protocol=illumina "
# this works: --print-snps  --show-refdiff  << print-snps lists snps: nn@idsnp,...
#......
# cgb2 reads: stranded??  --antistranded-penalty=1
# cgb2 reads are diff format.. not same qual, -stranded, lacking pair /1 /2 ..
# also has Ill chastity string : --filter-chastity= none,either,both << works
#3
# optsnp="--use-snps=$snpd --show-refdiff --print-snps"
# optreads=" --antistranded-penalty=1 --filter-chastity=both "
# snapopt1="-N 1 $optsnp $optreads"

## .. run thru gene.cds.fa only
## for cds mapping: no introns, no gmap, no failed reads
optsnp="--use-snps=$snpd --show-refdiff --print-snps"
optreads=" --antistranded-penalty=1 --filter-chastity=both "
optcds="-N 0 --nofails --gmap-mode=none "
snapopt1="$optcds $optsnp $optreads" 

#5 for genome map, .sam out, with introns, no snps print
# .. should add known intron database
# .. test for chasity: both, either, none
optsnp="--use-snps=$snpd "
optreads="--filter-chastity=both "
# optreads="--filter-chastity=none "
optgenome="-N 1 --nofails " 
  #? --gmap-mode=none "
snapopt1="$optgenome $optsnp $optreads" 

notef=$workd/rnas/gsnapset.$$.RUNNING
donef=`echo $notef | sed 's/RUNNING/DONE/'`

#fastdir=fa3uf18: uf18l1_GTCCGC_1.fastq.gz uf18l1_GTCCGC_2.fastq.gz

cd $workd/rnas/
if [ "X$fastdir" = "X" ]; then
  fastdir="fastq"
fi

## tie outdir to fastdir or not?
outdir=$workd/rnas/snp4o$gtag
outdir=`echo $fastdir | sed "s/^fa//; s/^/snp4o/;"`
#..

touch $notef
echo "START " >> $notef
echo `date`   >> $notef

mkdir $outdir
ls -l $fastdir >> $notef

#... START LOOP rnain : limit to suffix .fq, .fastq, sequence.txt ..

for rnain in $fastdir/*{_1.fastq,.1.fastq,_1_sequence.txt}.gz ;  do 
{ 
  cd $rund/
  # mkdir $rund/sams
  if ! test -f "$rnain" ; then continue; fi
  
  drna=`basename $rnain .gz | sed 's/.gz//; s/.fastq//; s/.fq//; s/_sequence.txt//;'`
  outna=$drna-$gtag
  echo "START gsnap : $rnain to $outna.bam" >> $notef
  echo `date`  >> $notef

  snapopts="--gunzip $snapopt1"
  
  rnain2=`echo $rnain | sed 's/1.fastq/2.fastq/; s/1_sequence/2_sequence/;'`
  if test -f "$rnain2" ; then
    inset="$rnain $rnain2"
  else 
    inset=$rnain
  fi

  # to narrow down break point, use i, j counters, 
  # i=1,ncpu, j=jstart,jend  subset of npart
  i=0;  
  j=$jstart;
  while [ $i -lt $ncpu ]; do { 
  
   $gmapd/bin/gsnap $snapopts  -A "sam"  --part=$j/$npart \
    -D $gdbpath/$gdb -d $dgenome $inset > $outdir/$drna.gsnap$i.samu &
  
   i=$(( $i + 1 ))
   j=$(( $j + 1 ))
  }
  done
  echo "wait gsnap : $drna" >> $notef
  
  wait

  echo "DONE gsnap : $drna" >> $notef
  du -h $outdir >> $notef
  ls -l $outdir >> $notef
  #..................
  ##...... sam sort/merge........

  samset=$drna
  #? cd $rund/sams/
  cd $outdir
  
  i=0; while [ $i -lt $ncpu ]; do {
    isamu="$samset.gsnap$i.samu"
    inam=`basename $isamu .samu`
    if ! test -f $isamu ; then  echo "missing $isamu"; continue;  fi
  
    ( $bindir/samtools view -u -t $gdbpath/$dgenosize  $isamu | $bindir/samtools sort - $inam.sort; \
      rm $isamu ) &
  
    i=$(( $i + 1 ))
  }
  done
  
  wait

  echo "DONE bam sort" >> $notef
  du -h >> $notef
  #..................
  
  $bindir/samtools merge $outna.bam $samset*.sort.bam
  echo "DONE merge to $outna.bam" >> $notef
  ls -l *.bam  >> $notef
 
} 
done

#... END LOOP rnain
echo "DONE " >> $notef
echo `date`  >> $notef
mv $notef $donef
