#! /bin/bash
### env  subd=bamq5treat qsub -q batch bam2strand.sh
#PBS -N bam2strand
#... -A ind114
#PBS -l nodes=1:ppn=16,walltime=22:55:00
#PBS -o bam2strand.$$.out
#PBS -e bam2strand.$$.err
#PBS -V

## bam2strandfix.sh : fix failed steps
ncpu=16
subd=bamq5treat

datad=$HOME/scratch
workd=$datad/chrs/daphmag/
rund=$workd/rnas/$subd
dgenomesize=$workd/genome/dmagna20100422assembly.fa.fai

sdir=$HOME/bio/evigene/scripts/rnaseq
# samtools=$HOME/bio/bin/samtools
module add samtools

cd $rund

# first, only tCO : control
# treats=`ls -d {HS,ND}tCOc[XI]`
# bams=`ls {HS,ND}tCOc[XI]/*dmag2.bam`
# run2: all ND data
treats=`ls -d NDt{CD,NC,P5,PB,UV}cX`
bams=`ls NDt{CD,NC,P5,PB,UV}cX/*dmag2.bam`

echo "start bam2strand.$subd : `date`"  

# steps 1,2,3  foreach bam in bams; 
#    samtools sort -n chr.bam > names.bam; 
#    samtools view -f 0x2 names.bam | perl tostrands : (fwd,rev,noo).sam
#    samtools view -1 -t xxx.size -o part.bam  part.sam : (fwd,rev,noo)

cat <<-'EOP' > tostrands.pl
 BEGIN{ $f=$ENV{fn} or die; open(F,">$f.fwd.sam"); open(R,">$f.rev.sam"); open(N,">$f.noo.sam");} 
 while(<>) { ($d)=split; ($o)= m/\tXS:A:(.)/; if($d eq $ld) { if($o or $lo) { $oc=$lo||$o; 
  if($oc eq "+") { print F $ll,$_; } elsif($oc eq "-") { print R $ll,$_; } } else { print N $ll,$_; } }
  ($ll,$ld,$lo)= ($_,$d,$o); }
 END{ close(F); close(R); close(N); }
EOP

### done this step 1-2
#. i=0; 
#. for bam in $bams; do {
#.   nam=`echo $bam | sed 's/.bam//; s/-dmag2//;'`
#.   ( samtools sort -n  $bam  $nam.names; 
#.     samtools view -f 0x2 $nam.names.bam | env fn=$nam perl tostrands.pl ) &  
#.   i=$(( $i + 1 ))
#.   if [ $i -ge $ncpu ]; then wait; i=0; fi
#. } done
#. wait
#.
#. echo "bam2strand done sortn, tostrands: `date` "

## FIX here..
i=0; 
for tdir in $treats; do {
  nam=`basename $tdir | sed 's/^/gr/;'`
  samstr=`ls $tdir/*.{fwd,rev,noo}.sam`
  for sam in $samstr; do {
    nam=`echo $sam | sed 's/.sam//;'`
    samtools view -1 -t $dgenomesize -o $nam.bam $sam &
    i=$(( $i + 1 ))
    if [ $i -ge $ncpu ]; then wait; i=0; fi
  } done
} done
wait

echo "bam2strand done sam2bam: `date` "

# steps 4,5  .. merge parts into treatment groups
# 4. samtools merge -1 -n groupj.fwd.bam  part1.fwd.bam part2.fwd.bam
# 5. samtools sort groupj.fwd.bam groupj.chrs.fwd

i=0; 
for tdir in $treats; do {
  nam=`basename $tdir | sed 's/^/gr/;'`
  for ord in fwd rev noo; do {
    parts=`ls $tdir/*.$ord.bam`
    groupf=$nam.$ord
    ( samtools merge -1 -n $groupf.bam  $parts;
      samtools sort  $groupf.bam  $groupf.chr ) &
    i=$(( $i + 1 ))
    if [ $i -ge $ncpu ]; then wait; i=0; fi
  } done  
} done
wait

# then bam2bed.sh for each groupf.chr.bam

echo "end bam2strand: `date` "

exit

#..... change this to split bam into stranded/spliced parts
#  for each bamfile:
#   1. samtools sort -n part1.chrs.bam > part1.names.bam
#   2. samtools view -f 0x2 part1.names.bam | perl .. split mated pairs into fwd.sam, rev.sam, noo.sam
#   3. samtools view -1 -t xxx.size part1.rev.sam -o part1.rev.bam
#   .. merge parts into treatment groups
# 4. samtools merge -1 -n groupj.fwd.bam  part1.fwd.bam part2.fwd.bam
#   .. sort again into chr order, keep name order for other uses?
# 5. samtools sort groupj.fwd.bam groupj.chrs.fwd

