#! /bin/bash
### qsub -q batch bam2strandv2.sh
#PBS -N bam2strandv2
#PBS -A ind114
#PBS -l nodes=1:ppn=16,walltime=22:55:00
#PBS -o bam2strandv2.$$.out
#PBS -e bam2strandv2.$$.err
#PBS -V

# simplified version for less complex data than daph.mag.
# this version: both proper pairs and imperfect read pairs
# bam2strandimp version: pull imperfect mates in gaps;

ncpu=16
groupnam="fungr1" # get from bams
subd=strasm1
# SPPSUF='-dmag2'
SPPSUF='-kfish2'

datad=$HOME/scratchn
workd=$datad/chrs/kfish/
rund=$workd/rnas/$subd
## BAD name here..
dgenomesize=$workd/genome/kfish2asm.fa.count

sdir=$HOME/bio/evigene/scripts/rnaseq
# samtools=$HOME/bio/bin/samtools
# module add samtools
export PATH=${PATH}:$HOME/bio/bin/

mkdir $rund
cd $rund
# temp file cleanup..
mkdir tmpfiles
# bams=`ls *.bam`
bams=`ls $workd/rnas/bam1fungr/SRR*.bam`

echo "start bam2strand.$subd : `date`"  

cat <<-'EOP' > tostrands.pl
 BEGIN{ $f=$ENV{fn} or die; open(F,">$f.fwd.sam"); open(R,">$f.rev.sam"); open(N,">$f.noo.sam");} 
 while(<>) { ($d)=split; ($o)= m/\tXS:A:(.)/; if($d eq $ld) { if($o or $lo) { $oc=$lo||$o; 
  if($oc eq "+") { print F $ll,$_; } elsif($oc eq "-") { print R $ll,$_; } } else { print N $ll,$_; } }
  ($ll,$ld,$lo)= ($_,$d,$o); }
 END{ close(F); close(R); close(N); }
EOP

cat <<-'EOX' > toimpairstrands.pl
 BEGIN{ $f=$ENV{fn} or die; open(F,">$f.fwdimp.sam"); open(R,">$f.revimp.sam"); open(N,">$f.nooimp.sam");} 
 while(<>) {
  my($d,@v)=split"\t"; next if(($v[0] & 12) == 12); ($o)= m/\tXS:A:(.)/;
  if($d eq $ld) { 
    if(bothmap($lv[0],$v[0])) { putimp( $lo||$o, $ld,\@lv,$d,\@v); }  
    if(swapmateloc(\@lv,\@v)) { putimp( $lo||$o, $ld,\@lv,$d,\@v); }  
  } 
  ($ld,$lo,@lv)=($d,$o,@v); 
 }
 sub bothmap { my($lf,$vf)=@_; return ($lf & 4 or $vf & 4)?0:1; } 
 sub putimp { my($oc,$ld,$lv,$d,$v)=@_; my $ll=join"\t", $ld,@$lv; my $vv= join"\t", $d,@$v;
  unless($oc){ print N $ll,$vv;} elsif($oc eq "+"){ print F $ll,$vv;} elsif($oc eq "-"){ print R $ll,$vv;} 
 }
 sub swapmateloc { my($l,$v)=@_; for my $r ($l,$v) { my($c,$cb,$mc,$mb)=@{$r}[1,2,5,6]; 
  if($mc eq "*") {} elsif($mc eq "=") { @{$r}[1,2,5,6]=($c,$mb,$mc,$cb); } 
  else { @{$r}[1,2,5,6]=($mc,$mb,$c,$cb); } 
  } return 1;} 
 END{ close(F); close(R); close(N); }
EOX

#......................
# 2. bam2strand; 2a/b merged
# combine steps 2,3 per input bam ? no, step3 has n-bam x 6 parts
# have enough cpu to do perf/imperf separately? but need sort -n for both

# 1. name sort : ~ 45 min for 4G bam / 80M reads
i=0; 
for bam in $bams; do {
  nam=`basename $bam .bam | sed "s/$SPPSUF//;"`
  samtools sort -n  $bam  $nam.names &
  i=$(( $i + 1 ))
  if [ $i -ge $ncpu ]; then wait; i=0; fi
} done
wait
echo "bam2strand done sortn: `date` "

# 2a,b. separate strand parts
i=0; 
for nabam in *.names.bam; do {
  nam=`basename $nabam .names.bam`
  ( samtools view -f 0x2 $nabam | env fn=$nam perl tostrands.pl ) &
  i=$(( $i + 1 ))
  ( samtools view -f 0x1 -F 0x2 $nabam | env fn=$nam perl toimpairstrands.pl ) &  
  i=$(( $i + 1 ))
  if [ $i -ge $ncpu ]; then wait; i=0; fi
} done
wait
echo "bam2strand done tostrands: `date` "

# 3. sam2bam
## FAIL here ??
i=0; 
for sam in *.{fwd,rev,noo}*.sam; do {
  nam=`basename $sam .sam`
  samtools view -1 -t $dgenomesize -o $nam.bam $sam &
  i=$(( $i + 1 ))
  if [ $i -ge $ncpu ]; then wait; i=0; fi
} done
wait
echo "bam2strand done sam2bam: `date` "

mv *.{fwd,rev,noo}*.sam tmpfiles/

# steps 4,5  .. merge parts into 6 strand groups; location sorted
# FIXME: add samtools index at end..
# FIX2: make sure ls *.ord.bam doesnt get groupnam.ord.bam; use parts prefix?
#       or above:  -o $nam.$ord.part.bam
i=0; 
for ord in fwd rev noo fwdimp revimp nooimp; do {
  parts=`ls *.$ord.bam`
  groupf=$groupnam.$ord
  ( samtools merge -1 -n $groupf.bam  $parts; mv $parts tmpfiles/ ; \
  samtools sort  $groupf.bam  $groupf.chr; mv $groupf.bam tmpfiles/ ; \
  samtools index $groupf.chr.bam ) &
  i=$(( $i + 1 ))
  if [ $i -ge $ncpu ]; then wait; i=0; fi
} done  
wait


# then bam2bed.sh for each groupf.chr.bam
echo "end bam2strand: `date` "

#........................................
# steps 1,2,3  foreach bam in bams; 
#    samtools sort -n chr.bam > names.bam; 
#    samtools view -f 0x2 names.bam | perl tostrands : (fwd,rev,noo).sam
#    samtools view -1 -t xxx.size -o part.bam  part.sam : (fwd,rev,noo)
#* Modify here for impaired: pull improperly paired mates in gaps;
#  step 2i: samtools view -F 0x2 names.bam | perl toimpairedstrands.pl ..

# 2a: samtools view -f 0x2 $nam.names.bam | env fn=$nam perl tostrands.pl
## see also below 2a. changes for 2i.
# i=0; 
# for bam in $bams; do {
#   nam=`echo $bam | sed 's/.bam//; s/-dmag2//;'`
#   ( samtools sort -n  $bam  $nam.names; 
#     samtools view -f 0x2 $nam.names.bam | env fn=$nam perl tostrands.pl ) &  
#   i=$(( $i + 1 ))
#   if [ $i -ge $ncpu ]; then wait; i=0; fi
# } done

# 2i: samtools view -F 0x2 $nam.names.bam | env fn=$nam perl toimpairstrands.pl
# 3. sam2bam: samtools view -1 -t $dgenomesize -o $nam.bam $sam
# 4. samtools merge -1 -n groupj.fwd.bam  part1.fwd.bam part2.fwd.bam
# 5. samtools sort groupj.fwd.bam groupj.chrs.fwd
