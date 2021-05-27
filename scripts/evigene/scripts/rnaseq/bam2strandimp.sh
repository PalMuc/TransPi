#! /bin/bash
### qsub -q batch bam2strandimp.sh
#PBS -N bam2strand
#... -A ind114
#PBS -l nodes=1:ppn=16,walltime=22:55:00
#PBS -o bam2strand.$$.out
#PBS -e bam2strand.$$.err
#PBS -V

ncpu=16
subd=bamq5treat
# bam2strandimp version: pull imperfect mates in gaps;

datad=$HOME/scratch
workd=$datad/chrs/daphmag/
rund=$workd/rnas/$subd
dgenomesize=$workd/genome/dmagna20100422assembly.fa.fai

sdir=$HOME/bio/evigene/scripts/rnaseq
# samtools=$HOME/bio/bin/samtools
module add samtools

cd $rund

# first, only tCO : control
treats=`ls -d {HS,ND}tCOc[XI]`
bams=`ls {HS,ND}tCOc[XI]/*dmag2.bam`
# run2: all ND data
# treats=`ls -d NDt{CD,NC,P5,PB,UV}cX`
# bams=`ls NDt{CD,NC,P5,PB,UV}cX/*dmag2.bam`

echo "start bam2strand.$subd : `date`"  

# steps 1,2,3  foreach bam in bams; 
#    samtools sort -n chr.bam > names.bam; 
#    samtools view -f 0x2 names.bam | perl tostrands : (fwd,rev,noo).sam
#    samtools view -1 -t xxx.size -o part.bam  part.sam : (fwd,rev,noo)
#* Modify here for impaired: pull improperly paired mates in gaps;
#  step 2i: samtools view -F 0x2 names.bam | perl toimpairedstrands.pl ..

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

i=0; 
for bam in $bams; do {
  nam=`echo $bam | sed 's/.bam//; s/-dmag2//;'`
  ( samtools sort -n  $bam  $nam.names; \
    samtools view -F 0x2 $nam.names.bam | env fn=$nam perl toimpairstrands.pl ) &  
  i=$(( $i + 1 ))
  if [ $i -ge $ncpu ]; then wait; i=0; fi
} done

wait

echo "bam2strand done sortn, tostrands: `date` "

i=0; 
for tdir in $treats; do {
  nam=`basename $tdir | sed 's/^/gr/;'`
  #2a: samstr=`ls $tdir/*.{fwd,rev,noo}.sam`
  samstr=`ls $tdir/*.{fwd,rev,noo}imp.sam`
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
  #2a: for ord in fwd rev noo; do 
  for ord in fwdimp revimp nooimp; do 
  {
    parts=`ls $tdir/*.$ord.bam`
    groupf=$nam.$ord
    ( samtools merge -1 -n $groupf.bam  $parts; \
      samtools sort  $groupf.bam  $groupf.chr ) &
    i=$(( $i + 1 ))
    if [ $i -ge $ncpu ]; then wait; i=0; fi
  } done  
} done
wait

# then bam2bed.sh for each groupf.chr.bam
echo "end bam2strand: `date` "

