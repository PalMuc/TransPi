#!/bin/tcsh
# bam2bw.sh
#  bam2bw.sh $dmag2/genome/dmagna20100422assembly.fa.fai [134]*.bam

set dolog=1
## fixme bpath
set bpath=/bio/bio-grid/mb/bin

if ( $?dgenomesize == 0 ) then 
  set dgenomesize=$1; shift; 
endif
if ( $?drna == 0 ) then 
  set drna=$1; shift; 
endif
if ( $drna == "" || $dgenomesize == "" ) then
  echo "usage: bam2bw.sh /path/to/genome.chr_size.tab rnaseqlib4{.bam}"; exit;
endif

while ($drna != "") 
  echo "# bam2bw chrsize=$dgenomesize bam=$drna"
  set brna=`basename $drna .bam`

  if( -f $brna.bw ) then
    echo "# exists $brna.bw"; 
  else

    $bpath/samtools pileup $drna | env dolog=$dolog gs=$dgenomesize perl -ne \
'BEGIN{$dolog=$ENV{dolog}; open(S,$ENV{gs});while(<S>){($r,$n)=split; $rs{$r}=$n; } close(S);}\
($r,$b,$xn,$c)=split"\t"; if($r eq $lr and $c == $lc and $le == $b-1 ) { $le=$b; } \
else { putb($lr,$lb-1,$le,$lc);  $le=$lb=$b; }  \
($lr,$lc)= ($r,$c); END{putb($lr,$lb-1,$le,$lc);}  \
sub putb{ $v=pop(@_); if($_[0] and $_[2] > $_[1] and $v >= 1){ \
my $s=$rs{$_[0]}||0; return if($_[2] > $s or $_[1] >= $s); \
if($dolog){$v=sprintf("%.2f",log($v));} print join("\t",@_,$v),"\n"; }  }' \
  > $brna.bed

    # nasty fails when read beyond end if dgenomesize.. filter
    $bpath/bedGraphToBigWig $brna.bed $dgenomesize  $brna.bw

    if( -f $brna.bw ) then
      echo /bin/rm $brna.bed
    endif

    # if( ! -f $brna.bam.bai ) then echo $bpath/samtools index $brna.bam endif
  endif

  set drna=""
  if( $# > 0 ) then
    set drna=$1; shift;
  endif

end

