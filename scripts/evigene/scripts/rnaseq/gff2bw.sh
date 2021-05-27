#!/bin/tcsh
# from bam2bw.sh
#  bam2bw.sh $dmag2/genome/dmagna20100422assembly.fa.fai [134]*.bam
#  input here chrsize.txt  *.gff3

set gsuf=".gff3"
set dolog=0
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
  set brna=`basename $drna $gsuf`

  if( -f $brna.bw ) then
    echo "# exists $brna.bw"; 
  else
 
    # cut -f1,4,5,6 $drna | grep -v '^#' > $brna.bed
    env gs=$dgenomesize dolog=$dolog perl -ne \
'BEGIN{ $dolog=$ENV{dolog};  open(S,$ENV{gs}); while(<S>){($r,$n)=split; $rs{$r}=$n; } close(S);} \
next unless(/^\w/); ($r,$b,$e,$c)=(split"\t")[0,3,4,5]; next unless($e >= $b); $e++ if($e == $b); \
$r="scaffold_10r" if($r eq "scaffold_10"); $s=$rs{$r}||0; \
next if($e>$s or $b>=$s); print join("\t",$r,$b,$e,$c),"\n"; ' \
  $drna > $brna.bed

# if($dolog){$v=sprintf("%.2f",log($v));} 
# nasty fails when read beyond end if dgenomesize.. filter
# line 10 of ccn51_Tcv1_coverage.bed: start and end coordinates the same; # They need to be at least one apart

    $bpath/bedGraphToBigWig $brna.bed $dgenomesize  $brna.bw
    if( -f $brna.bw ) then
      echo /bin/rm $brna.bed
    endif
  endif

  set drna=""
  if( $# > 0 ) then
    set drna=$1; shift;
  endif

end
