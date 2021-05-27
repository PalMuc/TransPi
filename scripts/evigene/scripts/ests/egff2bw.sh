#!/bin/tcsh
# from bam2bw.sh
#  bam2bw.sh $dmag2/genome/dmagna20100422assembly.fa.fai [134]*.bam
#  input here chrsize.txt  *.gff3
## modified for overlapped gff reads to non-over coverage bed
## in: est.{bean,leaf,pistil}.mars11.gff.gz   > est.bean.bed > .bw

set gsuf=".gff.gz"
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
  set brna=`basename $drna $gsuf | sed 's/.mars11//'`

  if( -f $brna.bw ) then
    echo "# exists $brna.bw"; 
  else
 
  if( ! -f $brna.bed ) then
    gzcat $drna | env gs=$dgenomesize dolog=$dolog perl -ne \
' ($r,$b,$e,$c,$o)=(split"\t")[0,3,4,5,6]; $w=1+$e-$b; $li=$b-$lb; \
if($lr ne $r) { putc($lr,$lb,$le) if($lr and $lb>0); $lie=$lb=$le=0;  } \
if($lr eq $r and $li>0) { putc($lr,$lb,$li); } \
for($i=0; $i<$w; $i++){ $c[$i]+=$c; }  ($lr,$lb,$le,$lc,$lo)=($r,$b,$e,$c,$o); \
sub putc{ my($r,$ib,$w)=@_;  my $ie=$ib; my $lc=0; \
for my $i (0..$w) { my $c= shift @c; $lc=$c if($i==0); if($c != $lc or $i==$w) { \
$ib=$lie if($ib<$lie); $ie=$ib+1 if($ie<=$ib); \
if($lc>0){ $lc=sprintf("%.2f",log(1+$lc)) if($dolog); print join("\t",$r,$ib,$ie,$lc)."\n"; } \
$lie=$ib=$ie; } $ie++; $lc=$c; } } BEGIN{$dolog=$ENV{dolog};} ' \
   > $brna.bed
  endif

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

