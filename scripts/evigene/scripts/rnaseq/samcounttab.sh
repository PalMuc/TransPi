#!/bin/tcsh
#/export/udisk3/work/daphplx/rnas/anreads  makegtab.sh

set dpx2=/export/udisk3/work/daphplx
set libs=(cCon cCha mCon mCad mArs mCop mNic mZin)
# cds utr
set onam=exonutr
set gnam=geneutr
set feats=exonutr19

## retest
#foreach lb ( 1 3 4 5 )
foreach lb ( 3 4 )
  set gp=$libs[$lb]
  
  samtools view ../bams/Manak_20101102_lib$lb.bam | \
  $dpx2/scripts/samcount.pl -over $feats.gff -group $gp -in stdin \
  > $onam.$gp.utab  &

end

wait
exit
#.....

foreach lb ( 2 6 7 8 )
  set gp=$libs[$lb]

  samtools view ../bams/Manak_20101102_lib$lb.bam | \
  $dpx2/scripts/samcount.pl -over $feats.gff -group $gp -in stdin \
  > $onam.$gp.utab  &

end

wait

# then sort exoncds.$gp.utab > exoncds.$gp.tab
# and 
#  cat exoncds.$gp.tab | perl -ne\
#  'if(/^feature/){print; next;} ($x,$c)=split; ($g=$x)=~s/[xcu]\d+$//; \
#  if($g ne $lg){ print "$lg\t$gc\n" if($lg); $gc=0;} $gc+=$c; $lg=$g; \
#  END{print "$lg\t$gc\n";} ' > genecds.$gp.tab


