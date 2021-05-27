#!/bin/tcsh
# makegtab.sh

set dpx2=/export/udisk3/work/daphplx/
set libs=(cCon cCha mCon mCad mArs mCop mNic mZin)
# cds utr
set onam=exonutr
set gnam=geneutr
set feats=exonutr19

foreach lb ( 1 2 3 4 5 6 7 8 )
  set gp=$libs[$lb]
  sort $onam.$gp.utab > $onam.$gp.tab

  cat $onam.$gp.tab | perl -ne\
  'if(/^feature/){print; next;} ($x,$c)=split; ($g=$x)=~s/[xcu]\d+$//; \
  if($g ne $lg){ print "$lg\t$gc\n" if($lg); $gc=0;} $gc+=$c; $lg=$g; \
  END{print "$lg\t$gc\n";} ' > $gnam.$gp.tab

end

paste $gnam.*.tab | cut -f 1,2,4,6,8,10,12,14,16 > allgrp.$gnam.tab
paste $onam.*.tab | cut -f 1,2,4,6,8,10,12,14,16 > allgrp.$onam.tab

