#!/bin/tcsh
# rrnafilter.sh : remove rrna hyperexpress reads from sams

set workd=/export/udisk2/work/aphid/

foreach bam (*.bam)
  set nam=`basename $bam .bam`
  if ( -f $nam.clean.sam ) continue
  if ( -f $nam.clean.sam.gz ) continue
  echo "# Filter rrna in $bam "
  samtools view -F 4 $bam | \
   $workd/scripts/samfiltergff.pl -hardclip -expandover=50 -mate \
    -in stdin -over $workd/misc/repeat_rrna.gff \
    -filter $nam.rrna.sam -notfilter $nam.clean.sam \
    >& log.samflt.$nam

end

