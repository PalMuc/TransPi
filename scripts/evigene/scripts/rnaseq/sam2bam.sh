#!/bin/tcsh

set rnas=/bio/bio-grid/mb/rnaseq
set workd=/export/udisk3/work/daphmag

# add -u -h for view

# foreach samz ( *.sam.gz )
foreach samz ( *.sam )

  #set samf=`basename $samz .sam.gz`
  set samf=`basename $samz .sam`
  if( -f $samf.sort.bam ) continue

  #was# $rnas/bin/samtools view -b -T $workd/genome/dmagna20100422assembly.fa $samz
  $rnas/bin/samtools view -u -h -b -T $workd/genome/dmagna20100422assembly.fa $samz \
    | $rnas/bin/samtools sort - $samf.sort

end

