#!/bin/tcsh
#env gmapset=*.gmap.out8.gz makegff.sh

set evigene=/bio/bio-grid/mb/evigene

foreach goz ( $gmapset )
 set ena=`basename $goz .gz | sed 's/.gmap.out.*//;'`
 set nn=`echo $ena | sed 's/.s[1-9]r.//; s/.mars11//; s/.cirad1c//; s/Assembly/a/; s/reads/er/;'`
 # echo gzcat $goz TO env src=$nn intron=-1 best=0 $evigene/scripts/gmap_to_gff.pl TO $ena.gmap.gff
 gzcat $goz | env src=$nn intron=-1 best=0 $evigene/scripts/gmap_to_gff.pl > $ena.gmap.gff
end

