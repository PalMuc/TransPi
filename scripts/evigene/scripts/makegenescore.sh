#!/bin/tcsh
# makegenescore.sh

foreach pred ( dmag2_epir3 dmag2_epir2 dmag_ep24 )
  foreach blp ($pred-*.blastp.gz)
    set runset=`basename $blp .blastp.gz`
    if( -f $runset.genescore ) continue;
    if( $blp =~ "$runset-$runset.blastp.gz" ) continue; # _self
    if( $blp =~ "$runset-self.blastp.gz" ) continue; # _self
  
    echo $runset
    gzcat $runset.blastp.gz | grep -v '^#' | perl -ne'($q,$t,$b)=(split)[0,1,-1]; print "$t\t$b\t$q\n";' |\
    sort -k1,1 -k2,2nr | perl -ne'($g,$b,$q)=split; print unless($did{$g}++);' > $runset.genescore

  end

  foreach blp ($pred-self.blastp.gz)
    set runset=`basename $blp .blastp.gz`
    if( -f $runset.genescore ) continue;

    echo $runset
    gzcat $blp | grep -v '^#' | perl -ne\
'($q,$t,$b)=(split)[0,1,-1]; $bs=0 if($lq ne $q); \
if($t eq $q) { $bs=$b unless($bs); if($v){ $v=~s/0$/$bs/; print $v; $v=""; } } \
else { unless($did{$q}++) { $v="$q\t$b\t$t\t$bs\n"; if($bs>0){ print $v; $v=""; } } } $lq=$q; ' \
 | sort -k1,1 > $runset.genescore

  end
end
