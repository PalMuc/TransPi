#! /bin/bash
### makeidbatr.sh : postprocess idba-trans outputs to evg trs

evigene=$HOME/bio/evigene
xbin=$HOME/bio/exonerate/bin

## datad == outdir
if [ "X" = "X$subd" ]; then 
  if [ "X" = "X$datad" ]; then echo "missing env datad=/path/to/data "; exit -1; fi
  subd=$datad
fi
nam=`basename $subd`
pt=$nam

# $idbin/idba_tran $traopts --read $inpe --num_threads $ncpu --out $outdir
#.........................

## each transcript-kmer.fa
cat $subd/transcript-??.fa | env tag=${nam}rid perl -ne \
'if(/^>/) { s/transcript-(\d+)_/$ENV{tag}k${1}Loc/;  $ok=1;  } print if($ok); ' \
 > $nam-trs.tr 
$xbin/fastanrdb $nam-trs.tr > $nam-trs.nrtr 
$evigene/scripts/cdna_bestorf.pl -nostop -minaa=30 -aa -cdna $nam-trs.nrtr

## final contig.fa ; set kmer tag == 01 for consistency.
## dang dang --no_correct gives wacky contig.fa IDs
# >_112631
# >_112632

if [ -f $subd/contig.fa ]; then

cat $subd/contig.fa | env tag=${nam}fid perl -ne \
'if(/^>/) { s/contig-\d+_//; s/>_/>/; s/>/>$ENV{tag}k01Loc/; ($w,$rc)=m/(?:length|read_count)_(\d+)/g;
s/(length|read_count)_(\d+)/$1=$2;/g; $ok=1; $maybeok=($w>=200 and $rc>0)?1:0; } print if($ok); ' \
  > $nam-fin.tr
$evigene/scripts/cdna_bestorf.pl -nostop -minaa=30 -aa -cdna $nam-fin.tr

fi

env stat=1 span=1 $evigene/scripts/prot/aaqual.sh $nam-*.aa

if [ -f $subd/contig.fa ]; then
# compress output set
gzip --fast $subd/*.fa
gzip --fast $subd/kmer
gzip --fast $subd/{transcript-path,connection,component,align}-*
fi

#................
