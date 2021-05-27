#!/bin/tcsh
# scripts/anaugmap.sh  : annotate prediction.gff

# 4.1. evidence annotate by base overlaps: protein, est, rnaseq, intron, tiletar 
#  redo scoring: -pct 10 when using markbase; otherwise lose partial real scores  
#  2010.10: add -strand overlapfilter options for most; note many ESTs not stranded;
#    PASA has bogus +strand for unknown

set pcto=10
set anv=an1

#NOT: | scripts/overlapfilter -strand -pass 'exon' -pct $pcto -act markidbase -mark ref \
#     -in stdin -over refseq/refseq-genes.gff.gz \

# if ( 0 == 1 ) then
# also genes/dmag_ep24augmap2.gff.gz 
foreach augmapgff (genes/dmag*augmap2.gff.gz)
#foreach augmapgff (augustus/*-augmap.gff.gz)

  set grp=`basename $augmapgff .gff.gz`
  if( -f genes/$grp.$anv.gff.gz ) continue
  if( -f genes/$grp.$anv.gff ) continue

  echo "$grp : $anv"
  gzcat $augmapgff \
| scripts/overlapfilter -strand -pass CDS -pct $pcto -act markidbase -mark pro \
    -in stdin -over prot/protein_uniq.gff.gz \
| scripts/overlapfilter -nostrand -pass 'exon,cDNA_match' -pct $pcto -act markidbase -mark pasa \
    -in stdin -over est/pasa_assemblies.gff.gz \
| scripts/overlapfilter -strand -pass 'exon,HSP' -pct $pcto -act markbase -mark est \
   -in stdin -over est/est_uniq.gff.gz \
| scripts/overlapfilter -strand -pass 'exon' -pct $pcto -act markidbase -mark rseq \
   -in stdin -over rnas/rnaseq_uniq.gff.gz \
| scripts/overlapfilter -intron2splice -pass 'exon,intron' -act markid -midtype scoresum \
    -mark intr -in stdin -over intron/intron.gff.gz \
| scripts/overlapfilter -strand -pass 'exon,transposon' -pct $pcto -act markbase -mark terepeat \
   -in stdin -over misc/transposon.gff.gz \
> genes/$grp.$anv.gff 

end
#endif

exit

#..... can we add equiv cDNAgene score: 100 for perfect exon matches; lower for too much/little
# epasa3/pasatrain_genes.best1.gff.gz : add above overfilt cgene score? then accum per mRNA

#..... add this homology genescore annot ; ho3= fixed in overbest; should be hbest=
#..... both plant_other and plant_arab scores/ids? : harabid=

cd genes/
foreach predf (*-augmap.$anv.gff)
  # set pred=`basename $predf -augmap.$anv.gff`
  set pred=`basename $predf .$anv.gff`
  set prednam=`echo $pred | sed -e 's/aphid2_/AUG/;'`
  setenv gpre $prednam 

  if( -f $pred-augmap.${anv}ho.gff ) continue
  
  if( -e aaeval/${pred}-*.genescore ) then
  cat aaeval/${pred}-*.genescore \
  | sort -k1,1 -k2,2nr | sed -e's/^/gscore /' | cat - ${predf} | perl -ne \
'if(/^gscore\s+(\S+)\t(\S+)\t(\S+)$/){ ($g,$h)=($1,"$2/$3"); $gb{$g}=$h  unless($gb{$g});\
if(/drosmel/){ $arb{$g}=$h unless($arb{$g}); } }\
else { if(/\tmRNA/) { ($g)=m/ID=(\w+)/; $b=$gb{$g}; s/$/;ho3=$b/ if $b; \
$b=$arb{$g}; s/$/;hdros=$b/ if $b; } print; }' \
  > $pred-augmap.${anv}ho.gff
  endif

end

exit
#...........................

