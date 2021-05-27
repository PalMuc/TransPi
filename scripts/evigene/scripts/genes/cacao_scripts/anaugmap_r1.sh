#!/bin/tcsh
# 4.1. evidence annotate by base overlaps: protein, est, rnaseq, intron, tiletar 
## .an3 add pasa= evidence
## FIXME, redo scoring for -pct 10 when using markbase; otherwise lose partial real scores  

set pcto=10

# if ( 0 == 1 ) then
# foreach augmapgff (augustus/*-augmap.gff.gz)
foreach augmapgff (genes/cacao9_fgenesh_fix.gff.gz)
# foreach augmapgff (augustus/{cacao9_epit2,cacao9_epi4}-augmap.gff.gz)

  set grp=`basename $augmapgff .gff.gz`
  echo "$grp : an3"
  gzcat $augmapgff \
| scripts/overlapfilter -pass CDS -pct $pcto -act markidbase -mark pro -in stdin -over prot/protein_uniq.gff.gz \
| scripts/overlapfilter -pass 'exon,cDNA_match' -pct $pcto -act markbase \
    -mark pasa -in stdin -over est/pasa_assemblies.gff.gz \
| scripts/overlapfilter -pass 'exon,HSP' -pct $pcto -act markbase -mark est -in stdin -over est/est_uniq.gff.gz \
| scripts/overlapfilter -pass 'exon' -pct $pcto -act markbase -mark rseq -in stdin -over rnas/rnaseq.gff.gz \
| scripts/overlapfilter -intron2splice -pass 'exon,intron' -act markid -midtype scoresum \
    -mark intr -in stdin -over intron/intron.gff.gz \
> genes/$grp.an3.gff 

end
# endif

exit

#..... add this homology genescore annot .an4
cd genes/
foreach predf (cacao9*-augmap.an3.gff)
set pred=`basename $predf -augmap.an3.gff`
set prednam=`echo $pred | sed -e 's/cacao9_/AUG/;'`
setenv gpre $prednam 

cat aaeval/${pred}-*.genescore | sort -k1,1 -k2,2nr | \
cat - ${pred}-augmap.an3.gff | perl -ne \
'if(/^($ENV{gpre}\w+)\t(\S+)\t(\S+)$/){ $gb{$1}="$2/$3" unless($gb{$1});}\
else { if(/\tmRNA/) { ($g)=m/ID=(\w+)/; $b=$gb{$g}; s/$/;ho3=$b/ if $b; } print; }' \
> $pred-augmap.an4.gff
end

#...........................
# special case for fgenesh w/o exons, just mRNA/CDS
## drop this; use genes/cacao9_fgenesh_fix.gff.gz with exon == CDS
set augmapgff=genes/cacao1_fgenesh_m_20100901.gff.gz

  set grp=`basename $augmapgff .gff.gz`
  echo "$grp : an2"
  gzcat $augmapgff \
| scripts/overlapfilter -pass CDS -pct $pcto -act markidbase -mark pro -in stdin -over prot/protein_uniq.gff.gz \
| scripts/overlapfilter -pass 'CDS,HSP' -pct $pcto -act markbase -mark est -in stdin -over est/est_uniq.gff.gz \
| scripts/overlapfilter -pass 'CDS' -pct $pcto -act markbase -mark rseq -in stdin -over rnas/rnaseq.gff.gz \
| scripts/overlapfilter -intron2splice -pass 'CDS,intron' -act markid -midtype scoresum -mark intr \
  -in stdin -over intron/intron.gff.gz \
> genes/$grp.an2.gff

