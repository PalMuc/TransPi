#!/bin/tcsh
# mustkeepgenes.sh
# for aphid2, mod for others?

##
##     ***  SEE Replacement in bestgenes_update.pl   ***
##

set lnew=genes/mustkeepdrop.new
set lorig=genes/mustkeepdrop.list

set cufdir=rnas/velmapt7/
set evigene=/bio/bio-grid/mb/evigene/

cat map/pickgene-aphid2x.results | cut -f4,5 | perl -pe \
's/gid=(m|m\w|x\w)\d+/gid=/; s/p[1-9]\t/\t/ if(/ACYPI/); s/(gpick|gid)=//g; ' \
> $lnew

# Fixme for aphid_cuf8r n=5; aphid2Trinity ids n=27;  
#   asmbl = pasa n=4; PASAgasmbl n=32; XM n=1; ACYPI n=4;
#   AUG n=232 .. all ok??  << check bestgenes set for missing AUG ids

# AUG inputs: 
#   aphid2_epi4.an7.gff.im  aphid2_epir16b.an7.gff.im  
#   aphid2_epir2.an7.gff.im  bestgenes_of4.an7g2.gff.im
# mustkeep AUG n=188, found n=132; 
#  16 missed = pasa joins: AUGepir16bs48g101t1_m7AUGepir9s48g91t1
#  40 other miss = AUGepi5, AUGepir2, AUGepir3, AUGepir9,   AUGepir10
# comm -32 mustkeep.augid mustkeep.haveaugid | grep -v '_' > mustkeep.getaugid
# gzcat aphid2_{epi5,epir3,epir9,epir10}-augmap.gff.gz | ggrep -F -f mustkeep.getaugid - > mustkeep.getaug.gff
#  ##  from aphid2_{epi5,epir3,epir9,epir10}-augmap.gff.gz
# : add to bestgenes_of4.an7g2.gff ?? or other

# $evigene/scripts/annotate_predictions.pl -noanno -vers an7s -conf=genes/evigene_aphid2.conf \
#  genes/*.gff.im rnagene/*.gff.im > & genes/bestgenes_of.an7s.log


# change this overlapfilter TO overgenedup.pl
# n same=274112 of in 344642 mRNA
if( ! -f $cufdir/aphid_rnaseq.cuf1to2.ids ) then
  $evigene/scripts/overgenedup.pl -pass 'mRNA,exon' -type exon -slopexon=8 \
  -in $cufdir/aphid_rnaseq.all27cuff8.gff.gz -over $cufdir/aphid_rnaseq.all17cuff8.gff.gz -act markid -mark cuf1 \
  | grep cuf1= | cut -f9 | sed 's/;;cuf1=/  /; s/ID=//;'  > $cufdir/aphid_rnaseq.cuf1to2.ids

endif

### convert cuff1to2 and pull those models
# cat $lnew | grep '^Gsc' | cut -f1 > $cufdir/aphid_rnaseq.all17cuff8.keepids
#
# gzgrep mRNA $cufdir/aphid_rnaseq.all17cuff8.gff.gz | ggrep -F -f $cufdir/aphid_rnaseq.all17cuff8.keepids - \
# > $cufdir/aphid_rnaseq.all17cuff8.keepmrna
#
# gzgrep mRNA $cufdir/aphid_rnaseq.all27cuff8.gff.gz | \
# $evigene/scripts/overlapfilter -strand -in stdin -over $cufdir/aphid_rnaseq.all17cuff8.keepmrna \
# -act markid -mark cuf1 -type samefeat | grep cuf1= | cut -f9 | sed 's/;;/  /; s/ID=//; s/cuf1=//;' \
# > $cufdir/aphid_rnaseq.cuf1to2.keepids

cat $cufdir/aphid_rnaseq.cuf1to2.ids $lnew | perl -ne\
'if(/^(aphid_cuf8r27\S+)\s(\S+)/) { ($n,$o)=($1,$2); @o=split",",$o; map{ $on{$_}=$n; } @o; } \
else { ($g,$q)=split; if($n=$on{$g}) { s/$g/$n/; s/$/\tUPDATE/; } print; }' > $lnew.2
## missing 43 Gsc old models : fixme

# mv $lnew $lnew.0
/bin/mv $lnew.2 $lnew

# ensure have models for gene mix
#? use this to pull any rnagene: aphid_rnaseq_cufftrinpasavelv.gff

# FIXME: use aphid_rnaseq_mustkeep.gff.old as 1st source of lnew genes,
#        then only grab missing.  But filter out any drops from old.

mv rnagene/aphid_rnaseq_mustkeep.gff rnagene/aphid_rnaseq_mustkeep.gff.old
touch rnagene/aphid_rnaseq_mustkeep.gff

cat  $lnew | grep aphid_cuf8 | grep -v drop | cut -f1  | sed 's/^/=/' | ggrep -F -f - rnagene/aphid_rnaseq.all27cuff8mecdso.an7.gff \
  >> rnagene/aphid_rnaseq_mustkeep.gff

## IDs end w/ p1 == same as .gff
## take from here to get all:  rnas/inch/aphid.trinity.gff.gz
# cat  $lnew | grep aphid2Trinity | grep -v drop | cut -f1 | ggrep -F -f - rnagene/aphid_trinity.ident1.an2.gff 

cat  $lnew | grep aphid2Trinity | grep -v drop | cut -f1 | sed 's/^/=/' > atrids.tmp
env GREP=/usr/sfw/bin/ggrep gzgrep -F -f atrids.tmp rnas/inch/aphid.trinity.gff.gz \
  >> rnagene/aphid_rnaseq_mustkeep.gff
/bin/rm atrids.tmp


cat  $lnew | grep PASAgasmbl_ | grep -v drop | cut -f1 | sed 's/^/=/' | ggrep -F -f - rnagene/pasa2_aphid3.asmbl_bestgenes.an7.gff  \
  >> rnagene/aphid_rnaseq_mustkeep.gff

# ** add velvet genes: aphid_vel (1 altmodel so far?)
cat  $lnew | grep aphid_vel | grep -v drop | cut -f1 | sed 's/^/=/' | ggrep -F -f - rnagene/aphid_vel7asm.ident1.an2.gff  \
  >> rnagene/aphid_rnaseq_mustkeep.gff

## fixme: add ?? best pasa asmbl_ from ../epasa3/pasa_outr/pasa2_aphid3.asmbl_genes.gff.gz
# asmbl_177712    best  < ok, no alternate == PASAgasmbl_177712
# asmbl_14110     best  < no good >> PASAgasmbl_14206 now
# asmbl_31633     altmodel
# asmbl_207199    altmodel
# asmbl_181346    altmodel

# ACYPI29102-RAp1 << need to cut that p1..n from mustkeepdrop.list
cat  $lnew | grep ACYPI | grep -v drop | cut -f1 | sed 's/^/=/'| ggrep -F -f -  genes/acyr1-ACYPImRNA.an7.gff \
  >> rnagene/aphid_rnaseq_mustkeep.gff


mv rnagene/aphid_rnaseq_cufftrin_kinfull.an7.gff rnagene/aphid_rnaseq_cufftrin_kinfull.an7.gff.old

cat rnagene/aphid_rnaseq_cufftrin_infull.an7.gff rnagene/aphid_rnaseq_mustkeep.gff | perl -ne\
'if(/^#/){print;} else { if(/\tmRNA/){ ($g)= m/ID=([^;\s]+)/; $p=($didg{$g}++) ? 0 : 1; } print if $p;} ' \
  > rnagene/aphid_rnaseq_cufftrin_kinfull.an7.gff

mv $lorig $lorig.old
mv $lnew $lorig


################

# look at ncbi genes, overlap but not same CDS as dgil, but w/ identical rnagene models
# .. could in most, all? cases require rnagene models for mustkeep
# .. some of mix7r are pasa-asml w/o CDS ** fix this, add in cds


#   cat acyr2_ncbirefgene.gff | sed 's/^#a.//' | \
#   $evigene/scripts/overgenedup.pl -over bestgenes.DGILmix8c.gff -in \
#   stdin -slopexon=8 -type CDSonly -act drop | \
#   $evigene/scripts/overlapfilter -in stdin -over \
#   bestgenes.DGILmix8c.gff -act keep -pct 50 -pass CDS | grep CDS | \
#   perl -ne'if(m/Parent=([^;\s]+)/){ print "$1\n";}' | sort -u \
#   > acyr2_ncbirefgene.mix8c.cdiff.ids
#   
#   cat acyr2_ncbirefgene.gff | sed 's/^#a.//' | \
#   $evigene/scripts/overgenedup.pl -in bestgenes.DGILmix8c.gff -over \
#   stdin -slopexon=8 -type CDSonly -act drop | \
#   $evigene/scripts/overlapfilter -in stdin -over \
#   acyr2_ncbirefgene.gff -act keep -pct 50 -pass CDS | grep CDS | \
#   perl -ne'if(m/Parent=([^;\s]+)/){ print "$1\n";}' | sort -u \
#   > bestgenes.DGILmix8c.ncbiref.cdiff.ids
  
#   ## fixme, this hits good mix7 genes, want only misses
#   ## fixme2, some of these not in mix7r are .notriv skips that do have expression evidence ; should have kept
#   #   .. bad evd scoring: see aug evd_ flags
#   cat acyr2_ncbirefgene.gff | sed 's/^#a.//' |
#   $evigene/scripts/overgenedup.pl -over bestgenes.DGILmix7r.gff 
#   -in stdin -slopexon=8 -type CDSonly -act drop |
#   $evigene/scripts/overlapfilter -in stdin -over
#   bestgenes.DGILmix7r.gff -act drop -pass CDS | grep CDS | perl
#   -ne'm/Parent=([^;\s]+)/; print "$1\n";' | sort -u >
#   acyr2_ncbirefgene.mix7r.cnone.id1
#   
#   comm -32 acyr2_ncbirefgene.mix7r.cnone.ids acyr2_ncbirefgene.mix7r.cdiff.id1 
#   > acyr2_ncbirefgene.mix7r.cdiff.ids
  
#   ggrep -F -f acyr2_ncbirefgene.mix8c.cdiff.ids \
#   acyr2_ncbirefgene.gff | $evigene/scripts/overgenedup.pl -over \
#   ../rnagene/aphid_rnaseq_cufftrinpasavelv.gff -in stdin \
#   -slopexon=8 -sloputr=149 -act markid -mark rnagene | grep mRNA | \
#   grep rnagene= > acyr2_ncbirefgene.cdiff8.rnagene.mrna
#    
#   ggrep -F -f bestgenes.DGILmix8c.ncbiref.cdiff.ids \
#   bestgenes.DGILmix8c.gff | $evigene/scripts/overgenedup.pl -over \
#   ../rnagene/aphid_rnaseq_cufftrinpasavelv.gff -in stdin \
#   -slopexon=8 -sloputr=149 -act markid -mark rnagene | grep mRNA | \
#   grep rnagene= > bestgenes.DGILmix8c.cdiff8.rnagene.mrna

##   n same (CDS_only):
#     8984 bestgenes.DGILmix7r.ncbiref.csame.ids
#     9567 bestgenes.DGILmix8c.ncbiref.csame.ids
#    11732 acyr2_gnomongene x DGILmix8c
#     6451 acyr1-ACYPImRNA x DGILmix8c
#    30020 DGILmix7r x DGILmix8c (of 36000)
#
##   n diff w/ ident rnagene:
#     1198 acyr2_ncbirefgene.cdiff.rnagene.mrna
#     3534 bestgenes.DGILmix7r.cdiff.rnagene.mrna
#      960 acyr2_ncbirefgene.cdiff8.rnagene.mrna
#     2878 bestgenes.DGILmix8c.cdiff8.rnagene.mrna

#   ggrep -F -f acyr2_ncbirefgene.mix7r.cdiff.ids
#   acyr2_ncbirefgene.gff | $evigene/scripts/overgenedup.pl -over
#   ../rnagene/aphid_rnaseq_cufftrinpasavelv.gff -in stdin
#   -slopexon=8 -pass='mRNA,exon' -typeover=exon -act markid -mark
#   rnaexon | grep mRNA | grep rnaexon= >
#   acyr2_ncbirefgene.cdiff.rnaexon.mrna
#   
#   # ^^ from both rnagene (cds ident), rnaexon (exon ident), pick best rnagene for mustkeep
#    
#   cat acyr2_ncbirefgene.cdiff.rnagene.mrna | perl
#   -ne'm/ID=([^;\s]+)/; print "$1\n";' >
#   acyr2_ncbirefgene.cdiff.rnagene.ids
#   
