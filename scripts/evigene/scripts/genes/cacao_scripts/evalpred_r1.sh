#!/bin/tcsh
# evalpred.sh

## FIXME: here or in annot, -strand should be used for EST, prot evidence check

set td=scripts

## up to .an3 now
## .. drop out cacao9_epit1 as useless
set OLDpredicts=(genes/cacao9_arb-augmap.an3 genes/cacao9_cacao0-augmap.an3  \
genes/cacao9_epit2-augmap.an3 \
genes/cacao9_epi3-augmap.an3 genes/cacao9_epi4-augmap.an3 \
genes/cacao9_mix1 genes/cacao9_mix2 )

#set predicts=( genes/cacao9_mix5 )
#set predicts=( genes/cacao9_mix6 )
set predicts=( genes/cacao9_mix7 genes/cacao9_fgenesh_fix.an4 )

## add intron/intron eval ? count splice matches
set evds=(est/est_uniq prot/protein_uniq rnas/rnaseq)
set evpass=('exon,HSP' CDS 'exon')

set ei=0
# if( 1 == 0 ) then 
foreach evd ($evds)
@ ei = $ei + 1
set passtype=$evpass[$ei]
echo "Evidence: $evd"

foreach pred ($predicts)
echo "Prediction: $pred"
$td/overlapfilter -pass "$passtype" -in $evd.gff.gz -over $pred.gff.gz -act keep -pct 50 -base > /dev/null
echo "-----------------------------"
echo
end
echo
end
# endif

# if( 1 == 0 ) then  
## all_evd_exons = est + rnaseq exons
echo "Evidence: all_evd_specif"
foreach pred ($predicts)
echo "Prediction: $pred"
$td/overlapfilter -pass "exon" -sumbase -over est/all_evd_exons.gff.gz -in $pred.gff.gz -act keep -pct 50 -base > /dev/null
echo "-----------------------------"
echo
end
# endif

## add full cDNA gene accuracy test
#if( 1 == 0 ) then  

echo "Evidence: cDNA_gene_accuracy"
foreach pred ($predicts)
echo "Prediction: $pred"

$td/overlapeval.pl -strand -pass exon -pct 50 \
 -in epasa3/pasatrain_genes.best1.gff.gz -over $pred.gff.gz

echo "-----------------------------"
echo

end
#endif

