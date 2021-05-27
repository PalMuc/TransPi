#!/bin/tcsh
# evalpred.sh

## FIXME: here or in annot, -strand should be used for EST, prot evidence check
## gene2/  .an6, new predictors, drop old _arb, epit1?
## add Prot Homology scores, pctfound + bitscore, from mRNA ho3= and harabid= 
## add transposon score 
## add overlapfilt -strand flag and redo all

set td=scripts

set predicts=`/bin/ls genes/dmag-*.an1.gff.gz`
# include dmag_ep24augmap2.gff.gz

## add intron/intron eval ? count splice matches

#..
set evds=(est/est_uniq prot/protein_uniq rnas/rnaseq_uniq misc/transposon)
set evpass=('exon,HSP' CDS 'exon', 'exon,transposon')

set ei=0
#if( 1 == 0 ) then 

foreach evd ($evds)
@ ei = $ei + 1
set passtype=$evpass[$ei]
echo "Evidence: $evd"
foreach pred ($predicts)
echo "Prediction: $pred"
$td/overlapfilter -strand -pass "$passtype" -in $evd.gff.gz -over $pred -act keep -pct 50 -base > /dev/null
echo "-----------------------------"
echo
end
echo
end

#endif

## not ready for -in intron way; -over introns only works now
# $td/overlapfilter -intron2splice -pass 'exon,intron' -act keep -pct 50 -base \
#   -over intron/intron.gff.gz -in $pred  > /dev/null

#if( 1 == 0 ) then  

## all_evd_exons = uniq of est + rnaseq exons + prot cds
## is -sumbase right or wrong here; got bad nums with it
echo "Evidence: all_evd_specif"
foreach pred ($predicts)
echo "Prediction: $pred"
$td/overlapfilter -strand -pass "exon" -over est/all_evd_exons.gff.gz -in $pred -act keep -pct 50 -base > /dev/null
echo "-----------------------------"
echo
end

#endif

## add full cDNA gene accuracy test
## ?? add refseq/refseq-genes.gff here ?
if( 1 == 0 ) then  

echo "Evidence: cDNA_gene_accuracy"
foreach pred ($predicts)
echo "Prediction: $pred"
$td/overlapeval.pl -strand -pass exon -pct 50 -in epasa/pasatrain_genes.best1.gff.gz -over $pred
echo "-----------------------------"
echo
end

endif

## add full protein gene accuracy: prot/protein_uniq.gff.gz
#if( 1 == 0 ) then

echo "Evidence: Prot_gene_accuracy"
foreach pred ($predicts)
echo "Prediction: $pred"
$td/overlapeval.pl -strand -pass CDS -pct 50 -in prot/protein_uniq.gff.gz -over $pred
echo "-----------------------------"
echo
end

#endif


##  Prot Homology scores, pctfound + bitscore, from mRNA ho3= and harabid= 
if( 1 == 0 ) then  

echo "Evidence: Protein_Homology"
foreach pred ($predicts)
echo "Prediction: $pred"
gzgrep '	mRNA' $pred | perl -ne\
'($aa)=m/aalen=(\d+)/; $saa+=$aa; $ng++; \
 ($ho)=m/;ho3=([\d\.]+)/; $sho+=($ho/$aa) if $ho; $no++ if $ho; \
 ($ha)=m/;horef=([\d\.]+)/; $sha+=($ha/$aa) if $ha; $na++ if $ha;  \
 END{ $na||=1; $no||=1; printf "protein_homol best n=%d, bits/aa=%.3f, arabid n=%d, bits/aa=%.3f \n", \
  $no, $sho/$no, $na, $sha/$na; }' 

echo "-----------------------------"
echo
end

endif

#if( 1 == 0 ) then

## fixme: remove #commented features ..
echo "Evidence: Gene_coverage"
foreach pred ($predicts)
echo "Prediction: $pred"
echo -n 'Coding bases: '
gzgrep -v '^#' $pred | grep '	CDS' | perl -ne\
 '($r,$b,$e)=(split)[0,3,4]; $sw+= 1+$e-$b; $n++; END{ $aw=int(10*$sw/(1024*1024))/10; print "$aw\n"; }'
echo -n 'Exon bases: '
gzgrep -v '^#' $pred | grep '	exon' | perl -ne\
 '($r,$b,$e)=(split)[0,3,4]; $sw+= 1+$e-$b; $n++; END{ $aw=int(10*$sw/(1024*1024))/10; print "$aw\n"; }'
echo -n 'Gene_count: '
gzgrep -v '^#' $pred | grep -c '	mRNA'
echo "-----------------------------"
echo
end

#endif

