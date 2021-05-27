#!/bin/bash
# env aa=xxx.aa blastp=xxx.blastp.gz omcl=$carp ./maketop.sh

evigene=/bio/bio-grid/mb/evigene
skipho=Thecc
skips=cacao

onam=`basename $blastp .blastp.gz | sed 's/.aa//; s/.blastp//;'`

env skipho=$skipho aa=$aa tall=1 $evigene/scripts/makeblastscore.pl $blastp > $onam.tall4 

$evigene/scripts/eval_orthogroup_genesets.pl -bitmed -mintaxa 3 -skips $skips \
 -out $onam.topout2 -tallscore $onam.tall4 \
 -groupgene $omcl/*_omclgn.tab -groupcount $omcl/*-orthomcl-count.tab -geneaa $omcl/*.aa.gz

