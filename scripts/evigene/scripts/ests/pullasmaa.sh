#!/bin/bash
# env gid=Nasvi2EG037195t1 asmid=nvit1v2Svelmid8Loc970t1 asmfa=nvit1_rnaseq.vel2rs13.fa.gz ./pullasmaa.sh

evigene=/bio/bio-grid/mb/evigene

gunzip -c $asmfa | perl -ne'if(/^>(\S+)/) {$t=$1; $p=($t eq $ENV{asmid})?1:0; } print if $p;' > $asmid.fa

$evigene/scripts/genefindcds.pl -act fasta -cdna $asmid.fa | perl -pe's/ / upd=$ENV{gid}; /;' > $asmid.aa
