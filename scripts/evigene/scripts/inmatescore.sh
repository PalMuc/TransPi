#!/bin/tcsh
# inmatescore.sh : intronscore tables for anno v7
# env DOAN=1 workd=$workd genes=(genes/aphid2_{epir2,epir16b,epi4}*.an7.gff .. )

set DOAN=0
set inref=$workd/intron/intron_good.gff.gz
# set inref=intron/intron_all.gff.gz

#.. skip mate scores for wasp, only have for est pairs
set materef=$workd/rnas/rnaseq.matejoin.gff
#set materef=rnas/aphid_rnaseq.all27cuff8.matejoin.gff

# best rna+est assemblies (28k? 10k w/o good homology)
set rnaref=$workd/rnas/rnaseq_best.gff.gz

set protref=$workd/prot/protein_best.gff.gz

set odir=""
if( -d score ) set odir="score/"

foreach geneset ( $genes )

  set gbase=`echo $odir$geneset | sed 's/.gz//; s/.gff//; '`
  echo "make $gbase.{ovintr,ovpro,ovrna} from $geneset"

##new introntab/ovintr from overintron.pl : problem -LONGINTRON makes  nintron=-err often; screws other calcs
#AUGepi6cp1s1g1t1        inqual=66,4/6  << 2nd is replacement for nintron (== validsplice/totalsplice)
#AUGepi6cp1s1g5t1        inqual=76,-3/29/34  << nintron=-error/valid/total

  if( ! -f $gbase.ovintr &&  -f $inref ) then
  $workd/scripts/overintron.pl -LONGINT=99999 -form=score -annotkey=inqual -intron $inref -gene $geneset \
      | perl -pe 's/,/\tnintron=/; s/,inerr:/\tinerr=/;' > $gbase.ovintr
  endif

  #if( ! -f $gbase.ovintr &&  -f $inref ) then
  if( 0 ) then
  $workd/scripts/overlapfilter.perl -intron2splice=error -pass 'exon,intron' -act markid -midtype scoresum \
    -mark intr -over $inref -in $geneset | $workd/scripts/intronscore.pl > $gbase.introntab
    ## cut revise introntab for genescore
  cut -f1-3 $gbase.introntab | egrep -v '^GeneID|^#'| perl -ne 'if(/^\w/){@v=split"\t"; $v[1]="inqual=".$v[1]; \
    $v[2]="intr=".$v[2]; $_=join("\t",@v); print; }' > $gbase.ovintr
  endif
  
  if( ! -f $gbase.ovpro && -f $protref ) then
  $workd/scripts/overgenedup.pl -self -exon=CDS -type CDSsimilar -slopexon=8 -act markid -mark ovpro \
    -over $protref -in $geneset | grep mRNA | perl -ne 'if(/^\w/ and /\tmRNA/) { ($id)=m/ID=([^;\s]+)/;  \
    $ov=(m/(ovpro=[^;\s]+)/)? $1: "na"; if($ov=~m,/[CI]?(\d+),){ $d=$1; $ov=~s/=/=$d,/;} \
    print "$id\t$ov\n";} ' > $gbase.ovpro
  # score format: ovpro=E9J921_SOLIN/C99.77,ID2/88.66,.. best 1st
  # .. modify to put number 1st: ovpro=99.77,E9J921_SOLIN/C99.77,ID2/88.66,..
  endif

  ## replace this inprot with model equiv from overgenedup, better
  # gunzip -c $protref | grep CDS |\
  # $workd/scripts/inmatescore.pl -exontype CDS -over stdin -in $geneset > $gbase.inprot

  if( ! -f $gbase.ovrna && -f $rnaref ) then
  $workd/scripts/overgenedup.pl -self -type CDSsimilar -slopexon=8 -act markid -mark ovrna \
    -over $rnaref -in $geneset | grep mRNA | perl -ne 'if(/^\w/ and /\tmRNA/) { ($id)=m/ID=([^;\s]+)/;  \
    $ov=(m/(ovrna=[^;\s]+)/)? $1: "na"; if($ov=~m,/[CI]?(\d+)[\.]?(\d*),){ $d=$1; $x=$2; $d=$x if($x>$d); $ov=~s/=/=$d,/; } \
    print "$id\t$ov\n";} ' > $gbase.ovrna
  # score  modify to put number 1st: ovrna=99.77,E9J921_SOLIN/C99.77,ID2/88.66,..
  endif

  if( -f $materef ) then
  $workd/scripts/inmatescore.pl -over $materef -in $geneset > $gbase.inmate
  endif

  # - quick score to an1.gff : but use genescoremerge.sh + evigene/annotate_genes instead
  if( $DOAN ) then
      # introntab score col3 = count valid - invalid introns
      cat $gbase.introntab | cut -f1,3 | sed 's,/.*,,; s/^/intr /;' > $gbase.intr.tmp
      
      # inprot change score col: col2=qual bad, use col7=bTrue
      #old# cat $gbase.inprot | cut -f1,7 | sed 's,/.*,,; s/^/ovpro /;' > $gbase.inpro.tmp
      #new: geneid  ovpro=E9J921_SOLIN/C99.77,
      cat $gbase.ovpro | sed 's,ovpro=,,; s/,.*//; s/^/ovpro /;' > $gbase.ovrna.tmp
      
      cat $gbase.ovrna | sed 's,ovrna=,,; s/,.*//; s/^/ovrna /;' > $gbase.ovrna.tmp

      set CAT="cat"
      if( $geneset =~ *.gz ) set CAT="gzcat"

      $CAT $geneset | cat $gbase.{intr,ovrna,ovpro}.tmp - | perl -ne \
      'if(/^(intr|ovrna|ovpro)\s/){ ($k,$d,$v)=split; $dv{$d}{$k}=$v; } \
      else{ if(/\tmRNA/){ ($id)=m/ID=([^;\s]+)/; foreach $k (qw(intr ovrna ovpro)) \
      { if( $v=$dv{$id}{$k}) { s/$/;$k=$v/; } } } print; }' \
      > $gbase.an1.gff

     /bin/rm $gbase.{intr,ovrna,ovpro}.tmp

  endif

end

## for .inmate:
#1  make gff of mated exons (maybe change cuff.exons to per geneset exons?)
#    gzcat velmapt7/aphid_rnaseq.all27cuff8.gff.gz | grep '  exon' | grep '\.1;' | $evigene/scripts/overlapjoins.pl -over stdin -in cleaned/allpe_matec.tab -format bed -strand -sorted > aphid_rnaseq.all27cuff8.matejoin.gff 
#
#2 make gene table of mated.exon scores
#   evigene/scripts/inmatescore.pl -over rnas/aphid_rnaseq.all27cuff8.matejoin.gff -in $geneset > $geneset.inmate
#
#3 add .inmate .introntab scores to genes.gff
#    env geneset=xxxx.gff  addinmate.sh

## old .inprot
# ==> /bio/bio-grid/nasv4/rnas/asmrna/awork3/nvit1_rnaseq.cuff8t13c.pinfix.inprot <==
# ModelID                         Quality         nExon   bpExon  nOvGene nTrue   bTrue   nFNeg   bFNeg   nFPos   bFPos   pFNeg   pFPos
# nvit1cuf83c_Gsc1g92t1           5/partial       1       130     1       1       62      8       665     0       68      0.915   0.523
# nvit1cuf83c_Gsc1g134t1          24/partial      1       183     1       1       183     8       1193    0       0       0.867   0

## new .inprot from overgenedup -exon=CDS -type CDSsimilar -slopexon=8 -act markid -over $prot -in $genes

#... add new intron scoring from overbestgene1 ?? rework overlapfilter, or other to do just new intron splice score

# gzcat $workd/intron/intron_good.gff.gz nvit_epi6c1-augmap.gff.gz | $evigene/scripts/overbestgene1.perl -typeover over -pct 10 -strand -score 'intron:20,CDS:1,UTR:2' -in stdin -skip | grep '   mRNA' > nvit_epi6c1-augmap.ovbest1.mrna &
# cat nvit_epi6c1-augmap.ovbest1.mrna | perl -ne'($d)=m/ID=(\w+)/; @v= grep/=/, m/((ints|cxlen|skip)=[^;\s]+)/g; ($cx)=m/svec=\d+.\d+.(\d+)/; map{ s/,i\d.*$// if(/ints/); s/$/,$cx\%/ if(/cxlen/); } @v; print join("\t",$d,@v),"\n";' 
# AUGepi6cp2s1g3t1        ints=100,6/6    cxlen=207/1170,17%
# AUGepi6cp2s1g10t1       ints=50,2/4     cxlen=603/3424,17%
# AUGepi6cp2s1g11t1       ints=66,4/6     cxlen=249/2707,9%

