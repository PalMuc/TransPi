#!/bin/bash
# asmrna2refalign.sh
# see also: evigene/scripts/rnaseq/trasmfixup.pl

## input is gmap output alignments, queryset= asmrna, dbset= refmrna
## gmap --nosplicing -n 9 -S -D $gdb -d $refmrna  $asmrna  > asmname.refname.gmap.out

evigene=/bio/bio-grid/mb/evigene/

# inmap=`/bin/ls {cacao3{tri,sop,cuf13a,v[1345],vel123},cgbAssembly.,rscuff8}*evg3itr.gmap.out8g.gz`
inmap=$*
outgff=gmapall.gff
outtab=align.tab
jointab=join.tab

# stranded "genome" here as ref mRNA.tr are known strands, joins/chimer often map reversed
# nobest=1 means collect all ref aligns/trquery

optgmap2gff="strand=1 noerrspan=1 intron=-1 nobest=1"

for gmapz in $inmap; do {

  pt=`basename $gmapz .gz | sed 's/\.gmap.*//; s/\.fa-/-/; s/\.tr-/-/; s/cgbAssembly./est/;'`
  # echo "$gmapz TO $nam.align.tab"
  
  # gmap_to_gff converts gmap.out to genes.gff with quality scores on mRNA lines
  gzcat $gmapz | env src=$pt $optgmap2gff $evigene/scripts/gmap_to_gff.pl  > $pt.$outgff

# pick ID=id or Target=tid ? ID has extra _Gnn, _Cn for multipath entries, Target same. use that.
  grep mRNA $pt.$outgff | perl -ne \
'BEGIN{ @k=qw(match qlen cov pid path indels cdsindel); 
print join("\t","GeneID","gespan","geor","QueryID","quspan",@k)."\n"; } 
next unless(/^\w/ and /\tmRNA/ and /ID=/ and /Target=/);
@v=split"\t"; ($gid,$gb,$ge,$go,$at)=@v[0,3,4,6,-1]; 
@at= map{ $v=($at=~m/\b$_=([^;\s]+)/)?$1:0; $v; } ("ID",@k);
($tid,$tb,$te)= m/Target=(\S+) (\d+) (\d+)/;  $id= shift @at; 
print join("\t",$gid,"$gb-$ge",$go,$tid,"$tb-$te",@at)."\n"; ' \
 | sort -k1,1 -k2,2n -k4,4 > $pt.$outtab

# calc join scores for multipath entries
cat $pt.$outtab | $evigene/scripts/rnaseq/asmrna2refjoins.pl > $pt.$jointab

# and recombine aligntab, jointab, aaqualtab for summary stats

} done

exit

#.... examples of multipath

# GeneID            gespan  geor  QueryID                     quspan  match   qlen    cov     pid     path    indels  cdsindel
# noj:alt; same locus, alttr
# Thecc1EG000014t1  44-534  +     cacao3sopcsc1k25loc1957t2   1-491   491     902     54.4    100.0   1/2     0       0
# Thecc1EG000014t2  594-879 -     cacao3sopcsc1k25loc1957t2   617-902 286     902     31.7    100.0   2/2     0       0

# noj:triv; trivial rev-td-join : reversed align for parts, but 1st is tiny (0.9% cov, 24 bp)
# Thecc1EG000026t1  4461-4484 -    cacao3sopcsc1k23loc3756t1  2596-2619   24    2619    0.9     100.0   3/3     0       1
# Thecc1EG000027t2  1-2482    +    cacao3sopcsc1k23loc3756t1  61-2542     2479  2619    94.8    99.9    1/3     0       

# noj:multimap; same span multipath, not join
# Thecc1EG045396t1   842-1486 +    cacao3sopcsc10rk23loc108t10  2251-2895 642   2895    22.3    99.5    1/2     0       0
# Thecc1EG045396t2   477-1121 +    cacao3sopcsc10rk23loc108t10  2251-2895 642   2895    22.3    99.5    2/2     0       0

# join:revnearloq;  partial rev-near-join (90t,94t); 85t and 94t are same span multipath; lowish pctid/indels suggests 2ndary aligns or other genes?
# Thecc1EG045390t1   115-667   +   cacao3sopcsc10rk23loc108t12  1-492      460   3323    14.8    92.7    2/3     0/4     -4
# Thecc1EG045394t1   2415-3015 -   cacao3sopcsc10rk23loc108t12  2620-3250  589   3323    19.0    93.3    1/3     30/0    30
# Thecc1EG045385t1   1172-1509 -   cacao3sopcsc10rk23loc108t12  2954-3291  289   3323    10.2    85.5    3/3     0       0


# sums ................
## , only ng= per refgene are useful?
## quality = total valid refgenes (onepath + altof)
# B...... only cdhit parsed prots
# cacao3cuf13aset ng=20163  join=21%,4432   onepath=50%,10124  samespan=3%,657  altof=19%,3980
# cacao3sopc11nr  ng=17881  join=14%,2576   onepath=59%,10662  samespan=4%,832  altof=16%,2939
# cacao3tri1asm   ng=28335  join=18%,5281   onepath=64%,18231  samespan=4%,1359 altof=16%,4805
# cacao3vel12     ng=31170  join=17%,5391   onepath=70%,21887  samespan=5%,1787 altof=22%,7077
##  per asmtr
# cacao3cuf13aset nt=18509  join=14%,2627   onepath=65%,12165  samespan=3%,689  altof=19%,3615
# cacao3sopc11nr  nt=18260  join=7%,1364    onepath=74%,13568  samespan=4%,848  altof=15%,2823
# cacao3tri1asm   nt=50476  join=11%,5821   onepath=69%,35109  samespan=4%,2437 altof=16%,8412
# cacao3vel12     nt=77116  join=8%,6506    onepath=73%,56731  samespan=3%,2705 altof=16%,12404
#
# B..... no prot parsing
# estbean         ng=18044   join=8%,1587    onepath=71%,12982    samespan=2%,435    altof=12%,2312
# estleaf         ng=15621   join=7%,1157    onepath=75%,11738    samespan=1%,273    altof=12%,1921
# cacao3v1asm     ng=36401   join=19%,6924   onepath=74%,27181    samespan=9%,3544   altof=24%,9003
# cacao3v3asm     ng=38532   join=21%,8359   onepath=75%,29229    samespan=11%,4555  altof=27%,10483
# cacao3v4asm     ng=36420   join=18%,6752   onepath=72%,26304    samespan=8%,2951   altof=26%,9624
# cacao3v5asm     ng=35357   join=18%,6647   onepath=71%,25283    samespan=7%,2697   altof=24%,8659
##  per asmtr
# estbean         nt=21204   join=4%,1057    onepath=81%,17309    samespan=2%,531    altof=11%,2536
# estleaf         nt=20982   join=4%,857     onepath=82%,17389    samespan=2%,460    altof=11%,2508
# cacao3v1asm     nt=129160  join=5%,6641    onepath=74%,96092    samespan=4%,6334   altof=17%,21974
# cacao3v3asm     nt=152034  join=4%,6714    onepath=74%,113089   samespan=5%,8219   altof=17%,26476
# cacao3v4asm     nt=128434  join=5%,6558    onepath=73%,94677    samespan=4%,5245   altof=18%,23603
# cacao3v5asm     nt=104424  join=6%,6828    onepath=73%,76789    samespan=3%,4144   altof=17%,18131

