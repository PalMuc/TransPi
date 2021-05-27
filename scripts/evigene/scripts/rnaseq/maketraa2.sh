#!/bin/bash
# env trset=xxx.fa.gz  maketraa.sh
# FIXME: use aa flag utrbad/utrpoor here, 
#  1. segregate: ok.aa, poor.aa , 
#  2. cd-hit -i ok.aa -o ok_cd.aa
#  3. cd-hit-2d -i ok_cd.aa -i2 poor.aa -o poor_cd.aa
## fixme2: filter poor.aa also by smallsize
## fixme3: cd-hit -c 0.95 not 0.90 default

## should these be in relation to trsize? ie very bad = pcds<10..20; bad=pcds<20..40, ..
AMINBAD=200
AMINPOO=100

trset=$*
evigene=/bio/bio-grid/mb/evigene
bindir=/bio/bio-grid/mb/bin/

for faz in $trset; do { 
  nam=`basename $faz .gz | sed 's/\.fa//; s/\.tr//;'`
  if test -f $nam.aa ; then  /bin/mv $nam.aa $nam.aa.old; fi

  $evigene/scripts/genefindcds.pl -nostop -mincds=60 -act fasta -cdna $faz > $nam.aa

  perl -pi -e'if(/^>/) { s/\si=\d+;//; s/\s+flag=.+$//; s/\s+path=.+$//; }' $nam.aa

  cat $nam.aa | env mpoo=$AMINPOO mbad=$AMINBAD nam=$nam perl -ne\
'BEGIN{ $mpoo=$ENV{mpoo}; $mbad=$ENV{mbad}; $nam=$ENV{nam}; 
open(OK,">$nam.ok.aa"); open(BAD,">$nam.poor.aa"); }
if(/^>/) { $bad=(/utr(poor|bad)/)?1:0; ($al)=m/aalen=(\d+)/; $skip=0; 
if($bad and $al>0){ $skip=1 if((/utrbad/ and $al<$mbad) or (/utrpoor/ and $al<$mpoo)); } }
if($skip){ } elsif($bad) { print BAD $_; } else { print OK $_; }'

  $bindir/cd-hit -c 0.95 -d 0 -i $nam.ok.aa -o ${nam}.ok_cd.aa >& log.cd1$nam 

  ## should 2d use ok_cd.aa or ok.aa ?
  $bindir/cd-hit-2d -d 0 -i $nam.ok.aa -i2 $nam.poor.aa -o ${nam}.poor_cd.aa >& log.cd2$nam 

} done

exit;

#... count utrqual per kmer.  large kmer have 66% ok, small have 33% ok ?
# 
# grep '^>' daphmag2vel9hallnr.aa | perl -ne\
# '($d)=m/>(\w+)/; ($al,$pc,$cl,$ut)=m/aalen=(\d+).(\d+)..(\w+)[-]?(\w*)/; ($cl)=m/clen=(\d+)/; \
# ($k)=$d=~m/k(\d+)[Ll]oc/; unless($k){ ($k)= $d =~ m/sub(\d+)[Ll]oc/; $k||=0; } $ut||="okutr"; \
# for $k ($k,"all") { $ku{$k}{Tot}++;  $ku{$k}{$ut}++; } $ut{$ut}++; \
# END{ @ut=sort keys %ut; print join("\t","kmer","total",@ut)."\n"; foreach $k (sort keys %ku) { \
# $tc=$ku{$k}{Tot}; print "$k\t$tc"; foreach $u (@ut) { $c=$ku{$k}{$u}; $pc=int(100*$c/$tc);  \
# print "\t$c,$pc"; } print "\n"; } }' 
#
#......... cacao soap  cacao3sopc11nr.aa 
# kmer    total   okutr        utrbad      utrpoor
# 23      61739   18437,29     30883,50    12419,20
# 25      67397   20304,30     33637,49    13456,19
# 29      63710   19532,30     31329,49    12849,20
# 35      64010   20936,32     29960,46    13114,20
# 49      63563   25038,39     25133,39    13392,21
# 69      2608    1387,53      715,27      506,19
# 89      27224   17895,65     4657,17     4672,17
# all     350251  123529,35    156314,44   70408,20
#
#........ cacao velvet velw/vel4g/cacao3vel4g3nr.aa (kset matching soap..)
# kmer    total   okutr        utrbad     utrpoor
# 21      151729  48985,32     76353,50   26391,17
# 25      122015  29023,23     72118,59   20874,17
# 29      109005  25059,22     64907,59   19039,17
# 35      94496   21656,22     55701,58   17139,18
# 47      71742   19918,27     38170,53   13654,19
# 91      19517   11786,60     4629,23    3102,15
# all     568504  156427,27    311878,54  100199,17
#...... all vel4g.kmer
# kmer    total   okutr        utrbad      utrpoor
# 21      151729  48985,32     76353,50    26391,17
# 23      1482    354,23       902,60      226,15
# 25      122015  29023,23     72118,59    20874,17
# 29      109005  25059,22     64907,59    19039,17
# 35      94496   21656,22     55701,58    17139,18
# 41      83193   21005,25     46629,56    15559,18
# 47      71742   19918,27     38170,53    13654,19
# 49      663     144,21       455,68      64,9
# 51      64680   19545,30     32354,50    12781,19
# 55      677     148,21       458,67      71,10
# 69      633     161,25       405,63      67,10
# 79      562     155,27       352,62      55,9
# 89      407     107,26       245,60      55,13
# 91      19517   11786,60     4629,23     3102,15
# 93      353     96,27        210,59      47,13
# all     721154  198142,27    393888,54   129124,17
#
#...... cacao trinity (no kmers) cacao3tri1asmnr.aa 
# kmer    total   okutr        utrbad     utrpoor
# all     95606   29644,31     46110,48   19852,20
#
#.....  cacao cufflinks13 (no kmer) cacao3cuf13aset.aa
# kmer    total   okutr        utrbad     utrpoor
# all     39244   15339,39     16350,41   7555,19
#
# ..........................................................
#  12/oct/23  .. cacao redo new vel12; maketraa2.sh w/ cdhit filt
#
#.......  vel4g/cacao3vel4g3nr.{ok,poor}_cd.aa ; drop few kmer: 93,89..55,49,23
# kmer    total   okutr           utrbad          utrpoor
# 21      29650   22835,77%,141aa 3133,10%,294aa  3682,12%,225aa
# 25      19564   10868,55%,187aa 4695,23%,308aa  4001,20%,249aa
# 29      17499   8117,46%,202aa  5321,30%,319aa  4061,23%,265aa
# 35      15749   6499,41%,248aa  5310,33%,316aa  3940,25%,268aa
# 41      14631   5670,38%,247aa  5031,34%,327aa  3930,26%,272aa
# 47      13502   5183,38%,288aa  4439,32%,320aa  3880,28%,284aa
# 51      13113   5078,38%,270aa  3904,29%,326aa  4131,31%,288aa
# 91      6188    4580,74%,299aa  352,5%,326aa    1256,20%,236aa
# all     130878  69451,53%,205aa 32418,24%,317aa 29009,22%,263aa
# kmer    total   complete        partial         partial3        partial5
# 21      29650   13105,44%,261aa 6239,21%,97aa   2172,7%,109aa   8134,27%,88aa
# 35      15749   11478,72%,329aa 1276,8%,116aa   590,3%,157aa    2405,15%,138aa
# 41      14631   10885,74%,331aa 1115,7%,114aa   496,3%,176aa    2135,14%,139aa
# 51      13113   9671,73%,339aa  967,7%,137aa    457,3%,195aa    2018,15%,166aa
# 91      6188    3437,55%,343aa  887,14%,220aa   509,8%,250aa    1355,21%,207aa
# all     130878  82750,63%,317aa 16049,12%,113aa 6647,5%,152aa   25432,19%,122aa

#....... cacao3/rnas/velw/vel12/cacao3vel12h.nrtr.{ok,poor}_cd.aa
#.. vel12: No 454est, only proper paired reads mapped to chrs (same data as soapc11?/trin)
#.. has improved %okutr over vel4g, but slightly shorter aa ave. * more partialaa *
# kmer    total   okutr           utrbad          utrpoor
# 21      19633   16714,85%,135aa 986,5%,283aa    1933,9%,184aa
# 25      28138   21269,75%,164aa 2980,10%,300aa  3889,13%,218aa
# 35      23526   14722,62%,238aa 4699,19%,314aa  4105,17%,251aa
# 45      18665   11537,61%,235aa 3449,18%,321aa  3679,19%,252aa
# 51      16952   10798,63%,220aa 2809,16%,323aa  3345,19%,251aa
# 91      7724    6373,82%,287aa  178,2%,321aa    1173,15%,211aa
# all     114638  81413,71%,199aa 15101,13%,313aa 18124,15%,234aa
# kmer    total   complete        partial         partial3        partial5
# 21      19633   6918,35%,228aa  4818,24%,110aa  1688,8%,133aa   6209,31%,90aa
# 25      28138   12355,43%,280aa 5812,20%,112aa  2251,7%,145aa   7720,27%,105aa
# 35      23526   13877,58%,335aa 3184,13%,129aa  1477,6%,185aa   4988,21%,138aa
# 45      18665   11056,59%,329aa 2346,12%,133aa  1264,6%,205aa   3999,21%,132aa
# 51      16952   9384,55%,324aa  2340,13%,134aa  1288,7%,203aa   3940,23%,128aa
# 91      7724    4691,60%,327aa  1154,14%,191aa  490,6%,238aa    1389,17%,188aa
# all     114638  58281,50%,307aa 19654,17%,124aa 8458,7%,173aa   28245,24%,119aa

#....  cacao3/rnas/trin/cacao3tri1asmnr.{ok,poor}_cd.aa.gz
# kmer    total   okutr           utrbad          utrpoor
# all     38077   25229,66%,272aa 5910,15%,377aa  6938,18%,318aa
# kmer    total   complete        partial         partial3        partial5
# all     38077   24450,64%,398aa 4956,13%,117aa  2141,5%,151aa   6530,17%,101aa

#....  cacao3/rnas/cuff/cuff13ars39/cacao3cuf13aset.{ok,poor}_cd.aa
# kmer    total   okutr   utrbad  utrpoor
# all     20942   12641,60%,391aa 3553,16%,389aa  4748,22%,307aa
# kmer    total   complete        partial         partial3        partial5
# all     20942   15790,75%,412aa 1368,6%,199aa   1811,8%,362aa   1973,9%,180aa

#..... cacao3/rnas/trsoap/cacao3sopc11nr.{ok,poor}_cd.aa.gz
# kmer    total   okutr   utrbad  utrpoor
# 23      17476   12861,73%,302aa 2116,12%,365aa  2499,14%,292aa
# 25      12670   7580,59%,239aa  2328,18%,367aa  2762,21%,290aa
# 29      11445   6626,57%,213aa  2244,19%,373aa  2575,22%,289aa
# 35      12556   7575,60%,209aa  2147,17%,388aa  2834,22%,294aa
# 49      16260   11466,70%,185aa 1705,10%,356aa  3089,18%,287aa
# 69      666     549,82%,159aa   18,2%,294aa     99,14%,227aa
# 89      8259    6053,73%,249aa  202,2%,312aa    2004,24%,214aa
# all     79332   52710,66%,235aa 10760,13%,369aa 15862,19%,280aa
# kmer    total   complete        partial         partial3        partial5
# 23      17476   11876,67%,404aa 1730,9%,104aa   782,4%,138aa    3088,17%,99aa
# 25      12670   8308,65%,359aa  1467,11%,105aa  579,4%,141aa    2316,18%,107aa
# 29      11445   7401,64%,345aa  1290,11%,106aa  556,4%,144aa    2198,19%,101aa
# 35      12556   7841,62%,351aa  1504,11%,99aa   582,4%,137aa    2629,20%,100aa
# 49      16260   8819,54%,333aa  2432,14%,93aa   839,5%,112aa    4170,25%,86aa
# 69      666     306,45%,242aa   140,21%,119aa   51,7%,114aa     169,25%,112aa
# 89      8259    5459,66%,286aa  941,11%,144aa   397,4%,168aa    1462,17%,161aa
# all     79332   50010,63%,353aa 9504,11%,105aa  3786,4%,136aa   16032,20%,103aa

#==============================
#........ daphnia magna soap7  daphmag3sopcsub1to16_nr.aa 
# kmer    okutr       utrbad        utrpoor
# 25      22455,51    13470,30      7917,18
# 29      22363,51    13085,30      7792,18
# 35      22077,53    11979,28      7261,17
# 45      25586,61    9634,23       6404,15
# 55      29545,69    7532,17       5651,13
# 65      27430,71    5938,15 	    4798,12
# 
#........ daphmag velvet9 daphmag2vel9hallnr.aa 
# kmer    total   okutr      utrbad     utrpoor
# 29      105916  64282,60   25583,24   16051,15
# 39      129029  83866,64   27817,21   17346,13
# 49      54811   30920,56   15023,27   8868,16
# 59      79913   53493,66   15633,19   10787,13
# 69      50685   35371,69   9242,18    6072,11
# all     420354  267932,63  93298,22   59124,14
#
#..... daphmag trinity7, no kmers; by subset, left out 10..16
# kmer    total   okutr   utrbad  utrpoor
# 1       7211    3650,50 2321,32 1240,17
# 2       9243    3953,42 3583,38 1707,18
# 3       7241    3313,45 2656,36 1272,17
# 4       9552    4565,47 3321,34 1666,17
# 5       6553    3086,47 2333,35 1134,17
# 6       6582    3297,50 2193,33 1092,16
# 7       7386    3456,46 2593,35 1337,18
# 8       5571    2627,47 1887,33 1057,18
# 9       8854    4403,49 2832,31 1619,18
# all   463099  273775,59  115254,24  74070,15  # includes nomap sub0 w/ mouse genes
#
#..... daphmag cufflinkth7, no kmers daphmag3cuf13th3.aa 
# kmer    total   okutr      utrbad     utrpoor
# all     28499   10700,37   13530,47   4269,14
#
#.... daphmag genes2001, no kmers: daphmagna_201104m8.aa 
# kmer    total   okutr      utrbad     utrpoor
# all     33484   18334,54   10354,30   4796,14
#
