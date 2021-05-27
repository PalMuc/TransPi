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
