#!/bin/bash
# evaltrx.sh : evaluate tr-assembly using est long reads 
# blastn map reads x trdb; count per read: nmap, avematch, bitscore; one quick/reliable eval. of trasm quality
# need similar quick blastp bestprot x aadb as est-eval isnt enough;
# ** add simpler eval CDS/exon ratio for quality : modify genefindcds for this? or just use header info.
# .. use this est-blastn, prot-blastp and CDS/exon stats to pick best of transcripts (for merge and for final)


# allow stdin for trs and/or est?  .gz **?
nbin=/bio/bio-grid/mb/ncbicpp/bin
dbdir=trdb
trs=$1
rds=$2
flags=$3

if ! test -f $trs ; then echo "missing trs; use: evaltrx.sh trseqs readseqs"; exit; fi
if ! test -f $rds ; then echo "missing reads; use: evaltrx.sh trseqs readseqs"; exit; fi
if ! test -x $nbin/makeblastdb ; then echo "please edit ncbibin path $nbin"; exit; fi


dbnam=`basename $trs .tr | sed 's/^cacao3//; s/.gz//; s/.fa//; s/.tr//;  s/.fasta//;'`
rdnam=`basename $rds .fa | sed 's/^reads.//; s/.gz//; s/.fa//; s/.tr//; s/.fasta//;'`
blout=eval.$rdnam-$dbnam.mblout

# if [ test -f $dbdir/$dbnam.nsq && ! $flag =~ /clean|keep/ ]; then
#   echo "have $dbdir/dbnam; use: evaltrx.sh trseqs readseqs flag=keep|clean"; exit
# fi  

if ! test -d $dbdir ; then mkdir $dbdir; fi
# need opt to make clean
if ! test -f $dbdir/$dbnam.nsq ; then
TCAT=cat; tz=`echo $trs | sed 's/.gz//'`; if [ $trs != $tz ]; then TCAT=gunzip -c; fi
$TCAT $trs | $nbin/makeblastdb -title $dbnam -out $dbdir/$dbnam -dbtype nucl -logfile $dbdir/log.m$dbnam
fi

if ! test -f $blout ; then
RCAT=cat; tz=`echo $rds | sed 's/.gz//'`; if [ $rds != $tz ]; then RCAT=gunzip -c; fi
$RCAT $rds | $nbin/blastn -evalue 1e-5 -db $dbdir/$dbnam -outfmt 7 -out $blout
fi

# add sort -k1,1 -k12,12nr to pick top match if have several blasts
echo -n "eval.$rdnam-$dbnam : " ; cat $blout | ggrep -A4 '^# Query' | egrep '^[a-zA-Z]' | perl -ne \
'next unless(/^\w/); ($q,$t,$pi,$ma,$mi,$x,$qb,$qe,$tb,$te,$ev,$bits)=split; next if($did{$q}++); \
$n++; $sm+=$ma; $spi+=$pi; $sbit+=$bits; \
END{ $abit=int($sbit/$n); $am=int($sm/$n); $ap=int(10*$spi/$n)/10; \
print "n=$n; match=$am; pident=$ap; bits=$abit; matchsum=$sm\n"; }'

#..... VelvetO retests  w/ various data slices, parameters ....
# .. problem of lower velo qual appears related to incorporation of read-mistakes/frameshifts/...
# .. but cannot so far trace to data causing this, and am using same data as before (or new slices)
# *add stats like gmap: pmapped=nmap/nest, pcover=ave match/estlen, mismat=ave mi
# *** CLEAN READS *** velv needs input of cleaned reads (remove lowqual, using sickle..)

# reads.tsh1188.s2r2-all: nt=599600  bp=383737523 alen=639;  573737 (96%) mapped to mars11asm;
# reads.tsh1188.s2r2-sc6: nt=45582   bp=29439521  alen=645;
#
# eval.tri1sc6    : n=45144; match=633; pident=98.8; bits=1126; matchsum=28613574
# eval.tri4sc6    : n=45121; match=634; pident=98.8; bits=1129; matchsum=28647559  # trimmed reads, same as vel4g3sc6
# eval.cuff08sc6  : n=44610; match=629; pident=98.8; bits=1121; matchsum=28091564
# eval.cuff13sc6  : n=44090; match=599; pident=98.8; bits=1068; matchsum=26448052
# eval.vel43sc6   : n=45098; match=610; pident=97.8; bits=1047; matchsum=27512025 ; best kmers: 12251 sc6k25, 32847 sc6k35
#... ^^ best of vel so far, near prior vel4.. vel43sc6  ** use higher kmers v44 kset="55 45 39 29"
# eval.vel44sc6   : n=45135; match=616; pident=97.5; bits=1048; matchsum=27808200 ; bestkm: 3417 44sc6k29, 4995 44sc6k39, 9115 44sc6k45, 27608 44sc6k55
# eval.velm4sc6a  : n=45109; match=638; pident=97.5; bits=1083; matchsum=28809281 ; merge vel4. works? not quite velall
#
# eval.vel45sc6   : n=45161; match=628; pident=98.4; bits=1101; matchsum=28400740 << best single run
# eval.vel46sc6   : n=45179; match=616; pident=97.8; bits=1058; matchsum=27846279
# eval.velm4sc6   : n=45147; match=641; pident=97.5; bits=1085; matchsum=28964145 # ok but not as good as should be
# eval.vel4allsc6 : n=45224; match=639; pident=98.5; bits=1124; matchsum=28910171 # ^all above part asm; ave.best but how to extract each best asm?
#  vel4all best kmer set for EST : do same for aaqual
#  10590 46sc6k89, 5174 46sc6k49, 6603 43sc6k35, 4575 45sc6k81, 3195 45sc6k71
#   3158 45sc6k91, 2080 43sc6k25, 1320 44sc6k55, 1002 46sc6k29
# vel47all: estbest kset="93 89 79 69 55 49 35 29 25 23"
#  most long aa from k49, but k23,29,35 are similarly high
#
#.. Cleaned reads; vel4t3=rmdup, t2=keepdup; also aa are longer.
# eval.vel4t3sc6  : n=44985; match=630; pident=98.5; bits=1108; matchsum=28351746   # in progress 5 kmers 89..41
# eval.vel4t3sc6  : n=45139; match=636; pident=98.5; bits=1117; matchsum=28711137   # finished 8 kmers: 89..25
# eval.vel4t2sc6  : n=45139; match=630; pident=98.3; bits=1101; matchsum=28480767   # in progress
#
# .. from gordon.sdsc, as vel4t2 but slighly diff kmer set (k89>k91, drop k69/add k21<poor choice?)
# eval.vel4g3sc6  : n=45156; match=635; pident=98.3; bits=1111; matchsum=28714978   # not quite above vel4t3sc6 .. more kmers is better for this test
# eval.velm4g3sc6 : n=45090; match=637; pident=97.9; bits=1097; matchsum=28758573   # merge, not quite as good, or for aa
#
# vel4t3 kmers: 1287 3sc6k25, 1848 3sc6k29, 2212 3sc6k35, 2799 3sc6k41, 3640 3sc6k47,
#     7875 3sc6k51, 12374 3sc6k69, 13104 3sc6k89
# vel4g3 kmers: 988 3sc6k21, 1152 3sc6k25, 1816 3sc6k29, 2822 3sc6k35, 3673 3sc6k41,
#     5648 3sc6k47, 15379 3sc6k51, -k69.., 13678 3sc6k91

# cat cdhits/cacao3vel47sc*_cd.aa.clstr | perl -ne'if(/^>Cluster (\d+)/){ $cl=$1; print $al{$ma}  if($lc
# and $al{$ma}); $ma=0; %al=(); next; } s/>//; s/aa,//; s/\.\.\.//; s/ at//; ($i,$al,$tr,$pc)=split; ($km)=$tr=~m/(
# k\d+)Loc/; next if($al<200); ($g=$tr)=~s/t\d+//;  if(/ \*/){ $ma=$al; print "$km\n"; $did{$g}++; } else { $al{$al}
# .="$km\n" unless($did{$g}++); } $lc=$cl; ' | sort | uniq -c | sort -k1,1nr
#
# 5281 k49, 5133 k35, 4994 k23, 4823 k29, 4674 k55, 4479 k25, 4157 k69, 3506 k79, 2501 k89, 2169 k93

#..
# v45 bestkm: 1350 45sc6k101, 2744 45sc6k31, 3081 45sc6k41, 4384 45sc6k51, 
#             5265 45sc6k61, 6548 45sc6k71, 8374 45sc6k81, 13415 45sc6k91**
# v46 bestkm: 1147 46sc6k21, 1497 46sc6k23, 2211 46sc6k25, 
#             5124 46sc6k29, 17188 46sc6k49*, 18012 46sc6k89**


# eval.vel94sc6   : n=44897; match=599; pident=97.7; bits=1031; matchsum=26937582 ; bestkm: 12486 94sc6k45, 32411 94sc6k63
# eval.velm9sc6   : n=45074; match=621; pident=97; bits=1036; matchsum=28026300  ; merge of vel9s; not so good

# .. best merge all cur vel : vel4/cacao3vel4[34]sc6.tr vel6/cacao3vel6sc6.tr vel8/cacao3vel8sc6.tr vel9/cacao3vel92sc6.tr
# eval.velasc6    : n=45206; match=637; pident=98; bits=1100; matchsum=28835014  << getting near best
# .. best parts: 4559 43sc6k25, 13753 43sc6k35, 2505 44sc6k45, 10821 44sc6k55
#                 295 6sc6k23,  515 6sc6k27, 2074 6sc6k39, 587 8sc6k23, 2123 8sc6k35
#                 792 92sc6k25, 5811 92sc6k35
# .. try highest kmer = 105; for cgb 106bp reads??

# 2kmer-eval.vel44sc6   : n=44828; match=595; pident=97.7; bits=1019; matchsum=26717086   # 13189 sc6k45 31639 sc6k55
#
# eval.cuff13sc6 : n=44090; match=599; pident=98.8; bits=1068; matchsum=26448052
#
#... vv why are newer velv runs so poor? and getting worse each new trial; vel6 has wider kmers, try highest kmer?
# eval.vel6sc6   : n=45075; match=576; pident=96.4; bits=943; matchsum=25967756
# eval.vel7sc6   : n=44999; match=551; pident=95.9; bits=891; matchsum=24798477
# eval.vel8sc6   : n=44989; match=564; pident=95.9; bits=910; matchsum=25383817
# eval.vel92sc6  : n=45031; match=556; pident=96.4; bits=912; matchsum=25081003
# eval.vel93sc6  : n=44747; match=548; pident=96.5; bits=907; matchsum=24546296 # changed opts not helpful v92
#.. using velvet10 old version, old data ; not same soft as vel4.0?
# eval.vel42sc6  : n=44753; match=485; pident=95.7; bits=775; matchsum=21738593


