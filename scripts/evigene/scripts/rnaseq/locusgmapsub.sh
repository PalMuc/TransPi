#! /bin/bash
## locusgmapsub.sh  was gmap2locsub.sh 
## rewrite as all perl; subscript of steps per ingmap
## version for trsplit, n=30 parts, input gmap.out, add gmap.gff step

##.......
## Strand/Sense problem, input tr is not stranded mRNA, map-strand, sense flags are ambibuous
## do we flip g-or for -strand to make same locus?  affects equalgene, gmaploci ..
# align.tab
# tscaffold2120 2462791-2490212 +       PitaTv1R022005t84       73-1672 1454    1674    95.6    90.3    0       25/10   2       0       395,70%,complete        0       99      0       3A18VOk40Loc2110t9o13   gmap
# tscaffold2120 2462744-2490249 -       PitaTv1R022005t101      1-1683  1580    1690    99.6    93.2    0       22/13   2       0       395,70%,complete        0       326     -1      6AB25VOk30Loc11414t3o10 gmap
# zGenomeID     gespan  geor    AQueryID        quspan  match   qlen    cov     pid     path    indels  nexon   splice  aalen   offs    aamap   sense   oid     tag
#.......

## .. should be caller env params:
evigene=$HOME/bio/evigene/
pubidtab=../evgtv2x/evgTv1.pubidx
keepat="aalen,oid,offs"

## loctabapp=./gmaploci.pl a caller param? or now use evigene/scripts/rnaseq/locustab.pl
loctabapp=$evigene/scripts/rnaseq/locustab_basic.pl
if [ ! -x $loctabapp ]; then echo "ERR: need locustab app: $loctabapp"; exit -1; fi

# input one gmap.out file name ; or may start with gmap.gff input
gmo=$1

# for gmo in $ingmap; do ... call this subscript
# xxx.gmap.gff input ok here
pt=`basename $gmo .gmap.out13 | sed 's/\.gmap.*//'`
gf=$pt.gmap.gff

if [ ! -f $pt.gmap.gff ]; then
  cat $gmo | env keepat=$keepat src=$pt skiplow=0 nopath=1 best=0 strand=0 noerrspan=1 intron=-1  \
  $evigene/scripts/gmap_to_gff.pl  > $pt.gmapoid.gff0

## add this oid2pubid mapping !!! change Target=oid also? below gff2align fix/bug
  # if [ -f $pubidtab ]; then ..
  cat $pubidtab $pt.gmapoid.gff0 | perl -ne \
'next if(/^\W/); if(/^PitaTv1/){ chomp;($pd,$od,$pg,$ti,$cl,$aq)=split"\t"; $pod{$od}=$pd; $aaq{$od}=$aq; } 
else { chomp; @v=split"\t"; ($td)=m/(?:Parent|ID)=([^\s;]+)/; $od=$td; $spl=""; 
if($od=~s/_([GC]\d+)$//) { $spl=";Split=$1"; }  $pd=$pod{$od}||$od; $aa=$aaq{$od}||0; s/=$td/=$pd$spl/; 
if(/\tmRNA/){s/;/;oid=$od;/; s/Target=$od/Target=$pd/; s/aalen=/aamap=/; s/$/;aalen=$aa/; }  
print $_."\n"; }' \
  > $pt.gmap.gff
  
  gf=$pt.gmap.gff
fi

if [ ! -f $pt.align.tab ]; then
  ## gff2aligntab.pl: $USEPUBID=$ENV{pubid}||0; # versus Target=id
  cat $gf | env pubid=1 $evigene/scripts/ests/gff2aligntab.pl | sed 's/^GenomeID/zGenomeID/;' | \
    sort -k1,1r -k2,2n -k4,4 > $pt.align.tab
fi

#.... equalgene
if [ ! -f $pt.eqgene ]; then
 $evigene/scripts/equalgene.pl -ov $gf -in $gf > $pt.eqgene
fi

## skip  for now: overlapfilter -pass exon -strand -in gf -over gf > $pt.overexon
## should add back?, as check on eqgene loci classing.
 
## add gmaploci.pl $pt.align.tab # makes loci tab, needs aalen not in this gff
## requires $pt.align.tab $pt.eqgene; writes $pt.aa80.loci

$loctabapp $pt.align.tab

##...... CLUSTER script for calling this w/ partitioned genome data inputs .........
#   #! /bin/bash
#   ## env ingmap=xxx.gmap.out datad=`pwd` qsub -q normal gmap2lociset.sh
#   #PBS -N gmap2lociset
#   #PBS -l nodes=1:ppn=32,walltime=9:55:00
#   #PBS -V
#   ## version for trsplit, n=30 parts, input gmap.out, add gmap.gff step
#   
#   ncpu=32
#   ## make these, other params part of env call plus defaults
#   export minaa=89
#   export mincov=10
#   export evigene=$HOME/bio/evigene/
#   
#   if [ "X" = "X$datad" ]; then echo "ERR: missing datad=pathtodata"; exit -1; fi
#   cd $datad/
#   ## works for both .gmap.out, gmap.gff .. skip .out if have .gff
#   if [ "X" = "X$ingmap" ]; then ingmap=`ls *.gmap.gff`; fi
#   if [ "X" = "X$ingmap" ]; then ingmap=`ls *.gmap.out*`; fi
#   
#   i=0; for igm in $ingmap; do { 
#   
#     $evigene/scripts/rnaseq/locusgmapsub.sh  $igm &
#     
#     i=$(( $i + 1 ));
#     if [ $i -ge $ncpu ]; then wait; i=0; fi
#   } done
#   wait
#   
