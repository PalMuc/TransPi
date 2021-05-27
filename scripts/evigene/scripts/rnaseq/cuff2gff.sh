#!/bin/tcsh

set evigene=/bio/bio-grid/mb/evigene
#set workd=/bio/bio-grid/cacao3
#set dgenome=$workd/genome/cacao11allasm.fa
#set fpre=cacao3
set workd=/bio/bio-grid/kfish2
set dgenome=$workd/genome/kfish2asm.fa
set fpre=kfish2
set MAKETR=1
set MIN_TR=180

set cufdir=$1
set src=$2
if ( -d $cufdir ) then

else
 echo "usage: $0  cufdir sourcename "; exit
endif

set cuffs=`find $cufdir/ -name transcripts.gtf  \! -size 0`

foreach ct ( $cuffs )
 set nam=`echo $ct | sed -e"s,/transcripts.gtf,,; s,$cufdir/,,; s,subset.*/,,;  s/subset.//; s/.cuff//; s,/,-,g; " `
  set trnam=$fpre$nam

#  #test name
 env nam=$nam perl -e\
'BEGIN{$in=$_=$ENV{nam}; s/.sam//;s/parts./p/; s/Scaffold/sc/i; s/combo_/c/i; s/cuff.//i; s/\W/_/g; $na=lc($_);}\
 print "in=$in name=$na\n"; exit; '
# continue;

## fix funky chr names for cacao11:  cacao1mito, cacao1chloroplast, ..
## transcript_id "CUFF.5.1";
## drop tiny mRNA / trsize < 40?; change name to Scaffold scNNN uniq

  cat $ct | env src=$src nam=$nam perl -ne\
'BEGIN{ $src=$ENV{src}; $in=$_=$ENV{nam}; s/.sam//; s/parts./p/; s/Scaffold/sc/i; s/combo_/c/i; s/cuff.//i; \
s/\W/_/g; $na=lc($_); $src ||="G$na"; print "##gff-version 3\n"; } \
($r,$b,$e,$v)=(split)[0,3,4,5]; ($sn=$r) =~ s/(Scaffold|super)[_]?/sc/i; $sn=~s/contig[_]?/ct/i; \
$sn=~s/mito/mt/; $sn=~s/chloroplast/cpl/; if($na=~m/(p\d+)/){ $sn.=$1; }\
s/\ttranscript/\tmRNA/; s/\tCufflinks/\t$src/; s/"//g; ($tid)=m/transcript_id ([^\s;]+)/; \
($gid=$tid) =~ s/\.(\d+)\.(\d+)$/g$1t$2/; $gid=~s/CUFF/${src}_G$sn/; \
$dt=(/\tmRNA/)?"ID":"Parent";  s/transcript_id $tid/$dt=$gid/; \
s/ frac .*$//; s/gene_id \S+ //;  \
s/ exon_number /xi=/; s/ FPKM (\S+)//; $k=$1; $k=($k>999)?int($k):sprintf"%.3g",$k;  \
s/\t$e\t$v/\t$e\t$k/; $w=1+$e-$b; if(/mRNA/){$skip=($w<40)?1:0;}\
print unless($skip); $lid=$tid; '\
  > $trnam.gff


# cat rscuff-*.gff > allcuff-rnaseq.gff
if( $MAKETR != 0 ) then
  $evigene/scripts/cdsgff2genbank.pl -MIN_CDS=$MIN_TR -v -pretty -t=exon -a dna \
    -gff $trnam.gff -fasta $dgenome > $trnam.tr

## add -nostop for cd-hit ??
##otherscript#  $evigene/scripts/genefindcds.pl -mincds=60 -act fasta -cdna $trnam.tr  > $trnam.aa 
 
endif

end

