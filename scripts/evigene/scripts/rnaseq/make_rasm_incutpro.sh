#!/bin/tcsh
# env workd=$workd genes=xxxx.gff  makeincutpro.sh

# ** obsolete, instead use
# $evigene/scripts/genefindcds.pl -dna $workd/genome/nasvit1asm.fa -intron $workd/intron/intron_good.gff.gz \
#  -genes $gset.gff.gz > $gset.pinfix.gff

set workd=/bio/bio-grid/nasv4
setenv PERL5LIB ${PERL5LIB}:${workd}/PerlLib

## dammmmmmmmmm
set path=(/usr/local/bin /bio/bio-grid/mb/bin /usr/sfw/bin /usr/bin  /usr/sbin  /usr/java/bin /usr/ccs/bin )

# env gset=nvit1_rnaseq.vel2rs13
set gset=`echo $genes | sed 's/.gz//; s/.gff//;' `

$workd/scripts/overlapfilter -in $genes  -over $workd/intron/intron_good.gff.gz \
 -pass 'exon,intron' -strand -pct 100 -typeover cut -act mark -mark inov \
 | grep -v CDS > $gset.incut.tmp
    
## bad env path; miss cdbyank
($workd/scripts/pa2dgg_gff3proteins.pl $gset.incut.tmp $workd/genome/nasvit1asm.fa addprotpart \
  > $gset.incutcd.tmp ) > & log.padd$gset

cat $gset.incutcd.tmp | perl -ne\
'if(/\tgene/){ next;} if(/\tmRNA/) { putg(); $dt="ID"; $u5=$u3=0; } \
elsif(/_prime_utr\t/) { $u5++ if(/five/); $u3++ if(/three/); next; } \
else { $dt="Parent"; s/ID=[^;\s]+(\d+);//; $xn=$1; } \
($d)=m/$dt=([^;\s]+)/; push @g,$_; END{putg();} \
sub putg {if(@g){ if($u5>2 or $u3>2) { $un=$u5+$u3; $g[0]=~s/$/;utrx=$u5,$u3/; } print @g;} @g=(); $u5=$u3=0; }'\
 > $gset.incutcd.gff

echo "/bin/rm $gset.incut*.tmp"

