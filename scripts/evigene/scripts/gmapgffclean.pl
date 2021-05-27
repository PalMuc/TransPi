#!/usr/bin/env perl
# gmapgffix.pl

## FIX gmap CDS bug: last 3' CDS end is -2 short, revcomp also ## gmap version? or all?

$notarget=1; #opt?
$src=$ENV{src};
$firstonly= $ENV{first} || 0; # only mrna1
$MINCP= $ENV{cov} || $ENV{Coverage} || 0; # filter out weak aligns?
$nopr=0;

while(<>) {
  s/Name=[^;]+;//; s/^###/#/; 
  unless(s/\.(path|mrna)1\b//g) { s/\.(path|mrna)(\d+)/p$2/g; next if($firstonly);  }
  s/\t\S+/\t$src/ if($src);
  if(/\tgene/){s/.*//;} 
  elsif(/\t(CDS|exon)/) { 
    s/(ID)=[^;\n]+;?//g; 
    s/Target=/targ=/ if($notarget);
    } 
  elsif(/\tmRNA/){ 
    s/Parent=/gene=/; ($cv,$pi)=m/Coverage=([^;\s]+);Identity=([^;\s]+)/; ($d)=m/ID=([^;\s]+)/; 
    ($gn,$dp)= ($d =~ m/^(.+)p(\d+)$/) ? ($1,$2) : ($d,0);
    $nopr=0; if($cv and $pi) {
    $cp=int(($cv+0.5) * ($pi+0.5)/100); $cp=100 if($cp>100); s/\t\.\t/\t$cp\t/; $d=~s/p\d+$//; 
    $nopr=($MINCP and $cp < $MINCP)?1:0;
    if(!$nopr and $did{$gn} > $cp) { $nopr=1; }  # drop 2ndaries if lesser
    $did{$gn}=$cp unless($nopr);
    }
     ## $s=$sc{$d}; s/$/;vscore=$s/ if($s); 
  }
  print unless($nopr);
}


=item filter junk models

## use CDS:1 score if .gff has CDS; dont use score: unless that is meaningful
#...  -mrna mRNA -exon 'CDS,exon'

 gzcat $rnagenes.gff.gz | $evigene/scripts/overbestgene2.perl -in stdin \
 -alttr -typeover exon -noOVEREXON2 \
 -scoretype='many.score:3,CDS:2,UTR:1' \
 -genescore  -trivial 10  -pctover 10  -summarize -noskip \
  > $rnagenes.best1.gff

=cut

