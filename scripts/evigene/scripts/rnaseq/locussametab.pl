#!/usr/bin/env perl
# sameloctab.pl

=item usage locussametab.pl  

  tabulates evgR diff loci that are same locus in genome map, using various inputs
  NEXT STEP: 
    evigene/scripts/rnaseq/locussamediff2pubid.pl 
    updates evigene publicset idclass table from genome-mapping inputs 
      diffloci.tab (evgR alts are not alts)  
      sameloci.tab (evgR diffloci are alts)  
  
  set for kfish2 evigeneR gmapped transcripts
  
  sameloctab.pl -eqgene $nam.all.eqgene -idtab $nam.newids \
    -diffloc $nam.diffloci.eqalt2main -alntab $nam.alnsense.tab \
      > $nam.sameloci.tab

  idtab : mrnaid, geneid, aasize, aaref optional
  eqgene: equal gene table,  $evigene/scripts/equalgene.pl -in nam.gff -ov nam.gff > nam.eqgene
  alntab: alignment table, $evigene/scripts/ests/gff2aligntab.pl  nam.gff > nam.align.tab
  diffloci: table of alts that map elsewhere from main tr (removed from sameloc tests)
      equalgene -in alts.gff -ov mains.gff > alt2main.eqgene  
      cat alt2main.eqgene  | sort -k1,1 -k3,3 | perl -ne \
'($td,$od,$ad,$loc)=split; if($ltd eq $td) { $ad=$lad if($ad eq "na"); } elsif($ltd and $ltd ne $td) { \
print join("\t",$ltd,$lod,$lad,$lloc)."\n" if($lad eq "na"); } ($ltd,$lod,$lad,$lloc)=($td,$od,$ad,$loc);' \
 > $nam.diffloci.eqalt2main

=item result 

  kfish2/rnas/kf2evgr/trevg367mixx/alt11/
  diffsame/kf2pub11_both.sameloci9.tab
    classes: 26782 total, 13168 notsamer, 9832 sameloc, 3782 samelocsw

  sameloc    Funhe2Exx11m016071t12   Funhe2Exx11m046585t2    220     124     CDD:220496,145  ,
  notsamer   Funhe2Exx11m009703t2    Funhe2Exx11m047099t1    562     132     UniRef50_P00740,1918    UniRef50_Q19UG5,307
  notsamer   Funhe2Exx11m009703t2    Funhe2Exx11m047099t2    562     125     UniRef50_P00740,1918    UniRef50_Q19UG5,307
  sameloc    Funhe2Exx11m016605t4    Funhe2Exx11m046202t1    100     133     UniRef50_A2VDW6,369     ,
  samelocsw  Funhe2Exx11m021998t1    Funhe2Exx11m021276t2    283     209     ,       ,
  samelocsw  Funhe2Exx11m021998t1    Funhe2Exx11m021276t3    283     162     ,       ,
  
  diffloci.eqalt2main : n=5739 merge of two ways (splits confuse this)
  
  
=item see also  

  kfish2/rnas/kf2evgr/trevg367mixx/alt11/pubgff.info 
    : describes this and several other steps to make clean, public evgRna mapped gene gff set
    
  evigene/scripts/rnaseq/locusgmapsub.sh 
    : simpler/similar steps, used for pine tr set
    
=cut

use strict;
use Getopt::Long;

my $MINCDS=20;
my $MINALN=85;

my($eqgene,$alntab,$diffloci,$idtable); # = @ARGV;
$eqgene="kf2pub11_both.all.eqgene";
$alntab="kf2pub11_both.alnsense.tab";
$diffloci="kf2pub11_both.diffloci.eqalt2main2";
$idtable="kf2pub11_both.newids";

my $optok= GetOptions(
  "eqgene=s", \$eqgene, 
  "alntab=s", \$alntab,  
  "diffloci=s", \$diffloci, 
  "idtable=s", \$idtable, 
);

my(%okmap,%poormap,%difflocus,%aasize,%aaref,%garef);

open(F,$alntab) or die $alntab;
while(<F>) {
  my($id,$cov,$pid)=(split"\t")[3,7,8]; my $coi=$cov*$pid/100;  
  $okmap{$id}++ if($coi>=$MINALN);  
  # $poormap{$id}++ if($coi<$MINALN and not $okmap{$id});
  # print "$id\t\n" if($coi<$MINALN and not $okmap{$id});
} close(F);

open(F,$diffloci) or die $diffloci;
while(<F>) {
  my($id,$oid,$othrd,$loc,$cla)=split"\t";
  $difflocus{$id}=$cla if($cla =~ /Df\d/);
} close(F);

open(F,$idtable) or die $idtable;
while(<F>) { next unless(/^\w/); 
  my @v=split"\t"; my($id,$gd,$aq)=@v[0,2,5]; 
  my ($aw)= $aq=~m/(\d+)/;  $aasize{$id}=$aw; 
  my ($arv,$ar)= m/aaref:(\d+),([:\w]+)/;  
  if($ar) { $aaref{$id}=$ar; $garef{$gd}{$ar}+=$arv; } 
} close(F);


open(F,$eqgene) or die $eqgene;
while(<F>) {
  my ($td,$od,$ad,$loc)=split"\t"; (my $tg=$td)=~s/t\d+$//; 
  next unless( $okmap{$td} and not $difflocus{$td} );
  my @ad= grep{ not /^$tg/ } split",",$ad; 
  my @sd= grep /\w/, map{ my($d,$p)=split"/"; (($p=~/^[CI]/ or $p>=$MINCDS))?$d:""; } @ad; 
  @sd= grep { $okmap{$_} and not $difflocus{$_} } @sd;
  foreach my $sd (@sd) { 
    my($md,$ad)= sort ($td,$sd); 
    my($clsame,@sameval)= aaorder($md,$ad);
    print join("\t",$clsame,@sameval)."\n"; # if($clsame =~ /^same/);
    } 
} close(F);

# warn "#stats: ...\n";

sub aaorder {
  my($td,$sd)= @_;
  my($tr,$sr,$trc,$src,$taw,$saw,$v,$diffr,$swap);
  my($tg,$sg)=map{ ($v=$_)=~s/t\d+$//; $v; } ($td,$sd); 
  $tr=$aaref{$td}; $sr=$aaref{$sd}; 
  $trc=$garef{$tg}{$tr}; $src=$garef{$sg}{$sr}; 
  $taw=$aasize{$td}; $saw=$aasize{$sd};
  $diffr=($tr and $sr and $tr ne $sr)?1:0;  
  $swap=($src>$trc)?1:($src<$trc)?0:($saw >$taw)?1:0; 
  
  my @val=($sd,$td,$saw,$taw,"$sr,$src","$tr,$trc");
  @val=@val[1,0,3,2,5,4] unless($swap);  
  if($diffr) { return ("notsamer",@val);  } 
  elsif($swap){ return ("samelocsw",@val); } 
  else { return("sameloc",@val); } 
}


__END__

cat $nam.all.eqgene | perl -ne 'BEGIN{$MINCDS=20;}\
($td,$od,$ad,$loc)=split"\t"; ($tg=$td)=~s/t\d+$//; @ad=grep{ not/^$tg/ } split",",$ad; \
@sd=grep /\w/, map{ ($d,$p)=split"/"; ($p=~/^[CI]/ or $p>=$MINCDS)?$d:""; } @ad; \
foreach $s (@sd) { @d=sort ($td,$s); print join("\t","samelocus",@d)."\n"; }  ' | sort -u \
 > $nam.sameloci.eqall

cat kf2pub11_both.alnsense.tab | perl -ne'($id,$cov,$pid)=(split)[3,7,8]; $coi=$cov*$pid/100; \
$ok{$id}++ if($coi>85);  print "$id\t\n" if($coi<85 and not $ok{$id});' | sort -u | \
 ggrep -v -F -f - kf2pub11_both.sameloci.eqall > kf2pub11_both.sameloci.eqall.goodmap

grep Df0 kf2pub11_both.diffloci.eqalt2main2 | cut -f1 | sed 's/^/samelocus  /; s/$/ /;' | \
 ggrep -v -F -f - kf2pub11_both.sameloci.eqall.goodmap > kf2pub11_both.sameloci.eqall.goodnodf

## FIXME correct sameloci table w/ aasize,aaref : dont join tr of diff aaref ?
cat kf2pub11_both.newids kf2pub11_both.sameloci.eqall.goodnodf | perl -ne\
'if(/^samelocus/) { ($ssl,$td,$sd)=split; ($tg,$sg)=map{ \
($v=$_)=~s/t\d+$//; $v; } ($td,$sd); $tr=$ar{$td}; $sr=$ar{$sd}; \
$trc=$gar{$tg}{$tr}; $src=$gar{$sg}{$sr}; $taw=$aw{$td}; $saw=$aw{$sd}; \
$diffr=($tr and $sr and $tr ne $sr)?1:0;  $swap=($src>$trc)?1:($src<$trc)?0: \
($saw >$taw)?1:0; @val=($sd,$td,$saw,$taw,"$sr,$src","$tr,$trc"); \
@val=@val[1,0,3,2,5,4] unless($swap);  if($diffr) { print \
join("\t","notsamer",@val)."\n";  } elsif($swap){ print \
join("\t","samelocsw",@val)."\n"; } else { print \
join("\t","sameloc",@val)."\n"; } } elsif(/^\w/) { @v=split"\t"; \
($id,$gd,$aq)=@v[0,2,5]; ($aw)=$aq=~m/(\d+)/; ($arv,$ar)=m/aaref:(\d+),(\w+)/; \
if($ar) { $ar{$id}=$ar; $gar{$gd}{$ar}+=$arv; }  $aw{$id}=$aw; }' \
> kf2pub11_both.sameloci.swaps.good