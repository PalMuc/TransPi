#!/usr/bin/env perl
# cuffaa2cds.pl

=item info

** instead of this use genefindcds.pl and regenerate proteins, gff+cds
   from input cdna=cuff.tr gff=cuff.gff

 $evigene/scripts/genefindcds.pl \
  -genes daphmag_rnaseq.cuff2.gff.gz \
  -cdna daphmag_rnaseq.cuff2.tr \
  -dna $dmag/genome/dmagna20100422assembly.fa \
  > daphmag_rnaseq.cuff2cds.gff
    
grep '^>' daphmag_rnaseq.cuff2.aa | perl -ne's/>//; @v=split /;?\s+/, $_; ($d,$al,$cl,$or,$ofs,$loc)=@v;
 print join("\t",@v[0..4])."\n";' > daphmag_rnaseq.cuff2aa.tab

# daphmag2cuf8Gcontig00001.15.1   aalen=33,complete       clen=315        strand=-        offs=126-227
# daphmag2cuf8Gcontig00001.15.2   aalen=39,complete       clen=826        strand=-        offs=23-142
# daphmag2cuf8Gcontig00003.5.1    aalen=62,partial        clen=189        strand=+        offs=2-187
# daphmag2cuf8Gcontig00116.16.1   aalen=121,complete      clen=891        strand=-        offs=331-696

# paste protein info and CDS offsets into cufflinks.gff (no CDS).

=cut

# cat daphmag_rnaseq.cuff2.gff.gz | cat daphmag_rnaseq.cuff2aa.tab - | perl -ne

while(<>) {
if(/^daphmag2cuf/){ 
  ($d,@v)=split; ($id=$d)=~s/daphmag2cuf8//; map{ ($k,$v)=split"="; $an{$id}{$k}=$v; } @v; }
elsif(/^\w/) { # gff
  ($d)=m/(?:ID|Parent)=([^;\s]+)/; %anv=%{$an{$d}}; @gf=split"\t"; ($r,$t,$b,$e,$o)=@gf[0,2,3,4,6];
  if($o eq ".") { $gf[6]=$anv{strand}; }  
  if(/\tmRNA/) { $al=$anv{aalen}; $x=$anv{clen}; $a3=3*$al; $cx="$a3/$x";
    ($ob,$oe)=$anv{offs}=~m/(\d+).(\d+)/; ($tl,$tb,$te)=(0,0,0); $gf[8]=~s/$/;aalen=$al;cxlen=$cx/; } 
  elsif(/\texon/) { 
    $w=1+$e-$b; $tl+=$w; 
    if($tl <= $oe and $tl >= $ob) { $bi=$ob-$tb; $ei=($cb,$ce)=($xb+$bi,$xe-$ei); } 
    }
  print @gf; 
}
