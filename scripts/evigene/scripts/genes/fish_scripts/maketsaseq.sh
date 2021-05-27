#!/bin/bash
# maketsaseq.sh ; add maketsannot.sh ..

evigene=/bio/bio-grid/mb/evigene 
# do for pt=mrna,aa,cds
cd pubgenes/  

tonam=kfish2rae5hx11; isx=0; 
# isgr=1; isgr=0;
for isgr in 0 1; do {

todir=../pubtsafunhe; if [ $isgr = 1 ]; then todir=../pubtsafungr; fi
mkdir $todir;

for ty in mrna aa cds ann; do {
  echo "#TO: " $todir/$tonam.$ty
gunzip -c kfish2rae5h.{main,alt}.pub.$ty.gz | cat kfish2rae5g.pubattr.txt - |\
  env isx=$isx isgr=$isgr ty=$ty perl -ne \
'BEGIN{ $TY=$ENV{ty}; $ISX=$ENV{isx}; $ISGR=$ENV{isgr}; }
if(/^Public/) { print if($TY eq "ann");  next; } 
elsif(/^Funhe/) { ($pd,$od,@an)=split"\t"; @od=split",",$od; 
($xd)=grep/Funhe2Exx11m/,@od; $xd=$pd if($ISX); $isgr=($od=~m/Fungr/)?1:0; 
if($xd and $isgr == $ISGR) { if($xok{$xd}) { $xdup{$xd}.="$pd,"; }
else { $pok{$pd}=$xd; $xok{$xd}=$pd; $od=~s/$xd/$pd/; 
print join("\t",$xd,$od,@an) if($TY eq "ann"); } } } 
else { if(/^>(\S+)/) { $pd=$1; $ok=$xd=$pok{$pd};  
if($xd and $xd ne $pd) { s/>$pd/>$xd/; s/oid=/pubid=$pd; oid=/;} } 
print if $ok; }' \
    > $todir/$tonam.$ty

} done
# ty
mv $todir/$tonam.ann $todir/$tonam.ann.txt

} done 
# isgr

cd ../
# -----------------

# extras  
cd pubmixx11/ 

tonam=kfish2evgxx11; isx=1; 
# isgr=0; isgr=1;
for isgr in 0 1; do {
todir=../pubtsafunhe; if [ $isgr = 1 ]; then todir=../pubtsafungr; fi

for ty in mrna aa cds ann; do {
  echo "#TO: " $todir/$tonam.$ty
gunzip -c kfish2evgxx11.$ty.gz | cat kfish2evgxx11.mann.notpub5h - |\
  env isx=$isx isgr=$isgr ty=$ty perl -ne \
 'BEGIN{ $TY=$ENV{ty}; $ISX=$ENV{isx}; $ISGR=$ENV{isgr}; }
if(/^Public/) { print if($TY eq "ann");  next; } 
elsif(/^Funhe/) { ($pd,$od,@an)=split"\t"; @od=split",",$od; 
($xd)=grep/Funhe2Exx11m/,@od; $xd=$pd if($ISX); $isgr=($od=~m/Fungr/)?1:0; 
if($xd and $isgr == $ISGR) { if($xok{$xd}) { $xdup{$xd}.="$pd,"; }
else { $pok{$pd}=$xd; $xok{$xd}=$pd; $od=~s/$xd/$pd/; 
print join("\t",$xd,$od,@an) if($TY eq "ann"); } } } 
else { if(/^>(\S+)/) { $pd=$1; $ok=$xd=$pok{$pd};  
if($xd and $xd ne $pd) { s/>$pd/>$xd/; s/oid=/pubid=$pd; oid=/;} } 
print if $ok; }' \
    > $todir/$tonam.$ty
 
} done
# ty
mv $todir/$tonam.ann $todir/$tonam.ann.txt
} done
# isgr

cd ../
# -----------------

# merge to submit set input
cd pubtsafunhe/
allnam=kfish2evg5hx11fh
for ty in mrna aa cds ann.txt; do {
  cat {kfish2rae5hx11,kfish2evgxx11}.$ty > $allnam.$ty
} done
cd ../

cd pubtsafungr/
allnam=kfish2evg5hx11fg
for ty in mrna aa cds ann.txt; do {
  cat {kfish2rae5hx11,kfish2evgxx11}.$ty > $allnam.$ty
} done
cd ../

# FIX fungr seq bad codes..
# FIXME2, need new aa cdna_bestorf.pl for changed mrna .. output only changed mrna, 
#  $evigene/scripts/cdna_bestorf.pl  -aa -cds -cdna ${allnam}fix.$ty
cd pubtsafungr/
allnam=kfish2evg5hx11fg
ty=mrna 
if [ ! -f $allnam.$ty.codefixed ]; then 
  cat $allnam.$ty | perl -ne \
  'if(/^>(\S+)/){ $d=$1; puts() if($fa); $hd=$_; $id=$d; $fa="";} else{ chomp; $fa.=$_;} END{puts();} 
  sub puts{ $fa=uc($fa); chomp($hd); $nx= $fa=~tr/ACGTN//c; 
  if($nx){ $fa=~tr/ACGTN/n/c; $hd.=" badcode=$nx;"; 
  $fa=~s/(.{60})/$1\n/g;  print "$hd\n$fa\n"; $fa=$hd=""; } }' \
   > ${allnam}fix.$ty
   
   #? -action=fwd ?, maybe add -outmrna to get proper direction
  $evigene/scripts/cdna_bestorf.pl -nostop -noutrorf -outaa -outcds -outmrna -cdna ${allnam}fix.$ty

  for ty in mrna aa cds; do {
    mv $allnam.$ty  $allnam.$ty.nofix
    perl -ne 'if(/^>(\S+)/){ $d=$1; $ok=($did{$d}++)?0:1; } print if $ok;' \
      ${allnam}fix.$ty $allnam.$ty.nofix > $allnam.$ty
  } done
  touch $allnam.$ty.codefixed ; gzip --fast $allnam.$ty.nofix 
  
fi
cd ../

## need also correct -strand mRNA/aa : handfull, do from id list

# ... old ..
# for ty in mrna cds; do {
#   if [ ! -s $allnam.$ty.nofix ]; then 
#   cat $allnam.$ty | perl -ne \
# 'if(/^>(\S+)/){ $d=$1; puts() if($fa); $hd=$_; $id=$d; $fa="";} else{ chomp; $fa.=$_;} END{puts();} 
# sub puts{ $fa=uc($fa); chomp($hd); $nx= $fa=~tr/ACGTN//c; 
# if($nx){ $fa=~tr/ACGTN/n/c; $hd.=" badcode=$nx;"; }
# $fa=~s/(.{60})/$1\n/g;  print "$hd\n$fa\n"; $fa=$hd=""; }' \
#  > ${allnam}fix.$ty
# 
#   mv $allnam.$ty  $allnam.$ty.nofix
#   mv ${allnam}fix.$ty $allnam.$ty 
#   touch $allnam.$ty.codefixed ; gzip --fast $allnam.$ty.nofix 
#   fi
# 
# } done # ty
