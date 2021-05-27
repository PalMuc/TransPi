#!/bin/bash
# asmdupfilt2.sh

clsuf=class4
alsuf=aln4
tinyaln=35

## add stats below ..

if [ "X" = "X$evigene" ]; then
evigene=/bio/bio-grid/mb/evigene/
fi

## name_cdenn.clstr input list
trset=$*

## in names not same:
# dmag4vel4xca_cd90.aa dmag4vel4xca_cd90.aa.clstr dmag4vel4xca_cd90.aa.qual   
# dmag4vel4xca_cde35.cds  dmag4vel4xca_cde35.cds.clstr
#........
## FIXME WRONG aa.qual from cdhit, need all.aa.qual from input.aa.gz
# ERR: missing daphmag5nalsoap_cde35.aa.qual
# ERR: missing daphmag5xc3soap_cde35.aa.qual

for cdscl in $trset; do {
  pt=`echo $cdscl | sed 's/.gz//; s/\.cds.*//; s/\.clstr//;'`
  #WRONG for aasize; WRIGHT for clstr# ptaacl=`echo $pt | sed 's/_cde.*//; s/$/_cd90/;'`
  ptaa=`echo $pt | sed 's/_cde.*//;'`

  if [ ! -f $pt.$clsuf ]; then 
    aaf=$ptaa.aa.qual
    aaclstr=${ptaa}_cd90.aa.clstr
    if [ ! -f $aaf ]; then
      aaf=$pt.aa.qual
      aaclstr=${pt}_cd90.aa.clstr
    fi

    if [ -f $aaf ]; then
      echo "# asmrna_dupfilter2 $pt"

      $evigene/scripts/rnaseq/asmrna_dupfilter2.pl -debug -CDSALIGN -tinyaln $tinyaln -aasize $aaf \
        -acdhit $aaclstr -bcdhit $cdscl -outeq $pt.$alsuf -outclass $pt.$clsuf 

    else
      echo "ERR: missing $pt.aa.qual";
    fi
    echo "#..."
  fi

  if [ -f $pt.$clsuf ]; then
echo "# Class Table for $pt.$clsuf "
cat $pt.$clsuf | cut -f2,3 | sort | uniq -c | perl -ne \
's/altmidfrag/amfrag/; s/maybeok/okay/; ($n,$ac,$cl)=split; for $ct ($cl,"total") { $tab{$ct}{$ac}+=$n; } $nt+=$n; 
END{ @ac=qw(okay drop); printf "%-9s\t","class"; print join("\t",@ac,@ac)."\n"; 
foreach $cl (sort keys %tab) { @pv=@nv=(); 
 foreach $ac (@ac) { $n=$tab{$cl}{$ac}||0; $p=int(1000*$n/$nt)/10; push @pv,$p; push @nv,$n; }
do{ print (("-") x 45); print"\n"; } if($cl eq "total");
printf "%-9s\t",$cl; print join("\t",@pv,@nv)."\n"; } }'
echo 
  fi

} done


