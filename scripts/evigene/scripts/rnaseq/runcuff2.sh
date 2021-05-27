#!/bin/tcsh

#? setenv LD_LIBRARY_PATH /usr/local/lib:/usr/sfw/lib:/usr/lib 

set bin_dir=/bio/bio-grid/mb/tophat/bin

## better to make subdir for each .sam, so can run parallel
## dang cuffl for fixed names; no stdin
foreach sam (`/bin/ls -t *.sam`)
  set gp=`echo $sam | sed -e's/.sam//;'`
  if( -f $gp.genes.expr ) continue

  $bin_dir/cufflinks $sam > & log.cuf.$gp

  cat transcripts.gtf | perl -ne\
  's/\ttranscript/\tmRNA/; ($r,$b,$e,$v)=(split)[0,3,4,5]; s/Cufflinks/Gsnap3/; s/CUFF/GSNA3/g; \
  s/ frac .*$//; $dt=(/\tmRNA/)?"ID":"Parent";  s/gene_id \S+ //; s/transcript_id /$dt=/; s/"//g; \
  s/ exon_number /xi=/; s/ FPKM (\S+)//; $k=$1; $k=($k>999)?int($k):sprintf"%.3g",$k;  \
  $w=1+$e-$b; s/\t$e\t$v/\t$e\t$k/;  print; BEGIN{print"##gff-version 3\n"; }'\
   > transcripts.gff

  mv genes.expr $gp.genes.expr
  mv transcripts.expr $gp.transcripts.expr
  mv transcripts.gtf $gp.transcripts.gtf
  mv transcripts.gff $gp.transcripts.gff

end

