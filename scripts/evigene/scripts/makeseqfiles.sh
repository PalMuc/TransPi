#!/bin/tcsh

# Make predicted protein and transcript files
set dgenome=../genome/aphid2asm.fa

foreach agf ( aphid2_mix3.gff.gz )

  set aan=`basename $agf .gff.gz`
  
  # if( 0 == 1 ) then
  # if (-f $aan.aa.gz ) continue
  gzcat $agf | grep '	mRNA' | perl -ne\
   '($d)=m/ID=([^;\s]+)/; ($aa)=m/protein=([^;\s]+)/; \
   ($al)=m/(aalen=[^;\s]+)/; ($ho)=m/ho3=[\d\.\/]+([^;\s]+)/; \
   ($pro)=m/pro1=[\d\.\/]+.([^;\s]+)/; \
   $ho.=","if($ho && $pro); $ho.=$pro; $ho=" Dbxref=$ho" if $ho; \
   if($aa){ $aa=~s/(.{0,60})/$1\n/g; chomp($aa); print ">$d $al$ho\n$aa\n";}' \
   > $aan.aa
  # endif

  ../scripts/cdsgff2genbank.pl -pretty -t=exon -a dna -gff $agf -fasta $dgenome > $aan.tr
  ## fixme:
  # perl -pi -e'if(/^>/){ s/;(est|rseq|intr|pasa|pro|pParent|ni)=[^;\s]+//g; }' $aan.tr

end

