#!/bin/bash -l
# gg-velrun.sh

velbin=$HOME/bio/velvet/bin
gmapbin=$HOME/bio/gmap/bin
workd=$HOME/scratch/chrs/aphid2
dgenome=aphid2asm

## filter out toobig sets : 830 Mb for now
sz=`du -sk | sed -e's/.\..*//;'`
if [ $sz -gt 830000 ]; then
  echo "ERROR: too big size = $sz"
  exit
fi

## dang bash; set nullglob for missing files to velveth: *.longfa 
shopt -s nullglob
echo Virt mem limit set 25000000 kbytes == 25 Gb
ulimit -v 25000000
ulimit -v
echo '#-------------------'

echo "#.. start rnavelv : `date`"
echo "#.. velveth"
$velbin/velveth vel 21 -fasta -shortPaired *.fa2 -short *.fa1 *.fa  -long *.longfa 
echo "#.. velvetg"
$velbin/velvetg vel -read_trkg yes 
echo "#.. oases"
$velbin/oases vel -ins_length 200 

if [ -f vel/transcripts.fa ]; then
echo "#.. gmap : `date`"
perl -pi -e'if(/^>/){ s,/(\d+)_Confidence_(\d\.\d+), nt=$1; score=$2;,g; s/Locus_/Loc/; s/_Transcript_/t/; }' \
  vel/transcripts.fa
$gmapbin/gmap -D $workd/genome/gmap -d $dgenome -n 4 -f 2 vel/transcripts.fa > vel/transcripts.gff 2> /dev/null

perl -pi.old -e's/Name=[^;]+;//; s/^###/#/; s/\.(path|mrna)(\d+)/p$2/g; if(/\tgene/){s/.*//;} 
elsif(/\t(CDS|exon)/) { s/(ID)=[^;\n]+;?//g; }
elsif(/\tmRNA/){ s/Parent=/gene=/; ($cv,$pi)=m/Coverage=([^;\s]+);Identity=([^;\s]+)/; ($d)=m/ID=([^;\s]+)/;
$cp=int($cv * $pi/100); s/\t\.\t/\t$cp\t/; $d=~s/p\d+$//; $s=$sc{$d}; s/$/;vscore=$s/ if($s); }              
BEGIN{open(S,"grep score= vel/transcripts.fa |");while(<S>){($d)=m/>(\S+)/; ($s)=m/score=([\d\.]+)/; $sc{$d}=$s;}}'\
  vel/transcripts.gff

else
  echo "ERROR: no vel/transcripts.fa"
fi

echo "#.. end rnavelv : `date`"

