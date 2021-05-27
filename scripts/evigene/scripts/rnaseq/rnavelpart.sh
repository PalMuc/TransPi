#! /bin/bash -l
### qsub -q normal rnavelv.sh
#PBS -N rnavel4
#PBS -A TG-MCB100147
## for -q normal ember.ncsa ncpus=6 6cpu/node mem=NNNmb
#PBS -l ncpus=6,mem=63gb,walltime=10:55:00
#PBS -o velv4.out
#PBS -e velv4.err
#PBS -V

# workd=$HOME/scratch/chrs/aphid2
workd=/export/udisk3/work/aphid
dgenome=aphid2asm
# datad=$workd/rnas/velv
# velbin=$HOME/bio/velvet/bin
# gmapbin=$HOME/bio/gmap/bin
velbin=/bio/bio-grid/mb/rnaseq/velvet/bin
gmapbin=/bio/bio-grid/mb/gmap10b/bin

echo "#.. start rnavelv : `date`"
# cd $datad

## FIXME suffixes are a mess now: need fa2 OR pairs, not both; ditto fa1/unpair
echo "#.. velveth"
$velbin/velveth vel 21 -fasta -shortPaired *.fa2 -short *.fa1 *.fa  -long *.longfa 
echo "#.. velvetg"
$velbin/velvetg vel -read_trkg yes 
echo "#.. oases"
$velbin/oases vel -ins_length 200 

echo "#.. end rnavelv : `date`"

## rename transcripts.fa first: want score in ID so gmap.gff can parse it??
perl -pi -e'if(/^>/){ s,/(\d+)_Confidence_(\d\.\d+), nt=$1; score=$2;,g; s/Locus_/Loc/; s/_Transcript_/t/; }' \
  vel/transcripts.fa

echo "#.. gmap"
$gmapbin/gmap -D $workd/genome/gmap -d $dgenome -n 4 -f 2 vel/transcripts.fa > vel/transcripts.gff 2> /dev/null

perl -pi.old -e's/Name=[^;]+;//; s/^###/#/; if(/\tgene/){s/.*//;} elsif(/\t(CDS|exon)/) { s/(ID)=[^;\n]+;?//g; } 
elsif(/\tmRNA/){ s/Parent=/gene=/; ($cv,$pi)=m/Coverage=([^;\s]+);Identity=([^;\s]+)/; ($d)=m/ID=([^;\s]+)/;
$cp=int($cv * $pi/100); s/\t\.\t/\t$cp\t/; $s=$sc{$d}; s/$/;vscore=$s/ if($s); } s/\.(path|mrna)(\d+)/p$2/g;
BEGIN{open(S,"grep '^>' vel/transcripts.fa |");while(<S>){($d)=m/>(\S+)/; ($s)=m/score=([\d\.]+)/; $sc{$d}=$s;}}'\
  vel/transcripts.gff

echo "#.. end rnavelv gmap : `date`"
# add final steps to collate rparts/.../vel/transcripts.{fa,gmap.gff} for genome.

