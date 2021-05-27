#! /bin/bash
### qsub -q normal trasmspan_kfish2sd2.sh
#PBS -N trasmspan
#PBS -A ind114
#PBS -l nodes=1:ppn=32,walltime=23:55:00
#PBS -o trasmspan.$$.out
#PBS -e trasmspan.$$.err
#PBS -V

ncpu=32

## /tmp space sort failes..
export TMPDIR=$HOME/scratchn/tmp
export LC_ALL=C 

# trasmspan_kfish2sd2.sh
# k1: fungr1 map to funhe2 with -grepread mismatch:2 read filter ; drops ~10% mapped reads
# k2: no grepread mismatch; trasmspan_kfish2sd2.config opts for velvetg,oases
dv=k2
groupnam="fungr1" # get from bams
config=trasmspan_kfish2sd2.config
strand=fwd
# strand=rev
#OLDoptrasm="-ver $dv -grepread mismatch:2 -allstrandspan"
optrasm="-ver $dv -allstrandspan"

spantable=kfish2_strandspan5m.tab
# lastids: ptf093       ptr083
firstpart=1; lastpart=94
subd=strasm2

datad=$HOME/scratchn
workd=$datad/chrs/kfish/
rund=$workd/rnas/$subd
sdir=$HOME/bio/evigene/scripts/rnaseq

cd $rund
## add *.chr for nooimp.chr.bam and noo.chr.bam
bamlist=`ls $groupnam.{$strand,noo}*.chr.bam`

echo "start trasmspan : `date`"  
i=0; # cpu counter
j=$firstpart; # data part counter
  j=$(( $j - 1 ))

while [ $j -lt $lastpart ]; do
{
  j=$(( $j + 1 ))
  nam=trasm${dv}p$j
  rlog=$nam.$strand.log
  # continue if bestx.gff exists.. NEED j++
  if test -f $nam/vel${dv}pt*.$strand.best1.gff; then echo "done $nam"; continue; fi

  echo $sdir/trasmstrandspan.pl -part=$j -strand=$strand -out $nam $optrasm \
    -config=$config -span=$spantable -bams $bamlist  
  $sdir/trasmstrandspan.pl -part=$j -strand=$strand -out $nam  $optrasm  \
    -config=$config -span=$spantable -bams $bamlist  >& $rlog  &

  i=$(( $i + 1 ))
  if [ $i -ge $ncpu ]; then wait; i=0; fi

} done
wait

echo "end trasmspan: `date` "

#........... span table  .................
# $evigene/scripts/rnaseq/strandspanintab.pl -partwin 5000000 \
#   -in $kfish2/intron/intron_good.gff.gz -chr ../bam1fungr/allfungr1-kfish2.chrcount -zerospan \
#    > kfish2_strandspan5m.tab
# #allfungr1-kfish2.chrcount is scaffold-size table with read counts/scaffold; read cnt not needed
# #kfish2_strandspan5m.tab
# #Chr    Beg     End     Nfwd    Nrev    Nno     Idfwd   Idrev   Or
# Scaffold1       1       28900   0       0       0       ptf001  0       0
# Scaffold1       28901   80200   6801    0       0       ptf001  0       f
# Scaffold1       80201   92800   0       208     0       0       ptr001  r
# Scaffold1       92801   103900  432     0       0       ptf001  0       f

#................ config file ..................
# trasmspan_kfish2sd2.config
# # trasmspan.config for sdsc
# 
# # kfish2/fungr1 rna: 
# # kmer1=55 45 35 25 21
# kmers=53 43 33 29 25 21
# optvelvet=-ins_length 200
# optoases=-min_pair_count 2 -scaffolding yes -min_trans_lgth 180 -ins_length 200
# #def.bestscore=inqual:5,cov:8,pid:5,CDS:2,UTR:1,cdsindel:-2
# bestscore=aalen:9,inqual:5,cov:4,pid:2,CDS:2,UTR:1
# 
# ## should enable ENV params below; $ENV{HOME} is good
# MINLEN=180 
# samtools=/home/ux455375/bio/bin/samtools 
# ## note velvet127s not older v123
# velbin=/home/ux455375/bio/velvet127s/bin2/ 
# gmapdb=/home/ux455375/scratchn/chrs/kfish/genome/gmap12
# gmap=/home/ux455375/bio/gmap1206/bin/gmap
# evigene=/home/ux455375/bio/evigene
# genomed=/home/ux455375/scratchn/chrs/kfish/genome/
# dgenome=killifish20121129asm
# introns=intron_good.gff  
