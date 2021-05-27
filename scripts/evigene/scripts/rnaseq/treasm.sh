#! /bin/bash
### env oid=dmag4vel4xfik55Loc5429t2 qsub -q shared pulltreads.sh
#PBS -N treasm.sh
#PBS -A ind114
#PBS -l vmem=8gb,nodes=1:ppn=2,walltime=2:55:00
#PBS -o treasm.$$.out
#PBS -e treasm.$$.err
#PBS -V

# pulltreads.sh

peset="concordant_mult|concordant_uniq|paired_mult|paired_uniq_long|paired_uniq_inv|paired_uniq_scr"
srset="halfmapping_mult|halfmapping_uniq|unpaired_mult|unpaired_uniq"
noset="nomapping"

workd=$HOME/scratch/chrs/daphmag/rnas
odir=$workd/dmag5trfine
cd $workd

# bamall=$workd/gsodm13nd_srt/nodamx*-gsnap-*s.bam
bamall="$workd/gsodm13nd_srt/nodamx*-gsnap-*s.bam $workd/gsodm13sd_srt/Dman*-gsnap-*s.bam"

if [ "X" = "X$oid" ]; then echo env oid=XXX; exit -1; fi

module add samtools

cd $odir
pesam=$oid.pe.sam
srsam=$oid.sr.sam
touch $pesam
touch $srsam
npe=0; nsr=0; nf=0;

inbams=`ls $bamall | egrep "$peset"`
osam=$pesam
for ibam in $inbams; do {
  nam=`echo $ibam | sed 's/\.bam//;'`
  nc=`grep $oid $nam.count | cut -f3`
  if [ $nc -gt 0 ]; then
    npe=$(( $npe + $nc ))
    samtools view $ibam $oid >>  $osam
  fi
} done


inbams=`ls $bamall | egrep "$srset"`
osam=$srsam
for ibam in $inbams; do {
  nam=`echo $ibam | sed 's/\.bam//;'`
  nc=`grep $oid $nam.count | cut -f3`
  if [ $nc -gt 0 ]; then
    nsr=$(( $nsr + $nc ))
    samtools view $ibam $oid >>  $osam
  fi
} done


# veltreasm.sh
#.................
#! /bin/bash
### env oid=dmag4vel4xfik55Loc5429t2 qsub -q shared veltreasm.sh
#PBS -N veltreasm
#PBS -A ind114
#PBS -l vmem=8,nodes=1:ppn=2,walltime=2:55:00
#PBS -o veltreasm.$$.out
#PBS -e veltreasm.$$.err
#PBS -V

##.. velv reasm takes less time than pulltreads  .. k31 single kmer ok? NO
## ** PROBLEMS w/ vel reasm .. can merge mistakes/junk

## FIXME: 1st test all t-reasm are much worse; indels added in cds, ..
## .. test 1. hi-kmer, 2. boost other vel quality options, 3. maybe drop -merge yes, skip lo-kmers?
## .. merge may force errors/indel reads into full.tr
## UPDATE: workd==rund==vel$oid on call; has oid.sam seq for this run...

if [ "X" = "X$oid" ]; then echo env oid=XXX; exit -1; fi

if [ "X" = "X$workd" ]; then 
  workd=`pwd`; # or fail?
  # workd=$HOME/scratch/chrs/daphmag/rnas/dmag5trfine/vel_$oid
fi

velpath=`which velvetg`
if [ "X" = "X$velpath" ]; then 
  velbin=$HOME/bio/velvet127s/bibinn2
  # velbin=/bio/bio-grid/mb/bin
  export PATH=$PATH:$velbin
fi
if [ "X" = "X$evigene" ]; then
  evigene=$HOME/bio/evigene/scripts
# evigene=/bio/bio-grid/mb/evigene/scripts
fi


## presume workd == trasm dir w/ pesam, trlong already inplace
cd $workd
#off# rund=$workd/vel$oid
#off# mkdir $rund

# NOTE: pesam should be name sorted: sort oidpe.samu > oidpe.sam
## FIXME: convert sam to fastq so velv uses QUAL scores..
sort $oid.pe.sam | cut -f1,10,11 | perl -ne'($id,$s,$q)=split; print "\@$id\n$s\n+\n$q\n";' > $oid.pe.fq
sort $oid.sr.sam | cut -f1,10,11 | perl -ne'($id,$s,$q)=split; print "\@$id\n$s\n+\n$q\n";' > $oid.sr.fq
pes=$oid.pe.fq
srs=$oid.sr.fq
trlong=$oid.tr

lopts="-ins_length 350"
vopts="-conserveLong yes $lopts" ## -scaffolding yes ?
# oopt1="-merge yes -min_pair_count 2 $lopts"
# kset1="31"
oopts="-scaffolding yes -min_pair_count 3 $lopts"
kset="97 87 77 65"

#off# cd $rund/

kseqdir=vel_seq
velveth $kseqdir 27 -fastq -shortPaired $pes -short $srs -fasta -long $trlong -noHash
# velveth $kseqdir 27 -sam -shortPaired $pesam -short $srsam -fasta -long $trlong -noHash

for k in $kset;  do { 
  ksubdir=vel_$k
  mkdir $ksubdir
  ln -s ../$kseqdir/Sequences $ksubdir/
  velveth $ksubdir $k  -reuse_Sequences
  velvetg $ksubdir $vopts -read_trkg yes 
  oases   $ksubdir $oopts 
  /bin/rm $ksubdir/{Graph2,LastGraph,PreGraph,Roadmaps,Sequences}

  $evigene/cdna_bestorf.pl -aa -cdna $ksubdir/transcripts.fa
  # use transcripts.aa to pick best kmer reasm; maybe stop kmer loop if got better asm?
  ## fixme: cdna_bestorf.pl : option for ngaps=NN in header
}
done

## basic test: compare each vel_kk/transcripts.aa aalen=nnn to $oid.tr aalen=nnn
## use transcripts.aa to pick best kmer reasm

## TEST: ../veltreasmfq.sh replace .sam w/ .fq : does vel use QUAL then? lots o bad TTTT/### seqs in sam vel
## YES, fastq is somewhat better; remaining problem: these mapped reads are not enough ..
##  full prior asm spans with NNN areas where this readset fails with bad seq: AAAAA.. or TTTT..
## **? can this be fixed using HALFMAPPING unmapped pair seq ?? all unmapped/crossmapped pair seqs?
##.. orig asm
# >dmag4vel4xcak45Loc2809t1 type=cdna; aalen=364,72%,complete; clen=1520;  strand=+; offs=330-1424;
# CCCGAAATAGCTGGCGGGCAACGAACAGCAGCAGATCTTGTCCCCCCCGCTTTTCACGTT
# ...
# GCAACCGATCTATTCTCTACAGTCGCAGCGCCAGACAAGAATATGAGAGGAGCTGCCCAT
# ATGGCAGATCATCTTCAAAGTTCATTGGGAATTTTGGGGTTTTTTTTTTTTTTTTTTNNN
# NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
# NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
# NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGGAGCAGCCGTTTATTTCGTGTTCTCTA
# AGGTTCACTCTGCATTGCGTTCTCCTATTCATTTTTCAGCAGAGACAAGCAGCAACATGC
##.. reasm, using sam > fastq
# qvel_97/transcripts.fa
# >Locus_1_Transcript_1/1_Confidence_1.000_Length_923
# AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA  << this is NNN above; missing prior gap-span seq
# AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGAGCAGCCGTTTATTTCGTGTTC  << GGAGCA.. as above
# TCTAAGGTTCACTCTGCATTGCGTTCTCCTATTCATTTTTCAGCAGAGACAAGCAGCAAC
#...............


