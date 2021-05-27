#! /bin/bash
### env name=ztickfem inpe=sraf/split/SRR*_1.fastq.gz qsub -q normal soaprun1.sh
#PBS -N soaptr7
#PBS -A ind114
#PBS -l nodes=1:ppn=16,walltime=23:55:00
#PBS -o soaptr7.$$.out
#PBS -e soaptr7.$$.err
#PBS -V

ncpu=6

## soap needs p=pair.fa or q1=read1.fq; q2=read2.fq
## ztick.v1: 64GB reads, soaptr fails (no error) -- outamem 64gb presumably
## v2: switch to fasta filtered; separate fem/male inpe= sets

dv=3; 
kset="31 25 29 85 75 65 55 45 39"
# F: gapfil; M: merge strength 1=def
sopt="-p $ncpu -F -t 99"

bindir=$HOME/bio/soaptrans
datad=$HOME/scratchg
workd=$datad/chrs/aabugs/tsa/ztick

# inpe=$workd/sraf/split/SRR*_1.fastq
if [ "X" = "X$inpe" ]; then inpe=sraf/allpe_1.fa; fi

## fixme: name needs to be in rund or else rdconfig$name ..
if [ "X" = "X$name" ]; then name=ztick; fi
oname="$name$dv"
rund=trsoap$name$dv

cd $workd/
mkdir $rund

## inpe presumed local path rel to workd
echo "max_rd_len=105" > $rund/rdconfig 
for fa1 in $inpe; do {
  fa2=`echo $fa1 | sed 's/_1/_2/;'`
  cat >> $rund/rdconfig <<EOT

[LIB]
rank=1
asm_flag=3
avg_ins=300
f1=$workd/$fa1
f2=$workd/$fa2

EOT

} done

cd $rund

for kmer in $kset; do {
 odir=sod$kmer
 outname=so${oname}.k$kmer
 mkdir $odir
if [ $kmer -lt 32 ]; then
  $bindir/SOAPdenovo-Trans-31kmer all -s rdconfig -o $outname -K $kmer $sopt
else 
  $bindir/SOAPdenovo-Trans-127mer all -s rdconfig -o $outname -K $kmer $sopt
fi

 mv $outname* $odir/
 if test -f $odir/$outname.contig ; then
  rm $odir/$outname.{readOnContig,readInGap,ctg2Read,preArc,vertex,edge.gz,Arc,updated.edge}
  gzip --fast $odir/$outname.{scafSeq,contig}
  ## drop read set...  gzip --fast $odir/$outname.{readOnContig,readInGap,ctg2Read,scafSeq,contig}
 fi

} done

