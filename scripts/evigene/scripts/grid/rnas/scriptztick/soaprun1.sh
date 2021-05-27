#! /bin/bash
### qsub -q normal soaprun1.sh
#PBS -N soaptr7
#PBS -A ind114
#PBS -l nodes=1:ppn=16,walltime=23:55:00
#PBS -o soaptr7.$$.out
#PBS -e soaptr7.$$.err
#PBS -V

ncpu=6

## soap needs p=pair.fa or q1=read1.fq; q2=read2.fq
## ztick.v1: 64GB reads, soaptr fails (no error) -- outamem 64gb presumably

dv=1; 
kset="31 25 29 85 75 65 55 45 39"
# F: gapfil; M: merge strength 1=def
sopt="-p $ncpu -F -t 99"

bindir=$HOME/bio/soaptrans
datad=$HOME/scratchg
workd=$datad/chrs/aabugs/tsa/ztick
# inpe=$workd/sraf/split/SRR*_1.fastq
inpe=$workd/sraf/allpe_1.fq

rund=trsoap$dv
oname="ztick$dv"

cd $workd/
mkdir $rund

echo "max_rd_len=105" > $rund/rdconfig 
for fq1 in $inpe; do {
  fq2=`echo $fq1 | sed 's/_1/_2/;'`
  cat >> $rund/rdconfig <<EOT

[LIB]
rank=1
asm_flag=3
avg_ins=300
q1=$fq1
q2=$fq2

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

