#! /bin/bash
### qsub -q normal soaprun1.sh
#PBS -N soaptr7
#PBS -A ind114
#PBS -l nodes=1:ppn=16,walltime=23:55:00
#PBS -o soaptr7.$$.out
#PBS -e soaptr7.$$.err
#PBS -V

ncpu=6

## kf_mdibl_rnaseq_reads_1.fq.gz  kf_whoi_rnaseq_reads_1.fq.gz
## soap needs p=pair.fa or q1=read1.fq; q2=read2.fq
# dv=1m; inpe1=fastq/kf_mdibl_rnaseq_reads_1.fq.gz
dv=1w; inpe1=fastq/kf_whoi_rnaseq_reads_1.fq.gz

inpe2=`echo $inpe1 | sed 's/1\./2./;'`;

## fail k75 : k too big for soap??
#kset1="31 25 29 21 75 65 55 45 39"
kset="31 25 29 21 37 45"
## F: gapfil; M: merge strength 1=def
sopt="-p $ncpu -F -t 99"

bindir=$HOME/bio/soaptrans
datad=$HOME/scratchg
workd=$datad/chrs/kfish/rnas
rund=trsoap$dv
oname="kf$dv"

cd $workd/
mkdir $rund

## check for existing rund/pereads.fa ; symlink doesnt work right for soap ??
if [ ! -f $rund/allpe_1.fq ]; then
  if [ -f $inpe1 ]; then 
  gunzip -c $inpe1 > $rund/allpe_1.fq; 
  gunzip -c $inpe2 > $rund/allpe_2.fq;
  fi
fi
if [ ! -f $rund/allpe_1.fq ]; then
 echo "missing $rund/allpe_1.fq "; exit -1
fi

cat > $rund/rdconfig  <<EOT
max_rd_len=110
[LIB]
rank=1
asm_flag=3
avg_ins=380
q1=allpe_1.fq
q2=allpe_2.fq
EOT
# p=pereads.fa

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
  ## drop also .readXXX
  rm $odir/$outname.{preArc,vertex,edge.gz,Arc,updated.edge}
  gzip --fast $odir/$outname.{readOnContig,readInGap,ctg2Read,scafSeq,contig}
 fi

} done

