#! /bin/bash
### env libset=runsoaptr_locust2sol.libs qsub -q normal runsoaptr_locust.sh
#PBS -N soaptr7
#PBS -A ind114
#PBS -l nodes=1:ppn=16,walltime=23:55:00
#PBS -o soaptr7.$$.out
#PBS -e soaptr7.$$.err
#PBS -V

dv=1
ncpu=12

kset="31 25 21 45 39"
# F: gapfil; M: merge strength 1=def
# sopt="-p $ncpu -F -M 2 -t 99"
sopt="-p $ncpu -F -t 99"

# libset=runsoaptr_banana1.libs
oname=`basename $libset .libs | sed "s/runsoaptr_//; s/$/v$dv/;"`
rdconfig=cf$oname

bindir=$HOME/bio/soaptrans
datad=$HOME/scratchn
workd=$datad/chrs/aabugs/tsa/plants

# link in trsoap: velfad=$workd/sraf
rund=$workd/trsoap

cd $rund

cat $libset | grep -v '^#' > $rdconfig

for kmer in $kset; do {
 odir=sod$kmer
 outname=so${oname}.k$kmer
 mkdir $odir
if [ $kmer -lt 32 ]; then
  $bindir/SOAPdenovo-Trans-31kmer all -s $rdconfig -o $outname -K $kmer $sopt
else 
  $bindir/SOAPdenovo-Trans-127mer all -s $rdconfig -o $outname -K $kmer $sopt
fi

 mv $outname* $odir/
 if test -f $odir/$outname.contig ; then
  rm $odir/$outname.{preArc,vertex,edge.gz,Arc,updated.edge}
  gzip --fast $odir/$outname.{readOnContig,readInGap,ctg2Read,scafSeq,contig}
 fi

} done

