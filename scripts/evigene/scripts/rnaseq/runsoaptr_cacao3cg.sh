#! /bin/bash
### env libset=confsoap_uf.libs qsub -q normal runsoaptr_cacao3cg.sh
#PBS -N soaptr7
#PBS -A ind114
#PBS -l nodes=1:ppn=32,walltime=18:55:00
#PBS -o soaptr7.$$.out
#PBS -e soaptr7.$$.err
#PBS -V

## run on cg3rna 4-cultivar sets, 75bp, pe reads, inserts 100,180,292; 
##  uf: 105,493,626 reads
## make config.libs ; insert size by lane# 1,2,3: 5,6:  7,8:
# ls -1 fquf/uf*_1.fastq.gz | perl -ne'chomp; $l=$_; ($r=$l)=~s/_1.fa/_2.fa/; 
#  ($k,$in)=($l=~/l[123]_/)?(2,180):($l=~/l[78]_/)?(3,292):(1,100); 
#  print "[LIB]\nrank=$k\nasm_flag=3\navg_ins=$in\nq1=$l\nq2=$r\n\n"; BEGIN{ print "max_rd_len=80\n"; }'\
#   > confsoap_uf.libs
#...

## segfaulting on trestles .. why? try no -a mem; -p cpu <= 16

subd=trsoap3cg
# ncpu=32
ncpu=16
kset="31 69 59 49 27 23 39"
#.. try lower -M 2 ; 1 is default
sopt="-p $ncpu -F -t 29"

## confsoap_uf.libs
oname=`echo $libset | sed 's/confsoap_//; s/.libs//;'`
rdconfig=cf$oname

bindir=$HOME/bio/soaptrans
datad=$HOME/scratchn
rund=$datad/chrs/cacao/rnas/$subd
cd $rund

cat $libset | grep -v '^#' > $rdconfig

for kmer in $kset; do {
 odir=sod${subs}k$kmer
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
