#! /bin/bash
### env name=ztickfem inpe=sraf/SRR*_1.fasta datad=`pwd` qsub -q normal soaprun1.sh
#PBS -N soaptr7
#PBS -A ind114
#PBS -l nodes=1:ppn=16,walltime=39:55:00
#PBS -V

ncpu=8
#which? if [ "X" = "X$ncpu" ]; then ncpu=6; fi

## soap needs p=pair.fa or q1=read1.fq; q2=read2.fq
## FIXME: no fa.gz, ungzip inpe ..

#bad-utr big tr: dv=4b; INSIZE=450;
# dv=4c; INSIZE=250; MAXRD=125
# dv=4d; INSIZE=200; MAXRD=105
dv=4m; INSIZE=250; MAXRD=152

#v3: kset="31 25 75 65 55 45"
#v4: best honbee k=27 31 45 55; but in evg3 final, k25 commonest; 
#v4b/c: kset="31 27 25 35 45 55 65"
#ve
kset="31 27 35 43 53 63 73 83 93 103"
## for soap ? *outer* insertsize, including read sizes, or not? : 450 is bad (250 inner + 2*100 rdlen)

# F: gapfil; M: merge strength 1=def; R: RPKM.stats
sopt="-p $ncpu -R -F -t 99"

## UPD: FIX 31kmer >> 31mer
bindir=$HOME/bio/soaptrans103

if [ "X" = "X$datad" ]; then echo "err missing datad=path/to/data"; exit -1; fi
if [ "X" = "X$inpe" ]; then echo "err missing inpe=sraf/SRR*_1.fasta"; exit -1; fi
# inpe=`echo $inpe | sed "s,sraf,$datad/sraf,g;"`

## fixme: name needs to be in rund or else rdconfig$name ..
if [ "X" = "X$name" ]; then name=ztick; fi
oname="$name$dv"
rund=trsoap$name$dv

cd $datad/
mkdir $rund

## inpe presumed local path rel to datad
echo "max_rd_len=$MAXRD" > $rund/rdconfig 
for fa1 in $inpe; do {
  fa2=`echo $fa1 | sed 's/_1/_2/;'`
  if [ $fa2 = $fa1 ]; then continue; fi

  nogz=`echo $fa1 | sed 's/\.gz//;'`
  if [ $nogz != $fa1 ]; then
    fa11=`basename $fa1 .gz`; fa11=$rund/$fa11;
    fa21=`basename $fa2 .gz`; fa21=$rund/$fa21;
    gunzip -c $fa1 > $fa11
    gunzip -c $fa2 > $fa21
    fa1=$fa11
    fa2=$fa21
  fi

  cat >> $rund/rdconfig <<EOT

[LIB]
rank=1
asm_flag=3
avg_ins=$INSIZE
f1=$datad/$fa1
f2=$datad/$fa2

EOT

} done

cd $rund

for kmer in $kset; do {
 odir=sod$kmer
 outname=so${oname}.k$kmer
 mkdir $odir
if [ $kmer -lt 32 ]; then
  $bindir/SOAPdenovo-Trans-31mer all -s rdconfig -o $outname -K $kmer $sopt
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

