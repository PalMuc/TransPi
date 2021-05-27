#! /bin/bash
### env trdb=pathtogmapdb readset=fastq/xxx*_1.fastq.gz rnam=yyy rund=path/to  qsub -q normal gsnap2genes.sh
#PBS -N gsnap2genes
#PBS -l nodes=2:ppn=16,walltime=39:55:00
#PBS -V

## new params: readset fails/last only ; need readseq="fq/xx*1.fq fq/yy*1.fq" not {curlies}
## gsnn2genes.sh: reset opts for notredam set: fastq5/nodamx1_1.fastq.gz ..
#  env readset=fq4x/{Dman_42,Dman_64,Dman_65}*_1.fastq.gz rnam=xfi trdb=trdb/xxx rund=`pwd` qsub ..
## takes 12h .. 15+hr for one readset, 48 cpu

ncpu=32; npart=32; jstart=0

if [ "X" = "X$trdb" ]; then echo "ERR: miss trdb=pathtogmapdb"; exit -1; fi
if [ "X" = "X$readset" ]; then echo "ERR: miss readset=fqs/xxx*_1.fastq.gz"; exit -1; fi
if [ "X" = "X$rund" ]; then echo "ERR: miss rund=pathto/data"; exit -1; fi

# keepset="concordant_mult,concordant_uniq,halfmapping_uniq,paired_mult,paired_uniq_long"
# poormap="halfmapping_mult,paired_uniq_inv,paired_uniq_scr,unpaired_mult,unpaired_uniq,nomapping"
# skipset="concordant_transloc,halfmapping_transloc,unpaired_transloc"  # these are empty probably
#-----------

bindir=$HOME/bio/bin
gmapd=$HOME/bio/gmap1309

## new opt 1303:  --maxsearch=1000 >= npaths ; smaller == speedup 
## new: SEGfaults .. no clue why..  NO GOOD cant get rid of segfaults

gdbpath=`dirname $trdb`
dgenome=`basename $trdb`
gtag=`echo $dgenome | sed 's/_.*//;'`

## add sam header controls:
#  --no-sam-headers               Do not print headers beginning with '@'
#  --sam-headers-batch=INT        Print headers only for this batch, as specified by -q

## defaults: px=200, pd=100
# optreads="--pairexpect=450 --pairdev=50"
optreads=""
## genome -N 1; transcript -N 0
## optgenome="-N 1  " 
## transcript
optgenome="-N 0 --gmap-mode=none " 
optsnp=""
snapopt1="--sam-headers-batch=0 $optgenome $optsnp $optreads" 

cd $rund

readfiles=`ls $readset`
if [ "X" = "X$rnam" ]; then
  reada=($readset)
  rnam=`dirname ${reada[0]}`
fi

outdir=gso$rnam
mkdir $outdir

for rnain in $readfiles ;  do 
{ 
  cd $rund/
  if ! test -f "$rnain" ; then continue; fi
  
  drna=`basename $rnain .gz | sed 's/\.fastq//; s/\.fq//; '`
  outna=$drna-$gtag
  snapopts="--gunzip $snapopt1"
  
  rnain2=`echo $rnain | sed 's/1.fq/2.fq/; s/1.fastq/2.fastq/; '`
  if test -f "$rnain2" ; then
    inset="$rnain $rnain2"
  else 
    echo "ERR: missing pair file2 $rnain2"; continue;
  fi

  i=0;  j=$jstart;

  echo "CMD: $gmapd/bin/gsnap $snapopts -A sam -D $gdbpath -d $dgenome $inset TO $outdir/$drna.gsnap0.samu "
  ##* change output to --split-output=$outdir/$drna.gsnap$i

  while [ $i -lt $ncpu ]; do { 
  
   $gmapd/bin/gsnap $snapopts  -A sam  --part=$j/$npart \
    -D $gdbpath -d $dgenome --split-output=$outdir/$drna.gsnap$i $inset &

   i=$(( $i + 1 ))
   j=$(( $j + 1 ))
  }
  done
  wait

} 
done

#... END LOOP rnain