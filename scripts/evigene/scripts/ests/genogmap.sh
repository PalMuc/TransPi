#! /bin/bash
### env trset=fastdir/*.tr  qsub -q normal genogmap.sh
#PBS -N genogmap
#PBS -A ind114
#PBS -l nodes=1:ppn=32,walltime=11:55:00
#PBS -o genogmap.$$.out
#PBS -e genogmap.$$.err
#PBS -V

ncpu=32
bindir=$HOME/bio/gmap1307/bin

# shoud be params:
datad=$HOME/scratchn/chrs/kfish
subd=$datad/genes
gdb=$datad/genome/gmap12
dgenome=killifish20130322asm; gtag=kfish2b

outdir=gmoutf

#for trmap: gopt="--nosplicing -n 9 -S"; suf=outns
#for gff, not good# gopt="-n 4 -f 2"; suf=gff
#12# gopt="--suboptimal-score=2 --min-intronlength=29 -n 9 -S"; suf=out9
#13a# gopt="--suboptimal-score=2 --min-intronlength=29 -n 9 -S"; suf=out137
## .. use this one, def minintr, but keep subopt
gopt="--suboptimal-score=2 -n 9 -S"; suf=out13

if [ "X$trset" = "X" ]; then echo "ERR: missing trset="; exit -1; fi

cd $subd/
echo "START `date`" 
mkdir $outdir

#... START LOOP estin
for estin in $trset;  do { 
  if ! test -f "$estin" ; then continue; fi

  nogz=`echo $estin | sed 's/\.gz//;'`
  dest=`basename $estin .gz | sed 's/\.tr//; s/\.fasta//; s/\.fa//;' `
  outf=$dest-$gtag.gmap.$suf
  echo "START gmap : $estin" >> $notef
  echo `date`  >> $notef
  ##gopt1=$gopt; ##NO good for gmap# if [ $estin != $nogz ]; then gopt1="--gunzip $gopt"; fi
 
  i=0; while [ $i -lt $ncpu ]; do { 

  if [ $estin != $nogz ]; then
    gunzip -c $estin | $bindir/gmap $gopt -D $gdb -d $dgenome --part=$i/$ncpu > $outdir/$dest.$gtag.gmap$i.$suf &
  else
    $bindir/gmap $gopt -D $gdb -d $dgenome --part=$i/$ncpu $estin > $outdir/$dest.$gtag.gmap$i.$suf & 
  fi

  i=$(( $i + 1 ))
  } done

  wait
  
  cat $outdir/$dest.$gtag.gmap*.$suf > $outf
  gzip --fast $outf

} done
#... END LOOP estin
echo "DONE `date` "
