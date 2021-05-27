#! /bin/bash -l
### env estin=myest.fa qsub -q normal estgmap8.sh
#PBS -N estgmap8
#PBS -A ddp138
#PBS -l nodes=1:ppn=32,walltime=23:55:00
#PBS -o estgmap8.$$.out
#PBS -e estgmap8.$$.err
#PBS -V

# make loop thru all est/fasta/*.fa.gz files.
ncpu=32


datad=/oasis/$USER
workd=$datad/chrs/cacao
bindir=$datad/bio/bin
gmapd=$datad/bio/gmap118
rund=/scratch/$USER/$PBS_JOBID
gdb=gmap118

dgenome=cacao11allasm ; gtag=mars11
# dgenome=cirad_cacao1_chrs; gtag=cirad1c

# dgenome=cacao10allasm ; gtag=mars10
# dgenome=cacao9asm ; gtag=mars9
# dgenome=tcacao_cirad1asm; gtag=cirad
# dgenome=cacao41asm ; gtag=mars41


## chimera: -x 99 < many parts failed ? segfault error in gmap
gopt="-n 4 -S"; suf=out
# gopt="-n 4 -f 2"; suf=gff

notef=$workd/est/estgmap.$$.RUNNING
donef=`echo $notef | sed 's/RUNNING/DONE/'`

touch $notef
echo "START " >> $notef
echo `date`   >> $notef

mkdir -p $rund
cd $rund/

# many estin here?
# cp -p $workd/est/fasta/$estin $rund/
cp -rp $workd/est/fasta $rund/
cp -rp $workd/genome/$gdb  $rund/

gunzip fasta/*.gz

du -h >> $notef
ls -l >> $notef
ls -l fasta/ >> $notef

#... START LOOP estin : limit to suffix .fa, .fasta ..
for estin in fasta/*.{fa,fasta};  do { 
  if ! test -f "$estin" ; then continue; 
  fi
  
  dest=`basename $estin .fa | sed 's/.fasta//; s/.fa//;' `
  outf=$dest-$gtag.gmap.$suf
  echo "START gmap : $estin" >> $notef
  echo `date`  >> $notef
  
  i=0; while [ $i -lt $ncpu ]; do { 
  
   $gmapd/bin/gmap $gopt -D ./$gdb -d $dgenome --part=$i/$ncpu $estin > $dest.$gtag.gmap$i.$suf &
  
   i=$(( $i + 1 ))
  }
  done
  
  wait
  
  cat $dest.$gtag.gmap*.$suf > $outf
  gzip $outf
  ls -l $outf.gz >> $notef
  cp -p $outf.gz $workd/est/
  /bin/rm $dest.$gtag.gmap*.$suf

  echo "DONE gmap : $outf.gz" >> $notef
  du -h >> $notef

} done
#... END LOOP estin

echo "DONE " >> $notef
echo `date`  >> $notef
mv $notef $donef

#..................
