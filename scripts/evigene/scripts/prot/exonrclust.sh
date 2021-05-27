#! /bin/bash
### env exontab=nasvit1-protexonregion.tab qsub -q normal exonrclust.sh
#PBS -N exonr1
#PBS -l nodes=1:ppn=32,walltime=13:55:00
#PBS -o exonr1.$$.out
#PBS -e exonr1.$$.err
#PBS -V

ncpu=32

exbin=$HOME/bio/exonerate/bin
export PATH=${exbin}:${PATH}
workd=$HOME/scratch/chrs/nasv1
dgenome=nasvit1asm

# xtype=prot
xtype=bestprot

cd $workd/prot/

if ! test -f $exontab ; then  echo "missing input $exontab"; exit;  fi
nam=`echo $exontab | sed 's/\..*//' `
echo "# start exonrclust: $nam `date`"

i=0 
while [ $i != $ncpu ]; do {
  out=$nam.x$xtype$i.gff ;  touch $out

  echo "# perl ./exonrclust.pl -exonerate $xtype -in=$exontab -out=$out -i=$i -n=$ncpu "
  #Fail2# perl exonrclust.pl -debug -in=$exontab -out=$out -i=$i -n=$ncpu 
  ( cat $exontab | perl exonrclust.pl -exonerate $xtype -debug -in=stdin -out=$out -i=$i -n=$ncpu ) &

  i=$(( $i + 1 ))
}
done

wait
echo "# done exonrclust: `date`"
