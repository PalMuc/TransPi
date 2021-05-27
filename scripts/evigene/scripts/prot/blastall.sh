#! /bin/bash
### env protin=prot.aa qsub -q normal blastall.sh
#PBS -N blastall 
#PBS -A ind114
#PBS -l nodes=6:ppn=16,walltime=26:55:00
#PBS -o blastall.$$.out
#PBS -e blastall.$$.err
#PBS -V

ncpu=96

nbin=$HOME/bio/ncbi/bin
evigene=$HOME/bio/evigene

#gordo# 
datad=$HOME/scratchg/
#tres# datad=$HOME/scratchn/
workd=$datad/chrs/nasv1
rund=$workd/prot

# prodb=hym11set12
prodb=$protin

pronam=`basename $prodb .aa`
odir=bloutg
bltag=sd; blopt="-evalue 1e-5"

cd $workd/prot/
if [ ! -f $prodb.psq ]; then $nbin/makeblastdb -dbtype prot -in $prodb -logfile $prodb.mbl.log; fi
mkdir $odir

if [ ! -f $protin.split.1.fa ]; then
 pindir=`dirname $protin`
 splitsize=`grep -v '^>' $protin | wc -c | sed 's/ .*//' `
 splitbp=$(( $splitsize / $ncpu ))
 $evigene/scripts/splitMfasta.pl --outputpath=$pindir --minsize=$splitbp $protin
fi

qset=`/bin/ls $protin.split.*.fa`

for qfile in $qset
{
  qnam=`basename $qfile .fa`
  onam=$odir/$bltag-$pronam-$qnam
  echo $nbin/blastp $blopt -outfmt 7 -db $prodb -query $qfile -out $onam.blastp
  $nbin/blastp $blopt -outfmt 7 -db $prodb -query $qfile -out $onam.blastp  &
}

wait

opack=`echo $bltag-$pronam-$qnam | sed 's/.split.*//'`
cat $odir/$opack.*.blastp > $opack.blastp
gzip --fast $opack.blastp


## .. blastall ran out of 19h time at 210k, 2/3 way thru  340K prots, on 96 cpu
## .. copy out last 1/3 from aa12b/hym*.split.nn.fa and rerun blastall
# hym12bsplitredo.tab from hym12baasplit.count - hym12bsplit17rblp.count (17hr of 18hr point)
# hym11set12b.	10.blastp	1897	10.fa	3179	1282
# hym11set12b.	11.blastp	1889	11.fa	3092	1203
# 
# cat hym12bsplitredo.tab | perl -ne \
# '($nfp,$ib,$nb,$ia,$na,$dab)=split; 
# open(F,"aa12b/${nfp}aa.split.$ia") or die"aa12b/$nfp$ia";
# open(O,">aa12re/${nfp}aa.split.$ia") or die "aa12re/$nfp$ia";
# $oa=$ja=0; while(<F>) { if(/^>(\S+)/) { $id=$1; $ja++; $ok=($ja>$nb)?1:0; $oa++ if($ok); } print O $_ if($ok); } 
# close(O); close(F);warn "copied $oa/$dab to aa12re/${nfp}aa.split.$ia \n"; ' \
