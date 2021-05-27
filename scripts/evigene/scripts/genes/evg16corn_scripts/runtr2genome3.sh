#! /bin/bash
# runtr2genome2b.sh
### env cdsseq=mycds.fa chrseq=mychrasm.fa datad=`pwd` qsub -q normal runtr2genome.sh
#PBS -N tr2genome2b
#PBS -l nodes=1:ppn=16,walltime=18:55:00
#PBS -V

## runtr2genome3 = update to cdsqual.pl, asmrna_dupfilter3.pl ?? Code/Noncode .. other
## runtr2genome2b = updated from mosquito/anoph/evg2anofunzh/merge24f,pub24set work
## revised: asmrna_dupfilter2b.pl, trclass2mainalt.pl, overeqcdsloc.pl
#
# $evigene/scripts/prot/overeqcdsloc.pl  -EQXONSHOW -alts -strand -minover 10 -debug \
#  -in evg24mergeanofun.btall -hotab evg24mergeanofun.aa.btall -out evg24mergeanofun.eqgene5
#
# $evigene/scripts/rnaseq/asmrna_dupfilter2b.pl -CDSALIGN -debug \
# -ablastab evg24mergeanofun.aa.btall  -aasize evg24mergeanofun.cds.qual -eqgene evg24mergeanofun.eqgene5 \
#  -outeqtab $pt.tgalntab -outclass $pt.tgclass5h
#
# $evigene/scripts/prot/trclass2mainalt.pl -cullx  -debug -idpre Anofunz4kEVm -trclass evg24m2banofun.tgclass5h
#

## test shell for evigene tr2genome.pl  
#----------------------------

vtag=evgm3
CDSBLAST_IDENT=98
CHRBLAST_IDENT=95
CDSBLAST_EVALUE=1e-19

evigene=$HOME/bio/evigene
export PATH=$HOME/bio/ncbi2230/bin:$PATH
# nbin=$HOME/bio/ncbi2230/bin

if [ "X" = "X$cdsseq" ]; then echo "ERR: missing cdsseq=pathto/cds.fa"; exit -1; fi
if [ "X" = "X$chrseq" ]; then echo "ERR: miss chrseq=pathto/chrasm.fa"; exit -1; fi
if [ "X" = "X$datad" ]; then echo "ERR: miss localdir datad=path/to/data"; exit -1; fi
if [ "X" = "X$gtag" ]; then gtag=`basename $chrseq | sed 's/\..*$//;'`; fi
if [ "X" = "X$ncpu" ]; then ncpu=16; fi
# if [ "X" != "X$aablast" ]; then echo see below; fi

## also need original cds-self blast for asmrna_dupfilter, or regenerate??
# cdsselfblast=
chrdb=`basename $chrseq | sed 's/\.fa.*//;'`
trname=`basename $cdsseq .cds | sed 's/\.cds//; s/\.fa.*//;'`

outdir=${vtag}f
onamp=$vtag$gtag
oname=$onamp-$trname

cd $datad/
mkdir $outdir

echo "# START tr2genome $cdsseq `date`" 

cdssize=$cdsseq.qual
if [ ! -f $cdssize ]; then
  # env ismrna=1 oid=1 off=1 $evigene/scripts/prot/aaqual.sh $cdsseq
  env ismrna=1 oid=1 off=1 $evigene/scripts/prot/cdsqual.pl $cdsseq
fi

# partition query.cds by ncpu
## FIXME: revise splitMfasta.pl to do all this: splitfasta.pl --nparts $ncpu
if [ ! -f $cdsseq.split.1.fa ]; then
 pindir=`dirname $cdsseq`
 splitsize=`grep -v '^>' $cdsseq | wc -c | sed 's/ .*//' `
 splitbp=$(( $splitsize / $ncpu ))
 $evigene/scripts/splitMfasta.pl --outputpath=$pindir --minsize=$splitbp $cdsseq
fi

# 6.1: blastn -query evgene.cds -db chrasm 
#----------------------------
blopt="-task megablast -perc_identity $CHRBLAST_IDENT -evalue $CDSBLAST_EVALUE -outfmt 7 "; opref="megi0";
echo "#t2g: makeblastdb -dbtype nucl -in $chrseq -out $chrdb"; 
makeblastdb -dbtype nucl -in $chrseq -out $chrdb;  

qset=`/bin/ls $cdsseq.split.*.fa`
i=0; for qfile in $qset
{
  qnam=`basename $qfile .fa`
  onam=$onamp-$qnam

  if [ $i -eq 1 ]; then 
  echo "#t2g.$i: blastn $blopt -db $chrdb  -query $qfile -out $outdir/$onam.blastn"; 
  fi
  blastn $blopt -db $chrdb  -query $qfile -out $outdir/$onam.blastn &
  i=$(( $i + 1 ))
}
wait

cdsblast=$oname.blastn
cat $outdir/$onamp-$trname*.blastn > $cdsblast

# 6.2: blastn to btall table
#----------------------------
# ? change -pmin=0.85 to -pmin=0.95 ? ie stringent for 2nd dup aligns
cdsbltab=$oname.btall
mbopts="-tall -spans=2 -onegenome=g1 -pctover=0.03 -pmin=0.95";
echo "#t2g: $evigene/scripts/makeblastscore3.pl  $mbopts -aasize $cdssize $cdsblast > $cdsbltab";
$evigene/scripts/makeblastscore3.pl  $mbopts -aasize $cdssize $cdsblast > $cdsbltab;

# 6.2b: aablast
#----------------------------
## UPDATE: add (require?) blastp evigene_pubset.aa x refprots, corresponding toinput cdsseq
## input may be refprots.aa, add logic to run blastp? convert to format for asmrna_dupfilter -ablastab or -anames
## use ref homol scores to retain good genes when filtering of junk loci/alts (ie many small, 1exon things)
##.. leave to caller for now..
#n if [ "X" != "X$aablast" ]; then
#n   if [ -f $aablast ]; then
#n     ##  need query,ref.aa.size table? or option -noNEEDLEN [default]
#n     aablbase=`basename $aablast | sed 's/.gz//;'`; 
#n     aabltab=`echo $aablbase | sed 's/\.blastp//;'`; 
#n     if [ $aablbase != $aabltab ]; then
#n       aabltab=$aabltab.btall
#n       mbaopts="-tall";
#n       if [ "X" != "X$aasize" ]; then  mbaopts="$mbaopts -aasize $aasize"; fi
#n       echo "#t2g: $evigene/scripts/makeblastscore3.pl $mbaopts $aablast TO $aabltab";
#n       $evigene/scripts/makeblastscore3.pl  $mbaopts $aablast > $aabltab;
#n       aablast=$aabltab
#n     fi
#n   fi
#n     
#n   if [ ! -s $aablast ]; then 
#n     echo "#WARN: missing -aablast=$aablast"; # continue
#n     aablast=""
#n   fi
#n fi

# 6.3 eqgene table
#----------------------------
cdseqgene=$oname.eqgene
#upd: 
# $evigene/scripts/prot/overeqcdsloc.pl  -EQXONSHOW -alts -strand -minover 10 -debug \
#  -hotab evg24mergeanofun.aa.btall -in evg24mergeanofun.btall  -out evg24mergeanofun.eqgene5

eqopts="-EQXONSHOW -alts -strand -minover 10 -debug";  
if [ "X" != "X$aablast" ]; then
  eqopts="$eqopts -hotab $aablast"
fi
echo "#t2g:" $evigene/scripts/prot/overeqcdsloc.pl $eqopts -in $cdsbltab -out $cdseqgene;
$evigene/scripts/prot/overeqcdsloc.pl $eqopts -in $cdsbltab -out $cdseqgene;


# 6.4 trclass
#----------------------------

## also need cds-self blast for asmrna_dupfilter, regenerate rather than reuse orig evigene set
# 6.4a: blastn -query evgene.cds -db evgcds 
oself=${vtag}self
cdsdb=$trname
selfblopt="-task megablast -ungapped -xdrop_ungap 4 -dust no -perc_identity $CDSBLAST_IDENT -evalue $CDSBLAST_EVALUE -outfmt 7 "; 

echo "#t2g: " makeblastdb -dbtype nucl -in $cdsseq -out $cdsdb; 
makeblastdb -dbtype nucl -in $cdsseq -out $cdsdb;  

qset=`/bin/ls $cdsseq.split.*.fa`
i=0; for qfile in $qset
{
  qnam=`basename $qfile .fa`
  onam=$oself-$qnam
  if [ $i -eq 1 ]; then 
  echo "#t2g.$i: blastn $selfblopt -db $cdsdb  -query $qfile -out $outdir/$onam.blastn"; 
  fi
  blastn $selfblopt -db $cdsdb  -query $qfile -out $outdir/$onam.blastn &
  i=$(( $i + 1 ))
}
wait

cdsselfblast=$oself-$trname.blastn
cat $outdir/$oself-$trname*.blastn > $cdsselfblast

# 6.4b: tr classing with  asmrna_dupfilter2b
#----------------------------
#upd:
# $evigene/scripts/rnaseq/asmrna_dupfilter2b.pl -CDSALIGN -debug \
# -ablastab evg24mergeanofun.aa.btall  -aasize evg24mergeanofun.cds.qual -eqgene evg24mergeanofun.eqgene5 \
#  -outeqtab $pt.tgalntab -outclass $pt.tgclass5h

outaln=$oname.tgalntab
outclass=$oname.tgclass
clopts="-CDSALIGN -debug";
# clopts="$clopts -acdhit $aacdhit" ; # option skip for now
if [ "X" != "X$aablast" ]; then
  clopts="$clopts -ablastab $aablast"
fi

echo "#t2g: " $evigene/scripts/rnaseq/asmrna_dupfilter3.pl $clopts -aasize $cdssize \
  -eqgene $cdseqgene -blastab $cdsselfblast -outeqtab $outaln -outclass $outclass ;
$evigene/scripts/rnaseq/asmrna_dupfilter3.pl $clopts -aasize $cdssize \
  -eqgene $cdseqgene -blastab $cdsselfblast -outeqtab $outaln -outclass $outclass ;


# 6.5: add make trclass to pubid table 
#----------------------------
#upd:
# $evigene/scripts/prot/trclass2mainalt.pl -cullx  -debug -idpre Anofunz4kEVm -trclass evg24m2banofun.tgclass5h

maopt="-cullx -debug"; 
if [ "X" != "X$idpre" ]; then maopt="$maopt -idpre $idpre"; fi
echo "#t2g: " $evigene/scripts/prot/trclass2mainalt.pl $maopt -trclass  $outclass ;
$evigene/scripts/prot/trclass2mainalt.pl $maopt -trclass  $outclass ;

# tidy up....
#----------------------------
echo "#t2g: tidyup tmpfiles"
rm $cdsseq.split.*.fa  
rm $cdsdb.{nsq,nin,nhr}
rm $chrdb.{nsq,nin,nhr}
#later: rm -f $outdir

echo "#t2g: gzip tmpfiles"
gzip --fast $cdsblast &
gzip --fast $cdsselfblast &
gzip --fast $cdsbltab &
gzip --fast $cdseqgene &
gzip --fast $outaln &

wait

echo "# DONE `date`" 
#----------------------------

#----------------------------
