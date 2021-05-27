#! /bin/bash
# runtr2genome.sh
### env cdsseq=mycds.fa chrseq=mychrasm.fa datad=`pwd` qsub -q normal runtr2genome.sh
#PBS -N tr2genome
#PBS -l nodes=1:ppn=16,walltime=18:55:00
#PBS -V

## test shell for evigene tr2genome.pl ; OK now
#   # 6.1. blastn okayset/my.cds to chrasm
#   my($cdsgmapblast)= blast2genome($cdsseq, $chrasm);
#   # 6.2. tabulate blast > cdsgmap.tall > cdsgmap.equalgene
#   my @xx= cdsgmap_maketables($cdsgmapblast);
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
  env ismrna=1 oid=1 off=1 $evigene/scripts/prot/aaqual.sh $cdsseq
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
## or leave to caller, expect aablast option to be == aabltab format
## input may be refprots.aa, add here logic to run blastp? convert to input for asmrna_dupfilter -ablastab or -anames (prefered)
## use ref homol scores to retain good genes when filtering of junk loci/alts (ie many small, 1exon things)
## leave off for now, expect user to provide makeblastscore table
# if [ "X" != "X$aablast" ]; then
#   if [ -f $aablast ]; then
#     ##  need query,ref.aa.size table? or option -noNEEDLEN [default]
#     aablbase=`basename $aablast .gz`; 
#     aabltab=`basename $aablbase .blastp`; 
#     if [ $aablbase != $aabltab ]; then
#       aabltab=$aabltab.btall
#       mbaopts="-tall";
#       if [ "X" != "X$aasize" ]; then  mbaopts="$mbaopts -aasize $aasize"; fi
#       echo "#t2g: $evigene/scripts/makeblastscore3.pl $mbaopts $aablast > $aabltab";
#       $evigene/scripts/makeblastscore3.pl  $mbaopts $aablast > $aabltab;
#       aablast=$aabltab
#     fi
#   fi
#     
#   if [ ! -s $aablast ]; then 
#     echo "#WARN: missing -aablast=$aablast"; # continue
#     aablast=""
#   fi
# fi

# 6.3 eqgene table
#----------------------------
cdseqgene=$oname.eqgene
# fixme: here or overeqcdsloc, -maxlist 10 too little, want all? overlap ids (noalts, minover)
eqopts="-noalts -strand -minover 10 -debug";  
echo "#t2g: $evigene/scripts/prot/overeqcdsloc.pl $eqopts -in $cdsbltab -out $cdseqgene";
$evigene/scripts/prot/overeqcdsloc.pl $eqopts -in $cdsbltab -out $cdseqgene;


# 6.4 trclass
#----------------------------

## also need cds-self blast for asmrna_dupfilter, regenerate rather than reuse orig evigene set
# 6.4a: blastn -query evgene.cds -db evgcds 
oself=${vtag}self
cdsdb=$trname
selfblopt="-task megablast -ungapped -xdrop_ungap 4 -dust no -perc_identity $CDSBLAST_IDENT -evalue $CDSBLAST_EVALUE -outfmt 7 "; 

echo "#t2g: makeblastdb -dbtype nucl -in $cdsseq -out $cdsdb"; 
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
outaln=$oname.tgalntab
outclass=$oname.tgclass
clopts="-CDSALIGN -debug";
# clopts="$clopts -acdhit $aacdhit" ; # option skip for now
if [ "X" != "X$aablast" ]; then
  clopts="$clopts -ablastab $aablast"
fi

echo "#t2g: $evigene/scripts/rnaseq/asmrna_dupfilter2b.pl $clopts -aasize $cdssize \
  -eqgene $cdseqgene -blastab $cdsselfblast -outeqtab $outaln -outclass $outclass" ;
$evigene/scripts/rnaseq/asmrna_dupfilter2b.pl $clopts -aasize $cdssize \
  -eqgene $cdseqgene -blastab $cdsselfblast -outeqtab $outaln -outclass $outclass ;

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

# test3 result:
# /bio/bio-grid/aabugs4/bugs/mosquito/anoph/evg2anofunzg
# evg2anofunz4g.cds@                          evgm2cchranofun-evg2anofunz4g.tgalntab.gz
# evg2anofunz4g.cds.qual                      evgm2cchranofun-evg2anofunz4g.tgclass
# evgm2cchranofun-evg2anofunz4g.blastn.gz     evgm2cself-evg2anofunz4g.blastn.gz
# evgm2cchranofun-evg2anofunz4g.btall         srun_comet.sh*
# evgm2cchranofun-evg2anofunz4g.eqgene        tr2chr.1624410.comet-26-72.out
# geno3f/tr2chr.1624410.comet-26-72.out | grep '^#t2g'
# START tr2genome evg2anofunz4g.cds Sun Feb 14 19:27:40 PST 2016
# DONE Sun Feb 14 19:55:12 PST 2016
#t2g: makeblastdb -dbtype nucl -in genoasm_AfunF1.fa -out genoasm_AfunF1
#t2g.1: blastn -task megablast -perc_identity 95 -evalue 1e-19 -outfmt 7  -db genoasm_AfunF1  -query evg2anofunz4g.cds.split.10.fa -out evgm2cf/evgm2cchranofun-evg2anofunz4g.cds.split.10.blastn
#t2g: /home/ux455375/bio/evigene/scripts/makeblastscore3.pl  -tall -spans=2 -onegenome=g1 -pctover=0.03 -pmin=0.85 -aa evg2anofunz4g.cds.qual evgm2cchranofun-evg2anofunz4g.blastn > evgm2cchranofun-evg2anofunz4g.btall
#t2g: /home/ux455375/bio/evigene/scripts/prot/overeqcdsloc.pl -noalts -strand -minover 10 -debug -in evgm2cchranofun-evg2anofunz4g.btall -out evgm2cchranofun-evg2anofunz4g.eqgene
#t2g: makeblastdb -dbtype nucl -in evg2anofunz4g.cds -out evg2anofunz4g
#t2g.1: blastn -task megablast -ungapped -xdrop_ungap 4 -dust no -perc_identity 98 -evalue 1e-19 -outfmt 7  -db evg2anofunz4g  -query evg2anofunz4g.cds.split.10.fa -out evgm2cf/evgm2cself-evg2anofunz4g.cds.split.10.blastn
#t2g: /home/ux455375/bio/evigene/scripts/rnaseq/asmrna_dupfilter2b.pl -CDSALIGN -debug -aasize evg2anofunz4g.cds.qual   -eqgene evgm2cchranofun-evg2anofunz4g.eqgene -blastab evgm2cself-evg2anofunz4g.blastn -outeqtab evgm2cchranofun-evg2anofunz4g.tgalntab -outclass evgm2cchranofun-evg2anofunz4g.tgclass
#t2g: tidyup tmpfiles
#t2g: gzip tmpfiles
#----------------------------
