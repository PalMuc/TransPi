#! /bin/bash
### env tag=t2 hits=xxx.mblastn genome=xxx mrna=mrna.fa datad=`pwd` qsub -q normal genosplign15.sh
#PBS -N genosplign 
#PBS -A ind114
#PBS -l nodes=1:ppn=16,walltime=18:55:00
#PBS -V

# rev 2015 for input of mrna_genome.mblast hit table, dont need other parts of splign
## fixme: need subset tag: t1, t2, t3, from hits name or param?
# need new method to split mblast table to ncpu parts
# NCBI splign align mRNA to genome
# MUST assume genome,mrna files exist in datad/ and DONT have full path

ncpu=15
#tgrid# 
nbin=$HOME/bio/ncbix/bin
#dgg# nbin=/bio/bio-grid/mb/ncbix/bin
evigene=$HOME/bio/evigene

## need also libpcre not in std set: BUT this aint on cluster nodes !
# export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$HOME/bio/lib/lib 

## >> opt for weaker aligns:  -min_compartment_idty below 0.7 default
# OLDsplopt="-type mrna -min_compartment_idty 0.5"
splopt="-type mrna -direction sense -min_compartment_idty 0.75 -compartment_penalty 0.20 -min_exon_idty 0.75 -min_singleton_idty_bps 150"

if [ "X" = "X$datad" ]; then echo "ERR: param datad=What?"; exit -1; fi
if [ "X" = "X$genome" ]; then echo "ERR: param genome=What?"; exit -1; fi
if [ "X" = "X$mrna" ]; then echo "ERR: param mrna=What?"; exit -1; fi
if [ "X" = "X$hits" ]; then echo "ERR: param hits=What?"; exit -1; fi

## dont allow gz here..
# genomez=$genome; genome=`echo $genome | sed 's/.gz//;'`
# mrnaz=$mrna; mrna=`echo $mrna | sed 's/.gz//;'`
mname=`basename $mrna .mrna | sed ' s/\.cdna//; s/\.tr//; s/\.fasta//; s/\.fa//;'`
gname=`basename $genome .fa | sed 's/\.fasta//; s/\.fa//;'`
## oname="$mname-$gname";

cd $datad/

## data-parallelize by splitting hits into ncpu parts; no good here need split by Funhe gene ids
#nogood# nhit=`wc -l $hits`; nsplit=$(( $nhit / $ncpu )) ; split -a 2 -l $nsplit $hits $hits.split.

##splnam=`basename $hits .mblastn`.split
splnam=`basename $hits`.split

## add tag to ispldir ..
spltag=`basename $hits .mblastn`
if [ "X" = "X$tag" ]; then tag=$spltag; fi

# fixme: input hits sort by -k1,1; output sort by -k2,2 -k1,1 for splign 
if [ ! -f $splnam.01 ]; then
  sort -k1,1 -k2,2 $hits | env ncpu=$ncpu snam=$splnam  perl -ne \
'BEGIN{ $ncpu=$ENV{ncpu}; $sn=$ENV{snam}; $nid=$lid=0;
for($i=0; $i<$ncpu; $i++) { $sni=$sn.sprintf ".%02d",(1+$i); open(my $fh,">",$sni); $fh[$i]=$fh; } }
($id)=split; if($id ne $lid) { unless($fh=$fhd{$id}) { 
$sid= $nid % $ncpu; $fhd{$id}=$fh= $fh[$sid]; $nid++; } } 
print $fh $_; $lid=$id;'
  
fi

echo "START `date` "
i=0;
qset=`/bin/ls $splnam.*`
## $i and split.i are diff now due to ls sort split.11 > split.2
for qfile in $qset
{
  j=$(( $i + 1 ))
  ispldir="spld$tag$j";
  qfile1=`basename $qfile`
  qnam=`basename $qfile`
  onam=$qnam

  ##sort fix# mv $qfile $ispldir/$qfile1
  ##buggers, need to run $nbin/splign -mklds before copy in qfile1
  mkdir $ispldir
  cd $ispldir ; 
  ln -s ../$genome .
  ln -s ../$mrna .
  $nbin/splign -mklds ./
  sort -k2,2 -k1,1 ../$qfile > $qfile1
  
  echo "# splign -hits $qfile1 in $ispldir";
  echo $nbin/splign $splopt -hits $qfile1 -log $onam.splog -ldsdir ./ TO $onam.splign > $onam.cmdline
  $nbin/splign $splopt -hits $qfile1 -log $onam.splog -ldsdir ./ > $onam.splign  &

  cd ../ ; 
  i=$(( $i + 1 ))
  if [ $i -ge $ncpu ]; then wait; i=0; fi
}

wait
echo "DONE `date` "

## tar it till get final output figured out? 
## need splog along w/ splign, separate files
# opack="$gname-$mname"
#? tar --exclude=$genome --exclude=$mrna --exclude=_SplignLDS2_ -cf $opack.tar spl*/
#no# cat spl*/*.splign > $opack.splign; gzip --fast $opack.splign

exit;
#................
## OLD
# if [ ! -f $genome.nsq ]; then 
#    if [ $genomez = "$genome.gz" ]; then
#     gunzip -c $genomez | $nbin/makeblastdb -parse_seqids -dbtype nucl -title $genome -out $genome
#    else
#     $nbin/makeblastdb -parse_seqids -dbtype nucl -in $genome; 
#    fi
# fi
## old
# if [ ! -f $mrna.split.1.fa ]; then
#  pindir=`dirname $mrna`
#  splitsize=`grep -v '^>' $mrna | wc -c | sed 's/ .*//' `
#  splitbp=$(( $splitsize / $ncpu ))
#  $evigene/scripts/splitMfasta.pl --outputpath=$pindir --minsize=$splitbp $mrna
# fi
#
# qset=`/bin/ls $mrna.split.*.fa`
#  # OLD this is fork set:
#   echo "# splign $qfile1 x $genome in $ispldir";
#   ( $nbin/splign -mklds ./; \
#     ln -s ../$genome.* ./ ; \
#     $nbin/makeblastdb -parse_seqids -dbtype nucl -in $qfile1; \  
#     $nbin/compart  -qdb $qfile1 -sdb $genome  > $onam.cpart; \  
#     $nbin/splign $splopt -comps $onam.cpart -blastdb $genome -ldsdir ./ -log $onam.splog > $onam.splign; ) &
