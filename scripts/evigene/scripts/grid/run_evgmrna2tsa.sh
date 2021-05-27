#! /bin/bash
### env idprefix=MysppEGm trclass=myspp.trclass  [vectrim=1 tbl2asn=1 species="Xxx yyy" ]  
###   datad=`pwd`  qsub -q shared run_evgmrna2tsa.sh
#PBS -N evgmrna2tsa 
#PBS -l nodes=1:ppn=8,walltime=5:55:00
#PBS -o egmrna2tsa.$$.out
#PBS -e egmrna2tsa.$$.err
#PBS -V

ncpu=8; # most 100K trsets run on 8cpu in 10m-20min; 300K set can take 1hr.
## common option: skip vecscreen and tbl2asn == -novectrim (not def)  -noruntbl2asn (default)
evigene=$HOME/bio/evigene/scripts

#old#export vecscreen=$HOME/bio/ncbic11/bin/vecscreen
## need to test ncbic++ version of vecscreen; uses makeblastdb ..
export vecscreen=$HOME/bio/ncbix/bin/vecscreen
export tbl2asn=$HOME/bio/ncbix/bin/tbl2asn

## FIXME: -TSADESC=tbl2asn flags for path/to/tsa.cmt files
## default $TSADESC="-w evgr_tsamethods.cmt -Y evgr_tsadesc.cmt -t evgr_tsasubmit.sbt";

## opts="-debug -runtbl2asn -NCPU $ncpu"
opts="-debug -NCPU $ncpu"
if [ "X" = "X$datad" ]; then echo "missing env datad=path/to/data"; exit -1; fi
if [ "X" = "X$trclass" ]; then "echo env trclass=path/to/name.trclass"; exit -1; fi
## .. these are now read via sra_result.csv, species => idprefix
## FIXME: \' bad, passed to perl opts.. remove ' ' for _
if [ "X" != "X$idprefix" ]; then opts="$opts -idprefix $idprefix"; fi
if [ "X" != "X$species" ]; then spp=`echo $species | sed 's/ /_/g;'`; opts="$opts -species=$spp"; fi
if [ "X" != "X$sra" ]; then sr=`echo $sra | sed 's/ /,/g;'`;  opts="$opts -sraids=$sr"; fi
#bad?# if [ "X" != "X$sra" ]; then opts="$opts -sraids=\'$sra\'"; fi

if [ "X" = "X$vectrim" ]; then opts="$opts -novectrim"; fi
if [ "X" != "X$tbl2asn" ]; then opts="$opts -runtbl2asn"; fi

cd $datad/
## add -outdir opt or : mkdir tsasubmit; cd tsasubmit
##FIXME: template files ; need these  or generate defaults?
if [ ! -f evgr_tsasubmit.sbt ]; then
	if [ -d ../tsasubmit ]; then cp ../tsasubmit/*.{cmt,sbt} ./; fi
fi
	
echo $evigene/evgmrna2tsa.pl  $opts -log -class $trclass
$evigene/evgmrna2tsa.pl  $opts -log -class $trclass
