#! /bin/bash
### env datad=path/to/data trset=myspecies_all.tr.gz qsub -q normal runtr2cds.sh
#PBS -N tr2cds
#PBS -A ind114
#PBS -l nodes=1:ppn=16,walltime=47:55:00
#PBS -V

## FIXME evg4corn : self.blastn fails after retry (2/21 parts, 1st fail 7/21), ncpu too high it appears
## retry ncpu=16 ?
#t2ac: getAaQual: naa=10013588 in evg4corn.aa.qual, val1 cornhi35p4msoapk25loc243751t1= 106,42,-1,complete-utrbad
#t2ac: nofragments_cds== evg4cornnrcd1.cds nrec= 2969415
#t2ac: ERR=blastn_cds_ncpu, nfail=7, nok=14, rerun cmd= /home/ux455375/bio/ncbi2230/bin/blastn -task megablast  -ungapped -xdrop_ungap 4 -dust no -perc_identity 98 -evalue 1e-19 -outfmt 7 -db evg4cornnrcd1_db 
#t2ac: ERR=blastn_cds_ncpu rerun nfail=2, nok=19, fatal cmd= /home/ux455375/bio/ncbi2230/bin/blastn -task megablast  -ungapped -xdrop_ungap 4 -dust no -perc_identity 98 -evalue 1e-19 -outfmt 7 -db evg4cornnrcd1_db 
#t2ac: blastn_cds= evg4cornnrcd1-self98.blastn
#...................................... 

## FIXMEd: upd to tr2aacds2.pl **
# comet.sdsc: caller srun_prog.sh sets ncpu, 24 cores/node, 128G mem
if [ "X" = "X$ncpu" ]; then ncpu=15; fi
if [ "X" = "X$maxmem" ]; then maxmem=50000; fi
evigene=$HOME/bio/evigene/scripts
evapp=$evigene/prot/tr2aacds2.pl

## 2016 new opts for asmrna_dupfilter: aminbad/poo should be same as aamin, as tiny aa usually have "normal" size utr
## new opt only in asmrna_dupfilter3.pl update: aagapmax = BAD_GAPS = 5% better try, was 25% too high
## this aaminbad change adds lots of new, excess fragment loci, need to be filtered out by other means after homol test

export asmrna_dupfilter2=$evigene/rnaseq/asmrna_dupfilter3.pl

export aagapmax=10
export aamin=30
export aapart=120
export aaminbad=30
export aaminpoo=30

traopts="-tidy -log -debug"
if [ "X" != "X$aablast" ]; then traopts="$traopts -ablastab $aablast"; fi

# lastz, etc in bio/bin
export PATH=$HOME/bio/bin:$PATH
# cd-hit
export PATH=$HOME/bio/cdhit/bin:$PATH
# fastanrdb
export fastanrdb=$HOME/bio/exonerate/bin/fastanrdb
# blastn:
## upd: ncbi2230, was ncbi2227 ?
export PATH=$HOME/bio/ncbi2230/bin:$PATH

if [ "X" = "X$trset" ]; then echo "missing env trset=xxxx.tr"; exit -1; fi
if [ "X" = "X$datad" ]; then echo "missing env datad=/path/to/data "; exit -1; fi
cd $datad/

echo "#START $trset : `date`"
echo $evapp -NCPU $ncpu -MAXMEM $maxmem $traopts -cdna $trset
$evapp -NCPU $ncpu -MAXMEM $maxmem $traopts -cdna $trset
echo "#DONE $trset : `date`"

