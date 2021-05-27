#! /bin/bash
### env genome=xxx mrna=mrna.fa datad=`pwd` qsub -q normal genosplign.sh
#PBS -N genosplign 
#PBS -A ind114
#PBS -l nodes=1:ppn=32,walltime=18:55:00
#PBS -o genosplign.$$.out
#PBS -e genosplign.$$.err
#PBS -V

# NCBI splign align mRNA to genome
# MUST assume genome,mrna files exist in datad/ and DONT have full path
# version, new options for compart,splign
dv=o2

## Tests show need weak idty opts to get valid split-gene parts in both compart,splign
## compart opts, need weaker ident to get true mrna aligns to draft/holey genome asm
cpopt="-min_idty 0.25"
splopt="-type mrna -min_compartment_idty 0.25"
## opt for weaker aligns:  -min_compartment_idty below 0.7 default
#v0# splopt="-type mrna -min_compartment_idty 0.5"

ncpu=32
nbin=$HOME/bio/ncbix/bin
evigene=$HOME/bio/evigene
## sdsc: need dang local/lib/libpcre replacement for libpcre.so* copy
export LD_LIBRARY_PATH=$HOME/bio/lib/lib

if [ "X" = "X$datad" ]; then echo "ERR: param datad=What?"; exit -1; fi
if [ "X" = "X$genome" ]; then echo "ERR: param genome=What?"; exit -1; fi
if [ "X" = "X$mrna" ]; then echo "ERR: param mrna=What?"; exit -1; fi

genomez=$genome; genome=`echo $genome | sed 's/.gz//;'`
mrnaz=$mrna; mrna=`echo $mrna | sed 's/.gz//;'`
mname=`basename $mrna .tr | sed 's/\.fasta//; s/\.fa//; s/\.cdna//;'`
gname=`basename $genome .fa | sed 's/\.fasta//; s/\.fa//;'`
## oname="$mname-$gname";

cd $datad/

if [ ! -f $genome.nsq ]; then 
   if [ $genomez = "$genome.gz" ]; then
    gunzip -c $genomez | $nbin/makeblastdb -parse_seqids -dbtype nucl -title $genome -out $genome
   else
    $nbin/makeblastdb -parse_seqids -dbtype nucl -in $genome; 
   fi
fi

## data-parallelize by splitting mrna into ncpu parts
if [ ! -f $mrna.split.1.fa ]; then
 pindir=`dirname $mrna`
 splitsize=`grep -v '^>' $mrna | wc -c | sed 's/ .*//' `
 splitbp=$(( $splitsize / $ncpu ))
 $evigene/scripts/splitMfasta.pl --outputpath=$pindir --minsize=$splitbp $mrna
fi

i=0;
qset=`/bin/ls $mrna.split.*.fa`
for qfile in $qset
{
  qfile1=`basename $qfile`
  qnam=`basename $qfile .fa`
  onam=$dv-$gname-$qnam
  # need subdir per part, use for splign-lds also
  # recall qfile has ".split.$i" sufx  NOTE: split.i IS NOT same as this i in for loop
  # change to this subdir naming?  kfish2p67vs.mrna.split.12.fa => spl12
  # ispldir=`echo $qnam | sed "s/$mrna.split./spl/;"`
  ispldir="spl$i";
  mkdir $ispldir

  ## what path on mrna,genome ?? ln needs to know *
  mv $qfile $ispldir/$qfile1
  
  ## need all of genome.blastdb in ispldir as symlinks
  #bad.wildcard.notthere# ln -s  ../$genome* $ispldir/

  cd $ispldir ;  
  ## worry:  $nbin/splign -mklds is picky: needs exist dir and above 2 files but nothing else??
  ## extra files bad for splign -mklds 
  ln -s ../$genome .
  
  # this is fork set:
  # later add evigene/rnaseq/splign2gff.pl after splign, as final output.
  echo "# splign $qfile1 x $genome in $ispldir";
  ( $nbin/splign -mklds ./; \
  	ln -s ../$genome.* ./ ; \
    $nbin/makeblastdb -parse_seqids -dbtype nucl -in $qfile1; \  
    $nbin/compart $cpopt -qdb $qfile1 -sdb $genome  > $onam.cpart; \  
    $nbin/splign $splopt -comps $onam.cpart -blastdb $genome -ldsdir ./ -log $onam.splog > $onam.splign; ) &

  cd ../ ; # or cd $datad/
  i=$(( $i + 1 ))
}

wait

## tar it till get final output figured  
opack="$mname-$gname"
tar -cf $opack-splign-$dv.tar spl*/$dv-$gname*
## package *.cpart  *.splog also ?? yes for now, want all to test
# cat spl*/*.splign > $opack.splign
# gzip --fast $opack.splign


#...... DANG DAMMIT ........
# make[4]: *** No rule to make target `/home/ux455375/bio/ncbi_cxx--12_0_0/
#   GCC412-Debug64/status/SQLITE3.enabled', needed by `requirements'.  Stop.
# make[4]: Leaving directory `/home/ux455375/bio/ncbi_cxx--12_0_0/GCC412-Debug64/build/app/splign'
# NOTE:  skipping project "splign" due to unmet requirements
## ** Systems trestles, mason have /usr/lib/libsqlite3 and /usr/include/sqlite3.h BUT config test fails
## version problem: sqlite3_pcache_methods not in  /usr/include/sqlite3.h 2009
## #define SQLITE_VERSION         "3.3.6" for 2009; 3.7 is curver
# http://www.sqlite.org/   .. make/install this HOME/bio/lib; cfg ncbix --with-sqlite
# http://www.sqlite.org/2013/sqlite-autoconf-3071700.tar.gz
## ncbix: ./configure --prefix=$HOME/bio/ncbix  --with-sqlite3=$HOME/bio/lib

#........ test run ........
# trestles.sdsc: /home/ux455375/scratchn/chrs/kfish/genes/gsplignf
#   kfish2asm.gdna -> ../../genome/killifish20130322asm.fa
#   kfish2p67vs.mrna -> ../trseq2/kfish2p67vs.mrna_pub.fa
#
#  env genome=kfish2asm.gdna mrna=kfish2p67vs.mrna datad=`pwd` qsub -q normal ../genosplign.sh
# 
# [ux455375@trestles-login2 gsplignf]$ ls
# kfish2asm.gdna	    kfish2asm.gdna.nsd	spl0   spl12  spl16  spl2   spl23  spl27  spl30  spl6
# kfish2asm.gdna.nhr  kfish2asm.gdna.nsi	spl1   spl13  spl17  spl20  spl24  spl28  spl31  spl7
# kfish2asm.gdna.nin  kfish2asm.gdna.nsq	spl10  spl14  spl18  spl21  spl25  spl29  spl4	 spl8
# kfish2asm.gdna.nog  kfish2p67vs.mrna	spl11  spl15  spl19  spl22  spl26  spl3   spl5	 spl9
#
# [ux455375@trestles-login2 gsplignf]$ ls spl10
# kfish2asm.gdna				       kfish2asm.gdna.nsd		kfish2p67vs.mrna.split.1.fa.nog
# kfish2asm.gdna.750219878.idc		       kfish2asm.gdna.nsi		kfish2p67vs.mrna.split.1.fa.nsd
# kfish2asm.gdna-kfish2p67vs.mrna.split.1.cpart  kfish2asm.gdna.nsq		kfish2p67vs.mrna.split.1.fa.nsi
# kfish2asm.gdna.nhr			       kfish2p67vs.mrna.split.1.fa	kfish2p67vs.mrna.split.1.fa.nsq
# kfish2asm.gdna.nin			       kfish2p67vs.mrna.split.1.fa.nhr
# kfish2asm.gdna.nog			       kfish2p67vs.mrna.split.1.fa.nin
# [ux455375@trestles-login2 gsplignf]$ ls -l spl10
#       17 Jul 19 13:41 kfish2asm.gdna -> ../kfish2asm.gdna
#   122160 Jul 19 13:42 kfish2asm.gdna.750219878.idc
#        0 Jul 19 13:41 kfish2asm.gdna-kfish2p67vs.mrna.split.1.cpart
#       21 Jul 19 13:41 kfish2asm.gdna.nhr -> ../kfish2asm.gdna.nhr
#       21 Jul 19 13:41 kfish2asm.gdna.nin -> ../kfish2asm.gdna.nin
#       21 Jul 19 13:41 kfish2asm.gdna.nog -> ../kfish2asm.gdna.nog
#       21 Jul 19 13:41 kfish2asm.gdna.nsd -> ../kfish2asm.gdna.nsd
#       21 Jul 19 13:41 kfish2asm.gdna.nsi -> ../kfish2asm.gdna.nsi
#       21 Jul 19 13:41 kfish2asm.gdna.nsq -> ../kfish2asm.gdna.nsq
#  7839876 Jul 19 13:41 kfish2p67vs.mrna.split.1.fa
#  1887312 Jul 19 13:41 kfish2p67vs.mrna.split.1.fa.nhr
#   127616 Jul 19 13:41 kfish2p67vs.mrna.split.1.fa.nin
#    42540 Jul 19 13:41 kfish2p67vs.mrna.split.1.fa.nog
#   530384 Jul 19 13:41 kfish2p67vs.mrna.split.1.fa.nsd
#    11019 Jul 19 13:41 kfish2p67vs.mrna.split.1.fa.nsi
#  1590403 Jul 19 13:41 kfish2p67vs.mrna.split.1.fa.nsq
#
# **** DID this fail? $nbin/splign -mklds ./;  : no _SplignLDS2_ subdir ...  
# *** yes, fail no /usr/local/lib/libpcre  << this is only /usr/local from ncbix configure.
#........ dammit local/lib lost ........
# /home/ux455375/bio/ncbix/bin/splign: error while loading shared libraries: 
#  libpcre.so.1: cannot open shared object file: No such file or directory
## ** this no good for cluster nodes:
##    export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
##  movit to $HOME/bio/lib/
#......
#  got lots of this also whine about NNN in trseq
# Warning: (1431.1) FASTA-Reader: First data line in seq is about 96% ambiguous nucleotides (shouldn't be over 40%) at line 4905
#...
#    --with-64  << need this for ncbix config / compile?? or not: checking build system type... x86_64-unknown-linux-gnu
# build is GCC412-Debug64 << only 64k
#.........  
# -mklds problems; cant do w/ extra files in subdir
#     $nbin/splign -mklds ./
# Error: (1441.9) Unrecognized file format: /oasis/projects/nsf/ind114/ux455375/chrs/kfish/genes/gsplignf/spl2/./kfish2asm.gdna-kfish2p67vs.mrna.split.12.cpart
# Error: (1441.9) Unrecognized file format: /oasis/projects/nsf/ind114/ux455375/chrs/kfish/genes/gsplignf/spl2/./kfish2asm.gdna-kfish2p67vs.mrna.split.12.splign
# Warning: (1441.6) Unrecognized top level object in /oasis/projects/nsf/ind114/ux455375/chrs/kfish/genes/gsplignf/spl2/./kfish2asm.gdna.nhr

