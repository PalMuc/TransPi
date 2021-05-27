#!/bin/bash
# fq2pairsr.sh
### env datad=`pwd` fastq=pairs_1.fq.gz qsub -q shared fq2pairsr.sh
#PBS -N fq2pairsr
#PBS -l vmem=16gb,nodes=1:ppn=2,walltime=28:55:00
#PBS -o fq2pairsr.$$.out
#PBS -e fq2pairsr.$$.err
#PBS -V

ncpu=2

if [ "X" = "X$fastq" ]; then echo "ERR: fastq= what?"; exit -1; fi
if [ "X" = "X$datad" ]; then echo "ERR: datad= what?"; exit -1; fi
cd $datad/
export TMPDIR=$datad/

for infq1 in $fastq; do {

nam=`echo $infq1 | sed 's/_[12]\..*//;'`
infq2=`echo $infq1 | sed 's/_1/_2/;'`
if [ ! -f $infq2 ]; then echo "ERR: missing fq2=$infq2"; continue; fi

CAT=cat;
innoz=`echo $infq1 | sed 's/.gz//;'`
if [ "$innoz.gz" = $infq1 ]; then CAT="gunzip -c"; fi

p=1; 
( $CAT $infq1 | lr=$p perl -e \
'BEGIN{$i=$ENV{lr};} while(<>){ ($fh)=split; $fs=<>; $qh=<>; $qs=<>; $fh=~s/^\@/>/; $fh.="/$i" unless($fh=~m,/[12]$,); print "$fh\t$fs";} ' \
 > $nam.fap$p ) &
# dang blank lines printed: chomp($fs) or no \n

p=2;
( $CAT $infq2 | lr=$p perl -e \
'BEGIN{$i=$ENV{lr};} while(<>){ ($fh)=split; $fs=<>; $qh=<>; $qs=<>; $fh=~s/^\@/>/; $fh.="/$i" unless($fh=~m,/[12]$,); print "$fh\t$fs";} ' \
 > $nam.fap$p ) &

wait

## this sort will take a while.. dang run1 had lots o blank lines
sort -T ./  $nam.fap[12] | env nam=$nam perl -ne\
'BEGIN{ $nam=$ENV{nam}; open(PE,">$nam.pe.fa2"); open(SR,">$nam.sr.fa"); }
next unless(/^>/); $nin++; ($pd,$sq)=split; ($id,$pn)=split"/",$pd;
if($pn == 2) { 
if($id eq $lid) { print PE "$lpd\n$lsq\n$pd\n$sq\n"; $npe++; $id=$pn=0;}
else { if($lpd) { print SR "$lpd\n$lsq\n"; $nsr++; }      
print SR "$pd\n$sq\n";  $nsr++; $id=$pn=0; } } 
elsif($pn == 1) { if($lpd) { print SR "$lpd\n$lsq\n"; $nsr++; } } 
else { $nerr++; } $lpd=$pd; $lsq=$sq; $lpn=$pn; $lid=$id; 
END{close(PE); close(SR); $no= $npe*2 + $nsr;
warn "#wrote nin=$nin, nout=$no, npe=$npe x 2, nsr=$nsr, nerr=$nerr to $nam.pe,sr.fa2\n";} '

} done

#wrote npe=211569956 x 2, nsr=1009184 to kf_whoi_rnaseq_reads.pe.fa2,sr.fa
#wrote npe=211080234 x 2, nsr=1926923 to kf_mdibl_rnaseq_reads.pe,sr.fa2

## bug: sr.fa has pairs:
# >FCD2AYNACXX:8:1201:11732:92871#ATTCCTTT/1
# GCCAGCTGCTTCTCGAACCTGGGGGCAACAAAAGACCAGGGACCCATGTTTTGAGGTTCCTCCTGACTCCAAACAAAGTCTTTGGCATTGGGGTACTTCC
# >FCD2AYNACXX:8:1201:11732:92871#ATTCCTTT/2
# GTGGTGCTGTGCTCCGGGAAGCATTATTACGCTCTGCTGAAACAGAGGGAGACATCAGCAGCCAACCAGAACACAGCGCTCATCCGTGTGGAGGAGCTGT
