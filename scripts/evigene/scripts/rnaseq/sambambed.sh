#! /bin/bash
### env rnain=xxx/xxx.sam qsub -q normal sambambed.sh
#PBS -N sambambed
#PBS -l nodes=1:ppn=24,walltime=9:55:00
#PBS -o sambambed.$$.out
#PBS -e sambambed.$$.err
#PBS -V

# for gsnap partitioned outputs.sam format
# peset="concordant_mult|concordant_uniq|paired_mult|paired_uniq_long|paired_uniq_inv|paired_uniq_scr"
# srset="halfmapping_mult|halfmapping_uniq|unpaired_mult|unpaired_uniq"
peuniq="concordant_uniq"
pemult="concordant_mult|paired_mult"
halfs="halfmapping_mult|halfmapping_uniq"
noset="nomapping"

bedopts="debug=1 log=1 ave=1 opt=-A"

module add samtools
evigene=$HOME/bio/evigene
# bam2bw=$evigene/scripts/rnaseq/bam2bw.pl

if [ "X" = "X$allsam" ]; then echo "ERR: miss allsam=gsodaplx/xxx*_.{concordant_uniq,concordant_mult,..}"; exit -1; fi
if [ "X" = "X$genoasm" ]; then echo "ERR: miss genoasm=genoasm.fa"; exit -1; fi
if [ "X" = "X$ncpu" ]; then ncpu=16; fi
if [ "X" = "X$datad" ]; then echo "ERR: datad=what?"; exit -1; fi

odir=sambambed$$

echo "start sambambed: `date` "
cd $datad
mkdir $odir

insam=`ls $allsam | egrep "$peuniq"`
insam=`ls $allsam | egrep "$pemult"`
insam=`ls $allsam | egrep "$halfs"`

if [ "X" = "X$samo" ]; then
  samo=($insam)
  samo=$samo.srt
fi

  # fai made from samtools
  # genosize=$genoasm.fai 
  
  cat $insam | samtools view -T $genoasm -S -u - | samtools sort -l 9 -o $samo.bam -
  env $bedopts $evigene/scripts/rnaseq/bam2bw.pl $genoasm.fai  $samo.bam > $samo.bed ;
  gzip --fast  $samo.bed; 


#   j=$(( $j + 1 ))
#   i=$(( $i + 1 ))
#   if [ $i -ge $ncpu ]; then wait; i=0; fi
# 
# } done
# wait

echo "end sambambed: `date` "

#.... test
# genosize=../../genome/gasm16ml/daphplx_gasm16ml.fa.count
# genoasm=../../genome/gasm16ml/daphplx_gasm16ml.fa
# insam=gsodaplx16gdnsurs/*gsnap.halfmapping_mult
# .. samtools sort - > $samo.bam OR samtools sort -l 9 -o $samo.bam -
##... new samtools ...
#  samtools view -T $genosize << changed to geno.fasta
# [ux455375@comet-ln3 rnamap]$ cat $insam | samtools view -T $genosize -S -u - | samtools sort - $samo
# [bam_sort] Use -T PREFIX / -o FILE to specify temporary and final output files
# Usage: samtools sort [options...] [in.bam]
# Options:
#   -l INT     Set compression level, from 0 (uncompressed) to 9 (best)
#   -m INT     Set maximum memory per thread; suffix K/M/G recognized [768M]
#   -n         Sort by read name
#   -o FILE    Write final output to FILE rather than standard output
#   -T PREFIX  Write temporary files to PREFIX.nnnn.bam
#   -@, --threads INT
#              Set number of sorting and compression threads [1]
#       --input-fmt-option OPT[=VAL]
#                Specify a single input file format option in the form
#                of OPTION or OPTION=VALUE
#   -O, --output-fmt FORMAT[,OPT[=VAL]]...
#                Specify output format (SAM, BAM, CRAM)
#       --output-fmt-option OPT[=VAL]
#                Specify a single output file format option in the form
#                of OPTION or OPTION=VALUE
#       --reference FILE
#                Reference sequence FASTA FILE [null]
