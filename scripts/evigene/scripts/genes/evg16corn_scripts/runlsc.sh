#! /bin/bash
### env inlong=longreads.fa inpe=shortreads.fa datad=`pwd` qsub -q normal runlongsc.sh
#PBS -N longsc
#PBS -A ind114
#PBS -l nodes=1:ppn=16,walltime=19:55:00
#PBS -V

module add samtools bowtie2

lbin=$HOME/bio/longsc/bin
evigene=$HOME/bio/evigene
bopts=""

if [ "X" = "X$outdir" ]; then outdir="lsco"; echo "WARN: env outdir=$outdir ... default"; fi
if [ "X" = "X$ncpu" ]; then ncpu=8; fi
if [ "X" != "X$lscopts" ]; then bopts=$lscopts; fi

if [ "X" = "X$datad" ]; then echo "err missing datad=path/to/data"; exit -1; fi
if [ "X" = "X$inlong" ]; then echo "err missing inlong=longreads.fa"; exit -1; fi
if [ "X" = "X$inshort" ]; then echo "err missing inshort=reads_1.fa"; exit -1; fi
# if [ "X" = "X$name" ]; then echo "err missing name=velsoapgroup"; exit -1; fi
# for large sort, /tmp not good
tmpdir=tmp$outdir

## caller sets OMP ncpu ?
## export OMP_NUM_THREADS=$ncpu

cd $datad/

# inrt=`echo $inpe | sed 's/_1\./_2./;'`
# mkdir $outdir
mkdir $tmpdir
bopts="$bopts --tempdir $tmpdir"

echo "START: `date`"
echo "$lbin/runLSC.py $bopts -o $outdir --long_reads $inlong --short_reads $inshort --threads $ncpu"
$lbin/runLSC.py $bopts -o $outdir --long_reads $inlong --short_reads $inshort --threads $ncpu

echo "DONE : `date`"

#................
#  parser.add_argument('--long_reads',help="FASTAFILE Long reads to correct. Required in mode 0 or 1.")
#  parser.add_argument('--short_reads',nargs='*',help="FASTA/FASTQ FILE Short reads used to correct the long reads. 
#  #fa: parser.add_argument('--short_read_file_type',default='fa',choices=['fa','fq','cps'],help="Short read file type")
#  parser.add_argument('--threads',type=int,default=0,help="Number of threads (Default = cpu_count)")
#  parser.add_argument('-o','--output',help="FOLDERNAME where output is to be written. Required in mode 0 or 3.")
#  parser.add_argument('--tempdir',default='/tmp',help="FOLDERNAME where temporary files can be placed")
#  parser.add_argument('--sort_mem_max',type=int,help="-S option for memory in unix sort")
#  parser.add_argument('--aligner',default='bowtie2',choices=['hisat','bowtie2'],help="Aligner choice. hisat parameters have not been optimized, so we recommend bowtie2.")
#  parser.add_argument('--sort_mem_max',type=int,help="-S option for memory in unix sort")
#  parser.add_argument('--minNumberofNonN',type=int,default=40,help="Minimum number of non-N characters in the compressed read")
#  parser.add_argument('--maxN',type=int,help="Maximum number of Ns in the compressed read")
#  parser.add_argument('--error_rate_threshold',type=int,default=12,help="Maximum percent of errors in a read to use the alignment")
#   parser.add_argument('--short_read_coverage_threshold',type=int,default=20,help="Minimum short read coverage to do correction")
#   parser.add_argument('--long_read_batch_size',default=500,type=int,help="INT number of long reads to work on at a time.  This is a key parameter to adjusting performance.  A smaller batch size keeps the sizes and runtimes of intermediate steps tractable on large datasets, but can slow down execution on small datasets.  The default value should be suitable for large datasets.")
#   parser.add_argument('--samtools_path',default='samtools',help="Path to samtools by default assumes its installed.  
   
 
if [ "X" = "X$datad" ]; then echo "err missing datad=path/to/data"; exit -1; fi
if [ "X" = "X$inlong" ]; then echo "err missing inlong=longreads.fa"; exit -1; fi
if [ "X" = "X$inshort" ]; then echo "err missing inshort=reads_1.fa"; exit -1; fi
# if [ "X" = "X$name" ]; then echo "err missing name=velsoapgroup"; exit -1; fi

cd $datad/
# inrt=`echo $inpe | sed 's/_1\./_2./;'`
# mkdir $outdir

echo "START: `date`"
echo "$lbin/runLSC.py $bopts --long_reads $inlong --short_reads $inshort --threads $ncpu"

$lbin/runLSC.py $bopts --long_reads $inlong --short_reads $inshort --threads $ncpu

echo "DONE : `date`"

