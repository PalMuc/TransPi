#! /bin/bash
### qsub -q batch bow2mag4trmap.sh
#PBS -N bowtrmap
#... -A ind114
#....  nodes=1:ppn=8,walltime=21:55:00
#PBS -l nodes=1:ppn=28,walltime=21:55:00
#PBS -o bowtrmap.$$.out
#PBS -e bowtrmap.$$.err
#PBS -V

# rsetn=7 x 4cpu/file = 28 total, add 1 file to get 32cpu?

# new bowtie2 may2012 inst on mason.iu; handles fq.gz; opts: --very-fast --fast --sensitive ..
module add bowtie/2.0.0_b6

#t1. trdb=trsets/goodalt3m8veltribest2
# odir=mag4trmapbt
trdb=trsets/daphmag3mtv3
odir=mag4tr3bt

ncpu=8
datad=$HOME/scratch
workd=$datad/chrs/daphmag/rnas/

# test case
lreads=$workd/fastq4/Dman21_TTAGGC_L006_R1_001.fastq.gz

## groups 2012.05.12 n=45
#   6 L002
#   3 L003
#   3 L004
#   8 L005
#  15 L006
#   4 L007
#   6 L008
ncpu=4
rset1="Dman_03_CGATGT_L005 Dman_13_ACAGTG_L006 Dman_23_TGACCA_L002 Dman_32_ACAGTG_L003 Dman_42_CAGATC_L004 Dman_60_TGACCA_L007 Dman_64_CAGATC_L008 "
rset2="Dman_05_TTAGGC_L005 Dman_15_CGATGT_L006 Dman_27_ACAGTG_L002 Dman_33_GCCAAT_L003 Dman_45_ACTTGA_L004 Dman_61_ACAGTG_L007 Dman_65_ACTTGA_L008 "
rset3="Dman_09_TGACCA_L005 Dman_21_TTAGGC_L006 Dman_29_GCCAAT_L002 Dman_53_ATCACG_L004 Dman_62_GCCAAT_L007 Dman_67_TGACCA_L008 Dman_74_ACTTGA_L003 "
rset4="Dman_31_CAGATC_L002 Dman_38_ATCACG_L005 Dman_46_ACAGTG_L006 Dman_68_ACAGTG_L008 Dman_78_TGACCA_L007 "
rset5="Dman_35_ACTTGA_L002 Dman_39_CGATGT_L005 Dman_47_GCCAAT_L006 Dman_70_GCCAAT_L008 "
rset6="Dman_36_GATCAG_L002 Dman_41_TTAGGC_L005 Dman_48_ATCACG_L006 Dman_73_CAGATC_L008 "
# rset7="Dman_44_TGACCA_L005 Dman_50_CAGATC_L006 "
# rset8="Dman_51_CGATGT_L006 Dman_58_TGACCA_L005 "
# rset9="Dman_52_ACTTGA_L006 "  && more to 15

if [ $iset = 1 ]; then rset=$rset1; fi
if [ $iset = 2 ]; then rset=$rset2; fi
if [ $iset = 3 ]; then rset=$rset3; fi
if [ $iset = 4 ]; then rset=$rset4; fi
if [ $iset = 5 ]; then rset=$rset5; fi
if [ $iset = 6 ]; then rset=$rset6; fi

#bowtie2: .. only sam output, -S outfile.sam?
#?? does --no-unal conflict with -un nomap.fq ? seems so; keep unal in .sam, drop un.fq
# opts="-q --threads $ncpu --very-fast -M 50 --no-head " # bad M opt, took ~14hr for 69M reads Dman21
# o2: 1.5hr for 58M reads Dman_70
opts="-q --mm --threads $ncpu --fast -k 10 --no-head"
# -M or -k 50 for  max aligns = need many for alt-tr
# --mm   use memory-mapped I/O for index; many 'bowtie's can share < use for multiruns?
# --no-mixed = no unpaired of pairs --no-discordant

cd $workd/
mkdir $odir
tnam=`basename $trdb .tr`

# for lreads in fastq4/*_R1_001.fastq.gz; do  
for rdnam in $rset; do {

lreads=$workd/fastq4/${rdnam}_R1_001.fastq.gz
rreads=`echo $lreads | sed 's/_R1/_R2/;'`
# rdnam=`basename $lreads _R1_001.fastq.gz `

echo bowtie2 $opts -x $trdb -1 $lreads -2 $rreads -S $odir/$tnam-$rdnam.sam 2>$odir/$tnam-$rdnam.blog
bowtie2 $opts -x $trdb -1 $lreads -2 $rreads -S $odir/$tnam-$rdnam.sam 2>$odir/$tnam-$rdnam.blog  &

} done

wait

#^^ NOT --mm;
# --mm : This facilitates memory-efficient parallelization of `bowtie` in
# situations where using `-p` is not possible or not preferable.

#... count reads/gene .. add here?

# pt=Dman21_TTAGGC_L006; cut -f3  mag4trmapbt/goodalt3*$pt.sam | perl -ne '$c{$_}++; $n++; \
# END{ print "#total\t$n\n"; foreach $t (sort keys %c) { $c=$c{$t}; chomp($t); print "$t\t$c\n"; } }'\
#  > mag4trmapbt/$pt.trcperl &

#.....................................................................
# # --threads 8 --very-fast -M 50 ; 14 hr
# ==> mag4trmapbt/goodalt3m8veltribest2-Dman21_TTAGGC_L006.blog <==
# 68881858 reads; of these:
#   68881858 (100.00%) were paired; of these:
#     27762051 (40.30%) aligned concordantly 0 times
#     22568258 (32.76%) aligned concordantly exactly 1 time
#     18551549 (26.93%) aligned concordantly >1 times
#     ----
#     27762051 pairs aligned concordantly 0 times; of these:
#       1341653 (4.83%) aligned discordantly 1 time
#     ----
#     26420398 pairs aligned 0 times concordantly or discordantly; of these:
#       52840796 mates make up the pairs; of these:
#         37094066 (70.20%) aligned 0 times
#         7142839 (13.52%) aligned exactly 1 time
#         8603891 (16.28%) aligned >1 times
# 73.07% overall alignment rate
# 
# 
# # --threads 8 --fast -k 10 ; 1.5 hr
# ==> mag4trmapbt/goodalt3m8veltribest2-Dman_70_GCCAAT_L008.blog <==
# 57601957 reads; of these:
#   57601957 (100.00%) were paired; of these:
#     20862572 (36.22%) aligned concordantly 0 times
#     20371323 (35.37%) aligned concordantly exactly 1 time
#     16368062 (28.42%) aligned concordantly >1 times
#     ----
#     20862572 pairs aligned concordantly 0 times; of these:
#       1112213 (5.33%) aligned discordantly 1 time
#     ----
#     19750359 pairs aligned 0 times concordantly or discordantly; of these:
#       39500718 mates make up the pairs; of these:
#         26000636 (65.82%) aligned 0 times
#         5660574 (14.33%) aligned exactly 1 time
#         7839508 (19.85%) aligned >1 times
# 77.43% overall alignment rate
#


# -----
# Bowtie 2 version 2.0.0-beta6 by Ben Langmead (blangmea@jhsph.edu)
# Usage: 
#   bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r>} [-S <sam>]
# 
#   <bt2-idx>  Index filename prefix (minus trailing .X.bt2).
#              NOTE: Bowtie 1 and Bowtie 2 indexes are not compatible.
#   <m1>       Files with #1 mates, paired with files in <m2>.
#              Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
#   <m2>       Files with #2 mates, paired with files in <m1>.
#              Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
#   <r>        Files with unpaired reads.
#              Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
#   <sam>      File for SAM output (default: stdout)
# 
#   <m1>, <m2>, <r> can be comma-separated lists (no whitespace) and can be
#   specified many times.  E.g. '-U file1.fq,file2.fq -U file3.fq'.
# 
# Options (defaults in parentheses):
# 
#  Input:
#   -q                 query input files are FASTQ .fq/.fastq (default)
#   --qseq             query input files are in Illumina's qseq format
#   -f                 query input files are (multi-)FASTA .fa/.mfa
#   -r                 query input files are raw one-sequence-per-line
#   -c                 <m1>, <m2>, <r> are sequences themselves, not files
#   -s/--skip <int>    skip the first <int> reads/pairs in the input (none)
#   -u/--upto <int>    stop after first <int> reads/pairs (no limit)
#   -5/--trim5 <int>   trim <int> bases from 5'/left end of reads (0)
#   -3/--trim3 <int>   trim <int> bases from 3'/right end of reads (0)
#   --phred33          qualities are Phred+33 (default)
#   --phred64          qualities are Phred+64
#   --int-quals        qualities encoded as space-delimited integers
# 
#  Presets:                 Same as:
#   For --end-to-end:
#    --very-fast            -D 5 -R 1 -N 0 -L 22 -i S,0,2.50
#    --fast                 -D 10 -R 2 -N 0 -L 22 -i S,0,2.50
#    --sensitive            -D 15 -R 2 -N 0 -L 22 -i S,1,1.25 (default)
#    --very-sensitive       -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
# 
#   For --local:
#    --very-fast-local      -D 5 -R 1 -N 0 -L 25 -i S,1,2.00
#    --fast-local           -D 10 -R 2 -N 0 -L 22 -i S,1,1.75
#    --sensitive-local      -D 15 -R 2 -N 0 -L 20 -i S,1,0.75 (default)
#    --very-sensitive-local -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
# 
#  Alignment:
#   -N <int>           max # mismatches in seed alignment; can be 0 or 1 (0)
#   -L <int>           length of seed substrings; must be >3, <32 (22)
#   -i <func>          interval between seed substrings w/r/t read len (S,1,1.25)
#   --n-ceil <func>    func for max # non-A/C/G/Ts permitted in aln (L,0,0.15)
#   --dpad <int>       include <int> extra ref chars on sides of DP table (15)
#   --gbar <int>       disallow gaps within <int> nucs of read extremes (4)
#   --ignore-quals     treat all quality values as 30 on Phred scale (off)
#   --nofw             do not align forward (original) version of read (off)
#   --norc             do not align reverse-complement version of read (off)
# 
#   --end-to-end       entire read must align; no clipping (on)
#    OR
#   --local            local alignment; ends might be soft clipped (off)
# 
#  Scoring:
#   --ma <int>         match bonus (0 for --end-to-end, 2 for --local) 
#   --mp <int>         max penalty for mismatch; lower qual = lower penalty (6)
#   --np <int>         penalty for non-A/C/G/Ts in read/ref (1)
#   --rdg <int>,<int>  read gap open, extend penalties (5,3)
#   --rfg <int>,<int>  reference gap open, extend penalties (5,3)
#   --score-min <func> min acceptable alignment score w/r/t read length
#                      (G,20,8 for local, L,-0.6,-0.6 for end-to-end)
# 
#  Reporting:
#   -M <int>           look for up to <int>+1 alns; report best, with MAPQ (5 for
#                      --end-to-end, 2 for --local)
#    OR
#   -k <int>           report up to <int> alns per read; MAPQ not meaningful (off)
#    OR
#   -a/--all           report all alignments; very slow (off)
# 
#  Effort:
#   -D <int>           give up extending after <int> failed extends in a row (15)
#   -R <int>           for reads w/ repetitive seeds, try <int> sets of seeds (2)
# 
#  Paired-end:
#   -I/--minins <int>  minimum fragment length (0)
#   -X/--maxins <int>  maximum fragment length (500)
#   --fr/--rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (--fr)
#   --no-mixed         suppress unpaired alignments for paired reads
#   --no-discordant    suppress discordant alignments for paired reads
#   --no-dovetail      not concordant when mates extend past each other
#   --no-contain       not concordant when one mate alignment contains other
#   --no-overlap       not concordant when mates overlap at all
# 
#  Output:
#   -t/--time          print wall-clock time taken by search phases
#   --un <path>           write unpaired reads that didn't align to <path>
#   --al <path>           write unpaired reads that aligned at least once to <path>
#   --un-conc <path>      write pairs that didn't align concordantly to <path>
#   --al-conc <path>      write pairs that aligned concordantly at least once to <path>
#   (Note: for --un, --al, --un-conc, or --al-conc, add '-gz' to the option name, e.g.
#   --un-gz <path>, to gzip compress output, or add '-bz2' to bzip2 compress output.)
#   --quiet            print nothing to stderr except serious errors
#   --met-file <path>  send metrics to file at <path> (off)
#   --met-stderr       send metrics to stderr (off)
#   --met <int>        report internal counters & metrics every <int> secs (1)
#   --no-unal          supppress SAM records for unaligned reads
#   --no-head          supppress header lines, i.e. lines starting with @
#   --no-sq            supppress @SQ header lines
#   --rg-id <text>     set read group id, reflected in @RG line and RG:Z: opt field
#   --rg <text>        add <text> ("lab:value") to @RG line of SAM header.
#                      Note: @RG line only printed when --rg-id is set.
# 
#  Performance:
#   -o/--offrate <int> override offrate of index; must be >= index's offrate
#   -p/--threads <int> number of alignment threads to launch (1)
#   --reorder          force SAM output order to match order of input reads
#   --mm               use memory-mapped I/O for index; many 'bowtie's can share
# 
#  Other:
#   --qc-filter        filter out reads that are bad according to QSEQ filter
#   --seed <int>       seed for random number generator (0)
#   --version          print version information and quit
#   -h/--help          print this usage message
# 
# 
# =====
# 
# 
# Usage: bowtie2-build [options]* <reference_in> <bt2_index_base>
#     reference_in            comma-separated list of files with ref sequences
#     bt2_index_base          write .bt2 data to files with this dir/basename
# *** Bowtie 2 indexes work only with v2 (not v1).  Likewise for v1 indexes. ***
# Options:
#     -f                      reference files are Fasta (default)
#     -c                      reference sequences given on cmd line (as <seq_in>)
#     -a/--noauto             disable automatic -p/--bmax/--dcv memory-fitting
#     -p/--packed             use packed strings internally; slower, uses less mem
#     --bmax <int>            max bucket sz for blockwise suffix-array builder
#     --bmaxdivn <int>        max bucket sz as divisor of ref len (default: 4)
#     --dcv <int>             diff-cover period for blockwise (default: 1024)
#     --nodc                  disable diff-cover (algorithm becomes quadratic)
#     -r/--noref              don't build .3/.4.bt2 (packed reference) portion
#     -3/--justref            just build .3/.4.bt2 (packed reference) portion
#     -o/--offrate <int>      SA is sampled every 2^offRate BWT chars (default: 5)
#     -t/--ftabchars <int>    # of chars consumed in initial lookup (default: 10)
#     --seed <int>            seed for random number generator
#     -q/--quiet              verbose output (for debugging)
#     -h/--help               print detailed description of tool and its options
#     --usage                 print this usage message
#     --version               print version information and quit
# 
