#! /bin/bash
### env name=outname inpe=readpairs.fa datad=path/to/data this.sh
#PBS -N tr2cds
#PBS -A ind114
#PBS -l nodes=1:ppn=16,walltime=47:55:00
#PBS -V

## runidbat2l == long_read len=151
## runidbat2r == bin2/ has short read maxlen 158 for longer shorts **
if [ "X" = "X$ncpu" ]; then ncpu=15; fi
if [ "X" = "X$maxmem" ]; then maxmem=50000; fi

# ov="o1"; traopts="--mink 27 --maxk 97 --step 10 --max_isoforms 21";
# ov="r2"; traopts="--mink 31 --maxk 91 --step 10 ";
ov="r3"; traopts="--mink 27 --maxk 117 --step 10 ";

evigene=$HOME/bio/evigene/scripts
idbin=$HOME/bio/idba/bin2
export PATH=$idbin:$PATH

if [ "X" = "X$inpe" ]; then echo "missing env inpe=readpairs.fa"; exit -1; fi
if [ "X" = "X$datad" ]; then echo "missing env datad=/path/to/data "; exit -1; fi
if [ "X" = "X$name" ]; then name=`basename $inpe .fa | sed 's/^/tridba/'`; fi
outdir=tidb$name$ov
cd $datad/

echo "#START  `date`"
echo "$idbin/idba_tran $traopts --read $inpe --num_threads $ncpu --out $outdir"
$idbin/idba_tran $traopts --read $inpe --num_threads $ncpu --out $outdir
echo "#DONE : `date`"

#................
# IDBA-Tran - Iterative de Bruijn Graph Assembler for next-generation transcriptome sequencing data.
# Usage: idba_tran -r read.fa -o output_dir
# Allowed Options: 
#   -o, --out arg (=out)                   output directory
#   -r, --read arg                         fasta read file (<=128)
#   -l, --long_read arg                    fasta long read file (>128)
#       --mink arg (=20)                   minimum k value (<=124)
#       --maxk arg (=60)                   maximum k value (<=124)
#       --step arg (=10)                   increment of k-mer of each iteration
#       --inner_mink arg (=10)             inner minimum k value
#       --inner_step arg (=5)              inner increment of k-mer
#       --prefix arg (=3)                  prefix length used to build sub k-mer table
#       --min_count arg (=2)               minimum multiplicity for filtering k-mer when building the graph
#       --min_support arg (=1)             minimum supoort in each iteration
#       --num_threads arg (=0)             number of threads
#       --seed_kmer arg (=30)              seed kmer size for alignment
#       --min_contig arg (=200)            minimum size of contig
#       --similar arg (=0.95)              similarity for alignment
#       --max_mismatch arg (=3)            max mismatch of error correction
#       --no_local                         do not use local assembly
#       --no_coverage                      do not iterate on coverage
#       --no_correct                       do not do correction
#       --pre_correction                   perform pre-correction before assembly
#       --max_isoforms arg (=3)            maximum number of isoforms
#       --max_component_size arg (=30)     maximum size of components
