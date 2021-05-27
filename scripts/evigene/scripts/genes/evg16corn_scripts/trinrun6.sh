#! /bin/bash
### env name=itick inpe=allpe_1.fa datad=`pwd` qsub -q normal trin1.sh
#... PBS -A ind114
#PBS -N trin4
#PBS -l vmem=255gb,nodes=1:ppn=16,walltime=27:55:00
#PBS -V

dv=6s
if [ "X" = "X$ncpu" ]; then ncpu=8; fi
if [ "X" = "X$mem" ]; then mem=120G; fi
if [ "X" = "X$datad" ]; then echo "err missing datad=path/to/data"; exit -1; fi
if [ "X" = "X$inpe" ]; then echo "err missing inpe=data_1.fa"; exit -1; fi
if [ "X" = "X$name" ]; then name=myspp; fi

# module at comet.sdsc
module add trinity

cd $datad

## newer trinity allows multifile --left infa1,infa1b --right infa2,infa2b ..
#o#inpe=`echo $inpe | sed "s,sraf,$datad/sraf,g;"`
inpe=`echo $inpe | sed "s,sraf,$datad/sraf,g; s/ /,/g;"`
infa1=$inpe
infa2=`echo $inpe | sed 's/_1/_2/g;'`
oname=trin$name$dv

#oldop="--kmer_method jellyfish --bflyHeapSpaceMax 10G --bflyCPU 2 --CPU $ncpu --max_memory $mem --SS_lib_type RF" 
#o# topts="--SS_lib_type RF --JM $mem --CPU $ncpu --bflyHeapSpaceMax 10G"
#new.old changed again# topts="--CPU $ncpu --max_memory $mem"
# #  --JM <string>            :(Jellyfish Memory) number of GB of system memory to use for 
#                            k-mer counting by jellyfish  (eg. 10G) *include the 'G' char 
## NOTE JM IS NOT maxmem, but mem/forked-job, eg. maxmem/ncpu or less, eg 2G
#. JM=8G
#. topts="--CPU $ncpu --JM $JM"
## dang mem param changed again !!
#  --max_memory <string>      :suggested max memory to use by Trinity where limiting can be enabled. (jellyfish, sorting, etc)
#                            provied in Gb of RAM, ie.  '--max_memory 10G'

JM=10G
topts="--CPU $ncpu --max_memory $JM"
# if [ $STRANDS ]; then
topts="$topts --SS_lib_type RF"

echo Trinity $topts --seqType fa --left $infa1 --right $infa2
Trinity $topts --seqType fa --left $infa1 --right $infa2

## sample: # Trinity --seqType fq --max_memory 50G --left reads_1.fq  --right reads_2.fq --CPU 6

