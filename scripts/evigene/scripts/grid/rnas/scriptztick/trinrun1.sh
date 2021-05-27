#! /bin/bash
### qsub -q normal trin1.sh
#PBS -A ind114
#PBS -N trin4d
#PBS -l nodes=1:ppn=16,walltime=29:55:00
#PBS -o trin4.$$.out
#PBS -e trin4.$$.err
#PBS -V

dv=1
pairmax=500
ncpu=7
mem=52G
## dv1 running ok in 50Gb, 22hr: 26k of 130k buttfly done (10hr for bf); likely outoftime at 29hr; need 39+hr
## is --bflyCPU 2 too low? maybe more, allow self-calc-cpu?
# gordon.sdsc
datad=$HOME/scratchg/
workd=$datad/chrs/aabugs/tsa/ztick
cd $workd

rund=$workd/trin
rdfile=$workd/sraf/allpe
oname=`basename $rdfile`
oname=trout$dv$oname

trindir=$HOME/bio/trinity
export PERL5LIB=$trindir/PerlLib:${PERL5LIB}

topts="--kmer_method jellyfish --bflyHeapSpaceMax 10G --bflyCPU 2 --CPU $ncpu --max_memory $mem --SS_lib_type RF " 
mkdir $rund
cd $rund

$trindir/Trinity.pl $topts --output $oname --seqType fq --left ${rdfile}_1.fq   --right ${rdfile}_2.fq


##... trin is chewing up mem again in Buttfly part, try dv=2 just buttfly w/ new opts
#    resources_used.cput = 67:08:36
#    resources_used.mem = 42609240kb     << but this is ok
#    resources_used.vmem = 172120548kb   << too much for 64G box?
#    resources_used.walltime = 12:52:19
##.. need -onlybuttfly option ..  reduce mem=36G  ??
## bfopt="--bflyHeapSpaceMax 10G --bflyCPU 2"  or --bflyCalculateCPU
##..
#  --max_number_of_paths_per_node <int>  :only most supported (N) paths are extended from node A->B,
#                                         mitigating combinatoric path explorations. (default: 10)
#  --group_pairs_distance <int>    :maximum length expected between fragment pairs (default: 500)
#                                   
#  --path_reinforcement_distance <int>   :minimum overlap of reads with growing transcript 
#                                        path (default: 75)
#
#  --lenient_path_extension        :require minimal read overlap to allow for path extensions. 
#                                   (equivalent to --path_reinforcement_distance=1)
#
#  --bflyHeapSpaceMax <string>     :java max heap space setting for butterfly
#                                   (default: 20G) => yields command
#                  'java -Xmx20G -jar Butterfly.jar ... $bfly_opts'
#  --bflyHeapSpaceInit <string>    :java initial hap space settings for
#                                   butterfly (default: 1G) => yields command
#                  'java -Xms1G -jar Butterfly.jar ... $bfly_opts'
#  --bflyGCThreads <int>           :threads for garbage collection
#                                   (default, not specified, so java decides)
#  --bflyCPU <int>                 :CPUs to use (default will be normal 
#                                   number of CPUs; e.g., 2)
#  --bflyCalculateCPU              :Calculate CPUs based on 80% of max_memory
#                                   divided by maxbflyHeapSpaceMax

