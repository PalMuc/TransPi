#! /bin/bash -l

# note: this doesnt work from gg_job.pl because it doesn't use a parts.list of genome parts
# fix somehow? or not

SHARED=$HOME/scratch
ACCOUNT=TG-xxx

QTYPE=PBS
#CLASS=NORMAL
CLASS=dque
CLOCKTIME=11:55:00

ncbibin=$HOME/bio/ncbi/bin
aug=$HOME/bio/augustus

# should be params:
out_dir=`pwd`
prot_dir=`pwd`
proteins=$prot_dir/MY_GENES.aa
protnam=queryaa
faparts=1000

db_dir=$SHARED/ncbidb
# use blast in parts on nrdb to save mem
dbset="nr.00 nr.01 nr.02"

#----------------

## assume blast db is formatted.
# for db1 in $dbset
# {
# if [ ! -f $db_dir/$db1.psq ]; then
#   $ncbibin/formatdb -p T -i $db_dir/$db1
# fi
# }

if [ ! -f $prot_dir/${protnam}1.fsa ]; then
  $aug/scripts/split_multifasta.pl --f=$protnam --seqs_per_file=$faparts \
    --output_dir=$prot_dir \
    --in $proteins
fi

qset=`/bin/ls -rt $prot_dir/${protnam}*.fsa`

for qfile in $qset
{
  for db1 in $dbset
  {
  query=`basename $qfile`
  jobn=bl$db1-$query
    
if [ $QTYPE == "PBS" ]; then
    cat > $jobn <<EOJ
#! /bin/bash -l
### qsub -q $CLASS $jobn
#PBS -N $jobn
#PBS -A $ACCOUNT
#PBS -l nodes=1:ppn=1,walltime=$CLOCKTIME
#PBS -o $out_dir/$jobn.out
#PBS -e $out_dir/$jobn.err
#PBS -V

$ncbibin/blastall -p blastp -m 9 -e 0.001 \
  -o $out_dir/$jobn.blastp \
  -d $db_dir/$db1  \
  -i $prot_dir/$query
EOJ
     # debug echo on/off
     echo "qsub -q $CLASS $jobn"

elif [ $QTYPE == "LL" ]; then
    cat > $jobn <<EOJ
#! /bin/bash -l
## llsubmit -q $jobn
# @ class = $CLASS
# @ account_no = $ACCOUNT
# @ wall_clock_limit = $CLOCKTIME
# @ error   = $out_dir/$jobn.err
# @ queue

$ncbibin/blastall -p blastp -m 9 -e 0.001 \
  -o $out_dir/$jobn.blastp \
  -d $db_dir/$db1  \
  -i $prot_dir/$query
EOJ
  
   echo "llsubmit -q $jobn"
fi

  }
}  


# eoj


