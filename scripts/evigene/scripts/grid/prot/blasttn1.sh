#! /bin/bash -l

SHARED=$HOME/scratch
ACCOUNT=TG-xxx

# ncbibin=$HOME/bio/ncbi/bin
ncbibin=/N/soft/linux-sles9-ppc64/ncbi-2.2.16/bin
aug=$HOME/bio/augustus
protdir=`pwd`

# should be params:
proteins=$protdir/dmoj_glean1_good.aa
protnam=dmojgln
faparts=1000

nadir=`pwd`/dgri3
naseq=dgri_caf060210.fa
# mkdir $nadir

if [ ! -f $nadir/$naseq.nsq ]; then
$ncbibin/formatdb -p F -i $nadir/$naseq
fi

if [ ! -f $protdir/${protnam}1.fsa ]; then
$aug/scripts/split_multifasta.pl --f=$protnam --seqs_per_file=$faparts \
  --output_dir=$protdir \
  --in $proteins
fi

qset=`/bin/ls -rt $protdir/${protnam}*.fsa`

for qfile in $qset
{
  query=`basename $qfile`
  jobn=runtbn-$query.sh
  
  cat > $jobn <<EOJ
#! /bin/bash -l
## IBM LoadLeveler job, llsubmit
# @ class = NORMAL
# @ account_no = $ACCOUNT
# @ wall_clock_limit = 04:50:00
# @ error   = $nadir/$jobn.err
# @ queue

$ncbibin/blastall -p tblastn -m 9 -e 0.00001 \
  -o $nadir/$query.tblastn \
  -d $nadir/$naseq  \
  -i $protdir/$query

EOJ
  
# llsubmit -q $jobn

}
