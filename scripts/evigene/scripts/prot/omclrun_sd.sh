#! /bin/bash
### qsub -q shared omclrun.sh
#PBS -N omclmake
#... ?? PBS -A ddp138
#PBS -A ind114
#PBS -l nodes=1:ppn=1,walltime=7:55:00
#PBS -o omcl1.$$.out
#PBS -e omcl1.$$.err
#PBS -V

ncpu=1

datad=/phase1/$USER
workd=$datad/chrs/nasv1
bindir=$datad/bio/bin
orlib=$datad/bio/orthomcl3
mcl=$datad/bio/mcl9/bin/mcl
evigene=$datad/bio/evigene

rund=/scratch/$USER/$PBS_JOBID

fasta=$workd/prot/arp11.aa
bldir=$workd/prot/arp11bl
onam=arp11u_mcl
indir=$workd/prot/omcl9in
outdir=$workd/prot/omcl9u4
notef=$workd/prot/omcl9.$$.RUNNING

touch $notef
echo "START " >> $notef
echo `date`   >> $notef
echo "rundir=$rund"   >> $notef

mkdir $outdir
mkdir -p $rund
mkdir $rund/inblast
cd $rund/

## indir has result of blast92orthomcl10, if done before
cp -rp $indir/* $rund/

# #make: takes ~1hr
cp -p $fasta $rund/infasta
cp -rp $bldir/* $rund/inblast/
$evigene/scripts/blast92orthomcl10.pl -fasta=infasta -in=inblast/ -out=$onam

if [ ! -s $onam.bpo  ]; then
  echo "FAIL: no $onam.bpo" >> $notef; exit -1
fi
if [ ! -s $onam.gg  ]; then
  echo "FAIL: no $onam.gg" >> $notef; exit -1
fi

# run:
export MCL=$mcl 
export ORTHOMCL=$rund/
perl -I$orlib $orlib/orthomcl.pl --mode 4 --bpo=$onam.bpo --gg=$onam.gg 

## problem : orthomcl makes own outdir name # exclude indir files... $onam.bpo $onam.gg inblast
du  >> $notef
ls -lR  >> $notef
rm -rf infasta inblast/
cp -rp $rund/* $outdir/

echo "DONE" >> $notef
echo `date`   >> $notef
