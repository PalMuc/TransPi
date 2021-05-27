#!/bin/bash
# velprepseq.sh  bamfile(s)
# .. inputs (a) mapped reads in .bam for now, later (b) unmapped fastq .. simpler
# ** add single reads (not properpairs);
# ** add duplseq reduction, maxdup=5? 8?; for pairs, lseq.rseq ; ?? use samtools rmdup insort.bam dedup.bam

dedup=1
minlen=35
# qualtype=sanger : for .sam/.bam
#.. sickle opts:
# -l, --length-threshold, Threshold to keep a read based on length after trimming. Default 20.
# -q, --qual-threshold, Threshold for trimming based on average quality in a window. Default 20.
# -x, --no-fiveprime, Don't do five prime trimming.
# -n, --discard-n, Discard sequences with any Ns in them.
# --quiet, do not output trimming info .. want

bindir=/bio/bio-grid/mb/rnaseq/bin
scripts=/bio/bio-grid/mb/evigene/scripts/rnaseq

if ! test -x ${bindir}/samtools; then echo "ERROR: missing bindir/samtools"; exit -1; fi
if ! test -x ${scripts}/velfq2fa.pl ; then echo "ERROR: missing scripts/velfq2fa.pl"; exit -1; fi

for bamf in $* ; do {

err=0
nam=`basename $bamf .bam`
echo "# making velvet inputs: trim.$nam.fa2 trim.$nam.fa1"

# properly-paired (0x2), no secondary 0x100 ; no dupls? 0x400; then sort by name -n
if ! test -s $nam-np.bam ; then 
  if [ $dedup = 1 ]; then
   ${bindir}/samtools view -u -f 0x2 -F 0x100 $bamf | ${bindir}/samtools rmdup - - |\
   ${bindir}/samtools sort -n - $nam-np ;
  else
   ${bindir}/samtools view -u -f 0x2 -F 0x100 $bamf | ${bindir}/samtools sort -n - $nam-np ;
  fi
  if ! test -s $nam-np.bam ; then err=1; exit $err; fi
fi

if ! test -s $nam.fq2 ; then 
# FIXME: this drops improper paired seqs now; save to $nam.fq0 ..
# FIXME: add location from sam loc=scaf:begin
  ${bindir}/samtools view $nam-np.bam | env nam=$nam perl -ne \
'sub revc{ my $c=shift; $c=reverse($c); $c=~tr/acgtACGT/tgcaTGCA/; return $c; }
BEGIN{ $nam=$ENV{nam}; open(F,">$nam.fq1"); open(R,">$nam.fq2"); open(S,">$nam.fq0");} 
($d,$fl,$chr,$cb,$s,$q)=(split"\t")[0,1,2,3,9,10]; $t=($fl & 0x80)?2:1;  
if($fl & 0x10) { $s=revc($s); $q=reverse($q); } 
if($d eq $ld) { if($lt == 2){ @sw=($ld,$lt,$lchr,$lcb,$ls,$lq); 
($ld,$lt,$lchr,$lcb,$ls,$lq)=($d,$t,$chr,$cb,$s,$q); ($d,$t,$chr,$cb,$s,$q)=@sw; }
print F join("\n","\@$ld/$lt loc=$lchr:$lcb",$ls,"+",$lq)."\n"; 
print R join("\n","\@$d/$t loc=$chr:$cb",$s,"+",$q)."\n"; 
$d=$s=$ld=$ls=""; } elsif($ld and $ls) { 
print S join("\n","\@$ld/$lt loc=$lchr:$lcb",$ls,"+",$lq)."\n"; $ld=$ls=""; }
($ld,$lt,$lchr,$lcb,$ls,$lq)=($d,$t,$chr,$cb,$s,$q);' 
  
  if ! test -s $nam.fq2 ; then err=2; exit $err; fi
fi

if ! test -s trim.$nam.fq2 ; then 
  ${bindir}/sickle pe -l $minlen -t sanger -f $nam.fq1 -r $nam.fq2 \
    -o trim.$nam.fq1 -p trim.$nam.fq2 -s trim.$nam.fq0 
  if ! test -s trim.$nam.fq2 ; then err=3; exit $err; fi
fi

# unpaired reads
if ! test -s $nam.fqs ; then 
  if [ $dedup = 1 ]; then
   ${bindir}/samtools view -u -F 0x102 $bamf | ${bindir}/samtools rmdup - - |\
   ${bindir}/samtools view - | env nam=$nam perl -ne \
'BEGIN{open(F,">$ENV{nam}.fqs")} ($d,$s,$q)=(split"\t")[0,9,10]; print F join("\n","\@$d",$s,"+",$q)."\n";'
  else
   ${bindir}/samtools view -F 0x102 $bamf | env nam=$nam perl -ne \
'BEGIN{open(F,">$ENV{nam}.fqs")} ($d,$s,$q)=(split"\t")[0,9,10]; print F join("\n","\@$d",$s,"+",$q)."\n";'
  fi
  # combine singles from above
  if test \( -s $nam.fq0 -a -s $nam.fqs \); then cat $nam.fq0 >> $nam.fqs; fi
  #NO# if ! test -s $nam.fqs ; then err=4; exit $err; fi
fi

if ! test -s trim.$nam.fqs ; then
  ${bindir}/sickle se -l $minlen -t sanger -f $nam.fqs -o trim.$nam.fqs
fi


# paste .fq1,fq2 to .fa2,.fa1 for velvet ; add -output name or pipe? add -pair/-single flags for input type
${scripts}/velfq2fa.pl trim.$nam.fq1 trim.$nam.fq2  > trim.$nam.fa2
${scripts}/velfq2fa.pl trim.$nam.fq0 > trim.$nam.fa1
${scripts}/velfq2fa.pl trim.$nam.fqs >> trim.$nam.fa1
if ! test -s trim.$nam.fa2 ; then err=5; exit $err; fi

echo "# velvet inputs: "; ls -l trim.$nam.fa[12]
gzip --fast trim.$nam.fa[12]
if [ $err == 0 ]; then echo "/bin/rm $nam-np.bam $nam.fq[012s] trim.$nam.fq[012s] "; fi


} done

#---------------------------
# velmakeseq4b.sh
# #!/bin/bash
# ##  qsub -q batch velmakeseq4b.sh
# #PBS -N velmakeseq
# #PBS -l nodes=1:ppn=16,walltime=11:55:00
# #PBS -o velmakeseq.$$.out
# #PBS -e velmakeseq.$$.err
# #PBS -V
# 
# export velbin=$HOME/bio/velvet/bin
# export scripts=$HOME/bio/evigene/scripts/rnaseq
# export bindir=$HOME/bio/bin
# 
# workd=$HOME/scratch/chrs/cacao/rnas
# bamdir=$workd/bamsc
# outdir=$workd/vel4bfa
# 
# scset="10r 9 8 7 6"
# # scset="5 4 3 2 1 x"
# sgrp="cgb nc1 nc2"
# 
# cd $workd
## .. tmpdir and outdir? for final files.
# mkdir $outdir
# cp -p longifa2.idlist $outdir/
# cd $outdir
# 
# for i in $scset; do {
#   sc="sc$i"
#   for grp in $sgrp; do { 
#     bamf=$bamdir/$sc.$grp.bam
#     nam=`basename $bamf .bam`
#     if [ $grp = "cgb" ]; then
#      ( $workd/velprepseq.sh $bamf; ) > log.velprep$nam 2>&1 &
#     else
#      ( $workd/velprepseq.sh $bamf; $workd/insplitfa.sh trim.$nam.fa2.gz ) > log.velprep$nam 2>&1  &
#     fi
#   } done
# } done
# 
# wait
#---------------------------
