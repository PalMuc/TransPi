#!/bin/bash
## env infa=xxx.fasta.keep datad=`pwd` qsub -q normal faspot2pair.sh
#PBS -N faspot2pair 
#PBS -A ind114
#PBS -l nodes=1:ppn=8,walltime=1:55:00
#PBS -V
## see also scriptdax/runsra2fa.sh

if [ "X" = "X$datad" ]; then echo "err missing datad=path/to/data"; exit -1; fi
if [ "X" = "X$infa" ]; then echo "err missing infa=xxx*.fasta"; exit -1; fi
if [ "X" = "X$nskip" ]; then nskip=0; fi
if [ "X" = "X$revcomp" ]; then revcomp=0; fi

spofa=$infa
# mkdir pairfa

i=0; for sfa in $spofa; do { 
  name=`basename $sfa .fasta | sed 's/\.keep//; s/\.fa.*//; s/$/dn/; '`
  if [ $revcomp -gt 0 ]; then name="${name}rc"; fi
  fa1=${name}_1.fa
  if [ -s $fa1 ]; then echo "exists: $fa1"; continue; fi

  env fn=$name rc=$revcomp nn=$nskip perl -ne \
'BEGIN{ $NN=$ENV{nn}; $RC=$ENV{rc}; $nok=$nsk=0;
 $f=$ENV{fn}; open(L,">${f}_1.fa"); open(R,">${f}_2.fa"); } 
if(/^>(\S+)/) { $d=$1; } else { chomp; $hl=int(length($_)/2); 
$sl=substr($_,0,$hl); $sr=substr($_,$hl);
if($NN>0) { $nl= $sl=~tr/Nn/Nn/; $nr= $sr=~tr/Nn/Nn/; if($nl>$NN or $nr>$NN) { $nsk++; next;} } 
if($RC) { $sl=revc($sl); $sr=revc($sr); }
print L ">$d/1\n",$sl,"\n"; print R ">$d/2\n",$sr,"\n"; $nok++; }
END{ warn "#spot2pair: $ENV{fn} nok=$nok, nskip=$nsk \n"; }
sub revc{ my($s)=$_[0]; $s=reverse($s); $s=~tr/ACGT/TGCA/; $s; }' $sfa   &

  i=$(( $i + 1 )); if [ $i -ge $ncpu ]; then wait; i=0; fi

} done
wait

exit;

#-----------
## add opt: revcomp for velvet -strand_specific? wants FR instead of RF orient of stranded pairs
## add opt: skip spot/pairs if NNN content of either part is high (NNN > cutoff?)

# gunzip -c spotfa/SRR514101.fasta.gz |  env fn=cornlo1rc  perl -ne  \
# 'BEGIN{ $f=$ENV{fn}; open(L,">${f}_1.fa"); open(R,">${f}_2.fa"); } \
# if(/^>(\S+)/) { $d=$1; } else { chomp; $hl=int(length($_)/2); \
# print L ">$d/1\n",revc(substr($_,0,$hl)),"\n"; print R ">$d/2\n",revc(substr($_,$hl)),"\n"; } \
# sub revc{ my($s)=$_[0]; $s=reverse($s); $s=~tr/ACGT/TGCA/; $s; } ' 

