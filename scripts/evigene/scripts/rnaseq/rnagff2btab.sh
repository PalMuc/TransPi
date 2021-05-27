#!/bin/tcsh
# env gff=aphid_rnaseq.all17cuff8  rnagff2btab.sh

# use env maketr=1 # set MAKETR=0
if ( ! $?maketr ) then
 set maketr=0
endif
if ( ! $?feature ) then
  set feature=exon
endif

set SCOREISIDENT=0
set workd=/export/udisk2/work/nasv4
set dgenome=nasvit1asm
## see gmapgff2btab.pl for gmap.gff

foreach trf ( $gff )
  # set nam=`basename $trf .gff.gz`
  set nam=`echo $trf | sed 's/\..*$//;'`
  set CAT=cat; 
  #FIXME: if($trf =~ .gz) set CAT=gzcat

  # .btab for PASA; note exon calcs for spliced
  $CAT $trf | grep "	$feature" | env piscore=$SCOREISIDENT perl -ne \
  '($r,$s,$t,$b,$e,$v,$o,$px,$at)= @v=split"\t";  \
  $w=1+$e-$b; $at=~/Parent=([^;\s]+)/; $id=$1;  \
  ($b,$e)=($e,$b) if($o eq "-");  \
  if($lid eq $id) { $xi++; } else { bput(); @x=(); $ci++; $xi=1; }  \
  push(@x,[$r,$id,$b,$e,$v,$o,$w,$ci,$xi]);  $lo=$o; $lid=$id;  \
  END{ bput();} sub bput { my $rev=($lo eq "-")?1:0; my $x0=0; my $nx=@x; my $ix=0; \
  if($rev) { @x=sort{$b->[2]<=>$a->[2]} @x;} else { @x=sort{$a->[2]<=>$b->[2]} @x;}  \
  foreach $xr (@x) { my($r,$id,$b,$e,$v,$o,$w,$ci,$xi)=@$xr; \
  $br=1 + $x0; $er= $x0 + $w; $x0 += $w; $ix++; $xi=$ix; $pi=$ENV{piscore}?$v:100; \
  print join("\t",$r,"","","gmap","",$id,$b,$e,$br,$er,$pi,"",1,$ci,$xi),"\n"; } } '\
    > $nam.btab
## NOT @x=reverse @x if($rev);  

  if( $maketr != 0 ) then
    $workd/scripts/cdsgff2genbank.pl -MIN_CDS_LEN=1 -t=$feature -a dna -gff $trf -fasta $workd/genome/$dgenome.fa > $nam.btab.fa 
  endif

end
