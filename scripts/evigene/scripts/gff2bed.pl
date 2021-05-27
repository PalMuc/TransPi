
# gff2bed.pl : overlapping gff locs to coverage intervals, non-overlapped, bed format

gzcat est.leaf.mars11.gff.gz | env dolog=1 perl -ne\
' ($r,$b,$e,$c,$o)=(split"\t")[0,3,4,5,6]; $w=1+$e-$b; $li=$b-$lb; \
if($lr ne $r) { putc($lr,$lb,$le) if($lr and $lb>0); $lie=$lb=$le=0;  } \
if($lr eq $r and $li>0) { putc($lr,$lb,$li); } \
for($i=0; $i<$w; $i++){ $c[$i]+=$c; } \
($lr,$lb,$le,$lc,$lo)=($r,$b,$e,$c,$o); \
sub putc{ my($r,$ib,$w)=@_;  my $ie=$ib; my $lc=0; \
for my $i (0..$w) { my $c= shift @c; $lc=$c if($i==0); if($c != $lc or $i==$w) { \
$ib=$lie if($ib<$lie); $ie=$ib+1 if($ie<=$ib); \
if($lc>0){ $lc=sprintf("%.2f",log(1+$lc)) if($dolog); print join("\t",$r,$ib,$ie,$lc)."\n"; } \
$lie=$ib=$ie; } $ie++; $lc=$c; } } BEGIN{$dolog=$ENV{dolog};} ' \
