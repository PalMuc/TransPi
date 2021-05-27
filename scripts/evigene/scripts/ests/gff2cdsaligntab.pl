#!/usr/bin/env perl
# gff2cdsaligntab.pl : tabulate alignment annots from gmap.gff and gsplign.gff
# special cases CDS exon spans > mRNA span fix
#..............
# $nam.align.tab columns
# GenomeID        gespan  geor    AQueryID        quspan  match   qlen    cov     pid     path    indels  nexon 
#	splice  aalen   offs    aamap   sense   oid     tag

# $dosplit=1; 
$OIDisID=$ENV{oid}||0;
$CDSt=$ENV{cds}||"CDS";
$SenseMRNA=1; $SenseMRNA=0 if($ENV{notmrna}); ## #?? default yes? turn off {notmrna} ??

# BEGIN{ 
@k=qw(match qlen cov pid path indels nexon splice aalen offs aamap sense Split oid tag); 
@gs=qw(gescore clen gaps chimera); @hd=grep{ not/Split/ } @k;
print join("\t","GenomeID","gespan","geor","AQueryID","quspan",@hd)."\n"; 
# } 

while(<>) {
chomp; s/^#s\.//; @v=split"\t"; 
($gid,$src,$typ,$gb,$ge,$gv,$go,$at)=@v[0,1,2,3,4,5,6,-1];
if($typ eq $CDSt) { $cw=1+$ge-$gb; $cdsw+=$cw; $cgb=$gb if($cgb==0 or $gb<$cgb); $cge=$ge if($ge>$cge); }
elsif($typ eq "intron"){ $inspl += ($gv>69)?2:($gv>44)?1:0; }
elsif($typ eq "mRNA") { 
$issplit=(/Split=|chimera=/)?1:0; 
($tid,$tb,$te)=m/(?:trg|Target)=(\S+) (\d+) (\d+)/; 
if($ltid){ putg(); } 
$inspl=$cdsw=0; $cgb=$cge=0; # new mRNA span
%at= map{ $v=($at=~m/\b$_=([^;\s]+)/)?$1:0; $_ => $v; } ("ID",@k);
#bad.here# ($tid,$tb,$te)=m/(?:trg|Target)=(\S+) (\d+) (\d+)/; 
$id=$at{ID}; $tid ||= $id;  # $at{Split}=~s/^C//;
%ags= map{ if($at=~m/\b$_=([^;\s]+)/) { $_=>$1; } } @gs;
if(/gescore=/) {  $tag="gspl";
  if($at{splice}) { $at{splice}= 2 + int($at{splice}/2); }
  $at{qlen}= $qlen = $ags{clen}||0; $gaps=$ags{gaps}||0;
  ($pcov,$match,$ql2)= $at{cov} =~ m/(\d+).,(\d+).(\d+)/;
  $at{pid}||=99; $at{cov} = $pcov; $at{match}= $match;
} else { $tag="gmap";
  $aamap= $at{aamap}||$at{aalen}||0;
  ($chi)= $id =~ m/_C(\d)$/; if($ags{chimera} or $chi) { $chi||=1; $at{Split}="$chi/2" unless($at{Split}); }
  if($at =~ /aalen=(\d+,\d+[^;\s]+)/) { $aaq=$1; $at{aalen}=$aaq; $at{aamap}=$aamap if($aamap ne $aaq); }
} 
$at{tag}=$tag; ($ltid,$ltb,$lte)=($tid,$tb,$te);
($lgid,$lgb,$lge,$lgv,$lgo,$lat)=($gid,$gb,$ge,$gv,$go,$at);
} 
}

#END{ 
$issplit=0; putg(); 
#}

sub putg { 
if($cge) { $lgb=$cgb; $lge=$cge; }
if($SenseMRNA and $at{sense} < 0) { $lgo=($lgo eq "-")?"+":($lgo eq "+")?"-":$lgo; }
if($issplit and $ltid eq $tid) { %lat=%at; $cov1=$at{cov}; 
  if($gid ne $lgid) { $lat{Split}="C3:$lgid:$lgb-$lge,$lgo,$cov1"; }
  elsif($lge < $gb-1000 or $lgb > $ge+1000) { $lat{Split}="C2:$lgb-$lge,$lgo,$cov1"; } 
  else { $lat{Split}="C1:$lgb-$lge,$lgo,$cov1"; } return; }
elsif($at{Split} and $lat{Split}) { %att=%lat; 
  map{ $att{$_}+=$at{$_}; } qw(cov match nexon aamap); %at=%att; %lat=(); } 
if($at{Split}) { $at{path}= $at{Split}; }
$at{splice}= 2+$inspl if($inspl>0 and not $at{splice});
$at{aamap}= int($cdsw/3) if($cdsw>0 and not $at{aamap});
$lpid=($OIDisID)?$ltid:$at{ID}||$ltid; $lpid=~s/_[CG]\d+$//;
@at= @at{@hd}; print join("\t",$lgid,"$lgb-$lge",$lgo,$lpid,"$ltb-$lte",@at)."\n"; 
} 

## sort not required ..
# | sort -k1,1 -k2,2n -k4,4 > $nam.align.tab


