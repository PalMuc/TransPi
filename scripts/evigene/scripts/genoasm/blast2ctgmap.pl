#!/usr/bin/env perl
# blast2ctgmap.pl
## variants for scarpa, sopra scaffolders
## in: contigs_scafs.blastns 

# ** ASSUMES input sorted by scaffold loc ..
# BEGIN
##$INBLAST=1; ## agp format.. $sr,$sb,$se,$cd,$cb,$ce,$co
$INAGP=$ENV{agp}||0;
$NORND=$ENV{norand}||0; 
$NJOIN=$ENV{join}||1; 
$NREAD=$ENV{nread}||1; $RDLEN=50;
$MO=$ENV{mo}||200; 
$NOSKIP=$ENV{noskip}||0; $BADSKIP= not $NOSKIP;
$OFORM=$ENV{form}||0; ## readmap, ..
$DOSAME=$ENV{csame}||0; # $DOSAME=1 if($OFORM=~/readmap/);
## fixme: scarpa contigs.map output needs many same-contig pairs to calculate lib support..

while(<>) {
  next if(/^\W/); my @v;
  if($INAGP) { @v=split; } else { @bv=split;  @v= @bv[1,8,9,0,6,7]; }
  my($sr,$sb,$se,$cd,$cb,$ce,$co)=@v;  $co ||="+"; 
  if($cb>$ce) { ($cb,$ce)=($ce,$cb); $co="-"; } 
  elsif($sb>$se) { ($sb,$se)=($se,$sb); $co="-"; }
  @v[1,2,4,5,6]=($sb,$se,$cb,$ce,$co); 
  $ok=1; 
  if($lcd and $lcd ne $cd) {
    if($NJOIN>1 and $llcd) { $ok=putc($llcd,$cd,\@llv,\@v); }
    putcsame($lcd,\@lv) if($DOSAME);
    $ok=putc($lcd,$cd,\@lv,\@v); 
  }
  unless($ok<0) { @llv=@lv; $llcd=$lcd; @lv=@v; $lcd=$cd; } 
}

# END{ $NOputc=($lcd,$cd); } 

sub putcsame { ## same-contig pairs for scarpa contigs.map
  my($lcd,$lv)=@_;
  #my($lcd,$cd,$lv,$vv)=@_; 
  my @lv=@$lv; # my @v=@$vv;
  my($lsr,$lsb,$lse,$lcd,$lcb,$lce,$lco)=@lv; 
  my($sr,$sb,$se,$cd,$cb,$ce,$co)=@lv; ## for balance
  # ($cd,$co)=($lcd,$lco); ## dont need also
  my $ds=0;
  my $soff=0; 
  $soff= $NREAD * $MO; 
  my $plc= $lcb + $MO; #  $soff; $plc=$lcb if($plc<$lcb);
  my $pcb= $plc + $soff;
  my $nreadsame=9999999; # cover all of contig for this? or $NREAD
  for(my $ir=0; $ir<$nreadsame; $ir++) {
    my $lm= moff($MO,$plc,1,$lce); 
    my $nm= moff($MO,$pcb,1,$ce); # ce == $lce
    my $dc= $ds + ($nm-$lm); # ($nm-0) + ($lce-$lm); 
    if($OFORM eq "readmap") { # for scarpa contigs.map
      my($ro,$cn);
      ($cn)= $lcd=~m/0*(\d+)/;  $ro=2; #? ($lco eq "-")?2:1; ## ro=2 or lco == co, use same as below?
      print join(" ",++$RID,$ro,$cn,$lm,$RDLEN,0)."\n";
      ($cn)= $cd=~m/0*(\d+)/;   $ro=1; #? ($co eq "-")?1:2; ## ro=1 ?? or as below
      print join(" ",++$RID,$ro,$cn,$nm,$RDLEN,0)."\n";
    } else {
      print join("\t",$lcd,$lm,$lco,$cd,$nm,$co,$dc)."\n";     
    }
    $plc+= 101; $pcb+= 101; # next pair offset
    last if($plc>=$lce or $pcb >= $ce);
  } return 1; 
}

sub putc { 
  my($lcd,$cd,$lv,$vv)=@_; my @lv=@$lv; my @v=@$vv;
  my($lsr,$lsb,$lse,$lcd,$lcb,$lce,$lco)=@lv; 
  my($sr,$sb,$se,$cd,$cb,$ce,$co)=@v; 
  my $ds=$sb - $lse; 
  if($BADSKIP and $ds < -$MO) { if($se-$sb < 2*$MO) { return -1; } else { return 0; } } 
  $ds=0 if($ds<0); #??
  my $soff=0; 
  $soff= $NREAD * $MO; 
  my $plc= $lce - $soff; $plc=$lcb if($plc<$lcb);
  my $pcb= $cb;  
  if(0) { ## dont ^^ change for rev lco, co ?
      if($lco eq "-") { $plc= $lcb; }
      if($co eq "-") { $pcb= $ce - $soff; $pcb=$cb if($cb > $pcb); }
  }
  for(my $ir=0; $ir<$NREAD; $ir++) {
    my $lm= moff($MO,$plc,1,$lce); # was moff($MO,$lce,-1,$lcb);
    my $nm= moff($MO,$pcb,1,$ce); 
    my $lmdist= $lce - $lm;
    my $nmdist= $nm;
    my $dc= $ds + $nmdist + $lmdist; # ($nm-0) + ($lce-$lm); 
    ## FIXME ^^ maybe offsets depend on lco/co orient !*??, not ro below
    if(1) { ## this helps, is it right way?
      if($lco eq "-") { $lm= moff($MO,$lce-$plc,1,$lce); $lmdist=$lm; }
      if($co eq "-") { $nm= moff($MO,$ce-$pcb,1,$ce); $nmdist=$ce - $nm; }
      $dc= $ds + $nmdist + $lmdist; # ($nm-0) + ($lce-$lm); 
    }
    
    if($OFORM eq "readmap") { # for scarpa contigs.map
      my($ro,$cn);
      ## ?? ro from lco/co may be bad .. should be 2 = first, 1 = 2nd in path ??
      ($cn)= $lcd=~m/0*(\d+)/; $ro=($lco eq "-")?1:2; # this is right; ;;;;v of read2 ; BAD ??
      print join(" ",++$RID,$ro,$cn,$lm,$RDLEN,0)."\n";
      ($cn)= $cd=~m/0*(\d+)/;  $ro=($co eq "-")?2:1; # NOTE: 1,2 vs read1 ?? this works
      print join(" ",++$RID,$ro,$cn,$nm,$RDLEN,0)."\n";
    } else {
      print join("\t",$lcd,$lm,$lco,$cd,$nm,$co,$dc)."\n";     
    }
    $plc+= 101; $pcb+= 101; # next pair offset
    last if($lm>=$lce or $nm>=$ce);
  } return 1; 
} 

sub moff{ my($MOFF,$m,$d,$e)=@_; 
  # my $r=($NORND)?1:rand(1); ## too much rand spread
  my $r=($NORND)?1: 1 - rand(0.15);   
  my $om= $m + $d*int($MOFF*$r); # 100 + 
  $om=$e if($d>0 and $om>$e); $om=$e if($d<0 and $om<$e); 
  return $om; 
} 


__END__

env norand=1 join=2  mo=500 skip=1 perl -ne \
'@bv=split; ($sr,$sb,$se,$cd,$cb,$ce,$co)=@v=@bv[1,8,9,0,6,7]; $co="+"; if($cb>$ce) { ($cb,$ce)=($ce,$cb); $co="-"; } @v[4,5,6]=($cb,$ce,$co); $ok=1; 
if($lcd and $lcd ne $cd) {
  if($NJOIN>1 and $llcd) { $ok=putc($llcd,$cd,\@llv,\@v); }
  $ok=putc($lcd,$cd,\@lv,\@v); 
}
unless($ok<0) { @llv=@lv; $llcd=$lcd; @lv=@v; $lcd=$cd; } 
END{ $NOputc=($lcd,$cd); } 
sub putc { my($lcd,$cd,$lv,$vv)=@_; my @lv=@$lv; my @v=@$vv;
my($lsr,$lsb,$lse,$lcd,$lcb,$lce,$lco)=@lv; $ds=$sb-$lse; 
if($BADSKIP and $ds < -$MO) { if($se-$sb < 2*$MO) { return -1; } else { return 0; } } 
$ds=0 if($ds<0); $lm=moff($lce,-1,$lcb); $nm= moff($cb,1,$ce); 
$dc=$ds + ($nm-0) + ($lce-$lm); print join("\t",$lcd,$lm,$lco,$cd,$nm,$co,$dc)."\n"; return $dc; } 
sub moff{ my($m,$d,$e)=@_; my $r=($NORND)?1:rand(1); my $om= $m + $d*int(100 + $MO*$r); $om=$e if($d>0 and $om>$e); $om=$e if($d<0 and $om<$e); return $om; } 
BEGIN{ $NORND=$ENV{norand}||0; $NJOIN=$ENV{join}||1; $MO=$ENV{mo}||200; $BADSKIP=$ENV{skip}||0; }' contigs_scafs.blastns | head

