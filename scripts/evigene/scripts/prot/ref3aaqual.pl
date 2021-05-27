#!/usr/bin/env perl
# ref3aaqual.pl
# /bio/bio-grid/aabugs4/bugs/omclbugs10/ref3aaqual
# rewrite sumqual2.sh to add stats: sd, sem, t-test for group compares (versions)

my $debug=$ENV{debug}||0;
my $REFAV=  $ENV{refaa}||560; #  ref ave aasize
# upd: let refaa=file sizetab get REFAV from ave in that
my $NT= $ENV{nt}||5173; # option
my $hd= $ENV{hd}||0;
my $USE2= $ENV{two} || 0;
my $SHOWTOTAL= $ENV{total} || 0;
my $SHOWSD= $ENV{sd} || 0;
my $SHOWSEM= $ENV{sem} || 0;
my $namwid= $ENV{wid} || 15;
my $TOPHIT= $ENV{tophit} || $ENV{top} ||0;
my $swapqt= $ENV{swap} || 0;
my $hasaaqual= 0;
my $hasrefid= 0;
my @inf= @ARGV;
my(%ns, %refid, %did, %svs, %ssq);

## Dang mixed ref3 idprefixes: 's/daphniaplx/daphniap/' or equivalence..
sub fixidp { local $_= shift; if(/daphniap/) { s/daphniaplx:/daphniap:/; } return $_; }

if($REFAV =~ /^\D/ and -f $REFAV) {
  my($na,$sa)=(0,0); 
  open(F,$REFAV); while(<F>) { my($id,$aw)=split; next if($id=~/^total|^\W/); 
    if($aw>0) { $id=fixidp($id); $refid{$id}=$aw;  $na++; $sa+=$aw; } } close(F);
  if($na>1) { $hasrefid=$na; $REFAV= int($sa/$na); $NT=$na unless($ENV{nt}); } #? reset NT or not
}

foreach my $inf (@inf) {
  open(F,$inf) or die $inf;
  $nam= $inf; 
  $nam =~ s,^.*/([^/]+)$,$1,; 
  $nam =~ s/ref\d*[_-]//;
  $nam =~ s/\.(tbest\w*|qual\w*)$//; # $nam =~ s/\..*$//; #?
  $nam=substr($nam,0,$namwid);  
  $hd++; %ns= %did= %svs= %ssq= ();
  while(<F>) {
    next unless(/^\w/); 
    # fixme, input aa.qual tab instead? sg, sl, gaps, aq, clen; no qg, bs, id, al, ql
    ($qg,$sg,@v)=split;  ## qg == ref gene expected, use swapqt=1 if needbe
    ($bs,$id,$al,$ql,$sl, $aq, $clen)=@v; 
    if($swapqt) { ($qg,$sg,$ql,$sl)=($sg,$qg,$sl,$ql); @v[3,4]=($ql,$sl); }
    $qg= fixidp($qg); # dang..
    ($qsp)= $qg=~m/^(\w+)/;   # qsq here is ref species, depending on options, use all or skip best match

    next if($TOPHIT and ($did{$qg} or $did{$qsp.$sg})); #?? OPTION for .tall4 input
    
    # next if($hasrefid and not $refid{$qg}); # TEST, want this or not?
    if($hasrefid) {
      next unless($refid{$qg}); # TEST, want this or not?
      $v[3]= $ql= $refid{$qg} if($ql<1);
    }

    # aq, clen: only for .qual set w/ prot quals added;
    ($al2,$ap,$ac)=split",",$aq; $cc=($ac=~/complete/)?4:($ac=~/partial(\d)/)?$1:0; ($uq)= $ac=~/(utr\w+)/;
    $hasaaqual++ if($ac);
    
    $n++; $ns{$qsp}++; $did{$qg}++; $did{$qsp.$sg}++;
    # add sum of squares for sd, ttest.  %svs == %sum; %ssq
    my $vaaf= ($cc == 4)? 100 : 0;
    my $vutr= ($uq =~ /utr(bad|poor)/) ? 0 :100;
    my @vx= (@v, $sl - $ql, $vaaf, $vutr, $clen); # last 3 for hasaaqual
    for $i (0..$#vx) { my $v=$vx[$i]; $svs{$qsp}[$i]+=$v; $ssq{$qsp}[$i]+= $v*$v; } # bit, idn, aln, lnq, lns
    
    #? want min,median,max values?
    
    # $svs{$qsp}[5] += $sl - $ql; # dln
    # $svs{$qsp}[6] += 100 if($cc == 4); # aafull
    # # $XXsvs{$qsp}[6] += 100 if($cc == 0);
    # $svs{$qsp}[7] += 100 unless($uq =~ /utr(bad|poor)/); #utrok
    # # $XXsvs{$qsp}[7] += 100 if($uq =~ /utrbad/);
    # $svs{$qsp}[8] += $clen; # trln
  } close(F);
  putrow(); # from %ns, %svs
}

sub putrow {

  my $K=6;  my $TK=6; 
  my @col=qw(bit idn aln lnq lns dln); # col[0..5] for Tcol
  if($hasaaqual) { $K += 3; push @col, qw( aafull utrok trln);  }
  
  my(@sv, @ssq, @av, @atv, @asd, @asdt);
  @ssq= @sv= @av= @avt= (0) x $K; 
  my @s= sort keys %ns; my $ns=0;  my $s1=$s[0]; my @s2=@s; 
  if($USE2) { 
    my %sbp; for my $s (@s) { $spb{$s} += $svs{$s}[2]; } # 2==aln OPTION here
    ($s1,@s2)= sort{$spb{$b}<=>$spb{$a}}keys %spb;
  }
  
  # my $NT3=$NT*scalar(@s2);# WRONG for hasrefid
  my $NT3=($hasrefid) ? $NT : $NT * scalar(@s2);

  for my $s (@s2) { $ns += $ns{$s}; 
    for my $i (0..$K-1) { $sv[$i] += $svs{$s}[$i]; $ssq[$i] += $ssq{$s}[$i]; } 
  } 
  
  for my $i (0..$K-1) { 
    my $sv= $sv[$i]; my $ssv= $ssq[$i];
    my($n,$mn,$sd)= meansd($ns,$sv,$ssv);
    $av[$i]= int(0.5 + $mn); 
    $asd[$i]= int(0.5 + $sd); #?
    # Totals, div by $NT3 not ns
    ## also Tlnq should == REFAV, so Tdln =~ Tlns - Tlnq, as dln = lns - lnq
    if($i == 5) { # Tdln, subtract ref ave aasize
      $sv= $sv[4]; $ssv= $ssq[4]; # wrong ssq for this calc; fix above?
      # $ssv -= $ns * $REFAV * $REFAV; #? is this right; probably not, may be neg ? remove  
      ($n,$mn,$sd)= meansd($NT3,$sv,$ssv);
      $mn -= $REFAV; # subtract ref ave aasize
    } else {
      ($n,$mn,$sd)= meansd($NT3,$sv,$ssv);
      if($i == 3) { $mn= $REFAV; } #? yes or no
    }
    $atv[$i]=int(0.5 + $mn); 
    $asdt[$i]= int(0.5 + $sd); #?
  } 
  ## $atv[5]= int(0.5 + $sv[4]/$NT3 - $REFAV); # Tdln, subtract ref ave aasize, FIXME opt
  
  if($SHOWTOTAL) {
    my @tcol=map{"T".$_} @col; # option:on or off : but not  aafull utrok trln
    if($SHOWTOTAL == 1) { @col=@av=@asd=(); } # show Total cols only; show both for total=2
    push @col, @tcol[0..$TK-1];
    push @av, @atv[0..$TK-1];
    push @asd, @asdt[0..$TK-1];
  }  
  my $nfmt="%-".$namwid."s";
  print join("\t",sprintf($nfmt,"group"),"ns",@col,"sppset")."\n" if($hd==1);  
  my $s2=join",",@s2; 
  print join("\t",sprintf($nfmt,$nam), $ns, @av, $s2)."\n"; 
  my $st=($SHOWSEM)?"SEM_":"SD_";
  print join("\t",sprintf($nfmt,$st.substr($nam,0,$namwid-4)), $ns, @asd, $s2)."\n" if($SHOWSEM or $SHOWSD); 
} 

sub meansd { # medmeansd
  my($n,$sm,$ss)= @_;
  return($n,$sm,0) if($n<2);
  my $mn= $sm/$n; 
  my $var= (($ss - $sm*$mn)/($n-1)); 
  my $sd= ($var>=1)?sqrt($var):0;
  $sd= $sd/sqrt($n) if($SHOWSEM); # SEM, option or only this?
  return($n,$mn,$sd); # $md,$min,$max);
}

# sub medmeansd {
#   my($arr)= @_;
#   my ($n,$sm,$ss,$md,$mn,$sd,$min,$max)=(0) x 9;
#   my @sl= sort{$b<=>$a} @$arr; # big first
#   $n=@sl;  ($max,$min)=  @sl[0,-1];
#   return($n,$md,$mn,$sd,$min,$max) if($n<1);
#   $md= @sl[int($n/2)]; 
#   return($n,$md,$mn,$sd,$min,$max) if($n<3);
#   for my $i (0..$n-1) { my $x= $sl[$i]; $sm+=$x; $ss += $x*$x; }
#   ##$mn= $sm/$n; $sd= sqrt(($ss - $sm*$mn)/($n-1));
#   $mn= $sm/$n; my $var= (($ss - $sm*$mn)/($n-1)); $sd=($var>=1)?sqrt($var):0;
#   return($n,$md,$mn,$sd,$min,$max);
# }

__END__

=item new output

flamingo2.% cat ref2ins3stats.tab1 | cut -f1-2,9-
group           ns      Tbit    Tidn    Taln    Tlnq    Tlns    Tdln    sppset
acypi1_2010     9540    353     224     402     542     507     -52     AMELL,daphniaplx
SEM_acypi1_2010 9540    4       2       4       6       5       5       AMELL,daphniaplx
acypi2_2011     9710    388     243     435     552     543     -16     AMELL,daphniaplx
SEM_acypi2_2011 9710    5       3       4       6       6       6       AMELL,daphniaplx
drosmel_r3_0    9779    378     240     429     554     597     37      AMELL,daphniaplx
SEM_drosmel_r3_ 9779    4       2       4       6       6       6       AMELL,daphniaplx
drosmel_r5_30   9887    391     248     444     558     627     67      AMELL,daphniaplx
SEM_drosmel_r5_ 9887    5       3       5       6       7       7       AMELL,daphniaplx
nasvit_ogs12    9765    418     256     434     538     574     14      TCAST,daphniaplx
SEM_nasvit_ogs1 9765    4       2       4       5       5       5       TCAST,daphniaplx
nvit2_evigenes_ 10080   435     269     458     554     584     24      TCAST,daphniaplx
SEM_nvit2_evige 10080   5       3       4       5       6       6       TCAST,daphniaplx

cat ref2ins3stats.tab2 | cut -f1-2,9-
group           ns      Tbit    Tidn    Taln    Tlnq    Tlns    Tdln    sppset
acypi1_2010     9540    353     224     402     542     507     -52     AMELL,daphniaplx
SD_acypi1_2010  9540    359     211     365     584     473     473     AMELL,daphniaplx
acypi2_2011     9710    388     243     435     552     543     -16     AMELL,daphniaplx
SD_acypi2_2011  9710    458     265     450     583     562     562     AMELL,daphniaplx
drosmel_r3_0    9779    378     240     429     554     597     37      AMELL,daphniaplx
SD_drosmel_r3_  9779    422     243     417     579     600     600     AMELL,daphniaplx
drosmel_r5_30   9887    391     248     444     558     627     67      AMELL,daphniaplx
SD_drosmel_r5_  9887    462     267     460     577     712     712     AMELL,daphniaplx
nasvit_ogs12    9765    418     256     434     538     574     14      TCAST,daphniaplx
SD_nasvit_ogs1  9765    432     246     403     550     530     530     TCAST,daphniaplx
nvit2_evigenes_ 10080   435     269     458     554     584     24      TCAST,daphniaplx
SD_nvit2_evige  10080   473     272     440     549     570     570     TCAST,daphniaplx

ref3ins3stats.tab1
group           ns      bit     idn     aln     lnq     lns     dln     Tbit    Tidn    Taln    Tlnq    Tlns    Tdln    sppset
acypi1_2010     14327   388     246     438     581     548     -32     358     227     404     536     506     -53     AMELL,TCAST,daphniaplx
SEM_acypi1_2010 14327   3       2       3       5       4       3       3       2       3       5       4       4       AMELL,TCAST,daphniaplx
acypi2_2011     14580   420     263     467     581     578     -2      395     247     438     546     543     -16     AMELL,TCAST,daphniaplx
SEM_acypi2_2011 14580   4       2       4       5       5       3       4       2       4       5       5       5       AMELL,TCAST,daphniaplx
drosmel_r3_0    14688   413     260     457     579     631     53      391     246     433     548     597     37      AMELL,TCAST,daphniaplx
SEM_drosmel_r3_ 14688   4       2       3       5       5       3       4       2       3       5       5       5       AMELL,TCAST,daphniaplx
drosmel_r5_30   14847   423     267     469     576     655     79      405     255     449     551     627     67      AMELL,TCAST,daphniaplx
SEM_drosmel_r5_ 14847   4       2       4       5       6       3       4       2       4       5       6       6       AMELL,TCAST,daphniaplx
nasvit_ogs12    14692   515     303     478     579     606     27      488     287     452     548     574     14      AMELL,TCAST,daphniaplx
SEM_nasvit_ogs1 14692   4       2       3       5       4       3       4       2       3       5       4       4       AMELL,TCAST,daphniaplx
nvit2_evigenes_ 15165   524     309     488     576     598     22      512     302     477     563     584     24      AMELL,TCAST,daphniaplx
SEM_nvit2_evige 15165   5       3       4       5       5       3       4       2       4       5       5       5       AMELL,TCAST,daphniaplx

=item old output

./sumqual2.sh ref3tabs/ref3-{acypi,drosmel_r[35],n}*.tbest3 << NO, wrong data, need .qual3 w/ prot qual added
group           ns      bit     idn     aln     lnq     lns     dln     aafull  utrok   trln    sppset
acypi1_2010     14327   388     246     438     581     548     -32     0       100     0       AMELL,TCAST,daphniaplx
acypi2_2011     14580   420     263     467     581     578     -2      0       100     0       AMELL,TCAST,daphniaplx
drosmel_r3_0    14688   413     260     457     579     631     53      0       100     0       AMELL,TCAST,daphniaplx
drosmel_r5_30   14847   423     267     469     576     655     79      0       100     0       AMELL,TCAST,daphniaplx
nasvit_ogs12    14692   515     303     478     579     606     27      0       100     0       AMELL,TCAST,daphniaplx
nvit2_evigenes_ 15165   524     309     488     576     598     22      0       100     0       AMELL,TCAST,daphniaplx

=cut

==> ref3sum.info <==
# Add this opt above? need tbest3 pair files to compare

Difference table for paired versions of gene sets : pair t-tests and graph

set pt=acypi; set one=acypi2
cat ref3tabs/ref3-$pt*.tbest3 | grep -v '^#' | sort -k1,1 -k2,2r | env one=$one perl -ne \
'($qd,$sd,@v)=split; ($bs,$id,$al,$ql,$sl)=@v; if($lq eq $qd) { @v2=@v; $s2=$sd;} \
else { puts($lq) ; @v1=@v; $s1=$sd; $s2=0;@v2=(); } $lq=$qd; END{ puts($lq); }\
sub puts{ return unless($lq and $s1); unless(@v2) { @v2=(0)x5; } \
unless($s1=~/$one/) { @t=@v1; @v1=@v2; @v2=@t; ($s1,$s2)=($s2,$s1); } @d=(); \
for $i (0..4) { $d[$i]=int($v1[$i]-$v2[$i]); } print join("\t",$lq,$s1,$s2,"D.$one",@d)."\n"; } \
BEGIN{ $one=$ENV{one}; }' \
  > ref3-$pt.diff3

==> sumqual2.sh <==
#!/bin/bash

## input: ref2tabs/*.qual3 = output of listtall + aaqual score (complete/partial, utrpoor/bad)
## drop @ tcol/avt out for qual

i=0;
for qual in $*; do {
  pt=`basename $qual .qual3| sed 's/^ref3-//; s/\..*//;'`
  i=$(($i + 1))
  cat $qual | env i=$i nam=$pt perl -ne \
'next unless(/^\w/); ($qg,$sg,@v)=split; ($bs,$id,$al,$ql,$sl, $aq, $clen)=@v; ($qsp)= $qg=~m/^(\w+)/; 
($al2,$ap,$ac)=split",",$aq; $cc=($ac=~/complete/)?4:($ac=~/partial(\d)/)?$1:0; ($uq)= $ac=~/(utr\w+)/;
$n++; $ns{$qsp}++; $did{$qg}++; $did{$qsp.$sg}++;
for $i (0..4) { $svs{$qsp}[$i]+=$v[$i];} 
$svs{$qsp}[5] += $sl - $ql; 
$svs{$qsp}[6] += 100 if($cc == 4);
$XXsvs{$qsp}[6] += 100 if($cc == 0);
$svs{$qsp}[7] += 100 unless($uq =~ /utr(bad|poor)/);
$XXsvs{$qsp}[7] += 100 if($uq =~ /utrbad/);
$svs{$qsp}[8] += $clen;
END{ $K=8; @s= sort keys %ns; $ns=0;  @s2=@s; $NT3=$NT*scalar(@s2);
for $s (@s2) { $ns += $ns{$s}; for $i (0..$K) { $sv[$i] += $svs{$s}[$i]; } } 
@av=@avt=(0)x6; 
for $i (0..$K) { $v=$sv[$i]; $av[$i]= int(0.5 + $v/$ns); $atv[$i]=int(0.5 + $v/$NT3); } 
$atv[5]= int(0.5 + $sv[4]/$NT3 - $REFAV); 
@col=qw(bit idn aln lnq lns dln aafull utrok trln); @tcol=map{"T".$_} @col; 
print join("\t",sprintf("%-15s","group"),"ns",@col,"sppset")."\n" if($hd==1);  
$s2=join",",@s2; print join("\t",sprintf("%-15s",$nam),$ns, @av, $s2)."\n"; } 
BEGIN{ $nam=$ENV{nam}; $nam=substr($nam,0,15); $hd=$ENV{i}; $REFAV=560; $NT=5173; }'

} done

## no print: @tcol; @atv
#?? for $s (@s) { $spb{$s} += $svs{$s}[2]; } ($s1,@s2)= sort{$spb{$b}<=>$spb{$a}}keys %spb; \
#..
# ave aa refset
# AMELL   5173    588
# TCAST   5173    557
# daphniaplx      5173    570
# amel+tcas: 572
# amel+dpx: 579
# tcas+dpx: 563

==> listtall2.sh <==
#!/bin/tcsh

## list version: write tall table of top genes picked as in sumtall2.sh
## with not qg or sg limit, need sort input tall3 by max bitscore?
## variant for lowest 2 ref spp; dropping ident-spp
## sort by max align or bits?  -k6,6nr = max align Ooops, 6=refsize, 5=align; redo?

set i=0;
foreach tal ($*)
  set pt=`basename $tal .tall3 | sed 's/^ref3-//; s/\..*//;'`
  set outf=`echo $tal | sed 's/\.tall/.tbest/;'`
  if( -f $outf ) continue; 
  @ i= $i + 1
  cat $tal | grep -v '^Query' | sort -k5,5nr -k3,3nr -k1,1 -k2,2 | env i=1 j=$i nam=$pt perl -ne \
'($qg,$sg,@v)=split; ($bs,$id,$al,$ql,$sl)=@v; ($qsp)= $qg=~m/^(\w+)/; \
unless($did{$qg} or $did{$qsp.$sg}){ $n++; $ns{$qsp}++; $did{$qg}++; $did{$qsp.$sg}++;\
$svs{$qsp}[5] += $sl - $ql; for $i (0..4) { $svs{$qsp}[$i]+=$v[$i];} \
print join("\t",$qg,$sg,@v)."\n"; } \
END{ @s=sort keys %ns; for $s (@s) { $spb{$s} += $svs{$s}[2]; } \
$ns=0; ($s1,@s2)= sort{$spb{$b}<=>$spb{$a}}keys %spb; \
for $s (@s2) { $ns += $ns{$s}; for $i (0..5) { $sv[$i] += $svs{$s}[$i]; } } \
@av=@avt=(0)x6; for $i (0..5) { $v=$sv[$i]; $av[$i]= int(0.5 + $v/$ns); \
$atv[$i]=int(0.5 + $v/$NT3); } $atv[5]= int(0.5 + $sv[4]/$NT3 - $REFAV); \
@col=qw(bit idn aln lnq lns dln); @tcol=map{"T".$_} @col; \
print "# ",join("\t",qw(group ns),@col,@tcol,"sppset")."\n" if($hd==1);  \
$s2=join",",@s2; print "# ",join("\t",$nam,$ns, @av, @atv, $s2)."\n"; } \
BEGIN{ $nam=$ENV{nam}; $hd=$ENV{i}; $REFAV=560; $NT=5173; $NT3=$NT*2; $NT3xxx=$NT*3; }' \
  > $outf

end

# ave aa refset : use instead of REFAV for all?
# AMELL   5173    588
# TCAST   5173    557
# daphniaplx      5173    570
# amel+tcas: 572
# amel+dpx: 579
# tcas+dpx: 563

==> maketall3.sh <==
#!/bin/tcsh

foreach blz ( sd-*blastp.gz )
  set pt=`basename $blz -ref3beebeetdaph.aa.blastp.gz | sed 's/^sd-//;'`
  env aa1=ref3.aa.count aa2=../aaset4/$pt.aa.count aagap=1 swap=0 tall=1 \
  ../makeblastscore2.pl $blz > ref3-$pt.tall3
end

==> sumtall2.sh <==
#!/bin/tcsh

## with not qg or sg limit, need sort input tall3 by max bitscore?
## variant for lowest 2 ref spp; dropping ident-spp
## sort by max align or bits?  -k6,6nr = max align

set i=0;
foreach tal ($*)
  set pt=`basename $tal .tall3 | sed 's/^ref3-//; s/\..*//;'`
  # echo -n "$pt : "
  @ i= $i + 1
  cat $tal | grep -v '^Query' | sort -k6,6nr -k3,3nr -k1,1 -k2,2 | env i=$i nam=$pt perl -ne \
'($qg,$sg,@v)=split; ($bs,$id,$al,$ql,$sl)=@v; ($qsp)= $qg=~m/^(\w+)/; \
unless($did{$qg} or $did{$qsp.$sg}){ $n++; $ns{$qsp}++; $did{$qg}++; $did{$qsp.$sg}++;\
$svs{$qsp}[5] += $sl - $ql; for $i (0..4) { $svs{$qsp}[$i]+=$v[$i];} } $lq=$qg; \
END{ @s=sort keys %ns; for $s (@s) { $spb{$s} += $svs{$s}[2]; } \
$ns=0; ($s1,@s2)= sort{$spb{$b}<=>$spb{$a}}keys %spb; \
for $s (@s2) { $ns += $ns{$s}; for $i (0..5) { $sv[$i] += $svs{$s}[$i]; } } \
@av=@avt=(0)x6; for $i (0..5) { $v=$sv[$i]; $av[$i]= int(0.5 + $v/$ns); \
$atv[$i]=int(0.5 + $v/$NT3); } \
$atv[5]= int(0.5 + $sv[4]/$NT3 - $REFAV); \
@col=qw(bit idn aln lnq lns dln); @tcol=map{"T".$_} @col; \
print join("\t",qw(group ns),@col,@tcol,"sppset")."\n" if($hd==1);  \
$s2=join",",@s2; print join("\t",$nam,$ns, @av, @atv, $s2)."\n"; } \
BEGIN{ $nam=$ENV{nam}; $hd=$ENV{i}; $REFAV=560; $NT=5173; $NT3=$NT*2; $NT3xxx=$NT*3; }'

end

# ave aa refset
# AMELL   5173    588
# TCAST   5173    557
# daphniaplx      5173    570
# amel+tcas: 572
# amel+dpx: 579
# tcas+dpx: 563

