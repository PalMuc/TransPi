#!/usr/bin/env perl
# estblastcover.pl 
#  transcripts coverage by ests, for blastn -query trs -db ests
#  cut from makeblastscore // annotate_predictions.pl

use strict;

my $OVSLOP=6; # FIXME: need overlapslop ~ < 0.01 of max part span ; or < 0.02..0.05 of min part span?
my $pctOVSLOP=$ENV{pctover} || 0.02; # need opt?
$pctOVSLOP=$pctOVSLOP/100 if($pctOVSLOP > 0.9);

# my $TALL= $ENV{tall}|| 1; # all targets >= pmin;  #? only this format?
# my $skipho= $ENV{skipho} || ""; # pattern to skip homol match, e.g. aphid x swiss _ACYPI = same spp

## pick 1 best target group: TBEST, given ID prefix patt TPRE eg:  '^....'
# my $TPRE= $ENV{tpre}||"";
# my $TBEST= $ENV{tbest}||""; # $TPRE="" unless($TBEST);
# my $pMINLOW= $ENV{pmin} || 0.3; # was 0.3?
# my $pMINBEST2 = 0.75; # dont use TBEST if below 0.5*topscore
# my $ADDALIGN= $ENV{align} || $TALL;

my $genesize= $ENV{'size'} ||"";
my $geneaa = $ENV{'tr'} || $ENV{aa} ||"";
my $CDSONLY= $ENV{'cds'} || 0;
my $SWAPQT= $ENV{'swap'} || 0;

# my( $bother, $bself)= @ARGV;
my( $bother)= @ARGV;
# ( $bother, $bself) = ( $bself, $bother) if($bother =~ /self/);

## fixme for alttr : idtag = t[2-n] : drop from paralog=
# my $ALTKEY= $ENV{alt} || 't'; ##'t\d+$';
# $ALTKEY .= '\d+$' if($ALTKEY =~ /\w/  and $ALTKEY !~ /\W/);

die "usage: env size=genesize.tab cds=1|0 estblastcover gene2est.blastn > genes.score \n"
  unless( -f $bother or ($bother =~ /^stdin|^-/i));
  
# my (%bself, %bparalog, %bother, %tother, %balt, %blen, $cat, $bsort, $cmd, @bspans, $lq, $bmax);
my (  %blen, %bcdsoff, %bcdslen, @bspans, %bmated, $lq, $bmax);
      
sub bint { local $_=shift; return (/e\+/) ? int($_) : $_; }

## replace geneaa/tr input w/ gene.sizetab 
# Query   Qlen    CDSlen  CDSoffs
# daphmag3tri7trimsub13loc1004c0t1        2213    1800    201-2003
# daphmag3tri7trimsub13loc1010c0t2        2216    1659    336-1997
# daphmag3tri7trimsub13loc1012c0t2        2642    1650    374-2026

my $haveqlen=0;
if($genesize) {
  open(AASIZE,$genesize) or die "FAIL: size=$genesize ...";
  while(<AASIZE>) { 
  my($id,$trsize,$cdssize,$cdsoff)=split; 
  my($b,$e)=split /[.-]+/,$cdsoff; 
  $blen{$id}=$trsize; $bcdslen{$id}=$cdssize; $bcdsoff{$id}=[$b,$e]; 
  } 
  close(AASIZE); $haveqlen=2;
} elsif($geneaa) { # warn:  not -f $geneaa
  open(AASIZE,"faCount $geneaa |") or die "FAIL: faCount  aa=$geneaa ...";
  while(<AASIZE>) { my($id,$al)=split; $blen{$id}=$al; } close(AASIZE); $haveqlen=1;
}
$CDSONLY=0 unless($haveqlen == 2);

my $inh;
if($bother =~ /^stdin|^-/i) {
  $inh= *STDIN
} elsif(not( $bother and -f $bother )) {
  die "# ERROR: Missing blastn input: $bother\n";
} elsif($bother =~ /\.gz/) {
  open(GSCORE,"gunzip -c $bother|")  or die "# ERROR: $bother\n"; $inh=*GSCORE;
} else {
  open(GSCORE,$bother)  or die "# ERROR: $bother\n"; $inh=*GSCORE;
}

## add maxgap ?
#  my $score= join("\t",$tbit,$ti,$ta,$tcov,"$gb-$ge",$tgap,$gaps); 
my @SCORES= (); # do below bestscore
# @SCORES= qw(Bits Iden Algn Ncov Covs Ngap Maxgap Gaps);
# print join("\t",qw(Query Qlen),@SCORES)."\n";

while(<$inh>) { 
  unless(/^\w/) { 
   #?? if(/^# Query: (\S+)/) { my $q=$1; my($al)=m/len=(\d+)/; $blen{$q}=$al if($al); $haveqlen++ if($al); }
   next;} 
  my @v=split; 
  my($q,$t,$bits,$aln,$mis,@bspan)= @v[0,1,-1,3,4, 6,7,8,9]; # 6-9 =  q. start, q. end, s. start, s. end, 
  # q == gene transcript, t == est 
  # add opt to swap q,t : need sort on t
  if($SWAPQT) { ($q,$t)= ($t,$q);  @bspan= @bspan[2,3, 0,1]; }
  
  ## %bother{query} not used now; see %tother{query}
  if($lq and $q ne $lq) {
      bestscore($lq); # now does putout/gene
#     my($lbits,$lt,$maxa,$maxi)= bestscore($lq);
#     my($tbits, $tt)= split",", $bother{$lq};
#     $bother{$lq}="$lbits,$lt,$maxi,$maxa" if($lt and ($lbits > $tbits) 
# 	    ); ##? or ($TPRE and $lbits > $pMINBEST2 * $tbits));  # dont change for same score
    @bspans=(); %bmated=(); $bmax=0;
  }

#  next if($skipho and $t =~ m/$skipho/); # not?

  $bits= bint($bits);
  if($q ne $t) { 
    my $aident= _max(0,$aln-$mis); # other way to calc: $aident = $pctident * $aln;
    sumscore( $q, $t, $bits,$aln,$aident, @bspan); ## if($bits > $pMINLOW * $bmax or @bspans); 
      # ?? limit # lowscore targets here? careful, 1/2 score * 2 can be best
    # $bother{$q}="$bits,$t,$aident,$aln" unless($bother{$q}); 
    $bmax= $bits if($bits > $bmax);
  }

  $lq= $q; 
 } close($inh);

bestscore($lq); # now does putout/gene

# my($lbits,$lt,$maxa,$maxi)= bestscore($lq);
# my($tbits, $tt)= split",", $bother{$lq};
# $bother{$lq}="$lbits,$lt,$maxi,$maxa" if($lt and ($lbits > $tbits) 
#   ); #?? or ($TPRE and $lbits > $pMINBEST2 * $tbits));  # dont change for same score


#.....................

sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }

# sub putout
# {
#   my($gid,$scores)= @_;
#   my $qlen= $blen{$gid}||0; # get elsewhere? OR include gid x gid scores in tother for self score?
#   # unless($qlen) { my @s= split"\t",$scores; $qlen=$s[3]; }
#   print join("\t",$gid,$qlen,$scores)."\n";
# }

sub matespan {
  my($pairid,$eid,$b,$e)=@_;
  $bmated{$pairid}{$eid}++;
  my($bmin,$emax)= ($b,$e);
  my ($mateid)= grep { $_ ne $eid } keys %{$bmated{$pairid}};
  return ($bmin,$emax) unless($mateid);
  #** Fixme, this reduces ident/align ratio by amount of untested insert span added here.. 
  foreach my $sp (@bspans) {
    my($xb,$xe,$tb,$te,$sbit,$sal,$sid,$starg)= @$sp;
    if($starg eq $mateid) { 
      $bmin= ($xb<$b)? $xb : $b;
      $emax= ($xe>$e)? $xe : $e;
      $sp->[0]= $bmin; $sp->[1]= $emax;
    }
  }
  return ($bmin,$emax);
}

# change sumscore, tESTid not needed, want union of all gene-bspans 

sub sumscore {
  my( $q, $t, $bits,$aln,$aident, $qb,$qe,$sb,$se) = @_;
  my $or=0;
  if($qb > $qe) { ($qb,$qe)= ($qe,$qb); $or--; }
  if($sb > $se) { ($sb,$se)= ($se,$sb); $or--; }

  if($CDSONLY) {
    my $cdsoff= $bcdsoff{$q} or return ;  
    my ($cdsb,$cdse)= @$cdsoff;
    return if($qe < $cdsb or $qb > $cdse);
    $qb=$cdsb if($qb < $cdsb);
    $qe=$cdse if($qe > $cdse);
  }
  
  # add mate-pair score, esp. add span between mates as covered..
  my $tmate=$t; 
  if( $tmate =~ s/\.(fwd|rev)$//i or $tmate =~ s,/([12])$,,) {
    ($qb,$qe)= matespan($tmate, $t, $qb,$qe); # $bmated{$tmate}{$t}++; 
    }
  
  unless(@bspans) { 
    push( @bspans, [$qb,$qe,$sb,$se,$bits,$aln,$aident,$t]); 
    return; }
    
  my $ov=0;
  
use constant OVEROK => 0; # test it, no this gathers too much junk, low qual aligns.  
# need more careful cutting of EST parts that fill gaps in align.
use constant OVERCUT => 1; # test it   

unless(OVEROK) {
  my $qlen=1+$qe-$qb; my $slen=1+$se-$sb;
  my $qslop= _max($OVSLOP, int($pctOVSLOP*$qlen));
  my $sslop= _max($OVSLOP, int($pctOVSLOP*$slen));
  my $cuts=0;
  # add case for this span exceeds last span?
  # ** should change this overlap filter so that part-overs are retained, part outside of last span.
  foreach my $sp (@bspans) {
    my($xb,$xe,$tb,$te,$sbit,$sal,$sid,$starg)= @$sp;
    unless(
      ($qe < $xb or $qb > $xe) or
      ($qe > $xe and $qb >= $xe - $qslop) or
      ($qb < $xb and $qe <= $xb + $qslop)) { 
      
      # change here, if qe > xe + qslop, or qb < xb-slop, then cut qb,qe to part outside of xb,xe..
      # need to cut sbit, sal, sid also .. messy
if(OVERCUT) {
      my ($cb,$ce)=(0,0);
      if($qe > $xe + $qslop) { ($cb,$ce)= ($xe, $qe); }
      elsif($qb < $xb - $qslop) { ($cb,$ce)= ($qb, $xb); }
      if($ce>0 and $cuts==0) { my $pc= ($ce-$cb)/($qe-$qb);
        $bits=int($pc*$bits); $aln=int($pc*$aln); $aident=int($pc*$aident); 
        ($qb,$qe)= ($cb,$ce); $cuts++;
      } else { $ov=1; last; }
} else {
      $ov=1; last; 
}      
      }

    if( $starg eq $t) {
    unless(  # target spans, skip this for est unless same target-id
      ($se < $tb or $sb > $te) or
      ($se > $te and $sb >= $te - $sslop) or
      ($sb < $tb and $se <= $tb + $sslop) ) { $ov=1; last; }
    }
  }
}    

  unless($ov) { push( @bspans, [$qb,$qe,$sb,$se,$bits,$aln,$aident,$t]); }
}


sub bestscore {
  my($gid)= @_;
  # my( $maxb,$maxt,$maxa,$maxi)= (0) x 5;
  my ($tbit,$taln,$tident,$t,$nid)= (0) x 9;
  my $eids="";
  
  @bspans= sort{ $a->[0] <=> $b->[0] } @bspans;
  my($gb,$ge)= ($bspans[0]->[0], $bspans[-1]->[1]);
  my $galignspan=1+$ge-$gb;  
  my $gsize= $blen{$gid} || $ge;  # ge = best guess...
  
  if($CDSONLY) { $gsize= $bcdslen{$gid} || $ge;  }
  my $cdsoff= $bcdsoff{$gid} || [$gb,$ge];  # ge = best guess...
  my ($cdsb,$cdse)= @$cdsoff;
  
  my @cov; $#cov= $ge; ## = (0) x $ge;
  my @cdscov; $#cdscov= $ge; ## = (0) x $ge;
  # .. change cdscov to include any hit to cds, with extensions past..
  # .. extend it also with utr-covers that also hit existing cdscov, so that
  # .. cov - cdscov becomes possible spurious extensions, report that, as tcov - ccov
  my @eids=();
  foreach my $sp (@bspans) {
    my($xb,$xe,$tb,$te,$xbit,$aln,$aident,$targid)= @$sp; 
    $tbit += $xbit; $taln+= $aln; $tident+= $aident;
    push @eids, $targid if(@eids < 3); ## $eids.="$targid," if($nid++ < 3);
    my $xmid= int( ($xb+$xe)/2);
    for my $i ($xb..$xe) { $cov[$i]++; } 
    # # $cdscov[$i]++ if($i>=$cdsb and $i<=$cdse); 
    if( ($xb >= $cdsb and $xb < $cdse) or ($xe > $cdsb and $xe <= $cdse)
      #? or $cdscov[ $xmid ] # which?
      or ($cdscov[ $xb+1 ] or $cdscov [ $xe-1 ]) # extend if already overlapping this xbe
      ) {
      for my $i ($xb..$xe) { $cdscov[$i]++; } 
      }
    }
  $eids= join(",",sort @eids);
    
  my $gaps=""; 
  my ($tcov,$tgap,$gapb,$gape,$gapmax, $ccov)= (0) x 10;
  for(my $i=$gb; $i<=$ge; $i++) { 
    if($cdscov[$i]) { $ccov++; } 
    if($cov[$i]) { $tcov++; 
      if($gape) { $gaps.="$gapb-$gape,"; my $gw= 1+$gape-$gapb; $tgap += $gw; $gapmax=$gw if($gw>$gapmax); } 
      $gapb=$gape=0; }
    else { if($gapb) { $gape=$i; } else { $gapb=$i; } }
  } 
  if($gape) { $gaps.="$gapb-$gape,"; my $gw= 1+$gape-$gapb; $tgap += $gw; $gapmax=$gw if($gw>$gapmax); } 
  $gaps=0 unless($gaps);
 
  my($tmated,$tpaired,$unpaired)=(0,0,0);
  if(scalar(%bmated)) {
    my @mates= sort keys %bmated; 
    foreach my $m (@mates) { $tmated++; my @two= sort keys %{$bmated{$m}}; 
      if(@two>1) { $tpaired++; } else { $unpaired++; }
    }
  }
  my $mates= ($tmated==0)? 0 : "$tpaired/$tmated";

  # dont need both taln, tcov .. should be same, near same but for slop above
  # FIXME added mate insert span, tcov larger than taln for mated spans; revert to put both taln,tcov
  #NO# $taln= $tcov; # best one for now .. 
  
  # FIXME: add top 1,2 EST ids for matching transcripts by best EST id. and/or use best per span
  # FOXME: change tr=tr.fa input to table of lengths w/ CDS len and span, UTR len and span
  # ...    modify above sums to provide CDS-only and CDS-UTR cover to hunt for stat of gene joins.
  
  ## Scores from trasm review paper:
  # Accuracy = Ident / Total-align (should include gaps)  
  # my $accuracy = int(0.5+100 * $tident / $galignspan);
  # Completeness = alignment/total-size > criterion (0.80) 
  # my $complete = int(0.5+100 * $taln/$gsize); ## > 80% / Total-align (should include gaps)  
  
  # Contiguity = ESTs fully (80%?) covered by transcript; not here..
  # Chimeric/joingene : use maxgap > xxx ? 
  
   # add $tcov/$blen == pct-align,  $tident/$taln = pct-ident ?
  my $cdsUncov= ($ccov==0)? 0 : _max(0, $tcov - $ccov); # == UTRgap? only for not CDSONLY
  # ^^ got bad scores for some, is ccov bad?  val == gsize or nearly so when tcov is large.
  
  my $scores= join("\t",$gsize,$tbit,$tident,$taln,$tcov,"$gb-$ge",$mates,$cdsUncov,$tgap,$gapmax,$gaps, $eids); 
  unless(@SCORES) { # header
  if($CDSONLY) { 
  @SCORES= qw(CDSlen Bits Iden Algn Ncov Covs Mated UTRgap Ngap Maxgap Gaps ESTids);
  } else {
  @SCORES= qw(Trlen Bits Iden Algn Ncov Covs Mated UTRgap Ngap Maxgap Gaps ESTids);
  }
  print join("\t","Query",@SCORES)."\n";
  }
  print join("\t",$gid,$scores)."\n";
  return($scores);
}
  
#  $tother{$gid}= $score; # dont need hash now
#  #.. drop this, tother now only output
#  if($tbit > $maxb) { $maxb=$tbit; $maxt= $t; $maxa=$taln; $maxi=$tident;}
#  return($maxb, $maxt, $maxa, $maxi);


__END__


TEST case with only arp11-bestarp genes (6k to 8k)
-- not sure the stats below are showing useful results; but this subset is selected for best proteins of each method.
-- still need good stat for joins/aberrant-UTRs .. CDS/trsize ratio?  cover gaps b/n cds and long utrs?

env size=best3meth.sizetab cds=0 $evigene/scripts/estblastcover.pl bl-daphmagplx_estlong-best3meth.tr.blastn.gz >  daphest-best3meth.covtab &
env size=best3meth.sizetab cds=1 $evigene/scripts/estblastcover.pl bl-daphmagplx_estlong-best3meth.tr.blastn.gz > daphest-best3meth.cdscov &


# fixme: aln vs ncov, later includes matepair insert. use ncov for comp=ncov/slen and gap, use aln for ident/aln

cat daphest-best3meth.covtab | egrep '^Query|^daphmag3tri' | perl -ne'($gd,$gl,$b,$id,$aln,$ncov,$cov,$mt,$n
gap,$mgap,$gaps)=split; ($mh,$mt)=split"/",$mt; $smh+=$mh;  $sm+=$mt;  $ng++; $slen+=$gl; $saln+=$aln; $scov+=$ncov; $sid+=$id; $sgap+
=$mgap; END{ $ac=int(100* $sid/$saln); $comp=int(100*$scov/$slen); $gap=int(100*$sgap/$scov); $mate=int(100*$smh/$sm);
print "ng=$ng; acc=$ac; compl=$comp; gap=$gap; mate=$mate; sid=$sid; saln=$saln; scov=$scov; smated=$smh; sgap=$sgap; slen=$slen\n";
 }'

trin: ng=6998; acc=92; compl=61; gap=19; mate=75; sid=8692575; saln=9372766; scov=10295401; smated=36716; sgap=2010479; slen=16700586
velo: ng=7130; acc=92; compl=62; gap=20; mate=74; sid=8609427; saln=9293655; scov=10197253; smated=38959; sgap=2049265; slen=16367370
cuff: ng=5388; acc=94; compl=58; gap=22; mate=74; sid=9944577; saln=10516752; scov=11520153; smated=31908; sgap=2596667; slen=19810007

#// only saln..
# trin: ng=6998; acc=92; compl=56; gap=21; mate=75; sid=8692575; saln=9372766; smated=36716; sgap=2010479; slen=16700586
# velo: ng=7130; acc=92; compl=56; gap=22; mate=74; sid=8609427; saln=9293655; smated=38959; sgap=2049265; slen=16367370
# cuff: ng=5388; acc=94; compl=53; gap=24; mate=74; sid=9944577; saln=10516752; smated=31908; sgap=2596667; slen=19810007


CDSONLY
cat daphest-best3meth.cdscov | egrep '^Query|^daphmag3tri' | perl -ne'($gd,$gl,$b,$id,$aln,$ncov,$cov,$mt,$n
gap,$mgap,$gaps)=split; ($mh,$mt)=split"/",$mt; $smh+=$mh;  $sm+=$mt;  $ng++; $slen+=$gl;  $scov+=$ncov; $saln+=$aln; $s
id+=$id; $sgap+=$mgap; END{ $ac=int(100* $sid/$saln); $comp=int(100*$scov/$slen); $gap=int(100*$sgap/$scov); $mate=int(1
00*$smh/$sm);  print "ng=$ng; acc=$ac; compl=$comp; gap=$gap; mate=$mate; sid=$sid; saln=$saln; scov=$scov; smated=$smh;
 sgap=$sgap; slen=$slen\n"; }'

trin: ng=6998; acc=92; compl=62; gap=21; mate=75; sid=7109282; saln=7710873; scov=7420270; smated=31071; sgap=1610737; slen=11799441
velo: ng=7130; acc=92; compl=62; gap=21; mate=75; sid=7264605; saln=7882379; scov=7636392; smated=32335; sgap=1630383; slen=12146154
cuff: ng=5388; acc=93; compl=61; gap=23; mate=74; sid=5317929; saln=5710743; scov=5493081; smated=20264; sgap=1305201; slen=8891049

