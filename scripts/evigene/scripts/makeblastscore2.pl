#!/usr/bin/env perl
# makeblastscore2.pl   predictor-ortho.blastp predictor-self.blastp > predictor.blscore
#  cut from annotate_predictions.pl

## FIXME: need to look for protein part scores to combine as one
# melon2.% gzgrep 'EG2ap020688t1' blastp/aphid2pub8d-uparphumbac.blastp.gz | grep E0VM80_PEDHC
# # Query:  EG2ap020688t1 aalen=2955,92%,complete
# EG2ap020688t1   E2BNB5_9HYME    54.16   1852    610     43      1258    2955    1070    2836    0.0     1780 << top score is only 1/2 prot
# EG2ap020688t1   E2BNB5_9HYME    47.74   1175    475     32      1       1126    1       1085    0.0      908 << rest
# ..
# EG2ap020688t1   E0VM80_PEDHC    53.72   1763    615     44      1262    2955    1108    2738    0.0     1692 << best score when summed
# EG2ap020688t1   E0VM80_PEDHC    50.25   1184    446     32      1       1122    1       1103    0.0     1027 <<
# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score

#** fixme, makeblastscore needs correction for small-overlap bitscore parts; allow overlap < 1% span?
# .. problem esp those long, repetitive prots: dyenin hvy chn, etc.
#gzgrep Nasvi2EG013019t1 aaeval/bp1-apis2-nvit2_evg11d_all.aa.blastp.gz  | head
#Nasvi2EG013019t1        apis2:XP_396548.4   66.41   3168    1033 14      27      3175    12      3167    0.0     4380
#... slight overlap in query (not targ): 10aa / (3175,800)
#Nasvi2EG013019t1        apis2:XP_396548.4   64.60   644     224  1       3164    3803    3269    3912    0.0      905
# bits=5285  == same full bitscore as below equiv, nearly identical model
#
#gzgrep XP_003426444.1 aaeval/bp1-apis2-nvit2_ncbiref.aa.blastp.gz | head
#ncbiref2:XP_003426444.1 apis2:XP_396548.4   66.46   3166    1036 14      24      3175    12      3165    0.0     4380
#ncbiref2:XP_003426444.1 apis2:XP_396548.4   64.30   647     227  1       3176    3818    3266    3912    0.0      905


use strict;

my $KeyHOMOLOG="homolog";  # make config OPTION
my $KeyPARALOG="paralog";
my $KeyINSPLIT="insplit";
my $OVSLOP=6; # FIXME: need overlapslop ~ < 0.01 of max part span ; or < 0.02..0.05 of min part span?
my $pctOVSLOP=$ENV{pctover} || 0.02; # need opt?
$pctOVSLOP=$pctOVSLOP/100 if($pctOVSLOP > 0.9);

my $swapqt= $ENV{swap} || 0;
my $geneaa= $ENV{aa1} || $ENV{aa} ||"";
my $tgeneaa= $ENV{aa2} ||"";
my $AAGAP = $ENV{aagap} || 0;

my $ONEGENOME=$ENV{oneref}||0; # 2015.02 ?? best spans ignoring scaffold splits .. require showSPAN ?

my $NEEDLEN= $ENV{skipnolen}||$ENV{needlen}||0; # skipnolen require id in aasize/blen or skip

## CDD fixme: need change blastp gnl|CDD|238737  to CDD:nnn

my $skipho= $ENV{skipho} || ""; # pattern to skip homol match, e.g. aphid x swiss _ACYPI = same spp
my $keepho= $ENV{keepho} || ""; # pattern to skip homol match, e.g. aphid x swiss _ACYPI = same spp
my $onlyquery= $ENV{onlyquery} || ""; # species query to keep, or all
$skipho =~ s/[,\s]+/\|/g;
$keepho =~ s/[,\s]+/\|/g;
$onlyquery =~ s/[,\s]+/\|/g;

## pick 1 best target group: TBEST, given ID prefix patt TPRE eg:  '^....'
my $TPRE= $ENV{tpre}||"";
my $TBEST= $ENV{tbest}||""; # $TPRE="" unless($TBEST);
my $pMINLOW= $ENV{pmin} || 0.3; # was 0.3?
my $TALL= $ENV{tall}||0; # all targets >= pmin;
my $showSPAN= $ENV{spans}||0; # 201402, TALL only options output targ,src align spans; make default?
my $pMINBEST2 = 0.75; # dont use TBEST if below 0.5*topscore
my $ADDALIGN= $ENV{align} || $TALL;

my( $bother, $bself)= @ARGV;
( $bother, $bself) = ( $bself, $bother) if($bother =~ /self/);

## fixme for alttr : idtag = t[2-n] : drop from paralog=
my $altids= $ENV{altids} ||""; # table of altid  mainid
my $ALTKEY= $ENV{alt} || 't'; ##'t\d+$';
$ALTKEY .= '\d+$' if($ALTKEY =~ /\w/  and $ALTKEY !~ /\W/);

die "usage: makeblastscore genes-ortho.blastp genes-self.blastp > genes.blscore \n"
  unless( -f $bother or -f $bself);
  
my (%bself, %bparalog, %bother, %tother, %balt, %blen, $cat, $bsort, $cmd, %bspans, %tspan, $lq, $bmax);
      
#o sub bint { local $_=shift; return (/e\+/) ? int($_) : $_; }
sub bint { my $b=shift; return ($b<0)?0 : ($b=~/e\+/)? int($b) : $b; }

my $haveqlen=0;

foreach my $aaset ($geneaa,$tgeneaa) {
if($aaset) {   # add tgeneaa = target set; use same blen hash?
  ## fixme for aagaps:   
  # gunzip -c $az |perl -pe'unless(/^>/){s/X/n/g; s/[A-Z]/a/g;}' |faCount stdin | cut -f1,2,7 |
  my $iscount=0;
  if($aaset =~ /count|qual/) { $iscount=1; open(AASIZE,$aaset) or die "FAIL: faCount  aa=$aaset..."; }
  else { open(AASIZE,"faCount $aaset |") or die "FAIL: faCount  aa=$aaset ..."; }
  my $hasgap=($AAGAP && $iscount)?1:0;
  while(<AASIZE>) { chomp; my($id,$al,@ac)=split"\t"; 
    if($hasgap) { if(@ac>1 or @ac==0){ $hasgap=0; } else { $al -= $ac[0]; }} 
    $blen{$id}=$al; } 
  close(AASIZE); $haveqlen=1;
}
}


my %altids=(); my ($faltid,$naltid)=(0,0);
if($altids) {
  open(F, $altids) or die "ERR: altids=$altids bad file of ids\n";
  while(<F>) { if(/^\w/) { chomp; my($aid,$mid)=split"\t"; $altids{$aid}=$mid; $naltid++; } } close(F);
}

if( $bself and -f $bself ) {
# $cat= ($bself =~ /\.gz/) ? "gunzip -c" : "cat";
# $cmd= join(" ",$cat, $bself, "|");
$cmd= ($bself =~ /\.gz/) ? "gunzip -c $bself |" : $bself;
open(GSCORE, $cmd) or warn "# ERROR: $cmd\n"; # is it file or list? 
while(<GSCORE>) { 
 next unless(/^\w/); chomp; my @v=split"\t"; 
 my($q,$t,$bits,$aln,$mis)= @v[0,1,-1,3,4]; $bits=bint($bits);
 my $qg=$q; my $tg=$t; 
 next if($onlyquery and $q !~ m/$onlyquery/);
 if($altids) { $qg=$altids{$q}||$q; $tg=$altids{$t}||$t; }
 elsif($ALTKEY) { $qg=~s/$ALTKEY//; $tg=~s/$ALTKEY//; }
 
 if($q eq $t) { $bself{$q}= $bits unless($bself{$q}); 
   $blen{$q} += $aln unless($haveqlen); # always == length for self? but for long aa broken blast
 } 
 elsif($qg eq $tg) { $balt{$q}{$t}=1; $faltid++; } # alttr, what? need to keep ids for other match
 else { $bparalog{$q}="$bits,$t" unless($bparalog{$q}); }
 
 } close(GSCORE);
}

my $blerr=0;
if( $bother and -f $bother ) {
# $cat= ($bother =~ /\.gz/) ? "gunzip -c" : "cat";
# $cmd= join(" ",$cat, $bother, "|");
$cmd= ($bother =~ /\.gz/) ? "gunzip -c $bother |" : $bother;
open(GSCORE, $cmd) or warn "# ERROR: $cmd\n"; # is it file or list? 
while(<GSCORE>) { 
  unless(/^\w/) { 
   #?? if(/^# Query: (\S+)/) { my $q=$1; my($al)=m/len=(\d+)/; $blen{$q}=$al if($al); $haveqlen++ if($al); }
   next;} 
  chomp; my @v=split"\t"; 
  unless(@v==12){ warn"ERR: blasttab not 12 cols:'@v'\n"; $blerr++; die if ($blerr>9); next; }
  # my($q,$t,$bits)= @v[0,1,-1];  #old
  my($q,$t,$bits,$aln,$mis,@bspan)=  @v[0,1,11,3,4, 6,7,8,9]; #was @v[0,1,-1,3,4,..]; # 6-9 =  q. start, q. end, s. start, s. end, 

  # bspan = ($qb,$qe,$sb,$se);
  if($swapqt) {
    ($q,$t)= ($t,$q); @bspan= @bspan[2,3,0,1];
  }
  next if($onlyquery and $q !~ m/$onlyquery/);
  $t=~s/^gnl\|CDD\|/CDD:/; # hack fix CDD:

  if($lq and $q ne $lq) {
    ##my($lbits, $lt)= bestscore($lq); 
    my($lbits,$lt,$maxa,$maxi)= bestscore($lq);
    my($tbits, $tt)= split",", $bother{$lq};
    $bother{$lq}="$lbits,$lt,$maxi,$maxa" if($lt and ($lbits > $tbits) 
	    or ($TPRE and $lbits > $pMINBEST2 * $tbits));  # dont change for same score
    %bspans=(); $bmax=0; $lq="";
  }

  next if($skipho and $t =~ m/$skipho/);
  next if($keepho and $t !~ m/$keepho/);
  next if($NEEDLEN and not( $blen{$q} and $blen{$t}));

  $bits= bint($bits);
  if($q ne $t) { 
    my $aident= _max(0,$aln-$mis); 
    # $aident = $pctident * $aln; # other way to calc, maybe better? test
    my $havespan=($ONEGENOME)?1:exists($bspans{$t});
    sumscore( $q, $t, $bits,$aln,$aident, @bspan) if($havespan or $bits > $pMINLOW * $bmax); 
      # ?? limit # lowscore targets here? careful, 1/2 score * 2 can be best
    $bother{$q}="$bits,$t,$aident,$aln" unless($bother{$q}); 
    $bmax= $bits if($bits > $bmax);
    }
  elsif($q eq $t and not $bself) {
    $bself{$q}= $bits unless($bself{$q}); 
    $blen{$q} += $aln unless($haveqlen); # always == length for self? but for long aa broken blast
    if($TALL) {
      my $aident= _max(0,$aln-$mis); # other way to calc: $aident = $pctident * $aln;
      sumscore( $q, $t, $bits,$aln,$aident, @bspan); 
      # $bother{$q}="$bits,$t,$aident,$aln" unless($bother{$q}); 
    }
  } 

  $lq= $q; 
 } close(GSCORE);

  my($lbits,$lt,$maxa,$maxi)= bestscore($lq);
  my($tbits, $tt)= split",", $bother{$lq};
  $bother{$lq}="$lbits,$lt,$maxi,$maxa" if($lt and ($lbits > $tbits) 
    or ($TPRE and $lbits > $pMINBEST2 * $tbits));  # dont change for same score

}      

sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }

## 2011.aug BUG here, need to test sb-se outside tb-te spans also
sub sumscore {
  my( $q, $t, $bits,$aln,$aident, $qb,$qe,$sb,$se) = @_;
  my $or=0;
  my $ttrue=$t; if($ONEGENOME) { $t=$ONEGENOME.$q; } # 2015.02 ?? best spans ignoring scaffold splits ?
  if($qb > $qe) { ($qb,$qe)= ($qe,$qb); $or--; } # FIXME record $or in bspans
  if($sb > $se) { ($sb,$se)= ($se,$sb); $or--; }
  unless($bspans{$t}) { 
    $bspans{$t}=[]; 
    $bspans{$ttrue}=[] if($ONEGENOME); # need for above test bspans{torig} for sumscore()
    push( @{$bspans{$t}}, [$qb,$qe,$sb,$se,$bits,$aln,$aident,$or,$ttrue]); 
    return; }
  my $ov=0;
  ## 2011oct overslop pct fix
  my $qlen=1+$qe-$qb; my $slen=1+$se-$sb;
  my $qslop= _max($OVSLOP, int($pctOVSLOP*$qlen));
  my $sslop= _max($OVSLOP, int($pctOVSLOP*$slen));
  
  ## FIXME? 2013mar: maybe this should take span of overlaps (diff from last), not skip 2nd overlap
  ## .. complex change tho, need to match all query,targ spans and extract non-over pair spans
  ## .. simple eg: q:1-100,50-150  t:200-300,250-350 .. shouldn't happen often, but does for repetitive seq
  
  foreach my $sp (@{$bspans{$t}}) {
    #o my($xb,$xe,$tb,$te,$xbit,$orsp)= @$sp; # BUG here 201407: missing aln,aident!
    my($xb,$xe,$tb,$te)= @$sp;  my $ttrue1= $$sp[8];
    if($qe < $xb or $qb > $xe) { }
    elsif($qe > $xe and $qb >= $xe - $qslop) { }
    elsif($qb < $xb and $qe <= $xb + $qslop) { }
    else { $ov=1; last; }
  ## add 2011.aug
    if($ONEGENOME and $ttrue ne $ttrue1) { } # 2015.02
    elsif($se < $tb or $sb > $te) { }
    elsif($se > $te and $sb >= $te - $sslop) { }
    elsif($sb < $tb and $se <= $tb + $sslop) { }
    else { $ov=1; last; }
  }  
  unless($ov) { push( @{$bspans{$t}}, [$qb,$qe,$sb,$se,$bits,$aln,$aident,$or,$ttrue]); }
}

sub bestscore {
  my($q)= @_;
  my($maxb,$maxt,$maxa,$maxi)= (0) x 9;
  ## 2011.11: best tprefix (pick best of prefixa, but prefixb if no prefixa)
  my(%maxb,%maxt);
  foreach my $t (sort keys %bspans) {
    my @sp= @{$bspans{$t}}; ##ok#delete $bspans{$t}; # 2015.02 fix?
    my ($tbit,$ta,$ti,$mxb,$mxe,$mtb,$mte,$mor,$axs,$ats,$ttrue)= (0) x 19; 
    my (@axs,@ats); $axs=$ats="";
    foreach my $sp (@sp) {
      my($xb,$xe,$tb,$te,$xbit,$aln,$aident,$or,$ttrue1)= @$sp;
      $tbit += $xbit; $ta+= $aln; $ti+= $aident;
      # $axs.="$xb-$xe,"; $ats.="$tb-$te,"; # $showSPAN == 2 all of hsps for detail lovers..
      $ttrue=$ttrue1 unless($ttrue);
      if($ONEGENOME) { 
        my $gr="$ttrue1:";
        ## unless($ttrue) { $ttrue=$ttrue1;} elsif($ttrue1 ne $ttrue) {  }
        push @axs,"$xb-$xe"; push @ats,"$gr$tb-$te"; # $showSPAN == 2 
      } else {
        push @axs,"$xb-$xe"; push @ats,"$tb-$te"; # $showSPAN == 2 
      }
      $mxe=$xe if($xe>$mxe); $mte=$te if($te>$mte);
      $mxb=$xb if($xb<$mxb or $mxb==0); $mtb=$tb if($tb<$mtb or $mtb==0);
      $mor+=$or; #? sum of 0/-1 ?
    }
   if($TALL) { 
    # if($ONEGENOME) { } # problem $t not ttrue val here, change? multi-t gets what? 1st or most common?
    $tother{$q}{$t}="$tbit,$ti,$ta,$ttrue" if(@sp);  ##  $bother{$q}="$bits,$t,$aident,$aln"
    }
    
   ## 20140317: fix spans output format for all hsps, need for some uses..
   ## table cols Qspan Sspan:  1-222,229-333 [tab] 1-222,329-433:+
   ## maybe want sorted hsps? join",", sort { $a <=> $b } split",",$axs; cur order is highest score
   ## but need to sort axs,ats same way as these are paired.
   if($showSPAN) { 
     my @ox= sort{ $axs[$a]<=>$axs[$b] }(0..$#axs);
     my $axs= join",",@axs[@ox]; 
     my $ats= join",",@ats[@ox]; # map{s/,$//} ($axs,$ats); 
     my $oc=($mor<0)?'-':'+';
     
     ##? does $ONEGENOME require listing all ttrue on span summary $mtb-$mte  ?
     ## dang2, span summary $mtb-$mte is not meaningful w/ multiple refs
     if($ONEGENOME) { my $lr=""; my $allr=""; 
      #o# $ats= join",",map{ my($r,$v)= split":",$_,2; $_=$v if($lr eq $r); $lr=$r; $_; } @ats[@ox]; 
      my @atsx=(); for my $ts (@ats[@ox]) { 
        my($r,$v)= split":",$ts,2; 
        if($lr eq $r) { push @atsx,$v; } else { push @atsx,$ts; $allr.="$r,"; }
        $lr=$r; 
        } 
      $ats= join",",@atsx; # if($showSPAN>1);
      $allr=~s/,$//; if($showSPAN<2 and $allr=~/,/) { $mtb="$allr:$mtb"; }
      }
     #o# if($showSPAN==2) { $tspan{$q}{$t}="$axs\t$ats:$oc"; } #? add $mxb-$mxe\t$mtb-$mte also?
     #o2# if($showSPAN==2) { $tspan{$q}{$t}="$mxb-$mxe/$axs\t$mtb-$mte/$ats:$oc"; } #? add $mxb-$mxe\t$mtb-$mte also?
     if($showSPAN==2) { $tspan{$q}{$t}="$mxb-$mxe/$axs\t$ats:$oc"; } #? add $mxb-$mxe\t$mtb-$mte also?
     else { $tspan{$q}{$t}="$mxb-$mxe\t$mtb-$mte:$oc"; } # need flip orient ? $oc has info
     #old# else { $tspan{$q}{$t}=($mor<0)?"$mxb-$mxe\t$mte-$mtb:$oc" : "$mxb-$mxe\t$mtb-$mte:$oc"; } # chg , to \t 
   }
   if($tbit > $maxb) { $maxb=$tbit; $maxt= $t; $maxa=$ta; $maxi=$ti;}
   if($TPRE) { my ($tp)= $t =~ m/^($TPRE)/;  $tp||="other";
     if($tbit>$maxb{$tp}) { $maxb{$tp}=$tbit; $maxt{$tp}=$t; } # add ta,ti
     }
  }

 if($TPRE) {
    my ($tp)= $maxt =~ m/^($TPRE)/;
    unless($tp =~ m/$TBEST/) { 
      my $maxb1=$maxb{$TBEST}; my $maxt1= $maxt{$TBEST}; 
      if($maxt1) { $maxt=$maxt1; $maxb=$maxb1; }# add ta,ti
    }
 }
 return($maxb, $maxt, $maxa, $maxi);
}

# OUTPUT.....................................

my ($ng,$nho)=(0,0);
if($TALL) { # just ho table for all q x t
  my @tcols= qw(Query Source Bits Ident Align Qlen Slen);
  push @tcols, qw(Qspan Sspan) if($showSPAN);
  print join("\t",@tcols)."\n";
  foreach my $gid (sort keys %tother) {
    foreach my $tg (sort{ $tother{$gid}{$b} <=> $tother{$gid}{$a} } keys %{$tother{$gid}} ){
      my $bia= $tother{$gid}{$tg}; 
      my($bit,$idn,$aln,$tgtrue)=split",",$bia; 
      my $qlen= $blen{$gid}||0; 
      my $tlen= $blen{$tg}||0; 
      my $tgval= ($ONEGENOME)? $tgtrue : $tg; # problem $tg not true val here, change? multi-t gets what? 1st or most common?
      my @vcols=($bit,$idn,$aln,$qlen,$tlen);
      if($showSPAN) { my $span=$tspan{$gid}{$tg}||"0\t0"; push @vcols, $span; } # was split",", span
      print join("\t",$gid,$tgval,@vcols)."\n";
    }
  }
  exit; # done

} else {

#?? report balt faltids ?  $balt{qid}{$tid}
my %ball= map{ $_ => 1 } keys %bother, keys %bself;
foreach my $gid (sort keys %ball) { 
  $ng++;
  my @out=();
  my $s=$bself{$gid};    my($saln)= ($s=~s/,(\d+,\d+)$//)?$1:"";
  my $b=$bother{$gid};   my($baln)= ($b=~s/,(\d+,\d+)$//)?$1:"";
  my $p=$bparalog{$gid}; my($paln)= ($p=~s/,(\d+,\d+)$//)?$1:"";
  if($s or $b or $p) {
    #fixabove: map{ if(/e\+/){ my($b,$t)=split","; $b= int($b); $_="$b,$t"; } } ($s,$b,$p); # 1.033e+04
    ## FIXME: these can be e-notation:  1.23e+4 > convert to digits
    ## 12jan: val=$maxbits,tag,$maxa,$maxi  now; cut maxa,maxi unless requested
    
    # overbestgenes expects format  "ho=score/maxscore,.."
    if($b) { $b =~ s=,=/$s,= if($s); push @out, "$KeyHOMOLOG=$b";  $nho++;}
    else { push @out,"na"; }
    if($ADDALIGN){ $baln=~s=,=/=; push @out, (($baln) ? "align=$baln" : "na"); }
      
    if($p) { $p =~ s=,=/$s,= if($s); push @out, "$KeyPARALOG=$p"; }
    else { push @out, "na"; }
    if($s and ($b or $p)) {  
      $b =~ s,/.*,,; $p =~ s,/.*,,; $s =~ s,/.*,,;
      my($bb,$bl)=($p>$b)? ($p,"pa") : ($b,"ho"); $bb=int(100*$bb/$s); 
      push @out, "pHOBEST=$bb\%$bl"; 
      } else { push @out, "na"; }
  }
    
  print join("\t",$gid,"na",@out),"\n";
}
} # not TALL

__END__

=item deltablast CDD special case

  -- collect all non-overlapped CDD hits per query gene;
  -- assume input is delta-blastp, sorted by best matches, only CDD?

set pt=locust1evgc
gunzip -c deblastzunip/sd-uniref*-$pt*.deblastp.gz 
 ..
pt=nasvit_ogs2  
gunzip -c deblastzref3/sd-ref3all-$pt.aa.deblastp.gz | grep 'gnl\|CDD' | perl -ne\
's/gnl.CDD./CDD:/; next unless(/\tCDD:/); ($td,$rd,@v)=split; putc() if($ltd and $ltd ne $td);  \
($pi,$al,$mi,$xm,$tb,$te,$rb,$re,$eval,$bits)=@v; $ov=0; foreach $be (@tbe) { \
($lb,$le)=split"-",$be; if($le>$tb and $lb<$te) { $ov=1; last; } } \
unless($ov) { $tbe="$tb-$te"; push @tbe,$tbe; push @cd,join("\t",$rd,$bits,$al-$mi,$al,$tbe); } $ltd=$td;  \
sub putc{ unless($did{$ltd}++) { foreach my $cd (@cd) { print "$ltd\t$cd\n";} } @cd=@tbe=();} ' \
  > cdd-$pt.domtab

# stats for domtab  : ave/sum ndom/gene, ave/sum aln/gene, num uniq cdd, 

pt=nasvit_ogs2  
cat cdd-$pt.domtab | env na=$pt perl -ne\
'($td,$rd,$bs,$idn,$al,$span)=split; $rds{$rd}++; $tda{$td}+=$al; $tdn{$td}++; $n++; $tdi{$td}+=$idn; \
END{ @td=sort keys %tdn; @rd=sort keys %rds; $nrd=@rd; $ntd=@td; foreach $td (@td) { \
$sa+=$tda{$td}; $si+=$tdi{$td}; $sn+=$tdn{$td}; } ($ava,$avi,$avn)= map{ int(10*$_/$ntd)/10; } ($sa,$si,$sn); \
printf "%-18s", $ENV{na}; \
print "n=$n; ntd=$ntd; nrf=$nrd; ave aln:$ava, idn:$avi, nrf.td=$avn; sum aln:$sa, idn:$si, nrf:$sn\n"; }'

nasvit_ogs12: n=27040; ntd=14251; nrf=5511; ave aln:285.3, idn:121.2, nrf.td=1.8; sum aln:4065946, idn:1727426, nrf:27040
nasvit_ogs2 : n=48279; ntd=26461; nrf=6209; ave aln:269.9, idn:113.3, nrf.td=1.8; sum aln:7142868, idn:2999455, nrf:48279
  ogs2: nrf.td=1.82; ogs1: nrf.td=1.89 < want 2nd digit?
  also need paired genes for test if v2 > v1 for same gene..

acypi2_2011       n=39027; ntd=23236; nrf=6643; ave aln:231.9, idn:98.2,  nrf.td=1.6; sum aln:5390711, idn:2283741, nrf:39027
acypi1_2010       n=29456; ntd=18731; nrf=5992; ave aln:215.5, idn:91.2,  nrf.td=1.5; sum aln:4037478, idn:1708535, nrf:29456

daphplx_evg10     n=45709; ntd=27266; nrf=9422; ave aln:246.2, idn:107.3, nrf.td=1.6; sum aln:6715592, idn:2926802, nrf:45709
daphplx_jgi06     n=28330; ntd=18250; nrf=6175; ave aln:223.5, idn:95.9,  nrf.td=1.5; sum aln:4078877, idn:1751823, nrf:28330

drosmel_r5        n=44633; ntd=19224; nrf=5470; ave aln:358.7, idn:151.4, nrf.td=2.3; sum aln:6895660, idn:2910774, nrf:44633
drosmel_r3        n=30047; ntd=14829; nrf=5423; ave aln:322.1, idn:137.8, nrf.td=2;   sum aln:4777722, idn:2043610, nrf:30047

nasvit_ogs2       n=48279; ntd=26461; nrf=6209; ave aln:269.9, idn:113.3, nrf.td=1.8; sum aln:7142868, idn:2999455, nrf:48279
nasvit_ogs1       n=27040; ntd=14251; nrf=5511; ave aln:285.3, idn:121.2, nrf.td=1.8; sum aln:4065946, idn:1727426, nrf:27040
  
=cut

=item bug align > aasize?

BUG? align > aasize, this is what blastp sez is right ... see below
- should I truncate align to max(aaq,aas) or keep as is? ident lacks this problem. use instead?

number of cases, about 10%:
ref3tabs/ref3-acypi1_2010.tbest3  NT=14327   naln>maxaa= 1363
ref3tabs/ref3-acypi2_2011.tbest3  NT=14580   naln>maxaa= 1564
ref3tabs/ref3-drosmel_r5_30.tbest3 NT=14847  naln>maxaa= 1341
ref3tabs/ref3-nvit2_evigene.tbest3 NT=15165  naln>maxaa= 1892
ref3tabs/ref3-nasvit_ogs12.tbest3  NT=14692  naln>maxaa= 1678


cat ref3tabs/ref3-acypi*.tbest3 | grep -v '^#' | sort -k1,1 -k2,2r | head -40

Query           Source                  Bits    Ident   Align   Qlen    Slen
AMELL:GB40063   acypi2:ACYPI006093      421     293     635     625     633
AMELL:GB40063   acypi1:ACYPI006093-PA   374     308     645     625     617  << align too big, has 337 mismat 

aasize: AMELL:GB40063   625     0
aasize: acypi1:ACYPI006093-PA   618     1

# Query: AMELL:GB40063 GB40063-PA H9KSZ7 Uncharacterized protein IPR001208,IPR003593
# qid           sid                     piden   alen    mism    gap     qb      qe      sb      se      eval     bits  
AMELL:GB40063   acypi1:ACYPI006093-PA   36.59   645     337     11      1       608     9       618     1e-103   374
..
AMELL:GB40063   acypi1:ACYPI003255-PA   32.20   531     322     16      76      583     263     778     8e-69    258
AMELL:GB40063   acypi1:ACYPI005644-PA   33.21   557     293     17      87      583     110     647     2e-61    233
AMELL:GB40063   acypi1:ACYPI004569-PA   29.30   604     362     18      2       576     226     793     4e-58    222
AMELL:GB40063   acypi1:ACYPI008927-PA   31.10   598     330     20      19      565     53      619     1e-57    221
AMELL:GB40063   acypi1:ACYPI007914-PA   39.94   318     161     7       261     565     334     634     1e-55    214
AMELL:GB40063   acypi1:ACYPI003302-PA   39.94   318     161     7       261     565     255     555     2e-55    213
AMELL:GB40063   acypi1:ACYPI002361-PA   30.75   504     322     9       96      589     162     648     9e-50    194
AMELL:GB40063   acypi1:ACYPI001075-PA   31.40   449     252     14      150     583     139     546     9e-46    181

=cut

