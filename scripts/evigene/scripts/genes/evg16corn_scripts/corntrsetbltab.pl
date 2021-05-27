#!/usr/bin/env perl
# corntrsetbltab.pl

=item usage
# sort by refsize-long, ref-id, ident, targid; change to -k5,5nr align vs ident k4?
# add option for diff3 tables, pairwise 2 gene type sources only

sort -k7,7nr -k2,2 -k4,4nr -k1,1  refSbicolor_313-{evg2corn,maize_[bj]*,tidb*,trsoap*,velv*,trin*}.aa.btall |\
 grep -v utrorf | cat Sbicolor_313.prime.ids - | \
 env merge=1 mina=75 idfilt=1 ref2=2 ../corntrsetbltab.pl > refsorghumprime-evg2best7aa75n.btall1

=cut

use strict;

my $REF2=$ENV{ref2}||0; # refid col 2 or 1?
my $MERGE=$ENV{merge}||0; # merge source types (btall)
my $MINA=$ENV{mina}||25;  ## only for btall now, add to hodiff ? but not used there before.
my $IDFILT=$ENV{idfilt}||0;  ## input id list first
my $USEIDN=$ENV{iden}||0; # use identity vs align score; add bitscore
my $USEBIT=$ENV{bits}||0;
my $EQSLOP=$ENV{eqslop} || 3; # opt?
my $nok=$ENV{nok} || $ENV{nref} || 0; ##for indiff

my $ISCORE=($USEIDN)?3:($USEBIT)?2:4; # idn/aln/bits
my @ISNAM=qw(Queryid Refid Bits Ident Align Quelen Reflen ); # input btall table cols

my($nrd, $eqn,$bbest,$abest,$bmiss,$amiss,$bcov95,$acov95,$sd)= (0) x 19;
my(%ok,%tt,%rv,%rw,%aln,%paln); # global for hodiff3 + others

if($ENV{indiff}) { hodiff3read(); }
elsif($ENV{hodiff}) { hodiff3(); }
elsif($ENV{btall}) { btallone(); } 
else { warn"# opts: env (indiff|hodiff|btall)=1\n"; }

sub btallone {
  my(%did);
  my(%rd,%td,%tdrd,%tspan,%ts,%pid,%pal,%aln,%sumtw); # btallone() tables
  while(<>) { 
    my @v=split; my($td,$rd,$bs,$idn,$aln,$tw,$rw,$tspa,$rspa)=(0) x 9;
    if($REF2) { ($td,$rd,$bs,$idn,$aln,$tw,$rw,$tspa,$rspa)=@v; }
    else { ($rd,$td,$bs,$idn,$aln,$rw,$tw,$rspa,$tspa)=@v; }
  
    if(@v==1){ $ok{$v[0]}++; next; } elsif(/^Query/) { } 
    elsif($IDFILT) { next unless($ok{$rd}); }

    my $ts; ($ts=$td)=~s/[Ll]oc.*//; 
    $ts=~s/(an|sm)velv/velv/;
    my($km)=$ts=~m/k(\d\d+)/; $ts=~s/k$km// if($km); $km||=0;
    if($ts=~s/^tidb//){ $ts=~s/[fr]id/idba/; } 
  
## Zm00001d046083_T001 == ensembl/gramene maize v4 gene set, aug 2016 ens16
    my($sm,$tso); 
if(1) {
    ($sm)= idtype($ts); if($sm eq "miss") { $ts.="miss"; }
} else {
    ($sm)= $ts=~m/(EVm|Zm0|zeamjgi|ncbig|velv|soap|trin|idba|brid)/; 
    unless($sm){ $sm="miss"; $ts.="miss"; } 
}
    $ts=~s/:.*//; # zeamjgi:id, others?
    if($MERGE) { 
       if($MERGE>1 and $sm=~/EVm/) { $ts=~s/EVm\d+t\d+/EVm/; $tso=$ts; } ## Evg dont merge subsets
       else { $tso=$ts; $ts=$sm; }
    }
    elsif($sm=~/EVm/) { $ts=~s/EVm\d+t\d+/EVm/; $tso=$ts; }  ## FIXME, merge overrides, merge 2+ EVm sets
    elsif($sm=~/ncbig/) { $ts=~s/ncbig\d+t\d+/ncbig/; $tso=$ts; } 
    #xx elsif($MERGE) { $ts=$sm; }
    
    s/$/\t$sm:$ts/; 
    my($pid,$pal)= map{ my $p= ($rw>0)?sprintf"%.0f",100*$_/$rw:0; $p=100 if($p>100); $p; } ($idn,$aln); 
    if(/^Query/) { ($pid,$pal)=("pIden","pAlgn"); $sm="Method"; $ts=""; }; 
    unless($did{$rd.$ts}++) { 
      if($pal>=$MINA) { $rd{$rd}++; $td{$ts}{$td}++;  
      if($tspa) { $tdrd{$ts}{$td}{$rd}++; $tspan{$td}{$rd}=$tspa; }
      $ts{$ts}++; $pid{$ts}+=$pid; $pal{$ts}+=$pal;
      $aln{$ts}+=$aln; $sumtw{$ts}+=$tw; } 
      print join("\t",@v,$pid,$pal,"$sm:$ts")."\n"; 
      } 
  }

  my @ts=sort keys %ts;
  my @rd=sort keys %rd; $nrd=@rd; 
  my $nok=($IDFILT)?scalar(keys %ok):$nrd;
  my $pts=sprintf"%.1f",100*$nrd/$nok;
  print "#\n#Summary per method, min $MINA %Align, for $nrd/$nok $pts% found ref transcripts\n"; 
  my @hd=qw(Methd nHits pHits nTids pIden pAlgn Align Tlen Joins);
  if($IDFILT) { @hd=qw(Methd nHits pHits nTids pIdenH pAlgnH AlignH TlenH pAlgnT AlignT Joins); }
  print "# ".join("\t",@hd)."\n";
  for my $ts (@ts) { 
    my @td= sort keys %{$td{$ts}}; my $ntid=scalar(@td);
    my $joins=0;
    
    if(%tspan) { 
      for my $td (@td) { 
        my @trd=sort keys %{$tdrd{$ts}{$td}}; my $ntrd=@trd; 
        for my $i (0..$ntrd-2) { my $ird=$trd[$i];
         for(my $j=$i+1; $j<$ntrd-1; $j++) { my $jrd=$trd[$j];
            my $nov= notover($tspan{$td}{$ird}, $tspan{$td}{$jrd});
            $joins++ if($nov);
            }
          }
      }
    } 
    my $nts=$ts{$ts};
    #o#@ave= map{ sprintf "%.1f",$_/$nrd; }($pid{$ts},$pal{$ts},$aln{$ts},$sumtw{$ts});  
    my @ave= map{ sprintf "%.1f",$_/$nts; }($pid{$ts},$pal{$ts},$aln{$ts},$sumtw{$ts});  
    if($IDFILT) {  #? nok or $nrd here
      my @avet= map{ sprintf "%.1f",$_/$nrd; }($pid{$ts},$pal{$ts},$aln{$ts},$sumtw{$ts});  
      @ave= (@ave,@avet[1,2]);
    }
    my $pts=sprintf"%.1f",100*$nts/$nok; # /$nrd; 
    map{ s/4ms/s/; s/9agv/v/; s/cornhi8m//; s/cornhi12mer3//; s/Zeamay(.EVm)/z$1/; } ($ts);
    print "#".join("\t",$ts,$nts,$pts,$ntid,@ave,$joins)."\n"; 
  } 

}

sub notover {
  my($ispa,$jspa)=@_; #328-600:- ..
  my($ib,$ie,$jb,$je)=map{ my($b,$e)= m/(\d+)\-(\d+)/; ($b,$e) } ($ispa,$jspa);
  my($io,$jo)=(0,0); 
  ($ib,$ie,$io)=($ie,$ib,-1) if($ib>$ie);
  ($jb,$je,$jo)=($je,$jb,-1) if($jb>$je);
  return ($ib < $je and $ie > $jb and $jo eq $io)?0:1
}


=item hodiff3 calc

## hodiff3 calc, better? has pairwise best score

cut -f1 refset/Sbicolor_313main.aa.qual | cat - btallsor/refSbicolor_313-{$pt,$bpt}.aa.btall | \
perl -ne 'if(/^(\S+)$/){$ok{$1}=1} else{($td,$rd)=split; print if($ok{$rd});}' | sort -k2,2 -k5,5nr -k1,1 | \
env idp=0 perl -ne 'BEGIN{$IDP=$ENV{idp}||0; } next if(/^Query|^\W/); @v=split; ($td,$rd,$bs,$idn,$aln,$tw,$rw)=@v;  \
$tt=($IDP)?substr($td,0,6):($td=~/(EVm)\d/)?$1:($td=~m/16tra/)?"trasm16aedes":($td=~/(trin)(c|Loc)/)?$1:substr($td,0,4); \
$rv{$rd}{$tt}=[@v] unless($rv{$rd}{$tt}); $tt{$tt}++; $rw{$rd}=$rw if($rw); END{ @t=reverse sort keys %tt; \
$zv=[0,0,0,0,0,0,0]; for $rd (sort keys %rw) { ($tav,$tbv)= map{ $rv{$rd}{$_} || $zv } @t; \
($aid,$av)=@{$tav}[0,4]; ($bid,$bv)=@{$tbv}[0,4]; $dbs= $av - $bv; $sd += $dbs; $nr++; $rw=$rw{$rd}; \
$paln{$t[0]} += 100*$av/$rw; $paln{$t[1]} += 100*$bv/$rw; if(abs($dbs)<3) { $eqn++; } elsif($dbs<0) { $bbest++; } \
else{ $abest++;} $bmiss++ if($bv ==0); $amiss++ if($av ==0); print join("\t",$rd,$rw,"$aid,$av",$dbs,"$bid,$bv")."\n"; } \
$ad=int($sd/$nr); @pv=map{ $p=int(1000*$_/$nr)/10; "$_,$p%" }($eqn,$abest,$bbest,$amiss,$bmiss); \
@pa=map{ int(10*$paln{$_}/$nr)/10;} @t; \
print "#stat a,b=@t, nr=$nr, aved(a-b)=$ad, refalign%(a,b)=@pa, sumd=$sd\n"; \
print "#best (equal,abest,bbest,amiss,bmiss)=@pv\n";}' > ${pt}_$bpt-sorghumprime.hodiff3 

=cut

sub pctof{ my($n,$d)=@_; my $p=($d>0)?100*$n/$d:0; $p=100 if($p>100); $p; }

sub hodiff3end {
  #? @t=reverse sort keys %tt;
  my @tt;
  my($ta,$tb)= @tt= sort{ $tt{$b}<=>$tt{$a} or $a cmp $b} keys %tt; 
  ($ta,$tb)= sort ($ta,$tb); ## err if @t>2
  if($tt{miss} or @tt>2) { warn "#hodiff: handles only 2 methods, using $ta,$tb of @tt\n"; }
  my @rd=(sort keys %rw); 
  my $zv=[0,0,0,0,0,0,0]; my $nr=0;
  for my $rd (@rd) { 
    $nr++; my $rw=$rw{$rd};
    my($tav,$tbv)= map{ $rv{$rd}{$_} || $zv } ($ta,$tb); # @t; 
    my($aid,$av)=@{$tav}[0,$ISCORE]; my($bid,$bv)=@{$tbv}[0,$ISCORE]; 
    my($pav,$pbv)= map{ pctof($_,$rw) } ($av,$bv);
    $aln{$ta} += $av; $aln{$tb} += $bv;
    $paln{$ta} += $pav; $paln{$tb} += $pbv;
    my $dbs= $av - $bv; $sd += $dbs; 
    if(abs($dbs)<=$EQSLOP) { $eqn++; } elsif($dbs<0) { $bbest++; } else{ $abest++;} 
    $bmiss++ if($bv == 0); $amiss++ if($av == 0); 
    $bcov95++ if($pbv >= 95); $acov95++ if($pav >= 95) ;  ## FIXME for ISCORE <> 3, always use palign here?
    print join("\t",$rd,$rw,"$aid,$av",$dbs,"$bid,$bv")."\n"; 
  } 

  $nrd=@rd; # global
  $nok=($IDFILT)?scalar(keys %ok):$nrd; # global
  hodiff3sum($ta,$tb,$nrd,$nok); # ($nrd,$sd,$eqn,$abest,$bbest,$afound,$bfound)
}

sub hodiff3sum {
  my($ta,$tb,$nrd,$nok)=@_;
    ## bad nrd, nok == bad %rw,%ok ??? from input btalls, hodiff3() > hodiff3end() > hodiff3sum()
  unless($nrd) { warn "#ERR: nrd=$nrd fail \%rw\n"; $nrd=1; }
  unless($nok) { warn "#ERR: nok=$nok fail \%ok\n"; $nok=1; }
  
  my $ad= int($sd/$nrd); 
  #bad: $nrd is total# ($afound,$bfound)=map{ $nok - $_ } ($amiss,$bmiss);
  my($afound,$bfound)=map{ $nrd - $_ } ($amiss,$bmiss);
  my @pv= map{ my $p=int(1000*$_/$nrd)/10; "$_,$p%" }($eqn,$abest,$bbest); # ,$afound,$bfound); 
  my @ph= map{ my $p=int(1000*$_/$nok)/10; "$_,$p%" }($afound,$bfound); 
  push @ph, map{ my $p=int(1000*$_/$nrd)/10; "$_,$p%" }($acov95,$bcov95); 
  my @pa=map{ int(10*$paln{$_}/$nrd)/10;} ($ta,$tb); # @t;  # nrd or nok ?
  my @aa=map{ int(10*$aln{$_}/$nrd)/10;} ($ta,$tb); # @t;  # nrd or nok ?
  #? if($tt{miss} or @tt>2) { print "#hodiff: handles only 2 methods, using $ta,$tb of @tt\n"; }
  print "#stat a,b=$ta,$tb, nref=$nok, refhit=$nrd, score=$ISNAM[$ISCORE],$ISCORE \n";
  print "#stat alndiff(a-b)=$ad, align%(a,b)=@pa, align(a,b)=@aa, sumd=$sd\n"; 
  print "#best (equal,abest,bbest)=@pv\n";
  print "#best (afound,bfound,acov95,bcov95)=@ph\n";
}

sub idtype {
  my($tid)=@_;
  $tid=~s/tidb/idba/;
  # Zm00001d046083_T001 == ens16 ensemble/gramene 2016 aug v4 corn genes
  my($sm)= $tid=~m/(EVm|Zm0|zeamjgi|ncbig|velv|soap|trin|idba|brid)/;
  $sm="miss" unless($sm);
  $tt{$sm}++; 
  return ($sm);
}

sub hodiff3read {
  my($tas,$tbs);
  while(<>) { 
    next if(/^\W/); chomp; my @v=split"\t";
    my($rd,$rw,$tda,$dab,$tdb)=@v;
    my($aid,$av,$bid,$bv)= map{ split",",$_ } ($tda,$tdb); 
    my($ta)= idtype($aid); $tas=$ta unless($ta eq "miss");
    my($tb)= idtype($bid); $tbs=$tb unless($tb eq "miss");
    $nrd++; $rw{$rd}=$rw if($rw);
    my($pav,$pbv)= map{ pctof($_,$rw) } ($av,$bv);
    $aln{$ta} += $av; $aln{$tb} += $bv;
    $paln{$ta} += $pav; $paln{$tb} += $pbv;
    my $dbs= $av - $bv; $sd += $dbs; 
    if(abs($dbs)<=$EQSLOP) { $eqn++; } elsif($dbs<0) { $bbest++; } else{ $abest++;}
    $bmiss++ if($bv == 0); $amiss++ if($av == 0);
    $bcov95++ if($pbv >= 95); $acov95++ if($pav >= 95) ;
  }
  
  #NO, need input order: ($ta,$tb)= @tt= sort{ $tt{$b}<=>$tt{$a} or $a cmp $b} keys %tt; 
  my @rd=sort keys %rw; $nrd=@rd; 
  $nok=$ENV{nok} || $ENV{nref} ||$nrd;  #not: ($IDFILT)?scalar(keys %ok):$nrd
  hodiff3sum($tas,$tbs,$nrd,$nok); # ($nrd,$sd,$eqn,$abest,$bbest,$afound,$bfound)
}


sub hodiff3 {
  # require? $REF2=1;  no, swap @v for rv[]
  my @iswap=(1,0,2,3,4,6,5,8,7);
  while(<>) { 
    next if(/^Query|^\W/); my @v=split;
    my($td,$rd,$bs,$idn,$aln,$tw,$rw,$tspa,$rspa);
    if(@v==1){ $ok{$v[0]}++; next; } elsif(/^Query/) { next; }
    if($REF2) { ($td,$rd,$bs,$idn,$aln,$tw,$rw,$tspa,$rspa)=@v; } #rv[] expects this way
    else { ($rd,$td,$bs,$idn,$aln,$rw,$tw,$rspa,$tspa)=@v; @v=@v[@iswap]; }
  
    #if(@v==1){ $ok{$v[0]}++; next; } elsif(/^Query/) { next; } # elsif()
    if($IDFILT) { next unless($ok{$rd}); }
  
    my($sm)= idtype($td);
    # $ts=$td;
    # #bad#($ts=$td)=~s/[Ll]oc.*//;
    # #skip#$ts=~s/(an|sm)velv/velv/;
    # #skip#($km)=$ts=~m/k(\d\d+)/; $ts=~s/k$km// if($km); $km||=0;
    # ## if($ts=~s/^tidb//){ unless($ts=~s/[fr]id/idba/) { $ts.="idba"; } }
    # $ts=~s/tidb/idba/;
    # ($sm)= $ts=~m/(EVm|zeamjgi|ncbig|velv|soap|trin|idba|brid)/;
    # unless($sm){ $sm="miss"; $ts.="miss"; }
  
    $tt{$sm}++; 
    $rv{$rd}{$sm}=[@v] unless($rv{$rd}{$sm}); 
    $rw{$rd}=$rw if($rw);
    }

 hodiff3end();
}

__END__

