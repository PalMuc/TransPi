#!/usr/bin/env perl
# evigene3compare.pl

=item about

  classify best gene model per locus using homology, expression scores.
  
  inputs: 
    1. evidence tables of gene id, h=homology, x=expression scores, a=alternate transcripts, per gene set
        equal/*.equal.ids
    2. equal.id pairwise tables of 2 gene sets with equivalent transcript ids
         transcript-*-hxscore.txt == rnas/tr*.hoex.tab
  
  steps:
    1.  read input tables
    2.  merge as compare3 table, row of equivalent models per locus with evidence scores
    3.  classify best model per 1 locus row, using scores
    4.  classify join/split loci over several rows

  cat equal/*.equal.ids | sed 's/^/eq /' | cat - rnas/tr*.hoex.tab | \
  $evigene/scripts/evigene3compare.pl > compare3-hoex7.tab
  
  see also below for old parts

=item author
  
  don gilbert, gilbertd near indiana edu, ca 2010-2011
  part of EvidentialGene, http://arthropods.eugenes.org/EvidentialGene/evigene/scripts/

=cut

use strict;

my $debug=1;

my $DX=9; # express base score same
my $DH=9; # bitscore same

# discrepancy action: proportions to take express diff as best over homology diff
my $XpH=0.76; my $XpHi= 1.24;  # express best for this divergence, pX-ratio < XpH or pX > XpHi
my $HpX=0.89; my $HpXi= 1.10;  # homology same if HpX < pH-ratio < HpXi

my $X_LARGE = 250; # orig: 500; #? too big, maybe 250?
  my $X_SMALL = 100; # min scores to use discrepancy rule
  my $X_SMALL3 = 400; 
my $H_LARGE = 60;  
  my $H_SMALL = 40;
  my $H_SMALL3 = 120;  # h,x small3 filter for lowqual, single-source gene set Acy1 = GNUM_NONE_SKIP

my $SPLIT_MAJOR = 1.25; # classify0 rule 6SJ
my $SPLIT_MINOR = 1.50; # score factors to accept split over join

my $SPLIT_NEAR_JOIN = 0.80; # classifyJoinSplit: decide when to test sum(split-scores)  
my $SPLIT_NEAR_JOIN_MAJOR = 0.66; # or 0.66 ??

my %gname = ( 0 => "error", 1 => "Evig", 2 => "Ncbi", 3 => "Acy1", 4 => "other" );
sub gnum { local $_=shift; return (/acyp2eg/)?1:(/ACYP/)?3:(/^na/)?4:(/^[XNLYG]/)?2:4; } 
sub gname { $gname{ gnum(@_) }; }
sub _gsort{ my $ai=gnum($a); my $bi=gnum($b); return($ai <=> $bi or $a cmp $b); } 
    ## idbest for Noscore: choose one of Ncbi, Evigene if there, ignore ACYPI1 
    ## some of the Noscore's are data mistakes, esp missing expression scores

my @GNUM_NONE_KEEP = (1,2);  # higher qual genesets
my @GNUM_NONE_SKIP = (3); # lower qual geneset; no BestID if Noscore

#----------------------------------------------------

my(%eq, %ec, %ha, %hx, %ta);
my(%did,$didhead);

# read data; could be improved, note required /^eq / tag on id equivalence rows
while(<>) {
  chomp;  
  next unless(/^\w/);
  if(s/^eq //){  # equal.ids tables
    my($d,$e,$pi,$lc)=split"\t"; 
    $ec{$d.$e}="$pi,$lc"; 
    $eq{$d}.="$e," unless($d eq "na"); 
    $eq{$e}.="$d," unless($e eq "na");  

   } elsif(/\t(h|x|a)=/){ # transcript-hxscore.tables
    s/\t(h|x|a)=/\t/g; 
    my($d,$h,$x,$a)=split"\t"; 
    my $gg=gnum($d); 
    $ha{$d}= $hx{$gg}{$d}="$h,$x"; 
    
    $ta{$d}= $a if($a); #?? no? late update tab20
    map{ s,/.*,,; $ta{$_}=$d if(/\w/); } split",",$a; #? save alt scores for tests? or use prime tr score?
    
   } 
}

sub equalhigh {
  my($id,$ed,$pct)= @_;
  $pct ||= 66; 
  my $ecx= $ec{$id.$ed} || $ec{$ed.$id} || 0;  # == "$pi,$lc"; pi == [CI]pc.px
  my ($ec,$ex)= $ecx =~ m/(\d+)\.(\d+)/;
  # $ecx =~ s/,.*//; my ($ec,$ex)= split( /[\.]/, $ecx);
  return (($ec >= $pct or $ex >= $pct) ? 1 : 0); # ($ec =~ /^[CI]/ or 
}

# sub classifyAll

foreach my $gg (sort keys %hx) { 
  foreach my $id (sort keys %{$hx{$gg}}) 
  { 
    next if($did{$id}); # FIXME below; defer $did{id}++ to after have idbest
    
    my @hx=(0) x 4; 
    my $manyto1=0; my %many= ();
    my %groups=();
    
    # problem here, when many eq score > primary eq d,e scores
    # .. solve in classifyJoinSplit
    my $i= gnum($id); # 1..3, 4=other
    $hx[$i]= "$id,".$hx{$gg}{$id}; 
    $groups{$i}++;
    
    my $eq= $eq{$id}; 
    my @eq= ($eq) ? (sort _gsort split",",$eq ) : (); 
    my %deq= (); 
    $deq{$id}= 1; $deq{0}= 1; # map{ $_,1 } ($id,@eq);
    my %eqOfId= map{ $_,1 } @eq;
    
    # id is best of its alts, but eq may not be best of its alts
    foreach my $e (@eq) { 
      my $ha=$ha{$e}; 
      
      ## problem w/ alts is that equal.ids have some very partial equivalences to alts that mess up use
      ## need to know if alt of eq is good or bad equiv of this id row; skip if bad, test/drop if good equiv
      
#      # not right adding alts; losing non-alt near gene overlaps lower scored alt
## this revision better, dropping mostly dupl true alts
## >> only these 2 mistakes? and removed maybe all duplicated alt sets?
## >> AND these 2 evig are added back at end as no-eq models
## but for cases of ncbi "alts" that are not alt-tr but different genes, overlapping utrs, eg
##   XM_001950706 and XM_003245462.1 marked as alt is 2nd gene == acyp2eg0021327t1
# < 1.Evig        acyp2eg0021327t1        acyp2eg0021327t1,246,1197       0       ACYPI001249-RA,246,1004,C99.79,434
# sc:88963-100867 0HX,5Z,5L,
# < 1.Same2       acyp2eg0021328t1        acyp2eg0021328t1,330,1720       XM_001950706.2,330,1720,C100.94,434sc:8896
# 3-90939 ACYPI005939-RA,330,1648,C100.95,434sc:88963-90939       0HX,3S,5L,
# ---
# > 2.Ncbi        XM_001950706.2  acyp2eg0021327t1,246,1197       XM_001950706.2,330,1720,nil     ACYPI001249-RA,246
# ,1004,C99.79,434sc:88963-100867 0HX,1H,5L,
# 30818,30819c30816
# .. same problem, 2 ncbi alts are diff genes
# < 1.Evig        acyp2eg0029519t1        acyp2eg0029519t1,279,1424       0       ACYPI005250-RA,277,1303,C94.87,124
# 3sc:127396-140035       0HX,5Z,5L,
# < 1.Same2       acyp2eg0029520t1        acyp2eg0029520t1,451,1562       XM_001949152.2,451,1562,I100,1243sc:136773
# -140035 ACYPI007136-RA,451,1533,C100.96,1243sc:136773-140035    0HX,3S,5L,
# ---
# > 2.Ncbi        XM_001949152.2  acyp2eg0029519t1,279,1424       XM_001949152.2,451,1562,nil     ACYPI005250-RA,277
# ,1303,C94.87,1243sc:127396-140035       0HX,1H,5L,

      if(not $ha and $eqOfId{$e} and (my $ealt=$ta{$e})) {  # ealt maybe id,id,id list : not now?
        # check goodness of alt equivalence
        # my $ec= $ec{$id.$e} || $ec{$e.$id} || 0;  # == "$pi,$lc"; pi == [CI]pc.px
        # my $eqgood= ($ec =~ /^[CI]/ or $ec > 66 or ($ec < 1 and $ec > 0.66)) ? 1 : 0;
        $ealt= 0 unless(equalhigh($id,$e));  # eq(id,e) exists here
        
        foreach my $ea (split",",$ealt) {
          next if($deq{$ea} or $did{$ea});
          if(my $halt= $ha{$ea}) { 
            ## warn "#DEBUG: alt match: id=$id, eq=$e (noval), eqalt=$ea, val=$halt \n" if $debug; # LOTs of
            #? if( equalhigh($id,$ea) ) ## eq($id,$ealt) may not exist
            { $ha=$halt; $e=$ealt; last; } 
          } } }
        
      #?? ^ wrong? alts are already lower scored alternates
      #.. maybe not, eq here are best match to id, not nesc best score among alts
      #.. do we need to check ta{e} always for best ha score? equal/*equal.ids has best agreement of ids ..
      #.. any cases of no ha{e}  but ha{ealt} ?
      
      next unless($ha); #..  we skip this model if no hxscore
      next if($deq{$e} or $did{$e}); # did already
      $deq{$e}= 1;
      
      $ha||= "0,0"; # not here now
      my $ec= $ec{$id.$e} || $ec{$e.$id} || "nil";  
      my $val= "$e,$ha,$ec";
      my $j=gnum($e);  
      
      ## joins: can be several j/e per d
      if($hx[$j]) { push(@{$many{$j}}, $val); $manyto1++; }
      else { $hx[$j]= $val; }
      $groups{$j}++;
      
      #OFF# $did{$e}++; # ? right, maybe not? 
      #  should defer did{id} to after classify choice? dont skip many-to-many links (eqe)
      #  when primary equivalence is not chosen?
      # instead use @idbest, mark all did and all direct links thru eq{idbest}
 
      ## Fixme1: when   e2 or e3 is join : 
      ##   classjoinsplit wont have d1..d2 data, and e2/e3 will show up again,
      ## answer: split out equivalent models and check, is this right?

      ## Fixme2: Oooopps, need to stop at 3rd level links, not traverse all overlaps in chain
      ## do eq{e} only for eq of eq{id} 
      
      ## Fixme3: add alts of eq also? eq may not be best score among alts, just best match to id
      ## alts causing problems: all are lower scoring/locus, should ignore all? but for cases
      ##  of equal.ids cross match to alts
      
      if( $eqOfId{$e} and ($eq{$e}) ) {  ##  or $ta{$e}
        my $eqa= $eq{$e}; ## join",",$eq{$e},$ta{$e}; 
        my @eqe= grep { not ($deq{$_} or $did{$_}) } split",",$eqa;
        push(@eq, @eqe) if(@eqe > 0); 
      }
    } 
    
    my($class,$ibest,$idbest,$rules,$scoremax);
    
    #? add singleton classify?, for separate rules?
    if(scalar(%groups) < 2) { ($class,$ibest,$idbest,$rules,$scoremax)= classify1(@hx); }
    elsif($manyto1) { ($class,$ibest,$idbest,$rules,$scoremax)= classifyJoinSplit(\@hx,\%many); }
    else { ($class,$ibest,$idbest,$rules,$scoremax)= classify3(@hx); }

    ## idbest for Noscore: choose one of Ncbi, Evigene if there, ignore ACYPI1 ?
    ## some of the Noscore's are data mistakes, esp missing expression scores
    
    my @idbest=split(/\|/, $idbest);
    my $idb= $idbest;
    if($manyto1) { 
      $idb= $idbest[0]; # better check $hx line for match
      foreach my $d (@idbest) { if($hx[$ibest] =~ /$d/) { $idb= $d; last; } } 
    }
    
    # update: did idbest, eq{idbest} here; this leaves undone tertiary links to best 
    # idbest == 0 ? use $id ?
    foreach my $d (@idbest) {
      $did{$d}++; 
      if($ta{$d}){ foreach my $a (split",",$ta{$d}) { $did{$a}++; } } # mark off alts of this best
      if($eq{$d}){ foreach my $e (split",",$eq{$d}) { $did{$e}++; # mark off direct equivalents
         ## .. this step has effects, good and bad, is where eq.id cross-linked alts can be removed,
         #     but also removes partial overlap 2nd locus to 2ndary alts
         ##  equal($a,$d) may not exist; test only equal(e,d) ? or test all ta{d},ta{e} ?
        # if($ta{$e}){ foreach my $a (split",",$ta{$e}) { $did{$a}++ if(equalhigh($a,$d)); } } # mark off alts of links?? is this ok? 

        if($ta{$e} and equalhigh($e,$d)){ foreach my $a (split",",$ta{$e}) { $did{$a}++; } } # mark off alts of links?? is this ok? 

         ## ^ IS problm when e is really low score alt, with part overlap of idbest, 
         ##   but e-alt is true high score at 2nd locus

        } }
    }
    
    # $hx[0]= $class; #? add $idbest to hx0 as new column2 ?
    $hx[0]= "$ibest.$class\t$idb"; # to pick col eg for Same2; oops Same ibest == 4
    
    ##push @hx, $rules if($debug);    
    print join("\t",qw(BestSet BestID Geneset1 Geneset2 Geneset3 RulesOfChoice)),"\n" unless($didhead++);
    print join("\t", @hx,$rules),"\n";  
    
    #* output for joinsplit:
    #* if join model is best, output all split models in same 1 row?
    #* if splits are best, one row per split model, join model in all?
    ##if($rules =~ /splitjoin=join|splitjoin=\w+join/) { #? also splitjoin=discrepancyjoin
    if($manyto1) { #? all? so we see all IDs
      my %dididb = ( $idb => 1);
      my $tag= ($rules =~ /splitjoin=\w*split/)?"s2":($rules=~/splitjoin=\w*join/)?"j1":"sj";
      for( my $k=0; $k < $manyto1; $k++) { 
        my @krow= @hx; my $kadd=0;
        for my $j (1,2,3) {
          if($many{$j} and $many{$j}->[$k]) { 
            $krow[$j]= $many{$j}->[$k]; 
            $kadd++;
            }
          }
        
        # FIXME: "#tag" problem now for splits with idb not same as last idb
        # my $idlast=$idb; # not good enough need hash of done
        #? need global hash? NO, 40819 ids, dups are '0': now have 44341 idbest rows (nocomment) but 40820 uniq ids
        
        $idb=$idbest[0];
        foreach my $id (@idbest) { if($krow[$ibest] =~ /$id/) { $idb= $id; last; } }
        $krow[0]= "$ibest.$class\t$idb";
        my $tag1= ($dididb{$idb}++) ?  "#$tag." : "";
        print $tag1,join("\t", @krow,$rules),"\n" if($kadd);  
        last unless($kadd);
      }
    }
    
  } 
} 



#-----------------------------------------------------------------------------


sub classify3
{
  my (@v)= @_;
  return classify0( "", $DH,$DX,$XpH,$XpHi,$HpX,$HpXi, @v);
}

sub classify1
{
  my (@v)= @_;  
  
  my $i=0; my $ng=0;
  for my $j (1,2,3) { unless($v[$j] eq "0") { $i=$j; $ng++; } } 
  return classify3( @v) if($ng>1);
  
  my($dc,$hm,$xm,$c,$s,$class,$rule,$Hsmall,$Xsmall)= (0) x 20; 
  my $rules="";
  
  ($dc,$hm,$xm)=split",",$v[$i]; 
  $c=$i;
  
  ## more options for this? for low-qual gene set
  if(grep { $_ == $c } @GNUM_NONE_SKIP) {
    $Hsmall= $H_SMALL3; $Xsmall= $X_SMALL3;
  } else {
    $Hsmall= $H_SMALL; $Xsmall= $X_SMALL;
  }
  
  $v[0] = $class = ($hm < $Hsmall and $xm < $Xsmall) ? "Noscore" : ($s > 1) ? "Same$s" : gname($dc); 

  if( $class =~ /^(Noscore|Same)/ and (grep { $_ == $c } @GNUM_NONE_SKIP)) {
    $dc = 0; $c= 0; 
  }
  return (wantarray) ? ($class,$c,$dc,$rules, "$hm.$xm") : $class;
}


sub classify0
{
  my( $splitsums, $DH,$DX,$XpH,$XpHi,$HpX,$HpXi,@v)= @_;
  my($dc,$hm,$xm,$c,$s,$class,$rule)= (0) x 10; 
  my $rules="";

# need another rule here for sum(splits) > join; require sum >> join, not just a bit >
# $SPLITGREATER = 1.75; ??  JOINGREATER = 1.25 ?? or 1.5?
# case XM_003240023.1: split 155,565 + 90.1,0 = 245.565  >= join 198,565 ; this join looks best despite splitsum
# also should do single model tests, if split(1) ~= join, sum-split should win

# case2 XM_001952703.2 where join is better than split parts, but sum-split > join:
# splits: 391,310 + 135,221 = 526.531 > 518,310 best model
  
  my %splitsums= ($splitsums) ? map{ $_,1 } split",",$splitsums : ();

  for my $i (1,2,3) { 
    my($d,$h,$x)=split",",$v[$i]; 
    if($d eq "0" or ($h <= 0 and $x <=0)) { $rule="5Z"; }
    elsif( $dc eq "0" and ($x>0 or $h>0)) { $hm=$h; $xm=$x; $c=$i; $dc=$d; $s=0; $rule="0HX"; } 
    else { 
      my $px= $x > $X_LARGE ? $xm/$x : ($xm > $X_LARGE) ? ($xm+$X_SMALL)/(abs($x)+$X_SMALL) : 1; 
      my $ph= $h > $H_LARGE ? $hm/$h : ($hm > $H_LARGE) ? ($hm+$H_SMALL)/(abs($h)+$H_SMALL) : 1;

  ## rule 6SJa: if ($i is sumsplit and $c is join) 
  #    if( iscore / cscore > SPLITGREATER ) ibest elsif( cscore / iscore < JOINGREATER ) ibest
  ##  and 6SJb: if ($c is sumsplit and $i is join) ... should this prefer join if sumsplit ~ same?
      if($c and $splitsums and ($splitsums{$i} or $splitsums{$c})) {

        my $SPLITBEST = (scalar(%splitsums) > 1) ? $SPLIT_MAJOR : $SPLIT_MINOR; 
        # lower criterion when splits are majority vs joins   

        my $sjrule=0; # 0 == Same
        if($splitsums{$i} and not $splitsums{$c}) { # 1 = isplit best; 3 = cjoin best
          if($h > $H_SMALL) { $sjrule= ($hm < $H_SMALL or $h/$hm > $SPLITBEST) ? 1 : 3; }  # ($hm/$h > 0.85) ?  1/1.25 = 0.80
          else { $sjrule= ($hm > $H_SMALL) ? 3 : 0; } 
          
        } elsif($splitsums{$c} and not $splitsums{$i}) { # 2 = ijoin best; 4 = csplit best
          if( $hm > $H_SMALL) { $sjrule= ($h < $H_SMALL or $hm/$h > $SPLITBEST)? 4 : 2; } # ($h/$hm > 0.85) ? 
          else { $sjrule= ($h > $H_SMALL) ? 2 : 0; } 
          
        } else { # both splits ..
          $sjrule= -1; # continue to other rules
        }
        if($sjrule >= 0) {
          $rule= "6SJ$sjrule"; ## ($sjrule==1)? "6SJa" : ($sjrule==2)?"6SJb" : "6SJ0";
          if($sjrule==1 or $sjrule == 2) { $hm=$h; $xm=$x; $c=$i; $dc=$d; $s=0; }
          $rules.= "$rule,";
          next;
        }
      }
  
      ## rule 4XH == rule1H + h small + x large = h/x discrepancy
      if($h > $hm + $DH and $px < $XpH and $ph > $HpX and $ph < $HpXi){ $rule="4XHa"; $hm=$h; $xm=$x; $c=$i; $dc=$d; $s=0;} # discrepency rule: xpress >> homol
      elsif($h > $hm + $DH and $px > $XpHi and $ph > $HpX and $ph < $HpXi) { $rule="4XHb"; } # no change

      ## rule1H = diff h large enough
      elsif( $h > $hm + $DH) { $rule="1H"; $hm=$h; $xm=$x; $c=$i; $dc=$d; $s=0; }  # best by homol

      ## rule2XH = diff x large enough, h not too small
      elsif( $x > $xm + $DX and $h > $hm - $DH) { $rule="2XH"; $hm=$h; $xm=$x; $c=$i; $dc=$d; $s=0;} # best by expression

      ## rule3S = diff x and h both small
      elsif( abs($x-$xm) <= $DX and abs($h - $hm) <= $DH ) { $rule="3S"; $s++;  } # same; $c=4 <NO ; make c == ibest list? 1,2
      ## rule5L = x and h too small
      elsif( ($h < $hm - $DH and $x <= $xm ) or ( $h <= $hm and $x < $xm - $DX) ) { $rule="5L"; }
      else { $rule="5e"; } ## $rule="5err"; # what?  $c=0 ?? maybe not error, but last is best
      } 
    $rules.= "$rule,";
  } 

#         # 4XHa isnt always discrepency; same as 2XH often == same ho, better x; put 2XH before this? or add $h > $hm + $DH
#       if($px < $XpH and $ph > $HpX and $ph < $HpXi){ $rule="4XHa"; $hm=$h; $xm=$x; $c=$i; $dc=$d; $s=0;} # discrepency rule: xpress >> homol
#       elsif( $px > $XpHi and $ph > $HpX and $ph < $HpXi) { $rule="4XHb"; } # no change

  
  $s++; 
  $v[0] = $class = ($hm < $H_SMALL and $xm < $X_SMALL) ? "Noscore" : ($s > 1) ? "Same$s" : gname($dc); 

  # drop this dc idbest if Noscore, then pick up other if there..
  # ?? avoid these also for Same score results
  if( $class =~ /^(Noscore|Same)/ and (grep { $_ == $c } @GNUM_NONE_SKIP)) {
    $dc = 0;  $c=0; #?
  }

    ## idbest for Noscore: choose one of Ncbi, Evigene if there, ignore ACYPI1 ?
    ## some of the Noscore's are data mistakes, esp missing expression scores
  if($dc eq "0") {
    # skip 3 (FIXME); always choose 2 if there; else choose 1 ?? choose/skip by gname() ?
    # there are problems w/ this, many Ncbi RefSeq w/ no score are really poor acypi1 models, should drop
    # 200 of 300 Noscore Ncbi-only are LOCnnnn
    for my $i (@GNUM_NONE_KEEP) { my($d,$h,$x)=split",",$v[$i]; $dc=$d unless($d eq "0"); }
  } 
  
  #? maybe return max scores as $hm/$xm ? or $hm.$xm ?
  return (wantarray) ? ($class,$c,$dc,$rules, "$hm.$xm") : $class;
}


# use classify3, but add 
# a. majority vote when 1=join vs 2,3=split  OR 1,2=join vs 3=split
# b. sum scores for splits > join scores

## see below, replace hxscore for joinsplit use
#  simple sum? weigthed?  h*2 + x ? or more like classify3 ?
sub hxscore { my $v=shift; my($d,$h,$x)=split",",$v; return $h*2 + $x; } 

sub classifyJoinSplit
{
  my ($vref,$manyref)= @_;
  my @vone= @$vref;
  my %many= %$manyref;
  my @ngene=(0) x 4; my $maxgene=0; my @hxval=(0) x 4; 
  my @hsum=(0) x 4;  my @xsum=(0) x 4;  my @dsum=(0) x 4;
  for my $i (1,2,3) { 
    
    # my $hxval= hxscore($vone[$i]);
    my ($dv,$hv,$xv)= split",",$vone[$i];  
    next if( $dv eq "0");    
    
    $ngene[$i] = ($vone[$i]) ? 1 : 0;
    $dsum[$i]= $dv; # one of many? make list?
    my $hvs = $hv; my $xvs = $xv; 
    my $hxval= $hv*2 + $xv;
    my ($hxtop,$ktop)= ($hxval,0);
    if($many{$i}) { 
      # problem of many score >> primary score, problem using classify3 w/ only primary score
      # if $hxval2 > $hxval1, swap many{i}[j] for vone[i] ? YES
      my @ev= @{$many{$i}};
      $ngene[$i] +=  scalar(@ev);
      my $k=0;
      foreach my $ev (@ev) { 
        # my $hxval2= hxscore($ev); 
        my ($dv2,$hv2,$xv2)= split",",$ev; 
        my $hxval2= $hv2*2 + $xv2;
        $dsum[$i] .='|'.$dv2; # will this work? fix caller to parse
        $hvs += $hv2; $xvs += $xv2;
        $hxval+= $hxval2; $k++;
        if($hxval2 > $hxtop) { $hxtop=$hxval2; $ktop=$k; }
        }
      
      # problem  hxscore() not same as classify3 scoring; use that to check for swap?
      if($ktop>0) {
        my $ev=$ev[$ktop-1]; my $ov=$vone[$i]; 
        $vone[$i]= $ev; $ev[$ktop-1]= $ov; 
        # this swap needs to propogate back to $vref,$manyref for output
        $manyref->{$i} = \@ev;
        $vref->[$i]= $ev;
      }
    }
    $hsum[$i]= $hvs; $xsum[$i]= $xvs;
    $hxval[$i]= $hxval;
    $maxgene= $ngene[$i] if($maxgene < $ngene[$i]);
  } 
  
  #ngene[1] always == 1? NO, now 1 can have splits;  what of ngene[2,3] > ngene[1] but not same? skip?
  my (@im); my $i1=0; # note @im is majority indices, $i1 minority, if this is simple
  my($class,$ibest,$idbest,$rules,$scoremax,$vmajor,$vminor,$ok)=(0) x 20;
  $rules="";


  for my $k (1,2,3) {
    if($ngene[$k] == 1) {
      my($k1,$k2)= grep{ $k != $_ } (1,2,3);
      if($ngene[$k1]>1 and $ngene[$k2] > 1) { $i1=$k; @im=($k1,$k2); }
      elsif($ngene[$k1] == 1 and $ngene[$k2] > 1) { $i1=$k2; @im=($k,$k1); }
      elsif($ngene[$k1] > 1 and $ngene[$k2] == 1) { $i1=$k1; @im=($k,$k2); }
      last if($i1 > 0);
    } 
  }

   ## change this from == to >= ?? problems otherwise 
  if($ngene[1] == $ngene[2] and $ngene[1] == $ngene[3]) { @im=(1,2,3); $i1= 0;} # is this possible? yes
  
#   elsif($ngene[1] == $ngene[2]) { @im=(1,2); $i1= 3;}
#   elsif($ngene[1] == $ngene[3]) { @im=(1,3); $i1= 2;}
#   elsif($ngene[2] == $ngene[3]) { @im=(2,3); $i1= 1;}
  
  if($i1 == 0 and @im == 0) {
    # test "majority" with only two forms? { @im=(2); $i1= 1;}
    # if($ngene[3] == 0 and $ngene[2]>0 and $ngene[1] > 0 and $ngene[2] != $ngene[1]) { $i1=1; @im=(2); }
    for my $k (1,2,3) {
      if($ngene[$k] == 0) {
        my($k1,$k2)= grep{ $k != $_ } (1,2,3);
        if($ngene[$k1]>0 and $ngene[$k2] > 0 ) { $i1=$k1; @im=($k2); } #? and $ngene[$k1] != $ngene[$k2]
        last;
      }
    }
  }

  my $hasmajority= (@im > 1 and $i1 > 0) ? 1 : 0;
  
use constant REDO_JCTEST => 1;
if(REDO_JCTEST) {  
  # redo using hsum, xsum and classify3
  # now some problems w/ messier acypi1[3] splits/overlaps turning up more when dont want
  # .. part of problem is many acypi1 models overlap a lot, if clean those out this would work better
  # e.g.  acyp2eg0000501t2 acyp2eg0000539t1
  # .. this is doing ok now after removing acypi1 overlaps, BUT, splits too often when sum(2) > 1 by a bit
  # ..  two short genes should need higher sum than 1? adjust classify3 sets
  
  if(@im > 0) {
  
    ## first test best single models; if ho(best-split) < 0.75 ho(join), skip sum-test
    ($class,$ibest,$idbest,$rules,$scoremax)= classify3(@vone);
    
    # my $spljn= ($ibest==0) ? "eithersj": ($ngene[$ibest] > 1) ? "split" : "join" ;  
    # if($spljn =~ /join/) 
    my $scoresplit=0; 
    if($ngene[$ibest] == 1) { ## join
      for my $i (1,2,3) { 
        next if($i==$ibest);
        my($d,$h,$x)=split",",$vone[$i]; 
        $scoresplit=$h if($h>$scoresplit);
      }    
    }
    
    ## my $sjfactor= ($ibest == 3) ? 0.66 : 0.75; # HACK fix for poor qual acypi1: 
      # dont want too many of it, either split or join unless it agrees w/ other or has >> value
      # or use ($ibest == $i1) == minority result
    
    my $sjfactor= ($ibest == $i1 and $hasmajority) ? $SPLIT_NEAR_JOIN_MAJOR : $SPLIT_NEAR_JOIN; # minority has smaller chance of being right  
      
    if($scoresplit == 0 or $scoresplit > $sjfactor * $scoremax ) {
      ## not quite right here; for homology, sum(split) should be >> join,
      ## two lesser parts should not beat large gene w/ > ho, unless sum(split) >> join
      my @vsums= (0) x 4; for my $i (1,2,3) { $vsums[$i]= join(",",$dsum[$i],$hsum[$i],$xsum[$i]); }
      my $isplits=""; for my $i (1,2,3) { $isplits .= "$i," if($ngene[$i] > 1); }
      # modify classify splits for hasmajority == join or split
      
      my $k= 1; ## ($maxgene>3) ? 2.5 : ($maxgene>2) ? 2.2 : 2;
      my($kDH,$kDX)= ($k*$DH, $k*3*$DX); # should depend on many count?
      
      ($class,$ibest,$idbest,$rules,$scoremax)= classify0($isplits,$kDH,$kDX,$XpH,$XpHi,$HpX,$HpXi,@vsums);
    }
    
    # if ibest has many, need idbest as list of ids : dsum = id1|id2|id3..

    # if(@im == 3) ??
    if(@im == 2 or @im == 1) {  
      if($ibest == $i1) { 
        # my $i2= $im[0]; # or 1
        my $i2= ($hsum[$im[1]]+$xsum[$im[1]] > $hsum[$im[0]]+$xsum[$im[0]]) ? 1 : 0;
        $vminor= $scoremax; 
        $vmajor= "$hsum[$i2].$xsum[$i2]";
        }
      elsif($im[0] == $ibest or $im[1] == $ibest) { 
        $vmajor= $scoremax; 
        $vminor= "$hsum[$i1].$xsum[$i1]";
        }
    }
   
  } else {

    ($class,$ibest,$idbest,$rules,$scoremax)= classify3(@vone);
    
    if(@im == 2 or @im == 1) {  # dont get here now.
      map{ $vmajor = $hxval[$_] if($hxval[$_]>$vmajor); } @im; 
      $vminor= $hxval[$i1]; 
    }
  
  }
  
} else {

  ($class,$ibest,$idbest,$rules,$scoremax)= classify3(@vone);
  if(@im == 2 or @im == 1) { # was (@im == 2)
    map{ $vmajor = $hxval[$_] if($hxval[$_]>$vmajor); } @im; 
    $vminor= $hxval[$i1]; 
  }
}

  $rules.="sjscore=$scoremax,sjbest=$idbest,";
  
  my $spljn= ($ibest==0) ? "eithersj": ($ngene[$ibest] > 1) ? "split" : "join" ; #? always ibest, can be 0 
  
  # ?? favor join over split if either/nodecision ? at least if split is to too-small genes
  # FIXME: discrepancies should be resolved, in favor if largest vmajor/vminor score?
  
  if(@im == 3) {
    $rules.="splitjoin=all$spljn,"; $ok=1; # use $ibest
  
  } elsif($vminor == 0 and $vmajor == 0) { 
    $rules.="splitjoin=nodecision$spljn,"; 
    
  } elsif($vminor > 1.01 * $vmajor) {  # always $i1 or simple score?
    if($i1 == $ibest) { $ok=1; $rules.="splitjoin=$spljn,"; }
    else { $rules.="splitjoin=discrepancy$spljn$i1,"; }
    
  } elsif($vmajor > 1.01 * $vminor) { # best of @im ?
    if($im[0] == $ibest or $im[1] == $ibest) { $ok=1;  $rules.="splitjoin=$spljn,"; }
    else { $rules.="splitjoin=discrepancy$spljn".join("/",@im).","; }
  
  } else { # split decision, use simple score to choose best, 
    $ok=1; $rules.="splitjoin=either$spljn,";
  }
  
  return (wantarray) ? ($class,$ibest,$idbest,$rules,$scoremax) : $class;
}

__END__

=item 2 choice join/spit error

problem in splitjoin: these NCBI XM's are true 1-exon gene cassette of 15+ genes, but 1 evigene model stringing
   them together was chosen, has same ho score, but > single x score (15x)

- part of error is equal.ids has not linked all 15 ncbigenes to 1 evigene, though cds overlaps
  http://192.168.2.2:8091/gbrowse/cgi-bin/gbrowse/aphid2x/?name=
  Scaffold1438:1..9000

- each ncbi cds-exon/gene has nearly same ho score = 100..105, 1 bit less than evigene 14-exon
- express is spotty, but is roughly 1/5 of exon span (504 bp)
   -- evigene model should have xpres=sum(15 x 100) but has > 2x that: this is main choice problem?
   -- however classifyJoinSplit() should recognize that 1 gene ho=100 > 15 gene x ho=100 is wrong choice
      -- even case here of 1 gene > 5 gene w/ same ho is bad
- need to redo hxscore() : calc hsum, xsum separately and use classify3 algorithm?


yeast.% egrep 'XM_003248061.1|XM_003248062.1|XM_003248063.1|XM_003248064.1|XM_003248065.1|XM_003248066.1|XM_003248067.
1|XM_003248068.1|XM_003248069.1|XM_003248070.1|XM_003248071.1|XM_003248072.1|XM_003248073.1|XM_003248074.1|acyp2eg0030
319t1' compare3-hoex7.tab12

]1.Evig acyp2eg0030319t1        acyp2eg0030319t1,106,3889       XM_003248067.1,102,486,17.17,1438sc:1-8635      0       0HX,5L,5Z,splitjoin=join,
#j1.1.Evig      acyp2eg0030319t1        acyp2eg0030319t1,106,3889       XM_003248061.1,101,165,5.50,1438sc:1-8635       0       0HX,5L,5Z,splitjoin=join,
#j1.1.Evig      acyp2eg0030319t1        acyp2eg0030319t1,106,3889       XM_003248063.1,104,172,29.29,1438sc:1-8635      0       0HX,5L,5Z,splitjoin=join,
#j1.1.Evig      acyp2eg0030319t1        acyp2eg0030319t1,106,3889       LOC100570635:unknown_transcript_1,0,24,0.92,1438sc:1-8635       0       0HX,5L,5Z,splitjoin=join,
#j1.1.Evig      acyp2eg0030319t1        acyp2eg0030319t1,106,3889       XM_003248072.1,105,144,10.10,1438sc:1-8635      0       0HX,5L,5Z,splitjoin=join,
2.Ncbi  XM_003248062.1  0       XM_003248062.1,101,52   0       5Z,0HX,5Z,
2.Ncbi  XM_003248064.1  0       XM_003248064.1,105,150  0       5Z,0HX,5Z,
2.Ncbi  XM_003248065.1  0       XM_003248065.1,104,105  0       5Z,0HX,5Z,
2.Ncbi  XM_003248066.1  0       XM_003248066.1,98.6,134 0       5Z,0HX,5Z,
2.Ncbi  XM_003248068.1  0       XM_003248068.1,105,55   0       5Z,0HX,5Z,
2.Ncbi  XM_003248069.1  0       XM_003248069.1,102,486  0       5Z,0HX,5Z,
2.Ncbi  XM_003248070.1  0       XM_003248070.1,100,149  0       5Z,0HX,5Z,
2.Ncbi  XM_003248071.1  0       XM_003248071.1,105,110  0       5Z,0HX,5Z,
2.Ncbi  XM_003248073.1  0       XM_003248073.1,102,-63  0       5Z,0HX,5Z,
2.Ncbi  XM_003248074.1  0       XM_003248074.1,89.0,91  0       5Z,0HX,5Z,

egrep 'XM_003248061.1|XM_003248062.1|XM_003248063.1|XM_003248064.1|XM_003248065.1|XM_003248066.1|XM_003248067.
1|XM_003248068.1|XM_003248069.1|XM_003248070.1|XM_003248071.1|XM_003248072.1|XM_003248073.1|XM_003248074.1|acyp2eg0030
319t1' equal/ncbigene2.overevigene8f.tab5
LOC100570635:unknown_transcript_1       noid    acyp2eg0030319t1/0.92   1438sc:9-434
XM_003248061.1  noid    acyp2eg0030319t1/5.50   1438sc:436-939
XM_003248062.1  noid    na      1438sc:941-1444
XM_003248063.1  noid    acyp2eg0030319t1/29.29  1438sc:1446-1949
XM_003248064.1  noid    na      1438sc:1951-2454
XM_003248065.1  noid    na      1438sc:2456-2959
XM_003248066.1  noid    na      1438sc:2961-3464
XM_003248067.1  noid    acyp2eg0030319t1/17.17  1438sc:3466-3969
XM_003248068.1  noid    na      1438sc:4472-4975
XM_003248069.1  noid    na      1438sc:4977-5480
XM_003248070.1  noid    na      1438sc:5482-5985
XM_003248071.1  noid    na      1438sc:5987-6490
XM_003248072.1  noid    acyp2eg0030319t1/10.10  1438sc:6492-6995
XM_003248073.1  noid    na      1438sc:6997-7500
XM_003248074.1  noid    na      1438sc:8008-8484

egrep 'XM_003248061.1|XM_003248062.1|XM_003248063.1|XM_003248064.1|XM_003248065.1|XM_003248066.1|XM_003248067.
1|XM_003248068.1|XM_003248069.1|XM_003248070.1|XM_003248071.1|XM_003248072.1|XM_003248073.1|XM_003248074.1|acyp2eg0030
319t1' equal/evigene8f.overncbigene2.tab5
acyp2eg0030319t1        AUGepir2s1438g21t1      XM_003248063.1/29.29,XM_003248067.1/17.17,XM_003248072.1/10.10,XM_003248061.1/5.50,LOC100570635:unknown_transcript_1/0.92       1438sc:1-8635


#..............................

=item compare3 work 22 july 2011
# compare3 work, 22 July 2011

1. resolve express scoring
    -- rnas/genescore/... from genome mapped reads, give usable span but poor at wrong-gene detection
        -- check on neg scores for end-exons == revstrand matches, promote these errors as false UTRs?
    -- need to check trgsnap scores of wrong-tr/perfect ratio for use to detect false UTRs
        -- need to remove alt-tr effects: rescore from gsnap.bam removing alt-tr ? (but need 2.best-scored alttr)
        -- can only use this as model err for NON-paralog, NON-alttr  **
        
2. need to mark all model alt-tr and pick best isoform from ho+ex scoring
   -- have a few problems for multiple models/locus not marked as alt-tr, e.g. userchoice:

3. redo eqid.equal.hx from 1,2, using only best isoform/predictor/locus

4. redo compare3-ho4.best.tab from above, using ho+ex best score algorithm 
    -- how weight ho vs ex? discrepency if disagree much?
    -- allow diff(ex) to trump ho if diff(ho) is small, diff(ex) larger ?

5. public tables:
	pubdir/genes-bestof3/
	5.1: locus equivalences:
		equal/acypi1-evigene8f.equal.ids      equal/ncbigene2-acypi1.equal.ids      equal/ncbigene2-evigene8f.equal.ids

	5.2: gene model scores/locus:
		rnas/tracypi1gene.hoex.tab   rnas/treg8gene.hoex.tab      rnas/trncbi2gene.hoex.tab
    # rename transcript-XXX-hxscore.txt
    
	5.3: final best model choice table:
		compare3-hoex6.tab2 == aphid2-best-of-evigene2-ncbi2-acypi1-genemodels.txt
	
	5.3.2: other summaries?  
		gbrowse-map.html for best vs alt models?
	
	5.4: method document, scripts for producing these
	
	5.5: any other intermediate tables?
	
#---
#.. redo here, merge express + ho scores and pick best alt-tr for each gene set
-- merge using equal.ids, print all genes in each group.hoex.tab, 3 x 3 column  as per compare3.
   in= tr*gene.hoex.tab, ../equal/*.equal.ids
   
   foreach g (genes{eg}) 
   	 ng=eqnc{$g}; ag=eqac{$g}; 
   	 print hoex{g}, hoex{ng}||NA, hoex{ag} || NA
   	 hoex{g} = hoex{ng} = hoex{ag} = done
   foreach g (genes{nc}) 
   	..
   foreach g (genes{ac})
   	..

# FIXME: next if($did{$d}++); << bad when gene d is join vs split of other 2; need did{1,2,3} ?

cat equal/*.equal.ids | sed 's/^/eq /' | cat - rnas/tr*.hoex.tab | perl -ne\
'chomp; if(s/^eq //){ ($d,$e,$pi,$lc)=split"\t"; $ec{$d.$e}="$pi,$lc"; 
$eq{$d}.="$e," unless($d eq "na"); $eq{$e}.="$d," unless($e eq "na");  } 
elsif(/\t(h|x|a)=/){ s/\t(h|x|a)=/\t/g; ($d,$h,$x,$a)=split"\t"; $gg=gnum($d); 
$ha{$d}= $hx{$gg}{$d}="$h,$x";   map{ s,/.*,,; $ta{$_}=$d if(/\w/); } split",",$a; } 
sub gnum { local $_=shift; return (/acyp2/)?1:(/ACYP/)?3:(/^na/)?4:(/^[XNLYG]/)?2:4; } 
sub _gsort{$ai=gnum($a); $bi=gnum($b); return($ai <=> $bi or $a cmp $b); } 
END { foreach $gg (sort keys %hx) { foreach $d (sort keys %{$hx{$gg}}) { 
next if($did{$d}++); @hx=(0) x 4; $i=gnum($d); $hx[$i]= "$d,".$hx{$gg}{$d}; 
$eq=$eq{$d}; @eq= ($eq) ? (sort _gsort split",",$eq ) : (); 
foreach $e (@eq) { $ha=$ha{$e}; if(not $ha and $ea=$ta{$e}) { $ha=$ha{$ea}; $e=$ea; } 
$ha||= "0,0"; $ec= $ec{$d.$e} || $ec{$e.$d} || "nil";  
$j=gnum($e); $hx[$j]= "$e,$ha,$ec"; $did{$e}++; } 
print join("\t", @hx),"\n"; } } }' 
 > compare3-hoex6.tab

#-----

=item 2 problems

2. problems  w/ join/split and multiple genes in 1 group at same locus (evigene esp).

  compare3-hoex6.tab2:
	Acy1	acyp2eg0000128t1,107,1182	XR_118890.1,na,1254,0.64,1sc:1417946-1421367	ACYPI47564-RA,195.9,561,25.45,1sc:1417946-1421367
					^^ superceded by acyp2eg0037009t1
	Evig	acyp2eg0037009t1,345,1755	XR_118890.1,na,1254,0.45,1sc:1417439-1422964	ACYPI47564-RA,195.9,561,43.37,1sc:1417439-1422964

  -- disagreement on join/split, use 2 vs 1 to decide?
  -- ? use rule to split where 2 types disagree w/ 3rd, AND sum of ho+ex over 2 loci > than 1 ho+ex of 3rd?
  -- adding intron scores would help some decisions

# >> need to sort by geneset ID for this to work right
# >> Evigene lacks all but few, Ncbi, Acyp1 have 800-1000 joins vs 2 splits: is this artifact of table maker?

cat compare3-hoex6.tab2 | sort -k3,3 | perl -ne\
'$gl=$_; ($c,$eg,$nc,$ac)=split"\t"; @eg=split",",$eg; @nc=split",",$nc; @ac=split",",$ac; 
@gg=($eg[0],$nc[0],$ac[0]); print $lg,$gl if($lg and ($lg =~ /^$C/ or $gl =~ /^$C/) 
and $gg[$i] eq $lgg[$i] and $gg[$j] ne $lgg[$j] and $gg[$k] ne $lgg[$k]) ; $lg=$gl; @lgg=@gg; 
BEGIN{$k=2;$j=0;$i=1; $C="Ncbi";}' 

# BEGIN{$i=2;$j=0;$k=1; $C="Acy1";}' 
# BEGIN{$k=2;$j=0;$i=1; $C="Ncbi";}' 
# BEGIN{$k=2;$i=0;$j=1; $C="Evig";}' 


#.... messy decision, staggered overlapping models

 cat equal/*.equal.ids  | egrep 'acyp2eg0000366t1|XM_003240054.1|ACYPI007583|XM_003240078.1|acyp2eg0000363t1|ACYPI001900'
acyp2eg0000363t1        ACYPI001900-RA  17.18   2sc:1747687-1777487
acyp2eg0000363t1        ACYPI007583-RA  53.43   2sc:1747687-1777487
acyp2eg0000364t1        ACYPI001900-RA  17.56   2sc:1747726-1765548
acyp2eg0000366t1        ACYPI007583-RA  9.49    2sc:1768296-1776786
ACYPI001900-RA  XM_003240078.1  29.35   2sc:1747740-1751884
ACYPI007583-RA  XM_003240054.1  18.18   2sc:1762845-1769660
ACYPI007583-RA  XM_003240078.1  63.57   2sc:1762845-1769660
acyp2eg0000363t1        XM_003240054.1  37.44   2sc:1747687-1777487
acyp2eg0000363t1        XM_003240078.1  57.52   2sc:1747687-1777487
acyp2eg0000364t1        XM_003240078.1  5.61    2sc:1747726-1765548
acyp2eg0000366t1        XM_003240054.1  0.54    2sc:1768296-1776786
acyp2eg0000366t1        XM_003240078.1  9.25    2sc:1768296-1776786

cat compare3-hoex7.tab | egrep 'acyp2eg0000366t1|XM_003240054.1|ACYPI007583|XM_003240078.1|acyp2eg0000363t1|ACYPI001900'

Evig    acyp2eg0000363t1,712,2846       XM_003240054.1,226,1158,37.44,2sc:1747687-1777487       ACYPI001900-RA,88.2,396,17.18,2sc:1747687-1777487       
      0HX,5L,5L,splitjoin=either,
 ^^ does look like single best in this span; 
    acyp2eg0000366t1 is bad model, part of acyp2eg0000363t1
    acyp2eg0000364t1 is bad model, part of acyp2eg0000363t1
    
Acy1    acyp2eg0000366t1,0,1300 XM_003240054.1,226,1158,0.54,2sc:1768296-1776786        ACYPI007583-RA,487,1145,9.49,2sc:1768296-1776786        
    0HX,1H,1H,splitjoin=discrepancyjoin2,
 ^^ overlaps acyp2eg0000363t1, as well as acyp2eg0000366t1
 
XM_003240078 +  XM_003240054 are parts of acyp2eg0000363t1,  all of same exons
ACYPI001900 + ACYPI007583 are parts of acyp2eg0000363t1, some of same exons


=cut
