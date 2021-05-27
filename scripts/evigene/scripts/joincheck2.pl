#!/usr/bin/env perl
# joincheck.pl : check gene annots for join evidence

=item about

  joincheck.pl : check gene annots for join evidence
  redo and/or add, ovrna=, ovpro=
  
  test using existing rseq, intr exon annotations.  
      intr=count,N1,N2,.. < intron chains b/n exons
      rseq=obase/nbase,rid1,rid2 < rna id chains b/n exons

=item annotate before

  change -mark ovrna,ovpro to ocrna,ocpro
  
  #... all evg11e for joincheck
  
  ingff=nvit2_evigenes.pub11e.gff
  rnaset=$workd/genes/alts/*.analt.gff
  proset=$workd/prot/nasvit1-{apisbestpro_exonr4j,uparp_hym_exonr5.best9q}.gff.gz  

  onam=`echo $ingff | sed 's/.gz//; s/.gff//;'`
  
  cat $rnaset | grep CDS |\
  $evigene/scripts/overlapfilter -pass CDS -strand -over stdin -in $ingff -act markidbase -mark ocrna -pct 80 \
  > $onam.ja1.gff
  
  # redo prot using only mrna alignx > 65% : reduce to >39
  
  gzcat $proset | perl -ne\
  'if(/\tmRNA/){($pi,$ps)= m/alignx=(\d+)..(\d*)/; $ps=$pi if($pi>$ps); $p=($ps>39)?1:0; } elsif(/\tCDS/ and $p) { print; }' |\
  $evigene/scripts/overlapfilter -pass CDS -strand -over stdin -in $onam.ja1.gff -act markidbase -mark ocpro -pct 80 \
  > $onam.ja2.gff
  
  # reannotate intr using only CDS, or fiddle joincheck to use exon intr= + CDS ocpro/ocrna
  $evigene/scripts/overlapfilter.perl -intron2splice=error -pass 'CDS,intron' -act markid -midtype scoresum \
  -mark intr -over $workd/intron/intron_all.gff.gz -in  $onam.ja2.gff \
  >  $onam.ja3.gff
  
  cat $onam.ja3.gff | $evigene/scripts/joincheck.pl -overtag 'ocrna,ocpro' -intag intr -exon CDS > $onam.joincheck

=item splits

  can this be improved to detect split errors also?
    - check for chain evidence b/n genes, simplest if they come in loc-sorted
    
=cut

use strict;
use warnings;
use Getopt::Long;

my $exontype='CDS'; # 'exon' ... both? or not

my %vcode=( 
    0 => 'N0',   # no chain evidence
    1 => 'C1',   # complete
    2 => 'J2',   # join
    3 => 'C1T1', # complete but for terminal exon
    4 => 'IJ4',  # join split by intron only ; EP3 weak evd
    -1 => 'OV1', # overlapped chain
    -2 => 'ER2', # other error
    -3 => 'EP3', # poor/weak evidence
    -4 => 'EM4', # endmiss, ok 
    );

my %scorecode=( # +100..0..-100 score for genescores; to be weighted to keep/remove models
    0 => 0, #'N0',   # no chain evidence
    1 => 99, #'C1',   # complete : good
    2 => -99, #'J2',   # join : bad
    3 => 69, #'C1T1', # complete but for terminal exon : good
    4 => -9, #'IJ4',  # join split by intron only ; EP3 weak evd
    -1 => -9, # 'OV1', # overlapped chain
    -2 => -2, #was -39, #'ER2', # other error; mostly poor models  << SCORE TOO neg strong? score as EP3? -2 ; kicking out good models on weak evidence
    -3 => -1, # 'EP3', # poor/weak evidence; worse than no evd?
    -4 => 9, #'EM4', # endmiss, ok 
    );
           
my @overtag= qw( rseq ocpro ocrna); # change tags ovpro/rna to ocpro/rna to distinguish CDS annots
my $overtag=""; # "rseq ovpro ovrna"
my $intag="intr";
my $OVFIRST_ONLY= 0; # "ovpro"; ## 1; # opt?

my($gid, $nx, %gr, %gx, %gloc, $InQual); #%gi, 

my $optok= GetOptions(
  "exontype=s", \$exontype, 
  "intag=s", \$intag, 
  "overtag=s", \$overtag, 
  "firstover=s", \$OVFIRST_ONLY, 
);

if($overtag) { @overtag= split /[\s,;]+/, $overtag; }
unshift(@overtag, $intag) if($intag);

sub printit {
  my( $gid, $code, $nexon, @splits)= @_;
  # (0, $nx)
  # (1, $nx, $icx1, chainids($cmax))
  # (-3, $nx, $icx1, chainids($cmax))
  # (2, $nx, $icx1, @icxx, chainids($cmax, @cmore))
  my $vcode= $vcode{$code}||$code;
  my $score= $scorecode{$code}||$code;  
  print join("\t",$gid, $vcode, $score, $nexon, @splits, "intr=$InQual"),"\n" if($gid);
}

while(<>) {
  next unless(/^\w/);
  
  # gene-sorted recs; dont need gid/pid hashing
  if(/\tmRNA/) {  # dont need now? use only exons
  
    if($gid and %gx) { 
      printit( $gid, hasSplit()); # globals gx,gi,gr
      }
      
    ($gid)= m/ID=([^;\s]+)/;
    $nx= 0; %gloc= %gr= %gx= ();  # %gi= 
    $InQual=""; 
    
#     if( my($cw,$xw)= m/cxlen=(\d+).(\d+)/) { 
#       my $pcds= ($xw>0) ? int(0.5 + 100*$cw/$xw) : 100; 
#       }

  } elsif(/\t($exontype)/) { # assumes now is gene record sorted
  
    my($pid)= m/Parent=([^;\s]+)/ or next; 
    if($pid ne $gid) { warn "# Gene sort error: $pid ne $gid\n"; die; } # error
    my($tb,$te)=(split"\t")[3,4];
    
    # change {nx} back to {tb} and keep {tb}=te ; then can mix CDS, exon annots
    ++$nx;
    my $ix= $nx; # or $tb;
    $gx{$ix}= $nx; #  use tb or nx as index?
    # $gx{$tb}=[$te,$nx];  
    
    # my($in)= m/;$intag=[^,]+,(N[^;\s]+)/; 
    # $gi{$ix}= $in||""; # change to gr{$intag}

    ## FIXME: partial, 2ndary ovpro is screwing this up by paralog partial matches
    # .. m6AUGytlarvae51hrmp14c4g51t1 = join of 2 paralogs, both w/ some ovprot=q9ABCCB_HUMAN but better other splits
    # my($ra)= m/;rseq=[^,]+,([^;\s]+)/; 
    
    my $ra="";
    foreach my $tag (@overtag) {  # now includes intag
      ##my($rv)= m/;$tag=[^,]+,([^;\s]+)/; 
      my($rv)= m/\b$tag=([^;\s]+)/; # include score
      $gr{$tag}{$ix}= $rv||""; #?? use this
      }
  }
}

printit( $gid, hasSplit()); # globals gx,gi,gr


#==========

=item hasSplit

  gene-sorted recs; dont need gid/pid hashing
  BUT maybe want to save chains b/n genes for needJoin check
  
  ** Revise: use some variant of majority vote so that single false join doesn't throw out strong split evidence
  major problem is that 1 bad evidence join drags all others with it
  -- intron chains: single inner miss that is longish, when all outer introns found, is split evid.
  -- ovpro especially, use major vote (or 1.split vs 2.join better decision)
     based on how many CDS are matched by 1.splits vs 2.join
  -- ovrna: split == partial asmrna often, weight lower than ovpro, dont use this alone to split gene?

  # new

      in: for all exons, if( $nexon > 5 and $ninsplice >= 2*$nexon-4 ) .. stronger intr evid
      pro: count #cds exons matched per pro id, if proid1 + proid2 == nexon, stronger pro evid
           conversely if id1 joins parts but matches only few exons, weaker evid

  # new algo
    for gene :
      for x exons : 
        for e evtype (intr,ovpro,ovrna..) :
          chain{e,id}, tchain{e,x}= join exonchains
      end x    
    
    join tchain{intr,x} : join intron chain links (exon-exon only 2 part chains)
    
    chainsizes{tchains}= count','
    
    if max(chainsize)= all-exons : return complete
    else
      if max[1,2,3..] = all-exons : test join or split of max-n chains
    
    ** classify each id-chain like  evigene/scripts/intronscore.pl
      complete, completeinner, missinner, etc.
       then classify gene from best chains per evdtype
       -- allow intron chain to join ovrna, not ovpro?
       -- use some form of majority vote for ovpro ids, but not intron, ovrna?
          ie 2 split ovpro (2 species?) overrides 1 join ovpro
=cut

sub hasSplit
{
  my @xi= sort{ $a<=>$b } keys %gx; # change from gx{xi} to gx{xb}; use values for xi ?
  my $nx= @xi;
  my ($ic,$ix)=(0,0); 
  my (%chain, %xchain, %ochain, %tchain, %chain2mainchain, %chainequal);
  our %ichain=();
  
    ## FIXME: partial, 2ndary ovpro is screwing this up by paralog partial matches
    # .. use only 1st ID of list? esp for ovpro
    # .. m6AUGytlarvae51hrmp14c4g51t1 = join of 2 paralogs, both w/ some ovprot=q9ABCCB_HUMAN but better other splits

  my ( $inidchain, $intchain, $inscore, $inqual)= chainIntrons(); # using globals %gx,%gr,...
  $InQual= $inqual; # global hack for test output
  
  # if($inqual =~ /complete/) .. likely complete
  # if($inqual =~ /miss_inner/) .. suspect split

# fixme: intron=complete chain overrides ovpro,ovrna err:
# Nasvi2EV000198t1        ER2     -39     7       1,2,3,4,5,7,    6,      apis2gno_185004,; == nvit1v4aSallLoc530t4,
#        intr=complete

  my @tagnoi= grep{ $_ ne $intag } @overtag;
  
  foreach my $x (@xi) {
    #my($xe,$xi)= @{$gx{$x}};
    $ix++; $ichain{$x}= $ix;
    
    #old# my $in= $gi{$x}; my @in= split",",$in; 
    # my $ra= $gr{$x}; my @ra=split",",$ra;
    
    my @ra=(); 
    foreach my $tag (@tagnoi) { 
      my $rt= $gr{$tag}{$x}; 
      my ($rval,@rt)= split",",$rt;
      if(@rt) { if($OVFIRST_ONLY and $tag =~ /$OVFIRST_ONLY/) { push @ra, $rt[0]; } else { push @ra, @rt; } }
    }
    
    #v2: separate parse of @in, @pro, @rna
    #   in: for all exons, if( $nexon > 5 and $ninsplice >= 2*$nexon-4 ) .. stronger in evid
    my %xc=();
    foreach my $id (@ra) { # @in, 
      my $c= $chain{$id};  
      # $c= ++$ic unless($c); # change to $c= $id ??
      unless($c){ $c= $id; } # chainid: prepend tag/type sort char? P, R, si
      $chain{$id}= $c; $xc{$c}++;
      } 
    
    if($intchain->{$x}) { map{ $xc{$_}++ } split",", $intchain->{$x}; } # intchain == \%tchain
    
    if(%xc) { $tchain{$x}= join",", sort keys %xc; } # only one tchain/x
  }

  map{ $chain{$_}= $inidchain->{$_} } keys %$inidchain;

  my %cnum;
  foreach my $xc (sort values %tchain) {
    map{ $cnum{$_}++ } split",", $xc;
  }
  
  # this looks right now; 1 major problem is that 1 bad evidence join drags all others with it.
  # almost right? need to add some lesser chains to main?
  #   if($v{$c} or $v{$chain2mainchain{$c}}) ??
  # .. still missing Intron joiners b/n 2 exons; all would be at tail of cmax
  # .. because 2 larger rnas split, and intron chains reset to those 2
  # .. chain2mainchain problem > need small chains to link 2 big ones
  
  
  # all possible chains ordered from big to small
  # for complete: cmax0 = all/most exons; for split: cmax0,cmax1 split exons, remainder may join
  # for equal cnum, prefer: ovpro > ovrna > intr ?
  my @cmax= sort{$cnum{$b} <=> $cnum{$a} or $a cmp $b} keys %cnum;
  
  foreach my $x (@xi) { #2ndary chaining after 1st pass
    my $xc= $tchain{$x} or next;
    
    my @xc= split",", $xc;
    my %xc= map{$_,1} @xc;
    # my %xmc= map{ my $m=$chain2mainchain{$_}||$_; $m => 1; } @xc;
    
    # all chains held by this exon; should merge all?
    my ($cm,@cother)= grep { $xc{$_} } @cmax; # $xmc{$_} or 
    
if(0) { # not now?  use only most complete simple chains, after chainIntrons()  
    foreach my $c (@cother) { 
      ## add test here of quality: dont join subchain to main unless passes qual test
      if($xchain{$c}) { $xchain{$cm}.= delete $xchain{$c};} # move to cm ?
      $chain2mainchain{$c}= $cm; #? reset now? # unless($chain2mainchain{$c});
    } 
}

    if($cm) {  
      $xchain{$cm}.="$x,"; $ochain{$x}.="$cm,"; 
      my $ca=$cm; foreach my $c (@cother) { $ca.="-$c" if($cnum{$c} == $cnum{$cm}); }#bad for cm
      $chainequal{$cm} .= "$ca,"; # note: only equal at this exon
    }  
  }


  # sub classify( \%chain, \%xchain, \%ochain)
  # ------------------------------
  
  my @c= (sort  keys %xchain); # {$a <=> $b}
  return (0, $nx) unless(@c>0);
  
  my %nchain=(); 
  map{ my $n= $xchain{$_} =~ tr/,/,/; $nchain{$_}=$n; } @c;

  my ($cmax,@cmore)= sort{$nchain{$b} <=> $nchain{$a}} @c;
  my $cx1= $xchain{$cmax};
  my $nx1= $nchain{$cmax};
  my $cxall= join",", @xi; $cxall.=",";
  
  
  #?? add return of $chainequal{$cmax} ?? in chainids ?
  
  our %chain2id=(); 
  # foreach my $id (sort keys %chain) { my $c=$chain{$id}; $c= $chain2mainchain{$c}||$c; $chain2id{$c}.="$id,"; }
  foreach my $id (sort keys %chain) { my $c=$chain{$id}; $chain2id{$c}.="$id,"; }
  sub chainids { my $cd= join"; == ", map{ $chain2id{$_}||$_ } @_;  return $cd; }

  #sub indexchain { my @ic=(); foreach (@_){ push @ic, join(",", map{ $ichain{$_}||$_ } split",",$_);} return @ic;}
  sub indexchain { return @_; }

  my ($icx1)= indexchain($cx1);
  # complete
  return (1, $nx, $icx1, chainids($cmax)) if($nx1 >= $nx or $nx == 1 
    or ($nx > 5 and $nx1 == $nx-1 and (! $ochain{$xi[0]} or ! $ochain{$xi[-1]}) ));
       # allow missing end

  return (-3, $nx, $icx1, chainids($cmax)) 
    if($nx < 3 or $nx1 == 1 or $nx1 < 0.2*$nx); # poor chain evidence; 2 exons dont qualify for join/split
       
  # handle cases of 1 terminal exon not in chain.
  # count if/until sum( $nchain{ccc} ) == $nx; ie look for all exons, may be 3+ chains

# if($InQual =~ /complete/ and J2) change to 3? 4? -4? -3?
#     3 => 'C1T1', # complete but for terminal exon
#     4 => 'IJ4',  # join split by intron only ; EP3 weak evd
#     -1 => 'OV1', # overlapped chain
#     -2 => 'ER2', # other error
#     -3 => 'EP3', # poor/weak evidence
#     -4 => 'EM4', # endmiss, ok 


  if(@cmore) {
    my %cx1= map{$_,1}split",",$cx1;
    my @cxx=(); 
    my ($ov,$nxx, $nintron, $nxxmorethanintron)= (0,0,0,0);
    foreach my $cm (@cmore) {
      my $cx2= $xchain{$cm};
      my $nx2= $nchain{$cm};
      my @cx2= split",",$cx2; 
      $ov=0; map{ $ov++ if($cx1{$_}); $cx1{$_}++; } @cx2;
      push @cxx, $cx2; 
      $nxx+=$nx2;
      
      my $cid=$chain2id{$cm} || $cm; # BUG in chain2id ??
      $cid=~s/\b(N\d+)\b//g; 
      $nxxmorethanintron += $nx2 if($cid=~/\w/);
      ## my @nin= ($chain2id{$cm} =~ /\b(N\d+)\b/g); # wrong way
      ## $nintron += scalar(@nin);
      
      last if($ov); # return (-1, $nx, $cx1, @cxx, chainids($cmax, @cmore)) if($ov); # confused, overlapped chains
      last if($nxx+$nx1 >= $nx);
      # last if($nx2 < 2); # not enough chain evid left
      # last if($nxx+$nx1 >= $nx-1);
    }
    
    my @icxx= indexchain(@cxx);
    # my $cmoreids= chainids(@cmore);
    # my $moreintron = ($cmoreids =~ /N\d+,/g);
    
    return (-1, $nx, $icx1, @icxx, chainids($cmax, @cmore)) if($ov); # confused, overlapped chains
      # also require nx2 > min size? $nx2 < 0.2*$nx
    return (-3, $nx, $icx1, @icxx, chainids($cmax, @cmore)) if($nx>3 and  ($nx1+$nxx)/$nx < 0.75 ); # poor chain evidence
      # 3/C1T1: complete but for terminal 1 BUT check this isnt inner exon:  and index($cxall,$cx1) >= 0
    return (3, $nx, $icx1, @icxx, chainids($cmax, @cmore)) 
         if($nxx < 2 and $nx1 > 1 and $nx > 3  and index($cxall,$cx1) >= 0);
    return (4, $nx, $icx1, @icxx, chainids($cmax, @cmore))
          if($nxx > 1 and $nxxmorethanintron <= 1); # $nxx - $nintron < 1 :NOT # recode intron only splits
    return (-4, $nx, $icx1, @icxx, chainids($cmax, @cmore))
          if($nxx > 1 and $InQual =~ /complete/); 
    return (2, $nx, $icx1, @icxx, chainids($cmax, @cmore))
          if($nxx > 1); 

  } else {
    return (1, $nx, $icx1, chainids($cmax))  
        if($nx > 3 and $nx1 == $nx-1 and (! $ochain{$xi[0]} or ! $ochain{$xi[-1]} ));
    return (-3, $nx, $icx1, chainids($cmax)) 
        if($nx > 3 and $nx1/$nx < 0.75);# poor chain evidence
    # is this error?:  AUGepi4p2s9g59t1 ER-2    8       1,3,4,5,6,7,8, << missing 2 near end
    # .. yes, model has false exon-2, ovpro, ovrna and introns all reject
  }
  
  my @icmore= indexchain(@xchain{@cmore});
  # check cx1 for ends-only missing, recode that vs innermiss : EI-4, EE-5, ER-2
  # endmiss:
  return (-4, $nx, $icx1, @icmore, chainids($cmax, @cmore) )
      if($nx > 3 and $nx1 > 0.74*$nx and index($cxall,$cx1) >= 0);
  
  return (-2, $nx, $icx1, @icmore, chainids($cmax, @cmore) ); # chain error of other kind
}


# from  evigene/scripts/intronscore.pl
sub chainIntrons
{
  # my($g, $flag, $nx, $valin, $insum, $insplit, @exons)= @_;
  my @xi= sort{ $a<=>$b } keys %gx; # change from gx{xi} to gx{xb}; use values for xi ?
  my $nx= @xi;

  my @join=(0) x $nx; 
  my %tchain=(); 
  my %idchain=(); 
  my %lid=(); my %id; 
  my($good,$bad,$join,$ginval,$inscore,$ix,$ic) = (0) x 10;

  foreach my $x (@xi) {
    $ix++;  
    %lid=%id; %id=(); 
    my $rt= $gr{$intag}{$x} or next; 
    my ($ival,@ids)= split",",$rt;

    $ival =~ s,/.*,,; 
    $good++ if($ival > 0); 
    $bad++  if($ival < 0); 
    
    my $j=0;  my %xc=();
    foreach my $id (@ids) { 
      my $c= $lid{$id};  
      # each simple,good,chained,inner exon has 2 Inids = intron splice ids
      #   condense 2 chain ids to 1 continuous over intron chain. ie reset id chain #
      # if($c) { $j= $ix; } else { $c= ++$ic; }
      if($c) { $j= $ix; } else { $c= $id; $c=~s/^N/si/; } #  or $j=$c unless($j)?
      $id{$id}= $c;  $idchain{$id}= $c;
      $xc{$c}++;
      } 

    if(%xc) { $tchain{$x}= join",", sort  keys %xc; } # {$a <=> $b}; only one tchain/x
    if($j) { $join[$ix-1]=$j;  $join++; }
  }

  # now for introns (only) join 2ndary chains 
  my %cnum=();
  foreach my $xc (sort values %tchain) {  map{ $cnum{$_}++ } split",", $xc; }
  my @cmax= sort{$cnum{$b} <=> $cnum{$a} or $a cmp $b} keys %cnum;
  my (%xchain, %ochain, %chain2mainchain);
  
  foreach my $x (@xi) { #2ndary chaining after 1st pass
    my $xc= $tchain{$x} or next;
    my @xc= split",", $xc;
    my %xc= map{$_,1} @xc;
    my %xmc= map{ my $m=$chain2mainchain{$_}||$_; $m => 1; } @xc;

    # all chains held by this exon 
    my ($cm,@cother)= grep { $xmc{$_} or $xc{$_} } @cmax;
    # merge to mainchain for introns (only) 
    foreach my $c (@cother) { 
      if($xchain{$c}) { $xchain{$cm}.= delete $xchain{$c};} # move to cm ?
      $chain2mainchain{$c}= $cm; #? reset now? # unless($chain2mainchain{$c});
    } 
    if($cm) { $xchain{$cm}.="$x,"; $ochain{$x}.="$cm,"; } #? reset ochain or not
  }
  
  # reset id>chain
  map { my $c=$idchain{$_}; $idchain{$_}= $chain2mainchain{$c} || $c; } keys %idchain;

  # not used# my $bjoin = join",",@join; 
  
  my $notjoin = ($nx - 1) - $join; $notjoin=0 if($notjoin < 0);
  $inscore += $join - $notjoin; #? 

  if($good==0 and $bad==0) { $ginval="none"; }
  # if($insum < 0) { $ginval="err.rev"; } # <0 includes retained intron (same strand)
  elsif($bad > 0) { $ginval="poor"; }
  elsif( $nx > 2 and $good == $nx and $join + 1 == $nx) { $ginval="complete"; } 
  elsif( $nx > 3 and $good >= $nx-2 and ($join >= $nx-3 or $join > 0.75 * $nx) ) { 
    my($j1,$j2)=(2,$nx-2); if($join > 8) { $j1++; $j2--; } 
    if( grep { $_ == 0 } @join[$j1..$j2]) { $ginval="miss_inner";} 
    else { $ginval="complete_inner"; } 
    } 
  elsif($good and $join) { $ginval= ($good/$nx > 0.66) ? "good" : "ok"; }  
  else { $ginval="poor"; } 

  return( \%idchain, \%ochain, $inscore, $ginval);  
}

__END__

=item cases

  fixme: this is join (maybe 3 part)
AUGepit5_xfemp3s5g48t1  EP-3    22      2,3,4,5,6,7,8,9,        13,14,15,16,17,18,
AUGepi4as42g133t4       EP-3    25      1,2,3,4,5,6,7,8,9,10,11,12,13,14,       15,16,17,19,
AUGepi6cp3s4g87t2       EP-3    39      20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,  5,6,7,8,9,10,11,16,

is this error or join?
AUGepi6bp4s1g120t1      EP-3    11      1,2,3,4,5,6,  [missing7,8]  9,10,

is this error or complete?
AUGepi4p1s24g96t1       ER-2    10      1,2,3,4,5,8,9,10,

this is error/weak evid, not nec. join: too many 1,2-exon parts
-- set min #parts for 2nd largest as pct total? or min exons?

AUGepit5_xpupas34g58t1  J2      15      10,11,12,13,14, 3,4,    6,7,    8,9,    2,      5,      15,
AUGepit5_xadults52g11t1 J2      8       1,2,3,4,        5,6,    7,8,

  not enough evidence:
    AUGepit3p1s9g153t1      J2      17      6,7,8,  2,3,    4,5,    9,10,   11,12,  13,14,  1,
                                  ^^^^^^^^^^^^^^^^ too little if nx1 < 0.n * nx
                                  
    m6AUGepit2p13c4g152t1   J2      5       3,4,    2,      5,

    AUGaf_S9g4037t1 J2      2       1,      2,
    AUGepi4ap2s9g219t1      J2      2       1,      2,
    m6AUGytwingmp14c4g163t1 J2      2       1,      2,
    
    AUGepit5_xembryop2s9g139t1      J2      3       2,      3,
    AUGepi6cp2s9g139t1      J2      3       1,      2,
  
  err or ok? middle 1 exon out of chain == overlap?
  AUGepi6cp1s9g16t1       C1T1    8       1,2,3,4,6,7,8,  5,

  AUGepit5_xwingp2s9g23t1 J2      17      4,5,6,7,8,9,10,11,12,13,17,     14,15,16,
  
  r8nvit1v4aSallLoc49491t6        J2      12      3,4,5,6,7,8,9,10,       1,2,11,12,
    ^^ this is complete, not sure why ends split from middle << annot mistake.

    # .. both these are not-joins, 1 has more evid + overlapped asmrna, 2nd is weak evid.
    # AUGepit3p3s9g5t1        J2      16      1,2,3,4,5,6,7,8,9,10,11,12,13,  14,15,
    # AUGepi6cp3s9g4t1        J2      20      4,5,6,7,8,      11,12,13,14,

=item sample

head equal/evg11e.overjoin2.notfull2.uids | sed 's/^/=/' | ggrep -F -f - bestgenes_of10.11e.gff | $evig
ene/scripts/joincheck.pl 
AUGaf_S1g539t1          0
AUGaf_S10g4193t1        2       1,2,3,4,5,6,7,8,        10,11,12,13,
AUGaf_S118g12925t1      2       5,6,7,8,9,10,11,        1,2,3,
AUGaf_S14g5183t1        2       1,2,3,4,5,6,7,8,9,10,   11,12,13,14,
AUGaf_S17g6204t1        2       3,4,5,6,7,      8,9,10,
AUGaf_S2g879t1          2       6,7,8,9,10,11,12,13,    2,3,4,
AUGaf_S28g8005t1        2       5,6,7,8,9,10,   1,2,3,4,
AUGaf_S28g8009t1        2       1,      2,
AUGaf_S28g8046t1        -1      3,4,5,6,7,      2,3,
AUGaf_S36g8997t1        2       1,2,3,4,5,6,    7,8,9,


cat equal/evg11e.overjoin2.notfull2.uids | sed 's/^/=/' | grep AUGepit5_xmalep5s1g2t1 | ggrep -F -f - b
estgenes_of10.11e.gff | $evigene/scripts/joincheck.pl | less
AUGepit5_xmalep5s1g2t1  2       1,2,3,4,5,6,7,8,9,10,   11,12,13,


SCAFFOLD1       AUGepit5_xmale  mRNA    7933246 7943723 1288,77,0,67,67,1005,2774,22,84,0,1301,0,892,2418       +       .       ID=AUGepit5_xmalep5s1g2t1;aalen=805,72%;homolog=1288/1667,apis2ref:XP_392140.3;inqual=84;
  nintron=22/26;ovpro=67,E2BVB4_9HYME/67.00,apis2gno_150384/66.00;ovrna=67,r8nvit1v3S2big0Loc2674t2/67.64,r8nvit1cuf83c_Gsc1g8664t2/33.31;cxlen=2418/3310,73%;inexon=11/14/13;scoresum=21437

>> chain1 by rseq,intr
SCAFFOLD1       AUGepit5_xmale  exon    7933246 7933666 0,0,0,0,0,0,421,0,0,0,16,0,0,0  +       .       Parent=AUGepit5_xmalep5s1g2t1;
	rseq=421/421,b5nvit1v3S2big0Loc2674t2;intr=16,N3181,N3182;ref=31/421,NV10441-RA
SCAFFOLD1       AUGepit5_xmale  exon    7934100 7934192 0,0,0,0,0,93,93,0,0,0,84,0,0,0  +       .       Parent=AUGepit5_xmalep5s1g2t1;
	rseq=93/93,b5nvit1v3S2big0Loc2674t2;intr=84,N3182,N3183;ref=36/93,NV10441-RA
SCAFFOLD1       AUGepit5_xmale  exon    7934266 7934552 0,0,0,0,0,287,287,0,0,0,146,0,0,0       +       .       Parent=AUGepit5_xmalep5s1g2t1;
	rseq=287/287,b5nvit1v3S2big0Loc2674t2;intr=146,N3183,N3184;ref=287/287,NV10441-RA
SCAFFOLD1       AUGepit5_xmale  exon    7934636 7934767 0,0,0,0,0,132,132,0,0,0,147,0,0,0       +       .       Parent=AUGepit5_xmalep5s1g2t1;
	rseq=132/132,b5nvit1v3S2big0Loc2674t2;intr=147,N3184,N3185,N3186;ref=132/132,NV10441-RA
SCAFFOLD1       AUGepit5_xmale  exon    7935685 7935852 0,0,0,0,0,168,168,0,0,0,188,0,0,0       +       .       Parent=AUGepit5_xmalep5s1g2t1;
	rseq=168/168,b5nvit1v3S2big0Loc2674t2;intr=188,N3186,N3187,N3188;ref=168/168,NV10441-RA
SCAFFOLD1       AUGepit5_xmale  exon    7935932 7936003 0,0,0,0,0,72,72,0,0,0,176,0,0,0 +       .       Parent=AUGepit5_xmalep5s1g2t1;
	rseq=72/72,b5nvit1v3S2big0Loc2674t2;intr=176,N3188,N3189;ref=72/72,NV10441-RA
SCAFFOLD1       AUGepit5_xmale  exon    7936085 7936337 0,0,0,0,0,253,253,0,0,0,138,0,0,0       +       .       Parent=AUGepit5_xmalep5s1g2t1;
	rseq=253/253,b5nvit1v3S2big0Loc2674t2;intr=138,N3189,N3190,N3191;ref=253/253,NV10441-RA
SCAFFOLD1       AUGepit5_xmale  exon    7936611 7936787 0,0,0,0,0,0,177,0,0,0,82,0,0,0  +       .       Parent=AUGepit5_xmalep5s1g2t1;
	rseq=177/177,b5nvit1v3S2big0Loc2674t2;intr=82,N3191,N3192;ref=177/177,NV10441-RA
SCAFFOLD1       AUGepit5_xmale  exon    7936872 7937215 0,0,0,0,0,0,344,0,0,0,143,0,0,0 +       .       Parent=AUGepit5_xmalep5s1g2t1;
	rseq=344/344,b5nvit1v3S2big0Loc2674t2;intr=143,N3192,N3193;ref=344/344,NV10441-RA
SCAFFOLD1       AUGepit5_xmale  exon    7937294 7937461 0,0,0,0,0,0,168,0,0,0,92,0,0,0  +       .       Parent=AUGepit5_xmalep5s1g2t1;
	rseq=168/168,b5nvit1v3S2big0Loc2674t2;intr=92,N3193;ref=168/168,NV10441-RA
<< chain1end
-- split --
>> chain2 by rseq,intr
SCAFFOLD1       AUGepit5_xmale  exon    7941626 7941763 0,0,0,0,0,0,138,0,0,0,28,0,0,0  +       .       Parent=AUGepit5_xmalep5s1g2t1;
	rseq=138/138,b5nvit1cuf83c_Gsc1g8664t2;intr=28,N3196;ref=138/138,NV10442-RA
SCAFFOLD1       AUGepit5_xmale  exon    7942003 7942410 0,0,0,0,0,0,408,0,0,0,19,0,0,0  +       .       Parent=AUGepit5_xmalep5s1g2t1;
	rseq=408/408,b5nvit1cuf83c_Gsc1g8664t2;intr=19,N3197;ref=408/408,NV10442-RA
SCAFFOLD1       AUGepit5_xmale  exon    7942611 7942723 0,0,0,0,0,0,113,0,0,0,30,0,0,0  +       .       Parent=AUGepit5_xmalep5s1g2t1;
	rseq=113/113,b5nvit1cuf83c_Gsc1g8664t2;intr=30,N3198,N3199;ref=113/113,NV10442-RA
SCAFFOLD1       AUGepit5_xmale  exon    7943188 7943723 0,0,0,0,0,0,0,0,0,0,12,0,0,0    +       .       Parent=AUGepit5_xmalep5s1g2t1;
  intr=12,N3199;ref=133/536,NV10442-RA


=cut