#!/usr/bin/env perl
# overlapfilter.perl

=item usage

  perl overlapfilter -act keep|drop|mark -mark=terepeat 
      -overlaps terepeats.gff -itype=gff|blast -input stdin|tandy.gff > tandyfilt.gff

  generalized to filter or mark GFF for any GFF overlap data (protein & EST matches, repeats)
  
  gzcat *exons_tandy6j.gff.gz | perl $td/overlapfilter.perl \
    -act mark -mark=proteinhsp -over *prot9-hsp.gff.gz -in stdin > dere_exons_tandy6jm.gff

  # ... cut retained introns from rna-assemblies :
  overlapfilter -in nvit1_rnaseq.cuff1rs13.gff.gz -over intron_good.gff.gz \
  -pass 'exon,intron' -strand -pct 100 -typeover cut -act mark -mark inov > nvit1_rnaseq.cuff1rs13.incut.gff

=cut

use strict;
use warnings;
use Getopt::Long;

use constant SAMEBASE => 3; # or 0; # for _sameloc, slop allowed in loca == locb
use constant SPLICE   => 3; # bp for intron splice span tests
use constant SPLICEX  => SPLICE - 1; # for exon
use constant { ACT_DROP=>1, ACT_KEEP=>2, ACT_MARK=>3, ACT_MARK_WITH_ID=>4,
      ACT_MARK_OVERBASE=> 5, ACT_MARK_OVERBASEID=> 6};
use constant { kOVERLAP=>1, kSAMELOC=>2, kNEARLOC=>3, kCUTOVER=>4 };

use constant { kINTRON2SPLICE_OVER=>1, kINTRON2SPLICE_QUERY=>2, kINTRONERROR_OVER=>3, };

our $BINSIZE  = 1000 ; #was# 5000;
our $NEARDIST =  500; # was 15k; needs to be < BINSIZE
our $SUMBASEOVER= 0; #which? see below
our $SPANBASEOVER= 0; # 2015.04, return overlap spans, along w/ other possible items: overID, overBASEsum
our $OVERisNEAR=  1; # 0 or 1 default ? use -type nearover to set; -type nearonly to unset over

my $debug=0;
my ($overlaps,$overlaplist,$markidtype,$input,$itype,$action,$actid,$ok,$mark,$baseover);
my ($overtype,$typeover,$sametypes,$pctover,$passtypes)= (kOVERLAP,"",1,0,"");
my ($save_pctover,$save_baseover)=(0,0);
my ($save_miss,$n_overlaps)=(0,0);
my ($intron2splice,$stranded,$negstrand,$orderedoverlap,%ordered_over);

my $optok= GetOptions(
  "overlaps=s", \$overlaps, 
  "typeover=s", \$typeover, 
  "pctover=s", \$pctover, 
  "NEARDIST=i", \$NEARDIST, "BINSIZE=i", \$BINSIZE, 
  "SUMBASEOVER!", \$SUMBASEOVER,
  "SPANBASEOVER!", \$SPANBASEOVER,
  "mark=s", \$mark, 
  "input|in=s", \$input,  
  "itype=s", \$itype,  
  "midtype=s", \$markidtype, # return not ID= but other attribute or score
  "action=s", \$action, 
  "passtypes=s", \$passtypes,  
  "orderedoverlap!", \$orderedoverlap, 
  "intron2splice:s", \$intron2splice, 
  "stranded!", \$stranded, "negstrand!", \$negstrand, 
  "debug!", \$debug, 
  "baseover!", \$baseover, # show base,pct overlap counts per item and total
  );

die "
usage:  perl overlapfilter -act keep|drop|mark|markid|markbase -mark=terepeat 
        -overlaps terepeats.gff -typeover overlap|inside|sameloc|samefeat|near
        -pctover=0 (for overlap type, min % overlap)
        -baseover : show overlap base count 
        -intron2splice  : convert -over introns.gff to exon splice end points +/- 2bp
        -orderedoverlap : keep overlap ids ordered by input
        -neardist=$NEARDIST (for near type, base distance)
        -passtypes='CDS' | 'mRNA,CDS' | 'gene.Gnomon,gene.Chainer' : act on these types only
        -midtype=ID|Name|score|source|... (for markid, attribute or score to mark, ID default)
        -itype=gff|blast -input stdin|tandy.gff > tandyfilt.gff
" unless($optok and $action and (-f $overlaps or -f $input));

if($action =~ /^cut/) { $typeover= "cut"; $action= "mark"; }
## for matching IDs b/n updates; need also match type: gene/mRNA/exon/CDS/...

$mark  ||= "over";
$itype ||= "gff";
$actid= ($action =~ /keep/) ? ACT_KEEP : ($action =~ /mark/) ? ACT_MARK : ACT_DROP;
$actid= ACT_MARK_OVERBASEID if($actid == ACT_MARK && $action =~ /id/i && $action =~ /base/i);
$actid= ACT_MARK_WITH_ID if($actid == ACT_MARK && $action =~ /id/i);
$actid= ACT_MARK_OVERBASE if($actid == ACT_MARK && $action =~ /base/i);

$overtype= ($typeover =~ /^over/) ? kOVERLAP : ($typeover =~ /inside/) ? kOVERLAP :
           ($typeover =~ /same/) ? kSAMELOC : ($typeover =~ /near/) ? kNEARLOC :
           ($typeover =~ /cut/) ? kCUTOVER : kOVERLAP;
$sametypes= ($typeover =~ /feat/) ? 1 : 0; #??

$pctover= $pctover/100.0 if($pctover);
# fix dependency; UPD15, was 0.10 min proportion; want pctover= 1/infinity for any gap overlaps..
$pctover= 0.0000001 if(!$pctover and ($actid == ACT_MARK_OVERBASE or $actid == ACT_MARK_OVERBASEID));

$passtypes =~ s/[,]/\|/g;
# $SUMBASEOVER= ($actid == ACT_MARK_OVERBASE) ? 0 : 1; #?? which?

# ** need another intron2splice option: 
#  1. match intron:any part of exon (i.e. exon mistakes) vs 
#  2. match intron:exon ends (i.e. exon accuracy)
#.. for stats want both 1,2 from same comparison. can report as summary or 
#.. mark both on exons?: intr=+FwdEnds,~FwdMid,-RevEnds,~RevMid  
#..   intr=+FwdEnd.FwdMid/-RevEnd.RevMid; or add other mark? inerr= splice inside exon
if(defined $intron2splice) {
  if($intron2splice =~ /^over/) { $intron2splice=kINTRON2SPLICE_OVER; }  
  elsif($intron2splice =~ /^in/) { $intron2splice=kINTRON2SPLICE_QUERY; }
  elsif($intron2splice =~ /^err/) { $intron2splice=kINTRONERROR_OVER; } # = over exons: splice and inside  
  elsif($intron2splice =~ /\d/ and $intron2splice>1) {
    $intron2splice= ($intron2splice == 3)? kINTRONERROR_OVER : kINTRON2SPLICE_QUERY;
  } else {
    $intron2splice= kINTRON2SPLICE_OVER;
  }
} else {
  $intron2splice=0;
}


my $sumscore= ($markidtype and $markidtype =~ /sum/)?1:0; # as in scoresum, pctidsum, ...; remove sum from marktype?
my $avescore= ($markidtype and $markidtype =~ /ave/)?1:0; 
my $intronsum= ($intron2splice and $markidtype 
  and $markidtype =~ /^score|orient|strand/ and $markidtype =~ /sum/)?1:0;

if($overtype == kNEARLOC) {
$BINSIZE= int($NEARDIST*2) if($NEARDIST > $BINSIZE);
$OVERisNEAR= 1 if($overtype =~ /over/);
$OVERisNEAR= 0 if($overtype =~ /only/);
}

# allow -overlap gff == stdin if -input is a file?
##local *R;
my $ovh; ##= *OVR;
#$ok = ($overlaps =~ /.gz$/) ? open(R,"gunzip -c $overlaps |") : open(R,$overlaps);
   if($overlaps =~ /.gz$/) { $ok= open(OVR,"gunzip -c $overlaps |");  $ovh= *OVR; }
elsif($overlaps =~ /^(stdin|-)/) { $ovh= *STDIN; $ok=1; }
else { $ok= open(OVR,$overlaps); $ovh= *OVR; }
die "bad -overlaps=$overlaps" unless($ok);
$overlaplist= collect_overlaps($ovh); close($ovh);

my $inh= *STDIN;
$ok = ($input =~ /.gz$/) ? open($inh,"gunzip -c $input |") 
      : ($input =~ /^(stdin|-)/) ? $inh= *STDIN
      : open($inh,$input);
die "bad -input=$input" unless($ok);

my($sum_baseover,$sum_pctover,$sum_basetotal,$sum_pcttotal,$sum_miss)=(0) x 9;

my $nin=0;
my $nr=0;
if($itype =~ /blast/i) { $nr= filter_blast($inh); }
else { $nr= filter_gff($inh); }
warn"#overlaps found=$nr\n" if $debug;

$sum_basetotal||=1;
$sum_pcttotal||=1;
warn "# base statistics: overlaps n=$nr , input n=$nin ,  overset n=$n_overlaps\n"
  if ($baseover || $debug);
if ($baseover) {  $nr||= 1;
my $ave=sprintf("ave_baseover=%.3f, ave_pctover=%.3f, ave_miss=%.3f",
        $sum_baseover/$sum_basetotal,$sum_pctover/$sum_pcttotal, $sum_miss/$nr);
warn "# $ave 
# sum_baseover=$sum_baseover, sum_pctover=$sum_pctover, sum_miss=$sum_miss 
# sum_basetotal=$sum_basetotal, sum_pcttotal=$sum_pcttotal\n" 
} ;

#..................

sub filter_gff
{
  my($inh)= @_;
  my $nr=0;
  my $printpass=1;
    # ($actid == ACT_DROP or $actid == ACT_MARK or $actid == ACT_MARK_WITH_ID or $actid == ACT_MARK_OVERBASE)? 1 : 0;
  $printpass=0 if($actid == ACT_DROP or $actid == ACT_KEEP);

  while(<$inh>){
    unless(/^\w/){ next if(/^(#n |$)/); print and next; }
    my $line=$_;
    my($ref,$src,$typ,$tb,$te,$tp,$to,$tx,$tattr)= split"\t";

    if($passtypes and "$typ.$src" !~ m/$passtypes/) { print if $printpass; next; } # pass other types
      # ^ do in both filter_gff and in collect_overlaps **??
      
    $nin++;
    $save_baseover= $save_pctover= 0;
    
# ** need another intron2splice option: 
#  1. match intron:any part of exon (i.e. exon mistakes) vs 
#  2. match intron:exon ends (i.e. exon accuracy)
#.. for stats want both 1,2 from same comparison. can report as summary or 
#.. mark both on exons?: intr=+FwdEnds,~FwdMid,-RevEnds,~RevMid  
#..   intr=+FwdEnd.FwdMid/-RevEnd.RevMid
## $intron2splice == kINTRONERROR_OVER << new
    #** should do this for both intron2splice cases, dont count
    #   intron splices matching middle of exons **
    #   but for exons here (kINTRON2SPLICE_OVER), change offset to tb..tb+1, te-1..te

    my( $tb2, $te2)=(0,0);
    if($intron2splice == kINTRON2SPLICE_QUERY) {
      ($tb,$te, $tb2, $te2)= 
        ($to eq "-") ? ($te+1,$te + SPLICE,$tb - SPLICE,$tb-1) : ($tb - SPLICE,$tb-1,$te+1,$te + SPLICE); #intron 3bp + 1shift
    } elsif($intron2splice == kINTRON2SPLICE_OVER) {
      ($tb,$te, $tb2, $te2)= 
        ($to eq "-") ? ($te - SPLICEX,$te,$tb,$tb + SPLICEX) : ($tb,$tb + SPLICEX,$te - SPLICEX,$te); # 3bp of exon end
    }
    
    my $isover= ($overtype == kOVERLAP) ? overlaps($ref,$tb,$te,$to)
              : ($overtype == kSAMELOC) ? sameloc($ref,$tb,$te,$to,$typ) 
              : ($overtype == kNEARLOC) ? nearloc($ref,$tb,$te) 
              : ($overtype == kCUTOVER) ? cut_overlaps($ref,$tb,$te,$to) 
              : overlaps($ref,$tb,$te,$to); # added to 2010.10

    if( ($intron2splice == kINTRON2SPLICE_QUERY 
      or $intron2splice == kINTRON2SPLICE_OVER)
      and $te2) {
      
      my $over2= ($overtype == kOVERLAP) ? overlaps($ref,$tb2,$te2,$to)
          : ($overtype == kSAMELOC) ? sameloc($ref,$tb2,$te2,$to,$typ) 
          : ($overtype == kNEARLOC) ? nearloc($ref,$tb2,$te2) 
          : ($overtype == kCUTOVER) ? cut_overlaps($ref,$tb2,$te2,$to) 
          : overlaps($ref,$tb2,$te2,$to); # added to 2010.10
      if($over2 and $over2 ne $isover) { $isover = ($isover) ? "$isover,$over2" : $over2; }
    }
    
    
    if($actid >= ACT_MARK and m/$mark=/) { s/;$mark=[^;\s]+//; } ## ** remove any old mark, always **
    if($isover) { # now $isover == ID of match
      $nr++;
      
      if($overtype == kCUTOVER) {  # $isover == \@nonoverlap locs
        my $ci=1; my $cn= ($tattr=~/\S/) ? ";ci" : "ci";
        foreach my $be (@$isover) {
          my($nb,$ne)= @$be;
          my $ns= 1+$ne-$nb;
          (my $nline= $line) =~ s/\t$tb\t$te/\t$nb\t$ne/;
          $nline =~ s/$/$cn=$nr.$ci;clen=$ns/; $ci++;
          print $nline;
        }
      next;
      }
      
      next if($actid == ACT_DROP);
      if($actid == ACT_MARK) { s/$/;$mark=1/; } 
      elsif($actid == ACT_MARK_WITH_ID) { s/$/;$mark=$isover/ ;  } 
      elsif($actid == ACT_MARK_OVERBASE) { my $w=1+$te-$tb; s,$,;$mark=$save_baseover/$w,;  }
      elsif($actid == ACT_MARK_OVERBASEID) { my $w=1+$te-$tb; s|$|;$mark=$save_baseover/$w,$isover|;  }
      
      $sum_baseover += $save_baseover;
      $sum_pctover  += $save_pctover;
      $sum_miss     += $save_miss;
    } else { 
      next if($actid == ACT_KEEP); 
    }
    print;
  }
  return $nr;
}




sub filter_blast # ncbi format=8,9 blast table
{
  my($inh)= @_;
  my $nr=0;
  while(<$inh>){
    unless(/^\w/){ print and next; }
    my($qid,$ref,$pid,$align,$xa,$xb,$qb,$qe,$tb,$te,@bmore)= split"\t";
    # fix blast loc swap for orient
    my $to='+';  ($tb,$te,$to)= ($te,$tb,'-') if($tb>$te);
    my $typ="HSP";
    
    $nin++;
    $save_baseover= $save_pctover= 0;
    my $isover= ($overtype == kSAMELOC) ? sameloc($ref,$tb,$te,$to,$typ) #($ref,$tb,$te,$to,$typ)
              : ($overtype == kNEARLOC) ? nearloc($ref,$tb,$te) 
              : overlaps($ref,$tb,$te,$to); # added to 2010.10
    if($isover) {
      $nr++;
      next if($actid == ACT_DROP);
      if($actid == ACT_MARK) {  s/$/\t$mark=1/; } 
      elsif($actid == ACT_MARK_WITH_ID) {  s/$/\t$mark=$isover/;} 
      elsif($actid == ACT_MARK_OVERBASE) { s/$/\t$mark=$save_baseover/ ;  } 
      $sum_baseover += $save_baseover;
      $sum_pctover  += $save_pctover;
      $sum_miss     += $save_miss;
    } else { 
      next if($actid == ACT_KEEP); 
    }
    print;
  }
  return $nr;
}

my $warns=0;

# convert to array handling? _max(@list) 
sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }

sub overlaps
{
  my($ref,$tb,$te,$to)= @_;    
  my (@lid, @oids);
  $save_baseover= $save_pctover= 0;
  my $twidth= 1 + $te - $tb;
  $sum_basetotal += $twidth;
  $sum_pcttotal  += 1;
  $save_miss = 0;
  my($missb, $misse)=(undef,undef);
  my %didid=();
  my @overs; # improved sumoverlap  
  
  return 0 unless($overlaplist->{$ref});

  my($tb1,$te1, $tb2, $te2)=(0) x 4;
  if($intron2splice == kINTRONERROR_OVER) {
    # input == exon, over= intron, test if overlap is splice end or internal
    ($tb1,$te1, $tb2, $te2)= 
      ($to eq "-") ? ($te - SPLICEX,$te,$tb,$tb + SPLICEX) : ($tb,$tb + SPLICEX,$te - SPLICEX,$te);  
  }


  my ($ib1, $ib2)= (int($tb/$BINSIZE) , int($te/$BINSIZE));
  for (my $ib = $ib1; $ib <= $ib2; $ib++) 
  {
    $overlaplist->{$ref}{$ib} or next;
    my @locs= @{$overlaplist->{$ref}{$ib}};
    foreach my $rloc (@locs) {
      my ($lb,$le,$lid,$lo,$oid)= @{$rloc}[0,1,3,4,6];

      next if($didid{$oid.$lb.$le}++);
      
      my $over= ($tb <= $le && $te >= $lb) ? 1 : 0;
      
      if($over and $intron2splice == kINTRONERROR_OVER) {
        my $ok=0; # note: lb,le here are one splice end span of intron: 3 bp?
        my $samestrand= ($to =~ /[+-]/ and $lo =~ /[+-]/ and $to ne $lo) ? 0 : 1;
        
        if($tb1 <= $le && $te1 >= $lb) { $ok= ($samestrand) ? 1 : -1; } # $errt="rev" unless $samestrand
        elsif($tb2 <= $le && $te2 >= $lb) { $ok= ($samestrand) ? 1 : -1; }
        elsif(($tb + 2*SPLICE <= $lb && $te - 2*SPLICE >= $le) and $samestrand) 
          { $ok= -1; } # $errt="retain" ; same strand inside; change val: -3 ? == retained intron err
        # else what?
        
        unless($ok == 0) {
        my $val= $ok;
        if($intronsum) { $val= $ok * abs($lid); } # lid == intron strand * score
        ## new version for intronsum: keep val, add list @lid == intron oids: 
        ##   intr=+nn/-mm,oid1,oid2,..
        push @oids, $oid; # just for intronsum? 
          # FIXME: these oids count intron ENDs (2 per intron) .. fix here or client use N11 == N12
        push @lid, $val; # lid == intron score always?? or not
        $save_baseover += $val; 
        $save_pctover  += ($val < 0) ? -1 : 1; 
        }
        next;
      } 
      
      $over=0 if($stranded and ($to =~ /[+-]/ and $lo =~ /[+-]/ and $to ne $lo)); # NEW 2010.10; fixme for to or lo eq "."
      $over=0 if($negstrand and ($to !~ /[+-]/ or $lo !~ /[+-]/ or $to eq $lo)); # NEW 2010.10
      if($over and $typeover =~ /inside/) { # overlap is all inside feat
        $over= ($tb <= $lb && $te >= $le) ? 1 : 0;
      }
      
      if($over and $pctover) {
        # my $maxo= _max( $le - $tb, $te - $lb); # wrong
        my ($bb,$be)= ( _max($tb,$lb), _min($te,$le) );
        my $maxo= abs(1+$be - $bb);

        my $leno= _min( abs($le - $lb), abs($te - $tb)) || 1;
        my $pover= $maxo/$leno;
        $over = 0 if $pover < $pctover;

        #?? add save_basemiss_in, save_basemiss_over        
        # $save_basemiss = _max($save_basemiss, $bb - $tb, $bb - $lb, $te - $be, $le - $be);
        if($over) {
          $missb= (defined $missb) ? _min($missb, abs($tb - $lb)) : abs($tb - $lb);
          $misse= (defined $misse) ? _min($misse, abs($te - $le)) : abs($te - $le);
          if($avescore) { $lid = $maxo * $lid; } # expand by nbases over; see below?
          
          # which?
          unless($SUMBASEOVER or $SPANBASEOVER) {
            $save_baseover = $maxo  if($maxo > $save_baseover); 
            $save_pctover  = $pover if($pover > $save_pctover); 
            #? push( @overs, [$bb,$be]);
          } else { # which?
            # $save_baseover += $maxo; # FIX, really want all non-duplic overlaps
            # $save_pctover  += $pover; # FIX
            push( @overs, [$bb,$be]);
            }
          }
        
        }
      push @lid, $lid if($over);
      # ^^ collect *all* overlap ids as list to return
      }
    }
  $missb ||= 0; $misse ||= 0;
  $save_miss = $missb + $misse;

  ## 2015.04.20 : add opt to return @overs spans matching this, as per return join",",@lid?
  my @overo=(); # add1504
  if(($SUMBASEOVER or $SPANBASEOVER) and @overs) {
    my($bb,$be)= ($tb,$te);
    my $baseo= 0;
    @overs= sort _sort_over @overs;
    foreach my $ab (@overs) {
      my($lb,$le)= @$ab;
      if($le < $bb) { next; }
      elsif($lb > $te) { last; } 
      elsif($lb < $bb) {
        my $e= ($be < $le) ? $be : $le; 
        $baseo += 1 + $e - $bb; 
        push @overo, "$bb-$e" if($SPANBASEOVER);
        $bb= $e+1;
        }
      elsif($lb <= $be) {  
        my $e= ($be < $le) ? $be : $le; 
        $baseo += 1 + $e - $lb; 
        push @overo, "$lb-$e" if($SPANBASEOVER);
        $bb= $e+1; 
        } 
      last if($bb > $te);
    }
    $save_baseover= $baseo;
    $save_pctover = ($twidth > 0) ? $baseo / $twidth : 0;
    if($SPANBASEOVER and @overo) {
      my $ret = 'ovspan:'.join(",",@overo); #? will this work w/ conflict to other returns?
      push @lid, $ret unless($intronsum or $sumscore or $avescore);
    }
  }
  
  if(@lid) {
    if($intronsum) {
       # sum of +score,-score 
        ## new version for intronsum: keep val, add list @lid == intron oids: 
        ##   intr=+nn/-ee,oid1,oid2,..
       my($pos,$neg, $ret)=(0,0,0); 
       foreach my $v (@lid) { if($v<0) { $neg += -$v; } elsif($v>0) { $pos += $v; } }
       if($pos>0 and $neg>0) { $ret= ($neg>$pos)? "-$neg/+$pos" : "+$pos/-$neg"; }
       elsif($neg>0) { $ret= -$neg; }
       else { $ret= $pos; }
       # FIXME: 2013sep: order oids by val (@lid) : hifreq first.
       if(@oids>2) {
          my %ord=(); for my $i (0..$#oids) { $ord{$oids[$i]}= $lid[$i]; }; 
          @oids= sort{ $ord{$b} <=> $ord{$a} or $a cmp $b } @oids;
       }
       $ret .= "," . join(",",@oids) if(@oids); # always or not?
       return $ret;
      }
    elsif($sumscore or $avescore) { 
      my $sum=0; map{$sum+=$_}@lid; 
      if($avescore) { $sum= sprintf "%.3g", $sum/($te - $tb + 1); } # sprintf "%.3g" ?
      return $sum; 
    }
    
    my %lid= map{$_,1}@lid; 
    if($orderedoverlap) {
      @lid= sort{ $ordered_over{$a} <=> $ordered_over{$b} } keys %lid;
    } else {
      @lid= sort keys %lid;
    }
    
    return join",", @lid; 
    }
  return 0;
}


sub _sort_over { # @[b,e] min-b, max-e
  return ($a->[0] <=> $b->[0]) || ($b->[1] <=> $a->[1]);
}

sub cut_overlaps
{
  my($ref,$tb,$te,$to)= @_;    
  my (@lid,%didid,@overs);
  $save_baseover= $save_pctover= 0;
#   $sum_basetotal += $te - $tb + 1;
#   $sum_pcttotal  += 1;
    
  return 0 unless($overlaplist->{$ref});
  #my @bins= (int($tb/$BINSIZE) .. int($te/$BINSIZE));
  #foreach my $ib (@bins)
  my ($ib1, $ib2)= (int($tb/$BINSIZE) , int($te/$BINSIZE));
  for (my $ib = $ib1; $ib <= $ib2; $ib++) 
  {
    $overlaplist->{$ref}{$ib} or next;
    my @locs= @{$overlaplist->{$ref}{$ib}};
    foreach my $rloc (@locs) {
      my ($lb,$le,$lid,$lo,$oid)= @{$rloc}[0,1,3,4,6];
      next if($didid{$oid.$lb.$le}++);
      my $over= ($tb <= $le && $te >= $lb) ? 1 : 0;      
      $over=0 if($stranded and ($to =~ /[+-]/ and $lo =~ /[+-]/ and $to ne $lo)); # NEW 2010.10; fixme for to or lo eq "."

      # inside fix: eg, cut only introns internal to exons
      if($over and ($typeover =~ /inside/ or $pctover >= 0.99)) {   
        $over= ($tb <= $lb && $te >= $le) ? 1 : 0;  
        if($over) { $save_baseover += 1+$le-$lb; $save_pctover+=1; }
      }
      elsif($over and $pctover) {
        my ($bb,$be)= ( _max($tb,$lb), _min($te,$le) );
        my $leno= _min( abs($le - $lb), abs($te - $tb)) || 1;
        my $maxo= abs(1+$be - $bb);
        my $pover= $maxo/$leno;
        $over = 0 if $pover < $pctover;
        $save_baseover = $maxo  if($maxo > $save_baseover); # FIX, really want all non-duplic overlaps
        $save_pctover  = $pover if($pover > $save_pctover); # FIX
      }
        
      push @overs, [$lb,$le] if ($over);
      }
    }
    
  my @opens=();
  if(@overs) {
    my($bb,$be)= ($tb,$te);
    @overs= sort _sort_over @overs;
    foreach my $ab (@overs) {
      my($lb,$le)= @$ab;
      if($le < $bb) { next; }
      elsif($lb <= $bb && $le > $bb) { $bb= $le+1; }
      elsif($lb < $be) {  #  && $le > $be
        my($b1,$e1)= ($bb,$lb-1);
        push @opens, [$b1,$e1];
        $bb= $le+1; 
        } 
      elsif($lb > $te) { last; } #?
      last if($bb >= $te);
    }
    # end point
    if($bb < $te) { push @opens, [$bb,$te]; }
    return \@opens;
  } else {
    return 0; ## [[$tb,$te]];
  }
}


sub nearloc
{
  my($ref,$tb,$te)= @_;    
  return 0 unless($overlaplist->{$ref});
  my %didid=();
  my @lid;

  my $tm= int(($tb+$te)/2); # change to min-distance from ends
  # my @bins= (int($tb/$BINSIZE) .. int($te/$BINSIZE));  foreach my $ib (@bins) 
  my ($ib1, $ib2)= (int($tb/$BINSIZE) , int($te/$BINSIZE));
  for (my $ib = $ib1; $ib <= $ib2; $ib++) 
  {
    $overlaplist->{$ref}{$ib} or next;
    my @locs= @{$overlaplist->{$ref}{$ib}};
    foreach my $rloc (@locs) {
      my ($lb,$le,$lid,$lo,$oid)= @{$rloc}[0,1,3,4,6];
      next if($didid{$oid.$lb.$le}++);
      
      my $over= ($tb <= $le && $te >= $lb) ? 1 : 0;
      # $over=0 if($stranded and ($to =~ /[+-]/ and $lo =~ /[+-]/ and $to ne $lo)); # NEW 2010.10; fixme for to or lo eq "."

#       if($over and $pctover) {
#         my ($bb,$be)= ( _max($tb,$lb), _min($te,$le) );
#         my $maxo= abs($be - $bb);
#         my $leno= _min( abs($le - $lb), abs($te - $tb)) || 1;
#         $over = 0 if $maxo/$leno < $pctover;
#         }

      if($over) { # skip overlaps ? need option: nearAndOver vs nearNotOver, or pipe overlapfilter 2 ways      
        push @lid, $lid if ($OVERisNEAR); 
      } else {
        ## my $bd= abs($tb - $le); my $ed= abs($lb - $te); my $mind= ($ed < $bd) ? $ed : $bd; 
        my $mind= _min( abs($tb - $le), abs($lb - $te));
        push @lid, $lid if ($mind < $NEARDIST); 
      }
      }
    }
    
  if(@lid) {
    if($sumscore or $avescore) { 
      my $sum=0; map{$sum+=$_}@lid; 
      if($avescore) { $sum= sprintf "%.3g", $sum/($te - $tb + 1); }
      return $sum; 
      }
    my %lid= map{$_,1}@lid; 
    if($orderedoverlap) {
      @lid= sort{ $ordered_over{$a} <=> $ordered_over{$b} } keys %lid;
    } else {
      @lid= sort keys %lid;
    }
    return join",", @lid; 
    
  }
  return 0;
}

sub sameloc
{
  my($ref,$tb,$te,$to,$typ)= @_;   
  my @lid;
  return 0 unless($overlaplist->{$ref});
#  warn join",",("sameloc",$ref,$tb,$te,$to,$typ),"\n" if $debug and $warns++<10;
  $to ||= '+';
  ($tb,$te,$to)= ($te,$tb,'-') if($tb>$te);
  my @bins= (int($tb/$BINSIZE) .. int($te/$BINSIZE));
  foreach my $ib (@bins) {
    $overlaplist->{$ref}{$ib} or next;
    my @locs= @{$overlaplist->{$ref}{$ib}};
    foreach my $rloc (@locs) {
      my ($lb,$le,$lid,$lo,$ltyp)= @{$rloc}[0,1,3,4,5];
      # [$tb,$te,$ref,$gid,$to,$typ]
      # also need match ft types?, orient,
      next if($sametypes and not($typ eq $ltyp and $to eq $lo)); #?
      next if($stranded and ($to =~ /[+-]/ and $lo =~ /[+-]/ and $to ne $lo)); # NEW 2010.10; 
      next if($negstrand and ($to !~ /[+-]/ or $lo !~ /[+-]/ or $to eq $lo)); # NEW 2010.10

      push @lid, $lid if(abs($tb-$lb) <= SAMEBASE and abs($te-$le) <= SAMEBASE);
      }
    }
  if(@lid) { 
    if($sumscore or $avescore) { my $sum=0; map{$sum+=$_}@lid; return $sum; }  
    my %lid= map{$_,1}@lid; return join",", sort keys %lid; }
  return 0;
}



sub collect_overlaps
{
  my($gff)= @_;
  my %overlaps=(); my $nr=0;
  while(<$gff>){
    next unless(/^\w/); chomp;
    my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr)= split"\t"; #,@gffmore
    $tattr ||="";  
    
    if($passtypes and "$typ.$src" !~ m/$passtypes/) { next; } # pass other types
      # ^ do in both filter_gff and in collect_overlaps **??
    
    $nr++;
    my $oid= "N".$nr; # save always for uniq feature test
    my($gid); # fixme: separate $gid from markid : want two fields, one for id tests
    if($markidtype) {
    
      if($intron2splice and $markidtype =~ /^score|orient|strand/i) {
         $gid = "$to$tp"; # strand.score = +nnn,-nnn,.nnn,0nnn possible results
      }
      elsif($markidtype =~ /^score/i) { $gid=$tp; }
      elsif($markidtype =~ /^source/i) { $gid=$src; }
      elsif($markidtype =~ /^type/i) { $gid=$typ; }
      elsif($markidtype =~ /^orient|strand/i) { 
          my $ov=($to eq "-")?"-1":($to eq "+")?"+1":($to eq ".")?"0":$to;
          $gid=$ov; } # change to +1,-1,0
      elsif($tattr =~ m/\b$markidtype=([^;]+)/) { $gid=$1; }
    } elsif($tattr) {
      if($tattr =~ m/\bID=([^;]+)/) {  $gid=$1; }
      elsif($tattr =~ m/\bParent=([^;]+)/) {  $gid=$1; }
      elsif($tattr =~ m/\bTarget=([^;\s]+)/) { $gid=$1; } # NEW 2010.10
    }
    unless(defined $gid) { $gid = $oid; }

    $ordered_over{$gid}= $nr; # if $orderedoverlap;

    if($intron2splice == kINTRON2SPLICE_OVER or $intron2splice == kINTRONERROR_OVER) { # 2010jul
      my($s1b,$s1e,$s2b,$s2e)= 
        ($to eq "-") ? ($te+1,$te + SPLICE,$tb - SPLICE,$tb-1) : ($tb - SPLICE,$tb-1,$te+1,$te + SPLICE); # 3bp + 1 shift
      ($tb,$te)= ($s1b,$s1e);
      my $rloc= [$tb,$te,$ref,$gid,$to,$typ,$oid]; 
      my @bins= (int($tb/$BINSIZE) .. int($te/$BINSIZE));
      foreach my $ib (@bins) { push( @{$overlaps{$ref}{$ib}}, $rloc); }
      ($tb,$te)= ($s2b,$s2e);  #? change gid, oid?
      #NOT NOW# $nr++; $oid= "N".$nr;  # dang, should oid be same for same intron, or +1 ?
      
    } elsif($intron2splice == kINTRON2SPLICE_QUERY) {
      # here is exon over, need to trim to exon endpoints to test overlap at ends?
      my($s1b,$s1e,$s2b,$s2e)= 
        ($to eq "-") ? ($te - SPLICEX,$te,$tb,$tb + SPLICEX) : ($tb,$tb + SPLICEX,$te - SPLICEX,$te); # 3 bases of exon
      ($tb,$te)= ($s1b,$s1e);
      my $rloc= [$tb,$te,$ref,$gid,$to,$typ,$oid]; 
      my @bins= (int($tb/$BINSIZE) .. int($te/$BINSIZE));
      foreach my $ib (@bins) { push( @{$overlaps{$ref}{$ib}}, $rloc); }
      ($tb,$te)= ($s2b,$s2e);  #? change gid, oid?
      #NOT NOW# $nr++; $oid= "N".$nr; # dang, should oid be same for same intron, or +1 ?
    }
      
    my $rloc= [$tb,$te,$ref,$gid,$to,$typ,$oid]; # change to string; save mem
    # my @bins= (int($tb/$BINSIZE) .. int($te/$BINSIZE));
    my ($ib1, $ib2)= (int($tb/$BINSIZE) , int($te/$BINSIZE));
    for(my $ib=$ib1; $ib<=$ib2; $ib++) { push( @{$overlaps{$ref}{$ib}}, $rloc); }
    }
  $n_overlaps=$nr; # global
  warn"#collect_overlaps=$nr\n" if $debug;
  return \%overlaps;
}


# sub _isoverlap {
#   my($gb,$ge, $qb,$qe)= @_; 
#   return ($gb <= $qe && $ge >= $qb) ? 1 : 0;
# }

