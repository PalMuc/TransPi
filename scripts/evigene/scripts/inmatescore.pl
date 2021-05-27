#!/usr/bin/env perl
# inmatescore.pl

=item usage

  inmatescore  -over validexons.gff -in genes.gff > genescore.table
 
=item update2:
  
  test this with other -over gene_evidence.gff
    where gene_evid has Parent/Target IDs linking exonic evidence,  e.g. 
    -- aligned protein (exonerate) CDS exons.
    -- read_intron-joined exons
    -- read_mated exons

  gzgrep exonr.CDS prot/protein.gff.gz |\
  $evigene/scripts/inmatescore.pl -exontype CDS -over stdin -in $geneset \
  > $geneset.inprot
   
  use as replacement for gene-level sensitivity (genome-total) stats  
  as well as quality scores per gene model.

=item fixme some stats maybe flaky

  summary stats need checking for accuracy
    
=item update : matedgenescore.pl ? inmatescore.pl ?
  
    * revise to handle also mate-pair scoring, using mated-exons.gff
      with Parent=mated-exons-id
      Update to ditto with introns? creating inpaired-exons.gff (? need exon evidence to fill introns?)
      with Parent=introned-exons-id
      
    * try new scoring: 
          true-pos  = pairexon bases matched
          false-neg = pairexon bases missed by predictor exons, bp and miss-bp/pair-exons
          false-pos = pairexon bp exceeded by pred exons, bp and extra-bp/pred-exons
                  
=cut

use strict;
use warnings;
use Getopt::Long;

my $debug=1;
our $BINSIZE  = 500 ; #was# 5000;

use constant SPLICE   => 3; # bp for intron splice span tests
use constant SPLICEX  => SPLICE - 1; # for exon
use constant{  jBEGIN => 0, jEND => 1, jCHR => 2, jGID => 3, jSTRAND => 4,
               jTYPE  => 5, jOID => 6, jATTR => 7, jSCORE => 8, jLINE => 9 }; ## gff overlap record
              
use constant KB => 1000; # or 1024?
use constant MB => KB * KB;
use constant GB => MB * MB;

my $CHECK_OVEROVER=1; # option?
my ($input, $overlaps, $passtypes) = ("") x 10;
## my ($gid, $insplit, $valin, $insum, $nex, $header, @exons);
my ($ok, $domateloc, $didhead, $stranded, $insorted, $norun, $pctover, $showhitid, $skipnohit, $n_overlaps)= (0) x 20;
my (%sumscore, %genesplits, %overhitids);
my $exontype="exon"; # allow CDS

$stranded= 1;  # defaults
$pctover= 50;
$skipnohit= 1;

my $optok= GetOptions(
  "overlaps=s", \$overlaps,
  "input=s", \$input,
  #drop#"passtypes=s", \$passtypes,  
  "exontype=s", \$exontype,  
  "stranded!", \$stranded,  #? should be default
  "pctover=i", \$pctover,   #  default? 
  "skipnohit!", \$skipnohit,  #? should be default
  "showhitid!", \$showhitid,
  #"verbose|v!", \$verbose, 
  #"norun|n", \$norun, 
  "debug", \$debug, 
  );

die "usage: inmatescore  -over validexons.gff -in genes.gff > genescore.table
  opts: -stranded $stranded -pctover $pctover -noskipnohit -exontype $exontype [exon|CDS]  .. \n"
  unless($optok and $input);

#drop# $passtypes =~ s/[,]/\|/g;
$pctover= $pctover/100.0 if($pctover);

# NOT HERE: my $nintr = add_intronexon_scores( $genes, $overlaps); # exon intr=, inids=, 
# DONT NEED NOW: my $nsplit= add_insplit($genes, $overlaps); # %insplit{geneid}

my $ovh;  
   if($overlaps =~ /.gz$/) { $ok= open(OVR,"gunzip -c $overlaps |");  $ovh= *OVR; }
elsif($overlaps =~ /^(stdin|-)/) { $ovh= *STDIN; $ok=1; }
else { $ok= open(OVR,$overlaps); $ovh= *OVR; }
die "bad -overlaps=$overlaps" unless($ok);

my($noverlist, $overlaplist, $featlist)  = collect_overlaps($ovh); close($ovh);

my $geneh= *STDIN;
if($input) {
$ok = ($input =~ /.gz$/) ? open($geneh,"gunzip -c $input |") 
      : ($input =~ /^(stdin|-)/) ? $geneh= *STDIN
      : open($geneh,$input);
die "bad -input=$input" unless($ok);
}

my( $nin, $nover )= score_genes($geneh);

# output_gene_scores();

sub score_genes
{
  my($geneh)= @_;
  my (%genehits, %genehitbase);
  my ($gid, $nex, $nexhit, $exwidth, $genequal) = (0) x 10;
  
  while(<$geneh>) {
    unless(/^\w/) { next; } #push(@savegene, $_);  # print if outformat = gene.an.gff
    my($ref,$src,$typ,$tb,$te,$tp,$to,$tx,$tattr);
    ($ref,$src,$typ,$tb,$te,$tp,$to,$tx,$tattr)= split"\t";
   
    if( $typ eq "mRNA" ) { 
      
      $genequal= geneScore( $gid, $nex, $exwidth, \%genehits, \%genehitbase,) if($gid);
      
      ($gid)=  m/ID=([^;\s]+)/ or die "ERROR: missing mRNA ID=\n"; 
      
      $exwidth= $nexhit= $nex= 0; ## $valin= $insum=0; @exons=();
      %genehitbase= %genehits=(); #? only need 1 at a time?
      
    } elsif( $typ eq $exontype ) {  #? passtypes ==  exontype?  
   
      my($pg)= m/Parent=([^;\s]+)/; 
      if($pg ne $gid) { warn "ERROR: exon Parent=$pg in gene=$gid out of order\n"; next; } ## assume input gff is gene-ordered  
      $nex++; 
      my $twidth= 1 + $te - $tb;
      $exwidth += $twidth;
      
      my ($isover, $bpover) =  overlaps( $ref,$tb,$te,$to); # this marks featlist items
    
      if($isover) {  
        $nover++; $nexhit++;
        foreach my $overid (@$isover) {
          # result is list of overgene ids;
          # record which overgenes this gene matches, and n-bases (== exon width)
          $genehitbase{$overid} += $bpover; ## $twidth;  #** Problem, want overlap bp here, not this model exon width
          $genehits{$overid} ++;
        }
      } else {
        $genehitbase{'MISS'} += $twidth;
        $genehits{'MISS'} ++;
      }

    }
    
  }
  
  $genequal= geneScore( $gid, $nex, $exwidth, \%genehits, \%genehitbase, "LAST") if($gid);
}


#....................
# convert to array handling? _max(@list) 
sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }


sub geneScore
{
  my($geneid, $nexon, $bpexon, $hits, $hitbase, $flag)= @_;
  $nexon ||= 1;
  
  my @second=();
  my($noverids, $maxhit, $maxbphit, $maxid, $overn, $overbp)= (0) x 10;
  my $nmiss = delete $hits->{'MISS'}; # n-exon
  my $bpmiss= delete $hitbase->{'MISS'}; # n-bases

  unless( not $skipnohit or scalar(keys %{$hits}) or defined($flag)) { 
    my %score = ( 
      modelid => $geneid,
      n_nohit => 1,
      nexon_nohit => $nexon, bpexon_nohit => $bpexon, 
      );
    my @sums= qw( n_nohit nexon_nohit bpexon_nohit );
    sumScores( $geneid, \%score, \@sums, {} );

    my $genequal = "0/none";
    return $genequal; # or print? see opt
  }

  foreach my $hitid (sort keys %{$hits}) {
    my $nhit  = delete $hits->{$hitid};  
    my $bphit = delete $hitbase->{$hitid}; 
    
    if($nhit > 0) { # always >0 ?
    $noverids++; # check for joins w/ 2+
    $overhitids{$hitid}++;
    
    if($nhit > $maxhit) #? or use $bphit > $maxbphit or both ?
    # if($bphit > $maxbphit)  
      { 
      if($maxhit > 0 and $maxbphit > 99) {
        @second= ($maxid, $maxhit, $maxbphit, $overn, $overbp); # keep 2nd-max stats
      }
      my ($noverex, $noverbp, $noverscore)= getOvergene( $hitid); #? check for splits
      $maxid= $hitid;
      $maxhit= $nhit;
      $maxbphit= $bphit; # keep same? if($bphit > $maxbphit);
      $overn= $noverex;
      $overbp= $noverbp;
      #? record noverscore ? a quality value
      }
    }
  }

  push( @{$genesplits{$maxid}}, $geneid);
  
  sub npscore { my($ntrue, $nhit)=@_; 
    if($ntrue <= $nhit or $ntrue < 1) { return(0,0); } else { 
    my $d= $ntrue - $nhit; 
    return ($d, sprintf("%.3f",$d/$ntrue)); 
    }
  }
  sub truescore { my($ntrue, $nhit)=@_; 
    if($ntrue < $nhit or $ntrue < 1) { return 1; } else { 
    return sprintf("%.3f",$nhit/$ntrue); 
    }
  }
  
  ## use $nmiss, $bpmiss here? same as from nexon, bpexon subtraction
  ## maybe count false pos only where have overgene, not for all genes (for summary?)
  
  my ($nfneg, $pnfneg)= npscore( $overn, $maxhit);
  my ($nfpos, $pnfpos)= npscore( $nexon, $maxhit);
  my ($ntpos, $pntpos)= ($maxhit, truescore( _max($nexon,$overn), $maxhit) );
  
  my ($bfneg, $pbfneg)= npscore( $overbp, $maxbphit); #? bad calc here?
    # ?? some double counting in summary of same gene $overbp from split genes
  my ($bfpos, $pbfpos)= npscore( $bpexon, $maxbphit);
  
    ## bpexon can reasonably be > than maxbphit, false-pos not nec problem esp. for prot homol.
  my ($btpos, $pbtpos)= ($maxbphit, truescore( _max($bpexon,$overbp), $maxbphit) );
  
  my($issplit, $nsplit, $bsplit)=(0,0,0); # from where? need OverGenes report
  
  my($isjoin, $njoin, $bjoin)=(0,0,0); ## isjoin == maybejoin
  if(@second) {
    # @second=($maxid, $maxhit, $maxbphit, $overn, $overbp)
    if( $second[1] > 1 and $second[2] > 99) { ($isjoin, $njoin, $bjoin)= (1,@second[1,2]); }
    # @second stats ? is it a join problem?
    ## fneg says if it is a split problem
    
  }
      # should these vals be OPTIONS?
      # NOTE this assumes overlap set is unique loci, otherwise noverids >> 1 often
  my $perfect = ($pntpos > 0.94 and $pbtpos > 0.94 and $noverids == 1) ? 1 : 0;
  my $mostly  = ($pntpos > 0.69 and $pbtpos > 0.65 and $noverids == 1) ? 1 : 0;
      #? is this right; seems we are missing perf cases
      
  ## gene quality score to add to predictor mRNA:  bpscore/flags 
  # bpscore = $maxbphit * $pbtpos ?  and/or  -falsneg -falspos
  # bTrue == btpos == maxbphit may be best single score
  
  my $genequal= int($maxbphit * $pbtpos);
  $genequal .= '/' . (($perfect) ? "perfect" : ($mostly) ? "mostly" : "partial");
  # ^^ problem output format here? , but want this format for gene.annot.gff; ?  num/qual to num\tqual?
  
  my %score = ( 
    modelid => $geneid, 
    nexon => $nexon, bpexon => $bpexon, 
    perfect => $perfect, mostly => $mostly,
    overgene  => $noverids, overbp => $overbp,
    ## exonshit => $maxhit, baseshit => $maxbphit,  # == ntrue, btrue
    hitid => $maxid, # optional print this
    
    ntrue => $ntpos, pntrue => $pntpos,
    btrue => $btpos, pbtrue => $pbtpos,

    nfneg => $nfneg, pnfneg => $pnfneg,
    bfneg => $bfneg, pbfneg => $pbfneg,

    nfpos => $nfpos, pnfpos => $pnfpos,
    bfpos => $bfpos, pbfpos => $pbfpos,
    
    isjoin => $isjoin, bjoin => $bjoin,
    issplit => $issplit, # bsplit => $bsplit,
    );
    
  ## print not return $score;  
  ## ?? change column order: all N ... all BP ... pct
  my @shows= qw( nexon bpexon overgene ntrue btrue nfneg bfneg nfpos bfpos  pbfneg pbfpos);
  my @showl= qw( nExon bpExon nOvGene nTrue bTrue nFNeg bFNeg nFPos bFPos  pFNeg pFPos);
  if($showhitid) { push(@shows, "hitid"); push(@showl, "HitID"); }
  
  do{ printf "%-25s\t%-12s\t","ModelID","Quality"; print join("\t", @showl),"\n"; } unless($didhead++);
  printf "%-25s\t%-12s\t",$geneid,$genequal; print join("\t", @score{@shows}), "\n";
  
  ## maybe count false pos only where have overgene, not for all genes (for summary?)
  ## BUT then should record nexon, bpexon from nohit models:
  # ... nNohitGene, nNohitExon, bpNohitExon
  # ... OR option to ignore all no-hit models (but for count)?
  
  my @sums= qw( nexon bpexon perfect mostly overgene overbp ntrue btrue nfneg bfneg nfpos bfpos isjoin issplit);
  ## my %hitsums= ( nfpos=>1, bfpos=>1 );
  sumScores( $geneid, \%score, \@sums, {}); ## \%hitsums);
  
  ## do sums over all genes: total at end
  if(defined($flag) and $flag =~ /LAST/) {
 
    # want totals from OverGenes for this
    #    my ($novergene, $noverex, $noverbp, $noverscore)= sumOvergene(); 
    foreach my $overid (sort keys %genesplits) {
      my @gids= @{$genesplits{$overid}};
      if(@gids > 1) { $sumscore{issplit}++; }
    }
    
    %score = %sumscore;
    $score{"n_nohit"} ||= 0;
    my $nhitgene= $score{"n"} - $score{"n_nohit"};
    
    ## redo percents from totals
    ## want also count of perfect genes
    ( $score{ 'pperfect' }) = truescore( $nhitgene, $score{ 'perfect' });
    ( $score{ 'pmostly' }) = truescore( $nhitgene, $score{ 'mostly' });
    
    ## one of these is bad calc (or both): bfneg is too high? or overbp too low?
    ## some summary have  btrue  > overbp: sum($maxbphit) > sum($overbp) shouldnt be possible?
    #  .. is possible because maxbphit == model exons, vs overbp == over exons
    ## and sumOvergenes gives different/lower counts, due to double over hit counting above?
    
    (undef, $score{ 'pbfneg' }) = npscore( $score{ 'overbp' }, $score{ 'btrue' }); #?? sometimes 0 bad
    ## ( $score{ 'pbfneg' }) = truescore( $score{ 'overbp' }, $score{ 'bfneg' });  
    
    ##(undef, $score{ 'pbfpos' }) = npscore( $score{ 'bpexon' }, $score{ 'btrue' }); #* use bfpos here
    ( $score{ 'pbfpos' }) = truescore( $score{ 'bpexon' }, $score{ 'bfpos' });  

    ## convert bases to Mb for summary
    foreach my $bs (qw(bpexon overbp btrue bfneg bfpos bpexon_nohit)) { $score{$bs} = prbase($score{$bs}); }
    
    $score{"n"} = $nhitgene;
    
#     ## ?? change column order: all N ... all BP ... pct
#     my @shows= qw(     n nexon ntrue nfneg nfpos  bpexon overbp btrue bfneg bfpos  pbfneg pbfpos  perfect pperfect mostly isjoin issplit);
#     my @showl= qw( nGene nExon nTrue nFNeg nFPos  bpExon bpOver bTrue bFNeg bFPos  pFNeg  pFPos   pPerf nPerf Mostly Joins Splits);

    my @shows= qw(     n nexon bpexon overbp ntrue btrue nfneg bfneg nfpos bfpos pbfneg pbfpos  perfect pperfect mostly isjoin issplit);
    my @showl= qw( nGene nExon bpExon bpOver nTrue bTrue nFNeg bFNeg nFPos bFPos  pFNeg pFPos   nPerf pPerf Mostly Joins Splits);
    print"\n";

    printf "%-10s\t","# TOTALS"; print join("\t", @showl),"\n"; 
    printf "%-10s\t","# ModelHit"; print join("\t", @score{@shows}), "\n";

    if($score{"n_nohit"}) {
    @shows= qw( n_nohit nexon_nohit bpexon_nohit  );
    printf "%-10s\t","# ModelMiss"; print join("\t", @score{@shows}), "\n";
    }
    
    ## summary of OvergeneMiss == EvidMiss?  : n_nohit, nexon_nohit, bpexon_nohit 
    ## use this for summary bFNeg? or like above, keep separate line
    
    my( $nmissov, $nexon_missov, $bpexon_missov, $evscore_missov,
        $nhitov,  $nexon_hitov,  $bpexon_hitov, $evscore_hitov)=
        sumOvergenes( \%overhitids );    
    printf "%-10s\t","# OverHit"; print join("\t", $nhitov, $nexon_hitov, prbase($bpexon_hitov) ), "\n";
    printf "%-10s\t","# OverMiss"; print join("\t", $nmissov, $nexon_missov, prbase($bpexon_missov) ), "\n";

    print"\n";
  }

  return $genequal;
}


sub sumScores {
  my($geneid, $score, $show, $hitsums)= @_;
  $sumscore{'n'}++;  ## ngene
  ## my $nohit= ($score->{overgene} > 0) ? 0 : 1;
  foreach my $sc (@$show) {
    ## next if($nohit and $hitsums->{$sc});
    $sumscore{$sc} += $score->{$sc};  # median, min, max ?
  }
}

sub sumOvergenes
{
  my( $hitidhash )= @_;
  # sum of both hit, nohit here?
  my( $nmissov, $nexon_missov, $bpexon_missov, $evscore_missov,
      $nhitov,  $nexon_hitov,  $bpexon_hitov, $evscore_hitov)= (0) x 10;
  foreach  my $gid (sort keys %{$featlist} ) {
    ## next if($hitidhash->{$gid});
    my ($noverex, $noverbp, $noverscore)= getOvergene( $gid);
    if($hitidhash->{$gid}) {
    $nhitov++;
    $nexon_hitov += $noverex;
    $bpexon_hitov += $noverbp;
    $evscore_hitov += $noverscore;
    } else {
    $nmissov++;
    $nexon_missov += $noverex;
    $bpexon_missov += $noverbp;
    $evscore_missov += $noverscore;
    }
  }
  return ( $nmissov, $nexon_missov, $bpexon_missov, $evscore_missov,
          $nhitov,  $nexon_hitov,  $bpexon_hitov, $evscore_hitov);
}

sub getOvergene
{
  my($overid) = @_;
  my($noverex, $noverbp, $noverscore)= (0,0, 0);
  if( my $gexons= $featlist->{$overid} ) {
    foreach my $rloc (@$gexons) {
      my ($lb,$le,$lid,$lo,$oid,$lval)= @{$rloc}[ jBEGIN,jEND,jGID,jSTRAND,jOID,jSCORE ]; # exon span
      my $xwidth= 1 + $le - $lb;
      $noverex++;
      $noverbp += $xwidth;
      $noverscore += $lval;
    }
  }
  
  return ($noverex, $noverbp, $noverscore);
}

sub prbase {
  my $nb= shift;
  $nb ||=0;
  my @lv=(GB, MB, KB); my @lb=("Gb","Mb","Kb");
  foreach my $i (0..$#lv) { 
    my $bv= $nb / $lv[$i];
    if($bv > 1) {
      my $dig= ($bv > 9.5) ? 0 : 1; # was $DIGITS
      return sprintf "%.${dig}f".$lb[$i], $bv; 
      }
    }
  my $dig= ($nb >= 2) ? 1 : 2; 
  return sprintf "%.${dig}f",$nb;
}


sub _sort_over { # @[b,e] min-b, max-e
  return ($a->[0] <=> $b->[0]) || ($b->[1] <=> $a->[1]);
}

sub overlaps
{
  my( $ref,$tb,$te,$to)= @_;    # == predictor exon here

  my ($nover, $maxbpover)= (0,0); 
  my @oid= (); my @overs=();
  my $twidth= 1 + $te - $tb;
  my %didid=();
  
  return 0 unless($overlaplist->{$ref});

  my ($ib1, $ib2)= (int($tb/$BINSIZE) , int($te/$BINSIZE));
  for (my $ib = $ib1; $ib <= $ib2; $ib++) 
  {
    $overlaplist->{$ref}{$ib} or next;
    my @locs= @{$overlaplist->{$ref}{$ib}};
    
    foreach my $rloc (@locs) {
      my ($lb,$le,$lid,$lo,$oid)= @{$rloc}[ jBEGIN,jEND,jGID,jSTRAND,jOID ]; # exon span
      next if($didid{$oid.$lb.$le}++);
      
      my $over= ($tb <= $le && $te >= $lb) ? 1 : 0;
      $over=0 if($stranded and ($to =~ /[+-]/ and $lo =~ /[+-]/ and $to ne $lo));  

      if($over and $pctover) {
        my ($bb,$be)= ( _max($tb,$lb), _min($te,$le) );
        my $maxo= abs(1 + $be - $bb); ## or is this overwidth: 1+ $le - $lb
        
        my $leno= _min( abs($le - $lb), abs($te - $tb)) || 1;
        my $pover= $maxo/$leno;
        $over = 0 if $pover < $pctover;
        
        #no, use all overs# $maxbpover= $maxo if($maxo > $maxbpover); # or sum non-overlapped sections ??
        push( @overs, [$bb,$be]);
        }

      push @oid, $lid if($over); # ? oid or lid == geneid
      }
  }

  if(@overs) {
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
        $bb= $e+1;
        }
      elsif($lb <= $be) {  
        my $e= ($be < $le) ? $be : $le; 
        $baseo += 1 + $e - $lb; 
        $bb= $e+1; 
        } 
      last if($bb > $te);
    }
    $maxbpover= $baseo;
  }
  
  if(@oid) { 
    my %oid= map{ $_,1 } @oid; 
    @oid= sort keys %oid;     
    return (\@oid, $maxbpover); ## ((@oid > 1) ? \@oid : 0); # only care about joins of 2+ oid
    }
  return 0;
}

sub overlap1
{
  my( $overlaplist, $ref,$tb,$te,$to)= @_;    # == predictor exon here

  my @oid= (); my @overs=();
  my $twidth= 1 + $te - $tb;
  my %didid=();
  
  return 0 unless(ref $overlaplist and $overlaplist->{$ref});

  my ($ib1, $ib2)= (int($tb/$BINSIZE) , int($te/$BINSIZE));
  for (my $ib = $ib1; $ib <= $ib2; $ib++) 
  {
    $overlaplist->{$ref}{$ib} or next;
    my @locs= @{$overlaplist->{$ref}{$ib}};
    
    foreach my $rloc (@locs) {
      my ($lb,$le,$lid,$lo,$oid)= @{$rloc}[ jBEGIN,jEND,jGID,jSTRAND,jOID ]; # exon span
      next if($didid{$oid.$lb.$le}++);
      
      my $over= ($tb <= $le && $te >= $lb) ? 1 : 0;
      $over=0 if($stranded and ($to =~ /[+-]/ and $lo =~ /[+-]/ and $to ne $lo));  

      push @oid, $lid if($over); # ? oid or lid == geneid
      }
  }

  if(@oid) { 
    my %oid= map{ $_,1 } @oid; 
    @oid= sort keys %oid;     
    return (\@oid);  
    }
  return 0;
}

sub collect_overlaps
{
  my($ingff)= @_;  
  my( %overlaps, %rlocs, %feats); # returns
  ## these are exons to join
  my %types;
  
  my ($nr,$noverover)=(0,0);
  while(<$ingff>) {   
    next unless(/^\w/); 
    my $inline= $_;
    chomp;
    my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr)= split"\t";  
    $tattr ||="";  
    #drop# if($passtypes and "$typ.$src" !~ m/$passtypes/) { next; } # pass other types
    
    $nr++;
    my $oid= "N".$nr; # save always for uniq feature test
    my($gid);  
    if($tattr) {
      if($tattr =~ m/\bParent=([^;]+)/) {  $gid=$1; }  #* expect this for exons
      elsif($tattr =~ m/\bTarget=([^;\s]+)/) { $gid=$1; } # NEW 2010.10
      elsif($tattr =~ m/\bID=([^;]+)/) {  $gid=$1; } ## avoid this one
    }
    unless(defined $gid) { $gid = $oid; }

    $types{$typ}++;
    # limit overlaps to exontype ??
    
    # ?? handle alt-tr exons here? compress to one exon per loc, append all parent gid, use same oid
    # -- skip exon ID=, use only gid=Parent or Target

    # ** FIXME do: check overlaps for self-overlap; exclude, otherwise no class perfect/mostly
    if($CHECK_OVEROVER) {
      my ($isover) =  overlap1( \%overlaps, $ref,$tb,$te,$to);  
      if ($isover) { $noverover++; next; }  # isover == \@ids, save that info?
    }
    
    # fixme: need numeric score; also problem: score=0,0,0
    $tp=1 unless($tp =~ /\d/ and $tp =~ /^[\d+-]/);
    $tp =~ s/[^\de\.+-]+//g; # or .. $tp =~ s/,.*//;
    
    # use constant{ jBEGIN => 0, jEND => 1, jCHR => 2, jGID => 3, jSTRAND => 4,
    #               jTYPE => 5, jOID => 6, jATTR => 7, jSCORE => 8, jLINE => 9 }; ## gff overlap record
              
    my $rloc= [$tb,$te,$ref,$gid,$to,$typ,$oid, $tattr, $tp, $inline, 0, 0]; # change to string; save mem

    # $rlocs{$oid} = $rloc; # unique hash
    push( @{$feats{$gid}},$rloc);  # want this or not?

    my ($ib1, $ib2)= (int($tb/$BINSIZE) , int($te/$BINSIZE));
    for(my $ib=$ib1; $ib<=$ib2; $ib++) { push( @{$overlaps{$ref}{$ib}}, $rloc); }
    }
    
  $n_overlaps=$nr; # global
  
  my @types= sort keys %types;
  unless( grep { $exontype eq $_ } @types ) {
    warn "# WARNING: No -exontype=$exontype in overlap types= @types\n"; # eg over=CDS, defaultin=exon
  }
  warn"#collect_overlaps=$nr, skip dup_loci=$noverover\n" if $debug;
  return ($nr, \%overlaps, \%feats); ## \%rlocs
}


__END__


=item old details

  assumes annotations in geneset: inexon, insplit, intr, inids
  gene annotation methods:  
  
  1a. intr = per exon, count introns with errors  
  1b. inexon = per mRNA, sum of intr / exons  
      inexon= now is added to mRNA by overbestgenes, from intr= exons scores
      .. dont need inexon here, count below from eoxn intr=
      
    overlapfilter -intron2splice=error -pass 'exon,intron' -act markid  -midtype scoresum \
      -mark intr -over intron.gff -in $pred
    1a result: exon annot  intr=10 ; intr=-9 ; intr=+10/-9 ...
      where intr score is from intron.gff score (count), i.e. can be large for 1 intron
      
    1b result: mRNA annot  inexon=validin/nexon, sum of exon intr=1, may be -validin
      validin here is # exons w/ +valid (or -invalid) intron sccore, not sum of intr scores
    
  2. insplit = flag mRNA with introns split across genes, 
      now added to mRNA by sub in annotate_prediction.pl

      overlapfilter -strand -intron2splice=input -pass 'exon,intron' -act markid \
        -mark insplit -over $pred -in intron.gff
     2 result: mRNA annot insplit=inid1,inid2,...
    >> dont need now, use exon joined/notjoined info for same;
       ie. if exon has inid=N22 at no other gene exon  == split error, unless
       exon has alternate inid=N23 that does join other exon.
     
  3. inids = intron splice ids per exon
     -- can this be merged w/ intr= scoring above ?
     
      overlapfilter -strand -intron2splice=over -pass 'exon,intron' -act markid \
        -mark inids  -over intron.gff -in $pred   
     3 result: exon annot inids=inid1,inid2 .. use to score joined exons from shared valid introns
       .. could append to exon intr=n,id1,id2 score
     
  4. merge all to give mRNA inqual= score, per below table, 
     qualitative scoring from best to worst:
  
      10. complete, complete_inner: introns good and join all/most exons (3+ exons)
      6. miss_inner: introns good and mostly joined, but misses 1 or more inner exon joins (could be 2 genes)
      4. good, ok  : introns match many/some exons, but not much joining 
      2. poor      : 1+ introns but errors on 1 or more exons, or no exon joins
      0. none      : no introns
      -2. err.rev (wrong intron strand hits exon end)
      -1. err.split: intron joins this and other gene

  problems applying to gene model selection
    - pred1 w/ 10 exons, 5 joined vs pred2 w/ only same 5 joined exons
      pred1 could be better, from other evid. 
      inqual should score both same here.
      
    - pred1 w/ 10 exons, 4 introns scattered among, could be good or bad model
      vs pred2a, pred2b w/ same 4 introns (2 each), joining fully 3 exons each
      pred2a,b probably better, but pred1 could be better from other evid. 
      inqual should score both same or not?

  inqual numeric score? = +validsplices -invalidsplices +joined exons -unjoined exons -other errors?
    - for large # exons/introns, intron-complete is better than incomplete (eg longer/more exons) model
      e.g. 10-exon complete > 15-exon incomplete at locus,
      but not so for few introns, 4-exon complete vs 10-exon incomplete
      .. maybe depends on average # exons/gene ?

=cut



=item results

   genes/bestgenes.DGILmix7p.inmate
ModelID                       Quality           nExon bpEx  nOv  nTrue bTrue nFNeg bFNeg nFPos bFPos pFNeg pFPos
mp7AUGepir9p1s1g20t1          77/partial        13    3046  1     3     486   0     0     10    2560  0     0.840
mp7AUGepi5p1s1g18t1           1549/perfect      3     1549  1     3     1549  0     0     0     0     0     0
mp7AUGepir9p1s1g16t1          1540/partial      4     2464  2     3     1948  0     461   1     516   0.191 0.209
mp7PASAgasmbl_43              281/partial       6     2269  2     3     799   0     0     3     1470  0     0.648
mp7AUGepir10p1s1g13t1         2537/mostly       7     3071  1     6     2791  0     0     1     280   0     0.091
mp7AUGepir10p1s1g8t1          3295/mostly       2     3850  1     2     3850  0     646   0     0     0.144 0
mp7AUGepir9p1s1g30t1          1151/perfect      4     1151  1     4     1151  0     0     0     0     0     0
..

mp7AUGepir10s999g18t1         1278/perfect      4     1354  1     4     1354  0     80    0     0     0.056 0
mp7AUGepir16bs999g24t1        1465/partial      8     1579  2     7     1551  0     90    1     28    0.055 0.018
mp7AUGepir9s999g17t1          2304/mostly       9     2683  1     8     2486  1     0     1     197   0     0.073
mp7AUGepir9s999g21t1          1081/perfect      5     1081  1     5     1081  0     0     0     0     0     0
mp7PASAgasmbl_218826          0/partial         2     501   0     0     0     0     0     2     501   0     1.000

# TOTALS      nGene nExon bpEx  nTrue bTrue nFNeg bFNeg nFPos bFPos pFNeg pFPos nPerf pPerf Most  Join Splits
# HitGenes  	14707	110148 37Mb	86026	29Mb	12941	6.3Mb	24122	8.1Mb	0.097	0.219	5049	0.343	8177	371	1300
# MissGenes 	21555	77595	 37Mb

## better calc here
TOTALS          nGene   nExon   bpExon  bpOver  nTrue   bTrue   nFNeg   bFNeg   nFPos   bFPos   pFNeg   pFPos   nPerf   pPerf   Mostly  Joins   Splits
ModelHit        14707   110148  37Mb    32Mb    86026   25Mb    12941   7.3Mb   24122   12Mb    0.227   0.331   2288    0.156   6943    371     1300
ModelMiss       21555   77595   37Mb
OverHit         15058   90397   29Mb
OverMiss        2964    6023    2.7Mb


  genes/bestgenes.DGILmix7q.gff

mq7AUGepi4s999g14t1           1797/perfect      4     1797  1     4     1797  0     0     0     0     0     0
mq7AUGepir9s999g17t1          2304/mostly       9     2683  1     8     2486  1     0     1     197   0     0.073
mq7PASAgasmbl_218810          1743/perfect      6     1743  1     6     1743  0     0     0     0     0     0
mq7AUGepi4s999g17t1           2215/perfect      5     2215  1     5     2215  0     0     0     0     0     0
mq7PASAgasmbl_218826          0/partial         2     501   0     0     0     0     0     2     501   0     1.000


                *** 7q is worse than 7p *** more exons but less perfect, FN reduced, FP increased
TOTALS          nGene   nExon   bpExon  bpOver  nTrue   bTrue   nFNeg   bFNeg   nFPos   bFPos   pFNeg   pFPos   nPerf   pPerf   Mostly  Joins   Splits
ModelHit        14826   108896  43Mb    32Mb    86277   27Mb    12114   5.5Mb   22619   16Mb    0.164   0.373   1830    0.123   6027    355     1274
ModelMiss       21491   77747   37Mb
OverHit         15204   90652   29Mb
OverMiss        2818    5768    2.6Mb


  rnagene/pasa2_aphid3.asmbl_bestgenes.an7.inmate
PASAgasmbl_218773             1580/perfect      4     1580  1     4     1580  0     0     0     0     0     0
PASAgasmbl_218782             965/perfect       5     965   1     5     965   0     0     0     0     0     0
PASAgasmbl_218826             0/partial         2     501   0     0     0     0     0     2     501   0     1.000

TOTALS          nGene   nExon   bpExon  bpOver  nTrue   bTrue   nFNeg   bFNeg   nFPos   bFPos   pFNeg   pFPos   nPerf   pPerf   Mostly  Joins   Splits
ModelHit        14826   99487   34Mb    31Mb    85791   27Mb    12407   4.6Mb   13696   7.3Mb   0.143   0.216   4522    0.305   8820    343     1023
ModelMiss       10579   28103   7.4Mb
OverHit         15281   90848   29Mb
OverMiss        2741    5572    2.7Mb


  genes/aphid2_epir16b.an7.inmate
AUGepir16bs1249g118t1         3700/mostly       16    3878  1     15    3788  0     0     1     90    0     0.023
AUGepir16bs1250g120t1         6789/mostly       14    7822  1     14    7822  0     1193  0     0     0.132 0
AUGepir16bs1250g121t1         775/partial       14    3454  1     12    1637  0     190   2     1817  0.104 0.526
AUGepir16bs1250g122t1         3684/partial      9     5077  1     6     4350  2     787   3     727   0.153 0.143
AUGepir16bs1250g123t1         648/partial       5     2423  1     2     1826  6     3311  3     597   0.645 0.246

TOTALS          nGene   nExon   bpExon  bpOver  nTrue   bTrue   nFNeg   bFNeg   nFPos   bFPos   pFNeg   pFPos   nPerf   pPerf   Mostly  Joins   Splits
ModelHit        14482   119931  41Mb    32Mb    84468   23Mb    16049   8.9Mb   35463   18Mb    0.276   0.441   339 <   0.023   4236    417     1271
ModelMiss       22849   88608   40Mb
OverHit         15134   90436   29Mb
OverMiss        2888    5984    2.8Mb

# TOTALS      nGene nExon bpEx  nTrue bTrue nFNeg bFNeg nFPos bFPos pFNeg pFPos nPerf pPerf Most  Join Splits
# HitGenes    14482 119931 41Mb  84468 29Mb  16049 7.4Mb 35463 12Mb  0.086 0.295 2276  0.157 6110  418   1271
# MissGenes   22849 88608  40Mb

=item more results

 4811821 May  8 13:11 genes/intabs/aphid2_epir2.an7.introntab
  986087 Apr 27 13:06 genes/intabs/aphid2_epir2.an7.inmate

==> genes/intabs/aphid2_epir2.an7.inmate <==
ModelID                         Quality         nExon   bpExon  nOvGene nTrue   bTrue   nFNeg   bFNeg   nFPos   bFPos   pFNeg   pFPos
AUGepir2p1s1g3t1                1456/partial    3       5938    1       1       2942    1       1554    2       2996    0.346   0.505
AUGepir2p1s1g5t1                626/partial     9       6643    1       6       2041    0       378     3       4602    0.156   0.693
AUGepir2p1s1g6t1                317/partial     7       1684    1       3       731     0       97      4       953     0.117   0.566

==> genes/intabs/aphid2_epir2.an7.introntab <==
GeneID  Score/IntronQual        Valin/Exons     Good/Bad/Join   Joins   Exons_Introns
AUGepir2p1s1g1t1        0/none  0/1     0/0/0           x1:
AUGepir2p1s1g2t1        -2/none 0/3     0/0/0           x1: x2: x3:
AUGepir2p1s1g3t1        -1/poor 1/3     1/0/0   000     x1: x2: x3:intr=+14/-2,N4,N6,N6


=cut