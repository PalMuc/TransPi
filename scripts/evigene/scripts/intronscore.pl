#!/usr/bin/env perl
# intronscore.pl

=item usage

  intronscore.pl < $geneset.gff >  $geneset.inscore.tab

  scripts/overlapfilter.perl -intron2splice=error -pass 'exon,intron' -act markid -midtype scoresum \
  -mark intr -over intron/intron_all.gff.gz -in $geneset | scripts/intronscore.pl > $gbase.introntab

=item update : see inmatescore.pl 
  
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
my($gid, $insplit, $valin, $insum, $nex, $header, @exons);
my %sumscore;

# my $optok= GetOptions(
#   ## "introns=s", \$introns,
#   ## "genes=s", \$genes,
#   #"verbose|v!", \$verbose, 
#   #"norun|n", \$norun, 
#   );
# 
# die "usage: intronscore < genes.gff > table-or-gene.an.gff
#   opts: .. \n"
#   unless($optok and $genes);


while(<>) {
  unless(/^\w/) { next; } #push(@savegene, $_);  # print if outformat = gene.an.gff

  if(/\tmRNA/) { 
    # offer output formats 1. annotate gene.gff, 2. table intron scores per gene row
    
    geneRow( $gid, "none", $nex, $valin, $insum, $insplit, @exons) if($gid);
    
    ($gid)=  m/ID=([^;\s]+)/ or die "ERROR: missing mRNA ID=\n"; 
    ($insplit)= (m/;insplit=([^;\s]+)/) ? $1 : ""; 
    # my($inx)= m/;inexon=([^;\s]+)/;  # NO need, count from exons, intr= 
    
    $nex= $valin= $insum=0; @exons=();
  
  } elsif(/\texon/) { 
    my($pg)= m/Parent=([^;\s]+)/; 
    if($pg ne $gid) { warn "ERROR: exon Parent=$pg in gene=$gid out of order\n"; next; } ## assume input gff is gene-ordered  

    my($intr)=  (m/;intr=([^;\s]+)/) ? $1 : "";     
    # my($inid)=  (m/;inids=([^;\s]+)/) ? $1 : ""; # get ids from intr= ?
    #  # overlapf -mark intr -act markid -midtype scoresum << need to change overlapf to return intr + inid
      
    $nex++; 
    my $xval= "x$nex:";
    ## $xval .= "intr=$intr,$inid" if($intr or $inid);
    $xval .= "intr=$intr" if($intr);
    push @exons, $xval;
    if($intr) {
      my($inscore,@inid)= split ",", $intr;
      my($in1,$in2) = split "/", $inscore;
      $inscore= $in1; $inscore += $in2 if($in2);
      if($inscore > 0) { $valin++; } elsif($inscore < 0) { $valin--; } #? <0 reduce?
      $insum += $inscore; # can be neg
    }
  }
  
}


geneRow( $gid, "LAST", $nex, $valin, $insum, $insplit, @exons) if($gid);

#....................


sub geneRow
{
  my($g, $flag, $nx, $valin, $insum, $insplit, @exons)= @_;

  my($good,$bad,$join,$ginval,$inscore) = (0) x 10;

  #** add summary at end: count classes, sum/ave scores
  
  # inscore = +validsplices -invalidsplices +joined exons -unjoined exons -other errors?
  # * weight inscore by exon-bases to make comparable to other gene scores??
  $inscore= $valin; # == +valid -invalid
  my $bjoin="";
  
  if($valin == 0) { $ginval="none"; $inscore =  1 - $nx; } 
  ##elsif($insum < 0) { $ginval="err.rev"; } 
  ##drop## elsif($insplit) { $ginval="err.split"; } # FIXME, some are alt-splices **
  else { 
    my %lid=(); my @join=(0) x $nx; my $ix=0; 
    foreach my $x (@exons) { 
      ## $x =~ m/x(\d+)/ or next; 
      $ix++;  
      my $j=0; my %id=(); 
      if( my($v)= $x=~/intr=(.+)/ ){ 
        my($iv,@id)=split",",$v; $iv =~ s,/.*,,; 
        $good++ if($iv > 0); 
        $bad++  if($iv < 0); 

        # fixme: new intr=nnn, oid1, oid2 are intron-endpoints; N11 == N12 same intron
        map{ $j=1 if($lid{$_}); $id{$_}=1; } @id; 
        ## also: record intron ids w/ 1 exon end, count cases w/o 2nd exon end == splits

        }
      $join[$ix-1]=$j; $join++ if($j); 
      %lid=%id;  
      }
    
    $bjoin = join"",@join; 
    
    my $notjoin = ($nx - 1) - $join; $notjoin=0 if($notjoin < 0);
    $inscore += $join - $notjoin; #? 
    #  -score for valid introns w/ some joins, but more no-join
    # max inscore= nexons + (nexons-1) = 2 * nexons - 1; for all exons joined w/ introns
    # min inscore= 1 - (nexons-1) =  2 - nexons; for 1 intron, no joins
    # inscore=0 for no introns << should be min
    ## $inscore += $good - $bad; # this is already input from $valin

    if($insum < 0) { $ginval="err.rev"; } # <0 includes retained intron (same strand)
    elsif($bad > 0) { $ginval="poor"; }
    elsif( $nx > 2 and $good == $nx and $join + 1 == $nx) { $ginval="complete"; } 
    elsif( $nx > 3 and $good >= $nx-2 and ($join >= $nx-3 or $join > 0.75 * $nx) ) { 
      my($j1,$j2)=(2,$nx-2); if($join > 8) { $j1++; $j2--; } 
      if( grep { $_ == 0 } @join[$j1..$j2]) { $ginval="miss_inner";} 
      else { $ginval="complete_inner"; } 
      } 
      
      # add? and inscore > 0  for good/ok
    elsif($good and $join) { $ginval= ($valin/$nx > 0.66) ? "good" : "ok"; }  
    else { $ginval="poor"; } 
    
   } 
  
  #? weight inscore by ginval class? 100  for complete, 10 for good/ok, other?
  my $outscore = ($ginval =~ /complete/) ? 100*$inscore 
    : ($ginval =~ /miss_inner|good/) ? 25 * $inscore
    : ($ginval =~ /ok/) ? 10 * $inscore
    : $inscore;
    
    # merge col2/3 for genescore compatibility: Score/IntronQual
  print join("\t", qw(GeneID  Score/IntronQual Valin/Exons Good/Bad/Join Joins Exons_Introns)),"\n"
    unless($header++);  
  print join("\t", $g, "$outscore/$ginval", "$valin/$nx", "$good/$bad/$join", $bjoin, "@exons"), "\n";

#  my($good,$bad,$join,$ginval,$inscore) = (0) x 10;
  $sumscore{'gene'}++;
  $sumscore{'exon'} += $nx;
  $sumscore{'valin'} += $valin;
  $sumscore{'inscore'} += $inscore;
  $sumscore{'good'} += $good;
  $sumscore{'bad'}  += $bad;
  $sumscore{'join'} += $join;
  $sumscore{'quals'}{$ginval}++; 

  if($flag and $flag =~ /LAST/) {
    my @keys= qw( gene inscore valin exon good bad join);
    my @hd  = qw( nGene Score Valin nExon Good Bad Join Qual_counts);
    my @quals= sort keys %{$sumscore{'quals'}};
    @quals = map{ my $c= $sumscore{'quals'}{$_}; "$_=$c"; } @quals;
    print "\n#TOTALS\n";
    print join("\t", "#T", @hd),"\n";    
    print join("\t", "#T", @sumscore{ @keys }, @quals),"\n";

  }

}



__END__

=item use

set genesets=(genes/aphid2_{epir2,epir16b,epi4}*.an7.gff 
  genes/bestgenes_of4.an7g2.gff 
  rnagene/aphid_rnaseq_cufftrin_kinfull.an7.gff 
  rnagene/pasa2_aphid3.asmbl_bestgenes.an7.gff)

foreach geneset ( $genesets )

set otab=`echo $geneset | sed 's/.gff/.introntab/'`

scripts/overlapfilter.perl -intron2splice=error -pass 'exon,intron' -act markid -midtype scoresum \
-mark intr -over intron/intron_all.gff.gz -in $geneset | $evigene/scripts/intronscore.pl > $otab

end

==> genes/aphid2_epi4.an7.introntab <==

#TOTALS
#T      nGene   Score   Valin   nExon   Good    Bad     Join    Qual_counts
#T      35961   98700   92043   187470  96993   4751    79083  
complete=4202   complete_inner=1845     err.rev=1908    good=2143
      miss_inner=1235 none=19818      ok=1487 poor=3323

==> genes/aphid2_epir16b.an7.introntab <==

#TOTALS
#T      nGene   Score   Valin   nExon   Good    Bad     Join    Qual_counts
#T      37331   116304  107876  208539  110903  2900    89818  
complete=3709   complete_inner=4949     err.rev=1110    good=2870
      miss_inner=1339 none=18780      ok=1815 poor=2759

==> genes/aphid2_epir2.an7.introntab <==

#TOTALS
#T      nGene   Score   Valin   nExon   Good    Bad     Join    Qual_counts
#T      32505   111514  93864   177753  98047   4001    81449  
complete=3903   complete_inner=3028     err.rev=1408    good=1998
      miss_inner=1054 none=16911      ok=1392 poor=2811

==> genes/bestgenes_of4.an7g2.introntab <==

#TOTALS
#T      nGene   Score   Valin   nExon   Good    Bad     Join    Qual_counts
#T      14313   157848  98141   118448  99944   1672    81921  
complete=4362   complete_inner=2897     err.rev=10      good=2343
      miss_inner=1469 none=131        ok=969  poor=2132

==> rnagene/aphid_rnaseq_cufftrin_kinfull.an7.introntab <==

#TOTALS
#T      nGene   Score   Valin   nExon   Good    Bad     Join    Qual_counts
#T      13058   124762  68940   69158   69001   31      55961  
complete=8748   complete_inner=119      err.rev=7       good=3671
      miss_inner=2    none=15 ok=13   poor=483

==> rnagene/pasa2_aphid3.asmbl_bestgenes.an7.introntab <==

#TOTALS
#T      nGene   Score   Valin   nExon   Good    Bad     Join    Qual_counts
#T      25405   203591  114654  127590  116401  1491    95561  
complete=11096  complete_inner=463      err.rev=354     good=7778
      miss_inner=468  none=3842       ok=180  poor=1224


=item details

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
