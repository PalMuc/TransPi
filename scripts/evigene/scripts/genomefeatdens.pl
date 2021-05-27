#!/usr/bin/env perl
# genomefeatdens.pl

=item about

genome feature densities
density count per sliding base window

# be careful of mixed types: mRNA, CDS sorted by loc will mess up overlap check

cat genes/cacao11_bestgenes.pub3b.gff | egrep '^scaffold_[1-9]  ' | env src=0 $evigene/scripts/genomefeatdens.pl > map/dens/bestgene3b.bed

gzcat prot/protein_uniq.gff.gz | egrep '^scaffold_[1-9] ' | env feat=CDS src=0 $evigene/scripts/genomefeatdens.pl > map/dens/prot.bed

gzcat est/est.*.mars11.gff.gz | egrep '^scaffold_[1-9]  ' | env src=1 $evigene/scripts/genomefeatdens.pl > ! map/dens/est.bed

gzcat misc/transposon_mrho.gff.gz | sort -k1,1 -k4,4n -k5,5nr | egrep '^scaffold_[1-9]  ' | env src=0 $evigene/scripts/genomefeatdens.pl > map/dens/transposon.bed

gzcat criollo13gdna-mars10asm.bed.gz | env chrs=$caca/genome/asm10/cacao10chrs.fa.count \
snp=base snpscore=10 bin=100000 slide=1 $evigene/scripts/genomefeatdens.pl > criollogdna10.bed8 

=cut

use strict;


my $bins  = $ENV{bin} || 50000;  # 5000; #which default?
my $slide = $ENV{slide} || 0; 
my $snpbed= $ENV{snp} || 0; 
my $inbed = $ENV{bed} || 0; 
my $inbed1 = $ENV{bed1} || 0; # 1-base location, score: scaffold_1      40      6
my $addid = $ENV{addid} || 0; 
my $SPANSCORE = $ENV{span} || 0; #? default on
my $ftype = $ENV{feat} || "exon|transposon|intron";
my $typesource= defined $ENV{src} ? $ENV{src} : 0;
my $SNP_ONLYBASE= ($snpbed == -1 || $snpbed eq "base")?1:0; 
my $SNPSCORE= $ENV{snpscore} || 0;
my $CHRFULL= $ENV{chrfull}||0; # DONT abbreviate (chr|scaffold)[123] as 123 ?
my $addscore= defined $ENV{score} ? $ENV{score} : 1; # exon.c; only gff?: always do this or not?
my $addmiss = defined $ENV{miss} ? $ENV{miss} : 1; # exon.ib; only gff?: always do this or not?

my $usage=
"\nenv bin=$bins bed=$inbed|gff slide=$slide feat='$ftype' chrs=genome/cacao11chrs.fa.count genomefeatdens.pl\n";

my $chrs  = $ENV{chrs} or die "Need chrs=chrsizesfile with name\tsize; $usage";

my $MORESCORES= $ENV{morescores} || ""; # for bed now, add gff?; use ENV{feat} ?
my @MORESCORES= $MORESCORES ? split( /[,;\s]+/, $MORESCORES) : ();

my $bin2= 2*$bins;
my $binoff= int($bins/(1+$slide));
my $binstep= ($slide>0) ? int($bins/(1+$slide)) : $bins;

my (%nc,%ids,%tp,%last);
my($lr,$lb,$le,$lt);

my (@chrs,%chrs,%csteps);
open(C,$chrs) or die "$chrs";
while(<C>){ next unless(/^\w/); my($c,$clen)=split; push @chrs, $c; $chrs{$c}= $clen; } close(C);

foreach my $c (@chrs) {
  my $clen= $chrs{$c};
  # my @stepbounds=(); # my($cstart, $cend)=(0,0);
  for(my $k=1; $k<$clen; ) {
    my($cb,$ce)= ($k-1, $k-1+$bins); $ce=$clen-1 if($ce>=$clen);

    $csteps{$c}{$cb}= $ce;
    # push(@stepbounds, $cb,$ce); # use 0-origin; need only start?
    $k= $k + $binstep;
  }
}

## outofmem: need to dump SnpBed, Bed at chr change
if($snpbed) { readSnpBed(); }
elsif($inbed1) { readBed1(); } 
elsif($inbed) { readBed(); }
else { readGff(); }

putout();

#.........................

sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }

=item readGff

  : input gff features sorted by loc, filter feature type (exon|intron|transposon default)
  : output density/window
    feat.nb = number of bases/window with feature
    feat.ib = number of bases/window between feat (not very good, not inverse)
    feat    = count of features in window
      
  head -20 map/dens/transposon.bed
  scaf    loc     transposon      transposon.ib   transposon.nb
  10r     0       15      19677   20585
  10r     25000   25      33907   14875
  10r     50000   18      41734   10727
  10r     75000   6       47787   1759

=item readGff annots
    
   add annot parsing 
   Parent=Thecc1EG021518t1;
    -- gdna coverage/exon   cvXXX=bases/xbase,sumcov
        score: sumcov > cvXXX.c ; bases > cvXXX.nb ; and score exon.nb
    cvM16_3=1102/1102,33320;cvCriollo_13=581/1102,3645;
    cvPound_7=1102/1102,46005;cvUF273_Type1=1102/1102,39225;
    cvccn51=1102/1102,39425;
      -- struvar needs parsing; ty:XX = type; cnn:1 = strain:nread -- is nread useful here?
        score: stvTY.c 
    stv=ty:BK,ccn51:1,pound7:1,,ty:BK,ccn51:2,,ty:D,ccn51:1,pound7:1,tsh1188:4,,ty:D,ccn51:4,pound7:1,tsh1188:2,,ty:D,pound7:1,uf273:2,ty:D,pound7:1,uf273:3
    
=cut

sub readGff
{
  my @cstarts=();
  # my $addscore= defined $ENV{score} ? $ENV{score} : 1; # exon.c; only gff?: always do this or not?
  # my $addmiss = defined $ENV{miss} ? $ENV{miss} : 1; # exon.ib; only gff?: always do this or not?
  my $morekeys= (@MORESCORES>0) ? join('|',@MORESCORES) : "";
  
while(<>){
  next unless(/^\w/);
  my($r,$src,$ty,$tb,$te,$score,$tor,$tph,$tat)=split"\t";  # FIXME, add score per bed, option?
  my ($id)= m/ID=([^;]+)/;  # addid

  next unless($chrs{$r});
  unless($ty and $te>0) { die "NOT GFF format? $_"; } #?
  next if($ftype and $ty !~ m/($ftype)/);
  
  $ty= ($typesource) ? "$ty.$src" : $ty;
  ($lr,$lb,$le)= @{ $last{$ty} || [-1,0,0] };
  my $ov=($lt eq $ty and $r eq $lr and $tb < $le and $te > $lb) ? 1: 0;
  ($lb,$le)=(0,0) if($lr ne $r);

  my $btype="$ty.nb"; # use feattype prefix?
  my $ctype="$ty.c";  # $addscore
  my $itype="$ty.ib"; # $addmiss
  $score= int($score);
  
  # for CDS, mRNA, ...
  unless($ov) { 
    my $width= 1+$te-$tb;  # was nb
    # my $ib= ($tb>$le)? $tb - 1 - $le : 0; # negative ok?
    my $scorespan= ($SPANSCORE) ? $width * $score : $score;
    
    unless($r eq $lr and @cstarts) {
      @cstarts= sort{$a <=> $b} keys %{$csteps{$r}}; # new
    }
    
    my %moresc=(); my @moresc=(); 
    if($morekeys) {
      my @kv= $tat =~ m/\b((?:$morekeys)=[^;\s]+)/g;
      map{ my($k,$v)=split"="; 
        ($v)= $v =~ m/([+-\d]+)/; #need some field parsing for this to work well
        if($k and defined $v) { $k="$k.$src" if($typesource); $moresc{$k}= int($v); }
      } @kv;
      @moresc= sort keys %moresc; 
    }
    
    ## my @ibins= (int($tb/$binstep) .. int( $te/$binstep));
    ## match bn to csteps{$c}{xxx}
    my $j=0; 
    my @ibins= grep { $_ <= $te && $csteps{$r}{$_} >= $tb } @cstarts;  
    foreach my $bn (@ibins) {
      my $be= $csteps{$r}{$bn}; $j++;
      if($j==1 and $le < $bn) { $le=$bn-1; } #? 
 
      my $b1= _max($bn,$tb-1); my $e1= _min($be,$te);
      my $nb= _max(0,$e1-$b1); # not +1 here, b1 is 0-base
      #bad#my $ib= _max(0, $b1 - 1 - $le);
      ## le also bad when missing; need b1 - _max(lastbin,le)
      my $ib= ($j>1)? 0 : _max(0, $b1 - 1 - $le);
      
      $nc{$r}{$bn}{$btype} += $nb;
      $nc{$r}{$bn}{$ctype} += $scorespan if($addscore);  ## add ctype for score count... FIXME here or user: true count is sc * span, per below bed1 diff
      $nc{$r}{$bn}{$itype} += $ib if($addmiss); # dont usually want this one; option?
      $nc{$r}{$bn}{$ty}++; 
      $ids{$r}{$bn}{$id}++ if($addid); # allow 2+
      
      # option to process annots $tat for more scoring types
      if(@moresc) {
        foreach my $ktype (@moresc) {
        $nc{$r}{$bn}{$ktype} += $moresc{$ktype};   
        $tp{$ktype}++; 
        }
      }
      
    }

    $tp{$ty}++;
    $tp{$btype}++;
    $tp{$ctype}++ if($addscore);
    $tp{$itype}++ if($addmiss);

    ($lr,$lb,$le,$lt)=($r,$tb,$te,$ty);
    $last{$ty}=[$r,$tb,$te];
  }

}

}


sub readBed
{
  # my($type)=@_;
  my $type = $ENV{feat} || "bed";
  my $btype="$type.nb"; # use feattype prefix?
  my $ctype="$type.c";  # $addscore
  my $itype="$type.ib";  # $addmiss
  my @cstarts=();

  $tp{$btype}++;
  $tp{$ctype}++;
  map { $tp{$_}++; } @MORESCORES;
  
  my($lr,$lb,$le,$cstart,$cend,$bn)= (0) x 10;
  my $ov=0; #?? not for bed? ( $r eq $lr and $tb < $le and $te > $lb) ? 1: 0;
  $lr="nonesuch";
 
  while(<>){
    next unless(/^\w/); chomp; 
    my($r,$tb,$te,$sc,@sc2)=split"\t";  # $rc,$ps # ps always > 0, always snp at this loc

    #below# next unless($chrs{$r});

    ## OPTION: process @sc2 other scores, but need option of type list for those.
    ## only one type, dont need last{t} hash
    #No# my $t= $type; ##($typesource) ? "$ty.$s" : $ty;
    #no# ($lr,$lb,$le)= @{ $last{$t} || [-1,0,0] };

    #b ($lb,$le)=(0,0) if($lr ne $r);
    ## if($te==$tb or $te==$tb+1) { } 1-base snp ?
    #b my $width=1+$te-$tb;
    #b my $scorespan= ($SPANSCORE) ? $width * $sc : $sc;
    
    if(1) { # unless($ov) 
      unless($r eq $lr and @cstarts) {
        next unless($chrs{$r});
        @cstarts= sort{$a <=> $b} keys %{$csteps{$r}}; # new
      }
      
      ($lb,$le)=(0,0) if($lr ne $r);
      my $width=1+$te-$tb;
      my $scorespan= ($SPANSCORE) ? $width * $sc : $sc;

      my @ibins= grep { $_ <= $te && $csteps{$r}{$_} >= $tb } @cstarts;  
      foreach my $bn (@ibins) {
        my $be= $csteps{$r}{$bn};
        my $b1= _max($bn,$tb); my $e1= _min($be,$te);
        my $nb= _max(0,1+$e1-$b1);
        my $ib= _max(0, $b1 - 1 - $le);
        
        $nc{$r}{$bn}{$itype} += $ib; 
        $nc{$r}{$bn}{$btype} += $nb;  #  $t.".nb" == btype
        $nc{$r}{$bn}{$ctype} += $scorespan;  ## add ctype for score count... FIXME here or user: true count is sc * span, per below bed1 diff
        $nc{$r}{$bn}{$type}++; 

      ## OPTION: process @sc2 other scores, but need option of type list for those.
        if(@MORESCORES and @sc2) {
          for(my $j=0; $j < @sc2; $j++) { 
            my $ktype=$MORESCORES[$j];
            my $kval=$sc2[$j]; 
            # if($SPANSCORE) { $kval= $width * $kval; }
            $nc{$r}{$bn}{$ktype} += $kval;   
            $tp{$ktype}++; # above?
          }
        }

      }
  
      $tp{$type}++;
      $tp{$btype}++;
      $tp{$itype}++;
    
      ($lr,$lb,$le)=($r,$tb,$te);
      #no# $last{$t}=[$r,$tb,$te]; # dont need here
    }
  }
  #///////////    
  
}



sub readBed1
{
  my @cstarts=();
  my $type= $ENV{feat} || 'cov';
  my $btype="$type.nb"; # use feattype prefix?
  my $ctype="$type.c";  # $addscore
  my $itype="$type.ib";

  $tp{$btype}++;
  $tp{$ctype}++; # unless($SNP_ONLYBASE); 
  my($lr,$lb,$cstart,$cend,$bn)= (0) x 10;
  
  #FIXME? is there an edge bug here, bn == 0 when it shouldnt be, Or bn == max / not== max when should?
  
  while(<>){
    next unless(/^\w/);
    my($r,$tb,$score)=split"\t";  #only: scaffold_1      40      6
    next if($tb<1);
    
    unless($r eq $lr and @cstarts) {
      next unless($chrs{$r}); # or last if sorted..
      @cstarts= sort{$a <=> $b} keys %{$csteps{$r}}; # new
      $lr= $r; $cend=$cstart=0; $lb=0;
    }
    
    unless($tb >= $cstart and $tb <= $cend) {
      ($bn)= grep { $_ <= $tb && $csteps{$r}{$_} >= $tb } @cstarts;  
      #? fail if bn missing ..
      $cstart= $bn;
      $cend= $csteps{$r}{$bn};
      next unless($tb >= $cstart and $tb <= $cend); # bugfix?
    }
  
    $nc{$r}{$bn}{$btype} ++; # count any snp at basepos 
    $nc{$r}{$bn}{$ctype} += $score; 
    if($lb>0 and $tb > $lb+1) {
      $nc{$r}{$bn}{$itype} += $tb - $lb; 
    }
    
    $lb= $tb;
  }
}


=item SnpBed

  : output: per window (50kb) count basepos w/ any snp; pct = count/windowsize
  : input value per base, where snp exists, is
    chr, basepos, refbase, pctOfLinesWithSnp, count of snp types ( -del, A,C,G,T)
    
  gzcat datastore/snps/unfdSNPs_10_sets.bed.gz | head
  scaffold_1      30      C       10      0-      0A      9C      0G      1T
  scaffold_1      33      A       10      0-      9A      1C      0G      0T
  scaffold_1      37      C       20      0-      0A      8C      0G      2T
  scaffold_1      40      A       10      0-      9A      1C      0G      0T
  scaffold_1      44      C       20      0-      0A      9C      0G      2T
  scaffold_1      49      C       10      0-      0A      10C     0G      1T
  scaffold_1      51      C       60      0-      0A      5C      0G      6T

  output:  .nb = bases in window with any snp;
           .c =  sum of non-ref bases in window? (max would be nb/window * nlines=10 * nhetero ?)
  scaf    loc     snps.c  snps.nb
  10r     0       15      19677 
  10r     25000   25      33907
  10r     50000   18      41734
  10r     75000   6       47787

=cut 

sub readSnpBed
{
  my @cstarts=();
  my $type= $ENV{feat} || 'snp';
  my $btype="$type.nb"; # use feattype prefix?
  my $ctype="$type.c";  # $addscore
  my $itype="$type.ib";  # addmiss

  $tp{$btype}++;
  $tp{$ctype}++ unless($SNP_ONLYBASE); 
  
  my($lr,$lb,$cstart,$cend,$bn)= (0) x 10;
  my %imap= ( '-' => 0, A=>1, C=>2, G=>3, T=>4); # data specific cols
  
  while(<>){
    next unless(/^\w/);
    my($r,$tb,$rc,$ps,@sc)=split"\t";  # ps always > 0, always snp at this loc
    
    unless($r eq $lr and @cstarts) {
      next unless($chrs{$r}); # or last if sorted..
      @cstarts= sort{$a <=> $b} keys %{$csteps{$r}}; # new
      $lr= $r; $cend=$cstart=0;
    }
    
    unless($tb >= $cstart and $tb <= $cend) {
      ($bn)= grep { $_ <= $tb && $csteps{$r}{$_} >= $tb } @cstarts;  
      $cstart= $bn;
      $cend= $csteps{$r}{$bn};
      next unless($tb >= $cstart and $tb <= $cend); # bugfix?
    }

    next if($SNPSCORE and $ps < $SNPSCORE);
    # $nc{$r}{$bn}{$btype} ++; # count any snp at basepos 
  
    my $nb=1;
    if($SNP_ONLYBASE and $rc > $tb and $rc < $tb+100) {
      my $te= _min($rc,$cend);
      $nb= 1 + $te - $tb;
    }
    $nc{$r}{$bn}{$btype} += $nb; # count any snp at basepos 
  
    unless($SNP_ONLYBASE) { # or @sc == 0 or !defined(@sc)
    my $nsnp=0; 
    my $ir= $imap{$rc} || 9;
    foreach my $i (0..$#sc) { $nsnp += int($sc[$i]) unless($ir==$i); }
    $nc{$r}{$bn}{$ctype} += $nsnp; 
    }
    
    if($lb>0 and $tb > $lb+1) {
      $nc{$r}{$bn}{$itype} += $tb - $lb; 
    }
    $lb= $tb;
    
  }

}



# sub readGff_OLD
# {
# 
# while(<>){
#   my($r,$s,$ty,$b,$e)=split"\t"; 
#   
#   # next if($ftype and $ty !~ m/($ftype)/);
#   
#   my $t= ($typesource) ? "$ty.$s" : $ty;
#   my $rorig=$r;
#   $r =~ s/(contig|chr|scaffold|super)\D*//i; # numeric part only?
#   $r= $rorig unless($r);
#   
#   ($lr,$lb,$le)= @{ $last{$t} || [-1,0,0] };
#   my $ov=($lt eq $t and $r eq $lr and $b < $le and $e > $lb) ? 1: 0;
#   ($lb,$le)=(0,0) if($lr ne $r);
#   
#   # for CDS, mRNA, ...
#   unless($ov)  { 
#     my $nb=1+$e-$b; 
#     my $ib=$b - $le - 1; # negative ok?
#     
#     # this is one window, add slide to bin2 = $binoff + $bn ?
#     my $bn= $bins * int(($b+$e)/$bin2);
#     # fixme: b,e spans 2 bins
#     
#     #? ($id)=m/(?:ID|Parent)=([^;]+)/;
#     for( my $i= 0; $i<=$slide; $i++) {
#       $bn= $bins * int( ( $i*$binoff + $b+$e)/$bin2) if($i>0);
#       $nc{$r}{$bn}{$t.".ib"} += $ib; # if CDS && id == lid; i.e. introns only
#       $nc{$r}{$bn}{$t.".nb"} += $nb; 
#       $nc{$r}{$bn}{$t}++; 
#       ## $bn += $binoff;
#       }
#     $tp{$t}++;
#     $tp{$t.".nb"}++;
#     $tp{$t.".ib"}++;
#   
#   ($lr,$lb,$le,$lt)=($r,$b,$e,$t);
#   $last{$t}=[$r,$b,$e];
#   }
# 
# }
# 
# }

sub putout 
{
  my @tp= sort keys %tp; # these are columns!
  my @id= ($addid)? qw(ids) : (); 
  print join("\t","scaf","loc",@tp,@id),"\n";
  
  my $slide1=1+$slide;
  # my @r= sort{$a<=>$b}keys %nc; # numeric scaffold# changed
  my @r= @chrs;
  foreach my $r (@r) {  
    my $cnum=$r; 
    $cnum =~ s/^(contig|chr|scaffold|super)\D*(\d+)/$2/i unless($CHRFULL); # numeric part only?
    # my @b= sort{$a<=>$b}keys %{$nc{$r}}; 
    my @b= sort{$a <=> $b} keys %{$csteps{$r}}; # new
    foreach my $b (@b) {
      print "$cnum\t$b";
      foreach my $t (@tp) { 
        my $nc= $nc{$r}{$b}{$t} || 0;  
        # $nc=int($nc/$slide1) if($slide>0); # is this wrong now?
        print "\t",$nc; }
      if($addid) { my @id=sort keys %{$ids{$r}{$b}}; @id=qw(na) unless(@id); print "\t",join(",",@id); }
      print "\n";
      } 
  } 
}
