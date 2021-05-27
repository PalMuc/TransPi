#!/usr/bin/env perl
# samtcount.pl

use strict;
use Getopt::Long;

use constant DOT => '.';

my $debug=1;
my $MIN_IDENT = 90;
my $insorted=1;  # fixme
my $informat=0;  # fixme
my $defaultReadGroup = "rgc"; # for feature tag
my $outformat= "table"; # or gff
my %didputid;
my $didhead=0;
my $domeancount=0;
my (%count,%number,%rdgroup); # result per gid,rdgroup
my($overstart,$overend);
my @groupcols;

our $BINSIZE  = 500 ; #was# 5000;

my ($overlaps, $input, $stranded, $orderedoverlap, %ordered_over);
my ($ok, $nover, $nmiss, $n_overlaps)= (0) x 10;

my $optok= GetOptions(
  "input=s", \$input, 
  "informat=s", \$informat,  
  "overlaps=s", \$overlaps, # gff features to annotate and output
  "identity|MIN_IDENT=i", \$MIN_IDENT, 
  "group=s", \$defaultReadGroup, 
  "format=s", \$outformat, 
  "sorted!", \$insorted, 
  "debug!", \$debug, 
  );

die "opt error" unless($optok and $overlaps);

$domeancount=1 if($informat =~ /tiledif|tiling|tilegff/i);

my $ovh;  
   if($overlaps =~ /.gz$/) { $ok= open(OVR,"gunzip -c $overlaps |");  $ovh= *OVR; }
elsif($overlaps =~ /^(stdin|-)/) { $ovh= *STDIN; $ok=1; }
else { $ok= open(OVR,$overlaps); $ovh= *OVR; }
die "bad -overlaps=$overlaps" unless($ok);

my($overlaplist, $featlist, $reflist)
  = collect_overlaps($ovh); close($ovh);

my $inh= *STDIN;
if($input) {
$ok = ($input =~ /.gz$/) ? open($inh,"gunzip -c $input |") 
      : ($input =~ /^(stdin|-)/) ? $inh= *STDIN
      : open($inh,$input);
die "bad -input=$input" unless($ok);
}

my( $nin, $nerr, $nomap );
if($informat =~ /tiledif/i) {
  ( $nin, $nerr, $nomap )= readTilediff($inh);
  
} elsif($informat =~ /tiling|tilegff/i) {
  ( $nin, $nerr, $nomap )= readTilegff($inh);
  
} elsif($informat =~ /sam/i) {
  ( $nin, $nerr, $nomap )= readSam($inh);
  
} else {
  die "need -informat=tilediff|tilegff|sam\n";
}

warn "#samcount: nin=$nin, nerr=$nerr, nomap=$nomap, nover=$nover, nmiss=$nmiss, nfeat=$n_overlaps\n";


#..................................................................

sub readSam 
{
	my($inh)= @_;
	
  my ( $nin, $nerr, $nomap ) = (0) x 10;
  my ( $isover, $atb, $ate, $overids);
  $isover=$atb= $ate= 0; $overids=undef;
  my $lchr="";
  
  while(<$inh>) {
  next unless(/^\w/); # comments?
  
  # ?? add inner loop here for dupl reads, $rdcount++
  # $rdkey= "$flag.$chr.$cloc.$cigar.$seq";
  # if($rdkey eq $lrdkey) { $rdcount++; }
  # else { putlast($lrdkey, $rdcount); }
  
  chomp;
  my $linein=$_;  $nin++;
  my($qid, $flag, $chr, $cloc, $mapq, $cigar, $matechr, $mateloc, $isize, $seq, $qual, @opt)
    = split"\t";
  
#   my $f_pair= ($flag & 0x0001);
#   my $f_pairok= ($flag & 0x0002);
  my $f_mismatch= ($flag & 0x0004);  # also set when chr eq '*'
#   my $f_mismate= ($flag & 0x0008);
  my $f_isrev= ($flag & 0x0010);
#   my $f_revmate= ($flag & 0x0020);
  
  if($chr ne $lchr) {
    put_feats($lchr) if($lchr and $insorted);
    $isover=$atb= $ate= 0; $overids=undef;
  }
  
  $lchr= $chr; 
  do{ $nomap++; next;} if($chr eq '*'); # no map
  
  my $len  = length($seq);
  my $cend = $cloc; # calc correct end from align/intron spans
  my $intype="intron";
  my(@aln,@intr,@itype);
  # fixme for: ^nSnM ... was m/^(\d+)M/
  if($cigar =~ m/^(\d+)M/ ) { @aln=($1); $cend += $1; }
  elsif($cigar =~ /^\d+[HSN]/) { @aln=(0); } # need for nSnM cigar
  while($cigar =~ m/(\d+)([DHNS])(\d+)M/g) { 
    my($bi,$bt,$bx)=($1,$2,$3);
    $bi=0 if($bt eq "H" or $bt eq "S");
    push(@intr,$bi);  push(@itype,$bt);
    push(@aln,$bx); $cend += $bi + $bx; 
    }


# ** add RG readgroup parsing, cmdline option, ...
  my $instrand= DOT; my $inmismat=0;
  my $score= DOT; 
  my $readgroup= $defaultReadGroup;
  
  my $nx= scalar(@aln);
  my $alen=0; map{ $alen+=$_ }@aln; 
  foreach (@opt) { 
    if(m/NM:i:(\d+)/){ $alen -= $1; } 
    elsif(m/XS:A:([\+\-])/){ $instrand=$1; } 
    elsif(m/NS:i:(\d+)/) { $inmismat=$1; } #new, if XS, = Mismatches within min_anchor_len of a splice junction, > 0 is poor
    elsif(m/RG:A:(\S+)/){ $readgroup=$1; } #new
    }
  $score=int(0.5 + 100*$alen/$len);

  do{ $nerr++;  next;} if($score < $MIN_IDENT); # 2err/37bp = 35/37 = 95%; 36/37 = 97%

  my ($sid, $matesid, $matestrand)=(0,0,DOT);
  my $xb= $cloc; my $lin=0;
  for (my $i=0; $i < $nx; $i++) {
    my $aln= $aln[$i]; 
    my $intr=$intr[$i]||0;
    my $isintron= ($itype[$i] eq "N")?1:0;
    my $xe= $xb + $aln - 1; # this is right
    
#     if( $xe <= $ate and $xb >= $atb and $overids) { # iff contained in last feature set
#       addRead2($overids, $chr, $xb, $xe, $readgroup, 1);
#     } else {

      ( $isover, $atb, $ate, $overids)= 
      addRead($chr, $xb, $xe, $instrand,  $readgroup, 1);

#     }
    
#     if($intr>0 and $isintron) {
#       my ($sb,$se);
#       $sb= $xe + 1;  $se= $xe + $intr; 
#       #addIntron( $chr, $sb, $se, $instrand, $intype);  
#       #print "# $chr:$sb-$se:$instrand/$f_isrev\n# $linein\n" 
#       #  if($debug and $testin{c}{$chr} and $testin{b}{$sb} and $testin{e}{$se});
#     }
    
    $xb += $aln + $intr ; # this is right, not offby1
    $lin= $intr;
 	  }
 	  
  }

  put_feats("ALL");  
  return( $nin, $nerr, $nomap );
}


sub addRead
{
  my( $rchr, $rb, $re, $ror, $rdgroup,$count)= @_;  
  my( $nupd, $atb, $ate )=(0,0);
  
  my @overids= overlaps($rchr,$rb,$re,$ror,$rdgroup,$count); # call addRead2 here?
  $atb = $overstart; # _max($atb,$overstart);
  $ate = $overend;   # ($ate==0) ? $overend : _min($ate,$overend);
  
  my $isover= @overids;
  if($isover) {  
    $nover++; 
    ## ( $nupd, $atb, $ate )= addRead2( \@overids, $rchr, $rb, $re, $rdgroup,$count);
  } else {
    $nmiss++;
  }  
 return( $isover, $atb, $ate, \@overids); 
}

# drop this
# sub addRead2
# {
#   my( $overids,$rchr,$rb,$re,$rdgroup,$count)= @_;  
#   my( $nupd, $atb, $ate )=(0,0,0);
#   if(@$overids) { 
#     $nover++; 
#     foreach my $gid (@$overids) {
#       my( $iupd, $itb, $ite )= addRead1($gid,$rchr,$rb,$re,$rdgroup,$count);
#       $nupd += $iupd;
#       $atb = _max($atb,$itb);
#       $ate = ($ate==0) ? $ite : _min($ate,$ite);
#     }
#   } else {
#     $nmiss++;
#   }  
#   return( $nupd, $atb, $ate );
# }
# 
# 
# # drop this
# sub addRead1
# {
#   my($gid,$rchr,$rb,$re,$rdgroup,$count)= @_;
#   my($nupd, $atb, $ate )=(0,0,0);
#   ($overstart,$overend)= (0,0);
# 
#   my $fta = $featlist->{$gid} or return;
#   my $n= scalar(@$fta);
#   # $rdgroup{$rdgroup}++;
#   foreach my $i (0..$n-1) {
#     my($tb,$te,$ref,$gid,$to,$typ,$oid,$inline)= @ { $fta->[$i] };
#     next unless($rb  < $te and $re > $tb and $ref eq $rchr);
#     
#     addRead0($gid,$tb,$te,$rdgroup,$count);
#     $nupd++;
#   }
#   
#   $atb = $overstart; # _max($atb,$overstart);
#   $ate = $overend; # ($ate==0) ? $overend : _min($ate,$overend);
#   return($nupd, $atb, $ate );
# }
# 
# # drop this
# sub addRead0
# {
#   my($gid,$tb,$te,$rdgroup,$count)= @_;
#   $count{$gid}{$rdgroup} += $count;
#   $rdgroup{$rdgroup}++;
#   $overstart= _max($overstart,$tb);
#   $overend= ($overend==0) ? $te : _min($overend,$te);
# }



# convert to array handling? _max(@list) 
sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }

sub overlaps
{
  my($ref,$tb,$te,$to, @rdcount)= @_; ## ($rdgroup, $count)     
  my @lid;
  my %didid=();
  ($overstart,$overend)= (0,0);
  # $count ||= 1;
  
  return unless($overlaplist->{$ref});
  my @bins= (int($tb/$BINSIZE) .. int($te/$BINSIZE)); # 1 bin only for short reads, except rarely at breaks
  foreach my $ib (@bins) {
    $overlaplist->{$ref}{$ib} or next;
    my @locs= @{$overlaplist->{$ref}{$ib}};
    foreach my $rloc (@locs) {
      my ($lb,$le,$lid,$lo,$oid)= @{$rloc}[0,1,3,4,6];
      next if($didid{$oid.$lb.$le}++);
      
      my $over= ($tb <= $le && $te >= $lb) ? 1 : 0;
      $over=0 if($stranded and ($to =~ /[+-]/ and $lo =~ /[+-]/ and $to ne $lo)); # NEW 2010.10; fixme for to or lo eq "."

      if($over) {
#        addRead0($lid,$lb,$le,$rdgroup,$count);
        for(my $k=0; $k<@rdcount; $k+=2) {  # for tilescores
          my($rdgroup,$count)= @rdcount[$k,$k+1];
          if($count ne "NA") {
          $count{$lid}{$rdgroup} += $count;
          $number{$lid}{$rdgroup} ++;
          $rdgroup{$rdgroup}++;
          }
          }
        $overstart= _max($overstart,$tb);
        $overend= ($overend==0) ? $te : _min($overend,$te);
        push @lid, $lid;
        }
      }
    }
  
  if(@lid) {    
    my %lid= map{$_,1}@lid; 
    @lid= sort keys %lid;    
    return @lid; ## return join",", @lid; 
    }
  return;
}


sub put_feats
{
  my( $ref )= @_;  # or put by chr/ref
  my @gid;
  
  # if(ref $ids) { @gid= @$ids; } # elsif
  if($ref eq "ALL") { @gid= sort keys %ordered_over; }  
  elsif($ref) { 
    return unless($reflist->{$ref}); # reads w/o feats
    @gid= @{$reflist->{$ref}};  
    }
  else { return; } # never?
  
  unless(%rdgroup) { $rdgroup{$defaultReadGroup}++; }
  my @rdgroup= sort keys %rdgroup;  
  unless($didhead++) {
    print join("\t","feature_id","Number",@rdgroup),"\n" if($outformat =~ /table/i);
    print "##gff-version 3\n" if($outformat =~ /gff/i);    
  }  
  
  @gid= sort{ $ordered_over{$a} <=> $ordered_over{$b} } @gid;
  foreach my $gid (@gid) {
    next if ( $didputid{$gid}++ );
    
    if($outformat =~ /table/i) {
      my @gcount=(); my $nmax=0;
      if($domeancount) {
      ## take average for tilediff/tilegff, not count sum.
      @gcount= map{ 
        my $sc=$count{$gid}{$_}||0; 
        my $nc=$number{$gid}{$_}||1; #? add this as output column?
        $nmax=$nc if($nc>$nmax);
        sprintf"%.2f", $sc/$nc; } @rdgroup;      
      } else {
      @gcount= map{ 
        my $nc=$number{$gid}{$_}; $nmax=$nc if($nc>$nmax);
        sprintf"%.2f", $count{$gid}{$_} || 0 } @rdgroup;
      }
      print join("\t",$gid,$nmax,@gcount),"\n";
      
    } elsif($outformat =~ /gff/i) {
      my $gcount= join ";", map{ my $c=sprintf"%.2f", $count{$gid}{$_} || 0; "$_=$c" } @rdgroup;
      my $fta= $featlist->{$gid} || [];
      foreach my $ft (@$fta) {
        my($tb,$te,$ref,$gid,$to,$typ,$oid,$inline)= @$ft;
        $inline =~ s/$/;$gcount/;
        print $inline;
      }
    }
  }
  # flush(STDOUT);
}


sub collect_overlaps
{
  my( $ingff)= @_;
  my( %overlaps, %feats, %reflist); # returns
  my $nr=0;
  while(<$ingff>){
    next unless(/^\w/); 
    my $inline= $_;
    chomp;
    my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr)= split"\t"; #,@gffmore
    $tattr ||="";  
    
    #FIX# if($passtypes and "$typ.$src" !~ m/$passtypes/) { print if $printpass; next; } # pass other types
    
    $nr++;
    my $oid= "N".$nr; # save always for uniq feature test
    my($gid); # fixme: separate $gid from markid : want two fields, one for id tests
    if($tattr) {
      if($tattr =~ m/\bID=([^;]+)/) {  $gid=$1; }
      elsif($tattr =~ m/\bParent=([^;]+)/) {  $gid=$1; }
      elsif($tattr =~ m/\bTarget=([^;\s]+)/) { $gid=$1; } # NEW 2010.10
    }
    unless(defined $gid) { $gid = $oid; }

    my $rloc= [$tb,$te,$ref,$gid,$to,$typ,$oid,$inline]; # change to string; save mem

    push( @{$feats{$gid}},$rloc); #? only 1 per id ?
    push( @{$reflist{$ref}},$gid); #? for output by ref?
    $ordered_over{$gid}= $nr; # if $orderedoverlap;
    
    my @bins= (int($tb/$BINSIZE) .. int($te/$BINSIZE));
    foreach my $ib (@bins) { push( @{$overlaps{$ref}{$ib}}, $rloc); }
    }
    
  $n_overlaps=$nr; # global
  warn"#collect_overlaps=$nr\n" if $debug;
  return (\%overlaps, \%feats, \%reflist);
}

=item add tile gff equivalent

melon2.% gzcat $bg/dpulex/tilex/dpxtiles0907/daphnia_metal-mean3.gff.gz | head
scaffold_10000  daphnia_metal   tile    1       50      475,758,-0.467  .       .       ID=s10000FS000000001
scaffold_10000  daphnia_metal   tile    31      80      407,583,-0.358  .       .       ID=s10000FS000000031
scaffold_10000  daphnia_metal   tile    56      105     538,552,-0.026  .       .       ID=s10000FS000000056
scaffold_10000  daphnia_metal   tile    86      136     2730,571,1.564  .       .       ID=s10000FS000000086
scaffold_10000  daphnia_metal   tile    101     150     535,825,-0.433  .       .       ID=s10000FS000000101
scaffold_10000  daphnia_metal   tile    126     175     353,647,-0.604  .       .       ID=s10000FS000000126

=cut

sub readTilegff
{
  my( $ingff)= @_;
  my ( $nin, $nerr, $nomap, $lchr ) = (0) x 10;
  my ( $isover, $atb, $ate, $overids);
  $isover=$atb= $ate= 0; $overids=undef;
  my $readgroup= $defaultReadGroup;

  while(<$ingff>){
    next unless(/^\w/); 
    my($chr,$src,$typ,$xb,$xe,$xp,$xor,$xph,$xattr)= split"\t";  
    my($vexp, $vcon, $vlrat)= split",",$xp; # 2730,571,1.564 = scores expt,cont,log(exp/con)
    $nin++;

    if($chr ne $lchr) {
      put_feats($lchr) if($lchr and $insorted);
      $isover=$atb= $ate= 0; $overids=undef;
    }
    $lchr= $chr; 
    
    ( $isover, $atb, $ate, $overids)= 
      addTileScores($chr, $xb, $xe, $xor,  
        $readgroup, $vexp, $vcon, $vlrat);
    
  }

  put_feats("ALL");  
  return( $nin, $nerr, $nomap );
}

sub addTileScores
{
  my( $rchr, $rb, $re, $ror, $rdgroup,$vx,$vc,$vr)= @_;    
  my @rdcount= ( $rdgroup."N", 1,  $rdgroup."X", $vx, 
    $rdgroup."C", $vc,  $rdgroup."R", $vr);   

  my @overids= overlaps($rchr,$rb,$re,$ror,@rdcount); 
  
  my $isover= @overids;
  if($isover) {  
    $nover++; 
  } else {
    $nmiss++;
  }  
 return( $isover, $overstart, $overend, \@overids); 
}

=item add tilerepdiff table

/export/udisk2/work/daphplx/tiles/signorm0810

SEQ_ID  start   end     Amean   Fem1    Fem2    Fem3    Cha1    Cha2    Cha3    Met1    Met2    Met3
1       5983    6034    9.035   0.1538  -0.6379 0.6633  0.8629  NA      -0.2288 0.9128  -0.213  0.2739
1       26953   27007   9.067   -0.595  NA      0.592   0.0602  0.0405  0.4964  -0.3565 0.5216  0.0921
1       27033   27096   9.139   -0.9009 -0.4583 0.2444  -0.0995 -0.3497 0.9572  -0.3189 0.3154  NA

=cut

sub readTilediff
{
  my( $ingff)= @_;
  my ( $nin, $nerr, $nomap, $lchr, $isover, $atb, $ate ) = (0) x 10;
  # my $readgroup= $defaultReadGroup;
  
  while(<$ingff>){
    chomp; my($chr,$xb,$xe,@vals)= split"\t";  
    if(/^SEQ/) {
      @groupcols=@vals;
    } else {
      next unless(/^\d/); 
    }      
    
    next if($vals[0] eq "NA"); # missing all
    $chr ="scaffold_".$chr; ##if($chr =~ /^\d/);
    $nin++;

    if($chr ne $lchr) {
      unless($lchr or @groupcols) { @groupcols= map{ "Group$_" } (0..$#vals); }
      put_feats($lchr) if($lchr and $insorted);
      $isover=$atb= $ate= 0;  
    }
    $lchr= $chr; 
    
    ( $isover, $atb, $ate)= 
      addTileDiff($chr, $xb, $xe, @vals);
    
  }

  put_feats("ALL");  
  return( $nin, $nerr, $nomap );
}

sub addTileDiff
{
  my( $rchr, $rb, $re, @vals)= @_; 
  
  my($i,$j,@rdcount);
  for($i=0,$j=0; $i<@vals; $i++, $j+=2) {
    $rdcount[$j]= $groupcols[$i];
    my $v= $vals[$i]; 
    # ok now?# $v=0 if($v eq "NA");
    $rdcount[$j+1]= $v;
  }
  
  my @overids= overlaps($rchr,$rb,$re,DOT,@rdcount); 
  
  my $isover= @overids;
  if($isover) {  
    $nover++; 
  } else {
    $nmiss++;
  }  
 return( $isover, $overstart, $overend); # \@overids 
}
