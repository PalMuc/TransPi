#!/usr/bin/env perl
# samcount.pl

use strict;
use Getopt::Long;

use constant DOT => '.';

my $debug=1;
my $MIN_IDENT = 0;  # off for now, use opt
my $insorted=1;  # fixme
my $defaultReadGroup = "rgc"; # for feature tag
my $outformat= "table"; # or gff
my %didputid;
my $didhead=0;
my (%count,%rdgroup); # result per gid,rdgroup
my($overstart,$overend);

our $BINSIZE  = 500 ; #was# 5000;

my ($overlaps, $input, $stranded, $orderedoverlap, %ordered_over);
my ($ok, $nover, $nmiss, $n_overlaps)= (0) x 10;
my (@overmore, $mustdo, $mustnot,$dropids,%dropids,$keepids,%keepids);

my $optok= GetOptions(
  "sam|hits|input=s", \$input, 
  "overlaps=s", \$overlaps, # gff features to annotate and output
  ##"overmore=s", \@overmore, # gff features to annotate and output
  "mustdo=s", \$mustdo, 
  "mustnot=s", \$mustnot, 
  "dropids=s", \$dropids, # file of read ids to skip (dupl matched)
  "keepids=s", \$keepids, # file of read ids to keep (uniques)
  ## add overlap2/overmust/... 2nd gff set which read must match (or must not), after 1st overset
  ## for exons + tiles (i.e. drop untiled exon parts)
  "identity|MIN_IDENT=i", \$MIN_IDENT, 
  "group=s", \$defaultReadGroup, 
  "format=s", \$outformat, 
#  "chr=s", \$CHR, 
#  "start=i", \$CSTART, 
#  "length=i", \$CLENGTH, 
  "sorted!", \$insorted, 
  "debug!", \$debug, 
  );

die "opt error" unless($optok and $overlaps);

die "error, one only: -dropids or -keepids\n" if($dropids and $keepids);
if($dropids) {
  open(F,$dropids) or die "bad -dropids=$dropids"; while(<F>){ chomp; $dropids{$_}=1; } close(F);
}
if($keepids) {
  open(F,$keepids) or die "bad -keepids=$keepids"; while(<F>){ chomp; $keepids{$_}=1; } close(F);
}

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

my( $nin, $nerr, $nomap, $nskipid )= 
readSam($inh);

warn "#samcount: nin=$nin, nerr=$nerr, nomap=$nomap, nover=$nover, nmiss=$nmiss, nskipid=$nskipid nfeat=$n_overlaps\n";


#..................................................................

sub readSam 
{
	my($inh)= @_;
	
  my ( $nin, $nerr, $nomap, $nskipid ) = (0) x 10;
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
  
  do{ $nomap++; next;} if($f_mismatch); ## $chr eq '*'); # no map
  if($keepids) { do{ $nskipid++; next; } unless($keepids{$qid});  }
  elsif($dropids) { do{ $nskipid++; next; } if($dropids{$qid});  }
  
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
    elsif(m/NS:i:(\d+)/) { $inmismat=$1; } #new, if XS, = Mismatches within min_anchor_len ..
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
  return( $nin, $nerr, $nomap,$nskipid );
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
  my($ref,$tb,$te,$to, $rdgroup, $count)= @_;    
  my @lid;
  my %didid=();
  my $mustgot=0;
  my $mustnotgot=0;
  ($overstart,$overend)= (0,0,0);
  $count ||= 1;
  
  return unless($overlaplist->{$ref});
  my @bins= (int($tb/$BINSIZE) .. int($te/$BINSIZE)); # 1 bin only for short reads, except rarely at breaks
  foreach my $ib (@bins) {
    $overlaplist->{$ref}{$ib} or next;
    my @locs= @{$overlaplist->{$ref}{$ib}};
    foreach my $rloc (@locs) {
##     my $rloc= [$tb,$te,$ref,$gid,$to,$typ,$oid,$inline]; # change to string; save mem
      my ($lb,$le,$lid,$lo,$ltyp,$oid)= @{$rloc}[0,1,3,4,5,6];
      next if($didid{$oid.$lb.$le}++);
      
      my $over= ($tb <= $le && $te >= $lb) ? 1 : 0;
      $over=0 if($stranded and ($to =~ /[+-]/ and $lo =~ /[+-]/ and $to ne $lo)); # NEW 2010.10; fixme for to or lo eq "."

      if($over) {
        
        if($mustnot and $ltyp =~ /$mustnot/) {
          $mustnotgot++; last;
        }
        if($mustdo and $ltyp =~ /$mustdo/) {
          $mustgot++; next; # no ids collected
        }  
        
#        addRead0($lid,$lb,$le,$rdgroup,$count);
        $count{$lid}{$rdgroup} += $count;
        $rdgroup{$rdgroup}++;
        $overstart= _max($overstart,$tb);
        $overend= ($overend==0) ? $te : _min($overend,$te);
        push @lid, $lid;
        }
      }
    }
    
  if($mustnotgot > 0 or ($mustdo and $mustgot == 0)) {
    # need to undo any count
    foreach my $lid (@lid) { $count{$lid}{$rdgroup} -= $count;  $rdgroup{$rdgroup}--; }
    return;
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
    print join("\t","feature_id",@rdgroup),"\n" if($outformat =~ /table/i);
    print "##gff-version 3\n" if($outformat =~ /gff/i);    
  }  
  
  @gid= sort{ $ordered_over{$a} <=> $ordered_over{$b} } @gid;
  foreach my $gid (@gid) {
    next if ( $didputid{$gid}++ );
    ## new: skip mustdo, mustnot features ... leave out of ordered_over ?
    
    if($outformat =~ /table/i) {
      my @gcount= map{ $count{$gid}{$_} || 0 } @rdgroup;
      print join("\t",$gid,@gcount),"\n";
      
    } elsif($outformat =~ /gff/i) {
      my $gcount= join ";", map{ my $c=$count{$gid}{$_} || 0; "$_=$c" } @rdgroup;
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
    
    my $skipord= ( ($mustdo and $typ =~ /$mustdo/) or ($mustnot and $typ =~ /$mustnot/) )?1:0;
    unless($skipord) {
      push( @{$reflist{$ref}},$gid); #? for output by ref?
      $ordered_over{$gid}= $nr;
    }
    
    my @bins= (int($tb/$BINSIZE) .. int($te/$BINSIZE));
    foreach my $ib (@bins) { push( @{$overlaps{$ref}{$ib}}, $rloc); }
    }
    
  $n_overlaps=$nr; # global
  warn"#collect_overlaps=$nr\n" if $debug;
  return (\%overlaps, \%feats, \%reflist);
}

