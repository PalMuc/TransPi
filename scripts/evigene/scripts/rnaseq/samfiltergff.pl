#!/usr/bin/env perl
# samfiltergff.pl

=item about

  filter reads.sam that overlap filter.gff (i.e. exons, rRNA high-express repeats, ...)
  optional output: either/both overlap set and nonover set
  
=item options

  -input      : input.sam
  -sorted     : input is sorted (assumed true?)
  -overlaps   : filter.gff
  -notfiltered : output file for nooverlap-set, default = STDOUT, null to skip
  -filtered   : output file for overlap-set
  -mateloc    : filter based on mate location, if not mapped but mate is mapped
  -pctover=50 : require percent base overlap to call overlap
  -expandover=50 :  bp to extend overlap ends; e.g. for partial maps of rRNA frags
  -hardclip    : convert cigar nS to nH, clipping seq, for cufflinks bug

=item usage

set workd=/export/udisk3/work/aphid/

foreach sam (*.bam.sam)
  set nam=`basename $sam .bam.sam`
  if ( -f $nam.norna.sam ) continue
  $workd/scripts/samfiltergff.pl  -mate -in $sam -over $workd/misc/repeat_rrna.gff \
     -filter $nam.sam.rrna -notfilter $nam.norna.sam

end

=cut

use strict;
use Getopt::Long;

use constant DOT => '.';

my $debug=1;
my $MIN_IDENT = 90;
# my $insorted=1;  # not used
# my $defaultReadGroup = "rgc"; # for feature tag
# my $outformat= "table"; # or gff
my $didhead=0;
# my %didputid;
# my (%count,%rdgroup); # result per gid,rdgroup
# my($overstart,$overend);

our $BINSIZE  = 500 ; #was# 5000;

my ($overlaps, $input, $outfilterfile, $outnotfile, $stranded, $orderedoverlap, %ordered_over);
my ($ok, $domateloc, $dohardclip, $expandover, $pctover, $nover, $nmiss, $n_overlaps)= (0) x 20;
my ($dropids,%dropids,$keepids,%keepids);

my $optok= GetOptions(
  "input|sam=s", \$input, 
  "overlaps=s", \$overlaps, # gff features to annotate and output
  "expandover=i", \$expandover, # bp to extend overlap ends; e.g. for partial maps of rRNA frags
  "hardclip!", \$dohardclip, # convert cigar nS to nH, clipping seq, for cufflinks bug
 # "identity|MIN_IDENT=i", \$MIN_IDENT, 
  # "group=s", \$defaultReadGroup, 
  # "format=s", \$outformat, 
  "filtered=s", \$outfilterfile, 
  "notfiltered=s", \$outnotfile, #? or stdout
  "dropids=s", \$dropids, # file of read ids to skip (dupl matched)
  "keepids=s", \$keepids, # file of read ids to keep (uniques)
  "pctover=i", \$pctover, 
  "mateloc!", \$domateloc, 
  # "sorted!", \$insorted, 
  "debug!", \$debug, 
  );

die "opt error" unless($optok and $overlaps);

$pctover= $pctover/100.0 if($pctover);

die "error, one only: -dropids or -keepids\n" if($dropids and $keepids);
if($dropids){open(F,$dropids) or die "bad -dropids=$dropids"; while(<F>){chomp; $dropids{$_}=1;} close(F);}
if($keepids){open(F,$keepids) or die "bad -keepids=$keepids"; while(<F>){chomp; $keepids{$_}=1;} close(F);}

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

my $outnot= *STDOUT; # need undef option
if($outnotfile and $outnotfile !~ m/^(stdout|-)/) {
  if($outnotfile =~ m/^(null|undef|no)$/i) { $outnot=undef; $ok=1; }
  else {
  $ok = open($outnot,">$outnotfile");
  die "bad -notfiltered=$outnotfile" unless($ok);
  }
}

my $outfilter= undef;
if($outfilterfile) {
  $ok=0;
  if(($outfilterfile =~ /^(stdout|-)/) and ! defined $outnot) { $outfilter= *STDOUT; $ok=1; }
  else {  $ok = open($outfilter,">$outfilterfile"); }
  die "bad -filtered=$outfilterfile" unless($ok);
}

my( $nin, $nerr, $nomap )= readSam($inh);

close($outnot) if ($outnot);
close($outfilter) if ($outfilter);
  
warn "#samcount: nin=$nin, nomap=$nomap, nover=$nover, nmiss=$nmiss, nfeat=$n_overlaps\n";
# loqual=$nerr, 

#..................................................................

sub readSam 
{
	my($inh)= @_;
	
  my ( $nin, $nerr, $nomap ) = (0) x 10;
  my ( @overat, $isover, $atb, $ate, $overids);
  my $lchr="";
  
  while(<$inh>) {
  next unless(/^\w/); # comments ; @SQ stuff... keep or not?
    
  # chomp;
  my $linein=$_;  $nin++;
  my($qid, $flag, $chr, $cloc, $mapq, $cigar, $matechr, $mateloc, $isize, $seq, $qual, @opt)
    = split"\t";
  
  my $f_pair= ($flag & 0x0001);
  my $f_pairok= ($flag & 0x0002);
  my $f_mismatch= ($flag & 0x0004);  # also set when chr eq '*'
  my $f_mismate= ($flag & 0x0008);
#   my $f_isrev= ($flag & 0x0010);
#   my $f_revmate= ($flag & 0x0020);
    
  if($keepids) { do{ $nomap++; next; } unless($keepids{$qid});  }
  elsif($dropids) { do{ $nomap++; next; } if($dropids{$qid});  }

  $isover= 0;
  $lchr= $chr;
  my $len  = length($seq); # not if ($seq eq "*");
  
  # do{ $nomap++; next;} if($chr eq '*'); # f_mismatch == no map
  if($f_mismatch) { # ($chr eq '*')
    $nomap++;
    
    #?? should domateloc check BOTH locations and filter if EITHER matches gff?
    if($domateloc and $f_pair and ! $f_mismate) {
      my $matee= $mateloc + $len; # ok?
      if( ! $overlaplist->{$matechr} ) {
        $isover=0; @overat=();
      } elsif(@overat and over1($matechr,$mateloc,$matee,DOT, @overat) ) {
        $isover=1;
      } else {
        @overat= overlaps($matechr,$mateloc,$matee,DOT); 
        $isover= @overat;
      }
    }
    
  } else {  # read mapped
    
    my $cend = $cloc; # calc correct end from align/intron spans
    my $intype="intron";
    my(@aln,@intr,@itype);
    
    # ** H is treated odd cuz seq is clipped to it
    # ** S softclip also odd, $cloc is at NON-clipped seq start.. and cend = cloc + seqlen - end S
    
    # dohardclip: convert end dS..dS to dH, AND clip seq  ALSO clip qual, if exists
    if($dohardclip and ($cigar =~ m/\dS/)) { 
      my ($oldcig,$cl,$cr)= ($cigar,0,0);
      if($cigar =~ s/(\d+)S$/$1H/) { $cr=$1; } # $seq=substr($seq,0,-$cr); 
      if($cigar =~ s/^(\d+)S/$1H/) { $cl=$1; } # $seq=substr($seq,$cl); 
      ##bad for qual# $linein =~ s/\t$oldcig\t/\t$cigar\t/; 
      my @linein=split"\t",$linein;
      $linein[5]= $cigar;
      if($cr or $cl) {
        if($seq ne "*") { 
         my $oldseq=$seq;
         $seq=substr($seq,0,-$cr) if $cr;  $seq=substr($seq,$cl) if $cl; 
         # $linein =~ s/\t$oldseq\t/\t$seq\t/; 
         $linein[9]= $seq;
         $len  = length($seq); 
         }
        if($qual ne "*") { 
         my $oldqual=$qual;
         $qual=substr($qual,0,-$cr) if $cr; $qual=substr($qual,$cl) if $cl; 
         #badddd# $linein =~ s/\t$oldqual\t/\t$qual\t/; 
         $linein[10]= $qual;
         }
       }
    $linein= join"\t", @linein; 
    }

# Op BAM Description : SAM Cigar ops, 2010
# M 0 alignment match (can be a sequence match or mismatch) 
# I 1 insertion to the reference 
# D 2 deletion from the reference 
# N 3 skipped region from the reference 
# S 4 soft clipping (clipped sequences present in SEQ) 
# H 5 hard clipping (clipped sequences NOT present in SEQ) 
# P 6 padding (silent deletion from padded reference) 
# = 7 sequence match    << M eqiv
# X 8 sequence mismatch << M eqiv
    
    my $ciglen=0;
    # $cigar =~ s/[X=]/M/g; # for now
    if($cigar =~ m/^(\d+)M/ ) { my $m=$1; @aln=($m); $cend += $m; $ciglen+=$m; }
    elsif($cigar =~ /^(\d+)[NHSDIP]/) { @aln=(0); } # padding, need for nSnM cigar
    while($cigar =~ m/(\d+)([NHSDIP])(\d+)M/g) { 
      my($bi,$bt,$bx)=($1,$2,$3);      
      $ciglen += $bx; $ciglen += $bi if($bt =~ /SI/);
      $bi=0 if($bt =~ /HSI/); # $bt eq "H" or $bt eq "S"
      push(@intr,$bi);  push(@itype,$bt);
      push(@aln,$bx); 
      $cend += $bi + $bx; 
      }
    $len= $ciglen if($seq eq "*"); #?? do this for matee above
    
    my $instrand= DOT; my $inmismat=0;
    my $score= 0; 
    # ?? add RG readgroup parsing, cmdline option, ...
    # my $readgroup= $defaultReadGroup;
    
    my $nx= scalar(@aln);
    my $alen=0; map{ $alen+=$_ }@aln; 
    foreach (@opt) { 
      if(m/NM:i:(\d+)/){ $alen -= $1; } 
      elsif(m/XS:A:([\+\-])/){ $instrand=$1; } 
      # elsif(m/NS:i:(\d+)/) { $inmismat=$1; } #new, if XS, = Mismatches within min_anchor_len of a splice junction, > 0 is poor
      # elsif(m/RG:A:(\S+)/){ $readgroup=$1; } #new
      }
      
    # $score=int(0.5 + 100*$alen/$len);
    # do{ $nerr++; } if($score < $MIN_IDENT); # next;  2err/37bp = 35/37 = 95%; 36/37 = 97%

    if( ! $overlaplist->{$chr} ) {
      $isover=0; @overat=();
    } else {
    
    my $xb= $cloc;  
    for (my $i=0; !$isover && $i < $nx; $i++) {
      my $aln= $aln[$i]; 
      my $intr=$intr[$i]||0;
      my $isintron= ($itype[$i] eq "N")?1:0;
      if($aln > 1) {
      my $xe= $xb + $aln - 1; # this is right; not if aln == 0
      
      if(@overat and over1($chr,$xb,$xe,$instrand, @overat) ) {
        $isover=1;
      } else {
        @overat= overlaps($chr,$xb,$xe,$instrand); 
        $isover= @overat;
      }
      last if $isover;
      }
      
      $xb += $aln + $intr ; # this is right, not offby1
      # $lin= $intr;
      }
     }
     
    } # end mapped read
 	 
 	  if($isover) { print $outfilter $linein if($outfilter); $nover++; }
 	  else { print $outnot $linein if($outnot); $nmiss++; }
  }

  return( $nin, $nerr, $nomap );
}



# convert to array handling? _max(@list) 
sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }

sub over1 {
  my($ref,$tb,$te,$to, $lref, $lb, $le, $lo, $lid)= @_;    # , $rdgroup, $count
  my $over= ($ref eq $lref && $tb <= $le && $te >= $lb) ? 1 : 0;
  if($over and $pctover) {
    my ($bb,$be)= ( _max($tb,$lb), _min($te,$le) );
    my $maxo= 1+abs($be - $bb);
    my $leno= _min( abs($le - $lb), abs($te - $tb)) || 1;
    my $pover= $maxo/$leno;
    $over = 0 if $pover < $pctover;
  }
  return $over;
}

sub overlaps
{
  my($ref,$tb,$te,$to)= @_;    # , $rdgroup, $count
  my @lid;
  my %didid=();
  # ($overstart,$overend)= (0,0);
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
      #? $over=0 if($stranded and ($to =~ /[+-]/ and $lo =~ /[+-]/ and $to ne $lo)); # NEW 2010.10; fixme for to or lo eq "."

      if($over and $pctover) {
        my ($bb,$be)= ( _max($tb,$lb), _min($te,$le) );
        my $maxo= 1+abs($be - $bb);
        my $leno= _min( abs($le - $lb), abs($te - $tb)) || 1;
        my $pover= $maxo/$leno;
        $over = 0 if $pover < $pctover;
      }
      
      if($over) {
        # return here:  current over span, id
        return( $ref, $lb, $le, $lo, $lid);
        
#         $count{$lid}{$rdgroup} += $count;
#         $rdgroup{$rdgroup}++;
#         $overstart= _max($overstart,$tb);
#         $overend= ($overend==0) ? $te : _min($overend,$te);
#         push @lid, $lid;
        }
      }
    }
  
#   if(@lid) {    
#     my %lid= map{$_,1}@lid; 
#     @lid= sort keys %lid;    
#     return @lid; ## return join",", @lid; 
#     }
    
  return;
}


# sub put_feats
# {
#   my( $ref )= @_;  # or put by chr/ref
#   my @gid;
#   
#   # if(ref $ids) { @gid= @$ids; } # elsif
#   if($ref eq "ALL") { @gid= sort keys %ordered_over; }  
#   elsif($ref) { 
#     return unless($reflist->{$ref}); # reads w/o feats
#     @gid= @{$reflist->{$ref}};  
#     }
#   else { return; } # never?
#   
#   my @rdgroup= sort keys %rdgroup;  
#   unless($didhead++) {
#     print join("\t","feature_id",@rdgroup),"\n" if($outformat =~ /table/i);
#     print "##gff-version 3\n" if($outformat =~ /gff/i);    
#   }  
#   
#   @gid= sort{ $ordered_over{$a} <=> $ordered_over{$b} } @gid;
#   foreach my $gid (@gid) {
#     next if ( $didputid{$gid}++ );
#     
#     if($outformat =~ /table/i) {
#       my @gcount= map{ $count{$gid}{$_} || 0 } @rdgroup;
#       print join("\t",$gid,@gcount),"\n";
#       
#     } elsif($outformat =~ /gff/i) {
#       my $gcount= join ";", map{ my $c=$count{$gid}{$_} || 0; "$_=$c" } @rdgroup;
#       my $fta= $featlist->{$gid} || [];
#       foreach my $ft (@$fta) {
#         my($tb,$te,$ref,$gid,$to,$typ,$oid,$inline)= @$ft;
#         $inline =~ s/$/;$gcount/;
#         print $inline;
#       }
#     }
#   }
#   # flush(STDOUT);
# }
# 

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

    if($expandover) { $tb -= $expandover; $te += $expandover; }
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

