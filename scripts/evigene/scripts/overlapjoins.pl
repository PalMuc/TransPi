#!/usr/bin/env perl

=item about overlapjoins.pl

  -input=introns.gff, mates.gff, any span
  -over=exons.gff, locations to be joined from input
  -output=exons-joined.gff ..

   cut from overlapfilter.perl

=item revise 

  - measure gene/exon span covered by mates (and dangling mates)
  
  - prefered result would be gene exons w/ mate cover/uncover spans indicated
    (including holes in exon cover that may be introns), plus some
    end-of-gene measure of dangling mates: how much of gene model is missing?
    
gzcat velmapt7/aphid_rnaseq.all27cuff8.gff.gz | grep '  exon' | grep '\.1;' |\
$evigene/scripts/overlapjoins.pl -over stdin -in rnas/cleaned/allpe_matec.tab -format bed -strand -sorted \
> aphid_rnaseq.all27cuff8.matejoin.gff

# rnas/cleaned/allpe_matec.tab = mate endpoints, +/-/. strand, .bed format
  
=cut


use strict;
use Getopt::Long;

use constant DOT => '.';

my $debug=1;
our $BINSIZE  = 500 ; #was# 5000;

## these should be options
my $MAX_JOIN_SPAN = 500000;  # we can have longer introns; want reliability though
my $LARGE_JOIN_SPAN = 20000; # require more joiners to count
my $MIN_NJOINER = 3;  
my $MINOR_NJOIN = 0.05;
my $MIN_LONGJOINER = 12;

use constant SPLICE   => 3; # bp for intron splice span tests
use constant SPLICEX  => SPLICE - 1; # for exon
use constant{ jBEGIN => 0, jEND => 1, jCHR => 2, jGID => 3, jSTRAND => 4,
              jTYPE => 5, jOID => 6, jATTR => 7, jLINE => 8, jJOINID => 9, jJOINS => 10 }; ## gff overlap record
              
use constant{ inGFF => 0, inBED => 1, outJOINS => 0, outGFF => 1 };

my ($overlaps, $input, $passtypes, ) = ("") x 10;
my ($ok, $domateloc, $didhead, $stranded, $insorted, $informat, $outformat, $intron2splice, $n_overlaps)= (0) x 20;
my $JGID= 0;

my $optok= GetOptions(
  "input=s", \$input, 
  "format=s", \$informat, 
  "outformat=s", \$outformat, 
  "overlaps=s", \$overlaps, # gff features to annotate and output
#  "mateloc!", \$domateloc, 
  "sorted!", \$insorted, 
  "passtypes=s", \$passtypes,  
  "stranded!", \$stranded,  
  "debug!", \$debug, 
  );

die "opt error" unless($optok and $overlaps);

$passtypes =~ s/[,]/\|/g;
$intron2splice=1 if($passtypes =~ /intron/);
$intron2splice=1 if($informat =~ /intron/);

if($informat =~ /bed|mate/) { $informat= inBED; }
elsif($informat =~ /gff/) { $informat= inGFF; }
else { $informat= inGFF; }

if($outformat =~ /gff/) { $outformat= outGFF; }
else { $outformat= outJOINS; }

my $ovh;  
   if($overlaps =~ /.gz$/) { $ok= open(OVR,"gunzip -c $overlaps |");  $ovh= *OVR; }
elsif($overlaps =~ /^(stdin|-)/) { $ovh= *STDIN; $ok=1; }
else { $ok= open(OVR,$overlaps); $ovh= *OVR; }
die "bad -overlaps=$overlaps" unless($ok);

my($noverlist, $overlaplist, $featlist)  = collect_overlaps($ovh); close($ovh);

my $inh= *STDIN;
if($input) {
$ok = ($input =~ /.gz$/) ? open($inh,"gunzip -c $input |") 
      : ($input =~ /^(stdin|-)/) ? $inh= *STDIN
      : open($inh,$input);
die "bad -input=$input" unless($ok);
}

my %joinlist= ();
my %joinscore=();

my( $nin, $nover )= input_gff($inh, $informat);

output_joinlist();

warn "#done: njoin= $nover; ninput= $nin; noverlist=$noverlist\n" if $debug;

#-------------------------------------------------------------------

sub input_gff
{
  my($inh,$informat)= @_;
  my ($nin,$nover, $lastref)= (0) x 10;
  my @chrjoin=();
  
  while(<$inh>){
    next unless(/^\w/);

    my($ref,$src,$typ,$tb,$te,$tp,$to,$tx,$tattr);

      # also read this frm sampairtab.pl == bed variants
      # ("\t",$chr,$mateloc,$cend,$instrand,"mate",$isize)
    if($informat == inBED) {
    ($ref,$tb,$te,$to,$typ,$tattr)= split"\t";
    $to ||= "."; $typ ||=""; $tp=1;
    
    } else {
    ($ref,$src,$typ,$tb,$te,$tp,$to,$tx,$tattr)= split"\t";
    if($passtypes and "$typ.$src" !~ m/$passtypes/) {  next; } # pass other types
    }
    $nin++; # joiner oid, use this as joinID to attach two exons

    if($insorted and $lastref ne $ref) {
      output_joinlist($lastref, \@chrjoin,) if(@chrjoin);
      @chrjoin= ();
    }
    $lastref= $ref;

    my $twidth= 1 + $te - $tb;
    next if($twidth > $MAX_JOIN_SPAN); # use graded criteria? if long require many joiners
    
    my $isjoin=  overlaps( $ref,$tb,$te,$to); # this marks featlist items
    
    if($isjoin) {  
      $nover++;
      my $jid= "J".$nover;
      my $jscore = ($intron2splice) ? $tp : 1;
      $joinscore{$jid}= $jscore;
      
      #? do here or in overlaps() ?
      foreach my $exonid (@$isjoin) {
        my $exloc= $featlist->{$exonid} or next;
        $exloc->[jJOINS] .= "$jid,"; #? or do in caller, but need return rloc
        # introns: need $tp score == number of joiners
        push( @{$joinlist{$jid}}, $exonid); 
        }
        
      push(@chrjoin, $jid) if($insorted);
    }
    
  }
  
  return ($nin,$nover);
}


sub output_joinlist
{
  my($chr, $jidref,) = @_;
  
#   if( $outformat == outGFF ) {
#     output_gff($chr, $jidref);
#   } else {
    output_joinlist1($chr, $jidref);
#   }
}


#? not useful yet
# sub output_gff
# {
#   my($ref,$jidref)= @_;
#  
#   if($ref) {
#     return unless($overlaplist->{$ref});
#     foreach my $ib (sort keys %{$overlaplist->{$ref}}) {
#       my @locs= @{$overlaplist->{$ref}{$ib}};
#       
#       foreach my $rloc1 (@locs) {
#         my $oid= $rloc1->[jOID];
#         my $rloc= delete $featlist->{$oid}; #?? or mark done
#         next unless(ref $rloc);
#        
#         my $oline = $rloc->[jLINE];
#         my $ojoin = $rloc->[jJOINS];  # ** LONG LIST, avoid if possible
#         if($ojoin) { $oline =~ s/$/;join=$ojoin/; }
#         print $oline;
#       }
#       
#     }
#     
#   } else { # dump all? mark done ?
#     for (my $o = 1; $o <= $noverlist; $o++) {
#       my $oid= "N".$o;
#       my $rloc= delete $featlist->{$oid};
#       next unless(ref $rloc);
#       
#       my $oline = $rloc->[jLINE];
#       my $ojoin = $rloc->[jJOINS];  # ** LONG LIST, avoid if possible
#       if($ojoin) { $oline =~ s/$/;join=$ojoin/; }
#       print $oline;
#     }   
#   }
# }


sub output_joinlist1
{
  my($chr, $jidref,)= @_;
  
  unless($jidref and ref($jidref)) {
    my @jid= sort {
        my($an,$bn)=($a,$b); map{s/J//}($an,$bn); $an <=> $bn 
        } keys %joinlist;
    
    $jidref= \@jid;
  }
  
  ## merge jid data: want to see only exons joined by *any* jid, and *count* of joins
  ## output row: exon-join-list(gid,oid)  n-joiners
  
  my %jointexon;
  foreach my $jid ( @$jidref ) 
  {
    next unless(ref($joinlist{$jid}) and  @{$joinlist{$jid}} > 1);
    my $exid = delete $joinlist{$jid}; #?? ok 
    # for introns, jid has score to use below
    my $jscore= $joinscore{$jid} || 1;

    my $nx= @$exid;
    for (my $i=0; $i<$nx; $i++) {
      my $ioid= $exid->[$i];
      for(my $j=$i+1; $j<$nx; $j++) {
        my $joid= $exid->[$j];
        $jointexon{$ioid}{$joid} += $jscore;
        $jointexon{$joid}{$ioid} += $jscore;
      }    
    }
  }
  
  foreach my $ioid (sort keys %jointexon) {
    my $rloc= $featlist->{$ioid} or next; #?
    my ($lr,$lt,$lb,$le,$lid,$lo,$oid,$jgid)= @{$rloc}
      [ jCHR,jTYPE,jBEGIN,jEND,jGID,jSTRAND,jOID,jJOINID ]; # exon span

    unless($jgid) { $rloc->[jJOINID]= $jgid= ++$JGID; }
    
    my $attr="jgid=JG$jgid;gid=$lid;oid=$oid";
  
    my $maxj= 0;
    # my @jex= sort keys %{$jointexon{$ioid}}; #/ sort by count
    my @jex= sort { $jointexon{$ioid}{$b} <=> $jointexon{$ioid}{$a} }
      keys %{$jointexon{$ioid}};
    
    my $nx= @jex; my $ix=0;
    $attr .= ";njoin=$nx;xjoin=";
    foreach my $jex (@jex) {
      my $nj= $jointexon{$ioid}{$jex} or next;
      my $jloc= $featlist->{$jex};
      my ($jref,$jb,$je,$jstrand,$jlid)= @{$jloc}[ jCHR,jBEGIN,jEND,jSTRAND,jGID ];
      if($jlid ne $lid and $jstrand ne $lo) { next; } # cufflinks noise
      $maxj= $nj if($nj > $maxj);

      my $trivial= ($nj < $MIN_NJOINER or $nj/$maxj < $MINOR_NJOIN) ? 1 : 0;
      next if($trivial); #??
      my $separation= _min( abs($jb - $le), abs($je - $lb));
      next if( $jref ne $lr or ($separation > $LARGE_JOIN_SPAN and $nj < $MIN_LONGJOINER));
      
      $jloc->[jJOINID]= $jgid; # unless?
      
      $ix++;
      ## if($jlid ne $lid) { $nj .= "/$jlid"; } # other gene here? problem w/ alt-tr don't want to see all that noise
      if($jlid ne $lid) {
        $nj .= "*"; 
      } # other gene here? problem w/ alt-tr don't want to see all that noise
        # ^^ most from cufflinks.exons are those stupid reversed noise things; drop
        
      $attr .= "$jex:$nj,";
    }


    if($ix < $nx) { $attr =~ s/njoin=$nx/njoin=$ix/; $nx=$ix; }
    print join("\t",$lr,"joined",$lt,$lb,$le,$maxj,$lo,".",$attr)."\n" if($nx>0 and $maxj > 0);
    
  }  
  
}



# convert to array handling? _max(@list) 
sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }

sub overlaps
{
  my( $ref,$tb,$te,$to)= @_;    # == joiner span (intron, mate, ..)

  my ($nover, @oid);
  my $twidth= 1 + $te - $tb;
  my %didid=();
  
  return 0 unless($overlaplist->{$ref});

  my($tb1,$te1, $tb2, $te2)=(0) x 4;
  if($intron2splice) {
    ($tb1,$te1, $tb2, $te2)= 
      ($to eq "-") ? ($te+1,$te + SPLICE,$tb - SPLICE,$tb-1) : ($tb - SPLICE,$tb-1,$te+1,$te + SPLICE); #intron 3bp + 1shift
    $tb= $tb1; $te= $te2;
    
  } else {
    ($tb1,$te1, $tb2, $te2)= 
      ($to eq "-") ? ($te - SPLICEX,$te,$tb,$tb + SPLICEX) : ($tb,$tb + SPLICEX,$te - SPLICEX,$te);  
  }

  my ($ib1, $ib2)= (int($tb/$BINSIZE) , int($te/$BINSIZE));
  for (my $ib = $ib1; $ib <= $ib2; $ib++) 
  {
    $overlaplist->{$ref}{$ib} or next;
    my @locs= @{$overlaplist->{$ref}{$ib}};
    
    foreach my $rloc (@locs) {
      my ($lb,$le,$lid,$lo,$oid)= @{$rloc}[ jBEGIN,jEND,jGID,jSTRAND,jOID ]; # exon span
              
      next if($didid{$oid.$lb.$le}++);
      
      my $over= ($tb <= $le && $te >= $lb) ? 1 : 0;

      if($over) {  # test only if end points hit, not middle of joiner
        my $ok=0; # note: lb,le here are one splice end span of intron: 3 bp?
        my $samestrand= ($stranded and $to =~ /[+-]/ and $lo =~ /[+-]/ and $to ne $lo) ? 0 : 1;
        
        if($tb1 <= $le && $te1 >= $lb) { $ok= ($samestrand) ? 1 : -1; }
        elsif($tb2 <= $le && $te2 >= $lb) { $ok= ($samestrand) ? 2 : -2; }
        else { $ok= 0; } ## ignore inside
        
        unless($ok == 0) {
          push @oid, $oid;  # return exon id
          # do this in caller
          # $rloc->[jJOINS] .= "$jid,"; #? or do in caller, but need return rloc
          # push( @{$joinlist{$jid}}, $oid); 
          
          $nover++;
          }

        }
             
      }
  }
    
  # return $nover;  
  
  if(@oid) { 
    my %oid= map{ $_,1 } @oid; 
    @oid= sort keys %oid;     
    return ((@oid > 1) ? \@oid : 0); # only care about joins of 2+ oid
    }
  return 0;
}


sub collect_overlaps
{
  my($ingff)= @_;  
  my( %overlaps, %rlocs, %feats); # returns
  ## these are exons to join
  
  my $nr=0;
  while(<$ingff>) {   
    next unless(/^\w/); 
    my $inline= $_;
    chomp;
    my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr)= split"\t";  
    $tattr ||="";  
    if($passtypes and "$typ.$src" !~ m/$passtypes/) { next; } # pass other types
    
    $nr++;
    my $oid= "N".$nr; # save always for uniq feature test
    my($gid);  
    if($tattr) {
      if($tattr =~ m/\bID=([^;]+)/) {  $gid=$1; }
      elsif($tattr =~ m/\bParent=([^;]+)/) {  $gid=$1; }
      elsif($tattr =~ m/\bTarget=([^;\s]+)/) { $gid=$1; } # NEW 2010.10
    }
    unless(defined $gid) { $gid = $oid; }


    # ?? handle alt-tr exons here? compress to one exon per loc, append all parent gid, use same oid
    # -- skip exon ID=, use only gid=Parent or Target
    
# use constant{ jBEGIN => 0, jEND => 1, jCHR => 2, jGID => 3, jSTRAND => 4,
#               jTYPE => 5, jOID => 6, jATTR => 7, jLINE => 8, jJOINS => 10 }; ## gff overlap record
              
    my $rloc= [$tb,$te,$ref,$gid,$to,$typ,$oid, $tattr, $inline, 0, 0]; # change to string; save mem

    $rlocs{$oid} = $rloc; # unique hash
    # push( @{$feats{$gid}},$rloc);  # want this or not?

    my ($ib1, $ib2)= (int($tb/$BINSIZE) , int($te/$BINSIZE));
    for(my $ib=$ib1; $ib<=$ib2; $ib++) { push( @{$overlaps{$ref}{$ib}}, $rloc); }
    }
    
  $n_overlaps=$nr; # global
  warn"#collect_overlaps=$nr\n" if $debug;
  return ($nr, \%overlaps, \%rlocs); ## \%feats
}

