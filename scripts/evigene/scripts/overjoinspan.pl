#!/usr/bin/env perl

=item overjoinspan.pl 
    from overlapjoins.pl

  - measure gene/exon span covered by mates (and dangling mates)
  
  - prefered result would be gene exons w/ mate cover/uncover spans indicated
    (including holes in exon cover that may be introns), plus some
    end-of-gene measure of dangling mates: how much of gene model is missing?
    
  usage:
  
  gzcat genes/XXX.gff.gz |  \
  $evigene/scripts/overjoinspan.pl -over stdin -in rnas/cleaned/allpe_matec.tab \
  -format bed -strand -sorted > genes/XXX.matespan.gff
  
  -- revise for sam input
  samtools view  -f2 $aphid2/rnas/bams/aphidpe_SRR098330.bam Scaffold25 | ...
  
  ** this works, maybe what is needed (if score is correct :)
  -- needs path with samtools and options (now -f2 = only perfect pairs)
  -- should parallelize with chrlist subsets; use velvet chr subsets?
  
  cat *.sc25.gff | $evigene/scripts/overjoinspan.pl -over stdin -chrlist Scaffold25 \
  -format sam -strand $aphid2/rnas/bams/aphidpe_*.bam > compare3.matespans.sc25.gff
  
=item memory overloads

  - when run parallel 30 parts x 10 chrs each using chrlist, 5 bam inputs, samtools
  - possibly due to single chr/scaffolds w/ lots of pe reads
  - somewhere input_joins is eating up mem, more than should, 
    only persistant addition I see is in markjoin() add to gene/exon annots.
  - savememory opt doesnt solve this: process 1 chr per samtools call, dump out/delete genes after that chr  
  
=cut


use strict;
use Getopt::Long;

use constant DOT => '.';

my $debug=1;
our $BINSIZE  = 500 ; #was# 5000;

my $samtoolsview= $ENV{samview} || "samtools view -f 2"; # proper pairs only? option

## these should be options
my $MAX_JOIN_SPAN = 500000;  # we can have longer introns; want reliability though
my $LARGE_JOIN_SPAN = 20000; # require more joiners to count?

use constant SPLICE   => 3; # bp for intron splice span tests
use constant SPLICEX  => 24; #was SPLICE - 1; # for exon x mate ends FIXME: for sam/bam use readsize

use constant{ jBEGIN => 0, jEND => 1, jCHR => 2, jGID => 3, jSTRAND => 4,
              jTYPE => 5, jOID => 6, jATTR => 7, jSRC => 8, jLINE => 9, 
              jJOINS => 10, jJOINERR => 11, jJOINID => 12, }; ## gff overlap record
              
use constant{ inGFF => 0, inBED => 1, inSAM => 2, outJOINS => 0, outGFF => 1 };

my $mrnaType="mRNA";
my $exonType="exon";
my ($overlaps, $passtypes, $chrlist ) = ("") x 10;
my ($ok, $savememory, $stranded, $insorted, $informat, $outformat, $intron2splice, $n_overlaps)= (0) x 20;
my $JGID= 0;
my @input= ();

my $optok= GetOptions(
  "input=s", \@input,       # required: exon joiner spans, eg. mate pairs, introns, in format bed,sam/bam,gff
  "overlaps=s", \$overlaps, # required: gene gff to annotate and output
  "format=s", \$informat, 
  #"outformat=s", \$outformat, # not used
  "chrlist=s", \$chrlist, # subset of chrs, for .bam input only?
  "samtoolsview=s", \$samtoolsview, # for .bam input only
  "sorted!", \$insorted,     # if input sorted by chr
  "mrnaType=s", \$mrnaType,  # gene mRNA types to annotate w/ exon join scores
  "exonType=s", \$exonType,  # gene exon types to join
  "stranded!", \$stranded,   # use stranded input tests : default?
  "memorysave!", \$savememory,   
  "debug!", \$debug, 
  );

die "opt error" unless($optok and $overlaps);

push @input, @ARGV;

$mrnaType =~ s/[,]/\|/g;
$exonType =~ s/[,]/\|/g;
##$intron2splice=1 if($passtypes =~ /intron/);
$intron2splice=1 if($informat =~ /intron/);

unless($informat) { 
  if($input[0] =~ /\.(bam|sam)/) { $informat= inSAM; }
  elsif($input[0] =~ /\.(bed|mate)/) { $informat= inBED; }
  elsif($input[0] =~ /\.(gff)/) { $informat= inGFF; }
}

if($informat =~ /bed|mate/) { $informat= inBED; }
elsif($informat =~ /sam|bam/) { $informat= inSAM; }
elsif($informat =~ /gff/) { $informat= inGFF; }
else { $informat= inGFF; }

# if($outformat =~ /gff/) { $outformat= outGFF; } # not used
# else { $outformat= outJOINS; }

$chrlist ||= "";
my %chrlist= ($chrlist) ? map{ $_,1 } split( /[,;|\s]+/, $chrlist) : ();
# maybe using too much mem with chrlist sets:
#  31 processes x 20 scaffolds each  on 64GB 32 core machine using all mem
# put in loop here per chr in list? then need to re-read gff, reads/chr
# or only for samtools, loop input_joins x bams per chr ?
# only expanding mem per read is per exon jJOINS,jJOINERR strings : compress those in loop?

my $ovh;  
   if($overlaps =~ /.gz$/) { $ok= open(OVR,"gunzip -c $overlaps |");  $ovh= *OVR; }
elsif($overlaps =~ /^(stdin|-)/) { $ovh= *STDIN; $ok=1; }
else { $ok= open(OVR,$overlaps); $ovh= *OVR; }
die "bad -overlaps=$overlaps" unless($ok);

# %$exonlocs is unused, drop

my($noverlist, $overlaplist, undef, $generecs)  = collect_overlaps($ovh); 
close($ovh);

my %havegene=();
# my %joinlist= ();
# my %joinscore=();
my( $nin, $nover )= (0,0);

$chrlist =~ s/[,;|\s]+/ /; # SPACE for samtools
my @chrs= split " ",$chrlist;
if($savememory and @chrs > 2 and $informat == inSAM) {

} elsif($chrlist) {
  @chrs=($chrlist); # may be empty; use dummy?
} else { 
  @chrs=("dummy"); #?
}

foreach my $onechr (@chrs) {
  $chrlist= ($onechr eq "dummy") ? "" : $onechr;
  %chrlist= ($chrlist) ? map{ $_,1 } split( /[,;|\s]+/, $chrlist) : ();
  
foreach my $input (@input) {
  # allow multiple @input; # but this presumes not insorted / -nosort
  # allow .bam using samtools view $bam[$i] $chrlist
  
  my $inh= undef; ## *STDIN;
  if($input =~ /\.bam$/) { $informat= inSAM; $ok= open($inh,"$samtoolsview $input $chrlist |"); }
  else { $ok = ($input =~ /.gz$/) ? open($inh,"gunzip -c $input |") 
        : ($input =~ /^(stdin|-)/) ? $inh= *STDIN : open($inh,$input); 
  }
  
  warn "#input.$informat: $input\n" if $debug;
  die "bad -input=$input" unless($ok);
  
  my( $nin1, $nover1 )= input_joins($inh, $informat);
  close($inh);  
  $nin += $nin1; $nover += $nover1;
  warn "#input.$informat:  nin=$nin1; nover=$nover1\n" if $debug;
}

## deal w/ generecs output for loop
my $outchr= ($onechr and $onechr !~ /\s/)? $onechr : "";
output_joinlist($outchr);  # this does mark/delete generecs as output .. so shouldnt do twice

warn "#done: njoin= $nover; ninput= $nin; noverlist=$noverlist\n" if $debug;
}

#-------------------------------------------------------------------

sub input_joins
{
  my($inh,$informat)= @_;
  my ($nin,$nover, $lastref, $lastb)= (0) x 10;
use constant USE_TLIST => 1;
  my @tlist;
  my $issplice= 0;
  
  while(<$inh>) {
    next unless(/^\w/);

    my($ref,$src,$typ,$tb,$te,$tp,$to,$tx,$tattr);
    my @span;
    
    if($informat == inBED) {
      ($ref,$tb,$te,$to,$typ,$tattr)= split"\t";
      $to ||= "."; $typ ||=""; $tp=1;
      @span=( $to, $tb, $te); $issplice=0;
      
    } elsif($informat == inSAM) {
      my($flg,$cref,$cb,$cigar,$mc,$mb,$pairspan)=(split"\t")[1,2,3,5,6,7,8]; 
      
      # update here to use cigar N introns, 
      # use spliced reads == mates and also count splice site aligns
      # also, use read length instead of fixed SPLICEX
      # check $flg flags?
      
      if($mc ne "=") { next; } # bad pair
      
      $ref= $cref;  
      $to= (/XS:A:([+-])/) ? $1 : ".";
      $typ="read"; $tp=1;
      # 36M76N18M .. parse cigar for introns that change total span, set strand
      my $ce= $cb; my $rlen=0;
      
      ## add intron splices to @span **
      @span= ($to, $cb );  $issplice=1;
      while( $cigar =~ m/(\d+)([A-Z])/g ) { 
        my($n,$t)=($1,$2);
        if($t=~/^[MNS]/) {
          $rlen+= $n if($t=~/^[MS]/); 
          my $cb1= $ce + $n;
          if($t=~/^[N]/) { push( @span, $ce, $cb1); }
          $ce = $cb1; 
          }
        }
      
      push( @span, $ce);

      my $me= $mb + $rlen; # allow missing pair, if have intron splice
      push( @span, $mb,$me) if($mb>0); 
      ##@span= ($to, $cb, $ce, $mb,$me);

      # large spans here may be problem: instead use quartet: $cb,$ce,$mb,$me
      if($mb < $cb) { $tb= $mb; $te= $ce; } else { $tb= $cb; $te= $me; }
      
    } else {
      ($ref,$src,$typ,$tb,$te,$tp,$to,$tx,$tattr)= split"\t";
      @span=( $to, $tb, $te); $issplice=0;
    }
    
    next if($chrlist and not $chrlist{$ref});    
    $nin++; # joiner oid, use this as joinID to attach two exons

    if($insorted and $lastref and $lastref ne $ref and %havegene) {
      my @chrjoin = sort keys %havegene; %havegene=();
      output_joinlist($lastref, \@chrjoin,) if(@chrjoin);
    }

    my $twidth= 1 + $te - $tb;
    next if($twidth > $MAX_JOIN_SPAN); # use graded criteria? if long require many joiners

    my $isjoin;
if(USE_TLIST) {    
    #? bundle some input_joins for efficient overlaps() ?
    push( @tlist, \@span); ##  [$tb,$te,$to] 
    unless($lastref eq $ref and $lastb + $BINSIZE > $tb) {  ## oops dont update lastb here, memory overload!
      $isjoin=  overlaps( $ref, $issplice, @tlist); # fixed
      @tlist=(); $lastb= $tb;
    }
} else {    
    $isjoin=  overlaps( $ref, $issplice, \@span); ## [$tb,$te,$to]  # this marks exonlocs items
}
    $lastref= $ref; 

   if($isjoin) {  # drop here, do $havegene{$gid}++; in overlaps()
      $nover += (USE_TLIST) ? $isjoin : 1; # count from @tlist
   }   
    
  }

  if(@tlist) {
    my $isjoin= overlaps( $lastref, $issplice, @tlist);
    $nover += (USE_TLIST) ? $isjoin : 1; # count from @tlist
  }

  if($insorted) {
    my @chrjoin = sort keys %havegene; %havegene=();
    output_joinlist($lastref, \@chrjoin,) if(@chrjoin);
  }
  
  return ($nin,$nover);
}




sub output_joinlist
{
  my($chr, $jidref,) = @_;
  output_genes($chr,$jidref); 
#    output_joinlist1($chr, $jidref);
}





# convert to array handling? _max(@list) 
sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }

sub overlaps
{
  my( $ref, $issplice, @tbeolist)= @_;    # == joiner span (intron, mate, ..)
  #old#my( $ref,$tb,$te,$to)= @_;    # == joiner span (intron, mate, ..)

  return 0 unless($overlaplist->{$ref});
  
  my ($nover, @oid);
  foreach my $tbeo (@tbeolist) 
  {
    # my($tb1,$te1, $tb2, $te2);
    my($tb,$te, $to)= (0) x 10;
    my @tspan=();
    
    if($issplice) {   
      ($to,@tspan)= @$tbeo;
    } else {  
      ($to,$tb,$te)= @$tbeo;
      @tspan= ($tb,$tb + SPLICEX,$te - SPLICEX,$te); # short cut for now ; do in caller?
    }
    
#     my $dorev= 0; ## NOT ($to eq "-");
#     if($intron2splice) {
#       ($tb1,$te1, $tb2, $te2)= 
#         ($dorev) ? ($te+1,$te + SPLICE,$tb - SPLICE,$tb-1) : ($tb - SPLICE,$tb-1,$te+1,$te + SPLICE); #intron 3bp + 1shift
#       $tb= $tb1; $te= $te2;
#       
#     } else {
#       ($tb1,$te1, $tb2, $te2)= 
#         ($dorev) ? ($te - SPLICEX,$te,$tb,$tb + SPLICEX) : ($tb,$tb + SPLICEX,$te - SPLICEX,$te);  
#     }
    
    my $checkstrand= ($stranded and $to =~ /[+-]/) ? 1 : 0;
    
    ## my @tl=($tb1,$te1,$tb2,$te2);  # handle any number of pairs tb,te here 
    # oops, now @tspan can have 3+ pairs from introns + mates.
    # end = 1,2,... 
    my $nend= int(@tspan / 2);
    foreach my $end (1..$nend) {
      my %didid=();
      my ($ob,$oe)= splice( @tspan,0,2); 
      next unless($oe>0);
      my ($ib1, $ib2)= (int($ob/$BINSIZE) , int($oe/$BINSIZE));
      for (my $ib = $ib1; $ib <= $ib2; $ib++) {
        $overlaplist->{$ref}{$ib} or next;
        my @locs= @{$overlaplist->{$ref}{$ib}};
        foreach my $rloc (@locs) {
          my ($lb,$le,$lid,$lo,$oid)= @{$rloc}[ jBEGIN,jEND,jGID,jSTRAND,jOID ]; # exon span
          next if($didid{$oid.$lb.$le}++); ## is this bad? need hit same exon 2 times: tb1,tb2
  
          my $over= ($ob <= $le && $oe >= $lb) ? 1 : 0;
          if($over) {  # test only if end points hit, not middle of joiner
            my $samestrand= ($checkstrand and $lo =~ /[+-]/ and $to ne $lo) ? 0 : 1;
            
            my $ok= ($samestrand) ? $end : -$end; # end=1..3 if intron, 1..2 for 2 mates only
              
            my $span="$ob-$oe"; 
            my $retval= [$lid, $ok, $oid, $span];
            push @oid, $retval;  # return exon id
            $nover++;
            }
          }
        }
      }
  }  
  
  if(@oid) { 
    @oid = sort { $a->[0] cmp $b->[0] or $a->[1] <=> $b->[1] or $a->[2] cmp $b->[2] } @oid;
    
    # reduce result here?  good if @oid == 2 and ok == 1,2 and lid1 == lid2
    #  bad if (@oid == 1 or ok < 0 or lid1 != lid2 but check thru all for lid1==lid2
    my ($qual,$gid1,$oid1,$lgid,$lov,$ov1)= (0) x 10;
    ## my @oid2;    
    foreach my $ov (@oid) {
      my($gid,$ok,$oid,$span)= @$ov;

      if($lgid and $gid ne $lgid) { 
        if($lov and $qual<2) { ##push(@oid2, $lov); 
          markjoin( $lgid, $qual, $lov->[2], $lov->[3]); 
        }
        $gid1=$qual=0; 
      }
      
      if($ok == 1) { $qual=1;  $gid1= $gid; $oid1=$oid; $ov1=$ov; }
      elsif($ok >= 2 and $gid eq $gid1) { 
        $qual=($oid eq $oid1) ? 2 : 3; ## push(@oid2, $ov1, $ov);
        markjoin( $gid, $qual, $ov1->[2], $ov1->[3], $oid, $span);
        $gid1= $gid; $oid1=$oid; $ov1=$ov; # case of >2 splice ends
        } 
      elsif($ok < 0 and $qual < 2) { $qual=$ok; }
      $lgid= $gid; $lov= $ov;
    }
    if($lov and $qual<2) {
      ## push(@oid2, $lov); 
      markjoin( $lgid, $qual, $lov->[2], $lov->[3]);
    }  

    return scalar(@oid);
    # return \@oid2; # dont need oid2 no
    }
    
  return 0;
}


sub markjoin {
  my($gid, $qual, $xid1, $xspan1, $xid2, $xspan2)= @_;

  my $gr= $generecs->{$gid} or return; ## is this error? warn?
  #my $mrna = $gr->[0]; # check type
  $havegene{$gid}++ if($insorted); # from input_ DO ONLY if $insorted ; otherwise eats memory for no use
  
  my $xspan1q= "$xspan1/$qual,";
  my $xspan2q= ($xspan2)?"$xspan2/$qual,":"";
  
  foreach my $rloc (@$gr) {
    my ($ltyp,$lb,$le,$lid,$lo,$oid)= @{$rloc}[ jTYPE,jBEGIN,jEND,jGID,jSTRAND,jOID ]; # exon span
    
    my $ismrna= ($ltyp =~ m/$mrnaType/)?1:0;
    
    ## more efficent store for many small spans ?  
    if($qual >= 2) {  # qual==3 means cross-exon, same gene : save this
      my $js= $rloc->[jJOINS] || ""; my $up=0;
      if($ismrna or $oid eq $xid1) { $js .=  $xspan1q; $up++; }
      if($xspan2 and ($ismrna or $oid eq $xid2)) { $js .= $xspan2q;  $up++; }
      if($up and $savememory and length($js) > 2000) { $js= spanof($js, 1); }
      $rloc->[jJOINS]= $js if($up);
      ## condense spans at some point .. here?
      
    } else {  # error spans
      if($ismrna or $oid eq $xid1 or $oid eq $xid2) { # or both?
      my $js= $rloc->[jJOINERR] || ""; my $up=0;
      $js .=  $xspan1q; $up++;
      if($xspan2) { $js .=  $xspan2q; $up++; }
      if($up and $savememory and length($js) > 2000) { $js= spanof($js, 1); }
      $rloc->[jJOINERR]= $js if($up);
      }
     
    }

  }

}



sub spanof {
  my($joins,$inexon)= @_;
  $joins =~ s/^0+//; # init bug
  unless($joins =~ /\d/) { return  wantarray ? (0,0,0,0,"") : ""; } 
  my @sp=();
  foreach my $sp (split",",$joins) {
    my $err= ($sp =~ s,/(.*),,) ? $1 : 0;  # err now == qual for all
    my($sb,$se)= split(/[\.\-]+/,$sp);
    push(@sp, [$sb,$se,$err]) if($se>0);    
  }
  
  my $ns= scalar(@sp); 
  my($max,$tot,$mb,$me,$lb,$le)=(0) x 10;
  @sp= sort { $a->[0] <=> $b->[0] } @sp;
  for(my $i=$ns-1; $i >= 0; $i--) {
    my($sb,$se,$err)= @{$sp[$i]}; 
    
    my $ends= ($inexon) ? $se+50 : $se;
    my $ov= ($le > $sb and $lb <= $ends)? 1 : 0;
    if($ov) { 
    # collapsing errs here? but called for separate jerr, joins (err=0)
      $se= $le;  $sb= $lb if($lb < $sb);
      splice(@sp,$i+1,1); $sp[$i]= [$sb,$se,$err]; 
      }
    $lb=$sb; $le=$se;
  }
  
  $ns= scalar(@sp); my $span="";
  for(my $i=0; $i < $ns; $i++) {
    my($sb,$se,$err)= @{$sp[$i]};
    my $w= 1+$se-$sb; $tot+= $w;
    $span .= ($err) ? "$sb-$se/$err," : "$sb-$se,"; 
    if($w>$max) { $max=$w; $mb=$sb; $me=$se; }
  }
  
 return wantarray ? ($tot,$max,$mb,$me,$span) : $span;
}


sub output_genes 
{
  my($chr, $jidref,)= @_;
  # @$jidref == now is list of gene ids, use this
  my $nout=0;
  my @geneids= (ref $jidref) ? @$jidref : (sort keys %$generecs);
  foreach my $gid (@geneids) {  
    my $gr= $generecs->{$gid} or next;
    next if($gr->[0]->[jSRC] eq "done" or ($chr and $gr->[0]->[jCHR] ne $chr));
    # delete $generecs->{$gid}; # delete or mark as output ??
    # $gr->[0]->[jSRC]= "done"; # oops, do After output
    $nout++;
    
    ## should build mrna score/join attr from exons here
    foreach my $rloc (@$gr) {
      my ($lr,$lsrc,$lt,$lb,$le,$lid,$lo,$oid,$attr,$joins,$jerrs)= @{$rloc}
      [ jCHR,jSRC,jTYPE,jBEGIN,jEND,jGID,jSTRAND,jOID,jATTR,jJOINS,jJOINERR ]; # exon span
      next unless($le>0); # bugs ??
      my $score= 0;
      my $ismrna= ($lt =~ m/$mrnaType/)?1:0;
    
      my($jtot,$jmax,$jb,$je,$jlist) = spanof($joins, 1); # !$ismrna);
      my($etot,$emax,$eb,$ee,$elist) = spanof($jerrs, 1); # !$ismrna);
      $score= $jtot - $etot; #?
      $attr= ($ismrna) ? "ID=$lid" : "Parent=$lid" if($debug);
      $attr .= ";joins=$jtot/$jlist" if($joins);
      $attr .= ";jerrs=$etot/$elist" if($jerrs);
      
      print join("\t",$lr,$lsrc,$lt,$lb,$le,$score,$lo,".",$attr)."\n";
      # if($savememory) { for my $i (0..jJOINID) { $rloc->[$i]=0; } }
    }
  $gr->[0]->[jSRC]= "done"; # oops, do After output
  }
  
  #?? $| = 1; # unbuffer
  print "#done.out chr=$chr nout=$nout\n";
  #?? $| = 0; # bad??
}


sub collect_overlaps
{
  my($ingff)= @_;  
  my( %overlaps, %rlocs, %feats); # returns
  ## these are exons to join
  # rlocs == %$exonlocs is unused, drop

  #** FIXME: add to/create global $chrlist and %chrlist from ingff
  my %inchrlist=();
  
  my $nr=0;
  while(<$ingff>) {   
    next unless(/^\w/); 
    my $inline= $_;
    chomp;
    my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr)= split"\t";  
    $tattr ||="";  
    #? if($passtypes and "$typ.$src" !~ m/$passtypes/) { next; } # pass other types
    next if($chrlist and not $chrlist{$ref});
    
    $nr++;
    my $oid= "N".$nr; # save always for uniq feature test
    my($gid);  
    if($tattr) {
      if($tattr =~ m/\bID=([^;]+)/) {  $gid=$1; }
      elsif($tattr =~ m/\bParent=([^;]+)/) {  $gid=$1; }
      elsif($tattr =~ m/\bTarget=([^;\s]+)/) { $gid=$1; } # NEW 2010.10
    }
    unless(defined $gid) { $gid = $oid; }
    $inchrlist{$ref}++;

# use constant{ jBEGIN => 0, jEND => 1, jCHR => 2, jGID => 3, jSTRAND => 4,
#               jTYPE => 5, jOID => 6, jATTR => 7, jSRC => 8, jLINE => 9, 
#               jJOINS => 10, jJOINERR => 11, jJOINID => 12, }; ## gff overlap record
              
    my $rloc= [$tb,$te,$ref,$gid,$to,$typ, $oid, $tattr, $src, $inline, 0, 0, 0];    

    # handle mRNA,exon separately?  use mRNA only for exon summary annot?
    if($typ =~ /($mrnaType)/) {
      push( @{$feats{$gid}},$rloc);  

    } elsif($typ =~ /($exonType)/) {
      # $rlocs{$oid} = $rloc; # unique hash
      push( @{$feats{$gid}},$rloc);  

      my ($ib1, $ib2)= (int($tb/$BINSIZE) , int($te/$BINSIZE));
      for(my $ib=$ib1; $ib<=$ib2; $ib++) { push( @{$overlaps{$ref}{$ib}}, $rloc); }

    } else { # save other parts for output ?
    
    }
  }
    
  ## %rlocs unused now ?  
  # reset global chrlist
  %chrlist= %inchrlist;
  $chrlist= join " ", sort keys %inchrlist; # SPACE sep for samtools
  
  $n_overlaps=$nr; # global
  warn"#collect_overlaps=$nr\n" if $debug;
  return ($nr, \%overlaps, \%rlocs, \%feats); ## 
}



__END__

=item trbamspan.pl

# trbamspan.pl: quic subset of trbamstats.sh  
# ** FIXME3: try measuring largest span covered by mated pairs, not count, using perfect matches?
#
$debug=1;
my($outf,$bams,$hd);  
my $gr=$ENV{gr} || 0;
if($gr == 1) { $outf="treg8pspan.tab2";  $bams="*-ap2pub8tr.bam"; }   
elsif($gr == 2) {  $outf="trncbi2pspan.tab2"; $bams="*-a2ncbireftr.bam"; }
elsif($gr == 3) { $outf="tracypi1pspan.tab2"; $bams="*-ap1acypitr.bam"; }
else { die "which group? env gr=1,2,3\n"; }
warn "# out=$outf, bin=$bams\n" if $debug; 

my @bams= `/bin/ls -1 $bams`;
foreach my $bm (@bams) { 
  chomp($bm);  warn "# samtools view -f 1 $bm \n" if $debug;
  my($np1,$pp1,$pno1,$px1,$poor1)=(0) x 10;
  open(S,"samtools view -f 1 $bm |"); 
  while(<S>){ 
    my($c,$cb,$cg,$mc,$mb,$pairspan)=(split"\t")[2,3,5,6,7,8]; $np1++; $np{$c}++; 
    unless($cg =~ /^\d+M$/) { $pim{$c}++; $poor1++; } 
    elsif($mc eq "=") { $pp{$c}++; $pp1++; 
      if($pairspan<0) { span($c, $mb, $mb - $pairspan); } 
      elsif($pairspan>0) { span($c, $cb, $cb + $pairspan); }
    } 
    elsif($mc eq "*") { $pno{$c}++; $pno1++;} else { $px{$c}++; $px1++; } 
  } close(S); 
  
  my $log= "# np=$np1, pp=$pp1, pno=$pno1, px=$px1, imperf=$poor1 in $bm\n";
  if($debug) { warn $log; $_= $log; 
   my @v=m/\w+=(\d+)/g; my $n=$v[0]||1; my @p=map{ int(0.5+100*$_/$n); }@v; 
   my $o=$_; for my $i (0..$#v){ $o=~s/=$v[$i]/% $p[$i]/; } warn $o; 
   }
}  

open(OUTH, ">$outf");
foreach my $c (sort keys %np) { 
  my @v= map{ $_||0 } ($np{$c},$pp{$c},$pno{$c},$px{$c},$pim{$c}); 
  my ($spanx,$spanb,$spane)= maxspan($c);
  print OUTH join("\t",qw(trID pspan spanb spane npair pperf pnone pother imperf)),"\n" unless($hd++);
  print OUTH join("\t",$c,$spanx,$spanb,$spane,@v),"\n"; 
}  close(OUTH);

sub maxspan { # save also maxspan b,e
  my($c)=@_;  
  my @sp= (ref $span{$c}) ? @{$span{$c}} : (); 
  my $ns= scalar(@sp); 
  my($max,$mb,$me,$lb,$le)=(0) x 10;
  @sp= sort { $$a[0] <=> $$b[0] } @sp;
  for(my $i=$ns-1; $i > 0; $i--) {
    my($sb,$se)= @{$sp[$i]}; 
    if($le > $sb and $lb <= $se) {
      $se= $le;  $sb= $lb if($lb < $sb);
      splice(@sp,$i+1,1); $sp[$i]= [$sb,$se]; 
      }
    $lb=$sb; $le=$se;
  }
  for(my $i=0; $i < $ns; $i++) {
    my($sb,$se)= @{$sp[$i]}; my $w= 1+$se-$sb; 
    do { $max=$w; $mb=$sb; $me=$se; } if($w>$max);
  }
 return (wantarray) ? ($max,$mb,$me) : $max;
}

sub span { 
  my($c,$cb,$ce)= @_; 
  my @sp= (ref $span{$c}) ? @{$span{$c}} : (); 
  my $over=0;  my $ns= scalar(@sp);
  for(my $i=0; $i < $ns; $i++) {
    my($sb,$se)= @{$sp[$i]};
    my $ov= ($cb <= $se and $ce >= $sb)?1:0;
    if($ov) { $sb=$cb if($cb < $sb); $se=$ce if($ce > $se); $sp[$i]= [$sb,$se]; $over++; }
  }
  unless($over) { push(@sp, [$cb,$ce]); } #? @sp= sort { $$a[0] <=> $$b[0] } @sp; 
  $span{$c}= \@sp;
}

=cut

=item parallel joinspans

  ** this works, maybe what is needed (if score is correct :)
  -- needs path with samtools and options (now -f2 = only perfect pairs)
  -- should parallelize with chrlist subsets; use velvet chr subsets?
  
  cat *.sc25.gff | $evigene/scripts/overjoinspan.pl -over stdin -chrlist Scaffold25 \
  -format sam -strand $aphid2/rnas/bams/aphidpe_*.bam > compare3.matespans.sc25.gff

grep -c exon a*.gff
acyr1_acypi.gff:154602
acyr2_ncbirefgene.gff:139137
aphid2_evigene8f.gff:236287

  cat genes/a*.gff > genes/compare3_genes.gff
  
subset lists: 42 parts
    12 rnas/velmapt7/bigsub.list
     >> 30 parts for one trestles node (32 cpu)
    10 rnas/velmapt7/midasub.list
    10 rnas/velmapt7/midbsub.list
    10 rnas/velmapt7/smallsub.list

midasub.list:
subset  mida0   Scaffold60,Scaffold8,Scaffold114,Scaffold136,Scaffold53,Scaffold222,Scaffold50,Scaffold149,Scaffold158,Scaffold92,Scaffold229,Scaffold193,Scaffold502,Scaffold115,Scaffold257,Scaffold221,Scaffold90,Scaffold74
subset  mida1   Scaffold12,Scaffold78,Scaffold225,Scaffold37,Scaffold195,Scaffold148,Scaffold70,Scaffold38,Scaffold279,Scaffold272,Scaffold164,Scaffold182,Scaffold56,Scaffold256,Scaffold65,Scaffold33,Scaffold67,Scaffold244

#.....................
## memory overloads w/ this even after splitting mid,small scaf lists
##  64 GB mem eaten up in 15 min, less than 1/2 thru, -memorysave no cure
## fixes: 1. merge bamset to bam/allpe.bam, 2. use -sorted -in bam/allpe.bam 
#! /bin/bash
### env rnain=xxxxx qsub -q normal genernaspan.sh
#PBS -N genernaspan
#PBS -A ddp138
#PBS -l nodes=1:ppn=32,walltime=3:55:00
#PBS -o genernaspan1.$$.out
#PBS -e genernaspan1.$$.err
#PBS -V

## avoid HOME, use  /scratch/$USER/$PBS_JOBID   /oasis/$USER; drop bash -l login flag
ncpu=32
scratch=/oasis/$USER

export PATH=$scratch/bio/bin:${PATH}

workd=$scratch/chrs/aphid2
bindir=$scratch/bio/bin

# all 3 gene sets at one go: evigene2f,ncbiref2,acypi1 : 500,000 exons
genes=$workd/genes/compare3_genes.gff

cd $workd/rnas/

#OLD# bamset=bams/aphidpe_*.bam
bamset=bams/allpe.bam
# bams/aphidpe_SRR071347.bam  bams/aphidpe_SRR075803.bam	bams/aphidpe_SRR098330.bam
# bams/aphidpe_SRR075802.bam  bams/aphidpe_SRR097896.bam


# bigs: 12
# nam=bigs; subsets=(`cut -f3 sublists/bigsub.list`)
#OLD# nam=mids; subsets=(`cut -f3 sublists/{mid,small}*sub.list`)
# mids: 30, valid for i=0 .. 29
# recut subsets to 21 parts, run on 1 node (32 cores) to save mem.

i=0; while [ $i -lt $ncpu ]; do {

  ichrset=${subsets[$i]}
  if [ $ichrset ] ; then
   ( ./overjoinspan.pl -strand -over $genes  -format sam -sorted -in $bamset -chrlist $ichrset \
      > compare3_genes.pespan.$nam$i.gff ) &
  fi
  
  i=$(( $i + 1 ))
}
done
wait

#..................

  
=cut
