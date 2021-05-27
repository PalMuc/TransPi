#!/usr/bin/env perl
# overintron.pl

=item replace old intronscore

  intronscore.pl < $geneset.gff >  $geneset.inscore.tab

  scripts/overlapfilter.perl -intron2splice=error -pass 'exon,intron' -act markid -midtype scoresum \
  -mark intr -over intron/intron_all.gff.gz -in $geneset | scripts/intronscore.pl > $gbase.introntab
  
  upd 2017jan: new outform tags: 
    inid = append intron id to exon, 
    ichain = add intron chain id to mRNA flag
  -FIXINID = "fixed" intron splice id = location "c7b1234567f" = chr7:base1234567forward/+ splice site

=cut

use strict;
use warnings;
use Getopt::Long;

use constant { kINTRON2SPLICE_OVER=>1, kINTRON2SPLICE_QUERY=>2, 
               kINTRONERROR_OVER=> -1, kINTRONERROR_INSIDE => -3}; # only want last 2?
              
use constant { kIntronOver => 0, kIntronInside => 2 };
use constant SPLICE   => 3; # bp for intron splice span tests
use constant SPLICEX  => SPLICE - 1; # for exon

my $debug=1;

my $BINSIZE   = 50000 ; 
our $LONG_INTRON = 19999; # longest accepted w/o valid intron;
my ($overlaps,$passtypes,$dnasequence,@input,$intron2splice,
    $itype,$action,$actid,$typeover,$ok,$mark,$nin);
my $mrnatypes='mRNA';
my $exontypes='exon|CDS';
my $introntypes='intron|ip';
my $annotkey='ints';
my $introns={}; # was $overlaplist= {};
my $stranded=1; # only this opt?
my $outform="";
my $FIXINID=0;
my $doRETAINED_INTRONS=0;

my $optok= GetOptions(
  "introns=s", \$overlaps, 
  "genes=s", \@input,  
  "exontypes=s", \$exontypes, 
  "annotkey=s", \$annotkey, 
  "format=s", \$outform, # upd 2017jan: new outform tags: inid = append intron id to exon, ichain = add intron chain id to mRNA flag
  "passtypes=s", \$passtypes,  #??
  "debug!", \$debug, 
  "FIXINID!", \$FIXINID, # upd 2017jan, make default? -nofix to turn off
  "keptintrons|retainedintrons!", \$doRETAINED_INTRONS, # upd 2017feb, add retained/kept/inner introns to inchain
  "LONGINTRON=i", \$LONG_INTRON, 
  );

die "usage:
  overintron.pl -genes genes.gff  -introns intron_good.gff > genes.annot.gff
" unless($optok and $overlaps); ##  and @input);

if($overlaps) {
my $ovh; 
   if($overlaps =~ /.gz$/) { $ok= open(OVR,"gunzip -c $overlaps |");  $ovh= *OVR; }
elsif($overlaps =~ /^(stdin|-)/) { $ovh= *STDIN; $ok=1; }
else { $ok= open(OVR,$overlaps); $ovh= *OVR; }
die "bad -overlaps=$overlaps" unless($ok);

# $overlaplist= 
$introns= collect_overlaps($ovh); close($ovh);
}

$intron2splice= kINTRONERROR_OVER; # ($overlaps) ? kINTRONERROR_OVER : 0; # only choice?
my $hasintrons= 1; #(ref($introns) and scalar(%$introns))?1:0;
my $testintrons= ($hasintrons and $LONG_INTRON > 0)?1:0;

push @input, @ARGV;
foreach my $input (@input) {
my $inh= *STDIN;
$ok = ($input =~ /.gz$/) ? open($inh,"gunzip -c $input |") 
      : ($input =~ /^(stdin|-)/) ? $inh= *STDIN
      : open($inh,$input);
die "bad -input=$input" unless($ok);

# set output=$input.ovintron.gff ??
my ($nchanged,$ngene)= filter_gff($inh); 
warn"#findcds changed=$nchanged, ngene=$ngene\n" if $debug;
}

#..................

sub filter_gff
{
  my($inh)= @_;
  my ($ng,$nr,$nchange)= (0,0,0);
  my $docomm=($outform =~ /score/i) ? 0 : 1; # else gff
  my $printpass=$docomm;
  # $printpass=0 if($actid == ACT_DROP or $actid == ACT_KEEP);

  my @generec=(); my @otherft=();
  
  while(<$inh>){
    unless(/^\w/){ next if(/^(#n |$)/); print if $docomm; next; }
    my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr)= split"\t";
    $nr++; chomp($tattr);

    if($passtypes and "$typ.$src" !~ m/$passtypes/) { print if $printpass; next; } 
    
    my($gid,$pid); 
    if($tattr =~ m/\bID=([^;]+)/) {  $gid=$1; }
    if($tattr =~ m/\bParent=([^;]+)/) {  $pid=$1; 
      $gid=$pid unless($gid and ($typ =~ /^($mrnatypes)$/)); 
      }
    unless(defined $gid) { $gid = "N".$ng; }

    my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid]; 

    if($typ =~ /^($mrnatypes)$/) {  
      # allow gene and mRNA types .. and/or keep other types in generec
      $nchange += testgene(\@generec, \@otherft) if(@generec);
      
      $ng++;
      @generec = ($rloc); # maybe best store as array of [gffcols], sorted by type-loc
      @otherft=();
      
    } elsif($typ =~ /^($exontypes)$/) {
      push @generec, $rloc;         
    } elsif(/\t/) {
      push @otherft, $rloc;         
    }
      
  }
  
  $nchange += testgene(\@generec, \@otherft) if(@generec);
  
  return ($nchange,$ng);
}

use constant REVFIX17 => 0;

sub testgene
{
  my($generecIN, $otherft)= @_;

  my($mrna)= grep{ $_->[2] eq "mRNA" } @$generecIN;
  my @exons= grep{ $_->[2] eq "exon" } @$generecIN;
  unless(@exons) {  @exons= grep{ $_->[2] eq "CDS" } @$generecIN; }

  unless($mrna and @exons > 0) {
    putgene(undef, $generecIN, $otherft, "");
    return 0;
  }
  
  my $gstrand= $mrna->[6];
  my $nexons = @exons;
  @exons = sort _sort_exons @exons;
if(REVFIX17) {  
  @exons= reverse @exons if($gstrand eq "-"); # NO: 17jan :want this for rev intron chains?
  #^rev : maybe ok if change sub overlaps() to ignore strand ?
}
  
  # my ($inchanged, $infixnote, @exoninfix) = intron_error_cut( \@exongff);
  # my($sovin,$sinerr,$lxe,$lovin)= (0) x 10;
  my $flags="";
  my $iflag="";
  my( $lxe,$lovin, $incode, $insum, $ierrsum, $intotal)= (0) x 10;
  my (@inid,@exid); # upd 17jan
  my $ix=0; 
  foreach my $ex (@exons) {  $ix++;
    my($ref, $xt,$xb,$xe, $to)= ( $ex->[0], $ex->[2], $ex->[3], $ex->[4], $ex->[6]);
    my($ov,$ovin,$inerr,$inid)= overlaps($ix==1, $ix==$nexons, $ref, $xt, $xb,$xe, $to, []);
    push @inid, $inid;
    my $ixw= 1+$xe-$xb; # width 1:99,2:89,.. or relative span, 1:1-99,2:100-189,.. width enough
    push @exid, "x$ix-$ixw";
    if( 1 ) { # $hasintrons
      my $evin = ($ovin & 1) + ($lovin & 2); # should be == 3 for both
      my $evmax= (($ix==1) ? 0 : 2) + (($ix==$nexons) ? 0 : 1);
      
      $insum++ if($ovin & 1); $insum++ if($ovin & 2); # == found splices
      $intotal += (($ix==1) ? 0 : 1) + (($ix==$nexons) ? 0 : 1); # == num splices
      
      $ierrsum += $inerr; # minus mistakes, always -inerr
      my $ints= ($inerr) ? "$inerr/$evin" : $evin;
      if($LONG_INTRON>0 and $lxe>0 and ($xb - $lxe > $LONG_INTRON) and ($evin < $evmax) ) {
        $ints .=",longerr:".($xb - $lxe);
        $iflag .="longerr:".($xb - $lxe).",";
        #?OFF# $ierrsum += -1; # NOT -3 here, longerr is only possible err
      }
      unless($ints eq "0") {
        $ints.=",$inid" if($outform =~ /inid/);
        ## $iflag .= "i$ix:$ints,"; # $evin/$inerr,"; 
        $ex->[8] =~ s/$/;$annotkey=$ints/;
      }
    }

    $lxe=$xe; $lovin= $ovin;
  }
  
  if($intotal>0) {
    ## bad score  < -100 h for longerr == -3; redo output format:
    ## inqual=INCODE,insum/intotal,err:$ierrsum
    $incode= int (100 * ($insum + $ierrsum) / $intotal); # can be -
    #old# $flags .= "$annotkey=$incode," . (($ierrsum) ? "$ierrsum/$insum/$intotal" : "$insum/$intotal");
    $flags .= "$annotkey=$incode,$insum/$intotal";
    $flags .= ",inerr:$ierrsum" if($ierrsum); # drop iflag list of introns?
    $flags .= ",$iflag" if($iflag); # drop iflag list of introns?
    if($outform=~/chain/){ 
      my $ichain=join",",@inid; 
      my $xchain=join",",@exid;
      #? add exon ix chain?: xc:1/w1,2/w2,..n/wn
      $flags .= ";ic:$ichain;xc:$xchain"; 
      }
  } else {
    $incode= 1;
  }
  
  # $mrna->[8] =~ s,$,;$flags, if($flags);
  putgene( $mrna, $generecIN, $otherft, $flags);
  return ($insum,$ierrsum,$intotal); #? 
}

sub putgene {
  my($mrna, $exons, $otherft, $flag)= @_;
  my $cskip= ""; #($flag && $flag =~ /skip=/) ? "#x.": "";
  $otherft ||= [];
  # my $svec= $gref->[jSCOREVEC];
  # $flag .= ";svec=$svec" if($svec);
  
  # if($outform =~ /scoreonly/)  print only mrna ID, intron flag
  my %didx=();
  foreach my $ex ($mrna, @$exons, @$otherft) { # $gref, 
    next unless($ex and ref $ex); 
    next if($didx{$ex}++); # for mrna in exons
    my($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid)= @$ex;
    if($outform =~ /score/i) {
      print join("\t",$gid,$flag),"\n"; return;
    }
    $tattr =~ s/;$//; 
    $tattr .=";$flag" if($flag); $flag=""; # only for gref
    print join("\t", $cskip.$ref, $src, $typ, $tb, $te, $tscore, $to, $tph, $tattr),"\n";
  }
}

sub _sort_exons
{
  # ($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid)
  # fixed: ?? with ref 1st this isn't working; big score not at front
  return ($a->[0] cmp $b->[0]) # ref
      || ($a->[2] cmp $b->[2]) # type: CDS > exon
      || ($a->[3] <=> $b->[3]) # begin
      || ($b->[4] <=> $a->[4]) # end; small>large ?
      ;
}


sub overlaps {
  my($xfirst, $xlast, $ref, $typ, $tb, $te, $to) = @_; # , $locs

  my ($overin,$inerr, $inid)=(0,0,"");
  my $dointrons= ($hasintrons and not($xfirst and $xlast)) ? 1 : 0;
  my $isrev=($to eq "-")?1:0;
if(REVFIX17) {  
  $isrev=0; # input reverse order exons
}      
  # always check tb,te for intron supt? dont do for each overlap test.
  if($dointrons) {
    my($ov1, $ov2, $do1, $do2, $tb1,$te1, $tb2, $te2)= 
      ($isrev) ? # dont need? REVFIX17
        ( 2, 1, !$xlast, !$xfirst, $te - SPLICEX,$te, $tb,$tb + SPLICEX) 
      : ( 1, 2, !$xfirst, !$xlast, $tb,$tb + SPLICEX,$te - SPLICEX,$te);  
    
    if($do1) { 
      my ($ovi,$iid)=overintron($ref, $tb1, $te1, $to, 0); 
      $inid="0";
      if($ovi<0) { $inerr += $ovi; } elsif($ovi>0) { $inid=$iid; $overin+= $ov1; }
      # $overin += $ov1 if($ovi>0); $inerr += $ovi if($ovi<0); 
    }
    if($do2) { 
      my ($ovi,$iid)=overintron($ref, $tb2, $te2, $to, 0); 
      my $idd="0";
      if($ovi<0) { $inerr += $ovi; } elsif($ovi>0) { $idd=$iid; $overin+= $ov2; }
      $inid= ($do1 and $isrev)?"$idd,$inid" : ($do1)?"$inid,$idd":$idd;    # upd 17jan: if($rev) change order ??
      #o.$inid=($do1)?"$inid,$idd":$idd;      # upd 17jan: if($rev) change order ??
      # $overin += $ov2 if($ovi>0);  $inerr += $ovi if($ovi<0); 
    }
  }

    ## test exon contains intron: kINTRONERROR_INSIDE
    ## 1702 upd: option to add these to inid as "retained_intron" locations
  if($hasintrons) {
    my ($ovi,@iid)= overintron( $ref, $tb + 2*SPLICE, $te - 2*SPLICE, $to, kIntronInside); 
    if($ovi == kINTRONERROR_INSIDE) { 
      $inerr += -1;   # NOT -3 here
      if($doRETAINED_INTRONS) {
        my $iid=join",", map{ $_."k" } @iid;  # need useful intronID flag k=kept, "r", "i" bad
        $inid= ($inid)? "$inid,$iid" : $iid;
      }
    }
  }

  return (0,$overin,$inerr, $inid);
}

sub fixinid {  # use in collect_overlaps() stead of overintron() 
  my($ref,$lb,$lo)=@_;
  if($ref=~s/^Chr/c/i) {}
  elsif($ref=~s/^(Scaffold|Supercontig|Super)[_]?/s/i) {}
  elsif($ref=~s/^Contig/ct/i) {}
  else { $ref="c$ref"; }
  my $lid= $ref . "i$lb" . (($lo eq "-" )?"r":"f");# or $lo < 0
  return $lid;
}  

sub overintron {
  my($ref, $tb, $te, $to, $testtype) = @_;  # tb,te == exon splice point here
  return 0 unless($introns->{$ref});
  my $ok= 0; my @inid=();
  my ($ib1, $ib2)= (int($tb/$BINSIZE) , int($te/$BINSIZE));
  for (my $ib = $ib1; $ib <= $ib2; $ib++) {
    $introns->{$ref}{$ib} or next;
    my @locs= @{$introns->{$ref}{$ib}};
    foreach my $rloc (@locs) {
      #data: my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid,$inb,$ine]; 
      my ($lb,$le,$lo,$lid,$inb,$ine)= @{$rloc}[3,4,6,9,10,11];
      # lb,le are splice ends of intron; inb,ine is full span
      
      ## upd 17jan, lid to location??: $lid from collect_overlaps is N1234 or gff ID=xxx
      ## moved to collect_overlaps (but there changes inid!): 
      # $lid= fixinid($ref,$lb,$lo) if($FIXINID);
      
      my $samestrand= ($to =~ /[+-]/ and $lo =~ /[+-]/ and $to ne $lo) ? 0 : 1;  

      # wait, dont return 1st hit, if error, check for valid and clear error if valid found
      if($testtype == kIntronInside) { 
        if($ok==0 or $ok == kINTRONERROR_INSIDE) { 
          if($tb <= $inb and $te >= $ine and $samestrand) {
            $ok= kINTRONERROR_INSIDE;
            push @inid, $lid; #* fixme double counting some in splice sites
            # c2i154500rk-c2i154611rk and c2i154500rk-c2i154677rk < same start ID
            # only want splice-pair IDs, two inIDs with same inb-ine span, and longest span if 2+ share splice
            } 
          else {
            # any other? ignore intron_inside revstrand ?
          }
        }
      } else { 
        if($tb <= $le && $te >= $lb) { if($ok<=0) { $ok= ($samestrand) ? 1 : -1; } }
      }
    if($ok>0) { return ($ok,$lid); }
    }
  }
  if($testtype == kIntronInside and $ok == kINTRONERROR_INSIDE) { return($ok,@inid); }
  return $ok;
}


sub collect_overlaps  # collect_introns
{
  my($gff)= @_;
  my ($nr,$nx)=(0,0);
  my %overlaps=();
  while(<$gff>){
    next unless(/^\w/); chomp;
    my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr)= split"\t";
    $tattr ||="";

    # ?? only Introns here? need splice-site calcs ; no ID= for introns...
    next unless($typ =~ /^($introntypes)$/);

    $nr++;
    my($gid,$pid,$lid); 
    if($tattr =~ m/\bID=([^;]+)/) {  $gid=$1; }
    if($tattr =~ m/\bParent=([^;]+)/) {  $pid=$1; 
      $gid=$pid unless($gid and ($typ =~ /^($mrnatypes)$/)); 
      }
    unless(defined $gid) { $gid = "N".$nr; }
    
    #? only kINTRONERROR_OVER and kINTRONERROR_INSIDE
    my($inb,$ine)= ($tb,$te); # full intron span, keep
    if(1) { ## $intron2splice == kINTRON2SPLICE_OVER or $intron2splice == kINTRONERROR_OVER) # 2010jul
      my($s1b,$s1e,$s2b,$s2e)= 
        ($to eq "-") ? ($te+1,$te + SPLICE,$tb - SPLICE,$tb-1) : ($tb - SPLICE,$tb-1,$te+1,$te + SPLICE); # 3bp + 1 shift
      ($tb,$te)= ($s1b,$s1e);
      
      $lid= ($FIXINID) ? fixinid($ref,(($to eq "-")?$tb:$te),$to) : $gid; # NOTE this is splice id, not intron(2splice) id
      my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$lid,$inb,$ine]; 
      my @bins= (int($tb/$BINSIZE) .. int($te/$BINSIZE));
      foreach my $ib (@bins) { push( @{$overlaps{$ref}{$ib}}, $rloc); }
      ($tb,$te)= ($s2b,$s2e);  #? change gid, oid?
    }  
    
    $lid= ($FIXINID) ? fixinid($ref,(($to eq "-")?$te:$tb),$to) : $gid; # NOTE this is splice id, not intron(2splice) id
    my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$lid,$inb,$ine]; 
    my @bins= (int($tb/$BINSIZE) .. int($te/$BINSIZE));
    foreach my $ib (@bins) { push( @{$overlaps{$ref}{$ib}}, $rloc); } # $generec
  }
  
  warn"#collect_overlaps n=$nr\n" if $debug;
  return \%overlaps;
}

__END__

# old parts

# sub scoregene {
#   my($gid,$genes,$exons,$introns,$otherft)= @_; # ,$flags
#   # $flags="CDS,exon" unless($flags);
#   
#   my $gref= $genes->{$gid} or return;
#     
#   #?? also score introns, add to mrna score as per overbestgene2? see above
#   ## iscore as weight: (valid splices - error splices) / total splices ??
#   my $flags="";
#   my( $incode, $insum, $ierrsum, $intotal)= (0) x 10;
#   if($hasintrons) {
#     my $iflag=""; 
#     my($nover, $cdsover, $xbad, $ix, $lxe, $lovin) = (0) x 10;
#     my $rexons= $exons->{$gid};
#     my @exons= grep { $_->[2] ne "CDS" } @$rexons; # CDS > exon
#     @exons= @$rexons unless(@exons);
#     @exons = sort _sort_exons @exons;
#     
#     my $nexons= @exons;  
#     $ix=0; foreach my $ex (@exons) {  $ix++;
#       my($ref, $xt,$xb,$xe, $to)= ( $ex->[0], $ex->[2], $ex->[3], $ex->[4], $ex->[6]);
#       my($ov,$ovin,$inerr)= overlaps($ix==1, $ix==$nexons, $ref, $xt, $xb,$xe, $to, []);
#       if( 1 ) { # $hasintrons
#         my $evin= ($ovin & 1) + ($lovin & 2); # should be == 3 for both
#         
#        $insum++ if($ovin & 1); $insum++ if($ovin & 2); # == found splices
#         $intotal += (($ix==1) ? 0 : 1) + (($ix==$nexons) ? 0 : 1); # == num splices
#         
#         $ierrsum += $inerr; # minus mistakes
#         my $ints= ($inerr) ? "$inerr/$evin" : $evin;
#         if($LONG_INTRON>0 and $lxe>0 and ($xb - $lxe > $LONG_INTRON) and ($evin < 3) ) {
#           # $xbad++; 
#           $ints .=",longerr:".($xb - $lxe);
#         }
#         unless($ints eq "0") {
#         $iflag .= "i$ix:$ints,"; # $evin/$inerr,"; 
#         $ex->[8] =~ s/$/;ints=$ints/;
#         }
#       }
# 
#       $lxe=$xe; $lovin= $ovin;
#     }
#     if($intotal>0) {
#       $incode= int (100 * ($insum + $ierrsum) / $intotal); # can be -
#       $flags .= "ints=$incode," . (($ierrsum) ? "$ierrsum/$insum/$intotal" : "$insum/$intotal") .",$iflag;";
#     } else {
#       $incode= 1;
#     }
#   }
# 
# }

# sub intronoverlaps
# {
#   my ( $ref,$tb,$te,$to)= @_;    # exon here
#   my ( %didid,@overs,@ovok);
#   return 0 unless($introns->{$ref});
#   
#   my($tb1,$te1, $tb2, $te2)=(0) x 4;
#   if($intron2splice == kINTRONERROR_OVER) {
#     # input == exon, over= intron, test if overlap is splice end or internal
#     ($tb1,$te1, $tb2, $te2)= 
#       ($to eq "-") ? ($te - SPLICEX,$te,$tb,$tb + SPLICEX) : ($tb,$tb + SPLICEX,$te - SPLICEX,$te);  
#   }
#   
#   my($linb,$line)= (0,0);
#   my ($ib1, $ib2)= (int($tb/$BINSIZE) , int($te/$BINSIZE));
#   for (my $ib = $ib1; $ib <= $ib2; $ib++) 
#   {
#     $introns->{$ref}{$ib} or next;
#     my @locs= @{$introns->{$ref}{$ib}};
#     foreach my $rloc (@locs) {
#       #new# my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid,$inb,$ine]; 
#       
#       my ($lb,$le,$lo,$oid,$inb,$ine)= @{$rloc}[3,4,6,9,10,11];  # intron here 
#        # FIXME now is intron ENDs
#        # FIXME: need full intron lb,le and also splice ends
#        
#       next if($didid{$oid.$lb.$le}++);
#       
#       my $over  = ($tb <= $le && $te >= $lb) ? 1 : 0;      
#       # $over=0 if($stranded and ($to =~ /[+-]/ and $lo =~ /[+-]/ and $to ne $lo)); 
#       # my $inside= ($over and $tb <= $lb && $te >= $le) ? 1 : 0; ## intron inside  
#       # $inside=0 if($stranded and ($to =~ /[+-]/ and $lo =~ /[+-]/ and $to ne $lo)); 
#  
#       # if($inside) {  push @overs, [$lb,$le]; } else
#       
#       if($over and $intron2splice == kINTRONERROR_OVER) { # test for splice reversed errors
#        
#         my $ok=0; # note: lb,le here are one splice end span of intron: 3 bp?
#         my $samestrand= ($to =~ /[+-]/ and $lo =~ /[+-]/ and $to ne $lo) ? 0 : 1;
#         
#         if($tb1 <= $le && $te1 >= $lb) { $ok= ($samestrand) ? 1 : kINTRONERROR_OVER; } # $errt="rev" unless $samestrand
#         elsif($tb2 <= $le && $te2 >= $lb) { $ok= ($samestrand) ? 2 : kINTRONERROR_OVER; }
#         #old# elsif(($tb + 2*SPLICE <= $lb && $te - 2*SPLICE >= $le) and $samestrand) { $ok= kINTRONERROR_INSIDE; } 
#         elsif(($tb + 2*SPLICE <= $inb && $te - 2*SPLICE >= $ine) and $samestrand) { $ok= kINTRONERROR_INSIDE; } 
#                   # ^^ same strand inside  == retained intron err
#         # else what?
#         
#         if($ok < 0) {
#           if($ok == kINTRONERROR_INSIDE) { push @overs, [$inb,$ine, $ok]; } # need error type also: retained vs strand-err
#           else { push @overs, [$lb,$le, $ok]; } # need error type also: retained vs strand-err
#           #? need ok1, ok2 for both ends of exon to say if 1 end is ok?
#         } elsif( $ok > 0) {
#           push @ovok, [$lb,$le, $ok]; #? return which end of exon is supported by intron (or both)
#           # push @spliceok, ($ok == 2) ? (($to eq "-") ? $tb : $te) : (($to eq "-") ? $te : $tb);
#         }
#         # next;
#       } 
#         
#       # push @overs, [$lb,$le] if ($over);
#       }
#   }
# 
# # ** FIXME: this makes bad cuts, likely where 2+ alt-introns are overlapped
# 
#   my @opens=();
#   if(@overs) {
#     my $okover= scalar(@ovok);
#     my($bb,$be)= ($tb,$te); my $errlast= 0; my($llb,$lle)=(0,0);
#     @overs= sort _sort_over @overs;
#     foreach my $ab (@overs) {
#       my($lb,$le,$errcode)= @$ab;
#       ## .. for $errcode == kINTRONERROR_OVER, strand err at exon splice; need to chop off entire exon ?
#       # next if($lb >= $llb and $le <= $lle); # skip inside alt-intron
#       next if($lb < $lle and $le > $llb); # skip overany alt intron
#       
#       ## .. for $errcode == kINTRONERROR_INSIDE
#       if($le < $bb) {  }
#       elsif($lb <= $bb && $le > $bb) { $bb= $le+1; }
#       elsif($lb < $be) {  #  && $le > $be
#         my($b1,$e1)= ($bb,$lb-1);
#         push @opens, [$b1,$e1,$errcode,$okover];
#         $bb= $le+1; 
#         } 
#       elsif($lb > $te) { last; } #?
#       ($llb,$lle)= ($lb,$le);
#       $errlast= $errcode;
#       last if($bb >= $te);
#     }
#     if($bb < $te) { push @opens, [$bb,$te,$errlast,$okover]; } # add end point
#     return \@opens;
#     # return (\@opens, \@spliceok);
#   } else {
#     return 0; ## [[$tb,$te]];
#     # return ( [], \@spliceok);
#   }
# }

# sub intronadd {
#   my($introns, $rloc)= @_;
#   my($ref,$src,$typ,$tb,$te,$score,$to,$tph,$tattr,$gid)= @$rloc;
#   
#   my($s1b,$s1e,$s2b,$s2e)= 
#     ($to eq "-") ? ($te+1,$te + SPLICE,$tb - SPLICE,$tb-1) 
#     : ($tb - SPLICE,$tb-1,$te+1,$te + SPLICE); # 3bp + 1 shift
# 
#   ($tb,$te)= ($s1b,$s1e);
#   my $rloc1= [$ref,$src,$typ,$tb,$te,$score,$to,$tph,$tattr,$gid]; 
#   for( my $ib= int($tb/$BINSIZE); $ib <= int($te/$BINSIZE); $ib++) { 
#     push( @{$introns->{$ref}{$ib}}, $rloc1); }
# 
#   ($tb,$te)= ($s2b,$s2e);   
#   my $rloc2= [$ref,$src,$typ,$tb,$te,$score,$to,$tph,$tattr,$gid]; 
#   for( my $ib= int($tb/$BINSIZE); $ib <= int($te/$BINSIZE); $ib++) { 
#     push( @{$introns->{$ref}{$ib}}, $rloc2); }
# }


