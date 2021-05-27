#!/usr/bin/perl
# overbestgene1.perl

=item overbestgene1
  
  special case of gene gff overlap detection:
   - input one gene/cds gff file with quality scores
      and many overlapping predictions (e.g. exonerate predicts)
   - output best by score non-overlapping subset genes + cds
     (some overlap possible, want all mostly distinct predictions)
     
 gzcat exonerate-gldmoj.gff3.gz | $td/overbestgene1.perl -in stdin > exonerate-gldmoj-best.gff &
 #collect_gff=1660943, ngene=763282   + 6735 geneids with terepeat > about same...
 #done kept=23000, skipped=744001
    vs kept=15038, skipped=748244 for simple filter (includes terepeat filter) 
    
=item  eg

gzcat $wasp/intron/intron_good.gff.gz | cat - nasvit1-uparp_hym_exonr5.gff | $evigene/scripts/overb
estgene1.perl -longint=29999 -major=2 -typeover mostover -pct 6 -strand -score 'alignx:20,intron:20,CDS:2' -in stdi
n > nasvit1-uparp_hym_exonr5.best9o.gff

gzcat $workd/intron/intron_good.gff.gz nvit_epi6c1-augmap.gff.gz | $evigene/scripts/overbestgene1.p
erl -typeover over -pct 10 -strand -score 'intron:20,CDS:1,UTR:2' -in stdin -skip | grep '      mRNA' > nvit_epi6c1
-augmap.ovbest1.mrna

see $caca/prot/exonrs/over.info
  alignx = exonerate align score
  using -majorityvote here to pick best from models with 2+species agreement, ie skip outlier/bad proteins.
  
gzcat $caca/intron/intron_good.gff.gz cacao11all-*-best4.exonr4.gff.gz | $evigene/scripts/overbestg
ene1.perl -major=2 -typeover over -pct 6 -strand -score 'alignx:30,intron:30,CDS:1' -in stdin > cacao11-plant8-mbes
t7.gff

    
=cut

use strict;
use warnings;
use Getopt::Long;

use constant VERSION => '2014.05.28'; # handle Split= genes;
# "2012.08.03"; # dropscores, .. 
# "2011.09.11";  "2012.07.20";# added alttr processing; fixups

use constant SAMEBASE => 3; # for _sameloc, slop allowed in loca == locb
# use constant { ACT_DROP=>1, ACT_KEEP=>2, ACT_MARK=>3, ACT_MARK_WITH_ID=>4 };
use constant { kOVERLAP=>1, kSAMELOC=>2, kNEARLOC=>3, kINSIDE => 4 }; ## kBESTGENE=>4, kBESTEXONS=>5, };

use constant SPLICE   => 3; # bp for intron splice span tests
use constant SPLICEX  => SPLICE - 1; # for exon
# modified gff gene [rloc] record fields
use constant{ jCHR => 0, jSRC => 1, jTYPE => 2, jBEGIN => 3, jEND => 4, 
              jSCORE => 5, jSTRAND => 6, jPHASE => 7, jATTR => 8, jGID => 9, 
              jCDSLEN => 10, jTRLEN => 11, jSCOREVEC => 12,  }; 
              # last 3 for mRNA only; jSCOREVEC replaces jSCORESUM 
use constant USE_DROPSCORE_ANDOR => 1;
use constant { kDROP_OR => 0, kDROP_AND => 1, kDROP_NOT => -1, kDROP_MUSTNOT => 2 }; 
use constant {  # from overbestgene2.perl
  kSKIP_OVERLAP => 1, 
  kSKIP_TRIVIAL_EVD => 2, # later skim out paralog-only
  kSKIP_BAD => 3, # if flagged in attrib : not same as skip required ?
  kDEMOTE_MAYBEJOIN => 56, # 
  kKEEP_REQUIRED => 77,
  kKEEP_ALTERNATE => 69,
  kSKIP_REQUIRED => 666,   # from user input list
  kSKIP_NOEVD => 999,  # fails to pass dropscore
  };
  

our $BINSIZE  = 50000 ; #was# 1000; want large steps for this case? bin edge problems
our $LONG_INTRON = 19999; # longest accepted w/o valid intron;
our $PMAJORSAME  = 0.75; #? 0.80; #?? 0.85;  base agreement to call to models similar enough for majority vote
my  $PCDS_OK  = 50; # min pct CDS/exon ratio for valid mRNA, pcode UTR

# our $NEARDIST =  500; # was 15k; needs to be < BINSIZE

our $debug=1;
my ($showskips,%geneexons,$input,$itype,$action,$actid,$ok,$mark);

my $exontypes="CDS,exon,match_part,HSP";
my $mrnatypes="mRNA,match"; # fixme ?? need this
my $introntypes = "intron";

my $SPECIALFIELD='aalen|utrx|inqual'; # several scores need parsing..

# my $add_codingscore= 1; # option
my $badattr="";
my $typeover="overlap";
my $stranded=1; # should be =1 default?
my $pctover= 0.10; # was 0; default = 0.10 ??
my $doalttr=0;
my $domajority= 0;
my @saveargs= @ARGV;
my ($scoretype, $dropscore, @scoreweight, @scorefield, @dropscore, @droptype);
$scoretype= "score"; $dropscore="";


my $optok= GetOptions(
  "input=s", \$input,  
  "typeover=s", \$typeover, 
  "stranded!", \$stranded, 
  "exontypes=s", \$exontypes, 
  "badattr=s", \$badattr, 
  "mrnatypes=s", \$mrnatypes, 
  "pctover=i", \$pctover, 
  "BINSIZE=i", \$BINSIZE, 
  "LONGINTRON=i", \$LONG_INTRON, 
  "majorityvote:1", \$domajority, #? default :i==1
  "scoretype=s", \$scoretype,  
  "dropscore=s", \$dropscore,  
   #^ need multiscores now like overbestgene2: xxx,CDS,introns,
   
  "alttr!", \$doalttr, # 2012.07 add; change this to typeover == altover 
  "skipshow!", \$showskips, 
  "debug!", \$debug, 
  );

die "
usage:  perl overbestgene1 -input genes.gff|stdin  > bestgenes.gff
  -exontypes $exontypes  -mrnatypes $mrnatypes
  -pctover 90 -scoretype aalen -[no]stranded -typeover over|inside|samefeat
" unless($optok);
# -act keep|mark|markid -mark=bestgene

if($typeover =~ /alt/) { $doalttr=1; }  # which???
elsif($doalttr) { $typeover= "alt".$typeover; }

$pctover= $pctover/100.0 if($pctover > 1);
$exontypes =~ s/[,; ]+/|/g;
$mrnatypes =~ s/[,; ]+/|/g;

my $manyscore= ($scoretype =~ /,\w/)?1:0; # commas mean  many
getScoretypes();

my $inh= *STDIN;
$ok = ($input =~ /.gz$/) ? open($inh,"gunzip -c $input |") 
      : ($input =~ /^(stdin|-)/) ? $inh= *STDIN
      : open($inh,$input);
die "bad -input=$input" unless($ok);

my %keepgene=(); # print
my %skipgene=(); # optionally list ..
my($ngenein, $nskip, $nkeep)= (0) x 10;

sub MAIN{}
# MAIN:
my ($genebins, $genes, $exons, $introns, $otherft) 
      = collect_gff($inh); 

my $hasintrons= (ref($introns) and scalar(%$introns))?1:0;
## my $testintrons= ($hasintrons and $LONG_INTRON > 0)?1:0;

  # calc cxlen scores here ? cds/exon total, .. or only if scoretype=CDS,UTR,trlen ?
# if($add_codingscore or $hasintrons)  # always now
foreach my $gid (keys %$genes) { 
  $ngenein++; 
  scoregene($gid, $genes, $exons, $introns, $otherft); 
  # enabled dropscore= ; sets $skipgene{$gid} for any scored below drops
  }

# print gffheader,  version, scores # FIXME
HEADER: {
  my $version= VERSION;
  print <<"EOGFF";
##gff-version 3
# appl: $0
# args: @saveargs
# vers: $version
# scores: '$scoretype'
# drops : '$dropscore'
# ingene: $ngenein

EOGFF
}

# memsave: if -sorted can collect_gff and process per bin w/o reading all gff 1st

foreach my $gref (sort keys %$genebins) {
  my @bins= sort{$a <=> $b} keys %{$genebins->{$gref}};
  my %exonpatt=(); # for doalttr
  foreach my $ib (@bins) {  
    my %keepgbin=();
    my %major=();
    my @didexons=();
    #not here?# my %exonpatt=(); # for doalttr
    
    # as overbestgene2, this sort puts high score ahead of low, then throws out lower, overlapped models per locus
    # proper score is critical to picking best
    # .. modify _sort for manyscore_sum() ?
    
    my @sgenes = sort _sort_refscoreloc @ { $genebins->{$gref}{$ib} };
      
    foreach my $gref (@sgenes) {
      my($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid)= @$gref;
      # FIXME: Split= gene
      
      next if( $skipgene{$gid} );
      
      my $rexons= $exons->{$gid} ;
      unless(ref $rexons) {
        # are we missing exon-gene id links?
        $skipgene{$gid} = "error:missing exons";
        next;
      }
      
      if($keepgene{$gid}) { push(@didexons, @$rexons); next; }  #? never here

    #...... replacement sub ................................
#         my ($nover, $nexons, $xbover, $xbtotal, $ncds, $cdsover, $xbad, \%overid)
#           = overanygene( $gid, \@didexons);
    #.......................................................
  
        my($nover, $xbover, $xbtotal, $cdsover, $cdsbover, $isum, $xbad, $ix, $lxe, $lovin) = (0) x 10;
        my %overid=();
        my $flags="";
        
        #my @exons = sort _sort_exons grep { $_->[2] ne "CDS" } @$rexons; # CDS > exon
        my @exons= grep { $_->[2] ne "CDS" } @$rexons; # CDS > exon
        @exons= @$rexons unless(@exons);
        @exons = sort _sort_exons @exons;
        my $nexons= @exons;  
        my $exonpatt="";
        
  # FIXME: for doalttr, keep if score high and exons have diff splice pattern (not containedin others)
  # OPTION: too many 1exon fragments, both strands overlapping good genes; figure option to remove some by score
        
        my $savestrand= $stranded; # 1exon fix for overlaps()
        #at cdsexon# $stranded=0 if($nexons==1); #? or not, unless cdsover
        
        $ix=0; foreach my $ex (@exons) {  $ix++;
          my($xt,$xb,$xe)= ( $ex->[2], $ex->[3], $ex->[4]);
          my($ov,$ovin,$inerr,@overid)= 
            overlaps($ix==1, $ix==$nexons, $ref, $xt, $xb,$xe, $to, \@didexons, "allover");

          if($ov) { $nover++; map{ $overid{$_}++ } @overid; }; # for major; isalt bug? need all overid, not one
          $xbover += $ov;  $xbtotal+= 1 + $xe - $xb;
          $xbad++  if($badattr && $ex->[8] =~ m/$badattr/);
          $exonpatt.= "$xb-$xe,";
          
          ## always add ovin score as annot to exon/mrna
          if( $hasintrons ) {
            my $evin= ($ovin & 1) + ($lovin & 2); # should be == 3 for both
            $isum += $evin;
            if($LONG_INTRON>0 and $lxe>0 and ($xb - $lxe > $LONG_INTRON) and ($evin < 3) ) {
              $xbad++; 
            }
          }
          $lxe=$xe; $lovin= $ovin;
        }
        $exonpatt=~s/,$//; 
        $exonpatt{$gid}= $exonpatt; # doalttr
        $exonpatt{$exonpatt}= $gid unless( $exonpatt{$exonpatt}); # test bug; dont need now?
        
        $stranded=0 if($nexons==1); #? or not, unless cdsover
        @exons =  sort _sort_exons grep { $_->[2] eq "CDS" } @$rexons; # CDS > exon
        my $ncds= @exons;  
        $ix=0; foreach my $ex (@exons) {  $ix++;
          my($ov,$ovin,$inerr,$overid)= 
            overlaps($ix==1, $ix==$ncds, $ref, $ex->[2], $ex->[3], $ex->[4], $to, \@didexons);
          $cdsover++ if($ov);
          $cdsbover += $ov;  # $cdstotal += 1+$xe-$xb;
        }
        
        
        ## FIXME: below, here:  nover vs cdsover : need to allow UTR-overlaps but not CDS-over
        # also fix 1exon unstranded cdsover only .. 
        # want to allow cases of 1exon CDS in UTR of other
        #  --CCCC..CCCC...----->
        #                 <-CCCC----  valid 1exon in  utr of other.
        #  but also brings in trivials:
        #   CCCC--->  trivial 1exon non-cdsover cases to avoid; nearly same tr w/ orfs called at each end
        #   <---CCCC 
        
        if($nexons <= 1) { $nover= $cdsover; $nexons= $ncds; $xbover=$cdsbover; }
        $stranded= $savestrand; # 1exon fix
        my $utronlyover= ($cdsover == 0 and $ncds > 0 and $nover > 0) ? 1 : 0;
        
  #.............overanygene..................................
       
# FIXME updates:
# x- score for introns, NOT inside if matches alt-intron 
# ** score cxlen and avoid aberrant UTR models, at least keep inside ones if valid c/x ratio
#    include both CDS,exon types here? distinguish cds-over from exon-over?
# -- proteins: add major vote on best model: 2+ near same; use _similargene( gid1, gid2) and score by pctident
# .. need to count base overlap spans, major == near same bases; need overgenedup similar exon scoring?
# .. mark alttr here? for same-locus, lower scoring but not allinside
# .. ditto for asmrna ?

        my $majorsame= 0; my $majoroverid="";
        if($domajority>0 and $nover>0) {
          my @ovid = sort{ $overid{$b} <=> $overid{$a} } keys %overid; # all should be keepgene
          # use all w/ ovid same score? YES, 1st can be bad model
          my $ovmax= 0;
          foreach my $ovid (@ovid) {
            my $ovscore= $overid{$ovid};
            last if($majorsame != 0 and $ovscore < $ovmax);
            $ovmax= $ovscore if($ovscore> $ovmax);
            if(my $ovref = $genes->{$ovid}) {
              my $glen = $gref->[jTRLEN]  || $gref->[jCDSLEN];  
              my $ovlen= $ovref->[jTRLEN] || $ovref->[jCDSLEN];  
              my $ovscore= ($glen>0 or $ovlen>0) ? $xbover / _max($glen, $ovlen) : 0;
              $majorsame= (($ovscore >= $PMAJORSAME) 
                       and ($nover >= $PMAJORSAME * $nexons)) ? 1 : 0; 
              $majoroverid= $ovid;
              }
  
            if($majorsame > 0) { 
              push(@{$major{$ovid}}, $gid); 
              $majorsame= scalar(@{$major{$ovid}}); 
              # last; # @ovid
            } else {
              if( $major{$ovid} and @{$major{$ovid}} >= $domajority) { $majorsame= -1; } # last;  # can skip this
            }
           
          }
          
        } elsif($nover>0) { # not majority but info
          my @ovid = sort{ $overid{$b} <=> $overid{$a} } keys %overid; # all should be keepgene
          $majoroverid= $ovid[0]; # use to flag skips w/ better gene
        }  
        
        
  # FIXME: here? for doalttr, keep if score high and exons have diff splice pattern (not containedin others)
        my $isalt=0;
        if($nover > 0 and $nexons > 1 and $cdsover>0 and $doalttr) { # and ! $utronlyover or $cdsover>0
          my $clen= $gref->[jCDSLEN];
          my $xlen= $gref->[jTRLEN];
          my $pcode= ($clen>0 and $xlen>0) ? int(100*$clen/$xlen) : 100;   
          my $iexonpatt=  $exonpatt{$gid};
          $isalt= ($pcode >= $PCDS_OK and $iexonpatt) ? 1 : 0;  #  or $cdsover == $ncds
          $isalt=0 if($isalt and $exonpatt{$exonpatt} ne $gid); # got some identicals; why miss below??

          if($isalt) {
          my $iexonpatt0= $iexonpatt;  
          $iexonpatt =~ s/^\d+//; $iexonpatt =~ s/\d+$//; # trim endpoints : 1Exon problem? dont allow 1exon alts
          foreach my $ovid (sort keys %overid) {
            #?? BUG HERE missing some subset alts.. fix above, @overids for allover
            my $ovpatt=  $exonpatt{$ovid} or next; # next if($skipgene{$ovid}); #?? or next;
            if( index($ovpatt,$iexonpatt) >=0 ) { $isalt=0; last; } # subset model
            else { 
              $ovpatt=~s/^\d+//; $ovpatt=~s/\d+$//;  # test both ways? 
              # maybe NOT GOOD IDEA, shorter may clobber longer, lower score but still valid? but what else?
              if( index($iexonpatt0,$ovpatt) >=0 ) { $isalt=0; last; }
              }
            }

          # annotate if isalt
          if($isalt) { $flags.="isalt=1;"; }  ## $nover=0 ; 
          }
        } 
        elsif($nover > 0 and $typeover =~ /all|most/) { ##  and $typeover =~ /inside/ # allow also non-inside here          
          my $clen= $gref->[jCDSLEN];
          my $xlen= $gref->[jTRLEN];
          # my($clen,$xlen)= $tattr =~ m/cxlen=(\d+).(\d+)/;
          my $pcode= ($clen>0 and $xlen>0) ? int(100*$clen/$xlen) : 100;   
          
          if( $pcode < $PCDS_OK and $cdsover == $ncds) { } # skip; add pcode flag?
          elsif($typeover =~ /all/) { $nover=0 if($nover < $nexons); }
          elsif($typeover =~ /most/) { $nover=0 if($nover < 0.75 * $nexons); } # >= 3/4 
        }
        
        $flags.="nover=$nover;cdsover=$cdsover;";
        $flags.="overid=$majoroverid;" if($majoroverid and not $domajority);
        
        if($xbad > 0) { 
          $flags.="xbad=$xbad;"; # nover=$nover;cdsover=$cdsover;";
          $skipgene{$gid} = $flags;
          $nskip++;  

        ## not majorsame here for majority: overgene can be wrong model (join)
        } elsif($domajority>0 and $nover>0) { # need to save but may not output
        
          # if($majorsame > 0 or $majorsame < 0) 
          unless($majorsame == 0) { 
            # 1 case where 1+ better models of same general size exist, 
            # -1 case of smaller/poor but overgene has majority vote
            $flags .="majorid=$majoroverid;"; # nover=$nover;cdsover=$cdsover;";
            $skipgene{$gid} = $flags;
            $nskip++;  
            
          } else { 
            # no better same-size gene, keep this one for now?
            # BUT check if overgene has majority vote, then skip this
            $flags .="majorid=$majoroverid;"; # nover=$nover;cdsover=$cdsover;";
            # no: $keepgene{$gid}++;  $nkeep++; #??
            push(@didexons, @$rexons);
            $keepgbin{$gid}= $flags; #?
          }

#         } elsif($nover > 0) { 
#           $skipgene{$gid} = $flags;
#           $nskip++;  
          
        } elsif($nover < 1 or $isalt or $utronlyover) {  
          push(@didexons, @$rexons);
          
          if($domajority>0) { $keepgbin{$gid}= $flags; } 
          else { 
            $keepgene{$gid}++; $nkeep++;
            putgene( $gref, $rexons, $otherft->{$gid}, $flags);   
            }
        } else {
          $skipgene{$gid} = $flags; ## add kSKIP codes ?
          $nskip++;  
        }
          
      }  # sgenes loop
      
    putgbin(\%keepgbin, \%major) if($domajority>0 and %keepgbin);  
    # ^^ update , \%keepgene, \%skipgene, nskip/nkeep
    }
}


if($showskips) {
  print "# skipped genes ",(".") x 50," \n";
  foreach my $gid (sort keys %skipgene) {  
    my $sflag= $skipgene{$gid};
    putgene( $genes->{$gid}, undef, undef,"skip=1;$sflag");   
    }
}

warn"#done kept=$nkeep, skipped=$nskip\n" if $debug;

#.......................................

sub putgbin {
  my($gbin, $major)= @_;  # , $rkeepgene, $rskipgene
  return unless(ref($gbin) and %$gbin);
  
  my %didkeep=();
  my %hasmajority= ();
  if(ref $major) {
    map{ my $ra= $major->{$_}; my $n=scalar(@$ra); $hasmajority{$_}= $n if($n >= $domajority); } 
       (keys %$major);
  }
  
  # need to print 2 steps: 
  #  1. print if hasmajority, 
  #  2. print others if no overlap to major gene. >> need to recheck overlaps?
  # ** need score sort first
  my (@sgenes, @didexons);
  @didexons=();
  
  @sgenes = map{ $genes->{$_} } grep{ $hasmajority{$_} } keys %$gbin;
  @sgenes = sort _sort_refscoreloc @sgenes; # need best 1st, for overgene test again
  foreach my $gref (@sgenes) {
    my $gid= $gref->[jGID];
    # next if($didkeep{$gid}); #? or $keepgene{$gid}
    my $flags= $gbin->{$gid};  
    my $rexons= $exons->{$gid};

    ## no cant assume this are all best, majority scored even for lesser genes
    my($nover)= $flags =~ m/nover=(\d+)/;
    my $skipit=0;
    if($nover>0) {
      my ($over)= overanygene( $gid, \@didexons);
      $skipit=1 if($over);
    }
    $flags .= "ismajor=1;";
    
    if($skipit) { 
      $nskip++; $skipgene{$gid}= $flags; 
      # delete $keepgene{$gid};  #? not there?
    } else {
      $didkeep{$gid}=1; $keepgene{$gid}++; $nkeep++; 
      putgene( $gref, $rexons, $otherft->{$gid}, $flags);  
      push(@didexons, @$rexons);
    }
  }
  
  @sgenes  = map{ $genes->{$_} } grep{ !$hasmajority{$_} } keys %$gbin;
  @sgenes = sort _sort_refscoreloc @sgenes; # need best 1st, for overgene test again
  foreach my $gref (@sgenes) {
    my $gid= $gref->[jGID];
    next if($didkeep{$gid}); #? or $keepgene{$gid}
    my $flags= $gbin->{$gid};  
    my $rexons= $exons->{$gid};
    my($nover)= $flags =~ m/nover=(\d+)/;
    
    my $skipit=0;
    if($flags =~ /majorid=([^;\s]+)/) { # not enough, also test overany
      my $ovid=$1; 
      $skipit=1 if($didkeep{$ovid});
    } ##else 
    unless($skipit) {
      my ($over)= overanygene( $gid, \@didexons);
      $skipit=1 if($over);
    }
    
    if($skipit) {
      $nskip++; $skipgene{$gid}= $flags; 
      #? delete $keepgene{$gid};
    } else {
      $didkeep{$gid}=1; $keepgene{$gid}++; $nkeep++; 
      putgene( $gref, $exons->{$gid}, $otherft->{$gid}, $flags);   
      push(@didexons, @$rexons);
    }
  }
  
}



sub overanygene {
  my($gid, $rdidexons)= @_;
  
  my $gref= $genes->{$gid} or return -1;
  my($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid1)= @$gref;
  my $rexons= $exons->{$gid} or return -1;
  ## my @didexons= @$rdidexons;

  my($nover, $xbover, $xbtotal, $cdsover, $cdsbover, $isum, $xbad, $ix, $lxe, $lovin) = (0) x 10;
  my %overid=();
  
  my @exons= grep { $_->[2] ne "CDS" } @$rexons; # CDS > exon
  @exons= @$rexons unless(@exons);
  @exons = sort _sort_exons @exons;
  my $nexons= @exons;  
  
  $ix=0; foreach my $ex (@exons) {  $ix++;
    my($xt,$xb,$xe)= ( $ex->[2], $ex->[3], $ex->[4]);
    my($ov,$ovin,$inerr,$overid)= 
      overlaps($ix==1, $ix==$nexons, $ref, $xt, $xb,$xe, $to, $rdidexons);
    $nover++ if($ov);
    $overid{$overid}++ if($ov); # for major
    $xbover += $ov;  $xbtotal+= 1 + $xe - $xb;
    $xbad++  if($badattr && $ex->[8] =~ m/$badattr/);
    
    ## always add ovin score as annot to exon/mrna
    if( $hasintrons ) {
      my $evin= ($ovin & 1) + ($lovin & 2); # should be == 3 for both
      $isum += $evin;
      if($LONG_INTRON>0 and $lxe>0 and ($xb - $lxe > $LONG_INTRON) and ($evin < 3) ) {
        $xbad++; 
      }
    }

    $lxe=$xe; $lovin= $ovin;
  }
  
  @exons =  sort _sort_exons grep { $_->[2] eq "CDS" } @$rexons; # CDS > exon
  my $ncds= @exons;  
  $ix=0; foreach my $ex (@exons) {  $ix++;
    my($ov,$ovin,$inerr,$overid)= 
      overlaps($ix==1, $ix==$ncds, $ref, $ex->[2], $ex->[3], $ex->[4], $to, $rdidexons);
    $cdsover++ if($ov);
    $cdsbover += $ov;  # $cdstotal += 1+$xe-$xb;
  }
  unless($nexons > 0) { $nover= $cdsover; $nexons= $ncds; $xbover=$cdsbover; }

  # if($nover > 0 and $typeover =~ /inside/) 
  if($nover > 0 and $typeover =~ /all|most/) { ##  and $typeover =~ /inside/ # allow also non-inside here          
    my $clen= $gref->[jCDSLEN];
    my $xlen= $gref->[jTRLEN];
    my $pcode= ($clen>0 and $xlen>0) ? int(100*$clen/$xlen) : 100;   
    
    if( $pcode < $PCDS_OK and $cdsover == $ncds) { } # skip; add pcode flag?
    elsif($typeover =~ /all/) { $nover=0 if($nover < $nexons); }
    elsif($typeover =~ /most/) { $nover=0 if($nover < 0.75 * $nexons); } # >= 3/4 
  }
 
  return($nover, $nexons, $xbover, $xbtotal, $ncds, $cdsover, $xbad, \%overid); # 
}


sub _sort_refscoreloc 
{
# ($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid)
# fixed: ?? with ref 1st this isn't working; big score not at front
  return ($a->[0] cmp $b->[0]) # ref
      || ($b->[5] <=> $a->[5]) # score
      || ($a->[3] <=> $b->[3]) # begin
      || ($b->[4] <=> $a->[4]); # end large>small
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


  # if($flags) { $gref->[8] = addreplattribs($gref->[8],$flag); }
sub dropattribs
{
  my($attr,$drop)= @_;
  return $attr unless($drop);
  my @key= grep /\w/, map{ my($k,$v)= split"=",$_,2; $k; } split ";", $drop;
  if(@key) { my $dkey= join('|',@key); $attr =~ s/;($dkey)=[^;\n]+//g; }
  return $attr;
}
  
sub addreplattribs
{
  my($attr,$add)= @_;
  return $attr unless($add);
  $attr= dropattribs($attr,$add);
  $attr =~ s/$/;$add/; # not .= may end in \n ??
  return $attr;
}


sub putgene {
  my($gref,$exons,$otherft,$flag)= @_;
  my $cskip= ($flag && $flag =~ /skip=/) ? "#x.": "";
  $otherft ||= [];

  # FIXME: strip old flag same-scores from mrna
  
  my $svec= $gref->[jSCOREVEC];
  if($svec and $flag !~ m/svec=/) {
    my($tscore,$scorefix) = manyscore_sum($svec); # update, add doesnt sum now
    $flag .= "svec=$scorefix"; 
    }
    
  foreach my $ex ($gref, @$exons, @$otherft) {
    next unless($ex and ref $ex);
    my($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid)= @$ex;
    $tattr =~ s/;$//; 
    if($flag) { 
      $tattr= addreplattribs($tattr,$flag); $flag="";  # only for gref
      }
    print join("\t", $cskip.$ref, $src, $typ, $tb, $te, $tscore, $to, $tph, $tattr),"\n";

    # FIXME Split= gene, has 2+ mRNA recs in @$ex; are split exons loc-sorted?
    # 2x13=26 array of  $rloc= [$ref,$src,$typ,$tb,$te,$score,$to,$tph,$tattr,$gid, 0,0,0];  
    if($typ eq 'mRNA' and $tattr=~/;Split=/ and scalar(@$ex) > 20 and $ex->[15] eq 'mRNA') {
      my @mrna2= @{$ex}[13..21];  
      print join("\t",@mrna2),"\n" if($mrna2[9] =~ /ID=/); # SHOULD place before same ref exon
      }
  }
}


sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }

sub overlaps {
  my($xfirst, $xlast, $ref, $typ, $tb, $te, $to, $locs, $oflags) = @_;
  $oflags||="";
  # change to return($overbases,...) == maxo ?
  
  my ($overin,$inerr)=(0,0);
  my $dointrons= ($hasintrons and not($xfirst and $xlast)) ? 1 : 0;
  # always check tb,te for intron supt? dont do for each overlap test.
  if($dointrons) {
    my($ov1, $ov2, $do1, $do2, $tb1,$te1, $tb2, $te2)= 
      ($to eq "-") ? # dont need?
        ( 2, 1, !$xlast, !$xfirst, $te - SPLICEX,$te, $tb,$tb + SPLICEX) 
      : ( 1, 2, !$xfirst, !$xlast, $tb,$tb + SPLICEX,$te - SPLICEX,$te);  
    
    if($do1) { 
      my $ovi=overintron($ref, $tb1, $te1, $to); 
      $overin += $ov1 if($ovi>0); $inerr += $ovi if($ovi<0); 
    }
    if($do2) { 
      my $ovi=overintron($ref, $tb2, $te2, $to); 
      $overin += $ov2 if($ovi>0); $inerr += $ovi if($ovi<0); 
    }
  }
  
  my $alloverid=($oflags =~ /allover/)?1:0;
  my($lovbases,$loverin,$linerr,@lid)= (0,0,0);
  foreach my $rloc (@$locs) {
    ## $rloc= [$ref,$src,$typ,$tb,$te,$score,$to,$tph,$tattr,$gid]; 
    my ($lref,$lty,$lb,$le,$lo,$lid)= @{$rloc}[0,2,3,4,6,9];    
    next if($lref ne $ref or ($typ and $lty ne $typ));
    
    my $ovbases= 0;
    my $over= ($tb <= $le && $te >= $lb) ? 1 : 0;
    $over=0 if($stranded and ($to =~ /[+-]/ and $lo =~ /[+-]/ and $to ne $lo)); # NEW 2010.10; fixme for to or lo eq "."
    
    if($over and $typeover =~ /inside/) { 
      # kINSIDE, add SAMEBASE ? allow few bases outside
      
      $over= ($tb >= $lb - SAMEBASE && $te <= $le + SAMEBASE) ? 1 : 0;
      ##$over= ($tb <= $lb && $te >= $le) ? 1 : 0;
      # ^^ reverse this caller is inside, locs are outside
      
      # add intron test here: if inside, not inside if tb, te match introns
      $over= 0 if($over and $overin>0 and $inerr==0 
            and ( (($overin & 1) && $tb > $lb + SAMEBASE ) 
              or  (($overin & 2) && $te < $le - SAMEBASE)) ); # match new intron
      $ovbases= 1+ $te - $tb if ($over); # == this.width

    } elsif($over and $pctover) {
      my ($bb,$be)= ( _max($tb,$lb), _min($te,$le) );
      my $maxo= 1+abs($be - $bb);
      my $leno= _max( 1, _min( abs($le - $lb), abs($te - $tb)) );
      $over = 0 if($maxo/$leno < $pctover);
      
      ## add intron test? * Careful, dont clear if exon is outside of range
      $over= 0 if($over and $overin>0 and $inerr==0 
            and ( (($overin & 1) && $tb > $lb + SAMEBASE) 
              or  (($overin & 2) && $te < $le - SAMEBASE)) ); # match new intron

      $ovbases= $maxo if($over);
    }
    
    #* change here to collect/return all @lid this is over ?
    if($over) {
      if($alloverid) { push @lid, $lid; 
        ($lovbases,$loverin,$linerr)= ($ovbases,$overin,$inerr) if($ovbases>$lovbases); 
      } else { return ($ovbases,$overin,$inerr,$lid); }
    }
  }
  
  return($lovbases,$loverin,$linerr,@lid) if(@lid);
  return (0,$overin,$inerr);
}


sub overintron {
  my($ref, $tb, $te, $to) = @_;  # tb,te == exon splice point here
  return 0 unless($introns->{$ref});
  my ($ib1, $ib2)= (int($tb/$BINSIZE) , int($te/$BINSIZE));
  for (my $ib = $ib1; $ib <= $ib2; $ib++) {
    $introns->{$ref}{$ib} or next;
    my @locs= @{$introns->{$ref}{$ib}};
    foreach my $rloc (@locs) {
      my ($lb,$le,$lo,$lid)= @{$rloc}[3,4,6,9];
      my $samestrand= ($to =~ /[+-]/ and $lo =~ /[+-]/ and $to ne $lo) ? 0 : 1;  
      #?? add -score for wrongstrand ?
      ## return 1 if($tb <= $le && $te >= $lb && $samestrand);
      if($tb <= $le && $te >= $lb) {
        return ($samestrand) ? 1 : -1;
      }
    }
  }
  return 0;
}


sub collect_gff
{
  my($gff)= @_;
  my (%genebins, %genes, %exons, %introns, %others);
  %genebins= %genes= %exons= %introns= %others= ();
  
  my ($ng,$nr)= (0,0);
  while(<$gff>){
    next unless(/^\w/); chomp;
    my $line= $_;
    my($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr)= split"\t";
    $nr++;
    $tscore=0 if($tscore eq ".");

    my $score= 0; # want this to select; w/ field choice?
    # should add scoreweights also: -score 'intron:3,coding:2,CDS:1'
    # my @scorename= split ",",$scoretype;    #multiscore .. see obestgene2
    my @scorevec= (0) x scalar(@scorefield);
    foreach my $sname (@scorefield) {
      my $sc=0;
      if($sname =~ /^score/i) { $sc=$tscore; }
      elsif($tattr =~ m/\b$sname=([^;\s]+)/) { $sc=$1; }
      ## if($sc =~ m/([\d\.\+\-]+)/) { $sc=$1; }  # handle this in addscore()
        # FIXME for some fields w/ multiple scores: utrx=0,4 == 0left, 4right bad utrs; inqual=100,9/... > 900 score
      # $score += $sc;
      addscore( $sc, $sname, \@scorevec);  
    }
    ($score)= manyscore_sum(\@scorevec); # update, add doesnt sum now
    $score= $tscore if($score == 0);
    
    my($gid); 
    if($tattr =~ m/\bID=([^;\s]+)/) {  $gid=$1; }
    elsif($tattr =~ m/\bParent=([^;\s]+)/) {  $gid=$1; }
    unless(defined $gid) { $gid = "N".$ng; }
    
    ## FIXME: 201405; add Split= gene handling.  ID= is same for splits, followed by ';Split=[12..];'
    ## .. cant replace genes{gid} w/ 2nd split; either make new ID.#split (for all recs) 
    ## .. or append exons to one gid record, but then mRNA rloc needs to be array in genes{gid}= rloc1,rloc2 ..
    my $issplit= ($tattr =~ m/;Split=(\d+)\b/)?$1:0;
    
    my $rloc= [$ref,$src,$typ,$tb,$te,$score,$to,$tph,$tattr,$gid]; # change to string; save mem
    
    if($typ =~ /^($mrnatypes)$/) { 
      $tattr= dropattribs($tattr,"isalt=1"); # others?
      $rloc= [$ref,$src,$typ,$tb,$te,$score,$to,$tph,$tattr,$gid, 0,0,0];  
      $rloc->[jSCOREVEC]= join(";",@scorevec);
      
      #ID=xxx;Split=[123..];
      if($issplit) { 
        if(my $rloc1=$genes{$gid}) { 
          my $rloc2= [ @$rloc1, @$rloc ]; # can we use append?
          $genes{$gid}= $rloc2;  
        } else {
         $genes{$gid}= $rloc; $ng++; # should be uniq/id ?
        }
      } else { 
        $genes{$gid}= $rloc; $ng++; # should be uniq/id ?
      }
      
      for( my $ib= int($tb/$BINSIZE); $ib <= int($te/$BINSIZE); $ib++) { 
        push( @{$genebins{$ref}{$ib}}, $rloc);  
        }
    } elsif($typ =~ /^($exontypes)$/) {
      push( @{$exons{$gid}}, $rloc);
      
    } elsif($typ =~ /^($introntypes)$/) {
      intronadd(\%introns,$rloc);
      
    } else {
      push( @{$others{$gid}}, $rloc);
    } 
  }

  warn"#collect_gff=$nr, ngene=$ng\n" if $debug;
  return (\%genebins, \%genes, \%exons, \%introns, \%others);
}


sub manyscorefields {
  my($keepwt)=@_;
  ## return() unless($manyscore);
  (my $sf= $scoretype) =~ s/^many.//; 
  if($keepwt) { return (split",",$sf); }
  else { return map{ s/:.+$//; $_ } split",",$sf; } # was s/:[\-\d]+//
}

#* update: dont sumscore here; use manyscore_sum
sub addscore {
  my( $val, $sname, $svec)= @_;
  #?? 2012jul BUG here, got all neg tscore ... all effect of high cdsindel counts, neg wt ?
  ## is val neg? doesnt look like it: svec=0,159,100.0,0,0,83; score=inqual:10,cdsindel:-6,pid:2,path:-1,CDS:2,UTR:2
  $svec=[] unless(ref $svec);
  
  ## dang bad for cdsindel vs CDS here.. both match \b$sname
  return $svec unless($scoretype =~ m/\b$sname\b/i); # was \b$sname; allow sname == 'UTR|coding'
  
  for(my $i=0; $i<@scorefield; $i++) {
    if($scorefield[$i] =~ m/\b$sname\b/i) {
      if($sname =~ /^($SPECIALFIELD)$/) {  } 
      elsif($val =~ m/([\d\.\+\-]+)/) { $val=$1; }  # handle this in addscore()
      $val =~ s=,=/=g; # cant have ,commas, here
      $svec->[$i]= $val; #as vector, no weigth
      last;
      
#       my $wt=  $scoreweight[$i];
#       # hack fix for UTR|coding : weight peak at ~80% with lower score for 95+, much lower/neg for <50%
#       if($scorefield[$i] =~ m/^(UTR|coding)/i) {
#         if($val > 94) { $score += -$wt * ($val-94);  }   
#         elsif($val < 65) { $score += -$wt * (65-$val); } 
#         else { $score += $wt * $val; }
#       } else {
#         $score += $wt * $val; #== manyscore_sum
#       }
    }
  }
  return $svec;
}

sub manyscore_sum { 
#? fixme for output svec, rewrite scorevec from transforms below
# FIXME, need to preserve orig score or know scorefix done
  my($svec)= @_;
  my $score= 0;  
  
  unless(ref $svec) {
    if($svec and $svec =~ m/[,;]/) {  #? split /[,;]/, $svec
      my @scorevec= split /[,;]/,, $svec; # jSCOREVEC packed
      $svec= \@scorevec;
    } else {
      return ($score,"");
    } 
  }
    
  my @scorefix= (0) x scalar(@scorefield); 
  my %svals; for(my $i=0; $i<@scorefield; $i++) { $svals{$scorefield[$i]}= $svec->[$i]||0; }
  
  for(my $i=0; $i<@scorefield; $i++) {
    my $sname=$scorefield[$i];
    my $val= $svec->[$i] or next; 
    my $wt=  $scoreweight[$i];
    my $val2= $val;
    
    # hack fix for UTR|coding : weight peak at ~80% with lower score for 95+, much lower/neg for <50%
    if($sname =~ m/^UTR$|^coding/) { # NOT utrx
      ## fixme2: set max score at UTR/pCDS == 80? score tail off both dir from that? see overbestgene2
if(1) {      
      my $cdsbest=80; my $cdstoosmall=60; my $cdstoobig=90;
      if($val > $cdstoobig) { $val=100 if($val>100); $val2= 100 - $val; } # pos, 0..10 
      elsif($val < $cdstoosmall) { $val=0 if($val<0); $val2= $val - $cdsbest } # neg -20 .. -80
      else { $val2= 100 - abs($cdsbest - $val); } # 100..80, +max at cdsbest, tails off to 80
      $score += $wt * $val2; 
} else {      
      if($val > 94) { $score += -$wt * ($val-94);  }   
      elsif($val < 65) { $score += -$wt * (65-$val); } 
      else { $score += $wt * $val; }
}
      
    } elsif($sname =~ m/^inqual/) {  #FIXME: SPECIALFIELD
      ## ;utrx=0,19;ocds=eqn;inqual=41,22/-8/4; ## inbad = negval here
      ## FIXME: high inqual for bad utrx should not be rewarded.. subtract utrx from inqual?
      ## .. velb9ptr013k37Loc41t2 is bad join: utrx=7,5;inqual=100,24/0/0;aalen=656,30%,complete
      ## my $v2=int($val); # damn noise: Argument "274/86%/partial3" isn't numeric  
      ($val2)= $val =~ m/([\d\.-]+)/;
      
      my $utrx = $svals{'utrx'} || 0;
      if($utrx =~ m/^(\d+).(\d+)/) { $utrx= $1 + $2; $utrx-- if($utrx>0); } else { $utrx=0; }
      if($val =~ m/^([\d-]+).([\d-]+).([\d-]+)/) { 
        my($ip,$ig,$ib)=($1,$2,$3); 
        my $inv= $ig + $ib; $inv -= $utrx if($inv>=$utrx);
        $val2= abs($ip) * $inv; 
        }
      $score += $wt * $val2; 
      
    } elsif($sname =~ m/^utrx/) {  #FIXME: SPECIALFIELD
      ($val2)= $val =~ m/([\d\.-]+)/;  # maybe multiply val x 100? or just use high neg wt.. utrx:-100
      if($val =~ m/^(\d+).(\d+)/) { $val2= $1 + $2; } # sum not _max
      $score += $wt * $val2; 
      
    } elsif($sname =~ m/^aalen/) {  #FIXME: SPECIALFIELD
      ##aalen=122,20%,complete
      ($val2)= $val =~ m/([\d\.-]+)/;
      # if($val =~ m/^(\d+).(\d+)/) { $val2= _max($1,$2); }
      $score += $wt * $val2; 
      
    } else {
      $score += $wt * $val;  
    }
    $scorefix[$i]= $val2;
  }
  my $scorefix= join",", @scorefix;  ## , for output  ; for jSCOREVEC ??
  return ($score, $scorefix);
}

=item overgene2 UTR scoring

  my $UTR_TOO_BIG = 0.60; # NOW, portion of transc len.; WAS 2, proportional to expected size, e.g 1000 bp for 400 bp expect

  UTR_USE_NEW_WEIGHT
      $utr= ($exw - $cdsw); $utr= 0 if($utr < 0);
      # change back to ratio utr/exw ?? or use that to moderate score
      
      my $UTR_BEST_SIZE= 0.5 * $UTR_TOO_BIG; ## if($UTR_BEST_SIZE >= $UTR_TOO_BIG);      
      my $utr_toobig = $UTR_TOO_BIG * $exw; # change $UTR_TOO_BIG to this 0.50
      my $utr_best   = $UTR_BEST_SIZE * $exw;  #? UTR_BEST_SIZE option
      my $nutr= $nexon - $ncds;
      my $utr_toomany = ($nutr > 4) ? 1 : 0; # need left,right utr counts : look for utrx=5,2 (l,r) field

      # if nutr very big, drop this model always ?? or only if no alt models at locus
      # but some are valid ncrna

      my $utrscore= 0;
      if($utr <= $utr_best) 
      {
        $utrscore = $utr; ## == $utr_expected_size * $utr/$utr_expected_size; # positive
        $utrscore -= $utr * ($nutr - 3) if($utr_toomany);
      } 
      elsif( $utr >= $utr_toobig ) 
      { # too large, neg weight?
        $utrscore = ($utr_toobig - $utr) * 0.5; # * negative
        $utrscore -= $utr * ($nutr - 3) if($utr_toomany);
      } 
      else # > best < toobig
      {  
        $utrscore = $utr_toobig - $utr; # positive
        $utrscore -= $utr * ($nutr - 3) if($utr_toomany); # make negative?
      }
      $utr= int($utrscore);

=cut

sub manyscore_sum0 {
  my($score)= @_;
  my @as= split",",$score;
  my $n= scalar(@as);
  my($as1)= (0);
  for( my $i=0; $i<$n; $i++) {
    my($as2)=($as[$i]);  
    $as2 =~ s,/.*,,;   # $as2 =~ s,[^\d+-e\.].*,,;
    $as1 += $scoreweight[$i] * $as2;
  }
  return int($as1);
}


# sub getScoretypes( $scoretype, \@scoreweight, \@dropscore, \@droptype, \@genegroup );
# 10/09/18: can we add -scoreweight, for terepeats? yes
 
sub getScoretypes {
  my $i=0; 
  my @sf= manyscorefields(1);
  @scoreweight = (1) x scalar(@sf);
  foreach my $sf (@sf) {  
    if($sf =~ m/(\w+):([\d\.\-]+)/) { $scoreweight[$i]=$2; } #? allow fractions? 0.5  
    $i++;
    } 
  @scorefield= manyscorefields(0);

  @dropscore= (0) x scalar(@scorefield);
  @droptype = (0) x scalar(@scorefield);
  if($dropscore) {
    $dropscore =~ s/^many\.//; 
    foreach my $i (0..$#scorefield) { 
      my $sf= $scorefield[$i]; 
      my($dcode,$dval)=(undef,undef); 
      if($dropscore =~ m/$sf[:=]([+])(\d+)/) { ($dcode,$dval)=($1,$2); } # user-error fix...
      elsif($dropscore =~ m/(\W?)$sf[:=](\d+)/) { ($dcode,$dval)=($1,$2); }
      if(defined $dval) {
        $dcode= (not defined $dcode) ? kDROP_OR 
          : ($dcode eq "+") ? kDROP_AND : ($dcode eq "-") ? kDROP_NOT
          : ($dcode eq "*") ? kDROP_MUSTNOT : kDROP_OR;
        $dropscore[$i]= $dval; 
        $droptype[$i]= $dcode; 
        }
      }
    }

}


sub scoregene {
  my($gid,$genes,$exons,$introns,$otherft)= @_; # ,$flags
  # $flags="CDS,exon" unless($flags);
  
  my $gref= $genes->{$gid} or return;
    
  #?? also score introns, add to mrna score as per overbestgene2? see above
  ## iscore as weight: (valid splices - error splices) / total splices ??
  my $flags="";
  my( $incode, $insum, $ierrsum, $intotal)= (0) x 10;
  if($hasintrons) {
    my $iflag=""; 
    my($nover, $cdsover, $xbad, $ix, $lxe, $lovin) = (0) x 10;
    my $rexons= $exons->{$gid};
    my @exons= grep { $_->[2] ne "CDS" } @$rexons; # CDS > exon
    @exons= @$rexons unless(@exons);
    @exons = sort _sort_exons @exons;
    
    my $nexons= @exons;  
    $ix=0; foreach my $ex (@exons) {  $ix++;
      my($ref, $xt,$xb,$xe, $to)= ( $ex->[0], $ex->[2], $ex->[3], $ex->[4], $ex->[6]);
      my($ov,$ovin,$inerr)= overlaps($ix==1, $ix==$nexons, $ref, $xt, $xb,$xe, $to, []);
      if( 1 ) { # $hasintrons
        my $evin= ($ovin & 1) + ($lovin & 2); # should be == 3 for both
        
        ##$insum   += $evin;  # 1+2 splices w/ valid intron ends
        ##$intotal += (($ix==1) ? 0 : 1) + (($ix<2) ? 0 : 2); # match 1+2 of evin, or count 1 for each splice
        #..
        $insum++ if($ovin & 1); $insum++ if($ovin & 2); # == found splices
        $intotal += (($ix==1) ? 0 : 1) + (($ix==$nexons) ? 0 : 1); # == num splices
        
        $ierrsum += $inerr; # minus mistakes
        my $ints= ($inerr) ? "$inerr/$evin" : $evin;
        if($LONG_INTRON>0 and $lxe>0 and ($xb - $lxe > $LONG_INTRON) and ($evin < 3) ) {
          # $xbad++; 
          $ints .=",longerr:".($xb - $lxe);
        }
        unless($ints eq "0") {
        $iflag .= "i$ix:$ints,"; # $evin/$inerr,"; 
        $ex->[8] =~ s/$/;ints=$ints/;
        }
      }

      $lxe=$xe; $lovin= $ovin;
    }
    if($intotal>0) {
      $incode= int (100 * ($insum + $ierrsum) / $intotal); # can be -
      $flags .= "ints=$incode," . (($ierrsum) ? "$ierrsum/$insum/$intotal" : "$insum/$intotal") .",$iflag;";
    } else {
      $incode= 1;
    }
  }
  
  my($cw,$cdsg,$xw,$xong,$xn,$cn)= (0) x 9;
  if($gref->[8] =~ m/cxlen=/) {  # FIXME: now genefindcds lists cxlen= from transcript; want also GFF cds/exon cxlen score?
    ($cw,$xw)= $gref->[8] =~ m,cxlen=(\d+).(\d+),;  
  } 
  { #now always count cdsg xong
    # my($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid)= @$gref;
    my $rexons= $exons->{$gid};
    my $rmore = $otherft->{$gid};
    ##return unless(ref $rexons); # are we missing exon-gene id links?
    if(ref $rexons) {
    foreach my $ex (@$rexons, @$rmore) {
      next unless(ref $ex);
      my($t,$b,$e)= ($ex->[2], $ex->[3], $ex->[4]);
      my $w=1 + $e - $b;
      if($t eq "CDS") { $cdsg+= $w; $cn++; }
      elsif($t eq "exon") { $xong+= $w; $xn++; }
      }
    # $xw=$cw if($xw==0); # or not? for cds-only genes
    if($xw or $cw) { 
      # have annot
    } else {
      ($cw,$xw)= ($cdsg,$xong);
      $flags .= "cxlen=$cw/$xw;";
      }
    }
  }
  
  # FIXME: strip old same-scores 
  # $gref->[8] =~ s,$,;$flags, if($flags);
  if($flags) { $gref->[8] = addreplattribs($gref->[8],$flags); }

  $xw= _max($cw,$xw);
  my $pcode= ($cw>0 and $xw>0) ? int(100*$cw/$xw) : 1;  # use as sort score? should merge w/ tscore at collect_gff

    ##? do like best2, make vector of score?  20,10,0,9 ; so know each value?
  my $tscore= $gref->[5]; # check numeric?
  my @scorevec;
  @scorevec= (0) x scalar(@scorefield);
  @scorevec= split";", $gref->[jSCOREVEC] if($gref->[jSCOREVEC]); # FIXME, need to preserve orig score or know scorefix done
  # $scorevec[$k] = $tscore; #what field?
  addscore( $cw, 'CDS', \@scorevec); #  if($scoretype =~ /\bCDS/i)
  addscore( $xw, 'exon', \@scorevec);
  addscore( $cdsg, 'cdsg', \@scorevec); # special case for CDS from transcript not genome CDS ..
  addscore( $pcode, 'UTR|coding', \@scorevec);  # pcode reverse of prior UTR score 
  addscore( $incode, 'intron', \@scorevec);    
  ## my $scorefix=[];
  ($tscore) = manyscore_sum(\@scorevec); # update, add doesnt sum now
  #?? 2012jul BUG? here, got all neg tscore ... no hi count of neg wt field

  $gref->[5]= $tscore; 
  $gref->[jCDSLEN]= $cw;  
  $gref->[jTRLEN]= $xw;  
  $gref->[jSCOREVEC]= join(";",@scorevec);# FIXME, need to preserve orig score or know scorefix done

  scoredrops($gid, \@scorevec) if($dropscore);
}


sub scoredrops
{
  my($gid,$scorevec)= @_;
  
  my $keep=0;
# if(USE_DROPSCORE_ANDOR) 
  {        
  for(my $i=0; $i < @scorefield; $i++) { 
    my($ds,$sv)= ($dropscore[$i],$scorevec->[$i]);
    my $droptype= $droptype[$i] || 0; # reorder @sfields to put require/AND, reject/NOT at end
    map{ s=[,;/].*== }($ds,$sv); # "$pv/$px" # Argument "72/232" isn't numeric in numeric
    
    if(not defined $ds or $ds == 0) {
    
    } elsif($droptype == kDROP_AND) {
      if ($sv < $ds) { $keep= -1; } 
      elsif ($sv >= $ds and $keep >= 0 ) { $keep=1; } 
    
    } elsif($droptype == kDROP_NOT) { #? is this useful
      $keep= -1 if ($sv >= $ds);           
     
    } elsif($droptype == kDROP_MUSTNOT) { # add
      do{ $keep=1; last} if ($sv >= $ds);           
     
    } else { # kDROP_OR
      $keep= 1 if ($sv >= $ds and $keep >= 0);           
    }
  }
  $skipgene{$gid}=kSKIP_NOEVD unless($keep > 0); # 999 = flag dont unskip this one
    
} 
# else {  # original, OR only
#   for(my $i=0; $i < @scorefield; $i++) { 
#     my($ds,$sv)= ($dropscore[$i],$scorevec->[$i]);
#     map{ s=[,;/].*== }($ds,$sv); # "$pv/$px" # Argument "72/232" isn't numeric in numeric
#     # want -dropscore?  terepeat:-100 means drop if >=100 TE bases 
#     $keep=1 if ($ds > 0 and $sv >= $ds); 
#     }
#   $skipgene{$gid}=kSKIP_NOEVD if($keep==0); # 999 = flag dont unskip this one
# }

}


sub intronadd {
  my($introns, $rloc)= @_;
  my($ref,$src,$typ,$tb,$te,$score,$to,$tph,$tattr,$gid)= @$rloc;
  
  my($s1b,$s1e,$s2b,$s2e)= 
    ($to eq "-") ? ($te+1,$te + SPLICE,$tb - SPLICE,$tb-1) 
    : ($tb - SPLICE,$tb-1,$te+1,$te + SPLICE); # 3bp + 1 shift

  ($tb,$te)= ($s1b,$s1e);
  my $rloc1= [$ref,$src,$typ,$tb,$te,$score,$to,$tph,$tattr,$gid]; 
  for( my $ib= int($tb/$BINSIZE); $ib <= int($te/$BINSIZE); $ib++) { 
    push( @{$introns->{$ref}{$ib}}, $rloc1); }

  ($tb,$te)= ($s2b,$s2e);   
  my $rloc2= [$ref,$src,$typ,$tb,$te,$score,$to,$tph,$tattr,$gid]; 
  for( my $ib= int($tb/$BINSIZE); $ib <= int($te/$BINSIZE); $ib++) { 
    push( @{$introns->{$ref}{$ib}}, $rloc2); }
}



