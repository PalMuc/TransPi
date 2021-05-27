#!/usr/bin/env perl
# gff2hints.pl
# convert a GFF evidence file 
# to augustus hints file 
# adapted from blat2hints.pl by d.g. gilbert, 2007.12

## add conversion of other gene predictions to hints: esp mRNA span > start,stop hints
## tho these are CDS start/stop more often than mRNA start,stop
## ?? add intergene 'irpart' hints from gene preds? augustus may be joining genes too much

use strict;
use Getopt::Long;

my $usage = <<"EOU";
$0 -- convert GFF v3 file with EST or Protein alignments, or Gene predictions to hints file for AUGUSTUS

Usage: $0 --in=gfffile|stdin --out=hintsfile
  PREREQUISITE: input GFF file must be sorted by target (=genome) sequence names
  and should be sorted within the sequences by begin coordinates for efficiency
  options:
  --type=s           GFF type(s) to use (default all; e.g. type=EST_match,match_part  protein_match,HSP ...)
  --class=s          EST|protein|prediction (hint-type: exon or CDS, use with type and source)
    --sortmodels       GFF must be sorted by location, gene model (prediction class mixed sources?)
  --priority=n       priority of hint group (default 4)
  --minintronlen=n   alignments with gaps shorter than this and longer than maxgaplen are discarded (default 41)
  --maxintronlen=n   alignments with longer gaps are discarded (default 350000
  --ep_cutoff=n      this many bp are cut off of each exonpart hint at end of alignment (default 10)
  --source=s         source identifier (default 'E', 'P' for protein)
  --remove_redundant only keep the strongest hint for a region (default false)
  --maxcoverage=n    maximal number of hints at a given position. A high value causes long running time of
                     AUGUSTUS in regions with thousands of cDNA alignments. (default 3000)
  --ssOn             include splice site (dss, ass) hints in output (default false)
EOU
#  --maxQgaplen=n     maximum length of gap in query (cDNA) sequence
#  --maxgaplen=n      gaps at most this length are simply closed (default 14)

my $gfffilename;
my $hintsfilename;
my $minintronlen = 25; #? 41;
my $maxintronlen = 350000; #? not used
my $start_stop_radius = 9; # was 15?

my $ip_cut = 4; # at least 1; inset for intron parts from exon end/start
my $ep_cutoff =  9;
my $MAX_IRPART = 600; # for gene models; skip big ir regions, want ability to predict inside these; note this does ir hint only at start of gene

#my $maxgaplen = 14;
#my $maxQgaplen = 5;

my $min_endblock_len = 6;
my $source="E";
my $priority = 4;
my $prgsrc_prefix = "g2h";
my $prgsrc= $prgsrc_prefix;
my $line=0;
my $coloffset=0;
my $remove_redundant=0;
my @coverage=();
my $maxcoverage = 3000;
my $ssOn=-1;
my $cdsOn=-1;
my $dosortmodels=0; # option

my $exonpartid = "ep"; # abbreviate to decrease file size
my $types="";
my $class="EST";

## use CDS, CDSpart for Protein input

if ($#ARGV < 1 ) {
    print "$usage";
    exit;
}


my $ok= GetOptions(
           'in=s'=>\$gfffilename,
           'out=s'=>\$hintsfilename,
           'types=s'=>\$types,
           'minintronlen:i'=>\$minintronlen, # not used; merge w/ maxQgaplen
           'maxintronlen:i'=>\$maxintronlen, # not used; change
           ##'maxgaplen:i'=>\$maxgaplen, # not used
           ##'maxQgaplen:i'=>\$maxQgaplen, # should be same as minintronlen
           'ep_cutoff:i'=>\$ep_cutoff,
           'ip_cutoff:i'=>\$ip_cut,
           'source=s'=>\$source, # 'P' means exon > CDS,CDSpart ? ..
           'class=s'=>\$class, # protein, EST, usetype
           'priority:i'=>\$priority,
           'remove_redundant!'=>\$remove_redundant,
           'maxcoverage:i'=>\$maxcoverage,
           'ssOn!'=>\$ssOn,
           'cdsOn!'=>\$cdsOn,
           'sortmodels!'=>\$dosortmodels,
           );
$ok or die $usage;

if($gfffilename =~ /.gz$/){ $ok= open(BLAT,"gunzip -c $gfffilename|"); }
elsif($gfffilename =~ /^(stdin|-)$/) { $ok=open(BLAT,"<&STDIN"); } #  reopen perl ..
else { $ok= open(BLAT, "<$gfffilename") ; }
$ok or die "Couldn't open $gfffilename\n";

if($ip_cut < 1) { warn "ip_cut must be >=1; reset to 1\n"; $ip_cut=1; }

# stdout ok here
if(!$hintsfilename or $hintsfilename =~ /^(stdout|-)$/) { $ok=open(HINTS,">&STDOUT"); } #  reopen perl ..
else { $ok= open(HINTS, ">$hintsfilename"); }
$ok or die "Could not open $hintsfilename";

## most of these shouldnt be globals
my (@dsshints, @asshints, @exonhints, @exonparthints, @intronhints);
# my ($i, $j, $mstart, $mend, $badalignment, $gaplen);
# my ($match,$TgapCount,$strand,$qname,$qsize,$blockSizes,$tStarts, $qStarts, $tstart, $tend);
# my (@f,@b,@t,@q);
# my (@blockbegins, @blockends);
# my $numBlocks;

# hint lists are sorted by by increasing begin position
my @hint; # (begin, end, strand, tname, qname)
my $hintref;
my ($targetname, $oldtargetname);
$oldtargetname = "no name yet";
my $skiplines=0;

# for gff
my $gvers=0;
# my( $tsource, $ttype, $score, $phase, $attr);
my ($lastqname);
my (@starthints, @stophints, @tsshints, @ttshints,  @intronparthints, 
    @CDShints, @CDSparthints, @UTRparthints, @irparthints);
my @GroupKeys = qw(Parent Target ID); #? need choice?
my @parts=();
my @lastparts=();

my $isgenemodel= ($class =~ /pred/i) ? 1 : 0; # global?
my $isprotalign= ($class =~ /prot/i) ? 1 : 0; # global?
  my $isprotfull= ($isprotalign && $class =~ /full/i) ? 1 : 0;
  #^^ need -sort for this
my $isestalign= ($class =~ /EST/i) ? 1 : 0; # global?

if($isgenemodel) {
  # if($ssOn<0) { $ssOn=1; } # not useful
  # if($cdsOn<0) { $cdsOn=1; } # want CDSparts/cp only !
  unless($types){ $types="exon,CDS"; } # no others?
}
if($isestalign) {
  # if($ssOn<0) { $ssOn=1; } # not useful?
}
$ssOn=0 if($ssOn<0);
$cdsOn=0 if($cdsOn<0);

my %types= map{ $_,1 } split/[,;\| ]/, $types;
my @allparts=();

 
# worry about sorting of gff: for gene models, need each gm grouped together,
# but also gm sorted by location of full models. unix sort wont do this for gff w/ mixed predictors
sub _sort_genemodels {
  # in: a,b short gff: [$ref, $tsource, $ttype, $tstart, $tend, $score, $strand, $qname]
  my $ana= $a->[7]; my $bna= $b->[7];
  if($ana eq $bna) { # same model
    my($at)= grep(/$a->[2]/, ("4gene","3mRNA","2exon","1CDS"));
    my($bt)= grep(/$b->[2]/, ("4gene","3mRNA","2exon","1CDS"));
    return $a->[0] cmp $b->[0] #ref
      or $bt cmp $at #type : gene > mRNA > exon > CDS > other
      or $a->[3] <=> $b->[3] #start
      or $b->[4] <=> $a->[4] #end
      ;
  } else {
    return $a->[0] cmp $b->[0] #ref
      or $a->[3] <=> $b->[3] #start
      or $b->[4] <=> $a->[4] #end
      or $ana cmp $bna #name
      or $a->[1] cmp $b->[1] #source
      ;
  }
}


while (<BLAT>) {

  unless(/^\w/) {
    if (/^##gff-version\s+(\d+)/){
      $gvers=$1;
      die "Cannot yet handle gff-version $gvers\n" if($gvers<3);
      }
    next;
    }

    $line++;
    chomp;
    my @f = split /\t/, $_; 
    unless (@f == 9) { warn "Not GFF format"; next }  
    
    my($targetname, $tsource, $ttype, $tstart, $tend, $score, $strand, $phase, $attr)= @f;
#     $match       = $f[0];
#     $TgapCount   = $f[6];
#     $strand      = $f[8];
#     $qname       = $f[9];
#     $qsize       = $f[10];
#     $targetname  = $f[13];
#     $tstart      = $f[15];
#     $tend        = $f[16];
#     $blockSizes  = $f[18];
#     $qStarts     = $f[19];
#     $tStarts     = $f[20];

    next if($types and not $types{$ttype});

    my %attr= map{ my($k,$v)=split"=",$_,2; $k=>$v; } split ";", $attr;
    my($qname)= grep /\S/, @attr{@GroupKeys}; 
    $qname =~ s/\s.*$//; # for Target=name start end jaz
    $prgsrc= "$prgsrc_prefix.$tsource";
    my $thispart= [$targetname, $tsource, $ttype, $tstart, $tend, $score, $strand, $qname];
    
    # for gene model predictions, preserve last, next genespan to get intergene (irpart) hints
    if($qname ne $lastqname) { 
        if($dosortmodels) { push(@allparts, @parts); @parts=(); }
        if(@parts>0) { processparts($lastqname, \@parts, \@lastparts, [$thispart] ); }
        @lastparts= @parts; @parts=();  $lastqname= $qname;
    }
    
    if ($targetname ne $oldtargetname) {
      if($dosortmodels) { 
        my @sparts= sort _sort_genemodels @allparts;
        $lastqname=""; @parts=(); @lastparts= @parts; 
        foreach my $i (0..$#sparts) {
          my $newname= $sparts[$i]->[7]; 
          if($i>0 and $newname ne $lastqname) {
            processparts($lastqname, \@parts, \@lastparts, []);
            @lastparts= @parts; @parts=();
            }
          push(@parts, $sparts[$i]); $lastqname= $newname; 
          }
        processparts($lastqname, \@parts, \@lastparts, []); 
        @lastparts= @parts=(); 
        }
        printHints();
        $#coverage = -1;
        @lastparts=(); @allparts=(); 
    }

    push @parts, $thispart;
    $oldtargetname = $targetname;
    $lastqname= $qname;
} # while input

# print "\n";

if($dosortmodels) { push(@allparts, @parts); @parts=(); }
if(@parts>0) { processparts($lastqname, \@parts, \@lastparts, []); }
if($dosortmodels) { 
  my @sparts= sort _sort_genemodels @allparts;
  $lastqname=""; @parts=(); @lastparts= @parts; 
  foreach my $i (0..$#sparts) {
    my $newname= $sparts[$i]->[7];
    if($i>0 and $newname ne $lastqname) {
      processparts($lastqname, \@parts, \@lastparts, []);
      @lastparts= @parts; @parts=();
      }
    push(@parts, $sparts[$i]); $lastqname= $newname; 
    }
  processparts($lastqname, \@parts, \@lastparts, []);
  @lastparts= @parts=(); 
  }
printHints();


###########################################################################################
#
# subroutines
#
###########################################################################################



sub addSignalHint; # = \&addSSHint;

sub processparts 
{
  my($callname, $parts, $lastparts, $nextparts)= @_;
  # we are mixing up parts names from callname ??
  
  # find span of parts ...
  my ($tname, $tstart, $tend, $tstrand);
  my $numBlocks= scalar(@$parts);
  
  foreach my $pt (@$parts) {
    my($tref, $tsource, $ttype, $pstart, $pend, $score, $strand, $pqname)= @$pt;
        #     my $thispart= [$targetname, $tsource, $ttype, $tstart, $tend, $score, $strand, $name];
    $tstart= $pstart if(!$tstart || $pstart<$tstart);
    $tend= $pend if($pend>$tend);
    $tstrand= $strand if($strand && $strand ne '.');
    $tname= $pqname unless($tname); # should be same for all
    
    my $filterout=0;
    for (my $i = int($pstart/10); $i <= int($pend/10) && !$filterout;$i++) {
        if ($coverage[$i] >= $maxcoverage) {
            $filterout=1;
        }
    }
    if ($filterout) {
        return;
    }
    for (my $i=int($pstart/10); $i <= int($pend/10); $i++) {
        if (defined $coverage[$i]) {
            $coverage[$i]++;
        } else {
            $coverage[$i]=1;
        }
    }
    }


    # now add the hints
    # $numBlocks = scalar @blockbegins;
    my($cdsStart, $cdsEnd)=(0,0);
    my($geneStart, $geneEnd)=(0,0);
    my($mref, $msrc, $mtype, $mstart, $mend, $score, $strand, $pname);

## FIXME: should always shrink size of _parts hints due to way augustus measures these:
## hint is valid if aug.predict *contains* part, invalid if part extends over either end
## versus for CDS,exon,intron valid must match exact start,stop aug.prediction

## rewrite below to handle CDS and exon in same loop, but distinguish @hints
    my $exoniscds=0;
    if($isgenemodel) {
      my @cds=();
      my @exon=();
      foreach my $pt (@$parts) {
        my $mtype= $pt->[2];
        if($mtype =~ /CDS/) { 
          push(@cds,$pt); 
          ($mref, $msrc, $mtype, $mstart, $mend, $score, $strand, $pname)= @$pt; 
          $cdsStart= $mstart if(!$cdsStart or $cdsStart>$mstart);
          $cdsEnd= $mend if($cdsEnd<$mend);
          #? add +- $ep_cutoff  only at ends?
          @hint = ($mstart+$ep_cutoff, $mend-$ep_cutoff , $strand, $pname, $score, $msrc);
          addIntervalHint(\@CDSparthints, [@hint]);
          # ^ CDSpart/cp not CDS
          } 
        else { push(@exon,$pt); }
        
      }
      if(@exon==0) { @exon= @cds; $exoniscds=1; } # ? change type to exon?
      $parts= \@exon;
      $numBlocks= scalar(@exon);
    }
    
    # ensure parts are loc sorted; where is this problem?
    my @spart= sort{$a->[3] <=> $b->[3]} @$parts;
    $parts= \@spart;
    
    for (my $i=0; $i<$numBlocks; $i++) {
    
        ($mref, $msrc, $mtype, $mstart, $mend, $score, $strand, $pname)= @{$parts->[$i]};  
        #     my $thispart= [$targetname, $tsource, $ttype, $tstart, $tend, $score, $strand, $name];
        ## dang; now have sorted by mtype before mstart ! FIX
        ## need to respect mtype inside genemodel block
        # ?? split out exon and CDS sub parts ?
        
        #my $bstart= $mstart; # $blockbegins[$i];
        #my $bend  = $mend;   # $blockends[$i];
        my $bnext = ($i+1<$numBlocks) ? $parts->[$i+1]->[3] : $mend; # $blockbegins[$i+1] : 0;
        my $bnend = ($i+1<$numBlocks) ? $parts->[$i+1]->[4] : $mend;  # $blockends[$i+1] : 0;
        
        $geneStart= $mstart if(!$geneStart or $geneStart>$mstart);
        $geneEnd= $mend if($geneEnd<$mend);

        if ($i==0 && $i==$numBlocks-1 && !$exoniscds) {
            # just one exonpart, should not happen when spliced EST
            if ($mstart + 2*$ep_cutoff <= $mend ){
                @hint = ($mstart+$ep_cutoff, $mend -$ep_cutoff, $strand, $pname, $score, $msrc);
                addExonpartHint([@hint]);
            }
            
        } elsif ($i==0) {
            # first block
            if ($mstart + $min_endblock_len-1 <= $mend ){
                if ($mstart + $ep_cutoff <= $mend && !$exoniscds){
                    @hint = ($mstart+$ep_cutoff, $mend , $strand, $pname, $score, $msrc);
                    addExonpartHint([@hint]);
                }
                
                if ($ssOn) {
                  if($strand ne '-') {
                    @hint = ($mend +1, $mend +1, $strand, $pname, '.', $msrc);
                    addSignalHint(\@dsshints, [@hint]);
                    @hint = ($bnext-1, $bnext-1, $strand, $pname, '.', $msrc);
                    addSignalHint(\@asshints, [@hint]);
                  } else {
                    @hint = ($mend +1, $mend +1, $strand, $pname, '.', $msrc);
                    addSignalHint(\@asshints, [@hint]);
                    @hint = ($bnext-1, $bnext-1, $strand, $pname, '.', $msrc);
                    addSignalHint(\@dsshints, [@hint]);
                    }
                }

                my $isintron= $mend + $minintronlen < $bnext;
                if ($isintron && ($i<$numBlocks-2 || $bnend-$bnext+1 > $min_endblock_len)) {
                    #? add +- $ep_cutoff  only at ends?
                    @hint = ($mend +$ip_cut, $bnext-$ip_cut, $strand, $pname, $score, $msrc);
                    addIntervalHint(\@intronhints, [@hint]);
                    # ^ turn into ip
                }
            }
            
        } elsif ($i==$numBlocks-1) {
            # last block
            if ($mend  - $min_endblock_len + 1 >= $mstart){
                if ($mstart <= $mend -$ep_cutoff && !$exoniscds){
                    @hint = ($mstart, $mend -$ep_cutoff, $strand, $pname, $score, $msrc);
                    addExonpartHint([@hint]);
                }
            }
            
        } else { 
            # internal block, add following intron hint
#             @hint = ($mstart, $mend , $strand, $pname, $score, $msrc);
#             addIntervalHint(\@exonhints, [@hint]);
#              # ^ turn into addExonpartHint
          #? add +- $ep_cutoff  only at ends?
            @hint = ($mstart, $mend , $strand, $pname, $score, $msrc);
            addExonpartHint([@hint]) unless($exoniscds);

            my $isintron= $mend + $minintronlen < $bnext;
            if ($isintron && ($i<$numBlocks-2 || $bnend-$bnext+1 > $min_endblock_len)) {
                #? add +- $ep_cutoff  only at ends?
                @hint = ($mend +$ip_cut, $bnext-$ip_cut, $strand, $pname, $score, $msrc);
                addIntervalHint(\@intronhints, [@hint]);
                # ^^ turn into intronpart/ip
                
                if ($ssOn){
                  if($strand ne '-') {
                    @hint = ($mend +1, $mend +1, $strand, $pname, '.', $msrc);
                    addSignalHint(\@dsshints, [@hint]);
                    @hint = ($bnext-1, $bnext-1, $strand, $pname, '.', $msrc);
                    addSignalHint(\@asshints, [@hint]);
                   } else {
                    @hint = ($bnext-1, $bnext-1, $strand, $pname, '.', $msrc);
                    addSignalHint(\@dsshints, [@hint]);
                    @hint = ($mend +1, $mend +1, $strand, $pname, '.', $msrc);
                    addSignalHint(\@asshints, [@hint]);
                   }
                }
            }
        }
    }
    
        # $q[$i] + $b[$i] >= $q[$i+1] - maxgap   >> q = query/Target starts; b = block/hsp sizes
        # gff: part[i].end >= part[i+1].start ??
        ## oops ; i not == parts now; due to above gap checks; drop that check?
        
  ## for gene preds, add start/stop signals ; this is for CDS not transcript start/stop
  ## these don't look useful in tests, not default
  if($cdsOn) {  ## $isgenemodel ...
  if ($tstrand eq '+') {
    if ( $cdsStart > 0) {
        @hint = ($cdsStart-$start_stop_radius, $cdsStart+2+$start_stop_radius, $tstrand, $tname, '.', $msrc);
        addSignalHint(\@starthints, [@hint]);
      }
    if ($cdsEnd > 0) {
        @hint = ($cdsEnd-2-$start_stop_radius, $cdsEnd+$start_stop_radius, $tstrand, $tname, '.', $msrc);
        addSignalHint(\@stophints, [@hint]);
      }
    } elsif($tstrand eq '-'){
    if ( $cdsStart > 0) {
        @hint = ($cdsStart-$start_stop_radius, $cdsStart+2+$start_stop_radius, $tstrand, $tname, '.', $msrc);
        addSignalHint(\@stophints, [@hint]);
      }
    if ($cdsEnd > 0) {
        @hint = ($cdsEnd-2-$start_stop_radius, $cdsEnd+$start_stop_radius, $tstrand, $tname, '.', $msrc);
        addSignalHint(\@starthints, [@hint]);
      }
    }
  }

    ## this is no good; augustus invalidates all these due I guess to overlap cds/exon
    ## need intron/intra-match regions only for -irpart ??
#   if($isprotfull && $numBlocks>1 && $geneStart>0 && $geneEnd>$geneStart) {
#     my $irname= "nir.".$tname;
#     @hint = ($geneStart, $geneEnd , '.', $irname, -9, $msrc); # neg score to make it not-irpart
#     addIntervalHint(\@irparthints, [@hint]); 
#   }
  
  if($isgenemodel && $geneStart>0 && $geneEnd>0 && $lastparts) {
    my($lastGeneStart, $lastGeneEnd)=(0,0);
    my($mref, $msrc, $mtype, $mstart, $mend, $score, $strand, $pname);
    foreach my $lastpart (@$lastparts) {
        ($mref, $msrc, $mtype, $mstart, $mend, $score, $strand, $pname)= @{$lastpart};
        # $lastGeneStart= $mstart if(!$lastGeneStart or $lastGeneStart>$mstart);
        $lastGeneEnd= $mend if($lastGeneEnd<$mend);
        }

## note: irpart hint should chop off possible UTR if not in prediction ...
    if($lastGeneEnd>0 && $lastGeneEnd+4 < $geneStart && !$exoniscds) {
      my $irname= "ir.".$pname.".".$tname;
      my $irstart= $lastGeneEnd+2;
      if($lastGeneEnd < $geneStart - $MAX_IRPART) { $irstart= $geneStart - $MAX_IRPART; } # ~ 600 bp
      @hint = ($irstart, $geneStart-2 , '.', $irname, ".", $msrc);
      addIntervalHint(\@irparthints, [@hint]); 
      if($lastGeneEnd < $geneStart - $MAX_IRPART) {
        my $irend= $lastGeneEnd + $MAX_IRPART;
        @hint = ($lastGeneEnd+2, $irend, '.', $irname.".1", ".", $msrc);
        addIntervalHint(\@irparthints, [@hint]); 
        }
    }
  }
    
}

#
# printHints
# print and delete the hints
#
sub printHints {
    # finished computing the hints for the old sequence. output them.
    if($remove_redundant){
        # delete all exonpart hints that are contained in an exon hint
        my $startidx=0;
        my $curidx;
        foreach my $exon (@exonhints){
            my $start = $exon->[0];
            my $end = $exon->[1];
            my $strand = $exon->[2];
            # increase $startidx until the exonpart does not start to the left of start
            while ($startidx <= $#exonparthints && @exonparthints[$startidx]->[0] < $start){
                $startidx++;
            }
            $curidx = $startidx;
            while ($curidx <= $#exonparthints && @exonparthints[$curidx]->[0] <= $end){
                if (@exonparthints[$curidx]->[0] >= $start && @exonparthints[$curidx]->[1] <= $end 
                 && @exonparthints[$curidx]->[2] eq $strand ) {
                    #redundant, delete it
                    #print "deleting " , (join " ", @{$exonparthints[$curidx]}), " as it is contained in " , (join " ", @{$exon}), "\n";
                    splice @exonparthints, $curidx, 1;
                } else {
                    $curidx++;
                }
            }
        }
    }

  #? add scores from gff for exonpart,exon ??
  
   my %hintsets=();
     
  if($isgenemodel) {
     # $hintsets{"exon"}= \@exonhints; # drop for ep
     $hintsets{"ep"}= \@exonparthints; # abbrev for exonpart
     $hintsets{"cp"}= \@CDSparthints; # was CDS; cp better results
     $hintsets{"ip"}= \@intronhints; # changed to ip, better result
     $hintsets{"irpart"}= \@irparthints;
     if( $cdsOn ) {
     $hintsets{"start"}= \@starthints; # not useful
     $hintsets{"stop"}= \@stophints;   # not useful
     }
     
   } elsif($isprotalign) {
     # $hintsets{"CDS"}= \@exonhints; #? drop this one?  blastp (part only) vs genewise/fasty good cds
     $hintsets{"cp"}= \@exonparthints; # abbrev for exonpart
     # $hintsets{"irpart"}= \@irparthints if $isprotfull; ## not working
     $hintsets{"ip"}= \@intronhints if $isprotfull; ## try this
     ## can we use prot w/ separated gene spans to hint intron, intergene ?

   } else { # EST/mRNA data default
     # $hintsets{"exon"}= \@exonhints; # drop?
     $hintsets{"ep"}= \@exonparthints; # abbrev for exonpart
     $hintsets{"ip"}= \@intronhints; # was intron
   }
   
   if ($ssOn) {
     $hintsets{"dss"}= \@dsshints;
     $hintsets{"ass"}= \@asshints;
     }
   
  foreach my $htype (sort keys %hintsets) {
    foreach my $hr (@{$hintsets{$htype}}) {
      my $score= $hr->[4] || 0;
      my $psrc= "$prgsrc_prefix.". ($hr->[5] || "");
      #my $psrc= $prgsrc;
      print HINTS join("\t", $oldtargetname, $psrc, $htype,
        $hr->[0], $hr->[1], $score, $hr->[2],".","grp=$hr->[3];pri=$priority;src=$source"),"\n";
    } 
   }
       
    # delete all hints as the new sequence starts
    @dsshints = ();
    @asshints = ();
    @exonhints = ();
    @exonparthints = ();
    @intronhints = ();
    
  @starthints= @stophints= @tsshints= @ttshints=  @intronparthints=
    @CDShints= @CDSparthints= @UTRparthints= @irparthints= ();

}

#
# addExonpartHint(hintref)
# search in the list of existing exonpart hints for an including or included one
# if no such hint exists, sort the parameter hint into this list
# if there is a stronger one, do nothing, if there are weaker ones, replace them with this one
#

sub addExonpartHint {
    my $href = shift;
    my $begin = $href->[0];
    my $end = $href->[1];
    my $strand = $href->[2];
    my $k;
    #print (join " ", @{$href});
    #print "\n";
    
    my $redundant = 0;
    if ($remove_redundant) {
        # check whether the exonpart hint is contained in one of the exon hints.
        $k = $#exonhints;
        # check the list of previous exonpart hints
        if ($#exonparthints>=0) {
            $k = $#exonparthints;
            while ($k>=0 && $exonparthints[$k]->[0]> $begin - 10000 && !$redundant) { #assume exonpart hints are less than 10000bp
                if ($exonparthints[$k]->[0]<= $begin && $exonparthints[$k]->[1] >= $end && $exonparthints[$k]->[2] eq $strand){
                    #print "found including hint: ", (join " ", @{$exonparthints[$k]}), "\n";
                    $redundant=1;
                } elsif ($exonparthints[$k]->[0] >= $begin && $exonparthints[$k]->[1] <= $end && $exonparthints[$k]->[2] eq $strand){
                    #print "found included hint: ", (join " ", @{$exonparthints[$k]}), " delete it now.\n";
                    splice @exonparthints, $k, 1; #delete k-th element
                }
                $k--;
            }
        }
    }
    if (!$redundant) {
        #insert hint at the right position
        #print "found no redundant hint\n";     
        $k = $#exonparthints;
        if ($remove_redundant) {
            while ($k>=0 && $exonparthints[$k]->[0] > $begin) {
                $k--;
            }
        }
        my @temparray = ($href);
        if ($k == $#exonparthints) {
            #print "insert at end\n";
            push @exonparthints, @temparray;
        } else {
            #print "*** splicing list ***\n";
            splice (@exonparthints, $k+1, 0, @temparray);
        }
    }

}

#
# addSSHint(hintref)  == addSignalHint()
# add the hint if it is not already there

sub addSignalHint {
    my $hintlistref = shift;
    my $href = shift;
    my $begin = $href->[0];
    my $strand = $href->[2];
    #print (join " ", @{$href});

    # add it if the same hint does not exist already
    if (@{$hintlistref}<1 || !$remove_redundant) {
        push @{$hintlistref}, $href;
    } else {
        my $k = @{$hintlistref}-1;
        while ($k>=0 && $hintlistref->[$k]->[0]>= $begin) {
            $k--;
        }
        my @temparray = ($href);
        if (!(($k+1 <= @{$hintlistref}-1 && $hintlistref->[$k+1]->[0]== $begin && $hintlistref->[$k+1]->[2] eq $strand) ||
            ($k+2 <= @{$hintlistref}-1 && $hintlistref->[$k+2]->[0]== $begin && $hintlistref->[$k+2]->[2] eq $strand))) {
            # hint does not previously exist, insert it
            if ($k== @{$hintlistref}-1) {
                push @{$hintlistref}, @temparray;
            } else {
                splice (@{$hintlistref}, $k+1, 0, @temparray);
            }
        }
    }

}

# *addSignalHint = \&addSSHint;

#
# addIntervalHint(hintref)
# for exon and intron hints (not exonpart)
# add the hint if it is not already there

sub addIntervalHint {
    my $hintlistref = shift;
    my $href = shift;
    my $begin = $href->[0];
    my $end = $href->[1];
    my $strand = $href->[2];
    #print (join " ", @{$href});
    #print "\n";

    # shortcut to add all hints, regardless whether they are mutiples
    #push @{$hintlistref}, $href;
    #return;

    # add it if the same hint does not exist already
    if (@{$hintlistref}<1) {
        push @{$hintlistref}, $href;
    } else {
        my $k = @{$hintlistref}-1;
        if ($remove_redundant) {
            while ($k>=0 && ($hintlistref->[$k]->[0]> $begin || ($hintlistref->[$k]->[0]== $begin && $hintlistref->[$k]->[1]>= $end))) {
                $k--;
            }
        }
        my @temparray = ($href);
        if ( !$remove_redundant || 
            !(($k+1 <= @{$hintlistref}-1 && $hintlistref->[$k+1]->[0]== $begin && $hintlistref->[$k+1]->[1]== $end && $hintlistref->[$k+1]->[2] eq $strand) ||
            ($k+2 <= @{$hintlistref}-1 && $hintlistref->[$k+2]->[0]== $begin && $hintlistref->[$k+2]->[1]== $end && $hintlistref->[$k+2]->[2] eq $strand))
            ) {
            # hint does not previously exist, insert it
            if ($k== @{$hintlistref}-1) {
                push @{$hintlistref}, @temparray;
            } else {
                splice (@{$hintlistref}, $k+1, 0, @temparray);
            }
        }
    }

}

