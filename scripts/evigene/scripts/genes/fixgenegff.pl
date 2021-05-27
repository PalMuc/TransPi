#!/usr/bin/env perl
# fixgenegff.pl

=item about

  correct genome-mapped gene gff problems
    .. cut long unsupported introns, esp utr
    .. cut mixed strands across exons, esp utr, from joined genes
    .. cds add phasing / add cds
    .. local-split genes (C1 form), drop part 1 or 2 if possible, utr-only

=item usage

    fixgenegff.pl -actions splits,badutr,cdsupdate -aasizes genes.aaqual [-introns ??] \
       genes.gff > genesfixed.gff

=item  inputs

   genes.gff from various (gmap, splign, evg-blast2gff, ..)
   aaqual table with cds-spans, or from evg-mRNA attr
   intronok validated table, ? or introns.gff
   option: chr.fasta, for cds2prot tests
   
=item updates

    * add cds2prot tests, given chr.fasta
    -- new action from cdsupdate ? cds2prot
    -- before & after makeCDS? (i.e. check input CDS if there)
    -- check both fwd/rev for some (?all), esp. any antisense marked, or valintron != cdsor cases
    
=item see also genefindcds2

  evigene/scripts/genefindcds2.pl
  
  has similar gene.gff corrections, updates using chrasm.fasta
  .. cds trim/extend, 
  .. split gene handling
  .. complex/messy as heck
      
  sooner not later:
    ** MERGE with genefindcds2 : it is really doing what cdsupdate/cdscheck here wants to do
    .. cds-offset calcs here are full of problems because mrna x chr align is often mixed up with
        indels, base changes, gaps, etc making mrna cds-offset wrong for chr-align cdsoffset
    .. flipstrand for antisense mappings is important, should flip sense=- to + strand to
    compare w/ mrna aatrans, unless have homology/other saying fwdsense protein is accurate.  
          
=cut

use FindBin; 
use lib ("$FindBin::Bin", "$FindBin::Bin/../../lib/"); # ugh: evigene/lib/Bio/ from evigene/scripts/genes/this.pl

use strict;
# use warnings;
use Getopt::Long;
# use Bio::DB::Fasta; # now as get_dna() require  Bio::DB::Fasta, so can use this w/o bioperl  

my $debug= $ENV{debug}||0;
# my $PRINTALL=1; # $ENV{printall}||0;
   
my $mrnatypes='mRNA';
my $exontypes='exon|CDS';
my $passtypes="";
my $actions="splits,badutr,cdsupdate"; # add cdstrans 

my $OUTH= *STDOUT; 
my $MININTRON= 10; # NCBI sez  SEQ_FEAT.ShortIntron and SEQ_FEAT.AbuttingIntervals, ie exons overlap
  # SEQ_FEAT.ShortIntron  Introns should be at least 10 nt long
my $MINEXON= 20; # SEQ_FEAT.ShortExon  Internal coding region exon is too short .. but what is "too"?
my($intronoktab,$genesizes,$outfile)=('') x 9;
my $LONGOK_MAX=9999; my $LONGOK_MAXin=2; 
my $MAX_UTRX= 3; # 9999; # i.e. unlimited ??
my $ANTISENSE_INTRON_IS_ERROR= 0; # option now, as many daphnia cds appear valid but reverse of rna-intron strand
my $DIFF_CDS=0;
my($chrfasta, $fasta_db); # for prot trans


my $USAGE =<<"USAGE";
usage: fixgenegff.pl -actions $actions -aasize genes.aaqual  \\
   -intronok genes.intronok genes.gff > genesfixed.gff
   -chrfasta chr.fasta for -act cdstrans
USAGE

my $optok= GetOptions( 
  "output=s" => \$outfile, # not ready; add opt tabular output of changes, same as mRNA annots
  "sizes|aasizes=s" => \$genesizes,
  "chrfasta=s" => \$chrfasta,
  "actions=s", \$actions,
  "intronoktab=s" => \$intronoktab,
  "MININTRON=i", \$MININTRON,  
  "MINEXON=i", \$MINEXON,
  "longintronmax=i", \$LONGOK_MAX, "longcountintrons=i", \$LONGOK_MAXin,
  "maxutrx=i", \$MAX_UTRX,  
  "antierror|errantisense!" => \$ANTISENSE_INTRON_IS_ERROR, # default off, add OFFBY1_IS_ERR, use same flag?
  "DIFF_CDS!" => \$DIFF_CDS,
  "debug!" => \$debug,
  );
die $USAGE unless($optok);

my(@changelog); # debug output?
$MAX_UTRX= 9999 if($MAX_UTRX<1);  
# handle $outfile option.
 
#off# my $ONLYCHANGES=$debug; # debug
#off#    $ONLYCHANGES=0 if($PRINTALL);
#----------------------------------------


sub MAINstub {}

if($actions =~ /cdstrans/) { ## requires chrfasta
  unless($chrfasta and -f $chrfasta) {
    die "option -action=cdstrans requires -chrfasta /path/to/chr.fasta\n", $USAGE;
  }
}  

my($nsizes,$aasizeh,$trsizeh,$cdspanh,$aaqualh)
    = ($genesizes) ? readSizes($genesizes) : (0,0,0,0);   

my($nintrons, $inokerr, $introns)= ($intronoktab) ? readIntronOk($intronoktab) : (0);

# my ($ok,$input);
# my $inh= *STDIN; # or perl (<>) == @ARGV ?
# if($input) {
#   $ok = ($input =~ /.gz$/) ? open($inh,"gunzip -c $input |") 
#         : ($input =~ /^(stdin|-)/) ? $inh= *STDIN
#         : open($inh,$input);
#   die "bad -input=$input" unless($ok);
# }

my($ng,$nx,$nfix)= filter_gff(); # ($inh)

# debug output/alt tabular output @changelog, other info? 
# maybe append @changelog (# geneid \t changes) to end of gff?
if($debug) {
  print $OUTH "#fixgenegff: ng=$ng, nx=$nx, nfix=$nfix\n"; 
  for my $cl (@changelog) { print $OUTH "#fixg\t$cl\n"; }
}

#----------------------

sub splitPart {
  my($xattr)=@_;
  my($xid) = $xattr =~ m/(?:ID|Parent)=([^;\s]+)/;  
  ## which split key has precedence? Split=\d or ID_C\d ?
  #x my($spl)= ($xattr =~ m/;Split=(\d+)/)? $1 : ($xid =~ /_C(\d+)$/)? $1 : 0;
  my($spl)= ($xid =~ /_C(\d+)$/)? $1 : ($xattr =~ m/;Split=(\d+)/)? $1 : 0;
  return($spl,$xid);
}

# fixme for "." strand compare
sub _eqstrand { return($_[0] eq $_[1] || $_[0] eq '.' || $_[1] eq '.' )?1:0; }
# convert to array handling? _max(@list) 
sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }

sub _sortgene  
{
  # ($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid)
  my($ta,$tb)= map{ (m/gene/)?1:(m/mRNA/)?2:(m/CDS|exon/)?3:4; } ($a->[2],$b->[2]);
  return ($ta <=> $tb)
      || ($a->[0] cmp $b->[0]) # ref : should all be same for gene record
      || ($a->[3] <=> $b->[3]) # begin
      || ($a->[4] <=> $b->[4]) # end
      || ($b->[2] cmp $a->[2]) # type: exon > CDS
      ;
}

sub _sortSplitgene {
  # ($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid)
  #FIXME? for split genes, Split=1,2,3/3 order by split-part? all mRNA at top or not?
  # need to apply split sort to each part: mrna, exon,cds .. better way?
  my($spla,$aid)= splitPart($a->[8]);
  my($splb,$bid)= splitPart($b->[8]);
  
  my($ta,$tb)= map{ (m/gene/)?1:(m/mRNA/)?2:(m/CDS|exon/)?3:4; } ($a->[2],$b->[2]);
  return ($spla <=> $splb) || ($ta <=> $tb)
      || ($a->[0] cmp $b->[0]) # ref : should all be same for gene record
      || ($a->[3] <=> $b->[3]) # begin
      || ($a->[4] <=> $b->[4]) # end
      || ($b->[2] cmp $a->[2]) # type: exon > CDS
      ;
}

sub __revSplitGene {  # rev-sort -strand exons per part
  # my($aid,$bid)= map{ m/(?:Parent|ID)=([^;\s]+)/; $1; } ($a->[8],$b->[8]);
  # my($spla)=($a->[8] =~ m/;Split=(\d+)/)? $1 : ($aid =~ /_C(\d+)$/)? $1 : 0;
  # my($splb)=($b->[8] =~ m/;Split=(\d+)/)? $1 : ($bid =~ /_C(\d+)$/)? $1 : 0;
  my($spla,$aid)= splitPart($a->[8]);
  my($splb,$bid)= splitPart($b->[8]);
  
  my $ora= ($a->[6] eq '-' or $a->[6] eq -1)?-1:1;
  my $orb= ($b->[6] eq '-' or $b->[6] eq -1)?-1:1;
  return ($spla <=> $splb) # split part
    || ($a->[0] cmp $b->[0]) # ref
    || ($ora * $a->[3] <=> $orb * $b->[3]); # strand*begin
}


sub cds_span {
  my($geneid,$mrnaattr)= @_;
  # my($cdsb,$cdse)=(0,0);
  my $cdsor= 0; # want Target orient, mRNA is +1 always, but antisense, or cDNA targ can be -1
  my $cdsoff= $cdspanh->{$geneid} if( ref $cdspanh );
  unless($cdsoff) { ($cdsoff)= $mrnaattr =~ m/offs=(\d+.\d+)/; }
  my($cdsb,$cdse)= ($cdsoff and $cdsoff=~m/(\d+).(\d+)/) ? ($1,$2) : (0,0);
  if($cdse>0 and $cdse < $cdsb) { ($cdsb,$cdse)= ($cdse,$cdsb); $cdsor= -1; }
  elsif($cdsoff =~ m/(\d+).(\d+):(.)/) { my $cor=$3; $cdsor=($cor eq '-')?-1:+1; } # :. is common, same as :+
  else { $cdsor= +1; }
  return ($cdsb,$cdse,$cdsor);
}


sub filter_gff
{
  my($inh)= @_;
  my ($ng,$nx,$nr,$nfix,$nhit,$errord)= (0) x 10;
  my $nocomm= 1; ##($actid == ACT_NULL) ? 1 : 0;
  my @generec=();
  my @geneother=();
  my $geneid=""; my $lgeneid="";
  my $generow=undef;
  
  # while(<$inh>) 
  while(<>) { # perl special, STDIN or @ARGV
    unless(/^\w/){ print unless($nocomm or /^(#n |$)/); next; }
    my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr)= split"\t";
    $nr++; chomp($tattr);
   
    my($gid,$pid); 
    if($tattr =~ m/\bID=([^;]+)/) {  $gid=$1; }
    if($tattr =~ m/\bParent=([^;]+)/) {  $pid=$1; 
      $gid=$pid unless($gid and ($typ =~ /^($mrnatypes)$/)); 
      }
    unless(defined $gid) { $gid = "N".$ng; }

    # my $issplit= ($gid =~ m/_C(\d+)/)?$1 : (m/;Split=(\d+)/)?$1:0;
    my($issplit,$gidsp)= splitPart($tattr);

    if($issplit) {
      $gid =~ s/_C\d+//; #?? maybe not here, want split-part ID to process each exon subset
      $pid =~ s/_C\d+//; 
    }
    
    my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid,$issplit,$gidsp]; 

    if($typ =~ /^gene$/) { $generow= $rloc; }
    elsif($typ =~ /^($mrnatypes)$/) {  
      # allow gene and mRNA types .. and/or keep other types in generec
      
      if($gid eq $geneid and $issplit) {
        push @generec, $rloc; $ng++; # check Parent == $geneid ?
      
      } else {
      
        $nfix += fixgene($actions, $geneid, \@generec, \@geneother) if(@generec);
        # ($changes, $changelog, $rowsput)=  fixgene($actions, $geneid, \@generec, \@geneother)
        
        $ng++;
        $geneid= $gid;  # parse for gene vs tr id/alttr num ?
        @generec = ($rloc); # maybe best store as array of [gffcols], sorted by type-loc
        @geneother= (); 
        # if(ref $generow){ unshift @generec, $generow; $generowGlobal= $generow; }
        $generow=undef;
      }
      
    } elsif($typ =~ /^($exontypes)$/) {
      if($pid ne $geneid) { warn "#ERR: Out-of-order GFF $typ:$pid in mRNA:$geneid\n"; $errord++; next; }
      push @generec, $rloc; $nx++; # check Parent == $geneid ?
      
    } elsif( ! $passtypes or "$typ.$src" =~ m/$passtypes/ ) {
      push @geneother, $rloc; # check Parent == $geneid ?
    }
    
  }
  
  $nfix += fixgene($actions, $geneid, \@generec, \@geneother) if(@generec);
  return ($ng,$nx,$nfix);
}


sub putgene {
  my($generec, $geneother, $nchanges) = @_;
  $geneother ||= [];  my $nput=0;
  # note generec contains gene, mrna, exons, cds, in orig order; updated locs ..
  # my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid]; 
  foreach my $rloc (@$generec, @$geneother) {
    next unless(ref $rloc); # gene row may be undef
    print $OUTH join("\t", @$rloc[0..8])."\n"; $nput++;
  }
  return $nput; # or nput?
}


use constant DEFER_MRNA_FIXSPAN => 1; # leave this update to last     

sub mrna_fixspan {
  my($mrna,$exons)= @_;
  my $changed=0;
  ## splitgene fix: check ref same.
  
  # FIXME lost 1st gene row before mrna; FIXME change gene span w/ mrna span changes
  # ALSO check/fix strand, mixed strand exons: most common.
  ##my($mrna) = grep{ $_->[2] eq "mRNA" } @generecfix;
  my($mref,$oldstart,$oldstop,$oldor)= ($mrna->[0],$mrna->[3],$mrna->[4],$mrna->[6]);
  my($newstart,$newstop,$nrev,$nfwd)= (0) x 9;
  my $nxok= @$exons;
  foreach my $ex (@$exons) {
    my($xr,$xb,$xe,$xo)= @{$ex}[0,3,4,6];
    next if($xr ne $mref); #?? bug
    if($xo eq '-') { $nrev++; } elsif($xo eq '+') { $nfwd++; }
    ## bug split part mrna get all parts span on same chr .. match Split=n from mrna/exon?
    $newstart= $xb if($newstart == 0 or $xb < $newstart);
    $newstop = $xe if($newstop == 0 or $xe > $newstop);
    }
  
  my $newor= ($nrev > $nfwd)?'-':'+';
  unless($newstop==0 or ($newstart == $oldstart and $newstop == $oldstop and $newor eq $oldor)) { 
    # ** Record changes in mrna attr
    $changed++; 
    $mrna->[3] = $newstart; $mrna->[4] = $newstop; 
    $mrna->[6] = $newor; # unless($newor eq $oldor);
    $mrna->[8] =~ s/$/;fixoldspan=$oldstart-$oldstop:$oldor/;
    $mrna->[8] =~ s/;nexon=\d+/;nexon=$nxok/; # add if missing?
    my $hasmix= ($mrna->[8]=~/strandmix=/)?1:0;
    if($nrev>0 and $nfwd>0) {  
      $mrna->[8] =~ s,$,;strandmix=+$nfwd/-$nrev, unless($hasmix);
    } elsif($hasmix) {
      $mrna->[8] =~ s,;strandmix=[^;\s]+,,;
    }
    }
  return $changed;
}

sub fixgene
{
  my($actions, $geneid, $generecIN, $geneother)= @_;
  my($changes, $didchange)=(0,0);
  my $changelog="";
  
  # my $generec= $generecIN;
  my @mrna = grep{ $_->[2] eq "mRNA" } @$generecIN; # fixme issplit/Splitgene has more mrna parts
  my $issplit= (@mrna>1)? @mrna : 0;
  my @generec;
  if($issplit) { @generec= sort _sortSplitgene  @$generecIN;  } 
  else {  @generec = sort _sortgene  @$generecIN;  }
  my $generec= \@generec;
  
  # my $actions="splits,badutr,cdsupdate";
  ## can ignore geneother here ; subs dont need
  
  if($actions =~ /splits/) {
  ($didchange, $generec, $geneother)= 
    fixgene_split($geneid, $generec, $geneother);
  #Not for each gene: warn "#fix splits: $didchange\n" if($debug);
  $changelog.="fixsplit$didchange," if($didchange);
  $changes+= $didchange;
  }
  
  if($actions =~ /badutr/) {
  ($didchange, $generec, $geneother)= 
    fixgene_cutbadutr($geneid, $generec, $geneother);
  #Not for each gene: warn "#fix cutbadutr: $didchange\n" if($debug);
  $changelog.="cutbadutr$didchange," if($didchange);
  $changes+= $didchange;
  }
  
  if($actions =~ /cdsupdate/) {
  ($didchange, $generec, $geneother)= 
    fixgene_cdsupdate($geneid, $generec, $geneother);
  #Not for each gene: warn "#fix cdsupdate: $didchange\n" if($debug);
  $changelog.="cdsupdate$didchange," if($didchange);
  $changes+= $didchange;
  }

  my($ctransflag, $phaseb)=(0,0);
  if($actions =~ /cdstrans/) { # this is test/compare only, not changer?; before/after cdsupdate?
    ($didchange, $generec, $ctransflag, $phaseb)= 
      fixgene_cdstrans($geneid, $generec);
    $changelog.="cdstrans$didchange:$ctransflag," if($didchange); #? no changes here? mrna.attr flag added
    # $changes+= $didchange;
  }
  
  if($changes) { #split?
    my @mrnanew = grep{ $_->[2] eq "mRNA" } @$generec; # must have
    my @exons = grep{ $_->[2] eq "exon" } @$generec;
    if(@mrnanew > 1) { # FIXME.. need to separate exons/part
    
    } else {
      $changes += mrna_fixspan($mrnanew[0], \@exons);
    }
  }
  
  $changelog ||= "none";
  push @changelog, "$geneid\t$changelog"; # if($debug) 
  warn "#fixg\t$geneid\t$changelog\n" if($debug>1);
  
  my $rowsput= putgene($generec, $geneother, $changes);
  return ($changes) ? 1 : 0;
}


sub overexons {
  my($ind,$exons)= @_;
  my($ir,$ib,$ie,$io)=split/:/,$ind; 
  my($ib1,$ib2,$ie1,$ie2)= ($ib-1,$ib-2,$ie+1,$ie+2); # exon splice is +/-1 of intron ends, allow off by 1?
  
  #? check for exon pair matching intron instead of splice point only?
  # this way we *should* get 2 exons, for each end of intron .. if not?
  my @xov=();
  foreach my $rloc (@$exons) {
    # my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid]; 
    my($xr,$xt,$xb,$xe,$xo,$gid)= @{$rloc}[0,2,3,4,6,9];
    next unless($xr eq $ir);
    # next unless(_eqstrand($xo,$io)); # ?? antisense errors skip this
    my $ov=0;
    $ov |=1 if($xb == $ie1 or $xb == $ie2 or $xb == $ie);  
    $ov |=2 if($xe == $ib1 or $xe == $ib2 or $xe == $ib);  
  
    push @xov, $rloc if($ov);  
  }
  
  my $nov=@xov;
  if($nov >= 2) { # >2 is mistake?, drop all but 2? iterate i=0.., j=i+1 exon pairs?
    for( my $i=0; $i<$nov-1; $i++) { my $j=$i+1;
      my($xb,$xe,$xo)= @{$xov[$i]}[3,4,6];
      my($yb,$ye,$yo)= @{$xov[$j]}[3,4,6];
      my $ov=0;
if(0) { $ov=1; } else { # debug this
      $ov = 2 if(($xe == $ib1 or $xe == $ib2 or $xe == $ib) and ($yb == $ie1 or $yb == $ie2 or $yb == $ie));  # should be this way
      $ov = 1 if(($xb == $ie1 or $xb == $ie2 or $xb == $ie) and ($ye == $ib1 or $ye == $ib2 or $ye == $ib));  #? rev intron
}      
      return @xov[$i,$j] if($ov); # otherwise what?
      }
  }
  return ();
}


sub fixgene_cutbadutr # mixed strands, long-invalid-introns, 
{
  my($geneid, $generec, $geneother)= @_;
  my($changed, $badspan, $cdscut)=(0) x 9; 

  my (@newgenerec, @oldrecparts);
  my @mrna = grep{ $_->[2] eq "mRNA" } @$generec; # must have
  # my $mrna = $mrna[0];
  # my @cds  = grep{ $_->[2] eq "CDS"  } @$generec;
  # my @exons= grep{ $_->[2] eq "exon" } @$generec;
  
  my $issplit= (@mrna > 1)? @mrna : 0;
  
  ## FIXME: splitgenes .. @mrna loop, each w/ uniq ID
  for my $mrna (@mrna) {
    my $didaddpart=0;
    my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid,$issplit1,$gidsp)= @$mrna;
    my @cds  = grep{ $_->[2] eq "CDS"  and $_->[8] =~ m/=$gidsp/ } @$generec;
    my @exons= grep{ $_->[2] eq "exon" and $_->[8] =~ m/=$gidsp/ } @$generec;
    
    my @inerrexons=();
    #?? does intronoktab have geneid or gid-splitpart? 2nd?
    if(my $inok= $inokerr->{$gidsp}) { # flags: CDSok,UTRok,longerr,antierr
      if($inok =~ /err/) {
            
        # figure out where errs are, if in utr, cut
        # $introns{$gidsp}{$ind}= join"\t",$iflag,$inw,$iv,$xtype;
        # $ind == scaffold_116:273041:273111:+
        my @in= sort keys %{$introns->{$gidsp}};
        my @inerr= grep{ $introns->{$gidsp}{$_} =~ /err/ } @in; 
        
        # ** Make antisense optional error, not always, not default?
        # UPD: antisense err should be ignored but for UTR ne CDS orient, 
        #   .. introns strand appears unreliable for CDS strand (what we want here)
        
        # $ind == scaffold_116:273041:273111:+
        for my $ind (@inerr) {
          my($ir,$ib,$ie,$io)=split/:/,$ind; 
          my($iflag)= split"\t",$introns->{$gidsp}{$ind};
          my @xo= overexons($ind, \@exons); # 2 exons per intron?
          #?? not both exons of intron, but farthest one from cds-middle, need to know utr/cds here?
          push @inerrexons, @xo if(@xo);
        }
      # my $nerr=@inerrexons;
      # warn "#xinerr\t$geneid\txerr=$nerr\n" if($debug and $nerr); # not seen..now seen
      }
    }
  
    my $ncds= @cds; #may be missing, either not added, or mapped exons may be UTR-only : cant cut those cases
    my($cdsb,$cdse,$cdsor)= cds_span($geneid,$tattr); # geneid or gidsp ?
    if($cdse > 0) {
      # my($revor)= $cdsor; # ($to eq '-')?-1:0; # ?? target/cds-or or genome-or? .. not used in cututrx()
      my($ncut,$newexons)= cututrx($mrna,$cdsb,$cdse,$ncds,$cdsor,\@exons, \@inerrexons);
      
      if($ncut) {
        $changed += $ncut; #??
        # my @newgenerec=( $mrna, @$newexons, @cds);
        #x return ( $changed, \@newgenerec, $geneother); ## 0 == no change ??
        push @newgenerec,  $mrna, @$newexons, @cds; 
        $didaddpart=1;
      }
    }
    
    unless($didaddpart) {
      push @oldrecparts, $mrna, @exons, @cds;
    }  
    
  } # @mrna parts
  
  if($changed) {
    push @newgenerec, @oldrecparts if(@oldrecparts);
    if($issplit) {  @newgenerec= sort _sortSplitgene @newgenerec;  } 
    else {  @newgenerec = sort _sortgene  @newgenerec;  }
    return ($changed , \@newgenerec, $geneother); ## 0 == no change ??
  }
  
  return ($changed , $generec, $geneother); ## 0 == no change ??
}


=item intron_short_fix
  
  * update to fix_intron_tooshortlong ?
  from evigene/scripts/genefindcds2.pl
  remove too-short introns, from gmap/gsplign/other, to MININTRON, changing exons as needed
  call this BEFORE bestorf so dont have to redo.
  
=cut

sub intron_short_fix  {  
  my($generec,$exontype) = @_; 
  #?? return($changed,$generec,$infixnote);  
  $exontype||="exon";
  my($infixnote,$changed)=("",0);
  
  #** FIXME for isssplit, generec has all parts, ref diff
  
  my @generec= sort _sortgene @$generec; # sorts by genome loc, not stranded
  my @exongff= grep{ $_->[2] eq $exontype } @generec; # must be _sortgene
  if(not @exongff and $exontype eq "exon") {
    $exontype="CDS";
    @exongff= grep{ $_->[2] eq $exontype } @generec; # must be _sortgene
  }
  
  my @exoninfix=();
  my($ix,$lr,$lb,$lt,$le,$lo,$lid,$lloc)= (0) x 10;  $le=-99999;
  foreach my $rloc (@exongff) {
    #new# my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid,$inb,$ine]; 
    $ix++; my($xr,$xt,$xb,$xe,$xo,$gid)= @{$rloc}[0,2,3,4,6,9];

    # my $over  = ($xb <= $le && $xe >= $lb) ? 1 : 0; 
    # my $near  = (abs($xb - 1 - $le) < $MININTRON)?1:0; # abs? sorted should put xb >= le? but overlaps may be xb<le
    # my $inside= ($over and $tb <= $lb && $te >= $le) ? 1 : 0; ## intron inside  
    # $over=0 if($stranded and ($to =~ /[+-]/ and $lo =~ /[+-]/ and $to ne $lo)); 
    # $inside=0 if($stranded and ($to =~ /[+-]/ and $lo =~ /[+-]/ and $to ne $lo));

    my $indist= ($lr ne $xr)?99999:($xb - 1 - $le); # dont care about strand, other, just intron distance < min
    ## ^^ Off-by+1 bug, was ($xb - $le); intron size = 1 + ($xb-1) - ($le+1) == ($xb - 1 - $le)
    #?? check/fix short exon sizes?  xe-xb < MINEXONW
    if(($lr eq $xr) and $le and $indist < $MININTRON) { # abs? sorted should put xb >= le? but overlaps may be xb<le
      $changed++;
      if($xe > $le) { $le=$xe; $lloc->[4]= $le; } # join 2 exons.
      $infixnote .= "xcut$ix:$indist,"
      #NOT# $lloc= $rloc; ($lr,$lt,$lb,$le,$lo,$lid)=($xr,$xt,$xb,$xe,$xo,$gid);
   } else { 
      push @exoninfix,$rloc; 
      $lloc= $rloc; ($lr,$lt,$lb,$le,$lo,$lid)=($xr,$xt,$xb,$xe,$xo,$gid);
    }
    
  }
  
  if($changed) { 
    my @oldexon= map { my @xnew= @$_; \@xnew; } @exongff; # clone all
    my @generecfix= grep{ $_->[2] ne $exontype } @$generec;
    push @generecfix, @exoninfix;
    if($debug>1) { foreach my $ex (@oldexon) {
      $ex->[2] = "oldexon"; $ex->[0] =~ s/^/#/; push @generecfix, $ex;
    } }
    #?? any newstart/stop changes ?? not for squeezing overlap exons?
    $generec= \@generecfix;
    if($exontype eq "exon") { # also fix CDS ..
      my($xchanged, $xgenerec, $xfixnote)= intron_short_fix($generec,"CDS"); 
      $generec=$xgenerec if($xchanged);
    }
    $infixnote=~s/,$//;
  }

  return($changed,$generec,$infixnote);  
}

sub cututrx {
  my($mrna,$cdsb,$cdse,$ncds,$revor,$xons, $inerrexons)= @_; 
  
  my $maxutrx= $MAX_UTRX; #2; # FIXME option, cut utr for errors: longint, splicemix..
  my $nx= @$xons;
  ## cant trust $ncds, may not have added CDS at this point, trust if $ncds > 0 ?
  
  #?? cant rely on inerrexons in UTR, case of CDS antisense, then UTR wrong strand .. chop those
  my %xerrs=(); my $nxerrs=0;
  if($inerrexons and ref($inerrexons)) {
    # inerrexons == exons matched to error introns, cut these if outside cds
    #?? not both exons of intron, but farthest one from cds-middle, need to know utr/cds here?
    for my $ex (@$inerrexons) { $nxerrs++; $xerrs{$ex}="inerr1"; } # flag err type?
  }

  my(@cx,@uxb,@uxe);
  my($xor,$ltb,$ncfwd,$ncrev,$nufwd,$nurev)=(0) x 9;
  for my $ex (@$xons) {
    my($xr,$xs,$xt,$xb,$xe,$xp,$xo,$tph,$xat)= @$ex;
    my($tid,$tb,$te)= $xat =~ m/(?:Target|trg)=(\S+)\s(\d+)\s(\d+)/;
    next unless($te); # warn/error
    
    # UPD: use introns valid splice to override xo
    # my($inor)= valintron('strand', $gidsp, $xr, $xb, $xe, $xo);
    # $xo= $inor if($inor and $inor ne $xo);
    
    my $cdx=($tb <= $cdse and $te >= $cdsb)? 1 :($te<$cdsb)?-1 :0; 
    # my $haserr= $xerrs{$ex}||0;
    if($cdx==1) { push @cx,$ex ; if($xo eq '-') { $ncrev++; } elsif($xo eq '+') { $ncfwd++; } } 
    else {
      if($xo eq '-') { $nurev++; } elsif($xo eq '+') { $nufwd++; }
      if($cdx<0) { push @uxb,$ex; } else { push @uxe, $ex; }
    }
    if($ltb) { $xor= ($ltb <= $tb)? 1 : -1; }
    $ltb=$tb;
  }
  
  my $ncdsx= @cx; $ncds= $ncdsx if($ncdsx > $ncds);
  my $cor= ($ncrev > $ncfwd)?'-':($ncfwd>0)?'+':0;
  my $uor= ($nurev and $nufwd)?'m':($nufwd>0)?'+':($nurev>0)?'-':0;
  my $strandmix=($cor ne $uor)?"ormix":"";  
  $nxerrs++ if($strandmix);

  if($ncdsx < 1) { # possible ncRNA? or utr split-part, should not be here
    $mrna->[8] =~ s/$/;cututrx=fail:nocdsmap/;  
    return(0, $xons) ; # add mRNA flag? or let drop down to below cututrx=fail ?
  }
  return(0, $xons) unless($nxerrs or @uxe > $maxutrx or @uxb > $maxutrx);   

  my $nxp= @cx + @uxb + @uxe;
  if($nxp == $nx) { # bug otherwise?
    my(@uxbok, @uxeok);
    my($bcut,$ecut)=("") x 9;
    my $nxu= @uxb;
    for(my $i=0; $i<$nxu; $i++) {
      my ($k,$k2)= ($xor<0) ? ($i,$i+1) : ($nxu - 1 - $i, $nxu - 2 - $i);
      my $ex= $uxb[$k]; ## k2 is 2nd exon of intronerr pair, cut at k2
      my $err= $xerrs{$ex}||0;
      if($strandmix and $ex->[6] ne $cor) { $err=$strandmix; }
      if($err) { if($xerrs{$uxb[$k2]}) { push @uxbok, $ex; }  $bcut= $err; last;}
      elsif($i >= $maxutrx) { $bcut="maxutr$i"; last; }
      else { push @uxbok, $ex; }
    }
    $nxu= @uxe;
    for(my $i=0; $i<$nxu; $i++) {
      my ($k,$k2)= ($xor>0) ? ($i,$i+1) : ($nxu - 1 - $i, $nxu - 2 - $i);
      my $ex= $uxe[$k];
      my $err= $xerrs{$ex}||0;
      if($strandmix and $ex->[6] ne $cor) { $err=$strandmix; }
      if($err) { if($xerrs{$uxe[$k2]}) { push @uxeok, $ex; }  $ecut= $err; last; } 
      elsif($i >= $maxutrx) { $ecut="maxutr$i"; last; }
      else { push @uxeok, $ex; }
    }
    
    # what of this: cututrx=0,b,einerr1 ; none cut, but set error?
    my @newexons= sort _sortgene (@cx,@uxbok,@uxeok); # ,@cds
    my $nxnew= @newexons;
    my $ncut= $nx - $nxnew; 
    # ?? allow nxnew == 0 for split part only utr, i.e. $nxnew == $ncds == 0 for utr-only part
    # BUT fixgene_split should be handling utr-only parts.. dont add here also
    if($nxnew < 1 or $nxnew < $ncds) {
      $mrna->[8] =~ s/$/;cututrx=fail:$ncut,b$bcut,e$ecut/;  
      return(0,$xons);   
    }
    
    if($bcut or $ecut) { # if($ncut != 0)  # any -ncut adding utrx ? .. add bcut/ecut
      $mrna->[8] =~ s/$/;cututrx=$ncut,b$bcut,e$ecut/; # FIXME better flags
      return($ncut, \@newexons);
    }
  }
  return(0,$xons);   
}


=item inerr: no CDS in exons
  .. should keep mapping but annotate as CDS out of map range ..
  .. this case likely gene join tho, mapped UTR far from unmapped CDS
grep =daplx6ml10dn9anvelvk31Loc26t14 velvdtest.infix.gff     
scaffold_293	velml10d	mRNA	121147	123303	99	+	.	ID=daplx6ml10dn9anvelvk31Loc26t14;Target=daplx6ml10dn9anvelvk31Loc26t14 4561 6510;qlen=6510;cov=30;pid=99;match=1950;
    nexon=4;intralign=4;strandmix=+3/-1;aalen=333,15%,complete-utrbad; <<
    offs=1078-77:-;path=1/14;inhit=m2,2;cututrx=4:fail,b,eormix
scaffold_293	velml10d	exon	121147	121259	95	-	.	Parent=daplx6ml10dn9anvelvk31Loc26t14;Target=daplx6ml10dn9anvelvk31Loc26t14 4561 4673;align=113;introna=0,0
scaffold_293	velml10d	exon	121356	121457	94	+	.	Parent=daplx6ml10dn9anvelvk31Loc26t14;Target=daplx6ml10dn9anvelvk31Loc26t14 4678 4783;align=106;introna=-2,-3
scaffold_293	velml10d	exon	121516	121628	97	+	.	Parent=daplx6ml10dn9anvelvk31Loc26t14;Target=daplx6ml10dn9anvelvk31Loc26t14 4783 4894;align=112;introna=-1,-1
scaffold_293	velml10d	exon	121688	123303	100	+	.	Parent=daplx6ml10dn9anvelvk31Loc26t14;Target=daplx6ml10dn9anvelvk31Loc26t14 4892 6510;align=1619;introna=-2,0
  >> alt path has CDS
scaffold_293	velml10d	mRNA	24714	36205	99	+	.	ID=daplx6ml10dn9anvelvk31Loc26t14_G2;Target=daplx6ml10dn9anvelvk31Loc26t14 1 5033;qlen=6510;cov=77;pid=99;match=5033;
  nexon=14;intralign=14;strandmix=+10/-4;aalen=333,15%,complete-utrbad;offs=1078-77:-;path=2/14
  

=item inerr2: cut all of exons; why? strandmix cases

grep =daplx6ml10dn9anvelvk31Loc285t2 velvdtest.infix.gff | head
scaffold_23	velml10d	mRNA	562931	567522	98	+	.	ID=daplx6ml10dn9anvelvk31Loc285t2;Target=daplx6ml10dn9anvelvk31Loc285t2 1 3707;
  qlen=3731;cov=99;pid=98;match=3707;nexon=10;intralign=10;
  strandmix=+6/-4;aalen=416,33%,complete-utrbad;offs=1340-90:-;inhit=m4,7;
  cututrx=10,binerr1,einerr1
  .. no exons
scaffold_23	velml10d	CDS	563020	563401	98	-	.	Parent=daplx6ml10dn9anvelvk31Loc285t2
scaffold_23	velml10d	CDS	563462	563602	100	-	.	Parent=daplx6ml10dn9anvelvk31Loc285t2
scaffold_23	velml10d	CDS	563661	563981	99	-	.	Parent=daplx6ml10dn9anvelvk31Loc285t2
scaffold_23	velml10d	CDS	564044	564437	99	-	.	Parent=daplx6ml10dn9anvelvk31Loc285t2
    int scaffold_23:564438:564501:+ here
scaffold_23	velml10d	CDS	564502	564509	99	+	.	Parent=daplx6ml10dn9anvelvk31Loc285t2 # bad tiny exon, should cut

grep daplx6ml10dn9anvelvk31Loc285t2 velvdaplx6ml10dn9_inlocjtr_lo60s.inhitok.intronok.tab
daplx6ml10dn9anvelvk31Loc285t2	scaffold_23:563402:563461:-	60	valintron=2056	CDS
daplx6ml10dn9anvelvk31Loc285t2	scaffold_23:563603:563660:-	58	valintron=2018	CDS
daplx6ml10dn9anvelvk31Loc285t2	scaffold_23:563982:564043:-	62	valintron=2780	CDS
daplx6ml10dn9anvelvk31Loc285t2	scaffold_23:564438:564501:+	64	valintron=0,antisense=1322	CDS << wrong strand, but in CDS
daplx6ml10dn9anvelvk31Loc285t2	scaffold_23:565625:565688:+	64	valintron=52	UTR << wrong strand UTR ..
daplx6ml10dn9anvelvk31Loc285t2	scaffold_23:565840:565910:+	71	valintron=108	UTR

..

=item inerr, no cut: need to fix?

scaffold_23	velml10d	mRNA	763878	772525	98	+	.	ID=daplx6ml10dn9anvelvk31Loc1616t38;
  Target=daplx6ml10dn9anvelvk31Loc1616t38 5 4985;qlen=5009;cov=99;pid=98;match=4981;nexon=4;intralign=4;
  >> strandmix=+3/-1;aalen=398,23%,complete-utrbad;offs=213-1409:+;path=1/2;inhit=m2,2;cututrx=0,b,einerr1
scaffold_23	velml10d	exon	763878	765877	99	-	.	Parent=daplx6ml10dn9anvelvk31Loc1616t38;Target=daplx6ml10dn9anvelvk31Loc1616t38 2987 4985;align=1999;introna=0,-1
                                                ^^ bad strand, hits ind scaffold_23:765878:765979:+
                                                .. probably should flip strand for this case
scaffold_23	velml10d	exon	765980	768696	98	+	.	Parent=daplx6ml10dn9anvelvk31Loc1616t38;Target=daplx6ml10dn9anvelvk31Loc1616t38 278 2987;align=2710;introna=0,0
scaffold_23	velml10d	CDS	767558	768696	98	+	.	Parent=daplx6ml10dn9anvelvk31Loc1616t38
scaffold_23	velml10d	CDS	772055	772115	97	+	.	Parent=daplx6ml10dn9anvelvk31Loc1616t38
scaffold_23	velml10d	exon	772055	772146	97	+	.	Parent=daplx6ml10dn9anvelvk31Loc1616t38;Target=daplx6ml10dn9anvelvk31Loc1616t38 182 279;align=98;introna=-2,-4
scaffold_23	velml10d	exon	772348	772525	92	+	.	Parent=daplx6ml10dn9anvelvk31Loc1616t38;Target=daplx6ml10dn9anvelvk31Loc1616t38 5 186;align=182;introna=-1,0

grep daplx6ml10dn9anvelvk31Loc1616t38 velvdaplx6ml10dn9_inlocjtr_lo60s.inhitok.intronok.tab
daplx6ml10dn9anvelvk31Loc1616t38	scaffold_23:765878:765979:+	102	valintron=0,antisense=250	UTR <<
daplx6ml10dn9anvelvk31Loc1616t38	scaffold_23:768697:772054:+	3358	valintron=17	CDS
daplx6ml10dn9anvelvk31Loc1616t38	scaffold_23:772147:772347:+	201	valintron=12	UTR
  
=cut    


sub fixgene_split
{
  my($geneid, $generec, $geneother)= @_;
  my($changed, $badspan,$cdscut)=(0) x 9; 
  
  # my @generec= sort _sortgene @$generec; #? sort by genostart or 5'start?
  my @mrna = grep{ $_->[2] eq "mRNA" } @$generec; # must have
  my $mrna = $mrna[0];
  my @cds  = grep{ $_->[2] eq "CDS"  } @$generec;
  my @exons= grep{ $_->[2] eq "exon" } @$generec;
  my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid,$issplit)= @$mrna;

  # (my $gid2=$gid) =~ s/_[GC]\d+$//; #?
  #? $idlist{$gid}= $idlist{$gid2}= -9; # done flag
  # my $issplit= ($mrna->[8] =~ /splitgene|Split=/)?1:0;
  # my ($genegap)= ($generow) ? $generow->[8] =~ /(gap=[^;\s]+)/ : ""; # dont always have this annot.. opt
  # ^^ test this vs nfix

  # resort by loc, ignoring split
  @exons= sort _sortgene @exons;
  @cds=   sort _sortgene @cds;
  
  my $nexon= @exons; 
  my($lex,$lxr,$lxo,$iexon)= (0) x 9;
  my(@exok,@cdsok);
  foreach my $ex (@exons) {  
    $iexon++;
    # my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid,$issplit]; 
    my($xr,$xs,$xt,$xb,$xe,$xp,$xo,$tph,$xat)= @$ex;
    if($lex and $lxr eq $xr and  $lex->[4] > $xb and $lex->[3] < $xe) {
      # check for xe > lex.xe, extend lex.xe; record as errspan in tr > geno map
      # push @ierr, $iexon;
      $lex->[8] =~ s/$/;errspan=$iexon:$xb-$xe/;
      my $xcut = 1 + $xe-$xb;
      if($xe > $lex->[4]) { my $d= $xe - $lex->[4]; $lex->[4]=$xe; $xcut -= $d; }  
      $badspan += $xcut;
    } else {
      push @exok, $ex;
      $lex= $ex;
    }
    $lxr=$xr; $lxo=$xo;
  }

  unless($badspan) {
  
    # check if @cds is all in 1 split part, drop utr parts
    if($issplit and @mrna > 1) {
      my(%isp,$ispid); my $i=0;
      foreach my $ex (@cds) { 
        my($xr,$xs,$xt,$xb,$xe,$xp,$xo,$tph,$xat,$xid)= @$ex;
        my($isp,$xidss)= splitPart($xat);
        $isp{$isp}= ++$i; $ispid=$xid;
      }
      my @isp=sort keys %isp;
      if(@isp == 1) {
        $changed++;
        $issplit= $isp[0];
        my($mrnai)= grep{ $_->[8] =~ /_C$issplit;/ or $_->[8] =~ m/;Split=$issplit\b/ } @mrna;
        my @exoni= grep{ $_->[8] =~ /_C$issplit;/ or $_->[8] =~ m/;Split=$issplit\b/ } @exons;
        if($mrnai and @exoni) { @mrna=($mrnai); @exons= @exoni; }
      }
    }
    
    if($issplit and @mrna == 1) { # unsplit IDs .. have no 2nd path
      ## defer mrna loc update to mrna_fixspan()
      $mrna= $mrna[0];
      my($id)= $mrna->[8] =~ m/ID=([^;\s]+)/;  
      (my $nd=$id)=~s/_(C\d+)$//;
      $mrna->[8] =~ s/ID=$id/ID=$nd/; 
      $mrna->[8] =~ s/;chimera=[^;\n]+//; # drop
      unless($mrna->[8] =~ s/;Split=/;unsplit=/) { $mrna->[8] =~ s/$/;unsplit=$issplit/; }
      $changed++;
      for my $ex ( @exons, @cds) {
        my($id)= $ex->[8] =~ m/Parent=([^;\s]+)/;  
        (my $nd=$id)=~s/_(C\d+)$/;part=$1/;
        $ex->[8] =~ s/Parent=$id/Parent=$nd/; 
      }
      
    }
    if($changed) {
      my @newgenerec;
      if(@mrna > 1) {  @newgenerec = sort _sortSplitgene (@mrna, @exons, @cds); } 
      else { @newgenerec = sort _sortgene (@mrna, @exons, @cds); }
      $generec= \@newgenerec;
    }
    return ( $changed , $generec, $geneother); ## 0 == no change ??
    #was return ($ONLYCHANGES) ? 0 : putgene( \@generec, $geneother); ## $generow, $mrna,

  } else {
    ($lex,$lxr,$lxo,$iexon)= (0) x 9;
    foreach my $ex (@cds) {  
      $iexon++;
      my($xr,$xs,$xt,$xb,$xe,$xp,$xo,$tph,$xat)= @$ex;
      if($lex and $lxr eq $xr and  $lex->[4] > $xb and $lex->[3] < $xe) {
        $lex->[8] =~ s/$/;errspan=$iexon:$xb-$xe/;
        my $ccut = 1 + $xe-$xb;
        if($xe > $lex->[4]) { my $d= $xe - $lex->[4]; $lex->[4]=$xe; $ccut -= $d; }  
        $cdscut+= $ccut;
      } else {
        push @cdsok, $ex;
        $lex= $ex;
      }
      $lxr=$xr; $lxo=$xo;
    }
    
    #   my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid)= @$mrna;
    if(1) { # if($badspan or $mb != $tb or $me != $te)
      $changed++;
      ## fixme: other mrna attr changes: add cov=  of @mrna split parts, nexon=@exok, nintron= ???, path=2/2?
      ##  remove? chim[12]=.. and chimera= or rename?
      $mrna->[8] =~ s/;(chim[12]|chimera)=[^;\n]+//g; # drop
      # $issplit==0 dont add this unsplit? but have badspan/joinsplit fixes
      unless($mrna->[8] =~ s/;Split=/;unsplit=/) { $mrna->[8] =~ s/$/;unsplit=$issplit/; }
      my $nxok= @exok; $mrna->[8] =~ s/;nexon=\d+/;nexon=$nxok/; 
      my $covt=0; for my $m (@mrna) { if(my($c)=$m->[8] =~ m/;cov=(\d+)/) { $covt+=$c; } }
      $covt=100 if($covt>100); 
      $mrna->[8] =~ s/;cov=\d+/;nexon=$covt/ if($covt); 
      ## change from  joinsplit= to unsplit= ??
      $mrna->[8] =~ s/$/;joinsplit=xcut:$badspan,ccut=$cdscut/; ## ,oldspan:$tb-$te << DEFER_MRNA_FIXSPAN does 
      #  #err = "ERROR.span:genome_span:$gspan,tr_span:$tspan,$ref:$tb-$te";
     
      ## ?? defer mrna loc update to mrna_fixspan()
      unless(DEFER_MRNA_FIXSPAN) {
      my($mb, $me)=(0) x 9;
      foreach my $ex (@exok) {  
        my($xr,$xs,$xt,$xb,$xe,$xp,$xo,$tph,$xat)= @$ex;
        $mb= $xb if( $mb == 0 or $xb < $mb); 
        $me= $xe if( $xe > $me );
        }
      $mrna->[3]= $mb; $mrna->[4]= $me; 
      }
      
      my($id)= $mrna->[8] =~ m/ID=([^;\s]+)/;  
      (my $nd=$id)=~s/_(C\d+)$//;
      $mrna->[8] =~ s/ID=$id/ID=$nd/; 
      }
    for my $ex ( @exok, @cdsok) {
      my($id)= $ex->[8] =~ m/Parent=([^;\s]+)/;  
      (my $nd=$id)=~s/_(C\d+)$/;part=$1/;
      $ex->[8] =~ s/Parent=$id/Parent=$nd/; 
    }
  
    my @newgenerec= sort _sortgene ( $mrna, @exok, @cdsok); #?? @mrna ?
    return ( $changed , \@newgenerec, $geneother); ## 0 == no change ??
    #was return putgene( \@newgenerec, $geneother); ## $generow, $mrna,
  }
  
  return ( $changed , $generec, $geneother); ## 0 == no change ??
}

sub fixgene_cdstrans {
  my($geneid, $generec, $testrev)= @_;
  my($didchange)= (0);
  $testrev ||=0;
  
  my @mrna = grep{ $_->[2] eq "mRNA" } @$generec; # must have
  my $mrna = $mrna[0];
  my @cds  = grep{ $_->[2] eq "CDS"  } @$generec;
  # my @exons= grep{ $_->[2] eq "exon" } @$generec;
  my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid,$issplit)= @$mrna;
  my $gstrand= $to;

  ## FIXME: problem caused by indels in 5'UTR exons before CDS .. cdsbegin invalid
  ## can estimate offset, exon span vs target span, check with cdstrans()
  ## .. really need base alignment of mRNA x chr dna to pick out mistakes
  
  # ?? need to resort @cds by 5' to 3', Target start
  my($cdsb,$cdse,$cdsor0)= cds_span($geneid,$tattr);
  @cds= sort _sortgene @cds;
  
if(0) {  
  ## dang orientTargetCDS() NOT ready for @cdsgff, lacks Target..
  my($cdsor,$targor,$xrevorder,$chrstrand)= orientTargetCDS($cdsb,$cdse,$cdsor0,\@cds);
  if($targor and $cdsor ne $targor) { # if($chrstrand eq '-') is same
    # this is where cds.fwd (1..n) but target exons n..1, or other way, need cds ~ target order
    # chrstrand  is '-' when chror ne targor, + when eq
    @cds= reverse @cds; #? is this right? # need 5'cds first
  }
} else {
  if($gstrand eq '-') {
    @cds= reverse @cds unless($cdsor0 < 0);
  } else {
    @cds= reverse @cds if($cdsor0 < 0);
  }
}
  
  ## not working right ; many more InnerStop than cdsgff2genbank.pl for same cds.gff
  
  my @ctrans= getCDStrans($geneid, \@cds);
  ## @ctrans == ($aaflag,$phaseb,$aascore,$protaa,$cdsdna)
  # if($phaseb ne $phase0) .. change cds[0] phase
  my($ctransflag, $phaseb, $ctscore)= @ctrans[0,1,2];
  
  if($testrev and $ctscore < 2) { # only care when $ctrans.aaflag ne Complete ? 
    my $rstrand= ($gstrand eq '-')?'+':'-';
    my @revcds=map{ $_->[6]= $rstrand } reverse @cds;
    my @rtrans= getCDStrans($geneid, \@revcds);
    my($rtflag, $rphaseb, $rtscore)= @rtrans[0,1,2];
    if($rtscore > $ctscore) {
      $ctransflag.=",revaa:$rtflag"; # do more? change orient?
    }
  }
  
  $mrna->[8] =~ s/$/;cdstrans=$ctransflag/;
  return($didchange, $generec, $ctransflag, $phaseb);
}  

sub diffcds {
  my($generec,$cdsnew)= @_;
  ## $DIFF_CDS OPTION simple compare old to new CDS: same? strand change? exons differ?
  my($ndiff,$ordiff, $locdiff)=(0) x 9;
  my @oldcds = sort _sortgene grep{ $_->[2] eq "CDS"  } @$generec;
  my @newcds = sort _sortgene @$cdsnew;
  my $nnew= @newcds; my $nold= @oldcds;
  $ndiff = $nnew - $nold;
  if($nnew>0 and $nold>0) {
    if($nnew eq $nold) {
      
    } else {
    
    }
  }
  return($ndiff,$ordiff, $locdiff); 
}

sub fixgene_cdsupdate   
{
  my($geneid, $generec, $geneother)= @_;
  my($changed, $badspan,$cdscut)=(0) x 9; 

  # my @generec= sort _sortgene @$generec; #? sort by genostart or 5'start?
  my @mrna = grep{ $_->[2] eq "mRNA" } @$generec; # must have
  my $mrna = $mrna[0];
  # my @cds  = grep{ $_->[2] eq "CDS"  } @$generec;
  my @exons= grep{ $_->[2] eq "exon" } @$generec;
  my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid,$issplit)= @$mrna;
  my $gstrand= $to;
    
=item rev sort exons
  
  genefindcds.pl does it this way..
  @exongff = sort _sortgene grep(/exon/,@generec); # fwd chrloc sort
  my $revgene=($gstrand eq "-")?1:0; # NOTE exongff are rev order here; 5' is first
  if($issplit) { # sort by parts, by exon?
    @exongff= sort __revSplitGene @exongff;
  } elsif($revgene) {
    @exongff= reverse @exongff;  
  }

=item cdna vs mrna special case

  for un-oriented cdna mappings, check/correct for antisense map x cdsor-rev
  .. that is input cdna should have been flipped to be mRNA
  .. do before makeCDSexons?
  
=item add cds2prot tests, given chr.fasta

    -- new action from cdsupdate ? cds2prot
    -- before & after makeCDS? (i.e. check input CDS if there)
    -- check both fwd/rev for some (?all), esp. any antisense marked, or valintron != cdsor cases
    
=cut

  my($cdsb,$cdse,$cdsor)= cds_span($geneid,$tattr);
  if($cdse>0) {
    if($issplit) { # sort by parts, by exon?
      @exons= sort __revSplitGene @exons;
    } elsif($gstrand eq "-") {
      @exons= reverse @exons;  
    }

    my($cdsstrand, $cdsnew)= makeCDSexons($mrna,$cdsb,$cdse, $cdsor, $gstrand, \@exons);
    
    if($cdsnew and @$cdsnew) {
      
      #?? compare old cds to new, skip change if no diff?
      my @diffc= ($DIFF_CDS)? diffcds($generec,$cdsnew) : (0);
      
      $changed++;  
      my @newgenerec= sort _sortgene ( @mrna, @exons, @$cdsnew);
      if($cdsstrand ne $gstrand) {
        map{ $_->[6]= $cdsstrand; } @newgenerec;
        $mrna[0]->[8] =~ s/$/;flipstrand=$cdsstrand,$gstrand/; # ??
        $mrna[0]->[8] =~ s/;sense=-1/;oldsense=-1/; # ?? not all have such, how many? leave as is? changeto ;oldsense=?
        $changed++; #??
      }
      return ( $changed, \@newgenerec, $geneother); ## 0 == no change ??
    }
  }
  return ($changed , $generec, $geneother); ## 0 == no change ??
}

=item CDS orient problems
  
  ##  BAD cds-end clips wrong side for some: sense=-1 and strand=-, as one case
  ##  .. revor differing from exon strand is likely problem; 
  ## ** revor was wrong, need Target/mRNA orient here
  ##  .. may have same problem in splign2gff, others that make CDSexons
  ## CHANGE>> input xons are gloc sorted, xb0 < xb1, but target tb,te may be reversed
  ##  .. reverse @exons when chrstrand == -, so mRNA/Target start at 1st of exons
  ## may want chror, always clip 5' end for chr-lowest, 3'end for chr-highest, regardless of strand?
  ## ** another problem: target orient reverse of chr strand, ie putative antisense, but good prot..
  ## .. happens much w/ daphnia genes, need to calc "chror" by target align to chr, 
  ##    not by chr strand (from rna introns possibly .. this is where mismatch is from)
  ## ** also CDS-one-exon things can be either or.

=item case1: CDS/Target or reverse of chr-or  
  a. cdsoff=+ to transcript
  b. transcript aligns to +chr
  c. rna-introna are -chr
      ^^ b/c disagree, need to set cds-exons for b, flip strand from c.rna-introns

  orig gff, 
    c. rev -chr from rna-introna, but 
    b. tr Target orient is fwd +chr (exon1: T=1-409, exon2: T=409-469), and a. fwd offs=43-399:+ 
    
scaffold_7	velml10d	mRNA	1090761	1091310	100	-	.	ID=daplx6ml10dn9anvelvk43Loc9t2232;Target=daplx6ml10dn9anvelvk43Loc9t2232 1 469;qlen=470;cov=100;pid=100;match=469;
   nexon=2;intralign=2;aalen=118,75%,complete;offs=43-399:+;path=1/2;inhit=m1,1
scaffold_7	velml10d	exon	1090761	1091167	100	-	.	Parent=daplx6ml10dn9anvelvk43Loc9t2232;Target=daplx6ml10dn9anvelvk43Loc9t2232 1 409;align=409;introna=0,-2
scaffold_7	velml10d	CDS	1090771	1091125	100	-	0	Parent=daplx6ml10dn9anvelvk43Loc9t2232
scaffold_7	velml10d	exon	1091249	1091310	100	-	.	Parent=daplx6ml10dn9anvelvk43Loc9t2232;Target=daplx6ml10dn9anvelvk43Loc9t2232 409 469;align=61;introna=-2,0

  .. bad chr2aa from above gff:
  grep =daplx6ml10dn9anvelvk43Loc9t2232 velvdtest.infix.gff | .. cdsgff2genbank.pl
>daplx6ml10dn9anvelvk43Loc9t2232 loc=scaffold_7:1090771-1091125:-;type=CDS.velml10d;aaflag=Internalstop,;nx=1;len=118;partial_gene=true
LDPPYKILHTPNLIRFVCV*CISDI*LFISRES*CFKLRKRSSESNPSVPYVVIVEESFG
CR*VATIMSPSHQVTRVSCEDEIFAIVLETETLTESVFQIIACCSKFITELRCQFTFS
  
  .. proper chr2aa after flip strand, adjust cds start/end (1 cds-exon)
  grep =daplx6ml10dn9anvelvk43Loc9t2232 velvdtest.infix.gff | head -5 | sed 's/     \-/     +/;s/1090771/1090803/; s/1091125/1091158/;' | \
  $evigene/scripts/cdsgff2genbank.pl -debug  -a aa -t CDS -pretty -gff stdin -fasta ../genome/gasm16ml/daphplx_gasm16ml.fa 
>daplx6ml10dn9anvelvk43Loc9t2232 loc=scaffold_7:1090803-1091158:+;type=CDS.velml10d;aaflag=Start,;nx=1;len=118
MNLEQQAMIWNTLSVSVSVSRTMAKISSSQDTRVTWWLGDMIVATYRQPKDSSTITTYGT
DGLLSELRFRSLKHHDSREMNNQISDIHHTQTNLMRLGVWSILYGGSRQKWTASFSMA

=item case 6s: gsplign, Antisense=-1, +cdsoff (mRNA), +chr-tr align, -chr-splices
  .. gff-trans is bad, pid=99 could be factor
  .. change to +chr or, cds b/e adjust
  
grep =Daplx7b3EVm000781t7  gsplign7merr.gff    
scaffold_140	splign	mRNA	253685	260234	100	-	.	ID=Daplx7b3EVm000781t7;cov=100%,1772/1780;pid=99.5;nexon=5;splice=12;
  sense=-1;Target=Daplx7b3EVm000781t7 1 1780;gescore=84;aalen=375,63%,complete;clen=1780;offs=45-1172;oid=Daplx6mlEVm000472t58,tidbdaplx6ml25mrno2lridk33Loc51739
scaffold_140	splign	exon	253685	254129	0.998	-	.	Parent=Daplx7b3EVm000781t7;Target=Daplx7b3EVm000781t7 1 445;splice=nnCC+
scaffold_140	splign	exon	254708	255034	0.997	-	.	Parent=Daplx7b3EVm000781t7;Target=Daplx7b3EVm000781t7 446 772;splice=TACT+
scaffold_140	splign	exon	255199	255534	0.991	-	.	Parent=Daplx7b3EVm000781t7;Target=Daplx7b3EVm000781t7 773 1108;splice=ACCT+
scaffold_140	splign	exon	255904	256139	1	-	.	Parent=Daplx7b3EVm000781t7;Target=Daplx7b3EVm000781t7 1109 1344;splice=ACCT+
scaffold_140	splign	exon	259797	260234	0.989	-	.	Parent=Daplx7b3EVm000781t7;Target=Daplx7b3EVm000781t7 1345 1780;splice=ACnn+
scaffold_140	splign	CDS	253685	254085	1	-	.	Parent=Daplx7b3EVm000781t7 << bad CDS-e adjust, should be CDS-b
        >> should be  CDS 253685+44  254129 ; = 253729
scaffold_140	splign	CDS	254708	255034	1	-	.	Parent=Daplx7b3EVm000781t7
scaffold_140	splign	CDS	255199	255534	1	-	.	Parent=Daplx7b3EVm000781t7
scaffold_140	splign	CDS	256076	256139	1	-	.	Parent=Daplx7b3EVm000781t7  << bad CDS-b adjust, should be CDS-e
        >> should be  CDS 255904  256139-172 (1344 - 1172) = 255967 or 255968 (255904 + 64, from Trg 1108 + 64)
        >> try 255904 + 44 (cds-b off) = 255948, nogo

........ splign of Daplx7b3EVm000781t7
gsplignt5nmerr/maperr2set_part1.splog
+10153	Daplx7b3EVm000781t7	scaffold_140	Ok	-78.754
zgrep Daplx7b3EVm000781t7  gsplignt5nmerr/maperr2set_part1.splign.gz 
+10153	Daplx7b3EVm000781t7	scaffold_140	0.998	445	1	445	253685	254129	  <exon>CC	M83RM361
+10153	Daplx7b3EVm000781t7	scaffold_140	0.997	327	446	772	254708	255034	TA<exon>CT	M190RM136
+10153	Daplx7b3EVm000781t7	scaffold_140	0.991	336	773	1108	255199	255534	AC<exon>CT	M175RM10R2M148
+10153	Daplx7b3EVm000781t7	scaffold_140	1	236	1109	1344	255904	256139	AC<exon>CT	M236
+10153	Daplx7b3EVm000781t7	scaffold_140	0.989	438	1345	1780	259797	260234	AC<exon>  	M78RM24IM60RM219RM43IM9
.......
        
  gff-trans
>Daplx7b3EVm000781t7 loc=scaffold_140:253685-256139:-;type=CDS.splign;aaflag=Internalstop,;nx=4;len=376;partial_gene=true
KRETIWD*NFNKMN*AKNVVDVCIGTNGRMSVPSNREHHYRNLRDRYTNCTYVDGNLELT
WLQDEHLDLTFLQHIREVTGYVLISHVDVRRIVLPALQIIRGRTLFKLNVHDEEFALVVT
LSKMHNLEMPALRDILAGSCGLFNNYNLCHMKTINWEEIITGSNSKIFNVYNFTEPEREC
PACHESCEAGCWGEGPENCQKFSKINCSPQCHQGRCFGPKPRECCHLFCAGGCTGPKQSD
CLAGHAALQPHHVLLGNQPRGQIRLRRHLRQNVSGAFAQGQRGLRPHLPAQEEGRQRRMR
RLRRPVPKDLPGRRYRPRRKHRELPRLHRHRGIHHDPGSFLRRIPADLSQFHLRPALPAT
PSGPAGNIRHPPGDYW

 gff-trans, fix CDS1-b/e, CDS4-e/b, not complete..
>Daplx7b3EVm000781t7 loc=scaffold_140:253729-255967:-;type=CDS.splign;nx=4;len=376;partial_gene=true
SNVTSHHSSPADKANEVFKGKICIGTNGRMSVPSNREHHYRNLRDRYTNCTYVDGNLELT
WLQDEHLDLTFLQHIREVTGYVLISHVDVRRIVLPALQIIRGRTLFKLNVHDEEFALVVT
LSKMHNLEMPALRDILAGSCGLFNNYNLCHMKTINWEEIITGSNSKIFNVYNFTEPEREC
PACHESCEAGCWGEGPENCQKFSKINCSPQCHQGRCFGPKPRECCHLFCAGGCTGPKQSD
CLACRNFFDDGVCKQECPAMQRYNPITYSWEINPEGKYAYGATCVKTCPEHLLKDNGACV
RICPPKKKAVNGECVVCDGPCPKTCQGVDIVHAGNIESFRDCTVIEGSITILDHSFAGFQ
QIYPNFTFGPRYPRLH

 gff-trans, +chr or, fix CDS1-b/e, CDS4-e/b, this is it
>Daplx7b3EVm000781t7 loc=scaffold_140:253729-255968:+;type=CDS.splign;aaflag=Complete;nx=4;len=376
MESRVARAEGEIGINLLESGEGMIQDRDGSLDDGAIAEALDVSGVDDIDALAGLWARAVA
DDAFAVDGLLLGRADADAGPVVLEQMLRTRFDAGGAVGVFALGVDFPRVRDGVVALHGRA
FLLAHAVVEEIPASETVALLGSGASAGAKQVAALARLGSEATPLVALRTAVDFGEFLAIL
GPLAPAAGLARLVTSRTFTLRLREIVNVENLGIGARDDLLPVDRLHVAQVVVVEQSATTR
QYISEGGHFEVVHFGEGDDEGEFLVVDVEFEQSPSADDLQGRQDDAADVHVRDEDVAGDF
ADVLQKSQVQVLVLEPSQFQVAVDVRAVGVAVAQIPVVVLPVGRHRHPPVRPDANFAFKD
FVGFIGGTRVMGRDV*
  
=item case 5: 2exon, +cdsoff, -chr-tr, +chr-intron
  -- expect bad trans
  -- change to -chr, cds b/e
   ** BUG from blast2evgff.pl, Target offsets are off by introna bits, in these velvdtest.infix.gff data
   
 grep =daplx6ml10dn9anvelvk43Loc7651t1_G2 velvdtest.infix.gff | head
scaffold_88	velml10d	mRNA	400053	402707	100	+	.	ID=daplx6ml10dn9anvelvk43Loc7651t1_G2;Target=daplx6ml10dn9anvelvk43Loc7651t1 10 1805;qlen=1844;cov=97;pid=100;match=1796;nexon=3;intralign=3;aalen=464,75%,complete;offs=93-1487:+;path=2/3;inhit=m2,2
scaffold_88	velml10d	exon	400053	400204	100	+	.	Parent=daplx6ml10dn9anvelvk43Loc7651t1_G2;Target=daplx6ml10dn9anvelvk43Loc7651t1 1649 1805;align=157;introna=0,-5
      ^^ Target=id 1654 1805 from introna -5
scaffold_88	velml10d	exon	400312	400484	100	+	.	Parent=daplx6ml10dn9anvelvk43Loc7651t1_G2;Target=daplx6ml10dn9anvelvk43Loc7651t1 1481 1656;align=176;introna=-3,0
      ^^ Target=id 1481 1653 from introna -3
scaffold_88	velml10d	CDS	400481	400484	100	+	0	Parent=daplx6ml10dn9anvelvk43Loc7651t1_G2
scaffold_88	velml10d	exon	401237	402707	100	+	.	Parent=daplx6ml10dn9anvelvk43Loc7651t1_G2;Target=daplx6ml10dn9anvelvk43Loc7651t1 10 1481;align=1472;introna=-1,0
      ^^ Target=id 10 1480 from introna -1
scaffold_88	velml10d	CDS	401320	402707	100	+	2	Parent=daplx6ml10dn9anvelvk43Loc7651t1_G2

   
=item case 4: 2exon, +cdsoff, +chr-tr, -chr-intron
  -- expect bad trans: yes
  -- change to +chr, but cds b/e insets are +2 for both ends, no change: yes
grep =daplx6ml10dn9anvelvk91Loc12111t1 velvdtest.infix.gff | head
scaffold_6	velml10d	mRNA	552529	552959	100	-	.	ID=daplx6ml10dn9anvelvk91Loc12111t1;Target=daplx6ml10dn9anvelvk91Loc12111t1 1 382;qlen=382;cov=100;pid=100;match=382;
  nexon=2;intralign=2;aalen=126,98%,partial;offs=3-380:+;inhit=m1,1
scaffold_6	velml10d	exon	552529	552815	100	-	.	Parent=daplx6ml10dn9anvelvk91Loc12111t1;Target=daplx6ml10dn9anvelvk91Loc12111t1 1 293;align=293;introna=-5,-1
scaffold_6	velml10d	CDS	552531	552815	100	-	2	Parent=daplx6ml10dn9anvelvk91Loc12111t1
scaffold_6	velml10d	CDS	552870	552957	100	-	0	Parent=daplx6ml10dn9anvelvk91Loc12111t1
scaffold_6	velml10d	exon	552870	552959	100	-	.	Parent=daplx6ml10dn9anvelvk91Loc12111t1;Target=daplx6ml10dn9anvelvk91Loc12111t1 293 382;align=90;introna=0,0

 -orig trans
>daplx6ml10dn9anvelvk91Loc12111t1 loc=scaffold_6:552531-552957:-;type=CDS.velml10d;aaflag=Internalstop,;nx=2;len=124;partial_gene=true
ETSSNITNNRLYIGTTSDKVGNPLPSIYLSSLYISFLCRSASLAISSMAGLSFTARLRK*
*VELHCDV*SKIGRRSRAITEPSGNAISSSSS*YSSSCNPHDSRTFLNF*ILRCEARIEL
MRSR
 
 +chr-gff trans
>daplx6ml10dn9anvelvk91Loc12111t1 loc=scaffold_6:552531-552957:+;type=CDS.velml10d;nx=2;len=124;partial_gene=true
RDLINSIRASHRRIQKFKNVLESCGLHDEEYQEEEEEMAFPEGSVIARDLLPILDHTSQW
SSTYYFLRRAVKLRPAIDEIAKEADLQRNEMYKEDKYIEGSGFPTLSLVVPMYNRLLVIL
EDVS
   
=item case 3: multi-exon, >> -cdsoff, +chr-tr align, -chr-intron 
  -cdsoff differs from case1 +off, match -chrintrons, so gff should translate ok
  
grep =daplx6ml10dn9anvelvk63Loc5218t1 velvdtest.infix.gff | head
scaffold_73	velml10d	mRNA	106825	108819	100	-	.	ID=daplx6ml10dn9anvelvk63Loc5218t1;Target=daplx6ml10dn9anvelvk63Loc5218t1 1 1543;qlen=1549;cov=100;pid=100;match=1543;
  nexon=5;intralign=5;aalen=170,33%,complete-utrbad;offs=1163-651:-;path=1/4;inhit=m3,3
scaffold_73	velml10d	exon	106825	107165	100	-	.	Parent=daplx6ml10dn9anvelvk63Loc5218t1;Target=daplx6ml10dn9anvelvk63Loc5218t1 1 343;align=343;introna=0,-2
scaffold_73	velml10d	exon	107318	107657	100	-	.	Parent=daplx6ml10dn9anvelvk63Loc5218t1;Target=daplx6ml10dn9anvelvk63Loc5218t1 345 686;align=342;introna=-1,-1
scaffold_73	velml10d	CDS	107624	107657	100	-	2	Parent=daplx6ml10dn9anvelvk63Loc5218t1
scaffold_73	velml10d	exon	107719	107784	100	-	.	Parent=daplx6ml10dn9anvelvk63Loc5218t1;Target=daplx6ml10dn9anvelvk63Loc5218t1 686 753;align=68;introna=-1,-2
scaffold_73	velml10d	CDS	107719	107784	100	-	1	Parent=daplx6ml10dn9anvelvk63Loc5218t1
scaffold_73	velml10d	exon	107861	107929	100	-	.	Parent=daplx6ml10dn9anvelvk63Loc5218t1;Target=daplx6ml10dn9anvelvk63Loc5218t1 753 822;align=70;introna=-2,-2
scaffold_73	velml10d	CDS	107861	107929	100	-	2	Parent=daplx6ml10dn9anvelvk63Loc5218t1
scaffold_73	velml10d	CDS	108220	108439	100	-	0	Parent=daplx6ml10dn9anvelvk63Loc5218t1
scaffold_73	velml10d	exon	108220	108819	100	-	.	Parent=daplx6ml10dn9anvelvk63Loc5218t1;Target=daplx6ml10dn9anvelvk63Loc5218t1 944 1543;align=600;introna=0,0

cdsgff2genbank.pl 
  -- not quite full aa, missing 40aa an 3'stop,
  -- problem due to missed coverage 120 nt (40aa) in exon2 (rev) == MGap of splign
>daplx6ml10dn9anvelvk63Loc5218t1 loc=scaffold_73:107624-108439:-;type=CDS.velml10d;aaflag=Start,;nx=4;len=129
MSSPYLYMEAFQITSSIKLEDIKVYQSFSCFSEEFEELLGHSDEIEHLSENRIVKIKCSS
QDPLKGQPSFHFEACFLQPTTQFVTSNKVSVVTVMEGLADDSTDLTVLNVNEETILNMGT
MIVTVEENI

=item case 2: CDS/Target or matches chr or

  a. cdsoff=+ to transcript
  b. transcript aligns to -chr
  c. rna-introna are -chr
  
grep =daplx6ml10dn9anvelvk59Loc23312t1 velvdtest.infix.gff | head
scaffold_81	velml10d	mRNA	560186	560836	100	-	.	ID=daplx6ml10dn9anvelvk59Loc23312t1;Target=daplx6ml10dn9anvelvk59Loc23312t1 1 584;qlen=584;cov=100;pid=100;match=584;
      nexon=3;intralign=3;aalen=171,87%,partial3;offs=70-582:+;inhit=m2,2
scaffold_81	velml10d	exon	560186	560491	100	-	.	Parent=daplx6ml10dn9anvelvk59Loc23312t1;Target=daplx6ml10dn9anvelvk59Loc23312t1 277 584;align=308;introna=-1,-1
scaffold_81	velml10d	CDS	560188	560491	100	-	2	Parent=daplx6ml10dn9anvelvk59Loc23312t1
scaffold_81	velml10d	exon	560565	560680	100	-	.	Parent=daplx6ml10dn9anvelvk59Loc23312t1;Target=daplx6ml10dn9anvelvk59Loc23312t1 94 211;align=118;introna=0,-2
scaffold_81	velml10d	CDS	560565	560680	100	-	0	Parent=daplx6ml10dn9anvelvk59Loc23312t1
scaffold_81	velml10d	CDS	560742	560767	100	-	0	Parent=daplx6ml10dn9anvelvk59Loc23312t1
scaffold_81	velml10d	exon	560742	560836	100	-	.	Parent=daplx6ml10dn9anvelvk59Loc23312t1;Target=daplx6ml10dn9anvelvk59Loc23312t1 1 96;align=96;introna=-1,0

=cut
  
## later: add chrasm.fasta seq check: does chr-cds-seq convert to valid prot, same size?

sub phaseCDS {
  my($tb,$te,$cb,$ce,$cdsor, $cdslen)=@_;
  # tb,te= target_begin, target_end; cb,ce,cdsor = cds_begin,_end,_orient(-1,1), running cdslen
  
  # phase calc; FIXME partial5 may need start phase
  my($xend5,$xend3)= ($cdsor < 0)? ($te,$tb) : ($tb,$te); # if input is mRNA, tb == xend5
  my $d5= ($tb >= $cb) ? 0 : $cb - $tb; # pos
  my $c5= ($xend5 > $xend3) ? $xend5 - $d5 : $xend5 + $d5;    				  
  my $d3= ($te <= $ce) ? 0 : $ce - $te; # neg
  my $c3= ($xend5 > $xend3) ? $xend3 - $d3 : $xend3 + $d3;

  my $xwide  = 1 + abs($c3 - $c5);
  $cdslen   += $xwide;
  my $inc3   = $cdslen % 3;
  my $inc5   = ($xwide - $inc3) % 3; # only care about this one
  if ($inc5 == -1) { $inc5 = 2; }
  return( $inc5, $cdslen); # is this right?
}


sub orientTargetCDS {
  my($cb,$ce,$cdsor,$xons)= @_; 
	my($targor, $xrevorder)= (0) x 9; 
	return($cdsor,0,0,0) unless(ref($xons) and $xons->[1] and $xons->[1]->[3]);
  if($ce > 0 and $ce < $cb) { 
    ($cb,$ce)= ($ce,$cb); $cdsor= -$cdsor; #?
  }
  $xrevorder= ($xons->[0]->[3] > $xons->[1]->[3]) ? -1 : 0; # exon sort order, fwd or rev?
  my($tb1)= $xons->[0]->[8] =~ m/(?:Target|trg)=\S+\s(\d+)/; # Target MISSING for CDS exons
  my($tb2)= $xons->[1]->[8] =~ m/(?:Target|trg)=\S+\s(\d+)/;
	if($xrevorder) { $targor=($tb2 and $tb1 > $tb2)? 1 : ($tb1 and $tb1 < $tb2) ? -1 : 0; }
  else { $targor=($tb2 and $tb1 > $tb2)? -1 : ($tb1 and $tb1 < $tb2) ? 1 : 0; }
	
	my $chrstrand= 
	  ($targor > 0 and $cdsor >= 0) ? '+' :
	  ($targor > 0 and $cdsor  < 0) ? '-' :  
	  ($targor < 0 and $cdsor  < 0) ? '+' :  # CDSb == chr-start
	  ($targor < 0 and $cdsor >= 0) ? '-' :  # CDSb == chr-end
	  ''; # ($chror < 0)?'-':'+' ; # cant tell, use input
	  	 
  return($cdsor,$targor,$xrevorder,$chrstrand); #??  
}


sub makeCDSexons {
  my($mrna,$cb,$ce,$cdsor,$gstrand,$xons,$AS_STRING)= @_; 
  
  ## fixme for $cb > $cd .. flip or
	return() unless($cb and $ce and ref($xons)); #$ce > $cb
  ## UPD17: this should check, return new chr-strand matching CDS-on-genome strand
  ## .. cdsor is needed, but input chror should be ignored. Calc proper chr-or from
  ## cdsor and chr-trans align orient (exon Target order vs exon chr order). See above notes
  
  ## FIXME5: 16.10.11 : sense=-1; flipor CDS can be bad, wrong calc from exons, offset

	my @cdsgff=(); 
	my ($cdslen,$phase, $cdsfirst, $targor, $xrevorder, $chrstrand)=(0) x 9; 
  $cdsfirst= 1;  
  if($ce > 0 and $ce < $cb) { 
    ($cb,$ce)= ($ce,$cb); $cdsor= -$cdsor; #?
  }
  
  if(1) {
    ($cdsor,$targor,$xrevorder,$chrstrand)= orientTargetCDS($cb,$ce,$cdsor,$xons);
    $chrstrand ||= $gstrand;
  } else {
    $xrevorder= ($xons->[0]->[3] > $xons->[1]->[3]) ? -1 : 0; # exon sort order, fwd or rev?
    my($tb1)= $xons->[0]->[8] =~ m/(?:Target|trg)=\S+\s(\d+)/;
    my($tb2)= $xons->[1]->[8] =~ m/(?:Target|trg)=\S+\s(\d+)/;
    if($xrevorder) { $targor=($tb2 and $tb1 > $tb2)? 1 : ($tb1 and $tb1 < $tb2) ? -1 : 0; }
    else { $targor=($tb2 and $tb1 > $tb2)? -1 : ($tb1 and $tb1 < $tb2) ? 1 : 0; }
    $chrstrand= 
      ($targor > 0 and $cdsor >= 0) ? '+' :
      ($targor > 0 and $cdsor  < 0) ? '-' :  
      ($targor < 0 and $cdsor  < 0) ? '+' :  # CDSb == chr-start
      ($targor < 0 and $cdsor >= 0) ? '-' :  # CDSb == chr-end
      $gstrand; # ($chror < 0)?'-':'+' ; # cant tell, use input
	}
	  	 
  foreach my $x (@$xons) {
    my @c= @$x; # only this way here
    # my @c; # 
    # if(ref($x) =~ /ARRAY/) { @c= @$x; } # ?? ASARRAY, return same as xons form?
    # elsif($x =~ /\texon\t/) { @c= split"\t",$x; chomp($c[-1]); }	    
    # else { next; } # error
    
    # my $xor= ($revor < 0) ? '-' : '+'; # or $c[6]; # +/- strand
    my($chrb,$chre,$chrorx)= @c[3,4,6];
    my($tid,$tb,$te)= $c[8] =~ m/(?:Target|trg)=(\S+)\s(\d+)\s(\d+)/;
    next unless($te); # warn/error
    
    if($te > $cb and $tb < $ce) {
      ## fixme again; -strand needs c3,c4 swap also
      # FIXME: still bad, targor < 0, cdsor > 0, +gstrand , from gmap antisense; not xrevorder
      # CDSe hits this:  if($te >= $ce) { .. NOT ($xrevorder and $cdsfirst) > $c[3] += $d; should be c[4] += d

      if($tb <= $cb) { my $d=$cb-$tb;  ## chror < 0 implies chre == tb, col[4]
        if($xrevorder) { if($cdsfirst) { $c[4] -= $d; } else { $c[3] += $d; } }
        else { if($cdsfirst) { $c[3] += $d; } else { $c[4] -= $d; } }
        $cdsfirst= 0;
        #y if($chror<0) { $c[4] -= $d; } else { $c[3] += $d; }  
        #x if($cdsor<0) { $c[4] -= $d; } else { $c[3] += $d; } # FIXME: - bad or offset; flip d sign, $or * $d
        }
      if($te >= $ce) { my $d=$te-$ce; ## chror < 0 implies chrb == te, col[3]
        if($xrevorder) { if($cdsfirst) { $c[4] -= $d; } else { $c[3] += $d; } }
        else { if($cdsfirst) { $c[3] += $d; } else { $c[4] -= $d; } }
        $cdsfirst= 0;
        #y if($chror<0) { $c[3] += $d; } else { $c[4] -= $d; }
        #x if($cdsor<0) { $c[3] += $d; } else { $c[4] -= $d; }
        }

       
      ($phase,$cdslen)= phaseCDS($tb,$te,$cb,$ce,$cdsor,$cdslen);
      # if(0) {
      #   # phase calc; FIXME partial5 may need start phase
      #   my($xend5,$xend3)= ($cdsor < 0)? ($te,$tb) : ($tb,$te); # if input is mRNA, tb == xend5
      #   my $d5= ($tb >= $cb) ? 0 : $cb - $tb; # pos
      #   my $c5= ($xend5 > $xend3) ? $xend5 - $d5 : $xend5 + $d5;    				  
      #   my $d3= ($te <= $ce) ? 0 : $ce - $te; # neg
      #   my $c3= ($xend5 > $xend3) ? $xend3 - $d3 : $xend3 + $d3;
      # 
      #   my $xwide  = 1 + abs($c3 - $c5);
      #   $cdslen   += $xwide;
      #   my $inc3   = $cdslen % 3;
      #   my $inc5   = ($xwide - $inc3) % 3; # only care about this one
      #   if ($inc5 == -1) { $inc5 = 2; }
      #   $phase= $inc5; # is this right?
      # }
      
      $c[2]= "CDS";
      # $c[5]=1; # score.. leave as exon score (alignment?)
      # $c[6]= $chrstrand; # do this here? should change exons also
      $c[7]= $phase;  
      $c[8] =~ s/;Target=.*//;
      
      my $xout= ($AS_STRING) ? join("\t",@c)."\n" : \@c;  # which? array or string?
      push @cdsgff, $xout;
    }
  }
 
	return ($chrstrand, \@cdsgff);
}


sub readIntronOk {
  my($inf)= @_;
  my(%inokerr, %introns); # == inokerr, introns
  my ($nin,$nerr,$nval)=(0,0,0);
  my $ok= open(F,$inf);
  while(<F>) {
    next if(/^\W|^nogene/); 
    my @v=split;
    my($td,$ind,$inw,$iv,$xtype)=@v; # intronok.tab
    $td=~s/Dpx6imEVm/Daplx7b3inewEVm/;  # ugh, idprefix mixup: idtab Daplx7b3inewEVm03716t1 == inoktab Dpx6imEVm03716t1
    $xtype ||= "other"; # bug?
    
    # my($tclass)= classof($td);
    my($vi)= m/valintron=(\d+)/; 
    my($as)= m/antisense=(\d+)/; 
    my($obe)= m/offby=(\d+)/; # upd1707, gsplign +1/-1 splice site when antisense
    my($ler)= ($vi < $LONGOK_MAXin and $inw > $LONGOK_MAX)?1:0;
    
    my $iflag="";
    $iflag.="${xtype}ok," if($vi>0 and not $ler); # exclusive of errs  
    # $iflag.="antierr," if($as>0); #** Make this OPTION : error or not
    if($as>0) { $iflag.= (($ANTISENSE_INTRON_IS_ERROR) ? "antierr," : "antisense,");  }
    if($obe>0) { $iflag.= (($ANTISENSE_INTRON_IS_ERROR) ? "offbyerr," : "offby,");  }
    $iflag.="longerr," if($ler>0);
    $iflag ||= "noevd";
    $inokerr{$td} .= "$iflag,";
    $introns{$td}{$ind}= join"\t",$iflag,$inw,$iv,$xtype;
    $nin++; $nerr++ if($iflag=~/err/); $nval++ if($iflag=~/ok/);
    # $tclass{$td}= $tclass;
    # $tderr{$td}++ if($iflag=~/err/); # unless($iflag=~/inok/);
    # $tdinok{$td}++ if($iflag=~/(CDS|UTR)ok/);
  } close(F);
  warn  "# readIntronOk n=$nin,ok=$nval,errs=$nerr from $inf\n" if $debug;
  return($nin, \%inokerr, \%introns);
}

sub readSizes {
  my(@inf)= @_; 
  my $nt=0;
  my (%alen,%trlen,%cdspan,%aaqual); %alen=(); 
  my $hasgap=0; # ($AAGAP)?1:0;#  && $iscount .. my $iscount=1; 
  
  # my $CDSSPAN=1; #  want these.. but check table
  my $hasspan=0; # test for it.. $CDSSPAN; #  collect %trlen,%cdspan ? ** TEST sizes input for this?
  my $testspan=1; #(defined $CDSSPAN)?0:1; #? always test, ignore $CDSSPAN ?
  
  foreach my $aaf (grep/\w/, map{ split(/[,\s]+/,$_) }@inf) {
    if(open(F,$aaf)) { my $n=0;
      while(<F>){ 
        next if(/^\W/); chomp; my($id,$aw,@ac)=split"\t"; 
        if($hasgap) { if(@ac>1 or @ac==0){ $hasgap=0; } else { $aw -= $ac[0]; }} 
        ## dang new cds.qual has Code/Noncode col before cdspan col ..
        if($hasspan or $testspan) { 
          if(@ac>=3 and $ac[2]=~/^\d/) { 
	          my($gp,$aq,$tw,$csp,$cspx)=@ac; # csp == span OR Code/Noncode col ..
	          my $isutrorf=(m/utrorf/)?1:0; # key in $oid may be missing
	          my $cspan= ($csp=~/^\d/) ? $csp : ($cspx=~/^\d/) ? $cspx : 0;
	          if($testspan) {
	            #x if(($csp=~/^\d/ or $cspx=~/^\d/)) { $hasspan=1; $testspan=0; }
	            if($cspan) { $hasspan=1; $testspan=0; }
	            else { if(++$testspan>9) { $hasspan=0; $testspan=0; } }
	          }
            $aaqual{$id}=$aq; $trlen{$id}=$tw;
	          #?? still buggy? $csp=$cspx=""  if($isutrorf); # data bug: bad cds-offset for utrorfs, fixme
	          ## offs=1340-90:- << reverse strand parts
            $cdspan{$id}= $cspan; # ($csp=~/^\d/) ? $csp : ($cspx=~/^\d/) ? $cspx : 0; # Code-col=3?  
            } 
          else { if(++$testspan>9) { $hasspan=0; $testspan=0; } } 
        } 
        $alen{$id}=$aw; $n++; 
      } close(F); 
      
      $nt+=$n; warn  "# readSizes n=$n from $aaf\n" if $debug;
    } else {
      warn "# cant read sizes from $aaf\n" ;# if $debug
    }
  }
  return($nt,\%alen,\%trlen,\%cdspan,\%aaqual); # change to 1 hash w/ fields?
}

## translation methods, from evigene/scripts/cdsgff2genbank.pl
## some of these methods belong in evgxxx.pm package

sub getCDStrans {
  my($gid,$cdsgff)= @_;
  # check/sort @cds from 5' to 3', need strand?
  my $trustphase= 0; #option
  my ($cdsdna,$phase0,$complete)=("",0,0); # dorev: reverse when -strand cds
  my $nx= @$cdsgff;
  use constant DOREVDNA=>1;
  
  for(my $ix=0; $ix<$nx; $ix++){
    my $exonft= $cdsgff->[$ix];
    # my($ref,$start,$stop,$strand,$phase,$xat,$xid)= @{$exonft}[0,3,4,6,7,8,9];
    # my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid]; 
    my ($xdna,$isrev,$xlen)= getexondna($exonft, DOREVDNA);
    $cdsdna .= $xdna;
    if($ix == 0) {
      $phase0= $exonft->[7];  
      my $ATG= uc(substr($xdna, $phase0, 3));
      $complete |= 1 if($ATG eq 'ATG');    
    }
  }
  
  my $phaset= ($complete & 1 || $trustphase)? $phase0 : undef;
  my($phaseb,$protaa,$aascore,$aaflag) = getBestFrame( $cdsdna, $gid, $complete & 1, $phaset); 
  ## aascore = 2/Complete, 1/StartOrStop, <0 Innerstops
  #* add flag for $phaseb ne $phase0 : $aaflag .=",Newphase:$phaseb" ?
  return($aaflag,$phaseb,$aascore,$protaa,$cdsdna);
}

sub getexondna {
  my($exonft,$dorevcomp)= @_; 
  my($ref,$start,$stop,$strand,$phase,$xgid)= @{$exonft}[0,3,4,6,7,9];
  my $isrev=($strand eq '-')?1:0;
  my $exondna  = get_dna( $chrfasta, $ref, $start, $stop);
  $exondna ||="";
  if($debug and not $exondna) {
    warn "#getexondna=0 gid=$xgid loc=$ref:$start-$stop\n"; # loc=NOPATH:1-69
  }
  $exondna = uc($exondna); # always?
  my $xlen = length($exondna); # 1+$stop-$start; # add even if no dna
  if($dorevcomp and $isrev) {  
    $exondna = reverse $exondna;
    $exondna =~ tr/gatcGATC/ctagCTAG/;
    }
  return ($exondna,$isrev,$xlen);
}

sub get_dna {
  my($fasta, $ref, $start, $stop)= @_; #, $fasta_db_ref
  unless( $ref && $stop>0) { warn "need ref-ID, start, stop location\n"; return; }
 
  require Bio::DB::Fasta;  # FIXME: not portable w/o parts of BioPerl !
  my $havedb= ref $fasta_db;
  unless($havedb) {
    my $db = eval { Bio::DB::Fasta->new($fasta); } or die "$@\nCan't open sequence file(s). "; # return ?
    $fasta_db= $db;  
    }
  
  my $seq = $fasta_db->seq($ref, $start => $stop)  or return;  ## die "cant locate seq of $ref:$start-$stop\n"; 
  $seq= $seq->seq if(ref $seq); #  weird bioperl change  
  return $seq;
}


my @s5CodonTable = ();
BEGIN {
 @s5CodonTable = (
 [ ['K','N','K','N','X',], ['T','T','T','T','T',], ['R','S','R','S','X',], ['I','I','M','I','X',], ['X','X','X','X','X',],],
 [ ['Q','H','Q','H','X',], ['P','P','P','P','P',], ['R','R','R','R','R',], ['L','L','L','L','L',], ['X','X','X','X','X',],],
 [ ['E','D','E','D','X',], ['A','A','A','A','A',], ['G','G','G','G','G',], ['V','V','V','V','V',], ['X','X','X','X','X',],],
 [ ['*','Y','*','Y','X',], ['S','S','S','S','S',], ['*','C','W','C','X',], ['L','F','L','F','X',], ['X','X','X','X','X',],],
 [ ['X','X','X','X','X',], ['X','X','X','X','X',], ['X','X','X','X','X',], ['X','X','X','X','X',], ['X','X','X','X','X',],],
 );
}

sub ibase {
  my $c= substr($_[0],$_[1],1);
  return ($c eq 'A')? 0: ($c eq 'C')? 1: ($c eq 'G')? 2: ($c eq 'T')? 3: 4;
}  
  
sub translate {
  my($cds, $offset)= @_; # assumes uc($cds);
  my $aa=""; my $aa_length = int((length($cds) - $offset) / 3);
	for (my $i = 0; $i < $aa_length; $i++) {
		my $idx = 3 * $i + $offset;
		$aa .= $s5CodonTable[ ibase($cds,$idx)][ ibase($cds,$idx+1) ][ ibase($cds,$idx+2) ];
	}
  return $aa; 
}

# aaflag == 'Complete|Start,Stop,Internalstop'
sub getBestFrame {
  my($cds, $id, $isfullcds, $usephase)= @_;
  my ($bestscore, $besti, $bestpro, $bestflag)= (-999999,0,"","");
  my $ph0= ($usephase =~ /\d/) ? $usephase : 0;
  my $nframe= ($isfullcds ? $ph0+1 : 3);
  $cds = uc($cds);  
  for (my $i= $ph0; $i<$nframe; $i++) {
    my $pro= translate( $cds, $i );
    my $flag=""; my $score=0;
    if (substr($pro,0,1) eq 'M') { $score += 1; $flag.="Start,"; }
    if (substr($pro,length($pro)-1,1) eq '*') { $score += 1; $flag.="Stop,";} # adj internal == end
    my $instop= substr($pro,0,length($pro)-1) =~ tr/*/*/;
    if($instop){ $score += $instop * -3; $flag.="Internalstop,";}
    elsif($score == 2) { $flag="Complete"; }
    #?? change score to |=1 for Start, |=2 for Stop, $instop * -3 ?
    if ($score > $bestscore) { $besti= $i; $bestscore=$score; $bestpro=$pro; $bestflag=$flag; }
    #ugh not this way# warn("# bestFrame[$i,$id]: $score ; $pro \n") if($debug); # debug  
    last if($score >= 2) # best possible; stop here
    }
  $bestflag ||="notrans"; #?
  return wantarray ? ($besti,$bestpro,$bestscore,$bestflag) : $besti;
}

1;

__END__

