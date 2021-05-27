#!/usr/bin/env perl
# evigene/scripts/genes/blast2evgff.pl

=item about

  convert table from 'blastn -query genes.fa -db genoasm -outfmt 7' to genes.gff  
  update of evigene/scripts/blast92gff3.pl
  
=item update 2017.06

  .. merge of blast92gff3.pl and gmap2evgff.pl (of GMAP -S output)
  .. added -introns introns.gff input for splices
  .. added -aasizes genes.aaqual for CDS-spans, CDS output in gff

  BUG: Target= offsets should be adjusted with intron adjusts to chr start/stop
  
=item fixme 2018.02
  ..  mrna cover is wrong,  need sum( exon align=\d+)/mrnalen
  ..  mix-strand introns near same exon end are problem (for some loci)
      .. need better way to choose which in matchIntronSpliceSites()
      
      
=head1 AUTHOR

  Don Gilbert, gilbertd@indiana.edu  June 2008 (1st blast92gff3.pl)
    
=cut

use strict;
use warnings;
use Getopt::Long;

my $VERSION = '2018.02.09'; # fix intron aligns
# '2017.06.25'; # "2011.09.02"; # new 2011.09
my $MAX_EXON_SEPARATION = 200000; 
my $MAX_ALIGN_PER_QUERY = 49; # new 2011.09; quickly drop low qual queries; speeds up tblastn prot-query x genome a lot; default? 
my $MININTRON=  20; # close exons w/ less separation == tr-align-gaps not introns.
  # NCBI sez SEQ_FEAT.ShortIntron  Introns should be at least 10 nt long
my $MINEXON= 20; # SEQ_FEAT.ShortExon  Internal coding region exon is too short .. but what is "too"?
my $MINIDENT=  0; # e.g. = 0.90, filter low qual exons to prevent them from bumping higher qual part.

use constant QBASEOVER => 11; # query hsp overlap slop in bases
use constant SBASEOVER => 9; # source hsp overlap slop in bases
use constant SPLICE   => 3; # bp for intron splice span tests
use constant SPLICEX  => SPLICE - 1; # for exon
use constant { kINTRON2SPLICE_OVER=>1, kINTRON2SPLICE_QUERY=>2, kINTRONERROR_OVER=>3, };
our $BINSIZE  = 1000 ; #was# 5000;

my $OVERLAP_SLOP_QUERY  = 0.15; # was 0.15; for protein query HSPs
my $OVERLAP_SLOP_GENOME = 0.05; # was 0.50; for genome source HSPs
my $BIN10000= 10000; # binsize
my $GAPSIZE = 400; # what ?? for nearover1() test

my $LOWSCORE_SKIP  = 0.80; #0.50 ; 0.33; # i.e. skip matchs < 50% of max match score

use vars qw(
$swap_querytarget
$faprefix $debug $introngff
$querySource $skipquery 
$exonType $geneType  $dotarget $doself $addgene
$max_bit_score $min_bit_score
$bitcutoff $stringent2ndary $onespecies
%sumhash @sourcehsps %moregenes 
$species $queryid $npart
$tophsp $lqid $lsid $nwarn $outfile $genesizes
$gffout %qlength
);

$geneType="mRNA"; # or "protein_match" ...; 
$exonType="exon"; # or "match_part"; 
$querySource="blast"; 
$faprefix= 'gnl\|'; # drop NCBI extra id prefix
$stringent2ndary=1;
$dotarget= 1; $doself= 1; $addgene= 1; # == add mRNA/match row, old syntax
$min_bit_score=0; $lqid=0;

my $DROP_UTR_CHIMERA=1; #  1703upd:
my $CUTUTRX= $ENV{cututr}||0; # chomp off excess utr exons, cutval=max 5/3 end UTRs; genejoin errs
my $dointron= (defined $ENV{intron}) ? $ENV{intron} : 2; # gff only
my $addcds  = (defined $ENV{cds}) ? $ENV{cds} : 1; # gff only
my $ADDBS = $ENV{bitscore}||0;  #UPD1811: ADDBS option add ann: bitscore, identlen to exons,mRNA
my $GSTRANDED= $ENV{strand}||0;  # when "genome" is stranded, like mRNA ref
my $KEEPATTR= $ENV{keepat}||"aalen,oid,offs,clen"; # 1703 default: keepat="aalen,oid,offs" * clen for nomap lengths
# my $MAX_SPAN_ERR=100; #? error if exon span genomic - transcript > this, too big an error
# my $MAX_SPAN_ERR_HIDE= 9999; # inner span error too big and cannot fix, hide this as #err
#   # 1703upd: change default noerrspan=1
# my $noerrspan= $ENV{noerrspan}||0;  # cancel badspan changes

use constant UPD19 => 1;   # see also addbitscore/ADDBS
my $annotfields=""; #UPD19 = option to parse # Fields:".. evalue, bit score, subject seq"
my %annotfields; # per query hash

#  do own sorting for score, loc : still need queryID sort ...
#  ** WARNING : input blast must be sorted by queryID, and best score (sort -k1,1 -k12,12nr) **
#  cat  modelproteins-mygenome.tblastn| sort -k1,1 -k12,12nr |\

my $USAGE =<<"USAGE";
   blast2evgff.pl [options]  geneseq_chrseq.blastn >  genes.gff
   NOTE : input blast table, sorted by transcript queryID:  'blastn -query genes.fa -db chrseq -outfmt 7'
 options: 
 -output genes.gff, instead of STDOUT
 -introns introns.gff, input to align exons to valid splice sites
 -inmax $MAX_EXON_SEPARATION : maximum intron gap for exons of one gene
 -[no]addCDS, make CDS from exons, given cds-spans from aasize table
 -aasizes protein-size/CDS-span table, as from evigene prots/aaqual.sh, need for add CDS output
 -source $querySource : gff.source field
 -exonType $exonType , -mrnaType $geneType  : gff feature types
 -minbitscore $min_bit_score : skip if bit_score below this
 -lowscore $LOWSCORE_SKIP : skip matches < $LOWSCORE_SKIP of max score
 -alignmax $MAX_ALIGN_PER_QUERY : skip excess aligns to same query
 -swap : swap query, target (as for blastn -query chrseq -db transcripts, but sort on col 2)

USAGE

my $USESIZE= $ENV{usesize}||0;

my $optok= GetOptions( 
  "output=s" => \$outfile,
  "sizes|aasizes=s" => \$genesizes,
  "USESIZE!" => \$USESIZE, 
  "qoverlap=s" => \$OVERLAP_SLOP_QUERY,  
  "overlap=s" => \$OVERLAP_SLOP_GENOME,  
  "introns=s" => \$introngff,
  "inmax=s" => \$MAX_EXON_SEPARATION, # was intronmax
  "alignmax=s" => \$MAX_ALIGN_PER_QUERY,
  "LOWSCORE_SKIP=s" => \$LOWSCORE_SKIP,  
  "minbitscore=s" => \$min_bit_score,
  "source|src=s", \$querySource, # also == GFF source ?
  "skipsource=s", \$skipquery,
  "mrnaType|geneType=s", \$geneType,  
  "exonType=s", \$exonType,  
  "swap_querytarget!" => \$swap_querytarget,
  'stringent2ndary!' => \$stringent2ndary, # unused now
  'onespecies!' => \$onespecies, # keep only 1 gene/species of target species:gene * UNUSED now
      ## requires some revision; easier to grep each spp from blast out 
  "faprefix=s", \$faprefix,
  "self!" => \$doself, # -noself or default
  "cdsadd|addcds!" => \$addcds,  
  "target!" => \$dotarget,
  "match!" => \$addgene,
  "debug!" => \$debug,
  "addbitscore!" => \$ADDBS,
  "annotfields=s", \$annotfields,

  "MININTRON=i", \$MININTRON,
  "MINEXON=i", \$MINEXON,
  "MINIDENT=s", \$MINIDENT, 
  "CUTUTRX=i", \$CUTUTRX,   # not ready
  "cututrchimera|droputrchimera!", \$DROP_UTR_CHIMERA,  # not ready, need cds-spans; -nodroputrchim to turn off
#?  "keepat=s", \$KEEPATTR, 
#?  "intron=i", \$dointron,  
#?  "nopath|unmapped!", \$NOPATHGFF,  #? if default=1, -nonopath or -nounmapped to turn off??
  );
die $USAGE unless($optok);

$LOWSCORE_SKIP=$LOWSCORE_SKIP/100 if($LOWSCORE_SKIP>1); # proportion
$MINIDENT= $MINIDENT*100 if($MINIDENT <= 1); # percent
#----------------------------------------

my($ok,$n_overlaps,$hascomm,$qcomm,$qqid,@evgattr);

sub reset_target_vars {

  printGeneLocations(@_); # turn hsps for 1 transcript query id into gene.gff, possibly several loci
  $npart= 0; 
  $max_bit_score= 1;
  $bitcutoff=0;
  $tophsp='';
  @sourcehsps=(); 
  %moregenes=();
  %annotfields=();
  %qlength=();
  #bad. if( my($at)= grep /clen=/, @evgattr ) { if($at =~ m/\bclen=(\d+)/) { $qlength{$qqid}= $1; } }
}

sub MAINstub {}

my $overlaplist=undef; # intronlist == overlaplist
if($introngff) {
  my $ovh; 
  if($introngff =~ /.gz$/) { $ok= open(OVR,"gunzip -c $introngff |");  $ovh= *OVR; }
  else { $ok= open(OVR,$introngff); $ovh= *OVR; }
  die "bad -overlaps=$introngff" unless($ok);
  $overlaplist= collect_overlaps($ovh); close($ovh);
  ## warn"#introns found=$n_overlaps\n" if $debug;
}

my($nsizes,$aasizeh,$trsizeh,$cdspanh,$aaqualh)
    = ($genesizes) ? readSizes($genesizes) : (0,0,0,0);   

#NO: my $inh = *STDIN; # option?
$gffout = *STDOUT;
if( $outfile && open(OUTH,">$outfile")) { $gffout= *OUTH; }
print $gffout "##gff-version 3\n#evigene genes/blast2evgff.pl $VERSION\n\n";

my $DEF_BLAST_FIELDS="# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score";
my @bfields= split(/, */, $DEF_BLAST_FIELDS); map{ s/ /_/g; s/%/p/g; s/\.//g; s/\W/_/g; } @bfields;

my @next_evgattr=(); # next Query info, keep till after next valid align row & reset_target_vars()
while(<>) {
  chomp;

  ## ? pull # Query comment for info
  ## Query: aedes_AAEL015287-PA desc=sodium/chloride dependent amino acid transporter; loc=supercont1.1712:8485:20120:1; gene=AAEL015287; MD5=0b6585f2b1568115937872f8fea07ced; length=638; dbxref=euGenes:ARP2_G144
  if(m/^# Query: (\S+)/) {  $qqid=$1;
    # upd17 for evigene headers aalen|clen|offs|oid
    # Query: Daplx7b3EVm013394t3 type=mRNA; aalen=293,76%,complete; clen=1159; offs=178-1059; oid=Daplx5cEVm012684t1,daplx5ad9c9tvelvk47Loc8t1; organism=Daphnia_pulex; 
    # also chrmap=95a,99i,1159l,6x,scaffold_30:206222-207774:+;
    @next_evgattr= m/\b((:?aalen|clen|offs|oid)=[^;\s]+)/g;
    if(m/\b(Name=[^;\n]+);/) { push @next_evgattr, $1;  }
    ## add namealn=$nap;Dbxref=$nadx;Name= ??
    if(@next_evgattr and not m/clen=/) { if( my($cl)= m/,(\d+)l,/ ) { push @next_evgattr, "clen=$cl"; } }
    }
  elsif(m/^#.Fields: (.+)/) { # UP19
    # Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score, subject seq
    my $fs=$1; 
    @bfields= split(/, */, $fs); map{ s/ /_/g; s/%/p/g; s/\.//g; s/\W/_/g; } @bfields;
    }
  elsif(/^\w/) {  # dont assume blast input has comments .. -m8
    my @brow= split "\t";
      #UPD19 want blast table some options? e.g. subject seq? parse # Fields header?
      # Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score, subject seq
    my ($qid, $sid, $pctident, $alignment_length, $mismatches, $gap_openings, 
       $q_start, $q_end, $s_start, $s_end, $prob, $bit_score ) = @brow; 
    next unless($bit_score); # some error logs mixed in.. ## warn, die ??
    next if($skipquery and $qid =~ m/$skipquery/); # patch to skip some parts
    
    cleanid($qid);  
    cleanid($sid);  
    if($swap_querytarget) {
      ($qid,$sid,$q_start,$q_end,$s_start,$s_end)= 
      ($sid,$qid,$s_start,$s_end,$q_start,$q_end);
      }
  
    if($qid eq $sid) {
      $qlength{$qid}= $alignment_length unless($qlength{$qid}); # if($pctident > 95);
      next unless($doself);
    }
    
    # query/gene batching, not target/scaffold, with sort -k1,1 input
    if($qid ne $lqid) {
      reset_target_vars($lqid);
      @evgattr= @next_evgattr; @next_evgattr=();
      if( my($at)= grep /clen=/, @evgattr ) { if($at =~ m/\bclen=(\d+)/) { $qlength{$qqid}= $1; } }
    }
    
    my($s_strand,$q_strand)= ('+','+');
    if ($s_start > $s_end) { $s_strand='-'; ($s_start,$s_end)= ($s_end,$s_start);  }
    if ($q_start > $q_end) { $q_strand='-'; ($q_start,$q_end)= ($q_end,$q_start);  }
    
    # count hsps same query loc, diff db loc
    $tophsp= $sid if ($npart==0);


    my $tkey= "$qid-HSP:$q_start,$q_end";  # query-loc key
    my $saved= 0;
    my $qoverlapped=0;  # Not used now: DROP from hspval
    my $soverlapped=0;  # DROP ; need only for genelocs.other ?

    # this is where we assume input is sorted by query-id and top bitscore : 
    # we drop overlapped hsps silently here ... should reinstitute stringent2ndary test
    #not here# next if( $qoverlapped && $soverlapped ); # this happens only for same query, same scaffold
    ##next if($stringent2ndary && $bit_score < $bitcutoff && ($sid ne $lsid));

    #o $alignment_length -= $mismatches; # FIX: 100207 ; UNFIX 2017.04, below uses alignlen as full span, identity tracks mismatches
    my $identlen = $alignment_length * $pctident / 100;
    $bit_score= int($bit_score); # has spacing
    
    # my $hspvalOLD= [$qid,$sid,$alignment_length, $q_start,$q_end,$s_start,$s_end, 
    #   $bit_score,  $prob, $tkey, $s_strand,
    #   $soverlapped, $qoverlapped, $pctident]; # 
    my $hspval= [ $qid, $sid, $alignment_length, $q_start,$q_end,$s_start,$s_end, 
                  $bit_score, $prob, $tkey, $s_strand, $identlen, $pctident ]; 
    my $atlockey="$qid.$sid.$q_start.$s_start.$s_end.$bit_score"; # == hsp-key, uniq to each pair/row, assignBestgene()
    push(@sourcehsps, $hspval); # unless($qoverlapped && $soverlapped); # use this as recycled hits/scaffold

if(UPD19){
    if($annotfields) {
      my $ann="";
      for(my $i=0; $i<=$#brow; $i++) { # slow, improve
        my $v=$brow[$i]; my $k=$bfields[$i];  
        if($k and $v and $annotfields =~ m/$k/) { $ann.= "$k=$v;"; }
        }
      $annotfields{$atlockey} = $ann if($ann);
    }
}    
    
    $npart++; $lqid= $qid; $lsid= $sid;
  } # blast line
} # while(<>)


reset_target_vars($lqid); ## printGeneLocations();
close($gffout) if($gffout);

if($debug) { # debug/verbose / stats of parts saved ...
  my($tpart,@tkeys,$v);
  my @sumkeys= sort keys %sumhash;
  warn "# Summary of HSPs saved\n";
  foreach $tpart (@sumkeys) {
    @tkeys= sort keys %{$sumhash{$tpart}};
    foreach (@tkeys) { $v= $sumhash{$tpart}{$_}; warn "# $tpart $_ = $v\n"; }
  }  
}

# end main // subs ...................................................

# hspvalOLD: ($qid,$sid,2:$alignment_length, 3:$q_start,$q_end, 5:$s_start,$s_end, 
#         7:$bit_score, 8: $prob, $tkey, 10:$s_strand,  $soverlapped, $qoverlapped, 13:$pctident) # drop: $soverlapped, $qoverlapped,
# hspvalNEW: ($qid, $sid, 2:$alignment_length, $q_start,$q_end, 5:$s_start,$s_end, 
#                   7:$bit_score, $prob, $tkey, 10:$s_strand, 11:$identlen, $pctident)

sub _sortHsp_Score {
  #old# return ($b->[7] <=> $a->[7]); # top bitscore 1st only? add align, pctident, ??
  return ($b->[11] <=> $a->[11] or $b->[7] <=> $a->[7] or $b->[2] <=> $a->[2]); # identlen>bits>alnlen   
}

sub _sortHsp_SrcLocation  {
  my($ar,$ab,$ae,$ao)= @{$a}[1,5,6,10]; # old/new hsp same
  my($br,$bb,$be,$bo)= @{$b}[1,5,6,10];
  my $ocmp= ($ao eq $bo)? 0 : ($ao eq "-") ? 1 : -1;
  return ($ar cmp $br || $ocmp || $ab <=> $bb || $be <=> $ae); # sort all -strand  last
}

sub bestlocHspset {
  my($hsploc, $tsid, $ts_start, $ts_end, $ts_strand)= @_;
  my @before=(); my @after=();
  my $skiphsp= 0;
  
  ## maxexonsep cuts out split-genes/chimera, fix?
  my($trange0, $trange1)= ($ts_start - $MAX_EXON_SEPARATION, $ts_end + $MAX_EXON_SEPARATION);
## FIXME: let valid introns override MAX_EXON_SEPARATION

  foreach my $hspval (@$hsploc) {
#     my($qid,$sid,$alignment_length, $q_start,$q_end,$s_start,$s_end, 
#         $bit_score,  $prob, $tkey, $s_strand,
#         $soverlapped, $qoverlapped, $pctident)
    my( $qid, $sid, $alignment_length, $q_start,$q_end,$s_start,$s_end, 
                  $bit_score, $prob, $tkey, $s_strand, $identlen, $pctident )        
        = @$hspval;
    
    #? mixed strand problems here?    
    next unless($sid eq $tsid && $s_strand eq $ts_strand);
    next unless($s_start > $trange0 && $s_end < $trange1);
 
    # filter out poor exons
    my $xw= 1 + abs($q_end - $q_start);
    my $xpoorident= ($MINIDENT and $pctident < $MINIDENT)?1:0; # dont have:  ($ixn>0 and $ixn < $#xon)
    my $xpoortiny=  ($MINEXON and $xw < $MINEXON)?1:0; # alignment_length or xw ?
    next if($xpoorident or $xpoortiny);
   
    if($s_end <= $ts_start + SBASEOVER) { 
      unshift(@before, $hspval); 
      $trange0 = _min($trange0, $s_start - $MAX_EXON_SEPARATION);
      }
    elsif($s_start >= $ts_end - SBASEOVER) { 
      push(@after, $hspval); 
      $trange1 = _max($trange1, $s_end + $MAX_EXON_SEPARATION);
      }
    }  
  return (\@before, \@after); # sorted around hspbest
}

sub assignBestgene {  # version 3
  #? my($theqid) = @_; # not used

  my($topsid,$saved,$genenum,$lastgenenum)=(0) x 9;
  $npart=0;
  
  ## for speed and for clarity, optional max genenum cutoff, ie. dont consider very low score single-hsp genes even if the hsps/gene sum up 
  my $lastexon= undef;
  my @allsaved=(); # dont need global
  
  my @hspbest = sort _sortHsp_Score @sourcehsps;
  my @hsploc  = sort _sortHsp_SrcLocation @sourcehsps;
  
  # input sourcehsp should all be for same query-gene, sorted by location (NOT/was bitscore)
  # need location-sort to match up exon parts of genes
  my %donehsp=();
#   my($qid,$sid,$alignment_length, $q_start,$q_end,$s_start,$s_end, 
#         $bit_score,  $prob, $tkey, $s_strand,
#         $soverlapped, $qoverlapped, $pctident);
  my( $qid, $sid, $alignment_length, $q_start,$q_end,$s_start,$s_end, 
                  $bit_score, $prob, $tkey, $s_strand, $identlen, $pctident );       

  my $max_bit_score= 1;
  my $genebest;
  
  foreach my $hspbest (@hspbest) {
  #   ($qid,$sid,$alignment_length, $q_start,$q_end,$s_start,$s_end, 
  #         $bit_score,  $prob, $tkey, $s_strand,
  #         $soverlapped, $qoverlapped, $pctident)
    ( $qid, $sid, $alignment_length, $q_start,$q_end,$s_start,$s_end, 
          $bit_score, $prob, $tkey, $s_strand, $identlen, $pctident )        
          = @$hspbest;
    my $topkey="$qid.$sid.$q_start.$s_start.$s_end.$bit_score";
  # warn "#top $topkey\n"; #DEBUG
    next if $donehsp{$topkey}++;

    ## FIXME: next overlapsome is problem if we want to keep all genes with ~same score at same loc
    ## next if overlapsome($s_start, $s_end, \@allsaved, $OVERLAP_SLOP_GENOME);
    my $isover= overlapsome($s_start, $s_end, \@allsaved, $OVERLAP_SLOP_GENOME);
    
    my $lowscore = ($bit_score < $max_bit_score * $LOWSCORE_SKIP); # 2010.3 patch
    next if($isover and $lowscore);
    $max_bit_score= $bit_score if ($bit_score > $max_bit_score);
    
    $genenum++;
    last if($MAX_ALIGN_PER_QUERY>0 and $genenum > $MAX_ALIGN_PER_QUERY);
    #NOT USED# $topsid= $sid; ## if ($npart==0);
     
    # if $onespecies
    # my($sdb,$sid1)= ($sid =~m/:/) ? split(/[:]/, $sid,2) : ("",$sid);
  
    my $keynum= $genenum;
    unless(ref $moregenes{$keynum}) { $moregenes{$keynum}= []; }
    push( @{$moregenes{$keynum}}, $hspbest); $saved=1;
    $sumhash{'other'}{($saved?'saved':'notsaved')}++; 
    if ($saved) {
      push(@allsaved, $hspbest);
      $sumhash{'ALL'}{($saved?'saved':'notsaved')}++; 
      }
    
    my($tsid, $tq_start, $tq_end, $ts_start, $ts_end, $ts_strand)= 
      ($sid, $q_start, $q_end, $s_start, $s_end, $s_strand);
    my($trange0, $trange1)= ($ts_start - $MAX_EXON_SEPARATION, $ts_end + $MAX_EXON_SEPARATION);
    my($srange0, $srange1)= ($ts_start, $ts_end);
    my($qrange0, $qrange1)= ($tq_start, $tq_end);
    
## FIXME: let valid introns override MAX_EXON_SEPARATION

      ## FIXME.3 this really needs to step thru @hsploc starting at $hspbest loc and go down,up from there
      ## otherwise qrange, srange are bad.  getting two genes made from very good pieces of 1 gene match
      ## due to interior hsp's skipped in first pass nearest to hspbest.

use constant HSPBEFOREAFTER => 1; #debug opt? is this helpful? check

    my($before, $after, @before,@after);
if(HSPBEFOREAFTER) {
    #d my @beh= reverse grep { $_->[3] < $tq_start } @hsploc; # before tq_start; bad??
    #d my @afh= grep{ $_->[3] >  $tq_start } @hsploc;

    my @beh= grep { $_->[5] < $ts_start } @hsploc; 
       @beh= reverse @beh; ## reverse @hsploc is bad ** need to UNREV after bestlocHspset
    my @afh= grep { $_->[5] > $ts_start } @hsploc;

    my($bbefore, $bafter)= bestlocHspset(\@beh, $tsid, $ts_start, $ts_end, $ts_strand);
    my($abefore, $aafter)= bestlocHspset(\@afh, $tsid, $ts_start, $ts_end, $ts_strand);
    
    @beh= reverse @$bbefore; push @beh, @$abefore; $before= \@beh;
    @afh= reverse @$bafter;  push @afh, @$aafter;  $after= \@afh;
    # $before= [ @$bbefore, @$abefore ]; #?? abefore should be empty,  for +strand / rev -strand
    # $after = [ @$bafter, @$aafter ]; #?? bafter should be empty..
    
} else {      
    ($before, $after)= bestlocHspset(\@hsploc, $tsid, $ts_start, $ts_end, $ts_strand);
}
    
    foreach my $hspval (@$after, @$before) { # instead of foreach can we hash-find nearby hsps?
      # next unless (ref $hspval); #?
      
  #     ($qid,$sid,$alignment_length, $q_start,$q_end,$s_start,$s_end, 
  #         $bit_score,  $prob, $tkey, $s_strand,
  #         $soverlapped, $qoverlapped, $pctident)
      ( $qid, $sid, $alignment_length, $q_start,$q_end,$s_start,$s_end, 
          $bit_score, $prob, $tkey, $s_strand, $identlen, $pctident )        
          = @$hspval;
  
      ## FIXME here; should not skip, but keep some of these to check; 
      ## should not look far away if nearby exon fits; it it is done already or overlaps, count in qrange and skip on
      ## FIXME.2 new problem with this; far-hsp already done can eat away query-range; need to skip those
  
      ## already done# next unless($sid eq $tsid && $s_strand eq $ts_strand && $s_start > $trange0 && $s_end < $trange1);
  
      my $skiphsp= 0;
      my $atlockey="$qid.$sid.$q_start.$s_start.$s_end.$bit_score";
      ## warn "#at $atlockey\n"; #DEBUG
      next if($atlockey eq $topkey);
      $skiphsp=1 if $donehsp{$atlockey};
      $skiphsp=1 if $skiphsp || overlapsome($s_start, $s_end, \@allsaved, $OVERLAP_SLOP_GENOME);
  
      my $qover= overlap1($q_start, $q_end, $qrange0, $qrange1, $OVERLAP_SLOP_QUERY);
      next if($qover); # last?
      # $skiphsp=1 if ($qover);
      ## FIXME: need to look at qloc vs top-qloc; if -strand, @before must be higher qloc, @after lower
      ## and v.v.
  
      # if($skiphsp) now check that s_start,s_end is *near* top hsp; skip if not
      if($skiphsp) {
        my $nearover= nearover1($s_start, $s_end, $srange0, $srange1, $OVERLAP_SLOP_GENOME, 5000); # GAPSIZE
        next if($nearover >= 0); ## { next if $skiphsp; }
      }
      
      if($s_end <= $ts_start + SBASEOVER) { # before
  #       if($s_strand eq '-') { $qrange1 = $q_end if($q_end> $qrange1);  }
  #       else {  $qrange0 = $q_start if($q_start< $qrange0); } # not before      
        if($s_strand eq '-') { next if($q_start + QBASEOVER < $qrange1); $qrange1 = $q_end if($q_end> $qrange1);  }
        else { next if($q_end - QBASEOVER > $qrange0);  $qrange0 = $q_start if($q_start< $qrange0); } # not before      
        unshift(@before, $hspval) unless $skiphsp; 
        $srange0 = $s_start unless $skiphsp;
        }
        
      elsif($s_start >= $ts_end - SBASEOVER) { # after; bug in next here
  #       if($s_strand eq '-') {  $qrange0 = $q_start if($q_start < $qrange0); }
  #       else {  $qrange1 = $q_end if($q_end > $qrange1);  } # not after
        if($s_strand eq '-') { next if($q_end - QBASEOVER > $qrange0); $qrange0 = $q_start if($q_start< $qrange0); }
        else { next if($q_start  + QBASEOVER < $qrange1);  $qrange1 = $q_end if($q_end> $qrange1);  } # not after
        push(@after, $hspval) unless $skiphsp; 
        $srange1 = $s_end unless $skiphsp;
        }
    ##warn "#skip=$skiphsp $atlockey\n"; #DEBUG
    }  
    
    #? limit before, after size? dont try to find what isn't there by too far a match
    foreach my $hspval (@before, @after) { # .. and after
    
  #     ($qid,$sid,$alignment_length, $q_start,$q_end,$s_start,$s_end, 
  #         $bit_score,  $prob, $tkey, $s_strand,
  #         $soverlapped, $qoverlapped, $pctident)
      ( $qid, $sid, $alignment_length, $q_start,$q_end,$s_start,$s_end, 
          $bit_score, $prob, $tkey, $s_strand, $identlen, $pctident )        
          = @$hspval;
      my $atlockey="$qid.$sid.$q_start.$s_start.$s_end.$bit_score";
      $saved= 0; ##no## $genenum= 0;
         
      if(1) {
        #my $keynum= "G" . $genenum;
        my $keynum=  $genenum;
        unless(ref $moregenes{$keynum}) { $moregenes{$keynum}= []; }
        push( @{$moregenes{$keynum}}, $hspval); $saved=1;
        $sumhash{'other'}{($saved?'saved':'notsaved')}++; 
      }
  
      if ($saved) {
        $lastexon= $hspval;
        $lastgenenum= $genenum;
        $donehsp{$atlockey}++;
        push(@allsaved, $hspval);
        $sumhash{'ALL'}{($saved?'saved':'notsaved')}++; 
        }
     $npart++;
    }
  }
  return($genenum); # not == ngene ? not saved
}

# sub parseCDSoff {
#   my($geneid, $cdsoff, $g_strand)= @_;
#   return (0,0) unless($cdsoff);
#   my($cdsb,$cdse)= $cdsoff =~ m/(\d+)\-(\d+)/; # 199-339:+; 143-3:- << note rev locs:strand
#   return (0,0) unless($cdse>0);
#   my $cdsrev=0;
#   if($cdse < $cdsb) { $cdsrev=-1; ($cdsb,$cdse)=($cdse,$cdsb); }
#   my $cdstrlen= 1+$cdse-$cdsb; # global for putalign!
#   return($cdsb,$cdse,$cdsrev,$cdstrlen);
# }

sub addCDS {
  my($geneid, $cdsoff, $blaststrand, $exongff)= @_;
  
  ## expect caller to find cdsoff..
  # ($cdsoff)= grep /offs=/, @evgattr unless($cdsoff); # or where from?
  # $cdsoff= $cdspanh->{$id} if($genesizes and ref($cdspanh));

  return () unless($cdsoff and ref($exongff));
  my($cdsb,$cdse)= $cdsoff =~ m/(\d+)\-(\d+)/; # 199-339:+; 143-3:- << note rev locs:strand
  my $cdsrev=0;
  return () unless($cdse);
  my @cdsgff=();
  if($cdse < $cdsb) { $cdsrev=-1; ($cdsb,$cdse)=($cdse,$cdsb); }
  my $cdstrlen= 1+$cdse-$cdsb; # global for putalign!
  for my $xn (@$exongff) {
    my @x=split"\t",$xn;
    my($xc,$xt,$xb,$xe,$xo,$at)= @x[0,2,3,4,6,8];
    next unless($xt eq $exonType);
    $at||=""; 
    my($id)= $at =~ m/Parent=([^;\s]+)/; $id=$geneid unless($id); # shouldnt be needed
    my($tid,$tb,$te)= $at=~m/Target=(\S+).(\d+).(\d+)/;
    # ($tb,$te)= ($te,$tb) if($tb > $te);
    
    my $orientaln= $xo; # IS THIS strand ok? need introns or cds orient to know.
    if($orientaln ne $blaststrand) { # swap? *** NEED TO CHECK this effect, seems to work, unclear
      $orientaln= $blaststrand;
    }
    
    if($tb < $cdse and $te > $cdsb) { 
      # my($xb,$xe)= @g[3,4]; # genomic start,end
      my @cds= @x;  
      $cds[2]= "CDS"; $cds[8]= "Parent=$id"; 
      unless($tb >= $cdsb and $te <= $cdse) {
        if($tb < $cdsb) { my $d= $cdsb - $tb; 
          if($orientaln eq "-") { $xe -= $d; } else { $xb += $d;} }
        if($te > $cdse) { my $d= $te - $cdse; 
          if($orientaln eq "-") { $xb += $d; } else { $xe -= $d;} }
          # DAMN cds ends calc bug ** DUE to exon mismap Target span > Genome span;
        #  need to ? reduce diff ; cancel this CDS if genome span too short? or partial
        if($xb >= $xe) { @cds=(); } 
        else { @cds[3,4]= ($xb,$xe); }
      }
     ## need to calc cdsglen from all cds exons in @aligno
     ## unless($cdsglen == $cdstrlen) { my $d= $cdstrlen - $cdsglen; $xout[0]=~s,$,;cdsindel=$d,; } # :$cdstrlen/$cdsglen
     if(@cds) { push @cdsgff, join("\t",@cds)."\n"; } # $nput++;
    }
  }

  return(@cdsgff); # also $cdsrev * orientaln ?
}


=item overlap_intronOLD

sub overlap_intronOLD
{
  my($ref,$tb,$te,$to,$tid)= @_;   # this is exon
  return 0 unless($overlaplist->{$ref});

  my $stranded= 0; # not for this data; make option..
  my(%didid, @hits);
  my ($ib1, $ib2)= (int($tb/$BINSIZE) , int($te/$BINSIZE));
  for (my $ib = $ib1; $ib <= $ib2; $ib++)  {
    $overlaplist->{$ref}{$ib} or next;
    my @locs= @{$overlaplist->{$ref}{$ib}};
    foreach my $rloc (@locs) {
      #  my $rloc= [$tb,$te,$ref,$gid,$to,$typ,$oid,$sat]; # change to string; save mem
      ## sat == splice at this point .. 1st/last exon base before/after intron splice
      my ($lb,$le,$lid,$lo,$oid,$sat)= @{$rloc}[0,1,3,4,6,7]; # overlap record, not hspval
      next if($didid{$oid.$lb.$le}++);
      
      my $over= ($tb <= $le && $te >= $lb) ? 1 : 0;
      # $over=0 if($stranded and ($to =~ /[+-]/ and $lo =~ /[+-]/ and $to ne $lo)); # NEW 2010.10; fixme for to or lo eq "."
      if($over) {
        # do we collect all or just 1st?
        my $mind= _min( abs($tb - $sat), abs($te - $sat));
        my $retval= [$sat,$sat,$ref,$lid,$lo,$oid,$mind]; # sat == splice at, exon-end point, 1bp before/after intron
        #x my $retval= [$lb,$le,$ref,$lid,$lo,$oid,$mind]; # NOTE these are splice points, lb-le,  3bp wide
        
        push @hits, $retval; # return $retval;
      }
    }
    last if (@hits);
  }
  if(@hits) {
    my($bhit)= sort{ $a->[6] <=> $b->[6] or $a->[0] cmp $b->[0] } @hits;
    return $bhit;
  }
  return 0;
}

=cut

      
# my $esplice= overlap_intronNEW($sid, $s_end, 1, $s_strand, $qid ); # TEST 18.02,
# this works, matches splice sites gmap.gff better
sub overlap_intron
{
  my($ref,$tpoint,$tIsEnd,$to,$tid)= @_;   # this is exon
  return 0 unless($overlaplist->{$ref});
  my($tb,$te)= ($tpoint - QBASEOVER, $tpoint + QBASEOVER);
  
  my $stranded= 0; # not for this data; make option..
  my(%didid, @hits);
  my ($ib1, $ib2)= (int($tb/$BINSIZE) , int($te/$BINSIZE));
  my $onepass= ($ib1 == $ib2);
  for (my $ib = $ib1; $ib <= $ib2; $ib++)  {
    $overlaplist->{$ref}{$ib} or next;
    my @locs= @{$overlaplist->{$ref}{$ib}};
    foreach my $rloc (@locs) {
      # $rloc= [$tb,$te,$ref,$gid,$to,$typ,$oid,$sat,$tp]; # 18.02 added score=tp; change to string; save mem
      ## sat == splice at this point .. 1st/last exon base before/after intron splice
      
      my ($lb,$le,$lid,$lo,$oid,$sat,$iscore)= @{$rloc}[0,1,3,4,6,7,8]; # overlap record, not hspval
      next if($didid{$oid.$lb.$le}++);
      
      my $over= ($tb <= $le && $te >= $lb) ? 1 : 0;
      if($over) {
        # do we collect all or just 1st? NOT 1st, closest to $tpoint, prefer tb side for tIsEnd, te side for not end
        #old my $mind= _min( abs($tb - $sat), abs($te - $sat));
        my $mind=  abs($tpoint - $sat); # this is key min distance from aligned exon end
        my $retval= [$sat,$sat,$ref,$lid,$lo,$oid,$mind,$iscore]; # sat == splice at, exon-end point, 1bp before/after intron

        push @hits, $retval; # return $retval;
      }
    }
    last if($onepass); # last if (@hits);
  }
  
  if(@hits) {
    ## fixme: mixstrand problem, use $to AND $iscore to help pick top hit?
    my($bhit,$chit)= sort{ $a->[6] <=> $b->[6] or $a->[0] cmp $b->[0] } @hits;
    if($chit and $bhit->[4] ne $to and $chit->[4]  eq $to) { $bhit=$chit; } # 18.02 eqstrand wins
    return $bhit;
  }
  return 0;
}


sub matchIntronSplices {
  my($exons,$blaststrand) = @_;
  my $nalign=0;
  #   my $blaststrand= $exons->[0]->[10]; # for CDS strand mixup; but see g_strand/intrstrand

  ## Ugh.. *should* revise this to match  adjacent exon ends to one intron, not to intron splice ends as here..
  my $lq_end=0;
  
  my @ex=  sort{ $$a[5] <=> $$b[5] } @$exons; # sort by $s_start
  for(my $i=0;  $i < @ex; $i++) {
    my $ex= $exons->[$i];
    my( $qid, $sid, $alignment_length, $q_start,$q_end,$s_start,$s_end, 
        $bit_score, $prob, $tkey, $s_strand, $identlen, $pctident )      
      = @$ex;
      
    my $ashift=0;
    my($new_start,$new_end,$new_strand,$inhit)=($s_start,$s_end,$s_strand,"");
    ## s_strand is bogus..
 
    ##18.02 ..  mix-strand introns near same exon end are problem (for some loci)
    ##18.02  .. need better way to choose: cds-strand best; else mrna strand or max of exons, but s_strand is/maybe invalid here
   
    # FIXME: start-slop, end + slop are not prefered (ie. expand align to unknown, other contract to known align)
    # .. ?? try first start+slop,end-slop, then +/- if missed
    # .. also overlap_intron() can hit 2+ splices w/i slop, pick closest to align point
    
    #orig# my $bsplice= overlap_intronOLD($sid, $s_start, $s_start + QBASEOVER,$s_strand,$qid );     
    my $bsplice= overlap_intron($sid, $s_start, 0, $blaststrand, $qid );  # TEST 18.02, missing some offby1,2,3

    if($bsplice) {
      my($lb,$le,$ref,$lid,$lo,$oid)= @$bsplice;
      $new_start= $lb; # ($lo eq "-")?$lb:$le; # need lo?
      $new_strand= $lo; 
      $inhit.="$lid,";
    }
     
    #orig# my $esplice= overlap_intronOLD($sid, $s_end - QBASEOVER, $s_end, $s_strand, $qid );
    my $esplice= overlap_intron($sid, $s_end, 1, $blaststrand, $qid ); # TEST 18.02,
    if($esplice) {
      my($lb,$le,$ref,$lid,$lo,$oid)= @$esplice;
      $new_end= $lb; # ($lo eq "-")?$le:$lb; # need lo?
      $new_strand= $lo; 
      $inhit.="$lid,";
    }
    
    if($bsplice or $esplice) {
      my $ds= $s_start - $new_start;
      my $de= $new_end - $s_end;
      my $dsb= $ds + $de;
      my $dsv= ($bsplice)?$ds:".";
      my $dev= ($esplice)?$de:".";
 
      # if($s_strand eq "-") {  ($ds,$de)=($de,$ds); } # what? change q sign?
      
      ##?? not this way?
      # my $new_qstart= $q_start + $ds;
      # my $new_qend= $q_end + $de;
      my($new_qstart,$new_qend)= ($q_start,$q_end);
      my $newaln = 1 + $q_end - $q_start;
      
      # $newaln = $dsb + $alignment_length; #?
      # $newaln = 1 + $new_end - $new_start;  # this *should* be okay now..
      
      if($lq_end and $q_start <= $lq_end) { #** THIS IS/was BAD ; problem below fixed it?
        my $oldaln = 1 + $q_end - $q_start;
        my $d= $lq_end - $q_start;
        if($d > 0 and $d <= QBASEOVER) {
          $new_qstart= $q_start + $d;
          # $new_qend= $new_qstart + $oldaln -1;  
          if($newaln > $oldaln - $d) { $newaln -= $d; }
          # $new_qend= $new_qstart + $newaln -1;  #? no
          }
        }
      # $new_qend= $new_qstart + $newaln -1;  # this *should* be okay now..
      # my $newaln = 1 + $new_qend - $new_qstart;
      
      $nalign++;
      my $xflag= ($debug)?"introna=$dsv,$dev,$inhit":"introna=$dsv,$dev";
      $ex= [ $qid, $sid, $newaln, $new_qstart, $new_qend, $new_start, $new_end, 
            $bit_score, $prob, $tkey, $new_strand, $identlen, $pctident, $xflag ]; # intron aln flag? i13
      $exons->[$i]= $ex;
    }
    
    $lq_end= $ex->[4];
  }
  return($nalign);
}

sub printOneLocation {
  my($exons, $igene, $ngene) = @_;
  return 0 unless(ref $exons && @$exons > 0);

  #UPD1811: ADDBS option add ann: bitscore, identlen to exons,mRNA
  my $qdb=".";
  my @locs= ();
  my($sum_bit_score,$sum_align,$sum_ident,$iex)=(0)x10;
#   my($qid1,$qid,$sid,$alignment_length,$q_start,$q_end,$s_start,$s_end,
#       $bit_score, $prob, $tkey,$s_strand, $soverlapped, $qoverlapped, $pctident);
  my($qid1);
  my( $qid, $sid, $alignment_length, $q_start,$q_end,$s_start,$s_end, 
        $bit_score, $prob, $tkey, $s_strand, $identlen, $pctident, $xflag );        

  my @gffout= ();
  my($g_start,$g_end, $g_strand,$tm_start,$tm_end,)=(undef,0,0,0,0,0,0);
  my($lq_start,$lq_end)=(0,0);
  my($ifwd,$irev,$inostrand)=(0,0,0);
  my $nexon= @$exons;
  
  my $blaststrand= $exons->[0]->[10]; # for CDS strand mixup; but see g_strand/intrstrand
  
  my($intralign)= ($overlaplist) ? matchIntronSplices($exons,$blaststrand) : 0;
  
  foreach my $ex ((sort{ $$a[5] <=> $$b[5] } @$exons)) # sort by $s_start
  {
#     ($qid,$sid,$alignment_length, $q_start,$q_end,$s_start,$s_end, 
#       $bit_score,  $prob, $tkey, $s_strand,
#       $soverlapped, $qoverlapped, $pctident)
    ( $qid, $sid, $alignment_length, $q_start,$q_end,$s_start,$s_end, 
        $bit_score, $prob, $tkey, $s_strand, $identlen, $pctident, $xflag )      
      = @$ex;
    
    my $atlockey="$qid.$sid.$q_start.$s_start.$s_end.$bit_score"; # UPD19 hsp/exon key
    $qdb= $querySource; # default
    ($qdb,$qid1)= ($qid =~m/:/) ? split(/[:]/, $qid,2) : ($qdb,$qid);
    my $inaligned= ($xflag and $xflag=~/introna/)?1:0;
    ## ** add mixstrand check/flag for  inaligned ..
    
    #upd17: trim exon overlap slop, *should* look for proper splice sites on chr.fasta . later
    #upd: need to know if exon aligned to introns here?
    if( not $inaligned # ** cancel shift if have intron align
      and $lq_end and $q_start <= $lq_end) {
      my $d= 1 + $lq_end - $q_start;
      if($d > 0 and $d <= QBASEOVER) {
        #? change alignlen, other
        my $qaln= 1 + $q_end - $q_start;
        $q_start += $d;
        $s_start += $d;
        if($alignment_length > $qaln - $d) { $alignment_length -= $d; }
      }
    }
    
    $g_start= $s_start unless(defined $g_start);
    $g_end= $s_end;
    
    if($inaligned) { if($s_strand eq '-'){ $irev++; } elsif($s_strand eq '+'){ $ifwd++; }  }
    else { $inostrand++; }
    if($g_strand) {
      $g_strand=$s_strand if($inaligned);
    } else {
      $g_strand= $s_strand ;
    }
      
    #??was# $tm_start= $q_start unless($tm_start or $tm_start > $q_start); 
    $tm_start= $q_start unless($tm_start > 0 and $q_start > $tm_start); 
    $tm_end  = $q_end   unless($q_end < $tm_end);

    $iex++;
    $sum_bit_score += $bit_score;
    $sum_align += $alignment_length; #getting > tlength, exon overlap slop?
    $sum_ident += $pctident * $alignment_length; ## pctident is 100% not 1.00, need pctident/100 * alen? no, see use below
    if(1) { # == GFF_LOC
      $tkey =~ s/^\w+://;  $tkey =~ s/,/-/g; # extra dbid and no commas in values
      #was# $attr= "Parent=${qid1}_${igene};tkey=$tkey;tloc=$q_start-$q_end;align=$alignment_length";
      my $tid=$qid1; $tid .= "_G$igene" if($igene>1);
      my $attr= "Parent=$tid";
      $attr.= ";Target=$qid $q_start $q_end;align=$alignment_length" if($dotarget);
      $attr.= ";idlen=".int($identlen).";bits=".int($bit_score) if($ADDBS);
      ## add intron align id/flag ?
      $attr.=";$xflag" if($xflag); # if(my $xflag=$ex->[13]){ $attr.=";$xflag"; }

# use constant UPD19 => 0;   
if(UPD19){
    if($annotfields) {
      my $anf= $annotfields{$atlockey}||"";
      if($anf){ $attr.=";$anf"; } #? formatted
    }
}    
      
      # print $gffout 
      my $xval= int(0.5 + $pctident); # upd17: was bit_score
      push @gffout, join("\t",
        ($sid, $qdb, $exonType, $s_start, $s_end, $xval, $s_strand,".",$attr)). "\n";
      }
    # else { push(@locs,"$s_start..$s_end"); }
    ($lq_start,$lq_end)=($q_start,$q_end);
  }
    
  $max_bit_score = $sum_bit_score if($sum_bit_score > $max_bit_score ); # assumes output by best match 1st
  return 0 if ($sum_bit_score < $max_bit_score * $LOWSCORE_SKIP);
  # return 0 if ($sum_bit_score < $min_bit_score); # min == 1 default
  
  ## FIX mixstrand and nostrand exons if have introns for some.
  if($intralign) {
    if($ifwd and $irev) { $g_strand= ($irev>$ifwd)?"-":"+"; } # ties? should use strand of cds-span as best
    $blaststrand=$g_strand if($blaststrand ne $g_strand);  #? is this right?
    if($inostrand) {
      for(my $i=0; $i<@gffout; $i++) {
      my @xv= split"\t",$gffout[$i];
      # my $inaligned= ($xflag and $xflag=~/introna/)?1:0;
      unless($xv[8]=~/introna=/) { $xv[6]= $g_strand; $gffout[$i]=join"\t",@xv; }
      }
    }
  }
  
  my $cdsoff ="";
  # my($cdsb,$cdse,$cdsref,$cdslen)= parseCDSoff($qid,$cdsoff,$blaststrand);
  
  if(1) { # check for regardless of addcds;
    ($cdsoff)= grep /offs=/, @evgattr; # or where from?
    if(not $cdsoff and $genesizes and ref($cdspanh)) {
      $cdsoff= $cdspanh->{$qid} || $cdspanh->{$qid1};
    }
    if($cdsoff and $addcds) {
      # STRAND MIXUP: blastn tr align has 1 orient, intronsplices have other, CDS needs to know both?
      my @cds= addCDS($qid1, $cdsoff, $blaststrand, \@gffout);
      # my ($revcds,$cds_gstrand,@cds)= addCDS($qid1, $cdsoff, $blaststrand, \@gffout);
      push @gffout, @cds if(@cds);
    }
  }
  
  if($addgene) { # should be always here
    my $tid=$qid1; $tid .= "_G$igene" if($igene>1);
    my $attr= "ID=$tid"; #? drop _G1 tag as for gmap.gff ?
    
    #b my $tspan=1 + $tm_end - $tm_start; #?? NOW same as sum_align which DID remove mismatches in align span
    #    #^1802 bad calc, need exon aligns + xtarget b e to find missing exon aligns
    my $tspan=1 + $tm_end - $tm_start; # use this for qlen, match=?
    #x my $talign= $sum_align; # 18.02 revert
          
    my $pident= ($sum_align<1)? 0 : int(0.5 + $sum_ident/$sum_align);
    my $tlen= $qlength{$qid} || 0;
    #?? if(not $tlen and $genesizes and ref($trsizeh)) { $tlen=  $trsizeh->{$qid} || $trsizeh->{$qid1}; } # trust this?
    if($genesizes and ref($trsizeh) and ($USESIZE or not $tlen)) { $tlen=  $trsizeh->{$qid} || $trsizeh->{$qid1}; } # trust this?
    
    $attr.= ";Target=$qid $tm_start $tm_end" if($dotarget);
    if($tlen) {
      my $cov=int(0.5 + 100*$sum_align/$tlen); $cov=100 if($cov>100);
      $attr.= ";qlen=$tlen;cov=$cov";# ;#upd17: tlen > qlen
    }
      
# OPTION: output filter by MIN_COV, MIN_PID ?

    $attr.= ";pid=$pident;match=$tspan;nexon=$nexon";#upd17: change pident= to pid=, add cov=? targ_span/targlen clen=
    $attr.= ";idlen=".int(0.5 + $sum_ident/100).";bits=".int($sum_bit_score) if($ADDBS);
    if($intralign) { 
      $attr.=";intralign=$intralign"; 
      if($ifwd and $irev) { 
        $g_strand=($irev>$ifwd)?"-":"+"; # done above now
        $attr.=";strandmix=+$ifwd/-$irev"; 
      }
    }
     
    #  @evgattr= m/\b((:?aalen|clen|offs|oid)=[^;\s]+)/g; # may be, or not, avail
    # aalen=105,55%,complete;offs=109-426;
    my($aaq)= grep /aalen=/, @evgattr;  
    if(not $aaq and $genesizes and ref($aaqualh)) { $aaq= $aaqualh->{$qid} || $aaqualh->{$qid1}; }
    if($aaq){ $aaq=~s/^\w+=//; $attr.= ";aalen=$aaq"; }
    if($cdsoff) { $cdsoff =~ s/^\w+=//; $attr.= ";offs=$cdsoff"; }
    my($oid)= grep /oid=/, @evgattr;  
    if($oid){ $oid =~ s/^\w+=//; $attr.= ";oid=$oid"; }
    if( my($name)= grep /Name=/, @evgattr ) { $attr.= ";$name";}
    
    if($ngene>1) { $attr.= ";path=$igene/$ngene"; }
    my $xval= $pident || 1; # upd17: was sum_bit_score
    unshift @gffout, join("\t",
      ($sid, $qdb, $geneType, $g_start, $g_end, $xval, $g_strand,".",$attr)). "\n";
  }
  
  print $gffout @gffout;
  return 1;
}

sub printGeneLocations {
  my($qid) = @_; # not used
  
  # ** mem overload problem hash grows to end:  %didid # this is check only for sort error; drop?
  # if(@sourcehsps && $qid && 0 < $didid{$qid}++) { die "ERROR: $qid already seen\n$USAGE"; }
  
  #** FIXME: problem here for next best ~= best score
  my($NOTngene)= assignBestgene(); # NEW; this could be buggy ?**
 
  #? add where? chimera checks; if($ischim and $DROP_UTR_CHIMERA) ..
  #?add here/where?  @exons= cututrx($CUTUTRX,@exons) if($CUTUTRX); # upd1703
 
  my $ng= 0;
  my @genes= sort{$a<=>$b} keys %moregenes;
  my $ngene= @genes;
  foreach my $keynum (@genes) { ## NOW NUMERIC KEY
    (my $gnum= $keynum) =~ s/^(.)[a-z]+/$1/;  ## S2..S9; o2..o9
    $ng += printOneLocation( $moregenes{$keynum}, $gnum, $ngene);
    }
  return $ng;
}


sub overlapsome { # need scaffold?
  my($s_start, $s_end, $exons, $OVERLAP_SLOP)= @_; 
  foreach my $ex (@$exons) {
    return 1 if overlap1($s_start, $s_end, $$ex[5], $$ex[6], $OVERLAP_SLOP);
  }
  return 0;
}

sub overlap1 {
  my($s_start, $s_end, $b_start, $b_end, $OVERLAP_SLOP)= @_; 
  return 0 unless ( $s_start < $b_end && $s_end > $b_start  ); ## no overlap
  return 1 if ( $s_start >= $b_start && $s_end <= $b_end ); # contained-in
  if ( $s_end > $b_start && $s_start < $b_end ) {  
    ## e.g.  s= 20,50  ; b = 10,30 : o= 10* or 40
    my $b_len= 1 + $b_end - $b_start;
    my $s_len= 1 + $s_end - $s_start; # choose ? biggest
    my $olp1 =  $b_end - $s_start; # 30-20 = 10
    my $olp2 =  $s_end - $b_start; # 50-10 = 40
    $olp1= $olp2 if ($olp2<$olp1);
    $s_len= $b_len if ($s_len<$b_len);
    return 1 if (($olp1 / $s_len) > $OVERLAP_SLOP);
    }
  return 0;
}  

sub nearover1 {  # near == -1; over == +1; other == 0
  my($s_start, $s_end, $b_start, $b_end, $OVERLAP_SLOP, $gapsize)= @_; 
  
  # need some overlap slop here; note locs here are QUERY == protein; should not be large
  unless ( $s_start < ($b_end - QBASEOVER) && $s_end > ($b_start + QBASEOVER)  ) {
    if($s_start >= $b_end) { return ($s_start - $b_end < $gapsize) ? -1 : 0; }
    elsif($b_start >= $s_end) { return ($b_start - $s_end < $gapsize) ? -1 : 0; }
    return 0;
  }
  
  return overlap1($s_start, $s_end, $b_start, $b_end, $OVERLAP_SLOP);
}  

sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }

sub minabs {
  my($a,$b)= @_;
  $a= abs($a); $b= abs($b);
  return ($a>$b)? $b : $a;
}

sub toofar1 {
  my($q_start, $q_end, $b_start, $b_end)= @_; 
  return 1 if minabs($b_end - $q_start,$q_end - $b_start) > $MAX_EXON_SEPARATION;
  return 0;
}

sub cleanid { 
  unless( $_[0] =~ s/$faprefix//) { $_[0] =~ s/^gi\|\d+\|(\S+)/$1/; }  # drop gi nums for ids
  $_[0] =~ s/\|$//; ## some ids have '|' at end of id; chomp
  $_[0] =~ s/\|/:/; #?? change pipe to db:id format
}

sub collect_overlaps
{
  my($inh)= @_;
  my %overlaps=(); my $nr=0;
  my $passtypes=""; # option, intron type
  #which? my $intron2splice= kINTRON2SPLICE_QUERY; #? 
  my $intron2splice= kINTRON2SPLICE_OVER; # option? or only these?
  
  while(<$inh>){
    next unless(/^\w/); chomp;
    my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr)= split"\t"; #,@gffmore
    $tattr ||="";  
    
    # intron here..
    if($passtypes and "$typ.$src" !~ m/$passtypes/) { next; } # pass other types
      # ^ do in both filter_gff and in collect_overlaps **??
    
    $nr++;
    my $oid= "N".$nr; # save always for uniq feature test
    my($gid); # fixme: separate $gid from markid : want two fields, one for id tests
    # if($markidtype) {
    #   if($intron2splice and $markidtype =~ /^score|orient|strand/i) {
    #      $gid = "$to$tp"; # strand.score = +nnn,-nnn,.nnn,0nnn possible results
    #   }
    #   elsif($markidtype =~ /^score/i) { $gid=$tp; }
    #   elsif($markidtype =~ /^source/i) { $gid=$src; }
    #   elsif($markidtype =~ /^type/i) { $gid=$typ; }
    #   elsif($markidtype =~ /^orient|strand/i) { 
    #       my $ov=($to eq "-")?"-1":($to eq "+")?"+1":($to eq ".")?"0":$to;
    #       $gid=$ov; } # change to +1,-1,0
    #   elsif($tattr =~ m/\b$markidtype=([^;]+)/) { $gid=$1; }
    # } else
    if($tattr) {
      if($tattr =~ m/\bID=([^;\s]+)/) {  $gid=$1; }
      elsif($tattr =~ m/\bParent=([^;\s]+)/) {  $gid=$1; }
      elsif($tattr =~ m/\bTarget=([^;\s]+)/) { $gid=$1; } # NEW 2010.10
    }
    unless(defined $gid) { $gid = $oid; }

    # $ordered_over{$gid}= $nr; # if $orderedoverlap;
    my $sat= $tb;
    my $rloc= [$tb,$te,$ref,$gid,$to,$typ,$oid,$sat,$tp]; # 18.02 added score=tp; change to string; save mem
    
    if($intron2splice == kINTRON2SPLICE_OVER or $intron2splice == kINTRONERROR_OVER) { # 2010jul
      my($s1b,$s1e,$s2b,$s2e)= 
        ($to eq "-") ? ($te+1,$te + SPLICE,$tb - SPLICE,$tb-1) : ($tb - SPLICE,$tb-1,$te+1,$te + SPLICE); # 3bp + 1 shift
      # ($tb,$te)= ($s1b,$s1e);
      $sat= ($to eq "-") ? $te+1 : $tb-1;
      $rloc= [$s1b,$s1e,$ref,$gid,$to,$typ,$oid,$sat]; 
      my @bins= (int($tb/$BINSIZE) .. int($te/$BINSIZE));
      foreach my $ib (@bins) { push( @{$overlaps{$ref}{$ib}}, $rloc); }
      
      # ($tb,$te)= ($s2b,$s2e);  #? change gid, oid?
      $sat= ($to eq "-") ? $tb-1 : $te+1;
      $rloc= [$s2b,$s2e,$ref,$gid,$to,$typ,$oid,$sat]; 
      #NOT NOW# $nr++; $oid= "N".$nr;  # dang, should oid be same for same intron, or +1 ?
      
    } elsif($intron2splice == kINTRON2SPLICE_QUERY) {
      # here is exon over, need to trim to exon endpoints to test overlap at ends?
      my($s1b,$s1e,$s2b,$s2e)= 
        ($to eq "-") ? ($te - SPLICEX,$te,$tb,$tb + SPLICEX) : ($tb,$tb + SPLICEX,$te - SPLICEX,$te); # 3 bases of exon
      # ($tb,$te)= ($s1b,$s1e);
      $sat= ($to eq "-") ? $te : $tb;
      $rloc= [$s1b,$s1e,$ref,$gid,$to,$typ,$oid,$sat]; 
      my @bins= (int($tb/$BINSIZE) .. int($te/$BINSIZE));
      foreach my $ib (@bins) { push( @{$overlaps{$ref}{$ib}}, $rloc); }
      # ($tb,$te)= ($s2b,$s2e);  #? change gid, oid?
      $sat= ($to eq "-") ? $tb : $te;
      $rloc= [$s2b,$s2e,$ref,$gid,$to,$typ,$oid,$sat]; 
    }
      
    my ($ib1, $ib2)= (int($tb/$BINSIZE) , int($te/$BINSIZE));
    for(my $ib=$ib1; $ib<=$ib2; $ib++) { push( @{$overlaps{$ref}{$ib}}, $rloc); }
    }
  $n_overlaps=$nr; # global
  warn"#collect_overlaps=$nr\n" if $debug;
  return \%overlaps;
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
        ##upd: tblastn need trlen==aalen w/o cdspan
        $trlen{$id}=$aw;

        if($hasspan or $testspan) { 
          if(@ac>=3 and $ac[2]=~/^\d/) { 
	          my($gp,$aq,$tw,$csp,$cspx)=@ac; # csp == span OR Code/Noncode col ..
		  $csp||=""; $cspx||="";
                  $tw= $aw unless($tw); #? or 3*aw if aa > cds size?
	          my $isutrorf=(m/utrorf/)?1:0; # key in $oid may be missing
	          my $cspan= ($csp=~/^\d/) ? $csp : ($cspx=~/^\d/) ? $cspx : 0;
	          if($testspan) {
	            #x if(($csp=~/^\d/ or $cspx=~/^\d/)) { $hasspan=1; $testspan=0; }
	            if($cspan) { $hasspan=1; $testspan=0; }
	            else { if(++$testspan>9) { $hasspan=0; $testspan=0; } }
	          }
            $aaqual{$id}=$aq; $trlen{$id}=$tw;
	          #?? still buggy? $csp=$cspx=""  if($isutrorf); # data bug: bad cds-offset for utrorfs, fixme
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


__END__

