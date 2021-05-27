#!/usr/bin/env perl
# genefindcds2.pl

=item about

  genefindcds.pl -genes trassembly_genes.gff -dnasequence genome.fasta [ -introns intron_good.gff ]
    > trassmbly_corrected_with_proteins.gff

  finds best orf, adds proteins and CDS to transcript assembly GFF
  when valid introns.gff are provided, used to check/correct for retained introns,
    and reversed-intron errors from assembly join of 2 genes.
    
  ** ASSUMES gene model is mRNA > exon (optionally CDS)
  ** ASSUMES input gff is ordered by gene records (mRNA/exon,CDS all together per ID)

=item genefindcds2 variant 2015.03

  - update for better testing of mRNA transcript > genome align.gff
    .. need more careful, base-level discrepancy of mRNA cds vs genome-map cds
    .. with possible fixes for genome mistakes (gaps, snps, indels) that damage cds/protein
  
  - add split gene handling, ie make one prot/cdna from all parts,
      expect input.gff has sorted together split parts?,
      expect mRNA Split=[123]/3 attribute to id part ordering.
      
  - add more cdna/cds sequence testing opts, input cdna or cds seq,
      output opts for cdna, cds and aa as file (instead/also in mrna.gff)
      
  - remove, move out intron test parts, for now (code clutter)

=item UPDATE 170724, fcds perl debugging/changes

  -chrcds is desired, default now, gives CDS offset,exons that best match to chr, not to mrna cdsoff
  -noshortinfix may help with better chrcds, but is foul of ncbi rule
  -completecds extends mrna align to complete a protein on chr (200/6000), make default? 
  -ratiocdnabest ? probably should make 1.25 default (now 1.50)
  -nostop 
  -full 0.75 : is this too low? ie cuts 25% of partial cds to make complete, can be many aminos
      .. replace with nAmino cutoff? e.g. partial if have > 20 .. 30 .. 50 aa partial end
      
  set pt=dpx7b3splignx.fixg; 
  $evigene/scripts/genefindcds2.pl -chrcds -noshortinfix -completecds=2 \
     -ratiocdnabest 1.25 -debug -nostopcodon -full=0.75 \
  -dna $dpxevg/genome/gasm16ml/daphplx_gasm16ml.fa -cdna $pt.mrna -genes $pt.aan.gff \
  -outgff $pt.fcdsg.gff -outaa -outcds -outcdna > & log.fcdsg$pt

=item author
  
  don gilbert, gilbertd near indiana edu, 2011
  drawn largely from Brian Haas's PASA scripts
  part of EvidentialGene, evigene/scripts/

=cut

##.. 12.07 update to share package

use constant VERSION  =>  '20170725'; # various updates, *computed cdsoff for cdnaInIsBest (mrna off vs chr off differ)
 # '20150415'; # bugfix xtrim=1 rev had bad start/stop
 # '20150330';# v2b updates for cds-exon-valid corrections; '20150325'; # v2 updates
 # '20150125'; # short intron fixes; better split gene handling from gmap, gsplign gff 
 # '20120809'; # updates to cdna_proteins.pm, other
 # '20120805'; # revised again cdnain vs gff stranding; 20120731 had bad best-strand
 # '20120731';  # utrorf, cdnain fixups
 # '20120706'; # '20120221'; 

use FindBin; 
use lib ("$FindBin::Bin", "$FindBin::Bin/../lib/"); # find Bio::DB::Fasta in evigene/lib/Bio/... from evigene/scripts/this.pl

use strict;
use warnings;
use Getopt::Long;
use cdna_proteins;
# use Bio::DB::Fasta; # now as get_dna() require  Bio::DB::Fasta, so can use this w/o bioperl  

use constant { kINTRON2SPLICE_OVER=>1, kINTRON2SPLICE_QUERY=>2, 
               kINTRONERROR_OVER=> -1, kINTRONERROR_INSIDE => -3 }; # only want last 2?
use constant SPLICE   => 3; # bp for intron splice span tests
use constant SPLICEX  => SPLICE - 1; # for exon

my $debug= 0; #pm# DEBUG

my $USE_CDSEXONS = 0;
my $USEGOODLEN=1;
my $DO_INFIX= 0; # which default?
my $DO_CDSFIX= 0; #ov2b 
my $NODIFCAN=0;
my $REANNOTATE=0;
my $allowalts= 0; # NOT 1 default? : this kills most changes, good ones..
my $pMAXSCORE = 0.05; # 1/20; 1/50?; dont consider intron cuts for this inscore < pMAX * maxvalidscore
my $DO_CDSCOMPLETE= 0; # 201511 upd to DO_CDSFIX
my $NOSHORTINFIX=0; # 201707

## cdna_proteins pack
#pm# my $AA_cdna_GT_genome= 1.50;  # default for prot(cdna) > prot(genome) test; change to 1.25 default?
#pm# my $MINAA= 30;  # used 
#pm# my $MINEXON= 60; # for cut
#pm# my $MINGOOD= 0.75; # filter out prots w/ fewer good aminos
# FIXME: adjust when to take partial vs complete, eg partial5 often is a few aa longer than complete M
#pm# my $ORF_FULLvPART = 0.85;
#pm# my $NoStopCodon=0;

my $MININTRON= 10; # NCBI sez  SEQ_FEAT.ShortIntron and SEQ_FEAT.AbuttingIntervals, ie exons overlap
  # SEQ_FEAT.ShortIntron  Introns should be at least 10 nt long
my $MINEXON= 20; # SEQ_FEAT.ShortExon  Internal coding region exon is too short .. but what is "too"?

my $BINSIZE   = 5000 ; 
my ($overlaps,$passtypes,$dnasequence,$cdnaseq,@input,$intron2splice,
    $itype,$action,$actid,$typeover,$ok,$mark,$nin);
my $mrnatypes='mRNA';
my $genetypes='gene';
my $exontypes='exon|CDS';
my $introntypes='intron'; # option?

my $overlaplist= {};
my $fasta_db= undef;
my $stranded=1; # only this opt?
my $samecds= 0; # cdna_proteins:KEEPSAMECDS; prefer keep same CDS exons but can extend/shorten protein bounds
my $debugin= undef;
my %cdnaseq; my %cdnahead;

=item updates/work

 add -output option to file updated genes
 x.add -cdnaseq option to pick best prot from asmrna transcript, compare to gff-cds
 FIXME: utrorf: works, but now only for same strand as bestorf; need to keep fwd/rev orfs, and retest against final best.
 ... add option here? to split transcript/gff for strong utrorf cases. also deal with 3+orfs/transcript: recursive?

  fix more: gene rows : need mrna span changes; lost 1st gene row on output.
  
=cut

# use constant USE_CHRMAP_CDS => 1; #? option? 1707 test upd
my $USE_CHRMAP_CDS = 1; #? option? 1707 test upd

my($outgff,$outaa,$outcds,$outcdna)=('') x 9;

my $optok= GetOptions(
  # "introns=s", \$overlaps, # ov2:skip for now
  "genes=s", \@input,  
  "dnasequence=s", \$dnasequence,  
  "cdnaseq=s", \$cdnaseq,  # input, alternate cdsseq
  "outgff:s", \$outgff, "outcdna:s", \$outcdna, "outcds:s", \$outcds, "outaa|outprot:s", \$outaa,

  "exontypes=s", \$exontypes, 
  #ov2 unused: "action=s", \$action, 
  "allowalts!", \$allowalts, 
  "passtypes=s", \$passtypes,  #??
  "minaa|minprot=i", \$MINAA,  
  "minintron=i", \$MININTRON,  
  "fullorf:s", \$ORF_FULLvPART,  
  "ratiocdnabest=s", \$AA_cdna_GT_genome, 
  "samecds!", \$samecds, 
  "nostopcodon", \$NoStopCodon, #was \$nostopcodon, 
  "goodlen!", \$USEGOODLEN, "goodmin=s", \$MINGOOD,  
  "Selenocysteine|selc!", \$USESelenocysteine,  
  "nodifcancel!", \$NODIFCAN, 
  # "fixintronerrors!", \$DO_INFIX,  # v2:skip for now
  "fixcds!", \$DO_CDSFIX, # ov2b update  
  "chrcds!", \$USE_CHRMAP_CDS, # 1707 update  
  "completecds=s", \$DO_CDSCOMPLETE, # 201511 upd to DO_CDSFIX
  "reannotate!", \$REANNOTATE, 
  "CDSEXONS!", \$USE_CDSEXONS, 
  "NOSHORTINFIX!", \$NOSHORTINFIX, 
  "debug:i", \$debugin, 
  );

die "usage:
  genefindcds.pl -genes trassembly_genes.gff -dnasequence genome.fasta [-cdna cdna.fa ]
    [ -outgff x.gff -outaa x.aa -outcds x.cds -outcdna x.cdna ]
    > trassmbly_corrected_with_proteins.gff
" unless($optok and ((@input and $dnasequence and -f $dnasequence )));  #  or $cdnaseq
#  -introns intron_good.gff 

if(defined $debugin) { $debug=($debugin>0)?$debugin:1; }

## FIXME? add option to filter out huge-span asmrna with poor qual : no intron evidence of long intron; poor prot
## protein=MLSHQLLEDSTMMQMKHGLRQGRENICQGSRLLLIGNVLVDNXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX..
$MINGOOD= $MINGOOD/100.0 if($MINGOOD > 1); 

##was $ORF_FULLvPART= $ORF_FULLvPART/100.0 if($ORF_FULLvPART > 1); 
 # ^^ -full=0 means only full orfs?; -full=1 means never?? change to -full=1 means full only; -full=0 longest part always
$ORF_FULLvPART=1 unless($ORF_FULLvPART =~ /\d/);
$ORF_FULLvPART= $ORF_FULLvPART/100.0 if($ORF_FULLvPART >= 1); 
$ORF_FULLvPART=1 if($ORF_FULLvPART == 0); # means always longest

$InnerStopToX= -1; #  unless(what?); cdna_proteins flag, turn off instops
useSelenocysteine() if($USESelenocysteine);
my $hasintrons=0;

sub MAINstub {}

if($cdnaseq) {
  my $ovh; 
  if($cdnaseq =~ /.gz$/) { $ok= open(OVR,"gunzip -c $cdnaseq |");  $ovh= *OVR; }
  elsif($cdnaseq =~ /stdin|^\-$/) { $ok=1; $ovh= *STDIN; }
  else { $ok= open(OVR,$cdnaseq); $ovh= *OVR; }
  die "bad -cdnaseq=$cdnaseq" unless($ok);
  my $id="none";
  while(<$ovh>) { if(/^>(\S+)/) { $id=$1; $cdnaseq{$id}=""; 
    if(m/ +(\S.+)$/){ my $h=$1; $h =~ s/\s*(len|cf|nt)=\S+//g; $cdnahead{$id}=$h; }
    } elsif(/\w/) { chomp; $cdnaseq{$id}.= uc($_); } } 
  close($ovh);

  #ov2 cdna_bestorf() unless(@input); #ov2: drop this, dont need here
}

my($houtgff,$houtaa,$houtcds,$houtcdna)= (undef) x 9; # global output handles

$houtgff= *STDOUT;
my $otemp= $outgff || $input[0]; $otemp=~s,\.\w+$,,;  $otemp.="_fcds" unless($outgff);
  $outgff="$otemp.gff" if(defined $outgff and not $outgff);
  $outaa="$otemp.aa" if(defined $outaa and not $outaa);
  $outcds="$otemp.cds" if(defined $outcds and not $outcds);
  $outcdna="$otemp.cdna" if(defined $outcdna and not $outcdna);

$ok=1;
$ok= open($houtgff,'>',$outgff) if($ok and $outgff);
$ok= open($houtaa,'>',$outaa) if($ok and $outaa);
$ok= open($houtcds,'>',$outcds) if($ok and $outcds);
$ok= open($houtcdna,'>',$outcdna) if($ok and $outcdna);

foreach my $input (@input) {
  my $inh= *STDIN;
  $ok = ($input =~ /.gz$/) ? open($inh,"gunzip -c $input |") 
        : ($input =~ /^(stdin|-)/) ? $inh= *STDIN
        : open($inh,$input);
  die "bad -input=$input" unless($ok);
  
  my ($nchanged,$ngene)= filter_gff($inh);
  warn "#findcds changed=$nchanged, ngene=$ngene\n" if $debug;
}

# close($houtgff); close(..);

# end MAINstub
#..................

# sub bestorf_test{ } # in cdna_proteins.pm

# 201207 ** REPLACED by cdna_bestorf.pl and cdna_proteins.pm




sub filter_gff
{
  my($inh)= @_;
  my ($ng,$nr,$nchange)= (0,0,0);
  my $printpass=1;
  # $printpass=0 if($actid == ACT_DROP or $actid == ACT_KEEP);
  my @generec=(); my @otherft=(); my @geneft;
  
  ## add own header ?? version?
  my $version=VERSION;
  print $houtgff <<"EOGFF";
##gff-version 3
# appl: genefindcds
# vers: $version

EOGFF
  
  my ($lasplit,$lgid)=(0,0);
  while(<$inh>){
    unless(/^\w/){ next if(/^(##gff-ver|#n |$)/);  print $houtgff $_ and next; }
    my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr)= split"\t";
    $nr++; chomp($tattr);
    if($passtypes and "$typ.$src" !~ m/$passtypes/) { print $houtgff $_ if $printpass; next; } 

# FIXME: ref == NOPATH .. need cdna-only results, modify gff to use cdna-seq as scaffold?
    
    my($gid,$pid); 
    if($tattr =~ m/\bID=([^;]+)/) {  $gid=$1; }
    if($tattr =~ m/\bParent=([^;]+)/) {  $pid=$1; 
      $gid=$pid unless($gid and ($typ =~ /^($mrnatypes)$/)); 
      }
    unless(defined $gid) { $gid = "N".$ng; }
    my $gidfix= $gid; $gidfix =~ s/_C(\d+)$//;
    
    my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid]; 

    if($typ =~ /^($mrnatypes)$/) {  
      # allow gene and mRNA types .. and/or keep other types in generec
      
      #ov2: Split=1/3 .. gene fix here, need to collect all parts (mRNA:1,2,.. exons..) for testgene()
      my $issplit=($tattr =~ m/;Split=([^;\s]+)/)? $1 : ($gid =~ /_C(\d+)$/)? $1 : 0;
      if($issplit and $gidfix eq $lgid) {
        push @generec, $rloc; # add mRNA part2..     
      
      } else {
        $nchange += testgene(\@generec, \@otherft) if(@generec);
        @generec = ($rloc); $ng++; # maybe best store as array of [gffcols], sorted by type-loc
        @otherft=(); # ** drops prior gene ft
        if(@geneft) { unshift @generec, @geneft; @geneft=(); } #?? put into generec ? will cause problems
      }
      
      $lgid= $gidfix; $lasplit= $issplit;
      
    } elsif($typ =~ /^($exontypes)$/) {
      push @generec, $rloc;         
    } elsif($typ =~ /^($genetypes)$/) {
      @geneft =( $rloc ); # keep only latest
      
    } elsif($tb>0 and $te>0) {
      push @otherft, $rloc;         
    }
      
  }
  
  # FIXME lost 1st gene row before mrna; FIXME change gene span w/ mrna span changes
  
  $nchange += testgene(\@generec, \@otherft) if(@generec);
  
  return ($nchange,$ng);
}


sub putseq {
  my($outh,$id,$seq,$attr)= @_;
  return -1 unless($id and $seq and $outh);
  # $seq =~ s/\*$//; #?? always remove stop codon for outfile? sometimes? best use -nostop opt
  if($NoStopCodon) { $seq =~ s/\*$//; }
  my $slen=length($seq); # none have \n or \s ?
  $seq =~ s/(.{60})/$1\n/g; $seq.="\n" unless($seq=~m/\n$/);
  $attr||=""; #  aalen=$aalen,$pcds%,$compl; clen=$clen; strand=$crev; offs=$prostart5-$proend3; $cdnah
  $attr.=" len=$slen;" unless($attr =~ m/len=$slen/); # qlen? clen?
  print $outh ">$id $attr\n", $seq;
}

sub putgene {
  my ($generec,$otherft,$flags)= @_; # add?? opt: $seqs = [seqrec for aa,cds,cdna ]
  ##  my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid]; 
  my $cc= ($flags and $flags =~ /skip=/)? "#x." : "";
  $otherft ||= [];

  my $checkGeneSpan= ($flags and $flags =~ /mrnachanged/)? 1:0;
  if($checkGeneSpan) { $flags =~ s/mrnachanged=.//;
    my($gene) = grep{ $_->[2] eq "gene" } @$generec;
    my($mrna) = grep{ $_->[2] eq "mRNA" } @$generec; # FIXME Split gene @mrna
    my $didfix= ($gene and $mrna) ? mrna_fixspan($gene,[$mrna]) : 0; 
  }
  
  foreach my $ft (@$generec, @$otherft) { 
    if(ref $ft) { my @v= @$ft; 
      $v[8]=~s/$/;$flags/ if($flags and $v[2] eq "mRNA");
      print $houtgff $cc.join("\t",@v[0..8])."\n" if(@v>4); 
      }
    }
  print $houtgff "\n"; #?
}



# in cdna_proteins.pm
# sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
# sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }

sub splitPart {
  my($xattr)=@_;  # gff attribute col 9
  my($xid) = $xattr =~ m/(?:ID|Parent)=([^;\s]+)/;  
  ## which split key has precedence? Split=\d or ID_C\d ?
  my($spl)= ($xid =~ /_C(\d+)$/)? $1 : ($xattr =~ m/;Split=(\d+)/)? $1 : 0;
  return($spl,$xid);
}

sub _sort_over { # @[b,e] min-b, max-e
  return ($a->[0] <=> $b->[0]) || ($b->[1] <=> $a->[1]);
}

sub _sortgene  
{
  # ($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid)
  #FIXME? for split genes, Split=1,2,3/3 order by split-part? all mRNA at top or not?
  my($ta,$tb)= map{ (m/gene/)?1:(m/mRNA/)?2:(m/CDS|exon/)?3:4; } ($a->[2],$b->[2]);
  return ($ta <=> $tb)
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

sub _sortSplitgene {
  # ($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid)
  #FIXME? for split genes, Split=1,2,3/3 order by split-part? all mRNA at top or not?
  # need to apply split sort to each part: mrna, exon,cds .. better way?
  # my($aid,$bid)= map{ m/(?:Parent|ID)=([^;\s]+)/; $1; } ($a->[8],$b->[8]);
  # my($spla)=($a->[8] =~ m/;Split=(\d+)/)? $1 : ($aid =~ /_C(\d+)$/)? $1 : 0;
  # my($splb)=($b->[8] =~ m/;Split=(\d+)/)? $1 : ($bid =~ /_C(\d+)$/)? $1 : 0;
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

sub _sortloc { # _sortgene
  #  my($ref,$start,$stop,$strand)= @{$ft}[0,3,4,6];
  return ($a->[0] cmp $b->[0]) # ref : should all be same for gene record
      || ($a->[3] <=> $b->[3]) # begin
      || ($a->[4] <=> $b->[4]) # end
      ;
}



# sub getBestProt ## in cdna_proteins.pm; revised
  
# getBestProt fix this to return longest full prot, and longest partial (if longer)
# .. test which is best.
# FIX2: add test intron overlaps : retained = intron inside exon; err = intron rev at splice/inside
#   my $longorf= $longest_orf_finder->get_longest_orf($cdna); # == hash
# FIXME3: option to mark/return 2ndary orf(s) in aberrant long-utr transcripts, that dont overlap 1st orf

sub getBestProtOfOrfs
{
  my($longorf,$longfull,$orfs,
     $ptype, $cdna, $exongff, $oldStart_b, $oldStart_e)= @_;
  my($orfprot,$prostart5,$proend3)=("",0,0,"");
  my ($utrorf,$utrosize)=(undef,0);
   
  if(ref($longorf)) {
    # $orfprot= $longest_orf_finder->get_peptide_sequence();
    # ($prostart5,$proend3)= $longest_orf_finder->get_end5_end3();

    my $lookmore= ($ptype =~ /long/)?0:1;
    if($samecds and $oldStart_b > 0) { # may not be right yet.
      my($samestartorf);
      if($ORF_FULLvPART <= 0.8) {
      ## bad test should be oldStart_b <= orf.start <= oldStart_e, NOT orf.end
      ##($samestartorf) = grep { $_->{complete} == 3 and $_->{start} >= $oldStart_b and  $_->{stop} <= $oldStart_e } @$orfs;
      ($samestartorf) = grep { $_->{complete} == 3 and $_->{start} >= $oldStart_b and  $_->{start} <= $oldStart_e } @$orfs;
      } else {
      ($samestartorf) = grep { $_->{start} >= $oldStart_b and  $_->{start} <= $oldStart_e } @$orfs;
      }
      if(ref $samestartorf) { $longorf= $samestartorf; $lookmore=0; } #NOT# else { return (undef); } # not found == no change, here
    } 
    
    if($lookmore and $longorf->{complete} < 3 and ref($longfull) ) {
      my $keylen=($USEGOODLEN)?"goodlen":"length";
      my $lsize= $longorf->{$keylen}; # was {length}
      my $fsize= $longfull->{$keylen};
      
      # FIXME: adjust when to take partial vs complete, eg partial5 often is a few aa longer than complete M
      if($fsize >= $ORF_FULLvPART * $lsize) {
        if(ref($exongff)) {
        my ($cdslong, $attrL, $pcodeL, $maxutrL)= getCDSgff2( $exongff, $longorf);
        my ($cdsfull, $attrF, $pcodeF, $maxutrF)= getCDSgff2( $exongff, $longfull);
        $longorf= $longfull if( $maxutrF < 3 or ($pcodeF >= $ORF_FULLvPART * $pcodeL));
        } else {
	      $longorf= $longfull;
	      }
      }
    }
    
    ($orfprot,$prostart5,$proend3)= orfParts($longorf);
    ## ($orfprot,$prostart5,$proend3)= ($longorf->{protein},$longorf->{start}, $longorf->{stop});

    ($utrorf,$utrosize)= getUtrOrf($longorf, $orfs, length($cdna));  
  }

  return($orfprot, $prostart5, $proend3, $longorf, $utrorf);
}


# sub getUtrOrf # in cdna_proteins.pm

# in cdna_proteins.pm; revised as getCDSgff2 << use that one?
# sub getCDSgff # _OLD
# {
#   my($exons,$orfprot,$prostart5,$proend3,$cdnalen) = @_;
#   # FIXME: need trlen= length(cdnain) for -cdna, and/or use gmap qlen= tag
#   # FIXMEs: exon Split= needs copy to CDS .. losing it, here???
#   $cdnalen ||= 0;
#     
#   my ($cds5,$cds3,$cds3last)=(0,0,0);
#   my @cds= ();
#   my @utr= ();
#   ## for phase; need reverse @exons
#   my ($cdna1,$inc5,$inc3,$nt_length, $nu5, $nu3)= (0) x 10;
#   
#   $cdna1= 0; # was 1; # offby1 at end?
#   $nt_length= 0; # $prostart5 % 3; #??
#   
#   # ** FIXME 2011Dec : stopcodon split intron >> CDS ends w/o final 1,2 bases ** WRONG
#   # .. must make next exon(if exists) part of CDS stop
#   
#   #?? rev bug here? got bad CDS for good prot/cdna after completeCDSb, revgene
#   my $addat="";
#   if($exons->[0]->[6] eq "-") {
#     my @xbeg= @{$exons->[0]}[3,4,6]; # b,e,o
#     my @xend= @{$exons->[-1]}[3,4,6];
#     if($xbeg[0] < $xend[0]) {
#       $addat.=",badrev"; 
#       my @xrev= reverse @$exons; $exons= \@xrev; 
#       }
#   }
#   
#   foreach my $exon (@$exons) {
#     my ($ref,$src,$xtyp,$xend5, $xend3,$xv,$xo,$xph,$xattr,$gid)= @{$exon};
#     
#     my $cdsattr=$xattr; 
#     # my $KEEPAN='Split|err|gaps'; #? |Target|trg |gapfill|gapfix ?
#     my $DROPAN='Target|trg|gapfill|gapfix|splice';
#     $cdsattr=~s/;($DROPAN)=[^;\n]+//g; #** Target|trg has spaces **
#     $cdsattr=~s/Parent=[^;\s]+[;]?//; #? or leave on should be same as $gid
#     $cdsattr="" unless($cdsattr=~/\w+=/);
#     
#     ($xend5, $xend3)= ($xend3,$xend5) if($xo eq "-"); #patch rev?
#     my $xd= abs($xend3 - $xend5); # ?? +1 for width
#     
#     $cdna1++; # add 1 here, not end loop 
#     my $cdna2= $cdna1 + $xd; # ?? +1 for width
#     # ** offby1 here ?? YES, need <=, >= to get full CDS stop,start split by intron; see below cdna1=cdna2+1
#     #OLD.if($cdna1 < $proend3 and $cdna2 > $prostart5) 
#     if($cdna1 <= $proend3 and $cdna2 >= $prostart5) 
#     { # overlap
#                 
#       my $d5= ($cdna1 >= $prostart5) ? 0 : $prostart5 - $cdna1; # pos
#       my $c5= ($xend5 > $xend3) ? $xend5 - $d5 : $xend5 + $d5;    				  
#       
#       my $d3= ($cdna2 <= $proend3) ? 0 : $proend3 - $cdna2; # neg
#       my $c3= ($xend5 > $xend3) ? $xend3 - $d3 : $xend3 + $d3;
#   
#       my $elength = 1 + abs($c3 - $c5);
#       $nt_length  += $elength;
#       $inc3        = $nt_length % 3;
#       $inc5        = ($elength - $inc3) % 3; # only care about this one
#       # $frame       = ($c5 + $inc5) % 3;
#       if ($inc5 == -1) { $inc5 = 2; }
#       
#       my $phase= $inc5; # is this right?
#       
#       my($cb,$ce)= ($c5 > $c3) ? ($c3,$c5): ($c5,$c3); #? rev patch
#       # exon Split=, other xattr need copy to CDS HERE ***
#       my $rloc= [$ref,$src,"CDS",$cb,$ce,".",$xo,$phase,"Parent=$gid;$cdsattr",$gid]; 
#       push @cds, $rloc;
#      
#       if($cdna1 <= $prostart5) { $cds5=$c5; }
#       if($cdna2 >= $proend3) { $cds3=$c3; } else { $cds3last= $c3; } # cdna2 < proend3 here
#  
#     } elsif(1) { # $addutr .. not used?
#       my $d5= ($cdna1 >= $proend3) ? 0 : $proend3 - $cdna1; # pos
#       my $u5= ($xend5 > $xend3) ? $xend5 - $d5 : $xend5 + $d5;    				  
#       my $d3= ($cdna2 <= $prostart5) ? 0 : $prostart5 - $cdna2; # neg
#       my $u3= ($xend5 > $xend3) ? $xend3 - $d3 : $xend3 + $d3;
# 
#       my($ub,$ue)= ($u5 > $u3) ? ($u3,$u5): ($u5,$u3); 
#       my $up= ($cdna1 < $prostart5) ? "five" : ($cdna2 > $proend3) ? "three" : "odd";
#       if($cdna1 < $prostart5) { $nu5++; } elsif($cdna2 > $proend3) { $nu3++; }
#       my $rloc= [$ref,$src, $up."_prime_utr",$ub,$ue,".",$xo,0,"Parent=$gid;$cdsattr",$gid]; 
#       push @utr, $rloc;
#     }
#     
#     ##$cdna1= $cdna2+1;  # is this off-by-1 now? yes, dont +1 here, do above
#     $cdna1= $cdna2;  
#   }       
#   $cds3||=$cds3last; # partial3
#   
#   my $trlen= ($cdnalen>$cdna1) ? $cdnalen : $cdna1; # if($cdnalen> $cdna1 or > 0) ??
#   my $aalen=length($orfprot); 
#   $aalen-- if(substr($orfprot,-1) eq '*');
#   my $clen= $aalen * 3; # can be off by -1,-2 here. 
#   my $ap=int(100 * $clen/$trlen);
# 
# ## .. add this prot qual test, also test orig prot.. esp for augustus X inner stops
# ##        if($istop < 1 and $id =~ /AUG/) { $istop= index($aa,'X'); }
# ##        if($istop > 0 and $istop < $al-1) { $astat="partialinner"; }
#   
#   my $mattr="cxlen=$clen/$trlen;aalen=$aalen,$ap%";
#   if($orfprot) {
#     my $p5= (substr($orfprot,0,1) eq 'M')?0:1;
#     my $p3= (substr($orfprot,-1,1) eq '*')?0:1;# FIXME: -nostop/NoStopCodon .. need aaqual
#     my $prostat = ($p5 and $p3) ? "partial": ($p5)? "partial5" : ($p3)? "partial3" :"complete";
#     $mattr.= ",$prostat;protein=$orfprot";
#   }
#   
#   $mattr.= ";cdsoff=$prostart5-$proend3"; #? as per ;utroff=$ustart-$uend
#   $mattr.= ";cdsspan=$cds5-$cds3$addat"; #?add for other uses? rev if needed? or not?
#   $mattr.= ";utrx=$nu5,$nu3" if($nu5 > 2 or $nu3 > 2); # ;utrx=$u5,$u3
#   ## mattr keys: cxlen,aalen,protein,utrx
#   
#   # FIXME: resort @cds by loc, not reversed. : let caller do
#   # @cds = sort _sortgene @cds;
# 
#   ## return also: $clen, $trlen or $utrlen or $ap, $nu5+$nu3, 
#   ## ($cdslong, $attrL, $pcodeL, $maxutrL)
#   return (\@cds, $mattr, $ap, _max($nu5,$nu3), \@utr); 
# }


    #?? move TRIM2 out to sub ($exondna,$exonft,$didtrim)= trimcdnaexon($jexon,$exondna,$exonft) ?
use constant CDNA_TRIMNNN => 1;  # need user opt to turn on/off; debug  
use constant { TRIM2 => 1, kEXONNOTRIM => 0, kEXONTRIM => 1, kEXONDROP => 2, };

=item trimcdnaexon trim exon NNN, cdna ends only

    FIXME for gaps, chomp off end gaps NNN of getcdna() ? or remove from orfprot
    .. need to change exongff also for cdna trimnnn .. messy
    ** THIS MAY BE BUGGY .. problematic at least **
    .. maybe move it to exon loop above, trim end exons, then into next ones if needed?
    .. *should not* need to trim 2nd+ exons otherwise they would not be mapped as exons onto all gaps
    .. if(j==0) checktrim 0>up;  if(j==nx1) checktrim end>down
    
   see testgene() does this ONLY w/ cdnain; do for all?
   
=cut

sub trimcdnaexon {
  my($jexon,$exondna,$exonft,$rev,$nx1,$whichj)= @_; 
  my $didtrim=kEXONNOTRIM;  $whichj||="";
  my($start,$stop)= @{$exonft}[3,4]; 
  my $j=$jexon;
  if(($j==0 or $j==$nx1 or $whichj=~/first|last/) and index($exondna,'N') >= 0) { # $trimNNN and 
    ## dangit .. $rev == reversed exon order
    my($first,$last);
    if($whichj) {
      $first=(($rev and $whichj eq "last") or (not $rev and $whichj eq "first")) ?1:0;
      $last =(($rev and $whichj eq "first") or (not $rev and $whichj eq "last")) ?1:0;
    } else {
      $first=(($rev and $j==$nx1) or (not $rev and $j==0)) ?1:0;
      $last =(($rev and $j==0) or (not $rev and $j==$nx1)) ?1:0;
    }
    
    if($first) {  ## ($j==0)  # first
      my $texon= $exondna;
      my $xl= length($texon); my $xb=-1; 
      my $xi= ($texon=~/^N/) ? 0 : index($texon,"NN"); 
      while($xi>$xb and $xi<$xb+21) {
        $xi++ while(substr($texon,$xi,1) eq 'N'); # do again?
        my $xn= index($texon,"NN",$xi); 
        if($xn > $xi and $xn < $xi+21) { $xb=$xi; $xi=$xn;  }
        else { $xb=$xi-1; $xi=-1; } ## xi>=0 here
      }
      if($xb>=0) {
        # FIXME kEXONDROP implies change first/last to next in line.. (or last for rev & first)
        if($xb >= $xl) { $texon=""; $start=1+$stop; $didtrim=kEXONDROP; } # drop
        else { 
          $texon= substr($texon,$xb+1); 
          $start += $xb+1;  $didtrim=kEXONTRIM;
          #* FIXME cloned, changed exons not seen
          my @xclone= @$exonft; $exonft= \@xclone;  
          #BUG**# if($rev) { $exonft->[4] -= $xb+1; } else { $exonft->[3]= $start;  } 
          ##^^BUG first == last for rev
          $exonft->[3]= $start; # for both rev/fwd
          }
        $exondna=$texon; # $xchange++;
      } else { 
        # push @exok, $exonft; 
      }
      
    } elsif($last) { ## ($j==$nx1)  # last; 1exon do both ??
      my $texon= reverse($exondna); # reuse above code, reversed..
      my $xl= length($texon); my $xb=-1; 
      my $xi= ($texon=~/^N/) ? 0 : index($texon,"NN"); 
      while($xi>$xb and $xi<$xb+21) {
        $xi++ while(substr($texon,$xi,1) eq 'N'); # do again?
        my $xn= index($texon,"NN",$xi); 
        if($xn > $xi and $xn < $xi+21) { $xb=$xi; $xi=$xn;  }
        else { $xb=$xi-1; $xi=-1; } ## xi>=0 here
      }
      if($xb>=0) {
        # UNreverse here..
        if($xb >= $xl) { $texon=""; $stop=$start-1; $didtrim=kEXONDROP;} # drop
        else { 
          $texon= substr($texon,$xb+1); 
          $texon= reverse($texon);
          $stop -= $xb+1;  $didtrim=kEXONTRIM;
          my @xclone= @$exonft; $exonft= \@xclone;  
          #BUG#if($rev) { $exonft->[3] += $xb+1; } else { $exonft->[4]= $stop;   } 
          ##^^BUG last == firs for rev
          $exonft->[4]= $stop; # for both rev/fwd
          }
        $exondna=$texon; # $xchange++;
      } else { 
        # push @exok, $exonft; 
      }
    }
  }
  return($didtrim,$exondna,$exonft,$start,$stop); #upd: start,stop
}

use constant DOREVCOMP => 1;

sub getexondna {
  my($exonft,$dorevcomp)= @_; 
  my($ref,$start,$stop,$strand,$phase,$xgid)= @{$exonft}[0,3,4,6,7,9];
  my $rev=($strand eq "-")?1:0;
  my $exondna  = get_dna( $dnasequence, $ref, $start, $stop);
  $exondna ||="";
  if($debug and not $exondna) {
    ## exonft == $ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid
    warn "#xdna=0 gid=$xgid loc=$ref:$start-$stop\n"; # loc=NOPATH:1-69
  }
  $exondna = uc($exondna); # always?
  my $xlen = length($exondna); # 1+$stop-$start; # add even if no dna
  $exondna = revcomp($exondna) if($rev and $dorevcomp);  # leave to caller? may need unrev
  return ($exondna,$rev,$xlen);
}

# sub getcdna2 { return getcdna2Test(@_); } #or getcdna 

sub getcdna2 { #  getcdna2Test  only for trimNNN now
  my($exons, $asexons, $expand, $expend, $trimNNN)= @_;  # add flag: trimNNNends. drop expand,expend ..
  my $xchange=0;
  my $nexons=  @$exons;
  my $nx1= $nexons - 1;
  # @$exons; _sortgene, always start<stop ?? maybe bug here in sort
  # NO:  my @sexon= sort _sortloc @$exons;
  
  my @j= (0..$nx1); 
  my $je= pop(@j); unshift @j, $je;
  my @exok= (0) x $nexons;
  my @exseq= ("") x $nexons;
  my $nx0= 0;
  
  for my $j (@j) {
    my $ft= $exons->[$j]; # $sexon[$j]; ??
    my($exondna,$rev,$xlen)= getexondna($ft,0); # DOES NOT revcomp(dna) if rev; DOREVCOMP does it
    $exok[$j]= $ft; $exseq[$j]=$exondna;  
    if($trimNNN) {
      my $whichj= ($j == $nx1) ? "last" : ($j == $nx0) ? "first" : "middle";
      # not bug: trimcdnaexon expects rev/fwd always w/ fwd(dna)
      my($didtrim,$texondna,$tft,$tstart,$tstop)= trimcdnaexon($j,$exondna,$ft,$rev,$nx1,$whichj);
      if($didtrim) {
        $exondna=$texondna; $ft=$tft; $xchange++;
        if($didtrim == kEXONDROP) { 
          $exok[$j]=0; $exseq[$j]="";
          if($j==$nx0) { $nx0++; }  
          elsif($j==$nx1) { }
          $nx1--;
          # if($j==$nx1) recheck exons->[j-1]; ?? or do j=nx1 out of order? and remember if dropped
          # if($j==0) update j,nx1 for next call of trimcdnaexon()
          # for j==1, trimcdnaexon($j-1,$exondna,$ft,$rev,$nx1-1);
        } else {      
          $exok[$j]= $tft; $exseq[$j]= $texondna;  
        }
      }
    }
  }

  my(@exonft,@exonseq);
  my $cdna= ""; my $cdnalen=0; 
  for my $j (0..$nexons-1) {
    my $exondna= $exseq[$j] or next;
    my $exft= $exok[$j];
    my $rev= ($exft->[6] eq "-")?1:0;
    #unless DOREVCOMP# 
    $exondna = revcomp($exondna) if($rev);  # $exondna = reverse $exondna;  $exondna =~ tr/ACGTacgt/TGCAtgca/;
    $cdna .= $exondna;  
    push @exonseq, $exondna;
    push @exonft, $exft;
  }
  
  #?? if($trimNNN and $xchange) { $exons= \@exonft; } # TRIM2: other parts updated: $cdna,$cdnalen,\@asexons
  $cdnalen= length($cdna);
  return($cdna,$cdnalen,\@exonseq,\@exonft,$xchange);
}

  
sub getcdna {
  my($exons, $asexons, $expand, $expend, $trimNNN)= @_;  # add flag: trimNNNends
  my $cdna= ""; my $cdnalen=0; 
  my @asexons=(); my @exok=(); my $xchange=0;
  my $lstop= 1;   
  my $doexp=(defined $expand && $expand>0)?1:0;
  if($doexp) { $expend=$expand unless($expend); }
  my $nx1= @$exons - 1;
  
  ## NOTE: rev == exons are reversed order, top = 1st in array
  foreach my $j (0..$nx1) {   # @$exons; how are these sorted? _sortgene, always start<stop
    my $ft= $exons->[$j];
    my($ref,$start,$stop,$strand,$phase,$fgid)= @{$ft}[0,3,4,6,7,9];
    my $rev=($strand eq "-")?1:0;
    # if($samecds and $phase>0 and $cdnalen==0) { if($rev) { $stop-=$phase; } else { $start+=$phase; } }
    # ^^ not here, use phase in get_orfs ..
    # $cdnalen += 1+$stop-$start; # add even if no dna
    my($xstart,$xstop)= ($start, $stop);
    if($doexp) {
       my $nstart= ($j<$nx1)? $exons->[$j+1]->[3] : 999999999;
       my $xp=($j==0)? $expend : $expand; $xstart= _max(1, _max($lstop,$xstart - $xp));
       $xp=($j==$nx1)? $expend : $expand; $xstop=  _min($nstart, $xstop + $xp);
    }
    my $exondna  = get_dna( $dnasequence, $ref, $xstart, $xstop);
    $exondna ||="";
    if($debug and not $exondna) {
      ## exonft == $ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid
      warn "#cdnaseq=0 gid=$fgid loc=$ref:$xstart-$xstop\n"; # loc=NOPATH:1-69
    }
    $exondna = uc($exondna); # always?

    if($trimNNN) {
        # FIXME kEXONDROP implies change first/last to next in line.. (or last for rev & first)
      my($didtrim,$texondna,$tft,$tstart,$tstop)= trimcdnaexon($j,$exondna,$ft,$rev,$nx1);
      
      if($didtrim) {
        $exondna=$texondna; $ft=$tft; 
        ($start,$stop)= ($tstart,$tstop); # @{$ft}[3,4]; # update; drop: start=stop
        push @exok, $ft unless($didtrim == kEXONDROP);
        $xchange++;

        if($didtrim == kEXONDROP) { 
          # if($j==$nx1) recheck exons->[j-1]; ?? or do j=nx1 out of order? and remember if dropped
          # if($j==0) update j,nx1 for next call of trimcdnaexon()
          # for j==1, trimcdnaexon($j-1,$exondna,$ft,$rev,$nx1-1);
        }
        
      } else {
        push @exok, $ft;
      }
    } else {
      push @exok, $ft;
    }
    
    $cdnalen += 1+$stop-$start; # add even if no dna
    $lstop= $stop;
    next unless($exondna);
    $exondna = revcomp($exondna) if($rev);  # $exondna = reverse $exondna;  $exondna =~ tr/ACGTacgt/TGCAtgca/;
    $cdna .= $exondna;  
    push @asexons, $exondna; ## if($asexons or $trimNNN);
  }
  
  $cdnalen= length($cdna) if($cdna);
  if($trimNNN and $xchange) { $exons= \@exok; } # TRIM2: other parts updated: $cdna,$cdnalen,\@asexons
  return($cdna,$cdnalen,\@asexons,$exons,$xchange);
}  

#----------------------
#    #?? move TRIM2 out to sub ($exondna,$exonft,$didtrim)= trimcdnaexon($jexon,$exondna,$exonft) ?
# if(TRIM2) {     
#     if($trimNNN and ($j==0 or $j==$nx1) and index($exondna,'N') >= 0) {
#       ## dangit .. $rev == reversed exon order
#       my $first=(($rev and $j==$nx1) or (not $rev and $j==0)) ?1:0;
#       my $last=(($rev and $j==0) or (not $rev and $j==$nx1)) ?1:0;
#       
#       if($first) {  ## ($j==0)  # first
#         my $texon= $exondna;
#         my $xl= length($texon); my $xb=-1; 
#         my $xi= ($texon=~/^N/) ? 0 : index($texon,"NN"); 
#         while($xi>$xb and $xi<$xb+21) {
#           $xi++ while(substr($texon,$xi,1) eq 'N'); # do again?
#           my $xn= index($texon,"NN",$xi); 
#           if($xn > $xi and $xn < $xi+21) { $xb=$xi; $xi=$xn;  }
#           else { $xb=$xi-1; $xi=-1; }
#         }
#         if($xb>=0) {
#           if($xb >= $xl) { $texon=""; $start=1+$stop; } # drop
#           else { 
#             $texon= substr($texon,$xb+1); $start += $xb+1; 
#             #* FIXME cloned, changed exons not seen
#             my @xclone= @$ft; $ft= \@xclone;  $ft->[3]= $start;  push @exok, $ft;
#             }
#           $exondna=$texon; $xchange++;
#         } else { push @exok, $ft; }
#         
#       } elsif($last) { ## ($j==$nx1)  # last; 1exon do both ??
#         my $texon= reverse($exondna); # reuse above code, reversed..
#         my $xl= length($texon); my $xb=-1; 
#         my $xi= ($texon=~/^N/) ? 0 : index($texon,"NN"); 
#         while($xi>$xb and $xi<$xb+21) {
#           $xi++ while(substr($texon,$xi,1) eq 'N'); # do again?
#           my $xn= index($texon,"NN",$xi); 
#           if($xn > $xi and $xn < $xi+21) { $xb=$xi; $xi=$xn;  }
#           else { $xb=$xi-1; $xi=-1; }
#         }
#         if($xb>=0) {
#           # UNreverse here..
#           if($xb >= $xl) { $texon=""; $stop=$start-1; } # drop
#           else { 
#             $texon= substr($texon,$xb+1); 
#             $texon= reverse($texon);
#             $stop -= $xb+1; # $stop = $stop - 1 - $xb 
#             my @xclone= @$ft; $ft= \@xclone;  $ft->[4]= $stop;  push @exok, $ft;
#             }
#           $exondna=$texon; $xchange++;
#         } else { push @exok, $ft; }
#       }
#     } else {
#       push @exok, $ft;
#     }
# }


# unless(TRIM2) { #ov2: drop,obsolete
#   if($trimNNN and index($cdna,'N') >= 0) { 
#     my $cdnagood= $cdna; my $clen=$cdnalen;
#     my @xn= @asexons; my @exok= @$exons;
#     
#     $cdnagood =~ s/^N+//;
#     my $xi= index($cdnagood,"NN"); 
#     if($xi>0 and $xi<20) { 
#       $xi++ while($xi<$cdnalen and substr($cdnagood,$xi,1) eq 'N');
#       $cdnagood= substr($cdnagood,$xi); 
#     }
#     my($i,$d,$clenb)=(0) x 9;
#     $clenb= length($cdnagood);
#     $d=$clen - $clenb; $i=0;
#     while($d>0) {
#       my $x= $asexons[$i]; my $ft= $exons->[$i];
#       my $l=length($x);
#       # my($ref,$start,$stop,$strand,$phase)= @{$ft}[0,3,4,6,7];
#       if($l <= $d) {
#         $d -= $l; $i++;
#         shift @exok; shift @xn;
#       } else { # $l > $d
#         $xn[$i]= substr($x,$d); $ft->[3] += $d; $d=0;
#       }
#     }
#     
#     $clen= length($cdnagood);
#     $cdnagood =~ s/N+$//; # fix exon3
#     $clenb= length($cdnagood);
#     my $xe= rindex($cdnagood,"NN");
#     if($xe > $clenb-20) {
#       $xe-- while($xe>0 and substr($cdnagood,$xe,1) eq 'N');
#       $cdnagood= substr($cdnagood,0,$xe); 
#     }
#     $clenb= length($cdnagood);
#     $d=$clen - $clenb; $i= @$exons - 1;
#     while($d>0) {
#       my $x= $asexons[$i]; my $ft= $exons->[$i];
#       my $l=length($x);
#       # my($ref,$start,$stop,$strand,$phase)= @{$ft}[0,3,4,6,7];
#       if($l <= $d) {
#         $d -= $l; $i--;
#         pop @exok; pop @xn;
#       } else { 
#         $xn[$i]= substr($x,0,$l-$d); $ft->[4] -= $d; $d=0;
#       }
#     }
#     
#     
#     $clen=length($cdnagood);
#     if($clen<$cdnalen) { 
#       $cdna=$cdnagood; $cdnalen=$clen; 
#       @asexons= @xn;
#       $exons= \@exok; $xchange=1;
#     }
#   } 
# }  ## TRIM1 not 2
  


sub stripOldAnnot {
  my ($mrna,$oldtags)= @_;
  return unless(ref $mrna and $mrna->[8]);
  $oldtags= [qw(cxlen aalen protein cdnabest cdnaorf aautrlen utroff utrprot ocds oaaln inqual intronfix xcut xdrop)]
    unless(ref $oldtags);
  my($an)= $mrna->[8]; my $oldan= $an;
  my $olds= join('|', @$oldtags);
  if( $an =~ s/;($olds)=[^;\n]+//g ) { $mrna->[8]= $an; } # or not?
  return ($an,$oldan); 
}

sub seqAnnot {
  my($mrna,$keeptags)= @_;
  return unless(ref $mrna and $mrna->[8]);
  $keeptags= [qw(aalen cdsoff cxlen clen qlen aaold oid)] unless(ref $keeptags);
  ## offs == cdsoff ?; clen != qlen != cxlen ? all of these or what?
  my($an)= $mrna->[8]; 
  my $keepan="";
  for my $t (@$keeptags) { if($an =~ m/\b($t=[^;\n]+)/) { $keepan.="$1; "; } }
  return $keepan; 
}

sub mrna_fixspan {
  my($mrna,$exons)= @_;
  my $changed=0;
  ## splitgene fix: check ref same.
  
  # FIXME lost 1st gene row before mrna; FIXME change gene span w/ mrna span changes
  ##my($mrna) = grep{ $_->[2] eq "mRNA" } @generecfix;
  my($mpart,$geneid)= splitPart($mrna->[8]);
  my($mref,$oldstart,$oldstop,$oldor)= ($mrna->[0],$mrna->[3],$mrna->[4],$mrna->[6]);
  my($newstart,$newstop,$nrev,$nfwd, $nxok)= (0) x 9;
  
  foreach my $ex (@$exons) {
    my($xr,$xb,$xe,$xo,$xat)= @{$ex}[0,3,4,6,8];
    next if($xr ne $mref);
    my($xpart)= splitPart($xat);
    next if($xpart ne $mpart); # Split fix
    ## bug split part mrna get all parts span on same chr .. match Split=n from mrna/exon?
    if($xo eq '-') { $nrev++; } elsif($xo eq '+') { $nfwd++; }
    $newstart= $xb if($newstart == 0 or $xb < $newstart);
    $newstop = $xe if($newstop == 0 or $xe > $newstop);
    $nxok++;
    }

  my $newor= ($nrev > $nfwd)?'-':'+';
  unless($newstop==0 or ($newstart == $oldstart and $newstop == $oldstop and $newor eq $oldor)) { 
    $changed++; 
    $mrna->[3] = $newstart; $mrna->[4] = $newstop; 
    $mrna->[6] = $newor; 
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




=item intron_short_fix
  
  remove too-short introns, from gmap/gsplign/other, to MININTRON, changing exons as needed
  call this BEFORE bestorf so dont have to redo.
  
  UPD1707: NOSHORTINFIX option to turn off, causes some problems ..

=cut

sub intron_short_fix  {  
  my($generec,$exontype) = @_; 
  #?? return($changed,$generec,$infixnote);  
  $exontype||="exon";
  my($infixnote,$changed)=("",0);
  return($changed,$generec,$infixnote) if($NOSHORTINFIX);  
  
  #** FIXME for isssplit, generec has all parts, ref diff
  
  #?? cut CDS for short in also ?? #ASSUME ?? NO # 
  #o# my @generec= sort _sortgene grep{ $_->[2] ne "CDS" } @$generec; # sorts by genome loc, not stranded
  my @generec= sort _sortgene @$generec; # sorts by genome loc, not stranded
  #FIXME: _sortSplitgene
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

  




=item generecfix
  # replace parts of full generec w/ udpates, ie exonfix
  # .. used at intronfix, cdna-trimnnn, elsewhere?
  # update mrnaspan .. also check CDSexons vs full exons.
  # update genespan?? need all mrna alts still for that
=cut
  
sub generecfix {
  my($generec,$exongff,$exonfix,)= @_;
  my $mrnachanged=0;
  
  my @oldexon= map { my @xnew= @$_; \@xnew; } @$exongff; # clone all
  foreach my $ex (@oldexon) { $ex->[2] = "oldexon"; $ex->[0] =~ s/^/#/; }

  my @generecfix= grep{ $_->[2] ne "exon" } @$generec;
  my @mrna = grep{ $_->[2] eq "mRNA" } @generecfix; # fixme issplit/Splitgene has more mrna parts
  my $issplit= (@mrna>1)? @mrna : 0;
  
  push @generecfix, @$exonfix;
  #d   push @generecfix, @oldexon if($debug>1); # do for check cds; if($debug>1); #??
  
  if($issplit) {
    @generecfix= sort _sortSplitgene  @generecfix; # sorts by genome loc, not stranded
  } else {
    @generecfix= sort _sortgene  @generecfix; # sorts by genome loc, not stranded
  }
    #FIXME: _sortSplitgene # ^ not quite right for CDS check: oldexon ..
        
  #o my($mrna) = grep{ $_->[2] eq "mRNA" } @generecfix; # fixme issplit/Splitgene has more mrna parts
  # FIXME Split:  BUG, sets all part-mrna to same enclosing span 
  foreach my $mrna (@mrna) {
    $mrnachanged++ if( mrna_fixspan($mrna,$exonfix) );
  }
  
  #d  @generecfix= grep{ $_->[2] ne "oldexon" } @generecfix unless($debug>1);
  $generec= \@generecfix;
  return( $generec, $exonfix, $mrnachanged);
}
  
sub phaseorf {
  my($xonseq,$cdnafull)=@_;
  return(0,0,"",$cdnafull,undef) unless($xonseq);
  my @lorfs;
  my @phase= ($cdnafull and length($xonseq)>3) ? (0,1,2) : (0);
  for my $j (@phase) { 
    my $cdp= ($j==0) ? $xonseq : substr($xonseq,$j);
    my $workseq= ($cdnafull) ? $cdnafull.$cdp : $cdp;
    # cdna_proteins.pm methods
    my @stops  = identify_putative_stops($workseq);
    my @starts = identify_putative_starts($workseq,\@stops);
    my @orfs = get_orfs (\@starts, \@stops, $workseq, '+');
    ## use goodlen not length, no gaps here
    my($longorf) = sort {$b->{goodlen} <=> $a->{goodlen} or $b->{complete} <=> $a->{complete}} @orfs;
    if($longorf) { push @lorfs, [$longorf->{goodlen},$j,$cdp,$workseq,$longorf]; }
    }
  # return (0,0,"",$cdnafull,undef) unless(@lorfs);
  my ($lorf)= sort{ $$b[0] <=> $$a[0] or $$a[1] <=> $$b[1] } @lorfs;
  return ($lorf) ? @$lorf : (0,0,"",$cdnafull,undef);
  # my ($llen,$lj,$lexon,$lcdna,$longorf)= @$lorf;
  # return($llen,$lj,$lexon,$lcdna,$longorf); # ==  @$lorf
} 

# 201511. add completeCDS extend-ends option, ix==0 and ix==$nexon1, look for completing cdsbases of partial cds
## 1707 FIXME: get computed proStart,Stop from  exons; dont use oldProstart,stop : or in completeCDSb()
## ** this is bad, causing shorter/partials vs not recalc
use constant RECALC_CDSOFF1707 => 0;

sub completeCDSb { 
  my($geneid, $cdnafull, $exongff, $revgene, $oldProstart5, $oldProend3, $oldprot)=@_;
    # drop params: $oldoffs,$oldStart_b,$oldStart_e; ADD: oldprot
    ## drop bad ,$EXTcdna return
    # ($extbest,$EXTexongff,$extaddattr,$ext5,$ext3)= 
    #   completeCDSb($geneid, $cdna, $cdnaexons, $revgene, $oldProstart5,$oldProend3);
  
  my $addattr="";
  my $XOFF=270;  # 180 too short ??
  my $UTRSLOP=27; # dont offset where UTR exists beyond partial CDS end ?    

  if(RECALC_CDSOFF1707) { #?? this may be bad, causing shorter/partials vs not??
    $KEEPSAMECDS=1; ## cdna_proteins:getBestProt2(), need for oldStart_b
    my($protCOMP, $prostart5COMP, $proend3COMP)= getBestProt("partial", $cdnafull, $exongff, $oldProstart5, $oldProend3);
    $KEEPSAMECDS=$samecds; ## cdna_proteins:getBestProt2();  
    $oldProstart5= $prostart5COMP;
    $oldProend3 = $proend3COMP;
  }
  
  ## bug? reduced complete count for -rev; check/reset?
  my @xbloc= @{$exongff->[ 0]}[0,3,4,6,8]; # ($xr,$xb,$xe,$xo,$xa)
  my @xeloc= @{$exongff->[-1]}[0,3,4,6,8];
  if($xbloc[3] eq "-") { $addattr.=",rerev" unless($revgene); $revgene=1; }
  else { $addattr.=",unrev" if($revgene); $revgene=0; }
  if( ($revgene and $xbloc[1] < $xeloc[1]) or (not $revgene and $xbloc[1] > $xeloc[1]) ) { 
    $addattr.=",resort";  my @xrev= reverse @$exongff; $exongff=\@xrev; }
  if( ($revgene and $oldProstart5<$oldProend3) or (not $revgene and $oldProstart5>$oldProend3) ) { 
    $addattr.=",revcdsbe"; ($oldProstart5,$oldProend3)= ($oldProend3,$oldProstart5); }
  
  my $nexon1= scalar(@$exongff) - 1;
  our $cdnaext= $cdnafull;  # only extend ends..
  my @extongff= @$exongff; # copy all
  my ($ext5try,$ext3try)=(0,0); 

  #NOTE exons here are rev-sorted, exon[0] = end5 of CDS
  sub exoffb{ my($ex,$ic,$xof)=@_; $ex->[$ic]+= $xof; $ex->[8]=~s/$/;cexoff=$xof/; }
  sub exoffadd{ 
    my($ex,$ic,$jc,$xof)=@_; 
    our $cdnaext;
    my($ob,$oe,$oo)= @{$ex}[3,4,6];
    my $xend = ($jc < $ic)?1:0; 
    my $xapend=($oo eq "-")? !$xend : $xend;
    $ex->[$ic]+= $xof; 
    $ex->[$jc]= ($xend)?$oe+1:$ob-1; # shift to end point
    my($addseq)= getexondna($ex, DOREVCOMP); 
    my $ret= length($addseq); $ret= -$ret if($xof<0);
    if($xapend) { $cdnaext= $cdnaext . $addseq;  } # bad for revgene.. needs opposite
    else { $cdnaext= $addseq . $cdnaext;  } 
    $ex->[$jc]= ($xend)?$ob:$oe; # revert non-end
    $ex->[8]=~s/$/;cexoff=$xof/; 
    return $ret;
    }

  # only these: if($ix==0 or $ix==$nexon1)  
  for my $ix (0,$nexon1) {
    my $xft= $exongff->[$ix];
    my($xb,$xe)= @{$xft}[3,4];  
    my $cdnaexoni= [@$xft]; # @xclone= @$xft; $cxi= \@xclone;  # will be changes to this here; clone?
    $extongff[$ix]= $cdnaexoni;
    if($revgene) {
      if($ix==0 and ($xe <= $oldProstart5 + $UTRSLOP)) {
        $ext3try= exoffadd($cdnaexoni,4,3,$XOFF);
      }
      if($ix==$nexon1 and ($xb >= $oldProend3 - $UTRSLOP)) { 
        $ext5try= exoffadd($cdnaexoni,3,4,-$XOFF);
      }
   
    } else {
      if($ix==0 and ($xb >= $oldProstart5 - $UTRSLOP)) {
        $ext5try= exoffadd($cdnaexoni,3,4,-$XOFF);
      }
      if($ix==$nexon1 and ($xe <= $oldProend3 + $UTRSLOP)) { 
        $ext3try= exoffadd($cdnaexoni,4,3,$XOFF);        
      }
    }
  }

  if($cdnaext ne $cdnafull) {
    # oldStart_b needs cdna_proteins:KEEPSAMECDS; prefer keep same CDS exons but can extend/shorten protein bounds
    $KEEPSAMECDS=1; ## cdna_proteins:getBestProt2();  
    my($orfprot, $prostart5, $proend3, $bestorf)= 
      getBestProt("partial", $cdnafull, $exongff); ##, $oldStart_b,$oldStart_e);
    my($orfproti, $prostart5i, $proend3i, $bestorfi)= 
      getBestProt("partial", $cdnaext, \@extongff );  
    $KEEPSAMECDS=$samecds; ## cdna_proteins:getBestProt2();  

    if($oldprot) {
      ## FIX for bestaa=pubaa, oldprot .. must test agains that; BUT problem pubaa: XXXX, other?
      my $sameold= index($orfproti,substr($oldprot,2));
      if($sameold<0 and $oldprot=~/XX/) { # check gaps?
        my @oldaa= split /X+/, $oldprot; my @oai;
        for my $oaa (@oldaa) { my $i=index($orfproti,$oaa); push @oai, $i if($i>=0); }
        if(@oai == @oldaa) {
	        my($ob,$oe)= @oai[0,-1]; my $olen=$oe + length($oldaa[-1]) - $ob;
          $oldprot= substr($orfproti,$ob,$olen);
	        }
      }
      $orfprot=$oldprot;
    }

    ## test stats
    my($aawi,$aaw)= (length($orfproti),length($orfprot));
    my $aadif= $aawi - $aaw; 
    my($trwi,$trw)= (length($cdnaext) , length($cdnafull)); 
    my($aaci,$aac)= ($bestorfi->{complete}, $bestorf->{complete});
    my($relstart5i,$relend3i)= ($prostart5i + $ext5try, $proend3i + $ext5try); # relative to -XOFF shift
    my $ooff= join "-", $prostart5, $proend3;
    my $noff= join "-", $relstart5i,$relend3i;
    if($relstart5i <= 0) { $noff="EXT$noff"; }
    
    my $extbest=0;
    $extbest=1 if($aawi>$aaw); # not enough, for part[53] need other offset same
    $extbest=0 if($relend3i < $proend3-$UTRSLOP); # diff locations, happens w/ for same prot
    $extbest=0 if($relstart5i > $prostart5+$UTRSLOP);
    my $sameprot= (index($orfproti,substr($orfprot,2))<0)?0:"AASAME";  # substr() cuts possible part5 codon
    $extbest=0 unless($sameprot);
    # if($extbest and $sameprot) { $extbest=0; } # substr() cuts possible part5 codon

    warn "#cds.complete $geneid fix=$extbest aaquals=$aaci/$aac, aadif=$aadif,$aawi/$aaw,$sameprot, "
        . "coff=$noff/$ooff, trw=$trwi/$trw, at=$addattr \n"; ## if $debug; , oldaaq=$oldaalen, oldoffs=$oldoffs; 

    if($extbest) { 
      ## trim off extra XOFF at ends, leave NO UTR or tiny UTR ? FIXME revor swap *
      my($ext5,$ext3)=(0,0); # use $ext5try,$ext3try above?
      my($cdsOfExt, $cdsattrOfExt)= (undef,"");
      ## try again:
      ($cdsOfExt, $cdsattrOfExt)= getCDSgff2( \@extongff, $bestorfi); # see above
      #o ($cdsOfExt, $cdsattrOfExt)= getCDSgff( \@extongff, orfParts($bestorfi)); # see above
      my($cdsb,$cdse)= $cdsattrOfExt=~m/cdsspan=(\d+).(\d+)/;
      if($cdse) {
        my($oxb,$oxe);
        my $crev=0; if($cdse < $cdsb) { $crev=-1; ($cdsb,$cdse)=($cdse,$cdsb); }
        if($revgene) {  
          $oxb= $exongff->[-1]->[3];  $oxe= $exongff->[0]->[4]; # rev sort
          if($cdsb<$oxb) { $ext3= -($cdsb - 3 - $oxb); } ## -() for rel..5i compat, -ext offset
          if($cdse>$oxe) { $ext5= -($cdse + 3 - $oxe); }
          }
        else { 
          $oxb= $exongff->[0]->[3];  $oxe= $exongff->[-1]->[4]; 
          if($cdsb<$oxb) { $ext5= $cdsb - 3 - $oxb; }
          if($cdse>$oxe) { $ext3= $cdse + 3 - $oxe; }
          }
      }
      
      if($ext5 == 0 and $ext3 == 0) {
        # NOTE exongff are rev order here; 5' is first
        if($relstart5i <= 0) { $ext5= $relstart5i - 3; } # 
        if($relend3i >= $trw){ $ext3= $relend3i + 3 - $trw; } # is this proper offset?
      }
      
      if($ext5 != 0) { if($revgene){ exoffb($exongff->[0],4,-$ext5); } else { exoffb($exongff->[0],3,$ext5); } }
      if($ext3 != 0) { if($revgene){ exoffb($exongff->[-1],3,-$ext3); } else { exoffb($exongff->[-1],4,$ext3); } }
 
      ## TEST bug, return @extongff  == _fcdsext12.gff, vs off for _fcdsext13.gff
      ## result no change with ocds=eqn,eqaa; BUT cdsfix has aanew > aaold
      ## ** BUG maybe happens when ext5 == 0 and ext3 == 0, but offsets changed. ??
      # my $savex= $exongff; $exongff= \@extongff;
      
      ## *** FIXME cdnaext is wrong here; repull from exongff   ***
      ##caller#($cdnaext)= getcdna2( $exongff, 1, 0, 0, CDNA_TRIMNNN);
      ## DROP bad return: ,$cdnaext
       
      $addattr = "cdsfix=ext:$aaci/$aac,$noff/$ooff,$aawi/$aaw$addattr;";
      return($extbest,$exongff,$addattr,$ext5,$ext3);
      }
  }
  #?? debug return cdsnofix= attr??
  
  return(0);
}

sub completeCDS { 
  my($geneid,$exongff,$revgene,$mrnaat,
     $oldaalen,$oldoffs,$oldProstart5,$oldProend3,$oldStart_b,$oldStart_e)=@_;

  ## FIXME: revise to work w/ getcdna2() output (cdnaseq and exon.gff),
  ##  .. only extend cdnaseq ends from exon ends where CDS abuts, then test for complete aa
  ##  .. insert completeCDS() AFTER splitgene/nosplit handling, and mod here to know we have CDS end point in part
  
  my $XOFF=270;  # 180 too short ??
  my $UTRSLOP=27; # dont offset where UTR exists beyond partial CDS end ?
  my ($cdnafull,$cdnaext,$ext5try,$ext3try)=("","",0,0); 
  my ($oldcdsb,$oldcdse)= $oldoffs=~m/(\d+).(\d+)/; # fail unless have oldoffs?
  my $nexon1= scalar(@$exongff) - 1;
  my $dok=($mrnaat=~/cdsfix=complete/)?1:0; # skip this opt?
  if(!$dok and $oldaalen=~/partial/ and $DO_CDSCOMPLETE>1) { $dok=1; } #? default
  return(0) unless($dok);
  # $nexon1= -1 unless($dok);
  my @extongff=();
  
  #NOTE exons here are rev-sorted, exon[0] = end5 of CDS
  sub exoff{ my($ex,$ic,$xof)=@_; $ex->[$ic]+= $xof; $ex->[8]=~s/$/;cexoff=$xof/; }
  
  for my $ix (0..$nexon1) {
    my $xft= $exongff->[$ix];
    my @xclone= @$xft; 
    my $cdnaexoni= \@xclone;  # will be changes to this here; clone?
    my($cdnai,$rev,$cleni)= getexondna($cdnaexoni, DOREVCOMP); 
    $cdnafull= uc($cdnafull . $cdnai); # unextended..
    
    #FIXME here: XOFF only when partial CDS abuts mRNA end
    # .. need $oldCDSb,$oldCDSe
    if($ix==0 or $ix==$nexon1) { 
      my($xb,$xe)= ($cdnaexoni->[3],$cdnaexoni->[4]);
      if($revgene) {
         if($ix==0 and ($xe <= $oldProstart5 + $UTRSLOP)) {
          $ext3try= $XOFF; 
          # $cdnaexoni->[4] += $ext3try; # dang rev needs [4] here
          exoff($cdnaexoni,4,$ext3try);  
        }
        if($ix==$nexon1 and ($xb >= $oldProend3 - $UTRSLOP)) { 
          $ext5try= -$XOFF; 
          # $cdnaexoni->[3] += $ext5try;  # rev needs [3] here
          exoff($cdnaexoni,3,$ext5try);  
        }
     
      } else {
        if($ix==0 and ($xb >= $oldProstart5 - $UTRSLOP)) {
          $ext5try= -$XOFF; 
          # $cdnaexoni->[3] += $ext5try; # dang rev needs [4] here
          exoff($cdnaexoni,3,$ext5try);  
        }
        if($ix==$nexon1 and ($xe <= $oldProend3 + $UTRSLOP)) { 
          $ext3try= $XOFF; 
          # $cdnaexoni->[4] += $ext3try;  # rev needs [3] here
          exoff($cdnaexoni,4,$ext3try);  
        }
      }
      ($cdnai,$rev,$cleni)= getexondna($cdnaexoni, DOREVCOMP); 
    }
    $cdnaext= uc($cdnaext . $cdnai); # extended..
    push @extongff, $cdnaexoni;
  }

  if($cdnaext and $cdnafull and $cdnaext ne $cdnafull) {
    # oldStart_b needs cdna_proteins:KEEPSAMECDS; prefer keep same CDS exons but can extend/shorten protein bounds
    # add $oldprot test ..
    $KEEPSAMECDS=1; ## cdna_proteins:getBestProt2();  
    my($orfprot, $prostart5, $proend3, $bestorf)= 
      getBestProt("partial", $cdnafull, $exongff, $oldStart_b,$oldStart_e);
    my($orfproti, $prostart5i, $proend3i, $bestorfi)= 
      getBestProt("partial", $cdnaext, \@extongff ); ## $oldStart_b,$oldStart_e;
    $KEEPSAMECDS=$samecds; ## cdna_proteins:getBestProt2(); reset

    ## test stats
    my($aawi,$aaw)= (length($orfproti),length($orfprot));
    my $aadif= $aawi - $aaw; 
    # my $gadif= $bestorfi->{goodlen} - $bestorf->{goodlen}; 
    my($trwi,$trw)= (length($cdnaext) , length($cdnafull)); 
    my($aaci,$aac)= ($bestorfi->{complete}, $bestorf->{complete});

    my($relstart5i,$relend3i)= ($prostart5i + $ext5try, $proend3i + $ext5try); # relative to -XOFF shift
    my $ooff= join "-", $prostart5, $proend3;
    my $noff= join "-", $relstart5i,$relend3i;
    if($relstart5i <= 0) { $noff="EXT$noff"; }
    
    my $extbest=0;
    $extbest=1 if($aawi>$aaw); # not enough, for part[53] need other offset same
    $extbest=0 if($relend3i < $proend3-$UTRSLOP); # diff locations, happens w/ for same prot
    $extbest=0 if($relstart5i > $prostart5+$UTRSLOP);
    my $sameprot= (index($orfproti,substr($orfprot,2))<0)?0:"AASAME";  # substr() cuts possible part5 codon
    $extbest=0 unless($sameprot);
    # if($extbest and $sameprot) { $extbest=0; } # substr() cuts possible part5 codon

    warn "#cds.complete $geneid fix=$extbest aaquals=$aaci/$aac, aadif=$aadif,$aawi/$aaw,$sameprot, coff=$noff/$ooff, "
        . "trw=$trwi/$trw, oldaaq=$oldaalen, oldoffs=$oldoffs;  \n"; ## if $debug;

    if($extbest) { 
      ## trim off extra XOFF at ends, leave NO UTR or tiny UTR ? FIXME revor swap *
      my($ext5,$ext3)=(0,0); # use $ext5try,$ext3try above?

      my($cdsOfExt, $cdsattrOfExt)= (undef,"");
      
      ## try again:
      ($cdsOfExt, $cdsattrOfExt)= getCDSgff2( \@extongff, $bestorfi); # see above
      #o ($cdsOfExt, $cdsattrOfExt)= getCDSgff( \@extongff, orfParts($bestorfi)); # see above
      my($cdsb,$cdse)= $cdsattrOfExt=~m/cdsspan=(\d+).(\d+)/;
      if($cdse) {
        my($oxb,$oxe);
        my $crev=0; if($cdse < $cdsb) { $crev=-1; ($cdsb,$cdse)=($cdse,$cdsb); }
        if($revgene) {  $oxb= $exongff->[-1]->[3];  $oxe= $exongff->[0]->[4];  }
        else { $oxb= $exongff->[0]->[3];  $oxe= $exongff->[-1]->[4]; }
        if($cdsb<$oxb) { $ext5= $cdsb - 3 - $oxb; }
        if($cdse>$oxe) { $ext3= $cdse + 3 - $oxe; }
      }
      #nohelp# ($cdsOfExt, $cdsattrOfExt)= getCDSgff2( \@extongff, $bestorfi); # cdna_proteins.pm
      ## cdsattrOfExt == "cxlen=$orflen/$trlen;aalen=$aalen,$pcds%,$compl;protein=$orfprot";
      ## see below:
      ## ($cdsgffOfExt, $cdsattrOfExt)= getCDSgff( \@extongff, $orfproti, $prostart5i, $proend3i, length($cdnaext));
      # if(ref($cdsgffOfExt)) { # test this way.
      #   my($cb,$ce,$oxb,$oxe)=(0) x 9;
      #   $cb= $cdsgffOfExt->[0]->[3];
      #   $ce= $cdsgffOfExt->[-1]->[4];
      #   ($cb,$ce)=($ce,$cb) if($cb>$ce);
      #   if($revgene) {  $oxb= $exongff->[-1]->[3];  $oxe= $exongff->[0]->[4];  }
      #   else { $oxb= $exongff->[0]->[3];  $oxe= $exongff->[-1]->[4]; }
      #   if($cb<$oxb) { $ext5= $cb - 3 - $oxb; }  ## utrexon -1 of below method..
      #   if($ce>$oxe) { $ext3= $ce + 3 - $oxe; }
      # }
      
      if($ext5 == 0 and $ext3 == 0) {
      # NOTE exongff are rev order here; 5' is first
      #.. messy calc, use instead @extongff[0,-1], $prostart5i,proend3i are offsets of those exons
      if($relstart5i <= 0) { $ext5= $relstart5i - 3; } # 
      if($relend3i >= $trw){ $ext3= $relend3i + 3 - $trw; } # is this proper offset?
      #? if($relend3i > $proend3){ $ext3= $relend3i + 3 - $proend3; } # is this proper offset?
      }
      
      # if($ext5 != 0) { if($revgene) { $exongff->[0]->[4] -= $ext5; } else { $exongff->[0]->[3] += $ext5; } }
      # if($ext3 != 0) { if($revgene) { $exongff->[-1]->[3] -= $ext3; } else { $exongff->[-1]->[4] += $ext3; } }
      if($ext5 != 0) { if($revgene){ exoff($exongff->[0],4,-$ext5); } else { exoff($exongff->[0],3,$ext5); } }
      if($ext3 != 0) { if($revgene){ exoff($exongff->[-1],3,-$ext3); } else { exoff($exongff->[-1],4,$ext3); } }
      
      ## TEST bug, return @extongff  == _fcdsext12.gff, vs off for _fcdsext13.gff
      ## result no change with ocds=eqn,eqaa; BUT cdsfix has aanew > aaold
      ## ** BUG maybe happens when ext5 == 0 and ext3 == 0, but offsets changed. ??
      # my $savex= $exongff; $exongff= \@extongff;
      
      my $addattr= "cdsfix=ext:$aaci/$aac,$noff/$ooff,$aawi/$aaw";
      return($extbest,$exongff,$addattr,$ext5,$ext3,$relstart5i,$relend3i);
## CALLER does
#       ($cdna, $cdnalen, $cdnaexonseq, $xfixfix, $xtrimmed)= getcdna2( \@exongff, 1, 0, 0, CDNA_TRIMNNN); 
#       $xtrimmed++; # * REGISTER prior exon change for  FIXMEd, CDS extended NOT exons ..
#       $addattr .= ";cdsfix=ext:$aaci/$aac,$noff/$ooff,$aawi/$aaw"; ## if($xfixed);

      }
  }
  return(0);
}  

  
=item testgene : big, complex, messy sub

  v2: remove intron test parts
  
  xFIXME : check/remove existing CDS; compare to new 
  FIXME2: exist CDS : allow for problem cases like end-of-scaffold partials
  FIXME3: -samecds not much use in test w/ augustus calls.
  FIXME4: stop changes that replace all CDS with completely new CDS, happens for transposon-spans
         .. not quite -samecds, but -nocompletelydifferentcds
  FIXME5: -samecds needs to use phase0 : partial5 from AUG uses this (eg. at scaf startpos=1)
  FIXME6: chimera, multimap IDs have 2ndary tag: _C[12] _G[2..n] ; for cdnain fix this

  FIXME7: 2015.01, Fix/remove tiny introns, exon-overlap => join exons? w/ some gap syntax?

=cut
 
=item testgene notes

  #  FIX2: add test intron overlaps : retained = intron inside exon; err = intron rev at splice/inside
  #  FIX3? add here option to filter out huge-span asmrna with poor qual : no intron evidence of long intron; poor prot
  #  FIX? for gaps, chomp off end gaps NNN of getcdna() ? or remove from orfprot
  #  .. need to change exongff also for cdna trimnnn
  
  #ov2: $issplit: require caller to collect all mRNA parts first, generecIN should contain all exons?
  #  .. but also need all mRNA parts, $mrna > @mrna for splits?

  #ov2: $issplit: add cdna paste of Split parts .. collect all parts, then process cdna..
  #ov2: FIXME split parts paste with 3 phases, ie cant assume break point is accurate; add N,NN shifts
  # .. need to test protrans for 3 diff phases for each part-join, ugh..
  # .. return @cdnaparts, do pastes+NN, protrans tests till find best prot ~= oldaalen
  # .. ie. 2 parts = 3 pro tests, 3 parts = 3*3 = 9 tests, 4 parts = 3*3*3 ..
  # .. cancel fixstrand loop if issplit?

  #ov2: FIXME split parts paste with 3 phases, ie cant assume break point is accurate; add N,NN shifts
  #.. need special processing to find phases for joining cdna parts.. ugh ugh
  #.. improved 200 split genes to >= trasm prot for kfish2nsplign15n_fcds4d (800 vs 600 of 3000 splits)
  #.. but also improved shorter ones some. ave split aa size now? vs before
  #.. Splits  fcds4d: n=3047, n90%=711,  aaw=331, aow=546; fcds4c: n=3028, n90%=231,  aaw=245, aow=546 .. improved
  #.. nosplit fcds4d: n=9462, n90%=2918, aaw=363, aow=506; fcds4c: n=9270, n90%=2726, aaw=364, aow=509 .. same

=cut

=item DO_CDSFIX work

  # $DO_CDSFIX works something like SPLITGENEPROTFIX, but per cds-exon test orf with added exondna
  #  .. trim or drop cds-exons that damage orig protein, ie inside  oldoffs span, ignore utr exons?
  #  .. where "damage" is fuzzy measure, 
  #  .. replace all of below SPLIT.FIX, return ($cdna,$cdnalen, $cdnaexonseq, $xfixfix, $xtrimmed)
  # .. should use this only on trasm genomap identified as bad cds2aa but w/ full-ish cdna map


  # * first try past bugs worked for 4 cases of long pubaa > short cds .. got middle way to pubaa
 test1.bestof6top5
 Funhe2EKm010574t1	-2433.da	1924,41%,complete	4423,84%,complete,gaps:66  
 Funhe2EKm013880t1	-3011.da	1908,51%,complete	4983,95%,complete,gaps:64  
 Funhe2EKm027743t1	-3194.da	1469,34%,complete	4663,97%,complete	89       
 Funhe2EKm038182t1	-2268.da	1516,40%,complete	3784,97%,partial3	95

 kfish2nsplign15n_test1fcds5.gff
 Funhe2EKm010574t1 rev
   trg=1 15650;gaps=1814,MGap:275-486,MGap:7660-7926,MGap:2234-2608,MGap:12338-12720,MGap:6714-7290,;
   offs=832-14103;aaold=4423,84%,complete
   cxlen=11355/12989;aalen=3785,87%,partial3;utrx=1,4;ocds=NEn,..;
   cdsoff=617-11971;cdsfix=xtrimb7:35,xtrimb7:5,xdrop19:5954-6121,xtrimb20:23,xtrimb57:131,xtrimb57:74,xtrimb57:26,xtrimb58:5,xtrime61:50,;xtrim=1;
  
 Funhe2EKm013880t1 rev
   trg=1 15577;gaps=522,RGap:15578-15587,MGap:4070-4193,MGap:10073-10240,MGap:3707-3901,MGap:10260-10284,;
   offs=163-15114;aaold=4983,95%,complete
   cxlen=12006/12517;aalen=4002,95%,partial;
   cdsoff=1-12006; cdsfix=xdrop5:3996-4062,xtrimb7:18,xdrop8:4265-4323,xdrop9:4265-4293,xdrop37:9360-9369,xtrimb38:33,xdrop38:9360-9609,xdrop39:9360-9426,xdrop40:9360-9753,xtrimb41:42,xdrop41:9360-9447,xtrimb42:9,xtrime42:18,xdrop43:9378-9423,xtrimb44:9,xdrop44:9378-9474,xtrimb62:92,xtrimb62:74,xtrimb62:14,;xtrim=1;
   
  Funhe2EKm027743t1 fwd
    trg=1 14283;gaps=1471,MGap:4477-4893,MGap:6427-6921,polyA:14284-14288,MGap:6933-7486,;gescore=86;
    offs=31-14022;aaold=4663,97%,complete
    cxlen=11295/11648;aalen=3765,96%,partial;
    cdsoff=1-11295;cdsfix=xtrime20:33,xdrop29:5793-5838,xtrime30:39,xtrimb31:12,xtrimb31:36,xdrop34:6279-6333,xdrop48:8398-8604,xdrop69:11296-11424,;
    
  Funhe2EKm038182t1 fwd
    trg=2 11657;gaps=502,LGap:1-1,MGap:5981-6203,MGap:4660-4937,;
    offs=305-11656;aaold=3784,97%,partial3
    cxlen=10578/10579;aalen=3526,99%,partial;
    cdsoff=1-10578;cdsfix=xdrop0:1-194,xtrimb33:161,xtrimb33:47,xtrimb35:17,;xtrim=1;

try1h, best so far
Funhe2EKm010574t1;cov=86%,13469/15650;pid=99.3;nexon=62;splice=214;splicemix=214,16;trg=Funhe2EKm010574t1 1 15650;
	gaps=1814,MGap:275-486,MGap:7660-7926,MGap:2234-2608,MGap:12338-12720,MGap:6714-7290,;
	clen=15650;offs=832-14103;cxlen=11382/13972;	
	aalen=3794,81%,complete;cdsoff=1330-12714;utrx=3,4;;cdsfix=xtrimb57:12014-13017/135,xtrimb57:11987-12882/108,xtrimb57:12014-12774/135,;xtrim=1;	
	aaold=4423,84%,complete
Funhe2EKm013880t1;cov=96%,14997/15587;pid=99.5;nexon=63;splice=240;trg=Funhe2EKm013880t1 1 15577;
	gaps=522,RGap:15578-15587,MGap:4070-4193,MGap:10073-10240,MGap:3707-3901,MGap:10260-10284,;	
	clen=15587;offs=163-15114;cxlen=14493/14493;	
	aalen=4831,100%,partial;cdsoff=1-14493;;cdsfix=xtrimb62:14591-15049/97,xtrimb62:14582-14952/88,xdrop62:14494-14864/14619,;xtrim=1;	
	aaold=4983,95%,complete
Funhe2EKm027743t1;cov=89%,12742/14288;pid=99.3;nexon=70;splice=266;trg=Funhe2EKm027743t1 1 14283;
	gaps=1471,MGap:4477-4893,MGap:6427-6921,polyA:14284-14288,MGap:6933-7486,;	
	clen=14288;offs=31-14022;cxlen=7500/12069;	
	aalen=2500,62%,partial3;cdsoff=4570-12069;utrx=21,0;;cdsfix=xtrime20:4212-4464/6,xdrop48:8923-9272/9129,xdrop69:11823-12207/12105,;xtrim=1;	
	aaold=4663,97%,complete
	#^^ problem case, cdsoff=31-4440 also has orf .. should be 1 long gene
Funhe2EKm038182t1;cov=95%,11093/11657;pid=99.4;nexon=82;splice=316;trg=Funhe2EKm038182t1 2 11657;
	gaps=502,LGap:1-1,MGap:5981-6203,MGap:4660-4937,;	
	clen=11657;offs=305-11656;cxlen=10905/11155;	
	aalen=3635,97%,partial3;cdsoff=245-11149;;	
	aaold=3784,97%,partial3

=cut


sub testgene  
{
  my($generecIN, $geneother)= @_;
  my($addattr,$changed,$mrnachanged)=("",0,0);
  my($cdnaoutg,$orfcdsg,$orfprotg)=("") x 3; #output seqs outside dang fixstrand loop
  
  my @oldCDS = grep{ $_->[2] eq "CDS" } @$generecIN;

  ## do this before bestorf...  call on generecIN, no call after USE_CDSEXONS 
  ## BUT ugh reusing generecIN below
  my($xchanged, $xgenerec, $xfixnote)= intron_short_fix($generecIN); 
  if($xchanged) {  
    #.. use generecfix here instead?
    # my($igenerec, $iexongff, $mrnachanged1)= generecfix( \@generec, \@exongff, $xfixfix); 
    # @generec= @$igenerec; @exongff= @$iexongff; $mrnachanged += $mrnachanged1;
 
    $generecIN= $xgenerec;
    $changed+= $xchanged; $addattr .=";" if($addattr); $addattr .= "inshort=$xfixnote";
  } 


  my @mrna= grep{ $_->[2] eq "mRNA" } @$generecIN;
  my $issplit= (@mrna>1) ? @mrna : 0;
  my @generec;
  if($issplit) {
    @generec= sort _sortSplitgene grep{ $_->[2] ne "CDS" } @$generecIN; # sorts by genome loc, not stranded
  } else {
    @generec= sort _sortgene grep{ $_->[2] ne "CDS" } @$generecIN; # sorts by genome loc, not stranded
  }
  
  my @exongff= grep{ $_->[2] eq "exon" } @generec;
  #ov2: my($mrna)  = grep{ $_->[2] eq "mRNA" } @generec;
  # my @mrna= grep{ $_->[2] eq "mRNA" } @generec;
  my $mrna= $mrna[0]; 
  my $mrnaat= $mrna->[8];
  ## my $issplit= (@mrna>1) ? @mrna : 0;
  #or? my $issplit=($mrnaat =~ m/;Split=([^;\s]+)/)? $1 : ($geneid =~ /_C(\d+)$/)? $1 : 0;
  my $geneid= ($mrnaat =~ m/ID=([^;\s]+)/)? $1 : ""; # make one up?
  
  unless($mrna and @exongff > 0) {
    if($USE_CDSEXONS and $mrna and @oldCDS > 0) { # or caller can s/CDS/exon/ easily, or duplicate like this
      foreach my $oc (@oldCDS) { my @doc= @$oc; $doc[2]="exon"; push @exongff, \@doc; }    
      push @generec, @exongff;
    } else {
      putgene($generecIN, $geneother,"err=Missing-mrna-exon"); # flag="err=Missing-mrna-exon"
      return 0;
    }
  }
  
  ## FIXME: ref == NOPATH .. need cdna-only results, modify gff to use cdna-seq as scaffold?
  if($mrna->[0] =~ /^NOPATH/) {
    # handleNoPath();
    putgene($generecIN, $geneother,"err=Missing-mrna-location");  
    return 0;
  }
  
  # my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid)= @$mrna;
  my $gstrand= $mrna->[6]; # FIXME for "." ; add CDS bestpro strand (always +??)
  my(@genefwd, @generev);
  my $fixstrand=0; # *** FIXME 2011Dec   
  $fixstrand=1 if($gstrand eq ".");# need to test bestpro both strands
  #? $fixstrand=1 if($tattr=~m/;sense=-1/); # gmap wrong for mRNA mapping; fixme 1605
  $fixstrand=0 if($issplit); # not both..
 
  my($addattrIN,$changedIN)=($addattr,$changed);
  my ($changed0,$changed1)=(0,0);
  for( my $ifix= 0; $ifix <= $fixstrand; $ifix++) { 
    # BIG loop testing both ways, save protfwd, protrev
    $changed=$changedIN; $addattr=$addattrIN; # clear for loop step
    
    if($fixstrand) {
      $gstrand=($ifix == 1) ? "-" : "+";
      ## must clone generec for save @genefwd, @generev
      @generec= sort _sortgene grep{ $_->[2] ne "CDS" } @$generecIN; # sorts by genome loc, not stranded
      @generec= map { my @xnew= @$_; \@xnew; } @generec; # clone all
      map{ $_->[6]= $gstrand } @generec;
      ($mrna) = grep{ $_->[2] eq "mRNA" } @generec;
      $mrnaat= $mrna->[8];
      @exongff= grep{ $_->[2] eq "exon" } @generec;

      if( $ifix == 1 ) {  @generev= @generec; } 
      else {  @genefwd= @generec; }
    }  
    
  
  $geneid= ($mrnaat =~ m/ID=([^;\s]+)/)? $1 : ""; # make one up?
  (my $geneidfix= $geneid) =~ s/_[CG]\d+$//; # chimera/splitgene _C[12] and multimap _Gnnn id tags
  #ov2: Split gene work here, geneidfix == true gene id.
  #no, see above# my $issplit=($mrnaat =~ m/;Split=([^;\s]+)/)? $1 : ($geneid =~ /_C(\d+)$/)? $1 : 0;
  
  my $oldprot= ($mrnaat =~ m/protein=([^;\s]+)/) ? $1 : "";
  my $oldisbest=($oldprot =~ /\w\w/ and $mrnaat =~ m/(bestaa=pubaa)/)?$1:"";
  my $oldaalen= ## evg variants and gmap.gff stutter
    ($mrnaat =~ m/aalen=(\d+,[^;\s]+)/) ? $1 :
    ($mrnaat =~ m/aaSize=(\d+)/) ? $1 :
    ($mrnaat =~ m/aalen=([^;\s]+)/) ? $1 : ""; # FIXME gmap dup aalen and aaSize<>aalen
  my $oldoffs= ($mrnaat =~ m/(?:cdsoff|offs)=([^;\s]+)/) ? $1 : "";
  
  my $revgene=($gstrand eq "-")?1:0; # NOTE exongff are rev order here; 5' is first
  if($issplit) { # sort by parts, by exon?
    @exongff= sort __revSplitGene @exongff;
  } elsif($revgene) {
    @exongff= reverse @exongff; #?? issplit
  }
  
  # FIXME5: -samecds needs to use phase0 : partial5 from AUG uses this (eg. at scaf startpos=1)
  my ($oldStart_b,$oldStart_e,$oldProstart5,$oldProend3)= (0,0,0,0); 
  if(@oldCDS){ #old:  and $samecds #?? drop samecds requirement?
    if($issplit) {
      @oldCDS= sort _sortSplitgene @oldCDS; # dang2    
    } else {
      @oldCDS= sort _sortgene @oldCDS; # dang
    }
    $oldProstart5= ($revgene) ? $oldCDS[-1]->[4] : $oldCDS[0]->[3]; #$oldprostart5 below
    $oldProend3  = ($revgene) ? $oldCDS[0]->[3] : $oldCDS[-1]->[4]; #$oldproend3
    foreach my $ex (@exongff) { my($b,$e)= ($ex->[3], $ex->[4]);
      if($b <= $oldProstart5 and $e >= $oldProstart5) { ($oldStart_b,$oldStart_e)=($b,$e); last;}
    }
  }

  ## these are updates by various tests, REGISTER $xtrimmed++ when other vals updated.
  my($cdna, $cdnalen, $cdnaexonseq, $xfixfix, $xtrimmed)= (0) x 9;

      # 201511. add extend-ends option, ix==0 and ix==$nexon1, look for completing cdsbases of partial cds
      # can do that inside DO_CDSFIX? NO, do separate... with what flags/opts?
      # ** DO_CDSCOMPLETE completeCDS() should follow SPLITGENEPROTFIX call?
  # if($DO_CDSCOMPLETE) {
  #   # add $oldprot param
  #   my($extbest,$EXTENDEDexongff,$extaddattr,$ext5,$ext3,$relstart5i,$relend3i)= 
  #     completeCDS($geneid,\@exongff,$revgene,$mrnaat,
  #        $oldaalen,$oldoffs,$oldProstart5,$oldProend3,$oldStart_b,$oldStart_e);
  #   if($extbest) {
  #     $oldStart_b= $oldStart_e= 0; #FIX: zero out $oldStart_b,$oldStart_e so wont force final CDS below
  #     @exongff= @$EXTENDEDexongff; # ? dont need; completeCDS() updates input \@exongff
  #     ($cdna, $cdnalen, $cdnaexonseq, $xfixfix, $xtrimmed)= getcdna2( \@exongff, 1, 0, 0, CDNA_TRIMNNN); 
  #     $xtrimmed++; # * REGISTER prior exon change for  FIXMEd, CDS extended NOT exons ..
  #     $addattr .= ";$extaddattr" if($extaddattr); 
  #   }
  # }
  
  # DO_CDSFIX problems on splitgenes .. so far no large improvement seen on those.
  if($DO_CDSFIX) {
    my (@cdnap,@cdnaxp);
    my ($cdnafull, $cdnaflen, $cdnafend)=("",0,0); 
    my ($oldcdsb,$oldcdse)= $oldoffs=~m/(\d+).(\d+)/; # fail unless have oldoffs?
    my $nexon1= $#exongff;
    my $xfixed="";
    $nexon1= -1 unless($oldcdse and $oldcdse>0);
    for my $ix (0..$nexon1) {
      my $xft= $exongff[$ix];
        # want Target=id xstart xend info, but may not be there..
      my @xclone= @$xft; 
      my $cdnaexoni= \@xclone;  # will be changes to this here; clone?

      my $xfixi="";
      my($xtgd,$xtgb,$xtge)=(0,0,0);
      if($$xft[8] =~ /(?:Target|trg)=(\S+).(\d+).(\d+);/) { ($xtgd,$xtgb,$xtge)=($1,$2,$3); }
      
      my($cdnai,$rev,$cleni)= getexondna($cdnaexoni, DOREVCOMP); 
      #DOREVCOMP# $cdnai = revcomp($cdnai) if($rev);  # getexondna leaves revcomp, trimNNN to caller now..
      # no trimNNN == CDNA_TRIMNNN 
      
      $cdnafend= length($cdnafull);
      my $cdnafbeg= 1 + $cdnafend;
      $cdnafend += $cleni;
      ($xtgb,$xtge)= ($cdnafbeg,$cdnafend) unless($xtge);
      
      #?? should use xtgb,xtge here .. mapping target spans of exon, not local cdna string span
      my $overcds=   ($xtge > $oldcdsb and $xtgb < $oldcdse) ? 1 : 0;
      # my $overcds=   ($cdnafend > $oldcdsb and $cdnafbeg < $oldcdse) ? 1 : 0;
      my $insidecds= ($xtgb >= $oldcdsb and $xtge <= $oldcdse) ? 1 : 0;
      my $utrbeg= ($cdnafbeg < $oldcdsb)?1:0;
      my $utrend= ($cdnafend > $oldcdse)?1:0;
      
      my $okcds=0; # need oldaalen, oldoffs cds offs / span to validate
      unless($overcds) { 
        $okcds=1;  # utr??, only append cdnai
        #NObelow# $cdnafull .= $cdnai;
        #Below#   push @cdnaxp, $cdnaexoni;
      }
      for(my $try=0; $okcds == 0 and $try < 3; $try++) {
        $okcds=1; # if(longorf spans orig cds offset ?)
        ## merge below get_orfs() here .. changing cdnaexoni as needed to maintain orf spanning cds.
        # my $nn=''; # or NN, N phase shifts?

        my($workseq,$longorf,$llen,$lj,$lexon);
        ($llen,$lj,$lexon,$workseq,$longorf)= phaseorf( $cdnai, $cdnafull);
        $cdnai= $lexon; $cleni= length($lexon);
        if($lj) { $cdnafbeg += $lj;
          $xtrimmed++; #<< need to register exon change for mrna update; also? flag in exon attr?
          if($rev) { $$cdnaexoni[4] -= $lj; } else { $$cdnaexoni[3] += $lj; } 
        }
        
    
        my $BEGCUT= 0.33; # 0.25; # : 0.15,0.25,0.33
        my $ENDCUT= 0.66; # 0.75; # : 0.85,0.75,0.66

        $cdnafend= length($workseq);
        # orf fields: length, goodlen, complete, start,stop, orient, innerstop, protein, sequence == cds
        my $lend= ($longorf) ? $longorf->{stop} : 0; $lend||=0; # stop too soon?
        if($lend < $cdnafbeg) { # this means exon is utr, after orf, treat like not overcds ?
          $lend= $cdnafbeg + int(0.10 * $cleni); # punt?
        }  
        
        my $xokspan= $lend - $cdnafbeg; #? +2? is lend/stop 1st of codon3 ?
        $xokspan=0 if($xokspan<0);  ## ^^ NEGative xokspan ?? cdnafbeg supposed to be exon start in workseq. lend orf stop in workseq 

        if(not $overcds) { $okcds=1; } #? insidecds 
        elsif($lend <= $cdnafbeg) { #bug where?
          $okcds=1; # this is utr exon, not overcds
          # .. bug dont understand yet, see lend<fbeg for some large, ok exon additions Funhe2EKm002468t1
          # .. neither of these gets right answer.
          
        } elsif($lend < $cdnafend-2) { #  - 3?
          # trim end cdnai to lend?
          # trim start cdnai if lend near exon start = cdnafbeg
          # drop exon if longorf has minimal exon cover and  $isinsidecds
          # $okcds=0;
          
          ## no certain xokspan cut val, try both cut and drop, but need next exon to test drop effect on orf
          
          if($cleni < 39) { #? bad?
            $okcds= -2; # drop, exon too short to cut
          
          } elsif($cdnafend > $oldcdse and $xokspan >= $BEGCUT*$cleni) { 
            $okcds= 1;   # last cds exon, dont cut unless start is poor
          
          #? } elsif( $cdnafbeg < $oldcdsb  and $xokspan <= $BEGCUT*$cleni) {
          #?  $okcds= 1;   # first cds exon, cut at end but not start...
            
          } elsif($xokspan <= $BEGCUT*$cleni) { #? $lend < 9 + $cdnafbeg : 0.15,0.25,0.33 give diff effects. what?
            $okcds=0; 
            my $trimb= $xokspan + 3; #?? need to do phase checks w/ trims; simple orf check here w/ trim?
            ## ^^ NEGative trimb ??
            # $cdnai= substr($cdnai,$trimb); 
            my($llen,$lj,$lexon,$lcdna,$longorf)= phaseorf( substr($cdnai,$trimb), $cdnafull);
            $cdnai= $lexon; $trimb += $lj; # lj=0,1,2
            
            # my($tb,$te)= ($cdnafbeg+$trimb,$cdnafend); # use xtgb,xtge if have
            $xtgb+= $trimb;
            if($rev) { $$cdnaexoni[4] -= $trimb; } else { $$cdnaexoni[3] += $trimb; }
            $xfixi.="xtrimb$ix:$xtgb-$xtge/$trimb,";
            
          } elsif($xokspan >= $ENDCUT*$cleni) {  ## 0.75,0.66 what?
            $okcds=0; 
            my $trime= $cleni - ($xokspan + 3); #?  
            $trime=0 if($trime < 0);
            # $cdnai= substr($cdnai,0,$trime); # wrong** 
            $cdnai= substr($cdnai,0,$xokspan-3);  
            
            ## no phasing when trim at end?
            # my($llen,$lj,$lexon,$lcdna,$longorf)= phaseorf( substr($cdnai,0,$trime), $cdnafull);
            # $cdnai= $lexon; # no# $trime += $lj; # lj=0,1,2
            # $cdnafbeg += $lj; #?? want this? phaseorf() is cutting at start of cdnai only
            
            # my($tb,$te)= ($cdnafbeg,$cdnafend-$trime);
            $xtge -= $trime;
            if($rev) { $$cdnaexoni[3] += $trime;} else { $$cdnaexoni[4] -= $trime; } # rev problems! # change exon.stop
            $xfixi.="xtrime$ix:$xtgb-$xtge/$trime,";
            
          } else { # drop exon if error in middle?
            $okcds= -1; 
            # $cdnai=""; $cdnaexoni= undef;  
            # $xfixi.="xdrop$ix:$cdnafbeg-$cdnafend/$lend,";
          }
          $cleni= length($cdnai);
        } else {
          $okcds=1;
        }
                
        if($okcds>0) { $cdnafull= $workseq; }
        elsif($okcds<0) {
         # FIXME here? keep $cdnaexoni if xdrop cdnai but is utr part.. treat end drop-exons as utr seq/exon?
         ## end-drop becomes  not $overcds; $okcds=1; 
         if($okcds == -1 and ($ix == 0 or $ix == $nexon1)) { $okcds=1; $cdnafull= $workseq; $xfixi.="xdropc$ix:$xtgb-$xtge/$lend,";}#??
         else { $cdnai=""; $cdnaexoni= undef; $xfixi.="xdrop$ix:$xtgb-$xtge/$lend,"; }
         }
      }
      
      #......
      # split bug: xdrop exon on other part, leaving only mRNA of split=2 .. should convert to utr-exon not drop
      
      $xfixed.=$xfixi if($xfixi);
      if($cdnai and ref $cdnaexoni) {  
        ##? append xfixi to cdnaexoni[8] attr
        ##?? append xtrimi == phase shift?
        $$cdnaexoni[8].= ";cdsfix=$xfixi" if($xfixi); 
        push @cdnap, $cdnai ; # use cdnafull only?
        push @cdnaxp, $cdnaexoni; # exon feat to keep, maybe changed start/stop
      }
      
    }  # each exon

    ## this works only for some cases, eg fails if map cover is low.
    ## check and cancel if neworfaalen << oldaalen ? or need prior/nonfix aalen to see if improved enough to keep
    ## and odd cases improve beyond orig pubaa, presumably removing trasm errors or possibly got utrorf * BAD result
# Funhe2EKm036243t1	1980,28%,complete	3329,45%,complete-utrpoor	92  << nofix: cds2aa << pubaa
# Funhe2EKm036243t1 aalen=5292,58%,complete; << cdsfix improved too much??, possibly picked utrorf; this is DSCAM **
#  ID=Funhe2EKm036243t1;cov=92%,20249/22083;pid=98.6;nexon=75;splice=268;trg=Funhe2EKm036243t1 1 22083;
#  gaps=1529,MGap:12449-12561,MGap:12294-12435,MGap:6230-6395,MGap:6181-6203,MGap:12043-12282,MGap:20542-20783,MGap:20904-21410,MGap:20796-20891,;gescore=85;
#  clen=22083;offs=3714-13703;Name=Down syndrome cell adhesion molecule (54%M);
#  cxlen=15876/27197;aalen=5292,58%,complete;cdsoff=491-16369;utrx=3,41;
#  cdsfix=xtrimb25:6194-6474/35,xtrimb49:12987-14485/191,xtrimb49:12836-14294/40,xtrimb49:13102-14254/306,xtrimb50:12983-13586/186,xtrimb50:12873-13400/76,xtrimb50:12891-13324/94,;xtrim=1;
#  aaold=3329,45%,complete-utrpoor
    
#   * fixme genefindcds2.pl update -fixcds needs gap check/trim
#   gap err for Funhe2EKm017647t1
#    /protein_id="Funhe2GD:Funhe2EKm017647p1"  /translation="XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXLLEKAYAKLNG
# from kfish2n3gmap_test4bfcds5.gff
# Funhe2EKm017647t1 Target=Funhe2EKm017647t1 733 2479;cov=70.3;;offs=42-2402;
#  cxlen=1581/1583;aalen=527,99%,partial;protein=XXXXXXXXXXXXXXXXXX
#  cdsoff=3-1583;inshort=xcut12:7;cdsfix=xdrop17:2376-2479/1646,;xtrim=3;aaold=786,95%,complete
# Funhe2EKm017647t1 kfish2n3gmap_test4bfcds5.cdna
# >Funhe2EKm017647t1 aalen=527,99%,partial; cdsoff=3-1583; cxlen=1581/1583; qlen=2485; aaold=786,95%,complete; oid=Funhe2Exx11m006520t2,Funhe2E6bm006529t3;  len=1583;
# TNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
# NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCTCCTGGAGAAAGCTTACG
#   BUT getcdna() below calls this.  
#     my($didtrim,$texondna,$tft,$tstart,$tstop)= trimcdnaexon($j,$exondna,$ft,$rev,$nx1);
# * This bug because of dropped terminal exon .. dropped here? ix=19 missing, ix=18 has TNNN gap, but trim() skipped as it was 2nd, not 1st (rev)

    
    $addattr .= ";cdsfix=$xfixed" if($xfixed); 
    # update call cdnaseq w/ new exons; seems to work best
    if(@cdnaxp>0) {
      my $xtrim1;
      ($cdna,$cdnalen, $cdnaexonseq, $xfixfix, $xtrim1)= getcdna2( \@cdnaxp, 1, 0, 0, CDNA_TRIMNNN);
      $xtrimmed += $xtrim1; # xtrimmed used above
      $xtrimmed++ if($xfixed);  
    }
  
  } # DO_CDSFIX
  
use constant SPLITGENEPROTFIX => 1;  
unless(SPLITGENEPROTFIX) {
  if($cdnalen>0) { } # done DO_CDSFIX
  else { ($cdna,$cdnalen, $cdnaexonseq, $xfixfix, $xtrimmed)= getcdna2( \@exongff, 1, 0, 0, CDNA_TRIMNNN); }

} else { ## if(SPLITGENEPROTFIX)  
  if($cdnalen>0) { } # done DO_CDSFIX now
  elsif(not $issplit) { # not split, needs cdnalen
    ($cdna,$cdnalen, $cdnaexonseq, $xfixfix, $xtrimmed)= getcdna2( \@exongff, 1, 0, 0, CDNA_TRIMNNN);

  } elsif($issplit) {
    my ($lp,$lid,@xp,@cdnap,@cdnaxp); $lid=$lp=0;
    my @partexons;
    for my $i (0..$#exongff) {
      my $x= $exongff[$i];
      # my($aid) = ($x->[8] =~ m/(?:Parent|ID)=([^;\s]+)/)? $1 :0;
      # my($spla)= ($x->[8] =~ m/;Split=(\d+)/)? $1 : ($aid =~ m/_C(\d+)$/)? $1 : 0;
      my($spla,$aid)= splitPart($x->[8]);

      if($spla ne $lp) { # or $aid ne $lid
        if(@xp) { 
          my($cdnai,$cleni,$cdnaexonseq,$xfixfix, $xtrimmed)= getcdna2( \@xp, 1, 0, 0, CDNA_TRIMNNN);
          ## buggers, need $xfixfix, $xtrimmed here? NO EFFECT of this seen
          my $cdnaexoni= [@xp];
          $cdnaexoni= $xfixfix if($xtrimmed); ## NEED mRNA/gene span update also ****
          push @cdnap,$cdnai; push @cdnaxp,@$cdnaexoni; @xp=(); 
          push @partexons, $cdnaexoni;
          }  
        }
      push @xp,$x;  $lp=$spla; $lid=$aid; 
      }
    if(@xp) { 
      my($cdnai,$cleni,$cdnaexonseq,$xfixfix, $xtrimmed)= getcdna2( \@xp, 1, 0, 0, CDNA_TRIMNNN);
          ##  NO EFFECT of xtrimmed seen here
      my $cdnaexoni= [@xp];
      $cdnaexoni= $xfixfix if($xtrimmed); ## NEED mRNA/gene span update also ****
      push @cdnap,$cdnai;  push @cdnaxp,@$cdnaexoni;  @xp=(); 
      push @partexons, $cdnaexoni;
      }  
      
    # assume forward mrna, find orfs adding each part w/ phase NNN tests

    # if(1) .. phaseorf is good now
    my $cdnafull="";  
    for my $ip (0..$#cdnap) {
      my $cdp= $cdnap[$ip];
      # my $cdnaexoni= $cdnaxp[$ix]; # not right, these are all exons, cdp is *part* of many exons
      my $partexons= $partexons[$ip]; # @$partexons
      my($llen,$lj,$lexon,$lcdna,$longorf)= phaseorf( $cdp, $cdnafull); 
      $cdnafull= $lcdna;   
      if($lj) {  $addattr .= ";sphase=$ip:$lj";
        $xtrimmed++; #<< need to register exon change for mrna update; also? flag in exon attr?
        my $rev= ($partexons->[0]->[6] eq '-')?1:0;
        if($rev) { $partexons->[0]->[4] -= $lj; } else { $partexons->[0]->[3] += $lj; } 
      }
    }
    $cdna= $cdnafull; $cdnalen= length($cdnafull);
    $cdnaexonseq= \@cdnap; # is this right? need ?
    $xfixfix= \@cdnaxp; 
    # $xtrimmed =0;  # $xfixfix, $xtrimmed = what?
    #......
  } # issplit, SPLITGENEPROTFIX
} # SPLITGENEPROTFIX
  
  
    # DO_CDSCOMPLETE completeCDS() should follow SPLITGENEPROTFIX call?
    # 201511. add extend-ends option, ix==0 and ix==$nexon1, look for completing cdsbases of partial cds
    # can do that inside DO_CDSFIX? NO, do separate... with what flags/opts?
  if($DO_CDSCOMPLETE) {
    my $docomp=($mrnaat=~/cdsfix=complete/ or ($oldaalen=~/partial/ and $DO_CDSCOMPLETE>1))?1:0;  
    $docomp=0 if($issplit); # too many problems, cdsb,cdse on diff scafs
    my($extbest,$EXTcdna,$EXTexongff,$extaddattr,$ext5,$ext3)=(0) x 9;
    if($docomp) {
    
      my $cdnaexons= ($xtrimmed)? $xfixfix : \@exongff;
      my $oldprotPar= ($oldisbest)?$oldprot:"";
      
      ## 1707 FIXME: get computed proStart,Stop from  exons .. in completeCDSb()
      # $KEEPSAMECDS=1; ## cdna_proteins:getBestProt2(), need for oldStart_b
      # my($orfprotCOMP, $prostart5COMP, $proend3COMP,)= 
      #  getBestProt("partial", $cdna, $cdnaexons, $oldProstart5,$oldProend3);
      
      #** bad EXTcdna now, reextract from exongff... drop ,$EXTcdna return
      ($extbest,$EXTexongff,$extaddattr,$ext5,$ext3)=
        completeCDSb($geneid, $cdna, $cdnaexons, $revgene, $oldProstart5,$oldProend3,$oldprotPar);  
          # drop par: $mrnaat, $oldaalen,$oldoffs,$oldStart_b,$oldStart_e
          # add  par: $cdna,  ?? $oldprot
          
      if($extbest) {
        ## DANG, this is bad for splits, changing exon start/end to wrong split part in EXTexongff
        ## bad EXTcdna; $cdna= $EXTcdna; $cdnalen= length($cdna);
        ($cdna,$cdnalen)= getcdna2( $EXTexongff, 1, 0, 0, CDNA_TRIMNNN);
        @exongff= @$EXTexongff; # gosh darn too many same vars
        $xfixfix= $EXTexongff; # gosh darn too many same vars
        $xtrimmed++; # * REGISTER prior exon change for  FIXMEd, CDS extended NOT exons ..
        $addattr .= ";$extaddattr" if($extaddattr); 
      }
    }
  }
  
  if($xtrimmed) { 
    #o# @exongff= @$xfixfix; ## NEED mRNA/gene span update also ****
    #** FIXME: ^^^ cloned xfix exons NOT registering in exon.gff output;
    # .. all this dang to-fro of generecIN, @generec fwd/rev is losing changes ***
    
    my($igenerec, $iexongff, $mrnachanged1)= generecfix( \@generec, \@exongff, $xfixfix); 
    @generec= @$igenerec; 
    @exongff= @$iexongff; ## == xfixfix usually
    $mrnachanged += $mrnachanged1;    
    $changed ++; $addattr .=";" if($addattr); $addattr .= "xtrim=$xtrimmed"; 
    $oldStart_b=$oldStart_e=0; # dont use now for prot/cds
  }

  #ov2: $issplit : collect+paste cdna parts before getBestProt
  #.........
  
  # add samecds opt using oldStart > need start,stop of 1st CDS exon
  $KEEPSAMECDS=1; ## cdna_proteins:getBestProt2(), need for oldStart_b
  my($orfprot, $prostart5, $proend3, $bestorf, $utrorf)= 
    getBestProt("partial", $cdna, \@exongff, $oldStart_b,$oldStart_e);
  $KEEPSAMECDS=$samecds; ## cdna_proteins:getBestProt2(); reset
  ## buggers: NoStopCodon ignored here.. -nostop has '*' stop end
  ## $bestorf->{innerstop}  << SOME HaVE this fixme..
  
#... loop here for intron overlaps >>>>>>>>>>>>>>>>>>>> ov2: GONE..

  ## UPD: use these computed CDS, not cdnain CDS, as best match to chr locs
  my ($prostart5COMP, $proend3COMP)= ($prostart5, $proend3); # maybe just these needed?
  # my ($cdsgffCOMP, $cdsattrCOMP)= getCDSgff( \@exongff, $orfprot, $prostart5, $proend3, $cdnalen);
  # my ($cdsgffCOMP, $cdsattrCOMP)= getCDSgff2( \@exongff, newOrf(protein=>$orfprot, start=>$prostart5, stop=>$proend3, cdnalen=>$cdnalen) );

#... test input cdna ..............  
  ## FIXME: cdnain is strand-less ; rev for fixstrand ?
  ## FIXME2: test cdnain before / after intronfix, intronfix can be bad.
  
  my ($cdnain,$cdnainhdr)  = ("",""); 
  my ($cdnainlen, $cdnainNotPrimaryId, $cdnainIsBest)= (0) x 9;
  my $cdnagfflen= $cdnalen;
  if($cdnaseq) {  # add 2011
    $cdnain= $cdnaseq{$geneid}; $cdnainhdr= $cdnahead{$geneid};
    if( not $cdnain and ($geneid ne $geneidfix) ) { 
      $cdnain= $cdnaseq{$geneidfix}; $cdnainhdr= $cdnahead{$geneidfix}; 
      $cdnainNotPrimaryId=1 if($cdnain); }
   } 
   
  if($cdnain) { 
    $cdnainlen= length($cdnain);

# FIXed: some bug here preventing Much Better cdnain orf from replacing mapped orf
# FIXED: ** userev index(cdna) test was bad; check both orfs; and/or check mrna sense=-1 annot from gmap >> trrev of gff
#   is bestorf_test() ok? yes
#   .. in cases of < 50% mapped cdnain; could be as simple as userev test, may need to do both strands

# FIXME2: longer cdnain orf not always better.. as usual.  Some of these are trasm mashups of 2+ genes,
#   mapped transcript is chimera at same locus giving clue to mashup.  happens more for velvet low kmer
#   where smaller fragments of 2+ genes are mashed together.   Need to check mRNA chimera tags for clues?

## Dang BUG2: chimera-only? cdnabest prot: protein=*MTGHYYYE..HFEY*  ^* is bogus stopcodon, from where ??
##   also mayb only with aautr ; also -nostopcodon removes * prefix, but doesnt remove all end *
## FIXED below at $oldprot.'*'; 

# FIXME3: chimera, _C1 and _C2 getting same protein for cdnainIsBest, including bad cases of C1=95%, C2 = 5% mapped
#  .. here or later allow only 1 of pair to have cdnainIsBest.. should be based on mapped CDS best match to protein
#     not just longest mapping portion (which often is bogus gene join as UTR).

## FIXME4: strand reversal, $userev?, of cdnainIsBest vs gff should not be allowed sometimes.. 
#      .. or need to reverse strand of output gff
#   i.e. changing  $prostart5, $proend3, to cdnain set is not valid for rev strand, w/o flip of exon strand
#   see fixstrand..

## FIXME5: new bug: Now gap-cdna with XXX protein winning over cdnain with full protein, same length. 
#   adjust goodlen test? problem may be index(cdnain,cdna) fails due to XXXX garbage

# new # need test both ** YES, this cures problems
# ** BUT need to know if best cdnain is same/rev of gff cdna .. FIXME4
# .. revert to old, bad; need cdnain strand same as input gff/genome cdna for fully mapped cdna
# .. but test both strands for split-maps (chimera), partial maps (coverage < 60?), and oneexon trs.

    use constant CIN_SAMESTRANDasMAP => 1;
    use constant CIN_BOTHSTRANDS => 0;
    
    use constant USE_EVG_MRNA_HEADER => 1;  #1707 upd
    my $UseMRNAhdr= (USE_EVG_MRNA_HEADER && ($cdnainhdr && $cdnainhdr =~ m/(?:cdsoff|offs)=/)) ? 1:0;

    my $cdnain_isrev=0; # doesnt mean samestrand-as-map, but just that cdnain was reversed.
    my $cstrandIsOpposite= 0; # 1 if antisense opposite, 0 == sense same as genome map; add -1 == dont know/cant align?
       ## FIXME6: bug in cstrandIsOpposite, reporting same Sense strand when opposite.
       ## .. probable failure of $cdnaindex not finding cdnagff inside cdnain.
       
    my($orfproti, $prostart5i, $proend3i,$bestorfi, $utrorfi);
    my $cdnaindex= -1;
    # my $cdnainSense=1; # or -1 if opposite of mapped strand
    # from gmap, "sense=-1" means cdnain orient was opposite of mapped strand
    # ?? what of utrorf on rev strand of bestorf ??
    
    ORFCDNAIN: { 
    my $cdnainrev= revcomp($cdnain);
    # replace exon count with valid intron count?
    # my $validgffstrand= ( $inok > $inerr) ? 1 : ($inok+$inerr == 0 and @exongff > 1) ? 1 : 0;
    my $validgffstrand= ( @exongff > 1 ) ? 1 : 0;
    # cancel if (abs($inerr) > $inok) 

    ## FIX: 160530.  # gmap wrong for mRNA mapping; fixme 1605
    my $orfaalen= length($orfprot);
    my ($oldaalend)= ($oldaalen=~m/(\d+)/)?$1:0;
    my $badgfforf= ($mrnaat =~ m/;sense=-1/ or ($oldaalend * 0.80 > $orfaalen))?1:0;
    if($badgfforf) { $validgffstrand=0; }

    # problem here, some chimera cdna parts can index cdnain both ways; ie cdnain is mashup of 1, or 2 near same.
    # .. Or, best of both is reverse of this cdna part, but on other cdna part.. case of rev-genes-join
    # .. cdna part does not map to bestorf part, so cant correct cdsgff from this
    
# ** FIXME here dang gmap.gff sense=-1 strand is WRONG dont use, validgffstrand=0 ; use oldaalen >> new; testboth

    my $cstrand_gstrand= 0;
    #drop# my $cdnagffgoodlen= $bestorf->{goodlen}; # not cdnagfflen; dang this is cds not cdna len
    my $cdnagffgoodlen= $cdnagfflen;
    my $cdnagood= $cdna; # only for index cstrand_gstrand test
    if(index($cdnagood,'N') >= 0) { ## No ($cdnagffgoodlen < $cdnagfflen)
      $cdnagood =~ s/^N+//; $cdnagood =~ s/N+$//; 
      $cdnagffgoodlen=length($cdnagood);
      my $xi= index($cdnagood,"NN"); 
      if($xi>0) { 
        my $xe= rindex($cdnagood,"NN");
        if($xi>$cdnagffgoodlen/2) { $cdnagood= substr($cdnagood,$xe+2); }
        else { $cdnagood= substr($cdnagood,0,$xi); }
      }
    } 
    
    if( ($cdnaindex= index($cdnain,$cdnagood)) >=0 ) {  $cstrand_gstrand= 1; }
    elsif( ($cdnaindex= index($cdnainrev,$cdnagood)) >=0 ) { $cstrand_gstrand= -1; }
    if($cdnaindex < 0) { 
      ## FIXME6: try harder to match strands; indel problems
      ## see below also, per exon-gff test for rev
      my ($xfwd,$xrev)=(0,0);
      foreach my $xs (@$cdnaexonseq) {
        if(index($cdnain,$xs) >=0 ) { $xfwd += length($xs); }
        elsif(index($cdnainrev,$xs) >=0 ) { $xrev  += length($xs); }
      }
      $cstrand_gstrand= ($xfwd>$xrev)? 1 : ($xrev>$xfwd)? -1 : 0;
    }
    
    ## here use goodlength tests for cdnain,cdnagff
    ## testboth is problem now: getting many Anti of gene-mashup (fwd+rev), want only cases of gapfill tr?
    my $testboth= 0; #  2=best of both 
    if($fixstrand) { } # DONT test both in this strand loop
    elsif( $validgffstrand and $cdnainlen < 1.75 * $cdnagffgoodlen) { $testboth=0; }
    elsif( $cdnainNotPrimaryId or ($cdnainlen >= 1.75 * $cdnagffgoodlen)) { $testboth= 2; } # or not $validgffstrand
    elsif( $cstrand_gstrand == 0 ) { $testboth= 2; } #? or leave 0
    
    # gmap wrong for mRNA mapping; fixme 1605
    if($badgfforf) { $testboth=2; }

    if(CIN_SAMESTRANDasMAP) { }   
    elsif(CIN_BOTHSTRANDS) { $testboth=2; }

    ## cdnainhdr from Evigene has protein offset on mRNA.seq .. use that instead of recalc orfs ??
    my ($longorf,$longfull,$orfs);
    
    if( $UseMRNAhdr ) {
      # UGH offsets are bad here, rev? lots of instops .. getOneOrf zero-origin
      my($cdsoff)= ($cdnainhdr =~ m/(?:cdsoff|offs)=([^;\s]+)/) ? $1 : "";  
      my $cdsor=0;
      my($cdsb,$cdse)= ($cdsoff and $cdsoff=~m/(\d+).(\d+)/) ? ($1,$2) : (0,0);
      if($cdse>0 and $cdse < $cdsb) {  $cdsor= -1; } #  ($cdsb,$cdse)= ($cdse,$cdsb);: DONT swap for OneOrf
      elsif($cdsoff =~ m/(\d+).(\d+):(.)/) { my $o=$3; $cdsor=($o eq '-')?-1:1; } # :. is common, same as :+
      else { $cdsor= 1; }
      if($cdsb>0 and $cdse>0) {
        #?  $cstrandIsOpposite=1 if( $cstrand_gstrand ne $cdsor ); #?? bad likely
 	$cstrandIsOpposite=0; 
        $cdnain= $cdnainrev if($cdsor<0); #? here?
        $cdnain_isrev= ($cdsor<0)?1:0;
        my $orf= getOneOrf( $cdsb-1, $cdse-1, $cdnain, '+'); # ** ARG, zero-origin cdsb,e here ***
        $longorf= $longfull= $orf; $orfs= [$orf]; # ($longest, $longfull, \@orfs) = getAllOrfs()
      } else {
        $UseMRNAhdr= 0; 
      }
    }
    
    if($UseMRNAhdr and ref($longorf)) { } # done
    elsif( $testboth == 2 ) {    
      ($longorf,$longfull,$orfs) = getAllOrfs($cdnain,"fwd"); # ,"dropnnn"
      my ($rlongorf,$rlongfull,$rorfs) = getAllOrfs($cdnainrev,"fwd"); # ,"dropnnn"
      ($cdnain_isrev)= bestorf_test($longorf,$rlongorf);
      if($cdnain_isrev) {
        # $cstrand_gstrand= -1 if($cstrand_gstrand == 0); 
        ## report if cstrand_gstrand unknown: cstrandIsOpposite=-1 ??
        $cstrandIsOpposite=1 if($cstrand_gstrand == 1);
        $cdnain= $cdnainrev;
        ($longorf,$longfull,$orfs)= ($rlongorf,$rlongfull,$rorfs);
      } else {
        # $cstrand_gstrand=  1 if($cstrand_gstrand == 0);
        $cstrandIsOpposite=1 if($cstrand_gstrand == -1);
      }
    } elsif( $cstrand_gstrand == -1 ) {
      $cdnain_isrev=1;
      $cdnain= $cdnainrev;
      ($longorf,$longfull,$orfs) = getAllOrfs($cdnain,"fwd"); # ,"dropnnn"
    } else {  ##  $cstrand_gstrand == 1
      $cdnain_isrev=0;
      ($longorf,$longfull,$orfs) = getAllOrfs($cdnain,"fwd"); # ,"dropnnn" 
    }
    
    ($orfproti, $prostart5i, $proend3i, $bestorfi, $utrorfi)= 
      getBestProtOfOrfs($longorf,$longfull,$orfs, 
        "partial", $cdnain, \@exongff, $oldStart_b,$oldStart_e);
    }
    

    if($orfproti ne $orfprot) { 
      ($cdnainIsBest,undef,undef) = bestorf_test($bestorf,$bestorfi);
       # FIXME bestorf_test: option too high -ratiocdna 1.25; at least should replace cdna XXXX gaps with cdnain perfect
       # test for near-sameprot but for XXX gaps? 
       ## if($orfprot =~ /XXX/) { (my $op=$orfprot) =~ s/X/./g; $cdnainIsBest=1 if($orfproti =~ /$op/); }
      
      # FIXME1707: cdnain may have been intentionally cut (utrx) in gff, but not in mrna.seq
      # .. how to fix that? need flags or table?
      
      # FIXME2: longer cdnain orf not always better.. as usual.  Some of these are trasm mashups of 2+ genes,
      # add other checks here if cdnain better than genome-mapped cdna .. esp if cdnain == chimeric mapping
      #  versus partly unmapped due to genome gap.
      
      ## reasons cdnain is best: gaps in genome mapping, NNNs or indels, versus cleaner transcript
      ## reasons not best: chimeric mapping over same locus or nearby; reversed introns over longer cdnain mapping, ..
      ## .. maybe not .. leave this off for now.
#             
#       # my $cdnabestExplained= bestorf->goodlen/length < cdnainorf->goodlen/len
      #ov2: my $cdnainProblems = $cdnainNotPrimaryId or ($hasintrons and $inerr != 0) or ($mrnaat =~ /chimera/);
      my $cdnainProblems = $cdnainNotPrimaryId or ($mrnaat =~ /chimera/);
      
      if($cdnainIsBest and $cdnainProblems) { 
      
        if($cdnainNotPrimaryId) { ## and $mrna->[8] =~ /chim\d=(\w+[^;\s]+)/
          # cancel if $pctmapped < 25, < 33? <50?
          my $ingood= $bestorfi->{goodlen} || 1;
          my $bgood = $bestorf->{goodlen};
          
          $cdnainIsBest= 0 if($bgood/$ingood < 0.25); # bad test.. really need other part's stats here
          ## and/or cancel if mapped CDS are missing or small fraction of cdnain?? need from getCDSgff() or @oldCDS?
          
        }
        
#         #?? this will cancel some good changes along w/ bad; add annot for scoring instead?
#         if($hasintrons and $inerr != 0) { $cdnainIsBest= 0; }
#         elsif($mrnaat =~ /chim\d=(\w+[^;\s]+)/) { 
#           # pick out chimera span and see if near this one.. chim1=scaffold00024:99418-99575:.
#           use constant CHIOK_MIN => 25000;
#           my $cloc= $1; my($cr,$cb,$ce)= $cloc =~ m/(\w+):(\d+).(\d+)/;  
#           my($thisr,$thisb,$thise)= ($mrna->[0], $mrna->[3], $mrna->[4]);
#           # .. maybe not always cancel; add flag instead? ..
#           $cdnainIsBest= 0 if($thisr eq $cr and (abs($thisb - $ce) < CHIOK_MIN or abs($thise - $cb) < CHIOK_MIN));
#         }
        # elsif($cdnainNotPrimaryId) {} #? what check?
      }
      
    }
      

    if($debug) {  #  and not $cdnainIsBest; debug set note; add orig cxlen;aalen here if cdnainIsBest?
      my $trd = ( $cstrandIsOpposite ) ? "trAnti": "trSense";
      $trd .= ($cdnain_isrev) ? "Rev" : "Fwd"; # less useful than anti/sense
      if(1) { # new
        my $uqual= proteinqual($bestorfi); #in cdna_proteins.pm
        $uqual=~s/=/:/g; $uqual=~s/;/,/g;
        $addattr .=";" if($addattr); 
        $addattr .= "cdnaorf=$uqual,trbest:$cdnainIsBest,$trd";  
      } else { # old
        my( $uprot, $ustart, $uend) = orfParts($bestorfi);
        my $ulen=  $bestorfi->{length};  
        my $ucdnalen=  $bestorfi->{cdnalen};  
        my $ualen= length($uprot); $ualen-- if(substr($uprot,-1) eq '*');
        my $upcds  = ($ucdnalen>0 && $ulen>0) ? int(100*$ulen/$ucdnalen) : 0; # was cdnalen != cdnainlen here; 
        my $ucompl= $bestorfi->{complete};
        $ucompl= ($ucompl==3)?"complete":($ucompl==2)?"partial5":($ucompl==1)?"partial3":"partial";
        $addattr .=";" if($addattr); 
        $addattr .= "cdnaorf=$ualen,$upcds%,$ucompl,trbest:$cdnainIsBest,$trd"; #?  
      }  
    } 
    
    # UPD1707: also test mrna-cds seq align to chr-cds, not cdnainIsBest (as option?)
    # .. need to know if chr-cds is covering at least parts of same cds as mrna .. if not, diff gene
    
    if($cdnainIsBest) {
      # replace?, notice, may need to change CDS like intronfix : depends on diff
          
      my $aadif= length($orfproti) - length($orfprot); 
      my $gadif= $bestorfi->{goodlen} - $bestorf->{goodlen}; 
      # my $trdif= length($cdnain) - length($cdna); 
      my $trdif= $cdnainlen - $cdnagfflen;  
      $cdnalen= $cdnainlen; # change it
      
      my $trdif1=$trdif;
      map{ $_="+".$_ if($_>0); $_= (($_==0)?"eq:":"NE:").$_; } ($aadif,$gadif,$trdif);
      my $aac= $bestorfi->{complete} - $bestorf->{complete}; # complete == 3,2,1,0
      
      # also compare cdna vs cdnain : size, where diff (inside, ends?)
      # mrna-attr: Coverage=40.0;Identity=99.7 < use this to see diff?
      my $nx0= @exongff; # no way to count exons in cdnain here.
      my $trin=0; my $nxeq=0;
      
      ## change cdnain_isrev to cstrandIsOpposite here?
      # my $trd= ($cdnain_isrev) ? "trREV" : "tr";
      my $trd= ($cstrandIsOpposite) ? "trdAnti" : "trdSense"; # "trSens" ?
      $trd .= ($cdnain_isrev) ? "Rev" : "Fwd"; # less useful than anti/sense
      
      ## FIXME4: strand reversal, $cdnain_isrev?, of cdnainIsBest vs gff should not be allowed sometimes.. 
      #      .. or need to reverse strand of output gff
      #   i.e. changing  $prostart5, $proend3, to cdnain set is not valid for rev strand, w/o flip of exon strand
      #   see fixstrand..
      
      ## fixme, use cstrandIsOpposite; also $prostart5, $proend3, are wrong or exongff needs strand change
      ## BUG, dont change this part strand when it is other chimera part that has bestorf w/ other strand..
      ## drop this for now, not right for cases seen.; use trdAnti note ; other flag for problems?      
#       if($cstrandIsOpposite and not $fixstrand) {
#         if($gstrand eq "-") { $gstrand="+"; } elsif($gstrand eq "+") { $gstrand="-"; }
#         map{ $_->[6]= $gstrand } @generec;
#         ($mrna) = grep{ $_->[2] eq "mRNA" } @generec;
#         @exongff= grep{ $_->[2] eq "exon" } @generec;
#         ## .. other problems w/ restrand..  oldCDS ? 
#       }
  

      if($cdna eq $cdnain) { $trd .= "eq"; $nxeq=$nx0; }
      elsif( $cdnaindex >=0 ) { # from index($cdnain|cdnainrev, $cdna)
        $trd .= "$trdif,in$trin";
        $nxeq= ($nx0==1)?1:$nx0-1;
        if($cstrandIsOpposite) {}
      }
      # elsif( ($trin=index($cdnain, $cdna)) >=0 ) { $trd .= "in$trdif,$trin"; $nxeq= ($nx0==1)?1:$nx0-1; }  #nxeq maybe; calc?
      # elsif( ($trin=index($cdnain, revcomp($cdna))) >=0 ) { $trd .="rc$trdif,$trin"; }
      else { 
        $trd .= $trdif; # * check each \@exongff > exseq index cdnain ? report if ends or inner x diff
        # my($cdna2, $cdnalen2, $cdnaexonseq)= getcdna( \@exongff, 1); # use from above cdnaexonseq
        my ($ix,$le,$xfwd,$xrev)=(0,1,0,0); 
        foreach my $xs (@$cdnaexonseq) {
          $ix++; my $xi= index($cdnain, $xs);  my $xr=""; 
          if($xi>=0) { $xfwd+= length($xs); }
          else { $xi= index($cdnain, revcomp($xs)); if($xi>=0) { $xr="c"; $xrev+=length($xs); }  }
          if($xi>=0) { my $xe= length($xs)+$xi; my $xb=$xi+1; $trd.=",xeq$ix:$xr$xb-$xe"; $nxeq++; 
             if( $le>1 and (my $g= $xb - $le) > 0 ) { $trd.= "g$g"; } $le=1+$xe; }
          else { $trd.=",xne$ix"; }
          my $xgap= $xs =~ tr/N/N/;  $trd .="N$xgap" if($xgap>0);
        }
        if($xrev > $xfwd and not $cstrandIsOpposite) { $trd =~ s/trdSense/trdAntix/; } # ERROR??
      }

      #check further if $nxeq > 0; xinner change?  count NNN gaps in cdna, cdnain; look for gaps at end if inner match
      if($nxeq > 0) {
        my $gaps = $cdna =~ tr/N/N/;
        my $gapi = $cdnain =~ tr/N/N/;
        if($trdif1>25) { 
          my($cdna2, $cdnalen2)= getcdna( \@exongff, 0, 100, _max(100,$trdif1));
          $gaps = $cdna2 =~ tr/N/N/;
        }
        if($gaps > $gapi) { $trd.= ",tN:$gaps"; }
      }       
      
      my $cov= ($mrnaat =~ m/(cov|Coverage)=(\d+)/) ? $2 : 0;
      $trd .= ",tcov:$cov" if($cov); # cov < 99 explains
      my $changenote= "cdnabest=aa$aadif,gd$gadif,dfull:$aac,$trd,nx0:$nx0";
      $addattr .=";" if($addattr); $addattr .= $changenote; $changed++;
      
      # FIXME HERE >> want CDS exons computed before cdnain . this changes to oldCDS prostart/end
      # .. KEEP orfprot= orfproti, NOT pro5,3, ??not utrorf
      if($USE_CHRMAP_CDS) {
      ($orfprot)= ($orfproti); # $prostart5,3 == $prostart5COMP,3 here
      } else {
      ($orfprot, $prostart5, $proend3, $utrorf)= ($orfproti, $prostart5i, $proend3i, $utrorfi);
      }
    }
  }

 
  if($orfprot) {
    my($annew,$anold)= stripOldAnnot($mrna) if($REANNOTATE); # also changes mrna

    my ($cdsgff, $cdsattr);    
    if($USE_CHRMAP_CDS) {
      # UPD: computed CDS not cdnainIsBest CDS
      # ($cdsgff, $cdsattr)= ($cdsgffCOMP, $cdsattrCOMP); # bad cdsattr for cdnainIsBest
      ## cdnalen == cdnainlen for cdnainIsBest; mismatch w/ cdnalenCOMP problems?  for using start/stop COMP
      ($cdsgff, $cdsattr)= getCDSgff2( \@exongff, newOrf(protein=>$orfprot, start=>$prostart5COMP, stop=>$proend3COMP, cdnalen=>$cdnalen) );
    } else {  
      ($cdsgff, $cdsattr)= getCDSgff2( \@exongff, newOrf(protein=>$orfprot, start=>$prostart5, stop=>$proend3, cdnalen=>$cdnalen) );
      # was getCDSgff(\@exongff, $orfprot, $prostart5, $proend3, $cdnalen);
    }
        
    # test $prostart5, $proend3 vs old
    
    my $diff=0; 
    $diff=1 if($REANNOTATE or $cdnainIsBest);
    
    my $diffcancel=0; my $oldstat="";  
    if(ref $cdsgff and @$cdsgff > 0) {
      my @newCDS; # = sort _sortgene @$cdsgff; # do here.
      if($issplit) {
      @newCDS= sort _sortSplitgene @$cdsgff; # do here.
      } else {
      @newCDS= sort _sortgene @$cdsgff; # do here.      
      }

      my @oldsave=();
      ## reversed CDS look wrong at end exons (always both ends? for nc>2)    
      if(@oldCDS > 0) {
        my($oldal, $newal, $oldprostart5, $oldproend3)= (0,0,0,0);
        $oldprostart5= $oldProstart5; $oldproend3= $oldProend3; # above, dont need both
        # compare, report in $addattr
        if($issplit) {
        @oldCDS= sort _sortSplitgene @oldCDS; # dang2    
        } else {
        @oldCDS= sort _sortgene @oldCDS; # dang
        }
        
        my $newprot= ($cdsattr =~ m/protein=([^;\s]+)/) ? $1 : "";
        if (@oldCDS == @newCDS) { $oldstat="ocds=eqn"; } else {  $oldstat="ocds=NEn"; $diff++; }
        
        if($oldprot) { 
          (my $op=$oldprot) =~ s/\*$//; (my $np=$newprot) =~ s/\*$//; 
          $oldal= length($op);
          $newal= length($np);
          my $oldap= ($cdnalen>0) ? int(0.5 + 300 * $oldal / $cdnalen) : 0;

          my $eq=($op eq $np)?1:0; # count STOP also
          if($eq and ($newprot =~ /\*$/) and not ($oldprot =~ /\*$/)) { $eq=0; }
          
          $diff++ unless($eq); 
          my $da= $newal - $oldal; $da="+$da" if($da>0);
          $oldstat .= ($eq) ? ",eqaa" : ",NEaa:$da"; 
          
          my $astat=0;
          if($oldprot =~ /^M/) { $astat |= 1; }
          if($oldprot =~ /\*$/) { $astat |= 2;} # ** AUG proteins lack '*' unless pasa-updated
          elsif($geneid =~ /AUG/) { $astat |= 2; } #  also check for inner X == augustus-fake for stop ?

          my $istop= index($oldprot,'*');
          ## find single 'X' inside prot, but allow this for NNN genome: SFXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXQP
          if($geneid =~ /AUG/ and (my $ix=index($oldprot,'X')) >0 ) { 
            if(substr($oldprot,$ix,3) eq "XXX") { } # ignore 2+; really should look at genome dna to decide
            else { $istop= $ix if($istop < 1 or $ix<$istop); }
            }
          # if($istop < 1 and $geneid =~ /AUG/) { $istop= index($oldprot,'X'); }

          # my $prostat = ($p5 and $p3) ? "partial": ($p5)? "partial5" : ($p3)? "partial3" :"complete";
          if($istop > 0 and $istop < length($oldprot)-1) { $astat= "partialinner"; } # innerstop ?
          elsif($astat == 3) { $astat="complete"; }
          elsif($astat == 1) { $astat="partial3"; }
          elsif($astat == 2) { $astat="partial5"; }
          else { $astat="partial"; }
          $astat =~ s/^(.)/o/; # dont confuse w/ new in scans
          # change oaaln= to aaold=  for consistancy
          $oldstat .= ";aaold=$oldal,$oldap%,$astat";
          #o $oldstat .= ";oaaln=$oldal,$oldap%,$astat";
          
        } else { $diff++; }

   # FIXME2: exist CDS : allow for problem cases like end-of-scaffold partials
        my $partialAtEndOk= 0;
        for(my $i=0; $i<@oldCDS; $i++) { 
          my $oc= $oldCDS[$i];
          $partialAtEndOk=1 if($oc->[3] < 450); # dont know high end to check
          unless($i<@newCDS and $oc->[3] == $newCDS[$i]->[3] and $oc->[4] == $newCDS[$i]->[4]) { 
            $diff++; $oldstat .= ",NE$i";  ## BUG: clone oc
            if($debug>1) { my @doc= @$oc; $doc[2]="oldCDS"; $doc[0] =~ s/^/#/; push @oldsave, \@doc; }
            }
          }
          
       
        do{ $diff=0; $diffcancel=1; } if(!$NODIFCAN and $partialAtEndOk and $oldstat =~ /artial/ and $newal < $oldal);

   # FIXME4: cancel changes that replace all CDS with completely new CDS, eg for transposon-spans
        $oldprostart5= ($gstrand eq "-") ? $oldCDS[-1]->[4] : $oldCDS[0]->[3];  #above# 
        $oldproend3  = ($gstrand eq "-") ? $oldCDS[0]->[3] : $oldCDS[-1]->[4];  #above# 
        $diff++ if($diff==0 and ($oldproend3 != $proend3 or $oldprostart5 !=  $prostart5));   # off by -1,-2 bugs   

        #  oldc:   [-----]
        #  newc:           [-----]  ; bad = shifted too much
        #  newc:        [-----]     ; bad "
        #  newc:     [---]          ; ok  = not shifted, same stop
        
        if($diff and $oldal>0) { # and $oldstat !~ /artialinner/
          ## .. bugs here ??
          #my($tb,$te)=($prostart5,$proend3); ($tb,$te)= ($te,$tb) if($tb>$te);
          #my($lb,$le)=($oldprostart5,$oldproend3); ($lb,$le)= ($le,$lb) if($lb>$le);
          my($tb,$te)=($newCDS[0]->[3],$newCDS[-1]->[4]); ($tb,$te)= ($te,$tb) if($tb>$te);
          my($lb,$le)=($oldCDS[0]->[3],$oldCDS[-1]->[4]); ($lb,$le)= ($le,$lb) if($lb>$le);
          if($lb==0 or $le==0 or $tb==0 or $te==0) { $tb=$te=$lb=$le=0; } # bug somewhere ...
          
          # argg; location stats not best here, introns have big effect.
          # .. this cancels most span increases...
          if(($tb < $lb and $te < $le) or ($tb > $lb and $te > $le)) { # CDS shifted along exons
            my ($bb,$be)= ( _max($tb,$lb), _min($te,$le) );
            my $maxo= 1 + $be - $bb; # ** not abs
            my $leno= _min( abs($le - $lb), abs($te - $tb)) || 1;
            my $pover= $maxo/$leno; # neg desired: cancel any case of maxo < 0 
            #??do{ $diff=0; $diffcancel=1; } if($pover < 0.33); # cancel change, too little cds similar
            do{ $diff=0; $diffcancel=1; } if(!$NODIFCAN and $pover < 0); # cancel change, too little cds similar
          }
        }

        $cdsattr .=";$oldstat";
        $changed++ if($diff>0);
      } else {
        $changed++; $diff=1;
      }

      if($diff) {
        push @generec, @newCDS; 
        push @generec, @oldsave; 
      } else {
        push @generec, @oldCDS;  
      }
    }
    # else { $diff++; } ## is this error? cdnainIsBest = cdnain best?
    #? else { push @generec, @newCDS;  } # no @$cdsgff no oldCDS
  
    ## add to addattr: count utrs (and span?) and add utrx= if >2
    # change getCDSgff() to return utr-exons?
    # if($u5>2 or $u3>2) { $un=$u5+$u3; $g[0]=~s/$/;utrx=$u5,$u3/; } 

    ## fixme: diff==0, dont remove old attr, add new or not? 
    #     or make option? sometimes always replace? use ocds=eqn,eqaa as flag for no change
    # .. except for some changes above result in NEaa but are cancelled.
    
    if($utrorf) {
      my( $uprot, $ustart, $uend) = orfParts($utrorf);
      my $ulen=  $utrorf->{length};  
      my $ualen= length($uprot);  
      my $upcds  = ($cdnalen>0 && $ulen>0) ? int(100*$ulen/$cdnalen) : 0;
      my $ucompl= $utrorf->{complete};
      $ucompl= ($ucompl==3)?"complete":($ucompl==2)?"partial5":($ucompl==1)?"partial3":"partial";
      $addattr .=";" if($addattr); 
      $addattr .= "aautrlen=$ualen,$upcds%,$ucompl;utroff=$ustart-$uend;utrprot=$uprot";  
    }
    
    if($cdsattr and $diff==0) {
      # check prot, add '*' if missing and complete
      # remove prot, new aalen if  $diffcancel
      unless( $mrnaat =~ m/;protein=/ or $diffcancel) {  ## should set diff=1 for this
        $diff=1; # $mrna->[8] =~ s/$/;$cdsattr/; 
        } ## add if missing
      if($diffcancel) { 
        my ($st)= $oldstat =~ m/ocds=([^;\s]+)/; 
        my ($nal)= $cdsattr =~ m/aalen=([^;\s]+)/;
        $mrna->[8] =~ s/$/;ocds=cancel:$st,nal:$nal/; 
        
      } elsif(not $NoStopCodon and $cdsattr =~ /complete/ and $oldprot =~ /\w/ and $oldprot !~ /\*$/) {
        my $ops= $oldprot.'*';  # this WAS bug for protein=*Maaaa* 
        $mrna->[8] =~ s/;protein=$oldprot/;protein=$ops/; # no end ;
      }
      
    } # elsif()
    
    my($didoldaa,$didoldoff)=(0,0);
    if($cdsattr and $diff>0) {  
        ## cdsattr keys: cxlen,aalen,protein,utrx
      my @keys= $cdsattr =~ m/(\w+)=/g;
      my $keyp= join('|',@keys);
      #x $mrna->[8] =~ s/[;]?($keyp)=([^;\n]+)//g; # OR/opt rename keyp_old ?
      $mrna->[8] =~ s/;($keyp)=[^;\n]+//g;  
      $mrna->[8] =~ s/$/;$cdsattr/;
      $didoldaa= ($cdsattr =~ /(aaold|oaaln)=/);
      $didoldoff= ($cdsattr =~ /cdsoff=$oldoffs/);
    }
    
    #o $addattr.=";aaold=$oldaalen" if($oldaalen and not($cdsattr =~ /oaaln=/));  
    #            #^^ aaold= duplicates oaaln=.. in $oldstat ; change to oaaln= ?
    #  $astat =~ s/^(.)/o/; # dont confuse w/ new in scans
    #  $addattr .= ";oaaln=$oldal,$oldap%,$astat";
    if($oldaalen and not($didoldaa)) { 
      $oldaalen=~s/,[cp]/,o/; $addattr.=";aaold=$oldaalen";  
    }
    $addattr.=";offold=$oldoffs" if($oldoffs and not($didoldoff)); ## add ??
    $addattr =~ s/;;/;/g; $addattr =~ s/,,/,/g; # $addattr=~s/;.;/;/; # bad compress patt? Name=Timeless;n;;
    
    $mrna->[8] =~ s/$/;$addattr/ if($addattr);
    $mrnaat= $mrna->[8]; # update
    
  } else { # no orfprot : can happen, but error if oldCDS << need to restore those here ??
    push @generec, @oldCDS if(@oldCDS);
    my $addattr="ocds=cancel:NEn,aatrans-failed";
    $mrna->[8] =~ s/$/;$addattr/;
    $mrnaat= $mrna->[8]; # update
  }

  if($ifix == 1) {  @generev= @generec; $changed1=$changed; }  # already cloned, now preserve
  else {  @genefwd= @generec;  $changed0=$changed; }
  
  ## output seq vars outside dang fixstrand loop
  $cdnaoutg= $cdna;  $orfprotg= $orfprot; $orfcdsg= $bestorf->{sequence}||"nnn";
  
  } # BIG loop for fixstrand 

  ## somewhere add nostopcodon opt:
  #  if($NoStopCodon and substr($orfprot,-1) eq '*') { $orfprot =~ s/\*$//; }
  
  my($generec,$revbest)= (\@generec,0);
  if($fixstrand) {
    ($generec,$revbest)= fixstrand( \@genefwd, \@generev);
    $changed= ($revbest) ? $changed1 : $changed0;
  }
  
  my $gflags=($mrnachanged)?"mrnachanged=$mrnachanged":"";
  # FIXME15: add opt output -cdna -cds sequences, separate files (also -aa seq sep file instead of mrna.gff attr)

  ## FIXME:  sort gene/splitgene @generec before putgene() ??
  if($issplit) {
    @generec= sort _sortSplitgene  @$generec; $generec= \@generec;
  } else {
    @generec= sort _sortgene  @$generec; $generec= \@generec;
  }

  putgene( $generec, $geneother, $gflags);
  
  my $seqattr= seqAnnot($mrna); #? add aaold
  putseq($houtcdna,$geneid,$cdnaoutg,$seqattr) if($houtcdna);
  putseq($houtaa,$geneid,$orfprotg,$seqattr) if($houtaa);
  putseq($houtcds,$geneid,$orfcdsg,$seqattr) if($houtcds);
  
  return ($changed) ? 1 : 0;
}
# end sub testgene() too long ......


sub fixstrand {
  my($genefwd, $generev)= @_;
  my $generec= $genefwd;
  my $revbest=0;
  
  my ($mrna)  = grep{ $_->[2] eq "mRNA" } @$genefwd;
  my $aafwd   = ($mrna->[8] =~ m/aalen=([^;\s]+)/) ? $1 : "";
  my $protfwd = ($mrna->[8] =~ m/protein=([^;\s]+)/) ? $1 : "";
  ($mrna)     = grep{ $_->[2] eq "mRNA" } @$generev;
  my $aarev   = ($mrna->[8] =~ m/aalen=([^;\s]+)/) ? $1 : "";
  my $protrev = ($mrna->[8] =~ m/protein=([^;\s]+)/) ? $1 : "";

  my $dlenrev= length($protrev) - length($protfwd);
  my $plenrev= abs($dlenrev) / _max(1, _max(length($protrev) , length($protfwd)));
  my $acfwd= ($aafwd =~ /complete/)?1:0; 
  my $acrev= ($aarev =~ /complete/)?1:0; 
  if($acfwd ne $acrev and $plenrev < 0.33)  {
    $revbest= ($acrev) ? 1 : 0;
  } else {
    $revbest=($dlenrev > 0)?1:0;
  }
  $generec= ($revbest) ? $generev : $genefwd;
  return( $generec, $revbest); # defer caller changed.01
}
  
sub get_dna {
  my($fasta, $ref, $start, $stop)= @_; #, $fasta_db_ref
  unless( $ref && $stop>0) { warn "need ref-ID, start, stop location\n"; return; }
 
  require Bio::DB::Fasta;  # FIXME: not portable w/o parts of BioPerl !
  ## need also (Bio::DB::SeqI Bio::Root::Root)
  ## (DB_File GDBM_File NDBM_File SDBM_File) : are these core Perl modules?
  my $havedb= ref $fasta_db;
  unless($havedb) {
    my $db = eval { Bio::DB::Fasta->new($fasta); }
      or die "$@\nCan't open sequence file(s). "; # and return;
    $fasta_db= $db;  
    }
  
  my $seq = $fasta_db->seq($ref, $start => $stop) 
      or return; ## die "cant locate seq of $ref:$start-$stop\n";#? and return;
  $seq= $seq->seq if(ref $seq); # is this weird bioperl change here or not
  #?? Die if no seq? likely bad genome fasta lots of errs
  return $seq;
}


1;

__END__

## revise longest_orf_finder
# sub getAllOrfs { }  # in cdna_proteins.pm;  revised for both strands
# sub orfParts # in cdna_proteins.pm;  revised 
# sub get_orfs {} # in cdna_proteins.pm;  revised 
# sub identify_putative_starts {} # in cdna_proteins.pm;   
# sub identify_putative_stops {} # in cdna_proteins.pm;   
# sub revcomp {} # in cdna_proteins.pm;   
# sub revcomp_coord {} # in cdna_proteins.pm;   
# use vars qw ($currentCode %codon_table);# in cdna_proteins.pm;   
# sub translate_sequence {} # in cdna_proteins.pm;   

## main old introns
#ov2 $intron2splice= ($overlaps) ? kINTRONERROR_OVER : 0; # only choice?
#ov2 my $nintron= 0;
#ov2 if($overlaps) {
#ov2   my $ovh; 
#ov2      if($overlaps =~ /.gz$/) { $ok= open(OVR,"gunzip -c $overlaps |");  $ovh= *OVR; }
#ov2   elsif($overlaps =~ /^(stdin|-)/) { $ovh= *STDIN; $ok=1; }
#ov2   else { $ok= open(OVR,$overlaps); $ovh= *OVR; }
#ov2   die "bad -overlaps=$overlaps" unless($ok);
#ov2   
#ov2   # $overlaplist= 
#ov2   $nintron= collect_Intron_overlaps($ovh); close($ovh);
#ov2 }
#ov2 $hasintrons= ($nintron>0 and scalar(%$overlaplist))?1:0;

# sub testgene_OLDv1
# {
#   my($generecIN, $geneother)= @_;
#   my($addattr,$changed,$mrnachanged)=("",0,0);
#   
#   my @oldCDS = grep{ $_->[2] eq "CDS" } @$generecIN;
# 
#   ## do this before bestorf...  call on generecIN, no call after USE_CDSEXONS 
#   ## BUT ugh reusing generecIN below
#   my($xchanged, $xgenerec, $xfixnote)= intron_short_fix($generecIN); 
#   if($xchanged) {  
#  
#     #.. use generecfix here instead?
#     # my($igenerec, $iexongff, $mrnachanged1)= generecfix( \@generec, \@exongff, $xfixfix); 
#     # @generec= @$igenerec; @exongff= @$iexongff;
#     # $mrnachanged += $mrnachanged1;
#  
#     $generecIN= $xgenerec;
#     $changed+= $xchanged; $addattr .=";" if($addattr); $addattr .= "inshort=$xfixnote";
#   } 
# 
#   my @generec= sort _sortgene grep{ $_->[2] ne "CDS" } @$generecIN; # sorts by genome loc, not stranded
#   my($mrna)  = grep{ $_->[2] eq "mRNA" } @generec;
#   my @exongff= grep{ $_->[2] eq "exon" } @generec;
# 
#   unless($mrna and @exongff > 0) {
#     if($USE_CDSEXONS and $mrna and @oldCDS > 0) { # or caller can s/CDS/exon/ easily, or duplicate like this
#       foreach my $oc (@oldCDS) { my @doc= @$oc; $doc[2]="exon"; push @exongff, \@doc; }    
#       push @generec, @exongff;
#     } else {
#       putgene($generecIN, $geneother,"err=Missing-mrna-exon"); # flag="err=Missing-mrna-exon"
#       return 0;
#     }
#   }
#   
# #  ## do this before bestorf...  call on generecIN, no call after USE_CDSEXONS 
# #  ## BUT ugh reusing generecIN below
# #   my($xchanged, $xgenerec, $xfixnote)= intron_short_fix(\@generec); 
# #   if($xchanged) {   
# #     @generec= @$xgenerec; ## $changed+=$xchanged; $addattr .= "$xfixnote,";
# #   } 
#   
#   
#   # my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid)= @$mrna;
#   my $gstrand= $mrna->[6]; # FIXME for "." ; add CDS bestpro strand (always +??)
#   # *** FIXME 2011Dec : 
#   my $fixstrand=0;
#   my(@genefwd, @generev);
#   if($gstrand eq ".") {
#     $fixstrand=1; # need to test bestpro both strands
#   }
#   
#   my($addattrIN,$changedIN)=($addattr,$changed);
#   my ($changed0,$changed1)=(0,0);
#   for( my $ifix= 0; $ifix <= $fixstrand; $ifix++) { 
#     # BIG loop testing both ways, save protfwd, protrev
#     $changed=$changedIN; $addattr=$addattrIN; # clear for loop step
#     
#     if($fixstrand) {
#       $gstrand=($ifix == 1) ? "-" : "+";
#       ## must clone generec for save @genefwd, @generev
#       @generec= sort _sortgene grep{ $_->[2] ne "CDS" } @$generecIN; # sorts by genome loc, not stranded
#       @generec= map { my @xnew= @$_; \@xnew; } @generec; # clone all
#       map{ $_->[6]= $gstrand } @generec;
#       ($mrna) = grep{ $_->[2] eq "mRNA" } @generec;
#       @exongff= grep{ $_->[2] eq "exon" } @generec;
# 
#       if( $ifix == 1 ) {  @generev= @generec; } 
#       else {  @genefwd= @generec; }
#     }  
#   
#   my $geneid= ($mrna->[8] =~ m/ID=([^;\s]+)/)? $1 : ""; # make one up?
#   (my $geneidfix= $geneid) =~ s/_[CG]\d+$//; # chimera/splitgene _C[12] and multimap _Gnnn id tags
#   
#   my $oldprot= ($mrna->[8] =~ m/protein=([^;\s]+)/) ? $1 : "";
#   @exongff= reverse @exongff if($gstrand eq "-");
#   
# # FIXME5: -samecds needs to use phase0 : partial5 from AUG uses this (eg. at scaf startpos=1)
#   my ($oldStart_b,$oldStart_e)= (0,0); # only for -samecds 
#   if(@oldCDS and $samecds){ @oldCDS= sort _sortgene @oldCDS; # dang
#     my $cstart= ($gstrand eq "-") ? $oldCDS[-1]->[4] : $oldCDS[0]->[3]; 
#     foreach my $ex (@exongff) { my($b,$e)= ($ex->[3], $ex->[4]);
#       if($b <= $cstart and $e >= $cstart) { ($oldStart_b,$oldStart_e)=($b,$e); last; }
#     }
#   }
#   
# #  FIX2: add test intron overlaps : retained = intron inside exon; err = intron rev at splice/inside
# #  FIX3? add here option to filter out huge-span asmrna with poor qual : no intron evidence of long intron; poor prot
# #  # my ($inchanged, $infixnote, $badintrons, @exoninfix) = intron_error_cut( \@exongff);
# #  # if($badintrons) { putgene( $generecIN,"skip=1;badintron=$badintrons"); return 1; } 
# 
# #  FIX? for gaps, chomp off end gaps NNN of getcdna() ? or remove from orfprot
# #  .. need to change exongff also for cdna trimnnn
#   
#   # my($cdna,$cdnalen, $cdnaexonseq)= getcdna( \@exongff, 1, 0, 0, CDNA_TRIMNNN);
#   my($cdna,$cdnalen, $cdnaexonseq, $xfixfix, $xtrimmed)= getcdna( \@exongff, 1, 0, 0, CDNA_TRIMNNN);
#   if($xtrimmed) { 
#   
#     #o# @exongff= @$xfixfix; ## NEED mRNA/gene span update also ****
#     #** FIXME: ^^^ cloned xfix exons NOT registering in exon.gff output;
#     # .. all this dang to-fro of generecIN, @generec fwd/rev is losing changes ***
#     
#     my($igenerec, $iexongff, $mrnachanged1)= generecfix( \@generec, \@exongff, $xfixfix); 
#     @generec= @$igenerec; @exongff= @$iexongff;
#     $mrnachanged += $mrnachanged1;
#     # $mrnachanged++ if( mrna_fixspan($mrna,\@exongff) ); # $mrna->[3] = $newstart;  $newstop;
#     
#     $changed ++; $addattr .=";" if($addattr); $addattr .= "xtrim=$xtrimmed"; 
#   }
# 
#   # add samecds opt using oldStart > need start,stop of 1st CDS exon
#   my($orfprot, $prostart5, $proend3, $bestorf, $utrorf)= 
#     getBestProt("partial", $cdna, \@exongff, $oldStart_b,$oldStart_e);
#   
#   
# #... loop here for intron overlaps >>>>>>>>>>>>>>>>>>>> 
#     ## FIXME2: test cdnain before / after intronfix, intronfix can be bad.
#     # revise : option to annotate errors, but no change to gff.
#     # revised new sub infix
#   my $inerr=0;
#   if($hasintrons) {
#     if($DO_INFIX) {
#       #o my($ichanged,$igenerec,$iexongff,$ibestorf,$infixnote)= 
#       #o     intron_error_fix(\@generec,\@exongff,$bestorf,$gstrand);
#       my($ichanged,$iexongff,$ibestorf,$infixnote)= 
#            intron_error_fix(\@exongff,$bestorf,$gstrand);
# 
#       if($ichanged) {
#         my($igenerec, $mrnachanged1);
#         ($igenerec, $iexongff, $mrnachanged1)= generecfix( \@generec, \@exongff, $iexongff); 
#         @generec= @$igenerec; @exongff= @$iexongff; 
#         $bestorf= $ibestorf; ($orfprot,$prostart5,$proend3)= orfParts($bestorf);
#         $mrnachanged += $mrnachanged1;
#         $changed += $ichanged; $addattr .=";" if($addattr); $addattr .= $infixnote; 
#       }
#     } else { 
#       my($haserr,$inoknote)= intron_error_nofix(\@exongff);
#       $inerr= $haserr;
#       $addattr .=";" if($addattr); $addattr .= $inoknote; # always add if hasintrons
#     }
#   }
#   
# #  if($hasintrons and !$DO_INFIX) {
# #     my($inoknote, $inok,$inerr,$inzip)=("",0,0,0); # use below with cdnain checks
# #     foreach my $ex (@exongff) {
# #       my ($okOrErr, $inannot) = intronoverlaps( @{$ex}[0,3,4,6], 0, 2 ); 
# #       if($okOrErr == 0) { $inzip++; }  
# #       else {
# #         if($okOrErr>0) { $inok++; } elsif($okOrErr<0) { $inerr--; } # should this be inerr++ ? or add key?
# #         # @{$ex}[8] =~ s/$/;inok=$inannot/;  # drop this?
# #         } #? ## need also add to mRNA summary inqual= 
# #       }
# #       # my $inoknote="inqual=$inok/$inerr/$inzip"; # or inqual=9ok,3bad,2none ? or put -err first if > ok; 
# #       my $incode= int (100 * ($inok + $inerr) / ($inok + abs($inerr) + $inzip)); # can be -
# #       $inoknote="inqual=$incode," . ((abs($inerr) > $inok) ? "$inerr/$inok/$inzip" : "$inok/$inerr/$inzip");
# #       $addattr .=";" if($addattr); $addattr .= $inoknote;   
# # 
# #     # FIXME: inqual : overbestgenes1 has better? version; switch to that?
# #     #   .. add inqual score for DO_INFIX also?
# #     #      $incode= int (100 * ($insum + $ierrsum) / $intotal); # can be -
# #     #      $flags .= "ints=$incode," . (($ierrsum) ? "$ierrsum/$insum/$intotal" : "$insum/$intotal") .",$iflag;";
# #            
# #  } elsif($hasintrons and $DO_INFIX) { # old sub infix
# #     my ($inchanged, $infixnote, @exoninfix) = intron_error_cut( \@exongff);
# #     if($inchanged and @exoninfix) {
# #   
# #       if($gstrand eq "-") {   # 2011Dec fix **
# #         my @xs= reverse sort _sortloc @exoninfix;
# #         @exoninfix= @xs;
# #       }
# #       
# #       my($cdnafix, $cdnafixlen)= getcdna( \@exoninfix);
# #       my($fixprot, $fixprostart5, $fixproend5)= getBestProt("partial", $cdnafix, \@exoninfix);
# #       # test also: poor=(utr5>2 or utr3 >2) from getCDSgff()
# #       
# #       ## **?? replace infix when xdrop= error even if fixprot == orfprot
# #       if( (length($fixprot) > 2 + length($orfprot))
# #          or ($infixnote =~ /xdebug/)
# #          or ($infixnote =~ /xdrop=/ and length($fixprot) >= length($orfprot)) ) {
# #         # annotate this
# #         $addattr .=";" if($addattr); $addattr .= $infixnote; $changed++;
# #         ($orfprot, $prostart5, $proend3)= ($fixprot, $fixprostart5, $fixproend5);
# #         # @generec : need to replace exons in @generec here
# #         my @oldexon= map { my @xnew= @$_; \@xnew; } @exongff; # clone all
# #         @generec= grep{ $_->[2] ne "exon" } @generec;
# #         push @generec, @exoninfix;
# #         if($debug>1) { foreach my $ex (@oldexon) {
# #           $ex->[2] = "oldexon"; $ex->[0] =~ s/^/#/; push @generec, $ex;
# #           } }
# #         @exongff = @exoninfix; # preserve oldexon if debug ??
# #         my($newstart,$newstop)= (0,0);
# #         foreach my $ex (@exoninfix) {
# #           my($b,$e)= @{$ex}[3,4];
# #           $newstart= $b if($newstart == 0 or $b < $newstart);
# #           $newstop= $e if($newstop == 0 or $e > $newstop);
# #           }
# #         $mrna->[3] = $newstart; $mrna->[4] = $newstop;
# #       }
# #     }
# # } # old infix sub
# 
# #... end loop for intron overlaps <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 
# 
# #... test input cdna ..............  
#   ## FIXME: cdnain is strand-less ; rev for fixstrand ?
#   ## FIXME2: test cdnain before / after intronfix, intronfix can be bad.
#   
#   my $cdnain  = ""; 
#   my $cdnainlen= 0; my $cdnainNotPrimaryId= 0; my $cdnainIsBest=0;
#   my $cdnagfflen= $cdnalen;
#   if($cdnaseq) {  # add 2011
#     # $cdnain= $cdnaseq{$geneid} || $cdnaseq{$geneidfix}; 
#     $cdnain= $cdnaseq{$geneid};
#     if( not $cdnain and ($geneid ne $geneidfix) ) { $cdnain= $cdnaseq{$geneidfix}; $cdnainNotPrimaryId=1 if($cdnain); }
#    } 
#    
#   if($cdnain) { 
#     $cdnainlen= length($cdnain);
# 
# # FIXed: some bug here preventing Much Better cdnain orf from replacing mapped orf
# # FIXED: ** userev index(cdna) test was bad; check both orfs; and/or check mrna sense=-1 annot from gmap >> trrev of gff
# #   is bestorf_test() ok? yes
# #   .. in cases of < 50% mapped cdnain; could be as simple as userev test, may need to do both strands
# 
# # FIXME2: longer cdnain orf not always better.. as usual.  Some of these are trasm mashups of 2+ genes,
# #   mapped transcript is chimera at same locus giving clue to mashup.  happens more for velvet low kmer
# #   where smaller fragments of 2+ genes are mashed together.   Need to check mRNA chimera tags for clues?
# 
# ## Dang BUG2: chimera-only? cdnabest prot: protein=*MTGHYYYE..HFEY*  ^* is bogus stopcodon, from where ??
# ##   also mayb only with aautr ; also -nostopcodon removes * prefix, but doesnt remove all end *
# ## FIXED below at $oldprot.'*'; 
# 
# # FIXME3: chimera, _C1 and _C2 getting same protein for cdnainIsBest, including bad cases of C1=95%, C2 = 5% mapped
# #  .. here or later allow only 1 of pair to have cdnainIsBest.. should be based on mapped CDS best match to protein
# #     not just longest mapping portion (which often is bogus gene join as UTR).
# 
# ## FIXME4: strand reversal, $userev?, of cdnainIsBest vs gff should not be allowed sometimes.. 
# #      .. or need to reverse strand of output gff
# #   i.e. changing  $prostart5, $proend3, to cdnain set is not valid for rev strand, w/o flip of exon strand
# #   see fixstrand..
# 
# ## FIXME5: new bug: Now gap-cdna with XXX protein winning over cdnain with full protein, same length. 
# #   adjust goodlen test? problem may be index(cdnain,cdna) fails due to XXXX garbage
# 
# # new # need test both ** YES, this cures problems
# # ** BUT need to know if best cdnain is same/rev of gff cdna .. FIXME4
# # .. revert to old, bad; need cdnain strand same as input gff/genome cdna for fully mapped cdna
# # .. but test both strands for split-maps (chimera), partial maps (coverage < 60?), and oneexon trs.
# 
#     use constant CIN_SAMESTRANDasMAP => 1;
#     use constant CIN_BOTHSTRANDS => 0;
# 
#     my $cdnain_isrev=0; # doesnt mean samestrand-as-map, but just that cdnain was reversed.
#     my $cstrandIsOpposite= 0; # 1 if antisense opposite, 0 == sense same as genome map; add -1 == dont know/cant align?
#        ## FIXME6: bug in cstrandIsOpposite, reporting same Sense strand when opposite.
#        ## .. probable failure of $cdnaindex not finding cdnagff inside cdnain.
#        
#     my($orfproti, $prostart5i, $proend3i,$bestorfi, $utrorfi);
#     my $cdnaindex= -1;
#     # my $cdnainSense=1; # or -1 if opposite of mapped strand
#     # from gmap, "sense=-1" means cdnain orient was opposite of mapped strand
#     # ?? what of utrorf on rev strand of bestorf ??
#     
#     ORFCDNAIN: { 
#     my $cdnainrev= revcomp($cdnain);
#     # replace exon count with valid intron count?
#     # my $validgffstrand= ( $inok > $inerr) ? 1 : ($inok+$inerr == 0 and @exongff > 1) ? 1 : 0;
#     my $validgffstrand= ( @exongff > 1 ) ? 1 : 0;
#     # cancel if (abs($inerr) > $inok) 
#     
#     # problem here, some chimera cdna parts can index cdnain both ways; ie cdnain is mashup of 1, or 2 near same.
#     # .. Or, best of both is reverse of this cdna part, but on other cdna part.. case of rev-genes-join
#     # .. cdna part does not map to bestorf part, so cant correct cdsgff from this
#     
#     my $cstrand_gstrand= 0;
#     #drop# my $cdnagffgoodlen= $bestorf->{goodlen}; # not cdnagfflen; dang this is cds not cdna len
#     my $cdnagffgoodlen= $cdnagfflen;
#     my $cdnagood= $cdna; # only for index cstrand_gstrand test
#     if(index($cdnagood,'N') >= 0) { ## No ($cdnagffgoodlen < $cdnagfflen)
#       $cdnagood =~ s/^N+//; $cdnagood =~ s/N+$//; 
#       $cdnagffgoodlen=length($cdnagood);
#       my $xi= index($cdnagood,"NN"); 
#       if($xi>0) { 
#         my $xe= rindex($cdnagood,"NN");
#         if($xi>$cdnagffgoodlen/2) { $cdnagood= substr($cdnagood,$xe+2); }
#         else { $cdnagood= substr($cdnagood,0,$xi); }
#       }
#     } 
#     
#     if( ($cdnaindex= index($cdnain,$cdnagood)) >=0 ) {  $cstrand_gstrand= 1; }
#     elsif( ($cdnaindex= index($cdnainrev,$cdnagood)) >=0 ) { $cstrand_gstrand= -1; }
#     if($cdnaindex < 0) { 
#       ## FIXME6: try harder to match strands; indel problems
#       ## see below also, per exon-gff test for rev
#       my ($xfwd,$xrev)=(0,0);
#       foreach my $xs (@$cdnaexonseq) {
#         if(index($cdnain,$xs) >=0 ) { $xfwd += length($xs); }
#         elsif(index($cdnainrev,$xs) >=0 ) { $xrev  += length($xs); }
#       }
#       $cstrand_gstrand= ($xfwd>$xrev)? 1 : ($xrev>$xfwd)? -1 : 0;
#     }
#     
#     ## here use goodlength tests for cdnain,cdnagff
#     ## testboth is problem now: getting many Anti of gene-mashup (fwd+rev), want only cases of gapfill tr?
#     my $testboth= 0; #  2=best of both 
#     if($fixstrand) { } # DONT test both in this strand loop
#     elsif( $validgffstrand and $cdnainlen < 1.75 * $cdnagffgoodlen) { $testboth=0; }
#     elsif( $cdnainNotPrimaryId or ($cdnainlen >= 1.75 * $cdnagffgoodlen)) { $testboth= 2; } # or not $validgffstrand
#     elsif( $cstrand_gstrand == 0 ) { $testboth= 2; } #? or leave 0
#     
#     if(CIN_SAMESTRANDasMAP) { }   
#     elsif(CIN_BOTHSTRANDS) { $testboth=2; }
# 
#     my ($longorf,$longfull,$orfs);
#     if( $testboth == 2 ) {    
#       ($longorf,$longfull,$orfs) = getAllOrfs($cdnain,"fwd"); # ,"dropnnn"
#       my ($rlongorf,$rlongfull,$rorfs) = getAllOrfs($cdnainrev,"fwd"); # ,"dropnnn"
#       ($cdnain_isrev)= bestorf_test($longorf,$rlongorf);
#       if($cdnain_isrev) {
#         # $cstrand_gstrand= -1 if($cstrand_gstrand == 0); 
#         ## report if cstrand_gstrand unknown: cstrandIsOpposite=-1 ??
#         $cstrandIsOpposite=1 if($cstrand_gstrand == 1);
#         $cdnain= $cdnainrev;
#         ($longorf,$longfull,$orfs)= ($rlongorf,$rlongfull,$rorfs);
#       } else {
#         # $cstrand_gstrand=  1 if($cstrand_gstrand == 0);
#         $cstrandIsOpposite=1 if($cstrand_gstrand == -1);
#       }
#     } elsif( $cstrand_gstrand == -1 ) {
#       $cdnain_isrev=1;
#       $cdnain= $cdnainrev;
#       ($longorf,$longfull,$orfs) = getAllOrfs($cdnain,"fwd"); # ,"dropnnn"
#     } else {  ##  $cstrand_gstrand == 1
#       $cdnain_isrev=0;
#       ($longorf,$longfull,$orfs) = getAllOrfs($cdnain,"fwd"); # ,"dropnnn" 
#     }
#     
#     ($orfproti, $prostart5i, $proend3i, $bestorfi, $utrorfi)= 
#       getBestProtOfOrfs($longorf,$longfull,$orfs, 
#         "partial", $cdnain, \@exongff, $oldStart_b,$oldStart_e);
#     }
#     
# 
#     if($orfproti ne $orfprot) { 
#       ($cdnainIsBest,undef,undef) = bestorf_test($bestorf,$bestorfi);
#        # FIXME bestorf_test: option too high -ratiocdna 1.25; at least should replace cdna XXXX gaps with cdnain perfect
#        # test for near-sameprot but for XXX gaps? 
#        ## if($orfprot =~ /XXX/) { (my $op=$orfprot) =~ s/X/./g; $cdnainIsBest=1 if($orfproti =~ /$op/); }
#        
#       # FIXME2: longer cdnain orf not always better.. as usual.  Some of these are trasm mashups of 2+ genes,
#       # add other checks here if cdnain better than genome-mapped cdna .. esp if cdnain == chimeric mapping
#       #  versus partly unmapped due to genome gap.
#       
#       ## reasons cdnain is best: gaps in genome mapping, NNNs or indels, versus cleaner transcript
#       ## reasons not best: chimeric mapping over same locus or nearby; reversed introns over longer cdnain mapping, ..
#       ## .. maybe not .. leave this off for now.
# #             
# #       # my $cdnabestExplained= bestorf->goodlen/length < cdnainorf->goodlen/len
#       my $cdnainProblems = $cdnainNotPrimaryId or ($hasintrons and $inerr != 0) or ($mrna->[8] =~ /chimera/);
#       
#       if($cdnainIsBest and $cdnainProblems) { 
#       
#         if($cdnainNotPrimaryId) { ## and $mrna->[8] =~ /chim\d=(\w+[^;\s]+)/
#           # cancel if $pctmapped < 25, < 33? <50?
#           my $ingood= $bestorfi->{goodlen} || 1;
#           my $bgood = $bestorf->{goodlen};
#           
#           $cdnainIsBest= 0 if($bgood/$ingood < 0.25); # bad test.. really need other part's stats here
#           ## and/or cancel if mapped CDS are missing or small fraction of cdnain?? need from getCDSgff() or @oldCDS?
#           
#         }
#         
# #         #?? this will cancel some good changes along w/ bad; add annot for scoring instead?
# #         if($hasintrons and $inerr != 0) { $cdnainIsBest= 0; }
# #         elsif($mrna->[8] =~ /chim\d=(\w+[^;\s]+)/) { 
# #           # pick out chimera span and see if near this one.. chim1=scaffold00024:99418-99575:.
# #           use constant CHIOK_MIN => 25000;
# #           my $cloc= $1; my($cr,$cb,$ce)= $cloc =~ m/(\w+):(\d+).(\d+)/;  
# #           my($thisr,$thisb,$thise)= ($mrna->[0], $mrna->[3], $mrna->[4]);
# #           # .. maybe not always cancel; add flag instead? ..
# #           $cdnainIsBest= 0 if($thisr eq $cr and (abs($thisb - $ce) < CHIOK_MIN or abs($thise - $cb) < CHIOK_MIN));
# #         }
#         # elsif($cdnainNotPrimaryId) {} #? what check?
#       }
#       
#     }
#       
#     if($debug) {  #  and not $cdnainIsBest; debug set note; add orig cxlen;aalen here if cdnainIsBest?
#       my $trd = ( $cstrandIsOpposite ) ? "trAnti": "trSense";
#       $trd .= ($cdnain_isrev) ? "Rev" : "Fwd"; # less useful than anti/sense
#       my( $uprot, $ustart, $uend) = orfParts($bestorfi);
#       my $ulen=  $bestorfi->{length};  
#       my $ualen= length($uprot); $ualen-- if(substr($uprot,-1) eq '*');
#       my $upcds  = ($cdnalen>0 && $ulen>0) ? int(100*$ulen/$cdnalen) : 0;
#       my $ucompl= $bestorfi->{complete};
#       $ucompl= ($ucompl==3)?"complete":($ucompl==2)?"partial5":($ucompl==1)?"partial3":"partial";
#       $addattr .=";" if($addattr); 
#       $addattr .= "cdnaorf=$ualen,$upcds%,$ucompl,trbest:$cdnainIsBest,$trd"; #? ;aacdnaoff=$ustart-$uend;aacdnaprot=$uprot";  
#     }
#     
#     if($cdnainIsBest) {
#       # replace?, notice, may need to change CDS like intronfix : depends on diff
#           
#       my $aadif= length($orfproti) - length($orfprot); 
#       my $gadif= $bestorfi->{goodlen} - $bestorf->{goodlen}; 
#       # my $trdif= length($cdnain) - length($cdna); 
#       my $trdif= $cdnainlen - $cdnagfflen;  
#       $cdnalen= $cdnainlen; # change it
#       
#       my $trdif1=$trdif;
#       map{ $_="+".$_ if($_>0); $_= (($_==0)?"eq:":"NE:").$_; } ($aadif,$gadif,$trdif);
#       my $aac= $bestorfi->{complete} - $bestorf->{complete}; # complete == 3,2,1,0
#       
#       # also compare cdna vs cdnain : size, where diff (inside, ends?)
#       # mrna-attr: Coverage=40.0;Identity=99.7 < use this to see diff?
#       my $nx0= @exongff; # no way to count exons in cdnain here.
#       my $trin=0; my $nxeq=0;
#       
#       ## change cdnain_isrev to cstrandIsOpposite here?
#       # my $trd= ($cdnain_isrev) ? "trREV" : "tr";
#       my $trd= ($cstrandIsOpposite) ? "trdAnti" : "trdSense"; # "trSens" ?
#       $trd .= ($cdnain_isrev) ? "Rev" : "Fwd"; # less useful than anti/sense
#       
#       ## FIXME4: strand reversal, $cdnain_isrev?, of cdnainIsBest vs gff should not be allowed sometimes.. 
#       #      .. or need to reverse strand of output gff
#       #   i.e. changing  $prostart5, $proend3, to cdnain set is not valid for rev strand, w/o flip of exon strand
#       #   see fixstrand..
#       
#       ## fixme, use cstrandIsOpposite; also $prostart5, $proend3, are wrong or exongff needs strand change
#       ## BUG, dont change this part strand when it is other chimera part that has bestorf w/ other strand..
#       ## drop this for now, not right for cases seen.; use trdAnti note ; other flag for problems?      
# #       if($cstrandIsOpposite and not $fixstrand) {
# #         if($gstrand eq "-") { $gstrand="+"; } elsif($gstrand eq "+") { $gstrand="-"; }
# #         map{ $_->[6]= $gstrand } @generec;
# #         ($mrna) = grep{ $_->[2] eq "mRNA" } @generec;
# #         @exongff= grep{ $_->[2] eq "exon" } @generec;
# #         ## .. other problems w/ restrand..  oldCDS ? 
# #       }
#   
# 
#       if($cdna eq $cdnain) { $trd .= "eq"; $nxeq=$nx0; }
#       elsif( $cdnaindex >=0 ) { # from index($cdnain|cdnainrev, $cdna)
#         $trd .= "$trdif,in$trin";
#         $nxeq= ($nx0==1)?1:$nx0-1;
#         if($cstrandIsOpposite) {}
#       }
#       # elsif( ($trin=index($cdnain, $cdna)) >=0 ) { $trd .= "in$trdif,$trin"; $nxeq= ($nx0==1)?1:$nx0-1; }  #nxeq maybe; calc?
#       # elsif( ($trin=index($cdnain, revcomp($cdna))) >=0 ) { $trd .="rc$trdif,$trin"; }
#       else { 
#         $trd .= $trdif; # * check each \@exongff > exseq index cdnain ? report if ends or inner x diff
#         # my($cdna2, $cdnalen2, $cdnaexonseq)= getcdna( \@exongff, 1); # use from above cdnaexonseq
#         my ($ix,$le,$xfwd,$xrev)=(0,1,0,0); 
#         foreach my $xs (@$cdnaexonseq) {
#           $ix++; my $xi= index($cdnain, $xs);  my $xr=""; 
#           if($xi>=0) { $xfwd+= length($xs); }
#           else { $xi= index($cdnain, revcomp($xs)); if($xi>=0) { $xr="c"; $xrev+=length($xs); }  }
#           if($xi>=0) { my $xe= length($xs)+$xi; my $xb=$xi+1; $trd.=",xeq$ix:$xr$xb-$xe"; $nxeq++; 
#              if( $le>1 and (my $g= $xb - $le) > 0 ) { $trd.= "g$g"; } $le=1+$xe; }
#           else { $trd.=",xne$ix"; }
#           my $xgap= $xs =~ tr/N/N/;  $trd .="N$xgap" if($xgap>0);
#         }
#         if($xrev > $xfwd and not $cstrandIsOpposite) { $trd =~ s/trdSense/trdAntix/; } # ERROR??
#       }
# 
#       #check further if $nxeq > 0; xinner change?  count NNN gaps in cdna, cdnain; look for gaps at end if inner match
#       if($nxeq > 0) {
#         my $gaps = $cdna =~ tr/N/N/;
#         my $gapi = $cdnain =~ tr/N/N/;
#         if($trdif1>25) { 
#           my($cdna2, $cdnalen2)= getcdna( \@exongff, 0, 100, _max(100,$trdif1));
#           $gaps = $cdna2 =~ tr/N/N/;
#         }
#         if($gaps > $gapi) { $trd.= ",tN:$gaps"; }
#       }       
#       
#       my $cov= ($mrna->[8] =~ m/(cov|Coverage)=(\d+)/) ? $2 : 0;
#       $trd .= ",tcov:$cov" if($cov); # cov < 99 explains
#       my $changenote= "cdnabest=aa$aadif,gd$gadif,dfull:$aac,$trd,nx0:$nx0";
#       $addattr .=";" if($addattr); $addattr .= $changenote; $changed++;
#       ($orfprot, $prostart5, $proend3, $utrorf)= ($orfproti, $prostart5i, $proend3i, $utrorfi);
#     }
#   }
# 
# 
#   if($orfprot) {
#     my($annew,$anold)= stripOldAnnot($mrna) if($REANNOTATE); # also changes mrna
#     
#     my ($cdsgff, $cdsattr)= getCDSgff( \@exongff, $orfprot, $prostart5, $proend3, $cdnalen);
#     # test $prostart5, $proend3 vs old
#     
#     my $diff=0; 
#     $diff=1 if($REANNOTATE or $cdnainIsBest);
#     
#     my $diffcancel=0; my $oldstat="";  
#     if(ref $cdsgff and @$cdsgff > 0) {
#       my @newCDS= sort _sortgene @$cdsgff; # do here.
#       my @oldsave=();
#       # push @generec, @newCDS; # wait till decide if new != old
#       #? @generec= sort _sortgene @generec; # sorts by genome loc, not stranded
# 
#       ## reversed CDS look wrong at end exons (always both ends? for nc>2)    
#       if(@oldCDS > 0) {
#         my($oldal, $newal, $oldprostart5, $oldproend3)= (0,0,0,0);
#         # compare, report in $addattr
#         @oldCDS= sort _sortgene @oldCDS;
#         ## my @newCDS= sort _sortgene @$cdsgff;
#         my $newprot= ($cdsattr =~ m/protein=([^;\s]+)/) ? $1 : "";
#         if (@oldCDS == @newCDS) { $oldstat="ocds=eqn"; } else {  $oldstat="ocds=NEn"; $diff++; }
#         
#         if($oldprot) { 
#           (my $op=$oldprot) =~ s/\*$//; (my $np=$newprot) =~ s/\*$//; 
#           $oldal= length($op);
#           $newal= length($np);
#           my $oldap= ($cdnalen>0) ? int(0.5 + 300 * $oldal / $cdnalen) : 0;
# 
#           my $eq=($op eq $np)?1:0; $diff++ unless($eq); 
#           my $da= $newal - $oldal; $da="+$da" if($da>0);
#           $oldstat .= ($eq) ? ",eqaa" : ",NEaa:$da"; 
#           
#           my $astat=0;
#           if($oldprot =~ /^M/) { $astat |= 1; }
#           if($oldprot =~ /\*$/) { $astat |= 2;} # ** AUG proteins lack '*' unless pasa-updated
#           elsif($geneid =~ /AUG/) { $astat |= 2; } #  also check for inner X == augustus-fake for stop ?
# 
#           my $istop= index($oldprot,'*');
#           ## find single 'X' inside prot, but allow this for NNN genome: SFXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXQP
#           if($geneid =~ /AUG/ and (my $ix=index($oldprot,'X')) >0 ) { 
#             if(substr($oldprot,$ix,3) eq "XXX") { } # ignore 2+; really should look at genome dna to decide
#             else { $istop= $ix if($istop < 1 or $ix<$istop); }
#             }
#           # if($istop < 1 and $geneid =~ /AUG/) { $istop= index($oldprot,'X'); }
# 
#           # my $prostat = ($p5 and $p3) ? "partial": ($p5)? "partial5" : ($p3)? "partial3" :"complete";
#           if($istop > 0 and $istop < length($oldprot)-1) { $astat= "partialinner"; } # innerstop ?
#           elsif($astat == 3) { $astat="complete"; }
#           elsif($astat == 1) { $astat="partial3"; }
#           elsif($astat == 2) { $astat="partial5"; }
#           else { $astat="partial"; }
#           $astat =~ s/^(.)/o/; # dont confuse w/ new in scans
#           $oldstat .= ";oaaln=$oldal,$oldap%,$astat";
#           
#         } else { $diff++; }
# 
#    # FIXME2: exist CDS : allow for problem cases like end-of-scaffold partials
#         my $partialAtEndOk= 0;
#         for(my $i=0; $i<@oldCDS; $i++) { 
#           my $oc= $oldCDS[$i];
#           $partialAtEndOk=1 if($oc->[3] < 450); # dont know high end to check
#           unless($i<@newCDS and $oc->[3] == $newCDS[$i]->[3] and $oc->[4] == $newCDS[$i]->[4]) { 
#             $diff++; $oldstat .= ",NE$i";  ## BUG: clone oc
#             if($debug>1) { my @doc= @$oc; $doc[2]="oldCDS"; $doc[0] =~ s/^/#/; push @oldsave, \@doc; }
#             }
#           }
#           
#        
#         do{ $diff=0; $diffcancel=1; } if(!$NODIFCAN and $partialAtEndOk and $oldstat =~ /artial/ and $newal < $oldal);
# 
#    # FIXME4: cancel changes that replace all CDS with completely new CDS, eg for transposon-spans
#         $oldprostart5= ($gstrand eq "-") ? $oldCDS[-1]->[4] : $oldCDS[0]->[3]; 
#         $oldproend3  = ($gstrand eq "-") ? $oldCDS[0]->[3] : $oldCDS[-1]->[4]; 
#         $diff++ if($diff==0 and ($oldproend3 != $proend3 or $oldprostart5 !=  $prostart5));   # off by -1,-2 bugs   
# 
#         #  oldc:   [-----]
#         #  newc:           [-----]  ; bad = shifted too much
#         #  newc:        [-----]     ; bad "
#         #  newc:     [---]          ; ok  = not shifted, same stop
#         
#         if($diff and $oldal>0) { # and $oldstat !~ /artialinner/
#           ## .. bugs here ??
#           #my($tb,$te)=($prostart5,$proend3); ($tb,$te)= ($te,$tb) if($tb>$te);
#           #my($lb,$le)=($oldprostart5,$oldproend3); ($lb,$le)= ($le,$lb) if($lb>$le);
#           my($tb,$te)=($newCDS[0]->[3],$newCDS[-1]->[4]); ($tb,$te)= ($te,$tb) if($tb>$te);
#           my($lb,$le)=($oldCDS[0]->[3],$oldCDS[-1]->[4]); ($lb,$le)= ($le,$lb) if($lb>$le);
#           if($lb==0 or $le==0 or $tb==0 or $te==0) { $tb=$te=$lb=$le=0; } # bug somewhere ...
#           
#           # argg; location stats not best here, introns have big effect.
#           # .. this cancels most span increases...
#           if(($tb < $lb and $te < $le) or ($tb > $lb and $te > $le)) { # CDS shifted along exons
#             my ($bb,$be)= ( _max($tb,$lb), _min($te,$le) );
#             my $maxo= 1 + $be - $bb; # ** not abs
#             my $leno= _min( abs($le - $lb), abs($te - $tb)) || 1;
#             my $pover= $maxo/$leno; # neg desired: cancel any case of maxo < 0 
#             #??do{ $diff=0; $diffcancel=1; } if($pover < 0.33); # cancel change, too little cds similar
#             do{ $diff=0; $diffcancel=1; } if(!$NODIFCAN and $pover < 0); # cancel change, too little cds similar
#           }
#         }
# 
#         $cdsattr .=";$oldstat";
#         $changed++ if($diff>0);
#       } else {
#         $changed++; $diff=1;
#       }
# 
#       if($diff) {
#         push @generec, @newCDS; 
#         push @generec, @oldsave; 
#       } else {
#         push @generec, @oldCDS;  
#       }
#     }
#     # else { $diff++; } ## is this error? cdnainIsBest = cdnain best?
#     #? else { push @generec, @newCDS;  } # no @$cdsgff no oldCDS
#   
#     ## add to addattr: count utrs (and span?) and add utrx= if >2
#     # change getCDSgff() to return utr-exons?
#     # if($u5>2 or $u3>2) { $un=$u5+$u3; $g[0]=~s/$/;utrx=$u5,$u3/; } 
# 
#     ## fixme: diff==0, dont remove old attr, add new or not? 
#     #     or make option? sometimes always replace? use ocds=eqn,eqaa as flag for no change
#     # .. except for some changes above result in NEaa but are cancelled.
#     
#     if($utrorf) {
#       my( $uprot, $ustart, $uend) = orfParts($utrorf);
#       my $ulen=  $utrorf->{length};  
#       my $ualen= length($uprot);  
#       my $upcds  = ($cdnalen>0 && $ulen>0) ? int(100*$ulen/$cdnalen) : 0;
#       my $ucompl= $utrorf->{complete};
#       $ucompl= ($ucompl==3)?"complete":($ucompl==2)?"partial5":($ucompl==1)?"partial3":"partial";
#       $addattr .=";" if($addattr); 
#       $addattr .= "aautrlen=$ualen,$upcds%,$ucompl;utroff=$ustart-$uend;utrprot=$uprot";  
#     }
#     
#     if($cdsattr and $diff==0) {
#       # check prot, add '*' if missing and complete
#       # remove prot, new aalen if  $diffcancel
#       unless( $mrna->[8] =~ m/;protein=/ or $diffcancel) {  ## should set diff=1 for this
#         $diff=1; # $mrna->[8] =~ s/$/;$cdsattr/; 
#         } ## add if missing
#       if($diffcancel) { 
#         my ($st)= $oldstat =~ m/ocds=([^;\s]+)/; 
#         my ($nal)= $cdsattr =~ m/aalen=([^;\s]+)/;
#         $mrna->[8] =~ s/$/;ocds=cancel:$st,nal:$nal/; 
#         
#       } elsif(not $NoStopCodon and $cdsattr =~ /complete/ and $oldprot =~ /\w/ and $oldprot !~ /\*$/) {
#         my $ops= $oldprot.'*';  # this WAS bug for protein=*Maaaa* 
#         $mrna->[8] =~ s/;protein=$oldprot/;protein=$ops/; # no end ;
#       }
#       
#     } # elsif()
#     
#     if($cdsattr and $diff>0) {  
#         ## cdsattr keys: cxlen,aalen,protein,utrx
#       my @keys= $cdsattr =~ m/(\w+)=/g;
#       my $keyp= join('|',@keys);
#       $mrna->[8] =~ s/[;]?($keyp)=([^;\n]+)//g;
#       $mrna->[8] =~ s/$/;$cdsattr/;
#     }
#     
#     $mrna->[8] =~ s/$/;$addattr/ if($addattr);
# 
#   } else { # no orfprot : can happen, but error if oldCDS << need to restore those here ??
#     push @generec, @oldCDS if(@oldCDS);
#     my $addattr="ocds=cancel:NEn,aatrans-failed";
#     $mrna->[8] =~ s/$/;$addattr/;
#   }
# 
#   if($ifix == 1) {  @generev= @generec; $changed1=$changed; }  # already cloned, now preserve
#   else {  @genefwd= @generec;  $changed0=$changed; }
#   
#   } # BIG loop for fixstrand 
# 
#   ## somewhere add nostopcodon opt:
#   #  if($NoStopCodon and substr($orfprot,-1) eq '*') { $orfprot =~ s/\*$//; }
#   
#   my($generec,$revbest)= (\@generec,0);
#   if($fixstrand) {
#     ($generec,$revbest)= fixstrand( \@genefwd, \@generev);
#     $changed= ($revbest) ? $changed1 : $changed0;
#   }
#   
#   my $gflags=($mrnachanged)?"mrnachanged=$mrnachanged":"";
#   # FIXME15: add opt output -cdna -cds sequences, separate files (also -aa seq sep file instead of mrna.gff attr)
#   putgene( $generec, $geneother, $gflags);
#   return ($changed) ? 1 : 0;
# }
# # end sub testgene() too long ......

# sub intron_error_nofix {
#   my($exongff) = @_;  
#   my($inoknote, $inok, $inerr, $inzip)=("",0,0,0); # use below with cdnain checks
#   
#   foreach my $ex (@$exongff) {
#     my ($okOrErr, $inannot) = intronoverlaps( @{$ex}[0,3,4,6], 0, 2 ); 
#     if($okOrErr == 0) { $inzip++; }  
#     else {
#       if($okOrErr>0) { $inok++; } elsif($okOrErr<0) { $inerr--; } # should this be inerr++ ? or add key?
#       ## @{$ex}[8] =~ s/$/;inok=$inannot/;  # drop this?
#       } #? ## need also add to mRNA summary inqual= 
#     }
#     
#   ## my $inoknote="inqual=$inok/$inerr/$inzip"; # or inqual=9ok,3bad,2none ? or put -err first if > ok; 
#   my $incode= int (100 * ($inok + $inerr) / ($inok + abs($inerr) + $inzip)); # can be -
#   $inoknote="inqual=$incode," . ((abs($inerr) > $inok) ? "$inerr/$inok/$inzip" : "$inok/$inerr/$inzip");
#   #o# $addattr .=";" if($addattr); $addattr .= $inoknote;   
# 
#   # FIXME: inqual : overbestgenes1 has better? version; switch to that?
#   #   .. add inqual score for DO_INFIX also?
#   #      $incode= int (100 * ($insum + $ierrsum) / $intotal); # can be -
#   #      $flags .= "ints=$incode," . (($ierrsum) ? "$ierrsum/$insum/$intotal" : "$insum/$intotal") .",$iflag;";
#   my $haserr= ($inzip>0 or $inerr<0)?$inerr:0;
#   return($haserr,$inoknote);          
# }
# 
# 
# sub _sortinloc1 { 
#   #  my($ref,$start,$stop,$strand)= @{$ft}[0,3,4,6];
#   return ($a->[0] cmp $b->[0]) # ref : should all be same for gene record
#       || ($a->[3] <=> $b->[3]) # begin
#       || ($a->[4] <=> $b->[4]) # end
#       ;
# }
# sub _sortinloc { 
#   #  my($ref,$start,$stop,$strand)= @{$ft}[0,3,4,6];
#   return ($a->[0] cmp $b->[0]) # ref : should all be same for gene record
#       || (abs($a->[3] - $b->[3])>90 && ($a->[3] <=> $b->[3])) # begin1 <
#       || ($b->[5] <=> $a->[5]) # score > ; should be start+=50 first then score
#       || ($a->[3] <=> $b->[3]) # begin <
#       || ($a->[4] <=> $b->[4]) # end <
#       ;
# }
# 
# 
# =item intron_error_fix
#   
#   using *valid* intron.gff overlaps, check if exons should be cut to match introns
#   
# =cut
# 
# sub intron_error_fix  { # infix
#   my($exongff,$bestorf,$gstrand) = @_; ## ,$addattr
#   #o my($generec,$exongff,$bestorf,$gstrand) = @_; ## ,$addattr
#   # return($changed,$generec,$exongff,$bestorf,$infixnote);
#   my $changed=0; my $mrnachanged=0;
#   
#   my($inchanged, $infixnote, @exoninfix) = intron_error_cut( $exongff);
#   if($inchanged and @exoninfix) {
# 
#     if($gstrand eq "-") {   # 2011Dec fix **
#       my @xs= reverse sort _sortloc @exoninfix;
#       @exoninfix= @xs;
#     }
# 
#     my($cdnafix, $cdnafixlen, undef, $xfixfix,$xtrimmed)= getcdna( \@exoninfix, 0, 0, 0, CDNA_TRIMNNN);
#     # $cdna,$cdnalen,\@asexons,$exons,$xchange
#     if($xtrimmed) { 
#       @exoninfix= @$xfixfix;
#       # see below: mrna_fixspan($mrna,\@exoninfix);  
#     }
# 
#     my($fixprot, $fixprostart5, $fixproend5, $fixorf)= getBestProt("partial", $cdnafix, \@exoninfix);
#     my($orfprot,$prostart5,$proend3)= orfParts($bestorf); ## use param bestorf not orfParts..
#     # test also: poor=(utr5>2 or utr3 >2) from getCDSgff()
# 
#     ## **?? replace infix when xdrop= error even if fixprot == orfprot
#     if( (length($fixprot) > 2 + length($orfprot))
#        or ($infixnote =~ /xdebug/)
#        or ($infixnote =~ /xdrop=/ and length($fixprot) >= length($orfprot)) ) 
#     {
#       $changed++; 
#       #o# $addattr .=";" if($addattr); $addattr .= $infixnote; # annotate this
#       #o# ($orfprot, $prostart5, $proend3)= ($fixprot, $fixprostart5, $fixproend5);
# 
#       return($changed,\@exoninfix,$fixorf,$infixnote); 
#       #nn $bestorf= $fixorf; # new return param
#       #nn $exongff= \@exoninfix; # changes here
#       
#       #.. defer to caller
#       # my $mrnachanged1= 0;
#       # ($generec, $exongff, $mrnachanged1)= generecfix($generec, $exongff, \@exoninfix); 
#       # $mrnachanged += $mrnachanged1;
#       
# ## old2.............
# #       my @oldexon= map { my @xnew= @$_; \@xnew; } @$exongff; # clone all
# #       my @generecfix= grep{ $_->[2] ne "exon" } @$generec;
# #       push @generecfix, @exoninfix;
# #       if($debug>1) { foreach my $ex (@oldexon) {
# #         $ex->[2] = "oldexon"; $ex->[0] =~ s/^/#/; push @generecfix, $ex;
# #         } }
# #         
# #       $exongff = \@exoninfix; # preserve oldexon if debug ??
# #       my($mrna) = grep{ $_->[2] eq "mRNA" } @generecfix;
# #       $mrnachanged++ if( mrna_fixspan($mrna,$exongff) );
# #       $generec= \@generecfix;
# ## old1      
# #       my($newstart,$newstop)= (0,0);
# #       foreach my $ex (@exoninfix) {
# #         my($b,$e)= @{$ex}[3,4];
# #         $newstart= $b if($newstart == 0 or $b < $newstart);
# #         $newstop= $e if($newstop == 0 or $e > $newstop);
# #         }
# #       $mrna->[3] = $newstart; $mrna->[4] = $newstop;
#       
#     }
#   }
#     
#   return($changed,$exongff,$bestorf,$infixnote);  
#   #o return($changed,$generec,$exongff,$bestorf,$infixnote); ## $addattr,$change
# }
# 
# sub intronoverlaps
# {
#   my ( $ref,$tb,$te,$to, $inmaxscore, $validonly)= @_;    # exon here
#   my ( %didid,@overs,@ovok);
#   return 0 unless($overlaplist->{$ref});
#   $validonly ||= 0;
#   
#   my($tb1,$te1, $tb2, $te2)=(0) x 4;
#   if(1) { # always ($intron2splice == kINTRONERROR_OVER)
#     # input == exon, over= intron, test if overlap is splice end or internal
#     ($tb1,$te1, $tb2, $te2)= 
#       ($to eq "-") ? ($te - SPLICEX,$te,$tb,$tb + SPLICEX) : ($tb,$tb + SPLICEX,$te - SPLICEX,$te);  
#   }
#   
#   $inmaxscore ||= 0; 
#   my $pmaxscore= $inmaxscore * $pMAXSCORE;
#   my($nover)= (0); ## NEED MAX inscore over all exons not just 1; global? reset per gene
#   
#   my($linb,$line, $llo,
#      $okleft, $okright, $errleft, $errright, $errinside)= (0) x 10;
#   my ($ib1, $ib2)= (int($tb/$BINSIZE) , int($te/$BINSIZE));
#   for (my $ib = $ib1; $ib <= $ib2; $ib++) 
#   {
#     $overlaplist->{$ref}{$ib} or next;
#     my @locs= @{$overlaplist->{$ref}{$ib}};
#     
#     # sort, sameorient 1st, so bidir introns are less problem : NOT working
#     @locs= ( (sort _sortinloc grep { $_->[6] eq $to } @locs),  (sort _sortinloc grep { $_->[6] ne $to } @locs) );
#     
#     foreach my $rloc (@locs) {
#       #new# my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid,$inb,$ine]; 
#       
#       my ($lb,$le, $inscore,$lo,$oid,$inb,$ine)= @{$rloc}[3,4,5,6,9,10,11];  # intron here : 10-11 == full span, 3-4 == splice
#        # FIXME now is intron ENDs
#        # FIXME: need full intron lb,le and also splice ends
#        # FIX-11Dec: use inscore/count to decide if valid cut: if validin.c=1000 and insidein.c=10, skipit.
#        
#       next if($didid{$oid.$lb.$le}++);
#       #NO, count overs# next if($inb < $line && $ine > $linb); # overlap last intron, rev: skip
#       # ($linb, $line, $llo)= ($inb, $ine, $lo);
#       
#       my $over  = ($tb <= $le && $te >= $lb) ? 1 : 0;      
#       # $over=0 if($stranded and ($to =~ /[+-]/ and $lo =~ /[+-]/ and $to ne $lo)); 
#       # my $inside= ($over and $tb <= $lb && $te >= $le) ? 1 : 0; ## intron inside  
#       # $inside=0 if($stranded and ($to =~ /[+-]/ and $lo =~ /[+-]/ and $to ne $lo)); 
#  
#       # if($inside) {  push @overs, [$lb,$le]; } else
#       
#       if($over) { #always($intron2splice == kINTRONERROR_OVER)   # test for splice reversed errors
#         ($linb, $line, $llo)= ($inb, $ine, $lo);
#         
#         my $ok=0; # note: lb,le here are one splice end span of intron: 3 bp?
#         my $samestrand= ($to =~ /[+-]/ and $lo =~ /[+-]/ and $to ne $lo) ? 0 : 1;
# 
#         if($tb1 <= $le && $te1 >= $lb) { if($samestrand) { $ok=1; $okleft+=$inscore; } else { $ok= kINTRONERROR_OVER; $errleft+=$inscore; } } # $errt="rev" unless $samestrand
#         elsif($tb2 <= $le && $te2 >= $lb) { if($samestrand) { $ok=2; $okright+=$inscore; } else { $ok=kINTRONERROR_OVER; $errright+=$inscore; } } # -1 or -2?
#         elsif(($tb + 2*SPLICE <= $inb && $te - 2*SPLICE >= $ine) and $samestrand) { $ok= kINTRONERROR_INSIDE; $errinside+=$inscore; }  #-3
#         
# #         if($tb1 <= $le && $te1 >= $lb) { $ok= ($samestrand) ? 1 : kINTRONERROR_OVER; } # $errt="rev" unless $samestrand
# #         elsif($tb2 <= $le && $te2 >= $lb) { $ok= ($samestrand) ? 2 : kINTRONERROR_OVER; } # -1 or -2?
# #         elsif(($tb + 2*SPLICE <= $inb && $te - 2*SPLICE >= $ine) and $samestrand) { $ok= kINTRONERROR_INSIDE; }  #-3
#                   # ^^ same strand inside  == retained intron err
#         #old# elsif(($tb + 2*SPLICE <= $lb && $te - 2*SPLICE >= $le) and $samestrand) { $ok= kINTRONERROR_INSIDE; } 
#         # else what?
#         
#         if($ok < 0) {
#           $nover++;
#           # $enderr |= $ok; # 1,2 or 3=both
#           if($inscore < $pmaxscore) { } ## $inmaxscore * $pMAXSCORE)  # skip weak cases
#           elsif($ok == kINTRONERROR_INSIDE) { push @overs, [$inb,$ine, $ok]; } # need error type also: retained vs strand-err
#           else { push @overs, [$lb,$le, $ok]; } # need error type also: retained vs strand-err
#           #? need ok1, ok2 for both ends of exon to say if 1 end is ok?
#         } elsif( $ok > 0) {
#           $nover++;
#           if($inscore > $inmaxscore)
#             { $inmaxscore= $inscore; $pmaxscore= $inmaxscore * $pMAXSCORE; }
#           # $endok |= $ok; # 1,2 or 3=both; want count of each end ok?
#           #? push @overs, [$lb,$le, $ok]; # UPDATE 11Dec, see below 
#           push @ovok, [$lb,$le, $ok, $inscore]; #? return which end of exon is supported by intron (or both)
#           # push @spliceok, ($ok == 2) ? (($to eq "-") ? $tb : $te) : (($to eq "-") ? $te : $tb);
#         }
#         # next;
#       } 
#         
#       # push @overs, [$lb,$le] if ($over);
#       }
#   }
#   
#   ## add for annotations, no changes
#   if($validonly >= 2) {  ## is allowalts synonym of this now?
#     # my $nok=@ovok; my $nerr=@overs; 
#      ## dont want count of introns here, want to know if either/both ends of exon are supported or not
#      ## and if there is internal intron kINTRONERROR_INSIDE
#      
#     my $okOrErr = ($nover==0) ? 0 : ($okleft >= $errleft and $okright >= $errright and $errinside < $inmaxscore) ? 1 : -1;
#     # my $inannot = ($nover==0) ? 0 : join ",", $okOrErr, $okleft,$okright, -$errleft,-$errright,-$errinside;
#     # exon annot:  ;inok=[1/-1],20,10,-2,-20,-3; want all this data?  or inok=ngood,-nbad | inok=-nbad,ngood
# #     my $inannot= ($okOrErr == 0) ? 0 
# #         : ($okOrErr>0) ? join ",", ($okleft+$okright), (-$errleft-$errright-$errinside) 
# #         : join ",", (-$errleft-$errright-$errinside),($okleft+$okright);
#     my $inannot= ($okOrErr == 0) ? 0 
#         : ($okOrErr>0) ? join ",", $okleft.'l',$okright.'r', -$errleft.'l',-$errright.'r',-$errinside.'i'
#         : join ",", -$errleft.'l',-$errright.'r',-$errinside.'i',$okleft.'l',$okright.'r';
#     return($okOrErr, $inannot); ##, $okleft, $okright, $errleft, $errright, $errinside); 
#     }
#     
#   return \@ovok if($validonly);
#   return 0 if(@overs > 1 and $allowalts); # dont force cut where we want alts and they exist
#    
#   # ** FIXME: this makes bad cuts, likely where 2+ alt-introns are overlapped
#   my @opens=();
#   if(@overs) {
#     my $okover= scalar(@ovok); #** Need ovok in overs to let bidir good superceed bad/rev introns
#     my($bb,$be)= ($tb,$te); my $errlast= 0; my($llb,$lle)=(0,0);
#     #OFF. @overs= sort _sort_over @overs; # NOT NOW, see above sort
#     
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
#         if( 1 or $errcode < 0) { push @opens, [$b1,$e1,$errcode,$okover]; }
#         $bb= $le+1; 
#         } 
#       elsif($lb > $te) { last; } #?
#       ($llb,$lle)= ($lb,$le);
#       $errlast= $errcode;  #?? if($errcode==kINTRONERROR_INSIDE or $errlast==0); # upd 11dec
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
# 
# 
# sub intron_error_cut
# {
#   my($exongff)= @_;
#   
#   my $xdebug=($debug>1)?"/xdebug":"";
#   
#   my @exnew=(); 
#   my $changenote=""; 
#   my $changed= 0; my $xi=0;  my $droperr=0;
#   my $inmaxscore=0;
#   foreach my $ex (@$exongff) {
#     my $valids = intronoverlaps( @{$ex}[0,3,4,6], 0, 1 ); 
#     if(ref $valids and @$valids > 0) {
#       foreach my $xloc (@$valids) {
#         my($xb,$xe,$errcode,$score)= @$xloc;
#         $inmaxscore= $score if($score>$inmaxscore);
#         }
#     }
#   }
#   
#   foreach my $ex (@$exongff) {
#   
#     # ** FIXME: this makes bad cuts, likely where 2+ alt-introns are overlapped
#     my $exonfix = intronoverlaps( @{$ex}[0,3,4,6], $inmaxscore ); 
#     $xi++;
#     # result= [array] of [$bb,$te,$errcode] ; errcode == kINTRONERROR_OVER/bad splice, kINTRONERROR_INSIDE/retained in
# 
#     # 2011Dec: DAMN, reverse(@$exonfix) for strand "-" : see below
#     # 2011Dec: DAMNN, got dupl xcut locations, all with ? xdrop=-1/1009/12.2
# 
# # orig:
# # scaffold_2      caca11r39cuf8   mRNA    6247942 6251164 225     +       .       ID=caca11r39cuf8_Gsc2g4225t1;
# # scaffold_2      caca11r39cuf8   exon    6247942 6248192 225     +       .       Parent=caca11r39cuf8_Gsc2g4225t1;xi=1;
# # scaffold_2      caca11r39cuf8   exon    6248434 6248780 225     +       .       Parent=caca11r39cuf8_Gsc2g4225t1;xi=2;
# # scaffold_2      caca11r39cuf8   exon    6248876 6251164 225     +       .       Parent=caca11r39cuf8_Gsc2g4225t1;xi=3;
# #.. after xcut ..  intronfix=1,xcut:-3/792/3.1,xdrop:-1/498/3.2,
# #  # Thecc1EG007301t1	rna8b:r8caca11r39cuf8_Gsc2g4225t1  << BAD xcut, DUPL exons,CDS; from overlapped introns
# # scaffold_2	caca11r39cuf8	exon	6247942	6248192	225	+	.	Parent=caca11r39cuf8_Gsc2g4225t1;xi=1;
# # scaffold_2	caca11r39cuf8	exon	6248434	6248780	225	+	.	Parent=caca11r39cuf8_Gsc2g4225t1;xi=2;
# # scaffold_2	caca11r39cuf8	exon	6248876	6249667	225	+	.	Parent=caca11r39cuf8_Gsc2g4225t1;xi=3;;xcut=-3/792/3.1
# #   ^^^^^^ bad exon1? or good
# # scaffold_2	caca11r39cuf8	exon	6248876	6251164	225	+	.	Parent=caca11r39cuf8_Gsc2g4225t1;xi=3;;xdrop=-1/498/3.2
# #   ^^^^^^ bad exon2 -- xb should be after 6249667
#     
#     if(ref $exonfix and @$exonfix > 0) {
#       my $fixerr=0;
#       my @exadd;
#       my $xj=0;
#       my($lxb,$lxe,$didx)=(0,0,0);
#       foreach my $xloc (@$exonfix) {
#         my($xb,$xe,$errcode,$hasok)= @$xloc;  $xj++;
#         my $xw= 1+$xe-$xb;
#         $fixerr++ if($xw < $MINEXON);
#         
#         if($errcode == kINTRONERROR_OVER) {
#           # drop exon span? but annotate what?
#           my $xan="xdrop=$errcode$xdebug/$xw/$xi.$xj";
#           $changenote .= "$xan,";
#           # careful: dont drop 1st w/ intron error .. keep as last exon before error?
#           # exons are 5' > 3' order and prot/CDS probably only in 5' before err
#           # .. but doesnt look wrong w/o this.
#          if($didx) { # upd 11dec
#             my @xx= @$ex;  $xx[3]= $xb; $xx[4]= $xe;
#             $xx[8] =~ s/$/;$xan/;
#             push @exadd, \@xx;
#             ($lxb,$lxe)=($xb,$xe); $didx=0;
#           } elsif($droperr==0 and $hasok) {
#             my @xx= @$ex;  #? $xx[3]= $xb; $xx[4]= $xe; #< change to this or not?
#             ## lxe >> $xx[3]= $xb; $xx[4]= $xe;
#             $xx[8] =~ s/$/;$xan/;
#             push @exadd, \@xx;
#            }
#           $droperr++;
#           
#         } elsif($errcode == kINTRONERROR_INSIDE) {
#           my $xan="xcut=$errcode$xdebug/$xw/$xi.$xj";
#           $changenote .= "$xan,";
#           my @xx= @$ex;  $xx[3]= $xb; $xx[4]= $xe;
#           $xx[8] =~ s/$/;$xan/;
#           push @exadd, \@xx;
#           ($lxb,$lxe)=($xb,$xe); $didx++;
#         }
#       }
#       
#       if($fixerr>0) { 
#         push @exnew, $ex;
#       } else {
#         $changed++; push @exnew, @exadd;
#       }
#       
#     } else {
#       push @exnew, $ex;
#     }
#   }
#   
#   $changenote =~ s/=/:/g;
#   $changenote= "intronfix=$changed,$changenote" if($changed);
#   # return $badintrons == long, no support
#   return ($changed, $changenote, @exnew);
# }
# 
# sub collect_Intron_overlaps
# {
#   my($gff)= @_;
#   my ($nr,$nx)=(0,0);
#   while(<$gff>){
#     next unless(/^\w/); chomp;
#     my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr)= split"\t";
#     $tattr ||="";
#     
#     # if($passtypes and "$typ.$src" !~ m/$passtypes/) { next; } # pass other types
#     #^ drop passtypes for mrnatypes,exontypes
#     next unless($typ =~ /^($introntypes)$/);
# 
#     # ?? only Introns here? need splice-site calcs ; no ID= for introns...
#     $nr++;
#     my($gid,$pid); 
#     if($tattr =~ m/\bID=([^;]+)/) {  $gid=$1; }
#     if($tattr =~ m/\bParent=([^;]+)/) {  $pid=$1; 
#       $gid=$pid unless($gid and ($typ =~ /^($mrnatypes)$/)); 
#       }
#     unless(defined $gid) { $gid = "N".$nr; }
# 
#     #? only kINTRONERROR_OVER and kINTRONERROR_INSIDE
#     my($inb,$ine)= ($tb,$te); # full intron span, keep
#     # if($intron2splice == kINTRON2SPLICE_OVER or $intron2splice == kINTRONERROR_OVER) 
#     if(1) { # always intron here
#       my($s1b,$s1e,$s2b,$s2e)= 
#         ($to eq "-") ? ($te+1,$te + SPLICE,$tb - SPLICE,$tb-1) : ($tb - SPLICE,$tb-1,$te+1,$te + SPLICE); # 3bp + 1 shift
#       ($tb,$te)= ($s1b,$s1e);
#       my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid,$inb,$ine]; 
#       my @bins= (int($tb/$BINSIZE) .. int($te/$BINSIZE));
#       foreach my $ib (@bins) { push( @{$overlaplist->{$ref}{$ib}}, $rloc); }
#       ($tb,$te)= ($s2b,$s2e);  #? change gid, oid?
#     }  
#     
#     my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid,$inb,$ine]; 
#     my @bins= (int($tb/$BINSIZE) .. int($te/$BINSIZE));
#     foreach my $ib (@bins) { push( @{$overlaplist->{$ref}{$ib}}, $rloc); } # $generec
#   }
#   
#   warn"#collect_overlaps n=$nr\n" if $debug;
#   return $nr;# return \%overlaps;
# }
# 

#ov2 sub cdna_bestorf
#ov2 {
#ov2   my @cdnaid= sort keys %cdnaseq; my $ncdna= @cdnaid;
#ov2   
#ov2   warn "#REPLACED by cdna_bestorf.pl, 201207\n";
#ov2   warn "#findcds from cdnaseq:$cdnaseq n=$ncdna\n" if $debug;
#ov2   my $ngood=0;
#ov2   ## my $MINAA= $MINID_CDS || 40; 
#ov2   my $MINTR= $MINAA * 3;
#ov2   my @rev=($action =~ /rev/)? (1) : ($action =~ /fwd/)? (0) : (0,1);
#ov2   my $fullpart=($action =~ /full/)?"full":($action =~ /long/)?"longpart":"bestpart";
#ov2 
#ov2 # ADD maybe: detect joined/fused genes using 
#ov2 #  1. cds span < ~ 1/2 trspan, leaving long utr, << annotate these cases, require aa complete/partial5
#ov2 #  2. check the utr for long cds
#ov2  
#ov2   foreach my $id (@cdnaid) {
#ov2     my $cdnain= $cdnaseq{$id}; 
#ov2     my $clen= length($cdnain); 
#ov2     my $cdnabest= $cdnain;
#ov2     my($orfprot, $prostart5, $proend3, $bestorf, $utrorf, $orflen, $isrev, $aalen, $pcds, $compl)= (0) x 20;
#ov2     if($clen < $MINTR) { 
#ov2       # warn?
#ov2     } else {
#ov2     for my $rev (@rev) {
#ov2       if($rev==1){  $cdnain= revcomp($cdnain); } 
#ov2       my($orfproti, $prostart5i, $proend3i, $bestorfi, $utrorfi)= getBestProt($fullpart, $cdnain); # , undef, $oldStart_b,$oldStart_e
#ov2 
#ov2       my $change=0;
#ov2       ($change,$orfprot,$bestorf) = bestorf_test($bestorf,$bestorfi);
#ov2       if($change) { 
#ov2         $isrev= $rev; $cdnabest= $cdnain; $utrorf= $utrorfi;
#ov2         ( undef, $prostart5, $proend3) = orfParts($bestorf);
#ov2       }
#ov2       
#ov2       }
#ov2     } # MINTR
#ov2  
#ov2     my $crev=($isrev==1)?"-":"+";  
#ov2     if($orfprot) {
#ov2       $orflen= $bestorf->{length};
#ov2       $pcds  = ($clen>0 && $orflen>0) ? int(100*$orflen/$clen) : 0;
#ov2       my $u1len= $prostart5 - 1; my $u2len= $clen - $proend3;
#ov2       
#ov2       if($NoStopCodon and substr($orfprot,-1) eq '*') { $orfprot =~ s/\*$//; }
#ov2       $aalen= length($orfprot); # $bestorf->{length}; # this is cds-len
#ov2       $aalen-- if(substr($orfprot,-1) eq '*'); # annoyance
#ov2       $compl= $bestorf->{complete};
#ov2       $compl= ($compl==3)?"complete":($compl==2)?"partial5":($compl==1)?"partial3":"partial";
#ov2       ##? not bad if partial? if u1len or u2len == 0
#ov2       if($pcds < 35 or $u1len > $orflen or $u2len > $orflen) { $compl.="-utrbad"; }  
#ov2       elsif($pcds < 60) { $compl.="-utrpoor";  } #?? maybe change to use prostart5 OR protend3 > 35%? 40% ?
#ov2       if($utrorf) { my $ol=int($utrorf->{length}/3); $compl.="-utrorf$ol"; }
#ov2       $ngood++;
#ov2     } else {
#ov2       # print join("\t", $id, $clen, 0, 0, 0, 0, "na")."\n" unless($action =~ /fasta/);
#ov2     }
#ov2     
#ov2     if($action =~ /fasta/) {
#ov2       my $cdnah= $cdnahead{$id}||"";
#ov2       if($action =~ /all/ or $orfprot) {
#ov2       $orfprot =~ s/(.{60})/$1\n/g;
#ov2       # clen=$orflen/$clen,$prostart5-$proend3; ??
#ov2       print ">$id aalen=$aalen,$pcds%,$compl; clen=$clen; strand=$crev; offs=$prostart5-$proend3; $cdnah\n";
#ov2       print $orfprot,"\n"; 
#ov2       
#ov2       # FIXME for utrorf: find way to split cdna b/n 1st, 2nd orf: midway?
#ov2       if($utrorf) {
#ov2         my( $uprot, $ustart, $uend) = orfParts($utrorf);
#ov2         my $ulen=  $utrorf->{length};  
#ov2         my $ualen= length($uprot);  
#ov2         my $upcds  = ($clen>0 && $ulen>0) ? int(100*$ulen/$clen) : 0;
#ov2         my $ucompl= $utrorf->{complete};
#ov2         $ucompl= ($ucompl==3)?"complete":($ucompl==2)?"partial5":($ucompl==1)?"partial3":"partial";
#ov2         $uprot =~ s/(.{60})/$1\n/g;
#ov2         print ">$id.utrorf aalen=$ualen,$upcds%,$ucompl; clen=$clen; strand=$crev; offs=$ustart-$uend; flag=UTRprotein\n";
#ov2         print $uprot,"\n"; 
#ov2         }
#ov2         
#ov2       if($action =~ /cds|cdna/i) {
#ov2         my $tag=($action=~/cds/i)?"cds":"cdna";
#ov2         $cdnabest= $bestorf->{sequence} if($action=~/cds/i); # not last
#ov2         $cdnabest  =~ s/(.{60})/$1\n/g;
#ov2         print ">$id.$tag aalen=$aalen,$pcds%,$compl; clen=$clen; strand=$crev; offs=$prostart5-$proend3; \n";
#ov2         print $cdnabest,"\n";
#ov2         }
#ov2       }
#ov2     } else { 
#ov2       print join("\t", $id, $aalen, $pcds, $compl, $clen, $crev, $prostart5, $proend3, $orfprot),"\n"; 
#ov2     }
#ov2 
#ov2   }
#ov2   warn "#findcds from cdnaseq found $ngood / $ncdna\n" if $debug;
#ov2 }
