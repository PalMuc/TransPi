#!/usr/bin/env perl
# asmrna_dupfilter.pl V4,2017.03

=item about
  
  Classifier of gene loci and alternates from transcript over-assembly.
  
 ..

=item usage variants

  1. basic, cds.blastn inputs cds aligns
  $evigene/scripts/rnaseq/asmrna_dupfilter.pl -debug \
    -aasize $pt.cds.qual \
    -blastab ${pt}.cds.blastn.gz \
    -outalign $pt.alntab  -outclass $pt.trclass
  
  2. add aa-homology input -ablast names, -aligntab from last -outalign has cds aligns
  $evigene/scripts/rnaseq/asmrna_dupfilter.pl -debug \
    -aasize $pt.cds.qual \
    -ablast $pt.names \
    -aligntab $pt.alntab  -outclass $pt.tr2class

  3. add input genome located cds-equivalences, -eqgene 
     from evigene/equalgene.pl -in genes.gff -over genes.gff
  $evigene/scripts/rnaseq/asmrna_dupfilter.pl -debug \
    -aasize $pt.cds.qual \
    -ablast $pt.names \
    -blastab evgm3self-${pt}.blastn.gz \
    -eqgene  evgm3chra-$pt.eqgene.gz \
    -outalign evgm3chra-$pt.tgalntab  -outclass evgm3chra-$pt.tgclass

  
  4. test, good for arabidopsis, using genome-CDS alignments (mustmap, musteqgene), 
   *after* 1st evg reduction, which left many genome-overlapped-cds loci (mostly paralogs)
    -outclass is only output (trclass table), -eqgene cds.ovself.eqgene input primary locus alignments
    also need input -ablast names table of refaa best hits, -acdhit clstr, to classify aabestho and AADUP
    input -aligntab is still required, but -mustmap/musteqgene make its aligns unused
    
  $evigene/scripts/rnaseq/asmrna_dupfilter4stub.pl -mustmap -musteqgene  -debug  \
    -aasize $pt.cds.qual  -ablast $pt.names -acdhit evg5arath_evgeqset.cd97aa.clstr \
    -eqgene evg5arath_cds.ovself.eqgenexokd  -aligntab evgm3chra-$pt.tgalntab  \
    -outclass evgm3chra-$pt.icJclass  >& evgm3chra-$pt.icJclass.log &

  convert icJclass to pubid class table (retained main, alt and no[alt]class genes)
  $evigene/scripts/prot/trclass2mainalt.pl -debug -cullx -idpre Arath5iEVm  -trclass evgm3chra-$pt.icJclass
  
    
=cut

use FindBin;
use lib ("$FindBin::Bin/../","$FindBin::Bin"); # 201405 add; assume evigene/scripts/ layout: main, prot/, rnaseq/, ...

use strict;
use warnings;
use Getopt::Long;

use constant VERSION => '2017.03.18'; # 2017 vers.4 update
  # v4: add locus reclassify from alternate input locus class table
  # v4: clean out debug, old code tests
  
# classifier options
our($AAMIN,$AAPART,$AAMINBAD, $AAMINPOO, $AADUP_IDENT,
    $MINCDS, $MINUTR, $OK_CDSUTR, $BAD_CDSUTR, $noRESET_CDSUTR, $BAD_GAPS, 
    $TINYALN, $IS_CDSALIGN,$ALTFRAG,$PHIALN,$NHIALN,$PHI,$PMID,$PLOW);

# inputs are for cds-align    
$IS_CDSALIGN= 1; # default now, was $ENV{cdsw}, using this only/mostly w/ cds-aligns, recheck for rna aligns

## protein quals
$AAMIN =$ENV{aamin}||30;   #  aacomplete, utrok
$AAPART=$ENV{aapart}||100; #  aapartial  
$AAMINBAD=$ENV{aaminbad}||60; # utrbad class
$AAMINPOO=$ENV{aaminpoo}||60; # utrpoor class
$BAD_GAPS= $ENV{aagapmax} || $ENV{BAD_GAPS} || 10; # was 25; # need opt? 'aagapmax' ?
$AADUP_IDENT=$ENV{aadupid}||98; # %ident, option for aacluster ident drops

$MINUTR = $ENV{minutr} ||300; # ~ 300b "fixed" average utr sizes, maybe too low, 
$MINCDS = $ENV{MINCDS} || 0; # reset 0 to  3*AAMIN ; drop this as separate opt?

## tr-self %identity levels to classify alt-tr
$PHI = $ENV{phi} ||99; 
$PMID= $ENV{pmid}||90; 
$PLOW= $ENV{plow}||80;  
$ALTFRAG= 0.5;
$PHIALN=$ENV{ahi} || 98; # was 99; was 65 !!;  DROP THIS; use mina?
$NHIALN=$ENV{nhi} || 9; # what? for short cds cancel few base diff < PHIALN
  
my $OK_CDSUTR_default= $ENV{pcdspoor} ||60; ##  # too high? changes w/ cds-size
my $BAD_CDSUTR_default= $ENV{pcdsbad} ||30; ## % CDS/trlen
$OK_CDSUTR= 0;
$BAD_CDSUTR= 0; # % CDS/trlen
$noRESET_CDSUTR= 1; # flag to ignore input tqual/aaqual utrbad

use constant UseTINYALNBASES => 1;   
$TINYALN = 20; # %align, was 25; NOT USED NOW ! UseTINYALNBASES
my $TINYALNBASES= $ENV{MINALIGN}||90; # was $MINALIGN= 90; # REUSE//NOT USED NOW; use just TINYALN; change to several levels
my $MIN_EQGENE= $ENV{MINEQGENE}||15; # %equalgene cds-align, was 33, 15% too high a cut? try <10% .. <5%

my $OUTSPANTAB= 1;  # make default until replace outclusters()
my $OVSLOP=6; # FIXME: need overlapslop ~ < 0.01 of max part span ; or < 0.02..0.05 of min part span?
my $pctOVSLOP=$ENV{pctover} || 0.02; # need opt?
   $pctOVSLOP=$pctOVSLOP/100 if($pctOVSLOP > 0.9);
my $SORTEDALIGNTAB=0; # debug input: sort -k2,2nr -k7,7nr -k6,6n evg2anofunz4c.alntab
my $OIDFIX=undef;

my $EQGENE_OVERRIDES_ALN = $ENV{musteqgene}||0;# identityclass
my $mustmapqual= 0; #v4:test
my $MAPMINAL= $ENV{mapminal}||66; #  eqgene mapqual: 85a,100i,2469l,4x,Spl:29%,KB668900 == %align,%ident,bases,exons,Spl:split
my $MAPMINID= $ENV{mapminid}||66;   # merge these?
my $MAPNOSPLIT= $ENV{mapnosplit}||0; #depends on inputs, have prefiltered bad splits, need to erase tinysplit flags
my $EQ2ALNTAB= 0; # make alntab from eqgene tab, only for -musteqgene -mustmap 
my $SQUISHIDS="";
my $KEEP_G= $ENV{keepg}||0; # ids_G2,3,,, are they to be dropped or kept?

my $debug= 0;
my ($aasizes,$trsizes,$blatpsl,$blastab,$lastz,$bcdhit,$aablast,$aanames,$aacdhit,$outeqtab,$outclass,
    $locuspairs, $dupids,$logfile,$head,$eqgene)= (0) x 20;
my @opts= @ARGV; #debug, does GetOptions remove @ARGV?

my $optok= GetOptions(

    # prot, transc/cds size tables
  "aasizes=s", \$aasizes, 
  "trsizes=s", \$trsizes, 
  
    # protein quality info
  "ablastab=s", \$aablast,   # this is traa-refaa.blastp
  "anames=s", \$aanames,   # variant aablast used also for naming
  "acdhit=s", \$aacdhit,    # this is traa-self.cdhit.clstr

    # V4:input variants of transcript-self-align data, only using blastn now .. drop others
  "blastab=s", \$blastab,    # this is tr-self.blastn, -CDSALIGN for cds-self.blastn
  "lastz=s", \$lastz, # lastz general format; was -blastz option -blast[ab] conflict **
  "blat=s", \$blatpsl, # tr-self.blat;  other input format opts here..?  -informat=xxx
  "bcdhit=s", \$bcdhit, #  tr-self.cdhits.clstr; 
  
  "dupids=s", \$dupids, #  
   
    # V4: eqlocus alternate table of re-classified transcripts (cures for alt-paralog problems, other)
  "locuspairs=s", \$locuspairs,    # related to eqgene input, but different, with pairwise ID locus equivalences, difference
     
    # V4: eqgene needs work,
  "eqmap|eqgene=s", \$eqgene,    # 130901: mapping equivalence table, of $evigene/equalgene.pl 
  "MINEQGENE=i", \$MIN_EQGENE,   # pct equal/align
  "EQGENEOVERRIDESALN|musteqgene!", \$EQGENE_OVERRIDES_ALN, #v4:test
  "mustmapqual!", \$mustmapqual, #v4:test
  "EQ2ALNTAB!", \$EQ2ALNTAB, #v4:test; impiles -eqgene, -musteq, -mustmap or NO -aligntab -blastab
 
  "aligntab|outalign|outeqtab=s", \$outeqtab,  ## intermediate output pairwise aligns, from blast/other align,  
  "outclass=s", \$outclass,   # primary output table of classified transcripts
  "sortedaligntab!", \$SORTEDALIGNTAB, 
  "KEEP_G|keepg!", \$KEEP_G, 


  "CDSALIGN!", \$IS_CDSALIGN, 
    # various cut-off parameters
  "AAMIN=i", \$AAMIN,   
  "AAPARTMIN=i", \$AAPART,   
  "AABADMIN=i", \$AAMINBAD,   
  "AAPOOMIN=i", \$AAMINPOO,   
  "aagapmax=i", \$BAD_GAPS,   
  "pCDSOK=i", \$OK_CDSUTR,   # CDSOKUTR
  "pCDSBAD=i", \$BAD_CDSUTR, # CDSBADUTR 
  "ALTFRAG|fragment=s", \$ALTFRAG,  # pALTFRAG ?
  "OUTSPANTAB!", \$OUTSPANTAB, # now fixed=1 ? drop opt
  "TINYALN|MINALIGN=i", \$TINYALNBASES,  # for UseTINYALNBASES
  "SQUISHIDS=s", \$SQUISHIDS,  # mem saver, squish all input IDs by pattern
  "debug!", \$debug, 
  );
  #  "logfile=s", \$logfile,  ## add??
  #old: "overlocus!", \$DOOVERLOCUS,  # V4:drop
  # "TINYALN|MINALIGN=i", \$TINYALN,  # obsolete as %align, but keep min align bases for UseTINYALNBASES

my $hasinalign= ($blastab) ? 1 : 0;
#v4o: my $hasinalign= ($blastab or $lastz or $blatpsl or $bcdhit) ? 1 : 0;
my $hasaligntab= ($outeqtab and -s $outeqtab);
$EQ2ALNTAB=0 unless($eqgene and -s $eqgene);
$hasinalign=1 if($EQ2ALNTAB and not ($hasinalign or $hasaligntab));

warn "# EvidentialGene asmrna_dupfilter.pl VERSION ",VERSION,"\n" if($debug); # change to loggit() ?
warn "# CMD: $0 @opts\n" if($debug);

die  "usage:  asmrna_dupfilter.pl -aasize=name.cds.qual 
 outputs: -outclass name.trclass -outalign name.aligntab
 input cds-align: -blast self-cds.blastn | -aligntab name.aligntab (output of blastn, or input)
 input rna-align: -noCDSALIGN -blast self-rna.blastn 
 input-opt: -ablastab refaa.names|refaa-blastp.table  -acdhit selfaa-cdhit.clstr 
          : -eqgene self-genomap.eqgene
 sets: -AAMIN=$AAMIN -AAPART=$AAPART -CDSUTR=$OK_CDSUTR percents  
 .. and more options, see source.
" unless($optok and $aasizes and ($hasinalign or $hasaligntab));  
# -trsize=name.tr.count ;  | -blat=name.blatpsl | -bcdhit=name.cdhit.clstr | 

my(%aasize, %sizeval, %codepot, %trsize, %aaqual, %oids, %oidsof, %cdsoff); # globals from readSizes();  
my(%aablast, %aablastref,$naabl);  # globals for readAAblast
my(%aacluster,%aaclustermain); # globals for readAAcdhit
my(%validids, %bspans);  # readblasttab globals
my(%dupids, %dupfirst, $ndupdrop); # readDupIds globals
my(%eqflag,$eqgenes,$nineqgene,$neqgene,$gmapqual,$neqalts,$neqexon,$eqexons,$negenes); # readEqualGene
$naabl=$ndupdrop=$neqgene=0;

$noRESET_CDSUTR=($OK_CDSUTR>0 or $BAD_CDSUTR>0)?0:1; # user option set, ignore input aaqual 201402 update
$OK_CDSUTR ||=$OK_CDSUTR_default;
$BAD_CDSUTR||=$BAD_CDSUTR_default;
$ALTFRAG= $ALTFRAG/100 if($ALTFRAG>1); # prop not percent
$MINCDS= 3 * $AAMIN if($MINCDS < 3);
$TINYALN= 100*$TINYALN if($TINYALN>0 and $TINYALN<1);  # pct not prop; obsolete for TINYALNBASES
my $TINYALNCDSLEN= 3*$TINYALNBASES;

my $OUTH= *STDOUT;
my $tmpdir=$ENV{TMPDIR}; $ENV{TMPDIR}= './' unless($tmpdir and $tmpdir ne '/tmp');  # sort TMPDIR
$hasinalign=0 if( ($ENV{useouteqtab} or $ENV{usealigntab} ) and -s $outeqtab ); # dont remake eqtab

## stubs, later update
# sub readlastz {}
# sub readcdhit {}
# sub readblatpsl {}
#--------------------------------------------------  


sub MAINstub {}

# sizes, id equivalences
my($naasizes,$ntrsizes)= readSizes($aasizes,$trsizes);
my($nfirstidsindupset,$ndupids)= ($dupids)? readDupIds($dupids): (0,0);

# protein qualities
if($aablast and $aablast =~ /\.names/ and not $aanames) { $aanames= $aablast; $aablast=""; }
($naabl)= readAAnametab($aanames) if($aanames); # prefer this now?
($naabl)= readAAblast($aablast)   if($aablast and not $naabl);  # one or other of aablast,aanames

readAAcdhit($aacdhit) if($aacdhit); # also correctAAcluster()

# V4:drop?/change EQGENE_OVERRIDES_ALN EQGENE_CHANGES_NOALN DOOVERLOCUS 

if($EQ2ALNTAB) {
  unless($outeqtab) { ($outeqtab=$eqgene) =~ s/\.\w+$//; $outeqtab.=".eqalntab"; } # make file name
  if($outeqtab) {
    rename($outeqtab,"$outeqtab.old") if( -f $outeqtab );
    open(OUTH,">$outeqtab") or die "write $outeqtab";
    $OUTH= *OUTH;
  } 
  $hasinalign=0;
}

($eqgenes,$nineqgene,$neqgene,$gmapqual,$neqalts,$neqexon,$eqexons,$negenes) # drop negenes for neqalts
    = ($eqgene) ? readEqualGene($eqgene) : (undef,0,0,0,0,0);
# $EQGENE_OVERRIDES_ALN=0 unless($neqgene);
if($EQ2ALNTAB and $outeqtab) { close($OUTH); $OUTH=undef; } # *STDOUT ?

# V4: add new reloc input option
#? readLocusPairs($locuspairs)   if($locuspairs); #? use EqualGene format, tables?

if($hasinalign) {

  if($outeqtab) {
    rename($outeqtab,"$outeqtab.old") if( -f $outeqtab );
    open(OUTH,">$outeqtab") or die "write $outeqtab";
    $OUTH= *OUTH;
  }

  my $nbalign=0;
  if($blastab) { ($nbalign)= readblasttab($blastab); }
  # elsif($lastz) { ($nbalign)= readlastz($lastz); }    #V4: maybe drop alt align inputs, not updated
  # elsif($bcdhit) { ($nbalign)= readcdhit($bcdhit); }
  # elsif($blatpsl) { ($nbalign)= readblatpsl($blatpsl); } # detect from input table ??
  if($outeqtab) { close($OUTH); $OUTH=undef; } # *STDOUT ?
  warn "# readalign: nids=$nbalign to $outeqtab\n" if $debug;
}

if($outeqtab) {
  my $infile=$outeqtab;    
  my $insorted=$SORTEDALIGNTAB; 
  unless($infile and -f $infile) { 
    die "ERR: missing input align table -aligntab $outeqtab";
  } elsif($infile =~ /\.gz/) { 
    die "ERR: cant use gzipped align table -aligntab $outeqtab";
  }

  if($aacdhit) {
    my $havevalid= scalar(%validids)?1:0;
    $havevalid= readIdsFromAlnTab($infile, $ndupids>0) unless( $havevalid); 
    correctAAcluster($havevalid); # update %aacluster,%aaclustermain
  }

  $OUTH= *STDOUT; 
  if($outclass) {
    rename($outclass,"$outclass.old") if( -f $outclass );
    open(OUTC,">$outclass") or die "write $outclass";
    $OUTH= *OUTC;
  }
  
  # primary function
  identityclass($OUTH,$infile,$insorted);  
 
  if($outclass) { close($OUTH); $OUTH=undef; }  
}

#===========================

use constant { kAATINY => 1, kAAUTRBAD => 2, kAAUTRPOOR => 4, kAADUP => 8, kAAGAPS => 16, 
              kNONCODE => 32,  kBADMAP => 64, }; 
              # kNONCODE == cds.qual Noncode flag ; kDUPEXONS => 64 ?? for cullExonEq
              # V4 add: kBADMAP => 64  
use constant NOTPOORBAD => kAATINY + kAADUP + kAAGAPS + kNONCODE + kBADMAP; # - kAAUTRPOOR - kAAUTRBAD - kAAGAPS
use constant NOTPOOR => kAATINY + kAADUP + kAAUTRBAD + kAAGAPS + kNONCODE + kBADMAP; # - kAAUTRPOOR

sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }
sub _bint { my $b=shift; return ($b<0)?0 : ($b=~/e\+/)? int($b) : $b; }
sub _squishID { if($SQUISHIDS) { map{ s/$SQUISHIDS//; } @_; } return @_; } 

sub openRead {  
  my($fna, $nostdin)= @_; my($ok,$hin)= (0,undef);
  $ok= ($fna =~ /\.gz$/) ? open($hin,"gunzip -c $fna|") 
  	 : ($fna =~ /stdin|^-/ and not $nostdin) ? *STDIN 
  	 : open($hin,$fna);  
	die "ERROR: openRead $fna" unless($ok); # loggit(1,"ERR: openRead $fna") 
  return ($ok,$hin);
}


=item identityclass

  classifier of transcript locus/alternate by alignments, other qualities

  OutaMem problems w/ largish gene sets .. on 8GB dgg.box :( 2017.03
  #OM  = reduce memory, block unneeded hashes 
  .. also try shortening all ID prefixes by pattern, on input
  
=cut

sub identityclass 
{
  my($outh, $infile, $insorted)= @_;

  #v4o# use constant TEST3 => 1; # 13aug09 test fixes to alt/main classing
  #v4o# use constant TEST1602 => 1; ## UPD 2016.02 problem w/ dup equal mains, equivalence after 1st see both as mains, 1 to be alt,
  #v4o# use constant TEST1603 => 1; ## test use here eqgenes/eqexons??
  
  #x use constant EQGENE_OVERRIDES_ALN => 0;
  #a my $EQGENE_OVERRIDES_ALN = $ENV{musteqgene} || $ENV{eqgenemust} ||0;# identityclass

  ### infile is putspans() alntab: Qid Qclen Qtlen Tid Tclen Ttlen Align Ident Bits
  ###  sort ^Qclen, ^Align, vTtlen: cat $atab | grep -v '^Qid' | sort -k2,2nr -k7,7nr -k6,6n 
  ### should it be ^Qclen, ^QID, ^Align, vTtlen : to keep qids together? no want top Aligns first

  my $ALNSORTORD='-k2,2nr -k7,7nr -k3,3n -k6,6n -k1,1 -k4,4';  #add vQtlen/3 before vTt/6 IDs to order ties
  unless($IS_CDSALIGN) { # sort tlen, not clen 
    # .. input blast align k7 is for tlen, not clen, but want to choose best by long clen > long tlen,talign
    $ALNSORTORD='-k2,2nr -k7,7nr -k3,3n -k6,6n -k1,1 -k4,4';  # test same as IS_CDSALIGN  
  }
  
  my($inh, %class,%bestmatch,%ismain, %havepair); 
  if(ref($infile)) { $inh= $infile; } # STDIN,.. sort ??
  else {
    if($insorted) { open(IN,$infile) or die "reading $infile"; $inh=*IN; }
    else { open(IN,"sort $ALNSORTORD $infile |") or die "reading $infile"; $inh=*IN; }
  }
  
  my($lastd)=("");
  while(<$inh>) { # maybe local table, sorted, or from file
    next if(/^Qid|^\W/); chomp; 
    my @v= split"\t"; 
    my($qd,$qc,$qw,$td,$tc,$tw,$aln,$iden,$bits)= @v; 
    my($isfrag,$aclass,$alnmax,$pid,$pal,$samesize)= (0) x 10;
    ($qd,$td)= _squishID($qd,$td);
    
    my $antiflag= ($IS_CDSALIGN and $bits < 0) ? "/-sense" : "/."; ## append to bestmatch ??
    		## bestmatch="id,99/89/-sense" for antisense?
    
    $isfrag= $aclass="";
    $samesize=($tc == $qc)?1:0; # OR tinyalndiff? abs($tc - $qc) < $NHIALN
    if($tc > $qc){ ($qd,$td)=($td,$qd); ($qc,$qw,$tc,$tw)=($tc,$tw,$qc,$qw); } # swap to greater cds?    		
		
		## V4.BUG: mustmapqual is dropping true main/noclass due to missing alt-mapqual (qd = alt?, td=main)
		## need it to set %ismain or bestmatch, to pull out at end .. use other means?
    #NOT HERE?? but may need skip qd to break badjoin loci .. cant use badmap as main
    if($mustmapqual) { next unless($gmapqual->{$qd}); } # V4:test
    
		$lastd=$qd; # here?
		
    if($td eq "self" or $qd eq $td) { # putspan fix for no-align cases
      $td=$qd; # $aclass="noclass"; << should be this, but ??
      $bestmatch{$qd}="$qd,100/100/self1" unless($bestmatch{$qd});
      next;   
    } 
		#? fixme: no longer have self align in align.tab, always add?
		# dont need this as bestm placeholder if below always sets $bestmatch{$qd} to something
    # unless($bestmatch{$qd}) { $bestmatch{$qd}="$qd,100/100"; } 

    my($qsize,$tsize)= ($IS_CDSALIGN) ? ($qc,$tc) : ($qw,$tw);
    # note: alnmax is min(qsize,tsize) not max
    $alnmax= ($qsize>$tsize and $tsize>0)?$tsize:($qsize>0)?$qsize:$tsize;  
    $isfrag= ($tsize < $ALTFRAG*$qsize)?"frag":"";
        
    $pid= ($aln<1)?0: int(0.5+ 100*$iden/$aln); $pid=100 if($pid>100);
    $pal= ($alnmax<1)?0 : int(0.5+ 100*$aln/$alnmax); $pal=100 if($pal>100);
    #o my $palq= ($qsize<1)?0 : int(0.5+ 100*$aln/$qsize); $palq=100 if($palq>100);
    #o my $palt= ($tsize<1)?0 : int(0.5+ 100*$aln/$tsize); $palt=100 if($palt>100);
    
    my $tinyalndiff= ((($alnmax - $aln) < $NHIALN)) ? 1 : 0; # TEST3 add, use w/ $pal

		if($neqgene>0) { 
			my $palmap= $eqgenes->{$qd}{$td}||0; # note this is pct of td aligned to qd
			my $nomap = $neqalts->{$qd}{$td}||0; # neqalts or negenes ??

      if($mustmapqual) { $nomap=1 unless($gmapqual->{$td}); } # V4:test ?? yes or no?
      #?? need an ignore-eqgene palmap/nomap case, when ichain/eqgene loci are bad?
      
      ##V4: for -mustmap -musteqgene, these become only valid alignments: 
      ## palmap<0 or nomap require no locus association of qd,td, i.e. no bestmatch, etc below
      
if($EQGENE_OVERRIDES_ALN) { ##? == musteqgene;  old test, redo
			if($nomap) { 
			  $antiflag .="/altpar0.no"; 
        #OM $eqflag{$td}{$qd}="altpar0.no.o$pal"; # for report
			  $pid=1; $pal=1; $aln=1; # $pid below PLOW == not same locus
## bug here?? not breaking nomap into 2 loci?? try as if no pair, next
        if($mustmapqual) {
        unless($bestmatch{$qd}) { $bestmatch{$qd}="$qd,100/100/self1"; } 
			  next;
			  }	
			  		  
			} elsif($palmap > 0) { 
			  $antiflag .="/altmap$palmap"; 
        $eqflag{$td}{$qd}="altmap$palmap.o$pal"; # for report
  	    $pid=$PHI; $pal=$palmap; #?? $aln= $alnmax * $palmap/100; # always replace regardless of pal?

			} else { # no palmap implies nomap for mustmap + musteqgene ??
## bug here?? not breaking nomap into 2 loci?? try as if no pair, next
        if($mustmapqual) {
        unless($bestmatch{$qd}) { $bestmatch{$qd}="$qd,100/100/self1"; } 
			  next;
			  }
			}
			
}	else {		
      # if($nomap) { } ..
      if($palmap>0) { 
			  my $xeq= $eqexons->{$qd}{$td}||0; # only for ($neqexon>0)
			  my $XPHI = 95; my $XPLO= 3;
			  # xeq care about (a) xeq>= identity == not alt but redundant, 
			  #   (b) xeq <= noxover, not alt but paralog maybe, (c) middle = usual alt
        
        if($neqexon>0 and $xeq >= $XPHI) { 
          if($palmap < $PHI) { } #?? no change or yes
          $antiflag .="/altmap$palmap"."xeq"; 
          $eqflag{$td}{$qd}="altmapxe$palmap.$xeq"; # for report, which way?
          $pal=$palmap if($palmap >= $PHI and $pal < $PMID);  # FIXME bad annot
          }
        elsif($neqexon>0 and $xeq < $XPLO) { # problems here.. dont reset pal yet
          $antiflag .="/altparx$pal";  #later? $pal=3; $aln=13;
          $eqflag{$td}{$qd}="altparx$palmap.$xeq"; # for report
          }
        else { 
          ## skip this case for tid =~ qid already alts of same locus, dont really need now?
          (my $qg=$qd)=~s/t\d+$//;  unless($td=~/^$qg/) {
          $antiflag .="/altmap$palmap";
          $eqflag{$td}{$qd}="altmap$palmap.$xeq"; # for report
          $pal=$palmap if($palmap >= $PHI and $pal < $PMID);  # FIXME bad annot
          }
          }
        }
      elsif($pal>0 and exists( $eqgenes->{$td}) and exists( $eqgenes->{$qd})) { # defined or exists??
        ## bad here for nomap alts ** why isnt defined/exists working? for ids not in eqgene
        if($gmapqual->{$qd} and $gmapqual->{$td}) {
          my $revmap= $eqgenes->{$td}{$qd}||0; 
          ## unmapped gene cases are problem here.. need what? defined ($eqgenes->{td}) 
          $antiflag .="/altpar$pal" if($revmap<3);  # $pal=3; $aln=13; #? change later?
          $eqflag{$td}{$qd}="altpar$pal.$palmap.$revmap"; # for report
        } else {
          #lots# warn "#DBG: bad eqgene $qd,$td\n" if($debug);
        }
        }
} # EQGENE_OVERRIDES_ALN        
		}

		#OM $havepair{$qd}{$td}= $pal;  #?? unused??
    my $pidalnval="$pid/$pal$antiflag"; # old: $pid/$pal
    my $tisaltpar=0; # treat like $skiptinyaln for now
       $tisaltpar=1 if($pidalnval=~/altpar\d/); #?  and $qclass // not altparx\d
     
		##V4? fixme: no longer have self align in align.tab, always add?
		## dont need this as bestm placeholder if always set $bestmatch{$qd} to something, ie no next above
    unless($tisaltpar) { #??   #? block bestmatch or not? affects final out class table
    $bestmatch{$qd}="$td,$pidalnval" unless($bestmatch{$qd} and not($bestmatch{$qd}=~/self|altpar\d/));
    $bestmatch{$td}="$qd,$pidalnval" unless($bestmatch{$td} and not($bestmatch{$td}=~/self|altpar\d/)); 
    }
    unless($bestmatch{$qd}) { $bestmatch{$qd}="$qd,100/100/self1"; } 
    
    my $qclass= $class{$qd}||"";
    if($samesize and $qclass =~ /alt/ and $ismain{$td}) { # check/stop alt1 <> alt2 assigns
      next if($bestmatch{$qd} =~ /$td,/); # problem here? for many equal bestmatch
    }

    my $skiptinyaln=0;
    # if(UseTINYALNBASES)  
    if($alnmax >= $MINCDS) { # recall alnmax is min(qsize,tsize); never skip tiny prots, MINCDS =~ 90 bases ?
      my $minbaseover= ($alnmax >= $TINYALNCDSLEN)? $TINYALNBASES : $alnmax * 0.10; ##$MIN_PCTOVERLAP;
      $skiptinyaln=($aln < $minbaseover)?1:0;
    } else {
      $isfrag= "frag" unless($isfrag); #or "tiny" ? # dont skip assign alt/frag to tiny cds, but maybe always set isfrag ?
    }
  
    # set alt aclass here, or not.
    if( $skiptinyaln or $tisaltpar) { } ## defer: $aclass="noalign$isfrag";  
    elsif( $tc == 0 and $qc > 0 and $pid >= $PLOW ) { $aclass="frag0aa"; } # tc==0 special case "frag0"
    elsif( $pid >= $PHI ) { 
      $aclass= ($isfrag)?"parthi":"althi"; # Primary classifier; change parthi to althipart ? althi$isfrag ?
      $aclass .= "1" if($pal >= $PHIALN or $tinyalndiff); # althi1 == hi-align + hi-ident, nearly same (identicals gone here)
      }
    elsif( $pid >= $PMID ) { $aclass="altmid$isfrag"; }
    elsif( $pid >= $PLOW ) { $aclass="altlo$isfrag"; }

    if(my $ocl= $class{$td}) {  
      if($ocl =~ m/^alt|^part/) { $class{$td}= $aclass if($aclass eq "althi1"); next; } # test3
      elsif($ocl =~ m/^main/) { 
        if($qclass =~ m/^main/ or $tisaltpar) { 
          # skip to below aclass set, test5, yes, do this also?
          } 
        else { next; } # * YES, need this : test4
      } 
    }

    if($aclass) { 
      my $qmain=$qd; my $attr= $pidalnval; ## "$pid/$pal$antiflag";

      use constant FINDMAINID => 1;  
      if($class{$qd}) {  #  =~ /alt|part/
        my $qnext= $qd; 
        my $more= ($bestmatch{$qnext})?1:0;
        my %qdid=( $qd => 1, $td => 1);
        while($more) { 
          $more=0;
          my($qbest,$qpal)=split",",$bestmatch{$qnext};
          
          $qbest=0 if($ismain{$qnext} and not $ismain{$qbest}); # want always? maybe bad?
          
          if($qbest and not $qdid{$qbest}) {
            $qnext= $qbest; $qdid{$qbest}=1; 
            $more=1 if($class{$qnext} and $bestmatch{$qnext});
          }
        }
        if($qnext ne $qd) { $qmain= $qnext; $attr="$pidalnval/$qmain"; } #was $attr.="/$qmain";
        ## WARNING: now attr has /-sense or /. before /qmain; 
        ##   evgmrna2tsa2.pl needs to know this field's structure, /qmain at end esp.
      }
      
      $ismain{$qmain}++;   
      $class{$td}=$aclass; #? problem here when class{td} == main set before, should not reset to alt
  
      # this way works.  if($bestmatch{$td}) must be true from above
      if($bestmatch{$td} =~ /^$qd/) { $bestmatch{$td}="$qd,$attr"; }
      elsif(not $bestmatch{$td}) { $bestmatch{$td}="$qd,$attr"; }

      ## add to prevent 2+ circular alta <> altb <> altc with no main ..
      $class{$qmain}="main" unless($class{$qmain}); # should always be:  $qc >= $tc//$samesize..

    }
      
  } close($inh);

  
  # END:  print trclass table; add more fields to output: aaqual, aablast, tr,aa sizes?
  # V4 FIXME: cant rely on bestmatch/ismain for cases with no best match, other alignment ?
  { my($q,$pal,$c,$d);
  foreach $d (sort keys %class) {
    ($q,$pal)=split",",$bestmatch{$d}; $c= $class{$d}; 
    my @cla= classifytr( $d, $c, $q, $pal);
    print $outh join("\t",@cla)."\n"; 
    }
  foreach $d (sort grep{ not($class{$_}) } keys %ismain) {
    ($q,$pal)=split",",$bestmatch{$d}; $c= "main";    
    my @cla= classifytr( $d, $c, $q, $pal);
    print $outh join("\t",@cla)."\n"; 
    }
  foreach $d (sort grep{ not($class{$_} or $ismain{$_}) } keys %bestmatch) {
    ($q,$pal)=split",",$bestmatch{$d}; $c= "noclass";  
    my @cla= classifytr( $d, $c, $q, $pal);
    print $outh join("\t",@cla)."\n"; 
    }
  }  
  
}


=item classifytr

  classifier of locus primary, alternate, fragment, redundant 
  using identityclass() collection of overlapping transcripts (CDS or full tr)  
  focused on CDS qualities, maybe should be classifyCDS() ..
  
  see also updated subs trevidence($trid,xxx) + classifyFullTr($tid,$cla,$qid,$pidal)
    -- classifyTrFull  adds tr overlaps outside(?) of CDS x CDS:

=cut


sub classifytr {
  my($tid,$cla,$qid,$pidal)= @_;

  #v4f use constant CODEPOT1607 => 1; # add from above codepot work
  #v4f use constant ALTPOOR1607 => 1;

  # use constant { kAATINY => 1, kAAUTRBAD => 2, kAAUTRPOOR => 4, kAADUP => 8, kAAGAPS => 16, };
  # use constant NOTPOORBAD => kAATINY + kAADUP + kAAGAPS; # - kAAUTRPOOR - kAAUTRBAD - kAAGAPS
  # use constant NOTPOOR => kAATINY + kAADUP + kAAUTRBAD + kAAGAPS; # - kAAUTRPOOR
  # add: kBADMAP => 64 ?
  
  $qid||="0"; $pidal||="0";
  my($pidi,$pali)= $pidal =~ m/^(\d+).(\d+)/;
  
    # seq sizes, cds/utr ratio
  my $aw= $aasize{$tid} || 0;
  my $tqual= $aaqual{$tid} || ""; #? parse for "aasize,pcds,aaqual"
  my $tw= $trsize{$tid} || 1; 
  my $pcds= int(300*$aw/$tw);  
  my $butr= ($tw - 3*$aw); # okutr if butr<=300p, modify BAD_CDSUTR for small aw

  #?? use this: my $qv= aaqualscore($aaqual{$id}); range is (partial/utrbad) -3..+2 (complete,good)
  my $ispoor= ($aw < $AAMIN or ($tqual =~ m/partial/ and $aw < $AAPART))?kAATINY:0;

    # aa homology ref, self 
  my $tbits= $aablast{$tid} || "0,0";
  my($tbscore,$tbref)=split",",$tbits;
  
  my($tbrscore,$tbrefbest)= ($tbref) ? (split",", $aablastref{$tbref}) : (0,0); 
  my $rbits2= $aablastref{$tid} || "0,0";
  my($tbrscore2,$tbref2)= split",", $rbits2; # this is revhash tid => score,ref refbest

  my $aaclus= $aacluster{$tid} || "0,0";  
  my($aamainid,$aamainpi)= split",", $aaclus; # ="$mainid,$pi";
  if($aamainid eq $tid) {  $aamainid=""; $aamainpi=0; }
  
    # genome map qual
  my $mapqloc= ($nineqgene) ? $gmapqual->{$tid} : 0;  
  my($mapqual,$maploc)= ($mapqloc)?(split"\t",$mapqloc):(0,0);

  if($tqual =~ m/gapbad/) { $ispoor |= kAAGAPS; }  # gaps discrepancy w/ aaqual utrbad and pcds utrbad 
  if($butr >= $MINUTR) { # ~ 300b "fixed" average utr sizes
    if($tqual =~ m/utrbad/ and $noRESET_CDSUTR) { $ispoor |= kAAUTRBAD; } 
    elsif($tqual =~ m/utrpoor/ and $noRESET_CDSUTR) { $ispoor |= kAAUTRPOOR; }
    else { $ispoor |= ($pcds < $BAD_CDSUTR)? kAAUTRBAD: ($pcds < $OK_CDSUTR)? kAAUTRPOOR : 0; }
  }

  #OM my $tcode= $sizeval{$tid}{'codepot'}||""; # $tcode =~ /^Noncode/ set what? kAATINY+kAAUTRBAD
  my $tcode= $codepot{$tid}||""; # $tcode =~ /^Noncode/ set what? kAATINY+kAAUTRBAD
  $tcode=0 if($tbscore>0);
  $ispoor |= kNONCODE if($tcode and $tcode =~ /^Noncode/);  
  #?? want to keep NONCODE now if has introns? other value?
  
    # # eqgene mapqual: 85a,100i,2469l,4x,Spl:29%,KB668900 == %align,%ident,bases,exons,Spl:split
    # V4: special location flags from icchains.eqgene: 
    #   l4630.t8.icsub == lower qual alt, keep only if other quals good
    #   l4630.t7.icalt.icdup == duplicate alt, keep only if other quals good
    # V4: use icchain alt nums to reclass althi,altother .. all uniq icalt chain should be kept
  $ispoor |= kBADMAP if($mustmapqual and not $maploc); # V4:test

  
  $ispoor |= kBADMAP if($MAPNOSPLIT and $mapqual =~ /Spl:/); # << FIX, respect pMAXSPLIT or optional ignore Spl
  # Arath5EVm006792t1	drop	main	Arath5EVm006792t4	99/35/./altmap35	495,84%,complete	
  #  aaref:899,AT4G06534.1,refok,chrmap:98a,100i,1488l,2x,Spl:2%,chr4:3356788-3356814:.,chr4:3356810-3361255:+,pflag:64

  if( my($mal,$mid)= ($mapqual =~ m/(\d+)a,(\d+)i/) ) { 
    $ispoor |= kBADMAP if($mal < $MAPMINAL or $mid < $MAPMINID); }
  ## not icsub yet, need fuller ichain alt parsing.
  
  if($aamainpi >= $AADUP_IDENT) { # TEST aacluster added info on dup prots; AADUP_IDENT=98 default
    $cla.="a2" unless($cla=~/hi[012]/); # hi1/a2 redundant info
    # $tbits.=",aadup:$aamainid"; # defer, AFTER refbest/good/..
    $ispoor |= kAADUP if($cla =~ /althi/); # only for /althi|parthi/
  }
  unless( $ispoor & kAADUP ) {
    if($pidal =~ m/altmap\d+xeq/ and $cla =~ /althi1|part/) {
      $ispoor |= kAADUP;  # add new bitflag kDUPEXONS ?
    }
    # V4: eqgene altmap100 fixup, these have identical CDS set..
    if($pidal =~ m/altmap(\d+)/) { my $am=$1; 
      if($am >= $AADUP_IDENT and $cla =~ /althi1|part/) {
        $cla.="a2" unless($cla=~/hi[012]/); # hi1/a2 redundant info
        #? $tbits.=",aadup:$aamainid";
        $ispoor |= kAADUP;  # add new bitflag kDUPEXONS ?
      }
    }
    if(($maploc =~ /\.icdup/) and ($cla =~ /alt|part/)) { $ispoor |= kAADUP; }  #? kDUPEXONS
  }
  my $aadupflag= ($ispoor & kAADUP and $aamainid) ? ",aadup:$aamainid" : ""; # defer, AFTER refbest/good/..
  # below add to  $tbits.= $aadupflag
  
  my $MINBLASTSCORE= 60; # aablast only? bitscore always?
  my $keepdrop="";
  unless( $tbscore == 0 or $tbits=~/^0,0/) {
    my $risbest= (($tbrefbest eq $tid) or ($tbscore >= $tbrscore))?1:0;# was ($tbrefbest eq $tid), allow 2+ same score
    my $risgood= ($tbrscore > $MINBLASTSCORE and $tbscore >= 0.90 * $tbrscore)?1:0;
    my $risok= ($tbref2 and $tbrscore2 >= $MINBLASTSCORE)?1:0;
    if($risok and $risgood) { $risgood=0 if($tbrscore2 > $tbrscore); } #?
    my $isaadup=($ispoor & kAADUP);
    my $isbadmap=($ispoor & kBADMAP);
    # in absense of aacluster, use tbrefbest/tbrscore to set isaadup ? need new hash to track aascores?
    
    if($risbest) { 
      $tbits.=",refbest"; $keepdrop.="okay:refbest,"; $ispoor=0 unless($isaadup); }
    elsif($risgood) { 
      $tbits.=",refgood"; unless($isaadup){ $keepdrop.="okay:refgood,"; $ispoor=0;} }  
    elsif($risok) {  
      $tbits.=",refok"; unless($isaadup or $isbadmap){ $keepdrop.="okay:refok,"; $ispoor=0;} }
  }

  if($ispoor > kAATINY and $cla !~ /althi|parthi/) {  
    if($aw >= $AAMINBAD) { $ispoor = $ispoor & NOTPOORBAD; } ##  ^ (kAAGAPS+kAAUTRPOOR+kAAUTRBAD);  ##?? ^ XOR bad
    elsif($aw >= $AAMINPOO) { $ispoor = $ispoor & NOTPOOR; } ##  ^ (0+kAAUTRPOOR);
  }

  if($cla =~ /parthi|frag0aa/) {
    $keepdrop.= "drop";
    
  } elsif($cla =~ /althi/) {
    $keepdrop.= (($ispoor & NOTPOORBAD) or ($ispoor and $pali > 75))?"drop":"okay";

  } elsif($cla =~ /main/) {
    $keepdrop.= ($ispoor)?"drop":"okay";  # ??
    
  } elsif($cla =~ /noclass/) {
    $keepdrop.= ($ispoor)?"drop":"okay"; #? dont drop? ditto to altmid
    
  } else { # other altmid/low 
    $keepdrop.= ($ispoor)?"drop":"okay"; #? dont drop this class, ispoor maybe uniq prot: maybeok?
  }
  
  my $okay;
  if($keepdrop =~ /drop/) {
    if($keepdrop =~ /okay/) { $okay= "maybeok"; } else { $okay= "drop"; }
  } else {
    $okay= "okay";
  }
  
  $tbits.= $aadupflag if($aadupflag); # defer, AFTER refbest/good/..
  $tbits.= ",chrmap:$mapqual,$maploc" if($maploc and $mapqual);
  $tbits.= ",pflag:$ispoor"; # DEBUG info, but keep
    $tbits="aaref:$tbits" unless($tbits=~/^0,0/);
  
  if(defined $eqflag{$tid}) { # TEST1603
    my $eqfl= $eqflag{$tid}{$qid}||""; 
    if($eqfl) { $eqfl="$qid/$eqfl,"; }
    my @q= grep{ $_ ne $qid } sort keys %{$eqflag{$tid}}; 
    $eqfl .= join",",map{ "$_/".$eqflag{$tid}{$_} }@q;  
    $tbits.= ",feq:$eqfl";
  }
  
  return (wantarray) ? ($tid,$okay,$cla,$qid,$pidal,$tqual,$tbits) : $okay;
}


#----------------------------------------------------------
# align inputs

sub putspans {
  my($lq)= @_;
  my $nmatch=0;
  foreach my $lt (sort keys %bspans) {
    next if($lt eq $lq); # is lt eq lq allowed here?
    
    my @bspans= @{$bspans{$lt}};
    my ($tbit,$taln,$tidn,$aq,$at,$wq,$wt,$torient,$xbm,$xem,$tbm,$tem)= (0) x 19;
    foreach my $sp (@bspans) {
      my($xb,$xe,$tb,$te,$xbit,$aln,$aident,$or)= @$sp; # 2013.aug: IS_CDSALIGN add $or
      $tbit += $xbit; $taln+= $aln; $tidn+= $aident; 
      $torient += $aln * $or;  
      $xbm=$xb if($xbm==0 or $xb<$xbm); $xem=$xe if($xe>$xem);
      $tbm=$tb if($tbm==0 or $tb<$tbm); $tem=$te if($te>$tem);
      }

    $aq= $aasize{$lq}||0; $aq *=3;  # 3*convert to cds-size
    $at= $aasize{$lt}||0; $at *=3;
    $wq= $trsize{$lq}||0; 
    $wt= $trsize{$lt}||0;  
    my $alnmax= _max(1, ($IS_CDSALIGN) ? _min($aq,$at) : _min($wq,$wt) );
    
    ## if UseTINYALNBASES
    if($alnmax >= $MINCDS) { # never skip tiny prots, MINCDS =~ 90 bases ?
      my $minbaseover= ($alnmax >= $TINYALNCDSLEN)? $TINYALNBASES : $alnmax * 0.10; ##$MIN_PCTOVERLAP;
      next if($taln < $minbaseover);
    } 
    
    $nmatch++;
    if($IS_CDSALIGN and $torient<0) { $tbit= -$tbit; } #? what followon problems does this cause?
    print $OUTH join("\t",qw(Qid Qclen Qtlen Tid Tclen Ttlen Align Ident Bits Qspan Tspan))."\n" unless($head++);
    print $OUTH join("\t",$lq,$aq,$wq,$lt,$at,$wt,$taln,$tidn,$tbit, "$xbm-$xem","$tbm-$tem")."\n"; 
    $validids{$lq}++; $validids{$lt}++;
  } 

  if($nmatch==0) {
    my ($tbit,$taln,$tidn,$aq,$at,$wq,$wt)= (0) x 19; 
    my $lt="self"; # or lq, or use blast-self scores?
    $aq= $aasize{$lq}||0; $aq *=3;  # 3*convert to cds-size
    $wq= $trsize{$lq}||0; 
    $at=$aq; $wt=$taln=$tidn=$tbit= $wq; # or aq
    print $OUTH join("\t",qw(Qid Qclen Qtlen Tid Tclen Ttlen Align Ident Bits Qspan Tspan))."\n" unless($head++);
    print $OUTH join("\t",$lq,$aq,$wq,$lt,$at,$wt,$taln,$tidn,$tbit, "1-$wq","1-$wt")."\n"; 
    $validids{$lq}++; $validids{$lt}++;
  }
  
  %bspans=();
  return($nmatch);  
}


sub readblasttab {
  my ($bother)= @_;
  
  my($lt,$lq,$sa,$sm,$ss,$qd,$wq,$nids)= (0) x 10;
  my %dupskipids=(); my $dupskipspan=0;
  %bspans=();
  
  my($ok,$fh)= openRead($bother);
  while(<$fh>) { 
    unless(/^\w/) { next; } 
    my @v= split; 
    my($q,$t,$bits,$pctident,$aln,$mis,$gap,@bspan)= @v[0,1,-1,2,3,4,5, 6,7,8,9]; # 6-9 =  q. start, q. end, s. start, s. end, 
    ($q,$t)= _squishID($q,$t);

    if($lq and $q ne $lq) {  
      putspans($lq) unless($dupskipspan); 
      $nids++; %bspans=(); $dupskipspan=0; 
    }

    my $dupskip=0;
    if( $dupids ) { 
      if( $dupids{$q} and $q ne $dupids{$q} ) { 
        $dupskipspan=$q; ## if($q eq $t); # flag to skip putspan
        $dupskipids{$q}++; $dupskip++;   
      } elsif( $q ne $t and $dupids{$t} and $t ne $dupids{$t}) { 
        $dupskipids{$t}++; $dupskip++; 
      }  
      if($dupskip) { $ndupdrop++; }  
    }
  
    if($t eq $q) { }  #?? want self-align in table?
    elsif($dupskip==0 and $dupskipspan == 0) { 
      $bits= _bint($bits);
      # my $aident= _max(0,$aln-$mis-$gap); # other way to calc: $aident = $pctident * $aln;
      my $aident= int(0.5 + $aln*$pctident/100); # better?
      sumblastpart( $q, $t, $bits,$aln,$aident, @bspan);  
    } 
    $lq=$q; $lt=$t;  
  } close($fh);
  
  putspans($lq) unless($dupskipspan); $nids++; $dupskipspan=0;

  if($ndupdrop>0) { $ndupdrop=scalar(keys %dupskipids); }
  warn "# readblasttab: nids=$nids; ndupdrop=$ndupdrop\n" if $debug;
  return($nids);
}

sub sumblastpart {
  my( $q, $t, $bits,$aln,$aident, $qb,$qe,$sb,$se) = @_;
  my $or=0; # 2013.aug : ORIENT problem.  if both here are reversed, or=0; if only 1, or=-1
  if($qb > $qe) { ($qb,$qe)= ($qe,$qb); $or=-1; }
  if($sb > $se) { ($sb,$se)= ($se,$sb); $or=($or<0)?0:-1; } #was $or--
  unless($bspans{$t}) { 
    $bspans{$t}=[]; push( @{$bspans{$t}}, [$qb,$qe,$sb,$se,$bits,$aln,$aident,$or]); 
    return; }
  my $ov=0;
  my $qlen=1+$qe-$qb; my $slen=1+$se-$sb;
  my $qslop= _max($OVSLOP, int($pctOVSLOP*$qlen));
  my $sslop= _max($OVSLOP, int($pctOVSLOP*$slen));
  
  foreach my $sp (@{$bspans{$t}}) {
    my($xb,$xe,$tb,$te,$xbit)= @$sp;
    if($qe < $xb or $qb > $xe) { }
    elsif($qe > $xe and $qb >= $xe - $qslop) { }
    elsif($qb < $xb and $qe <= $xb + $qslop) { }
    else { $ov=1; last; }
    if($se < $tb or $sb > $te) { }
    elsif($se > $te and $sb >= $te - $sslop) { }
    elsif($sb < $tb and $se <= $tb + $sslop) { }
    else { $ov=1; last; }
  }  
  unless($ov) { push( @{$bspans{$t}}, [$qb,$qe,$sb,$se,$bits,$aln,$aident,$or]); }
}


sub readAAnametab {
  my($aanametab)= @_;
  my($swapids,$naabl,$nblerr)= (0) x 9;
  my %bscore;
  open(F, $aanametab) or die "FAIL: read $aanametab";
  while(<F>) { next unless(/^\w/); 
    my($td,$name,$alnscore,$rd,@more)=split"\t"; 
    ($td)= _squishID($td);
    # alnscore format expected: 72%,3270/4555,3282 ;  may be '72' or '72%' only
    my($apct,$aln)= $alnscore =~ m=^(\d+)%,(\d+)\b=; # partial match ok
    unless($aln){ ($aln)= $alnscore =~ m/^(\d+)/; }  # only care about relative score to other td
    next unless($aln and $aln>0);
    my $bscore= $aln;
    $naabl++; 
    $rd =~ s/^RefID://; # or CDD: or other
    unless($OIDFIX or $naabl>9) { 
      $OIDFIX=1 if(%oids and ($oids{$td} or $oids{$rd}));
    }
    if($OIDFIX){ $td= $oids{$td}||$td; $rd= $oids{$rd}||$rd; }
    # local $^W = 0;   # no warnings; no help
    unless($aablast{$td} and $bscore{$td} > $bscore) {
      $aablast{$td}="$bscore,$rd";  $bscore{$td}= $bscore;
      unless( $aablastref{$rd} and $bscore{$rd} > $bscore ) {
        $aablastref{$rd}= "$bscore,$td";  $bscore{$rd}= $bscore;
        $aablastref{$td}= "$bscore,$rd";   # revhash?
      }  
    } else {
      unless( $aablastref{$rd} and $bscore{$rd} > $bscore ) {
        $aablastref{$rd}= "$bscore,$td";  $bscore{$rd}= $bscore;
        $aablastref{$td}= "$bscore,$rd";   # revhash?
      }  
    }
  } close(F); 
 warn "# readAAnametab: naabl=$naabl\n" if $debug;
 return($naabl); 
}

sub readAAblast  {  
  my($aablast)= @_;
  if($aablast =~ /\.names$/) { return readAAnametab($aablast); }
  my $swapids=0; my $naabl=0; my $nblerr=0;
  my %bscore;
  my $BSCORE_COL = 0; # 0=bits, 1=ident, 2=aln
  
  use constant nCHECK => 29; ## precheck format before sort of possibly very large file ..
	my($ok,$inh)= openRead($aablast,1);
  while(<$inh>) { 
    next unless(/^\w/); 
    my($td,$rd,@v)=split"\t"; $naabl++; # expect: trid refid bits ident align .. may be refid,trid instead
    ($td)= _squishID($td);

    if(%oids) { 
      $OIDFIX=1 if($oids{$td} or $oids{$rd});   # is aasize{oid} ok here
      if($OIDFIX){ $td= $oids{$td}||$td; $rd= $oids{$rd}||$rd; }
    }
    if($swapids==1) {
      ($td,$rd)=($rd,$td);
    } elsif($swapids==0) {
      if(@v >= 12 and $v[8] =~ /^\d/) { $nblerr++; } # blast.tab
      elsif($aasize{$td}) { $swapids= -1; }
      elsif($aasize{$rd}) { $swapids= 1; }
      else { $nblerr++; }
    }
    last if($naabl > nCHECK);
  } close($inh);
  if($nblerr>2 or $swapids==0) { 
    die "ERR: expect table of aablast scores: trid refid bitscore identity align ..\n"
      ." $nblerr trids from aasize not found in -aablast=$aablast\n";
  }
  
  $naabl=$nblerr= 0;  
	my $BLSORTORD=($swapids==1)?'-k3,3nr -k5,5nr -k2,2 -k1,1' : '-k3,3nr -k5,5nr -k1,1 -k2,2'; # bits > aln > refid > qid
  $ok= ($aablast =~ /\.gz$/) 
        ? open($inh,"gunzip -c $aablast | sort $BLSORTORD |") 
        : open($inh,"sort $BLSORTORD $aablast |");  
  
  while(<$inh>) { 
    next unless(/^\w/); 
    my($td,$rd,@v)=split"\t"; $naabl++; # expect: trid refid bits ident align .. may be refid,trid instead
    ($td)= _squishID($td);
    my $bscore= $v[$BSCORE_COL]; # bits, ident, algn; want choice of? ident or algn maybe better
    if($swapids==1) { ($td,$rd)=($rd,$td); }
    if($OIDFIX) { $td= $oids{$td}||$td; $rd= $oids{$rd}||$rd; }
    
    # local $^W = 0;   # no warnings; no help
    unless($aablast{$td} and $bscore{$td} >= $bscore) {
      $aablast{$td}="$bscore,$rd";  $bscore{$td}= $bscore;
      unless( $aablastref{$rd} and $bscore{$rd} >= $bscore ) {
        $aablastref{$rd}= "$bscore,$td";  $bscore{$rd}= $bscore;
        $aablastref{$td}= "$bscore,$rd";   # revhash?
      }  
    } else {
      unless( $aablastref{$rd} and $bscore{$rd} >= $bscore ) {
        $aablastref{$rd}= "$bscore,$td";  $bscore{$rd}= $bscore;
        $aablastref{$td}= "$bscore,$rd";   # revhash?
      }  
    }
    
  } close($inh); 
 warn "# readAAblast: naabl=$naabl\n" if $debug;
 return($naabl); 
}


sub readSizes {  # see also cdna_evigenesub.pm getAaQual()
  my($aasizes,$trsizes)= @_;
  my($naa,$ntr,$nerr,$ok,$inh)=(0) x 10;
  if($aasizes) {  
    my $iscds=0;
    ($ok,$inh)= openRead($aasizes,1);
    while(<$inh>) { 
      next if(/^\W/ or /^total/); 
      my @v=split; # aa.qual cols; gap is removed from alen
      my($id,$alen,$gap,$aqual,$tlen,$offs,$oids);
      my $codepot=0; #  Code/Nonc/Unknown col maybe there..
      if($v[5] =~ /^(Code|Noncode|Unknown)/) {
        ($id,$alen,$gap,$aqual,$tlen,$codepot,$offs,$oids)=@v; 
      } else {
        for my $i (0..$#v) { if($v[$i] =~ /^(Code|Noncode|Unknown)/) { $codepot=$v[$i]; splice(@v,$i,1); last; } }
        ($id,$alen,$gap,$aqual,$tlen,$offs,$oids)= @v;  # ,offs,oids may be missing. and Code/Nonc
      }
      $offs||=0; $oids||=0;
      ($id)= _squishID($id);
      
      my $cdlen=0;
      unless($alen =~ /^\d/) { $nerr++; next; }
      if($aqual) { 
        my($aqlen)= $aqual=~m/(\d+)/; 
        if($aqlen < $alen and $aqlen*3 >= $alen-3) { $cdlen=$alen; $alen= int($alen/3); $iscds++; } #cdslen not alen
        elsif($aqlen == 0 and $tlen == 0 and $iscds) { $cdlen=$alen; $alen= int($alen/3);  } # datatable bug
        else { $cdlen=$alen*3; }
        $aqual .= "-gapbad" if($gap>0 and (100*$gap/($alen+$gap) > $BAD_GAPS)); # add qual flag for gap/(alen+gap) > MAXGAP:  
        $aaqual{$id}= $aqual; 
      }
      if($tlen =~ /^\d/) { $trsize{$id}= $tlen; $ntr++; } 
      if($offs =~ /^\d/) { $offs=~s/:.$//; $cdsoff{$id}= $offs; }  
      if($oids) { # capture, for fixups?
        my @oids=grep{ /^\w/ and not($_ eq "na" or $_ eq $id) } split",",$oids;
        for my $d (@oids) { $oids{$d}=$id; $aasize{$d}=$alen unless($aasize{$d});} 
        $oidsof{$id}=$oids if(@oids);  
      }       
      $aasize{$id}=$alen; $naa++; 
      $codepot{$id}= $codepot;
      #OM : using only sizeval{id}{codepot} now, replace this hash?
      #OM @{$sizeval{$id}}{ qw(aasize cdsize gap aaqual trsize cdsoff oids codepot) }
      #OM                   = ($alen,$cdlen,$gap,$aqual,$tlen,$offs,$oids,$codepot);
      } close($inh); 
  }
  
  if($trsizes) {
  if($trsizes =~ /^aaqual/ or $ntr>0) { 
    # got above;now default
  } elsif($trsizes =~ /^aasize|^cdssize/) { # for tr == cds
    foreach my $id (keys %aasize) { $trsize{$id}= 3*$aasize{$id}; }
  } else {  
    ($ok,$inh)= openRead($trsizes,1);
    while(<$inh>) { next if(/^\W/ or /^total/); my($id,$al)=split; $trsize{$id}=$al; $ntr++; } close($inh); 
  }
  }
 
 warn "# readSizes: naa=$naa; ntr=$ntr\n" if $debug;
 return($naa,$ntr); 
}


# fix from fastanrdb all.cds > all_nr.cds >> hdr has cds-identical ids; should have prefiltered this
sub readDupIds {
  my($infile)= @_;
  my($nids,$ndups,$inh,$ok)=(0,0,undef);
  %dupids= %dupfirst=();
	($ok,$inh)= openRead($infile,1);
  # open(IN,$infile) or die "reading $infile"; $inh=*IN; 
  while(<$inh>) {  
    s/^>//; # if from fastanrdb
    my @dupids= split; # grep /\w/ or any other?
    next unless(@dupids>1); #? or record all ids?
    @dupids= _squishID(@dupids);
    my $firstid= $dupids[0];
    if($dupids{$firstid}) {
      my $nextdup= $firstid;
      $firstid= $dupids{$nextdup};
      shift @dupids;
      push @{$dupfirst{$firstid}}, @dupids;  
    } else {
      $dupfirst{$firstid}= \@dupids;  $nids++; 
    }
    map{ $dupids{$_}= $firstid } @dupids;
    $ndups += @dupids;
  } close($inh);  
  warn "# readDupIds: nfirst=$nids; ndups=$ndups\n" if $debug;
  return ($nids,$ndups);
}

# pre-read all of infile ids into validids, before idenityclass
sub readIdsFromAlnTab {
  my($infile, $hasdupids)= @_;
  my($nids,$ok,$inh)=(0,0,undef);
  ($ok,$inh)= openRead($infile,1);
  # open(IN,$infile) or die "reading $infile"; $inh=*IN; 
  while(<$inh>) {  
    next if(/^Qid|^\W/); chomp; 
    my($qd,$qc,$qw,$td,$tc,$tw,$aln,$iden,$bits)= split"\t"; 
    ($qd,$td)= _squishID($qd,$td);
    if( $hasdupids ) {
      $validids{$qd}++ unless( $dupids{$qd} and $qd ne $dupids{$qd});    
      $validids{$td}++ unless( $dupids{$td} and $td ne $dupids{$td});    
    } else {
      $validids{$qd}++; $validids{$td}++;
    }
    $nids++;
  } close($inh); # rewind if ref() ???
  return $nids;
}

# change mainaa/subaa in %aacluster, using other input info: aaqual/size-gaps, cds/tr align ids
sub aaqualscore  # in cdna_evigenesub.pm
{
  my($mqual)= @_;
  my $mqv=0; $mqual ||="missing";
  if($mqual =~ /complete/) { $mqv += 2; } elsif($mqual =~ /partial[35]/) { $mqv += 1; }
  if($mqual =~ /utrbad/) { $mqv -= 2; } elsif($mqual =~ /utrpoor/) { $mqv -= 1; }
  if($mqual =~ /gapbad/) { $mqv -= 1; } # or -2?
  return $mqv; # range is now -3..+2, was -2 .. +2
}


sub correctAAcluster {
  my($havevalidids)= @_;
  # % aaclustermain == hash{mainid}[subids]
  # % aacluster == hash{eachid} = "mainid,pctident"
  my $nreset=0; my %didmain;
  foreach my $mid (sort keys %aaclustermain) {
    next if($didmain{$mid}); # probably dont need. 
    my $mainid= $mid;
    my @cids= @ { $aaclustermain{$mid} }; # has all ids incl mainid?
    my $maw = $aasize{$mid} || 0;
    my $mqv = aaqualscore($aaqual{$mid});  #? parse for "aasize,pcds,aaqual"
    my $reset=0;  my @goodids=(); my $ninval=0; 
    if($havevalidids and not $validids{$mainid}) { $maw=0; $mqv=-9; $mainid=0;  $ninval++; }
    foreach my $id (@cids) {
      if($havevalidids and not $validids{$id}) { 
        $ninval++; next; # drop from aacluster
      }
      push @goodids, $id; 
      next if($id eq $mainid);
      my $aw= $aasize{$id} || 0; 
      my $qv= aaqualscore($aaqual{$id});
      if($aw > $maw and $qv >= $mqv - 1) {  # reset main; ignore qv if aw >>> maw ?
        ($mainid,$maw,$mqv)=($id,$aw,$qv); $reset++; 
      } elsif($qv > $mqv and ($aw > 0.98*$maw)) {
        ($mainid,$maw,$mqv)=($id,$aw,$qv); $reset++;       
      }
    }
    
    if( @goodids == 0 and $ninval >= @cids) {
      foreach my $id (@cids) { delete $aacluster{$id}; }
      $aaclustermain{$mid}= [];
      $nreset++; $reset=0; 
    }
    if($reset) {
      foreach my $id (@goodids) {
        my($oldmain,$pi)=split",",$aacluster{$id};
        $aacluster{$id}= "$mainid,$pi";
        }
      $aaclustermain{$mainid}= \@goodids;
      $nreset++;
    }
    $didmain{$mainid}++;
  }
  warn "# correctAAcluster nreset=$nreset\n" if($debug); # got nreset=212165 for nclust=232463 TOO HIGH?
  return($nreset);
}

sub readAAcdhit {
  my ($bother)= @_;
  my($lt,$lq,$sa,$sm,$ss,$qd,$wq)= (0) x 10;
  my($ok,$fh)= openRead($bother);
  %aacluster=(); # global
  our @cluster=();
  my $iscdhit=0;
  my($cli,$ncl,$nalt,$mainid,$mainlen,$nerr)=(0)x10;

  # FIXME: use aaqual utrbad/poor for main/alt with pi=100; reset main if utrbad and alt utrok
  sub putclus { 
    my($mainid,$cluster1)=@_; ## our @cluster;
    if($mainid) {  
      $aaclustermain{$mainid}=[] unless($aaclustermain{$mainid});
      foreach my $cl (@$cluster1) { 
      my($id,$pi)=split",",$cl; $aacluster{$id}="$mainid,$pi"; 
      push @{$aaclustermain{$mainid}}, $id;
      } 
    }
  }
  
  while(<$fh>) { 
    if(/^>/) { 
      putclus($mainid,\@cluster) if(@cluster); @cluster=(); $mainid=0;
      ($cli)= m/(\d+)/;  $ncl++; #   m/Cluster\s*(\d+)/;  fixed or not?
      
    } elsif(/^(\d+)/) {
      my $i=$1;
      m/(\d+)(nt|aa), >(.+)\.\.\. (.+)$/;  
      my($tlen,$typ,$tid,$pinfo)=($1,$2,$3,$4); 
      my @pmore; ($pinfo,@pmore)= split /\t/, $pinfo;
      my $pi;
      my $ismain=($pinfo =~ /\*/)?1:0;  # FIXME: drop /$/; new merge fnum at end of line now;
      ($tid)= _squishID($tid);
      unless($tid) {
        $nerr++;
      } elsif($ismain) {
        $mainid= $tid; $mainlen=$tlen; $pi=100;
        push @cluster, "$tid,$pi";
      } else {
        # pi for cdhit-est: at 1:492:892:1383/+/100.00%
        # pi for cdhit-aa : at 99.76%   <<<<<
        $pinfo =~ s/at //;  $pinfo =~ s/^\D+//; #? always \digit start
        ($pi)= $pinfo =~ m/([\d\.]+)%/; $pi=~s/\.00//; $pi||=0;
        push @cluster, "$tid,$pi";  $nalt++;
      }
    } elsif(/^\w/) {
      $nerr++; # warn/die not cdhit format ..
    }
  }
  
  putclus($mainid,\@cluster) if(@cluster); ### putspans($mainid) if($mainid); %bspans=();
  warn "# readAAcdhit: nclust=$ncl, nalt=$nalt, nerr=$nerr \n" if $debug;
  return($ncl); 
}

=item readEqualGene

  Table of map equivalences, add as adjunct to align tab of blast/lastz, which have omissionn mistakes ~10%?
	$evigene/scripts/equalgene.pl -in kf2mixx.main.gff -over kf2mixx.main.gff > kf2mixx.main.eqgene

  adjust eqgenes to account for poor gmapping, paralogs map to same locus often, poorly, 
  but blastn/any align says they are different loci
  also have perfect gmap eqgenes, not scored same by blastn, should be called same locus dups
	- add mapqual to eqgene.table ? need it per geneID

  - in readEqualGene, mark which overlap alts are bad, which ok
  - also need to regard mapqual align, Split values to decide
  - need separate classifier to handle various eqgene attributes, decide which tr/alts are bad/good
  
=cut

sub readLocusPairs { # reuse readEqualGene for now
  my($locuspairs)=@_;
}  

sub readEqualGene {
  my($infile)= @_;
  my($nids,$nov,$nxeq, $nne)=(0) x 9; 
  my (%eqgenes,%eqexons,%neqalts,%gmapqual,%negenes); # extended eqgene.tab cols

	my($ok,$inh)= openRead($infile,1);
  while(<$inh>) {  
  	next unless(/^\w/);
		chomp; my @v=split"\t";
  	my($mid,$oid,$overids,$loc,$mapqual,$altnotover)=@v;
		# extended eqgene; mapqual = below; alt_notover column for paralog classing
		
		next if(not $KEEP_G and $mid =~ /_G\d+$/); # dup map id syntax; FIXME OPT to KEEP_G
    ($mid)= _squishID($mid);
		my $issplit=0;
    if($mid =~ s/_C(\d+)$//) { $issplit=$1; } #what of split maps? change mapqual?
    my $pmid= $oids{$mid}||$mid;
    
    # V4: special flags from icchains.eqgene location: 
    #   l4630.t8.icsub == lower qual alt, keep only if other quals good
    #   l4630.t7.icalt.icdup == duplicate alt, keep only if other quals good
    # eqgene mapqual: 85a,100i,2469l,4x,Spl:29%,KB668900 == %align,%ident,bases,exons,Spl:split

    if($issplit and not ($mapqual =~ /Spl:/)) { # FIXME; not accurate fixup
      my($pc)= $mapqual=~m/(\d+)a,/; my $po=90-$pc; $mapqual.=",Spl:$po%";
    }
    
    if($issplit and $gmapqual{$mid}) {
      # dont replace unless this mapqual > other mapqual
      my ($pspl)= ($mapqual=~m/Spl:(\d+)/) ? $1 : 0;
      $issplit=99 if($pspl > 50); #flag lesser
      $gmapqual{$pmid}= $gmapqual{$mid}= "$mapqual\t$loc" if($pspl < 50);
    } else {
		  $gmapqual{$pmid}= $gmapqual{$mid}= "$mapqual\t$loc";  # mapqual eq "na" == missing
		}
		
		## only set validids for mustmapqual ?? otherwise use aligntab validids
    #NOT HERE: if($mustmapqual) { $validids{$mid}= $validids{$pmid}=1; } # V4:test
		
    if($altnotover and $altnotover ne "na") {
      my @nealt= split",",$altnotover; $nne++;
      map{ my $pa= $oids{$_}||$_; 
        $neqalts{$pmid}{$pa}= $neqalts{$mid}{$_}=1; 
        } @nealt;
    }
    
      # add negative val: not-overlocus, for untangling altpar 
      # Arath5EVm000004t1	Arath5EVm000004t2/93xe100,Arath5EVm000005t1/0.0no <?? 0-notover == ignore any align
		#? next if($overids eq "na" or not $overids);
		my @ovd= grep{ $_ ne "na" } split",",$overids;
    $nids++; 
    my $jov=0; my $putspan=0;
		foreach my $ov (@ovd) {
			my($ovd,$cx)=split"/", $ov; 
			next if(not $KEEP_G and $ovd =~ /_G\d+$/); # dup map id syntax
      ($ovd)= _squishID($ovd);
		  my $ovsplit=0;
      if($ovd =~ s/_C(\d+)$//) { $ovsplit=$1; } #what of split maps? change mapqual?
			$cx=~s/^[IC]//; #dont.care# $cx="$cx.$cx" unless($cx=~/\.\d/);

      ## FIXME: xeq only from prot/overeqcdsloc.pl, not from equalgene.pl genes.gff
      # $evigene/scripts/prot/overeqcdsloc.pl $eqopts -in $cdsbltab -out $cdseqgene;
      # Arath5EVm000004t1	noid	Arath5EVm000004t2/93xe100,Arath5EVm000004t3/60xe50	chr4:9613617-9635988:-	93a,100i,12609l,54x	0	lsd:Arath5EVm000004t2

			my($xeq)= ($cx=~s/xe(\d+)//)?$1:0; 
			my($ca,$xa)=split /[\.]/,$cx;  $ca||=0;  $xa||=0;
  	  my $povd= $oids{$ovd}||$ovd;
      
      # if($cx eq '0.0no') {  # USE altnotover instead; or $cs =~ m/not/ ; new, block locus join: EVm004t1  EVm0005t1/0.0no
      #   $negenes{$pmid}{$povd}= $negenes{$mid}{$ovd}= $cx; # new negenes hashtab? or use neqalts ?
      #   $neqalts{$pmid}{$povd}= $neqalts{$mid}{$ovd}=1; 
      # }  
			
			  ## need to deal with issplit/ovsplit parts ..
			my $okeq= ($ca >= $MIN_EQGENE)?1:0;
			if($okeq and $issplit) { 
			  my $cas= $eqgenes{$mid}{$ovd}||0;  
			  if($cas > $ca) { $okeq=0; $putspan++; } # dont add 2nd eq2aln tho 
			  }
			if($okeq) { 
			  $jov++;  
			  $eqgenes{$pmid}{$povd}= $eqgenes{$mid}{$ovd}= $ca; # pct of mid aligned to ovd
			  if($xeq) { $nxeq++; $eqexons{$pmid}{$povd}= $eqexons{$mid}{$ovd}= $xeq; } # new, test
			  
        if($EQ2ALNTAB) { # new, make alntab from eqgene
          my $aq= $aasize{$mid}||0; $aq *=3;  # 3*convert to cds-size
          # my $at= $aasize{$ovd}||0; $at *=3;
          my $aln= int(0.5 + $aq * $ca/100);
          my $aident= $aln; # fake
          if(  my($mapi)= $mapqual=~m/(\d+)i,/) { $aident= int(0.5 + $aq * ($ca/100) * ($mapi/100)); }
          my $bits= 2*$aln; # fake
          my $or= 1;
          my $bspan= [1,$aq,1,$aq,$bits,$aln,$aident,$or];
          %bspans=(); #global for putspans();
          $bspans{$ovd}= [$bspan];
          putspans($mid); $putspan++; %bspans=();
        }
		  }
		}
		$nov++ if($jov);
		if($EQ2ALNTAB and not $putspan) {
      %bspans=(); #global for putspans();
      putspans($mid); $putspan++; %bspans=();
		}
	}
	warn "# readEqualGene: equal=$nov/$nids, eqnot=$nne\n" if $debug;
	return (\%eqgenes,$nids,$nov,\%gmapqual,\%neqalts,$nxeq,\%eqexons); #,\%negenes
# ($eqgenes,$nineqgene,$neqgene,$gmapqual,$neqalts,$neqexon,$eqexons) # drop negenes for neqalts
}

