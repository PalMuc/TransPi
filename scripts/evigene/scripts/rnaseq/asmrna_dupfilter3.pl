#!/usr/bin/env perl
# asmrna_dupfilter.pl

=item about
  
  Filter out hi identity duplicate transcripts from transcript assembly,
  using megablast AFTER identifying best proteins.
  
  Principle is that over-assembling transcript reads with many assembly options
  produces a subset of accurate assemblies in a superset of crappy ones.
  
  Best-protein detection is done on superset (longest orf + cds/utr quality),
  then clustered (cd-hit) to reduce to "best" subset by protein quality.
  
  However cd-hit protein filtering retains nearly identical transcripts assembled
  from same reads, but differ in protein enough to avoid cd-hit cluster.
  
  This step identifies high identity transcripts from bestcd subset, from megablast align,
  then marks/removes the hi-id fraction, keeping longest-protein/longest-tr of each
  hi-id cluster.

=item usage

 $evigene/scripts/rnaseq/asmrna_dupfilter.pl -aa=$pt.aa.count -tr=$pt.tr.count -blast=$pt.mblast\
    > $pt.dupclust


=item inputs

  inputs: 
    protein sizes, transcript sizes, as faCount table (id, size); see aacount for gappy proteins

   tr.megablast output from:
   makeblastdb -in $pt.tr -dbtype nucl -out $pt
   blastn -db $pt -query $pt.tr -task megablast -evalue 1e-19 -perc_identity 99 -dust no -outfmt 7 -out $pt.mblast

  outputs: 
    table of same/subset transcripts, rough equiv of cd-hit-est
    bestset.tr, filtered: 1.bestsetaa.ids only, 2.remove trsame subset (like cd-hit-est, but diff methods)

=item FIXMEs

  2013.aug : IS_CDSALIGN ORIENT or is IMPORTANT : need to know when alt-cds are reversed
  	patch in sumblastpart: $or add to bspans
  	
	## FIXME asmrna_dupfilter.pl: check for ID mismatches == too many zeros in aa.count /tr.count / blastn
	## eg. tables have diff prefixes like 'litova:vel' vs 'litovavel'

	** -eqgene genemap.eqgene use is problem .. maybe fixed, need
	   converts all main class to altmap, loses main link, and
		 causes further samemap-locus loci to be created via NOMAIN patch.
	
	2013.sep:
#m2t: ERR: trclass inline:Funhe2Exy3m004019t1_G2        drop    altmap  Fungr1EG3m004332t1      99/48           0,0,pflag:3
#m2t: ERR: trclass inline:Funhe2Exy3m149508t1_G2        drop    altmap  Funhe2Eq7m074376t1      99/35           0,0,pflag:3
## these are from new asmrna_dupfilter using xxx.eqgene table, should not be in trclass ..
## kfish2evg367mixx.trclass6:296 kfish2evg367mixx.trclass8:821 kfish2evg367mixx.trclass9:821
	grep  Funhe2Exy3m149508t1_G2 kfish2evg367mixx.eqgene
	Funhe2Eq7m074376t1      Funhe2Eq7m074376t1      Funhe2Exy3m149508t1_G2/35.55    1095sc:48614-48822:.
	Funhe2Exy3m149508t1_G2  noid    Funhe2Eq7m074376t1/35.55,Funhe2E6bm091492t1/35.39       1095sc:48708-48907:.
	Funhe2E6bm091492t1      Funhe2E6bm091492t1      Funhe2Exy3m149508t1_G2/35.39    1095sc:48614-48785:.

=item FIXME 201402 drop.mains problem

  201402: 
  AAPART cutoff is problem, with UTRBAD/POOR, causing drop.mains of uniq orthologs 
  AAPART & AAMIN, AAMINBAD, AAMINPOO  need fix to prevent drop uniq/valid main/noclass genes, 
  including recover perfectfrag drops of complete/utrgood trasm for AAPART-bad drop.mains
  part of problem that shorter-AA are sometimes best trasm
  
  update for drop.mains, perfectfrag replacements, needs input list from tr2aacds 
    of perfectfrag/dups from fastanr/cdhit prior steps
      
=item revise classifier/output

  FIXME: make OUTSPANTAB=1 default, see below alnclass.sh as better/faster classifier.
  
  update: fold in tested methods for tr class assignment, to replace outclusters
  -- newer class assignment: using %ident (>98) and %align (>50%)
          from  aabugs4qual/tsaevg/cdsidentclass.sh

  shorter alignments of identity are tested as valid criteria for alternate-transcripts
  using genome-mapped transcripts, hi-ident + part-align gives *mostly* alt-tr, some paralogs,
  depending on species & freq of hi-ident paralogs.
              
  update: add refblast  input table for added classification, per
            aabugs4/aabugs4qual/tsaevg/classparts.pl
        tall4 format: Refid Trid  Bits Iden Algn Rlen Tlen, where Ref/Tr may be swapped columns

  update: add 2ndary alt-tr table, in aa.qual format?,
      from alt-tr called by trasm soft, but excluded by cdhit/other filtering as not distinct
      these all should have same geneid + tNNN suffix as primary input tr (defined in -aasize -blast trself.blast)
      ID matching only used, plus aa.qual, to decide whether to keep/drop
      - input to classifier may be 2nd OUTSPANTAB, generated from 1st pass of this and 2nd-alts.aa.qual table.

=item add UTR BAD/POOR filters

  problem now is that main class accumulates utrbad/poor genes, many w/ other utrbad alternates,
  so as more trasm sets are added to this filtering, main class (and some of others) increase w/ junk.
  
  use input aasizes/aaqual scores
  when reading cdhit/blast aligns (esp cdhit clusters), where feasible swap/replace UTRBAD top/first gene 
    if equivalent utrok gene is in cluster/align (ie. utrok/utrbad have same prot/cds size, or ok is nearly same),
    esp. if utrok is also aacomplete.    
  per evigene/scripts/prot/aabest3.sh

=item add ncRNA classing 2016.02
  
  -- asmrna_dupfilter can do 1st pass, likely want tr2genome mapping for intron measures
  -- putative ncRNA can be pulled from tr2aacds dropset:
        drop.(main|noclass), maybe drop.altmid and frag/part classes
  -- use only subset with no CDS overlap (ie not althi, maybe not altmid)
    .. but will need added full tr/cdna overlap test      
  -- maybe select from only utrbad/utrpoor subset
  -- firmest ncrna class, w/o special tests, would include introns, ie tr2genome exon parts alignments

=item add classifyFullTr 2016.05 


=item add cd-hit-est input alignments

  maybe replace cd-hit-aa with cd-hit-cds(est) for both aa-based and nt-based equivalence/reduction
  of trasm sets.  cd-hit-aa has problem of lumping paralogs w/ silent subsititutions, want to keep
  those in early filtering.  cd-hit-cds/est gives approx same classes of same-locus vs diff-locus as
  blastn or blat (depending on parameters).
    cd-hit-est -c 0.99 -G 0 -aS $pALIGN -l 150 -d 0 -i cacao11pub3ig.cds
        pIDENT=0.99 is good as w/ others; pALIGN=0.25 .. 0.50
        
  cacao11pub3ig.class3
  # alternate classifier, replace both cd-hit aa and mblast-tr ? 
  cacao11g: 44403 cds, 29451 loci, 7630 loci have alts,
  blastn (99% id, >=50%? align): alt class=16179 alts, 484 diff locus; diffclass: noclass=noalts; altmid=236 alts, 1156 diff
      false-alts: 484; false-loci: 236 of 44403
  cdhit25: 29893 clusters, 730 loci have alts in diff clusters; 295 diff loci in same cluster;
      false-alts: 295; false-loci: 730 of 44403
        -- about same as above self.mblast (same loci?)
  cdhit50: 30057 clusters, 841 loci alts diff cluster; 252 diff loci in same cluster
    -- cd-hit pALIGN=0.33 may be good choice;
     
=item updated for blat psl input

  blat -out=psl -t=dna -minIdentity=98 -maxIntron=0 kfish1cds.tr kfish1cds.tr kfish1cds.selfblat
  blat -out=psl -t=dna  -minScore=99  -minIdentity=99 -maxIntron=0 kfish1cds.tr kfish1cds.tr kfish1cds.selfblat99
  #  -minScore=99 or such to reduce fragment aligns, 30 default, nmatch - nmiss
  : tests show blat and megablast give nearly same results; blat can be very slow vs mblast for large tr set
  
  
=item see also
  
  evigene/scripts/prot/aaqual.sh
  evigene/scripts/makeblastscore.pl
  evigene/scripts/rnaseq/asmrna_equalgenes.pl
  
=item author
  
  don gilbert, gilbertd near indiana edu, 2012
  part of EvidentialGene, evigene/scripts/

=cut

use FindBin;
use lib ("$FindBin::Bin/../","$FindBin::Bin"); # 201405 add; assume evigene/scripts/ layout: main, prot/, rnaseq/, ...

use strict;
use warnings;
use Getopt::Long;
#add# use cdna_evigenesub; # replacing local subs 201405

# subs now in asmrna_dupfilter3_overlocsub.pm
# sub classifyFullTr;
# sub overlocusclass;
##bad?# use asmrna_dupfilter3_overlocsub;

my $debug= 0;
use constant VERSION => '2016.07.09'; # BAD_GAPS param added, 25% default is too high
# '2016.05.27'; # classifyFullTr() / overLocusclass()
# '2016.02.11'; # MINBASEOVER/TINYALN changes
# '2014.05.17'; #'2014.02.21'; # test/fix drop.main.utrbad problems missed orthogenes
# '2013.09.07'; # 01'; # 08.09'; #08.07; 03.25; 24; # '2013.02.21';

my $OUTSPANTAB= 1;  # make default until replace outclusters()
my $OVSLOP=6; # FIXME: need overlapslop ~ < 0.01 of max part span ; or < 0.02..0.05 of min part span?
my $pctOVSLOP=$ENV{pctover} || 0.02; # need opt?
   $pctOVSLOP=$pctOVSLOP/100 if($pctOVSLOP > 0.9);

# classifier options
our($AAMIN,$AAPART,$AAMINBAD, $AAMINPOO, $AADUP_IDENT, $OK_CDSUTR, $BAD_CDSUTR, $noRESET_CDSUTR, $BAD_GAPS, 
    $TINYALN, $IS_CDSALIGN,$ALTFRAG,$PHIALN,$NHIALN,$PHI,$PMID,$PLOW);

# 201402: these 4 AAsize cutoffs become options, need to reset to cure drop.main.utrbad bug
# also reset BAD_CDSUTR,OK_CDSUTR options to ignore input aaqual if user sets limits: noRESET_CDSUTR
# FIXME2: need adjust pCDS/utrpoor calc for aaqual "utrorf" cases that now have full mrna size but should be split
$AAMIN =$ENV{aamin}||30; #was 40;   # for aacomplete, utrok
$AAPART=$ENV{aapart}||100; # was 100; # for aapartial # 201402: THIS aapart cut is problem, drop.mains uniq orthologs here
      # 201402: update for drop.mains, perfectfrag replacements, need input list of perfectfrag/dups from fastanr/cdhit prior steps
      # 201402: AAPART & AAMIN, AAMINBAD, AAMINPOO  need fix to prevent drop uniq/valid main/noclass genes, 
      # including recover perfectfrag drops of complete/utrgood trasm for AAPART-bad drop.mains
      # part of problem that shorter-AA are sometimes best trasm
$AAMINBAD=$ENV{aaminbad}||60; #was 200;  # for utrbad class
$AAMINPOO=$ENV{aaminpoo}||60; #was 100;  # for utrpoor class
## 2014 new opts for asmrna_dupfilter2:
## tcas4/tribol beetle: 2014 ncbi has 10 prot under 40aa, 5 are NP curated/refseq; 
## export aamin=30;  export aapart=120; export aaminbad=90; export aaminpoo=60

my $OK_CDSUTR_default= $ENV{pcdspoor} ||60; ##  # too high? changes w/ cds-size
my $BAD_CDSUTR_default= $ENV{pcdsbad} ||30; ## % CDS/trlen
$OK_CDSUTR= 0;
$BAD_CDSUTR= 0; # % CDS/trlen
$noRESET_CDSUTR= 1; # flag to ignore input tqual/aaqual utrbad

my $MINUTR=$ENV{minutr}||300; # ~ 300b "fixed" average utr sizes, maybe too low, 
# see sample of org gene set pcds for small aa .. per species min-utr, ranging from 400b insect .. 700b fish .. 1000b mouse
# should adjust %cds bad for cds size, for large cds> 10k, pcds >=90%, for small cds < 300, pcds may be 33%
# should replace utrbad/poor flags with coding potential flags, adding some cds-seq cp calcs.

# classifier should also downweight qual by more than aaspan-aagaps as aagaps damage homol scores
# eg:  score = aanogapsize - aagaps ? reduce qual by ngaps
#old# $BAD_GAPS= 25;  # % gaps in AA; 160709: make this opt, 25% too high, 5% maybe best default?
## BUT review of homol sez gappy prots can still be best/only locus rep.. mostly want to reduce qual score of gappy prots.
$BAD_GAPS= $ENV{aagapmax} || $ENV{BAD_GAPS} || 10; # was 25; # need opt? 'aagapmax' ?

## add UTR BAD/POOR filters, per evigene/scripts/prot/aabest3.sh
$AADUP_IDENT=98; # %ident, option for aacluster ident drops

$TINYALN = 20; # %align, was 25; was $ENV{mina}||50; ignore less than this == MINAL
  ## UPD 2016.02: -tinyaln 35  bad (sometimes), need opts and diff test: tiny by bases not pct overlap
  ## ie. cds-exon overlap should be valid for >= tinybases, eg. 90? need to sample ref species alt exon sizes
  ## UPD 2016.02: reinstate MINALIGN, for tiny by bases not pct overlap, or reuse TINYALN for bases overlap 
  ## ?? UseTINYALNBASES only for IS_CDSALIGN ?
use constant UseTINYALNBASES => 1;   
my $TINYALNBASES= $ENV{MINALIGN}||90; # was $MINALIGN= 90; # REUSE//NOT USED NOW; use just TINYALN; change to several levels
my $MINCDS = $ENV{MINCDS} || 90; # or 3*MINAA  
my $TINYALNCDSLEN= 3*$TINYALNBASES;
my $MIN_EQGENE= $ENV{MINEQGENE}||15; # %equalgene cds-align, was 33; # global/opt..

$ALTFRAG= 0.5;
$PHIALN=$ENV{ahi} || 98; # was 99; was 65 !!;  DROP THIS; use mina?
$NHIALN=$ENV{nhi} || 9; # what? for short cds cancel few base diff < PHIALN
#  my $tinyalndiff= ($aln - $alnmax < $NHIALN) ? 1 : 0 # TEST3 add

## tr-self %identity levels to classify alt-tr
## <90% ident prob too low to class alt-tr; use pmid=95 plow=90 ??
$PHI = $ENV{phi} ||99; 
$PMID= $ENV{pmid}||90; 
$PLOW= $ENV{plow}||80;  

$IS_CDSALIGN= $ENV{cdsw}||0; # WHICH default?  tralign tests better than cds, but cds for other, eg. cds-cdhits
my $SORTEDALIGNTAB=0; # debug input: sort -k2,2nr -k7,7nr -k6,6n evg2anofunz4c.alntab
my $OIDFIX=undef;
my $DOOVERLOCUS=0;

my ($aasizes,$trsizes,$blatpsl,$blastab,$lastz,$bcdhit,$aablast,$aanames,$aacdhit,$outeqtab,$outclass,
    $dupids,$logfile,$head,$eqgene)= (0) x 20;

my $optok= GetOptions(
  "aasizes=s", \$aasizes, 
  "trsizes=s", \$trsizes, 
  "ablastab=s", \$aablast,   # this is traa-refaa.blastp
  "anames=s", \$aanames,   # variant aablast used also for naming
 		## FIXME: allow -ablastab to be .names table, for tr2aacds, ie decide by file name
  "acdhit=s", \$aacdhit,    # this is traa-self.cdhit.clstr

  "blastab=s", \$blastab,    # this is tr-self.blastn, -CDSALIGN for cds-self.blastn
  "lastz=s", \$lastz, # lastz general format; was -blastz option -blast[ab] conflict **
  "blat=s", \$blatpsl, # tr-self.blat;  other input format opts here..?  -informat=xxx
  "bcdhit=s", \$bcdhit, #  tr-self.cdhits.clstr; 
  # ntalign=s, ntformat=blastab|blatpsl|cdhit ..
  "dupids=s", \$dupids, #  

  "eqmap|eqgene=s", \$eqgene,    # 130901: mapping equivalence table, of $evigene/equalgene.pl 
  
  "aligntab|outeqtab=s", \$outeqtab,  ## this should have aliases: -outalntab and -inalntab, or just -aligntab ?  
  "outclass=s", \$outclass,   # 2nd out: outclasstab
  "sortedaligntab!", \$SORTEDALIGNTAB, 

  # 201402 option update: $AAMIN,$AAPART,$AAMINBAD, $AAMINPOO
  "AAMIN=i", \$AAMIN,   
  "AAPARTMIN=i", \$AAPART,   
  "AABADMIN=i", \$AAMINBAD,   
  "AAPOOMIN=i", \$AAMINPOO,   
  "aagapmax=i", \$BAD_GAPS,   

  "TINYALN|MINALIGN=i", \$TINYALN,  # pTINYALN ? FIXME .. both pctMINALN, basesMINALN
  "pCDSOK=i", \$OK_CDSUTR,   # CDSOKUTR
  "pCDSBAD=i", \$BAD_CDSUTR, # CDSBADUTR 
  "ALTFRAG|fragment=s", \$ALTFRAG,  # pALTFRAG ?
  "OUTSPANTAB!", \$OUTSPANTAB, # now fixed=1 ? drop opt
  "CDSALIGN!", \$IS_CDSALIGN, 
  "overlocus!", \$DOOVERLOCUS, 
  "debug!", \$debug, 
  #  "logfile=s", \$logfile,  ## add??
  );


warn "# EvidentialGene asmrna_dupfilter.pl VERSION ",VERSION,"\n" if($debug); # change to loggit() ?
my $hasbalign= ($blastab or $lastz or $blatpsl or $bcdhit) ? 1 : 0;

die "usage:  asmrna_dupfilter.pl -aasize=name.aa.count -trsize=name.tr.count 
 input rna-align: -blast=name.blastn | -blat=name.blatpsl | -bcdhit=name.cdhit.clstr | -aligntab=name.aligntab
 opts: -CDSUTR=$OK_CDSUTR percents  -ablastab=traa-refaa-blast.table  .. more options, see source.
" unless($optok and $aasizes and ($hasbalign or $outeqtab)); ## and $trsizes < dont need if aasizes=aa.qual


$noRESET_CDSUTR=($OK_CDSUTR>0 or $BAD_CDSUTR>0)?0:1; # user option set, ignore input aaqual 201402 update
$OK_CDSUTR ||=$OK_CDSUTR_default;
$BAD_CDSUTR||=$BAD_CDSUTR_default;

if($TINYALN>0) {
  $TINYALN= 100*$TINYALN if($TINYALN<1);
  ### $TINYALN=$MINALIGN if($TINYALN>$MINALIGN); # drop MINALIGN for TINYALN
}
$ALTFRAG= $ALTFRAG/100 if($ALTFRAG>1); # prop not percent

# should be our() shared vars ..
my(%aasize, %sizeval, %trsize, %aaqual, %oids, %oidsof, %cdsoff); # all from readSizes(); one hash? == sizeval
my(%aablast, %aablastref,$naabl); $naabl=0;
my(%aacluster,%aaclustermain); # globals for readAAcdhit
my(%better, %outrows, %validids, %bspans); ## %validids was %outids; 
# oids from aasize? if($OIDFIX and %oids) { $td= $oids{$td}||$td; $rd= $oids{$rd}||$rd; }

my(%dupids, %dupfirst, $ndupdrop); $ndupdrop=0;
use constant DUPFILTER1 => 1; # tests ok
use constant DUPFILTER2 => 0; # bad?

sub MAINstub {}

my $OUTH= *STDOUT;

## set TMPDIR to localdir always? dont use problematic perl mods
# unless($ENV{TMPDIR}) {  ## sort TMPDIR, cant use /tmp
#   require Cwd; if(my $cdir=Cwd::getcwd()) { $ENV{TMPDIR}=$cdir; } 
#   # 201506: perl fail cwd2 .. require Cwd::getcwd()
#   #o# if(my $cdir=FindBin::cwd2()) { $ENV{TMPDIR}=$cdir; } 
# }
# unless($ENV{TMPDIR}) { if(my $cdir=FindBin::cwd2()) { $ENV{TMPDIR}=$cdir; } } ## sort TMPDIR, cant use /tmp
## perl-lib-version bug: Undefined subroutine &FindBin::cwd2 .. use `pwd` ? other hack

readSizes();
readDupIds($dupids)   if($dupids);
 		## FIXME: allow -ablastab to be .names table, for tr2aacds, ie decide by file name; was '.names'
if($aablast and $aablast =~ /\.name/ and not $aanames) { $aanames= $aablast; $aablast=""; }
($naabl)= readAAnametab($aanames) if($aanames); # prefer this now?
($naabl)= readAAblast($aablast)   if($aablast and not $naabl);  # one or other of aablast,aanames

readAAcdhit($aacdhit) if($aacdhit); # also correctAAcluster()

use constant EQGENE_OVERRIDES_ALN => 0;
use constant EQGENE_CHANGES_NOALN => 0; # not quite right yet, 160217, tho new eqgene data useful
## $eqgenes->{$tid}{$sid} >= $MIN_EQGENE alignment
## FIXME: need to adjust eqgenes to account for poor gmapping, 
##   paralogs map to same locus, poorly, 
##   but blastn/any align says they are different loci
## BUT also have perfect gmap eqgenes, not scored same by blastn, should be called same locus dups

#old# my($eqgenes,$nineqgene,$neqgene)= ($eqgene) ? readEqualGene($eqgene) : (undef,0,0);
## \%eqgenes,$nids,$nov,\%gmapqual,\%neqalts,$nxeq,\%eqexons
my(%eqflag); # set in identityclass, output?
my($eqgenes,$nineqgene,$neqgene,$gmapqual,$neqalts,$neqexon,$eqexons)
    = ($eqgene) ? readEqualGene($eqgene) : (undef,0,0,0,0,0);


# FIXME 201405, allow reuse $outeqtab even if hasbalign.. for tr2aacds restarts
if( ($ENV{useouteqtab} or $ENV{usealigntab} ) and -s $outeqtab ) { $hasbalign=0; }

my $nbalign=0;
if($hasbalign) {
  if($outeqtab) {
    rename($outeqtab,"$outeqtab.old") if( -f $outeqtab );
    open(OUTH,">$outeqtab") or die "write $outeqtab";
    $OUTH= *OUTH;
  }
  
  # allow input of OUTSPANTAB instead of regenerate?
  if($blastab) { ($nbalign)= readblasttab($blastab); }
  elsif($lastz) { ($nbalign)= readlastz($lastz); }
  elsif($bcdhit) { ($nbalign)= readcdhit($bcdhit); }
  elsif($blatpsl) { ($nbalign)= readblatpsl($blatpsl); } # detect from input table ??
  if($outeqtab) { close($OUTH); $OUTH=undef; } # *STDOUT ?
  warn "# readalign: nids=$nbalign to $outeqtab\n" if $debug;
}

# unless($OUTSPANTAB) 
if($outeqtab) {
  my $infile=$outeqtab;   # infile == outeqtab ? STDIN?
  my $insorted=$SORTEDALIGNTAB; # add input opt for sorted table ?? test classing bugs
  unless($infile and -f $infile) { 
    die "ERR: missing input align table -aligntab $outeqtab";
  } elsif($infile =~ /\.gz/) { 
    die "ERR: cant use gzipped align table -aligntab $outeqtab";
  }

  # >> set this BEFORE correctAAcluster from read{blast|cd|blat} : $validids{$id}
  if($aacdhit) {
    my $havevalid= scalar(%validids)?1:0;
    $havevalid= readIdsFromAlnTab($infile) unless( $havevalid); 
    correctAAcluster($havevalid); # update %aacluster,%aaclustermain
  }

  $OUTH= *STDOUT; 
  if($outclass) {
    rename($outclass,"$outclass.old") if( -f $outclass );
    open(OUTC,">$outclass") or die "write $outclass";
    $OUTH= *OUTC;
  }
  
  if($DOOVERLOCUS) {
  $IS_CDSALIGN=0; # always?
  overlocusclass($OUTH,$infile,$insorted, $IS_CDSALIGN);
  } else {  
  identityclass($OUTH,$infile,$insorted); 
  }
  if($outclass) { close($OUTH); $OUTH=undef; } # *STDOUT ?
}

# if(0) { outclusters($OUTH) unless($OUTSPANTAB); } # DROP this old version

#.................................

sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }
sub bint { my $b=shift; return ($b<0)?0 : ($b=~/e\+/)? int($b) : $b; }

sub openRead { # in cdna_evigenesub.pm
  my($fna, $nostdin)= @_; my($ok,$hin)= (0,undef);
  $ok= ($fna =~ /\.gz$/) ? open($hin,"gunzip -c $fna|") 
  	 : ($fna =~ /stdin|^-/ and not $nostdin) ? *STDIN 
  	 : open($hin,$fna);  
  # loggit(1,"ERR: openRead $fna") unless($ok);
	die "ERROR: openRead $fna" unless($ok);
  return ($ok,$hin);
}



sub readSizes {  # see also cdna_evigenesub.pm getAaQual()

  my($naa,$ntr,$nerr,$ok,$inh)=(0) x 10;
  if($aasizes) {  
    ## fix for aacount gaps: id,size,gaps : NOT NOW, aa.qual: id,size-gaps,gaps,..
    ## drop faCount? require aa.qual here?

    # ($ok,$inh)= openRead($aasizes,1);
    open(F,$aasizes) or die "FAIL: read $aasizes ..."; 
    # if($aasizes =~ /count|qual/) { 
    # } else { open(F,"faCount $aasizes |") or die "FAIL: faCount $aasizes ..."; } # drop this..

## FIXME 201405: Maybe replace w/ cdna_evigenesub:getAaQual() -- somewhat different %AAQUALH data
## getAaQual id => $alen,$pctcds,$acv,$aqual1
## FIXMEd for cds.qual not aa.qual, col2 == cds-size, col4=aasize,qual, col5=trsize
## FIXME 2016.05, cds.qual now *may* have "Code/Noncode,.." column after tlen, before offs, 

    my $iscds=0;
    while(<F>) { 
      next if(/^#/ or /^total/); 
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
      ## all in one?
      @{$sizeval{$id}}{ qw(aasize cdsize gap aaqual trsize cdsoff oids codepot) }
        = ($alen,$cdlen,$gap,$aqual,$tlen,$offs,$oids,$codepot);

      
      } close(F); 
  }
  
  if($trsizes) {
  if($trsizes =~ /^aaqual/ or $ntr>0) { 
    # got above;now default
  } elsif($trsizes =~ /^aasize|^cdssize/) { # for tr == cds
    foreach my $id (keys %aasize) { $trsize{$id}= 3*$aasize{$id}; }
  } else {  
    ## drop faCount? expected use is aa.qual w/ trsize column
    # ($ok,$inh)= openRead($trsizes,1);
    open(F,$trsizes) or die "FAIL: read $trsizes ..."; 
    # if($trsizes =~ /count|qual/) { 
    # } else { open(F,"faCount $trsizes |") or die "FAIL: faCount  $trsizes ..."; }
    while(<F>) { next if(/^#/ or /^total/); my($id,$al)=split; $trsize{$id}=$al; $ntr++; } close(F); 
  }
  }
 
 warn "# readSizes: naa=$naa; ntr=$ntr\n" if $debug;
 return($naa,$ntr); 
}

# fix from fastanrdb all.cds > all_nr.cds >> hdr has cds-identical ids; should have prefiltered this
sub readDupIds {
  my($infile)= @_;
  my($nids,$ndups,$inh)=(0,undef);
  %dupids= %dupfirst=();
	# ($ok,$inh)= openRead($infile,1);
  open(IN,$infile) or die "reading $infile"; $inh=*IN; 
  while(<$inh>) {  
    s/^>//; # if from fastanrdb
    my  @dupids= split; # grep /\w/ or any other?
    next unless(@dupids>1); #? or record all ids?
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

=item readEqualGene

Table of map equivalences, add as adjunct to align tab of blast/lastz, which have omissionn mistakes ~10%?
	$evigene/scripts/equalgene.pl -in kf2mixx.main.gff -over kf2mixx.main.gff > kf2mixx.main.eqgene
* need to change main/over ids to input oids somewhere ..
Mainid                  Oid                     Overlapids
Funhe2Exy3m110884t1     Fungr1EG3m041003t1      Funhe2Exy3m129932t1/C97.72,Funhe2Exy3m110083t1/82.77    10098sc:638046-638339:.
Funhe2Exy3m125954t1     Fungr1EG3m043564t1      na      479sc:197743-197931:.
Funhe2Exy3m091394t1     Fungr1EG3m037315t1      Funhe2Exy3m076422t1/I100,Funhe2Exy3m083869t1/89.90,Funhe2Exy3m044035t1_C1/87.88,Funhe2Exy3m063866t1/87.87,Funhe2Exy3m054191t1_C1/87.83,Funhe2Exy3m059312t1/87.81,Funhe2Exy3m025662t1/87.79,Funhe2Exy3m060328t1_C1/86.87,Funhe2Exy3m050942t1/86.86,Funhe2Exy3m052474t1/86.86,Funhe2Exy3m085750t1/86.86,Funhe2Exy3m050513t1/86.84,Funhe2Exy3m040128t1/86.82,Funhe2Exy3m052598t1/86.78,Funhe2Exy3m052599t1/86.78,Funhe2Exy3m049962t1/85.85,Funhe2Exy3m054654t1/85.85,Funhe2Exy3m123553t1/57.83,Funhe2Exy3m026282t1_C1/53.87,Funhe2Exy3m090395t1_C2/0.97    1967sc:10172-10679:-

FIXME: need to adjust eqgenes to account for poor gmapping, 
  paralogs map to same locus, poorly, 
  but blastn/any align says they are different loci
BUT also have perfect gmap eqgenes, not scored same by blastn, should be called same locus dups
	-- add mapqual to eqgene.table ? need it per geneID
	-- or separate table?

Fixme: dont output these to trclass
#m2t: ERR: trclass inline:Funhe2Exy3m149508t1_G2        drop    altmap  Funhe2Eq7m074376t1      99/35           0,0,pflag:3

grep  Funhe2Exy3m149508t1_G2 kfish2evg367mixx.eqgene
Funhe2Eq7m074376t1      Funhe2Eq7m074376t1      Funhe2Exy3m149508t1_G2/35.55    1095sc:48614-48822:.
Funhe2Exy3m149508t1_G2  noid    Funhe2Eq7m074376t1/35.55,Funhe2E6bm091492t1/35.39       1095sc:48708-48907:.
Funhe2E6bm091492t1      Funhe2E6bm091492t1      Funhe2Exy3m149508t1_G2/35.39    1095sc:48614-48785:.
	
=cut

=item eqgene classifier

  - do this in readEqualGene, mark which overlap alts are bad, which ok
  - also need to regard mapqual align, Split values to decide
  
  need separate classifier to handle various eqgene attributes, decide which tr/alts are bad/good
  eg eqgene classifier for this case: good g131, bad g453t5
  .. for this case need to know that g453t5 <mismap> g453t1,2,3,4 .. count overlap/alt/locus? drop outliers? use Split info?
ok Anofunz4gEVm000131t1	noid	Anofunz4gEVm000452t5/19	KB668936:299087-306372:-	99a,99i,7287l,3x	0
ok Anofunz4gEVm000131t2	noid	Anofunz4gEVm000452t5/33	KB668936:299737-303981:-	98a,99i,4269l,2x	0
ok Anofunz4gEVm000131t3	noid	Anofunz4gEVm000452t5/58	KB668936:299632-301584:-	70a,99i,1809l,2x	0
ok Anofunz4gEVm000452t1	noid	na	KB668900:3696-9724:-	99a,100i,4410l,8x	0
ok Anofunz4gEVm000452t2	noid	na	KB668900:3696-9724:-	99a,100i,4410l,8x	0
ok Anofunz4gEVm000452t3	noid	na	KB668900:3696-9724:-	96a,100i,4329l,8x	0
ok Anofunz4gEVm000452t4	noid	na	KB668900:3696-8875:-	100a,100i,3372l,9x	0
bad Anofunz4gEVm000452t5	noid	Anofunz4gEVm000131t1/56,Anofunz4gEVm000131t2/56,Anofunz4gEVm000131t3/43	KB668936:300530-301921:-	85a,100i,2469l,4x,Spl:29%,KB668900	0
ok Anofunz4gEVm000452t6	noid	na	KB668900:6592-8636:-	79a,100i,948l,3x	0

  g131,3/3 alts are over 1/6 g453 alts
  -- should this indicate g453t5 is not-much-over other g453t alts?
  
=cut

sub readEqualGene {
  my($infile)= @_;
  my($nids,$nov,$nxeq, $nne)=(0) x 9; 
  my (%eqgenes,%eqexons,%neqalts,%gmapqual); # extended eqgene.tab cols
		  #* change this eqgene/altmap classing,
		  # NO cds-align b/n td,qd means something, poor mapping paralogs end up here
		  # .. need (a) gmap qual score (align,ident,split), and
		  # ..      (b) paralog flag reverse of eqgene, for classed-alts that *dont* gmap same locus
	
		  # TEST1602: ($eqgenes,$nineqgene,$neqgene,$gmapqual,$neqalts)= readEqualGene($eqgene) extended table

  # open(IN,$infile) or die "reading $infile"; $inh=*IN; 
	my($ok,$inh)= openRead($infile,1);
  while(<$inh>) {  
  	next unless(/^\w/);
		my @v=split"\t";
  	my($mid,$oid,$overids,$loc,$mapqual,$altnotover)=@v;
		# extended overeqgene from altclassed-cds-genome.blastn calc: add alt_notover column for paralog classing
		
		next if($mid =~ /_G\d+$/); # dup map id syntax
    $mid =~ s/_C\d+$//; #what of split maps? change mapqual?
    my $pmid= $oids{$mid}||$mid;
    
    # new bleqgene mapqual: 85a,100i,2469l,4x,Spl:29%,KB668900 == %align,%ident,bases,exons,Spl:split
		$gmapqual{$pmid}=$gmapqual{$mid}="$mapqual\t$loc"; # or $loc\t$mapqual
		
    if($altnotover and $altnotover ne "na") {
      my @nealt= split",",$altnotover; $nne++;
      map{ my $pa= $oids{$_}||$_; 
        $neqalts{$pmid}{$pa}= $neqalts{$mid}{$_}=1; 
        } @nealt;
    }
    
		next if($overids eq "na" or not $overids);
    $nids++; my $jov=0; 
		my @ovd= grep{$_ ne "na"} split",",$overids;
		foreach my $ov (@ovd) {
			my($ovd,$cx)=split"/", $ov; 
			next if($ovd =~ /_G\d+$/); # dup map id syntax
		  $ovd =~ s/_C\d+$//;  ## what of split maps: id _C[12] ? skip or chop _C?
			$cx=~s/^[IC]//; #dont.care# $cx="$cx.$cx" unless($cx=~/\.\d/);

			# overeqcdsloc.pl test syntax,may change, adds exon-equal:  
			#  Anofunz4jEVm000002t47		Anofunz4jEVm000002t1/100xe100,..
			# problem cases, mis-align due to missing exon align of true alt, should aln<100 reduce xe value?
			# Anofunz4hEVm000070t2	Anofunz4hEVm000070t1/87xe100  	KB668689:297566-312085:-	87a<<,99i,10536l,8x
			# Anofunz4hEVm000070t1	Anofunz4hEVm000070t2/89xe88   	KB668689:297566-312085:-	100a,99i,10314l,8x
			
			my($xeq)= ($cx=~s/xe(\d+)//)?$1:0; 
			
			#old# $cx.=".100" if($cx =~ /^I100/); $cx=~s/^[IC]//; 
			my($ca,$xa)=split /[\.]/,$cx;  $ca||=0;  $xa||=0;

## FIXME : input has oids, want pubids (also) from readSize..
#   # my $qoids= $oidsof{$qd};  my $toids= $oidsof{$qd};
# $td= $oids{$td}||$td; $qd= $oids{$qd}||$qd; <<<

			if($ca >= $MIN_EQGENE) { 
  		  my $povd= $oids{$ovd}||$ovd;
			  ## ca is not reciprocal equal now: mid-larger over ovd-smaller, ca is large portion of ovd
			  ## mid-smaller over ovd-larger, ca is small portion of ovd
			  $jov++;  $eqgenes{$pmid}{$povd}= $eqgenes{$mid}{$ovd}= $ca; # pct of mid aligned to ovd
			  if($xeq) { $nxeq++; $eqexons{$pmid}{$povd}= $eqexons{$mid}{$ovd}= $xeq; } # new, test
			  }
		}
		$nov++ if($jov);
	}
	warn "# readEqualGene: nequal=$nov/$nids\n" if $debug;
	return (\%eqgenes,$nids,$nov,\%gmapqual,\%neqalts,$nxeq,\%eqexons);
	#oold#return (\%eqgenes,$nids,$nov);
}

			## Argument "58^I1802sc:10833-11091:+\n" isn't numeric : split $ov not $_
			#x# if($xa > $ca) { $ca=$xa if($ca<$MIN_EQGENE and $xa>25+$MIN_EQGENE); } # what? adjust ca to xa? may be wrong. but CDS mapping can be wrong.
			## ^^ drop this, want only cds overlaps, otherwise utr effects such as 2 gene parts cause problems
			

=item readAAnametab / readAAblast

  maybe revise this to also use trasm.names table, computed for publicset naming,
  has essentially same info w/ best ref per tr.
  names format:
  TrID   Name   Align_score  RefID  RepID (uniprot)
  Funhe2Exx4m000455t12    CDD: Na_trans_assoc, Sodium ion transport-associated    100%,245/230,1997       CDD:219069      pfam06512
  Funhe2Exx4m000455t12    Sodium channel protein type 5 subunit alpha     100%,2057/2016,1997     RefID:UniRef50_Q14524   RepID:SCN5A_HUMAN
  Align_score = val%,nalign/nref,ntr
  
  make aablast.tab from query-ref.blastp
   a. (best)?
   cat evg2anofunz4g-bugsref.aa.blastp | \
     env blsum=1 aa="$query.aa.qual,$ref.aa.qual" ISCORE=4 nst=3 $evigene/scripts/blast2bestgenes.pl \
     > evg2anofunz4g-bugs6ref.aa.bltab
     # bltab: Query   Source     Bits    Iden    Algn    Qlen    Slen
   b. older, mostly same for cols 1-7
    $evigene/scripts/makeblastscore3.pl -tall -aasize "$query.aa.qual,$ref.aa.qual" evg2anofunz4g-bugsref.aa.blastp \
     > evg2anofunz4g-bugsref.aa.btall
   
     
=cut 

sub readAAnametab
{
  my($aanametab)= @_;
  my $swapids=0; my $naabl=0; my $nblerr=0;
  my %bscore;
  open(F, $aanametab) or die "FAIL: read $aanametab";
  while(<F>) { next unless(/^\w/); 
    my($td,$name,$alnscore,$rd,@more)=split"\t"; 
    # alnscore format expected: 72%,3270/4555,3282 ;  may be '72' or '72%' only
    #x my($apct,$aln,$refsize,$tsize)= $alnscore =~ m=^(\d+)%,(\d+)/(\d+),(\d+)=;
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


sub readAAblast 
{  
  my($aablast)= @_;
  my $swapids=0; my $naabl=0; my $nblerr=0;
  my %bscore;
  if($aablast =~ /\.names$/) { return readAAnametab($aablast); }
  
  ## precheck format before sort of possibly very large file ..
  use constant nCHECK => 29;
	my($ok,$inh)= openRead($aablast,1);
  # open(F,$aablast) or  die "FAIL: read $aablast ..."; 
  while(<$inh>) { next unless(/^\w/); 
    ## FIXME: allow orig blastp table format, not postprocess .tall4 version ?
    ## at least check for blastp table
    my($td,$rd,@v)=split"\t"; $naabl++; # expect: trid refid bits ident align .. may be refid,trid instead

    if(%oids) { 
      $OIDFIX=1 if($oids{$td} or $oids{$rd});   # is aasize{oid} ok here
      if($OIDFIX){ $td= $oids{$td}||$td; $rd= $oids{$rd}||$rd; }
    }
    if($swapids==1) {
      ($td,$rd)=($rd,$td);
    } elsif($swapids==0) {
      if(@v >= 12 and $v[8] =~ /^\d/) { $nblerr++; } # blast.tab
      else {
      if($aasize{$td}) { $swapids= -1; }
      elsif($aasize{$rd}) { $swapids= 1; }
      else { $nblerr++; }
      }
    }
    
    last if($naabl > nCHECK);
  } close($inh);
  if($nblerr>2 or $swapids==0) { 
    die "ERR: expect table of aablast scores: trid refid bitscore identity align ..\n"
      ." $nblerr trids from aasize not found in -aablast=$aablast\n";
  }
  
  $naabl=$nblerr= 0;  
  #?? add -aablastsorted option? sort wastes time if not needed.
  #old# open(F,"sort -k3,3nr $aablast |") or die "FAIL: read $aablast ..."; 
	# ($ok,$inh)= openRead($aablast,1);
	
	## FIXME sort for swapids, td first, refd 2nd when tied scores, prob not much effect
	## FIXME TEST1602: -k3 bitscore sort may be flaky, opt try -k4 ident or -k5 align
	my $sortord=($swapids==1)?'-k3,3nr -k5,5nr -k2,2 -k1,1' : '-k3,3nr -k5,5nr -k1,1 -k2,2';
  # open(F,"sort $sortord $aablast |") or die "FAIL: read $aablast ..."; 
  $ok= ($aablast =~ /\.gz$/) ? open($inh,"gunzip -c $aablast | sort $sortord |") 
        : open($inh,"sort $sortord $aablast |");  
  
  while(<$inh>) { next unless(/^\w/); 
    ## FIXME: allow orig blastp table format, not postprocess .tall4 version ?
    ## at least check for blastp table
    my($td,$rd,@v)=split"\t"; $naabl++; # expect: trid refid bits ident align .. may be refid,trid instead
    
    if($OIDFIX) { # and %oids %oid2pubid; aablast, others?
      $td= $oids{$td}||$td; $rd= $oids{$rd}||$rd;
    }
    
    if($swapids==1) {
      ($td,$rd)=($rd,$td);
    }
#    # done already ..
#     } elsif($swapids==0) {
#       if(@v >= 12 and $v[8] =~ /^\d/) { $nblerr++; } # blast.tab
#       else {
#       if($aasize{$td}) { $swapids= -1; }
#       elsif($aasize{$rd}) { $swapids= 1; }
#       else { $nblerr++; }
#       }
#       if($nblerr>0 and $naabl > 19) { 
#         die "ERR: expect table of aablast scores: trid refid bitscore identity align ..\n"
#           ." trids from aasize not found in -aablast=$aablast\n";
#       }
#     }
    
    my $bscore= $v[0]; # bits, ident, algn; want choice of? ident or algn maybe better
    # local $^W = 0;   # no warnings; no help
    ## maybe fixme: missing some uniq blastref{rd} here? pull aablastref out of aablast{} ?
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


use constant SPANSUM => 1;
use constant ADDSPANS1605 => 1;

sub putspans {
  my($lq)= @_;
  # fixme: save no-match lq ids:  lq, aq, wq, na.... or self-score ?
  my $nmatch=0;
  # fixme: sort output bspans by tidn/taln? $bspans{$b}->[0]->[4] = xbit
  # fixme? throw away dupids here? should be in same bspans align cluster
  
  # our(%dupids, %dupfirst);
  my $dupfirst=""; 
if(DUPFILTER2) {    
  if($dupids) { $dupfirst= $dupids{$lq} || ""; }
}
  
  foreach my $lt (sort keys %bspans) {
    next if($lt eq $lq); # is lt eq lq allowed here?
    
    ## move this dupid filter before into readblastab, readcdhit ? see active DUPFILTER1 
if(DUPFILTER2) {    
    if( $dupids and $dupids{$lt} ) {
      if($dupfirst eq $dupids{$lt}) { $ndupdrop++; next; }
      $dupfirst= $dupids{$lt};
    }
}
    
    my @bspans= @{$bspans{$lt}};
    my ($tbit,$taln,$tidn,$aq,$at,$wq,$wt,$torient,$xbm,$xem,$tbm,$tem)= (0) x 19;
    foreach my $sp (@bspans) {
      my($xb,$xe,$tb,$te,$xbit,$aln,$aident,$or)= @$sp; # 2013.aug: IS_CDSALIGN add $or
      $tbit += $xbit; $taln+= $aln; $tidn+= $aident; 
      ##$torient+=$or; ## weight by aln so tiny -or dont throw it off?
      $torient += $aln * $or; ## weight by aln so tiny -or dont throw it off?
if(ADDSPANS1605) {
      $xbm=$xb if($xbm==0 or $xb<$xbm); $xem=$xe if($xe>$xem);
      $tbm=$tb if($tbm==0 or $tb<$tbm); $tem=$te if($te>$tem);
}
      }
    # my $mis= $taln - $tidn; # dont need this; replace w/ tidn
    if($OUTSPANTAB) { 
      # add Qsize,Tsize  cds/trlen ?
      $aq= $aasize{$lq}||0; $aq *=3;  # 3*convert to cds-size
      $at= $aasize{$lt}||0; $at *=3;
      $wq= $trsize{$lq}||0; 
      $wt= $trsize{$lt}||0;  
      ## NO: add here?? $aaqual{$lq} ; $aaqual{$lt}; 
      ## FIXME taln may be rel CDS-len (aq,at) or TR-len (wq,wt)
      my $alnmax= _max(1, ($IS_CDSALIGN) ? _min($aq,$at) : _min($wq,$wt) );
      
      ## UPD 2016.02: -tinyaln 35  bad (sometimes), need opts and diff test: tiny by bases not pct overlap
      ## eg. TINYALNBASES= 90; MIN_PCTOVERLAP=10;
  if(UseTINYALNBASES) {
      if($alnmax >= $MINCDS) { # never skip tiny prots, MINCDS =~ 90 bases ?
        my $minbaseover= ($alnmax >= $TINYALNCDSLEN)? $TINYALNBASES : $alnmax * 0.10; ##$MIN_PCTOVERLAP;
        next if($taln < $minbaseover);
      } 
  } else { # old, pct-TINYALN problem    
      next if( (100 * $taln / $alnmax ) < $TINYALN); # skip trival matches
  }      
      $nmatch++;
      
      # 2013.aug: IS_CDSALIGN add $torient here? only care if $tor < 0; add sign to one of these cols?
      ## cant use taln, will screw up sort -n by maxalign; use tbit, unused now in scoring
      if($IS_CDSALIGN and $torient<0) { $tbit= -$tbit; } #? what followon problems does this cause?
if(ADDSPANS1605) {
      print $OUTH join("\t",qw(Qid Qclen Qtlen Tid Tclen Ttlen Align Ident Bits Qspan Tspan))."\n" unless($head++);
      print $OUTH join("\t",$lq,$aq,$wq,$lt,$at,$wt,$taln,$tidn,$tbit, "$xbm-$xem","$tbm-$tem")."\n"; 
} else {      
      print $OUTH join("\t",qw(Qid Qclen Qtlen Tid Tclen Ttlen Align Ident Bits))."\n" unless($head++);
      print $OUTH join("\t",$lq,$aq,$wq,$lt,$at,$wt,$taln,$tidn,$tbit)."\n"; 
}
      $validids{$lq}++; $validids{$lt}++;
      } 
    else { 
      puts($lq,$lt,$taln,$tidn);   $nmatch++; # DROP this old version
    }
  } 

if(DUPFILTER2) {  
  if($nmatch==0 and $dupfirst and $dupids{$lq} and $dupids{$lq} ne $dupfirst) { $nmatch=-1; $ndupdrop++; }
}
  if($nmatch==0) {
    if($OUTSPANTAB) { 
      my ($tbit,$taln,$tidn,$aq,$at,$wq,$wt)= (0) x 19; 
      my $lt="self"; # or lq, or use blast-self scores?
      $aq= $aasize{$lq}||0; $aq *=3;  # 3*convert to cds-size
      $wq= $trsize{$lq}||0; 
      $at=$aq; $wt=$taln=$tidn=$tbit= $wq; # or aq
if(ADDSPANS1605) {
      print $OUTH join("\t",qw(Qid Qclen Qtlen Tid Tclen Ttlen Align Ident Bits Qspan Tspan))."\n" unless($head++);
      print $OUTH join("\t",$lq,$aq,$wq,$lt,$at,$wt,$taln,$tidn,$tbit, "1-$wq","1-$wt")."\n"; 
} else {      
      print $OUTH join("\t",qw(Qid Qclen Qtlen Tid Tclen Ttlen Align Ident Bits))."\n" unless($head++);
      print $OUTH join("\t",$lq,$aq,$wq,$lt,$at,$wt,$taln,$tidn,$tbit)."\n"; 
}
      $validids{$lq}++; $validids{$lt}++;
    } else {
    
    }
  }
  
  %bspans=();
  return($nmatch); # ,$ndupdrop
}




=item sub classifyFullTr
  
  further classing, after -CDSALIGN with identityclass()/classifytr(),
  using -noCDS (ie. full self-mrna-blastn table of aligns, CDS-main/alt classified)
  
  add classifyTrFull/TrFrag/ that adds tr overlaps outside(?) of CDS x CDS:
  separate CDS v UTR overlaps, extra locus classes depend on how much of which type
    .. need CDS qualities: hoscore (aablast), aasize, %CDS (coding pot)
    .. trfrag overlaps CDS of trgood > drop/relcass as partof or altof trgood.
      .. unless trfrag.UTR only over trgood.CDS and trfrag.CDS has good code
    .. trfrag overlaps only UTR of trgood, depends on trfrag coding quals: other locus or partof trgood
      .. trfrag.CDS overlaps trgood.UTR, trfrag > partof:trgood unless has strong code quals
      .. trfrag.UTR only over trgood.UTR, trfrag > other locus (code/noncode?)
    .. trfrag overlaps trfrag, depends on CDS v UTR, code quals

=item sub overlocusclass

  temp in asmrna_dupfilter3_overlocsub.pm

  inputs (via what? overLocusclass instead of identityclass? ) : 
    1. cds.qual with sizes cdslen,aaqual,trlen,cds-offs, and maybe cds-codepot (Code/Noncode)
    2. self-blastn btall table for all mrna x mrna, excluding already classified main x alt set,
        i.e. only diff-locus overlaps
    3. aablast btall table of ref-homol scores [option?]
  
=cut

#not this way# 
# require "asmrna_dupfilter3_overlocsub.pm";
#.........
# asmrna_dupfilter3_overlocsub.pm
# new subs for  evigene/scripts/rnaseq/asmrna_dupfilter3.pl

# use strict;
# use warnings;
#?? package main;

=item classifyFullTr
  from asmrna_dupfilter3.pl:sub classifytr()

  further classing, after -CDSALIGN with identityclass()/classifytr(),
  using -noCDS (ie. full self-mrna-blastn table of aligns, CDS-main/alt classified)
  inputs from overLocusclass
    
  add classifyTrFull/TrFrag/ that adds tr overlaps outside(?) of CDS x CDS:
  separate CDS v UTR overlaps, extra locus classes depend on how much of which type
    .. need CDS qualities: hoscore (aablast), aasize, %CDS (coding pot)
      
      readSizes() collects from aa/cds/mrna.qual: 
        aasize{id}, aaqual{}, trsize{}, oids{}/oidsof{}, cdsoff{} also now
      
    .. trfrag overlaps CDS of trgood > drop/relcass as partof or altof trgood.
      .. unless trfrag.UTR only over trgood.CDS and trfrag.CDS has good code
    .. trfrag overlaps only UTR of trgood, depends on trfrag coding quals: other locus or partof trgood
      .. trfrag.CDS overlaps trgood.UTR, trfrag > partof:trgood unless has strong code quals
      .. trfrag.UTR only over trgood.UTR, trfrag > other locus (code/noncode?)
    .. trfrag overlaps trfrag, depends on CDS v UTR, code quals

=item overlocusclass

  from asmrna_dupfilter3.pl:sub identityclass()
  
  inputs (via what? overLocusclass instead of identityclass? ) : 
    1. cds.qual with sizes cdslen,aaqual,trlen,cds-offs, and maybe cds-codepot (Code/Noncode)
    2. self-blastn btall table for all mrna x mrna, excluding already classified main x alt set,
        i.e. only diff-locus overlaps
    3. aablast btall table of ref-homol scores [option?]

=cut

## debug add for asmrna_dupfilter3_overlocsub.pm, dont need these
sub _minb { my($x,$y)=@_; return ($y < $x) ? $y : $x; }
sub _maxb { my($x,$y)=@_; return ($y > $x) ? $y : $x; }

use constant { kAATINY => 1, kAAUTRBAD => 2, kAAUTRPOOR => 4, kAADUP => 8, kAAGAPS => 16, 
              kNONCODE => 32, }; # kNONCODE == cds.qual Noncode flag
              # kDUPEXONS => 64 ?? for cullExonEq
## eor flags to clear others
use constant NOTTINY => kAADUP + kAAGAPS + kNONCODE; # - kAATINY - kAAUTRPOOR - kAAUTRBAD - kAAGAPS
use constant NOTPOORBAD => NOTTINY + kAATINY; # kAATINY + kAADUP + kAAGAPS + kNONCODE; # - kAAUTRPOOR - kAAUTRBAD - kAAGAPS
use constant NOTPOOR => NOTPOORBAD + kAAUTRBAD; # kAATINY + kAADUP + kAAUTRBAD + kAAGAPS + kNONCODE; # - kAAUTRPOOR

sub trevidence {
  my($tid,$hiid)= @_; # ,$cla,$qid,$pal
  $hiid||=0; ## $hiid= ( $pid >= $PHI ); ## eg >= 99% ident
  
  my $aw= $aasize{$tid} || 0;
  my $tqual= $aaqual{$tid} || ""; #? parse for "aasize,pcds,aaqual"
  ## add if there:
  my $tcode= $sizeval{$tid}{'codepot'}||""; # $tcode =~ /^Noncode/ set what? kAATINY+kAAUTRBAD
  my $tbits= $aablast{$tid} || "0,0";
  my($tbscore,$tbref)=split",",$tbits;
  $tcode=0 if($tbscore>0);
  
  my $mapqloc= ($nineqgene) ? $gmapqual->{$tid} : 0;  
  my($mapqual,$maploc)= ($mapqloc)?(split"\t",$mapqloc):(0,0);

  my($tbrscore,$tbrefbest)= ($tbref) ? (split",", $aablastref{$tbref}) : (0,0); 
  my $rbits2= $aablastref{$tid} || "0,0";
  my($tbrscore2,$tbref2)= split",", $rbits2; # this is revhash tid => score,ref refbest

  my $ispoor= ($aw < $AAMIN or ($tqual =~ m/partial/ and $aw < $AAPART))?kAATINY:0;
  my $tw= $trsize{$tid} || 1; 
  my $pcds= int(300*$aw/$tw);  my $butr= ($tw - 3*$aw); # okutr if butr<=300p, modify BAD_CDSUTR for small aw
  
  if($tqual =~ m/gapbad/) { $ispoor |= kAAGAPS; }  # gaps discrepancy w/ aaqual utrbad and pcds utrbad 
  if($butr >= $MINUTR) { # ~ 300b "fixed" average utr sizes
    # should adjust %cds bad for cds size, for large cds> 10k, pcds >=90%, for small cds < 300, pcds may be 33%
    if($tqual =~ m/utrbad/ and $noRESET_CDSUTR) { $ispoor |= kAAUTRBAD; } 
    elsif($tqual =~ m/utrpoor/ and $noRESET_CDSUTR) { $ispoor |= kAAUTRPOOR; }
    else { $ispoor |= ($pcds < $BAD_CDSUTR)? kAAUTRBAD: ($pcds < $OK_CDSUTR)? kAAUTRPOOR:0; }
    # if($butr <= 300) {} elsif($pcds < $BAD_CDSUTR) { $ispoor |= kAAUTRBAD; } elsif($pcds < $OK_CDSUTR) { $ispoor |= kAAUTRPOOR; } 
  }
  if($tcode and not $tbscore) { 
    # $ispoor |= kAATINY+kAAUTRBAD+kNONCODE if($tcode =~ /^Noncode/); # leave off AATINY+AAUTRBAD?
    $ispoor |= kNONCODE if($tcode =~ /^Noncode/); #?
    }
  
  my $MINBLASTSCORE= 60; # aablast only? bitscore always?
  my $keepdrop=""; # not used here..
  unless( $tbscore == 0 or $tbits=~/^0,0/) {
    my $risbest= (($tbrefbest eq $tid) or ($tbscore >= $tbrscore))?1:0;# was ($tbrefbest eq $tid), allow 2+ same score
    my $risgood= ($tbrscore > $MINBLASTSCORE and $tbscore >= 0.90 * $tbrscore)?1:0;
    my $risok=   ($tbrscore > $MINBLASTSCORE and $tbscore >= 0.50 * $tbrscore)?1:0;
       $risok |= ($tbref2 and $tbrscore2 >= $MINBLASTSCORE)?1:0;
    if($risok and $risgood) { $risgood=0 if($tbrscore2 > $tbrscore); } #?
    $ispoor = $ispoor & NOTTINY if($risgood); ## FIXME1712: clear kAATINY for $risok/risgood
    my $isaadup=($ispoor & kAADUP);
    if($risbest) { 
      $tbits.=",refbest"; $keepdrop.="okay:refbest,"; $ispoor=0 unless($isaadup); }
    elsif($risgood) { # was 3rd
      $tbits.=",refgood"; unless($isaadup){ $keepdrop.="okay:refgood,"; $ispoor=0;} } # maybe2 ok ?
    elsif($risok) {  # was 2nd ? why superceed refgood? 2nd ref match?
      $tbits.=",refok"; unless($isaadup){ $keepdrop.="okay:refok,"; $ispoor=0;} }
  }
  
  if($ispoor > kAATINY and not $hiid) { ## $cla !~ /althi|parthi/ == $hiid
    if($aw >= $AAMINBAD) { $ispoor = $ispoor & NOTPOORBAD; } ##  ^ (kAAGAPS+kAAUTRPOOR+kAAUTRBAD);  ##?? ^ XOR bad
    elsif($aw >= $AAMINPOO) { $ispoor = $ispoor & NOTPOOR; } ##  ^ (0+kAAUTRPOOR);
  }

  $tbits.= ",chrmap:$mapqual,$maploc" if($maploc and $mapqual); #? skip her?
  $tbits.= ",pflag:$ispoor"; # DEBUG info, but keep
  $tbits="aaref:$tbits" unless($tbits=~/^0,0/);

  return($ispoor,$tqual,$tbits); # .. others
}


sub classifyFullTr {
  my($tid,$cla,$qid,$pidal)= @_;
  my($aw,$ispoor,$tqual,$tbits,$mapqual,$maploc)= (0) x 9;
  my $keepdrop="";

  #   use constant { kAATINY => 1, kAAUTRBAD => 2, kAAUTRPOOR => 4, kAADUP => 8, kAAGAPS => 16, };
  #   use constant NOTPOORBAD => kAATINY + kAADUP; # - kAAUTRPOOR - kAAUTRBAD - kAAGAPS
  #   use constant NOTPOOR => kAATINY + kAADUP + kAAUTRBAD + kAAGAPS; # - kAAUTRPOOR
  $qid||="0"; $pidal||="0";
  my ($aclass,$fclass) =  split"/",$cla; ## ??
  $cla= $aclass; #?

if(1) {
  ($ispoor,$tqual,$tbits)= trevidence($tid,$cla =~ m/althi|parthi/);

  if($cla =~ /cull/) { $ispoor |= kAADUP; } # # add new flag kDUPEXONS ?
  if($cla =~ /althi1|part/ and $pidal =~ m/altmap\d+xeq/) { $ispoor |= kAADUP; } # kDUPEXONS?

  my $aaclus= $aacluster{$tid};  
  if($aaclus) {
    my($aamainid,$aamainpi)= split",", $aaclus; # ="$mainid,$pi";
    if($aamainid eq $tid) {  $aamainid=""; $aamainpi=0; }
    if($aamainpi > $AADUP_IDENT) { # TEST aacluster added info on dup prots; AADUP_IDENT=98 default
      $cla.="a2" unless($cla=~/hi[012]/); # hi1/a2 redundant info
      $tbits.=",aadup:$aamainid";
      $ispoor |= kAADUP if($cla =~ /althi/); # only for /althi|parthi/
    }
  } else { $aaclus ||= "0,0"; }
  
  if($tbits =~ m/aaref:/) {
    my $isaadup=($ispoor & kAADUP);
    if($tbits=~/,refbest/) { $keepdrop.="okay:refbest,"; $ispoor=0 unless($isaadup); }
    elsif($tbits=~/,refgood/ and not $isaadup) { $keepdrop.="okay:refgood,";  $ispoor=0;}
    elsif($tbits=~/,refok/ and not $isaadup) { $keepdrop.="okay:refok,";  $ispoor=0; }
  }
  
} else {  
  
  $aw= $aasize{$tid} || 0;
  $tqual= $aaqual{$tid} || ""; #? parse for "aasize,pcds,aaqual"
  $tbits= $aablast{$tid} || "0,0";
  my($tbscore,$tbref)=split",",$tbits;
  
  my $mapqloc= ($nineqgene) ? $gmapqual->{$tid} : 0;  
  ($mapqual,$maploc)= ($mapqloc)?(split"\t",$mapqloc):(0,0);

  my($tbrscore,$tbrefbest)= ($tbref) ? (split",", $aablastref{$tbref}) : (0,0); 
  my $rbits2= $aablastref{$tid} || "0,0";
  my($tbrscore2,$tbref2)= split",", $rbits2; # this is revhash tid => score,ref refbest

  my $aaclus= $aacluster{$tid} || "0,0";  
  my($aamainid,$aamainpi)= split",", $aaclus; # ="$mainid,$pi";
  if($aamainid eq $tid) {  $aamainid=""; $aamainpi=0; }
  
  my $ispoor= ($aw < $AAMIN or ($tqual =~ m/partial/ and $aw < $AAPART))?kAATINY:0;
  my $tw= $trsize{$tid} || 1; 
  my $pcds= int(300*$aw/$tw);  my $butr= ($tw - 3*$aw); # okutr if butr<=300p, modify BAD_CDSUTR for small aw

  if($tqual =~ m/gapbad/) { $ispoor |= kAAGAPS; }  # gaps discrepancy w/ aaqual utrbad and pcds utrbad 
  if($butr >= $MINUTR) { # ~ 300b "fixed" average utr sizes
    # should adjust %cds bad for cds size, for large cds> 10k, pcds >=90%, for small cds < 300, pcds may be 33%
    if($tqual =~ m/utrbad/ and $noRESET_CDSUTR) { $ispoor |= kAAUTRBAD; } 
    elsif($tqual =~ m/utrpoor/ and $noRESET_CDSUTR) { $ispoor |= kAAUTRPOOR; }
    else { $ispoor |= ($pcds < $BAD_CDSUTR)? kAAUTRBAD: ($pcds < $OK_CDSUTR)? kAAUTRPOOR : 0; }
    # if($butr <= 300) {} elsif($pcds < $BAD_CDSUTR){ $ispoor |= kAAUTRBAD; } elsif($pcds < $OK_CDSUTR) { $ispoor |= kAAUTRPOOR; } 
  }
  
  if($pidal =~ m/altmap\d+xeq/ and $cla =~ /althi1|part/) {
    $ispoor |= kAADUP;  # add new bitflag kDUPEXONS ?
  }
  
  if($aamainpi > $AADUP_IDENT) { # TEST aacluster added info on dup prots; AADUP_IDENT=98 default
    $cla.="a2" unless($cla=~/hi[012]/); # hi1/a2 redundant info
    $tbits.=",aadup:$aamainid";
    $ispoor |= kAADUP if($cla =~ /althi/); # only for /althi|parthi/
  }

  my $MINBLASTSCORE= 60; # aablast only? bitscore always?
  unless( $tbscore == 0 or $tbits=~/^0,0/) {
    my $risbest= (($tbrefbest eq $tid) or ($tbscore >= $tbrscore))?1:0;# was ($tbrefbest eq $tid), allow 2+ same score
    my $risgood= ($tbrscore > $MINBLASTSCORE and $tbscore >= 0.90 * $tbrscore)?1:0;
    my $risok= ($tbref2 and $tbrscore2 >= $MINBLASTSCORE)?1:0;
    if($risok and $risgood) { $risgood=0 if($tbrscore2 > $tbrscore); } #?
    $ispoor = $ispoor & NOTTINY if($risgood); ## FIXME1712: clear kAATINY for $risok/risgood
    my $isaadup=($ispoor & kAADUP);
    if($risbest) { 
      $tbits.=",refbest"; $keepdrop.="okay:refbest,"; $ispoor=0 unless($isaadup); }
    elsif($risgood) { # was 3rd
      $tbits.=",refgood"; unless($isaadup){ $keepdrop.="okay:refgood,"; $ispoor=0;} } # maybe2 ok ?
    elsif($risok) {  # was 2nd ? why superceed refgood? 2nd ref match?
      $tbits.=",refok"; unless($isaadup){ $keepdrop.="okay:refok,"; $ispoor=0;} }
  }
  
  if($ispoor > kAATINY and $cla !~ /althi|parthi/) {  
    if($aw >= $AAMINBAD) { $ispoor = $ispoor & NOTPOORBAD; } ##  ^ (kAAGAPS+kAAUTRPOOR+kAAUTRBAD);  ##?? ^ XOR bad
    elsif($aw >= $AAMINPOO) { $ispoor = $ispoor & NOTPOOR; } ##  ^ (0+kAAUTRPOOR);
  }
}

      ## classifyFullTr() parts .. do AFTER find cross-locus alignments
  my $fclact="";
  if($fclass) {  
    
    ## use ispoor ...
    # my $thoval= $aablast{$tid} || 0;  # my($qhoscore,$qhoref)=split",",$thoval;
    # my $pcds= ($tw>0) ? 100*$tc/$tw : $OK_CDSUTR;
    # my $tcdsisbad=0;
    # if($thoval) { $tcdsisbad=0; }
    # elsif($tcode=~/^Noncode/) { $tcdsisbad=$qcode; }
    # elsif($tc < 2*$MINCDS and $pcds <= $BAD_CDSUTR) { $tcdsisbad=2; }
    # elsif($tc < $MINCDS and $pcds < $OK_CDSUTR) { $tcdsisbad=1; }
    
    ## oneloc == cdsovcds or (utrovcds and tcdsisbad)
    ## twoloc == (utrovutr or cdsovutr) and not tcdsisbad
    ## ambig  == others?
    if($fclass =~ /cdsovcds/ or ($fclass =~ /utrovcds/ and $ispoor)) { 
      if($fclass=~/part/ or $ispoor) { $fclact="drop.$fclass"; } # frag/part
      else { $fclact="alt.$fclass"; } ## large set here
    } elsif($fclass =~ /utrovutr|cdsovutr/ and not $ispoor) {
      $fclact="locus.$fclass"; # keep? use orig locus class ?
    } elsif($fclass =~ /altof|mainof/) { #  and not $ispoor
      $fclact=$fclass; # "altof"; #?  mainof > mainlocus ??
    } else {
      $fclact="locusmaybe.$fclass"; # ispoor? includes .notover1
    }
  }
  $fclact ||= "notover2"; #?
  
  #  all main+noclass/notover (no utr overlaps) that are cds-ispoort (Noncode) should be kept as ncRNA
  
  my $claret="$cla/$fclact"; # return this?
  if($cla =~ /parthi|frag0aa/ or $fclact =~ /^drop/) {
    $keepdrop.= "drop";
    # $cla=$fclact unless($cla =~ /parthi|frag/);
  
  } elsif($fclact =~ /^locus/) { # locusmaybe?
    $keepdrop.= "okay"; #? includes notover1 + ispoor == noncode?
    if($ispoor and $cla =~ /main|noclass/) { $claret="noncode/$fclact"; }
    # $cla=$fclact;
    
  } elsif($cla =~ /althi/ or $fclact =~ /^alt/) {
    $keepdrop.= ($ispoor)?"drop":"okay";  # ispoor vs main size?
    # $cla=$fclact unless($cla =~ /alt/);
    
  } elsif($cla =~ /main/) { # mainlocus? other term for best of overlocus aligns?
    if($ispoor and $fclact =~ /notover1/) { $ispoor=0; $claret="noncode/$fclact"; }
    $keepdrop.= ($ispoor)?"drop":"okay";  # ??
    
  } elsif($cla =~ /noclass/) {
    if($ispoor and $fclact =~ /notover1/) { $ispoor=0; $claret="noncode/$fclact"; }
    $keepdrop.= ($ispoor)?"drop":"okay"; #? dont drop? ditto to altmid
    
  } else { # other altmid/low ; part ?
    $keepdrop.= ($ispoor)?"drop":"okay"; #? dont drop this class, ispoor maybe uniq prot: maybeok?
  }
  
  
  my $okay;
  if($keepdrop =~ /drop/) {
    if($keepdrop =~ /okay/) { $okay= "maybeok"; } else { $okay= "drop"; }
  } else {
    $okay= "okay";
  }
  
  # now added in trevidence() .. flag to skip this?
  unless($tbits=~/chrmap:/) { $tbits.= ",chrmap:$mapqual,$maploc" if($maploc and $mapqual); }
  unless($tbits=~/pflag:/) { $tbits.= ",pflag:$ispoor"; } # DEBUG info, but keep
  unless($tbits=~/^aaref:/) { $tbits="aaref:$tbits" unless($tbits=~/^0,0/); } 
  
  ## add feq: back for overlocus cullExonEq() test??  
  if(defined $eqflag{$tid} and not($tbits =~ /feq:/)) { 
    my $eqfl= $eqflag{$tid}{$qid}||""; 
    if($eqfl) { $eqfl="$qid/$eqfl,"; }
    my @q= grep{ $_ ne $qid } sort keys %{$eqflag{$tid}}; 
    $eqfl .= join",",map{ "$_/".$eqflag{$tid}{$_} }@q;  
    $tbits.= ",feq:$eqfl";
  }
  
  ## maybe add for output cds.qual stats tid/qid (debug?)
  return (wantarray) ? ($tid,$okay,$claret,$qid,$pidal,$tqual,$tbits) : $okay;
}

=item example overlocusclass

  see asmrna_dupfilter3_overlocsub.pm

eg5:evg12aedesmpub.m4class
Aedesg12bEVm000320t19   okay    althi/alt.cdsovcds      Aedesg12bEVm000215t3    100/7/. 1171,73%,complete       aaref:aaref:2110,AAEL010606-PA,pflag:0,pflag:0
Aedesg12bEVm000320t20   okay    althi/alt.cdsovcds      Aedesg12bEVm000215t3    100/7/. 1171,72%,complete       aaref:aaref:2110,AAEL010606-PA,pflag:0,pflag:0
    100/7 << tiny cds over 7%, should it be called alt?
.. these are oid-alts, however, should use that
Aedesg12bEVm000215t1	Aedesg1EVm000222t1	Aedesg12bEVm000215	1	main	2232,90%,complete	100/100/./altmap100xeq	aaref:4532,AAEL018134-PA,refgood,chrmap:100a,100i,6699l,9x,supercont1.241:1295228-1328778:-,pflag:0,feq:Aedesg1EVm000222t2/altmapxe100.100,
Aedesg12bEVm000215t2	Aedesg1EVm000222t2	Aedesg12bEVm000215	2	althi1maybe	2232,88%,complete	100/100/./altmap100xeq	aaref:4533,AAEL018134-PA,refbest,chrmap:100a,100i,6699l,9x,supercont1.241:1295228-1328778:-,pflag:8,feq:Aedesg1EVm000222t1/altmapxe100.100,
Aedesg12bEVm000215t3	Aedesg1EVm000222t3	Aedesg12bEVm000215	3	althi	2084,85%,partial3	100/94/.	aaref:3944,AAEL018134-PA,chrmap:100a,100i,6160l,9x,Spl:5%,supercont1.490,supercont1.241:1295550-1320691:-,pflag:0
>> diff locus? m000215t3 is split-map problem case, 0320t19 may be accurate call, diff from  m000215t[12]
>> should drop any t19 problem cases, larger alts ok
Aedesg12bEVm000320t19	Aedesg1EVm000222t4	Aedesg12bEVm000320	19	cullalthi1	1171,73%,complete	100/92/./altmap52	aaref:2110,AAEL010606-PA,chrmap:99a,100i,3512l,9x,supercont1.490:478821-528287:-,pflag:0,feq:Aedesg1EVm000328t1/altmap52.36,Aedesg1EVm000222t3/altpar9.0.0,Aedesg1EVm000222t5/altmapxe99.100,Aedesg1EVm000328t10/altmap84.73,Aedesg1EVm0
Aedesg12bEVm000320t1	Aedesg1EVm000328t1	Aedesg12bEVm000320	1	main	2012,81%,complete	100/95/.	aaref:3698,AAEL010606-PA,refgood,chrmap:99a,100i,6038l,22x,Spl:2%,supercont1.1301,supercont1.490:463365-690982:-,pflag:0
Aedesg12bEVm000320t2	Aedesg1EVm000328t2	Aedesg12bEVm000320	2	althi	2008,81%,complete	100/93/.	aaref:3731,AAEL010606-PA,refgood,chrmap:97a,100i,6026l,21x,supercont1.490:463365-690982:-,pflag:0
Aedesg12bEVm000320t3	Aedesg1EVm000328t3	Aedesg12bEVm000320	3	althi	1999,81%,complete	100/95/.	aaref:3846,AAEL010606-PA,refbest,chrmap:100a,100i,5999l,22x,supercont1.490:463365-690982:-,pflag:0
Aedesg12bEVm000320t12	Aedesg1EVm000328t12	Aedesg12bEVm000320	12	althi	1251,74%,complete	100/93/.	aaref:2342,AAEL010606-PA,chrmap:99a,100i,3752l,11x,supercont1.490:478821-550496:-,pflag:0
.. cds.qual
Aedesg12bEVm000215t1	6699	0	2232,90%,complete	7387	686-7384:.	Aedesg1EVm000222t1,tidbaedes2sr6feo2ridk97Loc29099
Aedesg12bEVm000215t2	6699	0	2232,88%,complete	7602	632-7330:.	Aedesg1EVm000222t2,tidbaedes2sr4mao2ridk41Loc45077
Aedesg12bEVm000215t3	6160	92	2084,85%,partial3	7299	1047-7298:.	Aedesg1EVm000222t3,aedes2sr4ma4dsoapk25loc10166t3
Aedesg12bEVm000320t12	3752	4	1251,74%,complete	5025	31-3786:.	Aedesg1EVm000328t12,aedes2fe6nm4dsoapk27loc22685t20
Aedesg12bEVm000320t19	3512	4	1171,73%,complete	4770	16-3531:.	Aedesg1EVm000222t4,aedes2fe6nm4dsoapk27loc22685t32

  
=cut

=item defaults

$AAMIN =$ENV{aamin}||30; #was 40;   # for aacomplete, utrok

$TINYALN = 20; # %align, was 25; was $ENV{mina}||50; ignore less than this == MINAL
$TINYALNBASES= $ENV{MINALIGN}||90; # was $MINALIGN= 90; # REUSE//NOT USED NOW; use just TINYALN; change to several levels
$TINYALNCDSLEN= 3*$TINYALNBASES;
$MINCDS = $ENV{MINCDS} || 90; # or 3*MINAA  

$ALTFRAG= 0.5;

$PHIALN=$ENV{ahi} || 98; # was 99; was 65 !!;  DROP THIS; use mina?
$NHIALN=$ENV{nhi} || 9; # what? for short cds cancel few base diff < PHIALN

$PHI = $ENV{phi} ||99; 
$PMID= $ENV{pmid}||90; 
$PLOW= $ENV{plow}||80;  

$OK_CDSUTR_default= 60;
$BAD_CDSUTR_default= 30; # % CDS/trlen

=cut

sub overlocus_eqgene {
  my($qd,$td,$pal)=@_;
  my($palmap,$antiflag,$eqflag)=(0,"","");
  
  $palmap= $eqgenes->{$qd}{$td}||0; # note this is pct of td aligned to qd

  ## oids: $palmap= $eqgenes->{$qd}{$td}||0; # note this is pct of td aligned to qd

  ## problem here, when oidof(qid) sameloc oidof(td) wont have eqgene entry??
  ## cancel this, oid-alts do have eqgene entries.. but some dont overlap  
  #   if(0 and %oidsof) {
  #     # $td= $oids{$td}||$td; $qd= $oids{$qd}||$qd; 
  #     ## this isn't good enough test, paralts diff loci show up w/ same oid
  #     ## need glocation test of sameoid, add to eqgene tables?
  #     my ($qoid)= split",", $oidsof{$qd};
  #     my ($toid)= split",", $oidsof{$td};
  #     my($qog,$tog)=map{ my $g=$_; $g=~s/t\d+$//; $g; } ($qoid,$toid);
  #     if($qog and $qog eq $tog) {
  #       $antiflag .="/altoid";
  #       $palmap= $pal if($palmap<$pal); # if($palmap <= 0) ?
  #       # return( $pal,$palmap,$antiflag );
  #     }
  #   }
  
  ## if($neqgene>0) 

  ## simple way:
  # if($palmap>0) { $antiflag .="/altmap$palmap"; $eqflag="altmap$palmap";
  #	if(1 && $palmap > $pal) { $pal=$palmap; } ##? not this change EQGENE_OVERRIDES_ALN
  #}

  ## hard way:
  if($palmap>0) { 
    my $xeq= $eqexons->{$qd}{$td}||0; # only for ($neqexon>0)
    my $XPHI = 95; my $XPLO= 3;
    # xeq care about (a) xeq>= identity == not alt but redundant, 
    #   (b) xeq <= noxover, not alt but paralog maybe, (c) middle = usual alt
    
    #* set aln=0/TINYALNBASES when pal=0/TINYALN
    #*? change 'paralt' to 'altpar' to avoid other evg parse problems?
    if($neqexon>0 and $xeq >= $XPHI) { 
      #n# $antiflag .="/altmapxe$palmap.$xeq"; #?? ugh, want this way for other uses
      $antiflag .="/altmap$palmap"."xeq/$qd"; #orig..
      $eqflag="altmapxe$palmap.$xeq"; # for report? or use only antiflag?
      $pal=$palmap if($palmap >= $PHI and $pal < $PMID);  # FIXME bad annot
      }
    elsif($neqexon>0 and $xeq < $XPLO) { # problems here.. dont reset pal yet
      $antiflag .="/altparx$pal/$qd";  #later? $pal=3; $aln=13;
      $eqflag="altparx$palmap.$xeq"; # for report
      }
    else { 
      (my $qg=$qd)=~s/t\d+$//;  
      unless($td=~/^$qg/) {
        $antiflag .="/altmap$palmap/$qd"; #orig .. add .$xeq";
        $eqflag="altmap$palmap.$xeq"; # for report
        $pal=$palmap if($palmap >= $PHI and $pal < $PMID);  # FIXME bad annot
        }
      }
    }
    
  # elsif($pal>0 and exists( $eqgenes->{$td}) and exists( $eqgenes->{$qd})) { # defined or exists??
  #   ## bad here for nomap alts ** why isnt defined/exists working? for ids not in eqgene
  #   if($gmapqual->{$qd} and $gmapqual->{$td}) {
  #     my $revmap= $eqgenes->{$td}{$qd}||0; 
  #     ## unmapped gene cases are problem here.. need what? defined ($eqgenes->{td}) 
  #     $antiflag .="/altpar$pal" if($revmap<3);  # $pal=3; $aln=13; #? change later?
  #     $eqflag="altpar$pal.$palmap.$revmap"; # for report
  #   } else {
  #     #lots# warn "#DBG: bad eqgene $qd,$td\n" if($debug);
  #   }
  # }
    
  ##}

  $eqflag{$td}{$qd}=$eqflag if($eqflag);
  return( $pal,$palmap,$antiflag,$eqflag );
}


use constant minXEQ => 99; # exon-equal, high to ensure keep valid alts 
use constant minXCP => 95; # cds-align min for cullExonEq
use constant minXAL => 90; # chr-align min for cullExonEq

=item cull1ExonEq 
  from trclass2mainalt.pl, change for 1 compare

  bugs:
  FIXME: this is culling differing alts, at least by aasize,maploc, gmap-nexon, ..
  use note attrs more. chrmap:
  these are md-cullers, but not same alt form as longer alts
  need also note param of md culler 
   Anofunz4kEVm000028t7	okay	mainlocus/notover2	Anofunz4kEVm000028t1	100/68/.	2766,93%,complete	
     aaref:4249,AGAP002523-PA,refbest,chrmap:99a,99i,8297l,9x,KB669425:23470-33108:-,pflag:0
   Anofunz4kEVm000028t22	maybeok	parthi/altof	Anofunz4kEVm000028t1	100/70/.	2064,69%,complete	
     aaref:4022,AGAP002523-PA,refgood,chrmap:100a,99i,6195l,8x,KB669425:25506-33108:-,pflag:0
   
  Anofunz4kEVm000028t4	drop	althicull/altof	Anofunz4kEVm000028t1	100/59/.	3292,86%,complete	
    aaref:3773,AGAP002523-PA,refok,chrmap:97a,99i,9879l,15x,Spl:28%,KB668673,KB669425:24898-33108:-,pflag:0,feq:Anofunz4kEVm000028t22/altmapxe100.100
  Anofunz4kEVm000028t5	drop	althicull/altof	Anofunz4kEVm000028t1	100/76/.	2786,93%,complete	
   aaref:3764,AGAP002523-PA,refok,chrmap:99a,99i,8361l,10x,KB669425:23470-33108:-,pflag:0,feq:Anofunz4kEVm000028t22/altmapxe100.100,Anofunz4kEVm000028t6/altmapxe99.100,Anofunz4kEVm000028t7/altmapxe99.100
  Anofunz4kEVm000028t6	drop	althicull/altof	Anofunz4kEVm000028t1	100/68/.	2785,93%,complete	
    aaref:4238,AGAP002523-PA,refgood,chrmap:99a,99i,8354l,10x,KB669425:23470-33108:-,pflag:0,feq:Anofunz4kEVm000028t7/altmapxe99.100

  ** maybe two loci here, diff aarefs, splitmaps, diff scafs
  >> revised cull1x better, only 028t5 is drop.cull, 028t6 becomes new main (refgood),
     028t22 becomes maybeok	parthi/altof t1/t6
     t1,t2 longer aa, but smaller aaref:hoscore become alts of t6 ?? could be poor AGAP prot
     
=cut

sub cull1ExonEq {  
  my($md,$mdcl,$td,$tcl,$feq)=@_;  # ,$note
  my($culled,$galn)=("",0); 
  return "" if($mdcl=~/cull/);

  #now from caller: my($feq)= ($note=~/feq:([^;:\s]+)/)?$1:0;
  
  ## overlocus_eqgene now returns these
  ##    $antiflag .="/altmap$palmap.xeq"; #orig..
  ##    $eqflag{$td}{$qd}="altmapxe$palmap.$xeq"; # for report? or use only antiflag?

  ## my $note= $notes->{$td};
  my($tispoor,$tqual,$note) = trevidence($td,1);    
  my($mispoor,$mqual,$mnote)= trevidence($md,1);   

  return "" if( ($note=~/refbest/) or not ($tcl=~/althi|parthi/) 
        or ($mispoor>0 and $tispoor==0) );
  
  # pick out all of note parts, mnote also ??
  my($tgaln,$txn,$tspl,$mgaln,$mxn,$mspl,$cmap)=(0) x 9;
  ($cmap)= $note =~ m/chrmap:([^;:\s]+)/; # ,scaff[:]span[:]or .. not needed 
    ($tgaln)= ($cmap =~ m/(\d+)a,/)?$1:0;
    ($txn)  = ($cmap =~ m/(\d+)x,/)?$1:0;
    ($tspl) = ($cmap =~ m/Spl:(\d+)/)?$1:0;
  ($cmap)= $mnote =~ m/chrmap:([^;:\s]+)/; 
    ($mgaln)= ($cmap =~ m/(\d+)a,/)?$1:0;
    ($mxn)  = ($cmap =~ m/(\d+)x,/)?$1:0;
    ($mspl) = ($cmap =~ m/Spl:(\d+)/)?$1:0;
  
  my @xe= grep /altmapxe/, split",",$feq; # 1 only here? xd == md
  for my $xe (@xe) {
    my($xcp,$xev)= $xe =~ m/altmapxe(\d+).(\d+)/; # ugh, diff syntax now from overloc
    next unless($xev and $xev >= minXEQ and $xcp >= minXCP);
    next if($tgaln < minXAL or ($mxn and $txn > $mxn) or $tspl > 9);
    $culled= $tcl."cull"; # was "cull$tcl/$xe";
    return $culled; # or last;
  }

  return $culled;  
}


sub overlocusclass {
  my($outh, $infile, $insorted, $iscdsalign)= @_;
  $iscdsalign||=0; ## default NOT $IS_CDSALIGN

  ## ?? infile == blast table -outfmt 7, this sort won't work..
  ## redo alntab to include offsets $qab,$qae,$tab,$tae ?
  
  my $ALNSORTORD='-k2,2nr -k7,7nr -k3,3n -k6,6n -k1,1 -k4,4';  #add vQtlen/3 before vTt/6 IDs to order ties
  unless($iscdsalign) { # sort tlen, not clen 
    # .. not sure what is right, want cds to play role in best choice
    # .. input blast align k7 is for tlen, not clen, but want to choose best by long clen > long tlen,talign
    #t1.NO: $ALNSORTORD='-k3,3nr -k7,7nr -k2,2nr -k6,6nr -k1,1 -k4,4';  # ^Qtlen,^Align,^Qclen,vTtlen,Qid,Tid
    $ALNSORTORD='-k2,2nr -k7,7nr -k3,3n -k6,6n -k1,1 -k4,4';  # test same as iscdsalign  
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
    my($qd,$td,$qc,$qw,$qcoff,$qcode,$qcoffb,$qcoffe, $qab,$qae,$qspan,
       $tc,$tw,$tcoff,$tcode,$tcoffb,$tcoffe, $tab,$tae,$tspan, $tor, 
       $aln,$iden,$bits,$pidn,$mis,$gap) = (0) x 30;
    my($isfrag,$aclass)=("") x 4;
    my($alnmax,$pid,$pal,$samesize)= (0) x 10;

    ## -bits == reverse align : bits not used here for scoring.. best choice other than adding column
   
if(0) {    ## no good for sort input..
    ($qd,$td,$bits,$pidn,$aln,$mis,$gap,$qab,$qae,$tab,$tae)= @v[0,1,-1,2,3,4,5, 6,7,8,9]; #blast table ?
         # 6-9 =  q. start, q. end, t. start, t. end, 

    # $iden= _maxb(0,$aln-$mis-$gap); # better? $iden = $pidn/100 * $aln;
    $iden= int($pidn/100 * $aln); # or $aln-$mis-$gap
    if($tab>$tae) { $bits=-$bits; $tor=-1; ($tab,$tae)=($tae,$tab); }
    ($qc,$qw,$qcoff,$qcode) = map{ $sizeval{$qd}{$_}||0 } qw( cdsize  trsize cdsoff  codepot);
    ($tc,$tw,$tcoff,$tcode) = map{ $sizeval{$td}{$_}||0 } qw( cdsize  trsize cdsoff  codepot);
    
} else {
    ($qd,$qc,$qw,$td,$tc,$tw,$aln,$iden,$bits,$qspan,$tspan)= @v;  #*** need blast align start,end vals
    ($qab,$qae)=split "-",$qspan; 
    ($tab,$tae)=split "-",$tspan; 
    $pidn= ($aln<1)? 0: int(0.5+ 100*$iden/$aln); $pidn=100 if($pidn>100);

    ($qcoff,$qcode) = map{ $sizeval{$qd}{$_}||0 } qw( cdsoff  codepot);
    ($tcoff,$tcode) = map{ $sizeval{$td}{$_}||0 } qw( cdsoff  codepot);
}

    # skip to next early on these
    if($td eq "self" or $qd eq $td) { # putspan fix for no-align cases
      $td=$qd; $bestmatch{$qd}="$td,100/100" unless($bestmatch{$qd});
      next;   
    }
    
    $qcoff= $cdsoff{$qd} unless($qcoff);
    $tcoff= $cdsoff{$td} unless($tcoff);

=item tcoff (only?) reading bug
    td == self bug?
  Use of uninitialized value $tcoffe in subtraction (-) at evigene/scripts/rnaseq/asmrna_dupfilter3.pl line 1224, <IN> line 2541139.
  Use of uninitialized value $tcoffb in subtraction (-) at evigene/scripts/rnaseq/asmrna_dupfilter3.pl line 1224, <IN> line 2541139.

=cut

=item special case near identicals
  .. maybe best handled here with swap to best, but need all evidence 1st.

  ## also check if class{$td}/class{$qd} ??
  ## do AFTER 1st tc > qc swap and new tcaln calcs..
  my $nearlysame=(abs($tc-$qc) < 30 and $tcaln >= 0.90*$qc)?1:0;    
  if($nearlysame) { 
    my $doswap=0;
    my($tispoor,$tqual,$tbits)= trevidence($td); # ,$tcla =~ m/althi|parthi/
    my($qispoor,$qqual,$qbits)= trevidence($qd); #,$qcla =~ m/althi|parthi/
    my($tho,$qho)=map{ my($ho)= (m/aaref:(\d+)/)?$1:0; $ho; } ($tbits,$qbits);
    if($tho > 5+$qho) { $doswap=1; } elsif($qho > 5+$tho) { $doswap=0; }
    elsif($tispoor > $qispoor) { $doswap=1; } 
    else { $doswap=0; }
    if($doswap) {
      ($qd,$td)=($td,$qd); 
      ($qc,$qw,$tc,$tw)=($tc,$tw,$qc,$qw); 
      ($qcoff,$qcode,$tcoff,$tcode)=($tcoff,$tcode,$qcoff,$qcode);
      ($qab,$qae,$tab,$tae)= ($tab,$tae,$qab,$qae);
      ($qgd,$tgd)=($tgd,$qgd); 
      ($qcoffb,$qcoffe,$tcoffb,$tcoffe)=($tcoffb,$tcoffe,$qcoffb,$qcoffe);
      ($qcaln,$qualn,$tcaln,$tualn)= ($tcaln,$tualn,$qcaln,$qualn);
    }
  }
  
eg2:evg12aedesmpub.m3class 
Aedesg12bEVm007420t1	okay	main/	Aedesg12bEVm068781t1	100/92/.	473,92%,complete	aaref:969,AAEL011510-PA,refgood,pflag:0
  ^^vv m068781t1 (3 alts) ~= m007420t1 (1 alt), same size, 781t has bettr aaref score, should replace m007420t1 
Aedesg12bEVm068781t1	okay	althi/alt.cdsovcds	Aedesg12bEVm007420t1	100/92/.	473,92%,complete	aaref:975,AAEL011510-PA,refbest,pflag:0
eg3:
Aedesg12bEVm004960t1	okay	main/	Aedesg12bEVm068777t1	99/91/.	605,91%,partial5	aaref:975,AAEL014089-PA,refgood,pflag:0
  ^^vv m068777t1 better aaref than m004960t1, same size, same locus
Aedesg12bEVm068777t1	okay	althi/alt.cdsovcds	Aedesg12bEVm004960t1	99/91/.	605,90%,partial5	aaref:989,AAEL014089-PA,refbest,pflag:0


=cut

    ## Q presumed larger than T (sort), if need swap, do before all other calcs.
    use constant SWAPBIGCDS => 1; 
    if(SWAPBIGCDS && $tc > $qc) { 
      ($qd,$td)=($td,$qd); 
      ($qc,$qw,$tc,$tw)=($tc,$tw,$qc,$qw); 
      ($qcoff,$qcode,$tcoff,$tcode)=($tcoff,$tcode,$qcoff,$qcode);
      ($qab,$qae,$tab,$tae)= ($tab,$tae,$qab,$qae);
      }
 
    my($qgd,$tgd)=map{ my $g=$_; $g=~s/t\d+$//; $g; } ($qd,$td);
    my $isaltof = ($qgd eq $tgd)?1:0; # fixme, fclass == altof
    my $ismainof= ($isaltof and $td=~/t1$/)?1:0; # want this? avoid reclassing t1mains, for now
    
use constant doALTCULL => 1;

    if(not doALTCULL and $isaltof) { 
      # dont care here for known alts??
      ## DO this after below cds-over aclass = althi 
      ## add cullExonEq() test? but need input with gmap feq:ID/altmapxe87.100 .. 
      # if(0 and $neqgene>0) {
      #   my( $palg,$palmap,$antiflagg,$eqflagg ) = overlocus_eqgene($qd,$td,$pal) ; # this sets feq altmapxe..
      #   ## $antiflag .="/altmapxe$palmap.$xeq"."xeq"; 
      #   if($eqflagg =~ m/altmapxe\d/) {
      #     ## also need aaref:val,refbest/good/ok,  $note=~/feq:([^;:\s]+)/,  in antiflag == note
      #     ## aaref from trevidence
      #     my($ispoor,$tqual,$tbits)= trevidence($td,1); # $cla =~ m/althi|parthi/
      #     my($culls)= cull1ExonEq($qd,$td,"althi","$tbits,feq:$qd/$eqflagg"); 
      #     # my $culls= ($CULLXEQ) ? cullExonEq($md,\@ad,\%alt,\%notes) : {}; #?? here
      #   }
      # }
      
      next;
    }
    
    ($qcoffb,$qcoffe)=split "-",$qcoff;
    ($tcoffb,$tcoffe)=split "-",$tcoff;
 
    my($qcaln,$qualn,$tcaln,$tualn)=(0) x 4; # split align to cds/utr portions, for each seq
    if($qae < 1) { #  or $qcoffe < 1 ## fail??
      $qcaln= _minb($aln, $qc);  $qualn= _minb($aln, $qw - $qc);
    } else {
      $qcaln= _maxb(0, _minb($qae,$qcoffe) - _maxb($qab,$qcoffb));
      $qualn= _maxb(0, $qae - $qcoffe) + _maxb(0, $qcoffb - $qab);
    }
    if($tae < 1) { #  or $tcoffe < 1
      $tcaln= _minb($aln, $tc);  $tualn= _minb($aln, $tw - $tc);
    } else {
      $tcaln= _maxb(0, _minb($tae,$tcoffe) - _maxb($tab,$tcoffb));
      $tualn= _maxb(0, $tae - $tcoffe) + _maxb(0, $tcoffb - $tab);
    }

    ## Q presumed larger than T 
    ## do AFTER 1st tc > qc swap and new tcaln calcs..
    ## also check if class{$td}/class{$qd} ??
    ## problems for aaref best but not nearlysame, for cdsover cases. want to swap anyway
    
    my $nearlysame=(abs($tc-$qc) < 30 and $tcaln >= 0.90*$qc)?1:0;  #?? tcaln >= 0.50*qc   
    my($thoval)= split",", ($aablast{$td} || "0");
    my($qhoval)= split",", ($aablast{$qd} || "0");
    $nearlysame=1 if($thoval > 5+$qhoval); # force swap check .. small diff?
    
    if($nearlysame) { 
      my $doswap=0;
      my($tispoor,$tqual,$tbits)= trevidence($td); # ,$tcla =~ m/althi|parthi/
      my($qispoor,$qqual,$qbits)= trevidence($qd); #,$qcla =~ m/althi|parthi/
      my($tho,$qho)=map{ my($ho)= (m/aaref:(\d+)/)?$1:0; $ho; } ($tbits,$qbits);
      my $hodiff=($tispoor and not $qispoor)?15:5;
      $hodiff += 10 if($isaltof and $qd=~/t1$/); # avoid swap out t1main
      if($tho > $hodiff+$qho) { $doswap=1; } # any quals but hoval to check?
      elsif($qho > 5+$tho) { $doswap=0; }
      # elsif($tispoor < $qispoor) { $doswap=1; } 
      elsif( ($qispoor & (kAATINY+kAAUTRBAD)) > ($tispoor & (kAATINY+kAAUTRBAD)) ) { $doswap=1; } 
      else { $doswap=0; }
      if($doswap) {
        ($qd,$td)=($td,$qd); 
        ($qc,$qw,$tc,$tw)=($tc,$tw,$qc,$qw); 
        ($qcoff,$qcode,$tcoff,$tcode)=($tcoff,$tcode,$qcoff,$qcode);
        ($qab,$qae,$tab,$tae)= ($tab,$tae,$qab,$qae);
        ($qgd,$tgd)=($tgd,$qgd); 
        ($qcoffb,$qcoffe,$tcoffb,$tcoffe)=($tcoffb,$tcoffe,$qcoffb,$qcoffe);
        ($qcaln,$qualn,$tcaln,$tualn)= ($tcaln,$tualn,$qcaln,$qualn);
      }
    }
    
    #?? aln as only-cds-overlap? need another val: qcaln x tcaln
    # overloc classes from q,t CDS v UTR aligns
    #   .. trfrag overlaps CDS of trgood > drop/relcass as partof or altof trgood.
    #     .. unless trfrag.UTR only over trgood.CDS and trfrag.CDS has good code
    #   .. trfrag overlaps only UTR of trgood, depends on trfrag coding quals: other locus or partof trgood
    #     .. trfrag.CDS overlaps trgood.UTR, trfrag > partof:trgood unless has strong code quals
    #     .. trfrag.UTR only over trgood.UTR, trfrag > other locus (code/noncode?)
    #   .. trfrag overlaps trfrag, depends on CDS v UTR, code quals

    $pid= int(0.5+$pidn); # replaces: ($aln<1)?0: int(0.5+ 100*$iden/$aln); $pid=100 if($pid>100);
    my $alnfull=$aln;
    
    ## use "part" not "frag" qualifier?
    ## should this use only cdslen? $ispart=($tc < 0.50 * $qc)
    # my $ispart= ($tw < $ALTFRAG*$qw)?"part":""; # see below; altfrag == 0.5
    my $ispart= ($tc < $ALTFRAG*$qc)?"part":""; # see below; altfrag == 0.5
    
    #? $samesize=($tw == $qw and $tc == $qc)?2:($tc == $qc)?1:0;    
    # my($tispoor,$tqual,$tbits)= trevidence($td,$tcla =~ m/althi|parthi/);
    # my($qispoor,$qqual,$qbits)= trevidence($qd,$qcla =~ m/althi|parthi/);
    
    #============= overlocus classing ====================
    my $fclass="";
    # my $alnsig = ($alnfull >= 0.10 * _min($qw,$tw))?1:0;
    if($isaltof) { # ($qgd eq $tgd) 
      $fclass=($ismainof)?"mainof":"altof";  # handle these other way, but separate if NO/minor overlap? == paralogs??

    # } elsif($qcaln >= $TINYALNBASES) { # == 90b ** ADD pctalnmax == qc/tc; $TINYALN == 20 now, too big?
    } elsif($qcaln >= $TINYALNBASES and ($alnfull >= 0.10 * $tw)) { # == 90b ** ADD pctalnmax == qc/tc 
      if($tcaln < $NHIALN) { # == 9b
        $fclass="utrovcds$ispart"; # this must be >= TINYALNBASES ?
      }
      else { $fclass="cdsovcds$ispart"; } # == frag
      
    } elsif($qualn >= $TINYALNBASES and ($alnfull >= 0.10 * $tw)) {
      if($tcaln < $NHIALN) { $fclass="utrovutr$ispart";  } # this must be >= TINYALNBASES ?
      else { $fclass="cdsovutr$ispart"; } # == frag
    }
    unless($fclass) {
      if( $alnfull < $TINYALNBASES or ($alnfull < 0.10 * _min($qw,$tw)) ) {
        $fclass = "notover1"; 
      }
    }
    
    ## classifyFullTr() parts .. do AFTER find cross-locus alignments
    ## oneloc == cdsovcds or (utrovcds and tcdsisbad)
    ## twoloc == (utrovutr or cdsovutr) and not tcdsisbad
    ## ambig  == others?
    #   my $fclact="";
    #   if($fclass =~ /cdsovcds/ or ($fclass =~ /utrovcds/ and $ispoor)) { 
    #     if($fclass=~/part/ or $ispoor) { $fclact="drop.$fclass"; } # frag/part
    #     else { $fclact="alt.$fclass"; }
    #   } elsif($fclass =~ /utrovutr|cdsovutr/ and not $ispoor) {
    #     $fclact="locus.$fclass"; # keep? use orig locus class ?
    #   } else {
    #     $fclact="maybelocus.$fclass";
    #   }
 
   
    #============= cdsalign classing ====================
    #above# my $alnfull=$aln;
    $aln= _maxb($qcaln,$tcaln); # ?? which min/max or qcaln?
    
    my $antiflag= ($iscdsalign and $bits < 0) ? "/-sense" : "/."; ## bestmatch="id,99/89/-sense" for antisense?
    $samesize=($tc == $qc)?1:0; # OR tinyalndiff? abs($tc - $qc) < $NHIALN
    
    ## eqgene classifier, not for overlocus ?? or yes?
		
		$havepair{$qd}{$td}++; #?? change this to pal align val ?
		$lastd=$qd;

    ## ABOVE now; 
    # if($td eq "self" or $qd eq $td) { # putspan fix for no-align cases
    #   $td=$qd; $bestmatch{$qd}="$td,100/100" unless($bestmatch{$qd});
    #   next;   
    # }

    my $qcds= ($qw>0) ? 100*$qc/$qw : $OK_CDSUTR;
    my $tcds= ($tw>0) ? 100*$tc/$tw : $OK_CDSUTR;
    my $qutrbad= ($qcds >= $OK_CDSUTR)?0:1;  
    my $tutrbad= ($tcds >= $OK_CDSUTR)?0:1;  

    my($qsize,$tsize)= ($iscdsalign) ? ($qc,$tc) : ($qw,$tw);
    $alnmax= ($qsize>$tsize and $tsize>0)?$tsize:($qsize>0)?$qsize:$tsize;  

    #deflt: $ALTFRAG= 0.5; .. change this for overlocus?
    $isfrag= ($tsize < $ALTFRAG*$qsize)?"frag":"";
    
    $pal= ($alnmax<1)?0 : int(0.5+ 100*$aln/$alnmax); $pal=100 if($pal>100);
    # my $palq= ($qsize<1)?0 : int(0.5+ 100*$aln/$qsize); $palq=100 if($palq>100);
    # my $palt= ($tsize<1)?0 : int(0.5+ 100*$aln/$tsize); $palt=100 if($palt>100);
    
    my $skiptinyaln=0;
    my $tisaltpar=0; # treat like $skiptinyaln for now
    my $tinyalndiff= ((($alnmax - $aln) < $NHIALN)) ? 1 : 0; #yes: TEST3 && 

		$havepair{$qd}{$td}= $pal; #?? change this to pal align val for altpar ?

    # for both isaltof/notaltof, add antiflag to pidval.
    my( $eqpal,$eqpalmap,$eqantiflag,$eqflagg ) 
        = ($neqgene>0)? overlocus_eqgene($qd,$td,$pal) : (0,0,0,0);

    if($isaltof) { # $neqgene>0 and 
      $antiflag .= $eqantiflag if($eqpalmap>0);  # ismainof?? /$qd ?? wrong got self:
      # Aedesg12bEVm000004t3	drop	althicull/altof	Aedesg12bEVm000004t1	100/92/./altmap98xeq/Aedesg12bEVm000004t3
    } else {
      # my( $palg,$palmap,$antiflagg,$eqflagg ) = overlocus_eqgene($qd,$td,$pal) ;
      if($eqpalmap>0) { #always?
        $pal= $eqpal; $antiflag .= $eqantiflag;
      } else {
        # cancel == paralog/alt
        $tisaltpar=1;  $antiflag .= "/noteqgene"; # $skiptinyaln= 1;
        $fclass =~ s/cdsovcds/utrovutr/; #??
        # bug: Aedesg12bEVm002491t1	drop	noclass/drop.cdsovcdspart	Aedesg12bEVm000291t1	99/77/./noteqgene
      }
    }
    
    my $pidalnval="$pid/$pal$antiflag"; # old: $pid/$pal
     ## .. change these to use palq, palt : palign per q,t size ?
    $bestmatch{$qd}="$td,$pidalnval" unless($bestmatch{$qd});
    $bestmatch{$td}="$qd,$pidalnval" unless($bestmatch{$td});  
        
    my $qclass= $class{$qd}||"";
    if($samesize and $qclass =~ /alt/ and $ismain{$td}) { # check/stop alt1 <> alt2 assigns
      next if($bestmatch{$qd} =~ /$td,/); # problem here? for many equal bestmatch
    }

    if(UseTINYALNBASES) { ## == 1
      if($alnmax >= $MINCDS) { # recall alnmax is min(qsize,tsize); never skip tiny prots, MINCDS =~ 90 bases ?
        my $minbaseover= ($alnmax >= $TINYALNCDSLEN)? $TINYALNBASES : $alnmax * 0.10; ##$MIN_PCTOVERLAP;
        $skiptinyaln=($aln < $minbaseover)?1:0;
      } else {
        $isfrag= "frag" unless($isfrag); #or "tiny" ? # dont skip assign alt/frag to tiny cds, but maybe always set isfrag ?
      }
    } else { # old, pct-TINYALN is a problem for large cds  
      $skiptinyaln= ($pal < $TINYALN)?1:0;
    }
  
    # ** change for fclass, eg. isfrag not relevant for separate loci, no cds overlap..
    unless($fclass =~ /cdsovcds|altof/) { 
      $isfrag=""; 
      if($fclass =~ /utrovcds|utrovutr|cdsovutr/) {
        $skiptinyaln=1; # yes, means dont reclassify for no/tiny cds-overlap, only for cdsovcds
        }
      }

    # set alt aclass here, or not.
    if( $skiptinyaln or $tisaltpar) { } ## defer: $aclass="noalign$isfrag";  
    elsif( $tc == 0 and $qc > 0 and $pid >= $PLOW ) { $aclass="frag0aa"; } # tc==0 special case "frag0"
    elsif( $pid >= $PHI ) { 
      $aclass= ($isfrag)?"parthi":"althi"; # Primary classifier; change parthi to althipart ? althi$isfrag ?
      $aclass .= "1" if($pal >= $PHIALN or $tinyalndiff); # althi1 == hi-align + hi-ident, nearly same (identicals gone here)
      }
    elsif( $pid >= $PMID ) { $aclass="altmid$isfrag"; } # partmid for isfrag ?
    elsif( $pid >= $PLOW ) { $aclass="altlo$isfrag"; }  # partlo for isfrag ?

    if(doALTCULL and $isaltof and not $ismainof and $aclass =~ /althi|parthi/) { 
      ## DO this after below cds-over aclass = althi 
      ## add cullExonEq() test? but need input with gmap feq:ID/altmapxe87.100 .. 
      ##  $eqpal,$eqpalmap,$eqantiflag,$eqflagg 
      if($eqflagg) { # ($neqgene>0)
        # my( $palg,$palmap,$antiflagg,$eqflagg ) = overlocus_eqgene($qd,$td,$pal) ; # this sets feq altmapxe..
        if($eqflagg =~ m/altmapxe\d/) {  
          ## need aaref:val,refbest/good/ok from trevidence,  $note=~/feq:([^;:\s]+)/,  
          ## add trevidence(qd) to get chrmap mapqual,loc for compare .. here? or cull1ExonEq
          my $qclass= $class{$qd}||""; 
          ## now in cull1x: my($ispoor,$tqual,$tbits)= trevidence($td,1); # $aclass =~ m/althi|parthi/
          my($culls)= cull1ExonEq($qd,$qclass,$td,$aclass,$eqflagg); # "$tbits,feq:$qd/$eqflagg" 
          if($culls) { $aclass= $culls; } # done above:  $antiflag .= $eqantiflag; 
        }
      }

    }

    #yes: if(TEST3)
    if(my $ocl= $class{$td} and not $aclass =~ /cull/) {  
      #yes: if(TEST1602)   
      if($ocl =~ m/^alt|^part/) { $class{$td}= $aclass if($aclass eq "althi1"); next; } # test3
      elsif($ocl =~ m/^main/) { 
        #* problem of 2+mains w/ several essentially identicals,  break tie here **
        if($qclass =~ m/^main/ or $tisaltpar) { 
          # skip to below aclass set, test5, yes, do this also?
          } 
        else { next; } # * YES, need this : test4
      } 
    }

    # add/rep fclass to aclass
    if($aclass and $fclass) { $aclass.="/$fclass"; }
    elsif($fclass) { $aclass="noclass/$fclass"; } #?? main/ noclass/ ; change part/ to other
    ## m4class most okay set
    # 10941 okay	althi/alt   << these are presumed prior mistakes, cdsover for diff loci, qualify w/ last-oid equiv?
    # 13132 okay	main/notover
    # 44160 okay	noclass/notover
    # 23274 okay	part/locus  << noclass/locus ?

    if($aclass) { 
      my $qmain=$qd; my $attr= $pidalnval; ## "$pid/$pal$antiflag";

      if($class{$qd}) {  #  =~ /alt|part/
        my $qnext= $qd; 
        my $more= ($bestmatch{$qnext})?1:0;
        my %qdid=( $qd => 1, $td => 1);
        while($more) { 
          $more=0;
          my($qbest,$qpal)=split",",$bestmatch{$qnext};
          
          #yes: if(TEST1603)
          $qbest=0 if($ismain{$qnext} and not $ismain{$qbest}); # want always? maybe bad?
          
          if($qbest and not $qdid{$qbest}) {
            $qnext= $qbest; $qdid{$qbest}=1; 
            $more=1 if($class{$qnext} and $bestmatch{$qnext});
          }
        }
        if($qnext ne $qd) { $qmain= $qnext; $attr="$pidalnval/$qmain"; } #was $attr.="/$qmain";
        ## WARNING: now attr has /-sense or /. before /qmain;  evgmrna2tsa2.pl needs to know...
      }
      
      $ismain{$qmain}++;   
      $class{$td}=$aclass; #? problem here when class{td} == main set before, should not reset to alt
      #yes: if(TEST1602)  # this way works.  if($bestmatch{$td}) must be true from above
      if($bestmatch{$td} =~ /^$qd/) { $bestmatch{$td}="$qd,$attr"; }
      elsif(not $bestmatch{$td}) { $bestmatch{$td}="$qd,$attr"; }

      #yes: if(TEST3) #add to prevent 2+ circular alta <> altb <> altc with no main ..
      $class{$qmain}="mainlocus" unless($class{$qmain}); # should always be:  $qc >= $tc//$samesize..

    }
      
  } close($inh);

  # fclass: cdsovcds,utrovcds,cdsovutr,utrovutr + part/isfrag
  ## problem here using q/pal when those are tiny cds of mostly utr overlap (for main/noclass at least)
  ## bestmatch is no longer best cds match (with alts), but best cross-locus match
  # END:  print trclass table; add more fields to output: aaqual, aablast, tr,aa sizes?
  { my($q,$pal,$c,$d);
  foreach $d (sort keys %class) {
    ($q,$pal)=split",",$bestmatch{$d}; $c= $class{$d}||"altclass"; # any missing?
    my @cla= classifyFullTr( $d, $c, $q, $pal);
    print $outh join("\t",@cla)."\n"; 
    }
  foreach $d (sort grep{ not($class{$_}) } keys %ismain) {
    ## FIXME: "main" here doesn't mean that, these are the better compared to overlapped loci-tr
    ## change to "oklocus" ? okovlocus?  "drop" classing from other aspects like ispoor Noncode val
    ($q,$pal)=split",",$bestmatch{$d}; $c= "mainlocus"; # was "main";    
    my @cla= classifyFullTr( $d, $c, $q, $pal);
    print $outh join("\t",@cla)."\n"; 
    }
  foreach $d (sort grep{ not($class{$_} or $ismain{$_}) } keys %bestmatch) {
    ($q,$pal)=split",",$bestmatch{$d}; $c= "noclass";  
    my @cla= classifyFullTr( $d, $c, $q, $pal);
    print $outh join("\t",@cla)."\n"; 
    }
  }  
  
}

# 1;

 
  
  
=item classifytr

  classifier of locus primary, alternate, fragment, redundant 
  using identityclass() collection of overlapping transcripts (CDS or full tr)  
  focused on CDS qualities, maybe should be classifyCDS() ..
  
=cut


sub classifytr {
  my($tid,$cla,$qid,$pidal)= @_;

  # use constant { kAATINY => 1, kAAUTRBAD => 2, kAAUTRPOOR => 4, kAADUP => 8, kAAGAPS => 16, };
  # use constant NOTPOORBAD => kAATINY + kAADUP + kAAGAPS; # - kAAUTRPOOR - kAAUTRBAD - kAAGAPS
  # use constant NOTPOOR => kAATINY + kAADUP + kAAUTRBAD + kAAGAPS; # - kAAUTRPOOR
  $qid||="0"; $pidal||="0";
  
  ## pidal == pident/palign,otherattr ==  "$pid/$pal$antiflag";
  my($pidi,$pali)= $pidal =~ m/^(\d+).(\d+)/;
  
  my $aw= $aasize{$tid} || 0;
  my $tqual= $aaqual{$tid} || ""; #? parse for "aasize,pcds,aaqual"
  ## add these utr class quals : this should be in tqual string: aasize, pcds, aaqual
  my $tw= $trsize{$tid} || 1; 
  my $pcds= int(300*$aw/$tw);  
  my $butr= ($tw - 3*$aw); # okutr if butr<=300p, modify BAD_CDSUTR for small aw
  
  my $ispoor= ($aw < $AAMIN or ($tqual =~ m/partial/ and $aw < $AAPART))?kAATINY:0;

  my $tbits= $aablast{$tid} || "0,0";
  my($tbscore,$tbref)=split",",$tbits;
  #FIXME: add tag to tbits output:  $tbits="aaref:$tbits" unless($tbits=~/^0,0/);
  
  #TEST1602 add on tbits, like aablast ref info? equalgene mapqual info 
  my $mapqloc= ($nineqgene) ? $gmapqual->{$tid} : 0;  
  my($mapqual,$maploc)= ($mapqloc)?(split"\t",$mapqloc):(0,0);
  # $gmapqual{$mid} == "$mapqual\t$loc"; 
  # new bleqgene mapqual: 85a,100i,2469l,4x,Spl:29%,KB668900 == %align,%ident,bases,exons,Spl:split

  # change to aablastref{tid} ?
  my($tbrscore,$tbrefbest)= ($tbref) ? (split",", $aablastref{$tbref}) : (0,0); 
  my $rbits2= $aablastref{$tid} || "0,0";
  my($tbrscore2,$tbref2)= split",", $rbits2; # this is revhash tid => score,ref refbest

  my $aaclus= $aacluster{$tid} || "0,0";  
  my($aamainid,$aamainpi)= split",", $aaclus; # ="$mainid,$pi";
  if($aamainid eq $tid) {  $aamainid=""; $aamainpi=0; }
  ## FIXME: aamain can be bad/missing; fix before this, need cds/tr align info to know if aamainid is good.
  
  
  ## main class should use aasize, aapartial and utrbad/poor with AAMINBAD, AAMINPOO
  ## FIXME: pcds bad here? comes from aw - gaps effects; separate gapbad from utrbad
  ## 201402: maybe ignore tqual utrbad/poor and use only pcds with these cutoffs? dont know tqual cutoffs
  if($tqual =~ m/gapbad/) { $ispoor |= kAAGAPS; }  # gaps discrepancy w/ aaqual utrbad and pcds utrbad 
  if($butr >= $MINUTR) { # ~ 300b "fixed" average utr sizes
    # should adjust %cds bad for cds size, for large cds> 10k, pcds >=90%, for small cds < 300, pcds may be 33%
    if($tqual =~ m/utrbad/ and $noRESET_CDSUTR) { $ispoor |= kAAUTRBAD; } 
    elsif($tqual =~ m/utrpoor/ and $noRESET_CDSUTR) { $ispoor |= kAAUTRPOOR; }
    else { $ispoor |= ($pcds < $BAD_CDSUTR)? kAAUTRBAD: ($pcds < $OK_CDSUTR)? kAAUTRPOOR : 0; }
    # if($pcds < $BAD_CDSUTR) { $ispoor |= kAAUTRBAD; } elsif($pcds < $OK_CDSUTR) { $ispoor |= kAAUTRPOOR; } 
  }

use constant CODEPOT1607 => 1; # add from above codepot work
use constant ALTPOOR1607 => 1;
  my $tcode=0;
if(CODEPOT1607) {
  $tcode= $sizeval{$tid}{'codepot'}||""; # $tcode =~ /^Noncode/ set what? kAATINY+kAAUTRBAD
  $tcode=0 if($tbscore>0);
  $ispoor |= kNONCODE if($tcode and $tcode =~ /^Noncode/);  
}

  ## TEST1603: $pidal =~ m/altmap\d+xeq/ ; xeq for althi1 means redundant == AADUP
  # problem cases, xeq but not same due to poor align, diff aa size, etc. : want to keep both these, t2 is valid alt
  # Anofunz4hEVm000070t1	okay	althi	Anofunz4hEVm000070t2	100/52/./altmap87xeq	3437,61%,complete	aaref:4110,AGAP000222-PB,refgood,chrmap:100a,99i,10314l,8x,KB668689:297566-312085:-,pflag:0
  # Anofunz4hEVm000070t2	okay	main	Anofunz4hEVm000070t1	100/52/./altmap87xeq	3511,52%,complete-utrpoor	aaref:4064,AGAP000222-PB,refgood,chrmap:87a,99i,10536l,8x,KB668689:297566-312085:-,pflag:0
  
  if($aamainpi > $AADUP_IDENT) { # TEST aacluster added info on dup prots; AADUP_IDENT=98 default
    ## check that aamainid is main ??
    $cla.="a2" unless($cla=~/hi[012]/); # hi1/a2 redundant info
    #defer# $tbits.=",aadup:$aamainid";
    $ispoor |= kAADUP if($cla =~ /althi/); # only for /althi|parthi/
  }
  unless( $ispoor & kAADUP ) {
    if($pidal =~ m/altmap\d+xeq/ and $cla =~ /althi1|part/) {
      $ispoor |= kAADUP;  # add new bitflag kDUPEXONS ?
    }
    # V4: eqgene altmap100 fixup, these have identical CDS set..
    # if($pidal =~ m/altmap(\d+)/) { my $am=$1; 
    #   if($am >= $AADUP_IDENT and $cla =~ /althi1|part/) {
    #     $cla.="a2" unless($cla=~/hi[012]/); # hi1/a2 redundant info
    #     #? $tbits.=",aadup:$aamainid";
    #     $ispoor |= kAADUP;  # add new bitflag kDUPEXONS ?
    #   }
    # }
    # if(($maploc =~ /\.icdup/) and ($cla =~ /alt|part/)) { $ispoor |= kAADUP; }  #? kDUPEXONS
  }
  my $aadupflag= ($ispoor & kAADUP and $aamainid) ? ",aadup:$aamainid" : ""; # defer, AFTER refbest/good/..
  # below add to  $tbits.= $aadupflag

  # TEST3 : add? ispoor for cla == althi1 and pid/pal == 100/100  ie identicals?? NO, pal is align/shortlen
  # NO, not this, althi1 can be good alt, need align/mainlen stat also. or use only AADUP score as now.
  
  # CHECKME: adding aablast kept 40k more okay, all althi/ahia2 + 2k parthi
  # .. are these true uniq aablast or just althi aadups ?
  my $MINBLASTSCORE= 60; # aablast only? bitscore always?
  my $keepdrop="";
  unless( $tbscore == 0 or $tbits=~/^0,0/) {
    my $risbest= (($tbrefbest eq $tid) or ($tbscore >= $tbrscore))?1:0;# was ($tbrefbest eq $tid), allow 2+ same score
    my $risgood= ($tbrscore > $MINBLASTSCORE and $tbscore >= 0.90 * $tbrscore)?1:0;
    my $risok= ($tbref2 and $tbrscore2 >= $MINBLASTSCORE)?1:0;
    if($risok and $risgood) { $risgood=0 if($tbrscore2 > $tbrscore); } #?
    $ispoor = $ispoor & NOTTINY if($risgood); ## FIXME1712: clear kAATINY for $risok/risgood
    ## ispoor & kAADUP should not be cleared here???
    my $isaadup=($ispoor & kAADUP);
    if($risbest) { 
      $tbits.=",refbest"; $keepdrop.="okay:refbest,"; $ispoor=0 unless($isaadup); }
    elsif($risgood) { # was 3rd
      $tbits.=",refgood"; unless($isaadup){ $keepdrop.="okay:refgood,"; $ispoor=0;} } # maybe2 ok ?
    elsif($risok) {  # was 2nd ? why superceed refgood? 2nd ref match?
      $tbits.=",refok"; unless($isaadup){ $keepdrop.="okay:refok,"; $ispoor=0;} }
  }
  
      # 201402: update for drop.mains, perfectfrag replacements, need input list of perfectfrag/dups from fastanr/cdhit prior steps
      # 201402: AAPART & AAMIN, AAMINBAD, AAMINPOO  need fix to prevent drop uniq/valid main/noclass genes, 
  ## maybe ignore utrbad for main .. otherwise can miss true longutr/ncutr genes; but refbest/good will keep those w/ homology
  ## keep largeaa,utrpoor for main, noclass, altlow; but drop smallaa,utrpoor
  if($ispoor > kAATINY and $cla !~ /althi|parthi/) {  
    if($aw >= $AAMINBAD) { $ispoor = $ispoor & NOTPOORBAD; } ##  ^ (kAAGAPS+kAAUTRPOOR+kAAUTRBAD);  ##?? ^ XOR bad
    elsif($aw >= $AAMINPOO) { $ispoor = $ispoor & NOTPOOR; } ##  ^ (0+kAAUTRPOOR);
  }

  if($cla =~ /parthi|frag0aa/) {
    $keepdrop.= "drop";
    
  } elsif($cla =~ /althi/) {
if(ALTPOOR1607) {
    $keepdrop.= (($ispoor & NOTPOORBAD) or ($ispoor and $pali > 75))?"drop":"okay";
} else {
    $keepdrop.= ($ispoor)?"drop":"okay";  # ispoor vs main size?
}

## 9jul16: alt class bug: drop.althi when <30% align to main, this should be valid alt due to small align
## problem ispoor from utrpoor, but is valid unique aahoref, possible paralog not alt
# cornhi8mtrinLocDN35783c0g3t2    drop    althi   cornhi8m9agvelvk35Loc6006t1     100/32/.        262,54%,complete-utrpoor        0,0,pflag:4
## test solution: ignore ispoor when $pal < 75% ?
    
  } elsif($cla =~ /main/) {
    # $keepdrop.= "okay"; # always ok if has alts? NO, 
    # Fixme: keeping all "main" class gives ever expanding trsets with added trasm
    # should apply aaqual, aaref criteria to main + its althi; 


## 21feb update: dmag5icd35all_cde35.class4
## class7 maindrops: utrpoor/bad:48144, then partial:19493, gapbad:2866,  397 are aalong+complete
## class6 maindrops: utrpoor/bad:57636, then partial:19815, gapbad:3178,  1075 are aalong+complete
# dmag4vel4ibnk31Loc11665t1       drop    main    dmag4vel4ifik31Loc13364t1       99/99   102,71%,complete        0,0,pflag:4
# dmag4vel4ibnk31Loc1189t4        drop    main    dmag4vel4ipak31Loc21101t1       100/35  364,62%,complete        0,0,pflag:18
## pflag:18 = kAAGAPS + kAAUTRBAD
#..
## class4 maindrops: utrpoor/bad:60742, then partial:19737, but 3000 are aalong+complete, WHY drop?
## class5 maindrops: utrpoor/bad:56971, then partial:18807, but 1788 are aalong+complete (have gaps, so utrbad/poor now)
# dmag4vel4ibnk31Loc10846t1       drop    main    dmag4vel4ipak31Loc13986t1       99/94   116,73%,complete        0,0,pflag:4
# dmag4vel4ibnk31Loc11089t1       drop    main    dmag4vel4ifik31Loc18962t1       99/68   129,70%,complete        0,0,pflag:4
# >dmag4vel4ibnk31Loc10846t1 is aacluster:main, not aatiny, not utrpoor, not kAADUP, not althi;
# .. its alt is poor: partial5-utrbad, 37aa; did AAcluster fix mangle main?  pflag:4 == kAAUTRPOOR, miscalc from $pcds ???
# ** GAPS in aa above; 116aa,73%pcds is wrong if -gaps removed. BUT gaps not removed from nt size, so pcds off by that.
#  dmag4vel4ibnk31Loc10846t1 68aa,48gap,475nt,pcds=3*68/475=43%;  pcds= 3*68/(475-3*48) = 62%
# .. so drop.main/aacomplete is ~correct calc given aagaps : probably dont want utrbad flag here, but a gapbad may be useful
# .. 48gaps/116aa is not good prot; add flag gapbad == kAAGAPS above in readAAqual ?


    # $ispoor = $ispoor & kAATINY; # clear bits 2,4,.. # now above
    $keepdrop.= ($ispoor)?"drop":"okay";  # ??
    
  } elsif($cla =~ /noclass/) {
    $keepdrop.= ($ispoor)?"drop":"okay"; #? dont drop? ditto to altmid
    
  } else { # other altmid/low 
    $keepdrop.= ($ispoor)?"drop":"okay"; #? dont drop this class, ispoor maybe uniq prot: maybeok?
  }
  
  my $okay;
  if($keepdrop =~ /drop/) {
    # okaymaybe/okmaybe instead of maybeok ?
    if($keepdrop =~ /okay/) { $okay= "maybeok"; } else { $okay= "drop"; }
  } else {
    $okay= "okay";
  }
  
  $tbits.= $aadupflag if($aadupflag); # defer, AFTER refbest/good/..
  $tbits.= ",chrmap:$mapqual,$maploc" if($maploc and $mapqual);
  $tbits.= ",pflag:$ispoor"; # DEBUG info, but keep
  #FIXME: add tag to tbits output, here at end?
    $tbits="aaref:$tbits" unless($tbits=~/^0,0/);
  #TEST1602 add on tbits, like aablast ref info? equalgene mapqual info 
  # $gmapqual{$mid} == "$mapqual\t$loc"; 
  
  if(defined $eqflag{$tid}) { # TEST1603
    my $eqfl= $eqflag{$tid}{$qid}||""; 
    # debug output all eqflag{tid}{otherid}
    if($eqfl) { $eqfl="$qid/$eqfl,"; }
    my @q= grep{ $_ ne $qid } sort keys %{$eqflag{$tid}}; 
    $eqfl .= join",",map{ "$_/".$eqflag{$tid}{$_} }@q;  
    $tbits.= ",feq:$eqfl";
  }
  
  return (wantarray) ? ($tid,$okay,$cla,$qid,$pidal,$tqual,$tbits) : $okay;
}


=item identityclass

  from  aabugs4qual/tsaevg/cdsidentclass.sh
  input from above putspans, sorted by ^Qclen, ^Align, vTtlen:
     qw(Qid Qclen Qtlen Tid Tclen Ttlen Align Ident Bits) 

  update classes per aabugs4qual/tsaevg/classparts.pl
  .. need blastp/blast table input for this..
  .. need 3-4 final categories:  keep-main, keep-alts, trash-fragments (alt,noalt) 
     .. keep includes partial prot where good enough.
     .. trash-alts maybe (not fragment but not distinct enough)
     .. add poor cut to main,alt,noclass per below a,b qualities.
     
cl_main.alttr and cl_alts, keep if: 
  a. homology unique or best (near best for alt),  or
  b. aasize >= 45% of main and aacomplete? and not aasame-as-main?
  c. main == fragment if aasize < MINAAFULL or (aasize < MINAAPART and partial/utrpoor)

FIXME/Check : aablast score has dropped/poor a set of 653 uniq ref (Nref1) matches vs prior clidtab/class1
  see env swapids=1 refblast=aaeval/refdebl/refde-.tall3 ./classparts.pl trsets/$pt.clid2tab
Partition               Nclass  Nmatch  Nrefa   Nref1   Bits    Iden    Algn    Rlen    Tlen
locust1best5.cl_poor    68122   1331    3485    653     113     61      170     361     389

# # FIXME2: maybe ignore utrbad for main .. otherwise can miss true longutr/ncutr genes; 
# # but refbest/good will keep those w/ homology
# cat $pt.class3 | grep 'drop.main' | sort -k6,6nr | head
# fungrvel4ik53Loc89t53   drop    main    fungrvel3bk43Loc142t24  99/66   1408,27%,complete-utrbad        187.2,DRERI:ENSDARG00000086955
#.. ^only aa equiv is also drop:
#  fungrvel3bk43Loc142t24  drop    althi   fungrvel4ik53Loc89t53   99/66   1373,32%,complete-utrbad        181.4,DRERI:ENSDARG00000086955
#
# fungrvel4k25Loc3571t1   drop    main    fungrvel4k29Loc22240t1  99/100  668,28%,complete-utrbad 0,0
# fungrvel4k25Loc15515t1  drop    main    fungrvel2k35Loc34556t1  100/100 580,24%,complete-utrbad 261,DRERI:ENSDARG00000068192
# .. other DR92 genes:
#  fungrvel3bk35Loc16459t1 okay    main    fungrvel3bk29Loc13518t1 99/92   876,82%,complete        543,DRERI:ENSDARG00000068192,refbest
#  fungrvel3bk29Loc13518t1 okay    althi   fungrvel3bk35Loc16459t1 99/92   836,61%,complete        572,DRERI:ENSDARG00000068192,refbest
#  fungrvel4k25Loc15515t1  drop    main    fungrvel2k35Loc34556t1  100/100 580,24%,complete-utrbad 261,DRERI:ENSDARG00000068192
#  fungrvel4k35Loc25356t1  okay    althi   fungrvel4k25Loc15515t1  100/65  363,99%,partial3        256.3,DRERI:ENSDARG00000068192
#
# fungrvel3bk29Loc5288t1  drop    main    fungrvel4k25Loc5314t2   99/87   565,23%,partial5-utrbad 0,0
# fungrvel4k25Loc3268t8   drop    main    fungrvel4ik53Loc14576t1 100/58  526,24%,complete-utrbad 396,DRERI:ENSDARG00000038737
#...    
    
=cut

sub identityclass {
  my($outh, $infile, $insorted)= @_;

use constant TEST3 => 1; # 13aug09 test fixes to alt/main classing
	my $havevalidids= scalar(%validids)?1:0; ## need this?

use constant TEST1602 => 1; ## UPD 2016.02 problem w/ dup equal mains, equivalence after 1st see both as mains, 1 to be alt,

## 2016.feb : problems now classing  alts => main/noclass, tho blastab infile has proper hi identity scores
##   found w/ genom-map eqgene, some essentially identical tr classed as sep loci ..
##   repeated runs of this on new oksets reduces this alt>loc misclass, but not entirely gone.

## 2013.aug : IS_CDSALIGN ORIENT in alntab, add -sign/antisense to one field
### -aln == reverse align ?? for CDS antisense problem
### -bits == reverse align ?? bits not used here for scoring.. best choice other than adding column
### -aln affects here sort: -k7,7nr; use other field Ident?

#... see above now
#   $TINYALN=$ENV{mina}||50; 
#   $IS_CDSALIGN=$ENV{cdsw}||0; # default using tralign, tests better than cds
#   $PHIALN=$ENV{ahi}||98; # was 65; # use mina?
#   $PHI=$ENV{phi}||99; $PMID=$ENV{pmid}||90; $PLOW=$ENV{plow}||80; 
#   # $ALTFRAG: add isfrag pct; now 50 (0.5)
  
### infile is putspans() alntab: Qid Qclen Qtlen Tid Tclen Ttlen Align Ident Bits
###  sort ^Qclen, ^Align, vTtlen: cat $atab | grep -v '^Qid' | sort -k2,2nr -k7,7nr -k6,6n 
### should it be ^Qclen, ^QID, ^Align, vTtlen : to keep qids together? no want top Aligns first
### * add vQtlen before vTtlen, else get utrbad before good ***
  # my $ALNSORTORD='-k2,2nr -k7,7nr -k6,6n '; #orig
  # my $ALNSORTORD='-k2,2nr -k7,7nr -k6,6n -k1,1 -k4,4'; #add IDs to order ties

  my $ALNSORTORD='-k2,2nr -k7,7nr -k3,3n -k6,6n -k1,1 -k4,4';  #add vQtlen/3 before vTt/6 IDs to order ties
  unless($IS_CDSALIGN) { # sort tlen, not clen 
    # ..not sure what is right, want cds to play role in best choice
    # .. input blast align k7 is for tlen, not clen, but want to choose best by long clen > long tlen,talign
    #t1: $ALNSORTORD='-k3,3nr -k7,7nr -k2,2nr -k6,6nr -k1,1 -k4,4';  # ^Qtlen,^Align,^Qclen,vTtlen,Qid,Tid
    $ALNSORTORD='-k2,2nr -k7,7nr -k3,3n -k6,6n -k1,1 -k4,4';  # test same as IS_CDSALIGN  
  }
  
  my($inh, %class,%bestmatch,%ismain, %havepair); 
  if(ref($infile)) { $inh= $infile; } # STDIN,.. sort ??
  else {
    ## FIXME fail /tmp sort use, what bug? using TMPDIR.. use -T ./
    my $tmpdir=$ENV{TMPDIR};
    $tmpdir= './' unless($tmpdir and $tmpdir ne "/tmp"); ## { $tmpdir=`pwd`; chomp($tmpdir); }
    if($insorted) { open(IN,$infile) or die "reading $infile"; $inh=*IN; }
    else { open(IN,"sort -T $tmpdir $ALNSORTORD $infile |") or die "reading $infile"; $inh=*IN; }
    #orig#else { open(IN,"sort $ALNSORTORD $infile |") or die "reading $infile"; $inh=*IN; }
    
    ## UPD 1602: maybe revise to drop sort, expect per-qd, top-hit sort order of blastn output?
    ## then change class assignments as new qd x td warrants?
  }
  
  my($lastd)=("");
  while(<$inh>) { # maybe local table, sorted, or from file
    next if(/^Qid|^\W/); chomp; 
    my @v= split"\t"; 
    my($qd,$qc,$qw,$td,$tc,$tw,$aln,$iden,$bits)= @v; 
    my($isfrag,$aclass,$alnmax,$pid,$pal,$samesize)= (0) x 10;
    
		## 2013.aug : IS_CDSALIGN ORIENT in alntab, add -sign/antisense to one field: not aln, iden?, bits?
    ## FIXME: *** 98/100/-sense/PitaEaR000975t24 << buggers, mixes up users of this table.
    ## all piad entries should have pd and sense/asense fields, placeholder where needed '.' ?
    
    #old# my $antiflag= ($IS_CDSALIGN and $bits < 0) ? "/-sense" : ""; ## append to bestmatch ??
    my $antiflag= ($IS_CDSALIGN and $bits < 0) ? "/-sense" : "/."; ## append to bestmatch ??
    		## bestmatch="id,99/89/-sense" for antisense?
    		## ^^ is this new flag causing problems?  Besthit is appended /after, parsed where?
    
    $isfrag= $aclass="";
    $samesize=($tc == $qc)?1:0; # OR tinyalndiff? abs($tc - $qc) < $NHIALN
    if($tc > $qc){ ($qd,$td)=($td,$qd); ($qc,$qw,$tc,$tw)=($tc,$tw,$qc,$qw); } # swap to greater cds?

    if(EQGENE_CHANGES_NOALN && TEST1602) {
      ## neqalts data not yet clean enough for this reclassing
      ## ? just ignore cds-align row for those marked as non-alts? 
      ## set any vals? class{},ismain{},havepair{} ?>
      #?? next if(ref $neqalts and $neqalts->{$qd}{$td}	);  

      # if(TEST1602 && $lastd && $lastd ne $qd) {	# reverse, neqalts	OR below/above: skiptonext if($neqalts and $neqalts->{$qd}{$td}	);	
      # if(ref $neqalts) { 
      #   my @ted= sort{ $aasize{$b} <=> $aasize{$a} or $a cmp $b } keys %{$neqalts->{$lastd}}; 
      #   foreach my $te (@ted) { 
      #     if($havepair{$lastd}{$te}) { # break pairing alt class
      #       my $laclass="mainmap"; #? mainparalog?
      #       if( $class{$te} =~ /^alt/) {  $class{$te}= $laclass; }
      #       #?? elsif( $class{$lastd} =~ /^alt/) {  $class{$lastd}= $laclass; }
      #       $havepair{$lastd}{$te}=0; 
      #       $havepair{$te}{$lastd}=0;
      #     } 
      #   }
      # }
      # }				  
    }
    

=item eqgene classifier

  - do this in readEqualGene, mark which overlap alts are bad, which ok
  - also need to regard mapqual align, Split values to decide
  
  need separate classifier to handle various eqgene attributes, decide which tr/alts are bad/good
  eg eqgene classifier for this case: good g131, bad g453t5
  .. for this case need to know that g453t5 <mismap> g453t1,2,3,4 .. count overlap/alt/locus? drop outliers? use Split info?
ok Anofunz4gEVm000131t1	noid	Anofunz4gEVm000452t5/19	KB668936:299087-306372:-	99a,99i,7287l,3x	0
ok Anofunz4gEVm000131t2	noid	Anofunz4gEVm000452t5/33	KB668936:299737-303981:-	98a,99i,4269l,2x	0
ok Anofunz4gEVm000131t3	noid	Anofunz4gEVm000452t5/58	KB668936:299632-301584:-	70a,99i,1809l,2x	0
ok Anofunz4gEVm000452t1	noid	na	KB668900:3696-9724:-	99a,100i,4410l,8x	0
ok Anofunz4gEVm000452t2	noid	na	KB668900:3696-9724:-	99a,100i,4410l,8x	0
ok Anofunz4gEVm000452t3	noid	na	KB668900:3696-9724:-	96a,100i,4329l,8x	0
ok Anofunz4gEVm000452t4	noid	na	KB668900:3696-8875:-	100a,100i,3372l,9x	0
bad Anofunz4gEVm000452t5	noid	Anofunz4gEVm000131t1/56,Anofunz4gEVm000131t2/56,Anofunz4gEVm000131t3/43	KB668936:300530-301921:-	85a,100i,2469l,4x,Spl:29%,KB668900	0
ok Anofunz4gEVm000452t6	noid	na	KB668900:6592-8636:-	79a,100i,948l,3x	0

  g131,3/3 alts are over 1/6 g453 alts
  -- should this indicate g453t5 is not-much-over other g453t alts?
  
=cut
		  
		if(EQGENE_CHANGES_NOALN && $lastd && $lastd ne $qd) { 
		  ## revise lastd class if warranted, when have new qd
		  #* change this altmap class: NO cds-align b/n td,qd means something, poor mapping paralogs end up here
		  # .. need (a) gmap qual score (align,ident,split), and
		  # ..      (b) paralog flag reverse of eqgene, for classed-alts that *dont* gmap same locus
		  # TEST1602: ($eqgenes,$nineqgene,$neqgene,$gmapqual,$neqalts)= readEqualGene($eqgene) extended table

			if($neqgene>0 && $eqgenes->{$lastd}) {
				my @ted= sort{ $aasize{$b} <=> $aasize{$a} or $a cmp $b } keys %{$eqgenes->{$lastd}}; 
      
        ## FIXMEd remove dupl map ids _G2,3.., from eqgene table: Funhe2Exy3m149508t1_G2 ..
        ## FIXME2: this turns all main class into altmap, need to keep 1 as main or otherwise flag main/altmap
        ## check each lastd/@ted for ismain{id}, class{id}, and aasize{id} : pick one to keep as main
			  
			  ## FIXME3: Still bad, 160217, replaces combined locus groups' main w/ althi1/altmap w/o reciprocal alt to main
			  
			  ##?? decide here w/ mapquals to reclass? or should equalgene table be revised for lowqual cases?
			  ## mapqual now == gmapalign%,nnn,location
			  # my $ldmapqual= (ref $gmapqual)? $gmapqual->{$lastd} : ""; # TEST1602
			  # skip_reclass if($ldmapqual < $MINGMAPQUAL);
			  
				my $miss=0;
				foreach my $te (@ted) { 
          unless($havepair{$lastd}{$te}) { 
            next if($te =~ /_G\d+$/); # dup map id syntax, dropped from eqgene table now?
            next if($havevalidids and not $validids{$te}); #?? fix for _Gnnn, or not? 
            ## should do in readEqualGene but dont have %validids then

 			      # my $temapqual= (ref $gmapqual) ? $gmapqual->{$te} : ""; # TEST1602
 			      # next if($temapqual < $MINGMAPQUAL);
 			      
            my $palmap= $eqgenes->{$lastd}{$te}||0;  # val is te align to lastd, not reverse
            my $palmapt= $eqgenes->{$te}{$lastd}||0;   
            # swap ids for largest pal ??
            if($palmap > 0) {
              $miss++; 
              #bad?# my $pidalnval="99/$palmap";  # .$fakeantiflag ??
              my $pidalnval="$palmap";  
             
              ## this is wrong way, te here is (often) main gene, lastd is fragment overlapping 
              ## maybe NOT bestmatch{lastd}, val is te align to lastd
              ## maybe replace bestmatch{te} if palmap > bestmatch{te} ?

              $bestmatch{$lastd}="$te,$pidalnval" unless($bestmatch{$lastd});

              #BAD??# $bestmatch{$te}="$lastd,$pidalnval" unless($bestmatch{$te});
              
              my $laclass="altmap"; # this may be replacing main wrongly, need mainmap? test class{te} first?
              #bad# $class{$te}= $laclass; #?? is this going to work
              $class{$lastd}= $laclass; 
              $havepair{$lastd}{$te}++; $havepair{$te}{$lastd}++; #??
              }  
            } 
					}
				if($miss) { } #.. reclass $lastd ??
			} # eqgene
		}
		
		
		#x $havepair{$qd}{$td}++; #?? change this to pal align val ?
		
		$lastd=$qd; # here?
    if($td eq "self" or $qd eq $td) { # putspan fix for no-align cases
      $td=$qd; # $aclass="noclass"; << should be this, but ??
      $bestmatch{$qd}="$td,100/100/self1" unless($bestmatch{$qd});
      next;  # can lose eqgene if no further td/qd ..
    }

    #FIXME: have 2nd perfect matches, e.g qd1 >> td but qd2 == td; need that class to drop dups.
    # ** BETTER: Remove these before align/cluster; fastanrdb on .cds, .aa?
    #  asmrna5/cdsx/alldmag5x.cds n=4093166; alldmag5x.nrcds n=1782663; 501494 have dups (some many-many)
    # .. info is in blastn aligns, not in cdhit clusters (no 2ndary aligns).
    # .. use alntab ident column, if cdssize1 = cdssize1 = ident-3 (stop), identicals
    # eg dmag4vel4xbxk55Loc9866t1  dmag4vel4xfik55Loc11404t1 dmag4vel4xpak25Loc1083t9 : cds-ident alts of main dmag4vel4ibxk55Loc9119t1
    
    ## add here?? $aaqual{$qd} ; $aaqual{$td}; .. maybe do this below, END
    #  my($qqual,$tqual)= ($aaqual{$qd},$aaqual{$td});
    #  my($qbits,$tbits)= ($aablast{$qd},$aablast{$td});

    ## FIXME: check qw-main pCDS for UTRBAD/POOR; dont call tw frag if qc/qw is UTRBAD/POOR
    ## FIXME1405: sometimes? use input $aaqual{$qd}, ignore OK_CDSUTR/BAD_CDSUTR unless option..
    ## my $maw = $aasize{$qd} || 0;
    ## my $mqv = aaqualscore($aaqual{$qd});  #? parse for "aasize,pcds,aaqual"

    if(0) {
    ## NOT YET USED: add these utr class quals?
    my $qcds= ($qw>0) ? 100*$qc/$qw : $OK_CDSUTR;
    my $tcds= ($tw>0) ? 100*$tc/$tw : $OK_CDSUTR;
    my $qutrbad= ($qcds >= $OK_CDSUTR)?0:1; # ($qcds > $BAD_CDSUTR)?1: 2;
    my $tutrbad= ($tcds >= $OK_CDSUTR)?0:1; # ($tcds > $BAD_CDSUTR)?1: 2;
    }
    
# FIXME: when qc == tc, can assign althi1 = althi2 instead of main1 = althi2
    my($qsize,$tsize)= ($IS_CDSALIGN) ? ($qc,$tc) : ($qw,$tw);
    # note: alnmax is min(qsize,tsize) not max
    $alnmax= ($qsize>$tsize and $tsize>0)?$tsize:($qsize>0)?$qsize:$tsize;  
    $isfrag= ($tsize < $ALTFRAG*$qsize)?"frag":"";
    
    #o if($IS_CDSALIGN) { $alnmax=($qc>$tc and $tc>0)?$tc:($qc>0)?$qc:$tc;  $isfrag= ($tc < $ALTFRAG*$qc)?"frag":""; } 
    #o else { $alnmax=($qw>$tw and $tw>0)?$tw:($qw>0)?$qw:$tw;  $isfrag= ($qutrbad==0 and $tw < $ALTFRAG*$qw)?"frag":"";  }
    ## if($qc < $MINCDS) { $class{$qd}.="tiny"; } elsif($qc < $MINCDSPART and $ispartial<needqual){ ..}
    
    $pid= ($aln<1)?0: int(0.5+ 100*$iden/$aln); $pid=100 if($pid>100);
    $pal= ($alnmax<1)?0 : int(0.5+ 100*$aln/$alnmax); $pal=100 if($pal>100);
    #o my $palq= ($qsize<1)?0 : int(0.5+ 100*$aln/$qsize); $palq=100 if($palq>100);
    #o my $palt= ($tsize<1)?0 : int(0.5+ 100*$aln/$tsize); $palt=100 if($palt>100);
    #^ add $palq = $aln/$qc; $palt= $aln/$tc; or palmax = aln/alnmin 
    
    ##WrongWayPeachy#my $tinyalndiff= (TEST3 && (($aln - $alnmax) < $NHIALN)) ? 1 : 0; # TEST3 add, use w/ $pal
    my $tinyalndiff= (TEST3 && (($alnmax - $aln) < $NHIALN)) ? 1 : 0; # TEST3 add, use w/ $pal
    #NOT this,see samesize# my $samealn= ($tinyalndiff && abs($qc - $tc) < $NHIALN) ? 1 : 0; # TEST3 add2; move up to IS_CDS ??BAD calc?

## PROBLEM using eqgenes .. wont get here unless td x qd is in blastn.alntab ; some are NOT.
## check all $eqgenes->{$td} ?? should have self match.
##* Need to change this altmap class: NO cds-align b/n td,qd does mean something, poor mapping paralogs end up here
## 13Sep01 : add eqgene map info for missed alignments : pal adjustment
## 13sep06 : maybe not here .. trust align score if it exists for pair? at least needs flag if changes pal.

use constant TEST1603 => 1; ## test use here eqgenes/eqexons??

		$havepair{$qd}{$td}= $pal; #?? change this to pal align val for altpar TEST1603?

=item altpar problem case

-- may need align distance tree/cluster for altpars, and decision w/ other data on where to cut tree to separate loci
-- part of problem is at alt findmain, iteration should stop when palign drops a lot
   it t1/t4 is wrong main for many of these (t2,t3 and children)
   palign drops from >60% to <20% for t2,t3,..
   
-- there are 2 loci at least, diff gmap, diff aaref, alts of 2nd remain linked to 1st tho.
egrep '^Anofunz4hEVm004829t.    ok' evg24m2banofun.tgclass3     
>>locus1
>> t1 should be main not t4, aafull, refbest, 100align, 
>> is t4 artifact? may be mashup of loc1,loc2 as it aligns well to both sets
Anofunz4hEVm004829t1	okay	althi	Anofunz4hEVm004829t4	99/68/./altmap65xeq/Anofunz4hEVm004829t5	509,95%,complete	aaref:972,AGAP002866-PA,refbest,chrmap:100a,100i,1530l,2x,KB669169:1413398-1415008:-,pflag:0
Anofunz4hEVm004829t4	okay	main	Anofunz4hEVm004829t5	100/98/./altmap58	554,96%,partial3	aaref:897,AGAP002866-PA,refgood,chrmap:65a,98i,1662l,2x,KB669169:1413434-1414585:-,pflag:0
Anofunz4hEVm004829t5	okay	althi1	Anofunz4hEVm004829t4	100/98/./altmap58	490,84%,complete	aaref:885,AGAP002866-PA,refgood,chrmap:100a,100i,1430l,3x,KB669169:1413398-1415008:-,pflag:0

>>locus2
Anofunz4hEVm004829t2	okay	main	Anofunz4hEVm004829t4	99/3/./paralt15	509,89%,complete	aaref:928,AGAP002865-PA,refbest,chrmap:100a,100i,1530l,2x,KB669169:1375626-1377211:-,pflag:0
Anofunz4hEVm004829t8	okay	althi	Anofunz4hEVm004829t4	99/3/./paralt17	433,90%,complete	aaref:779,AGAP002865-PA,chrmap:100a,98i,1302l,2x,KB669169:1375626-1376983:-,pflag:0
>> ? locus2b w/ crossmap ?
Anofunz4hEVm004829t7	okay	althi1	Anofunz4hEVm004829t4	99/3/./paralt18	433,92%,complete	aaref:778,AGAP002865-PA,chrmap:100a,98i,1302l,2x,KB669169:1375626-1412196:-,pflag:0
Anofunz4hEVm004829t9	okay	althi1	Anofunz4hEVm004829t4	99/3/./paralt18	433,90%,complete	aaref:772,AGAP002865-PA,chrmap:100a,99i,1302l,2x,KB669169:1375626-1412196:-,pflag:0
Anofunz4hEVm004829t6	okay	altmid	Anofunz4hEVm004829t4	99/3/./paralt18	433,79%,complete	aaref:775,AGAP002865-PA,chrmap:100a,98i,1302l,2x,KB669169:1410838-1412196:-,pflag:0

>> ? locus3 : belongs w/ Anofunz4hEVm004760t23 Anofunz4hEVm005570t1 Anofunz4hEVm005571t1 and others n=29 same aaref, all over KB668962:33384-35764:-
>> but some have low align, may be 4th? paralog loc?
Anofunz4hEVm004829t3	okay	althi	Anofunz4hEVm004829t4	100/3/./paralt16	508,83%,complete	aaref:820,AGAP012957-PA,chrmap:84a,99i,1527l,6x,KB668962:33634-35764:-,pflag:0

=cut


		if($neqgene>0) { 
if(TEST1603) {
			my $palmap= $eqgenes->{$qd}{$td}||0; # note this is pct of td aligned to qd
      if($palmap>0) { 
			  my $xeq= $eqexons->{$qd}{$td}||0; # only for ($neqexon>0)
			  my $XPHI = 95; my $XPLO= 3;
			  # xeq care about (a) xeq>= identity == not alt but redundant, 
			  #   (b) xeq <= noxover, not alt but paralog maybe, (c) middle = usual alt
        
        #* set aln=0/TINYALNBASES when pal=0/TINYALN
        #*? change 'paralt' to 'altpar' to avoid other evg parse problems?
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
} else {			
			my $palmap= $eqgenes->{$qd}{$td}||0; 
			if($palmap>0) { $antiflag .="/altmap$palmap"; 
  			if(EQGENE_OVERRIDES_ALN && $palmap > $pal) { $pal=$palmap; } ##? not this change
			}
}			
		}

    my $pidalnval="$pid/$pal$antiflag"; # old: $pid/$pal
    my $tisaltpar=0; # treat like $skiptinyaln for now
    if(TEST1603) {
      $tisaltpar=1 if($pidalnval=~/altpar\d/); #?  and $qclass // not altparx\d
    }
    if(TEST1603) { # this may be sort of working right now..
    $bestmatch{$qd}="$td,$pidalnval" unless($bestmatch{$qd} and not($bestmatch{$qd}=~/altpar\d/));
    $bestmatch{$td}="$qd,$pidalnval" unless($bestmatch{$td} and not($bestmatch{$td}=~/altpar\d/)); 
     # prob here for altpar ?? need to replace td>qd
    } else {    
    ## .. change these to use palq, palt : palign per q,t size ?
    $bestmatch{$qd}="$td,$pidalnval" unless($bestmatch{$qd});
    $bestmatch{$td}="$qd,$pidalnval" unless($bestmatch{$td});  
        #^^ maybe should replace best{td}= new qd if new qd.pid >> oldqd.pid ?
    }
    #new# unless($bestmatch{$qd}) { $bestmatch{$qd}="$qd,100/100/self1"; } # always need bestmatch{qd} entry
        
# Maybe problem here making too many fragment alts: a,b near identical alts of c, b <= a subset,
# but both classed/bestmatch to c main.  b should be dropped as fragment/99equiv of a, but dont see that due
# to b<c precedence.

    unless(TEST3) {
    if($class{$td}) { next; }  # TEST3: yes, defer this to afer aclass, change for althi1 
    }    
    
    my $qclass= $class{$qd}||"";
    if($samesize and $qclass =~ /alt/ and $ismain{$td}) { # check/stop alt1 <> alt2 assigns
      next if($bestmatch{$qd} =~ /$td,/); # problem here? for many equal bestmatch
    }

    # if(TEST1603) {
    #   #x $bestmatch{$td}= "" if($tisaltpar); # prob here for altpar ?? need to replace td>qd
    #   #? maybe should be; causes problems : empty bestmatch
    #   if(0 and $tisaltpar) { # problematic; leave in, replace above ??
    #     $bestmatch{$td}= "" if($bestmatch{$td}=~/^$qd/); 
    #     $bestmatch{$qd}= "" if($bestmatch{$qd}=~/^$td/); 
    #   }
    # }
        
    #..........    
      ## reclassify? althi > althi100 for 2ndary alts; althi100 = unneeded duplicates : TEST
      ##if($class{$td} eq "althi" and $pid>99 and $pal>=99) { $class{$td}= "althi100b"; }
      # ^ this probably wont happen, 2ndary althi, 1st will also be 100 pi/pa
      # ** problem using cdhits clusters here, it doesnt give 2ndary aligns as blastn does
      # .. so partial align to big tr may also be perfect align to smaller, valid alt tr
      ## next;   # problem above for qd1,qd2 > td, need to keep qd2 with bestmatch{qd}
    
    # UPD 1602: change TINYALN use for MINALIGN bases
    my $skiptinyaln=0;
    if(UseTINYALNBASES) {
      if($alnmax >= $MINCDS) { # recall alnmax is min(qsize,tsize); never skip tiny prots, MINCDS =~ 90 bases ?
        my $minbaseover= ($alnmax >= $TINYALNCDSLEN)? $TINYALNBASES : $alnmax * 0.10; ##$MIN_PCTOVERLAP;
        $skiptinyaln=($aln < $minbaseover)?1:0;
      } else {
        $isfrag= "frag" unless($isfrag); #or "tiny" ? # dont skip assign alt/frag to tiny cds, but maybe always set isfrag ?
      }
    } else { # old, pct-TINYALN is a problem for large cds  
      $skiptinyaln= ($pal < $TINYALN)?1:0;
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
      ## old $pid >= $PHI       
      # old0:  $aclass= ($pal<$PHIALN or $isfrag)?"parthi":"althi";
    	# TEST3: add min-basediff here to PHIALN tests, so shorty tr dont take over w/ minor base diff
    	# NO, see samesize: add $samealn flag, althi0 or althi2 ?? for ~identical size and ~perfect align
      # if($pid >= 99) { ## PHI == 99 default
      #   #old2# if($samealn) {$aclass .= "2"; } elsif # bad calls here, pal lowish 80%
      #   #old3.use this# if($pal >= $PHIALN or $tinyalndiff) { $aclass .= "1"; }  # "100"; use or not? PHIALN reused: was 99 here
      #   #old1 # $aclass .= "1" if($pid > 99 and $pal >= 99);  # "100"; use or not? PHIALN reused: was 99 here
      #   }


=item samesize/samebestmatch 
  problem of 2+mains w/ several essentially identicals, 
  after 2 mains created, then row linking them needs to break tie

egrep '^Anofunz4hEVm00398[2345]t1       ' $pt.tgclass | egrep -v 't[23456789]|t..       ' 
Anofunz4hEVm003982t1	okay	main	Anofunz4hEVm003985t1	100/100/./altmap98	629,94%,complete	aaref:1108,AGAP001965-PA,refbest,chrmap:98a,100i,1890l,5x,KB669169:464434-466658:+,pflag:0
Anofunz4hEVm003983t1	okay	main	Anofunz4hEVm003984t1	100/100/./altmap98	629,94%,complete	aaref:1101,AGAP001965-PA,refgood,chrmap:98a,99i,1890l,5x,KB669169:464434-466658:+,pflag:0
Anofunz4hEVm003984t1	okay	althi1	Anofunz4hEVm003983t1	100/100/./altmap98	629,94%,complete	aaref:1100,AGAP001965-PA,refgood,chrmap:98a,98i,1890l,5x,KB669169:464434-466658:+,pflag:0
Anofunz4hEVm003985t1	okay	althi1	Anofunz4hEVm003982t1	100/100/./altmap98	629,94%,complete	aaref:1108,AGAP001965-PA,refbest,chrmap:98a,100i,1890l,5x,KB669169:464434-466658:+,pflag:0

egrep '^Anofunz4hEVm00398[2345]t1       ' evg24mergeanofun.tgalntab | egrep -v 't[23456789]|t.. ' | sort -k2,2nr -k7,7nr -k6,6n
Anofunz4hEVm003982t1	1890	1996	Anofunz4hEVm003985t1	1890	1996	1888	1883	3602  # main1
Anofunz4hEVm003985t1	1890	1996	Anofunz4hEVm003982t1	1890	1996	1888	1883	3602  # altof main1
Anofunz4hEVm003983t1	1890	1996	Anofunz4hEVm003984t1	1890	1996	1886	1885	3621  # main2
Anofunz4hEVm003984t1	1890	1996	Anofunz4hEVm003983t1	1890	1996	1886	1885	3621  # altof main2
Anofunz4hEVm003983t1	1890	1996	Anofunz4hEVm003985t1	1890	1996	1841	1840	3535  # main2 = altof main1
Anofunz4hEVm003985t1	1890	1996	Anofunz4hEVm003983t1	1890	1996	1841	1840	3535  # altof main1 = main2
Anofunz4hEVm003982t1	1890	1996	Anofunz4hEVm003983t1	1890	1996	1811	1811	3483  # main1 = main2, break tie
Anofunz4hEVm003983t1	1890	1996	Anofunz4hEVm003982t1	1890	1996	1811	1811	3483  .. more, dont redo break
Anofunz4hEVm003982t1	1890	1996	Anofunz4hEVm003984t1	1890	1996	1800	1800	3462
Anofunz4hEVm003984t1	1890	1996	Anofunz4hEVm003982t1	1890	1996	1800	1800	3462
Anofunz4hEVm003984t1	1890	1996	Anofunz4hEVm003985t1	1890	1996	1800	1799	3456
Anofunz4hEVm003985t1	1890	1996	Anofunz4hEVm003984t1	1890	1996	1800	1799	3456

=cut

  if(TEST3) {
    if(my $ocl= $class{$td}) {  
    
  if(TEST1602) {  
      # test2: next is bad here, for class{td} = class{qd} == main both
      #* this helps, reduces main: 10632/17263, +noclass: 3964/4141, but adds those to okay.alt: 90048/60193
      #* need to fiddle w/ alt drops now
      ## maybe should not NEXT here, but set ismain{qd},class{qdmain} ??
      
      #** this may be LOSING valid mains, set in prior row, now its alt is resetting main to althi1 ..
      if($ocl =~ m/^alt|^part/) { $class{$td}= $aclass if($aclass eq "althi1"); next; } # test3
      elsif($ocl =~ m/^main/) { 
        # **  problem of 2+mains w/ several essentially identicals,  break tie here **
        if($qclass =~ m/^main/ or $tisaltpar) { 
          # skip to below aclass set, test5, yes, do this also?
          } 
        else { next; } # * YES, need this : test4
      } 
  } else {
      $class{$td}= $aclass if($ocl =~ m/^alt/ and ($aclass eq "althi1")); # is this right?  or  $aclass eq "althi2"??
      next;   
  } 
    }
  }    

    ## 1602 update: $ismain{$qd} and $ismain{$td} already called, later row sez they are equal also, reclass 1.
    ## .. but below changes 2nd main to alt: $class{$td}=$aclass

    if($aclass) { 
      ## FIXME ismain{qd} >> qd may be alt of other main, .. follow chain here?  keep qd as bestmatch? both?
      my $qmain=$qd; my $attr= $pidalnval; ## "$pid/$pal$antiflag";

      ## FIXME: this still leaves some alts w/o final main id link; including circular alt1 <=> alt2 links
      ## fix for (TEST1603) altpars, need to stop findmain iter when palign drops
      # my $isaltpar=$tisaltpar; if(TEST1603) { $isaltpar= $tisaltpar;}  ## == ($pidalnval=~/altpar/)?1:0; # 
      
      use constant FINDMAINID => 1;  
      # this is right now: if(FINDMAINID) ..
      if($class{$qd}) {  #  =~ /alt|part/
        my $qnext= $qd; 
        my $more= ($bestmatch{$qnext})?1:0;
        my %qdid=( $qd => 1, $td => 1);
        while($more) { 
          $more=0;
          my($qbest,$qpal)=split",",$bestmatch{$qnext};
          
          if(TEST1603) {
            $qbest=0 if($ismain{$qnext} and not $ismain{$qbest}); # want always? maybe bad?
 		        #? my $qtpal= $havepair{$qnext}{$td} || $havepair{$td}{$qnext}; #?? change this to pal align val for altpar TEST1603?
 		        #? $qbest=0 if($tisaltpar and defined $qtpal and $qtpal < $PLOW); # not good?
          }
          
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
if(TEST1602) {   
      if(1) { # this way works.  if($bestmatch{$td}) must be true from above
      if($bestmatch{$td} =~ /^$qd/) { $bestmatch{$td}="$qd,$attr"; }
      elsif(not $bestmatch{$td}) { $bestmatch{$td}="$qd,$attr"; }
      #? elsif($qmain ne $qd and $bestmatch{$td} !~ m,/,) { $bestmatch{$td} .= "/$qmain"; }
      }    
      # if(0) { # test16: this is BAD, bestmatch{td} always set here, usually to qd w/o updated attr
      #   $bestmatch{$td}="$qd,$attr" unless($bestmatch{$td}); #?keep 1st bestmatch? NOT this test, always have bmatch
      # }
      #^ this change if already bestmatch has effect, is it adding /qmain above, or is it not qd ?  test this:
      # if($ismain{$td}) { 
      #  #? delete $ismain{$td}; #? yes or no? prevent both aligned tr as main, but class{td}=alt does that
      # }
} else {      
      $bestmatch{$td}="$qd,$attr"; # orig way 
}
if(TEST3) {
      ## add to prevent 2+ circular alta <> altb <> altc with no main ..
      $class{$qmain}="main" unless($class{$qmain}); # should always be:  $qc >= $tc//$samesize..
}
    }
      
  } close($inh);

  
  # END:  print trclass table; add more fields to output: aaqual, aablast, tr,aa sizes?
  # FIXME: here, elsewhere create ID main,alt table, with new numeric IDs, old/cur ids, main/alt num
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


# pre-read all of infile ids into validids, before idenityclass
sub readIdsFromAlnTab {
  my($infile, $insorted)= @_;
  my($nids,$inh)=(0,undef);
  # ($ok,$inh)= openRead($infile,1);
  open(IN,$infile) or die "reading $infile"; $inh=*IN; 
  while(<$inh>) {  
    next if(/^Qid|^\W/); chomp; 
    my($qd,$qc,$qw,$td,$tc,$tw,$aln,$iden,$bits)= split"\t"; 
    
    #?? use dupids here to mark as not validids ?
    if( $dupids ) {
      $validids{$qd}++ unless( $dupids{$qd} and $qd ne $dupids{$qd});    
      $validids{$td}++ unless( $dupids{$td} and $td ne $dupids{$td});    
    } else {
      $validids{$qd}++; $validids{$td}++;
    }
    $nids++;
  } close($inh); # rewind if ref() ???
  return $nids;
}


sub readblasttab {
  my ($bother)= @_;
  
  my($lt,$lq,$sa,$sm,$ss,$qd,$wq,$nids)= (0) x 10;
  my($ok,$fh)= openRead($bother);
  # my $fh;
  # if($bother =~ /\.gz/) { open(F,"gunzip -c $bother |") or die $bother; $fh=*F; }
  # elsif($bother =~ /stdin|^-/) { $fh= *STDIN; }
  # else { open(F,$bother) or die $bother; $fh=*F; }
  
  %bspans=();
  my %dupskipids=(); my $dupskipspan=0;
  
  ## FIXMEd: self-only matches need  recording via putspans; should be but dupdrop is droping firstid also

  while(<$fh>) { 
    unless(/^\w/) { next; } 
    #  if(/^# Query/) { ## Unused info
    #  #Notnow# puts($lq,$lt,$sa,$sm) if($lt); 
    #  ($qd)=m/Query:\s+(\S+)/; $wq=(m/len=(\d+)/)?$1:0; }
      
    my @v= split; 
    my($q,$t,$bits,$pctident,$aln,$mis,$gap,@bspan)= @v[0,1,-1,2,3,4,5, 6,7,8,9]; # 6-9 =  q. start, q. end, s. start, s. end, 

    # FIXME: need to parse align parts before summing;
    # only(SPANSUM)  now
    if($lq and $q ne $lq) {  
      putspans($lq) unless($dupskipspan); $nids++; %bspans=(); $dupskipspan=0; 
    }

  ## FIXME: self-only matches need  recording via putspans
    my $dupskip=0;
if(DUPFILTER1) {    ## DF1 active; see also inactive DUPFILTER2
    # checkme: is dupids mainid always in blast set? if not dupskip drops useful items
    # ** DUPID list has ids not in blast/cluster set...
    #   add dupskipids{id}=mainid and check later in validids
    if( $dupids ) { 
      if( $dupids{$q} and $q ne $dupids{$q} ) { 
        $dupskipspan=$q; ## if($q eq $t); # flag to skip putspan
        $dupskipids{$q}++; $dupskip++;   
      } elsif( $q ne $t and $dupids{$t} and $t ne $dupids{$t}) { 
        $dupskipids{$t}++; $dupskip++; 
      }  
      if($dupskip) { $ndupdrop++; }  # next; below, not here  WRONG now, need to skip self putspan unless this dup is main 
    }
}   
    ## UNUSED now: qd, wq : FIXME save ids of no-match
    if($t eq $q) { } ## { $qd=$q unless($qd); $wq=$aln unless($wq); }
    elsif($dupskip==0 and $dupskipspan == 0) { 
 
     $bits= bint($bits);
      # only(SPANSUM)  
      my $aident= _max(0,$aln-$mis-$gap); # other way to calc: $aident = $pctident * $aln;
      #new: my $aident= int(0.5 + $aln*$pctident/100); # better?
      sumblastpart( $q, $t, $bits,$aln,$aident, @bspan); # if($bits > $pMINLOW * $bmax or $bspans{$t}); 
    } 
    $lq=$q; $lt=$t;  
  } close($fh);
  
  #only(SPANSUM)
  putspans($lq) unless($dupskipspan); $nids++; $dupskipspan=0;

  if($ndupdrop>0) { $ndupdrop=scalar(keys %dupskipids); }
  warn "# readblasttab: nids=$nids; ndupdrop=$ndupdrop\n" if $debug;
  return($nids);
}

## 2011.aug BUG here, need to test sb-se outside tb-te spans also
## 2013.aug : IS_CDSALIGN ORIENT or is IMPORTANT : need to know when alt-cds are reversed
sub sumblastpart {
  my( $q, $t, $bits,$aln,$aident, $qb,$qe,$sb,$se) = @_;
  my $or=0; # 2013.aug : ORIENT problem.  if both here are reversed, or=0; if only 1, or=-1
  if($qb > $qe) { ($qb,$qe)= ($qe,$qb); $or=-1; }
  if($sb > $se) { ($sb,$se)= ($se,$sb); $or=($or<0)?0:-1; } #was $or--
  unless($bspans{$t}) { 
    $bspans{$t}=[]; push( @{$bspans{$t}}, [$qb,$qe,$sb,$se,$bits,$aln,$aident,$or]); 
    return; }
  my $ov=0;
  ## 2011oct overslop pct fix
  my $qlen=1+$qe-$qb; my $slen=1+$se-$sb;
  my $qslop= _max($OVSLOP, int($pctOVSLOP*$qlen));
  my $sslop= _max($OVSLOP, int($pctOVSLOP*$slen));
  
  foreach my $sp (@{$bspans{$t}}) {
    my($xb,$xe,$tb,$te,$xbit)= @$sp;
    if($qe < $xb or $qb > $xe) { }
    elsif($qe > $xe and $qb >= $xe - $qslop) { }
    elsif($qb < $xb and $qe <= $xb + $qslop) { }
    else { $ov=1; last; }
    ## add 2011.aug
    if($se < $tb or $sb > $te) { }
    elsif($se > $te and $sb >= $te - $sslop) { }
    elsif($sb < $tb and $se <= $tb + $sslop) { }
    else { $ov=1; last; }
  }  
  ## IS_CDSALIGN add $or to bspans, problems?
  unless($ov) { push( @{$bspans{$t}}, [$qb,$qe,$sb,$se,$bits,$aln,$aident,$or]); }
}


=item readlastz: lastz align general format

test case:
  aabugs4/tsaevgc/daphmag5xbest5/dmag5xau13c2011f
  pt=dmag5xau13c2011

$evigene/scripts/rnaseq/asmrna_dupfilter2.pl -debug -CDSALIGN -tinyaln 35 \
  -aasize inputset/$pt.aa.qual -acdhit tmpfiles/${pt}_cd90.aa.clstr \
  -lastz tmpfiles/$pt-self97.lastz.gz  \
  -outeq tmpfiles/$pt.alnlztab2 -outclass $pt.trclasslz2

$bg/mb/galn/bin/lastz \
 --identity=100 --coverage=25 --step=10 --seed=match15 --notransition --exact=20 --match=1,5 \
 --ambiguous=n --nochain --nogapped  --format=general 'altset1.okboth.cds[multiple]' altset1.okboth.cds

#score  name1              strand1 size1 zstart1 end1  name2               strand2 size2 zstart2 end2 identity idPct  coverage covPct

779     dmvel4xpak25Loc11378t4   +  1083    0    779   dmvel4xpak25Loc11378t10  +   1083  0    779  779/779 100.0%  779/1083  71.9%
1083    dmvel4xpak25Loc11378t10  +  1083    0    1083  dmvel4xpak25Loc11378t10  +   1083  0    1083 1083/1083 100.0%  1083/1083 100.0%
770     dmvel4xpak25Loc11378t5   +  978     208  978   dmvel4xpak25Loc11378t10  +   1083  313  1083 770/770 100.0%  770/978 78.7%
251     dmvel4xpak25Loc11378t1   +  684     0    251   dmvel4xpak25Loc11378t10  +   1083  399  650  251/251 100.0%  251/684 36.7%
330     dmvel4xpak25Loc11378t1   +  684     354  684   dmvel4xpak25Loc11378t10  +   1083  753  1083 330/330 100.0%  330/684 48.2%
324     dmvel4xpak25Loc11378t7   +  753     429  753   dmvel4xpak25Loc11378t10  +   1083  753  1077 324/324 100.0%  324/753 43.0%
    ...
455     dmvel4xpak25Loc11378t3   +  753     0    455   dmvel4xpak25Loc11378t7   +   753   0    455  455/455 100.0%  455/753 60.4%
753     dmvel4xpak25Loc11378t7   +  753     0    753   dmvel4xpak25Loc11378t7   +   753   0    753  753/753 100.0%  753/753 100.0%
456     dmvel4xpak25Loc11378t1   +  684     222  678   dmvel4xpak25Loc11378t7   +   753   297  753  456/456 100.0%  456/684 66.7%
324     dmvel4xpak25Loc11378t10  +  1083    753  1077  dmvel4xpak25Loc11378t7   +   753   429  753  324/324 100.0%  324/753 43.0%
324     dmvel4xpak25Loc11378t5   +  978     648  972   dmvel4xpak25Loc11378t7   +   753   429  753  324/324 100.0%  324/753 43.0%

=cut

sub readlastz {
  my ($bother)= @_;
  
  my $islzg=0;
  my($lt,$lq,$sa,$sm,$ss,$qd,$wq, $nids)= (0) x 10;
  my($ok,$fh)= openRead($bother);
  # my $fh;
  # if($bother =~ /\.gz/) { open( $fh,"gunzip -c $bother |") or die $bother; }
  # elsif($bother =~ /stdin|^-/) { $fh= *STDIN; }
  # else { open($fh,$bother) or die $bother;  }
  %bspans=();

  ## TEST: is lastz output file sorted properly? LOOKS OK //name2 should all be grouped, does split-run undo that?
  
  my @hd; # $islzg=1 ;
  while(<$fh>) { 
    if(!$islzg and m/^#score\tname1/ and /\tcoverage/) { $islzg=1; chomp; s/^#//; @hd=split"\t"; } # should do but want some slack here?
    next unless(/^\d/);     
    chomp; my @v= split"\t";
    unless(@v==15 and $islzg) { die "# ERR: doesnt look like lastz general format I know: hd=",@hd," val=",@v,"\n" ; }
    
    ## allow subset columns?
    #score  name1   strand1 size1 zstart1 end1 name2   strand2 size2   zstart2 end2  identity idPct   coverage   covPct
    my( $lzscore, ## is lz score = count of ident bases? no,
      $tid, $tor, $tsize, $tstart, $tend,
      $qid, $qor, $qsize, $qstart, $qend,
      $nida, $pid, $ncovb, $pcov)= @v;
    $qstart++; $tstart++; # move to 1-origin

    ## NOTE: my lastz out has name2 as first-order == qid, name1 == tid
    if($lq and $qid ne $lq) {  
      putspans($lq);  $nids++; %bspans=();
    }

    if($tid eq $qid) { 
      #?? add:  $trsize{$qid} = $qsize unless($trsize{$qid});
      # $qd=$qid unless($qd); $wq=$aln unless($wq); 
    } else { 
      my($aident,$na)= split "/",$nida;
      my($aln,$nb)= split "/",$ncovb;
      my $bits=  $lzscore; # what? maybe 2*aident? 
      
# if(1) { ## SPANSUM
      ## *?* Need this; lastz hsp as for blastn can overlap lots ..
      sumblastpart( $qid, $tid, $bits,$aln,$aident, $qstart,$qend,$tstart,$tend); 
# } else {      
#       $bspans{$tid}=[]  unless(defined $bspans{$tid}); 
#       push( @{$bspans{$tid}}, [$qstart,$qend,$tstart,$tend,$bits,$aln,$aident]);  
# }
    } 
    ($lq,$lt)= ($qid,$tid);
  } close($fh);
  
  putspans($lq);  $nids++;   
  warn "# readlastz: nids=$nids\n" if $debug;  # ; ndupdrop=$ndupdrop
  return($nids);
}



=item cd-hit(est) align cluster format

  has align-span, percent ident, enough for use.
  
  cacao11pub3ig_cde25.cds.clstr
  >Cluster 0
  0       16221nt, >Thecc1EG029098t1... *
  >Cluster 1
  0       15408nt, >Thecc1EG019010t1... at 1:15408:1:15411/+/99.98%
  1       15411nt, >Thecc1EG019010t2... *
  >Cluster 2
  0       14343nt, >Thecc1EG007738t1... *
  1       10578nt, >Thecc1EG007738t2... at 1:10574:1831:12404/+/100.00%
  >Cluster 4
  0       12732nt, >Thecc1EG015810t1... at 1:12732:304:13035/+/100.00%
  1       13035nt, >Thecc1EG015810t2... *
  2       12504nt, >Thecc1EG015810t3... at 74:12503:148:12584/+/99.94%
  3       12717nt, >Thecc1EG015810t4... at 1:12717:304:13035/+/99.88%
  >Cluster 15
  0       9804nt, >Thecc1EG034527t1... *
  1       9804nt, >Thecc1EG034527t2... at 1:9804:1:9804/+/100.00%
  2       9804nt, >Thecc1EG034527t3... at 1:9804:1:9804/+/100.00%
  3       7512nt, >Thecc1EG034527t4... at 1:7512:2287:9804/+/99.92%
  >Cluster 21
  0       8751nt, >Thecc1EG006991t1... *
  1       8304nt, >Thecc1EG006991t2... at 2894:8304:3336:8751/+/99.83%
  2       8178nt, >Thecc1EG006991t3... at 2894:8178:3336:8630/+/99.74%
  3       8391nt, >Thecc1EG006991t4... at 3239:8377:3356:8494/+/100.00%

  dmag4vel4xfi_cde60.cds.clstr
  >Cluster 8633
  0       1383nt, >dmag4vel4xfik65Loc51t4... *
  1       1383nt, >dmag4vel4xfik75Loc49t1... at 1:1383:1:1383/+/100.00%
  >Cluster 8634
  0       1383nt, >dmag4vel4xfik65Loc96t9... *
  1       1383nt, >dmag4vel4xfik65Loc96t15... at 1:1383:1:1383/+/99.42%
  2       1383nt, >dmag4vel4xfik65Loc96t18... at 1:1383:1:1383/+/99.86%
  3       684nt, >dmag4vel4xfik81Loc3813t5... at 1:684:463:1146/+/99.12%
  4       864nt, >dmag4vel4xfik85Loc1522t3... at 1:864:379:1242/+/100.00%
  5       864nt, >dmag4vel4xfik85Loc1522t5... at 1:864:379:1242/+/99.31%
  6       774nt, >dmag4vel4xfik91Loc1049t2... at 1:774:373:1146/+/99.22%
  >Cluster 8635
  0       873nt, >dmag4vel4xfik45Loc1t9452... at 1:873:70:942/+/99.89%
  1       873nt, >dmag4vel4xfik45Loc1t9453... at 1:873:70:942/+/99.89%
  2       483nt, >dmag4vel4xfik45Loc1t9455... at 1:483:901:1383/+/100.00%
  3       483nt, >dmag4vel4xfik45Loc1t9479... at 1:483:901:1383/+/100.00%
  4       492nt, >dmag4vel4xfik55Loc412t36... at 1:492:892:1383/+/100.00%
  5       1383nt, >dmag4vel4xfik65Loc127t19... *
  6       240nt, >dmag4vel4xfik65Loc127t21... at 1:240:1144:1383/+/100.00%
  7       1383nt, >dmag4vel4xfik75Loc138t9... at 1:1383:1:1383/+/99.86%

=cut

sub readcdhit {
  my ($bother)= @_;
  
  my($lt,$lq,$sa,$sm,$ss,$qd,$wq)= (0) x 10;
  my($ok,$fh)= openRead($bother);
  # my $fh;
  # if($bother =~ /\.gz/) { open(F,"gunzip -c $bother |") or die $bother; $fh=*F; }
  # elsif($bother =~ /stdin|^-/) { $fh= *STDIN; }
  # else { open(F,$bother) or die $bother;  $fh=*F; }
  %bspans=();
  
  my $iscdhit=0;
  my($cli,$ncl,$nalt,$mainid,$mainlen,$nerr,$hasdupid)=(0)x10;
  while(<$fh>) { 
    if(/^>/) { 

    ## move this dupid filter before into readblastab, readcdhit ?
if(DUPFILTER1) {    
    if( $hasdupid > 0 ) {
      foreach my $lt (sort keys %bspans) {
        if( my $dpmain= $dupids{$lt} ) {
          if($bspans{$dpmain} and $lt ne $dpmain) {  # ok to drop ..
            if($lt eq $mainid) { $mainid= $dpmain; }
            delete $bspans{$lt}; $ndupdrop++; 
          }
        }
      }
    }  
}
      putspans($mainid) if($mainid); %bspans=(); $hasdupid=0;
      ($cli)= m/(\d+)/;  $ncl++; #   m/Cluster\s*(\d+)/;  fixed or not?
      
    # elsif(/^(\d\w*)\s+(\d+)(..), >(.+)\.\.\. (.+)$/) # problem: >(ID\.1)\.\.\.
    } elsif(/^(\d+)/) { 
      my $i=$1;
      m/(\d+)(nt|aa), >(.+)\.\.\. (.+)$/; # problem: >(ID\.1)\.\.\.
      my($tlen,$typ,$tid,$pinfo)=($1,$2,$3,$4); 
      ## FIXME: merge.clstr:  1f1  9999aa, >id... at 100  << .f1,.f2 added // DROP this?
      #cdhit perls: /(aa|nt), >(.+)\.\.\./
      my $pi;
      my $ismain=($pinfo =~ /\*$/)?1:0;
      $hasdupid++ if($dupids and $dupids{$tid}); # need to read all cluster then drop dups
      unless($tid) {
        $nerr++;
      } elsif($ismain) {
        $mainid= $tid; $mainlen=$tlen; $pi=100;
        push( @{$bspans{$tid}}, [1,$tlen,1,$tlen,$tlen,$tlen,$tlen]); 
      } else {
        # pi for cdhit-est: at 1:492:892:1383/+/100.00%
        # pi for cdhit-aa : at 99.76%
        $pinfo =~ s/at //;  $pinfo =~ s/^\D+//; #? always \digit start
        ($pi)= $pinfo =~ m/([\d\.]+)%/; $pi=~s/\.00//; $pi||=0;
        my($qstart,$qend,$tstart,$tend)= (0) x 4;
        my @pinfo= split( m=/=, $pinfo);
        my @aln= split/:/,$pinfo[0]; #? multiple align segs or 1 only?
        my $aor=$pinfo[1]; # dont need?
        # 249nt, >dmag5vel5xco1k75Loc13984t1... at 249:1:421:669/-/100.00% << revalign
        if(@aln>3) { ($qstart,$qend,$tstart,$tend)= @aln; # what if @aln>4 ?
          ($qstart,$qend)= ($qend,$qstart) if($qstart>$qend); } 
        my $aln= 1+$qend-$qstart;
        my $aident= int($aln * $pi/100);
        my $bits= $aln; # aident? 
        push( @{$bspans{$tid}}, [$qstart,$qend,$tstart,$tend,$bits,$aln,$aident]); # add tlen?
        $nalt++;
      }
      # $iscdhit=2; 
    } elsif(/^\w/) {
      $nerr++; # warn/die not cdhit format ..
    }
  }
  
if(DUPFILTER1) {    
    if( $hasdupid > 0 ) {
      foreach my $lt (sort keys %bspans) {
        if( my $dpmain= $dupids{$lt} ) {
          if($bspans{$dpmain} and $lt ne $dpmain) {  # ok to drop ..
            if($lt eq $mainid) { $mainid= $dpmain; }
            delete $bspans{$lt}; $ndupdrop++; 
          }
        }
      }
    }  
}
  putspans($mainid) if($mainid); %bspans=();
  warn "# readcdhit: nclust=$ncl, nalt=$nalt, ndupdrop=$ndupdrop, nerr=$nerr \n" if $debug;
  return($ncl+$nalt); 
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


sub correctAAcluster
{
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
    ## sort cids by aasize? or need to check thru all?
    # @cids = sort{$aasize{$b} <=> $aasize{$a}} @cids;
    my $reset=0;  my @goodids=(); my $ninval=0; 
    if($havevalidids and not $validids{$mainid}) { $maw=0; $mqv=-9; $mainid=0;  $ninval++; }
    foreach my $id (@cids) {
      ## also check each id is valid for tr/cds align, drop invalids including current mainid
      if($havevalidids and not $validids{$id}) { 
        ## if($id eq $mainid) { $maw=0; $mqv=-9; $mainid=0; } # above now: main could be last id.. fixme
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
  #   my $fh;
  #   if($bother =~ /\.gz/) { open(F,"gunzip -c $bother |") or die $bother; $fh=*F; }
  #   elsif($bother =~ /stdin|^-/) { $fh= *STDIN; }
  #   else { open(F,$bother) or die $bother;  $fh=*F; }
  
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
      
    # elsif(/^(\d\w*)\s+(\d+)(..), >(.+)\.\.\. (.+)$/)  ## BADD patt ??
    } elsif(/^(\d+)/) {
      my $i=$1;
      m/(\d+)(nt|aa), >(.+)\.\.\. (.+)$/;  
      my($tlen,$typ,$tid,$pinfo)=($1,$2,$3,$4); 
      ## FIXME: merge.clstr:  1.f1  9999aa, >id... at 100  << .f1,.f2 added // DROP this?
      ## FIXME2: new merge fnum at end of line now;
      my @pmore; ($pinfo,@pmore)= split /\t/, $pinfo;
      #cdhit perls: /(aa|nt), >(.+)\.\.\./
      my $pi;
      my $ismain=($pinfo =~ /\*/)?1:0;  # FIXME: drop /$/; new merge fnum at end of line now;
      # FIXME: merge.aa.clster bad; lacks some main * ; skip those for now
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

sub readblatpsl {
  my ($bother)= @_;
  
  my $ispsl=0;
  my($lt,$lq,$sa,$sm,$ss,$qd,$wq, $nids)= (0) x 10;
  my($ok,$fh)= openRead($bother);
  # my $fh;
  # if($bother =~ /\.gz/) { open(F,"gunzip -c $bother |") or die $bother; $fh=*F; }
  # elsif($bother =~ /stdin|^-/) { $fh= *STDIN; }
  # else { open(F,$bother) or die $bother;  $fh=*F; }
  %bspans=();
  
  $ispsl=1 ;
  while(<$fh>) { 
    ## $ispsl=1 if m/^psLayout/; # should do but want some slack here?
    next unless(/^\d/);     
    chomp; my @v= split"\t";
    if(@v==22) { shift(@v); } # ucsc psl starts with a 'bin' field?
    unless(@v==21 and $ispsl){ die "# error: doesnt look like psl format I know" ; }
    my( $mat, $mis, $rep_matches, 
      $qgapw, $tgapw, $orient,
      $qid, $qsize, $qstart, $qend,
      $tid, $tsize, $tstart, $tend,
      $blocksizes, $qstarts, $tstarts,
      )= @v[0..2, 5,7, 8..16, 18..20];
    $qstart++; $tstart++; # move to 1-origin

    ## FIXME: allow for blat split rows per query, happens sometimes.., use bspans as for blast
    if($lq and $qid ne $lq) {  
      putspans($lq); $nids++; %bspans=();
    }

    if($tid eq $qid) { 
      # $qd=$qid unless($qd); $wq=$aln unless($wq); 
    } else { 
      # my $pmat= 100 * $mat / $qsize;
      # next if($pmat < $TINYALN); # NOT NOW with bspan parts
      my $aident= _max(0,$mat-$mis);  # not sure -gap right there, align should be +gap, ident = match
      my $aln= $mat + $qgapw; # is this right? align span includes query gaps, and mismatches
      my $bits=  $mat; # what? maybe 2*mat ? 2*aident? 
      #** should split blocks in qstarts, tstarts to sum up q/t start,end
      ## but dont need see putspans:       $tbit += $xbit; $taln+= $aln; $tidn+= $aident;
      
      $bspans{$tid}=[] unless(defined $bspans{$tid}); ## [$qstart,$qend,$tstart,$tend,$bits,$aln,$aident]; 
      # putspans($qid);  # only 1 pairmatch/span per psl row.
      push( @{$bspans{$tid}}, [$qstart,$qend,$tstart,$tend,$bits,$aln,$aident]); 
    } 
    ($lq,$lt)= ($qid,$tid);
  } close($fh);
  
  putspans($lq); $nids++;  
  return($nids);
}




__END__


#################### DROP, Obsolete classifier ###################

=item obsolete outclusters

sub outclusters {  # DROP this old version
  my($outh)= @_;
  
  # $better{$lq}{$lt}++;
  my @ids= sort keys %validids;
  my %sumscore;
  foreach my $id (@ids) { 
    foreach my $jd (@ids) { next if($id eq $jd); 
    my $v= $better{$id}{$jd}||0; $sumscore{$id} += $v;
    }
  }
  
  ## need header:
  # outrow: join("\t",$typ,$lq,"$aq/$wq",$lt,"$at/$wt",$sa,$diffaln,$sm)."\n"; 
  print $outh join("\t",qw(Cluster Score12 Type Qid Qsize Tid Tsize Align Daln Ident))."\n" unless($head++);
      
  my %didid=(); my $cluid=0;
  my @topid= sort{ $sumscore{$b} <=> $sumscore{$a} } @ids;
  foreach my $topid (@topid) { 
    next if($didid{$topid});
    $cluid++; $didid{$topid}++;  
    my $topscore= $sumscore{$topid} || 0; 
    my @nextid=  sort{ $better{$topid}{$b} <=> $better{$topid}{$a} } keys %{$better{$topid}};
    foreach my $nextid (@nextid) {
      next if($didid{$nextid});
      my $nextscore= $sumscore{$nextid} || 0; 
      # my $bscore= $better{$topid}{$nextid} || 0;
      my $outrow= $outrows{$topid}{$nextid};
      unless($outrow) { 
        if($outrow= $outrows{$nextid}{$topid}) { 
          my @r=split"\t",$outrow;
          @r[1,2,3,4]= @r[3,4,1,2]; 
          my $ty=$r[0]; unless($r[0] =~ s/qsub/tsub/) { $r[0] =~ s/tsub/qsub/; }
          $outrow= join"\t",@r;
        }
      }
      $outrow ||= "na\n"; # outrows has \n
      print $outh "cl$cluid\t$topscore,$nextscore\t",$outrow;
      $didid{$nextid}++;
    }
  }
}

=item obsolete puts

# need to do more before output: cluster all same/subset by ids
sub puts { # DROP this old version
  my($lq,$lt,$sa,$sident)= @_;
  return if($lq eq $lt);
  # $sa -= $sm #?? adjust align by mismatches? No
  my($qsame,$tsame,$diffaln,$adiffqt)= (0) x 10;
  my $typ=""; 
  my $aq= $aasize{$lq}||0; $aq *=3;  # 3*convert to cds-size
  my $at= $aasize{$lt}||0; $at *=3;
  my $wq= $trsize{$lq}||0; 
  my $wt= $trsize{$lt}||0;  
  
  map{ $typ="missingsize" unless($_ > 0); } ($aq,$at,$wq,$wt);
  unless($typ) {
    my $paq= 100*$sa/$wq;
    my $pat= 100*$sa/$wt;
    my $cdq= 100*$aq/$wq;
    my $cdt= 100*$at/$wt;
    
    # UPD 1602: change TINYALN use for MINALIGN bases
    if($paq < $TINYALN and $pat < $TINYALN) { 
    
    } else {
    
    # add more classes here for pALIGN: sim90, sim80, sim70, ..
    # esp. for tiny fragments that mostly fit in bigger tr of good aa qual
    # need also class crappy big-tr of bad utr
    
    # $qsame= $paq >= $MINALIGN; 
    # $tsame= $pat >= $MINALIGN;   
  
    foreach my $pa (90,80,70,60) { # TEST; 
      # this doubles ndups, most at 80% level, but dont know yet if these are useless subset trasm or valid alts/paralogs
      if($qsame==0 and $paq >= $pa) { $qsame=$pa; }
      if($tsame==0 and $pat >= $pa) { $tsame=$pa; }
    }
  
    $diffaln= ($wt > $wq) ? $wq - $sa : $wt - $sa;
    $adiffqt = $aq - $at;
  
    #? check also wcds/wtr for utrbad score?
    # add class: when pCDS < 10%-33% ?, q,t/crappylong, NOT better .. want to trash these cases
    my $qutrbad= ($cdq >= $OK_CDSUTR)?0: ($cdq > $BAD_CDSUTR)?1: 2;
    my $tutrbad= ($cdt >= $OK_CDSUTR)?0: ($cdt > $BAD_CDSUTR)?1: 2;
    my $utrdiffqt= ($qutrbad and not $tutrbad)? -1 : (not $qutrbad and $tutrbad)? +1 : 0;
    
    if($qsame and $tsame) { 
    
      $typ= ($adiffqt==0 and $wq==$wt) ? "same$tsame"
      : ($adiffqt>0) ? "tsubset$tsame" 
      : ($adiffqt<0) ? "qsubset$qsame" 
      : ($utrdiffqt<0)? "qsubset$qsame"  # qutrbad==2 : "qsubsetbad"
      : ($utrdiffqt>0)? "tsubset$tsame"  # tutrbad==2 : "tsubsetbad"
      : ($wq<$wt)? "qsubset$qsame"
      : ($wq>$wt)? "tsubset$tsame"
      : "same$tsame"; 
    }
    ## add bad == 2 class here, for qsame but tutrbad==2, t is bad; vv
    elsif($qsame) { $typ= ($adiffqt <= 0) ? "qsubset$qsame" : "tsubsetlong$qsame"; } #? tutrbad
    elsif($tsame) { $typ= ($adiffqt >= 0) ? "tsubset$tsame" : "qsubsetlong$tsame"; } #? qutrbad
    }
  }
  
  # unless($head) { print join("\t",qw(Type Qid Qsize Tid Tsize Align Daln Ident))."\n"; $head++; }
  if($typ) { 
    my $val= join("\t",$typ,$lq,"$aq/$wq",$lt,"$at/$wt",$sa,$diffaln,$sident)."\n"; 
    # need clustering hash here to pick out best for many q,t subsets
    if($typ =~ /missing/) { $better{$lq}{$lt}=-999999; }
    elsif($typ =~ /tsubset|same/) { $better{$lq}{$lt}++; }
    elsif($typ =~ /qsubset/) {  $better{$lt}{$lq}++; }
    $validids{$lq}++; $validids{$lt}++;
    $outrows{$lq}{$lt} = $val; # if exists $outrows{$lq}{$lt}, *should* be same val, check?
    #NO# print $val;
    } 
}


=cut

__END__

=item problem case 1602, 
 .. two locus misclassed with dupl each
 .. all have align=4374, cds=4371, ident=4368.. varies a bit, trlen varies some
 .. ?? cdhitest -c 1 didnt mark as same due to ident variance.

Loc6707t2	= Loc6159t5, and Loc1420t2 = Loc1380t6
egrep '^(anofun2srr9afvelvk57Loc6159t5|anofun2srr9afvelvk31Loc6707t2|anofun2srr9afvelvk51Loc1420t2|anofun2srr9afvelvk57Loc1380t6)' evg2anofunx1/evg2anofunz4c.trclass | head                                       
anofun2srr9afvelvk31Loc6707t2	okay	maina2	anofun2srro1trinc17807g1t1	99/100/.	1457,96%,complete	0,0,aadup:tidbanofun2srro1ridk67Loc64281,pflag:0
anofun2srr9afvelvk57Loc6159t5	okay	maina2	tranofun2srr4dsoapk27loc24292t1	99/100/.	1457,97%,complete	0,0,aadup:tidbanofun2srro1ridk67Loc64281,pflag:0
anofun2srr9afvelvk51Loc1420t2	okay	maina2	tranofun2srr4dsoapk55loc6655t2	99/100/.	1443,62%,complete	0,0,aadup:tidbanofun2srro1ridk57Loc28893,pflag:0
anofun2srr9afvelvk57Loc1380t6	okay	maina2	tranofun2srr4dsoapk55loc6655t1	99/100/.	1443,62%,complete	0,0,aadup:tidbanofun2srro1ridk57Loc28893,pflag:0
 
 zegrep '^(anofun2srr9afvelvk57Loc6159t5|anofun2srr9afvelvk31Loc6707t2|anofun2srr9afvelvk51Loc1420t2|anofun2srr9afvelvk57Loc1380t6)' \
    evg2anofunx1/tmpfiles/evg2anofunz4c.alntab.gz | sort -k2,2nr -k7,7nr -k6,6n | less
    n=265 for Loc6707t2+
    n=58 for Loc1420t2+

Loc1420t2 = Loc1380t6.. blast alntab    
anofun2srr9afvelvk51Loc1420t2   4329    6896    tranofun2srr4dsoapk55loc6655t2  4329    5733    4332    4329    8312
anofun2srr9afvelvk57Loc1380t6   4329    6896    tranofun2srr4dsoapk55loc6655t1  4329    5733    4332    4329    8312
anofun2srr9afvelvk51Loc1420t2   4329    6896    anofun2srr9afvelvk73Loc8492t1   4329    6615    4332    4327    8300
anofun2srr9afvelvk51Loc1420t2   4329    6896    tidbanofun2srro1ridk27Loc20783  4329    6880    4332    4326    8295
anofun2srr9afvelvk51Loc1420t2   4329    6896    tidbanofun2srro1ridk57Loc28893  4329    6880    4332    4325    8289
anofun2srr9afvelvk57Loc1380t6   4329    6896    tidbanofun2srro1ridk27Loc20781  4329    6880    4332    4328    8306
anofun2srr9afvelvk57Loc1380t6   4329    6896    tidbanofun2srro1ridk27Loc20784  4329    6880    4332    4329    8312
anofun2srr9afvelvk57Loc1380t6   4329    6896    tidbanofun2srro1ridk37Loc30137  4329    6880    4332    4330    8318
anofun2srr9afvelvk57Loc1380t6   4329    6896    tidbanofun2srro1ridk57Loc28895  4329    6880    4332    4331    8324
anofun2srr9afvelvk51Loc1420t2   4329    6896    anofun2srr9afvelvk51Loc1420t4   4329    6896    4332    4326    8295
anofun2srr9afvelvk51Loc1420t2   4329    6896    anofun2srr9afvelvk57Loc1380t1   4329    6896    4332    4324    8283
anofun2srr9afvelvk51Loc1420t2   4329    6896    anofun2srr9afvelvk57Loc1380t2   4329    6896    4332    4326    8295
anofun2srr9afvelvk51Loc1420t2   4329    6896    anofun2srr9afvelvk57Loc1380t3   4329    6896    4332    4329    8312
anofun2srr9afvelvk51Loc1420t2   4329    6896    anofun2srr9afvelvk57Loc1380t4   4329    6896    4332    4321    8266
anofun2srr9afvelvk51Loc1420t2   4329    6896    anofun2srr9afvelvk63Loc1436t2   4329    6896    4332    4325    8289
anofun2srr9afvelvk51Loc1420t2   4329    6896    anofun2srr9afvelvk63Loc1436t4   4329    6896    4332    4326    8295
anofun2srr9afvelvk51Loc1420t2   4329    6896    anofun2srr9afvelvk63Loc1436t5   4329    6896    4332    4324    8283
anofun2srr9afvelvk51Loc1420t2   4329    6896    tranofun2srr4dsoapk55loc6655t1  4329    5733    4330    4321    8273
anofun2srr9afvelvk57Loc1380t6   4329    6896    tranofun2srr4dsoapk55loc6655t2  4329    5733    4330    4321    8273
anofun2srr9afvelvk57Loc1380t6   4329    6896    anofun2srr9afvelvk73Loc8492t1   4329    6615    4330    4319    8262
anofun2srr9afvelvk51Loc1420t2   4329    6896    tidbanofun2srro1ridk27Loc20781  4329    6880    4330    4322    8279
anofun2srr9afvelvk51Loc1420t2   4329    6896    tidbanofun2srro1ridk27Loc20784  4329    6880    4330    4319    8262
anofun2srr9afvelvk51Loc1420t2   4329    6896    tidbanofun2srro1ridk37Loc30137  4329    6880    4330    4318    8256
anofun2srr9afvelvk51Loc1420t2   4329    6896    tidbanofun2srro1ridk57Loc28895  4329    6880    4330    4317    8250
anofun2srr9afvelvk57Loc1380t6   4329    6896    tidbanofun2srro1ridk27Loc20783  4329    6880    4330    4322    8279
anofun2srr9afvelvk57Loc1380t6   4329    6896    tidbanofun2srro1ridk57Loc28893  4329    6880    4330    4323    8285
anofun2srr9afvelvk51Loc1420t2   4329    6896    anofun2srr9afvelvk57Loc1380t6   4329    6896    4330    4318    8256
    ^^ dup main equivalence here..
anofun2srr9afvelvk57Loc1380t6   4329    6896    anofun2srr9afvelvk51Loc1420t2   4329    6896    4330    4318    8256
anofun2srr9afvelvk57Loc1380t6   4329    6896    anofun2srr9afvelvk51Loc1420t4   4329    6896    4330    4318    8256
anofun2srr9afvelvk57Loc1380t6   4329    6896    anofun2srr9afvelvk57Loc1380t1   4329    6896    4330    4322    8279
anofun2srr9afvelvk57Loc1380t6   4329    6896    anofun2srr9afvelvk57Loc1380t2   4329    6896    4330    4324    8291
anofun2srr9afvelvk57Loc1380t6   4329    6896    anofun2srr9afvelvk57Loc1380t3   4329    6896    4330    4321    8273
anofun2srr9afvelvk57Loc1380t6   4329    6896    anofun2srr9afvelvk57Loc1380t4   4329    6896    4330    4319    8262
anofun2srr9afvelvk57Loc1380t6   4329    6896    anofun2srr9afvelvk63Loc1436t2   4329    6896    4330    4317    8250
anofun2srr9afvelvk57Loc1380t6   4329    6896    anofun2srr9afvelvk63Loc1436t4   4329    6896    4330    4318    8256
anofun2srr9afvelvk57Loc1380t6   4329    6896    anofun2srr9afvelvk63Loc1436t5   4329    6896    4330    4316    8244
anofun2srr9afvelvk51Loc1420t2   4329    6896    tranofun2srr4dsoapk25loc36821t1 4248    6805    4253    4249    8155
anofun2srr9afvelvk57Loc1380t6   4329    6896    tranofun2srr4dsoapk25loc36821t1 4248    6805    4253    4250    8160
anofun2srr9afvelvk51Loc1420t2   4329    6896    tranofun2srr4dsoapk31loc38144t1 4242    6740    4247    4243    8143
anofun2srr9afvelvk57Loc1380t6   4329    6896    tranofun2srr4dsoapk31loc38144t1 4242    6740    4247    4244    8149
anofun2srr9afvelvk51Loc1420t2   4329    6896    anofun2srr9afvelvk31Loc1478t1   3942    6912    3945    3941    7562
anofun2srr9afvelvk51Loc1420t2   4329    6896    anofun2srr9afvelvk41Loc1369t1   3942    6990    3945    3936    7533
anofun2srr9afvelvk57Loc1380t6   4329    6896    anofun2srr9afvelvk31Loc1478t1   3942    6912    3943    3935    7535
anofun2srr9afvelvk57Loc1380t6   4329    6896    anofun2srr9afvelvk41Loc1369t1   3942    6990    3943    3934    7529
anofun2srr9afvelvk51Loc1420t2   4329    6896    tranofun2srr4dsoapk27loc36801t1 3840    6752    3834    3831    7354
anofun2srr9afvelvk57Loc1380t6   4329    6896    tranofun2srr4dsoapk27loc36801t1 3840    6752    3834    3831    7354
anofun2srr9afvelvk57Loc1380t6   4329    6896    tidbanofun2srro1ridk77Loc26471  2829    5049    2832    2827    5416
anofun2srr9afvelvk51Loc1420t2   4329    6896    tidbanofun2srro1ridk77Loc26471  2829    5049    2830    2822    5395
anofun2srr9afvelvk51Loc1420t2   4329    6896    tidbanofun2srro1ridk27Loc20782  2343    2703    2343    2339    4482
anofun2srr9afvelvk51Loc1420t2   4329    6896    tidbanofun2srro1ridk27Loc20785  2343    2703    2343    2336    4465
anofun2srr9afvelvk57Loc1380t6   4329    6896    tidbanofun2srro1ridk27Loc20782  2343    2703    2343    2339    4482
anofun2srr9afvelvk57Loc1380t6   4329    6896    tidbanofun2srro1ridk27Loc20785  2343    2703    2343    2340    4488
anofun2srr9afvelvk51Loc1420t2   4329    6896    anofun2srr9afvelvk83Loc1440t1   2250    6752    2253    2250    4315
anofun2srr9afvelvk57Loc1380t6   4329    6896    anofun2srr9afvelvk83Loc1440t1   2250    6752    2253    2250    4315
anofun2srr9afvelvk57Loc1380t6   4329    6896    anofun2srr9afvelvk93Loc894t1    582     1252    503     503     968
anofun2srr9afvelvk51Loc1420t2   4329    6896    anofun2srr9afvelvk93Loc894t1    582     1252    502     499     949
anofun2srr9afvelvk51Loc1420t2   4329    6896    tidbanofun2srro1fidk01Loc64631  258     258     258     256     -485
anofun2srr9afvelvk57Loc1380t6   4329    6896    tidbanofun2srro1fidk01Loc64631  258     258     258     256     -485
----


=item classifytr new result1

# /bio/bio-grid/aabugs4/bugs/locust/bestof5x/trsets
# locust1best5 primary alt.tr of cl_main; n=108607 - 8545 in cl_alts,poor; 5446 cl_alts have cl_main
# .. need to filter out trival/aberrant alttr .. use trsize, aasize, close to main.tr
  108607 locust1best5.cl_main.altids
   25968 locust1best5.cl_main.altok.tab;  >=45% aasize of main;
   14754 are aa-complete (incl utrbad), not in cl_alts,poor
    5798 are same size as main; same prot?

   55376 locust1best5.cl_alts.ids
   33284 locust1best5.cl_main.ids
   40761 locust1best5.cl_poor.ids

** need to add other alt-set: needs to classify self-blast + aablast
/bio/bio-grid/aabugs4/bugs/locust/bestof5x/trsets
   locust1best5.cl_main.altids

/bio/bio-grid/aabugs4/tsaevg/trsets
      29 Feb  6 20:48 locust1best5.aa.qual -> .
 33933545 Jan 17 13:19 locust1best5.alntab
 31561602 Feb 13 23:18 locust1best5.alntab2
 1455946 Feb  8 14:24 locust1best5.cl_alts.ids
  853266 Feb  8 14:24 locust1best5.cl_main.ids
 1055964 Feb  8 14:24 locust1best5.cl_poor.ids
 12266926 Feb 13 23:18 locust1best5.class2
 5581623 Jan 19 22:57 locust1best5.classtab
 5939924 Feb  6 11:28 locust1best5.clid2tab
 6166194 Jan 14 15:22 locust1best5.tr.count

# test new classifier
pt=locust1best5
$evigene/scripts/rnaseq/asmrna_dupfilter2.pl -outspan -outeq $pt.alntab2 \
-aa=$pt.aa.qual -tr=aaqual \
-blast=sdoutz/sd-slf95-$pt.blastn.gz  \
-ablast=../aaeval/refdebl/refde-$pt.tall3 > $pt.class2

cat $pt.class2 | cut -f2 | sort | uniq -c | head
71657 drop
1729 maybeok
54816 okay

cat $pt.class2 | cut -f2,3 | sort | uniq -c | head -30
5595 drop       althi
2642 drop       altmid
4919 drop       altmidfrag
27183 drop      noclass
31318 drop      parthi   : check drops for aablast, aaqual; altmid may be valid due to part align to keeper
1729 maybeok    parthi   : all rescued by aablast score; fix this for >=90% of best?
16731 okay      althi
3511 okay       altmid
3432 okay       altmidfrag
19763 okay      main
11379 okay      noclass


=item newer cdsidentclass.sh  2013.02.13

if [ "X" = "X$phi" ]; then phi=99; fi   #below 99 gives many more misclass diff locus
if [ "X" = "X$pmid" ]; then pmid=95; fi
if [ "X" = "X$plow" ]; then plow=80; fi
if [ "X" = "X$mina" ]; then mina=50; fi
if [ "X" = "X$cdsalign" ]; then cdsalign=1; fi
if [ "X" = "X$suf" ]; then suf=clid2tab; fi

  cat $atab | grep -v '^Qid' | sort -k2,2nr -k7,7nr -k6,6n |\
  env mina=$mina cdsw=$cdsalign phi=$phi pmid=$pmid plow=$plow  perl -ne \
'BEGIN{ 
$MINAL=$ENV{mina}||50; $CDSW=$ENV{cdsw}; $PHIALN=65; $PHI=$ENV{phi}||99; $PMID=$ENV{pmid}||90; $PLOW=$ENV{plow}||80; }
chomp; @v=split; ($qd,$qc,$qw,$td,$tc,$tw,$aln,$iden,$bits)=@v; $isfrag= $aclass="";
if($tc > $qc) { ($qd,$td)=($td,$qd); ($qc,$qw,$tc,$tw)=($tc,$tw,$qc,$qw); } #?? always swap to greater cds?
# if($class{$qd} or $class{$td}) { next; } # maybe not. td-noclass > qd-main other
if($class{$td}) { next; }  
if($CDSW) { $ww=($qc>$tc and $tc>0)?$tc:($qc>0)?$qc:$tc;  $isfrag= ($tc < 0.5*$qc)?"frag":""; } 
else { $ww=($qw>$tw and $tw>0)?$tw:($qw>0)?$qw:$tw;  $isfrag= ($tw < 0.5*$qw)?"frag":"";  }
$pid= ($aln<1)?0: int(100*(0.5+$iden)/$aln); $pid=100 if($pid>100);
$pal= ($ww<1)?0 : int(100*(0.5+$aln)/$ww); $pal=100 if($pal>100);
$bestmatch{$qd}="$td,$pid/$pal" unless($bestmatch{$qd});
$bestmatch{$td}="$qd,$pid/$pal" unless($bestmatch{$td});
if($pal < $MINAL) { } #defer: $aclass="noalign$isfrag";  
elsif($tc == 0 and $qc > 0 and $pid >= $PLOW) { $aclass="frag0aa"; } # tc==0 special case "frag0"
elsif( $pid >= $PHI ) { $aclass= ($pal<$PHIALN or $isfrag)?"parthi":"althi"; }
## elsif( $pid >= $PHI ) { $aclass="althi$isfrag"; }
elsif( $pid >= $PMID ) { $aclass="altmid$isfrag"; }
elsif( $pid >= $PLOW ) { $aclass="altlo$isfrag"; }
if($aclass) { $class{$td}=$aclass; $bestmatch{$td}="$qd,$pid/$pal"; $ismain{$qd}++; }
END { foreach $d (sort keys %class) {
 ($q,$pal)=split",",$bestmatch{$d}; $c= $class{$d}; print join("\t",$d,$c,$q,$pal)."\n"; }
foreach $d (sort grep{ not($class{$_}) } keys %ismain) {
 ($q,$pal)=split",",$bestmatch{$d}; print join("\t",$d,"main",$q,$pal)."\n"; }
foreach $d (sort grep{ not($class{$_} or $ismain{$_}) } keys %bestmatch) {
 ($q,$pal)=split",",$bestmatch{$d}; print join("\t",$d,"noclass",$q,$pal)."\n"; }
}' > $nam.$suf


=item redo using -outspantab (quick vs slow for all above)

$evigene/scripts/rnaseq/asmrna_dupfilter.pl -outspan \
 -aa=$pt.aa.count -tr=$pt.tr.count -blast=sd-slf95-$pt.blastn.gz > $pt-s95.alntab

head locust1best5-s95.alntab
Qid                             Qclen   Qtlen   Tid                        Tclen   Ttlen   Align   Ident   Bits
locust1sop4p3k31loc9483t1       369     370     locust1Svel1K31L7036t1     1857    1858    370     368     673
locust1sop4p3k31loc9483t1       369     370     locust1sop4p1k23loc1146t2  1152    4185    370     369     678
locust1sop4p3k31loc9483t1       369     370     locust1sop4p2k23loc30t2    867     3971    370     368     673
locust1sop4p3k31loc9495t1       1632    4114    locust1sop4p0k31loc11571t1 288     687     307     306     563
locust1sop4p3k31loc9495t1       1632    4114    locust1sop4p1k31loc7134t3  3231    5461    2412    2406    4355

# pick longest aa, class those aligning well as (a) alts, (f) fragments (add utrbad class?)
# alts : aa >=60% of 1st, align < 100%
# frags: aa < 60% of 1st, align > 95% min(slen,tlength) : min to drop utrbad class

# FIXME: .alntab from blastn is missing any No-blast-match tr (or just self-match, which should be there).
# .. add nomatch to alntab via asmrna_dupfilter ?

#.............................................
# TEST: need to validate blastn align classes w/ genome-mapped trasm,
#    how many of each class are true/false same-locus,
#    vs paralog vs artifacts (poor genomap, poor intron agreement, ..)
#  -- cacao, daphmag, maybe fungr test cases; 
# Alternate, test by read-trmap (bowtie), for multi-maps, proper pairs, homogeneous? read cover
# I.e., need way to validate false pos/false neg and true pos for trasm that have some blastn identity.
# See aaeffect data sets, cacao & daphmag
/Users/gilbertd/Desktop/dspp-work/cacao/cacao3d/rnas/genes/
cacao align6f:
  cacao3tri1asmcd.aligny.tab	cacao3estasm_cd.aligny.tab	cacao3sopc11nr.aligny.tab ..
  $caca/rnas/genes/  trgmap, genogmap, introns
dmagalign6f:  classify each trset_allcd w/ asmrna_dupfilter, then count agree/disagree for alts == same-locus, paralog, other
	daphmag2nwb_allcd.aligny.tab	daphmag2vel9h_allcd.aligny.tab daphmag3cuf13th3.aligny.tab  ..
  $dmag/rnas/asmrna3/genes/  trgmap, genomap, introns,
  
# ** PROBLEMS: cacao, daph show 1/2 of these subset aligns are on diff scaffolds.. 
#  hi-id align trs are hard to classify as same/diff locus.  Use blastn mismatch as key? no mismatch?
#  use pident >= 99% as strict same-locus trasm classifier?

#.............................................

# sort by big q-aa, hi-align
# fixme: add utrbad tests, pick utrgood vs bad for (near)same-size-aa ; sort -k6,6n to get smaller Ttlen for same clen/align

#test# cat locust1best5-s95.alntab | grep -v '^Qid' | sort -k2,2nr -k7,7nr -k6,6n | perl -ne 
#.............................................
#!/bin/bash
# alnclass.sh

alnset=$*
for atab in $alnset; do {
  nam=`basename $atab .alntab`
  if [ ! -f $nam.classtab ]; then 
  echo "# $atab TO $nam.classtab"

# fixme: if td/tw==0 should call subclass frag or missing, NOT set ww=al;
# below patch will class td/tw==0 as frag

  cat $atab | grep -v '^Qid' | sort -k2,2nr -k7,7nr -k6,6n | perl -ne \
'chomp; @v=split; ($qd,$qc,$qw,$td,$tc,$tw,$al)=@v;
if($class{$qd} or $class{$td}) { next; }
$ww=($qw>$tw and $tw>0)?$tw:($qw>0)?$qw:$tw; 
$pal=($ww<1)?0:int(100*(0.5+$al)/$ww);
$bestmatch{$qd}="$td,$pal" unless($bestmatch{$qd});
$bestmatch{$td}="$qd,$pal" unless($bestmatch{$td});
$aclass="";
if($tc < 0.5*$qc and $al >= 0.95*$ww) { $aclass="frag"; }
elsif($tc >= 0.8*$qc and $al >= 0.95*$ww ) { $aclass="althi"; }
elsif($tc >= 0.5*$qc and $al >= 0.95*$ww ) { $aclass="althifrag"; }
elsif($tc >= 0.5*$qc and $al >= 0.80*$tw and $al < 0.95*$ww) { $aclass="altmid"; }
elsif($tc >= 0.5*$qc and $al >= 0.80*$ww ) { $aclass="altlo"; }
elsif($tc < 0.5*$qc and $al >= 0.80*$ww ) { $aclass="altlofrag"; }
if($aclass) { $class{$td}=$aclass; $bestmatch{$td}="$qd,$pal"; $ismain{$qd}++; }
END { foreach $d (sort keys %class) {
 ($q,$pal)=split",",$bestmatch{$d}; $c= $class{$d}; print join("\t",$d,$c,$q,$pal)."\n"; }
foreach $d (sort keys %ismain) {
 ($q,$pal)=split",",$bestmatch{$d}; print join("\t",$d,"main",$q,$pal)."\n"; }
foreach $d (sort grep{ not($class{$_} or $ismain{$_}) } keys %bestmatch) {
 ($q,$pal)=split",",$bestmatch{$d}; print join("\t",$d,"noclass",$q,$pal)."\n"; }
}' > $nam.classtab

fi

tabsum=`cut -f2 $nam.classtab | sort | uniq -c | perl -pe 's/\n/, /;'`
echo "# classes $nam: $tabsum"

} done

#.................
# recalc v3:

cut -f2 locust1best5-s95.classtab2 | sort | uniq -c | head
4062 althi
4915 althifrag
1590 altlo
10244 altlofrag
10248 altmid
18662 frag
20007 main
23808 noclass
-----
49721 allsubset
93507 total uniqids (27 are classed as main + altlo/..)
drop: frag, segregate: althi/hifrag, maybe altmid
keep as primary trset: main,noclass,altlo, maybe altmid

# calc v1:
cut -f2 locust1best5-s95.classtab | sort | uniq -c 
  9570 alt
  4127 althi        # keep but mark or segregate alt/althi set
  3231 althifrag    # maybe drop these
   629 altlo        # combine altlo/altlofrag, not distinct enough
  8397 altlofrag    # could be various: paralogs, alts, ..
 20723 frag         # drop these
 -----
 46677 allsubset
 20382 mainclass (has alts + frags)
 62362 unclassed (unique-ish) : add to .classtab ? > n=47089 this way, others?
129421 totalaa

#..recalc
cut -f2 locust1best5-s95.classtab2 | grep -v noclass | sort | uniq -c
 10214 alt   : rename altmid
  4061 althi
  4909 althifrag
  1624 altlo
 10242 altlofrag
 18659 frag
 -----
 49709 allsubset
 29007 mainclass
 23809 unclassed (unique-ish)   
102525 total 

1719 alt,althi have same gene id
 376 frag have same gene id (true alts?)
  42 altlo same gene id

cut -f2 locust1tri1_allcd-s95.classtab2 | grep -v noclass | sort | uniq -c
2588 althi
1090 althifrag
 123 altlo
 662 altlofrag
1990 altmid
1048 frag
----------
7501 allsubset
6882 mainclass
7413 unclassed
21796 total ?? low, 50168 in aa

cat locust1tri1_allcd-s95.classtab2 | egrep '    (frag|althi|alt )' | cut -f1 | \
sed 's/$/        /' | ggrep -v -F -f - locust1tri1_allcd.aa.count > locust1tri1_noalt.aa.qual


#locust1best5.aa          n=1000; aw=2082; med=1722; min,max=1307,14259; sw=2082891; sn=24999,24.9
#locust1best5noalt.aa     n=1000; aw=1784; med=1472; min,max=1144,14259; sw=1784691; sn=23075,23
#... vs ...
#locust1tri1_allcd.aa     n=1000; aw=1367; med=1179; min,max=919,6074; sw=1367083; sn=0,0
#locust1tri1_noalt.aa     n=1000; aw=1331; med=1152; min,max=894,6074; sw=1331403; sn=0,0



Eg.,
grep locust1Svel1K23L10030t  locust1best5.aa.count
locust1Svel1K23L10030t1 64      0       64,99%,partial  193   << frag
locust1Svel1K23L10030t4 795     0       795,95%,partial3        2510  << mainclass, best prot
locust1Svel1K23L10030t8 707     0       707,99%,partial 2122   << 2nd mainclass
locust1Svel1K23L10030t11        421     0       421,40%,complete-utrbad 3102
  ^^ has lowish align to locust1Svel1K23L10030t4, 1083/(2510,3102); probably bad assembly.
  
grep locust1Svel1K23L10030t locust1best5-s95.classtab
locust1Svel1K23L10030t1 frag    locust1Svel1K23L10030t4
locust1Svel1K39L26321t1 frag    locust1Svel1K23L10030t4
locust1Svel1K39L37361t1 frag    locust1Svel1K23L10030t4
locust1Svel1K39L51443t1 frag    locust1Svel1K23L10030t8
locust1sop4p0k31loc57894t1      frag    locust1Svel1K23L10030t8
locust1sop4p4k31loc57872t1      frag    locust1Svel1K23L10030t8
locust1tri1loc1401455c0t1       frag    locust1Svel1K23L10030t8
locust1tri1loc240526c0t1        altlofrag       locust1Svel1K23L10030t8
locust1tri1loc892983c0t1        frag    locust1Svel1K23L10030t4
locust1tri1loc897709c0t1        frag    locust1Svel1K23L10030t4
locust1tri1loc904858c0t1        frag    locust1Svel1K23L10030t4
locust1vel5k35Loc6409t1 frag    locust1Svel1K23L10030t4


=item FIXME for crappy long utrbad big.tr

whitefly1vel5k35Loc1t58595 << "locus" with 100k's of alternates must be crappy bin of all repetitive junk
  tr is 112k long but only 6000 longest orf used ..
  -- should skip such crappy big tr, classify in output list
  -- look for pCDS < 10%?
  
grep whitefly1vel5k35Loc1t58595 ../../../tsaevg/trsets/whitefly*.{tr,aa}.count 
whitefly1best3.tr.count:whitefly1vel5k35Loc1t58595       112482  31960   23096   23417   32074   1935    4896
whitefly1best3.aa.count:whitefly1vel5k35Loc1t58595       1869    3       1872,4%,complete-utrbad 112482
flamingo2.% grep whitefly1vel5k35Loc1t58595  whitefly*.dupclust | wc -l
     110


=item prelim tests

# duptrfilter.sh
# messy now; somewhat akin to cd-hit-est but 
# need to use blast-align instead of clustering, to get same tr w/ minor breaks
# or would cd-hit-est do what is needed? prob not. test?
#   this is not right yet; DROPPED LARGEST prots ***
#   problem: locust1vel5k55Loc697t2 = locust1sop4p1k23loc4964t1 longer tr, much shorter 4000 aa ***

# inputs: bestset.aa or bestsetaa.ids from cd-hit alltr.aa; all-source.tr
# outputs: bestset.tr, filtered: 1.bestsetaa.ids only, 2.remove trsame subset (like cd-hit-est, but diff methods)

ncbin=$bg/mb/ncbicpp/bin

pt=whitefly1best3vel1k_cd90
pt=locust1best51k

# s1: collect tr from bestcds.ids, limit to long tr? or do all as prep for TSA submit (drop crap)
# $d=~s/whiteflyvel/whitefly1vel/; s/whiteflyvel/whitefly1vel/; << do id change elsewhere.

gzcat vel?trs/*.tr.gz | cat bestof3/whitefly1best3vel_cd90.ids - | env idp=whitefly min=999 perl -ne \
'BEGIN{$IDP=$ENV{idp}; $MINTR=$ENV{min}||999;} if(/^($IDP\S+)\s*$/) { $ok{$1}=1 } else { \
if(/^>(\S+)/) { $d=$1; $ok=$ok{$d}; if(m/len=(\d+)/) { $ok=0 if($1<$MINTR); } } print if($ok); }' \
 > whitefly1best3vel1k_cd90.tr

gzcat vel*trs/locust1vel{4,5,u.all2_cd}.tr.gz  tr{soap,in}/locust*all2_cd.tr.gz | \
cat bestof5/locust1vel15trisop_cd90.ids - |  env idp=locust min=999 perl -ne \
...
> locust1vel15trisop_cd90.tr

# s2: megablast tr x tr; want high ident only

$ncbin/makeblastdb -in $pt.tr -dbtype nucl -out $pt
$ncbin/blastn -db $pt -query $pt.tr -evalue 1e-19 -perc_identity 99 -dust no -outfmt 7 -out $pt.mblast

# s3: critical, identify "same" transcripts.
# this needs adjusting to make sametr.tab from mblast
# need sizes of both q,t tr to judge sameness : same if align(tq) ~= sizeof(t) or sizeof(q)

# orig
# cat $pt.mblast | perl -ne'if(/^# Query/) { ($wq)=m/len=(\d+)/; ($qd)=m/Query:\s+(\S+)/; } elsif(/^(\w+)/) { ($q,$t,$p,$al,$mi,$ga,@v)=split;  unless($t eq $q) { if($t eq $lt) { $sa+=$al; $ss.=$_; $sm+=$mi+$ga;} else { if($sa/$wq > 0.9 and $wq/$sa > 0.9) { print join("\t","same",$lq,$lqw,$lt,$sa,$lqw-$sa,$sm)."\n";  print $ss; } $ss=$_; $sa=$al; $sm=$mi+$ga; } } $lq=$q; $lqw=$wq; $lt=$t;}' > $pt.sametab

# add tr.count input; AND aa.count here, pick same/subset cases from aasize > trsize ..
cat $pt.mblast | perl -ne\
'if(/^# Query/) { puts($lq,$lt,$sa,$sm) if($lt); ($qd)=m/Query:\s+(\S+)/; $wq=(m/len=(\d+)/)?$1:0; } \
elsif(/^\w/) { ($q,$t,$pi,$al,$mi,$ga,@v)=split; \
if($t eq $q) { $qd=$q unless($qd); $wq=$al unless($wq); } \
else { if($t eq $lt) { $sa += $al; $sm += $mi+$ga; $ss.=$_; } \
 else { puts($lq,$lt,$sa,$sm) if($lt); $sa=$al; $sm=$mi+$ga; $ss=$_;  } } \
$lq=$q; $lt=$t; } END{ puts($lq,$lt,$sa,$sm) if($lt); }
sub puts { my($lq,$lt,$sa,$sm)= @_;
my $wq= $sizes{$lq}||-1; my $wt=$sizes{$lt}||-1;
$qsame=($wq<1)? 0: $sa/$wq > 0.9; 
$tsame=($wt<1)? 0: $sa/$wt > 0.9;   
$adif=($wt > $wq) ? $wq-$sa : $wt - $sa;
if($qsame and $tsame) { $typ="same"; }
elsif($qsame) { $typ="qsubset"; }
elsif($tsame) { $typ="tsubset"; } else { $typ=""; }
if($typ) { print join("\t",$typ,$lq,$wq,$lt,$wt,$sa,$adif,$sm)."\n"; } }'\
> $pt.sametab


# s4: critical, create list of drop-tr ids, decide from sametr.tab which to keep
# fix to pick largest aa,aaqual before largest tr, for same trs
grep same $pt.sametab | perl -ne'($ss,$d,$dw,$t,$tw,$ddt,$mi)=split; $dw{$d}=$dw; $ds{$d}{$t}++; \
END { @d=sort{$dw{$b}<=>$dw{$a}} keys %dw; foreach $d (@d) { $w=$dw{$d}; @s=sort keys %{$ds{$d}}; \
foreach $s (@s) { next if($s eq $d); $sw=$dw{$s}; if($sw<=$w) { $drop=1; $drop{$s}++ unless($keep{$s}); } } \
$keep{$d}++ unless($drop{$d}); } foreach $d (sort keys %drop) { print "drop\t$d\n"; } \
foreach $d (sort keys %keep) { print "keep\t$d\n"; } }' > $pt.samedrop

grep drop $pt.samedrop | cut -f2 | sed 's/$/       /' > $tp.dropids

cat $pt.aa.count | egrep -v '^#|^total' | ggrep -v -F -f $pt.dropids - | sort -k2,2nr | head -1000 | env nam=$pt perl -ne '($aw,$nn)=(split)[1,2]; $aw= $aw-$nn;  $n++; $sw+=$aw; $sn+=$nn; push @aw,$aw; END{ $aw=int($sw/$n); $an=int(10*$sn/$n)/10; @aw=sort{$b <=> $a}@aw ; $md=$aw[int($n/2)]; print "$ENV{nam}\t  n=$n; aw=$aw; med=$md; sw=$sw; sn=$sn,$an\n"; }'

=cut
