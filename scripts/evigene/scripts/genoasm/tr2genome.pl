#!/usr/bin/env perl
# evigene/scripts/genoasm/tr2genome.pl

=item add to evigene pipeline
  
  name: tr2genome ? tr2chr? trcds2genome? 
  
  for use with/after trclassing of evigene/scripts/prot/tr2aacds2.pl 
  re-class transcripts with genome-alignments: 
    a. overlapped mains > altmap
    b. unoverlapped alts > mainmap (mainparalog?)
    
  ## UPD 2016.02 option step 6 
  # 6. genome-map reclassing, *after* 1st outclass using okayset
  # make separate pipeline * use after evgmrna2tsa2.pl publicset, for proper alt matching
  my $chrasm=0; # new opt input file -chrasm mygenome.fasta

  sub blast2genome{};
  sub cdsgmap_maketables{}; includes this overeqcdsloc
 
  if($chrasm) {
    # 6.1. blastn okayset/my.cds to chrasm
    my($cdsgmapblast)= blastn_cds_genome($chrasm, @okayset);
 
    # 6.2. tabulate blast > cdsgmap.tall > cdsgmap.equalgene
    my @xx= cdsgmap_maketables($cdsgmapblast);
 
    # 6.4. classify main/alternate cds, okay & drop subsets, using evigene/rnaseq/asmrna_dupfilter2b.pl
    # combine both cdsblastab + cdsgmapblast.eqgene sources of overlap classing
    ($outclass,$outaln)= asmdupfilter_cds($cdsgmapblast,$aasize,$aaclstr);
 
    # 6.5: new asmdupclass_sum(), asmdupclass_fileset(), ..
  }

=item add ncRNA classing 2016.02

  -- asmrna_dupfilter can do 1st pass, likely want tr2genome mapping for intron measures
  -- putative ncRNA can be pulled from tr2aacds dropset:
        drop.(main|noclass), maybe drop.altmid and frag/part classes
  -- use only subset with no CDS overlap (ie not althi, maybe not altmid)
    .. but will need added full tr/cdna overlap test      
  -- maybe select from only dropset/utrbad+utrpoor subset
  -- can ignore okayset/ as having significant coding potential? (even utrbad/poor)
  -- firmest ncrna class, w/o special tests, would include introns, ie tr2genome exon parts alignments
  
=item add reclassing of poor loci 2016.02
  
  -- work on okayset/ aashort, utrbad/poor, a large, mostly uninteresting group
  -- tr2genome chrmap info to reclass fragments/short near-to better loci, as possible detached exons
  -- chrmap 1-exon things, lacking homology, reclass as possible fragments or TE proteins
  
=item part scripts

  0. prelim 
     run evgtr2aacds.pl  : trclass and okayset/
     run evgmrna2tsa2.pl : publicset/
     
  1.  blastn -query evigene_pubset.cds -db chrasm > evgpubcds_chrasm.blastn
  
  2.  $evigene/scripts/makeblastscore3.pl -tall -spans=2 -oneref=1 $mbsopts -aa $cdssizes \
        evgpubcds_chrasm.blastn > evgpubcds_chrasm.tall2
      mbsopts="-pctover=0.03 -pmin=0.85 "; #?
      
  3. $evigene/scripts/prot/overeqcdsloc.pl -noalt -stranded -minover 10  \ 
       -in evgpubcds_chrasm.tall2 -out evgpubcds_chrasm.eqgenebl

  4. reclassify trset (publicset/) with asmrna_dupfilter2b.pl
        .. main/noclass > alts given hiqual chralign overlaps
        .. alts > main/noclass given hiqual chralign nonoverlap
        .. add splitgene class
        .. include reclass aashort/utrbad subset as culls, possible ncrna
        
  maybe add:
    5.  blastn -query ncrnaset/evg.cdna -db chrasm to look for spliced exons
    
=item author
  
  don gilbert, gilbertd near indiana edu, 2016
  part of EvidentialGene, http://arthropods.eugenes.org/EvidentialGene/
 
=cut

use constant VERSION => '2016.02.14'; 

##-----------------------------------
use FindBin;
use lib ("$FindBin::Bin/../"); # assume evigene/scripts/ layout: main, prot/, rnaseq/, ...
use strict;
use Getopt::Long;
use cdna_evigenesub; # replacing local subs 201405

# cdna_evigenesub globals:
our $EVIGENES=$ENV{EVIGENES} || "$FindBin::Bin/..";  
our $EGAPP='tr2chr';  
our $EGLOG='t2g';
our $dryrun=0; ## $DRYRUN ?
our $DEBUG= $ENV{debug}|| 0;
our $EVGLOGH;

my $CHRBLAST_IDENT= 98; # for cds x cds, change maybe for cds x chr ?
my $CDSBLAST_EVALUE= 1e-19; # is this ok?
my $NCPU= 1;
my $MAXMEM= 1000; # in Mb always 

my ($logfile,$aaseq,$aasize,$cdnaseq,$cdsseq,$chrseq,$aacdseq,$aaclstr,$aablast)= (undef) x 20; 
my (@okayset,@dropset,@inputset,@tmpfiles,@erasefiles); # tidyup file sets

my $optok= GetOptions( 
  #"cdnaseq|mrnaseq|trinput=s", \$cdnaseq, ## not here?
  #"aaseq|aainput:s", \$aaseq, # not here?
  "cdsseq|cdsinput=s", \$cdsseq, ## << primary query input
  "chrseq|chrinput|genomeseq=s", \$chrseq, ## << chr.fasta db input
  "logfile:s", \$logfile,
  # "ablastab=s", \$aablast,   # option traa-refaa.blastp.tall4 table for asmrna_dupfilter2.pl classifier
  # "runsteps=s", \$runsteps,  ## general options here?  -runsteps=lastz,xxx,
  "MINCDS=i", \$MINCDS,  
  "CHRBLAST_IDENT=i", \$CHRBLAST_IDENT, "EVALUE=i", \$CDSBLAST_EVALUE,  
  "NCPU=i", \$NCPU, "MAXMEM=i", \$MAXMEM,  
  "tidyup!", \$tidyup, 
  "dryrun|n!", \$dryrun,
  "debug!", \$DEBUG, ## was $debug, ## $DEBUG now, dont need both
);

die "EvidentialGene tr2genome.pl VERSION ",VERSION,"
  map protein coding sequences to chromosome assembly, to re-classify transcript loci
  by quality of duplicates, fragments, and alternate transcripts.
  See http://eugenes.org/EvidentialGene/about/EvidentialGene_trassembly_pipe.html
Usage: tr2genome.pl -cds evigene.cds_pub.fa -chr mychrasm.fa
  opts: -NCPU=$NCPU  -logfile -tidyup -dryrun -debug 
" unless($optok and $cdsseq and $chrseq); 

openloggit( $logfile, $cdsseq); 
loggit(1, "EvidentialGene tr2genome.pl VERSION",VERSION);
loggit(1, "ERR: unused arguments:",@ARGV) if(@ARGV>0); # do something w/ remaining @ARGV .. warn

$debug= $DEBUG; # dont need both
$tidyup= 1 unless($dryrun||$debug); # default on unless debug|dryrun ?

my $APPblastn=    findapp("blastn"); 
my $APPmakeblastdb= findapp("makeblastdb");

# my $APPfastanrdb= findapp("fastanrdb");
# my $APPcdhitest=  findapp("cd-hit-est"); 
# my $APPcdhit=     findapp("cd-hit");  
# ## .. these should call findevigeneapp() to warn/die if missing
# my $APPcdnabest= findevigeneapp("cdna_bestorf.pl"); # allow ENV/path substitutions?
# my $APPtraa2cds= findevigeneapp("prot/traa2cds.pl");

my $APPmakeblastscore = findevigeneapp("makeblastscore3.pl");
my $APPequalgene = findevigeneapp("prot/overeqcdsloc.pl");
my $APPtrdupfilter= findevigeneapp("rnaseq/asmrna_dupfilter2.pl");  
my $APPaaqual=    findevigeneapp("prot/aaqual.sh"); # replace w/ cdna_evigenesub:makeAaQual()
## FAIL at this point if any apps missing?
#-------------------------------------

sub MAIN_start {}
MAIN: {

  # 6.1. blastn okayset/my.cds to chrasm
  my($cdsgmapblast)= blast2genome($cdsseq, $chrasm);

  # 6.2. tabulate blast > cdsgmap.tall > cdsgmap.equalgene
  my @xx= cdsgmap_maketables($cdsgmapblast);

  ##.. these are tr2aacds steps, can we reuse that instead ?
  # 6.4. classify main/alternate cds, okay & drop subsets, using evigene/rnaseq/asmrna_dupfilter2b.pl
  # combine both cdsblastab + cdsgmapblast.eqgene sources of overlap classing
  ($outclass,$outaln)= asmdupfilter_cds($cdsgmapblast,$aasize,$aaclstr);

#   my $OUTH= (defined $EVGLOGH) ? $EVGLOGH : *STDOUT;
#   asmdupclass_sum($OUTH,$outclass,$aasize); 
#   
#   #??? 5. make final output files from outclass: okay-main, okay-alts, drops 
#   # %perfect_dups = redundant_idset($cdsseqnr); # read headers after rebest, for dupclass_fileset
#   # %perfect_frags= fragment_idset("$cdsseqnrcd1.clstr"); # read headers after rebest, for dupclass_fileset
# 
#   my @outfiles= asmdupclass_fileset($outclass,$cdnaseq,$aaseq,$cdsseq,\%perfect_dups,\%perfect_frags); 
#   loggit(0, "asmdupfilter_fileset=", @outfiles);
# 
#   tidyup_files() if($tidyup);
  
} # MAIN


sub cdsgmap_maketables
{
  my($cdsblast,$cdssize)=@_;
  # $APPmakeblastscore == "makeblastscore3.pl"
  # $APPequalgene == "prot/overeqcdsloc.pl"
  my($cmd,$runerr);
  
  my $insuf='blastn';
  my $cdsbltab  = makename($cdsblast,".btall",$insuf);  
  my $cdseqgene  = makename($cdsblast,".eqgene",$insuf);  
  #my $bslog   = makename($cdsblast,".btall.log",$insuf);
  #my $eqlog   = makename($cdsblast,".eqgene.log",$insuf);

  my $mbopts="-tall -spans=2 -onegenome=g1 -pctover=0.03 -pmin=0.85";
  $cmd="$APPmakeblastscore $mbopts -aasize $cdssize $cdsblast > $cdsbltab"; #  2> $bslog
  $runerr= runcmd($cmd); # unless(-s $cdsbltab) 
  # push @tmpfiles, $bslog; ## outclass is main output file?
  
  ## now work on $cdsbltab, using overeqcdsloc.pl methods : call sep prog or tuck subs here?
  #3. $evigene/scripts/prot/overeqcdsloc.pl -noalt -stranded -minover 10  \ 
  #     -in evgpubcds_chrasm.tall2 -out evgpubcds_chrasm.eqgenebl
  # readTallAlnTab($cdsbltab);  # *STDIN
  # putEqualGeneTab($cdseqgene); # *STDOUT

  my $eqopts="-noalts -strand -minover 10"; # -debug
  $cmd="$APPequalgene $eqopts -in $cdsbltab -out $cdseqgene"; # 2> $eqlog
  $runerr= runcmd($cmd); # unless(-s $cdseqgene)  
  # push @tmpfiles, $eqlog; ## outclass is main output file?
  
  return($cdseqgene,$cdsbltab);
}


sub asmdupfilter_cds
{
  my($cdsblast,$aasize,$aacdhit)=@_;

  my $insuf='aa.qual|aa.size|aa|fasta';
  my $outclass= makename($aasize,".trclass",$insuf);
  my $outaln  = makename($aasize,".alntab",$insuf); # change to .cds.alntab ?
  my $aflog   = makename($aasize,".adupfilt.log",$insuf);
  ## my $aacdhit = makename($aasize,".aa.clstr");
  my $aaclstropt= ( $aacdhit and -s $aacdhit ) ? "-acdhit $aacdhit" : "";
  my $dbg=($debug)? "-debug" : "";

	## opt for cdsblast ==  lastz : $USE_LASTZ global, file=name-self97.lastz
	my $cdsblastop= ( $cdsblast =~ m/blast/ ) ? "-blastab $cdsblast"
			: ( $cdsblast =~ m/lastz/ or $USE_LASTZ ) ? "-lastz $cdsblast" : "-blastab $cdsblast";

  my $aablastopt= ( $aablast and -s $aablast ) ? "-ablastab $aablast" : "";

  my $cmd="$APPtrdupfilter $dbg  -aasize $aasize -CDSALIGN $cdsblastop $aablastopt $aaclstropt"
    ." -outeqtab $outaln -outclass $outclass >$aflog 2>&1";
    
  unless(-s $outclass) {
  my $runerr= runcmd($cmd);
  }
  push @tmpfiles, $outaln, $aflog; ## outclass is main output file?
  return($outclass,$outaln);
}

sub tidyup_files 
{
#   if($tidyup and -s $outclass) {
#     loggit(0, "tidyup output folders: okayset dropset inputset tmpfiles");
#     ## tidyup needs tod/basename($fn) ?? use perl:rename($fn,"$tod/$tf") ? log nfiles moved not each filename?
#     sub tidyup{ my($tod,@td)= @_; mkdir($tod); foreach my $fn (@td) { my $tf=basename($fn); runcmd("mv $fn $tod/$tf") if(-f $fn); } }
#     tidyup("okayset",@okayset);  
#     tidyup("dropset",@dropset);  
#     tidyup("inputset",@inputset);  
#     tidyup("tmpfiles",@tmpfiles);  
#     ## my $rmlist= join" ",grep{ -f $_ } @erasefiles;
#     ## loggit(0,"tidyup erase:",$rmlist) if($rmlist); #too long for log .. chop
#     my @rmlist;
#     foreach my $fn (@erasefiles) { if(-f $fn) { unlink($fn); push @rmlist,$fn; } } # too verbose: runcmd("rm $fn") 
#     if(@rmlist) { my $nrm=@rmlist; my $rml=join" ",@rmlist[0..4]; loggit(0,"tidyup erase: n=$nrm, $rml .."); } 
#   }
}
  
sub blast2genome_opts
{
  my $opts=""; 
  if(my $blopt=$ENV{BLASTCHROPT}) { # reuse BLASTNOPT or new BLASTCHROPT
    $opts .= " $blopt"; 
    $opts .= " -perc_identity $CHRBLAST_IDENT" unless($opts=~/perc_identity/);
    $opts .= " -evalue $CDSBLAST_EVALUE" unless($opts=~/evalue/);
    }
  else { 
    $opts .= " "; # gapped, xdrop_gap?? NOT these:  -ungapped -xdrop_ungap 4 -dust no
    $opts .= " -perc_identity $CHRBLAST_IDENT -evalue $CDSBLAST_EVALUE";
    }  
  return $opts; 
}

sub blast2genome
{
  my($cdsseq,$chrasm)=@_;
  if($NCPU>1 and not $dryrun) { return blast2genome_ncpu($NCPU,$cdsseq,$chrasm); }
 
  my $cdsdb= makename($chrasm,"_db");
  my $cdsbltab= makename($cdsseq,"-$cdsdb$CHRBLAST_IDENT.blastn");
  my $blog= makename($cdsdb,".log","xxxsuf");
  return($cdsbltab) unless(-s $cdsseq);

  my $fmtcmd="$APPmakeblastdb -in $chrasm -dbtype nucl -out $cdsdb -logfile $blog";
  
  my $opts= blast2genome_opts();
  
  $opts.=" -num_threads $NCPU" if($NCPU>1); # -num_threads isnt effective
  my $cmd="$APPblastn -task megablast $opts -outfmt 7 -db $cdsdb -query $cdsseq -out $cdsbltab";
    
  unless(-s $cdsbltab) {
  runcmd($fmtcmd);
  my $runerr= runcmd($cmd);
  }
  push @erasefiles, "$cdsdb.nsq","$cdsdb.nin","$cdsdb.nhr";
  push @tmpfiles, $cdsbltab,$blog;
  return($cdsbltab);    
}

sub blast2genome_ncpu
{
  my($npart,$cdsseq,$chrasm)=@_;
  
  my $cdsdb= makename($chrasm,"_db");
  my $cdsbltab= makename($cdsseq,"-$cdsdb$CHRBLAST_IDENT.blastn");
  my $blog= makename($cdsdb,".log","xxxsuf");
  return($cdsbltab) unless(-s $cdsseq);
  
  my $opts= blast2genome_opts();

  my $blcmd0="$APPblastn -task megablast $opts -outfmt 7 -db $cdsdb "; # add parts: -query $cdsseq -out $cdsbltab
  my $fmtcmd="$APPmakeblastdb -in $cdsseq -dbtype nucl -out $cdsdb -logfile $blog";

  unless(-s $cdsbltab) {
    runcmd($fmtcmd);
    push @erasefiles, "$cdsdb.nsq","$cdsdb.nin","$cdsdb.nhr";

    my $csize= fasize($cdsseq);
    my $spldir= makename($cdsseq,"_blsplit/");
    my $splsize= 1 + int($csize/$npart);

    mkdir($spldir); # dryrun?
    my @splset= fasplit( $cdsseq, $spldir, $npart, $splsize); 
    my (@bloset);
    my $icpu= 0; my $ipart=0;
    foreach my $cds1 (@splset) {
      $ipart++;
      my $cdsbltab1= makename($cds1,"-$cdsdb$CHRBLAST_IDENT.blastn$ipart");
      my $cmd1= $blcmd0 . " -query $cds1 -out $cdsbltab1";  # add parts: -query $cdsseq -out $cdsbltab
      push @bloset, $cdsbltab1;
      my $pid= forkcmd($cmd1);    
      if(++$icpu > $npart) { while (wait() != -1) { }; $icpu= 0; }
    }
    while (wait() != -1) { };
    
    my $cmd= "cat ".join(' ',@bloset)." > $cdsbltab"; 
    my $runerr= runcmd($cmd); #??

    if($debug||$dryrun) {
      push @erasefiles, @splset, @bloset;  
    } else {
      foreach my $fn (@splset, @bloset) {  unlink($fn) if(-f $fn); } 
      rmdir($spldir); 
    }
    
  }

  push @tmpfiles, $cdsbltab,$blog;
  return($cdsbltab);    
}
