#!/usr/bin/env perl
# splign2evgff.pl : convert ncbi splign (mrna x genome) exon alignment table to gene.gff

use strict;
use warnings;
use Getopt::Long;

=item about

  convert ncbi splign to gff
  http://www.ncbi.nlm.nih.gov/sutils/splign/splign.cgi

=item evigene
  evigene/scripts/genes/splign2evgff.pl
  relocated from  evigene/scripts/rnaseq/splign2gff2.pl 

=item FIXME2.ncbi.splign changes IDs :(((

  lcl|Zeamay4bEVm000009t1 ref|NC_024459.1| << added damn pipes and lcl/ref crap, 2016 ncbic++ vers

=item FIXME

  -- need to filter out low qual 2ndary aligns,
    use splog table ** for each trid, order by -score, remove 2nd scores below 99%? 95%? of top
  ** splog score NOT GOOD for complex cases
    -- need also score canonical/standard splice bases; often have fwd-fwd good-splice, rev-rev poor-splice, higher score
    -- AND for split-genes, need align-span/scaffold checks: must keep all align parts that cover mRNA span
      
=cut

use constant VERSION => '20170714'; # '20150301'; 
  # updates from 2013 v: close up part-exons, no/tiny separation NOT intron, but tr-align-gap
  # many more updates, mostly split-part fixes; 
  
  ## TooManyParts: ugh, got >99 parts for some now, need limit
  #  ? all have > MINPARTSPAN 99bp; problem comes from MGap inner gaps allow insert of otherpart exons.
  # megasplit eg: Funhe2EKm023981t1_C1;Split=1/210;Name=Integrin alpha-M;
  
my $debug=0;
my $IDPREFIX="evgs"; 
my $SRC="splign";
my $MRNA_ATTR='clen|offs|oid|aalen|Name';
my $DUPMIN=  0.90; # 0.99; # score percent of best align : change this, splog score is unreliable
my $MININTRON=  20; # close exons w/ less separation == tr-align-gaps not introns.
  # NCBI sez SEQ_FEAT.ShortIntron  Introns should be at least 10 nt long
my $MINEXON= 20; # SEQ_FEAT.ShortExon  Internal coding region exon is too short .. but what is "too"?
my $MINIDENT=  0; # e.g. = 0.90, filter low qual exons to prevent them from bumping higher qual part.
my $FULLCOV = 98; # percent, was MINCOV_PART1 = 98; # 
my $MINPARTSPAN = 99; # base pairs
my $MAXPARTS = 9; # dont know, but 200 is too many :(


use constant MIN_EXON_SPAN => 20; # like minintron, but const?
use constant { seSTART => 1, seEND => 2 }; # "start", "end"; # if($startend == seEND) or ($s eq seEND) ..

my($infile,$output,$inlog,$dryrun,$trinfo,$insorted)= (0) x 10;
$output= undef;

my @saveopt= grep /^\-/, @ARGV;
my $optok= GetOptions(
  "input|tblfile=s", \$infile,
  "logfile=s", \$inlog,
  "trinfo|mrna=s", \$trinfo, # mrna headers from mrna.fasta ?
  "output|gff:s",  \$output,
  #unused# "idprefix=s", \$IDPREFIX,  
  "SRC|source=s", \$SRC,  
  "MININTRON=i", \$MININTRON,
  "MINEXON=i", \$MINEXON,
  "MINIDENT=s", \$MINIDENT, 
  "MINPARTSPAN=s", \$MINPARTSPAN,
  "MAXPARTS=s", \$MAXPARTS,
  # "dryrun|n!", \$dryrun, 
  "sorted!", \$insorted, 
  "debug!", \$debug, 
  );

die "splign2evgff: convert NCBI splign alignment of mRNA x chr to gene location gff, Version",VERSION,
"\nusage: splign2evgff.pl -in evgmrna.splign  -log evgmrna.splign.log [ -mrna evigene.mrna ] [-out evgmrna.gff ]
  -mrna option is mRNA fasta, for Evigene headers with attributes: $MRNA_ATTR
  -source $SRC -sorted -debug 
" unless($optok);

if($MINIDENT) { # want proportion 0..1
  if($MINIDENT>9) { $MINIDENT=$MINIDENT/100; }
  elsif($MINIDENT>1) { $MINIDENT=0.99; }
}

my(%trinfo);
readTrinfo($trinfo) if($trinfo);

my($ntrkeep,$ntrdup,%trkeep,%trhasdup); $ntrkeep=0;
$ntrkeep= readSplignLog($inlog) if($inlog);

my $INH= *STDIN;  # require file for sort ?
# if($infile) { open($INH,$infile) or die "reading $infile"; }  # add gzip open
## sort -k2,2 -k1,1 -k3,3 -k8,8n -k9,9n  input.splign  # << by trid / ipart / scaff / sbeg / send

if($insorted) {
	if($infile and $infile !~ /stdin|^-/) { open($INH,$infile) or die "reading $infile"; }
} else {
	die "ERR: need splign file as input to sort" unless($infile and -f $infile);
  open($INH,"sort -k2,2 -k1,1 -k3,3 -k8,8n -k9,9n $infile |") or die "sort/read $infile";
  ## ^^ maybe wrong sort, try instead -k2 geneid, -k6n genestart, -k4nr ident
  ## .. but then need to collect all of gene, and scaf parts (k1 or k3) to decide which parts are best
}

my $OUTH= *STDOUT;
if(not $output and defined $output and $infile) {
  $output= $infile; $output =~ s/\.\w+$//; $output.=".gff";
}
if($output) { open($OUTH,'>',$output) or die "writing $output"; }  

my($ntrout,$ltid,$ltidi,$lidir,$gaps,$mindels)= (0) x 10; 
my(@xon,%trdup,%trprintdup,%trexons);
my $lgeneid=""; my @geneset=();


=item sort multi-compart/split problems

## FIXME here? sort input.splign by trID? get all idir/trid, but preserve exon order per idir?
## sort -k2,2 -k1,1 -k3,3 -k8,8n -k9,9n  input.splign  # << by trid / ipart / scaff / sbeg / send
#
# +875	FunhetEGm056896t1	Scaffold348	0.933	360	604	963	79954	80313	GG<exon>  	M16RM18RM28RM43RM16RM79RMRM12RM16RM2RM3R2MRM4RM2RM2RM20RM8RM5RM2RM5RM2RM12RM22RM17
#	# .. rev map is slightly different: 77082 .. 77306,77313 
#	-875	FunhetEGm056896t1	Scaffold348	0.932	353	963	611	80313	79961	  <exon>GT	M17RM22RM12RM2RM5RM2RM5RM8RM20RM2RM2RM4RMR2M3RM2RM16RM12RMRM79RM16RM43RM28RM18RM9

## problem w/ sorting splign input: need to keep compartments together, sort by compar span + ident,etc
## .. BUT that way (sort -k2,2 -k1,1) + current methods are missing good but split sub-comparts
## .. found when sort by k2 k6n transcript locations; good eg is DSCAM Funhe2EKm000057t1 

compart sort gives best 1 compart, with middle gap in CDS:
KN805525.1      splign  mRNA    1517847 1525094 65      -       .       ID=Funhe2EKm000057t1;cov=65%,1415/2173;pid=99;nexon=8;splice=26;Target=Funhe2EKm000057t1 1 2173;gaps=742,MGap:866-1607,;gescore=62;splinfo=+654,nd2,sc341.624;aalen=577,79%,complete;clen=2173;offs=95-1828

trlocation sort (k6n) fills gap but messes up sub-comparts
KN805525.1      splign  mRNA    1523822 1525094 37      -       .       ID=Funhe2EKm000057t1_C1;Split=1/6;cov=37%,794/2173;pid=100;nexon=4;splice=14;Target=Funhe2EKm000057t1 1 794;gescore=41;splinfo=+654,nd2,sc341.624;aalen=577,79%,complete;clen=2173;offs=95-1828;oid=Funhe2Exx11m007377t4,Funhe2Eq7m074797t15;Name=Down syndrome cell adhesion molecule protein Dscam2 (97%S)
  ^^ comp#654, 2nd in same KN region (alt exons?)
JXMV01065661.1  splign  mRNA    96      344     11      +       .       ID=Funhe2EKm000057t1_C2;Split=2/6;cov=11%,246/2173;pid=98.8;nexon=1;splice=2;Target=Funhe2EKm000057t1 864 1112;gescore=8;splinfo=+483,nd2,sc48.081;aalen=577,79%,complete;clen=2173;offs=95-1828;oid=Funhe2Exx11m007377t4,Funhe2Eq7m074797t15;Name=Down syndrome cell adhesion molecule protein Dscam2 (97%S)
JXMV01065661.1  splign  mRNA    419     678     8       +       .       ID=Funhe2EKm000057t1_C3;Split=3/6;cov=8%,176/2173;pid=99.4;nexon=2;splice=8;Target=Funhe2EKm000057t1 1113 1289;gescore=11;splinfo=+483,nd2,sc48.081;aalen=577,79%,complete;clen=2173;offs=95-1828;oid=Funhe2Exx11m007377t4,Funhe2Eq7m074797t15;Name=Down syndrome cell adhesion molecule protein Dscam2 (97%S)
    ^^ should be 1 part, same compart #483, and could span tr:864-1638 at pi>98%, but comp #653 has lower qual overlaps
KN805525.1      splign  mRNA    1511068 1511193 6       -       .       ID=Funhe2EKm000057t1_C4;Split=4/6;cov=6%,129/2173;pid=95.7;nexon=1;splice=2;Target=Funhe2EKm000057t1 1290 1418;gescore=4;splinfo=+653,nd2,sc375.179;aalen=577,79%,complete;clen=2173;offs=95-1828;oid=Funhe2Exx11m007377t4,Funhe2Eq7m074797t15;Name=Down syndrome cell adhesion molecule protein Dscam2 (97%S)
    ^^ comp#653 is lower qual than #654 and #483; should skip
JXMV01065661.1  splign  mRNA    1154    1672    10      +       .       ID=Funhe2EKm000057t1_C5;Split=5/6;cov=10%,218/2173;pid=99.1;nexon=2;splice=6;Target=Funhe2EKm000057t1 1419 1638;gescore=12;splinfo=+483,nd2,sc48.081;aalen=577,79%,complete;clen=2173;offs=95-1828;oid=Funhe2Exx11m007377t4,Funhe2Eq7m074797t15;Name=Down syndrome cell adhesion molecule protein Dscam2 (97%S)
KN805525.1      splign  mRNA    1517847 1518249 18      -       .       ID=Funhe2EKm000057t1_C6;Split=6/6;cov=18%,402/2173;pid=99.8;nexon=1;splice=2;Target=Funhe2EKm000057t1 1771 2173;gescore=13;splinfo=+654,nd2,sc341.624;aalen=577,79%,complete;clen=2173;offs=95-1828;oid=Funhe2Exx11m007377t4,Funhe2Eq7m074797t15;Name=Down syndrome cell adhesion molecule protein Dscam2 (97%S)
    ^ comp#654 best at end, need to trim 1st span:1608-1739 using cigar R..,D8 M102< keep good end 1639-1739, follows good c#483:864-1638

cat Funhe2EKm000057t1.splign | egrep -v '[LR]-Gap'
+483    Funhe2EKm000057t1       JXMV01065661.1  0.988   249     864     1112    96      344     TT<exon>GT      M8RM63RM8RM167
+483    Funhe2EKm000057t1       JXMV01065661.1  0.988   84      1113    1196    419     502     AG<exon>GT      M18RM65
+483    Funhe2EKm000057t1       JXMV01065661.1  1       93      1197    1289    586     678     AG<exon>GT      M93
+483    Funhe2EKm000057t1       JXMV01065661.1  0.992   129     1290    1418    917     1045    AG<exon>GT      M30RM98
+483    Funhe2EKm000057t1       JXMV01065661.1  0.988   165     1419    1583    1154    1318    AG<exon>GT      M148RM5RM10
+483    Funhe2EKm000057t1       JXMV01065661.1  1       55      1584    1638    1618    1672    AG<exon>AA      M55
+653    Funhe2EKm000057t1       KN805525.1      0.986   71      795     865     1512626 1512556 AG<exon>GT      M28RM42
+653    Funhe2EKm000057t1       KN805525.1      0.946   92      866     957     1511600 1511509 TA<exon>TA      MRM2RMRM63RM8RM12
+653    Funhe2EKm000057t1       KN805525.1      -       321     958     1278    -       -       <M-Gap> -
+653    Funhe2EKm000057t1       KN805525.1      0.957   140     1279    1418    1511204 1511068 CT<exon>GT      M4RMD3M18RM11RM100
+653    Funhe2EKm000057t1       KN805525.1      0.952   165     1419    1583    1510960 1510796 AG<exon>GT      MRM15RM2RM2RM80RM13RM29RM5RM10
+653    Funhe2EKm000057t1       KN805525.1      0.868   53      1584    1636    1510494 1510448 AG<exon>GT      M18D6M25RM3
+653    Funhe2EKm000057t1       KN805525.1      0.895   105     1637    1741    1505926 1505822 AG<exon>GA      RM24RMRM4RM2RM13RM31RM3R2M12R2M4
+653    Funhe2EKm000057t1       KN805525.1      -       46      1742    1787    -       -       <M-Gap> -
+653    Funhe2EKm000057t1       KN805525.1      0.929   198     1788    1984    1504800 1504603 CT<exon>GA      M20RM4R2M33RM34RM6RM16R2IM23RM3RM2RM20RM21RM2
+654    Funhe2EKm000057t1       KN805525.1      1       158     1       158     1525094 1524937   <exon>GT      M158
+654    Funhe2EKm000057t1       KN805525.1      1       339     159     497     1524856 1524518 AG<exon>GT      M339
+654    Funhe2EKm000057t1       KN805525.1      1       129     498     626     1524395 1524267 AG<exon>GT      M129
+654    Funhe2EKm000057t1       KN805525.1      1       168     627     794     1523989 1523822 AG<exon>GT      M168
+654    Funhe2EKm000057t1       KN805525.1      1       71      795     865     1523729 1523659 AG<exon>GT      M71
+654    Funhe2EKm000057t1       KN805525.1      -       742     866     1607    -       -       <M-Gap> -
+654    Funhe2EKm000057t1       KN805525.1      0.886   132     1608    1739    1520110 1519987 TA<exon>GT      M3RM2RM3R2M2RMRM2RM2D8M102
+654    Funhe2EKm000057t1       KN805525.1      1       31      1740    1770    1519905 1519875 AG<exon>GT      M31
+654    Funhe2EKm000057t1       KN805525.1      0.998   403     1771    2173    1518249 1517847 AG<exon>        M396RM6

=cut

=item FIXME: need more parts scoring

	splign.log score and trkeep not good enougn
 ** splog score NOT GOOD for complex cases
    -- need also score canonical/standard splice bases; often have fwd-fwd good-splice, rev-rev poor-splice, higher score
    -- AND for split-genes, need align-span/scaffold checks: must keep all align parts that cover mRNA span

  ** need 2-pass thru splign table, per trancript (are all parts of trID together?)
    -- accumulate all of trID idir parts: +1, -1, +2, -2 ... then decide which to keep
    -- count canon splices, fwd: AG<exon>GT  rev: AC<exon>CT
    -- tabulate trid:span per idir, 
    
  ** But there are those 1-exon things that map 100s of places, slight diff in align score, no splice score.
  
=cut
		
sub cleanid_old{ my $d=$_[0]; $d=~s/\W/_/g; $d; } # ugh, prfix:id problems
sub cleanid{ my $d=$_[0]; $d=~s/^\w+\|//; $d=~s/\|$//;  $d=~s/[^\w\.]/_/g; $d; } # ugh, prfix:id problems
## cleanid UPD 2016, remove pipe crap ..

sub MAINhere {}

use constant XPOORFILT_EARLY => 1; # filter poor exons; _early: adds poor exons to MGap, otherwise same?

while(<$INH>) {
	next unless(/^[+-]\d/); #? all ok?
	chomp;
	my($idir, $tid, $gref, $pident, $alen, $tb, $te, $gb, $ge, $xsplice, $mismats) = split"\t";
	$tid=cleanid($tid);
	$gref=cleanid($gref); # ncbi added "ncbicrap|ID|" crap; allow '.' in ids

## change this trkeep/skip: dont skip from splog score of trkeep just yet; calc gene score from exons, canon splices, ..
## then decide to keep or drop.  For fwd/rev at same locus, need to look at both first here to decide?
## .. OK now, trkeep has up to MAXTOP=10 dups, regardless of log score, but this filter removes some trash w/ 100s align
  my $keepit= ( $ntrkeep > 0 ) ? $trkeep{$tid}{$idir} || 0 : 1;
  next unless($keepit);
    
	my $tidi= $tid.$idir;  # Trid+100  Trid-100 ..
	## DEFER  putgff() : put into trid array, then decide..
	if($tidi ne $ltidi) { 
		if(@xon) {
			my($geneid,$genescore,$genomespan,$mrnaspan,$genearr)= 
				putgff($ltid,$lidir,$gaps,\@xon,$mindels);
			# returns 0 for skipit...
			if($geneid) {		
			if($geneid ne $lgeneid) {
				processGeneSet($lgeneid, @geneset ) if($lgeneid && @geneset); @geneset=();
				}	
			push @geneset, [$geneid,$genescore,$genomespan,$mrnaspan,$genearr]; #??
			$lgeneid= $geneid;
			}
		}
			
		$gaps=""; $mindels=0; @xon=(); 
		}

  ## fixme: orient calc; NOTE tb,te,gb,ge are '-' for Gaps, other non-exon,non-locations
	# idir: +/- orient, item number
  my($dor)= $idir =~ /^([+-])/; # arghhhh .. - means mRNA is reversed (and tb > te ??)
	## ^^ irrelevant for needs here, use tb>te and gb>ge
 	my($strand,$tor,$gor,$hasloc)= ('.',0,0,0);
	$tor=1; if($tb > $te) { $tor=-1; ($tb,$te)= ($te,$tb); } # dont need rev trspan?
  ## FIXME: need to keep gor,tor and pass on to putgff splice test, or do splice flip here?
  if($gb =~ /\d/ and $ge =~ /\d/) {
  	$hasloc=1;
		$gor=1; if($gb > $ge) { $gor=-1; ($gb,$ge)= ($ge,$gb); } # should match idir; orient?
		$strand= ($gor * $tor < 0) ? '-' : '+';
		
    ## 201707: add exonpoor filter here, moving poor align to M-Gap, not putgff	
if(XPOORFILT_EARLY) {    
    my $xw= 1 + $ge - $gb;
    my $xpoorident= ($MINIDENT and $pident < $MINIDENT)?1:0; # dont have:  ($ixn>0 and $ixn < $#xon)
    my $xpoortiny=  ($MINEXON and $xw < $MINEXON)?1:0;
		if($xpoorident or $xpoortiny) {
		  $xsplice= "M-Gap"; # tb,te ok as is
		  $hasloc= 0; # dont record/print as MGap ..
		  # my $xout= "#px.".join("\t",@$x)."\n"; 
			# push @xonpoor, $xout; 
		}
}

	}
	
	# xsplice:  lbp <exon|gap|..> rbp; gap: M-Gap, L-Gap, what else?
	## if gap, then missing vals: pident, gb, ge, mismats
	my $xtype="dunno"; my $splchar="nnnn"; 
	my ($splscore,$rsplscore)= (0,0);
	if($xsplice =~ /exon/) {
  	$xtype="exon"; 
  	# my($sl,$sr)= $xsplice =~ m/(..).exon.(..)]/;  # AG<exon>GT .. 
  	my $sl = ($xsplice=~m/(..).exon/)?$1:'nn';
  	my $sr = ($xsplice=~m/exon.(..)/)?$1:'nn';
  	map{ $_='nn' if($_ eq '  ') }($sl,$sr);
  	
  	## ?? do splice flip here if gor<0 and +splice, or gor>0 and -splice?
  	## not here..
  	
  	my $splor=($gor < 0)?'-': ($gor > 0)?'+' : ''; ## buggers gor not enough? gor*tor ?
  	# my $splor=$strand; ## buggers gor not enough? gor*tor ?
		$splchar=$sl.$sr.$splor;
		##$splchar.='-' if($gor < 0); # DOES THIS DO IT RIGHT?
		
		## this gor> gor< is probably wrong need tor also? ie splice eq '-' or '+'
		## eg. FunhetEGm000081t3, -strand, 20 exons, all 'AGGT' splices .. lots like this, and reverse:
		##   +strand, spice=ACCT
		## sum: splice=4,AGGT = 10258, splice=0,AGGT'  = 8933
		##     splice=4,ACCT = 56, 'splice=0,ACCT' = 29 << this is not good.
		# Scaffold200     splign  mRNA    215515  223375  100     -       .       ID=FunhetEGm000081t3;
		# cov=100%,3976/3991;nexon=20;splice=0;Target=FunhetEGm000081t3 1 3991;gescore=20;
		# splinfo=+14924,nd0,sc176.854;aalen=1329,99%,partial;clen=3991;offs=3-3989;oid=kfish2eg6velvk55Loc85687t1
		## gmap has Scaffold200:215515-223375:-,cov=100% == Titin.
		
		# right?  MOVE this to putgff scoring all exons, keep splice=splchar
		if($sl eq 'AG') { $splscore+=2; } elsif($sl eq 'AC') { $rsplscore+=2; } 
		if($sr eq 'GT') { $splscore+=2; } elsif($sr eq 'CT') { $rsplscore+=2; } 
		#wrong# if($gor > 0) { $splscore+=2 if($sl eq 'AG'); $splscore+=2 if($sr eq 'GT'); }
		#wrong# elsif($gor < 0) { $splscore+=2 if($sl eq 'AC'); $splscore+=2 if($sr eq 'CT'); }
	} else { 
		($xtype=$xsplice) =~ s/[<>-]//g; 
		$gaps.="$xtype:$tb-$te,"; 
	}	
	

      ## FIXME PROBLEM : skip dups if they cover same/nearly genomic locs, in reverse... need best prot strand
      ## save exons? check $trkeep{$tid}{$d} dup count?
      ## PROBLEM 2 	: split-genes, cant use single top score, need to look at align spans,
      ##   when align parts are disjunct on sep Scaffolds, must keep all.      
        
	if($hasloc) {  # only exons have location
 		# insert comment or what for no-location rows? 
    
    ## M==Match; count R, D, I, other? R=replace, D=del, I=ins;  pident maybe == mmsum/alen
    ##?? distinguish Inserts and Dels ?
    my($mmsum,$indel)=(0,0);
    if(1) {
      while($mismats =~ m/(M\d*)/g) { my $m=$1; $m=($m=~m/(\d+)/)?$1:1; $mmsum+=$m; }
      while($mismats =~ m/([ID]\d*)/g) { my $m=$1; $m=($m=~m/(\d+)/)?$1:1; $indel+=$m; $mindels+=$m; }
    } else { # bad regex??
      while($mismats =~ m/M(\d*)/g) { my $m=$1; $m=1 unless($m); $mmsum+=$m; }
      while($mismats =~ m/[ID](\d*)/g) { my $m=$1; $m=1 unless($m); $indel+=$m; $mindels+=$m; }
    }
    
# 		if($rsplscore > $splscore) { 
# 			# want flipped case? need genesum of rsplscore > splscore
# 			# splice=-$rsplscore ? or splice=$splscore,$rsplscore,$splchar ?
# 		}

    #o# my $attr="Parent=$tid;Target=$tid $tb $te;splice=$splscore,$splchar;align=$mmsum/$alen";
    my $ix=1 + @xon;
    my $attr="Parent=$tid;Target=$tid $tb $te;splice=$splchar;ix=$ix";
    #?? add ;mismat=$mismat string for compart trim/fit tests?
    $attr .=";indel=$indel" if($indel); ## dropped with align= at exon output .. keep?
    $attr .=";align=$mmsum/$alen"; # now drop only xon align=, leave indel=
    
    #NOT here# $attr .=";gor=$gor" if($gor < 0); #?? need only for putgff splice fixup?
    push @xon, [$gref,$SRC,$xtype,$gb,$ge,$pident,$strand,".",$attr];
	}
	
	$ltid=$tid; $lidir=$idir; $ltidi= $tidi;
}

if(@xon) {
	my($geneid,$genescore,$genomespan,$mrnaspan,$genearr)= 
		putgff($ltid,$lidir,$gaps,\@xon,$mindels);
	if($geneid) {
    if($geneid ne $lgeneid) { # upd17
      processGeneSet($lgeneid, @geneset ) if($lgeneid && @geneset); 
      @geneset=();
      }
	  push @geneset, [$geneid,$genescore,$genomespan,$mrnaspan,$genearr]; #??
	  $lgeneid= $geneid;
	}
}
processGeneSet($lgeneid, @geneset ) if($lgeneid && @geneset); @geneset=();

close($OUTH); close($INH);
warn "# splign.gff: ntrout=$ntrout, \n" if $debug;

# MAINend --------------------------------------------
	
=item exon rows

	record 'AG<exon>GT' splice bases around exon for other uses,
	canon splice: fwd: AG<exon>GT  rev: AC<exon>CT
	spaces '  ' for no spl char/gap
	
   my($spl,$spr)= $xsplice =~ m/(..).exon.(..)]/;  
		what are variants on pattern?  ??<exon>?? missing/blank chars?  
     'CC<exon>  ' : blank splice before/after Gap
   need to revcomp xsplices, if or>0, but tor < 0?
   ALSO, can use canonical/std splice bases as qual score; some of splog scores are higher for rev-rev case,
     but those had many fewer std splices .. shifted exon breaks?

=item xsplice frequency

  according to augustus, GT/AG or GC/AG are two valid splice patts. AG-GT,AG-GC here
  547720 AG<exon>GT     canon. splice
   88468 AC<exon>CT      rev canon : when does this occur?
        ^ is this valid when mrna-rev + genome-rev ?
          splign is making some backass map orients, is this the case?
          
  350265 <L-Gap>        gaps
  346391 <M-Gap>
  333549 <R-Gap>
  42532   <exon>  
  
  62924   <exon>GT      half-splices
  61712 AG<exon>  
  22567 AC<exon>  
  17501   <exon>CT
  
  44958 AC<exon>CC      which of these are valid alt spicing?
  26913 TA<exon>CT
  23511 AC<exon>AC
  20438 AG<exon>GC     << count as valid? but occur often w/ many mismatches
  19190 CA<exon>CT
  18062 TA<exon>CC
  17852 AG<exon>AA
  17385 AG<exon>AT
  16846 AG<exon>CT
  16775 AC<exon>GT

=item AC-CT cases

>> this is useful case: various mappers have diff exon set, drop-outs as gap where good exon expression shown.
   .. introns give alternate paths that are not captured by tr-mapping.
   >> take a look at all of cufflinks direct location alt-transcripts, not tr-gmap. does it give same exon drop-outs?
       kf2bcuf2_Gsc74g70000t1, or  fungrk2cuf2_Gsc413g51437t1,
       kf2bcuf2_Gsc74g70000t1 kf2bcuf2_Gsc74g70000t2 kf2bcuf2_Gsc74g70000t3 kf2bcuf2_Gsc74g70000t4
       ^^ these all have same mappings as velv,soap,trin trasm, missing some of valid alt exon/intron cases
       .. splign gets other exon mappings of same tr than gmap or cufl.
       .. cant tell if this is trasm or genomeasm mixup, both?
       well-known gene: human:UniRef50_O75165, Name DnaJ subfamily C member 13; 2200aa : size requires extra exons

-323    FunhetEGm027001t1       Scaffold74      1       185     4619    4435    539792  539608  AC<exon>CT      M185
-323    FunhetEGm027001t1       Scaffold74      1       116     4219    4104    539076  538961  AC<exon>CT      M116
-323    FunhetEGm027001t1       Scaffold74      1       105     3257    3153    533407  533303  AC<exon>CT      M105
-323    FunhetEGm027001t1       Scaffold74      1       102     3152    3051    532137  532036  AC<exon>CT      M102
-323    FunhetEGm027001t1       Scaffold74      1       63      3050    2988    531940  531878  AC<exon>CT      M63
-323    FunhetEGm027001t1       Scaffold74      1       160     2987    2828    531768  531609  AC<exon>CT      M160
-323    FunhetEGm027001t1       Scaffold74      1       83      2569    2487    530868  530786  AC<exon>CT      M83
-323    FunhetEGm027001t1       Scaffold74      1       18      2486    2469    528828  528811  AC<exon>CT      M18
-323    FunhetEGm027001t1       Scaffold74      1       201     710     510     518887  518687  AC<exon>CT      M201
-323    FunhetEGm027001t1       Scaffold74      1       42      509     468     517005  516964  AC<exon>CT      M42
-323    FunhetEGm027001t1       Scaffold74      1       150     467     318     516648  516499  AC<exon>CT      M150
-323    FunhetEGm027001t1       Scaffold74      1       76      317     242     515179  515104  AC<exon>CT      M76

-7557   FunhetEGm021830t1       Scaffold10163   1       106     1592    1487    603990  603885  AC<exon>CT      M106
-7557   FunhetEGm021830t1       Scaffold10163   1       96      1486    1391    602760  602665  AC<exon>CT      M96
-7557   FunhetEGm021830t1       Scaffold10163   1       202     1390    1189    602558  602357  AC<exon>CT      M202
-7557   FunhetEGm021830t1       Scaffold10163   1       128     1188    1061    602156  602029  AC<exon>CT      M128
-7557   FunhetEGm021830t1       Scaffold10163   1       125     1060    936     601905  601781  AC<exon>CT      M125
-7557   FunhetEGm021830t1       Scaffold10163   1       112     935     824     601660  601549  AC<exon>CT      M112
-7557   FunhetEGm021830t1       Scaffold10163   1       72      584     513     598861  598790  AC<exon>CT      M72
-7557   FunhetEGm021830t1       Scaffold10163   1       29      512     484     598669  598641  AC<exon>CT      M29

>> utrbad : gene join?
-7704   FunhetEGm027172t1       Scaffold2       1       132     3124    2993    2562226 2562357 AC<exon>CT      M132
-7704   FunhetEGm027172t1       Scaffold2       1       52      2992    2941    2562482 2562533 AC<exon>CT      M52
-7704   FunhetEGm027172t1       Scaffold2       1       151     2712    2562    2564094 2564244 AC<exon>CT      M151
-7704   FunhetEGm027172t1       Scaffold2       1       71      2561    2491    2564343 2564413 AC<exon>CT      M71
-7704   FunhetEGm027172t1       Scaffold2       1       194     2328    2135    2579844 2580037 AC<exon>CT      M194
-7704   FunhetEGm027172t1       Scaffold2       1       55      2134    2080    2580400 2580454 AC<exon>CT      M55
	
=item non-exon rows

4248 <L-Gap>  == left-end (5p?) gap
4476 <M-Gap>  == middle-gap ?
3795 <R-Gap>  == right-gap 
 284 <poly-A> == what?
 165 <poly-T> == what?

=cut

#-------------------------------

sub dupLoc
{
	my($trid,$exons)= @_;
	
	if($trdup{$trid}++) {
      ## PROBLEM : skip dups if they cover same/nearly genomic locs, in reverse... need best prot strand
      ## this trexons test for overlap is bad; use mrna span?
	  my $hasloc= 0;
	  foreach my $x (@$exons) {
	  	my $xbe=  join "", (split"\t",$x)[0,3,4];
	    # my $xbe= $x->[0] . $x->[3] . $x->[4]; # exact match or overlap test?
	    $hasloc++ if($trexons{$trid}{$xbe});	    
	  }
	  if($hasloc) { $trdup{$trid}--; return $hasloc; } #?? is return mistake? need score of dup before deciding
	}
	return 0;
}


sub overparts {
  my($mbi,$mei,$parts)= @_;
  my $nparts= scalar(@$parts);	
  my @ov=(); my $ovspan=0;
	foreach my $i (0..$nparts-1) { 	
		my($md,$mb,$me)= split/[:-]/, $parts->[$i]->[3]; # mrna span
		if($mbi < $me and $mei > $mb) { 
		  if($mbi >= $mb and $mei <= $me) {
		  	return('inside',$i,0,0); # BUT some inside parts will fit in MGap .. skip all those? BAD idea
      }
		  if($mei > $me and $mbi > $mb) {
		    my $ovs= $me - $mbi;
		    if($ovs > $ovspan) { $ovspan=$ovs;
		    @ov=('start',$i,$mbi,$me+1); # last 2 = trim loc change
		    }
		    #no.return('start',$i,$mbi,$me+1); # last 2 = trim loc change
		  } elsif($mbi < $mb and $mei < $me) {
		    my $ovs= $mei - $mb;
		    if($ovs > $ovspan) { $ovspan=$ovs;
		    @ov=('end',$i,$mei,$mb-1); # last 2 = trim loc change
		    }
		    #no.return('end',$i,$mei,$mb-1); # last 2 = trim loc change
		  } else {
		    return('inside',$i,0,0); # BUT some inside parts will fit in MGap .. skip all those?
		  }
		}
	}
	return @ov if(@ov>0);
  return ('outside',-1,$mbi,$mbi);
}
			

sub  processGeneSet
{
	my( $geneid, @geneset )= @_;
	# geneset row == [$geneid,$genescore,$genomespan,$mrnaspan,$genearr];

	use constant kOVERSLOP => 1; # was 19; # add end-exon trim;  need to check, splign does some align adjusting to spread true tr-align
	use constant kOVERMAXTRIM => 9999; # add end-exon trim ; was 299; no real max?;   
	use constant DROPDUPGENES => 0;
  ## my $MINPART= 99; # == MINPARTSPAN now
	
	if(@geneset == 1) {
		printgene(@{$geneset[0]});
	
	} elsif(@geneset>1) {

  ## TooManyParts BUG: MGap let many parts insert into gaps for low ident but many-dup cases
  ## ID=Funhe2EKm023981t1_C1;Split=1/210;Name=Integrin alpha-M
  ## C1: trg=Funhe2EKm023981t1 1 3060;gaps=3683,MGap:419-1521,MGap:2848-2981,MGap:1548-1780,MGap:1804-2087,MGap:3002-3031,MGap:3061-3457,RGap:3489-4174,MGap:2105-2835,MGap:29-113,
  ##  .. mrna span large, but have only end exons: t1 1 28, <BIG MGap> t1 1522 1547, t1 2088 2104, ..t1 3032 3060
  ## C2: trg=Funhe2EKm023981t1 2 418; C3: trg=Funhe2EKm023981t1 4 1039; C4: trg=Funhe2EKm023981t1 4 420
  ## .. require no mRNA target part overlap ? ie. drop this xstarts/lstarts in favor of mRNA mb,me part overlaps
  # grep -c Split=1  *.gff
  # kfish2n_splign15.gff:1594  Split=9: 0; no trimGenePart()
  # kfish2nsplign15f.gff:5687  Split=9: 382; Split=10: 284; TooManyParts
  # kfish2nsplign15g.gff:3852  Split=9:  35; Split=10: 0; FixedTooManyParts, MAXPART=9,..
  # BUT, still have parts w/ mostly overlapping mrna target spans, exon overlap filtering not right?
  # eg: Funhe2EKm017371t1_C1: t1 1 2349; C2: t1 74 459; C3: t1 127 1011; C4: t1 130 1290
  
    if($MAXPARTS and @geneset > $MAXPARTS) { 
      # sort by gescore; pick top scored MAXPARTS only ? 
      @geneset= sort { $b->[1] <=> $a->[1] or $a->[3] cmp $b->[3] } @geneset;
      @geneset= splice(@geneset,0,$MAXPARTS);
      }

    # 2015 add MINIDENT filter low qual exons  
		  # filter: drop or keep as 2nd set to check after 1st parts collection
		  # in putgff() << or processGeneSet() ?
		
		## revert to old sort by genescore? then add parts inside gap? eg unless overparts(mb,me,@otherpartspans)  
    use constant PARTSORT_SCORE => 1;
    if(PARTSORT_SCORE) {     #old#
      @geneset= sort { $b->[1] <=> $a->[1] or $a->[3] cmp $b->[3] } @geneset; # sort gscore,mspan # BAD??
    } else {
      @geneset= sort { my($ar,$ab,$ae)= split/[:-]/,$a->[3]; my($br,$bb,$be)= split/[:-]/,$b->[3];
                    $ab <=> $bb or $b->[1] <=> $a->[1] } @geneset; # sort mspan<< for parts below
    }
     
# if(DROPDUPGENES) { # this filt should be option, parts check keeps out dups...
# 		my @gsnodup=(); 
# 		for my $gn (@geneset) {  # remove dupl calls at same locus, usual is +dir,-dir for splign
# 			my @exons= grep /\texon/, @{$gn->[4]}; 
# 			my $hasdup= (@exons) ? dupLoc($geneid,\@exons) : 0;
# 			push @gsnodup, $gn unless($hasdup);
# 		}
# 		@geneset= @gsnodup;
# }

		# my($gid1,$gsc1,$gsp1,$msp1,$gff1)=@{$geneset[0]};
		my @gff1= @{$geneset[0]->[4]}; # too messy a structure.. drop all but @gff genearr?
		my($cov1)= $gff1[0] =~ m/cov=(\d+)/; # or align= ??

		## FIXME 15, need to trim end exons for parts to fit .. ie kOVERSLOP test not enough.
		## geneset needs mRNA-loc sort here..
		my @parts=(); # ($geneset[0]);
		if($cov1 < $FULLCOV) { #$cov1 > 1 and ; look for more parts?
			my($md1,$mb1,$me1)= (0,0,0); # split/[:-]/, $geneset[0]->[3]; # or gff1 Target=xxx b e
			foreach my $ig (0 .. $#geneset) {
			  my $gparti= $geneset[$ig];
			  # my($gidi,$gsci,$gspi,$mspi,$gffi)=@{$gparti};
				my($md,$mb,$me)= split/[:-]/, $gparti->[3]; # $mspi; # 
				## FIXME: drop tiny parts, me - mb < MINPARTSPAN
				next if($me - $mb < $MINPARTSPAN); #== MINPARTSPAN; trimGenePart() also checks this after trim

if(1) {				
				my ($ovstartend,$ovpart,$at,$newat)= (@parts>0) ? overparts($mb,$me,\@parts) : ('outside',-1,$mb,$mb);
				if($ovpart >= 0) {
				  my $lastpart= $parts[$ovpart];
				  if($ovstartend eq 'start' or $ovstartend eq 'end') {
  				  my($newstart,$newend, $newpart, $newlast)= trimGenePart( $ovstartend, $at, $newat, $gparti, $lastpart);
				    if(ref $newlast) { if($newlast->[1] == 0) { splice(@parts,$ovpart,1); } else { $parts[$ovpart]= $newlast ; } }
				    if(ref $newpart) { push @parts, $newpart;  }# may be missing..
				    $me1=$newend   if($newend>$me1); 
				    $mb1=$newstart if($newstart>0 and $newstart < $mb1);
				  }
				} else {
				  push @parts, $gparti;  
				  $me1=$me if($me>$me1); $mb1=$mb if($mb < $mb1);
				}
} else {				
#/////////old/////////
				if($me1 == 0 or @parts == 0) { push @parts, $gparti; $mb1=$mb; $me1= $me; }
				elsif($mb + kOVERSLOP > $me1) { push @parts, $gparti;  $me1=$me; }
				elsif($me - kOVERSLOP < $mb1) { push @parts, $gparti;  $mb1=$mb; }
				else {  # end-exon trims
				  # change test? find largest distance, me - me1 or mb1 - mb ?
				  # if($me > $me1 and $mb < $mb1) spans last .. messy
				  # elsif($me > $me1 and $mb > $mb1) trim start
				  # elsif($mb < $mb1 and $me < $me1) trim end
				  # else inside last .. messy, but maybe ok insert missing part .. need part spans check?
				  
				  if($mb + kOVERMAXTRIM > $me1) { # only this with mrna-sort?
				    ## need to change all types, exon, mrna, cds? target+genome starts
				    my ($newstart,$newend,$newpart, $newlast)= trimGenePart( "start", $mb, $me1+1, $gparti, $parts[-1]);
				    #no# $parts[-1]= $newlast if(ref $newlast); # need opt to remove last, replace..
				    if(ref $newlast) { if($newlast->[1] == 0) { pop @parts; } else { $parts[-1]= $newlast ; } }
				    if(ref $newpart) { push @parts, $newpart;  $me1=$newend; }# may be missing..
				  } elsif($me - kOVERMAXTRIM < $mb1) {
				    my($newstart,$newend,$newpart, $newlast)= trimGenePart( "end", $me, $mb1-1, $gparti, $parts[-1]);
				    if(ref $newlast) { if($newlast->[1] == 0) { pop @parts; } else { $parts[-1]= $newlast ; } }
				    if(ref $newpart) { push @parts, $newpart;  $mb1= $newstart; }
				  }
				}
} # TEST old
				
			}
		}

		if(@parts > 1) {
			# unshift @parts, $geneset[0]; # done now
			my @order=();
			my $npart=@parts;
			foreach my $i (0..$#parts) { 	
				my($md,$mb,$me)= split/[:-]/, $parts[$i]->[3]; # mrna span
				push @order, [$i,$mb];
			}
			@order= sort { $a->[1] <=> $b->[1] } @order; # sort by mrna b?
			
			my $pn= 0;
			foreach my $o (@order) { 
				my $ip= $o->[0];
				my $part= $parts[$ip];
				$pn++; #my $pn=$i+1;
				my $geneid= $part->[0]; # $spart[$i]->[0];
				my $newid= $geneid."_C$pn"; 
				## ugh: FunhetEGm000040t3_C0, FunhetEGm000040t3_G2_C1
				## have to replace also $spart[$i]->[0]
				$part->[0]= $newid; # $spart[$i]->[0]= $newid;
				foreach (@{$part->[4]}) {  # @{$spart[$i]->[4]}
					 s/(ID|Parent)=$geneid/$1=$newid/; 
					 if(/\tmRNA/) { s,;,;Split=$pn/$npart;,; } # oops, pn is wrong after reorder
				}
				printgene(@$part); 
			}
		
		} elsif(@parts == 1) {
			printgene(@{$parts[0]});     
		} else {
			printgene(@{$geneset[0]}); # wait, look for more @parts before print; update attr if found; qnd answer
		}
			
		if($debug) {
			foreach my $i (1 .. $#geneset) {
				my($did,$dgenescore,$dgenomespan,$dmrnaspan,$dgenearr)= @{$geneset[$i]};
				print $OUTH "#d$i.",$dgenearr->[0] if($dgenearr); # show mRNA 2..n
			}
		}
	} 
}

# $trinfo= readTrinfo($trinfo) if($trinfo);
# >FunhetEGm033615t1 type=mRNA; aalen=645,23%,complete-utrbad; clen=8340; offs=11-1948; oid=kfish2qf7soapk21loc3t1; organism=Fundulus heteroclitus;
# >FunhetEGm033619t1 type=mRNA; aalen=1366,37%,complete-utrpoor; clen=10794; offs=3242-7342; oid=kfish2qf7soapk21loc3t3; organism=Fundulus heteroclitus;

sub readTrinfo {
	my($trinfo)= @_;
  return {} unless($trinfo);
	## my $MRNA_ATTR='clen|offs|oid|aalen|Name';
  ## look for aalen=  offs=  clen= 
  open(my $INH,$trinfo) or die "reading $trinfo"; 
  while(<$INH>) {
    if(/^>(\S+)/) { my $tid=$1; $tid=cleanid($tid);
   		#o# my $at= join ";", m/((?:$MRNA_ATTR)=[^;\s]+)/g;
   		my $at= join ";", m/((?:$MRNA_ATTR)=[^;\n]+)/g;
   		$trinfo{$tid}= $at || "";
    }
  } close($INH);
	return \%trinfo; #?? global now
}

=item readSplignLog

# spl25/kfish2asm.gdna-kfish2p67vs.mrna.split.3.splog 
# +1	FunhetEGm046329t1	Scaffold0	Ok	-8.248		: has only 1 exon partly aligned, thus low -score
# +2	FunhetEGm046569t1	Scaffold0	Ok	-17.252
# ..
# +8	FunhetEGm057611t1	Scaffold0	Ok	-333.571
# +9	FunhetEGm056896t1	Scaffold0	Ok	-44.051		: + = forward
# -9	FunhetEGm056896t1	Scaffold0	Ok	-44.051		: - = reverse
#
## problem ?? from combining split-compute set: duplicate idir numbers, wont uniqely say which items are linked
## no, each trid is limited to one split-part.  

damn ncbi splign 2016 added these ID pipe bits
  lcl|Zeamay4bEVm000009t1 ref|NC_024459.1|

=cut

sub readSplignLog {
	my($inlog)= @_;
	%trkeep= %trhasdup= ();  # global now
  return 0 unless($inlog);
  my(%tscore,@dkeep);
  open(my $LOGH,$inlog) or die "reading $inlog"; 
  while(<$LOGH>) {
    next unless(/^[+-]\d/);
    chomp; my($idir, $tid, $gref, $okay, $score) = split"\t";
    $tid=cleanid($tid);
    $score = -$score; # largest -score is best.. NOT Always; have +/- often, some hiscore - are wrong dir
    $tscore{$tid}{$idir}= $score;
  } close($LOGH);
  
  my $MAXTOP=10;
  my @tid= sort keys %tscore;
  foreach my $tid (@tid) {
  	## FIXME: score from splign.log may not be best to use alone..
  	## 	also count canonical splices for same tr align fwd/rev at same locus..
  	## change to keep top 10? scores, so split-tr parts are kept (? 10 enough, only care for max 3,4 split parts)
  	
    my ($dtop,@d) = sort{ $tscore{$tid}{$b} <=> $tscore{$tid}{$a} } keys %{$tscore{$tid}};
    my $scoretop= $tscore{$tid}{$dtop};
    my $scut= $DUPMIN * $scoretop;
    $trkeep{$tid}{$dtop}=$scoretop; $trhasdup{$tid}=0;
    # push @dkeep, $dtop;
    my $ntop=1;
    foreach my $d (@d) {
      my $s=$tscore{$tid}{$d};
      if($s >= $scut or $ntop < $MAXTOP) { 
        $trkeep{$tid}{$d}=$s; # push @dkeep, $d;
        $trhasdup{$tid}++; $ntrdup++; $ntop++;
      } 
      ## PROBLEM : skip dups if they cover same/nearly exons, in reverse... need best prot strand
      } 
    $ntrkeep++;
    # $trkeep{$tid}= \@dkeep; #??
  }
  
  warn "# splign.log: ntr=$ntrkeep, ntrdup=$ntrdup\n" if $debug;
  return $ntrkeep;
}

sub printgene {
	my($trid,$genescore,$genomespan,$mrnaspan,$refgene)= @_;

	  ## FIXME: here: defered ID change til here.. from putgff; need to clear %trdup / use new hash
	my $id= $trid; my $changeid=0;
	## this also will be used for split-gene cases, where want same 1 ID but add attribs: Split=1/2,2/2 ..
	if($trprintdup{$trid}++) {  
	  my $idup= $trprintdup{$trid}; $changeid=1;
	  $id= $trid ."_G". $idup; # what ID format for 2ndary matches?  _G2..9 ?
	}
	
	foreach my $gff (@$refgene) {
		$gff =~ s/(ID|Parent)=$trid/$1=$id/ if($changeid);
		## 1707: add exon score fix here, problems elsewhere
		if($gff=~/\t(exon|CDS)\t/){ my @x=split"\t",$gff; $x[5]= pctOf($x[5]); $gff=join"\t",@x; }
		print $OUTH $gff; # has \n newline for now
	}
}

sub putgff {
	my($trid,$idir,$gaps,$exons,$mindels)= @_;

	my $id= $trid;
	my @xon= sort{ $a->[3] <=> $b->[3] or $a->[4] <=> $b->[4] } @$exons;
  my $trattr= $trinfo{$trid} || ""; # source mrna attributes
  # clen= offs= oid= aalen= .. Name= ??

  # 2015 add MININTRON check/close exons w/ tiny/no gap; splign min sep=1, begin = 1 + last end.
  
	## dup align check:
	## FIXME: need genescore of all such dups first, to decide which to keep.
	## NOT HERE, see sub dupLoc()	#	if($trdup{$trid}++) {}
	
	my $tscore= $trkeep{$trid}{$idir} || 0;  $tscore =~ s/^\+//;
	my $tsdup=  $trhasdup{$trid}||0; #debug ??
	
	# make mRNA from exons, print
	## FIXME: mrna length needs gap spans added..
	my($mref,$mgb,$mge,$mor,$maln,$mlen,$mrnab,$mrnae, $mident, $lgb,$lge, $ltb,$lte, $lattr,
		 $mnxon,$msplscore,$rsplscore,$gapspan,$ngaps)=(0) x 20;
		 
	if($gaps) {  #? change gapspan to M-gaps only? ignore or separate LR end gaps? esp for split genes
	  my @gs= split",",$gaps; $ngaps=@gs;
	  foreach (@gs) { my($gt,$gb,$ge)= split/[:-]/,$_; $gapspan += 1+$ge-$gb; } # does this include polyA,polyT ?	  
	  $mlen += $gapspan;
	}
	
	# recalc cov= align + mismatchs, not indels or gaps ?
	$mnxon= @xon; my $gor=0; my $ixn=0;
	my(@xonkeep,@xonpoor);
	foreach my $x (@xon) {
		my($gref,$src,$xtype,$gb,$ge,$pident,$strand,$phx,$attr)= @$x;

    # 2015 add MINIDENT filter low qual exons  
    # 2016: end-exons allow lower ident.. if ($ix==0 or $ix==$#xon) { $minid=$MINID_END; }
    # 2017: add MINEXON .. too many tiny things (cds included) with long introns to match unmatchable seq..
    # FIXME: need to add these xonpoor to MGap, gap annots
    # filter: drop or keep as 2nd set to check after 1st parts collection
    # in putgff() or processGeneSet() ?
unless(XPOORFILT_EARLY) {   # see above in splign reader
    my $xw= 1 + $ge - $gb;
    my $xpoorident= ($MINIDENT and $pident < $MINIDENT and ($ixn>0 and $ixn < $#xon))?1:0;
    my $xpoortiny=  ($MINEXON and $xw < $MINEXON)?1:0;
		if($xpoorident or $xpoortiny) {
		  my $xout= "#px.".join("\t",@$x)."\n"; 
			push @xonpoor, $xout; 
			next; # xon
		}
}

		unless($mref) { $mref=$gref; $mor=$strand; }
		$mgb=$gb if($mgb==0 or $gb<$mgb); $mge=$ge if($ge>$mge);
		
		#o# if( my($splscore)= $attr=~m/splice=(\d+)/ ) { $msplscore+=$splscore; }  # score is # canon-splice-bases both sides of exon: 4,2,0.  div by 4? adjust ends? 2-max
		## add rev-splice score sum? splice=9,-12,... add n-gap count to score? or just gapspan?
		## right?  MOVE this to putgff scoring all exons, keep splice=splchar
		## use rsplscore > msplscore to flip strand, flag as sense=-1 as per gmap.gff
		# if($attr=~m/gor=-1/) { $gor--; } else { $gor++; } #splice= fix or add to splice= flag?
		if( my($sl,$sr)= $attr=~m/splice=(..)(..)/ ) {
			my $gorx= ( $attr=~m/splice=$sl$sr([+-])/ ) ? $1 : 0;
			if($gorx eq '-') { $gor--; } elsif($gorx eq '+') { $gor++; }
			if($sl eq 'AG') { $msplscore+=2; } elsif($sl eq 'AC') { $rsplscore+=2; } 
			if($sr eq 'GT') { $msplscore+=2; } elsif($sr eq 'CT') { $rsplscore+=2; } 
		}
		
		my($tgd,$tgb,$tge)= $attr=~m/Target=(\S+) (\d+) (\d+)/; # my $len= 1+$tge-$tgb;
		$mrnab=$tgb if($mrnab==0 or $tgb<$mrnab);
		$mrnae=$tge if($mrnae==0 or $tge>$mrnae);
		
		my($aln,$len)= $attr=~m/align=(\d+).(\d+)/;
	  $maln+= $aln; $mlen+= $len;
		#?? maybe drop align=100/100 from exon outputs, only mismatches?
		$mident += $pident * $aln;   # pident = 0..1 range, convert to int(100*pident) for output.gff
		
		if($lge and $gb - $lge < $MININTRON) { # close up, how?
	    my $lxbe= $gref.$lgb.$lge;
	    $trexons{$trid}{$lxbe}--;
	    
	    # $gref,$src,$xtype,$gb,$ge,$pident,$strand,$phx,$attr
	    # .. update what else? pident? attr.parts?
	    # .. should merge splices: splice1=AGCC+,splice2=ATGT+ >  AG{x1.gap.x2}GT+
	    ## 	need to change $attr=~m/Target=(\S+) (\d+) (\d+)/; .. add 2nd Target ? not good
	    
      #bad?# $attr=~s/Target=$tgd $tgb $tge/Target=$tgd $ltb $lte,$tgb-$tge/; # will this do? probably not
      ## these are same as Mgap annots..
      my($ub,$ue,$ugb,$uge)= ($ltb > $tge) ? ($tgb,$lte,$tge+1,$ltb-1) : ($ltb,$tge,$lte+1,$tgb-1);
      $attr=~s/Target=$tgd $tgb $tge/Target=$tgd $ub $ue;mgap=$ugb-$uge/; 
      my($lsp)= $lattr=~m/splice=(\w\w)/;
      if($lsp) { $attr=~s/splice=(\w\w)/splice=$lsp/; }
	    $x->[8]= $attr;      
	    $x->[3]= $gb= $lgb;
	    pop(@xonkeep); push @xonkeep, $x; #  $xonkeep[-1]= $x; 
		} else {
		  push @xonkeep, $x;
		}
		
	  my $xbe= $gref.$gb.$ge;# my $xbe= $x->[0] . $x->[3] . $x->[4]; # exact match or overlap test?
    $trexons{$trid}{$xbe}++; # for trdup test, only need if($trhasdup{$trid})
    $ixn++;
    ($lgb,$lge)=($gb,$ge);  ($ltb,$lte)=($tgb,$tge); $lattr=$attr;
	}
	@xon= @xonkeep; $mnxon= @xon;
	$maln ||= 1; # div 0 err for no @xons
	if($mnxon < 1) {
	  return(0,0,0,0,[]);
		# return ($trid,$genescore,$genomespan,$mrnaspan,\@gene); 
	}
	
	## FIXME3: is mis-orient problem in part due to bad ORF-call, mRNA is not right?  check for utrbad+misorient
	
	## FIXME2: need to flipor for both ways, test is strand vs fwdspl,revspl:
	##	main: $strand= ($gor * $tor < 0) ? '-' : '+';
	##  need genome-or only here, if fwdspl, expect gor > 0, flip if gor<0; vs revspl && gor < 0
	##  flipor=1 if(($strand eq '-' and fwdspl) or ($strand eq '+' and revspl))
	
  ## FIXME4: if both +score and -score are high, report in addat confusion == rev gene join likely
  ## FIXME5: 16.10.11 : sense=-1; flipor CDS can be bad, wrong calc from exons, offset

	my $flipor=0; my $addat="";
  my $msplscoresave= $msplscore;
	if($msplscoresave > 10 and $rsplscore > 10) {
	  $addat.=";splicemix=$msplscoresave,$rsplscore"; #? or append to splice=$msplscore, ?
	}
	if($rsplscore > 2+$msplscore and $rsplscore > 2) {
		$msplscore= $rsplscore;
		if( $gor > 0 ) { 
		  if($mor eq '+') { $mor='-'; $flipor=1; $addat.=";sense=-1"; }
		  elsif($mor eq '-') { $addat.=";sense=-1";  }
		} elsif($gor < 0 and $mor eq '+') {
		  $mor='-'; $flipor=1; $addat.=";sense=-1"; 
		}
	} elsif( $msplscore > 2 ) { # messier, sense=-1 when mrna-flipped (tor <0); do flip both ways
		if( $gor < 0 ) {
		  if($mor eq '+') { $mor='-'; $flipor=1; $addat.=";sense=-1"; } 
		} elsif($gor > 0 and $mor eq '-') {
		  $mor='+'; $flipor=1; $addat.=";sense=-1"; 
		}
	}
	
	if($trattr =~ /clen=(\d+)/) { $mlen=$1;	}
	$mlen||=1; #warn ?? unless($mlen);  # dies mlen==0 if no -trinfo mrna.hdr
	my $pctalign= pctOf( $maln/$mlen ); # int(0.5 + 100 * $maln/$mlen);
	my $pctident= pctOf( $mident/$maln, 10); # int(0.5 + 1000 * $mident/$maln)/10; # 
	
## change to match gmap.gff attrib??
## aaln=479;cov=17.8;indels=43/0;nexon=10;pid=94.5;qlen=8340;
## qlen == clen; pid*coverage == malign/mlen ;  indels?  nexon?
## change exon score to percent from .proportion?
	
	##? genescore: multiply proportions for all score parts; max = 1; BUT not .90 aln x 0 splice .. need fudge
	my $genescore = $pctalign/100; #? add in pctident?
	
	## fix for no splice at tr ends: +4 splice score: $msplscore+=4 if($msplscore>0)
	## fix for 1-exon spct= (2 + 4)/(5*2) = .60 ?? ;  (3 + 4)/(5*3) = 0.47
	
	use constant ONEX_SPLSCORE => 0.70; ## .50;
	if($mnxon>1) { my $spct= ($mnxon+$msplscore+4)/(5*$mnxon); $genescore *= $spct; } 
	else { $genescore *= ONEX_SPLSCORE; } # single-exon adjust? not good enough for dup choice w/ other has 2+


	## ^^ min spct for 0 splice is $nx / 5*nx .. 2/10, 3/15, 4/20 = .20;   
	## eg. hiscore 19 exons, 68 splscore, spct= .92 
	## gapspan included in pctalign;  
	# * $tscore? dont know what splign.log-tscore represents.. leave out?
	# uck: gescore=1.01084210526316;
	$genescore= ($genescore>1)? 100 : pctOf($genescore); # int(0.5 + 100*$genescore);
	
	## gmap.gff mrna attr: Target=xxx;aalen=220;cov=100.0;match=678;nexon=3;pid=99.9;qlen=679;..
	
	## 1707: cov=$pctalign%,$maln/$mlen << drop?? NO KEEP, used below ,maln/mlen excess
	my $mattr="ID=$id;cov=$pctalign%,$maln/$mlen;pid=$pctident;nexon=$mnxon;splice=$msplscore$addat"; # $pctalign%, ?? dont need pctalign 2 places.. keep in col5
	$mattr.= ";Target=$trid $mrnab $mrnae"; # always add, so can measure partial aligns of mrna
	# $mattr.= ";trid=$trid" if($id ne $trid); # or Target=$trid 1 $mlen ?? add always??
	# $mattr.= ";splscore=$msplscore" if($mnxon>1); #? or always
	$mattr.= ";indels=$mindels" if($mindels); # equiv gmap.gff
	#  gaps: my($gapt)=m/(\d+) bases in cdna/; == Inserts; my($gapg)=m/(\d+) bases in genome/; == Dels
  #  gmap.gff "indels=$gapt/$gapg" if($gapt>0 or $gapg>0);  

	$mattr.= ";gaps=$gapspan,$gaps" if($gaps);	## prefix with gapspan and/or ngaps
	$mattr.= ";gescore=$genescore"; # insert in col5 instead of pctalign?

	# $mattr.= ";tscore=$tscore" if(1); # if($KEEP_SPLIGNSCORE);
	$mattr.= ";splinfo=$idir,nd$tsdup,sc$tscore" if($debug); # if($KEEP_SPLIGNSCORE);
	$mattr.= ";$trattr" if($trattr);
	
	## ?? replace print here with store in @gff, then filter out low-score 2ndary aligns per trid
	## use genescore AND mrnaspan to decide to keep, need mrnaspan parts for split-genes
	my $genomespan="$mref:$mgb-$mge:$mor"; # mor ?
	my $mrnaspan="$trid:$mrnab-$mrnae"; #?
	
	my $mrnagff= join("\t",$mref,$SRC,"mRNA",$mgb,$mge,$pctalign,$mor,".",$mattr)."\n";
	# print $OUTH join("\t",$mref,$SRC,"mRNA",$mgb,$mge,$pctalign,$mor,".",$mattr)."\n";
	$ntrout++;
	
	my @exongff=();
	foreach my $x (@xon) {
	  $x->[8] =~ s/Parent=$trid/Parent=$id/ if($id ne $trid);
	  $x->[8] =~ s/;align=.*$// unless($debug); # dont need output? or leave to drop later?
	  $x->[6] = $mor if($flipor);
	  #causesBUGinSplitparts# $x->[5] = pctOf( $x->[5]); # int(0.5 + 100*$x->[5]); # .prop to pct%
	  my $xout= join("\t",@$x)."\n"; 
	  push @exongff, $xout;
	  # print $OUTH join("\t",@$x)."\n";
	}
	
	# OPTION: output introns where adjacent exons have proper splices ..
	
	# FIXMEd: add CDS exons, given mRNA CDS offsets .. need more inputs
        # FIXME5: 16.10.11 : sense=-1; flipor CDS can be bad, wrong calc from exons, offset

	my @cdsgff=();
	if($trattr =~ /offs=(\d+).(\d+)/) { # add cds-exons : option?  check $trattr=~/aalen=/ also ?
	  my($cb,$ce)=($1,$2); 
	  my $cdsdir= 1; # ($mor eq '-')?-1:1; # not this now, makeCDSexons() finds orient
	  #d if($ce < $cb) { ($cb,$ce)= ($ce,$cb); $cdsdir= -$cdsdir; } #?
	  @cdsgff= makeCDSexons($cb,$ce,$cdsdir,\@xon); ## add: $flipor -sense flag
	}
	
	my @gene=($mrnagff,@exongff,@cdsgff);
	if($debug and @xonpoor) { push @gene, @xonpoor; }
	return ($trid,$genescore,$genomespan,$mrnaspan,\@gene); 
}

=item genescore problems
	
	.. problem deciding what part to weight more : align cover is most important, but
		 more valid splices w/ slightly lower cover is better (??)
		 
## bad genescore, splice=0 above splice=8 .. nxon=1 for bad choice, nxon=3 for good miss.
# Scaffold1810    splkf2eg7m1     mRNA    22562   23156   88      -       .       ID=Funhe2Eq7m000962t1;cov=88%,566/642;nex
#   on=1;splice=0;Target=Funhe2Eq7m000962t1 49 642;gaps=48,LGap:1-48,;gescore=79;splinfo=+8643,nd9,sc44.929;aalen=164,76%,par
#   tial3;clen=642;offs=149-640;oid=kfish2qf7soapk21loc1010032t1
# #d1.Scaffold474 splkf2eg7m1     mRNA    170935  187263  89      +       .       ID=Funhe2Eq7m000962t1;cov=89%,569/642;nex
#   on=3;splice=8;Target=Funhe2Eq7m000962t1 2 642;gaps=1,LGap:1-1,;gescore=65;splinfo=+1647,nd9,sc54.389;aalen=164,76%,partia
#   l3;clen=642;offs=149-640;oid=kfish2qf7soapk21loc1010032t1
# #d2.Scaffold271 splkf2eg7m1     mRNA    267767  268201  62      -       .       ID=Funhe2Eq7m000962t1;cov=62%,395/642;nex
#   on=1;splice=2;Target=Funhe2Eq7m000962t1 207 642;gaps=206,LGap:1-206,;gescore=56;splinfo=+7228,nd9,sc123.48;aalen=164,76%,
#   partial3;clen=642;offs=149-640;oid=kfish2qf7soapk21loc1010032t1

## more genescore quandry : is spliced low cov better/worse than unspliced hicov?
## >> all spliced, but 33% unaligned
# Scaffold264     splkf2eg7m1     mRNA    363123  388491  67      -       .       ID=Funhe2Eq7m000197t1;cov=67%,369/547;nex
#  on=4;splice=12;sense=-1;Target=Funhe2Eq7m000197t1 93 545;gaps=147,LGap:546-547,MGap:129-181,RGap:1-92,;gescore=67;splinfo
#  =-1116,nd4,sc83.058;aalen=127,69%,partial3;clen=547;offs=167-547;oid=kfish2qf7soapk21loc1001798t1
# #d1.Scaffold264 splkf2eg7m1     mRNA    358922  388491  70      +       .       ID=Funhe2Eq7m000197t1;cov=70%,385/547;nex
#  on=5;splice=10;Target=Funhe2Eq7m000197t1 5 545;gaps=125,RGap:546-547,MGap:182-202,LGap:1-4,MGap:127-168,MGap:34-89,;gesco
#  re=53;splinfo=+1116,nd4,sc100.074;aalen=127,69%,partial3;clen=547;offs=167-547;oid=kfish2qf7soapk21loc1001798t1
## >> 1 exon 98% aligned << GMAP picks this one; is  Transposon gene xpressed.
# #d2.Scaffold827 splkf2eg7m1     mRNA    75515   76058   98      -       .       ID=Funhe2Eq7m000197t1;cov=98%,534/547;nex
#  on=1;splice=0;Target=Funhe2Eq7m000197t1 1 546;gaps=1,RGap:547-547,;gescore=49;splinfo=+8227,nd4,sc19.502;aalen=127,69%,pa
#  rtial3;clen=547;offs=167-547;oid=kfish2qf7soapk21loc1001798t1

	# ** BAD here? 1-exon lower cover scores above many-exon,higher cover model.. due to missing canon splices
	# .. maybe not bad, depends how much to weight canon-splices. 5-exon part align could be trash
	## >> bad case, 400 bp fragment likely TE-related, many valid aligns.. pick any
	# Scaffold2012:1096-16540:- FunhetEGm000040t3_G9;cov=58%,228/392;nexon=5;splice=12;Target=1 354;gaps=139;gescore=39;splinfo=-20302,nd9,sc120.736
  # Scaffold9953:1773835-1774051:- FunhetEGm000040t3_G6;cov=49%,193/392;nexon=1;splice=0;Target=1 213;gaps=179;gescore=44;splinfo=+24504,nd9,sc116.334
  # gmap got: Scaffold9941:1139485-1139880, 100% align, 1 exon
	## another problem case; gmap Scaffold10169:368500-368925:., 1exon
	# Scaffold151     11294   11509   47      +       ID=FunhetEGm000261t2;cov=47%,200/423;nexon=1;splice=4;Target=Fu
	# nhetEGm000261t2 147 395;gaps=174,RGap:1-146,LGap:396-423,;gescore=42;splinfo=-1038,nd9,sc121.769
	# #d1.Scaffold10107       555234  736229  58      -       ID=FunhetEGm000261t2;cov=58%,244/423;nexon=6;splice=6;T
	# arget=FunhetEGm000261t2 69 340;gaps=151,RGap:1-68,LGap:341-423,;gescore=23;splinfo=-26355,nd9,sc133.703
	# #d2.Scaffold10057       119821  1114429 72      +       ID=FunhetEGm000261t2;cov=72%,306/423;nexon=7;splice=2;T
	# arget=FunhetEGm000261t2 1 328;gaps=95,RGap:329-423,;gescore=19;splinfo=+12384,nd9,sc89.493
	## another .. gmap at Scaffold609:125478-125752:-, this is tiny frag. 
	# Scaffold609     124944  125176  52      +       ID=FunhetEGm000532t5;cov=52%,193/373;nexon=2;splice=4;Target=Fu
	# nhetEGm000532t5 1 203;gaps=170,LGap:204-373,;gescore=31;splinfo=-3552,nd3,sc103.482
	# #d1.Scaffold609 105775  125841  66      -       ID=FunhetEGm000532t5;cov=66%,248/373;nexon=7;splice=8;Target=Fu
	# nhetEGm000532t5 1 321;gaps=107,MGap:207-221,MGap:252-291,RGap:322-373,;gescore=28;splinfo=+17201,nd3,sc103.174;
	# #d2.Scaffold609 105775  125841  56      -       ID=FunhetEGm000532t5;cov=56%,209/373;nexon=5;splice=6;Target=Fu
	# nhetEGm000532t5 1 321;gaps=158,MGap:269-300,LGap:322-373,MGap:183-256,;gescore=25;splinfo=-17201,nd3,sc95.246;a

=cut

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

sub makeCDSexons {
  my($cb,$ce,$cdsor,$xons)= @_; ## add $flipor -sense flag
  ## FIXME5: 16.10.11 : sense=-1; flipor CDS can be bad, wrong calc from exons, offset
  ## .. likely Target span needs flip for -sense

  ##  BAD cds-end clips wrong side for some: sense=-1 and strand=-, as one case
  ##  ** cdsor MAY be wrong here, need Target/mRNA orient here
	##  $cdsor= ($mor eq '-')?-1:1;
  return () unless($cb and $ce and ref($xons));
  
  use constant NEWCDSX => 1;
	my @cdsgff=();
	my ($cdslen, $phase, $cdsfirst, $targor, $xrevorder)=(0) x 9; 
  $cdsfirst= 1;  
  my $nexon= @$xons;
  if($ce > 0 and $ce < $cb) { 
    ($cb,$ce)= ($ce,$cb); $cdsor= -$cdsor; #?
  }
  
  my(@xa,@xb);
  for(my $i=0; $i<2; $i++) { # foreach my $x (@{$xons}[0,1]) 
    last if($i>=$nexon); 
    my @c;  my $x= $xons->[$i];
    if(ref($x) =~ /ARRAY/) { @c= @$x; }
    elsif($x =~ /\texon\t/) { @c= split"\t",$x; chomp($c[-1]); }	
    else { next; }    
    if(@xa) { @xb=@c; last; } else { @xa=@c; }
  }  
  @xb=@xa if($nexon==1);
  $xrevorder= ($xa[3] > $xb[3]) ? -1 : 0; # exon sort order, fwd or rev?
  my($tb1)= $xa[8] =~ m/(?:Target|trg)=\S+\s(\d+)/;
  my($tb2)= $xb[8] =~ m/(?:Target|trg)=\S+\s(\d+)/;
	if($xrevorder) { $targor=($tb2 and $tb1 > $tb2)? 1 : ($tb1 and $tb1 < $tb2) ? -1 : 0; }
  else { $targor=($tb2 and $tb1 > $tb2)? -1 : ($tb1 and $tb1 < $tb2) ? 1 : 0; }

	my $chrstrand= # return this? 
	  ($targor > 0 and $cdsor >= 0) ? '+' :
	  ($targor > 0 and $cdsor  < 0) ? '-' :  
	  ($targor < 0 and $cdsor  < 0) ? '+' :  # CDSb == chr-start
	  ($targor < 0 and $cdsor >= 0) ? '-' :  # CDSb == chr-end
	  ($targor < 0)?'-':'+'; # or '.'; # ($chror < 0)?'-':'+' ; # cant tell, use input

  foreach my $x (@$xons) {
    my @c; # my @c= @$x;
    if(ref($x) =~ /ARRAY/) { @c= @$x; }
    elsif($x =~ /\texon\t/) { @c= split"\t",$x; chomp($c[-1]); }	    
    else { next; }    
    my($tb,$te)= $c[8] =~ m/Target=\S+ (\d+) (\d+)/;
    next unless($te); # warn/error

    if($te > $cb and $tb < $ce) {
if(NEWCDSX) {
      if($tb <= $cb) { my $d=$cb-$tb;  ## chror < 0 implies chre == tb, col[4]
        if($xrevorder) { if($cdsfirst) { $c[4] -= $d; } else { $c[3] += $d; } }
        else { if($cdsfirst) { $c[3] += $d; } else { $c[4] -= $d; } }
        $cdsfirst= 0;
        }
      if($te >= $ce) { my $d=$te-$ce; ## chror < 0 implies chrb == te, col[3]
        if($xrevorder) { if($cdsfirst) { $c[4] -= $d; } else { $c[3] += $d; } }
        else { if($cdsfirst) { $c[3] += $d; } else { $c[4] -= $d; } }
        $cdsfirst= 0;
        }
} else { # old method, bad
      ## fixme again; -strand needs c3,c4 swap also
      if($tb < $cb) { my $d=$cb-$tb; 
        if($cdsor < 0) { $c[4] -= $d; } else { $c[3] += $d; } # FIXME: - bad or offset; flip d sign, $or * $d
        }
      if($te > $ce) { my $d=$te-$ce; 
        if($cdsor < 0) { $c[3] += $d; } else { $c[4] -= $d; }
        }
}

      ($phase,$cdslen)= phaseCDS($tb,$te,$cb,$ce,$cdsor,$cdslen);

      $c[2]= "CDS";
      # $c[5]=1; #? score; leave exon score??
      $c[7]= $phase; ## FIXME: need phase calcs
      my $CUT_XATTR='Target|align|splice';
      $c[8] =~ s/;($CUT_XATTR)=[^;\n]+//g; # leave some xattr: ix=, indel= 
      # $c[8] =~ s/;Target=.*//; # leave some xattr: ix=, indel= 
      my $xout= join("\t",@c)."\n"; 
      push @cdsgff, $xout;
    }
  }
	 
	return @cdsgff;
}

sub pctOf { 
  my($prop,$decimalfac)=@_;
  $decimalfac ||= 1; # 10,100,1000 .. for 1,2,3 decimals
  my $p= ($prop <= 1.0) ? int(0.45 + $decimalfac * 100 * $prop)/$decimalfac : $prop;
  $p=100 if($p > 100);
  return $p;
}

=item trimGenePart 
  my ($newstart,$newend,$newpart)= trimGenePart( "start", $mb, $me1+1, $geneset[$ig]);
  * need to change all types, exon, mrna, cds? target+genome starts
  input=output of putgff: ($trid,$genescore,$genomespan,$mrnaspan,\@gene);
  			  # my($gidi,$gsci,$gspi,$mspi,$gffi)=@{$geneset[$ig]};
=cut

sub exonRanges {
  my($xon)=@_;
  my($gref,$ggb,$gge,$trid,$ttb,$tte)=(0) x 9;
  for my $x (@$xon) {
    my @x = split"\t",$x;
    my($gr,$gb,$ge)= @x[0,3,4];
    my($td,$tb,$te)= $x[8] =~ m/Target=(\S+)\s(\d+)\s(\d+)/;
    $trid=$td; $gref=$gr;
    $ggb= $gb if($ggb==0 or $gb<$ggb);
    $gge= $ge if($gge==0 or $ge>$gge);
    $ttb= $tb if($ttb==0 or $tb<$ttb);
    $tte= $te if($tte==0 or $te>$tte);
  }
  return($gref,$ggb,$gge,$trid,$ttb,$tte);
}


sub relocateExon {
  my($startend, $xon, $xdif)= @_;
  my @x=split"\t",$xon;  # if(ref($xon) =~ /ARRAY/)
  my $xoff= 1+$xdif; my $isrev=($x[6] eq '-')?1:0;
  my($td,$tb,$te)= $x[8] =~ m/Target=(\S+)\s(\d+)\s(\d+)/; 
  if($startend eq "end") {   
    my $ted = $te - $xoff;
    if($ted >= $tb + MIN_EXON_SPAN) {
      $x[8] =~ s/Target=$td $tb $te/Target=$td $tb $ted/;
      #notsure# $x[4] -= $xoff;  # shift down FIXME??? revor change gstart not gend
      if($isrev) { $x[3] += $xoff; } else { $x[4] -= $xoff; } # shift down FIXME: revor change gstart not gend
      return join"\t",@x; 
    }
  } elsif($startend eq "start") {
    my $tbd = $tb + $xoff;
    if($tbd <= $te - MIN_EXON_SPAN) {
      $x[8] =~ s/Target=$td $tb $te/Target=$td $tbd $te/;
      #notsure# $x[3] += $xoff;  # shift up. FIXME?? revor change gend not gstart
      if($isrev) { $x[4] -= $xoff; } else { $x[3] += $xoff; } # shift up. FIXME: revor change gend not gstart
      return join"\t",@x; 
    }
  }
  return "";
}

sub relocateGenePart {
  my($gidi,$gsci,$gffi,$xkeep)= @_;
  
  return(0,0,undef) unless(@$xkeep > 0);
  my($mrna)= grep /\tmRNA/, @$gffi;
  my @mrna = split"\t",$mrna;
  # my($gref,$src,$xtype,$gb,$ge,$pident,$strand,$phx,$attr)= @mrna;

  my($gref,$ggb,$gge,$trid,$ttb,$tte)= exonRanges($xkeep); # genome and mrna target

  ## FIXME here? mrna at: cov=112% for sum of splits, too high, recalc after trims.
  ## see above: my $mattr="ID=$id;cov=$pctalign%,$maln/$mlen;pid=$pctident;nexon=$mnxon;splice=$msplscore$addat";
  # Funhe2EKm000191t1_C1 cov=67%,1190/1775 : ok, untrimmed.
  # Funhe2EKm000191t1_C2 cov=46%,822/1775; : too high, trimmed to trspan=1204 1775 = 572, cov=32%
  my($pcov,$maln,$mlen)= $mrna[8]=~m/cov=(\d+)%,(\d+).(\d+)/;
  my $maln2= 1 + $tte - $ttb;
  if($maln and $maln2 != $maln) {
    return(0,0,undef) if($maln2 < $MINPARTSPAN);
    my $pcta2= pctOf($maln2/$mlen); # int(0.5 + 100 * $maln2/$mlen);    
    $mrna[8]=~s/cov=$pcov%,$maln/cov=$pcta2%,$maln2/;
  }
   
  my($cb,$ce)= $mrna =~ /\boffs=(\d+).(\d+)/;
  my $cdsdir= 1; # ($mrna[6] eq '-')?-1:1; # not this now, makeCDSexons() finds orient
	# if($ce < $cb) { ($cb,$ce)= ($ce,$cb); $cdsdir= -$cdsdir; } #?
  my @cdsnew= makeCDSexons($cb,$ce,$cdsdir,$xkeep);

  my $genomespan="$gref:$ggb-$gge:$mrna[6]"; # mor ?
  my $mrnaspan="$trid:$ttb-$tte"; #?
  $mrna[3]= $ggb; $mrna[4]= $gge; $mrna[8]=~s/(Target=\S+)\s(\d+)\s(\d+)/$1 $ttb $tte/;
  $mrna=join"\t",@mrna;
  my @gene=($mrna,@$xkeep,@cdsnew);
  my $newpart= [ $gidi,$gsci,$genomespan,$mrnaspan, \@gene];
  my $newstart= $ttb; my $newend= $tte;
  return($newstart,$newend,$newpart);
}

sub trimGenePart {
  my( $startend, $atmrnaloc, $newmrnaloc, $genepart, $lastpart)= @_;
  my($gidi,$gsci,$gspi,$mspi,$gffi)=@$genepart;
  my($gidl,$gscl,$gspl,$mspl,$gffl)=(ref $lastpart) ? @$lastpart : ("",0,"","",[]);
  my($md,$mb,$me)= split/[:-]/, $mspi;  # atmrnaloc == mb or me
  my($mdl,$mbl,$mel)= split/[:-]/, $mspl;  # atmrnaloc == mb or me
  
  my $newlast=undef;
  my($newstart,$newend,$newpart)= ($mb,$me,$genepart); # default no changes
  $newstart=$newmrnaloc if($startend eq "start");
  $newend=  $newmrnaloc if($startend eq "end");
  
  # my($mrna)= grep /\tmRNA/, @$gffi;
  my @xon= grep /\texon\t/, @$gffi;
  my @lxon= grep /\texon\t/, @$gffl; #? assume sorted by target.begin ?

if(1) {  
  my (@xdrop,@ldrop);
  my @xstarts= map{ my($td,$tb,$te)= m/Target=(\S+)\s(\d+)\s(\d+)/; [$tb,$te];} @xon;
  my @lstarts= map{ my($td,$tb,$te)= m/Target=(\S+)\s(\d+)\s(\d+)/; [$tb,$te];} @lxon;
  my ($ixdif,$jxdif,$xdif)=(0) x 4;
  
  ## TooManyParts BUG: MGap let many parts insert into gaps for low ident but many-dup cases
  ## ID=Funhe2EKm023981t1_C1;Split=1/210;Name=Integrin alpha-M
  ## C1: trg=Funhe2EKm023981t1 1 3060;gaps=3683,MGap:419-1521,MGap:2848-2981,MGap:1548-1780,MGap:1804-2087,MGap:3002-3031,MGap:3061-3457,RGap:3489-4174,MGap:2105-2835,MGap:29-113,
  ##  .. mrna span large, but have only end exons: t1 1 28, <BIG MGap> t1 1522 1547, t1 2088 2104, ..t1 3032 3060
  ## C2: trg=Funhe2EKm023981t1 2 418; C3: trg=Funhe2EKm023981t1 4 1039; C4: trg=Funhe2EKm023981t1 4 420
  ## .. require no mRNA target part overlap ? ie. drop this xstarts/lstarts in favor of mRNA mb,me part overlaps
  # Fix1..BUT, still have parts w/ mostly overlapping mrna target spans, exon overlap filtering not right?
  # eg: Funhe2EKm017371t1_C1: t1 1 2349; C2: t1 74 459; C3: t1 127 1011; C4: t1 130 1290
  # Maybe need all prior part exons in @lstarts?
  
  my($ldropj0,$ldropj1)=(-1,-1);
  for(my $ix=0; $ix<@xstarts; $ix++) {
    my($tb,$te)= @{$xstarts[$ix]};
    my($xident)= (split"\t",$xon[$ix])[5]; #  hack..  
    my $lidsum= 0; my @lov=();
    for(my $jx=0; $jx<@lstarts; $jx++) {
      next if($jx <= $ldropj1 and $jx >= $ldropj0);
      my($lb,$le)= @{$lstarts[$jx]};
      if($tb < $le and $te > $lb) {
        my($lident)= (split"\t",$lxon[$jx])[5]; #  hack..  
        $lidsum = $lident; #last only? if($lident>$lidsum); # $lidsum*$lident; 
        push @lov,$jx;
        $xdif= $le - $tb; $ixdif= $ix; $jxdif=$jx;
      }
    }
    #? change here for TooManyParts? require @ldrop or @xdrop to be contiguous set of exons
    # .. i.e. if drops = 2,5,9 of 11 change to 2..9; 
    # .. Also? require drop from one end, eg. 1..9 of 11, or 2..11 of 11
    # DEBUG swap# if($xident < $lidsum) { push @ldrop,@lov; } elsif($lidsum>0) { push @xdrop, $ix; } # DEBUG
    if($xident > $lidsum) { 
      if(@lov) { push @ldrop,@lov; @ldrop= sort{$a <=> $b} @ldrop; ($ldropj0,$ldropj1)= @ldrop[0,-1]; }
    } else { push @xdrop, $ix; }
  }
  
  if(@xdrop) { my @xkeep=();
    @xdrop= sort{$a <=> $b} @xdrop; my($ix0,$ix1)= @xdrop[0,-1];
    for(my $ix=0; $ix<@xon; $ix++) { 
      push @xkeep, $xon[$ix] unless($ix >= $ix0 and $ix <= $ix1); # TooManyParts fix
      #old# push @xkeep, $xon[$ix] unless(grep{$_ == $ix} @xdrop); 
      }
    if($xdif > 1 and grep{$_==$ixdif}@xdrop) { 
      if(my $xcut= relocateExon("start",$xon[$ixdif], $xdif )) { unshift @xkeep,$xcut; $xdif=0; }
      #   my $xcut= $xon[$ixdif]; my @xcut=split"\t",$xcut; 
      #   my($td,$tb,$te)= $xcut[8]=~m/Target=(\S+)\s(\d+)\s(\d+)/; 
      #   my $tbd = $tb+1+$xdif;
      #   if($tbd < $te-2) {
      #   $xcut[8]=~s/Target=$td $tb $te/Target=$td $tbd $te/;
      #   $xcut[3] += 1+$xdif; $xcut=join"\t",@xcut; # shift up.
      #   unshift @xkeep,$xcut; $xdif=0;
      #   }
    }
    ($newstart,$newend,$newpart)= relocateGenePart($gidi,$gsci,$gffi,\@xkeep);
  }
  
  if(@ldrop) { my @lkeep=();
    @ldrop= sort{$a <=> $b} @ldrop; my($ix0,$ix1)= @ldrop[0,-1];
    for(my $ix=0; $ix<@lxon; $ix++) { 
      push @lkeep, $lxon[$ix] unless($ix >= $ix0 and $ix <= $ix1); # TooManyParts fix
      #old# push @lkeep, $lxon[$ix] unless(grep{$_ == $ix} @ldrop);
      }
    if($xdif > 1 and grep{$_==$jxdif}@ldrop) { 
      if(my $xcut= relocateExon("end", $lxon[$jxdif], $xdif )) { push @lkeep,$xcut; $xdif=0; }
      #   my $xcut= $lxon[$jxdif]; my @xcut=split"\t",$xcut; 
      #   my($td,$tb,$te)= $xcut[8]=~m/Target=(\S+)\s(\d+)\s(\d+)/; 
      #   my $ted = $te-1-$xdif;
      #   if($ted > $tb+2) {
      #   $xcut[8]=~s/Target=$td $tb $te/Target=$td $tb $ted/;
      #   $xcut[4] -= (1+$xdif); $xcut=join"\t",@xcut; # shift down.
      #   push @lkeep,$xcut; $xdif=0;
      #   }
    }
    if(@lkeep) {
      my($lstart,$lend,$lnewpart)= relocateGenePart($gidl,$gscl,$gffl,\@lkeep);
      $newlast= $lnewpart;  # $lastpart update to $lnewpart;  
    } else {
      $newlast= [ $gidl,0,0,0, undef]; # remove this last part
    }
  }
  return($newstart,$newend,$newpart, $newlast);

} else {  
  my($lastident)= (split"\t",$lxon[-1])[5]; #  hack..  
  my @xkeep=(); my $droplast=0;
  for my $x (@xon) {
    my($td,$tb,$te)=  $x =~ m/Target=(\S+)\s(\d+)\s(\d+)/;
    if($tb < $mel and $te > $mbl) {
      my @xv=split"\t",$x; 
      my $xident= $xv[5]; #pident, check vs lxon.last.pident      
      if($xident > $lastident) { push @xkeep,$x; $droplast++; }
    } else {
      push @xkeep, $x;
    }
  }
  
  if(@xkeep < @xon) { # relocate $newpart # update: $gspi,$mspi,$gffi
    ($newstart,$newend,$newpart)= relocateGenePart($gidi,$gsci,$gffi,\@xkeep);
  }
  if($droplast) { # relocate lastpart
    my($lastx)= pop(@lxon);
    my($lstart,$lend,$lnewpart)= relocateGenePart($gidl,$gscl,$gffl,\@lxon);
    $newlast= $lnewpart; # $lastpart update to $lnewpart;  
 }
  return($newstart,$newend,$newpart, $newlast);
}  
}


sub trimGenePartOLD {
  my( $startend, $atmrnaloc, $newmrnaloc, $genepart, $lastpart)= @_;
  my($gidi,$gsci,$gspi,$mspi,$gffi)=@$genepart;
  my($md,$mb,$me)= split/[:-]/, $mspi;  # atmrnaloc == mb or me
  
  my($newstart,$newend,$newpart)= ($mb,$me,$genepart); # default no changes
  $newstart=$newmrnaloc if($startend eq "start");
  $newend=  $newmrnaloc if($startend eq "end");
  
  my($mrna)= grep /\tmRNA/, @$gffi;
  my @xon= grep /\texon\t/, @$gffi;
  my($td,$tb,$te);
  ($td,$tb,$te)=  $xon[0] =~ m/Target=(\S+)\s(\d+)\s(\d+)/;
  if($tb == $atmrnaloc) { # and $startend eq "start"
  
  }
  ($td,$tb,$te)=  $xon[-1] =~ m/Target=(\S+)\s(\d+)\s(\d+)/;
  if($te == $atmrnaloc) { # and $startend eq "end"
    
  }
  
  return($newstart,$newend,$newpart);
}



__END__

=item splign vs gmap

  info for kfish2p67vs.mrna gene set
  	.. overall pctalign is higher for splign
  	.. gmap for this mrna set has high proportion of split genes (chimera),
  		 but they are split at same locus, ie have gaps or mismatches to genome,
  		 and full span should be alignable .. maybe other gmap option will do that, dont know ..
  		 
  gmap nopath:14950  splign nomap:30106 of ntr=137660 ids
  .. check 'good' gmap chimera on 2 scaffolds versus splign
	
	split-gene FunhetEGm067171t5 failure in compart:
		- need to test other opts to get such, or replace w/ blastn  > splign..
		- this opt should pull many low-cover cases: -min_compartment_idty 0.50 => 0.20
					 
					 
	eg. spl-nomap hiqual
		FunhetEGm028320t1 4274,91%,complete	# gm: Scaffold10072:6500-8633,cov=15%  
			^^ end of scaf problem, split also?
		  ^^ there is more rna-seq at scaf:1-6500 should be gene asm parts; geno-align problem?
			^^ good full-length fish gene:  sacsin-like, 80% aln to 4284 aa, cichlids, tetraodon,xenopus,takifugu,..
		FunhetEGm067171t5 3577,62%,complete # gm split: Scaffold10067:190888-209862,cov=36%,3'part2 + Scaffold9902:670788-785051,cov=64,5'part
				# ^ this is true split-gene; note t1-4,t6 are mapped and longer **
		FunhetEGm025643t1 3222,61%,complete	# gm: Scaffold755:62292-137840,cov=52%
			^^ spans 75kb, with 2 smallish intronic genome gaps; other fish prots map at 52% here.
			>> there is 30kb gap at 5'end, start of scaf; this gene covers all of 150kb Scaffold755 (but top end, rev gene), 
			>> introns span gene.
			-- what part of tr/prot is not aligned? aug models are 1/2 size 1400 aa but same mapped exons.
			^^ good full-length fish gene: SZT2, 80% aln to 3698 aa, cichlids, tilapia,  ..
			>> Maylandia zebra, african cichlid is turning up top-hit for these, other kfish genes.
			# ^ these 3 not in compart output table, no align to genome?
		
		FunhetEGm018322t1 2941,93%,complete	# gm split: Scaffold630:10474-27099,cov=30,5'end + Scaffold848:4429-30540,cov=70,3'end
		FunhetEGm073739t1 2819,96%,complete # gm: Scaffold10038:1175876-1194030,cov=30
			^^ blastp align fragmented but full span for Maylandia, Oryzias fish top hits 2313aa, 40% cov/40% ident,  CDD:MDN1 AAA ATPase containing von Willebrand factor type A
			-- genome asm several big/small gaps both ends of mapped exons; missing 70% likely in 20kb gap.
			splign several alts map: FunhetEGm073739t13 cov=82%, 2/3 size; FunhetEGm073739t4 cov=75%, 2k-inc. 
						FunhetEGm073739t2 cov=88%, 2780 aa but exons spred beyond expected gene span.
		FunhetEGm058834t1 2741,90%,complete # gm: Scaffold10038:1365317-1389369,cov=46
			^^ blastp fragmented but full span for  Maylandia, tilapia fish, 3790aa, serine-rich adhesin for platelets-like
			splign 3 alts map here: FunhetEGm058834t8, 96% 1650aa-inc; FunhetEGm058834t10 ditto.. short aa
			.. problem alts: FunhetEGm058834t6, FunhetEGm058834t14 map past 5'end of t1..8: are these paralogs? evidence unclear
				.. prots are repetitive, but FunhetEGm058834t1 x FunhetEGm058834t6 only slight align.
			-- no genome gaps? ah, a large gap 30kb past 3' mapped end; 
			-- no other fish prots mapped here.. 
			
		FunhetEGm018012t1 2693,83%,complete	# gm: cov=47% Scaffold521:144238-196402 
			- several alts gm maps, spl misses; possible tandem dup gene, or one w/ exon mixup,
			- no compart align for FunhetEGm018012t*, but gmap got ~50% of many.
			- 52kb span, 1 small intronic scaf gap, but Scaffold521 ends at 3' exon .. end of scaf problem
			= 'Telomerase-associated protein 1' Funhe5EG032579t1 (1193aa), Funhe5EG032580t1(460 aa)
		FunhetEGm022232t1 2625,96%,partial3	# gm: cov=32% Scaffold9898:106970-129598 
			- strong express 5' end with messed up intron mapping: many overlaping, bi-dir.
			- Neuroblast differentiation-associated protein AHNAK,
			- splign maps FunhetEGm022232t2, 1180aa-inc, 100% cov; 
			- confused locus: several other trloci/diff prots, map here also; 3kb gap past 3' intron mess
			- expression suggests 2-3 long exons of 2k-5k, but introns map into these express exon spans.
		FunhetEGm009379t2 2605,73%,partial3	# gm: cov=32%
		
			* retest compart,splign on these cases of poor map (due to genogaps?), not splitscaf:
				FunhetEGm025643t1/Scaffold755, 
				FunhetEGm073739t1/Scaffold10038, FunhetEGm058834t1/Scaffold10038, 
				FunhetEGm018012t1/Scaffold521
					
			* blastp check proteins to see which of these are "real" full kfish genes?
			'FunhetEGm028320t1|FunhetEGm067171t5|FunhetEGm018322t1|FunhetEGm025643t1|FunhetEGm073739t1|FunhetEGm058834t1'
			
			* blastp check cases of gappy align to genome: is problem commoner w/ mrna or genome?
			
			
		FunhetEGm016085t4 2587,85%,complete	# t1 is mapped, very short 101aa, others t2-5 are 1500-2500aa
			# ^ what gives w/ evg tr2aacds calling short,partial 100aa FunhetEGm016085t1=kfish2qf7soapk21loc945695t1
			# t4 = kfish2qf7velvk71Loc2507t16
			# this is a NOMAIN, algo should have picked longest aa as t1: trevg67pt/publicset/kfish2p67vs.mainalt.tab
			# kfish2eg6velvk55Loc4503t9 NOMAIN  kfish2qf7soapk21loc945695t1/althi,kfish2eg6velvk55Loc4503t1/althi,
			#     kfish2eg6velvk55Loc4503t4/althi,kfish2qf7velvk71Loc2507t16/althi,kfish2eg6velvk61Loc4732t4/althi
 
=item retest missed compart

# test compart new opts for missed aligns:
cpdefault="-penalty 0.55 -min_idty 0.70"
cpopts="-penalty 0.25 -min_idty 0.25";  
test2: "-penalty 0.25" << no align
test3: "-min_idty 0.25"  << full align as 1st; use this min_idty (or lower?)
test4: "-min_idty 0.15" << same as .25
test5: "-min_idty 0.45" << same as .25

try this splign opt to match compart opt:
	splopt="-type mrna -min_compartment_idty 0.45"


CASE5: split FunhetEGm018322t1/s08split = Scaffold630,Scaffold848  
  compart cpopt=-min_idty 0.45
    >> only Scaffold848:4427-30540:+, tr:2832-9451
    
  compart cpopt=-penalty 0.25 -min_idty 0.25
    ++ adds Scaffold630:10474-27099:-, tr:1-2833
  
  compart cpopt=-penalty 0.25 -min_idty 0.45
    >> NO, same as 1.
    
  compart cpopt= -min_idty 0.25  
    ++ both, same as #2,  ie. skip -penalty, use low -min_idty
    
	** need new splign opts for this split case?  2nd part is lower score..
	a: splopt="-min_compartment_idty 0.45"
	b: splopt="-min_compartment_idty 0.25"  << this one to get 2nd part of split **
	   what of -compartment_penalty 0.55 def; No
	c: splopt="-min_compartment_idty 0.45 -compartment_penalty 0.25"   << no good
		"Multiple compartments will only be identified if they have at least this level of coverage"
			-- is that cover w/in the compart, or total cover of the mrna ?
	  -min_compartment_idty <Real, 0..1> Minimal compartment identity to align. Default = `0.7'

	
	$nbin/splign $splopt -comps $qfile1.compart8 -blastdb $genome -ldsdir ldsdir \
	  -log $qfile1.splog8a -aln $qfile1.aln > $qfile1.splign8a
	#.....
		FunhetEGm018322t1.splog8a
		+1      FunhetEGm018322t1       Scaffold848     Ok      -504.258
		-1      FunhetEGm018322t1       Scaffold848     Ok      -805.335	

		FunhetEGm018322t1.splog8b
		+1      FunhetEGm018322t1       Scaffold848     Ok      -504.258
		-1      FunhetEGm018322t1       Scaffold848     Ok      -805.335
		+2      FunhetEGm018322t1       Scaffold630     Ok      -1193.53 << 2nd part here



CASE4: split  FunhetEGm067171t5/s27split = Scaffold10067,Scaffold9902 
	$nbin/compart $cpopts  -qdb $qfile1 -sdb $genome  > $qfile1.compart5  

  newspan: Scaffold9902:670788-801381/lend , Scaffold10067:108914-209862/rend
     159 FunhetEGm067171t5.mrna.compart5  : none orig run
  
	$nbin/splign $splopt -comps $qfile1.compart5 -blastdb $genome -ldsdir ldsdir \
	  -log $qfile1.splog -aln $qfile1.aln > $qfile1.splign

  FunhetEGm067171t5.mrna.splog
  +1	FunhetEGm067171t5	Scaffold10067	Ok	-1797.34  << best, 63/93 stdsplice,fwd-fwd
  -1	FunhetEGm067171t5	Scaffold10067	Ok	-1892.56  < worse, 23/89 stdsplice,rev-rev
  +2	FunhetEGm067171t5	Scaffold9902	Ok	-1920.22  << best, 69/92 stdsplice,fwd-fwd
  -2	FunhetEGm067171t5	Scaffold9902	Ok	-2175.28  <? 34/87 stdsplice, rev-rev


CASE3: FunhetEGm073739t1/Scaffold10038
	-min_idty 0.45; $nbin/compart $cpopts  -qdb $qfile1 -sdb $genome  > $qfile1.compart5  
	newspan: Scaffold10038:1175876-1216671:+  # gm is partial align here. this span just part of gene.. gaps
	
	$nbin/splign $splopt -comps $qfile1.compart5 -blastdb $genome -ldsdir ldsdir \
	  -log $qfile1.splog -aln $qfile1.aln > $qfile1.splign
		^^ splign result similar to t2+ alt mapping, ok, gaps due to geno gaps.
		
FunhetEGm073739t1.mrna.compart5
lcl|FunhetEGm073739t1   lcl|Scaffold10038       100     17      0       0       8631    8647    1216655 1216671 0       34
lcl|FunhetEGm073739t1   lcl|Scaffold10038       87.0968 31      4       0       8594    8624    1216621 1216651 0       62
lcl|FunhetEGm073739t1   lcl|Scaffold10038       96.2025 79      3       0       7109    7187    1214040 1214118 0       158
..
lcl|FunhetEGm073739t1   lcl|Scaffold10038       100     87      0       0       766     852     1176638 1176724 0       174
lcl|FunhetEGm073739t1   lcl|Scaffold10038       99.3141 717     4       0       34      750     1175909 1176625 0       1434
lcl|FunhetEGm073739t1   lcl|Scaffold10038       93.3333 30      2       0       1       30      1175876 1175905 0       60


CASE2: FunhetEGm018012t1/Scaffold521
	-min_idty 0.45; $nbin/compart $cpopts  -qdb $qfile1 -sdb $genome  > $qfile1.compart5  
	newspan:  Scaffold521:144238-188797:+ << this is gmap span
	$nbin/splign $splopt -comps $qfile1.compart5 -blastdb $genome -ldsdir ldsdir \
	  -log $qfile1.splog -aln $qfile1.aln > $qfile1.splign
	
FunhetEGm018012t1.mrna.compart5
lcl|FunhetEGm018012t1	lcl|Scaffold521	99	100	1	0	7410	7509	188698	188797	0	200
lcl|FunhetEGm018012t1	lcl|Scaffold521	100	146	0	0	7264	7409	188279	188424	0	292
lcl|FunhetEGm018012t1	lcl|Scaffold521	100	132	0	0	7010	7141	178766	178897	0	264
lcl|FunhetEGm018012t1	lcl|Scaffold521	100	124	0	0	6885	7008	178507	178630	0	248
lcl|FunhetEGm018012t1	lcl|Scaffold521	98.8764	83	0	0	6797	6879	178313	178395	0	166
lcl|FunhetEGm018012t1	lcl|Scaffold521	100	128	0	0	6668	6795	178097	178224	0	256
..
lcl|FunhetEGm018012t1	lcl|Scaffold521	99.4898	194	0	0	627	820	149559	149752	0	388
lcl|FunhetEGm018012t1	lcl|Scaffold521	100	42	0	0	582	623	148245	148286	0	84
lcl|FunhetEGm018012t1	lcl|Scaffold521	98.6755	446	4	0	131	576	147791	148236	0	892
lcl|FunhetEGm018012t1	lcl|Scaffold521	93.0769	128	8	0	1	128	144238	144365	0	256


CASE1: FunhetEGm025643t1/Scaffold755
	newspan: Scaffold755:37619-137803:+
	
Lower opts:
$nbin/compart $cpopt -qdb FunhetEGm025643t1.mrna -sdb Scaffold755.gdna
	new compart opts gets Scaffold755:37619-77998,102301-137803 / tr:4511-8860,8863-15043
	^^ this is correct, introns cover this span, extends lower than gmap, catfish/zfish match full span.
	
Remap data created for Scaffold755.gdna; max offset = 146636
...
 Generating index volume: Scaffold755.gdna.424932660.p.v1 ... Ok
 Matching (strand = plus) ... 
lcl|FunhetEGm025643t1   lcl|Scaffold755 90      49      4       0       14995   15043   137755  137803  0       98
lcl|FunhetEGm025643t1   lcl|Scaffold755 82.2581 61      10      0       14933   14993   137688  137748  0       122
lcl|FunhetEGm025643t1   lcl|Scaffold755 94.1176 51      3       0       14870   14920   137624  137674  0       102
lcl|FunhetEGm025643t1   lcl|Scaffold755 95.8425 457     19      0       14398   14854   137151  137607  0       914
lcl|FunhetEGm025643t1   lcl|Scaffold755 99.6974 1319    3       0       13077   14395   135830  137148  0       2638
lcl|FunhetEGm025643t1   lcl|Scaffold755 99.8821 1648    0       0       11424   13071   134179  135826  0       3296
lcl|FunhetEGm025643t1   lcl|Scaffold755 100     155     0       0       11165   11319   133979  134133  0       310
..
lcl|FunhetEGm025643t1   lcl|Scaffold755 99.6109 253     0       0       8863    9115    102301  102553  0       506
lcl|FunhetEGm025643t1   lcl|Scaffold755 100     185     0       0       8676    8860    77814   77998   0       370
..
lcl|FunhetEGm025643t1   lcl|Scaffold755 100     128     0       0       5130    5257    51280   51407   0       256
lcl|FunhetEGm025643t1   lcl|Scaffold755 100     284     0       0       4840    5123    37619   37902   0       568
lcl|FunhetEGm025643t1   lcl|Scaffold755 100     158     0       0       4511    4668    3174    3331    0       316

Defaults:
$nbin/compart -qdb FunhetEGm025643t1.mrna -sdb Scaffold755.gdna
	new compart opts gets Scaffold755:37619-77998,102301-137803 / tr:4511-8860,8863-15043
 	.. 100% align for 37619-77998, 82%..100 for rest.
  gm: Scaffold755:62292-137840,cov=52%
 
$nbin/compart  -qdb $qfile1 -sdb $genome  | less
 Remap data created for Scaffold755.gdna; max offset = 146636
 Scanning 1 genomic sequences ... Ok
 Constructing FV ... Ok
 Remap data created for sequences; max offset = 15611
 Scanning sequences for N-mers and their positions.
 Generating index volume: qq.1505824399.p.v1 ... Ok
 Scanning Scaffold755.gdna for N-mers and their positions.
 Generating index volume: Scaffold755.gdna.1505824399.p.v1 ... Ok
 Matching (strand = plus) ... Ok
 Reading/transforming FV ... Ok
 Scanning sequences for N-mers and their positions.
 Generating index volume: qq.1505824399.m.v1 ... Ok
 Scanning Scaffold755.gdna for N-mers and their positions.
 Generating index volume: Scaffold755.gdna.1505824399.m.v1 ... Ok
 Matching (strand = minus) ... Ok

=cut

=item splign table

	spl25/kfish2asm.gdna-kfish2p67vs.mrna.split.3.splign
	
	+1	FunhetEGm046329t1	Scaffold0	0.994	344	1	344	2942702	2943045	  <exon>  	M89RMRM252   # msum= 252+1+89 = 342/344
	+2	FunhetEGm046569t1	Scaffold0	1	88	1	88	5831260	5831347	  <exon>GT	M88
	+2	FunhetEGm046569t1	Scaffold0	1	63	89	151	5831489	5831551	AG<exon>GT	M63
	+2	FunhetEGm046569t1	Scaffold0	1	212	152	363	5832998	5833209	AG<exon>  	M212
	+3	FunhetEGm047943t1	Scaffold0	-	5	1	5	-	-	<L-Gap>	-
	+3	FunhetEGm047943t1	Scaffold0	1	207	6	212	541156	541362	TT<exon>GT	M207
	+3	FunhetEGm047943t1	Scaffold0	1	425	213	637	544630	545054	AG<exon>GT	M425
	+3	FunhetEGm047943t1	Scaffold0	1	24	638	661	618363	618386	AG<exon>  	M24

=item splign logfile : useful?
	# spl25/kfish2asm.gdna-kfish2p67vs.mrna.split.3.splog 
	
	+1	FunhetEGm046329t1	Scaffold0	Ok	-8.248		: has only 1 exon partly aligned, thus low -score
	+2	FunhetEGm046569t1	Scaffold0	Ok	-17.252
	+3	FunhetEGm047943t1	Scaffold0	Ok	-28.458
	+4	FunhetEGm047950t1	Scaffold0	Ok	-28.65
	+5	FunhetEGm047996t1	Scaffold0	Ok	-16.632
	+6	FunhetEGm048121t1	Scaffold0	Ok	-18.912
	+7	FunhetEGm048618t1	Scaffold0	Ok	-151.9
	+8	FunhetEGm057611t1	Scaffold0	Ok	-333.571
	+9	FunhetEGm056896t1	Scaffold0	Ok	-44.051		: + = forward
	-9	FunhetEGm056896t1	Scaffold0	Ok	-44.051		: - = reverse

=item eg FunhetEGm056896t1
	# spl25/kfish2asm.gdna-kfish2p67vs.mrna.split.3.splog 
	-- several mappings, probably take highest -score, skip any 2ndary below <99% of top?
	+9	FunhetEGm056896t1	Scaffold0	Ok	-44.051
	-9	FunhetEGm056896t1	Scaffold0	Ok	-44.051
	+875	FunhetEGm056896t1	Scaffold348	Ok	-60.039    < best score, 
	-875	FunhetEGm056896t1	Scaffold348	Ok	-56.914
	+3190	FunhetEGm056896t1	Scaffold9988	Ok	-49.873
	-3190	FunhetEGm056896t1	Scaffold9988	Ok	-49.873
	+5626	FunhetEGm056896t1	Scaffold936	Ok	-47.645
	-5626	FunhetEGm056896t1	Scaffold936	Ok	-47.645
	+8236	FunhetEGm056896t1	Scaffold10172	Ok	-37.791
	-8236	FunhetEGm056896t1	Scaffold10172	Ok	-37.791
	
	# are these horrid long mismatch strings useful?
	#   M10RM9RM9RM16RMRM12RM23RM11RM6RM6RM51RM2I8MRM3IMRM2RMRMRMRM2D2M2RM33RM8RM3RM37RM62
	+875	FunhetEGm056896t1	Scaffold348	0.905	346	1	337	76738	77081	  <exon>GG	M10RM9RM9RM16RMRM12RM23RM11RM6RM6RM51RM2I8MRM3IMRM2RMRMRMRM2D2M2RM33RM8RM3RM37RM62
	+875	FunhetEGm056896t1	Scaffold348	-      36	338	373	-	-	<M-Gap>	-
	+875	FunhetEGm056896t1	Scaffold348	0.905	231	374	603	77082	77306	GA<exon>GT	M6RM15RM28RM24RM42DM41RM6IM2RMR2MD5M3R3M3R3M14RM23
	+875	FunhetEGm056896t1	Scaffold348	0.933	360	604	963	79954	80313	GG<exon>  	M16RM18RM28RM43RM16RM79RMRM12RM16RM2RM3R2MRM4RM2RM2RM20RM8RM5RM2RM5RM2RM12RM22RM17
	
	# .. rev map is slightly different: 77082 .. 77306,77313 
	-875	FunhetEGm056896t1	Scaffold348	0.932	353	963	611	80313	79961	  <exon>GT	M17RM22RM12RM2RM5RM2RM5RM8RM20RM2RM2RM4RMR2M3RM2RM16RM12RMRM79RM16RM43RM28RM18RM9
	-875	FunhetEGm056896t1	Scaffold348	0.908	238	610	374	77313	77082	TG<exon>TC	M30RM14R3M3R3M3R2M2RD5M3IM5RM41DM42RM24RM28RM15RM6
	-875	FunhetEGm056896t1	Scaffold348	-	     36	373	338	-	-	<M-Gap>	-
	-875	FunhetEGm056896t1	Scaffold348	0.905	346	337	1	77081	76738	CC<exon>  	M62RM37RM3RM8RM33RM2D2M2RMRMRMRM2RMIM3RM2I8MRM51RM6RM6RM11RM23RM12RMRM16RM9RM9RM10

=item eg. +8	FunhetEGm057611t1	Scaffold0	Ok	-333.571

	+8	FunhetEGm057611t1	Scaffold0	-	198	1	198	-	-	<L-Gap>	-
	+8	FunhetEGm057611t1	Scaffold0	1	176	199	374	2463919	2464094	AA<exon>GT	M176
	+8	FunhetEGm057611t1	Scaffold0	1	131	375	505	2465439	2465569	AG<exon>GT	M131
	+8	FunhetEGm057611t1	Scaffold0	1	505	506	1010	2468703	2469207	AG<exon>GT	M505
	+8	FunhetEGm057611t1	Scaffold0	1	771	1011	1781	2472316	2473086	AG<exon>GT	M771
	+8	FunhetEGm057611t1	Scaffold0	1	195	1782	1976	2474583	2474777	AG<exon>GT	M195
	+8	FunhetEGm057611t1	Scaffold0	1	164	1977	2140	2476325	2476488	AG<exon>GT	M164
	+8	FunhetEGm057611t1	Scaffold0	1	162	2141	2302	2476610	2476771	AG<exon>GT	M162
	+8	FunhetEGm057611t1	Scaffold0	1	77	2303	2379	2476868	2476944	AG<exon>GT	M77
	+8	FunhetEGm057611t1	Scaffold0	1	223	2380	2602	2478033	2478255	AG<exon>GT	M223
	+8	FunhetEGm057611t1	Scaffold0	1	234	2603	2836	2478338	2478571	AG<exon>GT	M234
	+8	FunhetEGm057611t1	Scaffold0	1	92	2837	2928	2478671	2478762	AG<exon>GT	M92
	+8	FunhetEGm057611t1	Scaffold0	1	128	2929	3056	2478935	2479062	AG<exon>GT	M128
	+8	FunhetEGm057611t1	Scaffold0	1	155	3057	3211	2479950	2480104	AG<exon>GT	M155
	+8	FunhetEGm057611t1	Scaffold0	1	135	3212	3346	2480194	2480328	AG<exon>GT	M135
	+8	FunhetEGm057611t1	Scaffold0	1	195	3347	3541	2480501	2480695	AG<exon>GT	M195
	+8	FunhetEGm057611t1	Scaffold0	1	183	3542	3724	2480929	2481111	AG<exon>GT	M183
	+8	FunhetEGm057611t1	Scaffold0	1	205	3725	3929	2481201	2481405	AG<exon>GT	M205
	+8	FunhetEGm057611t1	Scaffold0	1	147	3930	4076	2481498	2481644	AG<exon>GT	M147
	+8	FunhetEGm057611t1	Scaffold0	1	143	4077	4219	2481731	2481873	AG<exon>GT	M143
	+8	FunhetEGm057611t1	Scaffold0	1	204	4220	4423	2482561	2482764	AG<exon>GT	M204
	+8	FunhetEGm057611t1	Scaffold0	1	118	4424	4541	2482958	2483075	AG<exon>GT	M118
	+8	FunhetEGm057611t1	Scaffold0	1	103	4542	4644	2483719	2483821	AG<exon>GT	M103
	+8	FunhetEGm057611t1	Scaffold0	0.994	464	4645	5108	2483953	2484415	AG<exon>AA	M425RM7RM4DM25

=cut

=item problem case, low qual part-align of 3rd part bumped out hi qual align of 2nd middle

grep  Funhe2EKm005482t1 spldt1p2/*splign | sort -k6,6n -k4,4nr | egrep -v '[LR]-Gap'
  .. part 1 ..
+17786	Funhe2EKm005482t1	KN811587.1	1	172	4	175	855049	855220	TC<exon>GT	M172
+17786	Funhe2EKm005482t1	KN811587.1	1	125	176	300	883931	884055	AG<exon>GT	M125
+17786	Funhe2EKm005482t1	KN811587.1	1	102	301	402	885810	885911	AG<exon>GT	M102
+17786	Funhe2EKm005482t1	KN811587.1	1	101	403	503	890511	890611	AG<exon>GT	M101
+17786	Funhe2EKm005482t1	KN811587.1	0.901	121	504	624	901845	901962	AG<exon>AA	M19RM74D3M4R2M2RM2R3M2RMRM5
  .. part 2 ..
+17207	Funhe2EKm005482t1	KN811549.1	0.996	248	596	843	284223	283976	CA<exon>GT	M168RM79
+17207	Funhe2EKm005482t1	KN811549.1	1	121	844	964	283745	283625	AG<exon>GT	M121
+17207	Funhe2EKm005482t1	KN811549.1	0.993	137	965	1101	282493	282357	AG<exon>GT	M32RM104
+17207	Funhe2EKm005482t1	KN811549.1	1	213	1102	1314	254498	254286	AG<exon>GT	M213
+17207	Funhe2EKm005482t1	KN811549.1	1	175	1315	1489	254063	253889	AG<exon>GT	M175
+17207	Funhe2EKm005482t1	KN811549.1	1	121	1490	1610	220301	220181	AG<exon>GT	M121
+17207	Funhe2EKm005482t1	KN811549.1	1	164	1611	1774	218788	218625	AG<exon>GT	M164
+17207	Funhe2EKm005482t1	KN811549.1	0.99	103	1775	1877	217196	217094	AG<exon>GT	M50RM52
+17207	Funhe2EKm005482t1	KN811549.1	1	147	1878	2024	211408	211262	AG<exon>GT	M147
+17207	Funhe2EKm005482t1	KN811549.1	1	135	2025	2159	209759	209625	AG<exon>GT	M135
  .. part 3/KN805720 *bad* starts here, b=2026, ident=.70-.80 but part 2/KN811549 has ident=1.00 here
+2903	Funhe2EKm005482t1	KN805720.1	0.806	134	2026	2159	221716	221849	GT<exon>GC	M7RMRM2RM8RMR2M12RM7RM5RM2RM5RMR2M2R2M7RM8RM2RMRM3RM2RM8RM6R3M8RM10
+17207	Funhe2EKm005482t1	KN811549.1	0.992	129	2160	2288	200528	200400	AG<exon>GT	M100RM28
+2903	Funhe2EKm005482t1	KN805720.1	0.791	129	2160	2288	228413	228541	AG<exon>GT	MRM5RM2RM2RMRM7R2M3RM3RMRM4RM9RM9RM14RMRM8RM2RM5RMRM3RM7RM4RM4RM5R2MR2
+17207	Funhe2EKm005482t1	KN811549.1	1	129	2289	2417	200111	199983	AG<exon>GT	M129
+2903	Funhe2EKm005482t1	KN805720.1	0.76	129	2289	2417	229846	229974	AG<exon>GT	MRM5R2MRM2R4M8RM2RM2R2M2RM2RM23RM7RM3R3M5RM5R4M8RM2RM5RM2RM2RMRM6RM4
+17207	Funhe2EKm005482t1	KN811549.1	1	129	2418	2546	199009	198881	AG<exon>GT	M129
+2903	Funhe2EKm005482t1	KN805720.1	0.76	129	2418	2546	230616	230744	AG<exon>GT	MRM5RM2RM2RM8R5MRM8RM2RM5RM11RM2RM12RMRM5RM5R4M3RMRM11RM3RMR2M4RM4R2M
  .. part 3 *should* start here, b=2547, ident=.96 ok  ..
+2903	Funhe2EKm005482t1	KN805720.1	0.961	129	2547	2675	230863	230991	AG<exon>GT	M7R2MRM2RM2RM112
+17207	Funhe2EKm005482t1	KN811549.1	0.86	129	2547	2675	198800	198672	AG<exon>GT	M26RM13RM9RM7RM2RM2RM8RM5RM5RM8RM3R2M3RMR2M5RM11RM2RM
+2903	Funhe2EKm005482t1	KN805720.1	1	237	2676	2912	231150	231386	AG<exon>GT	M237
+2903	Funhe2EKm005482t1	KN805720.1	1	170	2913	3082	232804	232973	AG<exon>GT	M170
+2903	Funhe2EKm005482t1	KN805720.1	-	62	3083	3144	-	-	<M-Gap>	-
+2903	Funhe2EKm005482t1	KN805720.1	1	47	3145	3191	238014	238060	TC<exon>GT	M47
+2903	Funhe2EKm005482t1	KN805720.1	1	207	3192	3398	238700	238906	AG<exon>GT	M207
+2903	Funhe2EKm005482t1	KN805720.1	0.989	1628	3399	5026	240036	241646	AG<exon>  	M309RM504D16M202DM595

>> got this gff, missing middle found above
grep  Funhe2EKm005482t1 spldt1p2/*.gff | grep -v CDS
KN811587.1	splkf3n15	mRNA	855049	901962	12	+	.	ID=Funhe2EKm005482t1_C1;Split=1/2;cov=12%,609/5026;pid=98.2;nexon=5;splice=16;Target=Funhe2EKm005482t1 4 624;gaps=4405,LGap:1-3,RGap:625-5026,;gescore=12;aalen=1126,67%,complete;clen=5026;offs=194-3574;oid=Funhe2Exx11m002408t2,Funhe2E6bm002431t6;Name=Multiple epidermal growth factor domains protein 11 (97%T)
KN811587.1	splkf3n15	exon	855049	855220	1	+	.	Parent=Funhe2EKm005482t1_C1;Target=Funhe2EKm005482t1 4 175;splice=TCGT+
KN811587.1	splkf3n15	exon	883931	884055	1	+	.	Parent=Funhe2EKm005482t1_C1;Target=Funhe2EKm005482t1 176 300;splice=AGGT+
KN811587.1	splkf3n15	exon	885810	885911	1	+	.	Parent=Funhe2EKm005482t1_C1;Target=Funhe2EKm005482t1 301 402;splice=AGGT+
KN811587.1	splkf3n15	exon	890511	890611	1	+	.	Parent=Funhe2EKm005482t1_C1;Target=Funhe2EKm005482t1 403 503;splice=AGGT+
KN811587.1	splkf3n15	exon	901845	901962	0.901	+	.	Parent=Funhe2EKm005482t1_C1;Target=Funhe2EKm005482t1 504 624;splice=AGAA+
  .. missing  KN811549.1 part ...
KN805720.1	splkf3n15	mRNA	221716	241646	56	+	.	ID=Funhe2EKm005482t1_C2;Split=2/2;cov=56%,2801/5026;pid=96;nexon=10;splice=32;Target=Funhe2EKm005482t1 2026 5026;gaps=2087,LGap:1-2025,MGap:3083-3144,;gescore=52;aalen=1126,67%,complete;clen=5026;offs=194-3574;oid=Funhe2Exx11m002408t2,Funhe2E6bm002431t6;Name=Multiple epidermal growth factor domains protein 11 (97%T)
KN805720.1	splkf3n15	exon	221716	221849	0.806	+	.	Parent=Funhe2EKm005482t1_C2;Target=Funhe2EKm005482t1 2026 2159;splice=GTGC+
KN805720.1	splkf3n15	exon	228413	228541	0.791	+	.	Parent=Funhe2EKm005482t1_C2;Target=Funhe2EKm005482t1 2160 2288;splice=AGGT+
KN805720.1	splkf3n15	exon	229846	229974	0.76	+	.	Parent=Funhe2EKm005482t1_C2;Target=Funhe2EKm005482t1 2289 2417;splice=AGGT+
KN805720.1	splkf3n15	exon	230616	230744	0.76	+	.	Parent=Funhe2EKm005482t1_C2;Target=Funhe2EKm005482t1 2418 2546;splice=AGGT+
KN805720.1	splkf3n15	exon	230863	230991	0.961	+	.	Parent=Funhe2EKm005482t1_C2;Target=Funhe2EKm005482t1 2547 2675;splice=AGGT+
KN805720.1	splkf3n15	exon	231150	231386	1	+	.	Parent=Funhe2EKm005482t1_C2;Target=Funhe2EKm005482t1 2676 2912;splice=AGGT+
KN805720.1	splkf3n15	exon	232804	232973	1	+	.	Parent=Funhe2EKm005482t1_C2;Target=Funhe2EKm005482t1 2913 3082;splice=AGGT+
KN805720.1	splkf3n15	exon	238014	238060	1	+	.	Parent=Funhe2EKm005482t1_C2;Target=Funhe2EKm005482t1 3145 3191;splice=TCGT+
KN805720.1	splkf3n15	exon	238700	238906	1	+	.	Parent=Funhe2EKm005482t1_C2;Target=Funhe2EKm005482t1 3192 3398;splice=AGGT+
KN805720.1	splkf3n15	exon	240036	241646	0.989	+	.	Parent=Funhe2EKm005482t1_C2;Target=Funhe2EKm005482t1 3399 5026;splice=AGnn+

=cut

=item genosplign.sh cluster script

	see evigene/scripts/rnaseq/genosplign.sh
	-------------
	#! /bin/bash
	### env genome=xxx mrna=mrna.fa datad=`pwd` qsub -q normal genosplign.sh
	# NCBI splign align mRNA to genome
	splopt="-type mrna -min_compartment_idty 0.5"
	
	# split mrna.fasta to ncpu parts, run ncpu cluster jobs, output to spl$i part folders
	qset=`/bin/ls $mrna.split.*.fa`
	i=0; for qfile in $qset; do {
		# setup inputs ..
	  # this is fork set:
  	echo "# splign $qfile1 x $genome in $ispldir";
  	( $nbin/splign -mklds ./; \
  	  ln -s ../$genome.* ./ ; \
      $nbin/makeblastdb -parse_seqids -dbtype nucl -in $qfile1; \  
      $nbin/compart  -qdb $qfile1 -sdb $genome  > $onam.cpart; \  
      $nbin/splign $splopt -comps $onam.cpart -blastdb $genome -ldsdir ./ -log $onam.splog > $onam.splign; ) &
		i=$(( $i + 1 ))
	} done
	wait
	#.................. 

=cut

=item compart options (blastn variant)


bin/compart -help
USAGE
  /bio/bio-grid/mb/ncbix/bin/compart [-h] [-help] [-xmlhelp]
    [-qdb qdb] [-sdb sdb] [-ho] [-penalty penalty] [-min_idty min_idty]
    [-min_singleton_idty min_singleton_idty]
    [-min_singleton_idty_bps min_singleton_idty_bps] [-max_intron max_intron]
    [-dropoff dropoff] [-min_query_len min_query_len]
    [-min_hit_len min_hit_len] [-maxvol maxvol] [-noxf] [-seqlens seqlens]
    [-N N] [-version-full] [-dryrun]

DESCRIPTION
   Compart v.1.35. Unless -qdb and -sdb are specified, the tool expects
   tabular blast hits at stdin collated by query and subject, e.g. with 'sort
   -k 1,1 -k 2,2'

OPTIONAL ARGUMENTS
 -h	 Print USAGE and DESCRIPTION;  ignore all other parameters
 -help	 Print USAGE, DESCRIPTION and ARGUMENTS; ignore all other parameters
 -xmlhelp	Print USAGE, DESCRIPTION and ARGUMENTS in XML format; ignore all other parameters
 -qdb <String>	 cDNA BLAST database
 -sdb <String>	 Genomic BLAST database
 -ho		Print raw hits only - no compartments
 -penalty <Real, 0..1>	Per-compartment penalty   Default = `0.55'
 -min_idty <Real, 0..1>
   Minimal overall identity. Note: in current implementation  there is no
   sense to set different 'min_idty' and 'min_singleton_idty' (minimum is used
   anyway).
   Default = `0.70'
 -min_singleton_idty <Real, 0..1>
   Minimal identity for singleton compartments. The actual parameter passed to
   the compartmentization procedure is least of this parameter multipled by
   the seq length, and min_singleton_idty_bps. Note: in current implementation
    there is no sense to set different 'min_idty' and 'min_singleton_idty'
   (minimum is used anyway).
   Default = `0.70'
 -min_singleton_idty_bps <Integer>
   Minimal identity for singleton compartments in base pairs. Default =
   parameter disabled.
   Default = `9999999'
 -max_intron <Integer>
   Maximum intron length (in base pairs)
   Default = `1200000'
 -dropoff <Integer>
   Max score drop-off during hit extension.
   Default = `5'
 -min_query_len <Integer, 21..99999>
   Minimum length for individual cDNA sequences.
   Default = `50'
 -min_hit_len <Integer, 1..99999>
   Minimum length for reported hits in hits-only mode. No effect in
   compartments mode.
   Default = `16'
 -maxvol <Integer, 128..1024>
   Maximum index volume size in MB (approximate)
   Default = `512'
 -noxf
   [With external hits] Suppress overlap x-filtering: print all compartment
   hits intact.
 -seqlens <File_In>
   [With external hits] Two-column file with sequence IDs and their lengths.
   If none supplied, the program will attempt fetching the lengths from
   GenBank. Cannot be used with -qdb.
 -N <Integer>
   [With external hits] Max number of compartments per query (0 = All).
   Default = `0'
 -version-full
   Print extended version data;  ignore other arguments
 -dryrun
   Dry run the application: do nothing, only test all preconditions


=item splign options

bin/splign -help
USAGE
  /bio/bio-grid/mb/ncbix/bin/splign [-hits hits]
    [-comps comps] [-mklds mklds] [-blastdb blastdb] [-ldsdir ldsdir]
    [-query query] [-subj subj] [-disc] [-W mbwordsize] [-type type]
    [-compartment_penalty compartment_penalty]
    [-min_compartment_idty min_compartment_identity]
    [-min_singleton_idty min_singleton_identity]
    [-min_singleton_idty_bps min_singleton_identity_bps]
    [-min_exon_idty identity] [-min_polya_ext_idty identity]
    [-min_polya_len min_polya_len] [-max_intron max_intron]
    [-max_space max_space] [-direction direction] [-log log] [-asn asn]
    [-aln aln]

DESCRIPTION
   Splign: 1.39.8
   Package: public 12.0.0, build Jul 18 2013 16:30:51
    LLVMGCC_421-Debug64--x86_64-apple-darwin12.4.0-c_67_173_144_85
    
OPTIONAL ARGUMENTS
 -hits <File_In>
   [Batch mode] Externally computed local alignments (such as blast hits), in
   blast tabular format. The file must be collated by subject and query (e.g.
   sort -k 2,2 -k 1,1).
 -comps <File_In>
   [Batch mode] Compartments computed with Compart utility.
 -mklds <String>
   [Batch mode] Make LDS DB under the specified directory with cDNA and
   genomic FASTA files or symlinks.
 -blastdb <String>
   [Batch mode] Blast DB.
 -ldsdir <String>
   [Batch mode] Directory holding LDS subdirectory.
 -query <File_In>
   [Pairwise mode] FASTA file with the spliced sequence.
 -subj <File_In>
   [Pairwise mode] FASTA file with the genomic sequence.
 -disc
   [Pairwise mode] Use discontiguous megablast to facilitate alignment of more
   divergent sequences such as those from different organisms (cross-species
   alignment).
 -W <Integer>
   [Pairwise mode] Megablast word size
   Default = `28'
 -type <String, `est', `mrna'>
   Query cDNA type: 'mrna' or 'est'
   Default = `mrna'
 -compartment_penalty <Real, 0..1>
   Penalty to open a new compartment (compartment identification parameter).
   Multiple compartments will only be identified if they have at least this
   level of coverage.
   Default = `0.55'
 -min_compartment_idty <Real, 0..1>
   Minimal compartment identity to align.
   Default = `0.7'
 -min_singleton_idty <Real>
   Minimal singleton compartment identity to use per subject and strand,
   expressed as a fraction of the query's length.
 -min_singleton_idty_bps <Integer>
   Minimal singleton compartment identity to use per subject and strand, in
   base pairs. The actual value passed to the compartmentization procedure is
   the least of (min_singleton_idty * query_length) and
   min_singleton_identity_bps.
   Default = `9999999'
 -min_exon_idty <Real, 0..1>
   Minimal exon identity. Segments with lower identity will be marked as gaps.
   Default = `0.75'
 -min_polya_ext_idty <Real, 0..1>
   Minimal identity to extend alignment into polya. Polya candidate region on
   mRNA is detected first. Alignment is produced without the polya candidate
   region After that alignment will be extended into the polya candidate
   region to deal with case when initial polya detection was wrong
   Default = `1'
 -min_polya_len <Integer, 1..1000000>
   Minimal length of polya.
   Default = `1'
 -max_intron <Integer, 7..2000000>
   The upper bound on intron length, in base pairs.
   Default = `1200000'
 -max_space <Real, 500..4096>
   The max space to allocate for a splice, in MB. Specify lower values to
   spend less time stitching over large genomic intervals.
   Default = `4096'
 -direction <String, `antisense', `auto', `both', `default', `sense'>
   Query sequence orientation. Auto orientation begins with the longest ORF
   direction (d1) and proceeds with the opposite direction (d2) if found a
   non-consensus splice in d1 or poly-a tail in d2. Default translates to
   'auto' in mRNA and 'both' in EST mode
   Default = `default'
 -log <File_Out>
   Splign log file
   Default = `splign.log'
 -asn <File_Out>
   ASN.1 output file name
 -aln <File_Out>
   Pairwise alignment output file name
   
=cut
