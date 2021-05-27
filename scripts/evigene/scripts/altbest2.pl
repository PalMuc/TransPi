#!/usr/bin/env perl
# altbest.pl : pull valid alt-transcript assemblies

=item about

 altbest.pl : pull valid alt-transcript assemblies
    input: bestgenes.gff + all_altasm.gff
    input exons must be marked with intr=..N1,N2, intron splice ids
    input gff must be gene-record ordered (mRNA,exons)..
    
 intron chain annotation:
   $workd/scripts/overlapfilter -intron2splice=error -pass 'exon,intron' -act markid -midtype scoresum \
   -mark intr -over $workd/intron/intron_good.gff.gz -in genes.gff \

 steps: 
   1. equivalence allalt.gff to bestgenes (CDS equiv >= 66 %?)
   2. pull intron chains of equivalent alts, filter only uniq, drop short chains

=item major update 1607

	-- various mods for splitgenes (ID suffix test)
	-- keepids=list, other corrections, 
  -- dropnomain = retain alts.gff not in main-alt equal table
  	>> this results in many "new" loci, begs new method to id new loci + alts from
  		shared intron chains
  		
=item fixmees

	FIXME3: Nopath (unmapped) alts are kicked out, should not do that .. but no introns for such
			.. however, qlen=bogus from gmap-nopath, needed below.
			.. also may want to keep alts w/ mapto -introns but not annotated valid intr=
			
	** FIXME2: 
	a. done? revise to handle valid alternate exon ends (e.g. 1 alt continues w/ intron, 2nd ends cds)
	
  b. option to use genes.gff introns for valid intron chains, instead/with evd-intron annots.
     ^^b. do here?, or let caller make intron.tab and annotate genes, per above, would be simpler.
          i.e. gmap.out, gmap.gff output has valid intron splices for trgenes, use that?
  
  b2. input introns.gff/table instead of requiring genes.gff exon intr= annotation
      .. need parts from overlapfilter -intron2splice=error
          
 ** FIXME altbest.pl : reduce for partial prot vs shorter full prot of ~same rna models
 .. problem from pasa_asm > longest-orf-partials not as good as full prot alt models.
   
=cut

use FindBin;
use lib ("$FindBin::Bin"); # assume evigene/scripts/ layout: main, prot/, rnaseq/, ...

use strict;
use warnings;
use Getopt::Long;

use constant VERSION => '2016.07.20'; # altbest v2 ?

our $EVIGENES="$FindBin::Bin";  # $evigene="/bio/bio-grid/mb/evigene/scripts";
my $MINID_CDS=33; #was 66;
my $MINID_UTR=33; #was 40; # this is really full transcript span percent, utr+cds
my $MIN_PCDS=40;  # 0=all, assume mRNA cxlen=cds,exon or calculate? use to filter falseutr models
my $MIN_CDS2MAIN= 0.20; #was 0.55; # relative to main pcds, ie drop false long utr aberrations
my $MIN_ALTCHAINSIZE = 0.30; # OPTION: min proportion of introns in short alt, of longest
my $TAG_ALT='t'; # append/replace to prime id
my $CHANGEALTID=0; # 1802: not default now, opt -changealtid needed 
my $IDPREFIX= $ENV{idpre}||"newlocA";
my $ADDMAINID=1; #def:on
my ($maingenes, $altgenes, $eqtab, $outbase, $debug, $dropnomain, $skipnewloci, $keepnointron, $keepids, $dropids)= (0) x 9;
$dropnomain = 1; # 1803upd default; -nodropnomain to turn off
# maybe: $keepnointron = 1; # 1803upd default

my $optok= GetOptions(
  "altgenes=s", \$altgenes, 
  "maingenes=s", \$maingenes,  
  "eqtab=s", \$eqtab, 
  "outbase=s", \$outbase, 
  "mincds=i", \$MINID_CDS,  
  "minexon|minutr=i", \$MINID_UTR,  
  "minintronchainsize=s", \$MIN_ALTCHAINSIZE,  
  "mincoding=i", \$MIN_PCDS,  
  "minpcds2main=s", \$MIN_CDS2MAIN,  
  "keepnointron!", \$keepnointron, # is this flag or value option? to say keep NOPATH/others w/o intr= annots
  "idprefix=s", \$IDPREFIX,  #1607, newloc only?
  "tagalt=s", \$TAG_ALT, 
  "keepids=s", \$keepids,  #1607
  "dropids=s", \$dropids,  #1703, drop unless in keep set?
  "dropnomain!", \$dropnomain,  # 1607 = MAKE DEFAULT, alts w/o main eqgene overlap, output.altbest.gff , as newloc?
  "skipnewloci!", \$skipnewloci,  # debug newlocustables first
  "CHANGEALTID!", \$CHANGEALTID,  # noCHANGEALTID turns off, default on?
  "ADDMAINID!", \$ADDMAINID,  # upd1802 : always add to alt mrna annot?
  "debug!", \$debug, 
  ## separate debug opt to write invalidchain to ATAB ?
  );

die "usage: altbest.pl -main bestgenes.gff -alt allgenes.gff
   options: -mincds=$MINID_CDS -minexon=$MINID_UTR -mincoding=$MIN_PCDS -nochangealtid
Note: input gff must have exon intr=..,N1,N2,N3 intron chain annots.
" unless($optok and (-f $maingenes and -f $altgenes));

$MIN_CDS2MAIN=$MIN_CDS2MAIN/100 if($MIN_CDS2MAIN>1);
$MIN_ALTCHAINSIZE=$MIN_ALTCHAINSIZE/100 if($MIN_ALTCHAINSIZE>=1);

# change to -type similarCDS / utr and filter input by min cds/utr
# my $eqcmd="$EVIGENES/overgenedup.pl -exon=CDS -type CDS -mincds=$MINID_CDS -slopexon=8 -act markid -mark overg "
# ." -in $maingenes -over $altgenes";

my $eqcmd="$EVIGENES/overgenedup.pl -type similarCDS -mincds=$MINID_CDS  -minutr=$MINID_UTR -act markid -mark overg "
         ." -in $maingenes -over $altgenes";

# output files...
unless($eqtab) {
	$eqtab=$maingenes; $eqtab=~s/.gz//; $eqtab =~ s/.gff//; $eqtab .= ".equalalt.tab";
}
my $basetab = $outbase||$eqtab; $basetab =~ s/\.equalalt//; $basetab =~ s/\.eqgene//; $basetab=~s/\.tab//; 
my $altab 	= $basetab.".altbest.tab";
my $lnewtab	= $basetab.".newloci.tab";
my $valtgff = $altab; $valtgff =~ s/.tab//; $valtgff .=".gff";
my $lnewgff = $lnewtab; $lnewgff =~ s/.tab//; $lnewgff .=".gff";

my (%gmain,%galt,%ichain,%allalt2main,%allaltseen,%mainpcds,%gcds,%exonends,%alteqval);
my (%gcdslen,%keepids); #201607 split-gene fix?
my ($nerrmrna,$ninmrna,$nxin,$nxintr)=(0) x 10;

# for matchIntronChains()/putAltTables()
my %valclass=( 
	 1 => "okay", 2 => "okaltend", 3 => "oknewlocus", 4 => "okaykeep",
	-1 => "dropdup", -2 => "dropshort", -3 => "dropsubchain", 0 => "drop" );
	
my (%altnum, %didid, %galtok, %galtsame, %gmaindup, %gmaindupid, 
	  %alt2gid, %altbadchain, %nclass, $ntabmain, $ntabalt);

# ** FIXME: in busy regions get some bizarre assemblies that cover several genes,
# .. with long false UTRs; filter somewhere here, using multi-maingene overlaps?
# -- measure cdslen, exonlen and drop out false utr where cds/exon < 0.3,0.4 ?

warn "#options: -mincds=$MINID_CDS -minexon=$MINID_UTR -minintronchain=$MIN_ALTCHAINSIZE"
	." -mincoding=$MIN_PCDS -codingOfMain=$MIN_CDS2MAIN -tag=$TAG_ALT -keepids=$keepids\n" if $debug;

sub MAINstub {}

if($keepids) {
	open(IN,$keepids) or die "#ERR: cant read keepids=$keepids\n";
	while(<IN>) { my($id)=split; $keepids{$id}=1 if($id=~/^\w/); } close(IN);
}

my %dropids=(); if($dropids) {
	open(IN,$dropids) or die "#ERR: cant read dropids=$dropids\n";
	while(<IN>) { my($id)=split; $dropids{$id}=-1 if($id=~/^\w/ and not $keepids{$id}); } close(IN); # %keepids or %dropids ?
}

makeEqualgeneTab($eqtab,$eqcmd) unless( -f $eqtab);
readEqualgeneTab($eqtab);

collectIntronChains($maingenes, $altgenes);

matchIntronChains(); # and putAltTables()

# FIXME: new table of main-fragments, from eqgene, ismain that are subset overlaps to altgenes
# mainFragOfOtherGene(); # not right yet..

my($newalt,$newloci)= newlociFromIchains(); # TEST, table only now

putAltGFF($altgenes, $valtgff, $newalt, $newloci, $lnewgff);


#------- subs ------

sub splitID {
	my($gid,$attr)=@_; $attr||="";
	return($gid,$gid,0) unless($gid=~m/_C/ or $attr=~/Split=/); # most cases
 	my $gsplit=0; 
 	my $gidc= $gid; # gidc lacks split-tag, gid has it now
	if($gidc=~s/_C(\d+)$//) { $gsplit=$1; } # gmap2gff split-gene flag
	elsif($attr=~/;Split=(\w+)/) { $gsplit=$1; $gsplit=~s/^C//; $gid .="_C$gsplit"; } # make uniq id for splits ??
	return($gid,$gidc,$gsplit);
}

sub makeEqualgeneTab {
	#my($eqtab,$eqcmd)= @_;
	# make equiv table unless( -f $eqtab)  
	warn "# equalalt $eqtab = $eqcmd\n" if $debug;
	if(-s $eqtab) { warn "# equalalt reusing $eqtab\n"; return 0; }
	open(EOUT,">$eqtab") or die "writing $eqtab";
	open(ECMD,"$eqcmd |") or die "error $eqcmd";
	while(<ECMD>) {
		if(/\tmRNA/) { 
		my($id)=m/ID=([^;\s]+)/; 
		my($gid,$gidc,$gsplit)= splitID($id,$_); #* SPLIT id tags needed here also
		my($ov)= (m/overg=([^;\s]+)/) ? $1 : ""; next unless($ov);
		my $og= (m/oid=([^;\s]+)/) ? $1:"noid"; 
		my($r,$b,$e)=(split"\t")[0,3,4]; $r=~s/^Scaffold[_]?(\S+)/${1}sc/i; 
		print EOUT join("\t",$gid,$og,$ov,"$r:$b-$e"),"\n" if($ov);
		}
	} 
	close(ECMD);  close(EOUT);
}


sub readEqualgeneTab {
	# my($eqtab)= @_;
	# read equiv table
	
	# our (%alteqval); # global? FIXME: alts map onto 2+ gids ; last not best **
	sub gidbest{ my($aid,$gid)=@_;  # our (%alteqval);
		my @g= keys %{$alteqval{$aid}}; 
		if(@g>1) { my($gb)= sort { $alteqval{$aid}{$b} <=> $alteqval{$aid}{$a} or $a cmp $b} @g; return $gb; } 
		return $gid;
	}
	
	open(EOUT,$eqtab) or die "open $eqtab";
	while(<EOUT>) {
		my($gid,$altids)=(split"\t")[0,2];
		my($gidXXX,$gidc,$gsplit)=splitID($gid,"");
		
		# this filter is problem now, vs "newlocalt": $MINID_CDS, $MINID_UTR
		#v2: parse  ID/[CI]99.88,
		##upd17: record all CDS over:  0 < $c <$MINID also .. allaltseen
		
		my @altids= grep /\w/, map{ 
			my($c,$x)=m,/[CI]?(\d+)[.]?(\d*),; s,/.*$,,;
			$c||=0; $x=$c unless($x);
			$alteqval{$_}{$gid}="$c.$x";
			$_ if($c>=$MINID_CDS and $x>=$MINID_UTR); 
			} split",",$altids;
		
		#? should this track/log altids not retained in galt ? YES
		map{ s,/.*$,,; $allaltseen{$_}= $gid; } split",",$altids;
		 
		if(@altids) {
			# SPLIT ids: %galt, %allalt2main now have ID_C split tags
			$galt{$gid}= \@altids; # ok? eqtab gid col is uniq
			#o map{$allalt2main{$_}= $gid; } @altids; # FIXME: alts map onto 2+ gids ; last not best **
			map{ my($g)= gidbest($_,$gid); $allalt2main{$_}= $g; } @altids;  
			if($gidc ne $gid) { # what?
				if($galt{$gidc}) { } # merge ??
				else { $galt{$gidc}= \@altids; }
			}
		}
	}
}


sub collectIntronChains {
	#my($maingenes, $altgenes)= @_;
	# collect exon intron chains, only for equiv genes
	foreach my $gf ($maingenes, $altgenes) {
		my $ismain= ($gf eq $maingenes)?1:0;
		my $gop = ($gf =~ /.gz/) ? "gunzip -c $gf" : "cat $gf"; 
		# drop exon filter, test mRNA attribs.
		
		my ($gid,$gsplit,$keep,%gx,$iexon,$nexon); 
		open(G,"$gop |") or die "bad $gop";
		while(<G>){
			next unless(/^\w/);
			if(/\tmRNA/) {
				if($keep and $gid) {  ## FIXME3 here? for NOPATH, no intron alts; other mapped+nointron alts?
					my @xin=();
					if($gx{$gid}) { @xin= map{ $gx{$gid}{$_} } sort{ $a<=>$b } keys %{$gx{$gid}}; }
					# if($keepnointron and not @xin) { } # fake it? but how so we get some useful alts, but not lots o crap      
					$ichain{$gid}= join "; ",@xin if @xin; 
				}
					
				($gid)= m/ID=([^;\s]+)/; $ninmrna++;
				$keep=1; %gx= ();   #NO cant reset here: %gcds=();
				$iexon=0; ($nexon)= m/nexon=(\d+)/; $nexon||=1;
				
				# * need to handle split genes, _C[12] or ;Split=, but how?
				my $gidc= $gid; 
				($gid,$gidc,$gsplit)= splitID($gid,$_); #* SPLIT id tags needed here also
				# $gsplit=0;
				# if($gidc=~s/_C(\d+)$//) { $gsplit=$1; } # gmap2gff split-gene flag
				# elsif(/;Split=(\w+)/) { $gsplit=$1; $gsplit=~s/^C//; $gid .="_C$gsplit"; } # make uniq id for splits ??
				
				## FIXME: mRNA fields
				## gmap: ID=Funhe2Exy3m110884t1;trg=Funhe2Exy3m110884t1 1 292;
				## aalen=51;cov=93.0;indels=3/5;match=278;nexon=1;pid=93.6;qlen=314;cdsindel=-2;
				## aalen=63,61%,complete;offs=17-208;oid=Fungr1EG3m041003t1
				## fix2: 900/153061 with chimera=breakpoint.. are missing aalen= from gmap
				## FIXME3: nopath alts should be kept, but have no introns, invalid clen=99 false
				# NOPATH  kf2xx9gmap  mRNA  1   69  .. ID=Funhe2Exx11m004617t1;altc=althi1;trg=Funhe2Exx11m004617t1 1 69;
				#	 cov=0;match=0;pid=1;qlen=99;path=0/0;cdsindel=1866;aalen=803,94%,complete;offs=4-2415;oid=Funhe2E6bm004599t6
				# NOPATH  kf2xx9gmap  mRNA  1   69  .. ID=Funhe2Exx11m004617t2;altc=main;trg=Funhe2Exx11m004617t2 1 69;
				#  cov=0;match=0;pid=1;qlen=99;path=0/0;cdsindel=691;aalen=847,92%,complete;offs=84-2627;oid=Funhe2E6bm004599t1
				
				my($cw,$xw)= m/cxlen=(\d+).(\d+)/; ## BAD expectation .. deal w/ missing field
				unless($xw) { ($xw)=m/(?:qlen|clen)=(\d+)/;  }
				unless($cw) { 
					## my($ofs)=m/offs=([\d-]+)/; if($ofs) { my($ob,$oe)=split/-/,$ofs; $cw=1+$oe-$ob; }
					my($aw,$pcds)=m/aalen=(\d+),(\d+)/; 
					($aw)=m/aalen=(\d+)/ unless($aw); 
					if($aw) { 
						$cw=3*$aw; 
						if(not $xw and ($pcds>1 and $pcds<=100) ) { $xw= int(0.5+ $cw/($pcds/100)); } # recalc
					}
				}
				if($cw and $xw) { $gcdslen{$gidc}="$cw,$xw"; }
				elsif(my $cxw=$gcdslen{$gidc}) { ($cw,$xw)=split",",$cxw; } # 201607 split gene fix? C2 missing aalen annot
				$nerrmrna++ unless($cw and $xw);
				die "ERR: Missing mRNA annots: cxlen= or clen/qlen= and aalen= for $nerrmrna/$ninmrna at ID=$gid in $gf"
					if($nerrmrna/$ninmrna > 0.2);
					
				#? FIXME here: check for alts mapped to diff locus than main ? should remove such from alt test.
					
				if( $cw and $xw ) { 
					my $pcds= ($xw>0) ? int(0.5 + 100*$cw/$xw) : 70; 
					if($ismain) { $mainpcds{$gidc}= $pcds; }
					else {
						my $mainid= $allalt2main{$gid} || $allalt2main{$gidc}; # NOT gidc ? or both?; using splittag;  FIXME: 2+ mainids possible?
						my $mcds= ($mainid) ? $mainpcds{$mainid} : 70; $mcds||=0; # error?
						$keep=0 if($pcds < $MIN_PCDS or $pcds < $MIN_CDS2MAIN * $mcds); # maybe should filter relative to maingene pcds< YES
						$keep=0 if($dropids and $dropids{$gid}); # _C1/2 in this set, not on gidc
						$keep=1 if($keepids and $keepids{$gidc});
						# record not-keeps IDs/reason
						}
					}
				
			} elsif(/\texon/ and $keep) { # assumes now is gene record sorted
				## for altexonend, may need CDS end points also: $gcds{$pid}="begin,end" ? or exons b,e?
				my($pid)= m/Parent=([^;\s]+)/ or next; 
				my $pidc; ($pid,$pidc)= splitID($pid,$_); #* SPLIT id tags needed here also
	
				my($gref,$b,$e)=(split"\t")[0,3,4];  $nxin++;  $iexon++;
				# my($in)= m/;intr=[^,]+,(N[^;\s]+)/; 
				my($inscore,$in)= m/;intr=([^,]+),(N[^;\s]+)/; ## inscore=123 ; =+33/-22 ; =whatelse?
				$gmain{$pid}++ if($ismain);
			
				## Change here? for altexonend: add exon b,e to each in found .. to the right one ! that we dont know from inID
				## ;intr=+2041/-3,N204822,N204851,N204852,N204853
				# odd intr bug: N25789,N25789 >> N25785; N25785,N25786; N25786,N25787; N25787,N25788; N25788,N25789,N25789
				if($in) { 
					my($intscore)= $inscore=~m/([+-]?\d+)/; $intscore ||= 1;
					my @in=split",",$in; my %din; @in=grep{ !$din{$_}++ } @in; 
					if(@in) { 
						my $inkey = join ",",@in; 
						$gx{$pid}{$b} = $inkey;
						$gcds{$pid}{$inkey}="$b,$e"; 
						$exonends{$inkey}{$b}+= $intscore; # not ++ count,
						$exonends{$inkey}{$e}+= $intscore;  # dont reset exonends for each mRNA
						$nxintr++; 
						}          
				} elsif($keepnointron) { # fake something ?
					#? Fix here for too many trivial alts of mapped tr : maybe exclude end exons (but for NOPATH)
					# .. seems to work, altbest count about 1/2 between -nokeepnoin and last -keepnoin
					## naltbest, nokeep=77131/nopath=0, lastkeep=143763/nopath=15130, thiskeep=102248/nopath=14703
					## FIX-gsplign has lots of false long,middle introns for poormap cases; should not allow those as keepers here..
					## my $localkeepnointron= ($keepnointron and not ($src =~ /$GSPLIGN_SOURCE/))?$keepnointron:0;
					## .. or find other way for this; add some okintron annots to input.gff ?
					
					my $inkey = "NOINT".$gref; ## ."xe$b$e"; # for NOPATH, b,e are all same
					my $intscore=1; my $keepnoin=1; 
					if($gref =~ /^NOPATH/) { $inkey = "NOINT".$pid; } else { $keepnoin=0 if($iexon==1 || $iexon>=$nexon); }
					if($keepnoin) {
					$inkey.= ($iexon==1)? "xb$e" : ($iexon>=$nexon)? "xe$b" : "x$b$e";
					$gx{$pid}{$b} = $inkey;
					$gcds{$pid}{$inkey}="$b,$e"; 
					$exonends{$inkey}{$b}+= $intscore; # not ++ count,
					$exonends{$inkey}{$e}+= $intscore;  # dont reset exonends for each mRNA
					}
				}
	# #intronchains for mainID=Funhe2Exx11m004617t1, mainIn=0 == NOPATH w/ valid alts
	# invalidchain    Funhe2Exx11m004617t1    z2      dropshort,0,    NOINTFunhe2Exx11m004617t1xe169
	# invalidchain    Funhe2Exx11m004617t2    z3      dropshort,0,    NOINTFunhe2Exx11m004617t2xe169
	#...      	
			}
		} close(G);
		
		if($keep and $gid) {  ## FIXME3 here? for NOPATH, no intron alts; other mapped+nointron alts?
			my @xin=();
			if($gx{$gid}) { @xin= map{ $gx{$gid}{$_} } sort{ $a<=>$b } keys %{$gx{$gid}}; }
			# if($keepnointron and not @xin) { }# fake it? but how so we get some useful alts, but not lots o crap		
			$ichain{$gid}= join "; ",@xin if @xin; 
		}
	}
	warn "# intron chains: $ninmrna mRNA input, $nerrmrna errors skipped, $nxintr/$nxin exons w/ introns, from main, alt.gff\n" if $debug;
}




sub matchIntronChains {
	# my($altab)= @_;
	
	## open ATAB is used outside sub matchIntronChains
	warn "# $altab: valid-alts for main=$maingenes from alts=$altgenes\n" if $debug;
	open(ATAB,">$altab") or die "$altab"; # STDOUT?
	print ATAB "#altbest input main=$maingenes, alts=$altgenes\n";

# ** SPLIT ids: with and without _C[12..] tag, ichains need split tag, others not
# $gcds{$pid}{$inkey}="$b,$e"; #? use ichain/inkey instead of eqgene %galt to detect main/alt sets?

# sort keepids{} first so others can drop?
	foreach my $gid (sort keys %galt) {
		my (%did,@did,%inset);
		
		## eqgene filtering causing loss of some alts, end up now in "newlocalt" tab, shouldn't do that.
		## .. keep all eqgene overlaps even small, let this ichain testing classify ?
		
		my @alt= @{$galt{$gid}} or next;
		my $dn1= $ichain{$gid} or next; # FIX here for no introns.
		my @dns= sort{ length($ichain{$b}) <=> length($ichain{$a}) or $a cmp $b } grep{ $ichain{$_} } @alt;
		my $anum=1; 
		my $invalnum=0;
		my $i1= $dn1=~tr/;/;/; 
		my $imincount = int($MIN_ALTCHAINSIZE*$i1); # if long = 10 introns, < 3 is too short; can be 0 for keepnointron 
		$imincount||=1 unless($keepnointron);
		$ntabmain++;
		%inset=(); 
		map{ $inset{$_}++ } split /\W+/,$dn1; ## add gid to inset to track old/new loci?
		#? map{ $inset{$_}= $gid } split /\W+/,$dn1; ## add gid to inset to track old/new loci?
		
		# resort by dns cds size? prevent stackup of nearlysame-cds off by few bases at ends
		my (%cdsw,%didcds); 
		#no help??
		#n for my $ad (@dns) { if(my $cxw=$gcdslen{$ad}) { my($cw,$xw)=split",",$cxw; $cdsw{$ad}=$cw; } else { $cdsw{$ad}=1; } }
		#n @dns= sort{ length($ichain{$b}) <=> length($ichain{$a}) or $cdsw{$b} <=> $cdsw{$a} or $a cmp $b } @dns;
		
		print ATAB "#intronchains for mainID=$gid, mainIn=$i1\n"; # needs to be parsable to link altids 
		# sort keepids{} first so others can drop?
		
		foreach my $aid ($gid, @dns) { 
			my $ichain=$ichain{$aid}; my $inum= $ichain=~tr/;/;/; 
			my $ismain= ($aid eq $gid)?1:0;
			
			#?? do gid, @dns have _C split added if need?
			my($aidXXX,$aidc,$asplit)=splitID($aid,"");
			# my $aidc= $aid; $aidc=~s/_C\d+$//; # splitgene tag
			
			my($aidm,$aidi)= ($aid=~m/(\w+)t(\d+)/)?($1,$2):($aid,1); # _Csplit now
			my $altsameloc= ($gid =~ /$aidm/)?1:0;
			next if($didid{$aid}++); # and not ismain
			
			## oknewlocus: handle like new gid here?  need to reset dn1, imincount, inset ..
			## gets too messy .. should separate out newloci from alt.gff before this script.
			# if(0) {    
			#     my $gmain= $gid;
			#     my %ins=(); map{ $ins{$_}++ } split /\W+/,$ichain; 
			#     my %imain=(); map{ my $imain=$inset{$_}; $imain{$imain}++ if($imain); } keys %ins;
			#     my($imainid)= sort{ $imain{$b} <=> $imain{$a} } keys %imain; # gmain replaces gid?
			#    	if($gmain ne $imainid) { # $mainset{$imainid} and 
			#    		$ismain=1; $gmain= $imainid; $dn1= $ichain; 
			#  			$imincount = int($MIN_ALTCHAINSIZE*$inum); # if long = 10 introns, < 3 is too short
			# 		}
			# }		
			
			## FIXME: here to keep alternate end exons: intrchain may be same as others, or may be short,
			##   but some are valid alts, esp. when end exon is longish.
			## need examples to check params. record ic-invalids ? check for alt-drops that have large cds/aalen
			
			my $icsameasmain= (!$ismain and $ichain eq $dn1)?1:0;
				#^ does icsameasmain imply icvalid=-1/-4 ? set  $did{$dn1}=1 ? 
				#^ ah, for aid in (gid,dns) does this.
				
			## these lowqual intronchain tests may need adjust: $inum < int(0.3*$i1) or ichain subset in @did
			my $icvalid = 1; 
			$icvalid=-1 if($icvalid==1 && $did{$ichain});  # ichain is already done
			$icvalid=-2 if($icvalid==1 && ($inum < $imincount));  # ichain is too short ?? is this good test? drop: $inum==0 for keepnointron
			$icvalid=-3 if($icvalid==1 && @did and scalar( grep{ $_ =~ m/$ichain/ } @did )); # ichain is subset of valids
	
			## NOPATH,nointron alts get dropped here: inum==0 (only 1 exon), icvalid == -1, all alts have diff trID ichain tho
			## eg: invalidchain    Funhe2Exx11m004617t1    z2      dropshort,0,    NOINTFunhe2Exx11m004617t1xe169
			
			## FIXME: non-alt diff loci here, check for common introns.  dropshort/-2 wrong unless shared introns
			## ? do this test for all? no, but icvalid==1 can use it, 
			## mark new locus if alt ichain large and completely diff from main ichain
			if(($icvalid == -2 || $icvalid == 1) and $inum>0 and not $ismain) {
				my %ins=(); map{ $ins{$_}++ } split /\W+/,$ichain; 
				my $icommon=0; map{ $icommon++ if($inset{$_}) } keys %ins;
				## add gid to inset to track old/new loci?
				
				$icvalid=3 unless($icommon>0); # icvalid=3/oknewlocus ??
			}
	
=item newloc bugs
		
		* need to handle Split locus genes, somehow.  Use alt-majority locus? ie which intron(s) are most shared?
		
		a. several alts of one loc each become newlocs, should refresh loc-ichain list so those become part of newloc
		
		b. other-loc alt becomes newloc, should not do that shuffle: 1.move other-loc alt to this, 2. move out again
	
			>> alts should become alts of 1 newloc? not quite, these are tandem dups, some splitloc,
				loc1: Zeamay4gEVm000269t3     100.0   98.5    8/8.3   NC_024466.1:160312806-160319657:+
						00269t2,3,6,7,8,9,10,12,13,14 and 0267t2,6,12,25,
				loc2: Zeamay4gEVm000267t4     100.0   99.9    9/9.9   NC_024466.1:160464167-160470957:+
						00267t4,5,7,9,10,11,14,15,16  and 0269t5,11,16
						
		Zeamay4gEVm000267t1newt1        Zeamay4gEVm000267t10    Zeamay4EVm000243t13,cornhi8m9agvelvk67Loc7476t6
		Zeamay4gEVm000267t1newt1        Zeamay4gEVm000267t12    Zeamay4EVm000243t15,cornhi8m9agvelvk67Loc7476t9
		Zeamay4gEVm000267t1newt1        Zeamay4gEVm000267t14    Zeamay4EVm000243t17,cornhi8mtrinLocDN44414c0g3t2
		Zeamay4gEVm000267t1newt1        Zeamay4gEVm000267t15    Zeamay4EVm000243t18,cornhi8m9agvelvk67Loc7476t4
		Zeamay4gEVm000267t1newt1        Zeamay4gEVm000267t16    Zeamay4EVm000243t20,cornhi8m9agvelvk67Loc7476t3
		Zeamay4gEVm000267t1newt1        Zeamay4gEVm000267t17    Zeamay4EVm000243t23,cornlo1rs9sgvelvk93Loc7263t1
		Zeamay4gEVm000267t1newt1        Zeamay4gEVm000267t20    Zeamay4EVm000243t27,cornhi8ms9skvelvk95Loc9462t1
			>> dont belong to EVm000267t
		Zeamay4gEVm000267t1newt1        Zeamay4gEVm000269t11    Zeamay4EVm000243t31,cornlo2m9slvelvk75Loc1060t2
		Zeamay4gEVm000267t1newt1        Zeamay4gEVm000269t6     Zeamay4EVm000243t19,corn1dn9alvelvk85Loc9197t4
		Zeamay4gEVm000267t1newt1        Zeamay4gEVm000269t7     Zeamay4EVm000243t24,cornhi8m9agvelvk77Loc7976t1
	
			>> dont belong to EVm000286t locus, dont move in/out
			>> bug here is Zeamay4gEVm000286t1 is Split over 0286t and 01392t loci, otherwise alts of two are distinct
		Zeamay4gEVm000286t1newt1        Zeamay4gEVm001392t27    Zeamay4EVm000257t76,tidbcornhilor2ridk81Loc28332
		Zeamay4gEVm000286t1newt1        Zeamay4gEVm001392t28    Zeamay4EVm000257t103,cornhi8m9agvelvk35Loc846t30utrorf
		Zeamay4gEVm000286t1newt1        Zeamay4gEVm001392t30    Zeamay4EVm000257t93,cornhi12me4msoapk31loc13388t1
		Zeamay4gEVm000286t1newt1        Zeamay4gEVm001401t7     Zeamay4EVm000257t127,corn1dn4msoapk31loc20649t3
		
=cut
			
	my $debugval="";
	use constant TEST_ALTEND => 1;
	if(TEST_ALTEND) {
			# FIXME? get some dupl subchains -3? No, looks like 2 same-altends w/ diff ichains .. should be -1 did ichain
			if($icvalid == -3 and $gcds{$aid}) { 
				## gcds == geneends : cds ends or exon ends?
				# my($cb,$ce)=split",",$gcds{$aid}; ## do what now?
	
				my $validaltend=0;
				my($icb,$ice)=(split"; ",$ichain)[0,-1];  # front, back intron chain ids
				# test both icb,ice end points of ichain intron chain for exon longer than common size for internal exon.
				for my $ici ($icb,$ice) { 
					next unless($exonends{$ici});
					
					my @ends= sort{$a<=>$b} keys %{$exonends{$ici}}; 
					next unless(@ends >=2);
					# exonends are counted, most freq = common inner exons << No, sum of intr=score gets common inner exons
					# this max test fails unless larger ichain has several alts to count .. count not # alts, but intr=score
					my($end1max,$end2max)= sort{ $exonends{$ici}{$b} <=> $exonends{$ici}{$a} or $a <=> $b } @ends;			  
					($end1max,$end2max)= ($end2max,$end1max) if($end1max > $end2max);
					my($end1left,$end2right)= @ends[0,-1]; # sorted @ends
					
					my($end1this,$end2this)= split ",", $gcds{$aid}{$ici}; # == $exonendsid{$aid}{$ici};
					$end1this||=$end1max; $end2this||=$end2max;
					if($end1this < $end1max) { 
						#n if($end1this == $end1left) { $validaltend++; $exonends{$ici}{$end1this}=$exonends{$ici}{$end1max}; } # RESET max; exon left-end longer; keep? only if end1this == end1left ?
						$validaltend++ if($end1this == $end1left); # exon left-end longer; keep? only if end1this == end1left ?
						} 
					if($end2this > $end2max) { 
						 #n if($end2this == $end2right) { $validaltend++; $exonends{$ici}{$end2this}=$exonends{$ici}{$end2max}; }  #? exon right-end longer; only if end2this == end2right?
						 $validaltend++ if($end2this == $end2right);  #? exon right-end longer; only if end2this == end2right?
						}
			
					$debugval.="$ici:$end1this/$end1max/$end1left-$end2this/$end2max/$end2right,"
						if($debug);
				## need some debug output here; why reject this one?
				# invalidchain    Funhe2Exx6m124896t6     -3,6    N72309,N72313;...
				## but keep this:
				# intrchain       Funhe2Exx6m124896t12    2,4     N72313,N72314,N72315;...
						
				}
				
				#upd17: icvalid=2 isnt problem, icvalid=4 from keepids is too many .. below, slim down keepids
				#... use 2 pass? all keepids, then for each locus, pick only top scored alts to keep
				$icvalid=2 if($validaltend); #? too many here
				#n if($validaltend) { if($didcds{$ichain}) { } else { $icvalid=2; $didcds{$ichain}= $aid; } } #?? is this right?
				
				# test2 icvalid counts:		 # report these?
				# 50071 -1, 7576 -2, 14550 -3, 59776 1, 9642 2  << retained 10k as altend .. but right ones?
	
				## retain Only if cb,ce exon is longer than bigger chain exon; shorter is trash.
				## if($cb < exon-end[i] in @mainic or $ce > exon-end[j] in @mainic) {} have new alt end
				## then icvalid=1 ; modify ichain any as flag?
				##  maybe solution is add exon end locs to each intron chain string? should be all same for inner exons
				
				## append to intron chain as new end pts? or use simple flag for alt exon end?
				## need to know if exon end differs from larger intron chain exons
			}
	}
			
			if($icvalid<=0 and $keepids) {
				$icvalid=4 if($keepids{$aidc});
			}
			
			# FIXME 1607: if $icvalid<0 AND altid not mainid (diff locus), then also test diffloc main (alt t2..t9 only)
			if($icvalid<=0) {
				#above now# my($aidm,$aidi)= $aid=~m/(\w+)t(\d+)$/;
				# altsameloc= ($gid =~ /$aidm/)
				if($aidi>1 and $aidi<9 and not $altsameloc) {
					my $gido= $aidm."t1"; # BUT do only for gido higher number, not lower (larger)?
					my $icvalo=0;
					if(my $ichaino= $ichain{$gido}) {
						$icvalo=-1 if($ichaino eq $dn1 or $ichaino eq $ichain); # icsameasmain or thisalt
						if($icvalo<0 and not ($gmaindup{$gido} or $gmaindupid{$gido} or $gmaindupid{$gid})) {
						push @{$gmaindup{$gid}},$gido; 
						$gmaindupid{$gido}=$icvalo;
						my $valclass= $valclass{$icvalo} || $icvalo;
						# $debugval.="dupof:$gid,";
						my $inumo= $ichaino=~tr/;/;/; 
						print ATAB join("\t","invalidchainm",$gido,'z'.($invalnum+$anum),"$valclass,$inumo,dupof:$gid,",$ichaino)."\n"
							if($debug);  
						}
					}
				}
			}
			if($ismain and $icvalid == 1 and $gmaindupid{$gid}) {
				$icvalid= $gmaindupid{$gid}; # make invalid, from above..
			}
				
			# galtsame opt: collect alt == main set, ichain eq dn1 //did{ichain}, ie same intron chain as maintr, for asmrna check of main    
			push @{$galtsame{$gid}},$aid if($icsameasmain); ## ($aid ne $gid and $ichain eq $dn1);
			
			# unless($did{$ichain} or $inum==0 or $inum < int(0.3*$i1) or scalar(grep{ $_ =~ m/$ichain/ } @did)) 
			## add new alt anum to ATAB
			##? output icvalid as %lclass: 'okay$icv', 'drop$icv' instead of number?
			my $valclass= $valclass{$icvalid} || $icvalid;
			$nclass{$icvalid}++; $ntabalt++ unless($ismain); # report
			if($icvalid>0) {  
				push @did, $ichain unless($did{$ichain}); 
				$did{$ichain}++;  
				## Maybe add valid ichain introns to inset{} ; maybe not: throws off above test
				#? if($icvalid == 3) { } # new gid?
				#  else { map{ $inset{$_}= $gid } split /\W+/,$ichain; }

				if($ismain) { 
					$galtok{$gid}= [] unless($galtok{$gid}); # UPD1802: add valid-gid no alts, ismain ok
					# anum == 1 should be here
				} else { # unless($ismain)
					push @{$galtok{$gid}},$aid; $alt2gid{$aid} = $gid; $altnum{$aid}= ++$anum; 
					
					# upd17 bug: saving/printing many alts of same ichain/cds equal.. need equivalent of main2alt test for alts
					# in alt2gid .. comes from galtsame ? ie alts same chain as main ?
					
					if($icvalid == 3) { $alt2gid{$aid} = "$gid;newlocus=1" ; } # oknewlocus handle how? 
					## may be several alts of gid go to new but same locus, 
					## what of newlocus2, newlocus3? leave that? or try to guess from intron equivs
				} 
				print ATAB join("\t","intrchain",$aid,$TAG_ALT.$anum,"$valclass,$inum,$debugval",$ichain)."\n"; # this is long table, changed cols.

			} elsif(1) {  # debug or always?
				++$invalnum; # track  invalid alt num?
				#? keep $icvalid qual for report?
				$altbadchain{$aid}= $gid; ## keep hash equiv of 	push @{$galtok{$gid}},$aid; $alt2gid{$aid} = $gid;
				#? push @{$galtbad{$gid}},$aid;
				
				print ATAB join("\t","invalidchain",$aid,'z'.($invalnum+$anum),"$valclass,$inum,$debugval",$ichain)."\n"; # this is long table, changed cols.
			}
			# $didid{$aid}++ above
		 }
	}

	# sub putAltTables {} // part of matchIntronChains
	
	my $nclass= join", ", map{ $valclass{$_}.'='.$nclass{$_} } sort{ $b <=> $a } keys %nclass;
	my $nmdup= scalar(keys %gmaindupid);
	$nclass.= ", maindup=$nmdup" if($nmdup);
	my $classrep= "# sum.intrchains classified $nclass for $ntabmain mains, $ntabalt alts as to $altab\n";
	warn $classrep if $debug;
	print ATAB $classrep;
	
	# UPD1802: add valid-gid no alts, here in altsof list? or separate? ; leave out gmaindup
	print ATAB "#valid-alts for main=$maingenes, alts=$altgenes\n";
	foreach my $gid (sort keys %galtok) { 
		my @v= @{$galtok{$gid}}; @v=("none") unless(@v>0);
		print ATAB "altsof\t$gid\t",join(",",@v),"\n";
	}  
	
	# upd main equiv
	print ATAB "#altsame as main=$maingenes, alts=$altgenes\n";
	foreach my $gid (sort keys %galtsame) { 
		my @v= @{$galtsame{$gid}};
		print ATAB "altsame\t$gid\t",join(",",@v),"\n";
	}  
	
	# upd main equiv
	print ATAB "#maindup of main=$maingenes, alts=$altgenes\n";
	foreach my $gid (sort keys %gmaindup) { 
		my @v= @{$gmaindup{$gid}};
		print ATAB "maindup\t$gid\t",join(",",@v),"\n";
	}  
	
	close(ATAB);
}

sub mainFragOfOtherGene {

	# are these maindup? or like those?
	# problem: this finds 2 valid mains joined by joingene alt ..
	# need to use eqgene $allaltseen{$_}= $gid; instead of galt ? to get weak aligns; but that skips dup gid aligns
	
	my (%fragmain);
	foreach my $gid (sort keys %galt) {
		my @altids= @{$galt{$gid}};
		my($gidXXX,$gidc,$gsplit)=splitID($gid,"");
 		for my $aid (@altids) {
 			my $gmain= $alt2gid{$aid} or next; 
			my($gmXXX,$gmainc,$gsplit)=splitID($gmain,"");
			$gmainc =~ s/;newlocus.*//;
 			if($gmainc ne $gidc) {
 				my $gmainval= $alteqval{$aid}{$gmain}||$alteqval{$aid}{$gmainc}||0; # need eqgene overlap vals
 				my $gidval= $alteqval{$aid}{$gid}||$alteqval{$aid}{$gidc}||0;
 				if($gidval > $gmainval) { #?? which way <>? both? > means gid has > align to alt, smaller than gmain?
 					$fragmain{$gid}{$aid}="$gmainc/$gmainval"; #?
 				}
 			}
 		}
	}
	
	my @fragmain= sort keys %fragmain; my $nfg=@fragmain;
	warn "# mainfrag: $nfg fragment mains in main=$maingenes from alts=$altgenes\n" if $debug;
	open(ATAB,">$altab.mainfrag") or die "$altab";  
	print ATAB "#mainfrags n=$nfg in main=$maingenes, alts=$altgenes\n";
	foreach my $gid (@fragmain) { 
		my @alts= sort keys %{$fragmain{$gid}};
		my @agval= map{ "$_/$fragmain{$gid}{$_}" } @alts;
		print ATAB "mainfrag\t$gid\t",join(",",@agval),"\n";
	}  
	
	close(ATAB);
}

# sub newIDformat { sprintf "$IDPREFIX%06d",$_[0]; }

=item newloci..

	case1: 
	newlocus	newl1	2
	nintrchain	Zeamay4gEVm000001t1	newl1	t1	0	N103833; N103833,N103834; N103834; N103835; ..
	nintrchain	Zeamay4gEVm000001t3	newl1	t2	0	N103833; N103833,N103834; N103834; N103835; ..
	
	evg4cornok2reb.main2alt.eqgene     
	Zeamay4gEVm000001t1	
		Zeamay4gEVm000001t2/72.69, << validalt
		Zeamay4gEVm000001t3/27.28	 << newloc because overlap below criteria, but otherwise valid altof t1
		NC_024464.1:154250538-154294086

	evg4corn_tgok2oid.map.attr     
	Zeamay4gEVm000001t1	100.0	99.9	64/65.62	NC_024464.1:154250538-154294086:+	0	Zeamay4EVm000001t1
	Zeamay4gEVm000001t2	100.0	99.3	45/45.45	NC_024464.1:154268669-154294090:+	0	Zeamay4EVm000001t2
	Zeamay4gEVm000001t3	52.0	97.9	20/21.18	NC_024464.1:154250975-154268670:+	0	Zeamay4EVm000001t3

	case2: 	
		.. should be regular valid locus/alt?
		.. both t2,t3 should be listed in altbest.tab, why is t2 not? .. look same-ish, but t3 is 50%cov
	newlocus	newl9998	2
	nintrchain	Zeamay4gEVm005637t1	newl9998	t1	0	N14355; N14355,N14356; N14356,N14357; N14357,N14358; N14358
	nintrchain	Zeamay4gEVm005637t2	newl9998	t2	0	N14356; N14356,N14357; N14357,N14358; N14358
					^^ t2 is like t3? dropsubchain? need that test for newloc sub?
	eqgene:
	Zeamay4gEVm005637t1	Zeamay4gEVm005637t3/46.46,Zeamay4gEVm005637t2/45.32	NC_024459.1:277206217-277217949
	altbest.tab:
	#intronchains for mainID=Zeamay4gEVm005637t1, mainIn=4
	intrchain	Zeamay4gEVm005637t1	t1	okay,4,	N14355; N14355,N14356; N14356,N14357; N14357,N14358; N14358
	invalidchain	Zeamay4gEVm005637t3	z2	dropsubchain,3,..	N14356; N14356,N14357; N14357,N14358; N14358
	map.attr:
	Zeamay4gEVm005637t1	100.0	97.4	5/5.4	NC_024459.1:277206217-277217949:+	0	Zeamay4EVm004942t1
	Zeamay4gEVm005637t2	80.3	100.0	4/4.4	NC_024459.1:277207790-277217949:+	0	Zeamay4EVm004942t3
	Zeamay4gEVm005637t3	53.2	100.0	4/4.4	NC_024459.1:277207815-277217615:+	0	Zeamay4EVm004942t2


=cut

sub newlociFromIchains # test 1607
{
  #	my($lnewtab)= @_;

	# $ichain{$mrnaid}= join "; ",@xin if @xin; 
	# $gmain= $alt2gid{$tid} || ""; ignore these good alts but for "newlocus"
	# $altbadchain{$aid}= $gid; # skip these bad alts
	
	warn "# new loci/alts to $lnewtab\n" if $debug;
	open(NTAB,">$lnewtab") or die "$lnewtab"; # STDOUT?
	print NTAB "#newloci from main=$maingenes, alts=$altgenes\n";

	my @allid= sort keys %ichain;
	# skip also main.gff IDs, from $gmain{$pid}, these all presumably have been tested against alt chains
	# 	collectIntronChains:	$gmain{$pid}++ if($ismain);
  # #? skip allaltseen from eqgene? maybe reduce validalt overlap cutoff?
  #  	my $altineqtab= $allaltseen{$gid} || $allaltseen{$gidc} || 0; ## $gid,$gidc

	my @testid= grep{ my $mid= $alt2gid{$_}||0; $mid=0 if($mid=~/newlocus/); 
							  not ($mid or $altbadchain{$_}) } @allid;
	# or $gmain{$_}  << keep gmains, list only when have altof
	
	my (%alt2new,%newlocids,%inset,%idid,@idid,@invalid);
	my ($invalnum)= (0);
	## should record newloc ids to inset
	my $NEWLOCID= 0;
	my $MINCOMINT= $ENV{mincomint}||2; # or 3?
	my $SKIPNEW1EXON= not ($ENV{keepnew1exon}||0);
	my $imincount = ($SKIPNEW1EXON) ? 1 : 0; # int($MIN_ALTCHAINSIZE * $i1); # if long = 10 introns, < 3 is too short; can be 0 for keepnointron 
	
	for my $aid (@testid) {
		my $ichain= $ichain{$aid}; # look for common subchains to collect new loci
		my ($nida,$icommon)= (0,0);
		my (%ins,%nids);
		my($aidXXX,$aidc,$asplit)=splitID($aid,"");
		
		# my $wasbad= $altbadchain{$aid} || 0; # presume bad alt chains, not nec new loci
		my $wasnewloc= ($alt2gid{$aid} and $alt2gid{$aid}=~/newlocus/)?1:0;
		my $wasovercds= (not $wasnewloc and ($altbadchain{$aid} or $allaltseen{$aid}))?1:0; ## $gid,$gidc
		
		my $ismain= $gmain{$aid}; # keep here but list only with altof
		my $inum= $ichain=~tr/;/;/; 
		map{ $ins{$_}++ } split /\W+/,$ichain; # per intron test, but want >1 match for valid subchain equiv

		if($ismain) {
			$imincount = int($MIN_ALTCHAINSIZE * $inum); 
			$imincount = 1 if($SKIPNEW1EXON and $imincount<1); # unless($keepnointron);
		}

use constant ICVALID_NEWLOC => 1;	
		## add other inchain qual tests as above?
		my $icvalid = 1; 
		if(ICVALID_NEWLOC and not $ismain) {
			$icvalid=-1 if($icvalid==1 && $idid{$ichain});  # ichain is already done
			$icvalid=-2 if($icvalid==1 && ($inum < $imincount));  # ichain is too short ?? is this good test? drop: $inum==0 for keepnointron
			$icvalid=-3 if($icvalid==1 && @idid and scalar( grep{ $_ =~ m/$ichain/ } @idid )); # ichain is subset of valids
			if(($icvalid == -2 || $icvalid == 1) and $inum>0 and not $ismain) {
				#have# my %ins=(); map{ $ins{$_}++ } split /\W+/,$ichain; 
				my $icommon=0; map{ $icommon++ if($inset{$_}) } keys %ins;
				$icvalid=3 unless($icommon>0); # icvalid=3/oknewlocus ??
			}
		}
		
		if($icvalid>0 and $inum<1 and $SKIPNEW1EXON) {  # dont need if using imincount test above; ?? want inum<1,  ";"  means 2 exons/1 intron
			$icvalid= -2;
		}
		if($icvalid<=0 and $keepids) {
			$icvalid=4 if($keepids{$aidc});
		}
		
		if($icvalid > 0 and $wasovercds) { $icvalid= 0; } #upd17
			
		# 1 exon loci a problem still, keep here or skip?
		#? report $aid as invalidchain/oneexon if debug?
		if($icvalid < 1) {
			++$invalnum; # track  invalid alt num?
			my $valclass= $valclass{$icvalid} || $icvalid;
			push @invalid, ["invalidchain",$aid,'z'.$invalnum,"$valclass,$inum",$ichain];
			# print NTAB join("\t","invalidchain",$aid,'z'.($invalnum+$anum),"$valclass,$inum",$ichain)."\n";   
			# next;
		} 

		if($icvalid > 0) {
			for my $k (sort keys %ins) { if(my $nid= $inset{$k}) { $nids{$nid}++; $icommon++; } } ;
	
			if($icommon>=$MINCOMINT) {
				($nida)= sort{ $nids{$b}<=>$nids{$a} or $a cmp $b } keys %nids;
			} else {
				$nida= ++$NEWLOCID;
			}
			
			my $gidnew= $IDPREFIX.sprintf "%06d",$nida;
			for my $k (keys %ins) { $inset{$k}=$nida unless($inset{$k}); } # 2+ newloc/intron?
			push @{$newlocids{$gidnew}},$aid;
			$alt2new{$aid}= $gidnew;
			push @idid, $ichain unless($idid{$ichain}); $idid{$ichain}++;
		}
	}
	
	my $nalt=0; # scalar(keys %alt2new);
	my $nloci=0; #scalar(keys %newlocids);

	# foreach newloc, alts, list ichains?
	for my $gidnew (sort keys %newlocids) {
		my @aid= @{$newlocids{$gidnew}}; # sort by ichain size $inum*
		@aid= sort{ length($ichain{$b}) <=> length($ichain{$a}) or $a cmp $b } @aid;
		my($gmain)= grep{ $gmain{$_} } @aid; $gmain||=""; # keep here but list only with altof
		if($gmain) { @aid= grep{ $_ ne $gmain } @aid; next unless(@aid); }
		
  	print NTAB join("\t","newlocus",$gidnew,($gmain||"nomain"),scalar(@aid),)."\n"; 
  	my $anum=0; $nloci++;
  	unshift(@aid,$gmain) if($gmain);
		for my $aid (@aid) {
			my $ismain= ($gmain and $aid eq $gmain)?1:0; # keep here but list only with altof
			my $dni= $ichain{$aid}||"noin"; # look for common subchains to collect new loci
			my $inum= $dni=~tr/;/;/; 
			++$anum; $nalt++ unless($ismain);
			my $valclass=($ismain)?"origmain":"oknewlocalt";
  		print NTAB join("\t","nintrchain",$aid,$gidnew,$TAG_ALT.$anum,"$valclass,$inum",$dni)."\n"; 

			## upd17 for putAlt, reset gid, alt2gid, ..
			$alt2gid{$aid} = $gidnew; 
			$altnum{$aid}  = $anum; 

  		}
  	}

=item outp
newlocus        newlocA005415   nomain  6
nintrchain      Daplx5cEVm013894t49     newlocA005415   t1      oknewlocalt,4   N168427; N168427,N168428; N168428,N168429; N168429; N168430
nintrchain      Daplx6suEVm010416t74    newlocA005415   t2      oknewlocalt,3   N168427; N168427,N168428; N168428,N168429; N168422,N168429
nintrchain      Daplx5cEVm013894t54     newlocA005415   t3      oknewlocalt,4   N168423; N168423; N168428; N168428,N168429; N168429
=cut

  if($debug) { # all at end?
  	for my $inv (@invalid) {
  		my $aid= $inv->[1]; $alt2gid{$aid}= 0; #?? dont putAlt
 			print NTAB join("\t",@$inv)."\n"; #? "invalidchain",$aid,'z'.($invalnum+$anum),"$valclass,$inum",$ichain)."\n";   
 		}
  }	
  
  close(NTAB);

	warn "# newlociFromIchains nloci=$nloci, nalt=$nalt\n" if $debug;
	return($nalt,$nloci);
}



=item ic -3 case

	#intron chains of alts for Funhe2Exx6m124896t1, inmain=10
	intrchain       Funhe2Exx6m124896t1     1,10,   N72304,N72305;...
	intrchain       Funhe2Exx6m124896t9     1,11,   N72303,N72304,N72305;...
	intrchain       Funhe2Exx6m124896t7     1,10,   N72303,N72304,N72305;...
	intrchain       Funhe2Exx6m124896t4     1,9,    N72304,N72305;...
	intrchain       Funhe2Exx6m124896t2     1,8,    N72295,N72306,N72307,N72308;...
	invalidchain    Funhe2Exx6m124896t14    -3,7,N72295,N72306,N72307,N72308:235875/235875/235551-236175/236175/236175,
		N72318:244848/244848/244848-245291/245291/245291, N72295,N72306,N72307,N72308;...
	intrchain       Funhe2Exx6m124896t3     1,9,    N72304,N72305;...
	intrchain       Funhe2Exx6m124896t8     1,7,    N72304,N72305;...
	intrchain       Funhe2Exx6m124896t5     1,8,    N72304,N72305;...
	intrchain       Funhe2Exx6m124896t15    1,6,    N72295,N72306,N72307,N72308;...
	invalidchain    Funhe2Exx6m124896t16    -1,6,   N72295,N72306,N72307,N72308;...
	invalidchain    Funhe2Exx6m124896t17    -1,6,   N72295,N72306,N72307,N72308;...
	
	** gmap/gsplign mapping differs for t6, ichains differ & affect this test; above is all? gspl, later pub9set has gmap t6
	>> t6 end 237589 isn't longest, t7 is
	invalidchain    Funhe2Exx6m124896t6     -3,6,N72309,N72313:237589/237694/237549-237813/237813/237813,N72319:246562/24
	6562/245302-246834/246834/246834,       N72309,N72313;...
	
	invalidchain    Funhe2Exx6m124896t10    -3,5,N72309,N72313:237589/237694/237549-237813/237813/237813,N72319:246562/24
	6562/245302-246834/246834/246834,       N72309,N72313;...
	
	t7/Funhe2Eq7m077478t4/Funhe2Exx9m124883t16 longest exon at N72309,N72313: .. but doesnt end there. 
	grep 'Parent=Funhe2Exx6m124896t' kf2mixx_dedup.analt.gff | grep exon | grep 237549
	Scaffold10150   splkf2eg37mxx   exon    237549  237813  0.978   -       .       Parent=Funhe2Exx6m124896t7;trg=Funhe2Exx6m124896t7 1374 1641;splice=AGTC-
		;intr=+77/-46,N72309,N72313
	
	t6/Funhe2E6bm008467t4/Funhe2Exx9m124883t6
	grep 'Parent=Funhe2Exx6m124896t' kf2mixx_dedup.analt.gff | grep exon | grep 237589
	Scaffold10150   splkf2eg37mxx   exon    237589  237813  0.982   -       .       Parent=Funhe2Exx6m124896t6;trg=Funhe2Exx6m124896t6 1344 1567;splice=AGAG-
		;intr=+77/-46,N72309,N72313
	Scaffold10150   kfish2evg367    exon    237589  237813  98      -       .       Parent=Funhe2Exx6m124896t10;trg=Funhe2Exx6m124896t10 1224 1447;ix=6
		;intr=+77/-46,N72309,N72313

=item icvalid/invalid work

    ## ^ic -3 is commonest invalid, also where alt exons likely.
    ## kf2mixx count: 45270 -1(did) 7331 -2(short) 29196 -3(subset)
    ## eg: invalidchain    Funhe2Exx6m124896t6/Funhe2E6bm008467t4/Funhe2Exx9m124883t6     -3      N72309,N72313;...
		## may well be case of valid alt 3' exon end, longer than spliced exon at that spot of larger alts.
		## t1,t2,t4 have same exon w/ splice in middle..

	> t6 3'end exon, 105b longer exon, cds 69b longer
	Scaffold10150   splkf2eg37mxx  exon    237589  237813  0.982   -       .       Parent=Funhe2Exx6m124896t6;trg=Funhe2Exx6m124896t6 1344 1567;splice=AGAG-;
		intr=+77/-46,N72309,N72313 << this is same inID set as t2, but -46 means one N72309? maps inside exon ..
	Scaffold10150   splkf2eg37mxx  CDS     237625  237813  1       -       .       Parent=Funhe2Exx6m124896t6
	> t2 same exon
	Scaffold10150   kf2evg367mixy  exon    237694  237813  100     -       .       Parent=Funhe2Exx6m124896t2;trg=Funhe2Exx6m124896t2 1351 1470;ix=8;
		intr=123,N72309,N72313
	Scaffold10150   kf2evg367mixy  CDS     237694  237813  100     -       .       Parent=Funhe2Exx6m124896t2
	Scaffold10150   kf2evg367mixy  exon    237420  237539  100     -       .       Parent=Funhe2Exx6m124896t2;trg=Funhe2Exx6m124896t2 1471 1590;ix=9;intr=386,N72308,N72309,N72310,N72311,N72312
	
	grep Funhe2Exx6m124896t  kf2mixx_dedup.altbest.tab | sed 's/;.*/;.../;' | less
	#intron chains of alts for Funhe2Exx6m124896t1
	intrchain       Funhe2Exx6m124896t1     N72304,N72305;...
			^^ exon 234272  234420 ; intr=+113/-3,N72304,N72305 3'end exon;..
	intrchain       Funhe2Exx6m124896t9     N72303,N72304,N72305;...
	intrchain       Funhe2Exx6m124896t7     N72303,N72304,N72305;...
			^^ exon 234128  234420 ; intr=+113/-38,N72303,N72304,N72305 : 150b longer 3'end, hits new intron
	intrchain       Funhe2Exx6m124896t4     N72304,N72305;...
	intrchain       Funhe2Exx6m124896t2     N72295,N72306,N72307,N72308;...
	invalidchain    Funhe2Exx6m124896t14    -3      N72295,N72306,N72307,N72308;...
	intrchain       Funhe2Exx6m124896t3     N72304,N72305;...
	intrchain       Funhe2Exx6m124896t8     N72304,N72305;...
	intrchain       Funhe2Exx6m124896t5     N72304,N72305;...
	intrchain       Funhe2Exx6m124896t15    N72295,N72306,N72307,N72308;...
	invalidchain    Funhe2Exx6m124896t16    -1      N72295,N72306,N72307,N72308;...
	invalidchain    Funhe2Exx6m124896t17    -1      N72295,N72306,N72307,N72308;...
	invalidchain    Funhe2Exx6m124896t6     -3      N72309,N72313;...							<<* likely valid altend
	invalidchain    Funhe2Exx6m124896t10    -3      N72309,N72313;...
	invalidchain    Funhe2Exx6m124896t12    -3      N72313,N72314,N72315;...
	altsof  Funhe2Exx6m124896t1     Funhe2Exx6m124896t9,Funhe2Exx6m124896t7,Funhe2Exx6m124896t4,Funhe2Exx6m124896t2,Funhe2
	Exx6m124896t3,Funhe2Exx6m124896t8,Funhe2Exx6m124896t5,Funhe2Exx6m124896t15
	
	# start-5': N72318,N72319; N72319 (but for shorter 5' start N72317,N72318; N72318 inside longer : is it valid? see t15)
	intrchain       Funhe2Exx6m124896t1     N72304,N72305; N72305,N72306; N72306,N72307,N72308; N72308,N72309,N72310,N72311,N72312; N72309,N72313; N72313,N72314,N72315; N72311,N72314,N72316; N72312,N72315,N72316,N72317; N72317,N72318; N72318,N72319; N72319
	intrchain       Funhe2Exx6m124896t9     N72303,N72304,N72305; N72305,N72306; N72306,N72307,N72308; N72308,N72309,N72310,N72311,N72312; N72309,N72313; N72313,N72314,N72315; N72311,N72314,N72316; N72312,N72315,N72316,N72317; N72317,N72318; N72318; N72319; N72319
	intrchain       Funhe2Exx6m124896t7     N72303,N72304,N72305; N72305,N72306; N72306,N72307,N72308; N72308; N72309,N72313; N72313,N72314,N72315; N72311,N72314,N72316; N72312,N72315,N72316,N72317; N72317,N72318; N72318,N72319; N72319
	intrchain       Funhe2Exx6m124896t4     N72304,N72305; N72305,N72306; N72306,N72307,N72308; N72308,N72309,N72310,N72311,N72312; N72309,N72313; N72313,N72314,N72315; N72312,N72315,N72316,N72317; N72317,N72318; N72318,N72319; N72319
	intrchain       Funhe2Exx6m124896t2     N72295,N72306,N72307,N72308; N72308,N72309,N72310,N72311,N72312; N72309,N72313; N72313,N72314,N72315; N72311,N72314,N72316; N72312,N72315,N72316,N72317; N72317,N72318; N72318,N72319; N72319
	invalidchain    Funhe2Exx6m124896t14    -3      N72295,N72306,N72307,N72308; N72308,N72309,N72310,N72311,N72312; N72309,N72313; N72313,N72314,N72315; N72311,N72314,N72316; N72312,N72315,N72316,N72317; N72317,N72318; N72318
	intrchain       Funhe2Exx6m124896t3     N72304,N72305; N72305,N72306; N72306,N72307,N72308; N72308; N72313,N72314,N72315; N72311,N72314,N72316; N72312,N72315,N72316,N72317; N72317,N72318; N72318,N72319; N72319
	intrchain       Funhe2Exx6m124896t8     N72304,N72305; N72305,N72306; N72306,N72307,N72308; N72308,N72309,N72310,N72311,N72312; N72312,N72315,N72316,N72317; N72317,N72318; N72318,N72319; N72319
	intrchain       Funhe2Exx6m124896t5     N72304,N72305; N72305,N72306; N72306,N72307,N72308; N72308; N72313,N72314,N72315; N72312,N72315,N72316,N72317; N72317,N72318; N72318,N72319; N72319
	intrchain       Funhe2Exx6m124896t15    N72295,N72306,N72307,N72308; N72308; N72313,N72314,N72315; N72311,N72314,N72316; N72312,N72315,N72316,N72317; N72317,N72318; N72318
	invalidchain    Funhe2Exx6m124896t16    -1      N72295,N72306,N72307,N72308; N72308; N72313,N72314,N72315; N72311,N72314,N72316; N72312,N72315,N72316,N72317; N72317,N72318; N72318
	invalidchain    Funhe2Exx6m124896t17    -1      N72295,N72306,N72307,N72308; N72308; N72313,N72314,N72315; N72311,N72314,N72316; N72312,N72315,N72316,N72317; N72317,N72318; N72318
	invalidchain    Funhe2Exx6m124896t6     -3      N72309,N72313; N72313,N72314,N72315; N72311,N72314,N72316; N72312,N72315,N72316,N72317; N72317,N72318; N72318,N72319; N72319
	invalidchain    Funhe2Exx6m124896t10    -3      N72309,N72313; N72313,N72314,N72315; N72312,N72315,N72316,N72317; N72317,N72318; N72318,N72319; N72319
	invalidchain    Funhe2Exx6m124896t12    -3      N72313,N72314,N72315; N72312,N72315,N72316,N72317; N72317,N72318; N72318,N72319; N72319

=cut
		

# also create altbest.gff from these.. changing alt ID to main.tN...
# ADD option to NOT change alt id.
# ADD option not to write this valtgff .. make sure ATAB is sufficient to regenerate alt.gff
# FIXME: oknewlocus class needs to be flagged, ID changed to non-alt, per CHANGEALTID option ..
#   .. If this is accurate reclass, which seems to be mostly.
# kf2pub11_both.altbest.tab: valid-alts for main=kf2pub11_both.anmain.gff from alts=kf2pub11_both.analt.gff
# sum.intrchains classified oknewlocus=6787, okaltend=15150, okay=63817, dropdup=63739, dropshort=22672, dropsubchain=26996 for 25364 mains, 173797 alts as to kf2pub11_both.altbest.tab
## newloc eg
#intronchains for mainID=Funhe2Exx11m000009t1, mainIn=94
# intrchain       Funhe2Exx11m000009t1    t1      okay,94,        N158164;..
# intrchain       Funhe2Exx11m000009t2    t2      okay,93,        N158164;..
# intrchain       Funhe2Exx11m000009t3    t3      okay,93,        N158164;..
# ..
# intrchain       Funhe2Exx11m000009t26   t16     oknewlocus,45,  N99571,N99577;..
# .. N99571,N99577; N99578,N99579,N99580; N99580,N99583,N99584,N99586,N99587; N99586,N99588; N99587,N99588,N99589; N99589,N99590; N99590; N99592; N99592,N99593; N99593,N99594; N99594,N99595; N99595,N99596; N99596,N99598; N99598,N99599; N99599,N99600; N99600,N99601; N99601,N99602; N99602,N99603; N99603,N99606; N99606,N99607; N99607; N99608,N99609; N99608,N99609,N99610,N99611; N99611,N99612,N99618; N99618,N99619; N99616,N99620; N99617,N99620,N99621,N99622; N99622,N99623; N99623,N99624; N99624,N99625,N99626; N99626,N99629; N99629,N99630; N99630; N99631; N99631,N99632; N99632,N99633; N99633,N99634; N99634,N99635; N99635,N99637; N99637,N99638; N99638,N99639; N99639,N99641; N99641,N99642; N99642,N99643; N99643,N99644; N99644
# invalidchain    Funhe2Exx11m000009t29   z24     dropshort,34,   N158241;..
# intrchain       Funhe2Exx11m000009t25   t17     oknewlocus,39,  N99580,N99583,N99
# intrchain       Funhe2Exx11m000009t28   t18     oknewlocus,39,  N99580,N99583,N99
# invalidchain    Funhe2Exx11m000009t27   z27     dropshort,30,   N158169,N158170,N
# intrchain       Funhe2Exx11m000009t19   t19     oknewlocus,23,  N99631;..
# intrchain       Funhe2Exx11m000009t18   t20     oknewlocus,21,  N99631;..
# invalidchain    Funhe2Exx11m000009t32   z32     dropshort,13,   N158259;..
# intrchain       Funhe2Exx11m000009t30   t21     oknewlocus,10,  N99646;.. << dropshort-newlocus
# intrchain       Funhe2Exx11m000009t35   t22     oknewlocus,1,   N99658;.. << dropshort-newlocus


sub putAltGFF {
	#  my($altgenes, $valtgff, $newalt, $newloci, $lnewgff)= @_;
	my $gf= $altgenes;
	my $gop = ($gf =~ /.gz/) ? "gunzip -c $gf" : "cat $gf"; 
	open(G,"$gop |") or die "bad $gop";
	
	warn "# valid-alts to $valtgff for main=$maingenes from alts=$altgenes\n" if $debug;
	open(AGFF,">$valtgff") or die "$valtgff"; 
	print AGFF "##gff-version 3\n";
	print AGFF "# valid-alts for main=$maingenes from alts=$altgenes\n";

	my $putnewloci= ($newalt > 0 and not $skipnewloci)?1:0;;
	## upd17 write "newlocus" to separate .gff, many are not new loci but mistakes, some maybe ok
	if($putnewloci) { # hasnewloci and not skipnewloci
		open(NTAB,">$lnewgff") or die "$lnewgff";  
		print NTAB "##gff-version 3\n";
		print NTAB "#newloci from main=$newloci, alts=$newalt\n";
	}
	
	my($altid,$isalt,$isnewmain,$isnewalt,$hasnomain,$galt,$gmain)= (0) x 9;
	my($nnewput,$naltput,$nnomainput)=(0) x 9;
	my(%gmains);
	
	while(<G>) {
		next unless(/^\w/);
		if(/\tmRNA/) { 
			my($tid)=m/ID=([^;\s]+)/;
			my($gid,$gidc,$gsplit)= splitID($tid,$_); #* SPLIT id tags needed here also
			$gmain= $alt2gid{$tid} || "";
			$gmains{$gmain}++ if($gmain);
			(my $mainid=$gmain) =~ s/;.*//;  # add to alts to keep proper link
		
			$isalt= ($gmain)? 1: 0;
			## isnewmain/alt should be option still .. not debugged
			
			$isnewalt= ($gmain =~ /^$IDPREFIX/ or $gmain =~ /newloc/)?1:0;  # rename NEWLOCIDPREFIX ..
			$isnewmain= ($isnewalt and $altnum{$tid} < 2)?1:0;  #??$gid =~ /t1$/
			#bad $isnewmain= ($gmain =~ /newloc/)?1:0; #?? is this bad; $putnewloci ??
			$hasnomain=0; # 1607: == not($isalt or $isnewmain), no overlap to main.gff in eqgene.tab
			
			#? my $altineqtab= $allaltseen{$gid} || $allaltseen{$gidc} || 0; ## $gid,$gidc
			
			my $ismain= $gmain{$gid}; # fixme 170110: alt-input.gff can contain all of main-input.gff .. skip
			if($ismain) {
				s/$/;mainid=$mainid/ if($ADDMAINID and $mainid); # upd1802: default? 
				$isalt=$isnewmain=$hasnomain=0; # move on
				# ismain NOT printed below;  should print ismain be option?
								
			} elsif($isalt or $isnewmain) {  # fixme: isalt of isnewmain
				$galt=$tid; 
				($altid=$gmain) =~ s/;.*//;  # change t1 to taltnum, but if gmain == t2..tn then problem.
				if($isnewmain or $isnewalt) { $altid.="new" unless($altid=~/new/); }
				
				my $ai= $altnum{$tid} || 99;  
				my($altids,$altidc,$asplit)= splitID($altid,$_); #* SPLIT id tags needed here also
				if($altidc =~ /${TAG_ALT}1$/) { $altidc =~ s/${TAG_ALT}1$//; } # else problem or not?
				$altidc.= $TAG_ALT.$ai;
				$altids= $altidc; $altids.="_C$asplit" if($asplit);
				$altid= $altids; #split-tagged
				
				my $ctag= ($isnewalt) ? "locus" : "altid";
				if($CHANGEALTID) { s/ID=$tid/ID=$altid;old$ctag=$galt/; } 
				else { s/$/;/ unless(m/;/); s/;/;new$ctag=$altid;/; }
				unless(m/;alttr=1/ or $ai == 1) { s/$/;alttr=1/; } # old tag, drop?
				s/$/;mainid=$mainid/ if($ADDMAINID and $mainid); # upd1802: default?
				if($isnewmain or $isnewalt) { $nnewput++; } else { $naltput++; }
					
			} else {
				#? maybe treat like isnewmain .. give new locus id?
				s/$/;mainid=$mainid/ if($ADDMAINID and $mainid); # upd1802: default?
				s/$/;altnomain=1/;
				my $altineqtab= $allaltseen{$gid} || $allaltseen{$gidc} || 0; ## $gid,$gidc
				$hasnomain=1 unless($dropnomain or $altineqtab); # this is wrong, puts ALL, need flag for NotInEqgene
				$nnomainput++ if($hasnomain);
		 }
		} elsif(/\t/ and $CHANGEALTID and ($isalt||$isnewmain)) { # not $hasnomain 
			s/Parent=$galt/Parent=$altid/;
		} 
		
		# print AGFF $_ if ($isalt||$isnewmain||$hasnomain);
		if($isnewmain||$isnewalt||$hasnomain) { print NTAB $_ if($putnewloci); }
		elsif($isalt) { print AGFF $_; }  ##  
	}
	close(G); close(AGFF);
	close(NTAB) if($putnewloci);
	my $nloc=scalar(keys %gmains);
	my $natot=$naltput+$nnewput+$nnomainput;
	warn "# valid-alts put nloci=$nloc, altvalid=$naltput, altnewloc=$nnewput, altnomain=$nnomainput, alttotal=$natot\n" if $debug;
}

__END__

=item intronchains parsing

/bio/bio-grid/kfish2/rnas/kf2evgr/trevg367mixx/alt11
nam=kf2pub11_both

$evigene/scripts/overlapfilter \
 -intron2splice=error -pass 'exon,intron' -act markid -midtype scoresum -mark intr \
 -over $kfish2/intron/intron_okay.gff.gz -in $nam.gff.gz > $nam.an2.gff

 wc -l $nam.an2.*chains
   33075 kf2pub11_both.an2.galtchains
   11772 kf2pub11_both.an2.gsamechains
  183032 kf2pub11_both.an2.ichains

## intron chains table:

grep intr= $nam.an2.gff | perl -ne \
'($d)=m/Parent=(\w+)/; ($in)=m/intr=([^;\s]+)/;
($iv,$in1,$in2)=split",", $in; $iv=~s,/.*,,; if($iv>0) { map{
$gin{$d}{$_}++; $ic{$_}{$d}++ } grep/\w/, ($in1,$in2); } END{ for $g
(sort keys %gin) { @in=sort keys %{$gin{$g}}; %lg=(); $ni=@in; $niv=0;
for $i (@in) { $niv+=$gin{$g}{$i}; map{ $lg{$_}++; }  keys %{$ic{$i}};
} ($gd=$g)=~s/t\d+$//; $minc=$ni*0.20; @lg= map{ my $c=$lg{$_};
$c.="NO" if($c<$minc);  s/$gd/id_/; "$_:$c"; } grep{ $_ ne $g }
sort{$lg{$b}<=>$lg{$a} or $a cmp $b} keys %lg;
$ng=join",",grep/NO/,@lg; $lg=join",",grep{ not /NO/} @lg; $ng||="na";
$lg||="na"; print join("\t",$g,$ni,$niv,$lg,$ng)."\n";  }}  ' \
  > $nam.an2.ichains

cat $nam.an2.ichains | cut -f1,2,4 | perl -ne \
'($d,$ni,$ad)=split; ($g,$t)= $d=~m/(\w+)t(\d+)$/; @at=map{ s/id_t//;
s/:\w+//;  $_; } grep /id_/,split",",$ad; @ot=sort{$a<=>$b}($t,@at);
$ot=$ot[0]; map{ $sg{$g}{$ot}{$_}++ } @ot; if($lg and $lg ne $g) {
@t=sort keys %{$sg{$lg}}; %did=(); for $t (@t) { @at=sort{$a<=>$b}
keys %{$sg{$lg}{$t}}; $new=1; map{$new=0 if($did{$_}++); }@at;
if($new) { $at=join",",@at; print join("\t",$lg,"t".$t,$at)."\n"; } }
} $lg=$g; ' \
 > $nam.an2.galtchains

cat $nam.an2.ichains | cut -f1,2,4 | perl -ne \
'($d,$ni,$ad)=split; next if($ad eq "na"); ($g,$t)=
$d=~m/(\w+)t(\d+)$/; @at=map{  s/:\w+//;  $_; } grep { not /id_/ }
split",",$ad; next unless(@at); @ot=sort{$a<=>$b}(@at); $ot=$t; map{
$sg{$g}{$ot}{$_}++ } @ot; if($lg and $lg ne $g) { @t=sort keys
%{$sg{$lg}}; %did=(); for $t (@t) { @at=sort keys %{$sg{$lg}{$t}};
$new=1; map{$new=0 if($did{$_}++); }@at; if($new) { $at=join",",@at;
print join("\t",$lg,"t".$t,$at)."\n"; } } } $lg=$g; ' \
> $nam.an2.gsamechains

==> kf2pub11_both.an2.galtchains <==
  n=4855 loci (20%) have diff locus alts (>t1), ng=29135 loci have some altchain, ng1=26624 loci have t1 chain (vs 27492 below)
Funhe2Exx11m000001      t1      1
Funhe2Exx11m000002      t1      1
Funhe2Exx11m000003      t1      1
Funhe2Exx11m000003      t3      3
Funhe2Exx11m000004      t1      1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31
..
Funhe2Exx11m000012      t1      1,2,3,4,5,9,10,13,14,15,16,18  << Scaffold10067
Funhe2Exx11m000012      t17     17   << Scaffold9902
Funhe2Exx11m000012      t6      6,7,8,11,12   << Scaffold9981
  
==> kf2pub11_both.an2.gsamechains <==
  n=9988 loci (33%) have other at same locus, ng1=7312 t1 have same-locus
Funhe2Exx11m000009      t18     Funhe2Exx11m000013t1,Funhe2Exx11m000013t2,Funhe2Exx11m000013t3,Funhe2Exx11m000013t4,Funhe2Exx11m000013t7
   ^^ Scaffold14 09t18,013t (other 09t1.. at Scaffold323)
Funhe2Exx11m000012      t11     Funhe2Exx11m000093t1,Funhe2Exx11m000093t2,Funhe2Exx11m000093t3
Funhe2Exx11m000012      t17     Funhe2Exx11m000027t1,Funhe2Exx11m000027t10,Funhe2Exx11m000027t11,Funhe2Exx11m000027t12,Funhe2Exx11m000027t13,Funhe2Exx11m000027t14,Funhe2Exx11m000027t15,Funhe2Exx11m000027t16,Funhe2Exx11m000027t17,Funhe2Exx11m000027t18,Funhe2Exx11m000027t19,Funhe2Exx11m000027t2,Funhe2Exx11m000027t22,Funhe2Exx11m000027t25,Funhe2Exx11m000027t3,Funhe2Exx11m000027t4,Funhe2Exx11m000027t5,Funhe2Exx11m000027t6,Funhe2Exx11m000027t7,Funhe2Exx11m000027t8,Funhe2Exx11m000027t9
Funhe2Exx11m000013      t1      Funhe2Exx11m000009t18,Funhe2Exx11m000009t19,Funhe2Exx11m000009t25,Funhe2Exx11m000009t26,Funhe2Exx11m000009t28
	gsamechains t1 vs kf2pub11_both.main.eqgene:
	  of samechain t1=7312, noeqgene n=3631; eqcds>=20% n=2190; eqgene > 0 n=4410;
		of main.eqgene w/o samechain, n=74894, nno=68374 no eqgene, 20%eqcds n=2451
	main.eqgene total n=79899, nno=69957,  20%eqcds n=4641, 
		
 >> main.eqgene w/ no samechain, 20%eqcds n=2451 should be marked same locus likely missing valid introns (but check some?)
 >> samechain w/ no main.eqgene, noeqgene n=3631: what effects? : many are t1 x t2+, not tested in main.eqgene
   * check qual of shared intronchains for same-loci: weak vs strong? some are utr overlap junk
   * probably should use only equalgene cds overlap to call same-locus, otherwise utr overlaps pile up.
   * redo equalgene for all main+alts.
   
egrep '(Funhe2Exx11m000052|Funhe2Exx11m040951)' kf2pub11_both.an2.ichains | cut -f1,2,4 
Funhe2Exx11m000052t1    14      id_t2:10,id_t3:10,id_t5:10,id_t6:10,id_t8:9,id_t9:9,id_t4:8,id_t7:8,
		> utrover? Funhe2Exx11m040951t1:3,Funhe2Exx11m040951t2:3,Funhe2Exx11m040951t3:3,Funhe2Exx11m040951t4:3,Funhe2Exx11m040951t5:3,Funhe2Exx11m040951t6:3
		> no main/alt.eqgene cds or utr overlap for g052 x g40951
Funhe2Exx11m040951t1    4       id_t2:4,id_t3:4,id_t4:4,id_t5:4,Funhe2Exx11m000052t1:3,id_t6:3
	> share UTR exons, no CDS, w/ g52t1: Scaffold10165:612509-612543,642313-642387

==> kf2pub11_both.an2.ichains <==
	nt=183032 tr have ichains, and ng=29136 evasm loci, of nt=322715,ng=92004 in kf2pub11 mRNA; ng1=27492 loci with t1 in ichain
Funhe2Exx11m000001t1    65      127     na      Funhe2Exx11m007299t12:1NO,Funhe2Exx11m035699t12:1NO,Funhe2Exx11m035699t7:1NO
Funhe2Exx11m000002t1    66      131     na      Funhe2Exx11m000007t1:3NO,Funhe2Exx11m000007t2:3NO,Funhe2Exx11m000007t3:3NO,Funhe2Exx11m000007t4:1NO
Funhe2Exx11m000003t1    19      35      na      na
Funhe2Exx11m000003t3    4       6       na      na
Funhe2Exx11m000004t1    131     255     id_t2:129,id_t3:126,id_t4:126,id_t5:124,id_t6:124,id_t8:124,id_t7:123,id_t9:122,id_t10:121,id_t11:117,id_t12:113,id_t13:106,id_t14:99,id_t15:96,id_t18:94,id_t16:91,id_t17:85,id_t25:83,id_t19:82,id_t20:82,id_t21:79,id_t22:77,id_t23:77,id_t27:76,id_t24:72,id_t26:71,id_t28:67       id_t29:12NO,id_t30:12NO,id_t31:8NO,Funhe2Exx11m028732t1:8NO,Funhe2Exx11m052498t1:5NO,Funhe2Exx11m067108t1:4NO,Funhe2Exx11m049230t1:1NO
Funhe2Exx11m000004t10   121     239     id_t1:121,id_t2:121,id_t3:121,id_t4:121,id_t5:121,id_t6:121,id_t8:121,id_t7:120,id_t9:119,id_t11:117,id_t12:111,id_t13:103,id_t14:96,id_t18:89,id_t15:88,id_t16:88,id_t19:79,id_t17:77,id_t23:77,id_t25:76,id_t22:75,id_t20:74,id_t21:74,id_t24:70,id_t26:68,id_t27:67,id_t28:64        id_t29:12NO,id_t30:12NO,id_t31:8NO,Funhe2Exx11m028732t1:8NO,Funhe2Exx11m052498t1:5NO,Funhe2Exx11m049230t1:1NO,Funhe2Exx11m067108t1:1NO


=cut
