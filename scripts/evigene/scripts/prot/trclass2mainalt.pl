#!/usr/bin/env perl
# trclass2mainalt.pl

=item trclass2mainalt

  $evigene/scripts/prot/trclass2mainalt.pl -idpre Anofunz4iEVm  -trclass  evg2anofunz4h.tgclass3 
  cut from evigene/scripts/evgmrna2tsa2.pl  
  
=cut

use FindBin;
use lib ("$FindBin::Bin","$FindBin::Bin/.."); # assume evigene/scripts/ layout: main, prot/, rnaseq/, ...

use constant VERSION => '2017.09.01'; # shiftAltEq cut fix; altnomainfix; '03.18'; # 2017 v4 update with asmrna_dupfilter4.pl

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname);
use cdna_evigenesub; #added # 

our $EVIGENES="$FindBin::Bin";  
our $dryrun=0; ## $DRYRUN ?
our $DEBUG= $ENV{debug}|| 0;

my $IDPREFIX= $ENV{idprefix} || 'EVGm'; 
my $SHOWDROPS=0; 
my $SIZESORT=1; 
my ($trclass,$output,$logfile);  
my ($trpath,$trname, $sradatah, %settings);
my ($pubid_format,$altid_format,$GBPROID);
my $pubidnum_start=0;
my $CULLXEQ=0;
my $preserveOldIds=0; # change to ok pattern, IDPREOK

my $optok= GetOptions(
  "class|trclass=s", \$trclass,
  "dropshow!", \$SHOWDROPS,  
  "sizesort!", \$SIZESORT,  # default:on 
  "idprefix=s", \$IDPREFIX,  # FIXME: idpre option  overwritten by spppref
  "dryrun|n!", \$dryrun, 
  #"keepoldids|preserveOldIds!", \$preserveOldIds,  
  "keepoldids|preserveOldIds=s", \$preserveOldIds,  
  "CULLXEQ!", \$CULLXEQ,  
  "debug!", \$DEBUG, 
  );

die "EvidentialGene trclass2mainalt -trclass evigenes.trclass [-idprefix $IDPREFIX ] 
  makes tables of public ids and main-alt linkage, from trclass table of pair-wise locus links of tr2aacds.
  Pipeline precedence: asmrna_dupfilter -aligntab evigenes.aligntab -outclass evigenes.trclass ; 
  version ", VERSION, "\n"
  unless($optok and $trclass);  

my $IDPREOK= ($preserveOldIds and $preserveOldIds=~/^[a-zA-Z]/) ? $preserveOldIds : $IDPREFIX;

MAIN();

sub MAIN
{
  my($upstatus,$upfiles,$uptemp,$upokids)= (0) x 9; 

  loggit(0, "BEGIN $0  input=",$trclass,"date=",`date`);
	# do_settings("restore",$trclass);  

	# ($cdnaseq,$trpath,$trname)= get_evgtrset($trclass,$cdnaseq,"publicset"); # subs.pm; cdnaseq => mrnaseq here
	#	loggit(0, "get_evgtrset=",$cdnaseq,$trpath,$trname); ## facount($cdnaseq)
	#	loggit(LOG_DIE, "Missing -mrna",$cdnaseq) unless($cdnaseq and -s $cdnaseq);

	($pubid_format,$altid_format)= make_IDPREFIX(); # default abbrev of $organism now

	my($maintab,$pubids,$nmaintr,$nalltr)= trclass2maintab($trclass,"publicset",$upokids);
	# loggit(0, "trclass2maintab primary n=",$nmaintr,"allntr=",$nalltr,$pubids); 

  # skip this?
  # evigene/scripts/rnaseq/asmrna_altreclass.pl -trclass $trclass -altrenum -out -debug
  # loggit(0, "DONE at date=",`date`);
}


=item cullExonEq

  # add for updated flags from genomap data, althi1.exonequal can be culled
  # here maybe cull short, 1xon alts also
   cullExonEq FIXME* altmapxev now == altmapxe99.100, altmapxeCC.XX, where CC = cds-align%, XX = exon-over%
   .. cull should measure both, skip where XX, or either? below ~90..97 levels to avoid cull valid alts
   eg. altmapxe85.96, 
   Aedesg12bEVm000004t2,t3 altmapxe98.96 * same nex=29 but different exons == valid alts
   supercont1.69	evg12aedes	exon	1193328	1200847	100	+	Parent=Aedesg12bEVm000004t2;ix=18   *e
   supercont1.69	evg12aedes	exon	1193328	1200851	99	+	Parent=Aedesg12bEVm000004t3;ix=18   +e
   supercont1.69	evg12aedes	exon	1201310	1203497	99	+	Parent=Aedesg12bEVm000004t3;ix=19   -b
   supercont1.69	evg12aedes	exon	1201318	1203497	100	+	Parent=Aedesg12bEVm000004t2;ix=19   *b
      ...
   supercont1.69	evg12aedes	exon	1208459	1208576	+	Parent=Aedesg12bEVm000004t2;ix=24
   supercont1.69	evg12aedes	exon	1208459	1208576	+	Parent=Aedesg12bEVm000004t3;ix=24
   supercont1.69	evg12aedes	exon	1217618	1217775	+	Parent=Aedesg12bEVm000004t3;ix=25 <<
   supercont1.69	evg12aedes	exon	1227691	1227866	+	Parent=Aedesg12bEVm000004t2;ix=25 <<
   supercont1.69	evg12aedes	exon	1228196	1228350	+	Parent=Aedesg12bEVm000004t2;ix=26
   supercont1.69	evg12aedes	exon	1228196	1228350	+	Parent=Aedesg12bEVm000004t3;ix=26
  
  Aedesg12bEVm000033t6/Aedesg2EVm000025t3 = Aedesg2EVm000025t2/altmapxe96.98
  supercont1.593  386400  386513  100     + Parent=Aedesg12bEVm000033t5;ix=18
  supercont1.593  386400  386513  100     + Parent=Aedesg12bEVm000033t6;ix=18
  supercont1.593  390634  390759  100     + Parent=Aedesg12bEVm000033t6;ix=19 << t6
  supercont1.593  391146  391217  100     + Parent=Aedesg12bEVm000033t5;ix=19
  supercont1.593  391146  391217  100     + Parent=Aedesg12bEVm000033t6;ix=20
  ..
  supercont1.593  407003  407168  100     +  Parent=Aedesg12bEVm000033t5;ix=30
  supercont1.593  407003  407168  100     +  Parent=Aedesg12bEVm000033t6;ix=31
  supercont1.593  407433  407600  100     +  Parent=Aedesg12bEVm000033t5;ix=31  <<t5
  supercont1.593  418239  418381  100     +  Parent=Aedesg12bEVm000033t5;ix=32
  supercont1.593  418239  418381  100     +  Parent=Aedesg12bEVm000033t6;ix=32
  
  Aedesg12bEVm000035t4/Aedesg1EVm000041t4 = Aedesg12bEVm000035t3/Aedesg1EVm000041t3/altmapxe87.100 
    >> diff here is due to low map cover; xev==100 is spurious.
    >> check also chrmap:87a,99i,6570l,2x,.. vs chrmap:88a,99i,6564l,2x,..
    
=cut

use constant minXEQ => 99; # exon-equal, high to ensure keep valid alts 
use constant minXCP => 95; # cds-align min for cullExonEq
use constant minXAL => 90; # chr-align min for cullExonEq

sub cullExonEq {
  my($md,$ads,$alt,$notes)=@_;
  my %culled; my $ncul=0; my $galn; 
  for my $td (@$ads) {
    my $cl= $alt->{$md}{$td};
    my $note= $notes->{$td};
    if($cl=~/^althi/ and not($note=~/refbest/) and $note=~/feq:([^;:\s]+)/ ) { 
      my $feq=$1;
      my @xe= grep /altmapxe/, split",",$feq;
      for my $xe (@xe) {
        my($xd,$xpev)= $xe=~m,(\w+)/altmapxe(.+),;  # Aedesg1EVm000007t3 x Aedesg1EVm000007t2/altmapxe98.96 
        my($xcp,$xev)= $xpev=~m/(\d+).(\d+)/;
        my $xcl= $alt->{$md}{$xd}||"";
        if($xev >= minXEQ and ($galn)= $note=~m/chrmap:(\d+)a/) { $xev=0 if($galn < minXAL); }
        if($xev >= minXEQ and $xcp >= minXCP and not($xcl =~ /^cull/)) { 
          $alt->{$md}{$td}="cull".$cl;
          $culled{$td}=$xe; $ncul++; last;
        }
      }
    }
  }
  return \%culled;
}

=item shiftAltEq

 .. shift altnomain to new main if feq: overlap statement says other main exists
 .. Arath5EVm000200t5 is drop/lower qual to t2
 .. Arath5EVm000200t3/altmap99.o85 is main, same locus .. shift t2 to altof t3
 Arath5EVm000200t2	maybeok	althi1	Arath5EVm000200t5	99/99/./altmap99	1509,90%,complete	
  aaref:4338,AT5G43900.2,refbest,chrmap:57a,100i,4530l,22x,chr5:17657241-17662149:-,pflag:72,
  feq:Arath5EVm000200t1/altpar0.no.o35,Arath5EVm000200t3/altmap99.o85,Arath5EVm000200t4/altmap98.o73

        my($othermain)= shiftAltEq($td,$aaqual{$td},$notes{$td},\%main);

=cut

use constant minShiftAltCP => 25; # cds-align min for shiftAltEq, low is ok, pick 1st=best of feq: list
use constant minShiftAltXP => 15; # x-align min for shiftAltEq, low is ok, pick 1st=best of feq: list

sub shiftAltEq {
  my($td,$cl,$aaqual,$note,$mainref,$backalt)=@_;
  my ($newmain,$galn)=(0); 
  if(1){
    # my $cl= $alt->{$md}{$td};
    # my $note= $notes->{$td};
    if($cl=~/^althi/  and $note=~/feq:([^;:\s]+)/ ) { # and not($note=~/refbest/)
      my $feq=$1;
      my @am= grep /altmap/, split",",$feq;
      for my $am (@am) {
        my($xd,$xpev)= $am=~m,(\w+)/altmap(.+),;  # Aedesg1EVm000007t3 x Aedesg1EVm000007t2/altmap98.o96 
        my($xcp,$xov)= $xpev=~m/(\d+)\D+(\d+)/;
        if($xcp >= minShiftAltCP and $xov > minShiftAltXP) { # xov is location cds-overlap%, require what min?
          my $mcl= $mainref->{$xd}; # $main{$md}= $cl;
          if($mcl) { $newmain=$xd; return ($newmain); }
          elsif(my $bmd= $backalt->{$xd}) { 
            # is this.td alt of other alt.xd that has main? this may return bad main
            # .. want to know overlap of this.td and bmd
            if($mainref->{$bmd}) {  $newmain=$bmd; return ($newmain); } 
          }
        }
      }
    }
  }
  return ($newmain);
}


sub trclass2maintab
{
  my($trclass,$pubdir, $okids)=@_;
  my $ntr=0;  my $nerr=0;
  my $mainindex= $pubidnum_start;
  ## okids = \%validids after merge filesets
  my $hasokids= ($okids and ref($okids) and scalar(%$okids))?1:0;
  
  my $maintab = makename($trclass,".mainalt.tab","trclass");  # > $pt.mainalt.tab
  my $pubidtab= makename($trclass,".pubids","trclass");   
  if(not -f $pubidtab and $pubdir and -d $pubdir) {
  	my($pubd,$ft);
  	($pubd,$ft)= getFileset($pubdir,'pubids',$pubd);  $pubidtab=$ft if($ft);  
  	($pubd,$ft)= getFileset($pubdir,'mainalt.tab',$pubd);  $maintab=$ft if($ft);  
  }
  return($maintab,$pubidtab,$mainindex,$ntr) if( -s $maintab and -s $pubidtab);# or dryrun ..

  my($ok,%main,%mainsize,%alt,%altsize,%maindrops,%altdrops,%didaltdrops,%balt,%drop,$outh,$outpubidh,$inh);
  my(%aaqual,%piad,%notes, %newmain, %altdefer);
  
  ($ok,$inh)= openRead($trclass);
  $ok= open($outh,'>',$maintab) if($ok);
  $ok= open($outpubidh,'>',$pubidtab) if($ok);
  unless($ok) { loggit(1,"ERR: parse $trclass TO $maintab"); return; }

  ## FIXME: only althi are reliably locus alternates; altmid .. are more likely paralogs
  while(<$inh>) {
    next unless(/^\w/); chomp;
    my($td,$ok,$cl,$md,$piad,$aq,$fl)=split "\t";
    
    ## new data bug, md == 0, piad == 0, from asmrna_dupfilter2b; keep? call these 'noclass' for now?
    if($DEBUG and not $md) { $md="$td.miss"; $nerr++; }
    else { $md ||= $td; }
 
    unless($cl and $md and $aq) { $nerr++; loggit(1,"ERR: trclass inline:$_"); next; } # what?? report?

		my $dropit=0;
		## maybeok needs handling, drop? or keep/mark/cull ??
		## should revise asm dupfliter to avoid these maybes, includes 'refbest', others some good/bad
		## now all are from exoneq flag, all 'althi1|part|frag' classes; 5702 refbest(dups?) of 13059
		## * maybe keep all maybes ; found some refbest/good being dropped here..
		## * but maybe is large subset, most althi1: 236768 drop, 195819 maybeok, 90660 okay in evg4anoalb.trclass 
    ## .. use hoscore if available, else keep til have hoscores?
		if($ok eq 'maybeok') { 
		   $ok="okay"; $cl.="maybe"; # if($fl=~/refbest/) { $ok="okay"; $cl.="maybe"; } 
		}
		
		# FIXME170503: drops may be false-main of others, need drop record td/md at DeferAltNoMain17 below
    if($ok ne 'okay') { $drop{$td}=$md; $cl=$ok.$cl; $dropit=1;  } # next unless($SHOWDROPS);  OPTION: include drops?
    elsif($hasokids and not $okids->{$td}) { $drop{$td}=$md;  $dropit=1; $cl='dropid'.$cl; } # next unless($SHOWDROPS);  unless $td =~ /utrorf/ ??
    #o if($ok ne 'okay') { $drop{$td}=1; $cl=$ok.$cl; $dropit=1;  } # next unless($SHOWDROPS);  OPTION: include drops?
    #o elsif($hasokids and not $okids->{$td}) { $drop{$td}=1;  $dropit=1; $cl='dropid'.$cl; } # next unless($SHOWDROPS);  unless $td =~ /utrorf/ ??

    ## 98/100/-sense/PitaEaR000975t24 << buggers.
    ## all piad entries should have pd and sense/asense fields, placeholder where needed '.' ?
    ## NEW SYNTAX tgclass for piad: 99/96/./altmap95/Anofunz4hEVm003223t1
    # need more sensible way to stick in genomap flags into trclass struct
    ## UPD 170901: misplaced pd/other high ident ids in notes:feq:xxx
    ##  Daplx7b3EVm000398t4, new main, should become md = pd = main2 for t1/t2 pair to break circular alt/alt link
    ## eg: Daplx7b3EVm000398t1	okay	althi	Daplx7b3EVm000398t2	99/95/./altmap95	1989..	
    ##     aaref:...,pflag:0,feq:Daplx7b3EVm000398t2/altmap95.o95,Daplx7b3EVm000398t4/altmap90.o90

    my $tgalt="";
    my($pi,$pa,$asense,$pd,@px)=split"/",$piad; # asense before pd ** NEED To revise asm dupfilter to clarify
    if($pd=~/^(alt|main|par|nocla)/) { $tgalt=$pd; $pd=(@px<1)?"":shift @px; } 
    if(@px) {
       # @px is what else?
    }
     
    if($asense =~ /sense/) { $pd="" unless($pd =~ /^\w/); }
    elsif($asense =~ /^\w/) { $pd=$asense; $asense=""; }
    if($pd=~/^self/) { $pd=""; }
    elsif($pd =~ /^\w/) { $md=$pd; } # V4: this is 'look-ahead' locus main id, do away with this usage
   	$cl=~s/a2$//;  #? dont need a2 AADUP qualifier?

 		$aaqual{$td}= $aq; # save for pubtable
 		$piad{$td}= $piad; # save for pubtable
 		$notes{$td}= $fl;
 		
 		my $aasize= ($aq =~ m/(\d+).(\d*)/) ? "$1.$2" : 0; # size.pCDS for best sort
 		## FIXME: have 2+ mains/$md from inexchain.tab/pubids merges ( main,reloc to locus w/ main)
 		## .. one needs to be reclassed as alt. ? rely on $td eq $md or not?

# V4:2017.03: FixAltNoMain16 not good enough, mis-ordered input where alts preceed their main
#    need to collect all of table before reassign alt => mainl
use constant FixAltNoMain16   => 1; # this may work; fixme2: part == alt here? or reclass input trclass?
use constant DeferAltNoMain17 => 1; # this may work; fixme2: part == alt here? or reclass input trclass?

    my $isalt=1;
    if($cl =~ /^(main|noclass)/) { 
      $isalt= ($main{$md})?1:0;
      ##x $isalt=0; if($td ne $md) { if($main{$md}) { $isalt=1; } }
      if($isalt) {
        $cl =~ s/^(main|noclass)/althim/;
      } else {
    	  $main{$td}=$cl; $balt{$td}=$td;
    	  $mainsize{$td}= $aasize;  
    	}
    } #o: else
    elsif($cl =~ /^alt/) {  # need this recip? to reset alt>main if best after reorder?
      $isalt= ($main{$md})?1:0;
      unless($isalt) {
if(DeferAltNoMain17) {
        my $mdfix= ($drop{$md} and $drop{$md} eq $td)? $td : $md; #upd1705
        # preserve this localt, fixup after input all of locus pairs
        # my @altrec=(); #($td,$ok,$cl,$md,$piad,$aq,$fl) == input row; aasize
        my @altrec=( $td, $cl, $mdfix, $aasize);
        push @{$altdefer{$mdfix}}, \@altrec;
} elsif(FixAltNoMain16) {  #no good unless out-of-order alts preceding mains have  md/mainid attached to piad value
          my $nmd= $newmain{$md} || $td;
          if($nmd eq $td) {
          $cl =~ s/^alt\w+/mainl/;
          $main{$nmd}= $cl; $balt{$nmd}=$td;
          $mainsize{$nmd}= $aasize;
          } else {
          $isalt=1;
          }
          $newmain{$md}=$nmd; $md=$nmd;
} else {
        $cl =~ s/^alt\w+/mainl/; # this is wrong also
    	  $main{$td}=$cl; $balt{$td}=$td;
    	  $mainsize{$td}= $aasize; # ($aq =~ m/(\d+).(\d*)/) ? "$1.$2" : 0;
}
      }
    }
    
    if($isalt) { 
    	$alt{$md}{$td}= $cl; $balt{$td}=$md; 
    	$altsize{$td}= $aasize; 
    }  
  }

if(DeferAltNoMain17) {

  # upd1705? fix for dropped false mains ?
  for my $md (grep { $drop{$_} } sort keys %altdefer) {
    my $newmd=0;
    my $otherid= $drop{$md};
    my @aset= sort{ $$b[3] <=> $$a[3] or $$a[0] cmp $$b[0] } @{$altdefer{$md}};
    for my $ar (@aset) {
       my( $td, $cl, $mdd, $aasize)= @$ar;
       if($td eq $otherid) { $newmd=$td; last; } # debug log this?
       }
    if($newmd) {
      my @mdalts= @{$altdefer{$md}};
      map{ $_->[2]= $newmd; } @mdalts;
      if(my $nalts= $altdefer{$newmd}) { unshift @mdalts, @$nalts; } # check dups?
      $altdefer{$newmd}= \@mdalts;
      my $nalt=@mdalts;
      loggit(0, "DBG fixaltmaindrop dropm:$md, newm:$newmd, nalt=$nalt") if $DEBUG;
    }
  }
  
  for my $md (sort keys %altdefer) {
    my @aset= sort{ $$b[3] <=> $$a[3] or $$a[0] cmp $$b[0] } @{$altdefer{$md}};
    for my $ar (@aset) {
      my( $td, $cl, $mdd, $aasize)= @$ar;
      my $isalt= ($main{$md})?1:0;
      my $thismd= $md;
      
      unless($isalt) {
        my($othermain)= shiftAltEq($td,$cl,$aaqual{$td},$notes{$td},\%main,\%balt);
        if($othermain) { $thismd=$othermain; $isalt=1; }
      }  
      
      unless($isalt) { ## FixAltNoMain16 moves here, out of input loop
        my $nmd= $newmain{$md} || $td;
        if($nmd eq $td) {
          $cl =~ s/^alt\w+/mainl/;
          $main{$nmd}= $cl; $balt{$nmd}=$td;
          $mainsize{$nmd}= $aasize;
        } else { $isalt=1; }
        $newmain{$md}=$nmd; $thismd=$nmd;
      }
    if($isalt) { 
    	$alt{$thismd}{$td}= $cl; $balt{$td}=$thismd; 
    	$altsize{$td}= $aasize; 
      }  
    }
  }  
}
  
  
  my %hasmain;
  my @amain= grep { not $main{$_} } sort keys %alt; # dropmain here now only for SHOWDROPS !
  foreach my $am (@amain) { 
    my $md= $balt{$am} || $am; ## $md=$am if($drop{$md});
  	if(!$main{$md} and $drop{$md}) { my $md1= $balt{$md}||""; if($md1 and $main{$md1}) { $md=$md1; } }
  	
    if($main{$md}) { my @at=keys %{$alt{$am}}; map{ $alt{$md}{$_}=$alt{$am}{$_}; $balt{$_}=$md; } @at; } 
    elsif($md) { 
     $main{$am}="NOMAIN";  # FIXME: get rid of these by finding alt main
     } 
  }
  
  foreach my $td (keys %balt) {
    my $md= $balt{$td} || $td; 
    $main{$md}="NOMAIN" unless($main{$md});# FIXME: get rid of these by finding alt main
  }


     
  ## add headers to these:
  #originalID     MainClass  Alternates
  #Public_mRNA_ID         originalID      PublicGeneID    AltNum
  ##FIXME: use extended realt format now?
  # #Public_mRNA_ID originalID      PublicGeneID    AltNum  Class   AAqual  pIdAln  Notes

  print $outh '#'.join("\t",qw(originalID MainClass Alternates))."\n";
  #oprint $outpubidh '#'.join("\t",qw(Public_mRNA_ID originalID PublicGeneID AltNum))."\n"
  print $outpubidh '#'.join("\t",qw(Public_mRNA_ID originalID PublicGeneID AltNum Class AAqual pIdAln Notes))."\n"
    if($outpubidh);

  my %doneid=();
  my @mainlist;
  if($SIZESORT) { # or NOT IDSORT ?
  @mainlist= sort{ $mainsize{$b} <=> $mainsize{$a} or $a cmp $b } keys %main;
  } else { # IDsort, only special cases, like merge old/update, keeping old order for most
  @mainlist= sort{  $a cmp $b } keys %main;
  }

#upd1708: preserveOldIds option ..  keep input gene ids as best feasible, including altnum?
#  need all alts per main, check evg gene ids for 1 x many, also need check across loci for 2+ new loci w/ old geneids
#  where cant preserve geneid for all of locus, do what? new geneid, or modify old?
#  UPD: want/need main = t1, cant preserve old alt nums if alt classing changes ..
  my ($nnewids,$newidh)= ($preserveOldIds) ? preserveOldIds(\@mainlist, \%alt, \%drop, \%altsize) : (0, undef);

  foreach my $md (@mainlist) { 

    my @ad= sort{$alt{$md}{$a} cmp $alt{$md}{$b}
      or $altsize{$b} <=> $altsize{$a} or $a cmp $b } keys %{$alt{$md}}; 

    my $culls= ($CULLXEQ) ? cullExonEq($md,\@ad,\%alt,\%notes) : {}; #?? here
    my $ad= join",",map{ "$_/".$alt{$md}{$_} } @ad; 
    my $mc= $main{$md}; 
    
    # FIXME: change $mc eq "noclass" to "main" if alts exist, and vversa
    if($mc =~ /^NOMAIN/) { 
      # do below
    } elsif(@ad>0 and $mc !~ /^main/) {
      $mc =~ s/^\w+/main/;
      $main{$md}= $mc; 
    } elsif( @ad==0 and $mc !~ /^noclass/) {
      $mc =~ s/^\w+/noclass/;
      $main{$md}= $mc; 
    }
    
    if($SHOWDROPS) {
    	my @add= sort{$altdrops{$md}{$a} cmp $altdrops{$md}{$b}} keys %{$altdrops{$md}};  
    	map{ $didaltdrops{$_}=1 } @add; # later dump not didaltdrops
    	my $add=join",",map{ "$_/".$altdrops{$md}{$_} } @add; 
    	$ad .= ",$add" if ($add);
    }
    
    print $outh join("\t",$md,$mc,$ad)."\n";  # mainalt.tab
    
    if($outpubidh) { # should be required ??
      my $ialt= 0; my $needmain=0;
      my($cla,$aaq,$pida,$nots);
      $cla= $mc;  # cla=$main{$td}=$cl;  ; changed above?
      $aaq= $aaqual{$md}||"noaa";
      $pida=$piad{$md}||0;  
      $nots=$notes{$md}||"nonote";  
      if($mc eq "NOMAIN") { $cla=(@ad>0)?"main":"noclass"; } ## needs to change, to main? to noclass?
      elsif($mc =~ /^alt/) { } # is this were nomain show up? or @ad?

      my @sad= sort{ $altsize{$b} <=> $altsize{$a} or $a cmp $b } @ad;      
      
      if($drop{$md} or $doneid{$md}) { $needmain=1; }
      else {
      	$mainindex++; $needmain=0; # BUG: move below drop{}
      	my ($pubmrnaid,$pubgeneid);      	     	   
      	# if($nnewids and (my $nd= $newidh->{$md})) { 
      	#  $pubmrnaid=$nd; ($pubgeneid= $pubmrnaid) =~ s/t\d+$//; ++$ialt; }
      	# unless($pubmrnaid) {
      	($pubmrnaid,$pubgeneid)= make_pubid($md, $mainindex, ++$ialt, $newidh);
      	# }
      	print $outpubidh join("\t",$pubmrnaid,$md,$pubgeneid,$ialt,$cla,$aaq,$pida,$nots)."\n";  #n
      	$ntr++; $doneid{$md}++;
      	}
      
      foreach my $ad (@sad) {
        unless($drop{$ad} or $doneid{$ad}) {
        $cla=  $alt{$md}{$ad}||"nocl"; 
        $aaq=  $aaqual{$ad}||"noaa";
        $pida= $piad{$ad}||0; # fixme == piad above
        $nots= $notes{$ad}||"nonote"; # fixme
        # $cull= $culls->{$ad}||""; # $cla == "cullalt..";
      	if($needmain) { 
      	  $mainindex++; $needmain=0; 
      	  if($cla=~/^alt/){ $cla=(@sad>1)?"main":"noclass"; } 
      	}  
        my ($altmrnaid,$altgeneid);
      	# if($nnewids and (my $nd= $newidh->{$ad})) { 
      	#  $altmrnaid=$nd; ($altgeneid= $altmrnaid) =~ s/t\d+$//; ++$ialt; }
        # unless($altmrnaid) {
        ($altmrnaid,$altgeneid)= make_pubid($ad, $mainindex, ++$ialt, $newidh);
        # }
        print $outpubidh join("\t",$altmrnaid,$ad,$altgeneid,$ialt,$cla,$aaq,$pida,$nots)."\n"; 
        $ntr++; $doneid{$ad}++;
        }
      }
    }
  }
  
  if($SHOWDROPS) { # UPD.160911
    my @mains= sort keys %altdrops;
    for my $md (@mains) {
      my @misdrop= grep{ not $didaltdrops{$_} } keys %{$altdrops{$md}};
      if(@misdrop) {
        @misdrop= sort{$altdrops{$md}{$a} cmp $altdrops{$md}{$b}} @misdrop;  
        map{ $didaltdrops{$_}=1 } @misdrop; # later dump not didaltdrops
        my $misdrop= join",",map{ "$_/".$altdrops{$md}{$_} } @misdrop; 
        my $mc= $main{$md}||"NOMAINd"; 
        print $outh join("\t",$md,$mc,$misdrop)."\n";  # mainalt.tab
      }
    }
  }
    
  close($inh); close($outh);

	# push @publicset, $maintab, $pubidtab; #?
  return($maintab,$pubidtab,$mainindex,$ntr);  # return main,alt id hashes ....
}

#upd1708: preserveOldIds option ..  keep input gene ids as best feasible, including altnum?
#  need all alts per main, check evg gene ids for 1 x many, also need check across loci for 2+ new loci w/ old geneids
#  where cant preserve geneid for all of locus, do what? new geneid, or modify old?

sub idparts {
  my($id)=@_;
## FIXME for _G2 _C2? tags, G2,n is diff locus id
  my($gpre,$gnum,$gd,$ti,$gdup)=(0) x 9;
  $gd=$id; 
  # ($gdup)= ($gd=~s/_([GC]\d+)$//)?$1:0;
  ($gdup)= ($gd=~s/_(G\d+)$//)?$1:0;
  ($ti)= ($gd=~s/t(\d+)$//)?$1:0;
  ($gnum)= ($gd=~m/(\d+)$/) ? $1:0;
  ($gpre=$gd) =~ s/$gnum$//;
  if($gdup) { $gd.=$gdup; $gnum=0; } # gnum invalid?
  #($gd,$ti)= ($id=~m/^(\w+)t(\d+)$/) ? ($1,$2) : ($id,0);
  #($gpre,$gnum)= ($gd =~ m/^(\w+[a-zA-Z])(0\d\d+)$/)? ($1,$2) : ($gd,0);
  return($gpre,$gnum,$gd,$ti);
}

sub preserveOldIds {
  my($mainlist, $altsOfMain, $drop, $altsize)= @_;
  my(%gids, %gnums, %gdone, %newids, $gprefix);
  # return \%newids
  my $nids=0;
  
  ## prefix, gnum problem: have mixed gprefices .. diff gnum set for each
  my $GNEXTNUM=0;
  foreach my $md (@$mainlist) {
    my @okd = grep{ not($drop->{$_}) } ($md, keys %{$altsOfMain->{$md}} );
    for my $id (@okd) {
      my($gpre,$gnum,$gd,$ti)= idparts($id);
      next unless($gnum);
      $gnums{$gnum}++; # $gids{$gd}{$ti}= $id;
      $gprefix= $gpre unless($gprefix);
      }
  }
  my($glast)= sort{ $b <=> $a } keys %gnums;
  $GNEXTNUM= 9 + $glast;
  
  #my $idformat= $pubid_format; $idformat=~s/$IDPREFIX/N$gprefix/; # debug N addition
  my $idformat= "n$gprefix".'%06d';
  my(%havepubg);

  foreach my $md (@$mainlist) {
  
    # my @ad= sort{ $altsOfMain->{$md}{$a} cmp $altsOfMain->{$md}{$b} # class sort
    #  or $altsize->{$b} <=> $altsize->{$a} or $a cmp $b } keys %{$altsOfMain->{$md}}; 
    my @ad= sort{ $a cmp $b } keys %{$altsOfMain->{$md}}; 
    my @okd = grep{ not($drop->{$_}) } ($md,@ad);
  
    %gnums=(); %gids=();
    my $gnumfirst=0; my $gprefat=0;
    for my $id (@okd) {
      my($gpre,$gnum,$gd,$ti)= idparts($id);
      next unless($gnum); # for _G2,n.. new loci
      $gnums{$gnum}++; $gids{$gnum}{$ti}= $id;
      $gnumfirst=$gnum unless($gnumfirst or $gdone{$gnum}); 
      $gprefat= $gpre unless($gprefat); #?? problems? yes
      # if($gpre ne $gprefix) { } # error ?? #$gprefix= $gpre unless($gprefix);
    }
    
    my @gnums= sort{ $gnums{$b}<=>$gnums{$a}  or $a <=> $b } keys %gnums;
    unless($gnumfirst) { # pick most or first in (main) <<
      for(my $i=0; $i<=$#gnums; $i++) {
        unless($gdone{$gnums[$i]}) { $gnumfirst=$gnums[$i]; last; }
      }
    }
    # usage assumption: -keepids=IdPreA -idpre IdPreA, so inval gprefat will go to ++GNEXTNUM 
    unless($gprefat and $gprefat =~ /$IDPREOK/){ $gprefat=$IDPREFIX;  $gnumfirst=0; }
    #^^ bug, dupids
    $idformat= "n$gprefat".'%06d'; # maybe new for each main id??
    my($pubgn);
    do { 
      unless($gnumfirst) { $gnumfirst= ++$GNEXTNUM; } #??
      $pubgn = sprintf( $idformat, $gnumfirst);
      if($havepubg{$pubgn}) { $gnumfirst=0; }
    } while ($havepubg{$pubgn});
    $havepubg{$pubgn}=1; 

#try1 bad .. first locus, all mixed up
# Daplx7mEVm000001t1	Daplx7b3EVm000002t1	Daplx7mEVm000001  1
# Daplx7mEVm000001t2	Daplx7b3EVm000002t2	Daplx7mEVm000001  2
# NDaplx7b3EVm000002t000000t2	Daplx7b3EVm000076t2	NDaplx7b3EVm000002t000000 3   

if(1) { # try2
    my %tidone=(); my $timax= @okd; #? not same using orig ti > @okd
    for my $id (@okd) {
      my($gpre,$gnum,$gd,$ti)= idparts($id);
      my $tinew=0;
      if($gnum eq $gnumfirst and not $tidone{$ti}) { $tinew=$ti; }
      if($tinew==0) { do { $tinew++; } while( $tidone{$tinew} ); } # can put lowqual alts at ti top
      $tidone{$tinew}++; 
      $gdone{$gnumfirst}++;
      #above# my $pubgn = sprintf( $idformat, $gnumfirst); 
      my $pubti = sprintf( $altid_format, $tinew); # same ti or new? cant use same ti w/o checks
      my $pubid=  $pubgn . $pubti;
      $newids{$id}= $pubid;
      $nids++;
    }
}
if(0) { # try1
    for my $gn (@gnums) {
      my @ti=sort{ $a <=> $b } keys %{$gids{$gn}};
      my $gsame= ($gn eq $gnumfirst);
      
      unless($gsame) { $gnumfirst= ++$GNEXTNUM; } # NO,
      
      for my $i (@ti) {
        my $iid= $gids{$gn}{$i};
        $gdone{$gnumfirst}++;
        my $pubgn = sprintf( $idformat, $gnumfirst); 
        
        my $pubti = sprintf( $altid_format, $i); # same ti or new? cant use same ti w/o checks
        my $pubid=  $pubgn . $pubti;
        
        if($gsame and $pubid ne $iid) { } # mistake ? flag it?
        $newids{$iid}= $pubid;
        $nids++;
      }
    }
}
    
  } # mainlist
  
  return($nids,\%newids);
} # sub preserveOldIds


sub make_pubid # add preseveOld %newids here? use altnum w/ it or preserve old altnum?
{
  my($oid, $pubidnum, $altnum, $preservedIds)= @_;
  $pubidnum_start= $pubidnum; #?

  my($pubid,$pubgene);
  if($preservedIds and ref($preservedIds)) {
    if(my $pid= $preservedIds->{$oid}) {
      if(my($pgene,$palt)= $pid=~m/(\w+)t(\d+)$/) {
        $pubgene=$pgene;
        # if($palt ne $altnum) ..
        $pubid   = $pubgene . sprintf( $altid_format, $altnum); # or palt ?
        return($pubid,$pubgene);
      }
    }
  }
  $pubgene = sprintf( $pubid_format, $pubidnum); 
  $pubid   = $pubgene . sprintf( $altid_format, $altnum);
  return($pubid,$pubgene);
}

sub make_IDPREFIX
{
  my($existingID)= @_;
  my $digits=6; # or 7?
  my $ChangeDefaultIdPrefix= 0; # ($IDPREFIX eq $DEFAULT_SETTINGS{'IDPREFIX'}) ? 1:0;
  
  if($existingID and $existingID !~ m/^$IDPREFIX/) {
    $existingID=~ s/t\d+$//;
    my($prefix,$nums)= $existingID =~ m/^(\w+)(\d\d\d+)$/; ## prefix may have numbers, how many? 
    if($prefix) { $IDPREFIX= $prefix; $digits= length($nums) || 6; }
    $ChangeDefaultIdPrefix=0;
  }

  my $nd= ( $IDPREFIX =~ s/(0\d+)$// ) ? length($1) : $digits;
  my $pubid_format = $IDPREFIX.'%0'.$nd.'d';  
  my $altid_format = 't%d'; 
  # my $GBPROID= $IDPREFIX."_".$DATE;  
  return($pubid_format,$altid_format); # ,$GBPROID);
}
