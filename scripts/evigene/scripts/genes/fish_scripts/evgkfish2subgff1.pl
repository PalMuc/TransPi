#!/usr/bin/env perl
# evgkfish2subgff.pl
# Genome annots with cds2aa checked, trasm pubaa translate exceptions

use strict;

my $debug=$ENV{debug}||0;
my $TESTSET=$ENV{test}||0;  # subset, drop most alts..
my $NOSPLITSET=$ENV{nosplit}||0;  # subset, drop split genes
my $VERS=$ENV{vers}||"fc5"; # fc4..

my $EVIGENE=$ENV{evigene}||"/bio/bio-grid/mb/evigene";

my $S2='kf2a'; my $SANN='kf2rae5|kf2rae5alt'; # == kf2a, kfish2rae5h_kfish2a_findcds3t1.gff
my $S3='kf3n'; my $SOUT='kf2rae5h'; # == kf3n, kfish2rae5h_kfish2n_findcds3.gff
my $SPUB='pubaa'; # pub.aa replaces cds2aa

## redo SANN == S2SRC; $SOUT > $S3SRC;  $SRCOUT="kf2rae5$VERS";
## kf2rae5h

my $TRAMSID="Funhe2Exx11m";
my $DROPOID='Funhe2E6bm|Funhe2Emap3m|Funhe2Eq7m|Fungr1EG3m|fungrvel|rfung|^xp|UniRef|Ictpun';  #? estnewbg
## drop most.. say instead what to keep: 
my $KEEPOID='\bFunhe2Exx11m|\bkf4bAUG|\bAUG|PASA'; #? ^Funhe5EG| .. other?

# kf4bAUGpia9cs55g ??
my $SplitGeneFixID= 1; # ncbi sub needs; make split gene parts have new ids w/ suffix _C1/2 ..
my $IDSplitSuffix =  "_C";  

my $AADIFFMIN=$ENV{aadiff}||75; 
    #^^ makebestalt call, pay attention to AADIFF qual, use pubaa where cds2aa is short/poor match
    # 44000/110164 are >= 80%; 60000/110000 are >= 50%
    # note this aadiff % may be flaky calc. subsitute cds2aa.sizeNoX / pubaa.sizeNoX ?
    
# my $DROPAN="Target|match|score"; 
my $DROPAN='Target|match|score|oname|groupname|ocds|inshort|utrx|xtrim'; 

## oid=Funhex11,Funhexother << remove 2nd troid
# exon annots: 13651 error from gmap; keep/report in gb.annot? Region=map-exception...
# Parent=Funhe2EKm000208t1;error=ERROR.span:genome_span:1070,tr_span:1241,KN805525.1:5989263-5988193171
# Parent=Funhe2EKm000500t1;error=ERROR.span:genome_span:3687,tr_span:3875,KN805527.1:3629399-3625712188;

# also drop,maybe: ocds=NEn,NE0,NE2,NE3,NE4,NE5;inshort=xcut3:-32;
# also maybe drop extras: oname=othername,  
#------ data ----------------------------------------------

my $ORIG_bestidtab="gmapn/kfish2rae5h_asmbest_aahogmap.tab2";
my $bestidtab="gmapn/kfish2rae5h_asmbest_aahogmap.tab3";
#^vv^ make 2nd table for others, alts, etc. best of kf3n, kf2a gmap and/or cds2aa size
my $bestalttab="gmapn/kfish2rae5h_altbest.tab";
  #v1: 26387 kf2a, 83776 kf3n; note this includes all trs, asmbest_aahogmap as well
  #v2:  6511 kf2a, 40476 kf3n, 63176 pubaa : use %aadiff to call pubaa .. can do any better w/o hoval?
  #v3: 20350 kf2a, 58198 kf3n, 31615 pubaa : use only aasize-cds/aasize-pub, better than aadiff
## FIXME: revise makealtbest/ altbest.tab to combine asmbest_aahogmap.tab
## .. and correct asmbest cases where t1 mismatch to talt and poor map, pubaa/kf3n should be pubaa/kf2a

# add new table of id.src errors: Scaf change for kf2a, split gene missing exons, other ..?
my $ErrorIdTab= "gmapn/kfish2rae5h_kfish2a_findcds3.cantuse.ids";

# * FIXME: NEW Annots from pubtsa6 : Name, .. supercede pubgenes.gff annots
# ?? from pubtsa6/kfish2rae5x11t6_tsasubmit.aa.gz or .mrna.gz ?
#  also need pubaa prot, ..
my $pubaain="pubtsa6/kfish2rae5x11t6_tsasubmit.aa.gz";

my $scafidtab="pubgenome/ncbifunhe3_kfish2asm.sametab";
my $SUBMITCONF="pubgenome/evigene_kfish2_gbsubmit.conf";

# my $onam="kfish2rae5h_kfish2n_";
my $otmpgff="gmapnsub1/kfish2rae5h_kfish2n_${VERS}.gff";
my $outsrtgff="gmapnsub1/kfish2rae5h_kfish2n_${VERS}ans.gff";
my $outidtab="gmapnsub1/kfish2rae5h_kfish2n_${VERS}.idtab";  # all input mRNA ID.src + action for debug

my @ingff= qw(
  gmapn/kfish2rae5h_kfish2a_findcds3t1.gff.gz 
  gmapn/kfish2rae5h_kfish2a_findcds3ta.gff.gz 
  gmapn/kfish2rae5h_kfish2n_findcds3.gff.gz
  );

my($didhead);
my(%scaf,%scafn); # dont need scafn ?
my(%pubaa,%puban,%an,%didput,%seenid); # pub ann from 1 gff input, use for next..
my(%mapq,%cdsq,%mapv,%cdsv);  # for makeAltBestTab
my %skipids=(); # ErrorIdTab

#------ methods ----------------------------------------------
sub MAINstub {}

if($ENV{makealt}) { 
  #?? my ($bestidh) = readAsmBestTab($bestidtab); # %bestid{pid}{id,src,hoval,gmap,..}
  makeAltBestTab($bestalttab,undef); #? use skipids below?
  exit; 
}

# add new table of id.src errors: Scaf change for kf2a, split gene missing exons, other ..?
# and/or adjust bestid/bestalt tabs for these errors ?

#above: $ErrorIdTab= "gmapn/kfish2rae5h_kfish2a_findcds3.cantuse.ids";
if( open(F,$ErrorIdTab) ) {
  my $skipsrc=$S2; # if($ErrorIdTab =~ /kfish2a/);
  warn "in.skipids: src=$skipsrc, $ErrorIdTab\n" if $debug;
  while(<F>){ my($id,$err)=split; $skipids{$id}{$skipsrc}= $err; } close(F);
}

my ($bestidh) = readAsmBestTab($bestidtab); # %bestid{pid}{id,src,hoval,gmap,..}
my ($bestalth)= readAsmBestTab($bestalttab, $bestidh); 
# %bestid= %$bestidh; # global, for readGff()

readScaffoldTab($scafidtab);  

readPubaa($pubaain); # # return(\%pubaa,\%puban);

my $outh=undef;
if($debug>1 or $otmpgff =~ /stdout/) { 
  $outh= *STDOUT;
} else {
  warn "#gff.out: $otmpgff\n" if $debug;
  rename($otmpgff,"$otmpgff.old") if(-s $otmpgff);
  open($outh,'>',$otmpgff) or die $otmpgff;
}

my($nint,$nputt,$ndropt)=(0) x 3;
for my $ingff (@ingff) {
  my $flags=""; # per ingff ?
  my($nin,$nput,$ndrop)= readGff($ingff,$outh,$flags); #  does putSubGff
  $nint+=$nin; $nputt+=$nput; $ndropt+=$ndrop;
}
warn "#gff.totout: in=$nint,put=$nputt,drop=$ndropt\n" if($debug);
close($outh) unless($otmpgff =~ /stdout/);

warn "#idactions: $outidtab\n" if($debug);
if(open(my $tabh,'>',$outidtab)) {
  # %didput{$id} == nput ,%seenid{$id}{$insrc} == $action ..
  print $tabh "#gff.totout: in=$nint,put=$nputt,drop=$ndropt\n";
  my @insrc=($S2,$S3);
  print $tabh join("\t","AQueryID","nPut","Src2","S2act","Src3","S3act")."\n";
  for my $id (sort keys %seenid) {
    my @s= sort keys %{$seenid{$id}}; @s= @insrc if(@s<2);
    my $np= $didput{$id}||0;
    print $tabh "$id\t$np";
    for my $s (@s){ my $sv= $seenid{$id}{$s} || "noval"; print $tabh "\t$s\t$sv"; } 
    print $tabh "\n";
    }
  close($tabh);
}

sortGff($otmpgff, $outsrtgff) if(-s $otmpgff);

#------------------------------------------------------

=item try1

  valid tr tables:
   35053 gmapn/kfish2rae5h_asmbest_aahogmap.tab2
  110164 gmapn/kfish2rae5h_altbest.tab
  101356  .. with map loc : should have mRNA for each of these

  try3#gff.totout: in=236338,put=101733,drop=8926 .. is this correct? put includes Splits

  grep mRNA gmapnsub1/kfish2rae5h_kfish2n_fc5.gff | wc -l
     23958  ** too few, what of alts?
  
  grep mRNA gmapnsub1/kfish2rae5h_kfish2n_fc5.gff | 
    grep -c insrc=kf3n =   114  ** way to few here.. table reading bug? * DROPped *
        asmbestaa n=16863; altbest n=83821
    grep -c insrc=kf2a = 23844
        asmbestaa n=1997; altbest n=26342

=item try2

 env  makealt=0 debug=1 ./evgkfish2subgff.pl       
  #asmbest.in: gmapn/kfish2rae5h_asmbest_aahogmap.tab2
  #asmbest.in: gmapn/kfish2rae5h_altbest.tab
  #scaftab.in: pubgenome/ncbifunhe3_kfish2asm.sametab
  #pubaa.in: pubtsa6/kfish2rae5x11t6_tsasubmit.aa.gz
  #gff.out: gmapnsub1/kfish2rae5h_kfish2n_fc5.gff
  #gff.in: gmapn/kfish2rae5h_kfish2a_findcds3t1.gff.gz
  #gff.out: in=37063,put=2899,drop=0
  #gff.in: gmapn/kfish2rae5h_kfish2a_findcds3ta.gff.gz
  #gff.out: in=78520,put=20945,drop=0
  #gff.in: gmapn/kfish2rae5h_kfish2n_findcds3.gff.gz
  #gff.out: in=120755,put=114,drop=86701   << ** DROPS why?
  #gff.totout: in=236338,put=23958,drop=86701  
  
=item try3
  env  makealt=0 debug=1 ./evgkfish2subgff.pl 
  #asmbest.in: gmapn/kfish2rae5h_asmbest_aahogmap.tab2
  #asmbest.in: gmapn/kfish2rae5h_altbest.tab
  #scaftab.in: pubgenome/ncbifunhe3_kfish2asm.sametab
  #pubaa.in: pubtsa6/kfish2rae5x11t6_tsasubmit.aa.gz
  #gff.out: gmapnsub1/kfish2rae5h_kfish2n_fc5.gff
  #gff.in: gmapn/kfish2rae5h_kfish2a_findcds3t1.gff.gz
  #gff.out: in=37063,put=2899,drop=0
  #gff.in: gmapn/kfish2rae5h_kfish2a_findcds3ta.gff.gz
  #gff.out: in=78520,put=20945,drop=0
  #gff.in: gmapn/kfish2rae5h_kfish2n_findcds3.gff.gz
  #gff.out: in=120755,put=77889,drop=8926       << * FIXED drop bug
  #gff.totout: in=236338,put=101733,drop=8926
  
  * Splits problem still , no Split= mRNA..
  cat gmapnsub1/kfish2rae5h_kfish2n_fc5.gff | grep mRNA | grep bestaa=pubaa | wc -l
   11385 ?? is this right count, includes Splits .. should be 13766 ?

  cat kfish2rae5h_asmbest_aahogmap.tab2 | grep -v '#' |  grep -v NOPATH | egrep -v 'kf2a|kf3n' | wc -l
   13766; but this counts 'nohoref', need horef for pubaa ? n=9714

=item try4

  env  makealt=0 debug=1 ./evgkfish2subgff.pl 
  #asmbest.in: gmapn/kfish2rae5h_asmbest_aahogmap.tab2
  #asmbest.in: gmapn/kfish2rae5h_altbest.tab
  #scaftab.in: pubgenome/ncbifunhe3_kfish2asm.sametab
  #pubaa.in: pubtsa6/kfish2rae5x11t6_tsasubmit.aa.gz
  #gff.out: gmapnsub1/kfish2rae5h_kfish2n_fc5.gff
  #gff.in: gmapn/kfish2rae5h_kfish2a_findcds3t1.gff.gz
  #gff.out: in=37063,put=2899,drop=0
  #gff.in: gmapn/kfish2rae5h_kfish2a_findcds3ta.gff.gz
  #gff.out: in=78520,put=20945,drop=0
  #gff.in: gmapn/kfish2rae5h_kfish2n_findcds3.gff.gz
  #gff.out: in=120755,put=84063,drop=8926
  #gff.totout: in=236338,put=107907,drop=8926 : a few more, Splits?
 
  bestaa=pubaa : 13778 minus Split=2 .. closer to above aahogmap count
  
  grep  mRNA gmapnsub1/kfish2rae5h_kfish2n_fc5.gff | grep Split= | wc -l
   11132;  kf2a=4958; kf3n=6174, pubaa for 4786
  grep  mRNA gmapnsub1/kfish2rae5h_kfish2n_fc5.gff.old | grep Split= | wc -l
    0

=item try5 : redo makealt, more pubaa for alt w/ small aa
  .. problem new bestalt, pubaa calls not used? same bestaa=pubaa as before ..
  
  env  makealt=0 debug=1 ./evgkfish2subgff.pl 
  #asmbest.in: gmapn/kfish2rae5h_asmbest_aahogmap.tab2
  #asmbest.in: gmapn/kfish2rae5h_altbest.tab
  #scaftab.in: pubgenome/ncbifunhe3_kfish2asm.sametab
  #pubaa.in: pubtsa6/kfish2rae5x11t6_tsasubmit.aa.gz
  #gff.out: gmapnsub1/kfish2rae5h_kfish2n_fc5.gff
  #gff.in: gmapn/kfish2rae5h_kfish2a_findcds3t1.gff.gz
  #gff.out: in=37063,put=2512,drop=0
  #gff.in: gmapn/kfish2rae5h_kfish2a_findcds3ta.gff.gz
  #gff.out: in=78520,put=16701,drop=0
  #gff.in: gmapn/kfish2rae5h_kfish2n_findcds3.gff.gz
  #gff.out: in=120755,put=71474,drop=2337
  #gff.totout: in=236338,put=90687,drop=2337  <<? why fewer put/drop? splits? no, bug .. try6

  $evigene/scripts/bestgenes_update.pl -act=sort -debug -cadd addgene=1 -vers kf2rae5h -conf pubgenome/evigene_kfish2_gbsubmit.conf -in gmapnsub1/kfish2rae5h_kfish2n_fc5.gff -out gmapnsub1/kfish2rae5h_kfish2n_fc5ans.gff
  # gene count for kf2rae5h: ngene=34340, nmrna=89669, nexon=1763741, nother=0

=item try6 : after redo makealt, fix for alts

  .. as try5 .. 
  #gff.in: gmapn/kfish2rae5h_kfish2n_findcds3.gff.gz
  #gff.out: in=120755,put=93140,drop=8965   # corrected
  #gff.totout: in=236338,put=112353,drop=8965

  pubaa in kfish2rae5h_kfish2n_fc5.pep.report 
    41524 mismatches
    52299 matches
  
  mRNA n=112353; Split=1/2 n=8819/8822, 17641 tot
  bestaa=pubaa n=45835 minus Split=2
  -- split problems for tbl2asn, need diff IDs for parts? at least for prots
  -- also sort splits independently, sort by scaffolds at top level
  -- split IDs, NCBI uses 'A','B' example but undefined; use  _C1/2 ? S1/S2 ?
    .. same split idsuffix for locus_tag, transcript_id, protein_id
    
  NCBI ft Split gene annot: 
  * annotate these as separate genes with unique locus_tags, plus separate CDS/mRNAs 
    with different protein_id's and transcript_id's. * In addition, link the features together with 
  * notes that refer to the other part of the gene.*
  ..  do not create extremely short features 
    split end cutoff 50bp? diff for cds/utr?, ie utr split end: 200bp or 100bp + valid intron? 
      cds split end: 120bp ? 60 bp?

=item try 7
  
  env  test=1 vers=fc6t debug=1 ./evgkfish2subgff.pl 
  #asmbest.in: gmapn/kfish2rae5h_asmbest_aahogmap.tab2
  #asmbest.in: gmapn/kfish2rae5h_altbest.tab
  #scaftab.in: pubgenome/ncbifunhe3_kfish2asm.sametab
  #pubaa.in: pubtsa6/kfish2rae5x11t6_tsasubmit.aa.gz
  #gff.out: gmapnsub1/kfish2rae5h_kfish2n_fc6t.gff
  #gff.in: gmapn/kfish2rae5h_kfish2a_findcds3t1.gff.gz
  #gff.out: in=37063,put=2512,drop=0,pubaa=0
  #gff.in: gmapn/kfish2rae5h_kfish2a_findcds3ta.gff.gz
  #gff.out: in=78520,put=6510,drop=10191,pubaa=0
  #gff.in: gmapn/kfish2rae5h_kfish2n_findcds3.gff.gz
  #gff.out: in=120755,put=57847,drop=44258,pubaa=30912
  #gff.totout: in=236338,put=66869,drop=54449
  
  #sortGff: 
  $evigene/scripts/bestgenes_update.pl -act=sort -debug -cadd addgene=1 -vers kf2rae5h \
  -conf pubgenome/evigene_kfish2_gbsubmit.conf \
  -in gmapnsub1/kfish2rae5h_kfish2n_fc6t.gff -out gmapnsub1/kfish2rae5h_kfish2n_fc6tans.gff \
   >& gmapnsub1/log.srt$pt

=item try8

env  test=1 vers=fc6t debug=1 ./evgkfish2subgff.pl 
#asmbest.in: gmapn/kfish2rae5h_asmbest_aahogmap.tab2
#asmbest.in: gmapn/kfish2rae5h_altbest.tab
#scaftab.in: pubgenome/ncbifunhe3_kfish2asm.sametab
#pubaa.in: pubtsa6/kfish2rae5x11t6_tsasubmit.aa.gz
#gff.out: gmapnsub1/kfish2rae5h_kfish2n_fc6t.gff
#gff.in: gmapn/kfish2rae5h_kfish2a_findcds3t1.gff.gz
#gff.out: in=37063,put=2512,drop=0,pubaa=0
#gff.in: gmapn/kfish2rae5h_kfish2a_findcds3ta.gff.gz
#gff.out: in=78520,put=6510,drop=10191,pubaa=0
#gff.in: gmapn/kfish2rae5h_kfish2n_findcds3.gff.gz
#gff.out: in=120755,put=57847,drop=44258,pubaa=30912
#gff.totout: in=236338,put=66869,drop=54449

=item more tries 
  in  submitf/evgsubmitkfish2fc7.info

=cut
    

#... subs ........................................................

=item scaffoldTab =  pubgenome/ncbifunhe3_kfish2asm.sametab

0	same	KN805525.1	Scaffold0	6356336	6356336	258514	0
1	same	KN805526.1	Scaffold1	6250690	6250690	296747	0
2	same	KN805527.1	Scaffold2	4057079	4057079	186851	0

13	same	KN805538.1	Scaffold13	790706	790706	25993	0
14	diff	KN805539.1	Scaffold14	808538	808517	45575	-21
15	same	KN805540.1	Scaffold15	822042	822042	60333	0
--
46	same	KN805571.1	Scaffold46	727099	727099	33157	0
47	diff	KN805572.1	Scaffold47	695601	695582	16356	-19
48	same	KN805573.1	Scaffold48	788298	788298	95945	0

=cut

# my(%scaf,%scafn);

sub readScaffoldTab {
  my( $intab)=@_;
  warn "#scaftab.in: $intab\n" if($debug);
  open(F,$intab) or die $intab; 
  while(<F>) {
    next if(/^\W/);
    my($i,$same,$kf3n,$kf2a,$w3,$w2,$nn,$diff)=split;
    $scafn{$kf2a}= $kf3n; $scafn{$kf3n}= $kf2a; #? dont need?
    $scaf{$kf2a}{newna}= $kf3n;
    $scaf{$kf2a}{same}= ($same eq "same")?1:0;
    $scaf{$kf2a}{diff}= $diff;
    $scaf{$kf3n}{size}= $w3;
    $scaf{$kf3n}{oldna}= $kf2a;
  } close(F);
  # return(%scaf,%scafn);
}

=item asmbest.tab

kfish2rae5h_asmbest_aahogmap.tab2
AQueryID	BestHoID	HoDiff	clAAD	cov	pid	GenomeID:gespan	path
Funhe2EKm000003t1	Funhe2EKm000003t1.kf3n	same,0.a,platyfish:ENSXMAP00000008562	OK	100.0	99.9	KN805525.1:5067-10366:-	0
Funhe2EKm000004t1	Funhe2EKm000004t1	top,14.a,mayzebr:XP_004571781.1	SMAL	93.7	97.3	KN805525.1:25802-58390:+	0
Funhe2EKm000005t1	Funhe2EKm000005t1	top,26.a,azmolly:XP_007548349	DIFF0	99.9	92.3	KN805525.1:62577-77945:+	0

kfish2rae5h_altbest.tab
AQueryID	BestRepID	noHoDiff	clAAD	cov	pid	GenomeID:gespan	path
Funhe2EKm000003t1	Funhe2EKm000003t1.kf3n	noho	OK,100%	100	100	KN805525.1:5067-10366:-	0
Funhe2EKm000004t1	Funhe2EKm000004t1.kf2a	noho	DIFF1,99%	91	99	Scaffold0:25799-58393:+	0
Funhe2EKm000004t2	Funhe2EKm000004t2.kf3n	noho	DIFF1,95%	90	95	KN805525.1:39579-58049:+	0

=cut

sub readAsmBestTab {
  my( $intab, $betterh )=@_;
  # do also w/ altbest.tab; same format? use %bestid or 2nd hash? BUT skip bestid entries..
  $betterh ||= undef;
  my %bestid= (); my($ndup,$nerr)=(0,0);
  warn "#asmbest.in: $intab\n" if($debug);
  open(F,$intab) or die $intab; 
  while(<F>) {
    next if(/^\W/ or /^AQuery/);
    chomp; my @v=split"\t";
    my($pubid, $vidsr, $hoval, $clad, $cov, $pid, $genloc, $npath)=@v;
    my($vid,$src)= split /\./, $vidsr,2;
    $src ||= $SPUB; #"pubaa";
    # my $mapsrc=$src; if($src eq $SPUB) { $mapsrc= ($genloc =~ /Scaffold/)? $S2: $S3; }
    my $mapsrc= ($src =~ /$S2|$S3/) ? $src :($genloc =~ /Scaffold/)? $S2: $S3;  
    
    if($skipids{$vid}{$mapsrc}) { $nerr++; next; } # mapsrc/src here?
       
    if(ref $betterh) {
      if(exists $betterh->{$pubid} or exists $betterh->{$vid}) { $ndup++; next; }
      # my $isdup=0; my $bp;
      # $isdup=1 if(($bp= $betterh->{$pubid})); ## and ($bp->{id} eq $vid));
      # $isdup=1 if(($bp= $betterh->{$vid})); ## and ($bp->{id} eq $vid));
      # if($isdup) { $ndup++; next; }
    }  
    
    $bestid{$pubid}{id}=$vid; 
    ## redefine: src == aasrc : S2, S3, SPUB=pubaa; use only mapsrc to match input.gff insrc
    $bestid{$pubid}{src}= $src;
    $bestid{$pubid}{mapsrc}= $mapsrc;
    $bestid{$pubid}{hoval}= $hoval||"noho"; # clad ??
    $bestid{$pubid}{gmap}= "$cov,$pid,$genloc,$npath";
    # includes keep sources per gene locus, need also alt keep rows
  } close(F);
  return (\%bestid); #?
}

# >Funhe2Exx11m011522t4 aalen=418,73%,complete; clen=1711; strand=+; offs=7-1263;
#   pubid=Funhe2EKm003823t1; oid=Fungr1EG3m006328t3; organism=Fundulus_heteroclitus; isoform=1;
#   Name=Arrestin domain-containing protein 3; 
#   Dbxref=CDD:215866,TrEMBL:UniRef50_Q96B67,TrEMBL:ARRD3_HUMAN,family:FISH11G_G1836; 
#   uvcut=Moderate1,17,1712-1728,end3;

# my(%pubaa,%puban);

sub readPubaa {
  my($inaa,$flags)=@_;
  warn "#pubaa.in: $inaa\n" if($debug);
  if($inaa=~/.gz$/) { open(F,"gunzip -c $inaa|") or die $inaa; }
  else { open(F,$inaa) or die $inaa; }
  my($id,$aa);
  while(<F>) {
    if(/^>(\S+)/) { my $xid=$1; # x11 tsa id here
      $pubaa{$id}=$aa if($aa and $id); 
      ($id)= (m/pubid=(\w+)/)?$1:0; $aa="";
      if($id) { # update Name, aalen, Dbxref, clen/strand/offs? oid?
        ## FIXME: Selc: Funhe2EKm035762t1 aalen=238,28%,complete-utrbad,selcstop; Selcstop=1331,1118,962;
        for my $k (qw(aalen offs Dbxref Name Selcstop)) {
          my($v)=m/\b$k=([^;\n]+)/; 
          $puban{$id}{$k}=$v if($v); #? or as string ";$k=$v" ??
        }
      }
    } elsif($id and /\w/) { chomp; $aa.=$_; }
  } close(F);
  
  $pubaa{$id}=$aa if($aa and $id); 
  ## fixme, want pubaa format as for cds2aa?
  ##    aalen=999,xxx,complete;protein=MAZxxxxxVW;cdsoff=111-999;

  # return(\%pubaa,\%puban);
}

=item gmapn/kfish2rae5h_asmbest_aahogmap.tab2

  AQueryID	BestHoID	HoDiff	clAAD	cov	pid	GenomeID:gespan	path
  Funhe2EKm000003t1	Funhe2EKm000003t1.kf3n	same,0.a,platyfish:ENSXMAP00000008562	OK	100.0	99.9	KN805525.1:5067-10366:-	0
  Funhe2EKm000004t1	Funhe2EKm000004t1	top,14.a,mayzebr:XP_004571781.1	SMAL	93.7	97.3	KN805525.1:25802-58390:+	0
  Funhe2EKm000005t1	Funhe2EKm000005t1	top,26.a,azmolly:XP_007548349	DIFF0	99.9	92.3	KN805525.1:62577-77945:+	0
  Funhe2EKm000006t1	Funhe2EKm000006t1.kf3n	top,4.a,platyfish:ENSXMAP00000008658	DIFF0	100.0	99.5	KN805525.1:85415-96983:-	0
=cut

# my(%an,%didput); # pub ann from 1 gff input, use for next..

sub readGff {
  my($ingff,$outh,$flags)=@_;  
  # my $outh=*STDOUT;
  
  warn "#gff.in: $ingff\n" if($debug);
  if($ingff=~/.gz$/) { open(F,"gunzip -c $ingff|") or die $ingff; }
  else { open(F,$ingff) or die $ingff; }
  
  my($drop,$issplit,$id,$gid,$pid,$origid,$swapid,$altnum,$anchange,$nput,$ndrop,$nin,
    $hasan,$insrc,$putout,$putpubaa,$npubaa, $XXputout)= (0) x 29;
  
  # $insrc= ...
  while(<F>) {
    if(/^\W/) { next; }
    
    #gff rows
    my $inrow= $_; # changes only to attr? or src?
    my @v=split"\t"; 
    my($gr,$gs,$gt,$at)=@v[0,1,2,8]; chomp($at);
    $anchange="";
    
    my $SCAFOK=1; my($rnew,$rold);
    if($rnew= $scaf{$gr}{newna}) { 
      $SCAFOK= $scaf{$gr}{same}; #  # problem .. try to shift? need agpLift !
      $rold= $gr; $gr= $rnew;  s/\b$rold\b/$rnew/g; 
      $anchange .= ";scold$SCAFOK=$rold"; 
    }
    
    # problems w/ Split genes..
    
    # problems w/ input "gene" rows .. drop and remake next step ??
    if($gt eq "gene") { next; }
    elsif($gt eq "mRNA") { 
      ($id)= $at=~/\bID=([^;\s]+)/; $nin++;
      $origid = $id; #? or after split/o fixes
      #below# $gid= $id; $gid =~ s/_C\d+$//; $gid=~ s/t\d+$/t1/;
      
      $issplit= $drop= $putout= $putpubaa= $swapid= 0;
      $drop++ if($id=~m/_G\d+$/); # next or not?

      # $issplit= ($id=~s/_C(\d+)$//)?$1:0; # kf3n only; missing this? ** NEED this HERE
      if($id =~ s/_C(\d+)$//){ $issplit=$1; } # kf3n only; missing this? ** NEED this HERE
      elsif($at =~ m/\bSplit=(\w+)/){ $issplit=$1; } # kf2a uses this..
        ## dont trust mapCover annot? new kf3n gmap may not be split
        ## elsif($at =~ m/mapCover=[\d\%]+,Split(\w+)/) { $issplit=$1; } # ;mapCover=100%,Split2/2
        
      ($altnum)= $id =~ m/t(\d+)$/; $altnum ||= 1; # or 0?
      $gid= $id; $gid=~ s/t\d+$/t1/;
      
      # $act= $intab{$id}{action}; # keep/drop/swap-main-alt/other..
      $hasan= ($gs =~ /^($SANN)$/) ?1:0;
      $insrc= ($gs =~ /^($SOUT)$/)? $S3 : ($gs =~ /^($SANN)$/) ? $S2 : "dunno";
      $seenid{$id}{$insrc}="$origid;"; # add actions; allow for dup ids here = Splits ?
     
      ## use both %bestid w/ hovals and %bestalt nohoval .. one hash?
      #** ignore bestaltr when have bestgidr{id}, but not when id ne gid
      # my $bestaltr= $bestalth->{$id}; 
      # if($gid eq $id) { $bestaltr= undef; }
      
      my $Selcstop= $puban{$id}{'Selcstop'} ||""; #? is this right place
      $putpubaa=1 if($Selcstop);

      my $bestgidr= $bestidh->{$gid}; 
      my $checkalt= ($gid eq $id and $bestgidr) ? 0 : 1;
      if($bestgidr) {
        if($bestgidr->{id} eq $id) { 
          $checkalt=0;
          if($insrc eq $bestgidr->{mapsrc}) {
            $putout=1;
            unless($SCAFOK) { $bestgidr->{mapsrc}=$S3; $XXputout=0; } # dont need this hack now, skipids
            if($putpubaa or $bestgidr->{src} eq $SPUB) {
              $putpubaa=1; $anchange.=";bestaa=pubaa"; 
            }
          }
       
          #..........
	        # FIXME src eq pubaa with $S3 is bad map ; use besthash->{mapsrc}  instead
          #     if(($putpubaa or $bestgidr->{src} eq $SPUB)) {
          #       my $mapsrc= $bestgidr->{mapsrc};
          #       if( $insrc eq $mapsrc ) { # was eq $S3 # allow insrc eq S3 or S2 ?
          #         $putout=1; $putpubaa=1;
          #         $anchange.=";bestaa=pubaa"; }
          #     } elsif($bestgidr->{src} eq $insrc) {
          #       $putout=1; 
          #       unless($SCAFOK) { $bestgidr->{mapsrc}=$S3; $XXputout=0; } # dont need this hack now, skipids
          #     }
          
          if($putout and $id ne $gid) { 
            $swapid="$id,$gid";
            $anchange.=";swapmain=$id,$gid";
            ## FIXME: need output id change here? hash like splitgene ?
          } 
        }
        # no else  
        # locus pubid not this id, but this may be valid alt .. need alts in bestalt table..
      }  
      
      # if($checkalt and $bestaltr) 
      if($checkalt) {
        my $bestaltr= $bestalth->{$id}; 
        if($bestaltr and $bestaltr->{id} eq $id) {
          if($insrc eq $bestaltr->{mapsrc}) {
            $putout=1;
            unless($SCAFOK) { $bestaltr->{mapsrc}=$S3; $XXputout=0; } # dont need this hack now, skipids
            if($putpubaa or $bestaltr->{src} eq $SPUB) {
              $putpubaa=1; $anchange.=";bestaa=pubaa"; 
            }
          }
        } # no else
      }
#          #....
#           if(($putpubaa or $bestaltr->{src} eq $SPUB)) {
#             my $mapsrc= $bestaltr->{mapsrc};
# 	          if($insrc eq $mapsrc) { # was eq $S3
#             $putout=1; $putpubaa=1;
#             $anchange.=";bestaa=pubaa"; }
#           } elsif($bestaltr->{src} eq $insrc) {
#             $putout=1; 
#             unless($SCAFOK) { $bestaltr->{src}=$S3; $XXputout=0; } # only for kf2a ?
#           } 
      
      $anchange.=";insrc=$insrc"; # if($putout);

    } else {
      ($pid)= $at=~/\bParent=([^;\s]+)/; 
    }
      
    if($hasan) { 
      # CHANGE: SANN also may be out.gff : kf2a source of intab
      if($gt eq "mRNA") { 
        # ($id)=$at=~/ID=([^;\s]+)/; 
        my $anat= $at;
        $anat=~s/ID=$id;//; # do above now
        $anat=~s/Split=[^;\s]+;//; 
        $anat=~s/;($DROPAN)=[^;\n]+//g; # Target|match|score|ocds|inshort now .. same for all src?
        $an{$id}= $anat; 
      } 
    }

    #after putout# $seenid{$id}{$insrc}="$origid,"; # add actions
    
    if($putout) { # $gs =~ /^($SOUT)$/   ## may be same input gff as SANN
     
      if($gt eq "mRNA") { 
        # ($id)=$at=~/ID=(\w+)/; # above
        #a $issplit= ($id=~s/_C(\d+)$//)?$1:0; #?? missing this? ** NEED to move this up
        #a unless($issplit) {
        #a  if($at =~ m/\bSplit=(\w+)/) {  $issplit=$1; }
        #a  elsif($at =~ m/mapCover=[\d\%]+,Split(\w+)/) {  $issplit=$1; } # ;mapCover=100%,Split2/2
        #a }
        
        # FIXME: add option debug output : input ids, output ids (action val:keep/drop/swapmain/..)
        $drop++ if($id=~m/_G\d+$/); # see above
        $drop++ if(/^NOPATH/); # CHANGE, keep splits;  or $issplit
        $drop++ if($TESTSET and $altnum > 3 and not $swapid=~/$id\b/); # and not $swapid=~/$id\b/
        # ^^ BUGGERS, cancels asmbest_aahogmap main/alt swaps for t>3, n=620/35053; eg: Funhe2EKm000021t1	Funhe2EKm000021t4.kf3n
        $drop++ if($NOSPLITSET and $issplit); # this cuts out large majority of tbl2asn ERROR 
        my $anadd="";

        if(my $serr= $skipids{$id}{$insrc}) { $anchange.=";errskip=$serr"; } # should be removed from bestid hashes.

        ## cds2aa protein: keep for many, replace for pubaa-valid
        if($putpubaa) { 
          ## fixme, pubaa{} now just aa seq; want pubaa format as for cds2aa?
          ##    aalen=999,xxx,complete;protein=MAZxxxxxVW;cdsoff=111-999;
          #  add qual Protein:curated for pubaa cases, for tbl2asn cdsexcept:"annotated by transcript"
          if(my $aa=$pubaa{$id}) {
            my($aalen)= $puban{$id}{aalen} || "missaalen"; #  =~m/aalen=([^;\s]+)/;
            my($offs)=  $puban{$id}{offs} || "missoff"; # =~m/offs=([^;\s]+)/;  $offs||="missoffs";
            $anadd .= "aalen=$aalen;protein=$aa;cdsoff=$offs";
          } else {
            $anadd .= "err=pubaa_missing";
          }
        } else { # cds2aa attr
         my($add1)= $at =~ m/(aalen[^;]+;protein=[^;\n]+;cdsoff=[^;\s]+)/;
         $add1=~s/\*;cdsoff=/;cdsoff=/; # drop cds2aa stop
         $anadd.=";$add1" if($add1);
        }
        
        ## need sub fixAnnot() to handle these various annot updates    
        my $an= $an{$id}||""; 
        ## *should always have an, warn if not..
        warn "#errnoan: $id, origid:$origid\n" if($debug and not $an);
        
        my($anoid)= $an =~ m/\boid=([^;\s]+)/;
        my($inoid)= $at =~ m/\boid=([^;\s]+)/;
        
        #?? need both anoid,inoid ?? skip inoid
        ## oid=Funhe2Exx11m002607t4,kf2rae5alt:Funhe2Exx11m002607t4 << stutter
        # my %oid= map{ $_ =>1 } split",", "$anoid,$inoid";
        my %oid= map{ $_ =>1 } split",", $anoid;
        
        my @oid= grep /\w/, sort keys %oid;
        my($trasmid)= grep /$TRAMSID/,@oid;        
        if($trasmid) { 
          $trasmid=~s/^\w+://; 
          $anadd .= ";tblan=trid:$trasmid";  # annot for tbl2asn TSA valid trasm
        }
        # remove some excess, 2ndary  KEEPOID
        @oid= grep /$KEEPOID/,  @oid;
        # @oid= grep{ not m/($DROPOID)/ } @oid;
        my $oid= join",",@oid;
        unless($an =~ s/;oid=[^;\s]+/;oid=$oid/) { $an.=";oid=$oid"; }
        
        ## below..         
        #   s/;($DROPAN)=[^;\n]+//g; 
        #   for my $ak (qw(aalen protein cdsoff)) {
        #     s/;$ak=[^;\s]+// if($anadd=~/$ak=/);
        #     $an =~ s/;$ak=[^;\s]+//;
        #   }
       
        ##below# $an.= ";$anchange" if($anchange);
        if($puban{$id}) { 
          ## puban subkeys= aalen offs Dbxref Name
          ## fix an/puban: genegroup=FISH11G and in Dbxref=..,family:FISH11G : drop one?
          my @pk= sort keys %{$puban{$id}};
          #already in $anadd# unless($putpubaa) { @pk= grep{ not m/aalen|offs/ } @pk; }  
          for my $pk (@pk) { $an =~ s/\b$pk=[^;\n]+[;]?//; $an.= ";$pk=".$puban{$id}{$pk}; } 

          my $Selcstop= $puban{$id}{'Selcstop'} ||""; #? is this right place
          if($Selcstop) { # add /transl_except=(pos:1002..1004,aa:Sec);
            my @ss=split",",$Selcstop;
            for my $sb (@ss) { my $se=$sb+2; $an .=";transl_except=(pos:$sb..$se,aa:Sec)" if($sb>0); }
            #NOT HERE# $putpubaa=1 if($Selcstop);
           }
        }
        
        # if($didput{$id}) { } # ** check dup recs, but splits have 2 mrna same id..        
        if($issplit) { 
        
          unless($an){ $an=$at; $an=~s/;($DROPAN)=[^;\n]+//g;  }  #? keep $anadd="";
          $an="Split=$issplit;$an"; 
 
          if($SplitGeneFixID) {  
            # give new ID=oldid_C$issplit 
            # * update attr ;gene=gid  with _C$issplit ..
            my $sid= $id.$IDSplitSuffix.$issplit;
            $origid= $id;
            # s/ID=$origid/ID=$sid/; # leave to below ??
            $an =~ s/ID=$origid/ID=$sid/; 

            my($mgid)= $an =~ m/;gene=([^;\s]+)/;
            # problem: only some alts are split, eg t3 split, t1,t2 not. gid for t1,t2 and t3.s1 should be same
            if($mgid) { 
               if($issplit > 1) { # keep geneid of s=1 unchanged?
                my $sg= $mgid.$IDSplitSuffix.$issplit;
                $an =~ s/;gene=$mgid/;gene=$sg/; 
		          }
             }

            $id= $sid; # change here? will affect id hashes below..
            if(my $np= $didput{$id}) { $anchange.= ";errdupid=$np"; $drop++; }
          }         
          elsif(my $np= $didput{$id} > 1) { $anchange.= ";errdupid=$np"; $drop++; }
          
        } else {
          if(my $np= $didput{$id}) { $anchange.= ";errdupid=$np"; $drop++; }
        }

        if($putpubaa) {  # qual Protein:curated > validated
          if($an =~ /quality=([^;\s]+)/) { my $q=$1; my $r=$q; 
            unless($r=~s/Protein:/Protein:validated-/) { $r.=",Protein:validated"; }
            $an=~s/quality=$q/quality=$r;/;
          } else {
            $anadd.=";quality=Protein:validated";
          }
        }
        
        # dupid bug: $didput{$id} drop dupl, but make sure drop correct one, dupid should not occur..
        if($anchange =~ /errdupid=/) { 
          $drop++;  #n=2393 for fc6t
          warn "#errdupid: $id, changes=$anchange\n" if($debug);
        }

        # s/;($DROPAN)=[^;\n]+//g; # dont care, dropping..
        $an =~ s/;offs=[^;\s]+// if($anadd=~/cdsoff=/);
        for my $ak (qw(aalen protein cdsoff)) {
          $an =~ s/;$ak=[^;\s]+// if($anadd=~/$ak=/);
          ## s/;$ak=[^;\s]+// if($anadd=~/$ak=/);
        }
       
        #fix an/puban: genegroup=FISH11G and in Dbxref=..,family:FISH11G : drop one?
        $an.= ";$anchange" if($anchange);
        $an.= ";$anadd" if($anadd);  
        
        ## attr=paralog=FISH11G_G24467,,,,,,,,; what is this stutter?
        if($an) { # always have ?
          $an=~s/;;+/;/g; $an=~s/,,+/,/g; $an=~s/^;//; 
          unless(s/\bID=$origid.*$/ID=$id;$an;/) { # SplitGeneFixID fix
            warn "#erranmiss: $id, origid:$origid, an=$an\n" if($debug);
          }
        } 
       
      } else { 
        ($pid)=m/Parent=([^;\s]+)/; 
        ## SPLIT ID problem wrong _C1 assign to exons _C2
        my $npid= $id;
        if($issplit) {
          if(m/;Split=(\w+)/) { my $sn=$1;
            if($SplitGeneFixID) {
              my $sid= $origid.$IDSplitSuffix.$sn;
              if($sid ne $npid and $pid=~/^$origid/) { $npid= $sid; }
            }
          } else { 
            s/Parent=$pid/Parent=$pid;Split=$issplit/; 
          }
        }
        s/;($DROPAN)=[^;\n]+//g; s/Parent=$pid/Parent=$npid/; 
      }
      
=item split id bug

  .. problem because input gff has re-sorted split genes, 2 mRNA first, then all exons ..
  .. cant rely on last mRNA of split having right id, but use exon Split=1/2 annot? not on CDS tho
KN805525.1	kf2rae5	mRNA	5883935	5890995	89	-	.	ID=Funhe2EKm000203t1_C2;Split=2; insrc=kf2a;
KN805525.1	kf2rae5	mRNA	5888005	5890994	10	-	.	ID=Funhe2EKm000203t1_C1;Split=1;
KN805525.1	kf2rae5	exon	5883935	5884183	98	-	.	Parent=Funhe2EKm000203t1_C1;Split=2 <<< bug
KN805525.1	kf2rae5	exon	5887883	5888017	100	-	.	Parent=Funhe2EKm000203t1_C1;Split=2
KN805525.1	kf2rae5	exon	5888147	5888215	100	-	.	Parent=Funhe2EKm000203t1_C1;Split=1
   CDS also bad
   
=cut
   
      if($drop) {
        $ndrop++ if($gt eq "mRNA");
      } else {
        print $outh "##gff-version 3\n" unless($didhead++);
        ## change col2=SRC to SRCOUTmain and  SRCOUTalt ?
        
        print $outh "\n" if($gt eq "mRNA");
        print $outh "#SErr." unless($SCAFOK); # got 1018 of these
        print $outh $_;  # CHANGE, refer to intab for keep/drop; ?? use $inrow instead of $_ ??
        if($gt eq "mRNA" and $SCAFOK) { $didput{$id}++; $nput++; $npubaa++ if($putpubaa);  }
      }
    } 
    
    if($gt eq "mRNA") { 
      $anchange =~ s/;insrc=$insrc//;
      $seenid{$id}{$insrc} .= "ok:$putout/$drop$anchange,";  # add actions
    }
  } close(F);
  
  warn "#gff.out: in=$nin,put=$nput,drop=$ndrop,pubaa=$npubaa\n" if($debug);
  return($nin,$nput,$ndrop);
}



sub sortGff {
  my($otmpgff, $outsrtgff)= @_;
  
  my $cmd="$EVIGENE/scripts/bestgenes_update.pl -act=sort -debug -cadd addgene=1 "
    ."-vers $SOUT -conf $SUBMITCONF -in $otmpgff -out $outsrtgff";
  warn "#sortGff: $cmd\n"; ## if($debug);
  #debug# system($cmd);
}

=item sortGff
 
  #sortGff: 
  $evigene/scripts/bestgenes_update.pl -act=sort -debug -cadd addgene=1 -vers kf2rae5h5 \
  -conf pubgenome/evigene_kfish2_gbsubmit.conf \
  -in gmapnsub1/kfish2rae5h_kfish2n_fc5.gff -out gmapnsub1/kfish2rae5h_kfish2n_fc5ans.gff \
   > gmapnsub1/log.fc5srt 2>&1

  # genesort: mRNA duplicate for Funhe2EKm012296t2 ; one has Split=2; 2nd not.
  
  # scaffold splits.. *?? need to re-sort by scaffold for ncbi tbl2asn ? or not?
  # .. non-scaff/local C1 splits .. problem having 2 mRNA for overlapped spans .. revise to 1 mrna?

    insrc=kf3n;err=pubaa_missing; Name=TE:;oid=kf1gene9f:Funhe5EG035707t1,AUGie2:AUGpie2s3903g35t1;
    ^^ this is TE mapping mess. should drop it.
  KN805525.1  kf2rae5h  gene    2888508 2889451 76   -   ID=Funhe2EKm007172;Split=1
  KN805803.1  kf2rae5h  gene    262221  273731  76   -   ID=Funhe2EKm007172;Split=2
  KN805525.1  kf2rae5h  mRNA    2888508 2889451 17   -   ID=Funhe2EKm007172t1;gene=Funhe2EKm007172;Split=1;quality=Class:Medium,Express:Medium,Homology:OrthologMedium,Intron:Weak,Map:Strong,Protein:complete;aaSize=178;cdsSize=22%,204/1858;Name=TE: Transposase (41%C);oname=na;ortholog=catfish:IctpunEGm005529t1,41%;genegroup=FISH11G_G16309;Dbxref=UniProt:E7F105_DANRE,na;intron=25%,1/4;express=28%,rx:grandis,2a/0.9e/3g;mapCover=100%;equiv1=Funhe5EG035707t1/100,Funhe5EG016800t1/32,Funhe5EG035708t1/36,Funhe5EG016263t1/43,;oid=kf1gene9f:Funhe5EG035707t1,AUGie2:AUGpie2s3903g35t1;cxlen=201/1858;bestaa=pubaa;insrc=kf3n;err=pubaa_missing;
  KN805803.1  kf2rae5h  mRNA    262221  273731  76   -   ID=Funhe2EKm007172t1;gene=Funhe2EKm007172;Split=2;quality=Class:Medium,Express:Medium,Homology:OrthologMedium,Intron:Weak,Map:Strong,Protein:complete;aaSize=178;cdsSize=22%,204/1858;Name=TE: Transposase (41%C);oname=na;ortholog=catfish:IctpunEGm005529t1,41%;genegroup=FISH11G_G16309;Dbxref=UniProt:E7F105_DANRE,na;intron=25%,1/4;express=28%,rx:grandis,2a/0.9e/3g;mapCover=100%;equiv1=Funhe5EG035707t1/100,Funhe5EG016800t1/32,Funhe5EG035708t1/36,Funhe5EG016263t1/43,;oid=kf1gene9f:Funhe5EG035707t1,AUGie2:AUGpie2s3903g35t1;cxlen=201/1858;bestaa=pubaa;insrc=kf3n;err=pubaa_missing;
  
  KN805525.1  kf2rae5h  exon    2888508 2888822 85   -   Parent=Funhe2EKm007172t1;Split=1;ix=3
  KN805525.1  kf2rae5h  CDS     2888509 2888822 .    -   Parent=Funhe2EKm007172t1;Split=1
  KN805525.1  kf2rae5h  exon    2889243 2889278 100  -   Parent=Funhe2EKm007172t1;Split=1;ix=2
  KN805525.1  kf2rae5h  CDS     2889243 2889278 .    -   Parent=Funhe2EKm007172t1;Split=1
  KN805525.1  kf2rae5h  exon    2889333 2889451 97   -   Parent=Funhe2EKm007172t1;Split=1;ix=1
  KN805525.1  kf2rae5h  CDS     2889333 2889450 .    -   Parent=Funhe2EKm007172t1;Split=1
  
  KN805803.1  kf2rae5h  exon    262221  262552  100  -   Parent=Funhe2EKm007172t1;Split=2;ix=3
  KN805803.1  kf2rae5h  exon    262972  264430  100  -   Parent=Funhe2EKm007172t1;Split=2;ix=2
  KN805803.1  kf2rae5h  CDS     263068  263271  .    -   Parent=Funhe2EKm007172t1;Split=2
  KN805803.1  kf2rae5h  exon    273667  273731  100  -   Parent=Funhe2EKm007172t1;Split=2;ix=1

=cut

# my(%mapq,%cdsq,%mapv,%cdsv);  # for makeAltBestTab

sub bests { # for makeAltBestTab
  my($outh,$td,$last)= @_; 
  my $did=0;
  my $m2= $mapq{$td}{$S2}; my $m3= $mapq{$td}{$S3};
  if(($m2 and $m3) or $last) {
    my $cd2=$cdsq{$td}{$S2}||0; my $cd3=$cdsq{$td}{$S3}||0;
    my $mv2=$mapv{$td}{$S2}||0; my $mv3=$mapv{$td}{$S3}||0;
    if($mv2=~/ERR/ and $mv3=~/ERR/) { }
    elsif($mv2=~/ERR/) { $m2=$cd2=-9; }
    elsif($mv3=~/ERR/) { $m3=$cd3=-9; }

    my $bsrc= ($m2 > $m3 and $cd2 > $cd3) ? $S2 
         : ($m3 > $m2 and $cd3 > $cd2) ? $S3
         : ($cd2 > $cd3 and $m2 > 0.8*$m3) ? $S2 
         : ($cd3 > $cd2 and $m3 > 0.8*$m2) ? $S3 
         : $S3;
    my $cdv= $cdsv{$td}{$bsrc}||"misscd"; my $cdq= $cdsq{$td}{$bsrc}||0; $cdv.=",$cdq%";
    my $mapv=$mapv{$td}{$bsrc}||"missmapv"; my $mq= $mapq{$td}{$bsrc}||0; 

    #^^ redo call, pay attention to AADIFF qual, use pubaa where cds2aa is short/poor match
    if($cdq < $AADIFFMIN) { $bsrc=$SPUB; }       
    my $p=""; $p="#e." if($mapv =~ /ERR/); #skip this one??
    print $outh join("\t",$p.$td,"$td.$bsrc","noho",$cdv,$mapv)."\n"; 
    $did=1; # $did{$td}++;
    }
  return $did;
}

=item makeAltBestTab ERRmap alts

n=856  try2, excepting C3 splits, '.' strands
n=2401 : how many need check for keeping, any?
  281 split .. how many are Split problems : 159 C3: scaf split, 94 C1: split, 28 C2:
  .. ignore split C3 cases, not maperr.
  2120 non split;
  
  ** for poor maps, strand/ERRrev can be wacky/nonsense. 
    n=1016 pubaa; if pubaa can we flip strand to main?
  -- are any of these case of main/t1 wrong, all alts right?
  -- t1 split cases
  15 Funhe2EKm017943 : t1 C3:split KN807384.1=Scaffold1861 KN807896.1=Scaffold2382 ;
      ^^ alts on many diff scafs KN807422.1/Scaffold1899  KN806695.1/Scaffold1170 .. should be several loci
      ** AND asmbest_aahogmap sets t1 = t2.kf3n on 3rd scaf ** bad
asmbest_aahogmap: 
Funhe2EKm017943t1	Funhe2EKm017943t2.kf3n	same,0.a,azmolly:XP_007564356	OK	97.2	80.1	KN807384.1:9382-22382:-	C3:KN807896.1:14687-15904,-,33.7
altbest.tab     
Funhe2EKm017943t1	Funhe2EKm017943t1.pubaa	noho	DIFF1,74%	92	99	Scaffold2382:14684-15835:-	C3:Scaffold1899:6401-11515,-,92%	noer
Funhe2EKm017943t2	Funhe2EKm017943t2.kf3n	noho	OK,100%	100	100	KN807422.1:1566-18266:-	0	ERRscaf,
 
  14 Funhe2EKm003721
  12 Funhe2EKm020972  : ERRrev because t1 has '.' strand, others same scaf '+' or '.'
  11 Funhe2EKm009219
  11 Funhe2EKm015579
  11 Funhe2EKm015898
  10 Funhe2EKm021051
  10 Funhe2EKm032345

grep ERR gmapn/kfish2rae5h_altbest.tab | grep -v ERRnomap | cut -f9 | sort | uniq -c     
1017 ERRrev,
  77 ERRrev,ERRscaf,
 635 ERRrev,ERRscaf,ERRspan,
 104 ERRrev,ERRspan,
  70 ERRscaf,   ?? should these be other loci
 379 ERRscaf,ERRspan,
 119 ERRspan,

=item rev makeAltBestTab

## FIXME: revise makealtbest/ altbest.tab to combine asmbest_aahogmap.tab
## .. and correct asmbest cases where t1 mismatch to talt and poor map, pubaa/kf3n should be pubaa/kf2a

  add input: my $bestidtab="gmapn/kfish2rae5h_asmbest_aahogmap.tab2";

perl -ne '@v=split"\t"; ($pd,$tds,$hov,$aad,$cov,$pid,$gloc,$np)=@v;
($td,$sr)=split/\./,$tds; $sr||="pubaa";
$maps=($gloc=~/Scaffold/)?"kf2a":"kf3n"; if($ARGV=~/altbest/) {
$td{$td}{src}=$sr; $td{$td}{maps}=$maps; $pd{$pd}=$td;
$td{$td}{val}=[@v]; } else { $amaps=$td{$td}{maps};
$aval=$td{$td}{val}; if($sr eq "pubaa" and $maps ne $amaps) { print $_;
print join("\t",@$aval)."\n"; } } ' \
kfish2rae5h_altbest.tab kfish2rae5h_asmbest_aahogmap.tab2 \

-- update bestaaho.tab3

perl -ne 'my @v=split"\t"; chomp($v[-1]); my($pd,$tds,$hov,$aad,$cov,$pid,$gloc,$np)=@v;
($td,$sr)=split/\./,$tds; $maps=($gloc=~/Scaffold/)?"kf2a":"kf3n"; 
if($ARGV=~/altbest/){ $td{$td}{src}=$sr; $td{$td}{maps}=$maps; $td{$td}{val}=[@v]; } 
elsif(/^AQueryID/) { print; }
else { $amaps=$td{$td}{maps}; $aval=$td{$td}{val}; 
if($maps ne $amaps) { @vn= @$aval; $vn[0]= $v[0]; $vn[2]= $v[2]; $vn[8].=",upd2"; @v=@vn; }
if(not $sr or $sr eq "pubaa") { unless($v[1]=~s/\.\w+/.pubaa/) { $v[1] .= ".pubaa"; } }
print join("\t",@v)."\n"; } ' \
kfish2rae5h_altbest.tab kfish2rae5h_asmbest_aahogmap.tab2 \
  >  kfish2rae5h_asmbest_aahogmap.tab3

 ----
 
=cut

sub makeAltBestTab {
  my($outf,$outh)= @_;

  # $inaadiff= kfish2rae5h{t1,ta}_asm2a2n.diffaa.tab   
  # $inmapatt= kfish2rae5hpub-kfish2n.map.attr kfish2rae5h_asm2.map.attr
  # kfish2rae5h{t1,ta}_asm2a2n.diffaa.tab kfish2rae5hpub-kfish2n.map.attr  kfish2rae5h_asm2.map.attr \
  #  > kfish2rae5h_altbest.tab
  my @inaadiff= qw(gmapn/kfish2rae5ht1_asm2a2n.diffaa.tab gmapn/kfish2rae5hta_asm2a2n.diffaa.tab);
  my @inmapattr=qw(gmapn/kfish2rae5hpub-kfish2n.map.attr gmapn/kfish2rae5h_asm2.map.attr);

# GeneXrefStrandProblem : alts rev or to main/t1 .. drop here ???
# -- need to check these some how, drop alts on rev strand, else call as new gene locus if evidence/homol good?
  
  ## replace cdsp qual score w/ simple diff in aaNoGap sizes
  my @inaaqual= qw(pubgenes/kfish2rae5h.main.pub.aa.qual pubgenes/kfish2rae5h.alt.pub.aa.qual
    gmapn/kfish2rae5h_kfish2a_findcds3t1.cds2aa.qual gmapn/kfish2rae5h_kfish2a_findcds3ta.cds2aa.qual 
    gmapn/kfish2rae5h_kfish2n_findcds3.cds2aa.qual
    );
  
  unless($outh and ref $outh) {
    $outf ||= "kfish2rae5h_altbest.tab"; #? $outh= *STDOUT;
    warn "#malt.out: $outf\n" if($debug);
    rename($outf,"$outf.old") if(-s $outf);
    open($outh,'>',$outf) or die $outf;
    }
    
  # $S2="kf2a"; $S3="kf3n";
  print $outh join("\t",qw(AQueryID BestRepID noHoDiff  clAAD cov pid GenomeID:gespan path))."\n";  
  
  for my $inf (@inaadiff) {
    warn "#malt.in: $inf\n" if($debug);
    open(F,$inf) or die $inf;
    while(<F>) {
      next if(/^\W/ or /^AQueryID/); chomp; my @v=split"\t";
      my($aacl,$td,$cdsq,$src)=@v; 
      my($cdsp)= $cdsq=~/^(\d+)/; $cdsp||=0; $cdsp=100 if($cdsp>100);
      $aacl=~s/DIFF.0/DIFF0/; $aacl=~s/DIFF..*/DIFF1/; $aacl=~s/\..*$//;
      $cdsq{$td}{$src}=$cdsp; $cdsv{$td}{$src}= $aacl; 
    } close(F);
  }
  
  my(%awv);
  for my $inf (@inaaqual) {
    warn "#malt.in: $inf\n" if($debug);
    open(F,$inf) or die $inf;
    my $src=($inf=~/pub.aa.qual/)? $SPUB
      : ($inf=~/kfish2n/)? $S3 : ($inf=~/kfish2a/)? $S2 : $S3;
    while(<F>) {
      next if(/^\W/ or /^AQueryID/); chomp; my @v=split"\t";
      my($td,$aw)= @v; $awv{$td}{$src}= $aw;
      if($src ne $SPUB and (my $awp= $awv{$td}{$SPUB}) ) {
         my $cdsp= int(0.5 + 100* $aw/$awp);  $cdsq{$td}{$src}= $cdsp;
         my $aacl= $cdsv{$td}{$src};
         $aacl=~s/^(\d+)%/$cdsp%/;
         $cdsv{$td}{$src}= $aacl; 
         }
    } close(F);
  }

  my (%tdin, %nomap, %did, %gloc);
  for my $inf (@inmapattr) {
    warn "#malt.in: $inf\n" if($debug);
    open(F,$inf) or die $inf;
    ## got some dup input ids, from multimap .. skip 2nd
    while(<F>) {
      next if(/^\W/ or /^AQueryID/); chomp; my @v=split"\t";
      my($td,$cov,$pid,$nexon,$gloc,$path)= @v; 
      map{ s/%//;  $_= int(0.5 + $_); } ($cov,$pid);
      $tdin{$td}++;
      next if($did{$td}); #?? check for better score or not?
      my $src="none";
      if($gloc =~/^NOPATH/) { # or $cov == 0
        $nomap{$td}= $gloc; 
        $src= ($pid eq "1" and $nexon eq "0/0")? $S3 : ($pid eq "0" and $nexon eq "0")? $S2:"nulls"; # and $cov eq "0%"
        $mapq{$td}{$src}= -1 unless($mapq{$td}{$src}); 
        $mapv{$td}{$src}= join("\t",0,0,$gloc,$path,"ERRnomap");
      } else {
        my $mapq=0;
        $src= ($gloc =~ /^Scaffold/i)? $S2 :($gloc =~ /^KN|^J/)? $S3 :"nulls";
        my $err="";
        
        ## add location checks: alt rev t1? alt way outside span of t1 ? are they sorted by altnum?
        # Scaffold0:25799-58393:+ ; KN805525.1:5067-10366:-
        my $gid=$td; $gid=~s/t\d+/t1/;
        my $mvmain= $mapv{$gid}{$src};
        if($mvmain) {
          # FIXME here, sometimes mvmain/t1 is wrong vs other alts .. esp from asmbestho.tab ?
          my($cov1,$pid1,$loc1,$path1,$err1) = split"\t",$mvmain;
          my($mr,$mb,$me,$mo)= split /[:-]/, $loc1, 4; $mo=0 if($mo eq '.');
          my($gr,$gb,$ge,$go)= split /[:-]/, $gloc, 4; $go=0 if($go eq '.');
          $err.="ERRscaf," if($mr ne $gr);
          $err.="ERRspan," unless($gb < $me and $ge > $mb);
          $err.="ERRrev," if($mo and $go and $mo ne $go);
          if($err and ($path =~ /C3:/ or $path1 =~ /C3:/)) { $err=~s/ERR/split/g; }
        }
        
        if($path =~ /^C(\d):/) { my $sc=$1;  ## BAD split values for S2
          my @sv= split/[:,]/,$path; my $svc= ($sc > 2) ? $sv[4] : $sv[3]; $svc=~s/%//;
          if( $src eq $S2) { $cov=$svc if($svc>0); $mapq=$cov; $path=~s/,Split.*//;  # fix bad
          } else { $mapq = $cov - $svc; }
        } else {  
          unless($path eq "0" or $path eq "1") { next if($mapq{$td}{$src}); }
          $mapq= $cov; ## * $pid/100; ## pid == 0 for all kf2a, missing data
        }
        $mapq=100 if($mapq>100); $mapq{$td}{$src}= $mapq;
        $pid=99 if($pid==0);  $cov=100 if($cov>100);
        ##map{ $_= int(0.5 + $_) } ($cov,$pid);
        $err ||= "noer";
        $mapv{$td}{$src}= join("\t",$cov,$pid,$gloc,$path,$err);
        # $gloc{$td}{$src}= $gloc;
      }
      my $osrc=($src eq $S2)?$S3:$S2;
      $did{$td}= bests($outh,$td,0) if($mapq{$td}{$osrc});
  
    } close(F);
  }

  #END  
  for my $td (grep{ not $did{$_} } sort keys %tdin) { 
    if(exists $mapq{$td}) { $did{$td}= bests($outh, $td,1); } # this gets all nomap ?
    else{ my $nom=$nomap{$td}||"missmap"; 
      print $outh join("\t",$td,$td,"noho",0,0,0,$nom,0)."\n"; $did{$td}++; }
  }  
  
  # close($outh);
} # makeAltBestTab


__END__
   
#   gunzip -c gmapn/kfish2rae5h_kfish2a_findcds3{t1,ta}.gff.gz \
#             gmapn/kfish2rae5h_kfish2n_findcds3.gff.gz | \
#     cat gmapn/kfish2rae5h_asmbest_aahogmap.tab other.tabs - | \
#   env sann='kf2rae5|kf2rae5alt' sout=kf2rae5h perl -ne \
#    .. perl below needs update ..
#    ..
#    > tempan.gff  
   
#    cat tempan.gff | $evigene/scripts/bestgenes_update.pl -act=sort -debug -cadd addgene=1 \
#    -vers kf2rae5h -conf pubgenome/evigene_kfish2_gbsubmit.conf \
#    -in stdin -out gmapnsub1/kfish2rae5h_kfish2n_fc4ans.gff



=item steps:

S1: prepare inputs, new tables
  
  which gff inputs here? gmapn/findcds.gff or gmapn/kfish2rae5hpub-kfish2n.gmap.gff
    .. findcds versions have some improvements to cds locs; source for cds2aa valid prots
    
  1. kfish2rae5h_kfish2a_findcds3{t1,ta}.gff from kfish2rae5h.*pubsort.gff
    - kfish2rae5h*.pubsort.gff for orig mapping, need for cases of kf3n mismap
    * this set has all the pub annots needed for output, both pubsort and findcds
    * use kf2a gff locs for tabled cases of kf2a map is best
    * need to remap kfish2a Scaffold to kfish2n NCBI scafIDs, see
      pubgenome/ncbifunhe3_kfish2asm.sametab
    
  2. kfish2rae5h_kfish2n_findcds3.gff from kfish2rae5hpub-kfish2n.gmap.gff
    * use these kf3n gff locs for tabled cases of kf3n map is best, and default
  
  ^^ merge 1,2 gff first?, using tables of best mapping + best cds2aa ??
    - where asmbest_aahogmap doesnt have best call, use pct map + longest cds2aa
    - or is that Step2

  3. kfish2rae5h_asmbest_aahogmap.tab : table of locus best rep of gmap kf3n, kf2a or poormap aa exception
    - n=35053 gene loci, most/all?; n=19974 are cds2aa ok; n=15079 need poormap pubaa exception
    - ?* alts need cds2aa ok table, but need new blastp for that? or fudge it?
    - see below, update asmbest_aahogmap.tab ?  replace poor pAAD,clAAD w/ refho align stats?
    - some of small DIFF pubaa best can be recalled as cds2aa ok .. calc from blast align/ident diff < 1-9aa ?
      
    AQueryID                BestHoID                pAAD    cov     pid     clAAD   GenomeID:gespan:geor    path
    Funhe2EKm000003t1       Funhe2EKm000003t1.kf3n  100     100.0   99.9    OK      KN805525.1:5067-10366:- 0
    Funhe2EKm000004t1       Funhe2EKm000004t1       43      93.7    97.3    SMAL.227,-227   KN805525.1:25802-58390:+        0
                                ^^ make translate, poormap exception
    Funhe2EKm000005t1       Funhe2EKm000005t1.kf3n  97      99.9    92.3    DIFF.0  KN805525.1:62577-77945:+        0
    Funhe2EKm000009t1       nohoref                 100     100.0   99.9    OK      KN805525.1:121590-124252:-      0
    Funhe2EKm000021t1       Funhe2EKm000021t4.kf3n  90      99.9    88.7    DIFF.0  KN805525.1:297161-306776:+      0
                                ^^^ swap t1,t4 for submitset
    Funhe2EKm000032t1       Funhe2EKm000032t1.kf2a  92      75.6    99.4    DIFF.0  KN805525.1:1133312-1136423:+    0
                                ^^^ use kf2a mapping
    Funhe2EKm043040t1       nohoref                 0       0       1       MISS    NOPATH:1-69:.   0/0
                                ^^^ genome MISS, fix 'nohoref' wrong now

S2:
  * input annot.gff + in tables? to output, annotated, sorted by gene/alts, ans.gff
  - in table remapping alt to main for best cds2aa gmap cases?
  - in table of pubaa best (cds2aa bad) exceptions? need also pub.aa inputs to stick into out.gff?
  
  gunzip -c gmapn/kfish2rae5h_kfish2a_findcds3{t1,ta}.gff.gz \
            gmapn/kfish2rae5h_kfish2n_findcds3.gff.gz | \
    cat gmapn/kfish2rae5h_asmbest_aahogmap.tab other.tabs - | \
  env sann='kf2rae5|kf2rae5alt' sout=kf2rae5h perl -ne \
   .. perl below needs update ..
   ..
   > tempan.gff  
   
   cat tempan.gff | $evigene/scripts/bestgenes_update.pl -act=sort -debug -cadd addgene=1 \
   -vers kf2rae5h -conf pubgenome/evigene_kfish2_gbsubmit.conf \
   -in stdin -out gmapnsub1/kfish2rae5h_kfish2n_fc4ans.gff

  #old.in gmapn/kfish2rae5h.{main,alt}.pubsort.gff.gz gmapn/kfish2rae5h_kfish2n_findcds3.gff.gz

 
S3: ans.gff to tbl annots, pep for exceptions only? or all to check cds2aa

 s3a: make goodtab ? or not? 
      or drop.tab to list genes, alts too poorly mapped to include?
 
  ./evigene2genotbl_kfish2.pl -notbl2asn -skipnoann -debug \
   -conf pubgenome/evigene_kfish2_gbsubmit.conf \
   -in   gmapnsub1/kfish2rae5h_kfish2n_fc4ans.gff  -proteins gff \
   -change gmapnsub1/kfish2rae5h_kfish2n_fc4ans.goodtab \
   -out  gmapnsub1/kfish2rae5h_kfish2n_fc4.tbl > gmapnsub1/log.evg2tblfc4 2>&1


S4: tbl2asn

  submitf/gmapnsub1/
  pt=kfish2rae5h_kfish2n_fc4
  ln -s ../pubgenome/ncbifunhe302scaf.fa $pt.fa
  rename ".all" "" $pt.all.{tbl,pep}
  
  $nbin/tbl2asn -Mn -Vbv -a r10k -XE -Hy -t ./kfish2rae5gsub.sbt  \
  -j '[moltype=DNA] [organism=Fundulus_heteroclitus]'   \
  -Z $pt.discr -i $pt.fa -f $pt.tbl > log.asn$pt 2>&1

#.....................

=item make kfish2rae5h_altbest.tab

  Make $bestalttab="gmapn/kfish2rae5h_altbest.tab";
  #^vv^ make 2nd table for others, alts, etc. best of kf3n, kf2a gmap and/or cds2aa size
  similar to $bestidtab="gmapn/kfish2rae5h_asmbest_aahogmap.tab2";
  
  inputs: ?
  gmap quals: 
    kfish2rae5h_asm2.align.tab.gz : incomplete, enough? from pubgenes/attr.tab                
    kfish2rae5hpub-kfish2n.align.tab.gz
  -- or map.attr subset: 
    kfish2rae5h_asm2.map.attr kfish2rae5hpub-kfish2n.map.attr
  cds2aa quals:
    kfish2rae5h{t1,ta}_asm2a2n.diffaa.tab   ? these may be complete for alts
  or sources:
    kfish2rae5h_kfish2a_findcds3 t1,ta, kfish2rae5h_kfish2n_findcds3 .diffaaok = cds2aa compare outputs
  
  add this to cancel kf2a for changed scafs:
   my $scafidtab="pubgenome/ncbifunhe3_kfish2asm.sametab";

  result1:  26387 kf2a, 83776 kf3n;

# path: if($sc > 2) { ($scc,$sr,$sbe,$so,$scv)=@sv; } else { ($scc,$sbe,$so,$scv)=@sv; }
## bsrc="pubaa" when? for bad cdsq ?
## fixme: both src have asm2a2n.diffaa.tab vals, pick best: OK > DIFF0 > DIFF|SMAL|BIG > MISS

perl -ne  \
'next if(/^\W/ or /^AQueryID/); chomp; my @v=split"\t";
if($v[1] =~ m/^Funhe/ and $v[3]=~/($S2|$S3)/) { 
 my($aacl,$td,$cdsq,$src)=@v; ($cdsp)= $cdsq=~/^(\d+)/; $cdsp||=0;
 $aacl=~s/DIFF.0/DIFF0/; $aacl=~s/DIFF..*/DIFF1/; $aacl=~s/\..*$//;
 $cdsq{$td}{$src}=$cdsp; $cdsv{$td}{$src}= $aacl; 
 } 
elsif(/^Funhe/) { 
  my($td,$cov,$pid,$nexon,$gloc,$path)= @v; map{ s/%// }($cov,$pid);
  $tdin{$td}++;
  if($gloc =~/^NOPATH/) { # or $cov == 0
    $nomap{$td}= $gloc; 
    $src= ($pid eq "1" and $nexon eq "0/0")? $S3 :($pid eq "0" and $cov eq "0%")? $S2:"nulls";
    $mapq{$td}{$src}= -1 unless($mapq{$td}{$src}); 
  } else {
    $src= ($gloc =~ /^Scaffold/i)? $S2 :($gloc =~ /^KN|^J/)? $S3:"nulls";
    if($path =~ /^C(\d):/) { $sc=$1;  ## BAD split values for S2
      @sv= split/[:,]/,$path;  $scv= ($sc > 2) ? $sv[4] : $sv[3]; $svc=~s/%//;
      if( $src eq $S2) { $cov=$svc if($svc>0); $mapq=$cov; $path=~s/,Split.*//;  # fix bad
      } else { $mapq = $cov - $scv; }
    } else {  
      unless($path eq "0" or $path eq "1") { next if($mapq{$td}{$src}); }
      $mapq= $cov; ## * $pid/100; ## pid == 0 for all kf2a, missing data
    }
    $mapq=100 if($mapq>100); $mapq{$td}{$src}= $mapq;
    $pid=99 if($pid==0);  $cov=100 if($cov>100);
    $mapv{$td}{$src}= join("\t",$cov,$pid,$gloc,$path);
  }
  $osrc=($src eq $S2)?$S3:$S2;
  bests($td,0) if($mapq{$td}{$osrc});
}
sub bests {
  my($td,$last)= @_;
  my $m2= $mapq{$td}{$S2}; my $m3= $mapq{$td}{$S3};
  if(($m2 and $m3) or $last) {
    $cd2=$cdsq{$td}{$S2}||0; $cd3=$cdsq{$td}{$S3}||0;
    $bsrc= ($m2 > $m3 and $cd2 > $cd3) ? $S2 
         : ($m3 > $m2 and $cd3 > $cd2) ? $S3
         : ($cd2 > $cd3 and $m2 > 0.8*$m3) ? $S2 
         : ($cd3 > $cd2 and $m3 > 0.8*$m2) ? $S3 
         : $S3;
    $cdv= $cdsv{$td}{$bsrc}||"misscd"; $cdq= $cdsq{$td}{$bsrc}; $cdv.=",$cdq%";
    $mapv=$mapv{$td}{$bsrc}||"missmapv"; $mq= $mapq{$td}{$bsrc}; 
    print join("\t",$td,"$td.$bsrc","noho",$cdv,$mapv)."\n"; $did{$td}++;
    }
}
END{ for $td (grep{ not $did{$_} } sort keys %tdin) { 
  if(exists $mapq{$td}) { bests($td,1); } 
  else{ $nom=$nomap{$td}||"missmap"; print join("\t",$td,$td,"noho",0,0,0,$nom,0)."\n"; $did{$td}++; }
} }
BEGIN{ $S2="kf2a"; $S3="kf3n";
print join("\t",qw(AQueryID BestRepID noHoDiff  clAAD cov pid GenomeID:gespan path))."\n"; }' \
 kfish2rae5h{t1,ta}_asm2a2n.diffaa.tab kfish2rae5hpub-kfish2n.map.attr  kfish2rae5h_asm2.map.attr \
  > kfish2rae5h_altbest.tab

#...
  my $m2= $mapq{$td}{$S2}; my $m3= $mapq{$td}{$S3};
  if(0 and ($m2 and $m3)) {
    $cd2=$cdsq{$td}{$S2}||0; $cd3=$cdsq{$td}{$S3}||0;
    $bsrc= ($m2 > $m3 and $cd2 > $cd3) ? $S2 
         : ($m3 > $m2 and $cd3 > $cd2) ? $S3
         : ($cd2 > $cd3 and $m2 > 0.8*$m3) ? $S2 
         : ($cd3 > $cd2 and $m3 > 0.8*$m2) ? $S3 
         : $S3;
    $cdv= $cdsv{$td}{$bsrc}||"misscd"; $cdq= $cdsq{$td}{$bsrc}; $cdv.=",$cdq%";
    $mapv=$mapv{$td}{$bsrc}||"missmapv"; $mq= $mapq{$td}{$bsrc}; 
    print join("\t",$td,"$td.$bsrc","noho",$cdv,$mapv)."\n"; $did{$td}++;
    }

#-------
    kfish2rae5h_asmbest_aahogmap.tab2
    AQueryID        BestHoID        HoDiff  clAAD   cov     pid     GenomeID:gespan path
    Funhe2EKm000003t1       Funhe2EKm000003t1.kf3n  same,0.a,platyfish:ENSXMAP00000008562   OK      100.0   99.9    KN805525.1:5067-10366:- 0
    Funhe2EKm000004t1       Funhe2EKm000004t1       top,14.a,mayzebr:XP_004571781.1 SMAL    93.7    97.3    KN805525.1:25802-58390:+        0

      ** BUG in asm2.map.attr Split paths:
      Funhe2EKm000050t4 cov=200     0       Scaffold0:1159818-1160538:.     
          C1:1159318-1159910,+,100%,Split2/2 << wrong syntax, redo asm2.map.attr??
          
    ==> kfish2rae5h_asm2.map.attr <==
    AQueryID        cov     pid     nexon   GenomeID:gespan:geor    path
    Funhe2EKm000003t1       100%    0       50%,2/4 Scaffold0:5067-10366:-  0
    Funhe2EKm000004t1       91%     0       91%,40/44       Scaffold0:25799-58393:+ 0
    Funhe2EKm000004t2       88%     0       83%,40/48       Scaffold0:25660-58049:+ 0
    
    ==> kfish2rae5hpub-kfish2n.map.attr <==
    AQueryID        cov     pid     splice/nexon    GenomeID:gespan:geor    path
    Funhe2EKm000003t1       100.0   99.9    3/3     KN805525.1:5067-10366:- 0
    Funhe2EKm000004t1       93.7    97.3    22/22   KN805525.1:25802-58390:+        0
    Funhe2EKm000004t2       89.8    94.9    21/21   KN805525.1:39579-58049:+        0

  
    ==> kfish2rae5ht1_asm2a2n.diffaa.tab <==
    OK      Funhe2EKm000003t1       100%,214/214,214;       kf3n
    OK      Funhe2EKm000003t1       100%,214/214,214;       kf2a
    DIFF.-134       Funhe2EKm000004t1       5%,60/1103,969; kf3n
    DIFF.-11        Funhe2EKm000004t1       6%,66/1103,1092;        kf2a
    ==> kfish2rae5hta_asm2a2n.diffaa.tab <==
    DIFF.-94        Funhe2EKm000004t2       6%,61/1063,969; kf3n
    DIFF.-114       Funhe2EKm000004t2       5%,57/1063,949; kf2a
    DIFF.-134       Funhe2EKm000004t3       5%,51/931,797;  kf3n
    DIFF.-174       Funhe2EKm000004t3       5%,46/931,757;  kf2a

        
=cut

=item old pl

'BEGIN{ $SANN=$ENV{sann}||"nan"; $SOUT=$ENV{sout}||"none"; $DROPAN="Target|match"; $drop=0; }
if(/^\W/) { $NOTprint=1 if($ing); next; }
@v=split"\t"; my($r,$s,$t,$at)=@v[0,1,2,8]; chomp($at);
if($s =~ /^($SANN)$/) { if($t eq "mRNA"){ 
  ($id)=$at=~/ID=([^;\s]+)/; $at=~s/ID=$id;//; $at=~s/Split=[^;\s]+;//; 
  $an{$id}=$at; } }
elsif($s =~ /^($SOUT)$/) { $ing=1; if($t eq "mRNA") { 
   ($id)=$at=~/ID=(\w+)/; $drop=($id=~m/_G\d+$/)?1:0; $split=($id=~s/_C(\d+)$//)?$1:0;
   $drop=1 if(/^NOPATH/ or $split);
   ($add)= $at =~ m/(aalen[^;]+;protein=[^;\n]+;cdsoff=[^;\s]+)/;
   s/;($DROPAN)=[^;\n]+//g; 
   $an=$an{$id}||""; 
   if($split) { unless($an){ $an=$at; $an=~s/;($DROPAN)=[^;\n]+//g; $add=""; }  $an="Split=$split;$an"; }
   if($an) { s/ID=$id.*$/ID=$id;$an;$add;/; } }
  else { ($pid)=m/Parent=([^;\s]+)/; s/;($DROPAN)=[^;\n]+//g; s/Parent=$pid/Parent=$id/; }
  print unless($drop); }  ' \
| $evigene/scripts/bestgenes_update.pl -act=sort -debug -cadd addgene=1 \
 -vers kf2rae5h -conf pubgenome/evigene_kfish2_gbsubmit.conf \
 -in stdin -out gmapnsub1/kfish2rae5h_kfish2n_fc4ans.gff

=cut
