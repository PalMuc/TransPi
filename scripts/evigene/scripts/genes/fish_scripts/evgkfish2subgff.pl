#!/usr/bin/env perl
# evgkfish2subgff.pl
# Genome annots with cds2aa checked, trasm pubaa translate exceptions

use strict;

use constant VERSION => 2; # == 2015.04.02 : update w/ various new map.gff inputs, 

my $debug=$ENV{debug}||0;
my $TESTSET=$ENV{test}||0;  # subset, drop most alts..
my $NOSPLITSET=$ENV{nosplit}||0;  # subset, drop split genes
my $VERS=$ENV{vers}||"fc6"; # V2 ? // V1 = fc4.. fc5.. fc9t
my $domakealt= $ENV{makealt} || 0; # obsolete V2
my $doputalt= $ENV{putalt} || $ENV{alt} || 0; # add V2; change to gff reader, gid=nnnt1 filter to nnntalt
my $dogapfix= $ENV{gapfix} || 0; 

my $SplitGeneFixID= $ENV{splitid} || 0; # ncbi sub needs; make split gene parts have new ids w/ suffix _C1/2 ..

my $HODIFFMIN= $ENV{hodiff} || 95; # percent of pubaa-hoscore
my $AAabsDIFFMIN= $ENV{aadiff} || -49; # abs-aasize not pct-aa ??
my $AAPercentDIFFMIN= $ENV{aapctdiff}||75; # OLD AADIFFMIN, not used now?
    #^^ pubaa vs cds2aa criteria, attention to AADIFF qual, use pubaa where cds2aa is short/poor match
    # 44000/110164 are >= 80%; 60000/110000 are >= 50%

my $EVIGENE=$ENV{evigene}||"/bio/bio-grid/mb/evigene";


=item version 2015.04.02

  2015.04.02 : update w/ various new map.gff inputs, 
    new output folder gmapnsub2
    all w/ updated genefindcds2 including split-fix, cds-exon end adjusts
    need new best-version tables
      $bestidtab ="gmapn/kfish2rae5h_asmbest_aahogmapfc6.tab" ??
      $bestalttab="gmapn/kfish2rae5h_altbestfc6.tab";
      ^^ sub makeAltBestTab() needs update or make altbestfc6.tab by hand?
    
    S2/S3 not enough? or is S3 for all inputs of ncbi asm3 map?
    needs SplitGeneFixID work
 
  /bio-grid/kfish2/submitf/evgsubmitkfish2.info
 
  @ingff= 
  gmapn/kfish2rae5h_kfish2a_findcds3t1.gff.gz # replace new findcds2
  gmapn/kfish2rae5h_kfish2a_findcds3ta.gff.gz # replace new findcds2
  gmapn/kfish2rae5h_kfish2n_findcds3.gff.gz   # replace w/ new findcds2
  #adds:
  gmapn/findcds4upd/ 
    kfish2rae5h_kfish2a_t1fcds4 replaces kfish2rae5h_kfish2a_findcds3t1; mRNA n=34886 both (- Split=2)
    kfish2rae5h_kfish2a_tafcds4 replaces kfish2rae5h_kfish2a_findcds3ta; mRNA n=75212 both (-Split=2)
    kfish2rae5h_asm3n_fcds4 ? replaces kfish2rae5h_kfish2n_findcds3 ; has alts+t1
  #new gsplign maps  
  gmapn/gsplign15b/kfish2nsplign15n_fcds4h == findcds4upd/kfish2nsplign15n_fcds4h  ; has alts? + t1
  gmapn/gsplign15/kfish2nsplign15h_fcds4  == findcds4upd/kfish2nsplign15h_fcds4 

  gmapn/findcds4upd/fixset/ 
    {kfish2nsplign15h,kfish2nsplign15h,kfish2rae5h_gmap3n}.fcds5.gff.gz
      ^^ ~8000 t1 ids, cdsexon adjusted, pull out better subset than orig map set
     findcds4upd/kfish2daat1fx.bestof7.tab : picks best source set for fixcds?

  ** make new aahogmap blastp scoring of new cds2aa sets?
  tables of best cds2aa map source sets (and aa diff from pubaa/pubgenes set)
  use for new kfish2rae5h_asmbest_aahogmapfc6.tab, kfish2rae5h_altbestfc6.tab, 
    findcds4upd/kfish2cds2aatab4t1.bestof6.tab
    findcds4upd/kfish2cds2aatab4ta.bestof6.tab
    
    AQueryID	diffa	cds2alen	pubalen	cdsoff	puboff	pubclen	cxlen	cov	pid	path	goodmap	poormap
    
    kfish2cds2aatab4t1.bestof6, n=32733:
    goodmap : 26350 gmap2pub5h, 25151 gmap3npub5h, 18180 gmap3nx11tsa, 8922 gspl3n5n,  5976 gspl3n5h
    gooduniq: 3526 gmap3npub5h, 2291 gmap2pub5h, 722 gspl3n5n, 322 gspl3n5h  
    poormap : 7417 gmap3npub5h, 6305 gmap2pub5h, 4759 gspl3n5n, 3550 gspl3n5h
      skip gmap3nx11tsa as almost same as gmap3npub5h (same trasm)

    kfish2daat1fx.bestof7.tab : n=32733, includes fixset fixcds: gmap3n5fx,gspl3n5nfx,gspl3n5hfx
      -- use only good-fxonly set
    goodmap: 24636 gmap2pub5h, 23576 gmap3npub5h, 7434 gspl3n5n, 4738 gspl3n5h
              3149 gmap3n5fx, 2732 gspl3n5nfx, 2176 gspl3n5hfx      
    good-fxonly: 1359 gspl3n5nfx, 1234 gmap3n5fx, 1168 gspl3n5hfx ; nloci=2674
    good-notfx : 23082 gmap2pub5h, 21629 gmap3npub5h, 6030 gspl3n5n, 3697 gspl3n5h ; nloci=27218
      .. are any of the poor-fx instead good ortho but shorter aa?
      
=cut


my $OUTDIR='gmapnsub3'; # gmapnsub1

my $SPUB='pubaa'; # pub.aa replaces cds2aa
my $S2='kf2a';  # cds2aa on kfish asm2 / Scaffolds
  my $SANN='kf2rae5|kf2rae5alt'; # == kf2a, kfish2rae5h_kfish2a_findcds3t1.gff
  # V2 update SANN: no change?  BAD ANNOT/SANN now from findcds4upd/kfish2rae5ht1asm2_fcds4h.gff
  #  ** need annots from prior gmapn/kfish2rae5h_kfish2a_findcds3t1.gff.gz ..
  #  ** read annots, %an{id} from separate table or gff
  # findcds4upd/kfish2rae5ht1asm2_fcds4h.gff.gz = kf2rae5
  # findcds4upd/kfish2rae5htaasm2_fcds4h.gff.gz = kf2rae5alt
  #... no pub ann??
  # findcds4upd/kfish2rae5h_asm3n_fcds4.gff.gz = kf2rae5h
  # findcds4upd/kfish2nsplign15h_fcds4.gff.gz = splkf3n15h
  # findcds4upd/kfish2nsplign15n_fcds4h.gff.gz = splkf3n15n
  # findcds4upd/fixset/kfish2nsplign15h.fcds5.gff.gz = splkf3n15h
  # findcds4upd/fixset/kfish2nsplign15n.fcds5.gff.gz = splkf3n15n
  # findcds4upd/fixset/kfish2rae5h_gmap3n.fcds5.gff.gz = kf2rae5h

my $S3='kf3n'; # cds2aa on ncbi pubasm KN scafs, need for submit set
  my $SOUT='kf2rae5h'; # == kf3n, kfish2rae5h_kfish2n_findcds3.gff
  # V2 update SOUT: no change?


my $TRAMSID="Funhe2Exx11m";
my $PUBID="Funhe2EK";
my $DROPOID='Funhe2E6bm|Funhe2Emap3m|Funhe2Eq7m|Fungr1EG3m|fungrvel|rfung|^xp|UniRef|Ictpun';  #? estnewbg
## drop most.. say instead what to keep: 
my $KEEPOID='\bFunhe2Exx11m|\bkf4bAUG|\bAUG|PASA'; #? ^Funhe5EG| .. other?

# kf4bAUGpia9cs55g ??
#above/off# my $SplitGeneFixID= 1; # ncbi sub needs; make split gene parts have new ids w/ suffix _C1/2 ..
my $IDSplitSuffix =  "_C";  

    
# DROPAN now used both for public ann input and for working gmap.gff inputs, need slight diff sets 
my $DROPANP='Target|err|aaSize|cdsSize|match|score|oname|groupname|ocds|inshort|utrx|xtrim'; #  pubann.gff
  #?? dropanp: aaSize=432;cdsSize=76%,1296/1705; << duplicate/old/bad? have aalen=,cdsoff= in cds2aa and  pubaa vals
my $DROPANW='match|score|oname|groupname|ocds|inshort|utrx|xtrim'; # mapsrc.gff;  xtrim?
my $DROPAN= $DROPANW;

## V2 update source gff special attr to keep?
my $KEEPAN='Target|trg|nexon|err|gaps|gapfill|gapfix';


## oid=Funhex11,Funhexother << remove 2nd troid
# exon annots: 13651 error from gmap; keep/report in gb.annot? Region=map-exception...
# Parent=Funhe2EKm000208t1;error=ERROR.span:genome_span:1070,tr_span:1241,KN805525.1:5989263-5988193171
# Parent=Funhe2EKm000500t1;error=ERROR.span:genome_span:3687,tr_span:3875,KN805527.1:3629399-3625712188;

# also drop,maybe: ocds=NEn,NE0,NE2,NE3,NE4,NE5;inshort=xcut3:-32;
# also maybe drop extras: oname=othername,  
#------ data ----------------------------------------------

# my $ORIG_bestidtab="gmapn/kfish2rae5h_asmbest_aahogmap.tab2";
# my $OLD2_bestidtab="gmapn/kfish2rae5h_asmbest_aahogmap.tab3";
# #^vv^ make 2nd table for others, alts, etc. best of kf3n, kf2a gmap and/or cds2aa size
# my $OLD2_bestalttab="gmapn/kfish2rae5h_altbest.tab";
#   #v1: 26387 kf2a, 83776 kf3n; note this includes all trs, asmbest_aahogmap as well
#   #v2:  6511 kf2a, 40476 kf3n, 63176 pubaa : use %aadiff to call pubaa .. can do any better w/o hoval?
#   #v3: 20350 kf2a, 58198 kf3n, 31615 pubaa : use only aasize-cds/aasize-pub, better than aadiff
# ## FIXME: revise makealtbest/ altbest.tab to combine asmbest_aahogmap.tab
# ## .. and correct asmbest cases where t1 mismatch to talt and poor map, pubaa/kf3n should be pubaa/kf2a

##  fc6 = 2015.04.02 : update w/ various new map.gff inputs, 
## >> merged now, 1 table gmapn/findcds4upd/kfish2all5_fcds45.besttab3

use constant V2bestidtab => 1;
my $bestidtab ="gmapn/findcds4upd/kfish2all5_fcds45.besttab7"; #  was besttab3 .. 4 .. 5.. 6..
    # .. "gmapn/kfish2rae5h_asmbest_aahogmapfc6.tab";
my $bestalttab=""; # "gmapn/kfish2rae5h_altbestfc6.tab";

my($tsaids,$nxpd,%xpd,%ncbitsa,%augmod,%updates);  
  ## xidtab use this also to make public NCBI ids for submit xref ? prefer bestalt w/ ncbi xref
  # * xidtab updated w/ col3= NCBI TSA pubid
my $xidtab="gmapn/kfish2x11tsa.pubidtab";
my $augidtab="gmapn/kfish2augmod.pubidtab"; # these should come from kf2.gff direct, no cds2aa
my $auggff="gmapn/kfish2rae5h.augmod.gff.gz";
my $uptab="gmapn/kfish2rae5h.update1.tab"; # changes ie. drop/droplocus/replace/... as per genes update

#** THIS IS BAD, drops altnum t1,t2,.. NCBI folk must give ID mapping 
my $TSAIDprefix="GCES01"; # GCES01000000 replace 000000 w/ 009453 of Funhe2Exx11m009453t7 for Funhe2EKm002555t1
  # |From: <gb-admin@ncbi.nlm.nih.gov>
  # |Date: Mon, 2 Feb 2015 13:16:28 -0500
  # |Subject: Accession Fundulus heteroclitus; GCES00000000; PID  PRJNA269174
  # TSA pubid == GCES01000000 + x11 idnum

# add new table of id.src errors: Scaf change for kf2a, split gene missing exons, other ..?
## V2 update this cantuse.ids
  # 5453 err=Missing-mrna-exons : drop or recheck
  # 5700 err=Scaf-change : keep?
my $ErrorIdTab= "gmapn/kfish2rae5h_kfish2a_findcds3.cantuse.ids";

# * FIXME: NEW Annots from pubtsa6 : Name, .. supercede pubgenes.gff annots
# * FIXME2: need pubaa from pubgenes/kfish2rae5h.main.pub.aa.gz to fill in Funhe5EG/AUG set

# ?? from pubtsa6/kfish2rae5x11t6_tsasubmit.aa.gz or .mrna.gz ?
#  also need pubaa prot, ..
# my $pubaain="pubtsa6/kfish2rae5x11t6_tsasubmit.aa.gz";
my @pubaain=qw(
  pubtsa6/kfish2rae5x11t6_tsasubmit.aa.gz
  pubgenes/kfish2rae5h.main.pub.aa.gz
  );

my $scafidtab="pubgenome/ncbifunhe3_kfish2asm.sametab";
my $SUBMITCONF="pubgenome/evigene_kfish2_gbsubmit.conf";
my $GAPGFF="pubgenome/ncbifunhe302scaf.gaps.gff";

# my $onam="kfish2rae5h_kfish2n_";
my $otmpgff="$OUTDIR/kfish2rae5h_asm3n_${VERS}.gff";
my $outsrtgff="$OUTDIR/kfish2rae5h_asm3n_${VERS}ans.gff";
my $outidtab="$OUTDIR/kfish2rae5h_asm3n_${VERS}.idtab";  # all input mRNA ID.src + action for debug
my $outfixgff="$OUTDIR/kfish2rae5h_fx_${VERS}.gff";

my $STAG="none"; my $SMAP="none"; my $HOTAG=$STAG; my $HOTAG2="";
  ($STAG,$SMAP,$HOTAG)= sourceTag($ENV{src}); # set STAG from $ingff name 

## input first: t1asm2 taasm2 == SANN/hasann
## SANN not right on new fcds4, use old gff
## FIXME: pubannotgff use instead pubgenes/kfish2rae5g.{main,alt}.pubsort.gff.gz for annots and for AUG model gff source

my @pubannotgff= qw( 
  gmapn/kfish2rae5h_kfish2a_findcds3t1.gff.gz 
  gmapn/kfish2rae5h_kfish2a_findcds3ta.gff.gz 
);


## 2015.04.22: findcds upd fixes > fcds7 set
## 2015.04.15: findcds upd fix rev xtrim=1 bugggg > fcds6 set
## 2015.04.11: add ingff.cdna.endgaps table addition to these, need to trimNNN end gaps missed by findcds

## 2015.04.20: insert gapfix.gff, replacing others it is drawn from .. annots are output of this script
## gmapnsub3/kfish2rae5h_fc14m_gapfix.gff.gz mRNA n=3121, 
##   mRNA should have either gapfix=xxx n=133 or gapfill=xxx n=2043 annots to replace others, some lack, skip?
##   eg: ID=Funhe2EKm007709t5 insrc=kf2a:gmap2a5u;osrc=gmap2a5u;tblan=trid:gmap2a5u:Funhe2Exx11m005839t5;ggap=nnn/xxx,yyy;gapfill=115589-115592
##  .. exon.gapfill= needs to become misc_feature
## .. maybe can do away with endgaps table work

## 2015.04.25: add asm2 gmap lifted to asm3n ..  kfish2rae5htasm3lft_fcds7 from agplift2gff
##  add in good asm2 gmap lifted to asm3n .. otherwise miss some good mappings
##  $evigene/scripts/agplift2gff.pl -olda ../../pubgenome/kfish4asm2ncbi.agp -newa ../../pubgenome/ncbifunhe302scaf.agp -gff kfish2rae5htasm2dif_fcds7.gff
## FIXME: outofdate fc14m gapfix.gff >> fc15c4 needs auto-update here? or where? lots of errs from outofdate gff
## xxx:     gmapnsub3/kfish2rae5h_fc14m_gapfix.gff.gz
## kfish2rae5h.augmod.gff problems missing aaqual, protein vals; replace w/ kfish2rae5h.augmodfc7, kfish2rae5ht1asm2_fcds7 equiv;  4468 same cds2aa, 232 diff cds2aa

my @ingff= qw(
    gmapn/kfish2rae5h.augmod.gff.gz
    gmapn/findcds4upd/kfish2rae5ht1asm2_fcds7.gff.gz
    gmapn/findcds4upd/kfish2rae5htaasm2_fcds7.gff.gz
    gmapn/findcds4upd/kfish2rae5htasm3lft_fcds7.gff.gz
    gmapn/findcds4upd/kfish2rae5h_asm3n_fcds7.gff.gz
    gmapn/findcds4upd/kfish2nsplign15n_fcds7.gff.gz
    gmapn/findcds4upd/kfish2nsplign15h_fcds7.gff.gz
    gmapn/findcds4upd/fixset/kfish2rae5h_gmap3n_fcds7x.gff.gz
    gmapn/findcds4upd/fixset/kfish2nsplign15h_fcds7x.gff.gz
    gmapn/findcds4upd/fixset/kfish2nsplign15n_fcds7x.gff.gz
    );

# my @OLD2ingff= qw(
#     gmapn/findcds4upd/kfish2rae5ht1asm2_fcds4h.gff.gz
#     gmapn/findcds4upd/kfish2rae5htaasm2_fcds4h.gff.gz
#     gmapn/findcds4upd/kfish2rae5h_asm3n_fcds4.gff.gz
#     gmapn/findcds4upd/kfish2nsplign15n_fcds4h.gff.gz
#     gmapn/findcds4upd/kfish2nsplign15h_fcds4.gff.gz
#     gmapn/findcds4upd/fixset/kfish2rae5h_gmap3n.fcds5.gff.gz
#     gmapn/findcds4upd/fixset/kfish2nsplign15h.fcds5.gff.gz
#     gmapn/findcds4upd/fixset/kfish2nsplign15n.fcds5.gff.gz
#     );
# 
# my @OLDingff= qw(
#   gmapn/kfish2rae5h_kfish2a_findcds3t1.gff.gz 
#   gmapn/kfish2rae5h_kfish2a_findcds3ta.gff.gz 
#   gmapn/kfish2rae5h_kfish2n_findcds3.gff.gz
#   );


my($didhead);
my(%scaf,%scafn); # dont need scafn ?
my(%pubaa,%puban,%didput,%mustid,%seenid,%progsrc,%locustag,%endtrim); # pub ann from 1 gff input, use for next..
my(%mapq,%cdsq,%mapv,%cdsv);  # for makeAltBestTab
my %skipids=(); # ErrorIdTab

#------ methods ----------------------------------------------
sub MAINstub {}

if($domakealt) { 
  if(V2bestidtab) { die "ERR: V2bestidtab, makealt should be done in $bestidtab\n"; }
  makeAltBestTab($bestalttab,undef); #? use skipids below?
  exit; 
}

if($dogapfix) { 
  my $outgff= $ENV{out} || $outfixgff; ## $otmpgff outfixgff
  my @tingff= grep/gff/, @ARGV;
  # rename($outfixgff,"$outfixgff.old") if(-s $outfixgff);
  gapfix($outgff,@tingff); #? use skipids below?
  exit; 
}


readUpdatetab($uptab);  # read early, first? may change other inputs

#** V2bestidtab NEW FORMAT bestidtab
my ($bestidh) = readAsmBestTab($bestidtab); # %bestid{pid}{id,src,hoval,gmap,..}
my ($bestalth)= ($bestalttab) ? readAsmBestTab($bestalttab, $bestidh) : {}; 

checkGffErrors($ErrorIdTab);
readTSAidtab($xidtab);
readAUGMODidtab($augidtab);
readScaffoldTab($scafidtab);  

readPubaa(@pubaain); # # return(\%pubaa,\%puban); ## puban vs %an global, dont need both

# reads to global %an; my($pubannot)= 
readAnnot(@pubannotgff); # reuse $puban{$id}{'all'} ; new V2 .. now have globals: %an, %puban, %$pubannot mix of same

my $outh=undef;
if($debug>1 or $otmpgff =~ /stdout/) { 
  $outh= *STDOUT;
} else {
  warn "#gff.out: $otmpgff\n" if $debug;
  rename($otmpgff,"$otmpgff.old") if(-s $otmpgff);
  open($outh,'>',$otmpgff) or die $otmpgff;
}

my @tingff= grep/gff/, @ARGV;
   @tingff= @ingff unless(@tingff);

my($nint,$nputt,$ndropt)=(0) x 3;
for my $ingff (@tingff) {
  my $flags=""; # per ingff ?
  ($STAG,$SMAP,$HOTAG)= sourceTag($ingff); # set STAG from $ingff name 
  my($nin,$nput,$ndrop)= readGff($ingff,$outh,$flags); #  does putSubGff
  $nint+=$nin; $nputt+=$nput; $ndropt+=$ndrop;
}
warn "#gff.totout: in=$nint,put=$nputt,drop=$ndropt\n" if($debug);
close($outh) unless($otmpgff =~ /stdout/);

warn "#idactions: $outidtab\n" if($debug);
if(open(my $tabh,'>',$outidtab)) {
  # %didput{$id} == nput ,%seenid{$id}{$insrc} == $action ..
  print $tabh "#gff.totout: in=$nint,put=$nputt,drop=$ndropt\n";
  
  ## FIXME outidtab : seenid{id}{progsrc} instead of {insrc}, ncols = ~8 map sources now?
  ## outidtab: drop Src2/Src3 cols, add ActS[12345678..] cols
  my @insrc=($S2,$S3);
  
if(V2bestidtab) { # or VERSION > 1
  my @progsrc= sort keys %progsrc; my $npr=@progsrc;
  my @hd=qw(AQueryID nPut); 
  # for my $i (1..$npr) { push @hd, "actS$i"; } #?? use @progsrc instead of actS$i ?
  for my $i (1..$npr) { push @hd, "S$i".$progsrc[$i-1]; }  
  print $tabh join("\t",@hd)."\n";  
  for my $id (sort keys %seenid) {
    my $np= $didput{$id}||0;
    print $tabh "$id\t$np";
    for my $s (@progsrc){ my $sv= $seenid{$id}{$s} || "noval"; print $tabh "\t$sv"; } 
    print $tabh "\n";
    }
} else {  
  print $tabh join("\t","AQueryID","nPut","Src2","S2act","Src3","S3act")."\n";
  for my $id (sort keys %seenid) {
    my @s= sort keys %{$seenid{$id}}; @s= @insrc if(@s<2);
    my $np= $didput{$id}||0;
    print $tabh "$id\t$np";
    for my $s (@s){ my $sv= $seenid{$id}{$s} || "noval"; print $tabh "\t$s\t$sv"; } 
    print $tabh "\n";
    }
}    
  close($tabh);
}

sortGff($otmpgff, $outsrtgff) if(-s $otmpgff);
#---------------------------------------------

#... subs ........................................................

sub sourcemapTag { my($stag)= @_;  
  return $S3 if($stag =~ m/(gmap3|gspl3)/);
  return $S2 if($stag =~ m/(augmod2|gmap2|gspl2)/);
  return ""; #? $S3; # else should be all others ..
}

sub sourceAnnot {
  my($osrc, $gat)= @_;
  #caller: my($osrc)= $gat =~ m/\bosrc=([^;\s]+)/; # this is only for S2/ $stag="gmap2a5u/h" ?
  if(not $osrc and $gat =~ m/gescore=/ and $STAG=~/gmap2/) { $osrc="gspl2x11"; }
  if($osrc) {
    if($osrc=~/kf2x11gmap/) { $osrc="gmap2x11"; } # == STAG
    elsif($osrc=~/kf2x11gspl/) { $osrc="gspl2x11"; } #  
    elsif($osrc=~/kf4best/) { $osrc="kf2AUG"; } # all these only 
    elsif($osrc=~/kf1gene9f/) { $osrc="kf1AUG"; } # ?? Funhe5EG.asm1 map to asm2
  }
  $osrc ||= $STAG; # add only to tblan=trid:? but what of those w/o tsaid?
}

sub sourceTag {
  my($pt)= @_;
  ## ?? let pt == stag from table, set smap, hotag
  
  ## buggers: kfish2rae5h_gmap3n.fcds5 no stag/smap
  ## FIXME: add tag for $augmod{pubid} set .. these 4,000 are always best t1 set; pt =~ /augmod/ ??
  ## fixme: add gmapnsub3/kfish2rae5h_fc14m_gapfix.gff.gz == mix of all others, out of this script, fixed for gaps
  
  my $stag="none"; my $smap=""; my $hotag="";
  $HOTAG2=""; # global
  
  if($pt=~/fixset/) {
       if($pt=~/splign15h/) { $stag="gspl3n5hfx"; $hotag="gspl5hfx"; }    ## FIXME: fixset.fcds5
    elsif($pt=~/splign15n/) { $stag="gspl3n5nfx"; $hotag="gspl5nfx"; } 
    elsif($pt=~/5h_gmap3n|5h_asm3n/)  { $stag="gmap3hfx";   $hotag="gmap3hfx"; } 
  }  
  elsif($pt=~/augmod/) { $stag="augmod2"; $hotag="c1kf2a"; $HOTAG2="c1kf3n"; }   # UPD hotag == pubaa also
  elsif($pt=~/splign15h/) { $stag="gspl3n5h"; $hotag="gspl5h"; }   
  elsif($pt=~/splign15n/) { $stag="gspl3n5n"; $hotag="gspl5n"; } 
  elsif($pt=~/5h_gmap3n|5h_asm3n/)  { $stag="gmap3n5h"; $hotag="c1kf3n"; $HOTAG2="gmap3h"; } #* HOTAG2=gmap3h ##was gmap3npub5h; FIXME: fixset.fcds5
  elsif($pt=~/5h_kfish2a/){ $stag="gmap2a5h"; $hotag="c1kf2a"; }  # was gmap2pub5h fcds5 update to kfish2rae5ht1asm2_fcds5.*
  elsif($pt=~/5ht.asm2/)  { $stag="gmap2a5u"; $hotag="c1kf2a"; $HOTAG2="gmap2h"; }  #* HOTAG2=gmap2h # fcds5 update to kfish2rae5ht1asm2_fcds5.*
      # ^^findcds4upd/kfish2rae5ht1asm2_fcds4h replaces kfish2rae5h_kfish2a_findcds3t1
      # ^^FIXME 5ht.asm2 includes gsplign as well as gmap methods, but not coded.
  elsif($pt=~/5htasm3lft/)  { $stag="gmap3lft"; $hotag="c1kf2a"; $HOTAG2="gmap2h"; }  #* HOTAG2=gmap2h # fcds5 update to kfish2rae5ht1asm2_fcds5.*
      # ^^ kfish2rae5htasm3lft_fcds7 agplift of 5ht.asm2 to asm3n, special handling?
  elsif($pt=~/x11tsa/) { $stag="gmap3nx11tsa"; $hotag="none";  } 

  if($pt=~/gapfix/) { $stag="gmap3updx"; $smap=$S3; $hotag="c1kf2a"; $HOTAG2="c1kf3n"; } # FIXME: old or new tag?
   #^^ drop this as input src, run gapfix on output.gff before further steps to submit tbl
   
  if($stag =~ m/^(gmap3|gspl3)/) { $smap=$S3; } # else should be all others ..
  elsif($stag =~ /^(augmod2|gmap2|gspl2)/) { $smap=$S2; } 
  # $smap= sourcemapTag($stag); 
  # $STAG= $stag; $HOTAG= $hotag; $SMAP=$smap; # set globals?
  
  return (wantarray) ? ($stag,$smap,$hotag) : $stag;
}


=item readUpdatetab

  updates{id}= act t oid t upd tcomm;  
  various acts now: droplocus,dropalt|drop, partdrop=Split2, changesrc=gspl3n5n, addback
  addback overrides error, other, per mrna id .. need mrna mapsrc also? or use Asmbesttab mapsrc

=cut

sub readUpdatetab {
  my($intab)=@_;  
  warn "#updates.in: $intab\n" if($debug); ## gmapn/kfish2rae5h.update1.tab
  open(F,$intab) or warn "#ERR: Missing $intab\n"; #?? die??
  while(<F>){ 
    next if(/^\W/); s/\#.*//; # tail comments?
    my @v= split" ",$_,5; # hand edited, dont trust \tabs;  
    my($pd,$oid,$act,$upd,$comm)=@v; 
    unless($act =~ /^drop|^partdrop=|^changesrc=|^bestsrc=|^swapbest=|^add|^must/) {
      warn "#ERR updates: unknown action=$act, $pd, $oid, $upd, $comm \n";
    }
    if(my $upold= $updates{$pd}) { # 2+ ok or not? merge actions
      warn "#ERR updates: dup old=$upold, new action=$act, $pd, $oid, $upd, $comm\n";    
    }
    $updates{$pd}=join"\t", $act,$oid,$upd,$comm;
  } close(F); 
}

my $augids=0; # "gmapn/kfish2augmod.pubidtab"; # these should come from kf2.gff direct, no cds2aa
sub readAUGMODidtab {
  my($augidtab)=@_; $augids=0;
  if( open(F,$augidtab) ){ while(<F>){ next if(/^\W/); my($pd,$oid)=split; $augmod{$pd}=$oid; $augids++; } close(F);}
}

=item readTSAidtab

  updated 20150421 w/ NCBI TSA pubids, col3, use those to replace xx11 for tsaid
  ==> gmapn/kfish2x11tsaGCES01.pubidtab <==
Funhe2Exx11m000001t1	Funhe2EKm003842t1	GCES01002296
Funhe2Exx11m000004t1	Funhe2EKm028308t1	GCES01013506
Funhe2Exx11m000004t10	Funhe2EKm028308t10	GCES01058032

=cut

sub readTSAidtab {
  $tsaids=0; # TSA pubid ??? GCES01000000 + x11 idnum Funhe2Exx11m009453t7 ; no altnum needs adding in
  if( open(F,$xidtab) ){ while(<F>){ next if(/^\W/);
    my($xd,$pid,$ncd)=split; 
    $xpd{$pid}=$xd; #  ($ncd=~/\w/)? $ncd : $xd; 
    $ncbitsa{$pid}=$ncd if($ncd);
    $tsaids++; 
  } close(F);
  }
}

sub TSAid { 
  my($pid,$xid)=@_; $xid||=""; $xid=~s/^\w+://;
  my $xidp= $ncbitsa{$pid} || $xpd{$pid} || $xid; #?? both
  #** THIS IS BAD, drops altnum t1,t2,.. NCBI folk must give ID mapping 
  # if($xidp and my($xidnum)= $xidp =~ m/$TRAMSID(\d+)/) {
  #   return $TSAIDprefix . sprintf("%06d",$xidnum); # dont need sprintf ?
  # } 
  return $xidp; #?? or not
}


=item checkGffErrors

add new table of id.src errors: Scaf change for kf2a, split gene missing exons, other ..?
and/or adjust bestid/bestalt tabs for these errors ?

gunzip -c gmapn/kfish2rae5h_kfish2a_findcds3t[1a].gff.gz | grep 'mRNA' | \
 grep  'err=Missing-mrna' | cut -f9 | sed 's/ID=//; s/;.*/	err=Missing-mrna-exons/;' | sort -u \
  > gmapn/kfish2rae5h_kfish2a_findcds3.cantuse.ids

zgrep -c 'err=Missing-mrna' gmapn/kfish2rae5h_kfish2a_findcds3t1.gff.gz  = 2143
zgrep -c 'err=Missing-mrna-exon' gmapn/findcds4upd/kfish2rae5ht1asm2_fcds4h.gff.gz   = 5

=cut

sub checkGffErrors {
  my($erridtab)= @_;
  my($nerr)=(0);
  if(V2bestidtab) {
    # regen Missing-mrna-exons, per gff source, not all S2..
  if( open(F,$erridtab) ) {
    while(<F>) { 
      next if(/^\W/); my($id,$err)=split; $err=~s/err=//; $err=~s/-/_/g;
      if($err =~ /Scaf.change/i) { $skipids{$id}{$S2}= $err; $nerr++;}
      if($err =~ /Missing.mrna.exons/i) { } # skip these.
      #?? BUGS# else { $skipids{$id}{$S3}=$skipids{$id}{$S2}= $err; $nerr++;} # presume skip all src; add new problem cases.
      } close(F);
  }
  } else {
    #above: $ErrorIdTab= "gmapn/kfish2rae5h_kfish2a_findcds3.cantuse.ids";
    if( open(F,$erridtab) ) {
      my $skipsrc=$S2; # if($ErrorIdTab =~ /kfish2a/); # FIXME not this V2
      while(<F>){ next if(/^\W/); my($id,$err)=split; $err=~s/-/_/g; $skipids{$id}{$skipsrc}= $err; $nerr++; } close(F);
    }
  }
  warn "#in.skipids: n=$nerr from $erridtab\n" if $debug;
}


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

** used fa.count files to avoid *SCAFID SCREWUP*, NCBI local2ncbi scaf id table is BAD.
** order of fa scaffolds is same in both, but 1 short Scaffold6427 removed?? index=6427
** NOTE the off-by-1 patch here "samel1": $dn=$sn[$iwn-1]; $dnw=$sn{$dn}; if($dnw==$w) ..
   see Jan 29 12:55 /bio/bio-grid/kfish2/work/kfish2xh.hist
    .. should have used all fa.count cols to test samel1, same nA,nC,nG,nT counts ..
   grep samel1 ncbifunhe3_kfish2asm.sametab     
   6429	samel1	JXMV01064038.1	Scaffold6429	1043	1043	0	0

perl -ne '@v=split; ($d,$w,$g)=@v[0,1,6]; if(/^Scaf/) {  
push @sw,$d; $sw{$d}=$w; $snd=$sn[$iwn]; $snw=$sn{$snd}; if($snw == $w) { $cla="same"; } 
else { $cla="diff"; $dn=$sn[$iwn-1]; $dnw=$sn{$dn}; if($dnw==$w) { $iwn--; $snd=$dn; 
$snw=$dnw; $cla="samel1"; } } if(1) { print join("\t",$iw,$cla,$snd,$d,$w,$snw,$g,$snw-$w)."\n"; } 
$iwn++; $iw++; } else { push @sn,$d; $sn{$d}=$w; $in++; } BEGIN{$iw=$in=0; } ' \
 pubgenome/ncbifunhe302scaf.fa.count pubgenome/killifish20130322asm.fa.count \
 > pubgenome/ncbifunhe3_kfish2asm.sametab

--
head -6430 killifish20130322asm.fa.count | tail
Scaffold6425	1048	319	190	189	350	0	17
Scaffold6426	1047	287	194	193	373	0	16
Scaffold6427	1043	292	134	269	348	0	23
Scaffold6428	1043	271	212	209	351	0	16
Scaffold6429	1043	340	216	175	312	0	13

head -6430 ncbifunhe302scaf.fa.count | tail
JXMV01064035.1	1048	319	190	189	350	0	17  = Scaffold6425
JXMV01064036.1	1047	287	194	193	373	0	16  = Scaffold6426
  >> missing Scaffold6427
JXMV01064037.1	1043	271	212	209	351	0	16  = Scaffold6428, local2ncbi.tab renames Scaffold6427
JXMV01064038.1	1043	340	216	175	312	0	13  = Scaffold6429, local2ncbi.tab renames Scaffold6428

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

=item NEW format bestidtab

gmapn/findcds4upd/kfish2all5_fcds45.besttab3 
AGeneID	AQueryID	BestMap	BestAa	hodiff	hoval	aadiff	cds2alen	cdsoff	cov	pid	maploc	path	othersrc
Funhe2EKm000003t1	Funhe2EKm000003t1	gmap3n5h	gmap3n5h	100	89%,358,platyfish:ENSXMAP00000008562	0	214,94%,complete	19-663	100	100	KN805525.1:5067-10366:-	0	gmap3n5h,gmap2a5u
Funhe2EKm000004t1	Funhe2EKm000004t1	gmap2a5u	gmap2a5u	99	99%,1882,azmolly:XP_007548354	-11	1092,67%,complete	237-3515	91	0	Scaffold0:25799-58393:+	0	gmap2a5u,gspl3n5n,gmap3n5h
Funhe2EKm000004t2	Funhe2EKm000004t2	gmap3n5h	pubaa	0	0	-51	969,67%,complete357-3266	90	95	KN805525.1:39579-58049:+	0	gmap3n5h

=cut


sub readAsmBestTab {
  my( $intab, $betterh )=@_;
  # do also w/ altbest.tab; same format? use %bestid or 2nd hash? BUT skip bestid entries..
  $betterh ||= undef;
  my %bestid= (); my($ndup,$nerr,$nok)=(0) x 9;
  my %idbest=(); # local only?
  #o# warn "#asmbest.in: $intab\n" if($debug);
  open(F,$intab) or die $intab; 
  while(<F>) {
    next if(/^\W/ or /^AQuery|^AGene/);
    chomp; my @v=split"\t";
    my($pubid, $vidsr, $hoval, $clad, $cov, $pid, $genloc, $npath);
    my($vid,$progsrc,$mapsrc,$maperr,$endgaps); ## progsrc == BestAa = aasrc vs mapsrc
    my($bestaa,$hodiff,$aadiff, $cdsaalen,$cdsoff,$swapv);
    
    if(V2bestidtab) {
      #upd:  new maploc flag: ,egap:gapb=106/7 .. egap:gape=5332/5431 >> TRIM UTR exons to these target end points
      #BUG: was missing cdsoff,
      #AGeneID	AQueryID	BestMap	BestAa	hodiff	hoval	aadiff	cds2alen	cdsoff	cov	pid	maploc	path	othersrc
      ($pubid, $vid, $progsrc, $bestaa, $hodiff, $hoval, $aadiff, $cdsaalen, $cdsoff, $cov, $pid, $genloc, $npath)= @v;
      
      ## ,egap: dont have both that and ,err ?? follows ,err if there; syntax messy
      if($genloc=~m/,(egap\S+)/i) { $endgaps=$1; } #??  
      if($genloc=~m/,(err\S+)/i) { $maperr=$1; } #?? hack add maperr flag here? or end of maploc?
      
      $mapsrc= sourcemapTag($progsrc);
      $mapsrc ||= ($genloc =~ /Scaffold/)? $S2: $S3; # no longer valid test, besttab converted to KNscafs
      if($npath =~ s/,(swap\S+)//) { $swapv=$1; } else { $swapv=0; } # hack place for swap:id2/id1
      # 0,swap1:Funhe2EKm000251t1/Funhe2EKm000251t3; Split=1,swap1:Funhe2EKm000412t1/Funhe2EKm000412t6
    } else {
      ($pubid, $vidsr, $hoval, $clad, $cov, $pid, $genloc, $npath)=@v;
      ($vid,$progsrc)= split /\./, $vidsr,2;
      #x $progsrc ||= $SPUB; #"pubaa";
      $bestaa= $progsrc || $SPUB;
      $mapsrc= ($progsrc =~ /($S2|$S3)/) ? $1 : ($genloc =~ /Scaffold/)? $S2: $S3;  
    }

    ## updates{id}= act t oid t upd tcomm; acts: droplocus,dropalt|drop, partdrop=Split2, changesrc=gspl3n5n, addback
    ## .. addback overrides error, other, per mrna id .. need mrna mapsrc also? or use Asmbesttab mapsrc
    ## .. does changesrc=newmapsrc affect readAsmBestTab ? maybe yes
    my $upd=$updates{$pubid}||"";
    if($upd =~ m/^changesrc=(\w+)/) {
      my $newprogsrc=$1; my $oprogsrc= $progsrc; $progsrc=$newprogsrc;
      $mapsrc= sourcemapTag($progsrc);
      ## other things may change w/ newprogsrc, ie genloc, but dont have that data handy
      ## reset endgaps, maperr?
      # dont need upd copy, unless oprogsrc: $bestid{$pubid}{'update'}= $upd;
    }

    $progsrc{$progsrc}++; # for idtab
 
    #?? not here?? or flag in table? 
    ## ** revise makebest table to exclude the skipids as best sources..
    my $err=""; 
    if($err= $skipids{$vid}{$mapsrc}) { $nerr++;} # NOT  next; mapsrc/src here?
    if($maperr) { $err.="," if($err); $err.=$maperr; $nerr++; } # should skip?? flag it?
       
    if(ref $betterh) { # not V2..
      if(exists $betterh->{$pubid} or exists $betterh->{$vid}) { $ndup++; next; }
    }  
 
    # readAsmBestTab problems swapmain not noted, but have rows for both main>altswap and alt=alt
    #? need backmap id hash: idbest{$vid}=$pubid; and dont replace?
    my $vid1= $vid; $vid1=~s/_C\d+$//;
    #? if($vid1 ne $pubid) { } # swapid flag? .. swapv now done in input table
    
    if(my $publast= $idbest{$vid}) {
      $err.=";dup.$vid:$publast/$pubid"; $ndup++;
    }     
    #   my $vid0= $bestid{$publast}{id};
    #   if($vid0 eq $vid) {
    #     #? $vid= $publast; #.. do swap here? shouldnt.. ie turn displaced pubid into this real vid
    #     # but other values are wrong: hoval, bestaa, map cov,..
    #     # no, publast keeps this pubid: $pubid= $vid0;   
    #   }
    #   warn "#asmbest.swapid: $publast/$vid0 <> $pubid/$vid\n" if($debug); # n=621 seen
    #   next if($publast =~ /t1$/); # only for main t1 set?

    $bestid{$pubid}{'id'}=$vid;  $idbest{$vid}= $pubid;
        ## redefine: src == aasrc : S2, S3, SPUB=pubaa; use only mapsrc to match input.gff insrc
    $bestid{$pubid}{'aasrc'}= $bestaa; # == bestaa now, was src
    $bestid{$pubid}{'mapsrc'}= $mapsrc;
    $bestid{$pubid}{'progsrc'}= $progsrc;
    $bestid{$pubid}{'err'}= $err;
    $bestid{$pubid}{'hoval'}= $hoval||"noho"; # clad ??
    $bestid{$pubid}{'gmap'}= "$cov,$pid,$genloc,$npath";
    $bestid{$pubid}{'swapmain'}= $swapv;
    
    ## link alts? 
    my $locid=$pubid; $locid=~s/t\d+$//; $bestid{$locid}{'alts'} .= "$pubid,";
    
    ## note genloc may include ",flags:vals" .. ",err:xxx" or ",egap:xxx" .. move to end after npath?
    
  } close(F);
  warn "#asmbest.in: ndup=$ndup, nerr=$nerr, $intab\n" if($debug);
  return (\%bestid); #?
}

sub addPubAnnot {
  my($gat,$dropannpatt)= @_;
  $dropannpatt ||= $DROPANP; # was global DROPAN
  my($id)= $gat=~/\bID=([^;\s]+)/;
  my $origid = $id; #? or after split/o fixes
  my $issplit=0; #  $drop= $putout= $putpubaa= $swapid= 0;
  return(0,$origid,0,0) if($id=~m/_G\d+$/); # next or not?
  if($id =~ s/_C(\d+)$//){ $issplit=$1; } # kf3n only; missing this? ** NEED this HERE
  my $anat= $gat;
  $anat=~s/ID=$origid;//; 
  $anat=~s/Split=[^;\s]+;//;  $anat=~s/,Split\d[\d\/]+//; ## mapCover=65%,Split2/2; << drop Split append
  ## $anat=~s/Target=/trg=/; ## rename Target= to trg= instead of drop?
  ## err=Missing-mrna-exon # drop all err= in anat ?? added err to $DROPAN
  $anat=~s/;($dropannpatt)=[^;\n]+//g; # 'Target|match|score|oname|groupname|ocds|inshort|utrx|xtrim'  
  
  ## update quality= : drop Split, other things that may change..
  if($anat =~ /(quality=[^;\s]+)/) { 
    my $q=$1; my $r=$q; $r=~s/\-Split//; 
    ## clean up qual:Protein 'unavailable', None
    my($p)= $r=~m/Protein:([^;\s,]+)/; my $pr=$p; 
    if($pr=~s/None|unavailable/partial/i) { $r=~s/Protein:$p/Protein:$pr/; }
    $anat=~s/$q/$r/; }
    
  $puban{$id}{'all'}= $anat  unless($puban{$id}{'all'}); # replaces global $an{$id}= $anat ..
  return($id,$origid,$issplit,$anat);
}

sub readAnnot { ## use %puban, dont replace pubaa annots tho
  my(@ingff)=@_;
  # from orig pub.gff or other?
  my ($nin); 
  # my(%an); # << local or use global %puban as for readGff ??
  for my $inaa (@ingff) {
    warn "#annot.in: $inaa\n" if($debug);
    if($inaa=~/.gz$/) { open(F,"gunzip -c $inaa|") or die $inaa; }
    else { open(F,$inaa) or die $inaa; }
    while(<F>) {
      if(/^\W/) { next; }
      my @v=split"\t"; my($gr,$gs,$gt,$gat)=@v[0,1,2,8]; chomp($gat);
      if($gt eq "mRNA") {  $nin++;
        my($id,$origid,$issplit,$anat)= addPubAnnot($gat);
        } 
    } close(F);
  }
  # return(\%an);
}

# >Funhe2Exx11m011522t4 aalen=418,73%,complete; clen=1711; strand=+; offs=7-1263;
#   pubid=Funhe2EKm003823t1; oid=Fungr1EG3m006328t3; organism=Fundulus_heteroclitus; isoform=1;
#   Name=Arrestin domain-containing protein 3; 
#   Dbxref=CDD:215866,TrEMBL:UniRef50_Q96B67,TrEMBL:ARRD3_HUMAN,family:FISH11G_G1836; 
#   uvcut=Moderate1,17,1712-1728,end3;

# readPubaa: my(%pubaa,%puban);

sub putPubaa {
  my($id,$aa)=@_;
  return unless($id and $aa);
  if(index($aa,'X')>=0) {
    my $oaa=$aa; my $astop= ($aa=~s/\*$//)?1:0;
    my $olen= length($aa);
    my ($cutb)= ($aa=~s/^(X+)//i)?length($1):0;
    my ($cute)= ($aa=~s/(X+)$//i)?length($1):0;
    my $newlen= length($aa); ## add back astop or not?
    if($newlen == $olen) { $aa=$oaa; }
    else {
      my $ncut=$olen-$newlen; 
      warn "#pubaa.trimXXX: $id cut=$ncut, $cutb-$cute\n" if($debug>1);
      if(my $oaalen= $puban{$id}{'aalen'}) {
        $oaalen=~s/^(\d+)/$newlen/; $puban{$id}{'aalen'}= $oaalen;
      }
    }
  }
  $pubaa{$id}=$aa;
}

sub readPubaa {
  my(@inaa)=@_;
  # fixme: add inaa2: pubgenes/kfish2rae5h.main.pub.aa.gz, dont replace inaa id
  ## pubaa puban: aalen offs Dbxref Name Selcstop << treat as most valid
  my $inaaNotX11=0;
  for my $inaa (@inaa) {
    warn "#pubaa.in: $inaa\n" if($debug);
    if($inaa=~/.gz$/) { open(F,"gunzip -c $inaa|") or die $inaa; }
    else { open(F,$inaa) or die $inaa; }
    my($id,$aa);
    ## FIXME: some of these, from Funhe5EG gmap asm2 have end gaps XXX
    while(<F>) {
      if(/^>(\S+)/) { 
        my $xid=$1; # x11 tsa id here
        putPubaa($id,$aa) if($aa); # $pubaa{$id}=$aa if($aa and $id); 
        $aa=""; ($id)= (m/pubid=(\w+)/)?$1:0; 
        if(not $id and $xid =~ /^$PUBID/) { #? $inaaNotX11
          $id=$xid; $id=0 if($pubaa{$id}); # dont replace
          }
        if($id) { # update Name, aalen, Dbxref, clen/strand/offs? oid?
          ## FIXME: Selc: Funhe2EKm035762t1 aalen=238,28%,complete-utrbad,selcstop; Selcstop=1331,1118,962;
          for my $k (qw(aalen offs Dbxref Name Selcstop)) {
            my($v)=m/\b$k=([^;\n]+)/; 
            $puban{$id}{$k}=$v if($v); #? or as string ";$k=$v" ??
          }
        }
      } elsif($id and /\w/) { chomp; $aa.=$_; }
    } close(F);
    
    putPubaa($id,$aa) if($aa); # $pubaa{$id}=$aa if($aa and $id); 
    }
  
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
  
  ($STAG,$SMAP,$HOTAG)= sourceTag($ingff); # set STAG from $ingff name 
  my $UPDATEGFF= ($ingff =~ /gapfix.gff/)?1:0; #?? what?? new STAG?
  
  warn "#gff.in: STAG=$STAG,SMAP=$SMAP $ingff\n" if($debug);
  if($ingff=~/.gz$/) { open(F,"gunzip -c $ingff|") or die $ingff; }
  else { open(F,$ingff) or die $ingff; }
  
  my($drop,$droprow,$issplit,$id,$gid,$pid,$origid,$swapid,$swapfromid,
    $altnum,$anchange,$nput,$ndrop,$nin, $mustkeep,
    $hasan,$insrc,$inprog,$putout,$putpubaa,$npubaa, $endgaps, $XXputout)= (0) x 29;
  my( %didsrc, %splitref);
  %endtrim=();

  ## 2015.04.20: insert gapfix.gff, replacing others it is drawn from .. annots are output of this script
  ## gmapnsub3/kfish2rae5h_fc14m_gapfix.gff.gz mRNA n=3121, 
  ##   mRNA should have either gapfix=xxx n=133 or gapfill=xxx n=2043 annots to replace others, some lack, skip?
  ##   eg: ID=Funhe2EKm007709t5 insrc=kf2a:gmap2a5u;osrc=gmap2a5u;tblan=trid:gmap2a5u:Funhe2Exx11m005839t5;ggap=nnn/xxx,yyy;gapfill=115589-115592
  ##  .. exon.gapfill= needs to become misc_feature
  ## .. maybe can do away with endgaps table work

  
  # $insrc= ...
  while(<F>) {
    if(/^\W/) { next; }
    
    #gff rows
    my $inrow= $_; # changes only to attr? or src?
    my @v=split"\t"; 
    my($gr,$gs,$gt,$gb,$ge,$go,$gat)=@v[0,1,2,3,4,6,8]; chomp($gat);
    $anchange="";
    $droprow=0; ## drop == all of gene/mrna record, droprow == just this inrow
    
    my $SCAFOK=1; my($rnew,$rold);
    if($rnew= $scaf{$gr}{newna}) { 
      $SCAFOK= $scaf{$gr}{same}; #  # problem .. try to shift? need agpLift !
      $rold= $gr; $gr= $rnew;  s/\b$rold\b/$rnew/g; 
      $anchange .= ";scold$SCAFOK=$rold"; 
    }
    
    # LOTS of problems w/ Split genes..
    
    # problems w/ input "gene" rows .. drop and remake next step ??
    if($gt eq "gene") { next; }
    elsif($gt eq "mRNA") { 
      ($id)= $gat=~/\bID=([^;\s]+)/; $nin++;
      $origid = $id; #? or after split/o fixes
      #below# $gid= $id; $gid =~ s/_C\d+$//; $gid=~ s/t\d+$/t1/;
      
      $issplit= $mustkeep= $drop= $putout= $putpubaa= $endgaps= $swapid= $swapfromid= 0;
      $drop++ if($id=~m/_G\d+$/); # next or not?
      $drop++ if($gr =~ m/^NOPATH/); # CHANGE, keep splits;  or $issplit
      
      # check both ID_C and at:Split= annots?
      if($id =~ s/_C(\d+)$//){ $issplit=$1; } # kf3n only; missing this? ** NEED this HERE
      if($gat =~ m/\bSplit=(\w+)/){ $issplit=$1; } # kf2a uses this..
        ## dont trust mapCover annot: m/mapCover=[\d\%]+,Split(\w+)/;  new kf3n gmap may not be split
      ## FIXME: missing most _C2 Split=2  mRNA for _C1 
      
      ## FIXME: doalts not ok here, gid is only t1 now ...
      ($altnum)= $id =~ m/t(\d+)$/; $altnum ||= 1; # or 0?
      ($gid= $id) =~ s/t\d+$/t1/;
 
      ## updates{id}= act t oid t upd tcomm; acts: droplocus,dropalt|drop, partdrop=Split2, changesrc=gspl3n5n, addback
      ## .. addback overrides error, other, per mrna id .. need mrna mapsrc also? or use Asmbesttab mapsrc
      
      
      my $upd= $updates{$id} || $updates{$gid} || "";
      if($upd =~ m/^drop/) {  #  droplocus implies drop all gid 
        $drop++; 
      
      } elsif($upd =~ m/^(must|addback)/) {  
        # mustkeep|addback .. affect what? turn off drop, errs, force putout
        # where does progsrc fit? this implies only pubid must be kept, which version?
        $mustkeep=1; $mustid{$id}=$mustkeep;
      
      } elsif($upd =~ m/^partdrop=(\w+)/) {
        my $pdrop=$1; my($ip)= $pdrop=~m/(\d+)/; 
        ##?? should also remove Split status of other part (if only 1 other)
        if($ip>0 and $ip == $issplit) {
          $drop++; #?? 
        }
      }
      
      if($UPDATEGFF) { $mustkeep=2; $mustid{$id}=$mustkeep; } # also check gapfix= gapfill= annots??
         
      my $dropalt= ($doputalt and $altnum<2) ? 1 : (not $doputalt and $altnum>1) ? 2: 0;
      if($doputalt) {
        if($altnum>1) { $gid=$id; } else { $gid="SKIPt1";  } # what?? $gid="SKIPt1";   # skip t1, but not $drop++; 
        ## FIXME: swapmain changes here? see below cant do here, below at put?
        ## swapmain=Funhe2EKm015804t4,Funhe2EKm015804t1; == alt>main, rev for main>alt
      }
      
      # $act= $intab{$id}{action}; # keep/drop/swap-main-alt/other..
      ## V2 change this hasan,insrc ??
      $hasan= ($gs =~ /^($SANN)$/) ?1:0;
      $insrc= $SMAP; # was ($gs =~ /^($SOUT)$/)? $S3 : ($gs =~ /^($SANN)$/) ? $S2 : "dunno";
      $inprog= $STAG; # ?? any problem w/ STAG == progsrc ? BUGS...
      #^^   ($STAG,$SMAP,$HOTAG)= sourceTag($ingff); # set STAG from $ingff name 
      
      # if(my $augoid=$augmod{$ids}) { } # require orig map-source.gff not findcds set
      my $isaug= ($STAG =~ /augmod/)?1:0; ##? or $augmod{$ids}

      # FIXME: seenid{id}{inprog} instead of insrc !
      $seenid{$id}{$inprog}="$origid;"; # add actions; allow for dup ids here = Splits ?
     
      ## use both %bestid w/ hovals and %bestalt nohoval .. one hash?
      #** ignore bestaltr when have bestgidr{id}, but not when id ne gid
      
      my $Selcstop= $puban{$id}{'Selcstop'} ||""; #? is this right place
      $putpubaa=1 if($Selcstop);

      ## V2 putout FIXME: $STAG of bestidh needs to match ingff $STAG
      ## insrc not enough now;  mapsrc > progsrc ??
      ## THIS IS BAD: $inprog eq $bestgidr->{progsrc} .. whats wrong?
         
      ## FIXME: doalts not ok here, gid is only t1 now ...
      ## doalts: use $bestidh->{$id}||{$origid}  OR change gid to altid
      ## FIXME: new flag on bgloc/gmap: ",egap:gapb=885/0" ",egap:gape=5332/5431"
      ## MUST trim UTR-exons to those ends if flagged
      ## END GAPS: ** need to trim exon start/end if have $endgaps flag from above; requires target>genome  **
      
      my @BESTK= qw(id mapsrc progsrc aasrc err gmap);
      my($bgid,$bgmaps,$bgprgs,$bgaas,$bgerr,$bgloc)= (0) x 9;
      my $bestgidr= $bestidh->{$gid}; # is bestidhash wacky?

      #ok# use constant TESTAGAIN => 1; # bestgidr{@BESTK} is ok
      #ok# if(TESTAGAIN)
      if($gid eq "SKIPt1") { # get bgid to test t1 swapid
        if( my $t1bestr= $bestidh->{$id}) { ## not right, need something else?
          ($bgid,$bgmaps,$bgprgs,$bgaas,$bgerr,$bgloc) = @{$t1bestr}{@BESTK};
          ## bgid here should be swapto: t2,t3 .. bgid is same but _C1,2 split part..
          $bgid=~s/_C\d$//;
          if($bgid ne $id and $bgid ne $origid and $inprog eq $bgprgs and $insrc eq $bgmaps) { 
            $gid= $bgid; # for swapid/swap main to alt
            $bestgidr= $bestidh->{$gid}; # maybe this is it
            # $bestgidr= $t1bestr; #?? wrong? 
            # $putout=1; #?? need this will fail below test isbgid
            # $swapid="$id,$gid"; $dropalt=0;
            # $anchange.=";swapmain=$id,$gid";
          }
        }
      }
      
      if($bestgidr) {
        ($bgid,$bgmaps,$bgprgs,$bgaas,$bgerr,$bgloc) = @{$bestgidr}{@BESTK};
      }
 
      # else { # know this works; dont need above is ok         
      #   if($bestidh->{$gid}) {
      #     ($bgid,$bgmaps,$bgprgs,$bgaas,$bgerr,$bgloc) = @{$bestidh->{$gid}}{@BESTK};
      #   }
      # }
      
      my $checkalt= ($gid eq $id and $bestgidr) ? 0 : 1;
      $checkalt=0 if(VERSION == 2);
      
      if($mustkeep) { 
        $drop=0; 
        if($mustkeep == 2) { #?? if($UPDATEGFF);
          $putout=1; ## input gat annot is supposed to be valid.
          if($putpubaa or $bgaas eq $SPUB or $gat =~ /bestaa=pubaa/) {
            $putpubaa=1; $anchange.=";bestaa=pubaa"; 
          }
          $bestgidr=undef; # dang need more parts below..
          # ignroe bestgidr for mustkeep from gapfix.gff, not other? ..
        } 
      }

      if($bestgidr and not $drop) { #  or $dropalt: defer
        ## FIXME? missing some bgid = id for t2,3 alts here, and split _C1 ..
        ## FIXME let _C2,Cn match _C1 here.
        my $isbgid= ($bgid eq $id or $bgid eq $origid)?1:0;
        $isbgid=1 if($issplit and ($bgid =~ m/^${id}_C/)); # ($bgid eq $id."_C1")
        
        if($isbgid) { #? splits, some id= _C1 ..
          $checkalt=0;
          if($insrc eq $bgmaps and $inprog eq $bgprgs) {
            $putout=1;  
            if($bgerr) { 
              my $flag=($doputalt)?"errskip":"problem";
              $anchange.=";$flag=best:$bgerr"; #?? some of besterr need to be skipped
            }
            
            if($bgloc =~ m/egap:([^;,\s]+)/) { 
              my $egap=$1; # ",egap:gapb=885/0" ",egap:gape=5332/5431"
              if( my($etyp,$eoff,$eend)= $egap =~ m/(.)=(\d+).(\d+)/) {  # "gapb=885/0" "gape=5332/5431"
                ($eoff,$eend)=($eend,$eoff) if($eoff>$eend); # Target coords
                my $trimb= 1 + $eend-$eoff; # ($eoff>$eend) ? $eoff-$eend : $eend-$eoff;
                #?? need adjust eoff,eend to target start, if>1 ??
                if( my($tgd,$tgb,$tge)= $gat =~ m/(?:Target|trg)=(\w+).(\d+).(\d+)/ ) {
                  if($tgb>1) { $tgb-=1; my($eb,$ee)=($eoff+$tgb,$eend+$tgb);
                    $egap=~s,=\d+.\d+,=$eb/$ee,;
                  }
                }
               
                ## need mRNA end-points added to know target ends
                ## maybe need mRNA Target=span to reset gap b/e, ie gapb=0 == mRNA Target start, add start
                my($mgr,$mgb,$mbe,$mgo)=($gr,$gb,$ge,$go);                
                $egap=~s/=/:/;
                $endgaps="trim$etyp:$trimb,mrna:$gb-$ge:$go,$egap"; # gape:5332/5431 gapb:885/0
              }
            }
            
            unless($SCAFOK) { $bestgidr->{mapsrc}=$S3; $XXputout=0; } # dont need this hack now, skipids
            if($putpubaa or $bgaas eq $SPUB) {
              $putpubaa=1; $anchange.=";bestaa=pubaa"; 
            }
          }
       
          #..........
	        # FIXME src eq pubaa with $S3 is bad map ; use besthash->{mapsrc}  instead
          
          if($putout and $id ne $gid) { 
            $swapid="$id,$gid"; $dropalt=0;  #want this swap for seenid for NOT putout 
            $anchange.=";swapmain=$id,$gid";
            ## need output id change .. below
          } 
        }
        # no else  
        # locus pubid not this id, but this may be valid alt .. need alts in bestalt table..
                
        if($debug and 4 > $didsrc{$STAG}++) { 
          warn "#inrow $gid,$id,$SMAP,$STAG, putout=$putout, bgid=$bgid, bgmaps=$bgmaps, bgprgs=$bgprgs, bgaas=$bgaas; anchange=$anchange\n";
        } 
       
      }  
      
      #?? checkalt off now?
      #above# $checkalt=1 if(VERSION == 2); # update for locustag, need all alts in each pass
      if($checkalt and not $drop) {
        my $bestaltr= $bestidh->{$id} || $bestalth->{$id}; 
        if($bestaltr) {
          # my @BESTK= qw(id mapsrc progsrc aasrc err gmap);
          my($bgid,$bgmaps,$bgprgs,$bgaas,$bgerr) = @{$bestaltr}{@BESTK}; # ?? does this fail?
          my $isbgid= ($bgid eq $id or $bgid eq $origid)?1:0;
          $isbgid=1 if($issplit and ($bgid =~ m/^${id}_C/)); # ($bgid eq $id."_C1")
          if($isbgid) {
            if($insrc eq $bgmaps and $inprog eq $bgprgs) {  
              $putout=1;
              unless($SCAFOK) { $bestaltr->{mapsrc}=$S3; $XXputout=0; } # dont need this hack now, skipids
              if(not $dropalt and $putpubaa or $bestaltr->{aasrc} eq $SPUB) {
                $putpubaa=1; $anchange.=";bestaa=pubaa"; 
              }
            }
          }
        } # no else
      }

#      # ** bestgenes_update.pl -cadd addlocus=1 is now doing this.  
#      ## this locustag hash is no good, put mRNA as seen, but need to see all gff before
#      ## deciding below tags. need some premade table of these.
      
      if($mustkeep and not $putout) { # problem ??? YES, problem from besttab not matching..
        warn "#errmust: put=false, $id, origid:$origid\n" if($debug); # dropped?
      }
      
      $putout=0 if($dropalt);
      $anchange.=";insrc=$insrc:$STAG"; # FIXME add $STAG
      $anchange.=";must=$mustkeep" if($mustkeep);

    } else { # gt type NOT mRNA|gene; do nothing?
      ($pid)= $gat=~/\bParent=([^;\s]+)/; 
    }
      
    if($hasan) { 
      # FIXME : transfer some annots of new ingff to puban: gsplign gaps=..MGap; cov=?; trg=?
      # CHANGE: SANN also may be out.gff : kf2a source of intab
      if($gt eq "mRNA") { 
        my($id1,$origid1,$issplit1,$anat1)= addPubAnnot($gat,$DROPANP);
        # my $anat= $gat; # see improved readAnnot()
        # $anat=~s/ID=$origid;//; # do above now
        # $anat=~s/Split=[^;\s]+;//; 
        # $anat=~s/;($DROPAN)=[^;\n]+//g; # Target|match|score|ocds|inshort now .. same for all src?
        # # what of $puban{$id} from readPubaa() : aalen offs Dbxref Name Selcstop
        # unless($puban{$id}{'all'}){ $puban{$id}{'all'}= $anat; }
        # #o# unless($an{$id}) { $an{$id}= $anat; } # what if? presume 1st is best from readAnnot() ??
      } 
    }

    #after putout# $seenid{$id}{$insrc}="$origid,"; # add actions
    
    if($putout) {  ## may be same input gff as SANN
     
      if($gt eq "mRNA") { 
        my $anadd="";
        ## *should always have puban, warn if not..
        my $puban= $puban{$id}{'all'} || "";  # V2 use readAnnot by preference
        warn "#errnoan: $id, origid:$origid\n" if($debug and not $puban); # dropped?
         
        # FIXME: add option debug output : input ids, output ids (action val:keep/drop/swapmain/..)
        #above# $drop++ if($id =~ m/_G\d+$/); # see above
        #above# $drop++ if($gr =~ m/^NOPATH/); # CHANGE, keep splits;  or $issplit
        $drop++ if($TESTSET and $altnum > 3 and not $swapid=~/$id\b/); # and not $swapid=~/$id\b/
        $drop++ if($NOSPLITSET and $issplit); # this cuts out large majority of tbl2asn ERROR 

        if($endgaps and not $drop) {
          my($ntrim,@newg)= trimEndGap($endgaps,$id,$_,$gr,$gt,$gb,$ge,$go,$gat);
          ($_,$gr,$gt,$gb,$ge,$go,$gat)=@newg if($ntrim>0); ## id not in @newg
          ## check if new gb,ge changes cdsoff= .. then cancel/change cds2aa attr; protein= aalen= cdsoff= on mRNA ??
        }

        ## merge errors: errskip= besterr= >> errskip=best:xxx above .. 
        ## if doputalt, drop all/some of errors
        my $serr= $skipids{$id}{$insrc} || "";       
        if($gat =~ /err=(Missing.mrna.exon)/) { my $e=$1; $e=~s/-/_/g; $serr.="," if($serr); $serr.=$e;} # some gat err= are mistakes?
        $anchange.=";errskip=$serr" if($serr);  # should be removed from bestid hashes ??
        
        ## cds2aa protein: keep for many, replace for pubaa-valid
        ## FIXME: splits should always have putpubaa if orig cds offs  spans split points (? unless ho-valid non-split aa best?)
        ## .. this should be set in besttab
      
        unless($putpubaa or $gat=~/protein=\w/) { # fix miss prot, esp Split parts :(
          $putpubaa=1; $anchange.=";bestaa=pubaa"; 
        }
 

	# FIXME: missoff is set of mostly AUG models, or kfish1 FunheEG mapped (some est/pasa/..)
	# these should *not* be putpubaa, but are for vs reasons. many are ERRORs of CDS 
        # 243 / 911 of missoff are ERRORS; may have wrong aa in pub.aa ..
	# grep mRNA kfish2rae5h_asm3n_fc15d3.gff | grep cdsoff=missoff | cut -f9 | \
	#  egrep  'kf[12]AUG|oid=AUG|oid=kf1gene|oid=pasa|oid=kf2rae5:Funhe5EG'  | wc -l = 905;  NOT egrep n= 6

        if($putpubaa and not $puban{$id}{offs}) { # can we cancel and resolve those bugs? create new ones?
	  $putpubaa=0; $anchange=~s/;bestaa=pubaa//g;
	}

        if($putpubaa) { 
          ## fixme, pubaa{} now just aa seq; want pubaa format as for cds2aa?
          ##    aalen=999,xxx,complete;protein=MAZxxxxxVW;cdsoff=111-999;
          #  add qual Protein:curated for pubaa cases, for tbl2asn cdsexcept:"annotated by transcript"
          if(my $aa=$pubaa{$id}) {
            my($aalen)= $puban{$id}{aalen} || "missaalen"; #  =~m/aalen=([^;\s]+)/;
            my($offs)=  $puban{$id}{offs} || "missoff"; # =~m/offs=([^;\s]+)/;  $offs||="missoffs";
            $anadd .= ";aalen=$aalen;protein=$aa;cdsoff=$offs";
          } else {
            $anadd .= ";err=pubaa_missing";
          }
        } else { # cds2aa attr : problem if missing here !??
         my($add1)= $gat =~ m/(aalen[^;]+;protein=[^;\n]+;cdsoff=[^;\s]+)/;
         $add1=~s/\*;cdsoff=/;cdsoff=/; # drop cds2aa stop
         $anadd.=";$add1" if($add1);
        }
        
        ## ?? sub fixAnnot() to handle these various annot updates    
        my($anoid)= $puban =~ m/\boid=([^;\s]+)/;
        my($inoid)= $gat =~ m/\boid=([^;\s]+)/;
        
        #?? need both anoid,inoid ?? skip inoid
        ## oid=Funhe2Exx11m002607t4,kf2rae5alt:Funhe2Exx11m002607t4 << stutter
        # my %oid= map{ $_ =>1 } split",", "$anoid,$inoid";
        my %oid= map{ $_ =>1 } split",", $anoid;
        
        ## add gat:osrc= map sources; 12777 kf2x11gmap, 11243 kf2x11gspl, 4879 kf4best=AUG?, 2057 kf1gene9f=Funhe5EG/AUG, ..
        my($osrc)= $gat =~ m/\bosrc=([^;\s]+)/; # this is only for S2/ $stag="gmap2a5u/h" ?
	      $osrc = sourceAnnot($osrc, $gat);
        $anadd .= ";osrc=$osrc"; # but drop old osrc=?

	      # if(not $osrc and $gat=~m/gescore=/ and $STAG=~/gmap2/) { $osrc="kf2x11gspl"; }
        # if($osrc) {
        #   if($osrc=~/kf2x11gmap/) { $osrc="gmap2x11"; } # == STAG
        #   elsif($osrc=~/kf2x11gspl/) { $osrc="gspl2x11"; } #  
        #   elsif($osrc=~/kf4best/) { $osrc="kf2AUG"; } # all these only 
        #   elsif($osrc=~/kf1gene9f/) { $osrc="kf1AUG"; } # ?? Funhe5EG.asm1 map to asm2
        # }
        # $osrc ||= $STAG; # add only to tblan=trid:? but what of those w/o tsaid?
        
        my @oid= grep /\w/, sort keys %oid;
        my($trasmid)= grep /$TRAMSID/,@oid; # $trasmid=~s/^\w+://; 
        my $tsaid= TSAid($id,$trasmid); # returns TRAMSID if not in tsaid table
        if($tsaid) { $anadd .= ";tblan=trid:$osrc:$tsaid"; } # annot for tbl2asn TSA valid trasm
        ## ^^ add map-prog mid-prefix? gmap,gspl,xxx for asn tablemaker? tblan=trid:gspl3fx:tsaid ??
        
        ## NO GOOD locustag here, need table of all locid splits/scaffold     
        # remove some excess, 2ndary  KEEPOID
        @oid= grep /$KEEPOID/,  @oid; ##  not m/($DROPOID)/ 
        if(@oid) {
          my $oid= join",",@oid; ## skip if missing ??
          unless($puban =~ s/;oid=[^;\s]+/;oid=$oid/) { $puban.=";oid=$oid"; }
        }
        
        ##below# $an.= ";$anchange" if($anchange);
        if($puban{$id}) { # now has {'all'} == pubanmrna above,  to skip
          ## pubaa an subkeys= aalen offs Dbxref Name
          my @pk= grep { $_ ne 'all' } sort keys %{$puban{$id}};
          #already in $anadd# unless($putpubaa) { @pk= grep{ not m/aalen|offs/ } @pk; }  
          for my $pk (@pk) { $puban =~ s/\b$pk=[^;\n]+[;]?//; $puban.= ";$pk=".$puban{$id}{$pk}; } 

          my $Selcstop= $puban{$id}{'Selcstop'} ||""; #? is this right place
          if($Selcstop) { # add /transl_except=(pos:1002..1004,aa:Sec);
            my @ss=split",",$Selcstop;
            for my $sb (@ss) { my $se=$sb+2; $puban .=";transl_except=(pos:$sb..$se,aa:Sec)" if($sb>0); }
            #NOT HERE# $putpubaa=1 if($Selcstop);
           }
        }

        ## FIXME: swapmain changes here? change id t1<>ta ??
        ## $swapid="$id,$gid"; $anchange.=";swapmain=$id,$gid";
        #above# my $swapfromid=0;
        if($swapid and $swapid=~m/$id,(\w+)/) { 
          my $toid=$1; 
          $puban =~ s/ID=$origid/ID=$toid/; # change id below; ID= not in $puban ??
          $swapfromid=$id; $id=$toid;   
        }
        
        if(my $np= $didput{$id}) { # check if mustid has done this one..
          if($np>0 and $mustid{$id} and not $mustkeep) {
            $anchange.= ";mustnot=$np"; $drop++;
          }
        }
        
        if($issplit) { 
          unless($puban){ $puban=$gat; $puban=~s/;($DROPANP)=[^;\n]+//g;  }  #? keep $anadd="";
          $puban="Split=$issplit;$puban"; 
          $splitref{$origid}{$gr}= $issplit; # for CDS/exon missing _C/Split tags !* shouldnt be
          
          #?? add check/drop useless split part here?? as per updates.partdrop, remove problems
          ##  need Target/trg= of split map, and cdsoff=
          
## SplitPartCull bugs, ie drops valid part other alts have; need list of all alt locs? from readAsmBestTab ?
#errpartdrop: Funhe2EKm038575t1 S1, outofcds:84-488/490-2073, insrc=kf3n:gspl3n5h;errpartdrop=split:1
##?? exclude must= from partcull ??
# #errpartdrop: Funhe2EKm029058t1 S1, outofcds:1-386/674-1090, insrc=kf3n:gmap3updx;must=2;errpartdrop=split:1

          my $SplitPartCull= 1; # debug 1st;    
          if($SplitPartCull) {
            my($tgd,$tgb,$tge)= $gat =~ m/(?:Target|trg)=(\w+).(\d+).(\d+)/; $tge||=0; $tgb||=0;
            ## my($add1)= $gat =~ m/(aalen[^;]+;protein=[^;\n]+;cdsoff=[^;\s]+)/;
            my($cdsb,$cdse)= $gat =~ m/;cdsoff=(\d+).(\d+)/; $cdsb||=0; $cdse||=0;
            if($tge>1 and $cdse>1) {
              #??** need check pubaa offs also 
              my $outside=($tgb > $cdse or $tge < $cdsb)?1:0;
              if($outside and $putpubaa) { 
                my($offs)= $puban{$id}{offs}; # =~m/offs=([^;\s]+)/;  $offs||="missoffs";
                ($cdsb,$cdse)= $offs =~ m/(\d+).(\d+)/;
                $outside=($tgb > $cdse or $tge < $cdsb)?1:0;
              }
              
              if($outside) { # check alts at this gr ref; BUG??
                my $locid= $id; $locid=~s/t\d+$//; 
                my @alts= split",", $bestidh->{$locid}{'alts'}; # is bestidhash wacky?
                my $aok=0; 
                foreach my $ad (grep {$_ ne $id} @alts) { $aok++ if($bestidh->{$ad}{'gmap'} =~ m/$gr\b/); }
                $outside=0 if($aok);
              }
              
              if($outside) {
                my $outs="outofcds:$tgb-$tge/$cdsb-$cdse";
                $anchange.= ";errpartdrop=split:$issplit"; $drop++;                 
                warn "#errpartdrop: $id S$issplit, $outs, changes=$anchange\n" if($debug);
              }
            }
          }
          
          
          if($SplitGeneFixID) {  
            # give new ID=oldid_C$issplit 
            # * update attr ;gene=gid  with _C$issplit ..
            my $sid= $id.$IDSplitSuffix.$issplit;
            $origid= $id; # oops? dont change orig origid?
            $puban =~ s/ID=$origid/ID=$sid/; 

            my($mgid)= $puban =~ m/;gene=([^;\s]+)/; # OFF, see locustag work
            # problem: only some alts are split, eg t3 split, t1,t2 not. gid for t1,t2 and t3.s1 should be same
            if($mgid) { 
               if($issplit > 1) { # keep geneid of s=1 unchanged?
                #OFF# my $sg= $mgid.$IDSplitSuffix.$issplit;
                #OFF# $puban =~ s/;gene=$mgid/;gene=$sg/; 
		          }
             }

            $id= $sid; # change here? will affect id hashes below..
            if(my $np= $didput{$id}) { $anchange.= ";errdupid=$np,split=$issplit"; $drop++; }
          }         
          elsif(my $np= $didput{$id} > 1) { $anchange.= ";splitdupid=$issplit/$np";  } # NOT $drop++;
          
        } else {
          if(my $np= $didput{$id}) { $anchange.= ";errdupid=$np"; $drop++; }
        }

        if($putpubaa) {  # qual Protein:curated > validated
          if($puban =~ /quality=([^;\s]+)/) { my $q=$1; my $r=$q; 
            unless($r=~s/Protein:/Protein:validated-/) { $r.=",Protein:validated"; }
            $puban=~s/quality=$q/quality=$r;/;
          } else {
            $anadd.=";quality=Protein:validated";
          }
        }
        
        # dupid bug: $didput{$id} drop dupl, but make sure drop correct one, dupid should not occur..
        # also anchange=~ /errskip=/   drop++ ??

        if($anchange =~ /(errdupid|errskip)=/) { 
          my $errt=$1; ## BUGS were here errskip S3/gmap3 for Scaf-change only in S2/gmap2
          $drop++;  #n=2393 for fc6t          
          warn "#$errt: $id, changes=$anchange\n" if($debug);
        }

        # s/;($DROPAN)=[^;\n]+//g; # dont care, dropping..
        ## maybe $KEEPAN for transfer of new/$gat source ann to final output, e.g Target=, special source atts?
        my @akeep= split /\|/, $KEEPAN;
        for my $ak (@akeep) {
          if(not $anadd =~ m/$ak=/ and $gat =~ /(;$ak=[^;\n]+)/) { $anadd.=$1; }
        }
        
        $puban =~ s/;offs=[^;\s]+// if($anadd=~/cdsoff=/);
        for my $ak (qw(aalen protein cdsoff)) {
          $puban =~ s/;$ak=[^;\s]+// if($anadd=~/$ak=/);
          ## s/;$ak=[^;\s]+// if($anadd=~/$ak=/);
        }
       
        #fix an/puban: genegroup=FISH11G and in Dbxref=..,family:FISH11G : drop one?
        $puban.= ";$anchange" if($anchange);
        $puban.= ";$anadd" if($anadd);  
        
        if($puban) { # always have ?
          $puban=~s/;;+/;/g; $puban=~s/,,+/,/g; $puban=~s/^;//; # clean up an change stutters paralog=FISH11G_G24467,,,,,,,
          unless(s/\bID=$origid.*$/ID=$id;$puban;/) { # SplitGeneFixID fix
            warn "#erranmiss: $id, origid:$origid, an=$puban\n" if($debug);
          }
        } else {
          warn "#erranmiss: $id, origid:$origid\n" if($debug);
          if($origid ne $id and not s/\bID=$origid/ID=$id/) { 
            warn "#erranmiss: $id, origid:$origid, an=$puban\n"; 
          }        
        }
       
      } else { ## gt eq exon,CDS,other?
        ($pid)=m/Parent=([^;\s]+)/; 
        my $npid= $id; # == mRNA id, possibly changed.. need to worry about mismatches?
        
        ## BUG for swapfromid/swapid/swapmain and SPLITs issplit only?
        ## got dupl. exons, 2x total, split stutter w/ id swap?
        ## npid == toid, but Parent is pid == fromid ..
        ## Funhe2EKm009845t1/Funhe2EKm009845t3 swap, both mRNA (t1>t3 issplit), but only 2* t1.exons
        
        ## END GAPS: ** need to trim exon start/end if have $endgaps flag from above; requires target>genome  **
        ## FIXME: only first/last exon(s) .. this is clipping ALL of them :((
        ## FIXME CDS trim endgaps, if exon trim intrudes into cds offset ; need to record exon trims?
        ## .. also need to correct protein=XXX for such cases.
        
        if($endgaps and not $drop) {
          my($ntrim,@newg)= trimEndGap($endgaps,$pid,$_,$gr,$gt,$gb,$ge,$go,$gat); # skips CDS now
          if($ntrim>0){ ($_,$gr,$gt,$gb,$ge,$go,$gat)=@newg; } ## pid not in @newg
          # elsif($gt eq 'CDS' and my $xbe=$endtrim{$pid})  # $endtrim{$pid}.="$gr:$gb-$ge,";  do this in trimEndGap()
          # ** if CDS ntrim, should update mRNA protein attr: trim aa or recalc, fix cdsoffset, aasize, ..
          # .. diff handling for bestaa=pubaa
          # unfortunately, mRNA is already output. should clear this mess up in genefindcds2.pl
          
          $droprow++ if($_ =~ /^#/ and not $debug); #?
        }
        
  ## SPLIT ID problem wrong _C1 assign to exons _C2 << CDS has this problem still, need to track gr/ref of split
  # * Bug in evgkfish2subgff.pl:readGff() issplit << parsed only for mRNA, 
  #   input.gff has mRNA(all splits)>exons>cds split parts out of order need to parse Parent=_C[12] and/or Split=
        my $pids= $pid;  my $csplit=0;
        if($pids=~s/_C(\d+)//) { $csplit=$1; }
        if(m/;Split=(\w+)/) { $csplit=$1; }
        if($csplit or $issplit) {
          my $sn=$csplit; 
          if($csplit==0) { $sn= $splitref{$pid}{$gr}||0; } #  and $gt eq 'CDS' ??
          if($sn==0) { $sn=$issplit; } # if 0 ?? CDS here, bad issplit nums due to sort after exons
          
          if($SplitGeneFixID) {
            my $sid= $origid.$IDSplitSuffix.$sn;
            if($sid ne $npid and $npid=~/^$origid/ and $pid=~/^$origid/) { $npid= $sid; }
          } else { 
            #x# if($pids ne $npid and $pid=~/^$origid/) { $npid= $pids; } # THIS IS BUG for swapid at least
            s/Parent=$pid/Parent=$pid;Split=$sn/ unless(/;Split=$sn/); 
          }
        }
        
        s/Target=/trg=/; ## rename Target= to trg= instead of drop?
        s/;($DROPAN)=[^;\n]+//g; 
        unless($pid eq $npid or s/Parent=$pid/Parent=$npid/) { 
          warn "#err: $gt Parent=$pid to $npid fail\n" ; }
      }
      
      if($drop or $droprow) {
        $ndrop++ if($gt eq "mRNA");
      } else {
        print $outh "##gff-version 3\n" unless($didhead++);
        ## change col2=SRC to SRCOUTmain and  SRCOUTalt ?
        
        ##? #SErr also cases of anchange=~ /errskip=|besterr=/  or drop++ ??
        print $outh "\n" if($gt eq "mRNA");
        print $outh "#SErr." unless($SCAFOK); # got 1018 of these
        print $outh $_;  # CHANGE, refer to intab for keep/drop; ?? use $inrow instead of $_ ??
        ## swapfromid note: didput collects newid, wanted w/ above checks for dups
        ## but results in idtab where didput <> putout disagrees
        if($gt eq "mRNA" and $SCAFOK) { $didput{$id}++; $nput++; $npubaa++ if($putpubaa); }
      }
    } # putout full record
    
    if($gt eq "mRNA") { 
      $anchange =~ s/;insrc=$insrc:$STAG//;
      ## swapfromid may preceed/follow swapto new id; update both or not? maybe both
      ## see above didput{newid}
      if($swapfromid) {  # eg: Funhe2EKm000143t1:from > Funhe2EKm000143t2:to ??
        $seenid{$swapfromid}{$inprog} .= "ok:$putout/$drop$anchange,"; 
        $seenid{$id}{$inprog} .= "ok:$putout/$drop$anchange,";  
      } else {
        $seenid{$id}{$inprog} .= "ok:$putout/$drop$anchange,";  # add actions
      }
    }
  } close(F);
  
  warn "#gff.out: in=$nin,put=$nput,drop=$ndrop,pubaa=$npubaa\n" if($debug);
  return($nin,$nput,$ndrop);
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
   
sub trimEndGap {
  my($endgaps,$pid,$grow,$gr,$gt,$gb,$ge,$go,$gat) = @_;
  my $inrow=$grow;
  ## END GAPS: ** need to trim exon start/end if have $endgaps flag from above; requires target>genome  **
  ## buggy... also need to trim mRNA end
  # $endgaps="trim$etyp:$trimb,mrna:$gb-$ge:$go,$egap"; # gape:5332/5431 gapb:885/0

  
  my $trim=0; 
  my $OFFEND = 9; # allow how much inset from target start,end for gap start,end?
        # .. does it depend on CDS offset? want to trim CDS XXX start,end but not inside
    
  return ($trim,$grow,$gr,$gt,$gb,$ge,$go,$gat) unless($endgaps);
  
  if($gt eq 'CDS' and my $xbes=$endtrim{$pid}) { # do this in trimEndGap()
    #?? Split problems? .. for CDS not using Target, need to check split part exon==CDS
    my($oldgb,$oldge)=($gb,$ge); my $ncdscut=0;
    for my $xbe (split",",$xbes) { 
      my($xgr,$xb,$xe)= split/[:-]/,$xbe; 
      if($xgr ne $gr) { next; } # split diff part
      elsif($ge > $xb and $gb < $xb) {  $trim=$xb-$gb; $gb=$xb; $ncdscut++; } 
      elsif($ge > $xe and $gb < $xe) { $trim=$ge-$xe;  $ge=$xe; $ncdscut++; } 
      ## should be only 1 CDS end trim, check/warn if more? $ncdscut
      ## any case of cds extend beyond both ends of exon?
    }
    if($trim>0) {
      $grow =~ s/\t$oldgb\t$oldge/\t$gb\t$ge/; # do better? # and update gff..
      $grow =~ s/$/;xtrim=$endgaps/;
      my($etyp,$etrimb)= $endgaps =~ m/trim(.):(\d+)/;  # changed to gapb:span, from "gapb=885/0" "gape=5332/5431"
      warn "#endgap: trim=$trim.$etyp, $pid $gt gap:$endgaps\n" if($debug);
      }
    
  } elsif($gt eq "exon" or $gt eq "mRNA") {
    my($oldgb,$oldge)=($gb,$ge); my $xwid=1 + $ge - $gb;
    my($tgd,$tgb,$tge)= $gat =~ m/(?:Target|trg)=(\w+).(\d+).(\d+)/; $tge||=0; $tgb||=0;
    ## *** ^^^ PROBLEM NO Target|trg on exons for some in.gff ***

    my($etyp,$etrimb)= $endgaps =~ m/trim(.):(\d+)/;  # changed to gapb:span, from "gapb=885/0" "gape=5332/5431"
    my($mrnab,$mrnae)= $endgaps =~ m/mrna:(\d+).(\d+)/;
    my($gapb,$gape)  = $endgaps =~ m/gap.:(\d+).(\d+)/;
    ($gapb,$gape)=($gape,$gapb) if($gapb>$gape); # Target coords
    
    if($etrimb<1) { 
      # skip, warn?
      
    #  elsif($etrimb >= $xwid) .. skip exon?
    
    } elsif($tge>0) {
      # BUT to use Target b,e need also adjust gapb,e to mRNA target-start !*
      #? trim == etrimb - 1,off by 1 bug? gapb,gape are 0-origin; tgb,tge are 1-origin
      $gapb+=1; $gape+=1; # change to 1-origin
      if($tgb > $gape or $tge < $gapb) {
        $trim=0;  
      } elsif($etyp eq "e" and $gapb<$tge and $gape>=$tge-$OFFEND) {
        $trim= 1 + $tge - $gapb;  #? $etrimb;  
        if($trim >= $xwid) { }  # what?
        elsif($go eq '-') { $gb += $trim; } else { $ge -= $trim; }
      } elsif($etyp eq "b" and $gape>$tgb and $gapb <= $tgb+$OFFEND) { # problem here requiring gapb<=tgb?
        $trim= 1 + $gape - $tgb;  #? $etrimb;  
        if($trim >= $xwid) { } # what?
        elsif($go eq '-') { $ge -= $trim; } else { $gb += $trim; }
      }
      
    } else { # no Target :(
      my($trimgb,$trimge)=(0,0);
      ## this way has problem w/ splits .. need target offsets to be precise
      if($go eq '-') {
        if($etyp eq 'e') { $trimgb=$mrnab; $trimge=$mrnab+$etrimb; }
        if($etyp eq 'b') { $trimgb=$mrnae-$etrimb; $trimge=$mrnae; }
      } else {
        if($etyp eq 'b') { $trimgb=$mrnab; $trimge=$mrnab+$etrimb; }
        if($etyp eq 'e') { $trimgb=$mrnae-$etrimb; $trimge=$mrnae; }
      }
      if($gb > $trimge or $ge < $trimgb) { # not end exon..
        $trim=0;  
      } elsif($etrimb >= $xwid) { # skip exon??
        warn "#endgap-err: gap > exon size, $pid Target:$tgb-$tge gap:$endgaps\n";   
        #? $trim= $etrimb;  
        #? $grow =~ s/$/;xtrim=$endgaps/;
        #? $grow = "#xdrop.$grow" if($gt eq "exon"); # is this answer??
        
      } elsif($etyp eq "e") { # ge > trimgb + ; gb < trimge -
        $trim= $etrimb;  
        if($go eq '-') { $gb += $trim; } else { $ge -= $trim; }
      } elsif($etyp eq "b") {  # gb < $trimge + ; ge > trimgb -
        $trim= $etrimb;  
        if($go eq '-') { $ge -= $trim; } else { $gb += $trim; }
      }
    }
    
    if($trim>0) {
      $grow =~ s/\t$oldgb\t$oldge/\t$gb\t$ge/; # do better? # and update gff..
      $grow =~ s/$/;xtrim=$endgaps/;
      $endtrim{$pid} .= "$gr:$gb-$ge," if($gt eq 'exon'); #?? splits need target tgb,tge? or gr?
      }
    if($debug and $trim>0) { # show work?
      warn "#endgap: trim=$trim.$etyp, $pid $gt Target:$tgb-$tge gap:$endgaps\n";
      warn "#endgap-before: $inrow";
      warn "#endgap-after : $grow";
    }
  }
  return($trim,$grow,$gr,$gt,$gb,$ge,$go,$gat);  ## * pid not returned
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
    if($cdq < $AAPercentDIFFMIN) { $bsrc=$SPUB; }     # AADIFFMIN  
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


sub exonMapErrStub {} # see sub, notes at end

sub gapfixMarkExon {
  my($ingap, $ix, $nx, $canfill,$pubaa, $gff)= @_;
  #from mrna: canfill, pubaa,
  my($rf,$tp,$rb,$re,$ro,$attr)= @$gff[0,2,3,4,6,8];
  unless($ingap) { ($ingap)= ($attr=~m/;ggap=([^;\s]+)/)?$1:0; }

  if($ix<1) { $ix=1; } elsif($nx>0 and $ix>$nx) { $ix=$nx; } #?? shouldnt see these, CDS off?
  my $isend= ($ix == 1 or ($nx>0 and $ix == $nx))?1:0; ## when nx == 0, what?
  my $isrev= ($ro eq "-")?1:0;
  my $isexon=($tp eq 'exon')?1:0;
    
  ## handle $nx == 0 ; dont know end exon till read all exons
 
  if( my($ob,$oe)= $ingap =~ m/ovspan:(\d+).(\d+)/ ) {
    my ($didfix,$fixnote)=(0,0);
    my $isinner=($ob > $rb+2 and $oe < $re-2)?1:0;
    if($isinner) { $isend=0; }
    ## test both ix=1,ix=nx when nx=1?
    elsif($nx==1) { $isend=0 if($isinner); } ## always isend otherwise for 1 exon only, eg split parts
    elsif($isend and $ix==1) { if($isrev) { $isend=0 if($oe<$re-2); } else { $isend=0 if($ob > $rb+2); } }
    elsif($isend and $ix==$nx) { if($isrev) { $isend=0 if($ob > $rb+2); } else { $isend=0 if($oe<$re-2); } }

    ## $isinner=($isinner or ($canfill and not $isend)); #?? nx==0 bug
    
    ## FIXME for CDS ends hit gap but not exon.. need to trim CDS, separate isend?
    
    if($isinner or ($canfill and not $isend)) { 
      if($canfill and $isexon) { $fixnote="gapfill=$ob-$oe"; } 
      else { $fixnote="gapnofix=inner" if($isexon); } 
    } elsif($isend or $pubaa) { 
      my $fix= ($ob <= $rb+2 and $oe <  $re-2)? "btrimto.".($oe+1)
          : ($ob  > $rb+2 and $oe >= $re-2)? "etrimto.".($ob-1)
          : ($ob <= $rb+2 and $oe >= $re-2)? "dropexon": "";
      if($fix) { $fixnote="gapfix=$fix"; } 
      elsif($isend) { $fixnote="gapnofix=end"; }
      }
    if($fixnote) { $attr=~s/$/;$fixnote/; $gff->[8]= $attr; }  
    # only attr fixnote is added for genefixgaps.pl input
    }        
    
  return($attr);
}

#  do 3-level sort: mRNA rec sorted by type; gene>mRNA sorted by alt-tr num?; chr>genes sorted by location
sub _sortgene  
{
  # rec == [$ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr]
  my($ta,$tb)= map{ (/gene/)?1:(/mRNA/)?2:(/exon/)?3:(/CDS/)?4:5; } ($a->[2],$b->[2]); 
  return ($a->[0] cmp $b->[0]) # ref  : splits separate
      || ($a->[3] <=> $b->[3]) # begin
      || ($b->[4] <=> $a->[4]) # end: b<>a wanted; a<>b will put shorter CDS > exon > mRNA,
      #not# ($a->[4] <=> $b->[4]) # end: a<>b will put shorter CDS > exon,
      || ($ta <=> $tb) # type mRNA/exon/CDS, now only when same begin/end
      ;
}
#   $OLDsort= ($ta <=> $tb)
#       || ($a->[0] cmp $b->[0]) # ref : should all be same for gene record
#       || ($a->[3] <=> $b->[3]) # begin
#       || ($a->[4] <=> $b->[4]) # end
#       ;

sub gapfixAllExons {
  my($canfill,$pubaa,$mrnaexons)= @_;
  # sort @mrna so exon>CDS, gapfix paired exon,cds
  my @mrna= sort _sortgene @$mrnaexons;

if(1) { # use this way
  # FIXME: test CDS separate gapfix, cure? 150 ERROR:   SEQ_FEAT.FeatureBeginsOrEn..
  my($nx,@exonc);
  @exonc= grep{ $$_[2] =~ /^exon/ } @mrna;
  $nx= @exonc;
  for(my $i=0; $i < $nx; $i++) {
    my $x= $exonc[$i];
    if( my($ingap)= $$x[8]=~m/;ggap=([^;\s]+)/) {
      gapfixMarkExon($ingap, $i+1, $nx, $canfill, $pubaa, $x); }
  }
  @exonc= grep{ $$_[2] =~ /^CDS/ } @mrna;
  $nx= @exonc;
  for(my $i=0; $i < $nx; $i++) {
    my $x= $exonc[$i];
    if( my($ingap)= $$x[8]=~m/;ggap=([^;\s]+)/) {
      gapfixMarkExon($ingap, $i+1, $nx, $canfill, $pubaa, $x); }
  }
} else {
  my @exonc= grep{ $$_[2] =~ /^(exon|CDS)/ } @mrna;
  my $ix=0;  my $nx= scalar( grep{ $$_[2] =~ /^exon/ } @exonc);
  for my $x (@exonc) {
    $ix++ if($$x[2]=~/^exon/);
    my($ingap)= ($$x[8]=~m/;ggap=([^;\s]+)/)?$1:0;
    gapfixMarkExon($ingap, $ix, $nx, $canfill, $pubaa, $x) if($ingap);
  }
}
} 
 
sub gapfixLastExon {
  my($xingap, $nx, $canfill,$pubaa, $mrna)= @_;
  ## ADD: fixlastexon($lastxingap) if($nx==0 and $lastxingap); # but do both last exon, CDS
  if($nx==0 and $xingap and @$mrna > 1) { 
    my($li,$lc)=(0,0);
    if($mrna->[-1]->[2] =~ /exon/) { $li=-1; $lc=-2 if($mrna->[-2]->[2] =~ /CDS/); }
    elsif($mrna->[-2]->[2] =~ /exon/) { $li=-2; $lc=-1 if($mrna->[-1]->[2] =~ /CDS/);}
    if($li) {
      my $xon= $mrna->[$li]; # lastxn
      my $cxon=0;
      if($lc) {
        $cxon= $mrna->[$lc];  
        $cxon=0 unless($$cxon[3] == $$xon[3] or $$cxon[4] == $$xon[4]); # test before fix
        }
      gapfixMarkExon(0, 999, 999, $canfill, $pubaa, $xon); 
      if(ref $cxon) { gapfixMarkExon(0, 999, 999, $canfill, $pubaa, $cxon); }
    }
  }
}

=item gapfix

env gapfix=1 vers=fc15c4 debug=1 ./evgkfish2subgff.pl gmapnsub3/kfish2rae5h_asm3n_{fc15c4,fc15c4alt}.gff

#gapfix.in: outgapx=gmapnsub3/kfish2rae5h_fx_fc15c4.gapx.gff,outgapfix=gmapnsub3/kfish2rae5h_fx_fc15c4.gapfix.gff,outnogap=gmapnsub3/kfish2rae5h_fx_fc15c4.gapnot.gff, 
  ingff=gmapnsub3/kfish2rae5h_asm3n_fc15c4.gff gmapnsub3/kfish2rae5h_asm3n_fc15c4alt.gff, gaps=pubgenome/ncbifunhe302scaf.gaps.gff

#gapfix.cmd: /bio/bio-grid/mb/evigene/scripts/overlapfilter -in gmapnsub3/kfish2rae5h_asm3n_fc15c4.gff -over pubgenome/ncbifunhe302scaf.gaps.gff -act markbaseid -SPANBASEOVER -mark=ggap -pct=1
#gapfix.cmd: /bio/bio-grid/mb/evigene/scripts/overlapfilter -in gmapnsub3/kfish2rae5h_asm3n_fc15c4alt.gff -over pubgenome/ncbifunhe302scaf.gaps.gff -act markbaseid -SPANBASEOVER -mark=ggap -pct=1

#gapfix.cmd: /bio/bio-grid/mb/evigene/scripts/genes/genefixgaps.pl -noexonfix -input gmapnsub3/kfish2rae5h_fx_fc15c4.gapx.gff > gmapnsub3/kfish2rae5h_fx_fc15c4.gapfix.gff

#gapfix.out: err=0,nfix=2868,nnofix=99959,outgapfix=gmapnsub3/kfish2rae5h_fx_fc15c4.gapfix.gff,outnogap=gmapnsub3/kfish2rae5h_fx_fc15c4.gapnot.gff

=cut

sub gapfix {
  my($outgff,@ingff)= @_;
  my($nfix,$nnofix, $outnogap,$houtnogap, $outgapx, $houtgapx, $outgapfix);
  $outgff="genes.gff" unless($outgff);
  $outgff.=".gff" unless($outgff =~ /\.gff/);
  ($outnogap= $outgff)  =~ s/\.gff/.gapnot.gff/; # bad if not .gff
  ($outgapx= $outgff)   =~ s/\.gff/.gapx.gff/;
  ($outgapfix= $outgff) =~ s/\.gff/.gapfix.gff/;

  warn "#gapfix.in: ingff=@ingff, gaps=$GAPGFF\n" if($debug);
  # warn "#gapfix.outgapx=$outgapx,outgapfix=$outgapfix,outnogap=$outnogap, \n" if($debug);
  die "missing $GAPGFF" unless(-f $GAPGFF);
  open($houtnogap, '>', $outnogap) or die "write $outnogap";
  open($houtgapx,  '>', $outgapx) or die "write $outgapx";
  
  ## write only gap-fixed? but want to replace those from input.gff to output.
  ## write two files: gapxgenes.gff, nogapgenes.gff (most), then 
  #   $evigene/scripts/genes/genefixgaps.pl -noexonfix < gapxgenes.gff > gapfixgenes.gff

  for my $ingff (@ingff) {
    
    my $cmd="$EVIGENE/scripts/overlapfilter -in $ingff -over $GAPGFF -act markbaseid -SPANBASEOVER -mark=ggap -pct=1";
    warn "#gapfix.cmd: $cmd\n" if($debug);
    open(GFPIPE,"$cmd |") or die "ERR: $cmd";  ##  > kfish2rae5h_fc14m_gapall.gff
    
    my($mingap,$xingap,$canfill,$pubaa,$nx,$jx,$ix, @mrna);
    while(<GFPIPE>) {
      if(/^\W/) { 
        print $houtnogap $_; 
      } else { # gff..
        chomp; my @v=split"\t";  
        my($rf,$tp,$rb,$re,$ro,$attr)= @v[0,2,3,4,6,8];
        my($id)= $attr=~m/(?:ID|Parent)=([^;\s]+)/;
        my($ingap)= ($attr=~m/;ggap=([^;\s]+)/)?$1:0;
        
        if($tp eq 'mRNA') {
        
          if(@mrna) {
            ## ADD: fixlastexon($lastxingap) if($nx==0 and $lastxingap); # but do both last exon, CDS
            #old# gapfixLastExon($xingap,$nx, $canfill,$pubaa,\@mrna) if($nx==0 and $xingap);                        
            gapfixAllExons($canfill,$pubaa, \@mrna) if($xingap); 
            if($xingap) { $nfix++; } else { $nnofix++; }
            my $hout= ($xingap)? $houtgapx : $houtnogap; 
            for my $r (@mrna) { print $hout join("\t",@$r),"\n"; }
          }
          
          @mrna=(); $xingap=0; 
          $mingap= $ingap;
          $canfill= ($attr=~m/;insrc=\w+:augmod/ or $attr=~m/;osrc=\w+AUG/ or $attr=~/;oid=AUG/ or $attr=~/;trg=Funhe5EG/)? 0 :1;
          ## nocanfill add: oid=AUGepir8:AUGepir8s338g50t1 ; trg=Funhe5EG009381t1
          $pubaa= (m/bestaa=pubaa/)?1:0; 
          $nx= (m/nexon=(\d+)/)?$1:0; $jx=0;
          if($ingap) { $v[8]=~s/;ggap=[^;\s]+//; } # drop ggap for output annot, not informative on mRNA ..
          
        } elsif($tp =~ /^(exon|CDS)/ and $mingap) {
          my $isexon=0; if($tp =~ /^exon/) { $jx++; $ix= (m/;ix=(\d+)/)?$1:$jx; $isexon=1; }
          
          $xingap++ if($ingap); # use gapfixAllExons() at end instead of below
  
          ## BUG exon gapfix w/o CDS gapfix .. CDS sorted before its exon.
          ## .. need to collect full @mrna, sort then gapfix() ? ie gapfixAllExons($canfill,$pubaa,\@mrna);
          if( 0 and $ingap ) { 
            my($attrfix)= gapfixMarkExon($ingap, $ix, $nx, $canfill, $pubaa, \@v);
            $v[8]= $attrfix if($attrfix);
            $xingap= $ingap; ##not ++;
          }
          
          #  if( 0 and $ingap ) {}
          
        } # no else for other gff type
        push @mrna,\@v;
      }
    }

    if(@mrna) {
      gapfixAllExons($canfill,$pubaa, \@mrna) if($xingap); 
      #old# gapfixLastExon($xingap, $nx, $canfill,$pubaa, \@mrna) if($nx==0 and $xingap);                        
      if($xingap) { $nfix++; } else { $nnofix++; }
      my $hout= ($xingap)? $houtgapx : $houtnogap;
      for my $r (@mrna) { print $hout join("\t",@$r),"\n"; }
      @mrna=(); $xingap=0; 
    }
  } # next ingff ..
  
  close(GFPIPE); 
  close($houtnogap); close($houtgapx);
  
  my $fxcmd="$EVIGENE/scripts/genes/genefixgaps.pl -noexonfix -input $outgapx > $outgapfix";
  warn "#gapfix.cmd: $fxcmd\n" if($debug);
  my $err= system($fxcmd);
  warn "#gapfix.out: err=$err,nfix=$nfix,nnofix=$nnofix,outgapfix=$outgapfix,outnogap=$outnogap\n" if($debug);
  return($err,$nfix,$nnofix,$outgapfix,$outnogap);
} 

=item gapfix work in here??

  /bio-grid/kfish2/submitf/gapgenefix_kfish2.info
  
  set dgenome=ncbifunhe302scaf
  cat $dgenome.fa | env src=$dgenome perl -ne 'if(/^>(\S+)/){ putg() if $d; $d=$1; $s=""; } \
  elsif(/^\w/){ chomp; $s.= uc($_); } END{putg()}\
  sub putg { $src=$ENV{src}||"assembly"; $i=0; $w=length($s); while($i<$w) { $b=index($s,"N",$i); if($b<$i) { $i=$w; } \
  else { $e=$b; $f=index($s,"A",$e); if($f>=$e) { $g=rindex($s,"N",$f); $e=$g if($g>$e); } \
  if($e==$b) { while(substr($s,$e,1) eq "N" and $e<$w){ $e++; } } $i=$e+1; \
  $gw=1+$e-$b; $b++; $e++ if($e<$w); print join("\t",$d,$src,"gap",$b,$e,$gw,"+","."),"\n"; } } }' \
    > $dgenome.gaps.gff
  
  ## old off-by -1 gap b,e
  
  #  new -SPANBASEOVER opt for overgap spans on exons; act markbaseid, gapid=gapKN805525_13142 == ref:startb
  $evigene/scripts/overlapfilter -in kfish2rae5h_fc14m_all.gff.gz -over ../pubgenome/ncbifunhe302scaf.gaps.gff \
    -act markbaseid -SPANBASEOVER -mark=ggap -pct=1 > kfish2rae5h_fc14m_gapall.gff
  
  grep ggap= kfish2rae5h_fc14m_gapall.gff | grep -v mRNA | grep exon | perl -ne \
   '($d)=m/Parent=(\w+)/; print "$d\n";' | sort -u > kfish2rae5h_fc14m_gapx.ids
  
  # * add gap-trim/fix code, trim end exons w/ gaps at ends, but not internal gaps 
  # * genefixgaps add gapfix=btrimto-33293441 ; gapfix=etrimto-33293441; gapfix=drop all gaps; codes for end exons/cds w/ gap ovspan at ends;
  # .. trim all exon-end gaps, even internal exons? shouldnt affect cds, other submit data IF bestaa=pubaa, but can affect cds2aa set
  # .. gapfill= gap span filled by trasm, cancel if AUGmodel, ok insrc/osrc= other, trasm x11
  # .. need to handle list of ovspan:55015-55114,55255-55256 ..
  
  perl -ne 'if(/^(\w+)$/){ $ok{$1}=1 } 
  else { if(/\tmRNA/){ ($id)= m/ID=([^;\s]+)/; $ok=$ok{$id}; 
    $canfill= (m/insrc=\w+:augmod/ or m/osrc=\w+AUG/)? 0 :1;
    $pubaa=(m/bestaa=pubaa/)?1:0; $nx= (m/nexon=(\d+)/)?$1:0; $jx=0; } 
  elsif($ok and /\t(exon|CDS)/) { @v=split"\t"; ($rf,$rb,$re,$ro)=@v[0,3,4,6];
    $isexon=0; if(/\texon/) { $jx++; $ix= (m/ix=(\d+)/)?$1:$jx; $isexon=1; }
    $isend=($ix == 1 or $ix == $nx)?1:0; $isrev=($ro eq "-")?1:0;
    if(($ga)= m/;ggap=([^;\s]+)/) { 
      if( ($ob,$oe)= $ga=~m/ovspan:(\d+).(\d+)/ ) {
        $isinner=($ob > $rb+2 and $oe < $re-2)?1:0;
        if($isinner) { $isend=0; }
        elsif($isend and $ix==1) { if($isrev) { $isend=0 if($oe<$re-2); } else { $isend=0 if($ob > $rb+2); } }
        elsif($isend and $ix==$nx) { if($isrev) { $isend=0 if($ob > $rb+2); } else { $isend=0 if($oe<$re-2); } }
        $isinner=($isinner or ($canfill and not $isend));
        if($isinner) { s/$/;gapfill=$ob-$oe/ if($canfill and $isexon); }
        elsif($isend or $pubaa) { 
          $fix= ($ob <= $rb+2 and $oe <  $re-2)? "btrimto.".($oe+1)
              : ($ob  > $rb+2 and $oe >= $re-2)? "etrimto.".($ob-1)
              : ($ob <= $rb+2 and $oe >= $re-2)? "dropexon": "";
          s/$/;gapfix=$fix/ if($fix); }
        }
    }
  } print if $ok; }' \
   kfish2rae5h_fc14m_gapx.ids kfish2rae5h_fc14m_gapall.gff  \
    > kfish2rae5h_fc14m_gapx.gff
  
  $evigene/scripts/genes/genefixgaps.pl -noexonfix < kfish2rae5h_fc14m_gapx.gff > kfish2rae5h_fc14m_gapfix.gff


=cut

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

=item try4+

   cat gmapn/findcds4upd/kfish2all5_fcds45.besttab3 | grep 't1      Funhe' | cut -f1 | sort -u |  wc -l
    32765 : should have this many maint1 in output.gff : NOT
    #gff.totout: in=283014,put=29607,drop=230 in log.evgsubfc11f
    >> missing ~3,800 ? are these the nomap set? would the be in besttab3 ?
    cat kfish2rae5h_asm3n_fc11f.idtab | grep  't1  1       ' | wc -l =   28945
    cat kfish2rae5h_asm3n_fc11f.idtab | grep  't1  0       ' | wc -l =    5980
      n=5980 put=0 ie not in kfish2rae5h_asm3n_fc11f.gff
    >> likely odd AQueryID link to AGeneID is problem, eg 019t3, 041t1_C1 ..
    eg.
  cat kfish2rae5h_asm3n_fc11f.idtab | grep 't1   0  ' | head
AQueryID       	nPut	Src2	S2act                     	Src3	S3act
Funhe2EKm000019t1	0	kf2a	Funhe2EKm000019t1;ok:0/0;scold1=Scaffold0,gmap2a5u,	kf3n	Funhe2EKm000019t1;ok:0/0,gspl3n5h,
Funhe2EKm000041t1	0	kf2a	Funhe2EKm000041t1;ok:0/0;scold1=Scaffold0,gmap2a5u,	kf3n	Funhe2EKm000041t1_C2;ok:0/0,gspl3n5nfx,
Funhe2EKm000090t1	0	kf2a	Funhe2EKm000090t1;ok:0/0;scold1=Scaffold0,gmap2a5u,	kf3n	Funhe2EKm000090t1;ok:0/0,gspl3n5nfx,
Funhe2EKm000092t1	0	kf2a	Funhe2EKm000092t1;ok:0/0;scold1=Scaffold0,gmap2a5u,	kf3n	Funhe2EKm000092t1_C2;ok:0/0,gspl3n5n,
Funhe2EKm000124t1	0	kf2a	Funhe2EKm000124t1;ok:0/0;scold1=Scaffold0,gmap2a5u,	kf3n	Funhe2EKm000124t1_C2;ok:0/0,gspl3n5h,
Funhe2EKm000137t1	0	kf2a	Funhe2EKm000137t1;ok:0/0;scold1=Scaffold0,gmap2a5u,	kf3n	Funhe2EKm000137t1;ok:0/0,gspl3n5h,
Funhe2EKm000143t1	0	kf2a	Funhe2EKm000143t1;ok:0/0;scold1=Scaffold0,gmap2a5u,	kf3n	Funhe2EKm000143t1;ok:0/0,gspl3n5hfx,
Funhe2EKm000182t1	0	kf2a	Funhe2EKm000182t1;ok:0/0;scold1=Scaffold0,gmap2a5u,	kf3n	Funhe2EKm000182t1;ok:0/0,gspl3n5h,
Funhe2EKm000191t1	0	kf2a	Funhe2EKm000191t1;ok:0/0;scold0=Scaffold95,gmap2a5u,	kf3n	Funhe2EKm000191t1_C2;ok:0/0,gspl3n5h,
Funhe2EKm000203t1	0	kf2a	Funhe2EKm000203t1;ok:0/0;scold1=Scaffold0,gmap2a5u,	kf3n	Funhe2EKm000203t1;ok:0/0,gspl3n5h,
  grep t1 1 :  S2act ok:1/ n=4687; S3act ok:1/ n=15701
---
AGeneID          	AQueryID	        BestMap	  BestAa	  maploc	 
Funhe2EKm000019t1	Funhe2EKm000019t3	gmap3n5h	gmap3n5h	KN805525.1:267309-281534:+
Funhe2EKm000041t1	Funhe2EKm000041t1_C1	gspl3n5h	pubaa	KN805525.1:941170-941979:+
Funhe2EKm000090t1	Funhe2EKm000090t1	gmap2a5u	pubaa	Scaffold0:2282686-2286860:.
Funhe2EKm000092t1	Funhe2EKm000092t1	gmap2a5u	gmap2a5u	Scaffold0:2348604-2351174:-
Funhe2EKm000124t1	Funhe2EKm000124t1_C1	gspl3n5h	gspl3n5h	JXMV01064242.1:7209-7396:+
Funhe2EKm000137t1	Funhe2EKm000137t1	gmap2a5u	gmap2a5u	Scaffold0:3169380-3169568:+
Funhe2EKm000143t1	Funhe2EKm000143t2	gmap3n5h	gmap3n5h	KN805525.1:3239634-3244459:+
Funhe2EKm000182t1	Funhe2EKm000182t2	gmap3n5h	gmap3n5h	KN805525.1:4918131-4956066:-
Funhe2EKm000191t1	Funhe2EKm000191t2	gmap3n5h	gmap3n5h	KN805525.1:5110917-5134065:+
Funhe2EKm000203t1	Funhe2EKm000203t1	gmap2a5u	gmap2a5u	Scaffold0:5888005-5890994:-

  cat  ../gmapn/findcds4upd/kfish2all5_fcds45.besttab3 | egrep 
AGeneID          	AQueryID	BestMap	BestAa	hodiff	hoval	aadiff	cds2alen	cdsoff	cov	pid	maploc	path	othersrc
Funhe2EKm000019t1	Funhe2EKm000019t3	gmap3n5h	gmap3n5h	100	98%,908,azmolly:XP_007548373	0	599,61%,complete	96-1895	100	99	KN805525.1:267309-281534:+	0	gmap3n5h,gmap2a5u,gspl3n5h
Funhe2EKm000041t1	Funhe2EKm000041t1_C1	gspl3n5h	pubaa	0	0	-174	70,25%,partial5	3-215	52	99	KN805525.1:941170-941979:+	Split=1/2	gspl3n5h,gspl3n5n,gmap3n5h,none,gspl3n5hfx,gspl3n5nfx,gmap2a5u
   
  ----     
   env  test=1 vers=fc11t debug=1 ./evgkfish2subgff.pl >& gmapnsub2/log.evgsubfc11t  # test set   
   env vers=fc11f debug=1 ./evgkfish2subgff.pl >& gmapnsub2/log.evgsubfc11f # all
   
    cat gmapnsub2/log.evgsubfc11f | egrep -v '#errnoan|inrow' 
#in.skipids: src=kf2a, gmapn/kfish2rae5h_kfish2a_findcds3.cantuse.ids
#asmbest.in: gmapn/findcds4upd/kfish2all5_fcds45.besttab3
#scaftab.in: pubgenome/ncbifunhe3_kfish2asm.sametab
#pubaa.in: pubtsa6/kfish2rae5x11t6_tsasubmit.aa.gz
#gff.out: gmapnsub2/kfish2rae5h_asm3n_fc11f.gff
#gff.in: STAG=gmap2a5u,SMAP=kf2a gmapn/findcds4upd/kfish2rae5ht1asm2_fcds4h.gff.gz
#gff.out: in=37063,put=4728,drop=0,pubaa=2136
#gff.in: STAG=gmap2a5u,SMAP=kf2a gmapn/findcds4upd/kfish2rae5htaasm2_fcds4h.gff.gz
#gff.out: in=78520,put=0,drop=0,pubaa=0
#gff.in: STAG=gmap3n5h,SMAP=kf3n gmapn/findcds4upd/kfish2rae5h_asm3n_fcds4.gff.gz
#gff.out: in=120755,put=23138,drop=230,pubaa=3272
#gff.in: STAG=gspl3n5n,SMAP=kf3n gmapn/findcds4upd/kfish2nsplign15n_fcds4h.gff.gz
#gff.out: in=18694,put=345,drop=0,pubaa=148
#gff.in: STAG=gspl3n5h,SMAP=kf3n gmapn/findcds4upd/kfish2nsplign15h_fcds4.gff.gz
#gff.out: in=13667,put=767,drop=0,pubaa=420
#gff.in: STAG=gmap3hfx,SMAP=kf3n gmapn/findcds4upd/fixset/kfish2rae5h_gmap3n.fcds5.gff.gz
#gff.out: in=5096,put=0,drop=0,pubaa=0
#gff.in: STAG=gspl3n5hfx,SMAP=kf3n gmapn/findcds4upd/fixset/kfish2nsplign15h.fcds5.gff.gz
#gff.out: in=4207,put=409,drop=0,pubaa=331
#gff.in: STAG=gspl3n5nfx,SMAP=kf3n gmapn/findcds4upd/fixset/kfish2nsplign15n.fcds5.gff.gz
#gff.out: in=5012,put=220,drop=0,pubaa=179
#gff.totout: in=283014,put=29607,drop=230
#idactions: gmapnsub2/kfish2rae5h_asm3n_fc11f.idtab
#sortGff: /bio/bio-grid/mb/evigene/scripts/bestgenes_update.pl -act=sort -debug -cadd addgene=1 -vers kf2rae5h -conf pubgenome/evigene_kfish2_gbsubmit.conf -in gmapnsub2/kfish2rae5h_asm3n_fc11f.gff -out gmapnsub2/kfish2rae5h_asm3n_fc11fans.gff


=cut

=item exonMapErr
  
  FIXME V2: add MGap errors from gsplign to this misc_feature, locus map error listing
  FIXME: stutter, have same error from multiple alt exons .. 
  # add here? as Region/misc_feature following mRNA/CDS ..
   exon annots: 13651 error from gmap; keep/report in gb.annot? Region=map-exception...
   Parent=Funhe2EKm000208t1;error=ERROR.span:genome_span:1070,tr_span:1241,KN805525.1:5989263-5988193171
   Parent=Funhe2EKm000500t1;error=ERROR.span:genome_span:3687,tr_span:3875,KN805527.1:3629399-3625712188;
   ERROR: Bad location on feature misc_feature (start 3629398, stop -669255109)

    # possible misc attr:  /inference="[CATEGORY:]TYPE[ (same species)][:EVIDENCE_BASIS]"  
    # gene == /locus_tag="text" 

  from gmap: Funhe2EKm009624t3
     misc_feature    19132..19282
                     /locus_tag="D326_009624"
                     /note="trmap_genome_error:genome_span:150,tr_span:363"
  from gsplign: Funhe2EKm009624t2 MGap:1145-1351 b/n exon3:469-1144,exon4:1352-1429
     misc_feature    20055..21252
                     /locus_tag="D326_009624"
                     /note="trmap_genome_error:tr_align_gap:1145-1351"
    
  fixme V2 update: exonMapErr: locus_tag wanted, diff from gene gid, get from mrna
  
  ?? cancel dup maperr locations?
  WARNING: valid [SEQ_FEAT.DuplicateFeat] Features have identical intervals, but labels differ FEATURE: 
    misc_feature: map_error:genome_span:413,tr_span:768 [lcl|KN811437.1:22476-22889] [lcl|KN811437.1: delta, dna len= 1503504]
  
    same genome locus for 3 alts of locustag: D326_033565, but diff tr spans
  grep  '^22476  22889'  kfish2rae5h_fc14_knset8.tbl
  22476	22889	misc_feature <<?? gff says 22478<< not 76 ..	22889
        note	map_error:genome_span:413,tr_span:579 for Funhe2EKm033565t1
  22476	22889   # mrna exon w/o map err?
  22476	22889	misc_feature
        note	map_error:genome_span:413,tr_span:768 for Funhe2EKm033565t2
  
=cut


=item sub exonMapErr from evigene2genotbl_kfish2.pl

#?? move from evigene2genotbl_kfish2 to here? 
# .. need read all mRNA/exons to bundle before this calc
# .. cant do here yet

my(%didxerr);  

sub exonMapErr {
  my($exons,$geneid,$mrna)= @_;
  ## add gsplign MGap handling: now in mRNA annot as Target spans, match to exons/between-exons
  ## use same tag? eg. miscfeat span = exon1e .. exon2b, Note=err:trmapgap:[MGap 200-300]
  my @xerr= grep{ $_->[8]=~ m/error=ERROR.span/ } @$exons;

use constant CHECK_SPLIGN_EXON_ERRS => 1; # TEST
if(CHECK_SPLIGN_EXON_ERRS) {  
  if(not @xerr and ref $mrna) { # gsplign MGap ?
    my($gp)= $mrna->[8] =~ /gaps=([^;\s]+)/;
    ## Also see exon annot internal mgap trg=Funhe2EKm000098t1 1 1081;mgap=633-754;
    ## matches mRNA gaps=459,MGap:633-754,...
    if($gp=~/MGap/) { 
      my @mg= $gp=~m/MGap:([^,;]+)/g; # gaps=580,MGap:1371-1487,MGap:1000-1148,MGap:398-547,MGap:4545-4708,;
      my %mg=(); map{ my($b,$e)= m/(\d+).(\d+)/; $mg{$b}=$e; } @mg;
      my($ltb,$lte,$ltrb,$ltre)=(0,0,0,0);
      for my $x (@$exons) {
        my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$tid)= @$x;
        if( my($trb,$tre)= $tattr=~m/(?:Target|trg)=\S+.(\d+).(\d+)/ ) {
          for my $mb (sort keys %mg) { 
            my $me=$mg{$mb}; my $twid=1+$me-$mb; my($gb,$ge)=(0,0);
            if($mb<$tre and $me>$trb) { # inside exon
              if($to eq '-') { }
              else { $gb=$tb+($mb-$trb); $ge=$tb+($me-$trb); } ## REV STRAND needs te - xxx
            } elsif($ltb and $me<$tre and $mb>$ltb) { 
              if($to eq '-') { }
              else { $gb=$ltb+($mb-$ltb); $ge=$ltb+($me-$ltb); }
            }
            if($gb>0 and $ge>$gb) {
              my $xer="genome_gap,tr_span:$twid,$ref:$gb-$ge";
              $tattr =~s/$/;error=ERROR.span:$xer/; $x->[8]=$tattr; push @xerr,$x; last; 
            }
          }
        ($ltb,$lte,$ltrb,$ltre)= ($tb,$te,$trb,$tre);
        }
      }
    }
    # hassle need to match MGap tr-spans to exon Target/trg spans;
    # should revise gsplign2gff to keep MGap on exons
    # mimic: error=ERROR.span:genome_span:2315,tr_span:2428,KN811437.1:691859-689544,
  }
}
  
  return(0) unless(@xerr);
  my $locustag= (ref $mrna and $mrna->[8] =~ /locustag=([^;\s]+)/) ? $1 : 0;
  
  my @xloc=();
  for my $x (@xerr) {
    my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$tid)= @$x;
    my($xerr)= $tattr =~ m/error=ERROR.span:([^;\s]+)/;
    if($xerr) {
      my($gs,$ts,$gloc)=split",",$xerr;
      my($gr,$gb,$ge)= $gloc =~ m/(\S+):(\d+)-(\d+)/; #** ERROR.span: gb,ge can be rev, ignore? and use to
      ## WHAT THE FK it is giving rev strand from exon to misfeat *!*!*&(!@
      my $go=$to; # which strand?  #my $go= ($gb>$ge)?'-':'+'; # which to believe?
      if($to eq '.') { $go= ($gb>$ge)?'-':'+'; }
      ($gb,$ge)=($ge,$gb) if($gb>$ge);
      
      my $gsp= $ge - $gb;  my $xsp= $te - $tb;
      #? check $gr eq $ref ; maybe bad span?
      if($ge>0 and $gr eq $ref and $gsp < 1.1*$xsp) {
        my $xat=""; #"Note=trmap_genome_error:$gs,$ts"; 
        if($locustag) {
          $xat .= "locustag=$locustag;"; #? or Parent=pid ??
        } elsif($geneid) {
          $xat .= "gene=$geneid;"; #? or Parent=pid ??
        }
        
        #old# $xat .= "Note=trmap_genome_error:$gs,$ts;"; 
        $xat .= "Note=map_error:$gs,$ts;"; 
        #o# my $xerrid=join",",$gloc,$geneid,$tb,$te; #? maybe skip tb,te part?
        my $xerrid=join",",$gloc,$geneid; #?  skip tb,te part?
        next if($didxerr{$xerrid}++); # dont stutter.. SEQ_FEAT.DuplicateFeat.Features have identical intervals, but labels differ
        
        ## FIX strand, "." no good needs to be same as mRNA : FIX is worse .. all -strand bad strand?
        my $xl= [$gr,$src,"misc_feature",$gb,$ge,1,$go,0,$xat,$tid];
        push @xloc,$xl;
        }
    }
  }
  my $nx=@xloc;
  return ($nx,@xloc);
}

=cut

