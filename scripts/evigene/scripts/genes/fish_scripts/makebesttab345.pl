#!/usr/bin/env perl
# makebesttab345.pl
# /bio/bio-grid/kfish2/submitf/ 
# from submitf/gmapn/findcds4upd/makecdstab.sh and kfish2/submitf/evgkfish2subgff.pl

use strict;

my $S2='kf2a'; my $SANN='kf2rae5|kf2rae5alt'; # == kf2a, kfish2rae5h_kfish2a_findcds3t1.gff
my $S3='kf3n'; my $SOUT='kf2rae5h'; # == kf3n, kfish2rae5h_kfish2n_findcds3.gff
my $SPUB='pubaa'; # pub.aa replaces cds2aa

my $debug= $ENV{debug}||0;
my $doalt= $ENV{alt}||0; # separate tables for main/alts ?? but some mixing for best t1/main == alt

my $HODIFFMIN= $ENV{hodiff} || 95; # percent of pubaa-hoscore
my $AAabsDIFFMIN= $ENV{aadiff} || -49; # abs-aasize not pct-aa ??
my $AAPercentDIFFMIN= $ENV{aapctdiff}||75; # OLD AADIFFMIN, not used now?
    #^^ makebestalt call, pay attention to AADIFF qual, use pubaa where cds2aa is short/poor match
    # 44000/110164 are >= 80%; 60000/110000 are >= 50%
    # note this aadiff % may be flaky calc. subsitute cds2aa.sizeNoX / pubaa.sizeNoX ?
    
my $xidtab="gmapn/kfish2x11tsa.pubidtab";
my $augidtab="gmapn/kfish2augmod.pubidtab"; # these should come from kf2.gff direct, no cds2aa
my $auggff="gmapn/kfish2rae5h.augmod.gff.gz";
my $uptab="gmapn/kfish2rae5h.update1.tab"; # changes ie. drop/droplocus/replace/... as per genes update

# FIXME: add erridtab to exclude sources as best, from evgkfish2subgff.pl
my $ErrorIdTab= "gmapn/kfish2rae5h_kfish2a_findcds3.cantuse.ids";
  
  # ?? create new cds2aa.qual tables, parse from cds2aa.gff inputs?
my @inaaqual= qw(pubgenes/kfish2rae5h.main.pub.aa.qual pubgenes/kfish2rae5h.alt.pub.aa.qual);
  ## skip these? get from mRNA.gff parse
  ##  gmapn/kfish2rae5h_kfish2a_findcds3t1.cds2aa.qual gmapn/kfish2rae5h_kfish2a_findcds3ta.cds2aa.qual 
  ##  gmapn/kfish2rae5h_kfish2n_findcds3.cds2aa.qual
  
  ## add inhoscore, which variant?  
  # fish4ref.bltop4 has best hoscore for each cds2aa/pubaa model tested, t1 + some alts
  # fish4ref.best1q has best single model, .. best1p, bltop2 ..
  ## FIXME: srctag differs from STAG : gspl3n5h > gspl5h .. etc; need HOTAG equiv
my @inhoscore=qw( gmapnbl/aaeval2/kfish2rae5h_bestfcds45pubaa-fish4ref.bltop4);

  # add scaf rename Scaf > KNcbi
my $scafidtab="pubgenome/ncbifunhe3_kfish2asm.sametab";
  
  ## drop inaadiff
# my @inaadiff= qw(gmapn/kfish2rae5ht1_asm2a2n.diffaa.tab gmapn/kfish2rae5hta_asm2a2n.diffaa.tab);

  ## replace in map.attr w/ in cds2aa.gff and parse out map.attr?
my @inmapattr=qw(gmapn/kfish2rae5hpub-kfish2n.map.attr gmapn/kfish2rae5h_asm2.map.attr);

## 2015.04.25: add asm2 gmap lifted to asm3n ..  kfish2rae5htasm3lft_fcds7 from agplift2gff
## 2015.04.22: findcds upd fixes > fcds7 set
## 2015.04.15: findcds upd fix rev xtrim=1 bugggg > fcds6 set
## 2015.04.11: add ingff.cdna.endgaps table addition to these, need to trimNNN end gaps missed by findcds
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
    
    
=item usage

  env  debug=1 ./makebesttab345.pl gmapn/findcds4upd/kfish2rae5h_asm3n_fcds4.gff.gz |\
     sort > gmapn/findcds4upd/kfish2rae5h_asm3n_fcds4.besttab345

  # redo, all gff input, one table output w/ cds2aa best src + scores + samebest,poor srcs?
    @row= ($aadiff,@at,$maploc,$issplit,$STAG); # all for cds2aa hash? or just ROWT?
    $cds2aa{$id}{$STAG}= \@row;


=item inputs

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

=cut
    

#--------

# BEGIN  
my( %cds2aa, %sources, %cds2aaswap, %maplocswap, # output hash,
    %awv, %aaq, %cdsq, %cdsv, %hoval, %hoscore, %skipids, %maploc);
my(%scaf,%scafn);
my($tsaids,$nxpd,%xpd,%augmod,%updates);  
my(@COLNAMES);

my $STAG="none"; my $SMAP="none"; my $HOTAG=$STAG; my $HOTAG2="";
  ($STAG,$SMAP,$HOTAG)= sourceTag($ENV{src} || $ENV{pt}); # set STAG from $ingff name 

readUpdatetab($uptab);

readAaQuals(@inaaqual);
readTSAidtab($xidtab);
readAUGMODidtab($augidtab);
# "gmapn/kfish2augmod.pubidtab"; # these should come from kf2.gff direct, no cds2aa

readScaffoldTab($scafidtab);  
readHoScore(@inhoscore);
checkGffErrors($ErrorIdTab);

my $outh= *STDOUT;

my @tingff= grep/gff/, @ARGV;
   @tingff= @ingff unless(@tingff);

for my $ingff (@tingff) {
  my($ok,$inh);
  warn "#malt.in: $ingff\n" if($debug);
  
  if($ingff =~ /.gz/) { $ok=open($inh,"gunzip -c $ingff |"); }
  else { $ok=open($inh,$ingff); }
  unless($ok){ my $err="#ERR: missing ingff $ingff\n"; if($debug) { warn $err; next; } else { die $err; } }
  ($STAG,$SMAP,$HOTAG)= sourceTag($ingff); # set STAG from $ingff name 
  my($ngapid,$endgaphash)= readEndGaps($ingff);    
  cds2aaGffToAaMapTable($inh,$outh,$endgaphash); #
}

cds2aaBestOut($outh); # if (%cds2aa);

# ?? makeAltBestTab($outf,$outh);
#-------------------------------------------------

# GeneID  alt.src		Bits	Iden	Algn	Qlen	Slen	Algb	Mism	Ibpt  RefID
# Funhe2EKm000003	t1.c1kf3n	358	174	209	234	214	214	31	2,1	platyfish:ENSXMAP00000008562
# Funhe2EKm000003	t1.pubaa	358	174	209	234	214	214	31	1,1	platyfish:ENSXMAP00000008562
# Funhe2EKm000004	t1.pubaa	1909	935	1102	1102	1103	1103	166	1,1	azmolly:XP_007548354

# @ingff try3:
#malt.in: gmapn/findcds4upd/kfish2nsplign15n_fcds4h.gff.gz
#malt.in: gmapn/findcds4upd/kfish2nsplign15h_fcds4.gff.gz
#malt.in: gmapn/findcds4upd/kfish2rae5h_asm3n_fcds4.gff.gz
#malt.in: gmapn/findcds4upd/kfish2rae5ht1asm2_fcds4h.gff.gz
#malt.in: gmapn/findcds4upd/kfish2rae5htaasm2_fcds4h.gff.gz
#malt.in: gmapn/findcds4upd/fixset/kfish2nsplign15h.fcds5.gff.gz
#malt.in: gmapn/findcds4upd/fixset/kfish2nsplign15n.fcds5.gff.gz
#malt.in: gmapn/findcds4upd/fixset/kfish2rae5h_gmap3n.fcds5.gff.gz
#HOTAGs 19756 c1kf3n  4594 c1kf2a     2663 pubaa
#  816 gspl5h  465 gmap3h
#  446 gmap3hfx  418 gspl5hfx
#  348 gspl5n  261 gmap2h  166 gspl5nfx

sub sourceTag {
  my($pt)= @_;
  ## ?? let pt == stag from table, set smap, hotag
  
  ## buggers: kfish2rae5h_gmap3n.fcds5 no stag/smap
  ## FIXME: add tag for $augmod{pubid} set .. these 4,000 are always best t1 set; pt =~ /augmod/ ??
  # gmapn/kfish2rae5h.augmod.gff.gz
  
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
  elsif($pt=~/x11tsa/) { $stag="gmap3nx11tsa"; $hotag="none"; } 
   
  if($pt=~/gapfix/) { $stag="gmap3updx"; $smap=$S3; $hotag="c1kf2a"; $HOTAG2="c1kf3n"; } # FIXME: old or new tag?

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
    my($pd,$oid,$act,$upd,$comm)=@v; next unless($act);
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
  open(F,$augidtab) or die "#MISS: $augidtab\n"; 
  while(<F>){ next if(/^\W/); my($pd,$oid)=split; $augmod{$pd}=$oid; $augids++; } close(F); 
}

sub readTSAidtab {
  $tsaids=0; # TSA pubid == GCES01000000 + x11 idnum
  open(F,$xidtab) or die "#MISS: $xidtab\n"; 
  while(<F>){ next if(/^\W/); my($xd,$pd)=split; $xpd{$xd}=$pd; $tsaids++; } close(F); 
}

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

sub readHoScore {
  for my $inf (@inhoscore) {
    warn "#malt.in: $inf\n" if($debug);
    open(F,$inf) or die $inf;
    while(<F>) {
      next if(/^\W/ or /^Query/); 
      chomp; my @v=split"\t";  
      my($gd,$tisrc,$bits,$iden,$algn,$reflen,$tlen,$alnb,$mism,$ibpt,$refid)= @v; 
      my($ti,$src)= split/\./,$tisrc;
      ## src == HOTAG not STAG, but $SPUB == hotag ?
      my $td=$gd.$ti; 
      if($td ne $gd.'t1') { } #what? add as t1 unless t1?
      # $hoscore{$td}{$src}=$bits; 
      $hoval{$td}{$src}= [$bits,$iden,$algn,$reflen,$refid]; # pctalgn= 100*$algn/$reflen
    } close(F);
  }
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
    # regen Missing-mrna-exons, per gff source, not all S2..
  open(F,$erridtab) or die $erridtab;
  while(<F>) { 
    next if(/^\W/); my($id,$err)=split; 
    if($err =~ /Scaf.change/) { $err=~s/err=//; $err=~s/-/_/g; $skipids{$id}{$S2}= $err;  $nerr++;}
    else { } # fixme
  } close(F);
  warn "#in.skipids: n=$nerr from $erridtab\n" if $debug;
}

=item readEndGaps

  my($ngapid,$endgaphash)= readEndGaps($ingff);

kfish2rae5ht1asm2_fcds4h.cdna.endgaps
AGeneID         	TrSize	Cdsoff	Gapbeg	Gapend
Funhe2EKm003867t1	620	252-521	0	520-619
Funhe2EKm004111t1	2130	584-2128	0	2025-2124
Funhe2EKm010445t1	512	231-512	0	403-502
Funhe2EKm013892t1	1878	141-1790	0	1788-1877

    findcds4upd/
    cds  cdna
     15     4 kfish2nsplign15h_fcds4.cdna.endgaps
     11    35 kfish2nsplign15n_fcds4h.cdna.endgaps
    128   213 kfish2rae5h_asm3n_fcds4.cdna.endgaps
     70   480 kfish2rae5ht1asm2_fcds4h.cdna.endgaps
     75   374 kfish2rae5htaasm2_fcds4h.cdna.endgaps
    findcds4upd/fixset/
    cds   cdna
      2      2 kfish2nsplign15h.fcds5.cdna.endgaps
      2      2 kfish2nsplign15n.fcds5.cdna.endgaps

=cut

sub readEndGaps {
  my($ingff)= @_;
  my(%gaps); my $ngap=0;
  my $intab=$ingff; $intab=~s/\.gff.*/.cdna.endgaps/;
  warn "#endgaps.in: $intab\n" if($debug);
  open(F,$intab) or warn "#endgaps.in: miss=$intab\n";  
  ##FIXME ?? look for .cdna if no cdna.endgaps and make it??
  while(<F>) { 
    next if(/^AGene|^\W/);
    my($id,$trlen,$off,$gb,$ge)=split;
    my($ob,$oe)=split"-",$off; my($gbb,$gbe)=split"-",$gb; my($geb,$gee)=split"-",$ge; 
    ## ignore gaps inside cds, should have been checked/corrected, and ncbi Cheka should pass
    ## CANT ignore cds-end-gaps at very end, only 1-codon (3 bp) inside
    if($gbe>0) { if($gbe>$ob and $gbb>$ob+2) { } else { $gaps{$id}.="gapb=$gbe/$gbb,"; $ngap++; } }
    if($gee>0) { if($geb<$oe and $gee<$oe-2) { } else { $gaps{$id}.="gape=$geb/$gee,"; $ngap++; } }
    ## $gaps{$id}{beg}=$gbe;  $gaps{$id}{end}=$geb;
    #or# $gaps{$id}.="beg=$gbe,beg0=$gbb,"; $gaps{$id}.="end=$geb,end2=$gee,";
  } close(F);
  return($ngap,\%gaps);
}

sub readAaQuals {

  for my $inf (@inaaqual) {
    warn "#malt.in: $inf\n" if($debug);
    open(F,$inf) or die $inf;
    
    ## FIXME: Selcstop=1 prots *** special handling ** lost now as not bestaa *****
    ## annot not in aaqual yet, in >aahdr 
    
    ## FIXME: $STAG, not global, parse from inf name
    my $src=($inf=~/pub.aa/)? $SPUB : ($inf=~/kfish2n/)? $S3 : ($inf=~/kfish2a/)? $S2 : $S3;
    # if($STAG =~ /gmap2pub5h|gmap2pub5u/) { $src=$S2; } 
    # elsif($STAG =~ m/gspl3n5h|gspl3n5n|gmap3npub5h|gmap3nx11tsa/) { $src=$S3; }
    my($stag,$smap)= sourceTag($inf); # set STAG from $ingff name 
    #?? $src=$smap if($smap);
    
    while(<F>) {
      next if(/^\W/ or /^AQueryID/); 
      chomp; my @v=split"\t";  
      my($td,$aw,$gp,$aq,$clen,$offs)= @v; 
      if($gp>0 and not($aq =~ /gap/)) { $aq.=",gaps:$gp"; $v[3]=$aq; } 
      $offs=~s/:.//; $v[5]=$offs;
      $awv{$td}{$src}= $aw; $aaq{$td}{$src}= [@v];
        # ^^ == $pubaw{$td}=$aw; $pubaq{$td}=[@v];  
      
      # if($src eq $SPUB and $aq=~/,selc/) { } # FIXME
      
      if($src ne $SPUB and (my $awp= $awv{$td}{$SPUB}) ) {
         my $cdsp= int(0.5 + 100* $aw/$awp);  
         my $aacl= $cdsv{$td}{$src};
         $aacl=~s/^(\d+)%/$cdsp%/;
         $cdsq{$td}{$src}= $cdsp;
         $cdsv{$td}{$src}= $aacl; 
         }
    } close(F);
  }
}


sub cds2aaGffToAaMapTable {
  my($inh,$outh,$endgaphash)= @_;
  my ($nskipdup,$nswapmain,$nput)=(0) x 9;
  while(<$inh>) {
  if(/\tmRNA/) { 
    my($ids,$aw,$ao,$issplit);
    chomp; my @v=split"\t"; 
    my($ref,$rb,$re,$ro,$at)=@v[0,3,4,6,8]; chomp($at);
    my $hasprot= ($at =~ m/;protein=/)?1:0;
    #? next unless($hasprot); # pick only cds2aa w/ this for splits? defer this and collect all split maplocs ?

    my $SCAFOK=1; my($rnew,$rold);
    if($rnew= $scaf{$ref}{newna}) { 
      $SCAFOK= $scaf{$ref}{same}; #  # problem .. try to shift? need agpLift !
      $rold= $ref; $ref= $rnew;  s/\b$rold\b/$rnew/g; 
      # $anchange .= ";scold$SCAFOK=$rold"; 
    }

    my $maploc="$ref:$rb-$re:$ro";
    
    my($id)=m/ID=([^;\s]+)/; ($ids=$id)=~s/_C\d+$//; 
    next if($id =~ /_G\d+$/); # skip all 2nd gmap locs
    
    if($STAG eq "gmap3nx11tsa" and $nxpd) { my $pd=$xpd{$ids} or next; $id=~s/$ids/$pd/; $ids=$pd; } #?? 
    unless(($issplit)=m/(Split=[^;\s]+)/){ $issplit=($id=~/_C(\d+)/)?"Split=$1":0; }  

    $maploc{$ids}{$STAG} .= "$maploc,"; #? use this for all splits, alts, 
    
    # if(my $augoid=$augmod{$ids}) { } # require orig map-source.gff not findcds set
    ## augmod2 == STAG; 
    my $isaug= ($STAG =~ /augmod/)?1:0;
    
    my $maperr= $skipids{$id}{$SMAP} || ""; 
    $maperr .=",err:Scaf_change" unless($maperr or $SCAFOK);
    if($at =~ /err=(Missing.mrna.exon)/) { my $e=$1; $e=~s/-/_/g; $maperr .=";$e";}
    # other maploc errs may go here: Split C1: same-loc, rev-strand wont go past ncbi checks: ERR.MixedStrand
    $maploc.=",err:$maperr" if($maperr); # paste serr onto maploc,err=xxx
    
    ## BUG: 100 Scaf_change are creeping in to BestMap even w/ alt lift to asm3n set.. 
    my $CANTUSE= ($SCAFOK and not $maperr=~/Scaf.change/)?0:1;
    
    my $endgaps= $endgaphash->{$id} || $endgaphash->{$ids}; #?? add to maploc, but not as err:, as gaps:?
    if($endgaps) { $maploc.=",egap:$endgaps"; }

    # augmod.ann: D=Funhe2EKm000022t1;oid=kf4bAUGpia7bp1s0g19t1;
    # aalen=259,60%,complete;osrc=kf4best;Dbxref=UniProt:Q76IM1_DANRE;Name=Pol-like protein;nameln=15%,203/1329,259;
    # homolog=89.1,UniRef50_R4G9D5;inqual=0;nintron=0/8;tere=40;cxlen=780/1278,61%;
    if($isaug) {
      $at .= ";cov=100;pid=100"; # no offs,cdsoff; clen == cxlen.2
    }
    
    ## add flags col? for errors, other?
    @COLNAMES= qw(AQueryID hodiff hoval aadiff cds2alen pubalen cdsoff puboff clen cxlen cov pid maploc path maptag); 
    my @at=(); 
    for my $k (qw(aalen aaold cdsoff offs clen cxlen cov pid)) { 
      my($v)= $at =~ m/\b$k=([^;\s]+)/; $v||=0; 
      if($k =~ m/cov|pid/) { $v=~s/,.*//; $v=~s/%//; $v= int(0.5 + $v); }  
      push @at,"$k=$v"; 
    } 
     
    unless( ($aw)= $at=~ m/\baalen=(\d+),\w/){ ($aw)=m/\baalen=(\d+)/; }
    if($hasprot) {
      my($aa)= $at=~m/protein=([^;\s]+)/; 
      my($xa)= $aa=~tr/Xx/Xx/; if($xa) { $aw -= $xa; $at[0].=",gaps:$xa"; }
    }
    
    my $Selcstop=0;
    if(my $pubv=$aaq{$ids}{$SPUB} ) { # my $pubv=$pubaq{$ids}
      my($d,$w,$g,$aq,$cl,$off)=@$pubv; 
      $ao=$w; $at[1]="aaold=$aq"; $at[4]="clen=$cl"; $at[3]="offp=$off"; 
      $Selcstop=$1 if($aq=~m/,(selc\w*)/i);
      if( $Selcstop) { $hasprot=1; $at[0].=",$Selcstop" unless($at[0]=~/selc/i); }
    } else { 
      ($ao)=(m/\baaold=(\d+)/)?$1:$aw; # fail instead of using aw?
    }

    my($hodiff,$hoscore)=(0,0); # missing
    ## BUG in HOTAG use, some ingff may be one of 2 HOTAGs, changed cds2aa data?
    ## ........
    my $hoval=$hoval{$ids}{$HOTAG};
    if(not $hoval and $HOTAG2) { $hoval=$hoval{$ids}{$HOTAG2}; }
    #?? seem to be missing hoval for valid alts/main .. why? ids,hotag mismatch?
    
    if($hoval) { 
      my($bits,$iden,$algn,$reflen,$refid)=@$hoval; # $bits,$iden,$algn,$reflen,$refid
      my $paln= int(0.5 + 100*$algn/$reflen); $paln=100 if($paln>=100);
      $hoscore="$paln%,$bits,$refid";
      if( my $hopub=$hoval{$ids}{$SPUB} ) {
        my($pbits,$piden,$palgn,$preflen,$prefid)=@$hopub; # $bits,$iden,$algn,$reflen,$refid
        $hodiff= int(0.5 + 100* $bits/$pbits); #? check prefid == refid or just use bits as top val?
        # ^^ bad score to compare best, want abs bits 1st
      }    
    }
## old asmbest HoDiff : duplicate this? use pct(hocds2aa/hopubaa) ?
# Funhe2EKm000005t1.pubaa top,26.a,azmolly:XP_007548349   
# Funhe2EKm000007t1.kf3n  samedi,-3.a,azmolly:XP_007548366
# Funhe2EKm000008t1.kf3n  same,0.a,mayzebr:XP_004571788.1 

    my $aadiff= $aw - $ao; #? change to percent? both?
    #my $paadiff= ($ao>0) ? int(0.5 + 100*$aw/$ao) : 99;

    ## convert this reader to @inmapattr and AaDiff for makeAltBestTab ?
    ## diffaa > AaDiff;
    ## add col: maploc=GenomeID:gespan:geor 
    ## add col: HoDiff
    ## drop cols? cdsoff,puboff  pubclen,cxlen ??
    ## ?? drop pub cols here, use separate table, ie pub.aa.qual has these, no valid cds maploc for pubaa
    #o# AQueryID diffa cds2alen pubalen cdsoff puboff pubclen cxlen cov pid path maptag

    if($CANTUSE or not($hasprot or $isaug)) {
      # pick only cds2aa w/ this for splits? 
      # ? defer this and collect all split maplocs ?
      # $cds2split{$ids}{$STAG}= what??
      next;
    }
    
    map{ s/^\w+=//; } @at; 
    ## $ids or $id ?
    ## ?? fix here change talt to t1 if has hoscore? or add t1 of talt?
    ## FIXME: squeeze in $maperr flags, reduce score
    ## FIXME: bestaa=pubaa for Selcstop ***
    my @row= ($id, $hodiff, $hoscore, $aadiff, @at, $maploc, $issplit, $STAG); # all for cds2aa hash? or just ROWT?
    
    if(my $valo= $cds2aa{$ids}{$STAG}) {
      ## BAD??, dont do swap/skip here ?
      if(not $hoscore and $valo->[2]) { 
        my $idso= $valo->[0];
        if($cds2aa{$idso}{$STAG}) { 
          # warn "#malt.in: skip duprec $ids,$idso\n" if($debug); # lots of these, 471 below?
          $nskipdup++; next; 
        } else { 
          # warn "#malt.in: swapmain2 $ids,$idso\n" if($debug); # none of these
          #? $row[-2].=",swap2:$ids";
          #? $ids= $idso;  #swap main
        } 
      }# has old, swapmain ??
    }
    
    $cds2aa{$ids}{$STAG}= \@row; $nput++;
    $sources{$STAG}++;
    
    if($hoscore and not $ids=~/t1$/) {
      my $idso= $ids; $idso=~s/t\d+$/t1/; # FIXME output needs to use ids, and id,src
      # ^^ FIXME2, need to track these altid > mainid changes, and add reverse swapid (t1 > talt),
      # and dont put rows of two ids w/ same origid
      
      my $valo= $cds2aa{$idso}{$STAG};
      
      ## problem better Selcstop t1 main pubaa  lost here? need hopub
      # swapmain=Funhe2EKm010847t7,Funhe2EKm010847t1/Selenoprotein N
      # swapmain=Funhe2EKm036789t3,Funhe2EKm036789t1/Selenocysteine insertion sequence
      my $hot1best=0;
      if( my $hopub=$hoval{$idso}{$SPUB} ) {
         my($bits,$iden,$algn,$reflen,$refid)=@$hoval; # $bits,$iden,$algn,$reflen,$refid
         my($pbits,$piden,$palgn,$preflen,$prefid)=@$hopub; # $bits,$iden,$algn,$reflen,$refid
         my $hodifft1= int(0.5 + 100* $bits/$pbits);
         $hot1best=1 if($hodifft1 < $HODIFFMIN); # 95% now
      }
 
      if(my $upd= $updates{$idso}) {
        if($upd=~m/^swapbest=(\w+)/) { my $sbid=$1;
        if($sbid eq $ids) { $hot1best=0; $valo= 0; } # force swap this way ??        
        }
      }
     
      unless($hot1best or ($valo and $valo->[2])) { # 2: has hoscore 
        $cds2aa{$idso}{$STAG}=\@row;
        # FIXME: upd maploc{$idso} also ?? but 2 may not be splits?
        my $mapo= $maploc{$idso}{$STAG};
        $maploc{$idso}{$STAG}=""; $maploc{$idso}{$STAG} .= "$maploc,"; #? replace? not append?
        if($valo) { 
          $nswapmain++;
          # warn "#malt.in: swapmain1 $ids,$idso\n" if($debug); #  
          $valo->[-2].=",swap1:$idso/$ids"; # n=471 ?? this is bad now; need BestOut best selection before swap
          $cds2aaswap{$idso}{$STAG}= $valo; # swapmain rows for t2 alt
          $maplocswap{$idso}{$STAG}= $mapo if($mapo); # swapmain rows for t2 alt
       } 
      }
      # BUT on output need to indicate swap t1,t2.. 
    }
    
  } # mRNA
  }
  warn "#malt.in: nput=$nput, skip duprec n=$nskipdup, putative nswapmain=$nswapmain\n" if($debug); # lots of these, 471 below?
}

=item swapmain eg

AGeneID	AQueryID	BestMap	BestAa	hodiff	hoval	aadiff	cds2alen	cdsoff	cov	pid	maploc	path	othersrc
Funhe2EKm000139t1	Funhe2EKm000139t2	gmap2a5u	gmap2a5u	101	100%,476,azmolly:XP_007548524	-8	255,77%,complete	63-830	95	0	Scaffold0:3223767-3230956:-	0	gmap2a5u,gmap3n5h
Funhe2EKm000139t2	Funhe2EKm000139t1	gmap2a5u,swap1:Funhe2EKm000139t1	gmap2a5u	0	0	0	245,70%,complete	101-838	99	0	Scaffold0:3197624-3230994:-	0	gmap2a5u,gmap3n5h

Funhe2EKm000455t1	Funhe2EKm000455t6	gmap2a5u	gmap2a5u	106	100%,1430,azmolly:XP_007564497	-61	826,90%,complete	37-2517	82	0	Scaffold2:1999049-2008002:-	0	gmap2a5u,gmap3n5h,gspl3n5h
Funhe2EKm000455t6	Funhe2EKm000455t1	gmap2a5u,swap1:Funhe2EKm000455t1	gmap2a5u	0	0	0	817,96%,complete	16-2469	100	0	Scaffold2:1999191-2007981:-	0	gmap2a5u,gmap3n5h

=cut

# my(%maploc);  #? same as $splitmap{$ids}{$STAG}{$issplit} = $maploc
## use $maploc{$ids}{$STAG} .= "$maploc," and @splits= split",",$maplocs ?

sub locusaltmaperr { # from makeAltBestTab for cds2aaBestOut
  my($td,$src,$bestgid,$bestsrc,$swapbestid)=@_;
  my $err= "";
  my $gid=$td; $gid=~s/t\d+/t1/;
  
  my $loctd= $maploc{$td}{$src}; #? must have
  my $loct1= $maploc{$gid}{$src}; # same src or best src here??
  if($bestgid and $bestsrc) {
    $gid= $bestgid;
    $loct1= $maploc{$bestgid}{$bestsrc}; # same src or best src here??
    if($swapbestid) { # no effect, wrong id? maploc?
      #no# $loctd= $maplocswap{$swapbestid}{$src} || $loctd; #?? swapmain fix here?? maploc{bid} is proper loc,
      #no# $loctd= $maploc{$swapbestid}{$src} || $loctd; #?? swapmain fix here?? maploc{bid} is proper loc,
    }
  }
  return "" if($td eq $gid);
  ## FIXME: maploc should use {gid}{altid}{src} struct, check all alts of gid.
  ## for best selection of alts also need to check cross sources, or use some score independent of src
  
  if($loct1 and $loctd) {
    my @loct1= split",",$loct1;
    my @loctd= split",",$loctd;
    # if(@loct1>1 or @loctd>1) { } # handle splits .. test all @t1 vs td[0] ?? set err only if none agree
    my($jn,$kn)=(1,1);
    if(@loctd>1) { $kn=@loctd; } # ugh: @loct1>1 and @td..
    elsif(@loct1>1) { $jn=@loct1; } #?? elsif or both @ jn,kn
    my($okjk,$errlast,$errfirst)=("","","");
    for(my $k=0; $k<$kn; $k++) {
      for(my $j=0; $j<$jn; $j++) {
        my $errjk="";
        my($mr,$mb,$me,$mo)= split /[:-]/, $loct1[$j], 4; $mo=0 if($mo eq '.');
        my($gr,$gb,$ge,$go)= split /[:-]/, $loctd[$k], 4; $go=0 if($go eq '.');
        $errjk.="ERRscaf," if($mr ne $gr);
        $errjk.="ERRspan," unless($gb < $me and $ge > $mb);
        $errjk.="ERRrev," if($mo and $go and $mo ne $go);
        $errfirst= $errjk unless($errfirst); $errlast=$errjk;
        $okjk="$j,$k" unless($errjk or $okjk);
      }
    }
    unless($okjk) { $err.= $errfirst; }
    #...........
    # FIXME here, sometimes mvmain/t1 is wrong vs other alts .. esp from asmbestho.tab ?
    # my($mr,$mb,$me,$mo)= split /[:-]/, $loct1[$j], 4; $mo=0 if($mo eq '.');
    # my($gr,$gb,$ge,$go)= split /[:-]/, $loctd[$j], 4; $go=0 if($go eq '.');
    # $err.="ERRscaf," if($mr ne $gr);
    # $err.="ERRspan," unless($gb < $me and $ge > $mb);
    # $err.="ERRrev," if($mo and $go and $mo ne $go);
    
    #? if($err and ($path =~ /C3:/ or $path1 =~ /C3:/)) { $err=~s/ERR/split/g; }
  }
  return($err);
}
 
=item _sortgscore

  * revise to allow slop for hobits, hodiff, aadif, covpi ..
  ie. sort best w/o letting tiny diff in 1st score displace large diff in others
  use eg round(b.hobits/10) <=> round(a.hobits/10)
  
  # FIXME: hodiff bad sort 1st, need hobitscore absolute value 1st
  # .. shift up 1 more score, 0 = hobits

=cut

sub _sortgscore  {  
  ## new: $hobits,$hodiff,$aadiff,$covpi,$nsplit,$src,$errv, << first errv == col[6]
  return (
    $$b[6] <=> $$a[6]  # new: errv < 0 or > 0
    or $$b[0] <=> $$a[0] # NOW hobits, hodiff ?? use SLOP to better a few lower hodiff for higher map qual?
    or $$b[1] <=> $$a[1] # hodiff ?? use SLOP to better a few lower aadiff w/ higher mapqual?
    or $$b[2] <=> $$a[2] # aadiff  | merge these? ie lower covpi >> hi cov+split
    or $$b[3] <=> $$a[3] # covpi  | merge these? ie lower covpi >> hi cov+split
    or $$a[4] <=> $$b[4] # nsplit |
    or $$a[5] cmp $$b[5] # src
    );
}

sub cds2aaBestOut {
  my($outh)= @_;
  warn "#malt.out: cds2aaBestOut\n" if($debug);

  use constant ROWT => 3;# output slim table
  use constant cMAPLOC => 12;# index maploc col input

  my $didhdr=0;
  #above# @COLNAMES= qw(AQueryID hodiff hoval aadiff cds2alen pubalen cdsoff puboff clen cxlen cov pid maploc path maptag); 
  # @row= ($id, $hodiff,$hoscore, $aadiff, @at, $maploc, $issplit, $STAG); # all for cds2aa hash? or just ROWT?
  # $cds2aa{$id}{$STAG}= \@row;
  my @SCORECOLS=(1,2,3,10,11,12,13); # hodiff,hoval,aadiff,cov,pid,maploc,path
  
  ## fixme, need both ids (can differ): $cds2aa{$id} and row[0] .. also swap STAG to col2
  my @OUTCOLS;
  if(ROWT == 1) {
    my @atcol= map{ 4 + $_ }(0,1,2,3,6,7); # @at: 0..7
    @OUTCOLS= (0,3,@atcol,12,13,14);  # == $aadiff,@at[0,1,2,3,6,7],$maploc,$issplit,$STAG
    # @hdr=qw(AQueryID aadiff cds2alen pubalen cdsoff puboff cov pid maploc path maptag); # drop: clen,cxlen 
  } 
  if(ROWT == 2) {
    my @atcol= map{ 4 + $_ }(0,2,6,7);
    @OUTCOLS= (0,3,@atcol,12,13,14);  # == $aadiff,@at[0,2,6,7],$maploc,$issplit,$STAG
    # @hdr=qw(AQueryID aadiff cds2alen cdsoff cov pid maploc path maptag); # keep only aalen, cdsoff, cov, pid ?
  }
  if(ROWT == 3) {
    my @atcol= map{ 4 + $_ }(0,2,6,7);
    @OUTCOLS= (0,1,2,3,@atcol,12,13,14);  # == $hodiff,$hoscore,$aadiff,@at[0,2,6,7],$maploc,$issplit,$STAG
    # @hdr=qw(AQueryID hodiff hoval aadiff cds2alen cdsoff cov pid maploc path maptag); # keep only aalen, cdsoff, cov, pid ?
  }
  
  #was# AQueryID	hodiff	hoval	aadiff .. maptag
  #new# AGeneID  AQueryID BestMap BestAa  hodiff	hoval	aadiff
  my @outc= @OUTCOLS;
  @OUTCOLS=($outc[0],$outc[-1]); # ,@outc[1..-2] swap STAG to col2
  shift @outc; pop @outc; push @OUTCOLS, @outc;
  
  my @hdr= @COLNAMES[@OUTCOLS];
  $hdr[1]=~s/maptag/BestMap/;
  splice(@hdr,2,0,"BestAa"); # BestAa
  unshift(@hdr,"AGeneID");   # AGeneID?? LocusID ?? vs AQueryID
  push(@hdr,"othersrc");    # othersrc

  print $outh join("\t", @hdr)."\n" unless($didhdr++);

  my(%didput,$nskipdup,$nswap,$nswapb,%swapbest);
  my @src= sort keys %sources; #? ok sort
  
  # insert alt<>t1 mismap check here? ie alts should be same-loc-strand but for split parts

  sub vround { my($v,$d)=@_; return $d*int( 0.5 + $v/$d); } #? is this right? any d>1 for +/- v?
  
  my($bestGid, $bestGidsrc)=(0,0);
  for my $ids (sort keys %cds2aa) {
    my (%score, %sval); 
    my $gid=$ids; $gid=~s/t\d+//;
   
    if($bestGid) {
      my $g=$bestGid; $g=~s/t\d+$//; 
      unless($ids=~/^$g/) { ($bestGid, $bestGidsrc)=(0,0); }
    }
      # pick best among @src ...
    for my $src (@src) {
      my $vals= $cds2aa{$ids}{$src} or next; # no src val
      my($hodiff,$hoval,$aadiff,$cov,$pid,$maploc,$path)= @{$vals}[@SCORECOLS];
      my ($hobits)= $hoval=~m/,(\d+)/; $hobits||=0; # palign%,bitscore,refid
       ##  hodiff bad score to compare best, want abs bits 1st
      my $bid= $vals->[0]; $bid=~s/_C\d+$//; 
      
      my ($maperr)= $maploc=~m/,(err.\w*)/; ## remove ,err for splits?
      # my ($endgaps)= $maploc=~m/,(egap:.\w*)/; ## not error, but change score?
            
      my @maplocs= split",",$maploc{$ids}{$src}; 
      # if($bestGid and $maplocswap{$bid}{$src}){ @maplocs= split",",$maplocswap{$bid}{$src}; } ##? no good
      if($path =~ /swap1:/ and $maplocswap{$bid}{$src}) { @maplocs= split",",$maplocswap{$bid}{$src}; }
      
      if(@maplocs > 1) { # calc C1: map errs
        my %hasloc; my $sperr="";
        for my $sloc (@maplocs) { 
          my($sr,$sb,$se,$so)= split /[:-]/, $sloc, 4; $so=0 if($so eq '.');
          if($hasloc{$sr}) { $sperr.="C1:$sr:$sb-$se:$so,"; }
             #^? dont bother more checks? eg my($lb,$le,$lo)=@{$hasloc{$sr}}; .. overlap lb,sb le,se lo<>so
          push(@{$hasloc{$sr}},[$sb,$se,$so]);
          }
        if($sperr) { $maperr.= ($maperr)?",":"err:"; $maperr .= $sperr; }
      }
      
      ## ** locusaltmaperr  problem : src-specific, can pick t1=srcA, t2=srcB with maperr **
      ## fixme swapmain 2 alt? my $newmainid= $swapbest{$ids}{$src}; NO good; try bid ??
      
      my $altmaperr=""; # = locusaltmaperr($ids,$src); # do even if maperr?
      if($bestGid) { 
        # bestGid is always 't1' but may be swapmain w/ alt; ids should be swap of
        
        #No-effect# if($path =~ /swap1:/ and $maplocswap{$bid}{$src}) { $altmaperr= locusaltmaperr($ids,$src, $bestGid, $bestGidsrc,$bid); } else
        #no# if($ids ne $bid) { $altmaperr= locusaltmaperr($bid,$src, $bestGid, $bestGidsrc);  } else
        #No-effect# if($swapbest{$ids}{$src} eq $bestGid) { $altmaperr= locusaltmaperr($ids,$src, $bestGid, $bestGidsrc,$bid); } else
        #No-effect# if($swapbest{$bestGid}{$bestGidsrc} eq $ids) { $altmaperr= locusaltmaperr($ids,$src, $bestGid, $bestGidsrc,$bid); } else
        
        $altmaperr= locusaltmaperr($ids,$src, $bestGid, $bestGidsrc); 
          ## try to catch swap1 ERRrev misses.. nothing else is working
        if(not $altmaperr and $path =~ /swap1:/) { 
          #no# my $bvals= $cds2aa{$bestGid}{$bestGidsrc};
          #no# my $bmap= $$bvals[cMAPLOC];
          my $bmap= $maplocswap{$bestGid}{$bestGidsrc};
          my ($tor,$bor)= map{ (m/:([+\.-])/)?$1:0 } ($maploc,$bmap);
          $altmaperr.="ERRrev," if($tor and $bor and $tor ne $bor);
        }

      } else { $altmaperr= locusaltmaperr($ids,$src); }
      
      if($altmaperr and ($debug or not $maperr)) { $maperr.= ($maperr)?",":"err:"; $maperr .= $altmaperr; }
      
      my $covpi= $cov * $pid/100;
      #o# my $nsplit= ($path eq "0") ? 0 : (m/^(\d+)/) ? $1 : 1; #?
      my $nsplit= ($path =~ /^0/) ? 0 : (m/^(\d+)/) ? $1 : 1; #?
      
      ## better compare by rounding hobits, hodiff, aadiff, covpi ?
      $hobits= vround($hobits,8);
      $hodiff= vround($hodiff,8);
      $aadiff= vround($aadiff,8);
      $covpi = vround($covpi,4);
      
      my $errv=0;

      if(my $upd= $updates{$ids}) { #  || $updates{$gid}
        if($upd=~m/^bestsrc=(\w+)/) {  ## other swapbest=id1/id2 ??
          my $bs=$1; $errv= ($bs eq $src) ? 199 : -1;
        }
      }
      
      if($maperr) { 
        $errv= -9 unless($errv>0); # rank err types? eg gmap2:Scaf-change vs gmap3:C1_split
        ## $hobits= "-999.$hobits" ; #BAD try other way  # sort last .. BUT other problems? other map errs?
        $maperr=~s/^,//; $maperr=~s/^err://;
        $maploc=~s/,err.\w*//; $maploc.=",err:$maperr";
        $$vals[cMAPLOC]= $maploc;
        }
      
      if($errv==0 and $src =~ /augmod/) { $errv=999; } # always best !??
        ## augmod has n=175 err:Scaf_change .. only this maperr type, cant ignore
                
      if($errv==0 and $swapbest{$ids}{$src}) { $nswapb++; $errv=99; } #was $hobits=999999;
      
      $score{$src}= [$hobits,$hodiff,$aadiff,$covpi,$nsplit,$src,$errv,]; # sort for best
      $sval{$src}= $vals;  
      # record maperrs ? need for bestsrc: cancel if maperr (and isalt)?
    }
    
    my @sbest= sort _sortgscore values(%score);
    my ($hobits,$hodiff,$aadiff,$covpi,$nsplit,$src,$errv)= @{$sbest[0]};
    my $bestsrc= $src; # == $sbest[0]->[5];
    my $othersrc= join",", map{ $_->[5] }@sbest;
    my $bestval= $sval{$bestsrc};
    my @outrow= @{$bestval}[@OUTCOLS];  
    
    ## need to flag output table for bestsrc with maperr, presumably any others also have maperr..
    ## add to BestMap src col=1 or maploc col=12
    ## outrow[xx] == maploc  has ,err:xxx now
    
    my $bid= $bestval->[0]; $bid=~s/_C\d+$//; 
    if($ids ne $bid) { #?? was:  and not $didput{$bid} .. need these for swapmain maploc tests
      # ^^ FIXME2, need to track these altid > mainid changes, and add reverse swapid (t1 > talt),
      my $oswap= $cds2aaswap{$ids}{$bestsrc};
      my $omap= $maplocswap{$ids}{$bestsrc};
      my $oval= $cds2aa{$bid}{$bestsrc};
      if($oswap) { # or $oval
        # $cds2aaswap{$ids}{$bestsrc} == val for bestval->[0] id ?
        my($obid)= $oswap->[0]; # presume obid == ids
        if($oswap) { ##? was:  and not $didput{$obid}
          # insert b2 vals for this, flag swap
          $swapbest{$bid}{$bestsrc}=$ids; $nswap++;          
          $cds2aa{$bid}{$bestsrc}= $oswap; # swap AQuery ids .. but what else?
          #no# $maploc{$bid}{$bestsrc}= $omap if($omap); ## bad, sets err:ERRspan to many alts, swapmain locusaltmaperr bug?
        }
      }
    }

      #?? is this right or not?    
    if(my $idold= $didput{$bid}) {
      # warn "#malt.out: didput $bid as $idold, not $ids\n" if($debug);
      $nskipdup++; next;
    }
    
    if($ids =~ /t1$/) {
      ($bestGid, $bestGidsrc)= ($bid,$bestsrc);
      #? $maploc{$bestGid}{$bestGidsrc}= $bestmaploc; # same src or best src here??
    }
    
    ## ?? option to output scoring for all sources? esp maperr, hoscore, 
    ## add out col aasource == pubaa if hodiff,aadiff are low (<95% hodiff?)

    my $Selcstop=0;
    if(my $pubv=$aaq{$bid}{$SPUB}) {  ## NOT ids here, bid
      my($d,$w,$g,$aq,$cl,$off)=@$pubv; 
      $Selcstop=1 if($aq=~m/,selc/i);  # Funhe2EKm010847t1	556	0	556,80%,complete,selcstop
    } 

    my $aasrc=$bestsrc;
    ## these are now vround() scores .. need unround here?
    ## .. main effect of vround() reduces SPUB calls, good?
        ## FIXME: bestaa=pubaa for Selcstop ***
    
    if($hodiff == 0) { if($aadiff < $AAabsDIFFMIN) { $aasrc=$SPUB; } }
    elsif($hodiff < $HODIFFMIN) { $aasrc=$SPUB; }
    if($Selcstop) { $aasrc=$SPUB; } # always? .. has exception for *==Selc

    splice(@outrow,2,0,$aasrc); # BestAa     
    unshift(@outrow,$ids);   # LocusID ?? vs AQueryID
    push(@outrow,$othersrc); # othersrc
    
    $didput{$bid}=$ids;
    print $outh join("\t",  @outrow)."\n"; 
  } # id
  warn "#malt.out: nswap=$nswap,$nswapb; nskipdup=$nskipdup\n" if($debug);
}



#======= OLD ================================

my(%mapq,%mapv);
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
    if($cdq < $AAPercentDIFFMIN) { $bsrc=$SPUB; }       
    my $p=""; $p="#e." if($mapv =~ /ERR/); #skip this one??
    print $outh join("\t",$p.$td,"$td.$bsrc","noho",$cdv,$mapv)."\n"; 
    $did=1; # $did{$td}++;
    }
  return $did;
}

# from /bio-grid/kfish2/submitf/evgkfish2subgff.pl
sub makeAltBestTab {
  my($outf,$outh)= @_;


# GeneXrefStrandProblem : alts rev or to main/t1 .. drop here ???
# -- need to check these some how, drop alts on rev strand, else call as new gene locus if evidence/homol good?

  unless($outh and ref $outh) {
    $outf ||= "kfish2rae5h_altbest.tab"; #? $outh= *STDOUT;
    warn "#malt.out: $outf\n" if($debug);
    rename($outf,"$outf.old") if(-s $outf);
    open($outh,'>',$outf) or die $outf;
    }
  
  #o# print $outh join("\t",qw(AQueryID BestRepID HoDiff clAAD cov pid GenomeID:gespan path))."\n";  
  print $outh join("\t",qw(AQueryID BestRepID HoDiff AaDiff cov pid GenomeID:gespan path))."\n";  
  
  ## drop inaadiff tables ?? see above inAaQual for cdsq,cdsv
  # for my $inf (@inaadiff) {
  #   warn "#malt.in: $inf\n" if($debug);
  #   open(F,$inf) or die $inf;
  #   while(<F>) {
  #     next if(/^\W/ or /^AQueryID/); chomp; my @v=split"\t";
  #     my($aacl,$td,$cdsq,$src)=@v; 
  #     my($cdsp)= $cdsq=~/^(\d+)/; $cdsp||=0; $cdsp=100 if($cdsp>100);
  #     $aacl=~s/DIFF.0/DIFF0/; $aacl=~s/DIFF..*/DIFF1/; $aacl=~s/\..*$//;
  #     $cdsq{$td}{$src}=$cdsp; $cdsv{$td}{$src}= $aacl; 
  #   } close(F);
  # }

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
        $err ||= "noer";
        $mapv{$td}{$src}= join("\t",$cov,$pid,$gloc,$path,$err);
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

#--------

## add source tag col: gmap2pub5h ..
#i ==> kfish2nsplign15h_fcds4.cds2aatab <== #i Funhe2EKm000014t1       gspl3n5h
#i ==> kfish2nsplign15n_fcds4.cds2aatab <== #i Funhe2EKm000004t1       gspl3n5n
#i ==> kfish2rae5h_asm3n_fcds4.cds2aatab <== #i Funhe2EKm000003t1       gmap3npub5h
#i ==> kfish2rae5h_kfish2a_t1fcds4.cds2aatab <== #i Funhe2EKm027243l2t1     gmap2pub5h
#i ==> kfish2x11tsa_pubid_asm3n_fcds4.cds2aatab <== #i Funhe2EKm003823t1       gmap3nx11tsa
## add cov,pid map stats cols
## handle errs: ocds=cancel:NEn,aatrans-failed, ie no cds2aa; 
## now get cds2alen = 0, but only 1 has protein= byhand, so most are missing from output
#.. aatrans-failed:
# kfish2rae5h_asm3n_fcds4.gff.gz:9018
# kfish2rae5h_kfish2a_t1fcds4.gff.gz:2268
# kfish2rae5h_kfish2a_tafcds4.gff.gz:6712
# kfish2x11tsa_asm3n_fcds4.gff.gz:19725
##
## drop field= tags for header col?
## Funhe2EKm008403l2t1  0.da    aalen=78,10%,complete   aaold=78,10%,complete   cdsoff=546-782  offp=0  clen=0  cxlen=234/2332  0       gmap2pub5h
## AQueryID ADiff cds2aalen  pubaalen  cdsoff  puboff  pubclen  cxlen cov pid path maptag
#.............



=item try1 out
    cat gmapn/findcds4upd/kfish2rae5h_asm3n_fcds4.besttab345 | egrep -v 't[2-9]      0' | head
  AQueryID              hodiff  hoval                                   aadiff  cds2alen        cdsoff  cov     pid     maploc  path    maptag
  Funhe2EKm000003t1     100     89%,358,platyfish:ENSXMAP00000008562    0       214,94%,complete        19-663  100     100     KN805525.1:5067-10366:- 0       gmap3npub5h
  Funhe2EKm000004t1     0       0                                       -134    969,58%,complete        650-3559        94      97      KN805525.1:25802-58390:+        0       gmap3npub5h
  Funhe2EKm000005t1     102     97%,1834,azmolly:XP_007548349           2       984,70%,complete        4-2958  100     92      KN805525.1:62577-77945:+        0       gmap3npub5h
  Funhe2EKm000006t1     100     100%,1979,azmolly:XP_007548361          0       1257,85%,complete       163-3936        100     100     KN805525.1:85415-96983:-        0       gmap3npub5h
  Funhe2EKm000007t1     98      100%,513,azmolly:XP_007548366           0       325,34%,complete        801-1778        100     96      KN805525.1:97595-107565:+       0       gmap3npub5h
  Funhe2EKm000008t1     100     89%,610,mayzebr:XP_004571788.1          0       543,83%,complete        118-1749        100     100     KN805525.1:107343-127504:-      0       gmap3npub5h
  Funhe2EKm000009t1     100     29%,68,zebrafish3:ENSDARP00000098613    0       228,99%,partial 1-684   100     100     KN805525.1:121590-124252:-      0       gmap3npub5h
  Funhe2EKm000010t1     100     100%,478,azmolly:XP_007548367           0       285,62%,complete        24-881  99      100     KN805525.1:131954-136869:-      0       gmap3npub5h

=item try2 out : all ingff, best out
  ?? add alts in rae5htaasm2_fcds4h

  env debug=1 ./makebesttab345.pl \
    gmapn/findcds4upd/kfish2{nsplign15n_fcds4h,nsplign15h_fcds4,rae5h_asm3n_fcds4,rae5ht1asm2_fcds4h}.gff.gz \
    gmapn/findcds4upd/fixset/kfish2*.fcds5.gff.gz \
    | less
  
#malt.in: pubgenes/kfish2rae5h.main.pub.aa.qual
#malt.in: pubgenes/kfish2rae5h.alt.pub.aa.qual
#malt.in: gmapnbl/aaeval2/kfish2rae5h_bestfcds45pubaa-fish4ref.bltop4
#malt.in: gmapn/findcds4upd/kfish2nsplign15n_fcds4h.gff.gz
#malt.in: gmapn/findcds4upd/kfish2nsplign15h_fcds4.gff.gz
#malt.in: gmapn/findcds4upd/kfish2rae5h_asm3n_fcds4.gff.gz
#malt.in: gmapn/findcds4upd/kfish2rae5ht1asm2_fcds4h.gff.gz
#malt.out: cds2aaBestOut
  out .. skipping alts
AQueryID           hodiff  hoval                           aadiff  cds2alen           cdsoff  cov    pid  maploc               path    maptag       othersrc
Funhe2EKm000003t1  100     89%,358,platyfish:ENSXMAP08562  0       214,94%,complete   19-663   100   100  KN805525.1:5067-10366:-  0   gmap3npub5h  gmap3npub5h,gmap2pub5u
Funhe2EKm000004t1  99      99%,1882,azmolly:XP_007548354   -11     1092,67%,complete  237-3515 91    0    Scaffold0:25799-58393:+  0   gmap2pub5u   gmap2pub5u,gspl3n5n,gmap3npub5h
Funhe2EKm000004t2  0       0                               -51     969,67%,complete   357-3266 90    95   KN805525.1:39579-58049:+ 0   gmap3npub5h  gmap3npub5h
Funhe2EKm000005t1  102     97%,1834,azmolly:XP_007548349   2       984,70%,complete   4-2958   100   92   KN805525.1:62577-77945:+ 0   gmap3npub5h  gmap3npub5h,gmap2pub5u
Funhe2EKm000006t1  100     100%,1979,azmolly:XP_007548361  0       1257,85%,complete  163-3936 100   100  KN805525.1:85415-96983:- 0   gmap3npub5h  gmap3npub5h,gmap2pub5u
Funhe2EKm000007t1  98      100%,513,azmolly:XP_007548366   0       325,34%,complete   801-1778 100   96   KN805525.1:97595-107565:+  0 gmap3npub5h  gmap3npub5h,gspl3n5n,gmap2pub5u
Funhe2EKm000008t1  100     89%,610,mayzebr:XP_004571788.1  0       543,83%,complete   118-1749 100   100  KN805525.1:107343-127504:- 0 gmap3npub5h  gmap3npub5h,gmap2pub5u
Funhe2EKm000009t1  100     29%,68,zebrafish3:ENSDARP098613 0       228,99%,partial    1-684    100   100  KN805525.1:121590-124252:- 0 gmap3npub5h  gmap3npub5h,gmap2pub5u
Funhe2EKm000010t1  100     100%,478,azmolly:XP_007548367   0       285,62%,complete   24-881   99    100  KN805525.1:131954-136869:- 0 gmap3npub5h  gmap3npub5h,gmap2pub5u
Funhe2EKm000011t1  100     100%,1201,azmolly:XP_007546764  -22     661,97%,complete   31-2016  98    100  KN805568.1:402782-424715:- 0 gmap3npub5h  gmap3npub5h,gspl3n5n,gmap2pub5u
Funhe2EKm000012t1  100     100%,1813,platyfish:ENSXMAP08696 0      889,69%,complete   95-2764  99    98   KN805525.1:143911-166069:- 0 gmap3npub5h  gmap3npub5h,gspl3n5n,gmap2pub5u
Funhe2EKm000013t1  100     100%,882,platyfish:ENSXMAP08710  0      427,73%,complete   15-1298  100   99   KN805525.1:215256-220716:+ 0 gmap3npub5h  gmap3npub5h,gmap2pub5u

Funhe2EKm000014t1       100     97%,474,platyfish:ENSXMAP00000008715    -56     251,64%,complete        62-817  69      98      KN805525.1:221521-229295:-      0       gmap3npub5h     gmap3npub5h,gspl3n5h,gspl3n5n,gmap2pub5u
> Funhe2EKm000015t1 missing, where at?
Funhe2EKm000016t1       0       0       0       203,28%,complete        96-707  100     100     KN805525.1:245257-248644:-      0       gmap3npub5h     gmap3npub5h,gmap2pub5u
Funhe2EKm000017t1       100     100%,374,platyfish:ENSXMAP00000008730   0       187,61%,complete        123-686 100     99      KN805525.1:250195-255518:-      0       gmap3npub5h     gmap3npub5h,gmap2pub5u
Funhe2EKm000018t1       131     32%,64,zebrafish3:ENSDARP00000126870    -46     193,65%,complete        69-650  100     86      KN805525.1:258087-264147:-      0       gmap3npub5h     gmap3npub5h,gmap2pub5u,gspl3n5n
Funhe2EKm000019t1       98      94%,841,azmolly:XP_007548374    -28     568,70%,complete        142-1848        85      0       Scaffold0:267346-281089:+       0       gmap2pub5u      gmap2pub5u,gspl3n5h,gmap3npub5h
Funhe2EKm000019t3       100     98%,908,azmolly:XP_007548373    0       599,61%,complete        96-1895 100     99      KN805525.1:267309-281534:+      0       gmap3npub5h     gmap3npub5h
               ^^ ?? pick 19t3 over t1
Funhe2EKm000020t1       0       0       0       35,39%,partial5 3-110   0       0       Scaffold0:281953-289713:+       0       gmap2pub5u      gmap2pub5u,gmap3npub5h
Funhe2EKm000021t1       100     92%,889,platyfish:ENSXMAP00000008745    3       602,82%,complete        5-1813  100     89      KN805525.1:297161-306776:+      0       gmap3npub5h     gmap3npub5h,gspl3n5n,gmap2pub5u
Funhe2EKm000022t1       100     18%,107,mayzebr:XP_004574478.1  0       259,60%,complete        167-946 100     100     KN805525.1:339016-371906:+      0       gmap3npub5h     gmap3npub5h,gmap2pub5u
Funhe2EKm000025t1       52      51%,590,azmolly:XP_007548382    -329    287,45%,complete        97-960  72      86      KN805525.1:419179-427765:+      0       gmap3npub5h     gmap3npub5h,gspl3n5h,gspl3n5n,gmap2pub5u
Funhe2EKm000025t2       100     57%,655,azmolly:XP_007548382    -22     319,100%,partial        1-957   98      88      KN805525.1:419176-426811:+      0       gmap3npub5h     gmap3npub5h
               ^^ ?? pick 25t2 over t1

=item try3 out : all ingff, best out

  env  debug=1 ./makebesttab345.pl \
    gmapn/findcds4upd/{kfish2nsplign15n_fcds4h,kfish2nsplign15h_fcds4,kfish2rae5h_asm3n_fcds4,kfish2rae5ht[1a]asm2_fcds4h}.gff.gz \
    gmapn/findcds4upd/fixset/kfish2*.fcds5.gff.gz \
      > gmapn/findcds4upd/kfish2all5_fcds45.besttab2
    n=101618

  * fixme: alt-best-hoscore to replace t1main, if other map quals good.. see above 0019t1, 0025t1
  
try3a best for hoscore cases
  24537 gmap3npub5h, 8989 gmap2pub5u, 1100 gspl3n5h, 440 gspl3n5n
      ** ?? fx missing, some have best ho; bad HOTAG ??
      
try3b best cds2aa sources for hoscore (better pubaa not a class)
  24222 gmap3n5h,  8476 gmap2a5u, 1111 gspl3n5h  
    657 gspl3n5hfx, 444 gspl3n5n, 238 gspl3n5nfx
     total hoscore n=35149 (includes 7880 alts, 27269 t1main)
     by sort order gmap2 is best for ties. change that to gmap3n5h? how many ties?
     t1-no-hoscore n=3606, talt n=63083
 best cds2aa sources for no-hoscore:
   2372 gmap3n5h, 475 gmap2a5u, 467 gspl3n5h, 166 gspl3n5n
     71 none,  39 gspl3n5hfx,  16 gspl3n5nfx

  #malt.in: pubgenes/kfish2rae5h.main.pub.aa.qual
  #malt.in: pubgenes/kfish2rae5h.alt.pub.aa.qual
  #malt.in: gmapnbl/aaeval2/kfish2rae5h_bestfcds45pubaa-fish4ref.bltop4
  #malt.in: gmapn/findcds4upd/kfish2nsplign15n_fcds4h.gff.gz
  #malt.in: gmapn/findcds4upd/kfish2nsplign15h_fcds4.gff.gz
  #malt.in: gmapn/findcds4upd/kfish2rae5h_asm3n_fcds4.gff.gz
  #malt.in: gmapn/findcds4upd/kfish2rae5ht1asm2_fcds4h.gff.gz
  #malt.in: gmapn/findcds4upd/kfish2rae5htaasm2_fcds4h.gff.gz
  #malt.in: gmapn/findcds4upd/fixset/kfish2nsplign15h.fcds5.gff.gz
  #malt.in: gmapn/findcds4upd/fixset/kfish2nsplign15n.fcds5.gff.gz
  #malt.in: gmapn/findcds4upd/fixset/kfish2rae5h_gmap3n.fcds5.gff.gz
  #malt.out: cds2aaBestOut

 cat gmapn/findcds4upd/kfish2all5_fcds45.besttab2 | cut -f1-5,7,9-11 | egrep -v ' 0       0' | head 
  AQueryID      hodiff  hoval   aadiff  cds2alen        cov     maploc  path    maptag
  Funhe2EKm000003t1     100     89%,358,platyfish:ENSXMAP00000008562    0       214,94%,complete        100     KN805525.1:5067-10366:- 0       gmap3npub5h
  Funhe2EKm000004t1     99      99%,1882,azmolly:XP_007548354   -11     1092,67%,complete       91      Scaffold0:25799-58393:+ 0       gmap2pub5u
  Funhe2EKm000005t1     102     97%,1834,azmolly:XP_007548349   2       984,70%,complete        100     KN805525.1:62577-77945:+        0       gmap3npub5h
  Funhe2EKm000006t1     100     100%,1979,azmolly:XP_007548361  0       1257,85%,complete       100     KN805525.1:85415-96983:-        0       gmap3npub5h
  Funhe2EKm000007t1     98      100%,513,azmolly:XP_007548366   0       325,34%,complete        100     KN805525.1:97595-107565:+       0       gmap3npub5h
  Funhe2EKm000008t1     100     89%,610,mayzebr:XP_004571788.1  0       543,83%,complete        100     KN805525.1:107343-127504:-      0       gmap3npub5h
  Funhe2EKm000009t1     100     29%,68,zebrafish3:ENSDARP00000098613    0       228,99%,partial 100     KN805525.1:121590-124252:-      0       gmap3npub5h
  Funhe2EKm000010t1     100     100%,478,azmolly:XP_007548367   0       285,62%,complete        99      KN805525.1:131954-136869:-      0       gmap3npub5h
  Funhe2EKm000011t1     100     100%,1201,azmolly:XP_007546764  -22     661,97%,complete        98      KN805568.1:402782-424715:-      0       gmap3npub5h

=cut
