#!/usr/bin/env perl
# arp_condef.pl

=item notes

  make consensus def for gene groups
  cut from arp-blastpgp.pl
  
  inputs:
    aabugs3_omclgn.tab : table of ARPid, geneid
    *.aa.deflines.gz (just grep ^> from .aa )

  gunzip -c pro3/*.aa.deflines.gz | cat aabugs3_omclgn.tab  - \
   | perl arp_condef.pl > aabugs3_omclgn.consensus_def.txt
  
 Fish cleans:
  -- see also evigene2genbank for more name cleans per NCBI Genbank requirements
perl -pe 's/; src=.*//; 
s/\s+homolog//; 
s/H.sapiens //i; s/human and mouse\s*//i; s/human\s*//i; 
s/Novel\s*//; s/Putative\s+//i; 
s/\s*\([^\)\n]+\) *$//;  # many (Species) trailers in these names
s/\s*[\dA-Z]+$//;  # leave this trailing alphanum, cut for viewx
  
  FIXME? no: cut off src= from names. in arp_condef3. for omcl_subgroups call

=item upd 2017 notes
 
  $TEST_CONSENSUSNAME = $ENV{testcon} ||  $ENV{test} || 0; 
  -- 2017nov try again, no good protein_names:consensusname() ; rewrite or drop
  this get_condesc() needs work, esp 1 long name scores >> many consens short name; 
   .. add term/nspecies weight
  * should add CDD conserved domain names for consensus checks

  BUG in arp_condef3 parse of goodname/badname, requires species prefix in ID to match good/poorname spp
  from ($spp)=  id =~ m/^([a-zA-Z]+)/) ;  not same as protein_names: handling
  
=item 2017.nov daph10_omcl consname tests

  bio-grid/daphplx/rnasm/aaeval/omcldap/daph10_omcl

  myspecies=dpx17evg
  clade=DaphniaInsectFish
  speciesmap=bemtab=Whitefly,dapmaevg14=Daphnia_magna,dpx17evg=Daphnia_pulex,dapgal16tsa=Daphnia_galeata,dapsim17evg=Daphnia_similoides,drosmel16nc=Fruitfly,tribcas16nc=Beetle,guppync1=Guppy,zebrafish3=Zebrafish,human=Human

  # some dapmag names bad, skip? are dpx17 names ok? zebrafish?   no names for dapgal,dapsim
  ## * need to use ($spp)=  id =~ m/^([a-zA-Z]+)/) here, see arp_condef3
  ## need to remove crap from input ../names/*names (spp ids, etc)
  #try1: most human
  #try2:
  # goodname='guppy|tribcas|bemtab|drosmel|dpx17evg|human'
  # poorname='dapmaevg|zebrafish'
  #try3: leave dpx17evg,zebrafish unclassed
  # goodname='guppy:2,tribcas:2,bemtab:1,drosmel:2,human:1' << problem syntax for arp_condef3
  # poorname='dapmaevg'
  #try4: fix names: bemtab,dapmaevg14,dpx17evg,drosmel16nc,guppync1,human,tribcas16nc,zebrafish3
  # goodname='bemtab,guppync,drosmel,human,tribcas'
  # poorname='dapmaevg'
  
  # naming steps: output omclgn.consensus_def.txt, ugp.txt, ugp_brief.txt, omsubgrp1.names
  otopts="-steps=condef,genegroup,subgroup"
  
  # try8: env test=1 for arp_condef3.pl == TEST_CONSENSUSNAME of protnames.pm
  # .. testcon not better, making dumb nonconsensus mistakes
  #  my $cmd="$EVIGENES/omcl/arp_condef3.pl -noput -nodigit -nolike -gtag=$GTAG -idprefix=$IDPRE"
  #    . " -good='$goodname' -poor='$poorname' $inlist > $outfile";
  #x export test=1; export testcon=1
  
  # try9: remove good/poor names, use all equal, probably best for this sppname set
  goodname='NOGOODNAME'
  poorname='NOPOORNAME'

  env myspecies=$myspecies date=$date clade=$clade speciesmap="$speciesmap" \
    goodname="$goodname" poorname="$poorname" \
    $evigene/scripts/omcl/orthomcl_tabulate.pl $otopts -debug -idprefix=$idprefix -omclpath ./ -namepath ../names

  orthomcl_tabulate sub cmds:
  evigene/scripts/omcl/arp_condef3.pl -noput -nodigit -nolike \
    -gtag=DAPHa -idprefix=DAPHa_G -good='NOGOODNAME' -poor='NOPOORNAME' \
      ../names/bemtab16nc.names ../names/daphmag.names .. ../names/zebrafish2.names \
      daph10omcla_omclgn.tab > daph10omcla_omclgn.consensus_def.txt
  
  env xml=0 date=20171114 clade=DaphniaInsectFish gtag=DAPHa idprefix=DAPHa_G \
    speciesmap=bemtab=Whitefly,dapmaevg14=Daphnia_magna,..,human=Human \
    $evigene/scripts/omcl/genegroupbpo.pl \
    daph10omcla_omclgn2sum.tab  daph10omcla_omclgn.consensus_def.txt \
    ../names/*.names  ./all_orthomcl.out  > daph10omcla_genes.ugp.txt
  
  evigene/scripts/omcl/arp_condef3.pl .. ditto .. \
    daph10omcla_omsubgrp1.list > daph10omcla_omsubgrp1.names
    
=cut

use FindBin;
use lib ("$FindBin::Bin/..","$FindBin::Bin/../prot/"); # assume evigene/scripts/omcl/ << this path
my $EVIGENES=$ENV{EVIGENES} || "$FindBin::Bin/..";  

use strict;
use Getopt::Long;
use protein_names; # has nameclean() recase() isTEname()

my $TEST_CONSENSUSNAME = $ENV{testcon} ||  $ENV{test} || 0; 
  # 2017nov try again, bad consensus.. rewrite or drop
  
my $debug=0;
my $poorannots = $ENV{poorname} || 'POORSPECIES'; ##was  ENV{poor} 'amel|apis';
my $goodannots = $ENV{goodname} || 'GOODSPECIES'; ## was ENV{good} 'culex|aedes|ixodes|pediculus';
my $GTAG = $ENV{gtag} || "ARP";
my $OGPRE= $ENV{idprefix}||"${GTAG}1_G";
my $NOCLEAN= $ENV{noclean} || 0;
my $recase= $ENV{recase} || 0;

#n now in pack protein_names; our ...
 $NAME_NOPUTATIVE= $ENV{noput} || 0;;
 $NAME_NOLIKE= $ENV{nolike} || 0;  
 $NAME_NODIGITS= $ENV{nodigits} || 0;;
 $USE_TENAME=$ENV{USE_TENAME}||0; # 2017upd NOT default, was 1
## our $NAME_UNK  = "Uncharacterized protein"; # uniprot  

my $optok= GetOptions( 
  "idprefix=s", \$OGPRE, "GTAG=s", \$GTAG, # dont need both, make GTAG from IDPREFIX
  "goodannots=s", \$goodannots, "poorannots=s", \$poorannots, 
  "noclean", \$NOCLEAN,
  "nolike", \$NAME_NOLIKE,
  "noputative", \$NAME_NOPUTATIVE,
  "nodigits", \$NAME_NODIGITS,
  "recase!", \$recase,
  "debug!", \$debug,
);

my (%gdef, %ar);
sub MAIN_stub{}

# FIX: consensusname(): our($goodannots,$poorannots) not pack global yet.. also %sppwt=fixed set
# $goodannots =~ s/,/|/g; $poorannots =~ s/,/|/g;
map{ s/[,;\s]+/|/g; } ($goodannots,$poorannots);
$NAME_GOODSPP=$goodannots;
$NAME_POORSPP=$poorannots;

# sub readIntables() {} ..
while(<>) {
  if(/^>(\S+)/) { #?? require >
    my $gn=$1; my $an=$_; 
    $gn=~s/_/:/ unless($gn=~/:/); #DAMMIT.need.for.all.now# 
    if($an=~/\t/) { my @an=split"\t",$an; $an=join("\t",@an[0,1]); }
    $an =~ s/(MD5|length|loc)=[^\s;]+;?//g;  #? leave in drosmel ID=FBpp; parent=FBgn; ?
    $gdef{$gn}= $an;
    
  } elsif(/^$GTAG(\d+)\s+(\S+)/) { ## xxx_omclgn.tab
    my $ar=$1; my $gn=$2; # dang: species_gene here, species:gene names
    $gn=~s/_/:/ unless($gn=~/:/); #DAMMIT.need.for.all.now# 
    push( @{$ar{$ar}}, $gn); # $gg{$gn}=$ar; 

## upd 2013nov: parse, name subgroup tables:
# FISH11D998.s1:  catfish_IctpunEGm005221t1,human_UniRef50_Q9Y519,zfish_ENSDARP00000090847,zfish_ENSDARP00000118506,zfish_ENSDARP00000121734,
# FISH11D998.s2:  catfish_IctpunEGm031540t1,human_UniRef50_Q6ZMB5,kfish2_Funhe2EKm004356t1,medaka_ENSORLP00000005890,platyfish_ENSXMAP00000009340,spotgar_ENSGACP00000018640_1,stickleback_ENSGACP00000018640,tetraodon_ENSTNIP00000007115,zfish_ENSDARP00000055153,
# FISH11D998.s3:  kfish2_Funhe2EKm001955t1,mayzebr_XP_004557607.1,medaka_ENSORLP00000003534,stickleback_ENSGACP00000006021,tetraodon_ENSTNIP00000015909,tilapia_ENSONIP00000024767,
  } elsif(/^$GTAG([\w\.:]+)\s+(\S+)/) { 
    my $ar=$1; my $gns=$2; $ar=~s/:$//;
    my @gns=split",",$gns;
    for my $gn (@gns) {
    $gn=~s/_/:/ unless($gn=~/:/); #DAMMIT.need.for.all.now# 
    push( @{$ar{$ar}}, $gn);  
    }
    
  } elsif(/^(\w\S+)\s+/) { # more names; may be \t table, split off \textra ?
    my $gn=$1; my $an=$_;
    $gn=~s/_/:/ unless($gn=~/:/); #DAMMIT.need.for.all.now# 
    if($an=~/\t/) { my @an=split"\t",$an; $an=join("\t",@an[0,1]); }
    $an =~ s/(MD5|length|loc)=[^\s;]+;?//g;  #? leave in drosmel ID=FBpp; parent=FBgn; ?
    $gdef{$gn}= $an;
  }
}

if($TEST_CONSENSUSNAME) {

# FIX: consensusname(): our($goodannots,$poorannots) not pack global yet.. also %sppwt=fixed set
#above# $NAME_GOODSPP=$goodannots;
#above# $NAME_POORSPP=$poorannots;

foreach my $ad (sort{$a <=> $b or $a cmp $b} keys %ar) {
  ## consname input: list of tabbed "$de,$palign,$id1,$id2"

  ## DAMMIT ids here are spp_ not spp: !! what gives, fixed at input above.. Ahh gdef{id}=an an NOT fixed
  my @gdef= map { 
    my($gn,$de)= split" ",$gdef{$_}, 2; chomp($de); 
    $gn=~s/^>//; $gn=~s/_/:/ unless($gn=~/:/); #DAMMIT.need.for.all.now# 
    join("\t",$de,99,$gn); 
    } @{$ar{$ad}};
  
  my $bestnameid= consensusname(\@gdef); 
  my ($bna,$bpi,$bid)=split"\t",$bestnameid;
  my $nac= nameclean($bna);  ## NOTE: consname does nameclean() but returns orig name
  print $OGPRE.$ad,".consensus\t",$nac,"\tsrc=",$bid,"\n";
}

} else {
foreach my $ad (sort{$a <=> $b or $a cmp $b} keys %ar) {
  my @gdef= map{ $gdef{$_} } @{$ar{$ad}};
  my $condesc= get_condesc( $OGPRE.$ad, \@gdef); #  $OGPRE.$ad not used in condesc()?
  print $condesc,"\n";
}
}

=item condesc / consensus description

  FIXMEd: now uses  evigene/scripts/prot/protein_names.pm
  
  see also protein_names:consensusname() 
    derived from this but wants align scores, maybe better
    where input is list of ($name,$palign,$id1,$id2)
    return is best row of list, (name,palign,id1,id2), unchanged.
      palign = percent align (1 val) or /^(\d+)%,(\d+).(\d+),(\d+)/ == namealign string
        can be missing/empty.
      id2 == optional 2nd dbxref for de name;
  
=cut

sub get_condesc {
  my($cluid, $gdef)= @_;
  
  my (%didspp,%de,%src);
  my $noid=0;
  foreach (@$gdef) {
    ##nameclean: my ($isunk,$isput,$islike)= (0,0,0);
        
    next if(/by Gnomon/); #?? or not #  Gene predicted by Gnomon |  Partial gene predicted by Gnomon ..
    s/^>//; chomp;
    my($id)= (s/^(\S+)\s+//) ? $1 : "noid".++$noid; 
    $id=~s/_/:/ unless($id=~/:/); #DAMMIT.need.for.all.now# 
    my $spp= ($id =~ m/^([a-zA-Z]+)/) ? $1 : "NULLSPP"; # $id; #?? dont use id for default.. use one dummy?
    my $src= (s/\b(src[=:]\S+)//i) ? $1 : ""; #? dont need id becomes src=$bestid
    my $dbx= (s/\b(dbxref[=:]\S+)//i) ? $1 : "";
    if(s/(Name|desc)[=:]//) { ; } ##old: s/(\w+)[=:](.+)$//;  # extra.. gone? desc= old Name= key; 
    ## nameclean() looks for @NAMETAG_CUT @NAMETAG_KEEP key=value things to remove .. add more?
 
    #What bug ?? Name="; Uncharacterized protein;"
    unless($NOCLEAN) {
      $_ = nameclean($_) if(/\S/); # FIXME: have similar sub in 3 places, use evigene/prot/protein_names.pm ?
      s/\-like\b// if($NAME_NOLIKE); # nameclean fails??????? check for bug, our $NAME_NOLIKE fail?
      ## now in nameclean $NAME_NODIGITS; drop trailing nums: family form num:  protein 1|2|3... ; or user option? can be useful.
    } 

    s/^\s+//; s/\s+$//; s/  +/ /g; s/;\s*;/;/g; # done in nameclean() do again?
    $de{$id}= $_ if(m/\w\w/); # and not $didspp{$spp} #?? change: keep all/spp but downweight 2+ below by 0.3?
    $didspp{$spp}++;
  }
  
  # FIXME: option for subgroup names: pick 1 best spp among all of subgroup/grp names
  # fixme: many word name is overweighted by word count; add kwspp count over species per word
    
  #now what? need one readable/valid line from these; not all agree
  my (%kw, %kwd, %kwcomm, %kwspp, %wtd);

if(1) { # new
  my $NAME_NONE_PLUS = $NAME_NONE.'|putative|like|protein';
  foreach my $id (sort keys %de) {
    my $spp= ($id =~ m/^([a-zA-Z]+)/) ? $1 : "NULLSPP";  
    my @w= split /\W+/, lc($de{$id}); # words
    my $wt= 1; # $palign{$id} || 1; # palign or balign ?
    my %wseen=(); foreach (@w) { 
      unless($wseen{$_}++ or m/^($NAME_NONE_PLUS)/i) { 
        $kwspp{$_}{$spp}++; $kwcomm{$_}++; $kw{$_}+=$wt; $kwd{$id}{$_}+=$wt; }
    }
  }
} else { # old
  foreach my $id (sort keys %de) {
    my $w= $de{$id};  # $w =~ s/[,;=].*$//;
    my @w= split /\W+/, $w;
    ## $kw{$w} += @w; $kwd{$id}{$w} += @w; # score whole phrase? 
    foreach (@w) { 
      $_= lc $_; $kw{$_}++; $kwd{$id}{$_}++; 
      }
  }
}

  my @spp= sort keys %didspp; ## this maybe IDs if no spp tag found .. skip?
  my $nspp= @spp; $nspp||=1; $nspp=99 if($nspp>99);
  my %sppord=(); my %spo; 
  my $isp=$nspp; 
  map{ $sppord{$_}=$isp-- if(m/$goodannots/); } split/\W/,$goodannots;
  map{ $sppord{$_}=$isp-- unless($sppord{$_} or m/$poorannots/); } @spp;
  map{ $sppord{$_}=$isp-- if(m/$poorannots/); } split/\W/,$poorannots;
  
  %didspp=(); my $bestcomm=0;
  my @wt= sort{ $kw{$b} <=> $kw{$a} } keys %kw;
  foreach my $id (sort keys %kwd) {
    my $spp= ($id =~ m/^([a-zA-Z]+)/) ? $1 : "NULLSPP"; # $id;
    my $wt=0; my $wi= @wt;
    my $nw= int( 0.7 * $wi);
    $wi *= 3;

    ## new...  common isn't necessarily best, when e.g. lots of ensembl-zfish-named fishes from 1 source
    ## vs uniprot names from a few other sources
#     my @kw= sort keys %{$kwd{$id}};
#     my ($tuniq,$tcomm)=(0,0);
#     foreach my $kw (@kw) { if($kwcomm{$kw} > 1) { $tcomm++; } else { $tuniq++; } }
#     $wt += $tcomm - $tuniq; # ??
#     my $isgood= ($id =~ /$goodannots/)?2:($id =~ /$poorannots/)?0:1; 
#     my $bok=0;
#     if($tcomm > 0 and $tuniq == 0) { $bok=1; }
#     elsif($isgood > 1 and $tcomm > $tuniq) { $bok=1; }
#     elsif($tcomm - $tuniq > $bestcomm) { $bok=1; }
#    #.. more ..
#    #....
    
    foreach my $w (@wt[0..$nw]) { 
      #o# $wt += $wi if( $kwd{$id}{$w} ); 
      if( $kwd{$id}{$w} ) {
        $wt += $wi ;
        my $nsp= scalar(keys %{$kwspp{$w}}); #17upd
        $wt += 3 * $nsp if($nsp>1);      
      }
      $wi -= 3;
      }
    my $spo= $sppord{$spp}||0;
    # $wt = $wt * 0.4 if($id =~ /$poorannots/); # down/up weight others?
    # $wt = $wt * 1.2 if($id =~ /$goodannots/); # down/up weight others? : add more good annots
    $wt = $wt * ($spo/$nspp) if($id =~ /$poorannots/);  
    $wt = $wt * (1+$spo/$nspp) if($id =~ /$goodannots/); 
    $wt = $wt * 0.3 if($didspp{$spp}++);  
    $wtd{$id}= $wt;
    $spo{$id}=$spo;
  }
  
  # should require at least 2 cases of agreement, otherwise decline consensus desc
  my($bestid)= sort{ $wtd{$b} <=> $wtd{$a} or $spo{$b} <=> $spo{$a} or $a cmp $b } keys %wtd;
  my $de= $de{$bestid};
  $de =~ s/;\s*$//; $de =~ s/^[;\s]+//; # s/^\W+//; ???
  if($de=~/\w\w/){ $de.="; src=$bestid" if($bestid); } else { $de=$NAME_UNK; }
  
  my $clukey="$cluid.consensus";
  return (wantarray) ? ($clukey,$de) : $clukey."\t".$de;
}


