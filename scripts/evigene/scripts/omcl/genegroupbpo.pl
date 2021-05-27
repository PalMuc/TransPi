#!/usr/bin/env perl
# genegroupbpo.pl

=item usage

  gunzip -c aabugs_omclgns2.tab.gz pro2/*.aa.deflines.gz \
   omclfilt1_Aug_10/all_orthomcl.out.gz |\
    perl genegroupbpo.pl \
    >  arthropod_ogenesum4.ugp.xml

  #  gunzip -c  omclfilt1_Aug_10/all_orthomcl.out.gz | \
  #  cat aabugs_omclgn2sum.tab  proteomes/*.aa.deflines - | perl genegroupbpo.pl \
  #  > arthropod_ogenesum2.ugp.xml 

microbe% ll omclfilt1_Aug_10/all_orthomcl.out.gz 
  931529 10 Aug 19:08 omclfilt1_Aug_10/all_orthomcl.out.gz

microbe% ll proteomes/*.aa.deflines
   826371 14 Aug 12:32 proteomes/aedes.aa.deflines
 1857526 17 Aug 14:43 proteomes/aphid.aa.deflines
 1309004 17 Aug 15:16 proteomes/apis.aa.deflines
 1678628 14 Aug 12:32 proteomes/culex.aa.deflines
 1493906 14 Aug 12:32 proteomes/dmel.aa.deflines
 1433190 17 Aug 14:32 proteomes/nasonia.aa.deflines
 1557879 17 Aug 15:19 proteomes/tribolium.aa.deflines


=item aabugs_omclgns2.tab:

 gunzip -c aabugs_filtered.bpo.gz | \
  cat  aabugs_omclgn.tab  dmel_protid_keep.tab - \
 | perl pastebpo.pl > aabugs_omclgns2.tab &


=item notes
  
  * add group consensus descriptions from arp-blastpgp.pl
      arthropod_ogene3.consensus_def.txt
ARP1_G11.consensus conserved hypothetical protein; src=culex_CPIJ005972
ARP1_G12.consensus zinc finger protein 91; src=nasonia_NCBI_hmm386594
ARP1_G13.consensus
ARP1_G16.consensus Cytochrome P450 CYPm3r10 (Fragment).; src=anopheles_AGAP008213-PA
      
  * DONE: add in prot descriptions: fasta deflines if only that
  
  * DONE: need occurrence: format so lucene search can easily find
      phylogenetic sets of on/off gene clusters
  
  * DONE: lucene searches: need results.xsl result table w/ field values
        including descrip, occur, 
=cut

use strict;
use Getopt::Long;
use constant VERSION => '2013.12.30'; # 11.15'; # 10.08'; # 

my $debug= $ENV{debug}||0;
my $doxml= $ENV{xml}||0;
my $date=  $ENV{date}||"20080824"; # today() ?
my $OGPRE= $ENV{idprefix}||"ARP1_G";
my $GTAG = $ENV{gtag} || ""; # "ARP"; # no default now see below
my $CLADE= $ENV{clade} || "Eukaryotes"; #was  "Arthropoda";
my $GENE_TITLE=$ENV{title} || "$CLADE gene group"; # was cluster
my $GENE_SOURCE="euGenes";
my $GENE_TYPE="gene_group"; # a SO type
my $SPPGN= $ENV{sppgn} || 0; # leave species_GeneID on accession
my $speciesmap= $ENV{speciesmap} || 0; # leave species_GeneID on accession

# make these GetOptions: doxml, date, OGPRE, TITLE, ...

my (%gdef,%gsum,%aalen,$nAAlen);
my (%sppmap, %sppshort); # NO defaults now, caller must set ..
my (%FISH11sppmap, %FISH11sppshort,%FISH12sppmap, %FISH12sppshort,
    %PLA9sppmap, %PLA9sppshort,
    %BUG1sppmap, %BUG1sppshort,%BUG3sppmap, %BUG3sppshort);

my $optok= GetOptions( 
  "speciesmap=s", \$speciesmap, 
  "date=s", \$date, "title=s", \$GENE_TITLE, 
  "IDPREFIX=s", \$OGPRE, "GTAG=s", \$GTAG, # dont need both, make GTAG from IDPREFIX
  "CLADE=s", \$CLADE,  
  "xml!", \$doxml,
  "SPPGN!", \$SPPGN,
  "debug!", \$debug,
);

die "EvidentialGene genegroupbpo.pl VERSION ",VERSION,"
  orthomcl gene groups document
Usage: genegroupbpo  < clade_omclgn2sum.tab  clade_omclgn.consensus_def.txt \$names all_orthomcl.out > clade_ugp.txt
  opts: -xml -date 20131230 -idprefix FISH11_G -CLADE Fish  ..   
  -speciesmap 'kfish=Fundulus_heteroclitus,cacao=Theobroma_cacao,wasp=Nasonia_vitripennis'
" unless($optok); 

if($OGPRE and not $GTAG) { ($GTAG=$OGPRE)=~s/_[^_]*$//; } # or keep same?

MAIN();
#..........................

   
BEGIN{
%BUG1sppmap = (
  acyr1 => 'Acyrthosiphon pisum',
  aedes => 'Aedes aegypti',
  amel4 => 'Apis mellifera',
  anopheles => 'Anopheles gambiae',
  culex => 'Culex pipiens',
  daphnia => 'Daphnia pulex',
  dmel => 'Drosophila melanogaster',
  dmoj => 'Drosophila mojavensis',
  dpse => 'Drosophila pseudoobscura',
  ixodes => 'Ixodes scapularis',
  nasonia => 'Nasonia vitripennis',
  pediculus => 'Pediculus humanus',
  tcas3 => 'Tribolium castaneum',  
# arp3
  aphid => 'Acyrthosiphon pisum',
  apis => 'Apis mellifera',
  drosmel => 'Drosophila melanogaster',
  drosmoj => 'Drosophila mojavensis',
  drospse => 'Drosophila pseudoobscura',
  tribolium => 'Tribolium castaneum',  
  bombyx => 'Bombyx mori',
);

%BUG1sppshort = (
  acyr1 => 'Aphid',
  aedes => 'Aedes',
  amel4 => 'Apis',
  anopheles => 'Anopheles',
  culex => 'Culex',
  daphnia => 'Daphnia',
  dmel => 'DrosMel',
  dmoj => 'DrosMoj',
  dpse => 'DrosPse',
  ixodes => 'Ixodes',
  nasonia => 'Nasonia',
  pediculus => 'Pediculus',
  tcas3 => 'Tribolium',  
# arp3
  aphid => 'Aphid',
  apis => 'Apis',
  drosmel => 'DrosMel',
  drosmoj => 'DrosMoj',
  drospse => 'DrosPse',
  tribolium => 'Tribolium',  
  bombyx => 'Bombyx',
);

# clade=HymenoInsectArp
# spl="antc,anth,aphid,apis2ref,bombusimp,bombusterr,daphnia,drosmel,human,trica,wasp"
# Camponotus x;  Harpegnathos x

%BUG3sppmap= (
  aphid => 'Acyrthosiphon pisum',
  antc => 'Camponotus ant',
  anth => 'Harpegnathos ant',
  apis => 'Apis mellifera',
  apis2ref => 'Apis mellifera',
  daphnia => 'Daphnia pulex',
  drosmel => 'Drosophila melanogaster',
  wasp => 'Nasonia vitripennis',
  trica => 'Tribolium castaneum',
  tribol => 'Tribolium castaneum',
  bombusimp => 'Bombus imp.',
  bombusterr => 'Bombus terr.',
  human => 'Homo sapiens',
);

%BUG3sppshort = (
  antc => 'AntCamp',
  anth => 'AntHarp',
  daphnia => 'Daphnia',
  wasp => 'Nasonia',
  aphid => 'Aphid',
  apis => 'Apis',
  apis2ref => 'Apis',
  drosmel => 'DrosMel',
  trica => 'Tribolium',
  tribol => 'Tribolium',
  bombusimp => 'BombusI',
  bombusterr => 'BombusT',
  human => 'Human',
);

%PLA9sppmap = ( # FIXME, used for occurrance
# plant9
arath => 'Arabidopsis thaliana',
cacao => 'Theobroma cacao',
frave => 'Fragaria vesca',
poptr => 'Populus trichocarpa',
ricco => 'Ricinus communis',
soltu => 'Solanum tuberosum',  
sorbi => 'Sorghum bicolor',
soybn => 'Glycine max',
vitvi => 'Vitis vinifera',
);

%PLA9sppshort = (
# plant9
arath => 'Arabidopsis',
cacao => 'Cacao',
frave => 'Strawberry',
poptr => 'Poplar',
ricco => 'Castorbean',
soltu => 'Potato',  
sorbi => 'Sorghum',
soybn => 'Soybean',
vitvi => 'Grape',
);

my @FISH11spplist= split ",",
 "catfish,human,kfish2,mayzebr,medaka,platyfish,spotgar,stickleback,tilapia,tetraodon,zebrafish";
#set spplist="catfish human kfish2 mayzebr medaka platyfish spotgar stickleback tetraodon tilapia zfish"
%FISH11sppshort = map{ $_ => ucfirst($_) } @FISH11spplist;
  $FISH11sppshort{kfish2} = 'Killifish';
  $FISH11sppshort{mayzebr} = 'Maylandia';
  $FISH11sppshort{zfish} = 'Zebrafish';

%FISH11sppmap = (
  human => 'Homo sapiens',
  zfish  => 'Danio rerio', # old: zebrafish
  zebrafish  => 'Danio rerio', # old: zebrafish
  kfish2 => 'Fundulus heteroclitus', # old: killifish
  catfish => 'Ictalurus punctatus',
  mayzebr => 'Maylandia zebra',
  medaka => 'Oryzias latipes',
  platyfish => 'Xiphophorus maculatus',
  spotgar => 'Lepisosteus oculatus',
  stickleback => 'Gasterosteus aculeatus',
  tilapia => 'Oreochromis niloticus',
  tetraodon => 'Tetraodon nigroviridis',
  astymex => 'Astyanax mexicanus',
  #old# xenopus => 'Xenopus tropicalis',
);

# 2014/fish12 adds astymex=cavefish, spotgar3
# astymex catfish human   kfish2  mayzebr medaka  platyfish spotgar3 stickleback  tetraodon  tilapia zebrafish

my @FISH12spplist= split ",",
 "astymex,catfish,human,kfish2,mayzebr,medaka,platyfish,spotgar3,stickleback,tilapia,tetraodon,zebrafish";
%FISH12sppshort = map{ $_ => ucfirst($_) } @FISH12spplist;
  $FISH12sppshort{kfish2} = 'Killifish';
  $FISH12sppshort{astymex} = 'Cavefish';
  ## $FISH12sppshort{spotgar} = 'Spottedgar'; #?? why all ok but 'spotgar' ?
  $FISH12sppshort{spotgar3} = 'Spottedgar';
  $FISH12sppshort{mayzebr} = 'Maylandia';
  $FISH12sppshort{zebrafish} = 'Zebrafish';
  
%FISH12sppmap = %FISH11sppmap;
  delete $FISH12sppmap{zfish}; #  
  delete $FISH12sppmap{spotgar}; # below keys conflict w/ spotgar3 ..
  $FISH12sppmap{spotgar3}  = 'Lepisosteus oculatus';
  
}

## DAMIT: make spp:id same as omcl spp_id; other dingbats here?
sub fix_sppid {
  my $sppid=shift;
  my($spp1,$gn1)= split /[:_]/, $sppid, 2; $sppid= $spp1.'_'.$gn1 if($gn1); # or : ??
  return (wantarray) ? ($sppid,$spp1,$gn1) : $sppid;
}

# was  occurrence: Aedes:1, Anopheles:1, Aphid:0, Apis:0, Culex:1, Daphnia:0, DrosMel:6, DrosMoj:5, DrosPse:8, Ixodes:0, Nasonia:0, Pediculus:1, Tribolium:1;
# now  occurrence: Aedes=1 Anopheles=1 Aphid=0 Apis=0 Culex=1 Daphnia=0 DrosMel=6 DrosMoj=5 DrosPse=8 Ixodes=0 Nasonia=0 Pediculus=1 Tribolium=1;

sub occurrence {
  my($gna)= @_;
  # rename spp nicely: acyr => aphid, amel => apis, dmel => drosmel, ...
  my %spc= map{ my $sp=$sppshort{$_}|| $_; $sp => 0 } keys %sppmap; 
  foreach my $gn (@$gna) { $gn =~ m/^([^\W_]+)[:_]/;  my $sp= $sppshort{$1}|| $1; $spc{$sp}++;  } 
  my $occ= join " ", map{ "$_=$spc{$_}" } sort keys %spc;
  return $occ;
}

sub avescore {
  my($gna)= @_;
  my($sev,$spi,$ns)= (0)x3;
  foreach my $gn1 (@$gna) {
    my($sppgn,$spp,$gn)= fix_sppid($gn1); # my($sp1,$gn)= split /[:_]/, $gn1, 2;
    next unless (defined $gsum{$gn1} or defined $gsum{$sppgn});
    my($og,$ng,$ev,$pi,$g2,$e2,$p2) = (defined $gsum{$gn1}) ? @{ $gsum{$gn1} } :
      (defined $gsum{$sppgn}) ? @{ $gsum{$sppgn} }: ('.') x 10; # err if none?
    $sev += $ev; $spi += $pi;  $ns++;
  }
  if($ns>0) {
    $sev = sprintf "%.2g", $sev/$ns;  
    $spi = int($spi/$ns);
    }
  return($sev, $spi, $ns);
}

sub print_genesum {
  if(1 || $doxml) { print_ugpxml(@_); return; } # this works now for both .txt and .xml
  
  my($og2,$ng,$nt,$gna)= @_;
  my $ogid= $OGPRE.$og2;
  my $occ= occurrence($gna);
  my($sev,$spi,$ns)= avescore($gna);

  print 
"GeneSummary:
  title: $GENE_TITLE $ogid;
  source: $GENE_SOURCE;
  type:  $GENE_TYPE;
  basic_information:
    GeneID: $ogid;
    species: $CLADE;
    ntaxa: $nt;
    ngene: $ng;
    group-evalue: $sev;
    group-identity: $spi;
    occurrence: $occ;
    date:  $date;
  similar_genes:
";
  foreach my $gn1 (@$gna) {
    my($sppgn,$spp,$gn)= fix_sppid($gn1); # my($spp,$gn)= split /[:_]/, $gn1, 2;
    my $sna= $sppmap{$spp} || $spp;
    my($og,$ng,$ev,$pi,$g2,$e2,$p2) = (defined $gsum{$gn1}) ? @{ $gsum{$gn1} } :
      (defined $gsum{$sppgn}) ? @{ $gsum{$sppgn} }: ('.') x 10; # err if none?
    my $gdef= $gdef{$gn1} || $gdef{$sppgn} || "";
    $ev= sprintf "%.1g",$ev; # drop extra digits ?
    print "    similarity: ", join(", ",$sna,$gn,$ev,$pi,$gdef),";\n";
  }
  print "\n";
}


my $in= 0;

sub xtab { 
  my $t= ($in>0) ? ("  ") x $in : ""; print $t;
}

sub ptag {
  my($t,$v,$ln)= @_;
  xtab() unless($ln);
  print ($doxml? "<$t>":"$t: ");
  print ($doxml? "$v</$t>" : "$v;") if(defined $v); ## escape xml ***
  print (($ln) ? " " : "\n");
  #print "\n" unless($ln);
}

sub btag {
  my($t, $ln)= @_;
  xtab(); $in++;
  print ($doxml?"<$t>":"$t:");
  print (($ln) ? " " : "\n");
}

sub etag {
  my($t, $ln)= @_;
  $in--;
  xtab() unless($ln or !$doxml);
  print ($doxml?"</$t>\n":"");
  print "\n" if($ln and !$doxml); #??
}

sub print_ugpdummy {
btag "GeneSummary id=\"header\""; 
# ptag "title", "$GENE_TITLE $ogid";
ptag "source","$GENE_SOURCE";
# ptag "type","$GENE_TYPE";
# btag "basic_information"; 
# ptag "date",  $date;
# etag "basic_information"; 
etag "GeneSummary";  print "\n";
}

sub print_ugpxml {
  my($og1,$ng,$nt,$gna)= @_;

  my $ogid= $OGPRE.$og1;
  my $occ= occurrence($gna);
  my($sev,$spi,$ns)= avescore($gna);

  my $ogdef= $gdef{$ogid} || "";    
  $ogdef =~ s/[,\s]*src=(\S+)//i;
  $ogdef =~ s/desc=//;

  $in= 0;
btag "GeneSummary id=\"euGenes:$ogid\"";  # fixme: id= tag?
ptag "title", "$GENE_TITLE $ogid";
ptag "source","$GENE_SOURCE";
ptag "type","$GENE_TYPE";
btag "basic_information"; 
ptag "GeneID", $ogid;
ptag "species", $CLADE;
ptag "ntaxa", $nt;
ptag "ngene", $ng;
ptag "occurrence", $occ; # turn into xml subtags?? need numeric index of each
ptag "group-evalue", $sev;
ptag "group-identity", $spi;
ptag "description", $ogdef if($ogdef); #? group-def or def, same tag as similarity.def ?
ptag "date",  $date;
etag "basic_information"; 

btag "similar_genes";  
  foreach my $gn1 (@$gna) {
    my($sppgn,$spp,$gn)= fix_sppid($gn1); # my($spp,$gn)= split /[:_]/, $gn1, 2;
    my $sna= $sppmap{$spp} || $spp;
    my($og,$ng,$ev,$pi,$g2,$e2,$p2) = (defined $gsum{$gn1}) ? @{ $gsum{$gn1} } :
      (defined $gsum{$sppgn}) ? @{ $gsum{$sppgn} }: ('.') x 10; # err if none?
    my $gdef= $gdef{$gn1} || $gdef{$sppgn} || "";
    # def: cleaned out dingbat chars below: <>&;
    
    my ($gdbx)= ""; if($gdef =~ s/Dbxref=(\S+)\s*//i) { $gdbx=$1; }
    $gdef =~ s/desc=//;
    # Dbxref=RefSeq:XM_394073,gi:110763653 desc=PREDICTED: Apis mellifera similar to CG3295-PA 
    
    # (my $pgn=$gn); #?? =~ s/_/:/; # db prefix style? or leave in for indexing 
    $ev= sprintf "%.1g",$ev; # drop extra digits ?
    btag "similarity", 1;
    ptag "spp",$sna,1;
    ptag "acc",($SPPGN?$sppgn:$gn),1; # sppgn was gn1
    ptag "eval",$ev,1;
    ptag "iden",$pi,1;
    # if(my $aalen= $aalen{$sppgn}) { ptag "len",$aalen,1; } # put for all or none, 0=missing? if %aalen exists
    if($nAAlen>0) { my $aalen= $aalen{$sppgn}||0; ptag "len",$aalen,1; }
    ptag "def",$gdef,1 if($gdef);# or null value?
    ptag "dbx",$gdbx,1 if($gdbx);# or null value?
    etag "similarity", 1;
  }
etag "similar_genes"; 

etag "GeneSummary";  print "\n";
}

# 2014 add: genegroupbpo.pl fish12a_omclgn2sum.tab .. > look for fish12a.aa.qual/size
sub lookForAAlen {
  my($inf)= @_;
  my $fnam=$inf; $fnam=~s/_[\w\.]+$//;
  my $fn= "$fnam.aa.qual"; # readdir ?
  $fn="$fnam.aa.count" unless(-f $fn); 
  $fn="$fnam.fa.count" unless(-f $fn); 
  unless(-f $fn) { warn "#WARN: missing $fnam.aa.qual to report aa lengths\n"; return 0; }
  my $n=0; open(F,$fn); 
  while(<F>) { my($id,$aw)=split; if($aw){ $id = fix_sppid($id); $aalen{$id}=$aw; $n++; } } close(F);
  return $n;
}



#..................... input ...............
sub MAIN {

my $needspp=1; %sppmap=%sppshort=();
if($speciesmap=~/^FISH11/i) {  %sppmap=%FISH11sppmap; %sppshort=%FISH11sppshort; }
elsif($speciesmap=~/^FISH12/i) {  %sppmap=%FISH12sppmap; %sppshort=%FISH12sppshort; }
elsif($speciesmap=~/^PLA9/i) {  %sppmap=%PLA9sppmap; %sppshort=%PLA9sppshort; }
elsif($speciesmap) {
  map{ my($k,$v)=split/[=:]/,$_,2; 
    $sppmap{$k}=$v; 
    $sppshort{$k}= ucfirst($k); # OR split "_",$v;
  } split/[,;|]+/, $speciesmap;
}
$needspp=0 if(scalar(%sppmap));

print "<?xml version=\"1.0\"?>\n" if $doxml;
btag "GeneSummaries" if $doxml; 

print_ugpdummy() if $doxml; # this is for lucegene reports
#  # because 1st record has extra <xml> tag that web browsers hate  
#  # lucegene returns this  1st record  (a:= added by reporter; s:this sourcefile)
#a: <?xml version="1.0" encoding="ISO-8859-1"?>
#a: <?xml-stylesheet href="/templates/xsl/arthropodxml.xsl" type="text/xsl"?>
#a: <GeneSummaries>
#a:<!-- docurl="arthropod_ogenesum3.ugp.xml,0-56892" -->
#s: <?xml version="1.0"?>  # extra
#s: <GeneSummaries>   # extra
 
# add debug stuf?
my ($linfile,$nin,%instat);
$nAAlen= lookForAAlen($ARGV[0]); # 2014 add: genegroupbpo.pl fish12a_omclgn2sum.tab .. > look for fish12a.aa.qual/size
  
while(<>){

  my $infile=$ARGV; 
  $instat{$infile}{open}++ if($infile ne $linfile);
  $nin++;
  
## dang, this list still includes alt-tr removed elsewhere ; require gsum;
if(/^ORTHOMCL(\d+).(\d+) genes,\s*(\d+) taxa.:\s+(.+)$/) { # omclfilt1_Aug_10/all_orthomcl.out.gz
  my($og,$ng,$nt,$gn)= ($1,$2,$3,$4);
  my @gn= split" ",$gn; # cut spp: prefix?
  my @gn2=();
  foreach my $gn1 (@gn) {
    my($sppgn,$spp,$gn)= fix_sppid($gn1); # my($spp,$gn)= split /[:_]/, $gn1, 2;
    push(@gn2,$sppgn) if (defined $gsum{$gn1} or defined $gsum{$sppgn});
###    Similarity: Drosophila melanogaster, dmel_NP_476765.1, ., .;
### ^^ these are removed alt-tr *** drop
  }
  
  $instat{$infile}{omclgroup}++;
  print_genesum($og,$ng,$nt,\@gn2);

} elsif(/^>/ or /^$GTAG/) { # aabugs.deflines
  chomp;  s/>//;   
  s/^gi\|\d+\|\w+\|//; # do we have any ncbi gi junk to cut? >gi.nnn.ref.ID_HERE
  my($g,$d)= split /[\s\|]+/,$_,2; $g=~s/\W+$//;
  # add consensus deflines:  >ARP1_G299.consensus Rho1 CG8416-PD; src=dmel_NP_725524.1
  $g =~ s/\.consensus$//; # dont need, ARP1 will be uniq pattern
  $g = fix_sppid($g);
  # def now may be name table w/ \t\t\t : drop
  $d =~ s/\t.*$//;
  # def: clean out dingbat chars here 
  $d =~ s/[<>\&]/ /g; $d =~ s/[;]/,/g;  
  $gdef{$g}= $d;
  $instat{$infile}{defline}++;
  
} elsif(/\t$GTAG/) {  # aabugs_omclgn1sum.tab aabugs_omclgn2sum.tab

# ------> aabugs_omclgns2.tab.gz <------
# aedes_AAEL015645-PA	culex_CPIJ009686	ARP21863	7e-99	78
# aedes_AAEL001757-PA	culex_CPIJ004938	ARP7181	6e-152	67
# aedes_AAEL001757-PA	anopheles_AGAP005778-PA	ARP7181	9e-129	61

# ------> aabugs_omclgn2sum.tab.gz <------
# aedes_AAEL015645-PA     ARP21863        1       7e-99   78      culex_CPIJ009686        7e-99   78
# aedes_AAEL001757-PA     ARP7181 9       2.222e-74       51      culex_CPIJ004938        6e-152  67
# aedes_AAEL001761-PA     ARP6403 10      9.01e-38        49      anopheles_AGAP012826-PA 6e-63   57

## DAMIT: omcl tables use/require 'species_id', others have species:id
## DAMIT: species_id NOT species:id crap again; (spp,id)= split /[:_]/,$id,2; << is this ok?
## DAMIT: names have spp:id; here need species_id, what else ??

  chomp;
  my($g1,$og,$ng,$ev,$pi,$g2,$e2,$p2)= split"\t";
  
  ## DAMIT: make spp:id same as omcl spp_id; other dingbats here?
  ##$g1= fix_sppid($g1);
  ##$g2= fix_sppid($g2);

  my($spp1,$spp2);
  ($g1,$spp1)= fix_sppid($g1);
  ($g2,$spp2)= fix_sppid($g2);
  if($needspp) { $sppmap{$spp1}=$spp1;  $sppmap{$spp2}=$spp2; } # must collect all before print_genesum()
  
  ## $ev= sprintf "%.1g",$ev; # drop extra digits ?
  $gsum{$g1}= [$og,$ng,$ev,$pi,$g2,$e2,$p2]; 
  $instat{$infile}{genepair}++;
 
 }
 $linfile= $infile;
}

etag "GeneSummaries" if $doxml; 

if($debug) {
  for my $inf (sort keys %instat) {
    my @t= sort keys %{$instat{$inf}};
    print STDERR "#IN\t",$inf;
    for my $t (@t) { my $c=$instat{$inf}{$t}||0; print STDERR "\t$t:$c"; } 
    print STDERR "\n";
    }
}

}

__END__


# dgg: sample .bpo
# 0;  1                ; 2 ;    3              ; 4 ;  5   ; 6 ;   7
# 2;aedes_AAEL015645-PA;264;aedes_AAEL015645-PA;264;1e-145;100;1:1-250:1-250.


# revised steps, 2009.12.2

# input tables for ugp summary output

set omclbpo=aabugs_filtered.bpo.gz
    # this is gene-pair reciprocal blastp summary input to orthomcl, as
    #  num;gid1;len1;gid2;len2;evalue;pct-ident;hspalign1.hspalign2...
    # 1;aedes_AAEL015645-PA;264;aedes_AAEL015645-PA;264;1e-145;100;1:1-250:1-250.
    # 499;aedes_AAEL009771-PA;569;culex_CPIJ003773;619;5e-49;29;1:18-528:34-571.
    # 54995;acyr1_ncbi_hmm538284;371;tcas3_GLEAN_08954;541;3e-24;35;1:109-309:69-272.2:150-307:319-473.

set genegrp=aabugs_omclgn.tab 
    # this is (ARPid geneid) list, from all_orthomcl.out ??

set protdef=pro2/*.aa.deflines.gz
    # annotations per-gene with dbxref= desc= 
    
set omclout=omclfilt1_Aug_10/all_orthomcl.out.gz
    # this is orthomcl groupings output, as
    # ORTHOMCL0(414 genes,2 taxa):     daphnia:daphnia_NCBI_GNO_226344 ...
    # ORTHOMCL499(25 genes,3 taxa):	 aedes:aedes_AAEL008967-PA ..

per-species id filter lists:
  dmel_protid_keep.tab : dmel NCBI gene ids to keep (others are poorer alt-tr)
    n=12717 (missing some?; should redo w/ newer dmel gene set)
  dmel_NP_524618.1	dmel_CG2168-PA	Ribosomal protein S3A 
  dmel_NP_651908.1	dmel_CG1587-PA	Crk 
  
  alttr_*.ids : other species alt-tr ids to dorp
     179 alttr_amel4_gnohmm.ids
     137 alttr_daphnia.ids
     342 alttr_dmoj.ids
     698 alttr_dpse.ids
     530 alttr_other.ids
    1886 total

  ARP1 groups post-identified as likely transposons, keep / mark or drop?
  617 sums/arp_teall.ids
  136 sums/arp_teaphid.ids

  ## paste together bpo and arp ids
  ##  revise this to filter alttr_*.ids as below
  gunzip -c $omclbpo | cat  aabugs_omclgn.tab  dmel_protid_keep.tab - \
   | perl pastebpo.pl > aabugs_omclgns2.tab 
 
  # ------> aabugs_omclgns2.tab.gz <------
  # aedes_AAEL015645-PA	culex_CPIJ009686	ARP21863	7e-99	78
  # aedes_AAEL001757-PA	culex_CPIJ004938	ARP7181	6e-152	67
  # aedes_AAEL001757-PA	anopheles_AGAP005778-PA	ARP7181	9e-129	61
  
  # ------> aabugs_omclgn2sum.tab.gz <------
  #   gene1                 arpid     ngene   ave-eval   ave-pi     gene2 .. gene2eval  gene2pi < ignore these 3
  # aedes_AAEL015645-PA     ARP21863        1       7e-99   78      culex_CPIJ009686        7e-99   78
  # aedes_AAEL001757-PA     ARP7181 9       2.222e-74       51      culex_CPIJ004938        6e-152  67
  # aedes_AAEL001761-PA     ARP6403 10      9.01e-38        49      anopheles_AGAP012826-PA 6e-63   57

  # ** input is aabugs_omclgns2.tab or  aabugs_omclgn2sum.tab ?? which
  # >> use aabugs_omclgn2sum, that is one row/gene; ignoring 2nd gene
  # >> only genes in aabugs_omclgn2sum.gene1 column are reported, 
  #    groupings come from $omclout (and aabugs_omclgns2 arpid?)
  
  ## produce gene group docs
  # inputs are: 1. per-gene id, arpid, eval, pi scores; 2. per-gene annots, 3. genes per omcl group
  #  need to add 4. outmcl080826/apr2anno.clusters.tab : mcl clusters of groups
  # gene ids must be in both 1 and 3 for output
  
  gunzip -c aabugs_omclgns2.tab.gz $protdef $omclout |\
    env xml=1 idprefix=ARP1_G date=20080824 perl genegroupbpo.pl \
     >  arthropod_ogenesum4.ugp.xml

#....... genegrp.idtab from omclout .........

gunzip -c $omclout | env idpre=ARP1_G perl -ne \
'BEGIN{ $idpre=$ENV{idpre}||"Omcl"; print join("\t",qw(gene group ngene ntaxa)),"\n"; }\
if(/^ORTHOMCL(\d+).(\d+) genes,(\d+) taxa.:\s+(.+)$/) {  
my($og,$ng,$nt,$gn)= ($1,$2,$3,$4);
map{ ($sp,$g)=split":"; print "$g\t$idpre$oid\t$sp\t$ng\t$nt\n"; } split" ",$gn; }' \
> gene_mclgroup.idtab

#... filter idtab
# .. Change: this should only filter alt-tr in same group, not in separate groups

# grep -v -F -f alttr.ids gene_mclgroup.idtab > gene_mclgroup-noalt.idtab

#.... merge .bpo, gene.deflines, group.idtab to gene_group.xml (update genegroupbtp.pl)


#------- insert to ugp.xml -------

cat apr2anno.newick | perl bugmcl8.info > apr2anno.clusters.tab

gunzip -c data/arthropod_ogenesum4.ugp.xml.gz | cat outmcl080826/apr2anno.clusters.tab - |\
perl -ne'if(/^ARC1_/){ chomp; ($c,$gc)=split"\t"; while( $gc =~ m/(ARP\w+)/g) { $gc{$1}= $gc; } } \
else { if(m/<GeneSummary/) { $gc= (m/euGenes:(\w+)/) ? $gc{$1} : ""; } \
elsif( $gc and m,</GeneSummary>,){ print "<related_gene_groups>\n$gc\n</related_gene_groups>\n"; }  \
print; }'\

#..................................
# update 200901 to correct/remove alltr from group counts

# alttr_dpse.ids alttr_other.ids ..

set spl="acyr,aede,amel,anop,cule,daph,dmel,dmoj,dpse,ixod,naso,pedi,tcas"

gunzip -c omclfilt1_Aug_10/all_orthomcl.out.gz | \
cat dmel_pa.ids alttr_*.ids - | \
env spl=$spl perl -ne \
'if(/^(dmel\S+)/) { $dmpa{"dmel:".$1}++; } \
elsif(/^(acyr1|amel4|dmoj|dpse|nasonia|tcas3)(\S+)/) { $alttr{"$1:$1$2"}++; } \
elsif(/ORTHOMCL/) { ($om,$gn)=split /:\s+/,$_,2;  @gn= split" ",$gn; \
@sp{@spl}= (0)x20; %didg=();  map{ $sg=$_; $spn=substr($_,0,4); \
$ok=1; if($spn=~/dmel/){ $ok= $dmpa{$_}||0; } \
elsif($spn=~/acyr|amel|dmoj|dpse|naso|tcas/){ $ok=($alttr{$_}) ? 0:1;} \
elsif($spn=~/aede|anop|pedi/){ my($g1,$tn)= split/\-/; $ok=($didg{$g1}++)? 0:1;} \
$sp{$spn}++ if ($ok); $drop{$spn}++ unless($ok); } @gn; \
$om =~ m/ORTHOMCL(\d+).(\d+) genes,(\d+) taxa/; ($og,$g,$t)=($1,$2,$3); \
@sc= @sp{@spl}; print join("\t",$og,$t,$g,@sc),"\n"; }\
END{ print join("\t", "#drop",0,0,@drop{@spl}),"\n"; } \
BEGIN{$spl=$ENV{spl}; @spl=split",",$spl; @drop{@spl}=(0)x20; print join("\t","OID","Nt","Ng",@spl),"\n";} '  \
> aabugs-orthomcl-flt1-count-pa3.tab



# run2: use arp ortho/paralog set for 10+ species;
#  selection of arp genes: 
#    1. remove alttr gene ids: $aa/ -alttr*.ids and +dmel_pa.ids
#    2. remove TE-groups (arp_teall.ids )
#    3. arp groups with 1+ species having 2+ paralogs

# reduce species set?  drop ixod,aede,dpse?
#set spl="acyr,aede,amel,anop,cule,daph,dmel,dmoj,dpse,ixod,naso,pedi,tcas"
set spl="acyr,amel,anop,cule,daph,dmel,dmoj,naso,pedi,tcas"

gunzip -c omclfilt1_Aug_10/all_orthomcl.out.gz | \
cat sums/arp_teall.ids sums/arp_teaphid.ids aabugs_filtered.aalen \
  dmel_pa.ids alttr_dmoj.ids alttr_dpse.ids alttr_other.ids - | \
env spl=$spl perl -ne \
'if(/^(\S+)\t(\d+)$/) { ($g,$na)=($1,$2); ($s)=$g=~/^([^_]+)/; $sz{"$s:$g"}=$na; }\
elsif(/^(dmel\S+)/) { $dmpa{"dmel:".$1}++; } \
elsif(/^(acyr1|amel4|daphnia|dmoj|dpse|nasonia|tcas3)(\S+)/) { $alttr{"$1:$1$2"}++; } \
elsif(/^ARP1_G(\d+)/){ $gdrop{$1}++; }\
elsif(/ORTHOMCL/) { ($om,$gn)=split /:\s+/,$_,2;  @gn= split" ",$gn; \
@gns=(); @sp{@spl}= (0)x20; %didg=();  map{ $sg=$_; $spn=substr($_,0,4); \
$ok=1; if($spn=~/dmel/){ $ok= $dmpa{$_}||0; } \
elsif($spn=~/acyr|amel|daph|dmoj|dpse|naso|tcas/){ $ok=($alttr{$_}) ? 0:1;} \
elsif($spn=~/aede|anop|pedi/){ my($g1,$tn)= split/\-/; $ok=($didg{$g1}++)? 0:1;} \
$sp{$spn}++ if ($ok);  push(@gns,$sg) if($ok and $skeep{$spn}); } @gn; \
$om =~ m/ORTHOMCL(\d+).(\d+) genes,(\d+) taxa/; ($og,$ng,$nt)=($1,$2,$3); \
next if($gdrop{$og}); \
@sc= @sp{@spl}; $smax=0; map{$smax=$_ if $_>$smax; }@sc; next unless ($smax>1); \
@sz= sort{$b <=> $a} map{ $sz{$_} } @gns; \
$msz= int( 0.7 * $sz[ 0] ); \
@gns= grep { $sz{$_} >= $msz } @gns; \
print join("\t","ARP1_G".$og,@gns),"\n"; }\
BEGIN{$spl=$ENV{spl}; @spl=split",",$spl; @skeep{@spl}=(1)x20; @drop{@spl}=(0)x20; }'  \
> aabugs-omcl2noaltte-pagenes.list


# index ugp.xml

/bio/argos/daphnia/webapps/lucegene
admin/lucegene-index.sh -lib arthropodxml -run > & log.arugp &

# build score table?
# add g2 ntaxa to g1 gene scores? add g1 prot-len

gunzip -c aabugs_filtered.bpo.gz | cat dmel_prota.ids - | perl -ne\
'if(/^(dmel_\S+)/) { $dmprime{$1}++; } else { \
($g1,$g2,$ev,$pi)= (split";")[1,3,5,6]; \
next if($g1 eq $g2 or ($g1 =~ /^dmel/ && not $dmprime{$g1}) or ($g2 =~ /^dmel/ && not $dmprime{$g2})); \
dumpg() unless($g1 eq $lg); $lg=$g1;  $ng++; $sev+= $ev; $spi += $pi; } \
sub dumpg{ if($ng>0) { $aev=sprintf"%.4g",$sev/$ng; $api=int($spi/$ng); \
print join("\t",$lg,$ng,$aev,$api),"\n"; } \
$lg=$ng=$sev=$spi=0; } \
END{dumpg(); }' \
| more
...
aedes_AAEL001877-PA     99      3.931e-09       27
aedes_AAEL009959-PA     12      0       94
aedes_AAEL009965-PA     19      3.158e-07       37   # not all of these belong in ortho group
aedes_AAEL009970-PA     76      1.068e-16       37
aedes_AAEL009974-PA     225     1.344e-28       43
...

microbe% head aabugs_omclgn.tab
ARP0	daphnia_NCBI_GNO_226344
ARP0	nasonia_NCBI_hmm10024
ARP0	nasonia_NCBI_hmm100784
 : swap cols to geneid,arpid
..

# gunzip -c aabugs_filtered.bpo.gz | \
#  cat omclfilt1_Aug_10/aabugs_omclgn.tab dmel_prota.ids - \
#  | perl pastebpo.pl > aabugs_omclgns.tab &

grep dmel all_orthomcl_dm1.out | perl -ne\
'($og)=m/ORTHOMCL(\d+)/; @g= m/(dmel_\S+)/g; print "ARP$og\t"; @d=(); \
foreach $g (@g) { ($c=$g)=~s/\-P\w+//; print "$g " unless($c eq $lc); \
push(@d,$g) if($c eq $lc);  $lc=$c; } print "\n"; print "# dropped: @d \n" if(@d);' \
> dmel_keepids.tab


grep -v '# dropped' dmel_keepids.tab | perl -ne'@v=split; shift(@v); print join("\n",@v),"\n";' |\
sort | uniq | cat - dmel_protid.tab | perl -ne\
'if(/^(dmel_CG\S+)/){$cg{$1}++;} elsif(/^dmel/){($n,$c,$d)=split"\t"; print if($cg{$c});} '\
> dmel_protid_keep.tab


gunzip -c aabugs_filtered.bpo.gz | \
 cat  aabugs_omclgn.tab  dmel_protid_keep.tab - \
 | perl pastebpo.pl > aabugs_omclgns2.tab &



# strip alttr -P[B-Z] here only? for aedes/anopheles/pediculus ?

cat aabugs_omclgns2.tab | perl -ne\
'($g1,$g2,$og,$ev,$pi)=split;  \
next if(not($og =~ /ARP/) or $g1 =~ /\-P[B-Z]/ or $g2 =~ /\-P[B-Z]/); \
dumpg() unless($g1 eq $lg); \
$lg=$g1; $log=$og; $ng++; $sev+= $ev; $spi += $pi; \
if($pi>$mpi){ $mpi=$pi; $mg=$g2; $mev=$ev; } \
sub dumpg{ if($ng>0) { $aev=sprintf"%.4g",$sev/$ng; $api=int($spi/$ng); \
print join("\t",$lg,$log,$ng,$aev,$api,$mg,$mev,$mpi),"\n"; } \
$mpi=$mg=$mev=$lg=$ng=$sev=$spi=0; } \
END{ dumpg(); }' \
> aabugs_omclgn2sum.tab

# ^^^^^^^^^^^^^ 
 this is input to genegroupbpo  with list of transcript/gene ids to keep (filter alt-tr)
 change dmel_prota.ids; need other than -PA; 

#............. text format

GeneSummary id="euGenes:ARP1_G458":
  title: Arthropod gene cluster ARP1_G458;
  source: euGenes;
  type: gene_group;
  basic_information:
    GeneID: ARP1_G458;
    species: Arthropoda;
    ntaxa: 8;
    ngene: 26;
    group-evalue: 1.4e-119;
    group-identity: 60;
    occurrence: Aedes:1, Anopheles:1, Aphid:0, Apis:0, Culex:1, Daphnia:0, DrosMel:6, DrosMoj:5, DrosPse:8, Ixodes:0, Nasonia:0, Pediculus:1, Tribolium:1;
    date: 20080814;
  similar_genes:
    similarity: spp: Aedes aegypti; acc: aedes_AAEL006975-PA; eval: 3e-137; iden: 58; 
    similarity: spp: Anopheles gambiae; acc: anopheles_AGAP007904-PA; eval: 3e-142; iden: 56; 
    similarity: spp: Culex pipiens; acc: culex_CPIJ003539; eval: 2e-140; iden: 57; 
    similarity: spp: Drosophila melanogaster; acc: dmel_NP_610892.1; eval: 9e-129; iden: 63; 
    similarity: spp: Drosophila melanogaster; acc: dmel_NP_729613.1; eval: 3e-141; iden: 61; 
    similarity: spp: Drosophila mojavensis; acc: dmoj_NCBI_GNO_32440339; eval: 9e-137; iden: 58; 
    similarity: spp: Drosophila mojavensis; acc: dmoj_NCBI_GNO_32695044; eval: 4e-141; iden: 64; 
    similarity: spp: Drosophila pseudoobscura; acc: dpse_NCBI_GNO_32037620; eval: 3e-118; iden: 63; 
    similarity: spp: Pediculus humanus; acc: pediculus_PHUM005199-PA; eval: 4e-128; iden: 48; 
    similarity: spp: Tribolium castaneum; acc: tcas3_GLEAN_01175; eval: 4e-136; iden: 52; 

GeneSummary id="euGenes:ARP1_G564":
  title: Arthropod gene cluster ARP1_G564;


....

pull similarity Description field from .aa defline where possible

proteomes/aedes.aa.gz :
>aedes_AAEL000012-PA Gustatory receptor 61a, putative

proteomes/culex.aa.gz :
>culex_CPIJ000007|chitinase|protein_coding|supercont3.1|419596|432057

proteomes/dmel.aa.gz :
>dmel_NP_524618.1| Ribosomal protein S3A CG2168-PA, isoform A [Drosophila melanogaster]


foreach aa ($aadef)
 set aan=`echo $aa | sed -e's/.gz//'`
 gunzip -c $aa | grep '^>'  > $aan.deflines
end 

 # do this in here: perl -ne's/^>//; ($g,$d)=split/[\s\|]/,$_,2; $g=~s/\W+$//; print "$g\t$d";'
 
 
#...... more deflines annot and protein.aa version 2 w/ updated ARP id

gunzip -c nasonia-refseq.deflines.gz nasonia_gnomon_refseq.txt.gz | perl -ne'chomp; if(s/^>(gi\S+) //){$gi=$1; $gi=~s/\W+$//; s/, mRNA//; $gi=~s/gi./gi:/; $gi=~s/.ref.(\w+)/,RefSeq:$1/; $r=$1; $rd{$r}="$gi desc=$_";} elsif(/^lcl.(\w+)\s+(\S+)/){$h=$1; $r=$2; $rd=$rd{$r}; print ">nasonia_NCBI_${h} Dbxref=$rd\n";}' > nasonia.deflines

gunzip -c aphid-refseq.deflines.gz aphid_gnomon_refseq.txt.gz | perl -ne'chomp; if(s/^>(gi\S+) //){$gi=$1; $gi=~s/\W+$//; s/, mRNA//; $gi=~s/gi./gi:/; $gi=~s/.ref.(\w+)/,RefSeq:$1/; $r=$1; $rd{$r}="$gi desc=$_";} elsif(/^7029/){@v=split"\t"; ($h,$r,$p)=@v[2,6,7]; $rd=$rd{$r};  print ">acyr1_ncbi_${h}"; print " Dbxref=RefProt:$p,$rd" if($p); print "\n";}' > aphid.aa.deflines

## Note Apis rna.asn Gnomon IDs all are 912345 whereas prot hmm ids drop 9 and leading 0 if < 10000

gunzip -c /bio/biomirror/ncbigenomes/Apis_mellifera/rna/rna.asn.gz | perl -ne'if(/str "ModelId"/){$m=1;}elsif($m==1 and /str "(\d+)"/){$hm=$1; $m=2;} elsif($m==2 and /accession "(\w+)/){ $r=$1; $m=3; } elsif($m==3 and /gi (\d+)/){$gi=$1; $m=4; } elsif($m==4 and /title "(.+)/) { $ti=$1; $m=5; $m=6 if(s/" ,//); } elsif($m==5){ s/^\s+//; $ti .= $_; $ti=~s/" ,//; $m=6;} if($m==6){ chomp($ti); $ti=~s/, mRNA//; $hm=~s/^9[0]*//; print ">amel4_ncbi_hmm$hm Dbxref=RefSeq:$r,gi:$gi desc=$ti\n"; $m=$hm=$r=$gi=$ti=""; }' > apis.aa.deflines

#...........  vers2.aa

gunzip -c aphid.aa.gz | cat aphid.aa.deflines aabugs_omclgn.tab - | perl -ne'if(/^(ARP\w+)\t(\S+)/){ $arp{$2}=$1; $fa=1; } elsif($fa and /^>(\S+)/){$g=$1; $fa=2; $om=$arp{$g};  $df=$df{$g}; $_=$df if($df); if($om){ $om=~s/ARP/ARP1_G/; unless(s/dbxref=/dbxref=euGenes:$om,/){ s/$/ dbxref=euGenes:$om/;} } } elsif(/^>(\S+)/){ $df{$1}=$_;} print if($fa>1);' > aphid2.aa

gunzip -c apis.aa.gz | cat apis.aa.deflines aabugs_omclgn.tab - | perl -ne'if(/^(ARP\w+)\t(\S+)/){ $arp{$2}=$1; $fa=1; } elsif($fa and /^>(\S+)/){$g=$1; $fa=2; $om=$arp{$g};  $df=$df{$g}; $_=$df if($df); if($om){ $om=~s/ARP/ARP1_G/; unless(s/dbxref=/dbxref=euGenes:$om,/){ s/$/ dbxref=euGenes:$om/;} } } elsif(/^>(\S+)/){ $df{$1}=$_;} print if($fa>1);' > apis2.aa

gunzip -c daphnia.aa.gz | cat daphnia.aa.deflines aabugs_omclgn.tab - | perl -ne'if(/^(ARP\w+)\t(\S+)/){ $arp{$2}=$1; $fa=1; } elsif($fa and /^>(\S+)/){$g=$1; $fa=2; $om=$arp{$g};  $df=$df{$g}; $_=$df if($df); if($om){ $om=~s/ARP/ARP1_G/; unless(s/dbxref=/dbxref=euGenes:$om,/){ s/$/ dbxref=euGenes:$om/;} } } elsif(/^>(\S+)/){ $df{$1}=$_;} print if($fa>1);' > daphnia2.aa

gunzip -c tribolium.aa.gz | cat tribolium.aa.deflines aabugs_omclgn.tab - | perl -ne'if(/^(ARP\w+)\t(\S+)/){ $arp{$2}=$1; $fa=1; } elsif($fa and /^>(\S+)/){$g=$1; $fa=2; $om=$arp{$g};  $df=$df{$g}; $_=$df if($df); if($om){ $om=~s/ARP/ARP1_G/; unless(s/dbxref=/dbxref=euGenes:$om,/){ s/$/ dbxref=euGenes:$om/;} } } elsif(/^>(\S+)/){ $df{$1}=$_;} print if($fa>1);' > tribolium2.aa

gunzip -c nasonia.aa.gz | cat nasonia.aa.deflines aabugs_omclgn.tab - | perl -ne'if(/^(ARP\w+)\t(\S+)/){ $arp{$2}=$1; $fa=1; } elsif($fa and /^>(\S+)/){$g=$1; $fa=2; $om=$arp{$g};  $df=$df{$g}; $_=$df if($df); if($om){ $om=~s/ARP/ARP1_G/; unless(s/dbxref=/dbxref=euGenes:$om,/){ s/$/ dbxref=euGenes:$om/;} } } elsif(/^>(\S+)/){ $df{$1}=$_;} print if($fa>1);' > nasonia2.aa
