# evgpubsetsum.pm

# package evgpubsetsum;
package main;

use FindBin;
use lib ("$FindBin::Bin","$FindBin::Bin/.."); # assume evigene/scripts/ layout: main, prot/, rnaseq/, ...
use strict;
use warnings; 
use evigene_pubsets; # now has some of below subs

#shared from evigene_pubsets.pm
our $IDPREFIX; # caller sets = $ENV{idprefix} || 'EVGm'; 
our $ORGANISM; # caller sets
our @ANNO_COLS;
our (%pubids,%pubidinfo,%puboidc,%puballoids); ## cdna_evigenesub globals 
our (@publicset,@submitset); # tidyup file sets

our($TABLE_G1,@TABLE_G1_KEYS); # local

sub geneset_summary {
  my($pubids, $annotab)=@_;
  my($npubid, $pubidh, $annoth);

  if(ref($pubids)) {
    $pubidh= $pubids;
  } else {
	  ($npubid, $pubidh)= read_pubids($pubids); 
	}
  # global: %pubidinfo has details; %$pubidh has only pubid,oid=>pubid
  # my($oid,$gid,$alti,$class,$aqual,$pia,$notes,$alloids)= split"\t", $pubidinfo{$pubid};
  # notes = aaref:nnn,,, chrmap:nnn,,, other:...
  
  if(ref($annotab)) {
    $annoth= $annotab;
  } else {
    my %annot= read_annotab($annotab); #** fill in annothash from pubidh if no table
    $annoth= \%annot;
  }
  
  # $annothash{$oid}= $annothash{$pubid}= ann.txt row line (tabbed) of cols:
  unless(@ANNO_COLS > 8) { # dang perl package mixup
    @ANNO_COLS= qw(PublicID OrigID TrLen CDSoff AAqual TrGaps Dbxref Namealign Product_Name CDD_Name Class Maploc Mapqual);
  }
  
  our(%sum);
  #x.my @pd=sort keys %$annoth; # wrong, annoth has both pubid,oid keys to same val
  my @pd=sort keys %pubidinfo; # this has only pubid keys

  # debug mapquals in _pubsets.pm
  #d sub mapquall{ my($mq)=@_; my @m= $mq=~m/(\d+)a,(\d+)i,(\d+)l,(\d+)x/; return (@m>3)?@m:(0,0,0,0); }
  
  for my $pd (@pd) {
    my @van=split"\t",$annoth->{$pd};
    my %van=(); for my $i (0..$#van){ $van{ $ANNO_COLS[$i] }= $van[$i]; } #??
    
    my($oid,$gid,$alti,$evgclass,$aqual,$pia,$notes,$alloids)= (0) x 9;
    my($aaref,$dbxref,$nap,$name,$cddname,$maploc,$mapqual)= (0) x 9;
    my($hasho,$hasin,$goodmap, $cov, $ident, $clen, $nexon, $aafull)= (0) x 9;
    
    ($oid,$gid,$alti,$evgclass,$aqual,$pia,$notes,$alloids)= split"\t", $pubidinfo{$pd};


    #debug: ($dbxref,$nap,$name,$evgclass, $maploc,$mapqual)= @van[6,7,8,10,11,12];
    if($van{'Class'}) { 
    ($dbxref,$nap,$name,$evgclass, $maploc,$mapqual)= map{ $van{$_} || 0 } 
       qw(Dbxref Namealign Product_Name Class Maploc Mapqual);  # is this bad?
    }
    
    ( $cov, $ident, $clen, $nexon)= mapquals($mapqual || $notes);
    
    if( my($arf)= $notes=~m/aaref:([^;\s]+)/) { # bug: 'aaref:0,0,' is bad, new bug from where?
      my($av)= $notes=~m/aaref:(\d+)/?$1:0;
      $arf=~s/,(chrmap|pflag|tscore).*//;
      $aaref=$arf if($av>0);
    }
    
    #UPD1905: fixme dbxref + Namealign valid but no ":" in dbxref ..
    #o $hasho= ($aaref or $dbxref =~ /:/)?1:0; # bad: nap = 'na' or 0 (w/ homol) or 100%,596/596,596
    $hasho= ($aaref or $dbxref =~ /:\w/ or ($dbxref=~/\w/ and $nap =~/%,/))?1:0; # bad: nap = 'na' or 0 (w/ homol) or 100%,596/596,596
    # hasho=1 if($nap=~/%,/); # dont bother. other two capture ho .. not if no ':' !!
    $hasin= ($nexon > 1)?1:0;
    
    my $pai= $cov * $ident/100; # cov only? or cov * ident
    # $goodmap= ($pai > 79)? "good" : ($pai < 10)? "nomap" : "poor";  # 3 levels: good > poor > nomap
    $goodmap= ($pai > 79)? 2 : ($pai < 10)? 0 : 1;  # 3 levels: good > poor > nomap
    
    # $aafull = ($aqual=~/(complete|partial)/)?$1:"other";
    $aafull = ($aqual=~/complete/)?1:0;

    #a my $ecla= ($evgclass =~ /^(drop|cull)/)?$1:"pub";  # cull1,cull2 ..
    #x my $eclb= ($evgclass =~ /^$ecla(main|noclass|alt|part|noncode)/)?$1:"other";
    #x my $eclb= ($evgclass =~ /^$ecla\d*(main|noclass|alt|part|noncode)/)?$1:"other";
    #x my $eclb= ($evgclass =~ /(main|noclass|alt|part|noncode)/)?$1:"other";
    
    my($ecla,$eclb)= $evgclass=~/^(\w*)(main|noclass|alt|part|noncode)/ ? ($1,$2) : ("pub",$evgclass);
    $ecla ||="pub";
    $eclb=~ s/part/alt/; 
    $eclb=~s/(main|noclass)/coding/;
    
    $sum{$ecla}{type}{$eclb}++;
    
    # need some intersections:  goodmap & hasin ..
    # per locus: skip alt, or bestof(main|alt) ?
    if($eclb =~ /coding/) {
      $sum{$ecla}{hasho}{$hasho}++; # 1/0 or "hasho/noho" ?
      $sum{$ecla}{aafull}{$aafull}++; # complete/partial or 1/0 ?
      $sum{$ecla}{goodmap}{$goodmap}++; # 2/1/0 or good/poor/nomap ?
      $sum{$ecla}{hasin}{$hasin}++ if($goodmap>1); # 1/0 or "hasint/1exon" ?
      $sum{$ecla}{hasho}{noint}++ if($goodmap>1 and $hasho>0 and $hasin==0); # 1/0 or "hasint/1exon" ?
    
    } elsif($eclb=~/alt/) {
 
      $sum{$ecla}{altof}{$gid}++;
   
    }
    
  }
  
  sub pct{ my($vyes,$vt)=@_; my $p=($vt>0)?int(0.5 + 100*$vyes/$vt):0; $p=100 if($p>100); $p; }
  sub pctof{ my($c,$t,$kyes,$kno,$kc)=@_; our(%sum);
    $kc= -1 unless(defined $kc); 
    my($vyes,$vno,$vc)= map{ $sum{$c}{$t}{$_}||0 } ($kyes,$kno,$kc); 
    my $vt=$vyes+$vno+$vc; my $p=($vt>0)?int(0.5 + 100*$vyes/$vt):0; $p=100 if($p>100); 
    return ($p,$vt,$vyes,$vno,$vc); }
  
  my %tot;  
  for my $c (qw(pub cull drop)) {
    my($p,$vt,$vyes,$vno,$vpoor,$nloci);
    ($p,$vt,$vyes,$vno)= pctof($c,'type','coding','noncode');
    $tot{$c}{n_loci}= $nloci= $vt;  
    $tot{$c}{p_pcloci}= $p;
    $tot{$c}{n_pcloci}= $vyes;
    $tot{$c}{n_ncloci}= $vno;  # NA_n_ncloci
      # n_teloci from what? type tecoding ?
      
    ($p,$vt,$vyes,$vpoor,$vno)= pctof($c,'goodmap', 2, 1, 0);
    $tot{$c}{n_anymap}= $vt; 
    $tot{$c}{p_goodmap}= $p;
    $tot{$c}{n_goodmap}= $vyes;
    $tot{$c}{n_partmap}= $vpoor;  # NA_n_partmap
    $tot{$c}{n_nomap}= $vno;

    ($p,$vt,$vyes,$vno)= pctof($c,'hasho',1,0);
    $tot{$c}{t_hasho}= $vt; 
    $tot{$c}{p_hasho}= $p;
    $tot{$c}{n_hasho}= $vyes;
    $tot{$c}{n_noho} = $vno;

    ($p,$vt,$vyes,$vno)= pctof($c,'hasin',1,0);
    $tot{$c}{t_hasint}= $vt; 
    $tot{$c}{p_hasint}= $p;
    $tot{$c}{n_hasint}= $vyes; # CULL NA_n_hasint, was n_hasin
    $tot{$c}{n_noint}= $vno;
    $tot{$c}{p_noint}= 100 - $p;
    $tot{$c}{n_noint_hasho}= $sum{$c}{hasho}{noint}||0;
 
    ($p,$vt,$vyes,$vno)= pctof($c,'aafull',1,0);
    $tot{$c}{n_aafull}= $vyes; 
    $tot{$c}{n_aapart}= $vno; 
    $tot{$c}{p_aafull}=$p;
   
    my @gid= sort keys %{$sum{$c}{altof}};
    my $naloc= @gid;
    my ($nalts,@nalt)=(0);
    # 1+altof count main tr also here
    for my $g (@gid) { my $na= 1 + $sum{$c}{altof}{$g}; $nalts+=$na; push @nalt,$na; }
    @nalt= sort {$b <=> $a} @nalt; 
    $tot{$c}{n_altloci}= $naloc;
    $tot{$c}{p_altloci}= pct($naloc,$nloci); # ($nloci>0)? int(100*$naloc/$nloci): 0;
    $tot{$c}{n_alts}= $nalts;
    $tot{$c}{max_alts}= $nalt[0];
    $tot{$c}{med_alts}= $nalt[ int($naloc/2) ];
    $tot{$c}{ave_alts}= ($naloc>0)? int(10*$nalts/$naloc)/10 : 0;
    $tot{$c}{n10_alts}= scalar( grep{ $_ >= 10 } @nalt);
    $tot{$c}{n50_alts}= scalar( grep{ $_ >= 50 } @nalt);
  }
  
  #upd1905: $IDPREFIX $ORGANISM and all from %settings
  my $tb=$TABLE_G1; 
  $tb=~s/ORGANISM/$ORGANISM/g if($ORGANISM);
  $tb=~s/IDPREFIX/$IDPREFIX/g if($IDPREFIX);
  my @intab= split"\n",$tb; # old: $TABLE_G1;
  my @otab=();
  for (@intab) {
    next if(/^#/);
    my @k= m/\b([A-Z][0-9A-Z_]+)\b/g;
    for my $k (@k) { 
      next unless($k=~/_/);
      my $kl=lc($k); 
      my $c= ($kl=~s/^(cull)//)?$1:"pub";
      my $v= $tot{$c}{$kl}; $v="NA_$kl" unless(defined $v); 
      s/\b$k\b/$v/; 
    }
    push @otab, $_;
  }
  my $otab= join("\n",@otab);
  return($otab);
}

BEGIN{

our @TABLE_G1_KEYS= qw(AVE_ALTS
MAX_ALTS MED_ALTS N_AAFULL N_AAPART N_ALTLOCI N_ALTS N_ANTISENSE N_GOODMAP N_HASHO N_LOCI
N_NCLOCI N_NOINT N_NOINT_HASHO N_NOMAP N_PARTMAP N_PCLOCI N_SPLITMAP N_TELOCI N10_ALTS
N50_ALTS NCULL_ALTS NCULL_HASHINT NCULL_HASHO NCULL_LOCI NCULL_NOHO NCULL_NOINT
NCULL_NOMAP NCULL_PARTMAP P_AAFULL P_ALTLOCI P_HASHO P_NOINT P_PCLOCI T_HASHO V_IDPREFIX
V_SPECIES);


our $TABLE_G1=<<"EOT";

Table G1. ORGANISM gene set numbers, version IDPREFIX 
---------------------------------------------------
N_LOCI gene loci, all supported by RNA-seq, most also have protein homology evidence
  N_PCLOCI (P_PCLOCI%) are protein coding, N_NCLOCI are non-coding
   N_TELOCI are protein coding, expressed, loci with transposon domains

All genes (100%) are assembled from RNA evidence, 0 are genome-modelled 
# NNNNN (PP%) of loci have uniquely mapped RNA-seq, 100% have full-coverage RNA expression evidence
# NNNNN of loci have >= 90% read coverage, NNNNN have >= 60%, NNNNN have < 20% read cover. 
#
# NNNNN (PP%) of loci, and NNNNNN (PP%) of transcripts have RNA-evidence introns, 
#   of NNNNNN with model introns and genome mapping.
# NNNNNN (PP.P%) of NNNNNN RNA-evidence introns are recovered in gene transcripts.

N_HASHO/T_HASHO (P_HASHO%) have protein homology to other species genes. 
#  NNNNN (PP%) have conserved coding seq across species (sig. Ka/Ks)
#  NNNNN (PP%) are orthologs to other species (OrthoMCL clustering), 
#   NNNN (PP%) are inparalogs of orthologs (NNNNN,PP% ortholog or inparalog, of pc loci),  
#   NNNN (PP%) have ortholog only in other clade spp.
#   NNNN (PP%) have non-clade ortholog
#  NNNNN (PP%) are species-unique (in OrthoMCL clusters)
#   NNNN of these have significant other-species homology
#   NNNN of these lack other-species homology, NNNN of these have XXX paralog

N_ALTS alternate transcripts are at N_ALTLOCI (P_ALTLOCI%) loci, with MED_ALTS median, AVE_ALTS ave, transcripts per locus,
  with MAX_ALTS alts maximum, N50_ALTS loci have 50+ alts, N10_ALTS have 10+ alts, 

N_AAFULL (P_AAFULL%) have complete proteins, N_AAPART have partial proteins, of N_PCLOCI coding genes

N_GOODMAP (P_GOODMAP%) are properly mapped to chromosome assembly (>=80% align), 
   N_PARTMAP partial-mapped coverage ( 10% < align <80%), 
   N_NOMAP are ~un-mapped genes ( align < 10% ), 
# [ N_SPLITMAP loci have split-scaffold mapping ]

N_NOINT/N_GOODMAP (P_NOINT%) are single-exon loci of those mapping >= 50% to genome,
  N_NOINT_HASHO of these have homology to other species genes.
# N_ANTISENSE loci have antisense or mixed sense mapping (intron strand is opposite of coding seq)

CULLN_LOCI are culled loci, not in public gene set, but with some unique sequences. 
  CULLN_HASINT culls are multi-exon, well aligned; CULLN_NOINT are single exon, well aligned, 
  CULLN_PARTMAP are parially mapped, and CULLN_NOMAP are poorly aligned to chromosomes.
  CULLN_HASHO culls have protein homology, CULLN_NOHO lack it. 
CULLN_ALTS are culled alternate transcripts, at both public and culled loci, redundant
  in splicing patterns to public alternates, or lacking in alignment or evidence, 
  though differing somewhat by sequence alignment.

Gene locus IDs: IDPREFIXm000001t1 ..IDPREFIXm0NNNNNt1,
Alternate transcripts have ID suffix t2 .. t100.
EVm000001 is the longest protein, larger ID numbers are mostly shorter, but not for all.
--------------------------------------------------------
Gene set classification notes:

Culled transcripts are  classified as unique by transcriptome alignments, but
re-classified as redundant, or lacking sufficient evidence, by chromosome alignments.
These are separated from the public gene set as redundant or low quality, 
but are available as valid evidence transcripts, for instance by reclassifying with an 
updated chromosome assembly.

EOT

}

1;

__END__

=item evg17pig4wc genesum try1

  evg17pig4wc/publicset/pigevg4wd.genesum.txt
  
  Table G1. ORGANISM gene set numbers, version IDPREFIX 
  ---------------------------------------------------
  39840 gene loci, all supported by RNA-seq, most also have protein homology evidence
    39840 (100%) are protein coding, NA_n_ncloci are non-coding     <*FIXED, n_ncloci =  0
     NA_n_teloci are transposon protein coding, expressed, loci     <?? drop for now?
  
  All genes (100%) are assembled from RNA evidence, 0 are genome-modelled 
  
  25374/39840 (64%) have protein homology to other species genes.   <*CHECK this, eg for locus t1..9
  
  361899 alternate transcripts are at 30468 (76%) loci, with 4 median, 11.8 ave, transcripts per locus,
    with 979 alts maximum, 1297 loci have 50+ alts, 9109 have 10+ alts,   
    ^^** NEED CULL alts w/ same splice patt, 1exons; 900 alts is out of bounds, altpars or/and 1exon/same-splice alts
  
  27442 (69%) have complete proteins, 12398 have partial proteins, of 39840 coding genes
  
  37819 (95%) are properly mapped to chromosome assembly (>=80% align), 
     1174 partial-mapped coverage ( 10% < align <80%), 
     847 are ~un-mapped genes ( align < 10% ), 
  
  6639/37819 (18%) are single-exon loci of those mapping >= 50% to genome,
    3249 of these have homology to other species genes.
  
  97733 are culled loci, not in public gene set, but with some unique sequences. 
    NA_n_hasint culls are multi-exon, well aligned; 89304 are single exon, well aligned,  <*FIXED NA_n_hasint
    1265 are parially mapped, and 3901 are poorly aligned to chromosomes.
    16500 culls have protein homology, 81233 lack it. 
  129718 are culled alternate transcripts, at both public and culled loci, redundant
    in splicing patterns to public alternates, or lacking in alignment or evidence, 
    though differing somewhat by sequence alignment.
  
  Gene locus IDs: IDPREFIXm000001t1 ..IDPREFIXm0NNNNNt1,
  Alternate transcripts have ID suffix t2 .. t100.
  EVm000001 is the longest protein, larger ID numbers are mostly shorter, but not for all.
  --------------------------------------------------------
  Gene set classification notes:
  
  Culled transcripts are  classified as unique by transcriptome alignments, but
  re-classified as redundant, or lacking sufficient evidence, by chromosome alignments.
  These are separated from the public gene set as redundant or low quality, 
  but are available as valid evidence transcripts, for instance by reclassifying with an 
  updated chromosome assembly.

=cut

=item table_g1

TABLE G1. SPECIES gene set numbers, version IDPREFIX 
---------------------------------------------------
NNNNN gene loci, all supported by RNA-seq, most also have protein homology evidence
  NNNNN (PP%) are protein coding, NNNNN are non-coding
   NNNN are transposon protein coding, expressed, loci 

All genes (100%) are assembled from RNA evidence, 0 are genome-modelled 
# NNNNN (PP%) of loci have uniquely mapped RNA-seq, 100% have full-coverage RNA expression evidence
# NNNNN of loci have >= 90% read coverage, NNNNN have >= 60%, NNNNN have < 20% read cover. 

NNNNN (PP%) of loci, and NNNNNN (PP%) of transcripts have RNA-evidence introns, 
   of NNNNNN with model introns and genome mapping.
NNNNNN (PP.P%) of NNNNNN RNA-evidence introns are recovered in gene transcripts.

NNNNN/NNNNN (PP%) have protein homology to other species genes. 
#  NNNNN (PP%) have conserved coding seq across species (sig. Ka/Ks)
#  NNNNN (PP%) are orthologs to other species (OrthoMCL clustering), 
#   NNNN (PP%) are inparalogs of orthologs (NNNNN,PP% ortholog or inparalog, of pc loci),  
#   NNNN (PP%) have ortholog only in other clade spp.
#   NNNN (PP%) have non-clade ortholog
#  NNNNN (PP%) are species-unique (in OrthoMCL clusters)
#   NNNN of these have significant other-species homology
#   NNNN of these lack other-species homology, NNNN of these have XXX paralog

NNNNN alternate transcripts are at NNNNN (PP%) loci, with N median, N.N ave, transcripts per locus,
  with NNN alts maximum, NN loci have NN+ alts, NNNN have NN+ alts, 

NNNNN (PP%) have complete proteins, NNNNN have partial proteins, of NNNNN PC genes

NNNNN (PP%) are properly mapped to chromosome assembly (>=80% align), 
   NNNN partial-mapped coverage ( 10% < align <80%), 
   NNNN are ~un-mapped genes ( align < 10% ), 
# [ NNN loci have split-scaffold mapping ]

NNNNN/NNNNN (PP%) are single-exon loci of those mapping >= 50% to genome,
  9539 of these have homology to other species genes.
# NNNN loci have antisense or mixed sense mapping (intron strand is opposite of coding seq)
# ^not many of these, leave out

NNNNN are culled loci, not in public gene set, but with some unique sequences. 
  NNNNN culls are multi-exon, well aligned; NNNN are single exon, well aligned, 
  NNNNN are parially mapped, and NNNN are poorly aligned to GCR10 chromosomes.
  NNNNN culls have protein homology, NNNNN lack it. 
NNNNNN are culled alternate transcripts, at both public and culled loci, redundant
  in splicing patterns to public alternates, or lacking in alignment or evidence, 
  though differing somewhat by sequence alignment.

Gene locus IDs: IDPREFIXm000001t1 ..IDPREFIXm0NNNNNt1,
Alternate transcripts have ID suffix t2 .. t100.
EVm000001 is the longest protein, larger ID numbers are mostly shorter, but not for all.
--------------------------------------------------------
Gene set classification notes:

Culled transcripts are  classified as unique by transcriptome alignments, but
re-classified as redundant, or lacking sufficient evidence, by chromosome alignments.
These are separated from the public gene set as redundant or low quality, 
but are available as valid evidence transcripts, for instance by reclassifying with an 
updated chromosome assembly.

var keys:
AVE_ALTS MAX_ALTS MED_ALTS N_AAFULL N_AAPART N_ALTLOCI N_ALTS N_ANTISENSE
N_GOODMAP N_HASHO N_LOCI N_NCLOCI N_NOINT N_NOINT_HASHO N_NOMAP N_PARTMAP
N_PCLOCI N_SPLITMAP N_TELOCI N10_ALTS N50_ALTS NCULL_ALTS NCULL_HASHINT
NCULL_HASHO NCULL_LOCI NCULL_NOHO NCULL_NOINT NCULL_NOMAP NCULL_PARTMAP
P_AAFULL P_ALTLOCI P_HASHO P_NOINT P_PCLOCI T_HASHO V_IDPREFIX V_SPECIES

=cut

=item table_g2

TABLE G2. SPECIES gene set equivalents, Evigene v NCBI
---------------------------------------------------
# NNNNN/NNNNN of loci are also  in one or both of ncbi and ensembl gene sets
#  NNNNN are in ncbi, NNNNN are in ensembl 
#  of NNNNN in neither, NNNNN have other clade gene homology, NNNN have human gene match,

#.. counts using eqgene gff overlap differ
NNNNN/NNNNN (PP%) of loci have ncbi equivalent, NNNNN have ensembl equiv.
  of ~NNNNN without  equivalent, NNNN are multi-exon/well-mapped, NNNN are one-exon/well-mapped, 
     and NNNN are poorly mapped to chr. assembly
  NNNN of NNNNN have validated RNA introns (most of multi-exon set)
 NNNNN of NNNNN have protein homology to other clade and/or human genes
   of NNNNN with homology, intron validation is missing (NNNN are 1-exon)
   of NNNN without homology, NNNN have valid introns 
   NNN lack either primary evidence type, but are RNA assemblies that include a mix of 
      secondary gene evidence (paralogy, weak homology)

=cut

=item table_g3

TABLE G3.  SPECIES gene sets compared for gene evidence recovery

Conserved genes in SPECIES gene sets (REFDATA)
Gene set    Align   Compl Frag Miss
------------------------------------
Evigene     443.1   2572    5    9    
NCBI        433.8   2554   13   19    
Ensembl     426.8   2510   47   29    
------------------------------------

  Reference Human (n=20191 loci, Homo_sapiens)
Gene set    Found   Align   Frag  Best
Evigene     87.5%   90.8%   0.5%   39%
NCBI        86.9%   89.3%   1.2%   10%  51%  equal
Ensembl     86.3%   88.4%   2.2%   10%
---------------------------------------

=cut

=item table_g4

Intron recovery for  gene sets (n=316113 of RNA-seq mapped to chrs)
             RNA-Introns
 Geneset      Found%  RI%/mRNA   
 DrEvigene17   81.6   94.5  
 DrNCBI16      76.4   84.7  
 DrEnsembl17   57.6   67.1  
 ----------------------------------

=cut


=item evgene_stat_methods

Protein Homology Methods: 
  BLASTP -query built_genes.aa -db reference_species.aa -evalue 1e-5
  BUSCO.py -i built_genes.aa -l vertebrata_odb9 -m prot
Protein Statistics: 
  Found  = percent of reference genes with signif. align in target gene set
  AlignF = average % alignment to found reference genes (align-aa/ref-aa)
  Frag   = fragment target genes with size < 50% of reference length, of found genes
  Best   = which target set has longest alignment per ref gene, of found genes

Intron Methods 
  map RNA-seq (Illumina) to chromosome assembly with GSNAP, 
  extract splice-mapped reads and their intron locations, 
  tabulate gene-exon x rna-intron matches.
Intron Statistics
  InFound% = percent of all valid introns recovered b/n gene exons
  RI%/mRNA = percent of recoverable introns, those found by any gene set, 
  as found in mRNA transcripts of each gene set
  GeneTr  = gene transcripts total in gene set
========================================================================

=cut