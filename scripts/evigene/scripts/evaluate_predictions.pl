#!/usr/bin/env perl
# evaluate_predictions.pl 

=item about

 evigene/evaluate_predictions.pl : from evalpred.sh

 4.1. evidence annotate by base overlaps: protein, est, rnaseq, intron, tiletar 
 redo scoring: -pct 10 when using markbase; otherwise lose partial real scores  
 2010.10: add -strand overlapfilter options for most; note many ESTs not stranded;
   PASA has bogus +strand for unknown

=item todo

  want new script to set up evidence, using/creating evigene.conf
  evigene/evidence_setup.pl
    - ask for evidence types (est, intron, prot, rnaseq, refgenes, misc/transposons, tile/tar)
    - ask for source files/urls 
    - maybe make intron.gff from est.gff, rnaseq.gff
    - make evd_uniq.gff from evd.gff
    - make all_evd_exons.gff from all evd_uniq.gff
    - write parts of evigene.conf
    
=item usage

  scripts/evaluate_predictions.pl -conf=genes/evigene_eval.conf -ho='ho2,ho3' \
    genes/ogs12.gff.gz genes/nvit2_mix?asm1*.gff.gz > & eval_ogs1mix67.out &
  
  scripts/evaluate_predictions.pl -conf=genes/evigene_eval.conf -ho='ho2,ho3' \
    genes/all.*.augmap.gff.gz > & eval_augrnatile7.out &
  
  cat eval_ogs1mix67.out eval_augrnatile7.out | scripts/evaluate_predictions.pl -tab

=item output

  Gene Evidence Summary      ----------- Gene Models --------------------------------
  Evid.   Nevd    Statistic  ogs12   rs009r1 rs016r1 rss3t3  ti.afem ti.amal nvit2x6 
  EST     34.8Mb  poverbase  0.441   0.753   0.750   0.735   0.755   0.763   0.860   
  Pro     56.7Mb  poverbase  0.413   0.435   0.431   0.426   0.426   0.426   0.518   
  RNA     50.5Mb  poverbase  0.316   0.531   0.517   0.512   0.487   0.492   0.532   

=cut

use constant VERSION => '2013.08.31'; # ... way back 

use FindBin;
use lib ("$FindBin::Bin"); # assume evigene/scripts/ layout: main, prot/, rnaseq/, ...
our $EVIGENES="$FindBin::Bin";  

use strict;
use Getopt::Long;

my $pcto= 50;  # use config  
my $norun = 1;   # OPTION
my $verbose=1;   # OPTION
my $DIGITS=0; # or 1
my $NCRNA_CDS_RATIO = 0.33;
my $NCRNA_CDS_MIN =  120;
my $ALT_TRID_PATT =  't([2-9]|\d\d+)$'; ## 't[2-9][0-9]*$'; # ugh bad pat

use constant KB => 1000; # or 1024?
use constant MB => KB * KB;
use constant GB => MB * MB;
use constant KEEPCOUNTS => 0; # table extra output $ENV{count}||0;
 
my $overlapfilter="$EVIGENES/overlapfilter";
my $overlapeval="$EVIGENES/overlapeval.pl";
my $overbestgenes="$EVIGENES/overbestgene2.perl"; # not used here
my $overgenedup="$EVIGENES/overgenedup.pl"; #  
my %programs= (); # add to other evigene set

# my $KeyHOMOLOGY="ho3";  #? merge w/ evalkeys
# my $KeyHO2="ho2";
my $KeyHOMOLOG="homolog";  # make config OPTION
my $KeyPARALOG="paralog";
my $KeyINSPLIT="insplit";

my @evalkeys= qw(est pro rseq ref tar terepeat); # OPTION
my @moreevalkeys= qw(allevd cdna_eval progene_eval homology_eval); 

my (@annotkeys,@moreannotkeys); # annotate_predictions.pl config
my $genescoredir;

## change to read these filepaths from evigene.conf file;
## one conf for both eval_pred and annot_pred

my %evidence= (
  est => "est/est_uniq.gff.gz",
  pro => "prot/protein_uniq.gff.gz",
  rseq => "rnas/rnaseq_uniq.gff.gz",
  #? intr => "intron/introns.gff.gz",
  terepeat => "misc/transposon.gff.gz",
  # pasa => "est/pasa_assemblies.gff.gz",
  tar => "evidence/tar.exons_uniq.gff",
  # tar => "tiles/tilemax.gff.gz",
  'ref' => "refseq/refseq-exons_uniq.gff.gz",
  allevd => "est/all_evd_exons.gff.gz",  # for evaluate only
  cdna_eval => "est/pasa_genes.gff.gz",  # for evaluate only
  progene_eval => "prot/protein_uniq.gff.gz", #? all or uniq here
);

my %evaluate_options = (
  est => "overlapfilter -strand -pass 'exon,HSP' -pct $pcto -act keep -base",
  pro => "overlapfilter -strand -pass CDS -pct $pcto  -act keep -base",
  rseq =>"overlapfilter -strand -pass exon -pct $pcto  -act keep -base",
  intr => "overlapfilter  -intron2splice -pass 'exon,intron' -act keep -midtype scoresum", #??
  terepeat => "overlapfilter  -strand -pass 'exon,transposon' -pct $pcto -act keep -base",
  pasa => "overlapfilter  -nostrand -pass 'exon,cDNA_match' -pct $pcto  -act keep -base",
  tar => "overlapfilter -pass 'exon,ep' -pct $pcto -sumbase -act keep -base",
  'ref' => "overlapfilter -strand -pass 'exon' -pct $pcto  -act keep -base",
  allevd => "-strand -pass exon -pct $pcto -act keep -base", # -over $allevdfile -in $pred 
  cdna_eval =>"-strand -pass exon -pct $pcto", # opts for overlapeval
  progene_eval =>"-strand -pass CDS -pct $pcto", # opts for overlapeval
);

my %annotate_options = ();  ## not  same as evaluate_options
my %general_config = ();
my %public_options = ();
my @pquantiles= (0, 0.001, 0.33, 0.66, 0.95); #? change 0.001 to 0.01 == any support?

my $tabonly=0;
my $modelstats=undef;
my $evidstats=undef;
my $config;
my @configadd;
my $keepgeneset="";

my $optok= GetOptions(
  "homology=s", \$KeyHOMOLOG,  # add to config
  "config=s", \$config,
  "DIGITS=i", \$DIGITS, # drop
  "tabulateonly!", \$tabonly, 
  "modelstats!", \$modelstats, 
  "evidstats!", \$evidstats, 
  "verbose!", \$verbose, 
  "norun|n", \$norun, 
  "cadd=s", \@configadd,
  "keepgeneset=s", \$keepgeneset, # out filter for -tab
  );
# ALT_TRID_PATT # put in config

die "usage: evaluate -conf=evigene.conf  geneset1.gff geneset2.gff ...
  opts: -verbose -tabulateonly (tabulate an input set of eval.out) \n"
  unless($optok and ($tabonly or scalar(@ARGV)));

my @predlist= @ARGV;  # remaining args

evigene_config($config, \@configadd); # always even if $config null

$keepgeneset =~ s/[\s,]+/\|/g if($keepgeneset);

my %evalok =  map{ $_ => 1 } map{ split/[,;\|\s]/; } ( @evalkeys, @moreevalkeys );

## move these to config:
$KeyHOMOLOG= $general_config{'KeyHOMOLOG'} if(($general_config{'KeyHOMOLOG'}));
# my $homolog_max= $general_config{'homolog_max'} || 0;
# my $homolog_db= $general_config{'homolog_db'} || "";

##($KeyHOMOLOGY,$KeyHO2)=split",",$KeyHOMOLOGY if($KeyHOMOLOGY =~ /,/);
($KeyHOMOLOG,$KeyPARALOG)=split",",$KeyHOMOLOG if($KeyHOMOLOG =~ /,/);

$evidstats=1 unless(defined $evidstats);
$modelstats=1 unless(defined $modelstats);

my %gene_cover; # change to 
my %gene_stats;
my @result; # result tab
sub resout { push @result, @_; print @_; } # put to @result
sub verbose { print @_,"\n" if $verbose; } # should use resout instead ??
# sub docommand { verbose(@_); my $err= ($norun) ? 0 : system(@_); return $err; }

### new
# genes_stats for all
# model_stats_out

if($tabonly) {
  while(<>) { push @result, $_; } 
  eval_table(@result) if($evidstats); # option input from stdin
  model_table(@result) if($modelstats);

# } elsif($modelstats) {
# 
#   gene_stats(); # includes calc for prot_homology gene_coverage
#   model_stats_out();    
#   model_table(@result);  #??
  
} else {

  # add intron checks: insplit, inrev

  if($evidstats) {
  evid_sensitivity();     # if ??
  all_evid_specificity()  if($evalok{allevd});

  my %evnametags=( cdna_eval => "cDNA", progene_eval => "Prot", rnagene_eval => "RNA");
  # cdna_precision()        if($evalok{cdna_eval});
  # progene_precision()     if($evalok{progene_eval});
  foreach my $tag (sort keys %evnametags) {
  gene_accuracy($tag, $evnametags{$tag}) if($evalok{$tag});
  }

  }
  
  gene_stats();    # includes calc for prot_homology gene_coverage
  prot_homology()         if($evalok{homology_eval}); # part of evidstats
  gene_coverage();         # part of evidstats
  if($modelstats) {
    model_stats_out();    
  }
  
  eval_table(@result) if($evidstats); # option input from stdin
  model_table(@result) if($modelstats); #? or not
}

#..................

sub evigene_config {
  my($cfile, $addoptions)= @_;
  my $ctype=0;

  use constant{ kEVFILE => 1, kEVOPT => 2, kANOPT => 3, kEVPROG => 4, kPUBOPT => 5, kEVGENES => 6, };
  
  if($cfile) { #  and -f $cfile
    open(F,$cfile) or die "ERROR reading config: $cfile";
    my @CONFIG= <F>;
    push @CONFIG, @$addoptions if(ref $addoptions);
  
    my ($lastkey, $lastval);
    # while(<F>)   # key => value
    foreach (@CONFIG) {
      s/^\s+//;  s/\#.*$//; s/\s+$//; # end of line comments 

    ## need now to handle continuation lines, end with \
    ## or prefix with "+" << use this

      my($key,$val);
      if(/^[+\.]/ and $lastkey) { $key= $lastkey; s/^.//; $val= $lastval.$_; }
      elsif($lastval =~ m,\\$, and $lastkey) { $key= $lastkey; $lastval=~s,\\$,,; $val= $lastval.$_; }
      else { ($key,$val)= split( /[=\s]+/, $_, 2); }
      # ($key,$val)= split" ",$_,2; # space in option values : or split /[=\s]+/, $_, 2
      
      next unless($key =~ /^\w/); 
      $val =~ s/\s*$//;
      # verbose "# config k=v: $key=$val";
 			# FIXME: need FindBin / sub findevigeneapp() for these now
      if($ctype == kEVPROG and $val =~ m,^scripts,) { $val =~ s,scripts,$EVIGENES,; }

      if($key =~ /^evidence/) { $ctype= kEVFILE; }
      elsif($key =~ /^evoption/) { $ctype= kEVOPT; }
      elsif($key =~ /^anoption/) { $ctype= kANOPT; }
      elsif($key =~ /^pubopt/) { $ctype= kPUBOPT; }
      elsif($key =~ /^geneset/) { $ctype= kEVGENES; }
      elsif($key =~ /^program/) { $ctype= kEVPROG; }
      elsif($key =~ /^end$/) { $ctype= 0; }
      ##elsif($ctype == kEVFILE or $val =~ /gff/) { $evidence{$key}= $val; } 
      elsif($ctype == kEVFILE ) { $evidence{$key}= $val; } # /gff/ was bad for geneset
      elsif($ctype == 0 and $val =~ /\.gff/) { $evidence{$key}= $val; } 
      elsif($ctype == kEVOPT) { $evaluate_options{$key}= $val; } 
      elsif($ctype == kANOPT ) { $annotate_options{$key}= $val; } 
      elsif($ctype == kPUBOPT ) { } ## $public_options{$key}= $val; 
      elsif($ctype == kEVGENES) {} #{ $geneset{$key}= $val; }
      
      elsif($key =~ /^evkey|^evalkey/) {  @evalkeys= split( /[\s,;]+/, $val); }
      elsif($key =~ /^ankey|^annotkey/) {  @annotkeys= split( /[\s,;]+/, $val); }
      elsif($key =~ /^evmore|^evalmore/) {  @moreevalkeys= split( /[\s,;]+/, $val); }
      elsif($key =~ /^anmore|^annotmore/) {  @moreannotkeys= split( /[\s,;]+/, $val); }

        ## these now in %evidence hash
      elsif($key =~ /^genescore/) { $genescoredir= $val; }  # for annotation add_genescore()

      #?? change to program hash list? use eval_options to set program?
      # elsif($ctype == kEVPROG  ) 
      elsif($val =~ /overlapfilter/ or $key =~ /^overlapfilter/) { $overlapfilter=$val; }
      elsif($val =~ /overlapeval/ or $key =~ /^overlapeval/) { $overlapeval=$val; }
      elsif($val =~ /overbestgenes/ or $key =~ /^overbestgenes/) { $overbestgenes=$val; }
      elsif($val =~ /overgenedup/ or $key =~ /^overgenedup/) { $overgenedup=$val; } # replace overlapeval

      # generic keys: name date genome .. other?
      elsif($key =~ /\w/ and $val =~ /\w/) { $general_config{$key}= $val; }
        # #ALT_TRID_PATT = $general_config{'alt_transcript_id'}

      # also for now : overlap, other progs
      if($ctype == kEVPROG) { $programs{$key}= $val; }

      ($lastkey, $lastval)=($key, $val);
      
    } close(F);
  }
  
  # add overgenedup
  die "ERROR: missing overlapfilter program: $overlapfilter"
      unless( -x $overlapfilter);
  die "ERROR: missing overlapeval program: $overlapeval"
      if($overlapeval and not -x $overlapeval); # not used by annot
}



sub getEvCommand
{
  my($ev)= @_;

  my $evfile  = $evidence{$ev}; #  epasa/pasatrain_genes.best1.gff.gz 
  my $evcmd   = $evaluate_options{$ev}; # -strand -pass exon -pct 50 
  if($evfile and not -f $evfile) { print "# MISSING evidence $ev: $evfile\n";  return;  }
  verbose "# evidence-file: $evfile"; # or resout ??
  
  # foreach program ...
  $evcmd =~ s/overlapfilter/$overlapfilter/;
  $evcmd =~ s/overlapeval/$overlapeval/;
  $evcmd =~ s/overgenedup/$overgenedup/; # replace overlapeval .. dont need 3rd prog key?
  #.. ^ replace with programs{}
  # my $progpath= $programs{$progname} || $progname;

  verbose "# ev-command $evcmd"; # or resout ??
  ## caller needs to add files (prediction, evidence) to command
  return ($evcmd, $evfile);
}


sub evid_sensitivity 
{

  foreach my $ev (@evalkeys) 
  {
    next unless($evalok{$ev});
  
    resout  "Evidence: $ev\n"; # result tab
    my($evopt, $evfile)= getEvCommand($ev);
    next unless($evfile);
  
    foreach my $pred ( @predlist )
    {
    
      (my $grp=$pred) =~ s,\.gz,,; $grp =~ s,\.gff,,;  $grp =~ s,^.*/([^/]+)$,$1,;
      # my $grp= `basename $pred .gff.gz`; chomp($grp);
      my $cmd= "$evopt -in $evfile -over $pred"; # no -mark $ev
      
      ## for evaluate: want only STDERR summary output of overlapeval, not annotated gff
      # my $cmd="$overlapfilter -strand -pass "$passtype" -in $evd.gff.gz -over $pred -act keep -pct 50 -base > /dev/null
      
      resout  "Prediction: $pred\n"; # result tab
      my @res= `$cmd 2>&1 1>/dev/null`; #  # result tab stderr?
      resout  @res;
      
      resout  "-----------------------------\n";
      resout  "\n";
      
    }
    resout  "\n";
  }

}


# sub intron_sensitivity 
# {
#
# $td/overlapfilter -over genes/aphid2_aphid0-augmap.gff.gz -in intron/intron.gff.gz \
# -intron2splice=2 -pass 'exon,intron' -strand -act mark -mark ipred -baseo > /dev/null
# aphid2_aphid0:
# base statistics: overlaps n=125080 , input n=228342 ,  overset n=175846
# aphid2_epir9:
# base statistics: overlaps n=137288 , input n=228342 ,  overset n=236192
#   ... sens and spec from this ^^
#
# }



## all_evd_exons = uniq of est + rnaseq exons + prot cds
## is -sumbase right or wrong here; got bad nums with it

sub all_evid_specificity
{
  resout  "Evidence: all_evd_specif\n";
  my $ev= "allevd";
  my($evopt, $evfile)= getEvCommand($ev);
  return unless($evfile);
  
  foreach my $pred ( @predlist )
  {
    resout  "Prediction: $pred\n";
    my $cmd="$evopt -over $evfile -in $pred"; # no -mark $ev 
    #"$overlapfilter -strand -pass exon -over $allevdfile -in $pred -act keep -pct 50 -base";
    my @res= `$cmd 2>&1 1>/dev/null`; #  # result tab stderr?
    resout @res;
    resout  "-----------------------------\n";
    resout  "\n";
  }

}


=item gene_accuracy

      ## replace these stats w/ overgenedup as gene_accuracy ?
Evidence: cDNA_gene_accuracy
# evidence-file: evidence/pasa_genes.gff.gz
# ev-command scripts/overlapeval.pl -strand -pass exon -pct 50
Prediction: genes/ogs12.gff.gz
#collect_overlaps=103404
#overlaps found=31932
ngene=10194; genehit=7020; gperfect=912; nexon=43016; exonhit=31932
Exon: Sn=60.6; Sp=88.3; Sp2=30.9; Ave=74.45; hit=9133562; miss=1209208; exonbase=15059788; allover=29483158
BestTr: Sn=60.2; Sp=55.2; Sp2=30.7; Ave=57.7; hit=9071608; miss=7349429; allbase=15059788;
Gene: Sn=60.2; Sp=55.2; Sp2=30.7; Ave=57.7; hit=9071608; miss=7349429; genebase=15059788;


  cdna_eval     evidence/pasa_genes.gff.gz        # merge asmrna? for evaluate only
  rnagene_eval  evidence/rnaseq.gff.gz        # asmrna? for evaluate only
  progene_eval  evidence/protein.gff.gz      # up to prot sep11; uniq or not here?

  cdna_eval   overlapeval -strand -pass exon -pct 50 # opts for overlapeval
  progene_eval  overlapeval -strand -pass CDS -pct 50 # opts for overlapeval
TO
  cdna_eval     overgenedup -type similarCDS -mincds=10 -minutr=33 -slopexon=8  -act null -sum
  progene_eval  overgenedup -type CDSonly  -mincds=10 -slopexon=8 -act null -sum

$evigene/scripts/overgenedup.pl  -type CDSsimilar -mincds=10 -slopexon=8 -act null -sum   
# for prot only add: -exon=CDS

flamingo2.% $evigene/scripts/overgenedup.pl -type CDSsimilar -mincds=10 -slopexon=8 -act null -sum -over evidence/pas
a_genes.gff.gz -in genes/ogs12.gff.gz

#collect_overlaps ngene=10194, nexon=81089
#overlaps over=evidence/pasa_genes.gff.gz in=genes/ogs12.gff.gz genes=18941 same=6263
#overgenedup over=evidence/pasa_genes.gff.gz in=genes/ogs12.gff.gz
# novergene=10194 noverexon=81089 
# ningene=18941 genehit=6263 gperfect=2737 nexon=204904 sens=61.4 spec=33 
# identity ave. CDS=67.03, Exon=54.22, 
# identity levels: C=2594, C33=4888, C50=4028, C66=3491, C80=3142, C90=2737, I=143, X33=4761, X50=3425, X66=2268, X80=1436, X90=821

flamingo2.% $evigene/scripts/overgenedup.pl -type CDSsimilar -mincds=10 -slopexon=8 -act null -sum -over evidence/pas
a_genes.gff.gz -in genes/nvit_epi4a-augmap.gff.gz

#collect_overlaps ngene=10194, nexon=81089
#overlaps over=evidence/pasa_genes.gff.gz in=genes/nvit_epi4a-augmap.gff.gz genes=31777 same=9940
#overgenedup over=evidence/pasa_genes.gff.gz in=genes/nvit_epi4a-augmap.gff.gz
# novergene=10194 noverexon=81089 
# ningene=31777 genehit=9940 gperfect=3996 nexon=425995 sens=97.5 spec=31.2 
# identity ave. CDS=63.46, Exon=50.69, 
# identity levels: C=3985, C33=7304, C50=5939, C66=5059, C80=4550, C90=3996, I=11, X33=7023, X50=4822, X66=3126, X80=1808, X90=846


=cut

sub gene_accuracy
{
  my($ev, $evname)= @_;
  $evname ||= $ev;
  resout "Evidence: ",$evname,"_gene_accuracy\n";

  my($evopt, $evfile)= getEvCommand($ev);
  return unless($evfile);

  foreach my $pred ( @predlist ) {
    resout "Prediction: $pred\n";
    
    my $cmd;
    if($evopt =~ /overlapeval/) {
    $cmd="$evopt -in $evfile -over $pred"; # no -mark $ev 
    } else {
    $cmd="$evopt -over $evfile -in $pred";      
    }
    
    # old: $overlapeval -strand -pass CDS -pct 50 -in $evid -over $pred;
    # new: $overgenedup -type CDSsimilar -mincds=10 -slopexon=8 -act null -sum -over $evid -in pred
    resout  `$cmd`; #  # result tab stderr?
    resout "-----------------------------\n";
    resout "\n";
  }

}

## full cDNA gene accuracy test
# sub cdna_precision
# {
#   resout "Evidence: cDNA_gene_accuracy\n";
#   my $ev= "cdna_eval";
#   my($evopt, $evfile)= getEvCommand($ev);
#   return unless($evfile);
#   
#   foreach my $pred ( @predlist ) {
#     resout "Prediction: $pred\n";
#     
#     my $cmd="$evopt -in $evfile -over $pred"; # no -mark $ev 
#     # $overlapeval -strand -pass exon -pct 50 -in epasa/pasatrain_genes.best1.gff.gz -over $pred;
#     # new: $overgenedup -type CDSsimilar -mincds=10 -slopexon=8 -act null -sum
#     resout  `$cmd`; #  # result tab stderr?
#     resout "-----------------------------\n";
#     resout "\n";
#   }
# 
# }
# 
# sub progene_precision
# {
#   resout "Evidence: Prot_gene_accuracy\n";
#   my $ev= "progene_eval";
#   my($evopt, $evfile)= getEvCommand($ev);
#   return unless($evfile);
# 
#   foreach my $pred ( @predlist ) {
#     resout "Prediction: $pred\n";
#     my $cmd="$evopt -in $evfile -over $pred"; # no -mark $ev 
#     # $overlapeval -strand -pass CDS -pct 50 -in prot/protein_uniq.gff.gz -over $pred;
#     resout  `$cmd`; #  # result tab stderr?
#     resout "-----------------------------\n";
#     resout "\n";
#   }
# 
# }

##  Prot Homology scores, pctfound + bitscore, from mRNA ho3= and harabid= 
##  extend this to parse other predict evd tags (only good for annotated genes): 
#     1. sum all evd tags/{exon,cds} > pct_evd/gene
#     2. homology
#     3. gene_coverage: counts, base coverage

# sub prot_homology_OLD
# {
#   resout "Evidence: Protein_Homology\n";
#   foreach my $pred ( @predlist ) {
#     
#     my($no,$na,$sho,$sha,$ng)=(0) x 10;
#     ## merge with gene_coverage: reads same pred.gff
#     my( $gcds, $gtr, $scds, $ncds, $sexon, $nexon, $saa, $smrna, $nmrna)= (0) x 10;
#     my %medlen; 
#     map{ $medlen{$_}=[] } qw(len_protein len_exon len_cds len_transcript len_gene); # medianzero();
#    
#     my $openpipe=($pred =~ /.gz$/) ? "gunzip -c $pred |" : $pred;
#     open(G, $openpipe) or warn "#ERROR: $openpipe\n";
#     while(<G>) {
#       next unless(/^\w/);
#       my($r,$b,$e)=(split"\t")[0,3,4]; 
#       my $w= 1 + $e - $b;
#       if($w < 0) { verbose "#bad location: $pred : $_"; }
#       
#       if(/\tmRNA/) {
#         $smrna+= $w; $nmrna++;
#         push( @{ $medlen{'len_gene'} }, $w);
#         # my($gid)= m/ID=([^;\s]+)/;
#         push( @{ $medlen{'len_cds'} }, $gcds) if($gcds);
#         push( @{ $medlen{'len_transcript'} }, $gtr) if($gtr);
#         $gcds=$gtr= 0;
#         
#         my($clen,$xlen)= m/cxlen=(\d+).(\d+)/;
#         my($aa)=m/aalen=(\d+)/;  #? collect min,med,max of aa ?
#         $aa= int($clen/3) if(!$aa and $clen);
#         $saa += $aa; # not using saa
#         push( @{ $medlen{'len_protein'} }, $aa) if($aa); # only good for aalen= annotation
# 
#         #? add count of mRNA that look like ncRNA ? : cds << exon and/or <40aa
#         # check for cxlen=cds/exon key?
#         $aa||=1; # problem missing aalen= but have ho3=, for pasagenes; need to add aalen=
#         
#         ## change this, dont /$aa but use scds/3 at end: bits/aa = sum(ho) / (sum(cds)/3)
#         my($ho)=m/;$KeyHOMOLOG=([\d\.]+)/; $sho+=($ho/$aa) if $ho; $no++ if $ho; 
#         my($ha)=m/;$KeyPARALOG=([\d\.]+)/; $sha+=($ha/$aa) if $ha; $na++ if $ha;  
#         
#       } 
#       elsif(/\tCDS/) { $scds+= $w; $gcds+= $w; $ncds++;  }
#       elsif(/\texon/) { $sexon+= $w; $gtr+= $w; $nexon++;  
#         push( @{ $medlen{'len_exon'} }, $w);
#       }
#       
#     } close(G);
#     
#     my $acds= int(10*$scds/MB)/10; 
#     my $aexon=int(10*$sexon/MB)/10; 
#     
#     ## return also mean, sum from minMedianMax() == minMedMeanMaxSum() ?
#     my $xmmm = join",", minMedianMax( $medlen{'len_exon'} );  
#     my $ammm = join",", minMedianMax( $medlen{'len_protein'} );  
#     my $cmmm = SumNAveMinMedianMax( $medlen{'len_cds'} ); 
#     my $tmmm = SumNAveMinMedianMax( $medlen{'len_transcript'} ); 
#     my $gmmm = join",", minMedianMax( $medlen{'len_gene'} ); #  median2median('len_gene', 1);
#     
#     $gene_stats{'gene_cover'}{$pred}= [
#       "Coding bases: $acds $ncds   size: $ammm min,median,max\n",
#       "Exon bases: $aexon  $nexon  size: $xmmm min,median,max\n",
#       "CDS bases: $cmmm\n",
#       "Transcript bases: $tmmm\n",
#       "Gene_count: $nmrna  size: $gmmm min,median,max\n" ,
#       ]; 
# 
#     next unless($no); #?? or print 0 ?
#     
#     $na||=1; $no||=1; 
#     my $res= sprintf "protein_homol best n=%d, bits/aa=%.3f, parag n=%d, bits/aa=%.3f \n",
#       $no, $sho/$no, $na, $sha/$na; 
# 
#     resout "Prediction: $pred\n";
#     resout $res,"\n";
#     resout "-----------------------------\n";
#   }
# 
# }

sub prot_homology
{
  resout "Evidence: Protein_Homology\n";
  foreach my $pred ( @predlist ) {
    my $res= $gene_stats{'homology'}{$pred} or next;
    resout "Prediction: $pred\n";
    resout $res,"\n"; #? print if empty
    resout "-----------------------------\n";
  }
}
 
sub gene_coverage
{
  resout "Evidence: Gene_coverage\n";
  foreach my $pred ( @predlist ) {
    resout "Prediction: $pred\n";

    my $gene_cover= $gene_stats{'gene_cover'}{$pred};
        
    resout @$gene_cover if(ref $gene_cover); 
    resout "-----------------------------\n";
  }
}

##  Prot Homology scores, pctfound + bitscore, from mRNA ho3= and harabid= 
##  extend this to parse other predict evd tags (only good for annotated genes): 
#     1. sum all evd tags/{exon,cds} > pct_evd/gene
#     2. homology
#     3. gene_coverage: counts, base coverage

sub gene_stats
{

## ** NEED MORE CONFIG here
#  express = est rseq tar pasa intr
#  homology = prot ortholog paralog ref?
#  transposon = terepeat
  my $keyho= $evaluate_options{keyhomolog};
  ($KeyHOMOLOG,$KeyPARALOG)=split",",$keyho if($keyho);  
  my $keyexon = $evaluate_options{keyexon} || "est rseq ref tar pasa intr";
  my $keycds  = $evaluate_options{keycds} || "pro";
  my $keyte   = $evaluate_options{keytransposon} || "terepeat";
  map{ s/[,\s]+/\|/g } ($keyexon, $keycds, $keyte);
  my @keyexon= split /[|]/,$keyexon;
  
  my @basekeys= qw(len_protein len_exon len_cds len_transcript len_gene len_alt_transcript);
  my @supkeys = qw(pct_support pct_express pct_homology pct_ortholog pct_paralog pct_transposon );
  # my @medkeys= (@basekeys, @supkeys);

  # * add Alt-transcripts ?  now evigene.gff has no gene lines or ids, uses mRNA ID=...t[1-n] format
  # -- segregate all alt-tr, or include in base stats?
  # -- basic handling: remove alt-tr (mRNA/CDS/exon) from main stats, count only number
  $ALT_TRID_PATT = $general_config{'alt_transcript_id'} if(defined($general_config{'alt_transcript_id'}));

  my $homolog_db= $general_config{'homolog_db'} || "";
  my $homolog_max= $general_config{'homolog_max'} || 0;
    # ^ if not given, count from found hoid over all predlist.
  ## my $homolog_maxfound= 0;
  my %allhoid=();
  my $lastpred="";
  use constant CutSPLIT_C => 1; # upd1804, dang Splits      
  
  foreach my $pred ( @predlist ) {
    
    my($no,$na,$sho,$sha,$ng)=(0) x 10;
    my( $gid, $lgid, $gcds, $gtr, $scds, $ncds, $sexon, $nexon, $nmrna,
      $npartialaa, $nncrna, $pho, $pha, $aalen, $xlen, $clen, $ispartialaa,
      $isalt, $nalt,
      )= (0) x 30;
      
    my (%medlen, %evsum, %evmax, %hoid); 
    map{ $medlen{$_}=[] } (@basekeys, @supkeys);
   
    my($ismrna)=0;
    my $openpipe=($pred =~ /.gz$/) ? "gunzip -c $pred |" : $pred;
    open(G, $openpipe) or warn "#ERROR: $openpipe\n";
    while(<G>) {  
      my $atend= eof(G); #??
      
      if(/^\w/) {
      my($r,$gtyp,$b,$e)=(split"\t")[0,2,3,4]; 
      my $w= 1 + $e - $b;
      if($w < 0) { verbose "#bad location: $pred : $_"; }
      my($pid)= m/Parent=([^;\s]+)/; # exon, CDS

      $ismrna= ($gtyp =~ /RNA$/)?1:0; # or option ($gtyp =~ m/^($mrnatypes)/)
      ## from annots (if there): transposon = exon:terepeat
      if($ismrna) { # was mRNA > ncRNA, mRNA ..
        $nmrna++;
        ($gid)= m/ID=([^;\s]+)/;
        if(CutSPLIT_C) { $gid=~s/_C\d+$//; }
        
        $isalt= ($ALT_TRID_PATT and $gid =~ m/$ALT_TRID_PATT/) ? 1 : 0;
        if($isalt) {
          push( @{ $medlen{'len_alt_transcript'} }, $gtr) if($gtr);  
          $nalt++;  $nmrna--;
          next; #?
        }  
        push( @{ $medlen{'len_gene'} }, $w);
        push( @{ $medlen{'len_cds'} }, $gcds) if($gcds);
        push( @{ $medlen{'len_transcript'} }, $gtr) if($gtr);
        # $gcds=$gtr= 0; # now below
        
        ($clen,$xlen)= (m/;cxlen=(\d+).(\d+)/)?($1,$2):(0,0);  # fixme if missing
        unless($xlen) { ($xlen)= (m/\b(?:clen|trlen)=(\d+)/)?$1:0; }
        # ($aalen)=m/;aalen=(\d+)/;  #? aalen=nnn,pctaa? >> aalen=1227,78p,complete
        my($pcds,$aacomp)=(0,"");
        if(m/;aalen=(\d+),(\d+).,(\w+)/) { $aalen=$1; $pcds=$2; $aacomp=$3; 
          unless($clen) { $clen= 3*$aalen; }
        }
        else { ($aalen)= m/;aalen=(\d+)/?$1:0; }
        my ($aaseq)= m/;protein=([A-Za-z][^;\s]+)/;
        if($aaseq =~ /\w/) { $aaseq =~ s/\*$//; $aalen=length($aaseq); }        
        $aalen= int($clen/3) if( !$aalen and $clen);
        $clen=  3*$aalen if( ! $clen and $aalen>0);
        push( @{ $medlen{'len_protein'} }, $aalen) if($aalen); # only good for aalen= annotation

        # count of mRNA that look like ncRNA ? : cds << exon and/or <40aa
        ## below
        ##  my $isncrna= ($xlen>0 and ($clen/$xlen <= $NCRNA_CDS_RATIO or $clen < $NCRNA_CDS_MIN)) ? 1 : 0;
        ## $nncrna++ if ($isncrna);
        
        $ispartialaa= -1;     
        if( $aacomp ) {  $ispartialaa=($aacomp=~/partial/)?1:0; }
        elsif(  $aaseq ) {  # dont count twice? isncrna
          my $istop= index($aaseq,'*');  my $alen= length($aaseq);
          $ispartialaa= (($aaseq !~ /^M/) or ($istop >=0 and $istop < $alen-1) ) ? 1 : 0;
          ##below# $npartialaa ++ if ($ispartialaa);
        }
        
        $aalen||=1; # problem missing aalen= but have ho3=, for pasagenes; need to add aalen=

        $pho= 0;
        # my($hot)= m/;$KeyHOMOLOG=([^;,\s]+)/; # add HOid count
        my($hotd)= m/;$KeyHOMOLOG=([^;\s]+)/; # add HOid count: num,id
        my($hot,$hod)= split",",$hotd; 
        my($ho,$homax)= split"/",$hot;
        $hoid{$hod}++ if($hod); 
        # if($ho and $aalen) { $pho= ($ho/$aalen); $sho+=$pho; $no++; }
        if($ho and $homax) { $pho= ($ho/$homax); $sho+=$pho; $no++; }
        # ^^ should not switch b/n homolog/max an aalen; confuses comparisons: use one or other only

        $pha= 0;
        my($hat)= m/;$KeyPARALOG=([^;,\s]+)/; my($ha,$hamax)= split"/",$hat;
        # if($ha and $aalen) { $pha= ($ha/$aalen); $sha += $pha; $na++;}
        if($ha and $hamax) { $pha= ($ha/$hamax); $sha += $pha; $na++;}

        ## oops, some have this e format: 535/1.033e+04,AUGepir2s7967g196t1

        # ^^ is this ho/homax what we want here, or always ho/aalen ?
        #  report says bits/aa, but overbest uses ho/homax
        # ** FIXME, this/aalen doesnt work for use w/ pct_support, as paralog bits/aa > 1 is common
        # .. need another way to get percent support from bitscores here, require /hmax? for that
        # .. but use bitscore/aalen for average sho, sha scores?
        
      } elsif(/\tCDS/) { 
        next if($isalt);
        $gcds+= $w; $scds+= $w; $ncds++;  
        foreach my $k (qw(pro)) {
          $evmax{$k}{$pid} +=  $w;  
          if(/;$k=(\d+)/) { $evsum{$k}{$pid} += $1; }
          # not: $evmax{"any"}{$pid} += $w; $evsum{"any"}{$pid} += $eany;
        }
         
      } elsif(/\texon/) { 
        next if($isalt);
        $gtr+= $w; # also isalt # NO, 
        $sexon+= $w;  $nexon++;  
        push( @{ $medlen{'len_exon'} }, $w);
        my $eany= 0;
        # $evaluate_options{exonkeys} = est rseq ref tar pasa terepeat intr
        
        foreach my $k (@keyexon, $keyte) { ## CONFIG add tar,...
          $evmax{$k}{$pid} += ($k eq "intr") ? 1 :  $w; # intr: skip 1-exon genes
          if(/;$k=(\d+)/) { my $v=$1;  
            $eany=$v if($v>$eany and $k !~ /$keyte/);
            $evsum{$k}{$pid} += ($k eq "intr") ? 1 : $v;
          }
        }
        $evmax{"any"}{$pid} += $w; $evsum{"any"}{$pid} += $eany;
      }
      } # is gff
      

      if( $atend or ($ismrna and $lgid and $gid ne $lgid) # was /\tmRNA/ >> $gtyp =~ /RNA$/
        and not $isalt #??
        ) {  # DO AT END also, // and $gid ne $lgid 
          $lgid= $gid if($atend);
          
          $xlen= $gtr if($xlen == 0);
          $clen= $gcds if($clen == 0);
          $gcds=$gtr= 0;

          my $isncrna= ($xlen>0 and ($clen/$xlen <= $NCRNA_CDS_RATIO or $clen < $NCRNA_CDS_MIN)) ? 1 : 0;
          # remove good prot  from this count:
          if($isncrna and ($pho >= 0.75 or ($clen >= $NCRNA_CDS_MIN and $pha >= 0.75))) { $isncrna= 0; }
          
          $nncrna++ if ($isncrna);
          ## change to is_ncrna_or_excess_utr ... 
          ## separate?  for pcds < ncrna_cut and aalen >= 100, check homology, .., to decide?
          
          if($ispartialaa >= 0 and  not $isncrna ) {  # dont count twice? isncrna
            $npartialaa ++ if ($ispartialaa > 0);
          }

          my @k= sort keys %evmax;
          my $epro= $evsum{"pro"}{$lgid};  
          my $mpro= $evmax{"pro"}{$lgid} || 1;
          $evsum{"any"}{$lgid} = $epro if ($epro > $evsum{"any"}{$lgid});
          my ($pmax, $prna, $promax)= (0) x 10;

          ## count genes with eprot OR homology OR paralogy
          if($epro > 0 or $pho > 0 or $pha > 0) { 
            $promax= $epro/$mpro;
            #my $pho= ($ho/$aalen);
            #my $ppa= ($ha/$aalen);
            if($pho > $promax) { $promax= $pho; } 
            if($pha > $promax) { $promax= $pha; }
            $pmax= $promax if($promax > $pmax);
            push( @{ $medlen{'pct_homology'} }, $promax); 
            
            push( @{ $medlen{'pct_ortholog'} }, $pho); ## (($pho > $pha) ? $pho : 0)); #? is this right
            push( @{ $medlen{'pct_paralog'} },  $pha); ## (($pho > $pha) ? 0 : $pha)); 
            push( @{ $medlen{'pct_inparalog'} },  $pha) if($pho < $pha); ## (($pho > $pha) ? 0 : $pha)); 
            
           }
    
          foreach my $k (@k) {
            # intr: skip 1-exon genes
            my $v = $evsum{$k}{$lgid} || 0;
            my $vm= $evmax{$k}{$lgid} || 1;
            my $p= ($v / $vm);
            if($k =~ /$keyte/) {  ## CONFIG ?
            push( @{ $medlen{'pct_transposon'} }, $p); #? want this? count p > 0
            # ^^ only if no express ??
            } else {
            push( @{ $medlen{'pct_'.$k} }, $p); #? want this? count p > 0
            $pmax= $p if($p > $pmax);
            $prna= $p if($k =~ /$keyexon/ and $p > $prna); # fixme CONFIG express keys
            }
            
          }
          push( @{ $medlen{'pct_support'} }, $pmax); # FIXME: lower than pct_homology
          push( @{ $medlen{'pct_express'} }, $prna); 
       }
       $lgid= $gid if($ismrna); ##if(/\tmRNA/);
        
    } close(G);

    
    my $acds= int(10*$scds/MB)/10; 
    my $aexon=int(10*$sexon/MB)/10; 
    
    ## return also mean, sum from minMedianMax() == minMedMeanMaxSum() ?
    ## change all to SumNAveMinMedianMax()

    my @genestats=();
    foreach my $mkey (grep !/^pct_/, sort keys %medlen) { #was (@basekeys)
      my $sum = SumNAveMinMedianMax( $medlen{$mkey} ); 
      push @genestats, "$mkey: $sum\n";
    }
    
    my $fullprot= $nmrna - $nncrna - $npartialaa;
    push(@genestats, 
      "count_noncoding_gene: $nncrna  (for pCDS <= $NCRNA_CDS_RATIO or lenCDS < $NCRNA_CDS_MIN)\n",
      "count_partial_protein: $npartialaa (for missing start, internal stops)\n",
      "count_full_protein: $fullprot (complete and protein coding)\n",
      "count_alt_transcript: $nalt \n",
      );
    $gene_stats{'gene_cover'}{$pred}= \@genestats;

    my @support=();
    foreach my $mkey (grep /^pct_/, sort keys %medlen) { # (@supkeys) { ## grep /pct_/, @medkeys )
      my @qcount= quantiles( $medlen{$mkey}, \@pquantiles );
      my $sum= join ",",@qcount;
      push @support, "$mkey: $sum  for ".join("%,",map{100*$_}@pquantiles)."% \n";
    }
    $gene_stats{'evidence_support'}{$pred}= \@support;
    
#     my $xmmm = join",", minMedianMax( $medlen{'len_exon'} );  
#     my $ammm = join",", minMedianMax( $medlen{'len_protein'} );  
#     my $cmmm = SumNAveMinMedianMax( $medlen{'len_cds'} ); 
#     my $tmmm = SumNAveMinMedianMax( $medlen{'len_transcript'} ); 
#     my $gmmm = join",", minMedianMax( $medlen{'len_gene'} ); #  median2median('len_gene', 1);
#     
#     $gene_stats{'gene_cover'}{$pred}= [
#       "Coding_bases: $acds $ncds   size: $ammm min,median,max\n",
#       "Exon_bases: $aexon  $nexon  size: $xmmm min,median,max\n",
#       "CDS_bases: $cmmm\n",
#       "Transcript_bases: $tmmm\n",
#       "Gene_count: $nmrna  size: $gmmm min,median,max\n" ,
#       "ncRNA_count: $nncrna\n",
#       "partial_protein_count: $npartialaa\n",
#       ]; 

    if($no>0 or $na > 0) {  #?? or print 0 ?
      $na||=1; $no||=1; 
      my @hoid= sort keys %hoid; 
      my $nhoid= scalar(@hoid);
      map { $allhoid{$_}++ } @hoid;
      
      ##my $homolog_max= $general_config{'homolog_max'} || 0;
      ##my $homolog_db= $general_config{'homolog_db'} || "";
      # ** $homolog_max = dbsize  >> than hoid found as only best are counted. find other homax?
      # .. maybe set as total ids found for all gene sets?
      my $phoid= ($homolog_max>0) ? 100*$nhoid/$homolog_max : 0;
      #^^ not here, calc below in output after get homaxfound, but need write that to resout !
      ##  pho=%.1f%% of %d,
      
      my $res= sprintf 
      "protein_homol $KeyHOMOLOG n=%d, hoid=%d, bits/aa=%.3f, $KeyPARALOG n=%d, bits/aa=%.3f \n",
                                  $no, $nhoid, $sho/$no,       $na, $sha/$na; # , $phoid, $homolog_max
      $gene_stats{'homology'}{$pred}= $res;  # ? make list [$res] like gene_cover
    }
    $lastpred= $pred;
    
  }

  ## my $homolog_max= $general_config{'homolog_max'} || 0;
    # ^ if not given, count from found hoid over all predlist.
  my $homolog_maxfound= scalar(keys %allhoid);
  $gene_stats{'homology'}{$lastpred} =~ s/$/\nhomolog_maxfound = $homolog_maxfound\n/; # for resout
  $general_config{'homolog_maxfound'}= $homolog_maxfound;
  $general_config{'homolog_max'}= $homolog_maxfound unless($homolog_max);
  
  # return %gene_stats
}


#..... can we add equiv cDNAgene score: 100 for perfect exon matches; lower for too much/little
# epasa3/pasatrain_genes.best1.gff.gz : add above overfilt cgene score? then accum per mRNA

#..... add this homology genescore annot ; ho3= fixed in overbest; should be hbest=


sub minMedianMax  # == minMedianMaxAveSum[2..4] 
{
  my($aref)= @_;
  $aref ||= [];
  my @ma = sort{ $a <=> $b } @$aref;
  my $mid= int ( scalar(@ma) / 2);
  return @ma[0, $mid, $#ma];
}

#? want Sum,N,Ave,Min,Median,Max  order for back compat ?
sub SumNAveMinMedianMax  # nAveMinMedianMaxSum
{
  my($aref)= @_;
  $aref ||= [];
  my @ma = sort{ $a <=> $b } @$aref;
  my $n  = scalar(@ma);
  my $mid= int ( scalar(@ma) / 2);
  my $sum=0; map{ $sum+=$_ } @ma;
  my $ave= ($n>0) ? $sum/$n : 0;
  if(wantarray) {
    return ($sum, $n, $ave, @ma[0, $mid, $#ma]);
  } else {
    my $lb=" sum,n,ave,min,median,max";
    my $val= join ",", ( prbase($sum), $n, map{ prbase($_) } ($ave, @ma[0, $mid, $#ma]) );
    return $val.$lb;
  }
}

sub quantiles
{
  my($aref, $quantr)= @_;
  my @quant= @$quantr; # for aref values, cut points to count values >= cut
  my @ma = sort{ $a <=> $b } @$aref;
  my $n  = scalar(@ma);
  my @countq=();
  foreach my $i (0..$#quant) {
    my $cut= $quant[$i];
    my $j= 0; while( $j<$n and $ma[$j] < $cut) { $j++; }
    $countq[$i]= $n - $j;
  }
  return @countq; #?
}


=item gene model summary table


Gene Models
===============

Daphnia magna consensus genes (Version):
  nnnnn genes

Evidence:
  nnnnn have evidence (Protein homology, EST or RNAseq)

  nnnnn have Protein homology (e-value <= 1e-5)
  			(Parolog vs Ortholog here?)
   nnnn have only protein evidence
   nnnn have Transposon match and only protein evidence

  nnnnn have EST or Rnaseq evidence (nnnnn EST, nnnnn Rnaseq)
  nnnnn have only EST/Rna evidence

  # Percent of model supported by evidence; add later w/ homology score also
  nnnnn  >= 75% support
  nnnnn  >= 50% support
  nnnnn  >= 25% support
  nnnnn  >= 10% support

  # Additional no-evidence models from best single predictor: epir3
  # ** should check these for uniprot homology
  nnnnn mRNA genes
  

Quality:
	nnnnn are good protein coding genes (protein start/stop, >= 40aa, CDS >= UTR)
	nnnnn are ncRNA by criteria ( >= 60% UTR ) [? or < 40 aa ]
	nnnnn are incomplete protein genes (missed start, stop or internal stop)    
        Transposon here?
   
-------------------------------------------------------------

=cut

sub model_table 
{
  # my @results=@_;
  ## Note this stuff is already printed thru resout ...
  ## here we just reformat to look nicer

  my $info= join(", ", $general_config{'name'}, $general_config{'date'});
  print "\nGene Models Summary";
  print " for $info" if($info);
  print "\n";
  
  my($ev, $pr, %epct, %quals, $putpart);
  foreach (@_) {
    chomp;
  
    if(/^Gene_models: (\S+)/) { 
      $ev=$1;  
      # push(@goteval, $ev) unless( grep { $_ eq $ev } @goteval);
      # unless($tno{$ev}){ $tno{$ev}= {}; }

    } elsif(/^Gene Models Summary/) { # skip ..
      #? $ev= $pr="";
      $putpart=1;       
  
    } elsif(/^---/) { # skip ..
      $putpart=1;       
      
    } elsif(m,^Models from: (\S+),) {
      my $porig=$1; $pr=clean_predictor($porig); 
      # $pr{$pr}= $porig; 
      print ("-" x 60); print "\n";
      print " Count of genes from $porig\n";
      
    } elsif($ev and m,^(pct_\w+): (\S+)\s+for (\S+),) {
      my ($key,$val,$cuts)=($1,$2,$3);
      $epct{$key}=$val;
      
    } elsif($ev and m,^(len_\w+): (\S+)\s+(\S+),) {
      my ($key,$val,$stats)=($1,$2,$3);
      $quals{$key}=$val;
      
    } elsif($ev and m,^(count_\w+): (\S+)\s*(.*),) {
      my ($key,$val,$info)=($1,$2,$3);
      $quals{$key}="$val,$info";
    }

  if($putpart) {
    $putpart=0;       
    my($val,@val,$t, $ngene);
      
    if($ev =~ /Evidence/) {
      $t="";
      @val= split",",$epct{pct_support};
      printf "$t%6d %s\n", $val[0],"Genes (version: $pr)";
      $ngene= $val[0];
      
      $t=(" ") x 6;
      print "$t Evidence support:\n";
      $t="";
      printf "$t%6d %s\n", $val[1],"have evidence (homology, EST or RNAseq)";
      
      $t=(" ") x 4;
      @val= split",",$epct{pct_homology};
      printf "$t%6d %s\n", $val[1],"have Protein homology";
      $t=(" ") x 8;
      @val= split",",$epct{pct_ortholog};
      printf "$t%6d %s\n", $val[2],"have orthologs (>=33%)";
      @val= split",",$epct{pct_paralog};
      printf "$t%6d %s\n", $val[2],"have in-paralogs (>=33%)";

      $t=(" ") x 4;
      @val= split",",$epct{pct_express};
      printf "$t%6d %s\n", $val[1],"have Expression (EST or RNAseq)";
      $t=(" ") x 8;
      @val= split",",$epct{pct_est};
      printf "$t%6d %s\n", $val[2],"have EST (>=33%)";
      @val= split",",$epct{pct_rseq};
      printf "$t%6d %s\n", $val[2],"have RNAseq (>=33%)";
      #   nnnnn have only EST/Rna evidence

      $t="";
      @val= split",",$epct{pct_support};
      printf "$t%6d %s\n", $val[4],"have >= 95% evidence coverage";
      printf "$t%6d %s\n", $val[3],"have >= 66% evidence coverage";
      printf "$t%6d %s\n", $val[2],"have >= 33% evidence coverage";
      print "\n";
      }
      
    if($ev =~ /Quality/) {
      $t=(" ") x 6;
      print "$t Quality of models:\n";
      $t="";
      @val= split",",$quals{count_full_protein},2;
      ## add here: count_alt_transcript

      printf "$t%6d %s\n", $val[0],"are full protein genes $val[1]";
      $t=(" ") x 4;
      @val= split",",$quals{len_protein};
      printf "$t%6s,%6s,%6s %s\n", $val[4],$val[5],$val[0],"protein size (median, maximum, sum)";
      @val= split",",$quals{len_transcript};
      printf "$t%6s,%6s,%6s %s\n", $val[4],$val[5],$val[0],"transcript size (median, maximum, sum)";
      my ($trsum)= unprbase($val[0]); $trsum=1 unless($trsum > 0);
      my ($cdsum)= unprbase( split",",$quals{len_cds} );
      my $cdstr = 100 * $cdsum / $trsum;
      printf "$t%6.0f%% %s\n", $cdstr,"coding/transcript ratio";


      @val= split",",$quals{count_alt_transcript},2;
      if($val[0] > 0) {
      $t="";
      printf "$t%6d %s\n", $val[0],"are alternate transcripts to $ngene genes "; 
      $t=(" ") x 4;
      @val= split",",$quals{len_alt_transcript};
      printf "$t%6s,%6s,%6s %s\n", $val[4],$val[5],$val[0],"alt-transcript size (median, maximum, sum)";
      }
      
      $t="";
      @val= split",",$quals{count_partial_protein},2;
      printf "$t%6d %s\n", $val[0],"are partial protein genes $val[1]";
      @val= split",",$quals{count_noncoding_gene},2;
      printf "$t%6d %s\n", $val[0],"may be noncoding or aberrant models $val[1]"; 
      ##printf "$t%6d %s\n", $val[0],"are noncoding genes $val[1]"; # bad class: noncoding here
      @val= split",",$epct{pct_transposon};
      printf "$t%6d %s\n", $val[2],"have transposon match >=33% ";

      }

    $ev=""; #?? here
  } 

  }

  print ("-" x 60); print "\n";
#  print "# Predictor names:\n";
#  foreach my $pr (@pr) { my $pf=$pr{$pr}; print " $pr=$pf,"; } print "\n";

}

=item modout

Models from: genes/bestgenes_of7.an7f.gff
Gene_models: Evidence
pct_any: 34103,32275,26709,16612,5483  for 0%,0.1%,33%,66%,95% 
pct_est: 34103,24218,14885,5015,509  for 0%,0.1%,33%,66%,95% 
pct_express: 34103,28185,20508,10241,1703  for 0%,0.1%,33%,66%,95% 
pct_homology: 26660,26660,26244,24211,19135  for 0%,0.1%,33%,66%,95% 
pct_intr: 34103,7740,7594,6833,4459  for 0%,0.1%,33%,66%,95% 
pct_ortholog: 26660,13016,12191,10082,7734  for 0%,0.1%,33%,66%,95% 
pct_paralog: 26660,10556,10199,8828,7055  for 0%,0.1%,33%,66%,95% 
pct_pro: 34103,18676,18083,15502,8889  for 0%,0.1%,33%,66%,95% 
pct_rseq: 34103,23496,15700,7170,1221  for 0%,0.1%,33%,66%,95% 
pct_support: 34103,32275,27893,21603,12780  for 0%,0.1%,33%,66%,95% 
pct_transposon: 34103,1728,371,148,52  for 0%,0.1%,33%,66%,95% 
-----------------------------
Gene_models: Quality
len_protein: 10Mb,34103,297.9,1.00,150.0,11Kb sum,n,ave,min,median,max
len_exon: 54Mb,148032,361.5,3.0,205.0,14Kb sum,n,ave,min,median,max
len_cds: 31Mb,34101,895.8,3.0,453.0,32Kb sum,n,ave,min,median,max
len_transcript: 54Mb,34102,1.6Kb,31.0,1.0Kb,33Kb sum,n,ave,min,median,max
len_gene: 101Mb,34103,3.0Kb,75.0,1.4Kb,70Kb sum,n,ave,min,median,max
count_noncoding_gene: 9717  (for pCDS <= 0.33 or lenCDS < 120)
count_partial_protein: 6343 (for missing start, internal stops)
count_full_protein: 18043 (complete and protein coding)
-----------------------------

=cut

sub model_stats_out
{
  # my @basekeys= qw(len_protein len_exon len_cds len_transcript len_gene);
  # my @supkeys= qw(pct_support pct_homology pct_express  pct_transposon );

  foreach my $pred ( @predlist ) {
    resout "Models from: $pred\n";

    resout "Gene_models: Evidence\n";  
    my $evidence_support= $gene_stats{'evidence_support'}{$pred};
    resout @$evidence_support if(ref $evidence_support); 
    resout "-----------------------------\n";
    
    resout "Gene_models: Quality\n";  
    my $gene_cover= $gene_stats{'gene_cover'}{$pred};
    resout @$gene_cover if(ref $gene_cover); # some quality stats here: ncRNA, partial aa
    resout "-----------------------------\n";
    
  }
}


## from evaltab.pl
our(%tno, %pr, @goteval);

sub clean_predictor {
  my ($pr, $prgot)= @_;
  my $oldpr= $pr;
  $pr =~ s,^.*/,,;
  $pr =~ s/\.gz//; $pr =~ s/\.gff//;
  $pr =~ s/\W?augmap\S*//; $pr=~s/\.gmap\S*//; 
  $pr =~ s/(acyr|aphid|cacao|dmag|dpx|nvit)[a-z\d]*[_-]?//;  
  ## $pr =~ s/\.an\w+//;  # not always ! never?
  # $pr=$oldpr if($pr eq $prgot);
  return $pr
}

sub eval_table
{
  # my @results=@_;
  my $use_overgenedup=0;
  my($ev, $pr);
  foreach (@_) {
    chomp;
  
  if(/^Evidence: (\S+)/) { 
    $ev=$1;  
    push(@goteval, $ev) unless( grep { $_ eq $ev } @goteval);
    unless($tno{$ev}){ $tno{$ev}= {}; }

  } elsif(/^Gene Evidence Summary/) { # skip ..
    $ev= $pr="";
    # Gene Evidence Summary for nasonia_vit_genes2, 2011mar
    if( not $general_config{'name'} and m/Summary for (.*)$/) { $general_config{'name'}= $1; }
    
  } elsif(m,^Prediction: (\S+),) {
    my $porig=$1; 
    $pr= clean_predictor($porig); 
    my $i=""; my $pold; while( $pold= $pr{ $pr.$i } and $pold ne $porig ) { ++$i; } $pr.=$i;
    $pr{$pr}= $porig; 
    
  } elsif(/^# base statistics: (.+)/) { my $v=$1; 
    my @v=split /\s*,\s*/,$v; map{ my($k,$v)=split" n=",$_; $tno{$ev}{$k}{$pr}=$v; } @v;
    # base statistics: overlaps n=137288 , input n=228342 ,  overset n=236192

  } elsif(/^# ave_baseover=([\d\.]+), ave_pctover=([\d\.]+)/) { 
    $tno{$ev}{poverbase}{$pr}=$1; $tno{$ev}{poverlap}{$pr}=$2;  

  } elsif(/^# sum_basetotal=([\d\.]+)/) {
   $tno{$ev}{basetotal}{$pr}=$1; 
   # fixme for allevd/specif. need other basetotal here, but dont have in overfilter result, yet
   # Evidence: all_evd_specif

#overgenedup over=evidence/pasa_genes.gff.gz in=genes/nvit_epi4a-augmap.gff.gz
# novergene=10194 noverexon=81089 
# ningene=44040 ninexon=258538 genehit=37408 gperfect=8959 sens=1.975 spec=0.849 
# identity ave. CDS=63.46 Exon=50.69  
# identity levels: C=3985 C33=7304 C50=5939 C66=5059 C80=4550 C90=3996 I=11 X33=7023 X50=4822 X66=3126 X80=1808 X90=846
# .. spec= ghit/ingene,  sens=ghit/overgene : not so great;  sens= max(C66,X66) / genehit ?
#....... update ...........
# >> g66 change to equalcut ?
#collect_overlaps ngene=10194, nexon=81089
#overlaps over=evidence/pasa_genes.gff.gz in=genes/nvit_epi6c1-augmap.gff.gz genes=34327 same=8766
#overgenedup over=evidence/pasa_genes.gff.gz in=genes/nvit_epi6c1-augmap.gff.gz
# novergene=10194 noverexon=81089 
# ningene=34327 ninexon=412131 genehit=8766 gperfect=3578 g66=4397 sens=0.559 spec=0.128 
# identity ave. CDS=63.1 Exon=46.8 
# identity levels: C=3574 C10=8380 C33=6359 C50=5133 C66=4397 C75=4137 C90=3578 I=4 X10=8477 X33=5669 X50=3686 X66=2305 X75=1631 X90=589
# identity ovgenes: C10=9110 C33=7825 C50=6556 C66=5699 C75=5280 C90=4425

  } elsif(/^#overgenedup /){
    $use_overgenedup=1;
  } elsif($use_overgenedup and /^# novergene=(\d+)/){ #overgenedup
   $tno{$ev}{ac_novergene}{$pr}=$1;   
    
  } elsif($use_overgenedup and /^# ningene=(\d+)/ ) {  
    ##/^# ningene=(\d+) ninexon=\d+ genehit=(\d+) gperfect=(\d+) sens=(\S+) spec=(\S+)/ 
    # $tno{$ev}{ac_ngene}{$pr}=$1;  $tno{$ev}{ac_genehit}{$pr}=$2; $tno{$ev}{ac_gperfect}{$pr}=$3;
    s/^#\s+//; my %idc= map{ my($k,$v)=split"="; $k => $v } split;
    foreach (keys %idc) { $tno{$ev}{"ac_".$_}{$pr}= $idc{$_}; }

  } elsif($use_overgenedup and /^# identity ave. CDS=(\S+) Exon=(\S+)/){ #overgenedup
   $tno{$ev}{ac_avcds}{$pr}=$1;  $tno{$ev}{ac_avexon}{$pr}=$2;  
   
  } elsif($use_overgenedup and /^# identity levels:/){ #overgenedup
    s/# identity levels:\s*//; my %idc= map{ my($k,$v)=split"="; $k => $v } split;
    #.. drop this for ac_sens, ac_spec ac_equalcut
    my $c66= $idc{C66}; my $x66= $idc{X66}; $c66= $x66 if($x66>$c66);
    my $ngene= $tno{$ev}{ac_ngene}{$pr} || 1;
    my $ovgene= $tno{$ev}{ac_overgene}{$pr} || 1;
    my $p66= $c66 / $ngene;
    my $o66= $c66 / $ovgene;
    $tno{$ev}{ac_c66}{$pr}=$c66;  
    $tno{$ev}{ac_p66}{$pr}=$p66;
    $tno{$ev}{ac_o66}{$pr}=$o66;

      ## replace these stats w/ overgenedup as gene_accuracy ?
  } elsif(/^ngene=(\d+); genehit=(\d+); gperfect=(\d+)/){
   $tno{$ev}{ngene}{$pr}=$1;  $tno{$ev}{genehit}{$pr}=$2; $tno{$ev}{gperfect}{$pr}=$3;

      ## replace these stats w/ overgenedup as gene_accuracy ?
  } elsif(/^BestTr: Sn=([\d\.]+); Sp=([\d\.]+)/){
   my($sn,$sp)=($1,$2); $sn=sprintf "%.3f",$sn/100; $sp=sprintf "%.3f",$sp/100;
   $tno{$ev}{gsens}{$pr}=$sn; $tno{$ev}{gspec}{$pr}=$sp;

  ## fixme: 2nd bits.aa is missing/option
  } 
  
  # elsif(/^protein_homol (\S+) n=(\d+), bits.aa=([\d\.]+)/) # , (\S+) n=(\d+), bits.aa=([\d\.]+)
  elsif(/^protein_homol (\S+)/)  
  {
   #new: protein_homol $KeyHOMOLOG n=%d, hoid=%d, pho=%f bits/aa=%.3f, $KeyPARALOG n=%d, bits/aa=%.3
   # my($t1,$n1,$b1,$t2,$n2,$b2)=($1,$2,$3, 0,0,0);
    my($t1,$n1,$b1,$t2,$n2,$b2,$nid,$pho,$idmax)= (0) x 10;
    $t1= $1;
    $n1=$1 if(m/ n=(\d+)/); $b1=$1 if(m, bits.aa=([\d\.]+),); 
    $nid=$1 if(m/ hoid=(\d+)/);  
    $pho=$1 if(m/ pho=([\d\.]+)/); $idmax=$1 if(m/ pho=[\d\.]+. of (\d+)/); ## gone
    if(/, (\S+) n=(\d+), bits.aa=([\d\.]+)/) { ($t2,$n2,$b2)=($1,$2,$3); }
    $tno{$ev}{"ho1_n"}{$pr}=$n1; $tno{$ev}{"ho1_b"}{$pr}=$b1; $tno{$ev}{"ho1_t"}{$pr}=$t1;
    $tno{$ev}{"ho1_id"}{$pr}=$nid; 
    # $tno{$ev}{"ho1_pi"}{$pr}=$pho; $tno{$ev}{"ho1_idmax"}{$pr}=$idmax; # << drop, calc on output
    $tno{$ev}{"ho2_n"}{$pr}=$n2; $tno{$ev}{"ho2_b"}{$pr}=$b2; $tno{$ev}{"ho2_t"}{$pr}=$t2;  
  } 
  
  ##  $gene_stats{'homology'}{$lastpred} =~ s/$/\nhomolog_maxfound = $homolog_maxfound\n/; # for resout
  elsif(/^homolog_maxfound = (\d+)/) { $general_config{'homolog_maxfound'}= $1; }

  elsif(/^Coding.bases: (\S+)/) { $tno{$ev}{cbases}{$pr}=$1; }
  elsif(/^Exon.bases: (\S+)/) { $tno{$ev}{ebases}{$pr}=$1; }
  elsif(/^Gene.count: (\S+)/) { $tno{$ev}{gcount}{$pr}=$1; }
  
    ## new format for bases, counts: 11apr08
  elsif($ev and m,^(len_\w+): (\S+)\s+(\S+),) {
    my ($key,$val,$stats)=($1,$2,$3);    
    my @val= split",",$val; # sum,n,ave,min,median,max
    
    my ($ctype,$cval)= ($key =~ /len_cds/)? ("cbases",$val[0])
        : ($key =~ /len_transcript/) ? ("ebases",$val[0])
        : ($key =~ /len_gene/) ? ("gcount", $val[1])
        : ($key,$val[0]);
    $tno{$ev}{$ctype}{$pr}=$cval;
  }
# Gene_models: Quality
# len_cds: 31Mb,34101,895.8,3.0,453.0,32Kb sum,n,ave,min,median,max
# len_transcript: 54Mb,34102,1.6Kb,31.0,1.0Kb,33Kb sum,n,ave,min,median,max
# len_gene: 101Mb,34103,3.0Kb,75.0,1.4Kb,70Kb sum,n,ave,min,median,max

  }

final_putr($ev); 
}

sub unprbase {
  my $nb= shift;
  $nb =~ /^(\d+)/ or return;
  my $n= $1;  
  my @lv=(GB, MB, KB); my @lb=("Gb","Mb","Kb");
  foreach my $i (0..$#lv) { 
    if($nb =~ /$lb[$i]/i) { return $n * $lv[$i]; }
  }
  return $n;
}

sub prbase {
  my $nb= shift;
  $nb ||=0;
  my @lv=(GB, MB, KB); my @lb=("Gb","Mb","Kb");
  foreach my $i (0..$#lv) { 
    my $bv= $nb / $lv[$i];
    if($bv > 1) {
      my $dig= ($bv > 9.5) ? 0 : 1; # was $DIGITS
      return sprintf "%.${dig}f".$lb[$i], $bv; 
      }
    }
  my $dig= ($nb >= 2) ? 1 : 2; 
  return sprintf "%.${dig}f",$nb;
}

# fractions only, round to $dig || 3
sub rou { my $f=shift; my $d = shift || $DIGITS || 3; return sprintf "%.${d}f",$f; }
# if frac, rou()
##sub rou1 { my ($f)=@_; return($f >= 1.0 or $f <=  0) ? $f : rou(@_); }
sub rou1 { my ($f)=@_; return($f >= 1.2 or $f <=  0) ? $f : rou(@_); }

sub final_putr {
  my($ev)= @_;

  my $info= join(", ", $general_config{'name'}, $general_config{'date'});
  print "\nGene Evidence Summary";
  print " for $info" if($info);
  print "\n";
  
  my ($hd, @pr);
  my @prkeep= sort keys %pr;
  if($keepgeneset) { @prkeep= grep m/$keepgeneset/, @prkeep; }
  @pr= @prkeep;
  
  foreach my $ev (@goteval) {
    $tno{$ev} or next;
    my %no= %{$tno{$ev}};
    unless($hd++) {
      @pr= map{s/^\d//;$_} sort map {
          s/^ogs/0ogs/; s/^best/1best/;
          s/^(aphid0|peaaphid0)/9$1/;
          $_ } @prkeep; ## keys %pr; 
      print join("\t","Evid.", "Nevd", "Statistic", @pr),"\n"; 
      my @dash= (("------") x 2, "-------------", ("------") x scalar(@pr));
      print join("\t",@dash),"\n";
      }

    my $ni= $no{input}{$pr[0]};
    my $nb= prbase( $no{basetotal}{$pr[0]} );

    my $evfull= $ev;
    $ev =~ s,^.*/,,;
    $ev =~ s,_uniq,,;
    $ev =~ s,tarexons,TAR,;
    $ev =~ s,^all_evd_specif,Specif,; 
    #$ev =~ s,^cDNA_gene_accuracy,ESTgene,;
    #$ev =~ s,^Prot_gene_accuracy,Progene,; 
    if($ev =~ s,_gene_accuracy,gene,) { $ev=~s/cDNA/EST/; $ev=~s/Prot/Pro/; }
    $ev =~ s,^Gene_coverage,Genome,;
    $ev =~ s,^Protein_Homology,Homolog,; 
    $ev =~ s,^transposon|terepeat,T'poson,; 
    $ev =~ s,^rseq,RNA,; 
    $ev =~ s,^intr,Intron,; 
    $ev =~ s,^(est|rna|tar),\U$1,; 
    $ev =~ s,^(pro|ref),\u$1,; 
        
    if($ev =~ /Homol/){
      ## my $nmax= $no{"ho1_idmax"}{$pr[0]} || "--"; #== $homolog_max= 
      ## my $nmax= $no{"ho1_maxfound"}{$pr[0]} || "--"; #== $homolog_max= $general_config{'homolog_max'} || 0;
      my $nmax= $general_config{'homolog_max'} || $general_config{'homolog_maxfound'} || 0;

      foreach my $k (qw(ho1_n ho1_id ho1_pi ho1_b ho2_n ho2_b)) {

        if($k eq "ho1_pi" and $nmax > 0) {  ## calc this from ho1_id / maxfound
        foreach my $p (@pr) {
          my $nhoid = $no{"ho1_id"}{$p} || 0;
          ###my $phoid = ($nmax>0) ? int(1000*$nhoid/$nmax)/10 : 0;
          my $phoid = ($nmax>0) ? $nhoid/$nmax : 0;
          $no{$k}{$p}= $phoid;
          }
        }

        next unless($no{$k});
        my $nv= ($k =~ /ho1_/) ? $nmax : "--"; ## $no{ngene}{$pr[0]};
        (my $kt=$k) =~ s/_\w+/_t/;
        my $tv=$no{$kt}{$pr[0]} || "any";
        my $kv=$k;
         
        $kv =~ s,\w+_b,$tv.bits/aa,;
        $kv =~ s,\w+_n,$tv.Nmatch,; # change 2011.09, match vs found==nid
        $kv =~ s,\w+_id,$tv.Nfound,; # add config w/ total Nhoprot for Nf/Ntot = sensitiv.
        $kv =~ s,\w+_pi,$tv.\%found,; # add config w/ total Nhoprot for Nf/Ntot = sensitiv.
        
        print "$ev\t$nv\t$kv";  
        foreach my $p (@pr) { my $v=rou1($no{$k}{$p}||0); print "\t$v"; }
        print "\n";
        }

    } elsif($ev =~ /GCover|^Genome/){
      foreach my $k (qw( cbases ebases gcount )) {
        my %kv=(cbases=>"Coding Mb", ebases=>"Exon Mbase", gcount=>"Gene count");
        my $kv= $kv{$k}|| $k; my $nv="--";
        print "$ev\t$nv\t$kv";
        foreach my $p (@pr) { my $v=$no{$k}{$p}||0; print "\t$v"; }
        print "\n";
        }


      ## replace these stats w/ overgenedup as gene_accuracy ?
    #  elsif($ev =~ /ESTgene|Progene/) 
    } elsif($evfull =~ /gene_accuracy/) {
    
#    $tno{$ev}{ac_ngene}{$pr}=$1;  $tno{$ev}{ac_genehit}{$pr}=$2; $tno{$ev}{ac_gperfect}{$pr}=$3;
#    $tno{$ev}{ac_overgene}{$pr}=$1;   
#    $tno{$ev}{ac_avcds}{$pr}=$1;  $tno{$ev}{ac_avexon}{$pr}=$2;
#    $tno{$ev}{ac_c66}{$pr}=$c66;  $tno{$ev}{ac_p66}{$pr}=$p66;
# novergene=10194 noverexon=81089 
# ningene=34327 ninexon=412131 genehit=8766 gperfect=3578 g66=4397 sens=0.559 spec=0.128 

      # ac_novergene
      if($no{ac_novergene}{$pr[0]}) { # 11sep overgenedup
      # (ac_gperfect ac_genehit ac_c66 ac_o66 ac_p66 ac_avcds )
      foreach my $k (qw(ac_gperfect ac_equalcut ac_genehit ac_sens ac_spec)) { # ac_o66 ac_p66 ac_avexon
        my $nv=$no{ac_novergene}{$pr[0]};
        my %kv=( 
            ac_gperfect =>"Perfect ", 
            ac_genehit  =>"Some    ", # 
            ac_equalcut  =>"Equal66%", #  
            # ac_o66      =>"Sensitiv.", ac_p66      =>"Specific.",
            ac_sens      =>"Sensitiv.", 
            ac_spec      =>"Specific.",
            ac_avcds    =>"aveEqual", 
            ac_avexon   =>"EqualExon", #? drop, change to max(CDS,exon)?
            );
        my $kv= $kv{$k} || $k;
        print "$ev\t$nv\t$kv";  
        foreach my $p (@pr) { my $v=rou1($no{$k}{$p} || 0); print "\t$v"; }
        print "\n";
        }
        
      } elsif($no{ngene}{$pr[0]}) {
      foreach my $k (qw(gperfect genehit gsens gspec)) {
        my $nv=$no{ngene}{$pr[0]};
        my %kv=( gperfect => "Perfect ", genehit=>"Mostly  ", gsens=>"Sensitv.", gspec=>"Specifc.");
        my $kv= $kv{$k} || $k;
        print "$ev\t$nv\t$kv";  
        foreach my $p (@pr) { my $v=rou1($no{$k}{$p} || 0); print "\t$v"; }
        print "\n";
        }
      }
      
    } elsif($ev =~ /Intron/) {
        my $k= "overlaps";  my $kv="SplicesHit";
        print "$ev\t$ni\t$kv"; $ni||=1;   
        foreach my $p (@pr) { my $v=$no{$k}{$p}; my $pv=rou($v/$ni); my $x=int(1000*$v/$ni)/1000; print "\t$pv"; }
        print "\n";
        
    } else {
      ##  drop poverlap overlaps output unless requested
      my @k= (KEEPCOUNTS) ? qw(poverlap poverbase overlaps) : qw(poverbase);
      foreach my $k (@k) { # sort keys %no
        # my %kv=( poverlap => "pOver", poverbase=>"baseOv", overlaps=>"nOver");
        my $kv= ($k=~/base/) ? "BaseOverlap" : $k;
        my $nv= ($k=~/base/) ? $nb : $ni;
        ## FIXME: for all_evd_specif: want other basetotal for nv/nb: not available
        
        print "$ev\t$nv\t$kv";  # changed ni to nb here; but only for k ~ base
        foreach my $p (@pr) { my $v=$no{$k}{$p} || 0; $v=rou($v); print "\t$v"; }
        print "\n" if($ev =~ /Specif/);
        }
    }
  print "\n";
  } 
  
  print ("-" x 60); print "\n";
  print "# Predictor names:\n";
  foreach my $pr (@pr) { my $pf=$pr{$pr}; print " $pr=$pf,"; } print "\n";
  
}

