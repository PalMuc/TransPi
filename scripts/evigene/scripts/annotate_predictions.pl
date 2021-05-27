#!/usr/bin/env perl
# annotate_predictions.pl 

=item about

 evigene/annotate_predictions.pl : from anaugmap.sh

 4.1. evidence annotate by base overlaps: protein, est, rnaseq, intron, tiletar 
 redo scoring: -pct 10 when using markbase; otherwise lose partial real scores  
 2010.10: add -strand overlapfilter options for most; note many ESTs not stranded;
   PASA has bogus +strand for unknown

=item usage

  # -n dryrun to list commands
  scripts/annotate_predictions.pl -n -verbose -vers an5 -conf genes/evigene_eval.conf \
   genes/all.*.augmap.gff.gz   

  scripts/annotate_predictions.pl -vers an4 -conf genes/evigene_eval.conf \
   genes/all.rnaseq.{009r1,xu004}.augmap.gff.gz genes/all.tiling.adult.male.augmap.gff.gz \
   > &  genes/log.an4a


=item configuration

  evidence files and options are now contained in a config file.
  This is used by annotate_predictions.pl and evaluate_predictions.pl

  evigene_eval.conf is an example:
    evidence:
      est 	=   evidence/est_uniq.gff.gz
      pro 	=   evidence/protein_uniq.gff.gz
    evkeys  =   est pro rseq ref tar terepeat
    evoption:
      est  =  overlapfilter -strand -pass 'exon,HSP' -pct 50 -act keep -base
      pro  =  overlapfilter -strand -pass CDS -pct 50  -act keep -base
    ankeys =  est pro rseq ref intr pasa tar terepeat
    anmorekeys = homology_annot bestgenes 
    anoption:
      est  =  overlapfilter -strand -pass 'exon,HSP' -pct 10 -act markbase
      pro  =  overlapfilter -strand -pass CDS -pct 10 -act markidbase

=item result

  # add evidence est : evidence/est_uniq.gff.gz
  # add evidence pro : evidence/protein_uniq.gff.gz
  # add evidence rseq : evidence/rnaseq_uniq.gff.gz
  # add evidence ref : evidence/ref.exons_uniq.gff.gz
  # add evidence intr : evidence/introns.gff.gz
  # add evidence pasa : evidence/pasa_assemblies.gff.gz
  # add evidence tar : evidence/tar.exons_uniq.gff.gz
  # add evidence terepeat : evidence/transposon.gff.gz
  # annotate all.rnaseq.009r1.augmap to genes/all.rnaseq.009r1.augmap.an4.gff
    gunzip -c genes/all.rnaseq.009r1.augmap.gff.gz |\
    .. long pipe of processes .. \
    > genes/all.rnaseq.009r1.augmap.an4.gff
    
=item intron evd tests
  
  add here and/or evaluate_ and/or overbestgenes
  -- need good qual. introns (count > min), maybe max_len < 2000

  1. intron orientation flips in one gene == spurious join evidence

      overlapfilter -nostrand -pass 'intron,exon' -intron2splice=1 -act markid -mark inor -midtype=strand 
      -input $pred -over $introns 
      .. do this part in overbestgene: grep exon | perl -ne\
'($d)=m/Parent=([^;\s]+)/; ($or)=(split)[6]; ($in)=m/inor=([^;\s]+)/; $in||=0;  \
if($ld ne $d) {$v=($gf>0 and $gr>0)?"+/-":($gf>0)?"+":($gr>0)?"-":0;  $gor{$v}++; $geneor{$ld}="$v/$gf.$gr"; $gf=$gr=0; } \
if($in eq "+") { $gf++; }elsif($in eq "-"){ $gr++;} $nor{"$or=$in"}++; $ld=$d;
      .. mark pred mRNA, inrev=
      
  2. intron joining 2 genes == spurious split evidence
    ... capture introns w/ 2 geneids, mark pred mRNA insplit=
  
  overlapfilter -strand -pass 'intron,exon' -intron2splice=2 -act markid -mark gene \
  -in ../intron/intron_good.gff.gz -over dmag_ab24*.an3.gff | grep ',' | wc 

  insplit dmag_ab24 = 246
  insplit dmag2_epir3,dmag2_epir6 = 40
  insplit dmag2_epir7 = 980
  insplit dmag2_api5 = 140


=item pick best genes

  added here script for next step: pick best genes

  match ankeys and overbestgene2 -scoretype
  evigene.conf: ankeys    est pro rseq ref intr pasa tar terepeat
  
  cat \
    genes/all.rnaseq.xu004.augmap.an4.gff \
    genes/all.rnaseq.009r1.augmap.an4.gff \
    genes/all.tiling.adult.male.augmap.an4.gff \
    genes/nvit2_mix6asm1.an3.gff \
  | grep -v skip= | scripts/overbestgene2.perl -in stdin \
  -genegroup='pro,ref,rseq' \
  -scoretype='many.ho2:6,ref:6,est:4,pro:4,rseq:4,tar:3,intr:3,pasa:2,terepeat:-3,UTR:1,CDS:3' \
  -dropscore='ho2:10,ref:30,est:20,pro:30,rseq:30,tar:30,intr:1,pasa:30,terepeat:0,UTR:0,CDS:0' \
  -typeover CDS  -OVEREXON2 -genescore -trivial 1 -pctover 10 -summarize -noskip \
  >  genes/nasvit_mix4a.gff

=cut

use constant VERSION => '2013.08.31'; # ... way back 

use FindBin;
use lib ("$FindBin::Bin"); # assume evigene/scripts/ layout: main, prot/, rnaseq/, ...
our $EVIGENES="$FindBin::Bin";  

use strict;
use Getopt::Long;

my $norun = 0;   # OPTION
my $verbose=1;   # OPTION
# my $add_genescore= 1;  # make config OPTION

my $anv="an1";   # OPTION
my $KeyHOMOLOG="homolog";  # make config OPTION
my $KeyPARALOG="paralog";
my $KeyINSPLIT="insplit"; # dropping this

## my $APPoverlapfilter= findevigeneapp("overlapfilter");

my $overlapfilter="$EVIGENES/overlapfilter";
my $overlapeval=""; # not here
my $overbestgenes="$EVIGENES/overbestgene2.perl";
my $overgenedup="$EVIGENES/overgenedup.pl"; #  
my %programs= (); # add to other evigene set

## my $genedir="genes/";  # use instead input prediction path
##my $genescoredir="genes/aaeval/";  # config
##my ($blastself,$blastother); # config only?
my %insplit; 

my @annotkeys= qw(est pro rseq ref intr pasa tar terepeat); # config OPTION
my @evalkeys=();
my (@moreannotkeys,@moreevalkeys);
my %geneset;

my %evidence= (
  est => "est/est_uniq.gff.gz",
  pro => "prot/protein_uniq.gff.gz",
  rseq => "rnas/rnaseq_uniq.gff.gz",
  intr => "intron/intron.gff.gz",
  terepeat => "misc/transposon.gff.gz",
  pasa => "est/pasa_assemblies.gff.gz",
  tar => "tiles/tilemax.gff.gz",
  'ref' => "refseq/refseq-genes.gff.gz",
);

my %evaluate_options = ();
my %general_config = ();

my %annotate_options = (  #filteroptions
  est => "overlapfilter -strand -pass 'exon,HSP' -pct 10 -act markbase",
  pro => "overlapfilter -strand -pass CDS -pct 10 -act markidbase",
  rseq =>"overlapfilter -strand -pass exon -pct 10 -act markidbase",
  intr => "overlapfilter -intron2splice -pass 'exon,intron' -act markid -midtype scoresum",
  terepeat => "overlapfilter -strand -pass 'exon,transposon' -pct 10 -act markbase",
  pasa => "overlapfilter -nostrand -pass 'exon,cDNA_match' -pct 10 -act markidbase",
  tar => "overlapfilter -pass 'exon,ep' -pct 10 -sumbase -act markbase -mark tar",
  'ref' => "overlapfilter -strand -pass 'exon' -pct 10 -act markidbase",
);


my $config;
my @configadd;
my $annotate_genes = 1; # what??
my $make_bestgenes= undef;
my $rescore_only= 0;# = overbestgenes -rescore

my $optok= GetOptions(
  "homology=s", \$KeyHOMOLOG,
  "version=s", \$anv,
  "config=s", \$config,
  "annotate!", \$annotate_genes, 
  ##"bestgenes!", \$make_bestgenes, 
  "bestgenes:s", \$make_bestgenes, # arghh -nobest fails
  "rescore!", \$rescore_only, 
  "verbose|v!", \$verbose, 
  "norun|n", \$norun, 
  "cadd=s", \@configadd,
  );

my $USAGE= "usage: annotate -conf=evigene.conf geneset1.gff geneset2.gff ...
  opts: -version $anv -norun -verbose -noannotate (for bestgenes only) \n";
die $USAGE  unless($optok); #  and scalar(@ARGV) >> now config:geneset

evigene_config($config, \@configadd); # always even if $config null

# set/clear annotate_genes and make_bestgenes from config and from options; options last
my $add_genescore = grep( /^homolog|^genescore/, @moreannotkeys);
my $add_insplit = grep( /^insplit/, @moreannotkeys);

if(defined $make_bestgenes) {
  if($make_bestgenes =~ /\w/) { } # is config tag name, eg. bestrna, bestgene2
  #elsif(make_bestgenes) make_bestgenes=1 else make_bestgenes=0; ## else { $make_bestgenes = 1; } 
} else {
  $make_bestgenes= grep( /^bestgene/, @moreannotkeys);
} 
$annotate_genes= grep( /^annotate/, @moreannotkeys) unless(defined $annotate_genes);

($KeyHOMOLOG,$KeyPARALOG)=split",",$KeyHOMOLOG if($KeyHOMOLOG =~ /,/);

## replace/option: @predlist= @ARGV  with config:genesets group ?
my @predlist= @ARGV;  # remaining args
if(!@predlist and %geneset) {
  @predlist= map{ $geneset{$_} } sort keys %geneset; # from config
}
die "MISSING gene set\n".$USAGE unless(@predlist);

sub verbose { print @_,"\n" if $verbose; }
sub docommand { verbose(@_); my $err= ($norun) ? 0 : system(@_); return $err; }

my @annotpredlist=();

if($annotate_genes) {
  @annotpredlist= annotate_genes( @predlist) ;
} else {
  @annotpredlist= @predlist;
}

$make_bestgenes=0 if(@annotpredlist < 1); # was <2 why?
if($rescore_only) { # = overbestgenes -rescore
  rescore_bestgenes(@annotpredlist) ;
  
} elsif($make_bestgenes) {
  make_bestgenes($make_bestgenes, @annotpredlist); # opt to name best method
}

#..................
my %pgroup;

sub annotate_genes 
{
  my ( @predlist )= @_;

  my $annotatepipe= make_annotatefilter(); # from annotkeys
  warn "# WARNING: missing annotation commands\n" unless $annotatepipe;
  
  my @annotpredlist=();
  foreach my $pred ( @predlist )
  {
    my $genedir="";
    my $grp= $pred;
    if($grp =~ s,^(.*/)([^/]+)$,$2,) { $genedir=$1;  }
    $grp =~ s,\.gz,,; $grp =~ s,\.gff,,; 
    # $grp =~ s/\W?augmap//;  #? fixme: -augmapcf.gff now
    $grp =~ s/\W?augmap\w*//;  #? fixme: -augmapcf.gff now
    $grp =~ s,\.an\w+.*,,; #? or not; need to make sure have uniq grp names
    
    my $i=""; my $pold; 
    while( $pold= $pgroup{ $grp.$i } and $pold ne $pred ) { ++$i; } $grp.=$i;
    $pgroup{$grp}= $pred;
    
    my $anpred="${genedir}$grp.$anv.gff";
    
    if( -f $anpred ) { verbose "# already exists: $anpred"; next; }
    verbose "# annotate $grp to $anpred";
    
    my $cmd= ($pred =~ /\.gz$/)? "gunzip -c $pred" :"cat $pred";
    $cmd .= " | $annotatepipe > $anpred";
      
    my $err= docommand($cmd);
    die "ERROR: annotate $grp to $anpred: $err\n" if($err);
  
    # warn "# WARNING: now add_insplit only if add_genescore\n"
    #  if($add_insplit and not $add_genescore);
    
    #** Change again, 1. pre-make table of predictor scores, 2. add that table to pred mRNA here
    # -- table:  predID predSource  key=val key2=val key3=val ...
  
    if( $add_genescore and $evidence{'genescore'} ) {

    add_genescores($anpred, $grp); 
    
    } else {
    my $nin=0;
    $nin= add_insplit($anpred, $grp) if($add_insplit);  
    add_blastscores($anpred, $grp) if($add_genescore or ($add_insplit and $nin>0));  
    }
    
    
    push @annotpredlist, $anpred;
  }
  return @annotpredlist;
}


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
      else { ($key,$val)= split(/[=\s]+/, $_, 2); }
      # ($key,$val)= split" ",$_,2; # space in option values : or split /[=\s]+/, $_, 2
      
      ## FIXME maybe: allow $ENV{} substitutions to val ? from '\$key' ?
      ## for $EVIGENES/path apps ?
      
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
      elsif($ctype == kEVFILE ) { $evidence{$key}= $val; } # /gff/ was bad for geneset
      elsif($ctype == 0 and $val =~ /\.gff/) { $evidence{$key}= $val; } 
      elsif($ctype == kEVOPT) { $evaluate_options{$key}= $val; } 
      elsif($ctype == kANOPT ) { $annotate_options{$key}= $val; } 
      elsif($ctype == kPUBOPT ) { } ## $public_options{$key}= $val; 
      elsif($ctype == kEVGENES) { $geneset{$key}= $val; }

      elsif($key =~ /^evkey|^evalkey/) {  @evalkeys= split( /[\s,;]+/, $val); }
      elsif($key =~ /^ankey|^annotkey/) {  @annotkeys= split( /[\s,;]+/, $val); }
      elsif($key =~ /^evmore|^evalmore/) {  @moreevalkeys= split( /[\s,;]+/, $val); }
      elsif($key =~ /^anmore|^annotmore/) {  @moreannotkeys= split( /[\s,;]+/, $val); }

        ## these now in %evidence hash
      #elsif($key =~ /^genescore/) { $genescoredir= $val; }  # for annotation add_genescore()
      #elsif($key =~ /^blastself/) { $blastself= $val; }  # for annotation add_blastscores()
      #elsif($key =~ /^blastother/) { $blastother= $val; }  # for annotation add_blastscores()

      #?? change to program hash list? use eval_options to set program?
      # elsif($ctype == kEVPROG  ) 
      elsif($val =~ /overlapfilter/ or $key =~ /^overlapfilter/) { $overlapfilter=$val; }
      elsif($val =~ /overlapeval/ or $key =~ /^overlapeval/) { $overlapeval=$val; }
      elsif($val =~ /overbestgenes/ or $key =~ /^overbestgenes/) { $overbestgenes=$val; }
      elsif($val =~ /overgenedup/ or $key =~ /^overgenedup/) { $overgenedup=$val; } # replace overlapeval

      # generic keys: name date genome .. other?
      elsif($key =~ /\w/ and $val =~ /\w/) { $general_config{$key}= $val; }

			# FIXME: need FindBin / sub findevigeneapp() for these now
			# fixme2: generic version of evigene_config() is in cdna_evigenesub.pm
      # also for now : overlap, other progs
      if($ctype == kEVPROG) { $programs{$key}= $val; }
      
      ($lastkey, $lastval)=($key, $val);
      
    } close(F);
  }
  
#   die "ERROR: missing overlapfilter program: $overlapfilter"
#       unless( -x $overlapfilter); # do now in setprogram
#   die "ERROR: missing overlapeval program: $overlapeval"
#       if($overlapeval and not -x $overlapeval); # not used by annot
}


#..... can we add equiv cDNAgene score: 100 for perfect exon matches; lower for too much/little
# epasa3/pasatrain_genes.best1.gff.gz : add above overfilt cgene score? then accum per mRNA

#..... add this homology genescore annot ; ho3= fixed in overbest; should be hbest=
# replace w/ add_blastscores
# sub add_genescore
# {
#   my ($pred,$grp)= @_;
# 
#   ## split this into grp-self.genescore, grp-others.genescore ?
#   ## change to read .blastp instead of .genescore?
#   
#   my $genescorepatt= "$genescoredir/${grp}-*.genescore";
#   
#   # my $existscore = -e $genescorepatt;  ## this is bad
#   my $existscore= `ls $genescorepatt`; chomp($existscore);
#   
#   unless($existscore) {   
#     verbose "# MISSING: genescore $pred : $genescorepatt";
#   
#   } elsif ( $norun ) {
#     verbose "# genescore $pred from $genescorepatt";
# 
#   } elsif( -f $pred and $existscore ) {
#     verbose "# genescore $pred";
#     my %gb;
#     open(GSCORE, "cat $genescorepatt | sort -k1,1 -k2,2nr |"); # is this sort right?
#     while(<GSCORE>) { if(/^(\S+)\s(\S+)\s(\S+)/){ my($g,$h)=($1,"$2/$3"); $gb{$g}=$h  unless($gb{$g}); } }
#     close(GSCORE);
#     
#     my $hotemp= "$pred.ho.tmp";
#     open(GFF,$pred) or warn $pred;
#     open(OUT,">$hotemp"); # tmpfile
#     while(<GFF>) {
#       if(/\tmRNA/) { my($g)=m/ID=([^;\s]+)/; my $b=$gb{$g}; s/$/;$KeyHOMOLOG=$b/ if $b; }
#       print OUT $_;
#     }
#     close(GFF); close(OUT);
#     rename($hotemp, $pred) if( -s $hotemp );  # or docommand(mv); 
#     
#   }
# 
# }


=item  add_intronsplit

  2. intron joining 2 genes == spurious split evidence
    ... capture introns w/ 2 geneids, mark pred mRNA insplit=
  
  overlapfilter -strand -pass 'intron,exon' -intron2splice=2 -act markid -mark gene \
  -in ../intron/intron_good.gff.gz -over dmag_ab24*.an3.gff | grep ',' | wc 

  insplit dmag_ab24 = 246
  insplit dmag2_epir3,dmag2_epir6 = 40
  insplit dmag2_epir7 = 980
  insplit dmag2_api5 = 140

grep -c insplit= genes/*.an4.gff  # 2x above score, for each gene
genes/dmag2_api5-augmap.an4.gff:269
genes/dmag2_epir2-augmap.an4.gff:81
genes/dmag2_epir3-augmap.an4.gff:79
genes/dmag2_epir6-augmap.an4.gff:81
genes/dmag2_epir7-augmap.an4.gff:1617
genes/dmag_ab24augmap0.an4.gff:283
genes/dmag_ep24augmap2.an4.gff:602
genes/pasa2dm6_bestgenes.an4.gff:260

=cut

sub add_insplit # add_intronsplit
{
  my ($pred,$grp)= @_;
  
  my $nadd=0;
  %insplit=(); # clear gene ids
  my $ev="insplit";
  my $evfile  = $evidence{$ev};
  my $evcmd   = $annotate_options{$ev} or return;
  $evcmd =~ s/overlapfilter/$overlapfilter/;
  verbose "# add evidence $ev : $evfile";
  unless(-f $evfile and $evcmd) { warn "# MISSING evidence $ev: $evfile\n";  return; }
    
  # insplit = overlapfilter -strand -pass 'intron,exon' -intron2splice=2 -act markid 
  my $cmd="$evcmd -mark insplit -over $pred -in $evfile |";
  # $cmd .=" grep insplit= |";
  verbose($cmd); 
  # my $err= ($norun) ? 0 : system(@_); return $err; }
  unless ( $norun ) {  
    open(P, $cmd);
    while(<P>) {
      my($ingenes)= m/insplit=([^;\s]+)/;
      next unless($ingenes =~ /,/); # score only if 2+ genes hit same intron
      ## NOTE: alt-tr will cause problems here, unless ID is clearly annotated so we can skip IDs
      
      my ($ir,$is,$it,$ib,$ie,$iv,$io)=split"\t";
      my $inid= join(":","$iv/in",$ir,$ib,$ie,$io);
      my @ingenes= split ",", $ingenes;
      foreach my $g (@ingenes) { $insplit{$g} .= "$inid,"; $nadd++; }
      ## ^^ change value? should be score/id1,id2,.. where score == 
      ## of insplit/gene ? or intron count? expect only 1 intron split/gene?
    } close(P);
  }
  return $nadd;
}

sub add_genescores
{
  my ($pred,$grp)= @_;

  my $nadd=0;  my %gkeys=();
  my $gscore = $evidence{'genescore'};
  my @gscore= `ls $gscore`; chomp(@gscore); @gscore = grep( /$grp/, @gscore); 
  my $existscore= @gscore;

  verbose "# add_genescores $pred : @gscore";
  verbose "# genescore columns:  gene_id  gene_source  key1=val1  key2=val2 ..";
  if ( $norun ) {

  } elsif( -f $pred ) {
  
    my (%gscore, $cat, $cmd);
    if( $existscore ) {      
      $cat= ($gscore[0] =~ /\.gz/) ? "gunzip -c" : "cat";
      $cmd= join(" ",$cat, @gscore, "|");
      open(GSCORE, $cmd) or warn "# ERROR: $cmd\n"; # is it file or list? 
      while(<GSCORE>) { 
       next unless(/^\w/); 
       my($gid, $gsrc, @scores)= split;
       if($gsrc =~ /=/) { unshift(@scores, $gsrc); $gsrc=""; }
       @scores = grep /=/, @scores; # should all be key=value
       next unless(@scores);       
       if( my $lastscore= $gscore{$gid.$gsrc}) { unshift(@scores, $lastscore); }
       $gscore{$gid.$gsrc}= join(";",@scores); # ready to add to mRNA
       } close(GSCORE);
            
    } else {
      verbose "# MISSING: add_genescores $pred ";
    }
    
    unless(%gscore) {
      warn "# ERROR: add_genescores $pred : no scores to add";
      return;
    }
    
    my ($ng,$nho)=(0,0);
    my $hotemp= "$pred.gs.tmp";
    open(GFF,$pred) or warn $pred;
    open(OUT,">$hotemp") or warn $hotemp; # tmpfile
    while(<GFF>) {
      if(/\tmRNA/) { 
        my($gsrc)=(split)[1];
        my($gid)=m/ID=([^;\s]+)/;  $ng++;     
        my $s= $gscore{$gid.$gsrc} || $gscore{$gid}; # allow missing gsrc
        if( $s ) {
          # also strip old keys
          my @k= $s =~ m/(\w+)=/g; map{ $gkeys{$_}++ } @k;
          my $k= join("|",@k); 
          s/;?($k)=[^;\s]+//g;
          s/$/;$s/; $nadd++; 
          }
        }
      print OUT $_;
    }
    close(GFF); close(OUT);
    rename($hotemp, $pred) if( -s $hotemp );  # or docommand(mv); 
    verbose "# add_genescores $pred found: $nadd/$ng";
    verbose "# add_genescores $pred keys: ", join(", ", map{ "$_=$gkeys{$_}" } sort keys %gkeys);
  }

  return $nadd;
}


sub add_blastscores
{
  my ($pred,$grp)= @_;

  my $nadd=0;
  my $blastself = $evidence{'blastself'};
  my $blastother= $evidence{'blastother'};
  # do any annotate_options{blastother} work here ??
  
  my @bself= `ls $blastself`; chomp(@bself); @bself = grep( /$grp/, @bself); 
  my %bself= map{ $_,1 } @bself;
  my @bother= `ls $blastother`; chomp(@bother);   
  @bother = grep { !$bself{$_} } grep /$grp/, @bother;
  my $existscore= @bother or @bself;
  
  if ( $norun ) {
    verbose "# blastscores $pred : $KeyPARALOG=@bself, $KeyHOMOLOG=@bother";
    verbose "# ASSUMES blast table (thisid[0], otherid[1], bitscore[-1]) sorted high bitscore first";
    ## assumes blast-table format, where query == this predictor
    # Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score

  } elsif( -f $pred ) {
  
    my (%bself, %bparalog, %bother, $cat, $bsort, $cmd);
    if( $existscore ) {
      ## ASSUMES Blast output sorted by best score ?? or do we do that?
      verbose "# blastscores $pred : $KeyPARALOG=@bself, $KeyHOMOLOG=@bother";
      verbose "# ASSUMES blasttable (-m 8,9) scores sorted best first";
      
      $bsort=(@bself > 1) ? " sort -k1,1 -k12,12nr |" : "";
      $cat= ($bself[0] =~ /\.gz/) ? "gunzip -c" : "cat";
      $cmd= join(" ",$cat, @bself, "|$bsort");
      open(GSCORE, $cmd) or warn "# ERROR: $cmd\n"; # is it file or list? 
      while(<GSCORE>) { 
       next unless(/^\w/); my @v=split; my($q,$t,$bits)= @v[0,1,-1];
       if($q eq $t) { $bself{$q}= $bits unless($bself{$q}); } 
       else { $bparalog{$q}="$bits,$t" unless($bparalog{$q}); }
       } close(GSCORE);
      
      # ** if list of files, need to sort by score **
      $bsort=(@bother > 1) ? " sort -k1,1 -k12,12nr |" : "";
      $cat= ($bother[0] =~ /\.gz/) ? "gunzip -c" : "cat";
      $cmd= join(" ",$cat, @bother, "|$bsort");
      open(GSCORE, $cmd) or warn "# ERROR: $cmd\n"; # is it file or list? 
      while(<GSCORE>) { 
       next unless(/^\w/); my @v=split; my($q,$t,$bits)= @v[0,1,-1];
       if($q ne $t) { $bother{$q}="$bits,$t" unless($bother{$q}); }
       } close(GSCORE);
      
    } else {
      verbose "# MISSING: blastscores $pred : $blastother, $blastself";
    }
    
    unless(%bother or %bparalog or %insplit) {
      warn "# ERROR: blast $pred : no scores to add";
      return;
    }
    
    my ($ng,$nho)=(0,0);
    my $hotemp= "$pred.ho.tmp";
    open(GFF,$pred) or warn $pred;
    open(OUT,">$hotemp") or warn $hotemp; # tmpfile
    while(<GFF>) {
      if(/\tmRNA/) { 
        my($g)=m/ID=([^;\s]+)/;  $ng++;
        
        my $s=$bself{$g};    
        my $b=$bother{$g};   
        my $p=$bparalog{$g}; 
        if($s or $b or $p) {
        map{ if(/e\+/){ $_= int($_) } } ($s,$b,$p); # 1.033e+04
        ## FIXME: these can be e-notation:  1.23e+4 > convert to digits
        
        # overbestgenes expects format  "ho=score/maxscore,.."
        if($b) { $b =~ s=,=/$s,= if($s); s,$,;$KeyHOMOLOG=$b,;  $nho++;}
        if($p) { $p =~ s=,=/$s,= if($s); s,$,;$KeyPARALOG=$p,; }
        if($s and ($b or $p)) {  
          $b =~ s,/.*,,; $p =~ s,/.*,,; $s =~ s,/.*,,;
          my($bb,$bl)=($p>$b)? ($p,"pa") : ($b,"ho"); $bb=int(100*$bb/$s); 
          s/$/;pHOBEST=$bb\%$bl/; }
        }
        
        my $insplit= $insplit{$g};  s/$/;$KeyINSPLIT=$insplit/ if $insplit; 
        $nadd++;
        }
      print OUT $_;
    }
    close(GFF); close(OUT);
    rename($hotemp, $pred) if( -s $hotemp );  # or docommand(mv); 
    verbose "# blastscores $pred found: $nho/$ng";
  }

  return $nadd;
}

sub setprogram {
  my($evcmd, $defname, $defprog)= @_;
  my $progna=(split " ",$evcmd)[0]; 
  my $prog= $defprog;
  if($progna and $prog= $programs{ $progna }) {  }
  elsif($defname and $defprog) { $progna=$defname; $prog= $defprog;  }
  $evcmd =~ s/$progna/$prog/;
  unless( -x $prog) { die "ERROR: missing program $progna: $prog\n"; } # or die
  return $evcmd;
}

sub make_annotatefilter {
  my @cmd;
  my @usekey; # make public
  foreach my $ev (@annotkeys) {
    my $evfile  = $evidence{$ev};
    my $evcmd   = $annotate_options{$ev} or next;
    
    $evcmd= setprogram($evcmd, "overlapfilter", $overlapfilter);
    # $evcmd =~ s/overlapfilter/$overlapfilter/;
    # not here# $evcmd =~ s/overlapeval/$overlapeval/;
    verbose "# add evidence $ev : $evfile";
    unless(-f $evfile) { warn "# MISSING evidence $ev: $evfile\n";  next; }
    
    my $cmd="$evcmd -mark $ev -in stdin -over $evfile";
    push @cmd, $cmd;
    push @usekey, $ev;
  }
  return unless(@cmd);
  
  # strip off old keys before adding again
  my $keypatt= join"|",@usekey;
  unshift @cmd, "perl -pe 's/;($keypatt)=[^;\\s]+//g;'";

  my $cmd= join(" | ",@cmd);
  return $cmd;
}

sub make_bestgenes {
  my ($evtag, @annotpreds)= @_; # call w/ @annotpredlist
  my ($evcmd);
  if($evtag =~ /\D/) { $evcmd = $annotate_options{$evtag}; } # FIXME: option to use other tags, bestrna, bestgenes2, ...
  else { $evtag= "bestgenes"; $evcmd = $annotate_options{$evtag}; }
  unless($evcmd) { warn "# MISSING bestgenes command for $evtag\n"; return; }
  
  $evcmd= setprogram($evcmd, "overbestgenes", $overbestgenes);
  # $evcmd =~ s/overbestgenes/$overbestgenes/;

  ## mustkeepdrop  genes/mustkeepdrop.list         # for bestgenes, expert selections
  my $mustfile  = $evidence{"mustkeepdrop"};
  if( $mustfile and -f $mustfile and not $evcmd =~ m/ \-must/) {
    $evcmd =~ s/$/ -must=$mustfile/; #? do here or leave to .conf bestgenes ?
  }
  
  # -- @predlist is wrong list: need output annot list 
  my $pred  = $annotpreds[0];
  my $npred = @annotpreds;
  # my $grp= `basename $pred .gff.gz`; chomp($grp);
  # my $genedir= `dirname $pred`; chomp($genedir);

  # my ($genedir) = $pred =~ m,^(.*)/[^/]+$,;
  my $genedir="";
  if($pred =~ m,^(.*/)[^/]+$,) {  $genedir=$1;  }
  my $bestgenes="${genedir}bestgenes_of$npred.$anv.gff"; # option to name.
  
  verbose "# bestgenes of $npred to $bestgenes";
  if( -f $bestgenes ) { verbose "# SKIP: already exists: $bestgenes"; return; }

  my $cmd="";
  $cmd= ($pred =~ /\.gz/) ? "gunzip -c " : "cat ";
  $cmd .= join " ",@annotpreds;
  $cmd .= " | $evcmd";
  # $cmd .= " -in stdin > $bestgenes";
  $cmd .= " -in stdin -out $bestgenes"; # log?

  my $err= docommand($cmd);
  warn "# ERROR: make_bestgenes $bestgenes.gff : $err\n" if($err);
}

sub rescore_bestgenes {
  my @annotpreds= @_;  
  my $ev= "rescoregenes"; # tag option??
  #? should use all bestgenes options but add -rescore flag ?
  my $evcmd = $annotate_options{$ev} or return;
  $evcmd= setprogram($evcmd, "rescoregenes", $overbestgenes);
  # $evcmd =~ s/rescoregenes/$overbestgenes/;

  my $npred = @annotpreds;
  verbose "# rescoregenes for $npred";
  
  foreach my $pred (@annotpreds) {
    my $genedir=""; my $pname=$pred;
    # if($pred =~ m,^(.*/)[^/]+$,) {  $genedir=$1;  }
    if($pred =~ m,^(.*/)([^/]+)$,) {  $genedir=$1; $pname=$2; $pname =~ s/.gff.*//; }
    $pname =~ s/.$anv//;
    my $bestgenes="${genedir}$pname.$anv-sc.gff";
    if( -f $bestgenes ) { verbose "# SKIP: already exists: $bestgenes"; return; }
    my $cmd = "$evcmd -in $pred > $bestgenes";
    my $err= docommand($cmd);
    warn "# ERROR: rescore_genes $bestgenes.gff : $err\n" if($err);
  }
}

=item bestgenes call

  cat \
    genes/all.rnaseq.xu004.augmap.an4.gff \
    genes/all.rnaseq.009r1.augmap.an4.gff \
    genes/all.tiling.adult.male.augmap.an4.gff \
    genes/nvit2_mix6asm1.an3.gff \
  | grep -v skip= | scripts/overbestgene2.perl -in stdin \
  -genegroup='pro,ref,rseq' \
  -scoretype='many.ho2:6,ref:6,est:4,pro:4,rseq:4,tar:3,intr:3,pasa:2,terepeat:-3,UTR:1,CDS:3' \
  -dropscore='ho2:10,ref:30,est:20,pro:30,rseq:30,tar:30,intr:1,pasa:30,terepeat:0,UTR:0,CDS:0' \
  -typeover CDS  -OVEREXON2 -genescore -trivial 1 -pctover 10 -summarize -noskip \
  >  genes/nasvit_mix4a.gff

=cut


=item intron-complete model scoring test

  see now at  evigene/scripts/intronscore.pl

  
## test new intron-complete scoring for gene models:
set gns=genes/bestgenes_cdsrna.mix7k

scripts/overlapfilter.perl -intron2splice=over -pass 'exon,intron' -act markid -mark inids \
-over intron/intron_good.gff.gz -in $gns.gff | grep -v '^#' | perl -ne \
'if(/\tmRNA/){ ($g)=m/ID=([^;]+)/; print "$lg\t$gin{$lg}\n" if($gin{$lg}); ($ins)=/;insplit=([^;\s]+)/; ($inx)=m/;inexon=([^;\s]+)/; \
$xi=0; $gin{$g}="inexon=$inx"; $gin{$g}.=";insplit=$ins" if($ins); $lg=$g; } \
elsif(/\texon/) { ($intr)=m/;intr=([^;\s]+)/; ($inid)=m/;inids=([^;\s]+)/; $xi++; $gin{$g}.=";x$xi"; \
$gin{$g}.=":intr=$intr,$inid" if($intr or $inid); } END{ print "$lg\t$gin{$lg}\n" if($gin{$lg}); }' \
 >  $gns.inset

cat $gns.inset | perl -ne\
'chomp; ($g,@x)=split/[\t;]/; ($in,$nx)=m/inexon=(\d+).(\d+)/; $good=$bad=$join=0; \
if($in == 0) { $gin{$g}="none"; } elsif($in < 0) { $gin{$g}="err.rev"; } \
else { %lid=(); @join=(0) x $nx; $ix=0; \
foreach $x (@x) { $x=~m/x(\d+)/ or next; $ix++;  $j=0; my %id=(); \
if( ($v)= $x=~/intr=(.*)/ ){ ($iv,@id)=split",",$v; $iv=~s,/.*,,; \
$good++ if($iv>0); $bad++ if($iv < 0); map{ $j=1 if($lid{$_}); $id{$_}=1; } @id; }\
$join[$ix-1]=$j; $join++ if($j); %lid=%id;  }\
if($bad == 0 and $good == $nx and $join + 1 == $nx) { $gin{$g}="complete"; } \
elsif( ($join > $nx-3 or $join > 0.75 * $nx) and $good > $nx-2 and $bad == 0 ) \
{ ($j1,$j2)=(2,$nx-2); if($join > 8) { $j1++; $j2--; } \
if( $nx < 4 or grep { $_ == 0 } @join[$j1..$j2]) { $gin{$g}="miss_inner";} else { $gin{$g}="complete_inner"; } } \
elsif($good and $join and $bad == 0) { $gin{$g}= ($in/$nx > 0.66) ? "good" : "ok"; }  else { $gin{$g}="poor"; } \
$join .= ",".join"",@join; } print "$g\t$gin{$g}\t$good/$bad/$join\t@x\n"; ' \



# if long run of joined introns allow miss 2 at ends for complete_inner ?
# e.g. mk7AUGepir3p1s10g4t1  good    23/0/22,0001111111111111111111111       inexon=23/25 x1 x2 x3:.. joined thru x25
#  change from 'good' to 'ok' if few joins
# e.g. mk7AUGepir3p1s10g9t1    good    6/0/3,000010000100001000        inexon=6/18  << ok not good, use 6/18 < 0.66 ?

mk7AUGepir10s102g9t1    good.2/1/0      inexon=2/4 x1 x2:intr=12,N31699 x3:intr=12,N31699 x4
mk7AUGepir9s102g15t1    err.1/0/0       inexon=1/3 x1:intr=10,N31701 x2 x3

                                          v-- flag not complete, inner intron missed
mk7AUGepi5s102g27t1     complete.inner011101111111      inexon=10/12 x1:intr=96,N31749 x2:intr=244,N31751,N31749 x3:int
r=258,N31753,N31751 x4:intr=175,N31755,N31753 x5:intr=-65/+34,N31757 x6:intr=133,N31759,N31757 x7:intr=178,N31761,N3175
9 x8:intr=136,N31763,N31761 x9:intr=123,N31765,N31763 x10:intr=164,N31767,N31765 x11:intr=146,N31769,N31767 x12:intr=48
,N31769

mk7AUGepi5s102g28t1     complete        inexon=8/8 x1:intr=68,N31771 x2:intr=102,N31771,N31773 x3:intr=92,N31773,N31775
 x4:intr=111,N31775,N31777 x5:intr=95,N31777,N31779 x6:intr=78,N31779,N31781 x7:intr=90,N31781,N31783 x8:intr=54,N31783
mk7PASAgasmbl_4608      complete        inexon=7/7 insplit=83/in:Scaffold102:585343:612584:+,83/in:Scaffold102:585343:6
12584:+,36/in:Scaffold102:609583:612584:+, x1:intr=56,N31785 x2:intr=118,N31785,N31787 x3:intr=210,N31787,N31871,N31873
,N31875 x4:intr=33,N31873,N31877 x5:intr=226,N31875,N31877,N31879 x6:intr=174,N31879,N31881 x7:intr=119,N31881,N31883
mk7AUGepir9s102g59t1    good.2/1/0      inexon=2/9 x1 x2 x3 x4 x5 x6 x7 x8:intr=74,N31885,N31887 x9:intr=43,N31887
mk7PASAgasmbl_4614      complete        inexon=2/2 insplit=83/in:Scaffold102:585343:612584:+,36/in:Scaffold102:609583:6
12584:+,36/in:Scaffold102:609583:612584:+, x1:intr=36,N31883 x2:intr=119,N31881,N31883
mk7AUGepir9s102g92t1    complete.inner01111111111       inexon=11/12 x1 x2:intr=200,N32015 x3:intr=905,N32015,N32017 x4
:intr=1112,N32017,N32019 x5:intr=876,N32019,N32021 x6:intr=764,N32021,N32023 x7:intr=828,N32023,N32025 x8:intr=1178,N32
025,N32027 x9:intr=1684,N32027,N32029 x10:intr=1547,N32029,N32031 x11:intr=1264,N32031,N32033 x12:intr=756,N32033
mk7AUGepir10s102g55t1   complete        inexon=16/16 x1:intr=136,N31927 x2:intr=262,N31929,N31927 x3:intr=252,N31931,N3
1929 x4:intr=204,N31933,N31931 x5:intr=175,N31935,N31933 x6:intr=205,N31937,N31935 x7:intr=172,N31939,N31937 x8:intr=10
6,N31941,N31939 x9:intr=104,N31943,N31941 x10:intr=96,N31945,N31943 x11:intr=80,N31947,N31945 x12:intr=88,N31949,N31947
 x13:intr=96,N31951,N31949 x14:intr=91,N31953,N31951 x15:intr=99,N31955,N31953 x16:intr=62,N31955

                                                     v--- not complete, missed inner join
mk7AUGepir3s102g52t1    complete.inner011111111111111011111111111       inexon=27/27 x1:intr=21,N31957 x2:intr=36,N3195
7,N31959 x3:intr=21,N31959,N31961 x4:intr=16,N31961,N31963 x5:intr=17,N31963,N31965 x6:intr=33,N31965,N31967,N31969 x7:
intr=14,N31967,N31971 x8:intr=9,N31971,N31973 x9:intr=36,N31969,N31973,N31975,N31977 x10:intr=20,N31975,N31979 x11:intr
=26,N31979,N31981 x12:intr=41,N31977,N31981,N31983 x13:intr=25,N31983,N31985 x14:intr=11,N31985,N31987 x15:intr=6,N3198
7 x16:intr=8,N31989 x17:intr=24,N31989,N31991 x18:intr=28,N31991,N31993 x19:intr=33,N31993,N31995 x20:intr=35,N31995,N3
1997 x21:intr=34,N31997,N31999 x22:intr=62,N31999,N32001 x23:intr=76,N32001,N32003 x24:intr=46,N32003,N32005 x25:intr=3
4,N32005,N32007 x26:intr=58,N32007,N32009 x27:intr=36,N32009

mk7AUGepi5s102g65t1     complete.inner011111111 inexon=9/10 x1 x2:intr=+62/-16,N31893 x3:intr=105,N31893,N31895 x4:intr
=91,N31895,N31897 x5:intr=94,N31897,N31899 x6:intr=117,N31899,N31901 x7:intr=106,N31901,N31903 x8:intr=112,N31903,N3190
5 x9:intr=155,N31905,N31907 x10:intr=78,N31907

mk7AUGepir16bs102g41t1  complete.inner0111111111        inexon=10/11 x1 x2:intr=16,N31909 x3:intr=34,N31909,N31911 x4:i
ntr=64,N31911,N31913 x5:intr=116,N31913,N31915 x6:intr=148,N31915,N31917 x7:intr=147,N31917,N31919 x8:intr=155,N31919,N


=cut
