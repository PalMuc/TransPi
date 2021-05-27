#!/usr/bin/perl
# overbestgene2.perl

=item overbestgene2
  
  special case of gene gff overlap detection:
   - input one gene/cds gff file with quality scores
      and many overlapping predictions (e.g. exonerate predicts)
   - output best by score non-overlapping subset genes + cds
     (some overlap possible, want all mostly distinct predictions)
     
 gzcat exonerate-gldmoj.gff3.gz | $td/overbestgene2.perl -in stdin > exonerate-gldmoj-best.gff &
 #collect_gff=1660943, ngene=763282   + 6735 geneids with terepeat > about same...
 #done kept=23000, skipped=744001
    vs kept=15038, skipped=748244 for simple filter (includes terepeat filter) 
    
=item utr update

  need better aberrant-utr scoring
  -- options: 1. average/median utr bases, 2. score utr:weight, 3? count introns/utr, >2 is aberrant ?
     score how close utr size is to median bases,
       +weight for == median, -weight for distance from median?
     
=cut

use strict;
use warnings;
use Getopt::Long;

use constant VERSION => '2018.01.08'; 
# '2013.09.29'; # 27
# '2012.08.20'; << sourceweight/weightsrc ignored sort bug; and pCDS/UTR adjusted so pCDS > 0.95 ok
# "2011.11.03"; # "2011.10.28"; # "2011.04.24"; #"2010.10.03";

use constant UPD1801 => 1; # fixes for intron annots: hasKeyNINTRON, nintron & intr, sort of same, check for both

my $BINSIZE  = 200000; # was 99000 ; #was# 1000; want large steps for this case?
my $ONEBIN_GENES = 1; # is this working or not? THIS IS a problem at BINSIZE
      # 2013sep: problems from huge-genespans picking exons that overlap later bins.
      #  switch from genebins{ref}{bin},rloc to exonbins{ref}{bin}{gid}  for each exon ?


my $CHECKJOINS = 1; #  works, default now but option
#old# my $JOINSCORES = "CDS,ho3,pro,xde"; # FIXME: option ; CDS must be 1st, need ordered
my $JOINSCORES = "CDS,homolog,ovpro,ovrna"; # FIXME: option ; CDS must be 1st, need ordered
my $KeyINSTRAND = "intr"; # special exon attr for evidence intron orients (+/-/0) change to +1,-1,0 ?        
my $KeyNINTRON = "nintron";  # score like inexon=validin/nexon; 
  #  niscore = validin/nexon * nexon * 100 == 100 * validin; ?? to make base-comparable
my $KeyJOINCHECK= "joinck"; # new special annot field for joinerr score, +100..0..-100, <<0 means likely join; <0 is poor; >0 is good

my $DO_OVEREXON2 = 1; #? make default?
my $UTR_OVERLAP_MAX = 0.20; # some portion 0.1 - 0.5 to allow two models w/ UTR overlap ; DO_OVEREXON2
my $UTR_CDS_TEST = 0.20; # some portion 0.1 - 0.5 to allow two models w/ UTR overlap ; DO_OVEREXON2
my $SKIP_REVCDS = 0; #upd18.02: not sure want this, to skip antisense with cds overlap. some invalid, some valid

my $RESCORE_ONLY = 0;
my $VECSCORE=0; #  upd1801: =0 default; 
my $doalttr = 0; # add "2011.02.11"
my $EXPAND_LOCUS  = 10000; # for locuslist; option?

# use constant UTR_USE_NEW_WEIGHT => 2;  # only way now
my $UTR_TOO_BIG = 0.60; # NOW, portion of transc len.; WAS 2, proportional to expected size, e.g 1000 bp for 400 bp expect
my $utr_expected_size = 0; # NOT USED NOW ? but set UTR_TOO_BIG
  # see http://arthropods.eugenes.org/arthropods/summaries/arthropod-genestruc-table.pdf

# modified gff gene [rloc] record fields
use constant{ jCHR => 0, jSRC => 1, jTYPE => 2, jBEGIN => 3, jEND => 4, 
              jSCORE => 5, jSTRAND => 6, jPHASE => 7, jATTR => 8, jGID => 9, 
              jCDSLEN => 10, jTRLEN => 11, jSCORESUM => 12 }; # last 3 for mRNA only

my $sumgenescore=1; # now default
use constant { SCORE_SUM => 1, SCORE_BEST => 0 };
my $COMBINE_SCORE = SCORE_SUM;  
  # using SCORE_SUM now as default; works, simple and consistant w/ scoresum 
  # best == weighted sum of evidence bases per gene

use constant oSTRANDED => 1; # FIXME1705: CDS-overlap1 should also use oSTRANDED, unless 1-cds-exon

use constant { 
kSKIP_OVERLAP => 1, 
kSKIP_TRIVIAL_EVD => 2, # later skim out paralog-only
kSKIP_BAD => 3, # if flagged in attrib : not same as skip required ?
kDEMOTE_MAYBEJOIN => 56, # 
kKEEP_REQUIRED => 77,
kKEEP_ALTERNATE => 69,
kSKIP_REQUIRED => 666,   # from user input list
kSKIP_NOEVD => 999,  # fails to pass dropscore
};

my $debug=1;
my ($showskips,$typeover,$overlaps,$overlaplist,%geneexons,
    $input,$output,$itype,$action,$actid,$ok,$mark,$ordersource,$weightsrc,
    $mustkeepdrop,$locuslist);
my @input;

my $exontypes="CDS,exon"; #,match_part,HSP";
my $mrnatypes="mRNA"; #,match"; # fixme ?? need this
$typeover= "CDS"; # make default?

my $badattr="";
my $pctover= 0.1; # should set sensible default 0.1 ?
my ($overbases,$overlen,$trivialover,@overids); # globals set in overlap1()
my ( $scoretype,$genegroup, @genegroup, $dropscore, @dropscore, @droptype); 
use constant USE_DROPSCORE_ANDOR => 1;
use constant { kDROP_OR => 0, kDROP_AND => 1, kDROP_NOT => -1, kDROP_MUSTNOT => 2 }; 

$scoretype= "score";
$genegroup= $dropscore="";

# special scoretypes calc here:  CDS (size), UTR (CDS/exon ratio), nintron, pro

my $summarize=0; #? make default
my %summary=();
## add input option: -mustkeepdrop = file of [ geneid keep/good/drop/skip ... ]

my %mustkeepgene=();
my %keepgene=(); # cannot preload this, see sub overbestgene
my %skipgene=(); # can preload this w/ must-drop ids
my($nskip, $nkeep)= (0,0);
my $OUT= *STDOUT;

my $optok= GetOptions(
  "input=s", \@input,  # or @ARGV ??
  "output=s", \$output,  
  "typeover=s", \$typeover, 
  "exontypes=s", \$exontypes, 
  "badattr=s", \$badattr, 
  "mrnatypes=s", \$mrnatypes, 
  "pctover=i", \$pctover, 
  "trivialover=i", \$trivialover, #? same param as pctover?

  "scoretype=s", \$scoretype, # e.g. 'many.pro:10,est:2,CDS:2'
  "genegroup=s", \$genegroup, # which of scoretype have genegroup ids
  "dropscore=s", \$dropscore,  # same format as scoretype many., min score to keep

  "utrsize=s", \$utr_expected_size, # change to proportion of trans. size
  "mustkeepdrop=s", \$mustkeepdrop,  

  "locuslist=s", \$locuslist,  # 2011.oct add to pick only these new loci +/- region
      # input file is gff? bed?  chr start end
  
  "alttr!", \$doalttr, 
  "sourceorder=s", \$ordersource, # srcA,srcB or srcC:99,srcD:98 == order: srcA>srcB>srcD>srcC .. other = last
  "sourceweight=s", \$weightsrc,  # 2011.9: srcC:0.9, srcD:1.1
  "skipshow!", \$showskips, 
  "genescore!", \$sumgenescore, 
  "summarize!", \$summarize, 
  "CHECKJOINS!", \$CHECKJOINS, 
  "OVEREXON2!", \$DO_OVEREXON2, 
  "REVCDSSKIP!", => \$SKIP_REVCDS, # bad name? noREVCDS != noREVCDSskip; using -skip[show]
  "JOINSCORES=s", \$JOINSCORES, 
  "COMBINE_SCORE=s", \$COMBINE_SCORE, 
  "ONEBIN_GENES!", \$ONEBIN_GENES, 
  "RESCORE_ONLY!", \$RESCORE_ONLY, 
  "vecscore!", \$VECSCORE, # PUT_SCORESUM, # or [no]vecscore ? ; needs COMBINE_SCORE/sumgenescore 
  "debug!", \$debug, 
  );

die "
usage:  perl overbestgene2 
  -exontypes $exontypes  -mrnatypes $mrnatypes  
  -scoretype=$scoretype [-genegroup=scores -dropscore=scores:val]
  -typeover=$typeover -pctover=$pctover -trivialover=$trivialover
  -skipshow -summarize -debug
  -input genes.gff|stdin > bestgenes.gff
" unless($optok);


$pctover= $trivialover if($trivialover and not $pctover);
$trivialover= $pctover if($pctover and not $trivialover); #? want or not ?
$trivialover= $trivialover/100.0 if($trivialover >= 1);
$pctover= $pctover/100.0 if($pctover >= 1);

if($utr_expected_size > 0) {
  if($utr_expected_size > 90) { $UTR_TOO_BIG= $utr_expected_size / 1000; } # bases / transcript-size
  elsif($utr_expected_size > 1) { $UTR_TOO_BIG= $utr_expected_size / 100; } # percent
  else {  $UTR_TOO_BIG= $utr_expected_size; }
}

if($typeover) { $typeover  =~ s/[,; ]+/|/g; }
$exontypes =~ s/[,; ]+/|/g;
$mrnatypes =~ s/[,; ]+/|/g;

my $dualscore=  0; ## DROP: ($scoretype =~ /(\w+)\.score$/i) ? 1: 0; #? drop this or not for many score
my $manyscore= ($scoretype =~ /^many/i) ? 1: 0;
   $manyscore= 1 if ($scoretype =~ /,\w/); # commas mean same as many
   # ^simplify, drop many prefix: -scoretype=est:5,pro:4,rna:3,...
my @scoreweight= (1) x 100; #? fixme same size as scorevec
my $overlap1_ovrev=0; my @orevids=(); # global flag for overlap1.oSTRANDED reverse overlap

# sub getScoretypes( $scoretype, \@scoreweight, \@dropscore, \@droptype, \@genegroup );
# 10/09/18: can we add -scoreweight, for terepeats? yes
if($manyscore) {
  my $i=0; 
  my @sf= manyscorefields(1);
  @scoreweight= (1) x scalar(@sf);
  #$scoretype =~ /^many.(.+)/i; 
  foreach my $sf (@sf) { ## split",",$1
    if($sf =~ m/(\w+):([\d\.\-]+)/) { $scoreweight[$i]=$2; } #? allow fractions? 0.5 ?
    $i++;
    } 
  
#  ** revise -dropscore to allow AND / OR / NOT?, e.g.
#    -dropscore= CDS:180 AND ( est:60 or pro:60 or rseq:60) AND NOT terepeat:60
## format: -dropscore=est:60,pro:60,rseq:60,+CDS:180,-terepeat:60

  @sf= manyscorefields(0);
  @dropscore= (0) x scalar(@sf);
  @droptype = (0) x scalar(@sf);
  if($dropscore) {
    $dropscore =~ s/^many\.//; 
    foreach my $i (0..$#sf) { 
      my $sf= $sf[$i]; 
      my($dcode,$dval)=(undef,undef); 
      if($dropscore =~ m/(\W?)$sf[:=](\d+)/) { ($dcode,$dval)=($1,$2); }
      elsif($dropscore =~ m/$sf[:=]([\+\*-])(\d+)/) { ($dcode,$dval)=($1,$2); } # user-error fix...
      if(defined $dval) {
        $dcode= (not defined $dcode) ? kDROP_OR 
          : ($dcode eq "+") ? kDROP_AND : ($dcode eq "-") ? kDROP_NOT
          : ($dcode eq "*") ? kDROP_MUSTNOT : kDROP_OR;
        $dropscore[$i]= $dval; 
        $droptype[$i]= $dcode; 
        }
    }
  }
  
  @genegroup= (0) x scalar(@sf);
  if($genegroup) {
    $genegroup =~ s/^many\.//; 
    foreach my $i (0..$#sf) { my $sf= $sf[$i]; 
      #? add flags to group type: priority, exon/cds span?
      if($genegroup =~ m/\b$sf[:=]([\w.]+)\b/) { $genegroup[$i]= $1; }
      elsif($genegroup =~ m/\b$sf\b/) { $genegroup[$i]= 1; }
    }
  }
  
}

my %ordersrc=();
  #? can we combine this opt as ordering (1..9) and scoreweight (0.90 = low, 1.10 = high)?
  # -sourceorder=srcA:1,srcB:2,.. 1 > 2 ; NO GOOD; use 2 opts
  #   OR srcA:0.9,srcB:1.1 1.1 weighs more than 0.9
  # .. then 1.1 > 0.9 implies same order, use highest val=1st, lowest/0 = last
  
if($ordersource) {
  my $ord=1; # need to know input set to get nsource max = first, use -value ? NO, use other opt for weight
  foreach my $os (split /[,; ]/, $ordersource) {
    if($os =~ m/(\S+):([\d\.]+)/) { $ordersrc{$1}=$2; } # 2011.9: must be numeric>0
    else { $ordersrc{$os}=$ord++; }  
    ##else { $ordersrc{$os}=$ord--; } # 2011.9 : NO GOOD; reversed ++: count from $nsource=first .. 1=last
  }
}

my %weightsrc=(); # implied weight=1 for unspecified cases; otherwise 0 < wt <= infinity?
if($weightsrc) { # 2011.9
  foreach my $os (split /[,; ]/, $weightsrc) {
    if($os =~ m/(\S+):([\d\.]+)/) { my($s,$p)=($1,$2); $weightsrc{$s}=$p if($p>0); }  #NO change, must input proportional weight 0.9,1.1,..
    ## else { $weightsrc{$os}=1; }  # probably useless, need :number input
  }
}

# sub readMustKeepDrop( $mustkeepdrop, \%mustkeepgene, \%skipgene);
if($mustkeepdrop) {
  die "ERROR: cannot read mustkeepdrop=$mustkeepdrop" unless( -f $mustkeepdrop ); 
  open(F, $mustkeepdrop);
  my $kdkeys='keep|best|good|true|drop|skip|bad|false|altmodel|alt|join';
  #** add other acts: 
  #   swapbest/swapmain=Funhe2EKm006301t6/Funhe2Exx11m012615t10 ; 
  #   partdrop=Split2/3 
  #   changesrc=gffsrc
  my %mustact;
  while(<F>){
    next unless(/^\w/);
    my($gid,$qual,@xt)=split;  
    $qual=$xt[0] if(@xt and $xt[0] =~ /^($kdkeys)/ and not($qual =~ /^($kdkeys)/)); # other table format: id oid action
    # last in only : drop gid before this?
    if($mustact{$gid}) {
      my $lasta=$mustact{$gid};
      warn "mustkeepdrop dupid: id=$gid, action=$qual, ignoring lastact=$lasta\n";
      delete $skipgene{$gid} ; delete $mustkeepgene{$gid} ;
    }
    $mustact{$gid}=$qual; 
    if($qual =~ /^(drop|skip|bad|false)/i) { $skipgene{$gid} = kSKIP_REQUIRED; }
    elsif($qual =~ /^(keep|best|good|true)/i) { $mustkeepgene{$gid} = kKEEP_REQUIRED; }
    elsif($qual =~ /^(altmodel|alt)/i) { $mustkeepgene{$gid} = kKEEP_ALTERNATE; }
    elsif($qual =~ /^(join)/i) { $mustkeepgene{$gid} = kDEMOTE_MAYBEJOIN; }
    
    ## handle altmodel kKEEP_ALTERNATE ? see doalttr:   $gref->[jATTR]=~s/$/;alttr=1/;   # flag; add ID of main tr?

    else {
      warn "mustkeepdrop ambiguous: id=$gid, action=$qual not in drop|keep\n";
    }
  } close(F); 
}


my %locuslist=(); # global
readlocuslist() if($locuslist);


#$OUT= *STDOUT;
if($output and $output =~ /\w/) {
  if(-e $output) { warn "# OUTPUT exists; moved to $output.old\n"; rename($output,"$output.old"); }
  open(OUTFILE,">$output") or die "error on out=$output\n";
  $OUT= *OUTFILE;
}  

$ok = 0;
my $inh= *STDIN;
if(@input > 1) { # change to @ARGV ?
  # do like perl: collect_gff(<>) ... how?
  die "FIXME for input list";
  
} elsif(@input == 1) {
  $input= shift @input;
  $ok = ($input =~ /.gz$/) ? open($inh,"gunzip -c $input |") 
      : ($input =~ /^(stdin|-)/) ? do { $inh= *STDIN; 1; }
      : open($inh,$input);
  die "bad -input=$input" unless($ok);
}


my ($genebins, $genes, $exons, $genesources, $geneorder) = collect_gff($inh); 

$genesources={} unless(ref $genesources);
my $ngenesources= scalar(keys %$genesources); # so that unspec wont all be last if srcB:99
my $gsources=join",", sort keys %$genesources;
$ngenesources ||=1;

#?? add any more info at top? sources? genegroups?
print $OUT "##gff-version 3\n".
"#program: overbestgenes, selection of best gene set by evidence scores\n".
"#version: ",VERSION,"\n".
"#author: d. g. gilbert, gilbertd at indiana edu\n".
"#scoretype: $scoretype\n". 
"#dropscore: $dropscore\n".
"#sources: n=$ngenesources: $gsources\n".
"\n"; 

# $geneorder valid only for $RESCORE_ONLY
if($RESCORE_ONLY) { 
  foreach my $gid (@$geneorder) {  
    putgene( $genes->{$gid}, $exons->{$gid});  
    }
    
} else { 
  overbestgenes(); 
  
  $nskip= scalar(keys %skipgene);
  if($showskips and $nskip>0) {
    print $OUT "# skipped genes ",(".") x 50," \n";
    foreach my $gid (sort keys %skipgene) {  
      my $skipv= $skipgene{$gid};
      putgene( $genes->{$gid}, undef, "skip=$skipv") if (ref $genes->{$gid});  # $exons->{$gid},
      }
  }
  
  # add report of any mustkeep not found/kept
  if( %mustkeepgene ) {
    foreach my $mkid (sort keys %mustkeepgene) {
      unless($keepgene{$mkid}) { 
      print $OUT "#Must_missing: $mkid\t",$mustkeepgene{$mkid},"\n"; 
      }
    }  
  } 
  
}

putsummary() if($summarize);

warn"#done kept=$nkeep, skipped=$nskip\n" if $debug;
# close($OUT);
# exit(0);

#...........

sub overbestgenes 
{

  foreach my $ref (sort keys %$genebins) {
    my @bins= sort{$a <=> $b} keys %{$genebins->{$ref}};
    
      ## problem w/ genes spanning bins, BINSIZE: 
      ## need to check in both bins? before keepgene/putgene
      ## .. simple answer, collect all genes per ref, est. max 5k .. forgo bins here
    
    my @sgenes=();
  
  unless($ONEBIN_GENES) {  
    my %sgenes=();
    foreach my $ib (@bins) {  
      # my @didexons=(); # should keep above bin loop? or preload w/ 2-bin kept genes
      foreach my $gref (@ { $genebins->{$ref}{$ib} }) {
        my $gid= $gref->[jGID];
        next if $skipgene{$gid}; # pre-dropped if no evidence ...
        $sgenes{$gid}= $gref;  
      }
    }
    @sgenes= values %sgenes;
    @bins=(1);
  }
  
    my @lastbinexons=(); # keep some near bin boundary so dont make duplicates at ound
    foreach my $ib (@bins)   
    { 
  if($ONEBIN_GENES) {  
      @sgenes= @ { $genebins->{$ref}{$ib} };
  }    
       # dropped @bins loop here for all genes/ref
       # this however slows processing a bunch, and eats memory
      my @didexons= @lastbinexons; # should keep above bin loop? or preload w/ 2-bin kept genes

      my @mustdemotejoin=();
      if( %mustkeepgene ) { 
        if($CHECKJOINS) {
        @mustdemotejoin= grep { $mustkeepgene{$_->[jGID]} == kDEMOTE_MAYBEJOIN }
                      grep { $mustkeepgene{$_->[jGID]} } @sgenes; 
        }
        my @mustkeep= grep { $mustkeepgene{$_->[jGID]} != kDEMOTE_MAYBEJOIN }
                      grep { $mustkeepgene{$_->[jGID]} } @sgenes; 
        if(@mustkeep or @mustdemotejoin) {
          foreach my $gref (@mustkeep,@mustdemotejoin) { 
            my($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid)= @$gref;
            my $tkeep= $mustkeepgene{$gid};
            ## upd1508: add kKEEP_SWAPMAIN : alt id > main gid ID change required ..
            if($tkeep == kKEEP_REQUIRED || $tkeep == kKEEP_ALTERNATE) {
              my $rexons= $exons->{$gid}; #? bug here
              push(@didexons, @$rexons) if($rexons);  
              unless( $keepgene{$gid} ) {
                $keepgene{$gid} = $tkeep; $nkeep++;
                $gref->[jATTR]=~s/$/;must=$tkeep/;
                $gref->[jATTR]=~s/$/;alttr=1/ if($tkeep == kKEEP_ALTERNATE);
                putgene( $gref, $rexons); 
                }
            } elsif($tkeep == kDEMOTE_MAYBEJOIN ) {
              $gref->[jATTR]=~s/$/;must=$tkeep/; #? unless($gref->[jATTR]=~/;must=/);
            } else { 
              # any other kKEEP codes?
            }
          }
        }
            
        my $nmust= @mustkeep; $nmust += @mustdemotejoin;
        my $nbefore= @sgenes;
        @sgenes= grep { ! $mustkeepgene{$_->[jGID]} } @sgenes;
        my $ncut= @sgenes;
        warn "bad mustgenes: nmust=$nmust; nbefore=$nbefore; nafter=$ncut" 
          if( $nmust + $ncut != $nbefore);
        @didexons= sort { $a->[jBEGIN] <=> $b->[jBEGIN] } @didexons; # new 2011.10 for overlap1 test
      }
       
      
      ## FIXME: 20111002 remove skips from sgenes
      @sgenes= grep { ! $skipgene{$_->[jGID]} } @sgenes;
      
      if($manyscore) { @sgenes= sort _sort_manyscoreloc @sgenes; }
      #DROP: elsif($dualscore) { @sgenes= sort _sort_refscore2loc @sgenes; }
      else { @sgenes= sort _sort_refscoreloc @sgenes; }
  
      # check gene joins/splits here?
      #  for given 1st gene at location, collect all mostly overlapped genes,
      #   if any has >=150% exons or <= 50% exons, check for join/split
      #   "best" here may be sum of evidence scores of mRNAs per source-type 
      #    over longest span window?
  
      if($CHECKJOINS) { # Hack for now, uses scores of ho2, CDS; but works
        # 2011.10: add must-joinmaybe? put these last in list
        # 2010.10: modify to use genegroup best_ attrib scores
        my($singenes, $joingenes);
        if($genegroup) {
          ($singenes, $joingenes)= findjoins_genegroup( \@sgenes);
        } else {
          ($singenes, $joingenes)= findjoins( \@sgenes);
        }
        
        if(@$joingenes > 0 or @mustdemotejoin > 0) { @sgenes=( @$singenes, @$joingenes, @mustdemotejoin); }
      }
 
      
      my %skipcheck;
      foreach my $gref (@sgenes) {
        my($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid)= @$gref;
        my $rexons= $exons->{$gid};
        # next unless($rexons); # are we missing exon-gene id links?
  
        # problem here with skip 1st then unoverlapped later? No, my mistake
        # which takes priority: keepgene or skipgene? keepgene implies already output
  
        if($keepgene{$gid}) { # should not be here now
          # push(@didexons, @$rexons); #? is this needed? 2-bin gene already
          warn "#WARN: already kept $gid\n" if $debug;
          
        } elsif($skipgene{$gid}) { # should not be here now?
          ## $skipcheck{$gid}++ unless($skipgene{$gid}>99);
          
        } else {
          my ($nover, $xbad, $nex, $baseo, $baselen)= (0) x 10;
          my @overid1=(); my @orevid1=();
          # we are losing some gene models, possibly UTR-only overlaps;
          # count CDS overlaps only?
          
  # FIXME1705: CDS-overlap1 should also use oSTRANDED, unless 1-cds-exon
          my @cexons; # FIXME1705
          if($typeover) { @cexons= grep{ $_->[jTYPE] =~ m/$typeover/ } @$rexons; }
          else { @cexons=  @$rexons; }
          my $ovflags= (@cexons > 1) ? oSTRANDED : 0; # FIXME1705 : stranded CDS test, for >1 exon
          foreach my $ex (@cexons) { # @$rexons
            #o next if($typeover and $ex->[jTYPE] !~ /$typeover/);
            $nex++;
            $xbad++  if ($badattr && $ex->[jATTR] =~ m/$badattr/);
            if (overlap1($ex, \@didexons, $typeover, $ovflags)) {
              $nover++; 
            } elsif($overlap1_ovrev) { 
              $ex->[jATTR]=~s/$/;revcdsov=1/;    # flag it or @overids? but not for @overid1
              push @orevid1, @orevids; # save for skip report
            } 
            
            push @overid1, @overids; # save for skip report
            $baseo  += $overbases; 
            $baselen+= $overlen; #** overlen here WAS wrong; instead  $xw= $ex->[4] - $ex->[3]
            }
          $baselen ||= 1;
          $nover=0 if($nover>0 and $baseo > 0 and  $trivialover > 0 and $baseo/$baselen <= $trivialover);
          
          # my $nrevover=0; # see below
          # if($SKIP_REVCDS and @orevid1) { 
          #   my $nrev= @orevid1; my $ncds= @cexons;
          #   $nrevover= $nrev if($nrev >= 0.50 * $ncds);
          # }
          
          # Insert here? typeover2=exon, trivial2=33 (over1=CDS)
          # * this is dropping genes w/ good evidence for ones w/o: ~3,000 with prots
          # * check strand, permit at least some rev strand genes overlapping UTR?

          # ** Problem here w/ mostly UTR models overlap many CDS, should not keep

        if($DO_OVEREXON2 and $nover == 0 and $xbad == 0 and $typeover !~ /exon/) {
          ## my $trivial2=$UTR_OVERLAP_MAX; #  ANOTHER OPTION: 0.50; # 0.33;  # make this higher; some real genes lost due to long utr overlaps
          ## OR chomp out all those AUG25 long utrs that are unreal, causing problems.
          ## can we just leave long UTR exons out of @didexons? probably not,
          ## give precedence to good CDS overlapping long utr?
          ## what this really should filter is small-CDS model overlaps not checked above
          
          # **? add other criteria to keep for UTR-only overlap:
          #  -- if this has long CDS and/or good homology score should keep
          
          my $typeover2="exon"; 
          my $cdsw = $gref->[jCDSLEN] || 1; # "hidden" gref attributes
          my $exonw= $gref->[jTRLEN] || 1; 
          my $utrw= $exonw - $cdsw; $utrw=1 if($utrw < 1);
          my $skiptest= ($UTR_OVERLAP_MAX >= 1 or $utrw < $UTR_CDS_TEST * $exonw )?1:0;  # ANOTHER OPTION
          
          my $CDSKEEP_FROM_UTROVER = 400; # 180; # $dropscore[xxxCDS]
          $skiptest=1 if($cdsw >= $CDSKEEP_FROM_UTROVER  and $utrw/$exonw < $UTR_TOO_BIG and $gref->[jSCORESUM] > 300);
          # other criteria to keep from utr-overlap: CDS introns, homology, ..
          
          unless($skiptest) {
          my ($baseo, $baselen)= (0) x 10;
          ## my $saveto= $typeover; $typeover= $typeover2;
          foreach my $ex (@$rexons) {
            next if($ex->[jTYPE] !~ /$typeover2/);
            #o $nover++ if (overlap1($ex, \@didexons, $typeover2, oSTRANDED));
            if (overlap1($ex, \@didexons, $typeover2, oSTRANDED)) {
              $nover++;
              $ex->[jATTR]=~s/$/;ovutr=1/;    # flag it or @overids? but not for @overid1
            }
            push @overid1, @overids; # save for skip report
            $baseo  += $overbases; 
            $baselen+= $overlen;  
            }
          $baselen -= $cdsw; # should == $exonw - $cdsw; test only UTR, cdsover == 0
          $utrw = $baselen if($baselen > $utrw); #??
          $nover=0 if($nover>0 and ($baseo < 30 or  $baseo/$utrw <= $UTR_OVERLAP_MAX));
          $gref->[jATTR]=~s/$/;utrover=$nover/ if($nover);   # flag mRNA also;  
          ##$typeover= $saveto; 
          }
          
        }        
          
          if($doalttr and $nover>0 and $nex > $nover) { # keep if 1+ (CDS)exon is unique
            $nover=0;
            $gref->[jATTR]=~s/$/;alttr=1/;   # flag; add ID of main tr?
          }
 
          #upd1802: record overid1, orevid1 regardless? of nover/xbad
          if(@overid1 or @orevid1) {
            my($ovkey,@ovid)= ($nover>0 and @overid1) ? ("overids",@overid1) 
                : (@orevid1)? ("orevids",@orevid1) : ();
            if(@ovid) {
            my %ovid= map{$_,1} @ovid; 
            my $ovid= join",", sort keys %ovid;
            $gref->[jATTR]=~s/$/;$ovkey=$ovid/;  
            
            if($SKIP_REVCDS and @orevid1) {  
              # $nover= @orevid1; 
              my $nrev= @orevid1; my $ncds= $nex; # @cexons;
              $nover= $nrev if($nrev >= 0.50 * $ncds);              
              } 
            ## SKIP_REVCDS is this one partial overlap rev cds? or most rev cds?  does mrna have -sense flag?
            ## is it ~same prot as rev-over, but antisense map orient? 
            ## prob should score/filter this problem case elsewhere
            }
          }
                    
          if($nover > 0 || $xbad > 0) { # maybe keep if $nover < $nx * 0.75
            $skipgene{$gid} = ($xbad > 0) ? kSKIP_BAD : kSKIP_OVERLAP; #was ++; 
            $nskip++;  
            # if($nover > 0 and @overid1) {
            #   my %overid1= map{$_,1} @overid1; 
            #   my $overid1= join",", sort keys %overid1;
            #   $gref->[jATTR]=~s/$/;overids=$overid1/; # if($overid1 and $showskips);
            #   }
            #? do we have any case of gid kept before, then skipped?, but putgene would have kept
            # if so, need to remove from @didexons, ...
            }  
          else { 
            $keepgene{$gid}++; $nkeep++;
            push(@didexons, @$rexons);
            putgene( $gref, $rexons); 
            
            @didexons= sort { $a->[jBEGIN] <=> $b->[jBEGIN] } @didexons; # new 2011.10 for overlap1 test
            }
          }
        
        }
  
      
      # rev sort by end loc  #$ex->[3], $ex->[4],  = tb,te,
      @lastbinexons=();
      foreach my $ex (sort { $b->[jEND] <=> $a->[jEND] } @didexons) { 
        push(@lastbinexons, $ex); last if(@lastbinexons > 50); 
        }
        
      }
      
  }
}


#.......................................

    #? can we check gene joins/splits here?
    #  for given 1st gene at location, collect all mostly overlapped genes,
    #   if any has >=150% exons or <= 50% exons, check for join/split
    #   "best" here may be sum of evidence scores of mRNAs per source-type 
    #    over longest span window?
    # FIXME for gene-joins ; often have ==homology, but >> CDS
    #   : find all mRNA overlaps to best, check for ~= ho, << CDS

    # my($singenes, $joingenes)= findjoins( \@sgenes);

sub over1gene {
  my( $ref,$tb,$te,$to,$rloc)= @_;
  my ($lref,$lty,$lb,$le,$lo,$lid)= @{$rloc}[0,2,3,4,6,9];
  my $over= ($ref eq $lref && $tb <= $le && $te >= $lb && _eqstrand($lo,$to)) ? 1 : 0;
  #** ADD EXON CHECK?  defer till we need it
  return $over;
}

sub overgenes {
  my( $ref,$tb,$te,$to, $sgenes, $ig, $ng)= @_;
  my @overs=();
  for(my $i=$ig; $i<$ng; $i++) {
    my $g2= $sgenes->[$i];
    push(@overs, $g2) if(over1gene($ref,$tb,$te,$to, $g2));
  }
  return @overs;
}

sub overexons {
  my( $gref, $oref)= @_;
  my $rexons =  $exons->{$gref->[jGID]}; 
  my $lexons =  $exons->{$oref->[jGID]}; 
  if(@$lexons > @$rexons) { my $sx= $lexons; $lexons= $rexons; $rexons= $sx; }
  # my $nex    = @$lexons; # oops, adjust for exon typeover skip
  ##  my $rloc= [$ref,$src,$typ,$tb,$te,$score,$to,$tph,$tattr,$gid]; # exon/mrna record
  my $savepctover = $pctover; $pctover=0; # skip this
  my $nover  = 0;
  my $nex    = 0; # oops, adjust for exon typeover skip
  foreach my $ex (@$lexons) {
    next if($typeover and $ex->[jTYPE] !~ /$typeover/);
    $nex++; $nover++ if (overlap1($ex, $rexons, $typeover, 0));
    }
  $pctover= $savepctover ;
  return 0 unless($nex>0);
  return $nover / $nex; # want to know ratio of nover / min(rexons, lexons) ?
}


#    my @gscores= pickscores_genegroup $grjoinscores,$tscore,$tattr;   
sub pickscores_genegroup {  
  my($flds, $tscores, $tattr)=@_;
  my @fl= grep { $_ ne "CDS" } map{ s,\W.*,,; $_} split",",$flds;

  my %fldmax=();
  foreach my $sf (@fl) { if($flds =~ m/\b$sf[:=]([\w.]+)\b/) { $fldmax{$sf}= $1; } }
## FIXME: use $genegroup[$i] or genegroup='homolog,ovpro:100' < 100 gives score max
##       if($genegroup =~ m/\b$sf[:=]([\w.]+)\b/) { $genegroup[$i]= $1; }


  #  $tattr contains best_$sf1= scores; $tscores = fld1,fld2,.. for sfields
  my @sc= split",",$tscores;
  my @sf= manyscorefields(0);
  my %sfi; for my $i (0..$#sf) { $sfi{$sf[$i]}= $i; }

  my ($clen,$xlen)= $tattr =~ m/cxlen=(\d+).(\d+)/; # ok, cxlen= set in collect_gff
  $clen||= 0; $xlen||=0; # some mistakes, missing tattr ?
  # push(@ret, $clen, $xlen); # add $xlen?
  ## new intron valid score here?  inexon == KeyNINTRON ??
  my ($validin,$nexon)=(0,0); 
  if($tattr =~ m/$KeyNINTRON=(\d+).(\d+)/) { ($validin,$nexon)= ($1,$2); } 
  elsif($tattr =~ m/inexon=(\d+).(\d+)/) { ($validin,$nexon)= ($1,$2); }
  $nexon ||= 1; $validin ||= 0;

  my ($joinck,$maybejoin)=(0,0); 
  if($tattr =~ m/$KeyJOINCHECK=([\d\.+-]+)/) { $joinck= $1; $maybejoin=($joinck < 0)?(-$joinck):0; } 
  
  my $havescore= 0;
  my @ret=();
  foreach my $fl (@fl) { 
    my $j= $sfi{$fl}; 
    my ($v)= $tattr =~ m/best_$fl=([\d\/.e+-]+)/;
    unless($v){ ($v)= $tattr =~ m/\b$fl=([\d\/.e+-]+)/; }
    if(not $v and defined $j) { $v= $sc[$j]; }
    if($v)  { 
      $havescore++; 
      # fixme: append group maxscore : $v/$vmax
      if($v !~ m,/, and (my $vmax=$fldmax{$fl})) { $v= "$v/$vmax\%"; }
    } else { $v=0; }
    push(@ret, $v); 
    }
  @ret=() unless($havescore); ## { return ($clen,$xlen,$validin,$nexon,$maybejoin); } # or empty?
  return ($clen, $xlen,$validin,$nexon,$maybejoin,@ret); #? check for all empty scores?
}




sub findjoins_genegroup { # _v2 
  my($sgenes)= @_;
  my(@singles, @joins);
  
    # should be options
#  my $HOSLOP  = 10; #? bitscore diff; use ratio? 0.95
#  my $HORAT   = 0.95;

  my $PEXONOVER= 0.50;
  my $PGENEOVER= 0.50;
  my $PDOUBLESIZE= 1.75; # test cds > PDOUBLESIZE * holen; or use cds > holen + CDSDIFF ?
  my $CDSDIFF = 300; #? what; ratio?  
  my $grjoinscores= $genegroup || $JOINSCORES; ## = "CDS,ho3,pro,xde"; # FIXME: option ; CDS must be 1st, need ordered
  
  # sgenes is ordered by hi score .. low
  my $ng= scalar(@$sgenes);
  
  foreach my $ig (0..$ng-1) {
    my $gref= $sgenes->[$ig];
    my $isjoin="";
    
    my($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid)= @$gref;

    my @og = overgenes( $ref, $tb,$te,$to, $sgenes, $ig+1, $ng);
    
    # DEFER below# unless(@og) { push(@singles, $gref); next; } # none following are over this gene
    
    my ($cds,$xlen,$valin,$nexon,$maybejoin,@gscores)= pickscores_genegroup( $grjoinscores,$tscore,$tattr);   
    
    # FIXME: 20apr: check INEXON scores to counteract these join scores
    # ... all exons w/ valid introns == true model more likely, still can be mistake tho for bad introns
    # my $valintronjoin= ($valin >= $nexon) ? 1 : 0;  # or valin >= $nexon-1 ?
    # too weak valin; 1 miss of 20 but is middle long false intron:
    # my $valintronjoin= ($valin > 0 and $valin / $nexon >= 0.83) ? 1 : 0;  #  3/4 = .75; 4/5 = .80; 5/6 = .83 6/7 =.85n
    my $valintronjoin= ($valin > 0 and $valin / $nexon >= 0.98) ? 1 : 0;  #  3/4 = .75; 4/5 = .80; 5/6 = .83 6/7 =.85n

    if($maybejoin!=0) {
      #? $gref->[jATTR]=~s/$/;join=${maybejoin}m/;  # dont need, maybejoin is joinck= attr
      push(@joins, $gref); next;
    }
    
    ##unless(@og) { push(@singles, $gref); next; } # none following are over this gene
    unless(@gscores and not $valintronjoin and $cds > 180 + $CDSDIFF) { push(@singles, $gref); next; } # no scores to test
    # ^with many genegroup scores, should these be sorted high first?
    # .. dont test w/ weak homolog scores.
    # .. test only long prots: $cds > $CDSDIF
    
    my $nog= scalar(@og);
    # test not join if all @og match 1 shortest part
    if($nog > 1) { 
      my $minw= 9999999; my $i0=0;
      foreach my $i (0..$#og) {
         my($lb,$le)= ($og[$i]->[jBEGIN], $og[$i]->[jEND]); 
         my $lw=$le-$lb;
         if($lw < $minw) { $i0= $i; $minw= $lw; }
      }
      ##my($lb,$le,$lty)= ($og[$i0]->[jBEGIN], $og[$i0]->[jEND], $og[$i0]->[jTYPE]);
      my @og2 = @og;
      my($og1)= splice(@og2,$i0,1);
      
      my $nog1= overlap2( $og[$i0], \@og2, 0, $PGENEOVER);
      $nog=0 if( $nog1 == scalar(@og2)); # not join if all @og match 1 shortest part
      }
    unless($nog > 1) { push(@singles, $gref); next; }
    
    my @orderscore=(0..$#gscores);
if(0) { #? want this or not? want to use homolog= if real before pro, rna base over scores    
    @orderscore= sort{ 
      (my $sa=$gscores[$a]) =~ s/\D.*//; 
      (my $sb=$gscores[$b]) =~ s/\D.*//; 
      $sb <=> $sa; } @orderscore;
}

    my $didtest= 0;
    foreach my $iscore (@orderscore) {
      my $ho2 = $gscores[$iscore] or next;
      my $holen= 0; my $hoIsPct=0;
      ## handle other num chars: ([\d.e+-]+)
      if($ho2 =~ m,([\d.e+-]+)/([\d.e+-]+)([%]?),) { $ho2=$1; $holen=$2; $hoIsPct=$3; } 
      ##if($ho2 =~ m,(\d+)/(\d+),) { $ho2=$1; $holen=$2; } 
      next unless($ho2>0 and $holen>0); # cannot score unless holen == max score ?      
      # * FIXME: 2011sep: ovpro,ovrna == 100% scores, no xxx/yyy, flag to use?
      # * FIXME2: using ovpro/100 patch fails for holen/cds below **
      
      # maybe in all cases take ratio: ho2/holen ? join == this ratio << other ratio
      # * better would be to have orthlog's selfscore here,
      #  use ho/orthomax as basis to say if is join.
      my $horat= $ho2/$holen;
      $didtest++ if($horat > 0.66);
      # if($hoIsPct) { $holen= $cds; } # FIXME2 fix, do after get horat
      
      my $join=0; 
      my $valincase= 0;
      my $likelyjoin= ($holen > 0 and not $hoIsPct and $cds - $holen > $CDSDIFF #< problem picking this
                       and $cds >= $PDOUBLESIZE * $holen)?1:0;
      # need PDOUBLESIZE here for long prots that mostly match long homologs      
      if($valintronjoin) { $valincase=1; $likelyjoin= 0; }
      
      foreach my $og (@og) {
        my ($ocds,$oxlen,$ovalin,$onexon,$omaybejoin,@ogscores)=
          pickscores_genegroup( $grjoinscores,$og->[jSCORE],$og->[jATTR]);
        my $oho= $ogscores[$iscore] or next;
        my $oholen= 0;  my $ohoIsPct=0;
        if($oho =~ m,([\d.e+-]+)/([\d.e+-]+)([%]?),) { $oho=$1; $oholen=$2; $ohoIsPct=$3; }
        next unless($omaybejoin==0 and $oho>0 and $oholen > 0 and ($cds - $ocds > $CDSDIFF)); #?
       
        # nexon, $onexon > 0 always
        my $notjoin= ($valin > $ovalin and $valin/$nexon > $ovalin/$onexon ) ? 1 : 0;
        $valincase += $notjoin;
        #? next if($notjoin); #? is this ok? problem introns; flag and check results
        
        my $ohorat= $oho/$oholen;
        next unless($ohorat > 0.25); # not good test otherwise
        $didtest++ if($ohorat > 0.66);
        
        # my $likelysplit= ($ocds - $cds > $CDSDIFF)?1:0;          
        # next if ($likelysplit);
        
        # if( ($oho / $ho2 >= $HORAT) and ($cds - $ocds > $CDSDIFF)) 
        if( ($ohorat > 1.75 * $horat) and ($cds - $ocds > $CDSDIFF) )
        {
          my $pover= overexons($gref, $og); # pover is ratio nover/min(nexons)
          $join++ if($pover >= $PEXONOVER); #? or more
        }         
        
        # last if($join>0);  #? test all og
        }
      
      ## debug all types of join info
      $isjoin .= "l1" if $likelyjoin; #? test this
      if ($join>0) { $isjoin .= "p".$join; }  # >0 or  >1 ?
      ##if($xjoin>0) { $isjoin .= "x".$xjoin; }
      if($isjoin and $valincase) { $isjoin .= "vi".$valincase; }
      
      last if ($didtest>0);      
    }

    if($isjoin) { 
      $gref->[jATTR]=~s/$/;join=$isjoin/;   # tattr
      push(@joins, $gref); 
    } else { 
      push(@singles, $gref); 
    }
  }
  
  return (\@singles, \@joins);
}


sub pickscores {  
  my($flds, $scores)=@_;
  my @ret=();
  my @fl= split",",$flds;
  my @sc= split",",$scores;
  my @sf= manyscorefields(0);
  my %sfi; for my $i (0..$#sf) { $sfi{$sf[$i]}= $i; }
  foreach my $fl (@fl) { my $j= $sfi{$fl}; my $v=(defined $j)?$sc[$j]:0; 
    push(@ret, $v); }
  return @ret;
}

sub findjoins {  # FIXME should drop this
  my($sgenes)= @_;
  my(@singles, @joins);
  
    # should be options
  my $HOSLOP  = 10; #? bitscore diff; use ratio? 0.95
  my $HORAT   = 0.95;
  my $XDERAT  = 1.2; # 1.5?; is this good to indicate loss of xde from join?
  my $XDERAT2 = 1.5;  
      ## one test case has 4.497/3.67 = 1.22 for single/join XDE ratio
      ## xde is problematic indicator; back off, use only when no homology
  my $PEXONOVER= 0.50;
  my $PGENEOVER= 0.50;
  my $PDOUBLESIZE= 1.75; # test cds > PDOUBLESIZE * holen; or use cds > holen + CDSDIFF ?
  my $CDSDIFF = 300; #? what; ratio?  
  
  # sgenes is ordered by hi score .. low
  my $ng= scalar(@$sgenes);
  
  foreach my $ig (0..$ng-1) {
    my $gref= $sgenes->[$ig];
    my $isjoin="";
    
    my($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid)= @$gref;
    
    # add pro if ho2 == 0 ? see below proscores fix; check proscores for multiple genes?
    # 2010.jul: add new join score from _self blastp analyses
    #   .. easy to pick gene1 == gene2a + gene2b putative join (or 2a,b split)
    
    # $JOINSCORES = "ho2,pro,CDS" # HACK: fixme
    # AUG25 has errs of long UTRs that are not UTRs; no CDS diff but long UTR 
    #  mistakenly joins  other genes : add UTR test?  utrs = exonbases - cds

# .. fix findjoins to require alternates to cover 2+ regions of putative join?
#      .. i.e. don't drop longer gene just from shorter/same-ho on one half.
# ** add another score for non-homol gene joins; what? self-prot-homol= paralogy scores?
#    para score would work like ho score: best bitscore per gene model; 
#      join score should not exceed singles scores : BUT for cases of bad ortho models(?many?)
#      .. now using ho3=ho2 + paralog bits if ho2==0 (tried para>ho2, problems?)
#    maybe favor pro if pro >> ho2 : * MAYBE NOT, got mistake w/ same pro for tandems
#    ** add ho aa length to aid join test: if this.cds >> ho.cdslen, likely join
#    .. can we use converse this.cds << ho.cdslen as split info?

    my ($cds,$ho2,$pro,$xde)= pickscores $JOINSCORES,$tscore;   
    # cds score field is same as   $cds= $gref->[10]; # "hidden" gref attributes
    # change pro, ho2 scores to add target cds length as "align/length"
    my ($holen,$prolen)= (0,0);
    if($ho2 =~ m,([\d.e+-]+)/([\d.e+-]+),) { $ho2=$1; $holen=$2; }
    if($pro =~ m,([\d.e+-]+)/([\d.e+-]+),) { $pro=$1; $prolen=$2; }
    
    if($ho2 > 0 or $pro > 0 or $xde > 0) 
    {
      # my $usepro= ($ho2 < 1)?1:0;
      my $usepro= ($ho2 < 1 or $pro > 2*$ho2)?1:0;  #?? yes or no
      if ($usepro) { $ho2= $pro; $holen=$prolen; }
      elsif($holen == 0) { $holen=$prolen; } # can use this instead?
      
      my @og = overgenes( $ref, $tb,$te,$to, $sgenes, $ig+1, $ng);
      ## here require @og > 1 AND @og cover >1 part of gref long gene
      my $nov= scalar(@og);
      
      if($nov > 1) { 
        # are any overlaps outside of 1st ? test pctover < 0.50
        # no not 1st og, need to pick shortest
        # my($lb,$le,$lty)= ($og[0]->[3], $og[0]->[4], $og[0]->[jTYPE]);
        # my @og2 = @og[1..$nov-1];
        my $minw= 9999999; my $i0=0;
        foreach my $i (0..$#og) {
           my($lb,$le)= ($og[$i]->[jBEGIN], $og[$i]->[jEND]); my $lw=$le-$lb;
           if($lw < $minw) { $i0= $i; $minw= $lw; }
        }
        ## my($lb,$le,$lty)= ($og[$i0]->[3], $og[$i0]->[4], $og[$i0]->[jTYPE]);
        my @og2 = @og;
        my($og1)= splice(@og2,$i0,1);
        
        my $nov1= overlap2( $og[$i0], \@og2, 0, $PGENEOVER);
        $nov=0 if( $nov1 == scalar(@og2)); # not join if all @og match 1 shortest part
        }
        
      if($nov > 1) {        
        my $join=0; my $xjoin=0;
        # my $likelyjoin= ($holen > 0 and $cds >= $PDOUBLESIZE * $holen)?1:0;
        my $likelyjoin= ($holen > 0 and $cds - $holen > $CDSDIFF
                         and $cds >= $PDOUBLESIZE * $holen)?1:0;
        # need PDOUBLESIZE here for long prots that mostly match long homologs
        
        foreach my $og (@og) {
          my ($ocds,$oho,$opro,$oxde)= pickscores( $JOINSCORES, $og->[jSCORE]);
          my ($oholen,$oprolen)= (0,0);
          if($oho =~ m,([\d.e+-]+)/([\d.e+-]+),) { $oho=$1; $oholen=$2; }
          if($opro =~ m,([\d.e+-]+)/([\d.e+-]+),) { $opro=$1; $oprolen=$2; }
          if ($usepro) { $oho= $opro; $oholen=$oprolen; }
          elsif($oholen == 0) { $oholen=$oprolen; } # can use this instead?
          my $likelysplit= ($oholen > 0 and $oholen - $cds > $CDSDIFF)?1:0;          
          next if ($likelysplit);
          
          #if(($oho + $HOSLOP > $ho2) and ($cds - $ocds > $CDSDIFF)) 
          if($ho2 > 0 and ($oho / $ho2 >= $HORAT) and ($cds - $ocds > $CDSDIFF)) {
            my $pover= overexons($gref, $og); # pover is ratio nover/min(nexons)
            $join++ if ($pover >= $PEXONOVER); #? or more
            }
            
          if($xde > 0) {
            my $dode= ($ho2 < 1 or ($oxde / $xde > $XDERAT2));
            if($dode and ($oxde / $xde > $XDERAT) and ($cds - $ocds > $CDSDIFF)) {
              my $pover= overexons($gref, $og); 
              $xjoin++ if ($pover >= $PEXONOVER);  
              }
            }
          
          # last if($join>0);  
          }
      
      ## debug all types of join info
      $isjoin .= "l1" if $likelyjoin; #? test this
      if ($join>0) { $isjoin .= "p".$join; }  # >0 or  >1 ?
      if($xjoin>0) { $isjoin .= "x".$xjoin; }
      }
    }
    
    if($isjoin) { 
      $gref->[jATTR]=~s/$/;join=$isjoin/;   # tattr
      push(@joins, $gref); 
    } else { push(@singles, $gref); }
  }
  return (\@singles, \@joins);
}


sub _sort_refscoreloc # $overlaps == mrnatypes
{
# ($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid)
# fixed: ?? with ref 1st this isn't working; big score not at front
  return ($a->[jCHR]   cmp $b->[jCHR]) # ref
      || ($b->[jSCORE] <=> $a->[jSCORE]) # score : not numeric err from what?
      || ($a->[jBEGIN] <=> $b->[jBEGIN]) # begin
      || ($b->[jEND]   <=> $a->[jEND]) # end
      || ($a->[jSRC]   cmp $b->[jSRC]); # src
}


sub manyscorefields {
  my($keepwt)=@_;
  return() unless($manyscore);
  (my $sf= $scoretype) =~ s/^many.//; 
  if($keepwt) { return (split",",$sf); }
  else { return map{ s/:.+$//; $_ } split",",$sf; } # was s/:[\-\d]+//
}

# see below gff
sub add_required_scorefield {
  my($sname, $sweight)=@_;  #, $sindex
  # update: $scoretype, @dropscore, @droptype, @genegroup ?
  (my $sf= $scoretype) =~ s/^many.//; 
  my @sf= (split",",$sf);
  for my $i (0..$#sf) { return $i if($sname eq $sf[$i]); }
  
  # need sindex > $#sf to add ; grep for $sname
  # ?? need sindex == @sf to add
  my $sindex = @sf; # always to add? unless($sindex >= @sf);
  $sweight ||= 1;

  $scoretype .= ",$sname:$sweight";
  $scoreweight[$sindex] = $sweight;
  $dropscore[$sindex] = 0; # default?
  $droptype[$sindex] = 0;
  return $sindex;
}





sub manyscore_sum {
  my($score)= @_;
  my @as= split",",$score;
  my $n = scalar(@as); # fixme: does @as always match @scoreweight? err if not? or have extra flds?
  my $nw= scalar(@scoreweight); 
  die "ERROR: manyscore_sum scorevec n=$n != scoreweight nw=$nw" unless($n == $nw);
  my($as1)= (0);
  for( my $i=0; $i<$n; $i++) {
    my($as2)=($as[$i]);  
    $as2 =~ s,/.*,,;   # FIXME: pro score now like aln/size
    # 10/09/18: can we add -scoreweight, for terepeats?
    $as1 += $scoreweight[$i] * $as2;
  }
  return int($as1);
}



sub _sort_manyscoreloc # manyscore: val1,val2,val3,...
{
  my($as1,$bs1)= (0,0); # score to compute
  
    ## $a == \($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid)
#   my @as= split",",$a->[jSCORE];
#   my @bs= split",",$b->[jSCORE];
#   my $n= scalar(@as);
  ## add more scores to sort on: genegroup best_xxx, but not? in tscore field [5]
  ## .. convert genegroup best_ scores to two: 
  #   1. CDS best = match all of homolog prot; 
  #   2. UTR best = match all of exon group (EST/Rnaseq/...) BUT NOT too much more
  
  ## fixme for UTR score: (exon-cds)/exon ; poor when > 0.5, good when <= 0.5
  ## maybe this should be sum of all scores, not best?

  #2011.9: mod ordersrc to be, for some src, a weight * sumscore?
  # my($asrc,$bsrc)= map{ $ordersrc{$_} || $_ } ($a->[jSRC], $b->[jSRC]);
  my($asrc,$bsrc)= ($a->[jSRC], $b->[jSRC]);
  #below#my($apsrc,$bpsrc)= (1,1);# if($weightsrc) { ($apsrc,$bpsrc)=map{ $weightsrc{$_} || 1 } ($asrc,$bsrc); }

if($COMBINE_SCORE == SCORE_SUM) { # default now; precompute this == manyscore_sum ?
  $as1 = $a->[jSCORESUM];  # added precompute; FIXME: opt to erase old precompute
  $bs1 = $b->[jSCORESUM];
  unless($as1 =~ /\d/ and $bs1 =~ /\d/) {  # drop if precompute ok; should not be zero but might be
    ## should use instead manyscore_sum() as below
    ($as1, $bs1)= (0,0);
    my @as= split",",$a->[jSCORE];
    my @bs= split",",$b->[jSCORE];
    my $n= scalar(@as);
    for( my $i=0; $i<$n; $i++) {
      my($as2,$bs2)=($as[$i],$bs[$i]);  
      $as2 =~ s,/.*,,; $bs2 =~ s,/.*,,;  # FIXME: pro score now like aln/size
      # 10/09/18: can we add -scoreweight, for terepeats? yes
      $as1 += $scoreweight[$i] * $as2;
      $bs1 += $scoreweight[$i] * $bs2;
    }
  #bug# if($weightsrc) { $as1= $apsrc * $as1; $bs1= $bpsrc * $bs1; } # FIXME: move out of unless()
  }
  
} else { # SCORETYPE == SCORE_BEST
  my $ibest=0;
  my @as= split",",$a->[jSCORE];
  my @bs= split",",$b->[jSCORE];
  my $n= scalar(@as);
  ($as1,$bs1)=($as[0],$bs[0]);
  $as1 =~ s,/.*,,; $bs1 =~ s,/.*,,;  # FIXME: pro score now like aln/size
  
  for( my $i=1; $i<$n; $i++) {
    my($as2,$bs2)=($as[$i],$bs[$i]);
    $as2 =~ s,/.*,,; $bs2 =~ s,/.*,,;  # FIXME: pro score now like aln/size
    if($as2 == $bs2) {  #next;
    } elsif( ($bs1+$as1) < 1 || ($bs1 == $as1 && $as2 != $bs2) ) {
      $as1= $as2; $bs1= $bs2; $ibest= $i;
    } elsif( ($bs2+$as2) > 0  ) {
      # 10/09/18: can we add -scoreweight, for terepeats?
      my $sw2= ($scoreweight[$i] * abs($bs2 - $as2))/($bs2+$as2) ;
      my $sw1= ($scoreweight[$ibest] * abs($bs1 - $as1))/($bs1+$as1);
      if($sw2 > $sw1) { $as1= $as2; $bs1= $bs2; $ibest=$i; }
    } 
  }
}

  if($weightsrc) {
    my($apsrc,$bpsrc)=map{ $weightsrc{$_} || 1 } ($asrc,$bsrc);  
    $as1= $apsrc * $as1; $bs1= $bpsrc * $bs1; 
  } 
    
  ## ordersource : use $genesources count from collect_gff for default order?
  my($aisrc,$bisrc)= (0,0);
  if($ordersource) { 
    my $nsrc= $ngenesources; # scalar(keys %$genesources); # so that unspec wont all be last if srcB:99
    ($aisrc,$bisrc)= map{ $ordersrc{$_} || $nsrc } ($asrc,$bsrc); 
  }
  
  return ($a->[jCHR] cmp $b->[jCHR]) # ref
      || ($bs1 <=> $as1) # score
      || ($aisrc <=> $bisrc) # move above begin/end ? for equal scores, this should take precedence 
      || ($a->[jBEGIN] <=> $b->[jBEGIN]) # begin
      || ($a->[jEND]   <=> $b->[jEND]) # end, shortest 1st here
      || ($asrc cmp $bsrc);  # src : should allow option how to sort otherwise ident models
}


# DROP dualscore
# sub _sort_refscore2loc # dualscore: special case 2 scores, need to sort on rel diff


sub putsummary {
  # output where? STDERR? end of STDOUT ?
  my $info=""; 
  $info .= "scoretype=$scoretype, ";
  print $OUT "\n# SUMMARY of scores/source, $info\n";
  ## print $OUT "# ",join("\t",qw(Source N  Ave Weight Sum)),"\n";
  
  my @src= sort keys %summary;
  my @scs= grep !/^(n|ndrop|sum|wt)$/, sort keys %{ $summary{ $src[0] } };
  my @scols=@scs; map { s/^s_// } @scols;
  
  print $OUT "# ",join("\t",qw(Source N  Ave), @scols),"\n";
  
  foreach my $src (@src) {
    my $sm= $summary{$src}{'sum'} || 0;
    my $n = $summary{$src}{'n'}   || 1;
    my $wt= $summary{$src}{'wt'}  || $n; # weight bad for neg weights
    # my $av= sprintf "%.0f", $sm / $wt; # or $n; ## $wt; 
    ###? drop wt, sm; add per scoretype ave?
    ## printf $OUT "# %-10s\t",$src; print $OUT join("\t",$n,$av,$wt,$sm),"\n";

    printf $OUT "# %-10s\t%6d\t%5.0f",$src, $n, $sm / $wt; 
    foreach my $sc (@scs) {
      my $sm= $summary{$src}{$sc} || 0;
      printf $OUT "\t%5.1f", $sm / $n;
      }
    print $OUT "\n";

  }
}


sub putgene {
  my($gref,$exons,$flag)= @_;
  #fixed above: warn now;  skipgenes mustdrop not found
  unless($gref and ref $gref) { warn "# Bad call to putgene($gref,$exons,$flag); \n"; return; }
  my @sf = manyscorefields(0);
  foreach my $ex ($gref, @$exons) {
    # new bug: not defined $ex : got here from put skipgene, with null skipgene from mustdrop list
    next unless($ex and ref $ex);
    my($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid)= @$ex;
    my $skipit= ($flag && $flag =~ /skip=/)?1:0;
    my $cskip= ($skipit)?"#x.":"";
    $tattr .=";$flag" if($flag);
    
    my $pscore=$tscore; 
    unless($VECSCORE) { 
      if($tscore =~ m/,/) { # allow for not manyscore, tscore == 'score' field only
        $pscore= $ex->[jSCORESUM] || 1; # only in mRNA/gref
        if($typ =~ /$mrnatypes/) { 
          unless($tattr =~ s/;scorevec=[^;\s]*/;scorevec=$tscore/) { $tattr =~ s/$/;scorevec=$tscore/; }
        }
      } 
    }
    print $OUT join("\t", $cskip.$ref, $src, $typ, $tb, $te, $pscore, $to, $tph, $tattr),"\n";

    if($summarize && $ex eq $gref) { # tscore summary / $src  
      my $skiptype = ($skipit) ? "2_Dropped" : "1_Kept";
      # summarize by $src (all), by kept/skip, want src-kept/skip also?
      my($sumv, $sumwt)= (0,0);
      if($manyscore) {
        my @sv = split",", $tscore;
        for my $i (0..$#sf) { 
          my $sf= $sf[$i];
          # 10/09/18: can we add -scoreweight, for terepeats?
          my $sw= $scoreweight[$i];
          my $sv= $sv[$i];  

          map{s,/.*,,}($sv); # "$pv/$px" # Argument "72/232" isn't numeric in numeric
          $sv ||= 0;
          
          my $svw= $sv * $sw; 
          $sumv  += $svw; 
          $sumwt += $sw if($sv > 0); # not right for ave score
          
          $summary{$src}{"s_".$sf} += $sv unless($skipit);
          $summary{"$src.all"}{"s_".$sf} += $sv;
          $summary{$skiptype}{"s_".$sf} += $sv;             
        } 
        
      } else {
        $sumv= $tscore;
      }
      
      $summary{$src}{'ndrop'}++ if($skipit); 
      foreach my $stype ($src, "$src.all", $skiptype) {
        # leave out of summary for src if skiptype =~ Dropped
        unless($stype eq $src and $skipit) {
        $summary{$stype}{'sum'} += $sumv;
        $summary{$stype}{'n'}++; 
        $summary{$stype}{'wt'} += $sumwt; # for weighted,
        }

      }
    }
  }
}


sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }
# fixme for "." strand compare
sub _eqstrand { return($_[0] eq $_[1] || $_[0] eq '.' || $_[1] eq '.' )?1:0; }

# sub overlaps {
#   my($tb, $te, $typ, $locs) = @_;
#   ##  my $rloc= [$ref,$src,$typ,$tb,$te,$score,$to,$tph,$tattr,$gid]; # exon/mrna record
#   $overbases = 0; $overlen= 0; # globals for base counts
#   return 0 if($typeover and $typ !~ /$typeover/);
#   $overlen= 1+$te-$tb;  
#   foreach my $rloc (@$locs) {
#     my ($lty,$lb,$le,$lid)= @{$rloc}[2,3,4,9];
#     next if($typeover and $lty !~ /$typeover/);
#     my $over= ($tb <= $le && $te >= $lb) ? 1 : 0;
#     # add option to screen out trival 5% UTR overlaps .. pctoverlap=i
#     # also option to use CDS-exons only ignoring UTR overlaps
#     # .. opt for diff pctover for CDS and exon typeover?
#     if($over and $pctover) {
#       my ($bb,$be)= ( _max($tb,$lb), _min($te,$le) );
#       my $maxo= 1+abs($be - $bb);
#       my $leno= _min( abs($le - $lb), abs($te - $tb)) || 1;
#       #BAD# $overlen= $leno; # globals for base counts
#       $over = 0 if $maxo/$leno < $pctover;
#       $overbases = $maxo if($overbases < $maxo); 
#       }
#     return 1 if($over);
#     }
#   return 0;
# }



sub overlap1 {
  my($rexon, $locs, $typeover, $flags) = @_;
  my($typ, $tb, $te, $to)= @{$rexon}[2,3,4,6];
  ##  my $rloc= [$ref,$src,$typ,$tb,$te,$score,$to,$tph,$tattr,$gid]; # exon/mrna record
  @overids=(); @orevids=();
  $overbases = 0; $overlen= 0; $overlap1_ovrev=0; # globals for base counts
  return 0 if($typeover and $typ !~ /$typeover/);
  my $nover=0; #upd1802.globa: my @orevids=();
  $overlen= 1+$te-$tb;  
  foreach my $rloc (@$locs) {
    my ($lty,$lb,$le,$lo,$lid)= @{$rloc}[2,3,4,6,9];
    next if($typeover and $lty !~ /$typeover/);
    my $over= ($tb <= $le && $te >= $lb) ? 1 : 0;
    if($over and $flags & oSTRANDED) {
      unless( _eqstrand($lo,$to) ) { $overlap1_ovrev++;  push @orevids, $lid; $over= 0; }
    }
    #o $over= 0 if($flags & oSTRANDED and ! _eqstrand($lo,$to) ); ## $lo ne $to
    if($over and $pctover) {
      my ($bb,$be)= ( _max($tb,$lb), _min($te,$le) );
      my $maxo= 1+abs($be - $bb);
      my $leno= _min( abs($le - $lb), abs($te - $tb)) || 1;
      $over = 0 if $maxo/$leno < $pctover;
      $overbases = $maxo if($overbases < $maxo); 
      }
    do{ $nover++; push @overids, $lid; } if($over);
    return $nover if($nover>2); # or test all? ##if($over);
    #return $nover if($nover>0 and !$over); # this should work if loc-sorted 
    # problem here: 1st over may be tiny, but want maximal overlap. input locs not now sorted (or by >score)
    # change: locsort $locs; count over, return if nover>2 < or not, dont sort locs, just test all?
    }
  #? if($overlap1_ovrev and not $nover) { @overids= @orevids; } # want these? YES, some strand=-1 are false antisense
  return $nover; ##($nover>0)?1:0;
}


##        my $nov1= overlap2( $og[$i0], \@og2, 0, $PGENEOVER);
sub overlap2 {
  my($rexon, $locs, $flags, $pctover) = @_;
  ##my($tb, $te, $typ, $locs, $pctover, $typeover) = @_;
  my($typ, $tb, $te, $to)= @{$rexon}[2,3,4,6];
  my $typeover= $typ;
  $overbases = 0; $overlen= 0; # globals for base counts
  ## return 0 if($typeover and $typ !~ /$typeover/);
  $overlen= 1+$te-$tb;  
  my $nover=0;
  foreach my $rloc (@$locs) {
    my ($lty,$lb,$le,$lid)= @{$rloc}[2,3,4,9];
    next if($typeover and $lty !~ /$typeover/);
    my $over= ($tb <= $le && $te >= $lb) ? 1 : 0;
    if($over and $pctover) {
      my ($bb,$be)= ( _max($tb,$lb), _min($te,$le) );
      my $maxo= abs($be - $bb);
      my $leno= _min( abs($le - $lb), abs($te - $tb)) || 1;
      $overbases = $maxo; 
      $over = 0 if $maxo/$leno < $pctover;
      }
    $nover++ if($over); # WAS return 1
    }
  return $nover;
}


sub overlocuslist {
  my($ref,$tb,$te) = @_;
  return 0 unless($locuslist{$ref});
  foreach my $ib (int($tb/$BINSIZE) .. int($te/$BINSIZE)) {
    $locuslist{$ref}{$ib} or next;
    foreach my $lref ( @{$locuslist{$ref}{$ib}} ) {
      my ($lr,$lb,$le)= @$lref; # [$ref,$tb,$te]; 
      ##my $over= ($tb <= $le && $te >= $lb) ? 1 : 0;
      return 1 if($tb <= $le && $te >= $lb);
    }
  }
  return 0;
}

# my %locuslist=(); # global
sub readlocuslist {
  open(F, $locuslist);
  my ($nlocus,$span)=(0,0);
  while(<F>){
    next unless(/^\w/);
    my @v=split;  my($ref,$tb,$te)= ($v[0],0,0);
    if($v[1]=~/^\d/ and $v[2]=~/^\d/) { ($tb,$te)= @v[1,2]; } # bed 
    elsif($v[3]=~/^\d/ and $v[4]=~/^\d/) { ($tb,$te)= @v[3,4]; } # gff 
    if($te>0) { 
      $tb= _max(1,$tb-$EXPAND_LOCUS); $te=$te+$EXPAND_LOCUS;
      $nlocus++;  $span += 1+$te-$tb;
      my $rloc= [$ref,$tb,$te]; 
      foreach my $ib (int($tb/$BINSIZE) .. int($te/$BINSIZE)) {
        push( @{$locuslist{$ref}{$ib}}, $rloc);
      }
    }
  } close(F);
  
  warn "#locuslist $locuslist restriction: n=$nlocus, span=$span\n"; #  if $debug;
  $locuslist="" unless($nlocus);
  return $nlocus;
}



sub collect_gff
{
  my($gff)= @_;
  
  my (%genebins, %genes, %exons, %sources, @geneorder); # return these
  ## drop some of these inner hashes:
  ## my (%cdslen, %exonlen, %inor, ); # dont need
  my (%promax, %proscore, %groupscore, %groupmax, %genescore, %idremap); # local
  
  my ($ng,$nr)= (0,0);
  my @sfields= manyscorefields(0);
  my @nulls = (0) x @sfields;
  my %binspan;
  my $skiplocus=0; # with locuslist
  
  ## special many.score fields: CDS, UTR=trsize/cds ratio, pro, nintron, trsize  

  ##my $nsfields = @sfields;
  my $icds= -2; #$nsfields;
  my $intr= -2; # $nsfields+1; # add ALWAYS? == KeyINSTRAND/intr 
  my $inexon= -2; # $nsfields+2; # add ALWAYS? diff from KeyINSTRAND/intr; see below
      # ^? add requred scores to sfields if not there? : CDS, intr, inex/nintron; want nintron/nexon for mrna
  my $ipro= -1;  my $iutr=-1; my $itrs=-1; my $iexons= -1;
  # my $inin= -1; # nintron >> now is inexon score
  # my $inrev= -1; #DROP  
  my $hasKeyNINTRON= 0; # FIXME; inexon clash user supplied and internal field
  my $ipcds= -2; #$nsfields; calculated PCDS = 100*c/xlen, but cut out >95
  
  # phom = homolog/selfho computed score; is this old ipro ?
  
  use constant NEWINSCORE => 1;  # for KeyINSTRAND
  
  for( my $i=0; $i < @sfields; $i++) { 
    $icds=$i if($sfields[$i] =~ /^CDS/i); 
    $ipcds=$i if($sfields[$i] =~ /^PCDS/i); 
    $ipro=$i if($sfields[$i] =~ /^pro/i); 
    $intr=$i if($sfields[$i] =~ /^$KeyINSTRAND/i);  
    #old# $inin=$i if($sfields[$i] =~ /^nint/i);  ## replace/use this for inexon score
    $inexon=$i if($sfields[$i] =~ /^inexon|^$KeyNINTRON/i);   # was /^nint|/
    # $inrev=$i if($sfields[$i] =~ /^inrev/i); #? DROP:expand to all intron errors: split w/ other gene; reversed; exon-internal ? 
    $itrs=$i if($sfields[$i] =~ /^trsize/i); 
    $iexons=$i if($sfields[$i] =~ /^nexon/i);   
    $iutr=$i if($sfields[$i] =~ /^UTR/i); # what weight? (exon-cds)/exon : minimize this? or what? 
        # utr >> cds is poor, but utr ~< cds is good; utr= (exon-cds)/cds ? (exon-cds)/exon ?
  }
      #  add requred scores to sfields if not there? : CDS, intr, inex/nintron; want nintron/nexon for mrna
  $icds   = add_required_scorefield("CDS",1) if($icds < 0);
  $intr   = add_required_scorefield($KeyINSTRAND,2) if($intr < 0);
  if($inexon>=0) { $hasKeyNINTRON=1; }
  else { $inexon = add_required_scorefield($KeyNINTRON,2); }
  # update @sfields, @nulls also  
  @sfields= manyscorefields(0);
  @nulls = (0) x @sfields;
  my $si=0; my %sfields= map{$_ => ++$si} @sfields;
  
  while(<$gff>){
    next unless(/^\w/); chomp;
    my $line= $_;
    my($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr)= split"\t";
    $nr++;
    my @proid=(); my $promax=0;
    
    $tscore=0 if($tscore eq ".");
    my $score= $tscore; # want this to select; w/ field choice?
    my $exonlen= 1 + $te - $tb;
    
    my($gid); 
    if($tattr =~ m/\bID=([^;]+)/) {  $gid=$1; } # fixme for exon ID=xxx;Parent=gene
    elsif($tattr =~ m/\bParent=([^;]+)/) {  $gid=$1; }
    unless(defined $gid) { $gid = "N".$ng; }
    
    if($scoretype =~ /^score/i) { $score=$tscore; } # also default
          
    elsif($manyscore or $scoretype =~ /^many.(.+)/i) {  # manyscore flag?
      my @scorev=(0) x @sfields;
      $score= ""; my $si= 0;
      foreach my $sf (@sfields) {
        my $v=0; my $sval ="";
        # my $sval = ($tattr =~ m/\b$sf=([^;\s]+)/) ? $1 : "";
        if($sf eq "score") { $v=$tscore; }
        elsif($tattr =~ m/\b$sf=([^;\s]+)/) { $sval = $1;
          # add computed pfield = 100*a/b  when sval= "a/b" ?
          if($sval =~ m/([\d\.e+-]+)/) { $v=$1; my $iv=$sfields{'p'.$sf}||0;
            #o# if($v =~ /e\+/ or $v>9 and $v=~/\./) { $v=int($v); } 
            if($v =~ /e\+/ or $v>9 and $v=~/\./) { my($dv)= $v=~m/(\d[\d\.\+e]+)/; $v=int($dv); } 
            if($iv>0 and $sval =~ m|[\d\.e+-]+/([\d\.e+-]+)|) {  
              my $vmax=$1; $scorev[$iv-1]= int( 0.5+100*$v/$vmax) if($vmax>0); }
          }
        } else {
          ## NOT next; # no $sf score : scorev[si] == 0;
        }
        
           # KeyINSTRAND = +fwd/-rev, or +fwd or -rev
           # NEWINSCORE meaning = +ok/-error where error score includes splice-reversed and exon-internal splice errs
           # value == sum of (ok - error), format is largest abs first (?)
        # sval == +ok/-err  or -err/+ok  largest abs first
        ## FIXME2: new intr= introntab inok/exons  NOT pos/neg strand
        if($sf =~ /^$KeyINSTRAND/ and $sval) { 
          my $in1 = ($sval =~ m,^([\d+-]+),)? $1 : 0;
          my $in2 = ($sval =~ m,/([\d+-]+),)? $1 : 0;
          if($in2) { $v= $in1 + $in2; } else { $v= $in1; }
        }
        
        #old# $score .= "$v,"; # pval where? replace this with scorev[]
        $scorev[$si]= $v if($sval); $si++; # new
        
        # this proid works but is a hack
        # pro=299/304,tribolium_TC013709,tribolium_TC014083
        if($sf =~ /^pro/i && $typ =~ /CDS/ && $sval =~ m/^([\d.e+-]+).([\d.e+-]+),([^;\s]+)/) {
          my($proaln,$prolen,$proids)=($1,$2,$3);
          #maybe: $scorev[$si-1]= $proaln; ## same as $v above?
          $promax= $exonlen; ## $prolen; # for this cds; this is just CDS length
          @proid = split",",$proids;
          }
          
      }
      $score = join(",",@scorev); # new
    }
    elsif($tattr =~ m/\b$scoretype=([^;]+)/) { $score=$1; }

# modified gene gff [rloc] record fields
# see above constant{ jCHR => 0, jSRC => 1, jTYPE => 2, jBEGIN => 3, jEND => 4, 
#               jSCORE => 5, jSTRAND => 6, jPHASE => 7, jATTR => 8, jGID => 9, 
#               jCDSLEN => 10, jTRLEN => 11, jSCORESUM => 12 };
    
    my $rloc= [$ref,$src,$typ,$tb,$te,$score,$to,$tph,$tattr,$gid]; # exon/mrna record
        # change to string and save mem ?
    
    if($typ =~ /^($mrnatypes)$/) { # $overtype == kBESTEXONS and 
      
      # check for nonuniq gene ids ? need to change for exon,cds also
      ## should drop not rename ??
      if(defined $genes{$gid}) {
        my $nid= "$gid.$ng";
        while(defined $genes{$nid}) { $ng++; $nid= "$gid.$ng"; }
        warn "#WARN: already have gene $gid => $nid\n" if $debug;
        $idremap{$gid}= $nid; # short term map for exons
        $gid= $nid;  $rloc->[jGID]= $gid;
      }

      $skiplocus=0;
      if($locuslist and not overlocuslist($ref,$tb,$te)) {
        $skiplocus=1;
        next;
      }
      
      $genes{$gid}= $rloc; $ng++; # should be uniq/id ?
      $sources{$src}++; # count all gene sources
      push(@geneorder, $gid) if($RESCORE_ONLY); # for rescore w/o removals?
      
      #* can we adjust bins so no genes overlap bounds? important, see above
      #.. need to collect/save adjustments, no assumption of sorted input yet
      # 2013sep: problems from huge-genespans picking exons that overlap later bins.
      #  switch from genebins{ref}{bin},rloc to exonbins{ref}{bin}{gid}  for each exon ?
      # bug here? undef $bins[0]      
      my @bins= (int($tb/$BINSIZE) .. int($te/$BINSIZE));
      
if($ONEBIN_GENES) {     

      # * this is a problem, not binning right
      my $ib= $bins[0] || 0;
      if(@bins>1) {
        my $ib1= $bins[1];
        #?? need to check all prior binned gene spans if change bins here?
        if(defined $genebins{$ref}{$ib1} ) {
          my @grefs= @{$genebins{$ref}{$ib1}}; 
          my @gref2= @grefs; my $upd=0;
          for( my $ir= $#grefs; $ir>=0; $ir--) { # need to go backwards for splice ?
            my $gb = $grefs[$ir]->[jBEGIN]; 
            if($gb > 0 and $gb < $te) { # move to this bin
              push( @{$genebins{$ref}{$ib}}, $grefs[$ir]); 
              splice(@gref2, $ir, 1); $upd++;
            }
          }
        $genebins{$ref}{$ib1}= \@gref2 if $upd; 
        }
      } else { # need to check if tb,te in lower bin; BUT unsorted input can change this
        if($ib > 0 and defined($binspan{$ref}{$ib - 1})
          and $binspan{$ref}{$ib - 1}{e} > $tb) { $ib = $ib - 1; } #??
        ##if($ib > 0 and $binspan{$ref}{$ib - 1}{e} >= $te) { $ib = $ib - 1; } #??
      }
      
      unless(defined($binspan{$ref}{$ib})) { 
        $binspan{$ref}{$ib}{b}= $tb; $binspan{$ref}{$ib}{e}= $te; 
      } else {
        $binspan{$ref}{$ib}{b}= $tb if($binspan{$ref}{$ib}{b} > $tb);
        $binspan{$ref}{$ib}{e}= $te if($binspan{$ref}{$ib}{e} < $te);
      }
      
      push( @{$genebins{$ref}{$ib}}, $rloc);
} else {      
      foreach my $ib (@bins) {
        push( @{$genebins{$ref}{$ib}}, $rloc); #? or store $gid
        }
}

    } elsif($typ =~ /^($exontypes)$/) { #$overtype == kBESTEXONS and 
      next if ($skiplocus);
      
      if(defined $idremap{$gid}) { 
        $gid= $idremap{$gid};  $rloc->[jGID]= $gid;
        }

      push( @{$exons{$gid}}, $rloc);
      ##if($typ =~ /^CDS/) {  $cdslen{$gid} += $exonlen; } ## 1 + $te-$tb;
      ##else { $exonlen{$gid} += $exonlen; } ##  1 + $te-$tb; 
      
      # FIXME for gene-joins of genescore for "pro"/protein w/ prot-ids, 
      #  separately score each protid, keep only best if many
      #  need more data from above ($tattr =~ m/$sf=(\d+)/)
      #?  -scoretype='many.pro:10/ID,est:5,...' 
      
      #* 2011/04/20: add all-exon support stats: count number of exons (cds too?)
      #  supported in this mRNA.  Use like genegroup, indicate which scores: intr esp.
      # add ?? -fullgene='intr:20,ref:1,est:1,pro:1,rseq:1' opt to score mRNA on how many exons are supt by those evd
      # .. esp. for introns, model w/ 5 of 5 exons w/ intron match scores >> 2 models same locus w/ 2,3 exon-intron matches
      # .. same logic doesn't apply to expression/prot scores as they dont join exons
      # .. should use this as balance to findjoins : high join score should not override high exon-intron score

      if(1) {
        my @sc= split",", $score; # in @sfields order
        
        if($sumgenescore) {
          $genescore{$gid}= [@nulls] unless(ref $genescore{$gid});
          foreach my $i (0..$#sc) { 
            # intr == KeyINSTRAND = +fwd/-rev, or +fwd or -rev
            #old# if($i == $inexon and not $hasKeyNINTRON) #?? but not user supplied $KeyNINTRON < wrong
            if($i == $inexon) # 201208
            {  
              my($vintr,$vinexon,$isnin)=(0,0,0);
              if( UPD1801 ) { # upd1801: try both inexon, intr from exon annot .. later is used now?
              if( $hasKeyNINTRON ) { $vintr= $sc[$inexon]; $isnin=1; }
              unless($vintr) { $vintr = $sc[$intr]; $isnin=0; }
              if($vintr) {
                $vintr =~ s,/.*,,;
                $vinexon= ( $isnin ) ? $vintr : ($vintr > 0) ? 1 : ($vintr < 0) ? -1 : 0;
              }
              
              } else { # old:2014..
                $vintr= ( $hasKeyNINTRON ) ? $sc[$inexon] : $sc[$intr]; 
                $vintr =~ s,/.*,,; 
                $vinexon= ( $hasKeyNINTRON ) ? $vintr : ($vintr > 0) ? 1 : ($vintr < 0) ? -1 : 0;
              }
              $genescore{$gid}->[$i] += $vinexon ;   # should use Weight of 100 * nin to make comparable to base scores
            } else {
              $genescore{$gid}->[$i] += $sc[$i];  
            }
          }
        }
        
        if($typ =~ /^CDS/ and $ipro >= 0) { 
          foreach my $proid (@proid) { 
            $proscore{$gid}{$proid} += $sc[$ipro]; 
            $promax{$gid}{$proid} += $promax; # do need? not same as $cdslen{$gid} ??
            }
          } 
      
        if($genegroup) {
          foreach my $i (0..$#sfields) { 
            if($genegroup[$i]) {
            my $sf= $sfields[$i];
            #no,fixme1809: genegroup == ipro, use @proid, proscore, promax per above CDS .. no effect? below gets same vals
            # if($i == $ipro) {
            #   foreach my $proid (@proid) { 
            #     $groupscore{$i}{$gid}{$proid}= $proscore{$gid}{$proid}; 
            #     $groupmax{$i}{$gid}{$proid}=  $promax{$gid}{$proid};
            #     }
            #   next; # i
            # }
            my $sval = ($tattr =~ m/\b$sf=([^;\s]+)/) ? $1 : "";
            if($sval =~ s/^\d[^,]+,//) { 
              my $pid=$sval; 
              my @pid= split/[,\s]/,$pid;
              map{ s,/.*$,, } @pid; # fixme2:  ID/score,ID/score syntax
              foreach $pid (@pid) { 
                $groupscore{$i}{$gid}{$pid} += $sc[$i]; 
                $groupmax{$i}{$gid}{$pid} += $exonlen; #?? need this
                }
              }
            }
          }
        }
      
      }
      
    }
  }  # end gff
  
  
  ## make sub set_manyscore
  #  manyscore_set( \%genes, \%exons, \%genescore, \%groupscore, \%groupmax, \%proscore, \%promax)

  if($manyscore) { # $scoretype =~ /^many/
    foreach my $gid (sort keys %genes) {
    
      my( $cdsw, $exw, $nexon, $ncds)= (0) x 10;
      foreach my $ex (@{$exons{$gid}}) {
        my $w= 1 + $ex->[jEND] - $ex->[jBEGIN];
        if($ex->[jTYPE] =~ /CDS/) { $ncds++; $cdsw += $w; }
        else { $nexon++; $exw += $w; }
      }      
      ## $nin-- if( $nin>0 );
      
      #  in $genes{$gid}= $rloc for use above ??
      # exon= [$ref,$src,$typ,$tb,$te,$score,$to,$tph,$tattr,$gid] n=10
             
      ## REVISE UTR: use option utr_expected_size (in bases, median value)
      ## score: (utr_size == utr_expected_size) size * +weight;  
      #   or (utr_size <<>> utr_expected_size)  (size - expected) * -weight
      
      my $utr = 0;
      if($iutr >= 0) { #if(UTR_USE_NEW_WEIGHT)  

      my $nutr= $nexon - $ncds;
      $utr= _max($exw - $cdsw, 0); # $utr= 0 if($utr < 0);
      # change back to ratio utr/exw ?? or use that to moderate score
      
      # utr_best = 400/1600 // $utr/$exw 
      my $UTR_BEST_SIZE= 0.5 * $UTR_TOO_BIG; ## if($UTR_BEST_SIZE >= $UTR_TOO_BIG);      
      my $utr_toobig = $UTR_TOO_BIG * $exw; # change $UTR_TOO_BIG to this 0.50
      
      my $utr_best   = $UTR_BEST_SIZE * $exw;  #? UTR_BEST_SIZE option
      ## 2013^ this is problem, small % utr, large % cds is okay, need bp size test for long trs.
      ## for average 2000b exw, pUTRBEST=0.30, utrbest=600b ? too big to expect? need some data
      if($exw >= 1000 and $nutr>1) { $utr_best= _min(150,$utr_best); } # fixed bp? only care if utr small in bp, or missing one end.
      
      my $utr_toomany = ($nutr > 4) ? 1 : 0; # need left,right utr counts : look for utrx=5,2 (l,r) field
      ## my $utr_toofew = ($nutr < 2)?1:0;
      
      # if nutr very big, drop this model ?? only if no alt models at locus
      # but some are valid ncrna
      #  680  bestgenes.DGILmix8c.utrtoomany.tab  : BEFORE utr_toomany score adjust
      #  353  bestgenes.DGILmix8d.utrtoomany.tab  : AFTER
      
      my $utrscore= 0;
      if($utr <= $utr_best) # change utr_best, okay up to 95%+ as long as some 100bp exist?
      {
        $utrscore = $utr; ## == $utr_expected_size * $utr/$utr_expected_size; # positive
        $utrscore -= $utr * ($nutr - 3) if($utr_toomany);
      } 
      elsif( $utr >= $utr_toobig ) 
      { # too large, neg weight?
        $utrscore = ($utr_toobig - $utr) * 0.5; # * negative
        $utrscore -= $utr * ($nutr - 3) if($utr_toomany);
      } 
      else # > best < toobig # near right; probably should set score to ~1 for all in okay-range
      {  
        $utrscore = $utr_toobig - $utr; # positive, smaller utr > +score
        $utrscore -= $utr * ($nutr - 3) if($utr_toomany); # make negative?
      }
      $utr= int($utrscore);
      
# } else { # old utr proportion of exon
#       # reduce long utr more; 40% is best?  
#       $utr = int(100 * ($exw - $cdsw) / $exw); # see also PCDS
#       $utr = ($utr >= 40) ?  60 - $utr : $utr; # gives -score for > 60% utr; ok?
# }
      } # iutr
      
      my $rv = $genes{$gid}; 
      $rv->[jCDSLEN]= $cdsw; # this is same as $rv->[5]{CDS} field
      $rv->[jTRLEN]= $exw; # "hidden" attrs; add for above use
      
      my @sc = split",", $rv->[jSCORE]; # score
      $sc[$icds]= $cdsw if ($icds >= 0); 
      $sc[$iutr]= $utr if ($iutr >= 0); 
      $sc[$itrs]= $exw if ($itrs >= 0);
      $sc[$iexons]= $nexon if ($iexons >= 0);

      #NOT NOW: $sc[$inexon]= $nin if ($inexon >= 0); # was inin; change this scoring: num exons w/ +intron evd
      # done below: mrna.$sc[$inexon] == $genescore{$gid}->[$inexon]
      # ? score inexon= -1 / 0 for case where intron hits middle of exon?
      my $validin = 0;
      $validin = $sc[$inexon] || $genescore{$gid}->[$inexon] || 0; # hasKeyNINTRON : have mRNA nintron= scored
      # if($hasKeyNINTRON and $tattr =~ m/$KeyNINTRON=(\d+).(\d+)/) { ($validin,$nexon)= ($1,$2); } 
      unless( UPD1801 ) { # now validin score should == number of valid introns
      if($hasKeyNINTRON and $validin) { $validin= int(5*$validin)/10; } # dang, 1/2 validin no longer right, not splices but ints
      }
      
      my $score_inexon = ($nexon >= 1) ? ";inexon=$validin/$nexon/$ncds" : ""; # should be n exon, not n intron
      # $sc[$inexon]= $validin if ($inexon >= 0);
      # ?? add ncds somewhere, inexon=nin/nex/ncds ?
      
      my $proid="";
      my $addat= "";
      
      # add scorable calc field, alt to UTR .. modify when cdsw == exw, no UTR ?
      my $pcds="";
      if ($exw>0) {
        my $PCDS=int(100*$cdsw/$exw); $pcds= ",".$PCDS."%";
        $PCDS= ($PCDS>100) ? 0 : ($PCDS>95) ? 100 - $PCDS : $PCDS; # cds > exw ?
        $sc[$ipcds]= $PCDS if($ipcds >= 0);
      }
      $addat .=";cxlen=$cdsw/$exw$pcds"; #want cdsw,exonw in addat?
      $addat .= $score_inexon;
      
      if($genegroup) {
        foreach my $i (0..$#sfields) { 
          if( $genegroup[$i] and defined $groupscore{$i}{$gid}) {
          my ($pid1)= sort{ $groupscore{$i}{$gid}{$b} <=> $groupscore{$i}{$gid}{$a}
                          } keys %{$groupscore{$i}{$gid}};
          if($pid1) {
          my $sf= $sfields[$i];
          my $gsc= $groupscore{$i}{$gid}{$pid1};
          my $gmax= $groupmax{$i}{$gid}{$pid1}; # $gmax=$gsc if($gmax < $gsc); #?? FIXME
          $addat .= ";best_$sf=$gsc/$gmax,$pid1";  #?? need score/size,pid ? size=cdsw or exon
          # best_rna, best_pasa, best_pro, best_ref, ...
          }
          }
        }
      }
      
      if($sumgenescore) { # make default ??
        my $gv= $genescore{$gid} ; 
        if(ref($gv) and @$gv>0) {
          foreach my $i (0..$#sc) { 
            # careful ; have -score now
            # genescore/gv == exon sum of scores;  sc here == mRNA scores
            # what of special scores above: icds, iutr, ..
            next if($i == $icds or $i == $iutr); ## or $i == $inexon
            my $gv= $gv->[$i] or next;  
            if($sc[$i]) { $sc[$i] += $gv;}  # mRNA or exon-sum ?? sum both?
            else { $sc[$i]= $gv; } ##NO if ($gv > $sc[$i]);  
            }
        }
        
        # generalize this to any evd with gene_span_ID: EST/pasa-asm, refseq, rnaseq/cuff
        # score as greatest cover/evdID, downweight preds >> evdID genespan as well as under evdID
        # i.e. best score = 100% match to evdID exons, reduce for over/under
        # designate in scoretype ? or add -genegroup='...' as per augustus.cfg:1group1gene
        # use 2 ways: 1. CDS join/split when evdID cds <<>> model cds
        #     2. UTR aberrant when evdID cds == model cds, evdID UTR > 0, but <<>> model UTR
        # use also w/ JOINSCORES / findjoins
        
        if(defined $proscore{$gid} and $ipro >= 0) {
          my @pk= sort{$proscore{$gid}{$b} <=> $proscore{$gid}{$a}} keys %{$proscore{$gid}};
          $proid= $pk[0]; # what if 2+ same score?
          my $pv= $proscore{$gid}{$pk[0]}; # max score per pro gene
          my $px= $promax{$gid}{$pk[0]}; # just sum CDS length
          $sc[$ipro]= "$pv/$px";
          $addat .= ";pro1=$pv/$px,$proid";
        }

###? expand to all intron errors: insplit (w/ other gene); proper+reversed = join; reversed-only; exon-internal
###  ;inerr=count/type,count2/type2,... 
#     unless(NEWINSCORE) {
#         if(defined $inor{$gid}) { #  and $inrev >= 0
#           my $fwd= $inor{$gid}{1} || 0; 
#           my $rev= $inor{$gid}{-1} || 0;
#           my $geneor= $rv->[6]; # + or -
#           my $add="";
#           if($geneor eq "+" and $rev) { $add=abs($rev); }
#           elsif($geneor eq "-" and $fwd) { $add=abs($fwd); }
#           
#           $sc[$inrev]= $add if($add and $inrev);
#           $addat .= ";inrev=$fwd/$rev" if($add); #?? add only if both or differ from mRNA strand? inrev=
#         }
#     }        

      }
      
# see above constant { jSCORE => 5, jATTR => 8, jCDSLEN => 10, jTRLEN => 11, jSCORESUM => 12 };
      ## add two attributes: combined score: scoresum=, best proscore ID: proid=
      my $mscore= $rv->[jSCORE]= join(",",@sc);
      my $msum= manyscore_sum($mscore);  #* also store in $rv->[NNN] for sort use
      $addat .= ";scoresum=$msum";
      $rv->[jATTR] =~ s/;(pro1|scoresum|cxlen|inexon)=[^;\s]+//g; # drop any old addat
      $rv->[jATTR] =~ s/$/$addat/; #  == tattr;
      $rv->[jSCORESUM]= $msum;


#  ** revise -dropscore to allow AND / OR / NOT?, e.g.
#    -dropscore= CDS:180 AND ( est:60 or pro:60 or rseq:60) AND NOT terepeat:60
#  .. but need nesting parenth parsing?  ( (a b c) AND (x y z) ) NOT (p OR q)
#  .. maybe this syntax:  a or b or c +CDS:nnn  -TE:nnn -inerr:nnn
#      + = always required, - = never allowed, otherwise require a OR b OR c
#  FIX2: mustkeep should trump dropscore

      if($dropscore) {
        # my $keep=0;
        my $tkeep= ($mustkeepdrop) ? $mustkeepgene{$gid}||0 : 0;
        my $keep= ($tkeep == kKEEP_REQUIRED || $tkeep == kKEEP_ALTERNATE)?1:0;
        if($keep) { } # continue
        elsif(USE_DROPSCORE_ANDOR) {        
        for(my $i=0; $i < @sfields; $i++) { 
          my($ds,$sv)= ($dropscore[$i],$sc[$i]);
          my $droptype= $droptype[$i] || 0; # reorder @sfields to put require/AND, reject/NOT at end
          
          #o#map{s,/.*,,}($ds,$sv); # "$pv/$px" # Argument "72/232" isn't numeric in numeric 
          map{s=[/,;].*==}($ds,$sv); # "$pv/$px" # Argument "72/232" isn't numeric in numeric 
          #^ add other syms: ,;
          
          if(not defined $ds or $ds == 0) {
          
          } elsif($droptype == kDROP_AND) {
            if ($sv < $ds) { $keep= -1; } 
            #No.off1208# elsif ($sv >= $ds and $keep >= 0 ) { $keep=1; } 
# ** FIXME: drops not dropping all with no evd.
# * problem is dropscore: +CDS:201 / kDROP_AND: 
# - need to say that CDS alone is not a keeper, but need at least CDS:201 to keep
# - ie dont set keep=1 for CDS>=201, only set keep= -1 if CDS<201
# - likely need both ways: +nintron:2 means must have >=2 but also is sufficient to keep; is this MUSTNOT ?
          
          } elsif($droptype == kDROP_NOT) { #? is this useful
            $keep= -1 if ($sv >= $ds);           
           
          } elsif($droptype == kDROP_MUSTNOT) { # add
            do{ $keep=1; last} if ($sv >= $ds);           
           
          } else { # kDROP_OR
            $keep= 1 if ($sv >= $ds and $keep >= 0);           
          }
          
        }
        $skipgene{$gid}=kSKIP_NOEVD unless($keep > 0); # 999 = flag dont unskip this one
          
      } else {  # original, OR only
        for(my $i=0; $i < @sfields; $i++) { 
          my($ds,$sv)= ($dropscore[$i],$sc[$i]);
          map{s,/.*,,}($ds,$sv); # "$pv/$px" # Argument "72/232" isn't numeric in numeric
          # want -dropscore?  terepeat:-100 means drop if >=100 TE bases 
          $keep=1 if ($ds > 0 and $sv >= $ds); 
          }
        $skipgene{$gid}=kSKIP_NOEVD if($keep==0); # 999 = flag dont unskip this one
      }
      }

      
    }
  }
  
  my $nsrc= scalar(keys %sources);
  warn"#collect_gff=$nr, ngene=$ng, nsource=$nsrc\n" if $debug;
  return (\%genebins, \%genes, \%exons, \%sources, \@geneorder);
}



__END__


