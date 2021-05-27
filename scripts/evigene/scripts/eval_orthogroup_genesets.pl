#!/usr/bin/env perl
# eval_orthogroup_genesets.pl
# rename? evorthogenes.pl

=item about

Evigene script to evaluate gene sets for match to orthology groups (orthomcl)
Inputs: geneset-orthospecies.blastp, orthogroup.tables (of orthomclsum.pl)

  -- use geneset-i x omcl-species blastp results;
     collect omcl-species ids w/ omcl group oid, skip 1-species groups (or require nt-min>2 species/group?)
     find max,min,mean scores for geneset-i x omclid  (scores=bits, ident?, align? == .tall table)
        - assign geneset-i to only best matching group
        - exclude geneset-i if < min-score of omcl group?

  genesets x orthology group stats
  ** bug now, errors where OMCL-lowqual-match to tg preceeds OMCL-himatch-tg ; need best-tg x best-arp match
  ** do this for all species genes in wasp-arp11 spp set? and cacao, aphid, ...
  ** add count/calc of too-long, too-short genes in groups, using median aa-size of group
    >> this works BUT extreme cases (e.g. frave) throw off excess min/max counts for other species; add new filter?
    >> add reverse of -skips=xxx,yyy : -keeps=aaa,bbb,ccc
  stats for geneset-i x omcl groups:
    1. count ngroups found in geneset-i (2 levels? any bestmatch, bestmatch >= median/min group score)
    2. ave. bitscore (w/ and w/o groups common to all genesets); ident, align scores close to bitscore
    3. count/ave too-long, too-short aa in geneset-i, using group median aa-len, too=+/- 2*stdev ?

=item inputs

  orthomcl tables:
    arp11u11-orthomcl-count.tab : group gene counts/species table : select groups w/ species >= mintaxa
    arp11u11_omclgn.tab         : group oid, geneid table
    faCount arp11u11noalt.aa.gz : geneid, aasize table for orthogenes
  
  genesets x orthogenes blastp
  
  
  evorthogenes.pl  -tallscore evg11u.arp11.tall4 -mintaxa 3 -skips wasp \
   -groupgene arp11u11_omclgn.tab -groupcount arp11u11-orthomcl-count.tab -geneaa arp11u11noalt.aa.gz 

=item usage2

  see evigene/scripts/omcl/omcltoptab.sh
  call these:
  env skipho=$skipho aa=$aa tall=1 $evigene/scripts/makeblastscore.pl $blastp > $onam.tall4 
  
  $evigene/scripts/eval_orthogroup_genesets.pl -bitmed -mintaxa 3 -skips $skips \
   -out $onam.topout2 -tallscore $onam.tall4 \
   -groupgene $omcl/*_omclgn.tab -groupcount $omcl/*-orthomcl-count.tab -geneaa $omcl/*.aa.gz
  
  #  summarize ortho scores for several species
  cat xxx*.topout2 | $evigene/scripts/eval_orthogroup_genesets.pl -intab stdin
  


=item run1

$evigene/scripts/eval_orthogroup_genesets.pl -mintaxa 3 -skips wasp -out evig11u.topout -tallscore ge
nes/pub11u/aaeval/evg11u.arp11.tall4 -groupgene $warp/arp11u11_omclgn.tab -groupcount $warp/arp11u11-orthomcl-cou
nt.tab -geneaa $warp/arp11u11noalt.aa.gz

#i orthogroups from prot/omcl9u11//arp11u11-orthomcl-count.tab
#i group n=10477 with ntaxa>=3
#i geneaa n=204423
#i groupgenes n=98129
#i top groups for genes/pub11u/aaeval/evg11u.arp11.tall4 to evig11u.topout
#i top groups n=9269 for ngenes=410164

Summary stats for evig11u.topout   << stats differ due to 1 geneset
   Nasvi: n=8938; nreal=8938; abits=674.1; dlen=-12.9; arealbits=674.1; xtralen:xl-2=269(3%);  xl2=103(1.1%); 
# orig
# NasviEg: n=9159; nreal=8940; abits=658.1; dlen=-11.9; arealbits=674.2; xtralen:xl-2=269(3.0%);  xl2=104(1.1%); 


set spp=arath ; set pt=plant9-$spp ;  $evigene/scripts/eval_orthogroup_genesets.pl -sppidmap "Thecc=cac
ao" -mintaxa 3 -skips $spp -out $pt.topout -tallscore $pt.tall4 -groupgene $carp/plant9_omclgn.tab -groupcount $ca
rp/plant9-orthomcl-count.tab -geneaa $carp/../plant9.aa.gz

=item input orthodb / evgeney zbod... ?

  #20121009: ? revise this to input orthodb groups, geneaa ?
  also check/compare w/ orthomcl scores for equivalent gene sets? esp. look for missed cases
    -- orthodb includes short aa that orthomcl kicks out before MCL clustering,
      ie fragments of 1/2 size/score.  orthodb groupings are wider than MCL defaults, but maybe similar
      after MCL group clustering?
  
      >> have this data: groupid, counts, aa
  my( $ogenegrouph, $ogenesizeh)= orthodb_groupgenes($groupcounts,$geneaa,$genegroups);
  
      >> missing genescore from orthodb; punt? sum only aa sizes/group, and size deviants?
      >> or compute blastp, per odb groupings?
  $topout= topmcl($genescore_tall, $ogenegrouph, $ogenesizeh,$topout); ##, $outname);

=cut

use strict;
use Getopt::Long;
use File::Basename;

my $DEBUG=1;
my $AA_TOOBIG = 1.95;    # instead of 2 or 3 sd; add 2 levels? 1.3, 1.9
my $AA_TOOBIG1 = 1.45;    # instead of 2 or 3 sd; add 2 levels? 1.3, 1.9
my $AA_TOOSMALL = 0.60;  # 0.65; for 1 level # 2 levels?  0.75?, 0.65?
my $AA_TOOSMALL1 = 0.75; # 2 levels?  0.75?, 0.65?
my $AA_OUTLIER1= 0;

# per run configs
# my $phyla="arp11u11";
# my $spl="antc,anth,aphid,apis2ref,bombusimp,bombusterr,daphnia,drosmel,human,trica,wasp";  
my $idprefix="OID";
my $skipspecies="";
my $keepspecies="";
my $FixSpeciesIDt2p="wasp"; # see fixids()
my $keeptarget=""; 
my $NAMWID=15; # opt    

my $MINTAXA=3;
my $USE_MEDBITS=0; # default?
my $DOALL= 0;
my ($genescore_tall, $genegroups, $groupcounts, $geneaa, $topout, $intab, $sppidmap);
my $digits=1;
my $SCORECOL= 0; # 0=bitscore, 1=ident, 2=align from tall4 tables
my $ALLTARG= 0;     #  -alltarg, put all tg, using keys $arbm{$ar}{$tg} ?? or keys  $btx{$tg}
my $DOSWAP=0;
my $DOCOM= 1; my $DOREFPCT= 1; my $DOTABTINY=0;

my $optok= GetOptions(
##  "config=s", \$config,
  "tallscores=s", \$genescore_tall,
  "intab=s", \$intab,
  "output=s", \$topout,
  "groupgenes=s", \$genegroups,
  "groupcounts=s", \$groupcounts,
  "geneaa=s", \$geneaa,
  "skipspecies=s", \$skipspecies,
  "keepspecies=s", \$keepspecies,
  "keeptarget=s", \$keeptarget,
  "sppidmap=s", \$sppidmap,
  "digits=s", \$digits,
  "idprefix=s", \$idprefix,
  "FixSpeciesIDt2=s", \$FixSpeciesIDt2p,
  "MINTAXA=i", \$MINTAXA,
  "SCORECOL=i", \$SCORECOL,
  "bitmedian!", \$USE_MEDBITS, # group median bitscore rather than highest pair score in group
  "swap!", \$DOSWAP, # swap targ,ref cols 1,2
  "common!", \$DOCOM, 
  "refpct!", \$DOREFPCT, "tabtiny!", \$DOTABTINY, 
  "tiny1|big1!", \$AA_OUTLIER1, # 
  "allgene!", \$DOALL, # all gene or top gene/spp/group for topstats only
  "alltargetgene!", \$ALLTARG, # all gene for topout table 
  "debug!", \$DEBUG, 
  );


die "usage: eval_orthogroup_genesets.pl  -tallscore xxxx.tall -groupgene xxxx -groupcount xxx -geneaa xxxx 
  opts: -mintaxa $MINTAXA -idprefix $idprefix ... \n"
  unless($optok and ($genescore_tall or $intab)); ## and -d $orun and -f "$orun/all_orthomcl.out"); #  and $orun

# evigene_config($config, \@configadd); # always even if $config null

$skipspecies =~ s/[,\s]+/|/g if($skipspecies);
$keepspecies =~ s/[,\s]+/|/g if($keepspecies);
$SCORECOL=0 if($SCORECOL>2); # Bits, Iden, Align choices 0,1,2 in tall4 tables

my(%sppidmap,@sppidkey);
if($sppidmap) {
  map{ my($k,$v)=split /[=:]/, $_; $sppidmap{$k}=$v if($v); } split /[, ]/, $sppidmap;
  @sppidkey= sort keys %sppidmap;
}

if($intab) {
  topstats($intab);
} else {
  #20121009: ? revise this to input orthodb groups, geneaa ? see eval_orthodb_groups.pl
  my( $ogenegrouph, $ogenesizeh)= orthogroupgenes($groupcounts,$geneaa,$genegroups);
  $topout= topmcl($genescore_tall, $ogenegrouph, $ogenesizeh,$topout); ##, $outname);
  topstats($topout);
}


#........................ subs .........................

=item FIXME

** add species prefix add/remove
  cacao_Thecc1EG002561t1 == cacao:Thecc1EG002561t1 == Thecc1EG002561t1
  omcl tables have cacao_Thecc or cacao:Thecc
  geneaa has >Thecc
  tallscore has Thecc ...
  
=cut  

sub fixids
{
  local $_= shift;
  s/_/:/ unless(/:/); #  $id=~s/:/_/; 
  
  if($sppidmap) {
    # if(/(\w+):/) { my $s=$1; if($sppidmap{$s}) { } }
    unless(/^\w+:/) { 
      ## my @sp= sort keys %sppidmap;
      my $id=$_; my $idt="";
      # my($pr)= m/^([^\d\W]+)/; # look for spp tag from id start
      foreach my $sp (@sppidkey) { if($id =~ m/^$sp/) { my $t=$sppidmap{$sp}; $idt="$t:$id" if($t); last; } }
      $_= $idt if($idt);
    }
   } 

  ## s/t(\d+)$/p$1/ if($FIXID_T2P); # per species??
  s/t(\d+)$/p$1/ if($FixSpeciesIDt2p and m/$FixSpeciesIDt2p/); # damn this fixup crap
  # if(/aphid:gi/) { s/gi\|\w+.ref./; s/\|//;  } 
  # if(wantarray) { return split(":",$_,2); }
  return $_;
}

sub medmeansd {
  my($arr)= @_;
  my ($n,$sm,$ss,$md,$mn,$sd,$min,$max)=(0) x 9;
  my @sl= sort{$b<=>$a} @$arr; # big first
  $n=@sl;  ($max,$min)=  @sl[0,-1];
  return($n,$md,$mn,$sd,$min,$max) if($n<1);
  $md= @sl[int($n/2)]; 
  return($n,$md,$mn,$sd,$min,$max) if($n<3);
  for my $i (0..$n-1) { my $x= $sl[$i]; $sm+=$x; $ss += $x*$x; }
  ##$mn= $sm/$n; $sd= sqrt(($ss - $sm*$mn)/($n-1));
  $mn= $sm/$n; my $var= (($ss - $sm*$mn)/($n-1)); $sd=($var>=1)?sqrt($var):0;
  return($n,$md,$mn,$sd,$min,$max);
}

sub median {  # screwy results with 2 entries, got smallest
  my($arr)= @_;
  my ($n,$sm,$ss,$md,$mn,$sd,$min,$max)=(0) x 9;
  my @sl= sort{$b<=>$a} @$arr; # big first
  $n=@sl; ($max,$min)=  @sl[0,-1];
  $md= $sl[int($n/2)]; # if($n>0); 
  return wantarray ? ($md,$n,$min,$max) : $md;
}


=item topmcl : tall2top problem, wrong bitscore : -bitmed is wrong to compare two genesets on same protdb

 -- problem shows comparing two topmcl tables, for same ogroup match, same best match id, diff bitscore..
 -- doesnt happen often, dont see pattern.. 
 -- ** USE_MEDBITS effect **
    median bitscore from all matches, which can differ b/n two gene sets
    even though protein in question is same for both gene sets...

using -bitmed
daphmag12r:r9vb9ptf062k87L17t1    locust:1Sv1K31L3776t1   ARP603  3251    2154    2155    12
daphmag11a:m8AUGepir2s02581g199t1 locust:1Sv1K31L3776t1   ARP603  1524    2154    2155    13

using -nobitmed
daphmag11a:m8AUGepir2s02581g199t1 locust:1Sv1K31L3776t1   ARP603  3378    2154    2155    13
daphmag12r:r9vb9ptf062k87L17t1    locust:1Sv1K31L3776t1   ARP603  3378    2154    2155    12
 
grep m8AUGepir2s02581g199t1 arpx3-daphmagna_201104m8.topout2d 
daphmag2011m8AUGepir2s02581g199t1       locust:1Svel1K31L3776t1 2144    ARP603  
  1524 << wrong bits    1666    2154    2155    2142    13      0
        .. other scores from locust:1Svel1K31L3776t1 are right.
        
from this tall4 input table:
 grep m8AUGepir2s02581g199t1 arpx3-daphmagna_201104m8.tall4 
daphmag2011m8AUGepir2s02581g199t1       daphniaplx:hxAUG26res47g77t1    4244    2107    2155    2155
daphmag2011m8AUGepir2s02581g199t1       locust:1Svel1K31L3776t1 3378    1666    2154    2155
                                                                ^^ should be
daphmag2011m8AUGepir2s02581g199t1       trica:D6WBG8_TRICA      3372    1670    2165    2155
daphmag2011m8AUGepir2s02581g199t1       aphid:acyp2eg0021368t1  3352    1629    2165    2155
daphmag2011m8AUGepir2s02581g199t1       wasp:Nasvi2EG022168p1   3314    1620    2161    2155
daphmag2011m8AUGepir2s02581g199t1       ixodes:ISCW018611-PA    3298    1591    2165    2155
daphmag2011m8AUGepir2s02581g199t1       human:NP_054733.2       3288    1589    2159    2155
daphmag2011m8AUGepir2s02581g199t1       drosmel:NP_648818.3     3287    1613    2169    2155
daphmag2011m8AUGepir2s02581g199t1       zebrafish:NP_001116729.1        3251    1573    2159    2155
daphmag2011m8AUGepir2s02581g199t1       daphniaplx:hxAUG26us80g22t1     1536.5  814     1711    2155
daphmag2011m8AUGepir2s02581g199t1       human:NP_006819.2       1524    861     1795    2155
                                                                ^^ got this val
daphmag2011m8AUGepir2s02581g199t1       locust:1Svel1K23L4089t2 1508    936     2009    2155
daphmag2011m8AUGepir2s02581g199t1       tetur:s09g05980 1504    881     1279    2155

=cut

sub topmcl
{
  my($genescores, $ogenegrouph, $ogenesizeh, $outname)= @_;

  my(%togtab, %arlm, %arbm, %bax, %baxt, %baxs, %btx, %btxa, %btxs, $nod, $ngn);
  unless($outname) { ($outname=$genescores) =~ s/\.\w+$//; $outname .=".topmcl3"; }

  sub putbx { # in topmcl, uses local vars ..
    my($tg,$sg,$sl,$ar,$bt)= @_;
    push @{$arlm{$ar}}, $sl;  
    push @{$arbm{$ar}{$tg}}, $bt;  # do $bits,$ident,$aln ? NEED {ar}{tg} here .. per target gene, all bits
    # fixme maybe : mark if b(tg.sg) same as max btx/bax ? dont need for group tests, but for best-species tests
    if($bt>$bax{$ar}) { $bax{$ar}=$bt; $baxt{$ar}=$tg; $baxs{$ar}=$sg; }  
    if($bt>$btx{$tg}) { $btx{$tg}=$bt; $btxa{$tg}=$ar; $btxs{$tg}=$sg; } 
  }


  # grep -v wasp prot/omcl9u11/arp11u11_omcl11.gr3.gids | cat - genes/pub11u/aaeval/evg11u.arp11.tall4 |
  # env gtag=ARP perl -ne 

=item old vs new .tall4
  # FIXME: 2013oct: .tall4 has same-species gene rows, same-gene rows : skip those here.. remove from tall4? no
old: kfish1/../fish8-medaka.tall4.gz 
Query   Source  Bits    Ident   Align   Qlen
medaka:ENSORLP00000000001       killifish:Funhe5EG035383t1      149.1   101     212     441
medaka:ENSORLP00000000001       killifish:Funhe5EG040602t1      140     80      153     441

new kfish2/../blastz/fish11-medaka.tall4
Query   Source  Bits    Ident   Align   Qlen    Slen
medaka:ENSORLP00000000001       medaka:ENSORLP00000000001       901     441     441     441     441
medaka:ENSORLP00000000001       medaka:ENSORLP00000025275       208     116     154     441     171
medaka:ENSORLP00000000001       kfish2:Funhe2EKm006518t1        171     155     394     441     466
medaka:ENSORLP00000000001       kfish2:Funhe2EKm017828t1        152     97      200     441     433

=cut
  
  warn "#i top groups for $genescores to $outname\n" if($DEBUG);
  if($genescores =~ /\.gz/) { 
    open(INTALL,"gunzip -c $genescores |") or die "$genescores"; # handle .gz
  } else {
    open(INTALL,"$genescores") or die "$genescores"; 
  }
  while(<INTALL>){
    if(/^Query/){ next; }
    elsif(/^\w/ and /\t/) { 
      chomp; my($tg,$og,@v)= split"\t"; 
        # @v SHOULD be Bits,Iden,Aln,Tlen ; Tlen may be 0/missing
        # tall4 adds Slen col: Query Source Bits Ident Align  Qlen Slen
      next if($tg eq $og);
      
      if($DOSWAP) {
        ($tg,$og)=($og,$tg); @v[3,4]= @v[4,3]; #tlen,olen
      }
      
      $og= fixids($og); # and??  $ts= fixids($tg);
      my($ts,$os)=map{ my($s)=split ":"; $s; } ($tg,$og);
      next if($ts eq $os);
      
      
      ## add keepspecies/skipspecies to ts? os filtered in orthogroup sub
      ## need only 1 spp in keeptarget to build proper tables..
      ## need eq targ !! damn
      #??#next if($keeptarget and ($ts ne $keeptarget)); # both keep,skip?
      next if($keeptarget and not($ts =~ m/^($keeptarget)$/)); # both keep,skip?
      # next if($skipspecies and $ts =~ m/$skipspecies/);

      my $oid= $ogenegrouph->{$og};
      next unless(defined $oid); # allow oid == 0

      my $tlen= $ogenesizeh->{$tg}||0; 
      $v[3]= $tlen if(@v<4 or $v[3] == 0);
      
      my $olen=$ogenesizeh->{$og}||0; # || $v[4] ?? if aasizes missing but tall4 has them? 
      my $poid= ($oid =~ /^\d/) ? $idprefix.$oid : $oid;
      # print OUT join("\t",$tg,$og,$olen,$poid,@v)."\n"; # option
      $togtab{$tg}{$og}=join"\t",$olen,$oid,@v; # bav{}
      # OPTION: swap bitscore/v[0] for ident/[1] or align/v[2]
      my $score=$v[$SCORECOL];
      putbx( $tg, $og, $olen, $oid, $score); #? do all @v
      $ngn++; 
    }
  } close(INTALL);
  # > genes/pub11u/aaeval/evg11u.arp11.tomcl
  $nod= scalar(keys %bax);
  warn "#i top groups n=$nod for ngenes=$ngn\n" if($DEBUG);
  
  #  > genes/pub11u/aaeval/evg11u.arp11.topmcl3
  open(OUT,">$outname") or die "$outname";
  print OUT join("\t",qw(TargGeneid SrcGeneid Slen Ogroup Bit Idn Aln Tlen Olen dTO xTO))."\n";  
  my (%did);

  sub puta { 
    my($ar,$tg,$sg)=@_;  
    my($tl,$dl,$xl); 
    my($n,$md,$mn,$sd,$min,$max);
    my($olen,$oid,@v)=split"\t",$togtab{$tg}{$sg}; 
    # $tl= @v[-1];  # WRONG now.. Slen not Tlen: Bits Ident Align  Qlen Slen
    $tl= $v[3]; # T/Qlen; use  $v[4] == olen/Slen ?

   ## arbm: change here bits > not highest score but median of orthogroup?
    #   * ave Bits = highest bitscore match, switch to average( median bits to orthogroup species)?
    #     so close relatives (2 daphnia) wont skew score.
    
    if($USE_MEDBITS) {
      my($mdaa,$naa,$minaa,$maxaa)= median( $arbm{$ar}{$tg} ); # screwy results with eg 2 entries
      $v[$SCORECOL]= $mdaa if($naa>=$MINTAXA); # was naa>0, require MINTAXA 
    }
    
   ## redo this size-outlier calc: use md, median size, and +/- 2stdev outlier?
    if(1) {   
        # ($md,$n,$min,$max)= median( $arlm{$ar} );  
        my @sl=sort{$b<=>$a} @{$arlm{$ar}};
        ## update for MINTAXA
        $n=@sl; 
        if($n < $MINTAXA) { $md=$sl[0]; $xl=0; }
        else { $md= @sl[int($n/2)];  $xl=($tl > 1.9*$md)?2:($tl < 0.65*$md)?-2:0; } # maybe this, 1.5 or .66 outside median
    } else {   
        # 3sd outlier may be less meaningful for this non-normal, small sample data than above
        ($n,$md,$mn,$sd,$min,$max)= medmeansd( $arlm{$ar} );
        my $sd2= 3 * $sd;  # 2 or 3?
        $xl=($tl > $md+$sd2)?2:($tl < $md-$sd2)?-2:0; # this gives higher counts than above 1.5 or .66 cuts
        # $xl=($tl > $mn+$sd2)?2:($tl < $mn-$sd2)?-2:0; # mean? no
    }

   $dl=$tl-$md;
   my $poid= ($oid =~ /^\d/) ? $idprefix.$oid : $oid;
   print OUT join("\t",$tg,$sg,$olen,$poid,@v,$md,$dl,$xl)."\n"; 
   $did{$tg}++; $did{$ar}++; 
   } # sub puta

  my @oids= sort{$a <=> $b or $a cmp $b} keys %bax; 
  foreach my $ar (@oids) {  
    my($bt,$tg,$sg); $bt=$bax{$ar}; $tg=$baxt{$ar};  $sg= $baxs{$ar};  
    puta($ar,$tg,$sg) unless($did{$tg} or $btx{$tg}>2+$bt); 
    
    # FIXME: -alltarg, put all tg, using keys $arbm{$ar}{$tg} ?? or keys  $btx{$tg}
    if($ALLTARG) {
      my @tga= sort{ $btx{$b} <=> $btx{$a} } grep{ $btx{$_} and $_ ne $tg } keys %{$arbm{$ar}};
      for my $tga (@tga) {
        my $ara= $btxa{$tga} || "nada"; # should be same as ar ?
        my $sga= $btxs{$tga} || "nosid";
        puta($ara,$tga,$sga) unless($did{$tga});
      }
     }
    } 
   
  foreach my $ar (@oids) { next if($did{$ar}); 
    my($bt,$tg,$sg); $bt=$bax{$ar}; $tg=$baxt{$ar};  $sg= $baxs{$ar}; 
    puta($ar,$tg,$sg) unless($did{$tg}); 
    } 

  close(OUT);
  return($outname);
}  


sub putarow {  # fixme above mess .. not used..
  my($oid, $tg, $sg, $tgsize, $sgsize, $mdsize, $values)= @_;  
  # my @vfake= (199, 99, 99); # Bits Iden Algn == $values
  #?? add 1, -1 levels
  my $dl= $tgsize - $mdsize;
  my $xl=($tgsize > $AA_TOOBIG*$mdsize) ? 2 : ($tgsize < $AA_TOOSMALL*$mdsize) ? -2 : 0;  
  if($AA_OUTLIER1 and $xl == 0) { #? move this to topstats as output filter?
    $xl= ($tgsize > $AA_TOOBIG1*$mdsize) ? 1 : ($tgsize < $AA_TOOSMALL1*$mdsize) ? -1: 0;  
  }
  print OUT join("\t",$tg,$sg,$sgsize,$oid, @$values, $tgsize ,$mdsize,$dl,$xl)."\n"; 
} 

=item toptab

  TargGeneid      SrcGeneid       Slen    Ogroup  Bit     Idn     Aln     Tlen    Olen    dTO     xTO
  killifish:Funhe5EG030598t1      zebrafish:E7FDA5_DANRE  743     FISH0   221     225     518     705     289     416     2
  killifish:Funhe5EG015351t1      medaka:ENSORLP00000008122       267     FISH10  162     185     281     434     267     167     0

=cut

sub topstats
{
  my($toptab)= @_;
  # let toptab == stdin
  
  my( %arn, %arng, %arx, %ard, %arb, %arpa, %ts, %tsx);
  my $SCORELAB= qw(Bits Ident Algn)[$SCORECOL];
  
  my $inh={};
  if($toptab =~ /^-|^stdin/) { $inh=*STDIN; }
  else { open(IN,$toptab) or die "$toptab"; $inh= *IN; }
  while(<$inh>) {
    
    next if(/^TargGeneid/); 
    my($tg,$sg,$sgsize,$aoid,@vals)=split; 
    # ts tag from geneid: need help:  acyp2ref, acyp2eg << not same...
    my $ts="miss"; 
    unless( ($ts)= $tg=~m/^(\w+):/) {
      unless( ($ts)= $tg=~m/^([^\d\W]+\d[^\d\W]+)/ ) { ($ts)= $tg=~m/^([^\d\W]+)/; }
    }
    
    my($mdsize,$dl,$xl)= @vals[-3,-2,-1];
    ##x $xd=$v[-1]; $d=$v[-2]; # note end of vec address, inner @v has extra cols
    my $bscore= $vals[$SCORECOL];
    my $paln  = 100*$vals[2]/$sgsize;
    
## input row now of aaset3arp7-apimel14nc.omcl.topalgn has extra end col from tall4.
# now row: tg,sg,slen,ar, @v=[Bit,Idn,Aln,Tlen,Slen,medOlen,xx1,xx2, xx3, dTO, xTO ]
#   $tg,$sg,$sgsize,$oid,@values,$mdsize,$dl,$xl
# add other stats?  paln = $v[2]/slen or $v[2]/
# TargGeneid	SrcGeneid	Slen	Ogroup	Bit	Idn	Aln	Tlen	Olen	dTO	xTO
# apimel14nc:XP_006560612	human:UniRef50_P58005	492	ARP7f_G1001	[v 502	244	426	948	492]	[vx 491	210	1,1] [vo 480	468	2]
# apimel14nc:XP_001121381	nasvit:Nasvi2EG011337t1	315	ARP7f_G1003	263	142	252	323	315	325	158	1,1	315	8	0
 
    # >> fix here for many gene/spp/group : +=$dl +=$bscore ?? xd
    $arb{$aoid}{$ts} += $bscore; # blast score (bits,idn,aln)
    $arpa{$aoid}{$ts} += $paln; #  
    $ard{$aoid}{$ts} += $dl; # diff from med-aasize
    #move above: $xl=0 if( !$AA_OUTLIER1 and ($xl==1 or $xl==-1));
    #old? $arx{$aoid}{$ts}=$xl if($xl); #? change to 3hash:  $arx{$aoid}{$ts}{$xl}++
    $arx{$aoid}{$ts}{$xl}++ if($xl); # aasize outlier score
    $arn{$aoid}++; $arng{$aoid}{$ts}++; $ts{$ts}++;   $tsx{$ts}{$xl}++; 

  } close($inh);

  my @ts=sort keys %ts; my $NG=@ts; my @ar=sort keys %arb; 
  my @ap=grep{ scalar(keys %{$arb{$_}})==$NG } @ar; 
  print "Summary stats for $toptab\n" unless($toptab=~/stdin/);  
  
   #printf "%-9s: ","Geneset"; print join("\t",qw(Ng Ngr Bits dSize rBits oddSize)),"\n";
   #print "n=$nt; nreal=$nr; abits=$ab; dlen=$ad; arealbits=$abr; xtralen:@xl\n";
  my %xlkey= ( 2=>"x2Big", -2 =>"x2Tiny", 1=>"x1big", -1 =>"x1tiny" );
  
  my @cset= ($DOCOM and @ts>1) ? (0,1) : (1); # Common not for 1 species tables
  for my $k (@cset) {
    my @aset= ($k==0) ? @ap : @ar;
    print "# ...  ",(($k==0)? "Common" : "All"), " groups ... \n";
    printf "%-${NAMWID}s ","Geneset"; 
    
    # print join("\t",qw(nGene Bits dSize rGene rBits outlierSize)),"\n";
    # change nGene/rGene > nGroup/rGroup  add nGene for DOALL
    my @mlab= qw( pAln dSize nGene);
    push @mlab, "pGene" if($DOREFPCT);
    my @xl= ($DOTABTINY) ? qw(nTiny	pTiny	nBig pBig) : qw(outlierSize); 
    if($DOALL) {
    print join("\t","nGroup", $SCORELAB, @mlab, "rGroup","r".$SCORELAB,@xl),"\n";    
    ## was qw( pAln dSize rGroup nGene)
    } else {
    print join("\t","nGroup", $SCORELAB, @mlab,"r".$SCORELAB,@xl),"\n";    
    }
   
    my (%sb,%sn,%sal,%sx,%xdk,%sd,%snt,%sng);    
    foreach my $aoid (@aset) { foreach my $t (@ts) {  
      my($bscore,$paln,$dl,$ng);
      $bscore=$arb{$aoid}{$t}; $dl=$ard{$aoid}{$t}; 
      $paln= $arpa{$aoid}{$t};
      #old# my $x=$arx{$aoid}{$t}; $sx{$t}{$x}++ if($x); # ok for 2 outlier levels here?  -2,+2 and -1,+1
      my @x= keys %{$arx{$aoid}{$t}}; map{ $sx{$t}{$_}++; $xdk{$_}++; } @x; # new
      # FIXME: DOALL, need:  sng{$t} += $ng;
      $ng= $arng{$aoid}{$t}; $sng{$t} += $ng;
      $sb{$t}+=$bscore; $sal{$t}+=$paln;
      $sn{$t}++; $sd{$t}+=$dl; $snt{$t}++ if($bscore); 
      } 
    } 

    my $decp=($digits>0)?10:1;
    my @xdk= sort keys %xdk;
    foreach my $t (@ts) { my($nt,$ab,$apa,$nr,$pref,$abr,$ad,@xl,$ng);
      $nt=$sn{$t}||1; $nr=$snt{$t}||1; $ng= $sng{$t}||1;
      $ab=int(0.5+$decp*$sb{$t}/$nt)/$decp;  
      $apa=int(0.5+$decp*$sal{$t}/$nt)/$decp;  
      $abr=int(0.5+$decp*$sb{$t}/$ng)/$decp;  
      $ad=int(0.5+$decp*$sd{$t}/$ng)/$decp; # ng == nr unless DOALL : should be
      $pref= int(0.5+$decp*100*$nr/$nt)/$decp;  
      ## fixme missing x2Tiny .. set to zero
      if($DOTABTINY) { # full table tab, not xx=nnn (pct)
        @xl= map{ my $xc=$sx{$t}{$_}||0; my $xp=int(1000*$xc/$nr)/10; ($xc,$xp); } sort{$a<=>$b} @xdk;       
      } else {
        @xl= map{ my $xc=$sx{$t}{$_}||0; my $xp=int(1000*$xc/$nr)/10; $xlkey{$_}."=$xc ($xp%)" } 
          sort{$a<=>$b} @xdk; #was misst: keys %{$sx{$t}};      
      }  
      my @mval= ($apa,$ad,$nr); # qw( pAln dSize nGene);
      push @mval, $pref if($DOREFPCT);
      printf "%-${NAMWID}s ",$t; 
      # print "n=$nt; nreal=$nr; abits=$ab; dlen=$ad; arealbits=$abr; xtralen:@xl\n";  
      if($DOALL) {
      print join("\t",$nt,$ab,@mval,$ng,$abr,@xl),"\n"; # add ng
      } else {
      print join("\t",$nt,$ab,@mval,$abr,@xl),"\n"; 
      }
    }
  }
  
}





sub orthogroupgenes
{
  my($groupcounts,$geneaa,$genegroups)= @_;
  
  my(%ogroup,%aasize,%god,%speciesaa);
  warn "#i orthogroups from $groupcounts\n" if($DEBUG);

  open(CIN,"$groupcounts") or die "$groupcounts";
  ## add idprefix.od or not ??
  while(<CIN>) { my($od,$nt,$ng,@c)=split; $ogroup{$od}++ if($nt>=$MINTAXA); } close(CIN);
  warn "#i group n=",scalar(keys %ogroup)," with ntaxa>=$MINTAXA\n" if($DEBUG);
  
  my $iscount=0;
  if($geneaa =~ /count|\.qual/) { # bad but ok..
    open(AASIZE,$geneaa) or die "aa count in $geneaa ..."; $iscount=1;
  } else {
    open(AASIZE,"faCount $geneaa | cut -f1,2 |") or die "faCount  $geneaa ...";
  }
  while(<AASIZE>) {
    my($id,$al)=split; 
    next unless($id and $al>0); # if($iscount and (not defined($al) or $al=~/\D/)) {}
    $id= fixids($id); # $id=~s/:/_/; 
    my($sp)=split ":",$id; 
    # $id=~s/p(\d)/t$1/ if(/wasp/); # damn this fixup crap
    ##if(/aphid:gi/) { ($rd)= m/ref\|([\w.]+)/; $id=~s/gi\|(\w+).*/gi:$1/; $id2{$id}="aphid_".$rd if($rd); } 
    $aasize{$id}=$al;
    $speciesaa{$sp}++;
  } close(AASIZE);
  my $sppcount= join", ", map{ "$_:$speciesaa{$_}" } sort keys %speciesaa;
  warn "#i geneaa n=",scalar(keys %aasize)," sppcount $sppcount\n" if($DEBUG);
 
  open(OGENES,"$genegroups") or die "$genegroups"; # arp11u11_omclgn.tab
  while(<OGENES>) {
    my($od,$gn)=split; 
    next if($keepspecies and not($gn =~ m/$keepspecies/)); # both keep,skip?
    next if($skipspecies and $gn =~ m/$skipspecies/);
    #?? next unless(  $speciesaa{$sp} );
    $gn= fixids($gn); ## DANG: $warp/arp11u11_omclgn.tab has wasp_000t1, aasize, tall4 has wasp_000p1
    my($odnum)= $od =~ m/(\d+)$/;  # was BAD patt for od= FISH11d123 : (\d+)
    $god{$gn}=$od if($ogroup{$odnum}); # save od or odnum ?
    # print OGENETAB join("\t",$od,$gn,$al)."\n" if($ogroup{$od}); # internal tab????
  } # > arp11u11_omcl11.gr3.gids
  close(OGENES); ##close(OGENETAB);
  warn "#i groupgenes n=",scalar(keys %god),"\n" if($DEBUG);

  return(\%god, \%aasize);
}



=item blast2tab ==  makeblastscore

 env aa=$geneaa tall=1 $evigene/scripts/makeblastscore.pl $geneblastp > ogs12.arp11.tall

 perl -pi -e 'if(/\taphid:gi/) { s/aphid:gi.\w+.ref./aphid:/; s/\|//g; } ' ogs12.arp11.tall

=cut

# sub blast2tab
# {
#   my($geneblastp,$geneaa)= @_;
#
#   ## instead of this, fix makeblastscore.pl to add aasize ...
#   open(AASIZE,"faCount $geneaa | cut -f1,2 |") or die "faCount  $geneaa ...";
#   while(<AASIZE>) {
#     my($id,$al)=split; $id = fixids($id); 
#     $aasize{$id}=$al;
#   } close(AASIZE);
# 
#   open(BLTAB,"env tall=1 $evigene/scripts/makeblastscore.pl $geneblastp |") or die "makeblastscore.pl $geneblastp";
#   while(<BLTAB>) {
#      ($id,$og,@v)=split; $al=$aasize{$id}||0; $v[3]=$al unless(/^Query/); print OUT join("\t",$id,$og,@v)."\n";
#   }
#}



=item scripts to convert
  
  evorthogenes.pl  -tallscore evg11u.arp11.tall4 -groupgene arp11u11_omclgn.tab -groupcount arp11u11-orthomcl-count.tab -geneaa arp11u11noalt.aa.gz 
  
  env tall=1 $evigene/scripts/makeblastscore.pl bp1-arp11u11noalt-ogs12.aa.blastp.gz > ogs12.arp11.tall
  perl -pi -e 'if(/\taphid:gi/) { s/aphid:gi.\w+.ref./aphid:/; s/\|//g; } ' ogs12.arp11.tall
  
  cat arp11u11-orthomcl-count.tab | env nt=3 gtag=ARP perl -ne\
  '($od,$nt,$ng,@c)=split; print "$ENV{gtag}$od\t\n" if($nt>=$ENV{nt});' \
    > arp11u11_omcl11.gr3.oids
  
  faCount ../nasvit2genes11s_ball/arp11u11noalt.aa.gz | cut -f1,2 | sed 's/^/len /' |\
  cat - arp11u11_omcl11.gr3.oids arp11u11_omclgn.tab | perl -ne\
  'if(s/^len //){ ($id,$al)=split; $id=~s/:/_/; $id=~s/p(\d)/t$1/ if(/wasp/); \
  if(/aphid:gi/) { ($rd)= m/ref\|([\w.]+)/; $id=~s/gi\|(\w+).*/gi:$1/; $id2{$id}="aphid_".$rd if($rd); }  $al{$id}=$al; }\
  elsif(/^(\w+)\s*$/) { $ok{$1}=1 } else { ($od,$gn)=split; \
  $al=$al{$gn}||0; $gn=$id2{$gn}||$gn; print join("\t",$od,$gn,$al)."\n" if($ok{$od}); }' \
  > arp11u11_omcl11.gr3.gids
  
  
  faCount pubdir/genes/nvit2_evigenes_pub11u.aa.gz | cut -f1,2 | \
  sed 's/^/len /' | cat -  genes/pub11u/aaeval/evg11u.arp11.tall |\
  perl -ne 'if(s/^len //){ ($id,$al)=split;  $id=~s/t(\d+)$/p$1/; $al{$id}=$al if($al>0); }\
  else { ($id,$og,@v)=split; $al=$al{$id}||0; $v[3]=$al unless(/^Query/); print join("\t",$id,$og,@v)."\n"; }'\
  > genes/pub11u/aaeval/evg11u.arp11.tall4
  
  grep -v wasp prot/omcl9u11/arp11u11_omcl11.gr3.gids | cat - genes/pub11u/aaeval/evg11u.arp11.tall4 |\
  env gtag=ARP perl -ne 'if(/^$ENV{gtag}/){ ($od,$gn,$al)=split; $gn=~s/_/:/; $god{$gn}=$od;  $glen{$gn}=$al; }\
  elsif(/^\w/ and /:/) { ($tg,$gn,@v)=split"\t"; $gn=~s/\|$// if($gn=~/aphid/); \
  $od=$god{$gn}||"na"; $l=$glen{$gn}||0; print join("\t",$tg,$gn,$l,$od,@v); }'\
  > genes/pub11u/aaeval/evg11u.arp11.tomcl
  
  # fixme, need arls or arlmed : omcl group median/mean length for best diff-length
  # fixme2: alttr problem?; no, did remove by did{$ar} ; topmcl2 n=751 are alts not p1
  
  cat genes/pub11u/aaeval/evg11u.arp11.tomcl  | perl -ne \
  '{ ($tg,$sg,$sl,$ar,$b,$i,$a,$tl)=split; \
  $arn{$ar}++; $arls{$ar}+=$sl; push @{$arlm{$ar}}, $sl; $dl=$tl-$sl;  \
  $bav{$tg}{$sg}=join"\t",$sl,$ar,$b,$i,$a,$tl; $bar{$ar}{$tg}=$b;  \
  if($b>$bax{$ar}) { $bax{$ar}=$b; $baxt{$ar}=$tg; $baxs{$ar}=$sg; }  \
  if($b>$btx{$tg}) { $btx{$tg}=$b; $btxa{$tg}=$ar; $btxs{$tg}=$sg; } } \
  END{  \
  foreach $ar (sort keys %bax) { \
   $b=$bax{$ar}; $tg=$baxt{$ar};  $sg= $baxs{$ar};  \
   puta($ar,$tg,$sg) unless($did{$tg} or $btx{$tg}>2+$b); } \
  foreach $ar (sort keys %bax) { next if($did{$ar}); \
   $b=$bax{$ar}; $tg=$baxt{$ar};  $sg= $baxs{$ar}; \
   puta($ar,$tg,$sg) unless($did{$tg}); } } \
  BEGIN{ print join("\t",qw(TargGeneid SrcGeneid Slen Ogroup Bit Idn Aln Tlen Olen dTO xTO))."\n"; }\
  sub puta { my($ar,$tg,$sg)=@_;  my($ml,$tl,$dl,$xl); \
   my @v=split"\t",$bav{$tg}{$sg}; my @sl=sort{$b<=>$a} @{$arlm{$ar}}; \
   $ml= @sl[int(@sl/2)]; $tl= @v[-1]; $dl=$tl-$ml; \
   $xl=($tl > 1.5*$sl[0])?2:($tl < 0.65*$sl[-1])?-2:0; \
   print join("\t",$tg,$sg,@v,$ml,$dl,$xl)."\n"; \
   $did{$tg}++; $did{$ar}++; } '\
  > genes/pub11u/aaeval/evg11u.arp11.topmcl3
  
  ## summary stats
  
  cat genes/pub11u/aaeval/*.arp11.topmcl3 | perl -ne\
  'next if(/^TargGeneid/); ($tg,$sg,$sl,$ar,@v)=split; ($ts)=$tg=~m/^([^\d\W]+)/; \
  $b=$v[0]; $xd=$v[-1]; $d=$v[-2]; $arb{$ar}{$ts}=$b; $ard{$ar}{$ts}=$d; \
  $arn{$ar}++; $ts{$ts}++;  $arx{$ar}{$ts}=$xd if($xd); $tsx{$ts}{$xd}++; \
  END{ @ts=sort keys %ts; $NG=@ts; @ar=sort keys %arb; \
  @ap=grep{ scalar(keys %{$arb{$_}})==$NG } @ar; foreach $ar (@ap) { foreach $t (@ts) { \
  $b=$arb{$ar}{$t}; $d=$ard{$ar}{$t}; $x=$arx{$ar}{$t}; $sx{$t}{$x}++ if($x); \
  $sb{$t}+=$b; $sn{$t}++; $sd{$t}+=$d; $snt{$t}++ if($b); } } \
  foreach $t (@ts) { $nt=$sn{$t}; $ab=int(10*$sb{$t}/$nt)/10;  \
  $nr=$snt{$t}; $abr=int(10*$sb{$t}/$nr)/10;  $ad=int(10*$sd{$t}/$nr)/10; \
  @xl=map{ $xc=$sx{$t}{$_}; $xp=int(1000*$xc/$nr)/10; "xl$_=$sx{$t}{$_}($xp%); " } sort keys %{$sx{$t}}; \
  print "$t: n=$nt; nreal=$nr; abits=$ab; dlen=$ad; arealbits=$abr; xtralen:@xl\n"; } }'

=cut

=item prelim summary stats

  # all groups
  NasviEg: n=9159; nreal=8940; abits=658.1; dlen=-11.9; arealbits=674.2; xtralen:xl-2=269(3.0%);  xl2=104(1.1%); 
  NcbiRef: n=9159; nreal=8193; abits=629.7; dlen= -3.3; arealbits=703.9; xtralen:xl-2=179(2.1%);  xl2= 63(0.7%); 
  NVogs12: n=9159; nreal=7696; abits=589.3; dlen=-11.7; arealbits=701.3; xtralen:xl-2=242(3.1%);  xl2=170(2.2%); 
  
  # common groups :
  NasviEg: n=7077; nreal=7077; abits=750.0; dlen=-8.4; arealbits=750.0; xtralen:xl-2=109(1.5%);  xl2=29(0.4%);   
  NcbiRef: n=7077; nreal=7077; abits=749.8; dlen=-4.0; arealbits=749.8; xtralen:xl-2=117(1.6%);  xl2=32(0.4%); 
  NVogs12: n=7077; nreal=7077; abits=728.7; dlen=-9.7; arealbits=728.7; xtralen:xl-2=173(2.4%);  xl2=122(1.7%); 
      dlen = ave. aa size difference from group median
      xl-2 = ngene aalen < 66% ogroup min aa size; xl2 = ngene aalen > 150% ogroup max aa size

  

=cut

=item data

  ==> ../pub11u/aaeval/ncbiref.arp11.tall <==
  Query   Source  Bits    Ident   Align   Qlen
  NcbiRef2rna6326 bombusterr:XP_003394015.1       509     473     693     
  NcbiRef2rna6326 apis2ref:NP_001107660.1 496     534     725     
  NcbiRef2rna6326 drosmel:NP_524042.2     328     237     339     
  
  prot/omcl9u11/arp11u11_omclgn.tab
  ARP0    aphid_ACYPI000809-PA
  ARP0    aphid_ACYPI001952-PA
  ARP0    aphid_ACYPI002992-PA

=cut
