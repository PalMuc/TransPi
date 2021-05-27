#!/usr/bin/env perl
# asmrna_equalgenes.pl
# gene2asmrna.pl

=item about
  
  create ref transcript assembly from reference genes/genome and superset of redundant tr assemblies.

=item inputs
  
  refgenes.over*asmrna.tab  : genome location equivalence, scored with percent cds-over.exon-over
    (evigene/equalgene.pl)
  
  cd-hits/refgenes.cdhit*asmrna{100,99,98,97,96,95,91}.aa.clstr
    : aa-equvalence as cd-hit clustering (evigene/cdallsame.info,pl)
    : should include refgenes equivalences (e.g. alt-trs w/ same aa)
    
   option: order preference for asmrna subsets (by source, presumably), where
      duplicate/equivalent asm among sources should be removed.

=item eg

  $evigene/scripts/rnaseq/asmrna_equalgenes.pl -geneid Thecc \
    -ingeneover pub3ig.allover.tab1 -inaaclust cdall100/pub3i_trall9_all.clstr > pub3ig.trasm.match.tab4

  * add -trokids trok.ids to promote these above others
  
  #NOTE: genes_overtrasm: ngene=44404, ntr=264389, ntrskip=1123504
  #NOTE: aacluster: nclust=33167, ngene=33167, ntr=117057, neqgene=14775, ndup=0
  #NOTE: gene_trasm_equiv: nkeep=205881, nskip=62497

Thecc1EG000001t1        na      nemap:0 aamiss:0,0,0
Thecc1EG000002t1        cacao3nwbL_g13025t00001 full:100.99     aaeq:205,100,100
Thecc1EG000002t1        cacao3cuf8_Gsc1g253t1   full:99.92      aaeq:205,100,100
Thecc1EG000002t2        cacao3nwbB_g11013t00001 full:96.98      aaeq:199,100,97
Thecc1EG000002t2        cacao5tri1sub3sc1loc106c0t1     full:96.95      aaeq:199,100,97
Thecc1EG000002t2        cacao3nwbP2_g03816t00001        full:96.93      aaeq:199,100,97 < sort above cacao5tri
Thecc1EG000002t2        cacao5sopcsc1k89loc5659t1       full:96.90      aaeq:199,99.50,97
Thecc1EG000003t1        cacao3cuf8_Gsc1g6t5     full:99.79      aaeq:281,100,100
Thecc1EG000003t2        cacao3cuf8_Gsc1g6t4     full:99.92      aaeq:285,100,100
Thecc1EG000003t2        cacao3vel5sc1Loc2319t5  full:99.92      aaeq:285,100,100

Thecc1EG000014t1        cacao3cuf8_Gsc1g63t2    full:99.92      aamiss:0,0,0


=item output(s)

  table of refgenes + asmrna [best] equivalences from both location and aa-match,
    one row per mRNAid x asmrnaid
    with equal scores indicating degree of equivalence:
      identical = same CDS locus, same aa;
      same-CDS, diff-aa; 
      diff-CDS, same-aa;
      partial = some CDS, some aa agreement
      and CDS/UTR score? utrpoor, utrbad for asmrna w/ same CDS/aa but likely extra gene parts in utr
      
    refgenes_overasm : output refgenes.overall.tab as before?
       
=item author
  
  don gilbert, gilbertd near indiana edu, 2012
  part of EvidentialGene, evigene/scripts/

=cut

use strict;
use warnings;
use Getopt::Long;

my $DIFFMAX_OVERASM = 5.0;  # percent difference from best-overasm to worst to keep
my $MINGENEIDENT= 95; #? 85; # for gene1aa == gene2aa identity ? other?
my $FILTER_OVRNA= 0.90;
my $FILTER_AARNA= 0.90;

my $MINFULLCDS= 90;
my $MINFULLEXON= 60;
my $MINPID = 95;
my $MINPDIF= 95;

# my $pMINLEN = $ENV{minlen} || 0.95;
# my $pMINID  = $ENV{minid} || 0.95;

my $IDPREFIX= $ENV{idprefix} || ""; # refgene id pattern; need for aaclusters at least,    

my $debug= 1;
my ($geneinfo,$ingeneover,$outgeneover,$inaaclusters,$outaaclust,$outeqtab,$trokids,$trpoorids,$logfile);

my $optok= GetOptions(
  "ingeneover=s", \$ingeneover,  #?? maybe file list
  "outgeneover=s", \$outgeneover,  # multiple outfiles? change suffix?
  "inaaclusters=s", \$inaaclusters,  #?? maybe file list
  "outaaclust=s", \$outaaclust,   
  "outeqtab=s", \$outeqtab,   
  "trokids|trgoodids=s", \$trokids,   
  "trpoorids|trbadids=s", \$trpoorids,  
  "logfile=s", \$logfile,  
  "geneidprefix=s", \$IDPREFIX,  
#  "geneinfo=s", \$geneinfo,  
#   "MINSIZE=i", \$MINSIZE,  
#   "MAXGAP=i", \$MAXGAP,  
  "MINGENEIDENT=i", \$MINGENEIDENT,  
  "DIFFMAX_OVERASM=i", \$DIFFMAX_OVERASM,  
  "debug!", \$debug, 
  );

die "usage: asmrna_equalgenes.pl  -geneidprefix Thecc -ingeneover xxx -inaaclusters yyy ...
 opts:
" unless($optok and $ingeneover and $inaaclusters);
#  -idprefix Thecc1RA_ -vecscreen=infile -geneinfo=infile -MINSIZE=$MINSIZE  -MAXGAP=$MAXGAP  -out=outfasta  -log=outlog


my $MINPID1 = $MINPID*0.90;
my $MINPDIF1= $MINPDIF*0.90;
#.......



my ($inh_geneoverasm, $outh_geneoverallasm);
  # ingeneover#?? maybe file list
  if($ingeneover =~ /stdin|^-/) { $inh_geneoverasm=*STDIN; }
  else { open(INO,$ingeneover)  or die "reading $ingeneover"; $inh_geneoverasm=*INO; }
  if($outgeneover) { open(OUTO,">$outgeneover")  or die "writing $outgeneover"; $outh_geneoverallasm=*OUTO; }

my $genesoverall= genes_overtrasm($inh_geneoverasm, $outh_geneoverallasm);
  close($inh_geneoverasm) if($ingeneover);
  close($outh_geneoverallasm) if($outgeneover);

my ($inh_aaclusters, $outh_aaclust);
  # inaaclusters#?? maybe file list
  if($inaaclusters =~ /stdin|^-/) { $inh_aaclusters=*STDIN; }
  else { open(INAA,$inaaclusters)  or die "reading $inaaclusters"; $inh_aaclusters=*INAA; }
  if($outaaclust) { open(OUTAA,">$outaaclust")  or die "writing $outaaclust"; $outh_aaclust=*OUTAA; }

my $aaclusters= aaclusters($inh_aaclusters, $outh_aaclust);
  close($inh_aaclusters) if($inaaclusters);
  close($outh_aaclust) if($outaaclust);

my ($outh_eqtab);
  if($outeqtab) { open(OUTET,">$outeqtab")  or die "writing $outeqtab"; $outh_eqtab=*OUTET; }
tabulate_gene_trasm_equiv($genesoverall, $aaclusters, $outh_eqtab);  
  close($outh_eqtab) if($outeqtab);

if($outeqtab) { 
  my(%oktrids,%poortrids);
  if($trokids){ open(F,$trokids) or die "reading $trokids"; while(<F>) { chomp; $oktrids{$_}=1 if(/^\w/); } close(F); }
  if($trpoorids){ open(F,$trpoorids) or die "reading $trpoorids"; while(<F>) { chomp; $poortrids{$_}=1 if(/^\w/); } close(F); }
    # map{ $oktrids{$_}=1 } grep /\w/, split;
  my ($inh_eqtab,$outh_besttab);
  open(INET,$outeqtab)  or die "reading $outeqtab"; $inh_eqtab=*INET;
  open(OUTBT,">$outeqtab.best1")  or die "writing $outeqtab.best1"; $outh_besttab=*OUTBT; 
  my $okh= ($trokids) ? \%oktrids : undef;
  my $poorh= ($trpoorids) ? \%poortrids : undef;
  pick_gene_trasm_best($inh_eqtab,$okh,$poorh,$outh_besttab);
}

#-----------------------------------

=item old genetrmatch.tab

  head subtsa2.cdmatch.allover.tab3b
GeneID          flag  asmsource   asmid               ovqual:overid/ovscore      aaqual:aasize,pid,pdif
Thecc1EG000002t1	TSA	Newbler:10	cgbaL_g13025t00001	full:L_g13025t00001/100.99	aaeq:205,100,100
Thecc1EG000002t2	TSA	Newbler:10	cgbaL_g13025t00001	full:B_g11013t00001/96.98	naaeq:205,100,97,B_g11013t00001
Thecc1EG000003t1	TSA	Cufflinks:08	cacao11r39cuf8_Gsc1g6t5	full:caca11r39cuf8_Gsc1g6t5/99.79	aaeq:281,100,100
Thecc1EG000003t2	TSA	Cufflinks:08	cacao11r39cuf8_Gsc1g6t5	full:cacao3vel5sc1Loc2319t5/99.92	naaeq:285,100,100,cacao3vel5sc1Loc2319t5
Thecc1EG000005t1	TSA	Cufflinks:08	cacao11r39cuf8_Gsc1g128t2	full:caca11r39cuf8_Gsc1g128t2/99.89	aaeq:1269,100,100
Thecc1EG000005t2	TSA	Velvet:11	cacao3vel4sc1Loc938t2	full:cacao3vel4sc1Loc938t2/100.99	aaeq:1156,100,100

=item new genetrmatch.tab

GeneID            asmid             ovqual:ovscore      aaqual:aasize,pid,pdif
Thecc1EG000002t1	L_g13025t00001	  full:100.99	  aaeq:205,100,100
Thecc1EG000005t2	cacao3vel4sc1Loc938t2	full:100.99	aaeq:1156,100,100

=cut

sub tabulate_gene_trasm_equiv
{
  my($geneover,$aaclust,$outh,)= @_; 
  my($nkeep,$nall)=(0,0);
  unless(defined $outh) { $outh= *STDOUT; }
  
  # geneover: $overall{$lgd}= $ovp; : ovp = list of trid/overcds.overexon,.. sorted best>worst
  # aaclust:  $aaclust{$refid}{$tid}= join"\t",$pboth,$pid,$pdif,$rsize,$tsize;
  my %didid; # also filter dup asmrna1 x gene1,2 matches
             # .. need score best match first, of rna1 x gene1 or gene2
             
  my @gid= sort keys %$geneover;
  foreach my $gid (@gid) {
    my $ovrna = $geneover->{$gid};
    unless(defined $aaclust->{$gid}) { $aaclust->{$gid}={}; }
    
#    # FIXME: for now, ignore any aarna ids not in ovrna list .. should report as nonover/other
#     my @aarna= ();
#     if(ref $aaclust) {
#       #? need sorted aaclust? 
#       @aarna= sort{ $aaclust->{$b} <=> $aaclust->{$a} } keys %{$aaclust};
#     }
    
    #  _bestrnasort internal to tabulate_gene_trasm_equiv : problems shared vars
    our $gidS= $gid;
    sub _bestrnasort {
      our $gidS;
      my($ta,$ca)= split"/",$a; # $a=~m,(.+)/([\d\.]+),;
      my($tb,$cb)= split"/",$b; # $b=~m,(.+)/([\d\.]+),; 
      $ca =~ s/\.(\d+)//; my($xa)= $1;
      $cb =~ s/\.(\d+)//; my($xb)= $1;
      return $cb <=> $ca 
        or $aaclusters->{$gidS}{$tb} <=> $aaclusters->{$gidS}{$ta}
        or $xb <=> $xa
        or $a cmp $b;  #? add order by trsource == tidprefix
    }
    
    my @ovrna = split",", $ovrna;  # re-sort, using also aaclust score?
    @ovrna= sort _bestrnasort @ovrna;
    
    my($lov, $lpid, $lpdif, $lboth, $nov)=(0) x 9;    
    foreach my $trov (@ovrna) {
      my($tid,$tov)= split"/",$trov;  $tov||= 0;
      $nall++;
      # already done# $tov =~ s,/I100,/100.99,; $tov =~ s,/C,/,; 
      # already done: next if($tov < $maxov-$DIFFMAX_OVERASM);
      
      my $aascore= $aaclust->{$gid}{$tid} || ""; # FIXME: any trid fixes here?
      my ($pboth,$pid,$pdif,$rsize,$tsize)= (0) x 9;
      ($pboth,$pid,$pdif,$rsize,$tsize)= split"\t",$aascore if($aascore);

      next if($FILTER_OVRNA and $tov < $lov * $FILTER_OVRNA);
      next if($FILTER_AARNA and $pboth < $lboth * $FILTER_AARNA);
  
      my ($ce,$xe)= split /\./,$tov;
      my $ok=($ce>$MINFULLCDS and $xe>$MINFULLEXON)?2:1; # only 1,2 here; no 0/part
      my $ovflag= qw(part nemap full)[$ok]; 

      my $aaflag=($aascore eq "")?"aamiss" 
          : ($pid>=$MINPID && $pdif>$MINPDIF)?"aaeq" # ($pid>99 and $pdif>98) 
          : ($pid>=$MINPID1 && $pdif>$MINPDIF1)?"aasim":"aapoor";
     
       ## filter out lower qual ovrna/aascore
      $lov= $tov if($lov < $tov);
      ($lpid, $lpdif, $lboth)= ($pid, $pdif, $pboth) if($lboth < $pboth);
                  
      print $outh join("\t",$gid,$tid,"$ovflag:$tov","$aaflag:$rsize,$pid,$pdif")."\n"; 
      $nkeep++;
    }
  }
  
  my $nskip= $nall - $nkeep;
  warn "#NOTE: gene_trasm_equiv: nkeep=$nkeep, nskip=$nskip\n" if $debug;
  return($nkeep,$nskip);
}

=item option to filter to best asm/mrna/gene

  * also option input trokset.ids, or trpoorset.ids, to keep/remove best trasm from other criteria 
      (eg. split _C12 or multi align _G234)
      
cat pub3ig.trasm.match.tab4c | perl -ne \
'($gd,$td,$me,$ae)=split; ($mf,$ce,$xe)=split/[:\.]/,$me; ($af,$aw,$ai,$al)=split/[:,]/,$ae; 
$sc=1; map{ $v=$_/100; $sc *=$v; } ($ce,$xe,$ai,$al); ($gg=$gd)=~s/t\d+//; $gtd="$gd\t$td"; 
$gv{$gtd}=$_; $gs{$gd}{$td}=$sc; $gg{$gg}{$gtd}=$sc; pbest() if($lg and $lg ne $gg);  
$lgd=$gd; $lg=$gg; END{ pbest(); }
sub pbest{ @td=sort{$gg{$lg}{$b} <=> $gg{$lg}{$a} or $a cmp $b} keys %{$gg{$lg}};  
for $gt (@td) { ($gd,$td)=split"\t",$gt; $gv=$gv{$gt}; print $gv unless($did{$gd}++);
} }  ' 
> pub3ig.trasm.match.tab1c

=cut

sub put_gene_trasm_best { 
  my($outh,$lg,$gglg,$gv)= @_;
  our(%gdid,%tdid);  ## %gv,%gg,
  ## my %gglg= %{$gg{$lg}}; # copy ok?
  my @td=sort{ my($ga,$ta)= split"\t",$a; my($gb,$tb)= split"\t",$b; 
    $gglg->{$b} <=> $gglg->{$a}   # now get uninit value warn here
    or $tdid{$ta} <=> $tdid{$tb}  # or here**
    or $a cmp $b } keys %$gglg;  
  for my $gt (@td) { 
    my($gd,$td)=split"\t",$gt; my $v=$gv->{$gt}; 
    unless($gdid{$gd}++) { print $outh $v,"\n"; $tdid{$td}++; }
  }
}

sub pick_gene_trasm_best
{
  my($inh,$oktrids,$poortrids,$outh,)= @_; 
  unless(defined $outh) { $outh= *STDOUT; }
  my($lgd,$lg);
  our(%gv,%gg,%gdid,%tdid,%gglg);  
  %gv=(); %gg=(); %gdid=(); %tdid=(); %gglg=();
  while(<$inh>) {
    next unless(/^\w/); chomp; 
    my($gd,$td,$me,$ae)=split"\t"; # assume sorted by gene id gd input
    (my $gg=$gd)=~s/t\d+$//; # OPTION trsuffix here
    if($lg and $lg ne $gg) { put_gene_trasm_best($outh,$lg,\%gglg,\%gv); %gglg=(); %gv=(); }
    my($mf,$ce,$xe)=split/[:\.]/,$me;  $ce||=0; $xe||=0; # mapflag,cds-align,exon-align
    my($af,$aw,$ai,$al)=split/[:,]/,$ae; $ai||=0; $al||=0; # aaflag,aawidth,pident,peqsize
    my $sc=1; map{ my $pv=$_/100; $sc *= $pv; } ($ce,$xe,$ai,$al); 
    if(ref $oktrids and not $oktrids->{$td}) { $sc *= 0.5; }
    if(ref $poortrids and $poortrids->{$td}) { $sc *= 0.5; }
    ##no# if($tdid{$td}) { $sc *= 0.98; } #? or use in put sort?
    $tdid{$td}= 0 unless(defined $tdid{$td}); # stop warns
    my $gtd="$gd\t$td"; 
    $gv{$gtd}=$_; $gglg{$gtd}=$sc; # $gg{$gg}{$gtd}=$sc; #  $tid{$gtd}=$td; $gs{$gd}{$td}=$sc;
    $lgd=$gd; $lg=$gg;
  }
  put_gene_trasm_best($outh,$lg,\%gglg,\%gv) if($lg);
}



## UCK: this sub putclstr fails when enclosed in sub aaclusters 
#  : probably local hash tab mixup; using our (%aacluster,) might fix

sub putclstr {
  my($outh,$rrow,$trow,$aacluster,$eqgene)= @_;
  our ($nclust, $ndup, $nref, $ntr, $neqgene);
  #? @rrow= sort @$rrow; # sort by id; sort trow?
  my ($ref1)= shift @$rrow;  
  $nref++; $nclust++;
  my($refid,$rsize,$rident)= split"\t",$ref1;
  my $i=0; my $nrrow= @$rrow;
  foreach my $tr (@$rrow,@$trow) { 
    $i++; my($tid,$tsize,$pid)= split"\t",$tr; 
    if($pid == 0) { $pid=$rident; } # 
    $pid=~s/100.00/100/; $pid=~s/(\.\d)\d/$1/; 
    my $pdif= $tsize/$rsize; $pdif= 1.0/$pdif if($pdif>1.0);  $pdif= int(1000*$pdif)/10;
    my $pboth= int($pid * $pdif)/100; # use for sort
     #? FILTER here low qual?
    if($i <= $nrrow) { 
      $nref++;
      #eqgene was this:  ($pid>99 and $pdif>98) 
      $eqgene->{$tid}{$refid}= $eqgene->{$refid}{$tid} = $pboth if($pboth >= $MINGENEIDENT);
    } else {
      if($aacluster->{$refid} and my $pold= $aacluster->{$refid}{$tid}) { #? is this possible
        if($pboth < $pold) { next; } #$ndup++; 
      }
      $ntr++;
      $aacluster->{$refid}{$tid}= join"\t",$pboth,$pid,$pdif,$rsize,$tsize; # any dup rid,tid? 
      print $outh join("\t",$refid,$rsize,$tid,$tsize, $pid, $pdif)."\n" if(defined $outh);
    }
  }
}

sub aaclusters
{
  my($inh,$outh)= @_; 
  my (@rrow,@trow,%aacluster,%eqgene);
  our ($nclust, $ndup, $nref, $ntr,$neqgene)= (0) x 9;
  
  # ** Problems here 
  #** MISSING Thecc1EG000021t1,2 from cd100.aa.clstr
  
  while(<$inh>) {
    if (/^>/) {
      putclstr($outh,\@rrow,\@trow,\%aacluster,\%eqgene) if(@rrow>0 and (@rrow>1 or @trow > 0)); # ignore only @rrow? NO
      @rrow= @trow=();
    } elsif(/^\d/) {
      ## chomp();
      my($pident, $tsize, $tid)= (0) x 3;
      if (m/(\d+)(?:aa|nt), >(.+)\.\.\./) { # cd-hit pattern
         $tsize = $1; $tid = $2;
         $pident=(m/at ([\d\.]+)/)? $1 :0; # 0 for \*
      } else {
        die "cd-hit format error: $_";
      }
  
      my $val= join("\t",$tid,$tsize,$pident);
      if($tid =~ m/$IDPREFIX/) { push(@rrow, $val); } else { push(@trow, $val); }
    }
  }
  putclstr($outh,\@rrow,\@trow,\%aacluster,\%eqgene) if(@rrow>0 and (@rrow>1 or @trow > 0));
  
  # FIXME: crossref aacluster{IDPREFIX1}{IDPREFIX2} .. so all aac{IDPREF}{trasm} work
  # .. but use only ident/hisimilarity ?
  
  foreach my $gid (sort keys %eqgene) {
    next if(defined $aacluster{$gid}); #?? need to check each eqgene for tids? aac{gid}{tid} ??
    my @geq= sort{ $eqgene{$gid}{$b} <=> $eqgene{$gid}{$a} } keys %{$eqgene{$gid}};
    foreach my $geq (@geq) {
      if($aacluster{$geq}) {
        $neqgene++;
        # $aacluster{$gid}= $aacluster{$geq}; # is this enough?
        my @tid= sort keys %{ $aacluster{$geq} };
        foreach my $tid (@tid) {
          my $acval= $aacluster{$geq}{$tid};
          # my($pboth,$pid,$pdif,$rsize,$tsize)= split"\t",$acval;
          $aacluster{$gid}{$tid}= $acval;
        }
        last;
      }
    }
  }
  warn "#NOTE: aacluster: nclust=$nclust, ngene=$nref, ntr=$ntr, neqgene=$neqgene, ndup=$ndup\n" if $debug;
  return \%aacluster;
}

=item genes.overasm table : premake?

geneid            geneoid                       overasm-list    locus
Thecc1EG000001t1	AUGpier8a:AUGpier8ap1s_1g1t1	na	1sc:1369-2836
Thecc1EG000002t1	rna8b:r8L_g13025t00001	L_g13025t00001/100.99,P1_g05178t00001_C2/99.98,cacao3sopcsc1k89loc5659t1/99.97,P2_g03816t00001/99.96,caca11r39cuf8_Gsc1g253t1/99.92,B_g11013t00001/99.91,cacao3tri1sub3sc1loc106c0t1/99.87	1sc:7897-10405
Thecc1EG000005t1	mar7g.mar11f:AUGepir7p1s1g7t1	cacao3tri1sub3sc1loc558c0t2/99.90,caca11r39cuf8_Gsc1g128t2/99.89,cacao3v3sc1Loc777t15/95.80	1sc:17413-27097
Thecc1EG000005t2	vel4ma11:cacao3vel4sc1Loc938t2	cacao3vel4sc1Loc938t2/100.99,cacao3vel5sc1Loc924t2/100.99,cacao3vel4sc1Loc938t1/99.99,cacao3vel4sc1Loc938t10/99.99,cacao3vel4sc1Loc938t8/99.99,cacao3vel5sc1Loc924t1/99.99,cacao3vel5sc1Loc924t10/99.99,cacao3vel5sc1Loc924t8/99.99,cacao3vel4sc1Loc938t6/97.96,cacao3vel5sc1Loc924t6/97.96	1sc:18443-25907

=cut


  sub putgeneover { 
    my($overall,$outh,$lgd,$loid,$lov,$lloc)= @_;
    our($nref,$ntr,$ntrall);
    if( my $ovo= $overall->{$lgd} ) {  $lov="$ovo,$lov"; } # remove dups !!
    my @ov=grep{ /\w/ and $_ ne "na" } split",",$lov; 
    map{ s,/I100,/100.99,; s,/C,/,; } @ov;
    @ov=sort{ my $ca=($a=~m,/([\d\.]+),)?$1:0; my $cb=($b=~m,/([\d\.]+),)?$1:0; $cb <=> $ca or $a cmp $b } @ov; 
    $ntrall += @ov;
    my $lc=0; my %didov=();
    my @ovp= grep /\w/, map{ 
      my $cb=(m,/([\d\.]+),)?$1:0; $lc=$cb if($cb>$lc); 
      $_ if($cb>$lc-$DIFFMAX_OVERASM and not $didov{$_}++); 
      } @ov;
    my $ovp= join",", @ovp; $ovp||="na";  
    $overall->{$lgd}= $ovp; # only this? oid? lloc?
    print $outh join("\t",$lgd,$loid,$ovp,$lloc) if(defined $outh); 
    $nref++; $ntr += @ovp; 
  }
  
sub genes_overtrasm
{
  my($ovh,$outh)= @_; 
  # ovh = refgenes.over*asmrna.tab | sort -k1,1{geneid} -k3,3{asmid/score?} < dont need sort-k3,3
  my($lgd,$loid,$lov,$lloc, %overall);
  our($nref,$ntr,$ntrall)=(0) x 9;
# use MINGENEIDENT here ?
  # cat pub3ig.over*.tab | sort -k1,1 -k3,3 | perl -ne refgenes_overasm > pub3ig.allover.tab13
  # output to internal hash? file? both?
  
  while(<$ovh>) {
    next unless(/^\w/);
    my($gd,$oid,$ov,$loc)=split"\t";  # loc can be missing; oid needs column but not used here
    if($lgd and $gd eq $lgd) { 
      $lov.=",$ov" unless($ov eq "na" or not $ov =~ /\w/); 
    } else { 
      putgeneover(\%overall,$outh,$lgd,$loid,$lov,$lloc) if($lgd); 
      ($lgd,$loid,$lov,$lloc)=($gd,$oid,$ov,$loc); 
    } 
  }
  putgeneover(\%overall,$outh,$lgd,$loid,$lov,$lloc) if($lgd);
  my $nskip= $ntrall - $ntr;
  warn "#NOTE: genes_overtrasm: ngene=$nref, ntr=$ntr, ntrskip=$nskip\n" if $debug;
  return \%overall;
}

__END__

=item prelim work

# revise matchtsa scripts mess to 1 perl:
#  .. ignore for now prior gene-tsa annotation list; use this to recreate in reusable way?
#
#  -- start w/ pubgenes.over*.tab as source of valid asmrna transcripts for gene evidence
#   .. ? need some allowance for asmrna that covers gaps, scaffold ends/split scaffolds
#   .. but mark/downgrade asmrna with bad/poor utr, esp. gene joins; 
#      e.g split mapping to same region _C1/2 often joins or mangled asmrna
#  -- add aaequal scores to pubgenes.allover.tab
#   .. all from cdhit "cdallsame" tables? 
#   .. want some allowance for asmrna alt-tr where aa differs at end (middle?);
#   .. blastp gives better "mostlysame/different end" scoring than cdhit clusters
#  -- 2 pass or multi pass, 2nd to collect alt-tr only for main valid asmrna?

# sub refgenes_overasm
cat pub3ig.over*.tab | sort -k1,1 -k3,3 | perl -ne\
'($gd,$oid,$ov,$loc)=split"\t"; if($gd eq $lgd) { $lov.=",$ov"
unless($ov eq "na"); } else { putg() if($lgd); $lov=$ov; $lloc=$loc;
$loid=$oid; $lgd=$gd; } END{ putg(); } sub putg { @ov=grep{ /\w/ and
$_ ne "na" } split",",$lov; map{ s,/I100,/100.99,; s,/C,/,; } @ov;
@ov=sort{ ($ca)= $a=~m,/([\d\.]+),; ($cb)= $b=~m,/([\d\.]+),;
$cb<=>$ca } @ov; $lc=0; $ovp=join",", grep /\w/, map{ ($cb)=
m,/([\d\.]+),; $lc=$cb if($cb>$lc); $_ if($cb>$lc-5);  } @ov;
$ovp||="na";  print join("\t",$lgd,$loid,$ovp,$lloc); }' >
pub3ig.allover.tab13

=item FIXME id mismatches

  Thecc cdhit.aa has renames src id prefixes
  vs thecc.overrna ids
  
  ino: cat $caca/rnas/asmequal/pub3ig.over*.tab | sort -u -k1,1 -k3,3 > pub3ig.allover.tab1
  inc: cat cdall100/pub3i_trall9_cd{100,99,98,97,96,95,91}.aa.clstr > cdall100/pub3i_trall9_all.clstr

  #remapids:
perl -pi -e '@v=split"\t"; $_=$v[2]; \
s/cacao3sop/cacao5sop/g; s/cacao3tri/cacao5tri/g; s/cacao3vel12/cacao4vel12/g; \
s/cacao3v([123])/cacao3vel$1/g; s/caca11r39cuf/cacao3cuf/g; s/\b(B_|L_|P1_|P2_)/cacao3nwb$1/g; \
$v[2]=$_; $_= join"\t",@v;'  pub3ig.allover.tab1


=cut
