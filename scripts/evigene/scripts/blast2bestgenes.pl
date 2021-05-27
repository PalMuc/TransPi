#!/usr/bin/env perl
# blast2bestgenes.pl

=item about blast2bestgenes 

  geneset orthology assessment method, part of evidentialgene project
  
  read blastp tables of ref-genes (with ortho family assignments) x target gene sets
  compare targ gene sets by ortho group stats: bitscore, identity, alignment
  tabulate per orthogroup, gene set hits as best/same-as-best/diff-poorer/missing
  
  don gilbert, 2015.03, part of evidentialgene
  gilbertd at indiana.edu

=item RENAME this
  was blastsum.pl
  bl2bestgenes  ?
  blast2bestgenes <<
  makeblast2bestgenes?
  
=item usage

  -- add all in one step?
  
  step1: blastp.output table to bltab
  gunzip -c ../outz/sd-beetle4eset-refarp7s8set2.aa.blastp.gz | grep -v '^#' | \
   env blsum=1 aa="`ls aasetf/*.aa.qual`" ISCORE=3 nst=2 ./blastsum.pl \
    > beetle4enoculset-refarp7s8set2.bltab3

   env ISCORE=4 .. > beetle4enoculset-refarp7s8set2.bidtab3
    ^^ ident score better than bitscore for balanced table

    same format as makeblastscore2 tall but 3-columns added
      Query	Source	Bits	Iden	Algn	Qlen	Slen	Algb	Mism	Ibpt
    Algn and Idn are re-calculated for part-aligns after trimming overlapping hits
    Algb is blastp original align size (longer usually); Mism is blastp mismatches (ignoring indels)
    Ibpt is blast rank Ib, hit parts pt (eg: 1,1 ; 1,3 ; 99,1 ; ..)
    
  step2:
  cat beetle4enoculset-refarp7s8set2.bidtab3 | \
  env bestarp=1 ISCORE=4 DIFFSCORE=pct SAMEV=3 SAMEVABS=19  \
  SKIP="tribcas" alts="aasetf/tribcas_ncbi1406.aa.qual" arpids=refarp7s10fset1.omclgn.tab ./blastsum.pl \
   > beetle4enoculset-refarp7s8set2-bltab.arpbest66iden2

  step3: summary table
  
  cat beetle4enoculset-refarp7s8set2-bidtab.arpbest66iden2 | \
    env itab=1 ISCORE=5 nam=bidtab.arpbest66iden2 ./blastsum.pl \
    > beetle4enoculset-refarp7s8set2-bidtab.arpitab66algn
    
  cat .. | env itab=1 ISCORE=4 nam=bidtab.arpbest66iden2 ./blastsum.pl \
    > beetle4enoculset-refarp7s8set2-bidtab.arpitab66iden
  
  step4: per-orgroup table for plots (minor reformat of step2 table)
  cat beetle4enoculset-refarp7s8set2-bidtab.arpbest66iden2 | \
    env ieach=1 COMBEST=1 ISCORE=4 comall=allspp-refarp7s10fset1.comgrpid2 ./blastsum.pl \
      > beetle4enoculset-refarp7s8set2-bidtab.arpeach66algn2
 
  see R plot subs in arthropod/aabugs5/beetlet/tcas4evg/aaeval/beetle4aaevalR.txt
     .. maybe add here?

=item blastp ref-query x targeneset db needs step1a,b

  * if blast -db is target-genesets, -query ref-genes, then have a 
  problem adding hit parts then filter by score-sorted blast tables.
  -- need step1a= add hit parts (together in blast table), then 
  -- step1b= re-sort by 1:ref-query, 2:topscore
  -- add here in blastab2blsum() for swap=1 case ?
  
  # step1a: redo, nosort, nofilter 2nd hits:
  gunzip -c ../aaeval2/blastz/aaset3arp7-{apimel,apis,amel}*.blastp.gz | grep -v '#' |\
     env blsum=1 ISCORE=4 nst=999999 \
     aa="aasetf/apismel4set.aa.qual aasetf/refarp7s10fset1.aa.qual" \
     swap=1 keepho="amelevg|apimel|apismel" skipquery="honbee" \
     ./blast2bestgenes.pl \
      > apismel4set-refarp7s10fset1.bidtab4nofilt
  
  #step1b: resort by 1:refid, 4:ident (or/and 3:bits, 5:align), 2:targid
  sort -k1,1 -k4,4nr -k3,3nr -k2,2 apismel4set-refarp7s10fset1.bidtab4nofilt | perl -ne\
  '($rd,$td,$bs,$idn,$aln,$rw,$tw,$alg,$msm,$bip)=@v=split"\t"; ($ts,$tdd)=split":",$td,2;\
  print unless(++$did{$rd.$ts}>$NST); BEGIN{$NST=3;}' \
    > apismel4set-refarp7s10fset1.bidtab4filt

     
=item blastsum replaces

    env pctover=0.06 pmin=0.90  tall=1 evigene/scripts/makeblastscore2.pl xxx.blastp.gz > xxx.tall4
      ^ main problem this blast reader isnt adding hit parts as carefully as needed.
      env blsum=1 .. sub blastab2blsum() replaces
      
    makebestarp3.sh *.tall4 > geneset.arpbest63algn, .arpinfo63algn [ignore]
      env bestarp=1 .. sub makebestarp() replaces, about same method, improved
      
    makebestarp3cg.sh geneset.arpbest63algn > geneset.arpitab63algn
      env itab=1 .. sub makeitab() replaces, same summary tab for simple plot
      
    makebestarp3cgeach.sh geneset.arpbest63algn > geneset.arpeach63algn
      env ieach=1 .. sub makeieachtab() replaces, same per ogroup tab for detail plot
  
=item updates
  added part-overlap checks
  added aasize=aasetf/*.aa.qual 
    ---    

  FIXME: replace STDIN with file args
  
=cut

my $ISCORE=$ENV{ISCORE}||3; # 3,4,5 for  makebestarp; 0,1,2 for makeblsum?
#makeblsum was: $ISC=$ENV{iscore}||0; # bits=0, iden=1, algn=2
my $NST=$ENV{nst}||$ENV{NST}||2; # n-secondary target hits; need =0 option, or infinity=99999

### DAMN this here not hidden below:
  my $NEEDLEN= $ENV{skipnolen}||$ENV{needlen}||0; # skipnolen require id in aasize/blen or skip

my $makeblsum= $ENV{'blsum'}||0;

my $makebestarp= $ENV{'bestarp'}||0;
my $alttab=$ENV{alts}||"aasetf/tribcas_ncbi1406.aa.qual"; #bestarp
my $arpids=$ENV{arpids}||"refarp7s10fset1.omclgn.tab"; #bestarp
#unused# my $refset=$ENV{refset}||"NADA"; # for sum tables?

my $comall= $ENV{comall} || "allspp-refarp7s10fset1.comgrpid2";  #?? makebestitab need this
my $makebestitab= $ENV{'itab'}||$ENV{'sumtab'}||0; # summary table of bestarp.tab
my $makebesteach= $ENV{'ieach'}||0;

my $debug= $ENV{'debug'}||0;
my($si,$sbs,$sid,$saw,$saln,$smis,$np) = (0) x 10; # makeblsum
my(@ovs);

# revise: opt in/out files, use perl @ARGV  like perl -n .. while(<>) { .. }
my $INH= *STDIN; 
my $OUTH= *STDOUT;

sub MAINstub {}

my %alen= readSizes($ENV{aa}||$ENV{sizes}); # only for which steps?

if($makebestitab) { # step3
  #   arthropod/aabugs5/beetlet/tcas4evg/aaeval/makebestarp3cg.sh
  # PSIZE=1  # percent targ-aasize / ref-aasize, use this always
  # ONETAB=1  fixed
  $COMBEST=1; $COMBEST=$ENV{COMBEST} if(defined $ENV{COMBEST}); # v7, best/same/diff/miss,pbest for only comgrp subset
  $COMHITS=1; $COMHITS=$ENV{COMHITS} if(defined $ENV{COMHITS});
  # $comall= $ENV{comall} || "allspp-refarp7s10fset1.comgrpid2";  #?? need this
  # tset="$trgspplist" : ignore for this
  # env combest=$COMBEST COMHITS=$COMHITS psize=$PSIZE onetab=$ONETAB iscore=$ISCORE \
  #   tset="$trgspplist" nam=$outf perl -ne \
  
  readarplist($comall);
  makeitab($INH,$OUTH); #? fixme inh stdin or perl ARGV file open?

} elsif($makebesteach) { # step4
  $COMBEST=1; $COMBEST=$ENV{COMBEST} if(defined $ENV{COMBEST}); # v7, best/same/diff/miss,pbest for only comgrp subset

  readarplist($comall);
  makeieachtab($INH,$OUTH);
  
} elsif($makebestarp) { # step2
  $SKIP=$ENV{SKIP}||""; ## or skipset=tribcas .. honbee ..
  $SKIP2ND=$ENV{SKIP2ND}||0; 
  $SAMEV=$ENV{SAMEV}||9; 
  $SAMEVABS=$ENV{SAMEVABS}||19; # default 0 to skip?
  $DIFF2=$ENV{DIFFSCORE}||0; 
  $pFRAG= $ENV{FRAG}||75;
  
  ## add swapqt opt for input bltab w/ reversed target,refid cols, rw,tw cols also.
  
   # opts now are 0,2,"pct" .. add mixed pctSAMEV and SAMEV absolute diff >> SAMEVABS
   # e.g. bitscores 2653/2702 = 98.1% diff but 50 bits, ~25aa diff
   # same3.tribcas4evg2      ARP7f_G2956,dromel:FBgn0038693  tribcas4evg2:Tribca2aEVm000120t1 2702 1348    2603
   # same.tribcas4a  ARP7f_G2956,dromel:FBgn0038693  tribcas4a:TC010767tx1   2653    1310    2521
  
  readalts($alttab) if($alttab);
  readarptab($arpids);
  makebestarp($INH,$OUTH); #? fixme inh stdin or perl ARGV file open?
  
} elsif($makeblsum) {  # step1

  ## blastab2blsum uses these opts:
  # $skipho= $ENV{skipho} || ""; # pattern to skip homol match, e.g. aphid x swiss _ACYPI = same spp
  # $keepho= $ENV{keepho} || ""; # pattern to skip homol match, e.g. aphid x swiss _ACYPI = same spp
  # $onlyquery= $ENV{onlyquery} || ""; # species query to keep, or all
  # $swapqt = $ENV{swap} || 0;

  blastab2blsum($INH,$OUTH);  #? fixme inh stdin or perl ARGV file open?
  
} else { warn 
  "#usage: env blsum=1 OR bestarp=1 OR sumtab/itab=1 OR ieach=1 [opts] blastsum.pl ..\n".
  "# blsum opts  : aa=file.aa.qual prot sizes, nst=$NST 2ndary hits\n".
  "# bestarp opts: ISCORE=3|4|5 (bits,idn,aln), DIFFSCORE=pct|0|2, SAMEV=9%, SKIP=skip_refspp\n".
  "#   alts=file.trid[col0]..geneid[col7], arpids=ARPid,refid,..\n";
}  

#------------------------

# makebestarp3.sh for blsum tab instead of tall tabs

#? export SAMEV=9;
#? export SKIP2ND=0;
#? export DIFF2SCORE=0 DIFF2SCORE=2 

sub readSizes {
  my($inf)= @_; my $n=0;
  my %alen=(); return %alen unless($inf);
  my @aaf=split/[,\s]+/,$inf;
  for my $aaf (@aaf) { open(F,$aaf); while(<F>){ my($id,$al)=split; $alen{$id}=$al; $n++; } close(F); }
  warn  "# read n=$n from $inf\n" if $debug;
  return %alen;
}

sub readalts {
  my($inf)= @_; my $n=0;
  # FIXME: problem table format case: cols 0,6 expected as trid,geneid
  my @aaf=split/[,\s]+/,$inf;
  for my $aaf (@aaf) { open(F,$aaf) or warn "# Missing altf: $aaf\n";
    while(<F>) { my($td,$gd)=(split)[0,6]; if($gd){ $altof{$td}=$gd; $n++;} } close(F);  
  }
  warn  "# read n=$n from $inf\n" if $debug;
  # return %altof;
}

# refarp7s10fset1.omclgn.tab  n=131814, 2014.oct.20
# ARP7f_G3	daphmag:Dapma7bEVm016249t1	7,406
# ARP7f_G3	daphmag:Dapma7bEVm018856t1	7,406
sub readarptab {
  my($inf)= @_; my $n=0;
  open(F,$inf) or warn "# Missing arptab: $inf\n";
  while(<F>) { if(/^ARP/) { my($ad,$rd,$txgn)=split; $rad{$rd}=$ad; $n++;} } close(F); 
  warn  "# read n=$n from $inf\n" if $debug;
  # return %rad;
}

sub readarplist {
  my($inf)= @_; my $n=0;
  open(F,$inf) or warn "# Missing arptab: $inf\n";
  while(<F>) { if(/^ARP/) { my($ad)=split; $arpid{$ad}=1; $n++; } } close(F); 
  warn  "# read n=$n from $inf\n" if $debug;
  # return %arpid;
}


=item makeieachtab 

   minor reformat of bestarp, same col format as from
    arthropod/aabugs5/beetlet/tcas4evg/aaeval/makebestarp3cgeach.sh
   to generate per-orgroup graphs of best/poor gene sets
    * prior output was sorted by tgeneset, now is input.tab order by best. reorder for plots?
    
  input bestarp.tab: beetle4enoculset-refarp7s8set2-bidtab.arpbest66iden2
best3.tribcas14nc	ARP7f_G2613,dromel:FBgn0053196	tribcas14nc:XP_008197897	20558	10622	17554	22949	18259	17796606	2,3
diff.tribcas4evg2	ARP7f_G2613,dromel:FBgn0053196	tribcas4evg2:Tribca2aEVm000001t1	19031	10241	18259	22949	182718624	7601	1,3
diff.tribcas1	ARP7f_G2613,dromel:FBgn0053196	tribcas1:TC011986	17631.7	9855	19264	22949	21117	21202	9098	5,5
diff.tribcas4a	ARP7f_G2613,dromel:FBgn0053196	tribcas4a:TC011986tx2	16921	9613	17838	22949	19064	19732	8010	6,3
  
  output format: beetle4enoculset-refarp7s8set2.arpeach63algn
tgeneset    	orgroupid	refgeneid	trgeneid	ng	bits	algn	pal	dlen	rsize	tsize	isbest
tribcas1    	ARP7f_G2613	dromel:FBgn0053196	tribcas1:TC011986	1	18057	21049	91.7	92.0	22949	2111best
tribcas14nc 	ARP7f_G2613	dromel:FBgn0053196	tribcas14nc:XP_008197897	1	20558	17795	77.5	79.6	229418259	diff
tribcas4a   	ARP7f_G2613	nasvit:Nasvi2EG006740t1	tribcas4a:TC011986tx1	1	17436	17647	76.9	76.2	16711	1748diff
  
=cut
 
sub makeieachtab {
  my($inh,$outh)= @_;
  my $hdr=0;
  $ISC=$ISCORE-1; $IS2=($ISC==4)?2:4; 
  $PSIZE=1;  # percent targ-aasize / ref-aasize, use this always
  my($larpid,$lrid,@tgv);
  while(<$inh>) {
    next if(/^\W/); chomp;
    my($cl,$rd,$td,$bs,$idn,$al,$rw,$tw)=split"\t"; 
    $bs=int($bs); $rw=1 if($rw<1); # ERROR if missing val
    my($arpid,$rsppid)= split",",$rd,2; 
    my($cln,$cls)=split/\./,$cl,2; $cln=~s/\d$//; 
    my $ts=$td; $ts=~s/:.*//; if($td eq "miss") { ($ts=$cl)=~s/^\w+.//; }
    my $iscom=$arpid{$arpid}||0; 
    next if($COMBEST and not $iscom); #?? need this;  
    # my $rs=$arpid; ## change hashes to each arpid
    # my $score=($ISC==4)?$al:($ISC==3)?$idn:$bs; 
    
    my $dw= ($PSIZE)? (100*$tw/$rw) : $tw - $rw; 
    my $ng= 1; # $ng{$arpid}{$ts}||0; 
    # my $pal=$pal{$arpid}{$ts}; # recalcd above, rw1
    my $pal= 100*$al/$rw; 
    map{ $_= sprintf "%.1f",$_; $_=100 if($_ >= 100); } ($pal,$dw); #>100 only for dw.PSIZE, pal
    
    $rsppid ||= $lrid||"ridmiss";
    ## add idn col ??  drop ng col? add iscom FLAG to output tab
    #orig# my @col= ($arpid, $rsppid, $td, $ng, $bs, $al, $pal, $dw, $rw, $tw, $cln);
    my @col= (sprintf("%-12s",$ts), $arpid, $rsppid, $td, $bs, $idn, $al, $pal, $dw, $rw, $tw, $cln, $iscom);
    unless($hdr++) { 
      #orig# my @hdr=qw(orgroupid refgeneid trgeneid ng bits algn pal dlen rsize tsize isbest);
      my @hdr=qw(orgroupid refgeneid trgeneid bits iden algn pal dlen rsize tsize isbest comorgrp);
      unshift @hdr, sprintf( "%-12s","tgeneset");
      print $outh join("\t",@hdr)."\n"; 
      } 
    
    ## fixme here collect all/arpid and reorder by ts
    # printf $outh "%-12s\t",$ts; print $outh join("\t",@col)."\n";  
    if($larpid ne $arpid) {
      map{ print $outh join("\t",@$_)."\n"; }sort{$a->[0] cmp $b->[0]}@tgv;
      @tgv=();
    }
    push @tgv, \@col;  
    $larpid=$arpid; $lrid=$rsppid;
  }
  map{ print $outh join("\t",@$_)."\n"; }sort{$a->[0] cmp $b->[0]}@tgv;
  
}

sub makeitab {
  my($inh,$outh)= @_;
  # < $comall $intabs > $info
  my $hdr=0;
  
  # BEGIN
  $ISC=$ISCORE-1; $IS2=($ISC==4)?2:4; 
  # $TSET=$ENV{tset}; @TSET=split/\W+/,$TSET if($TSET);
  # $refset=$ENV{refset}||"NADA"; $refset=~s/[\s,]+/|/g; 
  #above# $COMHITS=$ENV{COMHITS}||0;
  #above# $COMBEST=$ENV{combest}||0; 
  $ONETAB=1; # fixed $ENV{onetab}||0; 
  $PSIZE=1; # fixed $ENV{psize}||0;
  my $nam=$ENV{nam}||"noname"; 
  my(%dups,%frags);
  print $outh "# summary for $nam, iscore=$ISC\n"; 
  #x %comid= %arpid; #readarplist
  while(<$inh>) {
    if(/^(ARP\w+)$/){ $arpid{$1}=1; $ncom++; next; } # see above readarplist
    next if(/^(Class|\W)/); 
    chomp; my @v=split"\t"; $inline=$_;
    my($cl,$rd,$td,$bs,$idn,$al,$rw,$tw)=@v;
    my $dupfrag= $v[-1]; #  "9.dup,1.frag" ; 1512 upd flag, fixme hack
    
    my $score=($ISC==4)?$al:($ISC==3)?$idn:$bs; 
    if($ISC==3){ ($al,$idn)=($idn,$al); } 
    ($rid=$rd)=~s/,.*//; 
    $iscom=$arpid{$rid}||0; $grpid{$rid}++; $comidhit{$rid}++ if($iscom);
    $rs=$rd; unless($rs=~s/[_:].*//) { $rs="ARPgroup"; } 
    ($ts=$td)=~s/:.*//; if($td eq "miss") { ($ts=$cl)=~s/^\w+.//; }
    ($cln,$cls)=split/\./,$cl; 
    ## 1512 fixme: count cl==diff2,diff3 as tiny/partial/fragment hits
    ##not now.. $cln="dtiny" if($cln=~/diff[234]/);
    $cln=~s/\d$//; $cln{$rs}{$cln}{$cls}++;
    $hit=($td =~ /miss/ or $score < 1)?0:1; $rw||=1;
    my $dw= ($PSIZE)? (100*$tw/$rw) : $tw - $rw; 
    my $pal=($hit)?(100*$al/$rw): 0; $pal=100 if($pal>100);
    $rs{$rs}++; $nr{$rs}++ unless($did{$rid}++); $ts{$ts}++; $ng{$rs}{$ts}++ if($hit);
    $bs{$rs}{$ts}+=$bs; $al{$rs}{$ts}+=$al; $dw{$rs}{$ts} += $dw; $pal{$rs}{$ts} += $pal;
    $dups{$rs}{$ts}++ if($inline=~/[1-9]\.dup/); ## dupfrag
    $frags{$rs}{$ts}++ if($inline=~/1\.frag/);
    
    if($hit) { $ts.="hit";  
      $bs{$rs}{$ts}+=$bs; $al{$rs}{$ts}+=$al; $dw{$rs}{$ts} += $dw; $pal{$rs}{$ts} += $pal; }
    if($iscom){ $ts=~s/hit$//; $ts.="com"; $ng{$rs}{$ts}++ if($hit); $cln{$rs}{$cln}{$ts}++;
      $bs{$rs}{$ts}+=$bs; $al{$rs}{$ts}+=$al; $dw{$rs}{$ts} += $dw; $pal{$rs}{$ts} += $pal; } 
  }

  # END
  for $r (sort keys %nr) {
   my $hdr=0; $nr=$nr{$r}; 
   @ts=sort keys %ts; 
   #off# @ts=@TSET if(@TSET);
   @cln= qw(best same diff miss); #? 
   # @cln=sort keys %{$cln{$r}}; #?? 1512 add dtiny
   $ncom=scalar(keys %arpid); 
   $ngrp=scalar(keys %grpid);
   $ncomhit= scalar(keys %comidhit); 
   print $outh "# ref: $r, ngene=$nr, ngroup=$ngrp, ncomgrp=$ncom, ncomgrphit=$ncomhit [$COMHITS] ................\n" ;
   $ncomn= ($COMHITS)?$ncomhit:$ncom;
   for my $ts (@ts) { $ts0=$ts;
     $ng=$ng{$r}{$ts}; $bs=$bs{$r}{$ts}; $al=$al{$r}{$ts}; $dw=$dw{$r}{$ts}; $pal=$pal{$r}{$ts};
     @av= map{ sprintf "%.1f",$_/$nr; } (100*$ng, $bs, $al, $pal, $dw);
     $ts.="hit"; $bs=$bs{$r}{$ts}; $al=$al{$r}{$ts}; $dw=$dw{$r}{$ts}; $pal=$pal{$r}{$ts};
     @avhit=map{ sprintf "%.1f",$_/$ng; } ( $bs, $al, $pal, $dw);
     $ts=~s/hit/com/; $nc=$ng{$r}{$ts}; $bs=$bs{$r}{$ts}; $al=$al{$r}{$ts}; $dw=$dw{$r}{$ts}; $pal=$pal{$r}{$ts};
     @avcom=map{ sprintf "%.1f",$_/$ncomn; } ( 100*$nc, $bs, $al, $pal, $dw);
     $png=shift @av; $pncom=shift @avcom; 
     @col= ($ng,$png,$pncom,@av,@avhit,@avcom);
     $ts=$ts0;
     #always: if($ONETAB) 
      if($COMBEST) { $ts.="com"; }#FIXME: add/swap best/poor for all vs for commongrp: ts.com
      @cltab=(); $tc=$pbs=$pdf=0; 
      $pdup=$dups{$r}{$ts0}; $pfrag=$frags{$r}{$ts0}; # 1512 upd
      for $cl (@cln) { $c=$cln{$r}{$cl}{$ts}||0; $tc+=$c; push @cltab, $c; 
        if($cl=~/(diff|miss)/) { $pdf+=$c; } else { $pbs+=$c; } } 
      map{ $_= sprintf "%.1f",100*$_/$tc;} ($pbs,$pdf,$pdup,$pfrag);
      push @cltab, $pbs, $pdf, $pdup, $pfrag; 
      push @col, @cltab;     
     
     unless($hdr++) { 
       ## change from align cols to ident when ISC==3
       if($ISC==3) {  @hdr=qw(tng png pncom bits iden pid dlen bhit ihit phit dhit bcom icom pcom dcom); }
       else {  @hdr=qw(tng png pncom bits algn pal dlen bhit ahit phit dhit bcom acom pcom dcom); }
       push @hdr, @cln,"pbest","ppoor","pdup","pfrag" if($ONETAB);
       printf $outh "%-12s\t","tspp"; print $outh join("\t",@hdr)."\n"; 
       } 
     printf $outh "%-12s\t",$ts0; print $outh join("\t",@col)."\n"; # $ng,$png,$pncom,@av,@avhit,@avcom
   }
  }    
}
#    unless($ONETAB) { ## off
#    print "#----------------------\n"; 
#    @ts=sort keys %ts; @ts=@TSET if(@TSET);
#    printf "%-12s\t","tspp"; print join("\t",@cln,"%best","%poor")."\n";
#    for my $ts (@ts) { 
#       @cltab=(); $tc=$pbs=$pdf=0; 
#       for $cl (@cln) { $c=$cln{$r}{$cl}{$ts}||0; $tc+=$c; push @cltab, $c; 
#         if($cl=~/(diff|miss)/) { $pdf+=$c; } else { $pbs+=$c; } } 
#       map{ $_= sprintf "%.1f",100*$_/$tc;} ($pbs,$pdf);
#       push @cltab, $pbs, $pdf;      
#       printf "%-12s\t",$ts; print join("\t",@cltab)."\n"; 
#     }
#    print "#----------------------\n\n"; 
#    }

  

# orig: gunzip -c $tallz | sort -k${ISCORE},${ISCORE}nr -k2,2 -k1,1 | cat $arpids - | makebestarp.sh
# blastsum.tab same 7 cols as tall + others; use input sort order (k1=refid, besthits/ref), 
# collect all inputs then decide bestarp from refids/arp, using arp-ref w/ largest ISCORE (bits,iden,aln)

sub makebestarp {
  my($inh,$outh)= @_;
  $ISC=$ISCORE - 1; $ISC=4 if($ISC>4 or $ISC<2); 
  $IS2=($ISC==2)?4:2; # was ($ISC==4)?2:4; 
  my $swapqt = $ENV{swap} || 0;

  # globals: my(%rv,%rvscore,%tss,%rw);
  while(<$inh>) { 
    next if(/^\W/);
    chomp; my @v=split"\t"; 
    if($swapqt) { @v[0,1]= @v[1,0]; @v[5,6]= @v[6,5]; } # ($rd,$td)= ($td,$rd); ($rw,$tw)= ($tw,$rw);
    my($rd,$td,$bs,$idn,$al,$rw,$tw)=@v; 

    # @v=@v[0..6];  # chop extra blsum cols for output? or not?
    my $score=$v[$ISC]; # bs, idn or al
    my $gd= geneid($td); # $altof{$td}||$td; $gd=~s/t\d+$//; $gd=~s/tx\d+$//; $gd=~s/\-[RP]\w+$//;
    next if($SKIP and $rd=~/$SKIP/); 
    ## see also SKIP/skipquery/skipho in blastab2blsum
    
    ## blsum: remove did{gd} filter here, do in putarp
    next unless($td and $score>0); # ( and not $did{$gd});  #  $gd here? was td
    
    $rad=$rad{$rd} or next; 
    my $ts=$td; $ts=~s/:.*//; 

    $v[0]= "$rad,$v[0]"; 
    $rv2nd{$v[0]}{$td}=[@v]; # save all for below?

  # UPD1512: add dup-gene scoring (eg. busco 1copy dups..) and frag score (diff << criteria)
  # use rv2nd, dup added table, criterion is score < rvscore.best but score > 0.95? of best
  # .. AND td align/ident level close to ref to count as valid dup

    unless($rv{$rad}{$ts} and $score <= $rvscore{$rad}{$ts}) { # $rv{$rad}{$ts}->[$ISC]
      #above# $v[0]= "$rad,$v[0]"; # v[0] == rd == refid
      # $did{$td}++; $did{$gd}++;
      $rv{$rad}{$ts}=[@v]; $tss{$ts}++; $rw{$rad}=$rw if($rw > $rw{$rad});
      $rvscore{$rad}{$ts}= $score;
      $rvscore{$rad}{'best'}= $score if($score>$rvscore{$rad}{'best'});
    } 
    else { # 1512 dup counts
      if($score > 0.95 * $rvscore{$rad}{$ts} and $al/$rw > 0.90) { # dup criterion, fixme
	 # FIXME: altof
         my $bd=$rv{$rad}{$ts}->[1];
         my $bg=$altof{$bd}||$bd; my $tdg=$altof{$td}||$td;
         $rvscore{$rad}{$ts.".dup"} ++ unless($bg eq $tdg);
         #? $rv{$rad}{$ts}->[999]++; # =[@v]; flag in @v ?? or rvscore{}{'dup'} ?
      }
    }
    # elsif($rv{$rad}{$ts} and $score <= $rvscore{$rad}{$ts}) { # $rv{$rad}{$ts}->[$ISC]
    #   # mark as 2nd best, dont use for other group * maybe wrong; try not
    #   # FIXME: makebest? putarp? need to use same refgenes/arpid for all, though some have best score on diff refg
    #   # $v[0]= "$rad,$v[0]"; 
    #   # $rv2nd{$rad}{$td}=[@v]; # {$ts} # not used?
    #   # if($SKIP2ND) {$did{$td}++; $did{$gd}++;}
    # }
  } 
  putarp($outh); 
  print $outh "# ISCORE=$ENV{ISCORE}, SAMEV=$SAMEV, DIFF2SCORE=$DIFF2, SKIP2ND=$SKIP2ND\n" ;   
}

sub geneid {
  my($td)=@_;
  # BUG: gd needs geneset prefix from td ..
  my($ts,$tdd)=split":",$td,2;
  my $gd=$altof{$td}; if($gd){ $gd="$ts:$gd"; } else { $gd=$td; } 
  $gd=~s/t\d+$//; $gd=~s/tx\d+$//; $gd=~s/\-[RP]\w+$//; # various alt tags
  # add -P[ABCD..]; add genetag tribcas4a:TC000tx1,2; problem of NCBI alts no idtag
  return $gd;
}

# putarp bestarp.tab
# best3.tribcas4evg2	ARP7f_G1000,honbee:Apimel3aEVm003441t1	tribcas4evg2:Tribca2aEVm003531t1	699.9	430	549	842	611
# diff.tribcas1	ARP7f_G1000,honbee:Apimel3aEVm003441t1	tribcas1:TC014149	572.9	343	409	842	534
# diff.tribcas4a	ARP7f_G1000,honbee:Apimel3aEVm003441t1	tribcas4a:TC014149tx1	562.9	338	404	842	529
# diff.tribcas14nc	ARP7f_G1000,honbee:Apimel3aEVm003441t1	tribcas14nc:XP_008193319	597	348	395	842	393
# same3.tribcas1	ARP7f_G10002,nasvit:Nasvi2EG004171t1	tribcas1:TC015373	52	57	195	191	194
# same.tribcas14nc	ARP7f_G10002,nasvit:Nasvi2EG004171t1	tribcas14nc:XP_975784	52	57	195	191	194
# same.tribcas4a	ARP7f_G10002,nasvit:Nasvi2EG004171t1	tribcas4a:TC015373tx1	52	57	195	191	194
# same.tribcas4evg2	ARP7f_G10002,nasvit:Nasvi2EG004171t1	tribcas4evg2:Tribca2aEVm009815t13	50.8	57	195	191	200

sub putarp { 
  my($outh)= @_;
  # decide bestarp here from refids/arp, using arp-ref w/ largest ISCORE (bits,iden,aln)
  my (%did);
  my @tss=sort keys %tss; my $ntss=@tss;
  for my $rad (sort { $rvscore{$b}{'best'} <=> $rvscore{$a}{'best'} or $a cmp $b} keys %rvscore) {
  
    my @mis=(0)x7; @mis[0,1,5]=($rad,"miss",$rw{$rad}); my $mis= \@mis;
    my @rv; my $nmiss=0;
    
    
    for $ts (@tss) { 
      my $v= $rv{$rad}{$ts}; 
      my $gd= (ref $v) ? geneid($$v[1]) : 0;
      # $did{$td}++; $did{$gd}++; #<< add here skip gene if used in other rad
      if($gd) {
        if($did{$gd}) { $v= undef;
          #? my @td2= sort keys %{$rv2nd{$rad}}; # sort by score!
          #? for my $t (@td2) { $v= $rv2nd{$rad}{$t}; my $gd= geneid($$v[1]); }
        } else { $did{$gd}++; }
      }
      unless($v) { $v=[@mis]; $$v[1]="$ts:miss"; $nmiss++; } 
      else {
        my $dupv= $rvscore{$rad}{$ts.".dup"}||0;
        my $vflag= $dupv.".dup"; # this replaces '1,1' col not used now ?
        my($rw,$tw)= @$v[5,6];  # ($rd,$td,$bs,$idn,$al,$rw,$tw)=@v;
        my $fragv= ($tw < $rw*$pFRAG/100)?1:0;
        $vflag .=",$fragv.frag";
        $v->[-1]= $vflag;  ## "9.dup,1.frag
      }
      push @rv,$v; 
    }
    if($nmiss >= $ntss) {
      # print $outh "#missall: $rad\n"; 
      next;
    }   
    
    @rv=sort{ $$b[$ISC]<=>$$a[$ISC] or $$b[$IS2]<=>$$a[$IS2]} @rv;

  # UPD1512: add dup-gene scoring (eg. busco 1copy dups..) and frag score (diff << criteria)
  #
  # FIXME: makebest? putarp? need to put same refgenes/arpid for all tss, 
  # though some have best score on diff refg
  # fix here, after rv sort. 
  # BUT changes scores, more diff/best, fewer same >> ~100 more evg2trib best, 10% increase, from same using all-refid
  # .. not perfect fix, some cases where 2nd targets lack rad1 in blsum table, due to output limit NST=2
      
    my($df,$sm,$i,$ts)=(0,0);  my @sd=(0)x3; 
    my @bv=@{$rv[0]}; 
    my @outv=(\@bv);
    my $rad1= $bv[0]; my $rw=$bv[5];
    for $i (1..$#rv) { 
      my @v= @{$rv[$i]}; 
      if($v[0] ne $rad1){ my $td=$v[1]; my $v2= $rv2nd{$rad1}{$td}; @v=@$v2 if($v2); } # ensure match to same arpid/refgene
      push @outv, \@v; # possible change
      my $bvi=$bv[$ISC];
      my $sd= $bvi - $v[$ISC]; my $sdabs=$sd;
      my $howdif=0; my $psd= ($bvi>0) ? 100*$sd/$bvi : 0;
      if($DIFF2 =~ /^p/) { $sd= $psd; } # pct diff
      elsif($DIFF2==2) { $sd += ($bv[$IS2] - $v[$IS2]); }
      if($v[1] =~ /miss/ and $sd>0) { $df++; $sd[$i]="miss"; }
      elsif($sd > $SAMEV){ $df++; $sd[$i]="diff";  $howdif=1; } 
      elsif($SAMEVABS and $sdabs > $SAMEVABS) { $df++; $sd[$i]="diff"; $howdif=2; } # eg for largish abs diff, smallish pct for large genes
      else { $sd[$i]="same"; $sm++; } 
      #if($howdif>0 and $psd < $pFRAG) { # 1512 tiny flag
      #  $sd[$i]="diff2"; #? ($psd < 0.66 * pFRAG)?"diff3":"diff2"; # 3 levels of diff now, 1=near same, 2=75%, 3=50%
      #  ## urk need pFRAG flag independently of best/same/diff, this is relative to ref not to best-target
      #}
      }
    if($sm>0) { $sd[0]="same$sm"; } else { $sd[0]="best$df"; } 
    # some ogroups are all-miss, drop** approx 6,000 above 10,000 found
    # some ogroup scores very low (weak align), drop?
    for my $i (0..$#outv) { 
      my @v=@{$outv[$i]}; #was @{$rv[$i]}; 
      ($ts=$v[1])=~s/:.*//; ## $dt{$ts}=1; 
      print $outh join("\t","$sd[$i].$ts",@v)."\n"; 
      }
  }
}

# blastab2blsum
# in: blast table of ref-query x target-subject, outfmt=6,7, 
# out: sum table ref-query, target-subject, bits, iden, align, adding blast align parts/subject

sub blastab2blsum {
  my($inh,$outh)= @_;

  $IBLIDX=5; # global for scores sort, blast order index, 1st = best single hit
  $ISC=$ISCORE - 3; #?? use same 3,4,5 ISCORE range as for makebestarp?
  $ISC=0 if($ISC>2 or $ISC<0); 
  my $swapqt = $ENV{swap} || 0;
  ##ABOVE NOW: my $NEEDLEN= $ENV{skipnolen}||$ENV{needlen}||0; # skipnolen require id in aasize/blen or skip
  $pMINLOW   = $ENV{pmin} || 0.20; # ?? ADD this, or some low qual 2nd hit filter
  $pMINBEST2 = $ENV{pmin2nd} || 0.0; # ADD, skip 2nd hits if below 0.5*topscore; replace $NST ?
  # my $showSPAN= $ENV{spans}||0; # 201402, TALL only options output targ,src align spans; make default?
  
  ## also SKIP option == skipquery, see above: if($SKIP and $rd=~/$SKIP/); 
  my $skipho= $ENV{skipho} || $ENV{skiptarg} ||""; # pattern to skip homol match, e.g. aphid x swiss _ACYPI = same spp
  my $keepho= $ENV{keepho} || $ENV{keeptarg} || ""; # pattern to skip homol match, e.g. aphid x swiss _ACYPI = same spp
  my $onlyquery= $ENV{onlyquery} || ""; # species query to keep, or all
  my $skipquery= $ENV{skipquery} || $ENV{SKIP} || ""; # species query to keep, or all
  map{ s/[,\s]+/\|/g; } ($skipho,$keepho,$onlyquery,$skipquery);
  
  %did=(); # $did{$lq.$ts} # clear this global..
  my($sbs,$sid,$saw,$saln,$smis,$si,$np)= (0) x 9;
  my $ltval= [$sbs,$sid,$saw,$saln,$smis,$si,$np];
  my $blerr=0;
  
  while(<$inh>) { 
    next if(/^\W/);
    chomp; my @v=split"\t"; 
    unless(@v==12){ warn"ERR: blasttab not 12 cols:'@v'\n"; $blerr++; die if ($blerr>9); next; }

    my($q,$t,$pi,$aln,$mis,$indl,$rb,$re,$tb,$te,$ev,$bs)=@v; # blast table columns, outfmt=6/7
    #o# my($q,$t,$bits,$aln,$mis,@bspan)=  @v[0,1,11,3,4, 6,7,8,9]; ## 6-9 =  q. start, q. end, s. start, s. end, 
    
    ## SORT problems here using swap .. wrong order then for putq() filter by NST
    ## need input ordered by 1:refid, 2:top-bitscore + 3:trgid..
    
    if($swapqt) { ($q,$t)= ($t,$q); ($rb,$re,$tb,$te)= ($tb,$te,$rb,$re); }
    
    ## filter for query/target species, as per makeblasttab2
    next if($skipquery and $q =~ m/$skipquery/);  
    next if($onlyquery and $q !~ m/$onlyquery/); # need rev: skipquery
    next if($skipho and $t =~ m/$skipho/);
    next if($keepho and $t !~ m/$keepho/);
    next if($NEEDLEN and not( $alen{$q} and $alen{$t}));
    
    my $o=1; if($tb>$te) { $o=-1; ($tb,$te)=($te,$tb); } 
    if($q ne $lq) { putt($lt,$lq,$ltval); putq($lq,$outh); $lt=$si=0; } 
    if($t ne $lt) { putt($lt,$q,$ltval); @ovs=(); $np=$saln=$smis=$sbs=$sid=$saw=0; $si++; } 
    my($ok,$trimb)=(0,0);
    ($ok,$rb,$re,$tb,$te,$trimb)= ovtrim($rb,$re,$tb,$te);
    if($ok) {
      $np++; $rw=1+$re-$rb; $tw=1+$te-$tb; $aw=($rw<$tw)?$rw:$tw; 
      if($trimb>0) { $bs -= $trimb; $bs=0 if($bs<0); } #?? 2*trimb? or 1.5*trimb?
      $saw+=$aw; $sbs+=$bs; $saln+=$aln; $smis+= $mis; # or $mis+$indl ?
      $bi=int(0.5+ $aw*$pi/100); $sid+=$bi; 
    }
    ($lq,$lt)=($q,$t);
    $ltval= [$sbs,$sid,$saw,$saln,$smis,$si,$np];
  } 
  putt($lt,$lq,$ltval); putq($lq,$outh); # END
}

sub putt { 
  my($lt,$lq,$ltval)=@_;  
  $sum{$lq}{$lt}=$ltval if($lt and $lq); # [$sbs,$sid,$saw,$saln,$smis,$si,$np];  # == $ltval
  # my($ts,$td)=split":",$lt,2; 
  # if(++$did{$lq.$ts} <= $NST) 
} 

sub bint { my $b=shift; return ($b<0) ? 0 : ($b=~/e\+/)? int($b) : int(0.49+$b); } # drop .decimals and e-001 vals
#orig# sub bint { my $b=shift; return ($b<0) ? 0 : ($b=~/e\+/)? int($b) : $b; }

sub putq {
  my($lq,$outh)=@_; return unless($lq);
  my %sumq= %{$sum{$lq}}; delete $sum{$lq}; 
  my $topaln= 0;
  for my $lt (sort{$sumq{$b}->[$ISC] <=> $sumq{$a}->[$ISC] or
      $sumq{$a}->[$IBLIDX] <=> $sumq{$b}->[$IBLIDX] or $a cmp $b} keys %sumq) {
    my($ts,$td)=split":",$lt,2; 
    next if($NST and ($did{$lq.$ts} >= $NST));
    
    my($bs,$id,$aw,$aln,$mis,$i,$p)=@{$sumq{$lt}}; 
    my $qlen=$alen{$lq}||0; my $slen=$alen{$lt}||0; # add: Qlen  Slen from aa.qual tables?
    
    ## use these to filter low qual 2nd hits: $pMINLOW $pMINBEST2
    ## but uncertain here that lq is ref?
    my $paln= ($qlen>0)? $aln/$qlen : 0;
    # next if($pMINLOW and $paln>0 and $paln < $pMINLOW); # absolute: Algn/Reflen < pMINLOW
    
    ## as for NST and $did{$lq.$ts}
    # next if($pMINBEST2 and $did{$lq.$ts} and $paln>0 and $paln < $pMINBEST2 * $topaln); # rel to top hit
    ### next if($pMINBEST2 and $paln>0 and $paln < $pMINBEST2 * $topaln); # rel to top hit
    ## ^^ respect $ts species, apply pMINBEST2 only for same species? ie keep 1 of each at least
    
    $bs= bint($bs);
    print $outh join("\t",qw(Query Source Bits Iden Algn Qlen Slen Algb Mism Ibpt))."\n" unless($HDR++); 
    print $outh join("\t",$lq,$lt,$bs,$id,$aw,$qlen,$slen,$aln,$mis,"$i,$p")."\n"; 
    ++$did{$lq.$ts};
    $topaln= $paln if($paln>$topaln);
  } 
} 

sub ovtrim { 
  my($qb,$qe,$sb,$se)=@_;
  my @inv=($qb,$qe,$sb,$se); 
  my $ov=0; my $maxt=0; my $ok=1;
  for my $sp (@ovs) {
    my($qd,$sd)=(0,0);
    my($xb,$xe,$tb,$te)= @$sp;   
    if($qe < $xb or $qb > $xe) { }
    elsif($qe > $xe and $qb > $xb) { 
      $qd= $xe - $qb; $qb=$xe+1;  
      $ok=0 if($qd > $qe - $xe); # overlap > extension
    } elsif( $qb < $xb and $qe < $xe) {
      $qd= $qe - $xb; $qe=$xb-1;
      $ok=0 if($qd > $xb - $qb); # overlap > extension
    } else {
      $qd= 1+ $qe - $qb; $qb=$qe=$xe; $ok=0;
    }
    if($se < $tb or $sb > $te) { }
    elsif($se > $te and $sb > $tb) { 
      $sd= $te - $sb; $sb=$te+1;      
      $ok=0 if($sd > $se - $te); # overlap > extension
    } elsif( $sb < $tb and $se < $te) {
      $sd= $se - $tb; $se=$tb-1;
      $ok=0 if($sd > $tb - $sb); # overlap > extension
    } else {
      $sd= 1+ $se - $sb; $sb=$se=$te; $ok=0;
    }
    ## fixme, ok=0 if overlap > extension..
    if($qd or $sd) {
      $ov++; $maxt+= ($qd>$sd)?$qd:$sd;
    }
    last unless($ok);
  }  
  push @ovs, [$qb,$qe,$sb,$se] if($ok);
  return($ok, $qb,$qe,$sb,$se,$maxt);
}

## method from makeblastscore2.pl .. not quite what is needed here.
# sub sumscore {
#   my( $q, $t, $bits,$aln,$aident, $qb,$qe,$sb,$se) = @_;
#   my $or=0;
#   my $ttrue=$t; # 2015.02 ?? best spans ignoring scaffold splits ?
#   if($qb > $qe) { ($qb,$qe)= ($qe,$qb); $or--; } # FIXME record $or in bspans
#   if($sb > $se) { ($sb,$se)= ($se,$sb); $or--; }
#   unless($bspans{$t}) { 
#     $bspans{$t}=[]; 
#     push( @{$bspans{$t}}, [$qb,$qe,$sb,$se,$bits,$aln,$aident,$or,$ttrue]); 
#     return; }
#   my $ov=0;
#   my $qlen=1+$qe-$qb; my $slen=1+$se-$sb;
#   my $qslop= _max($OVSLOP, int($pctOVSLOP*$qlen));
#   my $sslop= _max($OVSLOP, int($pctOVSLOP*$slen));
#   
#   foreach my $sp (@{$bspans{$t}}) {
#     #o my($xb,$xe,$tb,$te,$xbit,$orsp)= @$sp; # BUG here 201407: missing aln,aident!
#     my($xb,$xe,$tb,$te)= @$sp;  my $ttrue1= $$sp[8];
#     if($qe < $xb or $qb > $xe) { }
#     elsif($qe > $xe and $qb >= $xe - $qslop) { }
#     elsif($qb < $xb and $qe <= $xb + $qslop) { }
#     else { $ov=1; last; }
#     if($se < $tb or $sb > $te) { }
#     elsif($se > $te and $sb >= $te - $sslop) { }
#     elsif($sb < $tb and $se <= $tb + $sslop) { }
#     else { $ov=1; last; }
#   }  
#   unless($ov) { push( @{$bspans{$t}}, [$qb,$qe,$sb,$se,$bits,$aln,$aident,$or,$ttrue]); }
# }
