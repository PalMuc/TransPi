#!/usr/bin/env perl
# orthomcl_tabulate.pl

=item about

  collection of cmdlines to post-process orthomcl results to summary tables
  and annotated gene family document (ugp.xml) for web search/report
  
=item author

  from EvidentialGene collection
  D. Gilbert, gilbertd at indiana edu  2010..2012

=item updates to add

  see kfish2/prot/fish11c/omclwork.info   2013nov 
  
=cut


use FindBin;
use lib ("$FindBin::Bin/.."); # assume evigene/scripts/omcl/ << this path
my $EVIGENES=$ENV{EVIGENES} || "$FindBin::Bin/..";  

use strict;   
use Getopt::Long;
use constant VERSION => '2014.08.22'; # 03.22'; #'2013.12.30' 11.15'; # 10.08'; # 

## user options
my $DEBUG=(defined $ENV{debug}) ? $ENV{debug} : 1;
my $IDPRE= $ENV{idprefix} || "ARP9_G";
my $GTAG= $ENV{gtag} || ""; # make from IDPRE,dont need both but for history .."ARP";
my $mcl= $ENV{MCL} || "/bio/bio-grid/mb/mcl9/"; ## findapp('mclcm') and 'mcxdump'
## omcl_mcl2cluster -I $INFLATE is run option .. granularity of grouping (small = narrow groups; 1.5 default input omcl)
my $mcl2INFLATE = 3;
my $MINCOMMON= $ENV{mincommon} || 0; #  || $nspecies - 2;  # $SMIN
  
my $orun="Jan_25";      # get from dirlist
my $phyla="arp11u11";   # from dirlist
my($bpofile, $nspecies, $specieslist, $sppgenes, $omclpath, $namepath, $steps, $logfile);

##2017nov: drop default goodname/poorname
my $goodname=$ENV{goodname}||''; # was 'mayzebr|human|kfish2|platyfish'; 
my $poorname=$ENV{poorname}||''; # was 'stickleback|medaka|tetraodon';

# my $spl="antc,anth,aphid,apis2ref,bombusimp,bombusterr,daphnia,drosmel,human,trica,wasp"
# my $spl='ACEPH,AECHI,AMELL,CFLOR,HSALT,LHUMI,PBARB,SINVI,bombusimp,bombusterr,wasp'
#^ parse from $orun/parameter.log  SPECIES ..

my $optok= GetOptions( 
  "omclpath|datapath=s", \$omclpath, 
  "names|namepath=s", \$namepath, 
  "steps=s", \$steps, 
  "IDPREFIX=s", \$IDPRE, "GTAG=s", \$GTAG, # dont need both, make GTAG from IDPREFIX
  "MINCOMMON=i", \$MINCOMMON,  
  "mcl2INFLATE=i", \$mcl2INFLATE,  
##  "NCPU=i", \$NCPU, "MAXMEM=i", \$MAXMEM,  
##  "tidyup!", \$tidyup, 
##  "dryrun|n!", \$dryrun,
  "goodname=s", \$goodname, "poorname=s", \$poorname, # for arp_condef3.pl
  "logfile:s", \$logfile,
  "debug!", \$DEBUG,
);

die "EvidentialGene orthomcl_tabulate.pl VERSION ",VERSION,"
  convert orthomcl outputs to tables
Usage: ...
  opts: -idprefix=FISH11_G -omclpath=path/to/omclouts/ -namepath=path/to/*.names ..  -debug 
" unless($optok); 

if($IDPRE and not $GTAG) { ($GTAG=$IDPRE)=~s/_[^_]*$//; } # or keep same?
##  $goodname =~ s/'//g; $poorname =~ s/'//g; 

## separate/ pretab option: 
##  omcl_renumber(old_GTAG,old_omclgn.tab,newall_orthomcl.out) > newall_orthomcl.renum

sub MAINstub{}
MAIN:{
  chdir($omclpath) if($omclpath and -d $omclpath);
  omclout_startup();  # set orun, phyla from dirlist
  ($nspecies, $specieslist, $sppgenes)= omcl_getspecies();
  $MINCOMMON= $nspecies - 2 if($MINCOMMON==0);
  die "#ERROR: nspecies=$nspecies, list=@$specieslist\n" if($nspecies<2);
  
  omcl_gntab() unless($steps and $steps !~ /gntab/);
  omcl_pastebpo() unless($steps and $steps !~ /pastebpo/);
  omcl_count() unless($steps and $steps !~ /count/);
  omcl_gclass() unless($steps and $steps !~ /class/);
  omcl_gcommon() unless($steps and $steps !~ /common/);

  omcl_avematch(); # eval_orthogroup_genesets.pl, need input blastp.tall4 tables

  ## FIXME:  omcl_genegroupdoc:lookForAAlen($ARGV[0]) in:fish12a_omclgn2sum.tab > fish12a.aa.qual/.aa.count
  if($namepath) {
  omcl_consensusdef()  unless($steps and $steps !~ /condef/);  # maybe ok with -namepath ../names/; needs prepared names/*
  omcl_genegroupdoc() unless($steps and $steps !~ /genegroup/);  # maybe ok with -namepath; needs condef, names, ..
  }
  #notused# omcl_groupidtab(); # fixme

  omcl_subgroups() unless($steps and $steps !~ /subgroup/); # fixmed
  omcl_1to1tab() unless($steps and $steps !~ /1to1/);   # fixme opts
  omcl_orpartab(1) unless($steps and $steps !~ /orpar/);   # 1,0 : 2 tables? fixme
  omcl_orpartab(0) unless($steps and $steps !~ /orpar/);    
  
  # omcl_mcl2cluster();  # option or later
        
  omclout_finish() unless($steps and $steps !~ /finish/); # gzip big files
}

#-----------------------------------------
sub sysrun
{
  my @CMD= @_;
  if($DEBUG) { warn "#sysrun: ",join(" ",@CMD),"\n"; }
  return system(@CMD);
}

## input gene ids have species_ prefix, omcl adds species: in front; remove dupl prefix  spp:spp_ID

sub omclout_startup
{
  #above# chdir($omclpath) if($omclpath and -d $omclpath);
  opendir(D,"./"); my @fs= readdir(D); close(D);
  ($bpofile)= grep /\.bpo$/, @fs;
  $phyla= $bpofile;  $phyla =~ s/.bpo//; $phyla =~ s/_omcl$//; #?
  ##^^fixme handle bpo.gz ; need openread($fname) ...
  
  # ($orun) = grep /(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)/i, @fs; # better way? set outdir in omcl run?
  ##^^ this Date grep is bad, user option or other check needed; dont pick 1st match, try all for dirs
  ($orun)=  grep { (-d $_ && -f "$_/all_orthomcl.out") } @fs;
  
  # loggit(0,"input fileset=$bpofile,$orun/all_orthomcl.out");
  die "ERR: missing orthomcl speciesblast.bpo, Date/all_orthomcl.out fileset\n" # loggit(LOG_DIE,...)
    unless( -f $bpofile and -d $orun and -f "$orun/all_orthomcl.out");
  sysrun('perl', '-pi.old', '-e', 's/ \w+:/ /g;', $orun.'/all_orthomcl.out')
    unless($ENV{'orthoutisokay'} or $steps); # cut dupl spp:spp_geneID
}

sub omclout_finish
{
  sysrun("echo rm ${phyla}_omcl_bpo.se ${phyla}_omcl_bpo.idx  $orun/all_orthomcl.out.old");
  sysrun("gzip --fast $orun/all_blast.bbh $orun/tmp/all_ortho.*"); # NOTnow:  $orun/all_orthomcl.out 
  (my $ggfile=$bpofile)=~s/.bpo/.gg/; 
  sysrun("gzip --fast $bpofile $ggfile ${phyla}_omclgns2.tab");
}
  
sub omcl_getspecies
{
  my @spp=(); my %sppgenes; my $at;
  open(I,"$orun/parameter.log") or die "open $orun/parameter.log";
  while(<I>) {
    if(/^#+SPECIES/) { $at=1; }
    elsif($at==1 and /^#/) { $at= 0;} #  last; 
    elsif($at==1 and /^\s+(\w+)\s+(\d+) genes/) { push @spp, $1; $sppgenes{$1}=$2; }
    # if(/^#+PARAMETERS/) { $at="param"; }
    # if(/^#+FILES/) { $at="files"; }
    # if(/^#+START TIME/) { $at="time0"; }
    # if(/^#+END TIME/) { $at="time1"; }
  } close(I);
  my $nspp= @spp;
  warn "#info: nspecies=$nspp, list=@spp\n" if($DEBUG);
  return ($nspp, \@spp, \%sppgenes);
}


sub  omcl_gntab
{
  open(I,"$orun/all_orthomcl.out") or die "open $orun/all_orthomcl.out";
  open(O,">${phyla}_omclgn.tab");
  my( %alts, %didgn);
  while(<I>) { 
    if(/^ORTHOMCL/) { 
    my($om,$gn)=split /:\s+/,$_,2;  my @gn= split" ",$gn; 
    my($og,$ng,$nt)= $om =~ m/ORTHOMCL(\d+).(\d+) genes,\s*(\d+) taxa/;  
    my %didg=(); 
    foreach my $gn (@gn) { 
      my($sp)=split"_",$gn,2; (my $g1=$gn)=~s/[pt\.]\d+$//;  
      $didgn{$g1}++; 
      # $isaltNOT=($didgn{$g1}>1)?1:0; $alts{$sp}++ if($isalt);
      print O "$GTAG$og\t$gn\n"; # unless($isalt); 
      }
    }
  }
  # if(%alts){ print O "#alttr-drops: "; foreach $s (sort keys %alts){ print O "$s=$alts{$s}, "; } print O "\n"; } 
  close(O); close(I); 
}

sub alnspan
{
  my($spans,$qlen,$slen)=@_;
  #qlen=875;slen=1252;spans=1:3-366:378-735.2:505-875:838-1252.3:366-494:46-142.
  my($qa,$sa)=(0,0);
  for my $s (split/\./,$spans) {
    my($i,$qs,$ss)=split":",$s; 
    my($qa1,$sa1)= map{ my($bw,$ew)=split/\-/; 1+$ew-$bw; } ($qs,$ss);
    $qa += $qa1; $sa += $sa1;
    # my($qb,$qe,$sb,$se)= map{ split/\-/ } ($qs,$ss); $qa += 1+$qe-$qb; $sa += 1+$se-$sb;
  }
  my $qap = int(1000*$qa/$qlen)/10;
  my $sap = int(1000*$sa/$slen)/10;
  my $ap= int(1000*($qa+$sa)/($qlen+$slen))/10; # int(($qap+$sap)/2);
  $ap=100 if($ap>100);  $ap=~s/\.0$//;
  return (wantarray)?($ap,$qap,$sap,$qa,$sa):$ap;
}

sub omcl_pastebpo
{
  # sysrun("cat ${phyla}_omclgn.tab ${phyla}*_omcl.bpo | env idprefix=$GTAG $EVIGENES/omcl/pastebpo.pl > ${phyla}_omclgns2.tab");
  # omcl/pastebpo.pl
  my $tag=$GTAG; # $ENV{idprefix} || "ARP";
  my( %gog, %ogtax, $nout,);

  ## ? add group->ntaxa hash also here? for inpar vs upar
  open(I,"${phyla}_omclgn.tab") or die "open ${phyla}_omclgn.tab";
  while(<I>){ if(/^($tag\d+)\t(\S+)/) { my($og,$gn)=($1,$2); my($sp)=split"_",$gn; 
    $gog{$gn}=$og; $ogtax{$og}{$sp}++; } 
  } close(I);

  open(I,$bpofile) or die "open bpofile:$bpofile";
  open(O,">${phyla}_omclgns2.tab");
  
  my($lg1,$sp1,$ordid,$pardid,$orpar,$hasorth);
  while(<I>){
    if(/^\d+;/) { # xxx.bpo
      my @v= split";"; 

    ## FIXME2: add align scores, from simspan, pct-align = int(100*align/min(qlen,slen)) : min or max or ave len?
    ## FIXME: revise here to calc inparalogs, outparalogs. expect bpo-query x subj ordered by besthit
    ##  inparalogs = samespp-subj before otherspp
    ## bpoline= "$qid;$qlen;$sid;$slen;$prob;$percentIdent;$simspan"; 
    ## simspan is list of pairalign spans: 1:q1-q2:s1-s2.2:q1-q2:s1-s2.3:...
    #19102;kfish2_Funhe2EKm017121t1;875;mayzebr_XP_004572373.1;1252;0e+00;64;1:3-366:378-735.2:505-875:838-1252.3:366-494:46-142.
    #19103;kfish2_Funhe2EKm017121t1;875;tilapia_ENSONIP00000018723;1041;0e+00;69;1:1-367:313-676.2:517-875:674-1041.3:418-496:2-77.
use constant NEWGNTAB => 1;
if(NEWGNTAB) {
      my($si,$g1,$g1w,$g2,$g2w,$ev,$pi,$spans)= @v;
      next if($g1 eq $g2);
      my $og1= $gog{$g1};
      my $og2= $gog{$g2};
      if($g1 ne $lg1) { 
        $ordid=$pardid=$orpar=0; ($sp1)=split"_",$g1,2; 
        $hasorth=(1 < scalar( keys %{$ogtax{$og1}} ))?1:0;
      }
      next unless($og1 and $og1 eq $og2);
      my($sp2)=split"_",$g2,2;  
      ## NOTE: inpar here includes UDup, unique-species groups, should use other class for that: upar
      my $samespp= ($sp1 eq $sp2)?1:0;
      if($samespp) { 
        my $clpar=($hasorth)?"inpar":"upar";
        $pardid++; $orpar=($ordid)?"opar$pardid":"$clpar$pardid"; 
      } else { $ordid++; $orpar="orlog$ordid"; }
      #^^FIXME: tie-name-sort effect here for orlog ordid numbering. check ev,pi,gw1/2 for ties?
      
      my $alnave= alnspan($spans,$g1w,$g2w); # pct align for each, or min/max/ave ?
      print O join("\t",$g1,$g2,$og1,$ev,$pi,$orpar,$alnave),"\n";
      $lg1=$g1; $nout++;
} else {
      my($g1,$g2,$ev,$pi)= @v[1,3,5,6];
      next if($g1 eq $g2);
      my $og1= $gog{$g1};
      my $og2= $gog{$g2};
      next unless($og1 and $og1 eq $og2);
      print O join("\t",$g1,$g2,$og1,$ev,$pi),"\n";
      $nout++;
}
      }
  }
  # warn "#pastebpo n=$nout\n" if($DEBUG);
  close(I);

  open(I,"${phyla}_omclgns2.tab") or die "open ${phyla}_omclgns2.tab"; 
  open(O,">${phyla}_omclgn2sum.tab");
  #FIX2? add align ave to sum.tab ?
  
  my($g1,$g2,$og,$ev,$pi);
  my($lg,$log,$ng,$sev,$spi,$mpi,$mg,$mev)= (0) x 10;
  sub dumpg{ if($ng>0){
    my $aev=sprintf"%.4g",$sev/$ng; my $api=int($spi/$ng); 
    print O join("\t",$lg,$log,$ng,$aev,$api,$mg,$mev,$mpi),"\n"; } 
    $mpi=$mg=$mev=$lg=$ng=$sev=$spi=0; 
    } 
    
  while(<I>) {
    ($g1,$g2,$og,$ev,$pi)=split; dumpg() unless($g1 eq $lg); 
    $lg=$g1; $log=$og; $ng++; $sev+= $ev; $spi += $pi; 
    if($pi>$mpi){ $mpi=$pi; $mg=$g2; $mev=$ev; } 
    }
  dumpg();
  close(O); close(I);
  warn "#info: pastebpo=$nout\n" if($DEBUG);
}


sub omcl_count
{  
  open(I,"$orun/all_orthomcl.out") or die "open $orun/all_orthomcl.out";
  open(O,">${phyla}-orthomcl-count.tab");

  my($sog,$snt,$sng)=(0) x 3;
  my @spl= @$specieslist;  # $spl=$ENV{spl}; @spl=split",",$spl;
  my %drop; @drop{@spl}= (0) x $nspecies; 
  print O join("\t","OID","Nt","Ng",@spl),"\n";
  while(<I>) {  
    if(/^ORTHOMCL/) { 
    my ($om,$gn)=split /:\s+/,$_,2; my @gn= split" ",$gn; 
    my %spc; @spc{@spl}= (0) x $nspecies; 
    my %didgn=(); 
    foreach my $sg (@gn) { 
      (my $g1=$sg)=~s/[pt\.]\d+$//; # $isaltNOT=($didgn{$g1}++ > 0)?1:0; 
      my($sp,$gn)=split"_",$sg,2; $spc{$sp}++; ## unless($isalt); 
    }
    my($og,$g,$t)= $om =~ m/ORTHOMCL(\d+).(\d+) genes,\s*(\d+) taxa/; 
    my @sc= @spc{@spl}; print O join("\t",$og,$t,$g,@sc),"\n"; 
    $sog++; $sng+=$g; $snt+=$t;
    } 
  }
  close(O); close(I);
  
  $snt= int(100*$snt/$sog)/100; $sng=int(100*$sng/$sog)/100;
  warn "#info: ${phyla}-orthomcl-count.tab  ngr=$sog, taxa/gr=$snt, spp/gr=$sng;\n" if($DEBUG);
}

=item omcl_gclass

# summary gene class table
#? add spp correlations in genes/group ? vxy, vxx, vyy table == distance mat
#? add OrGrp count: number ortho groups species belongs; Orth1 + OrDup-grps
# iskip=8 == human  DANG NO, index from col=0 OID, human skip=11

=cut

sub omcl_gclass
{  
  open(I,"${phyla}-orthomcl-count.tab") or die "open orthomcl-count.tab";
  open(O,">${phyla}-orthomcl-gclass.tab");

  ## FIXed: use %sppgenes to fill in Uniq1, inGene
  #cat ${phyla}-orthomcl-count.tab | env iskip=-1 perl -ne
  my $iskip= -1; # $ENV{iskip};
  my $miss= 1; # $ENV{miss};
  my $onlys= 0; # $ENV{onlys}; # onlyspecies

  my(@hd, @spp, $MIS1);
  my @osp=(); if($onlys) { @osp= split /[,\s]+/,$onlys; }

  print O join("\t",
    qw(species inGene oGene Uniq1 UDup Orth1 OrDup nGroup UniqGrp OrGrp OrMis1 Guniq Gmax Gmin)),"\n";

  my(%sc,%sg,%sgu,%sgh1,%sghp,%sgOrG,%sgMIS,%suni,%smax,%smin);
  while(<I>) {
    my @v=split; my $nv=$#v; 
    if(/^OID/){ @hd=@v; @spp=@v[3..$nv]; $MIS1= @spp - $miss; next; } 
    elsif(/^\D/){ next; } #? ok
    
    my $nt=0; 
    #count instead# my $ntbad=$v[1];  my $ngbad=$v[2]; # fixme: ntbad : count
    # my @s=map{ $hd[$_]."=".$v[$_]; } (3..$nv); # fixme
    my %csp=(); my @s= map{ my $c=$v[$_]; $nt++ if($c>0); $csp{ $hd[$_] }=$c; $hd[$_]."=".$c; } (3..$nv);
    if(@osp) { my $ok=0; map{ $ok++ if( $csp{$_}) } @osp; next unless($ok); }

    my($vm,$smax,$vm2,$smin)= (0) x 10; my $vn=99999; 
    map{ my($s,$v)=split"="; 
      $sc{$s}++ if($v>0); $sg{$s} += $v; 
      $sgu{$s}  += $v if($nt==1); 
      $sgh1{$s} += $v if($nt>1 and $v==1); 
      $sghp{$s} += $v if($nt>1 and $v>1); 
      $sgOrG{$s}++ if($nt>1 and $v>=1); 
      $sgMIS{$s}++ if($v==0 and ($nt>=$MIS1 or ($iskip>=0 and $nt>=$MIS1-1 and $v[$iskip]==0))); 
      if($v<=$vn){ $smin=($v==$vn)?"$smin,$s":$s; $vn=$v;} 
      if($v>=$vm){ $smax=($v==$vm)?"$smax,$s":$s; $vm2=$vm; $vm=$v;}  
    } @s; 
    $suni{$smax}++ if($nt==1); 
    $smax{$smax}++ unless($nt < 2 or $smax =~ /,/);  
    $smin{$smin}++ unless($nt < 4 or $smin =~ /,/);  
  }
  foreach my $s (sort keys %sc) { 
    my($ing,$ong,$u,$h,$d,$c,$ngu,$ngo1,$ngo2,$ngo0,$nog,$ug1)= (0) x 20;
    $ing= $sppgenes->{$s} || 0; 
    $u=$suni{$s}||0; $h=$smax{$s}||0; $d=$smin{$s}||0; $c=$sc{$s}||0; $ong=$sg{$s}; 
    $ngu=$sgu{$s}||0; $ngo1= $sgh1{$s}||0; $ngo2= $sghp{$s}||0; $ngo0=$sgMIS{$s}||0; $nog= $sgOrG{$s}||0;  
    $ug1= ($ing>=$ong) ? $ing - $ong : "na";
    print O join("\t",$s, $ing, $ong, $ug1, $ngu,$ngo1,$ngo2,  $c,$c-$nog,$nog,$ngo0,  $u,$h,$d),"\n"; 
  }      
  close(O); close(I);    
  warn "#info: ${phyla}-orthomcl-gclass.tab;\n" if($DEBUG);
}
 

=item omcl_gclass reformat table this way

Fish orthology gene groups summary (OrthoMCL).
            ---------- GENES -----------    -------- GROUPS --------------  
            nGene Orth1 OrDup Uniq1 UDup    nGroup  OrGrp OrMis1 UniqGrp
            ----------------------------    -----------------------------  
killifish   30000 13958  9397* 2837 3808    17448   16587    107    861
zebrafish   29004 10053 12973  3730 2248    15081   14513    284    568
tilapia     21442 12656  7463   858  465    15207   15086    201    121
stickleback 20875 12669  6327  1560  319    14866   14797    264     69
medaka      19732 12130  5456  1660  486    14145   14015    538    130
tetraodon   19646 11624  5544  2263  215    13725   13649    523     76
-----------------------------------------------------------------------
  Uniq1,UDup  = single-copy and duplicated species-unique genes
  Orth1,OrDup = single-copy and duplicated orthologous genes
  UniqGrp,OrGrp = species-unique and orthologous groups
  OrMis1  = groups missing in species that all other species have

=item omcl_gcommon cmdline

# sppindex="0,1,3,4,5,6,7,8,9,10,11" == fish, human=2; min=10
# fix for skipset .. min of others.
cat *-orthomcl-count.tab | env min=11 skip=3 perl -ne \
'BEGIN{ $MINC=$ENV{min}||2; $SKIP=$ENV{skip}||0; }\
if(/^OID/){ @hd=split; next; } next if(/^\D/); my @v=split; my $nv=$#v; my $nt=$v[1]; \
$nt++ if($SKIP and $v[$SKIP-1]==0); next if($nt<$MINC); \
$ng++; for my $i (3..$nv) { my $c=$v[$i]; my $s=$hd[$i]; \
my $have=($c>0)?"have":"miss"; $have{$s}{$have}++; } \
END{  print "Common gene families presence, ncommon=$ng, min taxa=$MINC\n"; \
  printf "%12s\t%5s\t%5s\n",qw(Species Have Miss); \
  foreach my $s (sort{ $have{$a}{miss} <=> $have{$b}{miss} or $a cmp $b } keys %have)\
  { my($h,$m)= @{$have{$s}}{qw(have miss)}; printf "%12s\t%5d\t%5d\n",$s,$h,$m; }\
}'\


=cut

sub omcl_gcommon
{
  open(I,"${phyla}-orthomcl-count.tab") or die "open orthomcl-count.tab";
  open(O,">${phyla}-orthomcl-gcommon.tab");
  # cat ${phyla}-orthomcl-count.tab | env min=6 perl -ne 
  # my $SMIN= $ENV{mincommon} || $nspecies - 2;
  my(%have, @hd, $ng);
  while(<I>) {
    if(/^OID/){ @hd=split;  next; }    
    elsif(/^\D/){ next } 
    my @v=split; my $nv=$#v; my $nt=$v[1]; 
    if($nt>=$MINCOMMON) {  $ng++; 
      for my $i (3..$nv) { my $c=$v[$i]; my $s=$hd[$i]; my $have=($c>0)?"have":"miss"; $have{$s}{$have}++; } 
      } 
  }
  print O "Common gene families presence, ncommon=$ng, min taxa=$MINCOMMON\n"; 
  print O join("\t",qw(Species Have Miss))."\n";
  foreach my $s (sort{ $have{$a}{miss} <=> $have{$b}{miss} or $a cmp $b } keys %have) {  
    my($h,$m)= @{$have{$s}}{qw(have miss)}; 
    print O join("\t",$s,$h,$m)."\n"; 
  }
  close(O); close(I);
  warn "#info: ${phyla}-orthomcl-gcommon.tab;\n" if($DEBUG);
}

sub omcl_avematch
{
  my $info= 
"#info omcl_avematch not ready, would run this but need blast to species.tall tables: 
#i  for eachspecies specieslist; do 
#i    env aa=${phyla}_omcl.aa.qual tall=1 $EVIGENES/makeblastscore2.pl \\
#i     ${phyla}-${phyla}.blastp.gz > ${phyla}_eachspecies.tall4 
#i
#i    $EVIGENES/eval_orthogroup_genesets.pl -nodebug -bitmed -mintaxa 3 \\
#i    -out ${phyla}-eachspecies.topout -tallscore ${phyla}-_eachspecies.tall4 \\
#i    -groupgene ${phyla}_omclgn.tab -groupcount ${phyla}-orthomcl-count.tab \\
#i    -geneaa ${phyla}_omcl.aa.qual
#i  done
#i  cat ${phyla}-*.topout | $EVIGENES/eval_orthogroup_genesets.pl -intab stdin 
#i
#\n";
  warn $info if($DEBUG);
}

=item omcl_avematch

  set karp=$kfish/prot/kfish1ball/omclkf2
  set spp=human
  
  # $EVIGENES/makeblastscore2.pl fish8-xxx.blastp.gz > fish8-all.tall4
  env aa1=fish10main.aa.qual keepho="human:" tall=1 $EVIGENES/makeblastscore2.pl \
    fish10main-fish10main.blastp.gz > fish10main_human.tall4 
  
  .. split fish8-all.tall4 to fish8-{eachspecies}.tall4
  .. OR split blastp by '^species:', run makeblast tall4 on each
  
  $EVIGENES/eval_orthogroup_genesets.pl -nodebug -bitmed -mintaxa 3 \
  -out fish8-$spp.topout2 -tallscore fish8-$spp.tall4 \
  -groupgene $karp/fish8_omclgn.tab -groupcount $karp/fish8-orthomcl-count.tab\
  -geneaa $karp/fish8_omcl.aa.gz 
  
  #........... omclkf2/Aug_25/  : removed alts
  -- slight changes only from alt removal
  -- adding best of kfish alts/main doesnt change these scores; check alt for missed omcl groups
  cat fish8-{kill,medaka,stick,tetra,tilapia,zebr}*.topout2  | $EVIGENES/eval_orthogroup_genesets.pl -intab stdin -digit 0
  
   ...  Common Fish groups ... 
  Geneset  : nGene    Bits    dSize   rGene   rBits   outlierSize
  killifish: 10172    718     9       10172   718     x2Tiny=244 (2.3%)   x2Big=92 (0.9%)
  medaka   : 10172    699     -29     10172   699     x2Tiny=600 (5.8%)   x2Big=20 (0.1%)
  sticklebk: 10172    721     -19     10172   721     x2Tiny=387 (3.8%)   x2Big=22 (0.2%)
  tetraodon: 10172    688     -43     10172   688     x2Tiny=576 (5.6%)   x2Big=19 (0.1%)
  tilapia  : 10172    753     19      10172   753     x2Tiny=60 (0.5%)    x2Big=67 (0.6%)
  zebrafish: 10172    685     11      10172   685     x2Tiny=195 (1.9%)   x2Big=105 (1%)
   ...  All groups ...
  Geneset  : nGene    Bits    dSize   rGene   rBits   outlierSize
  killifish: 16319    574     8       15547   602     x2Tiny=576 (3.7%)   x2Big=311 (2%)
  medaka   : 16319    525     -34     13909   616     x2Tiny=1058 (7.6%)  x2Big=44 (0.3%)
  sticklebk: 16319    560     -26     14624   625     x2Tiny=832 (5.6%)   x2Big=53 (0.3%)
  tetraodon: 16319    515     -42     13580   619     x2Tiny=867 (6.3%)   x2Big=56 (0.4%)
  tilapia  : 16319    597     18      14896   654     x2Tiny=146 (0.9%)   x2Big=167 (1.1%)
  zebrafish: 16319    524     18      14634   584     x2Tiny=311 (2.1%)   x2Big=318 (2.1%)
  human    : 16319    441     45      12752   564     x2Tiny=68 (0.5%)    x2Big=280 (2.1%)
  xenopus  : 16319    429     0       12595   556     x2Tiny=261 (2%)     x2Big=131 (1%)

=item avematch reformat table this way

  Fish species average match to orthology gene groups.
           Common groups(n=10172)        All groups 
  Geneset    cBits  dSize     tBits  rGene   Small outliers
  --------- --------------------------------------
  killifish  718     9         574   15547   3.7%
  medaka     699     -29       525   13909   7.6%
  sticklebk  721     -19       560   14624   5.6% 
  tetraodon  688     -43       515   13580   6.3% 
  tilapia    753     19        597   14896   0.9% 
  zebrafish  685     11        524   14634   2.1% 
  ------------------------------------------------
    cBits = bitscore average to common groups
    tBits = bitscore average to all groups
    dSize = average size difference from group median
    Small outliers : percent species genes < 2sd of median gene size in group 
    rGene = number of gene groups found in species
   
=cut



sub getnames {
  return () unless($namepath and -d $namepath);
  opendir(D,$namepath); my @fs= readdir(D); close(D);
  my @names= map{ s,^,$namepath/,; $_ } grep /\.names/, @fs; # global? reuse..
  return (wantarray)? @names : join(" ",@names);
}

sub omcl_consensusdef
{
  # my $names= getnames(); 
  my @names= getnames();
  my $outfile="${phyla}_omclgn.consensus_def.txt";
  do { warn "ERR:omcl_consensusdef no names in $namepath"; return -1; } unless(@names>0);
  
  my $inlist = join(" ",@names,"${phyla}_omclgn.tab"); # "$names ${phyla}_omclgn.tab";
  #above# my $goodname=$ENV{goodname}||'mayzebr|human|kfish2|platyfish';
  #above# my $poorname=$ENV{poorname}||'stickleback|medaka|tetraodon';
  $goodname =~ s/'//g; $poorname =~ s/'//g; 

  # in: open(I,"cat $names ${phyla}_omclgn.tab |") or die "pipe $names ${phyla}_omclgn.tab";
  # out: > ${phyla}_omclgn.consensus_def.txt

  my $cmd="$EVIGENES/omcl/arp_condef3.pl -noput -nodigit -nolike -gtag=$GTAG -idprefix=$IDPRE"
    . " -good='$goodname' -poor='$poorname' $inlist > $outfile";
    
  my $err= sysrun($cmd);
  return ($err);
}

=item omcl_consensusdef

 ## use revised arp_condef3.pl; add weights to good spp names?
# note: mayz,hum,kfish use UniProt std names, better than ensembl/zfin names in other fish
  
  pt=11F
  cat ../names/*.names ${phyla}_omclgn.tab | \
  $evigene/scripts/omcl/arp_condef3.pl -noput -nodigit -nolike -gtag="FISH$pt" -idprefix="FISH${pt}_G" \
  -good='mayzebr|human|kfish2|platyfish' -poor='stickleback|medaka|tetraodon' \
  > ${phyla}_omclgn.consensus_def.txt

 
  cat ${phyla}_omclgn.consensus_def.txt | env recase=1 count=withID debug=1 \
  $EVIGENES/bestgenes_puban_wasp.pl | sed 's/^$IDPRE//' | sort -k1,1n | sed 's/^/$IDPRE/; s/TE:TE:/TE:/' 
    > ${phyla}_omclgn.consensus_def.rename.txt
  mv ${phyla}_omclgn.consensus_def.txt  ${phyla}_omclgn.consensus_def.txt0
  mv ${phyla}_omclgn.consensus_def.rename.txt ${phyla}_omclgn.consensus_def.txt 

=item pull names

  orthodb6.aa:  gzgrep '^>' $pt.aa.gz | perl -ne\
  '($id,$pid,$uid,$na)=m/>(\S+)\s(\S+)\s+(\w*)\s*(.*)$/; $na=~s/IPR.*$//; print "$id\t$na\n" if($na =~ /[a-z]/);' \
  > ../names/$pt.names
  
  evigene.aa: gzgrep '^>' $pt.aa.gz | perl -ne\
  '($id)=m/^>(\S+)/; ($na)=m/Name=([^;\n]+)/; print "$id\t$na\n" if($na =~ /\w/);' \
  > ../names/$pt.names
  
  refseq.aa:
     
=cut


sub omcl_genegroupdoc
{
  ## fixme: GetOptions
  my $doxml=$ENV{xml}||0; 
  my $myspecies= $ENV{spt}||$ENV{myspecies}||"FIXME_species"; # fixme
  my $speciesmap= $ENV{speciesmap}||"";  
  my $date= $ENV{date}||"20140000"; # fixme
  my $clade=$ENV{clade}||"CladeName";
  ## my $title=$ENV{title}||"Title here";
  my $names= getnames(); 
  my $inlist="${phyla}_omclgn2sum.tab  ${phyla}_omclgn.consensus_def.txt $names $orun/all_orthomcl.out ";
  my $outfile="${phyla}_genes.ugp.".($doxml?"xml":"txt");
  $myspecies =~ s/,/|/g;

  ## FIXME:  genegroupbpo.pl: lookForAAlen($inlist[0]) in: fish12a_omclgn2sum.tab > fish12a.aa.qual/.aa.count
  ## FIXME2: want GetOpts opts for -myspecies -date -clade etc and help info
  
  ## title from clade ..  title='$title'
  $speciesmap="speciesmap=$speciesmap" if($speciesmap);
  $date="date=$date" if($date);
  my $cmd="env xml=$doxml $speciesmap $date clade=$clade gtag=$GTAG idprefix=$IDPRE "
    ."$EVIGENES/omcl/genegroupbpo.pl $inlist > $outfile";
  ## FIXME: need other fish genegroupbpo.pl or many options ..
  
  my $err= sysrun($cmd);
  
  ## brief.txt failed.. no myspecies IDs
  unless($err or $doxml) { # require myspecies
    # drop egrep, do in perl loop ..
    # my $inf= "cat $outfile | egrep '(GeneID|ntaxa|ngene|occurrence|description): |$myspecies' |";
    my $inf=$outfile;
    my $outf="${phyla}_genes.ugp_brief.txt";
    open(I,$inf) or die "ERR: read $inf";
    open(O,'>',$outf) or die "ERR: write $outf";
    my($id,$d,$s);
    while(<I>) {
      next unless(m/(GeneID|ntaxa|ngene|occurrence|description): |$myspecies/);
      s/^ +//; 
      if(m/similarity:/ and /$myspecies/){ 
        ($id)=m/iden: (\d+)/; ($d)=m/acc: (\w[^;\s]+)/; $_="$d/$id,"; s/^/$myspecies: / if $s; $s=0;
      } else { 
        $s=1; print O "\n" if s/\s*GeneID:\s//; s/\n/ /;
      } 
      print O $_;
      }
    print O "\n"; 
    close(O); close(I);
  } 
  
  return ($err);
}

=item omcl_genegroupdoc

  cat \
  ${phyla}_omclgn2sum.tab  \
  ${phyla}_omclgn.consensus_def.txt \
  names/*.names  \
  $orun/all_orthomcl.out | \
  env xml=1 date=20120125 clade=HymenoInsectArp title='Arthropod gene group' gtag=$GTAG idprefix=$IDPRE \
  $EVIGENES/omcl/genegroupbpo.pl  
  > ${phyla}_genes.ugp.xml

=item genebriefdoc

	set spt=Nasvi
	set spt=Acyrthosiphon
  cat ${phyla}_genes.ugp.txt | egrep "GeneID|ntaxa|ngene|occur|descript|$spt" | env spt=$spt perl -pe \
  's/^ +//; if(m/similarity:/ and /$spt/){ ($id)=m/iden: (\d+)/; ($d)=m/acc: (\w[^;\s]+)/; $_="$d/$id,"; s/^/$spt: / if $s; $s=0;} \
  else{ $s=1; print "\n" if s/\s*GeneID:\s//; s/\n/ /;} BEGIN{$spt=$ENV{spt};} END{print"\n";} ' \
  >  ${phyla}_genes.ugp_brief.txt

=cut



sub omcl_groupidtab
{
}
=item omcl_groupidtab 

  # do this not for 11 taxa (all) but drop human, daphnia?
  cat ${phyla}-orthomcl-count.tab | env nt=11 gtag=$GTAG perl -ne
  '($od,$nt,$ng,@c)=split; print "$ENV{gtag}$od\t\n" if($nt==$ENV{nt});' 
	  > ${phyla}_omcl11.allgr.oids

  ggrep -F -f ${phyla}_omcl11.allgr.oids ${phyla}_omclgn.tab | cut -f2 | sed 's/^/gid /' > ${phyla}_omcl11.allgr.gids
  # use all.gids for restricted distance matrix restricted to common genes.

=cut

sub omcl_subgroups
{
  local(*SO);
  our( %sg, );
  my( %okoid,$ngn,$lorid,%inp,);
  my $inpipe="sort -k3,3 -k1,1 -k2,2 ${phyla}_omclgns2.tab| cat ${phyla}-orthomcl-count.tab - |";
  my $outfile="${phyla}_omsubgrp1.list";
    
  sub putorg {  our(%sg, );
  my $nig=0; my %ig=(); my @ga=sort keys %sg; 
  foreach my $ga (@ga) { my @gb= sort keys %{$sg{$ga}}; my $ig=0; my @og=(); 
    for my $g ($ga,@gb) { if(my $i=$ig{$g}) { push @og, $i; $ig=$i if($ig==0 or $i<$ig); }} 
    if($ig){ for my $j (@og){ next if($j == $ig); for my $k (keys %ig){ $ig{$k}=$ig if($ig{$k}==$j); }}}  
    else { $ig= ++$nig; }  
    for my $g ($ga,@gb) { $ig{$g}= $ig; } 
  } 
  my $lig=0;  my $jg=0; my @sg=sort{$ig{$a}<=>$ig{$b} or $a cmp $b} keys %ig; 
  for my $g (@sg) { my $ig=$ig{$g};
    unless($ig == $lig){ $jg++; print SO "\n$lorid.s$jg:\t"; } 
    print SO "$g,"; $lig=$ig; }  
  print SO "\n"; 
  }
  
  open(I,$inpipe) or die "ERR: pipe $inpipe";
  open(SO,'>',$outfile) or die "ERR: write $outfile";
  while(<I>) {
  if(/^\d/) { my($oid,$nt,$ng,@c)=split; my $sd=0; for my $c (@c) { $sd++ if($c>1); } 
    $okoid{$oid}=$sd if($nt>4 and $sd > $nt/2); } 
  elsif(/\t$GTAG/) { 
    my($ga,$gb,$orid,$ev,$pi,$orpar)=split; $ngn++; 
    my($oid)= $orid=~m/(\d+)$/; next unless($okoid{$oid});
    if($lorid and $orid ne $lorid) { putorg(); %sg=(); %inp=();} 
    if($orpar =~ m/orlog[1]$/) { $sg{$ga}{$gb}=$pi; $sg{$gb}{$ga}=$pi;} 
    elsif($orpar=~/inpar/) { $inp{$ga}{$gb}=$pi; $inp{$gb}{$ga}=$pi;} 
    $lorid=$orid; } 
  }
  putorg();  
  close(SO); close(I);
  
  ## do nametab also
  if(my @names= getnames()) {
  my $inlist = join(" ",@names,$outfile); # "$names ${phyla}_omsubgrp1.list";
  my $outnames="${phyla}_omsubgrp1.names";
  #above# my $goodname=$ENV{goodname}||'mayzebr|human|kfish2|platyfish'; #2017nov: drop default goodname/poorname
  #above# my $poorname=$ENV{poorname}||'stickleback|medaka|tetraodon';
  $goodname =~ s/'//g; $poorname =~ s/'//g;
 
  # -nonodigit... keep digits here ; 
  #? FIXME - no, leave: cut off src= from names. in arp_condef3.
  my $cmd="$EVIGENES/omcl/arp_condef3.pl -noput -nolike -gtag=$GTAG -idprefix=$IDPRE"
    . " -good='$goodname' -poor='$poorname' $inlist > $outnames";
    
  my $err= sysrun($cmd);
  }
  
}

=item omcl_subgroups

 sgroup choices: 
   1level,more: if($orpar =~ m/orlog1$/)
   2level,fewer: if($orpar =~ m/orlog[12]$/)

GTAG=FISH11D
IDPRE=FISH11D_G
 
cat ${phyla}_omclgns2.tab | sort -k3,3 -k1,1 -k2,2 | cat ${phyla}-orthomcl-count.tab - | \
env gtag=$GTAG perl -ne 'BEGIN{$GTAG=$ENV{gtag};} \
if(/^\d/) { ($oid,$nt,$ng,@c)=split; $sd=0; for $c (@c) { $sd++ if($c>1); } \
$okoid{$oid}=$sd if($nt>4 and $sd > $nt/2); } \
elsif(/\t$GTAG/) { \
($ga,$gb,$orid,$ev,$pi,$orpar)=split; $ngn++; \
($oid)= $orid=~m/(\d+)$/; next unless($okoid{$oid});\
if($lorid and $orid ne $lorid) { putorg(); %sg=(); %inp=();} \
if($orpar =~ m/orlog[1]$/) { $sg{$ga}{$gb}=$pi; $sg{$gb}{$ga}=$pi;} \
elsif($orpar=~/inpar/) { $inp{$ga}{$gb}=$pi; $inp{$gb}{$ga}=$pi;} \
$lorid=$orid;} END{ putorg(); } \
sub putorg {  $nig=0; %ig=(); @ga=sort keys %sg; \
foreach $ga (@ga) { @gb= sort keys %{$sg{$ga}}; $ig=0; @og=(); \
  for $g ($ga,@gb) { if($i=$ig{$g}) { push @og, $i; $ig=$i if($ig==0 or $i<$ig); }} \
  if($ig){ for $j (@og){ next if($j == $ig); for $k (keys %ig){ $ig{$k}=$ig if($ig{$k}==$j); }}} \
  else { $ig= ++$nig; } \
  for $g ($ga,@gb) { $ig{$g}= $ig; } } \
$lig=0; $jg=0; @sg=sort{$ig{$a}<=>$ig{$b} or $a cmp $b} keys %ig; \
for $g (@sg) { $ig=$ig{$g}; unless($ig == $lig){ $jg++; print "\n$lorid.s$jg:\t"; } print "$g,"; $lig=$ig; } \
print "\n"; }' > ${phyla}_omsubgrp1.list

cat ../names/*.names ${phyla}_omsubgrp1.list | \
env good='human|kfish|mayzebr|platyfish' poor='stickle|medaka|tetraodon' gtag=$GTAG idprefix=$IDPRE \
noclean=0 nolike=1 $evigene/scripts/omcl/arp_condef2.pl \
> ${phyla}_omsubgrp1.names

=item inparalog tab

grep inpar ${phyla}_omclgns2.tab | grep kfish | perl -ne\
'($ga,$gb,$oid,$ev,$pi,$cla)=split; print if($pi>=95);' > ${phyla}_kfish_inpar95.tab

=item orinpar count matrix

# ~/Desktop/genowork/kfish/fish11domcl.orpar.count
# or-inpar counts: inpar of orgrp only

cat ${phyla}_omclgns2.tab | perl -ne\
'($ga,$gb,$orid,$ev,$pi,$orpar)=split; $gs=join",",sort ($ga,$gb); next if($did{$gs}++);\
($sa)=split"_",$ga; ($sb)=split"_",$gb; map{ s/yfish/y/; s/aodon/ad/; s/leback/lb/; }($sa,$sb); \
if($orpar eq "orlog1") { $inp{$sa}{$sb}{$orid}++; $inp{$sb}{$sa}{$orid}++; $org{$orid}++; } \
elsif($orpar=~/inpar/){ map{ $inp{$sa}{$sa}{$orid}++ unless($didp{$_}++); }($ga,$gb); } \
$ngn++; END{ @sp=sort keys %inp; @og=sort keys %org; print join("\t","species",@sp)."\n";  \
for $sp (@sp){ print $sp; for $sb (@sp){ $c=0; for $o (@og) { $c += $inp{$sp}{$sb}{$o}; }\
print "\t$c"; } print "\n";}}' \
 > ${phyla}-orinpar-count.tab
  
=cut



sub omcl_1to1tab
{
  ## FIXME: sppindex for 1to1tab .. is it needed? caller opt? default?
  my $findex=$ENV{sppindex}||""; # "0,2,3,4,5,7,8,9,10"; # FIXME opt
  my $TMIN=$ENV{tmin}||0; #fixme opt 
  my @fi=();
  if($findex) { @fi=split/[,\s]+/,$findex; } # check-err
  if(@fi) { $TMIN= @fi-2 if($TMIN<1); $TMIN= @fi if($TMIN>@fi);  }

  my $infile="${phyla}-orthomcl-count.tab";
  my $outfile="${phyla}-orthomcl-1to1n.tab";
  open(I,$infile) or die "ERR: $infile";
  open(O,'>',$outfile) or die "ERR: write $outfile";
  my(@hd,@spp,$ncomm);
  while(<I>) {
    my($oid,$nt,$ng,@c)=split; 
    if(/^OID/) { @hd=split; @spp=@c; print O $_; next; }
    unless(@fi) {
      $#fi=scalar(@c); $TMIN= @fi-2 if($TMIN<1); $TMIN= @fi if($TMIN>@fi); 
      for my $i (0..$#c) { $fi[$i]=$i; }; 
    }
    my($fc,$fg)=(0,0); for my $i (@fi) { if($c[$i]) { $fc++; $fg+=$c[$i]; } } 
    next if($fc<$TMIN); $ncomm++; print O $_ if($fg == $fc);  
  } close(I); close(O);
}

=item omcl_1to1tab

  # 1:1 orlogs for 7..9 fish, ignoring spotgar (bad genes) and human
  cat ${phyla}-orthomcl-count.tab | perl -ne'BEGIN{ $TMIN=7; @fi=(0,2,3,4,5,7,8,9,10); }\
   ($oid,$nt,$ng,@c)=split; if(/^OID/) { @hd=split; @spp=@c; print; next; }\
   $fc=$fg=0; for $i (@fi) { if($c[$i]) { $fc++; $fg+=$c[$i]; } } \
   next if($fc<$TMIN); $ncomm++; print if($fg == $fc); END{ } ' \
  > ${phyla}-orthomcl-1to1n.tab

=cut

sub omcl_orpartab
{
  my($ORSPP)=@_; # 0 or 1 ?
  
  my $infile="${phyla}_omclgns2.tab";
  my $outfile="${phyla}-" . (($ORSPP)?"ors":"or1") . "inpar-count.tab"; # ORSPP=1
  # $outfile="${phyla}-or1inpar-count.tab" unless($ORSPP); # ORSPP=0
  open(I,$infile) or die "ERR: $infile";
  open(O,'>',$outfile) or die "ERR: write $outfile";
  
  my( %didp,%did,$ngn,%org,%inp);
  my($evo1,$pio1,$pavo1,$sao1,$sbo1)= (0) x 10;
  while(<I>) {
    my($ga,$gb,$orid,$ev,$pi,$orpar,$pav)=split; 
    my $gs=join",",sort ($ga,$gb); next if($did{$gs}++);
    my($sa)=split"_",$ga; my($sb)=split"_",$gb; 
    
if(1) { #Keep this way.
    if($ORSPP and $orpar =~ /orlog/ and not $inp{$sa}{$sb}{$orid}) { 
      $inp{$sa}{$sb}{$orid}++; $inp{$sb}{$sa}{$orid}++; $org{$orid}++; } 
    elsif($orpar eq "orlog1" and not $ORSPP) { 
      $inp{$sa}{$sb}{$orid}++; $inp{$sb}{$sa}{$orid}++; $org{$orid}++; } 
    elsif($orpar=~/inpar/){ 
      map{ $inp{$sa}{$sa}{$orid}++ unless($didp{$_}++); }($ga,$gb); 
    } 
    
} else {
    #FIX??: tie-namesort effect for orlog1 counts, tested method 2, not much change.. skip for dups
  if($orpar=~/inpar/){  
    map{ $inp{$sa}{$sa}{$orid}++ unless($didp{$_}++); }($ga,$gb); 
  } elsif($ORSPP) {
    if($orpar =~ /orlog/ and not $inp{$sa}{$sb}{$orid}) { 
      $inp{$sa}{$sb}{$orid}++; $inp{$sb}{$sa}{$orid}++; $org{$orid}++; 
      } 
  } else { # not ORSPP or inpar
    if($orpar eq "orlog1") { 
      $inp{$sa}{$sb}{$orid}++; $inp{$sb}{$sa}{$orid}++; $org{$orid}++; 
      $evo1= $ev; $pio1= $pi; $pavo1=$pav; $sao1=$sa; $sbo1=$sb;
    } elsif($orpar =~ /orlog/ 
      and $ev == $evo1 and $pi == $pio1 and $pav == $pavo1
      and not($sa eq $sao1 and $sb eq $sbo1)) {  ## $ev,$pi
      $inp{$sa}{$sb}{$orid}++; $inp{$sb}{$sa}{$orid}++; $org{$orid}++; 
    }
  }
}    
    
    $ngn++; 
  } close(I);

  #fixme.use printf names; # map{ s/yfish/y/; s/aodon/ad/; s/leback/lb/; }($sa,$sb); 
  my @sp=sort keys %inp; my @og=sort keys %org; 
  printf O "%-12s","species"; print O join("\t",@sp)."\n";  
  for my $sp (@sp){ 
    printf O "%-12s",$sp; 
    for my $sb (@sp){ my $c=0; for my $o (@og) { $c += $inp{$sp}{$sb}{$o}; } print O "\t$c"; } 
    print O "\n"; 
  } close(O);
  
}

sub omcl_renumber
{
  my($GTAG,$gntab,$omclout)= @_;
  (my $omclnew= $omclout) =~ s/\.out//; $omclnew.=".renum";
  open(GN,$gntab) or die "open $gntab";
  open(OM,$omclout) or die "open $omclout";
  open(ONEW,'>',$omclnew) or die "write $omclnew";
  my(%gog,$logn);
  while(<GN>) {
    my($og,$gid)=split; $og=~s/$GTAG//; 
    my($ogn)= $og=~m/(\d+)$/; $ogn="0z" if($ogn eq "0"); 
    $gog{$gid}=$ogn; $logn=$ogn; 
  } close(GN);
  
  my(%onew);
  while(<OM>) {
    if(/^ORTHOMCL/) { 
      my($om,$gn)=split /:\s+/,$_,2;  my @gn= split" ",$gn; 
      my($og,$ng,$nt)= $om =~ m/ORTHOMCL(\d+).(\d+) genes,\s*(\d+) taxa/; 
      $og="0z" if($og eq "0"); 
      my %oog=(); for my $g (@gn) { my $oog=$gog{$g}; $oog{$oog}++; } 
      my($oom)=sort{$b<=>$a} keys %oog; 
      $oom= ++$logn unless($oom); while($onew{$oom}) { $oom= ++$logn; } 
      $onew{$oom}=$og; # $nold{$og}=$oom; 
      print ONEW "ORTHOMCL$oom($ng genes,$nt taxa):\t $gn"; 
    } else { print ONEW $_; } # what?  
  } 
  close(OM); close(ONEW);
}

=item omcl_renumber
  Renumber omcl groups to match old version (nearly same.. added 10+ kfish2 missed genes..)

  cat ../fish11gor1/fish11g_omclgn.tab fish11gDec25/all_orthomcl.out | env gtag=FISH11G perl -ne \
  'BEGIN{ $GTAG=$ENV{gtag}; } if(/^$GTAG/){ ($og,$gid)=split; $og=~s/$GTAG//; \
  ($ogn)= $og=~m/(\d+)$/; $ogn="0z" if($ogn eq "0"); $gog{$gid}=$ogn; $logn0=$logn=$ogn; } \
  elsif(/^ORTHOMCL/) { my($om,$gn)=split /:\s+/,$_,2;  my @gn= split" ",$gn; \
  my($og,$ng,$nt)= $om =~ m/ORTHOMCL(\d+).(\d+) genes,\s*(\d+) taxa/; $og="0z" if($og eq "0"); \
  %oog=(); for my $g (@gn) { $oog=$gog{$g}; $oog{$oog}++; } ($oom)=sort{$b<=>$a} keys %oog; \
  $oom= ++$logn unless($oom); while($onew{$oom}) { $oom= ++$logn; } \
  $nold{$og}=$oom; $onew{$oom}=$og; print "ORTHOMCL$oom($ng genes,$nt taxa):\t $gn"; }' \
  > fish11gDec25/all_orthomcl.gor1.renum

=cut

#............. 
# mcl gene group tree info

sub omcl_mcl2cluster
{

  chdir($orun);
  sysrun("$mcl/bin/mclcm", 'tmp/all_ortho.mtx', '-a', '-I $mcl2INFLATE --shadow=vl -te 2'); # >& log.mcm1 # fork
  sysrun("$mcl/bin/mcxdump", '-imx-tree mcl.cone', '--newick', '-o', "${phyla}.newick", '-tab', 'tmp/all_ortho.idx');
  chdir("../");
  
=item mcl newick tree

  # cat ${phyla}_omclgn.tab ${phyla}.newick | perl -ne .. > ${phyla}.newick4
  set newk=${phyla}.newicki2

  ## improve these perls, mcltreesplit.pl
  cat ${phyla}_omclgn.tab  ${phyla}_omclgn.consensus_def.txt  ${newk} |
  env gtag=$GTAG perl -ne 'BEGIN{$gtag=$ENV{gtag};} 
  $p=1; if(/^$gtag\d\D+(\d+)/){ $g=$gtag.$1; s/\(LOC\d+\)//;  s/src=\S+//; 
  $gd{$g}= (m/\s(\S.+)$/)? $1.";" : "";  $p=0; }  
  elsif(/^($gtag\d+)\s+(\S+)$/){$ag{$2}= $1; $p=0;} 
  elsif(m/\(([^\)]+)\)/) { $gs=$1; @g=split",",$gs; %ga=(); $de="";  
  map{ $ag=$ag{$_}; $ga{$ag}++ if($ag); } @g; 
  map{ $de .= $gd{$_} } sort keys %ga; $ga= join ",", sort keys %ga;
  s/\([^\)]+\)/\($ga\)/; s/$/ # $de/ if $de; } print if $p; '
   > ${newk}b
  
  cat ${newk}b | env gtag=$GTAG perl -ne
  '$s=$_; s/\s*\#.*$//; s/[\,\s]+$//; $s=~s/$GTAG/$IDPRE/g; print $s unless($_ eq $ll); $ll=$_;' 
   > ${newk}c
   
  cat ${newk}c | env tag=$GTAG  perl ../omclw/mcltreesplit.pl > ${phyla}.clusters.tab

=cut

=item add clusters to ugp.xml

  cat  ${phyla}.clusters.tab  ${phyla}_genes.ugp.xml |
  perl -ne'if(/^ARC1_/){ chomp; ($c,$gc)=split"\t"; while( $gc =~ m/($GTAG\w+)/g) { $gc{$1}= $gc; } } 
  elsif(/^ARDE_/){next;} else { if(m/<GeneSummary/) { $gc= (m/id=.\w+:(\w+)/) ? $gc{$1} : ""; } 
  elsif( $gc and m,</GeneSummary>,){ print "<related_gene_groups>\n$gc\n</related_gene_groups>\n"; }  
  print; }'
  >  ${phyla}_genes.ugp.xml2

=cut

}

=item add

genegroupvenn/
  waspinsect_genefam_venn.pdf
  venn diagram of genes for 5 species, shared and unique (wasp,bee,ant,beetle,aphid)
  from http://bioinformatics.psb.ugent.be/webtools/Venn/
 
overgroups/   
  table.overgroups.nasonia.txt = list of over/under abundant groups in nasonia with statistics
  table.overgroups.* = same for other species

=cut

__END__

