#!/usr/bin/env perl
# # eval_orthodb_groups.pl

=item about

Evigene script to evaluate gene sets for match to orthology groups.
This is replacement to work with OrthoDB public tables, from  eval_orthogroup_genesets for orthomcl

  Inputs: xxx (from orthodb)

  Output: toptab same columns as from eval_orthogroup_genesets.pl / orthomcl, but
          only protein sizes and gene groups are valid,
          blast/align scores are fake (not avail from orthodb), and "top" gene/species-group is random, or aasize-related

=item  revise for input of orthodb tabtext of groups/genes

  gzcat OrthoDB6_Arthropoda_tabtext.gz | cut -f2,4,5,6 | grep -v '^ODB6_OG_ID' | perl -ne \
  '($od,$pid,$spf)=split"\t"; ($gs,$spl)=split " ",$spf; $sp=uc(substr($gs,0,1).substr($spl,0,4));\
  push @{$gn{$od}{$sp}},"$sp:$pid"; $od{$od}{$sp}++; $sp{$sp}++; $n++; END{ @od=sort keys %od;\
  @sp=sort keys %sp; print join("\t",qw(Query Source Bits Ident Align Qlen OID))."\n"; \
  foreach $d (@od) { $ng=$nt=0; foreach $sq (@sp) { @sg= @{$gn{$d}{$sq}}; next unless(@sg); \
  foreach $tq (@sp) { next if($tq eq $sq); @tg=@{$gn{$d}{$tq}}; foreach $sg (@sg) { \
  foreach $tg (@tg) { print join("\t",$sg,$tg,100,100,100,100,$d)."\n"; } } } } } }' \
  >  arp45-odball.fake.tall4  
  
  gzcat OrthoDB6_Arthropoda_tabtext.gz | cut -f2,4,5,6 | head -4
  ODB6_OG_ID      Gene_ID Organism        UniProt_Species
  EOG600001       ISCW001939      Ixodes scapularis       IXOSC
  EOG600001       ISCW003878      Ixodes scapularis       IXOSC
  EOG600001       ISCW004679      Ixodes scapularis       IXOSC
  
  ODB6_Level      ODB6_OG_ID      Protein_ID      Gene_ID Organism        UniProt_Species UniProt_ACC     UniProt_Description     InterPro_domains
  Arthropoda      EOG600001       ISCW001939-PA   ISCW001939      Ixodes scapularis       IXOSC   B7P7K5  Putative uncharacterized protein        NULL

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

my $idprefix="OID";
my $skipspecies="";
my $keepspecies="";
# my $FixSpeciesIDt2p="wasp"; # see fixids()

my $digits=1;
my $MINTAXA=3;
my $USE_MEDBITS=1; # default?
my $DOALL= 0;
my ($genescore_tall, $genegroups, $groupcounts, $geneaa, $topout, $intab, $sppidmap);

my $optok= GetOptions(
  # "tallscores=s", \$genescore_tall, # GONE
  "intab=s", \$intab,
  "output=s", \$topout,
  "groupgenes=s", \$genegroups, # ** in odb
  # "groupcounts=s", \$groupcounts, # GONE
  "geneaa=s", \$geneaa,#** in odb
  "skipspecies=s", \$skipspecies,
  "keepspecies=s", \$keepspecies,
  "sppidmap=s", \$sppidmap,
  "digits=s", \$digits,
  "idprefix=s", \$idprefix,
#  "FixSpeciesIDt2=s", \$FixSpeciesIDt2p,
  "MINTAXA=i", \$MINTAXA,
  "bitmedian!", \$USE_MEDBITS, # group median bitscore rather than highest pair score in group
  "tiny1|big1!", \$AA_OUTLIER1, # 
  "allgene!", \$DOALL, # all gene or top gene/spp/group
  "debug!", \$DEBUG, 
  );


die "usage: eval_orthodb_groups.pl  -groupgenes OrthoDB6_Arthropoda_tabtext.gz  -geneaa odbarthrogene.aa.count 
  opts: -mintaxa $MINTAXA -idprefix $idprefix ... \n"
  unless($optok and ($genegroups or $intab)); ## and -d $orun and -f "$orun/all_orthomcl.out"); #  and $orun

# evigene_config($config, \@configadd); # always even if $config null

$skipspecies =~ s/[,\s]+/|/g if($skipspecies);
$keepspecies =~ s/[,\s]+/|/g if($keepspecies);

my(%sppidmap,@sppidkey);
if($sppidmap) {
  map{ my($k,$v)=split /[=:]/, $_; $sppidmap{$k}=$v if($v); } split /[, ]/, $sppidmap;
  @sppidkey= sort keys %sppidmap;
}

if($intab) {
  topstats($intab); #? drop use eval_orthogroup_genesets
  
} else {
#  #20121009: ? revise this to input orthodb groups, geneaa ?
#  my( $ogenegrouph, $ogenesizeh)= orthogroupgenes($groupcounts,$geneaa,$genegroups); # odb CHANGE
#  $topout= topmcl($genescore_tall, $ogenegrouph, $ogenesizeh,$topout); ## odb CHANGE
# >> should handle all species here, and geneaa may be file list? directory?

  # this is "best" top gene table, add stats for  all genes/species/group
  my($topout1, $ogenegrouph,  $ogenesizeh, $ogroup, $ogcount, $species)
    = topodb($genegroups, $geneaa, $topout); ## odb CHANGE
  # return($outname, \%god, \%aasize, \%ogroup, \%ogcount, \%species);

  topstats($topout1);   # odb CHANGE
}

#..........

sub fixids
{
  local $_= shift; # id
  #? s/_/:/ unless(/:/); 
  # return separate spp: prefix if found?
  
  if($sppidmap) {
    unless(/^\w+:/) { 
      my $id=$_; my $idt="";
      foreach my $sp (@sppidkey) { if($id =~ m/^$sp/) { my $t=$sppidmap{$sp}; $idt="$t:$id" if($t); last; } }
      $_= $idt if($idt);
    }
   } 

  # s/t(\d+)$/p$1/ if($FixSpeciesIDt2p and m/$FixSpeciesIDt2p/); # damn this fixup crap
  if(wantarray) { return split(":",$_,2); }
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

sub median {
  my($arr)= @_;
  my ($n,$sm,$ss,$md,$mn,$sd,$min,$max)=(0) x 9;
  my @sl= sort{$b<=>$a} @$arr; # big first
  $n=@sl; ($max,$min)=  @sl[0,-1];
  $md= @sl[int($n/2)] if($n>0); 
  return wantarray ? ($md,$n,$min,$max) : $md;
}


# odb: merge sub orthogroupgenes() and topmcl() for topodb()
sub topodb
{
  my($genegroups,$geneaa,$outname)= @_;
  
  my(%ogroup,%species,%ogcount,%aasize,%pid2gid,%god,%odg,$nwarn);
  warn "#i orthogroups from $genegroups\n" if($DEBUG);

#   open(CIN,"$groupcounts") or die "$groupcounts";
#   while(<CIN>) { my($od,$nt,$ng,@c)=split; $ogroup{$od}++ if($nt>=$MINTAXA); } close(CIN);
#   warn "#i group n=",scalar(keys %ogroup)," with ntaxa>=$MINTAXA\n" if($DEBUG);

## groupcounts from genegroups now? or use 2ndary table from same source?
# gzcat OrthoDB6_Arthropoda_tabtext.gz | cut -f2,4,5,6 | grep -v '^ODB6_OG_ID' | perl -ne\
# '($od,$pid,$spf)=split"\t"; ($gs,$spl)=split " ",$spf; $sp=uc(substr($gs,0,1).substr($spl,0,4)); \
# $od{$od}{$sp}++; $sp{$sp}++; $n++; END{ @od=sort keys %od; @sp=sort keys %sp; \
# print join("\t",qw(OID Nt Ng),@sp)."\n"; foreach $d (@od) { $ng=$nt=0; @c= map{ $c=$od{$d}{$_}||0; \
# $ng+=$c; $nt++ if($c>0); $c; }@sp; print join("\t",$d,$nt,$ng,@c)."\n"; } }' \
# > arp45-orthodb6-count.tab

  my $inh=undef;
  if($genegroups =~ /stdin|^\-/) { $inh= *STDIN; }
  elsif($genegroups =~ /\.gz/) { open(OGENES,"gunzip -c $genegroups|") or die "$genegroups"; $inh= *OGENES; }
  else { open(OGENES,"$genegroups") or die "$genegroups"; $inh= *OGENES; }
  while(<$inh>) {
    # ODB6_Level ODB6_OG_ID Protein_ID Gene_ID Organism  UniProt_Species UniProt_ACC UniProt_Description   InterPro_domains
    my($od,$pid,$gn,$spf)=(split"\t")[1,2,3,4]; 
    next if($gn eq "Gene_ID");
    my($gs,$spl)=split " ",$spf; 
    my $sp=uc(substr($gs,0,1).substr($spl,0,4));  
    next if($keepspecies and not($sp =~ m/$keepspecies/)); # both keep,skip?
    next if($skipspecies and $sp =~ m/$skipspecies/);
    $gn= fixids("$sp:$gn");
    $pid2gid{"$sp:$pid"}=$gn; # for aasize id mixups
    $ogroup{$od}{$sp}++;
    $species{$sp}++; # count genes/spp or groups/spp?
    $god{$gn}=$od; #?? if($ogroup{$od}); # save od or odnum ?
    push @{$odg{$od}}, $gn; # for delete bad ogroup
  } 
  close($inh);  

    ## FIXME: if no ogroup of group counts, need to post-process ogroup/god hash for MINTAXA groups
  my $nogdrop=0;
  my @ogroup=sort keys %ogroup; my @sp=sort keys %species;
  foreach my $od (@ogroup) { 
    my($ng,$nt)=(0,0); 
    my @c= map{ my $c=$ogroup{$od}{$_}||0; $ng+=$c; $nt++ if($c>0); $c; }@sp;
    if($nt<$MINTAXA) { 
      my @dgn= @{$odg{$od}}; map{ delete $god{$_}; } @dgn;
      delete $odg{$od}; delete $ogroup{$od};  $nogdrop++;
      }
    $ogcount{$od}="$nt\t$ng".join("\t",@c);  #? save nt, ng
  }
        
  warn "#i groups n=",scalar(keys %ogroup)," with ntaxa>=$MINTAXA; ndrop=$nogdrop\n" if($DEBUG);
  warn "#i groupgenes n=",scalar(keys %god),"\n" if($DEBUG);
  warn "#i species n=",scalar(keys %species),"\n" if($DEBUG);

# geneaa here may be file list? at least multispp data; require spp:geneid prefix for aa? aa.count?  
  my %speciesaa;
  my $iscount=0;
  if($geneaa =~ /count/) { # bad but ok..
    open(AASIZE,$geneaa) or die "aa count in $geneaa ..."; $iscount=1;
  } else {
    open(AASIZE,"faCount $geneaa | cut -f1,2 |") or die "faCount  $geneaa ...";
  }
  while(<AASIZE>) {
    my($id,$al)=split; 
    next unless($id and $al>0); # if($iscount and (not defined($al) or $al=~/\D/)) {}
    # FIXME: here? BIMPA BTERR bombus imp/terr have NCBI prot ids in aasize, "gene id" is ?ncbi locus num
    if($pid2gid{$id}) { $id= $pid2gid{$id}; } # pid has sp: prefix
    $id= fixids($id); # MUST have species prefix here
    my($sp)=split ":",$id; 
    unless($god{$id}) { $nwarn++; next; }  # warn "# aa nogroup $id,$al\n" if($nwarn<10);
    $aasize{$id}=$al; # need to check valid species, some missing aa ..
    $speciesaa{$sp}++;
  } close(AASIZE);
  
  my $sppcount= join", ", map{ "$_:$speciesaa{$_}" } sort keys %speciesaa;
  warn "#i geneaa n=",scalar(keys %aasize)," nogeneaa=$nwarn, sppcount $sppcount\n" if($DEBUG);
  
  # check species not in speciesaa, drop missing from odg, ogroup
  my @nospp= grep{ !$speciesaa{$_} } sort keys %species;
  if(@nospp) {
    warn "#i drop no-species-aasize n=",scalar(@nospp)," : @nospp\n" if($DEBUG);
    @ogroup=sort keys %ogroup; @sp=sort keys %speciesaa;
    foreach my $od (@ogroup) { 
      my($ng,$nt)=(0,0); 
      my @c= map{ my $c=$ogroup{$od}{$_}||0; $ng+=$c; $nt++ if($c>0); $c; }@sp;
      if($nt<$MINTAXA) { 
        my @dgn= @{$odg{$od}}; map{ delete $god{$_}; } @dgn;
        delete $odg{$od}; delete $ogroup{$od};  $nogdrop++;
      } else {
        map{ delete $ogroup{$od}{$_} } @nospp;
      }   
      $ogcount{$od}="$nt\t$ng".join("\t",@c);  #? save nt, ng
    }
  }
  
  # write topout table here  
  unless($outname) { 
    ($outname=$genegroups) =~ s/\.\S+$//; 
    $outname .=  ($DOALL) ? ".allordb" : ".topordb"; 
  }
  putToptab($outname, \%odg, \%aasize, \%ogroup);  # \%god , \%ogcount, \%species

  return($outname, \%god, \%aasize, \%ogroup, \%ogcount, \%species); # \%odg
}



sub puta {
  my($oid, $tg, $sg, $tgsize, $sgsize, $mdsize)= @_;  
  my @vfake= (199, 99, 99); # Bits Iden Algn
  #?? add 1, -1 levels
  my $dl= $tgsize - $mdsize;
  my $xl=($tgsize > $AA_TOOBIG*$mdsize) ? 2 : ($tgsize < $AA_TOOSMALL*$mdsize) ? -2 : 0;  
  if($AA_OUTLIER1 and $xl == 0) { #? move this to topstats as output filter?
    $xl= ($tgsize > $AA_TOOBIG1*$mdsize) ? 1 : ($tgsize < $AA_TOOSMALL1*$mdsize) ? -1: 0;  
  }
  print OUT join("\t",$tg,$sg,$sgsize,$oid, @vfake, $tgsize ,$mdsize,$dl,$xl)."\n"; 
} 


  # this is "best" top gene table, add stats for  all genes/species/group
sub putToptab { 
  my($outname, $odg, $aasize, $ogroup)= @_;  # \%god , \%ogcount, \%species
  # unless($outname) { ($outname=$genegroups) =~ s/\.\S+$//; $outname .=".topordb"; }
  open(OUT,">$outname") or die "$outname";
  print OUT join("\t",qw(TargGeneid SrcGeneid Slen Ogroup Bit Idn Aln Tlen Olen dTO xTO))."\n";  
  
  my @ogroup= sort keys %$ogroup; # recount
  my (%did);
  foreach my $od (@ogroup) {  
    my @dgn= @{$odg->{$od}}; # all species genes/group
    my @dgnsizes= grep /\d/, map{ $aasize->{$_} } @dgn; # REMOVE zeros from missing spp
    
    #? fixme: do this per species, median/mean excluding put-species?
    # my($naa,$mdaa,$mnaa,$sdaa,$minaa,$maxaa)= medmeansd( \@dgnsizes );
    my($mdaa,$naa,$minaa,$maxaa)= median( \@dgnsizes );
    
    # need some selection here for "best" gene/species/group
    # .. pick species median genesize ; or closest to group median?
    my(%bestaa,%bestaad,%allgn);
    foreach my $g (@dgn) {
      my($sp)= split":",$g;
      my $gsize= $aasize->{$g} or next; # maybe missing?
      my $daa;  

# OTHERSPP_MEDIAN little effect here, dont bother
use constant OTHERSPP_MEDIAN => 0; 
if(OTHERSPP_MEDIAN) {      
      my @otheraa= grep /\d/, map{ $aasize->{$_} } grep !/$sp:/, @dgn;
      my($mdaaOther,$naa,$minaa,$maxaa)=  median( \@otheraa );
      $daa= abs($gsize - $mdaaOther);  
} else {
      $daa= abs($gsize - $mdaa);  
}
      push @{$allgn{$sp}}, $g;
      if( (!defined $bestaad{$sp}) or ($bestaad{$sp} > $daa)) { $bestaa{$sp}= $g; $bestaad{$sp}= $daa; }
      # add stats for spp size range.. minaa, maxaa, naa/meanaa ? noutlier?
    }
      
    foreach my $sp (sort keys %bestaa) {
      my ($sp2fake)= grep { $_ ne $sp } sort keys %bestaa; # sort?
      my $sg= $bestaa{$sp2fake};
      my $sgsize= $aasize->{$sg};
      if($DOALL) {
        my @tg=  @{$allgn{$sp}};
        foreach my $tg (@tg) {
          my $tgsize= $aasize->{$tg};
          puta($od,$tg,$sg,$tgsize,$sgsize,$mdaa);
          $did{$tg}++; 
        }
      } else {
      # push other stats as new columns?  size min,max,n2tiny,n2big
        my $tg= $bestaa{$sp};
        my $tgsize= $aasize->{$tg};
        puta($od,$tg,$sg,$tgsize,$sgsize,$mdaa);
        $did{$tg}++; 
      }
    }
    $did{$od}++; 
  } 
  close(OUT);

}

  

=item toptab : keep this format, but Bit/Idn/Aln will be fake

  TargGeneid      SrcGeneid       Slen    Ogroup  Bit     Idn     Aln     Tlen    Olen    dTO     xTO
  killifish:Funhe5EG030598t1      zebrafish:E7FDA5_DANRE  743     FISH0   221     225     518     705     289     416     2
  killifish:Funhe5EG015351t1      medaka:ENSORLP00000008122       267     FISH10  162     185     281     434     267     167     0

=cut

# odb NO change here, but Bits score fake
sub topstats
{
  my($toptab)= @_;
  # let toptab == stdin
  
  my( %arn, %arng, %arx, %ard, %arb, %ts, %tsx);
  
  my $inh={};
  if($toptab =~ /^-|^stdin/) { $inh=*STDIN; }
  else { open(IN,$toptab) or die "$toptab"; $inh= *IN; }
  while(<$inh>) {
    my($bt,$xd,$d,$ts);
    next if(/^TargGeneid/); 
    my($tg,$sg,$sl,$ar,@v)=split; 
    # ts tag from geneid: need help:  acyp2ref, acyp2eg << not same...
    unless( ($ts)= $tg=~m/^(\w+):/) {
    unless( ($ts)= $tg=~m/^([^\d\W]+\d[^\d\W]+)/ ) { ($ts)= $tg=~m/^([^\d\W]+)/; }
    }
 
    # FIXME: bt is fake for orthodb table
    
    $bt=$v[0]; $xd=$v[-1]; $d=$v[-2]; 
    # >> fix here for many gene/spp/group : +=$d +=$bt ?? xd
    $arb{$ar}{$ts} += $bt; $ard{$ar}{$ts} += $d; 
    #move above: $xd=0 if( !$AA_OUTLIER1 and ($xd==1 or $xd==-1));
    #old? $arx{$ar}{$ts}=$xd if($xd); #? change to 3hash:  $arx{$ar}{$ts}{$xd}++
    $arx{$ar}{$ts}{$xd}++ if($xd);
    $arn{$ar}++; $arng{$ar}{$ts}++; $ts{$ts}++;   $tsx{$ts}{$xd}++; 
  } close($inh);

  my @ts=sort keys %ts; my $NG=@ts; my @ar=sort keys %arb; 
  my @ap=grep{ scalar(keys %{$arb{$_}})==$NG } @ar; 
  print "Summary stats for $toptab\n";  
  
   #printf "%-9s: ","Geneset"; print join("\t",qw(Ng Ngr Bits dSize rBits oddSize)),"\n";
   #print "n=$nt; nreal=$nr; abits=$ab; dlen=$ad; arealbits=$abr; xtralen:@xl\n";
  my %xlkey= ( 2=>"x2Big", -2 =>"x2Tiny", 1=>"x1big", -1 =>"x1tiny" );
  
  for my $k (0,1) {
    my @aset= ($k==0) ? @ap : @ar;
    print " ...  ",(($k==0)? "Common" : "All"), " groups ... \n";
    printf "%-9s: ","Geneset"; 
    # change nGene/rGene > nGroup/rGroup  add nGene for DOALL
    if($DOALL) {
    print join("\t",qw(nGroup Bits dSize rGroup nGene rBits outlierSize)),"\n";    
    } else {
    print join("\t",qw(nGene Bits dSize rGene rBits outlierSize)),"\n";
    }
    
    my (%sb,%sn,%sx,%sd,%snt,%sng);
    
    foreach my $ar (@aset) { foreach my $t (@ts) {  my($bt,$x,$d,$ng);
      $bt=$arb{$ar}{$t}; $d=$ard{$ar}{$t}; 
      #old# $x=$arx{$ar}{$t}; $sx{$t}{$x}++ if($x); # ok for 2 outlier levels here?  -2,+2 and -1,+1
      my @x= keys %{$arx{$ar}{$t}}; map{ $sx{$t}{$_}++ } @x; # new
      # FIXME: DOALL, need:  sng{$t} += $ng;
      $ng= $arng{$ar}{$t}; $sng{$t} += $ng;
      $sb{$t}+=$bt; $sn{$t}++; $sd{$t}+=$d; $snt{$t}++ if($bt); 
      } 
    } 
    
    my $decp=($digits>0)?10:1;

    foreach my $t (@ts) { my($nt,$ab,$nr,$abr,$ad,@xl,$ng);
      $nt=$sn{$t}||1; $nr=$snt{$t}||1; $ng= $sng{$t}||1;
      $ab=int(0.5+$decp*$sb{$t}/$nt)/$decp;  
      $abr=int(0.5+$decp*$sb{$t}/$ng)/$decp;  
      $ad=int(0.5+$decp*$sd{$t}/$ng)/$decp; # ng == nr unless DOALL : should be
      @xl= map{ my $xc=$sx{$t}{$_}; my $xp=int(1000*$xc/$nr)/10; $xlkey{$_}."=$sx{$t}{$_} ($xp%)" } sort{$a<=>$b} keys %{$sx{$t}}; 
      printf "%-9s: ",$t; 
      # print "n=$nt; nreal=$nr; abits=$ab; dlen=$ad; arealbits=$abr; xtralen:@xl\n";  
      if($DOALL) {
      print join("\t",$nt,$ab,$ad,$nr,$ng,$abr,@xl),"\n"; # add ng
      } else {
      print join("\t",$nt,$ab,$ad,$nr,$abr,@xl),"\n"; 
      }
    }
  }
  
}



__END__

=item test1

$evigene/scripts/eval_orthodb_groups.pl -group OrthoDB6_Hymenoptera_tabtext.gz -geneaa hym-orthodb6.a
a.count
#i orthogroups from OrthoDB6_Hymenoptera_tabtext.gz
#i groups n=13201 with ntaxa>=3; ndrop=3474
#i groupgenes n=138204
#i species n=14
# aa nogroup ACEPH:ACEP10001,177   # these are just genes not in ortho groups, 
# aa nogroup ACEPH:ACEP10002,149
# aa nogroup ACEPH:ACEP10005,238
# aa nogroup ACEPH:ACEP10014,99
# aa nogroup ACEPH:ACEP10017,101
# aa nogroup ACEPH:ACEP10022,312

#i geneaa n=102453
Summary stats for OrthoDB6_Hymenoptera_tabtext.topordb
 ...  Common groups ... 
 ** Fail here due to no-aa species BIMPA, BTERR, 
 cat: cannot open fasta/AFLOR.*.fa.count
 cat: cannot open fasta/BIMPA.*.fa.count  << add Ncbi gene.aa
 cat: cannot open fasta/BTERR.*.fa.count  << Ncbi gene.aa
 cat: cannot open fasta/MROTU.*.fa.count

Geneset  : nGene        Bits    dSize   rGene   rBits   outlierSize
ACEPH    : 1    99      209     1       99      x2Big=1 (100%)
AECHI    : 1    99      0       1       99      x2Big=1 (100%)
AFLOR    : 1    99      0       1       99      x2Tiny=1 (100%)
AMELL    : 1    99      462     1       99      x2Big=1 (100%)
BIMPA    : 1    99      -148.9  1       99      x2Tiny=1 (100%)
BTERR    : 1    99      -148.9  1       99      x2Tiny=1 (100%)
CFLOR    : 1    99      0       1       99      x2Big=1 (100%)
HSALT    : 1    99      0       1       99      x2Big=1 (100%)
LHUMI    : 1    99      0       1       99      x2Tiny=1 (100%)
MROTU    : 1    99      0       1       99      x2Tiny=1 (100%)
NVITR    : 1    99      4       1       99      x2Big=1 (100%)
PBARB    : 1    99      215     1       99      x2Big=1 (100%)
SINVI    : 1    99      89      1       99      x2Big=1 (100%)

 ...  All groups ... 
Geneset  : nGene        Bits    dSize   rGene   rBits   outlierSize
ACEPH    : 12965        162.8   19.7    10609   199     x2Tiny=601 (5.6%)       x2Big=346 (3.2%)
AECHI    : 12965        151.3   76.3    9857    199     x2Tiny=189 (1.9%)       x2Big=687 (6.9%)
AMELL    : 12965        142.1   83.6    9257    199     x2Tiny=141 (1.5%)       x2Big=669 (7.2%)
CFLOR    : 12965        39.3    41.9    2563    199     x2Tiny=183 (7.1%)       x2Big=201 (7.8%)
HSALT    : 12965        148.9   54.5    9703    199     x2Tiny=424 (4.3%)       x2Big=544 (5.6%)
LHUMI    : 12965        163.6   43.9    10661   199     x2Tiny=348 (3.2%)       x2Big=475 (4.4%)
NVITR    : 12965        137.6   84.9    8964    199     x2Tiny=219 (2.4%)       x2Big=702 (7.8%)
PBARB    : 12965        163.5   25.7    10652   199     x2Tiny=468 (4.3%)       x2Big=349 (3.2%)
SINVI    : 12965        150.4   1.1     9800    199     x2Tiny=942 (9.6%)       x2Big=280 (2.8%)

nodata# AFLOR    : 12965        124.3   -510.4  8100    199     x2Tiny=7958 (98.2%)
nodata# BIMPA    : 12965        124.7   -505.5  8126    199     x2Tiny=8059 (99.1%)
nodata# BTERR    : 12965        123.5   -489.4  8049    199     x2Tiny=7857 (97.6%)
nodata# MROTU    : 12965        130.3   -499.1  8489    199     x2Tiny=8447 (99.5%)

=item run2
 
 -- added BIMPA,BTERR aa.size, fixed missing spp aa size
 
$evigene/scripts/eval_orthodb_groups.pl -group OrthoDB6_Hymenoptera_tabtext.gz -geneaa hym-orthodb6.aa.count
#i orthogroups from OrthoDB6_Hymenoptera_tabtext.gz
#i groups n=13201 with ntaxa>=3; ndrop=3474
#i groupgenes n=138204
#i species n=14
# aa nogroup ACEPH:ACEP10001,177
# aa nogroup ACEPH:ACEP10002,149
...

#i geneaa n=102453 nogeneaa=78901, sppcount 
# ACEPH:11808, AECHI:10935, AMELL:10621, CFLOR:11151, HSALT:11622, LHUMI:11868, NVITR:11475, PBARB:11801, SINVI:11172
  ^^ What? missed BTERR, BIMPA -- missed spp: prefix above
#i drop no-species-aasize n=5 : AFLOR BIMPA BTERR MROTU O
Summary stats for OrthoDB6_Hymenoptera_tabtext.topordb
 ...  Common groups ... 
Geneset  : nGene  Bits    dSize   rGene   rBits   outlierSize
ACEPH    : 6134   199     22.4    6134    199     x2Tiny=285 (4.6%)       x2Big=142 (2.3%)
AECHI    : 6134   199     82.3    6134    199     x2Tiny=77 (1.2%)        x2Big=315 (5.1%)
AMELL    : 6134   199     83.1    6134    199     x2Tiny=84 (1.3%)        x2Big=343 (5.5%)
CFLOR    : 6134   199     78.4    6134    199     x2Tiny=145 (2.3%)       x2Big=334 (5.4%)
HSALT    : 6134   199     59.4    6134    199     x2Tiny=199 (3.2%)       x2Big=291 (4.7%)
LHUMI    : 6134   199     51.7    6134    199     x2Tiny=140 (2.2%)       x2Big=226 (3.6%)
NVITR    : 6134   199     88      6134    199     x2Tiny=126 (2%)         x2Big=372 (6%)
PBARB    : 6134   199     28      6134    199     x2Tiny=230 (3.7%)       x2Big=170 (2.7%)
SINVI    : 6134   199     -8.1    6134    199     x2Tiny=694 (11.3%)      x2Big=137 (2.2%)
 ...  All groups ... 
Geneset  : nGene  Bits    dSize   rGene   rBits   outlierSize
ACEPH    : 12544  174.5   20.6    11002   199     x2Tiny=613 (5.5%)       x2Big=367 (3.3%)
AECHI    : 12544  161     76.9    10149   199     x2Tiny=194 (1.9%)       x2Big=724 (7.1%)
AMELL    : 12544  153.2   85.2    9654    199     x2Tiny=148 (1.5%)       x2Big=739 (7.6%)
CFLOR    : 12544  159     73.1    10025   199     x2Tiny=370 (3.6%)       x2Big=692 (6.9%)
HSALT    : 12544  159     55.6    10020   199     x2Tiny=432 (4.3%)       x2Big=586 (5.8%)
LHUMI    : 12544  175     44.6    11034   199     x2Tiny=358 (3.2%)       x2Big=505 (4.5%)
NVITR    : 12544  147.2   85.5    9278    199     x2Tiny=225 (2.4%)       x2Big=751 (8%)
PBARB    : 12544  175.1   26.4    11035   199     x2Tiny=481 (4.3%)       x2Big=379 (3.4%)
SINVI    : 12544  160.5   2.8     10119   199     x2Tiny=953 (9.4%)       x2Big=303 (2.9%)


=item run4 // 4 outlier classes: -2,-1,1,2 : make -1,1 option?

$evigene/scripts/eval_orthodb_groups.pl -group OrthoDB6_Hymenoptera_tabtext.gz -geneaa hym-orthodb6.aa.count
#i orthogroups from OrthoDB6_Hymenoptera_tabtext.gz
#i groups n=13201 with ntaxa>=3; ndrop=3473
#i groupgenes n=138204
#i species n=13
#i geneaa n=120294 nogeneaa=61060, sppcount ACEPH:11808, AECHI:10935, AMELL:10621, BIMPA:9043, BTERR:8798, CFLOR:11151, HSALT:11622, LHUMI:11868, NVITR:11475, PBARB:11801, SINVI:11172
#i drop no-species-aasize n=2 : AFLOR MROTU

Summary stats for OrthoDB6_Hymenoptera_tabtext.topordb
 ...  Common groups ... 
Geneset  : nGene   Bits    dSize   rGene   rBits   outlierSize
ACEPH    : 5512    199     -29.4   5512    199     x2Tiny=333 (6%)      x1tiny=331 (6%)    x1big=72 (1.3%)   x2Big=23 (0.4%)
AECHI    : 5512    199     25.7    5512    199     x2Tiny=89 (1.6%)     x1tiny=113 (2%)    x1big=154 (2.7%)  x2Big=67 (1.2%)
AMELL    : 5512    199     30.4    5512    199     x2Tiny=63 (1.1%)     x1tiny=87 (1.5%)   x1big=154 (2.7%)  x2Big=85 (1.5%)
BIMPA    : 5512    199     55.3    5512    199     x2Tiny=23 (0.4%)     x1tiny=28 (0.5%)   x1big=232 (4.2%)  x2Big=122 (2.2%)
BTERR    : 5512    199     54.2    5512    199     x2Tiny=45 (0.8%)     x1tiny=30 (0.5%)   x1big=221 (4%)    x2Big=131 (2.3%)
CFLOR    : 5512    199     22.5    5512    199     x2Tiny=131 (2.3%)    x1tiny=142 (2.5%)  x1big=170 (3%)    x2Big=78 (1.4%)
HSALT    : 5512    199     6.5     5512    199     x2Tiny=191 (3.4%)    x1tiny=202 (3.6%)  x1big=108 (1.9%)  x2Big=52 (0.9%)
LHUMI    : 5512    199     2.1     5512    199     x2Tiny=154 (2.7%)    x1tiny=150 (2.7%)  x1big=115 (2%)    x2Big=48 (0.8%)
NVITR    : 5512    199     35.7    5512    199     x2Tiny=110 (1.9%)    x1tiny=101 (1.8%)  x1big=195 (3.5%)  x2Big=102 (1.8%)
PBARB    : 5512    199     -21.8   5512    199     x2Tiny=258 (4.6%)    x1tiny=279 (5%)    x1big=90 (1.6%)   x2Big=29 (0.5%)
SINVI    : 5512    199     -61.3   5512    199     x2Tiny=750 (13.6%)   x1tiny=254 (4.6%)  x1big=84 (1.5%)   x2Big=41 (0.7%)
 ...  All groups ... 
Geneset  : nGene   Bits    dSize   rGene   rBits   outlierSize
ACEPH    : 12955   169.1   -18.2   11010   199     x2Tiny=732 (6.6%)    x1tiny=701 (6.3%)  x1big=329 (2.9%)  x2Big=165 (1.4%)
AECHI    : 12955   156.1   34.2    10159   199     x2Tiny=224 (2.2%)    x1tiny=254 (2.5%)  x1big=463 (4.5%)  x2Big=354 (3.4%)
AMELL    : 12955   153.7   39.4    10006   199     x2Tiny=139 (1.3%)    x1tiny=187 (1.8%)  x1big=439 (4.3%)  x2Big=292 (2.9%)
BIMPA    : 12955   134.5   56.9    8758    199     x2Tiny=47 (0.5%)     x1tiny=59 (0.6%)   x1big=442 (5%)    x2Big=268 (3%)
BTERR    : 12955   131.2   53.9    8540    199     x2Tiny=79 (0.9%)     x1tiny=62 (0.7%)   x1big=410 (4.8%)  x2Big=265 (3.1%)
CFLOR    : 12955   154.2   28.5    10037   199     x2Tiny=396 (3.9%)    x1tiny=358 (3.5%)  x1big=407 (4%)    x2Big=294 (2.9%)
HSALT    : 12955   154.4   13.9    10054   199     x2Tiny=469 (4.6%)    x1tiny=432 (4.2%)  x1big=327 (3.2%)  x2Big=219 (2.1%)
LHUMI    : 12955   169.7   5.3     11047   199     x2Tiny=427 (3.8%)    x1tiny=378 (3.4%)  x1big=380 (3.4%)  x2Big=188 (1.7%)
NVITR    : 12955   144.5   41      9408    199     x2Tiny=236 (2.5%)    x1tiny=201 (2.1%)  x1big=466 (4.9%)  x2Big=336 (3.5%)
PBARB    : 12955   169.8   -12.2   11052   199     x2Tiny=572 (5.1%)    x1tiny=572 (5.1%)  x1big=303 (2.7%)  x2Big=133 (1.2%)
SINVI    : 12955   155.7   -37.1   10134   199     x2Tiny=1110 (10.9%)  x1tiny=470 (4.6%)  x1big=273 (2.6%)  x2Big=133 (1.3%)

=item run5

$evigene/scripts/eval_orthodb_groups.pl -group OrthoDB6_Hymenoptera_tabtext.gz -geneaa hym-orthodb6.aa.count
#i orthogroups from OrthoDB6_Hymenoptera_tabtext.gz
#i groups n=13201 with ntaxa>=3; ndrop=3473
#i groupgenes n=138204
#i species n=13
#i geneaa n=120294 nogeneaa=61060, sppcount ACEPH:11808, AECHI:10935, AMELL:10621, BIMPA:9043, BTERR:8798, CFLOR:11151, HSALT:11622, LHUMI:11868, NVITR:11475, PBARB:11801, SINVI:11172
#i drop no-species-aasize n=2 : AFLOR MROTU
Summary stats for OrthoDB6_Hymenoptera_tabtext.topordb
 ...  Common groups ... 
Geneset  : nGene  Bits    dSize   rGene   rBits   outlierSize
ACEPH    : 5512   199     -29.4   5512    199     x2Tiny=333 (6%)      x2Big=23 (0.4%)
AECHI    : 5512   199     25.7    5512    199     x2Tiny=89 (1.6%)     x2Big=67 (1.2%)
AMELL    : 5512   199     30.4    5512    199     x2Tiny=63 (1.1%)     x2Big=85 (1.5%)
BIMPA    : 5512   199     55.3    5512    199     x2Tiny=23 (0.4%)     x2Big=122 (2.2%)
BTERR    : 5512   199     54.2    5512    199     x2Tiny=45 (0.8%)     x2Big=131 (2.3%)
CFLOR    : 5512   199     22.5    5512    199     x2Tiny=131 (2.3%)    x2Big=78 (1.4%)
HSALT    : 5512   199     6.5     5512    199     x2Tiny=191 (3.4%)    x2Big=52 (0.9%)
LHUMI    : 5512   199     2.1     5512    199     x2Tiny=154 (2.7%)    x2Big=48 (0.8%)
NVITR    : 5512   199     35.7    5512    199     x2Tiny=110 (1.9%)    x2Big=102 (1.8%)
PBARB    : 5512   199     -21.8   5512    199     x2Tiny=258 (4.6%)    x2Big=29 (0.5%)
SINVI    : 5512   199     -61.3   5512    199     x2Tiny=750 (13.6%)   x2Big=41 (0.7%)
 ...  All groups ... 
Geneset  : nGene  Bits    dSize   rGene   rBits   outlierSize
ACEPH    : 12955  169.1   -18.2   11010   199     x2Tiny=732 (6.6%)    x2Big=165 (1.4%)
AECHI    : 12955  156.1   34.2    10159   199     x2Tiny=224 (2.2%)    x2Big=354 (3.4%)
AMELL    : 12955  153.7   39.4    10006   199     x2Tiny=139 (1.3%)    x2Big=292 (2.9%)
BIMPA    : 12955  134.5   56.9    8758    199     x2Tiny=47 (0.5%)     x2Big=268 (3%)
BTERR    : 12955  131.2   53.9    8540    199     x2Tiny=79 (0.9%)     x2Big=265 (3.1%)
CFLOR    : 12955  154.2   28.5    10037   199     x2Tiny=396 (3.9%)    x2Big=294 (2.9%)
HSALT    : 12955  154.4   13.9    10054   199     x2Tiny=469 (4.6%)    x2Big=219 (2.1%)
LHUMI    : 12955  169.7   5.3     11047   199     x2Tiny=427 (3.8%)    x2Big=188 (1.7%)
NVITR    : 12955  144.5   41      9408    199     x2Tiny=236 (2.5%)    x2Big=336 (3.5%)
PBARB    : 12955  169.8   -12.2   11052   199     x2Tiny=572 (5.1%)    x2Big=133 (1.2%)
SINVI    : 12955  155.7   -37.1   10134   199     x2Tiny=1110 (10.9%)  x2Big=133 (1.3%)


=item run6 : subset of OrthoDB6_Arthropoda_tabtext

$evigene/scripts/eval_orthodb_groups.pl -group OrthoDB6_Arthropoda_tabtext.gz  -geneaa arp45-orthodb6.aa.count
#i orthogroups from OrthoDB6_Arthropoda_tabtext.gz
#i groups n=25402 with ntaxa>=3; ndrop=7991
#i groupgenes n=513391
#i species n=45
#i geneaa n=421013 nogeneaa=202129, sppcount AAEGY:12852, ACEPH:12406, ADARL:8842, AECHI:11281, AGAMB:11042, AMELL:10908, APISU:15578, ASTEP:9520, BMORI:11366, CFLOR:11790, CQUIN:13543, DANAN:13223, DEREC:13674, DGRIM:13056, DMELA:13285, DMOJA:12710, DPERS:13428, DPSEU:13482, DPULE:11299, DSECH:14265, DSIMU:13374, DVIRI:12869, DWILL:13051, DYAKU:14151, HMELP:10821, HSALT:12395, ISCAP:9235, LHUMI:12285, NVITR:13058, PBARB:12250, PHUMA:9126, RPROL:9307, SINVI:11831, TCAST:11281, TURTI:8429
#i drop no-species-aasize n=10 : AFLOR BIMPA BTERR DPLEX MDEST MMOLD MROTU MSEXT SMARI ZNEVA

* subset of arthropod species.
Summary stats for OrthoDB6_Arthropoda_tabtext.topordb
 ...  Common groups ... 
Geneset  : nGene  Bits    dSize   rGene   rBits   outlierSize
AMELL    : 1328   199     46.1    1328    199     x2Tiny=18 (1.3%)    x2Big=19 (1.4%)
APISU    : 1328   199     36.9    1328    199     x2Tiny=13 (0.9%)    x2Big=13 (0.9%)  << ok here, rel PHUMA
DMELA    : 1328   199     109.4   1328    199     x2Tiny=2 (0.1%)     x2Big=41 (3%)
DPULE    : 1328   199     -1.4    1328    199     x2Tiny=66 (4.9%)    x2Big=10 (0.7%)   << smallish, old JGI gene set fragments
HSALT    : 1328   199     33.1    1328    199     x2Tiny=32 (2.4%)    x2Big=15 (1.1%)
ISCAP    : 1328   199     -99.8   1328    199     x2Tiny=219 (16.4%)  x2Big=6 (0.4%)    << gene fragments
NVITR    : 1328   199     47.7    1328    199     x2Tiny=22 (1.6%)    x2Big=18 (1.3%)   << ok here, rel AMELL, HSALT
PHUMA    : 1328   199     31.6    1328    199     x2Tiny=22 (1.6%)    x2Big=15 (1.1%)
TCAST    : 1328   199     29.3    1328    199     x2Tiny=10 (0.7%)    x2Big=12 (0.9%)
 ...  All groups ... 
Geneset  : nGene  Bits    dSize   rGene   rBits   outlierSize
AMELL    : 21920  82.6    76.8    9098    199     x2Tiny=162 (1.7%)       x2Big=481 (5.2%)
APISU    : 21920  74.5    51.6    8211    199     x2Tiny=191 (2.3%)       x2Big=485 (5.9%)
DMELA    : 21920  99.6    75.7    10973   199     x2Tiny=31 (0.2%)        x2Big=353 (3.2%)
DPULE    : 21920  68.6    -10.1   7559    199     x2Tiny=675 (8.9%)       x2Big=259 (3.4%)
HSALT    : 21920  85.8    41.7    9456    199     x2Tiny=505 (5.3%)       x2Big=380 (4%)
ISCAP    : 21920  63.1    -73.4   6948    199     x2Tiny=1203 (17.3%)     x2Big=184 (2.6%)
NVITR    : 21920  83.1    73      9149    199     x2Tiny=221 (2.4%)       x2Big=518 (5.6%)
PHUMA    : 21920  71.7    51.2    7901    199     x2Tiny=278 (3.5%)       x2Big=341 (4.3%)
TCAST    : 21920  77.7    49.6    8561    199     x2Tiny=207 (2.4%)       x2Big=464 (5.4%)

=item run7 : fish7 ordb set

$evigene/scripts/eval_orthodb_groups.pl -group OrthoDB6_Actinopterygii_tabtext.gz -geneaa fish-orthodb6.aa
.count
#i orthogroups from OrthoDB6_Actinopterygii_tabtext.gz
#i groups n=18537 with ntaxa>=3; ndrop=1841
#i groupgenes n=120084
#i species n=7
#i geneaa n=120084 nogeneaa=26252, sppcount DRERI:17809, GACUL:17674, GMORH:16319, OLATI:16322, ONILO:18479, TNIGR:16610, TRUBR:16871
Summary stats for OrthoDB6_Actinopterygii_tabtext.topordb
 ...  Common groups ... 
Geneset  : nGene  Bits    dSize   rGene   rBits   outlierSize
DRERI    : 9290   199     5       9290    199     x2Tiny=222 (2.3%)   x2Big=71 (0.7%)
GACUL    : 9290   199     -9.6    9290    199     x2Tiny=170 (1.8%)   x2Big=7 (0%)
GMORH    : 9290   199     -35.2   9290    199     x2Tiny=287 (3%)     x2Big=5 (0%)
OLATI    : 9290   199     -16.2   9290    199     x2Tiny=334 (3.5%)   x2Big=17 (0.1%)
ONILO    : 9290   199     24      9290    199     x2Tiny=17 (0.1%)    x2Big=66 (0.7%)
TNIGR    : 9290   199     -31.2   9290    199     x2Tiny=345 (3.7%)   x2Big=17 (0.1%)
TRUBR    : 9290   199     -5.3    9290    199     x2Tiny=170 (1.8%)   x2Big=35 (0.3%)
 ...  All groups ... 
Geneset  : nGene  Bits    dSize   rGene   rBits   outlierSize
DRERI    : 18537  165.4   23.4    15407   199     x2Tiny=310 (2%)     x2Big=354 (2.2%)
GACUL    : 18537  178.2   -6.3    16604   199     x2Tiny=307 (1.8%)   x2Big=30 (0.1%)
GMORH    : 18537  166.6   -34.5   15519   199     x2Tiny=616 (3.9%)   x2Big=9 (0%)
OLATI    : 18537  165.6   -12.4   15430   199     x2Tiny=550 (3.5%)   x2Big=52 (0.3%)
ONILO    : 18537  182.3   28.4    16977   199     x2Tiny=40 (0.2%)    x2Big=192 (1.1%)
TNIGR    : 18537  166.7   -23.9   15525   199     x2Tiny=516 (3.3%)   x2Big=62 (0.3%)
TRUBR    : 18537  170.3   -1.2    15862   199     x2Tiny=316 (1.9%)   x2Big=72 (0.4%)

=item run8 : fish7 ordb set + kfish

$evigene/scripts/eval_orthodb_groups.pl -group vert53kf-orthodb6_fake_tabfix -geneaa fish-orthodb6.aa.coun
t
#i orthogroups from vert53kf-orthodb6_fake_tabfix
#i groups n=21692 with ntaxa>=3; ndrop=3679
#i groupgenes n=888224
#i species n=53
#i geneaa n=143879 nogeneaa=46022, sppcount DRERI:20940, FHETE:14647, GACUL:18536, GMORH:17646, OLATI:17230, ONILO:19732, TNIGR:17585, TRUBR:17572
#i drop no-species-aasize n=45 : ACARO AMELA BTAUR CFAMI CHOFF CJACC CPORC DNOVE DORDI ECABA EEURO ETELF FCATU GGALL GGORI HSAPI LAFRI LCHAL MDOME MEUGE MGALL MLUCI MMULA MMURI MMUSC NLEUC OANAT OCUNI OGARN OPRIN PABEL PCAPE PTROG PVAMP RNORV SARAN SHARR SSCRO STRID TBELA TGUTT TSYRI TTRUN VPACO XTROP

Summary stats for vert53kf-orthodb6_fake_tabfix.topordb
 ...  Common groups ... 
Geneset  : nGene  Bits    dSize   rGene   rBits   outlierSize
DRERI    : 5821   199     12.2    5821    199     x2Tiny=125 (2.1%)   x2Big=64 (1%)
FHETE    : 5821   199     -11.4   5821    199     x2Tiny=297 (5.1%)   x2Big=41 (0.7%)  ?? too many tiny vs my own compute: MISSING Odb Funhe genes
GACUL    : 5821   199     0.6     5821    199     x2Tiny=108 (1.8%)   x2Big=16 (0.2%)
GMORH    : 5821   199     -20.4   5821    199     x2Tiny=153 (2.6%)   x2Big=12 (0.2%)
OLATI    : 5821   199     -5.7    5821    199     x2Tiny=200 (3.4%)   x2Big=22 (0.3%)
ONILO    : 5821   199     31.2    5821    199     x2Tiny=10 (0.1%)    x2Big=65 (1.1%)
TNIGR    : 5821   199     -16.3   5821    199     x2Tiny=182 (3.1%)   x2Big=24 (0.4%)
TRUBR    : 5821   199     6.9     5821    199     x2Tiny=87 (1.4%)    x2Big=38 (0.6%)
 ...  All groups ... 
Geneset  : nGene  Bits    dSize   rGene   rBits   outlierSize
DRERI    : 13597  178     26.9    12159   199     x2Tiny=286 (2.3%)   x2Big=319 (2.6%)
FHETE    : 13597  125.6   -1.2    8579    199     x2Tiny=464 (5.4%)   x2Big=157 (1.8%)
GACUL    : 13597  179     -3.4    12231   199     x2Tiny=308 (2.5%)   x2Big=44 (0.3%)
GMORH    : 13597  174.2   -28.7   11901   199     x2Tiny=561 (4.7%)   x2Big=23 (0.1%)
OLATI    : 13597  169.9   -11.1   11612   199     x2Tiny=537 (4.6%)   x2Big=63 (0.5%)
ONILO    : 13597  179.3   32.7    12249   199     x2Tiny=39 (0.3%)    x2Big=197 (1.6%)
TNIGR    : 13597  168.7   -19.4   11529   199     x2Tiny=448 (3.8%)   x2Big=65 (0.5%)
TRUBR    : 13597  169     1       11547   199     x2Tiny=316 (2.7%)   x2Big=80 (0.6%)


=cut


