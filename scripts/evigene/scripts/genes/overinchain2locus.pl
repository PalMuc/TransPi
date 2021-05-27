#!/usr/bin/env perl
# evigene/scripts/genes/overinchain2locus.pl

=item overintron chain.table to locus classed genes

 usage: overinchain2locus.pl .pl < mygenes.inchains > mygenes.intronloc.tab

 classify by longest, alternate, subchain for overlapped intron chains (loci), summarize inchains,
 input inchains table from overintron.pl -format scoreinidchain

  set pt=evg45merge3yx ; 
  $evigene/scripts/overintron.pl -exon exon -fixin -format scoreinidchain \
    -introns corn6rseq-v4genoasm.introns.gff -genes $pt-v4genoasm.gmap.gff.gz \
    > $pt-v4genoasm.inchains

=item FIXME bug in chain classing

  got separate inchain loci for chains containing same introns ..
  problem with inc{intronId}{ichainId} 
  .. maybe inId belongs to >1 ichain? dont replace ichainIds
  .. or/and hashmap must allow multiple td
  >> sort _sortichain needed to get longest per intron
  
=item ADDed alt trancript i column

  add ti column: ichain = t1 (longest chain, locus main), icdup = t1 dup only?
    icalt = t2,..,tn; icsub = t2sub,.. or is it dup of alt t2?; 

=item current input genes.ichains table

  Arath4EVm008353t10	ints=100,8/8,ic:c3i10052378r,c3i10052632r,c3i10052465r,c3i10052942r,c3i10052713r,c3i10053205r,c3i10053020r,c3i10053376r
  Arath4EVm013060t3	ints=100,4/4,ic:c1i26130730f,c1i26130807f,c1i26130869f,c1i26130961f
  Arath4EVm010132t8	ints=100,8/8,ic:c4i15017884r,c4i15018041r,c4i15017968r,c4i15018364r,c4i15018139r,c4i15018793r,c4i15018724r,c4i15019048r
  Arath4EVm010066t4	ints=100,8/8,ic:c1i28490873f,c1i28491171f,c1i28491301f,c1i28491391f,c1i28491624f,c1i28491744f,c1i28491999f,c1i28492105f

=item current output genes.ichainloc table

  * probably should merge ichain/icdup in 1 row, only 2nd+ geneids added col3, as with icalt/icsub 
  ichain	l4789.t1	AT2G30920.1	16i,
  icdup	l4789.t1	AT2G30920.2,Arath4EVm012153t1	0
  >> becomes
  ichain	l4789.t1	AT2G30920.1,AT2G30920.2,Arath4EVm012153t1	16i,
  
arath16ap_evg4cds.boinchainri.loc
ichain	l4789.t1	AT2G30920.1	16i,c2i13157528r,c2i13157792r,c2i13157712r,c2i13158054r,c2i13157905r,c2i13158231r,c2i13158182r,c2i13158399r,c2i13158351r,c2i13158763r,c2i13158686r,c2i13159088r,c2i13159010r,c2i13159248r,c2i13159159r,c2i13159553r
icdup	l4789.t1	AT2G30920.2,Arath4EVm012153t1	0

ichain	l4790.t1	AT3G53580.1	16i,c3i19865059f,c3i19865166f,c3i19865333f,c3i19865411f,c3i19865481f,c3i19865582f,c3i19865639f,c3i19865743f,c3i19865922f,c3i19866011f,c3i19866095f,c3i19866218f,c3i19866277f,c3i19866379f,c3i19866452f,c3i19866791f
icdup	l4790.t1	Arath4EVm010612t1,Arath4EVm010612t2	0
icsub	l4790.t1	Arath4EVm010612t3	8i,c3i19865059f,c3i19865166f,c3i19865333f,c3i19865411f,c3i19865481f,c3i19865582f,c3i19865639f,c3i19865743f

ichain	l4791.t1	AT1G51460.1	16i,c1i19077488r,c1i19077677r,c1i19077570r,c1i19078045r,c1i19077773r,c1i19078412r,c1i19078128r,c1i19078661r,c1i19078500r,c1i19078866r,c1i19078744r,c1i19079402r,c1i19078941r,c1i19079566r,c1i19079487r,c1i19081149r
icdup	l4791.t1	Arath4EVm002883t1,Arath4EVm002883t2	0
icalt	l4791.t2	Arath4EVm002883t3	12i,c1i19078045r,c1i19078412r,c1i19078128r,c1i19078661r,c1i19078500r,c1i19078866r,c1i19078744r,c1i19079402r,c1i19078941r,c1i19079566r,c1i19079487r,c1i19081149r
icalt	l4791.t3	Arath4EVm002883t5	12i,c1i19077488r,c1i19077677r,c1i19077570r,c1i19078045r,c1i19077773r,c1i19078412r,c1i19078128r,c1i19078866r,c1i19078500r,c1i19078744rk,c1i19078661rk,c1i19078941r
icalt	l4791.t4	Arath4EVm002883t6	10i,c1i19077488r,c1i19077677r,c1i19077570r,c1i19078045r,c1i19077773r,c1i19078412r,c1i19078128r,c1i19078661r,c1i19078500r,c1i19078744r
  example keptintron "k" : l4791.t3 c1i19078500r,[c1i19078744rk,c1i19078661rk],c1i19078941r 
  
=item problems from gene join/split errors

 this elevate mistakes of gene joins, split genes to longest chain 
  .. need to filter out problematic genes.gff, for creating inchain loci (a1),
     then add in problematic/remaining genes to existing inchain loci (a2)
  .. filter genes.gff: a1 no split genes _C1/2, or dup maps _G2.., 
      limit a1 to ~100% map coverage+identity, 
      limit a1 to longest transcripts t1..t9 ? (this may exclude some paralog misloci; longest uniq loc?)
  
=cut

use strict; 

my(%dic, %icn, %icl, %inc, %icd, %icloc, %icval, %locic, %dxc,
  %dicbad, %icnbad, %incbad);

my $debug=$ENV{debug}||0;
my $noBADGENES= $ENV{nobad}||0;     # a1. dont use flagged-bad gene models to create/extend ichains
my $addBADGENES= $ENV{addbad}||0;   # a2. input a1 ichain table + genes.gff, only add to existing ichain
my $BADPATT=$ENV{'badpatt'} || 'bad|poor';
my $MAX_SAMECHAIN= $ENV{maxsame}||5; ## putEqualTab: need to filter out excessive alts of same ichain 

sub _sortichain{ return $icn{$b}<=>$icn{$a} or $a cmp $b }

# a2: opt to read prior ichain.loctab of a1, 
# should read/retain %locid here?
# does this imply dont add new ichain below? all added below should be icdup/icalt/icsub?
# use new hashes? %dicIN %icnIN %incIN 

if(my $inf= $ENV{readloc}) { readInchainloc($inf); }

if($ENV{eqgene}) { putEqualTab(); exit; } # after readInchain, or readInchainloc

 ## update inchains iv?
 ## Arath5EVm001221t1	ints=100,4/4,ic:c5i21474427r,c5i21474512r,c5i21474703r,c5i21474826r
 ## >> 	ints=100,4/4[;\t]ic:i,i,i[;\t]xc:1:10,2:99,3:88

while(<>) {
  my @v=split; # only trans.id, intron.value compound field from overintron.pl  ; may change
  my $isbad= ($_ =~ m/$BADPATT/)?1:0; # match all of line, ID and IV; eg ID(_C[12]|_G\d)$ can be bad
  my($td,$iv,$ic,$xc);
  ($td,$iv)=@v;
  if(@v == 2) {
    ($ic)= $iv=~m/\bic:([^;]+)/; 
    ($xc)= $iv=~m/\bxc:([^;]+)/; 
  } else {
    ($ic,$xc)= @v[2,3]; $xc||=""; 
    $ic=~s/^ic://;  $xc=~s/^xc://; 
  }
  $ic=~s/,0//g; $ic=~s/^0[,]*//; 
  if($ic) { 
   my @in=split",",$ic; my $ni=@in; 
   if($noBADGENES and $isbad) { 
    $icnbad{$ic}=$ni; $dicbad{$ic}{$td}++; 
    map{ $incbad{$_}{$ic}=$td; } @in; 
   } else {
    $icn{$ic}=$ni; $dic{$ic}{$td}++; # ichainId of trid
    map{ $inc{$_}{$ic}=$td; } @in; # intronId > ichainId; maybe many td per inid/icid > inc{in}{ic}{td}=1 ?
   }
  }
  if($xc) {
    $dxc{$xc}{$td}++; #= x1-150,x2-99,x3-250,x4-500,x5-120 ? exon num-width
  }
}

# if($ENV{eqgene}) { putEqualTab(); exit; } # after readInchain, or readInchainloc

print_ichains();

#------------------------

sub altchainsOfChain {
  my($ic, $altchains)= @_;
  $altchains->{$ic}= 1;
  for my $i (split",", $ic) {
    my @oc= sort _sortichain grep{ not $altchains->{$_} } keys %{$inc{$i}};  
    map{ $altchains->{$_}= 1; } @oc;
    map{ altchainsOfChain($_, $altchains); } @oc; # recurse
    }
}

# convert intron locids to letters ABC..Zabc..z, local string rep of intron chains/inchain

use constant COMMINSTR  => 0; # TEST option; NOT useful; intron lets ordered by frequency among alts
use constant COMMINRUNS => 1; # TEST option; intron runs in common intron groups

sub chainsStringsOf {
  my($mainic,$otds)= @_;
  my(%istr, %lij, %ljc, %lic, %cin, %irun);
  
  my $isrev=($mainic =~ m/\d(r|rk)\b/)?1:0;
  my @oc= sort _sortichain keys %$otds;
  #? use insplice pairs (i.e. intron) as letters?
  for my $ic ($mainic,@oc) {
    for my $i (split",", $ic) { my($ib)= ($i=~m/i(\d+)/)?$1:$i; 
      $ljc{$ib}++; $lij{$ib}=$i; $lic{$ib}{$ic}++; }
    }
    
  my @in;
if(COMMINSTR) { 
  # test, or instead use loc order below, then compress isplices by common frequency?
  @in= sort{ $ljc{$b}<=>$ljc{$a} or $a<=>$b } keys %lij;
} else {  
  @in= sort{ $a<=>$b } keys %lij;
  # @in= reverse @in if($isrev); #? tried, not good
}
  
  my $nin0= $#in;
  my @krun = (0..$nin0);
  my @lrun = (1) x (1+$nin0);
  # my @frun = (1) x (1+$nin0);
if(COMMINRUNS) { #? not both COMMINSTR & COMMINRUNS
  my (@inc); my $krun=0;
  for(my $k=0; $k<=$nin0; $k++, $krun++) { 
    my $ib= $in[$k];    
    my $fi= $ljc{$ib};
    my $lic= join",", sort keys %{$lic{$ib}}; # need ichains of in to say if same run..
    $krun[$k]= $krun; $lrun[$krun]=1; # $frun[$k]= $fi;
    # my $j= $k; my $run=$ib;
    for(my $k2= $k+1; $k2 <= $nin0; $k2++) {
      my $jb=$in[$k2];
      my $ljc= join",", sort keys %{$lic{$jb}}; # need ichains of in to say if same run..
      #?? need ichains + freq here to say if run is same
      if($ljc{$jb} == $fi and $ljc eq $lic){ $k=$k2; $krun[$k]= $krun; $lrun[$krun]++; } else { last; } #$run.= ",$jb";  $frun[$k]= $fi;
      #x if($ljc{$jb} == $fi){ $k=$k2; $krun[$k]= $krun; $lrun[$krun]++; } else { last; } #$run.= ",$jb";  $frun[$k]= $fi;
      }
      
    #x push @inc, $ib; 
    # if($j>$k) { 
    #   my $ir=1+$j-$k; $irun{$ib}=$ir; 
    #   #x for(my $k2=$k+1; $k2<=$j; $k2++) { $krun[$k2]= $krun; }
    #   #? $runv{$run}= $ib; # diff  run patts same ib, length : need diff let
    #   $k=$j; 
    # }
  }
  #x @in= @inc; # keep @in as is, change char val/@in via krun[k]
}

  my (@frun);
  for(my $k=0; $k<=$nin0; $k++) {     
    ## chr(48+($k-52)) == 123; chr(97+($k-26)) == abc; chr(65+$k) == ABC
    #x my $c= ($k>61)? '_' : ($k>=52)? chr(48+($k-52)) : ($k>=26) ? chr(97+($k-26)) : chr(65+$k);
    my $kc= $krun[$k]; ##||$k; #<< bad for k/krun == 0 valid num
    my $c= ($kc>51)? '_' : ($kc>=26) ? chr(97+($kc-26)) : chr(65+$kc);
    my $ib= $in[$k];
    my $in= $lij{ $ib };
    my $fi= $ljc{ $ib }; my $lrun=$lrun[$kc];
    $frun[$kc]= "$c$kc"."f$fi"."l$lrun"; # $frun{$in}= "$c/$kc:$fi";
    #x $frun[$kc]= "$c:$kc:$fi:$lrun"; # $frun{$in}= "$c/$kc:$fi";
    #x if(COMMINRUNS) { if(my $ri= $irun{$ib}){ $c= join"", (($c) x $ri); } } #? maybe dont add run-count
    #x if(COMMINRUNS) { if(my $ri= $irun{$ib}){ $c="$c$ri"; } } #? maybe dont add run-count
    $cin{$in}= $c;
  }
  
  my $locusfreq= (($isrev)?"rev_":"fwd_"). join",",@frun; # for each ichain? or just locus sum?
  # my %fstr; $fstr{'locus'}= $locusfreq;
  
  for my $ic ($mainic,@oc) {
    my @cic;
    for my $in (split",", $ic) { 
      my $c= $cin{$in}; $c||="-"; 
      #x if(COMMINRUNS) { next unless($c); }  else { $c||="-"; } # this may be bad
      push @cic, $c;
      # push @fic, $frun{$in};
      }
    $istr{$ic}=  join"",@cic;
    #  $fstr{$ic}= join",",@fic;
    }
  return(\%istr, $locusfreq);
}


sub readInchainloc { 
  my($inf)=@_;
  my($nloc)= (0);
  open(F,$inf); 
  while(<F>){ 
  if(/^ic/) { # if(/^ichain/)  # only care of this one? or read all iclasses
    my($icl,$locid,$tds,$nic,$istr)=split; 
    my($ni,$ic)=split",",$nic,2; $ni=~s/i$//; 
    my @in=split",",$ic; $ni=@in if($ni < @in);
    $icn{$ic}=$ni; $icl{$ic}= $icl; $icloc{$ic}= $locid; 
    my($li,$ti)=split/\./,$locid; map{ s/^[lt]//; } ($li,$ti); 
    my $icval= join"\t",$icl,$locid,$tds,$ni,$ic,$istr;  # OOPS dup locid, diff ic
    $icval{$ic}= $icval;

    #? push @{ $locic{$li}{$ti} }, $ic; # many ic/locid
    $nloc++;
    for my $td (split",",$tds) { 
      $icd{$td}=$ic; $icloc{$td}=$locid; 
      $dic{$ic}{$td}++; map{ $inc{$_}{$ic}=$td; } @in; } 
    } 
  } close(F);
  return($nloc);
}

sub _sortevgloc {
  my($ag,$at)= $a=~m/(.+)[t\.](\d+)$/;
  my($bg,$bt)= $b=~m/(.+)[t\.](\d+)$/;
  return ($ag cmp $bg or $at <=> $bt or $a cmp $b);
}

sub putEqualTab {

  #for my $li (sort{$a <=> $b} keys %locic) {
  #  for my $ti (sort{$a <=> $b} keys %{$locic{$li}}) {
  #    my $icv= $locic{$li}{$ti};
  # ... } }
  
  ## need to filter out excessive alts of same ichain .. need other quals to be sure,
  ## but add limit of 2..9 max?  $MAX_SAMECHAIN= 5;
  
  my(%lig,%gil,%tds,@tds,@stds);
  for my $ic (sort keys %icval) {
    my $icv= $icval{$ic};    
    my($icl,$locid,$tds,$ni,$icxx,$istr)=split"\t",$icv;
    my($li,$ti)=split/\./,$locid; map{ s/^[lt]//; } ($li,$ti);
    my $cloc= "$locid.$icl"; my $it=0;     
    my @tdi= sort _sortevgloc split(",",$tds);
    for my $td (@tdi) {
      next if($td=~ m/$BADPATT/); # filter ids
      my $g=$td; $g=~s/_[CG]\d+$//; $g=~s/t\d+$//; 
      $tds{$td}= "$li.$ti"; 
      ++$it;  my $cloct= $cloc;
      if($it > $MAX_SAMECHAIN) {  $cloct.=".icdup"; }  
      push @tds, $td; # include icdup
      #x if($it > $MAX_SAMECHAIN) { push @stds, $td; $cloct.=".icdup"; } else {  push @tds, $td; }
      #x push @tds, $td unless($it > $MAX_SAMECHAIN); 
      $gil{$g}{$li}{$td}=$cloct;
      $lig{$li}{$g}{$td}=$cloct;
      #NOT here: $it++; last if($it>=$MAX_SAMECHAIN);
    } 
  } 

  ## keep stds in order w/ tds? and flag as icdup for followon filter? dont know other quals of icdups here
  for my $td (@tds) {
    my $g=$td; $g=~s/_[CG]\d+$//; $g=~s/t\d+$//; 
    my($li,$ti)=split/\./,$tds{$td};
    my @gil= grep{ $_ ne $li } sort keys %{$gil{$g}};
    my @lig= grep{ $_ ne $g } sort keys %{$lig{$li}};
    
    my($eqgenes,$maploc,$mapqual,$negenes)= ("") x 9;
    $maploc= $gil{$g}{$li}{$td}; #? or from map.attr
    $mapqual="noqual"; # read from map.attr 99a,99i,999l,14x, 
    # these are map quals, should move from maploc?
    #   l4630.t8.icsub == lower qual alt, keep only if other quals good
    #   l4630.t7.icalt.icdup == duplicate alt, keep only if other quals good
    
    my @og=(); for my $gil (@gil) { push @og, keys %{$gil{$g}{$gil}}; }
    $negenes = join",", sort @og;    
    @og=(); for my $lig (@lig) { push @og, keys %{$lig{$li}{$lig}}; }
    $eqgenes = join ",",  map{ "$_/98.ic" } sort @og;
    
    map{ $_||="na" } ($eqgenes, $negenes);
    print join("\t",$td,"noid",$eqgenes,$maploc,$mapqual,$negenes)."\n";
  }    
  
}

=item add mapqual to eqgene 

perl -ne '($td)=@v=split; if(/^AQueryID/){} 
elsif($v[1]=~/^\d+$/){ ($td,$cov,$pid,$nix,$cloc,$npa,$oid,$tag,$clen)=@v; $cloc{$td}=$cloc; 
map{ s/\..*//; } ($cov,$pid,$nix); $mapq=join",",$cov."a",$pid."i",$clen."l",$nix."x"; 
if($npa=~/^C\d:/){ @nsp=split",",$npa; $sp=pop(@nsp); $sl=unshift(@nsp);  $sp="$sp%,$sl"; $mapq.=",Spl:$sp"; } 
$mapq{$td}=$mapq; } else { if($mapq=$mapq{$td}) { s/\tnoqual/\t$mapq/; } 
if($cloc=$cloc{$td}) { $lo=$v[3]; s/\t$lo/\t$cloc,$lo/; }  print; } ' \
 evg5arath_cds.map.attr arath16apcds_evg5cds.berinchain.evgeqgene \
 > arath16apcds_evg5cds.berinchain.evgeqgene2

Arath5EVm006926t5	noid	Arath5EVm006376t1/98.ic,..	chr4:8098471-8100165:-,l4630.t7.icalt	100a,100i,1556l,5x	na
Arath5EVm006376t8	noid	Arath5EVm006926t1/98.ic,..	chr4:8098742-8102732:-,l4630.t8.icsub	93a,100i,1344l,5x	na

=cut

sub print_ichains {

my($iloc,$ialt,%locid,%otds, %sumc, %suma);
my @ic= sort _sortichain keys %icn; 
for my $ic (@ic) { 
  my ($tdc,@tdd)= grep{ not $locid{$_}} sort keys %{$dic{$ic}}; 
  next unless($tdc); 
  $iloc++; $ialt=1; 
  my $lid="l$iloc.t$ialt";
  map{ $locid{$_}=$lid } ($tdc,@tdd);  $suma{$ialt}++;
  my %otds=();   

  my $locchains="/$ialt:$ic/";
  my %ocall= ();
  altchainsOfChain($ic, \%ocall);
  my @oc= grep{$_ ne $ic} sort _sortichain keys %ocall;  
  
  if($addBADGENES) { 
    my %bocall= ();
    altchainsOfChain($ic, \%bocall, \%incbad);
    my @boc= grep{$_ ne $ic} sort _sortichain keys %bocall; 
    push @oc, @boc if(@boc);
  }
  
  my %oty=();
  for my $oc (@oc) { 
    my $isub= index($locchains,$oc);
    my $ot= ($isub>0) ? "sub" : "alt";
    ## idup test isn't useful, already have unique ichains: $idup= index($locchains,":$oc/");
    if($ot eq "alt") { 
      ++$ialt; $locchains .="$ialt:$oc/"; $lid="l$iloc.t$ialt"; $suma{$ialt}++;
    } else { 
      my $salt=$ialt; # subchain of which alt?
      my $j= rindex($locchains,":",$isub); my $k= rindex($locchains,"/",$j); 
      $salt= substr($locchains,$k+1,$j-$k-1) if($k>=0 and $j>$k); 
      $lid="l$iloc.t$salt"; 
      } 
    
    my @tds=grep{ not $locid{$_}} sort keys %{$dic{$oc}}; 
    if(@tds){ my $tds=join",",@tds; $otds{$oc}="$ot\t$lid\t$tds"; $oty{$ot}++; map{ $locid{$_}= $lid; } @tds; }
  }      
  
  my($strh,$locsum)=  chainsStringsOf($ic,\%otds);
  
  #__ print above before alts 
  my $ns= 1 + $ic =~ tr/,/,/;
  $oty{alt}++; # ichain/main
  my $salt= join",", map{ my $n=$oty{$_}; "ic$_:$n" } sort keys %oty;
  ##   $salt ||= "icalt:0"; # count ichain  as alt
     
  print "\n#icloc\t$locid{$tdc}\t$salt; imax:$ns; iruns:$locsum;\n";
  my $tdlist=$tdc; my $si= $strh->{$ic}||".";
  if(@tdd) { $tdlist=join",",$tdc,@tdd; }
  print join("\t","ichain",$locid{$tdc},$tdlist,$ns."i,".$ic,$si)."\n"; $sumc{"ichain"}++;
  #x print join("\t","ichain",$locid{$tdc},$tdc,$ns."i,".$ic)."\n"; $sumc{"ichain"}++;
  #x if(@tdd){ my $tda=join",",@tdd; print join("\t","icdup",$locid{$tdd[0]}, $tda, 0)."\n"; $sumc{"icdup"} += @tdd; } 

  #add?? $locic{$li}{$ti}= join"\t",$icl,$tds,$ni,$ic,$istr;  
  
  #__ print in alts loop
  for my $oc (sort _sortichain keys %otds) { 
    my($ot,$lid,$tds)= split"\t",$otds{$oc}; # tds can be list, all with same locid{}
    my $ns= 1 + $oc =~ tr/,/,/;  my $si=$strh->{$oc}||".";
    print join("\t","ic$ot",$lid,$tds,$ns."i,".$oc,$si)."\n";  $sumc{"ic$ot"}++; 
    } 
    
} 

# summary: class sum sumc,  alt hist? t1..t20 ?
print "#n_class:"; for my $c (qw(ichain icalt icsub icdup)){ print " $c=",$sumc{$c}||0; } print "\n";
print "#n_alts :"; for my $i (1..20){ print " t$i=",$suma{$i}||0; } print "\n";
} # print_ichains

__END__

=item altchainsOfChain change
  
  #? change here? yes .. result is nearly same as nonrecurs version, good.
  #   a. collect all @oc from @in, 
  #     ar. add also all @in/@oc ? recursively 
  #   b. then _sortichain @ocall, call alts/subs
  # ens16corn32sep-v4genoasm.inchains 
  # old: #n_class: ichain=23893  icalt=54788  icsub=16736 icdup=9969
  # nv1: #n_class: ichain=23800- icalt=54822+ icsub=16795+ icdup=9906
  # nv2: #n_class: ichain=23800- icalt=54822+ icsub=16795+ icdup=9906 .. no change in dup class

  # } else {  # before altchainsOfChain()
  # for my $i (split",", $ic) {  
  #   # add to locus all ichains containing any introns of main
  #   my @oc= grep{$_ ne $ic} sort _sortichain keys %{$inc{$i}};  
  #   for my $oc (@oc) { 
  #     # my $ot=( $locchains =~ m/$oc,/ )? "sub":"alt"; 
  #     # want subchain class for sub of other alt ?? yes?
  #     my $isub= index($locchains,$oc);
  #     my $ot=($isub>0) ? "sub" : "alt";
  #     if($ot eq "alt") { 
  #       ++$ialt; $locchains .="$ialt:$oc/"; $lid="l$iloc.t$ialt"; $suma{$ialt}++;
  #     } else { 
  #       #? reclass sub as dup if full match to locchains alt
  #       my $salt=$ialt; # subchain of which alt?
  #       my $j= rindex($locchains,":",$isub); my $k= rindex($locchains,"/",$j); 
  #       $salt= substr($locchains,$k+1,$j-$k) if($k>=0 and $j>$k); 
  #       $lid="l$iloc.t$salt"; 
  #       } 
  #     
  #     my @tds=grep{ not $locid{$_}} sort keys %{$dic{$oc}}; 
  #     if(@tds) { my $tds=join",",@tds; $otds{$oc}="$ot\t$lid\t$tds";  map{ $locid{$_}= $lid } @tds; }
  #   }  
  # } 
  # }

=item joins detect 

  problem if intron/splice sites are mixed up across loci, e.g. w/ some gene joins
  detect joins? i.e. 1/few long chains join 2 distinct groups of shorter chains, split midway
    .. one way: gene-locus + altnum counts per ichain, 
   if > 1 locus/ichain check for large diff in intron counts at top (often evg t99 type), 
   vs t1/2/3 of each locus
    .. work with mixed gene set data: AT0000.iii and EVm0000tiiii
    .. but some cases may be true longer gene of fragment loci
   
  eg. 
ichain  l39.t1  Arath5EVm000492t23      80i  << join of g1(40i),g2(30i),g3(10i) below
icalt   l39.t2  AT1G01790.2     40i         << ATg1
icalt   l39.t3  AT1G01790.1,Arath5EVm000492t1   40i << ATg1, true EVg1
icalt   l39.t4  Arath5EVm000492t24      39i   EVg1
icalt   l39.t5  Arath5EVm000492t25      39i   EVg1
icalt   l39.t6  Arath5EVm000492t22      38i   EVg1
icalt   l39.t7  AT1G01770.1     32i           << ATg2
icsub   l39.t7  Arath5EVm003856t1       30i   << EVg2
icalt   l39.t8  Arath5EVm000492t26      24i
icalt   l39.t9  Arath5EVm003856t2       15i   EVg2 
icalt   l39.t10 AT1G01780.1     10i     << ATg3
icalt   l39.t11 Arath5EVm017606t2       8i  << EVg3
icalt   l39.t12 Arath5EVm017606t1       8i  << EVg3
--

ichain  l18.t1  Arath5EVm000800t66      94i   << join of several tandem genes
icalt   l18.t2  Arath5EVm000800t154     94i   << join ditto
icalt   l18.t3  Arath5EVm000800t48      48i
icalt   l18.t4  Arath5EVm000800t30      48i
icalt   l18.t6  Arath5EVm000800t46      48i
icalt   l18.t5  Arath5EVm000800t47      48i
icalt   l18.t7  Arath5EVm000800t83      48i
icalt   l18.t8  AT1G56130.2     48i         << ATg1
icalt   l18.t9  Arath5EVm000800t108     48i
icalt   l18.t10 Arath5EVm000800t74      48i
icalt   l18.t11 AT1G56140.1,Arath5EVm000800t23,Arath5EVm000800t24       48i  << ATg2
icalt   l18.t12 Arath5EVm000800t44      48i
icalt   l18.t13 Arath5EVm000800t25,Arath5EVm000800t28,Arath5EVm000800t29,Arath5EVm000800t34     48i
icalt   l18.t14 Arath5EVm000800t51,Arath5EVm000800t52,Arath5EVm000800t54        48i
icalt   l18.t15 Arath5EVm000800t2       48i
icalt   l18.t16 Arath5EVm000800t50,Arath5EVm000800t53   48i
icalt   l18.t17 Arath5EVm000800t43,Arath5EVm000800t45,Arath5EVm000800t49        48i
icalt   l18.t18 Arath5EVm000887t10      46i
icalt   l18.t19 Arath5EVm000887t5       46i
icalt   l18.t20 Arath5EVm000887t6       46i
icsub   l18.t8  AT1G56130.1,Arath5EVm000800t76  46i  << ATg1
icalt   l18.t21 AT1G56120.1     46i           << ATg3
icalt   l18.t22 Arath5EVm000800t59      46i
icalt   l18.t23 AT1G56120.2,Arath5EVm000800t35,Arath5EVm000800t36,Arath5EVm000800t37,Arath5EVm000800t38,Arath5EVm000800t39,Arath5EVm000800t40,Arath5EVm000800t41,Arath5EVm000800t42,Arath5EVm000800t58,Arath5EVm000800t60,Arath5EVm000800t61,Arath5EVm000800t64,Arath5EVm000800t65      46i
icalt   l18.t24 AT1G56145.2     46i         << ATg4 ?
icalt   l18.t25 AT1G56120.3,Arath5EVm000800t72,Arath5EVm000800t73,Arath5EVm000800t77,Arath5EVm000800t78 46i
icalt   l18.t26 Arath5EVm000800t3       46i
icalt   l18.t27 AT1G56145.1,Arath5EVm000887t2,Arath5EVm000887t3 46i  << ATg4 
icalt   l18.t28 Arath5EVm000800t134     45i
icalt   l18.t29 Arath5EVm000887t1       45i
icalt   l18.t30 Arath5EVm000800t87      45i
icalt   l18.t31 Arath5EVm000800t110     45i
icalt   l18.t32 Arath5EVm000800t84      45i
icsub   l18.t27 AT1G56145.3     44i       << ATg4 ?
icalt   l18.t33 Arath5EVm000887t7       44i
icalt   l18.t34 Arath5EVm000887t12      44i
icalt   l18.t35 Arath5EVm000887t8       44i
icalt   l18.t36 Arath5EVm000800t106     42i
icalt   l18.t37 Arath5EVm000887t4       42i
icsub   l18.t8  AT1G56130.3     40i       << ATg1
icalt   l18.t38 Arath5EVm000800t79      40i
icalt   l18.t39 Arath5EVm000800t88      38i
icalt   l18.t40 Arath5EVm000800t89      38i
icalt   l18.t41 Arath5EVm000800t143     32i
icalt   l18.t42 Arath5EVm000800t153     30i
icalt   l18.t43 Arath5EVm000887t9       28i
icalt   l18.t44 Arath5EVm000800t158     21i
icalt   l18.t45 Arath5EVm000800t151     16i
icalt   l18.t46 Arath5EVm000800t156     10i


--  
ichain  l42.t1  Arath5EVm002117t4       78i   << join of 3? loci, too long ichain
icalt   l42.t2  Arath5EVm002117t6       76i   << join ditto
icalt   l42.t3  Arath5EVm005931t3       40i   << diff locus
icalt   l42.t4  Arath5EVm002117t2       36i   << likely true ichain
icalt   l42.t5  AT4G21530.1,Arath5EVm002117t1   36i  << ditto, ATg1
icsub   l42.t5  Arath5EVm002117t3       30i
icalt   l42.t6  Arath5EVm002117t5       24i
icalt   l42.t7  Arath5EVm005931t12,Arath5EVm005931t13   22i
icalt   l42.t8  AT4G21534.1,Arath5EVm005931t5,Arath5EVm005931t6 22i  << ATg2
icalt   l42.t9  Arath5EVm002117t7       21i
icalt   l42.t10 AT4G21534.2,Arath5EVm005931t1,Arath5EVm005931t7 20i
icsub   l42.t3  AT4G21540.1,Arath5EVm007082t1,Arath5EVm007082t2 18i  << ATg3
icalt   l42.t11 Arath5EVm005931t2       18i
icalt   l42.t12 AT4G21540.3     18i
   
=cut
  
=item corn genes ichain locus summary

==> cshlrnapb-v4genoasm.inchain.locd <==
#n_class: ichain=21468 icalt=62719 icsub=33617 icdup=170929
#n_alts : t1=21468 t2=12204 t3=8616 t4=6420 t5=4909 t6=3864 t7=3133 t8=2593 t9=2191 t10=1894 t11=1601 t12=1378 t13=1237 t14=1060 t15=939 t16=839 t17=732 t18=659 t19=592 t20=527

==> ens16corn32sep-v4genoasm.inchain.locd <==
#n_class: ichain=23800 icalt=54822 icsub=16795 icdup=9906
#n_alts : t1=23800 t2=11268 t3=7740 t4=5673 t5=4373 t6=3443 t7=2798 t8=2296 t9=1928 t10=1635 t11=1387 t12=1202 t13=1051 t14=896 t15=794 t16=687 t17=613 t18=544 t19=485 t20=440

==> evg45merge3yx-v4genoasm.inchain.locd <==
#n_class: ichain=25031 icalt=40578 icsub=23648 icdup=37491
#n_alts : t1=25031 t2=12925 t3=8572 t4=5716 t5=3951 t6=2753 t7=1912 t8=1347 t9=950 t10=654 t11=468 t12=341 t13=239 t14=174 t15=126 t16=90 t17=68 t18=57 t19=42 t20=35

==> maize_jgi14denovo-v4genoasm.inchain.locd <==
#n_class: ichain=25041 icalt=28406 icsub=13727 icdup=1360
#n_alts : t1=25041 t2=11120 t3=6515 t4=3943 t5=2407 t6=1504 t7=1008 t8=607 t9=398 t10=266 t11=186 t12=137 t13=101 t14=67 t15=38 t16=26 t17=21 t18=14 t19=8 t20=5

=item ichain loci of two gene sets

  cat {evg45merge3yx,ens16corn32sep}-v4genoasm.inchains | overinchain2locus.pl \
     > bothevg5ens6.inchain.loci 
#n_class: ichain=27074 icalt=92363 icsub=43798 icdup=55367
#n_alts : t1=27074 t2=16080 t3=12377 t4=9641 t5=7842 t6=6342 t7=5213 t8=4352 t9=3687 t10=3119 t11=2640 t12=2264 t13=1944 t14=1684 t15=1459 t16=1278 t17=1114 t18=971 t19=866 t20=770

ichain	l1.t1	Zeamay5fEVm000018t1	156i
icalt	l1.t2	Zm00001d040166_T012	148i
icalt	l1.t3	Zm00001d040166_T014	148i

=cut

=item arath case

  pt=evg1arathmrna_arath11ga
 
 $evigene/scripts/overintron.pl -exon exon -fixin -format scoreinidchain \
   -introns $weed/aweed16/arath11arap.introns.gff \
   -genes evg1arath.mrna.gmap.gff > $pt.inchains

  $evigene/scripts/genes/overinchain2locus.pl < $pt.inchains > $pt.inchain.loctab
  cat $pt.inchain.loctab | cut -f1 | sort | uniq -c | head
  16804 ichain = locus
  15657 icalt
   8611 icsub
   3815 icdup

  # simple alts/locus count distribution table
  egrep '^(ichain|icalt)' $pt.inchain.loctab | egrep -v '_G' | \
  cut -f2 | sort | uniq -c | sed 's/ l.*/na/;' | sort | uniq -c | head 

 10089    1na
  2980    2na
  1652    3na
   862    4na
   467    5na
   313    6na
   166    7na
    99    8na
    45    9na
    38   10na

        c2/c    ln(c2/c) count
1.alts  na      0.00    10089
2.alts  1       0.00    2980
3.alts  1.8     0.59    1652
4.alts  3.5     1.24    862
5.alts  6.4     1.85    467
6.alts  9.5     2.25    313
7.alts  18.0    2.89    166
8.alts  30.1    3.40    99
9.alts  66.2    4.19    45
10.alts 78.4    4.36    38

 egrep '^(ichain|icalt)' $pt.inchain.loctab | egrep -v '_G' |   cut -f2 | sort | uniq -c | sed 's/ l.*/ alts/;' | sort | uniq -c | perl -ne '($c,$n)=split; $p="na"; if($n==1) { $c1=$c; $p="na"; } elsif($n==2) { $c2=$c; $p=1; } elsif($n>2 and $c>0) { $p=sprintf "%.1f",$c2/$c; } $lp= sprintf"%.2f",($c2>0 and $c>0)?log($c2/$c):"na"; print "$n.alts\t$p\t$lp\t$c\n";' | head 

=cut
