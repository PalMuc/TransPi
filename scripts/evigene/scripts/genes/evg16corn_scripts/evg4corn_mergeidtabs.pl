#!/usr/bin/env perl
# evg4corn_mergeidtabs.pl

=item about evg mergeidtabs
  
  evigene update to merge intronchain locus class table w/ prior publicset/pubids locus table
  -- tries to make sensible new locus groupings from both sorts of locus align evidence
  -- intronchain table now from pubset.cds x chrasm blastn locations
  -- output.trclass table then is input to trclass2mainalt to make new pubids table, w/ new locus/alt groupings
  
  evigene/genes/inexchains_gcb.pl -showexon 1 -pubidtab evg4corn2g.pubids -MODLOCUS 5 -pMINLOW 95.0 \
    -debug  -sizes evg4corn2g.cds.qual evg4corn2g_cds-genoasm.mblastn.gz \
    > evg4corn2g_cds-genoasm.inex9chain.tab 
    
  env swapids=1  newloc=1 format=trclass \
   evigene/genes/evg4corn_mergeidtabs.pl  evg4corn2g.pubids evg4corn2g_cds-genoasm.inex9chain.tab \
      >& evg4corn2g_inex9chain.trclass

  evigene/scripts/prot/trclass2mainalt.pl -idpre Zeamay4cEVm -debug -trclass evg4corn2g_inex9chain.trclass

  -- format=pubids not final, need trclass2mainalt process

  -- test case
  env nochr=0 ./mergeidtabs.pl evgm3chra-evg4corn.tgclass.pubids evg4cornmrna-genoasm.inv6chainct.tab \
     > & evg4corn.tgclasspub_inv6chainmr.idxtab

=item bugs in Ig locus classing
  
  -- is splitting true loci into 2+, mainly for  shortish alts with uncertain intron/exon overlaps
  -- need to repair, cancel "relocmiss" and relocadd,  given some added data such as GFF cds overlaps
       
=cut

use strict;

my(%aaref,%tqual,%chrmap,%chrmapi,%ptd,%tdp,%tdig,%igtd,%gids,%gnv,%partof);
my(%igt,%igg,%tig,%gig,%tdv,%newloc, %eqcds); my $neqcds=0;
my($keepid,$dropid,$nkeepid,$ndropid)= (0) x 9;

my $NOCHRBE= $ENV{nochr}||0;
my $IDPRE= $ENV{idpre}||"Zeamay4EVm"; # or "Zeamay4bEVm"; ## or ..
my $SWAPIDS= $ENV{swapids}||0;
my $IGONLYSET= $ENV{igonly}||0;
my $IGMISS='Ig999999'; ## FIXME: Ig000000 also is miss, nochr loc,  no introns ..
my $IGZERO='Ig000000'; ##my $IGMISSPAT='Ig999999|Ig000000'; #? or change Ig000000 to Ig999999 ?
my $OUTFORM=$ENV{format}||0; # for newloc only? format == 0|pubid|1|trclass|2|..
# my $CDSOVERTAB= $ENV{cdsover}||"";

sub readIdList{ 
 my($intab,$isdrop)=@_; 
 my $nid=0; my %ids=(); 
 open(F,$intab) or die "#ERR: idlist $intab \n"; 
 while(<F>) { next if(/^\W/); my($id,$oid)=split; $nid++; $ids{$oid}=1 if($oid); $ids{$id}=1; } close(F); 
 return ($nid,\%ids); 
}

if($ENV{keepid}) { ($nkeepid,$keepid)= readIdList($ENV{keepid},0); }
##later# if($ENV{dropid}) { ($ndropid,$dropid)= readIdList($ENV{dropid},1); }

if($ENV{diff}) { 
  chaindiff(); 

} elsif($ENV{newloc}) {
  mergeidtabs(); # == 2 input tables, pubids 1st
  # readPubids(); 
  # readIntronchains(); 
  
  #>> change readCdsOverTab to use eqgene tab, eg from runtr2genome3.sh > geno1f/evgm3chra-evg4corn.eqgene
  #    evigene/scripts/prot/overeqcdsloc.pl $eqopts -in $cdsbltab -out $cdseqgene
  #>> this is/should be same cds x chr mblastn as for Ig locus classing, merge these methods?
  
  $neqcds= ($ENV{cdsover})? readCdsOverTab($ENV{cdsover}, \%eqcds) : 0;
  
  my(%tdefer,%tdid);
  my($ngid,$ntid,$nthis,$ntadd,$ndef,$ndleft,$ncover)=(0) x 9;
  
  ## Zeamay4EVm000016t15utrorf <<< problem w/ id parse...
  
  for my $gid (sort keys %gids) { 
    my @tid= grep{ not $tdid{$_} } sort _gidsort keys %{$gids{$gid}};
    if($IGONLYSET) { @tid= grep { $tdig{$_} } @tid; }
    $ngid++; $ntid+= @tid;
    next unless(@tid); # done?
    
    ## BUG? no IGMISS in lists .. want to keep assigned to orig gid + top? Ig locus
    ## .. output bug?
    
    my(@tidigthis, @tidigo);
    my($igthis, @iggo)= grep { $_ ne $IGMISS } sort keys %{$tig{$gid}}; 
      #o: sort{ $tig{$gid}{$b}<=>$tig{$gid}{$a} or $a cmp $b} keys %{$tig{$gid}}; 
      #^ sort by Ig0000 instead of most Ig hits? and/or dont count td.partof for tig{gd}{ig}++ score
      
    if(@iggo) { # more than 1 ig/gid; noig cases not handled in 1st way? or not $tdig{$td} ?
      #x @tidigthis= grep{ $igtd{$igthis}{$_} or $igtd{$IGMISS}{$_} or not $tdig{$_} } @tid;
      #x @tidigo= grep{ not( $igtd{$igthis}{$_} or $igtd{$IGMISS}{$_} or not $tdig{$_} ) } @tid;
      ## other way:
      for my $td (@tid) { my $tdid=0;
	      for my $ig (@iggo) { if($igtd{$ig}{$td}) { push @tidigo, $td; $tdid=$ig; last; } }
        push @tidigthis, $td unless($tdid);
      }
    } else {
      @tidigthis= @tid; @tidigo=();
    }
    my @tidadd= grep { not (m/^$gid/ or $tdid{$_}) } sort _gidbsort keys %{$igtd{$igthis}};
    
    ## FIXUP for cdsoverlap table to cancel @tidadd and @tidigo Ig mistakes
    # main: %eqcds; $eqcds{$da}{$db}=$eqcds{$db}{$da}=xxx
    # ?? not working right yet ??
    if($neqcds) { ## %eqcds
      my(%dropadd,%dropmiss); my($nda,$ndm)=(0,0);
      for my $td (@tidadd) {
        for my $sd (@tidigthis) { if($eqcds{$td}{$sd}) { $dropadd{$td}=$sd; $nda++; last; } }
      }
      for my $td (@tidigo) {
        for my $sd (@tidigthis) { if($eqcds{$td}{$sd}) { $dropmiss{$td}=$sd; $ndm++; last; } }
      }
      $ncover += $nda + $ndm;
      if($nda) { @tidadd = grep{ not $dropadd{$_} } @tidadd; } # dropadd returns to next tidigthis via dropmiss
      if($ndm) { 
        my @tidm= grep{ $dropmiss{$_} } @tidigo; push @tidigthis, @tidm;
        @tidigo = grep{ not $dropmiss{$_} } @tidigo; 
        }
    }
    
    map{ $tdefer{$_}++ } @tidigo; # TOO many hits here, all of tid it seems

    $nthis+= @tidigthis; $ntadd+= @tidadd;
    ## shouldnt need now:  unless($tdid{$tid})
    
    ## order this,add may not be best alt order, resort by tid quals?
if(1) {
    my @tdord; my @act; 
    if(@tidadd) {  # is this reorder desired?
      my @uact = ("same") x scalar(@tidigthis);
      push @uact, ("add") x scalar(@tidadd);
      my @utdord= (@tidigthis,@tidadd);
      my @ord= sort { my $da=$utdord[$a]; my $db=$utdord[$b]; 
         $tqual{$db}{score} <=> $tqual{$da}{score} or $da cmp $db } (0..$#utdord);
      for my $i (0..$#ord) { $tdord[$i]= $utdord[$ord[$i]]; $act[$i]= $uact[$ord[$i]]; }
    } else { 
      @act=  ("same") x scalar(@tidigthis);  @tdord=@tidigthis; 
    }

     for my $i (0..$#tdord) {  my $tid=$tdord[$i];
      outnewlocgn($act[$i],$igthis,$gid,$tid) unless($tdid{$tid}); $tdid{$tid}++;
      $tdefer{$tid}=0 if($tdefer{$tid}); }
   
} else {
    for my $tid (@tidigthis) { outnewlocgn("same",$igthis,$gid,$tid) unless($tdid{$tid}); $tdid{$tid}++; 
      $tdefer{$tid}=0 if($tdefer{$tid}); }
    for my $tid (@tidadd) { outnewlocgn("add",$igthis,$gid,$tid) unless($tdid{$tid});  $tdid{$tid}++; 
      $tdefer{$tid}=0 if($tdefer{$tid}); }
}
  }
  
  $ndef= scalar(keys %tdefer); #?? too many
  my (%missloc);
  for my $tid ( grep{ $tdefer{$_}>0 } sort _gidbsort keys %tdefer) { # by gid or tid?
    my $tdig= $tdig{$tid}||"noig";
    my($gid)=map{ my($g,$t)= m/^(\w+)t(\d+)/; $g; } ($tid); # ** utrorf tag bug..
    my $missno= $missloc{$gid}{$tdig}||0; # for >2 Ig splits/gid locus
    my $missnew=0;
    unless($missno){ $missnew=1; my $mg= scalar(keys %{$missloc{$gid}}); $missloc{$gid}{$tdig}=$missno=1+$mg; }
    outnewlocgn("miss$missno",$tdig,$gid,$tid,$missnew) unless($tdid{$tid});  $tdid{$tid}++; $ndleft++; 
  }
  print "#sum.newloc ngid=$ngid, ntid=$ntid, nthis=$nthis, ntadd=$ntadd, ndefer=$ndef, ndefleft=$ndleft, ncdsover=$ncover\n";

} else { # merged pubid/ichainid table
  mergeidtabs(); # == 2 input tables, pubids 1st
  # readPubids(); 
  # readIntronchains(); 

  for my $gid (sort keys %gids) { 
    my @tid= sort _gidsort keys %{$gids{$gid}};
    if($IGONLYSET) { @tid= grep { $tdig{$_} } @tid; }
    for my $tid (@tid) { outgn($gid,$tid); }
  }
}


#--------------------------------

sub _gidsort { my($ta)= $a=~m/t(\d+)/; my($tb)= $b=~m/t(\d+)/; return($ta <=> $tb or $a cmp $b); } 
sub _gidbsort { my($ga,$ta)= $a=~m/(\w+)t(\d+)/; my($gb,$tb)= $b=~m/(\w+)t(\d+)/; return($ga cmp $gb or $ta <=> $tb); } 
# sub old_sortgid { my($ta)= $$a[1]=~m/(\d+)$/; my($tb)= $$b[1]=~m/(\d+)$/; return($ta <=> $tb or $a cmp $b); } 

sub  readCdsOverTab {  # ($ENV{cdsover}, \%eqcds);
  my($cdsoverf,$eqcdsh)=@_;
  my $neq=0;
  open(F,$cdsoverf) or return 0;
  while(<F>) {
    next unless(/^\w/);
    my($aid,$bid,$eqv)= split; $eqv||=1; #??
    $eqcds{$aid}{$bid}= $eqcds{$bid}{$aid}= $eqv; $neq++;
  } close(F);
  return $neq;
}

sub mergeidtabs {
  while(<>) {
    next if(/^\W|^Query/); 
    chomp; my @v=split"\t"; 
    
     # FIXME: add keepid/dropid hash lists.

     # opt to swap pd,td of pubids, in %gids, .., for locus calls depends which is best
     # need opt to report only those w/ inchains entries, test subsets, not all of pubids : gids? tdig*
    if(@v==8) { # evg.pubids tab, input 1st
      my($pd,$td,$pgORIG,$ti,$cla,$aaw,$pid,$flag)=@v;  
      
      if($SWAPIDS) { ($pd,$td)=($td,$pd); }
      next if(ref $keepid and not ($keepid->{$pd} or $keepid->{$td})); # also dropid??
      $ptd{$pd}=$td; $tdp{$td}=$pd; # both?
      my($tg,$pg)=map{ my($g,$t)= m/^(\w+)t(\d+)/; $g; } ($td,$pd); # ** utrorf tag bug..
      
      $gids{$tg}{$td}++; 
      $gnv{$td}=[$pd,$cla,$aaw,$pid,$flag];
      my($cloc)= ($flag=~m/chrmap:([^;\s]+)/)?$1:"";  $cloc=~s/,pflag.*//; $cloc=~s/,\d+i,/,/;
      $chrmap{$td}=$cloc;
      my($aaref)= ($flag=~m/aaref:([^;\s]+)/)?$1:""; $aaref=~s/,chrmap.*//; $aaref=~s/,pflag.*//;
      $aaref=~s/,oldid:\w+//;
      $aaref{$td}=$aaref;

      ## add tqual score for sort of new locus alts
      my($vref)= $aaref=~m/(\d+)/;  my($vaa)= $aaw=~m/(\d+)/; my $vscore=$vaa+$vref;
      $tqual{$td}{aaref}= $vref; $tqual{$td}{aalen}= $vaa;
      $tqual{$td}{score}= $vscore; # wtsum(aaref,aalen,align,partof) .. use for sort?

      }
      
    elsif(@v>10) { # intronchain.tab, input 2nd
      #** FIXME new ichain.tab MAY have rexons field before flag field
      my($td,$rc,$bits,$idn,$aln,$tw,$rw,$nx,$ig,$tbex,$rbe,$rexons,$flag)= (0) x 19;
      if(@v>12) {
        ($td,$rc,$bits,$idn,$aln,$tw,$rw,$nx,$ig,$tbex,$rbe,$rexons,$flag)= @v;
      } else {
        ($td,$rc,$bits,$idn,$aln,$tw,$rw,$nx,$ig,$tbex,$rbe,$flag)= @v;
      }
	
      $ig=$IGMISS if($ig eq $IGZERO); # fixup for Ig000000

      my $pd=0; 
      if($pd= $ptd{$td}) { ($td,$pd)=($pd,$td); }
      elsif($pd= $tdp{$td}) { }
      else { $pd=$td; } # error, not in pubids
      next if(ref $keepid and not ($keepid->{$pd} or $keepid->{$td})); # also dropid??

      #?NOT: if($SWAPIDS) { ($pd,$td)=($td,$pd); } # already swapped vi ptd/tpd
       # $td=$ptd{$td}||$td; ## evg4corn2g: td => pd of above pubids, not orig id
  
      my($tg,$pg)=map{ my($g,$t)= m/^(\w+)t(\d+)/; $g; } ($td,$pd); # ** utrorf tag bug..
       
      $gids{$tg}{$td}++; 
      $tdig{$td}=$ig; # other way not same: $tig{$tg}{$ig}++;
      $igtd{$ig}{$td}++; # or push @$igtd{$ig}, $td

      $igt{$ig}{$tg}++; $igg{$ig}{$pg}++; 
      $tig{$tg}{$ig}++; $gig{$pg}{$ig}++; 
      
      my($tbe,$txn)=split"/",$tbex;  
      $rbe=~s,/.*,,; $rbe="$rc:$rbe" unless($rbe=~m/^$rc:/); # now has $rc: prefix
      my($tb,$te)=split"-",$tbe; my $twx=1+$te-$tb; 
      $tw=$twx if($twx>$tw); $tw||=1;
      my $pal=int(0.5+100*$aln/$tw); $pal=100 if($pal>100); 
      my $partof=0; if($flag=~m/(partof=\w+)/) { $partof{$td}= $partof= $1;  }
      my $spl= ($flag=~m/split=(\d+)/)? ",Spl:$1":""; # split=2,NC_024459.1,NC_024467.1
       # chrmap:87a,100i,4627l,6x,NC_024465.1:70661022-70668846:+
      $chrmapi{$td}="${pal}a,${tw}l,${nx}x$spl,$rbe"; 

      $tqual{$td}{align}= $pal; #? $tqual{$td}{partof}= $partof; # split ?
      ## $tqual{$td}{score}= wtsum(aaref,aalen,align,partof) .. use for sort?
    } 
  }
}


sub outnewlocgn {
  my($act,$igthis,$gd,$td,$missnew)=@_; # (act==add|miss,$igthis,$gid,$tid);
  my($pd,$cla,$aaw,$pid,$onotes)= (0,"nocla",0,0,0);
  my $greloc= ($act=~/^(add|miss)/)?$act:""; # ?? replace add w/ qual sorted merge of same/add ??
  my $gdnew= ($greloc=~/miss/)? $gd.$greloc : $gd; #not for add

  ## FIXME: add gnv aaw,pid to output for .trclass/pubids equiv info
  ## fixme: ?? pid has main ID tag, maybe messing up trmainalt proc
  if( my $gnv=$gnv{$td} ) { ($pd,$cla,$aaw,$pid,$onotes)= @$gnv; }
  my $idpre=substr($td,0,3); $pid=~ s,/$idpre.*,,;
  
  my $chrmap= $chrmapi{$td}||$chrmap{$td}||"nochrmap"; 
  my $aaref=  $aaref{$td}||"noref"; $aaref=~s/,$//;
  my $partof= $partof{$td}||""; 
  my $ig= $tdig{$td}||"noig"; ## $igthis; # or Ig999999 missing  ; add noig flag

  # my $tinum= ++$newloc{$gdnew}; # need ti main/alt opts from qual
  my ($maintd,$tinum,$didmain)=($td,0,0);
  if($newloc{$gdnew}) {
    $tinum= ++$newloc{$gdnew}{ti}; # need ti main/alt opts from qual
    $maintd= $newloc{$gdnew}{main}; 
    #x $didmain=1 if($act =~ /add/ and $cla =~ /main/);
    $didmain=1 if($newloc{$gdnew}{havemain});
  } else {
    $newloc{$gdnew}{ti}= $tinum= 1;
    $newloc{$gdnew}{main}= $maintd= $td; # need ti main/alt opts from qual
    $newloc{$gdnew}{havemain}=0; # this is from $cla
    # if($missnew or here & act == miss) .. fiddle class to main, maintd ok?
    if($greloc =~ /miss/) { $cla =~ s/^\w+/main/; } #?
  }
  
  #o# my $ismain=($cla =~ /^main/ and not ($act =~ /add/))?1:0; # also noclass?
  my $ismain=(($cla =~ /^main/) and not $didmain)?1:0; # also noclass?
  ## ismain ||=  $td eq $maintd << not here, need also hange class if maintd == alt?
  $cla.=",reloc$greloc" if($greloc);
  if($partof and not($td eq $maintd) ) { 
     # FIXME: temp dont drop main partof, until have way to reclass of-id as new main
     # OR dont drop if td eq maintd 
     $cla.=",$partof"; $cla="drop$cla" unless($ismain or $cla=~/^(drop|cull)/); 
   }
  $newloc{$gdnew}{havemain}=$td if($ismain and not $newloc{$gdnew}{havemain});
  
  if($OUTFORM =~ /^[12]|pubid|trclass/) { ## $OUTFORM =~ /pubid/ .. add trclass format
    #pubids format: move igthis after chrmap? pack Notes field = "aaref:$aaref,chrmap:$chrmap,other"
    #Public_mRNA_ID	originalID	PublicGeneID	AltNum	Class	AAqual	pIdAln	Notes
    #Zeamay4EVm000002t1	okay	main	Zeamay4EVm000002t2	100/100/./altmap70xeq	5066,97%,complete	
    #    aaref:9468,Sobic.007G098500.1.p,refgood,chrmap:70a,100i,15201l,14x,NC_024465.1:70639407-70669113:+,pflag:0 
    ## trclass:
    # Zeamay4EVm000001t1	okay	main	Zeamay4EVm000001t2	100/87/.	5425,97%,complete	aaref:10109,Sobic.009G173600.1.p,refbest,chrmap:97a,100i,16278l,58x,NC_024464.1:154250572-154293534:+,pflag:0
    # Zeamay4EVm000001t2	okay	althi	Zeamay4EVm000001t1	100/87/.	3707,95%,complete	aaref:6846,Sobic.009G173600.1.p,chrmap:96a,100i,11103l,41x,NC_024464.1:154268760-154293534:+,pflag:0
    # Zeamay4EVm000001t4	drop	parthi	Zeamay4EVm000001t1	100/73/.	117,44%,complete-utrpoor	aaref:64.3,Sobic.009G173600.1.p,chrmap:73a,100i,354l,1x,NC_024464.1:154290981-154291238:+,pflag:4,feq:Zeamay4EVm000001t1/altpar73.0.73,Zeamay4EVm000001t2/altpar73.0.73

    my $pflag= ($onotes=~m/pflag:(\d+)/)?$1:0;
    my $nnotes= ($aaref ne "noref")?"aaref:$aaref":"0,0"; #<< funky noref format
    $nnotes  .= ",igloc:$igthis";
    $nnotes  .= ",missIg" if($ig eq $IGMISS or $ig eq "noig");
    if($cla =~ s/,(partof=[\w\.-]+)/,partof/) { my $pof=$1; $nnotes .= ",$pof"; }
    $nnotes  .= ($chrmap ne "nochrmap")?",chrmap:$chrmap":"";
    #? move partof=ID from cla to nnotes?
    $nnotes  .= ",pflag:$pflag"; # last?
    
    $gdnew =~ s/miss(\d+)/m$1/; # shorten
    #xx# $gdnew =~ s/bEVm/cEVm/; # FIXME ID option
    if($OUTFORM =~ /^2|trclass/) {
      ## FIXME: trclass out needs gdnew missNNN info for new locus, in maintd ?
      my $okdrop= ($cla =~ s/^(drop|cull)//)?"drop":"okay";

      #** maintd gdnew BAD for trclass2mainalt.pl, makes blank gdnew.t1 main entry; need to class one of missNNN as new main
      #TESTOFF: $maintd= $gdnew."t1"; # not right? should be same prefix as td, main of, store in newloc{gdnew}{main}=$td ?

      print join("\t",$td,$okdrop,$cla,$maintd,$pid,$aaw,$nnotes)."\n";

    } else {
      $nnotes.=",oldid:$pd"; #?? pd > pubid = $gdnew."t$tinum" instead?
      $gdnew =~ s/bEVm/cEVm/; # FIXME ID option
      my $pubid= $gdnew."t$tinum"; #  want this?  pd == orig id now
      print join("\t",$pubid,$td,$gdnew,$tinum,$cla,$pid,$aaw,$nnotes)."\n";
    }
  } else { # if($OUTFORM == 0) old default
    # $ig  ="$igthis/$ig" if($ig ne $igthis); # all are $IGMISS .. leave off ?
    #OR: 
    $chrmap.=",missIg" if($ig eq $IGMISS or $ig eq "noig"); # $ig=$igthis;

    print join("\t",$td,$pd,$gdnew,$igthis,$cla,$aaref,$chrmap)."\n" ; 
  }
  
} 


sub outgn { 
  my($gd,$td)=@_; 
  my($pd,$cla,$aaw,$pid,$note)= (0,"nocla",0,0,0);
  if( my $gnv=$gnv{$td} ) { ($pd,$cla,$aaw,$pid,$note)= @$gnv; }
  my $ig=$tdig{$td}||"noig"; # or Ig999999 missing 

  my $chrmap=$chrmap{$td}; my $chrmapi=$chrmapi{$td};
  if($NOCHRBE) { map{ s/,N[CW]_.*//; } ($chrmap,$chrmapi); }
  if($chrmap and $chrmapi) { 
    my @cm=split",",$chrmap;
    my @im=split",",$chrmapi;
    $chrmap="$chrmapi/pc/$chrmap";
    ## dont do this way, messy, just  chrmapi/chrmap ?
    # if($im[0] > $cm[0] or $im[2] > $cm[2]) { $chrmap="$chrmapi/pc/$chrmap"; }
    # else { $chrmap="$chrmap/ic/$chrmapi"; }
  } 
  elsif($chrmapi) { $chrmap||=$chrmapi; }

  my $aaref= $aaref{$td}||"noref";
  my $partof=$partof{$td}||""; $cla.=",$partof" if($partof);  

  print join("\t",$td,$pd,$ig,$cla,$aaref,$chrmap)."\n"; 
} 


sub chaindiff {
  my @hd=qw(IchainID Ntd NIg Ngd Ngi Tdloci Gdloci); print join("\t",@hd)."\n"; 
  while(<>){
    chomp; my @v=split"\t"; 
    my($td,$gd,$ig,$cla,$aaref,$loc)=@v;
    next if($gd eq "0" or $ig =~ m/noig|Ig999999|Ig000000/); 
    my($tg,$gg)=map{ my($g,$t)=split"t"; $g; } ($td,$gd); 
    # my($tg,$pg)=map{ my($g,$t)= m/^(\w+)t(\d+)/; $g; } ($td,$pd); # ** utrorf tag bug..

    $igt{$ig}{$tg}++; $igg{$ig}{$gg}++; 
    $tig{$tg}{$ig}++; $gig{$gg}{$ig}++; 
    $tdv{$td}=[@v]; 
  }

  my($nt,$ng,$tg,$gg,$nig,$toi,$goi)=(0)x9; my(%oi,@sn);
  for my $ig (sort keys %igg) { 
    my @tg=sort keys %{$igt{$ig}}; my @gg=sort keys %{$igg{$ig}}; 
    %oi=(); map{ map{$oi{$_}++}keys %{$tig{$_}}} @tg; $toi=scalar(keys %oi);  
    %oi=(); map{ map{$oi{$_}++}keys %{$gig{$_}}} @gg; $goi=scalar(keys %oi);
    $nt=@tg; $ng=@gg; $tg=join",",@tg; $gg=join",",@gg;  
    print join("\t",$ig,$nt,$toi,$ng,$goi,$tg,$gg)."\n"; $nig++; 
    my $i=0; map{ $sn[$i++] += $_; } ($nt,$toi,$ng,$goi); 
    $sn[4]+=($nt==1 and $toi==1)?1:0; $sn[5]+=($ng==1 and $goi==1)?1:0; 
    } 
  my @av=map{ int(100*$_/$nig)/100; } @sn; 
  print "#sum nig=$nig, ave(ntd,nti,ngd,ngi,t11,g11): @av\n"; 
}

__END__

=item newloc tests

try1: 
  #sum.newloc ngid=103918, ntid=536292, nthis=463032, ntadd=131338, ndefer=536292, ndefleft=53214
try2: 
  #sum.newloc ngid=103955, ntid=506923, nthis=439352, ntadd=74359, ndefer=536292???, ndefleft=49222
>> have dup tid, from tidadd
try3:
  #sum.newloc ngid=103955, ntid=506923, nthis=439352, ntadd=47718*, ndefer=536292??, ndefleft=49222

try4: using swapids=1, fixed
  #sum.newloc ngid=118980, ntid=514189, nthis=451389, ntadd=43180, ndefer=62800, ndefleft=41723

try5: using swapids=1, adds Ig9999 list n=208948 (same class as try4)
  #sum.newloc ngid=118980, ntid=514189, nthis=451389, ntadd=43180, ndefer=62800, ndefleft=41723

newloci = 120097, new missN= 5900, orig+add= 114197, vs input ngid=118980
  counting each 'Zeamay4bEVm*missN' as 1 new locus
  .. some are junk loci, tiny detatched exons or mislocated things
  newloci w/ aaref = 31925, 4433 new missN; 88172 newloci w/o aaref
  cut -f3 evg4corn2g.inex9chain.newloctab | sed 's/add//;' | sort -u | wc -l = 120097
  
try5 Ig999 tag, change? all IgNNN/Ig are /Ig999999
Zeamay4bEVm000002t7	Zeamay4EVm000002t13	Zeamay4bEVm000002	Ig000002/Ig999999	altmidfrag	noref	93a,177l,2x,NC_024465.1:70668943-70669063:-
Zeamay4bEVm000007t6	Zeamay4EVm000007t7	Zeamay4bEVm000007	Ig000005/Ig999999	altmidfrag	noref	100a,255l,1x,NC_024465.1:170154965-170155217:+



try4 relocs
eg: Zeamay4bEVm008969t, 72 tr, 5 Ig locs, 34 tr w/o reloc, 38 tr with reloc.
   1 Ig010131  23 Ig010380  1 Ig010489 13 Ig010577  34 Ig019320

try4:misses, ndefleft=41723
many miss1=33545, with miss2=5388, miss3=1587,... miss6=105, miss9=33, miss10=30, miss20=7 .. ok? bug?

....

Zeamay4bEVm012031t11	Zeamay4EVm012536t23	Zeamay4bEVm012031miss6	Ig020559	althi1,relocmiss6	noref	100a,291l,3x,NC_024464.1:30860623-125259335:+
Zeamay4bEVm006494t15	Zeamay4EVm001308t76	Zeamay4bEVm006494miss6	Ig020330	althi,relocmiss6	noref	97a,516l,3x,NC_024468.1:3486555-128952756:-
Zeamay4bEVm002995t21	Zeamay4EVm003108t21	Zeamay4bEVm002995miss6	Ig005984	althi,relocmiss6	noref	90a,1560l,3x,NC_024466.1:102231961-102233594:+
Zeamay4bEVm008053t21	Zeamay4EVm120908t1	Zeamay4bEVm008053miss6	noig	dropalthi1,relocmiss6	noref	nochrmap
..
Zeamay4bEVm000970t147	Zeamay4EVm000998t155	Zeamay4bEVm000970miss20	Ig020155	althi1,relocmiss20	noref	100a,699l,2x,NC_024461.1:10188759-50005804:-
Zeamay4bEVm000970t165	Zeamay4EVm000998t222	Zeamay4bEVm000970miss20	Ig020155	altmidfrag,relocmiss20	noref	100a,726l,2x,NC_024461.1:10188640-50005797:-
Zeamay4bEVm000970t169	Zeamay4EVm000998t185	Zeamay4bEVm000970miss20	Ig020155	altmidfrag,relocmiss20	noref	100a,705l,2x,NC_024461.1:10188757-50005812:-
Zeamay4bEVm001421t229	Zeamay4EVm001477t190	Zeamay4bEVm001421miss20	Ig018552	altmidfrag,relocmiss20	124,arath15:ATMG00810.1,	91a,669l,3x,NC_024459.1:81406-300721744:+
Zeamay4bEVm000319t367	Zeamay4EVm000328t334	Zeamay4bEVm000319miss20	Ig020099	altmidfrag,relocmiss20	noref	72a,516l,3x,NC_024464.1:116844507-161757580:-
Zeamay4bEVm000798t448	Zeamay4EVm000826t586	Zeamay4bEVm000798miss20	Ig020138	altmidfrag,relocmiss20	noref	91a,366l,3x,NC_024467.1:27944677-45434268:-


try4: same mixed Ig, one pubid loc
Zeamay4bEVm004063t4	Zeamay4EVm004219t8	Zeamay4bEVm004063	Ig005332	althi	951,Sobic.001G147200.1.p,	100a,1848l,6x,NC_024463.1:15221337-15225508:+
Zeamay4bEVm004063t5	Zeamay4EVm139673t4	Zeamay4bEVm004063	Ig005332	althi1	960,Sobic.001G147200.1.p,	100a,1848l,6x,NC_024463.1:15221337-15225508:+
Zeamay4bEVm004063t6	Zeamay4EVm139673t2	Zeamay4bEVm004063	Ig005332	althi1	946,Sobic.001G147200.1.p,	96a,1929l,7x,NC_024463.1:15221337-15225508:+
Zeamay4bEVm004063t7	Zeamay4EVm139673t3	Zeamay4bEVm004063	Ig005332	althi1	960,Sobic.001G147200.1.p,refok,	100a,1848l,6x,NC_024463.1:15221337-15225508:+
Zeamay4bEVm004063t8	Zeamay4EVm004219t5	Zeamay4bEVm004063	Ig005332	althi1	961,Sobic.001G147200.1.p,	100a,1848l,6x,NC_024463.1:15221337-15225508:+
Zeamay4bEVm004063t9	Zeamay4EVm004219t6	Zeamay4bEVm004063	Ig005332	althi	947,Sobic.001G147200.1.p,	100a,1848l,6x,NC_024463.1:15221337-15225508:+
Zeamay4bEVm004063t10	Zeamay4EVm004219t7	Zeamay4bEVm004063	Ig005332	althi1	961,Sobic.001G147200.1.p,	100a,1848l,6x,NC_024463.1:15221337-15225508:+
Zeamay4bEVm004063t11	Zeamay4EVm139673t5	Zeamay4bEVm004063	Ig005332	althi1	904,Sobic.001G147200.1.p,refok,	100a,1772l,6x,NC_024463.1:15221418-15225508:+
Zeamay4bEVm004063t12	Zeamay4EVm139673t6	Zeamay4bEVm004063	Ig005332	althi	919,Sobic.001G147200.1.p,	100a,1890l,7x,NC_024463.1:15221337-15225439:+
Zeamay4bEVm004063t14	Zeamay4EVm004219t9	Zeamay4bEVm004063	Ig005332	althi	748,Sobic.001G147200.1.p,	100a,1518l,4x,NC_024463.1:15223568-15225508:+
Zeamay4bEVm004063t15	Zeamay4EVm004219t10	Zeamay4bEVm004063	Ig005332	althi	618,Sobic.001G147200.1.p,	100a,1356l,4x,NC_024463.1:15223916-15225508:+
Zeamay4bEVm004063t16	Zeamay4EVm139673t8	Zeamay4bEVm004063	Ig005332	althi	733,Sobic.001G147200.1.p,	100a,1479l,3x,NC_024463.1:15223703-15225508:+
Zeamay4bEVm004063t17	Zeamay4EVm139673t9	Zeamay4bEVm004063	Ig005332	althi1	688,Sobic.001G147200.1.p,	96a,1494l,4x,NC_024463.1:15223758-15225508:+
Zeamay4bEVm004063t19	Zeamay4EVm139673t10	Zeamay4bEVm004063	Ig005332	althi	588,Sobic.001G147200.1.p,	99a,1275l,6x,NC_024463.1:15221436-15224759:+
Zeamay4bEVm004063t20	Zeamay4EVm139673t11	Zeamay4bEVm004063	Ig005332	althi	516,Sobic.001G147200.1.p,	100a,1029l,4x,NC_024463.1:15221337-15224362:+

Zeamay4bEVm004063t1	Zeamay4EVm004219t1	Zeamay4bEVm004063miss1	Ig004030	main,relocmiss1	1075,Sobic.001G147200.1.p,refgood,	100a,2028l,6x,NC_024459.1:258931311-258935850:+
Zeamay4bEVm004063t2	Zeamay4EVm004219t2	Zeamay4bEVm004063miss1	Ig004030	althi,relocmiss1	1098,Sobic.001G147200.1.p,refbest,	100a,2028l,6x,NC_024459.1:258931311-258935850:+
Zeamay4bEVm004063t3	Zeamay4EVm139673t1	Zeamay4bEVm004063miss1	Ig004030	althi,relocmiss1	1085,Sobic.001G147200.1.p,refgood,	100a,2028l,6x,NC_024459.1:258931311-258935850:+
Zeamay4bEVm004063t13	Zeamay4EVm139673t7	Zeamay4bEVm004063miss1	Ig004030	althi,relocmiss1	737,Sobic.001G147200.1.p,	96a,1584l,5x,NC_024459.1:258933882-258935850:+
Zeamay4bEVm004063t18	Zeamay4EVm004219t11	Zeamay4bEVm004063miss1	Ig004030	althi1,relocmiss1	617,Sobic.001G147200.1.p,	97a,1266l,4x,NC_024459.1:258934107-258935850:+
Zeamay4bEVm004063t21	Zeamay4EVm004219t12	Zeamay4bEVm004063miss1	Ig004030	althi1,relocmiss1	443,Sobic.001G147200.1.p,	100a,786l,4x,NC_024459.1:258931311-258934082:+
Zeamay4bEVm004063t22	Zeamay4EVm004219t13	Zeamay4bEVm004063miss1	Ig004030	althi1,relocmiss1	443,Sobic.001G147200.1.p,	100a,786l,4x,NC_024459.1:258931311-258934082:+

try3:
 egrep  '^(Zeamay4EVm004219t|Zeamay4EVm139673t)' evg4corn2g.inex9chain.newloctab
 * note both belong to Zeamay4bEVm004063t, but 2 Ig locs
Zeamay4EVm004219t5	Zeamay4bEVm004063t8	Zeamay4EVm004219	Ig005332	althi1	961,Sobic.001G147200.1.p,	100a,1848l,6x,NC_024463.1:15221337-15225508:+
Zeamay4EVm004219t6	Zeamay4bEVm004063t9	Zeamay4EVm004219	Ig005332	althi	947,Sobic.001G147200.1.p,	100a,1848l,6x,NC_024463.1:15221337-15225508:+
Zeamay4EVm004219t7	Zeamay4bEVm004063t10	Zeamay4EVm004219	Ig005332	althi1	961,Sobic.001G147200.1.p,	100a,1848l,6x,NC_024463.1:15221337-15225508:+
Zeamay4EVm004219t8	Zeamay4bEVm004063t4	Zeamay4EVm004219	Ig005332	althi	951,Sobic.001G147200.1.p,	100a,1848l,6x,NC_024463.1:15221337-15225508:+
Zeamay4EVm004219t9	Zeamay4bEVm004063t14	Zeamay4EVm004219	Ig005332	althi	748,Sobic.001G147200.1.p,	100a,1518l,4x,NC_024463.1:15223568-15225508:+
Zeamay4EVm004219t10	Zeamay4bEVm004063t15	Zeamay4EVm004219	Ig005332	althi	618,Sobic.001G147200.1.p,	100a,1356l,4x,NC_024463.1:15223916-15225508:+
Zeamay4EVm139673t2	Zeamay4bEVm004063t6	Zeamay4EVm004219add	Ig005332	althi1,relocadd	946,Sobic.001G147200.1.p,	96a,1929l,7x,NC_024463.1:15221337-15225508:+
Zeamay4EVm139673t3	Zeamay4bEVm004063t7	Zeamay4EVm004219add	Ig005332	althi1,relocadd	960,Sobic.001G147200.1.p,refok,	100a,1848l,6x,NC_024463.1:15221337-15225508:+
Zeamay4EVm139673t4	Zeamay4bEVm004063t5	Zeamay4EVm004219add	Ig005332	althi1,relocadd	960,Sobic.001G147200.1.p,	100a,1848l,6x,NC_024463.1:15221337-15225508:+
Zeamay4EVm139673t5	Zeamay4bEVm004063t11	Zeamay4EVm004219add	Ig005332	althi1,relocadd	904,Sobic.001G147200.1.p,refok,	100a,1772l,6x,NC_024463.1:15221418-15225508:+
Zeamay4EVm139673t6	Zeamay4bEVm004063t12	Zeamay4EVm004219add	Ig005332	althi,relocadd	919,Sobic.001G147200.1.p,	100a,1890l,7x,NC_024463.1:15221337-15225439:+
Zeamay4EVm139673t8	Zeamay4bEVm004063t16	Zeamay4EVm004219add	Ig005332	althi,relocadd	733,Sobic.001G147200.1.p,	100a,1479l,3x,NC_024463.1:15223703-15225508:+
Zeamay4EVm139673t9	Zeamay4bEVm004063t17	Zeamay4EVm004219add	Ig005332	althi1,relocadd	688,Sobic.001G147200.1.p,	96a,1494l,4x,NC_024463.1:15223758-15225508:+
Zeamay4EVm139673t10	Zeamay4bEVm004063t19	Zeamay4EVm004219add	Ig005332	althi,relocadd	588,Sobic.001G147200.1.p,	99a,1275l,6x,NC_024463.1:15221436-15224759:+
Zeamay4EVm139673t11	Zeamay4bEVm004063t20	Zeamay4EVm004219add	Ig005332	althi,relocadd	516,Sobic.001G147200.1.p,	100a,1029l,4x,NC_024463.1:15221337-15224362:+

Zeamay4EVm004219t1	Zeamay4bEVm004063t1	Zeamay4EVm004219miss	Ig004030	main,relocmiss	1075,Sobic.001G147200.1.p,refgood,	100a,2028l,6x,NC_024459.1:258931311-258935850:+
Zeamay4EVm139673t1	Zeamay4bEVm004063t3	Zeamay4EVm139673miss	Ig004030	althi,relocmiss	1085,Sobic.001G147200.1.p,refgood,	100a,2028l,6x,NC_024459.1:258931311-258935850:+
Zeamay4EVm004219t2	Zeamay4bEVm004063t2	Zeamay4EVm004219miss	Ig004030	althi,relocmiss	1098,Sobic.001G147200.1.p,refbest,	100a,2028l,6x,NC_024459.1:258931311-258935850:+
Zeamay4EVm139673t7	Zeamay4bEVm004063t13	Zeamay4EVm139673miss	Ig004030	althi,relocmiss	737,Sobic.001G147200.1.p,	96a,1584l,5x,NC_024459.1:258933882-258935850:+
Zeamay4EVm004219t11	Zeamay4bEVm004063t18	Zeamay4EVm004219miss	Ig004030	althi1,relocmiss	617,Sobic.001G147200.1.p,	97a,1266l,4x,NC_024459.1:258934107-258935850:+
Zeamay4EVm004219t12	Zeamay4bEVm004063t21	Zeamay4EVm004219miss	Ig004030	althi1,relocmiss	443,Sobic.001G147200.1.p,	100a,786l,4x,NC_024459.1:258931311-258934082:+
Zeamay4EVm004219t13	Zeamay4bEVm004063t22	Zeamay4EVm004219miss	Ig004030	althi1,relocmiss	443,Sobic.001G147200.1.p,	100a,786l,4x,NC_024459.1:258931311-258934082:+

try2:
grep  '^Zeamay4EVm139673|Ig005332' evg4corn2g.inex9chain.newloctab

Zeamay4EVm004219t5	Zeamay4bEVm004063t8	Zeamay4EVm004219	Ig005332	althi1	961,Sobic.001G147200.1.p,	100a,1848l,6x,NC_024463.1:15221337-15225508:+
Zeamay4EVm004219t6	Zeamay4bEVm004063t9	Zeamay4EVm004219	Ig005332	althi	947,Sobic.001G147200.1.p,	100a,1848l,6x,NC_024463.1:15221337-15225508:+
Zeamay4EVm004219t7	Zeamay4bEVm004063t10	Zeamay4EVm004219	Ig005332	althi1	961,Sobic.001G147200.1.p,	100a,1848l,6x,NC_024463.1:15221337-15225508:+
Zeamay4EVm004219t8	Zeamay4bEVm004063t4	Zeamay4EVm004219	Ig005332	althi	951,Sobic.001G147200.1.p,	100a,1848l,6x,NC_024463.1:15221337-15225508:+
Zeamay4EVm004219t9	Zeamay4bEVm004063t14	Zeamay4EVm004219	Ig005332	althi	748,Sobic.001G147200.1.p,	100a,1518l,4x,NC_024463.1:15223568-15225508:+
Zeamay4EVm004219t10	Zeamay4bEVm004063t15	Zeamay4EVm004219	Ig005332	althi	618,Sobic.001G147200.1.p,	100a,1356l,4x,NC_024463.1:15223916-15225508:+

>> bad sort?
Zeamay4EVm139673t10	Zeamay4bEVm004063t19	Zeamay4EVm004219add	Ig005332	althi,relocadd	588,Sobic.001G147200.1.p,	99a,1275l,6x,NC_024463.1:15221436-15224759:+
Zeamay4EVm139673t11	Zeamay4bEVm004063t20	Zeamay4EVm004219add	Ig005332	althi,relocadd	516,Sobic.001G147200.1.p,	100a,1029l,4x,NC_024463.1:15221337-15224362:+
Zeamay4EVm139673t2	Zeamay4bEVm004063t6	Zeamay4EVm004219add	Ig005332	althi1,relocadd	946,Sobic.001G147200.1.p,	96a,1929l,7x,NC_024463.1:15221337-15225508:+
Zeamay4EVm139673t3	Zeamay4bEVm004063t7	Zeamay4EVm004219add	Ig005332	althi1,relocadd	960,Sobic.001G147200.1.p,refok,	100a,1848l,6x,NC_024463.1:15221337-15225508:+
Zeamay4EVm139673t4	Zeamay4bEVm004063t5	Zeamay4EVm004219add	Ig005332	althi1,relocadd	960,Sobic.001G147200.1.p,	100a,1848l,6x,NC_024463.1:15221337-15225508:+
Zeamay4EVm139673t5	Zeamay4bEVm004063t11	Zeamay4EVm004219add	Ig005332	althi1,relocadd	904,Sobic.001G147200.1.p,refok,	100a,1772l,6x,NC_024463.1:15221418-15225508:+
Zeamay4EVm139673t6	Zeamay4bEVm004063t12	Zeamay4EVm004219add	Ig005332	althi,relocadd	919,Sobic.001G147200.1.p,	100a,1890l,7x,NC_024463.1:15221337-15225439:+
Zeamay4EVm139673t8	Zeamay4bEVm004063t16	Zeamay4EVm004219add	Ig005332	althi,relocadd	733,Sobic.001G147200.1.p,	100a,1479l,3x,NC_024463.1:15223703-15225508:+
Zeamay4EVm139673t9	Zeamay4bEVm004063t17	Zeamay4EVm004219add	Ig005332	althi1,relocadd	688,Sobic.001G147200.1.p,	96a,1494l,4x,NC_024463.1:15223758-15225508:+
>> misses, probably 2nd Ig for 1 locus only Zeamay4EVm139673
Zeamay4EVm139673t1	Zeamay4bEVm004063t3	Zeamay4EVm139673miss	Ig004030	althi,relocmiss	1085,Sobic.001G147200.1.p,refgood,	100a,2028l,6x,NC_024459.1:258931311-258935850:+
Zeamay4EVm139673t7	Zeamay4bEVm004063t13	Zeamay4EVm139673miss	Ig004030	althi,relocmiss	737,Sobic.001G147200.1.p,	96a,1584l,5x,NC_024459.1:258933882-258935850:+
  
=cut

=item chaindiff.tab

perl -ne \
'($td,$gd,$ig,$cla,$aaref,$loc)=@v=split; next if($gd eq "0" or $ig =~
m/noig|Ig999999|Ig000000/); ($tg,$gg)=map{ ($g,$t)=split"t"; $g; }
($td,$gd); $igt{$ig}{$tg}++; $igg{$ig}{$gg}++; $tig{$tg}{$ig}++;
$gig{$gg}{$ig}++; $tdv{$td}=[@v]; 
END{ for $ig (sort keys %igg) { 
@tg=sort keys %{$igt{$ig}}; @gg=sort keys %{$igg{$ig}}; %oi=(); map{
map{$oi{$_}++}keys %{$tig{$_}}} @tg; $toi=scalar(keys %oi);  %oi=();
map{ map{$oi{$_}++}keys %{$gig{$_}}} @gg; $goi=scalar(keys %oi);
$nt=@tg; $ng=@gg; $tg=join",",@tg; $gg=join",",@gg;  print
join("\t",$ig,$nt,$toi,$ng,$goi,$tg,$gg)."\n"; $nig++; $i=0; map{
$sn[$i++] += $_; } ($nt,$toi,$ng,$goi); $sn[4]+=($nt==1 and
$toi==1)?1:0; $sn[5]+=($ng==1 and $goi==1)?1:0; } @av=map{
int(100*$_/$nig)/100; } @sn; print "#sum nig=$nig,
ave(ntd,nti,ngd,ngi,t11,g11): @av\n"; } 
BEGIN{ @hd=qw(IchainID Ntd NIg Ngd Ngi Tdloci Gdloci); print join("\t",@hd)."\n"; } ' \
 evg4corn2g.pubclass_inex1chaincds.idxtab > evg4corn2g.inex1chaindiff.tab

=cut

=item input data

  ?  bug now for ichain tabulator, mrna blastn input but -cdstrim flag, measures align outside of cds, using cdssize(tw)
  fix tabulator.pl, but patch here by using col9 = mrna span (1-3488/ below)


# ... data1
# evg4corn/geno1f/evg*.inv4chainct.tab
# fix for aln/tw above for this -cdstrim but mrna extended align
# now aln,tw = 3379,2231; update aln,twx = 3379,3488
# Query   Source  Bits    Ident   Align   Qlen    Slen    Nexon   ILocus  Qexons  Sintrons        Geneinfo
# Zeamay4EVm000071t89     NC_024467.1     6244    3379    3379    2231    0       8       Ig000742        1-3488/1-1157,1158-1327,1328-1748,1749-1962,1963-2038,2091-2205,2206-2278,2352-3488     115910756-115915224:-/N1032263-N1032269,N1032271-N1032277,N1032285-N1032289,N1032290-N1032293,N65661-N65662,N65657-N65660,N65653-N65655,N65630-N65652   split=2,NC_024459.1,NC_024467.1;
# evg4corn2g/pubnew/evg4corn2g_cds-genoasm.inv7chainct.tab
# ncols=11 or 12
# Query	Source	Bits	Ident	Align	Qlen	Slen	Nexon	ILocus	Qexons	Sintrons	Geneinfo
# Zeamay4bEVm000001t3	NC_024464.1	15886	8625	8643	8922	0	27	Ig000002	1-8922/1-85,123-240,241-391,493-622,623-697,698-823,824-924,925-1018,1019-1154,1155-1209,1210-1386,1387-1510,1511-1669,1716-3376,3377-3443,3444-3640,3641-3744,3745-4049,4050-4564,4565-6369,6370-7908,7909-8095,8096-8163,8236-8481,8482-8617,8618-8714,8781-8922	154252404-154274323:+/N618115-N618116,N618117-N618118,N618119-N618120,N618121-N618122,N618123-N618124,N618125-N618126,N618127-N618128,N618129-N618130,N618131-N618132,N618133-N618134,N618135-N618136,N618137-N618138,N618139-N618141,N618140-N618142,N618143-N618144,N618145-N618146,N618147-N618148,N618149-N618150,N618151-N618153,N618152-N618157,N618158-N618159,N618160-N618161,N618162-N618163,N618164-N618165,N618166-N618167,N618168-N618169,N618170-N618171

#........ data2
# ncols=8
# evg4corn/geno1f/tgclass.pubids
# Zeamay4gEVm000002t5     Zeamay4EVm000002t5      Zeamay4gEVm000002       5       althi   2619,96%,partial3       100/73/.        aaref:4841,Sobic.007G098500.1.p,chrmap:42a,100i,7821l,6x,NC_024465.1:70643354-70658828:+,pflag:0,feq:Zeamay4EVm000002t4/altpar11.0.0
# evg4corn2g/publicset/evg4corn2g.pubids
#Public_mRNA_ID	originalID	PublicGeneID	AltNum	Class	AAqual	pIdAln	Notes
# Zeamay4bEVm000001t1	Zeamay4EVm000001t1	Zeamay4bEVm000001	1	main	5425,97%,complete	100/87/.	aaref:10109,Sobic.009G173600.1.p,refbest,

=cut
