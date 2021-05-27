#!/usr/bin/env perl
# overeqcdsloc.pl evg2anofunz4b_cds.tall3d >  evg2anofunz4b_cds.eqgenebl

=item about

  evigene/scripts/prot/overeqcdsloc.pl
  equalgene.pl/overgenedup.pl equivalent for blastn/tall table inputs
  input : blastn/tall table w/ cds-exon spans
  output: table of (id, oid, overids, location) as per equalgene.pl

  env noalts=1 strand=1 ./overeqcdsloc.pl \
    < genoAfunF1_evg2anofunz4b_cds.tall3d > genoAfunF1_evg2anofunz4b_cds.eqgenebl
  
=cut


use strict;
use Getopt::Long;

my $BINSIZE= 1000;
my $MINOVERLAP= $ENV{minover}||10; # percent, maybe too high? 5 or less?
my $stranded= $ENV{strand}||0; #? default on
my $ALTSHOW=  $ENV{alts}||0; #? was noalts default on
# my $UNALTS=  $ENV{unalt}||0; # reverse over() test for alts only, ie. which may be paralogs
my $MAXLIST=  $ENV{maxlist}||0; #was 10; # list of overlap ids .. show all unless ask for less
my $SHORTALT= $ENV{shortalt}||75; # percent of longest
my $DEBUG= $ENV{debug}|| 0;

my $debugOVERGENE=$ENV{overgene}|| 0; # BUGGY dont use yet, debug or drop
my $SHOWEQXON= $DEBUG; 

my ($bltab, $hotab, $outtab);

my $optok= GetOptions( 
  "bltab|input=s", \$bltab, ## << primary input
  "eqgene|output=s", \$outtab, ## << primary output
  "hotab=s", \$hotab, ## aaref blast score tab
  "MINOVERLAP=i", \$MINOVERLAP,  
  "stranded!", \$stranded, 
  "ALTSHOW!", \$ALTSHOW, #was rev $NOALTS, 
  "EQXONSHOW!", \$SHOWEQXON, 
  "MAXLIST=i", \$MAXLIST,  # may need all ?
  "overgene!", \$debugOVERGENE, 
  "DEBUG!", \$DEBUG, ## was $debug, ## $DEBUG now, dont need both
);

#? $bltab= shift @ARGV unless($bltab); #? only input file, can be stdin|-

die "usage:  overeqcdsloc.pl -in evgcds2genome.blastscoretab -out evgcds2genome.eqgene
  opts: -[no]stranded -[no]altshow  -minover $MINOVERLAP -debug
  gene locus equivalences, from blastscoretab of evigene/makeblastscore3.pl"
  unless($optok and $bltab);
  

my(%afloc, %afid, %afidloc, %afscore, %alts, %afexons, %ovw, %eqxon); # ovw == overlapwidth{aid}{bid}

my $OUTH= *STDOUT; 
if($outtab and $outtab !~ /stdout|^-/) {
  rename($outtab,"$outtab.old") if( -f $outtab );
  open(OUTC,">$outtab") or die "write $outtab";
  $OUTH= *OUTC;
}

my($nho,$hoscore)= readAablastTab($hotab) if($hotab);
readTallAlnTab($bltab);   
putEqualGeneTab($outtab); # to $OUTH

#------------ Blast2EqualGeneTab subs --------------

=item bugs? recheck/revise

>> these identical spans were missed in this output .. why?
egrep '^Anofunz4hEVm00044[5678]t1' evgm3chranofun-evg2anofunz4h.btall | cut -f1,2,5-
Anofunz4hEVm000445t1	KB668855	5460	5655	0	1-5609/1-375,376-2085,2086-2962,2963-4559,4560-4656,4657-5144,5309-5474,5475-5609	KB668855:40958-41332,41413-43122,43186-44066,44147-45743,45803-45905,45997-46484,46806-46971,47042-47181:+	0d	0d
Anofunz4hEVm000446t1	KB668855	5460	5655	0	1-5609/1-375,376-2085,2086-2962,2963-4559,4560-4656,4657-5144,5309-5474,5475-5609	KB668855:40958-41332,41413-43122,43186-44066,44147-45743,45803-45905,45997-46484,46806-46971,47042-47181:+	0d	0d
Anofunz4hEVm000447t1	KB668855	5460	5658	0	1-5609/1-375,376-2085,2086-2962,2963-4559,4560-4656,4657-5144,5309-5474,5475-5609	KB668855:40958-41332,41413-43122,43186-44066,44147-45743,45803-45905,45997-46484,46806-46971,47042-47181:+	0d	0d
Anofunz4hEVm000448t1	KB668855	5460	5658	0	1-5609/1-375,376-2085,2086-2962,2963-4559,4560-4656,4657-5144,5309-5474,5475-5609	KB668855:40958-41332,41413-43122,43186-44066,44147-45743,45803-45905,45997-46484,46806-46971,47042-47181:+	0d	0d

[ux455375@comet-ln3 geno1f]$ egrep '^Anofunz4hEVm000448t1' evgm3chranofun-evg2anofunz4h.eqgene     
Anofunz4hEVm000448t1	noid	Anofunz4hEVm000447t7/16,Anofunz4hEVm000447t6/15,Anofunz4hEVm000447t5/13,Anofunz4hEVm000447t3/13,Anofunz4hEVm000445t4/12,Anofunz4hEVm000447t2/10,Anofunz4hEVm000447t4/10	KB668855:40958-47181:+	97a,99i,5658l,8x	0
Anofunz4hEVm000445t1	noid	Anofunz4hEVm000447t7/16,Anofunz4hEVm000447t6/15,Anofunz4hEVm000448t3/15,Anofunz4hEVm000447t5/13,Anofunz4hEVm000447t3/13,Anofunz4hEVm000447t2/10,Anofunz4hEVm000447t4/10	KB668855:40958-47181:+	97a,99i,5655l,8x	0
Anofunz4hEVm000446t1	noid	Anofunz4hEVm000447t7/16,Anofunz4hEVm000447t6/15,Anofunz4hEVm000448t3/15,Anofunz4hEVm000447t5/13,Anofunz4hEVm000447t3/13,Anofunz4hEVm000445t4/12,Anofunz4hEVm000447t2/10,Anofunz4hEVm000447t4/10	KB668855:40958-47181:+	97a,99i,5655l,8x	0
Anofunz4hEVm000447t1	noid	Anofunz4hEVm000448t3/15,Anofunz4hEVm000445t4/12	KB668855:40958-47181:+	97a,99i,5658l,8x	0

=cut

sub readTallAlnTab {
  my($intab)=@_;
  my($ok,$inh) = openRead($intab); 
  my($nid, $nxt)=(0) x 2;
  
  while(<$inh>) {
    next if(/^Query|^\W/);
    chomp; my @v=split"\t"; 
    
    # format of: $evigene/scripts/makeblastscore3.pl -tall -spans=2 -oneref=1  trasm-genome.blastn > table
    my($id, $grm, $bits,$idn,$aln,$tw,$gw, $xspans, $gspans, $tdups, $gdups)= @v; 
  
    my($trspan)= ($xspans=~s,^([\d\-]+)/,,)?$1:0;
    my($goo)= ($gspans=~s/:([\+\-])$//)?$1:"+";
    # my @tx=  split",",$xspans; # should be same n @tx == @gx
    my @gx=  split",",$gspans; # includes gr's chrs ids with splits
    my $ngx=$#gx;
    my $lgr=$grm; my @grs=();
    my($gstart,$gend,$gsplit,$galn)=(0,0,0,0); 
    my $gor= "."; my $gorm=$goo; 

    (my $gd=$id) =~ s/t\d+$//;
    $alts{$gd}{$id}++;
    
    ## efficiency change: store only gene ids in afloc{gr}{bins},
    ## and separate geneexons{id}=@xbe hash, for revised over() method
    # %afidloc{gr}{bn}[ids], %afexons{id}=[@gx] replace afloc
    $afexons{$id}= \@gx;
    
    for my $i (0..$ngx) {
      my $gx= $gx[$i];  # my $tx=$tx[$i];
      my($gr)= ($gx=~s/^([\w\.-]+)://)? $1: $lgr; 
      if($gr ne $lgr) { push @grs, $gr; $gsplit++; $gor="."; }
      ## Splits: strand $goo is wrong .. get from @gx order
      my($xb,$xe)= $gx=~m/(\d+).(\d+)/;
      if($gor eq ".") { 
        my($xbn)= $gx[$i+1]=~m/(\d+).(\d+)/;
        $gor= ($i<$ngx) ? (($xb > $xbn)?"-":"+") : $goo; # Splits strand bug
        $gorm=$gor if($gr eq $grm); 
        } 
      $lgr=$gr;
      if($xe) { 
        unless($gsplit and $gr ne $grm) { 
          $gstart=$xb if($gstart==0 or $xb<$gstart); 
          $gend= $xe if($xe > $gend); 
          $galn += 1+ $xe-$xb;
        }
        my $afloc=join"\t",$id,$gr,$gor,$xb,$xe; $nxt++;
        for my $bn (glocbins($xb,$xe)) { 
          push @{$afloc{$gr}{$bn}},$afloc; 
          $afidloc{$gr}{$bn}{$id}++; # replacement
          }
      }
    }
    
    ## add $goo, gspans == exon spans w/o other??
    my $gloc="$grm:$gstart-$gend:$gorm";
    if(@grs) { 
      my $spaln= $aln - $galn;
      if($spaln<20) { $gsplit=0; } 
      else { my $p= pctof($spaln,$tw); @grs=grep{$_ ne $grm}@grs; $gsplit = "$p%,".join",",@grs; } # int(0.5+ 100*$spaln/$tw); 
    }
    $afid{$id}= join"\t", $bits, $idn, $aln, $tw, $xspans, $grm, $goo, $gloc, $gsplit, $gspans; ## "$gr:$xl\t$cw\t$nam";  
    $nid++;
    
    my $hoval= $hoscore->{$id}||0;
    my $pidn= int(100*$idn/$tw);
    my $score= $hoval * 9 + $pidn * 5 + $ngx * 2 + $idn * 0.01;
    $score /= 2 if ($gsplit);
    $afscore{$id}= $score;
  }
  
  warn "# readTallAlnTab: nid=$nid, nexon=$nxt\n" if $DEBUG;
  return($nid);
}

sub readAAnametab
{
  my($aanametab)= @_;
  my($nid, $nxt)=(0) x 9; my %hoscore;
  my($ok,$inh) = openRead($aanametab); 
  while(<$inh>) {  
    next if(/^Query|^\W/);
    my($td,$name,$alnscore,$rd,@more)=split"\t"; 
    # alnscore format: 72%,3270/4555,3282 ; may be only pctaln score (\d+)
    #x my($apct,$aln,$refsize,$tsize)= $alnscore =~ m=(^\d+)%,(\d+)/(\d+),(\d+)=;
    my($apct,$aln)= $alnscore =~ m=^(\d+)%,(\d+)\b=;
    unless($aln){ ($aln)= $alnscore =~ m/^(\d+)/; }
    if($aln > 0 and $aln > $hoscore{$td}) {
      $hoscore{$td}=$aln; $nid++; }
  } close($inh);
  warn "# readAAnametab: nid=$nid \n" if $DEBUG;
  return($nid,\%hoscore);
}

=item readAablastTab or readAAnametab

  readAAnametab() alt aablast data format
Anofunz4hEVm000001t1  Name_of_protein   95%,11401/11410,15789    AGAP009554-PA
Anofunz4hEVm000001t1  Name_of_protein   93%,10657/11328,15789    AAEL007898-PA
  readAablastTab() format
Anofunz4hEVm000001t1	AGAP009554-PA	19843	10113	11401	15789	11410
Anofunz4hEVm000001t1	AAEL007898-PA	11434	6103	10657	15789	11328

=cut

sub readAablastTab { 
  my($intab)=@_;
  if($intab =~ /\.names$/) { return readAAnametab($intab); }
  my($ok,$inh) = openRead($intab); 
  my($nid, $nxt)=(0) x 9; my %hoscore;
  while(<$inh>) {
    next if(/^Query|^\W/);
    chomp; my @v=split"\t"; 
    my($td,$rd,$bs,$idn,$aln,$tw,$rw)=@v;
    $hoscore{$td}=$idn if($idn>$hoscore{$td});
    $nid++;
  } close($inh);
  warn "# readAablastTab: nid=$nid \n" if $DEBUG;
  return($nid,\%hoscore);
}


sub putEqualGeneTab {
  my($outtab)=@_;
  my($tnov, $tnmisalt)=(0,0);
  
  # my @afids= sort _genealts keys %afid; # bad?? no, just extra v6 gene set
  my @afids= sort {   
    my($ag,$at)= $a=~m/(\w+)t(\d+)/;  
    my($bg,$bt)= $b=~m/(\w+)t(\d+)/; 
    $ag=~s/Anofun6EV/anofun6EV/; $bg=~s/Anofun6EV/anofun6EV/; #hack data test
    $afscore{$b} <=> $afscore{$a}  or $ag cmp $bg or $at <=> $bt or $a cmp $b;
    } keys %afid; 

  for my $id (@afids) { 
    my($bits, $idn, $aln, $tw, $xspans, $gr, $goo, $gloc, $gsplit, $gspans)= split"\t",$afid{$id};
    
      # over() step appears to be slow .. is it efficient?
    my ($nov,$ovids);
    if($debugOVERGENE) { 
      ($nov,$ovids)= overgene($id,$gr,$goo,$gspans);
    } else {
      ($nov,$ovids)= over($id,$gr,$goo,$gspans);
    }
    $tnov+=$nov;
    
    #*? add col for alt loc <> main/other alt ?? ie NOALTS but invert over() test for alts
    # .. do per id, test gn000t2,3,4 vs gn000t1 only? or gn000 all x all? or gn000ti x majority of others?
    # .. alt_notover : add positive count of alts over, for NOALTS/ALTSHOW
    # FIXME: filter out uninteresting misaltids, ie tiny alts no overlap

    my($nmisalt,$misaltids,$naltover)= ($ALTSHOW) ? (0,0,0) : alt_notover($id,$ovids); 
    $tnmisalt+= $nmisalt;
    
    my $nxon= 1 + $xspans=~tr/,/,/;
    my $mapval= mapval($idn,$aln,$tw,$nxon,$gsplit); # add nexon, later valid splices

    putover( $id, "noid",$tw, $ovids, $gloc, $mapval, $misaltids); # add out cols? aln/idn/tw ?
  }
  warn "# putEqualGeneTab: nover=$tnov, nmisalt=$tnmisalt\n" if $DEBUG;
}

#---------

sub _min{ ($_[0]<$_[1])?$_[0]:$_[1]; } 
sub _max{ ($_[0]>$_[1])?$_[0]:$_[1]; }

sub glocbins{ my($tb,$te)=@_; return() if($te<1); return(int($tb/$BINSIZE) .. int($te/$BINSIZE)); }

sub _genealts {   # sort by geneid, alt num, for Evigene IDs  gene0000t123
  my($ag,$at)= $a=~m/(\w+)t(\d+)/; 
  my($bg,$bt)= $b=~m/(\w+)t(\d+)/; 
  return ($ag cmp $bg or $at <=> $bt or $a cmp $b);
}

sub pctof {
  my($nn,$nd)=@_;
  my $pn= ($nd>0)? int(0.5+ 100*$nn/$nd): 0;  $pn=100 if($pn>100);
  return $pn; 
}

sub mapval { 
  my($idn,$aln,$tw,$nxon,$gsplit)=@_; 
  $nxon||=0;
  my $paln= pctof($aln,$tw); # ($tw>0)?int(0.5+ 100*$aln/$tw):0;   
  my $pidn= pctof($idn,$aln); # ($aln>0)?int(0.5+ 100*$idn/$aln):0;  

  #formats? 95%a,3x,1200l; 95a,98i,1200l,3x;  ? + ,Sp:xxx splits
  my $mapval= "${paln}a,${pidn}i,${tw}l,${nxon}x";
  # my $mapval= "${paln}a,${tw}l,${nxon}x";
  # $mapval.= ",${pidn}i" if($pidn < 95);# what level to show? skip?
  $mapval.=",Spl:$gsplit"  if($gsplit);
  return $mapval;
}

# replace ineff over() : BUGGY still, bad calls
sub overgene {
  my($id,$ingr,$go,$gspans)=@_; 
  my $lgr=$ingr;
  # %afidloc{gr}{bn}{ids}, %afexons{id}=[@gx] replace afloc
  return (0,[]) unless($afidloc{$lgr}); 
  
  my($nov,$eqn)=(0,0); 
  my @gx= split",",$gspans; # includes gr's chrs ids with splits
  my %gr;  
  for my $x (@gx) { 
    my($r)= ($x=~m/^(\w[\w\.-]+):/)?$1:$lgr; $lgr=$r;
    my($xb,$xe)= $x=~m/(\d+)\-(\d+)/; 
    my $gb=$gr{$lgr}{gb}||0; $gr{$lgr}{gb}=$xb if($gb==0 or $xb<$gb);
    my $ge=$gr{$lgr}{ge}||0; $gr{$lgr}{ge}=$xe if($xe>$ge);
    } 

  my (%ovd,%lid);
  for my $gr (keys %gr) { 
    for my $bn (glocbins($gr{$gr}{gb}, $gr{$gr}{ge}) ) {
      map{ $ovd{$_}++ if($_ ne $id); } keys %{$afidloc{$gr}{$bn}};
    }
  }
  for my $od (sort keys %ovd) {
    my($ovn1,$eqn1)= overxn($id,$od,$ingr,\@gx,$afexons{$od});
    if($ovn1) { $nov++; $lid{$od}++; }
    ## UPD: opt add report of ovn1,eqn1 per id-od : num exons over, identical
  }
  if($nov) {
    my @lid= sort{ $ovw{$id}{$b}<=>$ovw{$id}{$a} or $a cmp $b} keys %lid; 
    return (scalar(@lid),\@lid); 
  }
  return (0,[]);
}

sub overxn {
  my($id,$aid,$ingr,$gx,$ox)=@_;
  my($nov,$neq)=(0,0);
  return ($nov,$neq) unless(ref $ox);
  ## only have 1 aid here; my %lid; 
  my %didid; my $lgr=$ingr;
  for my $xbe (@$gx) { 
    my($gr)= ($xbe=~m/^(\w[\w\.-]+):/)? $1: $lgr; $lgr= $gr;
    my($xb,$xe)= $xbe=~m/(\d+)\-(\d+)/;
    my $xw=1+$xe-$xb;
    my $algr=$ingr;  # bad default
    for my $ax (@$ox) {  # dont s/// change ax/ox
      my($agr)= ($ax=~m/^([\w\.-]+):/)? $1: $algr; $algr=$agr;
      next unless($algr eq $lgr); 
      my($axb,$axe)= $ax=~m/(\d+).(\d+)/;
      
      my $over= ($xb <= $axe && $xe >= $axb) ? 1 : 0;
      #FIXME# $over=0 if($stranded and ($go =~ /[\+\-]/ and $ao =~ /[\+\-]/ and $go ne $ao)); 
      ## FIXME: add flags for gene equivalence or subset: same or nearly same exons
      if($over) { 
        $nov++;  
        my $ovw= 1+ _min($xe,$axe) - _max($xb,$axb);
        $ovw{$id}{$aid} += $ovw if($ovw>0);
        #o# $neq++ if($ovw >= 0.95*$xw);
        if($ovw >= 0.95*$xw) { $neq++; $eqxon{$id}{$aid}++; } # should use $pEQEXON = 0.95;
        # $neq++ if($xb == $axb); $neq++ if($xe == $axe); # or, if ovw > 0.95* max width?
        #? $didid{$aid.$axb.$axe}++; 
        } 
    }
  }
  return ($nov,$neq);
}

sub over { 
  my($id,$lgr,$go,$gspans)=@_; 
  return (0,[]) unless($afloc{$lgr}); 
  my $nov=0; my $neq=0; my %lid;  my %didid=();
  my @gx= split",",$gspans; # includes gr's chrs ids with splits
  
  ## efficiency change: 
  ## get all gx/scaf, or just genespan, glocbins(genestart,geneend), rather than per-exon
  ## then per-overgene id, test all this.exons x overgene.exons 
  
  for my $xbe (@gx) { 
    my($gr)= ($xbe=~s/^([\w\.-]+)://)? $1: $lgr; $lgr= $gr;
    my($xb,$xe)= $xbe=~m/(\d+).(\d+)/; my $xw=1+$xe-$xb;
    for my $bn (glocbins($xb,$xe)) {
      if(my $aflx=$afloc{$gr}{$bn}) { 
      for my $ax (@$aflx) { 
        my($aid,$ar,$ao,$axb,$axe)=split"\t",$ax;
        next if($aid eq $id); # or  $didid{$ax}
        my $over= ($xb <= $axe && $xe >= $axb) ? 1 : 0;
        $over=0 if($stranded and ($go =~ /[\+\-]/ and $ao =~ /[\+\-]/ and $go ne $ao)); 
        $over=0 if($over and $didid{$ax}); ## {$aid.$ar.$axb.$axe});  
        if($over) { 
          $nov++; $lid{$aid}++; # push @lid,$aid; 
          my $ovw= 1+ _min($xe,$axe) - _max($xb,$axb);
          $ovw{$id}{$aid} += $ovw if($ovw>0);
          if($ovw >= 0.95*$xw) { $neq++; $eqxon{$id}{$aid}++; } # should use $pEQEXON = 0.95;
          $didid{$ax}++; # {$aid.$ar.$axb.$axe}++; # this maybe bug, skipping not xbe
          } 
        } 
      } 
    }
  } 
  if($nov) {
    # my %lid= map{$_,1}@lid; 
    my @lid= sort{ $ovw{$id}{$b}<=>$ovw{$id}{$a} or $a cmp $b} keys %lid; 
    $nov= @lid; return ($nov,\@lid);  # add neq
  }
  return (0,[]);
}

sub altsof { 
  my($td)=@_;
  (my $gd=$td) =~ s/t\d+$//;
  my @aid= grep{$_ ne $td} keys %{$alts{$gd}};
  #xx my @aid= grep{$_ ne $td} grep /^$gd/, keys %afid; ## BAD time sink
  return ($gd,@aid);
}

use constant ALTSPLIT_SKIP => 1; # option?

sub alt_notover_short {
  my($id,$allalts)=@_; 
  # get sizes of shortest/longest alts 
  # .. maybe also check gsplit causes of notover ..
  my @alts= sort keys %$allalts;
  my %aval;
  for my $aid (@alts) { 
    my @aval= split"\t",$afid{$aid};
    # my($bits, $idn, $aln, $tw, $xspans, $gr, $goo, $gloc, $gsplit, $gspans)= @aval;
    $aval{$aid}= \@aval;
  } 
  
  my @asort= sort{ $aval{$b}->[3] <=> $aval{$a}->[3] or $a cmp $b } @alts;
  my $wmax= $aval{$asort[0]}->[3];
  my $wminalt= $wmax * $SHORTALT/100; # what?
  # my $wid = $aval{$id}->[3];
  # my @along= grep{ $aval{$_}->[3] >= $wminalt } @asort;
  my @ashort;
  if(ALTSPLIT_SKIP) {
  @ashort= grep{ $aval{$_}->[3] < $wminalt or $aval{$_}->[8] } @asort;
  } else {
  @ashort= grep{ $aval{$_}->[3] < $wminalt } @asort;
  }
  return @ashort;
}

sub alt_notover { 
  my($id,$ovids)=@_; 
  my ($gd,@alts)= altsof($id);
  if (@alts) {
    my $isover=0;
    my %alts= map{$_ => 1} @alts;
      # FIXME: filter out uninteresting misaltids, ie tiny alts no overlap  
    my @shortalts= alt_notover_short($id,\%alts);
    if(ref $ovids) { map{ $isover++; delete $alts{$_};  } (grep/^$gd/,@$ovids); }
    for my $ad (@shortalts) { delete $alts{$ad}; } ## $isshort++; 
    my @misalt= sort _genealts keys %alts; 
    return (scalar(@misalt),\@misalt,$isover); #?? or what? pull loc,info from afid?
  }
  return (0,0,0);
}

sub maxlist { return join",", (($MAXLIST<1)? @_ :splice(@_,0,$MAXLIST)); }

our %locsame;

sub putover {
  my($td,$oid,$tw,$ovd,$loc, $alnmap, $misaltids)=@_; 
  $oid||="noid"; $alnmap||=0;
  (my $gd=$td) =~ s/t\d+$//;
  
  my $misalt="";
  if(ref $misaltids and @$misaltids > 0) {
    # FIXME: filter out uninteresting misaltids, here or caller, ie tiny alts no overlap
    $misalt="novalt:".maxlist(@$misaltids); # join(",",@$misaltids);
  }
  $misalt||="0"; #"noaltmiss";
  my($nxt)= ($alnmap=~m/,(\d+)x/)?$1:1; 

  my $locsame="";
  my @ovv=();
  if(ref($ovd)) { for my $od (@$ovd) { 
    # my($bits, $idn, $aln, $odw, $xspans, $gr, $gloc, $gsplit, $gspans)=split"\t",$afid{$od}; # dont need here
    if(my $lsid= $locsame{$od}) { $locsame=$lsid unless($locsame); } # Splits problem here
    next if(!$ALTSHOW and $od=~m/^$gd/); # NOALTS
    my $ovt= $ovw{$td}{$od}||0; # bases of od aligned to td
    #my $tov=$ovw{$od}{$td}||0; # bases of td aligned to od, probably same, denom tw or ow changes pctof
    my $pov= pctof($ovt,$tw); # pct of td aligned to od
    # my $pov= pctof($tov,$ow); # want this instead?
    next if($pov < $MINOVERLAP);
if($SHOWEQXON) {
    my $eqx= $eqxon{$td}{$od}||0; my $peqx= pctof($eqx,$nxt); 
    $pov.="xe$peqx" if($peqx>9); # add $peqx if>0 ??
    # $pov.="xeq" if($peqx>95); # add $peqx if>0 ??
}
    push @ovv,"$od/$pov"; # "$od/$pov.$pov"; dont add .pov dupl..
    }
  }
  
  $locsame||=$td; $locsame{$td}= $locsame; # many possible?
  if(ref($ovd)) { for my $od (@$ovd) { 
    $locsame{$od}=$locsame unless($locsame{$od}); # Splits problem here
   } } 
   
  @ovv=("na") unless(@ovv);
  my $oval= maxlist(@ovv); # join(",",@ovv); # limit list to first 10? 20? one alt per locus?
  print $OUTH join("\t",$td,$oid,$oval,$loc,$alnmap,$misalt,"lsd:$locsame")."\n"; # eqgene table
}


sub openRead { # in cdna_evigenesub.pm
  my($fna, $nostdin)= @_; my($ok,$hin)= (0,undef);
  $ok= ($fna =~ /\.gz$/) ? open($hin,"gunzip -c $fna|") 
  	 : ($fna =~ /stdin|^-/) ? *STDIN : open($hin,$fna);  
  # loggit(1,"ERR: openRead $fna") unless($ok);
	die "ERROR: openRead $fna" unless($ok);
  return ($ok,$hin);
}



__END__

=item add to evigene pipeline
  
  add to updated evigene/scripts/prot/tr2aacds2.pl 
  ## UPD 2016.02 option step 6 
  # 6. genome-map reclassing, *after* 1st outclass using okayset
  # make separate pipeline * use after evgmrna2tsa2.pl publicset, for proper alt matching
  my $chrasm=0; # new opt input file -chrasm mygenome.fasta

  sub blastn_cds_genome{};
  sub cdsgmap_maketables{}; includes this overeqcdsloc
 
  if($chrasm) {
    # 6.1. blastn okayset/my.cds to chrasm
    my($cdsgmapblast)= blastn_cds_genome($chrasm, @okayset);
 
    # 6.2. tabulate blast > cdsgmap.tall > cdsgmap.equalgene
    my @xx= cdsgmap_maketables($cdsgmapblast);
 
    # 6.4. classify main/alternate cds, okay & drop subsets, using evigene/rnaseq/asmrna_dupfilter2b.pl
    # combine both cdsblastab + cdsgmapblast.eqgene sources of overlap classing
    ($outclass,$outaln)= asmdupfilter_cds($cdsgmapblast,$aasize,$aaclstr);
 
    # 6.5: new asmdupclass_sum(), asmdupclass_fileset(), ..
  }
  
  
=item input table

ptsize=evg1anofun.cds.qual
$evigene/scripts/makeblastscore3.pl -pctover=0.03 -pmin=0.85 -spans=2 -tall -oneref=1 -aa1 $ptsize \
  genoAfunF1_evg1anofun_cds.blastn > genoAfunF1_evg1anofun_cds.tall3d

cds.tall3d:
trid                  chrid     bits  iden  algn  trlen chrlen(0)
Anofunz4bEVm000224t1	KB668806	12526	6846	6860	6807	0	
  trspans   chrspans  duptrspans (ignore?) dupchrspans
  1-6799/1-105,106-295,296-411,412-611,612-755,756-1023,1024-1113,1114-1326,1327-1546,1547-1766,1767-1946,1947-2204,2205-2332,2333-2479,2480-2891,2892-3504,3505-3693,3694-3804,3805-3958,3959-4290,4291-4485,4486-4594,4595-4914,4915-5069,5070-5258,5259-5428,5429-5617,5618-5766,5767-6799	
  KB668806:290486-290594,290123-290312,289941-290060,288865-289064,288645-288792,288108-288375,287945-288036,287636-287850,287347-287568,287052-287271,284221-284401,283884-284141,283532-283663,217388-217535,216890-217303,216102-216714,215838-216026,215645-215760,215394-215548,214920-215251,214640-214836,214455-214565,214042-214361,213714-213877,213431-213619,212133-212308,211868-212056,211613-211765,210488-211520:-	
    1-1766/1-105,106-295,296-411,412-611,612-755,756-1023,1024-1113,1114-1326,1327-1545,1547-1766	
    KB668874:4186-4294,3823-4012,3641-3760,2565-2764,2345-2492,1808-2075,1645-1736,1336-1550,1050-1268,752-971:-

=item output sample v1
  .. revise, limit overlap/novalt lists to first 10 or so ? or need all?
  .. Splits: is this enough info? need split scaf spans?
  .. Splits: ** strand is wrong? not str for 1st scaf, but last?
  
  egrep 't[123]   noid' genoAfunF1_evg2anofunz4b_cds.eqgenebl | head -30
Anofunz4bEVm000001t1	noid	na	KB668665:428970-516596:-	100%,47370	0
Anofunz4bEVm000001t2	noid	na	KB668665:428970-516596:-	91%,46776	0
Anofunz4bEVm000001t3	noid	na	KB668665:428970-516596:-	91%,46698	0

Anofunz4bEVm000002t1	noid	na	KB669236:155244-223031:-	97%,42531,Sp:1%,KB668471	novalt:Anofunz4bEVm000002t2,Anofunz4bEVm000002t3,Anofunz4bEVm000002t4,Anofunz4bEVm000002t5,Anofunz4bEVm000002t6,Anofunz4bEVm000002t7,Anofunz4bEVm000002t10,Anofunz4bEVm000002t11,Anofunz4bEVm000002t12
Anofunz4bEVm000002t2	noid	na	KB669236:122732-223659:+	100%,31662	novalt:Anofunz4bEVm000002t1,Anofunz4bEVm000002t5,Anofunz4bEVm000002t6,Anofunz4bEVm000002t7,Anofunz4bEVm000002t8,Anofunz4bEVm000002t9,Anofunz4bEVm000002t10,Anofunz4bEVm000002t12,Anofunz4bEVm000002t13
Anofunz4bEVm000002t3	noid	na	KB669236:155244-223031:+	95%,31329	novalt:Anofunz4bEVm000002t1,Anofunz4bEVm000002t6,Anofunz4bEVm000002t8,Anofunz4bEVm000002t9,Anofunz4bEVm000002t10,Anofunz4bEVm000002t12,Anofunz4bEVm000002t13

Anofunz4bEVm000003t1	noid	na	KB668960:1239-36268:-	98%,29910,Sp:61%,KB668589	0
Anofunz4bEVm000003t2	noid	na	KB668960:1239-36268:-	97%,28941,Sp:60%,KB668589	0
Anofunz4bEVm000003t3	noid	na	KB668960:1239-36268:-	97%,28227,Sp:58%,KB668589	0

Anofunz4bEVm000004t1	noid	na	KB668754:240769-344656:+	100%,26484	0
Anofunz4bEVm000004t2	noid	na	KB668754:282358-344656:+	100%,26055	0
Anofunz4bEVm000004t3	noid	na	KB668754:282313-344656:+	100%,26016	0

Anofunz4bEVm000005t1	noid	na	KB668221:1691781-1723022:+	98%,24913	0
Anofunz4bEVm000005t2	noid	na	KB668221:1691781-1723022:+	99%,24562	0
Anofunz4bEVm000005t3	noid	na	KB668221:1691781-1723022:+	99%,24555	0

>> overlapped pair m06,m07
Anofunz4bEVm000006t1	noid	Anofunz4bEVm000007t1/33,Anofunz4bEVm000007t2/29,Anofunz4bEVm000007t3/26	KB668803:784878-812149:-	99%,20457	0
Anofunz4bEVm000006t2	noid	Anofunz4bEVm000007t1/32,Anofunz4bEVm000007t2/29,Anofunz4bEVm000007t3/25	KB668803:784878-812149:-	99%,20241	0
Anofunz4bEVm000006t3	noid	Anofunz4bEVm000007t1/32,Anofunz4bEVm000007t2/29,Anofunz4bEVm000007t3/25	KB668803:784878-812149:-	99%,19947	0

Anofunz4bEVm000007t1	noid	Anofunz4bEVm000006t5/68,Anofunz4bEVm000006t7/52,Anofunz4bEVm000006t6/48,Anofunz4bEVm000006t4/42,Anofunz4bEVm000006t1/33,Anofunz4bEVm000006t2/32,Anofunz4bEVm000006t3/30	KB668803:784878-812147:-	98%,20457	0
Anofunz4bEVm000007t2	noid	Anofunz4bEVm000006t5/54,Anofunz4bEVm000006t7/41,Anofunz4bEVm000006t6/35,Anofunz4bEVm000006t4/28,Anofunz4bEVm000006t1/16,Anofunz4bEVm000006t2/15,Anofunz4bEVm000006t3/13	KB668803:789162-812149:-	83%,16124	0
Anofunz4bEVm000007t3	noid	Anofunz4bEVm000006t5/62,Anofunz4bEVm000006t7/48,Anofunz4bEVm000006t1/40,Anofunz4bEVm000006t6/40,Anofunz4bEVm000006t2/38,Anofunz4bEVm000006t3/36,Anofunz4bEVm000006t4/32	KB668803:789091-808717:-	97%,13317,Sp:2%,KB668781	0

=item paralog examples
 .. diff map loci
Anofunz4bEVm000395t1    noid    na      KB668781:1064699-1097736:+      100%,5676       novalt:Anofunz4bEVm000395t2
Anofunz4bEVm000395t2    noid    na      KB669203:527311-600355:-        99%,4869        novalt:Anofunz4bEVm000395t1

Anofunz4bEVm000398t1    noid    na      KB669469:484224-500531:-        96%,5658        novalt:Anofunz4bEVm000398t2,Anofunz4bEVm000398t3,Anofunz4bEVm000398t4,Anofunz4bEVm000398t5,Anofunz4bEVm000398t6
Anofunz4bEVm000398t2    noid    na      KB668835:81343-87771:-  97%,5421        novalt:Anofunz4bEVm000398t1,Anofunz4bEVm000398t5,Anofunz4bEVm000398t6,Anofunz4bEVm000398t7

-- problem case t3, is it over or not, span suggests yes
Anofunz4bEVm000457t1    noid    na      KB669181:367792-374108:-        87%,5355        novalt:Anofunz4bEVm000457t3
Anofunz4bEVm000457t2    noid    na      KB669181:367792-374108:-        100%,5208       novalt:Anofunz4bEVm000457t3
  KB669181:374027-374108,373370-373552,369944-373181,369732-369856,368181-369520,368008-368119,367792-367938:-  
Anofunz4bEVm000457t3    noid    na      KB669181:369160-373119:-        100%,3675       0
  KB669181:>369944-373119,369732-369856<,369160-369520<:- * should be called overlap to t1,t2; bug where?
  
=item split alt problem case, likely not paralog but may be

Anofunz4bEVm000436t1    noid    na      KB668871:74603-89656:-  100%,5457       novalt:Anofunz4bEVm000436t3,Anofunz4bEVm000436t4
Anofunz4bEVm000436t2    noid    na      KB668871:85670-89656:-  100%,3645       novalt:Anofunz4bEVm000436t3,Anofunz4bEVm000436t4
  >> t3 is split mapped, rev or maybe or split or bug; keep if no other gene equiv at KB668755 ,
Anofunz4bEVm000436t3    noid    na      KB668871:85405-85834:+  87%,1221,Sp:52%,KB668755        novalt:Anofunz4bEVm000436t4


=cut
