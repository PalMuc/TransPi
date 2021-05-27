#!/usr/bin/env perl
# locussamediff2pubid.pl

=item usage

 cat $nam.diffloci.eqalt2main2 | egrep ' (Df1|Df2) ' | sed 's/^/difflocus  /;' | \
  cat - $nam.sameloci.eqall $nam.newids | perl locussamediff2pubid.pl
 
 ## add new inputs: 
  #MISS id, trclass info, from last pass of splitmainalt.pl
  #VALID oid .. from input.gff
 
 Results:
  322715 kf2pub11_both.newids     input  (-166531 drops)
  322715 kf2pub11_both.renewid9b  output (excludes drops)
    annots added: dfmap (diffloci), samap (sameloci),
      reorder_oldid, missid, ..

Kfish2 diffsame alt reclass using genome map results: newids > renewid9b
  TRclass  all         nochange     diffmap    samemap
  althi1  110744     93%,103435    6%,6879    0.3%,430
  althi2  74352      88%,65608     11%,8593   0.2%,151
  altmid  27791      81%,22435     19%,5234   0.4%,122
  altmap  20989      77%,16168     21%,4351   2.2%,470
  main    52638      99%,52345     0%,0       0.5%,293
  noclass 36200      99%,35871     0%,0       0.9%,329
  total   322714     92%,295862    8%,25057   0.5%,1795

       
=cut

while(<>) {

## FIXME: make td,ad t-num sortable with t02 > t10 ..  _sortid
## add new /^missid/ from splitmainalt #MISS: records,
  if(/^diffloc/) { 
    my($xx,$td,$od,$otherd,$loc,$cla,$mloc)=split; 
    $otherd=~s/,.*//;
    $tdiff{$td}=join("\t",$loc,$otherd,$cla,$mloc); 
    }
  elsif(/^notsame/) { } # skip, from sameloci.tab
  elsif(/^sameloc/) { # new variant: ^(samelocus|samelocsw|notsamer)
    ## FIXME: problem here td,ad ordering is not accurate to say which is 1st gene to get 2nd as new alt
    ## .. need 1. aasize, 2. aaref, of other alts for td, ad to decide?
    my($xx,$td,$ad)=split; ($gd=$td)=~s/t\d+$//; 
    #x push @{$tsame{$td}},$ad; # push @{$gsame{$gd}},$ad;
    ## next if($tdiff{$td}); #?? is this right way to break circular assigns? assume diffloc filled 1st
    $tsame{$td}{$ad}++ unless($tdiff{$td}); ## tsame,gsame : hash ad so no dups
    } 
  elsif(/^#MISS|^MISS/i) {
    #MISS:  Funhe2Exx11m000617t1    Funhe2E6bm000599t1      Funhe2Exx11m000617      1       althi1  1905,73%,complete       99/100/./Funhe2Eq7m083795t3
    #       aaref:1905,UniRef50_A3KMH1,reorder_oldid:Funhe2Exx11m000617t2,
    s/^\S+\s+//; chomp; @v=split"\t"; my($td,$oid,$gd,$ti,$cla)=@v;
    $missid{$td}="$oid,$gd,$ti"; #? $missinfo{$td}= $_;
    } 
  elsif(/^#VALID|^VALID/i) {
    s/^\S+\s+//; my($okoid,@otherstuff)=split; $validoid{$okoid}=1; $nvalidoid++; 
    }
  elsif(/^\W/) { print; } 
  else { 
    chomp; @v=split"\t"; $td=$v[0]; ($gd=$td)=~s/t\d+$//; 
    next if($v[4]=~/^drop/); #? bad
    $tinfo{$td}=$_; 
    $glist{$gd}{$td}++; #o push @{$glist{$gd}},$td; 
    } 
}

sub _sortid {
  my($ida,$idb)= @_;
  my($ga,$ia)= $ida=~m/(\w+)t(\d+)$/; # split "t",$ida;
  my($gb,$ib)= $idb=~m/(\w+)t(\d+)$/; # split "t",$ida;
  return ($ga cmp $gb or $ia <=> $ib);
}

sub diffloc {
  my($aloc,$bloc)= @_;
  my($ar,$ab,$ae)= split/[:-]/, $aloc;
  my($br,$bb,$be)= split/[:-]/, $bloc;
  return ($ar ne $br)?2:($ae < $bb or $ab > $be)?1:0;
}
  
END {

  @gids=sort keys %glist; $newg=$gids[-1]; 
  ($idp,$newgn)= $newg =~ m/^(\w+\D)(\d+)$/;

  ## FIXME: too many newgn, where group of alts from oldgn move to new locus, 
  ##   each made new locus but map to same place.
  ##   .. need some locus compare for tdiff set, esp by gd groups.
  
  foreach $td (sort _sortid keys %tdiff) { ## got duplicates, same oid ??
    my($loc,$otherd,$cla,$mloc)=split"\t",$tdiff{$td};
    $fo=$tinfo{$td} or next; @fo=split"\t",$fo; 
    next if($missid{$td}); #??
    $oid=$fo[1]; $gdold=$fo[2];
    $otherd="na" unless($otherd eq "na" or $tinfo{$otherd});
    
    if($otherd eq "na") { 
      #o $newgn++; $nig=$idp.$newgn; $ti=0; ## $ti=1; $nid=$nig."t1"; 
      ## FIX here test loc and gene group of td for bundling at same newgn
      if($old2newgn{$gdold}) {  
        foreach my $otd ( @{$old2newgn{$gdold}} ) { 
          my($oloc)= split"\t",$tdiff{$otd}; 
          unless(diffloc($loc,$oloc)) { 
            $nig=$newgn{$otd}; $ti=scalar( keys %{$glist{$nig}});  
            last; }
          }
      } else { $newgn++; $nig=$idp.$newgn; $ti=0; }
      
      do { $ti++; $nid= $nig."t$ti"; } while($tinfo{$nid});   
      push @{$old2newgn{$gdold}}, $td;
      $newgn{$td}= $nig;
      
    } else { 
      ($og,$ot)= $otherd=~m/^(\w+)t(\d+)$/; 
      $nig=$og;  $ti=scalar( keys %{$glist{$nig}});
      do { $ti++; $nid= $nig."t$ti"; } while($tinfo{$nid});   
    }
    $fo[0]=$nid; $fo[2]=$nig; $fo[3]=$ti; $fo[4].="dfmap"; $fo[7].="mapdiff_oldid:$td,";
    $tinfo{$nid}=join("\t",@fo); $tinfo{$td}=""; 
    $glist{$nig}{$nid}++; #opush @{$glist{$nig}},$nid;
    delete $glist{$gdold}{$td}; #o @gl1= grep{ $_ ne $td } @{$glist{$gdold}}; $glist{$gdold}= [@gl1];
  }
  
  # FIXME: tsame id1, id2 ordering not accurate. Use tinfo? aasize,aaref to reorder
  # fix in sameloci input table
  ## ($td,$od,$gd,$ti,$cla,$aaq,$pia,$notes) = split"\t", tinfo
  foreach $td (sort _sortid keys %tsame) { 
    next unless($tinfo{$td});
    next if($missid{$td}); #??
    ($gd=$td)=~s/t\d+$//; 
    $ti= scalar( keys %{$glist{$gd}});
    @ad = sort _sortid keys %{$tsame{$td}}; ## @ad= @{$tsame{$td}}; 
     
    foreach $ad (@ad) { 
      $fo=$tinfo{$ad} or next; @fo=split"\t",$fo; 
      next if($missid{$ad}); #??
      do { $ti++; $nid= $gd."t$ti"; } while($tinfo{$nid});   
      $oid=$fo[1];
      $fo[0]=$nid; $fo[2]=$gd; $fo[3]=$ti; $fo[4].="samap"; $fo[7].="mapsame_oldid:$ad,";
      $tinfo{$nid}=join("\t",@fo); $tinfo{$ad}=""; 
      $glist{$gd}{$nid}++; #o push @{$glist{$gd}},$nid;
      ($ag=$ad)=~s/t\d+$//; 
      delete $glist{$ag}{$ad}; #o @gl1= grep{ $_ ne $ad } @{$glist{$ag}}; $glist{$ag}= [@gl1];
      }
  }

  my %dido=();
  foreach my $gd (sort keys %glist) {
  
    ## fixme: number sort @td by t-num
    @td= sort _sortid keys %{$glist{$gd}}; # @td=  @{$glist{$gd}};
  
   # FIX4: reorder td alts by aa/aaref
    my (%gtaw,%gtho);
    foreach my $td (@td) {
      $fo=$tinfo{$td} or next;
      my($tdd,$oid,$gd,$ti,$cla,$aaq,$pia,$nos)= @fo=split"\t",$fo; 
      my($aw)=$aaq=~m/^(\d+)/; my($ho)=$nos=~m/aaref:(\d+)/; 
      if($nvalidoid) { $aw=0 unless($validoid{$oid}); }
      else { $aw=0 if($missid{$td}); } #?? will this force to last of alts      
      $gtaw{$ti}=$aw;  $gtho{$ti}=$ho||0;
    }
 
    my($hb,$ha,$wb,$wa,$whb,$wha);
    my @ti= sort{ $a <=> $b } keys %gtaw;
    my @hi= sort{ 
      $hb=$gtho{$b}; $ha=$gtho{$a}; $wb=$gtaw{$b}; $wa=$gtaw{$a}; 
      $whb=$wb*int($hb/10); $wha=$wa*int($ha/10); 
      $whb <=> $wha  or $a <=> $b; } keys %gtaw;
    my $oord=($ti[0] == 1)? 0: 1; $ni=$#ti; $ni=5 if($ni>5); 
    for my $i (0..$ni) { $oord++ unless($ti[$i] == $hi[$i]); }
    if($oord) { @ti=@hi; } # flag it?
    
    my $iord=0; 
    foreach my $ti (@ti) {
      $iord++; 
      my $td= $gd."t$ti"; # need this? err?
      my $nid= $gd."t$iord";
      my $fo=$tinfo{$td} or next; 
      @fo=split"\t",$fo; my $otd=$fo[0];
      ## FIXME: renumber ti info record; need oldid: newid changes also
      if($oord and $nid ne $td) {  # or $nid ne $td ??
        $fo[0]=$nid; $fo[3]=$iord; 
        unless($fo[7] =~ /oldid:/) { $fo[7].="reorder_oldid:$td,"; }
        $fo=join("\t",@fo); $td=$nid;
      }
      
      $fo.="missid:$otd," if($missid{$otd});   
      if($dido{$td}++) { $fo="#ERR_DUP:$fo"; }
      print $fo,"\n"; 
    }
  }
    
}

=item fixmes

  FIXME: need to check diffloc vs sameloc for circular re-mappings ..
  i.e. diffloc for one case can be sameloc for other.
  e.g Funhe2Exx11m000009t / Funhe2Exx11m000013
  g009 got 4 sameloc from 013 including main:Funhe2E6bm000010t1, 
  g013 got 4 diffloc from 009, 3 are Funhe2E6bm000010t* .. are same set
  in this case, all should be diffloc move to Funhe2Exx11m000013t locus, as g009 set are smaller t18..37 alts
  .. maybe diffloc always superceed sameloc ?
  
  FIXME2: poor mapping mixups, esp sameloc from Split genes, low qual aligns .. drop those ?
  FIX: check for A < B ; B < C assign mixups

  FIX3: merge missing mains, via trclass1 changes
   -- input here missmain2a.trclass1 replacements for newids, or change .newids first

cat missmain2a.trclass1 kf2pub11_both.newids | perl -ne\
'if(/^Funhe2Exx|^#Pub/){ ($pd,$od)=split; if($ad=$oad{$od}) { s/[\.]$//; s/$/oidold=$od,/; s/\t$od/\t$ad/; } print; }\
elsif(/\tdrop/) { ($od,$drp,$cla,$ad)=split; $oad{$od}=$ad; }' \
  > kf2pub11_both.newid2


  cat missmain2a.trclass1 kf2pub11_both.renewid6 | perl -ne\
'if(/\tdrop/) { ($od,$drp,$cla,$ad)=split; $oad{$od}=$ad; } \
else { ($pd,$od,$gd)=split; if($ad=$oad{$od}) { s/[\.]$//; s/$/oidold=$od,/; s/\t$od/\t$ad/; }\
 print; }' > kf2pub11_both.renewid7

  FIX4: Renumber miss-main alt2.., Reorder alts by aasize/aaref weight
  cat kf2pub11_both.renewid7 | grep -v '^#' | perl -ne \
'($td,$od,$gd,$ti,$cla,$aa,$pia,$nos)=split"\t"; $gti=$gd."t$ti";
$err=""; if($gti ne $td) { $err.="$td ne $gd,$ti"; }
unless($did{$gd}++) { $err.="alt>main $gd,$ti" unless($ti==1); }
($aw)=$aa=~m/^(\d+)/; ($ho)=$nos=~m/aaref:(\d+)/; $gtaw{$gd}{$ti}=$aw;
$gtho{$gd}{$ti}=$ho||0; $gter{$gd}{$ti}=$err; $xprint= "$td\t$err\n"
if($err); END{ foreach $g (sort keys %gtaw) { @ti=sort{$a<=>$b}keys
%{$gtaw{$g}};  @hi=sort{ $hb=$gtho{$g}{$b}; $ha=$gtho{$g}{$a};
$wb=$gtaw{$g}{$b}; $wa=$gtaw{$g}{$a}; $whb=$wb*int($hb/10);
$wha=$wa*int($ha/10); $whb <=> $wha  or $a <=> $b; } keys
%{$gtaw{$g}}; $err=""; for $ti (@ti) { $e=$gter{$g}{$ti}; $err.="$e,"
if($e); } $oord=0; $ni=$#ti; $ni=5 if($ni>5); for $i (0..$ni) {
$oord++ unless($ti[$i] == $hi[$i]); } if($err or $oord) { print
"$g\terr:$err\tord:@hi;old:@ti\n"; } } } ' \
 >  kf2pub11_both.renewid7.reord 
 
 nmiss main=464
 nreord = 20556, n-new-t1=12622
 
=cut

