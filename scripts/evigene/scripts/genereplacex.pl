#!/usr/bin/env perl
# genereplacex.pl : replace fix.gff into pub.gff w/o changing mRNA, but need checks on changes, misses.

=item about

  genereplacex.pl : replace fix.gff into pub.gff w/o changing mRNA, but need checks on changes, misses.

 -- upd annots, quals: Protein:complete;aaSize=249;cdsSize=76%,750/981
 -- must.tab: update class if fixid class; drop; rename IDs
 ./genereplacex.pl -must mustgene.1211.fixids -fix pub3g.good.cds2aa_fix.gff2 -in pub3g.fix4pan.gff -out pub3h.gff
 
 cacao 2011dec:
 ./replacex.pl -must mustgene.1211.fixids -fix pub3g.good.cds2aa_fix.gff2 -in pub3g.fix4pan2.gff -out pub3h.gff

 wasp 2012jan: 
 $evigene/scripts/genereplacex.pl -dupok -ver 10u -must pub11t.mustfix.1209 -fix pub11t.mustfix.1209.gff \
     -in ../pub11t.pan2.gff -out pub11u.gff
 wasp/genes/
  pub11t_cdserr/pub11t.mustfix.notes
  pub11t_cdserr/pub11t.mustfix.1206
  pub11t_cdserr/pub11t.mustfix.1209.gff
  pub11t_cdserr/pub11t.mustfix.1209
  pub11t_cdserr/pub11t.mustfix.1204.removals
Nasvi2EG011427t1        fix     cdserr.augshort AUGepit5_xfemp2s24g2t1
Nasvi2EG021354t1        drop    cdserr.dupgene  AUGaf_S142g13629t1
Nasvi2EG001614t1        fix     cdserr.stopmiss AUGepit3p2s2g118t1
  ..


=cut

use strict;
use Getopt::Long;

use constant kALT_DEFER => 1;
use constant kALT_ONLY => 2;
my $ANNOTKEY="rx";

my (%opt,%rename,%replace,%drop,%fix,%reclass,%nochange,%unknown,%fixx,
    %strand, %spanfix, %addalt, %mustdup, %mustnote,%did,%mx,%xx,%gefix);

my $optok= GetOptions( \%opt,"must=s", "fix=s", "in=s", "out=s", "vers=s","dupok!");
die unless($optok and $opt{fix} and $opt{in});
my $vers= $opt{vers} || "1";
#my $amark="rx=$vers";

## instead of handling out-of-order records here, run thru 
#  bestgenes_update.pl -conf evigene_cacao3.conf -act=sort -in xxx.gff -out xxx.sgff

# add actions: spanfix, gene add, defer/deferto .. from gff2genbanktbl
## gene action for out-of-order mRNA, these 2 come before gene ID=Thecc1EG047060
# Thecc1EG047060t3  gene  Thecc1EG047060  
# Thecc1EG047060t2  gene  Thecc1EG047060
# .. handle as defer mRNA tID till see gene gID ?
## defer action: Thecc1EG047019t2 is out of order FOLLOWING Thecc1EG047019t1 (correct, renamed from Thecc1EG047020t2)
# Thecc1EG047019    defer 
# # Thecc1EG047020t1 in between t1,t2
# Thecc1EG047019t2  deferto Thecc1EG047019

## warn "#FIXME: must update mRNA spans to match new exons\n";

if($opt{must}) { open(M,$opt{must}); while(<M>) {
  next unless(/^\w/); s/#.*$//;
  my $ce = (s/^CDSerr\s+//) ? 1 : 0;
  my($id,$act,$nid,@mmust)=split;
  if(0 < $mustdup{$id}++) { } # warn or die? some ok, like drop/replace, xxx/nochange, ..
  if($act =~ /^rename/){ $rename{$id}=$nid; }
  elsif($act =~ /^replace/){ $replace{$id}=$nid; } # needs new gff in -fix gff ** same/nearly as fix
    # replaceall means dump entire old rec including mRNA
    # replace may require new mRNA line, fix does not. BUT keep id=Thecc, nid here is oid not new id
    # add this syntax? Nasvi2EG012349t1,Nasvi2EG012351t1,Nasvi2EG012352t1 replace AUGepi4s29g55t1 
    #  .. means replace id1, drop rest; ie genefrags > 1 gene; or id1,id2,.. drop = list 
  elsif($act =~ /^drop/) { $drop{$id}=1; }
  # elsif($act =~ /^ignoredrop/) { delete $drop{$id}; } # too tricky; remove drop action
  elsif($act =~ /^fix/) { $fix{$id}=1; }
  elsif($act =~ /^addalt/){ $addalt{$id}=$nid; } # needs new gff in -fix gff 
  elsif($act =~ /^addgene/){ $addalt{$id}=$nid; } # needs new gff in -fix gff 
    ## addgene = newgene, handle like addalt: start id to add after
    # addalt: mainid addalt nid=altid oid < need both nid, oid to rename?
    # Nasvi2EG001908t1  addalt Nasvi2EG001908t2		oid=nvit1v3S2big1Loc1699t5	mustpick1201	altof
    # handle altid addalt/altof  oid  = new alt in fix.gff; not replace as no equiv mRNA
    # handle newid newgene oid = new gene in fix.gff; not replace, no equiv mRNA
  elsif($act =~ /^nochange/) { $nochange{$id}=1; } # special case such as splitgene, original is correct 
  elsif($act =~ /^strand/ and $nid =~ /^[+-]$/) { $strand{$id}=$nid; }  
  
## bad mRNA spans after rx=3 update == MISSING_GENES
# Thecc1EG002023t1 spanfix scaffold_1:11049068-11058532:-
  elsif($act =~ /^spanfix/ and $nid =~ /^(\w+):(\d+)[.-]+(\d+)[:]?([+-]?)$/) {
    my($nr,$nb,$ne,$no)= ($1,$2,$3,$4);
    # fix span only, not strand/ref; presumably to match exons span
    if($ne>0 and $nb>0) {  $spanfix{$id}="$nr\t$nb\t$ne\t$no"; }
    } 

  elsif($act =~ /class=(\S+)/) { my $c=ucfirst($1); $reclass{$id}=$c; }
  # exonedit now unknown; all have also fix action
  else { $unknown{$id} .= $act.","; }
  $#mmust=2 if(@mmust>3); # chop excess notes
  my $n= join(".", $act,$nid,@mmust )."\t"; # record in output attr?
  unless($act =~ /class=/) { $n =~ s/=/:/g; $mustnote{$id} .= $n; }
}  close(M); }


open(F,$opt{fix}) or die;
while(<F>) { next unless(/^\w/);
  my $k=(/\t(mRNA|gene)/)?"ID":"Parent"; 
  my ($d)=m/$k=([^;\s]+)/; 
  next if($nochange{$d}); # leave alone, ignore any fixes
  
  ## FIXME: 1208: gene span updates w/ fix of mrna,exons ***  genex{$d} ??
  if(/\tmRNA/){ $fixx{$d}++; 
    # need to resolve, remove duplicate fixes or will get wrong one
    # also need to check dup replace/rename ID,    
    ## splitgene.fix here has 2+ mrna
    # check mustnote{$d} =~ /replaceall/ ?? otherwise die
    if($fixx{$d}>1) { 
      my $dupok= ($opt{dupok} and ($mustnote{$d} =~ /replaceall/))?1:0;
      ## $mx{$d} .= $_; # << bad for splitgene; stutter.
      if($dupok) { $xx{$d}.=$_; $did{$d}{dupfixERR}++; } # put 2nd mRNA w/ exons? keep in order
      else { die "# *** duplicate fix.gff $d\n"; } 
    } else { $mx{$d}=$_; } 
  } elsif(/\tgene/) {
    $gefix{$d}= $_; # 1208 fixme
  } elsif($d) { $xx{$d}.=$_ unless($fixx{$d}>1 and not $opt{dupok}); } 
} close(F);

my $out= *STDOUT;
if($opt{out}) { open(OUT,">$opt{out}") or die; $out= *OUT; } 
open(IN,$opt{in}) or die; # no .gff.gz

my($repid,$nid,$addout,$outrec,$outdefer,$xadded) = (0) x 10;  

while(<IN>) {
  unless(/^\w/){ print $out $_; next; }
  $repid=$nid=$addout=""; 
  my $k=(/\t(mRNA|gene)/)?"ID":"Parent"; 
  my ($id)=m/$k=([^;\s]+)/; 
  my $tid=$id; $tid.="t1" if(/\tgene/);
  my $nochange=0;
  if($nochange{$id}) { $nochange=1; $did{$id}{nochange}++; }
  else {
    if($drop{$id} and not $replace{$id}) {  $did{$id}{drop}++; next; } # FIXME: drop/replace actions occur .. 
    if($nid=$rename{$id}) {  
      if(s/$k=$id/$k=$nid/) {  # fixme in handleMRNA or here: rename mRNA gene=oldID > gene=newID
        $did{$id}{rename}++; $did{$id}{renameTO}=$nid; 
        if(/\tmRNA/) { my($ogid,$ngid)=($id,$nid); map{ s/t\d+$// }($ogid,$ngid); s/;gene=$ogid/;gene=$ngid/; }
      } else { $did{$id}{renameERR}++; }
    }
    if($repid=$replace{$id}) {  # NO reversed case, d==oid,oldid
      #NO, if(s/$k=$id/$k=$nid/) { $did{$id}{replace}++; $did{$id}{replaceTO}=$nid; } else { $did{$id}{replaceERR}++; }
    }
    if(m/\tmRNA/ and (my $span=$spanfix{$tid})) { # or do in handleMRNA
      my($nr,$nb,$ne,$no)=split"\t",$span; my($r,$b,$e,$o)=(split"\t")[0,3,4,6];
      s/\t$b\t$e\t/\t$nb\t$ne\t/; $did{$id}{spanfix}++;  #  change all gene/mRNA/exon/CDS to this 
    }
    if(my $xor=$strand{$tid}) { my($b,$e,$v)=(split"\t")[3,4,5]; s/($b\t$e\t$v\t)./$1$xor/; $did{$id}{strand}++; } #  change all gene/mRNA/exon/CDS to this 
    elsif(/\tgene/ and (my $xor=$strand{$id})) { my($b,$e,$v)=(split"\t")[3,4,5]; s/($b\t$e\t$v\t)./$1$xor/; $did{$id}{strand}++; } #  change all gene/mRNA/exon/CDS to this 
  }
  
  
  # my $altid= $addalt{$id}; #? assumes none of other changes apply?
  # my $storeid= ($nochange{$id}) ? 0 : ($repid) ? $repid : ($altid) ? $altid : $id;
  # my $storeid= ($nochange{$id}) ? 0 : $repid || $altid || $id;

  if(/\tmRNA/) { 
    print $out $outdefer if($outdefer); $outdefer="";
    
    if($addalt{$id} and ( $replace{$id} or $fix{$id}) ) {
      # cant do both with this handleMRNA..
      ($outrec, $outdefer,$xadded)= handleMRNA( $id, $_, kALT_DEFER); # for alts, want to return record here, defer print till AFTER this record
      print $out $outrec;
      my($Aoutrec, $Aoutdefer, $Axadded)= handleMRNA( $id, $_, kALT_ONLY); # for alts, want to return record here, defer print till AFTER this record
      $outdefer= $Aoutdefer;
    } else {
      ($outrec, $outdefer,$xadded)= handleMRNA( $id, $_); # for alts, want to return record here, defer print till AFTER this record
      print $out $outrec;
    }

  } elsif(/\tgene/) {  # 1208 fixme : 
    # FIXME2, this fails when only alt mRNA t2.. is updated w/ new gene span
    # $tid.="t1" if(/\tgene/); above is problem .. test all t1..t99 ? 
    # ?? ignore tid and replace if gene gefix{id} ? not all input gefix{id} need changes, may be partial/false?
    
    if( $gefix{$id} ) {
      my $ok=0; my $gefix= $gefix{$id};
      if($gefix eq $_) { $ok=1; } # ok
      else { 
        for my $t (1..49) { my $tid= $id."t$t"; if($replace{$tid} or $fix{$tid}) { $ok=1; last; } }
        if($ok) { $_ = $gefix; $did{$id}{fixgene}++; } # 3 cases of alt t2.. here
        else { print "#RX.warn.gefix not applied\n#RX.gefix:$gefix"; } 
      }
    }
    print $out $_;

  } elsif(/\t(exon|CDS)/) {
    # my $xadded= ($nochange) ? 0 : ($repid)? $xx{$repid} : $xx{$id};
    #>> this is bad, use handleMRNA() return to say if exons need print
    ##my $xadded= $xx{$storeid}; # ($repid)? $xx{$repid} : $xx{$id};
    
    if($xadded) { $did{$id}{fixexon0}++;  }
    else { print $out $_; } # NOT here: print $out $addout if($addout);
    
  } else { # anything else, print
    print $out $_;
    #NOT here# print $out $addout if($addout); 
  }
  
} close(IN);
print $out $outdefer if($outdefer);

# print stats of work..
foreach my $d (sort keys %did) {
  my @act=sort keys %{$did{$d}}; 
  my @v= map{ "$_:$did{$d}{$_}" } @act;
  print $out join("\t","#RX",$d,@v)."\n";
}

#  check misses in %fixx, %drop %rename %reclass?
sub pute { my($d,$ac,$er)=@_; $er||="MISS"; print $out join("\t","#RX",$d,$er,$ac,$mustnote{$d})."\n" unless($did{$d}{$ac}); }
foreach my $d (sort keys %fixx) { pute($d,"fix") unless($did{$d}{drop}); } 
foreach my $d (sort keys %drop) { pute($d,"drop"); } 
foreach my $d (sort keys %rename) { pute($d,"rename"); } 
foreach my $d (sort keys %addalt) { pute($d,"addalt"); } 
foreach my $d (sort keys %replace) { pute($d,"replace"); } 
foreach my $d (sort keys %reclass) { pute($d,"reclass"); } 
foreach my $d (sort keys %nochange) { pute($d,"nochange"); } 
foreach my $d (sort keys %unknown) { pute($d,"unknown") unless($did{$d}{fix}); }  # exonedit here

close($out);

#.........

sub handleMRNA
{
  my($id,$inmrna, $altaction)= @_;

  my $inid=$id; # for alt change
  $altaction ||= 0;
  my $nid   = $rename{$id};
  my $repid = $replace{$id};
  my $exonsadded= 0;
  
  # ** FIXME: have gene with replace + addalt; handle both here? need loop..
  my $nochange= ($nochange{$id}) ? 1 : 0; #  annotate this
  
  # for addalt, input id/mrna is record to put this AFTER, with its own alt id..
  my $altid  = ($altaction == kALT_DEFER) ? 0 : $addalt{$id}; #? assumes none of other changes apply?
  my $storeid= ($altaction == kALT_ONLY) ? $altid 
    : ($nochange) ? 0 : ($repid) ? $repid : ($altid) ? $altid : $id;
  
  my $mx= ($storeid) ? $mx{$storeid} : 0; # fix{id} here??
  my $addout="";
  
  local $_= $inmrna;
  if($altid eq $storeid) {
    my($ref,$msrc,$mb,$me,$mv,$mor)= (split"\t")[0,1,3,4,5,6]; #my @mv = split"\t"; 
    $_= $mx;
    s/^$ref\t\S+/$ref\t$msrc/; # retain old src field
    $id= $altid;
    $nid= $rename{$altid}; # $altnid;
    if($nid) { my $k="ID"; 
      if(s/$k=$id/$k=$nid/) { 
        unless(s/;oid=[^;\s]+/;oid=$id/) { s/$/;oid=$id/; } 
        $did{$id}{rename}++; $did{$id}{renameTO}=$nid;  
      } else { $did{$id}{renameERR}++; } 
    }
  }
  elsif($repid eq $storeid and $mustnote{$id} =~ m/replaceall/ and $mx) {
    $_= $mx; 
    my $didall= ++$did{$repid}{replaceall};
    # new method: 1st call here replaces all parts of splitgene from mx{repid}, xx{repid}
    # 2nd+ calls here return empty output, ie delete all subsequent parts
    if($didall>1) { return ("","",1); } # $exonsadded=1;  
    $did{$id}{replace}++ if($id ne $repid); 
    $id= $repid; $nid= $rename{$id};
    # split genes still big problem; each part calls here w/ same ID; gets/outputs same replacement
    #?? $id= $repid; $nid= $rename{$id};  # replaceall presumes all of new record is correct
    ## $mx=""; # no, need to do addout xx{storeid}
  }
  
  if($mx) {
    ## FIXME: update mRNA strand if needed
    ## FIXME: replace (not all), copy quality= if repmRNA has expertchoice ? add to reclass{} if missing?
    
    my @mv = split"\t"; 
    my($ref,$msrc,$mb,$me,$mv,$mor)= @mv[0,1,3,4,5,6];
    my($xref,$xor,$mxan)= (split "\t",$mx)[0,6,8];
    if($mor eq "." and $xor ne ".") { s/($me\t$mv\t)$mor/$1$xor/; $did{$id}{strand}++;}
    if($repid) { s/;oid=[^;\s]+/;oid=$repid/ unless($repid eq $id); $did{$repid}{fix}++; $did{$id}{replace}++; $did{$id}{drop}++ if($drop{$id});}
    
    ## FIXME add skip prot flag; only exon changes
    # NO: not all fix updates here have aalen= cxlen= protein= to replace old
    my $aanochange= 0;
    $aanochange=1 if($mustnote{$id} =~ m/keepaa/);
    
    my($aan,$as,$ap,$aqual,$cx);
    if($aanochange) {
      # ($aan)= $_ =~ m/protein=([^;\s]+)/;     
    } else {
      ($aan)= $mxan =~ m/protein=([^;\s]+)/; 
      s/protein=[^;\s]+//; s/$/;protein=$aan/ if($aan);
      
      # cxlen=909/1413;aalen=303,64%,complete;
      # my($aal)= $mxan =~ m/aalen=([^;\s]+)/; 
      ($as,$ap,$aqual)= $mxan =~ m/aalen=(\d+),(\d+%)([^;\s]*)/;
      ($cx) = $mxan =~ m/cxlen=([^;\s]+)/;
      # ^^ FIXme: aaSize=;cdsSize=; bad data
    }
    
    # FIXME: missing/wrong gene= for some; addalt/addgene, other?
    # use bestgenes_update.pl -sort -pubid?/-addgene
    
    $addout= $xx{$storeid}; # ($repid)? $xx{$repid} : $xx{$id};
    my @addout= split"\n",$addout;
    my $nxa= scalar( grep /\t(exon|CDS)/, @addout); # check new vs old exon count
    ## FIXME2: update mRNA span if needed to new xspan
    my($cw,$xw,$tb,$te)=(0) x 9;
    foreach my $x (@addout) { 
      my @x=split"\t",$x; my($xt,$xb,$xe)=@x[2,3,4];
      my $w=1+$xe-$xb;
      if($xt =~ /exon/){ $xw+=$w; $tb=$xb if($tb==0 or $tb > $xb); $te=$xe if($xe>$te); }
      elsif($xt =~ /CDS/){ $cw+=$w; }
      }
    unless($mb == $tb and $me == $te) {  s/\t$mb\t$me\t/\t$tb\t$te\t/; $did{$id}{span}++;} 
    
    unless($aanochange) {
    ## ?? need cw -= 3 for stop codon?
    if($aan) {       
      $as=length($aan); $as-- if($aan=~m/\*$/);  # even if $as recalc; add complete/partial flag
      $aqual="complete" if($aan=~m/\*$/ and $aan=~m/^M/);
    }
    $ap= ($xw>0) ? int(0.5 + 100*$cw/$xw) : "0";
    $cx="$ap%,$cw/$xw";
    
    if($mxan =~ m/quality=Class:([^;\s]+)/) { my $mq=$1;
      $reclass{$id}= $mq if($mq =~ /expert/ and not $reclass{$id});
    }
    
    # fixme quality=StrongPartial-expertchoice .. Protein:complete ..
    s/aaSize=[^;\s]*/aaSize=$as/; 
    s/cdsSize=[^;\s]*/cdsSize=$cx/; # or add?   
    if($aqual =~ /complete/ and /Protein:\w*partial/) {
      s/Protein:\w*partial\w*/Protein:complete/;
      s/(quality=Class:\w+)Partial/$1/;   # quality=Class:StrongPartial
    }    
    } # aanochange.................
    
    # fixme addalt/addgene: msrc is wrong for mRNA + exons
    map{ s/$ref\t\S+/$ref\t$msrc/; 
      s/Parent=$repid/Parent=$id/ if($repid);
      s/Parent=$id/Parent=$nid/ if($nid); 
    } @addout;
    $addout= join"\n", @addout; $addout .="\n";
    
    ## FIXME: need some other mx annots copied over: cdnabest= ; ocds=? other? fxs=augnocds/etcerr2/..
    # ocds=NEn,NEaa:+77;oaaln=343,43%,oomplete,NE0,NE1,NE2,NE3;
    # cdnabest=aaNE:+77,gdNE:+306,dfull:0,trNE:-8,xeq1:1-868,xeq2:869-1217,xeq3:1218-1401,xne4N71,xeq5:1739-2358g337,tN:71,tcov:97,nx0:5;
    # homolog=675,apis2ref:XP_623782.1;fxs=fix8cdbesteq2
    
    ## FIXME: mark annot change rx=1 in mRNA line; in exon lines too? no
    my @fan=();
    foreach my $fk (qw( fsx ocds cdnabest )) { # also? cdnabest homolog=
      my ($an)= $mxan =~ m/($fk=[^;\s]+)/; 
      if($an and $fk eq "cdnabest") { 
        if($mxan =~ m/$an.homolog=([^;\s]+)/) { my $h=$1; $an .= ",hoscore:$h"; } # missing still
        $an =~ s/,(nx0|dfull|tcov|xeq\w*|xne\w*):[^,;\s]+//g; $an =~ s/,(xne\w*)//g;
        $an =~ s/,tN:/,gapspan,tN:/;
        }
      push @fan, $an if($an); 
      }
    if(my $mnote= $mustnote{$id}) { $mnote =~ s/\t/,,/g; $mnote =~ s/,*$//; push @fan, $ANNOTKEY."note=$mnote"; }
    my $fan= (@fan) ? ";".join";",@fan :"";
    s/$/;$ANNOTKEY=$vers$fan/;  
    $did{$id}{fix}++; $did{$id}{fixexon1}= $nxa;

  } elsif($nochange) {
    s/$/;$ANNOTKEY=$vers;fxvalid=1/;  # want to mark in genbanksubmit as exceptions where CDS <> prot. eg splitgene, cdnabest, 
  }
     
  if(my $rc=$reclass{$id}) { $did{$id}{reclass}++; 
    unless(/quality=Class/) { s/$/;quality=Class:$rc/; }
    elsif($rc =~ /,Protein/) { s/Class:[^;\n]+/Class:$rc/; }
    else { s/Class:[^;,\n]+/Class:$rc/; }
  }
    
#     print $out $_;
#     print $out $addout if($addout); 
  my($outrec,$outdefer);
  if($altid eq $storeid) {
    $outdefer= $_ . $addout;
    $outrec  = $inmrna; # unchanged from input ??
    $did{$id}{addalt}++; $did{$inid}{addalt}++;  
  } else {
    $outrec= $_ . $addout;  $exonsadded=1 if($addout);
    $outdefer="";
  }
  return ($outrec,$outdefer,$exonsadded); # or print
}
