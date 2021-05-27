#!/usr/bin/env perl
# evigene/scripts/genes/overselfdroptab.pl

=item about

  $evigene/scripts/genes/overselfdrop.pl  geneset.ovself.eqgene >  geneset.ovself.keepdrop
  opts: 
    env gscore=genescore.table ( ID gscore=999 ... )
    env mina=20 .. 10 : cds %align limit
    env gdup=0.50  : multimap gene duplicate down-weight, score below uniq/1st gmap gene
    
  use after equalgene, note cmerge.gff needs IDs of main/t1, alt/tn pattern to resolve loci
  $evigene/scripts/equalgene.pl -selfsame -this -oneexonstrandless -in geneset.gff \
	  > geneset.ovself.eqgene

  from zfish/sra17zmerge6evg/map5overselfloc.pl
  
=item problems
  
  * paralog gene joins are still a problem, gscore doesn't resolve all
  * ovself.eqgene  no-alt overlap col 7-1 may help, but need to collect range of alts/over-loci 
  * _G2..n multimaps should be down-weight versus 1st maps, i.e. drop candidates
  * gscore lspan>1 can reduce scores for likely joins
  
  see  cases in zfish/sra17zmerge6evg/map5best/Fixme.r3work.info
  
=item refine gscore using flags

keep flags, some are bad alt joins
maybe use lspan < -1 as bad flag
compare keep/drop pair for flags, drop one with bad flags
  69 flags=altbad --score
  99 flags=anti   --score
  11 flags=badspan  --score
1093 flags=blmap   ignore
 895 flags=inlongerr  --score? maybe not, or -small reduce
   8 flags=joinerr  --score
 387 flags=split    maybe ignore

=cut

use strict;

my $badflags='altbad|anti|badspan|joinerr'; # was also |inlongerr
my $USEBADFLAG=$ENV{badflag}||0;  
if($USEBADFLAG and $USEBADFLAG=~/[a-z]+/) { $badflags=$USEBADFLAG; }

my $MINA=$ENV{mina}||15;  #was 10 .. 20,  10% cds-overlap better? or too many spurious over?
my $GDUPWT= $ENV{gdup} || 0.50; $GDUPWT= $GDUPWT/100 if($GDUPWT>1); # reduce gscore for _G multimaps
# my $SHOW_NONOVER= $ENV{nover}||0; # default off?
my $SHOW_NONOVER= not( $ENV{skipnonover}||0 ); # default on?

my(%gscore,%gstable,%gsval,%pod,%pdv,%ovd,%pov,%ovm,%nov);
my(%chr,%cloc); # _locsort
sub _locsort { ($chr{$a} cmp $chr{$b} or $cloc{$a} <=> $cloc{$b} or $a cmp $b); }  

if(my $fs=$ENV{gscore}) { 
  open(F,$fs); while(<F>){ 
    next unless(/gscore=/); chomp; my($pd)=split; my($gs)=m/gscore=(\S+)/; $gs=1 if($gs eq "0"); 
    $gstable{$pd}=$gs; $gsval{$pd}=$_; 
  } close(F); 
}

while(<>) {
  my($pd,$oid,$ovd,$loc,$mapq)=split"\t"; # equalgene.pl output table
  my $pg=$pd; $pg=~s/_C\d+$//; $pg=~s/t\d+$//; 
  $pod{$pd}=$oid;
  my($lc,$lb)= split(/[:-]/,$loc);# =  $loc=~m/([\w\.]+):(\d+)/;
  $chr{$pd}=$lc; $cloc{$pd}=$lb; # for _locsort
  my $gs= gscore($pd); 
  unless($gs=~/\w/){ if($mapq=~m/,([\d-]+)gs/){ $gs=gscore($pd,$1); } else { $gs="no"; } }
  s/$/\t$gs.gs/; s/^$pd\s//; $pdv{$pd}=$_;
  if($ovd eq "na"){ next; } # $nov{$pd}=1; 
  my $ovm=0;
  for my $dv (split",",$ovd) { 
    my($d,$v)=split"/",$dv;  
    next if($d=~/^$pg/); # skip altover
    $v=~s/[CI]//; if($v>$MINA) { $ovd{$pd}{$d}=$v; if($v > $ovm){ $pov{$d}=$pd; $ovm{$pd}=$d; $ovm=$v; } } 
  }
}

# this now doesnt report nonoverlaps, add that?
overkeepdrop();
# END .... 

sub badflags {
  my($pd)=@_;
  return 0 unless($USEBADFLAG);
  my $bad=0;
  # my $badflags='altbad|anti|badspan|inlongerr|joinerr';
  # inlongerr maybe not bad; anti|badspan|joinerr, altbad? are bad
  if(my $gv=$gsval{$pd}) {
    my($fl)= $gv=~m/flags=(\S+)/;
    my @fl=split",",$fl;
    @fl= grep /$badflags/,@fl; $bad=@fl;
  }
  return $bad;
}

sub gscore{ 
  # ID_C1/2 split tag problems
  my($pd,$gsin)=@_; my $gs=0; 
  if($gsin) { if($pd=~m/_G\d/){ $gs= $GDUPWT*$gs; } $gscore{$pd}= $gs= $gsin; } 
  else { $gs= $gscore{$pd}; }
  return $gs if($gs); my $oid=$pod{$pd}; 
  for my $d (split",","$pd,$oid") { if($gs=$gstable{$d}) { 
    if($pd=~m/_G\d/ or $d=~m/_G\d/){ $gs= $GDUPWT*$gs; }
    $gscore{$pd}=$gs; $gsval{$pd}=$gsval{$d}; return $gs; } 
    }
  return "";
}


sub overkeepdrop { 
 my(%drop,%keep,%topov);
 # sort best score last
 # %nov nonover not in %ovm list
 my @pd=sort{ $gscore{$a} <=> $gscore{$b} or $ovm{$b} <=> $ovm{$a} or $b cmp $a } keys %ovm;
 for my $pd (@pd) { 
   my($gsomd,$tpm)= (0) x 9;
   my $gsp=gscore($pd); my $gsom=0; my $gsod=""; my @bet=(); 
   my $pbad= badflags($pd);
   my @od= grep{ not $drop{$_} } sort keys %{$ovd{$pd}};
   for my $od (@od) { my $gso= gscore($od); 
     ## gscore() refine: if(main($pd) and alt($od)) downweight gso; 
     ## check flagsof($pd) vs flagsof($od) : join,split,longerr,..
     my $obet= ($gso>$gsp)?1:0;
     
     # badflags rescue not great, gscore still valid, but need context of loci/alts overlapping
     # eg.  loca/t1,2,3 (hiscore) + locb/t1,2,3 (loscore); only locat2 (longerr/join/altbad) over locbt1,t,3 : drop locat2
     # but  loca/t1,2,3 (hiscore+longerr) + locb/t1 (losc); loct2,3 over locbt1 : drop locbt1
     my $obad= badflags($od);
     if($obad > $pbad) { $obet=0; } elsif($pbad > $obad) { $obet=1;} 
     
     if($obet) { push @bet,$od; if($gso>$gsom) { $gsom=$gso; $gsomd=$od; } } 
     my $pm= $ovd{$pd}{$od}; #o# $om=$ovd{$od}{$pd}||0; 
     if($pm > $tpm) { $topov{$pd}="$od/$pm"; $tpm=$pm; } ## topov od == ovm
     #o# if($om < $pm) { push @bet,$od; ..} 
   }
   if(@bet){ $drop{$pd}++; 
     #? maybe: $keep{$gsomd}=$pd or if(my $pov=$pov{$pd}){ $keep{$pov}=$pd; } .. pov == $gsomd ?
     my $v=$pdv{$pd}; $v=~s/$/,$gsom,$gsomd/; $pdv{$pd}=$v; }
   else { $keep{$pd}++; } # keep{$pd}= $pov{$pd}||"na";
 }

 # undrop if over ovd is drop, esp drop alts where main is keep
 # sort best score 1st
  my @dpd=sort{ $gscore{$b} <=> $gscore{$a} or $ovm{$a} <=> $ovm{$b} or $a cmp $b } keys %drop;
  for my $pd (@dpd) { my $v=$pdv{$pd}; my($od)= $v=~m/,(\w+)$/; 
    if($od and $drop{$od}) { delete $drop{$pd}; $keep{$pd}++; }
  }
 
  # change eqgene output, add top overlap locus id column; sort by location,
  my($act,$pd,$v,%act);
  $act="drop"; for $pd (sort keys %drop) { $act{$pd}=$act; }
  $act="keep"; for $pd (grep{ not $drop{$_} } sort keys %keep) { $act{$pd}=$act; } 
  if($SHOW_NONOVER){ $act="nover"; for $pd (grep{ not $act{$_} } sort keys %pdv) { $act{$pd}=$act; } }
  for $pd (sort _locsort keys %act) { print join"\t",$pd,$act{$pd},$topov{$pd}||0,$pdv{$pd}; }
}
