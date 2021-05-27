#!/usr/bin/env perl
# locustab_basic.pl was gmaploci_basic.pl
# merge in locusgmapsub.sh here as perl..

$MINAA=$ENV{minaa}||99;
$MINCV=$ENV{mincov}||85; # min coverage on genome, as cov * pident
# $MINOV=$ENV{minov}||20; # exon overlap, treat cds/exon ov score same here.
$MINOVC=$ENV{minov}||10; # exon overlap, treat cds/exon ov score same here.
$MINOVX=$ENV{minovx}||25; # exon overlap, treat cds/exon ov score same here.
# ^^ split out cds-ov >=15 and exon-ov >= 33
$MINALN=$ENV{minaln}||0;
$OKALN=($MINALN>40)?1:0; # special filter, only good mapping trs

#..................
# $USEGID=1; # $ENV{usegid}  # parse ID = GeneID + tAltnum
## strand problem w/ g-or x -sense; flip g-or for -sense?
## SENSEFIX no help here, need in gff before eqgene
use constant SENSEFIX => 0;
use constant USEGID => 1;

@alns=grep /align.tab/, @ARGV;
#notnow# unless(@alns) { opendir(D,"./"); @alns= grep /align.tab/, readdir(D); closedir(D); }
die "usage: gmaploci.pl xxx.align.tab + xxx.eqgene to xxx.loci" unless(@alns);

foreach my $aln (@alns) {

$pt=`basename $aln .align.tab`; chomp($pt);
$outf="$pt.a${MINAA}c${MINCV}l${MINALN}.loci";
$nin=$nloc=$skipcov=$skipaa=0;

# aln.tab
# zGenomeID     gespan  geor    
# 3:AQueryID    quspan  match   
# 6:qlen        cov     pid     path    indels  
# 11:nexon      splice  aalen offs      
# 15:aamap      sense   oid     tag

## DAMN align.tab got oid not pubid from Target.gff
open(F,"$pt.align.tab"); 
while(<F>) { next unless(/^\w/); @v=split"\t"; 
  ($gr,$gbe,$go,$id,$matc,$ql,$cov,$pid,$path,$nx,$aq,$am,$sens)=@v[0,1,2,3,5,6,7,8,9,11,13,15,16];
  ($aw)=$aq=~m/(\d+)/; $aw||=$am;  $nin++;
  if(SENSEFIX) { if($sens<0 and $go ne ".") { $go=($go eq "-")?"+":"-"; } }
   ## maybe add other filters here: minmatch = 300nt; geor ne "."; nexon>1 ?
  if($OKALN and ($go eq "." or $matc < $MINALN or $nx < 2)) { $skipcov++; $taa{$id}=0; next; }
  unless($OKALN) { $copi=$cov*$pid/100; if($copi < $MINCV) { $skipcov++; $taa{$id}=0; next; } }
  $tloc{$id}="$gr:$gbe:$go"; $taa{$id}=$aw; $tfo{$id}="$ql,$cov,$pid,$nx,$path";
  push @alids, $id;
} close(F);

open(F,"$pt.eqgene"); 
while(<F>) { 
  next unless(/^\w/); ($id,$oid,$eq,$loc)=split"\t"; 
## Fix align oid>pubid here? 
# if($tloc{$oid} && ! $tloc{$id}){ $taa{$id}=$taa{$oid}; $tloc{$id}=$tloc{$oid}; $tfo{$id}=$tloc{$oid};}

  next unless($tloc{$id});
  if($taa{$id} < $MINAA) { $skipaa++; next; }
  @eq= grep{ $taa{$_} >= $MINAA } 
    map{ ($ed,$ec,$ex)=split /[\/\.]/; $ec=~s/[IC]//; ($ec >=$MINOVC or $ex >=$MINOVX)?$ed:""; } split",",$eq;

## FIXME: missing loc links, via eq, get same ID in 2+locs..
## .. find smallest lid > 0 in all of id,eq .. reset all id,eq to that lid
if(1) {
  $lid=0; 
  foreach $d ($id,@eq) { if($ll=$idloc{$d}) { $lid=$ll if($lid==0 or $ll<$lid); } }
  $lid= ++$nloc unless($lid);
} else { # old
  $lid=$idloc{$id}||0;
  unless($lid) { 
   foreach $d (@eq) { if($lid=$idloc{$d}) { last; } } 
   $lid= ++$nloc unless($lid);
   }
}
  foreach $d ($id,@eq) { $locs{$lid}{$d}++;  $idloc{$d}=$lid; }
} close(F);

%ts=%tb=%te=();
foreach $d (@alids) {
   my($ts,$tb,$te,$to)=split/[:-]/,$tloc{$d};
   if($te){ $tb{$d}=$tb; $te{$d}=$te; $ts{$d}=$ts; }
}
@alids= sort{ $ts{$b} cmp $ts{$a} or $tb{$a} <=> $tb{$b} or $te{$b} <=> $te{$a} } @alids;
$nids=@alids;

## add: count same-trgene in altids, as "ntr/nsame";
##  count also trgene-alts at diff loci?
my %gloc;
if(USEGID) {
foreach my $t (keys %idloc) { (my $g=$t) =~ s/t\d+$//; $gloc{$g}{$idloc{$t}}++; }

if(1) {
## try to correct eqgene locs mistakes? got lots of same gid, overlap span, diff locid.
## looks too messy.  "lots" = ~ 4000 same gid adjacent loci; fail1: bug @sp
## 1st try after bugfix: maybe works, need also delete $gloc{$gdl}{rel}

 $nreloc=$nrelocd=0;
 for my $g (sort keys %gloc) { 
   my @idl= sort { $gloc{$g}{$b} <=> $gloc{$g}{$a} or $a <=> $b }  keys %{$gloc{$g}}; 
   if(@idl>1) { @spa=(); @reloc=(); @relocd=(); $l1=0;
    for my $l (@idl) { 
       ($d1,@d)= grep{ $locs{$l}{$_} } grep(/$g/, @alids); next unless ($d1);
       @sp=($ts{$d1},$tb{$d1},$te{$d1}); 
       if(@spa) { $ov=($spa[0] eq $sp[0] and $spa[1] < $sp[2] and $spa[2] > $sp[1])?1:0; 
         if($ov) { push( @reloc,$l); push(@relocd, $d1, @d);  }
         } 
       else { @spa=@sp; $l1=$l;  }
    }
  if(@reloc) { 
    #o for $rel (@reloc) { delete $locs{$rel}; delete $gloc{$g}{$rel}; $nreloc++; }
    #o for $d (@relocd) { $idloc{$d}= $l1; $locs{$l1}{$d}++; $nrelocd++; } 
    for $rel (@reloc) { 
      for $d (@relocd) { delete $locs{$rel}{$d}; } # same?: delete @{$locs{$rel}}{@relocd};
      delete $gloc{$g}{$rel}; $nreloc++; }
    for $d (@relocd) { $idloc{$d}= $l1; $locs{$l1}{$d}++; $gloc{$g}{$l1}++; $nrelocd++; } 

    }
  }
 }
 $nloc= scalar(keys %locs);
}

}

open(OUT,">$outf");
print OUT "#loci from eqgene, align.tab, nin=$nin, nloc=$nloc,reloc=$nreloc, ntr=$nids, "
  ."MINAA=$MINAA, MINCOV=$MINCV, MINALN=$MINALN, MINOVERCDS=$MINOVC\n";
print OUT "#".join("\t",qw(locid location idmain aasize map_info ntr altids ))."\n";

foreach $aid (@alids) 
{
 $lid=$idloc{$aid} or next; # what of case of no eqgene idloc??
 next if($didloc{$lid}++);
 @ids= sort keys %{$locs{$lid}}; $ntr=@ids;
 %ts=%tb=%te=(); foreach $d (@ids) {
   my($ts,$tb,$te,$to)=split/[:-]/,$tloc{$d};
   if($te){ $tb{$d}=$tb; $te{$d}=$te; $ts{$d}=$ts; }
 } 
 my %gid=();  my %olocs=();
 if(USEGID) {
   for my $t (@ids) { (my $g=$t)=~s/t\d+$//; $gid{$g}++; map{ $olocs{$_}++ } keys %{$gloc{$g}}; }
   for my $t (@ids) { (my $g=$t)=~s/t\d+$//; $gid{$t}=$gid{$g}; }
 } else {
   for my $t (@ids) { $gid{$t}= $taa{$t}; } # for sort == aasize
 }
 ($idl,@xids)= sort{ $gid{$b}<=>$gid{$a} or $taa{$b}<=>$taa{$a} 
   or $ts{$b} cmp $ts{$a} or $tb{$a} <=> $tb{$b} or $te{$b} <=> $te{$a} } @ids;
 my $nsameg=$gid{$idl};  (my $gdl=$idl)=~s/t\d+$//;
 my $noloc=scalar(keys %{$gloc{$gdl}}); # scalar(keys %olocs); # which?
 $ntrv= (USEGID) ? "$ntr/$nsameg/$noloc" : $ntr;
 $tloc= $tloc{$idl}; $taa=$taa{$idl}; $tfo=$tfo{$idl};
 $tloc=~s/:/\t/; $xids=join",",@xids;
 print OUT join("\t",$lid,$tloc,$idl,$taa,$tfo,$ntrv,$xids)."\n";
} 
print OUT "#sum: ninalign=$nin; nloc=$nloc; nskipcov=$skipcov; nskipaa=$skipaa\n";
close(OUT);

} # aln


=item CLUSTER script gmap2lociset.sh

  #! /bin/bash
  ## env ingmap=xxx.gmap.out datad=`pwd` qsub -q normal gmap2lociset.sh
  ##...... CLUSTER script for calling this w/ partitioned genome data inputs .........
  #PBS -N gmap2lociset
  #PBS -l nodes=1:ppn=32,walltime=9:55:00
  #PBS -V
  ## version for trsplit, n=30 parts, input gmap.out, add gmap.gff step
  
  ncpu=32
  ## make these, other params part of env call plus defaults
  export minaa=89
  export mincov=10
  export evigene=$HOME/bio/evigene/
  
  if [ "X" = "X$datad" ]; then echo "ERR: missing datad=pathtodata"; exit -1; fi
  cd $datad/
  ## input either .gmap.out, gmap.gff .. skip .out if have .gff
  if [ "X" = "X$ingmap" ]; then ingmap=`ls *.gmap.gff`; fi
  if [ "X" = "X$ingmap" ]; then ingmap=`ls *.gmap.out*`; fi
  
  i=0; for igm in $ingmap; do { 
  
    $evigene/scripts/rnaseq/locusgmapsub.sh  $igm &
    
    i=$(( $i + 1 ));
    if [ $i -ge $ncpu ]; then wait; i=0; fi
  } done
  wait
  

=cut 