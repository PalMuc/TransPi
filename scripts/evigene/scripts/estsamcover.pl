#!/usr/bin/env perl
# estsamcover.pl 
#  transcripts coverage by est/rnaseq, for SAM/BAM input alignments
#  estblastcover.pl converted to SAM input
# expect alignment of reads to transcripts (not genome), and input is sorted by transcripts (= chromosomes)

use strict;

my $samflags=$ENV{samflags} || ""; 
my $OVSLOP=6; # FIXME: need overlapslop ~ < 0.01 of max part span ; or < 0.02..0.05 of min part span?
my $pctOVSLOP=$ENV{pctover} || 0.02; # need opt?
$pctOVSLOP=$pctOVSLOP/100 if($pctOVSLOP > 0.9);
my $pctCOVMIN = $ENV{pctcovmin} || 0.2; # coverage break points minimum from median cover
$pctCOVMIN=$pctCOVMIN/100 if($pctCOVMIN > 0.99);

my $genesize= $ENV{'size'} ||"";
my $geneaa = $ENV{'tr'} || $ENV{aa} ||"";
my $CDSONLY= $ENV{'cds'} || 0;
#not here# my $SWAPQT= $ENV{'swap'} || 0;

my( $bother)= @ARGV; # one or several .bam or .sam files, expect sorted by

die "usage: env size=genesize.tab estsamcover.pl rna2transcript.bam|stdin > genes.score \n"
  unless( -f $bother or ($bother =~ /^stdin|^-/i));
  
my @SCORES= (); # do below bestscore
my (  %blen, %bcdsoff, %bcdslen, @bspans, %bmated,); #  $lq, $bmax
      
# sub bint { local $_=shift; return (/e\+/) ? int($_) : $_; }

## replace geneaa/tr input w/ gene.sizetab 
# Query   Qlen    CDSlen  CDSoffs
# daphmag3tri7trimsub13loc1004c0t1        2213    1800    201-2003
# daphmag3tri7trimsub13loc1010c0t2        2216    1659    336-1997
# daphmag3tri7trimsub13loc1012c0t2        2642    1650    374-2026

my $haveqlen=0;
if($genesize) {
  open(AASIZE,$genesize) or die "FAIL: size=$genesize ...";
  while(<AASIZE>) { 
  my($id,$trsize,$cdssize,$cdsoff)=split; 
  my($b,$e)=split /[.-]+/,$cdsoff; 
  $blen{$id}=$trsize; $bcdslen{$id}=$cdssize; $bcdsoff{$id}=[$b,$e]; 
  } 
  close(AASIZE); $haveqlen=2;
} elsif($geneaa) { # warn:  not -f $geneaa
  open(AASIZE,"faCount $geneaa |") or die "FAIL: faCount  aa=$geneaa ...";
  while(<AASIZE>) { my($id,$al)=split; $blen{$id}=$al; } close(AASIZE); $haveqlen=1;
}
$CDSONLY=0 unless($haveqlen == 2);

my $inh;
if($bother =~ /^stdin|^-/i) {
  $inh= *STDIN
} elsif(not( $bother and -f $bother )) {
  die "# ERROR: Missing blastn input: $bother\n";
} elsif($bother =~ /\.bam/) {
  open(GSCORE,"samtools view $samflags $bother|")  or die "# ERROR: samtools view $samflags $bother\n"; $inh=*GSCORE;
} elsif($bother =~ /\.sam/) {
  open(GSCORE,"samtools view $samflags -S $bother|")  or die "# ERROR: samtools view $samflags -S $bother\n"; $inh=*GSCORE;
} else {
  die "# ERROR: $bother not .sam or .bam\n"; 
}

readSam(); # $inh
close($inh);



sub readSam 
{
#? Globals
# my($n_in, $n_read, $n_aln, $n_aln_mate1, $n_pair, $n_pairok, $n_pairnear, # $n_pairfar, 
#    $n_strand, $n_map0, $n_map1, $n_mapn, $n_loident,
#    $n_mate1, $n_mate2, $n_secondary, $n_duplicate,
#    $n_idbreaks, $n_intron, $len_intron,  $len_mate, $len_read, $len_chr) = (0) x 50;

my ($lchr) = (0) x 10;
while(<$inh>) { 
  next unless(/^\w/); # comments ; @SQ stuff... 

  #my @v=split;  my($q,$t,$bits,$aln,$mis,@bspan)= @v[0,1,-1,3,4, 6,7,8,9]; # 6-9 =  q. start, q. end, s. start, s. end, 
  my ($qid, $flag, $chr, $cloc, $mapq, $cigar, $matechr, $mateloc, $isize, $seq, $qual, @opt)
      = split"\t";
  # here chr == transcript ID, qid = read id    
  next unless($flag =~ /\d/ and /^\w/);    
  #x $n_in++;

  if($lchr and $chr ne $lchr) { # assuming sorted by chr/transcript id
    bestscore($lchr); # now does putout/gene
    @bspans=(); # %bmated=(); $bmax=0;
  }
  $lchr= $chr; 
  
  my $f_pair    = ($flag & 0x0001)?1:0;
  my $f_pairok  = ($flag & 0x0002)?1:0; #??
  my $f_mismatch= ($flag & 0x0004)?1:0; #  f_mismatch == $chr eq "*"
  my $f_mismate = ($flag & 0x0008)?1:0;
  my $f_isrev   = ($flag & 0x0010)?1:0; # 16
  my $f_revmate = ($flag & 0x0020)?1:0; # 32
  my $f_first   = ($flag & 0x0040)?1:0; # 64; mate ids
  my $f_second  = ($flag & 0x0080)?1:0; # 128
  my $f_secondary= ($flag & 0x100)?1:0;  # should we skip or count these?
  my $f_duplicate= ($flag & 0x400)?1:0;
  
  if($f_mismatch) { # == ($chr eq '*') == unaligned
    #x $nomap++;  
    next;
  }
  
  # Op BAM Description : SAM Cigar ops, 2010
  # M 0 alignment match (can be a sequence match or mismatch) 
  # I 1 insertion to the reference 
  # D 2 deletion from the reference 
  # N 3 skipped region from the reference 
  # S 4 soft clipping (clipped sequences present in SEQ) 
  # H 5 hard clipping (clipped sequences NOT present in SEQ) 
  # P 6 padding (silent deletion from padded reference) 
  # = 7 sequence match    << M eqiv
  # X 8 sequence mismatch << M eqiv
  # ** H is treated odd cuz seq is clipped to it
  # ** S softclip also odd, $cloc is at NON-clipped seq start.. and cend = cloc + seqlen - end S

  ## $qid.="/2" if($f_second);
  
  my $cend = $cloc; # calc correct end from align/intron spans
  my(@aln,@intr,@itype,$instrand,$inmismat);
   
  my $ciglen=0;
  if($cigar =~ m/^(\d+)M/ ) { my $m=$1; @aln=($m); $cend += $m; $ciglen+=$m; }
  elsif($cigar =~ /^(\d+)[NHSDIP]/) { @aln=(0); } # padding, need for nSnM cigar
  while($cigar =~ m/(\d+)([NHSDIP])(\d+)M/g) { 
    my($bi,$bt,$bx)=($1,$2,$3);      
    $ciglen += $bx; $ciglen += $bi if($bt =~ /[SI]/);
    $bi=0 if($bt =~ /[HSI]/); # $bt eq "H" or $bt eq "S"
    #not used here# push(@intr,$bi);  push(@itype,$bt); 
    push(@aln,$bx); 
    $cend += $bi + $bx; 
    }
  if($cigar =~ /(\d+)[SI]$/) { $ciglen+= $1; } # end, which to count??

  my $lenc= ($seq eq "*") ? $ciglen : length($seq); # can happen but most data has seq..
  # my $nx= scalar(@aln);
  my $alen=0; map{ $alen+=$_ }@aln; 
  my $amis= 0;
  foreach (@opt) { 
    if(m/NM:i:(\d+)/){ $amis += $1; } 
    elsif(m/XS:A:([\+\-])/){ $instrand=$1; } 
    elsif(m/NS:i:(\d+)/) {  $amis += $1; } #new, if XS, = Mismatches within min_anchor_len of a splice junction, > 0 is poor
    }

  # FIXME ??? f_isrev says if cend < cloc or cend > cloc; NO, isrev means revcomp read, but start == cloc
  my @bspan=($cloc, $cend);
  my $matelinked= $f_pair;
  
  #x $n_pair++  if($f_pair);  # should this count if no align?
  if($f_pair and not $f_mismate) {
    #x $n_pairok++ if($f_pairok);
    if($matechr eq "=") {
      $matelinked=2;  # should we distinguish mate linked to same chr and to other chr = transcript? vs mate-nomap
      #x $n_pairnear++;
      # FIXME??? use f_revmate flag Rev/fwd and Mate rev/fwd to pick direction to extend..
      # my $lenm= ($mateloc < $cloc) ? $cend - $mateloc : ($mateloc+$lenc) - $cloc; # outer
      # $len_mate += $lenm;       
      
      #?? my $mateend= ($f_revmate) ? $mateloc - $lenm : $mateloc + $lenm;
      # my $mb= _min($cloc, $mateloc); my $me= _max($cend, $mateend);
      # @bspan=($mb,$me);  # is this right, we are also counting mate when it comes thru
      
      #.. read + innerspan, leave out matespan, counted later
      @bspan= ($mateloc < $cloc) ? (_min($mateloc+$lenc,$cloc), $cend) : ($cloc,_max($mateloc,$cend));
    }
  }

  sumscore( $qid, $chr, $lenc, $alen, _max(0,$alen-$amis), $matelinked, @bspan);  
  
#   # blast stats
#   my $aident= _max(0,$aln-$mis); # other way to calc: $aident = $pctident * $aln;
#   sumscore( $q, $t, $bits,$aln,$aident, @bspan);   
#   $bmax= $bits if($bits > $bmax);

} 

bestscore($lchr); # now does putout/gene
}


#.....................

sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }


# sub matespan {
#   my($pairid,$eid,$b,$e)=@_;
#   $bmated{$pairid}{$eid}++;
#   my($bmin,$emax)= ($b,$e);
#   my ($mateid)= grep { $_ ne $eid } keys %{$bmated{$pairid}};
#   return ($bmin,$emax) unless($mateid);
#   #** Fixme, this reduces ident/align ratio by amount of untested insert span added here.. 
#   foreach my $sp (@bspans) {
#     my($xb,$xe,$tb,$te,$sbit,$sal,$sid,$starg)= @$sp;
#     if($starg eq $mateid) { 
#       $bmin= ($xb<$b)? $xb : $b;
#       $emax= ($xe>$e)? $xe : $e;
#       $sp->[0]= $bmin; $sp->[1]= $emax;
#     }
#   }
#   return ($bmin,$emax);
# }

# change sumscore, tESTid not needed, want union of all gene-bspans 

sub sumscore {
  ##my( $q, $t, $bits,$aln,$aident, $qb,$qe,$sb,$se) = @_;
  my( $readid, $trid, $tlen,$aln,$aident, $mated, $qb,$qe) = @_; # add mated/unmated score
  # reduce hasmate,matelinked to 1 var: mate=0,1,2
  # my $or=0;
  # if($qb > $qe) { ($qb,$qe)= ($qe,$qb); $or--; } # not here
  # if($sb > $se) { ($sb,$se)= ($se,$sb); $or--; }

  if($CDSONLY) {
    my $cdsoff= $bcdsoff{$trid} or return ;  
    my ($cdsb,$cdse)= @$cdsoff;
    return if($qe < $cdsb or $qb > $cdse);
    $qb=$cdsb if($qb < $cdsb);
    $qe=$cdse if($qe > $cdse);
  }
  
#  # add mate-pair score, esp. add span between mates as covered..
#   my $tmate=$t; 
#   if( $tmate =~ s/\.(fwd|rev)$//i or $tmate =~ s,/([12])$,,) {
#     ($qb,$qe)= matespan($tmate, $t, $qb,$qe); # $bmated{$tmate}{$t}++; 
#     }
  
  unless(@bspans) { 
    push( @bspans, [$qb,$qe,$tlen,$aln,$aident,$mated, $readid]); 
    # push( @bspans, [$qb,$qe,$sb,$se,$bits,$aln,$aident,$readid]); 
    return; }
    
  my $ov=0;
  
# need more careful cutting of EST parts that fill gaps in align.
use constant OVERCUT => 1; # test it   

  my $qlen=1+$qe-$qb; # my $slen=1+$se-$sb;
  my $qslop= _max($OVSLOP, int($pctOVSLOP*$qlen));
  # my $sslop= _max($OVSLOP, int($pctOVSLOP*$slen));
  my $cuts=0;
  # add case for this span exceeds last span?
  # ** should change this overlap filter so that part-overs are retained, part outside of last span.
  foreach my $sp (@bspans) {
    # my($xb,$xe,$tb,$te,$sbit,$sal,$sid,$starg)= @$sp;
    my($xb,$xe,$slen,$sal,$sid,$smated,$sread)= @$sp;
    unless(
      ($qe < $xb or $qb > $xe) or
      ($qe > $xe and $qb >= $xe - $qslop) or
      ($qb < $xb and $qe <= $xb + $qslop)) { 
      
      # change here, if qe > xe + qslop, or qb < xb-slop, then cut qb,qe to part outside of xb,xe..
      # need to cut sbit, sal, sid also .. messy
if(OVERCUT) {
      my ($cb,$ce)=(0,0);
      if($qe > $xe + $qslop) { ($cb,$ce)= ($xe, $qe); }
      elsif($qb < $xb - $qslop) { ($cb,$ce)= ($qb, $xb); }
      if($ce>0 and $cuts < 9) { # can we cut same align multiple times? 
        my $pc= ($ce-$cb)/($qe-$qb);
        $tlen=int($pc*$tlen); $aln=int($pc*$aln); $aident=int($pc*$aident); 
        ($qb,$qe)= ($cb,$ce); $cuts++;
      } else { $ov=1; last; }
} else {
      $ov=1; last; 
}      
      }

#     if( $starg eq $t) {
#     unless(  # target spans, skip this for est unless same target-id
#       ($se < $tb or $sb > $te) or
#       ($se > $te and $sb >= $te - $sslop) or
#       ($sb < $tb and $se <= $tb + $sslop) ) { $ov=1; last; }
#     }

  }
   

  unless($ov) { push( @bspans, [$qb,$qe,$tlen,$aln,$aident,$mated,$readid]); }
  # unless($ov) { push( @bspans, [$qb,$qe,$sb,$se,$bits,$aln,$aident,$t]); }
}

use constant SHOW_READID => 0;

sub bestscore {
  my($gid)= @_;
  my ($tbit,$taln,$tident,$t,$nid)= (0) x 9;
  my($tmated,$tpaired,$unpaired)=(0,0,0);
  my $eids="";
  
  @bspans= sort{ $a->[0] <=> $b->[0] } @bspans;
  my($gb,$ge)= ($bspans[0]->[0], $bspans[-1]->[1]);
  my $galignspan=1+$ge-$gb;  
  my $gsize= $blen{$gid} || $ge;  # ge = best guess...
  
  if($CDSONLY) { $gsize= $bcdslen{$gid} || $ge;  }
  my $cdsoff= $bcdsoff{$gid} || [$gb,$ge];  # ge = best guess...
  my ($cdsb,$cdse)= @$cdsoff;
  
  my @cov; $#cov= $ge; ## = (0) x $ge;
  my @cdscov; $#cdscov= $ge; ## = (0) x $ge;
  # .. change cdscov to include any hit to cds, with extensions past..
  # .. extend it also with utr-covers that also hit existing cdscov, so that
  # .. cov - cdscov becomes possible spurious extensions, report that, as tcov - ccov
  my @eids=();
  foreach my $sp (@bspans) {
    #my($xb,$xe,$tb,$te,$xbit,$aln,$aident,$targid)= @$sp; 
    my($xb,$xe,$slen,$aln,$aident,$mated,$targid)= @$sp;
    $tbit += $slen; $taln+= $aln; $tident+= $aident;
    # mated==0, no mate; ==1 unlinked mate, ==2 linked mate
    if($mated>0) { $tmated++; $tpaired++ if($mated==2); }
    
    if(SHOW_READID) { push @eids, $targid if(@eids < 3); }  # ** Not useful for rnaseq..
    
    my $xmid= int( ($xb+$xe)/2);
    for my $i ($xb..$xe) { $cov[$i]++; } 
    if( ($xb >= $cdsb and $xb < $cdse) or ($xe > $cdsb and $xe <= $cdse)
      or ($cdscov[ $xb+1 ] or $cdscov [ $xe-1 ]) # extend if already overlapping this xbe
      ) {
      for my $i ($xb..$xe) { $cdscov[$i]++; } 
      }
    }
  $eids= join(",",sort @eids) if(SHOW_READID);
    
    # revise cov stats to use cov.count << median count as gap; ie big dip in coverage = gap
    # $pctCOVMIN=0.2; @ccov= sort{$b<=>$a} @cov; $covmd=$ccov[int($ncov/2)]; $covmin=_max(1,int($covmd*$COVMINPCT)); 
    # if($cov[$i] < $mincov) ..
    
  my @covsort= sort { $b<=>$a } grep{ $_ > 0 } @cov;
  my $covmed= $covsort[ int(@covsort/2) ];
  my $covmin= _max(1, int($covmed * $pctCOVMIN));
  
  my $gaps=""; 
  my ($tcov,$tgap,$gapb,$gape,$gapmax, $ccov)= (0) x 10;
  for(my $i=$gb; $i<=$ge; $i++) { 
    if($cdscov[$i]) { $ccov++; } 
    if($cov[$i] >= $covmin) { $tcov++;   # if($cov[$i] >= $mincov) ..
      if($gape) { $gaps.="$gapb-$gape,"; my $gw= 1+$gape-$gapb; $tgap += $gw; $gapmax=$gw if($gw>$gapmax); } 
      $gapb=$gape=0; }
    else { if($gapb) { $gape=$i; } else { $gapb=$i; } }
  } 
  if($gape) { $gaps.="$gapb-$gape,"; my $gw= 1+$gape-$gapb; $tgap += $gw; $gapmax=$gw if($gw>$gapmax); } 
  $gaps=0 unless($gaps);
 
#   if(scalar(%bmated)) {
#     my @mates= sort keys %bmated; 
#     foreach my $m (@mates) { $tmated++; my @two= sort keys %{$bmated{$m}}; 
#       if(@two>1) { $tpaired++; } else { $unpaired++; }
#     }
#   }

  my $mates= ($tmated==0)? 0 : "$tpaired/$tmated";

  # dont need both taln, tcov .. should be same, near same but for slop above
  # FIXME added mate insert span, tcov larger than taln for mated spans; revert to put both taln,tcov
  #NO# $taln= $tcov; # best one for now .. 
  
  # FIXME: add top 1,2 EST ids for matching transcripts by best EST id. and/or use best per span
  # FOXME: change tr=tr.fa input to table of lengths w/ CDS len and span, UTR len and span
  # ...    modify above sums to provide CDS-only and CDS-UTR cover to hunt for stat of gene joins.
  
  ## Scores from trasm review paper:
  # Accuracy = Ident / Total-align (should include gaps)  
  # my $accuracy = int(0.5+100 * $tident / $galignspan);
  # Completeness = alignment/total-size > criterion (0.80) 
  # my $complete = int(0.5+100 * $taln/$gsize); ## > 80% / Total-align (should include gaps)  
  
  # Contiguity = ESTs fully (80%?) covered by transcript; not here..
  # Chimeric/joingene : use maxgap > xxx ? 
  
   # add $tcov/$blen == pct-align,  $tident/$taln = pct-ident ?
  my $cdsUncov= _max(0, $tcov - $ccov); # == UTRgap? only for not CDSONLY
  my $scores= join("\t",$gsize,$tbit,$tident,$taln,$tcov,"$gb-$ge",$mates,$cdsUncov,$tgap,$gapmax,$gaps, $eids); 
  unless(@SCORES) { # header
  if($CDSONLY) { 
  @SCORES= qw(CDSlen TLen Iden Algn Ncov Covs Mated UTRgap Ngap Maxgap Gaps ESTids);
  } else {
  @SCORES= qw(Trlen TLen Iden Algn Ncov Covs Mated UTRgap Ngap Maxgap Gaps ESTids);
  }
  print join("\t","Query",@SCORES)."\n";
  }
  print join("\t",$gid,$scores)."\n";
  return($scores);
}
  

__END__

cacao reads x trasm test

/bio/bio-grid/cacao3/rnas/asmcov/bamtr
env size=cacao3tri4asm_cd.sizetab cds=0 $evigene/scripts/estsamcover.pl tc07l1_CAGATC_1-cacao3tri4asm_goods.bam > tc07l1_cacao3tri4asm_goods.covtab
env size=cacao3vel47all_cd.sizetab cds=0 $evigene/scripts/estsamcover.pl tc07l1_CAGATC_1-cacao3vel47all_goods.bam > tc07l1_cacao3vel47all_goods.covtab


foreach ctab (*.covtab)
foreach? echo $ctab
foreach? cat $ctab | perl -ne '($gd,$gl,$b,$id,$aln,$ncov,$cov,$mt,$ugap,$ngap,$mgap,$gaps,$eids)=split; $eids{$
eids}++; \
$mt="0/1" unless($mt); ($mh,$mt)=split"/",$mt; $mu=$mt-$mh; $nmu++ if($mu>0); $smh+=$mh;  $sm+=$mt; \
$nmh++ if($mh>0); $ng++; $slen+=$gl;  $scov+=$ncov; $saln+=$aln; $sid+=$id; $sugap+=$ugap; $sgap+=$mgap;\
END{ @eid2=grep{$eids{$_}>1} keys %eids; $neid2=@eid2; $ac=int(100* $sid/$saln); \
$comp=int(100*$scov/$slen); $ugap=int(100*$sugap/$scov); $gap=int(100*$sgap/$scov);\
$mate=int(100*$nmh/$ng); $mateu=int(100*$nmu/$ng);  $mateh=int(100*$smh/$sm);  \
print "ng=$ng; acc=$ac; compl=$comp; UTRunlink=$ugap; gap=$gap; mate=$mate,$mateu,$mateh;". \
" sharedest=$neid2;  sid=$sid; saln=$saln; scov=$scov; smated=$smh; sugap=$sugap; slen=$slen\n"; }'


tc07l1_cacao3cuf13aset_goods.covtab
ng=15401; acc=99; compl=67; UTRunlink=34; gap=12; mate=97,63,95; sharedest=1;  sid=12626439; saln=12682220; scov=23493519; smated=535660; sugap=8157706; slen=34840927

tc07l1_cacao3tri4asm_goods.covtab
ng=19645; acc=99; compl=73; UTRunlink=26; gap=7; mate=95,48,97; sharedest=1;  sid=15026122; saln=15090561; scov=27000377; smated=635061; sugap=7178910; slen=36735098

tc07l1_cacao3vel47all_goods.covtab
ng=33873; acc=99; compl=65; UTRunlink=31; gap=12; mate=95,59,94; sharedest=1;  sid=19946988; saln=20051199; scov=36818934; smated=742924; sugap=11463807; slen=56418889



