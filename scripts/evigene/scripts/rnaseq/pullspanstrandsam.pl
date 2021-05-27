#!/usr/bin/env perl
# pullspanstrandsam.pl

=item about

  env header=1 compress=1 $evigene/scripts/rnaseq/pullspanstrandsam.pl \
    $pt.*.dmag3rna5n.bw $pt.*.chr.count > $pt.strandspan.tab

  strandspan.tab has columns of [ chr,start,end,...,fwdid,revid,.. ]; 
  cols 1-3 of partid are fed to 'samtools view -L parttab bam', to pull reads per location.
  option -pid ptf011 means grep out this partid from table, use those locations to pull reads
  from bam for assembly (part=11, strand=fwd == partid ptf011)
  
  Span table maker:  evigene/scripts/rnaseq/pullspanstrandsam.pl  
      this now requires input bigwig (.bw) files made from strand-segregated reads for
      created location segment tables with spans of strandedness marked.

=item revise to use introns
  
  should be able to build strand-span tables using introns.gff
  pulled from reads; simpler than this way via .bw data
  -- expand exon edge of each intron by fixed width, eg. mate distance (outer insertsize dist)
     and sum into span windows, joining short empty spans b/n two intron locs.
     
  -- options for min intron counts and for min proprotion rev/fwd
  where both strands found but minor strand is much less frequent.
  
=item input data  

  pt=grNDtCOcX
  
  -- bigwig sets of read coverage per strand (fwd, rev, noo = no orient)
  ls $pt.*.dmag3rna5n.bw
  grNDtCOcX.fwd.7.dmag3rna5n.bw   grNDtCOcX.noo.8.dmag3rna5n.bw   grNDtCOcX.rev.9.dmag3rna5n.bw

  -- tables of chr/scaffold sizes, and read counts, from samtools idxstats
  ls $pt.*.chr.count
  grNDtCOcX.fwd.chr.count   grNDtCOcX.noo.chr.count   grNDtCOcX.rev.chr.count

=cut

use strict;
my $WIND=100;
my $MINLEN=200; # what?
my $MINREAD=500; # what?

#? add partition ID: given span window, set id-column of about that size over whole genome, rev and fwd
my $PARTWIN= $ENV{partwin} || 250000;
my $compressrows= $ENV{compress} || 0; #default?
my $header= $ENV{header} || 0;
my $didhead= ($header) ? 0 : 1;

my @bwfwd= grep /\.fwd/, grep /\.bw/, @ARGV;
my @bwrev= grep /\.rev/, grep /\.bw/, @ARGV;
my @bwnoo= grep /\.noo/, grep /\.bw/, @ARGV;
my @countf= grep /\.count/, @ARGV; # several?

my(%rlen,%rcount,$havecount, %partwin, %partids);

# open bamcount .. count reads/chr; skip empties.
foreach my $cf (@countf) {
  open(F,$cf);  
  while(<F>) { my($r,$rlen,$rc)=split; $rlen{$r}=$rlen; $rcount{$r}+=$rc; $havecount+=$rc; } close(F);
}

# foreach my $bwf (@bwfwd)  
for( my $fi=0; $fi<@bwfwd; $fi++)
{
  #(my $bwr=$bwf) =~ s/fwd/rev/;  # not right; has own names
  #(my $bwn=$bwf) =~ s/fwd/noo/; # dont need?
  my $bwf= $bwfwd[$fi];  my $havebwf= ($bwf and -f $bwf);
  my $bwr= $bwrev[$fi];  my $havebwr= ($bwr and -f $bwr);
  my $bwn= $bwnoo[$fi];  my $havebwn= ($bwn and -f $bwn);
  
  unless($havecount) {
  open(F,"bigWigInfo -chroms $bwf |");  
  while(<F>) { if(/^\s+(scaf|contig)/) { my($r,$ri,$rlen)=split; $rlen{$r}=$rlen; $rcount{$r}||=$MINREAD; } } close(F);
  open(F,"bigWigInfo -chroms $bwn |"); #no need for rev?
  while(<F>) { if(/^\s+(scaf|contig)/) { my($r,$ri,$rlen)=split; $rlen{$r}=$rlen;  $rcount{$r}||=$MINREAD; } } close(F);
  }
  
  my $pid=0;  
  my @r= sort{$rlen{$b} <=> $rlen{$a}} keys %rlen;
  foreach my $r (@r) {
    my $rlen= $rlen{$r}; 
    my $rc= $rcount{$r}||0;
    my $skip= ($rc < $MINREAD or $rlen < $MINLEN)?1:0;
    my $parts= int(0.5 + $rlen/$WIND); ## my $wid= int(0.5 + $rlen/$parts); # or $WIND?
    if($rc>0) { print "#skip" if($skip); print "#size: $r\t$rlen.len\t$rc.reads\t$parts.parts\n"; }
    next if($skip); # want read count instead..
    
    # ** need to trap bloody no data noise...
    # .. check for bwf, bwr, bwn?
    my (@fvals, @rvals, @nvals);    
    if($havebwf) { @fvals= split " ", `bigWigSummary -type=max $bwf $r 1 $rlen $parts 2> /dev/null`; }
    if($havebwr) { @rvals= split " ", `bigWigSummary -type=max $bwr $r 1 $rlen $parts 2> /dev/null`; }
    if($havebwn) { @nvals= split " ", `bigWigSummary -type=max $bwn $r 1 $rlen $parts 2> /dev/null`; }
    
    ## bug: nodatainregion result from bigwig.. need bam.count ?
    ## get "no data" for fwd,rev when other, noo has reads..
    
    my $n= @fvals; $n=@rvals if(@rvals>$n); 
    next if($n==0);
    my $b=1; my $lt="n";
    my @tp; $#tp= $n;
    for my $i (0..$n-1) { 
      my $fv=$fvals[$i]= int(10*$fvals[$i]); 
      my $rv=$rvals[$i]= int(10*$rvals[$i]); 
      my $nv=$nvals[$i]= int(10*$nvals[$i]); 
      my $tp;
      if($fv == 0 && $rv == 0) { $tp="0"; } ## below now: $tp= ($nv>0) ? $lt : 0; 
      else { $tp= ($fv==0) ? "r" : ($rv==0) ? "f" : ($fv > 5*$rv)? "fr" : ($rv > 5*$fv)? "rf" : "b"; }
      $tp[$i]= $tp;
      $lt=$tp;
    }
    
    ## FIXME: type tp should split difference for noo type between last, next fwd/rev
    $lt=0;
    for my $i (0..$n-1) { 
      my $fv=$fvals[$i]; 
      my $rv=$rvals[$i]; 
      my $nv=$nvals[$i]; 
      my $tp=$tp[$i]; 
      if($nv>0 and $tp eq "0") {
        for(my $j=$i+1; $j<$n; $j++) { 
          my $tj=$tp[$j]; my $nj= $nvals[$j];
          if($nj == 0) { 
            if( $lt ne "0" and $j>$i+1) {
              $tp=$lt; for(my $k=$i; $k<$j; $k++) { $tp[$k]=$lt; } 
            }
            last;  # fime, use $lt if got here...
          } elsif($tj ne "0") { 
            if($lt ne $tj and $lt ne "0" and $j>$i+2) { 
              $tp= $lt; my $ij= int(($i+$j)/2); 
              for(my $k=$i; $k<$ij; $k++) { $tp[$k]=$lt; } 
              for(my $k=$ij; $k<=$j; $k++) { $tp[$k]=$tj; }               
            } else {
              $tp=$tj; for(my $k=$i; $k<=$j; $k++) { $tp[$k]=$tj; } 
            }
            last; 
            }
          }
        }
      $lt=$tp;
      }

          
      # compressed rows
    if($compressrows) {  
    my($lb,$le,$lfv,$lrv,$lnv)=(0) x 10; $lt=-1;
    for my $i (0..$n-1) { 
      my $fv=$fvals[$i]; 
      my $rv=$rvals[$i]; 
      my $nv=$nvals[$i]; 
      my $tp=$tp[$i]; 
      my $e=$b + $WIND-1;
      if($lt eq $tp) {
        $le= $e;  $lfv+=$fv; $lrv+=$rv; $lnv+=$nv;
      } else {
      
        # for partid change find next diff strand
        my $diffi=-1;
        if($lt =~ /^(f|b|rf)/) {
          for(my $j= $i; $j<$i+100; $j++) { unless($tp[$j] =~ /^(f|b|rf)/) { $diffi=$j - $i; last; } }
        } elsif($lt =~ /^(r|b|fr)/) {
          for(my $j= $i; $j<$i+100; $j++) { unless($tp[$j] =~ /^(r|b|fr)/) { $diffi=$j - $i; last; } }
        }
        
        putrow($r,$lb,$le,$lfv,$lrv,$lnv,$lt, $diffi);
        # print join("\t",$r,$lb,$le,$lfv,$lrv,$lnv,$lt)."\n" if($lb>0); 
        ($lb,$le,$lfv,$lrv,$lnv)=($b,$e,$fv,$rv,$nv); 
      }
      $b= $e+1;
      $lt=$tp;
    }
    
    putrow($r,$lb,$le,$lfv,$lrv,$lnv,$lt, 0);
    # print join("\t",$r,$lb,$le,$lfv,$lrv,$lnv,$lt)."\n"; 
    
    } else {
    for my $i (0..$n-1) { 
      my $fv=$fvals[$i]; 
      my $rv=$rvals[$i]; 
      my $nv=$nvals[$i]; 
      my $tp=$tp[$i]; 
      my $e=$b + $WIND-1;
      putrow($r,$b,$e,$fv,$rv,$nv,$tp, 0);
      # print join("\t",$r,$b,$e,$fv,$rv,$nv,$tp)."\n"; # or save for back-fix tp
      $b= $e+1;
      $lt=$tp;
      }
    }  
    
  }  
}

sub putrow
{
  my($r,$lb,$le,$lfv,$lrv,$lnv,$lt, $difftype) = @_;
  
  if($lb>0 and $le>$lb) {
    my $win= 1+$le-$lb;
    
    #bad# my $wor = ($lt =~ /^(f|b|rf)/) ? "ptf" : ($lt =~ /^(r|b|fr)/) ? "ptr" : "ptn";

#    ## FIXME: dont change ID till at 0,0,0 break
    my($pidfor,$pidrev)=(0,0);
    if($lt =~ /^(f|b|rf)/) {
      $pidfor= partid( "ptf", $win, $difftype);
    }
    if($lt =~ /^(r|b|fr)/) {
      $pidrev= partid( "ptr", $win, $difftype);
    }
    
    unless($didhead) { print "#". join("\t",qw(Chr Beg End Nfwd Nrev Nno Idfwd Idrev Or)) ."\n";  $didhead++; }
    print join("\t",$r,$lb,$le,$lfv,$lrv,$lnv,$pidfor,$pidrev,$lt)."\n";
  }
}

sub partid
{
  my($wor,$win,$difftype)= @_;
  $partwin{$wor} += $win;  
  my $pid= 1 + int( $partwin{$wor} / $PARTWIN);
  my $lastid= $partids{$wor}; 
  if($lastid and $pid == 1 + $lastid and $difftype > 0) { $pid= $lastid; }
  $partids{$wor}= $pid;
  return sprintf("%s%03d",$wor, $pid);
}






__END__

# bigWigInfo -chroms strandgroups/grNDtCOcX.rev.*.dmag3rna5n.bw 
# bigWigSummary -type=max strandgroups/grNDtCOcX.rev.*.dmag3rna5n.bw scaffold00008 1 9000 90

flamingo2.% bigWigInfo -chroms strandgroups/grNDtCOcX.rev.*.dmag3rna5n.bw | less
version: 4
isCompressed: yes
isSwapped: 0
primaryDataSize: 8,707,561
primaryIndexSize: 184,916
zoomLevels: 9
chromCount: 3785
        contig00031 1620 1977
        contig00119 1621 772
        ...
        scaffold00001 0 24341
        scaffold00002 1 107191
        scaffold00005 2 912674
        scaffold00006 3 5887
        scaffold00007 4 707263
        scaffold00008 5 8942
        scaffold00009 6 2342
        ...
        scaffold03402 1619 3856
basesCovered: 14,132,213
mean: 2.937357
min: 0.600000
max: 10.800000
std: 1.852405

.......
 bigWigSummary -type=max strandgroups/grNDtCOcX.rev.*.dmag3rna5n.bw scaffold00008 1 9000 90
n/a     n/a     n/a     n/a     n/a     n/a     n/a     n/a     n/a   
 n/a     n/a     1       n/a     n/a     n/a     n/a     1.6     n/a  
  n/a     1.9     1.6     n/a     n/a     2.5     4.8     5.7     6.3 
   6.6     6.4     6.4     6       6       6.3     6.3     6.2     6.1
    6       6.4     6.4     6.1     5.7     5.7     6       6.1    
5.9     5.2     6.1     6       6       5.8     5.7     5.9     5.8   
 5.5     5.7     2.9     5       5.1     4.9     2.3     n/a     n/a  
  n/a     1.3     1.3     1.3     0.6     0.6     1.3     1.3     n/a 
   1       1       1       n/a     n/a     n/a     n/a     n/a     n/a
    n/a     n/a     n/a     n/a     n/a     1       1       1      
n/a     n/a
