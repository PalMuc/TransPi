#!/usr/bin/env perl
# strandspanintab.pl

=item about

  $evigene/scripts/rnaseq/strandspanintab.pl -in introns.gff -chr chrsize.tab  > strandspan.tab

  strandspan.tab has columns of [ chr,start,end,...,fwdid,revid,.. ]; 
  cols 1-3 of partid are fed to 'samtools view -L parttab bam', to pull reads per location.
  
=item parameters

  required:
    -chrsize file : table of chromosome sizes
    -introns file : gff of intron locations (stranded)
  options:
    -PARTWIN=250000 : partition window, for part ids, bigger for less data, smaller for more.    
    -MININTRON=6    : minimum number of introns in window
    -minrevstrand=0.05 : proportion below which reverse overlapping strand in window is ignored
    -window=100  : impl detail, controls granularity of output table
    -xspan=300   : impl detail, ~exon size > window
  
=item output

  This prepares input table for evigene/scripts/rnaseq/trasmstrandspan.pl
  It replaces evigene/scripts/rnaseq/pullspanstrandsam.pl
    that uses special BigWig tables  made from strand-segregated reads 

  Spans will be no smaller than window, which occurs with bi-directional introns in window
  Otherwise spans become xspan size, with data, or bigger w/o data.
  Spans are compressed (joined) if they have no/little data (introns are the data)
  Partition IDs collect multiple spans to PARTWIN size, with separate forward (ptf001)
  and reverse (ptr001) IDs for intron direction (columns Idfwd, Idrev). Column Or indicates
  f,r,b for direction.  Partitions span small chrs/scaffolds.

=item EG

  $evigene/scripts/rnaseq/strandspanintab.pl -partwin 1000000 \
   -in intron_good.gff.gz -chr allfungr1-kfish2.chrcount > kfish2_strandspan1m.tab
  
  #Chr            Beg     End     Nfwd    Nrev    Nno     Idfwd   Idrev   Or
  Scaffold1       1       28900   0       0       0       0       0       0
  Scaffold1       28901   80200   6801    0       0       ptf001  0       f
  Scaffold1       80201   92800   0       208     0       0       ptr001  r
  Scaffold1       92801   103900  432     0       0       ptf001  0       f
  Scaffold1       103901  119900  0       180     0       0       ptr001  r
  Scaffold1       119901  122300  123     3       0       ptf001  0       f
  Scaffold1       122301  122400  41      7       0       ptf001  ptr001  b
  Scaffold1       122401  208700  340     199322  0       0       ptr001  r
  Scaffold1       208701  218400  5888    0       0       ptf001  0       f
  ...
  Scaffold173     1       130500  0       0       0       0       0       0
  Scaffold173     130501  163100  0       4216    0       0       ptr425  r
  Scaffold173     163101  169000  3640    0       0       ptf405  0       f
  Scaffold173     169001  175400  0       8652    0       0       ptr425  r
  ...
  Scaffold10009   1       1100    70      0       0       ptf409  0       f
  Scaffold6753    1       1100    1260    0       0       ptf409  0       f
  Scaffold10082   1       1100    344     0       0       ptf409  0       f
  Scaffold10179   1       1000    168     0       0       ptf409  0       f   
  # last part-ids: ptf409,ptr415
  
  $evigene/scripts/rnaseq/strandspanintab.pl -partwin 5000000
  # last part-ids: ptf082,ptr083
  
=cut

#   option -pid ptf011 means grep out this partid from table, use those locations to pull reads
#   from bam for assembly (part=11, strand=fwd == partid ptf011)

use strict;
use Getopt::Long;

my $WIND=100;
my $MINALTSTRAND=0.05; # maybe 0.01 when hi>1000; proportion below which alt-strand is ignored
# my $MINLEN=200; # what?
# my $MINREAD=500; # what?
my $XSPAN=300; # what? must be > WIND
my $MININTRON= 6; # what?

#? add partition ID: given span window, set id-column of about that size over whole genome, rev and fwd
my $PARTWIN= $ENV{partwin} || 250000;
my $compressrows= $ENV{compress} || 1; #default? NOT USED NOW
my $header= $ENV{header} || 1;
my $didhead= ($header) ? 0 : 1;
my $zerospans= 0; # this occurs only at start of scaffolds w/ 0 introns before data?

my %args=();
my $optok= &GetOptions( \%args,
  "introns=s" , "chrsize=s",
  "window=i",\$WIND,
  "xspan=i",\$XSPAN,
  "partwindow=i",\$PARTWIN,
  "zerospans!",\$zerospans,
  #unused# "compress!",\$compressrows,
  "minrevstrand=s", \$MINALTSTRAND,
  "MININTRON=s", \$MININTRON,
  # "config=s", "output=s" ,  "logfile=s" , "version=s", 
  "debug:s", "nodebug" );

my $introns= $args{introns};
my $countf=  $args{chrsize};
die "usage: $0 -introns introns.gff -chrsize genome.fa.size\nsee: perldoc $0\n" 
  unless($optok and -f $introns and -f $countf);
  
my $NOSTRANDSCAN = int(9900/$WIND);  # FIXME: options; max exon size=10kb
my $NOSTRANDX    = int(9900/$WIND);


my(%rows, %rlen,%rcount,$havecount, %partwin, %partids);

# open bamcount .. count reads/chr; skip empties.
# NOT USED now: rcount{}
# FIXME: check $rcount{chr} for small scaff readcount >> large scafs, like MT, rRNA, .., segregate w/ own part-id & mark
#foreach my $cf ($countf) {
  open(F,$countf) or die "reading -chr $countf";  
  while(<F>) { next if(/^\W/); my($r,$rlen,$rc)=split; $rlen{$r}=$rlen; $rcount{$r}+=$rc; $havecount+=$rc; } close(F);
#}

my($lr,%span);

readIntrons(); #now  puts into %row{chr}

my @chrs= sort{$rlen{$b} <=> $rlen{$a}} keys %rlen;
#^^ prefer this chr sort order for output; change intron.gff input order?

foreach my $chr (@chrs)
{
  $rows{$chr} or next;
  my @rows= @{$rows{$chr}};
  foreach my $r (@rows) { putrow2(@$r); }
}

#------------------------------------------------------

sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }

sub readIntrons
{
  # assume introns.gff sorted by location; at least chr
  if($introns =~ /\.gz$/) { open(F,"gunzip -c $introns |") or die "piping gunzip -c $introns"; } 
  else { open(F,$introns) or die "reading $introns"; }
  while(<F>) {
    next unless(/^\w/); 
    chomp; 
    ## need to handle huge intron spans?
    ## ?? handle low scores here? skip if $score < $MINSCORE ?
    my ($r,$src,$ty,$b,$e,$score,$or,$dt,$attr) = split"\t"; # gff: chr,src,intron,b,e,val,+/-,.,id..
    if($lr and $r ne $lr) { dumpspan($lr); %span=(); }
    
    my @xspan= (_max(1,$b-$XSPAN), _max(1,$b-1), $e+1, $e+$XSPAN);
    while( my($xb,$xe)=splice(@xspan,0,2) ) {
      my($ib,$ie)= (int($xb/$WIND),int($xe/$WIND));
      for( my $i=$ib; $i<=$ie; $i++) { $span{$or}[$i] += $score; } #?
    }
    $lr=$r;
  } close(F);
  
  dumpspan($lr);
}
  
sub dumpspan
{
  my($chr)= @_;
  my $chrend= $rlen{$chr}; 
  my $endw= int($chrend/$WIND);
  my ($lt,$csum,@tp);
  
  ## 1. collect strand type array
  for(my $i=0; $i<=$endw; $i++) {
    my $fv=  $span{'+'}[$i] || 0;
    my $rv=  $span{'-'}[$i] || 0;
    $fv=0 if($fv<$MININTRON);
    $rv=0 if($rv<$MININTRON);
    $csum += $fv + $rv;
    
    my $tp=0; # switch from 0, f, r, b to ints: 0, +1, -1, 2==b
    unless($fv == 0 && $rv == 0) { 
     $tp= ($fv==0) ? -1 # "r" 
        : ($rv==0) ? +1 #"f" 
        : ($fv < $MINALTSTRAND*$rv)? -1 #"r" 
        : ($rv < $MINALTSTRAND*$fv)? +1  #"f"
        : 2; #"b"; 
     }
    $tp[$i]= $tp;
    $lt=$tp;
  }
  return(0) if ($csum < $MININTRON); #  and not $zerospans  FIXME here? for no-intron part ids
  
  ## 2. fill in holes b/n stranded spans.. scan back/for for same type
  $lt=0;
  for(my $i=0; $i<=$endw; $i++) {
    my $tp= $tp[$i];
    if($tp == 0) {
      my $je=$i+$NOSTRANDSCAN; $je=$endw if($je>$endw); 
      for(my $j=$i+1; $j <= $je; $j++) {
        my $tj=$tp[$j];
        if( $tj != 0 ) { 
          if($lt != $tj and $lt != 0 and $j > $i+2) { # split diff
            $tp= $lt; my $ij= int(($i+$j)/2); 
            for(my $k=$i; $k<$ij; $k++) { $tp[$k]=$lt; } 
            for(my $k=$ij; $k<=$j; $k++) { $tp[$k]=$tj; }               
            } 
          else { 
            $tp=$tj; for(my $k=$i; $k<=$j; $k++) { $tp[$k]=$tp; } 
            }
          last; 
         } # below# elsif($j == $je and $lt != 0) {} 
        }
        
      if($tp == 0 and $lt != 0) {
           # all  tj == 0 ; fill in with last non0? 
           # try w/o Maybe dont want this.
        $je= _min($endw, $i+$NOSTRANDX); 
        $tp=$lt; for(my $k=$i; $k<=$je; $k++) { last unless($tp[$k]==0); $tp[$k]=$tp; }       
      }
    } elsif( $lt == 0 ) { # and $tp != 0
      # extend down no-strand some with tp?
        my $ib= _max(1, $i-$NOSTRANDX);
        for(my $k=$i-1; $k>=$ib; $k--) { last unless($tp[$k]==0);  $tp[$k]=$tp; }
    } else {
    }
    $lt= $tp;
  }
    
  ## 3. compress same-strand windows to largest span, and putrow()
  $lt=0;  my $b=1;
  my($lb,$le,$lfv,$lrv,$lnv)=(0) x 10; $lt=-1;
  for(my $i=0; $i<=$endw; $i++) {
    my $fv=  $span{'+'}[$i] || 0;
    my $rv=  $span{'-'}[$i] || 0;
    my $tp= $tp[$i];  my $nv=0;
    my $e=$b + $WIND-1;
    if($lt eq $tp) {
      $le= $e;  $lfv+=$fv; $lrv+=$rv; #gone# $lnv+=$nv;
      
    } else {
      # for partid change find next diff strand
      my $diffi=-1;
      if($lt == 1 or $lt == 2) { ##  =~ /^(f|b|rf)/
        for(my $j= $i; $j<$i+50; $j++) { unless($tp[$j] == 1 or $tp[$j] == 2) { $diffi=$j - $i; last; } } # =~ /^(f|b|rf)/
      } elsif($lt == -1 or $lt == 2) { ## =~ /^(r|b|fr)/
        for(my $j= $i; $j<$i+50; $j++) { unless($tp[$j] == -1 or $tp[$j] == 2) { $diffi=$j - $i; last; } }
      }
      
      putrow($chr,$lb,$le,$lfv,$lrv,$lnv,$lt, $diffi);
      ($lb,$le,$lfv,$lrv,$lnv)=($b,$e,$fv,$rv,$nv); 
    }
    $b= $e+1;
    $lt=$tp;
  }
  
  # $le= $chrend; #?? Probably want this, peg to end of chr
  # need also IDs for chr-begin to start of data.
  putrow($chr,$lb,$le,$lfv,$lrv,$lnv,$lt, 0);
}

sub strandc
{
  my($lt)=@_;
  my $ltc=($lt == 1)?'f' : ($lt == -1)?'r' : ($lt == 2)?'b' : ($lt==0) ? '0' : 'n';
  return $ltc;
}

sub putrow
{
  # my($r,$lb,$le,$lfv,$lrv,$lnv,$lt, $difftype) = @_;
  my(@v)= @_;  my $r=$v[0];
  push( @{$rows{$r}}, \@v);
}

sub putrow2
{
  my($r,$lb,$le,$lfv,$lrv,$lnv,$lt, $difftype) = @_;
  
  # FIXME: add part ids if have reads but no introns .. DONT know read counts in windows here..
  # make this option to add pids to zero-intron partitions..
  
  if($lb>0 and $le>$lb) {
    my $win= 1+$le-$lb;

    $lt=0 if(!$zerospans and $lt != 0 
      and ($lfv<$MININTRON and $lrv<$MININTRON)); #does this defeat above strand-extensions?

    ## FIXME: dont change ID till at 0,0,0 break
    my($pidfor,$pidrev)=(0,0);
    if($lt == 1 or $lt == 2 or ($zerospans and $lt == 0)) { ##  =~ /^(f|b|rf)/
      $pidfor= partid( "ptf", $win, $difftype);
    }
    if($lt == -1 or $lt == 2) { ## =~ /^(r|b|fr)/
      $pidrev= partid( "ptr", $win, $difftype);
    }
    
    my $ltc=($lt == 1)?'f' : ($lt == -1)?'r' : ($lt == 2)?'b' : ($lt==0) ? '0' : 'n';
    unless($didhead) { print "#". join("\t",qw(Chr Beg End Nfwd Nrev Nno Idfwd Idrev Or)) ."\n";  $didhead++; }
    print join("\t",$r,$lb,$le,$lfv,$lrv,$lnv,$pidfor,$pidrev,$ltc)."\n";
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
