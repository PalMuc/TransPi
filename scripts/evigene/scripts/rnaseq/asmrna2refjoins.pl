#!/usr/bin/env perl
# asmrna2refjoins.pl
# parse joins in align.tab of asmrna2refalign.sh

use strict;
use Getopt::Long;

## classes:  join, pooj == poor-join, noj > nojo == not join
# expect table cols: 
my @tcols= qw(GeneID  gespan  geor QueryID quspan match qlen cov pid path indels cdsindel);
my @hd= @tcols;

my $MINCOV  = 10.0;   # for cov, min 10% ? or by bases?
my $HI_IDENT= 95.0; # for pid, percent identity align score
my $MINMATCH= 60;   # bases aligned
my $TRIVIAL_SAMESPAN = 0.20; # what?
my $GNUMPAT='(\d+)t(\d+)';   # bad \d match for prefix Thecc(1)EG .. need idpattern input

my (%jclass, %qval, %qpath);
    

while(<>) {
  next if(/^\W/);
  # chomp; my @v=split"\t";
  my @v=split;
  if(/^GeneID/) { @hd=@v; next; } # chomp($hd[-1]); 
  else {
    my %v=(); for my $i (0..$#v) { $v{$hd[$i]}= $v[$i]; }
    my $tid= $v{'QueryID'}; # Target alias ?
    my $gtid= $v{'GeneID'}; # Source or Reference alias ?
    # $qval{$tid} .= "$_\n"; # \@v;
    $qval{$tid}{$gtid}= \%v;  
    $qpath{$tid}++;
  }
}

foreach my $tid (sort keys %qpath) {
  my $np= $qpath{$tid};
  %jclass=();
  if($np<2) {
    my($gtid)= sort keys %{$qval{$tid}};
    # $jclass{$gtid}="nojo:onepath"; 
    print join("\t",$tid,$np,$gtid,"nojo:onepath")."\n";
    
  } else {
    my $gval= $qval{$tid};
    my %gval= %$gval;
#     my @qval= grep /\w/, split"\n", $qval{$tid};
#     foreach my $v (@qval) {
#       my @v= split"\t",$v; 
#       my %v=(); for my $i (0..$#v) { $v{$hd[$i]}= $v[$i]; }
#       my $gtid= $v{'GeneID'}; # really ref-mRNA id, tNNN alt suffix
#       $gval{$gtid}= \%v;   
#     }

    my @gd= sort keys %gval;
    for(my $i=0; $i<$#gd; $i++) {
      my $gdi= $gd[$i];  
#Doff#      next if($jclass{$gdi}); #? always? or only for i x j ?
      my $gvi= $gval{$gdi};
      my($ibe,$ior,$iqbe,$imatch,$icov,$ipid)= @{$gvi}{qw(gespan geor quspan match cov pid)};
      if(locov($gdi,$icov,$imatch)) { next; } 
      
      for(my $j=$i+1; $j<=$#gd; $j++) {
        my $gdj= $gd[$j];  
#Doff#        next if($jclass{gdj}); #? always? or only for i x j ?
        my $gvj= $gval{$gd[$j]};
        my($jbe,$jor,$jqbe,$jmatch,$jcov,$jpid)= @{$gvj}{qw(gespan geor quspan match cov pid)};
        if(locov($gdj,$jcov,$jmatch)) { next; } 
        if(samelocus($gdi,$gdj)) { next; }
        if(samespan($gdi,$iqbe,$gdj,$jqbe)) { next; }
        # is join, have 2 paths, not same tr-span, not low cover, not same locus
        my $jqual="";
        $jqual .= ($ior eq $jor or $ior eq ".")?"Fwd":"Rev";
        $jqual .= neargenes($gdi,$gdj); # FIXME: bad unless IDs are numeric location order..
        $jqual .= matchqual($gdi,$gvi,$ipid,$jpid); # hi/lo qual from ipid,jpid, indels?
## FIXME: reclass Loq and Same as not-join, as pooj == poor-join ?
	my $jtype= ($jqual =~ /Loq|Same/) ? "pooj:$jqual" : "join:$jqual";
        $jclass{$gdj} .= $jtype.".$gdi,$gdj;"; 
        # also  
        $jclass{$gdi} .= $jtype.".$gdi,$gdj;"; 
        # print join("\t",$tid,$jtype,$gdi,$gdj)."\n";
      }
    }
    
    for(my $i=0; $i<=$#gd; $i++) {
      my $gdi= $gd[$i];  
      my $jtype= $jclass{$gdi};
      print join("\t",$tid,$np,$gdi,$jtype)."\n";
    }
  }
  
}

sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }



sub neargenes { my ($dg,$eg)=  @_;
  my($gnd,$and)= $dg =~ m/$GNUMPAT/;
  my($gne,$ane)= $eg =~ m/$GNUMPAT/;
  return "Same" if($gnd==$gne);
  return "Tand" if($gnd>0 and $gne>0 and abs($gnd-$gne)==1);
  return "Near" if($gnd>0 and $gne>0 and abs($gnd-$gne)<5);
  return "Dist";
}


sub matchqual {
  my($dg,$eg,$did,$eid)= @_; 
  return "Loq" if($did<$HI_IDENT or $eid<$HI_IDENT);
  return "Hiq";
}

sub locov { 
  my ($dg,$cov,$match)= @_; 
  unless($cov >= $MINCOV or $match >= $MINMATCH) { $jclass{$dg} .= "nojo:locov.$dg;"; return 1; } # or? and?
  return 0;
}  
  
sub samelocus { my ($dg,$eg)=  @_;
  # my($gnd,$and)= $dg =~ m/$GNUMPAT/;
  # my($gne,$ane)= $eg =~ m/$GNUMPAT/;
  (my $gd=$dg) =~ s/t(\d+)$//; 
  (my $ge=$eg) =~ s/t(\d+)$//; 
  if($ge eq $gd) { $jclass{$eg} .= "nojo:altof.$dg;" unless($jclass{$eg} =~ /nojo:altof.$gd/); return 1; } # ;
  return 0;
}
  
sub samespan { 
  my($dg,$dbe,$eg,$ebe)= @_; 
  my($tb,$te)=split /[-]/,$dbe; my($lb,$le)= split /[-]/,$ebe;
  my $over= ($tb <= $le && $te >= $lb) ? 1 : 0;
  if($over) { # trival?
    my ($bb,$be)= ( _max($tb,$lb), _min($te,$le) );
    my $maxo= abs(1+$be - $bb);
    my $leno= _min( abs($le - $lb), abs($te - $tb)) || 1;
    my $pover= $maxo/$leno;
    $over = 0 if $pover < $TRIVIAL_SAMESPAN;
  }
  $jclass{$eg} .= "nojo:samespan.$dg;" if($over); 
  return $over;
}


__END__

#.... examples of multipath

# GeneID            gespan  geor  QueryID                     quspan  match   qlen    cov     pid     path    indels  cdsindel
# nojo:alt; same locus, alttr
# Thecc1EG000014t1  44-534  +     cacao3sopcsc1k25loc1957t2   1-491   491     902     54.4    100.0   1/2     0       0
# Thecc1EG000014t2  594-879 -     cacao3sopcsc1k25loc1957t2   617-902 286     902     31.7    100.0   2/2     0       0

# nojo:triv; trivial rev-td-join : reversed align for parts, but 1st is tiny (0.9% cov, 24 bp)
# Thecc1EG000026t1  4461-4484 -    cacao3sopcsc1k23loc3756t1  2596-2619   24    2619    0.9     100.0   3/3     0       1
# Thecc1EG000027t2  1-2482    +    cacao3sopcsc1k23loc3756t1  61-2542     2479  2619    94.8    99.9    1/3     0       

# nojo:multimap; same span multipath, not join
# Thecc1EG045396t1   842-1486 +    cacao3sopcsc10rk23loc108t10  2251-2895 642   2895    22.3    99.5    1/2     0       0
# Thecc1EG045396t2   477-1121 +    cacao3sopcsc10rk23loc108t10  2251-2895 642   2895    22.3    99.5    2/2     0       0

# join:revnearloq;  partial rev-near-join (90t,94t); 85t and 94t are same span multipath; lowish pctid/indels suggests 2ndary aligns or other genes?
# Thecc1EG045390t1   115-667   +   cacao3sopcsc10rk23loc108t12  1-492      460   3323    14.8    92.7    2/3     0/4     -4
# Thecc1EG045394t1   2415-3015 -   cacao3sopcsc10rk23loc108t12  2620-3250  589   3323    19.0    93.3    1/3     30/0    30
# Thecc1EG045385t1   1172-1509 -   cacao3sopcsc10rk23loc108t12  2954-3291  289   3323    10.2    85.5    3/3     0       0
