#!/usr/bin/env perl
# mcltreesplit.pl
# revise this to do better job; missing many good groupings
# also try for 2nd level, higher groups (as ARC1,2,3.. id?)
# parse () groupings, cut out subgroups mostly at 2nd-(), count n ARP ids
# drop singletons; subdivide (1?) huge TE group; remove :n and (ID1) parens

use strict;

my $MAXGG= $ENV{max} || 100; #1000; # 1 skipped at 1000; change to inf ?
my $MINGG= $ENV{min} || 0;   # change to 0 to put all
my $LEV1= $ENV{lev} || 1; # lev is bad, dropping many > lev
my $ATAG= $ENV{tag} || "ARP";  
my $UNK='^(hypothetical|unknown|uncharacterized)';
my $UNKNAME="unknown";

my($in,$inarp,$nc,$nskip,$nsing)=(0) x 10;
my($gv)=("") x 10;
my %desc; my %arp; my %sum;
my @sav= ();

while(<>){
  chomp; 
  parse($LEV1,$_);
}
parse($LEV1,")"); ## make sure end

sub parse {
  my($LEV,$line)= @_;
  
  push(@sav, $line);
  local $_= $line;
  
  # my $v = $_; $v =~ s/:\d+//; 
  my $arp= (/($ATAG[\w,]+)/) ? $1 : "";
  if($arp) {
    my $de= (s/\#\s*(.*)$//) ? $1 : ""; 
    #s/^\s+//; s/[\(\)]//g; s/:\d+//; s/\s+$//;
    #$gv .= $_;
    my @arp= grep /\w/, split ",", $arp;
    foreach my $ar (@arp) { $arp{$ar}++; $desc{$ar}= $de; }
    $gv .= "," if($gv=~/[\w\)]$/); $gv .= join",",@arp;
    $inarp=$in;
    return;
  } 

  s/^\s+//; s/\(\)//; s/:\d+//; s/\s+$//;
  return unless(/\S/);
  
  # $gv .= $_;
  if(/\(/) { $in++;  $gv .= "," if($gv=~/[\w\)]$/); $gv.="("; }
  if(/\)/) { 
    $in--; $gv.=")";   
    
    # my $dodump= ($in == $LEV or $in < 1)?1:0;
    my $dodump= ($inarp - $in > 1 or $in == $LEV or $in < 1)?1:0;
    
    if($in > $LEV ) { # 
      my $ngc= scalar(keys(%arp));
      $dodump=2 if($ngc > $MAXGG);
    }
  
    if($dodump) { # dump 
      my $ngc= scalar(keys(%arp));
      my $ndu=0; foreach my $ar (keys %arp) { $ndu++ if($arp{$_}>1); } 

use constant CUTDUP => 1;
      
      # $gv = cleangv1($gv);
      # have some dupl ARP ids from diff genes; remove? ** yes
      if(CUTDUP && $ndu>0) {
        foreach my $ac (sort keys %arp) { 
          for (my $nd= $arp{$ac}; $nd>1; $nd--) { $gv =~ s/$ac\b//; } # cuts 1st, better if last
        }
      }
      
      $gv= cleangv2($gv);
      # $gv= cleangv2($gv);

      $sum{ngc}{$ngc}++;
      $sum{ndu}{$ndu}++;
      
      if($ngc > $MAXGG and $dodump<2 ) { # only 1 of these, but has 1/4 of groups, rather divergent
        $nskip++;

# if(0) {  # is this working? no .. make it work .. need to splice out of sav/sav1 lines processed.
#         my($in0)= ($in); my @sav1=@sav;
#         $in=0; $gv="";  %desc= %arp= (); @sav=();
#         foreach my $s (@sav1) { parse($LEV+2,$s); } # hangs?
#         $in= $in0;
# }        

      } # split more ..
      elsif($ngc < $MINGG) { $nsing++; } # ignore singletons
      else { 
        $nc++; 
        my $de= desc(sort keys %arp);
        print "ARDE_C$nc\t$de\n" if $de;
        print "ARC1_C$nc\t$gv\n"; 
      } 

      $gv=""; %desc= %arp= (); @sav=();
    }
  }
}

warn "# n clust=$nc; ntoofew=$nsing; ntoomany=$nskip\n";
my @ngc= map { "$_:".$sum{ngc}{$_}; } sort{$a<=>$b} keys %{$sum{ngc}}; 
warn "# n groups=", join(",", @ngc), "\n";
my @ndu= map { "$_:".$sum{ndu}{$_}; }  sort{$a<=>$b} keys %{$sum{ndu}}; 
warn "# n dups=", join(",", @ndu), "\n";

sub cleangv1 {
  local $_= shift;
  s/,\)/)/g; s/,\s*$//;
  s/\((\w+)\)/$1/g;
  #? if(s/^\)\,/,/) { s/^\(//;  }  
  return $_;
}

sub cleangv2 {
  local $_= shift;
  s/,\)/)/g;  s/\(,/(/g; 
  s/\(\)//g;  s/,,+/,/g;
  s/\((\w+)\)/$1/g;
  
  $_= "(".$_ unless(/^\(/);  $_ .=")" unless(/\)$/);
  my $no= tr/\(/\(/; my $nc= tr/\)/\)/;
  my $d= $no - $nc;
  if($d>0) { my @ac= (")") x $d; $_.= join"",@ac; }
  elsif($d < 0) { my @ac= ("(") x (-$d);  $_= join("",@ac)  . $_;}
  return $_;
}


sub desc {
  my (@gn)= @_;
  my (%de,%degn);
  foreach my $gn (@gn) {
    my $de= $desc{$gn} or next;
    my @de=split /\s*;\s*/,$de;
    foreach (@de) {
    $_= lc($_);
    s/parent=\S+//;  s/gene=\S+//;  s/name=//;  s/\-p\w//g;
    s/\s\d\S+//g; s/\b(ens|loc)\w+//g; 
    s/conserved //ig; s/similar //g; s/probable //ig; s/putative //ig;
    s/acyrthosiphon pisum //ig; s/acyrthosiphon //i; s/nasonia vitripennis //g;
    s/phum_\w+//ig;  s/anopheles gambiae str. pest\W?//i;
    s/(uncharacterized|predicted) protein/$UNKNAME/i;
    s/\bcg\d+//g; s/\bagap\d+//g; s/\bga\d+//g;
    s/\s*;*\s*$//; s/^\W+//;  s/ \s+/ /g;
    $de{$_}++ if(/\w/); 
    }   
    #? $degn{$_}= $gn if(/\w/ and not $degn{$_});
  }

  ## fixme poor sort: add weight of ARP# (small better), downweight 'hypothetical's
  map{ $de{$_}=1 if(/$UNK/)} keys %de;
  my $ld=""; foreach my $d (sort{$a cmp $b} keys %de) { (my $ldc=$ld) =~ s/\W/./g;
    if($ld and $d =~ /$ldc/) { $de{$ld}++; } $ld=$d; }
  my @de= sort{ $de{$b} <=> $de{$a} } keys %de; # or $a cmp $b
  push @de, $UNKNAME unless @de;
  my $nd= @de;
  #NO# $nd= 2 if $nd > 2;
  my $de= join(" ;; ",@de[0..$nd]);
  $de ||= "nada";
  return $de;
}


