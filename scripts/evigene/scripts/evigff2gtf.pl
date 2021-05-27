#!/usr/bin/env perl
# evigff2gtf.pl

my $keep=1;
while(<>) {
  if(/^\W/) { 
    if(/^##gff-version/){ print "##gff-version 2.5 (gtf)\n"; } else { print; } 
  } else { 
    my ($c,$s,$t,$b,$e,$p,$o,$xp,$at)=split"\t";
    next unless($t =~ /^(mRNA|exon|CDS)/); # no CDS ?
    if($p =~ /,/) { my @p=split",",$p; $p=0; map{ $p=$_ if($_ > $p); } @p; }#? $p=$1 if($at=~/scoresum=(\d+)/); 
    if($t =~ /mRNA/) {  $keep=1; } ## ($p < 0.05 or /overpoor/)? 0: 1;  #? any filter
    next unless($keep);
    $at =~ m/(ID|Parent)=([^;\s]+)/; 
    my $tid=$2; (my $gid = $tid) =~ s/t\d+$//; # tNNN is transcript num at end of ID
    $t=~s/mRNA/transcript/;
    print join("\t",$c,$s,$t,$b,$e,$p,$o,$xp,'gene_id "'.$gid.'"; transcript_id "'.$tid.'"')."\n";
  }  
}


