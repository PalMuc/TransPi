#!/usr/bin/env perl
# pastebpo.pl
# gunzip -c aabugs_filtered.bpo.gz | \
# cat omclfilt1_Aug_10/aabugs_omclgn.tab dmel_prota.ids - \
# | perl pastebpo.pl > aabugs_omclgns.tab &
# but, remove dmel dupl pa: dmel_prota.ids
# are any genes in 2+ groups? no
# update 2009.12: replace dmel.ids w/ alttr.ids (drop these)

my $tag=$ENV{idprefix} || "ARP";

while(<>){

if(/^($tag\d+)\t(\S+)/) { # aabugs_omclgn.tab: ARP5688 acyr1_ncbi_hmm82293
  $gog{$2}=$1;
  
# } elsif(/^([a-z][a-z0-9]+)(\_\S+)$/) { # alttr.ids
#  $alttr{"$1$2"}++; 

#} elsif(/^(dmel_\S+)/) {  # dmel_prota.ids # drop this
#  $dmprime{$1}++;
  
} elsif(/^\d+;/) { # xxx.bpo
  my @v= split";"; 
  ($g1,$g2,$ev,$pi)= @v[1,3,5,6];
  next if($g1 eq $g2);
  #? next if($alttr{$g1} or $alttr{$g2});
  
  $og1= $gog{$g1};
  $og2= $gog{$g2};
  next unless($og1 and $og1 eq $og2);
  
  print join("\t",$g1,$g2,$og1,$ev,$pi),"\n";
  $nout++;
  }

}


# dgg: sample .bpo
# 0;  1                ; 2 ;    3              ; 4 ;  5   ; 6 ;   7
# 2;aedes_AAEL015645-PA;264;aedes_AAEL015645-PA;264;1e-145;100;1:1-250:1-250.

