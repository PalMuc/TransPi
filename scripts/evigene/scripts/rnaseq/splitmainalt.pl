#!/usr/bin/env perl
# splitmainalt.pl
# env idtab=kf2pub11_both.renewids nam=kf2pub11_both2a perl splitmainalt.pl  kf2pub11_both.an2.gff

## replace annot altc=oldclass with newclass from IDtab

#BEGIN
$nam=$ENV{nam}||"noname"; 
$idtab=$ENV{idtab} || "noidtab";
open(MG,">$nam.main.gff"); open(AG,">$nam.alt.gff"); open(EQ,">$nam.equalalt.tab"); 

open(ID,$idtab) or die "reading idtab=$idtab\n"; 
while(<ID>) { 
  next unless(/^\w/); chomp; 
  ($td,$oid,$gd,$ti,$cla,$aaq,$pia,$nots)=split"\t"; 
  if($nots=~/oldid:(\w+)/) { $old=$1; $onewid{$old}=$td; } 
  if($nots=~m/oidold=(\w+)/) { $oidnew{$oid}=$td; } # funky record replacement, patch gff IDs below  
  $aclass{$td}=$cla; $info{$td}=$_;
  $o2pub{$oid}=$td; $oids{$td}=$oid; 
}

## DANG, some missing mRNA for idtab entries .. count/show those ..
my($nerr,$nmain,$nalt,$nmiss,$ti)=(0) x 10;
while(<>) {
	if(/^\W/) { next; }
	elsif(/\tmRNA/) { 
		($id)=m/ID=([^;\s]+)/; ($oid)=m/oid=(\w+)/; # can miss oid on Split=2
		$oldid=$newid=""; 
		if($newid= $onewid{$id}) { 
			$oldid=$id; $id=$newid; s/ID=$oldid/ID=$newid/; 
		} elsif($newid= $oidnew{$oid}) {
			$oldid=$id; $id=$newid; s/ID=$oldid/ID=$newid/; 
	  }
		
		if($aclass= $aclass{$id}) { s/altc=\w+/altc=$aclass/; }
		$taboid= $oids{$id};
		if($oid and $oid ne $taboid) { s/$/;ERRoidmiss=$taboid-$oid/; $nerr++; }
		s/aalen=(\d+);/aamap=$1;/;
		s/trg=\w+/trg=$id/;
		$gotid{$id}++; # some are Split
		# problem now w/ 2 input.gff, if already gotid and not Split ..
		
		($ti)= $id =~ m/t(\d+)$/; 
		$ti=0 if($gotid{$id}>1 and not m/Split=/);
		
		if($ti == 1) { 
		  $mid=$id; $moid=$oid||"na"; $mloc=join":",(split)[0,3,4,6]; $minfo{$mid}="$moid\t$mloc"; 
		  $nmain++; }
		elsif($ti>1) { 
		  $mid=$id; $mid=~s/t\d+$/t1/; push @{$maid{$mid}}, $id; 
		  $nalt++; } 
	} else { # only CDS|exon ??
		my ($pid)= m/Parent=(\w+)/;
		if($oldid and $pid eq $oldid) { s/Parent=$pid/Parent=$newid/; } # ERR unless pid eq oldid
		s/trg=\w+/trg=$id/;
	}
	if($ti == 1) { print MG $_; } elsif($ti>1) { print AG $_; }
}

# END

foreach my $td (sort keys %info) {
  next if($gotid{$td});
  $nmiss++; $tinfo=$info{$td};
  ($ti)= $td =~ m/t(\d+)$/; 
  if($ti == 1) { print MG "#MISS:\t$tinfo\n"; } elsif($ti>1) { print AG "#MISS:\t$tinfo\n"; } 
}

my $stat="# splitmainalt $nam: nmain=$nmain, nalt=$nalt, nmiss=$nmiss, nerr=$nerr;\n";
print MG $stat; warn $stat;
close(MG); close(AG);
putalleq(); 

sub putalleq { 
	foreach my $mid (sort keys %minfo) { 
		my @aid=@{$maid{$mid}}; 
		my $aid= (@aid) ? join ",", map{ "$_/C90.90" } @aid : "na";  
		my($moid,$mloc)= split"\t",$minfo{$mid}; 
		print EQ join("\t",$mid,$moid,$aid,$mloc)."\n"; 
		} 
}  

__END__

=item work
#-------------------------------------------
#/bio/bio-grid/kfish2/rnas/kf2evgr/trevg367mixx/alt11
# env idtab=kf2pub11_both.renewid6 nam=kf2pub11_both2a ../scripts/splitmainalt.pl  kf2pub11_both.an2.gff

## run2
# splitmainalt kf2pub11_both2a: nmain=101221, nalt=239572, nmiss=2119, nerr=;

## these are not in input gff from gmapz/kf2pub11.*.gff
## OIDs not in ../gmapz/kf2pub11.*3m.idtab but are in ../gmapz/kf2pub11.idtabxy3m
# grep -c "#MISS" kf2pub11_both2a.*.gff
# kf2pub11_both2a.alt.gff:1508  trclass1.drop=1370; utrorf=138
# kf2pub11_both2a.main.gff:611  trclass1.drop=374;  utrorf=237
#   alt2 of main misses: n=482 of 611; n=472 exist in alt.gff
#    drops are (some) Fungr replacements, not gmapped in pubset1
#  ? replace these misses, or let alts cover? for utrorf letting alts cover is probably best.

#... replace main.miss w/ alts
gzcat ../gmapz/kfish2evg367mixx-kf2.gmap.gff.gz | cat missmain2a.oids - | perl -ne\
'if(/^(oid=\w+)/) { $ok{$1}=1 } else { if(/\tmRNA/) {
($oid)=m/(oid=\w+)/; ($id)=m/ID=(\w+)/;  $ok=$ok{$oid}; if($ok and $id
=~ /_C[12]$/) { ($idc=$id)=~s/_C[12]//; } elsif($id =~ /^$idc.C[12]/ ) {
$ok=1; } } print if($ok); }' > missmain2a.gff

cat missmain2a.trclass1 kf2pub11_both.renewid6 | perl -ne\
'if(/\tdrop/) { ($od,$drp,$cla,$ad)=split; $oad{$od}=$ad; } else {
($pd,$od,$gd)=split; if($ad=$oad{$od}) { s/[\.]$//; s/$/oidold=$od,/;
s/\t$od/\t$ad/; } print; }' > kf2pub11_both.renewid7

# dang need to put renewid pubid into missmain2a.gff .. via splitmainalt.pl ??

cat kf2pub11_both.an2.gff missmain2a.gff |\
env idtab=kf2pub11_both.renewid7 nam=kf2pub11_both2b ../scripts/splitmainalt.pl 

# splitmainalt kf2pub11_both2c: nmain=101630, nalt=239565, nmiss=1766, nerr=0;
# splitmainalt kf2pub11_both2b: nmain=101638, nalt=239565, nmiss=1766, nerr=0;
grep -c '#MISS' kf2pub11_both2b.*.gff
kf2pub11_both2b.alt.gff:1516
kf2pub11_both2b.main.gff:250 ## got the 372 replacement missmain2a.gff ok
 .. but 9 dupl mRNA records, one each in main,alt.gff, from shuffling missmain set.

nam=kf2pub11_both2c
$evigene/scripts/altbest.pl $altopts -main $nam.main.gff -alt $nam.alt.gff -eqtab $nam.equalalt.tab >& log.ab7

# sum.intrchains classified oknewlocus=3391, okaltend=10362, okay=80718, 
#   dropdup=63412, dropshort=24051, dropsubchain=20133 for 27314 mains, 174753 alts as to kf2pub11_both2c.altbest.tab
 
#............

head main.misses
#MISS:  Funhe2Exx11m000317t1    Funhe2Emap3m000029t1utrorf      Funhe2Exx11m000317      1       main    2376,36%,complete-utrbad        99/100/.        aaref:547,UniRef50_D6W540,refok,
    alt2: 1803,91%,complete,72X   aaref:1175,UniRef50_Q6JAN0,
#MISS:  Funhe2Exx11m000332t1    Funhe2Eq7m070760t1utrorf        Funhe2Exx11m000332      1       maina2  2345,37%,complete-utrbad        99/100/.        aaref:527,UniRef50_H7BYX6,
    alt2: 1225,83%,complete,161X  aaref:527,UniRef50_H7BYX6,
#MISS:  Funhe2Exx11m000488t1    Funhe2E6bm000458t1      Funhe2Exx11m000488      1       main    2053,48%,complete-utrbad        99/100/.        aaref:1528,UniRef50_Q9Y3S1,
    alt2: 1987,83%,complete       aaref:1510,UniRef50_Q9Y3S1,
#MISS:  Funhe2Exx11m000564t1    Funhe2E6bm000532t1      Funhe2Exx11m000564      1       main    1958,80%,complete       99/100/.        .
#MISS:  Funhe2Exx11m000573t1    Funhe2E6bm000540t1      Funhe2Exx11m000573      1       main    1950,73%,complete       99/100/.        .
#MISS:  Funhe2Exx11m000617t1    Fungr1EG3m000201t3      Funhe2Exx11m000617      1       main    1906,90%,complete       99/100/.        .
#MISS:  Funhe2Exx11m000699t1    Funhe2E6bm000676t2      Funhe2Exx11m000699      1       main    1835,78%,complete       99/100/.        .
#MISS:  Funhe2Exx11m000706t1    Funhe2E6bm000684t1      Funhe2Exx11m000706      1       maina2  1827,73%,complete       99/100/.        .
#MISS:  Funhe2Exx11m000712t1    Funhe2E6bm000689t1      Funhe2Exx11m000712      1       main    1823,86%,complete       99/99/. .
#MISS:  Funhe2Exx11m000729t1    Funhe2Eq7m039160t1      Funhe2Exx11m000729      1       maina2  1807,83%,complete       99/100/.        aaref:1852,UniRef50_M3ZVQ7,

alt2 of main.misses
Funhe2Exx11m000317t2    Funhe2E6bm000633t1      Funhe2Exx11m000317      2       althi   1803,91%,complete,72X   99/94/./Funhe2E6bm000633t3      aaref:1175,UniRef50_Q6JAN0,
Funhe2Exx11m000332t2    Funhe2E6bm000305t2      Funhe2Exx11m000332      2       althi1  1225,83%,complete,161X  99/100/./Funhe2Eq7m070760t1utrorf       aaref:527,UniRef50_H7BYX6,
Funhe2Exx11m000488t2    Funhe2E6bm000458t39     Funhe2Exx11m000488      2       althi1  1987,83%,complete       99/100/.        aaref:1510,UniRef50_Q9Y3S1,
Funhe2Exx11m000564t2    Funhe2Eq7m075704t1      Funhe2Exx11m000564      2       althi1  1958,90%,complete       99/100/.        aaref:2055,UniRef50_Q4ZG55,
Funhe2Exx11m000573t2    Fungr1EG3m000186t4      Funhe2Exx11m000573      2       althi1  1898,82%,partial3,36X   99/99/./Funhe2E6bm000540t5      aaref:1295,UniRef50_A7E2V4,
Funhe2Exx11m000617t2    Funhe2E6bm000599t1      Funhe2Exx11m000617      2       althi1  1905,73%,complete       99/100/./Funhe2Eq7m083795t3     aaref:1905,UniRef50_A3KMH1,
Funhe2Exx11m000699t2    Fungr1EG3m000233t5      Funhe2Exx11m000699      2       althi1  844,95%,partial3        98/100/./Funhe2Eq7m069285t1     aaref:507,UniRef50_Q14781,
Funhe2Exx11m000706t2    Fungr1EG3m000238t4      Funhe2Exx11m000706      2       althi1  1827,91%,complete       99/100/.        aaref:1439,UniRef50_M7BCT0,
Funhe2Exx11m000712t2    Funhe2E6bm000688t16     Funhe2Exx11m000712      2       althi1  1804,85%,complete,71X   99/100/./Funhe2E6bm000689t1     aaref:668,UniRef50_Q5M7N9,refok,
Funhe2Exx11m000729t2    Funhe2E6bm000704t1      Funhe2Exx11m000729      2       althi1  1807,85%,complete       99/100/.        aaref:1852,UniRef50_M3ZVQ7,refbest,
 
#..
# eg main misses: Funhe2E6bm000458t1|Funhe2E6bm000540t1
## not in kfish2evg367mixx-kf2.gmap.gff.gz 
## not in kfish2evg367mixx.mrna_pub-kf2.gmap.out13.gz, but alts are
## .. these were drops in publicset1 : got re-added thru kfish2evg367mixx.fungr.dupids 
## ** what to do now? substitute Fungr equiv? gmap the Funhe alts?
# egrep '^(Funhe2E6bm000458t1|Funhe2E6bm000540t1)     ' ../kfish2evg367mixx.trclass1 
# Funhe2E6bm000458t1      drop    althi1  Fungr1EG3m000150t1      100/100 2053,48%,complete-utrbad        0,0,aadup:Fungr1EG3m000150t1,pflag:14
# Funhe2E6bm000540t1      drop    althi1  Fungr1EG3m000186t1      99/100/Funhe2E6bm000540t5       1950,73%,complete       0,0,aadup:Fungr1EG3m000186t1,pflag:8
#
## Funhe2E6bm000458t == Funhe2Exx3m000505t; Funhe2E6bm000540t == Funhe2Exx3m000592t
# egrep '(Funhe2E6bm000458t1|Funhe2E6bm000540t1)      ' ../gmapz/kf2pub11.idtabxy3m
#Funhe2E6bm000458t1      Funhe2E6bm000458t1      Funhe2Exx11m000488      1       main    2053,48%,complete-utrbad        99/100/.        aaref:1528,UniRef50_Q9Y3S1,
#Funhe2E6bm000540t1      Funhe2E6bm000540t1      Funhe2Exx11m000573      1       main    1950,73%,complete       99/100/.        .


## run1
# nmain=101221, nalt=239572, nerr=; 
# wc -l kf2pub11_both2a*
#  3661963 kf2pub11_both2a.alt.gff
#    96955 kf2pub11_both2a.equalalt.tab << MISSING mains
#   756786 kf2pub11_both2a.main.gff
# grep -c mRNA kf2pub11_both2a*.gff
# kf2pub11_both2a.alt.gff:239572
# kf2pub11_both2a.main.gff:101221

=cut
