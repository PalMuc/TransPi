#!/usr/bin/env perl	
# env suf=fa2 list=bigsub.list velsubsetfa.pl  *.mapt.fa2.gz

$suf= $ENV{suf}||"fa";
$list=$ENV{list} or die "env list=xxx";
open(F,$list) or die "open $list"; 
while(<F>){ 
	if(/^subset\t(\w+)\t(\S+)/){ ($n,$v)=($1,$2); push @n,$n; 
	map{ $sc{$n}{$_}=1; }split",",$v; }
} close(F);

foreach $fa (@ARGV) {
 if($fa =~ /.gz$/) { open(FA,"gunzip -c $fa |") or die "open $fa"; }
 else { open(FA,$fa) or die "open $fa"; }
 while(<FA>) {
	if(/^>/) {
		($s)=m/loc=(\w+)/; $p=0; 
		foreach $n (@n) { if($sc{$n}{$s}) { $p=$n; last; } } 
		if($p) { $fp= $fps{$p}; 
			unless($fp){ my $fh; open($fh,">subset.$p.$suf"); $fps{$p}=$fp= $fh; } 
			print $fp $_; 
			} 
	} elsif($p) { print $fp $_; } 
 }
 close(FA);
}

