#!/usr/bin/env perl
## fq2fa.pl fastq to fasta, 

foreach $fq ( @ARGV ) {

($fn=$fq) =~ s/.gz//; $fn=~s/\.(fq|fastq)//;
$fao=$fn.".fa";
my($ok,$nerr,$npair)= (0) x 9;

if($fq and -f $fq) {
  next if( -f $fao);
  warn "# reduce $fq to $fao\n";
  if($fq=~/.gz/) { $ok=open(F,"gunzip -c $fq |"); } else { $ok=open(F,$fq); }
  open(L,">$fao") or die "write $fao"; 
  while(<F>) { 
    $fh=$_; $fs=<F>; $qh=<F>; $qs=<F>; ##($lh)=$fh=~/^.(\S+)/;  
    $fh=~s/^\@//;
    print L ">$fh",$fs; $npair++;
    } 
  close(L); close(F);  
  warn "# done: nseq=$npair for $fq to $fao\n";  
} else {
  warn "usage: fq2fa.pl reads.fastq[.gz]\n"; 
  exit(1);
}

}

