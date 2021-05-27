#!/usr/bin/env perl
## fq2spotfa.pl fastq pair files _1/2.fq to spot (join 1/2) fasta, 

$SAMESIZE=$ENV{samesize}||0; # default on?
$MAXERR=$ENV{maxerr}||99;

foreach $fq ( @ARGV ) {
  ($fqr=$fq) =~ s/_1\./_2./;
  if($fqr eq $fq) { warn "#err: skip not _1/_2 pair files: $fq\n"; next; }
  unless( -f $fq and -f $fqr) { warn "#err: missing $fq or $fqr\n"; next; }

  ($fn=$fq) =~ s/.gz//; $fn=~s/\.(fq|fastq)//; $fn=~s/_[12]/_sp/;
  $fao=$fn.".fa";
  my($ok,$nerr,$npair)= (0) x 9;
  if( -f $fao) { warn "#err: have alread $fao\n"; next; }

  warn "# fastq pair $fq,$fqr to spotfa $fao\n";
  if($fq=~/.gz/) { $ok=open(F,"gunzip -c $fq |"); } else { $ok=open(F,$fq); } die "read $fq" unless($ok);
  if($fqr=~/.gz/) { $ok=open(R,"gunzip -c $fqr |"); } else { $ok=open(R,$fqr); } die "read $fqr" unless($ok);
  open(L,">$fao") or die "write $fao"; 
  while(<F>) { 
    $fh=$_; $fs=<F>; $qh=<F>; $qs=<F>; 
    $rfh=<R>; $rfs=<R>; $rqh=<R>; $rqs=<R>; 
    $fh=~s/^\@//; $fh=~s,/[12],,; $rfh=~s/^\@//; $rfh=~s,/[12],,;
    chomp($fs); chomp($rfs); #? test same size? cant unsplit if not..
    if($fh ne $rfh) { $nerr++; $fs=$rfs=""; }
    if($SAMESIZE) { unless(length($fs) eq length($rfs)) { $nerr++; $fs=$rfs=""; } }
    if($fs and $rfs) { print L ">$fh",$fs,$rfs,"\n"; $npair++; }
    else { die "#err: too many readpair errors: $nerr\n" if($nerr > $MAXERR); }
    } 
  close(L); close(F);  close(R);
  warn "# done: nseq=$npair, nerr=$nerr in $fq to $fao\n";  
}
#  warn "usage: fq2fa.pl reads.fastq[.gz]\n"; 

