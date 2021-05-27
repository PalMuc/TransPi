#!/usr/bin/env perl
# fa2pairfa.pl : paste _1.fa, _2.fa to one file .fa2 (assumes _1,_2 have identical entries, same order)

use strict;
use Getopt::Long;

# change to getOptions()?
my $isFQ= $ENV{fastq}||0;
my $addnum=$ENV{addnum}||0;
my $joinpair= $ENV{joinpair}||0; # paste mate seqs as one line >hdr\n{leftseq}{rightseq}\n; put mate-len on hdr
my $joinoneline= $ENV{join1}||0; 
my $outdir;

my $optok= GetOptions(
  "fastq|fq!", \$isFQ, 
  "joinpair!", \$joinpair, 
  "oneline!", \$joinoneline, 
  "addnum!", \$addnum, 
  "outdir=s", \$outdir, 
  );
unless($optok and @ARGV) {
  die "usage: fa2pairfa.pl [opts] reads_1.(fa|fq)[.gz]\n
  opts: -fastq input format; -oneline all one line with header, -joinpair as header line + one seq line, not 2,  -addnum /1 /2 to id
  convert input pair files reads_1.fasta, reads_2.fasta to one fa2 file\n";
}

$addnum=0 if($joinpair);
$joinpair=1 if($joinoneline);

$outdir="" if($outdir and -f $outdir and not -d $outdir);
mkdir($outdir) if($outdir and not -d $outdir);
  
foreach my $fain ( @ARGV ) {
  unless($fain and -f $fain) { warn "# missing $fain\n"; next; }
  
  # need this?
  if($fain =~ /fastq|\.fq/ and not $isFQ) { warn "# Filename $fain says fastq/fq; reset to isFQ=1\n"; }
  elsif($fain =~ /fasta/ and $isFQ) {  warn "# Filename $fain says fasta; reset to isFQ=0\n";  }
  
  my $paired=0; my $fa2in=""; 
  if($fain =~ /_1\.|\W1\.fa|\W1\.fq/) { 
    $fa2in= $fain;  $fa2in=~ s/_1\./_2./; $fa2in=~ s/1\.(f[aq])/2.$1/;
    $paired= (-f $fa2in)?1:0;
  } elsif($fain =~ /_2\.|\W2\.fa|\W2\.fq/) {
    warn "# want only _1.fa or .fq here; skipping input $fain\n"; next;
  }
  unless($paired) { warn "# SKIPPING unpaired input: $fain\n"; next; }
  ## REMOVE paired tests below, allow only this..
  
  (my $fname=$fain) =~ s/.gz//; 
  $fname=~s/\.(fast[aq]|f[aq])//;
  if($outdir) { $fname =~ s,^.*/([^/]+)$,$outdir/$1,; }
  
  my $outf= ($joinpair) ? "$fname.fa2join" : "$fname.fa2";
  if( -s $outf ) {  warn "# already have $outf, wont rewrite \n"; next; }
  my($ok,$nin,$nskip,$nerr,$npair)= (0) x 9;
  if($fain=~/.gz/) { $ok=open(F,"gunzip -c $fain |"); } else { $ok=open(F,$fain); }
  if($fa2in=~/.gz/) { $ok=open(F2,"gunzip -c $fa2in |"); } else { $ok=open(F2,$fa2in); }
  open(OUT,">$outf") or die "writing $outf";
  
  my($fh,$fs,$lh,$fhr,$fsr,$rh, $qh,$qs, $qhr,$qsr)=("") x 10;
  while(<F>) { 
    $fh=$_; $fs=<F>; ($lh)= $fh=~/^.(\S+)/; $nin++;   
    $fhr=<F2>; $fsr=<F2>; ($rh)= $fhr=~/^.(\S+)/;   
    if($isFQ) { 
      $fh=~s/^./>/;  $qh=<F>; $qs=<F>;  # read/drop  # fh=='@header'
      $fhr=~s/^./>/; $qhr=<F2>; $qsr=<F2>;  # read/drop
      }
    if($addnum) {
      unless($fh=~m,/1 ,) { $fh=~s, ,/1 ,; $fhr=~s, ,/2 ,, unless($fhr=~m,/2 ,); }
    }
    if($lh ne $rh) { 
      $lh=~s,/[12],,; $rh=~s,/[12],,; 
      if($lh ne $rh) { $nerr++;   warn "# mismatch pair: $lh ne $rh\n";  
        die "ERR: too many mismatch:$nerr" if($nerr>99); 
        next;
      } 
    }
    
    if($joinpair) {  # paired
      chomp($fs);  my $ls=length($fs);
      chomp($fsr); my $rs=length($fsr);
      unless($ls == $rs) { my($rs1,$rse)=($ls+1,$ls+$rs); 
         $fh =~ s/$/ reads12=1..$ls,$rs1..$rse len2=$rs/; }
      $fh =~ s/$/ len1=$ls/;   #? reads12=1..$ls,$ls+1..$ls+$rs
      ## opt to put fh,fs,fsr all 1 line for sort -u removal of duplicate read pairs
      if($joinoneline) { chomp($fh); print OUT join("\t",$fh,$fs,$fsr)."\n"; }
      else { print OUT $fh,$fs,$fsr,"\n"; } $npair++;
    } else {
      print OUT $fh,$fs; 
      print OUT $fhr,$fsr; $npair++;
    }
  } 

  close(OUT); close(F); close(F2) if($paired);
  warn "# done: nin=$nin; npair=$npair, nskip=$nskip, nerr=$nerr for $fain,$fa2in to $outf\n";  
}
  
__END__

$evg/fa2pairfa.pl -fastq -join ../fqx/Dman_32_ACAGTG_L003_1.fastq.gz
 # done: nin=71868056; npair=71868056, nskip=0, nerr=0 for ../fqx/Dman_32_ACAGTG_L003_1.fastq.gz,../fqx/Dman_32_ACAGTG_L003_2.fastq.gz to ../fqx/Dman_32_ACAGTG_L003_1stq.fa2join

sort -k2,2 Dman_32_ACAGTG_L003_1stq.fa2join | perl -ne \
'($h,$sq)=split"\t"; unless($sq eq $lsq) { $sl=substr($sq,0,101); $sr=substr($sq,101,101); 
print join("\n",$h,$sl,$h,$sr)."\n"; $nu++; } else { $nd++; } $lsq=$sq; END{ warn "nuniq=$nu; ndup=$nd\n"; } ' \
> Dman_32_ACAGTG_L003_1stq.fa2uniq

nuniq=71867906; ndup=150  ** not worth effort for so few dups **
