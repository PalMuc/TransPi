#!/usr/bin/env perl
# fagff2sam.pl

use strict;
use warnings;
use Getopt::Long ;

my($ok, $fasta, $gff, $skiptypes, $debug);
my $dropdit=1; # dang id mess \.1 at end on some

my $optok= &GetOptions (
  "gff=s"=>\$gff,  
  "fasta=s"=>\$fasta,
  "skiptypes=s"=>\$skiptypes,
  "n|debug!"=>\$debug, # now default
  );

$fasta= shift @ARGV unless($fasta);
$gff  = shift @ARGV unless($gff);
die "USAGE: $0 -fasta my.fasta[.gz] -gff stdin|my.gff[.gz] > my.sam\n"
  unless($optok and $fasta and $gff);

if($fasta =~ /\.gz/){ $ok= open(FA,"gunzip -c $fasta |"); }
else { $ok= open(FA,$fasta); }
die $fasta unless $ok;

my $gfh;
if($gff =~ /\.gz/){ $ok= open(GF,"gunzip -c $gff |"); $gfh=*GF; }
elsif($gff =~ /stdin|-/){ $ok= 1; $gfh=*STDIN; }
else { $ok= open(GF,$gff); $gfh=*GF; }
die $gff unless $ok;

# dang .fa has ID=xxx.1 ; .gff has ID=xxx
my($id, %fa);
while(<FA>){ chomp; if(/^>(\S+)/) { $id=$1; $id=~s/\.\d$// if($dropdit);} else { $fa{$id}.=$_; } }
close(FA);

# handle paired seqs? if sameid, but for  some id suffix like .fwd .rev

# sam line
my($qid, $flag, $chr, $cloc, $mapq, $cigar, $matechr, $mateloc, $isize, $seq, $qual, @opt);
    $mapq=255;
    $isize=0; #? pair insert size
    $qual="*"; @opt=();
my(%did);

while(<$gfh>){
  next if($skiptypes and m/\t($skiptypes)\t/);
  if(/\bID=([^;\s]+)/) { 
    $qid=$1;  $qid=~s/\.\d$// if($dropdit);
    next if($did{$qid}++);
    $seq= $fa{$qid} or next;
    my @gff=split"\t";
    $flag=0;   #1=pair, 2=pairok, 4=mismatch 8=mismate, 16=rev, 32=revmate, 64=firstmate, 128=secondmate
    $flag |= 16 if($gff[6] eq "-");
    $chr=$gff[0]; $cloc=$gff[3];
    $cigar= length($seq)."M"; # add introns ? need gff match_part/exon/...
    $matechr="*"; $mateloc=0;
    @opt=();
    push(@opt,"XS:A:".$gff[6]) if($gff[6] =~ m/[+\-]/); # strand tag
    # ... if($paired and $qid =~ m/$pairpatt/) { }
    print join("\t",($qid, $flag, $chr, $cloc, $mapq, $cigar, $matechr, $mateloc, $isize, $seq, $qual, @opt)),"\n";
  } 
}

# x.8244810       16      Scaffold2       867     255     36M     *       0       0       CTAGACACCAAAAAAATAGCAGCGGCTATATGATTT      *
# x.1350890       0       Scaffold2       868     255     15M1I20M        *       0       0       TAGACACCAAAAAAAATAGCAGCGGCTATATGATTT      *
