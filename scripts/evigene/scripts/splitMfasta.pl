#!/usr/bin/env perl
#
# splitMfasta.pl
# split a multiple fasta file in smaller multiple fasta files
# 
# 
# Mario Stanke, 2.04.2007
# d.g.gilbert small mods, 2016 : input.seq.gz ok, leave input.suffix on split names, add --nparts=ncpu
#

use strict;
use Getopt::Long;
use File::Basename qw(basename);
use File::Spec::Functions qw(rel2abs);

my $usage = "$0 -- split a multiple fasta file in smaller multiple fasta files.\n\n";
$usage .= "Usage: $0 input.fa\n\n";
$usage .= "parameters:\n";
$usage .= "--minsize=n       each split output fasta file total to at least this man base pairs.\n";
$usage .= "                  set this to 0 to split the input in single sequence files.\n";
$usage .= "--nparts=n        split into this many parts (calculates minsize).\n";
$usage .= "--outputpath=s    prefix to output path. Output files are\n";
$usage .= "                  s/input.split.1.fa\n";
$usage .= "                  s/input.split.2.fa\n";
$usage .= "                  s/input.split.3.fa\n";
$usage .= "                  ... \n";
$usage .= "\n";


if (@ARGV < 1) {
    print $usage;
    exit;
}

my ($seqfilename, $ok, $nparts, $minsize, $outputpath,$cursize, $idx, $basename);
my $minsize = 0;
my $outputpath = "";

GetOptions( 'minsize:i'=>\$minsize, 'nparts:i'=>\$nparts, 'outputpath:s'=>\$outputpath);

$seqfilename = $ARGV[0];
if($outputpath =~ /^~/){
    my $HOME = $ENV{'HOME'};
    $outputpath=~s/~/$HOME/;
}
 
$outputpath=rel2abs($outputpath);
if ($outputpath ne "" && $outputpath !~ /\/$/){
    $outputpath .= "/";
}

if($nparts > 0) { # dgg: calc minsize, from seqfile[.gz]
  # splitsize=`grep -v '^>' $protin | wc -c | sed 's/ .*//' `
  # splitbp=$(( $splitsize / $ncpu ))
  my($nb,$nseq)=(0,0);
  if($seqfilename =~ /\.gz/) { $ok= open(SEQ, "gunzip -c $seqfilename |"); } else { $ok=open(SEQ, "<$seqfilename")  }
  while(<SEQ>) { if(/^>/) { $nseq++; } else { $nb+= length($_); } } close(SEQ);
  $minsize= int($nb/$nparts);
}

die("missing: -minsize or -nparts") if($minsize == 0);  # dgg: fail

if($seqfilename =~ /\.gz/) { $ok= open(SEQ, "gunzip -c $seqfilename |"); } else { $ok=open(SEQ, "<$seqfilename")  }
die ("Could not open $seqfilename") unless($ok);

$idx=0;
$cursize=0;
my $basename = basename($seqfilename);
#dgg: leave# $basename =~ s/(\.fa|\.fna|\.fasta)$//;

while (<SEQ>) {
    if ($idx == 0 || $cursize > $minsize && /^>/){
	$idx++;
	if ($idx > 1) {
	    close(SPLIT);
	}
	open (SPLIT, ">$outputpath$basename.split.$idx.fa") or die ("Could not open $basename.split.$idx.fa");
	$cursize = 0;
    }
    print SPLIT;
    if (!/^>/) {
	$cursize += length $_;
    }
}
close(SEQ); close(SPLIT);

