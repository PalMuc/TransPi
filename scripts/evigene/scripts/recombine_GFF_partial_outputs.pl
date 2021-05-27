#!/usr/bin/env perl

## from recombine_EVM_partial_outputs.pl
## this one just reads gff output from partitions, and adjusts location to full genome/scaffold locs.
## really need option to rewrite IDs from partial files which can be same when combined
## for both gff and fasta, other?
## see aug2mapgff.pl
## gff needs to understand gene model types (mRNA>CDS,exon,etc get same id)

# note also 

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long ;
# qw(:config no_ignore_case bundling);
use File::Basename;

my ($SEE, $nogff, $help, $partitions_file, $output_file_name);

my $usage = "usage: $0 [--nogff] --partitions part_file --output_file_name file_name\n\n";

&GetOptions ("partitions=s" => \$partitions_file,
             "output_file_name|O=s" => \$output_file_name,
             "verbose|S" => \$SEE,
             "nogff!" => \$nogff,
             "help|h" => \$help,
             );


if ($help || ! ($partitions_file && $output_file_name) ) {
    die $usage;
}

$output_file_name = basename($output_file_name); #just to be sure...


my %base_directories_to_partitions;
open (my $fh, $partitions_file) or die "Error, cannot open $partitions_file";
while (<$fh>) {
    chomp;
    my ($accession, $base_dir, $partitioned, $partition_dir) = split (/\t/);
    if ($partitioned eq 'Y') {
        $base_directories_to_partitions{$base_dir}->{$partition_dir} = 1;
    }
}
close $fh;

foreach my $base_dir (sort keys %base_directories_to_partitions) {
    my $partition_dirs_href = $base_directories_to_partitions{$base_dir};

    open (my $out_fh, ">$base_dir/$output_file_name") or die $!;
    print STDERR "\twriting output $base_dir/$output_file_name\n\n";

    foreach my $partition_dir (sort keys %$partition_dirs_href) {
        $partition_dir =~ /(\d+)-(\d+)$/ or die "Error, cannot extract coords from $partition_dir";
        my ($range_lend, $range_rend) = ($1, $2);
        
        my $part_file = $partition_dir . "/$output_file_name";
        if($nogff) {
        &cat_nongff($out_fh, $part_file, $range_lend);
        } else {
        &relocate_gff($out_fh, $part_file, $range_lend);
        }
    }
    
    close $out_fh;
}

exit(0);


####

sub cat_nongff {
    my ($out_fh, $part_file, $partition_lend) = @_;
    print STDERR "Parsing $part_file\n";    
    open (my $fh, $part_file) or warn $!; #? die or warn/return
    while (<$fh>) {
      # may want to rewrite fasta ids ? other?
      print $out_fh $_;
    }
    close $fh;
}

sub relocate_gff {
    my ($out_fh, $part_file, $partition_lend) = @_;
    
    print STDERR "Parsing $part_file\n";    
    open (my $fh, $part_file) or warn $!; #? die or warn/return
    while (<$fh>) {
      my $g=$_;
      if(/^\w/ and /\t/){
        my @x = split (/\t/, $g);
        $x[3] += $partition_lend -1; # gff start
        $x[4] += $partition_lend -1; # gff end
        $g= join("\t",@x);
        }
      print $out_fh $g;
    }
    close $fh;
}


