#!/usr/bin/env perl

## from recombine_EVM_partial_outputs.pl
## this one just reads gff output from partitions, and adjusts location to full genome/scaffold locs.
## really need option to rewrite IDs from partial files which can be same when combined
## for both gff and fasta, other?

## see aug2mapgff.pl
## gff needs to understand gene model types (mRNA>CDS,exon,etc get same id)

# ** FIXME: add opts to handle overlapped regions: at least drop or keep 2nd overlap area;
#    best would to combine unique non-overlap features (GFF) and drop full/part identicals
# note also 

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long ;
# qw(:config no_ignore_case bundling);
use File::Basename;

my ($SEE, %seenid, %oldid, $type, $help, $partitions_file, $toplevel, $output_file_name);
my $uniqueid=1;
my $skipover=0; # should become flag how to handle overlapped parts
my $SKIPCOM= "#skip."; # opt to skip entirely or comment out overlapped extras

my $usage = "usage: $0 --type [fasta|gff|list] [--skipover --nouniqueid] --partitions part_file --output_file_name file_name\n\n";

my $optok = &GetOptions ("partitions=s" => \$partitions_file,
             "type=s" => \$type,
             "uniqueid!" => \$uniqueid,
             "skipover!" => \$skipover, # add sk opts: skip entirely or print commented gff
             "SKCOMMENT=s" => \$SKIPCOM, #  
             "output_file_name|O=s" => \$output_file_name,
             "toplevel:s" => \$toplevel,
             "verbose|S" => \$SEE,
             "help|h" => \$help,
             );

#?? check/set type from output_file_name suffix?

if ($help || ! ($optok && $type && $type =~ /list|gff|fa/i && $partitions_file && $output_file_name) ) {
    die $usage;
}

# if(defined $toplevel) {
#  unless($toplevel) { $toplevel= cwd(); }
# }

$output_file_name = basename($output_file_name); #just to be sure...

if($type =~ /gff/i && $output_file_name !~ /\.gff/i) {
  warn "WARN: type=$type doesn't match output suffix \n";
} elsif($type =~ /fa/i && $output_file_name !~ /\.fa/i) {
  warn "WARN: type=$type doesn't match output suffix \n";
}

my %skiprange=();
my $skiprange= []; #start,end of overlap to skip
my ($lsegref, $lrange_lend, $lrange_rend)=(0) x 3;
my %partlist=();

my %base_directories_to_partitions;
open (my $fh, $partitions_file) or die "Error, cannot open $partitions_file";
while (<$fh>) {
    chomp;
    my ($accession, $base_dir, $partitioned, $partition_dir) = split (/\t/);
    $partlist{$base_dir} = $accession; #? want accession?
    if ($partitioned eq 'Y') {
        $base_directories_to_partitions{$base_dir}->{$partition_dir} = 1;

      # for skipover
        if( $partition_dir =~ /(\w+)_(\d+)-(\d+)$/ ) {
        my ($segref, $range_lend, $range_rend) = ($1,$2,$3);
        if($segref eq $lsegref && $range_lend < $lrange_rend) {
          $skiprange{$partition_dir}= [ $range_lend, $lrange_rend-1 ];
          }
        ($lsegref, $lrange_lend, $lrange_rend) =($segref, $range_lend, $range_rend);
        }
    }
}
close $fh;

foreach my $base_dir (sort keys %base_directories_to_partitions) {
    my $partition_dirs_href = $base_directories_to_partitions{$base_dir};

    my $new_outfile= "$base_dir/$output_file_name";
    $new_outfile .= ".list" if($type =~ /list/i);
    open (my $out_fh, ">$new_outfile") or die $!;
    print STDERR "\twriting output $new_outfile\n\n" if $SEE;

    my $part=0;
    foreach my $partition_dir (sort keys %$partition_dirs_href) {
#         $partition_dir =~ /(\d+)-(\d+)$/ or die "Error, cannot extract coords from $partition_dir";
#         my ($range_lend, $range_rend) = ($1, $2);
        $partition_dir =~ /(\w+)_(\d+)-(\d+)$/ or die "Error, cannot extract coords from $partition_dir";
        my ($segref, $range_lend, $range_rend) = ($1,$2,$3);
        $part++;
        
        if($skipover) { $skiprange= $skiprange{$partition_dir } || []; }
        
        my $part_file = $partition_dir . "/$output_file_name";
        if($type =~ /gff/i) {
          &relocate_gff($out_fh, $part_file, $segref, $range_lend, $range_rend,  $part);
        } elsif($type =~ /fa/i) {
          &cat_fasta($out_fh, $part_file, $segref, $range_lend, $range_rend, $part);
        } elsif($type =~ /list/i) {
          &list_file($out_fh, $part_file, $segref, $range_lend, $range_rend, $part);
        } 
    }
    
    close $out_fh;
}

if( defined $toplevel ) {
  my $new_outfile= "$toplevel$output_file_name";
  $new_outfile .= ".list" if($type =~ /list/i);
  open (my $out_fh, ">$new_outfile") or die $!;
  print STDERR "\twriting output $new_outfile\n\n" if $SEE;
  foreach my $base_dir (sort keys %partlist) {
        my $part_file = $base_dir . "/$output_file_name";
        if($type =~ /list/i) {
          if(-f $part_file) { print $out_fh $part_file,"\n"; }
          if(-f "$part_file.list") {
          open (my $fh, "$part_file.list") or warn $!;
          while(<$fh>){ print $out_fh $_; }
          close($fh);
          }          
        } else {  # cat files
          open (my $fh, $part_file) or warn $!; #? die or warn/return
          while(<$fh>){ print $out_fh $_; }
          close($fh);
        }
  }
  
  # --strip-components=N to remove leading /full/path/to/maindir ... ? do in .list?
  print STDERR "# tar -cf $output_file_name.tar --files-from=$new_outfile\n" if ($type =~ /list/i); # SEE
}


exit(0);


####
### add option to tar parts w/o need to rewrite locations in all formats, e.g. blast output
##   1. create file part list; 2. run tar w/ list (outside this?)
##      tar -cf myoutput.tar --files-from=outparts.list

sub list_file {
    my ($out_fh, $part_file, $segref, $partition_lend, $partition_rend,  $part) = @_;
    # print STDERR "Listing $part_file\n";    
#     if($docomment) {
#     print $out_fh "#part_location: $segref:$partition_lend-$partition_rend\n";
#     print $out_fh "#part_file: $part $part_file\n\n";
#     }
    print $out_fh $part_file,"\n" if (-f $part_file); # non zero ? -s 
}


sub cat_fasta {
    my ($out_fh, $part_file, $segref, $partition_lend, $partition_rend,  $part) = @_;
    print STDERR "Parsing $part_file\n" if $SEE;    
    open (my $fh, $part_file) or warn $!; #? die or warn/return
    while (<$fh>) {
      # may want to rewrite fasta ids ? other?
      if($uniqueid && /^>(\S+)/) {
        my $id= $1;
        my $newid= uniqid_gff("ID=$id", $segref, $part);
        $newid =~ s/ID=//;
        s/^>(\S+)/>$newid/;
      }
      print $out_fh $_;
    }
    close $fh;
}

sub uniqid_gff {
  my($at, $segref, $part)= @_;

  #? just use $part and $segref to make uniq ?
  my($id)= $at =~ m/ID=([^;\s]+)/; 
  if($id) { 
    my $newid= "p$part.$id";
    $newid= $segref.$newid unless($id =~ /$segref/);
    $oldid{$id}{$segref}{$part}= $newid;
    $at =~ s/ID=([^;\s]+)/ID=$newid/;
    }
    
  # also handle Parent= ; Target= ?
  if($at =~ m/Parent=([^;\s]+)/) { 
    my $par=$1; # can have 2+ parents ... split ",", $par;
    my $newid= $oldid{$par}{$segref}{$part};
    if($newid) { $at =~ s/Parent=([^;\s]+)/Parent=$newid/; }
    }
   
  return $at;
}

sub relocate_gff {
    my ($out_fh, $part_file, $segref, $partition_lend, $partition_rend, $part) = @_;
    
    my $didline=0; my $didskip= 0;
    print STDERR "Parsing $part_file\n" if $SEE;    
    open (my $fh, $part_file) or do {warn "bad part: $part_file: ". $!; return; }; #? die or warn/return
    while (<$fh>) {
      my $g=$_;
      if(/^\w/ and /\t/) {
        my @x = split (/\t/, $g); $didskip= 0;
        $x[3] += $partition_lend -1; # gff start
        $x[4] += $partition_lend -1; # gff end
        $x[8]= uniqid_gff($x[8], $segref, $part) if $uniqueid;
        $g= join("\t",@x);
        if( $skipover && @{$skiprange}>0 && _inside( @x[3,4], @$skiprange) ) { 
          ## print notice, parsable
          # $g=""; $g="#skipover: ".join(",",@x)."\n" if($x[2] =~ /^(gene|mRNA)/);
          ## fixme2: print all of skipover, commented out.
          $g= $SKIPCOM.$g; $didskip=1;
          }
        }
      elsif($didskip) { $g= $SKIPCOM.$g; } #add #sk. even to comments
          
      print $out_fh $g if($g);
      unless($didline++) { # comment where we are from
        print $out_fh "#part_location: $segref:$partition_lend-$partition_rend\n";
        print $out_fh "#part_file: $part $part_file\n\n";
      }
    }
    close $fh;
}

sub _inside {
  my($x,$y,$min,$max)= @_;
  return ($y <= $max && $x >= $min && $y >= $min && $x <= $max) ? 1 : 0;
}

