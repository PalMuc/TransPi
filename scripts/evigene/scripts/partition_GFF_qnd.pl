#!/usr/bin/env perl
# partition_GFF_qnd.pl from partition_GFF_genome.pl
# d.gilbert: quick and dirty form of partition_GFF_add.pl
# ?? no fasta

# FIXME: 2011sep: not handling overlapped partitions; 2nd over (start) is skipped.

use warnings;
use strict;
use FindBin;
use Getopt::Long; # qw(:config no_ignore_case bundling);
use Cwd;
use Carp;
use lib ("$FindBin::Bin/../PerlLib", "$FindBin::Bin"); # use Bin/ as alternate loc for Fasta_reader 
#? use Fasta_reader; # this is lightweight PASA/EvmUtils module; keep
use File::Basename;

umask(0000);

my $param_string = "HELP: see partition_GFF_add.pl\n";

my ($genome, @GFF, $type, $debug,
    $partitions_file, $accession_list, $help,
    );

# $accession_list="accession.list";
$type= "GFF";

my $optok= &GetOptions (
             'GFF=s', \@GFF,
             'type=s', \$type,
             'debug!', \$debug,
             'genome=s', \$genome,
             'partitions_file=s' => \$partitions_file,
             'accession_list=s' => \$accession_list,
             'h' => \$help,
             );

push(@GFF, @ARGV); # any remainder;

# process required options
if ($help) { die $param_string; }
unless ($optok && @GFF>0 && $partitions_file) {
    die $param_string;
}


my %base_directories_to_partitions;
my %base_directories_info;
open (my $fh, $partitions_file) or die "Error, cannot open $partitions_file";
while (<$fh>) {
    chomp;
    my ($accession, $base_dir, $partitioned, $partition_dir) = split (/\t/);
    $base_directories_info{$base_dir}{haspart}= ($partitioned eq 'Y');
    $base_directories_info{$base_dir}{accession}= $accession;
    if ($partitioned eq 'Y') {
        $base_directories_to_partitions{$base_dir}->{$partition_dir} = 1;
    }
}
close $fh;


my @files_to_process = ();
if( $type eq 'fagff' ) {
  ## not yet
#   my $ingff= shift @GFF;
#   my $infasta= shift @GFF; # FIXME option
#   my $fasta_list= get_fasta_list($infasta);
#   push(@files_to_process, 
#     {  type => "fagff", 
#        basename => basename($infasta), # output part fasta, not GFF
#        file => $ingff, # input file
#        fasta_list => $fasta_list,
#    } ) ;
  
} else { 
  foreach my $gf (@GFF) {
    unless(-f $gf) { die "missing GFF=$gf\n"; next; }
    
    push(@files_to_process, 
      {  type => "GFF", # nice idea but not used
         basename => basename($gf),
         file => $gf,
     } ) ;
  }
  
}

# my %partinfo;
my $partinfo= getPartInfo(); # base_directories_info

foreach my $file_struct (@files_to_process) {
  my $in_filename = $file_struct->{file};
  my $basename = $file_struct->{basename};
  my $type = $file_struct->{type};
  my $in_fasta_list = $type eq 'fagff' ? $file_struct->{fasta_list} : undef;
  
  warn "# process $type input $in_filename\n" if $debug;
  if($type =~ /^gff/i) {
    gff_partitioner( $in_filename, $partinfo, $basename);
  }
  
}

# exit(0);


#...............


sub getPartInfo
{
  my $partinfo= {};
  
  foreach my $base_dir (sort keys %base_directories_info) {
  
    my $accession = $base_directories_info{$base_dir}{accession};      
    unless (-d $base_dir) {
      die "Error, missing $accession folder: $base_dir";
    }
  
    my $haspart= $base_directories_info{$base_dir}{haspart};
  
    if( $haspart ) {
      my $partition_dirs_href = $base_directories_to_partitions{$base_dir};
      my $max=0; my $lmax_rend=0; my $lmin_lend=0;
      my @cstarts=();
      my @cends=();
      foreach my $partition_dir (sort keys %$partition_dirs_href) {
        $partition_dir =~ /(\w+)_(\d+)-(\d+)$/ or die "Error, cannot extract coords from $partition_dir";
        my ($segref, $min_lend, $max_rend) = ($1,$2,$3);
        # urk can have overlap where next min_lend < last max_rend
        # for this to work sort above needs to be numeric
        my $isover= ($min_lend > $lmin_lend and $lmax_rend > $lmin_lend and $min_lend < $lmax_rend) ? 1:0;
        push(@cstarts, $min_lend);
        push(@cends, $max_rend);
        $max= $max_rend if($max_rend > $max);
        $partinfo->{$accession}{partition_dir}{ $min_lend } =  $partition_dir;
        $partinfo->{$accession}{overlapped}{ $min_lend } = $isover;
        $partinfo->{$accession}{overlapped}{ $lmin_lend } ++ if($isover);
        $lmax_rend= $max_rend; $lmin_lend= $min_lend;
        } 
        ## need to have @{parts} ordered by min_lend ... or sort below
      # @cstarts = sort { $a <=> $b } (@cstarts,$max);
      my @ord = sort { $cstarts[$a] <=> $cstarts[$b] } (0..$#cstarts);
      @cstarts= @cstarts[@ord];
      @cends  = @cends[@ord];
      
      $partinfo->{$accession}{cstarts} = \@cstarts;
      $partinfo->{$accession}{cends} = \@cends;
    }
  
          
    else {  # no parts
      ## my $part_file = "$base_dir/$basename";
      my %accs=();
      
      if($accession =~ /^combo_/) {
        my $lh;
        if( $accession_list && open( $lh, "$base_dir/$accession_list") ) {  
          while(<$lh>){  s/^>//; if(m/^(\S+)/){ $accs{$1}++;} } close($lh); 
          # $accession= \%accs;
        } elsif( $genome && open( $lh, "$base_dir/$genome") ) { 
          while(<$lh>){ if(m/^>(\S+)/){ $accs{$1}++;} } close($lh); 
          # $accession= \%accs;
        } else {
          die "missing --accession_list or --genome for combo";
        }
        
        foreach my $acc (sort keys %accs) {
          $partinfo->{$acc}{partition_dir}{ 1 }= $base_dir;
          $partinfo->{$acc}{overlapped}{ 1 }= 0;
          $partinfo->{$acc}{combo}= $accession;
        }
        
      } else {
        $partinfo->{$accession}{partition_dir}{ 1 }= $base_dir;
        $partinfo->{$accession}{overlapped}{ 1 }= 0;
      }
      
    }
    
  }
  
  return $partinfo;
}



sub gff_partitioner {
  my($ingff,$partinfo,$basename)= @_;
  
  my %outhands= ();
  my($lcontig);
  
  # my($outgff,$seq_id,$min_lend,$max_rend,$adjust_to_1)= (0) x 10;
  
  open(my $inh,"$ingff") or do { warn "cant read $ingff\n"; return; };  
  my $got_spacer = 0;
  my $comments="";
  
  while (<$inh>) {
  
      unless (/\w/) { 
          # print $outh $_ unless ($got_spacer);
          $got_spacer = 1;
          next;
      }
  
      if (/^\#/) {
        $comments .= $_;
        # print $outh $_;
        next;
      }
  
      my @x = split (/\t/);
      my ($contig, $lend, $rend) = ($x[0], $x[3], $x[4]);
            
      my $partition_dir="";
      my $adjust_to_1= 0;
      my $adjust_coord = 0;
      my $min_lend= 0;

# FIXME: not handling overlapped partitions; 2nd over (start) is skipped.
# .. last<< bad, need to look for 2 parts that overlap lend,rend or flag if overlapping span
      my $overlapped= 0;
      my $partition_dir2="";
      my $adjust_coord2 = 0;
      
      if($partinfo->{$contig}{cstarts}) { # {haspart}
        my @cstarts= @{ $partinfo->{$contig}{cstarts} };
        my @cends  = @{ $partinfo->{$contig}{cends} }; # bad
        for(my $i=0; $i <= $#cstarts; $i++) { 
          if( $lend >= $cstarts[$i] && $rend <= $cends[$i] ) {
           $min_lend= $cstarts[$i];
           $partition_dir= $partinfo->{$contig}{partition_dir}{ $min_lend }; 
           $overlapped= $partinfo->{$contig}{overlapped}{ $min_lend }; 
           #^ over next or last ??
           $adjust_to_1= 1;
           $adjust_coord = $min_lend - 1;
           if($overlapped and $i < $#cstarts) {
            my $j=$i+1;
            if( $lend >= $cstarts[$j] && $rend <= $cends[$j] ) {
              $adjust_coord2 = $cstarts[$j];
              $partition_dir2= $partinfo->{$contig}{partition_dir}{ $adjust_coord2 }; 
              $adjust_coord2--;
              }
            }
           last;
           }
          }
      } else {
        $min_lend = 1;
        $partition_dir= $partinfo->{$contig}{partition_dir}{ $min_lend };       
        $overlapped= $partinfo->{$contig}{overlapped}{ $min_lend }; 
        # $adjust_to_1= 0;
        # $adjust_coord = 0;
      }
      
      unless( -d $partition_dir ) { 
        # die "bad path for $contig:$lend-$rend: $partition_dir\n";
        # got some hints w/ huge spans; ignore if few
        warn "bad path for $x[2]:$x[1] $contig:$lend-$rend: $partition_dir\n";
        next; ## return; # next/return/die ?
      } 
 
# scaffold_10	/N/gpfsbr/gilbertd/chrs/aug/dpulex/genoparts/scaffold_10	Y	
#  /N/gpfsbr/gilbertd/chrs/aug/dpulex/genoparts/scaffold_10/scaffold_10_1-2000000
# scaffold_10	/N/gpfsbr/gilbertd/chrs/aug/dpulex/genoparts/scaffold_10	Y	
#  /N/gpfsbr/gilbertd/chrs/aug/dpulex/genoparts/scaffold_10/scaffold_10_1950001-2169786
# bad path for scaffold_10:1949995-1950004: 

 
      # close outh after pass contig/ref ? assume ingff sorted by ref?
      if($lcontig and $contig ne $lcontig) {
        my $comacc= $partinfo->{$lcontig}{combo} || $lcontig;
        foreach my $file (sort keys %{$outhands{$comacc}}) { 
          next if ($file =~ m,combo,);
          ## this is bad for combo's
          warn "# close $file\n" if $debug;
          my $outh= delete $outhands{$comacc}{$file};
          close($outh) if $outh; 
        }
      }
      
      $lcontig= $contig;
      
      my @savex= @x;
      my $kend= ($overlapped && $partition_dir2) ? 2 : 1;
      for (my $kover= 1; $kover <= $kend; $kover++) {
        if($kover == 2) {
          $partition_dir= $partition_dir2;
          $adjust_coord= $adjust_coord2;
          @x= @savex;
        }
        my $part_file = $partition_dir . "/$basename";
        my $comacc= $partinfo->{$contig}{combo} || $contig;
  
        ## FIXME: need contig -> combo_n map here for files
        my $outh= $outhands{$comacc}{$part_file};
        unless($outh) {
          warn "# at $comacc, open $part_file\n" if $debug;
          open( $outh,">$part_file") or do { die "ERR: at $comacc, cant write $part_file\n"; return; };
          $outhands{$comacc}{$part_file}= $outh;
          print $outh $comments if $comments; #? only at top
        }
  
        if ($adjust_to_1) {
          $x[3] -= $adjust_coord;
          $x[4] -= $adjust_coord;
        }
        
        print $outh join ("\t", @x);
      } # kover
      
      $got_spacer = 0;
  }
  
  # close($outh); 
  close($inh);
  foreach my $contig (sort keys %outhands) { 
    foreach my $file (sort keys %{$outhands{$contig}}) { 
      warn "# close $file\n" if $debug;
      my $outh= delete $outhands{$contig}{$file};
      close($outh) if $outh; 
      }
    }
}


