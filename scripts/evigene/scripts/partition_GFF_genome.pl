#!/usr/bin/env perl
# partition_GFF_genome.pl
# d.gilbert: generalized from partition_EVM_inputs
# input: genome.fasta + any GFF files
## add combine option: for scaffolds smaller than minSegmentSize, 
##  combine fasta, gff in one file (or n files of max segmentSize)

use warnings;
use strict;
use FindBin;
use Getopt::Long; # qw(:config no_ignore_case bundling);
use Cwd;
use Carp;
use lib ("$FindBin::Bin/../PerlLib", "$FindBin::Bin"); # use Bin/ as alternate loc for Fasta_reader 
use Fasta_reader; # this is lightweight PASA/EvmUtils module; keep
use Data::Dumper;
use File::Basename;

umask(0000);

my $param_string = <<__PARAMS__;

###############################################################################################

The evidence for each contig is partitioned into a separate contig directory.

If the contig is longer than the segmentSize value, the data will be further partitioned 
into segmentSize data sets.  

Contigs shorter than minSegmentSize will be skipped.  Alternately, --combineSmallSegments
will combine these into one fasta, GFF file set (or n files of <= segmentSize).

The output file (specified by --partition_listing) will contain a list of all directories containing
the partitioned data and the partitioned data ranges.


###############################################################################################
#
# Options:
#  --genome                * :fasta file containing all genome sequences
#  --GFF                   * :files containing gene predictions + evidence; repeat as needed
#                            : all extra arguments taken as GFF input data
#  --segmentSize           * :length of a single sequence for running EVM
#  --overlapSize           * :length of sequence overlap between segmented sequences
#  --minSegmentSize          :length minimum of sequence required (skip small scaffolds;dgg)
#  --skipSegmentList         file with list of segment IDs to skip (dgg)
#  --skipSegmentList         file with list of segment IDs to keep (dgg); overrides skiplist
#  --combineSmallSegments    flag to create file set combining small scaffolds (<minSegmentSize)
#
#  --partition_listing     * :name of output file to be created that contains the list of partitions
#  -h                        : help message with more detailed information.
#
#   (*required options)
#
################################################################################################

__PARAMS__

    ;

my ($genome, @GFF, $gene_predictions, $protein_alignments, $transcript_alignments,
    $pasaTerminalExons, $repeats, $segmentSize, $overlapSize,
    $minSegmentSize, $skipSegmentList, $combineSmallSegments,
    $keepSegmentList,
    $partition_listing, $help,
    );


my $optok= &GetOptions ('genome=s' => \$genome, 'GFF=s', \@GFF,
             'segmentSize=i' => \$segmentSize,
             'overlapSize=i' => \$overlapSize,
             'partition_listing=s' => \$partition_listing,
             'minSegmentSize=i' => \$minSegmentSize, ## dgg added
             'skipSegmentList=s' => \$skipSegmentList,
             'keepSegmentList=s' => \$keepSegmentList,
             'combineSmallSegments!' => \$combineSmallSegments, #dgg
             'h' => \$help,
             );

push(@GFF, @ARGV); # any remainder;

# process required options
if ($help) { die $param_string; }

#unless ($optok && $genome && @GFF>0 && $segmentSize && defined($overlapSize) && $partition_listing) {
unless ($optok && $genome && $segmentSize && defined($overlapSize) && $partition_listing) {
    die $param_string;
}

my $genome_basename = basename($genome);
my @files_to_process = ();
 
foreach my $gf (@GFF) {
  unless(-f $gf) { die "missing GFF=$gf\n"; next; }
  push(@files_to_process, 
    {  type => "GFF", # nice idea but not used
       basename => basename($gf),
       file => $gf,
   } ) ;
}


my %skipseglist=();
my %keepseglist=();
my %seqseglist=();
my %gffseglist=();
if($skipSegmentList) { #dgg
  open(F,$skipSegmentList) or die "Cannot process skiplist: $skipSegmentList";
  while(<F>){ chomp; $skipseglist{$_}=1; }
  close(F);
}
if($keepSegmentList) { #dgg
  open(F,$keepSegmentList) or die "Cannot process skiplist: $keepSegmentList";
  while(<F>){ chomp; $keepseglist{$_}=1; }
  close(F);
}

my $curr_dir = cwd();
my $util_dir = $FindBin::Bin;
my $gff_partitioner = "$util_dir/gff_range_retriever.pl";
my ($combo,$combosize,$comboacc)=(1,0,"");
my %did_part= ();

&partition_files_based_on_contig(@files_to_process);


open (my $ofh, ">$partition_listing") or die "Error, cannot write $partition_listing file";

my $fasta_reader = new Fasta_reader($genome);
while (my $seq_obj = $fasta_reader->next()) {
    
    my $accession = $seq_obj->get_accession();
    my $sequence = $seq_obj->get_sequence();
    my $sequence_length = length($sequence);
    
    ## create accession partition area
    my $acc_dir = "$curr_dir/$accession";

    my $toosmall= ($minSegmentSize && $sequence_length < $minSegmentSize);
    my $docombine=($combineSmallSegments && $toosmall && ! $skipseglist{$accession} );
    if($keepseglist{$accession}) {
      
    } elsif(scalar(%keepseglist)) {
      $toosmall=1; $docombine=0;
    }
    
    if($docombine) {
      $toosmall= 0;
      if($combosize + $sequence_length > $segmentSize) {
        $combo++; #?? only add new one if have reached segmentSize
        $combosize= 0;
      }
      $combosize= $combosize + $sequence_length;
      
      $comboacc="combo_".$combo;
      my $comboacc_dir = "$curr_dir/$comboacc";
      unless (-d $comboacc_dir) {
        mkdir ($comboacc_dir) or die "Error, cannot mkdir $comboacc_dir";
      }
      
      my $combogeno_filename="$comboacc_dir/$genome_basename";
      unless (-e $combogeno_filename) {
          system "touch $combogeno_filename";
        }
  
      if (-d $acc_dir) { 
        # copy to combo, then delete gff files in $acc_dir 
        foreach my $file_struct (@files_to_process) {
            my $basename = $file_struct->{basename};
            my $acc_filename = "$acc_dir/$basename";
            my $comboacc_filename = "$comboacc_dir/$basename";
            unless (-e $comboacc_filename) {
              system "touch $comboacc_filename";
            }
            open(IN,"<$acc_filename") or die "Error $acc_filename";
            open(OUT,">>$comboacc_filename") or die "Error $comboacc_filename";
            while(<IN>) { print OUT $_; }
            close(IN); close(OUT);
            unlink($acc_filename);
            }
        rmdir($acc_dir);
        }
      print "copied small $accession to $comboacc.\n";
      $acc_dir= $comboacc_dir;
    }
    
# dgg: opt to skip short or unwanted seqs 
    if($skipseglist{$accession} || $toosmall) {
      if (-d $acc_dir) { 
        # delete gff files in $acc_dir and move on..
        foreach my $file_struct (@files_to_process) {
            my $basename = $file_struct->{basename};
            my $acc_filename = "$acc_dir/$basename";
            unlink($acc_filename);
            }
        rmdir($acc_dir);
        }
      print "skipped $accession.\n";
      next;
    }

    unless (-d $acc_dir) {
        mkdir ($acc_dir) or die "Error, cannot mkdir $accession";
    }
    
    ## write the genome sequence
    my $fcreate= ($docombine) ? ">>" : ">";
    open (my $acc_ofh, "$fcreate$acc_dir/$genome_basename") or die $!;
    
    #print $acc_ofh $seq_obj->get_FASTA_format();
    ## fixme for augustus that mangles fasta reading: all header instead of >ID
    { my $seqfmt= $sequence; $seqfmt =~ s/(\S{60})/$1\n/g;
      print $acc_ofh ">$accession\n"; 
      print $acc_ofh $seqfmt,"\n";
    }
    close $acc_ofh;
    $seqseglist{$accession}++;
 
    ## pull out the relevant data:
    foreach my $file_struct (@files_to_process) {
        my $filename = $file_struct->{file};
        my $basename = $file_struct->{basename};
        my $acc_filename = "$acc_dir/$basename";
        $file_struct->{acc_file} = $acc_filename;

        unless (-e $acc_filename) {
            system "touch $acc_filename";
            if ($?) { die "Error, died running cmd: touch $acc_filename"; }
        }
        
        #my $cmd = "$gff_partitioner $accession 1 $sequence_length < $filename > $acc_filename";
        #&process_cmd($cmd);
    }
    
    my @ranges = &get_range_list($sequence_length);
    my $num_ranges = scalar (@ranges);
    print "$accession has $num_ranges partitions.\n";
    
    ## for $docombine, need to write $comboacc to $ofh/partlist >> only once!
    my $part_accession= ($docombine) ? $comboacc : $accession;
    
    if ($num_ranges == 1) {  ## ($docombine) should always be here; err if not?
        # do not partition further. (as indicated below with 'N')
        print $ofh "$part_accession\t$acc_dir\tN\n" 
          unless($did_part{$part_accession}++);
    }
    else {
        ## partition the data into segments:
        foreach my $range (@ranges) {
            my ($range_lend, $range_rend) = @$range;
            my $range_length = $range_rend - $range_lend + 1;
            
            my $partition_dir = "$acc_dir/$accession" . "_$range_lend" . "-$range_rend";
            mkdir $partition_dir or die "Error, cannot mkdir $partition_dir";
            my $genome_subseq = substr($sequence, $range_lend - 1, $range_length);
            $genome_subseq =~ s/(\S{60})/$1\n/g; #make fasta format
            
            open (my $part_ofh, ">$partition_dir/$genome_basename") or die "Error, cannot write to $partition_dir/$genome_basename";
            print $part_ofh ">$accession\n";
            print $part_ofh $genome_subseq;
            close $part_ofh;
            
            foreach my $file_struct (@files_to_process) {
                my $filename = $file_struct->{acc_file};
                my $basename = $file_struct->{basename};
                gff_partitioner( $filename, "$partition_dir/$basename",
                                  $accession, $range_lend, $range_rend, "ADJUST_TO_ONE");
                # my $cmd = "$gff_partitioner $accession $range_lend $range_rend ADJUST_TO_ONE < $filename > $partition_dir/$basename";
                # &process_cmd($cmd);
            }
            print $ofh "$part_accession\t$acc_dir\tY\t$partition_dir\n";
        }
    }
}

close $ofh;

clear_gffonly_dirs();

exit(0);



####

sub clear_gffonly_dirs {

  foreach my $acc (sort keys %gffseglist) {
    next if($seqseglist{$acc}); #? or $skipseglist{$acc} 
    my $acc_dir = "$curr_dir/$acc";
    if (-d $acc_dir) {
        # delete gff files in $acc_dir and move on..
        foreach my $file_struct (@files_to_process) {
            my $basename = $file_struct->{basename};
            my $acc_filename = "$acc_dir/$basename";
            unlink($acc_filename);
            }
        rmdir($acc_dir);
        }
      print "skipped gff-only $acc.\n";
   }
}

sub get_range_list {
    my ($sequence_length) = @_;

    my $range_lend = 1;
    my $range_rend = $segmentSize;
    my @ranges;
    
    while ($range_lend < $sequence_length - $overlapSize) {
        $range_rend = $range_lend + $segmentSize - 1;
        if ($range_rend > $sequence_length) {
            $range_rend = $sequence_length;
        }
        push (@ranges, [$range_lend, $range_rend]);

        $range_lend += $segmentSize - $overlapSize;
    }

    unless (@ranges) {
        ## must be a shorty such that the overlap length is greater than the sequence length.
        if ($sequence_length > $segmentSize) {
            die "Error, no ranges and sequence length ($sequence_length) does not exceed segmentSize ($segmentSize)";
        }
        push (@ranges, [1, $sequence_length]);
    }
    
    return (@ranges);
}



 

####
sub process_cmd {
    my $cmd = shift;
    print "CMD: $cmd\n";
    die "error, $cmd, $? " if system $cmd;
}

####
sub partition_files_based_on_contig {
    my @file_structs = @_;

    foreach my $file_struct (@file_structs) {

        my $file = $file_struct->{file};
        my $basename = $file_struct->{basename};
        
        my $ofh;
        my $curr_contig = "";
        open (my $fh, $file) or die "Error, cannot open $file";
        my ($newfile,$topline); # dgg
        while (<$fh>) {
            my $line = $_;
	    #unless (/\w/) { next; }
            unless (/^\w/) { $topline=$line unless($topline); next; } # dgg: drop gff comments; ? cant loose ##gff-version
            my @x = split (/\t/);
            my $contig_id = $x[0];
            if ($contig_id ne $curr_contig) {
                $curr_contig = $contig_id; $newfile=1;
                unless (-d "$curr_contig") {
                    #dgg-off# print STDERR "writing $file for $curr_contig\n";
                    mkdir $curr_contig or die "Error, cannot mkdir $curr_contig";
                }
                close $ofh if $ofh;
                open ($ofh, ">>$curr_contig/$basename") or die "Error, cannot write to $curr_contig/$basename";
	      $gffseglist{$curr_contig}++;
            }
            if($newfile && $topline) { print $ofh $topline; $newfile=0; }
            print $ofh $line;
        }
    }
    return;
}


# gff_range_retriever.pl
#  $gff_partitioner $accession $range_lend $range_rend ADJUST_TO_ONE < $filename > $partition_dir/$basename 
#  gff_partitioner( $ingff,$outgff, $accession $range_lend $range_rend ADJUST_TO_ONE);

sub gff_partitioner {
  my($ingff,$outgff,$seq_id,$min_lend,$max_rend,$adjust_to_1)= @_;

#   my $seq_id = $ARGV[0] or die $usage;
#   my $min_lend = $ARGV[1] or die $usage;
#   my $max_rend = $ARGV[2] or die $usage;
#   my $adjust_to_1 = $ARGV[3] || 0; # default is false
  $adjust_to_1 ||= 0;
  
  my $adjust_coord = $min_lend - 1;
  
  unless ($min_lend =~ /^\d+$/ && $max_rend =~ /^\d+$/) {
     warn "bad min_lend,max: $min_lend, $max_rend for $ingff\n"; return; # die $usage;
  }
  
  open(my $inh,"$ingff") or do { warn "cant read $ingff\n"; return; };
  open(my $outh,">$outgff") or do { warn "cant write $outgff\n"; return; };
  
  my $got_spacer = 0;
  while (<$inh>) {
      unless (/\w/) { 
          print $outh unless ($got_spacer);
          $got_spacer = 1;
          next;
      }
  
      if (/^\#/) {
          print $outh;
          next;
      }
  
      my @x = split (/\t/);
      my ($contig, $lend, $rend) = ($x[0], $x[3], $x[4]);
      
      if ($contig eq $seq_id &&
          $lend >= $min_lend &&
          $rend <= $max_rend) {
          
          if ($adjust_to_1) {
              $x[3] -= $adjust_coord;
              $x[4] -= $adjust_coord;
          }
          
          print $outh join ("\t", @x);
          $got_spacer = 0;
      }
  }
  close($outh); close($inh);
}
