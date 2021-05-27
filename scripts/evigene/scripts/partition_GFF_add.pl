#!/usr/bin/env perl
# partition_GFF_add.pl from partition_GFF_genome.pl
# d.gilbert: generalized from partition_EVM_inputs
# input: partition list + any GFF files

=item example use

 $aug/scripts/partition_GFF_add.pl -part partition.list  match-modDM.gff 

 $aug/scripts/partition_GFF_add.pl -genome dgri_caf060210.fa -part partition.list -type fagff match-modDM.gff modDM.fa
 $aug/scripts/partition_GFF_add.pl -part partition.list -type fagff match-ensAG.gff ensAG.fa

=item problems here w/ large part set

 this script died adding to acyr 17,000 parts; 
 was very slow; 
 also prints GLOB(X) x 10000 to err log; dont see why ...

  # instead use this (which also takes a while w/ 17k part folders)

  grep -v Y partitions.list | perl -ne\
 '($s,$p,$n)=split; $c="fgrep \"$s\t\" acyr-predictfix.gff > $p/acyr-predictfix.gff"; system($c);' \
  >& log.fg &

  plus this script for 'grep Y partitions.list > party.list'
  

=cut

use warnings;
use strict;
use FindBin;
use Getopt::Long; # qw(:config no_ignore_case bundling);
use Cwd;
use Carp;
use lib ("$FindBin::Bin/../PerlLib", "$FindBin::Bin"); # use Bin/ as alternate loc for Fasta_reader 
use Fasta_reader; # this is lightweight PASA/EvmUtils module; keep
use File::Basename;

umask(0000);

my $param_string = <<__PARAMS__;

###############################################################################################

The evidence for each contig is partitioned into a separate contig directory.
This script adds new GFF data to existing genome partitions. See partition_GFF_genome.pl


###############################################################################################
#
# Options:
#  --GFF                   * :files containing gene predictions + evidence; repeat as needed
#                              all extra arguments taken as GFF input data
#  --partitions_file       * :name of input file that contains the list of partitions
#  --genome                  :genome fasta for accession list in combo_ folders
#  --accession_list          :OR accession list in combo_ folders
#  -h                        : help message with more detailed information.
#
# Alt use: --type=fagff  input.gff input.fasta  : output fasta parts by ID in input.gff parts
#
#   (*required options)
#
################################################################################################

__PARAMS__

    ;

my ($genome, @GFF, $type,
    $partitions_file, $accession_list, $help,
    );

# $accession_list="accession.list";
$type= "GFF";

my $optok= &GetOptions (
             'GFF=s', \@GFF,
             'type=s', \$type,
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
    } else {
       #? $base_directories_noparts{$base_dir} = 1;
    }
}
close $fh;

my @files_to_process = ();
 
if( $type eq 'fagff' ) {
  my $ingff= shift @GFF;
  my $infasta= shift @GFF; # FIXME option
  my $fasta_list= get_fasta_list($infasta);
  push(@files_to_process, 
    {  type => "fagff", 
       basename => basename($infasta), # output part fasta, not GFF
       file => $ingff, # input file
       fasta_list => $fasta_list,
   } ) ;
  
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

# my $curr_dir = cwd();
# my $util_dir = $FindBin::Bin;

foreach my $base_dir (sort keys %base_directories_info) {
    my $partition_dirs_href= {};
    my $haspart= $base_directories_info{$base_dir}{haspart};
    if( $haspart ) {
       $partition_dirs_href = $base_directories_to_partitions{$base_dir};
    }
    
    my $accession = $base_directories_info{$base_dir}{accession};
        ## ^^ problem here 'combo_n' accessions; should/dont have accesssion list in part file
        
    unless (-d $base_dir) {
       die "Error, missing $accession folder: $base_dir";
    }
        
    foreach my $file_struct (@files_to_process) {
        my $in_filename = $file_struct->{file};
        my $basename = $file_struct->{basename};
        my $type = $file_struct->{type};
        my $in_fasta_list = $type eq 'fagff' ? $file_struct->{fasta_list} : undef;

      if (!$haspart) {  
        my $part_file = "$base_dir/$basename";
        ## problem 'combo_n' accession bad here
        ## need list in part file
        my %accs=();
        if($accession =~ /^combo_/) {
          my $lh;
          if( $accession_list && open( $lh, "$base_dir/$accession_list") ) {  
            while(<$lh>){  s/^>//; if(m/^(\S+)/){ $accs{$1}++;} } close($lh); $accession= \%accs;
          } elsif( $genome && open( $lh, "$base_dir/$genome") ) { 
            while(<$lh>){ if(m/^>(\S+)/){ $accs{$1}++;} } close($lh); $accession= \%accs;
          }
        }

        if($type =~ /^fagff/i) { # need both input gff, input fasta here
          fasta_from_idlist($in_filename, $in_fasta_list, $part_file,  $accession, 1, 999999999, 0);
        } elsif($type =~ /^gff/i) {
          gff_partitioner( $in_filename, $part_file, $accession, 1, 999999999, ""); # want maxint/undef range
        }                    
      }
      else {
        my $part=0;
        foreach my $partition_dir (sort keys %$partition_dirs_href) {
          $partition_dir =~ /(\w+)_(\d+)-(\d+)$/ or die "Error, cannot extract coords from $partition_dir";
          my ($segref, $range_lend, $range_rend) = ($1,$2,$3);
          $part++;
          
          my $part_file = $partition_dir . "/$basename";
          if($type =~ /^fagff/i) { # need both input gff, input fasta here
            fasta_from_idlist($in_filename, $in_fasta_list, $part_file,  $segref, $range_lend, $range_rend, $part);
          } elsif($type =~ /^gff/i) {
            gff_partitioner( $in_filename, $part_file, $segref, $range_lend, $range_rend, "ADJUST_TO_ONE");
          }
        }
     
      }
   
   }  # end files_to_process
    
  
} # end base_dir


exit(0);



####


=item fasta_from_idlist

cat ensag.idlist $bg/shared/prots/protfa2/ens{AG,AM}.fa | \
  perl -ne'if(/^>(\S+)/){$d=$1;$f=1;$d=~s/gnl.\w+.//; $d=~s/\-..$//; $p=($d{$d}) ?1:0;} \
  unless($f){ chomp; $d{$_}++;} print if($p);' > ensag.aa

=cut

sub get_fasta_list {
  my($filename)= @_;
  my %seqobjlist= ();
  my $fasta_reader = new Fasta_reader($filename);
  while (my $seq_obj = $fasta_reader->next()) {
      my $accession = $seq_obj->get_accession();
      # $accession =~ s/gnl\|\w+\|//; # drop ncbi prefix; fix in input data ! assume gff, fa IDs match
      $seqobjlist{$accession}= $seq_obj;
  }
  return \%seqobjlist;
}


sub fasta_from_idlist {
  my($ingff, $seqobjlist, $outfasta, $seq_id, $min_lend, $max_rend, $part)= @_;

  open(my $inh,"$ingff") or do { warn "cant read $ingff\n"; return; };
  # open(my $infa,"$infasta") or do { warn "cant read $infasta\n"; return; };
  open(my $outh,">$outfasta") or do { warn "cant write $outfasta\n"; return; };
  my %didid;
  
  while (<$inh>) {
      next unless (/^\w/);
      my @x = split (/\t/);
      my ($contig, $lend, $rend) = ($x[0], $x[3], $x[4]);
      
      my $keep=0;
      if(ref $seq_id) { $keep= ($seq_id->{$contig}) ? 1 : 0; }
      else { $keep= ($contig eq $seq_id &&
              $lend >= $min_lend && $rend <= $max_rend) ? 1 : 0;
            }   
            
      if ($keep) {
        my($id)= $x[8]=~ m/(?:ID|Parent)=([^;\s]+)/;
        next if(!$id || $didid{$id}++);
        my $seq_obj= $seqobjlist->{$id};
        next unless(ref $seq_obj);
        print $outh $seq_obj->get_FASTA_format();
      }
  }
  
  close $inh;  close $outh;
}


# gff_range_retriever.pl
#  $gff_partitioner $accession $range_lend $range_rend ADJUST_TO_ONE < $filename > $partition_dir/$basename 
#  gff_partitioner( $ingff,$outgff, $accession $range_lend $range_rend ADJUST_TO_ONE);

sub gff_partitioner {
  my($ingff,$outgff,$seq_id,$min_lend,$max_rend,$adjust_to_1)= @_;

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
          print $outh $_ unless ($got_spacer);
          $got_spacer = 1;
          next;
      }
  
      if (/^\#/) {
          print $outh $_;
          next;
      }
  
      my @x = split (/\t/);
      my ($contig, $lend, $rend) = ($x[0], $x[3], $x[4]);
      
      my $keep=0;
      if(ref $seq_id) { $keep= ($seq_id->{$contig}) ? 1 : 0; }
      else { $keep= ($contig eq $seq_id &&
              $lend >= $min_lend && $rend <= $max_rend) ? 1 : 0;
            }          
      if ($keep) {
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
