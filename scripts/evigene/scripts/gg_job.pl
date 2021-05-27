#!/usr/bin/env perl
# gg_job.pl

=head1 TITLE

gg_job.pl : run genome grid jobs on parallelized data paritions

=head1 SYNOPSIS

  genogrid/scripts/partition_GFF_genome.pl ... --partitions partitions.list
    : prepare genome data partitions

  genogrid/scripts/gg_job.pl -debug \
    --jobtype script --script genogrid/scripts/blast2exonrfine.sh \
    --QTYPE PBS --CLASS dque --CLOCKTIME "8:55:00" --ACCOUNT "TG-xxx" \
    --do_parts_batch 10 \
    --bin_dir ~/bio/exonerate/bin/  \
    --genome mygenome.fa \
    --proteins `pwd`/yourgenes.aa \
    --output mygenome-yourgenes.gff \
    --partitions partitions.list

  genogrid/scripts/recombine_partial_outputs.pl .. --part partitions.list  
    : collate partitioned results
    
    
=head1 ABOUT

This is a tool [gg_job.pl] to run an analysis script for each data part.  The basic analysis is one
that can be run in any given subdirectory as a shell script fed into the batch queue.  
PBS and LoadLeveler queue types are now supported.   There are a few basic genome analyses
internally supported (Blast, exonerate, augustus, ..) and the facility to turn a shell script
template into scripts for each data partition, with variable substitution (e.g. genome file
names).

Typically several different analyses are run on the same data partition tree.  Some are slower
(e.g. 10 hr/part) and some are faster (e.g. 1 minute/part).  To use the same data partition
effectively with this, an option to batch jobs is included (gg_job -do_parts_batch 10). THis
will write scripts for 10 parts into one job file for the batch queue.

=head1 AUTHOR

Don Gilbert, gilbertd@indiana.edu

=head1 NOTES

=cut

use strict;
use warnings;
use Getopt::Long ;
use File::Basename;

my $ACCOUNT="none";
my $CLOCKTIME="02:00:00";
my $CLASS="NORMAL";
my $QTYPE=""; # one of supported: PBS, LL, later: globusrun/gram
my @QTYPES= qw(LL PBS); # available types


my $snap_hmm="minimal.hmm"; #"dgri_caf1_DGIL_SNO.hmm";
my( $base_dir, $oldbase, $bin_dir, $idname, $proteins, $jobtype,
  $blastdata, $blastprog, $blastopts,
  $species, $hintsfile, $jobbatch,
  $genome, $output_file_name, $partitions, $part_dir, $debug,
  ) = ("") x 40;

## TODO: allow option to pack several quick part_dir scripts into one?
## e.g. 200 x 1 minute exonerate calls take longer to queue up and admin than run

my $parts_per_batch = 10; 
my $do_parts_batch  = 0;

$jobtype="";
$output_file_name="ggjob"; # add .gff, .aa, .tr
$partitions="partitions.list";


my $optok= &GetOptions (
  "jobtype=s"=>\$jobtype, # snap, exonerate, blast1, blastpart
  "hmm=s"=>\$snap_hmm,
  "idname=s"=>\$idname,
  "genome=s"=>\$genome,
  "proteins|blastdata=s"=>\$proteins,
  "hintsfile=s"=>\$hintsfile, # also extrinsicCfgFile
  "species=s"=>\$species,
  "blastprog|script=s"=>\$blastprog, # use jobtype ??
  "blastopt|options=s"=>\$blastopts,
  "bin_dir=s"=>\$bin_dir,
  "partitions=s"=>\$partitions,
  "do_parts_batch=s"=>\$do_parts_batch,
  "output_file_name=s"=>\$output_file_name,
  "ACCOUNT=s"=>\$ACCOUNT,
  "CLOCKTIME=s"=>\$CLOCKTIME,
  "CLASS=s"=>\$CLASS, # same as Q-name
  "QTYPE=s"=>\$QTYPE, # snap, exonerate, blast1, blastpart
  "base_dir=s"=>\$base_dir, # drop this
  "oldbase=s"=>\$oldbase,   # drop this
  "n|debug!"=>\$debug,
  );

$blastdata=$proteins;
#$scriptfile=$blastprog;

unless($idname){ $idname=$snap_hmm; $idname =~ s/\.[^\.]*$//; $idname =~ s/\W/_/g; }

my $USAGE= <<"USAGE";

:: gg_job : genome grid job ::
usage: gg_job \\
        --jobtype=augustus|exonerate|blast1|blastpart|snap|scriptfile \\
        --genome=$genome \\
        --partitions=partitions.list \\
        --bin_dir=$bin_dir   : /path/to/binary \\
        --output=$output_file_name \\
options:
        --hmm=$snap_hmm      : snap hmm *
        --proteins=modDM.fa  : exonerate input proteins *
        --species=species_name  : augustus config species* 
        --hintsfile=mygenome.hints  : augustus hints
        --blastdata=modDM.fa : blast1: blastdb; blastpart: part/fasta  *
        --blastprog=tblastn  : blast: blastn,blastx,tblastn,tblastx   *
        --script=xyz.shell   : shell script with variables part_dir, genome, output_file_name
        --do_parts_batch=n   : collect n genome parts into one batch script
        --idname=$idname  
        --n|debug            : dont run, but write grid scripts
grid admin:
        --QTYPE=$QTYPE       : queue system, one of supported: @QTYPES 
        --CLASS=$CLASS       : queue name
        --ACCOUNT=$ACCOUNT   --CLOCKTIME=$CLOCKTIME
USAGE

### need these only if partion on one path, move to another path     
##   --base_dir=/path/to/genome --oldbase=/old/path

die "missing --bin_dir=$bin_dir".$USAGE unless (-d $bin_dir);
$partitions= "$base_dir/partitions" unless (-f $partitions);
die "missing --partitions=$partitions".$USAGE unless (-f $partitions);

die $USAGE if($jobtype =~ /augustus/ and not ($species));
die $USAGE if($jobtype =~ /blast/ and not ($blastdata && $blastprog));
die $USAGE if($jobtype =~ /exonerate/ and not $proteins);
die $USAGE if($jobtype =~ /snap/ and not $snap_hmm);

die $USAGE unless($optok and $jobtype and $genome and $output_file_name);
die $USAGE unless(grep{$QTYPE eq $_} @QTYPES);
# want to check ACCOUNT, CLASS, ??
die $USAGE unless($ACCOUNT and $CLASS);

$parts_per_batch= $do_parts_batch if($do_parts_batch>1);

my $scriptstring="";
if($jobtype =~ /script/) {
  my $ok;
  $ok= open(F,$blastprog);
  $ok= open(F,"$bin_dir/$blastprog") unless $ok;
  die "cant open script $blastprog".$USAGE unless $ok;
  $scriptstring = join("", <F>); close(F);
  die "missing script $blastprog" unless($scriptstring =~ /\S+/);
}

 
$bin_dir =~ s,/$,,;
my $partn=0;

open (my $fh, $partitions) or die "Error, cannot open $partitions";
while (<$fh>) {
    next unless(/^\w/);
    chomp;
    my ($accession, $acc_dir, $is_partitioned, $partition_dir) = split (/\t/);
    next unless($acc_dir);
    $partn++;
    
    $part_dir = $acc_dir;
    if ($is_partitioned eq "Y") {
        $part_dir = $partition_dir;
    }
    
    $part_dir =~ s,$oldbase,$base_dir, if($oldbase); # die if fail?
    $part_dir =~ s,/$,,;
    
    my $jobna= basename("$output_file_name.$partn");
    
    unless (-d $part_dir) {
      die "Error, cannot locate data directory $part_dir \n";
      
    } elsif($jobtype =~ /script/) {
      # replace $genome, $output_file_name, part_dir in scriptfile ...
      # need better way to substitute variables script may want.
      submit_script_vars($jobna, $scriptstring, $part_dir);

    } elsif($jobtype =~ /exonerate/) {
      submit_exonerate($jobna, $proteins, $part_dir);

    } elsif($jobtype =~ /augustus/) {
      submit_augustus($jobna, $species, $hintsfile, $part_dir);

    } elsif($jobtype =~ /blast1/) {
      $blastopts= "-m 9 -e 0.001" unless($blastopts);
      submit_blast1($jobna, $blastdata, $blastprog, $blastopts, $part_dir);
      
    } elsif($jobtype =~ /blastpart/) {
      $blastopts= "-m 9 -e 0.001" unless($blastopts);
      submit_blastpart($jobna, $blastdata, $blastprog, $blastopts, $part_dir);
      
    } elsif($jobtype =~ /snap/) {
      submit_snap($jobna, $snap_hmm, $idname, $part_dir);
    
    } else {
      die "Cant handle jobtype: $jobtype".$USAGE;
    }
    # die if not submit  ?
 }

my $jobna= basename("$output_file_name.$partn");
submit_part_collection( 0, $jobna, "", $part_dir) # flush last batch if any
  if ($do_parts_batch); # and $jobbatch
  
# exit(0);

#............ genome apps ....................

sub submit_snap {
  my($jobn, $snap_hmm, $idname,  $part_dir)= @_;
  # my $jobn="snaprun.$partn";

  my $job = <<"EOF";
export ZOE=${bin_dir}

${bin_dir}/snap -quiet -name $idname -gff \\
  -aa ${part_dir}/$output_file_name.aa \\
  -tx ${part_dir}/$output_file_name.tr \\
  ${snap_hmm} ${part_dir}/$genome \\
  > ${part_dir}/$output_file_name.gff

EOF
  
  return submit_script($jobn, $job, $part_dir);
}


# blastn, blastx, tblastn some query .fa against genome partions dna
# which is .db, which query? 
# either formatdb each part/genome.fa or part/genome is query
#  tblastn -i prots -d genome should be ~ blastx -i genome -d prots
# case 1: one blastdb w/ fixed path, part/genome = query, pre-formatdb 
# case 2: part/query x part/genome : (pre?) formatdb each part/genome

sub submit_blast1 {
  my($jobn, $blastdb, $blastprog, $blastopts,  $part_dir)= @_;

  unless($blastprog =~ /blastn|blastx|tblastn|tblastx/) { #psitblastn,...?
    die "missing blastprog: $blastprog\n"; return -1;
    }
  unless(-f "$blastdb.nsq" or -f "$blastdb.psq") {
    die "missing blastdb: $blastdb.[np]sq\n"; return -1;
    }

  my $job = <<"EOF";
${bin_dir}/blastall -p $blastprog  \\
  $blastopts \\
  -d $blastdb \\
  -i ${part_dir}/$genome \\
  -o ${part_dir}/$output_file_name.$blastprog
EOF
  
  return submit_script($jobn, $job, $part_dir);
}

sub submit_blastpart {
  my($jobn, $blastquery, $blastprog, $blastopts,  $part_dir)= @_;

  unless($blastprog =~ /blastn|blastx|tblastn|tblastx/) { #psitblastn,...?
    die "missing blastprog: $blastprog\n"; return -1;
    }
  unless(-f "${part_dir}/$blastquery") {
    warn "missing query:${part_dir}/$blastquery\n"; return -1;
    }
    
  my $job = <<"EOF";
${bin_dir}/formatdb -p F -i ${part_dir}/$genome 

${bin_dir}/blastall -p $blastprog  \\
  $blastopts \\
  -i ${part_dir}/$blastquery \\
  -d ${part_dir}/$genome \\
  -o ${part_dir}/$output_file_name.$blastprog
EOF
  
  return submit_script($jobn, $job, $part_dir);
}


sub submit_exonerate {
  my($jobn, $proteins, $part_dir)= @_;

  ## reverse this: pick partdir/prots first
  # my $protpath= (-f $proteins) ? $proteins : "${part_dir}/$proteins" ;
  my $protpath= (-f "${part_dir}/$proteins") ? "${part_dir}/$proteins" : $proteins;
  unless(-f $protpath) { warn "missing proteins: $protpath\n"; return -1; }
  
  my $job = <<"EOF";
${bin_dir}/exonerate --model protein2genome --minintron 20 --maxintron 5000 \\
  --showtargetgff --showvulgar 0 --showalignment 0 \\
  --ryo '#qi %qi length=%ql alnlen=%qal\\n#ti %ti length=%tl alnlen=%tal\\n' \\
  --query $protpath\\
  --target ${part_dir}/$genome \\
  > ${part_dir}/$output_file_name.gff
EOF
  
  return submit_script($jobn, $job, $part_dir);
}

sub submit_augustus {
  my($jobn, $species, $hintsfile, $part_dir)= @_;

#  my $protpath= (-f $proteins) ? $proteins : "${part_dir}/$proteins" ;
#  unless(-f $protpath) { warn "missing proteins: $protpath\n"; return -1; }

  my $aug_dir=$bin_dir; #??
  my $aug_bin="augustus"; # fixme: param?
  my $extrinsicCfgFile="$aug_dir/config/augrun.cfg"; ## FIXME
  
  my $addhints= ($hintsfile) ?
    " --hintsfile=${part_dir}/$hintsfile --extrinsicCfgFile=$extrinsicCfgFile "
    : "";
  
  my $job = <<"EOF";
${aug_dir}/bin/${aug_bin} --gff3=on --uniqueGeneId=true \\
  --species=$species $addhints \\
  --AUGUSTUS_CONFIG_PATH=$aug_dir/config/ \\
  ${part_dir}/$genome  \\
  > ${part_dir}/$output_file_name.gff

EOF
  
  return submit_script($jobn, $job, $part_dir);
}

=head2 script vars?

 substitute any of these as \$varname in script?
 
?  "hmm=s"=>\$snap_hmm,
*  "idname=s"=>\$idname,
x  "genome=s"=>\$genome,
x  "proteins|blastdata=s"=>\$proteins,
?  "hintsfile=s"=>\$hintsfile, # also extrinsicCfgFile
*  "species=s"=>\$species,
?  "blastprog|script=s"=>\$blastprog, # use jobtype ??
?  "blastopt|options=s"=>\$blastopts,
*  "bin_dir=s"=>\$bin_dir,
  "partitions=s"=>\$partitions,
x  "output_file_name=s"=>\$output_file_name,

=cut

sub replace_script_vars {
  my($jobn, $jobscript, $part_dir)= @_;
  # need better way to substitute variables script may want.

  $jobscript =~ s/\$options/$blastopts/g; #? not always set; alias blastopts = options
  $jobscript =~ s/\$proteins/$proteins/g;      # == $blastdata; may be set 
  $jobscript =~ s/\$idname/$idname/g;     # always set to a default

  $jobscript =~ s/\$bin_dir/$bin_dir/g;   # always set
  $jobscript =~ s/\$part_dir/$part_dir/g;  # always set
  $jobscript =~ s/\$genome/$genome/g;      # should be set
  $jobscript =~ s/\$output_file_name/$output_file_name/g;  # alway set
  
  return $jobscript;
}

sub submit_script_vars {
  my($jobn, $jobscript, $part_dir)= @_;

  # need better way to substitute variables script may want.
  $jobscript= replace_script_vars($jobn, $jobscript, $part_dir);

  return submit_script($jobn, $jobscript, $part_dir);
}

#....... subs .........................

# generalize this in package? or instead use ergatis, other workflow tool?
# submit_job_augustus_llsubmit
# submit_job_snap_gridrun

      # replace $genome, $output_file_name, part_dir in scriptfile ...
sub submit_script_ll {
  my($jobn, $jobscript, $part_dir)= @_;
  
  unless($jobscript =~ m/\#\s*\@\s*account_no\s*=/) {
  my $jobhead = <<"EOS";
#! /bin/bash -l
### llsubmit $jobn: IBM LoadLeveler
# @ class = $CLASS
# @ account_no = $ACCOUNT
# @ wall_clock_limit = $CLOCKTIME
# @ error   = $part_dir/$jobn.err
# @ notification = always
# @ environment=COPY_ALL;
# @ queue

EOS
  $jobscript = $jobhead . $jobscript;
  }

#  ## but see submit_script_vars ; should we do subst on each?
#   $jobscript =~ s/\$part_dir/$part_dir/g;
#   $jobscript =~ s/\$genome/$genome/g; # longer?
#   $jobscript =~ s/\$output_file_name/$output_file_name/g;
  
  open(J,">$jobn") or die "$jobn"; print J $jobscript; close(J);
  if($debug) {
    print "llsubmit -q $jobn\n";
  } else {
    my $res=`llsubmit -q $jobn`; chomp($res);
    print "# llsubmit $jobn : $res\n";
  }  
  return 1;
}

sub submit_script_pbs {
  my($jobn, $jobscript, $part_dir)= @_;

  unless($jobscript =~ m/\#PBS\s-A\s*\S+/) {
  my $jobhead = <<"EOS";
#! /bin/bash -l
### qsub -q $CLASS $jobn
#PBS -N $jobn
#PBS -A $ACCOUNT
#PBS -l nodes=1:ppn=1,walltime=$CLOCKTIME
#PBS -o $part_dir/$jobn.out
#PBS -e $part_dir/$jobn.err
#PBS -V

EOS
#xxx PBS -m abe
#xxx PBS -M gilbertd@indiana.edu

  $jobscript = $jobhead . $jobscript;
  }

#   $jobscript =~ s/\$part_dir/$part_dir/g;
#   $jobscript =~ s/\$genome/$genome/g; # longer?
#   $jobscript =~ s/\$output_file_name/$output_file_name/g;
  
  open(J,">$jobn") or die "$jobn"; print J $jobscript; close(J);
  if($debug) {
    print "qsub -q $CLASS $jobn\n";
  } else {
    my $res=`qsub -q $CLASS $jobn`; chomp($res);
    print "# qsub -q $CLASS $jobn : $res\n";
  }  
  return 1;
}

## TODO: allow option to pack several quick part_dir scripts into one?
## e.g. 200 x 1 minute exonerate calls take longer to queue up and admin than run

# my ($jobbatch);
# my $parts_per_batch = 10; 
# my $do_parts_batch= 0;

sub submit_part_collection {
  my($partn, $jobn, $jobscript, $part_dir)= @_;
  my $ok=1;
 
  # $jobscript= replace_script_vars($jobn, $jobscript, $part_dir);
  $jobbatch .= $jobscript . "\n";
  
  # $ipart++; == partn?
  if ( ($partn % $parts_per_batch) == 0 && $jobbatch =~ /\S/) { 
    $ok= submit_script1( $jobn, $jobbatch, $part_dir); 
    $jobbatch="";
  }
  
  return $ok;
}

sub submit_script1 {
  my($jobn, $jobscript, $part_dir)= @_;

  return 0 unless($jobscript =~ /\S/);
  return submit_script_pbs($jobn, $jobscript, $part_dir) if($QTYPE eq "PBS");
  return submit_script_ll($jobn, $jobscript, $part_dir)  if($QTYPE eq "LL");
  # else error/die 
  die "Cant handle $QTYPE"; #already checked
}

sub submit_script {
  # my($jobn, $jobscript, $part_dir)= @_;

  #? always or not? see  submit_script_vars
  # $jobscript= replace_script_vars($jobn, $jobscript, $part_dir);

  if($do_parts_batch) { return submit_part_collection($partn, @_); }
  else { return submit_script1(@_); }
  
#   return submit_script_pbs($jobn, $jobscript, $part_dir) if($QTYPE eq "PBS");
#   return submit_script_ll($jobn, $jobscript, $part_dir)  if($QTYPE eq "LL");
#   # else error/die 
#   die "Cant handle $QTYPE"; #already checked
}

1;

__END__

