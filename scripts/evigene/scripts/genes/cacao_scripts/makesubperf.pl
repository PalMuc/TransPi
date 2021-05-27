#!/usr/bin/env perl

my $debug=1;
my $vprefix="sub1."; 
my $workd="/bio/bio-grid/cacao3";
my $rnas="/bio/bio-grid/mb/evigene/scripts/rnaseq";

# cacao2 skip this: only rdna_scaffold_xxx
# my $rrnafilt="$workd/misc/repeat_rrna.gff";

my $samtools="samtools";

## FIXME: what is -noduplic effect on velv/oases assemblies, esp for alt-tr sharing exons? redo as -noduplic=8 ?
# my $samfilter = "$rnas/samfiltergff.pl -expandover=50 -over $rrnafilt ";
my $samfilter = "cat";

# my $sam2velv  = "$rnas/sam2velv.pl -debug -mismatch=0.05 -endclip=0.06 -hardclip -samout";
my $sam2velv  = "$rnas/sam2velv.pl -debug -noduplic=5 -mismatch=0.05 -endclip=0.25 -hardclip ";
my @outsuf= qw( fa fa1 fa2 ); # sam

my @sublist= grep /\.list/, @ARGV;
my @bamlist= grep /\.bam$/, @ARGV;

foreach my $subf (@sublist) {
  open(L,$subf) or die "open $subf"; 
    warn("# make list $subf\n") if $debug;
  while(<L>) {
    my ($tp,$na,$sublist)=split;  next unless($tp and $na); # dang got blanks
    $sublist=~s/,/ /g;

    my $subname="$vprefix$na"; 
    warn("# write subset $sugname\n") if $debug;
    foreach my $suf (@outsuf) { system("touch $subname.$suf"); }

    my $ib=0;
    foreach my $bam (@bamlist) {
      warn("# read $bam\n") if $debug;
      my $tmpname="sub$na.tmp$ib"; $ib++;

     {
      my $err;
      my $tmpsam1="$tmpname.samp1";
      my $tmpsam2=$tmpsam1; #OFF# "$tmpname.samp2";
      my $cmd1="$samtools view $bam $sublist > $tmpsam1";
      my $cmd2="$samfilter < $tmpsam1 > $tmpsam2";
      my $cmd3="$sam2velv -name $tmpname < $tmpsam2";

      warn "cmd1: $cmd1\n" if $debug; $err= system($cmd1);
      if($err) { warn "ERR: cmd1 failed: $?\n"; }
      if(-s $tmpsam1) { 
         #OFF# warn "cmd2: $cmd2\n" if $debug; $err= system($cmd2);
         #OFF# if($err) {  warn "ERR: cmd2 failed: $?\n"; }
        if(-s $tmpsam2) { warn "cmd3: $cmd3\n" if $debug; $err= system($cmd3);
          if($err) {  warn "ERR: cmd3 failed: $?\n"; }
          }
        }
      unlink($tmpsam1) if(-f $tmpsam1);
      unlink($tmpsam2) if(-f $tmpsam2);
    } 

      foreach my $suf (@outsuf) { 
        my $subt="$tmpname.$suf";
        if( -f $subt ) { system("cat $subt >> $subname.$suf"); unlink($subt); }
	}
    } close(SUBSAM);

  } close(L);
}

