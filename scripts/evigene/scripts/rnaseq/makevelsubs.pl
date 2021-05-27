#!/usr/bin/env perl

my $debug=1;
my $vprefix="sub1."; # subperf7
my $workd="/export/udisk2/work/aphid";
my $workd="/bio/bio-grid/nasv4";
my $rnas="/bio/bio-grid/mb/evigene/scripts/rnaseq";
my $rrnafilt="$workd/misc/repeat_rrna.gff";

my $samfilter = "$rnas/samfiltergff.pl -expandover=50 -over $rrnafilt ";
# my $sam2velv  = "$rnas/sam2velv.pl -debug -mismatch=0.05 -endclip=0.06 -hardclip -samout";
my $sam2velv  = "$rnas/sam2velv.pl -debug -noduplic -mismatch=0.15 -endclip=0.66 -hardclip ";
my @outsuf= qw( fa fa1 fa2 ); # sam

my @sublist= grep /\.list/, @ARGV;
my @bamlist= grep /\.bam$/, @ARGV;

foreach my $subf (@sublist) {
  open(L,$subf) or die "open $subf"; 
    warn("# make list $subf\n") if $debug;
  while(<L>) {
    my ($tp,$na,$sublist)=split; 
    $sublist=~s/,/ /g;

    my $subname="$vprefix$na"; 
    warn("# write subset $sugname\n") if $debug;
    foreach my $suf (@outsuf) { system("touch $subname.$suf"); }

    my $ib=0;
    foreach my $bam (@bamlist) {
      warn("# read $bam\n") if $debug;
      my $tmpname="sub$na.tmp$ib"; $ib++;

if(1) {
      my $err;
      my $tmpsam1="$tmpname.samp1";
      my $tmpsam2="$tmpname.samp2";
      my $cmd1="samtools view $bam $sublist > $tmpsam1";
      my $cmd2="$samfilter < $tmpsam1 > $tmpsam2";
      my $cmd3="$sam2velv -name $tmpname < $tmpsam2";

      warn "cmd1: $cmd1\n" if $debug; $err= system($cmd1);
      if($err) { warn "ERR: cmd1 failed: $?\n"; }
      if(-s $tmpsam1) { warn "cmd2: $cmd2\n" if $debug; $err= system($cmd2);
        if($err) {  warn "ERR: cmd2 failed: $?\n"; }
        if(-s $tmpsam2) { warn "cmd3: $cmd3\n" if $debug; $err= system($cmd3);
          if($err) {  warn "ERR: cmd3 failed: $?\n"; }
          }
        }
      unlink($tmpsam1) if(-f $tmpsam1);
      unlink($tmpsam2) if(-f $tmpsam2);

} else {
	## this pipe is mem intens; split to parts w/ tmp .tmpsam files
      my $cmd="samtools view $bam $sublist | $samfilter | $sam2velv -name $tmpname";
      warn("$cmd\n") if $debug;
      system($cmd);   # makes new subname.{sam,fa,fa1,fa2}; 
}

      foreach my $suf (@outsuf) { 
        my $subt="$tmpname.$suf";
        if( -f $subt ) { system("cat $subt >> $subname.$suf"); unlink($subt); }
	}
    } close(SUBSAM);

  } close(L);
}

