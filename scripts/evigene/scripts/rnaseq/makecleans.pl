#!/usr/bin/env perl
# make clean.sam, rrna.sam with hardclip for cufflinks v0.6
# samtools -F 0x104 = drop unmapped(4), duplicates(x100)
# ** revise DO NOT drop dupls, for cufflinks, need for paralog finding **

my $debug=1;
# $rnas="/bio/bio-grid/mb/rnaseq";
my @sublist= grep /\.list/, @ARGV;
my @bamlist= grep /\.bam$/, @ARGV;

foreach my $subf (@sublist) {
  open(L,$subf) or die "open $subf"; 
    warn("# make list $subf\n") if $debug;
  while(<L>) {
    my ($tp,$na,$sublist)=split; 
    $sublist=~s/,/ /g;
    my $subsam="subs.$na.sam.unsort"; 
    if( -f $subsam ) { warn "# exists $subsam; next\n"; next; }
    system("touch $subsam");
    $ns=0; $nclip=0;
    warn("# write subset $subsam\n") if $debug;
    open(SUBSAM,">>$subsam") or die "write $subsam";
    foreach my $bam (@bamlist) {
       warn("# read $bam\n") if $debug;
      open(SAM,"samtools view -F 0x4 $bam $sublist |") or die "samtools $bam";
      while(<SAM>) { $_= hardclip($_); print SUBSAM $_; $ns++; } close(SAM);
    } close(SUBSAM);
    warn("# wrote n=$ns, nclip=$nclip to subset $subsam\n") if $debug;
  } close(L);
}

## need to sort all subsam...

sub hardclip {
  my($linein)= @_;
  my($qid, $flag, $chr, $cloc, $mapq, $cigar, $matechr, $mateloc, $isize, $seq, $qual, @opt)
    = split"\t", $linein;
  if($cigar =~ m/\dS/) { 
      my ($oldcig,$cl,$cr)= ($cigar,0,0);
      if($cigar =~ s/(\d+)S$/$1H/) { $cr=$1; } # $seq=substr($seq,0,-$cr); 
      if($cigar =~ s/^(\d+)S/$1H/) { $cl=$1; } # $seq=substr($seq,$cl); 
      my @linein=split"\t",$linein;
      $linein[5]= $cigar;  $nclip++;
      if($cr or $cl) {
        if($seq ne "*") { 
         $seq=substr($seq,0,-$cr) if $cr;  $seq=substr($seq,$cl) if $cl; 
         $linein[9]= $seq;
         $len  = length($seq); 
         }
        if($qual ne "*") { 
         $qual=substr($qual,0,-$cr) if $cr; $qual=substr($qual,$cl) if $cl; 
         $linein[10]= $qual;
         }
       }
    $linein= join"\t", @linein; 
   }
  return $linein;
}

