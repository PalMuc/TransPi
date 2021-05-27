#!/usr/bin/env perl
# protein2cds.pl

=item about
  
  backtranslate protein to coding sequence dna

=cut

use FindBin qw($Bin);
use lib ("$Bin","$Bin/../lib","$Bin/../");

use strict;
use warnings;
use Getopt::Long;
use cdna_proteins;

my ($debugin, $cdnaseq,$output,%didid,$useAmbiguousNucs,$useFixedNucs);

my $optok= GetOptions(
  "aa|protein|input=s", \$cdnaseq,  # use instead @input/@ARGV
  "output:s", \$output,
#   "minaa=i", \$MINAA,  # not used ?
#   "fullorf:s", \$ORF_FULLvPART,  
#   "nostopcodon", \$NoStopCodon, 
#   "Selenocysteine|selc!", \$USESelenocysteine,  
  "ambiguous!", \$useAmbiguousNucs, # -ambig | -noambig
  "fixed", \$useFixedNucs, # default?
  "debug:i", \$debugin, 
  );

$cdnaseq= shift @ARGV unless($cdnaseq);  

die "usage: protein2cds.pl [ -ambiguous ] myproteins.aa -out myproteins.cds or >myproteins.cds
" unless($optok and  $cdnaseq);

$DEBUG=$debugin if(defined $debugin);
if(defined $useFixedNucs and not (defined $useAmbiguousNucs)) { $useAmbiguousNucs= not $useFixedNucs; }

my $outh= *STDOUT;
if(defined $output and not $output) { # use cdnaseq name
  ( $output= $cdnaseq ) =~ s/\.(aa|faa|pep|fasta|fa).*//; $output.=".cds";
}
if($output and $output!~/stdout|^-/) { open(OUT, ">$output") or die $output; $outh= *OUT; }

# filter out dup ids; buggy input.aa ..

sub putcds {
  my($id,$tag,$hd,$aa)= @_;
  if($aa and not $didid{$id}) {
    my $cds= backtranslate_protein($aa,$useAmbiguousNucs); 
    my $cdslen=length($cds);
    if($cdslen>0) {
      $cds  =~ s/(.{60})/$1\n/g;
      print $outh ">",$id," type=$tag; len=$cdslen; $hd\n";
      print $outh $cds;
      print $outh "\n" unless($cds=~/\n$/);
      $didid{$id}++;
      }
    return $cdslen;
  }
  return 0;
}

MAIN: { 
  my $ovh; my $ok= 0; 
  my($nin,$tnin,$tngood,$tnsplit,$tnskip)=(0) x 9;
  if($cdnaseq =~ /.gz$/) { $ok= open(OVR,"gunzip -c $cdnaseq |");  $ovh= *OVR; }
  elsif($cdnaseq =~ /stdin|^\-$/) { $ok=1; $ovh= *STDIN; }
  elsif($cdnaseq) { $ok= open(OVR,$cdnaseq); $ovh= *OVR; }
  die "bad -aa=$cdnaseq" unless($ok);

  my($id,$tag,$hd,$aa);
  $tag="cds"; # dont add to id, add attr type=tag;
  while(<$ovh>) {
    if(/^>(\S+)\s*(.*)$/) { 
      my($inid,$inhd)=($1,$2);
      $tngood++ if(putcds($id,$tag,$hd,$aa));
      ($id,$hd,$aa)=($inid,$inhd,"");
      $nin++;  $tnin++;
    } elsif(/\w/) { chomp; $aa.= uc($_); } 
  } close($ovh);
  
  $tngood++ if(putcds($id,$tag,$hd,$aa)); 

  warn "#protein2cds from aaseq:$cdnaseq nin=$tnin nout=$tngood\n" if $DEBUG;
}  
