#!/usr/bin/env perl
# evigene/prot/namegenes.pl

use FindBin;
use lib ("$FindBin::Bin","$FindBin::Bin/prot/"); # == evigene/scripts/
use strict;
use Getopt::Long;
use protein_names; # Evigene nameclean() etc

our $DEBUG=0;
our $NAME_KEEPID; # protein_names
our $NAMED_MINALIGN; # protein_names
my ( $blasttab, $refnames, $sortblast, $outnames, $cddnames, $keepid, $blformat, $outformat)= ("") x 9;
my $minalign=0;
$sortblast= undef;
## 2017 outformat = NAMEOUT_TABLE default ?? changed protein_names.pm

my $optok= GetOptions(
 # "config=s", \$config,
 "blasttab|input=s", \$blasttab, # stdin ok; table of  QueryID  RefID  Bitscore Ident Align QueryLen RefLen
 "names|refnames=s", \$refnames, # file required; table of RefID \t name
 "cddnames=s", \$cddnames, # file required; table of RefID \t name
 "output:s", \$outnames, # stdout ok
 "informat=s", \$blformat,  
 "formatout=s", \$outformat,  
 "MINALIGN=s", \$minalign, # NAMED_MINALIGN,  # 0.60 default; allow 60% or 0.60 input
 "sortblast:s", \$sortblast,
 "keepid!", \$keepid, 
 "debug!", \$DEBUG, 
 );

die "usage: $0 -blasttab blast2ref.score.table -names refidname.tab 
  options: -sort Align|Ident|Bitscore -out outname
" unless($optok and $blasttab and $refnames);

if($minalign) {
$minalign= $minalign/100 if($minalign >= 1); $NAMED_MINALIGN=$minalign if($minalign >= 0.01);
}
$NAME_KEEPID= $keepid;
$sortblast= "align" if(not $sortblast and defined($sortblast));

warn "# $0 from $blasttab, $refnames, $sortblast\n" if($DEBUG);

# add call to:   if($outformat =~ /nameclean/) 
#  ($nout, $outname1)= protein_names::cleannamefile($refnames, "$refnames.clean", 0);
#  ($nout, $outname1)= protein_names::cleannamefile($cddnames, "$cddnames.clean", 1);

my ($nout, $outname1)= 
  blast2names( $blasttab, $refnames, $sortblast, $outnames, $cddnames, $blformat, $outformat );
warn "# $0 : nout=$nout to $outname1\n" if($DEBUG or $nout<1);

#-------------------------------
# sub blast2names( $blasttab, $refnames, $sortblast, $outnames, $cddinfo, $blformat, $outformat );  

