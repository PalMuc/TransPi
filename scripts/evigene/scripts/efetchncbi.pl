#!/usr/bin/env perl -w

# Author:  Oleg Khovayko
# File Description: eSearch/eFetch calling example


=item dros ests

dros.mel : txid 7227

http://www.ncbi.nlm.nih.gov/sites/entrez?\
cmd=Search&db=nucest&term=%28txid7222%5Borgn%5D%20AND%20gbdiv_est%5Bprop%5D%29

http://www.ncbi.nlm.nih.gov/sites/entrez?cmd=Search&db=nucest&term=(txid7260[orgn] AND gbdiv_est[prop])

=cut

use LWP::Simple;
use Getopt::Long;

my $ask=1;
my $debug = 0;
sub debug { warn @_ if($debug); }

# ---------------------------------------------------------------------------
# Subroutine to prompt user for variables in the next section

sub ask_user {
  print "$_[0] [$_[1]]: ";
  my $rc= $_[1];
  if($ask or not $rc) { $rc = <>; chomp $rc; }
  if($rc eq "") { $rc = $_[1]; }
  return $rc;
}


# ---------------------------------------------------------------------------
# Define library for the 'get' function used in the next section.
# $utils contains route for the utilities.
# $db, $query, and $report may be supplied by the user when prompted; 
# if not answered, default values, will be assigned as shown below.

my $utils = "http://www.ncbi.nlm.nih.gov/entrez/eutils";
my($taxon,$db,$query,$report,$retpage,$maxret,$outfile);
  $taxon= 0;
  $db="nucest";
  $maxret= -1;
  $retpage= ($debug) ? 5 : 2000;
  $outfile="stdout";
  $report="fasta";
  
my $optok= GetOptions(
  "db=s", \$db, 
  "query=s", \$query, 
  "report=s", \$report, 
  "output=s", \$outfile, 
  "taxon|taxid=i", \$taxon, 
  "retpage=i", \$retpage, 
  "maxret=i", \$maxret, 
  "debug!", \$debug, 
  "ask!", \$ask, 
  );

print "NCBI eFetch nucest/nucleotide/protein per taxon\n";

$taxon  = ask_user("Taxonid", $taxon); 
die "need Taxonid\n" unless ($taxon);

# add more db div opts: nucest protein nucleotide  
$db     = ask_user("Database", $db);

# biomol_mRNA[prop] NOT gbdiv_est[prop]
# (txid7425[orgn] AND biomol_mRNA[prop] NOT gbdiv_est[prop])
# txid7425 AND srcdb_refseq_known[prop] AND mrna[filter]
# txid7425 AND srcdb_refseq_inferred[prop] AND mRNA[filter]

my $q= "txid".$taxon."[orgn]";
if($db =~ /est/) { $q.=" AND gbdiv_est[prop]"; }
elsif($db =~ /nucleotide/) { $q.=" AND srcdb_refseq_known[prop] AND mrna[filter]"; }
elsif($db =~ /protein/) { $q.=" AND srcdb_refseq_known[prop] AND mrna[filter]"; }
$q = "($q)";

$query  = ask_user("Query",    $q); 
$report = ask_user("Report",   $report);  ## "brief"

$maxret   = ask_user("Result Max",  $maxret); 
$retpage= $maxret if ($maxret>0 and $retpage>$maxret);
$retpage  = ask_user("Result Page", $retpage); 
$outfile  = ask_user("Output",    $outfile); 


# ---------------------------------------------------------------------------
# $esearch cont¡ins the PATH & parameters for the ESearch call
# $esearch_result containts the result of the ESearch call
# the results are displayed ¡nd parsed into variables 
# $Count, $QueryKey, and $WebEnv for later use and then displayed.

my $esearch = "$utils/esearch.fcgi?" .
              "db=$db&retmax=1&usehistory=y&tool=eugenes.org&term=";

my $esearch_result = get($esearch . $query);

debug "\nESEARCH RESULT: $esearch_result\n";

$esearch_result =~ 
  m|<Count>(\d+)</Count>.*<QueryKey>(\d+)</QueryKey>.*<WebEnv>(\S+)</WebEnv>|s;

my $Count    = $1;
my $QueryKey = $2;
my $WebEnv   = $3;

debug "Count = $Count; QueryKey = $QueryKey; WebEnv = $WebEnv\n";

# ---------------------------------------------------------------------------
# this area defines a loop which will display $retmax citation results from 
# Efetch each time the the Enter Key is pressed, after a prompt.

if($outfile =~ /stdout/i) { $outh= *STDOUT; }
else { open(OUT,">$outfile") or die $outfile; $outh= *OUT; }

my $retstart;
$retpage= $maxret if ($maxret>0 and $retpage>$maxret);
$Count= $maxret if ($maxret>0 and $Count>$maxret);

for($retstart = 0; $retstart < $Count; $retstart += $retpage) {
  my $efetch = "$utils/efetch.fcgi?" .
               "rettype=$report&retmode=text&retstart=$retstart&retmax=$retpage&" .
               "db=$db&query_key=$QueryKey&WebEnv=$WebEnv";
	
  debug "\nEF_QUERY=$efetch\n";     

  my $efetch_result = get($efetch);
  
  debug "---------\nEFETCH RESULT(". ($retstart + 1) . ".." . ($retstart + $retpage) . "): \n";
  
  print $outh $efetch_result,"\n";
  
  if($debug) {
    print "-----PRESS ENTER!!!-------\n";  
    <>;
    }
  else {
    sleep(3);
    }
    
}
