#!/usr/bin/env perl
# Old /usr/local/bin/env perl -w

# Author:  Oleg Khovayko
# File Description: eSearch/eFetch calling example

use constant DEBUG => 0;

=item dros ests

~20K per spp

dros.mel : txid 7227

dros.sec : txid 7238  :: no est
dros.sim : txid 7240
dros.yak : txid 7245  ::? testes only, not whole orgn? 11k
dros.ere : txid 7220

dros.ana : txid 7217
dros.per : txid 7234  :: no est
dros.pse : txid 7237

dros.wil : txid 7260
dros.moj : txid 7230
dros.vir : txid 7244
dros.gri : txid 7222

http://www.ncbi.nlm.nih.gov/sites/entrez?\
cmd=Search&db=nucest&term=%28txid7222%5Borgn%5D%20AND%20gbdiv_est%5Bprop%5D%29

http://www.ncbi.nlm.nih.gov/sites/entrez?cmd=Search&db=nucest&term=(txid7260[orgn] AND gbdiv_est[prop])

daphnia magna: txid 35525; 13,257 EST on 081101
daphnia pulex: txid 6669 ; 152,659 EST on 081101  ; 570 nuc seqs, 13 genes, 606 prots ??

=cut

# ---------------------------------------------------------------------------
# Subroutine to prompt user for variables in the next section

sub ask_user {
  print "$_[0] [$_[1]]: ";
  my $rc = <>;
  chomp $rc;
  if($rc eq "") { $rc = $_[1]; }
  return $rc;
}

sub debug { warn @_ if(DEBUG); }

# ---------------------------------------------------------------------------
# Define library for the 'get' function used in the next section.
# $utils contains route for the utilities.
# $db, $query, and $report may be supplied by the user when prompted; 
# if not answered, default values, will be assigned as shown below.

use LWP::Simple;

my $utils = "http://www.ncbi.nlm.nih.gov/entrez/eutils";
my $retpage= (DEBUG) ? 5 : 2000;
my $maxret= 0;

print "NCBI eFetch ESTs per taxon\n";

my $db     = ask_user("Database", "nucest");
my $taxon  = ask_user("Taxonid",    "7260"); 
my $q= "(txid".$taxon."[orgn] AND gbdiv_est[prop])";
my $query  = ask_user("Query",    $q); 
my $report = ask_user("Report",   "fasta");  ## "brief"

$retpage  = ask_user("Result Page",    $retpage); 
$maxret   = ask_user("Result Max",    $maxret); 
my $outfile= ask_user("Output",    "stdout"); 


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
$Count= $maxret if ($maxret and $Count>$maxret);

for($retstart = 0; $retstart < $Count; $retstart += $retpage) {
  my $efetch = "$utils/efetch.fcgi?" .
               "rettype=$report&retmode=text&retstart=$retstart&retmax=$retpage&" .
               "db=$db&query_key=$QueryKey&WebEnv=$WebEnv";
	
  debug "\nEF_QUERY=$efetch\n";     

  my $efetch_result = get($efetch);
  
  debug "---------\nEFETCH RESULT(". ($retstart + 1) . ".." . ($retstart + $retpage) . "): \n";
  
  print $outh $efetch_result,"\n";
  
  if(DEBUG) {
    print "-----PRESS ENTER!!!-------\n";  
    <>;
    }
  else {
    sleep(3);
    }
    
}
