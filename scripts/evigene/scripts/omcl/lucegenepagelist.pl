#!/usr/bin/env perl
# lucegenepagelist.pl

=item usage

  perl lucegenepagelist.pl # add options

=item test  
  # query to get all records 'Source:eugenes' is common field:value to all
  curl -s 'http://insects.eugenes.org/aphid/lucegene/search?
  q=Source:eugenes&lib=aphidgenexml
  &start=100&pagesize=50
  &nobatchform=1&norefineform=1&nopagecontrol=1&nolinkstable=1'

  # perl for all reports to set of pages w/ ? pagesize=100; ng=33,000 ; np=3300
  # pagesize=500; np=660

=item notes

 * cut off useless header parts
 
 * parse 1st page for No. matches = ...(\d+)

<td align="right">No. matches = <b>33713</b> 
  of 33716 documents, in 0.0080 sec.</td>
</tr></table>
  
# cut this gbrowse href:
# <a href="http://insects.eugenes.org/cgi-bin/gbrowsenew/gbrowse/aphid/?name=SCAFFOLD511:44016-54290">SCAFFOLD511:44016-54290</a>

=cut

use strict;

my $DEBUG= $ENV{debug}||0;
my $PAGESIZE= 500;
my $startdoc= 0;
#my $startdoc= 5000;

my $uget  = 'curl -s '; # use LWP?

# segregate all the species-specific variables .. add daphnia gene pages
# my $myorg="aphid";
# my $myorg="daphnia";
my $myorg=$ENV{org}||"arthropod";

my $PAGEFILE="genepages_";
my $lurl = 'http://microbe.bio.indiana.edu:7182/lucegene_aphid/search?';
my $query = 'q=Source:eugenes&lib=aphidgenexml'; # all pages here
my $reporturl= 'lookup.lib=aphidgenexml.id=([^">]+)';
my $newrepurl= '/aphid/lucegene/report/'; # rewrite lookup url to this nicer one
#  <a href="lookup?lib=aphidgenexml&id=DGIL_AUG5s102g1t1">

my $format= '&nobatchform=1&norefineform=1&nopagecontrol=1&nolinkstable=1';
my $page  =  '&pagesize='.$PAGESIZE.'&start=';

my $foundpatt= '<a href=.lookup.lib=';
my @cutpatts=(
'\s+(<.DOCTYPE)',
'<script.*</script>(\s</head>)',
'<table class="header".*(\s<h2>Results)',
'<i><a href="index.jsp">Search Page</a></i>',
);

# dang changed tag case: Source > source, ...

if($myorg eq "aphid") {
  $PAGEFILE="aphid_genepages_";
  # $lurl  = 'http://insects.eugenes.org/aphid/lucegene/search?';
  $lurl = 'http://microbe.bio.indiana.edu:7182/lucegene_aphid/search?';
  $query = 'q=Source:eugenes&lib=aphidgenexml'; # all pages here
  $reporturl= 'lookup.lib=aphidgenexml.id=([^">]+)';
  $newrepurl= '/aphid/lucegene/report/'; # rewrite lookup url to this nicer one
  push(@cutpatts, 
  '<a href="http://insects.eugenes.org/cgi-bin/[^>]+>([^<]+)</a>',
  '<a href="/lucegene_aphid/results.jsp.sort=\w+">([^<]+)</a>',
  );
  
} elsif($myorg eq "daphnia") {
  $PAGEFILE="daphnia_genepages_";
  # $lurl  = 'http://wfleabase.org/aphid/lucegene/search?';
  $lurl = 'http://microbe.bio.indiana.edu:7182/lucegene/search?';
  $query = 'q=Source:wfleabase&lib=daphngenexml'; # all pages here
  $reporturl= 'lookup.id=([^">]+)';
  $newrepurl= '/lucegene/report/'; # rewrite lookup url to this nicer one
  push(@cutpatts, 
  '<a href="/cgi-bin/[^>]+>([^<]+)</a>',
  '<a href="/lucegene/results.jsp.sort=\w+">([^<]+)</a>',
  );
  
} elsif($myorg eq "arthropod") {
  $PAGEFILE="arthropod_genepages_";
  #$lurl = 'http://microbe.bio.indiana.edu:7182/lucegene_aphid/search?';
  #x $lurl = 'http://server7.eugenes.org:7151/lucegene_arthropod/search?'; # flamingo/melon2/...
  $lurl= 'http://arthropods.eugenes.org/lucegene_arthropod/search?'; # flaming-inet
  #$query = 'q=source:eugenes&lib=arthropodxml'; # all pages here
  $query = 'q=docid:eugenes&lib=arthropod2xml';
  $reporturl= '/genepage/arthropod/([^">]+)';
  $newrepurl= '/genepage/arthropod/'; # rewrite lookup url to this nicer one
  # this one is right now, but see below need
  # <a href="/genepage/arthropod/ARP1_G7963">
  push(@cutpatts, 
  '<a href="http://server7.eugenes.org:7151/cgi-bin/[^>]+>([^<]+)</a>',
  '<a href="http://insects.eugenes.org/cgi-bin/[^>]+>([^<]+)</a>',
  '<a href="/lucegene_arthropod/resultxsl.jsp.sort=[^"]+">([^<]+)</a>',
  'group-(identity)',
  );

} else {
  die "Dont know about myorg=$myorg";
}


my $maxdoc=0;
for( my $more=1; $more; ) {

  my $url= $lurl . $query . $format . $page . $startdoc; 
  warn "# $uget '$url'\n" if $DEBUG;
  my $out= `$uget '$url'`;
  
  $maxdoc= $1 if ($out =~ m/No. matches = ...(\d+)/s); # No. matches = <b>27493</b> 
  
  # $more= 0 unless($out =~ m/$foundpatt/s); ## need jsp output that marks data.
  $more= 0 unless($out =~ s/$reporturl/$newrepurl$1/sg); ## need jsp output that marks data.
  
  foreach my $cutp (@cutpatts) { $out =~ s/$cutp/$1/sg; }
  
  my $pgfile= $PAGEFILE . sprintf("%05i",$startdoc+1) .".html"; 
  warn "# atdoc=$startdoc; maxdoc=$maxdoc; write $pgfile\n" if $DEBUG;
  open(PG,">$pgfile"); print PG $out; close(PG);
  
  $startdoc += $PAGESIZE;
  $more=0 if($maxdoc > 0 and $startdoc > $maxdoc);
  $more=0 if($DEBUG and $startdoc > 10 * $PAGESIZE);
  sleep(3); #??
}

