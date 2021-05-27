#!/usr/bin/env perl
# uniprotfetch.pl

# ---------------------------------------------------------------------------
my $eutils_root  = "http://www.uniprot.org/";
my $MAPPING='mapping/';
my $BATCH='batch/';
my $UAGENT="eugenes.org";

##my $iproclasstab= "iproclass.tb";
my $iproclasstab= "uniprotgi.tab";
my $GREP="ggrep"; # gnu grep diff names ...

# ---------------------------------------------------------------------------

=item ID Mapping is SLOW; use file

for faster mapping : see sub updateUniprotGiTab() below
  Q. Can I download a complete mapping table between different databases?

  You can download the file iproclass.tb.gz from
  ftp://ftp.pir.georgetown.edu/databases/iproclass/ 
  This file is updated every two weeks.

=cut

use strict;

use LWP::UserAgent;
#use LWP::Simple;
use Getopt::Long;

my $DEBUG=0;
sub debug { warn "# ",@_ if($DEBUG); }

my $DOCTAB="\t";
my $IDTYPE="P_GI";  # == gi number key at uniprot/mapping/ # fixme
my $nhits= 5;
my $BITRANGE = 80; 
my $partsize= 1000; # new 0806
my $saverecords= 1;
my $output= ""; my $outprefix= "sw";
my $UIDkey="ID"; # use uniprot AC field not ID for main doc key ??
   # ^^ PROBLEM: multiple AC values, only one ID value tho
   # can we use ID instead?
   
##my @UNIPROT_tags= qw(ID AC DE GN OS DR); #UNIPROT text keys we want
my @UNIPROT_tags= qw(ID AC DE GN OS DR KW CC); #UNIPROT text keys we want
          # add this? OX   NCBI_TaxID=45351;  GN ORFNames
## option add SQ to keep prot in table
my @DR_DBS = qw(GO InterPro PANTHER Pfam RefSeq); # subkeys from DR == dbxref ; 
  # ^'DR   dbname'
  
my($swparse,$gilist,$blasttable,$saveprot,$updategitab);
my($keeptags, $keepline, %keeptags, @keeptags);

# add output file name/pattern for saverecords and table output?
my $optok= GetOptions(
  "blast=s", \$blasttable, 
  "gilist=s", \$gilist, 
  "keeptags=s", \$keeptags, 
  "keepline=s", \$keepline, # DAMN, this fails too much; faster patt match for full record, e.g. 'OS   Camponotus floridanus'
  "swparse=s", \$swparse, 
  "nhits=i", \$nhits, 
  "BITRANGE=i", \$BITRANGE, 
  "partsize=s", \$partsize, 
  "output=s", \$output, #?
  "tab=s", \$DOCTAB,  
  "iproclasstab=s", \$iproclasstab, #?
  "saverecords!", \$saverecords, 
  "proteinsave!", \$saveprot,
  "updategitab!", \$updategitab, 
  "DEBUG!", \$DEBUG, 
);

die "usage: uniprotfetch.pl -blast blasttable | -gi gilist 
 options: 
 -output outfile or stdout
 -nhits $nhits  : for blast, keep only n best hits; 
 -bitrange $BITRANGE (% of best hit); 
 -partsize $partsize : keep web-server happy for genome sizes; 
 -[no]saverecords : save uniprot full records ($saverecords); 
 -updategitab : update local uniprot-gi-map.table; (web mapping is SLOW)
 -debug
 e.g. 
  cat acyr-allgenes-ncbinr.blastp | grep -v '^#' | sort -k1,1 -k12,12nr \\
  | perl uniprotfetch.pl -blast stdin | sort -k1,1 -k3,3
" 
unless($optok);

# no, some other bug..
# warn "# NOTE: *** -keepline=pattern fails too much ***\n" if($keepline);
# .. this produces same 3,000 subset of 30,000 DAPPU genes
# gzcat uniprot_trembl_invertebrates.dat.gz | perl -ne if(/^ID/) { $p=(/_DAPPU/)?1:0; } print if $p; | /bio/bio-grid/mb/evigene/scripts/uniprotfetchbeta.pl -sw=stdin -out 
# --proteinsave gets all recs, WHY??



updateUniprotGiTab() if($updategitab);

my %blasttab=();
my $gistring="";
my (%doc, %gidocs, @tags); # getGiDocs
my $FAILED=0; # dont die in efetch; save any results 1st
$BITRANGE= $BITRANGE/100; # pct>ratio
my %DR_DBS = map{$_,1} @DR_DBS;

if($keeptags) {
  @keeptags = split(/[,;]+/, $keeptags);
  map{ my($k,$v)=split"=",$_; $keeptags{$k}=$v; $_=$k; } @keeptags;
}

my $outh;
if($output && $output ne "stdout") {
  $outprefix= $output; $outprefix =~ s/\.\w+$//;
  open(OUTPUT,">$output") or die "writing $output";
  $outh= *OUTPUT;
  debug "output to $output\n";
} else {
  $outh= *STDOUT;
}

if($swparse) {
  my( $inh, $stdin, $swdoc)= ("") x 20;

  @tags= @UNIPROT_tags; # DR == dbxref ; UNIPROT text keys we want
  push(@tags, "SQ") if($saveprot);
  push(@tags,@keeptags) if(@keeptags);
  %doc=(); %gidocs=(); 
  @doc{@tags}="";

  debug "swparse $swparse \n";
  local $/ = "//\n";  #read by record
  if($swparse =~ /^stdin|-/) { $inh= *STDIN; $stdin=1; }
  else { open(IF, $swparse) or die "Can't open for read: $!\n"; $inh=*IF; }
  while( $swdoc = <$inh>) {
## add opt here or parseUnip to keep only taxon subset
# OC   Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;
    next if($keepline and not $swdoc =~ /$keepline/);
    parseUnipText( $swdoc, {} ); # this prints/saves to %gidocs
  }
  close $inh unless($stdin);

} elsif($blasttable) {
  my( $inh, $stdin, $lqid,$lbits,$p,$ng)= ("") x 20;
  if($blasttable =~ /^stdin|-/) { $inh= *STDIN; $stdin=1; } 
  else { open(IN,$blasttable) or die "error for blast in: $blasttable"; $inh=*IN; }

  while(<$inh>) {
    next unless(/^\w/);
    my($qid,$sid,@v)=split; 
    my $bits=$v[-1]; 
    if($qid eq $lqid){ $ng++; $p=($ng < $nhits && $bits > $BITRANGE * $lbits)?1:0;} 
    else { $ng=0; $p=1; }
    if($p) {
      ($lbits)=($bits);
      my($gi)= $sid =~ m/gi\|(\d+)/;
      $blasttab{$gi} .= join("\t",$qid,$bits,"hit$ng")."\n" if($gi);#? dont need sid
      }
    ($lqid)=($qid);
  }
  close($inh) unless($stdin); #? err  
  
} else {
  $gilist= shift @ARGV unless($gilist);
  
  undef $/;  #for load whole file
  open IF, $gilist || die "Can't open for read: $!\n";
  $gistring = <IF>;
  close IF;
  debug "Loaded file: [$gistring]\n";
  # Prepare file - substitute all separators to comma
  ## $gistring =~ s/\s+/,/gs;
  # debug "Prepared file: [$gistring]\n";
}

# read in partsize chunks..
my $ok= 1;
if(scalar(%blasttab)) {
  my $ipart=0;
  my @gi= sort keys %blasttab; # may be 100,000+
  $partsize= scalar(@gi) if ($partsize<1);
  while($ok && @gi > 0) {
    $ipart++;
    my $np= min($partsize, scalar(@gi));
    my @giset= splice(@gi,0,$np);  
    $gistring= join " ", @giset;
    $ok= getGiDocs( $gistring, \@giset, $ipart); # partnum?
    # sleep(3) if(@gi > 0); 
  }

} elsif($gistring) {

  $ok= getGiDocs( $gistring);
}

close($outh) if($outh);
exit($FAILED) if($FAILED);

#...........

sub min{ return ($_[1] < $_[0])? $_[1] : $_[0]; }

sub getGiDocs {
  my($gistring, $giset, $ipart)= @_; # or ref?

# outputs: if blasttab, produce table similar to esumprotgi (geneid, bitscore, hit, result:organism, defline, ...)
# also save uniprot-records.txt w/ gi accession keys and from-to gi-unip ID table
# add from uniprot.txt: acc taxid title/DE/organism
  
  debug "getGiDocs part=$ipart\n";

  @tags= @UNIPROT_tags; # DR == dbxref ; UNIPROT text keys we want
  push @tags, "SQ" if($saveprot);
  push(@tags,@keeptags) if(@keeptags);

  %doc=(); %gidocs=(); #@tags=(); # getGiDocs
  @doc{@tags}="";
  
  # $IDTYPE = "P_GI"; # see above
  # FIXME: reuse if exists     my($fn)= $outprefix."_gimap$ipart.tab";
  my $idmap;
  my $savemap= $saverecords;
  my($mapfile)= $outprefix."_gimap$ipart.tab";
  if(-e $mapfile) {
    debug "reading $mapfile\n";
    if(open(F,"$mapfile")) { $idmap= join "", <F>; close(F); $savemap=0; }
  }
  $idmap= mapIds( $gistring, $IDTYPE) unless(defined $idmap); # may be empty
  
  # idmap here is From-To table of IDs
  if($savemap) { # save even if no result?
    ## my($fn)= $outprefix."_gimap$ipart.tab";
    debug "saving $mapfile\n";
    open(F,">$mapfile"); print F $idmap; close(F);
  }
  my(%fmap,%rmap);
  foreach (split "\n", $idmap) {
    next if(/^From/i);
    my($from,$to,$to2)=split; 
    # fixme: use UID first, fall back to UACC if ID missed; but need to know which?
    # getDocs should work on both; 
    $to= $to2 if(!$to && $to2);
    $fmap{$from}= $to; 
    $rmap{$to}= $from if($to and $to !~ m/^UPI0/); # fixme
    debug "$from => $to\n" unless $saverecords; 
    }

  # note, these AC's 'UPI0000D9B387_9544' throw out the batch service
  # drop from rmap
  # FIXME: reuse if exists     my($fn)= $outprefix."_records$ipart.txt";
  my $batchacc= join(" ",sort keys %rmap);
  
  my $docdata;
  my $savedoc= $saverecords;
  my($docfile)= $outprefix."_records$ipart.txt";
  if(-e $docfile) {
    debug "reading $docfile\n";
    if(open(F,"$docfile")) { $docdata= join "", <F>; close(F); $savedoc=0; }
  }
  $docdata= getDocs( $batchacc) unless(defined $docdata); 
  
  if($savedoc) { # save even if no result?
    debug "saving $docfile\n";
    open(F,">$docfile"); print F $docdata; close(F);
    }
  parseUnipText($docdata, \%rmap); # this prints/saves to %gidocs
  
  if(%blasttab) {
    unless(ref $giset) { $giset= [ sort keys %blasttab ]; }
    print $outh "#".join("\t",qw(query bits hit ncbigi), @tags),"\n";
    foreach my $gi (@$giset) {
    
      my $dbid = $fmap{$gi}; # can this be either ID or ACC? gidocs{} needs right key
      my $doc  = $gidocs{$gi} || $gidocs{$dbid} || "nodoc gi:$gi";

      my $blast= $blasttab{$gi} || "noblast";
      my @blast= split "\n", $blast;
      map{ print $outh "$_\t$gi\t$doc\n"; } @blast;
      ## qid, bits, hitnum, Gi, Caption, TaxId, Title
    }
  }

  return ($FAILED) ? 0 : 1;
}


=item UniProt fields

Now keeping: ID, AC, GN, DE, OS, DR (selected)

# DR   Pfam; PF00368; HMG-CoA_red; 1.
# DR   PANTHER; PTHR11621; UBQ-conjugat_E2; 1.
# DR   InterPro; IPR015943; WD40/YVTN_repeat-like.
# DR   GO; GO:0000159; C:protein phosphatase type 2A complex; IDA:UniProtKB.

KW : add keyword field but for trivial

What are these special comments? these are perhaps most useful
function info... capture some/all. How to summarize?

microbe% cat acyr-augmap5a-pt1-nra_record*.txt | grep 'CC   -!' | perl -ne's/:.*$//; print;' | sort | uniq -c | sort -k1,1nr | more
1028 CC   -!- SIMILARITY    **
 284 CC   -!- SUBCELLULAR LOCATION  ** 
 282 CC   -!- FUNCTION       **
 175 CC   -!- CAUTION
 170 CC   -!- SUBUNIT
 121 CC   -!- CATALYTIC ACTIVITY
  94 CC   -!- ALTERNATIVE PRODUCTS
  81 CC   -!- TISSUE SPECIFICITY    *?
  50 CC   -!- COFACTOR
  45 CC   -!- SEQUENCE CAUTION
  36 CC   -!- DOMAIN
  36 CC   -!- PATHWAY               ??**
  34 CC   -!- DEVELOPMENTAL STAGE    **
  34 CC   -!- PTM
  28 CC   -!- INTERACTION
  25 CC   -!- MISCELLANEOUS
  17 CC   -!- INDUCTION
  16 CC   -!- DISEASE            ?**
  14 CC   -!- WEB RESOURCE
  12 CC   -!- ENZYME REGULATION
   3 CC   -!- BIOPHYSICOCHEMICAL PROPERTIES
   2 CC   -!- MASS SPECTROMETRY
   2 CC   -!- RNA EDITING

Eg,
  46 CC   -!- SIMILARITY: Contains 1 RING-type zinc finger.
  11 CC   -!- SIMILARITY: Contains 1 helicase C-terminal domain.
  10 CC   -!- SIMILARITY: Contains 1 homeobox DNA-binding domain.
  22 CC   -!- SUBCELLULAR LOCATION: Nucleus.
  10 CC   -!- SUBCELLULAR LOCATION: Cytoplasm.
   7 CC   -!- PATHWAY: Protein modification; protein glycosylation.
   6 CC   -!- PATHWAY: Protein modification; protein ubiquitination.
   6 CC   -!- FUNCTION: May be involved in transcriptional regulation.
   1 CC   -!- FUNCTION: GTP-binding protein that functions as an allosteric
   2 CC   -!- DEVELOPMENTAL STAGE: At E15-E17, mainly in the developing retina,
   1 CC   -!- DEVELOPMENTAL STAGE: Expressed throughout embryonic development.
   2 CC   -!- TISSUE SPECIFICITY: In midgut and Malpighian tubules.
   2 CC   -!- TISSUE SPECIFICITY: Ubiquitous.
   1 CC   -!- DISEASE: Antigen in chronic rheumatoid arthritis and in the

DE: tag, drop Contains, promote RecName: Full= and Short=

DE   RecName: Full=14-3-3 protein beta/alpha;
DE   AltName: Full=Protein kinase C inhibitor protein 1;
DE            Short=KCIP-1;
DE   AltName: Full=Protein 1054;
DE   Contains:
DE     RecName: Full=14-3-3 protein beta/alpha, N-terminally processed;

=cut

my $parseHist="";

sub parseValue {
  my ($k,$v)= @_;
  $v =~ s/\.$//; $v =~ s/;\s*$//; 
  if($k eq "ID") { $v =~ s/^\s+//; $v =~ s/\s.*$//; }
  ## add SQ for protein seq only, not header.
# SQ   SEQUENCE   246 AA;  28082 MW;  6BE1A9BF97468017 CRC64;
#     MTMDKSELVQ KAKLAEQAER YDDMAAAMKA VTEQGHELSN EERNLLSVAY KNVVGARRSS
  elsif($k eq "SQ") { $v =~ s/\s*SEQUENCE.*//; $v =~ s/\s+//g; }
  elsif($k eq "AC") { $v =~ s/;\s*/,/g; }
  # ^^ FIXME for AC at least handle multiple values: AAA; BBB; CCC;
  elsif($k eq "GN") { $v =~ s/(Name|ORFNames|OrderedLocusNames)=//g; $v =~ s/;\s*/,/g; }
  elsif($k eq "OX") { $v =~ s/NCBI_TaxID=//; }
  elsif($k eq "DE") { 
    $parseHist="$k\t$v" if($v =~ /Contains:|AltName:/);
    return "" if($parseHist =~ /DE\tContains:/);
    return "" if($parseHist =~ /AltName:/); # unless saveAltName ?   $parseHist =~ /DE\tRecName:/);
    #NO?# $v =~ s/RecName://;
    }
  elsif($k eq "CC") { 
    return "" unless ($v =~ s/\-\!\-\s*//); # keep only 1st line
    my($sk,$sv)= split ": ",$v,2;
    if($sv && $sk =~ /(FUNCTION|SIMILARITY|TISSUE|LOCATION|PATHWAY|DEVELOPMENTAL|DISEASE)/) {
      $sk= substr($1,0,3); $sv =~ s/;\s*/,/g; 
      $v= "$sk:$sv";
    } else { return ""; }
    }
  elsif($k eq "KW") { 
    $v =~ s/\s*(Complete proteome|Alternative splicing|Polymorphism);?//g;
    $v =~ s/;\s*/,/g;
    }
  elsif($k eq "DR") { 
    my($db,$acc,$na,@dx)= split( /\;\s*/, $v );
    return "" unless($DR_DBS{$db});
    $v= ($na and $na ne "-") ? "$db:$acc/$na" : "$db:$acc";
    }
  $parseHist="$k\t$v";
  return $v;
}


sub parseUnipText {
  my($result, $idmap)= @_;
  my %tags= map{$_,1} @tags;
  my $lastk;
  ## my $saveresult= ($saverecords) ? 1 : 0;

  foreach (split "\n", $result) {    
    if(m,^(\w+)\s+(\S.+),) { 
      my ($k,$v)= ($1,$2);  
      if($tags{$k}) {
      $v= parseValue($k,$v);
      if($v) {
      $doc{$k} .= "," if $doc{$k};
      $doc{$k} .= $v; 
      if($k eq $UIDkey && $idmap) { 
        foreach my $v1 (split(",",$v)) { $doc{'ncbigi'}= $idmap->{$v1} if $idmap->{$v1}; } 
        } 
      }
      }
      $lastk= $k;
    } elsif(m,^//,) { 
      # saverecord( $doc{$UIDkey}, $result) if (0 < $saveresult--);
      printdoc();  
    } elsif(m,^\s+(\S.+),) {
      my ($k,$v)= ($lastk,$1); 
      if($tags{$k}) {
      $v= parseValue($k,$v);
      $doc{$k} .= " ".$v if($v); 
      }
    }
  }
  printdoc();  
}

# sub saverecord {
#   my($id,$result)= @_;
#   return unless $result;
#   $id= time() unless($id);
#   debug "saving sw$id.txt\n";
#   open(TXT,">>sw$id.txt"); print TXT $result; close(TXT);
# }


sub savedoc {
  my $val= join("\t", @doc{@tags});
  if($doc{$UIDkey}) { $gidocs{ $doc{$UIDkey} }= $val; }
  if($doc{'ncbigi'}) { $gidocs{ $doc{'ncbigi'} }= $val; }   # also use ncbigi as key?  
  @doc{@tags}="";
}

my $didhead=0;
sub printdoc {
  my $saved= ($swparse) ? 0 : scalar(%blasttab) ? 1 : 0;
  if($saved) { savedoc(); }
  elsif(scalar(%doc)) {  #? bad: @doc{@tags}
    my @showtags= @tags;

    # filtertag : regex per tag like OC to keep/remove entries
    # are keeptags part of showtags or not? 
    if(@keeptags) {
    @showtags = grep { !$keeptags{$_} } @showtags;
    my $skip=0;
    foreach my $kt (@keeptags) {
	      my $v= $doc{$kt}; 
        my $keep= $keeptags{$kt} or next;
        $skip=1 unless($v =~ m/$keep/);
    }
    next if($skip);
    }
 
    unshift(@showtags, qw(ncbigi))  unless($swparse);
    print $outh "#".join($DOCTAB,@showtags),"\n" unless ($didhead++);
    print $outh join($DOCTAB, @doc{@showtags}),"\n";
  }
  @doc{@tags}="";
}


=item iproclass.tab

 mapIds is SLOW ... use instead  iproclass.tb.gz from
 ftp://ftp.pir.georgetown.edu/databases/iproclass/ 
  7372343 records (jul08)

This table includes the following IDs (or ACs) delimited by tab:
1. UniProtKB accession (UniProtKB accession with taxon_id)
2. UniProtKB ID
3. EntrezGene
4. RefSeq
5. NCBI GI number     **
6. PDB
7. Pfam
8. GO
9. PIRSF
10. IPI
11. UniRef100
12. UniRef90
13. UniRef50
14. UniParc
15. PIR-PSD accession
16. NCBI taxonomy
17. MIM
18. UniGene
19. Ensembl
20. PubMed ID
21. EMBL/GenBank/DDBJ
22. EMBL protein_id


=cut

sub updateUniprotGiTab {
  # fixme; keep empty column for 2,3 (EG,RS)
  
  # check date on iproclasstab; update if > 2wk old or on request ...
  my $cmd= <<'EOC';
curl -s ftp://ftp.pir.georgetown.edu/databases/iproclass/iproclass.tb.gz |\
perl -ne'@v=split"\t"; print join("\t",$v[0],$v[1],".",".",$v[4]),"\n";' > 
EOC
  $cmd.= $iproclasstab."\n";
  print STDERR "# Update $iproclasstab :\n", $cmd;
  print STDERR "# You can use the full iproclass.tb also, named $iproclasstab\n";
  die "\nPlease run the above command to update local uniprot-gi-map table.\n";
  # system($cmd); # or debug print this
}

sub mapIds {
  my($gistring,$IDTYPE)= @_;  
  return unless($gistring =~ /\d/);

  my $result= mapIdsFile($gistring,$IDTYPE);
  $result= mapIdsEutil($gistring,$IDTYPE) unless($result);
  $result ||= "";
  return $result;
}

sub mapIdsFile {
  my($gistring,$IDTYPE)= @_;  
  my $result="";

  if($GREP && -f $iproclasstab) { # use indexed by something? ggrep -F -f gilist iproclass.tb
     $gistring =~ s/[,\s\+]+/\n/g; # joiner should be space now ...
     my %gihash= map{ $_,1 } split "\n", $gistring;
     my ($nfound,$ngot)= (0, scalar(keys %gihash));
     my $gifile= $outprefix."tmp.gilist";
     open(T,">$gifile"); print T $gistring,"\n"; close(T);
     my $cmd= "$GREP -F -f $gifile $iproclasstab";
     debug "mapIdsFile: $cmd\n";
     open(T,"$cmd |");
     while(<T>){
      my($uac,$uid,$eg,$rf,$gi,@more)=split"\t";
      my @gi=  split(/;\s*/,$gi);
      foreach my $gi (@gi) {
        if($gihash{$gi}) { 
          # $result .= join("\t",$gi,$uid,$uac)."\n"; 
          # $gihash{$gi}++; 
          # $gihash{$gi}= join("\t",$uac,$uid); # uac always there, uid not
          $gihash{$gi}= join("\t",$uid,$uac); # uac always there, uid not
          ## FIXME: may be many uac / uniprot record; only one uid ??
          $nfound++; last; 
          }
        }
      }
      close(T); unlink($gifile);
      debug "mapIdsFile: want: $ngot; found: $nfound\n";
      #?? warn if $nfound > 0 but $nfound/$ngot < 0.50 ??
      foreach my $gi (sort keys %gihash) {
        my $v= $gihash{$gi}; $v="" if ($v==1);
        $result.= join("\t",$gi,$v)."\n";
        }
  }
  return $result;
}

sub mapIdsEutil {
  my($gistring,$IDTYPE)= @_;  
  return unless($gistring =~ /\d/);
    
  # my $ngi= $gistring =~ tr/,/,/; $ngi++ if($gistring);
  # print "# ncbi esummary db=$db_name, ngi=$ngi, part=$ipart\n";
  
  my $form_data = {
    from => $IDTYPE, to => $UIDkey, ##? 'ACC' or 'ID ?
    query => $gistring, format => 'tab',
    };  

  my $result= callEutil( $eutils_root . $MAPPING, $form_data);
  return $result;
 
}

sub getDocs {
  my($idstring)= @_; 
  return "" unless($idstring =~ /\w/);

  my $form_data = {
    query => $idstring,
    format => 'txt',
    };  

  my $result= callEutil( $eutils_root . $BATCH, $form_data);
  return $result;
  
}

sub callEutil {
  my($url,$params)= @_;
  
  my $agent = LWP::UserAgent->new;
  $agent->agent($UAGENT);

  push @{$agent->requests_redirectable}, 'POST';
  debug "Submitting to $url ...\n";
  
  my $response = $agent->post("$url", $params);
  
  ## this delay is REALY BAD for mapping; almost nothing for batch fetch
  my $step=0; 
  while (my $wait = $response->header('Retry-After')) {
    $step++; debug "Waiting ($step; $wait s)...\n";
    sleep $wait;
    $response = $agent->get($response->base);
  }
  
  unless($response->is_success) { 
    $FAILED++; 
    warn 'Failed, got ' . $response->status_line . 
        ' for ' . $response->request->uri . "\n"
  }
  
  return  $response->content;
}

#--------------------------

__END__


#... add select DR records w/ info

my @DR_DBS = qw(GO InterPro PANTHER Pfam RefSeq);

melon.% cat *records*.txt | grep '^DR' | sort | uniq | perl -ne's/;.*$//; print; ' | sort | uniq -c | more
 192 DR   ArrayExpress
  18 DR   BioCyc
  36 DR   CleanEx
   9 DR   DIP
2006 DR   EMBL
 220 DR   Ensembl
 191 DR   FlyBase
  14 DR   GO         <<< *
 119 DR   Gene3D
   2 DR   GeneDB_Spombe
1081 DR   GeneID
  12 DR   GenomeReviews
  50 DR   GermOnline
   8 DR   Gramene
  12 DR   H-InvDB
  25 DR   HGNC
 159 DR   HOGENOM
 180 DR   HOVERGEN
   1 DR   HPA
  86 DR   HSSP
  40 DR   IntAct
 738 DR   InterPro    <<< *
 548 DR   KEGG        <<< no
  18 DR   LinkHub
  11 DR   MEROPS
  22 DR   MGI
  14 DR   MIM
  87 DR   NMPDR
 113 DR   PANTHER   <<< *
  14 DR   PDB
  14 DR   PDBsum
  64 DR   PIR
  40 DR   PIRSF
 109 DR   PRINTS
 426 DR   PROSITE
   2 DR   PeptideAtlas
   3 DR   PeroxiBase
 494 DR   Pfam          <<< *
  19 DR   PharmGKB
  23 DR   PhosphoSite
  61 DR   ProDom
   4 DR   RGD
  16 DR   Reactome
1183 DR   RefSeq      <<< *
 215 DR   SMART
  92 DR   SMR
   1 DR   SWISS-2DPAGE
   1 DR   SagaList
   4 DR   TAIR
   1 DR   TIGR
  49 DR   TIGRFAMs
 482 DR   UniGene
 388 DR   VectorBase
   3 DR   WormBase
   3 DR   WormPep
  36 DR   ZFIN
   2 DR   dictyBase


melon.% 

#---------------

=item beta.uniprot.org perl eg


"curl -s -d 'query=Q0III0+Q966Z9+Q967E3+Q967E2+Q66IU1+A7SWX2&format=text'   http://www.uniprot.org/batch/

curl -i 'http://www.uniprot.org/batch/?query=Q0III0+Q966Z9+Q967E3+Q967E2+Q66IU1+A7SWX2'

&format=text'
format=text&

# cant batch fetch in one step ; need to hassle w/ second url from unibeta
&format=text

http://www.uniprot.org/batch/?query=Q0III0+Q966Z9+Q967E3+Q967E2+Q66IU1+A7SWX2
>>>
http://www.uniprot.org/jobs/Z4XI.txt




#   use http://beta.uniprot.org/uniprot/ ; http://beta.uniprot.org/faq/28
# http://beta.uniprot.org/uniprot/?query=organism:9606+AND+antigen&format=tab&compress=yes&columns=id,reviewed,protein names
# formats:  text = full entry; tab = table of default/spec columns

columns  	comma-separated list of values, e.g. for UniProtKB: 
citation | clusters | comments | database | domains | domain | ec
| id | entry name | existence | families | features | genes | go
| go-id | interpro | interactor | keywords | keyword-id |
last-modified | length | organism | organism-id | pathway |
protein names | reviewed | score | sequence | 3d | subcellular
locations | taxon | tools | version | virus hosts

## tool_example.pl ##
use strict;
use warnings;
use LWP::UserAgent;

my $base = 'http://beta.uniprot.org';
my $tool = 'mapping';
my $params = {
  from => 'ACC', to => 'P_REFSEQ_AC', format => 'tab',
  query => 'P13368 P20806 Q9UM73 P97793 Q17192'
};

my $agent = LWP::UserAgent->new;
push @{$agent->requests_redirectable}, 'POST';
print STDERR "Submitting...\n";
my $response = $agent->post("$base/$tool/", $params);

while (my $wait = $response->header('Retry-After')) {
  print STDERR "Waiting ($wait)...\n";
  sleep $wait;
  print STDERR "Checking...\n";
  $response = $agent->get($response->base);
}

$response->is_success ?
  print $response->content :
  die 'Failed, got ' . $response->status_line . 
    ' for ' . $response->request->uri . "\n";
    
=cut

=item mapped GI ids

From      To
115496628 Q0III0
115620465	
118086617	
118094585	
125819864	
13928075	Q966Z9
13928077	Q967E3
13928079	Q967E2
148229713	Q66IU1
156356358	A7SWX2


Q0III0 Q966Z9 Q967E3 Q967E2 Q66IU1 A7SWX2

=cut
