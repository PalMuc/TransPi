#!/usr/bin/env perl
# evigene2genbanktbl.pl

=item about

    evigene2genbanktbl.pl  : convert genes.gff per Evigene structure/annotation (GFF.v3)
        to GenBank genome submit table.
  
=item usage

  $evigene/scripts/evigene2genbanktbl.pl -debug \
  -change genes/gbsubmit.changelist \
  -conf genes/evigene_gbsubmit.conf \
  -in genes/pub3h/genes_pub3h.gff.gz  -out submit/pub3h.tbl > & log.sub22 &
	
  input is genes.gff in ordered 3-level structure records (gene/mRNA+/exon,CDS,other)
  it can contain non-gene annotations as well (e.g. transposons)
  
  tbl2asn can be run from this script or separately on its results
    tbl2asn -p submit/ -t ./cacao3g_gbsubmit.sbt -a s -V vb -X E -Z submit/pub3h.discrep
  
=item configuration

  evigene_gbsubmit.conf contains key=value configurations for
    inputs (genome, qual, proteins, ...)
    feature annotation configurations (mapping from evigene to ncbi)
    evidence value cut-offs, etc.
  
  gbsubmit.changelist contains changes to genes.gff per mRNA ID,
  drops, rename (ID), strand +/-, spanfix, defer, rename and others.
  
=item Outputs as per GMODTools/Bulkfiles/ GenbankSubmitWriter.pm
  
  outputs  
   mygenome.tbl    : feature table  
   mygenome.fsa    : genome dna 
   mygenome.qvl    : assembly quality values
   mygenome.pep    : protein aa  
   mygenome.sbt              : stub template 
  
  outputs of NCBI tbl2asn:
  mygenome.sqn : ASN.1 record to submit to NCBI
  mygenome.val : errors & warnings
  mygenome.gbf : Genbank format for review

  see also
  gmod/schema/GMODTools/data/genomes/
    Anopheles_gambiae_str._PEST/anogam_20080511/genbanksubmit/
    Drosophila_melanogaster/current/genbanksubmit/

=item changelist actions

    drop
    notpartof
    MisMatchAA
    defer
    gene    geneid
    rename  newid
    strand  +/-
    spanfix location
    # retype  
    # relocate

=item evigene2genotbl_kfish2.pl temp

  -- test runs for kfish2 geno.gff to genbank submit tbl
  
  -- no need for nameclean(), done for TSA submit, use same gene names
  
  #Ko# turn off for kfish2 tests
  #Kf# addition
  
=cut

use FindBin;
#Ko# use lib ("$FindBin::Bin","$FindBin::Bin/prot/"); # == evigene/scripts/

use strict;
use Getopt::Long;
use File::Basename;
#Ko# use protein_names; # Evigene nameclean() etc

use constant VERSION  => '20130326';  

use constant CHR_RENAME   => 1; # on, works right
use constant TBL2ASNready => 0; # DONT run tbl2asn yet, leave to operator

#KF------------------------------------------
#KF: copy vars from protein_names.pm
use constant { NAMEDIFF_MINOR => 1, NAMEDIFF_MAJOR => 2, };
use constant MAXNAMELEN => 70;
use constant CDD_MAXNAMELEN =>  MAXNAMELEN - 10; # was 39;
our $NAME_NONE = "Unknown|Uncharacterized conserved|Uncharacterized|Hypothetical|noname";  
our $NAME_UNKUNIP= 'Uncharacterized protein'; # Uniprot std
our $NAME_UNKNCBI= 'hypothetical protein';  # NCBI prefers this to Uncharacterized .. special handling
our $NAME_UNK  = $NAME_UNKUNIP;
our $NAME_UNKCDD = 'Domain of unknown function'; # DUF nnn stuf
our $NAME_UNKADDLOCUS= 1; # policy, add/not the gene id to end of UNK names
our $MIN_NAMEIDENT = 35;  # min similar% for protein evid align/equivalence, for naming; JCVI uses 35%  
our $MIN_IDLIKE   = 15;  # low for now; 15..20 seems right;add Note with loqualname
our $MIN_CERTAIN  = 66;  # putative naming, what?
our $NAMED_MINALIGN = 0.60; # or 60% or $MIN_CERTAIN ??
our $NAME_IDPATT= 
      '(?:At|AT|Os)\d{1,2}[Gg]\d{3,}|DUF\d\d+|OSJNB\w\d+[\w.]+|' # plants 
      .'(:?clone|[cC]DNA) \w+\d\d+[,]?|C\d+orf\d+\w?|DKFZp\d+\w+|(?:ENS\w+\d|Zgc:|DDB_G|LOC|FLJ|KIAA|LINC)\d\d+|' # human/vert, may have several..
      .'(?:AGAP|AAEL|ACYPI|BcDNA.[CG][A-Z]|CG|G[A-Z]|GLEAN_|NM_)\d\d+'; # bugs;
#KF------------------------------------------

my $DEBUG=1;

my %config = (); # replace other config hashes w/ this 2-level: config{part}{key} = val
my %evidence= ();
my %public_options = ();
my %geneset=(); # from config
my %programs= (); # add to other evigene set

my $mrnatypes='mRNA';
my $exontypes='exon|CDS';
my $passtypes="";

my $MIN_IDENTITY = 35;  # min similar% for protein evid align/equivalence, for naming; JCVI uses 35%  
my $MIN_ESTIDENT = 10; # min align% to keep  EST/Rna evidence ; not used for naming

## /// move to protein_names  ///
# my $NAME_NONE = "Unknown|Uncharacterized|Hypothetical";  
# my $NAME_UNK  = "Uncharacterized protein"; # uniprot  
# my $NAME_UNKADDLOCUS= 1; # policy, add/not the gene id to end of UNK names
# my $MIN_NAMEIDENT = 35; # == $MIN_IDENTITY;  # min similar% for protein evid align/equivalence, for naming; JCVI uses 35%  
#   # -- note MIN_ID is used for both naming and keeping prot evidence, split this? use lower IDLIKE for evid?
#   # MIN_IDENTITY == MIN_NAMEIDENT now
# my $MIN_IDLIKE   = 15;  # low for now; 15..20 seems right;add Note with loqualname
# my $MIN_CERTAIN  = 66;  # putative naming, what?
## /////////

my $TINYAA = 30;  # minsize to turn into misc feature/fragment
my $ADDGENE = 1; # config
my $PEPfromGFF=0;# config
my $PutPEPfromFile=0;# config
my $FAWIDTH = 60;
my $CDS_OFFBY_FIX=1; # config 
my $SKIPNOANN= 0;

my $config;
my @configadd;
my $IgnoreOutOfOrder=0;
my $notbl2asn=0;

my $pubidnum_start=0;
my $ALTKEY="t"; # but see config altid_format
my $MID=""; # make-id = temp/prelim bestgenes id prefix
my $MySRC="Evigene"; ## DGILmix8d == anv ? .. replace with my DBID ?
my $MyIDPrefix='\w';
my ( $ingff, $changelist, $outgff, $pepout, $logout, $inscore, $genome, $proteins,$thischr )= ("") x 20;
$logout= undef;

use constant SPLITGENE_FIX => 1; # 20120829
my %splitgene=(); # hash mainid>partid ; partid == Thecc1EG015312t1a,b,c.. for transcid,protid

# changelist == input list of drop/change ids ; changes? retype (mrna > misc), ??

my $optok= GetOptions(
  "config=s", \$config,
  "input=s", \$ingff,
  "output=s", \$outgff,
  "changelist=s", \$changelist, # or from config
  "logout:s", \$logout,
  #... easy replace of config vals
  ## input for added annots: genescore tables
  "annot|score=s", \$inscore,
  "genome=s", \$genome,
  "proteins=s", \$proteins,
  "chromosome|ref=s", \$thischr, # file of ref id list option?
  #....
  "version|MSRC=s", \$MySRC,
  "cadd=s", \@configadd,
  "notbl2asn!", \$notbl2asn, 
  "skipnoann!", \$SKIPNOANN,
  "debug!", \$DEBUG, 
  );
  
$ingff= shift @ARGV unless($ingff);

die "usage: evigene2genbanktbl.pl -conf=evigene.conf  in.gff ...
  opts: -out out.tbl -mid prefixID ... \n"
  unless($optok and $ingff); #  and $action


evigene_config($config, \@configadd); # always even if $config null

$NAME_NONE    = $public_options{nameless} || $NAME_NONE; # same as Uniprot
$NAME_UNK     = $public_options{nameunknown} || $NAME_UNK;  
$NAME_UNKADDLOCUS = $public_options{nameunkaddid} if (defined $public_options{nameunkaddid});  
$NAME_IDPATT= $public_options{nameidpatt} || $NAME_IDPATT;

$MIN_IDENTITY = $public_options{pctunknown} || $MIN_IDENTITY;  # $config{general}->{MIN_PROIDENT}
  $MIN_NAMEIDENT=$MIN_IDENTITY; # for protein_names
$MIN_IDLIKE   = $public_options{pctlike} || $MIN_IDLIKE;   # < pctunknown or >??
$MIN_CERTAIN  = $public_options{pctuncertain} || $MIN_CERTAIN;  # $config{general}->{MIN_PROIDENT}
$MIN_ESTIDENT = $public_options{pctrnaevidence} || $MIN_ESTIDENT;  # $config{general}->{MIN_ESTIDENT} 

our %prot;
$PEPfromGFF = 1 if(($proteins||$evidence{'proteins'}) eq "gff"); #KF was w/o -proteins flag
my $dopepout= ($proteins)?1:($notbl2asn)?0:1; #KF output option
my $nprot= readprot( $proteins||$evidence{'proteins'} ) unless($PEPfromGFF);
# $PutPEPfromFile=1; #??? $PEPfromGFF or $nprot; # =1 for some cases not PEPfromGFF
$PutPEPfromFile= $public_options{putpep} || $config{general}->{putpep} || 0;  

#forget this# $ADDGENE      = (defined $config{attr_gene}->{skip}) ? ! $config{attr_gene}->{skip} : $ADDGENE;
my $DBID= $config{general}->{databaseid} || $config{general}->{dbid} || "MyDBID"; # make global
$MySRC= $DBID; #?? yes or no; now only used for myxref eg paralogs
my $LOCUSTAG= $config{general}->{locus_tag} || $DBID;  

$MyIDPrefix= $public_options{publicid} || $MyIDPrefix;  $MyIDPrefix=~s/\d+$//;
  ##?  Split gene ID format? "id.n"  "id_Cn"  "idSn: : config option?
my $IDSplitSuffix = $public_options{'idsplit_suffix'} || "_C";  
my $IDSplitIsDone = $public_options{'idsplit_is_done'} || 0; # input has already handled splitgene ids
  #?? maybe default yes?
  
$CDS_OFFBY_FIX = $config{general}->{cdsoffbyfix} || $CDS_OFFBY_FIX;   

#      1 ERROR:   SEQ_DESCR.BioSourceMissing          << removed Source info per NCBI request; put back ?
#      1 ERROR:   SEQ_DESCR.NoOrgFound                << ditto
my $TBL_DOSOURCE = (defined $config{general}->{genbanktbl_dosource}) ? $config{general}->{genbanktbl_dosource} : 1; 

my $minpro = $public_options{pctproevidence} || $MIN_IDENTITY; # min for evid if exists 
my $MIN_PROTIDENT= _min($minpro, $MIN_IDLIKE); # for evidence, not naming... maybe use even lower val? 5%, 10%?

our %progvers = ( 
  cuf1  => "Cufflinks:11", 
  cuf8  => "Cufflinks:08", 
  gmap  => "GMAP:11", # default
  gspl  => "Splign:139",
  mbl   => "blastn:22",
  puevd => "PASA:2",
  pasa  => "PASA:2",
  vel   => "Velvet:11" ,
  aug   => "AUGUSTUS:2",  
);
if( my $progvers= $programs{progvers}) {
  my @pg= split/[,\s]+/, $progvers; 
  map{ my($k,$v)= split/=/,$_; $progvers{$k}=$v if($v); } @pg;
}

my $DBXOTHER= $config{general}->{dbxref_other} || "DBXMISSING"; # special code or blank?

my $DBXEST= $config{general}->{dbxref_est} || "dbEST";
my %DBXREF_OK= ( taxon=>1, TAIR=>1, TrEMBL=>1, PGN => 1, ); 
if( my $dbxrefok= $config{general}->{dbxref_ok}) {
  my @pg= split/[,\s]+/, $dbxrefok; 
  map{ my($k,$v)=split/[=:]/,$_; $v=1 unless(defined $v); $DBXREF_OK{$k}=$v; } @pg; 
  #map{ my($k,$v)=split/=/,$_; $DBXREF_OK{$k}=$v if($v); } @pg; # THIS IS WRONG ????!!!
}

my %DBXREF_RECODE= (); 
# dbxref_recode = TAIR=arath Phytozome=poptr|vitvi|soybn|soltu|sorbi
# dbxref_recode = TrEMBL=UniProt  GeneID=mayzebr ENSEMBL=tilapia|zfish|platyfish|medaka|stickleback
if( my $dbxrefre= $config{general}->{dbxref_recode}) {
  my @pg= split/[,\s]+/, $dbxrefre; 
  map{ my($v,$k)=split/[=:]/,$_; my @k=split/[|]/,$k; foreach my $k1 (@k) { $DBXREF_RECODE{$k1}=$v; } } @pg; 
}

our %genescore;
my $addgenescore= $config{general}->{addgenescore} || "ovrna|ovpro";
$inscore = $evidence{genescore} unless($inscore);
if($inscore and open(S,$inscore)) {
  while(<S>){ chomp; my($g,@v)=split"\t"; $genescore{$g}= join";", grep /($addgenescore)=/, @v; } close(S);  
}  

our %changelist;
$changelist= $config{general}->{changelist} unless($changelist);
if($changelist and open(S,$changelist)) {
  # geneid  action:(drop,rename,retype..) newid ... grep valid actions?
  # are there any cases with 2+ actions? need changelist{$g} .= ..
  ## new MisMatchAA = skip my pep, let ncbi translate.
  while(<S>){ next unless(/^\w/); s/#.*$//; chomp; my($g,@v)=split" ";
    warn "#ERR: replaced $g action: $changelist{$g}\n" if($changelist{$g});
    $changelist{$g}= join"\t", @v; } close(S);  
}
  
## add thischr chromosome partition loop here
#..................

our ($inh, $outh, $logh, $peph);
our (%chrhandle,%chrmain,$outref);
my  (@chrlist, @countall);

#Kf: stub for nameclean() assume cleaned names
#   ($newna,$lowqualname,$nadiff)= nameclean( $newna, $pi, $locusadd );
sub nameclean { my($na,$napi,$locadd)=@_;  $na.=$locadd if($locadd); return($na,"",0); }
sub tblsuf { my($na,$pt,$suf)= @_; if($na) { $na =~ s,\.\w+$,,; return join(".",$na,$pt,$suf); } $na; }
sub warno { if(ref $logh){ print $logh @_; } else { warn @_; } }
#  my $c=($_[0] =~ /^#|^\d/) ? shift @_ : 3; $c=($c=~/#/)?$c :($c==1)?"#ERROR: " :($c==2)?"#WARN: " :"#NOTE: "; 

sub MAINstub {}
#..................

if( $evidence{genomesplit} ) { # handle other != chrmain
  my @csplit=  split /[,\s]+/, $evidence{genomesplit}; # now renamed chr
  foreach my $c (@csplit) { next if($c =~ /other/); $chrmain{$c}=1; }
}

$thischr= $evidence{genomesplit} unless($thischr);
#^^ # file of ref id list option?
if($thischr) {
  @chrlist=  split /[,\s]+/, $thischr; 
} else {
  @chrlist= ("all");
}  

# open all output handles per chrlist
# ** FIXME problems with renamed scaffold_10r CHR_RENAME; now require chrlist be renamed, not orig?

%chrhandle=();
foreach my $achr (@chrlist) {
  # handle achr = all, other
  # if($achr eq "other") {  warno "#DEBUG: skipping chr=$achr for now\n"; next; }
    
  my $outgff1= tblsuf( $outgff, $achr, "tbl"); # maybe "stdout|-"
  my $pepout1= tblsuf( $outgff, $achr, "pep"); # whatif stdout ?
  my $logout1= tblsuf( $logout || $outgff, $achr, "gff2tbl.log") if(defined $logout);
  
  if( $outgff1 and -f $outgff1 ) {
    warno "#WARN: already have chr=$achr in $outgff1\n";  next;
  }
  my($outh, $logh, $peph)= (undef) x 3;
  (undef, $outh) = openio( undef, $outgff1); 
  (undef, $logh) = openio( undef, $logout1) if($logout1);
  (undef, $peph) = openio( undef, $pepout1) if($dopepout); #was unless($notbl2asn); #KF : change option for pep output from mRNA prot=
  
  $chrhandle{$achr}{'outh'}= $outh;
  $chrhandle{$achr}{'logh'}= $logh;
  $chrhandle{$achr}{'peph'}= $peph;
}

#? allow 2+ ingff? e.g. genes + transposons, result will be unsorted.tbl
if(%chrhandle) {
  ($inh, undef)= openio($ingff, undef); 
  @countall= filter_gff($inh, "chrlist");
  putsource_noannotchrs() unless($SKIPNOANN); # FIX 2012feb; ** OPT turn this off
  warno "#DONE: gff2tbl chrs: @chrlist, total: @countall \n";
  close($inh);
  foreach my $achr (@chrlist) { setrefout($achr,"close"); }  
}

my $status1= $countall[0] + $countall[2]; # counts = chr, gene, mrna, cds, other

## FIXME: need .fsa, .qvl now; sub putgenome() should handle this.
if(not $notbl2asn) { # $status1 > 0 and 
  %chrhandle=();
  # my $fastachr  = $genome  || $evidence{'genome'};  # die unless(-f $fastachr);
  my $fastaqual = $evidence{'genomequal'} || "$genome.qual";  
  my $havequal= (-f $fastaqual)?1:0; 
  foreach my $achr (@chrlist) {
    # if($achr eq "other") { warno "#DEBUG: skipping chr=$achr for now\n"; next; }
    my $outfsa1= tblsuf( $outgff, $achr, "fsa"); # maybe "stdout|-"
    my $outqvl1= tblsuf( $outgff, $achr, "qvl"); 
    if( $outfsa1 and -f $outfsa1 ) {
      warno "#WARN: already have chr=$achr in $outfsa1\n"; next;
    }
    my(undef, $outh)= openio( undef, $outfsa1); 
    $chrhandle{$achr}{'outh'}= $outh;
    if($havequal) {
    my(undef, $peph)= openio( undef, $outqvl1); 
    $chrhandle{$achr}{'peph'}= $peph; # mess; need declared handle tho
    }
  }
  if(%chrhandle) {
    putgenome(\%config, $outgff,"chrlist");
    foreach my $achr (@chrlist) { setrefout($achr,"close"); }
  }
}

if($status1 > 0 and not $notbl2asn) {
  my $status = tbl2asn( \%config, $outgff); # add thischr option for partitions
  warno "#DONE: tbl2asn: $status \n";
}

#------ subs -------------------------------

sub readprot {
  my($proteins)=@_;
  my $nin=0;
  return $nin if($PEPfromGFF or not defined $proteins);
  # $proteins  ||= $evidence{'proteins'};  

  my($inh, $ok, $id,$aa)= (undef) x 5;
  if($proteins =~ /.gz/) { $ok= open($inh, "gunzip -c $proteins |"); }
  elsif($proteins =~ /stdin|^-/) { $inh= *STDIN; $ok=1; }
  else { $ok= open( $inh, $proteins); } # die "ERROR: read $proteins" unless($ok);
  while(<$inh>) {
    if(/^>(\S+)/){ my $d=$1; if($id and $aa){ $prot{$id}=$aa; $nin++; } $id=$d; $aa="";}
    else{ chomp; $aa.=$_; }
  } close($inh);
  if($id and $aa){ $prot{$id}=$aa; $nin++; } 
  warno "#readprot: $nin from $proteins\n";
  return $nin;
}

sub setrefout
{
  my($ref, $flag)= @_;
  my $ok= 0; $outref="none";
  return 0 if( chrname($ref) eq "skip"); #? here
  # our( $outh, $logh, $peph) all are chrlist dependent
  ## ref == all , ref == other ?
  $ref='all' if($chrhandle{'all'});
  if( %chrhandle ) { # dont change outhandles unless using this hash
    my $hasref=  $chrhandle{$ref} ? 1:0;
    if(not $hasref and $chrhandle{'other'} and not $chrmain{$ref}) {
      $ref='other'; $hasref=1;
    }
    if( $hasref ) {
    $outref= $ref;
    $outh= $chrhandle{$ref}{'outh'}; delete $chrhandle{$ref}{'outh'} if($flag =~ /close/);
    $logh= $chrhandle{$ref}{'logh'}; delete $chrhandle{$ref}{'logh'} if($flag =~ /close/);
    $peph= $chrhandle{$ref}{'peph'}; delete $chrhandle{$ref}{'peph'} if($flag =~ /close/);
    } else {
    $outh= $logh= $peph= undef; #? turn off
    }
  }
  if($flag =~ /close/) {
    close($outh) if(ref $outh);  $outh= undef;
    close($logh) if(ref $logh);  $logh= undef;
    close($peph) if(ref $peph);  $peph= undef;
  }
  $ok= (ref $outh)?1:0; # outh not changed unless %chrhandle or close
  return $ok;
}



sub openio {
  my($ingff, $outgff)= @_;
  my($inh, $outh, $ok)= (undef) x 3;
  
  if(defined $ingff) {
  if($ingff =~ /.gz/) { $ok= open($inh, "gunzip -c $ingff |"); }
  elsif($ingff =~ /stdin|^-/) { $inh= *STDIN; $ok=1; }
  else { $ok= open( $inh, $ingff); }
  die "ERROR: read $ingff" unless($ok);
  }
  
  if(defined $outgff) {
  if(!$outgff or $outgff =~ /^stdout|^-/) { $outh= *STDOUT; $ok=1; }
  else { $ok= open( $outh, ">$outgff"); }
  die "ERROR: write $outgff" unless($ok);
  }
  
  return($inh, $outh);
}

sub setprogram {
  my($evcmd, $defname, $defprog, $nodie)= @_;
  my $progna=(split " ",$evcmd)[0]; 
  my $prog= $defprog;
  if($progna and $prog= $programs{ $progna }) {  }
  elsif($defname and $defprog) { $progna=$defname; $prog= $defprog;  }
  $evcmd =~ s/$progna/$prog/;
  unless( -x $prog) { warn "ERROR: missing program $progna: $prog\n"; die unless($nodie); } # or die
  return $evcmd;
}


our ( %chrsize, %chrattr, %chrrename, %chrhasannot, $didchrrename);

sub chrsize {
  my($chr)= @_;
  $chr= chrname($chr);  
  return 0 if($chr eq "skip");
  
  my $species= $config{general}->{species} || "Unknown";
  my $taxid  = $config{general}->{taxid} || 0;
  
  unless( %chrsize ) {
    $chrsize{$chr}= 0;
    my $chrindex = $evidence{'chrindex'}; # sort chr if this exists
    open(F,$chrindex) or warno "#WARN: missing chrindex=$chrindex\n";
    while(<F>) {
      next unless(/^\w/);
      my($tchr,$w,$attr);
      my @v= split; 
      if( @v > 8 and $v[4]=~/^\d/) { #gff
        ($tchr,$w,$attr)= @v[0,4,8]; 
      } else {
        ($tchr,$w)= @v[0,1];
      }  

      if( $attr =~ /[=,]chr/i and $attr !~ /chromosome=/) { # use only from Name= or Alias=
        my ($chr)= $attr =~ m/[=,]chr[a-z_]+(\d+)/i; # bad for chrI, chrV ..
        unless( $chr) { ($chr)= $attr =~ m/\b(chr\w+)/i; } 
        $attr .= ";chromosome=$chr";
      }
      
      my $rchr= chrname($tchr);
      next if($rchr eq "skip");
      $chrsize{$rchr}= int($w); 
      if($rchr ne $tchr) { 
        # convert all attr tchr to rchr here? problems otherwise for Name=$tchr
        $attr =~ s/$tchr/$rchr/g;
        $attr .= ";Alias=$tchr";  
        }
        
      $chrattr{$rchr}= $attr if($attr);
    } close(F);
  }
  
# 	db_xref	NC_004353
# 	name	NC_004353
# 	organism	Drosophila melanogaster
  my $attr= $chrattr{$chr} || "";  
  $attr .= ";ID=$chr" unless($attr=~/ID=/); # |db_xref=
  $attr .= ";Name=$chr" unless($attr=~/Name=/i);
  $attr .= ";species=$species" unless($attr=~/species=|organism=/);
  $attr .= ";Dbxref=taxon:$taxid" if($taxid and not $attr=~/taxon:/);
  $attr .= ";mol_type=genomic DNA" unless($attr=~/mol_type=/);
  $attr =~ s/^;//;
  
  my $csize= $chrsize{$chr}||0;
  return ($csize, $attr);
}



sub chrname {
  my($r)= @_;
if(CHR_RENAME) {  
  unless( %chrrename ) {
    $chrrename{$r}= $r;
    $didchrrename= 0;
    my $chrrename= $public_options{'chrrename'};
    if($chrrename) {
      my @chrrename= split /[\s,]+/, $chrrename;
      map{ my($k,$v)= split/[=:]/,$_; $chrrename{$k}=$v if($v); } @chrrename;
      $didchrrename=1;
    }
  }
  $r= $chrrename{$r} || $r;  
}  
  return $r;
}

sub writeGenbankHeader
{
	my $self= shift; ## self == config hash
  my($seqid)= @_;

  my $tblname="";
  my $date = $self->{date};
  my $sourcetitle = $self->{sourcetitle};  
  $tblname= join("_",$sourcetitle,$date); $tblname =~ s/\s+/_/g;
  #? add $db  or not; 
  
  print $outh ">Features\t$seqid\t$tblname\n"; 
 # gbtbl: first line is >Features SeqID table_name == chromosome SeqID, same as fasta
}


  #FIXME:  putsource_noannotchrs FIX 2012feb; ** OPT turn this off, or change %chrattr ?
sub putsource_noannotchrs 
{
  foreach my $chr (sort keys %chrattr) {  ## unless($SKIPNOANN)
    next if($chrhasannot{$chr});
    if(setrefout($chr)) { # output handle for chrsplitting; checks chrhandle valid chr
      putsource($chr);
    }
  }
}



sub reformatProteinID { my $t=shift; $t=~s/t(\d+)$/p$1/; return $t; } # called before adding split suffix 

sub reformat_myid
{
  my($id,$formopt)= @_;
  my $db= $DBID;  $formopt||="";
  unless($db) { ($db)= $id =~ m/^(\w+[A-Za-z])\d\d/; }
  $db= 'gnl|'.$db.'|';
  return $id if($id=~m/^gnl\|/);
  
  ## splitgene: FIXME here?? add hash lookup for few splitgene IDs; protein_id|transcript_id need split part-id tags, not gene
  if(SPLITGENE_FIX) { 
    # (my $tid=$id) =~ s/p(\d+)/t$1/; # not needed, splitgene{} has both tid, pid
    # Add option that caller has already added splitgene ID tags, dont change here..
    # Add param here to pass pn/splitgene num of id, for linke other split half
    my $pn= $splitgene{$id}; # || $splitgene{$tid};
    if($formopt =~ /Split=(\d+)/) { $pn=$1; }
    if($pn and $id !~ m/$IDSplitSuffix/) {
      return $db.$id.$IDSplitSuffix.$pn;  
      ##? format? "id.n"  "id_Cn"  "idSn: : config option now
    }
  }
  return $db.$id;  # prefered in ncbi doc #? return "lcl|$id";  
}


sub reformatval
{
  my($typ, $id, $key, $oldkey, $val)= @_;
  # my $rekey=$key;
  my @val; #?? = undef; << BAD

  my $QUALDROPCLASS=1; # qual note fixes, option?
 
  if($key =~ /^(locus_tag)$/) {
    # FIXME15: locus_tag needs splitgene tag also BUT cant use same _Cnn tag !!
    # my $pn= $splitgene{$id}; if($pn) { return $db.$id.$IDSplitSuffix.$pn;  }
    $val =~ s/^\w*([A-Za-z])(\d\d+)/${LOCUSTAG}_${2}/;  ## use config; ncbi registered locus tag = TCM 
    
  } elsif($key =~ /^(gene|protein_id|transcript_id)$/) {
    $val= reformat_myid($val);
    
  } elsif($key =~ /^product$/) { 
    $val =~ s/\s*\((\d+)%.?\)//; # trailing pct ident; may be removed.. (73%P)
  }
  ## FIXME: product for Unknown/Unchar should change to ID! And/or do for all Names that are not uniq?
  ## add isoform to product :: fix in processgene() below

#  ## fixme: dont add chromosome here, see chrsize()
#   if($typ eq "source" and $key eq "name" and $val =~ /chr/i) { # fixme: use chrattr table
#     @val=("$key\t$val"); # keep name, add /chrom
#     if( my ($chr)= $val =~ m/chr[a-z_](\d+)/i) { push @val, "chromosome\t$chr"; } 
#     else {  push @val, "chromosome\t$val"; }
#   }  

    
# for CDS, quality attr:	/pseudo  /exception="xxx"  /note="xxx" for Protein=partial ..
#  Gene fragments: as gene with tags: note=nonfunctional due to frameshift ; note=frameshift ; /pseudo ; 
#   or as misc_feature with note=xxx
 
  ## old paralog:    
  #     my ($did)= $val =~ m/($MyIDPrefix\w+)/; # can be many...
  #     return () unless($did and $pi >= $MIN_IDENTITY); 
  #     return () if ($changelist{$did} =~ /^drop/);
  #     # check if $changelist{$tid} =~ /^drop/ 
  #     ## return () if(not $did or $val =~ /paralog.na|^na|paralog.\d|^\d/ or $pi < $MIN_IDENTITY); 
  #     # note    paralog=Thecc1EG034652t1,83% < drop % .. id:pct
  #     # fixme:    note    gnl|CacaoGD|paralog=Thecc1EG005594t1
  #     # fixme2:   note    paralog  Thecc1EG034062t1 gnl|CacaoGD|Evigene.of
  #     # fixme3:   note    paralog of 47% gnl|CacaoGD|
  #     
  #     if($key =~ /^note/) {    # only for CDS:note,db_xref but for mRNA:inference, do like ortholog
  #       $val =~ s/(\d+)%$//;   # drop :pid # $val =~ s/%//; $val =~ s/,/:/; 
  #       $val =~ s/\s*$did[,]?//; $val=~s/,\s*$//; 
  #       @val=("$key\t$val " . reformat_myid($did));
  #       ### fixme: recode as comment not db_xref
  #       # TURN OFF for now#  push @val, "db_xref\t$DBXOTHER:$did" unless($DBXOTHER =~ /MISSING/);
  #     }

  if($oldkey =~ /paralog/) { 
    my ($pi)= ( $val =~ s/(\d+)%// ) ? $1 : 69;  # also drop :pid
    return() unless($pi >= $MIN_IDENTITY);
    my @did= $val =~ m/($MyIDPrefix\w+)/g;
    @did= grep{ not $changelist{$_} =~ /^drop/ } @did;
    return () unless(@did);
    if($key =~ /^note/) {    # only for CDS:note,db_xref but for mRNA:inference, do like ortholog
      map{ $val =~ s/$_[\s,]*// } @did; $val=~s/\W+$//; 
      my $did= join", ", map{ reformat_myid($_) } @did;
      @val=("$key\t$val $did");
    }

  } elsif($oldkey =~ /ortholog/) {  # new db_xref
    my($pi)= ( $val =~ m/(\d+)%/ ) ? $1 : 69; 
    return () if($val =~ /^na|^\d/ or $pi < $MIN_PROTIDENT); 
    
  } elsif($oldkey =~ /isoform/) {
     $val .= "; alternatively spliced"; # always?
     
  } elsif($oldkey =~ /express/) {
    my($pi)= $val =~ m/(\d+)%/;  
    return () if($pi < $MIN_PROTIDENT); 
    ##  infer no good, need genbank-valid dbxref IDs 
    ## .. format opts; ncbi annot uses 'xx%' coverage syntax
    #o $pi= ($pi>99)?"1.00":"0.$pi"; ## sprintf "%.2f", $pi/100; # change from % to prop; 1.0 ..
    $pi= ($pi>99)?"100%":"$pi%"; 
    $val= "similar to RNA sequence, $pi coverage"; #as note not inference; gb hates %%% NO, is inconsistent
    ## return ("note\t$val");
    
  } elsif($oldkey =~ /Split/) { 
    my($sn,$nsplit)= isSplit("Split=".$val);
    unless($sn){ ($sn)= ( $val =~ m/(\d+)/ ) ? $1 : 0; }
    my($on)= ($sn>1)?1:2;   # need multi:  Split=2/4 ..
    # if($nsplit){ } # what?
    my $idb=$id; $idb=~s/$IDSplitSuffix.*$//;
    if($on) { 
      my $tpd= reformat_myid( $idb, "Split=$sn"); 
      my $opd= reformat_myid( $idb, "Split=$on"); 
      $val .= ", this $tpd link to $opd"; 
      }

  } elsif($oldkey =~ /quality/) {  # reformat for genbank
     ## quality=Class:Strong,Express:Strong,Homology:ParalogStrong,Intron:Strong,Protein:complete
     ## >> note    quality=Class:Strong,Express:Weak,Homology:ParalogStrong,Intron:None,Protein:complete
     #nogood, other parts to qual: my $newval="";
     while($val =~ m/(\w+):([\w-]+)/g) { 
      my($k,$v)=($1,$2);
      if($k eq 'Class' and $QUALDROPCLASS){ $val =~ s/$k:$v[,\s]*//; } 
      else { 
        my($k1,$v1)=($k,$v);
        ## $v1 chop off all '-XXX' extras? confusing/maybe wrong
        # protein is unavailable; protein is validated-None << bad
        # map is poor-Split << maybe wrong
        # intron is none-inerr << confused
        $v1 =~ s/(\-Split|-\inerr)//ig;
        $v1 =~ s/validated-//i if($k eq 'Protein'); # validated-partial5 .. not useful here
        $v1='' if($k eq 'Protein' and $v =~/None|unavailable/i);
        
        $k1  =~ s/\bexpress\b/expression/i;
        $v1  =~ s/(\w)(Strong|Medium|Weak|Poor)/$1 $2/;  
        #nogo# $newval .= "\L$k1 is \L$v1, "; 
        if($v1) { $val =~ s/$k:$v/\L$k1 is \L$v1/; } else { $val =~ s/$k:$v[,\s]*//;  }
        }
      } 
     #nogood# $newval=~s/,\s*$//; $val=$newval; # 
     $val =~ s/,/, /g; 
  }

  # dbxref: split "," into @val ?
  if($key =~ /db_xref/) { 
    return () if($val =~ /^na|^\d/); 
    if($oldkey =~ /ID|paralog/) { $val=  reformat_myid($val); } # others? check $val for MYIDpatt ?
    # ^ do below in id map{}
    #... Need to use ncbi db_xref.html where possible.  read from config?
    if($oldkey =~ /Dbxref/) { @val= split",",$val; @val= map { s,/.*,,;  $_ }  grep /^[a-zA-Z]/, @val; }
    else { $val=~s/,.*$//; @val=($val);  }
    my $dbskip= $config{general}->{dbxref_skip} ||"";  
    if($dbskip) { $dbskip =~ s/[\s,]+/|/g; @val = grep { not m/^($dbskip)/ } @val; }
    
    ## need to deal w/ any not in db_xref.html (most) : PGN is closest plant generic db
    # dbxref_ok  taxon DDBJ EMBL NCBI GeneID GI dbEST dbSNP TrEMBL Swiss-Prot  UNILIB InterPro ENSEMBL GO 
    #   + JGIDB ESTLIB APHIDBASE AntWeb BEETLEBASE dictyBase FLYBASE GDB GeneDB GRIN MGI PDB PFAM PGN SGD SGN TAIR VectorBase WormBase ZFIN
    # valid:  ncbi201107/api/asn2gnb6.c superset of http://www.ncbi.nlm.nih.gov/collab/db_xref.html
    # ** and NCBI is not a valid db; none for NCBI ids listed ..
    ## see also legalRefSeqDbXrefs in  api/asn2gnb6.c : NASONIABASE 
    
    my $NCBIdb="GeneID"; # "LocusID";#? PIR? LocusID>GeneID, ok
    @val= grep /\w/, map { 
      my $origval=$_;
      s/^UniRef/TrEMBL:UniRef/;  
      if(/^[A-Z0-9]+[_][A-Z0-9]{4,9}$/){ s/^/TrEMBL:/; } #catch all uniprot bare ids: [A-Z0-9]+_[A-Z]{6}
      # E2AL44_9HYME
      
      s/^(tr|pr).([XN][MRP]_)/$NCBIdb:$2/;
      s/^(apis2ref).([XN][MRP]_)/$NCBIdb:$2/;
      #   db_xref PGN:tr.XM_001600948.2 PGN:pr.XP_001600998.2  PGN:GeneID.100116528 # ncbi ids
      #   NcbiRef2rna2410 = localid; skip
      #   apis2ref.XP_003251084.1 = NCBI
      #bug:   db_xref PGN:altof.NcbiRef2rna2410 db_xref PGN:t2

      # plant genome tags: move to DBXREF_RECODE config.
      s/^(arath):/TAIR:/;
      s/^(poptr|vitvi|soybn|soltu|sorbi):/Phytozome:/; # fixme recode; frave, ricco are not Phytozome.. dig up other?
      ## s/^(arath|poptr|vitvi|frave|ricco|soybn|soltu|sorbi):/$DBXOTHER:\u${1}./; # hack to stop whining; PGN plant genome db in ncbi/db_xref.html
      
      my($d)= m/^(\w+):/; 
      if( $DBXREF_RECODE{$d} ) {  #KF FIXME, change $_ also
        my $rd= $DBXREF_RECODE{$d}; s/$d:/$rd:/; $d=$rd;
        } 
      unless($d) { $_= "$DBXOTHER:$_"; }
      elsif(not $DBXREF_OK{$d}) { s/$d:/$DBXOTHER:$d./; } # need config
      ## need option to move DBXOTHER = Missing to comment not db_xref
      if(/MISSING:/) { $_="note\tdbxref $origval"; } 
      else { $_="$key\t$_";  }
      } @val;
      
    return @val; # may be empty

  } elsif($key =~ /inference/) { 
    @val= map{ "$key\t$_" } evidencecode($key, $oldkey, $val);
    return @val; # may be empty
  }
  
  @val=("$key\t$val") unless(@val);
  return @val; # add rekey
}


sub evidencecode
{
  my($key, $oldkey, $val)= @_;

## dont get this error but reformat to make it go away: ONLY for alignment TSA:xxx
# [SEQ_FEAT.InvalidInferenceValue] Inference qualifier problem - spaces in inference 
# (alignment:Cufflinks:08:TSA:cacao11r39cuf8_Gsc8g3847t3:97.83)  << change last : to _97
# :$pi > _$pi doesnt solve 'spaces in inference', try '-'
# ... drop the trailing :pi for genbank ...

  my $DROPPI=1; 
  my $PSEP='i'; #no? '-'; no '_'
  
## add genepred:
#   /inference="ab initio prediction:Genscan:2.0"
#   /inference="ab initio prediction with evidence:AUGUSTUS:2"
  if($oldkey =~ /prediction/) {
    my @v = grep /^[a-zA-Z]/, split",", $val;
    my @pk= sort keys %progvers;
    my $st="ab initio prediction";
    my @ev= grep /\w/, map {
      my $vv=$_;
      my $prog= $progvers{'aug'}||"AUGUSTUS";
      foreach my $pk (@pk) { $prog= $progvers{$pk} if($vv =~ m/$pk/); }
      my $sv=$st; $sv.=" with evidence" if($prog=~/AUG/); # fixme
      "$sv:$prog" if($prog);   
    } @v;
    return @ev;
  }
  
  # recode idprefix: r8caca11 > cacao11; m7vitvi> vitvi
  # leave out '%' char
  # add key insrc=xxx:yyy insrc=kf2a:gmap2a5u; ? should say map progver : gspl/gmap/other
  if($oldkey =~ /ovrna|trasm/) {
    my @v = grep  /^[a-zA-Z]/, split",", $val;
      ## fixme: CacaoGD mar11f:CGD are bad ovrna: alignment:GMAP:11:TSA:mar1g.mar11f:CGD0022153
      ## bad also:  alignment:GMAP:11:TSA:mar7g.mar11f:AUGepir7p1s1g7t1
      ## bad also:  alignment:GMAP:11:TSA:AUGpie3_p1s_1g7t1
    @v= grep { ! /mar1g.mar11f|mar7g.mar11f|^AUG/ } @v;
    my @pk= sort keys %progvers;
    my @ev= grep /\w/, map {
      my $vv=$_;
      
      s/^\w+://; # keep this
      #o s/^(cgba).\w+:/$1/; 
      #o s/^r8//; s/caca11/cacao11/; # add to config; 2015 drop these project tags?
      #o $_= "cgba".$_ if(/^(B|L|P1|P2)_g\d/);
      ## TSA:cgba.bean:B_g00887t00002 << cgba_B_g00
      ## TSA:B_g00887t00002   << recode these all as cgba_B_g00
      ## TSA:cgba.leaf:L_g00772t00003
      
      my $pi=99;  # simplify 88.77 to pimax(cds,exon) ?
      if(s,/[CI]?(\S+),$PSEP$1,){ $pi=$1; } 
      elsif($oldkey =~ /trasm/) { $pi=100; $_.=$PSEP."100"; }  #$_.="%"; 
      s/$PSEP$pi// if ($DROPPI);
      my $vo=$_;
      my $prog= $progvers{'gmap'}||"GMAP";
      foreach my $pk (@pk) { $prog= $progvers{$pk} if($vv =~ m/$pk/ or $vo =~ m/$pk/); }
      "alignment:$prog:TSA:$vo" if($pi >= $MIN_ESTIDENT);   
    } @v;
    # ovrna=100,75,r8B_g00334t00001/I100  >>   "alignment:GMAP:TSA:$tid";
    return @ev;
  }

  elsif($oldkey =~ /inference/) {   
    return ($val);
  }

    #  inference       0%,eq:0 << bad    
  elsif($oldkey =~ /express/) {  # need some EST ids here.. redo as note, not inference
    my($pi)= $val =~ m/(\d+)%/;  
    ##my $estdb="CGBest"; # config DBXEST
    ## this inf no good, need genbank-valid dbxref IDs 
    my $v= "similar to RNA sequence, EST:$DBXEST:$pi"; #% need some est-db config
    ## my $v= "similar to RNA sequence, EST"; #% this fails muster
    # my $v= "alignment:GMAP:EST:$p"; or this?
    return ($pi >= $MIN_ESTIDENT) ? ($v) : ();
  }

  ## fixme: remap protdb to new prot db:
  ##   similar to AA sequence:poptr:POPTR_0007s08080.1  >> Phytozome:POPTR_0007s08080.1
    
  elsif($oldkey =~ /ovpro/) { #? use homolog= score also/instead?
    # ovpro=98,98,m7ricco:29682.m000590/C98.00
    my @v= grep  /^[a-zA-Z]/, split",", $val;  # add pi score here also? should be 1st num
    # my @ev= map { s,/.*,,; s/^m7//; "similar to AA sequence:$_"; } @v;
    # alt form: alignment:exonerate.2:$protid:$pi
    # my $prog="protein-exonerate:2"; # should be 
    my $prog= $progvers{'exonerate_prot'}||"protein-exonerate:2";

    my @ev= grep /\w/, map { 
      my $pi=69; if(s,/[CI]?(\d+).*,$PSEP$1,){ $pi=$1; s/$PSEP$pi// if ($DROPPI); } 
      s,/.*,,; s/^m7//; 
      my($db)= m/(\w+):/; if( $db and (my $dbr= $DBXREF_RECODE{$db}) ) { s/$db:/$dbr:/; }
      "alignment:$prog:$_" if($pi >= $MIN_PROTIDENT); 
      } @v;
    return @ev;
  }
  
  elsif($oldkey =~ /ortholog|paralog|homolog/) {  
    # homolog=441/581,poptr:POPTR_0011s14980.1
    # ortholog=soltu:PGSC0003DMP400002993,71%;
    # .. for paralog, check if changelist drop
    if($oldkey =~ /paralog/) {
      ##bad# my ($did)= $val =~ m/\b(\w[^\s,;:]+)/; 
      my ($did)= $val =~ m/($MyIDPrefix\w+)/;
      return () if ($changelist{$did} =~ /^drop/);
    }
    my $ispara=($oldkey =~ /paralog/)?1:0;
    my $ss= ($ispara)?" (same species)":"";
    my $pi= 69;
    if( $val =~ m,(\d+)%, ) { $pi=$1; }
    elsif($val =~ m,(\d+)/(\d+),) { my($bt,$bs)= ($1,$2); $pi= int(0.5 + 100*$bt/$bs) if($bs>0); }
    my @v= grep  /^[a-zA-Z]/, split",", $val;  # add pi score here also? should be 1st num
  ## fixme para format bug:
  ##  fail:  similar to AA sequence (same species):Thecc1EG016762t1 << needs :CacaoGD:Thecc 
  ##  oknew: similar to AA sequence:Phytozome:PGSC0003DMP400024258
    @v= map {
      s,/.*,,;  
      my($db)= m/(\w+):/; 
      if( $db and (my $dbr= $DBXREF_RECODE{$db}) ) { s/$db:/$dbr:/; }
      elsif($ispara) { s/($MyIDPrefix)/$DBID:$1/; } 
      "similar to AA sequence$ss:$_"; } @v; # need :$pi here for dang format
    # unless $DROPPI: @v= map { s,/.*,,;  "similar to AA sequence$ss:$_".":$pi"; } @v; # need :$pi here for dang format
    return @v;
  }
  
}


=item evidencecode

  inference=[CATEGORY:]TYPE[ (same species)][:EVIDENCE_BASIS]
  cat= COORDINATES or EXISTENCE
  type= similar to AA sequence , similar to RNA sequence, "ab initio prediction", "alignment" ..
        "non-experimental evidence, no additional details recorded"
  evid= dbxref,dbxref2,...
  ##  inference="similar to RNA sequence, mRNA:RefSeq:NM_000041.2"
  ##  inference="alignment:Splign:1.26p:RefSeq:NM_000041.2,INSD:BC003557.1"
  ##  inference="ab initio prediction:Genscan:2.0"
  ** use TSA? for transcript assembly dbxref ?
  TSA	shotgun assemblies of primary transcriptome data deposited in INSDC, the Trace Archive (TA) or the Short-Read Archive (SRA) 
   Transcriptome Shotgun Assembly 
    inference  "alignment:GMAP:TSA:xxxx"

=cut


sub validattr
{
  my($typ, $id, $attr)= @_;
  
  my $mapattr= $config{"attr_".$typ};
  # need retype option, after config: $gtype= $config{"attr_".$typ}->{type} || $typ;
  my $typo= $mapattr->{type} || $typ;

  if($typo eq "mobile_element") { # typ eq "transposon"
    # here? turn Name,class attr to mobile_element_type=classn:name
    my($cl)= $attr=~/class=([^;\n]+)/;
    my($na)= $attr=~/Name=([^;\n]+)/;
    my $newna=$na;
    
    if($na =~ /^(transposon|retrotransposon|non-LTR retrotransposon|integron)/) { 
      # LINE|SINE|MITE < valid vocab here but push down to trans/retro key
      my $met=$1; unless($na=~/^$met:/) { $newna="$met:$na"; }
    } elsif($cl) {
      if($cl =~ /^II\b/i) { $newna="transposon:$na"; }
      elsif($cl =~ /^I\b/i) { $newna="retrotransposon:$na"; }
    } else {
      if($na =~ /^TIR|Helitron|MITE/) { $newna="transposon:$na"; }  # class=II
      elsif($na =~ /^LTR|^LINE/) { $newna="retrotransposon:$na"; } # class=I
      else { $newna="other:$na"; } # punt? or other:?
    }
    $attr =~ s'Name=$na'Name=$newna' if($na ne $newna);
  }
  
  my ($issplit,$nsplit)= isSplit($attr); 
  $issplit ||= $splitgene{$id}; #??
  
  my @attr; my %didk;
  # 201504 patch for locus_tag/gene; new input attr locustag= replaces gene= for NCBI stupidness
  my $hasloctag= ($attr=~/locustag=/)?1:0;
  foreach my $at (split";",$attr) {
    my($k,$v)= split"=",$at;     
    
    #? if($hasloctag and ($k eq 'gene' or ($typo eq 'gene' and $k eq 'ID'))) { $didk{$k}++; next; }
    next if($didk{$k}++ and $k !~ m/Note|Dbxref|Alias/); #?? not always ! eg Note
    
    my $kn= $mapattr->{$k} || "skip";
    next if($kn =~ /^skip/); # make mapattr full set of allowed tags?
    next if($hasloctag and ($kn eq 'locus_tag' and $k ne 'locustag')); #??
    
    my $kval="";
    ($kn,$kval)= split(/[=\s]+/, $kn, 2); 
    
    if($kval) { $v= ($kval=~m/[:=]$/)?"$kval$v":"$kval $v"; } # more?
    # isoform   note=encoded by transcript variant $v
    
    ## ?? need this patch for all ids?
    # locus_tag|gene|protein_id|transcript_id << NO GOOD locus_tag what now?
    ## DAng them only allow A,B,C sort of split part suffix
    # FATAL: DiscRep_ALL:BAD_LOCUS_TAG_FORMAT::12 locus tags are incorrectly formatted.
    # kfish2gspln15nfcds5_knset8:Gene D326_017569_C1  lcl|KN805979.1:182188->182354   D326_017569_C1
    # kfish2gspln15nfcds5_knset8:Gene D326_017569_C2  lcl|KN807334.1:<683348-854416   D326_017569_C2

    if($kn eq "locus_tag" and $issplit > 0) {
      ## check if $v ends in A-Z first; may need to prepare gene ids this format for splits w/ alts same idtag/scaf
      ## NO, drop this ; mistake now..
      # unless($v=~m/[A-Z]$/) {
      #  my $lts= substr('ABCDEFGHIJKLMNOPQRSTUVWXYZ',$issplit-1,1); # dang u ncbi
      #  $v .= $lts if($lts); 
      # }
    }
    elsif($kn =~ /^(gene|protein_id|transcript_id)$/ and $issplit and not $v=~m/$IDSplitSuffix/) {
      $v .= $IDSplitSuffix.$issplit; #? ok for protid?
    }
    
    my @keyval= reformatval( $typo, $id, $kn, $k, $v);
    push @attr, grep /\w/, @keyval;    
    
    #NO?# if($kn eq "locus_tag" and $typ eq "gene") {
      #NO?# my @kv= reformatval( $typ, $id, "gene", $k, $v);
      #NO?# push @attr, $kv[0] ;
    #NO?# }
  }
  return @attr;
} 


sub putattr
{
  my($typ, $id, $attr)= @_;
  ## here or below, add id as Dbxref=Evigene:$id
  $attr .= ";myxref=$MySRC:$id" if($id =~ /$MyIDPrefix/ and not $attr =~ /=$MySRC/); # check if is my id: MyIDPrefix
  my @attr= validattr($typ, $id, $attr); # attr may have dups
 
  foreach my $at (@attr) {
    my ($k,$v)= split "\t",$at,2; 
    #was# print $outh join("\t","",$k,$v),"\n";
    print $outh join("\t","","","",$k,$v),"\n"; # empty cols 1-3 for attrs
  }
  print $outh "\n";
}


sub putloc
{
  my($typ, @loc)= @_;
  putCDSloc($typ, 0, @loc);
}

sub putCDSloc
{
  my($typ, $partial, @loc)= @_;
  my $splitpart= ($partial =~ s/\WSplit=(\S+)//i)? $1 : 0; #? wedge onto partial flag? eg.: complete,Split=1/3
  
  if($partial eq "partial") { $partial.="53"; } # both
  # need retype option: if($config{attr_$typ}{type})
  my $typo= $config{"attr_".$typ}->{type} || $typ;
  my $addattr=""; # return to caller
  
  my $gstrand= $loc[0]->[6]; ## my(undef,undef,$gstrand)= split("\t",$loc[0]);
  my $isrev=($gstrand eq "-");
  my @iter= ( 0 .. $#loc );
  @iter= reverse @iter if ($isrev);
  my($p5,$p3)= ("<",">"); #?? or this# ($gstrand eq "-") ? (">","<") : ("<",">");

  ## FIXME 1504: issplit *parts* need partial loc syntax (but protein can be complete)
  ## add opt $splitpart .. strand check??
  # "1/3" == $p3 first/start is not partial, 
  # "2/3" == $p5/p3 both ends partial, 
  # "3/3" == p5 last/end not partial
  my($split5,$split3)=(0,0);
  if($splitpart =~ /(\d+).(\d+)/){ my($p,$n)=($1,$2); 
    #? this is local rev, not all parts# if($isrev) { $split3=1 if($p==1); $split5=1 if($p==$n); } else 
    $split5=1 unless($p==1); $split3=1 unless($p==$n); 
  }
  
  my $first= $iter[0];
  my $last = $iter[-1];
  foreach my $i ( @iter ) { 
    my($start,$stop,$strand,$phase)= @{$loc[$i]}[3,4,6,7];  # split("\t",$loc[$i]);
    ($start,$stop) = ($stop,$start) if ($strand eq "-"); ## ($gstrand < 0);
    if(($partial =~ /3/ or $split3) and $i==$last ) { $stop="$p3$stop"; }
    if(($partial =~ /5/ or $split5) and $i==$first) {  
      $start="$p5$start"; 
      if($typ =~ /CDS/) { $addattr .=";codon_start=".($phase+1); }
    }
    
    my @v= ($i==$first) ? ($start,$stop,$typo) : ($start,$stop);
    print $outh join("\t",@v),"\n"; 
  }
  
  return ($addattr);
}


our $lastref;

sub putsource
{
  my($ref)= @_;
  if($ref ne $lastref) {  # should call from putloc(); see other above
    my ($chrsize, $chrattr)= chrsize($ref);
    my $type="source";
    writeGenbankHeader( $config{general}, $ref);  
    # option: dont write source info ...
    if($TBL_DOSOURCE) {
    putloc($type, [$ref,0,$type,1,$chrsize,0,"+",]);  
    putattr($type, $ref, $chrattr);
    }
    $lastref=$ref;
    $chrhasannot{$ref}++;
    return 1;
  } else {
    return 0;
  }  
}


sub putprot
{
  my($pid,$prot)= @_;  # $pout,
  return 0 unless($pid and $prot);
  $prot =~ s/(.{$FAWIDTH})/$1\n/g;
  $prot .= "\n" unless($prot =~ m/\n$/);
  $pid= reformat_myid( reformatProteinID( $pid)); # do in caller, with Split=1/2..
  print $peph ">$pid\n",$prot if(ref $peph);
  return 1;
}  
 


=item nameclean
  
  2013Mar: replace w/ evigene/scripts/prot/protein_names.pm : sub nameclean()
  
  see also longer method in evigene/scripts/bestgenes_puban_cacao.pl : nameclean()

#NOTE rename Thecc1EG026905t1: 'Uncharacterized protein locus Thecc1EG026905' < 'Prefoldin chaperone subunit family protein (19%T)'
#NOTE rename Thecc1EG026905t2: 'Uncharacterized protein locus Thecc1EG026905 isoform 2' < 'Prefoldin chaperone subunit family protein (24%T)'
    # -- should these be putative instead of Unk? or -like? or family?
    # Proteins of unknown function which exhibit significant sequence similarity
    # to a defined protein family have been named in accordance with other members
    # of that family... e.g. "Holliday junction resolvase family endonuclease".
    # It is also possible to use "-like" in the name, for outliers

See update in evigene/scripts/namecleangb.pl
   -- use for genes.gff naming before here;
   -- fix below to use above as package sub?
   
=cut

# sub nameclean_OLD
# {
#   my($namin, $pi, $locusname)= @_;
#   my $lowqualname=undef;
#   my ($isunk,$isput,$islike)= (0,0,0);
# 
#   #*** ?? IS THIS FAILING .. missing haem, Tumour
#   # local $_= $namin; 
#   #?? study();
#   my $_ = $namin;   # perl 5.9.1
# 
#   s/\s+\(\d+%.*\)//; # trailing pctident
#   s/\s+sym:\S+//;  # trailing symbol tags
#   
#   #  horribles die quickly...
#   if(/The protein encoded by this gene was identified/) { $_=$NAME_UNK; }
#   elsif(/Encodes a close of the Cauliflower OR/) { $_="Orange protein"; }
#   
#   # typos and britglish > amglish; need config list
#   s/Uncharacterised/Uncharacterized/; 
#   s/dimerisation/dimerization/;  s/luminium/luminum/g;  # Aluminium
#   s/signalling/signaling/; # Two Ls in British English, one in American English. 
#   s/onoxygenase/onooxygenase/; # [Mm]onox..
#   s/sulphide/sulfide/ig; s/sulphur/sulfur/ig;
#   s/Tumour/tumor/ig;  
#   s/haemoprotein/hemoprotein/i;  s/\bhaem/heme/ig; # Quinohaemoprotein>Quinohemoprotein
#   s/characteris/characters/;  
#   s/proteine/protein/;  s/\bcomponenet/\bcomponent/;
#   s/ protei$/ protein/; # prior nameclean bug
#   
# #  # UPDATE: add to species.config some of these nameclean patterns
# #   my $nameclean = ($config{nameclean} and ref($config{nameclean})) ? $config{nameclean} : {};
# #   foreach my $nk (sort keys %$nameclean) {
# #     my $npatt= $nameclean->{$nk};
# #     # if($npatt =~ m/^s(.)/) { my $nc=$1; my($ncut,$nto)= split/[$nc]/,$npatt; s/$ncut/$nto/; }
# #     my $eok= eval $npatt; # is this ok? prefer s/$nacut/$nato/ ; 
# #     # any use for m/$napatt/ ; need specific keys for this, m/$ISGENEID/ =  m/^(Os|At|AT)\d{1,2}g\d{3,}/
# #   }
#   
#   if(s/^TE://i) { }  # $istransposon=1; 
#   if(s/^(Predicted|Conserved|Expressed|Novel)\s+//i) { }  # maybe set $isunk?
#   if(s/^(Possible|potential|probable)\s+//i) { $isput=1; }   
#   s/(homology|similar|Similarity|related) to\s*//i;
#   s/\s\(Fragment\)//; s/\s[Ii]soform\s*.*$//;
#   unless(m/protein\-protein/) { s/^(ORF|protein)\s*//i; }
#   
#   # no-no species names: Arabidopsis yeast  human
#   # >> staphylococcal? = Staphylococcal nuclease ue, TAIR name
#   # and: genome complete  pseudogene? = TAIR name
#   # my $namedrops= $public_options{namedrops} || 'Arabidopsis|thaliana|yeast|complete sequence|complete|genome|pseudogene'; #plants;
#   # s/\b($namedrops)[,\.\s]*//ig;
#   s/\b(Arabidopsis|thaliana|yeast|human|Staphylococcal|complete sequence|complete|genome|pseudogene)[,\.\s]*//ig;
#   s/\s*\((?:InterPro|TAIR):[\w\.-]+\)//ig; #  (TAIR:AT1G22000.1); (InterPro:IPR010678)
#   if(s/paralog of //i) { $islike=1; }
#   
#   # horrible names:
#   # hydrolases, acting on acid anhydrides, in phosphorus-containing anhydrides,ATP-dependent helicases,nucleic acid binding,ATP bi...
#   # tRNA (guanine-N(1)-)-methyltransferase, metazoa , tRNA (guanine-N1-)-methyltransferase, eukaryotic == duplicated
#   # mannose-1-phosphate guanylyltransferase (GDP)s,GDP-galactose:mannose-1-phosphate guanylyltransferases,GDP-galactose:glucose-1-
#   # ATP-dependent peptidases,nucleotide binding,serine-type endopeptidases,DNA helicases,ATP binding,damaged DNA binding,nucleosid>
#   # serine/threonine kinases,protein kinases,ATP binding,sugar binding,kinases,carbohydrate binding 
#   # "The protein encoded by this gene was identified as a part of pollen proteome by mass spec analysis, It has weak LEA proteins, Encodes protein phosphatase 2A B'gamma subunit, Targeted to nucleus and cytosol"
#   
#   if(length($_) > 70) {
#      my $nc= tr/[.,]/[.,]/; 
#      if($nc > 1) { # keep only 1st two phrases.
#      my $i= _min(index($_,',',0),index($_,'.',0));  
#      while($i>0 and $i<30) { $i= _min(index($_,',',$i+1),index($_,'.',$i+1)); }
#      $_= substr($_,0,$i) if($i>20);
#      } 
#   }       
# #   if(length($_) > 70) {
# #      my $nc= tr/,/,/; 
# #      if($nc > 2) { # keep only 1st two phrases.
# #      my $i= index($_,",",0);  $i= index($_,",",$i+1);
# #      $_= substr($_,0,$i) if($i>20);
# #      } 
# #   }   
# 
#   ## should drop these geneid == name cases: At3g18210 for arabid genes, Os12g0628400, ...
#   ## replace w/ special UNK name? Unchar prot gene id
#   my $isid=0;
#   
#   # my $nameidpatt= $public_options{nameidpatt} || '(?:Os|At|AT)\d{1,2}g\d{3,}|GLEAN_\d+|G[A-Z]\d\d+'; #plants;
#   my $nameidpatt= $NAME_IDPATT; ## $public_options{nameidpatt} || $NAME_IDPATT;
#   if (/^($nameidpatt)/) { $isid=$1; $_= $NAME_UNK; } #? move to species.config
#   elsif(s/($nameidpatt)//) { $isid=$1; } # $hasid=$1;?
#     
#   $isunk= ($pi < $MIN_NAMEIDENT or  m/^($NAME_NONE)/i)?1:$isunk; # Unknown|Uncharacterized|Hypothetical
#   $isput= (!$isunk and ($pi < $MIN_CERTAIN or $isput))?1:0;
#   
#   ## ugh: TRIGALACTOSYLDIACYLGLYCEROL 1,
#   my $nuc= tr/[A-Z]/[A-Z]/;  #  uc($_) eq $_
#   if(m/[A-Z]\s+[A-Z0-9]/ and ($nuc > 19 or uc($_) eq $_ ) ) { $_= ucfirst(lc($_)); } # SHOUTING phrase ...  TRICHOME BIREFRINGENCE-LIKE 19
#   
#   s/\s*[Hh]omolo(gy|gue|g)\s+\d+//g; s/\b[Hh]omolo(gy|gue|g)\s*//g; # ? set $isput ? set -like? # add 'ortholog'
#   s/\s*[Oo]rtholo(gy|gue|g)\s+\d+//g; s/\b[Oo]rtholo(gy|gue|g)\s*//g; 
#   s/^(of|with)\s+//ig; # bad leading words
#   # Add protein to the end when ending in:  'binding|domain|like|related'
#   s/\b(binding|domain|like|related)\s(\W*)$/$1 protein $2/;  
#   if( s/[,\s]*putative//g ) { $isput=1; } #s/putative, putative/putative/;
# 
#   # punctuation
#   s/[\|]/:/g; s/#/n/g; 
#   # s/_/ /g; # or leave uscores ??
#   s/\s*$//; s/^\s$//; # lead/end spaces
#   s/\s*[\/\.,;_-]*$//; # trailing punc
#   s/[.] /, /g; # no sentences?
#   s/^\W+//; # no leading crap
#   # s/([,;:]){2,}/$1/g; #? double punc
# 
#   # plurals ?? anything ending in 's' ?
#     
#   # SEQ_FEAT.ProteinNameEndsInBracket:  Phosphoenolpyruvate carboxykinase [ATP] << change [] for ncbi
#   if(/\]$/) {  s/\[/\(/g; s/\]/\)/g;}
#   # unbalanced brackets; # add {} ? not used; <> ? not brackets
#   if(/[\(\)\[\]]/) {
#     my($nb,$ne,$d);
#     $nb= tr/\[/\[/; $ne= tr/\]/\]/; $d=$nb - $ne;
#     while($d>0) { s/\[//; $d--; } while($d<0) { s/\]//; $d++; }
#     $nb= tr/\(/\(/; $ne= tr/\)/\)/; $d=$nb - $ne;
#     while($d>0) { s/\(//; $d--; } while($d<0) { s/\)//; $d++; }
# #  ##..  not working
# #     my @bp=( "[","]", "(",")" ); 
# #     while ( my($b,$e)= splice( @bp, 0, 2) ) { 
# #     my $nb= tr/$b/$b/; my $ne= tr/$e/$e/; my $d= $nb - $ne; 
# #     while($d>0) { s/\Q$b/./; $d--; } while($d<0) { s/\Q$e/./; $d++; }
# #     } 
#   }   
#   
#   $_=$NAME_UNK unless(/\w\w/); # 'Conserved protein' becomes blank
#   if($isunk) { # regularize, but check/keep some additions
#     if(/^($NAME_NONE)$/i or /^($NAME_NONE) protein/i) { $_= $NAME_UNK; }
#     elsif(/^($NAME_NONE)\s*\w+/) { s/^($NAME_NONE)\s*/$NAME_UNK /;} # what?
#     else { 
#       # ?? save old,clean-name as Note in some/all cases.
#       if($pi >= $MIN_IDLIKE) { $islike=1; } ##  for > 10-15% ident? or any?
#       # Name=Mitogen-activated protein kinase kinase kinase 7-interacting protein 2
#       # >> Mitogen-activated kinase kinase kinase 7-interacting protein 2-like protein
#       else { $lowqualname=$_; $_= $NAME_UNK; }  # replace entirely? or not; NOT yet
#     }  
#   }
#   if($islike) { unless(/family/ or /\blike/ or /protein kinase/) { s/\s+protein//; $_ .= "-like protein"; }  }
#   s/protein protein/protein/ig; # other stutters?
#     
#   $_ .= ", putative"  if($isput and not $isunk);
#   $_ .= " $locusname" if($isunk and $locusname);
#   ## ncbi complains about this ^^ locusname addition; leave out at least w/ config.
#   return ($_, $lowqualname);
# }


=item processgene

  format structured model gene/(mRNA|ncRNA)/exon,CDS

  maybe add processfeature, for other things like transposons, ...
  processgene does this now if no mRNA.
  
  1719	638	repeat_region
    db_xref	baggins{}1471
    db_xref	FLYBASE:FBti0020395
    transposon	baggins{}1471

=cut



## revise isPartial .. use aaqual for stopcodon
# my($aasize,$partial,$nxxaa,$prothasstopc,$hasselcstop)= aaQuals($attr,$prot,)
sub aaQuals {
  my( $attr, $prot )= @_; # mrnaAttr

  my ($at, $aalen, $cdsp, $aaqual, $selcstop)= ("") x 9;
  unless($prot) {
    ($prot) = $attr =~ m/protein=([^\s;]+)/; # use for alt-product=same/isoform marking
  }
  ## attr: aalen=1543,54%,complete
  ## Name=Selenoprotein N;Selcstop=1412;aalen=556,80%,complete,selcstop; << look for selcstop also
  ## aalen= dups? aalen=123; aalen=145,80%,complete < check for both or expect cleaner input?
  ($at)= $attr =~ m/aalen=(\d+,[^\s;]+)/i; # aalen=\d+,pct,qual; from evigene
  if( $at ) { 
    ($aalen,$cdsp,$aaqual)= $at =~ m/(\d+),(\d+)%,(\w\S+)/;
    unless($aalen) { ($aalen,$cdsp,$aaqual)= $at =~ m/(\d+).(\d+)\W+(\w\S+)/; }
    ($selcstop)= ($at=~/selcstop/i)?1:0;
    }
  # aalen=\d+; from gmap2gff, other and aaSize=\d+   
  if( not $aalen  and ($at)= $attr =~ m/(?:aaSize|aalen)=([^\s;]+)/ ) { ($aalen)= $at =~ m/(\d+)/; }
  if( not $aaqual and ($at)= $attr =~ m/quality=([^\s;]+)/) { ($aaqual)= $at =~ m/Protein:([^\s,;])+/; }
    
  my($protc, $hasast, $nxx)=(0) x 9;
  if($prot) {  
    $protc |= 1 if($prot =~ /^M/);
    if($prot =~ /\*$/) { $protc |= 2; $hasast=1; } # FIXME dont assume * on complete prot. use quality complete.
    else {
      $protc |=2 if($aaqual =~ /complete/ or $aaqual =~ /partial5\b/); # even if have prot
    }
    $aalen= length($prot); # supercede size attr
    $aalen-- if($hasast); # use this standard, aalen excludes stop coden
    $nxx = $prot =~ tr/X/X/; # have 2 CDShasTooManyXs in other scaffs
  } else {  
    if($aaqual =~ /complete/) { $protc= 3; }
    else {
      $protc |=1 if($aaqual =~ /partial3\b/);  
      $protc |=2 if($aaqual =~ /partial5\b/);  
    }
  }
  
  my $partial  = ($protc == 3) ? "complete" : ($protc == 1)? "partial3" : ($protc == 2) ? "partial5" : "partial";
  # $partial.="53" if($partial eq "partial"); # option? only for isPartial() ?
  return($aalen,$partial,$nxx,$hasast,$selcstop);
}

sub isPartial {
  my( $mrna )= @_;
if(1) {  
  my($aasize,$partial,$nxxaa,$prothasstopc)= aaQuals($mrna->[8]);
  $partial.="53" if($partial eq "partial");
  return($partial);
} else {  
  my $attr= $mrna->[8];  
  my $partial= "";
  my ($prot) = $attr =~ m/protein=([^\s;]+)/; # use for alt-product=same/isoform marking
  if($prot) { # tinyaa not here.
    my $protc=0;
    $protc += 1 if($prot =~ /^M/);
    $protc += 2 if($prot =~ /\*$/); # FIXME dont assume * on complete prot. use quality complete.
    $partial = ($protc == 3) ? "complete" : ($protc == 1)? "partial3" : ($protc == 2) ? "partial5" : "partial";
  }
  unless($partial) {
    my ($qual) = $attr =~ m/quality=([^\s;]+)/;
    if($qual=~/Protein:\w*(complete|partial\w*)/) { $partial=$1; }
  }
  $partial.="53" if($partial eq "partial");
  return($partial);
}
}


# FIXME: add the new gene discrepancy checks here to bestgenes_update.pl
# .. see also genereplacex.pl where some of errors came from (i.e. replaced exons w/o adjusting mRNA span)

sub isSplit { 
  my($attr,$qual)= @_;
  my $issplit= 0; my $nsplit= 0; ## need also $nparts of Split=2/3 ..
  unless($qual){ ($qual)= $attr =~ m/quality=([^\s;]+)/; }
  $attr||=""; $qual||="";
  if($attr=~m/[;,]Split=(\w+)/) { $issplit=$1; ($nsplit)= $attr=~m/Split=$issplit.(\d+)/; }
  elsif($qual =~ /(split\w+)/i) { $issplit=$1; ($nsplit)= $attr=~m/split$issplit.(\d+)/i; }
  elsif( my($id)= $attr=~/\bID=([^\s;]+)/ ) { ## add ID=xxx_C1/2 check?; dont expect this here..
    if($id=~/_C(\d+)$/) { $issplit=$1; }
  }
  if($issplit and not $nsplit) { $nsplit=$issplit; $nsplit++ if($nsplit==1); } # look where else?
  return (wantarray) ? ($issplit,$nsplit) : $issplit;
}

sub force1strandgene {
  my( $geneall )= @_;
  my($gene,$mrna,$gc,$gb,$ge,$go, $gid, $ngerr, $nmrna)=(0) x 10;
  my $nga= @$geneall;
  
  my(%gor,%tor, $noor);
  for(my $j=0; $j < $nga ; $j+=2 ) {
    my $generec= $geneall->[$j];  # what is j+1 ??
    
    if(ref $generec and @$generec > 0) {
      my($tc,$tb,$te,$to,$attr,$id);
      
      ($mrna)= grep{ $_->[2] =~ /^($mrnatypes)$/ } @$generec; 
      next unless($mrna); #? any such bug?
      ($tc,$tb,$te,$to,$attr)= @{$mrna}[0,3,4,6,8];
      $nmrna++;
      my ($issplit,$nsplit)= isSplit($attr);
      my ($id)= $attr =~ m/ID=(\w+)/;
      $tor{$tc}{$id}=$to;
      $noor++ unless($to=~/^[+-]/);
      my %mxo; $mxo{$to}++ if($to=~/^[+-]/);
      my @exons= grep{ $_->[2] =~ /^exon/ } @$generec;  # not CDS
      # if($issplit) .. strand can be suspect ! add weight?
      if(@exons>1) { # only check these? skip 1x as maybe strandless
        map{ my $xo= $_->[6]; $mxo{$xo}++ if($xo=~/^[+-]/); } @exons;
      }
      my($mo)= sort{$mxo{$b}<=>$mxo{$a}} keys %mxo;
      unless($mo) { $noor++; $gor{$tc}{'.'}++; next; } # need gor{tc} place holder
      #?? $mo='+' unless($mo=~/^[+-]/); # any such bugs? yes some cases of '.' for split parts all on 1 chr

      my $cmo= $mxo{$mo}||0;
      $gor{$tc}{$mo} += $cmo; # or ++ 1 mrna?
      if($cmo < 2) { } # is suspect, mark changeable ..
      #? if($mo ne $to and $cmo>1) { }
      
      # ignore gene for now? or record or in tor{tc}{gid}? 
      if(my($gene)= grep{ $_->[2] =~ /^gene/ } @$generec) {
        ($tc,$tb,$te,$to,$attr)= @{$gene}[0,3,4,6,8];
        $noor++ unless($to=~/^[+-]/);
      }
      
    }
  }
  
  my($nfix,$nup,$nupexon,%gstrand,)= (0,0,0);
  $nfix+= $noor;
  my @tc= sort keys %gor;
  for my $tc (@tc) {
    my @ids= sort keys %{$tor{$tc}};
    my ($go,@mo)= sort{$gor{$tc}{$b}<=>$gor{$tc}{$a}} keys %{$gor{$tc}};
    unless($go and $go=~/^[+-]/) { # any such bugs? yes, splits all w/ .
     $go='+'; $gor{$tc}{$go}++;
    }
    my $cgo= $gor{$tc}{$go}||0;
    $gstrand{$tc}{'gene'}= $go; 
    if(@mo>0) {  } # what?
    for my $id (@ids) {
      if($tor{$tc}{$id} ne $go) {        
        $gstrand{$tc}{$id}= $go; $nfix++; # flag it
      }
    }
  }
  
  my %gid=(); my @fixid=();
  if($nfix) {
    for(my $j=0; $j < $nga ; $j+=2 ) {
      my $generec= $geneall->[$j];  # what is j+1 ??
      if(ref $generec and @$generec > 0) {
        my($tc,$tb,$te,$to,$attr,$id,$gid);
        
        my($mrna)= grep{ $_->[2] =~ /^($mrnatypes)$/ } @$generec; 
        next unless($mrna); #? any such bug?
        ($tc,$tb,$te,$to,$attr)= @{$mrna}[0,3,4,6,8];
        ($id)= $attr =~ m/ID=(\w+)/;
        ($gid=$id)=~s/_C\d+$//; $gid=~s/t\d+$//;

        my $go= $gstrand{$tc}{'gene'}; 
        $gid{$gid.':'.$go}++; # BUGGERS: - or + eaten by perl in some "" strings
        unless($go=~/^[+-]/) {
          warno "#WARN: force1strandgene BUG $id go=$go \n"; next;          
        }
        
        # my ($issplit,$nsplit)= isSplit($attr);
        #? need gstrand{$tc}{$id} or not?
        if($to ne $go) { 
          $mrna->[6]= $go; $nup++; push @fixid, $id.':'.$to;
          my @exons= grep{ $_->[2] =~ /^(exon|CDS)/ } @$generec;  # not CDS
          for my $x (@exons) {
            if($x->[6] ne $go) { $x->[6]= $go; $nupexon++; } # 
          }
        }
        if(my($gene)= grep{ $_->[2] =~ /^gene$/ } @$generec) {
          ($tc,$tb,$te,$to,$attr)= @{$gene}[0,3,4,6,8];
          if($to ne $go) { $gene->[6]= $go; $nup++; }
        }
        
      }
    }
    my @gid= sort keys %gid;
    warno "#WARN: force1strandgene nmrna=$nmrna, nfix=$nfix, noor=$noor, nup=$nup/$nupexon, gid=@gid, fixid=@fixid\n";
    # buggers: force1strandgene n=36824 is so wrong 
  }
  return($nfix,$nup,\@fixid);
}

  
sub processallgene
{
  my( @geneall )= @_;
  return unless(@geneall);
  my @res=(0) x 3;
  my($gene,$mrna,$gc,$gb,$ge,$go, $gid, $ngerr, $nmrna)=(0) x 10;
  
  force1strandgene(\@geneall); # if($public_options{force1strandgene});
  
  # check for parentgene and that it spans all mrna; pgene should be in geneall[0]..
  for(my $j=0; $j < @geneall ; $j+=2 )
  {
    my $generec= $geneall[$j];  
    if(ref $generec and @$generec > 0) { 
    
      # double dang: need to check all mRNA prots for partial; need protpart code below here
      #  set genepartial ONLY if partialspan == genespan
      # trap double gene records from outoforder mrna
      unless(ref $gene) { 
        ($gene)= grep{ $_->[2] =~ /^gene$/ } @$generec; 
        if(ref $gene) { # problem if we already set gb,gc from outoforder mRNA
          ($gc,$gb,$ge,$go)= @{$gene}[0,3,4,6] unless($gc and $ge); 
          $gid= $gene->[9];
          }
        }
        
      ($mrna)= grep{ $_->[2] =~ /^($mrnatypes)$/ } @$generec; 
      my($tc,$tb,$te,$to)= @{$mrna}[0,3,4,6];
      my $attr= $mrna->[8];  
      my ($issplit,$nsplit)= isSplit($attr);
      # my ($qual)= $attr =~ m/quality=([^\s;]+)/;
      # my $mrnaispartial= isPartial($mrna);
      $nmrna++;
      
      unless($gc and $ge) {
        my $dat= join(" ",@{$mrna}[0,3,4,6,8]); $dat=~s/;quality=.*$//;
        warno "#ERROR: mrna missing gene at ",$dat," $issplit\n"; 
        ($gc,$gb,$ge,$go)= ($tc,$tb,$te,$to); $ngerr++;
        
      } elsif( $go eq "." ) {
        $go= $to; $ngerr++; # $gene->[6]= $mrna->[6]; 
        warno "#WARN: fix missing strand: gene $gid:$go\n";  

      } elsif( ($gc and $tc ne $gc)
        or($go and $to ne $go)
        or($gb and $tb < $gb)
        or($ge and $te > $ge) # these mostly ???
        ) {  
          my $er=""; 
          $er.="to:$go/$to," if($go and $to ne $go); # problem left from swapmain in kfish2all5_fcds45.besttab5
          $er.="tb:$gb/$tb," if($gb and $tb < $gb); 
          $er.="te:$ge/$te," if($ge and $te > $ge);
          my $dat= join(" ",@{$mrna}[0,3,4,6,8]); $dat=~s/;quality=.*$//;
          warno "#ERROR: mrna ne gene =$er, at ",$dat," $issplit\n";   ## splitgene here, ok
          #NO# $gc= $tc; 
          #NO# $go=$to; #NO
          $ngerr++; $gb= $tb if($tb<$gb); $ge=$te if($te>$ge);
          }
      }
  }
  $ngerr++ unless(ref $gene);
  
  # force reorder gene into 1st geneall[0] or revise processgene
  if($ngerr) {
    my($src,$typ,$gv,$gph,$gattr,$gid);
    if(ref $gene) { ($src,$typ,$gv,$gph,$gattr,$gid)= @{$gene}[1,2,5,7,8,9]; }
    else { 
      ($src,$typ,$gv,$gph,$gattr,$gid)= @{$mrna}[1,2,5,7,8,9]; 
      $typ="gene"; $gv=1; $gid=~s/t\d+$//; $gattr="ID=$gid;gerr=$ngerr";
    }
    $gene= [$gc,$src,"gene",$gb,$ge,$gv,$go,$gph,$gattr,$gid]; 

    my $j=0; my $generec= $geneall[$j];  
    my @newrec= grep{ $_->[2] !~ /^gene$/ } @$generec; 
    push(@newrec, $gene);
    $generec= \@newrec; $geneall[$j]= $generec; 
    # .. something here is bad; dont get gene record output
  }
  

  my $geneispartial= undef;
  for(my $j=0; $j<@geneall; $j+=2 )
  {
    my $generec= $geneall[$j]; 
    my($mrna)= grep{ $_->[2] =~ /^($mrnatypes)$/ } @$generec; 
    my($tc,$tb,$te,$to)= @{$mrna}[0,3,4,6];
    my $mrnaispartial= isPartial($mrna);
    if( ($mrnaispartial or $geneispartial) and ($tb<=$gb or $te>=$ge)) {
      my $gp="";
      if($mrnaispartial =~ /partial/) { # turn on
        my $p5=( $mrnaispartial =~ m/5/ )?1:0;
        my $p3=( $mrnaispartial =~ m/3/ )?1:0;
        if($to eq "-") { if($p5 and $te>=$ge) { $gp.="5"; } if($p3 and $tb<=$gb) { $gp.="3"; }  }
        else { if($p5 and $tb<=$gb) { $gp.="5"; } if($p3 and $te>=$ge) { $gp.="3"; } }
        $geneispartial= ($gp =~ /[53]/) ? "partial$gp" : "complete"
        } 
      elsif($geneispartial =~ /partial/) { # turn off ?? or not
        $gp= $geneispartial;
        if($to eq "-") { if($te>=$ge) { $gp=~s/5//; } if($tb<=$gb) { $gp=~s/3//;  }  }
        else { if($tb<=$gb) {  $gp=~s/5//; } if($te>=$ge) { $gp=~s/3//;  } }
        $geneispartial= ($gp =~ /[53]/) ? $gp : "complete"
      } 
    } 
  }


  # while( my($generec, $geneother)= splice(@geneall,0,2) ) 
  for(my $j=0; $j<@geneall; $j+=2 )
  {
    my $generec= $geneall[$j]; 
    my $geneother= $geneall[$j+1];
    my @res1= processgene( $generec, $geneother, $gene, $geneispartial, $nmrna); 
    if($res1[1]>0) { $gene="done"; } # works, otherwise each mrna alt has preceding gene (dupl)
    for my $j (0..$#res1) { $res[$j]+=$res1[$j]; }
  }
  
  return @res;
}



sub processgene
{
  my($generecIN, $geneother, $geneIN, $genepartial, $nmrna)= @_;
  my($nref,$ngene,$ntr,$ncds,$nother)=(0) x 10;
  ##my $genepartial= undef; # FIXME;need all mRNA span tested
  
  # my @generec= @$generecIN; 
  my @generec   = sort _sortgene @$generecIN; # dont need?? better in case wrong.. NOTE: both strands sort fwd
  my @geneother = (ref($geneother)) ? @$geneother : ();
  
  if(not @generec and @geneother) {
    foreach my $x (@geneother) {
      my $ref = $x->[0];
      my $type= $x->[2];
      setrefout($ref) if($ref ne $lastref); # output handle for chrsplitting
      $nref+= putsource($ref) if($ref ne $lastref);
      putloc($type, $x);  $nother++;
      putattr($type, 0, $x->[8]); 
    }
    return($nref,$ngene,$ntr,$ncds,$nother);
  }
  
  my($mrna)= grep{ $_->[2] =~ /^($mrnatypes)$/ } @generec; # $mrnatypes
  my($gene)= grep{ $_->[2] =~ /^gene$/ } @generec; 
  my @exon = grep{ $_->[2] eq "exon" } @generec; # opt type?
  my @CDS  = grep{ $_->[2] eq "CDS" } @generec;  # opt type?
  
  #not working!??# 
  if($geneIN) { if(ref($geneIN)) { $gene=$geneIN; } else { $gene= undef; } }

  unless($mrna and @exon > 0) {
    my $id ="nada"; 
    if($mrna){ ($id)= $mrna->[8]=~/ID=(\w+)/; }
    elsif(@exon){ ($id)= $exon[0]->[8]=~/Parent=(\w+)/; } # KF: 2250 exon misses
    warno "#ERROR: Missing mRNA/exons: id=$id; mrna=$mrna; exons=@exon;\n";
    #no# putgene($generecIN, $geneother,"err=Missing-mrna-exon"); # flag="err=Missing-mrna-exon"
    return($nref,$ngene,$ntr,$ncds,$nother);
  }
 
  my $ref= $mrna->[0];
  setrefout($ref) if($ref ne $lastref); # output handle for chrsplitting
  $nref += putsource($ref) if($ref ne $lastref);
  
# .. special cases (quality=Poor,Partial,...)
#  Gene fragments: as gene with tags: note=nonfunctional due to frameshift ; note=frameshift ; /pseudo ; 
#   or as misc_feature with note=xxx
# CDS /exception="annotated by transcript or proteomic data" for quality=...Protein:curated_

  my $type= $mrna->[2];
  my $attr= $mrna->[8];
  my($id) = $attr =~ m/ID=([^\s;]+)/;
  my $gid= $id; #?? what of  gene=gid
  #?? ($gid)= $attr =~ m/;gene=([^\s;]+)/;

  ## FIXME split gene gid bad; here? 't1_C1' should be '_C1' 
  ## ID=Funhe2EKm009686t1_C1;gene=Funhe2EKm009686t1_C1;Split=1;isoform=1
  
  my ($qual)= $attr =~ m/quality=([^\s;]+)/;
  my ($issplit,$nsplit)= isSplit($attr, $qual);
  my $splitpartialtag= ($issplit)? ",Split=$issplit/$nsplit" : ""; # append to putCDSloc();
  
  ## IDSplitIsDone / IDSplitSuffix 
  ## PROBLEM: some uses here require IDSplitSuffix, some require no IDSplitSuffix ..
  ##  .. strip IDSplitSuffix then replace?
  
  my $evalt=0; #?? use this, w/ option instead of isoform=
  if($issplit) { if($gid =~ s/t(\d+)$IDSplitSuffix/$IDSplitSuffix/) { $evalt=$1;} 
    else { $gid =~ s/t(\d+)$//; $evalt=$1; } 
  } else { $gid =~ s/t(\d+)$//; $evalt=$1; }
  
  my $pid = reformatProteinID( $id); # (my $pid=$id)=~s/t(\d+)$/p$1/;
  
  # my ($aaqual)= $attr =~ m/aalen=\d+,\d+.,([^\s;]+)/;
  # if($qual =~ /Protein:([\s;]+)/) { $aaqual=$1; } #? both?
  # ^^now below: my($aasize,$aaqual,$nxxaa,$prothasstopc)= aaQuals($attr);
  
  ##  FIXME: Genbank tbl2asn needs new mRNA ID for these? same gene ID, not alt
  #KF: Check Split syntax .. look for Split=1|2 attr ??
  #Kf: Splits, diff attrs: 1. ";Split=1/2;", 2. ";mapCover=99%,Split=1/2;", 3. ";quality=...Map=xxx-[Ss]plit";
 
  ## quality=Class:Strong-expertchoice-splitgene:2,
  if(SPLITGENE_FIX) {
    # Add option that caller has already added splitgene ID tags, dont change here..
    ## isdone means 2 things: Split=1,2,3 valid part id AND ID=xxx is split-uniq id; 
    ## dont make sid new part ids when 1 true..
    # my $IDSplitIsDone = $config{general}->{'idsplit_is_done'} || 0; # input has already handled splitgene ids
    if($issplit) {
      #?? check for IDSplitSuffix ? IDSplitIsDone may be wrong
      # my $isdone=  $IDSplitIsDone; #? ignore this flag?? no, use 1st meaning
      # $isdone= ($id =~ m/$IDSplitSuffix/) ? 1 : 0; # 2nd meaning
      my $sid= $issplit; # this is sid if IDSplitIsDone
      unless($IDSplitIsDone) { $sid= 1 + ($splitgene{$id} || 0); }
      $splitgene{$id} = $sid; # increment
      $splitgene{$pid}= $sid; # equivalence prot
      $splitgene{$gid}= $sid;   ## gid also ??
    }
  }
  
  my $partof=0; #? rely on this partof: annot?
  if($qual =~ /\-partof:(\w+)/) { $partof= $1; } # reclass as misc_feature fragment

  my $cdsexcept=0;      ## add other key for TSA? trassembly ? validated?
  ## ^^vv add splitgene prots here as Protein:validsplit
  # if($issplit) { unless($qual =~ s/Protein:/Protein:validsplit/) { $qual=.",Protein:validsplit"; } }
  
  if($qual =~ /Protein:(curate\w*|valid\w*)/ or $issplit) {  my $ptag=$1;
    ## drop proteomic: $cdsexcept= "annotated by transcript or proteomic data";
    ## not with ncbi picky soft: [SEQ_FEAT.ExceptionProblem] annotated by transcript data is not a legal exception 
    $cdsexcept= "annotated by transcript or proteomic data";
    ## see below now
    # my $cid=""; ## move this after trasm work below...
    # if($attr =~ m/\btrid:([:\w-]+)/) { $cid= $1; } ## added mid-fix = mapper> trid:gspl2x11:ID
    # if(!$cid and my($rx)= $attr =~ m/rxnote=([^;\n]+)/){ ($cid)= $rx=~ m/cdnabest.(\w+)/; }
    # unless($cid) { ($cid)= $attr =~ m/cdnabest.(\w+)/; }
    # #? unless($cid =~ /[a-zA-Z]\d/) { $cid="evgasm"; } # FIXME
    # $cdsexcept .=";inference=similar to RNA sequence, mRNA:TSA:$cid" if($cid);
  }  
  # also CDS needs /inference=    my $v= "similar to RNA sequence, EST:$DBXEST:$pi"; #% need some est-db config
  # my $v= "similar to RNA sequence, mRNA:$cdnaid"; #% need some est-db config

  my $altnum=0;
  #KF isoform= is valid annot
  #KF: $add .= ";tblan=trid:$trasmid";  # annot for tbl2asn TSA valid trasm
  my($trasm)= $attr =~ m/tblan=trid:([^;,\s]+)/; ## now is gmapxx:id ?
  my($oid)= $attr =~ m/oid=([^\s;]+)/;  # oid=vel4ma11:cacao3vel4sc1Loc938t2
  my @oid = split",",$oid;
  my @pred= grep /AUG/, @oid; # messy hack to distinguish preds from trasm
  
  my @trmap=(); #KF? add osrc= for map/pred source?
  my($osrc)= $attr =~ m/osrc=([^;,\s]+)/;
  if($osrc) { 
    my($prs)= $osrc =~m/(\w*AUG\w*)/; 
    push @pred,$prs if($prs); #?? need ids
    my($trs)= $osrc =~m/((?:gspl|gmap)\w*)/; 
    my @trd= grep{ not /AUG/ } @oid; # add trasm ID here?
    if($trs and @trd) { map{ my $sd="$trs:$_"; push @trmap,$sd; }@trd; }
    ## add gmap/gspl osrc but need ID, from oids?
    ## not good, lots of /inference="similar to RNA sequence, mRNA:TSA:gmap3n5h" << no ID
  }
  
  if($attr =~ /isoform=(\w+)/) { 
    $altnum=$1;  # only if altnum > 1? # add oid == TSA inference:I100% ;
    ## FIXME: swapmain=xxxt1,xxxt2 changes ID altnum but not isoform= 
    if($evalt and $evalt ne $altnum and $attr =~ /swapmain=/) { 
      $attr  =~ s'isoform=$altnum'isoform=$evalt'; # add it, but also parse ID for t1,2,altnum otherwise all are isoform=1
      $altnum= $evalt; 
      ## perl s/// chokes on some attr Name values !***
      ## "Quantifier follows nothing in regex; marked by <-- HERE in m/isoform=ID=Funhe2EKm002512t1;
      ## Name=H(+ <-- HERE )/Cl(-) exchange transporter 3;...
      ## .. ;/ at ./evigene2genotbl_kfish2.pl line 1648, <$inh> line 167860.
      ## "s///o" = compile patt 1 time; no
      ##  "s'''", "m''" == no substitutes:  m'@pattern'; or s'@[]*(+/'bob';
    }        
  } elsif($nmrna>1) { # ? or $evalt>1 ??
    $altnum= $evalt || 1; 
    $attr .= ";isoform=$altnum"; # add it, but also parse ID for t1,2,altnum otherwise all are isoform=1
  }
    
  if($trasm) {
    $attr .= ";trasm=$trasm"; # special key for evidence inference parsing .. fixme 
  } elsif(@trmap) {
    my $tr=join",",@trmap; $attr .= ";trasm=$tr";  # no both trmap, pred ??
  } elsif(@pred) {
    my $pr=join",",@pred; $attr .= ";prediction=$pr"; #? was $oid
  } elsif($altnum and $oid) {  #  and not @pred
    $attr .= ";trasm=$oid";  # ?? need other tag for trasm evidence but NOT tblan= ; ovrna= ??
  }

  my $cdsexnote="";
  if($cdsexcept =~ m/annotated by transcript/) {
    my $cid=""; 
    if($trasm =~ m/(\w+)$/) { $cid=$1; } # messy now trasm=mapsrc:id ? shoud put mapsrc elsewhere!
    elsif($attr =~ m/\btrid:([:\w-]+)/) { $cid= $1; } ## added mid-fix = mapper> trid:gspl2x11:ID
    if(!$cid and my($rx)= $attr =~ m/rxnote=([^;\n]+)/){ ($cid)= $rx=~ m/cdnabest.(\w+)/; }
    unless($cid) { ($cid)= $attr =~ m/cdnabest.(\w+)/; }
    if($cid) {
      $cdsexcept .=";inference=similar to RNA sequence, mRNA:TSA:$cid"; # cancel cdsexcept unless have cdsexinf ..
    } else {
      $cdsexnote= "$cdsexcept, without TSA support"; $cdsexcept="";
    }
  }
  
  my ($altprotsame, $altfirstid,  $tinyaa)= (0) x 10; # $protpart,
  # 201 ERROR: [SEQ_FEAT.DuplicateFeat] 
  ## see isPartial($mrna)
  
  my ($prot) = $attr =~ m/protein=([^\s;]+)/; # use for alt-product=same/isoform marking
  $prot ||= "";
  #? here or below# $prot="" if( $changelist{$id} =~ /MisMatchAA/i );
  
  # FIXME issplit prot .. does each part have same prot, diff protid _C1,2 ? 
  # ?? now all protein= on Split=1/ but that may not be 1st input order
  # only 1 w/ prot? no, need protid for each w/ CDS, and that implies also prot value, needs to be same all parts
  # need cdsexception for all split parts to allow same prot > part-cds ?
  
  ## FIXME, allow prots from input.aa instead/replace these GFF prot attr
  ## $PEPfromGFF = 1 if(($proteins||$evidence{'proteins'}) eq "gff"); #KF was w/o -proteins flag
  if($prot{$id} and ($issplit or not $PEPfromGFF)) { $prot= $prot{$id} || $prot; }
  else { $prot{$id}= $prot if($prot); }
  
  # retest prot complete here, qual is out of date for rx=3 updates
  my($aasize,$protpart,$nxxaa,$prothasstopc)= aaQuals($attr,$prot);
  my $aasizeWithStop= $aasize; $aasizeWithStop++ if($protpart =~ /complete|partial5/);
  #o my($aasize)= $attr =~ m/aaSize=(\d+)/; # trap tiny aa<40, recode misc_feature all? 30 total; 20 exceptions  BadProteinStart

  if($prot) { # tinyaa not here.
    ## replace w/ common aaQuals($mrnaattr)
    #o my $protc=0;
    #o my $aastop= ($prot =~ /\*$/ or $aaqual=~/complete/)?1:0;
    #o $protc += 1 if($prot =~ /^M/);
    #o $protc += 2 if($aastop); ## FIXME, use aaQual annot values for complete/partial ..
    #o $protpart = ($protc == 3) ? "complete" : ($protc == 1)? "partial3" : ($protc == 2) ? "partial5" : "partial";
    #o $aasize= length($prot); 
    #o my $nxxaa = $prot =~ tr/X/X/; # have 2 CDShasTooManyXs in other scaffs

    $tinyaa= "$nxxaa XXX" if($nxxaa * 2 >= $aasize);
    if($prot{$gid}) {
      $altfirstid= $gid."t1"; # fixme ; bad for _C1/2 split genes
      $altprotsame= ($altnum>1 and $prot{$gid} eq $prot)?1:0;
      # DO need to check all alts..  ** FIXME: out of order alts, t2 > t2 w/ same prot
      my $maxalt= _max($altnum,15);
      for ( my $t=1; $t <= $maxalt and not $altprotsame; $t++ ) {
        next if($t eq $altnum); my $aid=$gid."t$t"; 
        if($aid ne $id and $prot{$aid} and $prot{$aid} eq $prot) { $altprotsame=$t; $altfirstid=$aid; }
      }

    # From murphyte@ncbi.nlm.nih.gov  Wed Feb 15 08:19:25 2012
    # SEQ_FEAT.CDSmRNAmismatch                1536 -- from a quick glance, it looks like you have multiple mRNAs paired 
    # with the same CDS feature (probably UTR variants). For these, you need to produce duplicate CDS features (with a s
    # eparate protein_id) so there is 1:1 correspondence
  
      if($altprotsame) { 
        #NO# $pid= reformatProteinID( $altfirstid); # $pid =~ s/p\d+$/p1/; 
        # ** tbl2asn has some hidden controls of ERROR:overlaps another CDS with the same product name
        # ..~/Desktop/dspp-work/genomesoft/ncbi201107/api/sqnutil3.c : HasOverlapComment: overlap frameshift .. key words
        $attr .= ";Note=alternateUTRof:$altfirstid"; # CDS overlap for 
        }
    } else { 
      $prot{$gid}= $prot; 
    }
  }
  if($aasize>0 and $aasize < $TINYAA and not $prot) { $tinyaa=$aasize; }

  # edit Name= here to change NO_NAME to NO_NAME id ??
  ## FIXME: Name/product for Unknown/Unchar should change to ID! And/or do for all Names that are not uniq?
  ## if isoform, add to product ?? SHOULD if alt-prots differ, may in any case to solve tbl2asn bugs
  #   product     CCC_04562 isoform B
	#		note alternatively spliced

  ## FIXME: need to rename on GENE basis, not per transcript... all tr should have/not have gid,isoform added.
  ## FIXME2: need lots of NCBI-name-police fixes; see discrep report
  my($name)= $attr =~ m/Name=([^;]+)/; 
  unless($name) { $name= $NAME_UNK;  $attr =~ s/$/;Name=$name/; }
  
  # badname need to set newna per below nameclean(), keep old name as is
  my ($newna,$lowqualname,$nadiff); 
  $newna= $name;
  if( $changelist{$id} =~ /^badname/i or $changelist{$gid} =~ /^badname/i ) { 
    my $nna="";
    if( $changelist{$id} =~ /^badname\s+(.+)/i) { $nna=$1; }
    elsif( $changelist{$gid} =~ /^badname\s+(.+)/i) { $nna=$1; }
    else { ($nna)= $attr =~ m/oname=([^;]+)/;  $nna="" if($nna=~/^same/); }
    $newna= $nna if($nna);  
  }
  
  my $pi= ( $newna =~ m/\((\d+)%.\)/ ) ? $1 : $MIN_PROTIDENT; 
  my $locusadd= ($NAME_UNKADDLOCUS) ? "locus $gid" : "";
  $newna =~ s/\s*\((\d+)%.\)//; $name =~ s/\s*\((\d+)%.\)//;  # (73%P); pull of this pct ident crap before clean.
  
  ($newna,$lowqualname,$nadiff)= nameclean( $newna, $pi, $locusadd );
  # $nadiff ==1 == NAMEDIFF_MINOR; ==2 == NAMEDIFF_MAJOR
  
  #NOTE rename Thecc1EG026898t1: 'Dynamin-2B' < 'Dynamin-2B (78%U)'
  #NOTE rename Thecc1EG026898t2: 'Dynamin-2A, putative isoform 2' < 'Dynamin-2A, putative (82%U)'
  #NOTE rename Thecc1EG026898t4: 'Dynamin-2A, putative isoform 4' < 'Dynamin-2A, putative (85%U)'
  #NOTE rename Thecc1EG026898t3: 'Dynamin-2A, putative' < 'Dynamin-2A, putative (82%U)'
  ## ^^ problem here re gene-level isoform naming ... but see 2B, 2A also.

# http://www.ncbi.nlm.nih.gov/genbank/eukaryotic_genome_submission_annotation/#Alternativelysplicedgenes
#        Alternatively spliced genes
#   In many cases a gene can be alternatively spliced, yielding
#   alternative transcripts. These transcripts may differ in the
#   coding region and produce different products, or they may differ
#   in the non-translated 5' or 3' UTR and produce the same protein.
#   To annotate alternatively spliced genes, include one mRNA and CDS
#   for each transcript, and include only one gene over all of the
#   features. Give the corresponding mRNA and CDS the same name, and
#   include a note "alternatively spliced" on each. If there are
#   multiple CDS with the same name, then add a note to each mRNA and
#   CDS to refer to each other, eg "transcript variant A" and
#   "encoded by transcript variant A" for one mRNA/CDS pair. If the
#   CDS have different translations, then they should have different
#   product names. Make sure that all the proteins have unique
#   protein_id's.
#
# test reformatting for ERROR: [SEQ_FEAT.DuplicateFeat] << match NCBI euk submit doc but no joy, all variants give same err
#  ... tbn2asn seems deficient in handling this:
#  mRNA	 note  transcript variant 3; alternatively spliced
#        product	D3-type cyclin, putative
#				 locus_tag	TCM_033764
#  CDS	 note  encoded by transcript variant 3; alternatively spliced;
#        product	D3-type cyclin, putative
       
       #Yes? and not $altprotsame
  # if($altnum > 0  and not $altprotsame) { $newna .= " isoform $altnum";  } # should do only when prot differs from maintr
  if($altnum > 0) { $newna .= ($altprotsame) ? " isoform $altprotsame" : " isoform $altnum";  }  
    # ^^ problem for altprotsame isoform 1 << main t1 (usually) lacks isoform 1 suffix. ie it isnt alternate
  # if($altnum > 0) { $newna .= ($altprotsame) ? "" : ($altnum > 1)? " isoform $altnum" : "";  }  
  # if($altnum > 0) { $newna .= " isoform $altnum";   }  #$attr .= ";Note=alternatively spliced"; 
    # ERROR: valid [SEQ_FEAT.DuplicateFeat] Features have identical intervals, but labels differ
  #above# $attr .= ";Note=alternatively spliced" if($altnum>0);
  
  if($name ne $newna) { ## $nadiff != 0 or nadiff == NAMEDIFF_MAJOR
    $attr =~ s/Name=\Q$name/Name=$newna/ ;  #?? add back pctident ?
    $attr =~ s/$/;lowqualname=$lowqualname/ if($lowqualname); # need other tag to put in CDS
    # fixme: sometimes, always add oldname as Note, esp if recoded to NONAME from function
    my $na= $name; $na =~ s/\s*\(\d+%.*\)//; $na =~ s/\s*sym:\S+//; # all/most have
    warno "#NOTE: rename $id: '$newna' < '$name'\n" if ($DEBUG and index($newna, $na)<0); # ignore also isoform..
  }
  
  $attr .= ";protid=$pid" if(@CDS);  # display in both mRNA and CDS records
    # NOTE: above may have changed pid to main pid if altprotsame
    #Wrong?? > ALSO: skip @CDS output if $altprotsame (but enter protid in mRNA)
    # terence murphy, wasp submit: need CDS for all mRNA (but probably use shared protid)
  
  # do here now  
  my $partial= $protpart; # as long as we have full updated prot here, that is enough
  my ($cdserr,$aabigger)=(0,0);
  ($cdserr,$aabigger)= checkCDSexons($id, $mrna, $aasizeWithStop, \@exon, \@CDS); # do above; add curated aa checks, partial annot
  if($aabigger > 0 and $cdsexcept and not $issplit) {  
    # FIXME This partial isn't always true, end points right but MGap internal genome-assembly gaps/mistakes
    #old# $partial="partial"; #?? always, or cdsexcept : BAD, drop
    #ok# warno "#WARN: $id prot=$protpart NOT partial but cdsexcept\n" if ($protpart=~/complete/ and $DEBUG);  
    $cdserr=0;    
  } elsif($aabigger > 0 and $cdserr and $issplit) {
    $cdserr=0 if($cdserr =~ /cdsoffby:/); # cancel annoying not-error from split cases
    ## spl2 exon ne mRNA: cdsoffby:1125,
  }
  
  if($qual=~/Protein:\w*(partial\w*)/) { # can also be wrong other way: NOT complete
    my $qpartial=$1; # always ignore qual value ?
    unless($protpart) { $partial= $qpartial; }
    elsif($protpart and $qpartial ne $protpart) { ## FIXME
      $attr =~ s/Protein:\w*(partial\w*)/Protein:$protpart/; 
      warno "#WARN: reclass $id Protein:$qpartial > $partial\n" if ($DEBUG);  
    } 
  }
  elsif($qual=~/Protein:\w*(complete\w*)/) {  
    my $qpartial=$1; # always ignore qual value ?
    if($protpart and $qpartial ne $protpart) { ## FIXME
      warno "#WARN: reclass $id Protein:$qpartial > $partial\n" if ($DEBUG);  
    } 
  }
  # ERROR: valid [SEQ_FEAT.PartialProblem] == out-of-date qual annot; see 'rx=' attr
 
  ## FIXME V2: add MGap errors from gsplign to this misc_feature, locus map error listing
  ## fixme V2 update: exonMapErr: locus_tag wanted, diff from gene gid
  use constant XMEend => 1;
  my($nxmaperr,@xmaperr)= exonMapErr(\@exon,$gid,$mrna);
  if(XMEend and $nxmaperr) { push @geneother, @xmaperr; } # problems? only inside put mRNA/@exon set?
  
  $attr .= ";".$genescore{$id} if($genescore{$id}); 
  
  # FIXME: need gene locus_id also (addgene); for partof, change gene=; for tinyaa, add new gene.
  # .. cancel this for Protein:curated == $cdsexcept
  if(($partof or $tinyaa) and not $cdsexcept) { # reclass as misc_feature fragment
    $type= "misc_feature";
    # FIXME: partof:notgoodgene, remove geneid, keep Note=partof
    if($partof) {
    my $partofid= reformat_myid($partof,$splitpartialtag); #? or as MySRC:id
    my $partgene= ";gene=".$partof; $partgene=~s/t\d+$//;
    if( $changelist{$id} =~ /notpartof/i ) { $partgene=""; }
    $attr =~ s/^/Note=Gene fragment of $partof;Dbxref=$partofid;/;
    $attr =~ s/;gene=\w+/$partgene/;
    $gene= undef;
    } elsif($tinyaa) {
    ## need to putgene 1st if doing so. BUT gene needs mRNA or other such parts; needs CDS unless /pseudo or nc
    ## "misc_feature without a corresponding gene feature" ..
    ## FIXME: cant put CDS with huge gap: CDShasTooManyXs in other scaffs; change to /pseudo or chop off gap ?
    ##   $tinyaa= "$nxxaa XXX" if($nxxaa * 2 >= $aasize);
    ## only 2 cases in cacaogd:

    $attr =~ s/^/Note=Gene fragment aasize:$tinyaa;/;
    $attr =~ s/;(locustag|gene)=\w+//g; # cant have locus_tag here
    $gene= undef;
    }
    my($attradd)= putCDSloc($type, $partial.$splitpartialtag, @exon);  $nother++;
    putattr($type, $id, $attr); 
    @CDS=(); @exon= (); 
  }

  if($ADDGENE and $gene) { # FIXME. Defer this to after last alt-tr (and all parts of gene)
    # .. fix cases where later alt-tr is partial including gene span (and check all inside gene span)
    # .. maybe change filter_gff to collect all altnum into generecIN, but then
    # .. need 2-level array: mRNA1/exons/cdss need single package.
    $genepartial= $partial unless(defined $genepartial);
    my $gtype= $gene->[2];
    my $gattr= $gene->[8];
    if($gene->[6] eq ".") { my $go=$mrna->[6]; $gene->[6]= $go; warno "#WARN: fix missing strand: gene $gid:$go\n"; }  
    putCDSloc($gtype, $genepartial.$splitpartialtag, $gene);  $ngene++;
    putattr($gtype, $gid, $gattr);
  }
  
  # * Split gene needs Note attribute .. add to config Split=1/2 .. needs ncbi-format id of splits 
  # [m]RNA, but may have retyped to @geneother...
  if(@exon) {
  putCDSloc($type, $partial.$splitpartialtag, @exon);  $ntr++;
  putattr($type, $id, $attr); 

  if(not XMEend and $nxmaperr) {
  foreach my $x (@xmaperr) {   
    # does this need to follow @exon before CDS for locusid match? need test subset
    my $type= $x->[2]; putloc($type, $x);  putattr($type, 0, $x->[8]); $nother++;
  } }

  }
  
 # terence murphy, wasp submit: need CDS for all mRNA (but probably use shared protid)
# From murphyte@ncbi.nlm.nih.gov  Wed Feb 15 08:19:25 2012
# SEQ_FEAT.CDSmRNAmismatch                1536 -- from a quick glance, it looks like you have multiple mRNAs paired 
# with the same CDS feature (probably UTR variants). For these, you need to produce duplicate CDS features (with a s
# eparate protein_id) so there is 1:1 correspondence
#.. BUT this leads to new error with splitgenes .. before no CDS on 2nd+ part, but should have it, now get this
#  64 REJECT: [SEQ_FEAT.MultipleCDSproducts]
#  64 REJECT: [SEQ_INST.IdOnMultipleBioseqs]


  if(@CDS) {  #  and not $altprotsame
    $type= "CDS";
    my $cdsattr= $attr; # $mrna->[8];
    # already in attr now: $cdsattr = "protid=$pid;$cdsattr";
    $cdsattr .= ";partial=$partial" if ($partial =~ /partial/); 
    $cdsattr .= ";exception=$cdsexcept" if($cdsexcept);
    $cdsattr .= ";Note=$cdsexnote" if($cdsexnote); # either or, not both

        # ** tbl2asn has some hidden controls of ERROR:overlaps another CDS with the same product name
        # ..~/Desktop/dspp-work/genomesoft/ncbi201107/api/sqnutil3.c : HasOverlapComment: overlap frameshift .. key words
        # RemoveCodingRegionsWithSuppressionWords() =     if (DoesStringContainPhrase (product, "ABC", TRUE, TRUE) .. transposon ..
        # OverlappingProductNameSimilar()
        # update .conf : attr_CDS cdsNote=
    # $cdsattr .= ";cdsNote=CDS overlap for alternateUTRof:$altfirstid" if($altprotsame); # << in mrna attr but NOT here in cdsattr
    #see isoform# 
    #above# $cdsattr .= ";cdsNote=alternatively spliced" if($altnum>0);
    
    ## handle some attr for this for valid CDS <> protein
    ## exception: 	"unclassified transcription discrepancy" 
    ## need regex for Dbxref > mRNA, eg. ncbi [XN]M_ ids vs [XN]P_ ids

    #above# my ($cdserr)= checkCDSexons();  
    my $putpep= $PutPEPfromFile; # use config, 0 default
    $putpep=1 if($cdsexcept);
    if($putpep) { 
      #was PutPEPfromFile was $PEPfromGFF  ##  and not $altprotsame
      ## PutPEPfromFile or PEPfromGFF isnt right, dont put all, ie only bestaa=pubaa, need new opt?
      ## do if cdsexcept ..
      ## new changelist MisMatchAA = skip my pep, let ncbi translate. see putCDS
      $prot="" if( $changelist{$id} =~ /MisMatchAA/i );
      
      #o# my $pepok = putprot( $pid, $prot);
      my $pepok = putprot( reformat_myid( $pid,$splitpartialtag), $prot);
      $cdserr .= "noaa" unless($pepok);
    }
    
#     if($partial eq "partial5" or $partial eq "partial") { # do this in putCDS
#       $cdsattr .=";codon_start=".$phase+1;
#     }
    my($attradd)= putCDSloc($type, $partial.$splitpartialtag, @CDS);  $ncds++;
    putattr($type, $id, $cdsattr.$attradd);
    
    #KF try8: 4192/4094 splits 1/2 exon ne; 18933 0split ne; are these the pubaa set? need cdsexception check
    warno "#ERROR: $id CDSexon, spl$issplit exon ne mRNA: $cdserr\n" if($cdserr); # can have cdsold,new lines appended
  }
 
  foreach my $x (@geneother) {
    my $type= $x->[2];
    putloc($type, $x);  $nother++;
    putattr($type, 0, $x->[8]); 
  }
  
  return($nref,$ngene,$ntr,$ncds,$nother);
} 

=item exonMapErr

  FIXME V2: add MGap errors from gsplign to this misc_feature, locus map error listing
  FIXME: stutter, have same error from multiple alt exons .. 
  # add here? as Region/misc_feature following mRNA/CDS ..
   exon annots: 13651 error from gmap; keep/report in gb.annot? Region=map-exception...
   Parent=Funhe2EKm000208t1;error=ERROR.span:genome_span:1070,tr_span:1241,KN805525.1:5989263-5988193171
   Parent=Funhe2EKm000500t1;error=ERROR.span:genome_span:3687,tr_span:3875,KN805527.1:3629399-3625712188;
   ERROR: Bad location on feature misc_feature (start 3629398, stop -669255109)

    # possible misc attr:  /inference="[CATEGORY:]TYPE[ (same species)][:EVIDENCE_BASIS]"  
    # gene == /locus_tag="text" 

  from gmap: Funhe2EKm009624t3
     misc_feature    19132..19282
                     /locus_tag="D326_009624"
                     /note="trmap_genome_error:genome_span:150,tr_span:363"
  from gsplign: Funhe2EKm009624t2 MGap:1145-1351 b/n exon3:469-1144,exon4:1352-1429
     misc_feature    20055..21252
                     /locus_tag="D326_009624"
                     /note="trmap_genome_error:tr_align_gap:1145-1351"
    
  fixme V2 update: exonMapErr: locus_tag wanted, diff from gene gid, get from mrna
  
  ?? cancel dup maperr locations?
  WARNING: valid [SEQ_FEAT.DuplicateFeat] Features have identical intervals, but labels differ FEATURE: 
    misc_feature: map_error:genome_span:413,tr_span:768 [lcl|KN811437.1:22476-22889] [lcl|KN811437.1: delta, dna len= 1503504]
  
    same genome locus for 3 alts of locustag: D326_033565, but diff tr spans
  grep  '^22476  22889'  kfish2rae5h_fc14_knset8.tbl
  22476	22889	misc_feature <<?? gff says 22478<< not 76 ..	22889
        note	map_error:genome_span:413,tr_span:579 for Funhe2EKm033565t1
  22476	22889   # mrna exon w/o map err?
  22476	22889	misc_feature
        note	map_error:genome_span:413,tr_span:768 for Funhe2EKm033565t2
  
=cut


my(%didxerr);  

sub exonMapErr {
  my($exons,$geneid,$mrna)= @_;
  ## add gsplign MGap handling: now in mRNA annot as Target spans, match to exons/between-exons
  ## use same tag? eg. miscfeat span = exon1e .. exon2b, Note=err:trmapgap:[MGap 200-300]

  ## ADD here gapfill= as misc_feature, on exons and tagged in mRNA attr
  ## 2015.04.20: insert gapfix.gff, replacing others it is drawn from .. annots are output of this script
  ## gmapnsub3/kfish2rae5h_fc14m_gapfix.gff.gz mRNA n=3121, 
  ##   mRNA should have either gapfix=xxx n=133 or gapfill=xxx n=2043 annots to replace others, some lack, skip?
  ##   eg: ID=Funhe2EKm007709t5 insrc=kf2a:gmap2a5u;osrc=gmap2a5u;tblan=trid:gmap2a5u:Funhe2Exx11m005839t5;ggap=nnn/xxx,yyy;gapfill=115589-115592
  ##  .. exon.gapfill= needs to become misc_feature
 
  ## FIXME: some hit gaps, trigger ERROR:   SEQ_FEAT.FeatureBeginsOrEndsInGap ..
  ##  xgapfill handles as gap feature, others need some fix
  
  ## 
  my @xloc=(); # output misc_feature annots
   # as: my $xl= [$gr,$src,"misc_feature",$gb,$ge,1,$go,0,$xat,$tid];
  my @xerr= grep{ $_->[8]=~ m/error=ERROR.span/ } @$exons;
  my @xgapfill= grep{ $_->[8]=~ m/gapfill=/ } @$exons;
  my $mattr= $mrna->[8];
  
use constant CHECK_SPLIGN_EXON_ERRS => 1; # TEST
if(CHECK_SPLIGN_EXON_ERRS) {  
  if(not @xerr and ref $mrna) { # gsplign MGap ?
    my($gp)= $mattr =~ /gaps=([^;\s]+)/;
    ## Also see exon annot internal mgap trg=Funhe2EKm000098t1 1 1081;mgap=633-754;
    ## matches mRNA gaps=459,MGap:633-754,...
    if($gp=~/MGap/) { 
      my @mg= $gp=~m/MGap:([^,;]+)/g; # gaps=580,MGap:1371-1487,MGap:1000-1148,MGap:398-547,MGap:4545-4708,;
      my %mg=(); map{ my($b,$e)= m/(\d+).(\d+)/; $mg{$b}=$e; } @mg;
      my($ltb,$lte,$ltrb,$ltre)=(0,0,0,0);
      for my $x (@$exons) {
        my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$tid)= @$x;
        my($ingap)= ($tattr=~m/;ggap=([^;\s]+)/)?$1:0; # cancel if ingap overlaps gb,ge
        if( my($trb,$tre)= $tattr=~m/(?:Target|trg)=\S+.(\d+).(\d+)/ ) {
          for my $mb (sort keys %mg) { 
            my $me=$mg{$mb}; my $twid=1+$me-$mb; my($gb,$ge)=(0,0);
            if($mb<$tre and $me>$trb) { # inside exon
              if($to eq '-') { }
              else { $gb=$tb+($mb-$trb); $ge=$tb+($me-$trb); } ## REV STRAND needs te - xxx
            } elsif($ltb and $me<$tre and $mb>$ltb) { 
              if($to eq '-') { }
              else { $gb=$ltb+($mb-$ltb); $ge=$ltb+($me-$ltb); }
            }
            if($gb>0 and $ge>$gb) {
              my $xer="genome_gap,tr_span:$twid,$ref:$gb-$ge";
              $tattr =~s/$/;error=ERROR.span:$xer/; $x->[8]=$tattr; push @xerr,$x; last; 
            }
          }
        ($ltb,$lte,$ltrb,$ltre)= ($tb,$te,$trb,$tre);
        }
      }
    }
    # hassle need to match MGap tr-spans to exon Target/trg spans;
    # should revise gsplign2gff to keep MGap on exons
    # mimic: error=ERROR.span:genome_span:2315,tr_span:2428,KN811437.1:691859-689544,
  }
}
  
  return(0) unless(@xerr or @xgapfill);
  my $locustag= (ref $mrna and $mattr =~ /locustag=([^;\s]+)/) ? $1 : 0;
  
  for my $x (@xgapfill) { ## 201504 add
    my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$tid)= @$x;
    my($xerr)= $tattr =~ m/gapfill=([^;\s]+)/; # gapfill=115589-115592 == genome spans filled by trasm
    if($xerr) {
      my($gb,$ge)= $xerr =~ m/(\d+)-(\d+)/;  
      my $gr=$ref;
      my $go=$to; # which strand?  #my $go= ($gb>$ge)?'-':'+'; # which to believe?
      if($to eq '.') { $go= ($gb>$ge)?'-':'+'; }
      ($gb,$ge)=($ge,$gb) if($gb>$ge);
      next if($ge < $gb+2); # tiny gap fills .. filter before this .. how tiny? 1 is bad for ncbi submit, 2?
      
      my $xerrid= join",",$xerr,$geneid;  
      next if($didxerr{$xerrid}++); # SEQ_FEAT.DuplicateFeat.Features have identical intervals, but labels differ
     
      my $tid4ncb= reformat_myid($tid);       #  gnl|Funhe2GD|$tid
      my $xat = "Note=transcript $tid4ncb fills genome gap;"; 
            
      ## ** NO GOOD, NCbi disallows any gap-spanning feature .. need to inset as 'gap' with linkage_evidence=align trnscpt
      #o if($locustag) { $xat .= "locustag=$locustag;"; } #? or Parent=pid ??
      #o elsif($geneid) { $xat .= "gene=$geneid;"; }#? or Parent=pid ??
      #o my $xl= [$gr,$src,"misc_feature",$gb,$ge,1,$go,0,$xat,$tid];
      #o push @xloc,$xl;
      
      ## two choices: 'assembly_gap' or 'gap': which? try both?
      ## gap: estimated_length=$gapsize;inference=xxx:mRNA:TSA:ID;note=xxx;
      ## assembly_gap: estimated_length=$gapsize;gap_type=within scaffold;linkage_evidence=align_trnscpt;
      
      my $gapsize=1+$ge-$gb;
      my $gtype="assembly_gap";
      $xat .= "estimated_length=$gapsize;";
      
      ## use agap if agap == gapsize ? .. where is that info?
if(0) {      
      ## this one fails at least when this gap abuts rest of asm gap not in trasm
      $gtype="assembly_gap";
      $xat .= "gap_type=within scaffold;linkage_evidence=align_trnscpt;";
      #not allowed# if($locustag) { $xat .= "locustag=$locustag;"; } #? or Parent=pid ??
      #not allowed# elsif($geneid) { $xat .= "gene=$geneid;"; }#? or Parent=pid ??
} else {
      ## ** Note: ncbi tbl2asn collapses asmgap + gap when same span, leaves out this form's inference, keeps note
      $gtype="gap";
      my $cid=""; 
      if($mattr =~ m/trasm:([:\w-]+)/) { $cid=$1; } # messy now trasm=mapsrc:id ? shoud put mapsrc elsewhere!
      elsif($mattr =~ m/\btrid:([:\w-]+)/) { $cid= $1; } ## added mid-fix = mapper> trid:gspl2x11:ID      
      if($cid) {
        #o $xat .="inference=similar to RNA sequence, mRNA:TSA:$cid;" if($cid); #?? this or align?
        my ($ec)= evidencecode('inference', 'trasm', $cid);
        $xat .= "inference=$ec;";
      }
      #not allowed# if($locustag) { $xat .= "locustag=$locustag;"; } #? or Parent=pid ??
}

      my $xl= [$gr,$src,$gtype,$gb,$ge,1,$go,0,$xat,$tid];
      push @xloc,$xl;
     
      }
  }
  
  for my $x (@xerr) {
    my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$tid)= @$x;
    my($xerr)= $tattr =~ m/error=ERROR.span:([^;\s]+)/;
    if($xerr) {
      my($gs,$ts,$gloc)=split",",$xerr;
      my($gr,$gb,$ge)= $gloc =~ m/(\S+):(\d+)-(\d+)/; #** ERROR.span: gb,ge can be rev, ignore? and use to
      ## WHAT THE FK it is giving rev strand from exon to misfeat: gb>ge cases
      my $go=$to; # which strand?  #my $go= ($gb>$ge)?'-':'+'; # which to believe?
      if($to eq '.') { $go= ($gb>$ge)?'-':'+'; }
      ($gb,$ge)=($ge,$gb) if($gb>$ge);
      
      if( my($ingap)= $tattr=~m/;ggap=([^;\s]+)/ ) { # cancel if ingap overlaps gb,ge
        my($ob,$oe)= $ingap =~ m/ovspan:(\d+).(\d+)/ ;
        next if($ob <= $ge and $oe >= $gb); # any overlap to gap cancels miscfeat for NCBI errors.
      }

      my $gsp= $ge - $gb;  my $xsp= $te - $tb;
      if($ge>0 and $gr eq $ref and $gsp < 1.1*$xsp) {
        my $xat=""; #"Note=trmap_genome_error:$gs,$ts"; 
        if($locustag) {
          $xat .= "locustag=$locustag;"; #? or Parent=pid ??
        } elsif($geneid) {
          $xat .= "gene=$geneid;"; #? or Parent=pid ??
        }
        $xat .= "Note=map_error:$gs,$ts;"; # was trmap_genome_error:
        #o# my $xerrid=join",",$gloc,$geneid,$tb,$te; #? maybe skip tb,te part?
        my $xerrid= join",",$gloc,$geneid;  
        next if($didxerr{$xerrid}++); # dont stutter.. SEQ_FEAT.DuplicateFeat.Features have identical intervals, but labels differ
        my $xl= [$gr,$src,"misc_feature",$gb,$ge,1,$go,0,$xat,$tid];
        push @xloc,$xl;
        }
    }
  }
  my $nx=@xloc;
  return ($nx,@xloc);
}

# [SEQ_FEAT.CDSwithNoMRNAOverlap] 11 out of 2729 CDSs overlapped by 0 mRNAs : check here.
# new err: TransLen: CDSlen off by 1,2 longer than aa*3; need to chop CDS end
sub checkCDSexons
{
  my($id, $mrna, $aasize, $exon, $CDS)= @_; # should be loc-sorted 
  my $err=""; my $aabigger=0;
  return($err,$aabigger) unless(ref($CDS) or ref($exon));
  my($cr,$cb,$ce,$co)= @{$mrna}[0,3,4,6];
  my($xr,$xb,$xe,$xo,$cdsl)=(0) x 5;
  my $ci=0;
  # add test for len(prot) > len(cds/3); mark curated aa > cds as partial cds/mrna/gene
  foreach my $cx ($exon,$CDS) {
    my $ct=($ci==0)?"x":"c"; $ci++;
    foreach my $c (@$cx) { # better would be check all @CDS inside @exon
      my($tr,$tb,$te,$to)= @{$c}[0,3,4,6];
      $err .= "r$ct:$tr," if($tr ne $cr);
      $err .= "o$ct:$to," if($to ne $co);
      $err .= "b$ct:$tb," if($tb < $cb or $tb > $ce);
      $err .= "e$ct:$te," if($te < $cb or $te > $ce);
      $xb= $tb if($xb==0 or $tb<$xb);
      $xe= $te if($te>$xe);  $xo=$to unless($xo); $xr=$tr unless($xr);
      $cdsl += 1+$te-$tb if($ct eq "c");
    }
  }
  
  my $errcds="";
  $aabigger=  $aasize*3 - $cdsl; # smaller? *stop -1?
  ## lots of 'cdsoffby:-3,' 30k/58k -- is this a stop codon off/on bug?
  ## cdsl includes stop codon; aasize may not.
  
  unless($aabigger == 0) {
    if($CDS_OFFBY_FIX and ($aabigger == -1 or $aabigger == -2) ) {
      my ($cold,$ce)=(0,0); 
      if($co eq "-") {  $ce= $CDS->[0]; $cold=[@$ce]; $ce->[3] -= $aabigger; }
      else { $ce= $CDS->[-1];  $cold=[@$ce]; $ce->[4] += $aabigger; }
      $err .="cdsoffbyfix:$aabigger,";
      #?? report changed CDS rows for gff update ?? need old + new?
      if($ce and $DEBUG) {
       $errcds ="\n";
       $errcds .= "#cdsold:" . join(",",@$cold).";cdsoffby=$aabigger\n";  
       $errcds .= "#cdsnew:" . join(",",@$ce)."\n";  
       }
    } else {
      $err .= "cdsoffby:$aabigger,";
    }
  }
  if($err) { $err="mrna $cr:$cb-$ce:$co ne exons $xr:$xb-$xe:$xo; $err$errcds"; }
  return ($err,$aabigger);
}


sub filter_gff
{
  my($inh, $thechr)= @_;
  
  my($ng,$nx,$nr,$nsame,$nhit,$nwarn,$ndrop)= (0) x 10;
  my @res=(0) x 3; my @res1;
  my @generec=(); my @geneall=(); my @defergene=();
  my @geneother=();
  my ($parentgene,$gid); # 
  
  while(<$inh>) {
    next unless(/^\w/);
    
    my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr)= split"\t";

    my $reforig= $ref;
    $ref= chrname($ref); #FIXME:  rename before test now
    next if($ref eq "skip");
    
    # next if($thechr and $ref ne $thechr); ## thechr opt should be BEFORE rename
    # ** FIXME; now thechr, chrlist AFTER renamed
    if($thechr) {
      if($thechr eq "chrlist") {
        next unless($chrhandle{$ref} or $chrhandle{'all'} or ($chrhandle{'other'} and not $chrmain{$ref}) );
      } else {
        next unless($thechr eq $ref);
      }  
    }
    
    $nr++; chomp($tattr);
    
    my($tid,$pid,$id); 
    if($tattr =~ m/\bID=([^;]+)/) { $id=$1; $tid=$id unless($typ =~ /^($exontypes)$/); }
             #was $tid=$id if($typ =~ /^($mrnatypes)$/);
    if($tattr =~ m/\bParent=([^;]+)/) {  $pid=$1; 
      $tid=$pid unless($tid and ($typ =~ /^($mrnatypes)$/));  }
    unless(defined $tid) { $ng++; $tid = "N".$ng; }

    my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$tid]; 

    if($typ =~ /^gene$/) {  
      ## $rloc->[9]= $id;
      if($gid and $tid ne $gid) { # defer to mrna : NO, bad, changed gid
        if(@generec) { push( @geneall, [@generec], [@geneother]); @generec=@geneother=(); }
        
        # single case mixup: t2 follows separated by other locus
        #   FIXME: Thecc1EG047019t2 is out of order FOLLOWING Thecc1EG047019t1 (correct, renamed from Thecc1EG047020t2)
        if($changelist{$gid} and $changelist{$gid} =~ /^defer/i ) {
          @defergene= @geneall; @geneall=();
        } 
        
        @res1= processallgene( @geneall) if(@geneall); @geneall=();
        if(@res1) { for my $j (0..$#res1) { $res[$j]+=$res1[$j]; } @res1=(); }
      }
      
      $parentgene= $rloc; $gid= $tid;
      
    }  elsif($typ =~ /^($mrnatypes)$/) {  
    
      # FIXME here? collect all altnum of same gene before processgene().
      # .. check all tr fit in genespan and partials include gene end partials.
      my ($pgene)= $tattr =~ m/gene=([^;\s]+)/;  #  problem cases renamed mRNA have old gene=
      my $pgenet= $tid; $pgenet =~ s/t\d+$//;    #  check agreement  
      
      if( $changelist{$tid} ) {
       if( $changelist{$tid} =~ /^gene\s(\w+)/i ) { my $ng=$1; 
        $pgenet= $pgene= $ng; if($gid ne $pgene) { $parentgene=undef; }  # out of order fix 2; is this ok?
        } 
       elsif( $changelist{$tid} =~ /^rename\s(\w+)/i ) { my $nt=$1; my $ng=$nt; $ng =~ s/t\d+$//;
        $pgenet= $pgene= $ng; if($gid ne $pgene) { $parentgene=undef; }  
        } 
      }        
     
      unless($pgene) { $pgene=$pgenet; } # add to tattr ?
      elsif($pgenet ne $pgene) {
        if($pgenet eq $gid) { 
          warno "#WARN: $tid has wrong gene=$pgene; should be $pgenet\n";
          $pgene= $pgenet; $tattr =~ s/gene=([^;\s]+)/gene=$pgenet/;
        } elsif($pgene eq $gid) { # should not be here .. but did get
          warno "#WARN: $tid is wrong? gene=$pgene; not $pgenet\n";
        }
      }
      
      if($pgene and $gid and $pgene ne $gid) { #  wrong, gid updated already ..
        # now this can be out-of-order mrna, where gene follows.
        @res1= processallgene( @geneall) if(@geneall); @geneall=();
        $gid=$pgene; #?
      } elsif($pgene eq $gid) {
        push( @geneall, [@generec], [@geneother]) if(@generec); @generec=@geneother=();
      }

      
      # allow gene and mRNA types .. and/or keep other types in generec
      if(@generec or @geneother) { @res1= processgene( \@generec, \@geneother);  }
      if(@res1) { for my $j (0..$#res1) { $res[$j]+=$res1[$j]; } @res1=(); }

      if( $changelist{$tid} ) { # handle here or where?
        # geneid  action:(drop,rename,retype..) newid ... grep valid actions?
        if( $changelist{$tid} =~ /^drop/i ) { $ndrop++; @generec = (); next; } # can include gene or altnum drop

        elsif( $changelist{$tid} =~ /^strand\s([+-])/i ) { my $onew=$1;
          $rloc->[6]=$onew; #[$ref,$src,$typ,$tb,$te,$tp,$onew,$tph,$tattr,$tid]; 
          warno "#WARN: $tid strandfix: $onew \n";
          } # reset strand .. should do in.gff
          
        elsif( $changelist{$tid} =~ /^spanfix\s(\S+)/i ) { my $loc=$1;
          my($nr,$nb,$ne,$no)= $loc =~ m/^(\w+):(\d+)[.-]+(\d+)[:]?([+-]?)/;
          # fix span only, not strand/ref; presumably to match exons span
          if($ne>0 and $nb>0) {  $rloc->[3]=$nb;  $rloc->[4]=$ne; }
          warno "#WARN: $tid spanfix: $loc \n";
          } # messier than reset strand .. should do in.gff

        elsif( $changelist{$tid} =~ /^rename\s(\w+)/i ) { my $ntid=$1;  # relocus?; newid ..
          #see above # my $ng=$nt; $ng =~ s/t\d+$//;
          # .. fixme, failed to make gene for Thecc1EG047111t1
          $gid= $pgenet; # is this answer? NO
          # make new parentgene?
          my $gattr="ID=$gid";
          $parentgene= [$ref,$src,"gene",$tb,$te,$tp,$to,$tph,$gattr,$gid];

          warno "#WARN: $tid rename $ntid\n";
          $tattr =~ s/ID=$tid/ID=$ntid/;
          $tattr =~ s/gene=([^;\s]+)/gene=$pgenet/;
          $rloc->[9]=$ntid; #[$ref,$src,$typ,$tb,$te,$tp,$onew,$tph,$tattr,$tid]; 
          $rloc->[8]=$tattr;
          } 
## new changelist MisMatchAA = skip my pep, let ncbi translate. see putCDS
##      elsif( $changelist{$tid} =~ /^gene\s(\w+)/i ) { } #see above
##      elsif( $changelist{$tid} =~ /^retype/i ) {   } # new type ..
##      elsif( $changelist{$tid} =~ /^relocate/i ) {   } # defer now; save out-of-order mrna/exons ; precedes gene.
        elsif( $changelist{$tid} =~ /^defer/i ) {  # deferto
          if(@defergene) {  @geneall= (@defergene,@geneall); @defergene=(); $changelist{$pgene} =~ s/^defer/diddefer/; }
        }
     }
      
      $ng++;
      @geneother= ();
      @generec = ($rloc); # maybe best store as array of [gffcols], sorted by type-loc
      push @generec, $parentgene if( $parentgene); # wrong gene here???
      $parentgene=undef; # only 1st
      
    } elsif($typ =~ /^($exontypes)$/) {

      if($changelist{$tid}) {
        if( $changelist{$tid} =~ /^drop/i ) { next; }
        elsif( $changelist{$tid} =~ /^strand\s([+-])/i ) { my $onew=$1;
          $rloc->[6]=$onew; # [$ref,$src,$typ,$tb,$te,$tp,$onew,$tph,$tattr,$tid]; 
          # warno "#WARN: fix strand: $onew $tid\n";
        } elsif( $changelist{$tid} =~ /^rename\s(\w+)/i ) { my $ntid=$1;   
          # warno "#WARN: $tid rename $ntid\n";
          $tattr =~ s/Parent=$tid/Parent=$ntid/; $pid= $ntid; 
          $rloc->[9]=$ntid; #[$ref,$src,$typ,$tb,$te,$tp,$onew,$tph,$tattr,$tid]; 
          $rloc->[8]=$tattr;
        } 
      } 

    #... all skipped exon are from drops in changelist. no warning?
    #WARN: skipped exon, no parent exon=Thecc1EG007388t2
    #WARN: skipped exon, no parent CDS=Thecc1EG011170t3
       
      unless(@generec) { warno "#NOTE: skipped exon, no parent $typ=$pid\n"; }  
      else {
      my $rid= $generec[0]->[9];
      if($rid ne $pid and not $IgnoreOutOfOrder) { warno "#ERROR: out of order gene record: mRNA=$rid, $typ=$pid\n" if($nwarn++ < 9); } # limit warns
      else { push @generec, $rloc; $nx++; } 
      }
      
    } elsif( ! $passtypes or "$typ.$src" =~ m/$passtypes/ ) {
      push @geneother, $rloc;  # may be independent of gene/transcript .. dont require generec
    }
  }
  
  @res1=();
  if(@generec and @geneall) { push( @geneall, [@generec], [@geneother]);  @generec=@geneother=(); }
  if(@geneall) { @res1= processallgene( @geneall); @geneall=(); }
  if(@generec or @geneother) { @res1= processgene( \@generec, \@geneother); }
  if(@res1) { for my $j (0..$#res1) { $res[$j]+=$res1[$j]; } @res1=(); }

  return (@res,$ndrop); # $ndrop   return ($ng,$nx,$nsame,$nhit);
}


=item putgenome

  put genome.fsa, genome.qvl (qual), genes.pep (unless PEPfromGFF)
  
** fixed BUG in scaf selection:
$evig/evigene/scripts/evigene2genbanktbl.pl -conf genes/evigene_cacao3_gbsubmit.conf -debug \
-in genes/pub3h/cacao11genes_pub3h.good.gff.gz -out submit/pub3h.tbl -chr 'scaffold_6,scaffold_8' -log

dgbook3% grep '^>' submit/pub3h.scaffold_8.fsa | wc
     703     703    9793
.. other mixup ..

=cut

sub putgenome 
{
	my $self= shift; # == \%config
	my( $tblpath, $thechr)= @_;
	warno "#NOTE: putgenome( $tblpath, $thechr )\n" if $DEBUG;

#  613 ERROR:   SEQ_DESCR.NoOrgFound					# 701 scaffolds, 613 w/o genes/annotation, 88 in .tbl with source/organism annot
#  see putsource_noannotchrs(), need to write .tbl entry for each of these w/o genes.
	
  my ($tblname, $subdir, $tblsuf) = File::Basename::fileparse($tblpath, qr/\.[^.]*/);
  #updated.............................
  #.. need to process genome per chr file; pre-separate? or do here
  my $fastachr  = $genome  || $evidence{'genome'};  # die unless(-f $fastachr);
  my $fastaqual = $evidence{'genomequal'} || "$genome.qual";  
  my $fastapep  = $proteins || $evidence{'proteins'};  

  my $gbsubfsa = ($thechr eq "chrlist") ? undef : "$subdir$tblname.fsa";
  my $gbsubqvl = ($thechr eq "chrlist") ? undef : "$subdir$tblname.qvl"; # qual
  my $gbsubpep = "$subdir$tblname.pep";

	warno "#NOTE: putgenome( $fastachr, $gbsubfsa )\n" if $DEBUG;
	warno "#NOTE: putgenome( $fastaqual, $gbsubqvl )\n" if $DEBUG;
	warno "#NOTE: putgenome( $fastapep, $gbsubpep )\n" if $DEBUG;

# change thechr to @chrlist and split to @tblname list
  if(1 or $didchrrename or $thechr) {
    my($pin,$pout)= openio($fastachr,$gbsubfsa);
    my $nop=0;
    while(<$pin>) { 
      if(/^>(\S+)/) { my $t=$1; 
        my $p= chrname( $t);  #FIXME:  rename before test now
        if($p eq "skip") { $nop=1; }
        elsif($thechr eq "chrlist") { $nop= setrefout( $p) ? 0 : 1; $pout=$outh; }
        else { $nop= ($thechr and $p ne $thechr) ? 1 : 0;  } ## thechr opt should be BEFORE rename
        #FIXME:  rename before test now# my $p= chrname( $t); 
        s/>$t/>$p/; 
      } 
      print $pout $_ unless($nop);
    } close($pout) if(ref $pout); close($pin);
    
    if( -f $fastaqual) {
      $nop=0;
      ($pin,$pout)= openio($fastaqual,$gbsubqvl); 
        # ** pout is ignored for chrlist; 
        # ** setrefout() peph=qvl is used instead of outh=fsa
      while(<$pin>) { 
        if(/^>(\S+)/) { my $t=$1; 
          my $p= chrname( $t);  #FIXME:  rename before test now
          if($p eq "skip") { $nop=1; }
          elsif($thechr eq "chrlist") { $nop= setrefout( $p) ? 0 : 1; $pout=$peph; }
          else { $nop= ($thechr and $p ne $thechr) ? 1 : 0;  } ## thechr opt should be BEFORE rename
          s/>$t/>$p/; 
        } 
        print $pout $_ unless($nop); # undefined value
      } close($pout) if(ref $pout); close($pin);
    }
     
  } else {
  	my $curdir= $ENV{'PWD'};  #?? not safe?
    $fastachr= "$curdir/$fastachr" if(-f "$curdir/$fastachr");
    symlink( $fastachr, $gbsubfsa); 
  }
  
  # redo here: option pull aa from gff.mRNA protein=
  unless( $PutPEPfromFile ) { #was $PEPfromGFF  # and -f $gbsubpep
    ## new changelist MisMatchAA = skip my pep, let ncbi translate. see putCDS
    my($pin,$pout)= openio($fastapep,$gbsubpep);
    while(<$pin>) { 
      if(/^>(\S+)/) { my $t=$1; 
        # my $p= reformatProteinID( $t);  # need lcl| or gnl| ??  either way all .pep not found.
        my $p= reformat_myid( reformatProteinID( $t));  # need >lcl|protid or >gnl|db|protid ??
        s/>$t.*$/>$p/; } 
      print $pout $_ if(/\S/); # are blank lines a problem for tbl2asn?
    } close($pout); close($pin);
  }

}  




=item tbl2asn

  foreach tblouts/chrs do  tbl2asn( $tblname);
  TBL2ASNready controls if tbl2asn is executed; defer
  
=cut

sub tbl2asn 
{
	my $self= shift; # == \%config
	my( $tblpath, $thechr)= @_;
	warno "#NOTE: tbl2asn( $tblpath, $thechr )\n" if $DEBUG;
	
  my ($tblname, $subdir, $tblsuf) = File::Basename::fileparse($tblpath, qr/\.[^.]*/);
  ## fileparse: subdir/ not subdir
  my $log="$subdir/log.tbl2asn.$tblname"; # see also logout
  
	my $opts   = $programs{'tbl2asnopts'};
  my $tbl2asn= setprogram("tbl2asn","","", not TBL2ASNready);
  # tbl2asn  = ncbi/bin/tbl2asn  
  # tbl2asnopts = -t template.sbt -V vb -p ./ 
  unless($opts =~ s/\-p\s*\S+/-p $subdir/) { $opts .= " -p $subdir"; }

  my $species= $config{general}->{species};
  $opts .= " -n '$species'" if($species and $opts !~ /\-n/);

  my $db= $DBID;
  $opts .= " -C $db" if($db and $opts !~ /\-C/);

  my $discrep="$subdir$tblname.discrep";  
  $opts =~ s/\-Z\s*\S+/-Z $discrep/;
  
  # tbl2asn 18.1   arguments:
  # -C  Genome Center Tag [String] 
  # -n  Organism Name [String]
  # -V  Verification : v=normal; b=genbank flat

  my $sbt= $self->{submit_template}->{doc};
  my $sbtpath= "$subdir$tblname.sbt"; ## WRONG:$self->{submit_template}->{path} || "template.sbt" ;
  if($opts =~ /\-t/ and $sbt) {
    unless($opts =~ s,\-t\s*(\S+),-t $sbtpath,) {} # NO: $opts.= " -t $sbtpath";
    ## die "ERROR: tbl2asn missing submit_template $sbtpath" unless($sbt and $sbtpath);
    open(F,">$sbtpath") or die $sbtpath; print F $sbt,"\n"; close(F);
  }
  
  my $ok=-99;
  warno("#NOTE: $tbl2asn $opts\n") if $DEBUG;
      ## "# using files: $tblpath, $gbsubfsa, $gbsubpep \n") if $DEBUG;
  # chdir($subdir);  # NOT? yes, tbl2asn dont work well for paths
  if(TBL2ASNready) { $ok= system("$tbl2asn $opts > $log 2>&1"); }
  else { warno "#NOTE: NOT RUN: $tbl2asn $opts\n"; }
  # chdir($curdir);
  return $ok;
}


#........................................



sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }

#   my @generec= sort _sortgene @$generecIN;
#  do 3-level sort: mRNA rec sorted by type; gene>mRNA sorted by alt-tr num?; chr>genes sorted by location
sub _sortgene  
{
  # rec == ($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr)
  # ($b->[2] cmp $a->[2]) # typ: reverse order by mRNA, exon, CDS  
  ##my($ta,$tb)= map{  (/gene/)?1:(/mRNA/)?2:(/CDS|exon/)?3:4; } ($a->[2],$b->[2]);
  my($ta,$tb)= map{ (/gene/)?1:(/mRNA/)?2:(/exon/)?3:(/CDS/)?4:5; } ($a->[2],$b->[2]); 
  return ($ta <=> $tb)
      || ($a->[0] cmp $b->[0]) # ref : should all be same for gene record
      || ($a->[3] <=> $b->[3]) # begin
      || ($a->[4] <=> $b->[4]) # end
      ;
}

sub _sortaltid {
  my($ta,$tb)= map{ (m/t(\d+)/) ? $1 : 0  } ($a,$b); 
  return ($ta <=> $tb || $a cmp $b);
}




sub evigene_config {
  my($cfile, $addoptions)= @_;
  my $ctype=0;

  use constant{ kEVFILE => 'evidence', kEVOPT => 'evoption', kANOPT => 'anoption', 
                kEVPROG => 'programs', kPUBOPT => 'pubopt', kEVGENES => 'geneset', };
  
  if($cfile) { #  and -f $cfile
    open(F,$cfile) or die "ERROR reading config: $cfile";
  
    my @CONFIG= <F>;
    push @CONFIG, @$addoptions if(ref $addoptions);
    
    my ($lastkey, $lastval);
    # while(<F>) 
    foreach (@CONFIG)
    {  # key => value
      s/^\s+//;  s/\#.*$//; s/\s+$//; # end of line comments 

    ## need now to handle continuation lines, end with \
    ## or prefix with "+" << use this

      my($key,$val);
      if(/^[+\.]/ and $lastkey) { $key= $lastkey; s/^.//; $val= $lastval.$_; }
      elsif($lastval =~ m,\\$, and $lastkey) { $key= $lastkey; $lastval=~s,\\$,,; $val= $lastval.$_; }
      else { ($key,$val)= split(/[=\s]+/, $_, 2); }
      
      next unless($key =~ /^\w/); 
      $val =~ s/\s*$//;  $val =~ s/\\n/\n/g; $val =~ s/\\t/\t/g; # dequote tabs,newlines
      # allow for val == regex in match/m or subs/s;  
      #  names:
      #   cutdbx = s/\((InterPro|TAIR):[\w\.-]+\)//g ?
      #   isgeneid = m/^(Os|At|AT)\d{1,2}g\d{3,}/
      # verbose "# config k=v: $key=$val";

# revise to parse '^section:' into separate hash
      if($key =~ s/:$//) { $ctype=$key; }
      
      if($key =~ /^evidence/) { $ctype= kEVFILE; }
#       elsif($key =~ /^evoption/) { $ctype= kEVOPT; }
#       elsif($key =~ /^anoption/) { $ctype= kANOPT; }
#       elsif($key =~ /^geneset/) { $ctype= kEVGENES; }
      elsif($key =~ /^pubopt/) { $ctype= kPUBOPT; }
      elsif($key =~ /^program/) { $ctype= kEVPROG; }
      elsif($key =~ /^end$/) { $ctype= 0; }
      
      elsif($ctype eq kEVFILE ) { $evidence{$key}= $val; } # /gff/ was bad for geneset
      elsif($ctype == 0 and $val =~ /\.gff/) { $evidence{$key}= $val; } 
#       elsif($ctype eq kEVOPT) { $evaluate_options{$key}= $val; } 
#       elsif($ctype eq kANOPT ) { $annotate_options{$key}= $val; } 
#      elsif($ctype eq kEVGENES) { $geneset{$key}= $val; }
      elsif($ctype eq kPUBOPT ) { $public_options{$key}= $val; } 
      elsif($ctype eq kEVPROG) { $programs{$key}= $val; }

      # generic keys: name date genome .. other?
      if($key =~ /\w/ and $val =~ /\w/) { 
        my $ogroup= $ctype || "general";
        #which?  $config{$ogroup}{$key}= $val;
        $config{$ogroup}->{$key}= $val; # this one
      }
      
      # also for now : overlap, other progs
      ($lastkey, $lastval)=($key, $val);
      
    } close(F);
  }
  
#   die "ERROR: missing overlapfilter program: $overlapfilter"
#       unless( -x $overlapfilter);

}

__END__

=item tbl2asn fixme: 1chr file set, or -a s option

The -a s flag tells tbl2asn to package the multiple FASTA components as a set
of unrelated sequences.  This accommodates users who create a single file
instead of one file per sequence.  

=item more opts tbl2asn

To process a set of chromosomes, sets of .fsa and .tbl files (along with
optional .src, .pep, .rna, and .qvl files) are placed into a source
directory. 
  chr01.fsa
  chr01.tbl
  chr02.fsa

The -g flag causes tbl2asn to generate a genomic product set.  Within the
set, the products of each related mRNA and CDS are packaged together in an
internal nuc-prot set.

=item finally some success

dgbook3% cat genes/pub3h/pub3h.gff | grep '^scaffold_10r' | \
 $evig/evigene/scripts/evigene2genbanktbl.pl -conf genes/evigene_cacao3_gbsubmit.conf -debug \
 -in stdin -out submit/pub3h_sc10.tbl

  -- dont use opts= -t name.sbt 
  -- do use tablename.sbt not template.sbt
 
# DONE gff2tbl: 1 3301 4204 4204 0 
# tbl2asn( submit/pub3h_sc10.tbl )
# ./mac.tbl2asn -g -V vb -p submit/ 
# using files: submit/pub3h_sc10.tbl, submit/pub3h_sc10.fsa, submit/pub3h_sc10.pep 
# DONE tbl2asn: 0 


=item tbl2asm errs

cat genes/pub3h/pub3h.gff | grep '^scaffold_10r' | \
$evig/evigene/scripts/evigene2genbanktbl.pl -conf genes/evigene_cacao3_gbsubmit.conf -debug \
-in stdin -out submit/pub3h_sc10.tbl
# DONE gff2tbl: 1 3301 4204 4204 0 

# tbl2asn( submit/pub3h_sc10.tbl )
# ./mac.tbl2asn -t template.sbt -V vb -p submit/ 
# using files: submit/pub3h_sc10.tbl, submit//pub3h_sc10.fsa, submit//pub3h_sc10.pep 
[tbl2asn 18.1] Unable to read required template file
# DONE tbl2asn: 256 

** mac.tbl2asn doesnt do directories : cant find -t submit/template.sbt, but can in same dir

.. running, got this at start:

>>> FIXME: these are true IDs, dang gb wants new ID for proteins:  t>p in .pep file
dgbook3% less log.t1
[tbl2asn 18.1] Unable to find protein sequence Thecc1EG000001t1
[tbl2asn 18.1] Unable to find protein sequence Thecc1EG000002t1
[tbl2asn 18.1] Unable to find protein sequence Thecc1EG000002t2
   ... all are listed; what? all are in .pep

[tbl2asn 18.1] Flatfile pub3h_sc10

[tbl2asn 18.1] Validating pub3h_sc10

 BAD CODE ...
mac.tbl2asn(6468,0xa0745540) malloc: *** error for object 0x7eb8340: incorrect checksum for freed object - object was probably modified after being freed.
*** set a breakpoint in malloc_error_break to debug
mac.tbl2asn(6468,0xa0745540) malloc: *** error for object 0x7f1ce60: incorrect checksum for freed object - object was probably modified after being freed.
*** set a breakpoint in malloc_error_break to debug
   
=cut

=item submit6 .val

1160 [SEQ_INST.InternalNsInSeqRaw]
1083 [SEQ_FEAT.IllegalDbXref]    # reduced with fake dbxref
 474 [SEQ_FEAT.UnnecessaryGeneXref]  # from where now
 352 [SEQ_FEAT.CDSwithMultipleMRNAs]
 309 [SEQ_FEAT.NotSpliceConsensusAcceptor]
 275 [SEQ_FEAT.NotSpliceConsensusDonor]
 202 [SEQ_FEAT.DuplicateFeat]
 156 [SEQ_FEAT.PartialProblem]
  42 [SEQ_FEAT.MultipleGeneOverlap]
  41 [SEQ_FEAT.InvalidInferenceValue]
  10 [SEQ_FEAT.CDSmRNArange]
   7 [SEQ_FEAT.UnqualifiedException]
   3 [SEQ_FEAT.GeneXrefStrandProblem]
   2 [SEQ_FEAT.CDSmRNAmismatch]
   2 [SEQ_FEAT.MissingGeneXref]
   2 [SEQ_FEAT.NoStop]
   1 [SEQ_FEAT.CDSwithNoMRNAOverlap]
   1 [SEQ_FEAT.ShortIntron]


=item submit1 .val report sum:

cat submit1/pub3h_sc10.val | perl -ne'@v=split; print "$v[2]\n";' | sort | uniq -c | sort -k1,1nr 10792 [SEQ_FEAT.IllegalDbXref]

8408 [SEQ_DESCR.NoOrgFound]        # what ? need -opt species ?

3322 [SEQ_FEAT.InvalidInferenceValue]   # db_xref critic
1219 [SEQ_INST.InternalNsInSeqRaw]      # genome asm problem

 774 [SEQ_FEAT.PartialProblem]          # bug in <> usage or end not at splice site
 587 [SEQ_FEAT.NotSpliceConsensusAcceptor]  # what?
 483 [SEQ_FEAT.NotSpliceConsensusDonor]
 
 379 [SEQ_FEAT.CDSwithMultipleMRNAs]   #?
 320 [SEQ_FEAT.StartCodon]    
 204 [SEQ_FEAT.DuplicateFeat]        #? most form alt-tr ; doesnt handle those?
 194 [SEQ_FEAT.PartialsInconsistent]
 161 [SEQ_FEAT.MultipleGeneOverlap]
 156 [SEQ_INST.BadProteinStart]
  41 [SEQ_FEAT.InternalStop]
  33 [SEQ_FEAT.NoStop]
  20 [SEQ_INST.TerminalNs]
  17 [SEQ_FEAT.GenesInconsistent]
  15 [SEQ_INST.ShortSeq]
  13 [SEQ_FEAT.CDSmRNArange]
  13 [SEQ_INST.StopInProtein]
   9 [SEQ_FEAT.GeneXrefStrandProblem]
   8 [SEQ_FEAT.SeqLocOrder]
   7 [SEQ_FEAT.TransLen]
   4 [SEQ_FEAT.CDSmRNAmismatch]
   4 [SEQ_FEAT.TranscriptLen]
   2 [SEQ_FEAT.CDShasTooManyXs]
   2 [SEQ_FEAT.MissingGeneXref]
   2 [SEQ_INST.HighNContentPercent]
   1 [SEQ_FEAT.CDSwithNoMRNAOverlap]
   1 [SEQ_FEAT.ShortIntron]     # asmrna needs fixing
   1 [SEQ_INST.TrailingX]


=item tbl example1

Sample output so far ....
  data/genomes/Anopheles_gambiae/anogam_20080511/genbanktbl/*tbl

>Features       NC_004818       AnoGambia_20080511
1       24393108        chromosome
        locus_tag       NC_004818
        species Anopheles_gambiae str. PEST


8966    4779    gene
        locus_tag       AgaP_AGAP000003
        db_xref VectorBase:AGAP000003
        locus_tag       AgaP_AGAP000003

5263    4779    mRNA
8966    8696
        transcript_id   AgaP_AGAP000003.t01
        note    AGAP000003-RA coding for AGAP000003-PA
        db_xref GI:158289427,VectorBase:AGAP000003,VectorBase:AGAP000003-RA
        locus_tag       AgaP_AGAP000003
        product AGAP000003-RA


=item tbl example2

>Features	NC_004353	DrosMelGb_20080512
1	1351857	source
	db_xref	NC_004353
	name	NC_004353
	organism	Drosophila melanogaster

...
64403   53434   gene
        locus_tag       plexB
        gene    plexB
        note    plexin B; synonyms: PlexB, plex, unnamed
        db_xref FLYBASE:FBgn0025740
        db_xref GeneID:43766
        old_gene        plexB
        old_locus_tag   Dmel_CG17245
        cyt_map 102A1-102A1

64403   63540   mRNA
61911   57142
57083   56500
53999   53817
53751   53434
        transcript_id   plexB.t01
        product plexB.t01
        note    plexin B; synonyms: PlexB, plex, unnamed
        locus_tag       plexB
        gene    plexB
        db_xref FLYBASE:FBgn0025740
        db_xref GeneID:43766
        old_gene        plexB
        old_locus_tag   Dmel_CG17245
        old_product     plexin B CG17245-RA
        old_transcript_id       NM_079877.2

64050   63540   CDS
61911   57142
57083   56500
53999   53817
53751   53644
        locus_tag       plexB
        gene    plexB
        protein_id      plexB.p01
        product plexB.p01
        note    CG17245 gene product from transcript CG17245-RA
        db_xref FLYBASE:FBgn0025740
        db_xref GeneID:43766
        old_gene        plexB
        old_locus_tag   Dmel_CG17245
        old_product     plexin B CG17245-PA
        old_protein_id  NP_524616.2

=cut
