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
      
=cut

use FindBin;
use lib ("$FindBin::Bin","$FindBin::Bin/prot/"); # == evigene/scripts/

use strict;
use Getopt::Long;
use File::Basename;
use protein_names; # Evigene nameclean() etc

use constant VERSION  => '20130326';  

use constant CHR_RENAME   => 1; # on, works right
use constant TBL2ASNready => 0; # DONT run tbl2asn yet, leave to operator

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
my $FAWIDTH = 60;
my $CDS_OFFBY_FIX=1; # config 

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
  "chromosome|ref=s", \$thischr,
  #....
  "version|MSRC=s", \$MySRC,
  "cadd=s", \@configadd,
  "notbl2asn!", \$notbl2asn, 
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
$PEPfromGFF = 1 if($evidence{proteins} eq "gff"); # || $PEPfromGFF;
#forget this# $ADDGENE      = (defined $config{attr_gene}->{skip}) ? ! $config{attr_gene}->{skip} : $ADDGENE;
my $DBID= $config{general}->{databaseid} || $config{general}->{dbid} || "MyDBID"; # make global
$MySRC= $DBID; #?? yes or no; now only used for myxref eg paralogs
my $LOCUSTAG= $config{general}->{locus_tag} || $DBID;  
$MyIDPrefix= $public_options{publicid} || $MyIDPrefix;  $MyIDPrefix=~s/\d+$//;

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
  puevd => "PASA:2",
  pasa  => "PASA:2",
  vel   => "Velvet:11" ,
);
if( my $progvers= $programs{progvers}) {
  my @pg= split/[,\s]+/, $progvers; 
  map{ my($k,$v)=split/=/,$_; $progvers{$k}=$v if($v); } @pg;
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

sub warno { 
  my $c=""; #($_[0] =~ /^#|^\d/) ? shift @_ : 3; 
  #$c=($c=~/#/)?$c : ($c==1)?"#ERROR: " : ($c==2)?"#WARN: " : "#NOTE: "; 
  if(ref $logh) { print $logh $c, @_; } else { warn $c, @_; } 
  }
sub tblsuf { my($na,$pt,$suf)= @_; if($na) { $na =~ s,\.\w+$,,; return join(".",$na,$pt,$suf); } $na; }

if( $evidence{genomesplit} ) { # handle other != chrmain
  my @csplit=  split /[,\s]+/, $evidence{genomesplit}; # now renamed chr
  foreach my $c (@csplit) { next if($c =~ /other/); $chrmain{$c}=1; }
}

$thischr= $evidence{genomesplit} unless($thischr);
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
  (undef, $peph) = openio( undef, $pepout1) unless($notbl2asn);
  
  $chrhandle{$achr}{'outh'}= $outh;
  $chrhandle{$achr}{'logh'}= $logh;
  $chrhandle{$achr}{'peph'}= $peph;
}

#? allow 2+ ingff? e.g. genes + transposons, result will be unsorted.tbl
if(%chrhandle) {
  ($inh, undef)= openio($ingff, undef); 
  @countall= filter_gff($inh, "chrlist");
  putsource_noannotchrs(); # FIX 2012feb
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


  #FIXME: 
sub putsource_noannotchrs 
{
  foreach my $chr (sort keys %chrattr) {
    next if($chrhasannot{$chr});
    if(setrefout($chr)) { # output handle for chrsplitting; checks chrhandle valid chr
      putsource($chr);
    }
  }
}



sub reformatProteinID { my $t=shift; $t=~s/t(\d+)$/p$1/; return $t; } 

sub reformat_myid
{
  my($id)= @_;
  my $db= $DBID;
  unless($db) { ($db)= $id =~ m/^(\w+[A-Za-z])\d\d/; }
  ## splitgene: FIXME here?? add hash lookup for few splitgene IDs; protein_id|transcript_id need split part-id tags, not gene
  if(SPLITGENE_FIX) { ## bad for protein_id: p1>t1
    ## FIXED BELOW# my $tid=$id; $tid =~ s/p(\d+)$/t$1/; # or do below for splitgene{pid} ?
    if($splitgene{$id}) { my $partnum= $splitgene{$id}; return "gnl|$db|$id.$partnum";}
  }
  return "gnl|$db|$id";  # which? << prefered in one ncbi doc
  #?# return "lcl|$id"; # which?
}


sub reformatval
{
  my($typ, $id, $key, $oldkey, $val)= @_;
  # my $rekey=$key;
  my @val; #?? = undef; << BAD
  
  if($key =~ /^(locus_tag)$/) {
    # $val =~ s/([A-Z])(\d\d+)/${1}_${2}/; ##messy# 
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

  ##use constant PARALOG_SPECIAL => 0;
  
  if($oldkey =~ /paralog/) { 
      # check if $changelist{$tid} =~ /^drop/ 
    my ($pi)= ($val =~ m/(\d+)%/ ) ? $1 : 69; 
    my ($did)= $val =~ m/($MyIDPrefix\w+)/;
    return () unless($did and $pi >= $MIN_IDENTITY); 
    ## return () if(not $did or $val =~ /paralog.na|^na|paralog.\d|^\d/ or $pi < $MIN_IDENTITY); 
    return () if ($changelist{$did} =~ /^drop/);
    # note    paralog=Thecc1EG034652t1,83% < drop % .. id:pct
    # fixme:    note    gnl|CacaoGD|paralog=Thecc1EG005594t1
    # fixme2:   note    paralog  Thecc1EG034062t1 gnl|CacaoGD|Evigene.of
    # fixme3:   note    paralog of 47% gnl|CacaoGD|
    
    if($key =~ /^note/) {    # only for CDS:note,db_xref but for mRNA:inference, do like ortholog
      $val =~ s/,(\d+)%$//;    # drop :pid # $val =~ s/%//; $val =~ s/,/:/; 
      $val =~ s/$did//;
      @val=("$key\t$val " . reformat_myid($did));
      ### fixme: recode as comment not db_xref
      # TURN OFF for now#  push @val, "db_xref\t$DBXOTHER:$did" unless($DBXOTHER =~ /MISSING/);
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
    $pi= ($pi>99)?"1.00":"0.$pi"; ## sprintf "%.2f", $pi/100; # change from % to prop; 1.0 ..
    $val= "similar to RNA sequence, EST $pi coverage"; #as note not inference; gb hates %%%
    ## return ("note\t$val");
    
  } elsif($oldkey =~ /quality/) {  # reformat for genbank
     ## quality=Class:Strong,Express:Strong,Homology:ParalogStrong,Intron:Strong,Protein:complete
     ## >> note    quality=Class:Strong,Express:Weak,Homology:ParalogStrong,Intron:None,Protein:complete
     while($val =~ m/(\w+):(\w+)/g) { 
      my($k,$v)=($1,$2); my $v1=$v; $v1=~s/(\w)(Strong|Medium|Weak)/$1 $2/;  
      $val =~ s/$k:$v/\L$k is \L$v1/; 
      } 
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
      if( $DBXREF_RECODE{$d} ) { $d= $DBXREF_RECODE{$d}; }
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
  
  # recode idprefix: r8caca11 > cacao11; m7vitvi> vitvi
  # leave out '%' char
  if($oldkey =~ /ovrna|trasm/) {
    my @v = grep  /^[a-zA-Z]/, split",", $val;
      ## fixme: CacaoGD mar11f:CGD are bad ovrna: alignment:GMAP:11:TSA:mar1g.mar11f:CGD0022153
      ## bad also:  alignment:GMAP:11:TSA:mar7g.mar11f:AUGepir7p1s1g7t1
      ## bad also:  alignment:GMAP:11:TSA:AUGpie3_p1s_1g7t1
    @v= grep { ! /mar1g.mar11f|mar7g.mar11f|^AUG/ } @v;
    my @pk= sort keys %progvers;
    my @ev= grep /\w/, map {
      my $vv=$_;
      s/^(cgba).\w+:/$1/; 
      s/^\w+://; s/^r8//; s/caca11/cacao11/; # add to config
      $_= "cgba".$_ if(/^(B|L|P1|P2)_g\d/);
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
      # (/vel/) ? "Velvet:11" : (/cuf8/) ? "Cufflinks.08" : (/cuf/) ? "Cufflinks:11" : (/puevd|pasa/) ? "PASA:2" : "GMAP:11";
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
    $attr=~s/Name=$na/Name=$newna/ if($na ne $newna);
  }
  
  my @attr; my %didk;
  foreach my $at (split";",$attr) {
    my($k,$v)= split"=",$at; 
    next if($didk{$k}++ and $k !~ m/Note|Dbxref|Alias/); #?? not always ! eg Note
    
    my $kn= $mapattr->{$k} || "skip";
    next if($kn =~ /^skip/); # make mapattr full set of allowed tags?
    
    my $kval="";
    ($kn,$kval)= split(/[=\s]+/, $kn, 2); 
    
    if($kval) { $v= ($kval=~m/[:=]$/)?"$kval$v":"$kval $v"; } # more?
    # isoform   note=encoded by transcript variant $v
    
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
  if($partial eq "partial") { $partial.="53"; } # both
  # need retype option: if($config{attr_$typ}{type})
  my $typo= $config{"attr_".$typ}->{type} || $typ;
  my $addattr=""; # return to caller
  
  my $gstrand= $loc[0]->[6]; ## my(undef,undef,$gstrand)= split("\t",$loc[0]);
  my @iter= ( 0 .. $#loc );
  @iter= reverse @iter if ($gstrand eq "-");
  my($p5,$p3)= ("<",">"); #?? or this# ($gstrand eq "-") ? (">","<") : ("<",">");
  
  my $first= $iter[0];
  my $last = $iter[-1];
  foreach my $i ( @iter ) { 
    my($start,$stop,$strand,$phase)= @{$loc[$i]}[3,4,6,7];  # split("\t",$loc[$i]);
    ($start,$stop) = ($stop,$start) if ($strand eq "-"); ## ($gstrand < 0);
    if($partial =~ /3/ and $i==$last ) { $stop="$p3$stop"; }
    # if($partial =~ /5/ and $i==$first) { $start="$p5$start"; }
    if($partial =~ /5/ and $i==$first) {  
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
  $pid= reformat_myid( reformatProteinID( $pid)); 
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

sub nameclean_OLD
{
  my($namin, $pi, $locusname)= @_;
  my $lowqualname=undef;
  my ($isunk,$isput,$islike)= (0,0,0);

  #*** ?? IS THIS FAILING .. missing haem, Tumour
  # local $_= $namin; 
  #?? study();
  my $_ = $namin;   # perl 5.9.1

  s/\s+\(\d+%.*\)//; # trailing pctident
  s/\s+sym:\S+//;  # trailing symbol tags
  
  #  horribles die quickly...
  if(/The protein encoded by this gene was identified/) { $_=$NAME_UNK; }
  elsif(/Encodes a close of the Cauliflower OR/) { $_="Orange protein"; }
  
  # typos and britglish > amglish; need config list
  s/Uncharacterised/Uncharacterized/; 
  s/dimerisation/dimerization/;  s/luminium/luminum/g;  # Aluminium
  s/signalling/signaling/; # Two Ls in British English, one in American English. 
  s/onoxygenase/onooxygenase/; # [Mm]onox..
  s/sulphide/sulfide/ig; s/sulphur/sulfur/ig;
  s/Tumour/tumor/ig;  
  s/haemoprotein/hemoprotein/i;  s/\bhaem/heme/ig; # Quinohaemoprotein>Quinohemoprotein
  s/characteris/characters/;  
  s/proteine/protein/;  s/\bcomponenet/\bcomponent/;
  s/ protei$/ protein/; # prior nameclean bug
  
#  # UPDATE: add to species.config some of these nameclean patterns
#   my $nameclean = ($config{nameclean} and ref($config{nameclean})) ? $config{nameclean} : {};
#   foreach my $nk (sort keys %$nameclean) {
#     my $npatt= $nameclean->{$nk};
#     # if($npatt =~ m/^s(.)/) { my $nc=$1; my($ncut,$nto)= split/[$nc]/,$npatt; s/$ncut/$nto/; }
#     my $eok= eval $npatt; # is this ok? prefer s/$nacut/$nato/ ; 
#     # any use for m/$napatt/ ; need specific keys for this, m/$ISGENEID/ =  m/^(Os|At|AT)\d{1,2}g\d{3,}/
#   }
  
  if(s/^TE://i) { }  # $istransposon=1; 
  if(s/^(Predicted|Conserved|Expressed|Novel)\s+//i) { }  # maybe set $isunk?
  if(s/^(Possible|potential|probable)\s+//i) { $isput=1; }   
  s/(homology|similar|Similarity|related) to\s*//i;
  s/\s\(Fragment\)//; s/\s[Ii]soform\s*.*$//;
  unless(m/protein\-protein/) { s/^(ORF|protein)\s*//i; }
  
  # no-no species names: Arabidopsis yeast  human
  # >> staphylococcal? = Staphylococcal nuclease ue, TAIR name
  # and: genome complete  pseudogene? = TAIR name
  # my $namedrops= $public_options{namedrops} || 'Arabidopsis|thaliana|yeast|complete sequence|complete|genome|pseudogene'; #plants;
  # s/\b($namedrops)[,\.\s]*//ig;
  s/\b(Arabidopsis|thaliana|yeast|human|Staphylococcal|complete sequence|complete|genome|pseudogene)[,\.\s]*//ig;
  s/\s*\((?:InterPro|TAIR):[\w\.-]+\)//ig; #  (TAIR:AT1G22000.1); (InterPro:IPR010678)
  if(s/paralog of //i) { $islike=1; }
  
  # horrible names:
  # hydrolases, acting on acid anhydrides, in phosphorus-containing anhydrides,ATP-dependent helicases,nucleic acid binding,ATP bi...
  # tRNA (guanine-N(1)-)-methyltransferase, metazoa , tRNA (guanine-N1-)-methyltransferase, eukaryotic == duplicated
  # mannose-1-phosphate guanylyltransferase (GDP)s,GDP-galactose:mannose-1-phosphate guanylyltransferases,GDP-galactose:glucose-1-
  # ATP-dependent peptidases,nucleotide binding,serine-type endopeptidases,DNA helicases,ATP binding,damaged DNA binding,nucleosid>
  # serine/threonine kinases,protein kinases,ATP binding,sugar binding,kinases,carbohydrate binding 
  # "The protein encoded by this gene was identified as a part of pollen proteome by mass spec analysis, It has weak LEA proteins, Encodes protein phosphatase 2A B'gamma subunit, Targeted to nucleus and cytosol"
  
  if(length($_) > 70) {
     my $nc= tr/[.,]/[.,]/; 
     if($nc > 1) { # keep only 1st two phrases.
     my $i= _min(index($_,',',0),index($_,'.',0));  
     while($i>0 and $i<30) { $i= _min(index($_,',',$i+1),index($_,'.',$i+1)); }
     $_= substr($_,0,$i) if($i>20);
     } 
  }       
#   if(length($_) > 70) {
#      my $nc= tr/,/,/; 
#      if($nc > 2) { # keep only 1st two phrases.
#      my $i= index($_,",",0);  $i= index($_,",",$i+1);
#      $_= substr($_,0,$i) if($i>20);
#      } 
#   }   

  ## should drop these geneid == name cases: At3g18210 for arabid genes, Os12g0628400, ...
  ## replace w/ special UNK name? Unchar prot gene id
  my $isid=0;
  
  # my $nameidpatt= $public_options{nameidpatt} || '(?:Os|At|AT)\d{1,2}g\d{3,}|GLEAN_\d+|G[A-Z]\d\d+'; #plants;
  my $nameidpatt= $NAME_IDPATT; ## $public_options{nameidpatt} || $NAME_IDPATT;
  if (/^($nameidpatt)/) { $isid=$1; $_= $NAME_UNK; } #? move to species.config
  elsif(s/($nameidpatt)//) { $isid=$1; } # $hasid=$1;?
    
  $isunk= ($pi < $MIN_NAMEIDENT or  m/^($NAME_NONE)/i)?1:$isunk; # Unknown|Uncharacterized|Hypothetical
  $isput= (!$isunk and ($pi < $MIN_CERTAIN or $isput))?1:0;
  
  ## ugh: TRIGALACTOSYLDIACYLGLYCEROL 1,
  my $nuc= tr/[A-Z]/[A-Z]/;  #  uc($_) eq $_
  if(m/[A-Z]\s+[A-Z0-9]/ and ($nuc > 19 or uc($_) eq $_ ) ) { $_= ucfirst(lc($_)); } # SHOUTING phrase ...  TRICHOME BIREFRINGENCE-LIKE 19
  
  s/\s*[Hh]omolo(gy|gue|g)\s+\d+//g; s/\b[Hh]omolo(gy|gue|g)\s*//g; # ? set $isput ? set -like? # add 'ortholog'
  s/\s*[Oo]rtholo(gy|gue|g)\s+\d+//g; s/\b[Oo]rtholo(gy|gue|g)\s*//g; 
  s/^(of|with)\s+//ig; # bad leading words
  # Add protein to the end when ending in:  'binding|domain|like|related'
  s/\b(binding|domain|like|related)\s(\W*)$/$1 protein $2/;  
  if( s/[,\s]*putative//g ) { $isput=1; } #s/putative, putative/putative/;

  # punctuation
  s/[\|]/:/g; s/#/n/g; 
  # s/_/ /g; # or leave uscores ??
  s/\s*$//; s/^\s$//; # lead/end spaces
  s/\s*[\/\.,;_-]*$//; # trailing punc
  s/[.] /, /g; # no sentences?
  s/^\W+//; # no leading crap
  # s/([,;:]){2,}/$1/g; #? double punc

  # plurals ?? anything ending in 's' ?
    
  # SEQ_FEAT.ProteinNameEndsInBracket:  Phosphoenolpyruvate carboxykinase [ATP] << change [] for ncbi
  if(/\]$/) {  s/\[/\(/g; s/\]/\)/g;}
  # unbalanced brackets; # add {} ? not used; <> ? not brackets
  if(/[\(\)\[\]]/) {
    my($nb,$ne,$d);
    $nb= tr/\[/\[/; $ne= tr/\]/\]/; $d=$nb - $ne;
    while($d>0) { s/\[//; $d--; } while($d<0) { s/\]//; $d++; }
    $nb= tr/\(/\(/; $ne= tr/\)/\)/; $d=$nb - $ne;
    while($d>0) { s/\(//; $d--; } while($d<0) { s/\)//; $d++; }
#  ##..  not working
#     my @bp=( "[","]", "(",")" ); 
#     while ( my($b,$e)= splice( @bp, 0, 2) ) { 
#     my $nb= tr/$b/$b/; my $ne= tr/$e/$e/; my $d= $nb - $ne; 
#     while($d>0) { s/\Q$b/./; $d--; } while($d<0) { s/\Q$e/./; $d++; }
#     } 
  }   
  
  $_=$NAME_UNK unless(/\w\w/); # 'Conserved protein' becomes blank
  if($isunk) { # regularize, but check/keep some additions
    if(/^($NAME_NONE)$/i or /^($NAME_NONE) protein/i) { $_= $NAME_UNK; }
    elsif(/^($NAME_NONE)\s*\w+/) { s/^($NAME_NONE)\s*/$NAME_UNK /;} # what?
    else { 
      # ?? save old,clean-name as Note in some/all cases.
      if($pi >= $MIN_IDLIKE) { $islike=1; } ##  for > 10-15% ident? or any?
      # Name=Mitogen-activated protein kinase kinase kinase 7-interacting protein 2
      # >> Mitogen-activated kinase kinase kinase 7-interacting protein 2-like protein
      else { $lowqualname=$_; $_= $NAME_UNK; }  # replace entirely? or not; NOT yet
    }  
  }
  if($islike) { unless(/family/ or /\blike/ or /protein kinase/) { s/\s+protein//; $_ .= "-like protein"; }  }
  s/protein protein/protein/ig; # other stutters?
    
  $_ .= ", putative"  if($isput and not $isunk);
  $_ .= " $locusname" if($isunk and $locusname);
  ## ncbi complains about this ^^ locusname addition; leave out at least w/ config.
  return ($_, $lowqualname);
}


=item processgene

  format structured model gene/(mRNA|ncRNA)/exon,CDS

  maybe add processfeature, for other things like transposons, ...
  processgene does this now if no mRNA.
  
  1719	638	repeat_region
    db_xref	baggins{}1471
    db_xref	FLYBASE:FBti0020395
    transposon	baggins{}1471

=cut

our %prot;

sub isPartial
{
  my( $mrna )= @_;
  
  my $partial= "";
  my $attr= $mrna->[8];  
  my ($prot) = $attr =~ m/protein=([^\s;]+)/; # use for alt-product=same/isoform marking
  if($prot) { # tinyaa not here.
    my $protc=0;
    $protc += 1 if($prot =~ /^M/);
    $protc += 2 if($prot =~ /\*$/);
    $partial = ($protc == 3) ? "complete" : ($protc == 1)? "partial3" : ($protc == 2) ? "partial5" : "partial";
  }
  unless($partial) {
    my ($qual) = $attr =~ m/quality=([^\s;]+)/;
    if($qual=~/Protein:\w*(complete|partial\w*)/) { $partial=$1; }
  }
  $partial.="53" if($partial eq "partial");
  return($partial);
}


# FIXME: add the new gene discrepancy checks here to bestgenes_update.pl
# .. see also genereplacex.pl where some of errors came from (i.e. replaced exons w/o adjusting mRNA span)

sub processallgene
{
  my( @geneall )= @_;
  return unless(@geneall);
  my @res=(0) x 3;
  my($gene,$mrna,$gc,$gb,$ge,$go, $gid, $ngerr, $nmrna)=(0) x 10;
  
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
      ## my ($issplit)= $mrna->[8] =~ /(split\w+)/;
      my ($qual)= $mrna->[8] =~ m/quality=([^\s;]+)/;
      my ($issplit)= $qual =~ /(split\w+)/;    # FIXME: Genbank tbl2asn needs new mRNA ID for these? same gene ID, not alt
      #my $mrnaispartial= isPartial($mrna);
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
          $er.="to:$go/$to," if($go and $to ne $go); 
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
    #my @res1= processgene( $generec, $geneother); 
    my @res1= processgene( $generec, $geneother, $gene, $geneispartial, $nmrna); 
    if($res1[1]>0) { $gene="done"; } # not working?
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
    warno "#ERROR: Missing mRNA/exons: @$generecIN \n";
    #putgene($generecIN, $geneother,"err=Missing-mrna-exon"); # flag="err=Missing-mrna-exon"
    return 0;
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
  my $gid= $id; $gid =~ s/t\d+$//;
  my $pid = reformatProteinID( $id); # (my $pid=$id)=~s/t(\d+)$/p$1/;
  my ($qual)= $attr =~ m/quality=([^\s;]+)/;
  ## my ($issplit)= $attr =~ /(split\w+)/; # not from entire attr, just from quality=
  my ($issplit)= $qual =~ /(split\w+)/;    #  FIXME: Genbank tbl2asn needs new mRNA ID for these? same gene ID, not alt
  ## quality=Class:Strong-expertchoice-splitgene:2,
  if(SPLITGENE_FIX) {
    if($issplit) {
      my $sid= 1 + ($splitgene{$id} || 0);
      $splitgene{$id}= $sid; # increment
      $splitgene{$pid}= $sid; # equivalence prot
      # if($splitgene{$id}) { my $partnum= $splitgene{$id}; return "gnl|$db|$id.$partnum";}
    }
  }
  
  my $partof=0; #? rely on this partof: annot?
  if($qual =~ /\-partof:(\w+)/) { $partof= $1; } # reclass as misc_feature fragment
  my $cdsexcept=0;
  if($qual =~ /Protein:curated/) { 
    $cdsexcept= "annotated by transcript or proteomic data";
    my $cid="";
    if($attr =~ m/\btrid:(\w+)/) { $cid= $1; }
    if(!$cid and my($rx)= $attr =~ m/rxnote=([^;\n]+)/){ ($cid)= $rx=~ m/cdnabest.(\w+)/; }
    unless($cid) { ($cid)= $attr =~ m/cdnabest.(\w+)/; }
    unless($cid =~ /[a-zA-Z]\d/) { $cid="evgasm"; } # FIXME
    $cdsexcept .=";inference=similar to RNA sequence, mRNA:TSA:$cid" if($cid);
    }  
  # also CDS needs /inference=    my $v= "similar to RNA sequence, EST:$DBXEST:$pi"; #% need some est-db config
  # my $v= "similar to RNA sequence, mRNA:$cdnaid"; #% need some est-db config

  
  my $alttr=0;
  if($attr =~ /isoform=(\w+)/) { # only if alttr > 1? # add oid == TSA inference:I100% ;
    $alttr=$1;
    my($oid)= $attr =~ m/oid=([^\s;]+)/;  # oid=vel4ma11:cacao3vel4sc1Loc938t2
    $attr .= ";trasm=$oid" if($oid and $oid !~ /^AUG/); # this is not always rna asm; need to remove AUG models
  } elsif($nmrna>1) {
    $alttr=1; $attr .= ";isoform=1"; # add it; 
  }

  my ($altprotsame, $altfirstid, $protpart, $tinyaa)= (0) x 10;
  # 201 ERROR: [SEQ_FEAT.DuplicateFeat] 
  ## see isPartial($mrna)
  my ($prot) = $attr =~ m/protein=([^\s;]+)/; # use for alt-product=same/isoform marking
  #? here or below# $prot="" if( $changelist{$id} =~ /MisMatchAA/i );
  $prot ||= "";
  $prot{$id}= $prot;
  my($aasize)= $attr =~ m/aaSize=(\d+)/; # trap tiny aa<40, recode misc_feature all? 30 total; 20 exceptions  BadProteinStart
  # retest prot complete here, qual is out of date for rx=3 updates
  if($prot) { # tinyaa not here.
    my $protc=0;
    $protc += 1 if($prot =~ /^M/);
    $protc += 2 if($prot =~ /\*$/);
    $protpart = ($protc == 3) ? "complete" : ($protc == 1)? "partial3" : ($protc == 2) ? "partial5" : "partial";
    $aasize= length($prot); 
    #? $aasize-- if($protc & 2 == 2); # -1 for * or not?
    my $nxx = $prot =~ tr/X/X/; # have 2 CDShasTooManyXs in other scaffs
    $tinyaa= "$nxx XXX" if($nxx * 2 >= $aasize);
    if($prot{$gid}) {
      $altfirstid= $gid."t1"; # fixme
      $altprotsame= ($prot{$gid} eq $prot)?1:0;
      # DO need to check all alts..  ** FIXME: out of order alts, t2 > t2 w/ same prot
      my $maxalt= _max($alttr,15);
      for ( my $t=1; $t <= $maxalt and not $altprotsame; $t++ ) {
        next if($t eq $alttr); my $aid=$gid."t$t"; 
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
  # if($alttr > 0  and not $altprotsame) { $newna .= " isoform $alttr";  } # should do only when prot differs from maintr
  if($alttr > 0) { $newna .= ($altprotsame) ? " isoform $altprotsame" : " isoform $alttr";  }  
    # ^^ problem for altprotsame isoform 1 << main t1 (usually) lacks isoform 1 suffix. ie it isnt alternate
  # if($alttr > 0) { $newna .= ($altprotsame) ? "" : ($alttr > 1)? " isoform $alttr" : "";  }  
  # if($alttr > 0) { $newna .= " isoform $alttr";   }  #$attr .= ";Note=alternatively spliced"; 
    # ERROR: valid [SEQ_FEAT.DuplicateFeat] Features have identical intervals, but labels differ
  #above# $attr .= ";Note=alternatively spliced" if($alttr>0);
  
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
  ($cdserr,$aabigger)= checkCDSexons($id, $mrna, $aasize, \@exon, \@CDS); # do above; add curated aa checks, partial annot
  if($aabigger > 0 and $cdsexcept) { #?? 
    $partial="partial"; #?? always, or cdsexcept 
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
  
  $attr .= ";".$genescore{$id} if($genescore{$id}); 

  # FIXME: need gene locus_id also (addgene); for partof, change gene=; for tinyaa, add new gene.
  # .. cancel this for Protein:curated == $cdsexcept
  if(($partof or $tinyaa) and not $cdsexcept) { # reclass as misc_feature fragment
    $type= "misc_feature";
    # FIXME: partof:notgoodgene, remove geneid, keep Note=partof
    if($partof) {
    my $partofid= reformat_myid($partof); #? or as MySRC:id
    my $partgene= ";gene=".$partof; $partgene=~s/t\d+$//;
    if( $changelist{$id} =~ /notpartof/i ) { $partgene=""; }
    $attr =~ s/^/Note=Gene fragment of $partof;Dbxref=$partofid;/;
    $attr =~ s/;gene=\w+/$partgene/;
    $gene= undef;
    } elsif($tinyaa) {
    ## need to putgene 1st if doing so. BUT gene needs mRNA or other such parts; needs CDS unless /pseudo or nc
    ## "misc_feature without a corresponding gene feature" ..
    ## FIXME: cant put CDS with huge gap: CDShasTooManyXs in other scaffs; change to /pseudo or chop off gap ?
    ##   $tinyaa= "$nxx XXX" if($nxx * 2 >= $aasize);
    ## only 2 cases in cacaogd:

    $attr =~ s/^/Note=Gene fragment aasize:$tinyaa;/;
    $attr =~ s/;gene=\w+//; # cant have locus_tag here
    $gene= undef;
    }
    my($attradd)= putCDSloc($type, $partial, @exon);  $nother++;
    putattr($type, $id, $attr); 
    @CDS=(); @exon= (); 
  }

  if($ADDGENE and $gene) { # FIXME. Defer this to after last alt-tr (and all parts of gene)
    # .. fix cases where later alt-tr is partial including gene span (and check all inside gene span)
    # .. maybe change filter_gff to collect all alttr into generecIN, but then
    # .. need 2-level array: mRNA1/exons/cdss need single package.
    $genepartial= $partial unless(defined $genepartial);
    my $gtype= $gene->[2];
    my $gattr= $gene->[8];
    if($gene->[6] eq ".") { my $go=$mrna->[6]; $gene->[6]= $go; warno "#WARN: fix missing strand: gene $gid:$go\n"; }  
    putCDSloc($gtype, $genepartial, $gene);  $ngene++;
    putattr($gtype, $gid, $gattr);
  }
  
  # [m]RNA, but may have retyped to @geneother...
  if(@exon) {
  putCDSloc($type, $partial, @exon);  $ntr++;
  putattr($type, $id, $attr); 
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

        # ** tbl2asn has some hidden controls of ERROR:overlaps another CDS with the same product name
        # ..~/Desktop/dspp-work/genomesoft/ncbi201107/api/sqnutil3.c : HasOverlapComment: overlap frameshift .. key words
        # RemoveCodingRegionsWithSuppressionWords() =     if (DoesStringContainPhrase (product, "ABC", TRUE, TRUE) .. transposon ..
        # OverlappingProductNameSimilar()
        # update .conf : attr_CDS cdsNote=
    # $cdsattr .= ";cdsNote=CDS overlap for alternateUTRof:$altfirstid" if($altprotsame); # << in mrna attr but NOT here in cdsattr
    #see isoform# 
    #above# $cdsattr .= ";cdsNote=alternatively spliced" if($alttr>0);
    
## handle some attr for this for valid CDS <> protein
## exception: 	"unclassified transcription discrepancy" 
## need regex for Dbxref > mRNA, eg. ncbi [XN]M_ ids vs [XN]P_ ids

    #above# my ($cdserr)= checkCDSexons($id, $mrna, $aasize, \@exon, \@CDS); # do above; add curated aa checks, partial annot
    if($PEPfromGFF) {  ##  and not $altprotsame
      ## new changelist MisMatchAA = skip my pep, let ncbi translate. see putCDS
      $prot="" if( $changelist{$id} =~ /MisMatchAA/i );
      my $pepok = putprot( $pid, $prot);
      $cdserr .= "noaa" unless($pepok);
    }
    
#     if($partial eq "partial5" or $partial eq "partial") { # do this in putCDS
#       $cdsattr .=";codon_start=".$phase+1;
#     }
    my($attradd)= putCDSloc($type, $partial, @CDS);  $ncds++;
    putattr($type, $id, $cdsattr.$attradd);
    warno "#ERROR: $id CDSexon, $issplit exon ne mRNA: $cdserr\n" if($cdserr); # can have cdsold,new lines appended
  }
 
  foreach my $x (@geneother) {
    my $type= $x->[2];
    putloc($type, $x);  $nother++;
    putattr($type, 0, $x->[8]); 
  }
  
  return($nref,$ngene,$ntr,$ncds,$nother);
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
  
  # my $cdsaa= int($cdsl/3); # includes stop codon
  my $errcds="";
  $aabigger=  $aasize*3 - $cdsl; # smaller? *stop -1?
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
    
      # FIXME here? collect all alttr of same gene before processgene().
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
#WARN: Thecc1EG047060t2 is wrong? gene=Thecc1EG028355; not Thecc1EG047060 .. < this one right
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
        if( $changelist{$tid} =~ /^drop/i ) { $ndrop++; @generec = (); next; } # can include gene or alttr drop
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
        if( $changelist{$tid} =~ /^strand\s([+-])/i ) { my $onew=$1;
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
  unless( $PEPfromGFF ) { # and -f $gbsubpep
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
