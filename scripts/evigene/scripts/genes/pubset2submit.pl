#!/usr/bin/env perl
# evgs/genes/pubset2submit.pl 

=item pubset2submit.pl

 revision of parts of  evgmrna2tsa2.pl
 roughly same as MAIN_submitonly() of evgmrna2tsa2
 from genes/trclass2pubset.pl MAIN_pubsetonly of evgmrna2tsa2
   
=item usage (maybe)

  set pt=pig4321ew; 
  $evigene/scripts/genes/pubset2submit.pl -debug -log -mrna publicset/$pt.mrna_pub.fa
 
  equivalent to this, which needs both -mrna -class for -onlysubmit
  $evigene/scripts/evgmrna2tsa2.pl -onlysubmit -log -debug -NCPU 1 -class pig4321ew.trclass -mrna publicset/pig4321ew.mrna

  -- opts -noruntbl2asn -NCPU 1  [more cpu for runtbl2asn, to split job]
  -- look for publicset/ files in current dir
  -- need sra.csv = pig4321ew.sra_result.csv = "SRAtable|datatable=s", \$sratable, 
   
=item FIXME tsancbigapck.pl

 damn ever changing ncbi gap policies;
 need to add changable gap policy scanner, for mrna/cds to fsa/tbl
 ncbi tsa submit now has different endgap policy 
   not same as their current tbl2asn checker
   and it takes hours and hours to get thru the ncbi web tsa checker to find this out..
   and days of rechecking to figure out computation of that text-defined gap policy 

 /bio/bio-grid/verts/sra2genes/pig18evgsub
 tsancbigapck.pl

=cut

use FindBin;
use lib ("$FindBin::Bin","$FindBin::Bin/.."); # assume evigene/scripts/ layout: main, prot/, rnaseq/, ...
use strict;
use warnings;
use Getopt::Long;
use cdna_evigenesub; 
use evigene_pubsets; # now has some of below subs

# replace common settings subs below w/ this
# use evigene_settings;

use constant VERSION => '2018.07.10'; # '2018.06.25'; #  from evgmrna2tsa2.pl & altbest2pubset

our $EVIGENES= "$FindBin::Bin/.."; # ="$FindBin::Bin";  # WRONG place ; this is in evgs/genes/ subdir ..
our $DEBUG= $ENV{debug}|| 0;
our $EGAPP='pub2submit';   
our $EGLOG='m2t';
our $dryrun=0; ## $DRYRUN ?

our $IDPREFIX= $ENV{idprefix} || 'EVGm'; 
our $FAHDRisValid=0; # export from pubset.pm for annots from mrna hdr
our ($pubid_format,$altid_format,$GBPROID);
our $DATE;  
our $RNATYPES; # our $RNATYPES='mRNA|ncRNA'; # in evigene_pubsets.pm; transcribed_RNA mrnatypes exontypes 
our ($KEEPDROPX,%keepdropxtab); # evigene_pubsets.pm

# main opt:
my $AAMIN_NOCLASS=$ENV{aaminnoclass}||60; # asmrna_altreclass -noclasscut=$AAMIN_NOCLASS; drop noclass (needs user opt), but rescued by aaref
my $NOALTDROPS=0; # turn on if have keepdrop table?
our($ORGANISM,$sraids,$BioProject)=("Noname","SRR000000","");  ## SRR 346404
our @sraids;
# opt for tbl2asn:
my $NCPU= 1; 
my $DOtbl2asn=0; 
my $DOpep4asn=0; # default was on, maybe turn off?
our $TSADESC="-w evgr_tsamethods.cmt -Y evgr_tsadesc.cmt -t evgr_tsasubmit.sbt"; # FIXME: where from?
our $TSAOPTS="-a r10k -l paired-ends -Vt -Mt -XE"; # for tbl2asn, read from config.;  
    ## -Vb = gen genbank gbf, not needed, option; -Vtb == genbank + tsa
our $DBRECODE= $ENV{dbrecode};


my ($trclass,$output,$tblfile,$logfile,$keepdropin);  
my ($trpath,$trname, $sradatah, $sratable, %settings);
my ($oname, @genes, @seqset, $genenames);
my %didpubid=(); #only for make_tblfsaset ??
# my ($vecscreenf,$trclass,$genenames,$cdnaseq,$output,$logfile,$tblfile,$runsteps); ## ,$dryrun above

our @evgdirs = qw(okayset dropset inputset trimset tmpfiles erasefiles publicset submitset);
our (@okayset,@dropset,@inputset,@tmpfiles,@erasefiles,@publicset,@submitset); # tidyup file sets
my $tidyup=1;


my $optok= GetOptions(
  # "trclass|input=s", \$trclass, # many gff inputs : use perl input pipe while(<>) ?
  "mrna|cdna|sequences=s", \@seqset, # many inputs? 
  # "oname|outbase|output=s", \$oname,  # FIXME: idpre option  overwritten by spppref
  "SRAtable|datatable=s", \$sratable, # see sra2genes, look for
  "idprefix=s", \$IDPREFIX,  
  # "tagidnew=s", \$NEWIDpre,   
  # "source=s", \$SRC,   
  # "names=s", \$genenames,   
  # "keepoldids|preserveOldIds=s", \$preserveOldIds,  # allow -idpre=XXX -keepold w/o 2nd param
  # "keepdropin=s", \$keepdropin,   # other opt key? reclassin?
  "DBRECODE=s", \$DBRECODE,  
  "DATE=s", \$DATE,  
  "NCPU=i", \$NCPU,## "MAXMEM=i", \$MAXMEM,  
  "logfile:s", \$logfile,
  "runtbl2asn!", \$DOtbl2asn, 
  "pep4asn!", \$DOpep4asn,   # now off by default; 
  "dryrun|n!", \$dryrun, 
  "tidyup!", \$tidyup, # default on ?
  "debug!", \$DEBUG, 
  );


die "EvidentialGene pubset2submit -mrna publicset/myspp.mrna_pub.fa
  makes fileset for submission of gene set to NCBI/EBI/DDBJ Transcript Archive 
  opts: -[no]runtbl2asn -sratable sra_result.csv -pep4asn -log -debug -NCPU=$NCPU
  version ", VERSION, "\n"
  unless($optok and @seqset); # (@genes or @seqset) 

# $NEWIDpre="" unless($NEWIDpre=~/^[A-Za-z]/); # allow same old id?
# my $IDPREOK=$IDPREFIX;
unless($DATE){ $DATE=`date '+%Y%m%d'`; chomp($DATE); }
my $DEFAULTidpre= 'NonameEVm';  ## this seems to fix -idprefix ignored bug

my %DEFAULT_SETTINGS= ( 
  IDPREFIX=>$DEFAULTidpre, 
  DATE=>$DATE, 
  organism => $ORGANISM, 
  sraids => $sraids, 
  BioProject => $BioProject, #fixme bioproj lost
  assemblers => 'Velvet/Oases; idba_trans; SOAPDenovoTrans; Trinity;',
  trclass => '', mrna => '', genenames=>'', 
  # genome => 'genome/chromosome_assembly.fasta.gz',
  ); 


our %DBXREF_OK=(); # ncbi/db_xref; should get from evigene_config();
{
	my $dbxref_ok="AFTOL AntWeb APHIDBASE ApiDB ASAP ATCC Axeldb BEETLEBASE BGD Bold CABRI CCAP CDD dbEST
dbProbe dbSNP dbSTS dictyBase EcoGene ENSEMBL EPD ESTLIB FANTOM_DB FBOL FLYBASE GABI GDB
GenBank GeneDB GeneID GI GO GOA Greengenes GRIN H-InvDB HGNC HMP HOMD HSSP IKMC IMGT
InterPro IntrepidBio IRD ISFinder JCM JGIDB LocusID MaizeGDB MGI MIM
miRBase MycoBank NBRC NextDB niaEST NMPDR NRESTdb OrthoMCL Osa1 Pathema PBmice PFAM PGN
Phytozome PIR PomBase PSEUDO PseudoCap RAP-DB RATMAP RFAM RGD RiceGenes RZPD SEED SGD SGN
SK-FST SoyBase SubtiList TAIR taxon TIGRFAM UNILIB UniProtKB/Swiss-Prot
Swiss-Prot SwissProt UniProtKB/TrEMBL TrEMBL UNITE VBASE2 ViPR WorfDB WormBase ZFIN";
	my @pg= split/[,\s]+/, $dbxref_ok; 
	map{ my($k,$v)=split/[=:]/,$_; $v=1 unless(defined $v); $DBXREF_OK{$k}=$v; } @pg; 
}       

our %DBXREF_RECODE=();
if($DBRECODE) {
  my @pg= split/[,\s]+/, $DBRECODE; 
  map{ my($v,$k)=split/[=:]/,$_; my @k=split/[|]/,$k; foreach my $k1 (@k) { $DBXREF_RECODE{$k1}=$v; } } @pg; 
}

my $SKIPTRIMSET=0; # unused here
my $skipTSAparts=($SKIPTRIMSET and not $DOtbl2asn)?1:0;
 
# data globals
# my($mainindex,$ntr)= (0) x 9;
# my(%annotes, %main, %mainsize,%alt,%oclass,%altsize,%altdrops,%didaltdrops,%balt,%drop);  # globals of trclass2maintab
#unused here: my(%aaqual, %piad, %newmain, %altdefer,%maindrops);

# my($nkeepdrop, $keepdroph, )= (0);
# my($nmapqual, $mapqualh, $alntabh)= (0); #? globals from map align.tab, if exists

my $mrnaseq= (@seqset>0)? $seqset[0] : ""; # many?? all input .mrna .cds and/or .aa

openloggit($logfile,$mrnaseq); # is this good log: publicset/pig4321ew.mrna_pub.pub2submit.log
loggit(1, "EvidentialGene pubset2submit VERSION",VERSION);

my $APPtbl2asn  =  findapp("tbl2asn") if($DOtbl2asn); # add this, needs configs + data doc files (.cmt, .sbt,..)
# my $APPtraa2cds= findevigeneapp("prot/traa2cds.pl");
# my $APPtrimvec= findevigeneapp("rnaseq/asmrna_trimvec.pl"); # NOT run here, see get_trimvecset()
# my $APPaltreclass= findevigeneapp("rnaseq/asmrna_altreclass.pl",1); # upd 201405 ;  NODIE

## FAIL at this point if any apps missing?
#-------------------------------------


MAIN();

sub MAIN
{
  ## default output to publicset/oname
  ## maybe default input seqs from okayset/ ; see mrna2tsa:get_evgtrset()
  my($upokgenes, $upokids, $upstatus)= (0) x 9; 
  my $mrnaseq= (@seqset>0) ? $seqset[0] : ""; # many?? all input .mrna .cds and/or .aa

  loggit(0, "BEGIN $0  input=",$mrnaseq,"date=",$DATE);
  %settings= %DEFAULT_SETTINGS;
  do_settings("restore",$mrnaseq);  
	
  ($mrnaseq,$trpath,$trname)= get_evgtrset("ignoretrclass",$mrnaseq,"publicset"); # subs.pm; cdnaseq => mrnaseq here
  loggit(0, "get_evgtrset=",$mrnaseq,$trpath,$trname);  
  loggit(LOG_DIE, "Missing -mrna",$mrnaseq) unless($mrnaseq and -s $mrnaseq);
  
  unless($oname){ # oname == trname ?
    $oname= $trname;
    # $oname=$mrnaseq; 
    # unless( $oname =~ s/\.mrna_pub\.fa|\.mrna// ) { $oname=~s/\.\w+$//; }
  }

	# sracvs == "SRAtable|datatable=s", \$sratable, 
	unless($sratable and -s $sratable) {
    # my($srat)= getFileset('.','sra_result.csv$|sra.csv$');  
    my $csv="$trpath/$trname.sra_result.csv";    
	  $csv="$trpath/sra_result.csv" unless(-f $csv);  
	  $sratable=$csv if(-f $csv);
	}

	#o# ($sradatah)= get_srainfo($sratable,$trpath,$oname); # trname ?
	($sradatah)= get_srainfo($sratable); # ,$sraids ?
	
	## these pulled from get_srainfo ..
	my($nupinfo,$tsamethf,$tsadescf,$tsasubf)= tsa_infotemplates($trpath, $oname, $sradatah);
	loggit(0,"info updated $nupinfo TSADESC=",$TSADESC);
	unless($genenames) { #? put in get_evgtrset
		my $gnt="$trpath/$trname.names";  $genenames=$gnt if(-s $gnt); 
		loggit(LOG_WARN, "Missing product -names",$gnt) unless($genenames);
	}
	

  ## publicset/pig4321ew.mrna_pub.fa  >> trname & oname= pig4321ew
  my($pubids,$annotab,$maltab);
  my($pubd)= getFileset("publicset",'pubids|ann.txt|mainalt.tab');  
  if(@$pubd){
    ($pubids)= grep/\.pubids$/, @$pubd; # grep trname ?
    ($annotab)= grep/\.ann.txt$/, @$pubd;
    ($maltab)= grep/\.mainalt.tab$/, @$pubd;
  } else {
    ($pubids,$annotab,$maltab)= map{ "$oname.$_" } qw(pubids ann.txt mainalt.tab);
  }
  
  ## change this mrnaseq to array= grep 'mrna|cdna|xxx' @seqset
  ($upokgenes)= makeSubmitset($pubids, $annotab, $oname, $mrnaseq); 
  $upstatus++ if($upokgenes>0);
 
  do_settings("log|save",$mrnaseq,); 

  if( $tidyup and $upstatus > 0 ) {  # move to sub package
    tidyupFileset("publicset",@publicset);  
  	tidyupFileset("submitset",@submitset);  #? after or before tsa_tbl2asn? need path in fileset
    tidyupFileset("tmpfiles",@tmpfiles);  
    my @rmlist;
    foreach my $fn (@erasefiles) { if(-f $fn) { unlink($fn); push @rmlist,$fn; } }  
    if(@rmlist) { my $nrm=@rmlist; my $rml=join" ",@rmlist[0..4]; loggit(0,"tidyup erase: n=$nrm, $rml .."); } 
  }

  loggit(0, "DONE at date=",`date`);
  loggit(0, "======================================"); # log may append runs
}
  
# ---------------------------



sub makeSubmitset {
  my($pubids, $annotab, $outname, $cdnaseq)= @_;

  my $skiptrrun=0; my $upstatus=0;
  our @publicset; #global in pubsets.pm
	my $pubdir="publicset";
	
	## dangit NOT FAHDRisValid for many mrna seq entries, when .pubids/ann.txt valid
	## set invalid if have both pubids and ann.txt?
  $FAHDRisValid=1; # export from pubset.pm for annots from mrna hdr

  loggit(0,"makeSubmitset($pubids,$outname)"); # LOG_NOTE,LOG_DEBUG // if($DEBUG);

  my(@mrna,@aa,@cds);
  ## change this cdnaseq/mrnaseq to array= grep 'mrna|cdna|xxx' @seqset
  if($cdnaseq and -f $cdnaseq) {  @mrna=($cdnaseq); } 
  else {
    @mrna= grep /\.(mrna|ncrna|rna|cdna|tr)$/, @seqset;
    @aa  = grep /\.(aa|pep)$/, @seqset; # if empty but @mrna try s/.mrna/.aa/ ?
    @cds = grep /\.(cds)$/, @seqset;
    if(@mrna) { $cdnaseq= $mrna[0]; } 
  }
  
  if(@mrna and not @aa) { @aa= map{ (my $f=$_) =~ s/mrna/aa/; $f; } @mrna; }
  if(@mrna and not @cds) { @cds= map{ (my $f=$_) =~ s/mrna/cds/; $f; } @mrna; }
  my $aaseq=$aa[0]||"";
  
  if($genenames and -f $genenames) { }
  elsif( -f "$outname.names") { $genenames= "$outname.names"; } #? need -option
  else { my $pn= makename($pubids,'.names'); $genenames=$pn if(-f $pn); }
  
	my($npubid, $pubidh)= read_pubids($pubids, $cdnaseq); 
  set_setting('mrna',$cdnaseq) if($npubid);

  ## this *should* replace annotation names for tbl2asn outputs .. want to fix problem cases
  if($genenames and -f $genenames) {
    my($nagot,$nain)= parse_genenames($genenames); 
    set_setting('genenames',$genenames) if($nagot);
    }
  
  # my %annothash= read_annotab($annotab); #** fill in annothash from pubidh if no table
  my($nann,$annothash)= read_annotab2($annotab);
  $FAHDRisValid=2 if($npubid and $nann > 0); ## upd1809; FLAG 2 == sensible merge only missing vals
  
  
  my($outfa, $tblout, $ntrout, $ncdsout, $pepout)=("") x 9;
	($outfa, $tblout, $ntrout, $ncdsout, $pepout) 
			= make_tblfsaset($cdnaseq,$annothash,$aaseq,$skiptrrun);  

	map{ push @submitset, $_ if($_); } ( $outfa, $tblout, $pepout);  #  added pep, may be missing
	## fixme: @submitset needs sbt,cmt of  TSADESC=>$TSADESC,
	$upstatus++ if($ntrout>0); # 0 if files already made
	loggit(0,"submitset: ntr=$ntrout, ncds=$ncdsout in $outfa, $tblout"); 

	if($DOtbl2asn) {
    my($npartgot,$spldir,$sqnoutlist)= tsa_tbl2asn($outfa,$tblout,$ORGANISM,$sraids);
    loggit(0,"tsa_tbl2asn nparts=$npartgot, submitset=$sqnoutlist"); #?? $spldir  
    
    # tbl2asn output in publicset/tsa > move to submitset/ with input files
    # pig4321ew.mrna_pub_tsasubmit
    # pig4321ew.mrna_pub.pub2submit.info
    # pig4321ew.mrna_pub.pub2submit.log
    # pig4321ew.mrna_pub.sumval
    
	} ## FIXME: submitset tbl2asn: name.val,.sqn,.discrep,.fixedproducts errorsummary.val
  
  my $nout= 1; ## $npm || $npa || $npc || 1;
  #fix: loggit(0,"submitset: ",$pubmrna,$minfo,$pubaa,$ainfo,$pubcds,$cinfo,$pubgff,$ginfo,$annotab); 
  return($nout);
}

my($AAGAPENDS, %AAGAPENDS); # global for putCDS..

sub make_tblfsaset
{
  my($mrnaseq,$annothash,$aaseq,$skiprun)=@_;
  my($notr,$nocds)=(0,0);
  my $outfa = ($output)  ? $output  : makename($mrnaseq,".fsa"); # was .fna ; tbl2asn wants fsa
  my $tblout= ($tblfile) ? $tblfile : makename($mrnaseq,".tbl");
  
	my $tsadir="submitset";
  if(not -f $tblout and $tsadir and -d $tsadir) {
  	 my($tsad,@files)= getFileset($tsadir,'tbl|fsa');   
   	 my($ft)= grep /\.tbl/, @$tsad; $tblout=$ft if($ft);
     ($ft)= grep /\.fsa/, @$tsad; $outfa=$ft if($ft);
  }

  # fixme: look in publicset/ for annot, submitset/ for outfa,tblout
  return($outfa,$tblout,$notr,$nocds) if( -s $outfa or $skiprun); # or dryrun ..

  my ($inh,$outh,$tblh,$annoth,$hd,$oid,$fa,$ok, $inaah, $peph);
  #FIXME: 1412: add .pep output, from mrnaseq>aaseq, using same protids in tbl
  $aaseq= makename($mrnaseq,".aa") unless($aaseq); 
  my $pepout= ($DOpep4asn and -s $aaseq) ? makename($tblout,".pep"):""; 

	($ok,$inh)= openRead($mrnaseq);
  $ok= open($outh,'>',$outfa) if($ok);
  $ok= open($tblh,'>',$tblout) if($ok);
  unless($ok) { loggit(1,"ERR: make_tblfsaset $mrnaseq TO $outfa,$tblout"); return; }

  # patch for nasty bug from uvcut: partial3 'N' > aa.X at end : trim last codon AND aa.X 
  my($hasendgap, %aaendgap)=(0);
  # in putCDS:     if($AAGAPENDS and my $agaps= $AAGAPENDS->{$oid}) { .. trim mrna }
  if(-s $aaseq) { 
    my($ok,$inh)= openRead($aaseq); 
    my ($id,$aa,$isuvcut)=(0) x 9;
    if($ok) { while(<$inh>) {
      if(/^>(\S+)/){ my $d=1; 
        if($id and $aa and $aa=~m/X$/) { $aaendgap{$id}=1; $hasendgap++; }
        $id=$d; $aa=""; # $isuvcut=(m/uvcut=([^;\s]+)/)?$1:0; 
      } else { chomp; $aa.=$_; }
    }
    if($id and $aa and $aa=~m/X$/) { $aaendgap{$id}=1; }
    close($inh); 
 	  loggit(0, "aaendgap n= $hasendgap");
   }
  }
  #?? bad glob: $AAGAPENDS= ($hasendgap)? \%aaendgap : undef;
  if($hasendgap) { %AAGAPENDS= %aaendgap; $AAGAPENDS=$hasendgap; }
  else {  %AAGAPENDS=(); $AAGAPENDS=0; }
     
  my %protids=();
  my($itr,$otr,$ocds,$oerr,$protid)= (0) x 10;
  my $pubidin= 0; # FIXME: is this always oid here?
  while(<$inh>) { 
    if(/^>(\S+)/) { my $d=$1; # is this oid or pubid : EITHER **
			if($fa and $oid and ($pubids{$oid})) #  or not $skipdropseqs
      {
				# my $tblinfo= parse_evgheader($oid,$hd,length($fa)); 
				#package global: our $FAHDRisValid
     		#was.bug:  $FAHDRisValid for pubidin
     		my $annorec= annotab2tblinfo( $oid, $annothash->{$oid}, $hd, $pubidin); ## = $annothash{$oid}||{}; 
     	  ($otr,$ocds,$oerr,$protid)= putTblFsa($outh,$tblh,$oid,$hd,$fa,$annorec,$itr); 

     	  if($oerr ne "OK") { }
     	  $protids{$oid}= $protid if($protid);
	      $notr+= $otr; $nocds+= $ocds;
      }
      $oid=$d; $hd=$_; chomp($hd); 
      $fa=""; $itr++;
      }
    elsif(/^\w/) { chomp; $fa.=$_; }
  } 
  
	if($fa and $oid and ($pubids{$oid})) #  or not $skipdropseqs
	{
    # my $tblinfo= parse_evgheader($oid,$hd,length($fa)); 
 		my $annorec= annotab2tblinfo( $oid, $annothash->{$oid},$hd,$pubidin); ## = $annothash{$oid}||{}; 
		($otr,$ocds,$oerr,$protid)= putTblFsa($outh,$tblh,$oid,$hd,$fa,$annorec,$itr); 
    $protids{$oid}= $protid if($protid);
	  $notr+= $otr; $nocds+= $ocds;
	}
  close($inh); close($outh);
  
  ## do this at end? collect oids, protids from putTblFsa()
  ## pepout should be option, at least into tbl2asn, 
  ## using its orf-caller may be less hassle
  if($pepout and scalar(%protids)) { 
    my($ok,$inh)= openRead($aaseq);  
    my $cutendgap=0; my %didp; # block dupids
    $ok= open($outh,'>',$pepout) if($ok); 
    unless($ok) { loggit(1,"ERR: make_tblfsaset $aaseq TO $pepout"); }
    else { while(<$inh>) {
      if(/^>(\S+)/){ my $d=$1; $ok=0; $cutendgap=0;
        if(my $p= $protids{$d}){ s/>$d.*$/>$p/; $ok=($didp{$p}++)?0:1; 
          if($ok and $hasendgap) { $cutendgap= ($AAGAPENDS{$d}||0) == 2; } } 
      }
      else { if($cutendgap) { s/X$//; } } # FIXME only at end of aa seq
      print $outh $_ if $ok;    
    }
    close($inh); close($outh);
    }
  }
  
  return($outfa,$tblout,$notr,$nocds,$pepout); # pepout ??
}

sub tblMaploc {
  my($tblinfo)= @_;
  my ($mnote, $mapl, $maln, $midn, $mlen, $mxon)=( 0) x 9;
  if($mapl= $tblinfo->{'maploc'}) {
    my $mapq= $tblinfo->{'mapqual'} || "";
    my $mapp="";
    $maln= $mapq =~ m/(\d+)a,/ ? $1 : 0;
    $mapp .= ($maln > 9 and $maln < 90)?"$maln pct aligned,":"";
    $midn= $mapq=~m/,(\d+)i/ ? $1 : 0; # skip, or $1-identity, 
    $mlen= $mapq=~m/,(\d+)l/ ? $1 : 0; # skip length
    $mxon= $mapq=~m/,(\d+)x/ ? $1 : 0; 
    $mapp .= "$mxon exons," if($mxon);
    $mapl=~s/:[fr\+\.\-]$//; # drop strand always?
    $mnote= ($maln < 10)? "unlocated on chromosome assembly" : "located at $mapl, $mapp"; 
    }
  return($mnote);
  #return ($mnote, $mapl, $maln, $midn, $mlen, $mxon);
}

#already.def: use constant FOR_NOTE => 2;
sub putTblFsa
{
  my ($outh, $tblh, $oid, $hdr, $faseq, $tblinfo, $itr)=@_; 

	# FIXME: annorec becomes tblinfo, not from faseq.hdr parsing here
  # my $tblinfo= parse_evgheader($oid,$hdr,length($faseq));
	# ADD global   $protids{$oid}= $protid; for .pep output
  # OPTION: isoform attrs, use only pubid t[1..n] syntax? 
  #   need also protein == and not= info (isoform A,B,C) vs same isoform/diff mrna
  #  pubgenes/kfish2rae5h.puballnr.isoforms .. has needed attr
  ## fixme: ann.txt cdsoff has dang :+/- screws up tbl and mrna<>tbl .. fix in ann.txt
  
  my($ntrout,$ncdsout,$ncrnaout,$nfatrim)=(0) x 10;
  my $pubid= $tblinfo->{'pubid'} || $oid; ## or what? err?
  my $falen= length($faseq);  
  my $oidin= $oid; if($oid eq $pubid) { $oid=$tblinfo->{'oid'} || $oidin; }  # which? dont want oid == pubid

  my($cdsoff,$gname,$dbxref,$aaqual,$trlen,$protid,$lotag,$namepct, $cddname)= 
      map{ $_ || "" }
      @{$tblinfo}{qw(cdsoff name dbxref aaqual trlen protid locustag namepct cdd)};
  my $Selcstop= $tblinfo->{'Selcstop'} ||"";
  my $seqtype= $tblinfo->{'seqtype'} || "mRNA"; # FIXME ncRNA
 
  $namepct =~ s/,.*//; $namepct =~ s/%//; # my($nap)=$namepct=~m/(\d+)/;
  my $gnameref= $tblinfo->{nameref} || $genedbxref{$oid} || $dbxref || ""; #?? $pubid  # also require dbx DBXREF_OK limited set
  # missing nameref, ann.txt dbxref contains now? some are XM_ not XP_ ?
  if( $gnameref =~ /,/ ){
    my($ndx)= split",",$gnameref;  $ndx=~s,/.*,,;  $gnameref=$ndx;
  }
  
  $cdsoff =~ s/:[+.-]$//; # drop strand if there.
  my($cdsb,$cdse)= split/[-]/,$cdsoff; # dang, :[+-.] strand
  if($cdse and $cdsb > $cdse) { ($cdsb,$cdse)=($cdse,$cdsb); $cdsoff="$cdsb-$cdse"; } # FIX in putCDSloc or here?
  ## dang: does this b>e mean mRNA seq is not reversed?
  
  my($cdsOK)= ($cdse and $cdse > 0 and $cdsb > 0)?1:0; # buggy data, eg cdsoff == 0
  my($ncRNAok)= ($seqtype ne "mRNA" and not $cdsOK)?1:0; #?
  my $cdsphase=0; # unused here?
  my($aafull)= $aaqual =~ m/(complete|partial\w*)/? $1 : "partial"; #? missing aaqual??
  my $aastop=  ($aafull =~ /partial3|partial$/)? 0 : 1;
  my $aastart= ($aafull =~ /partial5|partial$/)? 0 : 1;

  map{ $_='' if($_ eq "na"); } ($protid,$lotag,$gname,$cddname,$dbxref,$aaqual);  

  $aafull ||= $aaqual;
  my($cdnap,$cdsym);# prodname: $cdsym also ? add to comments?
  my $hascddname= ($cddname and $cddname ne $gname)?1:0;
  ($gname,$namepct)= productname($gname,$aaqual,$namepct, FOR_NCBI); #? add FOR_NCBI flag?
  ($cddname,$cdnap,$cdsym)= productname($cddname,"100,99%,complete",0, FOR_NCBI|FOR_NOTE) if($hascddname); #? add FOR_NCBI flag?
  
  if($didpubid{$pubid}++) { # catch dup seq ids !
    loggit(1,"ERR: dup seqid $pubid/$oid");
    return($ntrout,$ncdsout,"ERR",0);
  }
  #* FIXME: ncRNA
  unless($cdsOK or $ncRNAok) { # or $cdsoff =~ m/\d+\-/  # err, missing cds span for tbl; print fsa seq anyway?
    loggit(1,"ERR: missing cds-offset $pubid/$oid:",$cdsoff,$aaqual,$gname);
  }

=item fix shift partial cds at +1/+2 of mrna-ends, or mrna-gaps

    FIXME.1510: NCBI requires mrna-cds partial to abut ends of mrna,unless gaps, 
    ie. ofs=mrnastart+2,3 wrong and ofs=mrnaend-1,2 wrong; for start shift, also shift cdsphase/codonstart
    -- problem with this, of course, when mrna-*end* is +2 from partial3-cds end,
      cannot extend cds to mrna end w/o calling 2-base codon, adding to protein :((
    -- answer: trim mrna-end to cds-end (only for +2 case or +1 also?)
    
   double darn ncbi again, trim mrna end when +2 or +1 of cds-partial3 end
   .. need putCDSloc() logic before print faseq. 
   .. mrna length doesnt appear in outputs here to ncbi asn except in print faseq.

=item also fix inner cds +1,2 of NNN gaps
      
      inner gap cds-off by 1,2            ofs     trlen obeg         v'5  oend     v'3
  Dapma6tiEVm002532t2	  814,67%,partial5;808-3252	3646	808	NNNNNNNCA.GGA	3252	TGA.ACACTTTTG
  Dapma6tiEVm003622t3	  520,71%,partial5;344-1906	2183	344	NNNNNNNAT.CAT	1906	TGA.AATAAAAAT
  Dapma6tiEVm015968t20	130,98%,partial;7-396	     397	  7    AANNNA.AAA  396	AAA.
      -- extend inner partial cds'5 -1,-2 if Not N, but precedes NNN by 1,2
      
      3'end-N-bug                         ofs     trlen obeg         v'5  oend     v'3
  Dapma6tiEVm009340t53	264,71%,partial3;316-1107	1108	316	GTATTTTTA.ATG	1107	NAN.A
  Dapma6tiEVm001867t281	235,99%,partial3;1-705     706    1	         .ATG  705	AAN.C
  Dapma6tiEVm009340t60	203,92%,partial3;49-657    658   49	AAGAAGATA.ATG  657	NGN.G
      -- trim -1 cds'3 if == 1-N only
=cut
  
  $faseq= uc($faseq); # ugh, some fix gaps are 'nnn'.
  
  my($cdsline,$fatrim,$cdsoff2,$cdsphase2,$mrna2dif,$mrna2)= (0) x 9;
  if($cdsOK) {
    ($cdsline,$fatrim,$cdsoff2,$cdsphase2,$mrna2dif,$mrna2)= 
          putCDSloc($cdsoff,$aafull,$cdsphase,$falen,$faseq,$oidin); 
    if($cdsoff2 and $cdsoff2 ne $cdsoff) { $cdsoff= $cdsoff2; ($cdsb,$cdse)= split/[-]/,$cdsoff; } ## always reset?
    $cdsOK=0; if( $cdsoff =~ m/(\d+)\-(\d+)/ ) { my($cb,$ce)=($1,$2); $cdsOK=($ce>0 and $cb>0)?1:0; }
    
    $cdsphase= $cdsphase2;
      ## mrna2dif problem: faseq change inner5'part, leading 1,2 AGCT before NNN, change to NN also,
    if($mrna2dif) { $faseq=$mrna2; $falen= length($faseq); $nfatrim++; } 
    if($fatrim>0) {
      $nfatrim++; # *should* log any changed mrna/cds for later check, as aa.qual table? id/alen/gap/aqual/tlen/cdsoff/oid
      $faseq= substr($faseq, 0, $falen - $fatrim);
      $falen= length($faseq);  
    }
    ## if($trlen != $falen) { } #?? warn? use falen
  }
  # loggit(0,"CDSerr: $pubid") unless($cdsOK); #see above
  
  $faseq =~ s/(.{60})/$1\n/g; 
  print $outh ">$pubid\n$faseq\n";  $ntrout++;
  #1806: change here, print tblh '>Features\tpubid' for each fsa >pubid, to ensure same entries in both files
  print $tblh ">Features\t$pubid\t$GBPROID\n\n"; #?? is this used now
  my $ctab="\t\t\t";
  
  # FIXME: for ncRNA, need other kind of annot, $cdsline is 'begin end CDS' .. 
  if($ncRNAok){ 
    # valid annots: maploc, orig_id, maybe dbxref
    # chrmap:99a,90i,0l,4x,chr20:10736002-10740141:+,mol:ncRNA
    
    print $tblh join("\t",1,$falen,"$seqtype")."\n";
    my($mapnote)= tblMaploc($tblinfo);
    print $tblh $ctab,"map\t",$mapnote,"\n" if($mapnote); #? /note= or /map= qualifier?
    
    print $tblh $ctab,"note\toriginal_id:$oid\n" if($oid); # moved to end
    print $tblh "\n";
    $ncrnaout++; #? always same count as ntrout ? no, not cdsOK
  }  
  elsif($cdsOK) {  
    $ncdsout++; #? always same count as ntrout ? no, not cdsOK
    #old. print $tblh ">Features\t$pubid\t$GBPROID\n\n"; #?? is this used now
    # >Features       Thecc1ER_sopcsc10rk23loc140t1     cacao11evigene_20120827

		unless($protid) { ## dont need protid in tblinfo if it is only this transform
			$protid= $pubid;
			$protid=~s/t(\d+)$/p$1/; # genbank requires diff protid from mrnaid
			$protid= $GDB_PREFIX.$protid if($GDB_PREFIX); #  protein_id  gnl|CacaoGD|Thecc1EG016762p1
			}
    #no# $protids{$oid}= $protid;

    print $tblh $cdsline;
    #o.print $tblh $ctab, "protein_id\t$protid\n" if($protid); # move to end? or option?

    # print $tblh $ctab,"locus_tag\t$lotag\n" if($lotag); # problems, needs gene entry for this..
    if($lotag) {
      my($alt)= ($pubid=~m/t(\d+)$/)?$1:0; # check for >1 alt?
      my $v="locus_tag:$lotag"; $v.=", isoform $alt" if($alt);
      print $tblh $ctab,"note\t$v\n";
    }
    #o.print $tblh $ctab,"note\toriginal_id:$oid\n" if($oid); # move to end

    ## print productname near top of ann list
    ## productname does: my $gnn=($gname eq $NAME_UNK or not $gname)?$NAME_UNKNCBI:$gname;
    print $tblh $ctab,"product\t$gname\n"; 
    print $tblh $ctab,"note\tconserved domain $cddname\n" if($hascddname); # use ORIG cdd not prodname
      ## cdd note: XXX domain-containing protein >> XXX 

    if($gname and $namepct >= $MIN_IDLIKE) {  # MIN_IDLIKE MIN_NAMEIDENT
    
      #UPD108: change to "product alignment is 98 to human gene [GSYM?] NP_nnnnn" or similar, merge db_xrefs
      # "protein is 98 similar to human gene NP_nnnn"
      #x my $pnote="product similarity is $namepct"; # alignment/similarity?
      #x if($gnameref) { $pnote.=" to $gnameref"; }  
      
      my $pref= ($gnameref)? $gnameref : "reference"; # want user/config opts for label here, e.g.
        # protein is 99 pct similar to human RefSeq:NP_NNNNN gene; to pig RefSeq:NP_XXXX gene; ..
        # need some hash of keyNameRefID => valNameclass, tie to DBXREF mess? 
      
      my $pnote="protein is $namepct pct similar to $pref gene."; # alignment/similarity?
      print $tblh $ctab,"note\t$pnote\n"; 
      
      #old1805: print $tblh $ctab,"note\tproduct alignment:blastp is $namepct\n"; # NO: \%
      ## SUSPECT_PRODUCT_NAMES:  8 cds comments or protein descriptions contain '%'
      #see below# inference  similar to AA sequence:Phytozome:PGSC0003DMP400040846
    }

    ## add more notes from ann.txt: maploc, mapqual
    # # @ANNO_TBLINFO = qw(pubid oid trlen cdsoff aaqual trgaps dbxref namepct name cdd evgclass maploc mapqual);
    ## FIXME: NOPATH or mapq align < 10% >> note 'unlocated on chromosome assembly'
  if(1) {
    my($mapnote)= tblMaploc($tblinfo);
    print $tblh $ctab,"map\t",$mapnote,"\n" if($mapnote); #? /note= or /map= qualifier?
  
  } else {
    if(my $mapl= $tblinfo->{'maploc'}) {
      my $mapq= $tblinfo->{'mapqual'} || "";
      my ($maln, $midn, $mlen, $mxon)=( 0) x 9;
      if( $mapq=~m/(\d+)a,/ ) { 
        $maln=$1; my $mc=($maln > 89)?"":($maln < 10)?"":"$maln pct aligned,";
        $mapq=~s/(\d+)a,/$mc/; 
        }
      $mapq=~s/,(\d+)i//; $midn= $1; # drop, or $1-identity, 
      $mapq=~s/,(\d+)l//; $mlen= $1; # drop length
      $mapq=~s/,(\d+)x/, $1 exons/; $mxon=$1;
      $mapl=~s/:[fr\+\.\-]$//; # drop strand always?
      my $mnote= ($maln < 10)? "unlocated on chromosome assembly" : "located at $mapl, $mapq"; 
      print $tblh $ctab,"note\t",$mnote,"\n"; # map quality is $mapq
    }
  }
  
    if($Selcstop) { # add /transl_except=(pos:1002..1004,aa:Sec);
      my @ss=split",",$Selcstop;
      for my $sb (@ss) { my $se=$sb+2; print $tblh $ctab,"transl_except\t(pos:$sb..$se,aa:Sec)\n" if($sb>0); }
    }

=item DBXREF config

    see evigene2genbanktbl.pl
		evigene.conf:
			dbxref_recode = TAIR=arath Phytozome=poptr|vitvi|soybn|soltu|sorbi
			dbxref_other DBXMISSING  

   		dbxref_ok  CDD taxon DDBJ EMBL NCBI GeneID GI dbEST dbSNP TrEMBL Swiss-Prot  UNILIB InterPro ENSEMBL GO 
  	   + JGIDB ESTLIB APHIDBASE AntWeb BEETLEBASE dictyBase FLYBASE GDB GeneDB GRIN MGI PDB PFAM PGN SGD SGN TAIR VectorBase WormBase ZFIN

		dbxref_ok from http://www.ncbi.nlm.nih.gov/genbank/collab/db_xref
			-- this html not easily parsible.
    see also ncbicsrc/api/asn2gnb6.c superset of http://www.ncbi.nlm.nih.gov/collab/db_xref.html
			
		init db_xref:
			my $DBXOTHER= $configh->{general}->{dbxref_other} || "DBXMISSING"; # special code or blank?
			my %DBXREF_RECODE= (); # dbxref_recode = TAIR=arath Phytozome=poptr|vitvi|soybn|soltu|sorbi
			my %DBXREF_OK= ( taxon=>1, TAIR=>1, TrEMBL=>1, SwissProt=>1, ENSEMBL=>1, ); 
			if( my $dbxrefok= $configh->{general}->{dbxref_ok}) {
				my @pg= split/[,\s]+/, $dbxrefok; 
				map{ my($k,$v)=split/[=:]/,$_; $v=1 unless(defined $v); $DBXREF_OK{$k}=$v; } @pg; 
			}
			
			if( my $dbxrefre= $configh->{general}->{dbxref_recode}) {
				my @pg= split/[,\s]+/, $dbxrefre; 
				map{ my($v,$k)=split/[=:]/,$_; my @k=split/[|]/,$k; foreach my $k1 (@k) { $DBXREF_RECODE{$k1}=$v; } } @pg; 
			}

		process db_xref:
      my($d)= m/^(\w+):/; 
      if( $DBXREF_RECODE{$d} ) { $d= $DBXREF_RECODE{$d}; }
      unless($d) { $_= "$DBXOTHER:$_"; } elsif(not $DBXREF_OK{$d}) { s/$d:/$DBXOTHER:$d./; } 

=cut
	
    if($dbxref =~ /\w/) {
    	my @dbother=(); # dont pass NCBI valid db_xref:
      foreach my $dx (split",",$dbxref) { 
        # DBXREF_RECODE fix for ncbi's  dbx: restrictions
        # FIXME999: more problems w/ gene.names table having odd/local DBprefix:ID
        #   .. fix where? should have here list of valid NCBI/ISxxx db prefixes. from where?
        # .. isnameref not seen , namepct '%' problem? or gnameref wrong/miss?
        ## must be bad:  $dx eq $gnameref; $gnameref == "RepID:xxx,RefID:yyy"
        my $isnameref= ($namepct >= $MIN_NAMEIDENT and $gnameref =~ m/^$dx/)?1:0; # $dx eq $gnameref
        
        # ? is HUGO:symbol valid at ncbi? hugo HGNC:idnum is valid, not symbol?
        
        $dx =~ s/^UniRef/SwissProt:UniRef/; ## UniProtKB:UniRef/; # or /TrEMBL:UniRef/;  SwissProt
        $dx =~ s/UniProt:/SwissProt:/; ## TrEMBL: .. UniProtKB: # need list to check/replace as per 
        # is CDD:id ok? accepted by tbl2asn .. YES
        
        ## ==> catfish1all4cf/okayset/catfish1all4.mrna.tbl2asn.sumval <==
        ## 60800 WARNING: SEQ_FEAT.IllegalDbXref
        ##  db_xref DRERI:ENSDARG00000019365 ## >> ENSEMBL:  FIXME in where? catfish/zfish.names
        ## anything:ENS\w+\d\d+ should become ENSEMBL:
        
        if($dx =~ /:ENS/ and $dx !~ /^ENSEMBL:/) { $dx =~ s/^\w+:(ENS\w+\d\d+)$/ENSEMBL:$1/; }
        
        ## fixme: somedb:XP_ NP_ are ncbi prot ids;   db_xref_other:mayzebr:XP_004557409.1;
        if($dx =~ /:[XN]P_/ and $dx !~ /^GenBank:/) { $dx =~ s/^\w+:/GenBank:/; }
        #x if($dx =~ /:[XN]P_/ and $dx !~ /^GenBank:/) { $dx =~ s/^\w+:([XN]P_\w[\.\w]+)/GenBank:$1/; }

        #old:if($dx =~ /:[XN]P_/ and $dx !~ /^GeneID:/)..
	      # ^^ UPD201705: new tbl2asn whines GeneID for AA inference, but not for db_xref ??

        ## FIX2: arath:AT5G60040.1 << should be TAIR: dammit ; fix .names instead of here..
        # evigene2genbanktbl.pl sub reformatval()
 
        ## print $tblh $ctab,"db_xref\t$dx\n" if($dx=~/\w/ and $dx ne "na"); 
	      if($dx=~/\w/ and $dx ne "na") { 
		      my($db,$did)= ($dx =~ m/^(\w[^:]+):(.+)/)?($1,$2):("ref",$dx); 
		      #^^ FIXME no dbprefix: .. guess from dx ID?
	        if( $DBXREF_RECODE{$db} ) { $db= $DBXREF_RECODE{$db}; }
	        $did=~s,/\d+,,; # ugh, dromel:FBgn0036451/46, dappu1:E9GZG3_DAPPU/100 extra junk
	        my $dxr="$db:$did";
      	  if( $DBXREF_OK{$db} ) { 
	      	  # UPD201705: new tbl2asn whines about db_xref GenBank|NCBI|GeneID:NP_nnn prot dbxrefs
      	    if($db eq "GenBank") { push @dbother, $dxr;  } # ugh...GenBank: is bad for db_xref, okay for inference
	      	  else { print $tblh $ctab,"db_xref\t$dxr\n"; }
	      	  # .. change db_xref to note 
            print $tblh $ctab,"inference\tsimilar to AA sequence:$dxr\n" if($isnameref);
	      	} else { 
	      	  push @dbother, $dxr; 
	      	}
	      }
      }
      if(@dbother) { # print as notes
    		map{ print $tblh $ctab,"note\tdb_xref_other:$_\n"; } @dbother;
      }
    }
    
    print $tblh $ctab,"protein_id\t$protid\n" if($protid); # move to end? or option?
    print $tblh $ctab,"note\toriginal_id:$oid\n" if($oid); # moved to end
    print $tblh "\n";
  } #end cdsOK
  
  return ($ntrout,$ncdsout,"OK",$protid);  # $ncrnaout
}

=item cdsloc bug bad b>e

>Features       vDanrer6pEVm000001t22   vDanrer6pEVm_20181005
21521   >3      CDS
                        product Titin-like protein
                        note    protein is 23 pct similar to sacavefish16nc:XP_016301102.1 gene.
                        map     located at chr9:43153041-43218249, 99i, 91 exons
                        note    db_xref_other:GenBank:XP_016301102.1
                        note    db_xref_other:ncbi:XM_017357900.99
                        note    db_xref_other:zfin:ZDB-TSCRIPT-090929-9358.76
                        protein_id      gnl|Evigene|vDanrer6pEVm000001p22
                        note    original_id:tridba4b_sRn1l2bSRR1524240ridk45Loc7033

>Features       vDanrer6pEVm000002t3    vDanrer6pEVm_20181005
44936   >3      CDS

 ../publicset/zebrafish17evigene_m6pt.pubids
vDanrer6pEVm000002t3    Danrer6pEVm000002t1d57  vDanrer6pEVm000002      3       alt     14978,99p,partial3      0       0,0,chrmap:100a,99i,45126l,74x,chr9:42944188-42998028:- Danrer6pEVm000002t1d57,tridba1a_sNn12l1SRR1524240ridk71Loc32051

../zebrafish17evigene_m6pt.ann.txt >> CDSoff end-begin
vDanrer6pEVm000002t3    Danrer6pEVm000002t1d57  45126   44936-3 14978,99p,partial3      0       zfish16nc:XP_017213390.1,ncbi:XM_017357901/100.99,zfin:ZDB-TSCRIPT-090929-22794/99.99,  55%,14978/27144,14978   Titin, partial  na      alt

../publicset/zebrafish17evigene_m6pt.public.aa
>vDanrer6pEVm000002t3 type=protein; Dbxref=zfish16nc:XP_017213390.1,ncbi:XM_017357901/100.99,zfin:ZDB-TSCRIPT-090929-22794/99.99; 
  aalen=14978,99%,partial3; clen=45126; offs=44936-3; 
  evgclass=alt; oid=tridba1a_sNn12l1SRR1524240ridk71Loc32051; pubid=nopd,dropalthi1
 
=cut

sub putCDSloc 
{
  my($cdsoff,$partial,$cdsphase,$mrnalen,$mrna,$oid)= @_;
  $oid ||= "nada";
  ## is cdsoff == "<123->456" allowed here?
  #ALREADYdone: $cdsoff =~ s/:[+.-]$//; # drop strand if there.
  my($start,$stop)= split/[-]/,$cdsoff; # or $cdsoff =~ m/(\d+)-(\d+)/; # 
  ## ^^bug? start > length(mrna) below, bug from what? trim(mrna)?
  if($stop and $start > $stop) { ($start,$stop)=($stop,$start); $cdsoff="$start-$stop"; } # FIX in putCDSloc or here?
  
  my($pstart,$pstop)= ($start,$stop);
  my($p5,$p3,$codonstart,$mrnatrim,$mrnadif)= ("<",">",0,0,0);
  if($partial =~ /complete/) { 
    # need to check start-stop are inside mrna, some utrorf bugs..
    my $upd=0;
    # if($start < 1) { $upd=1; } # partial5
    if($stop > $mrnalen) { $stop= $mrnalen; $upd=3; } #? partial3 now?
    if($start > $stop) { $start=$stop=0; $upd=4; } # 0 will cancel putCDSloc(), some data bug? or set start=1
    if($upd){
      $cdsoff="$start-$stop";
      ($pstart,$pstop)= ($start,$stop);
      if($upd==3) { $partial="partial3"; }
      }
  } 
  
  unless($partial =~ /complete/) {
    ## FIXME.1510: NCBI requires mrna-cds partial to abut ends of mrna,unless gaps, 
    ## ie. ofs=mrnastart+2,3 wrong and ofs=mrnaend-1,2 wrong; for start shift, also shift cdsphase/codonstart
    ## -- problem with this, of course, when mrna-*end* is +2 from partial3-cds end,
    ##   cannot extend cds to mrna end w/o calling 2-base codon, adding to protein :((
    ## -- answer: trim mrna-end to cds-end (only for +2 case or +1 also?)
    
    unless($partial =~ /partial[53]/) { $partial.="53"; } # both
    
    if($partial =~ /3/) {  
      if($stop < $mrnalen and $stop > $mrnalen - 3) { # was stop > len-3
        $mrnatrim= $mrnalen - $stop;  ## FIXME.1510:
        # my $ix= index(substr($mrna,$mrnalen - 3), 'N');
        
        #upd1806: aaendgap fix, if aa.X at end, and mrna N at end, trim both
        # problem may be only for ni == 1 ? end gap
        if($AAGAPENDS ) { 
          #?? bug: and  $AAGAPENDS->{$oid} .. bug2: and $AAGAPENDS{$oid} < not found
          my $ni= (substr($mrna,$stop-2,1) eq 'N') ? 2 
                : (substr($mrna,$stop-3,1) eq 'N') ? 3  
                : (substr($mrna,$stop-1,1) eq 'N') ? 1 : 0; # stop-1=N bad for tbl2asn, stop-2/-3 ok
          if($ni == 1) { # any of 1,2,3? or only 1?
            $stop -= $ni; $mrnatrim += $ni; $cdsoff="$start-$stop"; 
            $AAGAPENDS{$oid}= 2; # this works at putpep, not incoming here
		        loggit(0, "aagaptrim", $oid, "aa-$ni, rna3end=", substr($mrna,$mrnalen-6)); # LOG_DEBUG
          } 
        }  
        if(0) { 
          # Causes more problems : Warn: SEQ_FEAT.TerminalXDiscrepancy mismatch prot,cds..
          # without get  Warn: SEQ_INST.TerminalNs for same prot/cds set; 
          # problem derives from prior fix, gap trims of AGCTnnnnnGnnnnnA left  AGCTnGnA un-cut single n's
          if(substr($mrna,$stop-1,1) eq 'N') { $stop--; $mrnatrim++; $cdsoff="$start-$stop"; } ## was below: if($faseq=~m/[ACGT]N$/)
          }
      }
      elsif($stop <= $mrnalen - 3 and $mrna) { 
        ##?? problem w/ this case, trim inner=mrna-NNN to cds-end may be large ??
        # if(substr($mrna,$stop,1) eq 'N') { $eshift=1; } ## == stop-1 + 1, zero-origin
        # elsif(substr($mrna,$stop+1,1) eq 'N') { $eshift=2;  }  
        # $stop += $eshift; $cdsoff="$start-$stop"; 
        }
      # $mrna= substr($mrna, 0, $mrnalen - $mrnatrim);  # leave to caller.
      $pstop="$p3$stop"; 
      }
      
    if($partial =~ /5/) { # no problem if called here already?
      if($start > $mrnalen) { $start=$stop=0; } # 0 will cancel putCDSloc(), some data bug? or set start=1
      elsif($start > 1 and $start < 4) { $cdsphase += $start-1; $start=1; } ## FIXME.1510:
      elsif($start >= 4 and $mrna and substr($mrna,$start-2,1) ne 'N') {
        ## buggggg this is too messy, still ncbi checha errs; dont shift mrna but change leading 1,2 to NN?
        ## bugg here, offby1, NNNN preceding, only shift for what?
        if(substr($mrna,$start-3,1) eq 'N') {  ## == start-1 -1, zero-origin
          #x# $start -= 1; $cdsphase= ($cdsphase-1) % 3; 
          $mrnadif=1; substr($mrna,$start-2,1)= 'N'; # try this way
          } 
        elsif(substr($mrna,$start-4,1) eq 'N') { # shift-2
          #x# $start -= 2; $cdsphase= ($cdsphase-2) % 3; 
          $mrnadif=2; substr($mrna,$start-2,1)= 'N'; # try this way
          substr($mrna,$start-3,1)= 'N'; # try this way
          } 
      }
      $cdsoff="$start-$stop";
      
#////// bad vvv /////////////
#       my $bshift=0; 
#       if($start > 1 and $start < 4) { $bshift=$start-1; } ## FIXME.1510:
#       elsif($start >= 4 and $mrna) {
#         if(substr($mrna,$start-2,1) eq 'N') { $bshift=1; } ## == start-1 -1, zero-origin
#         elsif(substr($mrna,$start-3,1) eq 'N') { $bshift=2;  } # shift-2
#       }
#       if($bshift>0) {
#         # $cdsphase += $bshift; # ugh wrong way? no; ugh ugh .. bshift wrong way for start>4 vs start 2,3
#         $cdsphase -= $bshift; $cdsphase = $cdsphase % 3; # convert -2 > 1, -1 > 2 << only for start>3
#         $start -= $bshift; $cdsoff="$start-$stop"; 
#       }
#--------------      
      $codonstart=$cdsphase+1; # is this right?
      $pstart="$p5$start";  
    }
  }
  my $tbl= join("\t",$pstart,$pstop,"CDS\n");
  $tbl .= "\t\t\tcodon_start\t$codonstart\n" if($codonstart>0);
  return ($tbl,$mrnatrim,$cdsoff,$cdsphase,$mrnadif,$mrna);## NOT: ,$mrna
}


# evgmrna2tsa2.pl get_evgtrset(); move to package?
sub get_evgtrset { 
  my($trclass,$cdnaseq,$pubdir)= @_;
  my($trpath,$trname,$nsra,$sradatah)=("","",0,undef);

  ## getmRNA() without make_mrna portion, find files only..
  # FIXME for mrna_cull.fa , other variants in publicset/
  unless($cdnaseq and -s $cdnaseq) { 
    if($pubdir and -d $pubdir) { 
      my($pubd)= getFileset($pubdir,'fa|mrna'); 
      my($trf)= grep /\.mrna_pub\.fa|\.mrna$|\.mrna.gz$/, (@$pubd);
      if($trf and -s $trf) { $cdnaseq= $trf; }
      }
  }

  ## publicset/pig4321ew.mrna_pub.fa  >> trname & oname= pig4321ew

  ## FIXME, return trpath = publicset/.. one above it to find data
  if($cdnaseq and -s $cdnaseq) { 
    ($trpath,$trname)= ($cdnaseq =~ m,(.*)/([^/]+)$,)? ($1,$2) : ("./",$cdnaseq);
    $trpath =~ s,[/]?$pubdir,, if($pubdir);
    $trpath="./" unless($trpath);
  }
  # want this? drop suffix
  unless( $trname =~ s/\.mrna_pub\.fa|\.mrna// ) {
    $trname=~s/\.\w+$//; 
  }
  
  # return many @mrna ? pub, cull, drop, xxx
  return($cdnaseq,$trpath,$trname,$sradatah); 
  
#   if($cdnaseq) { 
#     $notokay=1; # dont look in okayset/? look in $pubdir now?
#     $trclass= makename($cdnaseq,".trclass") unless($trclass); 
#   }
#   if($trclass) {
#     my $trpname= makename($trclass,"","trclass"); 
#     if($trpname =~ m,/,) { ($trpath,$trname)= $trpname =~ m,(.*)/([^/]+)$,; } # was BAD?
#     else { $trname=$trpname; }
#     $trpath ||= '.';  
#     my $okpath= ($notokay) ? $trpath :"$trpath/okayset"; # should I check if curdir has okayset files?
#     #OBSOLETE# my($mrnaOfUtrorf,$nutrorf)= getmRNA_utrorf($okpath,$trname);
#     ($cdnaseq)= getmRNA($okpath,$trname,$pubdir) if(!$cdnaseq and -d $okpath);
#   }
  
}


sub set_setting { my($key,$val)=@_; $settings{$key}= $val; }
sub set_newsetting { my($key,$val)=@_; $settings{$key}= $val unless($settings{$key}); }
sub get_setting {
  my($key)=@_;  my $v= $settings{$key}||undef; # add? || $ENV{$key}
  return $v;
}

## another common sub; see similar cdna_evigenesub.pm:evigene_config()
sub do_settings {  
  my($action,$pathname)= @_;
  ## write these to work dir; reread next go
  ## action == 'log|save|restore' ; restore called AFTER read new options, shouldnt replace

  my $PRES=$EGLOG; # == 's2g';
  my $trpname= makename($pathname,".$EGAPP.info"); # ".sra2genes.info"   ## dang
  my $runpath= dirname($pathname); # == '.'
  unless($runpath =~ /\w/) { $runpath= `pwd`; chomp($runpath); } # ok?
  
  #?? $sraids=~s/ +;/;/g; $sraids=~s/ +/;/g; # ?? need this

	## merge this and get_srainfo() from sra.cvs .. allow fields in both?
	## $sradatah->{'assemblers'} == "Velvet/Oases v; SOAPDenovoTrans v; Trinity v;"
  ## $settings{runpath} = $runpath == `pwd`;   ## add runpath == datad somewhere, here?
  
  my %mysettings= map{ $_ => $DEFAULT_SETTINGS{$_} } (keys %DEFAULT_SETTINGS); # ensure these are in??

  # globals, yuck .. do away and use only global %settings ?
  %mysettings=( 
		IDPREFIX=>$IDPREFIX,  DATE=>$DATE, 
		sraids => $sraids, 
		organism => $ORGANISM, BioProject => $BioProject,  
		TSADESC => $TSADESC, TSAOPTS => $TSAOPTS,
    sratable => $sratable, 
		# runpath => $runpath,
		# MINSIZE=>$MINSIZE, MAXGAP=>$MAXGAP,
		# trclass => $trclass, mrna => $cdnaseq, genenames=>$genenames,
		);
	
  if($action =~ /restore|read/ and -f $trpname) {
    open(my $inh, $trpname); # or loggit(warn ..);
    while(<$inh>) { chomp; if(s/^$PRES.//) { 
      my($k,$v)=split /\s*=\s*/,$_,2; 
      $v=~s/;*$//; $v=~s/\s*\#.*$//; # should I chomp trailing '#' comments? 
      my ($v1)= ($v=~/;/)? (split";",$v)[0] : $v; # ugh: organism=Bob_white;Bob_white_black;.. from sra.csv mixtures
      my $ov= $mysettings{$k}; my $dov= $DEFAULT_SETTINGS{$k}||"";
      unless($ov and $ov ne $dov) { 
        ## fixme: need to reset global defaults
        $ORGANISM=  $v1 if($k eq 'organism'); # ONE only, not orga;orgb;orglist;
        do { $v=~s/ +;/;/g; $v=~s/ +/;/g; $sraids= $v; } if($k eq 'sraids');
        $BioProject= $v if($k eq 'BioProject'); # FIXME lost
        $IDPREFIX=  $v if($k eq 'IDPREFIX');
        $DATE=      $v if($k eq 'DATE');
        $TSADESC=   $v if($k eq 'TSADESC');  
        $TSAOPTS=   $v if($k eq 'TSAOPTS'); # 
        $sratable=   $v if($k eq 'sratable');
        # $MINSIZE=   $v if($k eq 'MINSIZE'); # require int
        # $MAXGAP=    $v if($k eq 'MAXGAP'); # require int
        # $genenames=  $v if($k eq 'genenames');
        # $trclass=  $v if($k eq 'trclass');
        # $cdnaseq=  $v if($k eq 'mrna'); # fixme: key != varname
        
        $mysettings{$k}=$v; # after possible changes
        } 
      }
    } close($inh);
    
    %settings = %mysettings; # make global now; ONLY after restore|read ??
    #?? @sraids= grep /^[A-Z]+\d+/, split /\W+/, $sraids; # NOW GLOBAL; should be SRRnnnn ERRnnnn DRRnnn 

  } else { # action==save; copy new single vars to settings..
    for my $ks (sort keys %mysettings) { 
      my $ov= $mysettings{$ks};  my $dov= $DEFAULT_SETTINGS{$ks}||"";
      $settings{$ks}= $ov if($ov and $ov ne $dov);    
      } ;
  }

  my $settings= join "\n", map{ "$PRES.$_=".$settings{$_} } sort keys %settings;
  if($action =~ /log/) { loggit(0, "$EGAPP.info:\n$settings");  }
  if($action =~ /save/) { open(my $outh, '>', $trpname); print $outh $settings,"\n"; close($outh); }
}

sub get_srainfo {
  my($sracvs,$sraidlist)= @_; # new
	# my($sracvs,$trpath,$trname)= @_; # orig
	
  # unless($sracvs and -s $sracvs) {
  #   my $csv="$trpath/$trname.sra_result.csv";    
  #   $csv="$trpath/sra_result.csv" unless(-f $csv);  
  #   $sracvs=$csv if(-f $csv);
  # }
	
	my($nsra,$sradatah)=(0,0);
	if(-f $sracvs) {
	  ($nsra,$sradatah)= parse_sra_result_cvs($sracvs);
	} else {
    $sradatah= {};
    # make dummy sra.cvs for other components?? pubset2submit.pl wants it
    $sradatah->{cvsformat}= 0;
	  if($sraidlist) {
      ##$sraids= $sraidlist; #? or always use sids?
      my @sraids= grep /^[A-Z]+\d+/, split /\W+/, $sraidlist; # ?? GLOBAL; should be SRRnnnn ERRnnnn list
      $nsra= @sraids;
     }
	}
	loggit(0,"sra_result from",$sracvs,"nresult=",$nsra);
	
	my @SRAK=();
	if($sradatah->{cvsformat} == 1) {
  	@SRAK= ("Assemblers", "Instrument", "Total Size, Mb","Total Spots",'Total Assemblies');
	} else {
  	@SRAK= ("Assemblers", "Platform", "size_MB","spots","Total Assemblies", "BioProject");
	}
	for my $ks (@SRAK) {
  	if($settings{$ks} and not $sradatah->{$ks}) { $sradatah->{$ks} = $settings{$ks}; }
  	elsif($sradatah->{$ks}) { $settings{$ks} = $sradatah->{$ks}; }
	}
	
	##* moved out of srainfo
  # my($nupinfo,$tsamethf,$tsadescf,$tsasubf)= tsa_infotemplates($trpath, $trname, $sradatah);
  # loggit(0,"info updated $nupinfo TSADESC=",$TSADESC);
  # 
  # unless($genenames) { #? put in get_evgtrset
  #   my $gnt="$trpath/$trname.names";  $genenames=$gnt if(-s $gnt); 
  #   loggit(LOG_WARN, "Missing product -names",$gnt) unless($genenames);
  # }
	
	return($sradatah);
}

use constant FIXME_FTP => 1;

#* see evgpipe_sra2genes.pl:parse_sra_result_cvs() ; should be same sub now
sub parse_sra_result_cvs
{
  my($sracvs)= @_;
  my($ngot,$nin,$nerr,$cvsformat)=(0) x 9;
  my (%sradata, @hd);
  
  # revise for this format; change old format keys to new ??
  #? are these fields fixed now? need to check if no hdr in sracvs
  my $SRA2CVSHDR='Run,ReleaseDate,LoadDate,spots,bases,spots_with_mates,avgLength,size_MB,AssemblyName,download_path,Experiment,LibraryName,LibraryStrategy,LibrarySelection,LibrarySource,LibraryLayout,InsertSize,InsertDev,Platform,Model,SRAStudy,BioProject,Study_Pubmed_id,ProjectID,Sample,BioSample,SampleType,TaxID,ScientificName,SampleName,g1k_pop_code,source,g1k_analysis_group,Subject_ID,Sex,Disease,Tumor,Affection_Status,Analyte_Type,Histological_Type,Body_Site,CenterName,Submission,dbgap_study_accession,Consent,RunHash,ReadHash';
  
  open( my $inh, $sracvs) or loggit(1,"ERR: parse_sra_result_cvs reading $sracvs");
  while(<$inh>) { 
    next unless(/^\w/ and (/,/ or /\t/)); chomp; 
    if($cvsformat == 0) { # always/no header?
      if(/^Run,ReleaseDate/){  $cvsformat=2; }
      elsif(/^[DES]RR\d+,/){ $cvsformat=2; # SRR,DRR,ERR .. others?
        @hd=split",",$SRA2CVSHDR; # guess this col set
      }
      elsif(/Experiment Accession/){ $cvsformat=1; } # "Experiment..","xxx","yyy"
      elsif(/^\w+\t/){ $cvsformat=3; } # tabbed ?
      $sradata{cvsformat}= $cvsformat;
    }
    
    my @cols;
    if($cvsformat == 3) {
      @cols= split"\t",$_;
      if(/^Run/){ @hd=@cols; next; }
    } elsif($cvsformat == 2) {
      s/","/\t/g;  s/^"//; s/"$//; #quotes or not?
      s/,/\t/g;  @cols= split"\t",$_;
      if(/^Run/){ @hd=@cols; next; }
    } elsif($cvsformat == 1) {
      @cols= map{ s/^"//; s/"$//; $_; } split /\",\s*\"/, $_;
      if(/^Experiment Accession/){ @hd=@cols; next; }
    } else {
      $nerr++; # warn/die ??
      next;
    }
    
    $ngot++;
    for(my $i=0; $i<@cols; $i++) { 
      my $hd=$hd[$i]||$i; my $v=$cols[$i]; 
      ## redo sradata by SRRid ? sradata{col}{sid}=val ?
      $sradata{$hd}.="$v;" 
        unless($sradata{$hd} and $sradata{$hd} eq "$v;"); 
            # ^^^ problem for mix of same/diff vals, many entries
      }
  } close($inh);
  
  foreach my $k (sort keys %sradata) { $sradata{$k}=~s/;$//; }
  
  #reset defaults:  
  # my($deforg,$defsra)=("Noname","SRR000000"); 
  my($deforg,$defsra)= ($DEFAULT_SETTINGS{'organism'},$DEFAULT_SETTINGS{'sraids'});
  my $sorg= $sradata{"ScientificName"} || $sradata{"Organism Name"} ||"";
  my $sids= $sradata{"Run"} || $sradata{"Experiment Accession"} ||"";
  # globals: set here?
  $sorg=~s/;.*//; # fixme for multiple orgs: orga;orgb;orgc;...
  $sorg=~s/ /_/g; 
  $ORGANISM= $sorg if($sorg and ($ORGANISM eq $deforg or $ORGANISM !~ m/\w\w/));
  $sraids=   $sids if($sids and ($sraids eq $defsra or $sraids !~ m/SRR/)); #? or always use sids?
  @sraids= grep /^[A-Z]+\d+/, split /\W+/, $sraids; # NOW GLOBAL; should be SRRnnnn ERRnnnn list
  
  if($sorg) {
    my($gn,$sp)= split/[_\W]/,$ORGANISM,2;
    if($sp and $gn) {
      my $shortorg= ucfirst(substr($gn,0,3)) . lc(substr($sp,0,3));
      #? $shortorg .= 'EVd';
      # $settings{oname}= $shortorg;
      set_newsetting('oname',$shortorg);
      $IDPREFIX= $shortorg.'EVm' if($IDPREFIX eq $DEFAULTidpre);
    } else {
      set_newsetting('oname', substr($ORGANISM,0,8));
    }      
  }
  ## also use sradatah to edit template evgr_tsamethods.cmt, evgr_tsadesc.cmt, evgr_tsasubmit.sbt ??
  ## rewrite template .cmt, .sbt unless updated already.

  ## ALSO use 'FTP Path to Experiment' to turn SRX into SRR for picky ncbi tsa : wget or curl calls?
  my $HAVEsrrids= ($sraids =~ /[DES]RR\d/ and $sraids ne $defsra)?1:0; ## Need to save to info/config file
  unless(FIXME_FTP or $HAVEsrrids) {
    my @urls;
    @urls= grep /http|ftp/, split/;/, $sradata{"download_path"} || $sradata{"FTP Path to Experiment"};  
    my @srr= sra_ftp2srr(@urls);
    if(@srr>0) {
      $sraids= join(";",@srr);
      $sradata{'SRAids'}= $sraids;
      loggit(1,"sra_id=",$sraids);
      }
  }
   
  return($ngot, \%sradata, $sraids); # sids ?
}

sub sra_ftp2srr {
  my(@ftps)= @_;  return () unless(@ftps);
  my @srrs=(); 
unless(FIXME_FTP) {  
  my $APPcurl= findapp('curl'); return () if($APPcurl =~ /MISSING/);
  foreach my $ftppath (grep /^ftp:/, @ftps) {
  	my $cmd="$APPcurl -s -l $ftppath/"; 
  	loggit( ($dryrun) ? 1 : 0,"CMD=",$cmd);  
    my $srrs= `$cmd`; # or http: ?? ## runcmd() instead ?
    push @srrs, grep /SRR/, map{ m/(SRR\w+)/; $1; } split " ",$srrs;
    }
}    
  return @srrs;
}


#.......... tsa/tbl2asn subs : move out to separate app ?? ...................

=item tsa_infotemplates

  - Make these if needed; update for trasm name as needed

evgr_tsamethods.cmt:
  StructuredCommentPrefix ##Assembly-Data-START##
  Assembly Method EvidentialGene v2013.03.02; << from VERSION
     Velvet/Oases v1.2.03/o0.2.06 (2012.02); SOAPdenovo-Trans v2011.12.22; Trinity v2012.03.17  << USER inputs
  Assembly Name   evigeneR13     << from local myspecies.trclass name
  Sequencing Technology   Illumina   << from sra_result

evgr_tsadesc.cmt : insert parts of sra_result
  RNA-Seq data of [[this species]] are
  assembled de-novo with [[various RNA assemblers]], using multiple options
  for kmer fragmenting, insert sizes, read coverage, quality and
  abundance filtering. 
  EvidentialGene tr2aacds pipeline software is used
  to process the [[several million]] resulting assemblies by coding
  sequences, translate to proteins, score gene evidence including
  CDS/UTR quality, homology, and classify/reduce into a biologically
  informative transcriptome of primary and alternate transcripts.

evgr_tsasubmit.sbt
Submit-block ::= {
  contact { contact { name name { ... } }
  ..
  Seqdesc ::= user { type str "DBLink",
  data { { label str "BioProject", num 1, data strs { "PRJNA99999" } } } << need PRJNA here..
  }
  

=cut

sub tsa_infotemplates  # make if not found as files; edit w/ sra_result info, other
{
  my($trpath, $trname, $sradata)= @_; ## add sraresult hash info
  my($methtxt,$desctxt,$subtxt)=("") x 3;
  my($tsamethf,$tsadescf,$tsasubf)=("") x 3;
  my $update=0; my $nupdate=0;
  my $NEEDUPDATE = ! $skipTSAparts;
  
  # FIXME: must have \tabs in StructuredComment
  # FIXME paths: trpath may be submitset/ for 2nd runs..
  my $tpath=$trpath;
  my($okd,@files)= getFileset($tpath,'cmt|sbt');  # check for trname in file set
  unless(@files) {
    $tpath =~ s,\w+$,submitset,; 
    $tpath= $trpath."/submitset" unless(-d $tpath);
    ($okd,@files)= getFileset($tpath,'cmt|sbt');  
  }
  
  # ($tsamethf)= grep/$trname.tsamethods.cmt/,@files; $update=0;
  # unless($tsamethf) { ($tsamethf)= grep/evgr_tsamethods.cmt/,@files; $update=$NEEDUPDATE; }
  
  $tsamethf= "$tpath/$trname.tsamethods.cmt";  $update=0;
  #? unless( -s $tsamethf) { ($tsamethf)= grep /tsamethods.cmt/, @files; $update=$NEEDUPDATE; }
  unless( -s $tsamethf) {  $tsamethf= "$tpath/evgr_tsamethods.cmt"; $update=$NEEDUPDATE; }
  if( -s $tsamethf) {
    open(F, $tsamethf); $methtxt= join "", <F>; close(F);
  } else {
    $update=$NEEDUPDATE; $methtxt=""; # TEMPlATE
    my $Illumina= $sradata->{"Instrument"} || "Illumina"; # from sradata
    my $asmsoft= $sradata->{'Assemblers'}||""; $asmsoft="; $asmsoft" if($asmsoft); 
    my $asmname="$trname.evigene"; # from $trname
    my $egvers= 'v'.VERSION;
    $methtxt=<<"EOT";
StructuredCommentPrefix\t##Assembly-Data-START##
Assembly Method\tEvidentialGene $egvers$asmsoft
Assembly Name\t$asmname
Sequencing Technology\t$Illumina
EOT
  }
  if($update) {
    $tsamethf= "$tpath/$trname.tsamethods.cmt"; $nupdate++;
    open(O,'>',$tsamethf); print O $methtxt,"\n"; close(O);
	  push @submitset, $tsamethf;  
  }
  
  $tsadescf= "$tpath/$trname.tsadesc.cmt"; $update=0;
  #? unless( -s $tsadescf) { ($tsadescf)= grep /tsadesc.cmt/, @files; $update=$NEEDUPDATE; }
  unless( -s $tsadescf) {  $tsadescf= "$tpath/evgr_tsadesc.cmt"; $update=$NEEDUPDATE; }
  if( -s $tsadescf) {
    open(F, $tsadescf); $desctxt= join "", <F>; close(F);
  } else {
    $update=$NEEDUPDATE; $desctxt=""; # TEMPlATE
    # m2t.Assemblers=Velvet/Oases v1.2.07/o0.2.08 (2012.10); SOAPdenovo-Trans v2011.12.22; Trinity v2012.03.17 ; Cufflinks v1.0.3 (2012.07)
    my $asmsoft= $sradata->{'Assemblers'} || "various RNA assemblers"; #?
    if($asmsoft) { $asmsoft =~ s/\s+v\d[^;,]+//g; $asmsoft =~ s/;/,/g; } ## drop version info, kept above..
## sra2: @SRAK= ("Assemblers", "Platform", "size_MB","spots","Total Assemblies", "BioProject");

    my $datasize= $sradata->{size_MB}||$sradata->{"Total Size, Mb"} || "0"; ## "Total Size, Mb"
    my $nspots= $sradata->{spots} || $sradata->{"Total Spots"} || "0";
    my @ds=split";",$datasize;  my @ns=split";",$nspots; 
    if(@ds>1) { my $ds=0; my $ns=0; 
      for(my $i=0; $i<@ds; $i++) { my $n;
        ($n)= $ds[$i] =~ m/(\d+)/; $ds+=$n; 
        ($n)= $ns[$i] =~ m/(\d+)/; $ns+=$n; 
        } 
      $datasize=$ds; $nspots= $ns; }
    # calc mb from nspots x 100 x 2 ?
    if($nspots>1000 and not $datasize) { $datasize= int($nspots * 200/1000000); } # ok or not?
    if($datasize > 0) {  $nspots||=""; $datasize=", $datasize Mb in $nspots read-pairs,"; }
    else { $datasize=""; }
    my $asmcount= $sradata->{'Total Assemblies'} || "many"; #?
    # m2t.Total Assemblies=5583471

    #was: "EvidentialGene tr2aacds pipeline.."
    $desctxt=<<"EOT";
RNA-Seq data$datasize of $ORGANISM are assembled with
$asmsoft, using multiple options. 
EvidentialGene pipeline is used to process the $asmcount
resulting assemblies by coding sequences, translated proteins, and gene
evidence, then classify/reduce to a biologically informative
transcriptome of primary and alternate transcripts of gene loci.
EOT
  }
  if($update) {
    $tsadescf= "$tpath/$trname.tsadesc.cmt"; $nupdate++;
    open(O,'>',$tsadescf); print O $desctxt,"\n"; close(O);
	  push @submitset, $tsadescf;  
  }

  # check for $ENV{tsasubmit_sbt} ?
  $tsasubf= "$tpath/$trname.tsasubmit.sbt"; $update=0;
  unless( -s $tsasubf) { $tsasubf= $ENV{tsasubmit_sbt} || "$tpath/evgr_tsasubmit.sbt";  $update=$NEEDUPDATE; }
  if( -s $tsasubf) {
    open(F, $tsasubf); $subtxt= join "", <F>; close(F);
  } else {
    $update=$NEEDUPDATE; $subtxt=""; # TEMPlATE cant really guess this. but empty sbt isnt useful
    
    $subtxt=<<'EOT';
Submit-block ::= {
  contact {
    contact {
      name name { last "Gilbert", first "Donald" },
      affil std { affil "Indiana University", sub "Indiana", country "USA", email "gilbertd@indiana.edu" }
    }
  },
  cit {
    authors {
      names std { 
      {
      name name { last "Gilbert", first "Donald" }, 
      affil std { affil "Indiana University", sub "Indiana", country "USA" }
      } 
    }
    }
  },
  subtype new
}

Seqdesc ::= user {
  type str "DBLink",
  data { { label str "BioProject", num 1, data strs { "PRJNA99999" } } }
}
EOT
    
  }
  if(my $bioprj= $sradata->{"BioProject"}) { $update=$NEEDUPDATE if($subtxt=~s/PRJNA99999/$bioprj/); }
  
  if($update and $subtxt =~ /\w/) {
    $tsasubf= "$tpath/$trname.tsasubmit.sbt"; $nupdate++;
    open(O,'>',$tsasubf); print O $subtxt,"\n"; close(O);
	  push @submitset, $tsasubf;  
  }

  ## edit TSADESC if($nupdate or paths not in $TSADESC)
  ## my $TSADESC="-w evgr_tsamethods.cmt -Y evgr_tsadesc.cmt -t evgr_tsasubmit.sbt"; # FIXME: where from?
  unless($TSADESC =~ m/$tsamethf/ and $TSADESC =~ m/$tsadescf/ and $TSADESC =~ m/$tsasubf/ ) {
    #BAD# $TSADESC="-w $tsamethf -Y $tsadescf -t $tsasubf"; 
    $TSADESC.=" -w $tsamethf" unless($TSADESC =~ s/\-w \S+/\-w $tsamethf/);
    $TSADESC.=" -Y $tsadescf" unless($TSADESC =~ s/\-Y \S+/\-Y $tsadescf/);
    $TSADESC.=" -t $tsasubf" unless($TSADESC =~ s/\-t \S+/\-t $tsasubf/);
  }
  
  return($nupdate,$tsamethf,$tsadescf,$tsasubf);
}


sub tsa_tbl2asn
{
  my($cdnaseq,$cdnatbl,$ORGANISM,$sraids)=@_;
  return unless(-s $cdnaseq);
 
  ## FIXME 2015.01, missing .pep/.pep.report for NCPU case .. need to split those also
   
  ##.. this file part input works..
  # pt=litova1all3.mrnatop2000
  # $ncbin/tbl2asn -i tsa.$pt.fsa -f tsa.$pt.tbl -Z discrep.$pt.log \
  # -a r10k -l paired-ends -Vtb -Mt -XE \
  # -w evgr_tsamethods.cmt -Y evgr_tsadesc.cmt -t evgr_tsasubmit.sbt \
  # -j "[moltype=mRNA] [tech=TSA] $org $sra"
  ##  -r  Path for Results [String]  Optional

  ## fix2: test cdnaseq.gz cdnatbl.gz - gunzip, dont think tbl2asn can read gz
  ## FIXME FIXME FIXME damn paths again
  ## .info has basedir ./ for info files; but above moves them to ./submitset/
  ## should resave info after tidyup w/ new paths..
  ## m2t.TSADESC=-w ./ztick4cvel_pt3.tsamethods.cmt -Y ./ztick4cvel_pt3.tsadesc.cmt -t ./ztick4cvel_pt3.tsasubmit.sbt
  
  # my $wdir= dirname($cdnaseq);
  # my $cdnabase= basename($cdnaseq);
  
  my $spldir= makename($cdnaseq,"_tsasubmit/",'cdna|fasta|fsa|fa');  # ok for _split dir
  # ^^ UNUSED NOW? 
  
  my $tsaopts= $TSAOPTS || "-a r10k -l paired-ends -Vt -Mt -XE"; # read from config.; use TSADESC ??
    ## -Vb = gen genbank gbf, not needed, option; -Vtb == genbank + tsa
  my $tsadesc= $TSADESC; # "-w evgr_tsamethods.cmt -Y evgr_tsadesc.cmt -t evgr_tsasubmit.sbt"; # FIXME: where from?
    ## ^^ need config for these; generate some of this?

    ## FIXME: tbl2asn dies silently if missing files -w .. -Y .. or -t ..
    ## UGH: qsub gave this via env sra='xxx':
    ##   -j '[moltype=mRNA] [tech=TSA] [organism=Pogonus chalceus] [SRA=\'SRR424344;SRR424342;SRR424340\']'
    ## egmrna2tsa.19433.err sh: -c: line 0: unexpected EOF while looking for matching `''
    #er2g: forkCMD= /home/ux455375/bio/ncbi2227/bin/tbl2asn -a r10k -l paired-ends -Vt -Mt -XE \
    #  -w ./pogonus1all3.tsamethods.cmt -Y ./pogonus1all3.tsadesc.cmt -t ./pogonus1all3.tsasubmit.sbt \
    # -j '[moltype=mRNA] [tech=TSA] [organism=Pogonus chalceus] [SRA=\'SRR424344;SRR424342;SRR424340\']' \
    # -i ./okayset/pogonus1all3.mrna_tsasubmit/pogonus1all3.mrna.split1.fsa \
    # -f ./okayset/pogonus1all3.mrna_tsasubmit/pogonus1all3.mrna.split1.tbl \
    # -Z ./okayset/pogonus1all3.mrna_tsasubmit/pogonus1all3.mrna.split1.discrep

  map{ s/^\W+//; s/\W+$//; s/_/ /g; } ($ORGANISM, $sraids); ## $sraids=~s/ +;/;/g; $sraids=~s/ +/;/g; ??
  $sraids =~ s/;/,/g;
  my $tsaqual="-j \'[moltype=mRNA] [tech=TSA] [organism=$ORGANISM] [SRA=$sraids]\'";
  my $cmd0="$APPtbl2asn $tsaopts $tsadesc $tsaqual "; # add in loop: -i $pt.fsa -f $pt.tbl -Z discrep.$pt.log
  
  my ($cmddone, $nout, $sqnout)=(0) x 10;
  my $runok= ($APPtbl2asn =~ /MISSING/) ? 0 : 1;

  #?? check for tbl2asn file set?  
  ## FIXME: may be in "submitset" subdir : use instead makename($cdnaseq,"_tsasubmit/")
   
  $sqnout= makename($cdnaseq,".sqn",'cdna|fasta|fsa|fa');   
  return(1,"",$sqnout) if( -s $sqnout ); ## all sqnout have full path ?
  if( -d $spldir ) { ## ! -s $sqnout and 
    opendir(D,$spldir); 
    my @sqnout= map{ "$spldir/$_" } grep /\.sqn/, readdir(D); 
    closedir(D);
    $sqnout= join ", ", @sqnout; $nout= @sqnout;
    return($nout,$spldir,$sqnout) if($nout>0); 
  }
  
  $sqnout="";
  #only>1c# mkdir($spldir); # dryrun? only ncpu or use for all cases?
  if($runok and $NCPU > 1) {
    my $ccount= facount($cdnaseq); # use this, not fasize; ERR if ccount<??
    if($ccount >= 50*$NCPU) {
      mkdir($spldir); # dryrun? only ncpu or use for all cases?
      ($nout,$sqnout)= tbl2asn_ncpu($NCPU,$cmd0,$ccount,$spldir,$cdnaseq,$cdnatbl);
      $cmddone=1;
      
      push @submitset, $spldir;
    	#?? push @submitset, (split /,\s*/, $sqnout); # leave in name_tsasubmit ; 
    	# move all of name_tsasubmit to submitset/ ?? YES?
      # tbl2asn output in publicset/tsa > move to submitset/ with input files
      # pig4321ew.mrna_pub_tsasubmit
      # pig4321ew.mrna_pub.pub2submit.info
      # pig4321ew.mrna_pub.pub2submit.log
      # pig4321ew.mrna_pub.sumval

    }
  } 
  
    # put into $spldir ?
  if($runok and ! $cmddone) {
    $spldir="./"; # cwd?
    my $dlog1=$cdnatbl; $dlog1 =~ s/\.\w+$//; $dlog1.=".discrep";
    my $cmd1= $cmd0 . " -i $cdnaseq -f $cdnatbl -Z $dlog1"; # want -r $spldir or not ?
    my $err= runcmd($cmd1);    
    my $sqnt= makename($cdnaseq,".sqn",'cdna|fasta|fsa|fa'); # add $spldir ??
    if( -s $sqnt ) { 
      my $f; $nout=1; $sqnout=$sqnt; 
      ## rename errorsummary so not overwrittne by others
      my $esum0= dirname($sqnout)."/errorsummary.val";
      my $esumf= makename($sqnout,".sumval"); # or cdnaseq.errorsummary.val ?
      rename($esum0, $esumf); 
      push @submitset, $sqnout, $dlog1, $esumf; 
      for my $suf (qw(val gbf fixedproducts pep.report)) { 
        ($f=$sqnt) =~ s/\.sqn/.$suf/; push @submitset,$f if(-f $f);
      }
      # t2a outputs include discrep, xxx.sqn, xxx.val, errorsum.val ..
      # FIXME: submitset tbl2asn: name.val,.sqn,.discrep,.fixedproducts errorsummary.val
    }
  }
  
  return($nout,$spldir,$sqnout); # list all?
}

## tbl_extract id pattern is '>\w+\s+(\S+)' instead of fasta '>(\S+)'
sub tbl_extract # from faextract in cdna_evigenesub.pm
{
  my ($fa,$newfa,$faids)=@_; 
  my ($n,$hout,$hids,%ids)=(0);  
  my ($ok,$hin)= openRead($fa); return 0 unless($ok);
  %ids= %$faids; # hash input only
  # rename($newfa,"$newfa.old") if(-s $newfa);  
  $ok= open($hout,'>',$newfa);   
  if($ok) { $ok=0; 
    while(<$hin>) { if(/^>/){ my($id)= (m/>\w+\s+(\S+)/)?$1:0; $ok=$ids{$id}||0; } 
    print $hout $_ if($ok); } }
  close($hout); close($hin); 
  return $newfa;  
}

sub tbl2asn_ncpu
{
  my($npart,$cmd0,$ccount,$spldir,$cdnaseq,$cdnatbl)=@_;
  
  # my $ccount= facount($cdnaseq); # use this, not fasize; ERR if ccount<??
  my $splcount= int(0.99 + $ccount/$npart);
  
  ## NOTE: fasplit adds spldir to set path
  #* FIXME diff num entries in fsa, tbl, pep, due to missing CDS info
  #* fsa *should* have most entries (all pubset.mrna), with maybe fewer in tbl .. 
  # * should repair tbl output to produce >entry for each mrna/fsa
  # BUT should change here to match tblset by IDs in fsa split set,
  # ie. use idh=faidlist(cdna1) for tbl same as pep
  
  my @splset= fasplitcount( $cdnaseq, $spldir, $npart, $splcount,"fsa"); 
  my @tblset=();
  my @pepset=(); 

  use constant SPLIT_TBLBYID => 1;
  unless(SPLIT_TBLBYID) {
  @tblset= fasplitcount( $cdnatbl, $spldir, $npart, $splcount,"tbl");  # need to split into same parts.
  }
  
  ## FIXME 2015.01, missing .pep/.pep.report for NCPU case .. need to split those also
  # * PROBLEM for pep split, out of order ids ? got 2 in wrong split set..
  # .. fsa,tbl above are same id-order as tbl is made from cdna/fsa ..
  # .. need idlist from fsa, or each part, and new fasplitordered(...) to manage those.
  
  my $cdnapeps= makename($cdnaseq,".pep"); 
  my $haspep= (-s $cdnapeps);
  #  @pepset= fasplitcount( $cdnapeps, $spldir, $npart, $splcount,"pep"); 
  
  # split by idlist..
  my $ns= @splset;
  for(my $ip=0; $ip < $ns; $ip++) {
    my $cdna1= $splset[$ip];
    my $idh = faidlist($cdna1,{},"ashash");
    my $pep1= makename($cdna1,".pep"); 
    my $tbl1= makename($cdna1,".tbl");  # upd add tbl

    if(SPLIT_TBLBYID) {
    $tbl1= tbl_extract( $cdnatbl, $tbl1, $idh, ); # fix for >Features $ID ... not fasta >ID
    push @tblset, $tbl1;
    }
    
    if($haspep) {
    # Ah, damn diff ids for .fsa > .pep Id000t1 > gnl|Evigene|Id000p1  
    # $protid= $GDB_PREFIX.$protid if($GDB_PREFIX); #  protein_id  gnl|CacaoGD|Thecc1EG016762p1
    for my $m (keys %$idh) { my $n=$m; $n=~s/t(\d+)$/p$1/;
      $n=$GDB_PREFIX.$n if($GDB_PREFIX); $idh->{$n}=1; 
    }
    $pep1= faextract( $cdnapeps, $pep1, $idh, );  
    push @pepset, $pep1;
    }
  } 
  
  my $npartgot= @splset;
  my $icpu= 0;  
  # NO: chdir($spldir); # bad!, need files in starting path ... use -i spldir/in -r spldir/
  for(my $ip=0; $ip< $npartgot; $ip++) {
    my $cdna1= $splset[$ip];
    my $tbl1 = $tblset[$ip];
    ## dont need to add pep1, tbl2asn looks for  $pep1 = $pepset[$ip];
    (my $dlog1=$tbl1) =~ s/\.\w+$/.discrep/;
    my $cmd1= $cmd0 . " -i $cdna1 -f $tbl1 -Z $dlog1"; 
    my $pid= forkcmd($cmd1);    
    if(++$icpu > $npartgot) { while (wait() != -1) { }; $icpu= 0; }
  }
  while (wait() != -1) { };
  
  ## One part problem, errorsummary.val parts overwrite .. rebuild from $name.val ? or
  # my @valout= # my @gbfout= 
  # my @sqnout= map{ my $f=$_; $f=~s/\.fsa/.sqn/; $f; } @splset;
  
  ## deal with $name.split*.fixedproducts also?
  ## deal with $name.split*.discrep also?
  
  ## remake errorsummary.val from parts  .. maybe for 1cpu also? count ERROR: and loggit ?
  my %et;
  my $esumf= makename($cdnaseq,".sumval"); # or cdnaseq.errorsummary.val ?
  foreach my $sf (@splset) { 
    (my $valf=$sf) =~ s/\.fsa/.val/;
    if(open(VF,$valf)) { 
      while(<VF>) { my($err,$typ)= m/^(\S+)[\w\s]+.([\w\.]+)/; $et{$err}{$typ}++ if($typ); } close(VF); 
    }
  }
  open(my $esumh,'>',$esumf);
  foreach my $e (sort keys %et) { my @t=sort keys %{$et{$e}}; 
    foreach my $t (@t) { my $c=$et{$e}{$t}; printf $esumh "%6d %s %s\n",$c,$e,$t; } 
  } close($esumh); 
  push @submitset, $esumf; 

  ## check existance -s sqn     # or readdir(D) as above  
  my @sqnout= grep { -s $_ } map{ my $f=$_; $f=~s/\.fsa/.sqn/; $f; } @splset; 
  my $sqnout= join ", ", @sqnout;
  my $nsqn= @sqnout;
  return($nsqn,$sqnout); # list all?
}





1;

__END__
