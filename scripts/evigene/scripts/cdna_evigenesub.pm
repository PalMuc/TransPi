# cdna_evigenesub.pm

# package cdna_evigenesub;
package main;

use strict;
use warnings;
use FindBin;
use File::Basename qw(basename dirname fileparse);

use constant TRAA2CDS_2018 => 1; # prot/traa2cds.pl usage updated for utrorf handling

#? require Exporter;
# our @ISA = qw (Exporter);
# our @EXPORT = qw (xxxx);

use vars qw ( $EVIGENES $EGAPP $GDB_PREFIX $DEBUG $dryrun
  $MIN_NAMEIDENT $MIN_IDLIKE $EVGLOGH 
  %genenames %genedbxref %genenamepct %namedgenes %cddnames 
  %pubids %pubidinfo $APPtraa2cds
  $AAQUALF $AAQUALH $BAD_GAPS
  );

## add globals for getAaQual here? should be in main caller?
use constant { kAAQUAL_MAX => 3, kAAQUAL_MIN => -3, kAAQUAL_NONE => 0, }; #  aaqualscore() range
our ($AAQUALF,$AAQUALH) = ("",undef,undef);
#NOT YET# $AASIZEH 
#  $AAQUALH->{$id}="$alen,$pctcds,$acv,$aqual1";  $acv == numeric score of aqual1
#  our %AAQUALS = (); %AASIZES = (); # global hash, keep in sync w/ file name

our $DEBUG=0;
our $dryrun=0; ## $DRYRUN ?
our $EVIGENES="$FindBin::Bin"; #??
our $EGAPP='mrna2tsa'; # FIXME?
our $EGLOG='egr';
our $GDB_PREFIX='gnl|Evigene|';  #see below;  use IDPREFIX? No, ID has
our $APPtraa2cds= undef; #findevigeneapp("prot/traa2cds.pl"); # move to cdna_evigenesub for get_mRNA

## name stuff 
our $MIN_NAMEIDENT = 35;  # min similar% for protein evid align/equivalence, for naming; JCVI uses 35%  
our $MIN_IDLIKE   = 15;  # low for now; 15..20 seems right;add Note with loqualname
our $BAD_GAPS= 25;  # % gaps in AA

use constant { LOG_NOTE => 0, LOG_WARN => 1, LOG_DIE => -1, LOG_DEBUG => 2, };
our $EVGLOGH= undef; # renamed EVGLOGH from logh; package local? now exported

sub loggit{ 
	# my $dowarn=shift; my $s= join(' ',@_); # dang warn @_ empty join for 1st call here, from where ??
	my($dowarn,@msg)= @_; return unless($dowarn or @msg);
	return if($dowarn== LOG_DEBUG and not $DEBUG);
  #my $s= join(' ',@msg); ## Use of uninitialized value $msg[3] in join
  my $s= join ' ', map{ defined($_) ? $_ : '.' } @msg;
  chomp($s); $s="FATAL $s" if($dowarn == LOG_DIE);
  if($EVGLOGH){ print $EVGLOGH "#$EGLOG: $s\n"; } elsif($dowarn>0||$DEBUG){ warn "#$EGLOG: $s\n"; }
  if($dowarn == LOG_DIE) { die "#$EGLOG: $s\n" ; }
}

sub openloggit {
  my($logfile,$trname)= @_;
  if(not $logfile and defined $logfile) { # use output name
    $logfile= $trname || $EGLOG;
    $logfile= makename($logfile,".$EGAPP.log");  # need program suffix??
  }
  if($logfile) { 
    open($EVGLOGH, '>>', $logfile) or die $logfile; 
    ## EVGLOGH should have immediate write, $| = 1 ??, like STDERR
    my $lastio = select(STDOUT);
    select($EVGLOGH); $| = 1;
    select($lastio); 
  } 
}

sub openRead { # add to cdna_evigenesub.pm
  my($fna)= @_; my($ok,$hin)= (0,undef);
  $ok= ($fna =~ /\.gz$/) ? open($hin,"gunzip -c $fna|") : open($hin,$fna);  
  loggit(1,"ERR: openRead $fna") unless($ok);
  return ($ok,$hin);
}

## note these are in cdna_protein also; need more package local privacy.
sub _min1 { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max1 { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }



sub parse_genenames
{
  my($genenames, $noEdits)= @_;
  my($ngot,$nin)=(0,0);
  $noEdits ||=0; 
  # returns in globals: (%genenames,%genenamepct,%genedbxref,%namedgenes,%cddnames) 
  %genenames=%genenamepct=%genedbxref=%namedgenes=%cddnames=();
  return($ngot,$nin) unless($genenames and -f $genenames);
  
  ## FIXME2: ** use uniq names vs ERR: too short to keep valid tr, e.g. 
  #er2g: ERR: too short:183 #LitvaEG0018688t4     oid=litovavel2k35Loc15824t1     len=183 name=CDD: clpS, ATP-dependent Clp protease ada..-like    
  # grep  'CDD: clpS,' *.names = 1 only = litovavel2k35Loc15824t1
  
  # FIXME: need better reader; 2+ rows/id; pick best .. format may change..
  # names.tab ==  id, name, pctalign, refid, repid  : now
  #  trid1  C-ets-2 protein, putative       89%,103/116,197 RefID:UniRef50_E0VFI2   RepID:E0VFI2_PEDHC
  #  trid2  DBH-like monooxygenase protein 1        73%,445/613,516 RefID:UniRef50_Q6UVY6   RepID:MOXD1_HUMAN
  
  my($ok,$inh)= openRead($genenames);
  unless($ok) { loggit(1,"ERR: parse_genenames reading $genenames"); return; }
  
  while(<$inh>) { 
    next unless(/^\w/ and /\t/);
    chomp; $nin++;
    my($id,$name,$pctalign,$refid,$repid)=split"\t"; # may have only id, name

    my $xtra=""; 
    ($name,$xtra)=split";",$name,2 unless($noEdits); #??? xtra may be valid some have ;
    $name =~ s/\s+$//;

    ## BUG in data: missing pctalign fields ; dont know why.
    ## output of evigene/scripts/prot/namegenes.pl 
    ## whitefly1vel5k45Loc9888t2       Synaptojanin-1-like protein     RefID:UniRef50_B4E1Z3   UniProt:B4E1Z3_HUMAN
    if($refid and not defined $repid and $pctalign and $pctalign =~ /^\D/) {
    	$repid=$refid; $refid=$pctalign; $pctalign="";
    }
    $pctalign||=""; $refid||=""; $repid||=""; # missing ok
    
    # FIXME: 2 names/id maybe: CDD: xxx and gene xxx; keep both in ann.txt ? and pctalign?
    ## pctalign == 100%,450/450,446 : pct,naln/nref,ntrg
    ## old geneinfo:
    #  ($td,$gd,$mapqual,$aaqual,$trlen,$cdsor,$cdsoff,$lotag,$nam1,$dbxref)= @$rinfo;
    ## usage below in putseq
    # my %tblinfo= (pubid => $pubid, oid => $oid,  protid => $protid, locustag => 0,
    #  aaqual => "na", trlen => 0, cdsoff => "", cdsor => 1, 
    #  name => $genenames{$oid}||"", dbxref => $genedbxref{$oid}|| "na"); 

    my($pcta)= $pctalign=~m/(\d+)/; # warn Argument "87%,80/92,86" isn't numeric 
    my $haspct= ($pctalign =~/^\d/ and $pcta > 0)?1:0; # zero is missing value
    if(!$noEdits and $haspct and $pcta < $MIN_NAMEIDENT) { # use MIN_IDLIKE : name-like ? got some '0%' align
      ## bad: Uncharacterized protein-like ; Nuclease HARBI1, putative-like; Protein-like
      ## *should* leave these name changes to nameclean()
      if($pcta >= $MIN_IDLIKE) { 
       unless($name =~ /\blike|^Uncharacterized/) {
        $name =~ s/, putative//; 
        unless( $name =~ s/\s+protein$/-like protein/ ) { $name .= '-like'; } ## fixme: 'xxxx protein -like'
        }
      } else { next; } ## caller should decide? should we preserve for ann.txt table ? as Unchar ?
    }

        # DBXREF_RECODE fix for ncbi's  dbx: restrictions : evgmrna2tsa2.pl:putTblFsa()
        # FIXME999: more problems w/ gene.names table having odd/local DBprefix:ID
        #   .. fix where? should have here list of valid NCBI/ISxxx db prefixes. from where?
    
    ## fixme: CDD:206692,cd04107,RefID:UniRef50_Q9NX57,UniProt:RAB20_HUMAN, 
    ## drop  RefID:; drop? cd04107
    $refid =~ s/^RefID://;  $repid =~ s/^RepID://;  ## RepID: also 
    ## ?? try here add right DbPrefix: ? Uniprot/Uniref easy, others a mess.
    map { if(/:/) { } # asis
    	elsif(/^UniRef/i or /^[A-Z0-9]+_[A-Z][A-Z]+$/) { $_="TrEMBL:$_"; } 
    	elsif(/^ENS\w+\d\d+$/) { $_="ENSEMBL:$_"; }
      } ($refid,$repid);
    
    #old# $genedbxref{$id} .= "$refid," if($refid);
    #old# $genedbxref{$id} .= "$repid," if($repid and not $refid =~ /^CDD:/);
    
    $namedgenes{$name} .= "$id,"; #? if($pctalign >= $MIN_NAMEIDENT); # for uniq name retention
    
    ## FIXME: keep CDD names for .ann.txt, maybe .tbl submit as 2nd note
    $cddnames{$id}= $name if($name =~ /CDD:/ and not $cddnames{$id});
    
    ## FIXME: 1st in dont replace.. not just for CDD: ?
    unless($genenames{$id}) { #was and $name =~ /CDD:/ 
      $pctalign ||= 0; $refid ||= 0;
      $genedbxref{$id} .= "$refid," if($refid);
      $genedbxref{$id} .= "$repid," if($repid and not $refid =~ /^CDD:/);
      $genenames{$id}= $name;  $ngot++;
      $genenamepct{$id}= $pctalign;
       # repid Yes or No? this is by default RefID=UniRef50_xxxx and RepID=UniProt:xxxx_HUMAN
    }
  } close($inh);
  
  return($ngot,$nin);
}

sub parse_evgheader
{
  my($oid,$hdr,$trlen,$seqoid)= @_;
  $seqoid ||= $oid; # for split/dup gff gene IDs: Id_C1,2 Id_G2,3,.. names have seqoid
  $trlen  ||= 0;
    ## this becomes param or not?
    ## oid maybe pubid, check hdr for others
  my $pubid= $pubids{$oid} || $oid; # is it ERR if no pubid{oid} ?

  my $protid= $pubid;
  $protid=~s/t(\d+)$/p$1/; # genbank requires diff protid from mrnaid
  $protid= $GDB_PREFIX.$protid if($GDB_PREFIX);
  #  protein_id      gnl|CacaoGD|Thecc1EG016762p1

#     if(my $lotagpre= $settings{'LOCUSTAG'}) {
#       my($pubidnum)= $pubid =~ m/$IDPREFIX(\d+)/;
#       $locustag= "$lotagpre$pubidnum" if($pubidnum);
#     }

  my %tblinfo= (pubid => $pubid, oid => $oid,  protid => $protid, locustag => 0,
      aaqual => "", trlen => $trlen, cdsoff => "", cdsor => 1, 
      name => "", namepct => 0, dbxref => '', cdd => ''); ## not yet "na" 
  # tblinfo.Specials:  Selcstop
  # tblinfo.addkeys:   locus/location/maploc, mapqual
  
  use constant CHECK_NAMEOIDS => 1;
  my $nameoid= $seqoid;
  if(CHECK_NAMEOIDS and not $genenames{$nameoid} ) {
    if($hdr =~ m/\boid=([^\s;]+)/) { my $alloids=$1;
      for my $d (split",",$alloids){ if( $genenames{$d} ) { $nameoid=$d; last; }  }
    }
  }  

  if( $genenames{$nameoid} ) {
    $tblinfo{'name'}= $genenames{$nameoid};
    $tblinfo{'namepct'}=  $genenamepct{$nameoid} || 0;
    $tblinfo{'nameref'}= $tblinfo{'dbxref'}=  $genedbxref{$nameoid}; # ||"na" # should this be 'nameref' instead?
    $tblinfo{'cdd'}=  $cddnames{$nameoid}; # ||"na"
  }
        
  my($cdsb,$cdse,$aafull)=(0,0,0);
  if($hdr =~ m/\b(?:offs|cdsoff)=([\d-]+)/) { my $cdsoff=$1; $tblinfo{'cdsoff'}= $cdsoff; 
    ($cdsb,$cdse)= split/[-]/,$cdsoff;  } # do in putseq
  if($hdr =~ m/\b(?:aaqual|aalen)=([^\s;]+)/) { my $aq=$1; $tblinfo{'aaqual'}= $aq; 
    ($aafull)= $aq =~ m/(complete|partial\w*)/; }
  if($hdr =~ m/\bclen=(\d+)/) { my $ln=$1; unless($trlen){ $trlen=$ln; $tblinfo{'trlen'}= $trlen; } } #  and not $trlen ; skip? 
  if($hdr =~ m/\bstrand=([^\s;]+)/) { $tblinfo{'cdsor'}= $1; } # expect all '+' from traa2mrna
  if($hdr =~ m/\bSelcstop=([^\s;]+)/) { $tblinfo{'Selcstop'}= $1; }  # Selcstop update 14.12.30
  if($hdr =~ m/\boid=([^\s;]+)/) { my $od=$1; 
    if($od eq $oid){$od="";} elsif($od=~/$oid/){} elsif($oid ne $pubid){$od="$oid,$od";} #?
    $tblinfo{'oid'}= $od if($od); } ## oid param maybe pubid, check others
  #^ add (?:Target|trg)= as oid alternate of gff
  if($hdr =~ m/\b(?:Target|trg)=([^\s;]+)/) { my $trg=$1; my $od=$tblinfo{'oid'}||"";
    unless($od=~/$trg/){ $od.="," if($od); $od.= $trg; $tblinfo{'oid'}= $od; } }
  if($hdr =~ m/\b(?:[Nn]ame|[Pp]roduct|[Pp]roduct[Nn]ame)=([^=\n;]+)/) { # Product_Name= ?
     my $na=$1; $tblinfo{'name'}= $na unless($tblinfo{'name'}); }
  if($hdr =~ m/\b[Nn]amepct=([^\s;]+)/) { my $nap=$1; 
    $tblinfo{'namepct'}= $nap unless($tblinfo{'namepct'}); }
  elsif($hdr =~ m/\b[Nn]amealn=([^\s;]+)/) { my $naa=$1;
    ## namepct =>  namealn=58p,17578/30278,17622; .. keep align stats?
    my($nap)= $naa=~m/^(\d+)/?$1:0;  my($nal)= $naa=~m/,(\d+)/?$1:0;
    $tblinfo{'namepct'}= $nap if($nap and not $tblinfo{'namepct'});
    $tblinfo{'namealn'}= $nal if($nal);
    }
  if($hdr =~ m/\bevgclass=([^\s;,]+)/) { $tblinfo{'evgclass'}= $1; } # upd1803
  if($hdr =~ m/\btype=(\w+)/) { $tblinfo{'seqtype'}= $1; } # upd1803; type| seqtype| moltype ?
    
    #?? merge hdr dbxref and genenames dbxref?  
    ## problems w/ new/old dbxref near same, drop old, treat genenames as accurate
    # .. check dup uniprot _species now..
  if($hdr =~ m/\b(?:[Dd]bxref|db_xref)=([^\s;]+)/) { my $dbx=$1; 
    if(my $tdx=$tblinfo{'dbxref'}) { 
      for my $dx (split",",$dbx) { 
        my($d,$x)=split":",$dx; 
        if($x and not $tdx=~/$x/) {
          my $ok=1;
          if($d=~/SwissProt|TrEMBL|UniProt/) { 
            my($sp)= $x=~m/(_\w+)$/; $ok=0 if($sp and $tdx =~ /$sp/);
          }
        $tdx.=",$dx" if($ok);  
        }
      }
      $tblinfo{'dbxref'}= $tdx;
    } else { $tblinfo{'dbxref'}= $dbx; }
  }

  map{ $tblinfo{$_} ||= "na" } qw(aaqual dbxref cdd);
  
  return \%tblinfo;
}


sub getOkFileset
{
  my($okpath,$SUFFIX,$okfiles,$trname)= @_;
  
  #FIXME.1712: added trname opt, may have many projects in okayset/ 
  # see also getmRNA() .. too complex to use here, elsewhere
  #    ($trf)= grep /$trname\.(mrna$|mrna.gz$)/, (@pubd, @okd); # drop okd here?

  my($oktr,$alttr)=("","");
  $SUFFIX='tr|cdna|fasta|fna' unless($SUFFIX); # is this enough? NOT .mrna  
  my @okd=(); 
  if($okfiles and ref($okfiles) and @$okfiles > 0) { @okd= @$okfiles; }
  elsif(-d $okpath) { opendir(D,$okpath); @okd= map{ chomp; "$okpath/$_" } readdir(D);  closedir(D); }
  if($trname) { @okd= grep /$trname/, @okd; } # fix for different subprojects in okayset/ ..
      
  #FIXME: '($SUFFIX)$' so dont match trxxx.logfile
  ($oktr) = grep /.okay\.($SUFFIX)/, @okd;  
  ($alttr)= grep /.okalt\.($SUFFIX)/, @okd; 
  return($oktr,$alttr,\@okd); ## change to? (\@okd,$oktr,$alttr)  
}

sub getFileset
{
  my($okpath,$SUFFIX,$okfiles,$trname)= @_;
  $SUFFIX='tr|cdna|fasta|fna' unless($SUFFIX); # is this enough? NOT .mrna  
  #FIXME.1712: added trname opt, may have many projects in okayset/ 
  my @okd=(); 
  if($okfiles and ref($okfiles) and @$okfiles > 0) { @okd= @$okfiles; }
  elsif(-d $okpath) { opendir(D,$okpath); @okd= map{ chomp; "$okpath/$_" } readdir(D);  closedir(D); }
  if($trname) { @okd= grep /$trname/, @okd; } # fix for different subprojects in okayset/ ..
  #FIXME: '($SUFFIX)$' so dont match trxxx.logfile
  my @files = grep /\.($SUFFIX)/, @okd; my $nok=@okd;
  #ok here# warn "#getFileset($okpath,$SUFFIX,$okfiles)= dir:$nok, suf:@files \n" if($DEBUG); 
  return(\@okd, @files); 
}

sub tidyupFileset { 
	my($tod,@td)= @_;  
	return unless($tod and @td);
	my @tdlist;
	mkdir($tod) unless(-d $tod);
	foreach my $fn (@td) { 
  	if($fn and not ($fn =~ m,$tod/,)) { 
    #   if(-f $fn ) {  # old only -f
    #     (my $tfn=$fn) =~ s,^[\w\.]+/,,;   ## ASSUMES tod is subdir in current path !!
    #     rename($fn,"$tod/$tfn");  push @tdlist, "$tod/$tfn";
    #   }
      (my $tfn=$fn) =~ s,^[\w\.]+/(.+),$1,;   ## ASSUMES tod is subdir in current path !!
      my $todfn= "$tod/$tfn";
      if(-f $fn ) { 
        rename($fn, $todfn); push @tdlist, $todfn;
      } elsif( -d $fn ) { # move subdirs ok this way?
        rename($fn, $todfn); push @tdlist, $todfn;
      }
  	}
  } 
	my $n=@tdlist; my $n1= _min1($n-1,4); 
	loggit(0,"tidy: n=",$n, @tdlist[0..$n1]); 
}


=item AaQual evigene attribute
  
  AaQual is a transcript attribute extensively used by Evigene.
  Value is a tuple: "#aa-length,#coding-percent,Completeness"
  where aa-length is count of aa residues, including gaps (usually)
  coding-percent is %(CDS-length/mRNA-length)
  Completeness is controlled vocabulary: complete|partial3|partial5|partial 
    (partial=missing 5' and 3' ends, partial5=missing 5', ..)
    with other appended: -utrbad|-utrpoor|-gapbad|..
  It is calculated from proteins of mRNA transcripts following ORF translation.
  
  Evigene ORF sequences (.aa and .cds) and size table (.aa.qual) have this and
  other  ORF translation values, 
    offs=CDS-offset (b-e) in mRNA/cDNA, (e-b) for revcomp
    strand=+|- in cDNA/mRNA,
    clen=cDNA/mRNA length, 
    aalen=AaQual tuple or simple aa-length,

=item AaQual score
   
   This is integer value of Completeness vocabulary, with "complete" only as highest value.
   "partial", "utr" and "gap" attributes reduce score.
   Current range is -3..+3 (kAAQUAL_MIN..kAAQUAL_MAX) 
   
   A single numeric comparison of transcript Aa would include aa-size, coding% and completeness,
   for instance for selecting or sorting transcripts / proteins.
   Perhaps aascore = aa-length * codingpct/100 * (aaqual - kAAQUAL_MIN) / (kAAQUAL_MAX - kAAQUAL_MIN)
   
=cut

#above# use constant { kAAQUAL_MAX => 3, kAAQUAL_MIN => -3, kAAQUAL_NONE => 0, };
sub aaqualscore
{
  my($mqual)= @_;  $mqual ||="missing";
  my $mqv= kAAQUAL_NONE; 
  if($mqual =~ /complete/) { $mqv = kAAQUAL_MAX; } 
  elsif($mqual =~ /partial[35]/) { $mqv = kAAQUAL_MAX - 1; }
  if($mqual =~ /utrbad/) { $mqv = ($mqv >= kAAQUAL_MAX-1) ? kAAQUAL_NONE - 1 : $mqv - 2; } 
  elsif($mqual =~ /utrpoor/) { $mqv -= 1; }
  if($mqual =~ /gapbad/) { $mqv -= 1; } # or -2?
  return $mqv; # range is now -3..+3, was -2 .. +2
}


## replacing evigene/scripts/prot/aaqual.sh
sub makeAaQual {
  my($aaseq,$aaqualSetHash)= @_;
  my $doff=1; # $ENV{doff}; 
  my $doid=1; # no oids here? drop this? option?  $ENV{doid};
  our $makeAaQual_SETHASH= $aaqualSetHash||0;
  
  # minor update bug: this always makes new, should instead reuse old aasize..
  my $aasize= makename($aaseq,".aa.qual"); 
  if( -s $aasize ) {
    my $aqhash= getAaQual($aasize) if($makeAaQual_SETHASH); # yes or not?
    return($aasize);
  }
  
  my ($ok,$inh)= openRead($aaseq);  
  $ok= open(AAQ,'>',$aasize) if($ok);
  return unless($ok);
  
  # use constant  makeAaQual_SETHASH => 1;
  my($id,$aat,$aag,$al,$cl,$naa)= (0) x 9;
  
  sub puta { 
    my($id,$aat,$aag,$al,$cl)= @_; our $makeAaQual_SETHASH;
    $al=($aat+$aag).",na" if($al eq "na"); 
    print AAQ join("\t",$id,$aat,$aag,$al,$cl)."\n"; 
    ## option here: set global hashes
    if($makeAaQual_SETHASH) {
       my($aww,$pctcds,$aqual1)=split",",$al; 
       $pctcds =~ s/\%//;      
       $aqual1 .= "-gapbad" if($aag>0 and (100*$aag/($aat+$aag) > $BAD_GAPS)); # add qual flag for gap/(alen+gap) > MAXGAP:  
       my $acv= aaqualscore($aqual1);
       #NOT YET# $AASIZEH->{$id}=$aat; 
       $AAQUALH->{$id}="$aat,$pctcds,$acv,$aqual1"; 
    }
    return 1;
  }

  if($makeAaQual_SETHASH) {  $AAQUALF=$aasize; $AAQUALH={}; } # $AASIZEH={}; 
  while(<$inh>) {
    if(/^>(\S+)/) { my $td=$1; 
      $naa+= puta($id,$aat,$aag,$al,$cl) if($id);  $id=$td; $aat=$aag=0; 
      ($al)=m/aalen=([^;\s]+)/; $al||="na";
      ($cl)=m/clen=(\d+)/; $cl||=0; 
      if($doff){ 
        my($or)=m/strand=(.)/; $or||="."; 
        my($ofs)=m/offs=([\d-]+)/; $cl.= ($ofs)?"\t$ofs:$or":"\t0";  
        } 
      if($doid){ my $oid;
        unless(($oid)=m/oid=([^;\s]+)/) { ($oid)=m/gene[=:]([^;\s]+)/; } 
        $oid||="noid"; $cl.="\t$oid"; 
      } 
    } else { s/\*$//; $aat += tr/A-WYZa-wyz/A-WYZa-wyz/; $aag += tr/Xx\*/Xx\*/; }
  } 
  $naa+= puta($id,$aat,$aag,$al,$cl);  
  close(AAQ); close($inh);
  loggit(0, "makeAaQual: naa=$naa IN $aasize\n") if($DEBUG);
  ## optionally set global hashes: $AAQUALF,$AAQUALH,$AASIZEH
  return($aasize);
}


# getAaQual returns $aaqual->{$id}="$alen,$pctcds,$acv,$aqual1"; see also asmdupfilter:readSizes()
sub getAaQual {
  my($aaqualf)= @_;
  my($naa,$ntr,$nerr,$ok,$inh)=(0) x 10;
  
  # our ($AAQUALF,$AAQUALH,$AASIZEH) = ("",undef,undef);
  if($AAQUALF and $AAQUALF eq $aaqualf and ref($AAQUALH)) { 
    $naa= scalar(keys %$AAQUALH);
     if($DEBUG) { my($id1)= (keys %$AAQUALH)[0];
      loggit(0, "getAaQual: naa=$naa in $AAQUALF, val1 $id1=",$AAQUALH->{$id1});
      }
    return (wantarray) ? ($AAQUALH,$naa) : $AAQUALH; # %aasize also? NOT YET ,$AASIZEH
  }
  
  my %aaqual =(); # use global hash, keep in sync w/ file name
  #NOT YET# my %aasize =(); # use global hash, keep in sync w/ file name
  if($aaqualf) {  
    ## fix for aacount gaps: id,size,gaps : NOT NOW, aa.qual: id,size-gaps,gaps,..
    ## drop faCount? require aa.qual here?

    ($ok,$inh)= openRead($aaqualf); # $ok= open($inh,$aaqualf);
    unless($ok) { loggit(1,"ERR: getAaQual $aaqualf"); return 0; } # die "FAIL: read $aaqualf ..."; 

    while(<$inh>) { 
      next if(/^\W/ or /^total/); 
      my($id,$alen,$gap,$aqual,$tlen)=split; # aa.qual cols; gap is removed from alen
      unless($alen =~ /^\d/) { $nerr++; next; }
      
      if($aqual) { 
        $aqual .= "-gapbad" if($gap>0 and (100*$gap/($alen+$gap) > $BAD_GAPS)); # add qual flag for gap/(alen+gap) > MAXGAP:  
        my $acv= aaqualscore($aqual);
        my($aww,$pctcds,$aqual1)=split",",$aqual;  
        $pctcds =~ s/\%//;  
        $aaqual{$id}="$alen,$pctcds,$acv,$aqual1"; 
        # $aaqual{$id}= $aqual; # change this?
        #?? $aqv{$id}= $alen * $pctcds * $acv; # want something like this, better, for 1 num score
        #? want# if($tlen =~ /^\d/) { $trsize{$id}= $tlen; $ntr++; }
      } else {
        my($pctcds,$aqual1)=(99,"na"); my $acv= aaqualscore($aqual1);
        $aaqual{$id}="$alen,$pctcds,$acv,$aqual1"; 
      }
      #NOT YET# $aasize{$id}=$alen; #? need both aaqual and aasize hashes?
      $naa++; 
    } close($inh); 
    
    if($naa) { $AAQUALF=$aaqualf; $AAQUALH= \%aaqual; }
    ## dont need yet:  $AASIZEH= \%aasize;
  }
  
  if($DEBUG) { my($id1)= (keys %$AAQUALH)[0];
    loggit(0, "getAaQual: naa=$naa in $AAQUALF, val1 $id1=",$AAQUALH->{$id1});
    }
  ## AAQUALH isnt reset unless $naa..
  return (wantarray) ? ($AAQUALH,$naa) : $AAQUALH; # which ?
  #? return (wantarray)? (\%aaqual,$naa) : \%aaqual; # %aasize also? not yet
}



sub getmRNA   ## move to cdnasubs.pm ?
{
  my($okpath,$trname,$pubdir,$ADDutrorf)= @_;
  my($cdnaseq)=(""); # == mrna (oriented), not cdna.parts.tr
  
	use constant ALSOMAKE_AACDS => 1;
  $ADDutrorf= 1; # what? always check for okayset/*.utrorf.mrna ?
  
  #? FIXME? suffix .tr may not be used: .cdna? .fa? .fasta? ...
  my $TRSUFFIX='tr|cdna|fasta|fna'; # is this enough? NOT .mrna|cds|aa,  
   ## ?? add .fa ?? BUT conflict in now publicset/*.{mrna,cds,aa}_pub.fa **

  #old1712@my($oktr,$alttr,$okd)= getOkFileset($okpath,$TRSUFFIX);
  my($oktr,$alttr,$okd)= getOkFileset($okpath,$TRSUFFIX,undef,$trname);
  my @okd= @$okd;
  my @pubd=();
  if($pubdir and -d $pubdir) { my($pubd)= getFileset($pubdir,$TRSUFFIX); @pubd= @$pubd; }
  #upd1712
  @okd= grep /$trname\./, @okd;
 
  ## another bad call: #egr: get_evgtrset= publicset/locust1all5asm.p4.mrna . locust1all5asm
  ##   instead of publicset/locust1all5asm.mrna .. my test version mistake..
  ## Ugh! bad call: #er2g: get_evgtrset= ./okayset/pogonus1all3.mrna0.ann.txt . pogonus1all3
  ## need grep /\.mrna$|\.mrna.gz$|.mrna.fa$/ ??? 
  
  # UPD1807: allow other okayset/okXXX parts? or add to okay/okalt files?
  
  my ($trf);
  ($trf)= grep /$trname\.(mrna$|mrna.gz$)/, (@pubd, @okd); # drop okd here?
  # unless($trf) { ($trf)= grep /\.mrna$|\.mrna.gz$/, (@pubd); } # , @okd want this or not?
  ## FIXME? add \.mrna_pub.fa[.gz]
  unless($trf) { ($trf)= grep /\.mrna_pub\.fa|\.mrna$|\.mrna.gz$/, (@pubd); } # , @okd want this or not?
  if($trf and -s $trf) { $cdnaseq= $trf; }
  else { 
    my $okall=0;
    my $cdnatmp= ($pubdir) ? $pubdir : $okpath;
    $cdnatmp .= "/$trname.mrna";  
    mkdir($pubdir) if($pubdir and not -d $pubdir);

    ## FIXME: utrorf : made okayset/*.utrorf.{mrna,aa,cds} ; merge into update_mrna_fileset() or getmRNA/okayset ??

    loggit(0,"Make mRNA $cdnatmp from okayset transcripts");
    my($okaa) = grep /.okay\.aa$|.okay\.aa.gz$/, @okd;  
    my($altaa)= grep /.okalt\.aa$|.okalt\.aa.gz$/, @okd; 
    # #.. problem here .aa.gz not .aa; should makename() do optionally?
    
  ## FIXME: fail here if missing oktr ?
		$APPtraa2cds= findevigeneapp("prot/traa2cds.pl") unless($APPtraa2cds); #  
		my $tropts="-trout "; $tropts .= "-nomiss " if($ADDutrorf);
		
=item traa2cds new usage

  .. detect version of traa2cds? traa2cds.pl -help : VERSION >= 2018
  .. need all tr inputs to resolve utrorfs:  
  prot/traa2cds.pl -mrnaout -aa $pt.okay.aa -aa $pt.okalt.aa -cdna $pt.okay.cdna -cdna $pt.okalt.cdna -out $pt.okall.mrna
  .. OR can make mrna of okay,okalt separately using one -aa, but both -cdna (only aa-valid IDs written)

=cut
 		
if(TRAA2CDS_2018) {
    $ADDutrorf=0;
    my $cmd= "$APPtraa2cds -mrnaout -out $cdnatmp";
    $cmd.=" -aa $okaa -cdna $oktr" if($oktr and -f $oktr and -f $okaa);
    $cmd.=" -aa $altaa -cdna $alttr" if($alttr and -f $alttr and -f $altaa);
    my $err= runcmd($cmd);
    $okall +=2  unless($err);
     
} else {
    if($oktr and -f $oktr and -f $okaa) { 
      my $err= runcmd("$APPtraa2cds $tropts -cdna $oktr -aa $okaa -out stdout >> $cdnatmp");
      $okall++ unless($err);
    }
    if($alttr and -f $alttr and -f $altaa) { 
      my $err= runcmd("$APPtraa2cds $tropts -cdna $alttr -aa $altaa  -out stdout >> $cdnatmp");
      $okall++ unless($err);
    }
    
    if($ADDutrorf and $okall > 0) {
 			my($okin) = grep /.utrorf.mrna$|.utrorf.mrna.gz$/, @okd;  
   		if($okin) { runcmd("cat $okin >> $cdnatmp"); loggit(0,"add $okin to $cdnatmp"); } # err check? loggit?
   		else { $ADDutrorf=0; } # dont do .aa,cds
    }
}
    
    $cdnaseq= $cdnatmp if(-s $cdnatmp);
    loggit(LOG_DIE,"FATAL: make mRNA $cdnatmp: missing files oktr=$oktr,okaa=$okaa,alttr=$alttr,altaa=$altaa,
    	from (oktr,alttr)=getOkFileset(path=$okpath,suffix=$TRSUFFIX)") unless($okall>1);
    
    # FIXmaybe: also make pubdir/.aa,.cds along with .mrna ? see hassle in update_mrna_fileset 
    # FIXME2: ok*.aa and utrorf.aa have dup entries, use utrorf if there.
    if(ALSOMAKE_AACDS and $pubdir) {
  		my($ok,$hin,$hout,$fout,$okin,$altin,$utrin,); my %ids=();
  		foreach my $suf (".aa",".cds") {
				$fout= makename($cdnatmp,$suf);  
				($okin) = grep /.okay$suf$|.okay$suf.gz$/, @okd;  
				($altin)= grep /.okalt$suf$|.okalt$suf.gz$/, @okd; 
				($utrin)= ($ADDutrorf) ? grep(/.utrorf$suf$|.utrorf$suf.gz$/, @okd) : (); 
				if($okin and $altin and not -f $fout) {
					$ok= open($hout,'>',$fout);	my(%did,$did);
					if($utrin) { ($ok,$hin)= openRead($utrin); 
						while(<$hin>){ if(/^>(\S+)/) { $did=($did{$1}++)?1:0; } print $hout $_ unless($did); } close($hin); }
					($ok,$hin)= openRead($okin);  
						while(<$hin>){ if(/^>(\S+)/) { $did=($did{$1}++)?1:0; } print $hout $_ unless($did); } close($hin);
					($ok,$hin)= openRead($altin); 
						while(<$hin>){ if(/^>(\S+)/) { $did=($did{$1}++)?1:0; } print $hout $_ unless($did); } close($hin);
					close($hout);
					map{ $ids{$_}{$suf}=$did{$_} } keys %did;
					}
				}
				
  	  ## should check for ID agreement among output .mrna, .aa, .cds  
			my $nidok=0;
			($ok,$hin)= openRead($cdnaseq); while(<$hin>){ if(/^>(\S+)/) { $ids{$1}{'.mrna'}++; } } close($hin);
			foreach my $id (sort keys %ids) {
				my @suf= sort keys %{$ids{$id}}; 
				if(@suf==3) { $nidok++; } 
				else { loggit(1,"ERR:getmRNA-misspart $id:",@suf); }
			}	
			loggit(0,"getmRNA $pubdir/mrna,aa,cds nid=",$nidok);
    }
  }
  
  return($cdnaseq);  # , $aaseq, $cdsseq    
}


## from asmrna_trimvec : generalize mrna_update_fileset ..
## see tr2aacds.pl:asmdupclass_fileset
sub update_mrna_fileset
{
  my($trpath, $inmrna, $updatename, $trimids, @trimfiles)= @_; 
  my $upstatus=0;
  # outputs now should go to pubdir, but use inmrna path.
  my($trimmrna, $trimaa, $trimcds, $trimidfile)= @trimfiles; # hash instead? for below %fset
	map{ $_||="" } ($trimmrna, $trimaa, $trimcds, $trimidfile);
	my %okids=(); ## return hash of valid ids * AND add oids if found, for cross checks
	
	my $ntrim=0;
	if($trimids and ref($trimids)) {
		$ntrim= scalar(keys %$trimids);

	} elsif($trimidfile and -f $trimidfile) { # hash of ids in trimfiles, regenerate from trimidfile
   	$trimids={};  
   	my($ok,$hin)= openRead($trimidfile); 
   	while(<$hin>) { if(/^\w+/) { 
   		my($id,$action,@info)=split; # parse more of trimidfile? expect: id, action=DROP/OKCUT/PROBLEM, trimnotes, aanewhdr 
   		$trimids->{$id}=$action; $ntrim++; } 
   	} close($hin);

	} else { # ERROR unless -f flagtrimvec
 		#below# loggit(1, "ERR update_fileset  empty trim ids= $trimidfile"); 
	}
	
## ........ **#@&@!*% this file name wrangling is a big waste of time ..........
## ........ use fewer naming choices ???  no okdir for pubdir set
  ##upd1807: work on mrna_pub.fa => aa_pub.fa, cds_pub.fa names
  
	my $flagtrimvec= makename($inmrna,'.'.$updatename); 
	my($outaa,$outcds);
	if($inmrna =~ m/mrna_\w+\./) { # publicset/mrna_(pub|cull|xxx).fa
	  ($outaa = $inmrna) =~ s/mrna_/aa_/;
	  ($outcds = $inmrna) =~ s/mrna_/cds_/;
	} else {
    $outaa = makename($inmrna,".aa");  
    $outcds= makename($inmrna,".cds"); 
  }
  
  my($pubdir);
  #badlocal#(my $mrnapath= $inmrna) =~ s,/[^/]+$,,; ## THIS MAY BE WRONG NOW .. publicset/ vs okayset/
  ## ^^ FIXME localdir ./name.mrna == name.mrna 
  (my $mrnapath= $inmrna) =~ s,[^/]+$,,; 
  unless($mrnapath) { $mrnapath="./"; } else { $mrnapath =~ s,/$,,; }
  ($pubdir)= getFileset($mrnapath);
  
  ## drop okalt here? see above getmRNA/ALSOMAKE_AACDS, expect/require pubdir/aa,cds?
  my $aapatt= basename($outaa);
  my($inaaseq) = grep /$aapatt$|$aapatt.gz$/, (@$pubdir);  
  if($inaaseq) { $outaa=$inaaseq;  } else { $inaaseq=""; } # fail? ignore miss
  
  my $cdspatt=basename($outcds);
  my($incdsseq) = grep /$cdspatt$|$cdspatt.gz$/, (@$pubdir);  
  if($incdsseq) { $outcds=$incdsseq; } else { $incdsseq=""; } # fail?
  
  ## FIXME:for -novectrim set updatename == 
  if( $updatename =~ /SKIP/ or -f $flagtrimvec) {
  	my($ok,$hin)= openRead($inmrna); 
  	while(<$hin>) { if(/^>(\S+)/) { my $oid= (m/\boid=([^,;\s]+)/)?$1:1;  $okids{$1}=$oid; } } 
  	close($hin);
  	return ( 0, $inmrna, $outaa, $outcds, \%okids);   	# what if outaa,outcds missing?
  }
  if($ntrim < 1 or ! -s $trimmrna) { # fixme for -novectrim
		loggit(1, "ERR update_fileset  empty trim files= $trimmrna, $trimidfile"); 
		return (-1, $inmrna, $outaa, $outcds, \%okids);
  	}

  my $upmrna  = makename($inmrna,".mrna_upd"); 
  my $upaaseq = makename($inmrna,".aa_upd"); # failed empty outaa from above
  my $upcdsseq= makename($inmrna,".cds_upd");  # failed empty outcds from above
  
  $outaa =~ s/.gz$//; $outcds =~ s/.gz$//; # output fupname
  (my $outmrna= $inmrna)=~ s/.gz$//;
  my %fset= (
            # $nup,$nsame,$fin,$ftrim,$fup,$fupname
    mrna => [ 0, 0, $inmrna, $trimmrna, $upmrna, $outmrna ],
    aa   => [ 0, 0, $inaaseq, $trimaa, $upaaseq, $outaa],  #?? fix here for .okay.aa + .okalt.aa ?
    cds  => [ 0, 0, $incdsseq, $trimcds, $upcdsseq, $outcds],
  );
  
  foreach my $suf (sort keys %fset) { # is inmrna always .mrna ?
    my($ok,$hin,$hup,%keptids);
    my($nup,$nsame,$fin,$ftrim,$fup,$fupname)= @{$fset{$suf}};

      # is one only allowed here? -s ftrim but no fin?
    if(-s $fin and -s $ftrim) {
      %keptids=();   
      $ok= open($hup,'>',$fup); # unless($ok)...
      ($ok,$hin)= openRead($fin);  
      $ok=0; while(<$hin>) { 
      	if(/^>(\S+)/) { my $d=$1; my $oid= (m/\boid=([^,;\s]+)/)?$1:1; 
      	$ok=($trimids->{$d})?0:1; if($ok){ $keptids{$d}=$oid; $nsame++; } } 
      	print $hup $_ if($ok); }
      close($hin);
      
    	## pull trimset/uvcut.$suf, check? for trimids{id} and/or collect above kept ids
      ($ok,$hin)= openRead($ftrim); 
      $ok=0; while(<$hin>) { 
      	if(/^>(\S+)/) { 
      	  my $d=$1; my $oid= (m/\boid=([^,;\s]+)/)?$1:1; 
      	  $ok=($trimids->{$d} and not $keptids{$d})?1:0; 
      		if($ok){ $keptids{$d}=$oid; $nup++; } }
      	print $hup $_ if($ok); } 
      close($hin);
      close($hup);
      $fset{$suf}->[0]= $nup; $fset{$suf}->[1]= $nsame;  # $fset{$suf}= [$nup,$nsame,$fin,$fin2,$ftrim,$fup,$fupname];
      $upstatus++ if($nup>0);
    	%okids= %keptids if($suf eq 'mrna');
    } else {
      ## error, warn/loggit, skip?
      loggit(1, "ERR update_fileset.$suf empty $fin or $ftrim"); 
    } 
  }
  ## end merge loop

  my (@outfiles,@tmpfiles);
  if($upstatus == 3)  { # $upstatus == 3 or > 0?
    ## rename input files to input.old, better: input.untrim .notrim? .pretrim? .old?
    ## rename newmerge files to input
    foreach my $suf (sort keys %fset) { 
      my($nup,$nsame,$fin,$ftrim,$fup,$fupname)= @{$fset{$suf}};  
      if(-s $fupname) { rename($fupname,"$fupname.untrim"); push @tmpfiles, "$fupname.untrim"; }
      rename($fup,$fupname); push @outfiles, $fupname;
      ## push @tmpfiles, $ftrim; # from @trimfiles. dont call tmpfile? 
      loggit(0, "update_fileset.$suf upd=$nup, same=$nsame, $fin + $ftrim > $fupname"); 
    }
    runcmd("touch $flagtrimvec"); ## touch flag-file that new input has uvtrim results ..
  } else {
    loggit(1, "ERR update_fileset missing status=$upstatus/3");  # list fupnames?
    foreach my $suf (sort keys %fset) { # sort: aa,cds,mrna
      my($nup,$nsame,$fin,$ftrim,$fup,$fupname)= @{$fset{$suf}};  
      push @outfiles, $fup;
      loggit(1, "ERR update_fileset.$suf upd=$nup, same=$nsame, $fin + $ftrim >$fup"); 
    }
  }
  
  return ($upstatus, \@outfiles, \@tmpfiles, \%okids); # $upaaseq, $upcdsseq,$upmrna, 
}


sub facount {
  my $fa=shift; my $n=0;
  my($ok,$hin)= openRead($fa);
  # my $ok= ($fa =~ /\.gz$/) ? open(F,"gunzip -c $fa|") : open(F,$fa); # openRead()
  if($ok) { while(<$hin>) { $n++ if(/^>/); } close($hin); }
  return $n;  
}

sub fasize {  
  my $fa=shift; my $n=0;
  my($ok,$hin)= openRead($fa);
  # my $ok= ($fa =~ /\.gz$/) ? open(F,"gunzip -c $fa|") : open(F,$fa);# openRead()
  if($ok) { while(<$hin>) { $n += (length($_) - 1) unless(/^>/); } close($hin); }
  return $n; 
}

sub faidlist 
{
  my ($fa,$faids,$options)=@_; # other opt to return hash not file.
  my ($n,$noupdate,$ashash,$hout)=(0) x 9;  
  $options||=""; 
  $noupdate=($options=~/noupdate/)?1:0;
  $ashash=($options=~/hash/)?1:0;
  if($ashash) {
    my %faids=(); $faids= \%faids;
    my ($ok,$hin)= openRead($fa); 
    if($ok) { while(<$hin>) { if(/^>(\S+)/){ my $id=$1; $n++; $faids{$id}=1; } } }
    close($hin); 
  } else {
    $faids= makename($fa,".ids") unless($faids); 
    return $faids if($noupdate and -s $faids);
    my ($ok,$hin)= openRead($fa);  return 0 unless($ok);
    rename($faids,"$faids.old") if(-s $faids); #? OR option to return faids as is if exists...
    $ok= open($hout,'>',$faids);   
    if($ok) { while(<$hin>) { if(/^>(\S+)/){ my $id=$1; $n++; print $hout "$id\n"; } } }
    close($hout); close($hin); 
  }
  return $faids;  
}

sub faheaders {  # use with parse_evgheader($id,$hdr)..
  my($fa,$headhash)=@_; 
  my $n=0; $headhash={} unless(ref $headhash); 
  my($ok,$hin)= openRead($fa);
  if($ok){ while(<$hin>){ if(/^>(\S+)/){ my $id=$1; $headhash->{$id}= $_; $n++; }} close($hin); }
  return ($headhash,$n); 
}


sub faextract #  in cdna_evigenesub.pm
{
  my ($fa,$newfa,$faids,$noupdate)=@_; 
  my ($n,$hout,$hids,%ids)=(0);  $noupdate||=0; 
  $newfa= makename($fa,".extract") unless($newfa); 
  return $newfa if($noupdate and -s $newfa);
  my ($ok,$hin)= openRead($fa); return 0 unless($ok);
  
  ## faids may be ref() HASH or filename
  if(ref($faids) =~ /HASH/) {
    %ids= %$faids;
  } else {
    $faids= makename($fa,".ids") unless($faids); 
    ($ok,$hids)= openRead($faids); return 0 unless($ok); 
    while(<$hids>) { if(/^\w/) { my($id)=split; $ids{$id}=1; } } close($hids);
  }
  rename($newfa,"$newfa.old") if(-s $newfa); #? OR option to return faids as is if exists...
  $ok= open($hout,'>',$newfa);   
  if($ok){ $ok=0; while(<$hin>) { if(/^>(\S+)/){ my $id=$1; $ok=$ids{$id}||0; } print $hout $_ if($ok); } }
  close($hout); close($hin); 
  return $newfa;  
}

sub fadupids #  in cdna_evigenesub.pm
{ 
  my $fa=shift; my($ndup)=(0,0); my %ids; 
  my($ok,$hin)= openRead($fa);
  if($ok) { while(<$hin>) { if(/^>(\S+)/) { my $id=$1; my $ni= ++$ids{$id}; 
    if($ni>1) { $ndup++; loggit(1,"ERR: dup id:$id in $fa"); } 
    } } close($hin); }
  # die "ERR: $ndup duplicate ids in $fa\n" if($ndup>0); # leave to caller ?
  return (wantarray) ? ($ndup,\%ids) : $ndup;
}

sub fasizes_nogap {  
  my($fa,$okayc,$gapc)= @_;
  $okayc ||= 'ACGTagct'; $gapc ||= 'Nn';
  if($okayc =~ /^aa|amino|prot/) { $okayc='A-WYZa-wyz'; $gapc='Xx\*'; }
  my ($nokay,$ngap,$nt,$id)=(0,0,0,0); # NOT total, per record
  my %fasizes=();
  my($ok,$hin)= openRead($fa);
  if($ok) { while(<$hin>) { 
    if(/^>(\S+)/) {
      if($id) { $fasizes{$id}= join"\t",$nokay,$nt,$ngap; }
      $id=$1; ($nokay,$ngap,$nt)=(0,0,0);
    } else {
      chomp; s/\*$//;
      $nt += length($_);
      $nokay += tr/$okayc/$okayc/; # tr/A-WYZa-wyz/A-WYZa-wyz/; ## aa chars; na chars=/ACGTagct/ gaps=/Nn/
      $ngap  += tr/$gapc/$gapc/;  # tr/Xx\*/Xx\*/;
    }
  } close($hin); }
  if($id) { $fasizes{$id}= join"\t",$nokay,$nt,$ngap; }
  return \%fasizes; # ($nokay,$nt,$ngap);
}


sub fasplit { #  this splits by size; need alternate split by count: facount/ncpu
  my($fa, $spldir, $npart, $splsize, $fasuf)=@_;
  $fasuf ||= "fa";
  my @splist= (); my ($ok,$hin)= (0,undef);
  my($atsize,$atfile)=(0,0);
  my $fabase= basename($fa);
  # if($fa =~ /\.gz$/) { $ok= open(F,"gunzip -c $fa|"); } else { $ok= open(F,$fa); }
  ($ok,$hin)= openRead($fa);
  unless($ok) { loggit(LOG_DIE,"ERR: fasplit $fa $spldir $npart $splsize"); return @splist; }
  mkdir($spldir) unless(-d $spldir);
  while(<$hin>) {
    if($atfile==0 or ($atsize > $splsize && /^>/)) {
      $atfile++; if($atfile>$npart) { } # finishup???
      close(SPL) if($atfile>1);
      my $spname= $spldir . makename($fabase,".split$atfile.$fasuf",$fasuf.'|tbl|fasta|fsa|fa'); ##  choose suffix
      $ok= open(SPL,'>',$spname);  $atsize= 0;
      unless($ok) { loggit(LOG_DIE,"ERR: fasplit $atfile $spname"); return @splist; }
      push @splist, $spname;
      }
    print SPL;
    $atsize+= (length($_) - 1) unless(/^>/);
  } 
  close(SPL); close($hin);
  return @splist;
}

sub fasplitcount { # alternate split by count: facount/ncpu
  my($fa, $spldir, $npart, $splcount, $fasuf)=@_;
  $fasuf ||= "fa";
  my @splist= ();
  my($ok,$hin, $irec,$atfile)=(0) x 10; 
  my $fabase= basename($fa); #? is this bad
  # if($fa =~ /\.gz$/) { $ok= open(F,"gunzip -c $fa|"); } else { $ok= open(F,$fa); }
  ($ok,$hin)= openRead($fa);
  unless($ok) { loggit(LOG_DIE,"ERR: fasplit $fa $spldir $npart $splcount"); return @splist; }
  mkdir($spldir) unless(-d $spldir);
  
  while(<$hin>) {
    if($atfile==0 or ($irec >= $splcount && /^>/)) {
      $atfile++;  $irec=0;
      close(SPL) if($atfile>1);
      my $spname= $spldir . makename($fabase,".split$atfile.$fasuf",$fasuf.'|tbl|fasta|fsa|fa'); ##  choose suffix
      $ok= open(SPL,'>',$spname); 
      unless($ok) { loggit(LOG_DIE,"ERR: fasplit $atfile $spname"); return @splist; }
      push @splist, $spname;
      }
    print SPL;
    $irec++ if(/^>/);
    ## $atsize+= length($_) unless(/^>/);
  } 
  close(SPL); close($hin);
  return @splist;
}


##upd1804 fasort
sub fasort {
  my($infa,$outf)=@_;
  our $LOGT=($dryrun||$DEBUG)? LOG_WARN : LOG_DIE; 

  sub fasort_put {
    my($sout,$id,$hd,$fa)=@_;
    my $sid=$id;
    my($ni,$ti)= ($id=~m/(\d+)t(\d+)/)?($1,$2):(0,0); #evigene id form
    if($ni) {
      (my $gp=$id)=~s/${ni}t${ti}.*//;
      $sid="$gp\t$ni\t$ti"; 
    } else { # alphasort id
      $sid="$id\t0\t0";
    }
    map{ chomp($_); s/\n/\t/g; } ($hd,$fa);
    print $sout "$sid\t$hd\t$fa\n";
  }

  sub fasort_sort {
   my($tmpsrt,$outfa)=@_; our($LOGT);
   my $ns=0;
   open(OUT,">",$outfa) or loggit( $LOGT, "write $outfa"); 
   open(SI,"sort -T ./ -k1,1 -k2,2n -k3,3n $tmpsrt |") or loggit( $LOGT, "fail sort $tmpsrt");  
   while(<SI>) {
     my ($sa,$sb,$sc,$val)= split"\t",$_,4; 
     $val=~s/\t/\n/g;  chomp($val);
     print OUT $val,"\n"; $ns++
   }
   close(SI); close(OUT);
   return($ns);
  }
  
  my $tmpsrt= makename($infa,"$$.so"); 
  $outf= makename($infa,".sorted") unless($outf); 
  my($nt,$id,$hd,$fa)=("") x 9;  
  my($ok,$hin)= openRead($infa);
  unless($ok) { loggit( LOG_WARN, "cant sort $infa"); return -1; } 
  open(SO,">",$tmpsrt) or loggit( $LOGT, "write $tmpsrt");  
  while(<$hin>) {
    if(/^>(\S+)/){ my $d=$1; 
      fasort_put(*SO,$id,$hd,$fa) if($fa); 
      $id=$d; $hd=$_; $fa=""; }
    else { $fa.=$_; }
  } close($hin);
  fasort_put(*SO,$id,$hd,$fa) if($hd);
  close(SO); 
  
  my($ns)= fasort_sort($tmpsrt,$outf);
  unlink($tmpsrt);
  return($ns);
}

## in cdna_evigenesub.pm
sub findapp
{
  my($aname, $nodie)=@_;
  my $app=""; $nodie ||= 0;
  $app=$ENV{$aname} if(not $app and $ENV{$aname});  
  $app=$ENV{uc($aname)} if(not $app and $ENV{uc($aname)});  
  $app=`which $aname` unless($app); 
  chomp($app);
  ## #tr2aacds: app=blastn, path=no blastn in 
  my $dol=0; if(not $app or $app =~ /^no $aname/) { 
    $app="echo MISSING_$aname"; 
    $dol=($dryrun||$DEBUG||$nodie)? LOG_WARN : LOG_DIE; 
    }
  loggit( $dol, "app=$aname, path=$app");
  return (wantarray) ? ($app,dirname($app)) : $app;
}
#edit for containers
sub findapp2
{
  my($aname, $nodie)=@_;
  my $app=""; $nodie ||= 0;
  $app=`which -a $aname | grep conda` unless($app);
  chomp($app);
  ## #tr2aacds: app=blastn, path=no blastn in 
  my $dol=0; if(not $app or $app =~ /^no $aname/) {
    $app="echo MISSING_$aname";
    $dol=($dryrun||$DEBUG||$nodie)? LOG_WARN : LOG_DIE;
    }
  loggit( $dol, "app=$aname, path=$app");
  return (wantarray) ? ($app,dirname($app)) : $app;
}
sub finddata # new 1712
{
  my($aname, $nodie)=@_;
  my $dfile="";  $aname||="finddata.missing";
  $nodie=1; # for now; $nodie ||= 0;
  $dfile=$ENV{$aname} if(not $dfile and $ENV{$aname});  
  $dfile=$ENV{uc($aname)} if(not $dfile and $ENV{uc($aname)});  
  #NO# $dfile=`which $aname` unless($dfile); 
  ## check file exists:  $isfile= -f $dfile
  ## BUT may be ncbi db prefix: UniVec == UniVec.{nsq,..} or directory
  ## also check file wildcards: *?
  my $dexists=0;
  if($dfile and -e $dfile){ $dexists=1; }  
  elsif(-f $aname){ $dfile=$aname; $dexists=1; }  
  elsif($aname =~ m/[\*\?]/) {
    my $dls=`/bin/ls $aname`;  # trap sys err msg 'no such file'
    chomp($dls); if($dls and -e $dls) { $dfile=$dls; $dexists=1; }
  } 
  my $dol=0; if(not $dfile or $dfile =~ /^no $aname/) { 
    #NO# $dfile="echo MISSING_DATA_$aname"; 
    $dol=($dryrun||$DEBUG||$nodie)? LOG_WARN : LOG_DIE; 
    }
  loggit( $dol, "data=$aname, exists=$dexists, path=$dfile");
  return (wantarray) ? ($dfile,dirname($dfile)) : $dfile;
}

## in cdna_evigenesub.pm
sub findevigeneapp
{
  my($aname, $nodie)=@_;
  my $app=""; 
  my $ename= basename($aname,'.pl');
  $app= $ENV{$ename} if($ENV{$ename});  
  $app= $aname if(not $app); $nodie ||= 0;
  # my $EVIGENES="$FindBin::Bin/.."; # ok?
  # my $APPcdnabest= findevigeneapp("$EVIGENES/cdna_bestorf.pl"); # allow ENV/path substitutions?
  $app="$EVIGENES/$aname" unless(-x $app);
  my $dol=0; 
  unless( -x $app) { 
    $app="echo MISSING_$aname"; 
    $dol=($dryrun||$DEBUG||$nodie)? LOG_WARN : LOG_DIE; 
    }
  loggit( $dol, "evigeneapp=$aname, path=$app");
  return($app);
}

sub cat_splitset
{
  my($tofile, $splitset)=@_;
  my $nin= (ref($splitset))? @$splitset : 0; 
  loggit( ($dryrun) ? 1 : 0,"CMD=","cat_splitset to $tofile from n=$nin,",$$splitset[0],".. ");  
  return 0 if($dryrun); return -1 if($nin<1);
  my (@qfail,@ofail);
  if(-s $tofile) { } # save?  
  my $nok= 0;
  my $ok= open(O,'>', $tofile);
  if($ok) { for my $i (0..$nin-1) {
    my $outf= $splitset->[$i];
    my $isend= 0;
    if($outf and -s $outf) {
      my $lastl="";      
      $ok= open(I,$outf); while(<I>) { $lastl=$_; print O $_; } close(I);
      $isend=$ok; #not blast# $isend++ if( $lastl =~ m/^# BLAST processed/); 
    }
    if($isend>0) { $nok++; } else { push @ofail, $outf; }
  }
  close(O); }
  my $runerr= $nin - $nok;
  if($runerr) { loggit(1,"ERR=$runerr ","cat_splitset to $tofile"); } 
  return ($runerr, $nok, \@ofail);
}

sub runcmd
{
  my @cmd= @_;
  ## fail if $cmd[0] =~ /MISSING_/
  my $isdry=($dryrun or $cmd[0] =~ m/^echo\b/)?1:0;
  loggit( ($isdry) ? 1 : 0,"CMD=",@cmd);  
  my $err= ($isdry) ? 0 : system(@cmd);  ## ?? add option run: @result=`$cmd`;
  if($err) { loggit(1,"ERR=$err ",$cmd[0]); } # ..
  return $err;
}

sub forkcmd
{
  my @cmd= @_;
  loggit( ($dryrun) ? 1 : 0,"forkCMD=",@cmd);  
  unless($dryrun) {
    my $pid= fork();
    if($pid) { # parent
      return $pid;
    } else { # child
      my $err= system(@cmd);
      exit($err);
    }
  }
  # if($err) { loggit(1,"ERR=$err ",$cmd[0]); } # ..
}

sub makename
{
  my($infile,$osuf,$insuf)=@_;
  #BAD?#  $insuf ||= 'aa|blast|cdna|mrna|cds|tr|trclass|tbl|fasta|faa|fsa|fa';  ## fixme need insuf: tr|fasta|fa
  ## in: dmagset56ri.tr.split4.fa : OUT dmagset56ri_split/dmagset56ri.aa  LOST split4; all splits to one.aa
  unless($infile) { warn "cdna_evigenesub:makename MISSING infile"; $infile="Noname$$"; } # or die / buggy soft
  my $outfile= $infile; $outfile =~ s/\.gz$//; # bad for infile empty/undef .. do what?

  #BADold# $outfile =~ s,\.($insuf)[^\/\s]*$,,; ## ADD \. cut only final suffix !!!
  ## or instead use: smallest at of $at=rindex($outfile,'.'.$insuf[i]); $outfile=substr($outfile,0,$at)
  #o2# $insuf ||= 'aa|blast|cdna|mrna|cds|fasta|fa';  ## fixme need insuf: tr|fasta|fa
  #o2# $outfile =~ s,\.($insuf)[^\.\/\s]*$,,; ## FIXMEd, or s,\.($insuf)\w*$,,; or s,\.($insuf)$,,;
  
  if($insuf) { $outfile =~ s,\.($insuf)[^\.\/\s]*$,,; } # insuf can have \.: 'aa.qual|aa.size|...'
  else { $outfile =~ s,\.\w*$,,;  } ## default insuf maybe chop any \.\w+$ instead of guess list?
  
  $outfile.= $osuf if($osuf); 
  $outfile.= "_out" if($outfile eq $infile);
  return $outfile;
}

=item evigene_config

	GetOptions(.. "config=s", \$configfile, "cadd=s", \@configadd,);
	$configh= evigene_config($configfile, \@configadd); # always even if $config null

	$DBID= $configh->{general}->{dbid} || "MyDBID"; 
	$LOCUSTAG= $configh->{general}->{locus_tag} || $DBID; 
	$IDPrefix= $configh->{pubopt}->{publicid} || "Evigene"; 
	
 	from  evigene/scripts/evigene2genbanktbl.pl
	but drop global special hashes: %evidence; %public_options; %geneset; %programs; 
	change to: return \%config # replace other config hashes w/ this 2-level: config{part}{key} = val

=cut

sub evigene_config {
  my($cfile, $addoptions)= @_;
  my %config=(); my $ctype=0;

  use constant{ kGENERAL => 'general', kEVFILE => 'evidence', kEVPROG => 'programs', kPUBOPT => 'pubopt',  };
  	#old: kEVOPT => 'evoption', kANOPT => 'anoption',  kEVGENES => 'geneset',
  	
  if($cfile) { #  and -f $cfile
    open(F,$cfile) or die "ERROR reading config: $cfile";
  
    my @CONFIG= <F>;
    push @CONFIG, @$addoptions if(ref $addoptions);
    
    my ($lastkey, $lastval);
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

# revise to parse '^section:' into separate hash
      if($key =~ s/:$//) { $ctype=$key; }
      
      if($key =~ /^evidence/) { $ctype= kEVFILE; } # old/fixed set
      elsif($key =~ /^pubopt/) { $ctype= kPUBOPT; }
      elsif($key =~ /^program/) { $ctype= kEVPROG; }
      elsif($key =~ /^end$/) { $ctype= 0; }
#       elsif($key =~ /^evoption/) { $ctype= kEVOPT; }
#       elsif($key =~ /^anoption/) { $ctype= kANOPT; }
#       elsif($key =~ /^geneset/) { $ctype= kEVGENES; }
#.. drop special config hashes...
#       elsif($ctype eq kEVFILE ) { $evidence{$key}= $val; } # /gff/ was bad for geneset
#       elsif($ctype == 0 and $val =~ /\.gff/) { $evidence{$key}= $val; } 
# #       elsif($ctype eq kEVOPT) { $evaluate_options{$key}= $val; } 
# #       elsif($ctype eq kANOPT ) { $annotate_options{$key}= $val; } 
# #      elsif($ctype eq kEVGENES) { $geneset{$key}= $val; }
#       elsif($ctype eq kPUBOPT ) { $public_options{$key}= $val; } 
#       elsif($ctype eq kEVPROG) { $programs{$key}= $val; }

      # generic keys: name date genome .. other?
      if($key =~ /\w/ and $val =~ /\w/) { 
        my $ogroup= $ctype || "general";
        $config{$ogroup}->{$key}= $val; # this one #which?  $config{$ogroup}{$key}= $val;
      }
      
      # also for now : overlap, other progs
      ($lastkey, $lastval)=($key, $val);
    } close(F);
  }
  return \%config;
}


1;
__END__
