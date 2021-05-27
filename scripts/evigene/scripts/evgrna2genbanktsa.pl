#!/usr/bin/env perl
# evgrna2genbanktsa.pl or mrna2tsa.pl

=item notes

  EvidentialGene mrna2tsa.pl
  process tr2aacds.pl outputs for ncbi tsa  submit
  
  -- main/alternate id table
  -- pubids version of main/alt ids
  -- vecscreen mrna.tr
  -- asmrna2ncbitsa.pl process vecscreen data, annot table?

  parts from
# evigene2genbanktbl.pl
# asmrna2ncbitsa.pl
# bestgenes_update.pl 
# bestgenes_puban_kfish.pl ??

  $evigene/scripts/rnaseq/asmrna2ncbitsa.pl -GAPSOK -idpre Thecc1ER_ \
  -cdna ../tr5parts/pub3ig.$pt.tab4g.tr.gz -vec ../tr5parts/pub3ig.trasm.tab4g.vector.tab \
  -geneinfo ../tr5parts/pub3ig.trasm.tab4g.geneinfo1.tab  -log tr4g.$pt.log \
  -out $pt/TCM01.tsa_rasm.$pt.fsa -tbl $pt/TCM01.tsa_rasm.$pt.tbl

=cut

use FindBin;
use lib ("$FindBin::Bin"); # assume evigene/scripts/ layout: main, prot/, rnaseq/, ...

use strict;
use Getopt::Long;
#maybe# use cdna_proteins;

use constant VERSION => '2013.04.16'; # 03.20

warn "OBSOLETE: evgrna2genbanktsa.pl .. USE: evigene/scripts/evgmrna2tsa.pl 
http://arthropods.eugenes.org/EvidentialGene/about/EvidentialGene_trassembly_pipe.html\n\n";

## evigene path = self path
my $EVIGENES="$FindBin::Bin"; #??

my $debug= 0;
my $GENEINFO_VERS=2; # constant
my $MAXGAP=15; # NCBI 
my $ENDGAP=10; # trim ends if gaps w/in this of ends; NCBI, was 10
my $MINSIZE=200; # NCBI
my $GAPSOK=1; # default on? new policy 2012Dec for TSA tbl2asn: -a r10u -l paired-ends
my $MINGENEIDENT=85; # for asm == gene identity
my $IDPREFIX= $ENV{idprefix} || "evgr"; #  opt
my $GDB_PREFIX='gnl|Evigene|';  #see below; 'gnl|Evigene|'; # use IDPREFIX ? or not, since ID has this also
my $DATE=0;
my $pubidnum_start=0;

## namegenes messed up, didn't screen out loqualnames:
##  1%,309/22971,716        Dumpy
##  0%,126/34350,248        Titin
our $MIN_NAMEIDENT = 35;  # min similar% for protein evid align/equivalence, for naming; JCVI uses 35%  
our $MIN_IDLIKE   = 15;  # low for now; 15..20 seems right;add Note with loqualname

my ($vecscreenf,$trclass,$geneinfo,$genenames,$cdnaseq,$output,$logfile,$tblfile,$dryrun);

my $optok= GetOptions(
  # "config=s", \$config,
  "mrna|cdna=s", \$cdnaseq,
  "class|trclass=s", \$trclass,
  "vectors|vecscreen=s", \$vecscreenf,  
  "geneinfo=s", \$geneinfo,  #?? separate naming table?
  "names|genenames=s", \$genenames,   
  
  "output:s",  \$output,
  "tblfile:s", \$tblfile,  
  "logfile:s", \$logfile,
  # "annot|score=s", \$inscore,
  # "genome=s", \$genome,
  "idprefix=s", \$IDPREFIX,  
  "DATE=s", \$DATE,  
  "gdbprefix=s", \$GDB_PREFIX, #?? IDPREFIX default 
  # "proteins=s", \$proteins, # maybe
  # "version|MSRC=s", \$MySRC,
  # "cadd=s", \@configadd,
  "MINSIZE=i", \$MINSIZE,  
  "MAXGAP=i", \$MAXGAP,  
  # "MINGENEIDENT=i", \$MINGENEIDENT,  ## not used?
  "GAPSOK!", \$GAPSOK, 
  "dryrun|n!", \$dryrun, 
  "debug!", \$debug, 
  );


# OPTION: -class evg_tr2aacds.trclass only requirement, gives paths to mrna, names.tab ?
# OPTION: here or caller? make mrna/cdnaseq  from evg traa2cds.pl okayset/name.{okay,okalt}.{tr,aa}  
# $evigene/scripts/prot/traa2cds.pl -trout -cdna $pt.tr.gz -aa $pt.aa.gz -out -log
##  -trout : cdna output, not cds, as name.mrna.tr

#?No# $cdnaseq= shift @ARGV unless($cdnaseq);

die "usage: evgrna2genbanktsa.pl -mrna inmrna.fasta -class inmrna.trclass ...
 opts: -out=outfasta -idprefix Thecc1RA_ -vecscreen=infile -geneinfo=infile  -conf=evigene.conf
    -MINSIZE=$MINSIZE  -MAXGAP=$MAXGAP  -log=outlog\n"
  unless($optok and ($cdnaseq or $trclass));  

use constant { LOG_NOTE => 0, LOG_WARN => 1, LOG_DIE => -1, };
my $logh= undef;
sub loggit{ my $dowarn=shift; 
  my $s= join(' ',@_); chomp($s); $s="FATAL $s" if($dowarn == LOG_DIE);
  if($logh){ print $logh "#er2g: $s\n"; } elsif($dowarn>0||$debug){ warn "#er2g: $s\n"; }
  if($dowarn == LOG_DIE) { die "#er2g: $s\n" ; }
}

if(not $logfile and defined $logfile) { # use output name
  $logfile= $cdnaseq || $trclass;
  $logfile= makename($logfile,".mrna2tsa.log");  
}
if($logfile) { open(LOG, ">>$logfile") or die $logfile; $logh= *LOG; }

my $GAPSMAX = ('N') x $MAXGAP;
# evigene_config($config, \@configadd); # always even if $config null

my $nd= ( $IDPREFIX =~ s/(0\d+)$// ) ? length($1) : 6;
my $pubid_format = $IDPREFIX.'%0'.$nd.'d'; # $public_options{'publicid'} || "evgr000000";
my $altid_format = 't%d'; # t%d or t%02d # $public_options{'altid'} || "t00";
#? No?# unless($GDB_PREFIX) { $GDB_PREFIX= "gnl|$IDPREFIX|"; } #?

unless($DATE) { $DATE=`date '+%Y%m%d'`; chomp($DATE); } # fixme
my $GBPROID= $IDPREFIX."_".$DATE; # "cacao11evigene_20120827";

loggit(1, "EvidentialGene mrna2tsa.pl VERSION",VERSION);
loggit(1, "ERR: unused arguments:",@ARGV) if(@ARGV>0);

my $APPvecscreen=    findapp("vecscreen"); # ncbi/bin/... also need ncbi/data/UniVec
# my $APPmakeblastdb= findapp("makeblastdb");
# my $APPcdnabest= "$EVIGENES/cdna_bestorf.pl"; # allow ENV/path substitutions?
# my $APPtraa2cds= "$EVIGENES/prot/traa2cds.pl";
# my $APPtrdupfilter= "$EVIGENES/rnaseq/asmrna_dupfilter2.pl"; # complex call, inputs
# my $APPaaqual=    "$EVIGENES/prot/aaqual.sh";
my $APPtraa2cds= findevigeneapp("prot/traa2cds.pl");

## FAIL at this point if any apps missing?
#-------------------------------------


my( %vecscreen, %gene, %genenames, %genedbxref, %genenamepct, %namedgenes, %pubids, %geneinfo, $trpath,$trname);

($cdnaseq,$trpath,$trname)= get_evgtrset($trclass,$cdnaseq,);
  loggit(0, "get_evgtrset=",$cdnaseq,$trpath,$trname); ## facount($cdsseqnr)
  loggit(LOG_DIE, "Missing -cdnaseq",$cdnaseq) unless($cdnaseq and -s $cdnaseq);

unless($genenames) { 
  my $gnt="$trpath/$trname.names";  $genenames=$gnt if(-s $gnt); 
  loggit(LOG_WARN, "Missing product -names",$gnt) unless($genenames);
}


my($maintab,$pubids,$nmaintr,$nalltr)= trclass2maintab($trclass);
loggit(0, "trclass2maintab primary n=",$nmaintr,"allntr=",$nalltr,$pubids); 
if($pubids) { # require?
  open(F,$pubids) or loggit(1,"ERR: reading $pubids");
  while(<F>) { chomp; my($id,$oid,$gid,$alti)=split"\t"; 
     $pubids{$oid}= $id; $pubids{$id}="$oid\t$gid\t$alti";
  } close(F);
}

($vecscreenf)= vecscreen($cdnaseq,$vecscreenf);
if($vecscreenf) { 
  my %vids=();
  open(F,$vecscreenf) or loggit(1,"ERR: reading $vecscreenf");
  while(<F>) { chomp; my($vd,$vb,$ve,$vt)=split"\t"; $vecscreen{$vd}.="$vb\t$ve\t$vt\n"; $vids{$vd}++; } close(F);
  my $nvid=scalar(keys %vids); 
  loggit(0, "vectors found in ntr=",$nvid,$vecscreenf); 
}

my($outfa, $tblout, $annot,$notr,$nocds) = trprocess($cdnaseq,$trclass); # add main/alt pub ids, other geneinfo 
loggit(0,"DONE output ntr=$notr, ncds=$nocds in files $maintab, $pubids, $outfa, $tblout, $annot"); 

#---------------------------------------------------  

## make this subrtn: 
# my($cdnaseq,$trpath,$trname)= get_evgtrset($trclass,$cdnaseq,);

sub get_evgtrset {
  my($trclass,$cdnaseq)= @_;
  my($trpath,$trname)=("","");
  
  ## add gene name file here?
  # unless($genenames) { my $gnt="$trpath/$trname.names"; $genenames=$gnt if(-f $gnt); } # look for it

  if($trclass) {
    my $trpname= makename($trclass,"","trclass"); 
    if($trpname =~ m,/,) { ($trpath,$trname)= $trpname =~ m,(.*)/([^/]+)$,; } # BADDDDD
    else { $trname=$trpname; }
    $trpath ||= '.';  
    return($cdnaseq,$trpath,$trname) if($cdnaseq);
    
    my $okpath="$trpath/okayset";
    if(-d $okpath) { # got tr2aacds outset
      opendir(D,$okpath); 
      my @okd= map{ "$okpath/$_" } readdir(D); 
      closedir(D);
 
 ## Ugh! bad call: #er2g: get_evgtrset= ./okayset/pogonus1all3.mrna0.ann.txt . pogonus1all3
 ## need grep /\.mrna$|\.mrna.gz$|.mrna.fa$/ ??? 
      my ($trf)= grep /\.mrna$|\.mrna.gz$/, @okd;
      if($trf and -s $trf) { $cdnaseq= $trf; }
      else { 
        my $cdnatmp="$okpath/$trname.mrna";
        loggit(0,"Make cdna $cdnatmp from okayset transcripts");
        my($oktr) = grep /.okay\.tr/, @okd;  (my $okaa=$oktr) =~ s/\.tr/.aa/;
        my($alttr)= grep /.okalt\.tr/, @okd; (my $altaa=$alttr) =~ s/\.tr/.aa/;
        if($oktr and -f $oktr and -f $okaa) { 
          runcmd("$APPtraa2cds -trout -cdna $oktr -aa $okaa -out stdout >> $cdnatmp");
        }
        if($alttr and -f $alttr and -f $altaa) { 
          runcmd("$APPtraa2cds -trout -cdna $alttr -aa $altaa  -out stdout >> $cdnatmp");
        }
        $cdnaseq= $cdnatmp if(-s $cdnatmp);
        }
      }
  }
  
  return($cdnaseq,$trpath,$trname);
}


sub parse_evgheader
{
  my($oid,$hdr,$trlen)= @_;

  my $pubid= $pubids{$oid} || $oid; # is it ERR if no pubid{oid} ?

  my $protid= $pubid;
  $protid=~s/t(\d+)$/p$1/; # genbank requires diff protid from mrnaid
  $protid= $GDB_PREFIX.$protid if($GDB_PREFIX);
  #  protein_id      gnl|CacaoGD|Thecc1EG016762p1

  my %tblinfo= (pubid => $pubid, oid => $oid,  protid => $protid, locustag => 0,
      aaqual => "na", trlen => $trlen, cdsoff => "", cdsor => 1, 
      name => "", namepct => 0, dbxref => "na" ); 

  if( $genenames and $genenames{$oid} ) {
    $tblinfo{'name'}= $genenames{$oid};
    $tblinfo{'namepct'}=  $genenamepct{$oid} || 0;
    $tblinfo{'dbxref'}=  $genedbxref{$oid}||"na";
  }
        
  my($cdsb,$cdse,$aafull)=(0,0,0);
  if($hdr =~ m/\boffs=([\d-]+)/) { my $cdsoff=$1; $tblinfo{'cdsoff'}= $cdsoff; 
    ($cdsb,$cdse)= split/[-]/,$cdsoff;  } # do in putseq
  if($hdr =~ m/\b(?:aaqual|aalen)=([^\s;]+)/) { my $aq=$1; $tblinfo{'aaqual'}= $aq; 
    ($aafull)= $aq =~ m/(complete|partial\w*)/; }
  if($hdr =~ m/\bclen=(\d+)/ and not $trlen) { $trlen=$1; $tblinfo{'trlen'}= $trlen; } # skip? 
  if($hdr =~ m/\bstrand=([^\s;]+)/) { $tblinfo{'cdsor'}= $1; } # expect all '+' from traa2mrna
  if($hdr =~ m/\b(?:[Nn]ame|[Pp]roduct)=([^=\n;]+)/) { my $na=$1; 
     $tblinfo{'name'}= $na unless($tblinfo{'name'}); }

  return \%tblinfo;
}

  # my($outfa,$tblout,$annot,$notr,$nocds) = trprocess($cdnaseq,$trclass);    
sub trprocess
{
  my($cdnaseq,$trclass)=@_;
  my($notr,$nocds)=(0,0);
  my $outfa = ($output)  ? $output  : makename($cdnaseq,".fsa"); # was .fna ; tbl2asn wants fsa
  my $tblout= ($tblfile) ? $tblfile : makename($cdnaseq,".tbl");
  my $annot =  makename($cdnaseq,".ann.txt"); ## FIXME: change suffix: was .annotab
  return($outfa,$tblout,$annot,$notr,$nocds) if( -s $outfa); # or dryrun ..
  ## maybe option: $REGEN_AA
  (my $outaa=$outfa) = s/\.fsa/.faa/;
  
  if($genenames) {
    # FIXME2: ** use uniq names vs ERR: too short to keep valid tr, e.g. 
#er2g: ERR: too short:183 #LitvaEG0018688t4     oid=litovavel2k35Loc15824t1     len=183 name=CDD: clpS, ATP-dependent Clp protease ada..-like    
# grep  'CDD: clpS,' *.names = 1 only = litovavel2k35Loc15824t1
    
    # FIXME: need better reader; 2+ rows/id; pick best .. format may change..
    # names.tab ==  id, name, pctalign, refid, repid  : now
    #  trid1  C-ets-2 protein, putative       89%,103/116,197 RefID:UniRef50_E0VFI2   RepID:E0VFI2_PEDHC
    #  trid2  DBH-like monooxygenase protein 1        73%,445/613,516 RefID:UniRef50_Q6UVY6   RepID:MOXD1_HUMAN
    open(F,$genenames) or loggit(1,"ERR: reading $genenames");
    while(<F>) { 
      chomp; 
      my($id,$name,$pctalign,$refid,$repid)=split"\t"; # may have only id, name
      my $xtra; ($name,$xtra)=split";",$name,2; 
      
      # FIXME: 2 names/id maybe: CDD: xxx and gene xxx; keep both in ann.txt ? and pctalign?
      ## old geneinfo:
      #  ($td,$gd,$mapqual,$aaqual,$trlen,$cdsor,$cdsoff,$lotag,$nam1,$dbxref)= @$rinfo;
      ## usage below in putseq
      # my %tblinfo= (pubid => $pubid, oid => $oid,  protid => $protid, locustag => 0,
      #  aaqual => "na", trlen => 0, cdsoff => "", cdsor => 1, 
      #  name => $genenames{$oid}||"", dbxref => $genedbxref{$oid}|| "na"); 

      if($pctalign =~/^\d/ and $pctalign < $MIN_NAMEIDENT) { # use MIN_IDLIKE : name-like ? got some '0%' align
        ## bad: Uncharacterized protein-like ; Nuclease HARBI1, putative-like; Protein-like
        if($pctalign >= $MIN_IDLIKE) { unless($name =~ /\blike|^Uncharacterized/) {
          $name =~ s/, putative//; $name .= '-like'; }
        } else { next; } ## should we preserve for ann.txt table ?
      }
      
      ## fixme: CDD:206692,cd04107,RefID:UniRef50_Q9NX57,UniProt:RAB20_HUMAN, 
      ## drop  RefID:; drop? cd04107
      $refid =~ s/RefID://; 
      $genedbxref{$id} .= "$refid," if($refid);
      $genedbxref{$id} .= "$repid," if($repid and not $refid =~ /^CDD:/);
      $namedgenes{$name} .= "$id,"; #? if($pctalign >= $MIN_NAMEIDENT); # for uniq name retention
      
      unless($genenames{$id} and $name =~ /CDD:/) { # or refid =~ /CDD:/
        $pctalign ||= 0; $refid ||= 0;
        $genenames{$id}= $name;
        $genenamepct{$id}= $pctalign;
        ## $genedbxref{$id}= ($repid) ? $repid : $refid;  # do list here for all dbxref
         # repid Yes or No? this is by default RefID=UniRef50_xxxx and RepID=UniProt:xxxx_HUMAN
      }
    } close(F);
  }

  my ($inh,$outh,$tblh,$annoth,$hd,$oid,$fa,$ok);
  ## cdnaseq eq stdin also?
  if($cdnaseq =~ /\.gz$/) { $ok= open($inh,"gunzip -c $cdnaseq|"); }
  else { $ok= open($inh,$cdnaseq); }  
  $ok= open($outh,'>',$outfa) if($ok);
  $ok= open($tblh,'>',$tblout) if($ok);
  $ok= open($annoth,'>',$annot) if($ok);
  # $ok= open($outaah,'>',$outaa) if($outaa and $ok);
  unless($ok) { loggit(1,"ERR: trprocess $cdnaseq TO $outfa,$tblout"); return; }
  
  # replace >oid >pubid; but need to keep oid for other input tables.
  # add gene names?
  my($itr,$otr,$ocds,$oerr)= (0) x 10;
  while(<$inh>) { 
    if(/^>(\S+)/) { my $d=$1; 
      ## add: $outaah
      ($otr,$ocds,$oerr)= putseq($outh,$tblh,$annoth,$oid,$hd,$fa,$itr) if($fa); 
      $notr+= $otr; $nocds+= $ocds;
      $oid=$d; $hd=$_; chomp($hd); $fa=""; $itr++;
      }
    elsif(/^\w/) { chomp; $fa.=$_; }
  } 
  
  ## add: $outaah
  ($otr,$ocds,$oerr)= putseq($outh,$tblh,$annoth,$oid,$hd,$fa,$itr) if($fa); # last
  $notr+= $otr; $nocds+= $ocds;
  close($inh); close($outh);
  return($outfa,$tblout,$annot,$notr,$nocds); 
}


=item tr2main

 aabugs4/aabugs4qual/tsaevg/daphmag5_evgt2c_2013mar7.txt
 ## .. this main{md} isnt right either, as alt-2ndid can be other alt.

  grep okay $pt.trclass | perl -ne'($td,$ok,$cl,$md,$piad,$aq,$fl)=split; $cl=~s/a2$//; \
  ($pi,$pa,$pd)=split"/",$piad; $md=$pd if($pd); if($cl =~ /main|noclass/) { $main{$td}=$cl; } \
  else { $alt{$md}{$td}= $cl; $balt{$td}=$md; } $n++; \
  END { @amain= grep { not $main{$_} } sort keys %alt; \
  for $am (@amain) { $md= $balt{$am}; \
  if($main{$md}) { @at=keys %{$alt{$am}}; map{ $alt{$md}{$_}=$alt{$am}{$_} } @at; } \
  elsif($md) {  $main{$am}="NOMAIN";  } } \
  foreach $md (sort keys %main) { @ad= sort{$alt{$md}{$a} cmp $alt{$md}{$b}} keys %{$alt{$md}};\
  $ad=join",",map{ "$_/".$alt{$md}{$_} } @ad; $mc=$main{$md}; print join("\t",$md,$mc,$ad)."\n"; } }'\
    > okayset/$pt.mainalt.tab

=cut

sub trclass2maintab
{
  my($trclass)=@_;
  my $maintab = makename($trclass,".mainalt.tab","trclass");  # > $pt.mainalt.tab
  my $pubidtab= makename($trclass,".pubids","trclass");   
  return($maintab,$pubidtab) if( -s $maintab and -s $pubidtab); # or dryrun ..
  
  my(%main,%alt,%balt,%drop,$outh,$outpubidh,$inh);
  my $ntr=0;
  ## my $ok= open($inh,"grep okay $trclass |");
  my $ok= open($inh,$trclass);
  $ok= open($outh,'>',$maintab) if($ok);
  $ok= open($outpubidh,'>',$pubidtab) if($ok);
  unless($ok) { loggit(1,"ERR: parse $trclass TO $maintab"); return; }

  ## FIXME: only althi are reliably locus alternates; altmid .. are more likely paralogs
  while(<$inh>) {
    my($td,$ok,$cl,$md,$piad,$aq,$fl)=split;
    unless($ok eq 'okay') { $drop{$td}=1; next; } # OPTION: include drops?
    $cl=~s/a2$//;  #$n++; 
    my($pi,$pa,$pd)=split"/",$piad; $md=$pd if($pd); 
    if($cl =~ /^main|^noclass/) { $main{$td}=$cl;  $balt{$td}=$td; } 
    else { $alt{$md}{$td}= $cl; $balt{$td}=$md; }  
  }

  ## Fix MISSING main links, from alt to other alts ..  
  ## FIXME2: some of these NOMAIN are drops *** dont retain;
  ## FIXME3: adding drop{xxx} has screwed up alt-links somewhere; ** STILL MESSED UP
  # ... use drop{xx} only on output?
  
  my %hasmain;
  my @amain= grep { not $main{$_} } sort keys %alt; 
  foreach my $am (@amain) { 
    # next if($drop{$am});
    my $md= $balt{$am} || $am; ## $md=$am if($drop{$md});
    # next if($drop{$am} or $drop{$md}); #??
    if($main{$md}) { my @at=keys %{$alt{$am}}; map{ $alt{$md}{$_}=$alt{$am}{$_}; $balt{$_}=$md; } @at; } 
    elsif($md) {  $main{$am}="NOMAIN";  } # "nomain" ?
  }
  foreach my $td (keys %balt) {
    # next if($drop{$td});
    my $md= $balt{$td} || $td; 
    # $md=$td if($drop{$md}); # what ??
    $main{$md}="NOMAIN" unless($main{$md});
  }
     
  my $mainindex= $pubidnum_start;
  foreach my $md (sort keys %main) { 
    my @ad= sort{$alt{$md}{$a} cmp $alt{$md}{$b}} keys %{$alt{$md}}; # not here?  grep { ! $drop{$_} } 
    my $ad=join",",map{ "$_/".$alt{$md}{$_} } @ad; 
    my $mc=$main{$md}; 
    
    print $outh join("\t",$md,$mc,$ad)."\n"; 
    if($outpubidh) { # should be required ??
      $mainindex++;
      my $ialt= 0;
      unless($drop{$md}) {
      my ($pubmrnaid,$pubgeneid)= make_pubid($md, $mainindex, ++$ialt);
      print $outpubidh join("\t",$pubmrnaid,$md,$pubgeneid,$ialt)."\n"; $ntr++;
      }
      foreach my $ad (@ad) {
        unless($drop{$ad}) {
        my ($altmrnaid,$altgeneid)= make_pubid($ad, $mainindex, ++$ialt);
        print $outpubidh join("\t",$altmrnaid,$ad,$altgeneid,$ialt)."\n"; $ntr++;
        }
      }
    }
  }
  close($inh); close($outh);
  return($maintab,$pubidtab,$mainindex,$ntr);  # return main,alt id hashes ....
}

sub make_pubid
{
  my($oid, $mainindex, $altnum)= @_;

  ## use/check oid? keep global hash for pubid <=> oid ?
  ## my $pubid_format = $IDPREFIX.'%06d'; # $public_options{'publicid'} || "evgr000000";
  ## my $altid_format = 't%d'; # t%d or t%02d # $public_options{'altid'} || "t00";

  # $pubidnum_start++ if($altnum == 1); ## or $mainindex == last index?
  # my $pubidnum= $pubidnum_start;  # ONLY if altnum == 1 ? or if not seen this case..
  my $pubidnum= $mainindex;
  $pubidnum_start= $pubidnum; #?
  
  my $pubgene = sprintf( $pubid_format, $pubidnum); 
  my $pubid   = $pubgene . sprintf( $altid_format, $altnum);
  return($pubid,$pubgene);
}


=item putseq

  ## putseq() from asmrna2ncbitsa.pl
  Note: evg mrna.tr may have basic cds annot info
  okayset/dmag5xau13c2011_okmrna.tr.gz
  >sodmag4nalk25loc591t24 type=cdna; aalen=906,76%,complete; clen=3562;  strand=+; offs=93-2813;
  >sodmag4nalk25loc2322t3 type=cdna; aalen=812,86%,complete; clen=2827;  strand=+; offs=203-2641;

=cut

  ## see evigene2genbanktbl.pl:putCDSloc()
sub putCDSloc 
{
  my($cdsoff,$partial,$cdsphase)= @_;
  
  ## is cdsoff == "<123->456" allowed here?
  my($start,$stop)= split/[-]/,$cdsoff; # or $cdsoff =~ m/(\d+)-(\d+)/; # 
  
  my($p5,$p3,$codonstart)= ("<",">",0);
  if($partial =~ /complete/) { }
  else {
    unless($partial =~ /partial[53]/) { $partial.="53"; } # both
    if($partial =~ /3/) { $stop="$p3$stop"; }
    if($partial =~ /5/) { $start="$p5$start";  
      $codonstart=$cdsphase+1; # is this right?
    }
  }
  my $tbl= join("\t",$start,$stop,"CDS\n");
  $tbl .= "\t\t\tcodon_start\t$codonstart\n" if($codonstart>0);
  return $tbl;
}

sub putseq 
{
  ## add param: $outaah
  my ($outh, $tblh, $annoth, $oid, $hdr, $fa, $itr)=@_; 
  my $tlog="";  
  
  # tblinfo: my($td,$gd,$mapqual,$aaqual,$trlen,$cdsor,$cdsoff,$lotag,$nam,$dbxref)= @$rinfo;
  # maybe write all tblinfo as separate table, more info than .tbl .. put in logfile instead of below crap?

  my $ol= length($fa); 
  my($pubid,$def,$tblinfo);
  my($cdsb,$cdse,$cdsphase,$aafull,$ntrout,$ncdsout)=(0) x 10;

# if(1) {
  $tblinfo= parse_evgheader($oid,$hdr,length($fa));
  
  $pubid= $tblinfo->{'pubid'};
  $def= $pubid; # only this for ncbisubmit.fsa ?
  
  my $cdsoff= $tblinfo->{'cdsoff'}; ($cdsb,$cdse)= split/[-]/,$cdsoff; 
  my $aq= $tblinfo->{'aaqual'};    ($aafull)= $aq =~ m/(complete|partial\w*)/; 
  
#  #..........................
# } else {
#   
#   $pubid= $pubids{$oid} || $oid; # HERE or caller? is it ERR if no pubid{oid} ?
#   $def= $pubid; # only this for ncbisubmit.fsa ?
#   #NO? if(my $gname= $genenames{$oid}) { $def.=" [product=$gname]"; } # not prefered, goes in tblh
#   
#   my $protid= $pubid; #  $gd ? gbasn requires protein_id here, not same as mrna id
#   $protid=~s/t(\d+)$/p$1/; 
#   $protid= $GDB_PREFIX.$protid if($GDB_PREFIX);
#   #  protein_id      gnl|CacaoGD|Thecc1EG016762p1
# 
#   $tblinfo= { pubid => $pubid, oid => $oid,  protid => $protid, locustag => 0,
#       aaqual => "na", trlen => 0, cdsoff => "", cdsor => 1, 
#       name => $genenames{$oid}||"", dbxref => $genedbxref{$oid}|| "na" }; 
#       
#   if($hdr =~ m/\boffs=([\d-]+)/) { my $cdsoff=$1; $tblinfo->{'cdsoff'}= $cdsoff; 
#     ($cdsb,$cdse)= split/[-]/,$cdsoff;  }
#   if($hdr =~ m/\b(?:aaqual|aalen)=([^\s;]+)/) { my $aq=$1; $tblinfo->{'aaqual'}= $aq; 
#     ($aafull)= $aq =~ m/(complete|partial\w*)/; }
#   if($hdr =~ m/\bclen=(\d+)/) { $tblinfo->{'trlen'}= $1; }
#   if($hdr =~ m/\bstrand=([^\s;]+)/) { $tblinfo->{'cdsor'}= $1; }
#   if($hdr =~ m/\b(?:[Nn]ame|[Pp]roduct)=([^=\n;]+)/) { my $na=$1; 
#      $tblinfo->{'name'}= $na unless($tblinfo->{'name'}); }
# }  
  
## FIXME: vectrim in CDS should be disallowed *@#($ causing loss of stop codons in valid/hi-homology cds
  my($vectrimw)=(0,0);
  if(my $vec=$vecscreen{$oid}) { 
    # FIXMEd: multi locs per vd; sep by \n
    my @vec=split"\n",$vec; my $nv= @vec;
    foreach $vec (@vec) {
      my($vb,$ve,$vt)=split"\t",$vec; 
      if($ve > $cdsb and $vb <= $cdse and $aafull !~ /partial3|partial$/) { $vb=1+$cdse; } # dont trim stop codon !!!
      next if($vb >= $ol);
      my $trimw=1+$ve-$vb; $vectrimw += $trimw;
      substr($fa,$vb-1,$trimw)= ("N") x $trimw; 
      }
    } 


  # geneinfo for ncbi tsa submit: move this to feature.tbl not cdna.fsa; add cds offsets
  # >Feature Thecc1RA_L_g13025t00001
  #  bg  eg  gene  locus_tag TCM_000002
  #  bc  ec  CDS   product Cystathionine beta-synthase (CBS) family protein isoform 1
  #                protein_id  gnl|CacaoGD|Thecc1EG000002p1
  
  
  my $nn= $fa =~ tr/Nn/Nn/; 
  my $ncut=0;
  if($nn>0) { 
    # ** FIXME: cdsb,cdse adjust for inner gaps
    # fixme2: must adjust cdsb,e when cut BEFORE cdsb
    # fixme3: this is a mess; better to a. cut NNN, b. rerun cdna_bestorf for new cds offset?
    my ($lcdsb,$lcdse)= ($cdsb,$cdse);
    $fa=~s/n/N/g;
    ## SEQ_INST.HighNContentStretch: stretch of at least 5 Ns within the last 10 bases
    my $ne= rindex($fa,'N'); if($ne >= $ol - $ENDGAP) { 
      $fa=substr($fa,0,$ne); if($fa=~s/(N+)$//) {  my $ncut=length($1); $ne-=$ncut; }
      if($ne < $cdse) { $cdse = $ne; } #??
      }
    ## FIXME: cds-phase/codon_start changes w/ mod 3 of n1   
    my $n1= index($fa,'N'); if($n1 <= $ENDGAP) { 
      $n1++; $fa= substr($fa,$n1);  if($fa=~s/^(N+)//) { my $ncut=length($1); $n1+=$ncut; }
      if($cdsb>0) { $cdsb -= $n1; $cdse -= $n1;  } # $cdsphase = $n1 % 3; ??
      }
      
    $ncut=0; my $gapw= length( $GAPSMAX); #== MAXGAP
    unless($GAPSOK) {
    for (my $in= index($fa,$GAPSMAX); $in >= 0; ) {
      my $w=length($fa); my $en=$in+$gapw; 
      $en++ while($en<$w and substr($fa,$en,1) eq "N"); 
      my $wn= $en-$in; my $keep= 3 + ($wn % 3); my $cut= $wn-$keep; $ncut+=$cut; 
      my $facut= substr($fa,0,$in).substr("NNNNNN",0,$keep).substr($fa,$en); 
      $fa=$facut; 
      if($cdse>0) {
        if($en < $cdsb) { $cdsb -= $cut; $cdse -= $cut; } ##  $cdsphase = $cut % 3;
        elsif($in < $cdse and $en > $cdsb) {
          if($in <= $cdsb) { $cdsb -= $cut; $cdse -= $cut; } #??  $cdsphase = $n1 % 3; ??
          else { $cdse -= $cut; }
        }
      }
      $in=index($fa,$GAPSMAX); 
    } 
    }

    ## FIXME: cant have neg cdsb ..
    #er2g: >dmag5xevgr001990t1      oid=dmag4vel4xbxk75Loc4074t1    len=607; olen=615; cut=0; nnn=12/13;    cutcds=-5-501,oldcds=3-509,
    unless($cdse == $lcdse and $cdsb == $lcdsb) {
      if($cdsb <= 0) { 
        $cdsphase = (2 + $cdsb) % 3; # ?? this seems right, dont know why !!
        $cdsb=1;  my $aafull0= $aafull;
        $aafull="partial5" if($aafull=~/complete/); # 
        $aafull="partial"  if($aafull=~/partial3/); # 
        $tblinfo->{'aaqual'} =~ s/$aafull0/$aafull/ if($aafull0 ne $aafull);
        #FIX: $aaqual =~ s/complete|partial3/partial/; 
        }
      $tblinfo->{'cdsold'}= $tblinfo->{'cdsoff'}; # for annotab?
      $tblinfo->{'cdsoff'}= "$cdsb-$cdse";
      $tlog.="cutcds=$cdsb-$cdse,oldcds=$lcdsb-$lcdse,"; 
    }
  } 
  
  my $nl= length($fa); 
  $tblinfo->{'trlen'}= $nl;
  my $nn1= $fa=~tr/N/N/; 
  my $lendelta= ($nl==$ol and $nn==0) ? $nl : "$nl; olen=$ol; nnn=$nn1/$nn;";
  $lendelta .= " cut=$ncut;" if($ncut>0); 
  $tlog.="vectrim=$vectrimw," if($vectrimw); 
  ##er2g: >PogonEG0028883t1        oid=sobeetlepogo1ak31loc13234t4 
  ## len=413; olen=436; cut=0; nnn=18/41;    cutcds=1-314,oldcds=2-337,vectrim=23,
  ## for annotab:  cdsoff = cutcds; oldcds,olen not needed;  add: nnn=xxx,vectrim=xxx
  my $annogaps= ($nl==$ol and $nn==0) ? "0" : "gaps=$nn1/$nn";
  $annogaps.= ",oldcds=".$tblinfo->{'cdsold'} if($tblinfo->{'cdsold'});
  $annogaps.= ",cut=$ncut" if($ncut>0); 
  $annogaps.= ",vectrim=$vectrimw" if($vectrimw); 
  $annogaps = "TOOSHORT=$nl,$annogaps" if($nl<$MINSIZE);
  
  my($cdsoff,$gname,$dbxref,$aaqual,$trlen,$protid,$lotag,$namepct)= 
      @{$tblinfo}{qw(cdsoff name dbxref aaqual trlen protid locustag namepct)};

#   if($nl<$MINSIZE) #? do AFTER annoth so table has full record of skips?

  if($annoth) { 
    # usable output table ; FIXMEd: add header at top
    # NOTE: add vectrim info, other?; will need to regen .aa, .cds from .fsa output to be accurate..
    # .. user choice: ncbi submit restrictions may not be desired.
    my $aname= $gname || "hypothetical protein"; # || "na";  # which ??
    print $annoth join("\t",qw(PublicID OrigID TrLen CDSoff AAqual TrGaps Dbxref Namealign Product_Name))."\n" if($itr==1);
    print $annoth join("\t",$pubid,$oid,$trlen,$cdsoff,$aaqual,$annogaps,$dbxref,$namepct,$aname)."\n";
  }
  
  if($nl<$MINSIZE) { #? do AFTER annoth so table has full record of skips?
    ## ALSO retain shorties w/ uniq name **
    my $uniqname=0; 
    if($gname) { 
      my $allids= $namedgenes{$gname}; 
      $uniqname=1 if($allids eq "$oid,");
    }
    my $logt="#$def\toid=$oid\tlen=$lendelta";
    $logt.= "\tname=$gname" if($gname); $logt.= "\t$tlog" if($tlog);
    if($uniqname) {
      loggit(0,"PROBLEM: keep unique name but too short:$nl",$logt);
    } else {
      loggit(0,"ERR: skip too short:$nl",$logt);
      return (0,0,"ERR: skip too short:$nl");
    }
  }
        
  map{ $_="" if($_ eq "na"); } ($lotag,$gname,$dbxref,$aaqual);  

  $aafull ||= $aaqual;
  if($aafull and $aafull =~ /partial/) {
    $gname .= " fragment" if($gname);
#     $cdsoff =~ s/^(\d)/<$1/   unless($aafull =~ /partial3/);
#     $cdsoff =~ s/\-(\d)/>$1/  unless($aafull =~ /partial5/);
#     ##? drop cds if not complete? else have to deal w/ partial cds: "<cdsb >cdse CDS"
#     # $cdsoff=$lotag=""; #? unless($aafull =~ /complete/); #?? drop CDS or not
  }

  if($cdsoff =~ m/\d+\-/) { # allow for <123->456
    $ncdsout++; #? always same count as ntrout ?
    print $tblh ">Features\t$pubid\t$GBPROID\n\n"; #?? is this used now
    # >Features       Thecc1ER_sopcsc10rk23loc140t1     cacao11evigene_20120827
    
    my $cdsline= putCDSloc($cdsoff,$aafull,$cdsphase);
    print $tblh $cdsline;
    # $cdsoff =~ s/\-/\t/; print $tblh "$cdsoff\tCDS\n";
    
    print $tblh "\t\t\tprotein_id\t$protid\n" if($protid);
    print $tblh "\t\t\tlocus_tag\t$lotag\n" if($lotag);
    (my $npct=$namepct) =~ s/,.*//;
    print $tblh "\t\t\tnote\talignment:blastp is $npct\n" if($gname and $npct>0);
    print $tblh "\t\t\tnote\toriginal_id:$oid\n" if($oid);
    my $pname= $gname || "hypothetical protein";
    print $tblh "\t\t\tproduct\t$pname\n";
#         # ^^ FIXME: SEQ_FEAT.MissingCDSproduct; need some name

    if($dbxref =~ /\w/) {
      foreach my $dx (split",",$dbxref) { 
        # DBXREF_RECODE fix for ncbi's  dbx: restrictions
        $dx =~ s/^UniRef/SwissProt:UniRef/; ## UniProtKB:UniRef/; # or /TrEMBL:UniRef/;  SwissProt
        $dx =~ s/UniProt:/SwissProt:/; ## UniProtKB:/; # need list to check/replace as per 
        # is CDD:id ok? accepted by tbl2asn ..
        # evigene2genbanktbl.pl sub reformatval()
        #  my($d)= m/^(\w+):/; if( $DBXREF_RECODE{$d} ) { $d= $DBXREF_RECODE{$d}; }
        print $tblh "\t\t\tdb_xref\t$dx\n" if($dx=~/\w/ and $dx ne "na"); 
      }
    }
    print $tblh "\n";
  }
  
  $fa =~ s/(.{60})/$1\n/g; 
  print $outh ">$def\n$fa\n";  $ntrout++;
  loggit(0,">$def\toid=$oid\tlen=$lendelta\t$tlog"); # \tgene=$gn
  
  ## Maybe add option regen .aa at least from updated fa, cdsoff
  # if($outaah) { # == $REGEN_AA 
  #   my($XXaalen,$XXpcds,$XXcompl,$XXorflen,$fahead,$faprot) = cdna_proteins:translate1($fa,$cdsb,$cdse,AS_FASTA=1,...);
  #   print $outaah ">",$def," $fahead\n",$faprot,"\n";
  # }
  
  return ($ntrout,$ncdsout,"OK"); # no, return ($ntrout,$ncdsout)
}

#  ## FIXME here for new geneinfo
#   if($GENEINFO_VERS == 2) {
# 	  #	trid	geneid	mapqual	aaqual	trlen	or	cdsoff	locustag	product	dbxref
# 	  ## do revcomp and rev-offset for cdsor=-
#     $rinfo= $geneinfo{$oid} || [];
#     my($td,$gd,$mapqual,$aaqual,$trlen,$cdsor,$cdsoff,$lotag,$nam1,$dbxref)= @$rinfo;
#     $gn=$gd; $gname=$nam1;
#     
#     ($cdsb,$cdse)= split/[-]/,$cdsoff; # this SHOULD be offset in non-rev tr; any problems?
#     if($cdsor eq "-") {
#       my $farev= revcomp($fa); 
#       $fa= $farev; #? save orig
#       # my($cb,$ce)= split/[-]/,$cdsoff; # this SHOULD be offset in non-rev tr; any problems?
#       ## ** FIXME dammmit; if cb > ce, then this is offset BEFORE rev; if cb < ce is offset AFTER rev **
#       my($cbr,$cer)=(0,0);
#       if($cdsb > $cdse) {  ## this is new vers
#         $cbr= 1+$trlen - $cdsb;
#         $cer= 1+$trlen - $cdse;
#       } else { # DAMMIT mess w/ cdna_org version changed offs= meaning; this is old vers
#         $cer= $cdse; ## 1+$trlen - $cdsb;
#         $cbr= $cdsb; ## 1+$trlen - $cdse;
#       }
#       ($cdsb,$cdse)=($cbr,$cer);    
#       $tlog.="revcds:$cdsb-$cdse,";  
#       $cdsoff="$cdsb-$cdse"; # update !!!
#       $rinfo->[6]= $cdsoff;
#       $rinfo->[5]= $cdsor= "+";
#     }
#   } 
#...................
#  ## FIXME here for new geneinfo
#   if($GENEINFO_VERS == 2 and $tblh) {
# 	  #	trid	geneid	mapqual	aaqual	trlen	or	cdsoff	locustag	product	dbxref
# 	  ## do revcomp and rev-offset for cdsor=-
# 	  # use from above# my $rinfo= $geneinfo{$oid} || [];
#     my($td,$gd,$mapqual,$aaqual,$trlen,$cdsor,$cdsoff,$lotag,$nam,$dbxref)= @$rinfo;
# 
#     $tlog.="mapq=$mapqual,";
#     map{ $_="" if($_ eq "na"); } ($gd,$lotag,$nam,$dbxref);
#     # $gd="" if($gd eq "na"); # not all have this
#     # $lotag="" if($lotag eq "na"); # not all have this
#     my $protid= $id; # $td; # $gd ? gbasn requires protein_id here, use genes.tbl id? or tr id?
#     $protid=~s/t(\d+)$/p$1/; 
#     $protid= $GDB_PREFIX.$protid if($protid);
#     #  protein_id      gnl|CacaoGD|Thecc1EG016762p1
#     
#     # require aaqual =~ /complete/ and $mapqual =~ /mapfull:/ and $mapqual =~ /aaeq|aasim/
#     # FIXME: for genome-less asmrna: no mapqual, flag as 'nomap' or such
#     unless($aaqual =~ /complete/ and $mapqual =~ /mapfull/ and $mapqual =~ /aaeq|aasim/) {
#       $dbxref=""; # dont claim this annot unless hi qual gene match
#       ## FIXME: SEQ_FEAT.MissingCDSproduct
#       if($nam =~ /\w/ and $nam ne "na") { $nam .= " fragment"; } 
#       else { $nam="hypothetical protein"; }
#       #? drop cds if not complete? else have to deal w/ partial cds: "<cdsb >cdse CDS"
#       $cdsoff=$lotag="" unless($aaqual =~ /complete/); #?? drop CDS or not
#     }
#     
#     ## require cdsoff, skip if have only locustag
#     if($cdsoff =~ m/\d+\-\d+/ and $lotag) { ##  or $lotag
#       print $tblh ">Features\t$id\tcacao11evigene_20120827\n\n";
#       # >Features       Thecc1ER_sopcsc10rk23loc140t1     cacao11evigene_20120827
#       
#       if(0 and $lotag) {
#         print $tblh "1\t$trlen\tgene\n";
#         print $tblh "\t\t\tlocus_tag\t$lotag\n"; #? no gene row? put in CDS?
#         print $tblh "\n";
#         }
#       if($cdsoff =~ m/\d+\-\d+/) {
#         $cdsoff =~ s/\-/\t/;
#         print $tblh "$cdsoff\tCDS\n";
#         print $tblh "\t\t\tprotein_id\t$protid\n" if($protid);
#         print $tblh "\t\t\tlocus_tag\t$lotag\n" if($lotag);
#         print $tblh "\t\t\tproduct\t$gname\n" if($gname);
#         # ^^ FIXME: SEQ_FEAT.MissingCDSproduct; need some nam: hypothetical??? or orig-name + similarto/
#         if($dbxref =~ /\w/) {
#         foreach my $dx (split",",$dbxref) { print $tblh "\t\t\tdb_xref\t$dx\n" if($dx=~/\w/ and $dx ne "na"); }
#         }
#         print $tblh "\n";
#       }
#     }
#   }



sub vecscreen
{
  my($cdnaseq,$vectab)=@_;
  
  $vectab= makename($cdnaseq,".vector.tab") unless($vectab);
  return($vectab) if( -s $vectab); # or dryrun ..
  
  our($id,$vb,$ve,$ty,$vd,$outh,$inh);  
  sub putv { our($id,$vb,$ve,$ty,$vd,$outh,$inh); 
    print $outh join("\t",$id,$vb,$ve,$ty,$vd)."\n" if($id and $ty and not $ty=~/Weak/); 
  }
  
  #  $ncbi/bin/vecscreen
  ## FIXME: can we use ncbic++ instead? output not same as c-vecscreen... need to check curr ncbi source
  my $univecdb="";
  (my $ncbid=$APPvecscreen) =~ s,/vecscreen,/..,; ## want option for db UniVec path
  if( -f "$ncbid/data/UniVec.nsq") {
    $univecdb= "$ncbid/data/UniVec"; 
  } else {
    loggit(1,"ERR: $APPvecscreen missing ../data/UniVec.nsq"); return; 
  }
  
  ## lots of this warn: [vecscreen] WARNING:  [000.000]  Blast: No valid letters to be indexed on context 0
  my $vectmp= makename($cdnaseq,".vecscreen.tmp");
  my $veclog= makename($cdnaseq,".vecscreen.log");
  runcmd("$APPvecscreen -i $cdnaseq -d $univecdb -f3 -o $vectmp 2> $veclog");
  my $ok= open($inh,$vectmp); 
  
  ## use runcmd() here? write vec.tmpfile and read that instead of pipe.
  # loggit(0,"$APPvecscreen -i $cdnaseq -d $univecdb -f3");   
  # my $ok= open($inh,"$APPvecscreen -i $cdnaseq -d $univecdb -f3 |");
  
  $ok= open($outh,'>',$vectab) if($ok);
  unless($ok) { loggit(1,"ERR: $APPvecscreen -i $cdnaseq -d $univecdb TO $vectab"); return; }
  my($b,$e,$so);
  while(<$inh>) {
    chomp; 
    if(/^>/) { putv() if($id); ($id)=m/>Vector (\S+)/; ($vd)=m/Database: (\S+)/; $vb=$ve=$so=$ty=0; } 
    elsif(/^No hits/) { $id=0; } 
    elsif(/^\w/) {  
      if(/match/) { $ty=$_; } elsif(/^Suspect origin/) { $so=1; } 
      elsif(/^(\d+)\s+(\d+)/) { ($b,$e)=($1,$2); 
        if($so and $ve) { $vb=$b if($b<$vb); $ve=$e if($e>$ve); }  
        elsif($ty) { ($vb,$ve)=($b,$e); }  } 
      }
  } 
  putv() if($id); close($outh); close($inh);
  return($vectab);
}

sub revcomp {
  my ($seq) = @_;
  my $reversed_seq = reverse ($seq);
  $reversed_seq =~ tr/ACGTacgtyrkm/TGCAtgcarymk/;
  return ($reversed_seq);
}

sub findapp
{
  my($aname)=@_;
  my $app="";
  $app=$ENV{uc($aname)} if(not $app and $ENV{uc($aname)});  
  $app=$ENV{$aname} if(not $app and $ENV{$aname});  
  $app=`which $aname` unless($app); 
  chomp($app);
  ## #tr2aacds: app=blastn, path=no blastn in 
  my $dol=0; if(not $app or $app =~ /^no $aname/) { 
    $app="echo MISSING_$aname"; $dol=($dryrun||$debug)? LOG_WARN : LOG_DIE; 
    }
  loggit( $dol, "app=$aname, path=$app");
  return($app);
}

sub findevigeneapp
{
  my($aname)=@_;
  my $app= $aname;
  # my $EVIGENES="$FindBin::Bin/.."; # ok?
  # my $APPcdnabest= findevigeneapp("$EVIGENES/cdna_bestorf.pl"); # allow ENV/path substitutions?
  $app="$EVIGENES/$aname" unless(-x $app);
  my $dol=0; 
  unless( -x $app) { 
    $app="echo MISSING_$aname"; 
    $dol=($dryrun||$debug)? LOG_WARN : LOG_DIE; 
    }
  loggit( $dol, "evigeneapp=$aname, path=$app");
  return($app);
}



sub runcmd
{
  my @cmd= @_;
  loggit( ($dryrun) ? 1 : 0,"CMD=",@cmd);  
  ## fail if $cmd[0] =~ /MISSING_/
  my $err= ($dryrun) ? 0 : system(@cmd);  
  if($err) { loggit(1,"ERR=$err ",$cmd[0]); } # ..
  return $err;
}

sub makename
{
  my($infile,$osuf,$insuf)=@_;
  $insuf ||= 'aa|blast|cdna|cds|tr|trclass|fasta|fa';  ## fixme need insuf: tr|fasta|fa
  my $outfile= $infile; $outfile =~ s/\.gz$//;
  $outfile =~ s,\.($insuf)[^\/\s]*$,,; 
  $outfile.= $osuf if($osuf); 
  $outfile.= "_out" if($outfile eq $infile);
  return $outfile;
}


__END__
 
 
=item tsasub11 cacao

  evigene=/Users/gilbertd/Desktop/dspp-work/genes2/evigene
  ncbin=/Users/gilbertd/Desktop/dspp-work/genomesoft//ncbi/bin
  nwbsra='[SRA=SRR531454; SRR531455; SRA058778; SRA058779]'
  vel3sra='[SRA=SRA058777; SRA058780; SRA058781; SRA058782; SRR531454; SRR531455; SRA058778; SRA058779]'
  rnasra='[SRA=SRA058777; SRA058780; SRA058781; SRA058782]'
    ## redo asmrna2ncbitsa.pl -GAPSOK
    
  cd tsasub11
  
  for pt in cacao[345]{cuf[28],nwb,vel,sop,tri}; do {
  
  if [ -f $pt/TCM01.tsa_rasm.$pt.tbl ]; then continue; fi
  echo asmrna2ncbitsa $pt/TCM01.tsa_rasm.$pt
  
  $evigene/scripts/rnaseq/asmrna2ncbitsa.pl -GAPSOK -idpre Thecc1ER_ \
  -cdna ../tr5parts/pub3ig.$pt.tab4g.tr.gz -vec ../tr5parts/pub3ig.trasm.tab4g.vector.tab \
  -geneinfo ../tr5parts/pub3ig.trasm.tab4g.geneinfo1.tab  -log tr4g.$pt.log \
  -out $pt/TCM01.tsa_rasm.$pt.fsa -tbl $pt/TCM01.tsa_rasm.$pt.tbl
  
  } done
  
  
    # was -a r10u; unknown gap len; YES, r10k got rid of warning
  
  for pt in cacao[345]{cuf[28],nwb,vel,sop,tri}; do {
   if [ $pt == cacao3nwb ]; then sra=$nwbsra; 
   elif [ $pt == cacao3vel ]; then sra=$vel3sra; 
   else sra=$rnasra; fi
   
   if [ -f $pt/TCM01.tsa_rasm.$pt.sqn ]; then continue; fi
  
   echo $pt : $sra
   
   $ncbin/tbl2asn -p $pt/ -Z $pt/$pt-discrep.log  -w tsamethods.$pt.cmt \
   -Y tsadesc.$pt.cmt -t cacao3i_tsasubmit.sbt \
   -a r10k  -l paired-ends -Vtb -Mt -XE \
   -j "$sra [bioproject=PRJNA51633] [moltype=mRNA] [tech=TSA] [organism=Theobroma cacao]"
  
  } done
  
=cut
 