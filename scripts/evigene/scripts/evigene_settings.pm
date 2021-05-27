# evigene_settings.pm projectinfo or settings?

# package evigene_settings;
package main;

=item about package

  
=cut

use constant evigene_settings_VERSION => '2019.05.12';  
use constant FIXME => 1;

use FindBin;
use lib ("$FindBin::Bin"); # assume evigene/scripts/ layout: main, prot/, rnaseq/, ...

use strict;
use warnings;
use File::Basename qw(basename dirname fileparse);

use cdna_evigenesub;  
use evigene_pubsets; 

# pipeline global vars
our ($IDPREFIX, $ORGANISM, $DATE, $EGLOG, $EGAPP, $DEBUG, $TSADESC, $TSAOPTS);  
our ($BioProject, $sraids, $sratable, @sraids);
our (%settings); # rename for shared global?

our $DEFAULTidpre= 'NonameEVm';  ## this seems to fix -idprefix ignored bug
$IDPREFIX= $ENV{idprefix} || $DEFAULTidpre; ## "evgr"; #  opt
$ORGANISM= $ENV{organism} ||$ENV{species} || "Noname";
$DATE= `date '+%Y%m%d'`; chomp($DATE); # default is today; use perl func?

our %DEFAULT_SETTINGS= ( 
  IDPREFIX=>$DEFAULTidpre, 
  DATE=>$DATE, 
  organism => $ORGANISM, 
  #?? sraids => $sraids, 
  BioProject => $BioProject, #fixme bioproj lost
  # assemblers => 'Velvet/Oases; idba_trans; SOAPDenovoTrans; Trinity;',
  trclass => '', mrna => '', genenames=>'', 
  # genome => 'genome/chromosome_assembly.fasta.gz',
  ); 

#---------------------------------------------------  

sub add_setting {
  my($key,$val)=@_;
  if(my $ov= $settings{$key}) { 
    return if($ov =~ m/\b$val\b/);
    $val="$ov;$val";
  }
  $settings{$key}= $val; # return $val; ?
}

sub set_setting { my($key,$val)=@_; $settings{$key}= $val; }
sub set_newsetting { my($key,$val)=@_; $settings{$key}= $val unless($settings{$key}); }

sub get_setting {
  my($key)=@_;  my $v= $settings{$key}||undef; # add? || $ENV{$key}
  return $v;
}

=item do_settings

  do_settings = readwrite_settings 
  =~ cdna_evigenesub.pm: sub evigene_config() ?? use that?
  -- split into auto-updated.settings and user-supplied settings? ie dont replace latter
  
=cut

sub do_settings {  
  my($action,$pathname)= @_;
  ## write these to work dir; reread next go
  ## action == 'log|save|restore' ; restore called AFTER read new options, shouldnt replace


  my $PRES=$EGLOG; # == 's2g';
  # my $trpname= makename($pathname,".sra2genes.info"); # FIXME, CALLER sets pname.info
  my $trpname= makename($pathname,".$EGAPP.info"); 
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
		sraids => $sraids, sratable => $sratable, 
		organism => $ORGANISM, BioProject => $BioProject,  
		runpath => $runpath,
		TSADESC => $TSADESC, TSAOPTS => $TSAOPTS,
		# MINSIZE=>$MINSIZE, MAXGAP=>$MAXGAP,
		# trclass => $trclass, mrna => $cdnaseq, genenames=>$genenames,
		);
	
  if($action =~ /restore|read/ and -f $trpname) {
    open(my $inh, $trpname); # or loggit(warn ..);
    while(<$inh>) { chomp; if(s/^$PRES.//) { 
      my($k,$v)=split /\s*=\s*/,$_,2; 
      $v=~s/;*$//; $v=~s/\s*\#.*$//; # should I chomp trailing '#' comments? 
      my ($v1)= ($v=~/;/)? (split";",$v)[0] : $v; # ugh: organism=Bob_white;Bob_white_black;.. from sra.csv mixtures
      my $ov= $mysettings{$k}; 
      unless($ov and $ov ne $DEFAULT_SETTINGS{$k}) { 
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
    @sraids= grep /^[A-Z]+\d+/, split /\W+/, $sraids; # NOW GLOBAL; should be SRRnnnn ERRnnnn DRRnnn 

  } else { # action==save; copy new single vars to settings..
    for my $ks (sort keys %mysettings) { 
      my $ov= $mysettings{$ks};  
      $settings{$ks}= $ov if($ov and $ov ne $DEFAULT_SETTINGS{$ks});    
      # set_setting($ks,$ov) if($ov and $ov ne $DEFAULT_SETTINGS{$ks});

      } ;
  }

  my $settings= join "\n", map{ "$PRES.$_=".$settings{$_} } sort keys %settings;
  if($action =~ /log/) { loggit(0, "$EGAPP.info:\n$settings");  }
  if($action =~ /save/) { open(my $outh, '>', $trpname); print $outh $settings,"\n"; close($outh); }
}


=item SRA data queries

  http://www.ncbi.nlm.nih.gov/sra
  query=
  (("biomol transcript"[Properties]) AND "platform illumina"[Properties]) AND "library layout paired"[Properties] 
  * 2016.06: changed vocab: "biomol transcript"[Properties]
  
  RNA illumina paired : Public(164,016) 
  (("biomol transcript"[Properties]) AND "platform illumina"[Properties]) AND "library layout paired"[Properties] 
  RNA pacbio smrt     : Public(303) 
  ("biomol transcript"[Properties]) AND "platform pacbio smrt"[Properties]

=item parse_sra_result_cvs

  add: parse sra_result.csv if exists, for species, sraids .. other metadata
  evgr2tsa/litova1all3f/
    whiteshrimp_sra_result.csv # rename litova1all3.sra_result.csv

  "Experiment Accession","Experiment Title","Organism Name","Instrument",
    "Submitter","Study Accession","Study Title","Sample Accession","Sample Title",
    "Total Size, Mb","Total RUNs","Total Spots","Total Bases","FTP Path to Experiment",
    "Library Name","Library Strategy","Library Source","Library Selection"
    
  "SRX098246","Litopenaeus vannamei  transcriptome","Litopenaeus vannamei","Illumina HiSeq 2000","BGI",
    "SRP008317","BGI Litopenaeus vannamei transcriptome sequencing","SRS265043","Pacific white shrimp",
    "1515.37","1","13697473","2465545140","ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX098/SRX098246",
    "Ex-zai-2_l1","RNA-Seq","TRANSCRIPTOMIC","cDNA"

=item sra_cvs 2017 format

Run,ReleaseDate,LoadDate,spots,bases,spots_with_mates,avgLength,size_MB,AssemblyName,download_path,
  Experiment,LibraryName,LibraryStrategy,LibrarySelection,LibrarySource,LibraryLayout,InsertSize,InsertDev,
  Platform,Model,SRAStudy,BioProject,Study_Pubmed_id,ProjectID,Sample,BioSample,SampleType,
  TaxID,ScientificName,SampleName,g1k_pop_code,source,g1k_analysis_group,Subject_ID,
  Sex,Disease,Tumor,Affection_Status,Analyte_Type,Histological_Type,Body_Site,CenterName,Submission,
  dbgap_study_accession,Consent,RunHash,ReadHash

SRR3247180,2016-10-04,2016-06-13,26018777,6504694250,26018777,250,3026,,
  https://sra-download.ncbi.nlm.nih.gov/srapub/SRR3247180,SRX1645097,,RNA-Seq,
  PCR,TRANSCRIPTOMIC,PAIRED,0,0,ILLUMINA,Illumina HiSeq 2500,
  SRP072048,PRJNA315720,2,315720,SRS1350203,SAMN04569208,simple,743457,
  Daphnia similoides,Parthenogenetic female,,,,,female,,no,,,,,HUAIBEI NORMAL UNIVERSITY,
  SRA388003,,public,49CAA8290AC1FD00A6A5F0C61EDD392E,1EFA3DFE6E3818211435373BFC191AE2
...

  Run,ReleaseDate,LoadDate,spots,bases,spots_with_mates,avgLength,size_MB,AssemblyName,download_path,Experiment,LibraryName,LibraryStrategy,LibrarySelection,LibrarySource,LibraryLayout,InsertSize,InsertDev,Platform,Model,SRAStudy,BioProject,Study_Pubmed_id,ProjectID,Sample,BioSample,SampleType,TaxID,ScientificName,SampleName,g1k_pop_code,source,g1k_analysis_group,Subject_ID,Sex,Disease,Tumor,Affection_Status,Analyte_Type,Histological_Type,Body_Site,CenterName,Submission,dbgap_study_accession,Consent,RunHash,ReadHash
  ERR1249549,2017-01-20,2017-01-20,0,0,0,0,0,,https://sra-download.ncbi.nlm.nih.gov/srapub/ERR1249549,ERX1321488,Es1,WGS,cDNA,TRANSCRIPTOMIC,PAIRED,500,0,ILLUMINA,Illumina MiSeq,ERP014147,,,0,ERS1052381,SAMEA3865247,simple,1793945,Eremophila serrulata,Es,,,,,,,no,,,,,CEBITEC,ERA560968,,public,,
  SRR1558510,2015-06-09,2014-09-08,53574,7687869,0,143,3,,https://sra-download.ncbi.nlm.nih.gov/srapub/SRR1558510,SRX687195,,RNA-Seq,cDNA,TRANSCRIPTOMIC,PAIRED,0,0,ILLUMINA,Illumina MiSeq,SRP045791,PRJNA259519,2,259519,SRS689854,SAMN03008982,simple,10090,Mus musculus,GSM1488161,,,,,,,no,,,,,GEO,SRA180334,,public,19848BFFA4A7DAF75B80D1D108FBEC6B,F299F1EBB22A74DF5483D1FAE8BDBFDE
  ..
  case daphnia sim:
  SRR3247516,2016-10-04,2016-06-13,31197702,7799425500,31197702,250,3645,,https://sra-download.ncbi.nlm.nih.gov/srapub/SRR3247516,SRX1645182,,RNA-Seq,PCR,TRANSCRIPTOMIC,PAIRED,0,0,ILLUMINA,Illumina HiSeq 2500,SRP072048,PRJNA315720,2,315720,SRS1350211,SAMN04569209,simple,743457,Daphnia similoides,Sexual female,,,,,,,no,,,,,HUAIBEI NORMAL UNIVERSITY,SRA388003,,public,BDB5ADFE4321990CD6C79171FB645086,DF05C2947FBF60CBD911B697C274AC67

  * use LibraryLayout = PAIRED vs SINGLE, and spots_with_mates vs spots
  
=cut


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

  make_IDPREFIX_4org() if($sorg); # rep w/  pm make_IDPREFIX

  ## also use sradatah to edit template evgr_tsamethods.cmt, evgr_tsadesc.cmt, evgr_tsasubmit.sbt ??
  ## rewrite template .cmt, .sbt unless updated already.

  ## ALSO use 'FTP Path to Experiment' to turn SRX into SRR for picky ncbi tsa : wget or curl calls?
  my $HAVEsrrids= ($sraids =~ /[DES]RR\d/ and $sraids ne $defsra)?1:0; ## Need to save to info/config file
  unless(FIXME or $HAVEsrrids) {
    my @urls;
    @urls= grep /http|ftp/, split/;/, $sradata{"download_path"} || $sradata{"FTP Path to Experiment"};  
    #o @ftps= grep /ftp/, split/;/, $sradata{"FTP Path to Experiment"};  
    my @srr= sra_ftp2srr(@urls);
    if(@srr>0) {
      $sraids= join(";",@srr);
      $sradata{'SRAids'}= $sraids;
      loggit(1,"sra_id=",$sraids);
      }
  }
   
  return($ngot, \%sradata, $sraids); # sids ?
}

sub make_IDPREFIX_4org # see evigene_pubsets.pm:make_IDPREFIX()
{
  if($ORGANISM and $ORGANISM =~ m/\w\w/ and $ORGANISM ne "Noname") {
    my $shortorg="";
    my($gen,$spp)= split /[_\W]/, $ORGANISM,2;
    if($spp and length($gen)>1) { $shortorg= ucfirst(substr($gen,0,3)) . lc(substr($spp,0,3)) ; }
    else { $shortorg= ucfirst(substr($ORGANISM,0,6)); }
    $IDPREFIX= $shortorg.'EVm' if($IDPREFIX eq $DEFAULTidpre);
    set_newsetting('oname',$shortorg);
    return 1;
  }
  return 0;
}

sub get_srainfo {
  my($sracvs,$sraidlist)= @_;
	#o# my($trpath,$trname)= @_;
	
	unless($sracvs) {
	  # my($sracvs)= getFileset('.','sra_result.csv$|sra.csv$');  
	  # my $sracvs="$trpath/$trname.sra_result.csv";    
	  # $sracvs="$trpath/sra_result.csv" unless(-f $sracvs);  
  }
  
	my @SRAK=();
	my($nsra,$sradatah,$gotsraids)=(0,0,0);
	if(-f $sracvs) {
	  ($nsra,$sradatah,$gotsraids)= parse_sra_result_cvs($sracvs);
	} else { 
    $sradatah= {};
    # make dummy sra.cvs for other components?? pubset2submit.pl wants it
    $sradatah->{cvsformat}= 0;
	  make_IDPREFIX_4org(); # if $ORGANISM
	  if($sraidlist) {
      $sraids= $sraidlist; #? or always use sids?
      @sraids= grep /^[A-Z]+\d+/, split /\W+/, $sraids; # NOW GLOBAL; should be SRRnnnn ERRnnnn list
      $nsra= @sraids;
     }
	}
	
	# $gotsraids == global $sraids
	
	loggit(0,"sra_result from",$sracvs,"nresult=",$nsra);
  if($nsra>0) {
	  if ($sradatah->{cvsformat} == 1) {
      @SRAK=("Experiment Accession","Organism Name","Instrument",
              "Submitter", "Total Size, Mb","Total Spots");
    } else {
      @SRAK=("Run","ScientificName","Platform", # Instrument == Platform,Model
             "CenterName", "size_MB","spots", "BioProject");
    }
    my @v= map{ $_ .'='. $sradatah->{$_}.';'; } @SRAK; 
    loggit(0,"sra_info:",@v); 
    if(0 && $DEBUG){ map{ my $v=$sradatah->{$_}; loggit(0,"sra.",$_,'=',$v) if($v=~/\w/); } (sort keys %$sradatah);  }
  }
  
	if($sradatah->{cvsformat} == 1) {
  	@SRAK= ("Assemblers", "Instrument", "Total Size, Mb","Total Spots",'Total Assemblies');
	} else {
  	@SRAK= ("Assemblers", "Platform", "size_MB","spots","Total Assemblies", "BioProject");
	}
	for my $ks (@SRAK) {
  	if($settings{$ks} and not $sradatah->{$ks}) { $sradatah->{$ks} = $settings{$ks}; }
  	elsif($sradatah->{$ks}) { $settings{$ks} = $sradatah->{$ks}; }
	}
	
	return($sradatah,$nsra);
}


=item sra_ftp2srr from FTPpath

  curl  'ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX124/SRX124618/'
  dr-xr-xr-x 1073741824 ftp      anonymous        0 Aug 12  2012 SRR424344
  
  >> best:
  curl -s -l 'ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX160/SRX160070/'
  SRR521835
  SRR521944
  
  wget -A 'sra' -r -l 2 -nv -nd -nH -np --spider \
   'ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX124/SRX124618/'
  14:25:43 URL: ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX124/SRX124618/ [74] -> ".listing" [1]
  14:25:44 URL: ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX124/SRX124618/SRR424344/ [74] -> ".listing" [1]
  14:25:44 URL: ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX124/SRX124618/SRR424344/SRR424344.sra [0] -> "SRR424344.sra" [1]

=cut

sub sra_ftp2srr {   # FIXME
  my(@ftps)= @_;  return () unless(@ftps);
  my $APPcurl= findapp('curl',1); return () if($APPcurl =~ /MISSING/);
  # my $APPwget= findapp('wget'); return () if($APPwget =~ /MISSING/);
  ## dang new sra.csv lacks ftp: url, add it:
  my $baseu="ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra";
  my @srrs=(); 
  
  foreach my $ftppath (grep /^ftp:/, @ftps) {
  	my $cmd="$APPcurl -s -l $ftppath/"; 
  	loggit( ($dryrun) ? 1 : 0,"CMD=",$cmd);  
    my $srrs= `$cmd`; # or http: ?? ## runcmd() instead ?
    push @srrs, grep /SRR/, map{ m/(SRR\w+)/; $1; } split " ",$srrs;
    }
  return @srrs;
}

sub sraget { # FIXME
  my $srrh=undef;

  # mirror fetch from ncbi sra_result.csv table listing ftp address of data
  # # $mr ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX151/SRX151669
  ## Need option for fork/background.  Ncbi/local gets sick w/ too many calls at once.. fails some.
  ## 2014.11: ncbi ftp allowing only 2-same time; add sleep(5); wait?
  ## dang new sra.csv lacks ftp: url, add it:

  my $baseu="ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra";
  # $DEBUG=$ENV{debug}||0;
  my $dofork=$ENV{fork}||0;  # == wget -b
  # my $SNOOZE=$dofork;
  my $ndone=0;
  #my $mr='lftp -c mirror ';
  my $mr='wget -m -nv -np -nH --cut-dirs=7 ';
  $mr.=" -b " if($dofork); # wget -b will fork to background..
  my $pid=0;
  
  while(<$srrh>) {
    next if(/^#/);  
    chomp; my @v= map{ s/^"//; s/"$//; $_; } split",";
    my ($u)= grep /ftp:/, @v;
    unless($u) { my $sx=$v[0]; if($sx=~m/^[A-Z]{3}\d{6}/) { 
      my($sa)=substr($sx,0,3); my($sb)=substr($sx,0,6); $u="$baseu/$sa/$sb/$sx"; } 
    }
    my $ok=0;
    if($u) { 
      $ok= ($DEBUG)? "debug" : system("$mr $u "); 
     }
    warn "#fork=$ok,$pid $mr $u\n"; 
    $ndone++;
    sleep(5) if($dofork and $ndone % 2 == 0);
  }
}

