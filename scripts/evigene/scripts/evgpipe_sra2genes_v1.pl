#!/usr/bin/env perl
# evgpipe_sra2genes.pl

=item about

  ombnibus pipeline to reconstruct genes from SRA data (RNA-seq), with evigene methods
  outlined in evigene/docs/EvigeneAssemblyFromSRAdata.txt
  take 1: 2017.nov.04
  d.g.gilbert

=item method documents

  evigene/docs/evigene_rnapipe_methods1708.txt
  evigene/docs/perfect-mrna-assembly-2013jan.txt  
  evigene/docs/evgmrna2tsa_help.txt
  evigene/docs/evigene_goals2015.txt
  evigene/docs/evigene_goals2015b.txt

=item FIXME

  x Fixed: BUG STEP3_selectrna
  - add opt: -NCPUforClusterScripts=24  -NCPUforPipe=8 ... separate amounts
  - add opt: -MAXMEM
    .. add calc for sensible NCPU/asm script, dont want 100 cpu x 50 GB mem calls
  - change opt: maxmb|MAX_SIZE_MB >> -maxdataMB | MAX_DATASIZE_MB ..

  - need ref data options, path-to, methods : 
    S8. refblast, S9. names (ref.aa names)
    S7..S9: vecscreen UniVec, contamscreen NCBI contamseqs/rRNA, 
    S8b: NCBI CDD or other cons domain data, rpsblast or hmmer

  
=item pipe design
  
  goal is automated RNA-data to annotated, publishable gene set
  basic fetch rna-data, prepare, assemble several ways, combine assemblies, evigene reduce,
    reference protein/domain blast annotate genes, contaminant/vector screen, public file set
  
  input requirements
    a. SRA data ID/table of metadata (species, contents)
  input requirements/options
    b. reference species protein set (fetch? have ready? option?)   
    c. conserved prot domains (CDD)
    d. reference contaminant data (vectors, rRNA db, dfam? transposons)
    
  serial step-wise pipeline
    -- each step should be run with checks
    -- design pipe to be called for each step, or run continuous waiting for each step to finish
    -- step run time vaguely known, variable (minutes, hours, day+)

  some steps should write script to run async, under various cluster batch control systems
    then restart this pipe at next step (checking results of async step)
          
  options will be added
  
  See also evigene/scripts/evgmrna2tsa2.pl (pipe script that is now too messy for customers..)
  
=item components

  1. get data
    sraget.pl : small script, enclose
    = wget to ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/INSERTDATAPATH.sra
    env fork=1 ./sraget.pl daphsim16huau_srarna.csv  >& get1.log &
  * Merge 1+2, sratoolkit/fastq-dump now does web-data-fetch (http?), given SRRid
  
  2. sra to fasta, uses NCBI sratoolkit/fastq-dump
    env sra=SRX*/SRR*/SRR*.sra datad=`pwd` prog=./runsra2fa.sh sbatch ../srun_shared.sh

  3. subset for data size; max data size for basic assembly : 15-20GB for given maxmem=120 GB
    3a. option: clean/trim rna.fa
    3b. option: digital normalize rna.fa
  
  4.  run assemblers, with kmer size options, other opts
    4a. velvet/oases, ~10 kmer steps
    4b. idba_tran, ~10 kmer steps
    4c. soap_trans, ~10 kmer steps
    4d. trinity / other / user choices
    -- est. time at 12hr x 8cpu x 120 GB mem per assembler multi-kmer set (same for 1kmer trinity)
    
  5. post process assembly sets (trformat.pl)
    $evigene/scripts/rnaseq/trformat.pl -pre $nap -out $subd.tr -log -in $subd/vel*/transcripts.fa 
    .. etc per assembler/output set
    -- collect into evgrundir/subsets/
    
  6. quick qual assessment: cdna_bestorf > aastats per assembly, report
    $evigene/scripts/cdna_bestorf.pl -nostop -minaa=30 -aa -cdna $subd.tr.gz 
    $evigene/scripts/prot/aaqual.sh $subd.aa
    $evigene/scripts/prot/aastat.sh $subd.aa.qual >> assembly.aastat
    
  7. run evg over-assembly reduction, tr2aacds.pl
    env trset=$pt.tr datad=`pwd` prog=./runtr2cds.sh sbatch srun_comet.sh

  8. ref protein blastp x evg okayset
    env aaset=okayset/$pt.aa refaa=refset/$refaa ncpu=20 datad=`pwd` prog=./run_evgaablast.sh sbatch srun_comet.sh 

  9. namegenes from ref names, for annotation
    $evigene/scripts/prot/namegenes.pl -blast $aabltab -refnames refset/$refaa.names -out $pt.names
  
  add options:
    9b. vecscreen -db UniVec vectors == SEE evgmrna2tsa2.pl
    9c. blastn -db rRNA contam 
        ^^^ 9b,c extract from evgmrna2tsa2 into separate pipe-part script
    9d. rpsblastp -db CDD conserved domains
    9e. hmmscan -db dbfam transposons
    
  10. annotated publicset 
    env idprefix=$idp trclass=$pt.trclass names=$pt.names  species=$spp datad=`pwd` \
       prog=./run_evgmrna2tsa.sh sbatch srun_comet.sh
  
=item runme script

  #! /bin/bash
  ### env sratable=sraset.csv datad=`pwd` ncpu=16 qsub -q normal run_evgsra2genes.sh
  #PBS -N evgsra2genes
  #PBS -A PutAccountIdHere
  #PBS -l nodes=1:ppn=16,walltime=39:55:00
  #PBS -V
  
  if [ "X" = "X$ncpu" ]; then ncpu=8; fi
  if [ "X" = "X$maxmem" ]; then maxmem=64000; fi
  if [ "X" = "X$datad" ]; then echo "ERROR: missing datad=/path/to/data"; exit -1; fi
  if [ "X" = "X$sratable" ]; then echo "ERROR: missing sratable=/path/to/data"; exit -1; fi
  # opt name=testsra2evg
  
  # XSEDE .sdsc.edu
  if [ 1 = 1 ]; then
    bioapps=$HOME/bio
    evigenes=$bioapps/evigene/scripts
    # NOTE need current sratoolkit281 for web fetch by SRR id
    srabin=$bioapps/sratoolkit/sratoolkit281/bin
    #  velvet: fixme multi kmer binaries, bin4 = 151mer; bin2 = 99mer
    velobin=$bioapps/velvet1210/bin4  
    idbabin=$bioapps/idba/bin
    soapbin=$bioapps/soaptrans103 
    trinbin=$bioapps/trinity
    xnrbin=$bioapps/exonerate/bin
    cdhitbin=$bioapps/cdhit466/bin   
    ncbibin=$bioapps/ncbi/bin
  fi
  
  export PATH=$srabin:$velobin:$idbabin:$soapbin:$trinbin:$xnrbin:$cdhitbin:$ncbibin:$evigenes:$PATH
  
  evopts="-NCPU $ncpu -log -debug"
  if [ "X" != "X$name" ]; then evopts="$evopts -runname $name"; fi
  
  cd $datad
  echo $evigenes/evgpipe_sra2genes.pl $evopts -SRAtable $sratable
  $evigenes/evgpipe_sra2genes.pl $evopts -SRAtable $sratable


=item test settings

  srabin=$gs/sratools/sratoolkit.2.8.1-2-mac64/bin
  velobin=$gs/rnaseq/velvs/velo120/velbin1  # fixme multi kmer binaries
  idbabin=$gs/rnaseq/idba/idba-1.1.1/bin
  soapbin=$gs/rnaseq/soaptr/SOAPdenovo-Trans-r104/bin  # fixme binaries
  trinbin=$gs/rnaseq/trin12/trinityrnaseq_r2012-10-05/ # fixme binaries
  xnrbin=$gs/exonerate220/bin
  cdhitbin=$gs/cdhit17/bin   
  nbin=/bio/bio-grid/mb/ncbix/bin

  .. see findapp() for needed paths
  findapp('fastq-dump');  
  findapp('idba_tran', 1);
  findapp('velveth', 1); # revise velbin compile version names: vel{h,g},oases_NNNmer for max kmer
  findapp('SOAPdenovo-Trans-127mer', 1);
  findapp('Trinity', 1);
  findapp('fastanrdb', 1); 
  findapp('cd-hit-est', 1); 
  findapp('blastn', 1); # also vecscreen tbl2asn maybe
    
=cut


use constant VERSION => '2017.11.04'; # prelim
use constant FIXME => 1;

use FindBin;
use lib ("$FindBin::Bin"); # assume evigene/scripts/genes; layout: main, prot/, rnaseq/, ...

use strict;
use Getopt::Long;
use File::Basename qw(basename dirname);
use cdna_evigenesub;  
#maybe# use cdna_proteins;
#maybe# use protein_names;

our $EVIGENES="$FindBin::Bin";  
our $EGAPP='sra2genes';  
our $EGLOG='s2g';
our $dryrun=0; ## $DRYRUN ?
our $DEBUG= $ENV{debug}|| 0;

# global vars duplicated in %settings, DEFAULT_SETTINGS : do away with?
my($organism,$sraids,$BioProject)=("Noname",0,"");  ## SRR 346404; "SRR000000"
my $DEFAULTidpre= 'NonameEVm';  ## this seems to fix -idprefix ignored bug
my $IDPREFIX= $ENV{idprefix} || $DEFAULTidpre; ## "evgr"; #  opt
my $DATE=`date '+%Y%m%d'`; chomp($DATE); # default is today; use perl func?
my $MAX_SIZE_MB = $ENV{maxdatamb}||10000; # 10 Gb for testing

##default/opt $settings{'assemblers'} == "Velvet/Oases; idba_trans; SOAPDenovoTrans; Trinity;"

my %S2G_SUBDIRS = ( # informational for now,
  spotfa  => " 1. SRA spot (joined read pairs) files, from fastq-dump of SRAids",
  pairfa  => " 2. unjoined read pair files, _1.fa and _2.fa",
  rnasets => " 3. read pair rna sets, input to assemblers, various pairfa data slices",
  tra_XXX=>  " 4. subfolders per assembler/data slice",
  trsets  => " 5. assembled transcripts from several assembly runs", # also called subsets

  inputset => " 6. all transcripts/cds/aa from trsets as input to tr2aacds reduction",
  okayset  => " 7. non-redundant transcripts of tr2aacds, as gene locus primary (okay) and alternates (okalt)",
  dropset  => " 8. redundant transcripts of tr2aacds",

  refset   => " 9. reference sequences for annotation, eg refgenes.aa for homology, vector/contam screen",
  publicset => "10. public transcript/cds/aa sequences, annotations of evgmrna2tsa",
  submitset => "11. submission set for TSA database,  of evgmrna2tsa",

  genome  => "20. chromosome assembly, where available, for EvigeneH methods",
  aaeval  => "21. protein homology annotations, comparisons",
  geneval => "22. mRNA/CDS sequence annotations, comparisons",
);

my %DEFAULT_SETTINGS= ( 
  IDPREFIX=>$DEFAULTidpre, 
  DATE=>$DATE, 
  organism => $organism, 
  sraids => $sraids, 
  BioProject => $BioProject, #fixme bioproj lost
  assemblers => 'Velvet/Oases; idba_trans; SOAPDenovoTrans; Trinity;',
  trclass => '', mrna => '', genenames=>'', 
  ); 
  # vecscreen => '',
  # MINSIZE=>$MINSIZE, MAXGAP=>$MAXGAP,
  # LOCUSTAG => $DEFAULTidpre,
  # TSADESC=>$TSADESC, 

# our @evgdirs = qw(okayset dropset inputset trimset tmpfiles erasefiles publicset submitset);
my($runsteps,$runname,$sratable,$output,$logfile,$NCPU,$MAXMEM,$tidyup,$helpme) = (0) x 9;
my %settings= %DEFAULT_SETTINGS;
my ($sradatah, @sraids);

my @saveopt= @ARGV; shift @saveopt; # grep /^\-/, @ARGV;
my $optok= GetOptions(
  "runsteps=s", \$runsteps,   
  "runname=s", \$runname,   
  "organism|species=s", \$organism,   
  "SRAtable|datatable=s", \$sratable, "SRAids=s", \$sraids,   # one only
  "output:s",  \$output,
  "logfile:s", \$logfile,
  "idprefix=s", \$IDPREFIX,  # FIXME: idpre option  overwritten by spppref
  "DATE=s", \$DATE,  
  "NCPU=i", \$NCPU, "MAXMEM=i", \$MAXMEM,  
  "maxdatamb|MAX_DATASIZE_MB=i", \$MAX_SIZE_MB, 
  "dryrun|n!", \$dryrun, 
  "tidyup!", \$tidyup, # always on only?
  "debug!", \$DEBUG, 
  "help!", \$helpme, 
 );

  #   - need ref data options, path-to, methods : S8. refblast, S9. names (ref.aa names)
  # REFAA = refgenes.aa ; for now:  env REFAA=refset/refprots.aa VECDB=refset/UniVec  CONTAMDB=refset/contams .. 

  # "config=s", \$config, "cadd=s", \@configadd, #for evigene_config() .. want this?
  # "mrna|cdna=s", \$cdnaseq,
  # "class|trclass=s", \$trclass,
  # "names|genenames=s", \$genenames, ## ? allow for 2 files: myspecies.namerefids + allrefprot.names
  # "tblfile:s", \$tblfile,  
  # "MINSIZE=i", \$MINSIZE,  
  # "MAXGAP=i", \$MAXGAP,  
  # "runtbl2asn!", \$DOtbl2asn, 
  # "notrimvec|novectrim", \$SKIPTRIMSET, # * Change to -[no]vectrim .. only run if asked, as STEPnum

sub USAGE { 
  my $u= join"", "EvidentialGene sra2genes VERSION ",VERSION,"
  ombibus pipe for evigene methods, from SRA RNA-seq data to annotated public gene set
  
Usage: evgpipe_sra2genes.pl -SRAtable=myspecies_sra.csv | -SRAids=SRRnnnn,SRRmmmm 
opts: -help -runname MyProjectXXX -nCPU=$NCPU -idprefix Thecc1EG ..  
      -runstep 1,2,3,4..10  -log -dryrun -debug \n"; 

  if($_[0]) {
    my $ifo="\n";
    $ifo.="  *** EARLY DRAFT VERSION, Expect problems ***\n";
    $ifo.="---------------------------------------------------------------\n";
    
    $ifo.="Current pipleline design: 
    Process SRA RNA-seq data to a finished, annotated gene set, in steps, using existing, 
    tested Evigene methods.   Compute-intensive steps are run
    asynchronously, by generating cluster-ready shell scripts that you then submit to your
    cluster batch queue.  These steps include runassemblers, reduceassemblies, refblastgenes.
    
    See 'run_evgsra2genes.sh' an example cluster script to call this ombibus pipe. It sets
    paths to component software (assemblers, NCBI tools, others) that you must adjust.
    
    After these cluster runs, rerun this pipeline to proceed to next steps.  E.g.

evgpipe STEPs 1..4:
  env sratable=daphsim16huau_srarna.csv  name=dapsim_sra2evg  datad=`pwd` prog=./run_evgsra2genes.sh sbatch srun_shared.sh

ASYNC run assemblers (~ 8 hr each)
  env ncpu=8  datad=`pwd` prog=./runvelo.sh sbatch srun_comet.sh
  env ncpu=12  datad=`pwd` prog=./runidba.sh sbatch srun_comet.sh

evgpipe STEPs 5..7:
  env sratable=daphsim16huau_srarna.csv  name=dapsim_sra2evg  datad=`pwd` prog=./run_evgsra2genes.sh sbatch srun_shared.sh

ASYNC run assembly reduction to genes (~ 2 hr)
  env ncpu=20 maxmem=120000 prog=./run_tr2aacds.sh datad=`pwd` sbatch srun_comet.sh

evgpipe STEPs 8..9:
  env REFAA=refset/refarp7s10fset1.aa  sratable=daphsim16huau_srarna.csv  name=dapsim_sra2evg  datad=`pwd` prog=./run_evgsra2genes.sh sbatch srun_shared.sh

ASYNC run blastp (~ 8 hr)
  env ncpu=20 maxmem=120000 prog=./run_evgblastp.sh datad=`pwd` sbatch srun_comet.sh

evgpipe STEPs 10:
  env sratable=daphsim16huau_srarna.csv  name=dapsim_sra2evg  datad=`pwd` prog=./run_evgsra2genes.sh sbatch srun_shared.sh

STEPS in pipeline (will change)
  STEP1_sraget
  STEP2_sra2fasta
  STEP2a_sra2spot
  STEP2b_pairfa
  STEP3_selectrna
  STEP4_runassemblers
  STEP5_collectassemblies
  STEP6_qualassemblies
  STEP7_reduceassemblies
  STEP8_refblastgenes
  STEP9_namegenes
  STEP10_publicgenes
  More details via 'pod2man evgpipe_sra2genes.pl | nroff -man |less' \n";
    $ifo.="---------------------------------------------------------------\n";
    
    $ifo.="INPUT  -SRAtable=myspecies_sra.csv format is NCBI SraRunInfo.csv, 2017 format
    from https://www.ncbi.nlm.nih.gov/sra/ ( Send TO File, Format RunInfo)
Run,ReleaseDate,LoadDate,spots,bases,spots_with_mates,avgLength,size_MB,AssemblyName,download_path,
  Experiment,LibraryName,LibraryStrategy,LibrarySelection,LibrarySource,LibraryLayout,InsertSize,InsertDev,
  Platform,Model,SRAStudy,BioProject,Study_Pubmed_id,ProjectID,Sample,BioSample,SampleType,
  TaxID,ScientificName,SampleName,g1k_pop_code,source,g1k_analysis_group,Subject_ID,
  Sex,Disease,Tumor,Affection_Status,Analyte_Type,Histological_Type,Body_Site,CenterName,Submission,
  dbgap_study_accession,Consent,RunHash,ReadHash
Expected sra.csv format input may change; use of only -SRAids to be enabled.
Now requires NCBI sratoolkit/fastq-dump that has enabled web-fetch of data by SRAid.
That will become one option, others you fetch SRA/ENA data, or supply RNA-read-pairs.fasta/fastq
";   
    $ifo.="---------------------------------------------------------------\n";

    $ifo.="Layout of project directory:\n";
    map{ $ifo.="  $_:\t".$S2G_SUBDIRS{$_}."\n" } 
      sort{$S2G_SUBDIRS{$a} cmp $S2G_SUBDIRS{$b}} keys %S2G_SUBDIRS;
    $ifo.="\n";
    $ifo.="==============================================================\n";
    
    $u.= $ifo;
  }
  return $u;
}

die USAGE($helpme)  unless($optok and ($sratable or $sraids) );  

$runname ||= $organism."_SRA2Evigene"; # add $$? need constant name for run steps
#^^ reset after read sratable with organism?? but need for logfile; use org shortname,
# my $runtag= $settings{oname} || $runname;

$tidyup= 1 unless($dryrun); # ||$DEBUG default on unless debug|dryrun ?
# evigene_config($config, \@configadd); # always even if $config null
unless($DATE) { $DATE=`date '+%Y%m%d'`; chomp($DATE); } # fixme
$NCPU=1 unless($NCPU);
$MAXMEM= 50000 unless($MAXMEM);

openloggit($logfile,$runname); # which first? publicset/logfile or ./logfile ?
loggit(1, "EvidentialGene evgpipe_sra2genes.pl (-help for info), VERSION",VERSION);
loggit(1, "CMD: evgpipe_sra2genes.pl ",@saveopt);
loggit(1, "ERR: unused arguments:",@ARGV) if(@ARGV>0);

MAIN_sub(); 

=item Async STEP notes

  STEP4_runassemblers step should write script to run async, under various cluster batch control systems
      .. then restart at next step.  Ditto some other steps: 
  STEP7_reduceassemblies  runtr2cds.sh, 
  STEP8_refblastgenes     run_evgaablast
  
  add async alt steps, like
  STEP3b_selectrna:diginorm
  
=cut

use constant { STEPnotdone => 0, STEPok => 1, STEPerr => -1, STEPdone => 2,  };
 
sub MAIN_sub {  
  loggit(0, "BEGIN with input=",$sratable||$sraids,"date=",`date`);

	do_settings("restore",$runname,); # ("log|restore|save");

  # $donestep= .. # check for step output filesets
  # if($runsteps =~ /.1,/) ...
  my @stepok=(0);
  
  # ... prepare data
  @stepok= STEP1_sraget($sratable,$sraids);     # combine 1/2, sra fastq-dump does web-get now  
  loggit(LOG_DEBUG,"done STEP1_sraget", @stepok);
  
  @stepok= STEP2_sra2fasta();  # prep data
  loggit(LOG_DEBUG,"done STEP2_sra2fasta", @stepok);
  
  ##>> need option alternates here, several data slices to be assembled, with reruns for new slices
  ##   also ensure user-supplied selrna data slices can be used: subfolder rnasets/ ?
  ## add to step2/3: SCRIPT_diginorm(), works on spotfa/ fasta, takes time to run
  
  @stepok= STEP3_selectrna();  # prep data
  loggit(LOG_DEBUG,"done STEP3_selectrna", @stepok);
  # -------
	do_settings("save",$runname,); #  ("log|restore|save");
  
  @stepok= STEP4_runassemblers(); # writes scripts, async step to run..
  loggit(LOG_DEBUG,"done STEP4_runassemblers", @stepok);
  
  @stepok= STEP5_collectassemblies();  # postproc asm
  loggit(LOG_DEBUG,"done STEP5_collectassemblies", @stepok);
  
  @stepok= ($runsteps=~/6/)? STEP6_qualassemblies() : (STEPnotdone);     # postproc asm .. OPTIONAL
  loggit(LOG_DEBUG,"done STEP6_qualassemblies", @stepok);
  
  @stepok= STEP7_reduceassemblies();   # tr2aacds
  loggit(LOG_DEBUG,"done STEP7_reduceassemblies", @stepok);
  
  # ------------------
  # annotate .. various steps, substeps
  
  #x @stepok= ($runsteps=~/8/) ?  STEP8_refblastgenes() : (STEPnotdone);      # annotate .. various
  @stepok= STEP8_refblastgenes();  
  loggit(LOG_DEBUG,"done STEP8_refblastgenes", @stepok);

  @stepok= STEP9_namegenes();
  loggit(LOG_DEBUG,"done STEP9_namegenes", @stepok);
  # ADD: vecscreen(); contamcheck(); consdomains();  transposons();
  
  @stepok= STEP10_publicgenes();       # finish
  loggit(LOG_DEBUG,"done STEP10_publicgenes", @stepok);

	do_settings("save",$runname,); # was "log|save";   ("log|restore|save");
  loggit(LOG_DEBUG,"settings saved to $runname.info");

  loggit(0, "DONE at date=",`date`);
  loggit(0, "======================================"); # log may append runs
} # end MAIN


sub STEP1_sraget {
  my($sratable,$sraidlist)= @_;
  loggit(LOG_DEBUG,"STEP1_sraget");
  ## NOTE: current sratoolkit/fastq-dump will fetch from NCBI-sra by 'SRRnnn' accession.
  ##  .. skip this wget, use fastq-dump ..
  #   1. get data
  #     sraget.pl : small script, enclose
  #     = wget to ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/INSERTDATAPATH.sra
  #     env fork=1 ./sraget.pl daphsim16huau_srarna.csv  >& get1.log &

  my $nsra=0;
  ($sradatah,$nsra)= get_srainfo($sratable,$sraidlist); 
	# $gotsraids == global $sraids == new $sraidlist
  # sraget($sradatah); # skip this one
  
  return ($nsra>0) ? STEPok : STEPerr;
}


sub STEP2_sra2fasta {
  loggit(LOG_DEBUG,"STEP2_sra2fasta");

  #   2. sra to fasta, uses NCBI sratoolkit/fastq-dump
  #     env sra=SRX*/SRR*/SRR*.sra datad=`pwd` prog=./runsra2fa.sh sbatch ../srun_shared.sh
  #   uses fqdump=$HOME/bio/sratoolkit/fastq-dump

  my $ok=0;
  $ok= STEP2a_sra2spot();
  $ok= STEP2b_pairfa() if($ok>0);
  
  if($ok) {
    my($oks,$runapp,$runfile)= SCRIPT_diginorm(); # add here? can use in STEP3 selectrna
    my($aok,$aname)= add_script($oks,$runapp,$runfile,'diginorm','STEP2_sra2fasta'); 
  }
  
  return $ok;
}


sub STEP2a_sra2spot {

  ## DONT need to web-fetch w/ fastq-dump, use SRRids
  # FIXME: @srafiles == URL now, from %sradata{download_path}
  # my @srafiles= grep /http|ftp/, split/;/, 
  #   $sradatah->{"download_path"} || $sradatah->{"FTP Path to Experiment"};  

  # @sraids now global;
  @sraids= grep /^[A-Z]+\d+/, split /\W+/, $sraids; # NOW GLOBAL; should be SRRnnnn ERRnnnn list

  loggit(LOG_DIE,"missing SRA IDs: SRRnnn") unless(@sraids);
  loggit(0,"sra2fasta ids:",@sraids);
  
  # cd $datad
  my $havespots= (-d "spotfa" and -s "spotfa/$sraids[0].fasta")?1:0;
  if($havespots) {
    loggit(0,"sra2fasta using","spotfa/$sraids[0].fasta");
    return STEPdone;
  } else {
    my $fqdump= findapp('fastq-dump'); # does loggit(LOG_DIE,"missing ...") if($findapp =~ /MISSING/);
    mkdir("spotfa"); 
    my $icpu=0; 
    foreach my $sr (@sraids)  { 
      ## should trap STDERR from fqdump.. in openloggit() ?? in cmd line "2>&1" ?
      my $cmd="$fqdump -O spotfa --qual-filter --fasta 0 $sr";
      my $pid= forkcmd($cmd);  
      if(++$icpu >= $NCPU) { while (wait() != -1) { }; $icpu= 0; }
    } 
    while (wait() != -1) { };
    $havespots= (-d "spotfa" and -s "spotfa/$sraids[0].fasta")?1:0;
  }
  return ($havespots) ? STEPok : STEPerr;
}


sub STEP2b_pairfa {
  ## FIXME: dont replace existing pairfa/*.fa **
  # now global# my @sraids= grep /[A-Z]\d+/, split /\W+/, $sraids; # should be SRRnnnn ERRnnnn list
  my $havepairs= (-d "pairfa" and -s "pairfa/$sraids[0]_2.fa")?1:0;
  my %PAIRINFO=();
  
  if($havepairs) {
    return STEPdone if($settings{readsize} and $settings{nreads});
  }
  
  unless($havepairs) {
    my ($spotall,@spofa)= getFileset('spotfa','fasta$');  
    return STEPerr unless( @spofa); 
    mkdir("pairfa") unless(-d "pairfa");
    my $icpu=0; 
    foreach my $sfa (@spofa)  { 
      chomp($sfa); 
      my $name=$sfa; 
      $name =~ s/\.fasta.*//; $name =~ s/spotfa/pairfa/;  
      #dont need# $PAIRINFO{$name}{nreads}= $PAIRINFO{$name}{maxlen}= 0; # set global
      my $fa1="${name}_1.fa";
      if(-s $fa1) { warn "exists: $fa1\n"; next; }
      
      # use this fork
      my $pid= fork();
      if($pid) { # parent
        # doesn't get global var sets back from child
        
      } else { # child
        my($cerr,$nreads,$maxlen,$totlen,$lfn,$rfn)= splitspots($name,$sfa);
        if($cerr == 0) {
        open(FI,'>',"$name.fa.info"); 
        print FI "nreads=$nreads;maxlen=$maxlen;totlen=$totlen;lfn=$lfn;rfn=$rfn\n"; 
        close(FI);
        }
        ## child cant set parent global vars, use info file
        ## if($cerr != 0) { $PAIRINFO{$name}{error}=$cerr; ... }
        exit($cerr)
      }
      if(++$icpu >= $NCPU) { while (wait() != -1) { }; $icpu= 0; }
      } 
    while (wait() != -1) { };
    $havepairs= (-d "pairfa" and -s "pairfa/$sraids[0]_2.fa")?1:0;
  }
  
  if($havepairs) {
    my($nreads,$readsize)=(0) x 9;
    %PAIRINFO= get_readpair_info('pairfa'); # UPD: keys are SID not pairfa/SID
    for my $sid (sort keys %PAIRINFO) {
      my $pn= "pairfa/$sid"; # UPD
      my $nr= $PAIRINFO{$sid}{nreads}||0;
      my $ml= $PAIRINFO{$sid}{maxlen}||0; # totlen ?  
      $readsize= $ml if($ml>$readsize);
      $nreads += $nr;
      my $vals= join";", map{ "$_=".$PAIRINFO{$sid}{$_} } sort keys %{$PAIRINFO{$sid}}; 
      $sradatah->{$pn}= $vals; # pn == 'pairfa/SRRnnnn' okay?
    }
    
    # return stats from splitspots(): max readsize, nreads/npairs ..
    $settings{readsize} = $sradatah->{readsize} = $readsize;
    $settings{nreads}   = $sradatah->{nreads} = $nreads;
  }
  return ($havepairs) ? STEPok : STEPerr;
}



sub STEP3_selectrna_OLD {
  loggit(LOG_DEBUG,"STEP3_selectrna");

=item STEP3_selectrna info

  3. subset for data size (max data size for basic assembly : 15-20GB for given maxmem=120 GB
    3a. option: clean/trim rna.fa
    3b. option: digital normalize rna.fa
    3c. option: make selrna/trasm set for each SRR.fa >= MAX_SIZE_MB 
          .. get more genes fully asm from data that way, along w/ diginorm,

  FIXME: for several selectrna data slices
        need data/asm file tags for each slice: sb; sp; sn ?
        sb=basic 1 from all/nonorm, sr=1 per run/pairfa/SRR, sn=diginorm of all
        
  Data slice naming,location:
     try1: $outna="srr$samplen". ($sradatah->{BioProject} || $IDPREFIX ); 
     try2: $outna="(sb|sp|sn)[num]$PRID" ; and keep "(sb|sp|sn)[num]" with trsets/asms.tr
        
  MAX_SIZE_MB =  10Gb default, changeable, is max sized for all assemblers,
      depends on system rez, MAXMEM, NCPU; may want to calc some data size options
      
=cut

  my $SUBT='basic'; # | run=pair | norm
  
  ## fixed: BUG mixed order of _1/2
  ## also STEP2 sets : NOT YET, fork child vars.
  #  $settings{readsize} = $sradatah->{readsize} = $readsize;
  #  $settings{nreads}   = $sradatah->{nreads} = $nreads;

  # srainfo: spots,bases,spots_with_mates,avgLength,size_MB .. use spots,aveLength == bases
  my($tbases,$tspots,$lenspot)=(0) x 9;
  #now global# my @sraids= grep /[A-Z]\d+/, split /\W+/, $sraids; # should be SRRnnnn ERRnnnn list
  # FIXME: change cvs1 hdr to cvs2 so dont need this 
  # FIXME2: match $sradatah->{sraids} to each id attr
  if ($sradatah->{cvsformat} == 1) { 
    map{ $tbases += $_ } split";",$sradatah->{"Total Bases"};
    map{ $tspots += $_ } split";",$sradatah->{"Total Spots"};
  } else {  
    map{ $tbases += $_ } split";",$sradatah->{"bases"};
    map{ $tspots += $_ } split";",$sradatah->{"spots"};
  }
  
  my $tmb= int($tbases/1000000);
  my $maxread= $tspots;
  my $samplen= 1;
  if( $tmb > $MAX_SIZE_MB) {
    $maxread= int($tspots * $MAX_SIZE_MB/$tmb);
    $samplen= 1 + int($tmb/$MAX_SIZE_MB); # 2,3,.. read steps
  }
  ##bad# my $maxread= int($tspots/$samplen); ## wrong.. too few; need fraction
  
  my ($np1,$np2)=(0) x 9;
  ## should  outna have subdir?
  
  mkdir('rnasets');
  #bad# my $outna="srr$samplen". ($sradatah->{BioProject} || $IDPREFIX ); # BAD BioProj many;vals;...
  my $sdid= $sraids[0];
  my $st= ($SUBT =~ /norm/)?'sN':($SUBT =~ /basic/)?'sB':($SUBT =~ /run|pair/)?'sR':'sO';
  my $sn= 'n'.scalar(@sraids);
  my $sl= 'l'.$samplen;
  my $outna= "rnasets/$st$sn$sl$sdid";
  # $settings{readset} .= "$outna;" unless($settings{readset} =~ /$outna/);
  
  my $readset1=$outna."_1.fa";
  my $readset2=$outna."_2.fa";
  
  my ($pairall,@pairs)= getFileset('pairfa','fa$'); # cant do _1. in getFileset()

  # FIXME: grep @sraids, maybe more in pairfa than want..
  my $sraidset= join '|',@sraids;
  my @mypairs= grep /$sraidset/, @pairs;
  
  ## BUG here, sort @pair1,2
  my @pair1= sort grep /_1.fa/, @mypairs;
  my @pair2= sort grep /_2.fa/, @mypairs;
  my $pok= (@pair1 and @pair2 == @pair1);

  my $havereadset= ( -s $readset1 and -s $readset2);
  if($havereadset) { loggit(0,"readsample, reusing $readset2"); $pok=0; }
  
  my $pselreads= $maxread/$tspots; $pselreads=1 if($pselreads>1);
  my $pmax= int(100*$pselreads);  
  loggit(0,"readsample nreads=$maxread ($pmax% for $MAX_SIZE_MB maxMB), to $outna, pok=$pok");
  # report nread/pair?
  if($pok) {

# FIXME: slow, use instead unix: "head -n$NLINES >> $outfile"  ; 
#  .. need nreads & %keep = %pmax per pairfa
use constant USE_UHEAD => 1;
    my %PAIRINFO= get_readpair_info('pairfa'); # nreads per pairfa
    ##  my $nr= $PAIRINFO{$pn}{nreads}||0;

    open(O,'>',$readset1) or loggit(LOG_DIE,"write $readset1");
if(USE_UHEAD) {
    close(O);
}    
    for my $pf (@pair1) {
    
if(USE_UHEAD) {
      (my $pn=$pf)=~s/_[12].fa.*//;
      my ($sid)= $pf=~m/(\w+)_[12].fa/; # *should be ok, maybe not..
      my $nr= $PAIRINFO{$sid}{nreads}||0; #? fail if missing
      $np1= int($pselreads * $nr); # should count file
      my $nrkeep= 2 * $np1;
      # BUGGERS, 2 lines per read: 2*nrkeep for head
      my $cmd= "head -n $nrkeep $pf >> $readset1";
      my $err= runcmd($cmd);
      
} else {      
    
      open(P,$pf) or loggit(LOG_DIE,"read $pf");
      while(<P>) { 
        if(/^>/){ last if($np1 >= $maxread); $np1++; }
        #try# if($samplen>1) { next if($np1 % $samplen == 1); } 
        print O $_; 
      } close(P);
}
      
    } 
unless(USE_UHEAD) {
    close(O);
}
    open(O,'>',$readset2) or loggit(LOG_DIE,"write $$readset2");
if(USE_UHEAD) {
    close(O);
}    
    for my $pf (@pair2) {
if(USE_UHEAD) {
      (my $pn=$pf)=~s/_[12].fa.*//;
      my ($sid)= $pf=~m/(\w+)_[12].fa/; # *should be ok, maybe not..
      my $nr= $PAIRINFO{$sid}{nreads}||0;
      $np2= int($pselreads * $nr); # should count file
      my $nrkeep= 2 * $np2;
      my $cmd= "head -n $nrkeep $pf >> $readset2";
      my $err= runcmd($cmd);
} else {      
      open(P,$pf) or loggit(LOG_DIE,"read $pf");
      while(<P>){ 
        if(/^>/){ last if($np2 >= $maxread); $np2++; }
        #try# if($samplen>1) { next if($np2 % $samplen == 1); } 
        print O $_; 
      } close(P);
    } 
  }    
unless(USE_UHEAD) {
    close(O);
}
    
    # open all pair1/2, write to topdata1,2, subsampling
    loggit(0,"readsampled to $outna pair1,2=$np1,$np2");
    $havereadset= ( -s $readset1 and -s $readset2);
  }
  
  if($havereadset) { 
    $settings{readset}  = $sradatah->{readset} = $outna;
    ## allow many readsets ..  seta;setb;setc;.. different set key: rnasets?
    $settings{rnasets} .= "$outna;" unless($settings{rnasets} =~ /$outna/);
    #>> check np1, missing if done once..
    $settings{nreadset} = $sradatah->{nreadset} = $np1 if($np1>0); # ~ nreads == maxread
  }
  return ($havereadset) ? STEPok : STEPerr;
}


sub STEP3_selectrna {
  loggit(LOG_DEBUG,"STEP3_selectrna new");

  my %PAIRINFO= get_readpair_info('pairfa'); # nreads per pairfa
  my ($havereadset,$nreadsets)=(0,0);
  my (@perrun, @mixrun, @mixset, $sumtmb);
  
  #old# my @pid= grep /pairfa/, sort keys %PAIRINFO; # dang, has both pathto/SID and SID keys, drop pathto/
  my @sid= sort keys %PAIRINFO; # dang, has both pathto/SID and SID keys, drop pathto/
  my $sraidset= join '|', @sraids;
  @sid= grep /$sraidset/, @sid;
  
  ## sort pid by totlen ??
  for my $sd (@sid) {
    #old# my($sd)= (split"/",$pd)[-1];
    my $nspots= $PAIRINFO{$sd}{nreads}; 
    my $nbases= $PAIRINFO{$sd}{totlen}; 

    my $tmb= int($nbases/1000000);
    my $maxread= $nspots;
    if( $tmb > 0.90 * $MAX_SIZE_MB) {
      push @perrun, $sd;
    } else {
      push @mixrun, $sd;
      $sumtmb += $tmb;
      if( $sumtmb > 0.90 * $MAX_SIZE_MB) {
        push @mixset, join(";",@mixrun); 
        # push @mixset, [@mixrun]; 
        @mixrun=(); $sumtmb=0;
      }
    }    
  }
  
  for my $sd (@perrun) {
    my($sok)= sample_reads('pairs',$sd,\%PAIRINFO);
    $nreadsets++ if($sok>0);
  }
  for my $set (@mixset) {
    # my $sds= join";", @$set;
    my($sok)= sample_reads('mixset',$set,\%PAIRINFO);
    $nreadsets++ if($sok>0);
  }
  return ($nreadsets) ? STEPok : STEPerr;
}


sub sample_reads {
  my($subtype,$sraids,$pairinfo)=@_;
  ## loggit(LOG_DEBUG,"STEP3_selectrna");

  my %PAIRINFO;
  if(ref($pairinfo)) { %PAIRINFO= %$pairinfo; }
  else {  %PAIRINFO= get_readpair_info('pairfa'); } # nreads per pairfa
  
  my $SUBT=$subtype; # 'basic'; # | run=pair | norm
  
  my @sraids= grep/\w/, split/\W+/, $sraids;  
  
  my ($tspots,$tbases)=(0,0); 
  for my $sd (@sraids) { $tspots+=$PAIRINFO{$sd}{nreads}; $tbases+=$PAIRINFO{$sd}{totlen}; }  
  my $tmb= int($tbases/1000000);
  my $maxread= $tspots;
  my $samplen= 1;
  if( $tmb > $MAX_SIZE_MB) {
    $maxread= int($tspots * $MAX_SIZE_MB/$tmb);
    $samplen= 1 + int($tmb/$MAX_SIZE_MB); # 2,3,.. read steps
  }
  
  my ($np1,$np2)=(0) x 9;
  mkdir('rnasets');
  my $sdid= $sraids[0];
  my $st= ($SUBT =~ /norm/)?'sN':($SUBT =~ /basic|mix/)?'sB':($SUBT =~ /run|pair/)?'sR':'sO';
  my $sn= 'n'.scalar(@sraids);
  my $sl= 'l'.$samplen;
  my $outna= "rnasets/$st$sn$sl$sdid";
  my $readset1=$outna."_1.fa";
  my $readset2=$outna."_2.fa";
  
  my ($pairall,@pairs)= getFileset('pairfa','fa$'); # cant do _1. in getFileset()
  # FIXME: grep @sraids, maybe more in pairfa than want..
  my $sraidset= join '|',@sraids;
  my @mypairs= grep /$sraidset/, @pairs;
  
  ## BUG here, sort @pair1,2
  my @pair1= sort grep /_1.fa/, @mypairs;
  my @pair2= sort grep /_2.fa/, @mypairs;
  my $pok= (@pair1 and @pair2 == @pair1);

  my $havereadset= ( -s $readset1 and -s $readset2);
  if($havereadset) { loggit(0,"readsample, reusing $readset2"); $pok=0; }
  
  my $pselreads= $maxread/$tspots; $pselreads=1 if($pselreads>1);
  my $pmax= int(100*$pselreads);  
  loggit(0,"sample_reads type=$subtype, nreads=$maxread/$tspots ($pmax% for $MAX_SIZE_MB maxMB) of ids=$sraids to $outna, pok=$pok");
  # report nread/pair?
  if($pok) {
    # FIXME: slow, use instead unix: "head -n$NLINES >> $outfile"  ; 
    # use constant USE_UHEAD => 1;
    #have# my %PAIRINFO= get_readpair_info('pairfa'); # nreads per pairfa
    ##  my $nr= $PAIRINFO{$pn}{nreads}||0;

    open(O,'>',$readset1) or loggit(LOG_DIE,"write $readset1");
    close(O); # USE_UHEAD   
    for my $pf (@pair1) {
      # if(USE_UHEAD) 
      (my $pn=$pf)=~s/_[12].fa.*//;
      my ($sid)= $pf=~m/(\w+)_[12].fa/;
      my $nr= $PAIRINFO{$sid}{nreads}||0; #? fail if missing
      $np1= int($pselreads * $nr); # should count file
      my $nrkeep= 2 * $np1;
      # BUGGERS, 2 lines per read: 2*nrkeep for head
      my $cmd= "head -n $nrkeep $pf >> $readset1";
      my $err= runcmd($cmd);
    } 

    open(O,'>',$readset2) or loggit(LOG_DIE,"write $$readset2");
    close(O);
    for my $pf (@pair2) {
      # if(USE_UHEAD) 
      (my $pn=$pf)=~s/_[12].fa.*//;
      my ($sid)= $pf=~m/(\w+)_[12].fa/; # *should be ok, maybe not..
      my $nr= $PAIRINFO{$sid}{nreads}||0;
      $np2= int($pselreads * $nr); # should count file
      my $nrkeep= 2 * $np2;
      my $cmd= "head -n $nrkeep $pf >> $readset2";
      my $err= runcmd($cmd);
    }    
    
    # open all pair1/2, write to topdata1,2, subsampling
    loggit(0,"readsampled to $outna pair1,2=$np1,$np2");
    $havereadset= ( -s $readset1 and -s $readset2);
  }
  
  if($havereadset) { 
    $settings{readset}  = $sradatah->{readset} = $outna;
    ## allow many readsets ..  seta;setb;setc;.. different set key: rnasets?
    $settings{rnasets} .= "$outna;" unless($settings{rnasets} =~ /$outna/);
    #>> check np1, missing if done once..
    $settings{nreadset} = $sradatah->{nreadset} = $np1 if($np1>0); # ~ nreads == maxread
  }
  return ($havereadset) ? STEPok : STEPerr;
}


sub STEP4_runassemblers {
  loggit(LOG_DEBUG,"STEP4_runassemblers");
  my($whichasm,$aname,$nok)=("",0,0);
  
#   4.  run assemblers, with kmer size options, other opts
#     4a. velvet/oases, ~10 kmer steps
#     4b. idba_tran, ~10 kmer steps
#     4c. soap_trans, ~10 kmer steps
#     4d. trinity / other / user choices
#     -- est. time at 12hr x 8cpu x 120 GB mem per assembler multi-kmer set (same for 1kmer trinity)

  if(my $needsrun= $settings{STEP4_runassemblers}) {
    # clear?  return w/ notice?  
    # OR check for outputs of these.. move on to STEP5 if have
    loggit(LOG_WARN,"STEP4_runassemblers have scripts: $needsrun");
    return(STEPdone); # STEPerr?
  }

  if($settings{assemblers} =~ /velv/i) {
    my($ok,$sh,$outsh)= SCRIPT_velvetoases();
    ($ok,$aname)= add_script($ok,$sh,$outsh,'velveto','STEP4_runassemblers'); # velvetoases
    if($ok>0) { $nok++; $whichasm.= "$aname,"; }
    # if($ok>0) {
    #   open(S,'>',$outsh); print S $sh; close(S);
    #   system("chmod +x $outsh"); #? or not
    #   $settings{STEP4_runassemblers} .="$outsh;";
    # }
  }

  if($settings{assemblers} =~ /idba/i) {
    my($ok,$sh,$outsh)= SCRIPT_idba_trans();
    ($ok,$aname)= add_script($ok,$sh,$outsh,'idba','STEP4_runassemblers');
    if($ok>0) { $nok++; $whichasm.= "$aname,"; }
  }

  if($settings{assemblers} =~ /soap/i) {
    my($ok,$sh,$outsh)= SCRIPT_soap_trans();
    ($ok,$aname)= add_script($ok,$sh,$outsh,'soap','STEP4_runassemblers');
    if($ok>0) { $nok++; $whichasm.= "$aname,"; }
  }
  
  if($settings{assemblers} =~ /trinity/i) {
    my($ok,$sh,$outsh)= SCRIPT_trinity();
    ($ok,$aname)= add_script($ok,$sh,$outsh,'trinity','STEP4_runassemblers');
    if($ok>0) { $nok++; $whichasm.= "$aname,"; }
  }
  
   # return ($nok) ? STEPok : STEPerr;
  return ($whichasm) ? (STEPok,$whichasm) : (STEPerr,$whichasm);
}


sub STEP5_collectassemblies {
  loggit(LOG_DEBUG,"STEP5_collectassemblies");

#   5. post process assembly sets (trformat.pl)
#     $evigene/scripts/rnaseq/trformat.pl -pre $nap -out $subd.tr -log -in $subd/vel*/transcripts.fa 
#     .. etc per assembler/output set
#     -- collect into evgrundir/subsets/

# see daphplx/rnasm/evg1dapsim/evigene_methods/runscripts_dpx17/trformgroup.sh

  my $needsrun= $settings{STEP4_runassemblers};
  return (STEPerr) unless($needsrun);
  
  mkdir("trsets");
  my $runtag= $settings{oname} || $runname;
  loggit(LOG_DEBUG,"runtag=$runtag","runassemblers=",$needsrun);
  #s2g.STEP4_runassemblers=runvelo.sh;runidba.sh;runsoap.sh;runtrin.sh

  # HERE check for outputs of these.. 
  my @runs= split";",$needsrun; 
  my %trsets=();
  for my $rs (@runs) {
    #? use settings to hold runscript outdir/outfiles?
    ## now have 'trvelv', 'trsoap', 'tridba', 'trtrin' ...
    ## FIXME: trvelv > trvelo needs to match runvelo.sh
    
    my($rundir)= "nada";
    (my $rd=$rs)=~s/\.\w+$//; $rd=~s/run//;  # guess at std outdir name
    ($rundir)=`/bin/ls -d tr$rd*`; chomp($rundir);

## trap errs..
# ls: trvelo*: No such file or directory
# ls: tridba*: No such file or directory

    my @rundir= split" ",$rundir; # many maybe, but should be 1 per run .. 
    # need better file names: run_name.sh == tra_name outdir ?
    #bad,params# if(open(F,$rs)) { while(<F>){ if(m/^outdir=(\S+)/) { $rundir=$1; last; } } close(F); }
    
    if( $rundir and -d $rundir ) {
      # per method: 
      # velo  = rundir/vel_kset../transcripts.fa 
      # idbat = rundir/transcript-kmer.fa and rundir/contig.fa (final merge)
      #  ?? this also: $xbin/fastanrdb idba-trs.tr > idba-trs.nrtr ; many dups
      # soap = rundir/soap_kset../xxx.scafSeq  
      # trin = rundir/Trinity.fasta ??
      my @trs; my $dummy;
      if($rundir =~ /vel/) { ### Dang name bugs: velv velo velveto ...
        @trs=`ls $rundir/*/transcripts.fa`;  # gz also?
      } elsif($rundir =~ /idba/) {
        @trs=`ls $rundir/transcript-*.fa`;  #  gz also? also contig.fa ?
      } elsif($rundir =~ /soap/) {
        @trs=`ls $rundir/sod*/so*.scafSeq.gz`; #? gz or not
      } elsif($rundir =~ /trin/) {
        @trs=`ls $rundir/Trinity*.fasta`; #  gz also?
      } else {
        ($dummy,@trs)= getFileset( $rundir,'fa$|fasta$'); # uck, skip? warn?
      }

        # check if trsets/$rundir.tr exists ..
        # may have cases of partial finish, count @trs each call, check trformat.log count
      my $trout= "trsets/$rundir.tr";
      if( -s $trout) { 
        $trsets{$rundir}= STEPdone;
      } elsif(@trs) { 
        chomp(@trs);   
        my @cmd= ("$EVIGENES/rnaseq/trformat.pl", "-pre", $runtag,
           "-out", $trout, "-log", "-in", @trs);
        my $err= runcmd(@cmd); 
        $trsets{$rundir}= ($err) ? STEPerr : STEPok;
      } else {
        $trsets{$rundir}= STEPnotdone;
      }
    } else {
      $trsets{$rundir}= STEPnotdone; # missing?
    }
  }

  my($dummy,@trs)= getFileset( "trsets",'tr$'); 
  # return (@trs) ? STEPok : STEPerr;
  
  my ($trok,$trerr,$trwait)=(0,0,0);
  my @trsets= sort keys %trsets;
  for my $ts (@trsets) { my $v= $trsets{$ts}; if($v>0){ $trok++; } elsif($v<0) { $trerr--; } else { $trwait++; } }
  if(@trs > $trok) { $trok= @trs; } #?
  if($trok > 0 and $trerr == 0 and @trs) { return STEPok; }
  elsif($trerr > 0) { return STEPerr; }
  elsif($trwait > 0) { return STEPnotdone; }
  elsif(@trs > 0) { return STEPok; } # maybe
  else { return STEPnotdone; }
}

sub STEP6_qualassemblies {
  loggit(LOG_DEBUG,"STEP6_qualassemblies");

## OPTION: do later ? do after tr2aacds ?
#   6. quick qual assessment: cdna_bestorf > aastats per assembly, report
#     $evigene/scripts/cdna_bestorf.pl -nostop -minaa=30 -aa -cdna $subd.tr.gz 
#     $evigene/scripts/prot/aaqual.sh $subd.aa
#     $evigene/scripts/prot/aastat.sh $subd.aa.qual >> assembly.aastat

  my($trall,@trs)= getFileset( "trsets",'tr$'); 
  return (STEPerr) unless(@trs);
  
  my $aastat="trsets/assembly_aastat.txt";
  my($tralla,@aa)= getFileset( "trsets",'aa$'); 
  return STEPdone if(@aa and @aa == @trs and -f $aastat);

  my $MINCDS = $ENV{MINCDS} || 90; # tr2aacds set
  my $MINAA= int($MINCDS/3);

  my $icpu=0; my $haveaa=0;
  for my $cdna (@trs) {
    ## if we do this step, do for input to tr2aacds .. tr2aacds is more efficient with many NCPU
    #  ($cdsseq,$aaseq)= make_bestorf_ncpu($NCPU,$cdnaseq,$cdsseq,$aaseq); 
    #  my $cmd1="$APPcdnabest -nostop -minaa=$MINAA -cdna $cdna1 -aaseq $aa1 -cdsseq $cds1";
    
    (my $aaf=$cdna)=~s/\.(tr|fa|cdna)\w*$//; $aaf.=".aa";
    if( -s $aaf) { $haveaa++; next; }
    my @cmd= ("$EVIGENES/cdna_bestorf.pl", "-nostop","-minaa=$MINAA", "-aaseq", "-cdsseq", "-cdna", $cdna);
    my $pid= forkcmd(@cmd);  
    if(++$icpu >= $NCPU) { while (wait() != -1) { }; $icpu= 0; }
  } 
  while (wait() != -1) { };

  ($tralla,@aa)= getFileset( "trsets",'aa$'); 
  for my $aa (@aa) {
    # env stat=1 span=1 $evigene/scripts/prot/aaqual.sh $aa == subs:makeAaQual()
    my($aaqual)= makeAaQual($aa);
    if(-f $aaqual) {
      my $cmd="$EVIGENES/prot/aastat.sh $aaqual >> $aastat";
      runcmd($cmd);
    }
  } 

  return (@trs) ? STEPok : STEPerr;
}


sub STEP7_reduceassemblies {
  loggit(LOG_DEBUG,"STEP7_reduceassemblies");

#   7. run evg over-assembly reduction, tr2aacds.pl
#     env trset=$pt.tr datad=`pwd` prog=./runtr2cds.sh sbatch srun_comet.sh

  if(my $needsrun= $settings{STEP7_reduceassemblies}) {
    loggit(LOG_WARN,"STEP7_reduceassemblies have scripts: $needsrun");
    return(STEPdone); # STEPerr?
  }

  ## Note SCRIPT sub makes input.tr data file from trsets/*.tr
  my($ok,$runapp,$runfile)= SCRIPT_tr2aacds();
  my($aok,$aname)= add_script($ok,$runapp,$runfile,'tr2aacds','STEP7_reduceassemblies'); 

  return ($aok>0) ? STEPok : STEPerr;
}

sub STEP8_refblastgenes {
  loggit(LOG_DEBUG,"STEP8_refblastgenes");

#   8. ref protein blastp x evg okayset
#     env aaset=okayset/$pt.aa refaa=refset/$refaa ncpu=20 datad=`pwd` prog=./run_evgaablast.sh sbatch srun_comet.sh 
  
  # FIXME: need ref data options, path-to, methods : 
  # S8. refblast, S9. names (ref.aa names)
  # S7..S9: vecscreen UniVec, contamscreen NCBI contamseqs,rRNA, 
  # S8b: NCBI CDD or other cons domain data, rpsblast or hmmer
   
  if(my $needsrun= $settings{STEP8_refblastgenes}) {
    loggit(LOG_WARN,"STEP8_refblastgenes have scripts: $needsrun");
    return(STEPdone); # STEPerr ?
  }

  my($ok,$runapp,$runfile)= SCRIPT_evgblastp();
  my($aok,$aname)= add_script($ok,$runapp,$runfile,'evgblastp','STEP8_refblastgenes'); 

  return ($aok>0) ? STEPok : STEPerr;
}

sub STEP9_namegenes {
  loggit(LOG_DEBUG,"STEP9_namegenes");

#   9. namegenes from ref names, for annotation
#     $evigene/scripts/prot/namegenes.pl -blast $aabltab -refnames refset/$refaa.names -out $pt.names
# * ADD this to STEP8_refblastgenes , SCRIPT_evgblastp
# .. but need ref.names collection (if refgenes.aa but not refgenes.names)

  return (0) ? STEPok : STEPerr;
}

sub STEP10_publicgenes {
  loggit(LOG_DEBUG,"STEP10_publicgenes");

#   10. annotated publicset 
#     env idprefix=$idp trclass=$pt.trclass names=$pt.names  species=$spp datad=`pwd` \
#        prog=./run_evgmrna2tsa.sh sbatch srun_comet.sh

  #   my $P_TRCLASS= $settings{trclass} || "$trname.trclass";

  if(my $needsrun= $settings{STEP10_publicgenes}) {
    loggit(LOG_WARN,"STEP10_publicgenes have scripts: $needsrun");
    return(STEPdone); # STEPerr?
  }

  my($ok,$runapp,$runfile)= SCRIPT_evgpubset();
  my($aok,$aname)= add_script($ok,$runapp,$runfile,'evgpubset','STEP10_publicgenes'); 

  return ($aok>0) ? STEPok : STEPerr;
}


#---------------------------------------------------  

sub do_settings {  ## cdna_evigenesub.pm: sub evigene_config() ??
  my($action,$pathname)= @_;
  ## write these to work dir; reread next go
  ## action == 'log|save|restore' ; restore called AFTER read new options, shouldnt replace

  my $PRES=$EGLOG; # == 's2g';
  my $trpname= makename($pathname,".sra2genes.info"); 
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
		organism => $organism, BioProject => $BioProject,  
		runpath => $runpath,
		# TSADESC=>$TSADESC, MINSIZE=>$MINSIZE, MAXGAP=>$MAXGAP,
		# trclass => $trclass, mrna => $cdnaseq, genenames=>$genenames,
		);
	
  if($action =~ /restore|read/ and -f $trpname) {
    open(my $inh, $trpname); # or loggit(warn ..);
    while(<$inh>) { chomp; if(s/^$PRES.//) { 
      my($k,$v)=split /\s*=\s*/,$_,2; 
      $v=~s/;*$//; $v=~s/\s*\#.*$//; # should I chomp trailing '#' comments? 
      my $ov= $mysettings{$k}; 
      unless($ov and $ov ne $DEFAULT_SETTINGS{$k}) { 
        ## fixme: need to reset global defaults
        $organism=  $v if($k eq 'organism');
        do { $v=~s/ +;/;/g; $v=~s/ +/;/g; $sraids= $v; } if($k eq 'sraids');
        $BioProject= $v if($k eq 'BioProject'); # FIXME lost
        $IDPREFIX=  $v if($k eq 'IDPREFIX');
        $DATE=      $v if($k eq 'DATE');
        # $TSADESC=   $v if($k eq 'TSADESC');
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
    @sraids= grep /^[A-Z]+\d+/, split /\W+/, $sraids; # NOW GLOBAL; should be SRRnnnn ERRnnnn list

  } else { # action==save; copy new single vars to settings..
    for my $ks (sort keys %mysettings) { 
      my $ov= $mysettings{$ks};  
      $settings{$ks}= $ov if($ov and $ov ne $DEFAULT_SETTINGS{$ks});    
      } ;
  }

  my $settings= join "\n", map{ "$PRES.$_=".$settings{$_} } sort keys %settings;
  if($action =~ /log/) { loggit(0, "sra2genes.info:\n$settings");  }
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
      elsif(/^(SRR|ERR)\d+,/){ $cvsformat=2; 
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
  $sorg=~s/ /_/g;
  $organism= $sorg if($sorg and ($organism eq $deforg or $organism !~ m/\w\w/));
  $sraids=   $sids if($sids and ($sraids eq $defsra or $sraids !~ m/SRR/)); #? or always use sids?
  @sraids= grep /^[A-Z]+\d+/, split /\W+/, $sraids; # NOW GLOBAL; should be SRRnnnn ERRnnnn list
  
  if($sorg) {
    my($gn,$sp)= split/[_\W]/,$organism,2;
    $settings{oname}= $organism . 'EVd';
    if($sp and $gn) {
      my $shortorg= ucfirst(substr($gn,0,3)) . lc(substr($sp,0,3));
      $shortorg .= 'EVd';
      $settings{oname}= $shortorg;
      $IDPREFIX= $shortorg if($IDPREFIX eq $DEFAULTidpre);
      }
  }
  ## also use sradatah to edit template evgr_tsamethods.cmt, evgr_tsadesc.cmt, evgr_tsasubmit.sbt ??
  ## rewrite template .cmt, .sbt unless updated already.

  ## ALSO use 'FTP Path to Experiment' to turn SRX into SRR for picky ncbi tsa : wget or curl calls?
  my $HAVEsrrids= ($sraids =~ /SRR|ERR/ and $sraids ne $defsra)?1:0; ## Need to save to info/config file
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

sub get_srainfo {
  my($sracvs,$sraidlist)= @_;
	#o# my($trpath,$trname)= @_;
	
	# my $sracvs="$trpath/$trname.sra_result.csv";    
	# $sracvs="$trpath/sra_result.csv" unless(-f $sracvs);  
	my @SRAK=();
	
	# if($sraidlist and not $sracvs) { ... } # get sracvs from ncbi given ids?
	my($nsra,$sradatah,$gotsraids)= parse_sra_result_cvs($sracvs) if(-f $sracvs);
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
  	@SRAK= ("Assemblers", "Platform", "size_MB","spots",'Total Assemblies', 'BioProject');
	}
	for my $ks (@SRAK) {
  	if($settings{$ks} and not $sradatah->{$ks}) { $sradatah->{$ks} = $settings{$ks}; }
  	elsif($sradatah->{$ks}) { $settings{$ks} = $sradatah->{$ks}; }
	}
	
	 ## also use sradatah to edit template evgr_tsamethods.cmt, evgr_tsadesc.cmt, evgr_tsasubmit.sbt ??
	# my($nupinfo,$tsamethf,$tsadescf,$tsasubf)= tsa_infotemplates($trpath, $trname, $sradatah);
	# loggit(0,"info updated $nupinfo TSADESC=",$TSADESC);
 
  # unless($genenames) { #? put in get_evgtrset
  #   my $gnt="$trpath/$trname.names";  $genenames=$gnt if(-s $gnt); 
  #   loggit(LOG_WARN, "Missing product -names",$gnt) unless($genenames);
  # }
	
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


sub splitspots {
  my($fn,$sfa)=@_; 
  my($id,$err,$mlen,$nr,$tlen)=(0,0,0,0,0);
  my $lfn=$fn."_1.fa"; my $rfn=$fn.="_2.fa";
  open(S,$sfa) or return -1;
  open(L,'>',$lfn) or return -1; 
  open(R,'>',$rfn) or return -1;
  while(<S>) {
    if(/^>(\S+)/) { $id=$1;  $nr++; } else { chomp; 
    my $len=length($_); my $hl=int($len/2); 
    $tlen+=$len; $mlen=$hl if($hl>$mlen);
    print L ">$id/1\n",substr($_,0,$hl),"\n"; 
    print R ">$id/2\n",substr($_,$hl),"\n"; 
    } 
  } close(L); close(R); close(S);
  return($err,$nr,$mlen,$tlen,$lfn,$rfn);
  # return 0; # return max read size?
}

sub get_readpair_info {
  my($pairdir)=@_; # == pairfa
  my %PAIRINFO=();
  my ($pall,@pinfo)= getFileset('pairfa','fa.info$');  

  # FIXME: grep @sraids, maybe more in pairfa than want..
  # my $sraidset= join '|',@sraids;
  # @pinfo= grep /$sraidset/, @pinfo;

  for my $pinfo (@pinfo) { 
    (my $pn=$pinfo)=~s/\.fa.info//;
    my ($sid)= $pinfo=~m/(\w+)\.fa.info/; # *should be ok, maybe not..
    # messy, both pathto/SID and SID keys, drop pathto/ keys
    
    $PAIRINFO{$sid}{nreads}=0; 
    if(open(FI,$pinfo)) { my $finfo=<FI>; close(FI); chomp($finfo); 
      map{ my($k,$v)= split"="; $PAIRINFO{$sid}{$k}= $v; } split";",$finfo;
    }
  }
  return %PAIRINFO;
}
  
sub fa2pairfa { 
  # evigene/scripts/rnaseq/fa2pairfa.pl == interleave pair1/2 files as one.fa2 file
  my($fain,$fa2in)= @_;
  my($ok,$nin,$nskip,$nerr,$npair)= (0) x 9;
  my $paired=1; my $addnum=1; my $isFQ=0; # fixed opts
  (my $fname=$fain) =~ s/.gz//; 
  $fname=~s/\.(fast[aq]|f[aq])//;
  $fname=~s/_\d$//; 
  my $outf= "$fname.fa2";
  return($outf,STEPdone) if( -s $outf );
  loggit(0,"interleave $outf from $fain,$fa2in");
  
  open(OUT,">$outf") or loggit(LOG_DIE,"writing $outf"); # loggit(LOG_DIE,..)
  if($fain=~/.gz/) { $ok=open(F,"gunzip -c $fain |"); } else { $ok=open(F,$fain); }
  if($fa2in=~/.gz/) { $ok=open(F2,"gunzip -c $fa2in |"); } else { $ok=open(F2,$fa2in); }
  
  my($fh,$fs,$lh,$fhr,$fsr,$rh, $qh,$qs, $qhr,$qsr)=("") x 10;
  while(<F>) { 
    $fh=$_; $fs=<F>; ($lh)= $fh=~/^.(\S+)/; $nin++;   
    $fhr=<F2>; $fsr=<F2>; ($rh)= $fhr=~/^.(\S+)/;   
    # if($isFQ) { # not yet
    #   $fh=~s/^./>/;  $qh=<F>; $qs=<F>;  # read/drop  # fh=='@header'
    #   $fhr=~s/^./>/; $qhr=<F2>; $qsr=<F2>;  # read/drop
    #   }
    if($addnum) {
      unless($fh=~m,/1 ,) { $fh=~s, ,/1 ,; $fhr=~s, ,/2 ,, unless($fhr=~m,/2 ,); }
    }
    if($lh ne $rh) { 
      $lh=~s,/[12],,; $rh=~s,/[12],,; 
      if($lh ne $rh) { $nerr++;   
        warn "# mismatch pair: $lh ne $rh\n";  
        loggit(LOG_DIE, "ERR: too many mismatch:$nerr") if($nerr>99); 
        # die "ERR: too many mismatch:$nerr" if($nerr>99); 
        #^^ loggit(LOG_DIE,...);
        next;
# mismatch pair: SRR3247516.99 ne SRR3247180.99
# mismatch pair: SRR3247516.100 ne SRR3247180.100
# ERR: too many mismatch:100 at /home/ux455375/bio/evigene/scripts/evgpipe_sra2genes.pl line 1158, <F2> line 200.
# [comet-02-09.sdsc.edu:mpispawn_0][child_handler] MPI process (rank: 0, pid: 26403) exited with status 25

      } 
    }
  print OUT $fh,$fs; 
  print OUT $fhr,$fsr; $npair++;
  }
  close(OUT); close(F); close(F2) if($paired);
  return ($npair>0)?($outf,STEPok):($outf,STEPerr);
}


#============= Gene Assembler Template Scripts =======================

sub add_script {
  my($ok,$sh,$outsh,$asmNAME,$stepname)= @_;
  # my($ok,$sh,$outsh)= SCRIPT_velvetoases();
  if($ok>0) {
    open(S,'>',$outsh); print S $sh; close(S);
    system("chmod +x $outsh"); #? or not
    $settings{$stepname} .="$outsh;"; # STEP4_runassemblers
    loggit(0,"wrote $asmNAME to $outsh"); # save outsh in settings ...
    return(1,$asmNAME);
    # $whichasm.="$asmNAME,"; #? "velvetoases,";
  } else {
    loggit(LOG_WARN,"error $asmNAME"); 
    return(0);
  }
}

sub SCRIPT_common {
  my($runapp)=@_;
  
  my $P_DATAD = $settings{DATAD} || $settings{runpath} || `pwd`; # FIXME; use settings{DATAD} instead?
  my $P_NAME  = $settings{readset} || $runname; # asm run name
  $P_NAME= basename($P_NAME); # $P_NAME =~ s,^.*/,,g;
  
  ## use/require bioproject or equiv as asm run name
  ##   my $outna="srr$samplen". ($sradatah->{BioProject} || $IDPREFIX );
  ##  $settings{readset}  = $sradatah->{readset} = $outna;

  chomp($P_DATAD);
  # replace script vars here, maybe all vars have P_ prefix in template scripts?
  map{ 
    s/P_EVIGENES/$EVIGENES/g;  s/P_NCPU/$NCPU/g; s/P_MAXMEM/$MAXMEM/g; # globals
    s/P_DATAD/$P_DATAD/g;   s/P_NAME/$P_NAME/g; 
    #maybe# s/P_APP/$P_APP/g; s/P_ABIN/$P_ABIN/g; 
  } ($runapp);
  
  # check, return all P_ vars in script for caller
  my @pvar= $runapp =~ m/\b(P_\w+)/g; 
  my %pvar= map{ $_ => 1} @pvar;
  @pvar= sort keys %pvar;
  ## check in settings{var} ??
  # return($runapp, @pvar);
  
  return($runapp);
}


sub SCRIPT_idba_trans {
my $SCRIPT_idba_trans = <<'EOS'; # dont allow perl $var subs
#! /bin/bash
## env name=outname inpe=readpairs.fa datad=path/to/data this.sh
#PBS -N idba_trans
#PBS -A PutAccountIdHere
#PBS -l nodes=1:ppn=16,walltime=39:55:00
#PBS -V

if [ "X" = "X$ncpu" ]; then ncpu=P_NCPU; fi
if [ "X" = "X$maxmem" ]; then maxmem=P_MAXMEM; fi
if [ "X" = "X$datad" ]; then datad=P_DATAD; fi
if [ "X" = "X$inpe" ]; then inpe='P_INPE'; fi
if [ "X" = "X$name" ]; then name=P_NAME; fi

dv="1a"; # versnum, add opt
traopts="--mink P_KMIN --maxk P_KMAX --step P_STEPS --max_isoforms 9";

evigenes=P_EVIGENES
idba_tran=P_APP
export PATH=P_ABIN:$PATH

#o# outdir=tridba$name$dv
outdir=tridba${dv}_$name
cd $datad/

echo "#START  `date`"
echo "$idba_tran $traopts --read $inpe --num_threads $ncpu --out $outdir"
$idba_tran $traopts --read $inpe --num_threads $ncpu --out $outdir

## add cleanups, where contigs.fa  is ok .. all graph files?
## gzip transcripts,contigs,other keepers?

echo "#DONE : `date`"
EOS

  my $runapp=$SCRIPT_idba_trans;
  my $P_DATAD = $settings{runpath} || `pwd`; # in SCRIPT_common but use below
  my $P_INPE  = $settings{readset}; # lacks _[12].fa suffix ??

  # need better file names: 
  #my $runfile="runidba_$P_INPE.sh" ; my $outdir= "traidba_$P_INPE";
  my $runfile= "runidba.sh";
  $runapp= SCRIPT_common($runapp);
  
  ## FIXME: idba needs fa2, interleaved, not 12.fa pair files
  ## my $P_INPE  = $settings{readset_ileave}; # reads.fa2 ??
  my $okf=0;
  my $fa1= $P_INPE."_1.fa"; # assume ok
  my $fa2= $P_INPE."_2.fa";
  my $fa12= $P_INPE.".fa2"; # maybe not there, could be from velveth seq/xxx.fa
  $okf= (-s $fa12)?1:0;
  unless($okf>0) { ($fa12,$okf)= fa2pairfa($fa1,$fa2); }
  unless($fa12 =~ m,/,) { my $ft="$P_DATAD/$fa12"; $fa12=$ft if(-s $ft); }   
  unless($okf>0){ loggit(1,"data error $fa12"); return(STEPerr) ; }
  $P_INPE= $fa12;

  my $readsize= $settings{readsize}||100; # $sradatah->?
  my $P_KMIN=27;
  my $P_STEPS=10;
  my $P_KMAX=$readsize - 9;
  for(my $k=$P_KMIN+$P_STEPS; $k<$readsize; $k+=$P_STEPS) { $P_KMAX=$k; }
  
  # idba, velv, need alternate paths to diff compiled configs (ie max read size)
  my($P_APP,$P_ABIN)= findapp('idba_tran', 1);
  return(STEPerr) if($P_APP =~ /MISSING/); # { loggit(1,"app error $P_APP");  ; }

  ##nono# my $rundir= $P_DATAD."/tridba";
  ##nono# mkdir($rundir) or return(STEPerr);
  ##Nono# $P_DATAD= $rundir; #?? this may mess up data path to P_INPE
  ## INPE should have full path anyway?
  
  # replace script vars here
  map{ 
    #common# s/EVIGENES/$EVIGENES/g;  s/NCPU/$NCPU/g; s/MAXMEM/$MAXMEM/g; 
    #common# s/P_DATAD/$P_DATAD/g; 
    s/P_APP/$P_APP/g; s/P_ABIN/$P_ABIN/g; 
    s/P_INPE/$P_INPE/g; 
    s/P_KMIN/$P_KMIN/g; s/P_KMAX/$P_KMAX/g; s/P_STEPS/$P_STEPS/g; 
  } ($runapp);
  
  return(STEPok,$runapp,$runfile);# ,$rundir
}


=item IDBA-Tran help

 IDBA-Tran - Iterative de Bruijn Graph Assembler for next-generation transcriptome sequencing data.
 Usage: idba_tran -r read.fa -o output_dir
 Allowed Options: 
   -o, --out arg (=out)                   output directory
   -r, --read arg                         fasta read file (<=128)
   -l, --long_read arg                    fasta long read file (>128)
       --mink arg (=20)                   minimum k value (<=124)
       --maxk arg (=60)                   maximum k value (<=124)
       --step arg (=10)                   increment of k-mer of each iteration
       --inner_mink arg (=10)             inner minimum k value
       --inner_step arg (=5)              inner increment of k-mer
       --prefix arg (=3)                  prefix length used to build sub k-mer table
       --min_count arg (=2)               minimum multiplicity for filtering k-mer when building the graph
       --min_support arg (=1)             minimum supoort in each iteration
       --num_threads arg (=0)             number of threads
       --seed_kmer arg (=30)              seed kmer size for alignment
       --min_contig arg (=200)            minimum size of contig
       --similar arg (=0.95)              similarity for alignment
       --max_mismatch arg (=3)            max mismatch of error correction
       --no_local                         do not use local assembly
       --no_coverage                      do not iterate on coverage
       --no_correct                       do not do correction
       --pre_correction                   perform pre-correction before assembly
       --max_isoforms arg (=3)            maximum number of isoforms
       --max_component_size arg (=30)     maximum size of components
       
=cut

sub SCRIPT_velvetoases {
my $SCRIPT_velvetoases = <<'EOS'; # dont allow perl $var subs
#! /bin/bash
## env name=outname inpe=readpairs.fa datad=path/to/data this.sh
#PBS -N velvetoases
#PBS -A PutAccountIdHere
#PBS -l nodes=1:ppn=16,walltime=39:55:00
#PBS -V

# more than 8 cpu bad, overuses memory ..
if [ "X" = "X$ncpu" ]; then ncpu=P_NCPU; fi
if [ "X" = "X$maxmem" ]; then maxmem=P_MAXMEM; fi
if [ "X" = "X$datad" ]; then datad=P_DATAD; fi
if [ "X" = "X$inpe" ]; then inpe='P_INPE'; fi
if [ "X" = "X$name" ]; then name=P_NAME; fi

dv="1a"; 
outdir=trvelo${dv}_$name
INSIZE=P_INSIZE #or 200
INSIZESD=50 # or P_INSIZESD  

# kset_example="95 85 75 71 65 63 55 45 35"
kset="P_KSET"

evigenes=P_EVIGENES
velbin=P_ABIN
# FIXME appbin variants...
# #k>99 bin4
# velbin4=$HOME/bio/velvet1210/bin4
# #k>61# 
# velbin2=$HOME/bio/velvet1210/bin2
# #k<=61# 
# velbin1=$HOME/bio/velvet1210/bin
# # velbin=$velbin4; 
##......................

export OMP_NUM_THREADS=$ncpu

vopth=""
iopts="-ins_length $INSIZE -ins_length_sd $INSIZESD"
vopts="$iopts -max_gap_count 5 " 
oopts="-scaffolding yes -min_pair_count 3 -edgeFractionCutoff 0.03 -min_trans_lgth 200 $iopts"

echo "#START `date` " 
cd $datad
mkdir $outdir
cd $outdir

#..... run loop for vel kmer steps
shopt -s nullglob

kseqdir=vel${dv}_seq
if [ ! -f $kseqdir/Sequences ]; then
 $velbin/velveth $kseqdir 27 $vopth -noHash -shortPaired -fmtAuto -separate $inpe
fi

for k in $kset;  do { 
  ksubdir=vel${dv}_$k
  echo "#.. start velrun $ksubdir : `date`"
  #FIXME# velbin=$velbin4; 
  #FIXME# if [ $k -le 99 ]; then velbin=$velbin2; fi
  #FIXME# if [ $k -le 61 ]; then velbin=$velbin1; fi
  mkdir $ksubdir
  ln -s ../$kseqdir/Sequences $ksubdir/
  $velbin/velveth $ksubdir $k $vopth -reuse_Sequences
  $velbin/velvetg $ksubdir $vopts -read_trkg yes 

  echo "#.. end velvetg $ksubdir : `date`"
} done
## end loop 1

## loop 2 oases, 1cpu; save memory, oases doesnt use multicpu
export OMP_NUM_THREADS=2
export OMP_THREAD_LIMIT=2
i=0; 
for k in $kset;  do {
  ksubdir=vel${dv}_$k
  if [ -f $ksubdir/transcripts.fa ]; then continue; fi
  if [ -f $ksubdir/contigs.fa ]; then
    ## DONT fork here save mem ..
    #FIXME# velbin=$velbin4;
    #FIXME# if [ $k -le 99 ]; then velbin=$velbin2; fi
    #FIXME# if [ $k -le 61 ]; then velbin=$velbin1; fi
    $velbin/oases   $ksubdir $oopts 
    i=$(( $i + 1 ))
  fi
} done

wait

## add cleanups, where transcripts.fa  is ok .. all Graph files?
## gzip transcripts,contigs,other keepers?
for k in $kset;  do {
  ksubdir=vel${dv}_$k
  if [ -f $ksubdir/transcripts.fa ]; then
    /bin/rm $ksubdir/{Graph2,LastGraph,PreGraph,Roadmaps}
    gzip --fast $ksubdir/{contig-ordering.txt,contigs.fa,stats.txt}
    # gzip $ksubdir/transcripts.fa later
  fi
} done


echo "#DONE : `date`"
EOS

  my $runapp=$SCRIPT_velvetoases;
  my $P_DATAD = $settings{runpath} || `pwd`; # in SCRIPT_common but use below
  my $P_INPE  = $settings{readset}; # lacks _[12].fa suffix ??

  # need better file names: 
  # my $runfile="runvelo_$P_INPE.sh" ; my $outdir= "travelo_$P_INPE";
  my $runfile= "runvelo.sh";
  $runapp= SCRIPT_common($runapp);
  
  my $P_INSIZE= $settings{insertsize} || 200; # P_INSIZESD ?

  my $okf=0;
  my $fa1= $P_INPE."_1.fa";  
  my $fa2= $P_INPE."_2.fa";
  $okf= (-s $fa1 and -s $fa2)?1:0;
  unless($okf>0){ loggit(1,"data error $fa1"); return(STEPerr) ; }
  #x if($runapp =~ /\bcd\b|\bchdir/){ .. }
  map { $_="../$_" unless(m,/,) } ($fa1,$fa2); # find from run subdir this way?
  # bug fixed: inpe='../srr2PRJNA315720_1.fa srr2PRJNA315720_2.fa'; # <<< MISSING ../_2.fa
  
  $P_INPE= "$fa1 $fa2"; #? okay both ; need quotes?

  my $readsize= $settings{readsize}||100; # $sradatah->?
  my $P_KMIN= 31;
  my $P_STEPS=10;
  ## my $P_KMAX= $readsize - 9;
  my $P_KMAX= (int($readsize/10) -1 )* 10 + 5;
  my $halfr1= int($readsize/2) - 6; my $halfr2= $halfr1 + 12;
  my $P_KSET=""; for(my $k=$P_KMAX; $k>=$P_KMIN; $k-=$P_STEPS) {  
    $P_KSET.="$k "; 
    if($k > $halfr1 and $k < $halfr2) { my $k3=$k+4; $P_KSET.="$k3 "; } #my $k2=$k-4;  $P_KSET.="$k3 $k2 "; 
    }
  #result: kset=95 85 75 65 55 59 51 45 49 41 35 : 
  # halfsteps at midpt: 59<55>51,49<45>41 ? maybe want just 65, 59,55, 49,45, 35
  
  # idba, velv, need alternate paths to diff compiled configs (ie max read size)
  my($P_APP,$P_ABIN)= findapp('velveth', 1);
  return(STEPerr) if($P_APP =~ /MISSING/); #{ loggit(1,"app error $P_APP");   ; }

  ##nono# my $rundir= $P_DATAD."/tridba";
  ##nono# mkdir($rundir) or return(STEPerr);
  ##Nono# $P_DATAD= $rundir; #?? this may mess up data path to P_INPE
  ## INPE should have full path anyway?
  
  # replace script vars here
  map{ 
    #common# s/EVIGENES/$EVIGENES/g;  s/NCPU/$NCPU/g; s/MAXMEM/$MAXMEM/g; 
    #common# s/P_DATAD/$P_DATAD/g; 
    s/P_APP/$P_APP/g; s/P_ABIN/$P_ABIN/g; 
    s/P_INPE/$P_INPE/g; 
    s/P_INSIZE/$P_INSIZE/g; 
    s/P_KSET/$P_KSET/g;  
  } ($runapp);
  
  return(STEPok,$runapp,$runfile); 
}

sub SCRIPT_soap_trans {
my $SCRIPT_soap_trans = << 'EOS';
#! /bin/bash
## env name=outname inpe=readpairs.fa datad=path/to/data this.sh
#PBS -N soap_trans
#PBS -A PutAccountIdHere
#PBS -l nodes=1:ppn=16,walltime=39:55:00
#PBS -V

if [ "X" = "X$ncpu" ]; then ncpu=P_NCPU; fi
if [ "X" = "X$maxmem" ]; then maxmem=P_MAXMEM; fi
if [ "X" = "X$datad" ]; then datad=P_DATAD; fi
if [ "X" = "X$inpe" ]; then inpe='P_INPE'; fi
if [ "X" = "X$name" ]; then name=P_NAME; fi

dv="1a"; 
outdir=trsoap${dv}_$name
onamep=$name$dv

INSIZE=P_INSIZE #or 200
INSIZESD=50 # or P_INSIZESD  
# F: gapfil; M: merge strength 1=def
sopt="-p $ncpu -F -t 99"

kset_example="31 27 25 35 45 55"
kset="P_KSET"

evigenes=P_EVIGENES
export PATH=P_ABIN:$PATH

cd $datad/
mkdir $outdir

# write soap.config .. inpe presumed local path rel to datad
echo "max_rd_len=P_READSIZE" > $outdir/rdconfig 
for fa1 in $inpe; do {
  fa2=`echo $fa1 | sed 's/_1/_2/;'`
  if [ $fa2 = $fa1 ]; then continue; fi
  # nogz=`echo $fa1 | sed 's/\.gz//;'` # skip this, require ungzipped data
  cat >> $outdir/rdconfig <<EOT

[LIB]
rank=1
asm_flag=3
avg_ins=$INSIZE
f1=$datad/$fa1
f2=$datad/$fa2

EOT

} done

cd $outdir
echo "#START `date` " 

for kmer in $kset; do {
  odir=sod$kmer
  outname=so${onamep}.k$kmer
  mkdir $odir
  if [ $kmer -lt 32 ]; then
    SOAPdenovo-Trans-31mer all -s rdconfig -o $outname -K $kmer $sopt
  else 
    SOAPdenovo-Trans-127mer all -s rdconfig -o $outname -K $kmer $sopt
  fi

  mv $outname* $odir/
  if test -f $odir/$outname.contig ; then
    rm $odir/$outname.{readOnContig,readInGap,ctg2Read,preArc,vertex,edge.gz,Arc,updated.edge}
    gzip --fast $odir/$outname.{scafSeq,contig}
  fi

} done

echo "#DONE : `date`"
EOS

  my $runapp=$SCRIPT_soap_trans;
  my $P_DATAD = $settings{runpath} || `pwd`; # in SCRIPT_common but use below
  my $P_INPE  = $settings{readset}; # lacks _[12].fa suffix ??

  # need better file names: 
  # my $runfile="runsoap_$P_INPE.sh" ; my $outdir= "trasoap_$P_INPE";
  my $runfile= "runsoap.sh"; # change to match outdir; runs= run_$rname.sh, outs=tra_$rname/
  $runapp= SCRIPT_common($runapp);

  my $P_DATAD = $settings{runpath} || `pwd`; # FIXME
  my $P_INPE  = $settings{readset}; # lacks _[12].fa suffix ??
  my $P_INSIZE= $settings{insertsize} || 200; # P_INSIZESD ?

  # idba, velv, need alternate paths to diff compiled configs (ie max read size)
  ## I have: SOAPdenovo-Trans-127mer SOAPdenovo-Trans-31kmer ... use those?
  my($P_APP,$P_ABIN)= findapp('SOAPdenovo-Trans-127mer', 1);
  return(STEPerr) if($P_APP =~ /MISSING/); # dup{ loggit(1,"app error $P_APP");  ; }

  my $okf=0;
  my $fa1= $P_INPE."_1.fa";  
  my $fa2= $P_INPE."_2.fa";
  $okf= (-s $fa1 and -s $fa2)?1:0;
  unless($okf>0){ loggit(1,"data error $fa1"); return(STEPerr) ; }
  map { $_="../$_" unless(m,/,) } ($fa1,$fa2); # find from run subdir this way?
  $P_INPE= $fa1; #? soap script adds _2.fa .. change?

  # kset_example="31 27 25 35 45 55"
  my $readsize= $settings{readsize}||100; # $sradatah->?
  my $P_KMIN= 25;
  my $P_STEPS= 6;
  my $P_KMAX= _min1(75, (int($readsize/10) -1 )* 10 + 5);
  #x my $halfr1= int($readsize/2) - 6; my $halfr2= $halfr1 + 12;
  my $P_KSET=""; for(my $k=$P_KMIN; $k<=$P_KMAX; $k+=$P_STEPS) {  
    $P_KSET.="$k "; 
    #x if($k > $halfr1 and $k < $halfr2) { my $k3=$k+4; $P_KSET.="$k3 "; }   
    }

  # replace script vars here
  map{ 
    s/P_APP/$P_APP/g; s/P_ABIN/$P_ABIN/g; 
    s/P_INPE/$P_INPE/g; 
    s/P_INSIZE/$P_INSIZE/g; 
    s/P_KSET/$P_KSET/g;  
    s/P_READSIZE/$readsize/g; 
  } ($runapp);
  
  return(STEPok,$runapp,$runfile); 
}

sub SCRIPT_trinity {
my $SCRIPT_trinity = << 'EOS';
#! /bin/bash
## env name=outname inpe=readpairs.fa datad=path/to/data this.sh
#PBS -N trinity
#PBS -A PutAccountIdHere
#PBS -l nodes=1:ppn=16,walltime=39:55:00
#PBS -V

if [ "X" = "X$ncpu" ]; then ncpu=P_NCPU; fi
if [ "X" = "X$maxmem" ]; then maxmem=P_MAXMEM; fi
if [ "X" = "X$datad" ]; then datad=P_DATAD; fi
if [ "X" = "X$inpe" ]; then inpe='P_INPE'; fi
if [ "X" = "X$name" ]; then name=`basename $inpe .fa2 | sed 's/\.fa.*$//'`; fi

dv="1a"; 

evigenes=P_EVIGENES
export PATH=P_ABIN:$PATH

outdir=trtrin${dv}_$name
cd $datad/
mkdir $outdir
cd $outdir

echo "#START `date` " 

echo "#DONE : `date`"
EOS

  my $runapp=$SCRIPT_trinity;
  my $P_DATAD = $settings{runpath} || `pwd`; # in SCRIPT_common but use below
  my $P_INPE  = $settings{readset}; # lacks _[12].fa suffix ??

  # need better file names: 
  #my $runfile="runtrin_$P_INPE.sh" ; my $outdir= "tratrin_$P_INPE";
  my $runfile= "runtrin.sh";
  $runapp= SCRIPT_common($runapp);
  
  # idba, velv, need alternate paths to diff compiled configs (ie max read size)
  my($P_APP,$P_ABIN)= findapp('Trinity', 1);
  return(STEPerr) if($P_APP =~ /MISSING/); # { loggit(1,"app error $P_APP");   ; }

  # replace script vars here
  map{ 
    s/P_APP/$P_APP/g; s/P_ABIN/$P_ABIN/g; 
    s/P_INPE/$P_INPE/g; 
    # s/P_INSIZE/$P_INSIZE/g; 
    # s/P_KSET/$P_KSET/g;  
  } ($runapp);
  
  return(STEPok,$runapp,$runfile); 
}

sub SCRIPT_diginorm {
my $SCRIPT_diginorm = << 'EOS';
#! /bin/bash
## env name=outname inpe=readpairs.fa datad=path/to/data this.sh
#PBS -N diginorm
#PBS -A PutAccountIdHere
#PBS -l nodes=1:ppn=16,walltime=39:55:00
#PBS -V

if [ "X" = "X$ncpu" ]; then ncpu=P_NCPU; fi
if [ "X" = "X$maxmem" ]; then maxmem=P_MAXMEM; fi
if [ "X" = "X$datad" ]; then datad=P_DATAD; fi
if [ "X" = "X$inspot" ]; then inspot='P_INSPOT'; fi
# options: what defaults?
if [ "X" = "X$kmer" ]; then kmer=25; fi
if [ "X" = "X$keep" ]; then keep=20; fi

dv="1a"; 

#** requires python2.6+ .. need this?
module add python

evigenes=P_EVIGENES

## P_ABIN == $HOME/bio/khmer ? or $HOME/bio/khmer/scripts
export PATH=P_ABIN:$PATH
export PYTHONPATH=P_ABIN/../python
normapp=P_ABIN/normalize-by-median.py

#  ncpu=2 # this is 1-cpu app, I think, may have changed.
kopt=""
outdir=rnasets

# FIXME: calc xhash from MAXMEM, for data size, along with $kmer and $keep
# 8e9 x 4 = 32Gb ; 12e9 x 4 = 48Gb
xhash=8e9

# output to rnasets/ folder of assembler inputs? or need tmp folder?
cd $datad/
mkdir $outdir
cd $outdir

echo "#START `date` " 

## inspot can be list..
echo $normapp $kopt -C $keep -x $xhash  -k $kmer  $inspot
$normapp $kopt -C $keep -x $xhash  -k $kmer  $inspot

echo "#DONE : `date`"
EOS

  my $runapp=$SCRIPT_diginorm;
  my $P_DATAD = $settings{runpath} || `pwd`; # in SCRIPT_common but use below
  my ($spotall,@spofa)= getFileset('spotfa','fasta$');  

  # FIXME: grep @sraids, maybe more in pairfa than want..
  my $sraidset= join '|',@sraids;
  @spofa= grep /$sraidset/, @spofa;
  
  map{ s,^,../, } @spofa; # subdir to subdir
  my $P_INSPOT=join" ",@spofa;

=item diginorm
  
  #^^ run on spotfa/*.fasta not pairfa/..; can do all spotfa/*.fasta
  #  need splitspots on dnorm result after rundignorm.sh
# sub dignorm_fini { my($name,$sfa)= ($name,$dnormout)...
#    my($cerr,$nreads,$maxlen,$totlen,$lfn,$rfn)= splitspots($name,$sfa);
#   if($cerr == 0) {
#   open(FI,'>',"$name.fa.info"); 
#   print FI "nreads=$nreads;maxlen=$maxlen;totlen=$totlen;lfn=$lfn;rfn=$rfn\n"; 
#   close(FI);

=cut
  
  # need better file names: 
  my $runfile="rundiginorm_$runname.sh" ; my $outdir= "rnasets";
  
  $runapp= SCRIPT_common($runapp); # let common() set runfile name from P_NAME?
  
  my($P_APP,$P_ABIN)= findapp('normalize-by-median.py', 1);
  return(STEPerr) if($P_APP =~ /MISSING/); 

  # replace script vars here
  map{ 
    s/P_APP/$P_APP/g; s/P_ABIN/$P_ABIN/g; 
    s/P_INSPOT/$P_INSPOT/g; 
    # s/P_KMER/$P_KMER/g;  s/P_KEEP/$P_KEEP/g; 
  } ($runapp);
  
  return(STEPok,$runapp,$runfile); 
}

sub SCRIPT_evgpubset {
my $SCRIPT_evgpubset = << 'EOS';
#! /bin/bash
### env idprefix=MysppEGm trclass=myspp.trclass datad=`pwd` qsub -q shared run_evgmrna2tsa.sh
#PBS -N evgpubset
#PBS -A PutAccountIdHere
#PBS -l nodes=1:ppn=16,walltime=39:55:00
#PBS -V

if [ "X" = "X$ncpu" ]; then ncpu=P_NCPU; fi
if [ "X" = "X$maxmem" ]; then maxmem=P_MAXMEM; fi
if [ "X" = "X$datad" ]; then datad=P_DATAD; fi
if [ "X" = "X$trclass" ]; then trclass=P_TRCLASS; fi
if [ "X" = "X$species" ]; then species=P_ORGANISM; fi
if [ "X" = "X$idprefix" ]; then idprefix=P_IDPREFIX; fi

export evigenes=P_EVIGENES
export EVIGENES=P_EVIGENES
# export vecscreen=P_NCBIBIN/vecscreen
# export tbl2asn=P_NCBIBIN/tbl2asn

opts="-debug  -dropshow -skipdropseq -NCPU $ncpu"

## .. these are now read via sra_result.csv, species => idprefix
if [ "X" != "X$idprefix" ]; then opts="$opts -idprefix $idprefix"; fi
if [ "X" = "X$vectrim" ]; then opts="$opts -novectrim"; fi
if [ "X" != "X$tbl2asn" ]; then opts="$opts -runtbl2asn"; fi
if [ "X" != "X$species" ]; then spp=`echo $species | sed 's/ /_/g;'`; opts="$opts -species=$spp"; fi

if [ "X" = "X$names" ]; then names=`echo $trclass | sed 's/.gz//; s/\.trclass/.names/;'`; fi
if [ -s $names ]; then opts="$opts -names $names"; fi

cd $datad/

echo "#START `date` " 
echo $evigenes/evgmrna2tsa2.pl  $opts -log -class $trclass
if [ ! -f $trclass ]; then echo "ERR: missing -class $trclass"; exit -1; fi
$evigenes/evgmrna2tsa2.pl  $opts -log -class $trclass
echo "#DONE : `date`"

EOS

  my $runapp=$SCRIPT_evgpubset;
  my $runfile= "run_evgpubset.sh";
  $runapp= SCRIPT_common($runapp);
  
  # param: P_TRCLASS  P_NCBIBIN P_ORGANISM P_IDPREFIX
  my $trname= $runname; # 
  my $P_TRCLASS= $settings{trclass} || "$trname.trclass";
  unless(-s $P_TRCLASS) {
    loggit(1,"missing project.trclass:",$P_TRCLASS);
    return(STEPerr);
  }
  $settings{trclass} = $P_TRCLASS;
  
  my $P_ORGANISM=$organism;
  my $P_IDPREFIX=$IDPREFIX;
  
  # opt: -novectrim or -vectrim : need ncbi/vecscreen AND UniVec.db
  my($P_BLASTN,$P_NCBIBIN)= findapp('blastn', 1); #? dont need for pubset, but want vecscreen
  #x return(STEPerr) if($P_BLASTN =~ /MISSING/);

  # replace script vars here
  map{ 
    s/P_TRCLASS/$P_TRCLASS/g; 
    s/P_ORGANISM/$P_ORGANISM/g; 
    s/P_IDPREFIX/$P_IDPREFIX/g;  
    s/P_NCBIBIN/$P_NCBIBIN/g; 
  } ($runapp);
  
  return(STEPok,$runapp,$runfile); 
}

sub SCRIPT_tr2aacds {
my $SCRIPT_tr2aacds = << 'EOS';
#! /bin/bash
### env datad=path/to/data trset=myspecies_all.tr.gz qsub -q normal runtr2cds.sh
#PBS -N tr2aacds
#PBS -A PutAccountIdHere
#PBS -l nodes=1:ppn=16,walltime=39:55:00
#PBS -V

if [ "X" = "X$ncpu" ]; then ncpu=P_NCPU; fi
if [ "X" = "X$maxmem" ]; then maxmem=P_MAXMEM; fi
if [ "X" = "X$datad" ]; then datad=P_DATAD; fi
if [ "X" = "X$trset" ]; then trset='P_TRSET'; fi
if [ "X" = "X$name" ]; then name=`basename $trset .tr | sed 's/\.fa.*$//'`; fi

evigenes=P_EVIGENES
evapp=$evigenes/prot/tr2aacds2.pl
traopts="-tidy -log -debug"

# cd-hit-est/aa
# export PATH=$HOME/bio/cdhit/bin:$PATH
export PATH=P_CDHITBIN:$PATH
# export fastanrdb=$HOME/bio/exonerate/bin/fastanrdb
export fastanrdb=P_FASTANRDB
# blastn:
export PATH=P_NCBIBIN:$PATH

cd $datad/

echo "#START `date` " 
echo $evapp -NCPU $ncpu -MAXMEM $maxmem $traopts -cdna $trset
$evapp -NCPU $ncpu -MAXMEM $maxmem $traopts -cdna $trset
echo "#DONE : `date`"
EOS

  my $runapp=$SCRIPT_tr2aacds;
  my $runfile= "run_tr2aacds.sh";
  $runapp= SCRIPT_common($runapp);
  
  my $trname= $runname; # 
  #? my $trname= $settings{oname} || $runname;
  # $settings{trclass} = "$trname.trclass";

  ## BUG: this got trsets/tridbasrr2PRJNA3157201a.trformat.log logfile ** need 'suffix$' or make that default, ie 'suf.*' for otherwise
  # my($trall,@trs)= getFileset( "trsets",'tr|cdna|fasta'); 
  my($trall,@trs)= getFileset( "trsets",'tr$|cdna$|fasta$'); 
  return (STEPerr) unless(@trs);
  ## cdna_bestorf .aa,.cds in trsets/ can be used as input to tr2aacds ..
  
  # my $P_DATAD = $settings{runpath} || `pwd`; # FIXME
  my $P_TRSET = "$trname.tr"; # make all? or use trsets/*.tr ?
  unless(-s $P_TRSET) {
    my($runerr, $nok, $ofail)= cat_splitset($P_TRSET, \@trs);
    if($runerr) {  } # return(STEPerr)
  }
  
  #need paths: P_NCBIBIN P_FASTANRDB P_CDHITBIN  P_ABIN == bioappbin, may have others...
  my($P_FASTANRDB,$P_EXONRBIN)= findapp('fastanrdb', 1); 
  return(STEPerr) if($P_FASTANRDB =~ /MISSING/);
  my($P_CDHITEST,$P_CDHITBIN)= findapp('cd-hit-est', 1); 
  return(STEPerr) if($P_CDHITEST =~ /MISSING/);
  my($P_BLASTN,$P_NCBIBIN)= findapp('blastn', 1); 
  return(STEPerr) if($P_BLASTN =~ /MISSING/);

  # replace script vars here
  map{ 
    s/P_FASTANRDB/$P_FASTANRDB/g; s/P_CDHITBIN/$P_CDHITBIN/g; 
    s/P_NCBIBIN/$P_NCBIBIN/g;  
    s/P_TRSET/$P_TRSET/g;  s/P_TRNAME/$trname/g; # not P_NAME?
  } ($runapp);
  
  return(STEPok,$runapp,$runfile); 
}

sub SCRIPT_evgblastp {
my $SCRIPT_evgblastp = << 'EOS';
#! /bin/bash
### env aaset=my.aa refaa=ref.aa datad=path/to/data qsub -q normal evgblastp.sh
#PBS -N evgblast
#PBS -A PutAccountIdHere
#PBS -l nodes=1:ppn=16,walltime=39:55:00
#PBS -V

if [ "X" = "X$ncpu" ]; then ncpu=P_NCPU; fi
if [ "X" = "X$maxmem" ]; then maxmem=P_MAXMEM; fi
if [ "X" = "X$datad" ]; then datad=P_DATAD; fi
if [ "X" = "X$aaset" ]; then aaset=P_AASET; fi
if [ "X" = "X$refaa" ]; then refaa=P_REFAA; fi

evigenes=P_EVIGENES
# blastp/n/..:
export PATH=P_NCBIBIN:$PATH

cd $datad/

qname=`basename $aaset .aa`
refnam=`basename $refaa .aa`

blopt="-evalue 1e-5"
odir=blout1$qname
mkdir $odir

echo "#START `date` " 

if [ ! -f $refaa.psq ]; then
  echo makeblastdb -dbtype prot -in $refaa -logfile $refaa.mblog
  makeblastdb -dbtype prot -in $refaa -logfile $refaa.mblog
fi

if [ ! -f $aaset.split.1.fa ]; then
 pindir=`dirname $aaset`
 splitsize=`grep -v '^>' $aaset | wc -c | sed 's/ .*//' `
 splitbp=$(( $splitsize / $ncpu ))
 $evigenes/splitMfasta.pl --outputpath=$pindir --nparts $ncpu --minsize=$splitbp $aaset
fi

qset=`/bin/ls $aaset.split.*.fa`

for qfile in $qset
{
  qnamspl=`basename $qfile .fa`
  onam=$odir/$refnam-$qnamspl
  echo blastp $blopt -outfmt 7 -db $refaa -query $qfile -out $onam.blastp
  blastp $blopt -outfmt 7 -db $refaa -query $qfile -out $onam.blastp  &
}

wait

rqname=$refnam-$qname
aablast=$rqname.blastp
aabltab=$rqname.btall

cat $odir/$rqname.*.blastp > $aablast
/bin/rm $qset
echo $odir is temp blast output, check for completeness then erase
echo /bin/rm -r $odir

env oid=1 off=1 $evigenes/prot/aaqual.sh $aaset
env oid=1 off=1 $evigenes/prot/aaqual.sh $refaa
mbaopts="-tall -aasize $aaset.qual,$refaa.qual"

$evigenes/makeblastscore3.pl $mbaopts $aablast > $aabltab 

# * ADD here SCRIPT_evgblastp
#   9. namegenes from ref names, for annotation
#     $evigenes/prot/namegenes.pl -blast $aabltab -refnames $refaa.names -out $qname.names
# .. but need ref.names collection (if refgenes.aa but not refgenes.names)

gzip --fast $aablast

echo "#DONE : `date`"
EOS

  my $runapp=$SCRIPT_evgblastp;
  my $runfile= "run_evgblastp.sh";
  $runapp= SCRIPT_common($runapp);
  
  # need data: P_AASET = okayset/*.aa;  P_REFAA == user supplied
  # my $P_DATAD = $settings{runpath} || `pwd`; # FIXME

  my $trname= $runname; # 
  # my($okall,@aas)= getFileset( "okayset",'aa'); 
  my($okmain,$okalt)= getOkFileset( "okayset",'aa$');
  return (STEPerr,'okayset.aa') unless($okmain);

  my $P_REFAA = $settings{REFAA} || finddata('REFAA') || finddata('refgenes.aa'); # default 'refgenes.aa' ?
  return (STEPerr,'REFAA') unless($P_REFAA and -f $P_REFAA);
  $settings{REFAA}= $P_REFAA;
  
  my $P_AASET = $trname."_okall.aa"; # make all? or use trsets/*.tr ?
  unless(-s $P_AASET) {
    my($runerr, $nok, $ofail)= cat_splitset($P_AASET, [$okmain,$okalt]);
    if($runerr) {  } # return(STEPerr)
  }
  
  my($P_BLASTN,$P_NCBIBIN)= findapp('blastp', 1); 
  return(STEPerr,'blastp') if($P_BLASTN =~ /MISSING/);

  # replace script vars here
  map{ 
    s/P_NCBIBIN/$P_NCBIBIN/g;  
    s/P_REFAA/$P_REFAA/g;  s/P_AASET/$P_AASET/g; 
  } ($runapp);
  
  return(STEPok,$runapp,$runfile); 
}

1;

__END__
