#!/usr/bin/env perl
# trasmstrandspan.pl

=item about trasmstrandspan

  trasmstrandspan.pl  -pid=ptf011 [ OR -part=1..n -strand=fwd|rev ]  
  require: -bams=bamdir|bamlist -spantable=strandspan.tab 
  options: -version=v1 -[no]debug -log=log.xxx
  
Assembles short/long rna at partition aligned to genome, using now velvet/oases,
pulling reads from bam files using a table of strandedness spans, 
with long read bams named as 'estsr' or 'estpe'.

Short read bam files are separated by strand as fwd, rev, and no-strand
subsets.  No-strand segments are matched and pulled with 
neighbor fwd / rev segments, to allocate full read sets for stranded assembly.

This is designed for use with paired-end rna-seq that is not collected with 
strand-specific methods, but has been mapped to genome and introns detected that
confer strandedness.  Long reads (EST/454) are also used where stranded (pair end optional).

=item ADD OPTION: trasm for user span

  modify this to process (small) user span input + bamfile list.
  alternately take input as stdin of sam rows, from samtools output stream.
  return trcds.gff (with options)
  rely on common config file, maybe including bamlist

=item strand span table

  -spantable=strandspan.tab has columns of [ chr,start,end,...,fwdid,revid,.. ]; 
  cols 1-3 of partid are fed to 'samtools view -L parttab bam', to pull reads per location.
  option -pid ptf011 means grep out this partid from table, use those locations to pull reads
  from bam for assembly (part=11, strand=fwd == partid ptf011)
  
  Span table maker:  evigene/scripts/rnaseq/pullspanstrandsam.pl  (v1)
      Uses to input bigwig (.bw) files made from strand-segregated reads for
      created location segment tables with spans of strandedness marked.
  Replaced with evigene/scripts/rnaseq/strandspanintab.pl  (v2)
      Uses input introns.gff to produce strand spans table
      

=item data partition table

  #Chr    Beg     End     Nfwd    Nrev    Nno     Idfwd   Idrev   Or
  scaffold00512   799201  800200  164     0       127     ptf001  0       f
  scaffold00512   800201  801000  0       0       0       0       0       0
  scaffold00512   801001  801600  67      0       91      ptf001  0       f
  scaffold00512   801601  802800  0       0       33      0       0       0
  scaffold00512   802801  804900  42      0       536     ptf002  0       f
  scaffold00512   804901  813600  0       0       20      0       0       0
  scaffold00512   813601  813800  34      0       58      ptf002  0       f
  scaffold00512   813801  814100  123     48      90      ptf002  ptr002  b
  scaffold00512   814101  814300  89      12      61      ptf002  ptr002  fr

=item data processing

  1. pullreads per id in spantable and bam list, to new data-id subdir
     creating fixed name files:
        $datas/rna.sam    rna2.sam  rna3.sam 
            file numbers match velvet insert size options for -shortPaired{,2,3..}
        $datas/estpe.sam  estsr.sam  for -longPaired and -long reads          

  2. use velvet/oases multi-kmers to assemble transcripts
      (with some fault tolerance to failed kmer runs)
  3. collate velvkmers/transcripts.fa to trs.fa collection  
  4. map trs.fa to trs.gff on genome with gmap
  5. find orfs and annotate/score trs.gff, using trs.fa, genome.fa, and introns.gff when available.
  6. select best of cdna transcripts, scored for CDS size, CDS/UTR, introns matched, gmapping coverage/ident.
  
=item requires
  
  samtools (now in config); 
    "samtools view -L location.bed bam  region:1-9999" must work
  velveth, velvetg, oases for transcript assembly from reads
  gmap, to align transcripts.tr to genome

  evigene/scripts/ ..
    $evigene/scripts/gmap_to_gff.pl
    $evigene/scripts/genefindcds.pl  
    BioPerl Bio/DB/Fasta.pm used by genefindcds; now in $evigene/lib/Bio/ 
       genefindcds.pl: use lib ("$FindBin::Bin/../lib/"); # $Bin for genefindcds is evigene/scripts/ ; ../lib is right

  genomedir/
    genome.fasta
    gmapdb/
    
  introns.gff for validation of mapped transcripts
    - pulled from aligned reads.bam
    
  prepare input.stranded.chr.bam (see below eg)
      evigene/scripts/rnaseq/bam2stranded.sh
      evigene/scripts/rnaseq/bam2strandimp.sh
      
=item tests

  sbams=grNDtCOcX.fwd.chr.bam,grNDtCOcX.noo.chr.bam,daphmagna09_estsr.bam,daphpulex05_estpe.bam
  $evigene/scripts/rnaseq/trasmstrandspan.pl -pid=ptf002 -spantab=grNDtCOcX.strandspan.tab -bams="$sbams" -debug=p

  $evigene/scripts/rnaseq/trasmstrandspan.pl -debug=p -pid ptf002 -span grNDtCOcX.strandspan.tab \
    grNDtCOcX.{fwd,noo}.chr.bam daph*_est{sr,pe}.bam

=cut

use strict;
use Getopt::Long;

my $DOSPLITJOINS= 0; # needs testing
my $MINLEN= 180; # what?

my $samtools="samtools"; # on path or specify path

# .. CONFIG for velvet
my $velbin="/bio/bio-grid/mb/bin/"; # can be blank for PATH

# FIXME: k89 has highest proportion of quality models (50%); k29 lowest (10%); add k97? for 101bp reads; add k79?
my $kset="89 69 49 41 33 29 23"; # orig runs
#config: kmers= 97 89 79 69 57 47 37 29 # try this set
#        optoases= -edgeFractionCutoff 0.3 -paired_cutoff 0.3 -min_trans_lgth 180 -ins_length 475 -ins_length_long 800

#? maybe split into 2 kset parts; w/ diff options: 
#  hik= fewer reads minpair=4; -edge 0.10 ; lowk = -min_pair_count 9, -edgeFractionCutoff 0.2
# run hik 1st, remove bestmodel reads, then lowk on smaller data set?
my $vopts="-ins_length 200  -ins_length_long 500";
my $oopts="-edgeFractionCutoff 0.10 -min_trans_lgth $MINLEN -ins_length 200 -ins_length_long 500";
# my $vopts="-ins_length 475 -ins_length2 380 -ins_length_long 800";
# my $oopts="-edgeFractionCutoff 0.10 -min_trans_lgth $MINLEN -ins_length 475  -ins_length2 380 -ins_length_long 800";
my $pairinsize= ""; #  "1=grHS,2=grND"; # for -ins_length{,2,3..} and for data input -shortPaired{,2,3..}
  ## pairinsize=grND:380,grHS:475

# .. CONFIG for trmap : fixme ..
my $gmap="/bio/bio-grid/mb/gmap/bin/gmap";
my $evigene="/bio/bio-grid/mb/evigene";
my $genomed="/bio/bio-grid/daphmag/genome/";
my $dgenome="dmagna20100422assembly";
my $gmapdb="$genomed/gmap";
my $introns="intronall12_good.gff"; ##  # FIXME ... caller Path Changes **

# .. CONFIG for trfilterjunk; fixme these score weights need some test, adjust to avoid poor choices
# .. choosing by UTR/%CDS still a hassle; need also cov[erage] to avoid chimera when other alts as good.
## my $bestscores= "inqual:3,cdsindel:-5,pid:2,path:-5,CDS:1,UTR:2";
my $bestscores= "inqual:5,cov:8,pid:5,CDS:2,UTR:1,cdsindel:-2";


# partition table: has partiion ids, scaffold spans w/ fwd/rev reads
#  grNDtCOcX.strandspan.tab

# input bams:
# grNDtCOcX.fwd.chr.bam grNDtCOcX.rev.chr.bam grNDtCOcX.noo.chr.bam
# daphmagna_est09.bam  daphpulex_est05pe.bam
# .. pull part data: samtools view -L part.tab bam | sort > asm.sam

my $optin= join " ", @ARGV;
my %args=();
my $optok= &GetOptions( \%args,
  "spantable|input=s" , "bams=s",
  "pid=s", "i|part=i" , "strand=s",  
  # maybe add istart, istep, imax options to iterate over several part ids for cluster call.
  # trasmtrspan -istart=$i -istep=$ncpu -imax=$lastpart 
  "allstrandspan!", "splitjoins!",
  "samfilter=s", "grepreads=s",  ##  options of read filtering
  "config=s", "output=s" , # outputdir for all files created 
  "logfile=s", "version=s", "debug:s", "nodebug" );


readconfig($args{config}) if($args{config}); # moved 1st, adds to %args now...

my $datatab = $args{spantable} || shift(@ARGV);  ## data table .. now can be (single?) span, -span=scaffold7:110000-210000
my $bams    = $args{bams} || "";
my @morebams= grep /\.bam$/, @ARGV;
if(@morebams) { $bams.="," if($bams); $bams .= join(",",@morebams); }
elsif(@ARGV) { $optok=0; } # no remaining args? OR append to $bams list

my $debug= $args{debug} || 1;
$debug=0 if($args{nodebug});

$DOSPLITJOINS= $args{splitjoins} || 0; # needs testing

my $partid= $args{pid} || 0; 
my $ipart= $args{i} || 0;
my $strand= $args{strand} || ""; #? allow n nostrand? or b/both
my $allstrandspan= $args{allstrandspan} || 0;   # test1 says this is better than using strandspans for noor data..
my $dv= $args{version} || "r".$ipart;  # not: $partid << duplicates subset.

=item grepgreads option

  $evigene/scripts/rnaseq/trasmstrandspan.pl -out tes3hn7e -strand fwd -span scaffold00512:1427411-1477411
   -grepread "| egrep 'NM:i:0|NM:i:1      '" \
   -config=trflam.config -ver n9 -pid=ptf002 gr{HS,ND}tCOcX.{fwd,noo}*.chr.bam > & tes3hn7e.log

  alt: -grepread mismatch:1
  
=cut 

my $samfilter= $args{samfilter} || ""; 
  if($samfilter) { unless($samfilter =~ m/\-[fF] \d/) { warn "#BAD samfilter=$samfilter\n"; $samfilter="";} }

# add to config; allow same config & args, do these after readconfig()
my $grepreads= $args{grepreads} || ""; 
  if($grepreads) { 
    if($grepreads =~ /^mismat\w*\W+(\d+)(.*)/) {
      my($mis,$gx)=($1,$2); # gx useless?
      if($mis>0) { 
        my $nm="0"; for(my $j=1; $j<=$mis; $j++) { $nm.= '|'.$j; }
        $grepreads="| egrep 'NM:i:($nm)\t'"; 
      } else {
        $grepreads="| grep NM:i:0";
      }
    }
    unless($grepreads =~ /grep/) { $grepreads='grep '. $grepreads; }
    unless($grepreads =~ /\|/) { $grepreads='| '. $grepreads; }
  }
  
  
my $outputdir= $args{output} || "";
# my $procid=$$;  # process id for temp file names
my $LOGFILE= $args{logfile} || "log.trasm$$.$dv$partid";  

##above## readconfig($args{config}) if($args{config});

## FIXME: allow input of both strands (mixed)
## strandc should be '' for both/nostrand, but need check below; pass blank only if $strand =~ /both|^b|none|noo|^n/
sub strandc { my $s=shift; return ($s =~ /fwd|for|^f|\+/) ? '+' : ($s =~ /rev|^r|\-/) ? '-' : ''; }
sub strandl { my $s=shift; return ($s =~ /fwd|for|^f|\+/) ? 'f' : ($s =~ /rev|^r|\-/) ? 'r' : 'n'; }
sub strandsuf { my $s=shift; return ($s =~ /fwd|for|^f|\+/) ? 'fwd' : ($s =~ /rev|^r|\-/) ? 'rev' : 'noo'; }


#......... Maybe add for cluster computes
##   # trasmtrspan -istart=$i -istep=$ncpu -imax=$lastpart 
## if($istep and $imax) # loop here for partid,
##   for( my $ipart=$istart; $ipart<=$imax; $ipart += $istep )
##      $partid= sprintf "pt%s%03d",$or,$ipart;
##      MAIN($partid);
##   .. next ipart
#..............



#need strand for this below:  
if($partid and not $strand) {
  if($partid =~ m/([fr]|fwd|for|rev)\d/) { $strand=$1; }  
}

#DEBUG: pullreads( ptn010, kfish2_strandspan5m.tab, , allfungr1-kfish2.bam, )
#CMD: grep ptn010 kfish2_strandspan5m.tab > datak3ptn010/ptn010.sam.tab
#WARN: empty pullreads()

if($ipart and $strand and not $partid) {
  my $or = strandl($strand);
  $partid= sprintf "pt%s%03d",$or,$ipart;
}

my $strandc= strandc($strand); 
## strandc should be '' for both/nostrand, but need check below; pass blank only if 
my $strandok= ($strandc =~ /[+-]/ or ($strandc eq '' and $strand =~ /both|^b|none|noo|^n/))?1:0;
$allstrandspan=1 if($strandc eq '' and $strand =~ /both|^b|none|noo|^n/);

## bug here for strand=both; strandl > partid = ptn<< but strandspan*tab has ptf<< and ptr<< 
## .. need new kfish2_strandspan5m.tab w/ new ptn ids
#DEBUG: pullreads( ptn010, kfish2_strandspan5m.tab, , allfungr1-kfish2.bam, )
#CMD: grep ptn010 kfish2_strandspan5m.tab > datak3ptn010/ptn010.sam.tab
#WARN: empty pullreads()

die "called: $optin
usage: trasmstrandspan.pl  -pid=ptf011 [ OR -part=1..n -strand=fwd|rev ]  
  require: -bams=bamdir|bamlist -spantable=strandspan.tab 
  options: -version=v1 -[no]debug -log=log.xxx  -allstrandspan
  
Assembles short/long rna at partition aligned to genome, using now velvet/oases,
pulling reads from bam files, with long read bams named as 'estsr' or 'estpe',
short read bam files are separated by strand as fwd, rev, and noo subsets.
strandspan.table has columns of [ chr,start,end,...,fwdid,revid,.. ]; 
  cols 1-3 of partid are fed to 'samtools view -L parttab bam', to pull reads per location.
  grep partid is used to pull partition subset (part=11, strand=fwd == partid ptf011)
Span table maker:  evigene/scripts/rnaseq/strandspanintab.pl (v2)
" unless($optok and $datatab and $bams and $partid and $strandok);

# need to trap stderr from velvet, gmap, .. in runcmd() ?
system("touch $LOGFILE"); #? always
my @cleanup;

## ** FIXME: some 1-exon genes have no strand; look for/mark ptnoorID columns in strand.table ?

## ** FIXME: Problem of non-overlapped partitions, solve maybe by adding a few if id+1, id-1 rows to pullreads,
## .. then remove tr outside of span traasm in trmap()  ??

## * OPTION: do only 1 step below, eg redo trmap(), need input of trfile
##     .. pullreads() and trasm_velvet return params if already done.

## FIXME: save partspans in result files for checking, along w/ bamlist? ; end of trasm.gff now

## add version $dv to file/dir names; partid not enough: pullreads?
sub MAIN {}

## FIXME: check in $outputdir for existing data
  warn "#DEBUG: pullreads( $partid, $datatab, $strandc, $bams, )\n" if($debug); ## to LOGFILE ?
my($datadir,$subset,$partspans,$bamlist)= pullreads( $partid, $datatab, $strandc, $bams, ); #... this creates datapart.sams for velasm.
  warn  "#DEBUG: done: ($datadir,$subset,$partspans,$bamlist) = pullreads( $partid, $datatab, $strandc, $bams, )\n" if($debug);

# exit if no datadir/subset : handle last partid outofrange
unless(-d $datadir or  -d "$outputdir/$datadir") { 
 warn "#FAIL: No such ($datadir,$subset) = pullreads( $partid, $datatab, $strandc, $bams, )\n";
 exit -1;
}

if($outputdir) {  #  .. can do rest in subdir ok?
  # problem: makes 2 LOGFILEs, one outside, other inside... unless fullpath
  warn "#DEBUG: chdir($outputdir)\n" if($debug);
  mkdir $outputdir unless(-d  $outputdir); 
  system("mv $datadir $outputdir/") if (-d $datadir); 
  system("mv $LOGFILE $outputdir/") unless($LOGFILE =~ m,/,);
  chdir($outputdir);
}

  warn "#DEBUG: trasm_velvet( $dv, $subset, $datadir, )\n" if($debug);
my @veldirs= trasm_velvet( $dv, $subset, $datadir ); #..

  warn "#DEBUG: trcollate_velvet( $subset, @veldirs)\n" if($debug);
my($trfile, $ntr)= trcollate_velvet( $subset, @veldirs);

if($DOSPLITJOINS) {
  warn "#DEBUG: trsplitjoins( $trfile, $strandc, )\n" if($debug);
  ($trfile)= trsplitjoins( $trfile, $strandc, ); # test, optional?
}

## FIXME: separate out ($trcdsgff)= trfindcds($trgff,$strandc) from trmap()
  warn "#DEBUG: trmap( $trfile, $strandc, )\n" if($debug);
my($trgff)= trmap( $trfile, $strandc, ); #...

## FIXME need to filter both strands in overbestgenes
  warn "#DEBUG: trfilterjunk( $trgff, $trfile, $strandc, )\n" if($debug);
my($trgffclean)= trfilterjunk( $trgff, $trfile, $strandc, ); #...

  warn "#DEBUG: cleanup()\n" if($debug);
cleanup();


#.....................

sub runcmd_sub
{
  my ($dieOrNot,@cmd)= @_;  # options?
  my $cmd= join(" ",@cmd);
  warn "#CMD: ",$cmd,"\n" if($debug);

  # if(@cmd > 1) {  push @cmd, ">> $LOGFILE";  }  # bother; what shell syntax for stderr redir?; use perl STDERR capture?
  # else {  $cmd[0].= " >>$LOGFILE"; }
  # $cmd.=" 1>>$LOGFILE 2>>$LOGFILE";  # NO; still not working dangit.
  #.. sigh .. this works
  open OLDOUT, ">&", \*STDOUT;
  open OLDERR, ">&", \*STDERR;
  open STDOUT, ">>$LOGFILE" ;
  open STDERR, ">>$LOGFILE" ;
  warn "\n#CMD: ",$cmd,"\n";

  my $err= system(@cmd);
  
  open STDOUT, ">&OLDOUT";
  open STDERR, ">&OLDERR";
  
  # FIXME here: dont die here when velvet pukes out on 1 kmer.. put this into LOGFILE ?
  if($err) { 
    warn "#CMD: ",$cmd,"\n" unless($debug);  
    warn "#FAIL: $cmd[0]: $?";
    die if($dieOrNot);
    }
  return $err;
}

sub runCmdOrDie { return runcmd_sub(1, @_); }
sub runcmd { return runcmd_sub(0, @_); }

sub cleanup
{
  # write this as separate cleanup.sh script til done debugging..
  if(1 or $debug) { # always for now
    my $clf= "cleanup.".$LOGFILE.".sh"; $clf =~ s/.log//;
    warn "#FINISH: source $clf\n" if($debug);
    open(O,">$clf");
    foreach my $cup (@cleanup) {
      if( -f $cup) { print O "/bin/rm $cup\n"; }
      elsif( -d $cup) { print O "/bin/rm $cup/*; rmdir $cup\n"; }
    } close(O);
  } else {
    foreach my $cup (@cleanup) {
      if( -f $cup) { unlink($cup); }
      elsif( -d $cup) { system("/bin/rm $cup/*"); rmdir($cup); }
    }
  }  
}

sub readconfig {
  my($cf)= @_;
  open(F, $cf) or die "Missing config: $cf";
  # FIXME: also add these to global %args (dont replace)
  while(<F>) { 
    next unless(/^\w/); chomp;
    s/#.*$//;  s/\s*$//;   # s/=/ /;
    my($k,$v)= split /[=\s]+/,$_,2;
    $k=lc($k); 
    $args{$k}= $v unless($args{$k});
    CFG: for ($k) {
      /^velbin$/  && do { $velbin=$v; last CFG; };
      /^kmers/  && do { $kset=$v; last CFG; };
      /^optvelvet/ && do { $vopts=$v; last CFG; };
      /^optoases/ && do { $oopts=$v; last CFG; };
      /^gmap$/    && do { $gmap=$v; last CFG; };   # was bad for gmapdb xxx
      /^samtools$/ && do { $samtools=$v; last CFG; };  
      /^evigene$/ && do { $evigene=$v; last CFG; };
      /^genomed$/ && do { $genomed=$v; last CFG; };
      /^dgenome$/ && do { $dgenome=$v; last CFG; };
      /^gmapdb$/  && do { $gmapdb=$v; last CFG; };
      /^intron/   && do { $introns=$v; last CFG; };
      /^log$/     && do { $LOGFILE=$v; last CFG; };
      /^minlen$/  && do { $MINLEN=$v; last CFG; };
      /^bestscore/  && do { $bestscores=$v; last CFG; };
      ## more configs: input bam name patt >> velvet insertsize config
      ## pairinsize=1=grHS,2=grND   ... grND:380,grHS:475
      /^pairinsize/ && do { $pairinsize=$v; last CFG; };
      
    }
  }
}

=item bam2strand.sh 

  # split pairs.bam into stranded fwd,rev,noo unstranded parts, including mates w/ stranded reads
  # Note: this now only uses proper pairs -f 0x2, see impaired below.

  #! /bin/bash
  ### env  subd=bamq5treat qsub -q batch bam2strand.sh
  #PBS -N bam2strand
  #PBS -l nodes=1:ppn=24,walltime=22:55:00
  #PBS -o bam2strand.$$.out
  #PBS -e bam2strand.$$.err
  #PBS -V
  
  ncpu=24
  subd=bamq5treat
  datad=$HOME/scratch
  workd=$datad/chrs/daphmag/
  rund=$workd/rnas/$subd
  dgenomesize=$workd/genome/dmagna20100422assembly.fa.fai
  
  sdir=$HOME/bio/evigene/scripts/rnaseq
  # samtools=$HOME/bio/bin/samtools
  module add samtools
  
  cd $rund
  # first, only tCO : control
  treats=`ls -d {HS,ND}tCOc[XI]`
  bams=`ls {HS,ND}tCOc[XI]/*dmag2.bam`
  
  echo "start bam2strand.$subd : `date`"  
  # steps 1,2,3  foreach bam in bams; 
  #    samtools sort -n chr.bam > names.bam; 
  #    samtools view -f 0x2 names.bam | perl tostrands : (fwd,rev,noo).sam
  #    samtools view -1 -t xxx.size -o part.bam  part.sam : (fwd,rev,noo)

  #* Modify here for impaired: pull improperly paired mates in gaps;
  #  step 2i: samtools view -F 0x2 names.bam | perl toimpairedstrands.pl ..
     
  cat <<-'EOP' > tostrands.pl
   BEGIN{ $f=$ENV{fn} or die; open(F,">$f.fwd.sam"); open(R,">$f.rev.sam"); open(N,">$f.noo.sam");} 
   while(<>) { ($d)=split; ($o)= m/\tXS:A:(.)/; if($d eq $ld) { if($o or $lo) { $oc=$lo||$o; 
    if($oc eq "+") { print F $ll,$_; } elsif($oc eq "-") { print R $ll,$_; } } else { print N $ll,$_; } }
    ($ll,$ld,$lo)= ($_,$d,$o); }
   END{ close(F); close(R); close(N); }
  EOP
  
cat <<-'EOX' > toimpairstrands.pl
BEGIN{ $f=$ENV{fn} or die; open(F,">$f.fwdimp.sam"); open(R,">$f.revimp.sam"); open(N,">$f.nooimp.sam");} 
while(<>) {
  my($d,@v)=split"\t"; next if(($v[0] & 12) == 12); ($o)= m/\tXS:A:(.)/;
  if($d eq $ld) { 
    if(bothmap($lv[0],$v[0])) { putimp( $lo||$o, $ld,\@lv,$d,\@v); }  
    if(swapmateloc(\@lv,\@v)) { putimp( $lo||$o, $ld,\@lv,$d,\@v); }  
  } 
  ($ld,$lo,@lv)=($d,$o,@v); 
}
sub bothmap { my($lf,$vf)=@_; return ($lf & 4 or $vf & 4)?0:1; } 
sub putimp { my($oc,$ld,$lv,$d,$v)=@_; my $ll=join"\t", $ld,@$lv; my $vv= join"\t", $d,@$v;
  unless($oc) { print N $ll,$vv; } elsif($oc eq "+") { print F $ll,$vv; } elsif($oc eq "-") { print R $ll,$vv; } 
}
sub swapmateloc { my($l,$v)=@_; for my $r ($l,$v) { my($c,$cb,$mc,$mb)=@{$r}[1,2,5,6]; 
  if($mc eq "*") {} elsif($mc eq "=") { @{$r}[1,2,5,6]=($c,$mb,$mc,$cb); } 
  else { @{$r}[1,2,5,6]=($mc,$mb,$c,$cb); } 
  } return 1;} 
EOX

#................

  #2i: samtools view -F 0x2 $nam.names.bam | env fn=$nam perl toimpairstrands.pl 
      
  i=0; 
  for bam in $bams; do {
    nam=`echo $bam | sed 's/.bam//; s/-dmag2//;'`
    ( samtools sort -n  $bam  $nam.names; 
      samtools view -f 0x2 $nam.names.bam | env fn=$nam perl tostrands.pl ) &  
    i=$(( $i + 1 ))
    if [ $i -ge $ncpu ]; then wait; i=0; fi
  } done
  wait
  
  samstr=`ls {HS,ND}tCOc[XI]/*.{fwd,rev,noo}.sam`
  #2i: samstr=`ls {HS,ND}tCOc[XI]/*.{fwd,rev,noo}imp.sam`
  i=0; 
  for sam in $samstr; do {
    nam=`echo $sam | sed 's/.sam//;'`
    samtools view -1 -t $dgenomesize -o $nam.bam $sam &
    i=$(( $i + 1 ))
    if [ $i -ge $ncpu ]; then wait; i=0; fi
  } done
  wait
  
  # steps 4,5  .. merge parts into groups, location sorted
  # 4. samtools merge -1 -n groupj.fwd.bam  part1.fwd.bam part2.fwd.bam
  # 5. samtools sort groupj.fwd.bam groupj.chrs.fwd
  
  i=0; 
  for tdir in $treats; do {
    nam=`basename $tdir | sed 's/^/gr/;'`
    #2i: for ord in fwdimp revimp nooimp; do
    for ord in fwd rev noo; do
    {
      parts=`ls $tdir/*.$ord.bam`
      groupf=$nam.$ord
      ( samtools merge -1 -n $groupf.bam  $parts;
        samtools sort  $groupf.bam  $groupf.chr ) &
      i=$(( $i + 1 ))
      if [ $i -ge $ncpu ]; then wait; i=0; fi
    } done  
  } done
  wait
  
  # then bam2bed.sh for each groupf.chr.bam
  echo "end bam2strand: `date` "

=cut 


=item impaired: pull mates in gaps

 FIXME: pullreads BUG, input reads now from sam -f 0x2 proper pairs are missing mates that fall in gaps ** 
 i.e. mapped on small scafs/contigs that should fill gaps
 for trasmstrandspan.pl pullreads, need way to pull reads + improperpair mates from location select.
  ..  need way other than location span to pull unmapped/othermapped mates.
      small subset of current dmag bamq5 data, 5% of reads, 5Mill/100Mil (2.5% unmapped, 3% improper pairing)
      but pull by read-id is tedious; create mateloc bam of impaired reads? 
          samtools view -f 0x2 $bam $span >> rna.samu   # proper pairs
          samtools view $bam.impair $span >> rna.samu # improp pairs; locations here are of pairs in bam
          sort rna.samu > rna.sam
          
      make bam.impaired:
        # hex: 8, 9, 10=A, 11=B, 12=C, 13=D, 14=E, 15=F
        # -F 0x2 = skip proper pairs; flag 12 = 0x4+0x8 = ? 0xC = both unmapped
        # ?? need print both orig locs and swapped ?
        
    nam=`echo $bam | sed 's/.bam//;'`
    samtools view -F 0x2 $bam | sort | perl -ne \
    'my($d,@v)=split"\t"; next if(($v[0] & 12) == 12); \
    if($d eq $ld) { if(bothmap($lv[0],$v[0])) { print join"\t",$ld,@lv; print join"\t",$d,@v; } \
    if(swapmateloc(\@lv,\@v)) { print join"\t",$ld,@lv; print join"\t",$d,@v; } } ($ld,@lv)=($d,@v); \
    sub bothmap { my($lf,$vf)=@_; return ($lf & 4 or $vf & 4)?0:1; } \
    sub swapmateloc { my($l,$v)=@_; for my $r ($l,$v) { my($c,$cb,$mc,$mb)=@{$r}[1,2,5,6]; \
    if($mc eq "*") { } elsif($mc eq "=") { @{$r}[1,2,5,6]=($c,$mb,$mc,$cb); } \
    else { @{$r}[1,2,5,6]=($mc,$mb,$c,$cb); } \
    } return 1;} '\
    > $nam.impair.sam
    
    samtools view -u -t $chrsize $nam.impair.sam | samtools sort - $nam.impair
      # side effect of^ dupl locs, pairs are name-sorted as well as location sorted.


bam f2-9:
p1      scaffold00512   106160  40      6S95M   scaffold02494   2341    0
p1      scaffold00512   106160  40      8S93M   scaffold02494   2352    0
p1      scaffold00512   106160  40      1S100M  scaffold02494   2277    0
p1      scaffold00512   106160  40      8S93M   scaffold02494   2315    0
         
 eg. scaffold00512:105000-110600  gap in myosin heavy chain,  imp mates are on scaffold02494
   .. could fill gaps from rna mates of improp pairs, analyzing mate map to small scaf, and
      insert direction effects, .. all -inserts below gap are proper pairs, all +inserts missing
      for imp pair from scaffold00512 to scaffold02494; vice versa above gap.



=item gapfill tables

  input assembly gap.gff locations, paired reads
  output table of improper paired reads near gaps with a consensus location of mates, eg other scaffold or unmapped.
  
  set bam=Dman_62.bam ; 
  set bam=Dman_76.bam ; grep scaffold00512 dmagna20100422assembly.gaps.gff | env bam=$bam perl -ne \
  'my($r,$s,$t,$bb,$e,$w)=split; $bl=$bb-$GAPOFS; $er=$e+$GAPOFS; print join("\t",$r,$bb,$e,"gap",$w)."\n"; \
  pullg($r,$bl,$bb,"gapl"); pullg($r,$e,$er,"gapr"); \
  sub pullg{ my($r,$rb,$re,$tag)=@_; open(F,"samtools view -X $bam $r:$rb-$re | cut -f2-9 |"); my(%mc,%mb,$mn); \
  while(<F>){  next if(/\t=\t/); my($fl,$c,$cb,$mq,$cg,$mc,$mb,$in)=split"\t"; $mc{$mc}++; $mb{$mc}.="$mb,"; $mn++;} \
  my @mc= sort{ $mc{$b}<=>$mc{$a} } keys %mc; my($n1,$n2)= map{ $mc{$_} } @mc; if($n1>9 and $n1 > $n2*2) {\
  my $mc=$mc[0]; my @mb=sort{$a<=>$b} split",",$mb{$mc}; ($mb,$me)=@mb[0,-1]; \
  print join("\t",$r,$rb,$re,$tag,"$n1/$mn",$mc,"$mb-$me")."\n"; } } \
  BEGIN{ $bam=$ENV{bam}||"Dman_62.bam"; $GAPOFS=$ENV{ofs}||1000; } ' > gapfill.sc512.dman76.tab

=cut


sub pullreads
{
  my( $partid, $datatab, $strandc, $bams, )= @_;

  ##opt# my $allstrandspan= 0; # option to use all noo reads for stranded span
  ##FIXed: new opt, allow datatab to be span: "$ptchr:$ptb-$pte "
  ##FIXME: new opt, filter reads by criteria like NM:i:0/no mismatches.. no cigar =~ \dS softclips, ..
  
  ## my $partid= sprintf "pt%s%03d",$strand,$ipart;
  ## add version $dv to file/dir names; partid not enough: pullreads?
  my $subset= $partid; #??
  my $ptspans=""; # FIXME: save ptspans in result files for checking, along w/ bamlist?.
  my $cmd;

  my($datadir) = "data$dv$partid";  ## version here dv
  ## return if exist datadir/rna.sam ..  ?
  ## FIXME: check in $outputdir for existing data

  if( -d $datadir and -s "$datadir/rna.sam" ) {
    warn "#NOTE: reusing pullreads()\n";
    return($datadir,$subset,$ptspans); # return nonzero files.sam ?
  } 
  elsif( -d $outputdir and -d "$outputdir/$datadir" and -s "$outputdir/$datadir/rna.sam" ) {
    warn "#NOTE: reusing pullreads()\n";
    return($datadir,$subset,$ptspans); # return nonzero files.sam ?
  }
  
  mkdir($datadir); push @cleanup, $datadir;

use constant NEWTABREAD => 1;  
    ## change this, drop grep, open & read file; assume datatab is location sorted
  my $ptab="$datadir/$partid.sam.tab";
  my($ptchr,$ptb,$pte,$ngot,$callerspan)= (0) x 10; # may be many ptchr.. do what?
  
  if($datatab =~ /^(\w+):(\d+)[\.\-]*(\d+)/) { #is span, allow 2+?  chr1:1-9,chr1:19-29,.. ?
    ($ptchr,$ptb,$pte)= ($1,$2,$3);  $callerspan=1;
    $ptspans.= "$ptchr:$ptb-$pte ";
    open(O,">$ptab") or die "writing $ptab";
    # $partid="xxx1" unless($partid);
    print O join("\t",$ptchr,$ptb,$pte,$partid)."\n"; $ngot++;
    close(O);
   
  } elsif(  -f $datatab) {
  
  # FIXME: ptf004 includes scaffold00024 after scaffold00512; failed to grep to ptf004.sam.tab
  # if(NEWTABREAD) {
  warn "#CMD: grep $partid $datatab > $ptab\n" if($debug);
  open(T,$datatab) or die "no $datatab"; 
  open(O,">$ptab") or die "writing $ptab";
  while(<T>) { 
    next unless(/^\w/); 
    my @v=split; 
    my($c,$b,$e)= @v[0,1,2];

    if(/\t$partid\t/) {  # FIXME: both == ptn000, strandids == ptf000,ptr000
      print O $_; $ngot++;
      
      if($ptchr and $c ne $ptchr) { # new chr span
        $ptspans.= "$ptchr:$ptb-$pte ";
        $ptb=$pte=0; $ptchr=$c;
      } else { 
        $ptchr=$c; 
      }
      $ptb= $b if($ptb==0 or $ptb>$b);
      $pte= $e if($pte<$e);
      # keep a few rows of last partition?
      
    } else {
      # keep a few of next partition? need to parse partition ids ..
      # last if($ngot); # NO, need to know if it is last partid .. last if($ngot and $c ne $ptchr) ?
      # .. need more info to decide to stop..
      #BAD# last if($ngot and $c ne $ptchr);  #ptf004 scaffold00024 FAILs here..
    }

  } 
  close(O); close(T);
  $ptspans.= "$ptchr:$ptb-$pte ";  # FIXME: preserve ptspans for gff filtering : separate inspan/outspan trmap?
}

# } else {  
#   $cmd="grep $partid $datatab > $ptab";  # add "cut -f 1-3" for samtools?
#   runCmdOrDie($cmd); 
#   $ngot=`wc -l $ptab`; $ngot=int($ngot);
# }

  if( $ngot < 1 ) {
    warn "#WARN: empty pullreads()\n";
    return("",0,0); #   return nonzero files.sam ?
  }

  my @bams;
  if( -d $bams) {
    @bams= `ls -1 $bams/*.bam`;
  } else {
    @bams= grep /\.bam/, split /[,\s]+/, $bams;
  }
  
  my $locustab= ($allstrandspan||$callerspan) ? "" : "-L $ptab";   # allstrandspan works better .. more noor data
    ## FIXed: allow both-strands w/o fwd,rev name when strandc eq '', ok for allstrandspan
  
  my @estpe= grep /estpe/, @bams; # need estpe, estsr names ??
  my @estsr= grep /estsr/, @bams; # need estpe, estsr names ??
  @bams = grep !/est/, @bams; #??

  ## need option for -shortPaired2 $datas/rna2.sam .. other insertsizes
  ## ASSUME vopts, oopts have needed -ins2 -ins3 .. values
  # my $vopts="-ins_length 475 -ins_length2 380 -ins_length_long 800";
  # my $oopts="-edgeFractionCutoff 0.10 -min_trans_lgth $MINLEN -ins_length 475  -ins_length2 380 -ins_length_long 800";
  ## .. handle this in pullreads .. make rna[1234..].sam
  # my $pairinsize= "1=grHS,2=grND"; << this way  "grHS=1,grND=2,grXX=3"
  
  my $havereadins=0;
  my (@rins,%rinpatt);
  if($pairinsize) {
    @rins= split /[, ]+/, $pairinsize;
    if(@rins>1) { 
      map{ my($k,$v)= split /[:=]/, $_; $k= int($k)||1; chomp($v); $rinpatt{$k}=$v; } @rins;
      @rins= sort keys %rinpatt;
      $havereadins= (@rins>1)?1:0;
    }
  }

  system("touch $datadir/rna.sam");
  foreach my $bam (@bams) {
    ## change this pull: for fwd,rev.bam, pull entire span from ptab start-end,
    ## only for noo strandless use the ptab segments near stranded spans.
    
    my $pullto="$datadir/rna.sam";  # for pairinsert config by filename
    if($havereadins) {
      foreach my $rin (@rins) { 
        if($bam =~ m/$rinpatt{$rin}/) { 
          $pullto="$datadir/rna$rin.sam" if($rin>1);
          last;
        }
      }  
    }
      
    # FIXME? add strand filter here now? reverse of est: grep -v 'XS:A:$otherstrandc'; 
    # no good this ignores unstranded mates, set mate strand before ?
    # FIXed: add filter opts for samtools:  -f xx -F yy and attribs: 'NM:i:0' and cigar flags? like \dS softclip
    ##my $sopt = ($samfilter) ? $samfilter:'';
    ##my $grepr= ($grepreads) ? $grepreads:''; ## "| grep 'NM:i:0'";
    
    if($bam =~ /(fwd|rev)/) { $cmd= "$samtools view $samfilter $bam $ptspans "; }
    else { $cmd= "$samtools view $samfilter $locustab $bam $ptspans "; }  # FIXed use both -L tab and region
    $cmd.= $grepreads." | sort >> $pullto";

    runCmdOrDie($cmd);
  }
  
  if(@estsr) {
  system("touch $datadir/estsr.sam");
  foreach my $bam (@estsr) {  
    # allow nostrand, all orient?
    if($strandc) { $cmd="$samtools view $bam $ptspans | grep 'XS:A:$strandc' | sort >> $datadir/estsr.sam"; }
    else { $cmd="$samtools view  $locustab $bam $ptspans | sort >> $datadir/estsr.sam"; }
    runCmdOrDie($cmd);
    }
  }
  if(@estpe) {
  system("touch $datadir/estpe.sam");
  foreach my $bam (@estpe) {  
    if($strandc) { $cmd="$samtools view $bam $ptspans | grep 'XS:A:$strandc' | sort >> $datadir/estpe.sam"; }
    else { $cmd="$samtools view  $locustab $bam $ptspans  | sort >> $datadir/estpe.sam"; }
    runCmdOrDie($cmd);
    }
  }

  my $bamlist= "rna:". join(",",@bams) . ",est:" .join(",",@estpe,@estsr);
  # fixme maybe: collect linecount from datadir/*.sam for outputs?
  unless( -s "$datadir/rna.sam" ) {
    warn "#WARN: empty pullreads()\n";
    return("",0,0); # $datadir,$subset,$ptspans # return nonzero files.sam ?
  }

  return($datadir,$subset,$ptspans,$bamlist); # return nonzero files.sam ?
}



sub trasm_velvet
{
  my($dv, $subset, $datas, )= @_;  
  
  my @subdirs;
  #global now# my $kset="89 69 49 41 33 29"; # config param ... CHANGE see above

  my $kmiss="";
  foreach my $k (split" ",$kset) 
  {
    my $ksubdir="vel$dv${subset}_$k"; ## check for extra ks ?
    if(-d $ksubdir) { push @subdirs, $ksubdir; }
    else { $kmiss.="$k "; }
  }
  if(@subdirs) { 
    if($kmiss) { $kset=$kmiss; } 
    else { warn "#NOTE: reusing trasm_velvet()\n"; return (@subdirs); } 
  }
    
  # ins: Hels=475 Ndame=350
  # longpaired estdpx05: median insert = 800+
  #global# my $vopts="-ins_length 475  -ins_length_long 800";
  #global# my $oopts="-edgeFractionCutoff 0.10 -min_trans_lgth $MINLEN -ins_length 475 -ins_length_long 800";
  # .. maybe increase edge to reduce run-ons/joins
  
  # .. extra oopts: edgecut, paircut maybe/maybe not useful.
  # velg:  -conserveLong yes   oases: -merge yes 
  
  my $kseqdir="vel$dv${subset}_seq";
  my $cmd="";
  my @cmd;
  
  if($velbin and $velbin !~ m,/$,) { $velbin.="/"; }
  
  ## need option for -shortPaired2 $datas/rna2.sam .. other insertsizes
  ## ASSUME vopts, oopts have needed -ins2 -ins3 .. values
  # my $vopts="-ins_length 475 -ins_length2 380 -ins_length_long 800";
  # my $oopts="-edgeFractionCutoff 0.10 -min_trans_lgth $MINLEN -ins_length 475  -ins_length2 380 -ins_length_long 800";
  # my $pairinsize= "grHS=in1,grND=in2"; # for -ins_length{,2,3..} and for data input -shortPaired{,2,3..}
  ## .. handle this in pullreads .. make rna[1234..].sam
  # my @rins= split /[, ]+/, $pairinsize;
  # if(@rins>1) { }
  
  @cmd = ("${velbin}velveth", $kseqdir, 27, "-sam", "-shortPaired",  "$datas/rna.sam");
  # .. iterate i=2..n here ?  
  push @cmd, "-shortPaired2",  "$datas/rna2.sam" if(-s "$datas/rna2.sam");
  push @cmd, "-shortPaired3",  "$datas/rna3.sam" if(-s "$datas/rna3.sam");
  #?unpaired# push @cmd, "-short",  "$datas/rna1.sam" if(-s "$datas/rna1.sam");

  push @cmd, "-long", "$datas/estsr.sam" if( -s "$datas/estsr.sam");
  push @cmd, "-longPaired", "$datas/estpe.sam" if( -s "$datas/estpe.sam");
  push @cmd, "-noHash";
  runCmdOrDie(@cmd);  
  
# velveth: Could not open 2>>log.trasm16024: No such file or directory
# #FAILED: /bio/bio-grid/mb/bin/velveth: 256 at /bio/bio-grid/mb/evigene/scripts/rnaseq/trasmstrandspan.pl line 132.
  my $err; my $kfails;
  my @kset= split" ",$kset;
  while( my $k= shift @kset ) 
  {
    my $ksubdir="vel$dv${subset}_$k";
    mkdir $ksubdir;
    system("ln -s ../$kseqdir/Sequences $ksubdir/"); #?
    
    @cmd= ("${velbin}velveth", $ksubdir, $k, "-reuse_Sequences");
    $err= runcmd(@cmd); 
    if($err and not($kfails =~ m/$k/)) { my $k2=$k+2; push @kset, $k2; $kfails .= "$k,$k2,"; }  
    next if $err;
    
    @cmd= ("${velbin}velvetg", $ksubdir, split(" ",$vopts), "-read_trkg","yes");
    $err= runcmd(@cmd); 
    if($err and not($kfails =~ m/$k/)) { my $k2=$k+2; push @kset, $k2; $kfails .= "$k,$k2,"; }  
    next if $err;

    @cmd= ("${velbin}oases", $ksubdir, split(" ",$oopts));
    $err= runcmd(@cmd); # fails common here
    #.. if oases fails on this kmer, try another k+2 or k-2 ? FIXME: ditto velvetg, velveth?
    if($err and not($kfails =~ m/$k/)) { my $k2=$k+2; push @kset, $k2; $kfails .= "$k,$k2,"; }  
    # but return ksubdir as having data, can use contigs.fa from velvetg
    
    #system("/bin/rm $ksubdir/{Graph2,LastGraph,PreGraph,Roadmaps,Sequences}");
    foreach my $vtf (qw(Graph2 LastGraph PreGraph Roadmaps Sequences)) { unlink "$ksubdir/$vtf"; }
    #? veln9ptf001_27/{Graph2,LastGraph,PreGraph,Roadmaps,Sequences}: No such file or directory << solaris, what err?
    # dont need; later erase full ksubdir; but do anyway?
    
    push @subdirs, $ksubdir; # ok or fail?
  }
  
  # wait(); #?? waitall() ??  NOT HERE.. caller has cpu loop
  push @cleanup, $kseqdir;
  return(@subdirs); # what? list of ok ksubdir?
}



sub trcollate_velvet # veltrmake.sh
{
  my($subset,@subdirs)= @_;

  my $rnam=$subdirs[0];
  $rnam =~ s,^.*/,,;  $rnam =~ s/_[0-9][0-9]/s/; 
  # $rnam =~ s/veldm/daphmag2vel/;
  
  my $ntr=0;
  rename("$rnam.tr","$rnam.tr.old") if(-f "$rnam.tr");
  system("touch $rnam.tr");
  foreach my $sd (@subdirs) {
    my $tr="$sd/transcripts.fa";
    my $ct="$sd/contigs.fa";
    push @cleanup, $sd;
    my $nam=$sd; $nam=~s,^.*/,,; $nam=~s/_/k/; 
    
    if( -s $tr ) { # Fixed: was -f ; need nonzero
      open(F,$tr); open(O,">>$rnam.tr");
      while(<F>) { 
        if(/^>/) { $ntr++;
          s,/(\d+)_Confidence_, nt=$1; cf=,; s/_Length_/; len=/; s/Locus_/Loc/; s/_Transcript_/t/; s/>/>$nam/; 
          }
        print O $_;
      } close(O); close(F);
        
    } elsif( -s $ct ) { # Fixed: was -f ; need nonzero
      open(F,$ct); open(O,">>$rnam.tr"); my $p=0;
      while(<F>) { 
        if(/^>/) { $ntr++; s/_length_(\d+)_cov_(\d+).*/; len=$1; cf=$2;/; my $w=$1; 
          $p=($w>=$MINLEN)?1:0; s/>NODE_(\d+)/>${nam}Loc${1}t0/; }
        print O $_ if($p);
      } close(O); close(F);    
    }
  }
  
  return("$rnam.tr", $ntr);
}

sub gffsuf{ my($f,$s,$d)= @_; my $fn=$f; $fn=~s/\.$d// if($d); unless($fn =~ s/\.gff/.$s.gff/) { $fn.=".$s"; } return $fn;}

# FIXME: insert new step b/n trcollate(), trmap() : generic split joins using cdna_bestorf tests of utrorf
sub trsplitjoins   
{
  my($trfa,$strandc)= @_;
  
  my $trfasplit= $trfa.".split";
  my $cmd= "$evigene/scripts/cdna_bestorf.pl -splitutrorf -act cdnaonlyfasta -cdna $trfa > $trfasplit";
  my $err= runcmd($cmd);
  if($err) { return ($trfa); }
  else { push @cleanup, $trfa; return($trfasplit); }
}

sub trmap  # trgmap.sh
{
  my($trfa,$strandc)= @_;

  my $cmd;
  my $name= $trfa; $name=~ s/\.(tr|fasta|fa)$//; # `basename $trfa .tr`;
  my $trgff= "$name.gff";
  
#  # .. CONFIG for trmap ..
#   my $genomed="/bio/bio-grid/daphmag/genome/";
#   my $dgenome="dmagna20100422assembly";
#   my $introns="intronsc389.gff";
#   my $gmapdb="$genomed/gmap";
#   my $gmap="/bio/bio-grid/mb/gmap/bin/gmap";
#   my $evigene="/bio/bio-grid/mb/evigene";
  ## gmap opt?  -c, --chrsubset=string         Chromosome subset to search

  $cmd="$gmap --npaths=0 --min-intronlength=30 -S -D $gmapdb -d $dgenome $trfa | ".
  "env src=$name noerrspan=1 intron=-1 best=0 $evigene/scripts/gmap_to_gff.pl > $trgff";
  runCmdOrDie($cmd); 
  
  ## FIXME here: allow for unstranded 1-exon genes, use genefindcds prot to set strand
  my $strf= strandsuf($strandc);
  if($strandc) {
    my $err=0;
    my $trgffor= gffsuf($trgff,$strf); # (my $trgffor=$trgff) =~ s/.gff/.$strf.gff/; #  "$trgff.$strf";
if(1) {
    my $p=1; my $nok=0; 
    open(F,"$trgff") or $err++; 
    open(O,">$trgffor") or $err++;
    while(<F>) {
      if(/^\w/ and /\tmRNA/) { 
        my @v=split"\t"; my $or=$v[6]; 
        if($or eq $strandc) { $p=1; }
        elsif($or eq "." and m/nexon=1;/) { 
          my %at=(); 
          while(m/(aalen|cov|pid|cdsindel)=(\d+)/g) { $at{$1}=$2; }
          my $pc= $at{cov} * $at{pid} / 100;
          $p= ($pc>97 and $at{aalen}>=40 and $at{cdsindel} < 1)? 1: 0; # FIXME: config these opts; only for 1exon cases
        } else { $p=0; }
        $nok++ if($p);
      }
      print O $_ if($p);
    }
    
if(1) { # add some docs to end of gff; better put elsewhere; losing this later steps
    print O "\n";
    print O "#n dataspans=$partspans\n";
    print O "#n datainput=$bamlist\n";
}      
    close(O); close(F);
    $err++ unless($nok);
} else {  
    ## dont grep this, read in gff, decide to keep no-strand if : nexon=1, not chimera, aalen>minaa, cov*pid>97%, cdsindex=0
    $cmd= "cat $trgff | grep -v '#' | grep '\t$strandc\t\.' > $name.gff1";
    $err= runcmd($cmd); 
}    
    unless($err) {
    push @cleanup, $trgff; $trgff= $trgffor; 
    # rename( "$trgff", "$name.gff0");  push @cleanup,  "$name.gff0";
    # rename( "$name.gff1", "$trgff");
    }
  }
  
  #? allow for no intron.gff here?
  my $trgffan= gffsuf($trgff,"an","gmap"); # $trgff; $trgffan =~ s/.gff/.an.gff/;
  my $dna= "$genomed/$dgenome.fa";
  unless( -f $dna ) { $dna="$genomed/$dgenome.fasta"; }
  unless( -f $dna ) { die "#FAIL: MISSING genome dna $dna"; }
  
  # FIXME: failed here via Fasta.pm index try on dirty dgenome.fa ..
  $cmd= "$evigene/scripts/genefindcds.pl -ratiocdna 1.25 -genes $trgff -cdna $trfa -dna $dna";
  if($introns) {
    if(-f $introns) { $cmd.= " -nofixintron -intron $introns "; }
    elsif( -f "../$introns") {  $cmd.= " -nofixintron -intron ../$introns "; }
    else { warn "#ERROR: MISSING introns $introns\n"; }
  }
  $cmd.= " > $trgffan";
  runCmdOrDie($cmd); #? nodie; can we bypass failed genefindcds? can run overbestgene1 after failed genefindcds
  #?keep# push @cleanup, $trgff; 
  $trgff= $trgffan;
  
  # LOG result : 
  my $ntrmap= `grep -c mRNA $trgff`; chomp($ntrmap);
  warn "#NOTE: trmap n=$ntrmap mRNA in $trgff\n";
  return( $trgff);  
}

=item fail no BioPerl

#FAIL: /N/u/gilbertd/Mason/bio/evigene/scripts/genefindcds.pl 

/N/u/gilbertd/Mason/bio/evigene/scripts/genefindcds.pl -genes veln4ptr005s.gff -cdna veln4ptr005s.tr \
-dna /N/u/gilbertd/Mason/scratch/chrs/daphmag/genome//dmagna20100422assembly.fa \
# > veln4ptr005s.an.gff: 512 at /N/u/gilbertd/Mason/bio/evigene/scripts/rnaseq/trasmstrandspan.pl line 234.

Can't locate Bio/DB/Fasta.pm in @INC () at /N/u/gilbertd/Mason/bio/evigene/scripts/genefindcds.pl line 1292, <STDIN> line 18.

=cut


sub trfilterjunk
{
  my($trgff, $trfile, $strandc, )= @_;
  ## my $trgffbest= $trgff; $trgffbest =~ s/.gff//; $trgffbest =~ s/\.an//; $trgffbest.=".best1.gff";
  my $trgffbest= gffsuf($trgff,"best1","an"); # drop .fwd.an.gff
  
  my $cmd;
  
## FIXME need to filter both strands in overbestgenes
  
  ## fixed: inqual=$inok/$inerr/$inzip ; need to score -inerr > +inok << redid genefindcds inqual
  ##  .. inqual=nn%,inok/inerr/inzip now
  ## fixme: scale weights for score values: 
  #    inqual=nn%; cdsindel=1..10s also -val; pid=nn%, path=1,2 CDS=bases; UTR=nn%
  # score: inqual:10  from genefindcds. should score aalen=size,pCDS,complete/partial
  #    .. inqual=pct << waight by nintrons, i.e inqual=90,7/0/0 > 90 * 7; better than inqual=100,2/0/0
  # score: cdsindel:-9 from genefindcds << BUGgers; cdsindel neg ; need some scorewt fix: abs-cdsindel:-5 ??
  # ^^ may not be good score ..
  # score: chimera:-1 reduce? use path=n/m as proxy
  # score: cov=, pid= from gmap, coverage, percentidentity; score both?
  # score? utrx= extra utr exons from genefindcds, neg weight 
  
  # ** ADD TO CONFIG **
  # .. choosing by UTR/%CDS still a hassle; need also cov[erage] to avoid chimera when other alts as good.
  ##old# my $bestscores= "inqual:3,cdsindel:-5,pid:2,path:-5,CDS:1,UTR:2";
  ##new# my $bestscores= "inqual:5,cov:8,pid:5,CDS:2,UTR:1,cdsindel:-2";
  my $tscore= $bestscores || "inqual:3,cov:5,pid:2,CDS:1,UTR:1";
  $cmd="$evigene/scripts/overbestgene1.perl -alttr -pct 10 -strand -score '$tscore' -in $trgff > $trgffbest";
  my $err= runcmd($cmd);

  my $ntrmap= `grep -c mRNA $trgffbest`; chomp($ntrmap);
  warn "#NOTE: trbest n=$ntrmap mRNA in $trgffbest\n";
  return($trgffbest);
}


__END__

=item trasmspan.config for dgg

MINLEN=180 
samtools=samtools 
# .. CONFIG for velvet
velbin=/bio/bio-grid/mb/bin/ 
# later: 
# kset=89 69 49 41 33 29 # param...
# vopts=-ins_length 475  -ins_length_long 800
# oopts=-edgeFractionCutoff 0.10 -min_trans_lgth $MINLEN -ins_length 475 -ins_length_long 800

# .. CONFIG for trmap : fixme ..
gmap=/bio/bio-grid/mb/gmap/bin/gmap
evigene=/bio/bio-grid/mb/evigene
genomed=/bio/bio-grid/daphmag/genome/
gmapdb=/bio/bio-grid/daphmag/genome/gmap
dgenome=dmagna20100422assembly
introns=intronall12_good.gff  

=item trasmspan.config for mason.iu

MINLEN=180 
samtools=samtools 
velbin=/N/u/gilbertd/Mason/bio/velvet127s/bin2/ 
gmap=/N/u/gilbertd/Mason/bio/gmap/bin/gmap
evigene=/N/u/gilbertd/Mason/bio/evigene
genomed=/N/u/gilbertd/Mason/scratch/chrs/daphmag/genome/
gmapdb=/N/u/gilbertd/Mason/scratch/chrs/daphmag/genome/gmap
dgenome=dmagna20100422assembly
introns=intronall12_good.gff  

=item trasm cluster

#! /bin/bash
### env  spantab=xxx bamlist=xxx strand=fwd|rev qsub -q batch trasmspan.sh
#PBS -N trasmspan
#PBS -l nodes=1:ppn=16,walltime=22:55:00
#PBS -o trasmspan.$$.out
#PBS -e trasmspan.$$.err
#PBS -V

ncpu=16
dv=n4

strand=fwd
spantable=grNDtCOcX.strandspan.tab
bamlist=`ls grNDt*cX.{$strand,noo}.chr.bam daph*_est{sr,pe}.bam`

firstpart=1
lastpart=32
# lastpart=115
## fwd max=108 rev max=112 for grNDtCOcX.strandspan.tab
subd=bamq5treat/strandgroups
config=trasmspan.config

## example runs: 
## replace -pid=xxx$i with -strand=fwd -part=$i
## fwd strand
# $evigene/scripts/rnaseq/trasmstrandspan.pl -out rungnd$i -ver n3 -pid=ptf00$i \
#  -span grNDtCOcX.strandspan.tab \
#  -bams grNDt*cX.{fwd,noo}.chr.bam daph*_est{sr,pe}.bam > & rungnd$i.log &
## rev strand
# $evigene/scripts/rnaseq/trasmstrandspan.pl -out rungnd$i -ver n3 -pid=ptr00$i 
# -span grNDtCOcX.strandspan.tab \
# -bams grNDt*cX.{rev,noo}.chr.bam daph*_est{sr,pe}.bam > & runrnd$i.log &

##.. need trasmspan.config file for several paths: evigene, gmap, velvet, samtools; data
datad=$HOME/scratch
workd=$datad/chrs/daphmag/
rund=$workd/rnas/$subd
sdir=$HOME/bio/evigene/scripts/rnaseq

# export samtools or put on path
module add samtools

cd $rund

echo "start trasmspan : `date`"  
i=0; # cpu counter
j=$firstpart; # data part counter
while [ $j -le $lastpart ]; do
{
  nam=trasm${dv}p$j
  if test -d $nam; then echo "done $nam"; continue; fi

  # will exit -1 if jpart > parts found in spantable; use other signal?
  $sdir/trasmstrandspan.pl -part=$j -strand=$strand -out $nam  -ver $dv \
    -config=$config -span=$spantable -bams $bamlist  >& $nam.log &

  j=$(( $j + 1 ))
  i=$(( $i + 1 ))
  if [ $i -ge $ncpu ]; then wait; i=0; fi

} done
wait

echo "end trasmspan: `date` "


=cut

=item method for strand span transcript assembly

     .. dont need to do fwd,rev together; one strand steps 1-3  or 1-4b, as script?
     
  1.  split strandgroups/grNDtCOcX.strandspan.tab into partitions, ~ 200k wide, forward/reverse read locations
       .. do this for fixed window? and/or slide wind so have some overlap?  break at 0,0,0 no reads?
       grep '^scaffold00389' grNDtCOcX.strandspan.tab | egrep '        (f|b|rf)$' | head -92  >  grNDtCOcX.sc389f1.tab
       grep '^scaffold00389' grNDtCOcX.strandspan.tab | egrep '        (r|b|fr)$' | head -105 >  grNDtCOcX.sc389r1.tab
  
  2.  select sam partition, fwd/rev, sorting for pair names
    samtools view -L grNDtCOcX.sc389r1.tab grNDtCOcX.rev.chr.bam | sort > grNDtCOcX.sc389r1.rev.sams &
    samtools view -L grNDtCOcX.sc389r1.tab grNDtCOcX.noo.chr.bam | sort >> grNDtCOcX.sc389r1.rev.sams &
    
  2b. and est
    samtools view -L grNDtCOcX.sc389f1.tab ../../bamest/daphmagna_est09.bam | grep 'XS:A:+' > estmag09.sc389f1.fwd.sam
    samtools view -L grNDtCOcX.sc389f1.tab ../../bamest/daphpulex_est05pe.bam | grep 'XS:A:+' | sort > estplx05.sc389f1.fwd.sams
    
    samtools view -L grNDtCOcX.sc389r1.tab ../../bamest/daphmagna_est09.bam | grep 'XS:A:-' > estmag09.sc389r1.rev.sam
    samtools view -L grNDtCOcX.sc389r1.tab ../../bamest/daphpulex_est05pe.bam | grep 'XS:A:-' | sort > estplx05.sc389r1.rev.sams
       
  3. velvet asm rna .. takes only 30 sec on 1 cpu, 100k part, 5,6 kmers
    env strand=fwd ./runvelspan10sg.sh > & log.v10f &
    env strand=rev ./runvelspan10sg.sh > & log.v10r &
  
  4.  collate transcripts, map to genome, filter out junk
        -- improve filtering out subset tr; see nosub=1 gff2ggb.pl
  
  4a. trmake
    ./veltrmake.sh vel10ss389rev_??/transcripts.fa
    ./veltrmake.sh vel10ss389fwd_??/transcripts.fa
  
  4b. trmap
    -- note uses intronsc389.gff, genefindcds, other to score best trs
   env strand=+ ./trgmap.sh vel10ss389fwds.tr
   env strand=- ./trgmap.sh vel10ss389revs.tr
  
    .. improve filtering out subset tr; see nosub=1 gff2ggb.pl
    .. ?? want filter by quality?
    cat vel10ss389fwds.an.gff | grep mRNA | perl -ne \
    '($al,$ap,$ac)=m/aalen=(\d+).(\d+)..(\w+)/; ($ig,$ib,$iz)= m/inqual=(\d+).([-\d]+).(\d+)/;  \
    ($chi)=m/chimera=(\w+)/; ($d)=m/ID=(\w+)/;  print "$d\n"  unless( $ap<40 or abs($ib)>$ig);' \
     | wc  # > good.ids
          
  4c. view map
    cat vel10ss389revs.an.gff | grep '      exon' | env nosub=1 s=v10r $evigene/scripts/gff2ggb.pl > dmag3vel10s389p1.ggb
    cat vel10ss389fwds.an.gff | grep '      exon' | env nosub=1 s=v10f $evigene/scripts/gff2ggb.pl >> dmag3vel10s389p1.ggb
  
  #................................................................................

=cut

