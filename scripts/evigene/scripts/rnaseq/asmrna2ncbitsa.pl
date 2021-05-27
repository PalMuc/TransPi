#!/usr/bin/env perl
# asmrna2ncbitsa.pl

=item about
  
  clean transcript assembly cdna for NCBI TSA submission

  -Sequences containing vector contamination 
    -- use ncbi/blast/vecscreen, parse results for Strong|Medium vec spans, turn vec spans to NNN, then chomp as below
  * -Sequences less than 200bp
  * -Sequences containing more than 10% n's or greater than 14 n's in a row
    -- also new gap restriction: SEQ_INST.HighNContentStretch: stretch of at least 5 Ns within the last 10 bases
  -Annotation errors if annotation is provided(internal stops codons, etc.)
  -General formatting errors (missing tissue type, BioProjectID, SRR numbers, assembly data, or comment)
  -Product name formatting errors if annotation is provided

=item author
  
  don gilbert, gilbertd near indiana edu, 2012
  part of EvidentialGene, evigene/scripts/

=cut

use strict;
use Getopt::Long;

my $debug= 0;
my $GENEINFO_VERS=2; # constant
my $MAXGAP=15; # NCBI 
my $ENDGAP=10; # trim ends if gaps w/in this of ends; NCBI, was 10
my $MINSIZE=200; # NCBI
my $GAPSOK=0; # default on? new policy 2012Dec for TSA tbl2asn: -a r10u -l paired-ends
my $MINGENEIDENT=85; # for asm == gene identity
my $IDPREFIX= $ENV{idprefix} || ""; #  "Thecc1RA_"; # opt
my $GDB_PREFIX= $ENV{gbprefix} || 'gnl|CacaoGD|'; # .tbl annot, need option; no default?
my ($vecscreen,$geneinfo,$cdnaseq,$output,$logfile,$tblfile);

my $optok= GetOptions(
  "cdnaseq=s", \$cdnaseq,  
  "output=s", \$output,  
  "tblfile=s", \$tblfile,  
  "logfile=s", \$logfile,  
  "idprefix=s", \$IDPREFIX,  
  "gdbprefix=s", \$GDB_PREFIX,  
  "vectors|vecscreen=s", \$vecscreen,  
  "geneinfo=s", \$geneinfo,  
  "MINSIZE=i", \$MINSIZE,  
  "MAXGAP=i", \$MAXGAP,  
  "MINGENEIDENT=i", \$MINGENEIDENT,  
  "GAPSOK!", \$GAPSOK, 
  "debug!", \$debug, 
  );

die "usage: asmrna2ncbitsa.pl -cdna cdna.fa|stdin > cdna4tsa.fsa
 opts: -idprefix Thecc1RA_ -vecscreen=infile -geneinfo=infile -MINSIZE=$MINSIZE  -MAXGAP=$MAXGAP  -out=outfasta  -log=outlog
" unless($optok and  $cdnaseq);

my( %vec, %gene, %nam, %geneinfo);
my $GAPSMAX = ('N') x $MAXGAP;

my $logh= undef;
if($logfile) { open(LOG, ">$logfile") or die $logfile; $logh= *LOG; }

my $outh= *STDOUT;
if($output and $output!~/stdout|^-/) { open(OUT, ">$output") or die $output; $outh= *OUT; }

my $tblh= undef; # genbank submit annot
if($tblfile) { open(TBL, ">$tblfile") or die $tblfile; $tblh= *TBL; }

  
MAIN: {

=item geneinfo.tab v1

  asmid, percentident, geneid, genename (other info? need cds-start,stop for genbank submit)
  B_g00001c00019	100	Thecc1EG019703t1	Kinase superfamily protein
  B_g00001t00001	100	Thecc1EG036832t1	WRKY protein
  B_g00002c00004	36	Thecc1EG030081t1	
  B_g00002c00030	76	Thecc1EG021938t1	
  
=item geneinfo.tab v2:
  use or/strand=- to revcmp trfa,
  tsa.annotate.tbl uses locustag, product, dbxref, any other?
  # FIXME: for genome-less asmrna: no mapqual, flag as 'nomap' or 'na' or what?

trid                            geneid                  mapqual                         aaqual                  trlen  or       cdsoff          locustag        product                                                 dbxref
cacao5sopcsc10rk23loc1038t1     Thecc1EG043360t1        mapfull:100.97,aaeq:590,100,100 590,61%,complete        2861    -       228-2000        TCM_043360      NAC domain protein, IPR003441, putative isoform 1       TAIR:AT5G04410.1,TrEMBL:UniRef50_B9HUM1,Phytozome:POPTR_0008s03170.1
cacao5sopcsc10rk23loc1116t1     Thecc1EG042722t1        mapfull:100.97,aaeq:367,100,100 367,69%,complete        1585    -       185-1288        TCM_042722      2-oxoglutarate and Fe(II)-dependent oxygenase superfamily protein, putative     TAIR:AT3G21420.1,TrEMBL:UniRef50_D7TLQ8,Phytozome:GSVIVT01016505001
cacao5sopcsc10rk23loc1166t1     Thecc1EG042983t1        mapfull:99.99,aaeq:324,100,100  324,65%,complete        1486    -       274-1248        TCM_042983      Alpha/beta-Hydrolases superfamily protein       TAIR:AT2G39420.1,TrEMBL:UniRef50_Q8RXN7,Phytozome:POPTR_0010s21900.1
cacao5sopcsc10rk23loc1194t2     Thecc1EG045490t1        mapfull:94.78,aaeq:409,97.9,95.3        390,68%,complete        1712    -       67-1239 TCM_045490      OB-fold nucleic acid binding domain-containing protein isoform 1        TAIR:AT2G40660.1,TrEMBL:UniRef50_Q93VB0,Phytozome:GSVIVT01032669001
cacao5sopcsc10rk23loc1301t2     Thecc1EG043107t1        mapfull:99.86,aaeq:574,100,100  574,77%,complete        2233    +       212-1936        TCM_043107      Seven transmembrane MLO family protein  TAIR:AT1G61560.1,TrEMBL:UniRef50_Q94KB7,Phytozome:POPTR_0008s04100.1
  
=cut

if($geneinfo) { 
  my $nerr=0;
  open(F,$geneinfo) or die $geneinfo;
  while(<F>) { 
    chomp; 
    
    if($GENEINFO_VERS == 2) {
	  #	trid	geneid	mapqual	aaqual	trlen	or	cdsoff	locustag	product	dbxref
    my @row= split"\t";
    my($td,$gd,$mapqual,$aaqual,$trlen,$cdsor,$cdsoff,$lotag,$nam,$dbxref)= @row;
    ##    my($td,$gd,$mapqual,$aaqual,$trlen,$cdsor,$cdsoff,$lotag,$nam,$dbxref)=@{$geneinfo{$oid}};

    unless($trlen =~ /\d/ or $trlen !~ /\D/ or $cdsor =~ /^[+-]$/) {
      warn "#ERR.geneinfo.v2 format:$_\n"; $nerr++; die if($nerr>9);
    }
    $row[8]= $nam="" if($nam eq "na"); # missing
    $geneinfo{$td}= \@row;
    $nam{$td}=$nam; $gene{$td}=$gd; #skip this?
    
    } else {
    my($td,$pi,$gd,$nam)=split"\t";  #** add CDS offset cb,ce in tr
    if($pi >= $MINGENEIDENT) { $nam{$td}=$nam; $gene{$td}=$gd; } 
    }
  } close(F);
  # if(/^GeneInfo/) { chomp; ($gg,$td,$pi,$gd,$nam)=split"\t"; $nam{$td}=$nam; $gene{$td}="CacaoGD:$gd"; next; } 
}

=item cdna_vector.tab

  asmid, start, end, type
  B_g02526t00001	450	475	UniVec, Moderate match
  B_g02526c00001	263	299	UniVec, Strong match
  B_g03099t00001	1	27	UniVec, Moderate match
  
  ** bug here, can have 2+ vector lines; join or output as is?

 $ncbi/bin/vecscreen -i cdna.fa -d $ncbi/data/UniVec -f3 | perl -ne \
'chomp; if(/^>/) { putv() if($id); ($id)=m/>Vector (\S+)/; ($vd)=m/Database: (\S+)/; $vb=$ve=$so=$ty=0; } \
elsif(/^No hits/) { $id=0; } elsif(/^\w/) {  \
  if(/match/) { $ty=$_; } elsif(/^Suspect origin/) { $so=1; } \
  elsif(/^(\d+)\s+(\d+)/) { ($b,$e)=($1,$2); \
    if($so and $ve) { $vb=$b if($b<$vb); $ve=$e if($e>$ve); }  \
    elsif($ty) { ($vb,$ve)=($b,$e); }  } } END{ putv() if($id);}\
sub putv { print join("\t",$id,$vb,$ve,$ty,$vd)."\n" if($id and $ty and not $ty=~/Weak/); }'  \
 > cdna_vector.tab
 
  ** bug here, can have 2+ vector lines; join or output as is?
  ** also revcomp here created one apparent vector match.
  
>Vector cacao4vel12sc4k51Loc434t1 Screen: VecScreen Database: UniVec (build 6.0)
Strong match
8       41
Moderate match
1       20

>Vector cacao5sopcsc7k49loc670t1 Screen: VecScreen Database: UniVec (build 6.0)
Strong match
2834    2867
Moderate match
2855    2877

revcomp new
>Vector Thecc1ER_nwbL_g06984t00001 Screen: VecScreen Database: UniVec (build 6.0)
Moderate match
1382    1400
Suspect origin
1401    1414

  
=cut

if($vecscreen) { 
  open(F,$vecscreen) or die $vecscreen;
  # FIXME: multi locs per vd; .= val\n
  while(<F>) { chomp; my($vd,$vb,$ve,$vt)=split"\t"; $vec{$vd}.="$vb\t$ve\t$vt\n"; } close(F);
}

my($fa, $hd, $fh);
if($cdnaseq =~ /stdin|^-/) { $fh= *STDIN; }
elsif($cdnaseq =~ /\.gz/) { open(F,"gunzip -c $cdnaseq|") or die $cdnaseq; $fh=*F; }
else {open(F,$cdnaseq) or die $cdnaseq; $fh=*F; }
while(<$fh>) { 
  if(/^>(\S+)/) { my $d=$1; putseq($hd,$fa) if($fa); $hd=$d; $fa=""; }
  elsif(/^\w/) { chomp; $fa.=$_; }
} close($fh);

putseq($hd,$fa) if($fa); # last

}


#-----------------------------

sub cleanid { 
  local $_= $_[0];
  ## project specific cleaning : read from where? cmdline or properties?
  s/^cacao3(vel|v)(\d+)/vel$2/; # do we still need this?
  s/^caca11r39cuf(\d+)/cuf$1/;  
  s/^(B|L|P1|P2)/nwb$1/;
  if($IDPREFIX and !/^$IDPREFIX/) {
    s/^cacao[345]//; # FIXME; option?
    $_= $IDPREFIX.$_;
  }
  return $_;
}
  
sub putseq 
{
  my ($oid, $fa)=@_; 
  my $id= cleanid($oid);
  my $def= $id;
  my($ncut, $vw)=(0,0);

  my($vec,$gn,$nam,$rinfo,$cdsb,$cdse,$tlog);
  if($vec=$vec{$oid}) { 
    # FIXME: multi locs per vd; sep by \n
    my @vec=split"\n",$vec; my $nv= @vec;
    foreach $vec (@vec) {
      my($vb,$ve,$vt)=split"\t",$vec; $vw=1+$ve-$vb; substr($fa,$vb-1,$vw)= ("N") x $vw; 
      }
    } 
  $tlog="";  ($cdsb,$cdse)=(0,0);
  # geneinfo for ncbi tsa submit: move this to feature.tbl not cdna.fsa; add cds offsets
  # >Feature Thecc1RA_L_g13025t00001
  #  bg  eg  gene  locus_tag TCM_000002
  #  bc  ec  CDS   product Cystathionine beta-synthase (CBS) family protein isoform 1
  #                protein_id  gnl|CacaoGD|Thecc1EG000002p1
  
  if($GENEINFO_VERS == 2) {
	  #	trid	geneid	mapqual	aaqual	trlen	or	cdsoff	locustag	product	dbxref
	  ## do revcomp and rev-offset for cdsor=-
    $rinfo= $geneinfo{$oid} || [];
    my($td,$gd,$mapqual,$aaqual,$trlen,$cdsor,$cdsoff,$lotag,$nam1,$dbxref)= @$rinfo;
    $gn=$gd; $nam=$nam1;
    
    ($cdsb,$cdse)= split/[-]/,$cdsoff; # this SHOULD be offset in non-rev tr; any problems?
    if($cdsor eq "-") {
      my $farev= revcomp($fa); 
      $fa= $farev; #? save orig
      # my($cb,$ce)= split/[-]/,$cdsoff; # this SHOULD be offset in non-rev tr; any problems?
      ## ** FIXME dammmit; if cb > ce, then this is offset BEFORE rev; if cb < ce is offset AFTER rev **
      my($cbr,$cer)=(0,0);
      if($cdsb > $cdse) {  ## this is new vers
        $cbr= 1+$trlen - $cdsb;
        $cer= 1+$trlen - $cdse;
      } else { # DAMMIT mess w/ cdna_org version changed offs= meaning; this is old vers
        $cer= $cdse; ## 1+$trlen - $cdsb;
        $cbr= $cdsb; ## 1+$trlen - $cdse;
      }
      ($cdsb,$cdse)=($cbr,$cer);    
      $tlog.="revcds:$cdsb-$cdse,";  
      $cdsoff="$cdsb-$cdse"; # update !!!
      $rinfo->[6]= $cdsoff;
      $rinfo->[5]= $cdsor= "+";
    }
        
  } else {
    if($nam=$nam{$oid}) { $def.=" [product=$nam]"; }
    if($gn=$gene{$oid}) { $def.=" [protein_id=$gn]"; } 
  }
  
  my $ol= length($fa); 
  my $nn= $fa =~ tr/Nn/Nn/; 
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
    my $n1= index($fa,'N'); if($n1 <= $ENDGAP) { 
      $n1++; $fa= substr($fa,$n1);  if($fa=~s/^(N+)//) { my $ncut=length($1); $n1+=$ncut; }
      if($cdsb>0) { $cdsb -= $n1; $cdse -= $n1; }
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
        if($en < $cdsb) { $cdsb -= $cut; $cdse -= $cut; }
        elsif($in < $cdse and $en > $cdsb) {
          if($in <= $cdsb) { $cdsb -= $cut; $cdse -= $cut; } #??
          else { $cdse -= $cut; }
        }
      }
      $in=index($fa,$GAPSMAX); 
    } 
    }
    
    unless($cdse == $lcdse and $cdsb == $lcdsb) {
      $rinfo->[6]= "$cdsb-$cdse"; # update !!!
      $tlog.="cutcds=$cdsb-$cdse,oldcds=$lcdsb-$lcdse,"; 
    }
  } 
  my $nl= length($fa); 
  my $nn1= $fa=~tr/N/N/; 
  my $len=($nl==$ol and $nn==0)?$nl:"$nl; olen=$ol; cut=$ncut; nnn=$nn1/$nn;"; 
  $tlog.="vectrim=$vw," if($vw); 
  if($nl<$MINSIZE) {
    print $logh "#$def\toid=$oid\tgene=$gn\tlen=$len\t$tlog\tERR: too short:$nl\n" if($logfile);
    return "ERR: too short:$nl" 
  }

  if($GENEINFO_VERS == 2 and $tblh) {
	  #	trid	geneid	mapqual	aaqual	trlen	or	cdsoff	locustag	product	dbxref
	  ## do revcomp and rev-offset for cdsor=-
	  # use from above# my $rinfo= $geneinfo{$oid} || [];
    my($td,$gd,$mapqual,$aaqual,$trlen,$cdsor,$cdsoff,$lotag,$nam,$dbxref)= @$rinfo;

    $tlog.="mapq=$mapqual,";
    map{ $_="" if($_ eq "na"); } ($gd,$lotag,$nam,$dbxref);
    # $gd="" if($gd eq "na"); # not all have this
    # $lotag="" if($lotag eq "na"); # not all have this
    my $protid= $id; # $td; # $gd ? gbasn requires protein_id here, use genes.tbl id? or tr id?
    $protid=~s/t(\d+)$/p$1/; 
    $protid= $GDB_PREFIX.$protid if($protid);
    #  protein_id      gnl|CacaoGD|Thecc1EG016762p1
    
    # require aaqual =~ /complete/ and $mapqual =~ /mapfull:/ and $mapqual =~ /aaeq|aasim/
    # FIXME: for genome-less asmrna: no mapqual, flag as 'nomap' or such
    unless($aaqual =~ /complete/ and $mapqual =~ /mapfull/ and $mapqual =~ /aaeq|aasim/) {
      $dbxref=""; # dont claim this annot unless hi qual gene match
      ## FIXME: SEQ_FEAT.MissingCDSproduct
      if($nam =~ /\w/ and $nam ne "na") { $nam .= " fragment"; } 
      else { $nam="hypothetical protein"; }
      #? drop cds if not complete? else have to deal w/ partial cds: "<cdsb >cdse CDS"
      $cdsoff=$lotag="" unless($aaqual =~ /complete/); #?? drop CDS or not
    }
    
    ## require cdsoff, skip if have only locustag
    if($cdsoff =~ m/\d+\-\d+/ and $lotag) { ##  or $lotag
      print $tblh ">Features\t$id\tcacao11evigene_20120827\n\n";
      # >Features       Thecc1ER_sopcsc10rk23loc140t1     cacao11evigene_20120827
      
      if(0 and $lotag) {
        print $tblh "1\t$trlen\tgene\n";
        print $tblh "\t\t\tlocus_tag\t$lotag\n"; #? no gene row? put in CDS?
        print $tblh "\n";
        }
      if($cdsoff =~ m/\d+\-\d+/) {
        $cdsoff =~ s/\-/\t/;
        print $tblh "$cdsoff\tCDS\n";
        print $tblh "\t\t\tprotein_id\t$protid\n" if($protid);
        print $tblh "\t\t\tlocus_tag\t$lotag\n" if($lotag);
        print $tblh "\t\t\tproduct\t$nam\n" if($nam);
        # ^^ FIXME: SEQ_FEAT.MissingCDSproduct; need some nam: hypothetical??? or orig-name + similarto/
        if($dbxref =~ /\w/) {
        foreach my $dx (split",",$dbxref) { print $tblh "\t\t\tdb_xref\t$dx\n" if($dx=~/\w/ and $dx ne "na"); }
        }
        print $tblh "\n";
      }
    }
  }
  
  $fa =~ s/(.{60})/$1\n/g; 
  # defline option? strip  oid=..; len=.. from fahead before tbl2ans
  # optional 2nd outfile for  newid, oldid, len, changes..
  if($logfile) {
    print $logh ">$def\toid=$oid\tgene=$gn\tlen=$len\t$tlog\n";
    print $outh ">$def\n$fa\n";
  } else {
    if($debug) { print $outh ">$def oid=$oid; len=$len\n$fa\n"; }
    else { print $outh ">$def\n$fa\n"; }
  }
  return 0;
}

sub revcomp {
    my ($seq) = @_;
    my $reversed_seq = reverse ($seq);
    $reversed_seq =~ tr/ACGTacgtyrkm/TGCAtgcarymk/;
    return ($reversed_seq);
}

__END__

=item tsa annot

  -- revise geneinfo input for fuller tsa annots (as .tbl output per gene.tbl)
     -- use instead evigene2genbanktbl.pl with genes.gff input for tsa.tbl annots? will need revisions
     -- or copy/paste valid CDS annots from genome/genes.tbl?
  -- for CDS annots, option to revcomp flip tr when aa is on -strand
  tsasub5/TCM01.tsa_rasm.val 
  ERROR: valid [SEQ_FEAT.CDSonMinusStrandMRNA] CDS should not be on minus strand of mRNA molecule FEATURE: CDS: hypothetical protein [lcl|Thecc1RA_cacao5sopcsc10rk23loc28t1:c2486-696] [lcl|Thecc1RA_cacao5sopcsc10rk23loc28>: raw, rna len= 2554] -> [lcl|Thecc1RA_cacao5sopcsc10rk23loc28t1_1]

=item tsasub11 cacao

  evigene=/Users/gilbertd/Desktop/dspp-work/genes2/evigene
  ncbin=/Users/gilbertd/Desktop/dspp-work/genomesoft//ncbi/bin
  nwbsra='[SRA=SRR531454; SRR531455; SRA058778; SRA058779]'
  vel3sra='[SRA=SRA058777; SRA058780; SRA058781; SRA058782; SRR531454; SRR531455; SRA058778; SRA058779]'
  rnasra='[SRA=SRA058777; SRA058780; SRA058781; SRA058782]'
    ## fixme, vel3sra= both sets: ill + 454
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
  

  
  
=item example tsa annot

  ** maybe change asmid to pub3i id, w/ other prefix ?? 
  --- but only for full matches ... otherwise need new Thecc1ER ids ... easier to keep oid

  pub3ig.trasm.match.tab4e.best1:
  Thecc1EG045554t1 cacao5sopcsc10rk23loc140t1      full:100.93     aaeq:505,100,100
  aa>cacao3sopcsc10rk23loc140t1 aalen=505,73%,complete; clen=2054; strand=-; offs=323-1840; << diff offs after rev.
        revoffs= 1+2054 - 1840=215  .. 1+2054 - 323=1732
  over>Thecc1EG045554t1        AUGepir1a cacao5sopcsc10rk23loc140t1/C100.93      10rsc:25151751-25154909:-
  gmap> scaffold_10r:25151806-25154807:-  ID=cacao3sopcsc10rk23loc140t1;aalen=505;cov=100.0;match=2054;nexon=3;pid=100.0;qlen=2054;sense=-1

  ... input geneinfo table for tsa annot:
  trid                        geneid            mapequal    aaequal           aaqual           cdsoff   cdstrand  trlen
  cacao5sopcsc10rk23loc140t1  Thecc1EG045554t1  full:100.93 aaeq:505,100,100  505,73%,complete 323-1840  -        2054
  locustag    aadbxref                                    aaproduct    << from genome/genes.tbl for Thecc1EG045554t1 CDS
  TCM_045554  TAIR:AT2G30490.1,TrEMBL:UniRef50_P92994     Cinnamate-4-hydroxylase
  --  no CDS tbl annot unless mapequal==full and aaequal==aaeq ??
  ..  or maybe add CDS as "hypothetical protein" when aaqual=complete ??
  
  #   [2] Since we will use PID 51633 for both submissions this 
  #   will allow you to use the same locus_tag in both your genome 
  #   and TSA submissions.  Therefore, you are not required to 
  #   edit your locus_tags.  
  
  ... tsa.an.tbl
  >Features       Thecc1ER_sopcsc10rk23loc140t1     cacao11evigene_20120827
  1   2054  gene
          locus_tag       TCM_045554    # use this only as link w/ genome gene annot?
  323 1840  CDS
          ##no# transcript_id   gnl|CacaoGD|Thecc1EG045554t1    # Maybe; probably not; this is source ID
          ##no# protein_id      gnl|CacaoGD|Thecc1EG045554p1    # Maybe crossref genes.tbl, but not as protid, as dbxref? 
          ##? protein_id      gnl|CacaoGD|Thecc1ER_sopcsc10rk23loc140p1  # local protid for any CDS feat w/ product? or leave out?
          product Cinnamate-4-hydroxylase   # YES
          db_xref TAIR:AT2G30490.1   # Maybe
          db_xref TrEMBL:UniRef50_P92994   # Maybe

  ... gene.an.tbl
  tsa: NO mRNA or gene tbl, but source == mRNA, CDS only useful annot?
  # 25154909        25151751        gene   #YES
                        locus_tag       TCM_045554    # Maybe

  # 25154909        25153701        mRNA : 
                        transcript_id   gnl|CacaoGD|Thecc1EG045554t1
                        locus_tag       TCM_045554
                        note    quality class is strong, express is strong, homology is ortholog strong, intron is strong, protein is complete
                        product Cinnamate-4-hydroxylase
                        inference       similar to AA sequence:frave:gene28093
                        inference       similar to AA sequence (same species):CacaoGD:Thecc1EG042294t1
                        note    similar to RNA sequence, EST 0.41 coverage
                        protein_id      gnl|CacaoGD|Thecc1EG045554p1
                        inference       alignment:exonerate-protein:2:Phytozome:POPTR_0013s15380.1
                        note    inference alignment:Newbler:10:TSA:cgbaP2_g03566t00001
                        # ^^ this TSA will replace this genes.tbl mRNA inference

  25154485        25153701        CDS
                        transcript_id   gnl|CacaoGD|Thecc1EG045554t1    # Maybe
                        product Cinnamate-4-hydroxylase   # YES
                        db_xref TAIR:AT2G30490.1   # Maybe
                        db_xref TrEMBL:UniRef50_P92994   # Maybe
                        #NO note    dbxref frave:gene28093
                        #NO note    paralog of  gnl|CacaoGD|Thecc1EG042294t1
                        #no crossref genes.tbl# protein_id      gnl|CacaoGD|Thecc1EG045554p1    # Maybe
  
... tsa.gbf from this;  -kc makes CDS, want own annot
 $ncbin/tbl2asn -kc -Vtb -as -Mt -XE -t cacao3i_tsasubmit.sbt -p tsasub5/ -Z tsasub5/discrep.log \
 -j '[tissue_type=leaf,pistil,bean] [moltype=mRNA] [tech=TSA] [organism=Theobroma cacao] [cultivar=Matina 1-6]' \
          
LOCUS       Thecc1RA_cacao5sopcsc10rk23loc140t1 2054 bp   mRNA linear TSA 06-NOV-2012
DEFINITION  TSA: Theobroma cacao cultivar Matina 1-6, mRNA sequence.
ACCESSION   
VERSION
DBLINK      BioProject: PRJNA51633
KEYWORDS    TSA; Transcriptome Shotgun Assembly.
SOURCE      Theobroma cacao
  ORGANISM  Theobroma cacao
            Unclassified.
FEATURES             Location/Qualifiers
     source          1..2054
                     /organism="Theobroma cacao"
                     /mol_type="mRNA"
                     /cultivar="Matina 1-6"
                     /tissue_type="leaf,pistil,bean" << FIXME: "leaf; pistil; bean", or drop ttype
     CDS             complement(215..1732)
                     /codon_start=1
                     /product="hypothetical protein"
                     # Thecc1EG045554t1        Cinnamate-4-hydroxylase
                     /translation="MDLLLLEKALVSLFITVILAILISKLRSKRFRLPPGPIPVPVFG...

=cut


=item prelim script

  cat genes.tab vectors.tab cacao3i_tsaok.fsa | perl -ne
  'if(/^GeneInfo/) { chomp; ($gg,$td,$pi,$gd,$nam)=split"\t"; $nam{$td}=$nam; $gene{$td}="CacaoGD:$gd"; next; } 
  if(/\tUniVec/) { chomp; ($vd,$vb,$ve,$vt)=split"\t"; $vec{$vd}="$vb\t$ve\t$vt"; next; } 
  if(/^>(\S+)/) { $d=$1; putf() if($fa); $hd=$d; $fa=""; } elsif(/^\w/) { chomp; $fa.=$_; }
  END{ putf()} sub putf{ $id=$hd; $id=~s/^cacao3(vel|v)(\d+)/vel$2/; 
  $id=~s/^caca11r39cuf(\d+)/cuf$1/; $id=~s/^(B|L|P1|P2)/nwb$1/; $def="Thecc1RA_$id"; 
  if($vec=$vec{$hd}) { ($vb,$ve,$vt)=split"\t",$vec; $vw=1+$ve-$vb; substr($fa,$vb-1,$vw)=("N") x $vw; } 
  if($nam=$nam{$hd}) { $def.=" $nam,"; } if($gn=$gene{$hd}) { $def.=" $gn"; } 
  $ol=length($fa); $nn= $fa=~tr/N/N/; if($nn>0) { 
  $fa=~s/^N+//; $fa=~s/N+$//; $ncut=0; $in=index($fa,"NNNNNNNNNNNNNNN"); 
  while($in > 0) { $w=length($fa); $en=$in+15; $en++ while($en<$w and substr($fa,$en,1) eq "N"); 
  $wn= $en-$in; $keep= 3 + ($wn % 3); $cut= $wn-$keep; $ncut+=$cut; 
  $nb=substr($fa,0,$in).substr("NNNNNN",0,$keep).substr($fa,$en); 
  $fa=$nb; $in=index($fa,"NNNNNNNNNNNNNNN"); } } 
  $nl= length($fa); $nn1= $fa=~tr/N/N/; $len=($nl==$ol)?$nl:"$nl; olen=$ol; cut=$ncut; nnn=$nn1/$nn;"; 
  $len.=" vectrim=$vw" if($vec); 
  $fa =~ s/(.{60})/$1\n/g; print ">$def oid=$hd; len=$len\n$fa\n" if($nl>199); }' 
  > cacao3i_tsa_devecgap.fsa

=cut

