#!/usr/bin/env perl
# cdna_bestorf.pl cut from genefindcds.pl

=item about
  
  finds best orf, adds proteins and CDS to transcript assembly cdna
    cdna_bestorf.pl  -cdna transcript.fasta  > transcript.aa

=item todo
   
   FIXMEd: read all tr 1st is bad idea, huge.tr set chokes
   fixed: new option: noutrorf=turn off utrorf output, see splitutrorfs

	 ** add hdr report of longest revaa (2nd longest?) unless fwd/rev spec.
	 -- some (rare?) cases of shorter revaa is true prot; homology &/or gmap intron splicing validate this
	 
   clean up "action" options; separate opts
   add frameshift detect option (cdna_proteins.pm)
   
=item FIXME 201505 annot mismatch bugs

  * mismatching annots for -outaa -outmrna for cases of revc=1 and for offset=high-low;
    ie., where in.cdna is not mrna +strand to bestorf, needs revcomp for -outmrna,
      outaa annot has pre-rev values: offs=456-123; strand=-; 
      but -outmrna has revc=1; offs=123-456; strand=+ and is revcomp(incdna)
    also sometimes not above, but offs=456-123 problem, due to 2 bestorf (utrorf), but -noutrorf output,
      .. 1st outmrna has 2nd utrorf offs=bad (only that?), outaa has proper 1st offs ?
  
=item FIXME 201402 orthobest prots
  
  orthoprot refs x all asmtr finds ~1% .. 5% cases where cdna_bestorf fails to pick best ortho-orf
  can be parsed into table of diff.aa showing where 
   a. revaa is orthobest, replace with shorter+rev than longest-full-orf (sometimes much shorter, 50%..)
   b. overlap-extended orthobest, replace best aa is extended over longest-full-orf, sometimes with inner-stop, or full>partial
   c. overlap-reduced orthobest, replace best aa is shorter than longest{partial} orf
   d. utrorf(rev/fwd) has orthobest, main-longest-orf also still valid, use to ensure utrorf is annotated output.
   
  steps 
    0. cdna_bestorf first run, allasm.tr > allasm_bestorf.aa,cds
    1. tblastn refprot x allasm.tr, have orthobest table of alignments, cds-span:orient
    2. diff table of first bestorf-aa-spans x orthobest-align-spans 
    3. cdna_bestorf rerun, adding input of orthobest-diffaa.tab, output new aa,cds,mrna seq files w/ diffaa updates
        - need to call orf given orthobest-cds-span:orient, including inner-stops, extend out from part-align ends.
        - annotate orthobest types (e.g. ordiff:a,b,c,d above)

  input2: diff table like this?
  head catfish1evg8tr-refkfish.diffaa.tab2 | cut -f1,2,5,8,9,10 | head
  # RefID             TrID                       nAlign   refaaSpan  trAlignSpan  Class:aaOrfSpan
  Funhe2EKm000007t1   socatfishv2k31loc62810t1    326     1-325      2341-1364:-  revutr:3589-7470:+  # utrorf
  Funhe2EKm000108t1   socatfishv1k25loc49492t1    923     23-919     1-2706:+     over:22-2709:+   # small over, 20n
  Funhe2EKm000123t1   socatfishv1k21loc199140t1   500     1077-1575  22-1488:+    over:172-1491:+  # middle ovr, 150n
  Funhe2EKm000125t1   socatfishv1k25loc84696t1    525     51-739     177-1989:+   over:198-1919:+  # small ovr, 20n
  Funhe2EKm000125t1   socatfishv1k21loc76166t1    549     51-739     177-1985:+   over:198-1988:+  #  "
  Funhe2EKm000154t1   catfishtrin1loc50544c1t1    108     18-124     983-660:-    revaa:509-850:+  # revaa ~50n shorter
  Funhe2EKm000198t1   socatfishv1k31loc529495t1   129     1-117      95-478:+     revaa:559-8:-    # rev ~80n shorter
    -- use refaaSpan to guide bestorf to ~max extention from trAlignSpan ? need refaa-len?

  # aa.spans from bestorf1.aa hdr info.
   grep '^>' catfish1evg8.aa | perl -ne \
  '($id)=m/>(\S+)/; ($aw,$cl,$or,$ofs)=m/(?:aalen|clen|strand|offs)=([^;\s]+)/g; print join("\t",$id,$aw,$cl,$ofs,$or)."\n"; ' 
  > catfish1evg8.aa.spans
  
  cat catfish1evg8.aa.spans | sed 's/^/span /;' | cat - catfish1evg8tr-refkfish.tall1 | perl -ne \
  'if(/^span /) { ($ss,$id,$aq,$cl,$spa,$or)=split; $spa{$id}=$spa; $or{$id}=$or; $aq{$id}=$aq; } 
  else { ($rd,$td,@v)=split; ($bs,$idn,$aln,$rl,$tl,$rspa,$tspa)=@v; ($tbe,$tor)=split":",$tspa; 
  ($tb,$te)=split"-",$tbe; $aspa=$spa{$td}; ($ab,$ae)=split/[:-]/,$aspa; $aor=$or{$td};
  ($ab,$ae)=($ae,$ab) if($ae<$ab); ($tb,$te)=($te,$tb) if($tb>$te); $tw=$te-$tb; $aw=$ae-$ab; $ovt=0; 
  if($tb < $ae and $te > $ab){ $td=abs($tb-$ae); $ad=abs($ab-$te); $d=($td<$ad)?$td:$ad; $ovt=$d if($d/$aw >0.1); }
  $fl=""; if($aor ne $tor){ $fl=($ovt)?"revaa":"revutr"; } elsif($tw > $aw+9) { $fl=($ovt)?"over":"utr"; } 
  if($fl) { s/$/\t$fl:$aspa:$aor/; print; } }' \
   > catfish1evg8tr-refkfish.diffaa.tab2


          
=item author
  
  don gilbert, gilbertd near indiana edu, 2011
  part of EvidentialGene, evigene/scripts/
  cut from genefindcds.pl, which is drawn largely from Brian Haas's PASA scripts

=cut

use FindBin;
use lib ("$FindBin::Bin");

use strict;
use warnings;
use Getopt::Long;
use cdna_proteins;

use constant VERSION  => '20150526'; # TRIMCDSENDGAP added, bugs fixed; outmrna <> outaa annot bug fixes
  # '20141231'; # Selc annots added; mRNA revc=1 annot fixes;
  # '20140317';  # '20131124' 20130818 0228'20120721'; 

my ($cdnaseq, $nin, $skipsubsetseq, $splitutrorfs, $noutrorf, $revaa_report)= (0) x 20;

my $action="fasta";  # or cdnafasta  cdsfasta ..
my $debugin= undef;
my %cdnaseq; my %cdnahead; my %orthobest;

# add -output option to file updated genes
my($oformat,$ostrand,$opart,$orthobest)=("") x 10;
my($aaout,$cdsout,$cdnaout,$mrnaout)= (undef) x 10;
my $FIXUPMRNA=0; # hack or update old to new annots: Selc and revc=1

$InnerStopToX=0;  # bug 1505 getting best w/ innerstops now, why?
$USESelenocysteine=0;
# $TRIMCDSENDGAP=1; # bug 1505 need debug opt = 0 

my $optok= GetOptions(
  "cdnaseq=s", \$cdnaseq,  
  "aaseq|outaa|output:s", \$aaout, #? make default -out=infile.cds unless -out stdout
  "cdsseq|outcds:s", \$cdsout, #? make default -out=infile.cds unless -out stdout
  "outcdna:s", \$cdnaout, #? make default -out=infile.cds unless -out stdout
  "outmrna:s", \$mrnaout, #? same as outcdna but +strand from bestorf (cds is always +strand)
  "action=s", \$action,    # need to specify actions: 
      # now packing on:  "fwd|rev" strand + cdna|cds|noaa,fasta|notfa + full|longpart|bestpart
   "oformat=s", \$oformat,    # replace action: 
   "ostrand=s", \$ostrand,    # replace action: 
   "opart=s", \$opart,    # replace action: 
     
  "orthobest=s", \$orthobest,    # 201402 input table
  "minaa=i", \$MINAA,  # not used ?
  "fullorf:s", \$ORF_FULLvPART,  
  "nostopcodon", \$NoStopCodon, 
  "trimcdsgap!", \$TRIMCDSENDGAP,  # default = 1, use notrim
  "splitutrorfs!", \$splitutrorfs, 
  "noutrorf", \$noutrorf, 
  "revaareport", \$revaa_report, # dont need special opt? debug/verbose/action or default?
  "goodlen!", \$USEGOODLEN,  # add MINGOOD
  "Selenocysteine|selc!", \$USESelenocysteine,  
  "fixmrna|FIXUPMRNA!", \$FIXUPMRNA,  #  
  "goodmin=s", \$MINGOOD,  
  "debug:i", \$debugin, 
  );

die "usage: cdna_bestorf.pl -cdna cdna.fa > protein.aa
 opts: -minaa=$MINAA  -nostopcodon -goodmin=$MINGOOD 
      -outaa -outcds -outmrna -outcdna : output files (input.aa,.cds,.mrna,.cdna default) 
      -noutrorf : ignore long utr-orfs, -splitutrorfs : make cds/cdna from utr for utrorf
      -action=fasta|table, [cdna|cds][only]fasta, [fwd|rev|both]strand, [full|long|best]part
" unless($optok and  $cdnaseq);

if(defined $debugin) { $DEBUG=($debugin>0)?$debugin:1; } else { $DEBUG=0; }

if($action =~ /reportrevaa|reportrev|report/) { # messy place to add this; revfasta conflict
	$revaa_report=1; $action =~ s/reportrevaa|reportrev|report//;
}
$action=$ostrand."strand$action" if($ostrand);
$action=$opart."part$action" if($opart);
$action=$oformat.$action if($oformat);

# "action=s", # now packing on:  "fwd|rev" strand + cdna|cds|noaa,fasta|notfa + full|longpart|bestpart
# .. split action into format=fasta|table + strand=fwd|rev|both + part=full|long|best
# .. -act cdnanoaafasta  or -act cdnaonlyfasta
#  my $whichstrand= ($action =~ /rev/)? "rev" : ($action =~ /fwd/) ? "fwd" : "both";
#  my $fullpart= ($action =~ /full/)?"full":($action =~ /long/)?"longpart":"bestpart";
#  my $outtype = ($action =~ /fasta/)?"fasta":"table";
my $noprotaction= ($action=~/cdna|cds/ and $action=~m/noaa|only/)?1:0;  # and $outtype eq "fasta"

# FIXed: add option to filter out huge-span asmrna with poor qual : no intron evidence of long intron; poor prot
## protein=MLSHQLLEDSTMMQMKHGLRQGRENICQGSRLLLIGNVLVDNXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX..
$MINGOOD= $MINGOOD/100.0 if($MINGOOD > 1); 

##was $ORF_FULLvPART= $ORF_FULLvPART/100.0 if($ORF_FULLvPART > 1); 
 # ^^ -full=0 means only full orfs?; -full=1 means never?? change to -full=1 means full only; -full=0 longest part always
$ORF_FULLvPART=1 unless($ORF_FULLvPART =~ /\d/);
$ORF_FULLvPART= $ORF_FULLvPART/100.0 if($ORF_FULLvPART >= 1); 
$ORF_FULLvPART=1 if($ORF_FULLvPART == 0); # means always longest

useSelenocysteine($USESelenocysteine) if($USESelenocysteine);

## for tr2aacds, check for file orthobest == cdnaseq=~s/.suffix/.orthobestorf/
unless($orthobest) {
  my $obf= $cdnaseq; $obf=~s/\.gz//; $obf=~s/\.\w+$//; $obf.=".orthobestorf";
  $orthobest= $obf if(-s $obf); 
}
if($orthobest) {
  # our $InnerStopToX= -1; # option!
  if($orthobest =~ /^(\w+)=(\d+\-\d+:.)/) { 
    # debug input:  -orthobest socatfishv1k21loc463048t1=109-1889:+
    my($tid,$tralnspan)=($1,$2);
    $orthobest{$tid}= join"\t",$tralnspan,"1-99","debug","1-33","ref:id"; #?    
  } else {
  open(OB,$orthobest) or die "ERR reading $orthobest";
  while(<OB>) {
    next if(/^\W/); my($rid,$tid,$naln,$refspan,$tralnspan,$claaspan)=split; #??
    $claaspan||="clmiss:1-99:.";
    my($cla,$aaspan)=split":",$claaspan,2;
    $orthobest{$tid}= join"\t",$tralnspan,$refspan,$cla,$aaspan,"ref:$rid"; #?
  } close(OB);
  }
}
#   orthobest: diff table like this?
#   head catfish1evg8tr-refkfish.diffaa.tab2 | cut -f1,2,5,8,9,10 | head
#   # RefID             TrID                       nAlign   refaaSpan  trAlignSpan  Class:aaOrfSpan
#   Funhe2EKm000007t1   socatfishv2k31loc62810t1    326     1-325      2341-1364:-  revutr:3589-7470:+  # utrorf
#   Funhe2EKm000108t1   socatfishv1k25loc49492t1    923     23-919     1-2706:+     over:22-2709:+   # small over, 20n

sub openout {
  my($oname,$osuf,$ohandle)= @_;
  my $OUTH= undef;
  if(defined $oname and not $oname) { # use cdnaseq name
    $oname= ($cdnaseq =~ /stdin|^\-$/)? "cdna_bestorf_out" : $cdnaseq;
    $oname=~s/\.gz//; $oname=~ s/\.\w+$//; $oname.= $osuf;
  }
  if($oname) { 
    $oname.= "out$osuf" if($oname eq $cdnaseq);
    open($OUTH, ">$oname") or die "ERR writing $oname"; $ohandle= $OUTH; # *OUTH; 
  }
  return($ohandle,$oname);
}

# outputs:
## drop use of "-action=xxxcds -action=xxxcdna for outputs
$cdsout="" if($action =~ /cds/i and not defined $cdsout); # define output
$cdnaout="" if($action =~ /cdna/i and not defined $cdnaout); # define output

my($OUTAA,$onameaa)  = openout($aaout,".aa",*STDOUT);
my($OUTCDS,$onamecds)= openout($cdsout,".cds",undef);
my($OUTCDNA,$onamecdna);
if(defined $mrnaout) { 
  ($OUTCDNA,$onamecdna)= openout($mrnaout,".mrna",undef);
} else {
  ($OUTCDNA,$onamecdna)= openout($cdnaout,".cdna",undef);
}

MAIN: { 
  my $ovh; my $ok= 0; 
  my($nin,$tnin,$tngood,$tnsplit,$tnskip)=(0) x 9;
  if($cdnaseq =~ /.gz$/) { $ok= open(OVR,"gunzip -c $cdnaseq |");  $ovh= *OVR; }
  elsif($cdnaseq =~ /stdin|^\-$/) { $ok=1; $ovh= *STDIN; }
  elsif($cdnaseq) { $ok= open(OVR,$cdnaseq); $ovh= *OVR; }
  die "bad -cdnaseq=$cdnaseq" unless($ok);
  my $id="none";
  while(<$ovh>) {
    if(/^>(\S+)/) { $id=$1; 
      do{ my($ngo,$nspl,$nsk)= cdna_bestorf(); $tngood+=$ngo;$tnsplit+=$nspl;$tnskip+=$nsk;
         %cdnaseq= %cdnahead=(); $nin=0; } if($nin>49);
      $cdnaseq{$id}="";  $nin++;  $tnin++;
      if(m/ +(\S.+)$/){ 
      	my $h=$1; $h =~ s/\s*\b(len|cf|nt)=\S+//g; 
      	$h =~ s/\baa c offs=/offs=/;	## fixup old mess: 'aa c offs=xxx' from bad (aa|c)len= cut
      	$cdnahead{$id}=$h; }
      ## FIXME: cut (len|cf|nt)= is bad? need option/spec of atts to drop, esp. replace new atts
    } elsif(/\w/) { chomp; $cdnaseq{$id}.= uc($_); } } 
  close($ovh);

  my($ngo,$nspl,$nsk)= cdna_bestorf(); 
  $tngood+=$ngo;$tnsplit+=$nspl;$tnskip+=$nsk;
  warn"#cdna_bestorf from cdnaseq:$cdnaseq n=$tnin, ngood=$tngood , nsplit=$tnsplit, nskip=$tnskip\n" 
     if $DEBUG;
  $onameaa||=""; $onamecds||=""; $onamecdna||=""; #nogood# $^W=0; # $WARNING=0;
  warn"#cdna_bestorf output to aa:$onameaa  cds:$onamecds cdna:$onamecdna \n" if $DEBUG;

}  

#..................

# ## this is slow.. not best way to do this. cd-hit-454 -c 1.00 will be faster
# my %containedin;
# sub cdna_containedin  


sub cdna_bestorf
{
  # FIXmaybe: add other option before output: skip subset cdna contained in *good* longer cdna
  # .. change sorting of @cdnain to long > short-containedin
  my @cdnaid; 
     @cdnaid= sort keys %cdnaseq; 
  
  my $saveSelcOn= $USESelenocysteine;
  my %isgoodseq;
  my $ncdna= @cdnaid;
  # warn"#cdna_bestorf from cdnaseq:$cdnaseq n=$ncdna\n" if $DEBUG;
  # 1505 subtle bug, only for -outmrna, changes strand=- to strand=+, new offs, mrna=revcomp(cdna)
  # .. so that outaa, outcds strand/offs need changing to match;
  # .. change here to put mRNA 1st, CDS, AA follow?
  my @outadd=(); 
  # push @outadd, "CDS" if(defined $OUTCDS); 
  if(defined $OUTCDNA) {
    if(defined $mrnaout) { push @outadd, "mRNA"; } else{ push @outadd, "cDNA"; }
  }
  push @outadd, "CDS" if(defined $OUTCDS); 
  push @outadd, "AA" if(defined $OUTAA); # FIX1505 

  my ($ngood,$nsplit,$nskip)= (0) x 10; 
  my $MINTR= $MINAA * 3;
  my @rev=($action =~ /rev/)? (1) : ($action =~ /fwd/)? (0) : (0,1);
  my $whichstrand= ($action =~ /rev/)? "rev" : ($action =~ /fwd/) ? "fwd" : "both";
  my $fullpart= ($action =~ /full/)?"full":($action =~ /long/)?"longpart":"bestpart";
  $fullpart= $whichstrand."strand".$fullpart;
  
# ADD maybe: detect joined/fused genes using 
#  1. cds span < ~ 1/2 trspan, leaving long utr, << annotate these cases, require aa complete/partial5
#  2. check the utr for long cds
 
  foreach my $id (@cdnaid) {
    my $cdnain= $cdnaseq{$id}; 
    my $cdnainsize= length($cdnain); 
    my $cdnabest= $cdnain;
    my $cdnasize= $cdnainsize; 
    my($orfprot, $prostart5, $proend3, $bestorf, $utrorf, $orflen, $strand, $isrev, $aalen, $pcds, $compl)= (0) x 20;
    my $doskip= 0;
    my($fahead,$orfprotfa,$revaalen,$revinfo)=("") x 9;
    my $bestproflags= $fullpart;
     
    if($cdnasize < $MINTR) {  # and not $orthobest{$id} ??
      $doskip=1; $nskip++; # warn??
    } else {

    #...... new : both dirs ...
    if($FIXUPMRNA) {
      # assume called as
      # cdna_bestorf.pl -fixmrna -nostop -outaa fixed.aa -outcds fixed.cds -outmrna fixed.mrna -cdna input.err.mrna
      # .. does this require same cds offs span as input.err.mrna ?
      # parse header for selcstop or revc=1;  'offs=1472-756:-;' is rare alias of revc=1 ??
      my $cdnah= $cdnahead{$id}||"";
      if($cdnah =~ /\bselcstop/) { useSelenocysteine(1); } else { useSelenocysteine($saveSelcOn); } # ($USESelenocysteine) ??
      if($cdnah =~ /\brevc=1/) { # input mrna is forward, but offs=rev strand=-; need to change or not?
        # set fwd only? require cds offs ?
        my $fullpart1= ($action =~ /full/)?"full":($action =~ /long/)?"longpart":"bestpart";
        $bestproflags="fwdstrand$fullpart1";
      }
    }
    # FIXed: respect action =~ rev or fwd here: fullpart fixed
    my $allorfs;
    
    ## orthobest: flag for getBestProt2 as /bestspan=([\d\:\+\-]+)/ 
    if($orthobest and $orthobest{$id}) { 
      my($tralnspan,$refaaspan,$cla,$aaspan,$rid)= split"\t",$orthobest{$id};
      # FIXME-not-here: spans here are 1origin, subs want 0origin
      # $bestproflags .=",frameck" if($XXXwhat);
      $bestproflags .=",bestspan=$tralnspan:$cla:$rid:$refaaspan";
      our $InnerStopToX= -1 unless($DEBUG); # option!
    }
    
    ( $bestorf, $allorfs)= getBestProt2( $bestproflags, $cdnain); # , undef, $oldStart_b,$oldStart_e
 
#    my $lbestorf=undef;
#    # orthobest: in addition to getBestProt2() ; merge w/ that so we have @starts,@stops array ??
#     if($orthobest and $orthobest{$id}) {  
#       $lbestorf= $bestorf;
#       my($tralnspan,$refaaspan,$cla,$aaspan,$rid)= split"\t",$orthobest{$id}; #? aaspan should be above lbestorf span
#       my($cbe,$co)=split":", $tralnspan; ## 2341-1364:-
#       my($cb,$ce)= split /[\-]/, $cbe; 
#       if($cb>$ce) { $co='-' unless($co eq '-'); ($cb,$ce)=($ce,$cb); }
#       my($rab,$rae)= split /[\-]/, $refaaspan; # 5-99 : aa-span, not cds, no orient
#       ## need to extend cb,ce by refspan inset .. or fix input orthobest table.
#       ## ** orthobest tab doesnt know where start/stop codons may be, nor frame of this cb unless refspan says.
#       $bestorf = getOneOrf($cb,$ce,$cdnain,$co);  
#     }
       
    # Ooops, orf->{start,stop} are reversed for -strand; start>stop : ($lb,$le)=($le,$lb) if($lb>$le);
    ( $orfprot, $prostart5, $proend3, $orflen, $strand)= orfParts($bestorf, [qw(protein start stop length orient)]);
    $isrev= ($strand eq '-') ? 1 : 0;
    #NO, always input strand now# $cdnabest= revcomp($cdnain) if($isrev); #?? maybe not, only for output here.
    
    ($utrorf)= utrorf_test( $bestorf, $allorfs, $cdnasize); # noutrorf option here or below?

    ## 201402 fixme for shorter bestspan  ; and $bestorf->{flags}=~/bestspan/  ??
    ## is this best choice? some bestorf are complete vs partial, longer utrorf
use constant FIXME201402utrorfbest => 0; ## 1505 CANCEL, getBestProt2() already has bestorf from goodlen, other quals.
    if(FIXME201402utrorfbest and $utrorf and ($utrorf->{goodlen} > $bestorf->{goodlen})) { 
      ($bestorf,$utrorf)=($utrorf,$bestorf); 
      ## FIXME.1505 ^^ annot bugs here from swap; need to swap orfParts()
      ( $orfprot, $prostart5, $proend3, $orflen, $strand)= orfParts($bestorf, [qw(protein start stop length orient)]);
      $isrev= ($strand eq '-') ? 1 : 0;
    }
    
		if($revaa_report) {
			($revaalen,$revinfo)= revorf_report( $bestorf, $allorfs, $cdnasize); # return info only for longest reversed orf
			# revinfo == "aarev=55%[rev2fwd],aa99,90%,complete,s-,o9-309"
		}

    } # MINTR
 
     
    #above# my($fahead,$orfprotfa,$revaalen,$revinfo)=("") x 9;
    if($orfprot) {
      if($bestorf->{mrnaseq}) { $cdnabest= $bestorf->{mrnaseq}; } # 201402: is this ok now?
      ($aalen,$pcds,$compl,$orflen,$fahead,$orfprotfa)
          = proteindoc($bestorf,$cdnasize,$strand,($action =~ /fasta/));

		  ## note: fahead= "aalen=$aalen,$pcds%,$compl; clen=$cdnalen; strand=$cdnarev; offs=$prostart-$proend;";
      #  if( $skipsubsetseq ) {}   #.. drop this, use other methods to remove subset seq
      $isgoodseq{$id}++ unless($doskip);
      $ngood++ unless($doskip); # unless($utrorf) ??

#      ## change this, use bestorf->flags to add to fahead, in proteindoc now *
# 			if($bestproflags=~/bestspan/) {  
# 			  # my $ns= $orfprot =~ tr/*/*/; # or $ns= $bestorf->{innerstop};
# 			  # $ns-- if(substr($orfprot,-1,1) eq '*'); # stopatend, or ($bestorf->{complete} & 2)
# 			  # $fahead.=" innerstop=$ns;" if($ns>0); # moved this to proteindoc()
#         $fahead.=" bestflag=$bestproflags;";
#       }
      
			$fahead.=" $revinfo;" if($revaa_report and $revaalen >= $MINAA and $revinfo);
				# revinfo == "aarev=55%[rev2fwd],aa99,90%,complete,s-,o9-309"
        ## revaa_report : do we want option to also print OUTAA >id.rev\nrevprotfa ??
        ## OR option to replace fwdaa with revaa output (aa and cds), IF %revaa/fwdaa >= minrevaa
    }
    
    if($doskip) {
      # warn? any output?  $nskip++ above
      
    } elsif($action =~ /fasta/) {  # FIXME, do fasta or table here, format at print...
      my $cdnah= $cdnahead{$id}||"";
			$cdnah =~ s/\s*\b(aalen|clen|strand|offs|type|revc)=[^;\s]+[;]?//g;
			$cdnah =~ s/\s*\b(path)=[^;\n]+[;]?//g; # trin crap
			## above fixup old mess: 'aa c offs=xxx'
			# 2013feb: add output files: aaout, cdsout, cdnaout?
      
      #x if($action =~ /all/ or $orfprot) #?? no longer useful test
      { 
        # FIXME.1505 fahead mismatch aa<>mrna after revc(mrna) below, 
        # .. move all out{aa,cds,mrna} to loop here w/ same fahead annots? move revc(mrna) above this
        # .. dont need noprotaction now?, using @outadd list for all
        use constant FIX1505outaa => 1;
        
        unless(FIX1505outaa or $noprotaction) {
        print $OUTAA ">$id $fahead $cdnah\n";
        print $OUTAA $orfprotfa,"\n";
        }
        ## revaa_report : do we want option to also print OUTAA >id.rev\nrevprotfa ??
        ## OR option to replace fwdaa with revaa output (aa and cds), IF %revaa/fwdaa >= minrevaa
        
        my $utrorfseq=""; my $utrcut=0;
        if($utrorf and not $noutrorf) {
            # main:  $orfprot, $prostart5, $proend3, $orflen, $strand
          my($uprostart5, $uproend3, $uorflen, $ustrand)= orfParts($utrorf, [qw( start stop length orient)]);
          my($aalen,$pcds,$compl,$orflen,$fahead,$orfprotfa)
              = proteindoc($utrorf,$cdnasize,$ustrand, 1);
          unless(FIX1505outaa or $noprotaction) {
            print $OUTAA ">",$id,"utrorf $fahead $cdnah\n";
            print $OUTAA $orfprotfa,"\n";
            }
            
          # FIXMEd for utrorf: find way to split cdna b/n 1st, 2nd orf: midway?
          unless($splitutrorfs) { 
            # part of splitutrorfs, need valid utrorfseq for calcs.. do that part
            if($prostart5 > $uprostart5) {
              my $uend= _max($uprostart5,$uproend3) + 9;
              $utrorfseq= substr($cdnain, 0, $uend);
            } else {
              my $ustart= _min($uprostart5,$uproend3) - 9;  
              $utrorfseq= substr($cdnain, $ustart - 1);
              ($uprostart5,$uproend3)= map{ $_ - ($ustart-1) } ($uprostart5,$uproend3);  
            }
            if($uprostart5 > $uproend3) {
              $utrorfseq= revcomp($utrorfseq);
              $ustrand='+';
              my $usize= length($utrorfseq);
              ($uprostart5,$uproend3)=($usize +1 - $uprostart5, $usize +1 - $uproend3);   
            }
           $utrorf->{start}= $uprostart5;
           $utrorf->{stop} = $uproend3;
           $utrorf->{orient} = $ustrand;
           # $utrorf->{mrnaseq}= $utrorfseq;#?
          } 
          
          if($splitutrorfs) {
            # Ooops, orf->{start,stop} are reversed for -strand; start>stop : ($lb,$le)=($le,$lb) if($lb>$le);
            # *** matters here if cdnabest= revcomp(cdnain) ***

            my $bestlast= ($prostart5 > $uprostart5) ? 1 : 0;
            my $revbest= ($prostart5 > $proend3)?1:0;
            my $revuorf= ($uprostart5 > $uproend3)?1:0; 
            my($b1,$e1)= ($revbest) ? ($proend3,$prostart5) : ($prostart5,$proend3);
            my($b2,$e2)= ($revuorf) ? ($uproend3,$uprostart5) : ($uprostart5,$uproend3);
            
            my $cut= ($b2 >= $e1) ? $e1 + int(($b2 - $e1)/2) 
                   : ($e2 <= $b1) ? $e2 + int(($b1 - $e2)/2) : 0;
                   
            if($cut>0) {
              $utrcut  = $cut; # THIS utrcut is flag for annot, add bestfirst/last
              my $seq1 = substr($cdnain,0,$cut); # $cdnabest $cdnain same here??
              my $seq2 = substr($cdnain, $cut);
              ($cdnabest, $utrorfseq)=  ($bestlast) ? ($seq2, $seq1) : ($seq1, $seq2);
              $cdnabest = revcomp($cdnabest) if($revbest); # probably should also swap $prostart5 > $proend3
              $utrorfseq= revcomp($utrorfseq) if($revuorf);
              
use constant SPLITUP18 => 1;  

=item   SPLITUP18 splitutrorfs

   seems to work now  
   NO: got many ** inner stops, getOneOrf() problem
   may be offby1 from utrcut?, bug for fwd+rev 2nd orfs
   okaa: >Susscr4EVm008595t2utrorf aalen=150,62%,partial3; uorfcut=1718/2441; clen=723; strand=+; offs=272-721; 
   ofby: >Susscr4EVm008595t2utrorf aalen=150,62%,partial;  uorfcut=1718/2441; clen=723; strand=+; offs=273-722;
  
    
   still some bug.. maybe -strand bug, offs=HI-low
    >> getOneOrf()  Assumes seq == revcomp(seq) if dir eq '-'
  grep -c '*' old  pig4321ew.tsaerr_pt1.aa:330  
  grep -c '*' new  pig4321ew.tsaerr_pt1.aa:200
  grep -c '*' revc pig4321ew.tsaerr_pt1.aa:187 .. whats left?
    .. may be revcomp + cutat interaction >> cut1 = cut -1 for revcomp? .. NO, need 
    revcomp_coord(start|stop, seqlength) == $seq_length - $coord + 1
    
>Susscrtridba1a_sBn1l1SRR6236889idbaidbtk97Loc3631utrorf  uorfcut=456/9918; clen=456; strand=-; offs=455-57; 
RRRRRRRRRQRWRRRRRRRRRGG*AVVV
>Susscrtrsoap1a_sBn2l1SRR1519321soapk73Loc22064t1utrorf  uorfcut=2636/3948; clen=1312; strand=-; offs=843-1; 
LYIRKKDFNGXXXXXXXVVI*GIG*FQRIFIPYC*HYLQYLNYMAGFISVFQTVPVFPFL
>Susscrtrsoap1a_sBn2l1SRR5027062soapk25Loc88192t1utrorf  uorfcut=2219/3235; clen=1016; strand=-; offs=785-3; 
LPGLPRLAEGFFQLVSVGCVALLQVLQLLIFFLVQDPQEIL*LWQAEGFPLSNFMFTKNK
    
=cut

if(SPLITUP18) {
              #notnow: my $cut1= $cut + 1; # YES +1 to fix offby1? NOT -1 .. 0-origin needed for getoneorf, skip cut+1
              
              #x my ($bb,$be,$bs)= ($revbest) ?  ($proend3,$prostart5,'-') : ($prostart5,$proend3,'+');
              #y my ($bb,$be,$bs)= ($revbest) ?  ($cdnainsize-$prostart5+1, $cdnainsize-$proend3+1,'-') : ($prostart5,$proend3,'+');
              #x if($bestlast) { my $bcut=($revbest) ? $cut : $cut;  ($bb,$be)=map{ $_ - $bcut } ($bb,$be); }
              #x if($revbest) { $bs='-'; }
              
              $cdnasize= length($cdnabest);
              my ($bb,$be,$bs)= ($prostart5-1,$proend3-1,'+'); #  -1 for 0-origin
              if($bestlast) { ($bb,$be)=map{ $_ - $cut } ($bb,$be); }
              if($revbest) { $bs='+'; ($bb,$be)=($cdnasize-($bb+1),$cdnasize-($be+1)); } # add revcomp_coord AFTER cut down bb,be
              
              $bestorf= getOneOrf( $bb, $be, $cdnabest, $bs); 
              # my $bstops= $bestorf->{innerstop}; 

              # remake bestorf vars .. revise to avoid this
              ( $orfprot, $prostart5, $proend3, $orflen, $strand)
                  = orfParts($bestorf, [qw(protein start stop length orient)]);
              ($aalen,$pcds,$compl,$orflen,$fahead,$orfprotfa)
                 = proteindoc($bestorf,$cdnasize,$strand,($action =~ /fasta/));
              
              #x my ($ub,$ue,$us)= ($revuorf) ? ($uproend3,$uprostart5,'-') : ($uprostart5,$uproend3,'+');
              #y my ($ub,$ue,$us)= ($revuorf) ? ($cdnainsize-$uprostart5+1, $cdnainsize-$uproend3+1,'-') : ($uprostart5,$uproend3,'+');
              #x unless($bestlast) { my $ucut=($revuorf) ? $cut : $cut; ($ub,$ue)=map{ $_ - $ucut } ($ub,$ue); }
              #z if($revuorf) { $us='-'; } # getOneOrf handles rev ??: revcomp(seq),revcomp_coord

              my $ucdnasize=length($utrorfseq);
              my ($ub,$ue,$us)= ($uprostart5-1,$uproend3-1,'+'); # -1 for 0-origin
              unless($bestlast) { ($ub,$ue)=map{ $_ - $cut } ($ub,$ue); }
              if($revuorf) { $us='+'; ($ub,$ue)=($ucdnasize - ($ub+1), $ucdnasize - ($ue+1)); } # THIS IS IT: +1 to undo 0origin .. add revcomp_coord AFTER cut down bb,be
              
              $utrorf= getOneOrf( $ub, $ue, $utrorfseq, $us); # uses 0-origin, not 1
              # my $ustops= $utrorf->{innerstop}; 
              
} else {              
              $bestorf->{mrnaseq} = $cdnabest; # FIXME move this splitutr to package..
              $utrorf->{mrnaseq}  = $utrorfseq; # FIXME move this splitutr to package..
}              
              ## bugs below with old cdnalen, new uorflen, and cds offsets and revcomp
              # $cdnasize=length($cdnabest); #?? change or not; need change prostart/end also
              ## was BAD cut: rev, cant get 651aa from 857cds << reverse cseq, orfseq
              }
            }
            
          }
  

        # FIXME: remove tag from ID, add as type=tag
        # FIXME.1505, move AA into this out loop so mRNA annot changes will match in AA
        foreach my $tag (@outadd) {
          ## 201402: look for bestorf->{mrnaseq} before using cdnabest/cdnain seq
          ## .. maybe always add orf->mrnaseq with proper revcomp and cuts of utrorf, introns, ..
          my $mrnaseq= $bestorf->{mrnaseq} || $cdnabest;
          #o my $seq=  ($tag eq "CDS") ? $bestorf->{sequence} : $mrnaseq; 
          #o my $outh= ($tag eq "CDS") ? $OUTCDS : ($tag =~ /cDNA|mRNA/) ? $OUTCDNA : undef;

          my $tagshow= ($tag eq "AA") ? "protein" : $tag;
          my $seq=  ($tag eq "AA") ? $orfprotfa : ($tag eq "CDS") ? $bestorf->{sequence} : $mrnaseq; 
          my $outh= ($tag eq "AA") ? $OUTAA : ($tag eq "CDS") ? $OUTCDS : ($tag =~ /cDNA|mRNA/) ? $OUTCDNA : undef;

          ## option: mRNA (+strand always) instead of cDNA output:
          next unless($seq =~ /\w/ and defined $outh); # report?
            
          ## 201402: replace cutadd with bestorf->mrnacuts/mrnasegments info ?
          ## 1505: change output clen= to cutlen? should match output mrna length.
          #o: my $cutadd=($utrcut>0)? " cutlen=".length($mrnaseq).",at=$utrcut/$cdnainsize;" : "";
          my $cutadd=($utrcut>0)? " uorfcut=$utrcut/$cdnainsize;" : "";
          # my $cutadd=""; if($utrcut<0) { $cutadd.=" uorfcut=1..$utrcut/$cdnainsize;" } 
          # elsif($utrcut>0) {  $cutadd.=" uorfcut=$utrcut..$cdnainsize;"  }
          
          ## FIXME.1505, want to match outaa, outcds annots (offs,strand) to outmrna AFTER revcomp changes
          ##  .. but only if -outmrna requested 
          ## also have cases of prostart5 > proend3 but strand == +, this seems to be utrorf / mainorf mixup
          ## ^^ fixed above, from swap(bestorf,utrorf) but missed redo orfParts()
          ## dang2: should add revcom(utrorfseq) for utrorf.mrna strand=-
          
          if($tag eq "mRNA" and $strand eq "-") {
            my $flipflag="";
            ($flipflag,$seq,$strand,$prostart5,$proend3)= 
              mRNAflip($seq,$strand,$prostart5,$proend3, $cdnasize, $bestorf);
            $cutadd.= $flipflag if($flipflag);
          }
          # if($tag eq "mRNA" and $strand eq "-") { ## do for mRNA, not for cDNA, others
          #   $seq= revcomp($seq); $cutadd.=" revc=1;";  # need extra flag? YES
          #   my $rbeg= revcomp_coord($prostart5, $cdnasize); $prostart5= $rbeg;
          #   my $rend= revcomp_coord($proend3, $cdnasize); $proend3= $rend;
          #   ## !:( also do revcomp on Selc positions..
          #   if($bestorf->{Selc}) { 
          #     my @ss= split",",$bestorf->{Selc}; 
          #     $bestorf->{Selc}= join",", map{ revcomp_coord($_, $cdnasize) } @ss; 
          #     }
          #   $strand="+";
          #   #^^ *?? should also change offs=b>e to offs=len-e<len-b (revcomploc); strand=- to =+ **??
          #   #?? but also should do on aa,cds ?? or not, want that info..
          #   #^^ dang, this is bad now for $bestorf->{mrnaseq}, probably :((
          # }
            
            ## Dang, need aa fahead info: Selcstop=nnn,nnn, other?
          if( my $isel= $bestorf->{Selc}){ $cutadd.=" Selcstop=$isel;"; }

          $seq  =~ s/(.{60})/$1\n/g unless($seq eq $orfprotfa);
          # append $cdnah to header ??
          my $taghead= "aalen=$aalen,$pcds%,$compl;$cutadd clen=$cdnasize; strand=$strand; offs=$prostart5-$proend3;";
          print $outh ">",$id, " type=$tagshow; $taghead $cdnah\n";
          print $outh $seq,"\n";
 
          if($utrorf and not $noutrorf) {
            my $umrnaseq= $utrorf->{mrnaseq} || $utrorfseq;
            my $ucdnasize=length($umrnaseq);
            #o my $uoseq=  ($tag eq "CDS") ? $utrorf->{sequence} : $umrnaseq;  
            ## dang2: should add revcom(utrorfseq) for utrorf.mrna strand=-
            my $uoseq=  ($tag eq "AA") ? $utrorf->{protein} : ($tag eq "CDS") ? $utrorf->{sequence} : $umrnaseq;  
            if($uoseq) {
              my( $uprostart5, $uproend3, $uorflen, $ustrand)= orfParts($utrorf, [qw( start stop length orient)]);
              #o my($aalen,$pcds,$compl) = proteindoc($utrorf,$cdnasize,$ustrand, 1);
              my($aalen,$pcds,$compl,$uorflen1,$ufahead,$uorfprotfa) = proteindoc($utrorf,$ucdnasize,$ustrand, 1);
              ## 201402: replace cutadd with utrorf->mrnacuts/mrnasegments info ?
              
              #o my $cutadd=($utrcut>0)? " cutlen=$uorflen,at=$utrcut/$cdnainsize;" : "";
              my $cutadd=($utrcut>0)? " uorfcut=$utrcut/$cdnainsize;" : "";
              
              ## -splitutrorf bug: offs=-3385--1949; after flip, due to prostart,end in cdnasize span, not uorflen span
              if($tag eq "mRNA" and $ustrand eq "-") {
                my $flipflag="";
                ($flipflag,$uoseq,$ustrand,$uprostart5,$uproend3)= 
                  mRNAflip($uoseq,$ustrand,$uprostart5,$uproend3, $ucdnasize, $utrorf); #WRONG for rev-off: $uorflen,
                $cutadd.= $flipflag;
              }
              
              if($tag eq "AA"){ $uoseq= $uorfprotfa; } else { $uoseq  =~ s/(.{60})/$1\n/g; }
              print $outh  ">",$id,"utrorf", " type=$tagshow; aalen=$aalen,$pcds%,$compl;$cutadd clen=$ucdnasize; strand=$ustrand; offs=$uprostart5-$uproend3; \n";
              print $outh $uoseq,"\n";
              $nsplit += 2; #not now# $ngood--;
            }
          }
        }
          
      }
      
    } else { # move this outformat above with >fasta

			if($revaa_report and $revaalen >= $MINAA and $revinfo) {
    		$compl.=":$revinfo" ; # HACK, not a good place, add Notes column before orfprot ?
    		# revinfo == "aarev=99,90%,complete,-,9-309"
			}
			
      print $OUTAA join("\t", $id, $aalen, $pcds, $compl, $cdnasize, $strand, $prostart5, $proend3, $orfprot),"\n"; 
      if($utrorf and not $noutrorf) {
        my( $orfprot, $prostart5, $proend3, $orflen, $strand)= orfParts($utrorf, [qw(protein start stop length orient)]);
        my($aalen,$pcds,$compl) = proteindoc($utrorf,$cdnasize,$strand, 0);
        print $OUTAA join("\t", $id."utrorf", $aalen, $pcds, $compl, $cdnasize, $strand, $prostart5, $proend3, $orfprot),"\n"; 
      }
    }

  }
  # warn"#cdna_bestorf found $ngood good, $nsplit split, $nskip skip of $ncdna\n" if $DEBUG; 
  # fixme, write split count ngood
  return($ngood,$nsplit,$nskip);
}

sub mRNAflip {
  my($seq,$strand,$prostart5,$proend3, $cdnasize, $bestorf)= @_;
  my $flipflag="";
  if($strand eq "-") { ## do for mRNA, not for cDNA, others
    $seq= revcomp($seq); $flipflag=" revc=1;";  # need extra flag? YES
    my $rbeg= revcomp_coord($prostart5, $cdnasize); $prostart5= $rbeg;
    my $rend= revcomp_coord($proend3, $cdnasize); $proend3= $rend;
    ## !:( also do revcomp on Selc positions..
    if($bestorf->{Selc}) { 
      my @ss= split",",$bestorf->{Selc}; 
      $bestorf->{Selc}= join",", map{ revcomp_coord($_, $cdnasize) } @ss; 
      }
    $strand="+";
    #^^ *?? should also change offs=b>e to offs=len-e<len-b (revcomploc); strand=- to =+ **??
    #?? but also should do on aa,cds ?? or not, want that info..
    #^^ dang, this is bad now for $bestorf->{mrnaseq}, probably :((
  }
  return($flipflag,$seq,$strand,$prostart5,$proend3, $cdnasize, $bestorf);  
}

__END__

### MOVED to cdna_proteins.pm ........................................................
# convert to array handling? _max(@list) 
# sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
# sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }
#
# sub getBestProt {}
# 
# ## revise longest_orf_finder
# 
# sub getAllOrfs {}
# 
# sub get_orfs {}
# 
# sub identify_putative_starts {}
# 
# sub identify_putative_stops {}
# 
# sub revcomp {}
# 
# sub revcomp_coord {}
# 
# # parts from PASA/PasaLib/Nuc_translater.pm 
# use vars qw ($currentCode %codon_table);
# 
# # sub codon1 {}
# 
# sub translate_sequence {}
# 
# BEGIN {}
# 
# 1;
