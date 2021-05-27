#!/usr/bin/env perl
# traa2cds.pl

=item about
  
  extract coding seq from transcripts, using offset from protein.aa headers
  as produced by cdna_bestorf.
  
  traa2cds.pl -cdna catfish1best2.tr.gz -aa catfish1best2.aa.gz > catfish1best2.cds
  
  eg.  offs=start-end (before revcomp, strand=-)
  aa>socatfishv1k25loc147014t1 aalen=238,99%,partial; clen=718; strand=-; offs=716-3; 
  tr>socatfishv1k39loc386507t1 
  CCTCCTTTCGCTCTCCTTCAGCCACGACGAAACTCTCTCTCTTGTTCGCCGAGTGGCACAAGGGGCTTTT
  cds> revcomp( substring(tr,3-1,1+716-3) )

=item utrorf update 2018.06
  
  traa2cds -utrorf -mrnaout : will now generate full mrna seq of utrorf(cdna) AND mainorf(cdna)
    ie. two+ outputs from one input cdna seq, depending on aa seq IDs ( ID1 & ID1utrorf )
    .. ?? make default -utrorf for -mrnaout, -noutrorf to cancel
    
  traa2cds -utrorf -mrnaout -cdna catfish.tr -aa catfish.aa  > catfish.mrna
  
=item author
  
  don gilbert, gilbertd near indiana edu, 2012
  part of EvidentialGene, evigene/scripts/

=cut

use strict;
use Getopt::Long;

my (%argh);
my $optok= GetOptions( \%argh,
  "cdnaseq|trseq|input=s", ## @ARGV also see below
  "aaseq=s", 
  "output|cdsseq:s", #? make default -out=infile.cds unless -out stdout
  "logfile:s", 
  "mrnaout|trout!", #? combine w/ output|outcdna opt ?
  "utrorf!", "nomiss", 
  "debug!", 
  );

## FIXME here?: add utrorf options, need cdna_bestorf -splitutrorfs
##  for input -aa only >IDutrorf  and -cdna all.cdna, 
##   outputs .cds, .mrna for just -aa IDutrorf with mrna as -splitutrorf
my $cdnaseq= $argh{'cdnaseq'} || '';
my $aaseq= $argh{'aaseq'} || '';
my $debug= $argh{'debug'}  || 0;
my $TROUT= $argh{'trout'} || 0; # mrnaout alias
my $UTRORF= $argh{'utrorf'}  || 0;
my $MISSNOLOG= $argh{'nomiss'}  || $UTRORF; # or what?

die "usage: traa2cds.pl -cdna catfish1best2.tr.gz -aa catfish1best2.aa.gz [-out xxx.cds] > catfish1best2.cds
  opts: -trout : cdna output, not cds (revcomp input.tr to match aa strand)  -log : logfile  -debug
" unless($optok and $cdnaseq and $aaseq);

my( %vec, %gene, %nam, %geneinfo, %idput);
# my $GAPSMAX = ('N') x $MAXGAP;

## want append-out option ?
my $outh= *STDOUT;
my $output= $argh{'output'}; #? make default -out=infile.cds unless -out stdout
unless($output) { # use cdnaseq name; NOT  and exists($argh{'output'})
  my $osuf=($TROUT) ? ".mrna.tr" : ".cds";
  $output= makename( $output||$cdnaseq, $osuf);  
}
if($output and $output!~/stdout|^-/) { open(OUT, ">$output") or die $output; $outh= *OUT; }

my $logh= undef;
sub loggit{ my $dowarn=shift; my $s= join(' ',@_); chomp($s);
  if($logh){ print $logh "#ta2c: $s\n"; } elsif($dowarn||$debug){ warn "#ta2c: $s\n"; }}
# sub loggit { my $dowarn=shift; if($logh) { print $logh @_; } elsif($dowarn) { warn  @_; } }

my $logfile= $argh{'logfile'};
if(not $logfile and exists($argh{'logfile'})) { # use output name
  $logfile= makename( $output||$cdnaseq, ".traa2cds.log");  
}
if($logfile) { open(LOG, ">$logfile") or die $logfile; $logh= *LOG; }


# my $tblh= undef; # genbank submit annot
# if($tblfile) { open(TBL, ">$tblfile") or die $tblfile; $tblh= *TBL; }

  
MAIN: {
loggit(0, "out=$output from $cdnaseq, $aaseq");

my($fa, $hd, $fh)=("","",undef);
if($aaseq =~ /stdin|^-/) { $fh= *STDIN; }
elsif($aaseq =~ /\.gz/) { open(F,"gunzip -c $aaseq|") or die $aaseq; $fh=*F; }
else {open(F,$aaseq) or die $aaseq; $fh=*F; }
while(<$fh>) { 
  if(/^>(\S+)/) {  # only header from aa
    my $oid=$1;  
    #>socatfishv1k25loc147014t1 aalen=238,99%,partial; clen=718; strand=-; offs=716-3; 
    my($mapqual,$aaqual,$trlen,$cdsor,$cdsoff);
    unless( ($cdsoff)= m/offs=([^;\s]+)/ ) {  ## FAIL if missing cdsoff
      die "ERR: $oid missing cds offset info, off=$cdsoff, in $aaseq\n#ERRline:$_";
    }
    ($cdsor)= (m/strand=(.)/)?$1:".";
    ($aaqual)= (m/aalen=([^;\s]+)/)?$1:"";
    ($trlen)= (m/clen=([^;\s]+)/)?$1:0;
    $mapqual="na";
    
    ###  my($td,$gd,$mapqual,$aaqual,$trlen,$cdsor,$cdsoff)= @$rinfo; # ,$lotag,$nam1,$dbxref
	  #	trid	geneid	mapqual	aaqual	trlen	or	cdsoff	locustag	product	dbxref
    if( $geneinfo{$oid} ) {
      loggit(1, "ERR:dup id $oid in $aaseq"); next;
    }
    $geneinfo{$oid}= [ $oid, $oid, $mapqual, $aaqual, $trlen, $cdsor, $cdsoff ];
  }
  ## elsif(/^\w/) { chomp; $fa.=$_; }
} close($fh);


($fa, $hd, $fh)=("","",undef); 
my ($ok,$bpout,$nout)=(0,0);

## allow list of cdnaseq here? for utrorf: -aa name.utrorf.aa -cdna name.{okay,okalt,drop}.tr.gz

foreach my $cdna1 ($cdnaseq,@ARGV) {
  if($cdna1 =~ /stdin|^-/) { $fh= *STDIN; }
  elsif($cdna1 =~ /\.gz/) { open(F,"gunzip -c $cdna1|") or die $cdna1; $fh=*F; }
  else { open(F,$cdna1) or die $cdna1; $fh=*F; }
  my $okid=0;
  while(<$fh>) { 
    if(/^>(\S+)/) { my $d=$1; if($fa and $hd){ $nout++ if putseq($hd,$fa); }
      ## BAD, no utrorf on id # $hd=($UTRORF and not $geneinfo{$d})?"":$d;  ##  $oidfix=$oid.'utrorf';
      if($UTRORF) { $d="" unless($geneinfo{$d.'utrorf'} or $geneinfo{$d}); }
      $hd=$d; $fa=""; }
    elsif(/^\w/ and $hd) { chomp; $fa.= $_; }
  } close($fh);
  if($fa and $hd){ $nout++ if putseq($hd,$fa); } # last
}

loggit(0, "nout=",$nout);

## FIXME: check for aa.oid not in cdna.putseq >> utrorf set, others? .. list all idmiss?
## Fixme2: skip this if caller knows , or give extra.idfile?
my @idmiss= sort grep { not $idput{$_} } keys %geneinfo;
if(@idmiss) { ## MISSNOLOG here ??
  my $nmiss=@idmiss;  ##my @some=grep /\w/, (@idmiss[0..4],@idmiss[-5..-1]);
  loggit(1,"MISSaa: nmiss=$nmiss");## , eg=",@some
  unless($MISSNOLOG) {
  for(my $i=0; $i<$nmiss; $i+=20) { # may be 100s to 1000s
    my $j=$i+20; $j=$nmiss-1 if($j>=$nmiss); my @some= @idmiss[$i..$j];
    loggit(1,"MISSaa: ids=",@some);
  } }
}

}


#-----------------------------

sub makename
{
  my($infile,$osuf,$insuf)=@_;
  $insuf ||= 'aa|cds|tr|fasta|fa';  ## fixme need insuf: tr|fasta|fa
  ( my $outfile= $infile ) =~ s,\.($insuf)[^\/\s]*$,,; 
  $outfile.= $osuf if($osuf); 
  $outfile.= "_out" if($outfile eq $infile);
  return $outfile;
}



sub cleanid { 
  local $_= $_[0];
#   ## project specific cleaning : read from where? cmdline or properties?
#   s/^cacao3(vel|v)(\d+)/vel$2/; # do we still need this?
#   s/^caca11r39cuf(\d+)/cuf$1/;  
#   s/^(B|L|P1|P2)/nwb$1/;
#   if($IDPREFIX and !/^$IDPREFIX/) {
#     s/^cacao[345]//; # FIXME; option?
#     $_= $IDPREFIX.$_;
#   }
  return $_;
}


sub revcomp {
  my ($seq) = @_;
  my $reversed_seq = reverse ($seq);
  $reversed_seq =~ tr/ACGTacgtyrkm/TGCAtgcarymk/;
  return ($reversed_seq);
}

  
sub putseq 
{
  my ($oid, $fa)=@_; 
  my $id= $oid; # cleanid($oid);
  my($ncut, $vw)=(0,0);
  my $falen= length($fa);
  
  my($vec,$gn,$nam,$rinfo,$cdsb,$cdse,$tlog);
  $tlog="";  ($cdsb,$cdse)=(0,0);

  #  if($GENEINFO_VERS == 2)  
  #	trid	geneid	mapqual	aaqual	trlen	or	cdsoff	locustag	product	dbxref
  ## do revcomp and rev-offset for cdsor=-
  my $oidfix= $oid;
  $rinfo= $geneinfo{$oid};
  
  ## FIXME: utrorf need to do this also if have oid, but also oid.utrorf .. make special utrorf.aa input subset
  unless($rinfo) {
    $oidfix=$oid.'utrorf';
    $rinfo= $geneinfo{$oidfix};
    $id=$oidfix if($rinfo);
  }
  $rinfo ||= [];
  my($td,$gd,$mapqual,$aaqual,$trlen,$cdsor,$cdsoff)= @$rinfo; # ,$lotag,$nam1,$dbxref
  $trlen= $falen; # use what we have as input
  ## SKIP cds output unless cdsoff
  unless($cdsoff) { loggit(1,"MISS: $oid\t$falen,$aaqual") unless($MISSNOLOG); return 0; } # option nolog here
  
  ## MISS err neg cdsoff:  offs=-6355--800
  ## .. bug in cdna_bestorf.pl? or here? no DATA bug, got dup ids, diff seqs, in .aa,.cds and .tr *!*!*!*!*!
  ## .. all from one trasm set mixup has same id prefix: daphmag/rnas/soap5xun3set soap5xinb3_co, soap5xinb3_cop2, 
  ## MISS: dmag5xun3cosoapk21loc758t5 type=cdna; aalen=1851,79%,complete; clen=592;  strand=+; offs=-6355--800;   revcds:-6355--800,
  ## MISS: dmag5xun3cosoapk21loc2222t1 type=cdna; aalen=304,44%,complete-utrbad; clen=455;  strand=+; offs=-1530--616;    revcds:-1530--616,
  #bad >dmag5xun3cosoapk21loc2195t2 type=cdna; aalen=569,81%,partial3; clen=618;  strand=+; offs=-1088-618;  revcds:-1088-618,
  #bad >dmag5xun3cosoapk29loc38286t1 type=cdna; aalen=696,60%,complete; clen=3467;  strand=+; offs=-2261--171;  

  ($cdsb,$cdse)= split/[-,]/,$cdsoff; # this SHOULD be offset in non-rev tr; any problems?

  if($UTRORF and $oidfix =~ /utrorf/) {
    # see cdna_bestorf.pl splitutrorfs : need to split fa == cdna
    my($cbr,$cer)= ($cdsb > $cdse) ? ($cdse,$cdsb) : ($cdsb,$cdse);
    my $u5= $cbr; my $u3= 1+$trlen - $cer;
    my($cut,$cutseq,$cutfl);
      ## also adjust cdsb,cdse
    if($u5 >= $u3) { # longer utr is one to remove
      $cut=$cbr-3; $cut=1 if($cut<1); $cutseq=substr($fa,$cut-1);
      $cbr -= $cut-1; $cer -= $cut-1; $cutfl="u5";
      ($cdsb,$cdse)= ($cdsb > $cdse) ? ($cer,$cbr) : ($cbr,$cer); # update cds
    } else {  
      $cut=$cer+3; $cut=$trlen if($cut>$trlen); $cutseq=substr($fa,0,$cut);  $cutfl="u3";
      ## no change to cbr,cer as cut is above them
    }
    $tlog.="utrorfcut:$cut.$cutfl,";  
    $fa= $cutseq;
    $falen= length($fa);  $trlen= $falen;  
    $cdsoff="$cdsb-$cdse"; # update !!!
  }
      
  if($cdsor eq "-") {
    my $farev= revcomp($fa); 
    $fa= $farev; #? save orig
    ## ** FIXME dammmit; if cb > ce, then this is offset BEFORE rev; if cb < ce is offset AFTER rev **
    my($cbr,$cer)=(0,0);
    if($cdsb > $cdse) {  ## this is new vers
      $cbr= 1+$trlen - $cdsb;
      $cer= 1+$trlen - $cdse;
    } else { # DAMMIT mess w/ cdna_org version changed offs= meaning; this is old vers
      $cer= $cdse; ## 1+$trlen - $cdsb;
      $cbr= $cdsb; ## 1+$trlen - $cdse;
    }
    $tlog.="revcds_of:$cdsoff:$cdsor,";  
    ($cdsb,$cdse)=($cbr,$cer);    
    $cdsoff="$cdsb-$cdse"; # update !!!
    $rinfo->[6]= $cdsoff;
    $rinfo->[5]= $cdsor= "+";
  }
  
  $cdsb=1 if($cdsb < 1);
  
  my($typ,$cdsfa,$len,$atlen)=(0) x 10;
  if($TROUT) {  # cdnaout
    $typ="mRNA"; # was "cdna"; 
    $cdsfa=$fa; # either same or revcom(input)
    $len= length($cdsfa);
    $atlen="clen=$len; ";
  } else { # CDSOUT
    $typ="CDS"; # was 'cds'
    $cdsfa= substr($fa, $cdsb-1, 1 + $cdse - $cdsb);
    $len= length($cdsfa);
    $atlen= "cdslen=$len; clen=$trlen; ";
  }
  # my $cdsfa= substr($fa, $cdsb-1, 1 + $cdse - $cdsb);
  #keep#>socatfishv1k25loc147014t1 aalen=238,99%,partial; clen=718; strand=-; offs=716-3; 
  my $def= "$id type=$typ; aalen=$aaqual; $atlen strand=$cdsor; offs=$cdsoff;";
  if($cdse < 1 or $len<1) { loggit(1,"MISS: $def  $tlog"); return 0; }

  $cdsfa =~ s/(.{60})/$1\n/g; 
  print $outh ">$def\n$cdsfa\n";  $idput{$oidfix}++;
  loggit(0,">$def  $tlog"); 
  return $len;
}


__END__

=item evgrutrorf.sh

  #!/bin/bash
  ## evgrutrorf.sh  : fix make okayset/*.utrorf.mrna to go w/ utrorf.aa,cds
  
  evigene=/bio/bio-grid/mb/evigene/
  
  ptnames=`/bin/ls -1 {banana,catfish,litova,locust,pogonus,shrimp,whitefly,ztick}*/*.names`
  # ptnames=`/bin/ls -1 locust*/*.names`
  ptdirs=`echo $ptnames | sed 's,/.*,,g;'`
  thisdir=`pwd`
  
  for ptn in $ptnames; do {
   ptd=`echo $ptn | sed 's,/.*,,g;'`
   pt=`basename $ptn .names`
   ptar=outz/$ptd.tar
   if [ ! -f $ptar ]; then 
     if [ -f outz/$ptd.tar.gz ]; then ptar=outz/$ptd.tar.gz; fi
   fi
  
   if [ -f $ptd/okayset/$pt.utrorf.mrna ]; then continue; fi
   echo "# make $ptd/okayset/$pt.utrorf.mrna";
  
   if [ ! -f $ptd/dropset/$pt.drop.tr.gz ]; then
    gtar -xvf $ptar $ptd/dropset/$pt.drop.tr.gz
   fi
  
   cd $ptd/okayset/
   gunzip -c $pt*{okay,okalt}.aa.gz | perl -ne'if(/^>/) { $ok=(m/utrorf /)?1:0; } print if($ok);' > $pt.utrorf.aain
  
   $evigene/scripts/prot/traa2cds.pl -utrorf -trout -aa $pt.utrorf.aain -out $pt.utrorf.mrna  -log \
       -cdna $pt.{okay,okalt}.tr.gz ../dropset/$pt.drop.tr.gz
  
     ## and remake .aa,.cds for hopefully small changes..
   $evigene/scripts/cdna_bestorf.pl -nostop -noutrorf -act fwdfasta -cdna $pt.utrorf.mrna \
        -outaa $pt.utrorf.aa -outcds $pt.utrorf.cds
    # /bin/rm $pt.utrorf.aain # old
  
    cd $thisdir
  } 
  done

=cut

#... old parts from  evigene/scripts/rnaseq/asmrna2ncbitsa.pl
#  my $ol= length($fa); 
#  my $nn= $fa =~ tr/Nn/Nn/; 
#  ## this is for ncbi submit, gap handling...
#   if($nn>0) { 
#     # ** FIXME: cdsb,cdse adjust for inner gaps
#     # fixme2: must adjust cdsb,e when cut BEFORE cdsb
#     # fixme3: this is a mess; better to a. cut NNN, b. rerun cdna_bestorf for new cds offset?
#     my ($lcdsb,$lcdse)= ($cdsb,$cdse);
#     $fa=~s/n/N/g;
#     ## SEQ_INST.HighNContentStretch: stretch of at least 5 Ns within the last 10 bases
#     my $ne= rindex($fa,'N'); if($ne >= $ol - $ENDGAP) { 
#       $fa=substr($fa,0,$ne); if($fa=~s/(N+)$//) {  my $ncut=length($1); $ne-=$ncut; }
#       if($ne < $cdse) { $cdse = $ne; } #??
#       }
#     my $n1= index($fa,'N'); if($n1 <= $ENDGAP) { 
#       $n1++; $fa= substr($fa,$n1);  if($fa=~s/^(N+)//) { my $ncut=length($1); $n1+=$ncut; }
#       if($cdsb>0) { $cdsb -= $n1; $cdse -= $n1; }
#       }
#       
#     $ncut=0; my $gapw= length( $GAPSMAX); #== MAXGAP
#     unless($GAPSOK) {
#     for (my $in= index($fa,$GAPSMAX); $in >= 0; ) {
#       my $w=length($fa); my $en=$in+$gapw; 
#       $en++ while($en<$w and substr($fa,$en,1) eq "N"); 
#       my $wn= $en-$in; my $keep= 3 + ($wn % 3); my $cut= $wn-$keep; $ncut+=$cut; 
#       my $facut= substr($fa,0,$in).substr("NNNNNN",0,$keep).substr($fa,$en); 
#       $fa=$facut; 
#       if($cdse>0) {
#         if($en < $cdsb) { $cdsb -= $cut; $cdse -= $cut; }
#         elsif($in < $cdse and $en > $cdsb) {
#           if($in <= $cdsb) { $cdsb -= $cut; $cdse -= $cut; } #??
#           else { $cdse -= $cut; }
#         }
#       }
#       $in=index($fa,$GAPSMAX); 
#     } 
#     }
#     
#     unless($cdse == $lcdse and $cdsb == $lcdsb) {
#       $rinfo->[6]= "$cdsb-$cdse"; # update !!!
#       $tlog.="cutcds=$cdsb-$cdse,oldcds=$lcdsb-$lcdse,"; 
#     }
#   } 
  
#   my $nl= length($fa); 
#   my $nn1= $fa=~tr/N/N/; 
#   my $len=($nl==$ol and $nn==0)?$nl:"$nl; olen=$ol; cut=$ncut; nnn=$nn1/$nn;"; 

#   $tlog.="vectrim=$vw," if($vw); 
#   if($nl<$MINSIZE) {
#     print $logh "#$def\toid=$oid\tgene=$gn\tlen=$len\t$tlog\tERR: too short:$nl\n" if($logfile);
#     return "ERR: too short:$nl" 
#   }
#
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
#         print $tblh "\t\t\tproduct\t$nam\n" if($nam);
#         # ^^ FIXME: SEQ_FEAT.MissingCDSproduct; need some nam: hypothetical??? or orig-name + similarto/
#         if($dbxref =~ /\w/) {
#         foreach my $dx (split",",$dbxref) { print $tblh "\t\t\tdb_xref\t$dx\n" if($dx=~/\w/ and $dx ne "na"); }
#         }
#         print $tblh "\n";
#       }
#     }
#   }
#   
#   $fa =~ s/(.{60})/$1\n/g; 
#   # defline option? strip  oid=..; len=.. from fahead before tbl2ans
#   # optional 2nd outfile for  newid, oldid, len, changes..
#   if($logfile) {
#     print $logh ">$def\toid=$oid\tgene=$gn\tlen=$len\t$tlog\n"; #??
#     print $outh ">$def\n$fa\n";
#   } else {
#     if($debug) { print $outh ">$def oid=$oid; len=$len\n$fa\n"; }
#     else { print $outh ">$def\n$fa\n"; }
#   }
#   return 0;
 


