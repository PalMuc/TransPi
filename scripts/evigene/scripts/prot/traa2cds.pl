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
  
  * update2: allow/need multiple inputs to resolv utrorf, as aa/cds utrorf may be classed
    in other file from cdsorf
   -aa okay.aa -aa okalt.aa -cdna okay.tr -cdna okalt.tr -out okall.mrna

=item new usage

  .. need all cdna inputs to resolve utrorfs: -aa okay.aa -aa okalt.aa -cdna okay.tr -cdna okalt.tr -out okall.mrna
  .. but can do one -aa at a time, only aa.ids output
  
  $evigene/scripts/prot/traa2cds.pl -mrnaout -aa $pt.okay.aa -aa $pt.okalt.aa -cdna $pt.okay.cdna -cdna $pt.okalt.cdna -log -out $pt.okall.mrna
  
=item author
  
  don gilbert, gilbertd near indiana edu, 2012
  part of EvidentialGene, evigene/scripts/

=cut

use strict;
use Getopt::Long;
use constant VERSION => '2018.06.18'; 
use constant DEFAULT_UTRORF => 1; 

# my (%argh);
# my $optok= GetOptions( \%argh,
#   "cdnaseq|trseq|input=s", ## @ARGV also see below
#   "aaseq=s", 
#   "output|cdsseq:s", #? make default -out=infile.cds unless -out stdout
#   "logfile:s", 
#   "trout|mrnaout!", #? combine w/ output|outcdna opt ?
#   "utrorf!", "nomiss", 
#   "debug!", 
#   );
# 
# ## FIXME here?: add utrorf options, need cdna_bestorf -splitutrorfs
# ##  for input -aa only >IDutrorf  and -cdna all.cdna, 
# ##   outputs .cds, .mrna for just -aa IDutrorf with mrna as -splitutrorf
# my $cdnaseq= $argh{'cdnaseq'} || '';
# my $aaseq= $argh{'aaseq'} || '';
# my $debug= $argh{'debug'}  || 0;
# my $output= $argh{'output'}; #? make default -out=infile.cds unless -out stdout
# my $logfile= $argh{'logfile'};
# my $TROUT= $argh{'trout'}  ||0; # mrnaout alias || $argh{'mrnaout'}
# my $UTRORF= $argh{'utrorf'}  || 0;
# my $MISSNOLOG= $argh{'nomiss'}  || $UTRORF; # or what?
# if(DEFAULT_UTRORF) { # should be default for both mRNA & CDS out
#   $UTRORF=1 unless(defined $argh{'utrorf'}); # UPD -noutrorf cancels default
# }

#upd1806 options:
my(@cdna,@aaseq);
my($output,$logfile)=(undef,undef);
my($TROUT,$UTRORF,$MISSNOLOG,$debug)=(0) x 9;
$UTRORF= undef;

my $optok= GetOptions( 
  "cdnaseq|trseq|input=s", \@cdna, 
  "aaseq=s", \@aaseq, 
  "output|cdsseq:s", \$output, #? make default -out=infile.cds unless -out stdout
  "logfile:s", \$logfile,
  "trout|mrnaout!", \$TROUT, #? combine w/ output|outcdna opt ?
  "utrorf!", \$UTRORF,
  "nomiss", \$MISSNOLOG, 
  "debug!", \$debug,
  );

my $cdnaseq=$cdna[0];
my $aaseq= $aaseq[0];
if(DEFAULT_UTRORF) { # should be default for both mRNA & CDS out
  $UTRORF=1 unless(defined $UTRORF); # UPD -noutrorf cancels default
}

die "EvidentialGene traa2cds VERSION ",VERSION,
"\ncreate CDS or mRNA sequences from input cDNA & protein (aa) sequences.
Usage: traa2cds.pl -cdna catfish.tr -aa catfish.aa [-out xxx.cds] > catfish.cds
opts : -mrnaout, mRNA output instead of CDS (revcomp input.cdna to match aa strand); -log, log to file\n"
  unless($optok and $cdnaseq and $aaseq);

my( %geneinfo, %idput);
# my @cdna=($cdnaseq); 
# push @cdna, @ARGV if(@ARGV > 0); #? want @ARGS here, need equiv aa?

## want append-out option ?
my $outh= *STDOUT;
unless($output) { # use cdnaseq name; NOT  and exists($argh{'output'})
  #o# my $osuf=($TROUT) ? ".mrna.tr" : ".cds";
  my $osuf=".cds";  
  if($TROUT) { $osuf=($cdnaseq=~/\.mrna/)?".outmrna":".mrna";  } 
  $output= makename( $output||$cdnaseq, $osuf);  
}
if($output and $output!~/stdout|^-/) { 
  if( -s $output) { rename($output,"$output.old"); } #preserve old?
  open(OUT, ">$output") or die $output; $outh= *OUT; 
}

my $logh= undef;
sub loggit{ my $dowarn=shift; my $s= join(' ',@_); chomp($s);
  if($logh){ print $logh "#ta2c: $s\n"; } elsif($dowarn||$debug){ warn "#ta2c: $s\n"; }}

if(defined $logfile and not $logfile) { #  and exists($argh{'logfile'}); use output name
  $logfile= makename( $output||$cdnaseq, ".traa2cds.log");  
}
if($logfile) { open(LOG, ">$logfile") or die $logfile; $logh= *LOG; }
  
MAIN: {
loggit(0, "out=$output from $cdnaseq, $aaseq");

my($fa, $fid, $hd, $fh)=("","","",undef);

foreach $aaseq (@aaseq) {
  if($aaseq =~ /stdin|^-/) { $fh= *STDIN; }
  elsif($aaseq =~ /\.gz/) { open(F,"gunzip -c $aaseq|") or die $aaseq; $fh=*F; }
  else { open(F,$aaseq) or die $aaseq; $fh=*F; }
  while(<$fh>) { 
    if(/^>(\S+)/) {  # only header from aa
      my $oid=$1;  
      #>socatfishv1k25loc147014t1 aalen=238,99%,partial; clen=718; strand=-; offs=716-3; 
      my($mapqual,$aaqual,$trlen,$cdsor,$cdsoff);
      unless( ($cdsoff)= m/(?:offs|cdsoff)=([^;\s]+)/ ) {  ## FAIL if missing cdsoff; any aliases? cdsoff=
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
}


($fa, $fid, $hd, $fh)=("","","",undef); 
my ($ok,$bpout,$nout)=(0,0);

## allow list of cdnaseq here? for utrorf: -aa name.utrorf.aa -cdna name.{okay,okalt,drop}.tr.gz

foreach my $cdna1 (@cdna) { #was ($cdnaseq,@ARGV)
  if($cdna1 =~ /stdin|^-/) { $fh= *STDIN; }
  elsif($cdna1 =~ /\.gz/) { open(F,"gunzip -c $cdna1|") or die $cdna1; $fh=*F; }
  else { open(F,$cdna1) or die $cdna1; $fh=*F; }
  my $okid=0;
  while(<$fh>) { 
    if(/^>(\S+)/) { my $d=$1; if($fa and $fid){ $nout += putseq($fid,$hd,$fa); }
      ## BAD, no utrorf on id # $hd=($UTRORF and not $geneinfo{$d})?"":$d;  ##  $oidfix=$oid.'utrorf';
      if($UTRORF) { $d="" unless($geneinfo{$d.'utrorf'} or $geneinfo{$d}); }
      $fid=$d; $hd=$_; $fa=""; 
      }
    elsif(/^\w/ and $hd) { chomp; $fa.= $_; }
  } close($fh);
  if($fa and $fid){ $nout += putseq($fid,$hd,$fa); } # last
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
  $insuf ||= 'aa|cds|cdna|tr|fasta|fa';  ## fixme need insuf: tr|fasta|fa
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
  my ($oid, $hdr, $fa)=@_; 
  my $id= $oid; # cleanid($oid);
  my $falen= length($fa);
  my $tlog="";  
  #UPD.drop# my $oidfix= $oid;
  my $rinfo= $geneinfo{$oid};
  $hdr=~s/>\S+\s*/ /; chomp($hdr);
  $hdr=~s/\b(type|strand|offs|aalen|clen|cdslen|len)=[^\s;]+[;]?//g;
  
  #UPD1806 for two+ outputs (cdsorf & utrorf)
  my $nput= 0;
  my $utroid=0; my $utrinfo=[]; 
  if($UTRORF and ($utrinfo= $geneinfo{$oid.'utrorf'})) { $utroid= $oid.'utrorf'; }
  
  ## FIXME: utrorf need to do this also if have oid, but also oid.utrorf .. make special utrorf.aa input subset
  unless($rinfo or $UTRORF) {
    if( $rinfo= $geneinfo{$oid.'utrorf'} ) { $id= $oid.'utrorf'; }  # $oidfix= 
  }
  
  #UPD1806, loop here for 2+ seq outputs
  my $cdsrinfo= $rinfo; my $fasave= $fa; my $falensave= $falen;
  #------- outloop -------------
  for $rinfo ($cdsrinfo, $utrinfo) { 
  $rinfo ||= [];
  next unless(@$rinfo > 6);
  $fa= $fasave; $falen=$falensave; # UPD1806
  
  my($td,$gd,$mapqual,$aaqual,$trlen,$cdsor,$cdsoff)= @$rinfo; # ,$lotag,$nam1,$dbxref
  $trlen= $falen; # use what we have as input
  $id= $td; # UPD1806
  
  ## SKIP cds output unless cdsoff
  unless($cdsoff) { loggit(1,"MISS: $oid\t$falen,$aaqual") unless($MISSNOLOG); return 0; } # option nolog here
  
  ## MISS err neg cdsoff:  offs=-6355--800
  ## .. bug in cdna_bestorf.pl? or here? no DATA bug, got dup ids, diff seqs, in .aa,.cds and .tr *!*!*!*!*!
  ## .. all from one trasm set mixup has same id prefix: daphmag/rnas/soap5xun3set soap5xinb3_co, soap5xinb3_cop2, 
  ## MISS: dmag5xun3cosoapk21loc758t5 type=cdna; aalen=1851,79%,complete; clen=592;  strand=+; offs=-6355--800;   revcds:-6355--800,
  ## MISS: dmag5xun3cosoapk21loc2222t1 type=cdna; aalen=304,44%,complete-utrbad; clen=455;  strand=+; offs=-1530--616;    revcds:-1530--616,
  #bad >dmag5xun3cosoapk21loc2195t2 type=cdna; aalen=569,81%,partial3; clen=618;  strand=+; offs=-1088-618;  revcds:-1088-618,
  #bad >dmag5xun3cosoapk29loc38286t1 type=cdna; aalen=696,60%,complete; clen=3467;  strand=+; offs=-2261--171;  

  my($cdsb,$cdse)= split/[-,]/,$cdsoff; # this SHOULD be offset in non-rev tr; any problems?
  if($cdsb > $cdse and $cdsor ne "-") {
    $cdsor="-";
  }
  
  if($UTRORF and $td =~ /utrorf/) { # was $oidfix =~ /utrorf/
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
  my $def= "$id type=$typ; aalen=$aaqual; $atlen strand=$cdsor; offs=$cdsoff;$hdr";
  if($cdse < 1 or $len<1) { loggit(1,"MISS: $def  $tlog"); return 0; }
  
  $cdsfa =~ s/(.{60})/$1\n/g; 
  print $outh ">$def\n$cdsfa\n"; 
  $nput++; $idput{$id}++; #was $oidfix
  loggit(0,">$def  $tlog"); 
  }  #------- end outloop -------------

  return $nput; # return($len,$ulen); ?
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

  #... old parts from  evigene/scripts/rnaseq/asmrna2ncbitsa.pl

=cut


