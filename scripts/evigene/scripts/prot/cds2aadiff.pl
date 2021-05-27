#!/usr/bin/env perl
# cds2aadiff.pl : calc genes.cds2aa <> genes.aa 
# cat $genes.cds2aa $genes.aa  | cds2aadiff.pl  > $genes.diffa
# fixme: cds2aadiff.pl -cds $genes.cds2aa -aa $genes.aa > $genes.diffa

=item how to make cds2aa

  faTrans pub3i.cds stdout |\
  perl -ne 'if(/^>(\S+)/){$d=$1; puta() if($g); $g=$_; $aa=""; } \
  else { chomp; $aa.=$_; } END{puta();} sub puta{ $aa=~s/Z$/\*/; $aa=~s/(.{60})/$1\n/g; \
  $aa.="\n" unless($aa=~m/\n$/); print $g,$aa; }'\
   > pub3i.cds2aa     

=cut

use strict;
use Getopt::Long;

my $SHOWOK=$ENV{ok}||0;
my $CUTP= $ENV{cut} || 'aalen|Name|oid';
my @ERRS=qw(OK DIFF SMAL BIGG MISS MISC);
my $STOPCHECK= $ENV{stop}||0;
my $GAPOFFBYONE= 1;
my $GAPTEST= 0;

my($incdsaa,$incds,$inaa);

my $optok= GetOptions(
  "aa=s", \$inaa,
  "cds=s", \$incds,
  "c2aa=s", \$incdsaa, # add output default: inaa.diffa
  "stopcheck!", \$STOPCHECK, 
  "gapoffbyone!", \$GAPOFFBYONE, # gap3codon ?
  "testgaps!", \$GAPTEST, 
  "ok!", \$SHOWOK,
  "cut=s", \$CUTP,
  );
die "usage: cat genes.cds2aa genes.aa  | cds2aadiff.pl  > genes.diffa ...
  opts: 
  -aa=genes.aa input , -c2aa=genes.cds2aa input
  -cds=genes.cds, not cds2aa, uses JKent faTrans
  -cut='aalen|Name|oid' : header flag in genes.aa to signal input switch
  -[no]gapoffbyone [$GAPOFFBYONE] : policy on unambigous codon with 3-position gap char: is it gap or valid aa?
  -[no]stopcodon [$STOPCHECK] : include stopcodon at end of aa in check
  -[no]ok [$SHOWOK] : show OK set;  
  \n"  unless($optok);
  
our($ina,$haveca,$nao)=(0) x 10;
our(%ca,%gotca,%nerr);
our($g,$aa)=('','');

=item gap X or * stop diff

#diffaa: nOK=34949, nDIFF=82, nSMAL=7, nBIGG=1, nMISS=7, nMISC=0
nMISS=7 : all but one are tiny aa/cds left out by MINCDS. these may be drop genes
nSMAL=7 : 3 are human: prot2gene, as with catfish need kf2protgene.cds fix, 4 are kf4bAUG : what? 
nBIGG=1 : bad catfish.cds from kf2protgene.gff,mrna,cds,aa; redo kf2protgene.cds from mrna
nDIFF=82: most are either X/* or X/validaa single base/aa diffs

-- cases of X or *, are these soft policy diffs or seq diffs?
-- most remaining nDIFF=82 are either X/* or X/validaa single base/aa diffs

  n=13 with X/* are all AUG models, 6 from kf1, 7 from kf2, should have been aa-fixed for this problem.
  3 of 7 kf2 are TE genes, others Unchar mabye also TE.  Funhe2EKm017681t1 has useful Ortho FISH11E_G17285; drop others as needed.
  
DIFF.0  Funhe2EKm031595t1
Da: MAWCMYKGKCKFHCTNGNCLRLGSLICNQLNNCGDNSDEENCPVVTQHPPPGIFNXAGRRMKDAEVNLIHFSILEL
Dc: .......................................................*....................
Funhe2EKm031595t1 oid=Funhe5EG033461t1,AUGepir1s3208g61t1 aalen=76,7%  << X may be Aug special gap-fix, find other model? drop?
MAWCMYKGKCKFHCTNGNCLRLGSLICNQLNNCGDNSDEENCPVVTQHPPPGIFNX<<AGRR

=item off-by-one gap X 

** diff aa processors have diff policy on when to start gap if 2 of 3 codon bases exist before NNN.
-- which is right?  cdna_bestorf is using 2codon unambig = amino, others faTrans, ncbi?? require 3 base codon.

#diffaa: nOK=34949, nDIFF=82, nSMAL=7, nBIGG=1, nMISS=7, nMISC=0
-- after GAPOFFBYONE runs fix2, ie all? offby1 gaps gone

#diffaa: nOK=34630, nDIFF=401, nSMAL=7, nBIGG=1, nMISS=7, nMISC=0
-- after GAPOFFBYONE terminal X fix

DIFF.0  Funhe2EKm033821t1 : remaining gap-run bug?
Da: LESPRQTFTPALALLAEQRRRARSTLSWRFESAGTKVPRRLKCRQKLTRIQPRRIPALRSTAFSSYLHTHTXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXHTHTHTHIHT
YLPSVLLFPLQAFTYHSPESRSTHKLSSWFSVSACREKQRPKCTCLISRSGQAVTHSPASNVAPSSPAVFFPPLFRHYSYSAXXXXXXXXXXXXXXXXXXXXXEALFPPTFQALQLFSV
LCPNPLVCVCQISIGALCAWAARTGNCSLDRVGPLFLLLNTPCLLHRVILWIYGGGEEKKNCM
Dc: ...................................................................................................................
.................................................................................X.....................................
...............................................................

ggrep -A9 Funhe2EKm033821t1 kfish2rae5d.main4.cds2aa
>Funhe2EKm033821t1
LESPRQTFTPALALLAEQRRRARSTLSWRFESAGTKVPRRLKCRQKLTRIQPRRIPALRS
TAFSSYLHTHTXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXHTHTHTHIHTYLPSV
LLFPLQAFTYHSPESRSTHKLSSWFSVSACREKQRPKCTCLISRSGQAVTHSPASNVAPS
SPAVFFPPLFRHYSYSXXXXXXXXXXXXXXXXXXXXXXEALFPPTFQALQLFSVLCPNPL
VCVCQISIGALCAWAARTGNCSLDRVGPLFLLLNTPCLLHRVILWIYGGGEEKKNCM*

ggrep -A9 Funhe2EKm033821t1 kfish2rae5d.main4.aa
>Funhe2EKm033821t1 type=protein; aalen=297,78%,partial5; clen=1138; offs=1-894; oid=Funhe2Exx11m021246t1,Fungr1EG3m013257t1; organism=Fundulus_heteroclitus; evgclass=main;
LESPRQTFTPALALLAEQRRRARSTLSWRFESAGTKVPRRLKCRQKLTRIQPRRIPALRS
TAFSSYLHTHTXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXHTHTHTHIHTYLPSV
LLFPLQAFTYHSPESRSTHKLSSWFSVSACREKQRPKCTCLISRSGQAVTHSPASNVAPS
SPAVFFPPLFRHYSYSAXXXXXXXXXXXXXXXXXXXXXEALFPPTFQALQLFSVLCPNPL
VCVCQISIGALCAWAARTGNCSLDRVGPLFLLLNTPCLLHRVILWIYGGGEEKKNCM

#diffaa: nOK=34617, nDIFF=414, nSMAL=7, nBIGG=1, nMISS=7, nMISC=0
-- after GAPOFFBYONE iterate gap runs .. still buggy here got some offby1 X cases yet
   also terminal codon is sometimes X from faTrans .. is this offby1 cds length? or gap/end of mrna?
   .. treat end X like inner offby1 gap?
   
#diffaa: nOK=34561, nDIFF=470, nSMAL=7, nBIGG=1, nMISS=7, nMISC=0
-- after GAPOFFBYONE=1  .. remaining nDIFF due to several gap runs .. need smarter gapoffbyone 
DIFF.0  Funhe2EKm027867t1
Da: MASFQLLLFLGYFAAARAGVVLNCTNDYEKISCQLSAKQCSEYAVNIRNDQGYGAYNCSLKRCNSELCCCSVKITFITGETYTVETFNGSEKVDTKTLDLWESFKPKTPTIVSVK
EYEGIFTVNWKTNMGADFVGTLXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXIVSVKEYEGIFTVNWKTNMGADFVGTLTAEVTVSKKGEQGKVYGNIRPATDEGLQSFD
INGHDLEPGTPYVVSVRSYTDRSQMFSDRSQECEFKTXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXPGTAYVVSVRSYTDRSQMFSDRSQECEFKTPPSKSYWILGIIIAISAFGIIL
                                    ^gap2off1
Dc: ...................................................................................................................
.......................................................................................................................
....................................X..................................................................................

///
#diffaa: nOK=33014, nDIFF=2017, nSMAL=7, nBIGG=1, nMISS=7, nMISC=0
-- most of these nDIFF are gap X policy diff.

DIFF.0  Funhe2EKm007301t1
Da: MATYITELDVSLNDAEATYLRTQGFKQIYVNLNMNGSPEVYLWYKEGSSGPITRIQFSFTNSMAEGLNASGFEKVKKNLNTGTCGDDVLLWYFKGTTESHVPITKIEVSRDDEAR
NKVYLWYKKEDPGHKPIQAXXXXXXXXXXXXXXXXXXXXXXXXXXXXXHKPIQAITLLSNQKLIKDYEKAEVTVLQENLNSGCGCPCFPATPLYLGYYN
Dc: ...................................................................................................................
..................X................................................................................

DIFF.0  Funhe2EKm008138t1
Da: KSIFRLFCCCCRHVLALNEPTVSRGLFAQTTTRGFISVRKPLIFSRLNRGRDAAXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
Dc: .....................................................X.............................................................

DIFF.0  Funhe2EKm027610t1
Da: MRLREEERVGKRVCFVQTPPSWRRALSAPSLQLPPAFRVTEPCGSKLRGRVPCSRSKPARPISISPXXXXXXXXXXXXXXXXXXXXXXXXSSRSKPSRPISSSPSGYPSADMDFS
Dc: .................................................................X.................................................

=item stopcodon mixup

genes.aa lacks '*' end stop by design, faTrans adds. this now screws difa.  
parse genes.aa hdr for 'complete' and add * or skipit
BIGG.0,1        Funhe2EKm005219t1
Da: MNFLKK MDKAGE
Dc: ...... ......*

=cut

# FIXME: $inaa, $incdsaa instead of STDIN
if($incds) { 
  ($incdsaa,$nao)= makecds2aa($incds,$incdsaa); # set $ca{$g}=$aa also;  
  $haveca=1 if($nao);
} elsif($incdsaa) {
  ($g,$aa)=('','');
  open(INA,$incdsaa) or die "ERR: $incdsaa";  
  while(<INA>) {
    if(/^>(\S+)/){ $ca{$g}=$aa if($g); $g=$1; $aa=""; } elsif(/^\S/) { chomp; $aa.=$_; } 
  }
  $ca{$g}=$aa if($g); $haveca=1;
}

($g,$aa)=('','');
$ina=($haveca)?1:0;
my $fin= *STDIN;
if($inaa) { open(INA,$inaa) or die "ERR: $inaa"; $fin= *INA; $ina=1; }
while(<$fin>) {
  if(/^>(\S+)/){ my $d=$1; puta($g,$aa) if($g); $g=$d; $ina=1 if(/$CUTP/); $aa=""; } 
  elsif(/^\S/) { chomp; $aa.=$_; } 
}
puta($g,$aa) if($g); 
cmiss();
print "#diffaa: ", (join", ",map{ "n$_=".($nerr{$_}||0) } @ERRS),"\n"; 


sub makecds2aa # drop this for evigene/bestorf ??
{
  my($cdsin,$cdsaaout)=@_;
  unless($cdsaaout) { ($cdsaaout=$cdsin) =~ s/\.\w+$//; $cdsaaout.=".cds2aa"; }
  open(FAT,"faTrans $cdsin stdout |") or die "ERR: faTrans $cdsin stdout";
  open(FO,'>',$cdsaaout) or die "ERR: write $cdsaaout";
  my($ha,$ga,$ca)=("","",""); # %ca is global store for cds2aa
  our($no)=(0);
  sub putca{ my($ha,$ga,$ca)=@_; our($no); return unless($ha);
    $ca=~s/Z$/\*/; $ca{$ga}=$ca; # URK store ca before prettify
    $ca=~s/(.{60})/$1\n/g; $ca.="\n" unless($ca=~m/\n$/); 
    print FO $ha,$ca; $no++; }
  while(<FAT>) { 
    if(/^>(\S+)/){ my $d=$1; putca($ha,$ga,$ca) if($ha); $ga=$d; $ha=$_; $ca=""; } 
    else { chomp; $ca.=$_; } 
  }
  putca($ha,$ga,$ca); close(FO);
  return ($cdsaaout,$no);  
}

sub cmiss {
  foreach my $g (sort keys %ca) { 
    next if($gotca{$g});
    my $ca=$ca{$g};  my $lca=length($ca);
    my $et="MISC"; my $er="$et"; $nerr{$et}++; 
    print "$er\t$g\t$lca\n" if($lca>9);   
  }
}

sub puta { 
  my($g,$aa)= @_;
  my($i,$et,$er,$egap);
  return unless($g);
  unless($ina){ $ca{$g}=$aa; } # cds2aa here
  else { 
    # upcase aa,ca ? expected here.
    my $ca= uc($ca{$g}); $aa=uc($aa); # MISS handled
    unless($STOPCHECK){ map{ s/\*$//; } ($ca,$aa); }
    my $isdiff=($ca ne $aa)?1:0;

    my $cgap= index($ca,'X');
    my $agap= index($aa,'X');
    $egap="";
    
    ## option here to also check removing XXX runs form aa, pick best of 2 ?
    if($GAPTEST and $isdiff and ($cgap>=0 or $agap>=0)) {
      my($aax,$cax)= ($aa,$ca);
      map{ s/XX+/X/g } ($aax,$cax);
      my($wca,$waa)= map{ length($_) } ($ca,$aa);
      my($wcax,$waax)= map{ length($_) } ($cax,$aax);
      my($dap,$neq,$nne)= difa(0,$aa,$ca); 
      my($dapx,$neqx,$nnex)= difa(0,$aax,$cax); 
      my $peq = ($waa >0) ? 100*$neq/$waa : 0;
      my $peqx= ($waax>0) ? 100*$neqx/$waax : 0;
      if($peqx > 1.10*$peq) {
        # flag it and use aax,cax instead ??
        map{ $_= int($_) } ($peq,$peqx);
        $egap="Degap: $peqx%,$neqx/$waax,$wcax vs $peq%,$neq/$waa,$wca";
        $aa=$aax; $ca=$cax;
        $cgap= index($ca,'X');
        $agap= index($aa,'X');
        $isdiff=($ca ne $aa)?1:0;
        }
      }
        
    if($GAPOFFBYONE and $isdiff and ($cgap>=0 or $agap>=0)) { # need to iterate thru gapruns here ..
      my($wca,$waa)= map{ length($_) } ($ca,$aa);
      my $wmin=($wca>$waa)?$waa:$wca;
      my $ioff=0; my $did1=0; 
      while( $ioff < $wmin ) {
        # another case is single codon X vs OK .. how to handle? is this single N in CDS?
        my $xc=index($ca,'X',$ioff);  
        my $xa=index($aa,'X',$ioff); 
        if($xa>=0 and abs($xc-$xa)==1) { # fix here check for next run if xc==xa
          if($xa>$xc) { substr($aa,$xc,1)='X'; } 
          elsif($xa<$xc) { substr($ca,$xa,1)='X'; } 
          $isdiff=($ca ne $aa)?1:0; $did1++;  # should OK here add note?
          }
        # scan for next gaprun even if xc==xa last run
        if($isdiff and $xa>=0) { $ioff= $xa; $ioff++ while($ioff < $wmin and substr($aa,$ioff,1) eq 'X'); }
        unless($isdiff) { $ioff=$wmin; }
        elsif($xa < 0 or $xc < 0 or abs($xc-$xa) > 1) { $ioff=$wmin; } # cant fix, stop
      }

      if($isdiff and $wca==$waa and $ca =~ m/X$/ and $aa !~ m/X$/) { # terminal X/ok check, both ways?
        if(substr($ca,0,$wca-1) eq substr($aa,0,$waa-1)) { 
        substr($aa,$waa-1,1)='X'; $isdiff=($ca ne $aa)?1:0; 
        }
      }

    }
    
    
    if($isdiff) {
      $ca=~s/Z/\*/g; 
      my($wca,$waa)= map{ length($_) } ($ca,$aa); my($dap,$neq,$nne)=("",0,0);
      my $da= $wca - $waa; 
      my($xca,$xaa)=($ca,$aa); if($STOPCHECK) { map{ s/\*$//; } ($xca,$xaa); }
      unless($ca){  $et="MISS"; $er="$et"; $nerr{$et}++; }
      elsif( $da>0 and ($i=index($xca,$xaa)) >= 0) { ($dap,$neq,$nne)=difa(0,$aa,$ca); 
        $et="BIGG";  $er="$et.$i,$da"; $nerr{$et}++;  } 
      elsif( ($i=index($xaa,$xca))>=0) { ($dap,$neq,$nne)= difa($i,$aa,$ca); 
        $et="SMAL";  $er="$et.$i,$da"; $nerr{$et}++; } 
      else { ($dap,$neq,$nne)=difa(0,$aa,$ca);
        $et="DIFF"; $er="$et.$da"; $nerr{$et}++; 
      } 
      my $peq= ($waa>0) ? int(0.5 + 100*$neq/$waa) : 0; 
      my $eqnote= "$peq%,$neq/$waa,$wca;$egap";
      print "$er\t$g\t$eqnote\n$dap";  $gotca{$g}= $et;
    } else { my($wca,$waa)= map{ length($_) } ($ca,$aa);
      $et="OK"; print "$et\t$g\t100%,$waa/$waa,$wca;$egap\n" if($SHOWOK); $gotca{$g}= $et; $nerr{$et}++; 
    } 
  } 
} 

## add count of difa residue diff/same, eg. '.'
sub difa { 
  my($io,$aa,$ca)=@_; my @ca=split"",$ca; my @aa=split"",$aa;
  my $co=($io>0)?$io:0; unshift(@ca, ('-') x $co) if($co); my @ba=@ca; 
  my($nsam,$ndif)=(0,0);
  for my $i (0..$#aa) { my($e,$c)=($aa[$i],$ca[$i]); if($e eq $c){ $ba[$i]='.'; $nsam++; } else { $ndif++; } } 
  my $ba=join"", @ba; 
  return ("Da: $aa\nDc: $ba\n",$nsam,$ndif); 
} 


