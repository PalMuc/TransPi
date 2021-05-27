#!/usr/bin/env perl
# cdsqual.pl : mod from evigene/scripts/cdna_evigenesub.pm and prot/aaqual.sh
# adds coding potential calcs on CDS seq: Fickett + codon-usage hexamers, coding-class (code,nocode,unkn)
# tested vs cpat calcs, same fickett score, similar hex (diff calc), similar class call
# 2016.05.04, d.gilbert

use strict;

use constant VERSION => 'alpha:2016.05.04';

## need global params for utrbad/poor; see below .. opt to reduce/drop use of aabad/poor call of Noncode
our $pCDSbad = $ENV{pcdsbad}  ||30; # adjust DOWN *?
our $pCDSpoor= $ENV{pcdspoor} ||60;
our $BAD_GAPS= $ENV{aagapmax} || $ENV{BAD_GAPS} || 15;  # % gaps in AA
our ($AAQUALF,$AAQUALH) = ("",undef,undef);
our $NoStopCodon= 0;
our $USESelenocysteine= 0;

our $MINUTR=$ENV{minutr}||300; # ~ 300b "fixed" average utr sizes, maybe too low, 
# see sample of org gene set pcds for small aa .. per species min-utr, ranging from 400b insect .. 700b fish .. 1000b mouse
# should adjust %cds bad for cds size, for large cds> 10k, pcds >=90%, for small cds < 300, pcds may be 33%
# should replace utrbad/poor flags with coding potential flags, adding some cds-seq cp calcs.

use constant { kAAQUAL_MAX => 3, kAAQUAL_MIN => -3, kAAQUAL_NONE => 0, };
use constant { FNoncode => 0.74, FCoding => 0.95 };  
use constant { HNoncode => -1e-4, HCoding => 1e-4 };  #?? hexamer score, dunno what range to call
use constant { CPATNoncode => 0.20, CPATCode => 0.40 };  # paper picks 0.36 from human nc/cds, allow slop here

# GetOptions( in,out,other..);
my $debug=$ENV{debug}||0;
my $doCompareCpat=$ENV{cpat}||0;
my $hexcode={};
my ($hexfile,$cdsin)=("","");
our %hexcodencode; # == inline table, default?
 
MAIN: {
  ($hexfile)= shift @ARGV;  ## figure out how to embed this, 4096 rows, 2+ cols, 
    # as  hash{codonhex}=loglike_coding_nocode    .. see below 800 line %hexcodencode data.
  ($cdsin)= shift @ARGV; # foreach cdsin ..
  
  if($doCompareCpat) {
    my $cpattab = $hexfile;
    my $cqualtab= $cdsin;
    my($outfile)= shift @ARGV; $outfile||="";
    die "usage: env cpat=1 $0 cpat.out cds.qual [stdout|outfile] \n" unless(-f $cpattab and -f $cqualtab);
    cpat_compare($cpattab,$cqualtab,$outfile);
    # env cpat=1 cdsqual.pl -comparecpat evg24m2banofun.cpat.out evg24m2banofun5ht1.cds.qual  > evg24m2banofun5ht1.both.tab
    
  } else {
    my $usehexinline=0;
    if($hexfile and -f $hexfile and not -f $cdsin) {
      $cdsin=$hexfile; $hexcode=\%hexcodencode; $usehexinline=1;
    } else {
      die "usage: $0 [hexamers] input.cds \n" unless($cdsin and -f $cdsin); # and $hexfile and -f $hexfile 
    }
    my $qualout= ($debug)?"stdout":"$cdsin.qual"; # need in/out other GetOptions
    $qualout =~ s/\.gz//;
    $hexcode= read_hexcode($hexfile) unless($usehexinline);
    makeCdsQual($cdsin,$qualout);
  }
}


sub cpat_compare {
  my($cpattab,$cqualtab,$outfile)=@_;
  # use constant { CPATNoncode => 0.20, CPATCode => 0.40 };  # paper picks 0.36 from human nc/cds, allow slop here
  my(%cpat);
  open(F,$cpattab) or die "Missing cpat.out table:$cpattab";
  while(<F>) {
    my($id,$tw,$cw,$pfk,$phx,$cp)=split; next unless($tw>0);
    my $cl=($cp >= CPATCode)?"Code":($cp <= CPATNoncode)?"Nonc":"Unkn";  
    map{ s/(\.\d\d\d\d\d)\d+/$1/; } ($cp,$phx); # ($cp,$phx)=sprintf "%.5g",$cp,$phx;
    $cpat{uc($id)}="CY:$cl/$cp,$pfk,$phx"; # cpat upcases all ids :(
  } close(F);
  
  open(F,$cqualtab) or die "Missing cds.qual table:$cqualtab";
  my $outh= *STDOUT;
  if($outfile) { open(O,'>',$outfile) or die "writing table $outfile"; $outh=*O; }
  my($ncall,$sagree)=(0,0);
  while(<F>) {
    # my($id,$cw,$gp,$aw,$tw,$cf,$ch,$off,$oid)=split; # changed col format: cf,ch => "allc,cf,ch" one col
    my($id,$cw,$gp,$aw,$tw,$call,$off,$oid)=split;  
    my($ctot,$cf,$ch)=split",",$call;
    my $cpat=$cpat{uc($id)}||0; 
    map{ s/Noncode/Nonc/g; s/Unknown/Unkn/g; } ($call,$ctot,$cf,$ch); 
    ## count agree/disagree and report.. add cdsqual overall Code/Nocode/Unkn call/score
    ## ctot == combined cdsqual class
    if($cpat) { my $agree= ($cpat=~/$ctot/)?1:0; $sagree+= $agree; $ncall++; }
    print $outh join("\t",$id,$cw,$tw,$call,$cpat,$aw)."\n";  #  $cf,$ch  
  } 
  my $da=$ncall-$sagree;
  my($pa,$pd)=map{ int(100*$_/$ncall) } ($sagree,$da);
  print $outh "#sum: agree:$sagree,$pa%  disagree:$da,$pd%  total:$ncall\n";
  close(F); close($outh);
}

=item cpat_compare_noho class

  - called Code by cdsq, Noncode by cpat, lack homology 
  
perl -ne 'if(/^Ano/ and /\t(Code|Nonc|Unk)/) { @cp=split;
($id,$cw,$tw,$cclass,$yclass,$aaq)=@cp; $cc{$id}=$cclass; 
$cy{$id}=$yclass;  } elsif(/^>(\w+)/) { $id=$1;
($nt)=m/notes=([^;\s]+)/; $ho=(m/aaref:/)?"hoho":"noho";
if($cc=$cc{$id}) { ($cl)=split",",$cc; $cy=$cy{$id}; $cy=~s,/.*,,;
($aaq)=m/aalen=([^;\s]+)/; ($oid)=m/oid=([^;\s]+)/; print
join("\t",$id,$ho,$cc,$cy,$aaq,$nt)."\n" if($ho eq "noho" and
$cl=~/Code/ and $cy=~/Nonc/);  } } ' \
 evg24m2banofun5ht1.both.tab7  evg24m2banofun5ht1.cds \
  > evg24m2banofun5ht1.noho-disagree.tab7

  disagree n=3515, 1434 are both CF+CH:Nonc, but 227 of these are sigKaKs 
  - which can be shifted to Noncode reliably? or use Unknown as unreliable class, need ho test
    a. shortest (<60aa) without CF|CH Code call (and lack homolog): either Noncode or Unknown
    b. ignore CY:Nonc for longish partial, due to mis-orf call
    c. any CF+CH:Nonc also partial = Nonc ? n=308
    d. any CF+CH:Nonc other??
    
   also of 9235 cdqual:Nonc evg24m2banofun5ht1.both.tab7, 202 have sigKaKs, 
     of nnn cdqual:Code, 648 have sigKaKs (these are all noho?)
    
.. should also check sig KaKs scores, 925 ids ( culicifacies+mimimus )
cat  $anon/evg2anofunzh/annot/calc6a/anofunevg24m_*.noho*.kaks | cut -f1,4-8 | egrep -v 'Value|NA|nan' |\
  sort -k6,6g | head -1350 | cut -f1 | sort -u | grep 't1$' > anofunevg24m.sigkaks.ids
  
  n=135 of these are in evg24m2banofun5ht1.noho-disagree.tab7
grep -F -f anofunevg24m.sigkaks.ids evg24m2banofun5ht1.noho-disagree.tab7 |  head
Anofunz4kEVm013616t1	noho	Code,CF:Nonc/0.5545,CH:Nonc/-0.006156,CA:Code	CY:Nonc	156,76%,partial5	0,0,chrmap:100a,99i,471l,2x,KB668706:122506-124509:+,pflag:0,phi:100/100
Anofunz4kEVm013776t1	noho	Code,CF:Nonc/0.4147,CH:Nonc/-0.04469,CA:Code	CY:Nonc	154,80%,complete	0,0,chrmap:100a,100i,465l,1x,KB668759:733998-734462:+,pflag:0,phi:100/100
Anofunz4kEVm015750t1	noho	Code,CF:Nonc/0.4392,CH:Nonc/-0.008275,CA:Code	CY:Nonc	136,73%,complete	0,0,chrmap:100a,100i,411l,1x,KB668289:371403-371817:-,pflag:0,phi:100/100
Anofunz4kEVm022957t1	noho	Code,CF:Code/0.9511,CH:Code/0.00316,CA:Nonc	CY:Nonc	87,48%,complete-utrpoor	0,0,chrmap:100a,100i,264l,2x,KB668221:1283114-1283464:+,pflag:0,phi:100/100
Anofunz4kEVm021364t1	noho	Code,CF:Unkn/0.8344,CH:Nonc/-0.006755,CA:Code	CY:Nonc	98,65%,complete	0,0,chrmap:100a,100i,297l,1x,KB668803:841062-841358:+,pflag:0,phi:100/100
Anofunz4kEVm018043t1	noho	Code,CF:Code/1.0346,CH:Code/0.009592,CA:Code	CY:Nonc	122,82%,partial5	0,0,chrmap:100a,99i,369l,3x,KB668681:861746-864769:+,pflag:0,phi:100/100
Anofunz4kEVm017473t1	noho	Code,CF:Nonc/0.4203,CH:Nonc/-0.001472,CA:Code	CY:Nonc	125,84%,complete	0,0,chrmap:100a,100i,378l,1x,KB668737:278442-278819:+,pflag:0,phi:100/100
Anofunz4kEVm019118t1	noho	Code,CF:Nonc/0.6008,CH:Nonc/-0.0154,CA:Code	CY:Nonc	114,77%,complete	0,0,chrmap:96a,100i,345l,1x,KB668892:23706-24037:-,pflag:0,phi:100/100
Anofunz4kEVm022185t1	noho	Code,CF:Nonc/0.489,CH:Nonc/-0.0006118,CA:Code	CY:Nonc	92,74%,complete	0,0,chrmap:100a,100i,279l,1x,KB668589:106524-106804:-,pflag:0,phi:99/82/-sense/altpar82
Anofunz4kEVm023280t1	noho	Code,CF:Unkn/0.7424,CH:Nonc/-0.001904,CA:Code	CY:Nonc	83,72%,complete	0,0,chrmap:100a,100i,252l,1x,KB668670:494485-494736:+,pflag:0,phi:100/100

  
eg. problem case, 171,72%,complete is goodish, but -CF+CH sasy to skip
Anofunz4kEVm012599t1	noho	Code,CF:Nonc/0.4303,CH:Nonc/-0.006036,CA:Code	CY:Nonc	171,72%,complete	0,0,chrmap:100a,100i,516l,1x,KB668836:1918516-1919031:+,pflag:0,phi:100/100
    
sort evg24m2banofun5ht1.noho-disagree.tab7 | head
Anofunz4kEVm006633t1	noho	Code,CF:Unkn/0.8303,CH:Code/0.01794,CA:Code	CY:Nonc	399,76%,partial5	0,0,pflag:0,phi:100/100/.
Anofunz4kEVm007082t1	noho	Code,CF:Code/1.1428,CH:Code/0.03735,CA:Code	CY:Nonc	375,93%,partial5	0,0,pflag:0,phi:100/80/.
Anofunz4kEVm007603t1	noho	Code,CF:Unkn/0.9122,CH:Code/0.04503,CA:Code	CY:Nonc	346,8%,complete-utrbad	0,0,chrmap:17a,99i,1041l,1x,KB669074:1-179:-,pflag:0,phi:100/100
Anofunz4kEVm008364t1	noho	Code,CF:Nonc/0.6214,CH:Code/0.03786,CA:Code	CY:Nonc	312,12%,complete-utrbad	0,0,chrmap:7a,100i,939l,1x,KB669070:169062-169125:-,pflag:0,phi:100/100
Anofunz4kEVm008383t1	noho	Code,CF:Nonc/0.6272,CH:Code/0.01514,CA:Code	CY:Nonc	311,18%,complete-utrbad	0,0,pflag:0,phi:98/73/.
Anofunz4kEVm009110t1	noho	Code,CF:Code/1.1168,CH:Code/0.03845,CA:Code	CY:Nonc	276,67%,partial5	0,0,chrmap:100a,99i,801l,3x,KB668737:1010827-1014489:+,pflag:0,feq:Anofunz4hEVm007857t1/altpar20.0.0,phi:99/68/.
Anofunz4kEVm009160t1	noho	Code,CF:Code/1.0352,CH:Code/0.01871,CA:Nonc	CY:Nonc	273,17%,complete-utrbad	0,0,chrmap:100a,99i,822l,3x,KB669292:110584-115812:+,pflag:0,phi:98/81/./altparx81
Anofunz4kEVm009218t1	noho	Code,CF:Nonc/0.5088,CH:Nonc/-0.01577,CA:Code	CY:Nonc	270,82%,partial5	0,0,chrmap:100a,100i,813l,1x,KB668389:12072-12884:+,pflag:0,phi:100/100
Anofunz4kEVm009413t1	noho	Code,CF:Nonc/0.4499,CH:Nonc/-0.01076,CA:Code	CY:Nonc	261,74%,partial5	0,0,chrmap:87a,100i,786l,1x,KB668670:879810-880493:+,pflag:0,phi:100/100
Anofunz4kEVm009533t1	noho	Code,CF:Code/1.1211,CH:Code/0.01029,CA:Nonc	CY:Nonc	257,52%,partial5-utrpoor	0,0,chrmap:100a,98i,774l,4x,KB668680:367374-430062:+,pflag:0,phi:100/100
..
sort evg24m2banofun5ht1.noho-disagree.tab7 | tail
Anofunz4kEVm027460t1	noho	Code,CF:Code/0.986,CH:Nonc/-0.00202,CA:Code	CY:Nonc	37,63%,complete	0,0,chrmap:100a,100i,114l,1x,KB668289:323907-324020:-,pflag:0,phi:100/100
Anofunz4kEVm027461t1	noho	Code,CF:Nonc/0.5514,CH:Code/0.0007878,CA:Code	CY:Nonc	37,63%,complete	0,0,chrmap:100a,98i,114l,1x,KB668478:531907-532020:+,pflag:0,phi:100/100
Anofunz4kEVm027462t1	noho	Code,CF:Unkn/0.8692,CH:Nonc/-0.003086,CA:Code	CY:Nonc	37,62%,complete	0,0,chrmap:97a,99i,114l,1x,KB669480:445080-445190:-,pflag:0,phi:100/100
Anofunz4kEVm027463t1	noho	Code,CF:Unkn/0.8689,CH:Nonc/-0.001934,CA:Code	CY:Nonc	37,61%,complete	0,0,chrmap:89a,100i,114l,1x,KB668661:404690-404790:+,pflag:0,phi:100/100
Anofunz4kEVm027464t1	noho	Code,CF:Unkn/0.8984,CH:Nonc/-0.004554,CA:Code	CY:Nonc	36,61%,complete	0,0,chrmap:100a,99i,111l,1x,KB668780:343498-343608:-,pflag:0,phi:100/100
Anofunz4kEVm027837t1	noho	Code,CF:Nonc/0.4935,CH:Nonc/-0.005373,CA:Code	CY:Nonc	102,71%,complete	0,0,chrmap:100a,100i,309l,3x,KB668868:69541-69998:-,pflag:0,feq:Anofunz4hEVm020512t5/altmapxe100.100,phi:100/97/.
Anofunz4kEVm027839t1	noho	Code,CF:Unkn/0.8481,CH:Code/0.0001879,CA:Code	CY:Nonc	90,83%,complete	0,0,chrmap:100a,100i,273l,1x,KB668816:21891-22163:+,pflag:0,phi:100/100/-sense
Anofunz4kEVm027842t1	noho	Code,CF:Nonc/0.6917,CH:Nonc/-0.01061,CA:Code	CY:Nonc	65,71%,complete	0,0,chrmap:100a,100i,198l,1x,KB669181:728380-728577:+,pflag:0,feq:Anofunz4hEVm021875t1/altparx30.0,phi:100/40/./altparx40
Anofunz4kEVm027844t1	noho	Code,CF:Nonc/0.452,CH:Nonc/-0.00286,CA:Code	CY:Nonc	56,73%,complete	0,0,pflag:0,phi:100/93/.
Anofunz4kEVm027845t1	noho	Code,CF:Nonc/0.3983,CH:Nonc/-0.01611,CA:Code	CY:Nonc	68,71%,complete	0,0,chrmap:100a,100i,207l,1x,KB669292:712600-712806:+,pflag:0,feq:Anofunz4hEVm023477t2/altmapxe69.100,phi:100/76/./altmap69xeq

=item cpat_compare  

  v7 compare .. force aaclass == fclass=hclass = noncode
  #sum: agree:21520,77%  disagree:6256,22%  total:27776
    17762 Code, 9235 Noncode cdq, moved 3200 to Noncode
      3893 of cdq:Code are CY:Noncode, 
      1434 of these also CF:Nonc+CH:Nonc,
            all have >=70% pCDS, >100aa, but some partial, etc.
  
  ** check calls against homology score.
perl -ne 'if(/^Ano/ and /\t(Code|Nonc|Unk)/) { @cp=split;
($id,$cw,$tw,$cclass,$yclass,$aaq)=@cp; $cc{$id}=$cclass; 
$cy{$id}=$yclass;  } elsif(/^>(\w+)/) { $id=$1;
($nt)=m/notes=([^;\s]+)/; $ho=(m/aaref:/)?"hoho":"noho";
if($cc=$cc{$id}) { ($cl)=split",",$cc; $cy=$cy{$id}; $cy=~s,/.*,,; 
print join("\t",$id,$ho,$cl,$cy)."\n";  } } ' \
  evg24m2banofun5ht1.both.tab7  evg24m2banofun5ht1.cds | cut -f2,3 | sort | uniq -c 
cdqual class:
12608 hoho	Code
  287 hoho	Nonc  * bad, but not too many
  351 hoho	Unkn
 5154 noho	Code  * too many? check cdqual vals
 8948 noho	Nonc
  500 noho	Unkn
cpat class: cut -f2,4
12146 hoho	CY:Code
  594 hoho	CY:Nonc * more bad
  498 hoho	CY:Unkn
 1245 noho	CY:Code + not too many
12223 noho	CY:Nonc
 1070 noho	CY:Unkn
    
  v6 compare .. utrpoor change, aaclass == fclass when aascore<max and fclass == hclass == nocode
  #sum: agree:20166,72%  disagree:7610,27%  total:27776
   .. of cdqual, 21113 Code, 6033 Noncode 
    (7143 CY:Nonc, 3908 also CF+CH:Nonc, 2043 are 100-200aa/complete/utrok), << biased too much to Code?
  
  ## v5 compare
  #sum: agree:18449,66%  disagree:9399,33%  total:27848
   7000 are cqual Code, cpat Nonc, 2400 of these utrpoor (40%..), 
    .. 2000 have CF:Code or CH:Code
    .. 4000 have CF:Nonc and CH:Nonc .. 1443 utrpoor is problem set? 
    .. others include partial/complete, utrok, may be bad orf calls by cpat?

=cut
  
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

use constant { kAABIG => 300, kAATINY => 100 };

sub aaqualscore
{
  my($mqual,$aaw)= @_;  $mqual ||="missing";

  # $mvq++ if($aaw >= 300); # biggies  
  # $mvq-- if($aaw < 90); # shorties

  ## FIXME: option to ignore or reduce utrpoor/utrbad class as those appear to be unreliable
  use constant { qUTRBAD => 1.5, qUTRPOOR => 0.6 };
  
  my $mqv= kAAQUAL_NONE; 
  if($mqual =~ /complete/) { $mqv = kAAQUAL_MAX; } 
  elsif($mqual =~ /partial[35]/) { $mqv = kAAQUAL_MAX - 1; }
  if($mqual =~ /utrbad/) { 
    $mqv = ($mqv >= kAAQUAL_MAX-1 and $aaw>= kAABIG) ? kAAQUAL_NONE
      : $mqv - qUTRBAD; # change, short aa + utrbad, ignore complete..
    } 
  elsif($mqual =~ /utrpoor/) { $mqv -= qUTRPOOR; } #?oldpoor# $mqv-- if($aaw<kAATINY);  # $mqv-- if($aaw<kAABIG);
  if($mqual =~ /gapbad/) { $mqv -= 1; } # or -2?
  
  if($aaw >= kAABIG and $mqv < kAAQUAL_MAX-1) { 
    $mqv= ($mqv >= kAAQUAL_NONE-1) ? kAAQUAL_MAX-1 : $mqv+1; 
  }
  elsif($aaw>0 and $aaw<kAATINY) { $mqv = kAAQUAL_MIN if($mqv < kAAQUAL_MAX-1); } # call these short-bad
  else { $mqv = kAAQUAL_MIN+1 if($mqv < kAAQUAL_MAX-1); } # call these short-bad also?
  # ?? mqv-- if aasmall and utrpoor?
  
  return $mqv; # range is now -3..+3, was -2 .. +2
}

=item bad calls

Anofunz4kEVm017181t1	384	127,3%,complete-utrbad	666	Coding,CF:Noncode/0.723,CH:Noncode/-0.004062
 .. aaqual should be Unknown or Noncode .. kAABIG doesnt apply, 3% utrbad is very bad..
 
=cut

# decide Coding/Noncoding/Dunno from mix of aasize, cds/mrna ratio, fickett, hexamer scores
#     my($cpqual)= codingpotQual($al,$fclass,$fval,$hclass,$hval); ## add new 
sub codingpotClass {
  my($aaqual,$fclass,$fval,$hclass,$hval,$isutrorf)= @_;
  my($aaw,$pcds,$mqual)=split",",$aaqual;  $pcds=~s/%//;
  # my $paaw = ($aaw > 120)? 1 : $aaw/120; # heuristic cutoff or ratio?

  ## oid= utrorf is bad qual ..
  
  use constant qMINPCDS => 40;
  # AAQUAL Current range is -3..+3 (kAAQUAL_MIN..kAAQUAL_MAX) 
  my $mqv= aaqualscore($mqual,$aaw); # kAAQUAL_NONE;
  
  # if($mqual =~ /complete/) { $mqv = kAAQUAL_MAX; } 
  # elsif($mqual =~ /partial[35]/) { $mqv = kAAQUAL_MAX - 1; }
  # if($mqual =~ /utrbad/) { $mqv = ($mqv >= kAAQUAL_MAX-1) ? kAAQUAL_NONE - 1 : $mqv - 2; } 
  # elsif($mqual =~ /utrpoor/) { $mqv -= 1; }
  # $mqv= kAAQUAL_MAX if($aaw >= 300); #   is this long enough to always call best?
  # $mqv++ if($aaw >= 300); # biggies  
  # $mqv-- if($aaw < 90); # shorties
  
  ## change to Noncode for mqv < max, and fclass == hclass == Noncode?
  my $aclass=($mqv > 0)?"Code":($mqv<0)?"Noncode":"Unknown"; # too simple as yet
  
  if($fclass eq "Noncode" and $aclass ne $fclass and $fclass eq $hclass) {
    if($mqv < kAAQUAL_MAX-2 or $isutrorf or $pcds < qMINPCDS ) { $aclass=$fclass; }
  }
    #? add agreement level? 3,2,1?
  my $tclass= 
    (($aclass eq $fclass) and ($fclass eq $hclass))?$aclass:
    #n#($fclass eq $hclass)?$fclass:
    #n#($aclass eq "Code")?$aclass: # this one superceeds other two by inspection
    (($aclass eq $fclass) or ($aclass eq $hclass))?$aclass: #?? problematic classing if aclass/utrbad suspect
    ($aclass eq "Code")?$aclass: # this one superceeds other two by inspection
    ($fclass eq $hclass)?$fclass:
    "Unknown"; # no consensus
  return($tclass,$aclass);
}


sub makeCdsQual {
  my($aaseq,$aasize,$aaqualSetHash)= @_;
  my $ismrna=1; # only this, really iscds=1 : $ENV{ismrna}||$ENV{mrna}||0; 
  my $doff=1; # $ENV{doff}; 
  my $doid=1; # no oids here? drop this? option?  $ENV{doid};
  our $makeAaQual_SETHASH= $aaqualSetHash||0;
  
  # my $aasize= makename($aaseq,".cds.qual"); 
  # if( -s $aasize ) {
  #   my $aqhash= getAaQual($aasize) if($makeAaQual_SETHASH); # yes or not?
  #   return($aasize);
  # }
  # my ($ok,$inh)= openRead($aaseq);  
  # my $aasize= "$aaseq.qual";

  our($ok,$inh,$outh);
  if($aaseq =~ /\.gz/) { $ok= open($inh,"gunzip -c $aaseq |"); }
  else { $ok= open($inh,$aaseq); }
  if($ok) {
    if($aasize=~/stdout|-$/) { $outh=*STDOUT; $ok=1; }
    else { $ok= open(AAQ,'>',$aasize); $outh=*AAQ; }
  }
  return unless($ok);
  
  my($id,$seq,$aat,$aag,$al,$cl,$naa)= (0) x 9;
  my @xval;
  
  sub puta { 
    my($id,$seq,$aat,$aag,$al,$cl,@xval)= @_; 
    our ($outh,$makeAaQual_SETHASH);
    $seq=uc($seq); # dang need upcase for iscode_ tests..
    unless($aat) {
     $aat = $seq=~tr/ACGTacgt/ACGTacgt/; # nogap count 
     $aag = $seq=~tr/Nn/Nn/; # gap count
    }
    
    # $al=($aat+$aag).",na" if($al eq "na"); 
    # if($al eq "na" or not $al) 
    if($al =~ /na/ or not $al) 
    {
      my $cdsw= $aat + $aag; 
      my $aaw= int( $cdsw / 3);
      my $pcds= ($cl>0 and $cdsw>0)?int(0.5 + 100*$cdsw/$cl):0;
      my $ac=(substr($seq,0,3) eq 'ATG') ? 1 : 0;
      $ac |= 2; # need codon table to check stop codons :( .. see also below
      my $compl= ($ac==3)?"complete":($ac==2)?"partial5":($ac==1)?"partial3":"partial";
      if($cl - $cdsw <= $MINUTR) { } # ignore pcds if utr small
      elsif($pcds <= $pCDSpoor) { $compl.= ($pcds <= $pCDSbad)?"-utrbad":"-utrpoor"; }
      $al="$aaw,$pcds%,$compl";
    }
    my($fval,$fclass,$fn)= iscode_fickett($id,$seq); #fail now, 000 only 
    my($hval,$hclass,$hn)= iscode_hexamer($id,$seq);

    my $isutrorf = (grep/utrorf/,@xval)?1:0; # in oids of some data.. not a good test
    my($cpclass,$aaclass)= codingpotClass($al,$fclass,$fval,$hclass,$hval,$isutrorf); ## add new 

    $hval= sprintf("%.4g",$hval);
    my @val=($cl);
    push @val,"$cpclass,CF:$fclass/$fval,CH:$hclass/$hval,CA:$aaclass"; # ,$fn .. merge CF:,CH: in 1 col?
    #push @val,"CF:$fclass/$fval"; # ,$fn .. merge CF:,CH: in 1 col?
    #push @val,"CH:$hclass/$hval";
    push @val,@xval if(@xval);
    
    print $outh join("\t",$id,$aat,$aag,$al,@val)."\n"; 
    
    ## drop this from cdsqual, option here: set global hashes
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
      $naa+= puta($id,$seq,$aat,$aag,$al,$cl,@xval) if($id); 
      $id=$td; $aat=$aag=0; $seq=""; @xval=();
      ($al)=m/aalen=([^;\s]+)/; $al||="na";
      ($cl)=m/clen=(\d+)/; $cl||=0; 
      if($doff){ 
        my($or)=m/strand=(.)/; $or||="."; 
        my($ofs)=m/offs=([\d-]+)/; 
        push @xval, (($ofs)?"$ofs:$or":0);  
        } 
      if($doid){ my $oid;
        unless(($oid)=m/oid=([^;\s]+)/) { ($oid)=m/gene[=:]([^;\s]+)/; } 
        $oid||="noid"; 
        push @xval, $oid; # $cl.="\t$oid"; 
      } 
    } elsif(/\w/) {
       chomp(); $seq.=$_; # only for for aa.seq: s/\*$//; 
       #calc from seq now: if($ismrna) { $aat += tr/ACGTacgt/ACGTacgt/; $aag+= tr/Nn/Nn/; } # replace these from seq
       #calc from seq now: else { $aat += tr/A-WYZa-wyz/A-WYZa-wyz/; $aag += tr/Xx\*/Xx\*/; }
    }
  } 
  $naa += puta($id,$seq,$aat,$aag,$al,$cl,@xval);  
  
  close($outh); close($inh);
  # loggit(0, "makeAaQual: naa=$naa IN $aasize\n") if($DEBUG);
  return($aasize);
}

=item read_hexcode

mod3hexamer.tsv -- change format? want only hex, logratio
hexamer	coding	noncoding	ratio	logr
AAAAAA	0.000489501	0.00228415	0.99821	-0.001792
AAAAAC	0.000466459	0.000870139	0.999597	-0.000403369
AAAAAG	0.000845298	0.000780122	1.00007	6.51325e-05
AAAAAT	0.000402803	0.00119556	0.999208	-0.000792043
  nrow=4096
  
=cut

sub read_hexcode {
  my($hexfile)= @_;
  my %hexcode=();
  open(H,$hexfile); 
  while(<H>) { 
    my($hex,$pcd,$pnc,$rpc,$lrpc)=split;
    unless($lrpc) { $lrpc= log( (1+$pcd)/(1+$pnc) ); }
    $hexcode{$hex}= $lrpc;
    }
  close(H);
  return(\%hexcode);
}

sub iscode_hexamer {
  my($id, $cds)=@_;
  my $n=length($cds);
  my($val,$nx)=(0,0);
  for(my $i=0; $i<$n-6; $i+=3) {
    my $hex=substr($cds,$i,6);
    if(my $h= $hexcode->{$hex}) { $val+=$h; $nx++; }  # loglike value of coding+1 /noncode -1,    
    }

  my $ftype= ($val <= HNoncode) ? "Noncode" : ($val >= HCoding)? "Code" : "Unknown";
  return($val,$ftype,$nx);
}


sub _mina { my @v=@_; my $m=shift @v; for(@v){ $m=$_ if($_<$m); }; return $m; }
sub _maxa { my @v=@_; my $m=shift @v; for(@v){ $m=$_ if($_>$m); }; return $m; }

our (%position_prob, %position_weight, @position_para, %content_prob, %content_weight, @content_para);
BEGIN{ # dangit perl..
# Fickett TESTCODE data; NAR 10(17) 5303-531
our %position_prob =(
'A'=>[0.94,0.68,0.84,0.93,0.58,0.68,0.45,0.34,0.20,0.22],
'C'=>[0.80,0.70,0.70,0.81,0.66,0.48,0.51,0.33,0.30,0.23],
'G'=>[0.90,0.88,0.74,0.64,0.53,0.48,0.27,0.16,0.08,0.08],
'T'=>[0.97,0.97,0.91,0.68,0.69,0.44,0.54,0.20,0.09,0.09]
);
our %position_weight=('A'=>0.26,'C'=>0.18,'G'=>0.31,'T'=>0.33);
our @position_para  = (1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 0.0);

our %content_prob=(
'A'=>[0.28,0.49,0.44,0.55,0.62,0.49,0.67,0.65,0.81,0.21],
'C'=>[0.82,0.64,0.51,0.64,0.59,0.59,0.43,0.44,0.39,0.31],
'G'=>[0.40,0.54,0.47,0.64,0.64,0.73,0.41,0.41,0.33,0.29],
'T'=>[0.28,0.24,0.39,0.40,0.55,0.75,0.56,0.69,0.51,0.58]
);
our %content_weight=('A'=>0.11,'C'=>0.12,'G'=>0.15,'T'=>0.14);
our @content_para  =(0.33,0.31,0.29,0.27,0.25,0.23,0.21,0.17,0);
}

sub fickett_position_prob {
  my($value, $base, $debug)= @_;

if(0) { #DEBUG  
  # ok now with BEGIN eval..
  my ($ix)= grep { $value >= $position_para[$_] } (0..9);
  my $p= $position_prob{$base}->[$ix] || 0;
  my $w= $position_weight{$base} || 0; 
  print "#fpp: $value,$base i$ix p=$p, w=$w\n" if($debug);
  return $p * $w;
} else {  
  return 0 if($value<0);
  for my $i (0..$#position_para) {
    if ($value >= $position_para[$i]) { return $position_prob{$base}->[$i] * $position_weight{$base}; }
  }
}  
  return 0;
}

sub fickett_comp_prob {
  my($value, $base)=@_;
  return 0 if($value<0);
  for my $i (0..$#content_para) {
     if ($value >= $content_para[$i]) { return $content_prob{$base}->[$i] * $content_weight{$base}; }
  }
  return 0;
}



sub iscode_fickett {
  my($id,$cds)=@_;
  my @cds= split"",$cds;
  my %ACGT=(A=>1,C=>2,G=>3,T=>4);
  my @kACGT=qw(A C G T); # (0=>'A',1=>'C',2=>'G',3=>'T');

  my $ndebug=0;
  my $ncd=0;
  my($c,$f,$i,$j,$k,$v);
  my(@scores);
  for($i=0; $i<@cds; $i++) {
    if( $k=$ACGT{$cds[$i]} ) { $scores[$k-1][ $i % 3 ]++; $ncd++; }
  }
  my ($tsum,@bsum);
  for $k (0,1,2,3) {
    for $j (0,1,2) { my $s= $scores[$k][$j]; $bsum[$k]+=$s; $tsum+=$s; }
  }
  my @comp= map{ $_/$tsum } @bsum;
  my @posi;
  for $k (0,1,2,3) {
    my $vmax= _maxa(@{$scores[$k]});
    my $vmin= _mina(@{$scores[$k]});
    $posi[$k]= $vmax / (1+$vmin);
  }  

  my $fscore=0;
  for $k (0,1,2,3) {  
    my $base=$kACGT[$k];
    my ($pk,$ck);
    my $tdebug=0; # ($debug and $ndebug < 9)?1:0;
    my $pp= fickett_position_prob(($pk=$posi[$k]), $base, $tdebug); 
    my $cp= fickett_comp_prob(($ck=$comp[$k]), $base, $tdebug); 
    #if($debug and $ndebug++ < 9) {
    #  print "#dbg.fick.$id, $k.$base, pos:$pk=$pp, com:$ck=$cp\n";
    #}
    $fscore += $pp + $cp;
  }
  my $ftype= ($fscore < FNoncode) ? "Noncode" : ($fscore >= FCoding)? "Code" : "Unknown";
  
  return($fscore,$ftype,$ncd);
}



BEGIN { # average codon-pair stat log((1+code)/(1+noncode)) from cpat/dat/{human,zfish,fly}hexamer.tsv
our %hexcodencode= (
AAAAAA => -0.001792, AAAAAC => -0.000403369, AAAAAG => 6.51325e-05, AAAAAT => -0.000792043, AAAACA => -0.000679693, 
AAAACC => -4.12391e-05, AAAACG => -0.000153407, AAAACT => -0.000328468, AAAAGA => -0.000418644, AAAAGC => -0.000171301, 
AAAAGG => -0.000207088, AAAAGT => -0.000354192, AAAATA => -0.000892698, AAAATC => -0.000127767, AAAATG => -0.000309964, 
AAAATT => -0.000543834, AAACAA => -0.000657578, AAACAC => -0.000228344, AAACAG => 0.000102252, AAACAT => -0.000393765, 
AAACCA => -7.101e-05, AAACCC => 9.78338e-05, AAACCG => 1.56289e-05, AAACCT => -3.0542e-06, AAACGA => -9.51859e-05, 
AAACGC => 4.98873e-05, AAACGG => 3.72378e-06, AAACGT => -9.84276e-05, AAACTA => -0.000219944, AAACTC => 6.36571e-05, 
AAACTG => 0.000256327, AAACTT => -0.000186055, AAAGAA => -3.68559e-05, AAAGAC => 0.000346301, AAAGAG => 0.000570842, 
AAAGAT => 0.000192557, AAAGCA => -0.000119841, AAAGCC => 0.000234036, AAAGCG => -2.51453e-05, AAAGCT => 5.50655e-05, 
AAAGGA => 5.3764e-05, AAAGGC => 0.000123013, AAAGGG => -7.45839e-05, AAAGGT => 1.07997e-05, AAAGTA => -0.00020482, 
AAAGTC => 1.85231e-05, AAAGTG => 0.000134443, AAAGTT => -0.000193516, AAATAA => -0.00121122, AAATAC => -8.63985e-05, 
AAATAG => -0.000363919, AAATAT => -0.000658324, AAATCA => -0.000324167, AAATCC => 2.41329e-05, AAATCG => -8.13856e-05, 
AAATCT => -0.000176975, AAATGA => -0.000639329, AAATGC => -0.000200671, AAATGG => -0.000164389, AAATGT => -0.000556736, 
AAATTA => -0.000664623, AAATTC => -0.000117341, AAATTG => -0.000264469, AAATTT => -0.000542744, AACAAA => -0.000345911, 
AACAAC => 0.000482256, AACAAG => 0.000588018, AACAAT => 6.8214e-05, AACACA => -0.000113513, AACACC => 0.000243386, 
AACACG => 0.00010869, AACACT => -2.53991e-05, AACAGA => -0.000169587, AACAGC => 0.00037576, AACAGG => -1.87326e-05, 
AACAGT => 5.02378e-05, AACATA => -0.000136357, AACATC => 0.000523415, AACATG => 0.000320944, AACATT => -9.33998e-05, 
AACCAA => -0.000130475, AACCAC => 0.000111641, AACCAG => 0.000433559, AACCAT => -4.2085e-05, AACCCA => 2.98097e-05, 
AACCCC => 0.000173508, AACCCG => 6.18311e-05, AACCCT => 6.99701e-05, AACCGA => 1.86601e-05, AACCGC => 0.000180407, 
AACCGG => 9.30835e-05, AACCGT => 3.97605e-05, AACCTA => -1.21928e-05, AACCTC => 0.000213068, AACCTG => 0.000589102, 
AACCTT => -9.79311e-06, AACGAA => 5.58091e-05, AACGAC => 0.000281988, AACGAG => 0.000499462, AACGAT => 0.000234427, 
AACGCA => -1.30733e-05, AACGCC => 0.000277441, AACGCG => 1.76053e-05, AACGCT => 7.27397e-05, AACGGA => 0.000193937, 
AACGGC => 0.000289311, AACGGG => 8.59625e-05, AACGGT => 8.22511e-05, AACGTA => -4.63083e-05, AACGTC => 0.00010802, 
AACGTG => 0.000244386, AACGTT => -5.19572e-05, AACTAA => -0.00037507, AACTAC => 0.000347453, AACTAG => -0.000169084, 
AACTAT => 3.39023e-05, AACTCA => -6.50413e-05, AACTCC => 0.000246574, AACTCG => 0.000154714, AACTCT => 1.41511e-05, 
AACTGA => -0.000403368, AACTGC => 9.70073e-05, AACTGG => 0.000104724, AACTGT => -9.58292e-05, AACTTA => -0.000216277, 
AACTTC => 0.000351034, AACTTG => 1.25416e-05, AACTTT => -3.27902e-05, AAGAAA => 0.000341338, AAGAAC => 0.000615088, 
AAGAAG => 0.00151823, AAGAAT => 0.000132819, AAGACA => 8.2981e-07, AAGACC => 0.000452398, AAGACG => 0.000225469, 
AAGACT => 0.000111169, AAGAGA => -2.22786e-05, AAGAGC => 0.000354331, AAGAGG => 0.000155664, AAGAGT => 9.02773e-05, 
AAGATA => -6.79285e-05, AAGATC => 0.000676026, AAGATG => 0.000465837, AAGATT => 0.000171341, AAGCAA => -5.49358e-05, 
AAGCAC => 0.000221675, AAGCAG => 0.000715833, AAGCAT => -1.62093e-05, AAGCCA => 0.000132368, AAGCCC => 0.000367553, 
AAGCCG => 0.000230031, AAGCCT => 0.000155475, AAGCGA => 0.000103833, AAGCGC => 0.00033696, AAGCGG => 0.000212628, 
AAGCGT => 0.000130379, AAGCTA => 3.29744e-05, AAGCTC => 0.000252751, AAGCTG => 0.000957191, AAGCTT => 6.07675e-05, 
AAGGAA => 0.000487179, AAGGAC => 0.000735124, AAGGAG => 0.00151314, AAGGAT => 0.000621452, AAGGCA => 0.000137531, 
AAGGCC => 0.000576149, AAGGCG => 0.000234049, AAGGCT => 0.000347217, AAGGGA => 0.000152827, AAGGGC => 0.000421101, 
AAGGGG => 2.19883e-05, AAGGGT => 0.000195664, AAGGTA => -1.39826e-06, AAGGTC => 0.000276675, AAGGTG => 0.000671716, 
AAGGTT => 9.47386e-05, AAGTAA => -0.000359124, AAGTAC => 0.000428054, AAGTAG => -0.000186818, AAGTAT => 3.68549e-05, 
AAGTCA => -6.59604e-05, AAGTCC => 0.000266966, AAGTCG => 0.000159835, AAGTCT => 7.29119e-05, AAGTGA => -0.000355903, 
AAGTGC => 0.000113112, AAGTGG => 5.83962e-05, AAGTGT => -0.000100127, AAGTTA => -0.00017872, AAGTTC => 0.000472912, 
AAGTTG => 6.75091e-05, AAGTTT => 2.90168e-05, AATAAA => -0.0010012, AATAAC => -3.40026e-05, AATAAG => 1.31342e-05, 
AATAAT => -0.000480847, AATACA => -0.000319435, AATACC => 1.40456e-05, AATACG => -6.6721e-06, AATACT => -0.000172979, 
AATAGA => -0.000197812, AATAGC => -2.99809e-05, AATAGG => -0.000113792, AATAGT => -0.000133756, AATATA => -0.000472595, 
AATATC => -4.56661e-05, AATATG => -8.40714e-05, AATATT => -0.00054377, AATCAA => -0.000246964, AATCAC => -3.62761e-05, 
AATCAG => 0.000151426, AATCAT => -0.000219756, AATCCA => -6.44179e-06, AATCCC => 0.000178658, AATCCG => 0.000117374, 
AATCCT => 5.21977e-05, AATCGA => -6.69585e-05, AATCGC => 9.37876e-05, AATCGG => 8.9461e-06, AATCGT => -2.34491e-05, 
AATCTA => -9.28307e-05, AATCTC => 3.01312e-05, AATCTG => 0.000287011, AATCTT => -0.00010194, AATGAA => 3.24777e-05, 
AATGAC => 0.00027331, AATGAG => 0.000450746, AATGAT => 7.58247e-05, AATGCA => -6.39009e-06, AATGCC => 0.000440653, 
AATGCG => 8.03806e-05, AATGCT => 5.11349e-05, AATGGA => 0.000235456, AATGGC => 0.000376738, AATGGG => 6.29088e-05, 
AATGGT => 8.95189e-05, AATGTA => -0.000248139, AATGTC => 5.09047e-05, AATGTG => 0.000269032, AATGTT => -0.000271452, 
AATTAA => -0.00078489, AATTAC => -9.40919e-05, AATTAG => -0.00028659, AATTAT => -0.000460772, AATTCA => -0.000246197, 
AATTCC => -3.92827e-05, AATTCG => -1.0601e-05, AATTCT => -0.0002087, AATTGA => -0.000392659, AATTGC => -0.000138305, 
AATTGG => -0.00015223, AATTGT => -0.000374697, AATTTA => -0.000568338, AATTTC => -0.000197362, AATTTG => -0.000227357, 
AATTTT => -0.000672928, ACAAAA => -0.000627524, ACAAAC => -0.000238352, ACAAAG => -9.12142e-05, ACAAAT => -0.000359372, 
ACAACA => -0.000285519, ACAACC => -1.7555e-05, ACAACG => 4.70775e-06, ACAACT => -0.000139129, ACAAGA => -0.00021941, 
ACAAGC => -9.64185e-05, ACAAGG => -0.000113957, ACAAGT => -0.000119293, ACAATA => -0.000275333, ACAATC => -7.60247e-05, 
ACAATG => -7.52734e-05, ACAATT => -0.000246231, ACACAA => -0.000301282, ACACAC => -0.00077191, ACACAG => -7.58208e-06, 
ACACAT => -0.000273379, ACACCA => 7.07658e-06, ACACCC => 6.04589e-05, ACACCG => 6.90201e-05, ACACCT => -4.03591e-05, 
ACACGA => -4.96748e-05, ACACGC => -2.66212e-05, ACACGG => -2.26986e-05, ACACGT => -5.11818e-05, ACACTA => -0.000109802, 
ACACTC => -4.31659e-05, ACACTG => 6.5593e-05, ACACTT => -0.000200163, ACAGAA => -1.97934e-05, ACAGAC => 0.000151922, 
ACAGAG => 0.000325465, ACAGAT => 7.90016e-05, ACAGCA => -0.00010102, ACAGCC => 0.000138456, ACAGCG => 1.14322e-05, 
ACAGCT => 1.24976e-05, ACAGGA => 5.02974e-05, ACAGGC => 6.1108e-05, ACAGGG => -4.81957e-06, ACAGGT => -8.47915e-07, 
ACAGTA => -0.000103936, ACAGTC => 2.70316e-06, ACAGTG => 0.000128321, ACAGTT => -0.000112726, ACATAA => -0.000394971, 
ACATAC => -0.000116758, ACATAG => -0.00016257, ACATAT => -0.000300268, ACATCA => -0.000152591, ACATCC => 2.25304e-06, 
ACATCG => 2.95217e-05, ACATCT => -5.86605e-05, ACATGA => -0.000304458, ACATGC => -0.000140864, ACATGG => -8.55596e-05, 
ACATGT => -0.00021229, ACATTA => -0.000259346, ACATTC => -8.28704e-05, ACATTG => -0.000120627, ACATTT => -0.00050758, 
ACCAAA => 0.000115971, ACCAAC => 0.000350437, ACCAAG => 0.000626839, ACCAAT => 0.000199571, ACCACA => 0.000144405, 
ACCACC => 0.000485351, ACCACG => 0.000157723, ACCACT => 0.000151768, ACCAGA => -9.62676e-05, ACCAGC => 0.00033296, 
ACCAGG => 3.41386e-06, ACCAGT => 0.000145597, ACCATA => 1.35154e-05, ACCATC => 0.000543384, ACCATG => 0.000316201, 
ACCATT => 0.000200983, ACCCAA => -6.58477e-05, ACCCAC => 3.80057e-05, ACCCAG => 0.000240954, ACCCAT => -2.16944e-05, 
ACCCCA => 1.54622e-05, ACCCCC => 1.93212e-05, ACCCCG => 5.43563e-05, ACCCCT => 4.08756e-05, ACCCGA => 1.01199e-05, 
ACCCGC => 9.48005e-05, ACCCGG => 4.70201e-05, ACCCGT => 2.43469e-05, ACCCTA => -1.6996e-05, ACCCTC => 0.0001048, 
ACCCTG => 0.000367477, ACCCTT => -3.86114e-05, ACCGAA => 7.93476e-05, ACCGAC => 0.000160072, ACCGAG => 0.000294769, 
ACCGAT => 0.000167555, ACCGCA => 2.61344e-05, ACCGCC => 0.000243806, ACCGCG => 6.58909e-06, ACCGCT => 4.62381e-05, 
ACCGGA => 0.000100923, ACCGGC => 0.000183514, ACCGGG => 3.39327e-05, ACCGGT => 6.29577e-05, ACCGTA => 2.51739e-05, 
ACCGTC => 0.000153217, ACCGTG => 0.000235526, ACCGTT => 3.15495e-05, ACCTAA => -0.00016974, ACCTAC => 0.000290934, 
ACCTAG => -0.000100179, ACCTAT => 8.90984e-05, ACCTCA => 1.83814e-05, ACCTCC => 0.000180223, ACCTCG => 0.000125076, 
ACCTCT => -8.61168e-06, ACCTGA => -0.000255886, ACCTGC => 0.000158812, ACCTGG => 4.03339e-05, ACCTGT => 2.12326e-05, 
ACCTTA => -8.24342e-05, ACCTTC => 0.000376336, ACCTTG => 3.71011e-05, ACCTTT => 5.32419e-06, ACGAAA => -0.000136355, 
ACGAAC => 1.45158e-05, ACGAAG => 3.50814e-05, ACGAAT => -3.68697e-05, ACGACA => -2.87117e-05, ACGACC => 1.4055e-05, 
ACGACG => 3.50257e-05, ACGACT => -2.0332e-05, ACGAGA => -7.48104e-05, ACGAGC => 1.63604e-05, ACGAGG => -4.89141e-05, 
ACGAGT => -4.20485e-05, ACGATA => -3.95565e-05, ACGATC => 1.50024e-05, ACGATG => 2.18581e-05, ACGATT => -5.08624e-05, 
ACGCAA => -2.95647e-05, ACGCAC => 1.01324e-05, ACGCAG => 0.000131431, ACGCAT => -3.65331e-05, ACGCCA => 8.98531e-05, 
ACGCCC => 0.000219863, ACGCCG => 0.000122345, ACGCCT => 3.93488e-05, ACGCGA => -3.411e-06, ACGCGC => 4.71265e-05, 
ACGCGG => 1.65406e-05, ACGCGT => 5.84341e-07, ACGCTA => 1.14042e-05, ACGCTC => 0.000108246, ACGCTG => 0.000321431, 
ACGCTT => -4.08218e-06, ACGGAA => 7.37309e-05, ACGGAC => 0.000191595, ACGGAG => 0.000365439, ACGGAT => 0.000157741, 
ACGGCA => 4.76792e-05, ACGGCC => 0.000226374, ACGGCG => 0.00011596, ACGGCT => 0.000102799, ACGGGA => 8.91109e-05, 
ACGGGC => 0.000193873, ACGGGG => 1.82653e-05, ACGGGT => 5.48352e-05, ACGGTA => -2.24275e-05, ACGGTC => 7.44574e-05, 
ACGGTG => 0.000262814, ACGGTT => 3.27148e-05, ACGTAA => -0.00010865, ACGTAC => 6.37655e-05, ACGTAG => -6.53236e-05, 
ACGTAT => -2.1624e-05, ACGTCA => -2.67305e-05, ACGTCC => 5.42654e-05, ACGTCG => 2.98821e-05, ACGTCT => -1.12199e-05, 
ACGTGA => -0.000123144, ACGTGC => 4.23535e-07, ACGTGG => -2.15627e-05, ACGTGT => -7.42601e-05, ACGTTA => -6.61455e-05, 
ACGTTC => 4.66399e-05, ACGTTG => -5.66458e-06, ACGTTT => -7.41618e-05, ACTAAA => -0.000155874, ACTAAC => -4.9798e-05, 
ACTAAG => 2.32672e-06, ACTAAT => -0.00012117, ACTACA => -9.70977e-05, ACTACC => 1.41842e-05, ACTACG => -2.56624e-06, 
ACTACT => -7.29701e-05, ACTAGA => -0.000110095, ACTAGC => -2.06871e-05, ACTAGG => -6.77459e-05, ACTAGT => -5.09481e-05, 
ACTATA => -0.000154306, ACTATC => 9.04492e-07, ACTATG => -3.59488e-05, ACTATT => -0.000118705, ACTCAA => -0.000127701, 
ACTCAC => -5.56756e-05, ACTCAG => 0.000115359, ACTCAT => -9.54106e-05, ACTCCA => 4.19969e-05, ACTCCC => 4.15986e-05, 
ACTCCG => 8.36756e-05, ACTCCT => 6.26565e-05, ACTCGA => -3.27393e-05, ACTCGC => 6.29985e-06, ACTCGG => 1.26736e-05, 
ACTCGT => -1.22399e-05, ACTCTA => -5.45235e-05, ACTCTC => -3.62489e-05, ACTCTG => 0.00019222, ACTCTT => -9.32151e-05, 
ACTGAA => -3.77652e-06, ACTGAC => 8.37719e-05, ACTGAG => 0.000157609, ACTGAT => 6.55127e-05, ACTGCA => -6.12053e-05, 
ACTGCC => 0.000161274, ACTGCG => -4.88468e-06, ACTGCT => 3.19937e-05, ACTGGA => 9.77626e-05, ACTGGC => 0.00012235, 
ACTGGG => 1.15827e-05, ACTGGT => 5.41386e-05, ACTGTA => -0.000135735, ACTGTC => 4.6301e-05, ACTGTG => 0.000173433, 
ACTGTT => -9.58246e-05, ACTTAA => -0.000393015, ACTTAC => -1.48264e-05, ACTTAG => -0.000174199, ACTTAT => -0.0001512, 
ACTTCA => -0.000143756, ACTTCC => -3.28334e-05, ACTTCG => -1.94715e-05, ACTTCT => -0.000129933, ACTTGA => -0.000294947, 
ACTTGC => -8.22734e-05, ACTTGG => -0.000152518, ACTTGT => -0.000202201, ACTTTA => -0.000266511, ACTTTC => -6.55027e-05, 
ACTTTG => -0.000127529, ACTTTT => -0.00043392, AGAAAA => -0.000623279, AGAAAC => -0.000178754, AGAAAG => -0.000176775, 
AGAAAT => -0.000353043, AGAACA => -0.000260544, AGAACC => -4.85123e-05, AGAACG => -7.50749e-05, AGAACT => -0.000153548, 
AGAAGA => -0.000333199, AGAAGC => -0.00018175, AGAAGG => -0.000205808, AGAAGT => -0.000170649, AGAATA => -0.00024129, 
AGAATC => -8.13757e-05, AGAATG => -0.000143449, AGAATT => -0.000243041, AGACAA => -0.000187789, AGACAC => -0.000128895, 
AGACAG => -2.89102e-05, AGACAT => -0.000165036, AGACCA => -9.30433e-05, AGACCC => -3.85683e-06, AGACCG => -2.23505e-05, 
AGACCT => -4.31862e-05, AGACGA => -5.80657e-05, AGACGC => -5.81531e-06, AGACGG => -3.30492e-05, AGACGT => -4.29465e-05, 
AGACTA => -8.19503e-05, AGACTC => -5.41979e-05, AGACTG => -3.56853e-05, AGACTT => -0.000156667, AGAGAA => -0.000168991, 
AGAGAC => 6.65167e-05, AGAGAG => -7.87718e-05, AGAGAT => -1.14671e-05, AGAGCA => -0.000180569, AGAGCC => -3.62514e-06, 
AGAGCG => -9.0237e-05, AGAGCT => -0.000101978, AGAGGA => -0.000100301, AGAGGC => -6.47925e-05, AGAGGG => -0.000154603, 
AGAGGT => -7.79076e-05, AGAGTA => -0.000115087, AGAGTC => -5.77063e-05, AGAGTG => -3.57581e-05, AGAGTT => -0.000145062, 
AGATAA => -0.000300424, AGATAC => -4.07949e-05, AGATAG => -0.00015873, AGATAT => -0.000197391, AGATCA => -0.000154371, 
AGATCC => -7.57069e-05, AGATCG => -5.36365e-05, AGATCT => -0.000144247, AGATGA => -0.000359529, AGATGC => -0.000190384, 
AGATGG => -0.000238309, AGATGT => -0.000243954, AGATTA => -0.000202305, AGATTC => -0.000114603, AGATTG => -0.000160189, 
AGATTT => -0.000332365, AGCAAA => 7.50724e-07, AGCAAC => 0.000313104, AGCAAG => 0.000421068, AGCAAT => 0.000141963, 
AGCACA => 4.09781e-05, AGCACC => 0.000338976, AGCACG => 8.02404e-05, AGCACT => 7.76577e-05, AGCAGA => -0.000188318, 
AGCAGC => 0.000536891, AGCAGG => -0.000100609, AGCAGT => 0.000197123, AGCATA => -0.000100983, AGCATC => 0.00030411, 
AGCATG => 0.000212209, AGCATT => -1.83775e-05, AGCCAA => -8.11719e-05, AGCCAC => 3.7595e-05, AGCCAG => 0.000307788, 
AGCCAT => -0.000110864, AGCCCA => 5.19418e-06, AGCCCC => 8.95571e-05, AGCCCG => 1.89544e-05, AGCCCT => 2.86208e-05, 
AGCCGA => 4.45911e-06, AGCCGC => 9.42738e-05, AGCCGG => 5.69465e-05, AGCCGT => -1.3979e-06, AGCCTA => -3.56493e-05, 
AGCCTC => 7.84555e-05, AGCCTG => 0.000303914, AGCCTT => -4.55916e-05, AGCGAA => 8.45131e-05, AGCGAC => 0.000272333, 
AGCGAG => 0.000356578, AGCGAT => 0.000208526, AGCGCA => -5.70264e-06, AGCGCC => 0.000202451, AGCGCG => -5.23837e-06, 
AGCGCT => 2.82019e-05, AGCGGA => 9.44981e-05, AGCGGC => 0.000255208, AGCGGG => 3.4722e-05, AGCGGT => 6.24852e-05, 
AGCGTA => -1.91106e-05, AGCGTC => 5.4968e-05, AGCGTG => 0.000135828, AGCGTT => 5.57298e-06, AGCTAA => -0.000244559, 
AGCTAC => 0.000195686, AGCTAG => -0.000164544, AGCTAT => 3.95694e-05, AGCTCA => -4.8138e-05, AGCTCC => 0.000265596, 
AGCTCG => 7.66443e-05, AGCTCT => 2.30764e-05, AGCTGA => -0.000392598, AGCTGC => -3.30934e-05, AGCTGG => -8.71395e-05, 
AGCTGT => -0.000153286, AGCTTA => -0.000126795, AGCTTC => 0.000218158, AGCTTG => 6.35035e-06, AGCTTT => -5.71981e-05, 
AGGAAA => -0.000146416, AGGAAC => 6.44307e-05, AGGAAG => 0.000121868, AGGAAT => -6.96034e-05, AGGACA => -0.000170619, 
AGGACC => -1.8793e-05, AGGACG => -2.66172e-05, AGGACT => -9.49631e-05, AGGAGA => -0.000190619, AGGAGC => -7.91468e-05, 
AGGAGG => -0.000189635, AGGAGT => -8.73049e-05, AGGATA => -8.73936e-05, AGGATC => 2.17722e-05, AGGATG => -3.72754e-05, 
AGGATT => -8.50491e-05, AGGCAA => -0.000134518, AGGCAC => -4.57499e-05, AGGCAG => -3.38975e-05, AGGCAT => -0.000107865, 
AGGCCA => -9.9681e-05, AGGCCC => -1.87789e-05, AGGCCG => 9.81149e-07, AGGCCT => -6.89106e-05, AGGCGA => -4.64129e-05, 
AGGCGC => 2.55255e-05, AGGCGG => -2.27416e-05, AGGCGT => -3.22238e-05, AGGCTA => -8.94324e-05, AGGCTC => -4.33441e-05, 
AGGCTG => 3.62676e-06, AGGCTT => -0.000135888, AGGGAA => -3.39638e-05, AGGGAC => 6.76451e-05, AGGGAG => 0.000131385, 
AGGGAT => 5.63133e-05, AGGGCA => -0.000104518, AGGGCC => 2.12954e-05, AGGGCG => -3.65799e-06, AGGGCT => -5.68721e-05, 
AGGGGA => -0.000104501, AGGGGC => -1.4064e-05, AGGGGG => -0.000109105, AGGGGT => -5.7394e-05, AGGGTA => -6.68826e-05, 
AGGGTC => -2.38209e-05, AGGGTG => 1.54098e-06, AGGGTT => -0.000109853, AGGTAA => -0.000179123, AGGTAC => -1.55746e-05, 
AGGTAG => -0.000156058, AGGTAT => -6.13451e-05, AGGTCA => -0.000133581, AGGTCC => -4.42055e-05, AGGTCG => -3.11843e-05, 
AGGTCT => -9.00337e-05, AGGTGA => -0.000227545, AGGTGC => -9.92764e-05, AGGTGG => -0.000156534, AGGTGT => -0.000180889, 
AGGTTA => -0.00010826, AGGTTC => -3.30062e-05, AGGTTG => -9.50087e-05, AGGTTT => -0.000135504, AGTAAA => -0.000182878, 
AGTAAC => 4.98148e-06, AGTAAG => -1.51795e-05, AGTAAT => -9.38945e-05, AGTACA => -0.000104686, AGTACC => 6.33347e-06, 
AGTACG => -1.8086e-05, AGTACT => -8.84418e-05, AGTAGA => -0.000128155, AGTAGC => 2.34821e-05, AGTAGG => -9.23872e-05, 
AGTAGT => -5.36592e-05, AGTATA => -0.000147859, AGTATC => -2.93059e-05, AGTATG => -3.80717e-05, AGTATT => -0.000195069, 
AGTCAA => -0.000138481, AGTCAC => -5.3666e-05, AGTCAG => 9.6814e-05, AGTCAT => -0.00010696, AGTCCA => 1.45978e-05, 
AGTCCC => 0.000114934, AGTCCG => 8.55754e-05, AGTCCT => 5.13174e-05, AGTCGA => -1.16463e-05, AGTCGC => 6.8818e-05, 
AGTCGG => -3.27613e-06, AGTCGT => 1.21089e-05, AGTCTA => -6.85313e-05, AGTCTC => -3.88274e-05, AGTCTG => 0.000152575, 
AGTCTT => -9.62775e-05, AGTGAA => 2.05647e-05, AGTGAC => 0.000222856, AGTGAG => 0.00024382, AGTGAT => 0.000148858, 
AGTGCA => -3.27795e-05, AGTGCC => 0.000286001, AGTGCG => 4.99998e-06, AGTGCT => 1.96126e-05, AGTGGA => 9.68193e-05, 
AGTGGC => 0.000209187, AGTGGG => 7.84417e-06, AGTGGT => 4.18192e-05, AGTGTA => -0.000140708, AGTGTC => 3.2395e-05, 
AGTGTG => 5.90453e-05, AGTGTT => -0.000185042, AGTTAA => -0.000336699, AGTTAC => -1.62927e-05, AGTTAG => -0.000165979, 
AGTTAT => -0.000165481, AGTTCA => -0.000114593, AGTTCC => -6.39091e-06, AGTTCG => 1.43236e-05, AGTTCT => -9.44179e-05, 
AGTTGA => -0.000311268, AGTTGC => -0.000125863, AGTTGG => -0.000151008, AGTTGT => -0.000208261, AGTTTA => -0.000269478, 
AGTTTC => -9.80031e-05, AGTTTG => -0.000136448, AGTTTT => -0.000443637, ATAAAA => -0.000886435, ATAAAC => -0.000264927, 
ATAAAG => -0.000186786, ATAAAT => -0.000704642, ATAACA => -0.000242406, ATAACC => -3.85765e-05, ATAACG => -3.36233e-05, 
ATAACT => -0.000169212, ATAAGA => -0.000202092, ATAAGC => -8.54943e-05, ATAAGG => -0.000110779, ATAAGT => -0.000181276, 
ATAATA => -0.000571351, ATAATC => -0.000105858, ATAATG => -0.000190453, ATAATT => -0.000446384, ATACAA => -0.000305069, 
ATACAC => -0.000229127, ATACAG => -4.50174e-05, ATACAT => -0.00041364, ATACCA => -5.7882e-05, ATACCC => 4.56575e-05, 
ATACCG => 1.56439e-05, ATACCT => -6.49919e-05, ATACGA => -5.03367e-05, ATACGC => 3.75561e-05, ATACGG => -2.6032e-05, 
ATACGT => -6.07081e-05, ATACTA => -0.000113869, ATACTC => -2.79438e-05, ATACTG => -1.98737e-05, ATACTT => -0.00021102, 
ATAGAA => -9.27511e-05, ATAGAC => 4.76664e-05, ATAGAG => 5.93165e-05, ATAGAT => -3.74122e-05, ATAGCA => -7.5633e-05, 
ATAGCC => 6.64803e-05, ATAGCG => 1.34342e-05, ATAGCT => -6.82658e-05, ATAGGA => -5.82569e-05, ATAGGC => -1.54315e-05, 
ATAGGG => -5.19939e-05, ATAGGT => -7.81999e-05, ATAGTA => -0.000133983, ATAGTC => -3.75014e-05, ATAGTG => 1.80986e-05, 
ATAGTT => -0.000175115, ATATAA => -0.00060653, ATATAC => -0.000254446, ATATAG => -0.000263693, ATATAT => -0.000928164, 
ATATCA => -0.000177887, ATATCC => -3.71763e-05, ATATCG => -1.91985e-05, ATATCT => -0.00020289, ATATGA => -0.000374818, 
ATATGC => -0.000159945, ATATGG => -0.000139625, ATATGT => -0.000368546, ATATTA => -0.00044117, ATATTC => -0.000155007, 
ATATTG => -0.000271092, ATATTT => -0.000780754, ATCAAA => 2.86071e-05, ATCAAC => 0.000495195, ATCAAG => 0.000676097, 
ATCAAT => 0.000150793, ATCACA => 0.000112592, ATCACC => 0.000412479, ATCACG => 0.00012713, ATCACT => 0.000125809, 
ATCAGA => -0.000116756, ATCAGC => 0.000310652, ATCAGG => -6.3755e-06, ATCAGT => 8.12715e-05, ATCATA => -5.46999e-05, 
ATCATC => 0.000485842, ATCATG => 0.000301269, ATCATT => 8.00897e-05, ATCCAA => -4.25602e-05, ATCCAC => 0.000188253, 
ATCCAG => 0.000520742, ATCCAT => -3.71842e-05, ATCCCA => 3.52466e-05, ATCCCC => 8.59007e-05, ATCCCG => 4.8646e-05, 
ATCCCT => 5.41873e-05, ATCCGA => 6.38088e-05, ATCCGC => 0.000205322, ATCCGG => 0.000123024, ATCCGT => 8.06539e-05, 
ATCCTA => 4.82362e-06, ATCCTC => 0.000212613, ATCCTG => 0.000624357, ATCCTT => 5.32162e-06, ATCGAA => 0.000122661, 
ATCGAC => 0.000343347, ATCGAG => 0.000623966, ATCGAT => 0.000294399, ATCGCA => 4.09955e-05, ATCGCC => 0.000327958, 
ATCGCG => 2.2029e-05, ATCGCT => 0.000126892, ATCGGA => 0.000104307, ATCGGC => 0.000223161, ATCGGG => 6.8498e-05, 
ATCGGT => 6.42189e-05, ATCGTA => -3.0141e-06, ATCGTC => 0.000162151, ATCGTG => 0.000280743, ATCGTT => 1.23377e-05, 
ATCTAA => -0.00025551, ATCTAC => 0.000340367, ATCTAG => -0.000129806, ATCTAT => 1.07275e-05, ATCTCA => -3.49001e-05, 
ATCTCC => 0.000203879, ATCTCG => 0.000115488, ATCTCT => 6.20058e-06, ATCTGA => -0.000306486, ATCTGC => 0.000119989, 
ATCTGG => 3.99685e-05, ATCTGT => -8.43784e-05, ATCTTA => -0.000138225, ATCTTC => 0.00037279, ATCTTG => 1.67519e-05, 
ATCTTT => -5.43122e-05, ATGAAA => -0.000135349, ATGAAC => 0.000332401, ATGAAG => 0.000571625, ATGAAT => -4.38031e-05, 
ATGACA => -3.23719e-05, ATGACC => 0.000289386, ATGACG => 8.99586e-05, ATGACT => 1.29295e-06, ATGAGA => -8.17302e-05, 
ATGAGC => 0.000228817, ATGAGG => 6.10319e-05, ATGAGT => -3.89939e-07, ATGATA => -0.000139556, ATGATC => 0.000301332, 
ATGATG => 0.00026102, ATGATT => -8.88129e-05, ATGCAA => -0.000115535, ATGCAC => 9.6709e-05, ATGCAG => 0.000482524, 
ATGCAT => -0.000152242, ATGCCA => 3.30139e-05, ATGCCC => 0.000234901, ATGCCG => 0.00013828, ATGCCT => 7.62863e-05, 
ATGCGA => 3.8808e-05, ATGCGC => 0.000180861, ATGCGG => 8.42071e-05, ATGCGT => 4.80001e-05, ATGCTA => -3.6628e-05, 
ATGCTC => 0.000125548, ATGCTG => 0.000555371, ATGCTT => -8.41814e-05, ATGGAA => 0.000217024, ATGGAC => 0.000574446, 
ATGGAG => 0.000960225, ATGGAT => 0.000348092, ATGGCA => 0.000142089, ATGGCC => 0.000527468, ATGGCG => 0.000217471, 
ATGGCT => 0.000259719, ATGGGA => 0.000127176, ATGGGC => 0.000383479, ATGGGG => 1.19841e-05, ATGGGT => 0.000127164, 
ATGGTA => -2.96423e-05, ATGGTC => 0.000159483, ATGGTG => 0.000442682, ATGGTT => -1.58575e-05, ATGTAA => -0.000402806, 
ATGTAC => 0.000201434, ATGTAG => -0.000229424, ATGTAT => -0.000258466, ATGTCA => -5.84022e-05, ATGTCC => 0.000223506, 
ATGTCG => 0.000116867, ATGTCT => -2.8073e-08, ATGTGA => -0.000355041, ATGTGC => -4.30666e-05, ATGTGG => -2.84231e-05, 
ATGTGT => -0.000283321, ATGTTA => -0.000241779, ATGTTC => 0.000118903, ATGTTG => -3.11259e-05, ATGTTT => -0.000263775, 
ATTAAA => -0.000499369, ATTAAC => -8.02839e-05, ATTAAG => 5.28717e-05, ATTAAT => -0.000371223, ATTACA => -0.00019783, 
ATTACC => -3.21166e-05, ATTACG => -1.95474e-05, ATTACT => -0.000147607, ATTAGA => -0.000183971, ATTAGC => -6.23895e-05, 
ATTAGG => -8.50433e-05, ATTAGT => -0.000162101, ATTATA => -0.000442408, ATTATC => 1.48651e-06, ATTATG => -0.000120734, 
ATTATT => -0.000549957, ATTCAA => -0.000204405, ATTCAC => -1.79656e-05, ATTCAG => 0.000185965, ATTCAT => -0.00017094, 
ATTCCA => -4.20359e-05, ATTCCC => 4.8309e-05, ATTCCG => 5.16953e-05, ATTCCT => 2.81999e-05, ATTCGA => -1.57399e-05, 
ATTCGC => 9.26699e-05, ATTCGG => 2.31915e-05, ATTCGT => 1.5317e-05, ATTCTA => -0.000121853, ATTCTC => -2.99997e-05, 
ATTCTG => 0.000180241, ATTCTT => -0.000169748, ATTGAA => 9.22352e-05, ATTGAC => 0.00026839, ATTGAG => 0.000460362, 
ATTGAT => 0.000143331, ATTGCA => -2.27044e-05, ATTGCC => 0.00031564, ATTGCG => 4.4727e-05, ATTGCT => 9.0019e-05, 
ATTGGA => 0.000134429, ATTGGC => 0.000214457, ATTGGG => 6.98486e-06, ATTGGT => 4.9666e-05, ATTGTA => -0.000227291, 
ATTGTC => 9.2808e-05, ATTGTG => 0.000230999, ATTGTT => -0.000231136, ATTTAA => -0.000823887, ATTTAC => -0.000174034, 
ATTTAG => -0.000345496, ATTTAT => -0.000681583, ATTTCA => -0.000343364, ATTTCC => -0.000125704, ATTTCG => -6.5602e-05, 
ATTTCT => -0.000345898, ATTTGA => -0.000481644, ATTTGC => -0.000215854, ATTTGG => -0.000209763, ATTTGT => -0.000495766, 
ATTTTA => -0.000841808, ATTTTC => -0.000352415, ATTTTG => -0.000383953, ATTTTT => -0.00107457, CAAAAA => -0.000540161, 
CAAAAC => -0.000196169, CAAAAG => -1.36768e-05, CAAAAT => -0.000329413, CAAACA => -0.000394418, CAAACC => -6.16688e-05, 
CAAACG => -4.2671e-05, CAAACT => -0.000196759, CAAAGA => -0.00026722, CAAAGC => -0.000121677, CAAAGG => -0.000129184, 
CAAAGT => -0.000151326, CAAATA => -0.000383694, CAAATC => -0.000111545, CAAATG => -9.57977e-05, CAAATT => -0.000315391, 
CAACAA => -0.000175061, CAACAC => -9.24905e-05, CAACAG => 0.000396422, CAACAT => -0.000169788, CAACCA => -7.42989e-05, 
CAACCC => -1.05124e-05, CAACCG => 4.17716e-05, CAACCT => -4.89019e-05, CAACGA => -9.9194e-05, CAACGC => 3.01436e-06, 
CAACGG => -9.20128e-06, CAACGT => -3.22859e-05, CAACTA => -7.75221e-05, CAACTC => -3.49832e-05, CAACTG => 0.000135746, 
CAACTT => -0.000205972, CAAGAA => 2.28733e-05, CAAGAC => 9.93438e-05, CAAGAG => 0.0002258, CAAGAT => 4.70756e-05, 
CAAGCA => -6.55766e-05, CAAGCC => 6.84662e-05, CAAGCG => 1.21691e-05, CAAGCT => 1.14972e-05, CAAGGA => 4.77518e-05, 
CAAGGC => 5.84996e-05, CAAGGG => -3.30108e-05, CAAGGT => 1.29828e-05, CAAGTA => -8.75939e-05, CAAGTC => -4.7567e-05, 
CAAGTG => 8.88377e-05, CAAGTT => -0.000108396, CAATAA => -0.000463988, CAATAC => -3.79096e-05, CAATAG => -0.000176621, 
CAATAT => -0.000181897, CAATCA => -0.000170867, CAATCC => -3.9072e-05, CAATCG => 1.49348e-05, CAATCT => -9.00687e-05, 
CAATGA => -0.000265173, CAATGC => -0.000115902, CAATGG => -0.000119033, CAATGT => -0.000177226, CAATTA => -0.000279503, 
CAATTC => -0.000118201, CAATTG => -8.71379e-05, CAATTT => -0.000366628, CACAAA => -0.000182976, CACAAC => 0.000163262, 
CACAAG => 0.000322203, CACAAT => 1.35942e-05, CACACA => -0.00079428, CACACC => 7.05792e-05, CACACG => 6.73363e-05, 
CACACT => -8.67254e-05, CACAGA => -0.000203864, CACAGC => 7.62076e-05, CACAGG => -8.06086e-05, CACAGT => -6.94934e-06, 
CACATA => -0.000144157, CACATC => 0.000225229, CACATG => 0.000172908, CACATT => -0.000108407, CACCAA => -0.000111482, 
CACCAC => 5.012e-05, CACCAG => 0.000270369, CACCAT => -2.08737e-05, CACCCA => -0.000105717, CACCCC => -2.20682e-05, 
CACCCG => 2.38861e-05, CACCCT => -3.15629e-05, CACCGA => -3.99085e-05, CACCGC => 2.60074e-05, CACCGG => 4.89492e-05, 
CACCGT => -6.54128e-06, CACCTA => -5.69467e-05, CACCTC => -1.11053e-05, CACCTG => 0.000318147, CACCTT => -9.26157e-05, 
CACGAA => 5.20765e-05, CACGAC => 0.000170712, CACGAG => 0.000298834, CACGAT => 0.00012037, CACGCA => -5.44531e-05, 
CACGCC => 0.000122219, CACGCG => -4.0589e-06, CACGCT => -2.33682e-05, CACGGA => 7.04096e-05, CACGGC => 0.000137017, 
CACGGG => 2.37517e-05, CACGGT => 1.62516e-07, CACGTA => -2.46221e-05, CACGTC => 5.2626e-05, CACGTG => 0.000111595, 
CACGTT => -3.13662e-05, CACTAA => -0.000213182, CACTAC => 0.000194617, CACTAG => -0.000139977, CACTAT => 2.44887e-05, 
CACTCA => -0.000130044, CACTCC => 8.56188e-06, CACTCG => 9.02021e-05, CACTCT => -5.92061e-05, CACTGA => -0.000349832, 
CACTGC => -1.57746e-05, CACTGG => -6.81884e-05, CACTGT => -0.000116613, CACTTA => -0.0001692, CACTTC => 0.000135402, 
CACTTG => -1.05381e-05, CACTTT => -0.000184422, CAGAAA => 0.000128594, CAGAAC => 0.000469581, CAGAAG => 0.000884534, 
CAGAAT => 0.00021176, CAGACA => -9.71663e-06, CAGACC => 0.000308032, CAGACG => 0.000204944, CAGACT => 5.40718e-05, 
CAGAGA => 1.15224e-06, CAGAGC => 0.000308767, CAGAGG => 9.36534e-05, CAGAGT => 0.000133628, CAGATA => -2.90217e-05, 
CAGATC => 0.000470393, CAGATG => 0.000425743, CAGATT => 0.000119663, CAGCAA => 0.000377754, CAGCAC => 0.000320902, 
CAGCAG => 0.00191664, CAGCAT => 0.000111901, CAGCCA => 0.000165647, CAGCCC => 0.000285272, CAGCCG => 0.000296348, 
CAGCCT => 0.000141688, CAGCGA => 0.000113104, CAGCGC => 0.000283738, CAGCGG => 0.000195954, CAGCGT => 9.57192e-05, 
CAGCTA => 5.81263e-05, CAGCTC => 0.000272059, CAGCTG => 0.0010586, CAGCTT => 7.0026e-05, CAGGAA => 0.000380481, 
CAGGAC => 0.000565439, CAGGAG => 0.00128325, CAGGAT => 0.000449627, CAGGCA => 0.000198553, CAGGCC => 0.000598724, 
CAGGCG => 0.00034704, CAGGCT => 0.000304317, CAGGGA => 0.000130947, CAGGGC => 0.000387347, CAGGGG => 2.16307e-05, 
CAGGGT => 0.00019191, CAGGTA => -9.28496e-06, CAGGTC => 0.000183998, CAGGTG => 0.000681635, CAGGTT => 6.09523e-05, 
CAGTAA => -0.000259723, CAGTAC => 0.000388581, CAGTAG => -0.000159051, CAGTAT => 0.000149979, CAGTCA => -3.35018e-05, 
CAGTCC => 0.000167659, CAGTCG => 0.000169622, CAGTCT => 8.53161e-05, CAGTGA => -0.000325116, CAGTGC => 9.17101e-05, 
CAGTGG => 5.46728e-05, CAGTGT => -4.2241e-05, CAGTTA => -7.59225e-05, CAGTTC => 0.000374816, CAGTTG => 0.000131572, 
CAGTTT => 5.8839e-05, CATAAA => -0.000292011, CATAAC => -8.22036e-05, CATAAG => -1.12098e-05, CATAAT => -0.000203008, 
CATACA => -0.000254139, CATACC => -8.12936e-05, CATACG => -2.46849e-05, CATACT => -8.89716e-05, CATAGA => -0.000129107, 
CATAGC => -6.36102e-05, CATAGG => -7.44691e-05, CATAGT => -0.00010789, CATATA => -0.000308862, CATATC => -6.61912e-05, 
CATATG => -0.000168093, CATATT => -0.00031299, CATCAA => -0.000158567, CATCAC => -1.41644e-05, CATCAG => 0.00011437, 
CATCAT => -0.000153307, CATCCA => -7.17961e-05, CATCCC => 2.17381e-05, CATCCG => 0.000111289, CATCCT => -7.50438e-05, 
CATCGA => -4.01336e-05, CATCGC => 7.29333e-05, CATCGG => -7.81389e-06, CATCGT => -1.53365e-05, CATCTA => -9.96671e-05, 
CATCTC => -7.16946e-05, CATCTG => 9.16175e-05, CATCTT => -0.000186322, CATGAA => -6.72756e-06, CATGAC => 6.98196e-05, 
CATGAG => 0.000208264, CATGAT => 2.05518e-05, CATGCA => -0.000104815, CATGCC => 0.000125619, CATGCG => 1.43491e-05, 
CATGCT => -2.88717e-05, CATGGA => 2.12962e-05, CATGGC => 9.72885e-05, CATGGG => -2.92172e-05, CATGGT => 3.90032e-06, 
CATGTA => -0.000144976, CATGTC => -2.65536e-05, CATGTG => 7.94124e-05, CATGTT => -0.000162831, CATTAA => -0.000370387, 
CATTAC => -4.54482e-05, CATTAG => -0.000194144, CATTAT => -0.000216736, CATTCA => -0.000243669, CATTCC => -8.76321e-05, 
CATTCG => -2.38285e-05, CATTCT => -0.0001687, CATTGA => -0.00028165, CATTGC => -0.000161054, CATTGG => -0.000149354, 
CATTGT => -0.000236787, CATTTA => -0.000431168, CATTTC => -0.000235539, CATTTG => -0.000236739, CATTTT => -0.000698105, 
CCAAAA => -0.000204764, CCAAAC => -2.74765e-05, CCAAAG => 0.00010438, CCAAAT => -0.000123842, CCAACA => -3.63987e-05, 
CCAACC => -1.77105e-05, CCAACG => 6.75834e-06, CCAACT => -3.63099e-05, CCAAGA => -0.000116992, CCAAGC => -3.1471e-05, 
CCAAGG => -8.92718e-05, CCAAGT => -5.93533e-05, CCAATA => -0.00013241, CCAATC => -3.49447e-05, CCAATG => 6.34105e-05, 
CCAATT => -0.000112678, CCACAA => -8.10811e-05, CCACAC => -0.000120124, CCACAG => 0.000215038, CCACAT => -8.6699e-05, 
CCACCA => 0.00028411, CCACCC => 0.00014435, CCACCG => 0.000201388, CCACCT => 0.000192254, CCACGA => 1.59411e-05, 
CCACGC => 4.11673e-07, CCACGG => -8.07645e-06, CCACGT => 1.70112e-05, CCACTA => -5.8623e-05, CCACTC => -2.72578e-05, 
CCACTG => 0.000155507, CCACTT => -0.000117664, CCAGAA => 0.000166387, CCAGAC => 0.000196068, CCAGAG => 0.000460608, 
CCAGAT => 0.000204903, CCAGCA => -9.88195e-06, CCAGCC => 0.000111431, CCAGCG => 5.32332e-05, CCAGCT => 6.64067e-05, 
CCAGGA => 0.000167694, CCAGGC => 9.08464e-05, CCAGGG => -1.7575e-06, CCAGGT => 8.7979e-05, CCAGTA => -2.96688e-05, 
CCAGTC => 3.52003e-05, CCAGTG => 0.000197134, CCAGTT => -2.94897e-05, CCATAA => -0.000220079, CCATAC => 9.04744e-06, 
CCATAG => -0.000147674, CCATAT => -7.40413e-05, CCATCA => -6.28447e-05, CCATCC => 3.73574e-05, CCATCG => 8.36227e-05, 
CCATCT => -4.06602e-05, CCATGA => -0.000213976, CCATGC => -9.79055e-05, CCATGG => -8.90286e-05, CCATGT => -0.0001446, 
CCATTA => -0.000116506, CCATTC => -5.93707e-06, CCATTG => -4.84735e-05, CCATTT => -0.000228015, CCCAAA => 0.000232864, 
CCCAAC => 0.00038695, CCCAAG => 0.000563361, CCCAAT => 0.000228737, CCCACA => 0.000129716, CCCACC => 0.0002757, 
CCCACG => 0.000153212, CCCACT => 8.50519e-05, CCCAGA => -4.96359e-05, CCCAGC => 0.000302093, CCCAGG => -2.72935e-05, 
CCCAGT => 0.000193585, CCCATA => -2.09816e-05, CCCATC => 0.000346012, CCCATG => 0.000309023, CCCATT => 0.000141384, 
CCCCAA => -6.38495e-05, CCCCAC => -5.9185e-05, CCCCAG => 0.000199536, CCCCAT => -3.36712e-05, CCCCCA => 4.73477e-05, 
CCCCCC => -8.50233e-05, CCCCCG => 4.48494e-05, CCCCCT => 3.87239e-05, CCCCGA => 6.63129e-06, CCCCGC => 4.95807e-06, 
CCCCGG => 1.6926e-05, CCCCGT => 5.41099e-06, CCCCTA => -3.43417e-05, CCCCTC => -4.77251e-05, CCCCTG => 0.00016729, 
CCCCTT => -8.50648e-05, CCCGAA => 0.000118847, CCCGAC => 0.000171187, CCCGAG => 0.000333154, CCCGAT => 0.000163845, 
CCCGCA => 1.65729e-05, CCCGCC => 0.00014529, CCCGCG => 2.99344e-05, CCCGCT => 6.53491e-05, CCCGGA => 0.000100857, 
CCCGGC => 0.000157512, CCCGGG => 3.27901e-06, CCCGGT => 6.16688e-05, CCCGTA => 1.20201e-05, CCCGTC => 0.000106985, 
CCCGTG => 0.000188184, CCCGTT => 5.10343e-05, CCCTAA => -0.000150319, CCCTAC => 0.000238755, CCCTAG => -0.000101098, 
CCCTAT => 9.10103e-05, CCCTCA => -4.65153e-06, CCCTCC => 2.01478e-05, CCCTCG => 0.000104367, CCCTCT => -2.50587e-05, 
CCCTGA => -0.000219822, CCCTGC => -2.27005e-06, CCCTGG => -2.90609e-05, CCCTGT => -7.961e-05, CCCTTA => -6.60869e-05, 
CCCTTC => 0.000212268, CCCTTG => 2.73288e-05, CCCTTT => -3.63364e-05, CCGAAA => -5.87941e-05, CCGAAC => 2.85252e-05, 
CCGAAG => 7.43835e-05, CCGAAT => 3.98901e-05, CCGACA => 5.44765e-06, CCGACC => 3.41403e-05, CCGACG => 4.33406e-05, 
CCGACT => 6.2797e-07, CCGAGA => -6.9342e-05, CCGAGC => 4.62432e-06, CCGAGG => -4.32212e-05, CCGAGT => 2.13902e-06, 
CCGATA => -4.61813e-05, CCGATC => 8.1904e-06, CCGATG => 8.04965e-05, CCGATT => -2.09001e-05, CCGCAA => 3.41265e-05, 
CCGCAC => 8.27651e-05, CCGCAG => 0.000317088, CCGCAT => 3.69274e-05, CCGCCA => 0.000197244, CCGCCC => 0.000184245, 
CCGCCG => 0.0002591, CCGCCT => 8.55202e-05, CCGCGA => 1.09461e-05, CCGCGC => 4.78631e-05, CCGCGG => -1.41384e-05, 
CCGCGT => 3.91329e-06, CCGCTA => 9.01392e-06, CCGCTC => 7.27468e-05, CCGCTG => 0.000330843, CCGCTT => -3.67649e-05, 
CCGGAA => 0.000149126, CCGGAC => 0.000180023, CCGGAG => 0.000398216, CCGGAT => 0.000179276, CCGGCA => 5.40716e-05, 
CCGGCC => 0.000146687, CCGGCG => 0.000112454, CCGGCT => 7.24385e-05, CCGGGA => 7.50307e-05, CCGGGC => 0.000130055, 
CCGGGG => -3.37608e-06, CCGGGT => 5.39996e-05, CCGGTA => 5.98725e-06, CCGGTC => 5.16649e-05, CCGGTG => 0.000241864, 
CCGGTT => 3.28317e-05, CCGTAA => -7.68989e-05, CCGTAC => 0.000114672, CCGTAG => -4.70232e-05, CCGTAT => 7.25366e-05, 
CCGTCA => 1.17123e-05, CCGTCC => 8.13339e-05, CCGTCG => 5.54517e-05, CCGTCT => 3.18559e-05, CCGTGA => -9.31469e-05, 
CCGTGC => -2.28281e-05, CCGTGG => -1.63239e-05, CCGTGT => -3.94688e-05, CCGTTA => -3.80398e-05, CCGTTC => 7.88055e-05, 
CCGTTG => 2.15966e-05, CCGTTT => -2.74066e-05, CCTAAA => 1.0894e-05, CCTAAC => 2.01289e-05, CCTAAG => 2.0608e-05, 
CCTAAT => -3.16795e-05, CCTACA => -3.21855e-05, CCTACC => 7.52192e-06, CCTACG => 9.653e-06, CCTACT => -4.43417e-06, 
CCTAGA => -6.32842e-05, CCTAGC => 6.02081e-06, CCTAGG => -4.43889e-05, CCTAGT => -9.50691e-06, CCTATA => -8.10501e-05, 
CCTATC => -1.46788e-06, CCTATG => -3.04532e-06, CCTATT => -1.92564e-05, CCTCAA => 4.62515e-06, CCTCAC => -1.57316e-05, 
CCTCAG => 0.000197412, CCTCAT => -5.79433e-05, CCTCCA => 0.000240398, CCTCCC => 2.67853e-06, CCTCCG => 9.40138e-05, 
CCTCCT => 0.000163977, CCTCGA => 2.39687e-05, CCTCGC => 2.53053e-05, CCTCGG => 3.76678e-06, CCTCGT => -2.36874e-05, 
CCTCTA => -3.84089e-05, CCTCTC => -0.000105846, CCTCTG => 0.000142331, CCTCTT => -0.000112629, CCTGAA => 0.000125581, 
CCTGAC => 0.000140268, CCTGAG => 0.000322395, CCTGAT => 0.000174972, CCTGCA => 4.45838e-05, CCTGCC => 8.73278e-05, 
CCTGCG => 3.36052e-05, CCTGCT => 0.000143568, CCTGGA => 0.000219531, CCTGGC => 0.00014476, CCTGGG => 1.07631e-06, 
CCTGGT => 0.000134296, CCTGTA => -6.40358e-05, CCTGTC => 6.32346e-06, CCTGTG => 0.000213437, CCTGTT => 2.24933e-06, 
CCTTAA => -0.000226472, CCTTAC => 5.21054e-05, CCTTAG => -0.000118815, CCTTAT => -1.25586e-05, CCTTCA => -6.33626e-05, 
CCTTCC => -2.6913e-05, CCTTCG => 1.4054e-05, CCTTCT => -3.13826e-05, CCTTGA => -0.000216935, CCTTGC => -7.68285e-05, 
CCTTGG => -0.000126419, CCTTGT => -0.000137189, CCTTTA => -0.000172372, CCTTTC => -5.47552e-05, CCTTTG => -8.38693e-05, 
CCTTTT => -0.000276866, CGAAAA => -0.000151661, CGAAAC => -1.40306e-05, CGAAAG => 6.0328e-05, CGAAAT => -9.83644e-05, 
CGAACA => -5.27122e-05, CGAACC => 1.19276e-05, CGAACG => -3.49889e-05, CGAACT => -3.20629e-05, CGAAGA => -4.89805e-05, 
CGAAGC => -1.33055e-05, CGAAGG => -8.74865e-06, CGAAGT => -2.73989e-05, CGAATA => -7.85653e-05, CGAATC => -3.16649e-06, 
CGAATG => 3.82102e-05, CGAATT => -6.7604e-05, CGACAA => -3.51078e-05, CGACAC => 4.92668e-05, CGACAG => 0.000116776, 
CGACAT => -1.60364e-05, CGACCA => -6.34795e-06, CGACCC => 7.2746e-05, CGACCG => 2.0572e-05, CGACCT => -4.67717e-06, 
CGACGA => -2.33268e-06, CGACGC => 5.88961e-05, CGACGG => 5.34854e-06, CGACGT => 1.90302e-05, CGACTA => -8.35084e-06, 
CGACTC => 2.48568e-05, CGACTG => 0.000156462, CGACTT => -3.24035e-05, CGAGAA => 3.45339e-05, CGAGAC => 0.000113545, 
CGAGAG => 0.000190383, CGAGAT => 0.000113632, CGAGCA => -2.84175e-07, CGAGCC => 8.67617e-05, CGAGCG => -3.31859e-06, 
CGAGCT => 4.78056e-05, CGAGGA => 2.99755e-05, CGAGGC => 3.59842e-05, CGAGGG => 1.01993e-05, CGAGGT => 2.97222e-05, 
CGAGTA => -1.80429e-05, CGAGTC => 4.76766e-06, CGAGTG => 0.000104809, CGAGTT => -2.59164e-05, CGATAA => -0.000117926, 
CGATAC => 4.87199e-05, CGATAG => -7.0738e-05, CGATAT => -1.74192e-05, CGATCA => -5.31837e-05, CGATCC => 2.83201e-05, 
CGATCG => 8.4269e-06, CGATCT => 1.02461e-05, CGATGA => -0.000123142, CGATGC => -3.36107e-05, CGATGG => -2.97952e-05, 
CGATGT => -5.12247e-05, CGATTA => -8.66595e-05, CGATTC => 3.27753e-05, CGATTG => 7.34217e-06, CGATTT => -8.76625e-05, 
CGCAAA => 0.000155497, CGCAAC => 0.000254329, CGCAAG => 0.000525016, CGCAAT => 0.000118504, CGCACA => 6.26072e-07, 
CGCACC => 0.000262074, CGCACG => 4.93064e-05, CGCACT => 5.67279e-05, CGCAGA => 8.60703e-06, CGCAGC => 0.000192731, 
CGCAGG => 8.50504e-05, CGCAGT => 0.000116246, CGCATA => -5.08414e-06, CGCATC => 0.000340301, CGCATG => 0.000253843, 
CGCATT => 0.000135586, CGCCAA => 3.39239e-05, CGCCAC => 8.30081e-05, CGCCAG => 0.000309961, CGCCAT => 1.34851e-05, 
CGCCCA => -2.89085e-05, CGCCCC => -2.61211e-05, CGCCCG => -1.05849e-05, CGCCCT => -3.39566e-05, CGCCGA => -1.0474e-05, 
CGCCGC => 8.23471e-05, CGCCGG => 3.87672e-05, CGCCGT => 2.85666e-05, CGCCTA => 1.68069e-05, CGCCTC => 7.47275e-05, 
CGCCTG => 0.00035029, CGCCTT => 1.80767e-05, CGCGAA => 9.37664e-05, CGCGAC => 0.000147977, CGCGAG => 0.000323827, 
CGCGAT => 0.000135506, CGCGCA => -2.55727e-05, CGCGCC => 7.88255e-05, CGCGCG => -3.81409e-05, CGCGCT => -1.50232e-05, 
CGCGGA => 2.77288e-05, CGCGGC => 0.000114114, CGCGGG => -1.8392e-05, CGCGGT => 8.84432e-06, CGCGTA => -1.17255e-05, 
CGCGTC => 6.22059e-05, CGCGTG => 7.2625e-05, CGCGTT => -7.84228e-07, CGCTAA => -7.83432e-05, CGCTAC => 0.000272811, 
CGCTAG => -4.53071e-05, CGCTAT => 0.000119508, CGCTCA => 3.76621e-06, CGCTCC => 0.000210232, CGCTCG => 2.97046e-05, 
CGCTCT => 3.09027e-05, CGCTGA => -0.000149825, CGCTGC => 0.000108096, CGCTGG => 0.000125721, CGCTGT => 3.43141e-05, 
CGCTTA => -4.36561e-05, CGCTTC => 0.000314487, CGCTTG => 2.88948e-05, CGCTTT => 5.31836e-05, CGGAAA => 4.64976e-05, 
CGGAAC => 7.02805e-05, CGGAAG => 0.000149802, CGGAAT => 3.19392e-05, CGGACA => -1.27362e-05, CGGACC => 4.05079e-05, 
CGGACG => 4.76415e-06, CGGACT => -1.57071e-06, CGGAGA => -1.08767e-05, CGGAGC => 2.1169e-05, CGGAGG => 2.46322e-06, 
CGGAGT => -7.42178e-06, CGGATA => -1.101e-05, CGGATC => 5.28126e-05, CGGATG => 6.48907e-05, CGGATT => 8.15292e-06, 
CGGCAA => -1.28973e-05, CGGCAC => 6.09478e-05, CGGCAG => 0.000226732, CGGCAT => 9.11735e-08, CGGCCA => -1.98156e-06, 
CGGCCC => 8.66544e-05, CGGCCG => 1.78803e-05, CGGCCT => 2.48391e-05, CGGCGA => 8.95501e-06, CGGCGC => 6.01918e-05, 
CGGCGG => -6.73715e-06, CGGCGT => 1.659e-05, CGGCTA => 7.09751e-07, CGGCTC => 4.37459e-05, CGGCTG => 0.000222078, 
CGGCTT => -3.42048e-05, CGGGAA => 0.000144542, CGGGAC => 0.000201475, CGGGAG => 0.000369815, CGGGAT => 0.000167901, 
CGGGCA => 2.8854e-05, CGGGCC => 0.000143692, CGGGCG => 4.10456e-06, CGGGCT => 4.80392e-05, CGGGGA => -9.95212e-06, 
CGGGGC => 5.44415e-05, CGGGGG => -1.86642e-05, CGGGGT => 1.2676e-05, CGGGTA => -5.19302e-07, CGGGTC => 5.85296e-05, 
CGGGTG => 0.000115507, CGGGTT => 3.14353e-07, CGGTAA => -8.63405e-05, CGGTAC => 6.00979e-05, CGGTAG => -6.38425e-05, 
CGGTAT => -6.41426e-06, CGGTCA => -2.37661e-05, CGGTCC => 1.48058e-05, CGGTCG => -1.92244e-05, CGGTCT => -1.80968e-05, 
CGGTGA => -8.99409e-05, CGGTGC => -4.71636e-05, CGGTGG => -7.47592e-05, CGGTGT => -5.21253e-05, CGGTTA => -5.26757e-05, 
CGGTTC => -2.83439e-07, CGGTTG => -1.53964e-05, CGGTTT => -3.64029e-05, CGTAAA => -5.8711e-05, CGTAAC => 4.73402e-06, 
CGTAAG => 3.45949e-05, CGTAAT => -3.39055e-05, CGTACA => -4.00403e-05, CGTACC => 1.84484e-05, CGTACG => -2.63201e-06, 
CGTACT => -1.80839e-05, CGTAGA => -2.19834e-05, CGTAGC => -4.07415e-06, CGTAGG => -1.04376e-05, CGTAGT => -2.55893e-05, 
CGTATA => -6.31941e-05, CGTATC => 4.0105e-05, CGTATG => -9.54733e-06, CGTATT => -3.8591e-05, CGTCAA => 6.91353e-06, 
CGTCAC => 1.90923e-05, CGTCAG => 0.000166436, CGTCAT => -3.19671e-05, CGTCCA => -1.24315e-05, CGTCCC => 5.04258e-05, 
CGTCCG => 1.89006e-05, CGTCCT => 2.42971e-05, CGTCGA => -3.4636e-06, CGTCGC => 9.62927e-05, CGTCGG => -1.2488e-05, 
CGTCGT => 3.89634e-05, CGTCTA => -1.27932e-05, CGTCTC => 4.34326e-05, CGTCTG => 0.000209524, CGTCTT => -2.14153e-05, 
CGTGAA => 3.82163e-05, CGTGAC => 8.36641e-05, CGTGAG => 0.000161571, CGTGAT => 4.85861e-05, CGTGCA => -1.27277e-05, 
CGTGCC => 0.000129966, CGTGCG => -2.67435e-05, CGTGCT => 3.45522e-05, CGTGGA => 5.32417e-05, CGTGGC => 8.33943e-05, 
CGTGGG => -9.29329e-06, CGTGGT => 4.26383e-05, CGTGTA => -3.19339e-05, CGTGTC => 1.18409e-05, CGTGTG => 9.60907e-05, 
CGTGTT => -7.46649e-05, CGTTAA => -0.000127916, CGTTAC => 2.9956e-05, CGTTAG => -5.9175e-05, CGTTAT => 4.35424e-06, 
CGTTCA => -3.93473e-05, CGTTCC => 4.4876e-05, CGTTCG => 9.59565e-06, CGTTCT => -3.26272e-05, CGTTGA => -0.000104829, 
CGTTGC => -2.29871e-05, CGTTGG => -3.47896e-05, CGTTGT => -6.00943e-05, CGTTTA => -0.000108235, CGTTTC => 1.05076e-05, 
CGTTTG => -2.38927e-05, CGTTTT => -0.000170894, CTAAAA => -0.000257119, CTAAAC => -4.16738e-05, CTAAAG => 0.000130581, 
CTAAAT => -0.00012944, CTAACA => -7.73878e-05, CTAACC => 1.25369e-05, CTAACG => 2.87034e-06, CTAACT => -9.44609e-05, 
CTAAGA => -8.19264e-05, CTAAGC => -3.11969e-05, CTAAGG => -5.766e-05, CTAAGT => -7.02591e-05, CTAATA => -0.000167933, 
CTAATC => -2.77825e-05, CTAATG => 1.34793e-05, CTAATT => -0.000192644, CTACAA => -0.000121106, CTACAC => -9.3167e-05, 
CTACAG => 0.000124876, CTACAT => -0.000152892, CTACCA => -4.53833e-05, CTACCC => 2.16543e-05, CTACCG => 1.23917e-06, 
CTACCT => -6.54116e-05, CTACGA => -8.55918e-06, CTACGC => 2.16986e-05, CTACGG => 1.41773e-05, CTACGT => -2.60575e-05, 
CTACTA => -7.98336e-05, CTACTC => -9.93149e-06, CTACTG => 3.18818e-05, CTACTT => -0.000125024, CTAGAA => 2.15429e-05, 
CTAGAC => 6.13797e-05, CTAGAG => 0.000181, CTAGAT => 4.91073e-05, CTAGCA => -5.02152e-05, CTAGCC => 4.40389e-05, 
CTAGCG => 1.37e-05, CTAGCT => -2.17602e-05, CTAGGA => -3.54979e-05, CTAGGC => 2.26481e-05, CTAGGG => -3.00292e-05, 
CTAGGT => -3.71271e-05, CTAGTA => -5.99707e-05, CTAGTC => -9.78324e-06, CTAGTG => 4.25723e-05, CTAGTT => -0.000108024, 
CTATAA => -0.000236155, CTATAC => -5.60781e-05, CTATAG => -0.000154781, CTATAT => -0.000183243, CTATCA => -8.9049e-05, 
CTATCC => 8.3169e-06, CTATCG => 1.15526e-05, CTATCT => -9.99187e-05, CTATGA => -0.000168062, CTATGC => -8.73865e-05, 
CTATGG => -7.1792e-05, CTATGT => -0.000142697, CTATTA => -0.000181577, CTATTC => -5.66933e-05, CTATTG => -8.68055e-05, 
CTATTT => -0.000293543, CTCAAA => 9.67589e-05, CTCAAC => 0.000443188, CTCAAG => 0.000595457, CTCAAT => 0.000193551, 
CTCACA => 1.2878e-05, CTCACC => 0.000332322, CTCACG => 8.83842e-05, CTCACT => 4.83682e-05, CTCAGA => -0.000129645, 
CTCAGC => 0.000234017, CTCAGG => -2.46604e-05, CTCAGT => 6.02479e-05, CTCATA => -5.77171e-05, CTCATC => 0.000451722, 
CTCATG => 0.00024834, CTCATT => 3.85041e-05, CTCCAA => -6.85492e-05, CTCCAC => -5.9673e-05, CTCCAG => 0.000158877, 
CTCCAT => -0.000114185, CTCCCA => -9.50217e-05, CTCCCC => -9.6128e-05, CTCCCG => -1.0603e-05, CTCCCT => -0.000117927, 
CTCCGA => 1.06083e-05, CTCCGC => 4.48369e-05, CTCCGG => 2.20214e-05, CTCCGT => -2.5052e-05, CTCCTA => -4.04067e-05, 
CTCCTC => -4.38652e-05, CTCCTG => 0.000299477, CTCCTT => -0.000101536, CTCGAA => -1.69711e-05, CTCGAC => 0.000100761, 
CTCGAG => 0.000149731, CTCGAT => 7.79609e-05, CTCGCA => -4.04684e-05, CTCGCC => 7.6069e-05, CTCGCG => -1.27274e-05, 
CTCGCT => -3.12928e-05, CTCGGA => 2.1328e-05, CTCGGC => 9.14024e-05, CTCGGG => 4.54438e-06, CTCGGT => 1.0018e-05, 
CTCGTA => -2.89942e-05, CTCGTC => 1.60048e-06, CTCGTG => 6.33581e-05, CTCGTT => -3.4738e-05, CTCTAA => -0.000215485, 
CTCTAC => 0.000303813, CTCTAG => -0.000156178, CTCTAT => 7.7911e-05, CTCTCA => -8.67847e-05, CTCTCC => 8.83546e-05, 
CTCTCG => 3.97203e-05, CTCTCT => -0.000253777, CTCTGA => -0.000328754, CTCTGC => -5.0895e-05, CTCTGG => -0.000120536, 
CTCTGT => -0.000194064, CTCTTA => -0.000133208, CTCTTC => 0.000350287, CTCTTG => -6.27986e-05, CTCTTT => -0.000101639, 
CTGAAA => 0.000168605, CTGAAC => 0.000413762, CTGAAG => 0.00109897, CTGAAT => 0.000125503, CTGACA => 7.19634e-05, 
CTGACC => 0.000439968, CTGACG => 0.000199239, CTGACT => 4.72541e-05, CTGAGA => -4.49988e-05, CTGAGC => 0.00031318, 
CTGAGG => 0.00011908, CTGAGT => 0.000101541, CTGATA => -5.24816e-05, CTGATC => 0.000423585, CTGATG => 0.000417184, 
CTGATT => 8.44354e-05, CTGCAA => 0.000117652, CTGCAC => 0.000344087, CTGCAG => 0.00129408, CTGCAT => 2.9528e-05, 
CTGCCA => 0.000120976, CTGCCC => 0.000563287, CTGCCG => 0.00024179, CTGCCT => 7.04616e-05, CTGCGA => 0.000156241, 
CTGCGC => 0.000529788, CTGCGG => 0.000255348, CTGCGT => 0.000155881, CTGCTA => 7.20921e-05, CTGCTC => 0.000520638, 
CTGCTG => 0.0013831, CTGCTT => -2.04077e-06, CTGGAA => 0.000491731, CTGGAC => 0.00106884, CTGGAG => 0.00189904, 
CTGGAT => 0.000650189, CTGGCA => 0.00024049, CTGGCC => 0.00109219, CTGGCG => 0.00029366, CTGGCT => 0.000417077, 
CTGGGA => 0.000240644, CTGGGC => 0.000712056, CTGGGG => 6.82058e-05, CTGGGT => 0.000191104, CTGGTA => 4.49309e-05, 
CTGGTC => 0.000408922, CTGGTG => 0.000983086, CTGGTT => 0.000114719, CTGTAA => -0.00035471, CTGTAC => 0.0002895, 
CTGTAG => -0.000212863, CTGTAT => 2.9122e-05, CTGTCA => -2.75288e-05, CTGTCC => 0.000335757, CTGTCG => 0.000136714, 
CTGTCT => 2.08037e-05, CTGTGA => -0.000360156, CTGTGC => 0.000106608, CTGTGG => -4.35873e-06, CTGTGT => -0.000153318, 
CTGTTA => -0.000113624, CTGTTC => 0.000196639, CTGTTG => 6.41933e-05, CTGTTT => -9.35927e-05, CTTAAA => -0.00026807, 
CTTAAC => -2.32245e-05, CTTAAG => 2.82999e-05, CTTAAT => -0.000126766, CTTACA => -0.000109848, CTTACC => 1.26638e-05, 
CTTACG => -1.55396e-05, CTTACT => -7.85973e-05, CTTAGA => -0.00013103, CTTAGC => -2.02907e-05, CTTAGG => -6.92318e-05, 
CTTAGT => -6.78053e-05, CTTATA => -0.000124362, CTTATC => 1.66947e-06, CTTATG => -2.77762e-05, CTTATT => -0.00016741, 
CTTCAA => -0.000112487, CTTCAC => -6.37767e-05, CTTCAG => 0.000178393, CTTCAT => -0.000155668, CTTCCA => -4.02872e-05, 
CTTCCC => -5.62822e-05, CTTCCG => -6.95278e-06, CTTCCT => -7.68024e-05, CTTCGA => 1.35583e-05, CTTCGC => 4.65713e-05, 
CTTCGG => 2.81346e-05, CTTCGT => -1.28713e-05, CTTCTA => -8.8826e-05, CTTCTC => -7.34821e-05, CTTCTG => 9.75977e-05, 
CTTCTT => -0.000164396, CTTGAA => -7.10681e-05, CTTGAC => 4.61288e-05, CTTGAG => 0.000102583, CTTGAT => 4.47763e-06, 
CTTGCA => -5.552e-05, CTTGCC => 2.10658e-05, CTTGCG => -1.65071e-05, CTTGCT => -6.95107e-05, CTTGGA => -1.0634e-06, 
CTTGGC => 3.71911e-05, CTTGGG => -9.20742e-05, CTTGGT => -3.69884e-05, CTTGTA => -0.000151256, CTTGTC => -4.68228e-05, 
CTTGTG => -1.73457e-05, CTTGTT => -0.000199906, CTTTAA => -0.000520416, CTTTAC => -5.46884e-05, CTTTAG => -0.000231591, 
CTTTAT => -0.000238653, CTTTCA => -0.000206222, CTTTCC => -9.10689e-05, CTTTCG => -3.42651e-05, CTTTCT => -0.000243729, 
CTTTGA => -0.000359641, CTTTGC => -0.000195437, CTTTGG => -0.000197517, CTTTGT => -0.000303521, CTTTTA => -0.000405872, 
CTTTTC => -0.000208987, CTTTTG => -0.000196338, CTTTTT => -0.000608554, GAAAAA => -0.000229384, GAAAAC => 0.000142795, 
GAAAAG => 0.000444498, GAAAAT => -9.24188e-05, GAAACA => -0.000129625, GAAACC => 9.64168e-05, GAAACG => 5.67823e-05, 
GAAACT => -2.73604e-05, GAAAGA => -0.000142904, GAAAGC => 2.52471e-05, GAAAGG => -4.11937e-05, GAAAGT => -1.02506e-05, 
GAAATA => -0.00023412, GAAATC => 0.000167261, GAAATG => 0.000132982, GAAATT => -0.000100873, GAACAA => 7.73723e-06, 
GAACAC => 1.05697e-05, GAACAG => 0.000377422, GAACAT => -9.17839e-05, GAACCA => 7.50654e-05, GAACCC => 0.000133776, 
GAACCG => 7.43509e-05, GAACCT => 7.0643e-05, GAACGA => 3.6238e-05, GAACGC => 0.000133514, GAACGG => 8.10633e-05, 
GAACGT => 3.16696e-05, GAACTA => 7.50422e-06, GAACTC => 9.29795e-05, GAACTG => 0.000503253, GAACTT => 3.75393e-05, 
GAAGAA => 0.000583366, GAAGAC => 0.000527982, GAAGAG => 0.0010178, GAAGAT => 0.000554538, GAAGCA => 0.000112986, 
GAAGCC => 0.000382759, GAAGCG => 9.5662e-05, GAAGCT => 0.000283181, GAAGGA => 0.00016631, GAAGGC => 0.000225011, 
GAAGGG => 5.40367e-05, GAAGGT => 0.000114478, GAAGTA => 3.71922e-05, GAAGTC => 0.000108664, GAAGTG => 0.000396868, 
GAAGTT => 6.99883e-05, GAATAA => -0.00035, GAATAC => 0.000148287, GAATAG => -0.000149038, GAATAT => -6.09468e-05, 
GAATCA => -4.4464e-05, GAATCC => 8.90138e-05, GAATCG => 6.50913e-05, GAATCT => 2.97169e-05, GAATGA => -0.000328451, 
GAATGC => 5.94476e-06, GAATGG => 1.32728e-05, GAATGT => -6.82192e-05, GAATTA => -0.000156113, GAATTC => 8.87291e-05, 
GAATTG => 5.61582e-06, GAATTT => -0.000180116, GACAAA => 0.000303252, GACAAC => 0.000559032, GACAAG => 0.000898461, 
GACAAT => 0.000282599, GACACA => 0.000143355, GACACC => 0.000443269, GACACG => 0.000214395, GACACT => 0.000204644, 
GACAGA => 1.06072e-05, GACAGC => 0.000536881, GACAGG => 9.18979e-05, GACAGT => 0.000268371, GACATA => 9.45696e-05, 
GACATC => 0.000700351, GACATG => 0.000607904, GACATT => 0.00030801, GACCAA => 5.36429e-05, GACCAC => 0.000163914, 
GACCAG => 0.000515651, GACCAT => 4.59019e-05, GACCCA => 0.000128625, GACCCC => 0.000249863, GACCCG => 0.000109397, 
GACCCT => 0.000177229, GACCGA => 8.93862e-05, GACCGC => 0.000173496, GACCGG => 0.000130769, GACCGT => 9.19401e-05, 
GACCTA => 5.54096e-05, GACCTC => 0.000296027, GACCTG => 0.000781253, GACCTT => 9.02991e-05, GACGAA => 0.000223792, 
GACGAC => 0.000516432, GACGAG => 0.000847494, GACGAT => 0.000386669, GACGCA => 5.38048e-05, GACGCC => 0.000262167, 
GACGCG => 7.75571e-05, GACGCT => 0.000105898, GACGGA => 0.000123591, GACGGC => 0.000250913, GACGGG => 8.89843e-05, 
GACGGT => 9.70682e-05, GACGTA => 2.39932e-05, GACGTC => 0.000155031, GACGTG => 0.00033361, GACGTT => 5.03591e-05, 
GACTAA => -0.000150841, GACTAC => 0.000505614, GACTAG => -0.000104266, GACTAT => 0.000256283, GACTCA => 5.62626e-05, 
GACTCC => 0.000343093, GACTCG => 0.000200631, GACTCT => 0.000200639, GACTGA => -0.000265021, GACTGC => 0.000208642, 
GACTGG => 0.000208927, GACTGT => 6.01426e-05, GACTTA => -3.36125e-05, GACTTC => 0.000546199, GACTTG => 0.00020924, 
GACTTT => 0.00025064, GAGAAA => 0.000581848, GAGAAC => 0.000847066, GAGAAG => 0.00177189, GAGAAT => 0.000432935, 
GAGACA => 0.000119783, GAGACC => 0.000519255, GAGACG => 0.000307951, GAGACT => 0.000200104, GAGAGA => -5.34869e-05, 
GAGAGC => 0.000485875, GAGAGG => 0.000155467, GAGAGT => 0.000243718, GAGATA => 5.58269e-05, GAGATC => 0.000866082, 
GAGATG => 0.000645974, GAGATT => 0.000352043, GAGCAA => 0.000254355, GAGCAC => 0.000394662, GAGCAG => 0.00125597, 
GAGCAT => 0.000135249, GAGCCA => 0.000213151, GAGCCC => 0.000398617, GAGCCG => 0.000246118, GAGCCT => 0.000238523, 
GAGCGA => 0.000203522, GAGCGC => 0.000462072, GAGCGG => 0.000268116, GAGCGT => 0.000163621, GAGCTA => 0.000129972, 
GAGCTC => 0.000437274, GAGCTG => 0.00164119, GAGCTT => 0.000229567, GAGGAA => 0.000946121, GAGGAC => 0.00114683, 
GAGGAG => 0.00261521, GAGGAT => 0.00103885, GAGGCA => 0.000279535, GAGGCC => 0.00079218, GAGGCG => 0.000383783, 
GAGGCT => 0.000468906, GAGGGA => 0.000195897, GAGGGC => 0.000593674, GAGGGG => 7.0164e-05, GAGGGT => 0.000279942, 
GAGGTA => 7.32227e-05, GAGGTC => 0.000340814, GAGGTG => 0.000968542, GAGGTT => 0.000222965, GAGTAA => -0.000182616, 
GAGTAC => 0.00051369, GAGTAG => -0.000127203, GAGTAT => 0.000226404, GAGTCA => 8.29465e-05, GAGTCC => 0.000318649, 
GAGTCG => 0.000197151, GAGTCT => 0.000198944, GAGTGA => -0.000226678, GAGTGC => 0.000255607, GAGTGG => 0.000170794, 
GAGTGT => 0.000100631, GAGTTA => -2.12677e-05, GAGTTC => 0.000580567, GAGTTG => 0.000245878, GAGTTT => 0.00033816, 
GATAAA => 9.38568e-05, GATAAC => 0.000156153, GATAAG => 0.000309956, GATAAT => 7.01005e-05, GATACA => -3.93589e-05, 
GATACC => 9.97943e-05, GATACG => 7.69961e-05, GATACT => 3.14424e-06, GATAGA => -5.24028e-05, GATAGC => 8.70291e-05, 
GATAGG => -9.29687e-06, GATAGT => 1.38154e-05, GATATA => -5.84907e-05, GATATC => 0.000171673, GATATG => 0.000151407, 
GATATT => 2.81092e-05, GATCAA => 5.36163e-05, GATCAC => 0.000131431, GATCAG => 0.000435907, GATCAT => 4.66654e-05, 
GATCCA => 0.000145535, GATCCC => 0.000295584, GATCCG => 0.000144497, GATCCT => 0.00015746, GATCGA => 7.36209e-05, 
GATCGC => 0.00024975, GATCGG => 8.81948e-05, GATCGT => 9.67497e-05, GATCTA => 5.06683e-05, GATCTC => 0.000195253, 
GATCTG => 0.000672519, GATCTT => 0.000107661, GATGAA => 0.000742498, GATGAC => 0.000779827, GATGAG => 0.00121602, 
GATGAT => 0.000712567, GATGCA => 0.000240367, GATGCC => 0.000655019, GATGCG => 0.000186939, GATGCT => 0.000345839, 
GATGGA => 0.000387328, GATGGC => 0.000559604, GATGGG => 0.000197325, GATGGT => 0.000253972, GATGTA => 3.28468e-05, 
GATGTC => 0.000311136, GATGTG => 0.000783413, GATGTT => 0.000165881, GATTAA => -0.000275218, GATTAC => 0.000166734, 
GATTAG => -0.000145629, GATTAT => -7.95705e-06, GATTCA => 2.50587e-06, GATTCC => 8.22633e-05, GATTCG => 0.000100699, 
GATTCT => 1.72015e-05, GATTGA => -0.000225681, GATTGC => -1.41078e-05, GATTGG => -1.15936e-05, GATTGT => -0.000114068, 
GATTTA => -0.000170522, GATTTC => 0.000107763, GATTTG => 7.76486e-05, GATTTT => -0.000245757, GCAAAA => -0.00019472, 
GCAAAC => -7.95089e-05, GCAAAG => 5.25044e-05, GCAAAT => -0.000158223, GCAACA => -4.62505e-05, GCAACC => 1.34423e-05, 
GCAACG => -7.90554e-06, GCAACT => -3.08661e-05, GCAAGA => -0.00010866, GCAAGC => -6.43625e-05, GCAAGG => -6.32464e-05, 
GCAAGT => -7.19939e-05, GCAATA => -0.000152706, GCAATC => -3.06912e-05, GCAATG => 3.05989e-05, GCAATT => -0.0001517, 
GCACAA => -7.30117e-05, GCACAC => -0.000122136, GCACAG => 0.00012064, GCACAT => -0.000105176, GCACCA => 7.52217e-05, 
GCACCC => 0.000104984, GCACCG => 6.75143e-05, GCACCT => 6.25485e-05, GCACGA => 1.23513e-05, GCACGC => 1.20948e-05, 
GCACGG => 1.39882e-05, GCACGT => -1.8177e-05, GCACTA => -3.12293e-05, GCACTC => -3.58523e-06, GCACTG => 0.000209255, 
GCACTT => -7.71106e-05, GCAGAA => 0.000138189, GCAGAC => 0.000205871, GCAGAG => 0.000489559, GCAGAT => 0.000203554, 
GCAGCA => 0.000145918, GCAGCC => 0.000369973, GCAGCG => 7.94116e-05, GCAGCT => 0.000232505, GCAGGA => 0.000122179, 
GCAGGC => 7.36104e-05, GCAGGG => -1.59251e-05, GCAGGT => 1.68132e-05, GCAGTA => -3.16505e-05, GCAGTC => 5.41675e-06, 
GCAGTG => 0.00015295, GCAGTT => -4.57692e-05, GCATAA => -0.00023402, GCATAC => -3.42352e-05, GCATAG => -0.000133445, 
GCATAT => -0.000123281, GCATCA => -2.85486e-05, GCATCC => 9.63383e-05, GCATCG => 3.30674e-05, GCATCT => 1.37511e-05, 
GCATGA => -0.000196179, GCATGC => -0.000111488, GCATGG => -8.6207e-05, GCATGT => -0.000150658, GCATTA => -0.000106415, 
GCATTC => -3.38975e-05, GCATTG => -4.58996e-05, GCATTT => -0.000263025, GCCAAA => 0.000410102, GCCAAC => 0.000601397, 
GCCAAG => 0.00110288, GCCAAT => 0.000417984, GCCACA => 0.00026105, GCCACC => 0.000650306, GCCACG => 0.000234292, 
GCCACT => 0.000226488, GCCAGA => 4.76651e-05, GCCAGC => 0.000493572, GCCAGG => 6.55028e-05, GCCAGT => 0.000319134, 
GCCATA => 7.43136e-05, GCCATC => 0.000731422, GCCATG => 0.000647763, GCCATT => 0.000351241, GCCCAA => 8.55251e-05, 
GCCCAC => 0.000151658, GCCCAG => 0.000638221, GCCCAT => 6.58755e-05, GCCCCA => 4.62548e-05, GCCCCC => 8.31956e-05, 
GCCCCG => 5.61042e-05, GCCCCT => 7.62368e-05, GCCCGA => 7.40155e-05, GCCCGC => 0.000174368, GCCCGG => 9.76034e-05, 
GCCCGT => 7.74113e-05, GCCCTA => 3.52487e-05, GCCCTC => 0.000207793, GCCCTG => 0.000645725, GCCCTT => 6.00276e-05, 
GCCGAA => 0.000193559, GCCGAC => 0.000232416, GCCGAG => 0.000565085, GCCGAT => 0.000265231, GCCGCA => 9.97642e-05, 
GCCGCC => 0.000506725, GCCGCG => 5.78372e-05, GCCGCT => 0.000182052, GCCGGA => 0.000191364, GCCGGC => 0.000274569, 
GCCGGG => 2.32687e-05, GCCGGT => 0.000127641, GCCGTA => 4.20939e-05, GCCGTC => 0.000201981, GCCGTG => 0.000315612, 
GCCGTT => 6.76083e-05, GCCTAA => -0.000121403, GCCTAC => 0.000391486, GCCTAG => -0.000108611, GCCTAT => 0.000200608, 
GCCTCA => 4.89578e-05, GCCTCC => 0.000349302, GCCTCG => 0.000138929, GCCTCT => 8.18831e-05, GCCTGA => -0.0002032, 
GCCTGC => 0.000147509, GCCTGG => 8.22954e-05, GCCTGT => 2.32253e-05, GCCTTA => -2.73618e-05, GCCTTC => 0.000531865, 
GCCTTG => 0.000172492, GCCTTT => 0.000147695, GCGAAA => -6.32395e-05, GCGAAC => 2.90523e-05, GCGAAG => 8.48293e-05, 
GCGAAT => -2.80083e-07, GCGACA => -9.35874e-06, GCGACC => 3.19935e-05, GCGACG => 2.64089e-05, GCGACT => -2.04189e-05, 
GCGAGA => -6.39475e-05, GCGAGC => 3.15996e-06, GCGAGG => -4.20439e-05, GCGAGT => -1.85839e-05, GCGATA => -3.36306e-05, 
GCGATC => 1.99982e-05, GCGATG => 7.82665e-05, GCGATT => -4.34018e-06, GCGCAA => -8.96939e-06, GCGCAC => 3.24502e-05, 
GCGCAG => 0.000174533, GCGCAT => -2.04835e-05, GCGCCA => 2.54085e-05, GCGCCC => 0.000127997, GCGCCG => 5.24662e-05, 
GCGCCT => 1.85304e-05, GCGCGA => -4.4772e-06, GCGCGC => -6.16604e-06, GCGCGG => -2.72514e-05, GCGCGT => -2.19923e-05, 
GCGCTA => 2.75936e-06, GCGCTC => 6.19001e-05, GCGCTG => 0.000314542, GCGCTT => -3.00483e-05, GCGGAA => 0.00010578, 
GCGGAC => 0.000197143, GCGGAG => 0.00041569, GCGGAT => 0.000184888, GCGGCA => 0.000111785, GCGGCC => 0.000353628, 
GCGGCG => 0.000171258, GCGGCT => 0.000151208, GCGGGA => 4.59326e-05, GCGGGC => 0.00015201, GCGGGG => -6.54534e-06, 
GCGGGT => 6.15657e-05, GCGGTA => -7.5947e-06, GCGGTC => 7.44804e-05, GCGGTG => 0.000271708, GCGGTT => 2.51495e-05, 
GCGTAA => -9.19563e-05, GCGTAC => 7.65966e-05, GCGTAG => -5.66894e-05, GCGTAT => 1.03984e-05, GCGTCA => -2.64254e-05, 
GCGTCC => 7.52413e-05, GCGTCG => -2.14136e-05, GCGTCT => 2.33e-05, GCGTGA => -9.76735e-05, GCGTGC => -4.43574e-05, 
GCGTGG => -5.05043e-05, GCGTGT => -9.45169e-05, GCGTTA => -3.17675e-05, GCGTTC => 3.19381e-05, GCGTTG => 1.57719e-05, 
GCGTTT => -3.75375e-05, GCTAAA => 2.31438e-05, GCTAAC => 5.05294e-05, GCTAAG => 0.000120201, GCTAAT => -2.72217e-05, 
GCTACA => -2.99054e-05, GCTACC => 4.03777e-05, GCTACG => 2.21294e-05, GCTACT => -4.58987e-06, GCTAGA => -6.38013e-05, 
GCTAGC => 2.37878e-05, GCTAGG => -4.94404e-05, GCTAGT => -3.44239e-06, GCTATA => -5.19427e-05, GCTATC => 5.37191e-05, 
GCTATG => 8.16532e-05, GCTATT => -1.85969e-06, GCTCAA => 4.52767e-05, GCTCAC => 4.79408e-05, GCTCAG => 0.000310511, 
GCTCAT => -5.43135e-06, GCTCCA => 0.000165528, GCTCCC => 0.000119917, GCTCCG => 7.4294e-05, GCTCCT => 0.000145631, 
GCTCGA => 5.8083e-05, GCTCGC => 5.97504e-05, GCTCGG => 3.18335e-05, GCTCGT => 4.35806e-05, GCTCTA => 1.10698e-05, 
GCTCTC => 3.70661e-05, GCTCTG => 0.000484888, GCTCTT => 7.39668e-06, GCTGAA => 0.000179246, GCTGAC => 0.000248028, 
GCTGAG => 0.000470509, GCTGAT => 0.000228597, GCTGCA => 0.000172449, GCTGCC => 0.000494405, GCTGCG => 8.30694e-05, 
GCTGCT => 0.000304893, GCTGGA => 0.000240042, GCTGGC => 0.000227972, GCTGGG => 1.13923e-05, GCTGGT => 0.000131883, 
GCTGTA => -3.003e-05, GCTGTC => 9.65104e-05, GCTGTG => 0.000473822, GCTGTT => -3.6162e-05, GCTTAA => -0.000221766, 
GCTTAC => 4.08989e-05, GCTTAG => -0.000141674, GCTTAT => -3.84326e-05, GCTTCA => -5.97243e-05, GCTTCC => 4.31613e-05, 
GCTTCG => 2.42002e-05, GCTTCT => -5.37947e-05, GCTTGA => -0.000195798, GCTTGC => -7.63841e-05, GCTTGG => -0.000114592, 
GCTTGT => -0.000140308, GCTTTA => -0.000133748, GCTTTC => 1.07846e-05, GCTTTG => 3.83756e-05, GCTTTT => -0.000188132, 
GGAAAA => -0.000143173, GGAAAC => 0.000122062, GGAAAG => 0.000233238, GGAAAT => -3.09526e-05, GGAACA => -1.21268e-05, 
GGAACC => 0.00012506, GGAACG => 4.71674e-05, GGAACT => 1.89242e-05, GGAAGA => -0.000171614, GGAAGC => 3.34905e-05, 
GGAAGG => -0.000138344, GGAAGT => -4.9806e-06, GGAATA => -6.83552e-05, GGAATC => 0.000106311, GGAATG => 0.000112936, 
GGAATT => -1.44632e-05, GGACAA => -2.78528e-05, GGACAC => 5.55675e-05, GGACAG => 0.000233307, GGACAT => -3.02103e-05, 
GGACCA => 0.000172116, GGACCC => 0.000247713, GGACCG => 5.93351e-05, GGACCT => 0.000142228, GGACGA => 3.98059e-05, 
GGACGC => 0.000103285, GGACGG => 3.10314e-05, GGACGT => 2.6715e-05, GGACTA => -3.318e-05, GGACTC => 7.22254e-05, 
GGACTG => 0.000215518, GGACTT => 3.19598e-06, GGAGAA => 0.000263603, GGAGAC => 0.000358921, GGAGAG => 0.000488793, 
GGAGAT => 0.00030515, GGAGCA => 0.000145244, GGAGCC => 0.000329839, GGAGCG => 7.31409e-05, GGAGCT => 0.000199295, 
GGAGGA => 0.000283855, GGAGGC => 0.000246936, GGAGGG => -4.65873e-05, GGAGGT => 0.000151029, GGAGTA => 3.19953e-06, 
GGAGTC => 9.64226e-05, GGAGTG => 0.00023458, GGAGTT => 5.49454e-05, GGATAA => -0.000177135, GGATAC => 0.000136084, 
GGATAG => -0.000114422, GGATAT => -1.61001e-05, GGATCA => 3.45791e-06, GGATCC => 0.000156925, GGATCG => 0.000104112, 
GGATCT => 7.35342e-05, GGATGA => -0.000244345, GGATGC => -4.19567e-05, GGATGG => -6.9253e-05, GGATGT => -0.00012374, 
GGATTA => -9.90596e-05, GGATTC => 0.000153804, GGATTG => 4.48278e-05, GGATTT => -2.67441e-05, GGCAAA => 0.000245712, 
GGCAAC => 0.000425183, GGCAAG => 0.000704037, GGCAAT => 0.000238093, GGCACA => 0.000134873, GGCACC => 0.000425607, 
GGCACG => 0.000105573, GGCACT => 0.000175392, GGCAGA => -7.39188e-05, GGCAGC => 0.000501277, GGCAGG => -1.88216e-05, 
GGCAGT => 0.000240067, GGCATA => 1.39175e-05, GGCATC => 0.000552022, GGCATG => 0.000415996, GGCATT => 0.000217051, 
GGCCAA => 7.42355e-06, GGCCAC => 0.000133709, GGCCAG => 0.000355228, GGCCAT => 9.94756e-06, GGCCCA => 4.65877e-05, 
GGCCCC => 0.000105652, GGCCCG => 2.28946e-05, GGCCCT => 8.58351e-05, GGCCGA => 6.0287e-05, GGCCGC => 0.000154715, 
GGCCGG => 7.34308e-05, GGCCGT => 7.33636e-05, GGCCTA => 1.41087e-05, GGCCTC => 0.000107411, GGCCTG => 0.000416535, 
GGCCTT => 1.27238e-05, GGCGAA => 0.000147189, GGCGAC => 0.00025345, GGCGAG => 0.000488633, GGCGAT => 0.000285069, 
GGCGCA => 5.06078e-05, GGCGCC => 0.000227155, GGCGCG => 8.94504e-06, GGCGCT => 9.11531e-05, GGCGGA => 0.000235482, 
GGCGGC => 0.000412868, GGCGGG => 2.93474e-05, GGCGGT => 0.000229327, GGCGTA => 2.8877e-05, GGCGTC => 0.000169619, 
GGCGTG => 0.000224718, GGCGTT => 7.54095e-05, GGCTAA => -0.000145374, GGCTAC => 0.000431475, GGCTAG => -0.000119731, 
GGCTAT => 0.000205415, GGCTCA => 4.71781e-05, GGCTCC => 0.000346728, GGCTCG => 7.53375e-05, GGCTCT => 0.000128294, 
GGCTGA => -0.000253421, GGCTGC => 0.000113848, GGCTGG => 0.000101859, GGCTGT => 1.60696e-05, GGCTTA => -3.25522e-05, 
GGCTTC => 0.000449261, GGCTTG => 8.25029e-05, GGCTTT => 0.000109794, GGGAAA => -1.89623e-05, GGGAAC => 7.84636e-05, 
GGGAAG => 0.000199249, GGGAAT => -7.04473e-06, GGGACA => -3.01351e-05, GGGACC => 1.54761e-05, GGGACG => 1.50761e-05, 
GGGACT => -2.75426e-05, GGGAGA => -0.000123497, GGGAGC => -1.04162e-05, GGGAGG => -0.000176159, GGGAGT => -6.14338e-05, 
GGGATA => -5.29934e-05, GGGATC => 5.74772e-05, GGGATG => 4.15466e-05, GGGATT => -1.47941e-05, GGGCAA => -6.05091e-05, 
GGGCAC => -1.67632e-05, GGGCAG => 5.52412e-05, GGGCAT => -6.51674e-05, GGGCCA => 7.43104e-06, GGGCCC => 7.80137e-05, 
GGGCCG => 1.42262e-05, GGGCCT => 2.28174e-06, GGGCGA => -1.41085e-05, GGGCGC => -3.90048e-06, GGGCGG => -5.53788e-05, 
GGGCGT => -2.29178e-05, GGGCTA => -4.23602e-05, GGGCTC => 2.05236e-06, GGGCTG => 9.81063e-05, GGGCTT => -6.73187e-05, 
GGGGAA => 2.05822e-05, GGGGAC => 0.000152483, GGGGAG => 0.000200851, GGGGAT => 0.000134688, GGGGCA => 1.71441e-05, 
GGGGCC => 0.000129633, GGGGCG => 1.27731e-05, GGGGCT => 4.53076e-06, GGGGGA => -6.31744e-05, GGGGGC => 0.000109418, 
GGGGGG => -0.000138924, GGGGGT => 1.98871e-05, GGGGTA => -3.78888e-05, GGGGTC => 4.37226e-05, GGGGTG => 3.67938e-05, 
GGGGTT => -4.36925e-05, GGGTAA => -0.000126349, GGGTAC => 2.8444e-05, GGGTAG => -9.87761e-05, GGGTAT => -3.03341e-05, 
GGGTCA => -5.72009e-05, GGGTCC => 8.44204e-06, GGGTCG => -1.54638e-05, GGGTCT => -3.63256e-05, GGGTGA => -0.000150333, 
GGGTGC => -8.5968e-05, GGGTGG => -0.000172796, GGGTGT => -0.0001265, GGGTTA => -7.71318e-05, GGGTTC => -5.13445e-05, 
GGGTTG => -6.43275e-05, GGGTTT => -0.000129359, GGTAAA => -3.20913e-05, GGTAAC => 3.33212e-05, GGTAAG => -1.05366e-05, 
GGTAAT => -2.93074e-05, GGTACA => -2.40091e-05, GGTACC => 3.43005e-05, GGTACG => 5.6381e-06, GGTACT => -6.43816e-06, 
GGTAGA => -7.90377e-05, GGTAGC => 3.48038e-05, GGTAGG => -9.5398e-05, GGTAGT => -8.80265e-06, GGTATA => -8.78079e-05, 
GGTATC => 4.00117e-05, GGTATG => 6.63882e-06, GGTATT => -1.91166e-05, GGTCAA => 1.6761e-05, GGTCAC => 6.18018e-05, 
GGTCAG => 0.000179799, GGTCAT => 3.44342e-05, GGTCCA => 9.61231e-05, GGTCCC => 0.000104068, GGTCCG => 6.0831e-05, 
GGTCCT => 0.000130503, GGTCGA => 4.63453e-05, GGTCGC => 0.000139889, GGTCGG => 7.48712e-06, GGTCGT => 8.33953e-05, 
GGTCTA => -2.00943e-05, GGTCTC => 5.89509e-05, GGTCTG => 0.000298158, GGTCTT => 1.17874e-05, GGTGAA => 0.00016749, 
GGTGAC => 0.000201682, GGTGAG => 0.000226751, GGTGAT => 0.000215456, GGTGCA => 7.96452e-05, GGTGCC => 0.000267481, 
GGTGCG => 4.61294e-05, GGTGCT => 0.000163659, GGTGGA => 0.000273192, GGTGGC => 0.000385855, GGTGGG => -1.90666e-05, 
GGTGGT => 0.000226901, GGTGTA => -1.67472e-05, GGTGTC => 9.88804e-05, GGTGTG => 0.000185194, GGTGTT => 2.6391e-05, 
GGTTAA => -0.000197324, GGTTAC => 9.34548e-05, GGTTAG => -0.000147116, GGTTAT => 9.88067e-06, GGTTCA => -6.63248e-05, 
GGTTCC => 6.17997e-05, GGTTCG => 2.06238e-05, GGTTCT => -2.90442e-05, GGTTGA => -0.00014785, GGTTGC => -6.2847e-05, 
GGTTGG => -8.08971e-05, GGTTGT => -0.000128108, GGTTTA => -0.000115527, GGTTTC => 3.93631e-05, GGTTTG => -3.33523e-05, 
GGTTTT => -0.000237121, GTAAAA => -0.000333219, GTAAAC => -7.64174e-05, GTAAAG => 2.13946e-05, GTAAAT => -0.000291227, 
GTAACA => -0.000109358, GTAACC => -1.0598e-05, GTAACG => -3.45867e-05, GTAACT => -8.32767e-05, GTAAGA => -0.000130575, 
GTAAGC => -7.18705e-05, GTAAGG => -6.93555e-05, GTAAGT => -0.000123206, GTAATA => -0.000196692, GTAATC => -5.63005e-05, 
GTAATG => -6.79978e-05, GTAATT => -0.000235481, GTACAA => -0.000109769, GTACAC => -6.25354e-05, GTACAG => 4.28066e-05, 
GTACAT => -0.000165477, GTACCA => -9.90741e-06, GTACCC => 3.8192e-05, GTACCG => 1.339e-05, GTACCT => -1.70726e-05, 
GTACGA => -6.19456e-06, GTACGC => 2.1185e-05, GTACGG => -3.19492e-06, GTACGT => -3.36989e-05, GTACTA => -6.44142e-05, 
GTACTC => -3.20376e-05, GTACTG => 3.32388e-06, GTACTT => -0.000140341, GTAGAA => -1.06176e-05, GTAGAC => 6.5995e-05, 
GTAGAG => 0.000129391, GTAGAT => 3.31883e-05, GTAGCA => -2.29536e-05, GTAGCC => 6.3708e-05, GTAGCG => -1.92626e-06, 
GTAGCT => -3.23812e-05, GTAGGA => -3.87693e-05, GTAGGC => -1.77931e-05, GTAGGG => -3.76133e-05, GTAGGT => -4.76084e-05, 
GTAGTA => -6.41968e-05, GTAGTC => -2.05808e-05, GTAGTG => 1.63629e-05, GTAGTT => -0.000111571, GTATAA => -0.000253129, 
GTATAC => -5.83551e-05, GTATAG => -0.000127985, GTATAT => -0.000278981, GTATCA => -5.49452e-05, GTATCC => 2.37413e-05, 
GTATCG => -2.11628e-06, GTATCT => -8.22392e-05, GTATGA => -0.000162921, GTATGC => -7.24124e-05, GTATGG => -8.00405e-05, 
GTATGT => -0.000237569, GTATTA => -0.000215524, GTATTC => -5.73751e-05, GTATTG => -0.000117441, GTATTT => -0.000365108, 
GTCAAA => 0.000124486, GTCAAC => 0.000331845, GTCAAG => 0.000468772, GTCAAT => 0.0001273, GTCACA => 4.98468e-05, 
GTCACC => 0.000350664, GTCACG => 8.22556e-05, GTCACT => 8.82787e-05, GTCAGA => -6.61188e-05, GTCAGC => 0.00020993, 
GTCAGG => -4.13907e-05, GTCAGT => 3.15791e-05, GTCATA => -3.20264e-05, GTCATC => 0.000455194, GTCATG => 0.000219703, 
GTCATT => 0.000104758, GTCCAA => -1.93071e-05, GTCCAC => 4.89323e-05, GTCCAG => 0.000282163, GTCCAT => -2.99386e-05, 
GTCCCA => -3.28768e-05, GTCCCC => -2.4495e-06, GTCCCG => 1.54129e-05, GTCCCT => -2.09007e-05, GTCCGA => 8.68646e-06, 
GTCCGC => 6.53827e-05, GTCCGG => 3.75672e-05, GTCCGT => 9.68408e-06, GTCCTA => -2.2761e-05, GTCCTC => 6.42987e-05, 
GTCCTG => 0.000302312, GTCCTT => -2.57652e-05, GTCGAA => -5.93147e-06, GTCGAC => 9.17598e-05, GTCGAG => 0.000152776, 
GTCGAT => 8.02066e-05, GTCGCA => -1.29665e-05, GTCGCC => 8.08468e-05, GTCGCG => -2.14168e-05, GTCGCT => -5.71202e-06, 
GTCGGA => 2.66547e-05, GTCGGC => 8.8814e-05, GTCGGG => -2.55475e-05, GTCGGT => 4.11623e-05, GTCGTA => -7.44989e-06, 
GTCGTC => 5.76868e-05, GTCGTG => 7.81386e-05, GTCGTT => -1.05518e-05, GTCTAA => -0.000166885, GTCTAC => 0.000234934, 
GTCTAG => -0.000120261, GTCTAT => 4.3485e-05, GTCTCA => -5.25117e-05, GTCTCC => 0.000129579, GTCTCG => 4.14967e-05, 
GTCTCT => -5.9327e-05, GTCTGA => -0.000249599, GTCTGC => 2.44871e-05, GTCTGG => -4.38067e-05, GTCTGT => -0.00014699, 
GTCTTA => -7.89592e-05, GTCTTC => 0.000236959, GTCTTG => -1.45364e-05, GTCTTT => -3.87103e-05, GTGAAA => 9.11677e-06, 
GTGAAC => 0.000314223, GTGAAG => 0.000702709, GTGAAT => 9.70859e-05, GTGACA => 6.68173e-05, GTGACC => 0.000354459, 
GTGACG => 0.000144858, GTGACT => 7.76639e-05, GTGAGA => -7.44517e-05, GTGAGC => 0.000143587, GTGAGG => 3.45255e-06, 
GTGAGT => -2.8897e-05, GTGATA => -2.5463e-05, GTGATC => 0.000359812, GTGATG => 0.000321501, GTGATT => 0.00010217, 
GTGCAA => 2.98924e-05, GTGCAC => 0.000183483, GTGCAG => 0.000560087, GTGCAT => -6.38111e-05, GTGCCA => 0.000135206, 
GTGCCC => 0.000410807, GTGCCG => 0.000187992, GTGCCT => 9.83114e-05, GTGCGA => 8.94268e-05, GTGCGC => 0.000273788, 
GTGCGG => 0.000133246, GTGCGT => 4.41695e-05, GTGCTA => 2.91522e-05, GTGCTC => 0.000268422, GTGCTG => 0.000724732, 
GTGCTT => -9.17058e-06, GTGGAA => 0.000386814, GTGGAC => 0.000856684, GTGGAG => 0.00135974, GTGGAT => 0.000566357, 
GTGGCA => 0.000206165, GTGGCC => 0.000897033, GTGGCG => 0.000214046, GTGGCT => 0.000386774, GTGGGA => 0.000183608, 
GTGGGC => 0.000566722, GTGGGG => 2.04874e-05, GTGGGT => 0.000166876, GTGGTA => 4.37022e-05, GTGGTC => 0.000412508, 
GTGGTG => 0.000843607, GTGGTT => 0.000143537, GTGTAA => -0.000274851, GTGTAC => 0.00016894, GTGTAG => -0.000176778, 
GTGTAT => -5.62158e-05, GTGTCA => -2.8292e-06, GTGTCC => 0.000258644, GTGTCG => 0.000102855, GTGTCT => 7.80741e-05, 
GTGTGA => -0.000331506, GTGTGC => 1.26846e-05, GTGTGG => -1.49331e-05, GTGTGT => -0.000646837, GTGTTA => -0.000131924, 
GTGTTC => 0.000175265, GTGTTG => -2.86701e-05, GTGTTT => -8.68872e-05, GTTAAA => -0.000189627, GTTAAC => -1.61166e-05, 
GTTAAG => 2.99606e-05, GTTAAT => -0.000141102, GTTACA => -9.19628e-05, GTTACC => 4.22531e-05, GTTACG => -2.1933e-06, 
GTTACT => -4.46829e-05, GTTAGA => -9.1203e-05, GTTAGC => -1.34354e-05, GTTAGG => -6.7112e-05, GTTAGT => -0.000110283, 
GTTATA => -0.000117747, GTTATC => 1.86127e-05, GTTATG => -4.56391e-05, GTTATT => -0.000194845, GTTCAA => -0.000143194, 
GTTCAC => 1.3693e-05, GTTCAG => 0.000157874, GTTCAT => -8.91187e-05, GTTCCA => 2.39708e-05, GTTCCC => 5.52398e-05, 
GTTCCG => 3.99086e-05, GTTCCT => 4.13224e-05, GTTCGA => -3.8754e-05, GTTCGC => 4.88132e-05, GTTCGG => -3.07977e-06, 
GTTCGT => 9.91165e-06, GTTCTA => -8.01227e-05, GTTCTC => -2.97724e-05, GTTCTG => 0.000150633, GTTCTT => -0.00012419, 
GTTGAA => 1.64498e-05, GTTGAC => 2.46986e-05, GTTGAG => 0.000109122, GTTGAT => 2.95042e-05, GTTGCA => -5.3734e-05, 
GTTGCC => 6.60542e-05, GTTGCG => -2.24079e-05, GTTGCT => -1.51722e-05, GTTGGA => 7.00598e-05, GTTGGC => 9.72732e-05, 
GTTGGG => -4.63353e-05, GTTGGT => 1.37986e-05, GTTGTA => -0.000117969, GTTGTC => 1.29851e-05, GTTGTG => 6.54516e-05, 
GTTGTT => -0.00021826, GTTTAA => -0.000436309, GTTTAC => -2.39183e-05, GTTTAG => -0.000228904, GTTTAT => -0.000264469, 
GTTTCA => -0.00014213, GTTTCC => -7.00106e-05, GTTTCG => -2.80746e-05, GTTTCT => -0.000205455, GTTTGA => -0.00035369, 
GTTTGC => -0.000152757, GTTTGG => -0.00016515, GTTTGT => -0.000367138, GTTTTA => -0.000433697, GTTTTC => -0.000225473, 
GTTTTG => -0.000281709, GTTTTT => -0.000630574, TAAAAA => -0.0011481, TAAAAC => -0.000604573, TAAAAG => -0.000512056, 
TAAAAT => -0.000959441, TAAACA => -0.000581552, TAAACC => -0.000254276, TAAACG => -0.000175908, TAAACT => -0.000389711, 
TAAAGA => -0.00047006, TAAAGC => -0.000327949, TAAAGG => -0.00031342, TAAAGT => -0.000447195, TAAATA => -0.000873454, 
TAAATC => -0.000391638, TAAATG => -0.000579752, TAAATT => -0.000706588, TAACAA => -0.000396772, TAACAC => -0.000207995, 
TAACAG => -0.000258375, TAACAT => -0.000330125, TAACCA => -0.000246227, TAACCC => -0.000147587, TAACCG => -8.56285e-05, 
TAACCT => -0.000179049, TAACGA => -0.000114441, TAACGC => -8.41322e-05, TAACGG => -9.21856e-05, TAACGT => -0.00011373, 
TAACTA => -0.000247755, TAACTC => -0.000187333, TAACTG => -0.000239015, TAACTT => -0.000324706, TAAGAA => -0.0003679, 
TAAGAC => -0.000166526, TAAGAG => -0.000217918, TAAGAT => -0.000249055, TAAGCA => -0.000280598, TAAGCC => -0.000165071, 
TAAGCG => -9.06445e-05, TAAGCT => -0.00023884, TAAGGA => -0.000213576, TAAGGC => -0.000161442, TAAGGG => -0.000145974, 
TAAGGT => -0.000173189, TAAGTA => -0.000271458, TAAGTC => -0.000159634, TAAGTG => -0.000265377, TAAGTT => -0.000324412, 
TAATAA => -0.00077139, TAATAC => -0.000226557, TAATAG => -0.000232198, TAATAT => -0.000522392, TAATCA => -0.00031699, 
TAATCC => -0.00017345, TAATCG => -0.000102599, TAATCT => -0.000262118, TAATGA => -0.000363929, TAATGC => -0.000231246, 
TAATGG => -0.000217316, TAATGT => -0.000401749, TAATTA => -0.000546353, TAATTC => -0.000295254, TAATTG => -0.000341884, 
TAATTT => -0.000777646, TACAAA => -0.000124069, TACAAC => 0.000343481, TACAAG => 0.000515161, TACAAT => 4.53507e-05, 
TACACA => -8.90264e-05, TACACC => 0.000268349, TACACG => 0.00012286, TACACT => -3.80998e-05, TACAGA => -2.10121e-05, 
TACAGC => 0.000251852, TACAGG => 5.33998e-05, TACAGT => -1.67974e-05, TACATA => -0.00029015, TACATC => 0.000334732, 
TACATG => 0.000250899, TACATT => -0.000164645, TACCAA => -3.94923e-05, TACCAC => 0.000120644, TACCAG => 0.000407338, 
TACCAT => -5.04406e-06, TACCCA => 3.47041e-05, TACCCC => 6.27699e-05, TACCCG => 4.8133e-05, TACCCT => 2.53527e-05, 
TACCGA => 5.29931e-05, TACCGC => 0.000194578, TACCGG => 9.15748e-05, TACCGT => 4.85854e-05, TACCTA => -2.15064e-05, 
TACCTC => 0.000116805, TACCTG => 0.000490372, TACCTT => -5.32403e-05, TACGAA => 7.80086e-05, TACGAC => 0.000305592, 
TACGAG => 0.000458991, TACGAT => 0.000220868, TACGCA => 1.39151e-05, TACGCC => 0.000259108, TACGCG => 4.93362e-05, 
TACGCT => 6.03554e-05, TACGGA => 0.000119632, TACGGC => 0.000277775, TACGGG => 6.47936e-05, TACGGT => 6.39383e-05, 
TACGTA => -2.53321e-05, TACGTC => 0.000100784, TACGTG => 0.000205988, TACGTT => -2.37196e-05, TACTAA => -0.000196672, 
TACTAC => 0.000275363, TACTAG => -9.0168e-05, TACTAT => 2.92971e-05, TACTCA => -5.09235e-05, TACTCC => 0.000188398, 
TACTCG => 9.9043e-05, TACTCT => -1.15212e-05, TACTGA => -0.000217508, TACTGC => 0.000126062, TACTGG => 7.73504e-05, 
TACTGT => -9.05088e-05, TACTTA => -0.00018724, TACTTC => 0.000298002, TACTTG => 1.27846e-05, TACTTT => -0.000107467, 
TAGAAA => -0.00042038, TAGAAC => -0.000160145, TAGAAG => -0.000248878, TAGAAT => -0.000259483, TAGACA => -0.000197936, 
TAGACC => -0.000115806, TAGACG => -8.08734e-05, TAGACT => -0.000162684, TAGAGA => -0.000248196, TAGAGC => -0.000183919, 
TAGAGG => -0.000170103, TAGAGT => -0.000184303, TAGATA => -0.000204579, TAGATC => -0.000142337, TAGATG => -0.000235708, 
TAGATT => -0.000251978, TAGCAA => -0.000243319, TAGCAC => -0.000153588, TAGCAG => -0.000202265, TAGCAT => -0.000213987, 
TAGCCA => -0.000203811, TAGCCC => -0.000128733, TAGCCG => -7.38621e-05, TAGCCT => -0.000166769, TAGCGA => -7.85955e-05, 
TAGCGC => -6.8378e-05, TAGCGG => -6.10378e-05, TAGCGT => -6.72352e-05, TAGCTA => -0.000168057, TAGCTC => -0.000172649, 
TAGCTG => -0.000241179, TAGCTT => -0.000222566, TAGGAA => -0.000235839, TAGGAC => -0.000120169, TAGGAG => -0.000157705, 
TAGGAT => -0.000159797, TAGGCA => -0.00015718, TAGGCC => -9.4954e-05, TAGGCG => -5.1879e-05, TAGGCT => -0.000163353, 
TAGGGA => -0.000130989, TAGGGC => -0.000102229, TAGGGG => -0.000101461, TAGGGT => -0.000123664, TAGGTA => -0.000144503, 
TAGGTC => -9.78277e-05, TAGGTG => -0.000138542, TAGGTT => -0.000156941, TAGTAA => -0.000220111, TAGTAC => -0.000123848, 
TAGTAG => -0.000140082, TAGTAT => -0.000188834, TAGTCA => -0.000165219, TAGTCC => -0.000124003, TAGTCG => -7.08957e-05, 
TAGTCT => -0.000166313, TAGTGA => -0.000195044, TAGTGC => -0.000146244, TAGTGG => -0.000151867, TAGTGT => -0.000233486, 
TAGTTA => -0.000240822, TAGTTC => -0.000192573, TAGTTG => -0.000214521, TAGTTT => -0.000445836, TATAAA => -0.00048274, 
TATAAC => -5.79559e-05, TATAAG => -2.79573e-05, TATAAT => -0.000300216, TATACA => -0.000362149, TATACC => -3.71778e-05, 
TATACG => -1.07637e-05, TATACT => -0.000126661, TATAGA => -0.000205672, TATAGC => -7.37756e-05, TATAGG => -8.93141e-05, 
TATAGT => -0.000170984, TATATA => -0.000881797, TATATC => -8.45327e-05, TATATG => -0.00024502, TATATT => -0.000516037, 
TATCAA => -0.0001561, TATCAC => 7.03662e-06, TATCAG => 0.000125405, TATCAT => -0.000119546, TATCCA => -5.37202e-06, 
TATCCC => 6.69293e-05, TATCCG => 8.26635e-05, TATCCT => -5.20416e-06, TATCGA => -2.425e-05, TATCGC => 8.68769e-05, 
TATCGG => 2.05533e-05, TATCGT => -1.14884e-05, TATCTA => -0.000146544, TATCTC => -4.25932e-05, TATCTG => 0.000191154, 
TATCTT => -0.000151982, TATGAA => 2.46945e-05, TATGAC => 0.000262081, TATGAG => 0.000346626, TATGAT => 7.39004e-05, 
TATGCA => -5.87041e-05, TATGCC => 0.000251253, TATGCG => 5.46035e-05, TATGCT => 4.06246e-06, TATGGA => 0.000119558, 
TATGGC => 0.000178013, TATGGG => 2.65072e-05, TATGGT => -1.38781e-06, TATGTA => -0.000353905, TATGTC => 3.02733e-05, 
TATGTG => 0.00019158, TATGTT => -0.000188348, TATTAA => -0.000590815, TATTAC => -5.71516e-05, TATTAG => -0.000232931, 
TATTAT => -0.000484761, TATTCA => -0.000259319, TATTCC => -7.20131e-05, TATTCG => -5.55569e-05, TATTCT => -0.000215336, 
TATTGA => -0.000302551, TATTGC => -0.000144537, TATTGG => -0.000124093, TATTGT => -0.000389523, TATTTA => -0.000754087, 
TATTTC => -0.000192238, TATTTG => -0.000309108, TATTTT => -0.00107385, TCAAAA => -0.000417248, TCAAAC => -0.000123175, 
TCAAAG => -8.13273e-05, TCAAAT => -0.000277558, TCAACA => -0.000170793, TCAACC => -4.35716e-05, TCAACG => -1.38809e-05, 
TCAACT => -0.000158614, TCAAGA => -0.00019482, TCAAGC => -8.09077e-05, TCAAGG => -0.000108014, TCAAGT => -0.000134064, 
TCAATA => -0.00022162, TCAATC => -0.000112225, TCAATG => -0.000111797, TCAATT => -0.000273849, TCACAA => -0.000170242, 
TCACAC => -0.000171063, TCACAG => -7.04347e-06, TCACAT => -0.000210241, TCACCA => -5.35808e-05, TCACCC => 1.98423e-05, 
TCACCG => 2.46664e-05, TCACCT => -6.02171e-05, TCACGA => -2.20179e-05, TCACGC => 6.78412e-06, TCACGG => -2.34778e-05, 
TCACGT => -4.72282e-05, TCACTA => -9.58809e-05, TCACTC => -9.73226e-05, TCACTG => -1.01892e-05, TCACTT => -0.000199433, 
TCAGAA => -4.56272e-05, TCAGAC => 0.000124959, TCAGAG => 0.00020341, TCAGAT => 4.14427e-05, TCAGCA => -8.78785e-05, 
TCAGCC => 5.70377e-05, TCAGCG => -1.60104e-05, TCAGCT => -3.32468e-05, TCAGGA => -5.63882e-06, TCAGGC => 9.67855e-06, 
TCAGGG => -2.01285e-05, TCAGGT => -6.15477e-05, TCAGTA => -0.000109863, TCAGTC => -9.54829e-05, TCAGTG => -9.97484e-06, 
TCAGTT => -0.000236532, TCATAA => -0.000299431, TCATAC => -6.72068e-05, TCATAG => -0.000167159, TCATAT => -0.000252894, 
TCATCA => -0.000156403, TCATCC => 2.17017e-05, TCATCG => 2.83939e-05, TCATCT => -0.000112205, TCATGA => -0.000235333, 
TCATGC => -0.000129374, TCATGG => -0.000112266, TCATGT => -0.000236967, TCATTA => -0.000267527, TCATTC => -0.000127606, 
TCATTG => -0.000169094, TCATTT => -0.000462824, TCCAAA => 8.01892e-05, TCCAAC => 0.000288348, TCCAAG => 0.000489848, 
TCCAAT => 0.000155829, TCCACA => 1.54818e-05, TCCACC => 0.000305554, TCCACG => 0.000183445, TCCACT => 3.74709e-05, 
TCCAGA => -8.76936e-05, TCCAGC => 0.000317402, TCCAGG => -5.19919e-05, TCCAGT => 0.000132257, TCCATA => -6.72737e-05, 
TCCATC => 0.000288308, TCCATG => 0.000263071, TCCATT => 5.04668e-05, TCCCAA => -7.94072e-05, TCCCAC => -1.29033e-05, 
TCCCAG => 0.000184777, TCCCAT => -5.95514e-05, TCCCCA => -2.10119e-05, TCCCCC => -3.86313e-05, TCCCCG => 4.25794e-05, 
TCCCCT => -2.80884e-05, TCCCGA => -1.93599e-06, TCCCGC => 0.000112528, TCCCGG => 3.09699e-05, TCCCGT => 2.96364e-05, 
TCCCTA => -4.70218e-05, TCCCTC => -3.3555e-05, TCCCTG => 0.000246204, TCCCTT => -9.79779e-05, TCCGAA => 2.8867e-05, 
TCCGAC => 0.000145697, TCCGAG => 0.000253568, TCCGAT => 0.000120566, TCCGCA => 2.50371e-05, TCCGCC => 0.00018606, 
TCCGCG => 1.47433e-05, TCCGCT => 4.04058e-05, TCCGGA => 8.2179e-05, TCCGGC => 0.000159162, TCCGGG => 9.73286e-06, 
TCCGGT => 6.70843e-05, TCCGTA => -1.23078e-05, TCCGTC => 5.46123e-05, TCCGTG => 0.000169271, TCCGTT => 5.45148e-06, 
TCCTAA => -0.000191522, TCCTAC => 0.000219441, TCCTAG => -0.000120453, TCCTAT => 4.7695e-05, TCCTCA => 1.46929e-05, 
TCCTCC => 0.000170671, TCCTCG => 0.000145038, TCCTCT => -4.17052e-05, TCCTGA => -0.00031437, TCCTGC => -3.33424e-05, 
TCCTGG => -4.57162e-05, TCCTGT => -0.000120741, TCCTTA => -8.48427e-05, TCCTTC => 0.000137152, TCCTTG => 2.81612e-05, 
TCCTTT => -0.000135329, TCGAAA => -0.000118439, TCGAAC => 1.04981e-05, TCGAAG => 6.71148e-05, TCGAAT => -2.39043e-05, 
TCGACA => -2.73611e-05, TCGACC => 2.99597e-05, TCGACG => 3.42891e-05, TCGACT => -3.23301e-05, TCGAGA => -5.95502e-05, 
TCGAGC => 1.72373e-05, TCGAGG => -2.81984e-05, TCGAGT => -4.57029e-05, TCGATA => -5.10946e-05, TCGATC => 1.04092e-07, 
TCGATG => 5.42264e-05, TCGATT => -6.42019e-05, TCGCAA => 1.02194e-06, TCGCAC => 5.98611e-05, TCGCAG => 0.000230923, 
TCGCAT => -5.27236e-06, TCGCCA => 8.14453e-05, TCGCCC => 0.000182616, TCGCCG => 0.000138147, TCGCCT => 1.01637e-05, 
TCGCGA => 1.33131e-05, TCGCGC => 6.60996e-05, TCGCGG => 5.07743e-06, TCGCGT => -1.30104e-05, TCGCTA => 2.48645e-06, 
TCGCTC => 3.26436e-05, TCGCTG => 0.000334187, TCGCTT => -7.03929e-05, TCGGAA => 8.28781e-05, TCGGAC => 0.000164115, 
TCGGAG => 0.000309343, TCGGAT => 0.000142888, TCGGCA => 5.93125e-05, TCGGCC => 0.000178643, TCGGCG => 0.000133709, 
TCGGCT => 6.26422e-05, TCGGGA => 8.01594e-05, TCGGGC => 0.000102002, TCGGGG => 6.29257e-06, TCGGGT => 3.87791e-05, 
TCGGTA => -1.17002e-05, TCGGTC => 2.68537e-05, TCGGTG => 0.000242874, TCGGTT => -3.52386e-05, TCGTAA => -0.000107807, 
TCGTAC => 6.71646e-05, TCGTAG => -6.38478e-05, TCGTAT => -1.37903e-05, TCGTCA => -6.20279e-06, TCGTCC => 4.70324e-05, 
TCGTCG => 7.02809e-05, TCGTCT => -1.5391e-05, TCGTGA => -8.27723e-05, TCGTGC => -1.79561e-05, TCGTGG => -5.2292e-06, 
TCGTGT => -6.63648e-05, TCGTTA => -7.54275e-05, TCGTTC => 1.50101e-05, TCGTTG => 4.38961e-06, TCGTTT => -0.000125537, 
TCTAAA => -0.000170166, TCTAAC => -3.27999e-05, TCTAAG => 1.34651e-05, TCTAAT => -0.00012597, TCTACA => -8.68025e-05, 
TCTACC => -1.90401e-05, TCTACG => -1.7533e-05, TCTACT => -6.34634e-05, TCTAGA => -0.000109533, TCTAGC => -1.23431e-05, 
TCTAGG => -8.41739e-05, TCTAGT => -7.00044e-05, TCTATA => -0.000154232, TCTATC => -5.74309e-05, TCTATG => -4.38933e-05, 
TCTATT => -0.000170404, TCTCAA => -0.00010817, TCTCAC => -5.66325e-05, TCTCAG => 0.000127249, TCTCAT => -0.00014629, 
TCTCCA => 2.27284e-05, TCTCCC => -3.75147e-05, TCTCCG => 1.60986e-05, TCTCCT => -3.35391e-05, TCTCGA => 2.41962e-06, 
TCTCGC => 2.4302e-05, TCTCGG => -1.79585e-06, TCTCGT => -3.74178e-05, TCTCTA => -0.0001104, TCTCTC => -0.000268858, 
TCTCTG => 0.000136566, TCTCTT => -0.000225515, TCTGAA => 1.68982e-06, TCTGAC => 0.000109354, TCTGAG => 0.000221285, 
TCTGAT => 6.90957e-05, TCTGCA => -2.43353e-05, TCTGCC => 8.72656e-05, TCTGCG => -1.03792e-05, TCTGCT => -2.6239e-05, 
TCTGGA => 6.34431e-05, TCTGGC => 8.23082e-05, TCTGGG => -5.34102e-05, TCTGGT => -6.49242e-06, TCTGTA => -0.000183578, 
TCTGTC => -9.11728e-05, TCTGTG => 9.02706e-05, TCTGTT => -0.000196284, TCTTAA => -0.000374593, TCTTAC => -4.31607e-05, 
TCTTAG => -0.000175044, TCTTAT => -0.000153721, TCTTCA => -0.000174349, TCTTCC => -7.58406e-05, TCTTCG => 3.2735e-06, 
TCTTCT => -0.000173979, TCTTGA => -0.000272809, TCTTGC => -0.000138858, TCTTGG => -0.000146373, TCTTGT => -0.000255389, 
TCTTTA => -0.000292076, TCTTTC => -0.00020122, TCTTTG => -0.0001809, TCTTTT => -0.00048536, TGAAAA => -0.000774813, 
TGAAAC => -0.000375077, TGAAAG => -0.000416048, TGAAAT => -0.000624474, TGAACA => -0.000408411, TGAACC => -0.000203399, 
TGAACG => -0.000135254, TGAACT => -0.0003679, TGAAGA => -0.00053, TGAAGC => -0.000329961, TGAAGG => -0.000328623, 
TGAAGT => -0.000361765, TGAATA => -0.000400464, TGAATC => -0.000263842, TGAATG => -0.000448006, TGAATT => -0.000468346, 
TGACAA => -0.000336145, TGACAC => -0.000205897, TGACAG => -0.000318215, TGACAT => -0.000310308, TGACCA => -0.00025749, 
TGACCC => -0.00018145, TGACCG => -8.44655e-05, TGACCT => -0.000252302, TGACGA => -0.000122569, TGACGC => -0.000108305, 
TGACGG => -9.34465e-05, TGACGT => -0.000130741, TGACTA => -0.000197098, TGACTC => -0.000225622, TGACTG => -0.000327001, 
TGACTT => -0.000347053, TGAGAA => -0.000400046, TGAGAC => -0.000239719, TGAGAG => -0.000323382, TGAGAT => -0.000294572, 
TGAGCA => -0.000320974, TGAGCC => -0.000225737, TGAGCG => -0.000125196, TGAGCT => -0.000315173, TGAGGA => -0.000343576, 
TGAGGC => -0.000220574, TGAGGG => -0.00021761, TGAGGT => -0.000223878, TGAGTA => -0.000219829, TGAGTC => -0.000196905, 
TGAGTG => -0.000299309, TGAGTT => -0.000335735, TGATAA => -0.000340134, TGATAC => -0.000162537, TGATAG => -0.000169307, 
TGATAT => -0.000331349, TGATCA => -0.000259396, TGATCC => -0.000181226, TGATCG => -7.94892e-05, TGATCT => -0.000255539, 
TGATGA => -0.000428417, TGATGC => -0.000245069, TGATGG => -0.000294508, TGATGT => -0.000359595, TGATTA => -0.000326044, 
TGATTC => -0.000259638, TGATTG => -0.000295562, TGATTT => -0.000611694, TGCAAA => -0.00021605, TGCAAC => 0.00012052, 
TGCAAG => 0.000193011, TGCAAT => -9.3841e-05, TGCACA => -0.000163827, TGCACC => 0.000103767, TGCACG => 1.89431e-05, 
TGCACT => -0.000103974, TGCAGA => -0.000239523, TGCAGC => 2.83438e-06, TGCAGG => -0.000150441, TGCAGT => -0.000144575, 
TGCATA => -0.000187336, TGCATC => 0.000115496, TGCATG => -1.19308e-05, TGCATT => -0.000233314, TGCCAA => -0.000171576, 
TGCCAC => -3.29679e-05, TGCCAG => 0.00015345, TGCCAT => -0.00012379, TGCCCA => -0.00010249, TGCCCC => -3.86125e-10, 
TGCCCG => 6.38308e-06, TGCCCT => -0.000107092, TGCCGA => -3.78196e-05, TGCCGC => 6.48938e-05, TGCCGG => 4.5253e-05, 
TGCCGT => -3.47619e-05, TGCCTA => -7.91184e-05, TGCCTC => -6.40645e-05, TGCCTG => 0.000107668, TGCCTT => -0.000192758, 
TGCGAA => -7.551e-07, TGCGAC => 0.000104163, TGCGAG => 0.000210825, TGCGAT => 7.11899e-05, TGCGCA => -7.67342e-05, 
TGCGCC => 7.66989e-05, TGCGCG => -4.44331e-05, TGCGCT => -4.14445e-05, TGCGGA => 2.92597e-05, TGCGGC => 0.000108452, 
TGCGGG => -4.08249e-06, TGCGGT => -9.19358e-06, TGCGTA => -5.30702e-05, TGCGTC => 2.43955e-05, TGCGTG => 4.15397e-05, 
TGCGTT => -8.37267e-05, TGCTAA => -0.000268922, TGCTAC => 8.18359e-05, TGCTAG => -0.000150085, TGCTAT => -5.20216e-05, 
TGCTCA => -0.00014919, TGCTCC => 3.79722e-05, TGCTCG => -1.57325e-05, TGCTCT => -0.000184399, TGCTGA => -0.000397024, 
TGCTGC => -0.000131556, TGCTGG => -0.000195227, TGCTGT => -0.000260752, TGCTTA => -0.000198556, TGCTTC => 3.75933e-05, 
TGCTTG => -0.000120753, TGCTTT => -0.000328639, TGGAAA => -0.000269923, TGGAAC => 3.25089e-05, TGGAAG => 7.82256e-05, 
TGGAAT => -0.000103977, TGGACA => -0.000118073, TGGACC => 5.95572e-05, TGGACG => -1.01754e-05, TGGACT => -0.000128976, 
TGGAGA => -0.00028236, TGGAGC => -0.000102726, TGGAGG => -0.000170414, TGGAGT => -0.000142083, TGGATA => -0.000167287, 
TGGATC => 2.64823e-05, TGGATG => -1.38459e-05, TGGATT => -0.000159355, TGGCAA => -0.00020519, TGGCAC => -4.91819e-05, 
TGGCAG => 2.71645e-05, TGGCAT => -0.000126801, TGGCCA => -0.000196863, TGGCCC => -7.84206e-05, TGGCCG => -5.253e-05, 
TGGCCT => -0.000169586, TGGCGA => -4.30267e-05, TGGCGC => 2.6754e-05, TGGCGG => -3.54201e-05, TGGCGT => -3.6645e-05, 
TGGCTA => -9.0309e-05, TGGCTC => -3.84778e-05, TGGCTG => 0.000111915, TGGCTT => -0.000214331, TGGGAA => -0.000108281, 
TGGGAC => 0.000108151, TGGGAG => 0.000164221, TGGGAT => 4.5359e-05, TGGGCA => -7.98183e-05, TGGGCC => 4.10159e-05, 
TGGGCG => -1.79015e-05, TGGGCT => -6.91615e-05, TGGGGA => -9.81962e-05, TGGGGC => 3.3397e-05, TGGGGG => -0.000140839, 
TGGGGT => -6.49273e-05, TGGGTA => -9.39546e-05, TGGGTC => -2.9025e-05, TGGGTG => 4.63609e-05, TGGGTT => -0.000161233, 
TGGTAA => -0.000220336, TGGTAC => 6.90141e-05, TGGTAG => -0.000149001, TGGTAT => -6.27032e-05, TGGTCA => -0.000132107, 
TGGTCC => -3.58107e-05, TGGTCG => -2.60446e-06, TGGTCT => -0.000115892, TGGTGA => -0.000263281, TGGTGC => -8.77855e-05, 
TGGTGG => -0.000108713, TGGTGT => -0.000209359, TGGTTA => -0.00019032, TGGTTC => 2.68775e-05, TGGTTG => -0.000154505, 
TGGTTT => -0.000238385, TGTAAA => -0.000447255, TGTAAC => -0.000100863, TGTAAG => -9.60747e-05, TGTAAT => -0.000322315, 
TGTACA => -0.000265305, TGTACC => -7.53704e-05, TGTACG => -6.00994e-05, TGTACT => -0.000182881, TGTAGA => -0.000217397, 
TGTAGC => -0.000144405, TGTAGG => -0.000139173, TGTAGT => -0.000181352, TGTATA => -0.000390098, TGTATC => -0.00010106, 
TGTATG => -0.000228356, TGTATT => -0.000429714, TGTCAA => -0.000185885, TGTCAC => -0.000115205, TGTCAG => -3.70635e-05, 
TGTCAT => -0.000222614, TGTCCA => -0.000122604, TGTCCC => -4.23575e-05, TGTCCG => -2.38359e-06, TGTCCT => -0.000136073, 
TGTCGA => -4.09819e-05, TGTCGC => 8.27457e-06, TGTCGG => -1.14563e-05, TGTCGT => -6.18752e-05, TGTCTA => -0.000156294, 
TGTCTC => -0.000160958, TGTCTG => -0.000101767, TGTCTT => -0.000282994, TGTGAA => -0.000157643, TGTGAC => 7.18128e-05, 
TGTGAG => 0.000106349, TGTGAT => -7.74688e-05, TGTGCA => -0.000203628, TGTGCC => 7.67355e-05, TGTGCG => -8.98324e-05, 
TGTGCT => -0.00013049, TGTGGA => -7.53823e-05, TGTGGC => 3.30029e-05, TGTGGG => -6.013e-05, TGTGGT => -0.00011976, 
TGTGTA => -0.000337102, TGTGTC => -0.000116412, TGTGTG => -0.000619128, TGTGTT => -0.000432541, TGTTAA => -0.000424988, 
TGTTAC => -0.000108898, TGTTAG => -0.000208566, TGTTAT => -0.00029594, TGTTCA => -0.000319366, TGTTCC => -0.000113407, 
TGTTCG => -7.31095e-05, TGTTCT => -0.000304023, TGTTGA => -0.000394536, TGTTGC => -0.00022765, TGTTGG => -0.000266066, 
TGTTGT => -0.00044194, TGTTTA => -0.000549353, TGTTTC => -0.000278248, TGTTTG => -0.000435644, TGTTTT => -0.000997099, 
TTAAAA => -0.000962354, TTAAAC => -0.000311848, TTAAAG => -0.000315728, TTAAAT => -0.000719407, TTAACA => -0.000299078, 
TTAACC => -0.000116621, TTAACG => -9.20596e-05, TTAACT => -0.000260134, TTAAGA => -0.000263598, TTAAGC => -0.000184134, 
TTAAGG => -0.000164565, TTAAGT => -0.000304991, TTAATA => -0.000517799, TTAATC => -0.000209063, TTAATG => -0.000283923, 
TTAATT => -0.000640395, TTACAA => -0.000302917, TTACAC => -0.000170293, TTACAG => -0.000102928, TTACAT => -0.000335646, 
TTACCA => -0.000120091, TTACCC => -6.12418e-05, TTACCG => -3.15079e-05, TTACCT => -0.000119036, TTACGA => -7.2786e-05, 
TTACGC => -2.85294e-05, TTACGG => -2.48037e-05, TTACGT => -9.03999e-05, TTACTA => -0.00017366, TTACTC => -0.000123615, 
TTACTG => -9.62898e-05, TTACTT => -0.000293648, TTAGAA => -0.000126898, TTAGAC => -1.78703e-05, TTAGAG => -1.2857e-05, 
TTAGAT => -9.75671e-05, TTAGCA => -0.000138058, TTAGCC => -3.85536e-05, TTAGCG => -4.75322e-05, TTAGCT => -0.000135134, 
TTAGGA => -0.000106056, TTAGGC => -6.5707e-05, TTAGGG => -8.50402e-05, TTAGGT => -9.794e-05, TTAGTA => -0.000175438, 
TTAGTC => -9.9281e-05, TTAGTG => -0.000110599, TTAGTT => -0.00033463, TTATAA => -0.000515752, TTATAC => -0.000204304, 
TTATAG => -0.000256421, TTATAT => -0.000581799, TTATCA => -0.000242155, TTATCC => -9.32498e-05, TTATCG => -6.05512e-05, 
TTATCT => -0.000201804, TTATGA => -0.000328501, TTATGC => -0.000166461, TTATGG => -0.000169531, TTATGT => -0.000375785, 
TTATTA => -0.000651655, TTATTC => -0.00025605, TTATTG => -0.000345171, TTATTT => -0.00110642, TTCAAA => -0.000225372, 
TTCAAC => 0.000340553, TTCAAG => 0.000457948, TTCAAT => -1.82502e-05, TTCACA => -5.93191e-05, TTCACC => 0.000354114, 
TTCACG => 8.40084e-05, TTCACT => -6.8487e-06, TTCAGA => -0.000161808, TTCAGC => 0.000224878, TTCAGG => -5.73069e-05, 
TTCAGT => -5.54392e-05, TTCATA => -0.000182929, TTCATC => 0.000411793, TTCATG => 0.000198841, TTCATT => -0.000166874, 
TTCCAA => -0.000113239, TTCCAC => 0.000100925, TTCCAG => 0.000466536, TTCCAT => -0.000109089, TTCCCA => -0.000110276, 
TTCCCC => 6.79226e-05, TTCCCG => 1.73902e-05, TTCCCT => -5.64151e-05, TTCCGA => 3.42271e-05, TTCCGC => 0.000215897, 
TTCCGG => 0.000133726, TTCCGT => 4.84525e-05, TTCCTA => -6.10379e-05, TTCCTC => 0.000179975, TTCCTG => 0.000556138, 
TTCCTT => -0.000155634, TTCGAA => -3.86751e-05, TTCGAC => 0.000246744, TTCGAG => 0.000359076, TTCGAT => 0.000200427, 
TTCGCA => -1.61589e-05, TTCGCC => 0.000277021, TTCGCG => -3.95066e-06, TTCGCT => 1.19598e-05, TTCGGA => 7.17819e-05, 
TTCGGC => 0.000188185, TTCGGG => 5.80451e-05, TTCGGT => 1.6055e-05, TTCGTA => -3.81999e-05, TTCGTC => 9.92942e-05, 
TTCGTG => 0.000209482, TTCGTT => -8.00322e-05, TTCTAA => -0.000327763, TTCTAC => 0.000316486, TTCTAG => -0.000200968, 
TTCTAT => -1.44242e-05, TTCTCA => -0.000135348, TTCTCC => 0.000193389, TTCTCG => 3.99166e-05, TTCTCT => -0.000160398, 
TTCTGA => -0.000376307, TTCTGC => 1.16288e-05, TTCTGG => -2.24567e-05, TTCTGT => -0.000274875, TTCTTA => -0.000244594, 
TTCTTC => 0.00031091, TTCTTG => -9.61561e-05, TTCTTT => -0.000318164, TTGAAA => -0.000303033, TTGAAC => -4.07765e-05, 
TTGAAG => 7.82151e-05, TTGAAT => -0.000221235, TTGACA => -0.000148336, TTGACC => 2.8195e-05, TTGACG => -1.41402e-05, 
TTGACT => -0.000134317, TTGAGA => -0.000188419, TTGAGC => -3.66686e-06, TTGAGG => -8.92947e-05, TTGAGT => -0.000153805, 
TTGATA => -0.000204076, TTGATC => -2.00314e-06, TTGATG => -0.000107192, TTGATT => -0.000297989, TTGCAA => -0.00015138, 
TTGCAC => -4.86148e-05, TTGCAG => 9.10904e-05, TTGCAT => -0.00021708, TTGCCA => -5.05284e-05, TTGCCC => 7.84557e-05, 
TTGCCG => 2.79249e-05, TTGCCT => -9.71673e-05, TTGCGA => -2.33633e-05, TTGCGC => 4.46235e-05, TTGCGG => 1.01294e-05, 
TTGCGT => 3.79248e-06, TTGCTA => -0.000126108, TTGCTC => -5.3236e-05, TTGCTG => 0.000100689, TTGCTT => -0.00025034, 
TTGGAA => 6.71948e-05, TTGGAC => 0.000233717, TTGGAG => 0.000375223, TTGGAT => 0.000211712, TTGGCA => -3.99173e-05, 
TTGGCC => 0.000268012, TTGGCG => 5.52117e-05, TTGGCT => 5.94365e-05, TTGGGA => -2.65113e-05, TTGGGC => 0.000156788, 
TTGGGG => -8.66696e-05, TTGGGT => 1.23692e-05, TTGGTA => -8.92782e-05, TTGGTC => 2.06953e-05, TTGGTG => 0.0001415, 
TTGGTT => -0.000191078, TTGTAA => -0.000485, TTGTAC => -6.83946e-05, TTGTAG => -0.000241649, TTGTAT => -0.000319742, 
TTGTCA => -0.000181526, TTGTCC => -1.71664e-05, TTGTCG => -4.65661e-06, TTGTCT => -0.000150903, TTGTGA => -0.000356992, 
TTGTGC => -0.000168348, TTGTGG => -0.000168933, TTGTGT => -0.000418629, TTGTTA => -0.000326058, TTGTTC => -0.000185844, 
TTGTTG => -0.00032268, TTGTTT => -0.000845791, TTTAAA => -0.000931805, TTTAAC => -0.000167313, TTTAAG => -0.000132264, 
TTTAAT => -0.000614623, TTTACA => -0.000404296, TTTACC => -8.4655e-05, TTTACG => -2.97584e-05, TTTACT => -0.000310679, 
TTTAGA => -0.000273057, TTTAGC => -0.000131107, TTTAGG => -0.000187128, TTTAGT => -0.000298336, TTTATA => -0.000592386, 
TTTATC => -0.000133852, TTTATG => -0.000294011, TTTATT => -0.000984052, TTTCAA => -0.000404833, TTTCAC => -0.000165284, 
TTTCAG => -4.7353e-05, TTTCAT => -0.000370435, TTTCCA => -0.000275636, TTTCCC => -0.000139086, TTTCCG => -3.54067e-05, 
TTTCCT => -0.000258722, TTTCGA => -8.81537e-05, TTTCGC => -3.11539e-05, TTTCGG => -5.31082e-05, TTTCGT => -0.000103604, 
TTTCTA => -0.000277738, TTTCTC => -0.000236986, TTTCTG => -0.000123285, TTTCTT => -0.000613839, TTTGAA => -8.34812e-05, 
TTTGAC => 0.000224637, TTTGAG => 0.000374008, TTTGAT => 7.3855e-05, TTTGCA => -0.000144449, TTTGCC => 0.000213853, 
TTTGCG => 1.23822e-05, TTTGCT => -6.45814e-05, TTTGGA => 6.06571e-05, TTTGGC => 0.000169539, TTTGGG => -5.03031e-05, 
TTTGGT => -5.67304e-05, TTTGTA => -0.000416051, TTTGTC => -4.42401e-05, TTTGTG => 0.000194442, TTTGTT => -0.000665815, 
TTTTAA => -0.00121922, TTTTAC => -0.000260428, TTTTAG => -0.000485603, TTTTAT => -0.00090895, TTTTCA => -0.000542466, 
TTTTCC => -0.000341289, TTTTCG => -0.000136383, TTTTCT => -0.000668613, TTTTGA => -0.000643246, TTTTGC => -0.00032112, 
TTTTGG => -0.000377112, TTTTGT => -0.000835562, TTTTTA => -0.00110572, TTTTTC => -0.00055282, TTTTTG => -0.000662522, 
TTTTTT => -0.0021966, ); #hexcodencode
}


__END__


=item about tests


../cdsqual.pl ../fly_Hexamer.tsv evg24m2banofun5ht1.cds

** add cdsqual utrbad/poor, aasize/complete + CF:class + CH:class for consensus code/nocode/unkn coding class
-- how much agreement?
-- cpat hex score differs, but ~ same dir, size; fickett scores agree, same calc

>> add in as debug/test routine to compare this vs cpat results; note cpat upcases all ids, fix also
>> what are cpat coding-prob cutoffs for Code/Unkn/Nonc calls? Unkn now too broad? (0.10..0.50)
  .. paper says cp.cut=0.36, from sens/spec intersection for human code/noncode data; use 0.20..0.40 as Unkn range?
  
perl -ne \
'@v=split; if(@v==6) { ($id,$tw,$cw,$pfk,$phx,$cp)=@v;
$cl=($cp>0.5)?"Code":($cp<0.1)?"Nonc":"Unkn";  
map{s/(\.\d\d\d\d\d)\d+/$1/; } ($cp,$phx);
$cpat{uc($id)}="CY:$cl/$cp,$pfk,$phx"; } 
else { ($id,$cw,$gp,$aw,$tw,$cf,$ch,$off,$oid)=@v; $cpat=$cpat{uc($id)}||0; 
map{ s/Noncode/Nonc/; s/Unknown/Unkn/; } ($cf,$ch); 
print join("\t",$id,$cw,$tw,$cf,$ch,$cpat,$aw)."\n"; } ' \
 evg24m2banofun.cpat.out evg24m2banofun5ht1.cds.qual  > evg24m2banofun5ht1.both.tab

mRNA_size       ORF_size        coding_prob                     0       Hexamer_score
Anofunz4kEVm006490t1    408     0       CF:Unkn/0.7862  CH:Nonc/-0.004488       0       408,na
Anofunz4kEVm000439t1    5427    13367   CF:Unkn/0.8559  CH:Code/0.3914  CY:Code/0.99999,0.8559,0.47519  1808,40%,complete-utrbad
Anofunz4kEVm000135t1    8100    11820   CF:Code/1.1192  CH:Code/0.7038  CY:Code/1,1.1192,0.55240        2699,68%,complete
Anofunz4kEVm000019t1    14598   17031   CF:Unkn/0.7945  CH:Code/0.5584  CY:Code/1,0.7945,0.21707        4865,85%,complete
Anofunz4kEVm000875t1    4149    10626   CF:Code/1.1277  CH:Code/0.3677  CY:Code/0.99999,1.1277,0.62489  1382,39%,complete-utrpoor
Anofunz4kEVm000120t1    8430    10498   CF:Code/0.9779  CH:Code/0.5646  CY:Code/1,0.9779,0.45377        2809,80%,complete
Anofunz4kEVm000476t1    5265    9975    CF:Unkn/0.9025  CH:Code/0.1629  CY:Code/1,0.9025,0.13082        1754,52%,complete-utrpoor
Anofunz4kEVm017181t1    384     666     CF:Nonc/0.723   CH:Nonc/-0.007536       CY:Unkn/0.12030,0.723,-0.33895  127,3%,complete-utrbad
Anofunz4kEVm002749t1    2214    9792    CF:Code/0.9778  CH:Code/0.07114 CY:Code/0.97670,0.9778,0.10541  737,22%,complete-utrbad
Anofunz4kEVm013135t1    492     2167    CF:Nonc/0.6609  CH:Nonc/-0.05008        CY:Nonc/0.00979,0.6609,-0.90751 163,5%,complete-utrbad
Anofunz4kEVm012229t1    537     1619    CF:Nonc/0.6719  CH:Nonc/-0.01677        CY:Nonc/0.09721,0.6719,-0.40097 178,5%,complete-utrbad
Anofunz4kEVm002138t1    2646    9633    CF:Unkn/0.8413  CH:Code/0.09633 CY:Code/0.99953,0.8413,0.23581  881,27%,complete-utrbad
Anofunz4kEVm016335t1    399     3763    CF:Nonc/0.4253  CH:Nonc/-0.0184 CY:Nonc/0.00054,0.4253,-0.57007 132,4%,complete-utrbad
Anofunz4kEVm000206t1    7065    9235    CF:Code/0.9582  CH:Code/0.5645  CY:Code/1,0.9582,0.52768        2354,76%,complete
Anofunz4kEVm016672t1    393     1222    CF:Nonc/0.5295  CH:Nonc/-0.03788        CY:Nonc/0.01442,0.5295,-0.74207 130,4%,complete-utrbad

=cut

# evigene/scripts/cdna_proteins.pm
sub proteindoc
{
  my($orf, $cdnalen, $cdnarev, $forFASTA) = @_;
  my($aalen, $compl, $pcds, $istop, $Selc) = (0) x 10;
  my( $orfprot, $prostart, $proend, $orflen, $orient) = orfParts($orf);
  $cdnarev ||= $orient; # fill in blank; use always? shouldnt need to pass cdnarev as param.
    # Ooops, orf->{start,stop} are reversed for -strand; start>stop : ($lb,$le)=($le,$lb) if($lb>$le);
  if($orfprot) {
    $pcds  = ($cdnalen>0 && $orflen>0) ? int(100*$orflen/$cdnalen) : 0;
    my $urev= ($prostart>$proend)?1:0;
    my $u1len= ($urev) ? $cdnalen - $prostart : $prostart - 1; 
    my $u2len= ($urev) ? $proend - 1 : $cdnalen - $proend;
    
    $aalen= length($orfprot); # $bestorf->{length}; # this is cds-len
    if(substr($orfprot,-1,1) eq '*') { $aalen--; if($NoStopCodon) { $orfprot =~ s/\*$//; } }
    $istop= $orf->{innerstop} || 0; # add 201403:
    $compl= $orf->{complete};
    $compl= ($compl==3)?"complete":($compl==2)?"partial5":($compl==1)?"partial3":"partial";
    
    ##? not bad if partial? if u1len or u2len == 0
    ## need global params for utrbad/poor
    
    if($cdnalen - $orflen <= $MINUTR) { } # ignore pcds if utr small
    elsif($pcds < $pCDSbad or $u1len > $orflen or $u2len > $orflen) { $compl.="-utrbad"; }  
    elsif($pcds < $pCDSpoor) { $compl.="-utrpoor";  } #?? maybe change to use prostart OR protend3 > 35%? 40% ?
    ##? add istop flag to compl ??

    # 2014.12 add Selc flags 
    # Funhe2Exx11m009903t6 aalen=556,80%,complete,selcstop; Selcstop=index; .. Name=Selenoprotein N 
    if($USESelenocysteine and index($orfprot,'u')>0) { 
      $compl.=",selcstop"; 
      my $isel= $orf->{Selc}||1; $Selc="Selcstop=$isel"; # want mRNA/CDS index * 1-origin,may be list: 123,456,888
    }
    
    if($forFASTA) { $orfprot =~ s/(.{60})/$1\n/g; } #? only for fasta output
  } else {
    $orfprot="X"; # ? make dummy orfprot?
  }
  
  my $fahead= "aalen=$aalen,$pcds%,$compl; clen=$cdnalen; strand=$cdnarev; offs=$prostart-$proend;";
  $fahead .= " $Selc;" if($Selc);
  $fahead .= " innerstop=$istop;" if($istop>0);
  if(my $orflags= $orf->{flags}) { $fahead .= " orflags=$orflags;"; }
  return($aalen,$pcds,$compl,$orflen,$fahead,$orfprot);
}
