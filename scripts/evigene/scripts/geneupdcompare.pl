#!/usr/bin/env perl
# geneupdcompare.pl

=item about

  geneupdcompare.pl -in geneafter.gff -over genebefore.gff
  
  - compare updated gene models by shared ID, e.g. pasa updates
  - report exon,CDS changes in length, exon number/location
  
#   count code    desc   Pasa updates to handle
#    7020 13      EST assembly extends UTRs.     << change exons
#    1922 14      EST assembly alters protein sequence, passes validation. < change CDS,exons,prot
#    1140 16      EST assembly properly stitched into gene structure.   < can change CDS,exons

# .. maybe this or not: joins, esp a.all same CDS models at one locus vs b. join of separate loci
#     424 36      EST-assembly found capable of merging multiple genes.  < replace mRNA,CDS,exons old gene recs



=item author
  
  don gilbert, gilbertd near indiana edu, ca 2010-2011
  part of EvidentialGene, evigene/scripts/
  see similar.  evigene/scripts/overgenedup.pl
  
=cut

use strict;
use warnings;
use Getopt::Long;

#my $UTRSLOP = 49; # <=  when samecds, UTR difference < this is same gene
#  my $UTRSLOPN = 0; # UTR exon count diff
#my $TINYCDS = 40; # ignore tiny CDS tests for typeover == CDS
my $EXONSLOP = 0;
my $SAME_CDS= 90; # ** option
#my $SHOW_CDSUTR= 0; # similar only ?
my $MINID_CDS= 0;
my $MINID_UTR= 0;
#my $EQUAL_CUT= 66; # ** option for summary

my $mrnatypes='mRNA';
my $exontypes='exon|CDS';
my $passtypes="";

#my $overlaplist= {};
my %overgenes= (); ## my $overgeneid= 0;

my($input,$overlaps,$debug,$idequalpatt,$skippatt,$ok);

$idequalpatt='ev\w+'; # match puevd3b024410t1 to evd3b024410t1 in hash

my $optok= GetOptions(
  "overlaps=s", \$overlaps, 
  "input=s", \$input,  
  "idequal=s", \$idequalpatt,
  "skippatt=s", \$skippatt,
#   "typeover=s", \$typeover, 
#   "exontypes=s", \$exontypes, 
#   "action=s", \$action, 
#   "mark=s", \$mark, 
#   "passtypes=s", \$passtypes,  #??
#   "slopexon=i", \$EXONSLOP,  
#   "sloputr=i", \$UTRSLOP,  
  "samecds=i", \$SAME_CDS,  
  "mincds=i", \$MINID_CDS,  
  "minutr=i", \$MINID_UTR,  
#   "SELFKEEP!", \$SELFKEEP,
#   "summary:s", \$dosum,
  "debug!", \$debug, 
  );

die "usage:
  geneupdcompare.pl -in geneafter.gff -over genebefore.gff

  ** ASSUMES gene model is mRNA > exon > CDS
  ** ASSUMES input gff is ordered by gene records (mRNA/exon/CDS all together per ID)
" unless($optok and $input and $overlaps);


my $ovh; 
   if($overlaps =~ /.gz$/) { $ok= open(OVR,"gunzip -c $overlaps |");  $ovh= *OVR; }
elsif($overlaps =~ /^(stdin|-)/) { $ovh= *STDIN; $ok=1; }
else { $ok= open(OVR,$overlaps); $ovh= *OVR; }
die "bad -overlaps=$overlaps" unless($ok);

my($ngov, $nxov)= collect_overlaps($ovh); close($ovh);

my $inh= *STDIN;
$ok = ($input =~ /.gz$/) ? open($inh,"gunzip -c $input |") 
      : ($input =~ /^(stdin|-)/) ? $inh= *STDIN
      : open($inh,$input);
die "bad -input=$input" unless($ok);

# %sums=(); %sumid=(); 
my($ngin, $nxin, $ngsame) = filter_gff($inh);

warn"#overlaps over=$overlaps in=$input genes=$ngin same=$ngsame\n" if $debug;

# summary($input, $ngov, $nxov, $ngin, $nxin, $ngsame) if ($dosum);

#------------------------------------


my $didhead= 0;

sub testgene
{
  my($geneid, $generecIN, $geneother)= @_;
  
  my @generec= sort _sortgene @$generecIN;
  my $generec= \@generec;

  # my $mrna= $generec->[0]; #?? or check type?
  my($mrna)= grep{ $_->[2] eq "mRNA" } @generec;
  my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid)= @$mrna;
 
  return 0 if($skippatt and $tattr =~ m/$skippatt/);
  
  my $gidp=$geneid; #same as $gid;
  if($idequalpatt and $geneid =~ /($idequalpatt)/) { $gidp=$1; }
  my $ovgene= $overgenes{$gidp};
  my $issame= 0;
  
  unless($ovgene) {
    print join("\t",$geneid,0,0,0,"No_equiv"),"\n";
    
  } else {
    my($rref,$rsrc,$rtyp,$rb,$re,$rp,$ro,$rph,$rattr,$lid)= @{$ovgene->[0]};  
    # next unless($ref eq $rref and $rb < $te and $re > $tb); # no gene overlap

    my @HD=qw( EQgene EQcds Similar CXident 
               nca ncb cD cdsa cdsb nxa nxb xD xwa xwb);

    my @res= _similargene($generec, $ovgene);  
    my($issame1, $cdssame, $issimilar, $cxident,
      $nca, $ncb, $cd, $acdsw, $bcdsw, $nxa, $nxb, $xd, $axw, $bxw) = @res;
      # $nca, $acdsw, $nxa, $axw, $ncb, $bcdsw, $nxb, $bxw)= @res;
    print join("\t","Gene_ID",@HD),"\n" unless($didhead++);
    print join("\t",$geneid,@res),"\n";
    $issame= $issame1;
  }
  

  return ($issame)?1:0;
}

sub _samefeat 
{
  my($a,$b)= @_;
  return 0 unless( ref($a) and ref($b));
# ($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid)

  my($ab,$ae)=($a->[3],$a->[4]);
  my($bb,$be)=($b->[3],$b->[4]);

# if($MRNA_SPAN_WRONG) { # ncbi gff mRNA == gene span, not mRNA, for alt-tr; fix source
#   ($bb,$be)= ($ab,$ae)
#     if($a->[2] eq "mRNA" and $b->[2] eq "mRNA" and
#       ((abs($ab - $bb) <= $EXONSLOP) or (abs($ae - $be) <= $EXONSLOP)) );
# }
  if($EXONSLOP) {
  return ($a->[0] eq $b->[0]) # ref
      && ($a->[2] eq $b->[2]) # typ
      && (abs($ab - $bb) <= $EXONSLOP) # begin
      && (abs($ae - $be) <= $EXONSLOP) # end
      && _eqstrand($a->[6], $b->[6]) # ($a->[6] eq $b->[6]) # orient; allow "." to match + or - ?? ** FIXME yes
      ;
  } else {
  return ($a->[0] eq $b->[0]) # ref
      && ($a->[2] eq $b->[2]) # typ
      && ($ab == $bb) # begin
      && ($ae == $be) # end
      &&  _eqstrand($a->[6], $b->[6]) # ($a->[6] eq $b->[6]) # orient; allow "." to match + or - ??
      ;
  }
}

sub _samebases 
{
  my($a,$b)= @_;
  return (0,0) unless( ref($a) and ref($b));
# ($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid)
 
  my $wmax= 1 + _max( $b->[4] - $b->[3], $a->[4] - $a->[3]);
  if ( ($a->[0] eq $b->[0]) # ref
      && ($a->[2] eq $b->[2]) # typ
      &&  _eqstrand($a->[6], $b->[6]) # ($a->[6] eq $b->[6]) # orient; allow "." to match + or - ??
      && ($a->[3] < $b->[4]) # begin - end
      && ($a->[4] > $b->[3]) # end - beg
      )
    {
      my $omin= 1 + _min( $b->[4], $a->[4]) - _max( $b->[3], $a->[3]);
      # my $wmax= 1 + _max( $b->[4] - $b->[3], $a->[4] - $a->[3]);
      return( $omin, $wmax); #?
    }
  return (0,$wmax);
}


sub _similargene  
{
  my($agene,$bgene)= @_;  # agene=input; bgene=over
  my $na= @$agene;
  my $nb= @$bgene;
# ($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid)
  
  # option: classify _sameCDS, do both?
  my $gdiff= 0;
  my $nn= _max($na,$nb);
  my(@acds,@bcds, @aexon, @bexon);
  my($cdiff, $csame, $cba,$cbb, $cea,$ceb, $acdsw, $bcdsw,$lasti,$lastj)= (0) x 10;
  my($xdiff, $xsame, $axw, $bxw, $cdir,$xdir)= (0) x 10;
  
  for(my $i=0; $i<$nn; $i++) {
    my $issame= ($i>=$na or $i>=$nb)? 0 : _samefeat( $agene->[$i], $bgene->[$i]);
    $gdiff++ unless($issame);  
    
    if($i < $na and $agene->[$i]->[2] eq "CDS") { push(@acds,$i); }
    if($i < $nb and $bgene->[$i]->[2] eq "CDS") { push(@bcds,$i); }
    if($i < $na and $agene->[$i]->[2] eq "exon") { push(@aexon,$i); }
    if($i < $nb and $bgene->[$i]->[2] eq "exon") { push(@bexon,$i); }
  }

  my($nca,$ncb)= (scalar(@acds),scalar(@bcds));
  my($nxa,$nxb)= (scalar(@aexon),scalar(@bexon));
  
  # return (1,1,1,"I100") if($gdiff==0);
  return ( 1,1,1,"I100.100", 
    $nca, $ncb, $cdir, $acdsw, $bcdsw, $nxa, $nxb, $xdir, $axw, $bxw) if($gdiff==0);
  

# FIXME2: Probably bug, when overlap exon sizes differ greatly,
#         Aexon1tiny - Bexon1long = 1% over ; Aexon2long - Bexon2tiny = 1%over, final score is 1-2%over
#         but total base overlap is 99%, from exon split disagreement (e.g. in TE genes, but others too)
#          A:    --  ----------------------------------
#          B:   ---------------------------   ---     -      xsame/xmax is tiny here from exon shift
#
# >> redo exon overlap as max align transcript base overlap regardless of exon breaks
#    but keep cds overlap method same? diff introns mean diff protein

# FIXME: bug here, similar* : need to look inside acds,bcds for matching subset, not start from 0

  $nn= _max($nca,$ncb);
  ($lasti,$lastj)= (-1,-1);
  for ( my ($i,$j)=(0,0); $i < $nca or $j < $ncb; ) {  # USE i or j; need full acdsw,bcdsw
    my $ia= $acds[$i];  
    my $ib= $bcds[$j];
    my ($aok,$bok)= (defined $ia, defined $ib); # ? ($i<$nca, $b<$ncb)
    # my $aok= ($i < $nca)?1:0; my $bok= ($j < $ncb)?1:0;
    my ($sameb)= ($aok and $bok) ? _samebases( $agene->[$ia], $bgene->[$ib]) : 0;
    $csame += $sameb;

    if($aok) { 
      $acdsw += 1 + $agene->[$ia]->[4] - $agene->[$ia]->[3] if($i > $lasti);
      $cba= $agene->[$ia]->[3] if($i==0) ;
      $cea= $agene->[$ia]->[4] if($i==$nca-1);
      }
    if($bok) { 
      $bcdsw += 1 + $bgene->[$ib]->[4] - $bgene->[$ib]->[3] if($j > $lastj); 
      $cbb= $bgene->[$ib]->[3] if($j==0) ;
      $ceb= $bgene->[$ib]->[4] if($j==$ncb-1);
      }
    ($lasti,$lastj)= ($i,$j);
    
    unless($sameb) { 
      $cdiff++;
      if(($aok and $bok) and $agene->[$ia]->[4] < $bgene->[$ib]->[3]) { $i++; }
      elsif(($aok and $bok) and $bgene->[$ib]->[4] < $agene->[$ia]->[3]) { $j++; }
      else { $i++; $j++; } #??
    } else { 
      $i++; $j++; 
    }
    
  }

  my $cdsmax= _max( $acdsw, $bcdsw);
  my $cdsmin= _min( $acdsw, $bcdsw);
  my $cdsident= ($cdsmax>0) ? int( 0.5 + 100 * $csame / $cdsmax) : 0;
  # add: issimilar= cds_aencloses if $csame == $bcdsw;  cds_ainside if $csame == $acdsw
  
# FIXME: find align subset, same as cds
# FIXME2: >> redo exon overlap as max align transcript base overlap regardless of exon breaks; not same as cds

# FIXME ODD: 30jul wacky new results from same data, esp exon = 0 when CDS = 99;
  # die unless(checkgene($bgene, "_similargene.ov"));
  # ok here, still get wacky  C100.00  = no exon over when all CDS exon over

  $nn= _max($nxa,$nxb);
  ($lasti,$lastj)= (-1,-1);
  for ( my ($i,$j)=(0,0); $i < $nxa or $j < $nxb; ) { # USE i< OR j <
    my $ia= $aexon[$i];  
    my $ib= $bexon[$j];
    #?? my ($aok,$bok)= (defined $ia, defined $ib); # ? ($i<$nca, $b<$ncb)
    my $aok= ($i < $nxa)?1:0;
    my $bok= ($j < $nxb)?1:0;
    my $abok= ($aok and $bok)?1:0; ## this may have been bug: need 1:0 ??
    
    my ($sameb, $maxb)= ($abok) ? _samebases( $agene->[$ia], $bgene->[$ib]) : (0,0);
    $xsame += $sameb;
    
    my($ab,$ae,$bb,$be)= (0) x 10;
    if($aok) { 
      ($ab,$ae)= ($agene->[$ia]->[3], $agene->[$ia]->[4]);
      $axw += 1 + $ae - $ab if($i > $lasti); # WHOO, only if $i is NEW
      }
    if($bok) { 
      ($bb,$be)= ($bgene->[$ib]->[3], $bgene->[$ib]->[4]);
      $bxw += 1 + $be - $bb if($j > $lastj); 
      }
    ($lasti,$lastj)= ($i,$j);
            
    # new: 
    if($maxb > 0 and $sameb >= $maxb) {
      $i++; $j++;
      
    } elsif($sameb>0 ) { # and $sameb < $maxb .. look for more align to maxb
      if($abok and $ae < $be) { $i++; } # see if next aexon hits this bexon
      elsif($abok and $be < $ae) { $j++; } # other way      
      else { $i++; $j++; }

    } else { # if($sameb < 1)
      $xdiff++;
      if($abok and $ae < $bb) { $i++; }
      elsif($abok and $be < $ab) { $j++; }
      else { $i++; $j++; }
    }
      
  }

  my $xmax  = _max( $axw, $bxw);
  my $trident = ($xmax>0)   ? int( 0.5 + 100 * $xsame / $xmax) : 0;

  my $usame = _max( 0, $xsame - $csame);
  my $umax  = _max( $usame, _max( $axw - $acdsw, $bxw - $bcdsw));
  
  my $genesize= int( $cdsmax + $umax); # not same as xmax
  my $genesame= ($genesize>0) ? int( 0.5 + 100 * ($usame + $csame) / $genesize) : 0;
  my $issimilar= $genesame > 0 ? 1 : 0;

  # FIXME: need to format pad .trident w/ 0 for numeric compare
  if($trident < 10) { $trident = "0$trident"; } elsif($trident >= 100) { $trident="99"; } # do we ever get 100.100 ? should be all I100
  $genesame= ($cdsident>0 or $trident>0) ? "$cdsident.$trident" : 0;

  my $flag="";
  if($cdsident >= $SAME_CDS) { $flag = "C"; } ## $csame/$cdsmax ;   $cdiff == 0 and 
  elsif($MINID_CDS>0 and $cdsident < $MINID_CDS and ($MINID_UTR == 0 or $trident < $MINID_UTR))
  {
  # * use trident here not utrident, where exon overlap can be noncds x cds
  # check for large-cds over most of small-cds gene, but pct ID to cdsmax is small:
  $issimilar=0 if($cdsmin < 1 or ($MINID_CDS > 100 * $csame / $cdsmin) );
  }

  $cdir= ($acdsw>$bcdsw)?"+":($acdsw<$bcdsw)?"-":0;
  $xdir= ($axw>$bxw)?"+":($axw<$bxw)?"-":0;
  
  return ($gdiff?0:1, $cdiff?0:1, $issimilar, $flag.$genesame, 
    $nca, $ncb, $cdir, $acdsw, $bcdsw, $nxa, $nxb, $xdir, $axw, $bxw);

#   return ($gdiff?0:1, $cdiff?0:1, $issimilar, $flag.$genesame, 
#     $nca, $acdsw, $nxa, $axw, $ncb, $bcdsw, $nxb, $bxw);
 
}



# sub putgene
# {
#   my ($generec, $geneother)= @_;
#   ##  my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid]; 
#   foreach my $ft (@$generec, @$geneother) { 
#     if(ref $ft) { my @v= @$ft; print join("\t",@v[0..8]),"\n" if(@v>4); }
#     }
# }


sub _sortgene  
{
  # ($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid)
  ## return ($b->[2] cmp $a->[2]) # typ: reverse order by mRNA, exon, CDS : add something to ensure this?
  my($ta,$tb)= map{  (/gene/)?1:(/mRNA/)?2:(/CDS|exon/)?3:4; } ($a->[2],$b->[2]);
  return ($ta <=> $tb)
      || ($a->[0] cmp $b->[0]) # ref : should all be same for gene record
      || ($a->[3] <=> $b->[3]) # begin
      || ($a->[4] <=> $b->[4]) # end
      ;
}

# fixme for "." strand compare
sub _eqstrand { return($_[0] eq $_[1] || $_[0] eq '.' || $_[1] eq '.' )?1:0; }

# convert to array handling? _max(@list) 
sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }


sub filter_gff
{
  my($inh)= @_;
  my ($ng,$nx,$nr,$nsame,$nhit)= (0) x 10;
  my $nocomm= 1; ##($actid == ACT_NULL) ? 1 : 0;
  my @generec=();
  my @geneother=();
  my $geneid="";
  
  while(<$inh>) {
    unless(/^\w/){ print unless($nocomm or /^(#n |$)/); next; }
    my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr)= split"\t";
    $nr++; chomp($tattr);
   
    my($gid,$pid); 
    if($tattr =~ m/\bID=([^;]+)/) {  $gid=$1; }
    if($tattr =~ m/\bParent=([^;]+)/) {  $pid=$1; 
      $gid=$pid unless($gid and ($typ =~ /^($mrnatypes)$/)); 
      }
    unless(defined $gid) { $gid = "N".$ng; }

    my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid]; 

    if($typ =~ /^($mrnatypes)$/) {  
      # allow gene and mRNA types .. and/or keep other types in generec
      $nsame += testgene($geneid, \@generec, \@geneother) if(@generec);
      
      $ng++;
      @generec = ($rloc); # maybe best store as array of [gffcols], sorted by type-loc
      @geneother= ();
      $geneid= $gid;
      
    } elsif($typ =~ /^($exontypes)$/) {
      push @generec, $rloc; $nx++; # check Parent == $geneid ?
      
    } elsif( ! $passtypes or "$typ.$src" =~ m/$passtypes/ ) {
      push @geneother, $rloc; # check Parent == $geneid ?
    }
  }
  
  $nsame += testgene($geneid, \@generec, \@geneother) if(@generec);
  return ($ng,$nx,$nsame,$nhit);
}





sub addgene
{
  my($geneid, $generecIN)= @_;
  my @generec= sort _sortgene @$generecIN;
  my $generec= \@generec;

  if($overgenes{$geneid}) { # ERROR ?
    warn "# ERROR: duplicate overgene id=$geneid\n";
  }
  $overgenes{$geneid}= $generec;
  
#   my $mrna= $generec->[0]; #?? or check type?
#   my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid)= @$mrna;
#   my @bins= (int($tb/$BINSIZE) .. int($te/$BINSIZE));
#   foreach my $ib (@bins) { push( @{$overlaplist->{$ref}{$ib}}, $geneid); } # $generec
}

sub collect_overlaps
{
  my($gff)= @_;
  my ($ng,$nx)=(0,0);
  my @generec=(); my $geneid="";
  while(<$gff>){
    next unless(/^\w/); chomp;
    my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr)= split"\t";
    if($passtypes and "$typ.$src" !~ m/$passtypes/) { next; } # pass other types
    #^ drop passtypes for mrnatypes,exontypes
    $tattr ||="";
    
    my($gid,$pid); 
    if($tattr =~ m/\bID=([^;]+)/) {  $gid=$1; }
    if($tattr =~ m/\bParent=([^;]+)/) {  $pid=$1; 
      $gid=$pid unless($gid and ($typ =~ /^($mrnatypes)$/)); 
      }
    unless(defined $gid) { $gid = "N".$ng; }

    my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid]; 

    if($typ =~ /^($mrnatypes)$/) {  
      addgene($geneid, \@generec) if(@generec);
      $ng++; @generec = ($rloc); # maybe best store as array of [gffcols], sorted by type-loc
      $geneid= $gid;
      
    } elsif($typ =~ /^($exontypes)$/) {
      push @generec, $rloc;  $nx++;
    }
  }
  
  addgene($geneid, \@generec) if(@generec);
  
  warn"#collect_overlaps ngene=$ng, nexon=$nx\n" if $debug;
  return ($ng, $nx);
}

__END__

=item sample

>> evd3b024410t1 pasa change is alt-tr ok, but current best is valid longer prot (other alt intron used)

grep =evd3b024410t1 cacao11_bestgenes.pub3b.*gff | sed 's/;trg=.*//' 
scaffold_5      evd3b   mRNA    31383316        31387411        15979   -       .       
ID=evd3b024410t1;oid=AUGpier8ap13s_5g36t1;gene=evd3b024410;aalen=649,64%;
protein=
MKHKKLLLSLLEKCPNLKNLKKIHAHATTLGLLQNHNQALSCKILTTYANLNNPDDANRT
FNQIQRPDIVSWTCLIKLFVQYEDPFKSVLAFSQLIRNGLRPDTYSVVAALSACGKNKDL
DNGKLIHGVVSKYELGFENPIVGNALIDMYSRNGDILASELVFEWMFVKDIASWNSLLNG
FLLCNDLEASRRVFDKMPSRNAVSWTAMISGYVKGKEPLVGLKLFKEMKSEGKVDPTMVT
IVAVLAGCADSGGLYFGVSVHGYVKKVNLNEKNVVLSNALMDMYSKCGYLDVTAKIFNDM
VVRDVFSWTTMITGYAFHGKGKQALELFFDMLESRVAPNEVTFLSALSACSHEGLLVQGQ
RLFRIMVQRYGFKPKIEHYGCVLDLLGRAGLLEEAKMFIEEMPILPDAVVWRSLLCACLV
HGKLDLAEVAGKKVIELEPDDDGAYLLLWHMYSSTNRPEGAVKIRKLMRNQKVRKRPGCS
WLEVNGIVREFLAEYRPHYAGSDSCCILEGISEQSKLNEEFLWGRGGEEGFAKDFHKQTC
SLSVRSPYIFCFCSEELAQRIMRKTCPLPVYCLSISRRRTITGQFKAKANRWQRQLNFQA
IDCSKKKEKGSPAAQEEKADMASLGKPFFEHRISQGASLILLISARLLL;
homolog=646/1346,vitvi:GSVIVT01012938001;paralog=359/1346,
AUGpier8ap12s_5g169t1;pHOBEST=47%ho;inqual=66;nintron=4/6;
ovpro=80,m7vitvi:GSVIVT01012938001/80.00;ovrna=80,
r8cacao3vel4sc5Loc654t3/80.71,r8cacao3v3sc5Loc571t3_C1/7.12;
cxlen=1950/3040,64%;scoresum=15979;osrc=AUGpier8a;scorevec=
646,359,80,80,746,2624,4,66,0,42,0,734,1950
scaffold_5      evd3b   exon    31383316        31384429        1091    -       .       Parent=evd3b024410t1;est=329/1114;rseq=749/1114,cacao3v1Svsc5Loc644t2,cacao3vel4sc5Loc3385t1;intr=13,N99683
scaffold_5      evd3b   CDS     31384307        31384429        0       -       0       Parent=evd3b024410t1;
scaffold_5      evd3b   CDS     31385146        31385264        0       -       2       Parent=evd3b024410t1;
scaffold_5      evd3b   exon    31385146        31385264        189     -       .       Parent=evd3b024410t1;est=100/119;rseq=119/119,cacao3v1Svsc5Loc644t4;intr=-51/+21,N99683,N99684,N99685,N99686
scaffold_5      evd3b   CDS     31385541        31385680        0       -       1       Parent=evd3b024410t1;
scaffold_5      evd3b   exon    31385541        31385680        339     -       .       Parent=evd3b024410t1;est=140/140;rseq=140/140,cacao3v1Svsc5Loc644t3;intr=59,N99684,N99685,N99686
scaffold_5      evd3b   exon    31385745        31387411        1793    -       .       Parent=evd3b024410t1;est=177/1667;rseq=1616/1667,cacao3v1Svsc5Loc644t3
scaffold_5      evd3b   CDS     31385745        31387312        0       -       0       Parent=evd3b024410t1;


cacao11_bestgenes3b.update.ap3.gff # has several alttr
scaffold_5      puevd3b mRNA    31383316        31387411        .       -       .       
ID=puevd3b024410t1;pstatus=16;aalen=589;
protein=
MKHKKLLLSLLEKCPNLKNLKKIHAHATTLGLLQNHNQALSCKILTTYANLNNPDDANRT
FNQIQRPDIVSWTCLIKLFVQYEDPFKSVLAFSQLIRNGLRPDTYSVVAALSACGKNKDL
DNGKLIHGVVSKYELGFENPIVGNALIDMYSRNGDILASELVFEWMFVKDIASWNSLLNG
FLLCNDLEASRRVFDKMPSRNAVSWTAMISGYVKGKEPLVGLKLFKEMKSEGKVDPTMVT
IVAVLAGCADSGGLYFGVSVHGYVKKVNLNEKNVVLSNALMDMYSKCGYLDVTAKIFNDM
VVRDVFSWTTMITGYAFHGKGKQALELFFDMLESRVAPNEVTFLSALSACSHEGLLVQGQ
RLFRIMVQRYGFKPKIEHYGCVLDLLGRAGLLEEAKMFIEEMPILPDAVVWRSLLCACLV
HGKLDLAEVAGKKVIELEPDDDGAYLLLWHMYSSTNRPEGAVKIRKLMRNQKVRKRPGCS
WLEVNGIVREFLAEYRPHYAGSDSCCILEGISEQSKLNEEFLWGRGGEEGFAKDFHKQTC
SLSVRSPYIFCFCSEELAQRIMRKTCPLPAADEQSLVNSRLRLIDGRDN*;inqual=66
;nintron=4/6;ovpro=88,m7vitvi:GSVIVT01012938001/88.00;ovrna=
89,r8cacao3vel4sc5Loc654t3/89.72,r8cacao3v3sc5Loc571t3_C1/0.
12
scaffold_5      puevd3b exon    31385745        31387411 =      .       -       .       Parent=puevd3b024410t1;est=177/1667;rseq=1616/1667,cacao3v1Svsc5Loc644t3
scaffold_5      puevd3b CDS     31385745        31387312 =      .       -       0       Parent=puevd3b024410t1
scaffold_5      puevd3b exon    31385541        31385680 =     .       -       .       Parent=puevd3b024410t1;est=140/140;rseq=140/140,cacao3v1Svsc5Loc644t3;intr=59,N99684,N99685,N99686
scaffold_5      puevd3b CDS     31385541        31385680 =      .       -       1       Parent=puevd3b024410t1
scaffold_5      puevd3b exon    31385146        31385245 x      .       -       .       Parent=puevd3b024410t1;est=100/100;rseq=100/100,cacao3v1Svsc5Loc644t4;intr=-44/+20,N99683,N99684,N99685
scaffold_5      puevd3b CDS     31385184        31385245 x      .       -       2       Parent=puevd3b024410t1
scaffold_5      puevd3b exon    31383316        31384429 =      .       -       .       Parent=puevd3b024410t1;est=329/1114;rseq=749/1114,cacao3v1Svsc5Loc644t2,cacao3vel4sc5Loc3385t1;intr=13,N99683

//////////
cacao11_bestgenes3b.update.ap3.gff 
scaffold_5      puevd3b mRNA    30753500        30758143        .       +       .       ID=puevd3b024294t1;pstatus=16;aalen=729;
protein=MASMNNWLAFSLSPQELPSQTVDQDHHSQTAVSRLGFNSDDISGADVSGECFDLTSDSSAPSLNLPPPFGILEAFNRNNQSQDWNMKGLGMNSDGNYKTSSELSMLMGSSCNGQSLDQSNQ
EPKLENFLGNHSFSNHQQNKLHGCNTMYNTTTGEYMFPNCSLQLPSEDTTNARTSNGGDDNDNNNNKNNNNNTNINTGNGSSSIGLSMIKTWLRNQPAPPQPEAKNNGGASQSLSLSMSTGSQTGSPLP
LLTSSTGGGSGGESSSSDNNKQQKTPTGMDSESGAIEAMPRKSIDTFGQRTSIYRGVTRHRWTGRYEAHLWDNSCRREGQTRKGRQGGYDKEEKAARAYDLAALKYWGTTTTTNFPISNYEKELEEMKH
MTRQEYVASLRRKSSGFSRGASIYRGVTRHHQHGRWQARIGRVAGNKDLYLGTFSTQEEAAEAYDIAAIKFRGLNAVTNFDMSRYDVKSILESSTLPIGGAAKRLKDVEQAEMALDVQRVDDDNMSSQLT
DGINNYGAAHHGWPTIAFQQAQPFSMHYPYGQRVWCKQEQDSDANHTFQDLHQLQLGSTHNFFQPSVLHNLMAMDSSSMEHSSGSNSVIYCNGGGGDAAGSNGASGSYQAVGYGGNGGYVIPMGTVVASD
SNQNQGNGFGDNEVKTLGYETMYGSADPYHPRNLYYLSQQSSTGGVKASSYDQASACNNWVPTAVPTIAQRSSNMAVCHGAPTFTVWNDT*;
inqual=0;nintron=0/14;ovpro=79,m7ricco:29929.m004537/79.00;ovrna=41,r8B_g10091t00001/31.41,r8B_g16428t00001/21.21

scaffold_5      puevd3b exon    30753500        30753873        .       +       .       Parent=puevd3b024294t1
scaffold_5      puevd3b CDS     30753627        30753873        .       +       0       Parent=puevd3b024294t1
scaffold_5      puevd3b exon    30754068        30754746        .       +       .       Parent=puevd3b024294t1;
scaffold_5      puevd3b CDS     30754068        30754746        .       +       2       Parent=puevd3b024294t1
scaffold_5      puevd3b exon    30754869        30754951        .       +       .       Parent=puevd3b024294t1
scaffold_5      puevd3b CDS     30754869        30754951        .       +       1       Parent=puevd3b024294t1
scaffold_5      puevd3b exon    30755312        30755400        .       +       .       Parent=puevd3b024294t1
scaffold_5      puevd3b CDS     30755312        30755400        .       +       2       Parent=puevd3b024294t1
scaffold_5      puevd3b exon    30755819        30755892        .       +       .       Parent=puevd3b024294t1
scaffold_5      puevd3b CDS     30755819        30755892        .       +       0       Parent=puevd3b024294t1
scaffold_5      puevd3b exon    30755994        30756044        .       +       .       Parent=puevd3b024294t1
scaffold_5      puevd3b CDS     30755994        30756044        .       +       1       Parent=puevd3b024294t1
scaffold_5      puevd3b exon    30756367        30756443        .       +       .       Parent=puevd3b024294t1
scaffold_5      puevd3b CDS     30756367        30756443        .       +       1       Parent=puevd3b024294t1
scaffold_5      puevd3b exon    30756702        30758143        .       +       .       Parent=puevd3b024294t1;
scaffold_5      puevd3b CDS     30756702        30757591        .       +       2       Parent=puevd3b024294t1

cacao11_bestgenes.pub3b.gff 
scaffold_5      evd3b   mRNA    30751675        30758143        19815   +       .       ID=evd3b024294t1;oid=AUGepir1ap12s_5g224t1;aalen=729,68%;

scaffold_5      evd3b   exon    30751675        30752065        291     +       .       Parent=evd3b024294t1
scaffold_5      evd3b   exon    30753582        30753873        506     +       .       Parent=evd3b024294t1
                                ^^^^^^^^^^^^^^^^^^^^^^^^ == UTR change; ignore
scaffold_5      evd3b   CDS     30753627        30753873        0       +       0       Parent=evd3b024294t1
scaffold_5      evd3b   CDS     30754068        30754746        0       +       2       Parent=evd3b024294t1
scaffold_5      evd3b   exon    30754068        30754746        438     +       .       Parent=evd3b024294t1
scaffold_5      evd3b   CDS     30754869        30754951        0       +       1       Parent=evd3b024294t1
scaffold_5      evd3b   exon    30754869        30754951        0       +       .       Parent=evd3b024294t1
scaffold_5      evd3b   CDS     30755312        30755400        0       +       2       Parent=evd3b024294t1
scaffold_5      evd3b   exon    30755312        30755400        0       +       .       Parent=evd3b024294t1
scaffold_5      evd3b   CDS     30755819        30755892        0       +       0       Parent=evd3b024294t1
scaffold_5      evd3b   exon    30755819        30755892        0       +       .       Parent=evd3b024294t1
scaffold_5      evd3b   CDS     30755994        30756044        0       +       1       Parent=evd3b024294t1
scaffold_5      evd3b   exon    30755994        30756044        0       +       .       Parent=evd3b024294t1
scaffold_5      evd3b   CDS     30756367        30756443        0       +       1       Parent=evd3b024294t1
scaffold_5      evd3b   exon    30756367        30756443        0       +       .       Parent=evd3b024294t1
scaffold_5      evd3b   exon    30756702        30758143        1683    +       .       Parent=evd3b024294t1
scaffold_5      evd3b   CDS     30756702        30757591        0       +       2       Parent=evd3b024294t1;

=cut
