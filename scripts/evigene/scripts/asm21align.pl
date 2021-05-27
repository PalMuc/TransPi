#!/usr/bin/env perl
# asm21align.pl

=item about

  -in   nvit2_ncbiref.gff  : genes built on asm2
  -over nvit2_ncbiref-nvit1g8.gff.gz : transcripts aligned to asm1
  
  -out  nvit2_ncbiref.map2asm1.gff : genes moved to asm1 locs (problems marked)
  
  simple, identical offset of all exons means gene can be moved ok
    - includes shift to new scaffold and reversed orient
  
  match in,over gene by same id, not location overlap
  -id1 -id2 need flag which id tags to equate

  over (from gmap):
    idoverpatt= 'ID=ncbiref2:([\w\.]+)\.mrna\d'
    'ID=ncbiref2:([\w\.]+)\.mrna1;'  << only mrna1;
    
SCAFFOLD16      nasvit1asm      mRNA    5889    11002   .       -       .       
  ID=ncbiref2:XM_003425452.1.mrna1;Name=ncbiref2:XM_003425452.1;Parent=ncbiref2:XM_003425452.1.path1;Coverage=100.0;Identity=100.0
SCAFFOLD16      nasvit1asm      exon    10314   11002   100     -       .       
  ID=ncbiref2:XM_003425452.1.mrna1.exon1;Name=ncbiref2:XM_003425452.1;Parent=ncbiref2:XM_003425452.1.mrna1;Target=ncbiref2:XM_003425452.1 1 689 +
 
 in gene:
    idinpatt= 'transcript_id=([\w\.]+)'
SCAFFOLD16      RefSeq  mRNA    5889    11002   .       -       .       
  ID=NcbiRef_rna0;gid=gene0;gene=LOC100117425;product=XM_003425452.1;transcript_id=XM_003425452.1;Dbxref=tr:XM_003425452.1,pr:XP_003425500.1,GeneID:100117425,NASONIABASE:NV13865;Name=disintegrin and metalloproteinase domain-containing protein 10-like
SCAFFOLD16      RefSeq  exon    10314   11002   .       -       .       
  Parent=NcbiRef_rna0

=cut

use strict;
use warnings;
use Getopt::Long;

my $EXONSLOP = 0;
my $BINSIZE   = 900000 ; 
my $debug=1;
my $pctover= 0;

my ($overlaps,$dosum,$passtypes,$input,$itype,$typeover,$ok,$nin);
my $mrnatypes='mRNA|ncRNA|tRNA|transcript'; # ncRNA ??
my $exontypes='exon'; ## leave out CDS here?
my $overlaplist= {};
my %offsets= (); my $offsets;
# my @overgenes= (); my $overgeneid= 0;
my (%sums,%sumid);

my $idinpatt  ='ID=([^;\s]+)';
my $idoverpatt='ID=([^;\s]+)';


my $optok= GetOptions(
  "overlaps=s", \$overlaps, 
  "offsets=s", \$offsets, 
  "input=s", \$input,  
  "idinpatt=s", \$idinpatt, 
  "idoverpatt=s", \$idoverpatt, 
#   "typeover=s", \$typeover, 
#   "exontypes=s", \$exontypes, 
#   "passtypes=s", \$passtypes,  #??
#   "slopexon=i", \$EXONSLOP,  
#   "summary:s", \$dosum,
  "debug!", \$debug, 
  );

die "usage:
  asm12align.pl -in asm2gene.gff -idin $idinpatt -over asm1aligntr.gff -idover $idoverpatt > asm2gene_toasm1.gff

  ** ASSUMES gene model is mRNA > exon > CDS
  ** ASSUMES input gff is ordered by gene records (mRNA/exon/CDS all together per ID)
" unless($optok and $input and $overlaps);

my $ovh; 
   if($overlaps =~ /.gz$/) { $ok= open(OVR,"gunzip -c $overlaps |");  $ovh= *OVR; }
elsif($overlaps =~ /^(stdin|-)/) { $ovh= *STDIN; $ok=1; }
else { $ok= open(OVR,$overlaps); $ovh= *OVR; }
die "bad -overlaps=$overlaps" unless($ok);
my($ngov, $nxov)= collect_overlaps($ovh); close($ovh);

my $noff=0;
if($offsets) {
   if($offsets =~ /.gz$/) { $ok= open(OVR,"gunzip -c $offsets |");  $ovh= *OVR; }
elsif($offsets =~ /^(stdin|-)/) { $ovh= *STDIN; $ok=1; }
else { $ok= open(OVR,$offsets); $ovh= *OVR; }
die "bad -offsets=$offsets" unless($ok);
($noff)= collect_offsets($ovh); close($ovh);
}

my $inh= *STDIN;
$ok = ($input =~ /.gz$/) ? open($inh,"gunzip -c $input |") 
      : ($input =~ /^(stdin|-)/) ? $inh= *STDIN
      : open($inh,$input);
die "bad -input=$input" unless($ok);

%sums=(); %sumid=(); 
my($ngin, $nxin, $ngsame) = filter_gff($inh);

warn"#overlaps over=$overlaps nover=$ngov offsets=$offsets noff=$noff in=$input genes=$ngin same=$ngsame\n" if $debug;

#-----------


sub filter_gff
{
  my($inh)= @_;
  die "Err: need -idinpatt like 'ID=([^;]+)'" unless($idinpatt =~ m/[\(\)]/);  
  my ($ng,$nx,$nr,$nsame,$nhit)= (0) x 10;
  # my $printpass=1;
  # $printpass=0 if($actid == ACT_DROP or $actid == ACT_KEEP or $actid == ACT_NULL);
  my $nocomm= 0; # ($actid == ACT_NULL) ? 1 : 0;
  my @generec=();
  my @geneother=();
  my($lgid, $gid,$pid); 
  
  while(<$inh>) {
    unless(/^\w/){ print unless($nocomm or /^(#n |$)/); next; }
    my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr)= split"\t";
    $nr++; chomp($tattr);
    
    ##if($typ =~ /^($mrnatypes)$/ and $tattr =~ m/$idinpatt/) { $gid=$1; }
    if($typ =~ /^($mrnatypes)$/) { ($gid)= $tattr =~ m/$idinpatt/; }
    #not used# elsif($typ =~ /^($exontypes)$/ and $tattr =~ m/\bParent=([^;]+)/) {  $pid=$1; }

    my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid]; 

    if($typ =~ /^($mrnatypes)$/) {  
      # allow gene and mRNA types .. and/or keep other types in generec
      $nsame += testgene($lgid, \@generec, \@geneother) if(@generec);
      
      $ng++;
      @generec = ($rloc); # maybe best store as array of [gffcols], sorted by type-loc
      @geneother= ();
      ## warn "# Missing geneid: $_" unless($gid);
      $lgid= $gid;
      
    } elsif($typ =~ /^($exontypes)$/) {
      push @generec, $rloc; $nx++;       
      
    } elsif( ! $passtypes or "$typ.$src" =~ m/$passtypes/ ) {
      push @geneother, $rloc;
    }
  }
  
  $nsame += testgene($lgid, \@generec, \@geneother) if(@generec); # and lgid ??
  
  return ($ng,$nx,$nsame,$nhit);
}


sub testgene
{
  my($geneid, $generecIN, $geneother)= @_;
  
  my @generec= sort _sortgene @$generecIN;
  my $generec= \@generec;

  # my $mrna= $generec->[0]; #?? or grep mRNA
  # my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid)= @$mrna;
 
  $geneid ||= "nada";
  my $overgene= $overlaplist->{$geneid};
  if($overgene and ref $overgene) {
    my ($sameloc,$newgene,$newother,$matchtype,$oldloc)= 
      matchgene($generec, $overgene,$geneother) ;

    if($sameloc) { 
      putgene($newgene,$newother, "OK:$matchtype"); 
    } else {
      putgene($newgene,$newother, "ERROR:$matchtype");
    }
  
  } else {
    putgene($generec, $geneother, "ERROR:NOT_FOUND_overgene");
  }
}


sub putgene
{
  my ($generec, $geneother, $flags)= @_;
  ##  my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid]; 
  my $com=""; # ($flags =~ /NOT_FOUND/)?"#m." : ""; # ERROR?
  
  ##my $mrna= $generec->[0];
  my($mrna)= grep { $_->[2] =~ /^($mrnatypes)$/ } @$generec;
  if($mrna and $flags) { $mrna->[8]=~s/$/;aalign=$flags/; }
  
  foreach my $ft (@$generec, @$geneother) { 
    if(ref $ft) { my @v= @$ft; print $com, join("\t",@v[0..8]),"\n" if(@v>4); }
    }
}



sub matchgene
{
  my($generec,$ovgene,$geneother)= @_;   
  my $generecOLD= $generec;
  
  ## my $mrna= $generec->[0]; #?? or check type?
  my($mrna)= grep { $_->[2] =~ /^($mrnatypes)$/ } @$generec;
  unless($mrna) { $mrna= $generec->[0]; } # warn "# gene missing mrnatype: $gid\n";

  my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid)= @$mrna;
   
  my($rref,$rsrc,$rtyp,$rb,$re,$rp,$ro,$rph,$rattr,$lid)= @{$ovgene->[0]}; 
  
  my $origloc="$ref:$tb-$te:$to";
  my ($sameloc,  $error) = (0) x 10;
  if($ref eq $rref and $tb == $rb and $te == $re and _eqstrand($to,$ro) ) {
    $sameloc= _samegene($generec, $ovgene);
    ## option to return almostsameloc if here .. ie if mRNA sameloc, nochange
    return ($sameloc, $generecOLD, $geneother, "NO_CHANGE") if($sameloc);
  }
  
  my $shift="";
  my($nref,$offb,$offe,$no)= ($ref,0,0,$to);
  
  unless($ref eq $rref) { $nref= $rref; $shift .="r:$ref>$nref,"; }
  
  unless( _eqstrand($to,$ro) ) { 
    $no= $ro; $shift .="o:$to>$ro,";
    # need to reverse?  
    my $newgene = reverseft( $nref, $tb, $te, $no, $generec);
    $generec= $newgene;
    ##my $newother= reverseft( $nref, $tb, $te, $no, $geneother);
    ##$geneother= $newother;
  }
  
  #  # test all exons for same shift?
  $offb= $rb - $tb;
  $offe= $re - $te;
  if($offe == $offb) { } # ok or what?
  # my($shifterr,$offb,$offe,$errxi)= _exonshift($generec, $ovgene);
  
  # if($offe != $offb and $offb==0 or $offe==0) {} # probably not error but gmap slop
  
  my($roffb,$roffe)= shiftrange($rref,$rb,$re); # wrong one in offs file: rb,re not tb,te
  if( $roffb ) {
    $offb=$roffb; $offe=$roffe; # always?
    $shift .= "srange.";
  }
  
  if($offb or $offe) {
    $shift .="b:$offb,"; # $nb= $tb + $offb; 
    $shift .="e:$offe,"; # $ne= $te + $offe; 
  }

  # if(not $shift and not $sameloc) { } # recheck with _almostsameloc() ?
  return(1, $generecOLD,$geneother, "NOSHIFT:notsame") unless($shift); # error?  
  
  my $newgene = moveto( $nref, $offb, $offe, $no, $generec);
  $sameloc= _samegene($newgene, $ovgene);

  unless($sameloc or $shift =~ /srange/) {
    my($offb0,$offe0,$nerr,$nsame,$erri)= _exonshift($generec, $ovgene);
    if($nerr==0 or $nsame>$nerr) {
      $newgene = moveto( $nref, $offb0, $offe0, $no, $generec);
      $sameloc= _samegene($newgene, $ovgene);
      $shift .="xs$sameloc:$offb0,";
    }
  }

  my $newloc="";
  my $shiftok= ($sameloc or $shift =~ /r:/ or $shift =~ /srange/)?1:0;
  unless($shiftok) {
    my $newrna= $newgene->[0];  
    my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid)= @$newrna;
    $newloc="$ref:$tb-$te:$to";
  }
  
  # also write newgene but not sameloc if shifted r: o: ?
  if($shiftok) {
    $geneother = reverseft( $nref, $tb, $te, $no, $geneother) if($shift =~ /o:/);
    my $newother= moveto( $nref, $offb, $offe, $no, $geneother);
    return ($sameloc, $newgene, $newother, "SHIFTED:$shift;oldloc=$origloc");

  } elsif(($offb==0 or $offe==0) and $ref eq $rref and  _eqstrand($to,$ro)) {
    return(1, $generecOLD,$geneother, "ODDSHIFT:$shift;badnewloc=$newloc"); # error?  
  }
  
  return(0, $generecOLD, $geneother, "ALIGN_ERROR:$shift;badnewloc=$newloc"); # error? cannot match
}


sub shiftrange
{
  my( $ref, $tb, $te, $to)= @_;
  return unless($offsets{$ref});
  my @bins= (int($tb/$BINSIZE) .. int($te/$BINSIZE));
  foreach my $ib (@bins) {
    $offsets{$ref}{$ib} or next;
    my @rr= @{$offsets{$ref}{$ib}};
    foreach my $rr (@rr) {
      my($rb,$re,$roff)= @$rr;
      return ($roff,$roff) if($tb >= $rb and $te <= $re);
     }
  }
  return; 
}


sub moveto
{
  my( $nref, $offb, $offe, $no, $generec)= @_;
  my @agene= (); # @$generec;
  my $na= @$generec;
  for(my $i=0; $i<$na; $i++) {
    my @ft= @{$generec->[$i]}; # clone
    $ft[0]= $nref; $ft[6]= $no;
    $ft[3] += $offb; $ft[4] += $offe;
    push(@agene, \@ft);
  }
  return \@agene;
}

sub reverseft
{
  my( $nref, $starfof, $endof, $no, $generec)= @_;
  my @agene= (); 
  my $na= @$generec;
  for(my $i=$na-1; $i>=0; $i--) {
    my @ft= @{$generec->[$i]}; # clone
    $ft[0] = $nref; 
    $ft[6] = $no;
    $ft[3] = $starfof + $endof - $ft[4]; # $offb; 
    $ft[4] = $starfof + $endof - $ft[3]; # $offe;
    push(@agene, \@ft);
  }
  return \@agene;
}

sub _exonshift
{
  my($agene,$bgene)= @_;
  my $na= @$agene;
  my $nb= @$bgene; 
  # what if($na != $nb);
  my %offs;
  my($offb0,$offe0,$nsame,$nerr,$erri)= (0) x 10;
  for(my $i=0; $i<$na; $i++) {
    my $ax= $agene->[$i];
    my $bx= $bgene->[$i];
    next unless($ax->[2] =~ /^($exontypes)$/ and $bx);
    
    my($ab,$ae)=($ax->[3],$ax->[4]);
    my($bb,$be)=($bx->[3],$bx->[4]);
    my $offb= $bb - $ab;
    my $offe= $be - $ae;
    $offs{$offb}++; $offs{$offe}++; 
    if($i==0) { $offb0= $offb; $offe0= $offe; }
    elsif( $offb != $offb0 or $offe != $offe0) {
      $nerr++; $erri.="$i,";
      # $offb0= $offb; $offe0= $offe; #? reset
    } else {
      $nsame++;
    }
  }
  
  if($nerr>0 and $nsame > $nerr) { 
    my ($offmax)= sort{$offs{$b}<=>$offs{$a}} keys %offs;  #ok?? no?
    return($offmax,$offmax,$nerr,$nsame,$erri);
    }
  return($offb0,$offe0,$nerr,$nsame,$erri);
}


sub _samegene  
{
  my($agene,$bgene)= @_;
  my $na= @$agene;
  my $nb= @$bgene;
  return 0 unless($na == $nb); ## fixme: allow some end-exon slop
  for(my $i=0; $i<$na; $i++) {
    my $issame= _samefeat( $agene->[$i], $bgene->[$i]);
    return 0 unless $issame;
  }
  return 1;
}

sub _samefeat 
{
  my($a,$b)= @_;
  return 0 unless( ref($a) and ref($b));
# ($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid)

  my($ab,$ae)=($a->[3],$a->[4]);
  my($bb,$be)=($b->[3],$b->[4]);

  return ($a->[0] eq $b->[0]) # ref
      && ($a->[2] eq $b->[2]) # typ 
      && ($ab == $bb) # begin
      && ($ae == $be) # end
      &&  _eqstrand($a->[6], $b->[6]) # ($a->[6] eq $b->[6]) # orient; allow "." to match + or - ??
      ;
}

# convert to array handling? _max(@list) 
sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }

# fixme for "." strand compare
sub _eqstrand { return($_[0] eq $_[1] || $_[0] eq '.' || $_[1] eq '.' )?1:0; }



sub _sortgene  
{
  # ($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid)
  my($ta,$tb)= map{ (/gene/)?1:(/^($mrnatypes)$/)?2:(/^($exontypes)$/)?3:4; } ($a->[2],$b->[2]);
  return ($a->[0] cmp $b->[0]) # ref : should all be same for gene record
      || ($ta <=> $tb)
      || ($a->[3] <=> $b->[3]) # begin
      || ($a->[4] <=> $b->[4]) # end
      ;
}

sub addover
{
  my($geneid, $generecIN)= @_;
  my @generec= sort _sortgene @$generecIN;
  $overlaplist->{$geneid}= \@generec;
}

sub collect_offsets
{
  my($gff)= @_;
  my $noff=0;
  while(<$gff>){
    next unless(/^\w/); 
    my($ref,$tb,$te,$toff)= split;
    if($toff and $te>0) {
      my @bins= (int($tb/$BINSIZE) .. int($te/$BINSIZE));
      foreach my $ib (@bins) { push( @{$offsets{$ref}{$ib}}, [$tb,$te,$toff]); }  
      $noff++;
      }
  }
  return($noff); 
}

sub collect_overlaps
{
  my($gff)= @_;
  die "Err: need -idoverpatt like 'ID=([^;]+)'" unless($idoverpatt =~ m/[\(\)]/);  
  
  my ($ng,$nx, $lgid)=(0,0);
  my($gid,$pid); 
  my @generec=();
  while(<$gff>){
    next unless(/^\w/); 
    my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr)= split"\t";
    if($passtypes and "$typ.$src" !~ m/$passtypes/) { next; } # pass other types
    chomp($tattr); $tattr ||="";
    
    if($typ =~ /^($mrnatypes)$/) { ($gid)= $tattr =~ m/$idoverpatt/; }
    #not used# elsif($typ =~ /^($exontypes)$/ and $tattr =~ m/\bParent=([^;]+)/) {  $pid=$1; }
      ##no? $gid=$pid unless($gid and ($typ =~ /^($mrnatypes)$/)); 
    # unless(defined $gid) { $gid = "N".$ng; } # should be error

    my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid]; 

    if($typ =~ /^($mrnatypes)$/) {  
      addover($lgid, \@generec) if($lgid and @generec);
      $ng++; @generec = ($rloc); # maybe best store as array of [gffcols], sorted by type-loc
      ## warn "# Missing geneid: $_" unless($gid);
      $lgid= $gid;
    
    } elsif($lgid and $typ =~ /^($exontypes)$/) {
      push @generec, $rloc;  $nx++;
    }
  }
  
  addover($lgid, \@generec) if($lgid and @generec);
  
  warn "#collect_overlaps ngene=$ng, nexon=$nx\n" if $debug;
  return ($ng, $nx);
  # return \%overlaps;
}


=item problems

  some/many gmap problems .. how to decide?  if shift small, ignore?
  
    n=152 errors with o: strand change, most with scaf change, some histone/1-exon w/no strand or multiple identical loci
        -- ignore? only strand change if 1-exon
    n=1086 other errors, including trivial gmap mistakes
    
    n=96 are "unclassified transcription discrepancy"
    
 wc -l nvit2_ncbiref.remap1.*
     207 nvit2_ncbiref.remap1.err
    2199 nvit2_ncbiref.remap1.ok ; 1284 SHIFTED

** need to shift genes in shifted spans even if not SHIFTED (err,notsame,..)


SCAFFOLD7       432770  3885267 432039
SCAFFOLD1       5858735 9412973 5833309
SCAFFOLD6       1070217 4394652 135
SCAFFOLD8       953176  3674950 938744
SCAFFOLD20      783310  2033528 66
SCAFFOLD24      2091652 2447922 2049453
SCAFFOLD37      951325  1292287 365
SCAFFOLD5       4524461 4542011 392
SCAFFOLD98      407482  476323  110
SCAFFOLD41      711138  875719  92


 cat nvit2_ncbiref.remap1.ok | cut -f1,9 | sed 's/ID=.*;aalign=//; s/ODDSHIFT:.*/ODDSHIFT/; s/r:SCAFFOL
D.*/r:SCAF/;' | sort | uniq -c | sort -k1,1nr | less
 368 SCAFFOLD7  OK:SHIFTED:b:432039,e:432039,   span=432770-3885267
 307 SCAFFOLD1  OK:SHIFTED:b:5833309,e:5833309, span=5858735-9412973
 238 SCAFFOLD6  OK:SHIFTED:b:135,e:135,         span=1070217-4394652
 129 SCAFFOLD8  OK:SHIFTED:b:938744,e:938744,   span=953176-3674950
  84 SCAFFOLD20 OK:SHIFTED:b:66,e:66,           span=783310-2033528
  40 SCAFFOLD24 OK:SHIFTED:b:2049453,e:2049453, span=2091652-2447922
  29 SCAFFOLD2  OK:NOSHIFT:notsame
  26 SCAFFOLD37 OK:SHIFTED:b:365,e:365,         span=951325-1292287
   6 SCAFFOLD5  OK:SHIFTED:b:392,e:392,         span=4524461-4542011
   3 SCAFFOLD98 OK:SHIFTED:b:110,e:110,         span=407482-476323
   2 SCAFFOLD41 OK:SHIFTED:b:92,e:92,           span=711138-875719

   7 SCAFFOLD4793       OK:SHIFTED:r:SCAFFOLD4784>SCAFFOLD4793,
   6 SCAFFOLD4559       OK:SHIFTED:r:SCAFFOLD4550>SCAFFOLD4559,
   2 SCAFFOLD16 OK:SHIFTED:r:SCAFFOLD1610>SCAFFOLD16,b:3898950,e:3898950,
   2 SCAFFOLD5126       OK:SHIFTED:r:SCAFFOLD5115>SCAFFOLD5126,
   2 SCAFFOLD5285       OK:SHIFTED:r:SCAFFOLD5274>SCAFFOLD5285,
   2 SCAFFOLD5350       OK:SHIFTED:r:SCAFFOLD5339>SCAFFOLD5350,
   1 SCAFFOLD1028       OK:SHIFTED:r:SCAFFOLD1658>SCAFFOLD1028,b:14755,e:14781,
   1 SCAFFOLD7  OK:SHIFTED:b:441132,e:441132,

  many: SCAFFOLD6       ERROR:ALIGN_ERROR:b:135,e:135,
  many: SCAFFOLD7       ERROR:ALIGN_ERROR:b:432039,e:432039,

cut -f4 *badmap99 | sed 's/:.*//' | sort | uniq -c | sort -k1,1nr | less
 413 7sc
 352 1sc
 264 6sc
 141 8sc
  92 20sc
  50 24sc
  33 37sc
  22 5sc
  21 2sc
  18 4sc
  17 16sc
  16 10sc
  16 12sc
  16 18sc
  15 28sc
  13 9sc
  12 11sc
  11 3sc
  10 40sc

====    

  XM_001601976.2 ERROR:ALIGN_ERROR:b:90,e:0,
  XM_001599663.2 ERROR:ALIGN_ERROR:
  XM_003425420.1 ALIGN_ERROR:b:10408,e:0,

  XM_001608184.1 ERROR:ALIGN_ERROR:o:+>-,b:-1731028,e:-1731028  histone

  NM_001193319.1  ALIGN_ERROR:r:SCAFFOLD401>SCAFFOLD757,o:+>-,b:-145735,e:-143311 transcription discrepancy
  
  XM_001600739.2  ERROR:ALIGN_ERROR:r:SCAFFOLD5144>SCAFFOLD50,o:->+,b:770908,e:770908


  XM_003425500.1  ERROR:ALIGN_ERROR:o:->+,
  XM_001607187.2  ERROR:ALIGN_ERROR:r:SCAFFOLD16>SCAFFOLD5285,o:->+,b:-2751670,e:-2753147,
    ^^ has 2 mappings; this is mrna2, should have skipped
  XM_001608188.1  ERROR:ALIGN_ERROR:r:SCAFFOLD16>SCAFFOLD1610,o:->+,b:-3904474,e:-3904474,
  XM_003423751.1  ERROR:ALIGN_ERROR:b:2049453,e:2049453, << 1 exon P450 ; gmap alt map, 
  XM_001605083.2  OK:SHIFTED:b:2049453,e:2049453,
  
  XM_001607141.1  aalign=OK:SHIFTED:r:SCAFFOLD16>SCAFFOLD1,b:3883972,e:3883972
    SCAFFOLD1       RefSeq  mRNA    6050198 6050572 .       -

-- SCAFFOLD6 should all be simple +135 shift
SCAFFOLD6       RefSeq  mRNA    1051597 1057057 .       +       .       ID=NcbiRef_rna6192;gid=gene5471;gene=LOC10
0678678;product=XM_003424151.1;transcript_id=XM_003424151.1;Dbxref=tr:XM_003424151.1,pr:XP_003424199.1,GeneID:1006
78678;Name=tetratricopeptide repeat protein 17-like;aalign=ERROR:ALIGN_ERROR:b:135,e:135,
SCAFFOLD6       RefSeq  mRNA    1364249 1366474 .       -       .       ID=NcbiRef_rna6203;gid=gene5482;gene=LOC10
0121160;product=XM_001604696.2;transcript_id=XM_001604696.2;Dbxref=tr:XM_001604696.2,pr:XP_001604746.2,GeneID:1001
21160,NASONIABASE:NV11915;Name=hypothetical protein LOC100121160;aalign=ERROR:ALIGN_ERROR:b:135,e:135,
SCAFFOLD6       RefSeq  mRNA    1437517 2457583 .       +       .       ID=NcbiRef_rna6204;gid=gene5484;gene=LOC10
0677936;product=XM_003424174.1;transcript_id=XM_003424174.1;Dbxref=tr:XM_003424174.1,pr:XP_003424222.1,GeneID:1006
77936;Name=hypothetical protein LOC100677936;aalign=ERROR:ALIGN_ERROR:b:959650,e:135,
SCAFFOLD6       RefSeq  mRNA    1662934 1680532 .       +       .       ID=NcbiRef_rna6207;gid=gene5487;gene=LOC10
0117269;product=XM_001604324.2;transcript_id=XM_001604324.2;Dbxref=tr:XM_001604324.2,pr:XP_001604374.2,GeneID:1001
17269,NASONIABASE:NV11918;Name=zinc transporter ZIP11-like;aalign=ERROR:ALIGN_ERROR:b:135,e:135,
SCAFFOLD6       RefSeq  mRNA    1728205 1732891 .       +       .       ID=NcbiRef_rna6211;gid=gene5491;gene=LOC10
0121326;product=XM_001604882.1;transcript_id=XM_001604882.1;Dbxref=tr:XM_001604882.1,pr:XP_001604932.1,GeneID:1001
21326,NASONIABASE:NV11921;Name=transcription factor Sox-7-like;aalign=ERROR:ALIGN_ERROR:b:222,e:135,
SCAFFOLD6       RefSeq  mRNA    1736813 1743085 .       -       .       ID=NcbiRef_rna6212;gid=gene5493;gene=LOC10
0121363;product=XM_003424161.1;transcript_id=XM_003424161.1;Dbxref=tr:XM_003424161.1,pr:XP_003424209.1,GeneID:1001
21363,NASONIABASE:NV11923;Name=85 kDa calcium-independent phospholipase A2-like;aalign=ERROR:ALIGN_ERROR:b:135,e:1
35,
SCAFFOLD6       RefSeq  mRNA    1783373 1784759 .       +       .       ID=NcbiRef_rna6224;gid=gene5504;gene=LOC10
0121567;product=XM_001605126.1;transcript_id=XM_001605126.1;Dbxref=tr:XM_001605126.1,pr:XP_001605176.1,GeneID:1001
21567,NASONIABASE:NV11935;Name=calcium and integrin-binding family member 3-like;aalign=ERROR:ALIGN_ERROR:b:135,e:
135,
SCAFFOLD6       RefSeq  mRNA    1792629 1804349 .       +       .       ID=NcbiRef_rna6227;gid=gene5507;gene=LOC10
0117500;product=XM_001605234.2;transcript_id=XM_001605234.2;Dbxref=tr:XM_001605234.2,pr:XP_001605284.1,GeneID:1001
17500,NASONIABASE:NV11939;Name=60S ribosomal protein L18a-like;aalign=ERROR:ALIGN_ERROR:b:10627,e:135,
SCAFFOLD6       RefSeq  mRNA    2052707 2054514 .       -       .       ID=NcbiRef_rna6264;gid=gene5536;gene=LOC10
0122030;product=XM_003424237.1;transcript_id=XM_003424237.1;Dbxref=tr:XM_003424237.1,pr:XP_003424285.1,GeneID:1001
22030,NASONIABASE:NV50086;Name=probable multidrug resistance-associated protein lethal(2)03659-like;aalign=ERROR:A
LIGN_ERROR:b:135,e:49,
SCAFFOLD6       RefSeq  mRNA    2059668 2061162 .       +       .       ID=NcbiRef_rna6265;gid=gene5537;gene=LOC10
0680417;product=XM_003424170.1;transcript_id=XM_003424170.1;Dbxref=tr:XM_003424170.1,pr:XP_003424218.1,GeneID:1006
80417;Name=hypothetical protein LOC100680417;aalign=ERROR:ALIGN_ERROR:b:135,e:135,
SCAFFOLD6       RefSeq  mRNA    2566945 2603676 .       -       .       ID=NcbiRef_rna6327;gid=gene5587;gene=LOC10
0122742;product=XM_001606290.2;transcript_id=XM_001606290.2;Dbxref=tr:XM_001606290.2,pr:XP_001606340.2,GeneID:1001
22742,NASONIABASE:NV12007;Name=transcriptional repressor p66-beta-like;aalign=ERROR:ALIGN_ERROR:b:31551,e:135,
SCAFFOLD6       RefSeq  mRNA    2636519 2656928 .       +       .       ID=NcbiRef_rna6336;gid=gene5594;gene=LOC10
0122788;product=XM_001606341.2;transcript_id=XM_001606341.2;Dbxref=tr:XM_001606341.2,pr:XP_001606391.2,GeneID:1001
22788,NASONIABASE:NV12011;Name=spectrin beta chain, brain 1-like;aalign=ERROR:ALIGN_ERROR:b:3129,e:135,
SCAFFOLD6       RefSeq  mRNA    2988257 2990718 .       -       .       ID=NcbiRef_rna6361;gid=gene5620;gene=LOC10
0123080;product=XM_001606636.1;transcript_id=XM_001606636.1;Dbxref=tr:XM_001606636.1,pr:XP_001606686.1,GeneID:1001
23080,NASONIABASE:NV12024;Name=hypothetical protein LOC100123080;aalign=ERROR:ALIGN_ERROR:b:135,e:135,
SCAFFOLD6       RefSeq  mRNA    3221405 3228822 .       -       .       ID=NcbiRef_rna6378;gid=gene5634;gene=LOC10
0678887;product=XM_003424248.1;transcript_id=XM_003424248.1;Dbxref=tr:XM_003424248.1,pr:XP_003424296.1,GeneID:1006
78887;Name=hypothetical protein LOC100678887;aalign=ERROR:ALIGN_ERROR:b:135,e:135,
SCAFFOLD6       RefSeq  mRNA    3310979 3312662 .       +       .       ID=NcbiRef_rna6386;gid=gene5642;gene=Or69;
product=NM_001190580.1;transcript_id=NM_001190580.1;Dbxref=tr:NM_001190580.1,pr:NP_001177509.1,GeneID:100463076;Na
me=odorant receptor 69;aalign=ERROR:ALIGN_ERROR:b:135,e:135,
SCAFFOLD6       RefSeq  mRNA    3811362 3812479 .       +       .       ID=NcbiRef_rna6414;gid=gene5666;gene=LOC10
0123409;product=XM_001607712.1;transcript_id=XM_001607712.1;Dbxref=tr:XM_001607712.1,pr:XP_001607762.1,GeneID:1001
23409,NASONIABASE:NV12046;Name=hypothetical protein LOC100123409;aalign=ERROR:ALIGN_ERROR:b:135,e:135,
SCAFFOLD6       RefSeq  mRNA    3980229 3982306 .       -       .       ID=NcbiRef_rna6420;gid=gene5671;gene=LOC10
0123463;product=XM_003424253.1;transcript_id=XM_003424253.1;Dbxref=tr:XM_003424253.1,pr:XP_003424301.1,GeneID:1001
23463,NASONIABASE:NV50089;Name=sodium channel protein Nach-like;aalign=ERROR:ALIGN_ERROR:b:135,e:135,
SCAFFOLD6       RefSeq  mRNA    4018367 4026044 .       -       .       ID=NcbiRef_rna6424;gid=gene5675;gene=LOC10
0123500;product=XM_001607094.2;transcript_id=XM_001607094.2;Dbxref=tr:XM_001607094.2,pr:XP_001607144.2,GeneID:1001
23500,NASONIABASE:NV12054;Name=hypothetical protein LOC100123500;aalign=ERROR:ALIGN_ERROR:b:213,e:135,
SCAFFOLD6       RefSeq  mRNA    4058086 4081733 .       +       .       ID=NcbiRef_rna6435;gid=gene5685;gene=LOC10
0118263;product=XM_001607808.2;transcript_id=XM_001607808.2;Dbxref=tr:XM_001607808.2,pr:XP_001607858.1,GeneID:1001


-- SCAFFOLD24 all +2049453  shift
SCAFFOLD24      RefSeq  mRNA    57345   73512   .       +       .       ID=NcbiRef_rna807;gid=gene697;gene=LOC1001
21412;product=XM_001604975.2;transcript_id=XM_001604975.2;Dbxref=tr:XM_001604975.2,pr:XP_001605025.1,GeneID:100121
412,NASONIABASE:NV15329;Name=duodenase-1-like;aalign=ERROR:ALIGN_ERROR:b:2049453,e:2034575,
SCAFFOLD24      RefSeq  mRNA    80824   82968   .       +       .       ID=NcbiRef_rna813;gid=gene704;gene=CYP9AH5
;product=XM_003423751.1;transcript_id=XM_003423751.1;Dbxref=tr:XM_003423751.1,pr:XP_003423799.1,GeneID:100121504,N
ASONIABASE:NV15333;Name=cytochrome P450 9AH5;aalign=ERROR:ALIGN_ERROR:b:2049453,e:2049453,
SCAFFOLD24      RefSeq  mRNA    107633  108861  .       -       .       ID=NcbiRef_rna819;gid=gene709;gene=SPH103;
product=NM_001172785.1;transcript_id=NM_001172785.1;Dbxref=tr:NM_001172785.1,pr:NP_001166256.1,GeneID:100121577,NA
SONIABASE:NV15336;Name=serine protease homolog 103;aalign=ERROR:ALIGN_ERROR:b:2049453,e:2049453,
SCAFFOLD24      RefSeq  mRNA    324665  327593  .       +       .       ID=NcbiRef_rna844;gid=gene731;gene=LOC1001
15951;product=XM_001606735.2;transcript_id=XM_001606735.2;Dbxref=tr:XM_001606735.2,pr:XP_001606785.1,GeneID:100115
951,NASONIABASE:NV15352;Name=ubiquitin-conjugating enzyme E2 N-like;aalign=ERROR:ALIGN_ERROR:b:2049453,e:2049453,

-- SCAFFOLD20 all +66 shift
SCAFFOLD20      RefSeq  mRNA    1191751 1209143 .       -       .       ID=NcbiRef_rna6972;gid=gene6169;gene=LOC10
0118137;product=XM_001602132.2;transcript_id=XM_001602132.2;Dbxref=tr:XM_001602132.2,pr:XP_001602182.2,GeneID:1001
18137,NASONIABASE:NV14708;Name=heterogeneous nuclear ribonucleoprotein Q-like;aalign=ERROR:ALIGN_ERROR:b:66,e:66,
SCAFFOLD20      RefSeq  mRNA    1688563 1689531 .       -       .       ID=NcbiRef_rna7023;gid=gene6218;gene=LOC10
0115908;product=NM_001161678.1;transcript_id=NM_001161678.1;Dbxref=tr:NM_001161678.1,pr:NP_001155150.1,GeneID:1001
15908,NASONIABASE:NV14748;Name=GOBP-like venom protein;aalign=ERROR:ALIGN_ERROR:b:66,e:66,
...

-- SCAFFOLD7 all +432039 shift
SCAFFOLD7       RefSeq  mRNA    43849   48782   .       +       .       ID=NcbiRef_rna10709;gid=gene9477;gene=LOC1
00117239;product=XM_001601489.2;transcript_id=XM_001601489.2;Dbxref=tr:XM_001601489.2,pr:XP_001601539.2,GeneID:100
117239,NASONIABASE:NV12124;Name=hypothetical protein LOC100117239;aalign=ERROR:ALIGN_ERROR:b:432039,e:432039,
SCAFFOLD7       RefSeq  mRNA    61353   63533   .       +       .       ID=NcbiRef_rna10713;gid=gene9480;gene=LOC1
00117272;product=XM_001601514.2;transcript_id=XM_001601514.2;Dbxref=tr:XM_001601514.2,pr:XP_001601564.1,GeneID:100
117272,NASONIABASE:NV12128;Name=CD151 antigen-like;aalign=ERROR:ALIGN_ERROR:b:432039,e:432039,
SCAFFOLD7       RefSeq  mRNA    68840   226907  .       +       .       ID=NcbiRef_rna10716;gid=gene9481;gene=LOC1
00678286;product=XM_003424397.1;transcript_id=XM_003424397.1;Dbxref=tr:XM_003424397.1,pr:XP_003424445.1,GeneID:100
678286;Name=hypothetical protein LOC100678286;aalign=ERROR:ALIGN_ERROR:b:432039,e:432039,
...

-- Sc1 all +5833309 shift
SCAFFOLD1       RefSeq  mRNA    540277  547455  .       -       .       ID=NcbiRef_rna11149;gid=gene9843;gene=LOC1
00124010;product=XM_001607770.2;transcript_id=XM_001607770.2;Dbxref=tr:XM_001607770.2,pr:XP_001607820.2,GeneID:100
124010,NASONIABASE:NV10354;Name=hypothetical protein LOC100124010;aalign=ERROR:ALIGN_ERROR:b:5833309,e:5833309,
SCAFFOLD1       RefSeq  mRNA    594448  597391  .       +       .       ID=NcbiRef_rna11158;gid=gene9852;gene=LOC1
00115111;product=XM_001608191.2;transcript_id=XM_001608191.2;Dbxref=tr:XM_001608191.2,pr:XP_001608241.1,GeneID:100
115111,NASONIABASE:NV10361;Name=26S proteasome non-ATPase regulatory subunit 12-like;aalign=ERROR:ALIGN_ERROR:b:58
33309,e:5833309,


-- these are same, but got ERROR:ALIGN_ERROR:
SCAFFOLD16      RefSeq  mRNA    163739  167490  .       -       .       ID=NcbiRef_rna30;gid=gene20;gene=LOC100118350;product=XM_001602288.2;transcript_id=XM_001602288.2;Dbxref=tr:XM_001602288.2,pr:XP_001602338.2,GeneID:100118350,NASONIABASE:NV13885;Name=hypothetical protein LOC100118350
SCAFFOLD16      RefSeq  exon    165936  167490  .       -       .       Parent=NcbiRef_rna30
SCAFFOLD16      RefSeq  exon    165041  165862  .       -       .       Parent=NcbiRef_rna30
SCAFFOLD16      RefSeq  exon    163739  164848  .       -       .       Parent=NcbiRef_rna30


SCAFFOLD16      nasvit1asm      mRNA    163739  167490  .       -       .       ID=ncbiref2:XM_001602288.2.mrna1;Name=ncbiref2:XM_001602288.2;Parent=ncbiref2:XM_001602288.2.path1;Coverage=100.0;Identity=100.0
SCAFFOLD16      nasvit1asm      exon    165936  167490  100     -       .       ID=ncbiref2:XM_001602288.2.mrna1.exon1;Name=ncbiref2:XM_001602288.2;Parent=ncbiref2:XM_001602288.2.mrna1;Target=ncbiref2:XM_001602288.2 1 1555 +
SCAFFOLD16      nasvit1asm      exon    165041  165862  100     -       .       ID=ncbiref2:XM_001602288.2.mrna1.exon2;Name=ncbiref2:XM_001602288.2;Parent=ncbiref2:XM_001602288.2.mrna1;Target=ncbiref2:XM_001602288.2 1556 2377 +
SCAFFOLD16      nasvit1asm      exon    163739  164848  100     -       .       ID=ncbiref2:XM_001602288.2.mrna1.exon3;Name=ncbiref2:XM_001602288.2;Parent=ncbiref2:XM_001602288.2.mrna1;Target=ncbiref2:XM_001602288.2 2378 3487 +

    
=cut

