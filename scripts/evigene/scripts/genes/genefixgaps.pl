#!/usr/bin/env perl
# genefixgaps.pl from geneattr.pl

=item about

  genefixgaps.pl -in pub3i.gapgenes.fix.gff >
  tedious hack handler for gene.gapfix.gff, 
  exon/CDS w/ asm gaps have been marked, with gapfix= action to handle
  -- most gapfix in utr-exon are chop of end loc, mRNA and gene locs need adjust to match
  -- a few gapfix=dropexon or gapfix=cds-split-gene
  
=item gapfix annots

    gapfix=b+29 ; gapfix=e-11; 
        EXON_LOC_IS_ADJUSTED ** input exon end has been adjusted: -exonfix , or -noexonfix
    gapfix=reclass:Poor == drop this record ..
      gapfix=reclass:Poor-partof == drop also
    gapfix=keep-gene-as-cds-split
    gapfix=dropexon  == utr exon to drop

** FIXME:    
# gap in utr of internal exon, drop 5'utr exon and trim two utr gaps; trim at 33293441  ; LLR disease gene must keep
scaffold_5	evg3	exon	33293441	33294661	0	+	.	Parent=Thecc1EG025537t1;;gap=100/1416;gapfix=btrimto-33293441
    
    
=item author
  
  don gilbert, gilbertd near indiana edu, ca 2011
  part of EvidentialGene, evigene/scripts/
  
=cut


use strict;
use warnings;
use Getopt::Long;

my $EXON_LOC_IS_ADJUSTED= 1; #?? off by default or on?
#my $TINYCDS = 40; # ignore tiny CDS tests for typeover == CDS
# my $EXONSLOP = 0;
# my $SAME_CDS= 90; # ** option
# my $MINID_CDS= 0;
# my $MINID_UTR= 0;

my $mrnatypes='mRNA';
my $exontypes='exon|CDS';
my $passtypes="";

my %overgenes= (); ## my $overgeneid= 0;
my($input,$didhead,$overlaps,$debug,$aadef,$idlist,$ok,$offby1list,$generowGlobal)= (0) x 20;

my $optok= GetOptions(
  # "overlaps=s", \$overlaps, 
  "input=s", \$input,  
  # "aadef=s", \$aadef,  
  "offby1list=s", \$offby1list,  
  "idlist=s", \$idlist,  
  "exonfixed!", \$EXON_LOC_IS_ADJUSTED, 
  "debug!", \$debug, 
  );

die "usage: genefixgaps.pl < genes.gff  > genes.fixed.gff
  ** ASSUMES gene model is mRNA > exon > CDS
  ** ASSUMES input gff is ordered by gene records (mRNA/exon/CDS all together per ID)
" unless($optok); ##  and $input and $overlaps

my %aadef=();
my %idlist=();
if($idlist) {
  open(F,$idlist) or die "missing idlist=$idlist";
  while(<F>) { my($id)=split; $id=~s/>//; $idlist{$id}=1; } close(F);
}

# NCBI tbl2asn 2012.Jun.14,v200, has offby1 bug for gap detection; Gap at gene begin-1 is called gene-gap error.
# need offby end indicator even though all these are begin, use base value?  +1 +10 = begin+n, -1 -10 = end-n ?
my %offby1list;
if($offby1list) { 
  open(F,$offby1list) or die "missing offby1list=$offby1list"; 
  while(<F>) { my($id,$offend)=split; $id=~s/>//; $offby1list{$id}=$offend; } close(F);
}


my $inh= *STDIN;
if($input) {
$ok = ($input =~ /.gz$/) ? open($inh,"gunzip -c $input |") 
      : ($input =~ /^(stdin|-)/) ? $inh= *STDIN
      : open($inh,$input);
die "bad -input=$input" unless($ok);
}

# my($ngin, $nxin, $ngsame) = 
  filter_gff($inh);

# if(0) {
# foreach my $gid (grep{ $idlist{$_} > 0 } sort keys %idlist)
# {
#   # (my $gid2=$gid) =~ s/_[GC]\d+$//;
#   my($src, $nintron, $aalen, $cxlen, $location)=("na",0,"na","na","unmapped"); # unmapped?
#   ($aalen,$cxlen)= aadef( $gid, $aalen, $cxlen);
#   print join("\t", $gid, $src, $aalen, $cxlen, $nintron, $location)."\n";
#   # $idlist{$gid}= -9; # done flag
# }
# }
# warn"#overlaps over=$overlaps in=$input genes=$ngin same=$ngsame\n" if $debug;


#------------------------------------


sub putgene {
  my($generec, $geneother) = @_;
  $geneother ||= [];  my $nput=0;
  # note generec contains gene, mrna, exons, cds, in orig order; updated locs ..
  # my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid]; 
  foreach my $rloc (@$generec, @$geneother) {
    next unless(ref $rloc); # gene row may be undef
    # next if($rloc->[0] =~ /^#drop./); # is this ok?
    print join("\t", @$rloc[0..8])."\n"; $nput++;
  }
  return $nput; # or nput?
}

sub testgene
{
  my($geneid, $generecIN, $geneother)= @_;
  
  my ($dropgene,$nfix,)= (0) x 10;
  my $fixann="";
  
  ## problem? @generec[0] may be undef gene slot
  
  my @generec=  @$generecIN; # need to sort or not?
  #? my @generec= sort _sortgene @$generecIN; #? sort by genostart or 5'start?

  my($generow)= grep{ $_->[2] eq "gene" } @generec;  
  my($mrna)= grep{ $_->[2] eq "mRNA" } @generec; # must have
  my @cds  = grep{ $_->[2] eq "CDS"  } @generec;
  my @exons= grep{ $_->[2] eq "exon" } @generec;
  my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid)= @$mrna;
  #return 0 if($skippatt and $tattr =~ m/$skippatt/);

  my $issplit= ($mrna->[8] =~ /splitgene|Split=/)?1:0;
  my ($genegap)= ($generow) ? $generow->[8] =~ /(gap=[^;\s]+)/ : ""; # dont always have this annot.. opt
  # ^^ test this vs nfix
  
  (my $gid2=$gid) =~ s/_[GC]\d+$//;
  $idlist{$gid}= $idlist{$gid2}= -9; # done flag

  ## my($gmb,$gme,$gmb0,$gme0,$gmoff,$gmid)=(0) x 10;
  my ($gmid,$gmb0,$gme0,)= (0) x 10;
  if($generow) {
    ($gmid)= $generow->[8] =~ /ID=([^;\s]+)/;
    ($gmb0,$gme0)= ($generow->[3], $generow->[4]);
  } elsif($generowGlobal) { # for alt mrna .. hack fix
    ($gmid)= $generowGlobal->[8] =~ /ID=([^;\s]+)/;
    ($gmb0,$gme0)= ($generowGlobal->[3], $generowGlobal->[4]);
  } else {
    ($gmid = $gid2) =~ s/t\d+$//;
  }
  my $gmoff= $offby1list{$gmid} || 0;
  
#     gapfix=reclass:Poor == drop this record ..
#       gapfix=reclass:Poor-partof == drop also
#     gapfix=keep-gene-as-cds-split
#     gapfix=dropexon  == utr exon to drop
#     gapfix=btrimto-33293441  << is this ok? almost same as b+n, e-n
#     gapfix=b+29 ; gapfix=e-11; 
## upd 20150420: add annot gapfill=gapb-gape .. on exons from trasm to genome, indicates tr exon fills genome gap
##     .. dont change xons, but add annots to mrna
## eg  KN805956.1	kf2rae5	exon	195439	195481	87	+	.	Parent=Funhe2EKm009615t1;
#   ggap=2/43,gapKN805956_195464,ovspan:195464-195465;gapfill=195464-195465
  my @gapfill=();
  
  my ($cmin,$cmax,$clen,$cdiff,$cdsb,$cdse)= (-1,0,0,0,0,0);
  foreach my $ex (@cds) {
    my($b,$e,$xat)= @{$ex}[3,4,8]; 
    my($b0,$e0)= ($b,$e);
    my $xd=0; 
    my($fx)= $xat =~ m/gapfix=([^;\s]+)/;
    if($fx) { 
      $nfix++; my $xto=0;
      if($fx =~ /^drop/) { $ex->[0] =~ s/^/#drop./; $fixann.="dropcds,"; next; } # is this ok?
      elsif($fx =~ /^reclass.Poor/) { $dropgene=1; $fixann.="dropgene,"; }
      elsif($fx =~ /^b([+-]\d+)/) { $xd=$1; $b += $xd unless($EXON_LOC_IS_ADJUSTED); $fixann.="$fx,";}
      elsif($fx =~ /^e([+-]\d+)/) { $xd=$1; $e += $xd unless($EXON_LOC_IS_ADJUSTED); $fixann.="$fx,";}
      elsif($fx =~ /^btrimto.(\d+)/) { $xto=$1; $xd=$xto-$b; $b = $xto unless($EXON_LOC_IS_ADJUSTED); $fixann.="$fx,";}
      elsif($fx =~ /^etrimto.(\d+)/) { $xto=$1; $xd=$xto-$e; $e = $xto unless($EXON_LOC_IS_ADJUSTED); $fixann.="$fx,";}
      elsif($fx =~ /^keep-gene-as-cds-split/ and $issplit) { $fixann.="splitgene,"; }
      else { $fixann.="fail:$fx,"; }
      }

    # *** FIXME: not EXON_LOC_IS_ADJUSTED or gmoff by1 change to exon b,e needs push back to $ex->[3,4] ***
    $ex->[3]= $b if($b0 != $b); $ex->[4]= $e if($e0 != $e);
    
    my $w= 1+$e-$b; $clen += $w;  $cdiff += $xd;
    foreach my $v ($b,$e) { 
      $cmax=$v if($v>$cmax); 
      $cmin=$v if($cmin<0 || $v<$cmin); 
      } 
  }

  my ($xmin,$xmax,$xlen,$xdiff)= (-1,0,0,0,0,0);
  my $nexon= @exons; my $iexon=-1;
  foreach my $ex (@exons) {  $iexon++;
    # my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid]; 
    my($xr,$xs,$xt,$b,$e,$p,$o,$tph,$xat)= @$ex;
    my($b0,$e0)= ($b,$e);
    my $xd=0; 
    my($fx)= $xat =~ m/gapfix=([^;\s]+)/;
    if($fx) { 
      $nfix++; my $xto=0;
      if($fx =~ /^drop/) { $ex->[0] =~ s/^/#drop./; $fixann.="dropexon,";next; } # is this ok?
      elsif($fx =~ /^reclass.Poor/) { $dropgene=1; $fixann.="dropgene,"; }
      elsif($fx =~ /^b([+-]\d+)/) { $xd=$1; $b += $xd unless($EXON_LOC_IS_ADJUSTED); $fixann.="$fx,";}
      elsif($fx =~ /^e([+-]\d+)/) { $xd=$1; $e += $xd unless($EXON_LOC_IS_ADJUSTED); $fixann.="$fx,";}
      elsif($fx =~ /^btrimto.(\d+)/) { $xto=$1; $xd=$xto-$b; $b = $xto unless($EXON_LOC_IS_ADJUSTED); $fixann.="$fx,";}
      elsif($fx =~ /^etrimto.(\d+)/) { $xto=$1; $xd=$xto-$e; $e = $xto unless($EXON_LOC_IS_ADJUSTED); $fixann.="$fx,";}
      elsif($fx =~ /^keep-gene-as-cds-split/ and $issplit) { $fixann.="splitgene,"; }
      else { $fixann.="fail:$fx,"; }
      }
    if( my($fl)= $xat =~ m/gapfill=([^;\s]+)/ ) { push @gapfill,$fl; }## upd 20150420
  
    #?? ONLY for exons, not CDS; apply gmoff only for b,e+-n or btrimto/etrimto ? where xd != 0
    # if($gmoff) { $b=$gmb if($b == $gbm0); $e= $gme if($e==$gme0); } # fixme, need this even when no generow for alts..
    if($gmoff>0 and $iexon==0) { if($b-$xd <= $gmb0){ $b += $gmoff; $fixann.="offby$gmoff,";} }
    elsif($gmoff<0 and $iexon==$nexon-1) { if($e - $xd >=$gme0){ $e += $gmoff; $fixann.="offby$gmoff,"; } }

    # *** FIXME: not EXON_LOC_IS_ADJUSTED or gmoff by1 change to exon b,e needs push back to $ex->[3,4] ***
    $ex->[3]= $b if($b0 != $b); $ex->[4]= $e if($e0 != $e);
    
    my $w= 1+$e-$b; $xlen += $w; $xdiff += $xd;
    foreach my $v ($b,$e) { 
      $xmax=$v if($v>$xmax); 
      $xmin=$v if($xmin<0 || $v<$xmin); 
      } 

      ## is this useful?
    if($b >= $cmin and $e <= $cmax) { # inside cds
      $cdse += $w; ## 1+$e-$b; # next;
    } elsif($b < $cmin and $e > $cmin) { # split utr
      $cdsb += 1 + $cmin - $b; 
      if($e>$cmax) {  $cdse += 1 + $cmax - $b; } else { $cdse += $w; } ##$cdse += 1 + $e-$cmin;
    } elsif($b < $cmax and $e > $cmax) { # split utr
      $cdse += 1 + $cmax - $b;  ##$utr3 += 1 + $e - $cmax;
    } else { # outside cds
      if($e <= $cmin) { $cdsb += $w; }
      elsif($b >= $cmax) { }
    }
    
    # push( @{$utr{$gid}}, [$b,$e,$o,$p,$gid,$r]);
  }
  
  if($to eq "-") { ($cdsb,$cdse)= (1+$xlen-$cdse,1+$xlen-$cdsb); }

  $nfix++ if(@gapfill);
  unless($nfix) {    
    if($genegap) { print "#ERROR: $gid genegap=$genegap but nothing  fixed\n";  } 
    elsif($debug) { print "#warn: $gid nothing  fixed\n"; }
    return putgene( $generecIN, $geneother);  
  }
  elsif($nfix and $generow and not $genegap) { # ok mistake? missing gene gap= flag
    print "#warn: $gid fixed but no gap= flag\n" if($debug);
  }
  
  # my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid)= @$mrna;
  if($cmin>0 and $cmin<$xmin) { } # problem need exon change
  if($cmax>$xmax) { } # problem need exon change

  if($xmin != $tb or $xmax != $te) {  # need mRNA, gene loc change
    my($tbo,$teo)= ($tb,$te);
    if($xmin != $tb) { $tb=$xmin; } 
    if($xmax != $te) { $te=$xmax; } 
 
    ## if($gmoff) { $tb=$gmb if($tb == $gbm0); $te= $gme if($te==$gme0); } # not needed if exons fixed

    $mrna->[3]= $tb; $mrna->[4]= $te;
    my $td= ($tb - $tbo) .",". ($te - $teo);
    $mrna->[8] =~ s/;gap=[^;\s]+//; $mrna->[8] =~ s/$/;gapfix=$td/;
    if($generow) {
      my($gb,$ge)= ($generow->[3], $generow->[4]); ## see above: gmb,gme
      $gb=$tb if($gb == $tbo and $gb != $tb);
      $ge=$te if($ge == $teo and $ge != $te);
      ($generow->[3], $generow->[4])= ($gb,$ge);
      $generow->[8] =~ s/;gap=[^;\s]+//; $generow->[8] =~ s/$/;gapfix=$td/;
    }
  }
  
#    ## not needed/bad if end exons fixed; from above; apply offby1 fix after other exon gapfix=
#   if($generow and $gmoff) { # and/or test gid = mRNA ID?
#     my($gmb,$gme)= ($generow->[3], $generow->[4]); ## see above: gmb,gme
#     if($gmoff>0) { $gmb += $gmoff; } elsif($gmoff<0) { $gme += $gmoff; }
#     ($generow->[3], $generow->[4])= ($gmb,$gme);
#   }
  
  $mrna->[8] =~ s/$/;gapfixan=$fixann/ if($fixann); # debug only?
  if(@gapfill) {
    my $gf=join",",@gapfill; 
    $mrna->[8] =~ s/$/;gapfill=$gf/;
  }
     
  if($dropgene) { 
    $mrna->[0] =~ s/^/#drop./;
    return putgene( [$mrna] );
  } else {
    return putgene( \@generec, $geneother); ## $generow, $mrna,
  }
}

#   my $location="$ref:$tb-$te:$to";
#   my $aalen="";
#   my $cxlen="";
#    
#   ($aalen,$cxlen)= aadef( $gid, $aalen, $cxlen);
#   if(my $ad= $aadef{$gid2}) {
#     # ad == "aalen=$al\txlen=$xl\toffs=$ofs";
#     ($aalen)= $ad=~m/aalen=(\S+)/; # this format is aasize,complete no %cds !!
#     my ($xlen)= $ad=~m/xlen=(\d+)/;
#     my ($coff)= $ad=~m/offs=(\S+)/;
#     my ($clen)= $ad=~m/aalen=(\d+)/; $clen= 3*(1+$clen);
#     $cxlen="$clen/$xlen,$coff";
#     my $pcds  = int(0.5 + 100*$clen/$xlen);
#     $aalen =~ s/^(\d+)/$1,$pcds%/ unless($aalen=~/,\d+%/);
#   }  
#
#   if(!$aalen and $tattr=~/aalen=([^;\s]+)/) { ## 249,76%,complete
#     my $aal=$1; # my @aal=split",",$aal;
#     # if(m/;aalen=(\d+),(\d+%),(\w+)/) 
#     $aalen=$aal if($aal =~ m/^\d+,\d+\%/);
#   } 
#   
#   if(!$cxlen and $tattr=~/cxlen=([^;\s]+)/) { ## 750/981,100-850
#     my $cxl=$1; # my @cxl=split /[,\/]/, $cxl;
#     $cxlen=$cxl if($cxl =~ m=^\d+[,/]\d+=);
#   }
#.......
#   $xlen||=1;
#   # these exon counts maybe bad for partial mapped genes..
#   my $aasize= int($clen/3); # -1 for stopcodon ?
#   my $pcds  = int(0.5 + 100*$clen/$xlen);
#   
#   unless($cxlen) { $cxlen= "$clen/$xlen,$cdsb-$cdse"; }  # is pcds on cxlen= or on aalen= now ??
#   else { 
#     my($c1,$x1)= $cxlen=~m=(\d+)/(\d+)=; 
#     if($x1-3 > $xlen) { 
#       my $px=int(0.5+100*$xlen/$x1); $px="0$px" if($px<10);
#       $cxlen.=",partmap$px" if($px<98);
#     } elsif($cdse and not($cxlen =~ /,\d/)) { $cxlen .= ",$cdsb-$cdse"; }
#   }
#   
#   unless($aalen) { $aalen= "$aasize,$pcds%"; } # complete or partial?
#   if($aalen !~ /complete|partial/ and $tattr=~/protein=([^;\s]+)/) { 
#     my $aa=$1; my $ac=0; # test for aacomplete
#     $ac|=1 if($aa =~ /^M/);  $ac|=2 if($aa =~ /\*$/ or $gid =~ /AUG/); # still missing stop for aug genes
#     $aalen .= ($ac == 3)?",complete":($ac == 1)?",partial3":($ac==2)?",partial5":",partial";
#   } 
#  
#   my $nintron= @exons - 1;
#   # full gene-attr.tab: ID osrc oid gene alttr aalen cxlen .. score location
#   print join("\t", qw(ID osrc aalen cxlen nintron location))."\n" unless($didhead++);
#   print join("\t", $gid, $src, $aalen, $cxlen, $nintron, $location)."\n";
  


sub _sortgene  
{
  # ($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid)
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
  my ($ng,$nx,$nr,$nsame,$nhit,$errord)= (0) x 10;
  my $nocomm= 1; ##($actid == ACT_NULL) ? 1 : 0;
  my @generec=();
  my @geneother=();
  my $geneid="";
  my $generow=undef;

# problems w/ input of partial records eg mRNA only, no gap : remove from input.gff
#  warn here if no gap=/gapfix= annot found?
  
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

    ## FIXME? need to collect all alt mRNA per gene record and fix all together, for genespan ? use generowGlobal

    if($typ =~ /^gene$/) { $generow= $rloc; }
    elsif($typ =~ /^($mrnatypes)$/) {  
      # allow gene and mRNA types .. and/or keep other types in generec
      $nsame += testgene($geneid, \@generec, \@geneother) if(@generec);
      
      $ng++;
      @generec = ($rloc); # maybe best store as array of [gffcols], sorted by type-loc
      @geneother= (); 
      if(ref $generow){ unshift @generec, $generow; $generowGlobal= $generow; }
      $generow=undef;
      $geneid= $gid;  # parse for gene vs tr id/alttr num ?
      
    } elsif($typ =~ /^($exontypes)$/) {
      if($pid ne $geneid) { warn "#ERR: Out-of-order GFF $typ:$pid in mRNA:$geneid\n"; $errord++; next; }
      push @generec, $rloc; $nx++; # check Parent == $geneid ?
      
    } elsif( ! $passtypes or "$typ.$src" =~ m/$passtypes/ ) {
      push @geneother, $rloc; # check Parent == $geneid ?
    }
  }
  
  $nsame += testgene($geneid, \@generec, \@geneother) if(@generec);
  return ($ng,$nx,$nsame,$nhit);
}



