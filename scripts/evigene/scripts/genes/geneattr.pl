#!/usr/bin/env perl
# geneattr.pl

=item about

  geneattr.pl -in gene.gff [-aadef genes.aadef]
  see also evigene annotation table collection of scripts, plus combiner bestgenes_puban_{xxx}.pl

=item example

  cat $pt.gmap.gff | $evigene/scripts/geneattr.pl -aadef=$pt.aadef > annot/$pt-attr.tab
  
  -aadef may be def header of aa seq file with added gene attributes: aalen, clen, offs
  
  $dmag/rnas/asmrna3/goodset/annot/daphmagna_201104m8-attr.tab
  ID      osrc    aalen   cxlen   nintron location
  m8PASAgasmblc_1 DGILmix8        37,18%,complete 111/616,258-144 4       contig00001:26-931:+
  m8AUGepir7_contig00003g154t1    DGILmix8        68,99%,partial5 206/206,0-206   1       contig00003:1-309:+
  m8AUGapi5_contig00031g159t1     DGILmix8        574,99%,partial5        1725/1725,1-1726        4       contig00031:1-1977:-

=item updates

  add input/output of attribs of unmapped genes (no gene.gff) as ID, unmapped .. here? or elsewhere?

  add further project attribute options, merging other gene-row tables (ID, attr1=xxx, attr2=yyy,..)
  eg.
    -attrs 'ID osrc oid gene alttr aalen cxlen nintron inqual ovrna ovtilex tere join must scorevec'
      ^^ input/output fields to use? as required output column field order?
    
  see also genes/score/geneset.ovxxx and score/geneset.genescore, scoremerge.sh
  for merging tables of (geneid, key=value, key2=value2, ...)
  -- allow mixed input of genes.gff + geneid,attrib tables ?
     cat genes.gff genes.*attr.tab | geneattr.pl -in stdin -idkey Nasvi -attrkeys 'aalen,cxlen,homolog,paralog,xxx,...'

=item extended geneattr

  FIXME to handle this sort of project-specific table set.
  
  cat  pub11us.gff | grep mRNA | cat pub11u.genescore - | perl -ne\
  'chomp; if(/^(Nasvi2\w+)/) { $gd=$1; ($g,$ta)= $gd=~/^(\w+)t(\d+)$/; $gat{$g}{ta}=$ta unless($gat{$g}{ta}>$ta);\
  map{ ($k,$v)=split"="; $gat{$gd}{$k}=$v; } grep /(nintron|inqual|tere|ovrna|ovpro|ovtilex)=/, split"\t"; next; } \
  ($gd)=m/ID=(\w+)/; ($g,$ta)= $gd=~/^(\w+)t(\d+)$/; $ta=0 unless($ta>1 or $gat{$g}{ta} > 1); \
  if(/;aaSize=\d/) { s/;(cxlen|aalen)=[^;\s]*//g; ($pp)=m/,Protein:(\w+)/; ($ap)=m/cdsSize=(\d+%)/; \
  s/cdsSize=$ap,/cxlen=/; $pp||="complete"; $ap||="50%"; s/aaSize=(\d+)/aalen=$1,$ap,$pp/; } \
  elsif(m/;aalen=(\d+),(\d+%),(\w+)/) { s/;aalen=\d+;/;/g; } else { warn "#NOaalen: $gd\n"; } \
  @v=split"\t"; @aa= split";",$v[-1]; @a= grep(/($KT)=/, @aa); \
  push @a,"alttr=$ta" if($ta); push @a,"score=$v[5]"; push @a,"location=$v[0]:$v[3]-$v[4]:$v[6]"; \
  %a=(); map{ ($k,$v)=split"=",$_; $a{$k}=$v unless($a{$k}); } @a; \
  $gid=$a{ID}; if($gat{$gid}) { map{ $v= $gat{$gid}{$_}; $a{$_}=$v if($v); } @KT; } \
  @a= map{ $a{$_}||"." } @KT; print join("\t",@a),"\n"; \
  BEGIN{ @KT=qw(ID osrc oid gene alttr aalen cxlen nintron inqual ovrna ovtilex tere join must scorevec); \
  $KT=join"|",@KT; push(@KT,"score","location"); print join("\t",@KT),"\n"; } ' \
  > pub11u-attr.tab
    
=item author
  
  don gilbert, gilbertd near indiana edu, ca 2011
  part of EvidentialGene, evigene/scripts/
  
=cut


use strict;
use warnings;
use Getopt::Long;

#my $TINYCDS = 40; # ignore tiny CDS tests for typeover == CDS
my $EXONSLOP = 0;
my $SAME_CDS= 90; # ** option
my $MINID_CDS= 0;
my $MINID_UTR= 0;
my $MISSING= "na"; # or "." ?

my $mrnatypes='mRNA';
my $exontypes='exon|CDS';
my $passtypes="";

my($input,$didhead,$overlaps,$debug,$aadef,$idlist,$ok,$idkey,$attrkeys)= (0) x 20;
## FIXME3: add input/output of unmapped (no gene.gff) as ID, unmapped .. here? or elsewhere?

my $optok= GetOptions(
  "gff|input=s", \$input,  
  "aadef=s", \$aadef,  
  "idlist=s", \$idlist,  
  "idkey=s", \$idkey, # for mixed input genes.gff, tables(idkey,attr)
  "attributes|attrkeys|columns=s", \$attrkeys,  
  "debug!", \$debug, 
  );

die "usage: geneattr.pl < genes.gff  > genes.table
  ** ASSUMES gene model is mRNA > exon > CDS
  ** ASSUMES input gff is ordered by gene records (mRNA/exon/CDS all together per ID)
" unless($optok); ##  and $input and $overlaps

my (%attrib, %aadef, %idlist);

my @ATCOL= qw(osrc aalen cxlen nintron location); # default/orig ; leave ID  out
my %ATDEF= (osrc=>$MISSING, nintron=>0, aalen=>$MISSING, cxlen=>$MISSING, location=>"unmapped");
my $ATDEFAULTS=1;

## extended genes-attr.tab used for pubgenes..
# @KT=qw(ID osrc oid gene alttr aalen cxlen nintron inqual ovrna ovtilex tere join must scorevec);
# nintron has mixed use: count exons-1 from genes.gff simple way, other use: validintron/-invalid/modelintron

#?? attrkeys eq "all" == collect from input?
if($attrkeys) { # output columns;
  # fixed keys: location, others? aalen, cxlen, nintron ?
  my @at= grep { /^\w/ && !/^(ID|location)$/ } split /[,|\s]+/, $attrkeys;
  @at= (@at,"location");  # leave ID out..
  $ATDEFAULTS= (@at eq @ATCOL)?1:0;
  if($ATDEFAULTS){ for my $i (0..$#ATCOL) { $ATDEFAULTS=0 unless($ATCOL[$i] eq $at[$i]); } }
  @ATCOL= @at;
}
my %ATCOL= map{ $_=>1 } @ATCOL;

if($idkey) { $idkey=~s/\d*$//; }

if($idlist) {
  open(F,$idlist) or die "missing idlist=$idlist";
  while(<F>) { my($id)=split; $id=~s/>//; $idlist{$id}=1; } close(F);
}

if($aadef) {
  open(F,$aadef) or die "missing aadef=$aadef";
  # >daphmag3tri7trimsub13loc1004c0t1 aalen=600,complete; clen=2213; strand=+; offs=201-2003;
  while(<F>) {
    my($id)=split; $id=~s/>//;
    my($al)=m/aalen=([^;\s]+)/;
    my($xl)=m/clen=(\d+)/; my($ofs)=m/offs=([^;\s]+)/;
    if($al) { $aadef{$id}="aalen=$al\txlen=$xl\toffs=$ofs"; }
  } close(F);
  map{ $idlist{$_}=2 } sort keys %aadef;
}

my $inh= *STDIN;
# allow for stdin mixing gff and non-gff tables of geneid,attr.key=value ?
if($input) {
$ok = ($input =~ /.gz$/) ? open($inh,"gunzip -c $input |") 
      : ($input =~ /^(stdin|-)/) ? $inh= *STDIN
      : open($inh,$input);
die "bad -input=$input" unless($ok);
}

## need separate output handler when:
##  NOT $ATDEFAULTS or  $attrib{$gid}{$k}= $v;

print join("\t", "ID", @ATCOL)."\n" unless($didhead++); # do outside
filter_gff($inh);

# add unmapped gene ids; look for @attr= @{ $attr{$id} } from non-gff inputs ??
foreach my $gid (grep{ $idlist{$_} > 0 } sort keys %idlist)
{
  outattr($gid);
}
#   my($src, $nintron, $aalen, $cxlen, $location)=("na",0,"na","na","unmapped"); # unmapped?
#   ($aalen,$cxlen)= aadef( $gid, $aalen, $cxlen);
#   print join("\t", $gid, $src, $aalen, $cxlen, $nintron, $location)."\n";
#  # $idlist{$gid}= -9; # done flag




#------------------------------------

sub outattr
{
  my($gid)= @_;
  
  my $atref= $attrib{$gid};
  unless(ref $atref) {  my %at= %ATDEF; $atref= \%at; }
  my($aagot,$aalen,$cxlen)= aadef( $gid);
  if($aagot) {
    $atref->{aalen}= $aalen;
    $atref->{cxlen}= $cxlen;
  }
  foreach my $k (keys %ATDEF) { $atref->{$k}= $ATDEF{$k} unless(defined $atref->{$k}); }
  
  my @attr= map{ defined $atref->{$_} ? $atref->{$_} : $MISSING } @ATCOL;
  print join("\t",$gid,@attr)."\n";
  #ATDEFAULTS: print join("\t", $gid, $src, $aalen, $cxlen, $nintron, $location)."\n";
}

sub aadef
{
  my($gid,$aalen,$cxlen)= @_;
  (my $gid2=$gid) =~ s/_[GC]\d+$//;
  my $aagot=0;
  my $ad= $aadef{$gid} ||  $aadef{$gid2};
  if($ad) {
    # ad == "aalen=$al\txlen=$xl\toffs=$ofs";
    ($aalen)= $ad=~m/aalen=(\S+)/; # this format is aasize,complete no %cds !!
    my ($xlen)= $ad=~m/xlen=(\d+)/;
    my ($coff)= $ad=~m/offs=(\S+)/;
    my ($clen)= $ad=~m/aalen=(\d+)/; $clen= 3*(1+$clen);
    $cxlen="$clen/$xlen,$coff";
    my $pcds  = int(0.5 + 100*$clen/$xlen);
    $aalen =~ s/^(\d+)/$1,$pcds%/ unless($aalen=~/,\d+%/);
    $aagot= ($aalen =~ /\d/)?1:0;
  }  
  return($aagot,$aalen,$cxlen);
}

sub testgene
{
  my($geneid, $generecIN, $geneother)= @_;
  
  my @generec= sort _sortgene @$generecIN; #? sort by genostart or 5'start?
  # my $generec= \@generec;
  
  my($mrna)= grep{ $_->[2] eq "mRNA" } @generec;
  my($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid)= @$mrna;
  my $location="$ref:$tb-$te:$to";
  my $aalen="";
  my $cxlen="";
    
  #return 0 if($skippatt and $tattr =~ m/$skippatt/);
  my @cds  = grep{ $_->[2] eq "CDS" } @generec;
  my @exons= grep{ $_->[2] eq "exon" } @generec;

## FIXME: using exons to measure CDSsize, trsize is bad for poormap genes, eg 30% mapped.
## use aadef in preference (from findcds from tr)
## FIXME2:  gid = nnn_C1,2 nnn_G2,3.. chomp map flag

  (my $gid2=$gid) =~ s/_[GC]\d+$//;
  my $aagot=0;
  ($aagot, $aalen,$cxlen)= aadef( $gid, $aalen, $cxlen);

  if(!$aalen and $tattr=~/aalen=([^;\s]+)/) { ## 249,76%,complete
    my $aal=$1; # my @aal=split",",$aal;
    # if(m/;aalen=(\d+),(\d+%),(\w+)/) 
    $aalen=$aal if($aal =~ m/^\d+,\d+\%/);
  } 
  
  if(!$cxlen and $tattr=~/cxlen=([^;\s]+)/) { ## 750/981,100-850
    my $cxl=$1; # my @cxl=split /[,\/]/, $cxl;
    $cxlen=$cxl if($cxl =~ m=^\d+[,/]\d+=);
  }

  my ($cmin,$cmax,$clen,$xlen,$cdsb,$cdse)= (-1,0,0,0,0,0);
  foreach my $ex (@cds) {
    my($b,$e)=@{$ex}[3,4]; $clen+= 1+$e-$b;
    foreach my $v ($b,$e) { 
    $cmax=$v if($v>$cmax); 
    $cmin=$v if($cmin<0 || $v<$cmin); 
    } }

  foreach my $ex (@exons) {
    my($xr,$xs,$xt,$b,$e,$p,$o)= @$ex;
    my $w= 1+$e-$b; $xlen += $w;
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
    # push( @{$utr{$tid}}, [$b,$e,$o,$p,$gid,$r]);
  }
  
  if($to eq "-") { ($cdsb,$cdse)= (1+$xlen-$cdse,1+$xlen-$cdsb); }
    
  my $nintron= @exons - 1;

  $xlen||=1;
  # these exon counts maybe bad for partial mapped genes..
  my $aasize= int($clen/3); # -1 for stopcodon ?
  my $pcds  = int(0.5 + 100*$clen/$xlen);
  
  # my $cxlen= "$clen/$xlen,$pcds%,$cdsb-$cdse";  # is pcds on cxlen= or on aalen= now ??
  #?? add flag to cxlen when exon sizes << aadef xlen eg: ,partmap30% ?
  unless($cxlen) { $cxlen= "$clen/$xlen,$cdsb-$cdse"; }  # is pcds on cxlen= or on aalen= now ??
  else { 
    my($c1,$x1)= $cxlen=~m=(\d+)/(\d+)=; 
    if($x1-3 > $xlen) { 
      my $px=int(0.5+100*$xlen/$x1); $px="0$px" if($px<10);
      $cxlen.=",partmap$px" if($px<98);
    } elsif($cdse and not($cxlen =~ /,\d/)) { $cxlen .= ",$cdsb-$cdse"; }
  }
  
  unless($aalen) { $aalen= "$aasize,$pcds%"; } # complete or partial?
  if($aalen !~ /complete|partial/ and $tattr=~/protein=([^;\s]+)/) { 
    my $aa=$1; my $ac=0; # test for aacomplete
    $ac|=1 if($aa =~ /^M/);  $ac|=2 if($aa =~ /\*$/ or $gid =~ /AUG/); # still missing stop for aug genes
    $aalen .= ($ac == 3)?",complete":($ac == 1)?",partial3":($ac==2)?",partial5":",partial";
  } 
  
  # ? defer this output? maybe hash{gid} = attrib{key=val} ?
  # gene-attr.tab: ID osrc oid gene alttr aalen cxlen .. score location
  if($ATDEFAULTS and not $idkey) {
    ## print join("\t", qw(ID osrc aalen cxlen nintron location))."\n" unless($didhead++);
    print join("\t", $gid, $src, $aalen, $cxlen, $nintron, $location)."\n";
    $idlist{$gid}= $idlist{$gid2}= -9; # done flag
    
  } else {
    $idlist{$gid}= 3; # NOT done flag; idlist val = 1,2  above
    $attrib{$gid}{osrc}=$src;
    $attrib{$gid}{aalen}=$aalen;
    $attrib{$gid}{cxlen}=$cxlen;
    $attrib{$gid}{nintron}=$nintron;
    #after# $attrib{$gid}{location}=$location;
    $attrib{$gid}{score}= $tscore;
    
    ## map other tattr values per @ATCOL ; possibly replace above vals ?
    foreach my $at (grep /=/, split /;/, $tattr) {
      my($ak,$av)= split"=",$at,2;
      $attrib{$gid}{$ak}= $av if($ATCOL{$ak});
    }
    
    $attrib{$gid}{location}=$location; # dont allow replacement
  }
  
  return 1;
}


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
  my ($ng,$nx,$nr,$nsame,$nhit,$errord,$inerr)= (0) x 10;
  my $nocomm= 1; ##($actid == ACT_NULL) ? 1 : 0;
  my @generec=();
  my @geneother=();
  my $geneid="";
  
  while(<$inh>) {
    unless(/^\w/){ print unless($nocomm or /^(#n |$)/); next; }
    chomp; my @v= split"\t";
    my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr);
    
    # here allow for stdin mixing gff and non-gff tables of geneid,attr.key=value ?
    if($idkey and /^$idkey/) {
      my $gid= shift(@v);
      foreach my $kv (grep /\w=./, @v) {
        my($k,$v)= split "=", $kv, 2;
        $attrib{$gid}{$k}= $v; # allow multiples or not? .= "$v,"
      } 
      next;   
    } elsif(@v == 9) {
      ($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr)= @v;
    } else { # bad gff? skip?
      warn "# BAD input, not gene gff or attr: $_\n" unless($inerr++ > 10); next;
    }
    
    $nr++;    
    my($gid,$pid); 
    if($tattr =~ m/\bID=([^;]+)/) {  $gid=$1; }
    if($tattr =~ m/\bParent=([^;]+)/) {  $pid=$1; 
      $gid=$pid unless($gid and ($typ =~ /^($mrnatypes)$/)); 
      }
    unless(defined $gid) { $gid = "N".$ng; }

    my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid]; 

    # if($typ =~ /^gene$/) {}
    if($typ =~ /^($mrnatypes)$/) {  
      # allow gene and mRNA types .. and/or keep other types in generec
      $nsame += testgene($geneid, \@generec, \@geneother) if(@generec);
      
      $ng++;
      @generec = ($rloc); # maybe best store as array of [gffcols], sorted by type-loc
      @geneother= ();
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



