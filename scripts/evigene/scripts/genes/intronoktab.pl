#!/usr/bin/env perl
# evigene/scripts/genes/intronoktab.pl  introns.gff genes.gff > intronok.tab

# ** FIXME need location sort exons/gene for some genes.gff 
# .. still assume all exons follow mRNA, but may be unordered
# UPD1707: add errlongintron: invalid long introns; OFFBYCHECK splice-site off-by-1 = antisense err from gsplign

# use strict;
# use warnings;
use Getopt::Long;

my $IMIN=$ENV{imin}||1; 
my $debug=$ENV{debug}||0; 
my $FT=$ENV{ft}||"exon";  
my $ids=$ENV{ids}||"";
my $ANTICHECK= 1; # $ENV{anti}||0; # make default
my $OFFBYCHECK=0; #example of 'splign -direction wrongmrnadir' produces offby1 splice sites
my $LONGOK_MAX=9999; my $LONGOK_MAXin=2; #? not using MAXin here
my $SHOWERRIDS=0;
my($genegff,$introngff);
my @input; # == genes.gff
my $flagCDS=1;
my $ValidANTIOFF=0;  # count antisense= offby1= as valid intron hits (valid in some senses )
# my $IDFILTER=0; # id pattern, dont need with ids= 

my $optok= GetOptions(
  "ids=s", \$ids, # list of keep gene ids 
  "introns=s", \$introngff,
  "genes=s", \@input,
  "exontype=s", \$FT,
  "imin|MININWIDTH=i", \$IMIN,
  "ANTISENSE!", \$ANTICHECK, # -noanti turn off default
  "OFFBYCHECK!", \$OFFBYCHECK, #  
  "SHOWERRIDS!", \$SHOWERRIDS, #  
  "ValidANTIOFF!", \$ValidANTIOFF, #  
  "longintronmax=i", \$LONGOK_MAX, "longcountintrons=i", \$LONGOK_MAXin,
  "debug!", \$debug, 
);

unshift(@ARGV, @input) if(@input);
#x unshift(@ARGV, $introngff) if($introngff); # UNDO, read separately 1707

my $nkeepids= ($ids) ? readids($ids) : 0;

# push @input, $introngff if($introngff);
# push @input, @ARGV;
# my $inh=*STDIN;
# my $outh=*STDOUT;

if($introngff) {
  open(F,$introngff); 
  while(<F>){ 
    next if(/^\W/); chomp; @v=split"\t"; 
    my($rc,$s,$t,$rb,$re,$v,$ro,$x,$at)=@v; 
    #? next unless($t eq "intron" or $t =~ /$introntype/); # assume data only introns
    my $loc=join":",$rc,$rb,$re,$ro; 
    $in{$loc}=1+$v;
  } close(F);

}

# redo, collect gene record, sort, then calc introns
my($td,$pd,$nskip)=(0) x 9;
my(@xons,@cds); 
while(<>){
  next if(/^\W/); chomp; my @v=split"\t"; 
  my($rc,$s,$t,$rb,$re,$v,$ro,$x,$at)=@v; 
  if($t eq "intron"){  ## fixme need $introntype, or separate read -in; ? skip here if have above?
    my $loc=join":",$rc,$rb,$re,$ro; $in{$loc}=1+$v; 
  } else { 
    my($id)=m/(?:ID|Parent)=([^;\s]+)/; 
    ## mRNA .. add ncRNA , other transcript things w/ exons/introns
    ## other types trigger exon dump, new ID trigger exon dump
    
    ##drop this, use ids= file instead
    # if($IDFILTER) { # genes.gff subsets, eg t1 main .. could prefilter, or use idlist input?
    #  my $idok= ($id =~ m/$IDFILTER/)?1:0;
    #  unless($idok) { $nskip++ if($t eq "mRNA"); next; }
    #}
    
    if($t eq "mRNA"){ 
      $td=$id; #bad# $pd=0; @lx=(); 
      $exgene{$td}{0}++;  # count all ids here for no-intron genes 
    } elsif($t eq $FT){ 
    
      if($id ne $pd) { # exon ids only needed
        putxons($pd, \@xons, \@cds) if(@xons);
        @xons=(); @cds=();
      }
      $pd= $id; 
      push @xons, [@v];
      # == putxons
      # add uniq exon loc/geneid counts, to filter dupl mRNA as w/ corn-pacbio 
      # $xloc=join":",$rc,$rb,$re,$ro; 
      # $xloc{$xloc}=$pd unless($xloc{$xloc});
      # if(@lx) { 
      #   if($lx[3]>$rb) { ($ib,$ie)=($re+1,$lx[3]-1); } 
      #   else { ($ib,$ie)=($lx[4]+1,$rb-1); } 
      #   putx($pd,$rc,$ib,$ie,$ro); 
      #   }  
      # @lx=@v; 
    } elsif($t eq 'CDS' and $flagCDS) { 
      push @cds, [@v]; # keep to flag UTR/CDS for valintrons
    } 
    
    $lid=$id;
    } 
}
putxons($pd, \@xons, \@cds) if(@xons);

# END 
my @nogenehit;
if($ValidANTIOFF) {
  my %indidao=%indid; # indid
  map{ $indidao{$_}++; } keys %ingeneoff; # == indid for antisense,offby1
  @nogenehit= (grep{ not $indidao{$_} } sort keys %in);
} else { 
  @nogenehit= (grep{ not $indid{$_} } sort keys %in);
}

for $in (@nogenehit) { ## grep{ not $indid{$_} } sort keys %in
  $inv=$in{$in}; ($rc,$rb,$re,$ro)=split":",$in; $iw=1+$re-$rb; 
  print join("\t","nogenehit",$in,$iw,"cintron=$inv\n"); 
} 

sumstat();
#-----

sub readids {
  my($idf)=@_;
  ## note: these can be gene ids, where genes.gff has mRNA ids, diff suffix
  my $nk=0; open(F,$idf) or die "reading $idf";
  while(<F>){ next unless(/^\w/); my($id)=split; $nk++; $keepids{$id}=1; } close(F);
  return($nk); # return(\%keepids); #? or global
}

sub _sortgene  {
  # ($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr)
  my($ta,$tb)= map{  (/gene/)?1:(/mRNA/)?2:(/CDS|exon/)?3:4; } ($a->[2],$b->[2]);
  return ($ta <=> $tb)
      || ($a->[0] cmp $b->[0]) # ref : should all be same for gene record
      || ($a->[3] <=> $b->[3]) # begin
      || ($a->[4] <=> $b->[4]) # end 
      || ($a->[2] cmp $b->[2]) # CDS > exon
      ;
      ## [6] == orient, ignore? all same for generec
}

sub putxons {
  my($pd,$exons,$cds)=@_;
  if($nkeepids) {
    (my $gd=$pd)=~s/[\.t]\d+$//;
    return unless($keepids{$pd} or $keepids{$gd});
  }
  
  my (@lx,@sexons,@cloc);
  if($flagCDS and ref($cds)) {
    @sexons= sort _sortgene (@$exons,@$cds); # sort together, CDS>exon if same loc
  
  } else {
    @sexons= sort _sortgene @$exons;
  }

  for(my $i=0; $i<=$#sexons; $i++) {
    my @v= @{$sexons[$i]};
    my($rc,$s,$t,$rb,$re,$v,$ro,$x,$at)= @v; 
    
    if($flagCDS and $t eq 'CDS') {
      @cloc=($rc,$rb,$re,$ro);
    }
    next unless($t eq $FT);
    # add uniq exon loc/geneid counts, to filter dupl mRNA as w/ corn-pacbio 
    my $xloc=join":",$rc,$rb,$re,$ro; 
    $xloc{$xloc}=$pd unless($xloc{$xloc});
    
    if(@lx) { 
      if($lx[3]>$rb) { ($ib,$ie)=($re+1,$lx[3]-1); } 
      else { ($ib,$ie)=($lx[4]+1,$rb-1); }       
      my $iscds= ($flagCDS) ? sameloc( \@cloc, [$rc,$rb,$re,$ro]) : -1;
      putx($pd,$rc,$ib,$ie,$ro,$iscds); 
      }  
    @lx=@v; 
  }
}

sub sameloc{ my($c,$x)=@_;  return ($$c[0] eq $$x[0] and ($$c[1] eq $$x[1] or $$c[2] eq $$x[2])) ? 1:0;  }  

# #intronok.sum: nintron=189710, ingenevalid=101599,  inmiss=88111, hit%=53.5, genenotvalid=44680, inval%=nv/nin or nv/ing
# add ngenes w/ valid introns, ngenes w/o valid introns; rename genevalid=>exonvalid, genenotval=>exonnotval

sub sumstat {
  my @in= sort keys %in;
  my @ingene= sort keys %ingene; ## indid; = valid only
  my @gids= sort keys %exgene; 
  my @gidval= sort keys %exgeneval;
  my @goneexon= grep{ $exgene{$_}{0} } @gids; # == no-intron genes
  
  ## is this right fix?
  my %gxuniq=(); for $xl (sort keys %xloc) { my $gid=$xloc{$xl}; $gxuniq{$gid}++; }
  my @gxuniq= sort keys %gxuniq;
  @gids= grep{ $gxuniq{$_} } @gids;
  @gidval= grep{ $gxuniq{$_} } @gidval;
  @goneexon=grep{ $gxuniq{$_} } @goneexon;
  
#  # gene counts by ID dont correspond to location counts, due to many dupl, eg corn-pacbio set
#  # count by what? uniq exon/intron locations?
# cshlrnapb intronok: n intron=198622, exon=194967, genetr=477415, gene1exon=127260, inmiss=63156, hit%=68.2 
# cshlrnapb intronok: valid exon=135466, genetr=464342, noval exon=59501,30.5%, genetr=13073,2.7% / valgene= 97.3
#    ^^ gene counts wacky, due likely to many identical genetr, noval genetr is higher, same count but div dupl genetr*
# still wacky counts, more genetr (id) than exons (location)
# .. should not have more genetr (by id) than exons (by location)
# cshlrnapb intronok: n intron=198622, exon=194967, genetr=392946, gene1exon=89527, inmiss=63156, hit%=68.2 
# cshlrnapb intronok: valid exon=135466, genetr=380895, noval exon=59501,30.5%, genetr=12051,3% / valgene = 97
  
  my $nin= @in; my $ingene= @ingene;
  my $inmiss= scalar grep{ not $indid{$_} } @in;
  my $valing= scalar grep { $in{$_} } @ingene; # should == $nin - $inmiss

  if($ValidANTIOFF) {
    my %indidao=%indid;
    map{ $indidao{$_}++; } keys %ingeneoff;
    $inmiss= scalar grep{ not $indidao{$_} } @in;
    $valing= scalar grep{ $indidao{$_} } @in;
  }

  my $ngeneall= @gids; 
  my $geneval= @gidval; 
  my $ngeneonex= @goneexon;
  my $ngene= $ngeneall - $ngeneonex;
  my $geneinval= $ngene - $geneval;
  my $inval= $ingene - $valing;
  my $phit= int(1000 * $valing / $nin)/10;
  my $pnovalx= int(1000 * $inval / $ingene)/10; #? inval/nin or inval/ingene # this is exons now, change to gene counts?
  my $pnovalg= int(1000 * $geneinval / $ngene)/10; # genes w/o valid introns (includes 1x genes? or not)

   # gene-anti,amix=$naso,$nasmx; add OPTION to print anti,offby gids, type only,mix
  my $asense="";
  if($ANTICHECK) { 
    my @asgid= sort keys %asgeneval; # got wacky counts, >> antisense= intron rows; bad counter below: added non asense/novalint
    my @asmix= grep{ $exgeneval{$_} } @asgid;
    my @asonly= grep{ not $exgeneval{$_} } @asgid;
    my $nonly=@asonly; my $nmix=@asmix;
    $asense=" antisense,amix=$nonly,$nmix"; 
    if($SHOWERRIDS) {
      map{ print "#antisense_err\t$_\n" } (@asonly);
      map{ print "#antimix_err\t$_\n" } (@asmix);
      }
    }
  my $offby="";
  if($OFFBYCHECK) { 
    my @obgid= sort keys %offgeneval;
    my @obmix= grep{ $exgeneval{$_} } @obgid;
    my @obonly= grep{ not $exgeneval{$_} } @obgid;
    my $nonly=@obonly; my $nmix=@obmix;
    $offby=" offby1,omix=$nonly,$nmix"; 
    if($SHOWERRIDS) {
      map{ print "#offby1_err\t$_\n" } (@obonly);
      map{ print "#offmix_err\t$_\n" } (@obmix);
      }
    }
  my $longerr="";
  if($ANTICHECK or $OFFBYCHECK) {
    ##   $longerrgeneval{$pd}{$iloc}++ ; # count invalid long ints
    my @lergid= sort keys %longerrgeneval; # valid exgeneval not of interest
    my $nonly=@lergid;  
    $longerr=" longinterr=$nonly"; 
    if($SHOWERRIDS) {
      map{ print "#longint_err\t$_\n" } (@lergid);
      }
  }
   
  print "#intronok: n intron=$nin, exon=$ingene, genetr=$ngene (keep=$nkeepids), gene1exon=$ngeneonex, inmiss=$inmiss, hit%=$phit \n";
  print "#intronok: valid exon=$valing, genetr=$geneval, noval exon=$inval,$pnovalx%, genetr=$geneinval,$pnovalg%, $asense$offby$longerr  \n";
  # print "#intronok.sum: nintron=$nin, genevalid=$valing, inmiss=$inmiss, hit%=$phit, genenotval=$inval, notval%=$pnoval\n";
}


=item  UPD170713: OFFBYCHECK

  .. add offby1,2 test of intron/exon splices .. 
   cases in daphnia of such, may be biology or may be align soft problem: gsplign?
   eg: Daplx7b3EVm012731t1, revaa case, 1/2 exons have 1nt diff end from gmap vs gsplign
   .. aatrans same for both..
  # gmap (revmrna)
scaffold_720    dpx7b5  CDS     4983    5046    95      +
scaffold_720    dpx7b5  CDS     5114    5271    98      +
scaffold_720    dpx7b5  CDS     5348    5573    96      +
scaffold_720    dpx7b5  CDS     5636    5760    96      +
scaffold_720    dpx7b5  CDS     5824    5946    96      +
scaffold_720    dpx7b5  CDS     6017    6242    99      +
scaffold_720    dpx7b5  CDS     6318    6485    81      +
  # gsplign .. this offby *could* be splign method error, or antisense bug, or could be wobbly(?) splice site
scaffold_720    splign  CDS     4932    5046    0.954   -  : b-diff,e-same
scaffold_720    splign  CDS     5114    5271    0.987   -  = b,e-same
scaffold_720    splign  CDS     5348    5574    0.96    -  * -,e+1
scaffold_720    splign  CDS     5637    5761    0.96    -  * b+1,e+1
scaffold_720    splign  CDS     5825    5947    0.967   -  * b+1,e+1
scaffold_720    splign  CDS     6018    6243    0.991   -  * b+1,e+1
scaffold_720    splign  CDS     6319    6456    0.966   -  * b+1
  # valintrons of intron tab: gmap
scaffold_720:5047:5113:+        67      valintron=479   CDS:1,2
scaffold_720:5272:5347:+        76      valintron=1355  CDS:2,3
scaffold_720:5574:5635:+        62      valintron=1655  CDS:3,4
scaffold_720:5761:5823:+        63      valintron=572   CDS:4,5
scaffold_720:5947:6016:+        70      valintron=2098  CDS:5,6
scaffold_720:6243:6317:+        75      valintron=1219  CDS:6,7  
  # rnaintrons/intron3okids.gff | grep 'scaffold_720_5..._' | head
scaffold_720    intron  5047    5113    478     +   ID=scaffold_720_5047_5113_f;w=67
scaffold_720    intron  5231    7683    720     +   ID=scaffold_720_5231_7683_f;w=2453
scaffold_720    intron  5272    5347    1354    +   ID=scaffold_720_5272_5347_f;w=76
scaffold_720    intron  5574    5635    1654    +   ID=scaffold_720_5574_5635_f;w=62
scaffold_720    intron  5761    5823    571     +   ID=scaffold_720_5761_5823_f;w=63
scaffold_720    intron  5947    6016    2097    +   ID=scaffold_720_5947_6016_f;w=70
scaffold_720    intron  5951    6016    3       +   ID=scaffold_720_5951_6016_f;w=66
  # redo splign of revmrna, matches gmap splice ends .. problem is splign option -direction sense
  # .. redo all or subset of spligns w/o -direction?
cat Daplx7b3EVm012731t1revc.splign | grep '^+' | cut -f1,3-
+1      scaffold_720    0.954   175     1       175     4874    5046      <exon>GT      M65RM3RM11RM6RMDM9RMRMDM70
+1      scaffold_720    0.987   158     176     333     5114    5271    AG<exon>GT      M83RM5RM68
+1      scaffold_720    0.96    226     334     559     5348    5573    AG<exon>GT      M27RM2RM17RM18RM31RM38RM12RM13RM14RM45
+1      scaffold_720    0.96    125     560     684     5636    5760    AG<exon>GT      M5RM5RM41RM36RM4RM29
+1      scaffold_720    0.967   123     685     807     5824    5946    AG<exon>GT      M36RM5RM2RM44RM32
+1      scaffold_720    0.991   226     808     1033    6017    6242    AG<exon>GT      M75RM134RM15
+1      scaffold_720    0.966   266     1034    1298    6318    6583    AG<exon>GT      M27RM19RM38RM33IM97R2M10RM3RM17RM13
+1      scaffold_720    0.82    50      1299    1343    6629    6678    AT<exon>CT      M2I4M10RM3RM4IM16RM3RM3
.. rev dir has offby1: cat Daplx7b3EVm012731t1revc.splign | grep -v '^+' | cut -f1,3- | sort -k7,7n
-1      scaffold_720    0.954   175     175     1       5046    4874    AC<exon>        M71RDMRM9DMRM6RM11RM3RM65
-1      scaffold_720    0.987   158     333     176     5271    5114    AC<exon>CT      M68RM5RM83
-1      scaffold_720    0.96    227     560     334     5574*   5348    TA<exon>CT      M46RM14RM13RM12RM38RM31RM18RM17RM2RM27
-1      scaffold_720    0.96    125     685     561     5761*   5637*   CA<exon>CC      M30RM4RM36RM41RM5RM4
-1      scaffold_720    0.967   123     808     686     5947*   5825*   GA<exon>CC      M33RM44RM2RM5RM35
-1      scaffold_720    0.991   226     1034    809     6243*   6018*   CA<exon>CC      M16RM134RM74
-1      scaffold_720    0.966   265     1298    1035    6583    6319*   AC<exon>CC      M13RM17RM3RM10R2M97IM33RM38RM19RM26
-1      scaffold_720    0.82    50      1343    1299    6678    6629    AG<exon>AT      M3RM3RM16IM4RM3RM10I4M2  

=cut
  
sub putx { 
  my($pd,$rc,$rb,$re,$ro, $iscds)=@_; 
  my $iloc=join":",$rc,$rb,$re,$ro; 
  my $iw=1+$re-$rb; 
  my $inv= $in{$iloc}||0; 
  $ingene{$iloc}++; 
  $exgene{$pd}{$iloc}++; delete $exgene{$pd}{0}; # no longer no-intron gene

  if($inv) {
    $indid{$iloc}++ ;  ## indid == ingeneval
    $exgeneval{$pd}{$iloc}++ ; # count geneids, exons/gene and genehasvalidintrons
  } else {
    my($isanti,$isoffby,$islongerr)=(0,0,0);
    if($ANTICHECK) {
      my $or=($ro eq '-')?'+' : ($ro eq '+')?'-' :0;
      if($or) { 
        my $aloc=join":",$rc,$rb,$re,$or; 
        my $ainv= $in{$aloc}||0; 
        if($ainv) {
          ##?? opt to count as indid{} ? inerrdid{} ?
          $inv= "0,antisense=$ainv"; $isanti=1;
          $asgeneval{$pd}{$iloc}++ ; # count antisense
          $ingeneoff{$aloc}++;
         }
      }
    }
    if($OFFBYCHECK and not $isanti) {
      # example cases are offby1 both ends of exon vs intron; only shift both rb,re by +1, -1 ? or all permutes?
      # include antisense : usual case w/ gsplign d+1,o- or d-1,o+ ie revor
      # more offby? for $d (-1,1,-2,2){ for $xbe in ([$rb+$d,$re], [$rb,$re+$d], [$rb+$d,$re+$d] ) { .. } }
      my $or=($ro eq '-')?'+' : ($ro eq '+')?'-' :0;
      my ($oloc,$oinv);
      OFFBYLOOP: 
      for my $d (-1,+1) {
        my($orb,$ore)= map{ $_ + $d } ($rb,$re);
        for my $odir ($ro,$or) {
          $oloc=join":",$rc,$orb,$ore,$odir; 
          if($oinv= $in{$oloc}) {
          $inv= "0,offby=$oinv,d$d,o$odir"; $isoffby=1;
          $offgeneval{$pd}{$iloc}++ ; # count antisense
          $ingeneoff{$oloc}++;
          last OFFBYLOOP;
          }
        }
      }
    }
    if(not ($isanti or $isoffby) and $iw > $LONGOK_MAX) {
      $inv= "0,longerr=$iw"; $islongerr=1;
      $longerrgeneval{$pd}{$iloc}++ ; # count invalid long ints
      # my($ler)= ($vi < $LONGOK_MAXin and $iw > $LONGOK_MAX)?1:0;
    }
  }
  
  my $cpart=($iscds<0)?"0":($iscds>0)?"CDS":"UTR";
  print join("\t",$pd,$iloc,$iw,"valintron=$inv",$cpart)."\n" if($iw>=$IMIN); #? print if inv but not IMIN?
} 


# fixme for "." strand compare
sub _eqstrand { return($_[0] eq $_[1] || $_[0] eq '.' || $_[1] eq '.' )?1:0; }

# convert to array handling? _max(@list) 
sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }


=item calc
#...... intronok.tab calc : general tabulation of gene-exons.gff vs valid.introns.gff ......
#...... ? as evigene/scripts/genes/intronoktab.pl
#...... ? add sum at end, n,% found valid, nogene, invalids (val=0) ?

cat ../daplx16_erintron3okay.gff gsplign7{sets,merr}.gff |  env ft=exon perl -ne \
'next if(/^\W/); s/daphplx16ml_//g; @v=split"\t"; ($rc,$s,$t,$rb,$re,$v,$ro,$x,$at)=@v; if($t eq "intron"){ $loc=join":",$rc,$rb,$re,$ro; $in{$loc}=1+$v; } else { ($id)=m/(?:ID|Parent)=(\w+)/; if($t eq "mRNA"){ $td=$id; $pd=0; @lx=(); } elsif($t eq $FT){ $pd=$id; if(@lx) { if($lx[3]>$rb) { ($ib,$ie)=($re+1,$lx[3]-1); } else { ($ib,$ie)=($lx[4]+1,$rb-1); } putx($pd,$rc,$ib,$ie,$ro); }  @lx=@v; } } sub putx{ my($pd,$rc,$rb,$re,$ro)=@_; $iloc=join":",$rc,$rb,$re,$ro; $iw=1+$re-$rb; $inv=$in{$iloc}||0; $indid{$iloc}++ if($inv);  print join("\t",$pd,$iloc,$iw,"valintron=$inv\n") if($iw>=$IMIN); } END{ for $in (grep{ not $indid{$_} } sort keys %in) { $inv=$in{$in}; ($rc,$rb,$re,$ro)=split":",$in; $iw=1+$re-$rb; print join("\t","nogenehit",$in,$iw,"cintron=$inv\n"); } }  BEGIN{ $IMIN=$ENV{imin}||1; $FT=$ENV{ft}||"exon"; } ' \
  > gsplign7sm.intronok.tab

=cut

