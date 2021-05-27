#!/usr/bin/env perl
# joinsplitgff.pl

use strict;
# use warnings;
# use Getopt::Long;

my $debug=$ENV{debug}||0;
my $PRINTALL=$ENV{printall}||0;
my $ONLYCHANGES=$debug; # debug
   $ONLYCHANGES=0 if($PRINTALL);
   
my $mrnatypes='mRNA';
my $exontypes='exon|CDS';
my $passtypes="";
my ($ok,$input);

my $inh= *STDIN;
if($input) {
$ok = ($input =~ /.gz$/) ? open($inh,"gunzip -c $input |") 
      : ($input =~ /^(stdin|-)/) ? $inh= *STDIN
      : open($inh,$input);
die "bad -input=$input" unless($ok);
}

filter_gff($inh);

#----------------------

sub _sortgene  
{
  # ($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid)
  my($ta,$tb)= map{  (/gene/)?1:(/mRNA/)?2:(/CDS|exon/)?3:4; } ($a->[2],$b->[2]);
  return ($ta <=> $tb)
      || ($a->[0] cmp $b->[0]) # ref : should all be same for gene record
      || ($a->[3] <=> $b->[3]) # begin
      || ($b->[4] <=> $a->[4]) # end .. reverse larger > smaller
      # || ($a->[4] <=> $b->[4]) # end .. 
      # strand?
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
  my $geneid=""; my $lgeneid="";
  my $generow=undef;
  
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

    my $issplit= ($gid =~ m/_C(\d+)/)?$1 : (m/;Split=(\d+)/)?$1:0;
    my $gids= $gid;
    if($issplit) {
      $gid =~ s/_C\d+//; $pid =~ s/_C\d+//; 
    } else {
      print if($PRINTALL);
      next;
    }
    
    my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid,$issplit]; 

    if($typ =~ /^gene$/) { $generow= $rloc; }
    elsif($typ =~ /^($mrnatypes)$/) {  
      # allow gene and mRNA types .. and/or keep other types in generec
      
      if($gid eq $geneid and $issplit) {
        push @generec, $rloc; $ng++; # check Parent == $geneid ?
      
      } else {
      $nsame += testgene($geneid, \@generec, \@geneother) if(@generec);
      
      $ng++;
      @generec = ($rloc); # maybe best store as array of [gffcols], sorted by type-loc
      @geneother= (); 
      # if(ref $generow){ unshift @generec, $generow; $generowGlobal= $generow; }
      $generow=undef;
      $geneid= $gid;  # parse for gene vs tr id/alttr num ?
      }
      
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



sub putgene {
  my($generec, $geneother) = @_;
  $geneother ||= [];  my $nput=0;
  # note generec contains gene, mrna, exons, cds, in orig order; updated locs ..
  # my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid]; 
  foreach my $rloc (@$generec, @$geneother) {
    next unless(ref $rloc); # gene row may be undef
    print join("\t", @$rloc[0..8])."\n"; $nput++;
  }
  return $nput; # or nput?
}

sub testgene
{
  my($geneid, $generecIN, $geneother)= @_;
  
  my @generec= sort _sortgene @$generecIN; #? sort by genostart or 5'start?
  my($generow)= grep{ $_->[2] eq "gene" } @generec;  
  my @mrna = grep{ $_->[2] eq "mRNA" } @generec; # must have
  my $mrna = $mrna[0];
  my @cds  = grep{ $_->[2] eq "CDS"  } @generec;
  my @exons= grep{ $_->[2] eq "exon" } @generec;
  my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid,$issplit)= @$mrna;

  # my $issplit= ($mrna->[8] =~ /splitgene|Split=/)?1:0;
  # my ($genegap)= ($generow) ? $generow->[8] =~ /(gap=[^;\s]+)/ : ""; # dont always have this annot.. opt
  # ^^ test this vs nfix
  
  (my $gid2=$gid) =~ s/_[GC]\d+$//; #?
  #? $idlist{$gid}= $idlist{$gid2}= -9; # done flag

  my $nexon= @exons; 
  my($lex,$lxr,$lxo,$iexon)= (0) x 9;
  my($badspan,$cdscut)=(0) x 9; 
  my(@exok,@cdsok);
  foreach my $ex (@exons) {  
    $iexon++;
    # my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid,$issplit]; 
    my($xr,$xs,$xt,$xb,$xe,$xp,$xo,$tph,$xat)= @$ex;
    # if($lxo and $xo ne $lxo) { } #...
    if($lex and $lxr eq $xr and  $lex->[4] > $xb and $lex->[3] < $xe) {
      # check for xe > lex.xe, extend lex.xe; record as errspan in tr > geno map
      # push @ierr, $iexon;
      $lex->[8] =~ s/$/;errspan=$iexon:$xb-$xe/;
      my $xcut = 1 + $xe-$xb;
      if($xe > $lex->[4]) { my $d= $xe - $lex->[4]; $lex->[4]=$xe; $xcut -= $d; }  
      $badspan += $xcut;
    } else {
      push @exok, $ex;
      $lex= $ex;
    }
    $lxr=$xr; $lxo=$xo;
  }

  unless($badspan) {
  
    if($issplit and @mrna == 1) { # unsplit IDs .. have no 2nd path
      my($id)= $mrna->[8] =~ m/ID=([^;\s]+)/;  
      (my $nd=$id)=~s/_(C\d+)$//;
      $mrna->[8] =~ s/ID=$id/ID=$nd/; 
      # $mrna->[8] =~ s/(Split|chimera)=/${1}old=/g; # ugh.. want better attrib
      $mrna->[8] =~ s/;chimera=[^;\n]+//; # drop
      $mrna->[8] =~ s/;Split=/;unsplit=/; 
      for my $ex ( @exons, @cds) {
        my($id)= $ex->[8] =~ m/Parent=([^;\s]+)/;  
        (my $nd=$id)=~s/_(C\d+)$/;part=$1/;
        $ex->[8] =~ s/Parent=$id/Parent=$nd/; 
      }
      
    }
    return ($ONLYCHANGES) ? 0 : putgene( \@generec, $geneother); ## $generow, $mrna,

  } else {
    ($lex,$lxr,$lxo,$iexon)= (0) x 9;
    foreach my $ex (@cds) {  
      $iexon++;
      my($xr,$xs,$xt,$xb,$xe,$xp,$xo,$tph,$xat)= @$ex;
      if($lex and $lxr eq $xr and  $lex->[4] > $xb and $lex->[3] < $xe) {
        $lex->[8] =~ s/$/;errspan=$iexon:$xb-$xe/;
        my $ccut = 1 + $xe-$xb;
        if($xe > $lex->[4]) { my $d= $xe - $lex->[4]; $lex->[4]=$xe; $ccut -= $d; }  
        $cdscut+= $ccut;
      } else {
        push @cdsok, $ex;
        $lex= $ex;
      }
      $lxr=$xr; $lxo=$xo;
    }
    
    my($mb, $me)=(0) x 9;
    foreach my $ex (@exok) {  
      my($xr,$xs,$xt,$xb,$xe,$xp,$xo,$tph,$xat)= @$ex;
      $mb= $xb if( $mb == 0 or $xb < $mb); 
      $me= $xe if( $xe > $me );
      }
  
    #   my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid)= @$mrna;
    if(1) { # if($badspan or $mb != $tb or $me != $te)
      ## fixme: other mrna attr changes: add cov=  of @mrna split parts, nexon=@exok, nintron= ???, path=2/2?
      ##  remove? chim[12]=.. and chimera= or rename?
      ##$mrna->[8] =~ s/(Split|chimera)=/${1}old=/g; # ugh.. want better attrib
      $mrna->[8] =~ s/;(chim[12]|chimera)=[^;\n]+//g; # drop
      $mrna->[8] =~ s/;Split=/;unsplit=/; 
      my $nxok= @exok; $mrna->[8] =~ s/;nexon=\d+/;nexon=$nxok/; 
      my $covt=0; for my $m (@mrna) { if(my($c)=$m->[8] =~ m/;cov=(\d+)/) { $covt+=$c; } }
      $covt=100 if($covt>100); $mrna->[8] =~ s/;cov=\d+/;nexon=$covt/ if($covt); 
      ## change from  joinsplit= to unsplit= ??
      $mrna->[8] =~ s/$/;joinsplit=xcut:$badspan,ccut=$cdscut,oldspan:$tb-$te/;
      #  #err = "ERROR.span:genome_span:$gspan,tr_span:$tspan,$ref:$tb-$te";
     
      $mrna->[3]= $mb; $mrna->[4]= $me; 
      my($id)= $mrna->[8] =~ m/ID=([^;\s]+)/;  
      (my $nd=$id)=~s/_(C\d+)$//;
      $mrna->[8] =~ s/ID=$id/ID=$nd/; 
      }
    for my $ex ( @exok, @cdsok) {
      my($id)= $ex->[8] =~ m/Parent=([^;\s]+)/;  
      (my $nd=$id)=~s/_(C\d+)$/;part=$1/;
      $ex->[8] =~ s/Parent=$id/Parent=$nd/; 
    }
  
    my @newgenerec=( $mrna, @exok, @cdsok);
    return putgene( \@newgenerec, $geneother); ## $generow, $mrna,
  }
  
    
}