#!/usr/bin/env perl
# intronokfilter.pl

=item about

  select missed coding loci from intronoktab.pl tables
  inputs: id table, intronok.tab for current best gene set + alt models
  
  idtable: $pd,$td,$src,$loc,$oids : src includes env bestsrc=best, loc=chr|scaffold:bbb-eee:or
  
  inok tabs now have CDS/UTR attribute
  collect valintrons for CDS exons, count alt models that 
    (a) add all-new cds, (b) extend best cds;
  use also antisense, long err flags to flag poor models

  select new loci by 'CDSok' now (valid CDS introns),  UTRok option ?
  
  ?? parse best/alt sets by ID?

=item evigene
  
  should become evigene/scripts/genes/intronokfilter.pl ? 
  -- next to intronoktab.pl
  -- needs options for classing best set vs alt set, other? or use --idtable bestgenes.idtab
      maybe: --idtab bestset.ids bestmodels*.intronok.tab  othermodels*.intronok.tab
    
=item usage
  
  intronokfilter.pl best3gs_alliup.idtab  *.intronok.tab > best_addloci.tab
  
=cut

use strict;
use Getopt::Long;

# my $BESTIDPRE=$ENV{bestid}||"Daplx7b3EV"; # not used..
my $BESTSRC= $ENV{bestsrc}||'beo3AUGcex2|beo3dpxevg5a|beo3alt'; # use 
my $LONGOK_MAX=9999; my $LONGOK_MAXin=2; 
my $PUTBEST= $ENV{showbest}||0;
my $PUTX= $ENV{showext}||0;
my ($ANTIOK,$idtab,$debug)=(0) x 9;

my $optok= GetOptions(
  "idtable=s", \$idtab, # list of best gene ids & source tags
  "bestsource=s", \$BESTSRC,  
  "longintronmax=i", \$LONGOK_MAX,
  "longcountintrons=i", \$LONGOK_MAXin,
  "showbest!", \$PUTBEST, 
  "showext!", \$PUTX, 
  "ANTIOK!", \$ANTIOK, 
  "debug!", \$debug, 
);

die "usage: intronokfilter.pl -idtable best3gs_alliup.idtab  *.intronok.tab > best_addloci.tab\n"
  unless($optok);
  
my(%pod,%tdsrc,%tdloc,%tdoid,%tclass,%tderr,%tdinok,%ind); 

#  class from best3gs_alliup.idtab
#  sources (col3) mains: 13845 beo3AUGcex2 14770 beo3dpxevg5a;  89218 beo3alt  ; 6091 evg7inew  2532 nomapsrc
sub classof { 
  my($td)=@_; 
  my $src= $tdsrc{$td} || $td; # use id pattern if no source field
  my $cl= ($src=~/$BESTSRC/)?"best":"other";
  return $cl; 
}

sub readIdrow {
  my($row)= @_;
  return unless($row=~/^\w/);
  my @v=split" ",$row; map{ $_||="" } @v[1..4];
  my($pd,$td,$src,$loc,$oids)=@v;
  $td=$pd unless($td); # list of ids only?
  $pod{$td}=$pd;
  $tdsrc{$pd}=$tdsrc{$td}= $src;
  $tdloc{$pd}=$loc;  # not used
  $tdoid{$pd}=$oids; # not used
}

if($idtab){ open(F,$idtab); while(<F>){ readIdrow($_); } close(F); }

while(<>) {
  next if(/^\W|^nogene/); 
  my @v=split;
  if($v[3]=~m/^(chr|scaffold|noloc)/) { # id table, 
    readIdrow($_); next;
  }
  
  # intronok.tab
  my($td,$ind,$inw,$iv,$xtype)=@v;
  $td=~s/Dpx6imEVm/Daplx7b3inewEVm/;  # ugh, idprefix mixup: idtab Daplx7b3inewEVm03716t1 == inoktab Dpx6imEVm03716t1
  $xtype ||= "other"; # bug?
  
  my($tclass)= classof($td);
  my($vi)= m/valintron=(\d+)/; 
  my($as)= m/antisense=(\d+)/; 
  my($ler)= ($vi<$LONGOK_MAXin and $inw > $LONGOK_MAX)?1:0;
  my $revind = 0;
  if($ANTIOK and $as>0 and $vi == 0) { 
    $vi= $as;  
    ##($or)=substr($revind,-1,1); $revind=substr($revind,0,-1); #??
    $revind=$ind; $revind=~s/(.)$//; my $or=$1; my $ror=($or eq "-")?"+":"-"; $revind.=$ror;
    # also need flipstrand on ind: s/:$or/:$revor/;
  } 
  my $iflag="";
  # $iflag.="inok," if($vi>0); # exclusive of errs # change to CDSok, UTRok
  $iflag.="${xtype}ok," if($vi>0); # exclusive of errs  
  $iflag.="antierr," if($as>0);
  $iflag.="longerr," if($ler>0);
  $iflag ||= "noevd";
  
  $tclass{$td}= $tclass;
  $tderr{$td}++ if($iflag=~/err/); # unless($iflag=~/inok/);
  $tdinok{$td}++ if($iflag=~/(CDS|UTR)ok/);
  
  # $ind{$td}{$ind}= [$iflag,$xtype,$inw,$iv]; # 
  $ind{$td}{$ind}= join"\t",$iflag,$inw,$iv,$xtype;
  $ind{$td}{$revind}= join"\t",$iflag,$inw,$iv,$xtype if($revind); #both? or only revind for ANTIOK?
}

my @ids= sort keys %tclass;
my @idbest= grep{ $tclass{$_} eq 'best' } @ids;
my @idalts= grep{ $tclass{$_} ne 'best' } @ids;
my (%incdsbest,%incdsnew,%incdshave); 
my ($nadd,$nextend)=(0,0);

for my $td (@idbest) { 
  my @in= sort keys %{$ind{$td}}; 
  my @cin= grep{ $ind{$td}{$_} =~ /CDSok/ } @in; 
  map{ $incdsbest{$_}++ } @cin;
}  
 
for my $td (@idalts) { 
  my @in= sort keys %{$ind{$td}}; 
  my @cin= grep{ $ind{$td}{$_} =~ /CDSok/ } @in; 
  for my $in (@cin) { my $ib= $incdsbest{$in}; if($ib){ $incdshave{$td}{$in}++; } else{ $incdsnew{$td}{$in}++; } }
}  

for my $td (sort keys  %incdsnew) {
  my @inew= sort keys %{$incdsnew{$td}};
  my @iold= sort keys %{$incdshave{$td}};
  if(@inew > @iold) {
    $nadd++;
    putgene($td, "locadd", \@inew, \@iold);
  } else {
    $nextend++;
    putgene($td, "locextend", \@inew, \@iold) if($PUTX);
  } 
}

if($PUTBEST) {
for my $td (@idbest) { 
  my @iold= sort keys %{$ind{$td}}; 
  putgene($td, "locbest", [], \@iold);
}
}

sub geneloc {
  my($ina)=@_;
  my($gr,$gb,$ge,$go)=(0) x 9;
  for my $in (@$ina) {
    my($r,$rb,$re,$ro)=split":",$in;
    unless($gr){ ($gr,$gb,$ge,$go)=($r,$rb,$re,$ro); }
    else { $gb=$rb if($rb<$gb); $ge=$re if($re>$ge); }
    }
  return join":",$gr,$gb,$ge,$go;
}

sub putgene {
  my($td, $group, $inew, $iold)=@_;
  my @in= sort keys %{$ind{$td}};  my $ni=@in;
  my $tclass= $tclass{$td};
  my $pd=$pod{$td} || $td;
  my $tdok = $tdinok{$td}||0;
  my $tderr= $tderr{$td}||0;
  my $nnew= @$inew;
  my $nold= @$iold;
  my $gloc= geneloc(\@in);
  print join("\t",$td,$group,$tclass,$pd,"cdsnewold=$nnew,$nold/$ni","okerr=$tdok,$tderr",$gloc)."\n";
  ## note: inew, iold only have valid-CDS subset, put UTR/invals..
  my %idid;
  for my $in (@$inew) {
    print join("\t",$td,"inew",$in, $ind{$td}{$in})."\n"; $idid{$in}=1;
  }
  for my $in (@$iold) { # include bestset id here?
    print join("\t",$td,"ihave",$in, $ind{$td}{$in})."\n"; $idid{$in}=1;
  }
  for my $in (grep{ not $idid{$_} } @in) {
    print join("\t",$td,"iother",$in, $ind{$td}{$in})."\n"; 
  }
}

__END__

# $nas=scalar(keys %inas); $nok=scalar(keys %inok); @td=sort keys %ind; $nt=@td; 
# print "# nt=$nt,nok=$nok, nanti=$nas\n"; 
# for $td (grep{ $inas{$_}} @td) { 
#   $as=$inas{$td}||0; $ok=$inok{$td}||0; print join("\t",$td,"$as.anti","$ok.ok")."\n";
# }  


=item inputs

id table
  best3gs_alliup.idtab
  sources (col3) mains: 13845 beo3AUGcex2 14770 beo3dpxevg5a;  89218 beo3alt  ; 6091 evg7inew  2532 nomapsrc
  ugh: idtab Daplx7b3inewEVm03716t1 == inoktab Dpx6imEVm03716t1
  
Daplx7b3EVm025199t1	Daplx6cgEVm032947t1	beo3dpxevg5a	scaffold_1:1216807-1217381:.	Daplx6suEVm044958t1,daplx6de3su9amvelvk43Loc14039t1
Daplx7b3EVm000051t1	AUGcex2s1g262t1	beo3AUGcex2	scaffold_1:1366957-1395701:+	AUGcex2s1g262t1,
Daplx7b3EVm000216t1	AUGcex2s1g218t1	beo3AUGcex2	scaffold_1:1150580-1182293:+	AUGcex2s1g218t1,
Daplx7b3EVm000316t1	AUGcex2s1g163t1	beo3AUGcex2	scaffold_1:905456-913699:-	AUGcex2s1g163t1,
..
Daplx7b3inewEVm03716t1	Daplx5cEVm035340t1_C2	evg7inew	scaffold_446:5680-33559:-	Daplx5cEVm035340t1_C2
Daplx7b3inewEVm00168t1	Daplx5cEVm003850t1	evg7inew	scaffold_116:272330-275128:+	Daplx5cEVm003850t1,daplx5ad9etu9nvelvk81Loc8829t1
Daplx7b3inewEVm00463t5	Daplx5cEVm035661t1	evg7inew	scaffold_42:646895-647418:+	Daplx5cEVm035661t1,daplx5ao11c9uvelvk57Loc17247t1
..
Daplx7b3nomapEVm01242t1	nomapid	nomapsrc	noloc	Daplx6suEVm051972t2

best inok
  
alt model inok

  ==> dpx6cex2aug_ann.intronok.tab <==
  AUGcex2s20g2t1	scaffold_20:13474:16069:-	2596	valintron=0	CDS
  AUGcex2s20g2t1	scaffold_20:16329:17160:-	832	valintron=0	CDS
  AUGcex2s20g2t1	scaffold_20:17238:18980:-	1743	valintron=0	CDS
  AUGcex2s20g4t1	scaffold_20:46127:46209:+	83	valintron=0	CDS
  
  ==> dpx7b3splign.intronok.tab <==
  sDaplx7b3EVm000021t1	scaffold_45:567124:567272:+	149	valintron=137	UTR
  sDaplx7b3EVm000021t1	scaffold_45:567348:567985:+	638	valintron=33	UTR
  sDaplx7b3EVm000021t1	scaffold_45:568092:568161:+	70	valintron=30	UTR
  sDaplx7b3EVm000021t1	scaffold_45:568368:568431:+	64	valintron=114	UTR
  
  ==> evg34okinmiss.intronok.tab <==
  Dpx6imEVm03716t1	scaffold_446:5727:5799:-	73	valintron=21	UTR
  Dpx6imEVm03716t1	scaffold_446:5900:33364:-	27465	valintron=0	UTR
  Dpx6imEVm00168t1	scaffold_116:272556:272620:+	65	valintron=205	CDS
  Dpx6imEVm00168t1	scaffold_116:273041:273111:+	71	valintron=262	CDS
  
  ==> velvdaplx6ml10dn9_inlocjtr_lo60s.inhitok.intronok.tab <==
  daplx6ml10dn9anvelvk31Loc1t5	scaffold_120:379885:379940:-	56	valintron=170	CDS
  daplx6ml10dn9anvelvk31Loc1t22_G2	scaffold_120:379885:379940:-	56	valintron=170	CDS
  daplx6ml10dn9anvelvk31Loc1t71	scaffold_42:765165:765231:-	67	valintron=55	CDS
  daplx6ml10dn9anvelvk31Loc1t123	scaffold_146:41258:41309:-	52	valintron=11	UTR


=cut

