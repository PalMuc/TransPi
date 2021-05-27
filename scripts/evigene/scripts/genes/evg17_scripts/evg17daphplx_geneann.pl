#!/usr/bin/env perl
# genescore_dpx7b.pl
 

=item
 
 for Daphnia_pulex 2017 Evigene set annotation tables
 merge gene annotation tables for public gene set annots
 - per table parsing
 
 .. add here? gff, seq files w/ updated pubid, annots
 
 See plants/arabidopsis//evg5weed/annob/evg5arathann.pl 2017.may
 .. same methods
 
=cut

use strict;

my $anpath="evg5daplx/map5ann";


my @ANTABS= qw(
best3gs_alliup_mrna.attr 
best3gs_alliup_homol.sumtab1
best3gs_alliup_hoval.tab 
dpx7b7fcdsn.names 
dpx7n7pub.intronsum.tab dpx7n7pub_newing.intronsum.tab
dpx7n7pub.iste.tab 
dpx7n7pub.noncode.rereclass.tab
dpx7n7pub_cullnc.align.tab dpx7n7pub_cullnc.masource.uptab
dpx7n7pub_newing.align.tab dpx7n7pub_newing.pubids 
dpx7n7pub_updbyhand.tab 
);

my $NEWANNTAB = $ENV{"anntab"} || "dpx7n7pub_upd.genescore";
my $ORGANISM= $ENV{organism}||"Daphnia_pulex";
my $IDPRE="nDaplx7b3EVm"; ##current nDaplx7b3EVm; change to Daplx7b3EVm (oldpub) or Daplx7pEVm ?
        #  "Arath5nEV" used as $newid=~s/^Arath5\w+EV/$IDPRE/; 
        
my $debug=$ENV{debug} || 1;
my $oclass=$ENV{class}||""; # keep these
my $noclass=$ENV{noclass}||""; # skip these
my $GFFSRC=(defined $ENV{gffsrc}) ? $ENV{gffsrc} : undef;
my $CUTOIDS  = 'AUGcex\w+(t[2-9]+|t1\d+)|Daplx6cgEVm\w+'; # drop these oids by pattern .. others?
my $SEQTYPE=$ENV{seqtype}||""; #  seq type=protein|CDS|mRNA|ncRNA ..

my(%annot, %pubids, %oldpubid, %dop, %pod, 
  %anntab, %anncols, %tod, %nod); # read_newanntab

sub MAINstub {}

if(my $seq=$ENV{seq}) { # seq=inputs.aa .. use stdin?
  read_newanntab("forseq",$NEWANNTAB); # *STDIN
  putseq($seq); 
  
} elsif(my $gff=$ENV{gff}) {
  read_newanntab("forgff",$NEWANNTAB); # *STDIN
  putgff($gff); 

} elsif(1) {

  readMrnaAttr("best3gs_alliup_mrna.attr"); # read_anntxt()

  readTable("dpx7n7pub.noncode.rereclass.tab","Noncode|Code|Unknown");
  readTable("dpx7n7pub.iste.tab", "iste=");
  readTable("dpx7n7pub.intronsum.tab","ints=");
  readTable("dpx7n7pub_newing.intronsum.tab","ints=");
  
  readHoval("best3gs_alliup_hoval.tab","xxx=");
  readTable("best3gs_alliup_homol.sumtab1","hospp="); # skip hoval= if have hoval.tab ?
  readNames("dpx7b7fcdsn.names","xxx="); # read_names(); use as primary name source? ie replace mrna.attr ?
  
  readAligntab("dpx7n7pub_cullnc.align.tab"); # read_mapattr()
  readAligntab("dpx7n7pub_newing.align.tab"); # read_mapattr()
  
  readIdtable("dpx7n7pub_cullnc.masource.uptab"); # read_newidtab()
  readIdtable("dpx7n7pub_newing.pubids"); # read_newidtab()

# ** FIXME have some dups b/n id tables, same oid, some are cases of diff mapping locus, others are overlapped/same locus
  
  readReclasstab("dpx7n7pub_updbyhand.tab");
  
  putGenescore($NEWANNTAB); # == write_newanntab()
}

=item main() from evg5arathann.pl

if($seq=$ENV{seq}) { 
  read_newanntab("forseq"); # *STDIN
  putseq($seq); 
  
} elsif($gff=$ENV{gff}) {
  read_newanntab("forgff"); # *STDIN
  putgff($gff); 

} elsif(1) { # $ENV{anntab} ?

  read_newidtab("evg5arath_icnannbest.idtab");
  read_keeps("evg5arath.keepids");
  read_anntxt("evg5arath.ann.txt");
  read_puboids("evg5arath.puboids");
  read_names("evg5arath.names");
  read_alntab("evg5arath.cacaorange.alntab");
  read_mapattr("evg5arath_mrna_input.map.attr");
  
  write_newanntab();
}

=cut

#------------------------------------------

# arath > dpx: mapclass => Class
# my $mapclass = 'Class';
# my $origid   = 'obID'; # oid or obID
# my $locus    = 'locus'; #or 'maploc'; # change back to locus? YES
# my $namealn  = 'napct'; # this AND Dbxref ? revert napct => namealn ?

sub _sortid { my($ag,$ai)= $a=~m/EVm(\d+)t(\d+)/;  my($bg,$bi)= $b=~m/EVm(\d+)t(\d+)/; 
  return ($ag <=> $bg or $ai <=> $bi or $a cmp $b); }

sub putseq {
  my($seq)=@_;
  warn "#annot of $seq\n";
  
  #* drop all in seq attr but for select set, add from %anntab
  # keepattr: aalen  srcID 
  # rewrite from annot: oid= (subset? 1st, last only?) ; cutoids='AUGcex\w+(t[2-9]+|t1\d+)|Daplx6cgEVm\w+'
  # rewrite : organism=Daphnia_pulex  type=protein/CDS/mRNA/ncRNA .. other?
  # my $SEQTYPE=$ENV{seqtype}||"na"; # for seq type=protein|CDS|mRNA|ncRNA ..
  
  # special case seq tags: asrc=beo3alt,mapaa; srcID=nDaplx7b3EVm000001t16;srcFn=dpx7n7pub.gff2aa; srcOK=ok;
  # .. keep 'mapaa' somewhere
  
  my @keepseqann= qw( type aalen  ); # srcID ??
  my @addann= qw( aalen locus Dbxref Name ); # Class ; seq type=protein/CDS/mRNA/ncRNA 
  ## keep: type, aalen, clen, offs, organism
  # my @seqann= ( 'locus', 'Class', 'oid','Name' );# $origid, $namealn 
  # my @keepann=qw(type aalen clen offs organism);
  # my @dropann=qw(uvcut uvfa evgclass ); my $dropann=join('|',@dropann); # others? oid 
  #above#$oclass=$ENV{class}||"";

  my ($ok,$no,$inf)=(0,0,undef);
  # outfile instead of stdout?
  if($seq =~ /^stdin|^-/){ $inf= *STDIN; }
  else { open(S,$seq) or die $seq; $inf= *S; }
  
  while(<$inf>){
    if(/^>(\S+)/) { 
      my $td=$1; my $nd=0;
      if(my $d=$nod{$td}) { $nd=$td; $td=$d; }
      unless($nd) { $nd= $tod{$td}; } # || $tod{$tdc}; 
      
      $ok=($nd)?1:0; 
      ## dropdup is noclass; keep _frag or not? segregate
      if($ok and $noclass) { $ok=($anntab{$nd}{'Class'} =~ /$noclass/) ? 0: 1; }
      elsif($ok and $oclass) { $ok=($anntab{$nd}{'Class'} =~ /$oclass/) ? 1: 0; }
      if($ok) {
        my $hdr=$_; $hdr=~s/>\S+//; chomp($hdr);
        # my $oid=$anntab{$nd}{'obID'}||"noid";
        my $oid= $anntab{$nd}{'oid'}||"noid"; $oid=~s/\b($CUTOIDS),//g;
        unless($oid=~m/\b$td\b/ or $td eq $nd){ $oid="$td,$oid"; }
        
        my $cla= $anntab{$nd}{'Class'}; # must exist?
        my $ncla= reClassify($cla); 
        my $stype= ($hdr =~ m/\b(type=[^;\n]+)/)? $1 : "";
        my $aalen= ($hdr =~ m/\b(aalen=[^;\n]+)/)? $1 : "";
        my $asrc = ($hdr =~ m/\b(asrc=[^;\n]+)/)? $1 : "";
        
        # fixme: bad type=exon... or type=input RNA or Class=noncode type=mRNA
        #  >nDaplx7b3EVm000002t17 type=exon.dpx7n7_ar; aalen=6139,84%,complete; aafrommap=1; 
        #  evgclass=alt,dpx7n7_ar; locus=scaffold_29:324312-352460:+; Dbxref=daphsim:Dapsim1EVm000004t1; Name=Nesprin-1; 
        # oid=Daplx7b3EVm000076t17,AUGcex2s29g213t20,Daplx6mlEVm000002t13,daplx6ml634a9amvelvk77Loc838t33; organism=Daphnia_pulex
        #   nDaplx7b3EVm000079t6 type=RNA; aalen=1768,86%,complete; evgclass=alt,dpx7n7_ar; locus=scaffold_12:1091947-1098225:+; oid=daplx6ml10dn9anvelvk35Loc3395t8;
        # nDaplx7b3EVm058775t1 type=mRNA; aalen=98,5%,complete-utrbad; evgclass=noncode,dpx7n7_nc; ..
        
        if($stype) {
          #  type=exon.dpx7n7_ar; type=CDS.dpx7n7_ar; with aafrommap
          $stype=~s/=cds/=CDS/;
          if($stype=~m/=exon/) { $stype="type=mRNA"; } # for aafrommap=1;
          elsif($stype=~m/=CDS./) { $stype="type=CDS"; } # for aafrommap=1;
          if($stype=~m/=mRNA/ and $ncla =~ /^noncode|^altnc/){ $stype="type=ncRNA"; } # altnc or noncodealt ?
        } elsif($SEQTYPE) { 
          $stype= "type=$SEQTYPE"; 
        }
        ##no # my @keeps= map{ ($hdr =~ m/\b($_=[^;\n]+)/)? $1 : ""; } @keepseqann;
        my @keeps;
        push @keeps, $stype, $aalen, "evgclass=$ncla,$cla";
        push @keeps,"aafrommap=1" if($asrc =~ m/,mapaa/);
        my $hasaa= ($aalen =~ m/=/);
        my @adds = map{ my $v=$anntab{$nd}{$_}; $v="" if($hasaa and $_ eq "aalen"); ($v)?"$_=$v":"" } @addann;
        my $addann = join "; ", grep(/=/, @keeps,  @adds, "oid=$oid", "organism=$ORGANISM");
        $nd=~s/_[CG]\d+$//; 
        s/^>\S+.*$/>$nd $addann/; $no++;
      }
    }
    print if($ok); # print $outh $_ if($ok);
  } close($inf);
  warn "#seq output n=$no\n";
  return($no);
}


sub putgff {
  my($gffin)=@_;
  warn "#annot of $gffin\n";

  # FIXME arath>dpx attr keys  
  # FIXME: drop most gff annots, add most from %anntab
  # namealn => Dbxref, napct .. napct => namealn ?

  my @keepmrnan= qw( Split Target trg msrc  );  # replace cov pid  with mapqual;  gaps from splign only?
  my @keepexonn= qw( Target trg ix mgap error );  # fix mgap,error > aligngap=
  my @addann= qw( aalen intron intronerr mapqual homolog hokaks Dbxref napct Name ); #  iste? ; evgclass=

  #... old ...  
  my @gffann= grep{ not m/origid|locus|aaqual|cdsoff|mapqual/ } 
    sort { $anncols{$a} <=> $anncols{$b} or $a cmp $b } keys %anncols;
  my @dropann= qw(aamap match qlen clen chimera cdsindel); # qw( aasize offs ) ? 
  push @dropann, qw(pid chim1 chim2); # mapqual; leave cov,nexon for other uses
  my $dropann= join('|',@dropann);
    
  my ($ok,$no,$inf)=(0,0,undef);
  if($gffin =~ /^stdin|^-/){ $inf= *STDIN; }
  else { open(S,$gffin) or die $gffin; $inf= *S; }

  # FIXME: nDaplx7b3EVm064001t1	obID=daplx6ml10dn9anvelvk31Loc1t350	Class=vel10in_nc
  #   .. change mapclass  vel10in to dpx7n7, same as others
  #   gff ID='u' + obID
  
  while(<$inf>){
    if(/^\w/){ 
      chomp; my @gff=split"\t";
      my($gsrc,$gtype,$gann)= @gff[1,2,8]; chomp($gann);
      my($IDkey,$IDval)= m/(ID|Parent)=([^;\s]+)/; # _C1,2 split and _G2,3.. dup map ids?

      my $commout="";
      my($ok,$nd,$tdc)= (0) x 9;
      my $td=$IDval; #changed
      $td =~ s/^udaplx/daplx/; # vel10in_nc dpx7n7pub_newing.gff fixup .. NOT WORKING
      ($tdc=$td)=~s/_[CG]\d+$//;
      # $tdc=~s/^udaplx/daplx/; # vel10in_nc dpx7n7pub_newing.gff fixup .. NOT WORKING
      
      if(my $d=$nod{$td}) { $nd=$td; $td=$d; }
      elsif($tdc ne $td and my $d=$nod{$tdc}) { $nd=$tdc; $td=$d; }
      ($tdc=$td)=~s/_[CG]\d+$//;
      unless($nd) { $nd= $tod{$td} || $tod{$tdc}; }
      
      $ok=($nd)?1:0; 
      if($ok and $noclass) { $ok=($anntab{$nd}{'Class'} =~ /$noclass/) ? 0: 1; }
      elsif($ok and $oclass) { $ok=($anntab{$nd}{'Class'} =~ /$oclass/) ? 1: 0; }
      
      if($ok) {
        my $mapclass= $anntab{$nd}{'Class'}; # must exist?
        my $reclass= reClassify($mapclass); 
        
        ## dang hack; fix this in dpx7n7pub_upd.genescore table
        my $gannadd=$gann;
        if($mapclass =~ /vel10in_/) { $gannadd .= ";msrc=$mapclass"; $mapclass =~ s/vel10in_/dpx7n7_/; } 
        
        if($gtype eq "mRNA") { 
          my $obid= $anntab{$nd}{'obID'}||""; # same as td ??
          my $oid = $anntab{$nd}{'oid'}||""; # check, remove dup oids..
	        my ($oidg)= (m/;oid=([^;\s]+)/)?$1:""; # maybe ignore gff oid=?
	        my %oid; my $i=0; 
	        map{ if(/\w/ and not m/\b($CUTOIDS)/) { $oid{$_}= ++$i; } } split",","$tdc,$obid,$oid,$oidg"; 
          $oid= join",",  sort{ $oid{$a}<=>$oid{$b} } keys %oid;
          
           ## have mapqual in genescore table now
	        my $mapq= $anntab{$nd}{mapqual}||0; # mapqual=97a,99i,14787l,62x
	        $mapq=~s/a,/c,/; ##  sub mapqual() 'c' not 'a': 97c,99i,62x .. no 'l' for clen; change back?
	        unless($mapq and ($mapq =~ tr/,/,/)>1) { $mapq=mapqual(@gff[8,0,3,4]); }
          
          $gff[2]= "ncRNA" if($reclass =~ /noncode/); # and ($gtype eq "mRNA")
          
	        ## shorten this: napct=96p,1517/1580,1472 ; to namealn=96p ?
	        if(1) { # drop most ann
	          my @keep= grep/=/, map{ my($kv)= $gannadd =~ m/\b($_=[^;\n]+)/;  $kv; } @keepmrnan;
            my @add = grep/=/, map{ my $v= $anntab{$nd}{$_}; 
                if($_ eq "napct" and $v) { $v=~s/,.*//; "namepct=$v"; }
                else { ($v and $v ne "na")?"$_=$v":""; } } @addann;
            my $newann= "$IDkey=$nd;" . join ";", 
              map{ s/;/,/g; s/\%/p/g; s/,$//; s/Target=/trg=/; $_; } 
              ("evgclass=$reclass", @keep, @add,  "oid=$oid");    
            $gff[8]= $newann; # s/\t$gann/\t$newann/;
	        } 
	        if(0) { # old, keepann
            my $addann = join "", 
              map { my $v=$anntab{$nd}{$_}; $v=~s/;/,/g; 
                if($_ eq "namealn" and $v=~/,\d/) { my($rd,$nv)=split",",$v,2; $v="$nv;Dbxref=$rd"; }
                ($v and $v ne "na")?"$_=$v;":"" } @gffann;
            $addann =~ s/name=/Name=/;
            $addann .= ";mapqual=$mapq";
            s/;($dropann)=[^;\n]+//g; #? drop existing Name, other?
            s/$/;$addann/;
            s/Target=/trg=/; # do above?
            unless( s/;oid=[^;\s]+/;oid=$oid/ ) { s/$/;oid=$oid/; } 
            s/;;/;/g;
          }
        } else {
          # strip excess exon/CDS tags: ix, nx, xxx ... keep Target/trg, what else?
          my @keep= grep/=/, map{ 
            my($kv)= $gann =~ m/\b($_=[^;\n]+)/;  
            if($kv=~/error=/){ my($gp)= $kv=~m/notfixed=([-\d]+)/; $kv="aligngap=".abs($gp); }
            elsif($kv=~/mgap=(\d+).(\d+)/){ $kv="aligngap=".abs($2-$1); }
            $kv; } @keepexonn;
          ## my @add = grep/=/, map{ my $v= $anntab{$nd}{$_}; ($v and $v ne "na")?"$_=$v":"" } @addann;
          my $newann= "$IDkey=$nd;" . join";", 
             map{ s/;/,/g; s/\%/p/g; s/Target=/trg=/; $_; } ( @keep );    
          $gff[8]=$newann; # s/\t$gann/\t$newann/;
          if($gtype eq "CDS" and $reclass =~ /noncode/){ $commout="#no."; } #?? $ok=0;
        }
        
        my $gsrcnew= (defined $GFFSRC) ? $GFFSRC : $mapclass;  $gsrcnew=~s/\W.*//; 
        $gff[1]= $gsrcnew||$gsrc; # s/\t$gsrc/\t$gsrcnew/ if($gsrcnew);
        print $commout . join("\t",@gff)."\n" if($ok);  $no++;
       
        # s/$IDkey=$IDval/$IDkey=$nd/;  # ID|Parent=newid .. done above?
        # print if($ok);  $no++;
      }
      
    } else {
      print if(/^##gff/ or (/^#/ and /evigene/)); #  gff-version, other comm?
    }
  } close($inf);
  warn "#gff output n=$no\n";
  return($no);
}

sub mapqual { ## have mapqual in genescore table now
  my($gat,$chr,$cb,$ce)=@_;
  my($cov)= $gat=~/;cov=([\d+\.]+)/;
  my($pid)= $gat=~/;pid=([\d+\.]+)/;
  my($nx)=  $gat=~/;nexon=(\d+)/;
  my $npath= "";
  if($gat =~ m/;chim[12]=([^;\s]+)/) {  # chim1=chr1:3760586-3764044:-
     $npath=$1; 
     my($nc,$nb,$ne)= $npath=~m/([\w\.]+):(\d+).(\d+)/;
     my $ci=($nc ne $chr)?3:($ne > $cb - 1000 and $nb < $ce + 1000)?1:2;
     $npath="C$ci:$npath";
  }
  ($cov,$pid,$nx)= map{ int(0.5 + $_) } ($cov,$pid,$nx);
  my $mapq= join",", $cov."c", $pid."i", $nx."x"; # use locus above, not gloc?
  $mapq.= ",Spl:$npath" if($npath=~/^C/); # add other paths?
  return $mapq;
}

sub reClassify {
  # change tag names to sensible Class : dpx7n7_ar,_frag,_loq,.. 
  my($cla)=@_;
  $cla=~s/^(dpx7\w+|vel\w+)_//; #others?
  $cla=~s/^ar/alt/; 
  $cla=~s/^nca/noncodealt/; # noncodealt or altnc ?
  $cla=~s/^mr/main/; 
  $cla=~s/^loq/lowqual/; 
  $cla=~s/^nc/noncode/; 
  $cla=~s/^te/transposon/;
  $cla=~s/^frag/fragment/; 
  # dropdup, newfrom..
  return $cla;
}


=item genescore Class types

76511 _alt  # includes _mr and _nc alts
28310 _mr   # includes 1800 _loq main cds
 1593 _nc   # keep with _mr,_alt as 'ncRNA'
 1333 _te   # keep for checks, segregate?
 2117 _frag # drop or segregate
 1843 dropdup # drop
    2 newfromaug  # fixme
  total valid loci: 29,900 (_mr,_nc)

26317 dpx7n7_mr  76197 dpx7n7_ar
1870 dpx7n7_loq
 638 dpx7n7_nc      50 dpx7n7_nca
 123 vel10in_mr     28 vel10in_ar
 955 vel10in_nc    236 vel10in_nca
1333 dpx7n7_te    2117 dpx7n7_frag   : keep segregated?
1843 dropdup      : drop these from pub files
   2 newfromaug  # fixme, others?
  **? still missing some nomap/badmap sequence loci?
  
=item genescore eg  

nDaplx7b3EVm000021t2	obID=Daplx7b3EVm000021t3	
  Class=dpx7n7_ar	Dbxref=daphmag:Dapma7bEVm000593t5	Name=Mucin 5AC, oligomeric mucus/gel-forming	
  aalen=4497,91%,complete	hobest=dapmag,	hokaks=dapmag,dapsim,dapgal,	
  homolog=daphmag:Dapma7bEVm000593t5/4595a,drosmel:NP_524060.2/3708a,tribcas:XP_015839051.1/3994a,	
  intron=59/61	intronerr=-1,antisense	
  maploc=scaffold_45:576525-609532:+	mapqual=97a,99i,14787l,62x	
  napct=98%,4595/4688,4497	# prior namealn= .. which?
  oid=Daplx7b3EVm000021t3,Daplx6cgEVm000019t5,Daplx6mlEVm000027t1,daplx6ml10dn9anvelvk53Loc542t7

=item genescore keys

111709 Class  == mapclass of arath
111709 obID
111707 oid
111599 aalen   #? aaqual= or aalen= 
   # want cdsoff also? ; mapqual has trlen as '999l', revert to trlen=999 ?
111599 maploc > locus
111599 mapqual
108803 intron   # found missing, should be ok now
98830 homolog
90543 hokaks
19412 hobest      # useful or not?
83778 Dbxref
83778 Name
83778 napct
10932 intronerr
4105 rnatype    # mostly noncode, some code, add rnatype=code to all mRNA ?
1343 iste
  26 updid

//old//
111714 Class
111606 aalen  # want cdsoff also? (mapqual has trlen) 
..

=cut

sub putGenescore { # write_newanntab()
  my($outf)= @_;  warn "#putGenescore $outf\n" if($debug);
  if(-f $outf){ rename($outf,"$outf.old"); }
  open(OUTF,'>',$outf) or die $outf; my $nmap=0;

  ## .. table w/ headers, no key=val, or key=val format ?
  # my @hdr= grep{ not m/^oid/ } sort keys %anhdr;
  # print join("\t","PublicID","OrigID",map{ ucfirst $_ } @hdr)."\n";
  
  #  for my $pd (sort _sortid keys %pubids)
  for my $pd (sort keys %pubids) 
  {
    # my $td= $dop{$pd}||$pd; # oid1 of pubid; check all oids if missing annot val?
    ## handle also newid{pd}, some Class=drop should be dropped here?
    my $pdo= $oldpubid{$pd} || $pd;
    my $td= $dop{$pdo} || $pd; 

    my %ank= map{ $_ => 1 } (keys %{$annot{$pdo}}, keys %{$annot{$td}});
    my @ank= sort keys %ank;
    my $obid= $td; $obid.=",$pdo" if($pdo ne $pd);
    print OUTF join("\t",$pd,"obID=$obid");  # add pdo to obID ?
    for my $k (@ank) { 
      my $v= $annot{$pdo}{$k} || $annot{$td}{$k};
      print OUTF "\t$k=$v" if(defined $v);
    } print OUTF "\n";
  }
  close(OUTF);
  return $nmap;
}

sub read_newanntab {
  my($forwhat,$inf)=@_;
  my ($nan,$nskip,$nrep)=(0) x 9;
  open(F,$inf) or die $inf;   # stdin? or filename
  while(<F>) { next unless(/^\w/);
    chomp; my @v=split"\t"; 
    
    my $nd= shift @v; #no ID=xxx; not $v[0];
    my %ann= map{ my($k,$v)=split"=",$_,2; $k=>$v; } grep/=/,@v;
    my $obd= $ann{'obID'}; $obd=~s/,.*//; # changed pubid, keep old?
    my $oids=$ann{'oid'};

    $nd=~s/_[CG]\d+$// if($forwhat =~ /seq/); # have any _CG tags in this ann table??
    
    if($nod{$nd}) {
      my $replace=0;
#       my $imc=$anncols{'Class'}||0; # CHANGE
#       $replace=1 if(($noclass and $anntab{$nd}{'Class'} =~ /$noclass/) or ($oclass and $v[$imc] =~ /$oclass/));
#       $replace=0 if(($noclass and $v[$imc] =~ /$noclass/) or ($oclass and $anntab{$nd}{'Class'} =~ /$noclass/));
      if($replace) { delete $anntab{$nd}; delete $nod{$nd}; $nrep++; }
      else { $nskip++; next; }
    }    

    my @oids= split",",$oids;
    my $tda= $obd || $oids[0]; # $tda=$oids[0];
    # my $tdv= $obd; # ($tdv)= grep /Arath5EV/, @oids;
    # $tdv=$tda unless($tdv);
    for my $td ($tda) { # was ($tda,$tdv)
      $tod{$td}= $nd unless($tod{$td});  # nd dup _C1,C2,G2,..
      $nod{$nd}= $td; 
    }
    # for my $i (1..$#v) { $anntab{$nd}{$hd[$i]}=$v[$i]; } 
    for my $k (keys %ann) { $anntab{$nd}{$k}= $ann{$k}; $anncols{$k}++; }
    $nan++;
    
    #------- arath version .. headers, not key=val
    # if(/^Public/){ 
    #   @hd= map{ lc($_) } @v; for $i (1..$#hd) { $anncols{$hd[$i]}=$i; } 
    # } else { 
    #   my $nd=$v[0]; my $oids=$v[1]; 
    #   
    #   my $ndin=$nd; 
    #   $nd=~s/_[CG]\d+$// if($forwhat =~ /seq/); # unless($forwhat =~ /gff/); # ..$ENV{gff} # for seq out,not gff??
    #   if($nod{$nd}) {
    #     my $replace=0;
    #     my $imc=$anncols{mapclass}||0;
    #     $replace=1 if(($noclass and $anntab{$nd}{mapclass} =~ /$noclass/) or ($oclass and $v[$imc] =~ /$oclass/));
    #     $replace=0 if(($noclass and $v[$imc] =~ /$noclass/) or ($oclass and $anntab{$nd}{mapclass} =~ /$noclass/));
    #     if($replace) { delete $anntab{$nd}; delete $nod{$nd}; $nrep++; }
    #     else { $nskip++; next; }
    #   }
    # 
    #   my @oids= split",",$oids;
    #   $tda=$oids[0];
    #   ($tdv)= grep /Arath5EV/, @oids;
    #   $tdv=$tda unless($tdv);
    #   for my $td ($tda,$tdv) { 
    #     $tod{$td}= $nd unless($tod{$td});  # nd dup _C1,C2,G2,..
    #     $nod{$nd}= $td; 
    #   }
    #   for my $i (1..$#v) { $anntab{$nd}{$hd[$i]}=$v[$i]; } 
    #   $nan++;
    # }
    
  }
  warn "#read_newanntab=$nan, replace=$nrep, skip=$nskip\n"; 
}



sub readNames {
  my($inf)= @_;  warn "#readNames $inf\n" if($debug);
  open(F,$inf) or die $inf; my $nmap=0;
  while(<F>){  next if(/^\W/); chomp; my @v=split"\t"; 
    my($td,$name,$napct,$nadbx)=@v; 
    next if($name eq "noname");
    $annot{$td}{Name}=$name; # replacing mrna.attr
    $annot{$td}{napct}= $napct;
    $annot{$td}{Dbxref}= cleanRefspp($nadbx);
  } close(F);
  return($nmap); # or %annot
}

# best3gs_alliup_hoval.tab = 4 ref spp columns: Dapma, drosmel, daphother, insectother
# Daplx7b3EVm000001t1	12593a,62%l,dapmaevg14:Dapma7bEVm009447t1	3019a,168%l,drosmel16nc:NP_001261019.1	same	same
# Daplx7b3EVm000001t10	3328a,15%l,dapmaevg14:Dapma7bEVm009447t2	893a,176%l,drosmel16nc:NP_725511.2	same	933a,182%l,bemtab:bem6nc109032496t2
sub cleanRefspp {
  my ($ho)= @_;
  $ho =~ s/drosmel16nc:/drosmel:/;
  $ho =~ s/tribcas16nc:/tribcas:/;
  $ho =~ s/dapmaevg14:/daphmag:/;
  $ho =~ s/dapgal6tsa:/daphgal:/;
  $ho =~ s/dapsim:/daphsim:/;
  # bemtab: ok, other?
  return($ho);      
}

sub readHoval {
  my($inf)= @_;  warn "#readHoval $inf\n" if($debug);
  open(F,$inf) or die $inf; my $nmap=0;
  while(<F>){  next if(/^\W/); chomp; my @v=split"\t"; 
    # my($td,$hodmag,$hodmel,$hodaph,$hoinsect)=@v; 
    my $td=shift @v; my $oldho="";
    if($annot{$td}{homolog} and grep(/:/,@v)) { $oldho= delete $annot{$td}{homolog}; } # from mrna.attr .. drop?
    for my $ho (@v) { 
      next if($ho =~ /^(na|same)/); # next unless($ho =~ /^\d/);
      $ho= cleanRefspp($ho);
      # $ho =~ s/,\d+%l,/,/; # 893a,180%l,drosmel16nc:NP_725511.2 > 893a,drosmel:NP_725511.2 ? or drosmel:NP_725511.2/893a
      my($ha,$hl,$hd)=split",",$ho; $ho="$hd/$ha";
      $annot{$td}{homolog}.= "$ho," unless($ho =~ /^same/); 
      }
  } close(F);
  return($nmap); # or %annot
}

sub readReclasstab {
  my($inf)= @_;  warn "#readReclasstab $inf\n" if($debug);
  open(F,$inf) or die $inf; my $nmap=0;
  while(<F>){  next if(/^\W/); chomp; my @v=split"\t";
    # mod of dpx7n7pub_cullnc.masource.uptab, col2 has act: drop,reloc,replace,..
    # dpx7n7pub_updbyhand.tab 
    # drop,jointe # reloc,at1 # reloc,at2

    my($pd,$bd,$cd,$aaw) = @v;
    if($bd =~ /^drop/) { 
      $bd=~s/\W/_/g; $annot{$pd}{Class}=$bd; #? or remove from pubids{}
      # my $oldi= delete $pubids{$pd};
      my $newid="drop_".$pd; 
      $oldpubid{$newid}=$pd;   my $i= delete $pubids{$pd};

    } elsif($bd =~ /^replace/) { my($bc,$bv)=split",",$bd;  
       ## replace pd annots w/ new data? need replace gff,seq records also, new ID?

    } elsif($bd =~ /^reloc/) { my($bc,$bv)=split",",$bd;  
       ## bv = new pubid template  '[abcd]t[123]' .. other forms?
       ## other: nDaplx7b3EVm011223t999	reloc,at1 newfromaug	gff:scaffold_202 .. =AUGcex2s202g98t1;
       ## out: nDaplx7b3EVm011223at1	obID=nDaplx7b3EVm011223at1,nDaplx7b3EVm011223t999	updid=nDaplx7b3EVm011223at1/nDaplx7b3EVm011223t999
       ## ^^ need more annots for "newfrom.."

       my $newid=$pd; $newid=~s/t\d+$/$bv/;
       $oldpubid{$newid}=$pd;  #? dont add old pd to oids..
       my $i= delete $pubids{$pd}; $pubids{$newid}= $i;
       $annot{$pd}{updid}="$newid/$pd"; # idchange/reid/newid/updid/..
       my $oclass= $annot{$pd}{Class};
       if($cd =~ /^new/ and not $oclass) { $annot{$pd}{Class}= $cd; } ## more annots
       elsif($bv=~/t1$/ and $oclass=~/_ar/) { $oclass=~s/_ar/_mr/; $annot{$pd}{Class}= $oclass; }
     }

  } close(F);
  return($nmap); # or %annot
}

sub readIdtable {
  my($inf)= @_;  warn "#readIdtable $inf\n" if($debug);
  open(F,$inf) or die $inf; my $nmap=0;
  while(<F>){  next if(/^\W/); chomp; my @v=split"\t"; 
    my($pd,$bd,$cd,$aaw) = @v;

    ## there are some dup rows in id tables, seem to have same data.. not quite, 3 are dropdup in 2nd
# < nDaplx7b3EVm007612t6	obID=Daplx7b3EVm021726t2	Class=dropdup
# > nDaplx7b3EVm007612t6	obID=Daplx7b3EVm021726t2	Class=dpx7n7_ar
# < nDaplx7b3EVm012442t2	obID=Daplx7b3inewEVm04204t1	Class=dropdup
# < nDaplx7b3EVm024420t1	obID=Daplx7b3EVm024420t1	Class=dropdup

    if($pubids{$pd} and not($bd =~ /^drop/)) { next; } # all 3 drop dup rows are valid
    $pubids{$pd}= ++$nmap;
    my($oid)= grep/oid=/,@v;
    my $obid=0;

    ## valid aaqual is in these tables?? col4; also ,mapaa tag 
    if($aaw=~m/^\d+,\d+%,\w/){ $annot{$pd}{aalen}= $aaw; } #? aaqual= or aalen=

    # my $CLASSTAG= 'Class'; #genetype= or Class= ?    
    # fixme alts of vel10in_mr,_nc new class tag
    my $cla="noclass";
    if($cd=~/^vel10in_(\w+)/){ 
      $cla=$cd;  # $annot{$pd}{Class}=$cd;
      my $rt=$1; $annot{$pd}{rnatype}=($rt=~/nc/)?"noncode":"code";  # rnatype=noncode/code
    } elsif($bd=~/^(dpx7|drop)/) { 
      $bd=~s/.okdup//; # bd includes dropdup
      $bd=~s/(drop\w*)\..*/$1/; #  dropdup
      $cla=$bd; # $annot{$pd}{Class}=$bd;
    }
    if($cla =~ /_(mr|nc)/ and $pd=~/t[2-9]$/) { $cla =~ s/_mr/_ar/; $cla=~s/_nc/_nca/; }
    $annot{$pd}{Class}=$cla;
    # Class/genetype =  dpx7n7_ar _frag _loq _mr _nc _te , dropdup, vel10in_nc/mr

    unless($oid){ if($bd=~/^\w/ and $bd=~/t\d+/) { $oid=$bd; } } # for vel10in_ table
    if($oid){      
      my $i=0; $oid=~s/oid=//; # pod == pubid of oids; dop == oid1 of pubid
      for my $d (split",",$oid) { 
        (my $dc=$d)=~s/_[CG]\d+$//;
        $pod{$d}= $pod{$dc}= $pd; 
        if($i++ < 1){ $dop{$pd}=$d; $obid=$d; } 
        }
      my $ooid= $annot{$pd}{oid} || $annot{$obid}{oid} || 0;
      if(length($oid) > length($ooid)) { $annot{$pd}{oid}= $oid; } # or count tr/,/,/
    }
    # ..
  } close(F);
  return($nmap); # or %annot
}

=item idtabs
==> dpx7n7pub_newing.pubids <==
nDaplx7b3EVm064001t1	daplx6ml10dn9anvelvk31Loc1t350	vel10in_nc	35,50%,partial5	noho	98a,100i,212l,2x,1in,anti	scaffold_68:420605:420905:+
nDaplx7b3EVm064001t2	daplx6ml10dn9anvelvk31Loc1t462	vel10in_nc	35,50%,partial5	noho	98a,100i,212l,2x,1in	scaffold_68:420605:420905:+
nDaplx7b3EVm064002t1	daplx6ml10dn9anvelvk31Loc1t1200	vel10in_nc	80,23%,complete-utrbad	noho	100a,100i,1038l,4x,3in,anti	scaffold_381:39195:40606:+

==> dpx7n7pub_cullnc.masource.uptab <==
nDaplx7b3EVm000018t2	dpx7n7_ar	beo3dpxevg5a	5051,91%,complete	aaref:100%,5100/5033,5051,dapmaevg14:Dapma7bEVm004989t1	scaffold_296:115195:132549:-	oid=Daplx7b3EVm000018t2,Daplx6mlEVm000015t1,Daplx6cgEVm000017t1,Daplx6mlEVm000015t1,daplx6ml25mrn9alvelvk35Loc3339t2	best3gs_alliup.aa	dupid:nDaplx7b3EVm000018t2/Daplx6mlEVm000015t1
nDaplx7b3EVm000018t3	dropdup.dpx7n7_ar	beo3alt	5051,91%,complete	aaref:100%,5100/5033,5051,dapmaevg14:Dapma7bEVm004989t1	scaffold_296:115195:132549:-	oid=Daplx7b3EVm000018t2,Daplx6cgEVm000017t5,Daplx6mlEVm000015t1,daplx6ml25mrn9alvelvk35Loc3339t2	best3gs_alliup.aa	dupid:nDaplx7b3EVm000018t2/Daplx6mlEVm000015t1
nDaplx7b3EVm000023t1	dpx7n7_mr	beo3dpxevg5a	4755,78%,complete	aaref:100%,4765/4757,4755,dapmaevg14:Dapma7bEVm019116t1	scaffold_6:715470:735935:+	oid=Daplx7b3EVm000023t2,Daplx6mlEVm000018t1,Daplx6cgEVm000021t1,Daplx6mlEVm000018t1,daplx6ml10dn9anvelvk31Loc2177t10	best3gs_alliup.aa	dupid:nDaplx7b3EVm000023t1/Daplx6mlEVm000018t1

=cut

sub readMrnaAttr { # best3gs_alliup_mrna.attr 
  my($inf)= @_;  warn "#readMrnaAttr $inf\n" if($debug);
  open(F,$inf) or die $inf; my $nmap=0;
  while(<F>){  next if(/^\W/); chomp; my @v=split"\t"; 
    my($td,$mapid,$src,$loc)=@v; 
    my($oid)=grep/^oid=/,@v; 
    ## add? aalen= , ann key as aaqual= or aalen=
    my @at=grep /(Name|homolog|aalen)=/, @v; # is this homolog ok? hoval may be better
    if($oid){ $oid =~ s/=/=$mapid,/ unless($oid=~m/$mapid/); } else { $oid="oid=$mapid"; }
    map{ 
      my($k,$av)= m/^(\w+)=(.+)/; 
      if($k eq "homolog" and $av=~/^\d/) {  # homolog=4916,dapmaevg14:xxx
       my($ha,$hd)=split",",$av; $hd= cleanRefspp($hd); $av="$hd/$ha,"; } 
      $annot{$td}{$k}=$av; } ($oid,@at); 
  } close(F);
  return($nmap); # or %annot
}

sub readTable { # dpx7n7pub.noncode.rereclass.tab
  my($inf, $key)= @_;  warn "#readTable(key $key): $inf\n" if($debug);
  open(F,$inf) or die $inf; my $nmap=0;
  while(<F>){  next if(/^\W/); chomp; my @v=split"\t"; 
    my($td,$bv,$cv)=@v;
    if($bv=~/^(Noncode|Code|Unknown)/){  $annot{$td}{rnatype}= lc($bv); } # rnatype/class?
    elsif($bv=~m/iste=(\d+)/){  $annot{$td}{iste}= $1; }  
    elsif($bv=~m/^ints=(\d+\S+)/){ 
      my $ints=$1; $ints=~s/^\d+,//; $annot{$td}{intron}= $ints; 
      if($cv=~m/inerr:(\S+)/){  my $ine=$1; # fix shorthand a/o/l only antisense/longerr?
       # nDaplx7b3EVm000001t1	ints=94,16/17	inerr:-1,aol:1/0/0
       # nDaplx7b3EVm000001t14	ints=90,19/21	inerr:-1,aol:0/0/1,longerr:15218
        if( my($ae,$oe,$le)= $ine=~m/aol:(\d+).(\d+).(\d+)/ ) {
          $ine=~s/,aol:.*//; 
          $ine.=",antisense" if($ae>0); 
          $ine.=",longerr" if($le>0);
          }
        $annot{$td}{intronerr}=$ine; 
        } 
      }  
    elsif($bv=~m/^hospp=(\w+\S+)/){ # $annot{$td}{homols}= $1; 
      foreach my $kv (@v) { 
        my($k,$vv)=split"=",$kv; $k=~s/^kas/hokaks/; 
        my $hashoval= $annot{$td}{homolog};
        $annot{$td}{$k}=$vv unless($hashoval and $k=~/hoval|hospp/); # skip hospp also? should be in homolog=
        } # hospp,hoval,hobest,kaks
      }  
    $nmap++;
  } close(F);
  return($nmap); # or %ann
}

sub readAligntab { # xxx.align.tab from gene.gff
  my($inf)= @_; warn "#readAligntab $inf\n" if($debug);
  open(F,$inf) or die $inf; my $nmap=0;
  while(<F>){  next if(/^\W/); chomp; my @v=split"\t";
    # my($td,$cov,$pi,$nx,$mloc,$np,$oid,$tag,$clen,$cdspan,$aalen)= @v; # map.attr table
    # align.tab
    # GenomeID	gespan	geor	AQueryID	quspan	match	qlen	cov	pid	path	indels	nexon	splice	aalen	offs	aamap	sense	oid	tag	cdspan

    my($mr,$mbe,$mo,
      $td,$qusp,$mat,$clen,$cov,$pid,$np,$indel,$nx,$spli,
      $aalen,$cdoff,$aamap,$sense,$oid,$tag,$cdspan)= @v;  

    my $issplit= ($np=~/C[123]:/)?1:0; # cdspan vs unmapped aalen, Splits and other effects are problems
    # make mapqual tuple: 99a,98i,1234l,4x,[Spl:] ; eg 98a,100i,212l,2x,1in ** change 'a' to 'c' for prior key
    my $mloc= join":",$mr,$mbe,$mo;
    my $mapq= join(",", int($cov)."c", int($pid)."i", $clen."l", int($nx)."x");
    # arath: join",", $cov."c", $pid."i", $nx."x";
    if($issplit) { my($ps)=$np=~m/,(\d+)$/; $ps||=9; $mapq.=",Spl:$ps%"; }    
    $annot{$td}{'locus'}= $mloc;
    $annot{$td}{'mapqual'}= $mapq;
    $annot{$td}{'aalen'}= $aalen; #? aaqual= or aalen=
    
    # $mapattr{$td}= \@v;
    # my($mr,$mb,$me,$mo)= split/[:-]/,$mloc,4; #? ..
    # my($mb,$me)= split/[:-]/,$mbe,2;
    # $maploc{$td}{mr}=$mr; $maploc{$td}{mo}=$mo; $maploc{$td}{mb}=$mb; $maploc{$td}{me}=$me; 
    # my($cb,$ce)= split /\W+/,$cdspan; # 99-199 or 99..199
    # $maploc{$td}{cb}= $cb; $maploc{$td}{ce}=$ce; # $cspan; # or split?

    # my $poormap= ($cov < $MINCOVMAP)?1:0; # separate out split val
    # my($aaw)= ($aalen=~m/^(\d+)/)?$1:0;
    # $maploc{$td}{aalen}= $aaw; #? or from mapattr{}  
    # $maploc{$td}{issplit}= $issplit;   # keep this, update w/ eqgene id_C tag
    # $maploc{$td}{poormap}= $poormap;    # separate out split val
    
    $nmap++;
  } close(F);
  return $nmap;
}

# sub readIsTE { #  nDaplx7b3EVm016125t1	drop,iste=100	obID=Daplx7b3EVm016125t1,
#   my($inf)= @_;  warn "#readNoncodeclass $inf\n" if($debug);
#   open(F,$inf) or die $inf; my $nmap=0;
#   while(<F>){  next if(/^\W/); chomp; my @v=split"\t"; 
#     if($v[1]=~m/iste=(\d+)/){  $annot{$td}{iste}= $1; }  
#   } close(F);
#   return($nmap); # or %ann
# }
# 
# sub readIntronsum { # dpx7n7pub.intronsum.tab
#   my($inf)= @_;  warn "#readIntronsum $inf\n" if($debug);
#   open(F,$inf) or die $inf; my $nmap=0;
#   while(<F>){  next if(/^\W/); chomp; my @v=split"\t"; 
#   # nDaplx7b3EVm000001t13	ints=100,21/21	
#   # nDaplx7b3EVm000001t14	ints=90,19/21	inerr:-1,aol:0/0/1,longerr:15218
#     my($td,$inv,$inerr)=@v;
#     $inv=~s/^\w+=//; $annot{$td}{inval}=$inv;
#     $inerr=~s/^\w+=//; $annot{$td}{inerr}=$inerr if($inerr);
#   } close(F);
#   return($nmap);
# }

__END__

=item dpxevg17 gff out

** MISSING vel10in, ID=udaplx buggers
  also, change src=vel10in_xx to dpx7n7_xx
  scaffold_68	vel10in_nc	mRNA	..	ID=udaplx6ml10dn9anvelvk31Loc1t350;trg=daplx6ml10dn9anvelvk31Loc1t350 1 207;aalen=34;cov=97.6;match=207;nexon=2;pid=100.0;qlen=212;sense=-1

gff ID counts (mRNA|ncRNA):
  uniq id: 109822; total id: 109968 ; 146 dups, vs 230 Split annots, ie. not all Split are valid
    Split= : 103 uniq id, 230 total id, 127 extra parts, 19 unaccounted for..
    Nosplit dup ids (19)
nDaplx7b3EVm019113t4
nDaplx7b3EVm025868t2
nDaplx7b3EVm029163t1
nDaplx7b3EVm060962t3
nDaplx7b3EVm061369t1
nDaplx7b3EVm061542t1
nDaplx7b3EVm062229t3
nDaplx7b3EVm062851t1
nDaplx7b3EVm062873t1
nDaplx7b3EVm062975t1
nDaplx7b3EVm063214t1
nDaplx7b3EVm063238t1
nDaplx7b3EVm063257t1
nDaplx7b3EVm063272t1
nDaplx7b3EVm063301t1
nDaplx7b3EVm063380t1
nDaplx7b3EVm063478t1
nDaplx7b3EVm063516t1
nDaplx7b3EVm063650t1

genescore entries with oid sharing 2+ pubids, may be diff map loci
.. hopefully most are already checked as valid diff map loci (or alt dups should be fixed, in gff? )
  Daplx5cEVm005536t18	nDaplx7b3EVm027780t1	nDaplx7b3EVm062279t1
  Daplx5cEVm010341t10	nDaplx7b3EVm012809t1	nDaplx7b3EVm012809t7
  Daplx5cEVm011874t2	nDaplx7b3EVm059216t2	nDaplx7b3EVm059232t3
  Daplx5cEVm013283t12	nDaplx7b3EVm025311t2	nDaplx7b3EVm059347t3
  Daplx5cEVm016381t5	nDaplx7b3EVm027766t1	nDaplx7b3EVm063096t1
  Daplx5cEVm016958t1	nDaplx7b3EVm019214t1	nDaplx7b3EVm060499t1
    .... n=397 (excluding dropdup), n=287 are velml10d/dpx7n7pub_newing
    
    
nDaplx7b3EVm019113t4 = C1: split type, 2 overlapping parts
  p1: cov=40.6;nexon=1;pid=100.0;path=1/2;chim2=sc76:144656-152949:+;unsplit=1;
  p2: cov=59.0;nexon=2;pid=98.7;path=2/2;chim1=sc76:151543-152011:.;
nDaplx7b3EVm029163t1 = TE, C3: split type, 2 scaffolds  Pol Polyprotein
  p1: cov=58.3;indels=3/1;nexon=2;pid=98.8;path=2/5;chim1=sc392:54959-56828:.;unsplit=2
  p2: cov=40.9;nexon=1;pid=99.0;path=1/5;chim2=sc131:5501-8413:+;unsplit=1
nDaplx7b3EVm061369t1 = data sutter, same locus 
  has 1 entry in dpx7n7pub_cullnc.gff, ID=nDaplx7b3EVm061369t1; msrc=velml10d
  has 1 entry in dpx7n7pub_newing.gff, ID=udaplx6ml10dn9anvelvk53Loc1887t3
  2 entries output to dpx7fin.gff
nDaplx7b3EVm063238t1 = data stutter, same locus, like above, entry in each input gff

** FIXME: dpx7n7pub_upd.genescore has 2+ rows, diff pubIDs for same oid=daplx6ml10dn9anvelvk77Loc19489t1, others
  .. same obID=daplx6ml10dn9anvelvk77Loc19489t1 for some, same oid=
  .. from dpx7n7_xxx and vel10in_nc
  ** BUT also have cases of same seq oid, diff map locus, need to distinguish,
  also pick better of two "dups", e.g. mapqual diffs
  grep daplx6ml10dn9anvelvk77Loc19489t1 dpx7n7pub_upd.genescore
poor>nDaplx7b3EVm063238t1	obID=daplx6ml10dn9anvelvk77Loc19489t1	Class=dpx7n7_loq	aalen=44,56%,complete	
  intron=1/1	locus=scaffold_120:244445-244727:-	mapqual=87c,100i,238l,2x	oid=daplx6ml10dn9anvelvk77Loc19489t1
best>nDaplx7b3EVm064964t1	obID=daplx6ml10dn9anvelvk77Loc19489t1	Class=vel10in_nc	aalen=44,56%,complete	
  intron=2/2	locus=scaffold_120:244445-244880:-	mapqual=100c,100i,238l,3x	oid=daplx6ml10dn9anvelvk77Loc19489t1	rnatype=noncode
  
  
  
attr mRNA/ncRNA
try2:
109968 ID
109968 evgclass
109968 oid
109951 msrc
109860 aalen
109860 mapqual
107116 intron
97371 homolog
95801 trg
89512 hokaks
82290 Name
82290 namepct
81953 Dbxref
10765 intronerr
 230 Split
 
try1:
109968 ID
109968 evgclass
109968 oid
109860 aalen
109860 mapqual
108627 msrc
107116 intron
97371 homolog
95801 trg
89512 hokaks
82290 Name
82290 namepct
81953 Dbxref
27678 namepctNUL << what?
10765 intronerr
  230 Split

types, try3: counts Split parts
  -- vel10in_ > dpx7n7_ 
  26457 dpx7n7_mr	mRNA
  76314 dpx7n7_ar	mRNA
   1883 dpx7n7_loq	mRNA
   1578 dpx7n7_nc	ncRNA
    284 dpx7n7_nca	ncRNA
   1335 dpx7n7_te	mRNA
   2117 dpx7n7_frag	mRNA

try2:
  76286 dpx7n7_ar	mRNA
  26335 dpx7n7_mr	mRNA
   638 dpx7n7_nc	ncRNA
    50 dpx7n7_nca	ncRNA
  1883 dpx7n7_loq	mRNA
  1335 dpx7n7_te	mRNA
  2117 dpx7n7_frag	mRNA
    28 vel10in_ar	mRNA
   122 vel10in_mr	mRNA
   940 vel10in_nc	ncRNA
   234 vel10in_nca	ncRNA

  644466 dpx7n7_ar	CDS
  717449 dpx7n7_ar	exon
  4432 dpx7n7_frag	CDS
  5783 dpx7n7_frag	exon
  3629 dpx7n7_loq	CDS
  4342 dpx7n7_loq	exon
  149867 dpx7n7_mr	CDS
  162692 dpx7n7_mr	exon
  2066 dpx7n7_nc	exon
   232 dpx7n7_nca	exon
  6050 dpx7n7_te	CDS
  6372 dpx7n7_te	exon
    43 vel10in_ar	CDS
    88 vel10in_ar	exon
   153 vel10in_mr	CDS
   317 vel10in_mr	exon
  2781 vel10in_nc	exon
   817 vel10in_nca	exon

cat dpx7fin.gff | grep -v '^#' | cut -f2,3 | sort | uniq -c       
  76283 dpx7n7_ar	mRNA
  2117 dpx7n7_frag	mRNA
  1870 dpx7n7_loq	mRNA
  26334 dpx7n7_mr	mRNA
   638 dpx7n7_nc	mRNA      << ncRNA
    50 dpx7n7_nca	mRNA     << ncRNA
  1335 dpx7n7_te	mRNA


=item dpxevg17 gff attr
  grep -h '        mRNA'  dpx7n7pub_{cullnc,newing}.gff | grep -v '^#' | perl -ne '@v=split"\t"; $an=$v[8]; map{ ($k,$v)=split"=",$_,2; $k.="NUL" unless($v); print "$k\n"; } split";",$an;' | sort | uniq -c | sort -k1,1nr 
  
  mRNA:
  112798 aalen  *     5454 xtrim
  111844 ID           5408 cdsindel
  110904 nexon        5360 ovrna
  110502 msrc         5243 path
  110502 obID         4976 cdsfix
  106439 oid          4362 inmissadd
  98296 NUL == ;;     3819 utrx
  97666 pid           3791 cututrchimer
  97666 trg           3139 flipstrand   : maybe keep
  97663 cov           3109 oldsense     : ^^ ditto, one or other, and sense=?
  97054 cxlen         2531 inhit
  95495 offs          2531 intralign
  93736 protein       2281 sense
  93598 ocds          2125 ioldid
  93499 aaold         1389 chim2
  92927 cdsoff        1307 nNUL
  92927 cdsspan       1112 chimera
  91835 cdnaorf       943 exons
  83286 Name          717 chim1
  81393 clen          390 spliceNUL  == splice=;
  71130 qlen          372 unsplit
  67059 alttr         362 oaaln
  67059 oldaltid      352 dgcds
  64814 offold        271 dmcds
  63881 indels        232 Split
  63784 Dbxref        186 dscds
  60291 nintron       168 splicemix
  32544 cdnabest      155 strandmix
  28831 homolog       129 iste
  18038 pro1          128 must
  13457 fixoldspan    109 evg5ovcds
  13127 DbxrefNUL     24 sphase
  12000 gescore       12 nocds
  11610 splice        8 aamap
  11084 gaps          7 cutNUL
  10675 cututrx       3 joinsplit
  9677 match              
  6159 evgreclass
  -----------------------

  exon only: 
    914946 Parent
    823431 trg
    709185 ix   : maybe keep
    105053 splice == splign splice=AGGT-, drop
     6045 mgap   : splign align error (gap in exon align)
      trg=Daplx7b3EVm000048t1 2 8317;mgap=5892-6056;splice=CGTC-  >> aligngap=164
     2040 error  : gmap align err, keep these? for exon, not CDS
      trg=Daplx6suEVm010692t1 149 713;ix=3;error=ERROR.span:genome_span:428,tr_span:564,sc10:773093-772665,notfixed=-136
      ^^ rewrite to same format? aligngap=136 (from notfixed)

  exon/CDS :
    1736967 Parent
    1260268 ix
    824644 trg
    158617 
    158617 NUL
    105053 splice
    11194 mgap
    11056 cexoff
    9193 align
    8804 introna
    5989 part
    3538 error
       6 errspan
  
=item dpxevg17 seq attr

keepattr: aalen  srcID?
rewrite from annot: oid= (subset? 1st, last only?) ; cutoids='AUGcex\w+(t[2-9]+|t1\d+)|Daplx6cgEVm\w+'
rewrite : organism=Daphnia_pulex  type=protein/CDS/mRNA/ncRNA .. other?
dropattr: aaold codepot cov cxlen exons homolog inewid length locus mapq nexon read_count sense src strand
   asrc?? = beo3alt beo3dpxevg5a .. etc, drop
srcFn/ID/OK : useful or not? keep srcID not other srcNN
  srcID=Daplx7b3EVm000118t5;srcFn=best3gs_alliup.mrna; srcOK=ok;

fixup: cdsoff/offs/   clen/cxlen/length/qlen

grep -h '^>'  *.aa | perl -ne 's/>\S+//; @v=grep/=/, split/[; ]+/; map{ ($k,$v)=split"="; print "$k\n";} @v; ' | sort | uniq -c | head -40
110700 aalen
 702 aaold
108532 asrc
 803 cdsoff
101001 clen
2411 codepot
2411 cov
 803 cxlen
 938 exons
26356 homolog
2411 inewid
 212 length
 938 locus
1342 mapq
15497 nexon
1342 offs
109874 oid
82587 organism
 615 qlen
 212 read_count
1060 sense
1342 src
108532 srcFn
108532 srcID
108532 srcOK
3776 strand
108133 type

*.cds *.mrna add
8628 cdsflag
  48 error
 882 cexoff
3776 cdslen
6870 ix
8684 len
8684 loc
 938 locus
  28 part
 344 partial_gene
 615 qlen
.. more

=cut 


