#!/usr/bin/env perl
# evg5arathann.pl

my $IDPRE="Arath5nEV";
my $debug=$ENV{debug} || 1;
my $oclass=$ENV{class}||""; # keep these
my $noclass=$ENV{noclass}||""; # skip these
my $GFFSRC=(defined $ENV{gffsrc}) ? $ENV{gffsrc} : undef;

# output updates: env seq=evg5arath.aa perl evg5arathann.pl evg5arath_icnannbest.anntab > evg5arath_icn.aa

# >Arath5EVm000112t1 type=protein; aalen=1821,94%,complete; clen=5784; offs=127-5592; oid=Arath3EVm000114t1; organism=Arabidopsis_thaliana; 
# uvcut=siam_begin:0,52,1-52,end5; uvfa=CCCATTTACCGTATCCTCTA..ATTGCTCAATGAGAGATTCG; evgclass=main,okay,match:Arath4EVm000111t1,pct:100/100/.;

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

#------------------------------------------

sub _sortid { my($ag,$ai)= $a=~m/EVm(\d+)t(\d+)/;  my($bg,$bi)= $b=~m/EVm(\d+)t(\d+)/; 
  return ($ag <=> $bg or $ai <=> $bi or $a cmp $b); }

# PublicID	OrigID	Aaqual	Cdsoff	Locus	Mapclass	Mapqual	Mustkeep	Name	Namealn	Orlog_cacao	Orlog_orange	Trlen
# Arath5nEVm000001t1	Arath5EVm000001t1,Arath3EVm000001t1,arab6roo3dnta6ahvelvk65Loc1184t5	5438,98%,complete	61-16377	chr1:25069487:25095586:-	AtEvg5nmain,Arath5nEVm000001t1	100c,100i,74x	na	midasin-like protein;	AT1G67120.1,98%,16170	cacao:Thecc7EVm000001t2,47%,10672	orange:ncbig18036614t1,43%,9406	16565
sub read_newanntab {
  my($forwhat)=@_;
  my ($nan,$nskip,$nrep)=(0) x 9;
  while(<>) {
    chomp; my @v=split"\t";
    if(/^Public/){ 
      @hd= map{ lc($_) } @v; for $i (1..$#hd) { $anncols{$hd[$i]}=$i; } 
    } else { 
      my $nd=$v[0]; my $oids=$v[1]; 
      # ($td)= grep /Arath5EV/, split",",$oids; ## BUG, need other 3EV, 4EV ..

      # FIXME2: have _C1/2 in table, same seq, diff mapclass, keep good, skip bad noclass/oclass
      # if($anntab{$nd}{mapclass} =~ /$noclass/)
      # my $havend= $nod{$nd};
      
      my $ndin=$nd; 
      $nd=~s/_[CG]\d+$// if($forwhat =~ /seq/); # unless($forwhat =~ /gff/); # ..$ENV{gff} # for seq out,not gff??
      if($nod{$nd}) {
        my $replace=0;
        my $imc=$anncols{mapclass}||0;
        $replace=1 if(($noclass and $anntab{$nd}{mapclass} =~ /$noclass/) or ($oclass and $v[$imc] =~ /$oclass/));
        $replace=0 if(($noclass and $v[$imc] =~ /$noclass/) or ($oclass and $anntab{$nd}{mapclass} =~ /$noclass/));
        if($replace) { delete $anntab{$nd}; delete $nod{$nd}; $nrep++; }
        else { $nskip++; next; }
      }

      my @oids= split",",$oids;
      $tda=$oids[0];
      ($tdv)= grep /Arath5EV/, @oids;
      $tdv=$tda unless($tdv);
      for my $td ($tda,$tdv) { 
        $tod{$td}= $nd unless($tod{$td});  # nd dup _C1,C2,G2,..
        $nod{$nd}= $td; 
      }
      for my $i (1..$#v) { $anntab{$nd}{$hd[$i]}=$v[$i]; } 
      $nan++;
    }
  }
  warn "#read_newanntab=$nan, replace=$nrep, skip=$nskip\n"; 
}

sub write_newanntab {
  @hdr= grep{ not m/^oid/ } sort keys %anhdr;
  print join("\t","PublicID","OrigID",map{ ucfirst $_ } @hdr)."\n";
  for $nd (sort _sortid keys %nid) {
    $td = $nid{$nd} || "noid"; # fixme
    $oid= $anntab{$td}{oid} || "na";
    @an= map{ $nidtab{$nd}{$_} || $anntab{$td}{$_}||"na" } @hdr;
    #x @an= map{ $anntab{$td}{$_}||"na" } @hdr;
    print join("\t",$nd,"$td,$oid",@an)."\n";
  }
}

sub read_newidtab {
  $inf= $_[0] || "evg5arath_icnannbest.idtab";
  # $inf="evg5arath_icnannbest.idtab";
  open(F,$inf) or die $inf; warn "reading $inf\n";
  while(<F>){ next if(/^\W/); chomp; @v=split"\t"; 
    my($nd,$oids,$md,$cla,$loc)=@v; 
    ($td)= grep /Arath5EV/, split",",$oids; $td||="noid";
    ## new ids have gmap tags _C1/2 _G2,3,... need chop that ?
    ## Arath5nEVm000007t8_C1	Arath5EVm000006t37,Arath4EVm000006t38	Arath5nEVm000007t1	AtEvg5nalt	chr4:17067407:17068349:.
    ## Arath5nEVm000007t8_C2	noid	Arath5nEVm000006t1	AtEvg5nfrag	chr2:7773132:7791386:-
    ## Arath5nEVm000041t1_G2	noid	Arath5nEVm000041t1	AtEvg5nmain	chrchloro:143390:152918:.
    
    $nnd=$nd; $nnd=~s/_[CG]\d+//;## skip the _C/G but for locus, mapclass ? use nd or td ?
    # ($nnd,$nnt) = ($nd =~ m/^(\w+)(\_[CG]\d+)$/) ? ($1,$2) : ($nd,0);
    if($td eq "noid") {
      if($tdd= $nidx{$nnd}) { $td=$tdd; } # $tdd.$nnt
      else { $td=$nd; } #? or not
    } 
    
    $nid{$nd}= $td; $nidx{$nnd}= $td;
    #? $nidtab{$nd}=[@v]; 
    $nidtab{$nd}{locus}= $loc; 
    $nidtab{$nd}{mapclass}="$cla,$md";
    if($td ne "noid") {
    $anntab{$td}{locus}=$loc; $anhdr{locus}++;
    $anntab{$td}{mapclass}="$cla,$md"; $anhdr{mapclass}++;
    }
    
  } close(F);
}

sub mapqual {
  my($gat,$chr,$cb,$ce)=@_;
  my($cov)= $gat=~/;cov=([\d+\.]+)/;
  my($pid)= $gat=~/;pid=([\d+\.]+)/;
  my($nx)=  $gat=~/;nexon=(\d+)/;
  my $npath= "";
  if($gat =~ m/;chim[12]=([^;\s]+)/) {
     # chim1=chr1:3760586-3764044:-
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

sub putgff {
  my($gffin)=@_;
  warn "#annot of $gffin\n";

  #above# read_newanntab();
  
  my @gffann= grep{ not m/origid|locus|aaqual|cdsoff|mapqual/ } 
    sort { $anncols{$a} <=> $anncols{$b} or $a cmp $b } keys %anncols;
  #	OrigID	Aaqual	Cdsoff	Locus	Mapclass	Mapqual	Mustkeep	Name	Namealn	Orlog_cacao	Orlog_orange	Trlen

  ## FIXME: change namealn=DbxrefID,alnval to Dbxref=id;namealn=alnval
  ## FIXME: gmap2evg input.gff needs annot drops
# Target=Arath5oEVm033402t1 12 216;aamap=69;cov=54.8;match=177;nexon=2;pid=86.3;qlen=374;path=1/2;chimera=exon-exon boundary at 216;chim2=chr1:3763748-3763905:.;cdsindel=140;aalen=124,99%,partial;clen=374;offs=3-374;
# oid=Arath5EVm022860t1,Arath3EVm021214t1,tidbarab6roo3dno1sridk117Loc103323,Arath5EVm022860t1,Arath3EVm021214t1,tidbarab6roo3dno1sridk117Loc103323 locus=chr1:3763748:3763905:.; mapclass=AtEvg5nmain,Arath5oEVm033402t1; mustkeep=arath16alnbest; Name=transmembrane protein%2C putative (DUF761);; namealn=AT1G11230.2,30%,304; orlog_cacao=na; orlog_orange=na; trlen=374;

  my @dropann= qw(aamap match qlen clen chimera cdsindel); # qw( aasize offs ) ? 
  push @dropann, qw(pid chim1 chim2); # mapqual; leave cov,nexon for other uses
  my $dropann= join('|',@dropann);
    
  my ($ok,$no)=(0,0);
  open(S,$gffin) or die $gffin;
  while(<S>){
    if(/^\w/){ 
      my @gff=split"\t";
      my($tk,$td)= m/(ID|Parent)=([^;\s]+)/; # _C1,2 split and _G2,3.. dup map ids?
      my($ok,$nd,$tdc);
      ($tdc=$td)=~s/_[CG]\d+$//;

      #FIXME2: input.gff may have splits _C1,2 not in ann table, 

      #FIXME: input.gff may have newids, nd, not td of annot table; test
      if(my $d=$nod{$td}) { $nd=$td; $td=$d; }
      elsif($tdc ne $td and my $d=$nod{$tdc}) { $nd=$tdc; $td=$d; }

      ($tdc=$td)=~s/_[CG]\d+$//;
      $nd=$tod{$td}||$tod{$tdc}; 
      $ok=($nd)?1:0; 
      if($ok and $noclass) { $ok=($anntab{$nd}{mapclass} =~ /$noclass/) ? 0: 1; }
      elsif($ok and $oclass) { $ok=($anntab{$nd}{mapclass} =~ /$oclass/) ? 1: 0; }
      if($ok) {
        if(/\tmRNA/) { 
          my $oid= $anntab{$nd}{origid}||""; # check, remove dup oids..
	  my ($oidg)= (m/;oid=([^;\s]+)/)?$1:""; # both?
          # $oid= "$tdc,$oid" unless($oid=~/\b$tdc\b/);
	  my %oid; my $i=0; map{ $oid{$_}= ++$i; } split",","$tdc,$oid,$oidg"; 
          $oid=join",", sort{$oid{$a}<=>$oid{$b}} keys %oid;

	  my $mapq= mapqual(@gff[8,0,3,4]); ##$dropann.='|cov|pid|nexon|chim1|chim2';
          my $addann = join "", 
	    map { my $v=$anntab{$nd}{$_}; $v=~s/;/,/g; 
              if($_ eq "namealn" and $v=~/,\d/) { my($rd,$nv)=split",",$v,2; $v="$nv;Dbxref=$rd"; }
              ($v and $v ne "na")?"$_=$v;":"" } @gffann;
          $addann =~ s/name=/Name=/;
          $addann .= ";mapqual=$mapq";
          s/;($dropann)=[^;\n]+//g; #? drop existing Name, other?
          s/$/;$addann/;
          unless( s/;oid=[^;\s]+/;oid=$oid/ ) { s/$/;oid=$oid/; } 
	  s/;;/;/g;
        }
        my $gsrc= (defined $GFFSRC) ? $GFFSRC : $anntab{$nd}{mapclass};
        $gsrc=~s/\W.*//;
        s/\t\S+/\t$gsrc/ if($gsrc);
	s/Target=/trg=/;
        s/$tk=$td/$tk=$nd/; 
        print ;  $no++;
      }
      
    } else {
      print if(/^##gff/ or (/^#/ and /evigene/)); #  gff-version, other comm?
    }
  } close(S);
  warn "#annot output n=$no\n";
  return($no);
}


## fixme misses in to out
# ?? missing tod{od} = nd here?  >Arath3EVm032560t2 nd=Arath5nEVm033814t1
# anntab = 
# Arath5nEVm033814t1      Arath3EVm032560t2,arab6roo46aivelvk87Loc18100t1 39,22%,complete-utrbad  184-303 chr5:16714604:16715931:+        AtEvg5nmain     100c,100i,5x    intmiss na      na      na      na      532

sub putseq {
  my($seq)=@_;
  warn "#annot of $seq\n";
  ## keep: type, aalen, clen, offs, organism
  @seqann=qw( locus mapclass origid namealn name );
  @keepann=qw(type aalen clen offs organism);
  @dropann=qw(oid uvcut uvfa evgclass ); $dropann=join('|',@dropann);
  #above#$oclass=$ENV{class}||"";
  
  #above read_newanntab();
  ## split output files by mapclass: main, alt, frag
  my ($ok,$no)=(0,0);
  open(S,$seq) or die $seq;
  while(<S>){
    if(/^>(\S+)/) { 
      $td=$1; $nd=$tod{$td}; $ok=($nd)?1:0; 
      if($ok and $noclass) { $ok=($anntab{$nd}{mapclass} =~ /$noclass/) ? 0: 1; }
      elsif($ok and $oclass) { $ok=($anntab{$nd}{mapclass} =~ /$oclass/) ? 1: 0; }
      if($ok) {
        $oid=$anntab{$nd}{origid}||"noid";
        $addann = join "", map { $v=$anntab{$nd}{$_}; ($v)?" $_=$v;":"" } @seqann;
        $addann =~ s/origid=/oid=/;
        $addann =~ s/name=/Name=/;
        s/\s*($dropann)=\S+//g;
        s/$/$addann/;
        $nd=~s/_[CG]\d+$//; 
        s/^>$td/>$nd/; $no++;
      }
    }
    print if($ok);
  } close(S);
  warn "#annot output n=$no\n";
  return($no);
}

sub read_keeps {
  $inf= $_[0] || "evg5arath.keepids";
  open(F,$inf) or die $inf; warn "reading $inf\n";
  while(<F>){ next if(/^\W/); chomp; @v=split"\t"; 
    my($td,$keep,$nd,$why)=@v; $nd=~s/^Arath5\w+EV/$IDPRE/; 
    $nid{$nd}=$td unless($nid{$nd} or $nidtab{$td}); #warn??
    $keepid{$td}=$nd; 
    $keeptab{$td}="$nd,$keep,$why";
    $anntab{$td}{mustkeep}=$why; $anhdr{mustkeep}++;
  } close(F);
}

sub read_anntxt {
  $inf= $_[0] || "evg5arath.ann.txt";
  open(F,$inf) or die $inf; warn "reading $inf\n";
  while(<F>){ next if(/^\W/); chomp; @v=split"\t"; 
    my($td,$od,$clen,$cdsoff,$aaq)=@v; 
    $anntab{$td}{trlen}=$clen;  $anhdr{trlen}++;
    $anntab{$td}{cdsoff}=$cdsoff; $anhdr{cdsoff}++;
  } close(F);
}

sub read_puboids {
  $inf= $_[0] || "evg5arath.puboids";
  open(F,$inf) or die $inf; warn "reading $inf\n";
  while(<F>){ next if(/^\W/); chomp; @v=split"\t"; 
    my($td,$oid,$gdx,$tix,$clax,$aaq)=@v; 
    $anntab{$td}{oid}=$oid; $anhdr{oid}++;
    $anntab{$td}{aaqual}=$aaq; $anhdr{aaqual}++;
  } close(F);
}

sub read_names {
  $inf= $_[0] || "evg5arath.names";
  open(F,$inf) or die $inf; warn "reading $inf\n";
  while(<F>){ next if(/^\W/); chomp; @v=split"\t"; 
    my($td,$name,$nap,$ndx,)=@v; $nap=~s,/.*,,;
    $anntab{$td}{name}=$name;   $anhdr{name}++;
    $anntab{$td}{namealn}="$ndx,$nap"; $anhdr{namealn}++;
    # $anntab{$td}{namealn}=$nap; $anhdr{namealn}++;
    # $anntab{$td}{nameref}=$ndx;  $anhdr{nameref}++; #? or "$ndx,$nap" ?
  } close(F);
}

sub read_alntab {
  $inf= $_[0] || "evg5arath.cacaorange.alntab";
  open(F,$inf) or die $inf; warn "reading $inf\n";
  while(<F>){ next if(/^\W/); chomp; @v=split"\t"; 
    my($td,$name,$nap,$ndx,)=@v; $nap=~s,/.*,,;
    ($s)=split":",$ndx;
    $anntab{$td}{"orlog_$s"}="$ndx,$nap";  $anhdr{"orlog_$s"}++;
  } close(F);
}

sub read_mapattr {
  $inf= $_[0] || "evg5arath_mrna_input.map.attr";
  open(F,$inf) or die $inf; warn "reading $inf\n";
  while(<F>){ next if(/^\W/); chomp; @v=split"\t"; 
    my($td,$cov,$pid,$nx,$gloc,$npath)=@v; 
    $nx=int($nx); # $nx=~s,\.*,,;
    $mapq= join",", $cov."c", $pid."i", $nx."x"; # use locus above, not gloc?
    $mapq.= ",Spl:$npath" if($npath=~/^C/); # add other paths?
    $anntab{$td}{mapqual}= $mapq; $anhdr{mapqual}++;
  } close(F);
}

__END__

=item about
  evg5arathann.pl : quick perl paste evg5weed annotations together

/bio/bio-grid/plants/arabid/evg5weed/annob

ls evg5arath_icnannbest.idtab evg5arath.keepids evg5arath.ann.txt evg5arath.puboids evg5arath.names evg5arath.cacaorange.alntab evg5arath_mrna_input.map.attr 
evg5arath.ann.txt@		evg5arath.names@		evg5arath_mrna_input.map.attr@
evg5arath.cacaorange.alntab@	evg5arath.puboids@
evg5arath.keepids		evg5arath_icnannbest.idtab

>> new ids, old ids, mainid, idclass (gff), locus
==> evg5arath_icnannbest.idtab <==
Arath5nEVm000001t1	Arath5EVm000001t1,Arath3EVm000001t1	Arath5nEVm000001t1	AtEvg5nmain	chr1:25069487:25095586:-
Arath5nEVm000001t2	Arath5EVm000001t2,Arath4EVm000002t4	Arath5nEVm000001t1	AtEvg5nalt	chr1:25069368:25093094:-
Arath5nEVm000002t1	Arath5EVm000002t1,Arath4EVm000001t1	Arath5nEVm000002t1	AtEvg5nmain	chr3:430769:448665:-
Arath5nEVm000002t2	Arath5EVm000002t2,Arath3EVm000002t2	Arath5nEVm000002t1	AtEvg5nalt	chr3:430982:444938:-

>> mustkeep(nodrop) table for extras (nomaps..)
==> evg5arath.keepids <==
Arath5EVm000068t21	keep	Arath5lEVm033001t1	namedrop
Arath5EVm000179t26	keep	Arath5lEVm033002t1	namedrop
Arath5EVm000181t141	keep	Arath5lEVm033003t1	namedrop
Arath5EVm000209t28	keep	Arath5lEVm033004t1	namedrop

>> orig publicset/ann.txt (no names), keep trlen, cdsoff only
==> evg5arath.ann.txt <==
PublicID	OrigID	TrLen	CDSoff	AAqual	TrGaps	Dbxref	Namealign	Product_Name	CDD_Name
Arath5EVm022860t1	Arath3EVm021214t1	374	3-374	124,99%,partial	0	na	0	Uncharacterized protein	na
Arath5EVm011787t1	Arath3EVm010263t1	1239	24-1064	346,84%,complete	0	na	0	Uncharacterized protein	na
Arath5EVm007665t1	Arath3EVm007299t1	2507	55-1455	466,55%,complete-utrpoor	0	na	0	Uncharacterized protein	na

>> orig publicset/pubids + orig oids of trasm, keep oid, aaqual
==> evg5arath.puboids <==
#Public_mRNA_ID	originalID	PublicGeneID	AltNum	Class	AAqual	pIdAln	Notes
Arath5EVm000001t1	Arath3EVm000001t1,arab6roo3dnta6ahvelvk65Loc1184t5	Arath5EVm000001	1	main	5438,98%,complete	100/100/.	pflag:0,tscore:5438,
Arath5EVm000001t2	Arath4EVm000002t4,arab6cmix3dn4csoapk27loc11152t1	Arath5EVm000001	2	althi1	3274,64%,complete,2X	100/100/.	pflag:0,tscore:3274,
Arath5EVm000001t3	Arath4EVm000002t6,arab6jmix6amvelvk67Loc5372t1	Arath5EVm000001	3	althi1	1576,92%,complete	100/99/.	pflag:0,tscore:1576,

>> use arath align,names
==> evg5arath.names <==
Arath5EVm000001t1	midasin-like protein;	98%,16170/16203,16317	AT1G67120.1	gmap
Arath5EVm000001t2	midasin-like protein;	60%,9850/16203,9828	AT1G67120.1	gmap
Arath5EVm000002t1	auxin transport protein (BIG);	100%,15297/15297,15234	AT3G02260.2	gmap
Arath5EVm000002t2	auxin transport protein (BIG);	62%,9527/15297,9474	AT3G02260.1	gmap

>> non-arath homol vals
==> evg5arath.cacaorange.alntab <==
Arath5EVm000001t1	noname;	47%,10672/16371,16317	cacao:Thecc7EVm000001t2
Arath5EVm000001t1	noname;	43%,9406/15537,16317	orange:ncbig18036614t1
Arath5EVm000001t2	noname;	34%,7926/16371,9828	cacao:Thecc7EVm000001t1
Arath5EVm000001t2	noname;	24%,5116/15537,9828	orange:ncbig18036614t1

>> map qual, 
==> evg5arath_mrna_input.map.attr <==
AQueryID	cov	pid	nexon	GenomeID:gespan:geor	path	oid	tag	qlen
Arath5EVm000001t1	100	100	74.74	chr1:25069487-25095586:-	0	Arath3EVm000001t1	gmap	16565
Arath5EVm000001t2	100	100	66.48	chr1:25069368-25093094:-	0	Arath4EVm000002t4	gmap	15145
Arath5EVm000001t3	100	100	21.18	chr1:25069549-25076931:-	0	Arath4EVm000002t6	gmap	5142


=cut
