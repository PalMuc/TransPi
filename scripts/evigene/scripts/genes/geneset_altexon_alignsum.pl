#!/usr/bin/env perl
# cdsxref.pl : compare ref-exons found, align from blastn tables
# env byhit=1 allg=1 ./cdsxref.pl cacao17uniqcdsx-arath3trsets.blastn

=item about cdsxref
  
  this tabulates results of blastn -db genesets_cdna/cds -query ref_gene_unique_exons.fa
    where ref uniq exons include alternate transcript uniq exons, 
  this summary helps dissect gene set qualities of alternate reconstruction, 
    transcript full/part length (as recovery of start/end exons, vs overall alignment)
  
  relocate as evigene/scripts/genes/geneset_altexon_alignsum.pl
  
  current code is specific to gene sets tested (id patterns): 
    corn x ref(arabid,sorghum), arabid x ref(arabid,cacao,orange)
  
  
=item full usage

  a. tabulate uniq ref exons, eg. from refgenes.gff (CDS) exon locations 
    .. this step is simple exon equal filter
    .. ref gene ids need to have parsable alternate trans suffix: '\.[123]' or 't[123]' now
    .. exons need to be order numberd by target start .. end (e.g. x1=5' to xnnn=3')
        ^^ "ix=NNN" exon tag from gmap2gff.perl, is supposed to be correct transcript order
       
  perl -ne 'if(/^\w/ and /\texon\t/) {
    ($rc,$sg,$tp,$rb,$re,$v,$ro,$xp,$at)=@v=split"\t";
    ($td,$tb,$te)=m/Target=(\S+) (\d+) (\d+)/; ($ix)=m/ix=(\w+)/;
    $xloc=join":",$rc,$rb,$re,$ro; unless($didx{$xloc}) {
    $xid="$td.x$ix"; $xidv{$xid}=$xv="$td:$tb:$te"; $didx{$xloc}=$xv;
    print join("\t",$xid,$xv,$xloc)."\n";  } }' \
      evg7cacao.cds-caca11asm.gmap.gff > evg7cacao_cds.uniqexon.tab
    
  Thecc7EVm024998t1.x1    Thecc7EVm024998t1:1:246 scaffold_1:1508:1753:+
  Thecc7EVm024998t1.x2    Thecc7EVm024998t1:247:319       scaffold_1:1895:1967:+
  Thecc7EVm024998t1.x3    Thecc7EVm024998t1:320:502       scaffold_1:2055:2237:+
  Thecc7EVm024998t1.x4    Thecc7EVm024998t1:503:750       scaffold_1:2497:2744:+
  Thecc7EVm045540t1.x1    Thecc7EVm045540t1:1:61  scaffold_1:155690:155750:+
  Thecc7EVm045540t1.x2    Thecc7EVm045540t1:62:126        scaffold_1:155960:156024:+
  
  a2. improved exon-equal filter, drop any span overlaps
    cat evg7cacao_cds.uniqexon.tab | cut -f1,3 | sed 's/:/  /g;' | sort -k2,2 -k3,3n -k4,4nr | \
      grep '       scaff' | perl -ne '($xd,$rs,$rb,$re,$ro)=split; next if($rs eq $lrs and $rb<$lre and $re > $lrb);
        print; ($lrs,$lrb,$lre)=($rs,$rb,$re);' \
      > evg7cacao_cds.uniqexon.unilocs
    perl -ne '($xd)=@v=split; if(@v==5) { $ok{$xd}=1 } else { print if($ok{$xd}); }' \
      evg7cacao_cds.uniqexon.unilocs evg7cacao_cds.uniqexon.tab > evg7cacao_cds.uniqexonok.tab
  
  b. pull uniqexon.cds.fasta (with or without gene id selections)
  
  env allx=1 minw=50 perl -ne 'if(/^(Thec\w+)$/) { $gok{$1}=1; $GOK=1; } 
  elsif(/^Thec\S+\tThec/){ ($xid,$tbe,$gloc)=split; ($td,$tb,$te)=split":",$tbe; ($gd=$td)=~s/t\d+$//; 
  next if($GOK and not $gok{$gd});  $tok{$td}{$xid}=[$tb,$te,$gloc]; } else{ if(/^>(\S+)/) { 
  $ltd=$1; putx($td,$fa) if($td and $fa); $td=$ltd; $sok=0; $fa=""; $sok=($tok{$td})?1:0;  } 
  elsif($sok) { chomp; $fa.=$_; }  } sub putx{my($td,$fa)=@_; @xd=sort keys %{$tok{$td}}; 
  for $xd (@xd) { my($tb,$te,$gl)=@{$tok{$td}{$xd}}; $w=1+$te-$tb; $xs=substr($fa,$tb-1,$w);
   print ">$xd xlen=$w; xloc=$tb-$te/$gl;\n$xs\n" if($w>=$MINW);  } } BEGIN{ $MINW=$ENV{minw}||20; }' \
     evg7cacao_cds.uniqexon.hasalt.gids evg7cacao_cds.uniqexon.tab evg7cacao.cds \
     > evg7cacao_cds_hasalt.uniqexons.fa
  
  c. blastn -query uniqexons.fa -db genesets_cdna
  
    geneset_db=arath4pb2genes.cdna
    $nbin/blastn -evalue 1e-9 -outfmt 7 -task blastn \
      -query evg7cacao_cds_hasalt.uniqexons.fa -db $geneset_db -out cacao17uniqcdsxya-$geneset_db.blastn 
  
  d. cdsxref tabulations by gene set of uniqcds.blastn
    
    d1. -- exons hit per alternate transcripts, restrict to ref genes found in all gene sets:
    env byhit=1 allg=1 ./cdsxref.pl arath16uniqcdsx-arath4pb2genes.cdna.blastn  
  ref exons n=111649, AllSrcHitGene=1, aveByHit=1, isTrID=0, AlnRef=170691
  Src     nGene   nAlt    nXhit   pX      pRefAln AlnX    IdnX    pIdnX
  EVm.Ar  14170   15460   88393   79.2    94.1    267.0   265.3   99.5
  atap16  14170   16821   97596   87.4    100.0   301.9   301.9   100.0
  pacbc   14170   14701   73891   66.2    91.4    260.7   255.7   98.2
  pacbm   14170   14759   74475   66.7    91.5    260.6   255.7   98.2
  ..........

    d2a. -- exons hit per alternate, not restricted ref-genes in all gene sets
    env byhit=1 allg=0 ./cdsxref.pl arath16cds-arath4pb2genes.cdna.blastn  
    d2b. -- transcripts hit,  blastn -query refgenes_cdna instead of refgenes_uniqexons
    env istr=1 byhit=1 allg=0 ./cdsxref.pl arath16cds-arath4pb2genes.cdna.blastn 
    
    d3. -- exons hit per exon position (5', 3' and each found 1..max)
    .. this is restricted to ref transcripts where all gene sets have 1+ exon hit
    env exons=15 ./cdsxref.pl cacao17uniqcdsxya-arath4pb2genes.cdna.blastn

      geneset_db = EVm.Ar  atap16  pacbc   pacbm
  arath16uniqcdsx-arath4pb2genes.cdna.blastn
  Exons Hit by position in transcript
  Totals  gn=20706        tn=30993        xn=140746       gnt=27643       trt=34122
  Exon    nHit    EVm.Ar  atap16  pacbc   pacbm   EVm.Ar  atap16  pacbc   pacbm
  -5      18959   92.2    100     77.9    78.2    17497   18959   14783   14843
  -3      18959   97.2    100     93.7    93.8    18435   18959   17767   17798
  1       24933   91      100     75.2    75.5    22699   24933   18757   18834
  2       20413   95.7    100     81.8    82.1    19546   20413   16716   16778
  3       16306   96.1    100     83.1    83.4    15673   16306   13562   13614
  4       13303   96.1    100     84.6    84.9    12792   13303   11263   11299
  ..........
    
=cut

# fixme: collect ref-exon sizes from # QUery xlen=nnn if there, calc %align/refsize
# # Query: Thecc7EVm036769t1.x1 xlen=121; xloc=1-121/scaffold_1:27673:27793:+;

# BEGIN  
$ISTRID=$ENV{istr}||0; $ALLG=$ENV{allg}||0; $BYHIT=$ENV{byhit}||0; 
$SHOWDIFF=$ENV{diff}||0;
$REFLEN=$ENV{reflen}||0;
$EXONHITS=$ENV{exons}||0;
$DEBUG=$ENV{debug}||0;

if($refexontab= $ENV{reftab}) { readRefExonTab($refexontab); } # %xord set

## fixme : corn pacbio cshl16= csh6roo2: csh6ea  csh6em  csh6en  csh6po  csh6ro  csh6ta << all one source
sub tsrc { 
 my($td)=@_; 
 $td=~s/(hdfpac|hd5fpac)/pac/;
 my($ts)= ($td=~m/^(\w+)(EVm)\d/) ? "$2.$1" 
 	  : ($td=~m/^(\w+):/) ? $1 
  	: ($td=~m/^\w+(pacbc|pacbm)/) ? $1
  	: ($td=~m/^(SRR\d+)([A-Za-z]+)/) ? "$2.$1"
  	: substr($td,0,3);
  $ts=~s/^(csh6)\w+/${1}pb/; # corn pacbio cshl16= csh6ear1: csh6roo2: csh6endu:  ..
  if(1) { ($ts)= substr($ts,0,6); }
  return($ts);
}

=item readRefExonTab, check/fix exon order

#  current exon order may be wrong, check here w/ exon.tab transcript locations
Thecc7EVm033748t1.x2    Thecc7EVm033748t1:247:462       scaffold_1:1712152:1712367:-
Thecc7EVm033748t1.x3    Thecc7EVm033748t1:463:543       scaffold_1:1709916:1709996:-
                   ^-----^ this is proper xord
Thecc7EVm009179t1.x1    Thecc7EVm009179t1:1:123 scaffold_1:2154699:2154821:-
Thecc7EVm009179t1.x2    Thecc7EVm009179t1:124:562       scaffold_1:2154164:2154602:-
                  ^-----^ this is proper xord
evg7cacao_cds.uniqexon.tab
Thecc7EVm024998t1.x1    Thecc7EVm024998t1:1:246 scaffold_1:1508:1753:+
Thecc7EVm024998t1.x2    Thecc7EVm024998t1:247:319       scaffold_1:1895:1967:+
Thecc7EVm024998t1.x3    Thecc7EVm024998t1:320:502       scaffold_1:2055:2237:+
Thecc7EVm024998t1.x4    Thecc7EVm024998t1:503:750       scaffold_1:2497:2744:+
Thecc7EVm045540t1.x1    Thecc7EVm045540t1:1:61  scaffold_1:155690:155750:+
Thecc7EVm045540t1.x2    Thecc7EVm045540t1:62:126        scaffold_1:155960:156024:+
Thecc7EVm045540t1.x3    Thecc7EVm045540t1:127:203       scaffold_1:156122:156198:+

=cut

sub readRefExonTab {
  my($inf)=@_;
  %xord=(); # global
  my(%tstart,%xid); 
  my($nt,$tok,$tomis)=(0) x 9; 
  open(F,$inf) or return {};
  while(<F>) { next if(/^\W/);
    my($xd,$tloc,$gloc)=split;
    my($td,$tb,$te)=split":",$tloc; 
    my($xi)= $xd=~m/x(\d+)$/;
    $tstart{$td}{$tb}=$xi; 
    $xid{$td}{$xi}= $xd;
  } close(F);
  for my $td (sort keys %tstart) { 
    my @tb=sort{$a <=> $b} keys %{$tstart{$td}}; 
    my($xj,$ok,$omis)=(0) x 9;
    for my $i (0..$#tb) { 
      my $xi= $tstart{$td}{$tb[$i]}; 
      my $xid= $xid{$td}{$xi};
      $xord{$xid}= 1+$i; 
      if($xi>$xj){$ok++;} else{ $omis++; } 
      $xj=$xi; 
    } 
    $nt++; if($omis) { $tomis++; } else { $tok++; }
  } 
  warn "#dbg: ntrans=$nt, xorder err:$tomis, ok:$tok\n"; # if $DEBUG;
  return \%xord;
}

while(<>) {
  if(/^\W/) { 
    if(/^# Query:/ and /(length|len)=\d/i){ 
      ($rd)=m/Query: (\S+)/; ($rw)=m/(?:length|len)=(\d+)/i; 
      $rd=~s/Sobic\./Sobic/; $rlen{$rd}=$rw; $REFLEN++ if($rw); 
      #nogo# if(($tb,$te)= m/xloc=(\d+)-(\d+)/) { $xord{$rd}=$tb; }
    } next; }
  ($rd,$td,$pi,$al)=@v=split;  
  # full blasttab: ($rd,$td,$pi,$al,$mi,$ind,$rb,$re,$tb,$te,$ev,$bits)=@v;
  $rd=~s/Sobic\./Sobic/;
  if(@v==1) { $gok{$rd}=1; $GOK= 1;  next; } # geneID okay list
  ($ts)= tsrc($td);  # substr($td,0,3)
  @rd= split/\./,$rd;
  push @rd,"x1" if($ISTRID);
  if(@rd==3){ ($rg,$rt,$rx)= @rd; }
  elsif(@rd==4){ ($rss,$rg,$rt,$rx)= @rd; $rg="$rss.$rg"; }
  elsif(@rd==2 and $rd[0]=~m/t\d+$/){ # Evigene EVm0001t2.x3 format instead of AT1G0001.2.x3
    ($rg,$rx)= @rd; ($rt)= ($rg=~s/t(\d+)$//)?$1:0; 
    $rdorig=$rd; $rd="$rg.$rt.$rx"; # FIXME this and above rlen{rd} need to match
    $rdorig{$rd}=$rdorig;
  }
  else { warn"#bad refid: $rd, $td\n"; next; }
  next if($GOK and not $gok{$rg}); 
  $rgt="$rg.$rt"; $rx=~s/x//;  
  if(%xord) { $rxo= $xord{$rd} || $xord{$rdorig{$rd}}; # have input uniqexon.tab 
    if($rxo) { $rx=$rxo; } else { $nmisxord++; warn "# miss xord $rd\n" if $DEBUG; } # bug blastn has refexons not in uniqexon.tab
    }  
  #nogo# if( $xo=  $xord{$rd} || $xord{$rdorig{$rd}} ) { $rx= $xo; } # bad: not 1,2,3 but tpos 1,99,299,..
  if($EXONHITS) {
    $ts{$ts}++; # == below
    ## merge w/ below; replace ghit/thit with trhit{$rt}
    $trhit{$rgt}{$ts}++;
    # $ghit{$rg}{$ts}++; # == $rg{$rg}{$ts}++;
    # $thit{$rg}{$rt}{$ts}++;  # or thit{$rg.$rt}{$ts} ?
    $xhit{$rgt}{$rx}{$ts}++; # ~= val{$rd}{$ts}{n} or didx{$rd}{$ts}
  } else {
    $xg{$rd}= $rg;  $ts{$ts}++;
    unless($didx{$rd}{$ts}) {
      $xok=1; if($xdid=$didg{$td} and $xdid ne $rgt) { $xok=0; }
      if($xok) { $didg{$td}=$rgt; $rg{$rg}{$ts}++; $idn= $al*$pi/100;
      $val{$rd}{$ts}{n}=1; $val{$rd}{$ts}{al}=$al; $val{$rd}{$ts}{idn}=$idn; $val{$rd}{$ts}{pid}=$pi;
      print if($ENV{dopr}); } $rds{$rd}++; $didx{$rd}{$ts}++; 
      }
  }
}

if($EXONHITS) { printexonhits(); }
else {
  @ts=sort keys %ts; $nts=@ts; @rx=sort keys %didx;
  if($ALLG){ my %okg; for $g (keys %rg){ @t=keys %{$rg{$g}}; $okg{$g}=(@t == $nts); }
    @rx=grep{ $okg{$xg{$_}} }@rx; }
  $nrx=@rx; print "ref exons n=$nrx, AllSrcHitGene=$ALLG, aveByHit=$BYHIT, isTrID=$ISTRID, AlnRef=$REFLEN\n";
  %sum=(); %gts=(); for $rx (@rx){ @t=keys %{$val{$rx}}; for $t (@t) {
   for $k (qw(n al idn pid)) { $sum{$t}{$k}+=$val{$rx}{$t}{$k}; }
   if($REFLEN) { 
    my $aln= $val{$rx}{$t}{al}||0;
    my $rlen= $rlen{ $rdorig{$rx} }||$rlen{$rx}; 
    my $alnr= ($rlen>0) ? 100*$aln/$rlen : 0; $alnr=100 if($alnr>100); 
    $sum{$t}{'alr'}+= $alnr;
    }

   $hit=$val{$rx}{$t}{n}||0;
   my @idp= split/\./,$rx; shift @idp if(@idp==4); ($gi,$ti,$xi)=@idp;
   if($hit) { $gts{$t}{gn}{$gi}++; $gts{$t}{tn}{$gi.$ti}++; $gts{$t}{xn}{$gi.$ti.$xi}++; }
   } 
   if($SHOWDIFF and @ts==2) {
    ## @t == found both vs @ts == can be missing
    my($t0,$t1)=@ts; my($v0,$v1)= map{ $val{$rx}{$_}{al}||0 } @ts;
    my $da= $v0 - $v1; my $vm=($da<0)?$v1:$v0;
    if(abs($da) > 0.5*$vm) {
      print join("\t","diffaln",$da,$t0,$v0,$t1,$v1,$rx)."\n";
    }
   }
  }
  @var= qw(al idn pid); @vlab=qw(AlnX IdnX pIdnX);
  if($REFLEN) { unshift @var, 'alr'; unshift @vlab, 'pRefAln'; }
  print join("\t",qw(Src nGene nAlt nXhit pX),@vlab)."\n";
  for $ts (@ts) {
    $n=$c=$sum{$ts}{n}; $n||=1; $pc= sprintf"%.1f",100*$c/$nrx;
    $nav=($BYHIT)?$n:$nrx;
    @av= map{ sprintf"%.1f", $sum{$ts}{$_}/$nav; } @var;
    ($ng,$nt,$nx)=map{ scalar(keys %{$gts{$ts}{$_}}); } qw(gn tn xn);
    print join("\t",$ts,$ng,$nt,$c,$pc,@av)."\n";
  } 
} 

sub readexonhits {
  $EXONHITS=1; 
  while(<>) { # see above reader
    ($rd,$td)=split; ($rg,$rt,$rx)=split/\./,$rd;
    ($ts)= tsrc($td);  # substr($td,0,3)
    #.. ^^ as above..
    if($EXONHITS) {
      $rx=~s/x//;  
      $ghit{$rg}{$ts}++;
      $thit{$rg}{$rt}{$ts}++; 
      $xhit{$rg}{$rt}{$rx}{$ts}++; $ts{$ts}++; 
    }
  }
}

sub printexonhits {
  ## limit output to ? nexons=20? or those w/ >nn hits?
  ## modify all @ts for @rga= ghit, @tia= thit : foreach @ts pick any thit, count completeness of exons hit
  my $EXONMAX= ($EXONHITS>1)? $EXONHITS : 20;
  @ts=sort keys %ts;  @rg=sort keys %trhit; 
  $sum{gnt}=@rg; 
  if(1 or $ALLG){ @rg= grep{ @t=keys %{$trhit{$_}}; (@t==@ts) } @rg; } # keep always ?
  for $rg (@rg) {
    $sum{gn}++;  
    ($ti)= $rg=~m/(\d+)$/;
    @xi=sort{$a <=> $b} keys %{$xhit{$rg}};   
    $xe=($ti == 1)?$xi[-1]:99999; $xb=($ti == 1)?$xi[0]:999991; 
    for $xi (@xi) { $sum{xn}++; ($ix=$xi)=~s/x//; 
      $hitx{-5}{n}++ if($xi eq $xb);
      $hitx{-3}{n}++ if($xi eq $xe); $hitx{$ix}{n}++; 
      for $ts (@ts) {
        $h=$xhit{$rg}{$xi}{$ts}; if($h){ $hitx{-3}{$ts}++ if($xi eq $xe);
        $hitx{-5}{$ts}++ if($xi eq $xb); $hitx{$ix}{$ts}++; } 
      } 
    } 
  } 
  
  @sm= map{ $c=$sum{$_}; "$_=$c" } qw(gn xn gnt); 
  @xi=sort{$a <=> $b} keys %hitx; $nxi=@xi;
  print "Exons Hit by position in transcript (nxi=$nxi, showmax=$EXONMAX)\n";
  print join("\t","Totals",@sm)."\n"; 
  print join("\t",qw(Exon nHit),@ts)."\n"; 
  for $xi (@xi) { 
    last if($xi > $EXONMAX);
    $nx=$hitx{$xi}{n}; @c=map{ $hitx{$xi}{$_} } @ts; 
    @p=map{ int(1000*$_/$nx)/10 } @c; print join("\t",$xi,$nx,@p,@c)."\n"; 
  }
}

sub OLD_printexonhits {
  ## limit output to ? nexons=20? or those w/ >nn hits?
  ## modify all @ts for @rga= ghit, @tia= thit : foreach @ts pick any thit, count completeness of exons hit
  my $EXONMAX= ($EXONHITS>1)? $EXONHITS : 20;
  @ts=sort keys %ts;  @rg=sort keys %ghit; 
  @rga= grep{ @t=keys %{$ghit{$_}}; (@t==@ts) } @rg; $sum{gnt}=@rg; 
  for $rg (@rga) {
    @ti=sort{$a<=>$b} keys %{$thit{$rg}}; 
    @tia= grep{ @t=keys %{$thit{$rg}{$_}}; (@t==@ts); } @ti; 
    $sum{trt}+= @ti; $sum{gn}++; 
    for $ti (@tia) { 
      $sum{tn}++; @xi=sort{$a <=> $b} keys %{$xhit{$rg}{$ti}};   
      $xe=($ti == 1)?$xi[-1]:99999; $xb=($ti == 1)?$xi[0]:999991; 
      for $xi (@xi) { $sum{xn}++; ($ix=$xi)=~s/x//; 
        $hitx{-5}{n}++ if($xi eq $xb);
        $hitx{-3}{n}++ if($xi eq $xe); $hitx{$ix}{n}++; 
        for $ts (@ts) {
          $h=$xhit{$rg}{$ti}{$xi}{$ts}; if($h){ $hitx{-3}{$ts}++ if($xi eq $xe);
          $hitx{-5}{$ts}++ if($xi eq $xb); $hitx{$ix}{$ts}++; } 
        } 
      } 
    } 
  } 
  
  @sm= map{ $c=$sum{$_}; "$_=$c" } qw(gn tn xn gnt trt); 
  @xi=sort{$a <=> $b} keys %hitx; $nxi=@xi;
  print "Exons Hit by position in transcript (nxi=$nxi, showmax=$EXONMAX)\n";
  print join("\t","Totals",@sm)."\n"; 
  print join("\t",qw(Exon nHit),@ts)."\n"; 
  for $xi (@xi) { 
    last if($xi > $EXONMAX);
    $nx=$hitx{$xi}{n}; @c=map{ $hitx{$xi}{$_} } @ts; 
    @p=map{ int(1000*$_/$nx)/10 } @c; print join("\t",$xi,$nx,@p,@c)."\n"; 
  }
}

__END__


=item exon first,last counts
  add here? same input blastn?
  
grep -v '^#' ./corn/geneval/trsetevg5/arath16uniqcdsx-corn4gsetcdsb.blastn | \
grep -v '^#' ./corn/geneval/trsetevg5/sorghumuniqcdsx-corn4gsetcds.dcblastn | \
perl -ne 's/Sobic./Sobic/; ($rd,$td)=split; ($rg,$rt,$rx)=split/\./,$rd;
$rx=~s/x//;  ($ts)=substr($td,0,3); $ghit{$rg}{$ts}++;
$thit{$rg}{$rt}{$ts}++; $xhit{$rg}{$rt}{$rx}{$ts}++; $ts{$ts}++; END{
@ts=sort keys %ts;  @rg=sort keys %ghit; @rga= grep{ @t=keys
%{$ghit{$_}}; (@t==@ts) } @rg; $sum{gnt}=@rg; for $rg (@rga) {
@ti=sort{$a<=>$b} keys %{$thit{$rg}}; @tia= grep{ @t=keys
%{$thit{$rg}{$_}}; (@t==@ts); } @ti; $sum{trt}+= @ti; $sum{gn}++; for
$ti (@tia) { $sum{tr}++; @xi=sort{$a <=> $b} keys %{$xhit{$rg}{$ti}};   
$xe=($ti == 1)?$xi[-1]:99999; $xb=($ti == 1)?$xi[0]:999991; for $xi
(@xi) { $sum{xn}++; ($ix=$xi)=~s/x//; $hitx{-1}{n}++ if($xi eq $xb);
$hitx{0}{n}++ if($xi eq $xe); $hitx{$ix}{n}++; for $ts (@ts) {
$h=$xhit{$rg}{$ti}{$xi}{$ts}; if($h){ $hitx{0}{$ts}++ if($xi eq $xe);
$hitx{-1}{$ts}++ if($xi eq $xb); $hitx{$ix}{$ts}++; } } } } } @sm= map{
$c=$sum{$_}; "$_=$c" } qw(gn tr xn gnt trt); print
join("\t","Totals",@sm)."\n";  @xi=sort{$a <=> $b} keys %hitx; print
join("\t",qw(Exon nHit),@ts)."\n"; for $xi (@xi) { $nx=$hitx{$xi}{n};
@c=map{ $hitx{$xi}{$_} } @ts; @p=map{ int(1000*$_/$nx)/10 } @c; print
join("\t",$xi,$nx,@p,@c)."\n"; } }  ' | head -25
 
 -- supposed to count only genes, tr/gene where all sources hit .. check

=cut


=item fixme targ ids need opts to pick tsource

==> arath16ap.mrna <==
>atap16:AT1G01010.1 | NAC domain containing protein 1 | Chr1:3631-5899 FORWARD LENGTH=1688 | 201606

==> evg1arath.mrna <==
>Arath1EVm019457t1 type=mRNA; aalen=144,81%,partial5; clen=537; offs=3-437; oid=tidbarab6roo3dno1sridk107Loc104763; organism=Arabidopsis_thaliana;

==> pbclass_set1.mrna <==
>SRR3655759hdfpacbm150514100_41_1121_CCS type=mRNA; aalen=84,23%,complete-utrbad; clen=1080; strand=+; offs=118-372; fiveseen=1;polyAseen=1;threeseen=1;fiveend=41;polyAend=1121;threeend=1154;primer=1;chimera=0

==> pbclust_set2.mrna <==
>SRR3655759hdfpacbc1023_21_1009 type=mRNA; aalen=266,79%,complete; clen=1009; strand=+; offs=73-873; isoform=c1023;full_length_coverage=21;isoform_length=1009

=cut


