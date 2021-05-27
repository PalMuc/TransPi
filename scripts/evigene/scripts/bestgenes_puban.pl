#!/usr/bin/env perl
# bestgenes_puban.pl
# ?? merge into  evigene/scripts/bestgenes_update.pl
# cat aphid2pub8d-{attr,homolog,express,gequal,name}.tab | env table=1 ../bestgenes_puban.pl > aphid2pub8d.attr.tbl

#scoretype: homolog:9,paralog:1,ref:2,est:3,pro:3,rseq:2,intr:20,nintron:40,inqual:20,maqual:5,terepeat:-3,UTR:3,CDS:1

use strict;
my $jTEscore= 10; # index from scorevec in gene.gff ??

# ** FIXME-maybe: keep protein= in mRNA attribute, some now are curated: add that quality=Protein:Curated
my @vkey=qw(ID osrc oid gene alttr aalen cxlen inexon inqual maqual  join must scorevec
            chimera chim1 chim2 egover ); # 4 from chimera
my @atkey=qw(ID gene isoform  quality aaSize cdsSize Name Dbxref express ortholog paralog intron oid chimera score scorevec);
    #  change col header for Microstupid ExcelNOT:
my %recode_key = ( ID => "transcriptID", gene => "geneID");
my $MISSING= $ENV{missing} || 0; #"."; # or na, for TABLE only

# should documnet fields, in comments at table output top?

## fixme, new attr: 
# desc= from chimera > Name= ; egover= from chimera: save where?

## add Dbxref = equivalent acypi1,ncbiref2 genes
# .. relabel:  Ref1: == ACYPI or APHIDBASE: ?  Ref2=RefSeq: or RefSeq:
## add location to table output: chr:b-e:o
## fixme? missing Name, Dbxref should have null value?
## fixme? Express:Strong with few reads but spanning tr : should be :Medium unless reads/est > min?

my $formatTABLE = $ENV{table} || 0;
my $debug= $ENV{debug} || 0;
my $recase= $ENV{recase} || 0;

my $GENE_NONAME = "hypothetical protein";

# need to match aphid2pub8d-gequal.tab columns: use col header as REFtag ??
my $REF1tag = "APHIDBASE";  # ref1gene REF1tag = APHIDBASE
my $REF2tag = "RefSeq";     # ref2gene REF2tag = RefSeq

my ($tb, %gattr, @ids, @ha, @hn, @hg, @hh, @hx, @hd);

splice(@atkey,12,0,"location") if($formatTABLE);

my @TEnames= getTEnames();

if($ENV{count}) { namecount(); }
else { processtables(); }

#...................................................

sub getattr { 
  my $id= shift; 
  my $attr= $gattr{$id};  
  unless(ref $attr) { 
    $attr={}; $gattr{$id}= $attr; 
    if($formatTABLE) { $$attr{Name}= $$attr{Dbxref}= "na"; } #??
  }
  return $attr;
}


sub processtables 
{
while(<>) {
  chomp; my @v=split"\t"; 
  
  if(/^\W/) { next; }
  elsif(/^ID\tosrc/){ $tb="a"; @ha=@v; push @hd,@v; } 
  elsif(/\tgene_id\tname\t/) { $tb="n"; splice(@v,0,2); @hn=@v; push @hd, @v; } 
  elsif(/\tgene_id\tortholog\t/) { $tb="h"; splice(@v,0,2); @hh=@v; push @hd, @v; } 
  elsif(/\tgene_id\txcover\t/) { $tb="x"; splice(@v,0,2); @hx=@v; push @hd, @v; } 
  elsif(/\tgene_id\tref1gene\t/) { $tb="g"; splice(@v,0,2); @hg=@v; push @hd, @v; } 
  elsif(/^gene/i or /\tgene_id\t/i) { $tb="HUH?"; warn "# SKIPPING. Dont know table input: $_\n"; }
  
  elsif($tb eq "a") {
    my %v= vmap(\@ha,\@v); 
    
    my $id= $v{ID}; push(@ids, $id);
    my $attr= getattr($id); # $gattr{$id};  unless($attr) { $attr={}; $gattr{$id}= $attr; }
    $$attr{ID}= $id;   # RECODE for Table: Damn Microstupid Excel barfs on ID column name
    $$attr{gene}= $v{gene}; # can be null; fix?  # recode for table : gene_id:
    $$attr{aaSize}= $v{aalen};  
    $$attr{aaSize} =~ s/,[a-z][\w-]+//; # aalen,pctcode,quality : 145,35%,curated-complete
    
    my $cxlen= $$attr{cdsSize}= $v{cxlen};  # cds/tr
    $$attr{intron}= $v{inexon}; # in/ex
    $$attr{oid}="$v{osrc}:$v{oid}"; 
    $$attr{score}= $v{score}; 
    $$attr{location}= $v{location} if($v{location} and $formatTABLE); 
    
    my($clen,$trlen)= $cxlen=~ m/(\d+).(\d+)/;  $trlen||=1;
    $$attr{trlen}= $trlen;
    $$attr{isoform} = 0;
    if($v{alttr}) { my($t)=$v{ID} =~ m/t(\d+)$/; $$attr{isoform} = $t||$v{alttr}; }
    ($$attr{qualIntron} = $v{inqual}) =~ s,^.*/,, if $v{inqual};
    ($$attr{qualMated}  = $v{maqual}) =~ s,^.*/,, if $v{maqual};
    
    # ($$attr{qualProtein} = $v{aalen}) =~ s/^\d+,\d+.,// if $v{aalen};
    # ^^ add qualProtein=Poor  for ncRNA-like: utr too long or too many utr exons, an aalen < minaa
    # FIXME: add qualProtein=curated-complete?
    if( $v{aalen} =~ /,/ ) {
      my($aal,$aap,$aaq)= split",",$v{aalen};
      $aap =~ s/\%//;  
      if( $aap =~ /\d/ and $$attr{qualHomology} !~ /Homology:Ortholog(Strong|Medium)/) {
        if($aap < 16) { $aaq= "noncode_$aaq"; } elsif($aap < 34) { $aaq= "poor_$aaq"; }
        }      
      $aaq =~ s/curated-complete/curated_complete/;
      $$attr{qualProtein}= $aaq;
    }
     
    $$attr{qualExpertchoice} = 1 if($v{must} > 0);  # must= 77 or 69 (alt)
    $$attr{qualFusionMaybe} = 1 if($v{'join'} =~ /\w/); # join=l1/p2vi1vi1l1/p1p2/..

    $$attr{chimera}=0;
    if($v{chimera}) {
      (my $chi = $v{chimera} ) =~ s/,.*//; #chi=1 or 2, place chim1/2 here
      $chi .= ",".$v{chim1} if($v{chim1});
      $chi .= ",".$v{chim2} if($v{chim2});
      my $eg= ($v{egover}); $eg=~s/:\d+//g; # keep this or drop?
      $chi .= ",".$eg if ($eg); 
      $$attr{chimera}= $chi; 
    }
    
    my $scorevec= $v{scorevec}; 
    my @scorevec=split ",", $scorevec;
    $$attr{scorevec}= join",", map{ s,/.*,,; $_ } @scorevec;
    
    ## fixme here? set if >= minTEbases, but include other criteria : Express:Strong, ..
    my $pte= $scorevec[$jTEscore] / $trlen;
    $$attr{qualTransposon} = 1 if($pte > 0.33); 

  } elsif($tb eq "n") {
    my($g0,$id,@v)= @v;
    my %v= vmap(\@hn,\@v); 
    my $attr= getattr($id); 
    my $na= $v{name};  # "na" is no name
    
    $na= nameclean($na); # use qualities to decide for nonames: 
    my $iste= isTEname($na);
    $$attr{qualTransposon} = 1 if($iste); 
    
    if(my $pna= $v{namepct}) { 
      $pna =~ s/[CI]//;
      $pna =~ s/$/\%/ unless($pna =~ m/\%/);
      $na .=  " ($pna)"; 
      }
    
    $$attr{Name} = $na;

  } elsif($tb eq "g") {
    my($g0,$id,@v)= @v;
    my %v= vmap(\@hg,\@v); 
    my $attr= getattr($id); 
    my $na = $v{ref1gene}; # maybe na REF1tag = APHIDBASE
    my $na2= $v{ref2gene}; # maybe na REF2tag = RefSeq

    ## FIXME: list id1,id2,..  : chop off or apply to all 
    ## FIXME2: separate GNOMON: and RefSeq: i.e. keep best of both
    ## FIXME3: /I and/C flags indicate I/perfect or C >= 90% coding identity: keep this flag
if(1) {
    # na == ref1 = aphidbase now
    my(@na2,$r2,$r3); # refseq,gnomon now, preserve best of each, in order (1st == best)
    foreach (split",",$na2) {
      if(/:/) { push @na2,$_ unless($r2++); }
      elsif(/\w/ and $_ ne "na") { push @na2, "$REF2tag:$_" unless($r3++); }
    }
    
    map{ s/,.*//; my($n,$t,$p)=m,(.+)/([CI]?)(\d+),; $_="$n,$p\%$t" if($p); } ($na, @na2);

    $na  = ($na eq "na") ? "" : ($na =~ /\w/) ? "$REF1tag:$na" : "";    
    $na2 = join ",", @na2;
    
} else {
    map{ s/,.*//; my($n,$t,$p)=m,(.+)/([CI]?)(\d+),; $_="$n,$p\%$t" if($p); } ($na,$na2);  # ?? chop off "/pscore" ?
    $na  = ($na eq "na") ? "" :($na =~ /\w/) ? "$REF1tag:$na" : "";    
    $na2 = ($na2 eq "na") ? "" : ($na2 =~ /^GNOMON:/) ? $na2 : ($na2 =~ /\w/) ? "$REF2tag:$na2" : "";
}    
    $$attr{Dbxref}= ($na and $na2) ? "$na,$na2" : $na.$na2;

  } elsif($tb eq "h") {
    my ($g0,$id,@v)= @v;
    my %v= vmap(\@hh,\@v); 
    my $attr= getattr($id); 
    
    ## fixme: ortholog=0%,0,UniProt:na,Arp:na << drop if 0 bits, and/or "na"
    $$attr{ortholog} = ($v{ortholog} =~ /0%,0,/) ? "" : $v{ortholog};
    $$attr{paralog}  = ($v{paralog} =~ /0%,0,/ ) ? "" : $v{paralog};

    my $pho= $v{pbest};    
    my $qho= ($pho =~ /O/)?"Ortholog":($pho =~ /P/)?"Paralog":"";
    $pho =~ s/\D+//;
    $qho .= pqual($pho);
    $$attr{qualHomology} = $qho;
    if($qho =~ /Ortholog(Strong|Medium)/ and $$attr{qualProtein} =~ /poor|noncod/i) {
      $$attr{qualProtein} =~ s/(poor|noncode)[_]?//i;
      }

  } elsif($tb eq "x") {
    my($g0,$id,@v)= @v;
    my %v= vmap(\@hx,\@v);
    my $attr= getattr($id);
    my $trlen= $$attr{trlen}||1;
    
    # express: this is null value : 1       0       1,0r    
    # gives: express=0%,0r; .. is this right or should be no value?
    my $px= pct($v{xcover} / $trlen),
    my $ec= $v{est_cover}; # cov,%ident
    (my $nr= $v{rna_cover}) =~ s/^\d+,//;    
    $$attr{express} = "$px\%,$nr";
    
    my $q= pqual($px); $q="Medium" if($q eq "Strong" and ($nr*70)/$trlen < 1.5);
    $$attr{qualExpress} = $q;
  }

}

#? as table or gff attr?
foreach my $id (@ids) {
  putt( $id);
}

}

my $didhead=0;

sub putt {
  my($g)= @_;
  my @v; my %didk;
  my $attr= $gattr{$g};
  my @ak= sort keys %$attr;
  my @qk= grep /^qual/, @ak;
  my $qual= join ",", map{ $didk{$_}++; (my $k=$_ ) =~ s/^qual//; "$k:$$attr{$_}" } @qk;
  
  #? here set/reset Name = hypothetical protein depending on quality of evidence
  #    hypoth = weak evid;  uncharacterized = strong evd;  ? expressed uncharact = strong express
  my $na= $$attr{'Name'};
  $na= "hypothetical protein" if($na eq "na" or $na !~ /\w/);
  if($na =~ /hypothetical protein/i) {
    $na =~ s/hypothetical protein/uncharacterized protein/i
      if($qual =~ /(Express|Homology):\w*(Strong|Medium)/i);
    #?? on or off# $na .=", expressed" if($qual =~ /Express:(Strong|Medium)/i);
    $$attr{'Name'}= $na;  
  }
  
  foreach my $k (@atkey) {
    my $v=$$attr{$k}; $didk{$k}++;
    my $kl= $k;
    if($k eq "quality") { $v=$qual; }
    $v= $MISSING if($formatTABLE and not defined $v);
    push @v, (($formatTABLE) ? $v : "$kl=$v") if($formatTABLE or $v); 
    }

# if(0) {    
#   foreach my $k (sort keys %$attr) {
#     next if $didk{$k};  my $v=$$attr{$k}; push @v, (($formatTABLE) ? $v : "$k=$v") if($v); 
#     }
# }

  my $tab= ($formatTABLE) ? "\t" : ";";
  if($formatTABLE and 1 > $didhead++) {  
    my @colhead = map { $recode_key{$_} || $_ } @atkey;
    print join( $tab, @colhead),"\n"; 
    }
  print join( $tab, @v),"\n";
}

## add TE protein name classing; also quality flag Transposon


sub nameclean {
  local $_= shift;
  ## FIXME:  3&apos,-5&apos, exoribonuclease
  study;
  return $GENE_NONAME if($_ eq "na");
  
  s/\&apos[;,]/'/g;
  s/;/,/g; s/=/:/g; 
  s/PREDICTED:\s*//i;  s/\(predicted\)//i;
  s/,GN:[^\s,;]+//;  # ? drop/keep EC?
  s/,EC:[^\s,;]+//;
  s/Isoform [\w-]+ of //i; s/\s*(isoform|isozyme.*) [\w-]+//i; 
  s/\s*transcript variant \w+//i;
  s/,\s*partial\s*\w*//;  
  s/(Acyrthosiphon pisum|Nasonia vitripennis) similar\s*//i; 
  s/similar to\s*//i;
  s/conserved expressed\s*//i; s/conserved\s*//i;  
  
  s/GROwth/growth/; # odd problems
  $_= lc($_) if(/PROTEIN/);
  # s/COIL-COILED (MYOSIN-LIKE) PROTEIN/Coiled-coil (myosin-like) protein/i;
  
  if( s/\]\]/\]/) { s/\[//; }# [Pyruvate dehydrogenase [lipoamide]] kinase
  
  # CG9915 CG9915-PB < drop duplicate
  my($nb,$nc)= m/(\w+) (\w+)\-\w+/; if($nc and $nb eq $nc) { s/ $nc\-\w+//; }  
  
  s/[,.]?\s+$//; s/^\s+//;
  s/predicted protein/$GENE_NONAME/i; 
  s/Putative uncharacterized protein/uncharacterized protein/i;
  s/hypothetical protein LOC.*/$GENE_NONAME/i; 
  if(/^(probable|putative)\s+\w+/i) { my $p=$1; s/^$p\s*//i; s/$/, putative/; }
  s/[,.]?\s+$//; s/^\s+//;

  # dang arp2 has all lowercased names; re-uppecase?  g-protein-coupled > G-protein; class [a-z] > Class [A-Z]
  ## reset all uncharacterized to "hypothetical protein" for now, change at output from quality
  ## but what of specific genes: uncharacterized protein KIAA1530-like; uncharacterized protein C4orf34 homolog
  unless(m/uncharacterized protein\s+\w+/i) { s/uncharacterized protein/$GENE_NONAME/i; }

  $_= $GENE_NONAME if( m/^(na|protein|unknown)$/i or ! m/\w/);
  return $_;
}

# Class A rhodopsinG-protein coupled receptor GPRnna8
# class a rhodopsing-protein coupled receptor gproar4

my (%casewords, %casefirst);

sub recase {
  my($na,$namehash)= @_;  # wont always work, case is context dependent
  if(ref $namehash) {
    #my $D=($debug)?'#':"";
    foreach my $n (sort keys %$namehash) {
      my @w= grep /[A-Z]/, split(/\W+/, $n); 
      if($n =~ /^[A-Z][a-z]+/) { my $w= shift @w; $casefirst{ lc($w) }= $w; }
      map{ $casewords{ lc($_) }= $_ unless(/^Protein$/i); } @w;  # $D.
      
      #? COIL-COILED PROTEIN problem
    }
  }
  
  my @uw; my $i=0;
  my @w= split( /(\W+)/, $na);
  foreach my $w (@w) { 
    my $uw= $casewords{ lc($w) };
    unless($uw or $i>0) { $uw= $casefirst{ lc($w) }; }   
    push(@uw, $uw || $w);  $i++;
    }
  return join("",@uw);
}



# add sort/count names, removing extras: -like putative homolog
sub namecount {

  my (%nac, %nav, %ncase, $n);
  while(<>) {  
    s/\s*\(\d+[CI%]*\)$//; # namepct
    s/\s*\[[^\]]+\]\s*$//; # ncbi [species]
    $_ = nameclean($_);
    s/(homolog|putative|probable|partial)\s*//g;  s/-like//g;
    # maybe drop trailing nums: family form num:  protein 1|2|3...
    s/\s+\d+$//;
    s/\s*\([^\)]+.//g;  s/\s*\[[^\]]+\]\s*$//;
    s/[,;.]*\s*$//; s/^\s+//;
    
    #? collapse ids:  CGnnnn ; GA27495; AGAP009532-PA ; ACYPI001775 ; apis_ncbi_hmm3391 ; daphnia_hxAUG25s13g261t1
    
    s/^CG\d+\b/CG-drosmel/; # -PA sometimes
    s/^G[A-K]\d+\b/drosophila gene/; # 
    s/^AGAP\d+-P.\b/anopheles gene/;
    s/^ACYPI\d+\b/aphid gene/;
    s/^(tribolium|apis|daphnia|pediculus|drosmel)_[\w-]+/$1 gene/;
    
    my $iste= isTEname($_);
    my $lna= lc($_);  
    $lna = "zzte_".$lna if($iste);
    $nac{$lna} ++;  $n++;
    $nav{$lna}= $_;
    $ncase{$_}++;
  }
  
  if($recase) { 
    recase("xxx", \%ncase);
    while( my($k,$v) = each(%nav) ) { $nav{$k}= recase($v); }
  }
  
  foreach my $na (sort keys %nac) { 
    my $tep= ($na =~ /^zzte_/) ? "TE:" : "";
    print  $nac{$na}, "\t", $tep.$nav{$na}, "\n";
  }
  
}



sub vmap{ 
  my($h,$v)=@_; my %r=(); 
  for my $i (0..$#$h){ 
    # some cleanup should be elsewhere: no [=;] in fields; drop null='.'
    my $w=$$v[$i]; $w="" if($w eq ".");  $w =~ s/=/:/g; $w =~ s/;/,/g;
    $r{$$h[$i]}= $w; 
    } 
  return %r; 
}

sub pct { my $p=shift; return int(0.5 + 100*$p); }
sub pqual { my $p=shift; return(($p >= 66)? "Strong" : ($p >= 33) ? "Medium" : ($p >= 5) ? "Weak" : "None"); }


sub isTEname {
  my $na= shift;
  foreach my $te (@TEnames) { return 1 if($na =~ /$te/i); }
  return 0;
}

# BEGIN { @TEnames= map{ s/^TE\s+\d+\s*//; } split"\n", TEnames(); }
sub getTEnames
{
# likely transposon gene

my @TEnames= map{ s/^TE\s+\d+\s*//; $_; } 
  split "\n", <<"EOTE";
TE  212 gag-pol polyprotein precursor
TE  136 Transposase
TE  133 Endonuclease-reverse transcriptase
TE  131 transposase domain-containing protein
TE  127 reverse transcriptase/RNaseH
TE  114 polyprotein of retroviral origin
TE  105 ORF2 for putative reverse transcriptase
TE  102 has homology with reverse transcriptase, 5 end of coding region for ORF2 undetermined
TE  100 Reverse transcriptase
TE  100 gag-pol protein, putative
TE   92 F-element protein
TE   87 120.7 kDa protein in NOF-FB transposable element
TE   76 endonuclease and reverse transcriptase-like protein
TE   73 Retrovirus-related Pol polyprotein from transposon TNT 1-94
TE   63 transposase
TE   57 Transposable element P transposase
TE   56 Transposable element Tc3 transposase
TE   50 Endonuclease and reverse transcriptase-like protein
TE   44 PiggyBac transposable element-derived protein=
TE   43 RNA-directed DNA polymerase homolog (R1)
TE   42 AC9 transposase, putative
TE   41 Pol protein
TE   39 Probable RNA-directed DNA polymerase from transposon X-element
TE   34 Gag protein
TE   32 Polyprotein
TE   31 pol-like protein
TE   29 pogo family transposase
TE   27 Enzymatic polyprotein, putative
TE   25 BEL12_AG transposon polyprotein
TE   25 Tigger transposable element-derived protein
TE   24 Tigger transposable element-derived protein
TE   24 non-LTR retrotransposon R1Bmks ORF2 protein
TE   22 gag-pol polyprotein precursor
TE   22 mariner transposase
TE   21 Yabusame-2
TE   20 AC transposase
TE   20 Nucleic-acid-binding protein from transposon X-element
TE   16 ORF2, putative
TE   16 Pol-like protein
TE   16 Tigger transposable element-derived protein
TE   16 reverse transcriptase, putative
TE   15 RNA-directed DNA polymerase from mobile element jockey
TE   15 Transposase
TE   15 enzymatic polyprotein
TE   14 Reverse transcriptase-like protein
TE   12 Reverse transcriptase
TE   11 Blastopia polyprotein
TE   11 Polyprotein of retroviral origin, putative
TE   11 uncharacterized transposon-derived protein F52C9.6, putative
TE   10 Copia protein (Gag-int-pol protein)
TE   10 DNA Pol B2 domain-containing protein
TE   10 Gag-pol polyprotein
TE   10 Transposase Tn3
TE    9 Copia protein
TE    8 Probable RNA-directed DNA polymerase from transposon BS
TE    8 reverse transcriptase
TE    7 PiggyBac transposable element-derived protein
TE    6 PiggyBac transposable element-derived protein
TE    5 Tigger transposable element-derived protein
TE    4 Retrovirus-related Pol polyprotein from transposon 297
TE    4 Reverse transcriptase homolog protein
TE    4 transposase, putative
TE    3 Protease, reverse transcriptase, ribonuclease H, integrase
TE    3 Tigger transposable element-derived protein
TE    3 Transposon-derived Buster3 transposase-like protein
TE    3 polyprotein of viral origin
TE    3 reverse transcriptase homolog (LOC662692)
TE    2 115 kDa protein in type-1 retrotransposable element R1DM, putative
TE    2 Retrovirus-related Pol polyprotein from transposon 412
TE    2 Transposase insG for insertion sequence element IS4
TE    2 piggybac transposable element, putative
TE    1 120.7 kDa protein in NOF-FB transposable element-like
TE    1 Endonuclease/reverse transcriptase
TE    1 Pol polyprotein
TE    1 Protease and reverse transcriptase-like protein
TE    1 Retrotransposable element Tf2 155 kDa protein type
TE    1 Retrovirus-related Pol polyprotein from type-1 retrotransposable element R1
TE    1 Reverse transcriptase/RNaseH
TE    1 Transposon, putative
TE    1 Uncharacterized transposase-like protein HI1328.1
TE    1 bel12_ag transposon polyprotein
TE    1 enzymatic polyprotein, Endonuclease, Reverse transcriptase
TE    1 gag-pol polyprotein precursor, hypothetical protein
TE    1 pol polyprotein
TE    1 polyprotein
EOTE

warn "# TEnames: @TEnames[0..9,-1] \n" if($debug);
return @TEnames;
}

__END__


cat aphid2pub8d-{attr,homolog,express,gequal,name}.tab | env table=0 ../bestgenes_puban.pl\
> aphid2pub8d.attr.gfi

cat aphid2pub8d.attr.gfi ../bestgenes.DGILpub8d3.gff | perl -ne\
'if(/^ID=([^;]+)/){ $d=$1; chomp; $at{$d}=$_; } \
else{ if(/\tmRNA/ and /ID=([^;]+)/){ $at=$at{$1}; s/ID=.*$/$at/ if($at); }  print; }' \
> aphid2pub8dan1.gff 


# this gene attr redo needs fairly complex perl... put in bestgenes_update or one-off ?

.. pub8d annot table contents
ID, geneidm oid, osrc [keep/drop EG2ap?] : change to oid=osrc:oid?
equivalence: pub8d.ident.tab for acypi1, ncbi2ref, gnomon2?, valrefseq?, pasa2_aphid3??
Name : from unipName/acypiName/other?
       use best from unipNAME/ID,Arp2NAME,acypi1NAME ; dont keep all these names but keep IDs, Dbxref= for equivalents, homolog= for blastp
  
homolog=,paralog=,pHOBEST (improved) : from aphid2pub8d-uparphumbac.bltab,aphid2pub8d-arp5hum.bltab
  ^^ now have 2 homolog: arp5 and uniparp; choose which? want unip id, desc: use both?
express= or pEXPRESS= (tr-rna/est) : from aphid2pub8d-rnastats.tab: nmap1, npairok, len_chr=len_cover,
      maybe value express=len_cover/len_tr %,nread [rna + est?]

ARP2name/id ?
flags= or reame qualities=
   Express/Express+ = hasexpress, goodexpress,
   Ortholog=hasho OR Paralog=hasparalog (from pBEST)
   Fusion_maybe=join, Userchoice=must, 
   >> these have values, not flags= ?
   Protqual=complete/partial (aaflag)/poor_coding=ncRNA?, Transposon?
   Intronqual=, Matequal=
   
pub8d annots (- protein, -best_rseq, -best_xxx)
   : keep cxlen, aalen, inexon, osrc, scoresum=col5, scorevec? : convert inqual,matequal to Flags
   : keep as flags:  must, join,  inqual, maqual,
   : use scorevec for TEgene, calc ncRNA flag from cxlen low and/or aalen low


gff mrna annot fields:
   ID=xxx;
   Name=
   oid=osrc:oid;
   aasize=  [aalen]
   cdssize= [cxlen]
   introns= [inexon]
   alttr=1 (altnum?)
   express=pctcov,nreads(rna+est[ER])
   ortholog=pctbits,bits,unip:gene,arp:gene
   paralog=bits,pctbits,selfid
   quality=   Evidence:Strong|Medium|Partial, < compute from all else how?
              Homology:Ortholog|Paralog|None, # strong/weak modifier?
              Express:Strong|Weak|None,
              Protein:Complete|Partial|Poor/Noncode
              Intron:inqual
              Mated:maqual
              Expertchoice,Fusion_maybe,Transposon, (others?)
   scorevec=
   protein= ?? keep in pub.gff or not
   
where Strong == > 66/75/90% depending? ; Medium == 33%..66/75%,  Weak == 5/10% .. 33%,  None < 5/10%
#-------------------

Guide to pea aphid Evigene annot.txt columns and GFF mRNA  attributes:
  transcriptID    (ID in gff mRNA)
  geneID          (gene in gff mRNA, is Parent= to mRNA)
  isoform   : alternate transcript number if > 1, matches ID suffix (t2,t3...)
  quality   : list of quality values for Expression Homology Intron Mate-pairing, Protein,         
  aaSize    : protein aa length, percent of transcript
  cdsSize   : cds length / transcript length
  Name      : homology-derived gene name, UniProt arthropods and related databases
  Dbxref    : cross reference gene IDs to AphidBase v1, NCBI RefSeq v2
  express   : expressed span as percet of transcript, and read count for EST, RNA-seq
  ortholog  : protein orthology percent identity, bit score, and protein IDs
  paralog   : protein paralogy percent identity, bit score and gene ID
  intron    : evidence introns from expression / model introns
  location  : genome location
  oid       : original model ID
  chimera   : validated split or chimeric model from ACYPI v1 gene, has 2 locations (and 2 transcript IDs)
  score     : evidence score sum
  scorevec  : evidence score vector

Quality notes:
  Values are generally Strong/Medium/Weak/None
  Homology:  Ortholog if best match is other species, Paralog for this species
  Protein:  curated_complete indicates curated by expert, including chimera split ACYPI genes
            and that protein cannot be computed from genome sequence.
  Intron: and Mated: (mate pairing) qualities include perfect/complete for all exons supported in gene,
          good, poor, none : levels of intron, mate pair quality

Other field notes:
  Dbxref  = gene cross reference, includes percent equivalence, and "I" or "C" flag.
            I = identical model, C = >= 90% coding sequence identity
  
  chimera = includes location of other split part, and computed gene model that matches part
    ID=acyp2eg0037508t1 chimera=1,Scaffold298:481226-487632:+,acyp2eg0018229t1,complete
    ID=acyp2eg0037509t1 chimera=2,Scaffold298:618846-621340:+,acyp2eg0018215t1,complete
    These should/will have gene records added to show part equivalence.
    
  
  scorevec fields are defined in top of GFF file and used to make total gene model score, using weighted values
  ##gff-version 3
  #program: overbestgenes, selection of best gene set by evidence scores
  #scoretype: homolog:9,paralog:1,ref:2,est:3,pro:3,rseq:2,intr:20,nintron:40,inqual:20,maqual:5,terepeat:-3,UTR:3,CDS:1
              
Sample
  transcriptID=acyp2eg0000002t1
  gene=acyp2eg0000002
  quality=Express:Strong,Homology:None,Intron:good,Mated:perfect,Protein:complete
  aaSize=982,66%
  cdsSize=2949/4496
  Name=uncharacterized protein (66%)
  Dbxref=APHIDBASE:ACYPI52644-RA,54%C,RefSeq:XM_003239978.1,83%C
  express=100%,12226r
  ortholog=4%,71.6,UniProt:E0VP94_PEDHC,Arp:pediculus_PHUM354810-PA
  intron=2/2
  oid=ars27cuf8:aphid_cuf8r27Gsc1.248.1
  score=31792
  scorevec=71,0,0,0,0,67,28,2,75,4496,0,1150,2949


Alternate transcript indicator in ID and isoform field:
  transcriptID=acyp2eg0000002t2
  gene=acyp2eg0000002
  isoform=2



gene annot table fields: as above, more details?

# fixes: change alttr=1 to ($alttr)= $id =~ m/t(\d+)/;
#  clean off inqual,maqual nums? leave only qual value
#  aalen: split out aaqual: complete/parial/.. ?
#  gene= problem for PASAgasmbl_ = gcluster_798 : 2nd gene= value  ; should be acyp2eg0000165
# add location for attr.tsv

# Not all EG2AP now; have EG2APu .. fixme ..
cat aphid2pub8e.gff | grep 'mRNA' | perl -ne'chomp; @v=split"\t"; @aa=split";",$v[-1]; \
@a= grep( /($kt)=/, @aa); push @a,"score=$v[5]","location=$v[0]:$v[3]-$v[4]:$v[6]"; \
%a=(); map{ ($k,$v)=split"=",$_; $a{$k}=$v unless($a{$k}); } @a; \
@a= map{ $a{$_}||"." } @kt; print join("\t",@a),"\n"; \
BEGIN{ @kt=qw(ID osrc oid gene alttr aalen cxlen inexon inqual maqual join must scorevec ); \
$kt=join"|",@kt; push(@kt,"score","location"); print join("\t",@kt),"\n"; } ' \
> aphid2pub8e-attr.tab

wc -l  aphid2pub8d-{attr,homolog,express,name}.tab
   41269 aphid2pub8d-attr.tab
   41530 aphid2pub8d-homolog.tab
   39271 aphid2pub8d-express.tab
   27719 aphid2pub8d-name.tab

wc -l aphid2pub8e-{attr,homolog,express,gequal,name}.tab
   41588 aphid2pub8e-attr.tab      << new
   41530 aphid2pub8e-homolog.tab   << needs update from aphid2pub8mustup-uparphumbac.bltab
   39271 aphid2pub8e-express.tab   << update for mustup ??
   25232 aphid2pub8e-gequal.tab    << update chimera add, mustup
   27720 aphid2pub8e-name.tab      << update mustup,chimera
   
cat aphid2pub8e-{attr,homolog,express,gequal,name}.tab | \
env debug=1 table=0 $evigene/scripts/bestgenes_puban.pl > aphid2pub8e.attr.gfi

cut -f2 aphid2pub8e-attr.tab | sort | uniq -c
   8053 AUGepi4
    966 AUGepi5
    974 AUGepir10
   7059 AUGepir16b
   5487 AUGepir2
    752 AUGepir3
   1239 AUGepir9
    171 DGILmix8do
   4979 DGILmix8du
    190 acypi
    180 acypi1split2
    963 ars17trinity
   5520 ars27cuf8
     87 mustaltmodel
    193 mustbest            << 11 are velvet
   4741 pasa2_aphid3
     33 ref_aphid2          << from where? is this ncbi refseq1 map to asm2? NO, ncbi2
                ref:NM_001161949.1   ref:NM_001163077.1   

#.......


Problems:

  Long introns, > 20kb, a few >100kb, > 35 genes span over 250kb (more than bee, but same ballpark)

  False UTRs   : extended into next gene, or includes introns, sometimes many utr-exons.
                These are areas of high expression, joined to gene ends when should not be, 
                or coding section broken artifactually to non-coding (artifactually);
                e.g. commonest in est/rna-assemblies by PASA, cufflinks

  Chimera/split genes from version 1: problem to process; 1000 computed but <100 validated,
                  some match alternate models.
            include a few well known genes like dicer-1, maleless, sex-determining fem-1
            ACYPI000122|ACYPI006952|ACYPI005290|ACYPI003167|ACYPI006652
            
  Check/compare genes in ncbi2, acypi1 not in evigene/acyp2eg,
    recover any w/ strong evidence


Finalization:
  assemble ESTs and rna-part assemblies together to fuller expressed gene models (PASA), w/ checks
  check/validate proteins
  match to acypi1, ncbi refseq2 (and gnomon2) genes
  annotate for public use
  quality checks: prot homology, expression via read mapping, 
          gene to read-intron and mate-pair agreement
  
Validation work: 
  check transcript expression by re-map rna-reads, ests to model transcripts;
  check protein orthology/paralogy with current uniprot arthropod genes + bacteria, human
  gene quality measures
  transposon gene assessment (bact/TE gene homology).
  [ clean rna-est assemblies for use : protein section problems ]
  

#.......... aphid2pub8e attr

cat aphid2pub8e.gff | grep mRNA | cut -f9 | perl -ne'@a=split ";"; map{ ( $k,$v ) =split"=" ; print "$k\n" ; }
 @a ; ' | sort | uniq -c
  41587 ID
    194 Name     <<< from osrc=acypi, many GLEAN, ignore have other name
    180 Target    # from chimera
  41587 aalen
   5001 alttr
   1805 best_ref
  26670 best_rseq
  13653 checkprot
     90 chim1     # chimera, split loc : keep?
     90 chim2     # ditto
    180 chimera   # keep
    180 cov       # chimara gmap, drop
  36204 cxlen
      6 db_xref   #?
    180 desc      << from osrc=acypi1split2 chimera, acypi1, save non-GLEAN
    158 egover    << keep, chimera
  47415 gene
  26586 homolog
  35981 inexon
  36001 inqual
    152 join
  14305 maqual
    180 match     # chimera, drop
    656 must
    280 mustord
    180 nexon     # chimera drop
      3 note     << what?  these are only? ncbi2 from mustbest
  41587 oid
  41576 osrc
  32992 pHOBEST
  30469 paralog
  24667 pct_support
    180 pid       # drop, chimera
  41497 protein
    180 qlen      # drop chimera
  35981 scoresum
  35981 scorevec

#......

# name FIXME : no [=;] in field; drop name= parent=\w+

# change conserved hypothetical protein > Uncharacterized; or dont change hypoth.
# Naming changes: see daph magna: 
# conserved prot = homology but no func, expressed prot = no sig homol, but express
# .. use "similar to" or "-like" for weak homolog (or paralog stronger)
# .. see uniprot clusters: cluster name better for cases of weaker homol or poor species gene name

cat aphid2pub8d-*.detab | sort -k1,1 -k3,3 | perl -ne\
'chomp; ($g,$q,$na,$dx)=split"\t"; ($t)= $na=~m/^(\w+:)/; $dx=$t.$dx; if($lg ne $g) { \
putg() if($lg); $lna=$na; $lna=~s/^\w+://; $lq=$q; } else { $dx="$lx,$dx" if($lx); } \
($lg,$xq,$lx)=($g,$q,$dx); sub putg { $lx=~s/^,+//; $lx=~s/.Uni/Uni/; $lna||="na"; \
$lna =~ s/conserved hypothetical/conserved/;  $lna =~ s/;/,/g; $lna=~s/=/:/g; \
$lna =~ s/Acyrthosiphon pisum similar //; \
$lna =~ s/predicted protein/conserved protein/; \
$lna =~ s/^hypothetical protein\s*$/Uncharacterized protein/; \
$lq=~s/\w+=//; $lq=~s/,.*//; $lq=~s/([CI])$/\%${1}/; print join("\t",$lg,$lna,$lq,$lx),"\n"; }\
END{putg();} BEGIN{print join("\t",qw(gene_id name namepct nameref )),"\n"; }' \
>  aphid2pub8d-name.tab

#.........
UniP cluster names, from
ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref50/uniref50.xml.gz

gzcat uniref50.xml.gz | egrep '^<(entry id=|name.Cluster|dbReference type=)' |\
perl -ne'END{putc()} if(/^<entry id=.(\w+)/){ putc(); $cid=$1; } \
elsif(/^<name>Cluster:\s*([^<\n]+)/) { $cna=$1; } \
elsif(/^<dbReference/){ ($dx)=m/id=.([^">]+)/; push @dx,$dx if($dx); } \
sub putc { print join("\t",$cid,$cna,@dx),"\n" if($cid); $cid=$cna=""; @dx=(); }' \

<entry id="UniRef50_Q9LGK8" updated="2011-05-03">
<name>Cluster: Os01g0160700 protein</name>
..
<dbReference type="UniProtKB ID" id="Q9LGK8_ORYSJ">
............


# id update:
set xtab=express
set xtab=homolog
set xtab=name
cat ../bestgenes.DGILpub8d3.oids aphid2pub8d-$xtab.tab | perl -ne\
'if(/^(acyp2eg\w+)\t(EG2ap\S+)/) {($a,$e)=($1,$2); map{ $eg{$_}=$a; } split ",",$e; } \
else { ($eg)=split; $d=$eg{$eg} || "na"; s/\t/\t$d\t/; \
while(($ep)=m/,(EG2ap\w+)/g){ $d=$eg{$ep}||"na"; s/,$ep/,$d/; } print; }'\
> aphid2pub8d-$xtab.tab2

 

head aphid2pub8d-{attr,name,homolog,express}.tab

==> aphid2pub8d-attr.tab <==
ID      osrc    oid     gene    alttr   aalen   cxlen   inexon  inqual  maqual  join    must    scorevec        score
acyp2eg0000001t1        AUGepi4 AUGepi4p1s1g2t1 acyp2eg0000001  .       80,20%,complete 243/1196        0/2     -1/none .       .       .       0,72,0,0,0,190,0,0,-1,0,0,-117,243      324
acyp2eg0000002t1        ars27cuf8       aphid_cuf8r27Gsc1.248.1 acyp2eg0000002  .       982,66%,complete        2949/4496       2/2     75/good 4496/perfect    .       .       71,0,0,0,0,67,28,2,75,4496,0,1150,2949  31792
acyp2eg0000002t2        DGILmix8du      aphid_cuf8r27Gsc1.248.1t2       acyp2eg0000002  1       853,47%,complete        .       .       .       .       .       .       .       .
acyp2eg0000003t1        AUGepir16b      AUGepir16bp1s1g6t1      acyp2eg0000003  .       232,70%,complete        699/988 0/2     -1/none .       .       .       0,0,0,0,0,128,0,0,-1,0,0,289,699        1802

==> aphid2pub8d-name.tab <== #? add col head
gene_id0  gene_id   name  namepct   nameref
EG2ap000001t1   acyp2eg0000014t1        hypothetical protein LOC100570313       91%C    Ref2:XM_003239983.1
EG2ap000002t1   acyp2eg0000012t1        Peptidyl-prolyl cis-trans isomerase,EC=5.2.1.8, GN:GM13410      72%     UniProt:B4IF57_DROSE,Arp2:ARP2_G2333,drosmel_CG9916-PA,Ref1:ACYPI003541-RA,XM_001945068,Ref2:XM_001945068.2
EG2ap000003t1   acyp2eg0000010t1        Peptidyl-prolyl cis-trans isomerase,EC=5.2.1.8  88%     UniProt:A2I3U7_MACHI,Arp2:ARP2_G2333,apis_ncbi_hmm3417,Ref1:ACYPI001656-RA,XM_001945331
EG2ap000004t1   acyp2eg0000009t1        prefoldin subunit       44%     Arp2:ARP2_G3677,tribolium_TC010530,Ref1:ACYPI009869-RA,XM_001945448,Ref2:XM_001945448.2,UniProt:D6WDY6_TRICA
EG2ap000005t1   acyp2eg0000007t1        Replication factor C subunit 2, GN:RFC2 60%     UniProt:C1BUQ0_9MAXI,Arp2:ARP2_G4095,apis_ncbi_hmm20068,Ref1:ACYPI007996-RA,XM_001945736,Ref2:XM_001945736.2,XM_003239981.1
EG2ap000006t1   acyp2eg0000008t1        Uncharacterized protein         16%     Arp2:ARP2_G2601,apis_ncbi_hmm7633,Ref2:XM_003239982.1,UniProt:E2A7L8_9HYME
EG2ap000007t1   acyp2eg0000006t1        na      52      Ref1:ACYPI52641-RA,Ref2:XM_003239980.1,UniProt:E9IQJ0_SOLIN
EG2ap000008t1   acyp2eg0000002t1        na      66%C    Ref1:ACYPI52644-RA,Ref2:XM_003239978.1
EG2ap000015t1   acyp2eg0000015t1        Decapentaplegic protein, GN:dpp 31%     UniProt:Q75WK6_9HYME,Arp2:ARP2_G3682,pediculus_PHUM346320-PA,Ref2:XM_001944112.2

==> aphid2pub8d-homolog.tab <==

gene_id ortholog        paralog pbest
EG2ap000001t1   4%,68.6,UniProt:B0W788_CULQU,Arp:drosmel_CG13127-PA     23%,384,EG2ap000019t1   23%P
EG2ap000002t1   72%,314,UniProt:B4IF57_DROSE,Arp:drosmel_CG9916-PA      66%,289,EG2ap000003t1   72%O
EG2ap000003t1   88%,296,UniProt:A2I3U7_MACHI,Arp:apis_ncbi_hmm3417      86%,289,EG2ap000002t1   88%O
EG2ap000004t1   44%,117,UniProt:D6WDY6_TRICA,Arp:tribolium_TC010530     0%,0,na 44%O
EG2ap000005t1   66%,499,UniProt:C1BUQ0_9MAXI,Arp:apis_ncbi_hmm20068     28%,213,EG2ap036107t1   66%O
EG2ap000006t1   16%,120,UniProt:E2A7L8_9HYME,Arp:apis_ncbi_hmm7633      0%,0,na 16%O
EG2ap000007t1   12%,124,UniProt:E9IQJ0_SOLIN,Arp:sp|SNR48_HUMAN 0%,0,na 12%O
EG2ap000008t1   4%,71.6,UniProt:E0VP94_PEDHC,Arp:pediculus_PHUM354810-PA        0%,0,na none

==> aphid2pub8d-express.tab <==
# attr: express= %xcov/xlen, xcover, rrna,iest;  quality=Express:Strong|Weak|None
#
gene_id xcover  est_cover       rna_cover
EG2ap000001t1   2895    965,97i 2895,2812r
EG2ap000002t1   1359    890,100i        1359,364793r
EG2ap000003t1   2366    1586,99i        2366,46165r
EG2ap000004t1   1857    1062,99i        1857,7407r


