#!/usr/bin/env perl
# bestgenes_puban.pl

# cacao3 version 2011.oct, from aphid2
# cat aphid2pub8d-{attr,homolog,express,gequal,name}.tab | env table=1 ../bestgenes_puban.pl > aphid2pub8d.attr.tbl

##aphid#scoretype: homolog:9,paralog:1,ref:2,est:3,pro:3,rseq:2,intr:20,nintron:40,inqual:20,maqual:5,terepeat:-3,UTR:3,CDS:1

# 1022.Oct.10 genes/bestgenes_of11.ba3ec.gff
# 2011.Oct.10 genes/cacao11_bestgenes.pub3e.gff
#scoretype: homolog:9,paralog:1,ovpro:1,ovrna:4,est:1,rseq:1,nintron:50,inqual:2,inerr:20,intr:0,terepeat:-2,UTR:5,CDS:1
#dropscore: *homolog:80,*paralog:149,ovpro:50,*nintron:2,inqual:20,+CDS:201
#sources: n=11: AUGepir1,AUGepir1a,AUGepir3,AUGie3,AUGpier6,AUGpier8,AUGpier8a,AUGpiern7,mar1g.mar11f,mar7g.mar11f,rna8b

## FIXMEx: formatTABLE : split all? multi-value columns with ',' and / into separate cols (for excel)

use strict;



# unused:
my $jTEscore= 10; #  index from scorevec in gene.gff  # FIXME, read from .gff:#scoretype
# unused:
my @VKEY= qw(ID osrc oid gene alttr aalen cxlen nintron join must scorevec);
# old: inexon inqual maqual  chimera chim1 chim2 egover   # 4 from chimera

# output columns
# my @atkey_aphid=qw(ID gene isoform  quality aaSize cdsSize Name Dbxref express ortholog paralog intron oid chimera score scorevec);

# Dbxref here can get long with all possible;
# .. add Uniref, or Dbxref==Uniref? OR Uniref goes into ortholog; 
# .. othergene=cons9,cirad equiv; << was Dbxref, use that for now?

# equiv1,2 // eqgene1,2 == two equiv gene sets; was Dbxref
# Dbxref = TAIR gene ID for cacao
# add class to quality; estgroup
# FIXME: isoform == 0 for no alts, == 1 if has alts?

my @HOKEY= qw( Dbxref ortholog paralog uniprot genegroup);
my @ATKEY= ( qw(ID gene isoform quality aaSize cdsSize Name oname groupname),
   @HOKEY, qw( equiv1 equiv2 intron express estgroup location oid score scorevec) );

##my @naKEYS= qw(Dbxref Name oname); # table output, empty = "na"
my @naKEYS= qw(Dbxref Name oname genegroup groupname); # table output, empty = "na"

my @QUAL_CLASS=qw(qualExpress qualHomology qualIntron qualProtein); # ensure these tags exist for all
# before qualClass

use constant { pSTRONG => 66, pMEDIUM => 33, pWEAK => 5, pTEWEAK => 19 };

# gequal.tab columns: use col header as REFtag ??
my $REF1tag = "cacaoGD09";  # CGDnnnn   //"APHIDBASE";  # ref1gene REF1tag = APHIDBASE
my $REF2tag = "cacaoTCR1";  # Tc01_t000030 .. "RefSeq";     # ref2gene REF2tag = RefSeq

    #  change col header for Microstupid excelNOT:
my %recode_key = ( ID => "transcriptID", gene => "geneID",
                   equiv1 => $REF1tag, equiv2 => $REF2tag, );

my %subcols=(); # for simpletab

my $NO_NAMECLEAN= $ENV{noclean} || 0;
my $MISSING= $ENV{missing} || 0; #"."; # or na, for TABLE only
my $formatTABLE = $ENV{table} || 0;
my $simpletab= $ENV{simtable} || 0; # convert ',' in cols to subcols
$formatTABLE=1 if($simpletab); # $simpletab=0 unless($formatTABLE);
my $debug= $ENV{debug} || 0;
my $recase= $ENV{recase} || 0;
my $namecount= $ENV{count} || 0;
my $USE_TENAME= $ENV{tename} || 0; # FIXME?

my $NAME_NONE  = "Unknown"; ##"hypothetical protein";
my $NAME_NOFUN = "Uncharacterized protein"; # same as Uniprot

# all evigene.conf options
my $NAME_NONE2 = "Unknown|Uncharacterized|Hypothetical";  
## my $NAME_UNK  = "Uncharacterized protein"; # uniprot  
my $MIN_IDENTITY = 35;  # min similar% for protein evid align/equivalence, for naming; JCVI uses 35%  
my $MIN_IDLIKE   = 15;  # low for now; 15..20 seems right;add Note with loqualname
my $MIN_CERTAIN  = 60;  # putative naming, what?



@ATKEY = grep !/^location/, @ATKEY unless($formatTABLE);
@ATKEY = grep !/^scorevec/, @ATKEY if($simpletab && $formatTABLE);

my @TEnames= ($USE_TENAME) ? getTEnames() : (); # drop this ; use annot table

# my ($tb, %gattr, @ids, @ha, @hn, @hc, @hg, @hh, @hx, @hd);
my ( %gattr, @ids);
my $didhead=0;

if($namecount) { namecount(); }
else { processtables(); }

#...................................................

sub getattr { 
  my $id= shift; 
  my $attr= $gattr{$id};  
  unless(ref $attr) { 
    $attr={}; $gattr{$id}= $attr; 
    if($formatTABLE) { map{ $$attr{$_}= "na"; } @naKEYS; }
      # $$attr{Name}= $$attr{Dbxref}= "na"; 
  }
  return $attr;
}


# ==> cacao11pub3e-class.tab <==
# gene_id class   estgroup
# Thecc1EG000001t1        good    
# Thecc1EG000002t1        good    BLP

sub processtables 
{
  my ($tb,@ha, @hn, @hc, @hg, @hh, @hx, @hd, @heg);

  while(<>) {
  chomp; my @v=split"\t"; 
  
  if(/^\W/) { next; }
  elsif(/^ID\tosrc/){ $tb="a"; @ha=@v; push @hd,@v; }  # genes-attr.tab # should have all mRNA ids
  elsif(/^gene_id\tname\t/) { $tb="n"; splice(@v,0,1); @hn=@v; push @hd, @v; } # genes-name.tab
  elsif(/^gene_id\tclass\t/) { $tb="c"; splice(@v,0,1); @hc=@v; push @hd, @v; } # genes-class.tab
  elsif(/^gene_id\tortholog\t/) { $tb="h"; splice(@v,0,1); @hh=@v; push @hd, @v; } # genes-homolog.tab
  elsif(/^gene_id\t(xcover|express)\t/) { $tb="x"; splice(@v,0,1); @hx=@v; push @hd, @v; } # genes-express.tab
  elsif(/^gene_id\testgroup/) { $tb="eg"; splice(@v,0,1); @heg=@v; push @hd, @v; } # genes-estgroup.tab
  elsif(/^gene_id\tref1gene\t/) { $tb="g"; splice(@v,0,1); @hg=@v; push @hd, @v; } # genes-gequal.tab
  elsif(/^gene/i or /gene_id\t/i) { $tb="HUH?"; warn "# SKIPPING. Dont know table input: $_\n"; }
  
  elsif($tb eq "a") {
    my %v= vmap(\@ha,\@v); 
    
    my $id= $v{ID}; push(@ids, $id);
    my $attr= getattr($id); # $gattr{$id};  unless($attr) { $attr={}; $gattr{$id}= $attr; }
    $$attr{ID}= $id;   # RECODE for Table: Damn Microstupid Excel barfs on ID column name
    $$attr{gene}= $v{gene}; # can be null; fix?  # recode for table : gene_id:

# transcriptID            aaSize  cdsSize
# Thecc1EG000001t1        249,76% 750/981
# Thecc1EG000002t1        205,62% 618/977
## aaSize cdsSize : shift %cds to cdsSize ??
# Thecc1EG000001t1        249     76%,750/981

    my $aalen= $v{aalen};  # aalen,pctcode,quality : 145,35%,curated-complete
    $aalen =~ s/,[a-z][\w-]+//; 

    my $cxlen= $v{cxlen};  # cds/tr,%pctcode  ** MISSING in some
    
    my($clen,$trlen)= $cxlen=~ m/(\d+).(\d+)/;  
    my $pcds=0;
    if($cxlen =~ s/,(\d+)%$//) { $pcds=$1; }
    if($aalen =~ s/,(\d+)%//) { $pcds=$1; }
    unless($pcds>0 or $trlen==0) { $pcds= int(0.5 + 100*$clen/$trlen); }
    $cxlen = "$pcds%,$cxlen"; #was# $aalen = "$aalen,$pcds%";
    
    $$attr{aaSize} = $aalen;
    $$attr{cdsSize}= $cxlen;
    
    # FIXME : missinc trlen, cxlen: get from where?  aalen*3 + utrfudge ?
    if($trlen<1 and $aalen>0) { $trlen= 3 * $aalen; $trlen += int(0.1*$trlen); }
    $$attr{trlen}= $trlen; #$trlen||=1; 
    
    $$attr{oid}="$v{osrc}:$v{oid}"; 
    $$attr{score}= $v{score}; 
    $$attr{location}= $v{location} if($v{location} and $formatTABLE); 
    
    # FIXME: isoform == 0 for no alts, == 1 if has alts? change main-attr.tab input
    $$attr{isoform} = 0;
    if($v{alttr}) { my($t)=$v{ID} =~ m/t(\d+)$/; $$attr{isoform} = $t||$v{alttr}; }
    
    #..replace qualIntron value from nintron (nfound/nsplices) ?
    # .. ie, not 100 for 2 found of 2 splices, but OK
    # .. Strong if nsplice >= 10 and nfound/nsplice > 0.75  or nsplice >= 4 and nf/ns == 1
    # .. else Medium if ns>=10 and nf/ns > 0.65 or nsplice >= 4 and nf/ns>0.84 or ns>1 and nf/ns == 1
    # .. else Ok/weak if (nf > 0)

#     ($$attr{qualIntron} = $v{inqual}) =~ s,^.*/,, if $v{inqual};    # old
#     ($$attr{qualMated}  = $v{maqual}) =~ s,^.*/,, if $v{maqual};    # old
    
      # note: missing nintron has no introns; should have that flagged.
    my $inqual="None";  # use none, not poor for no evidence introns
    my $nin= $v{nintron} || 0;
    my $inqual1= $v{inqual} || 0;
    
    ## fix: Thecc1EG005832t3        2/12,longerr:332434,
    if($simpletab) { $nin =~ s/,(\w*err)/-$1/g;  }
    $nin=~s/,$//;
    
    $$attr{intron}= $nin; # was inexon=in/ex
    if( $nin =~ m,(\d+)/(\d+), ) { 
      my($ni,$ns)=($1,$2); ## split"/",$nin;
      if($ns>0) {
        my $pin= int(0.5+100*$ni/$ns);
        if($ns>9) { $inqual=($pin>74)?"Strong":($pin>66)?"Medium":($pin>0)?"Weak":"None"; }
        elsif($ns>3) { $inqual=($pin>99)?"Strong":($pin>84)?"Medium":($pin>0)?"Weak":"None"; }
        elsif($ns>1) { $inqual=($pin>99)?"Medium":($pin>0)?"Weak":"None"; }          
        $$attr{intron}= "$pin%,$nin";
      } 
    }
    $$attr{qualIntron} = $inqual;
    
    ($$attr{terepeat} = $v{tere}) if $v{tere}; # transposon bases; keep or not?
    
    # ($$attr{qualProtein} = $v{aalen}) =~ s/^\d+,\d+.,// if $v{aalen};
    # ^^ add qualProtein=Poor  for ncRNA-like: utr too long or too many utr exons, an aalen < minaa
    # FIXME: add qualProtein=curated-complete?
    
    if( $v{aalen} =~ /,/ ) { # ** FIXME: do after all tables, qualHomology not read yet?
      my($aal,$aap,$aaq)= split",",$v{aalen};
      $aap =~ s/\%//;  
      
      ## FIXME: have some StrongOrtho + StrongerParalog genes. MOVE to fixClass
      if( $aap =~ /\d/ ) {#  and $$attr{qualHomology} !~ /Homology:Ortholog(Strong|Medium)/
        if($aap < 16) { $aaq= "noncode_$aaq"; } 
        elsif($aap < 34) { $aaq= "utrpoor_$aaq"; } #?? << utrpoor instead; FIXME: use -class.tab instead
        }      
      $aaq =~ s/curated-complete/curated_complete/;
      $$attr{qualProtein}= $aaq;
    }
     
    $$attr{'qualExpertchoice'} = 1 if($v{must} > 0);  # must= 77 or 69 (alt)
    # $$attr{qualFusionMaybe} = 1 if($v{'join'} =~ /\w/); # join=l1/p2vi1vi1l1/p1p2/..

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
    unless($simpletab) { #bug get this or 0 0 0 0 in place in simtab
    $$attr{scorevec}= join",", map{ s,/.*,,; $_ } @scorevec;
    }
    
    ## fixme here? set if >= minTEbases, but include other criteria : Express:Strong, ..
    #x. my $pte= $scorevec[$jTEscore] / $trlen;
    #x. $$attr{qualTransposon} = 1 if($pte > 0.33); 

    my $iste= ($v{class} and $v{class} =~ /transposon/) ? 1 : 0;  # NOT here ?? or {class} transposon?
    $$attr{qualTransposon} = 1 if($iste); 

  } elsif($tb eq "c") {
    my($id,@v)= @v;
    my %v= vmap(\@hc,\@v); 
    my $attr= getattr($id); 
    
    # @{$attr}{@hc} = @v{@hc};
    
    my $cla = $v{class}; # goes to qual.. ??  class=good,transposon,poorcoding/partof/...
    # $$attr{class}= $cla;
    
    # allow for asis class.tab: quality= from updated pub.gff
    if($cla =~ /^quality[:=]/) {
      $cla =~ s/^quality[:=]//;
      $$attr{'quality'}= $cla;  # see fixClass()
      $$attr{'fixed'}=2;    
      my($cc)= $cla=~m/^([^;,]+)/; $$attr{qualClass}=$cc||"None";
      
    } else {
     ## modify good class for evidence: strong Ho / Ex == Good, good+middle/weak == Ok, poor/part == Poor
      ## or Strong/Medium/Weak/None
      $$attr{qualClass}= $cla || "None"; # class default is blank ? == good/ok ?
    }  
    $$attr{qualTransposon} = 1 if($cla =~ /transposon/i); 

    if(my $eg= $v{estgroup}) { # may be in express.tab
    $eg=~s/B/Bean/; $eg=~s/L/Leaf/; $eg=~s/P/Pistil/; # FIXME input tab
    $$attr{estgroup}= $eg;
    }
    
  } elsif($tb eq "n") {
    my($id,@v)= @v;
    my %v= vmap(\@hn,\@v); 
    my $attr= getattr($id); 
    
    # updated cols: name    namepct nameref oname   onamepct        onameref
    # preserve oname/opct/oref if found as new field?
    
    my $na= $v{name};  # "na" is no name  ; 
    # problem for genes w/o name.tab entry; change na to NAME_NONE < fixClass does this but pub3i input bypasses it
    $na= nameclean($na);   
    
    # FIXME: input -name.tbl has bad namepct; replace with ortholog or uniprot pct
    # ** FIXME3: need to add putative from namepct here. 
    
    my $pna= $v{namepct}; $pna =~ s/[CI]//;
    if($pna) { 
      
      # name add: putative, like, Unknown? from namepct  
      my $islike = ($na =~ /\-like/)?2:0;
      my $isput  = ($na =~ / putative/)?2:0;
      my $isunk  = ($na =~ m/^($NAME_NONE2)/i) ? 2 : ($pna < $MIN_IDENTITY)? 1:0;  
      unless($isput or $isunk) { $isput= ($pna < $MIN_CERTAIN)?1:0; }
      if($isunk == 1 and $pna >= $MIN_IDLIKE) { $islike=1; $isput=0; } ##  for > 10-15% ident? or any?
      if($islike == 1 and not ($isput==2 or $na =~ /family|\blike|protein kinase/)) { $na =~ s/\s+protein//; $na .= "-like protein"; }  
      if($isput==1 and not $isunk) { $na  .= ", putative"; }
      
      $na =~ s=[,]=.=g if($simpletab); # convert (33%T) to Name2 column for simpletab ??
      
      $pna =~ s/$/\%/ unless($pna =~ m/\%/);
      my $naref= $v{nameref};   # include naref ? abbrev? should already have ID
      $pna .= substr($naref,0,1) if($naref);
      $na .= ($simpletab)? ",$pna" : " ($pna)"; #  (55%U) or (66%T) ?
      }
      
    $$attr{'Name'} = $na; 
      
    ## skip tename when namepct < 20%
    my $iste= 0;
    $iste= isTEname($na) if($pna > pTEWEAK); # FIXME: from input table ..
   
    my $ona= $v{oname};  # "na" is no name  
    $ona= nameclean($ona) if($ona); # use qualities to decide for nonames: 

    if($ona) {
      my $ndiff= namediff($ona,$na);
      unless($ndiff) { $ona="same"; } # same/ditto/?
      else {
        if(my $opct= $v{onamepct}) {
          my $oref= $v{onameref};   # include naref ? abbrev? should already have ID
          $opct .= substr($oref,0,1) if($oref);
          $ona  .= ($simpletab)? ",$opct" : " ($opct)"; #  (55%U) or (66%T) ?
          $iste= isTEname($ona) if($opct > pTEWEAK and !$iste);
          }
        }
      $$attr{'oname'} = $ona;
    }

    my $gna=$v{'ggname'};  # add 4dec2011
    if( $gna and $gna ne "na") { # dont add if same as Name, oname?
      my $ndiff= namediff($gna,$na);
      if($ndiff) { $ndiff= namediff($gna,$ona); }
      $gna =~ s=[,]=.=g if($simpletab);
      if($ndiff or $formatTABLE) { $$attr{'groupname'} = $gna;  } # homol has genegroup source
    }
    
    $$attr{'qualTransposon'} = 1 if($iste); 
    
    
  } elsif($tb eq "g") {
    my($id,@v)= @v;
    my %v= vmap(\@hg,\@v); 
    my $attr= getattr($id); 
    my $na = $v{ref1gene}; # maybe na REF1tag = APHIDBASE
    my $na2= $v{ref2gene}; # maybe na REF2tag = RefSeq
# gene_id ref1gene        ref2gene
# Thecc1EG000001t1        CGD0000036/78.83        Tc01_t000030/10.23,Tc01_t000020/10.15

    ## fixme for simpletab
    if($simpletab) { $na =~ s/,/+/g; $na2 =~ s/,/+/g; }
    
    $$attr{equiv1}= $na;   # print colname = REF1tag
    $$attr{equiv2}= $na2;  # print colname = REF2tag
    
#   if(0) { # NOT for cacao...
#     # na == ref1 = aphidbase now
#     my(@na2,$r2,$r3); # refseq,gnomon now, preserve best of each, in order (1st == best)
#     foreach (split",",$na2) {
#       if(/:/) { push @na2,$_ unless($r2++); }
#       elsif(/\w/ and $_ ne "na") { push @na2, "$REF2tag:$_" unless($r3++); }
#     }
#     
#     # new equiv format: CGD0000016/C99.77 ; CGD0000034/I100 ; Tc01_t000070/55.44 == /cds.exon percent
#     # .. use max(cds,exon)? both? only cds?
#     map{ s/,.*//; my($n,$t,$p)=m,(.+)/([CI]?)([\d\.]+),; $_="$n,$p\%$t" if($p); } ($na, @na2);
#     $na  = ($na eq "na") ? "" : ($na =~ /\w/) ? "$REF1tag:$na" : "";    
#     $na2 = join ",", @na2;
#     $$attr{Dbxref}= ($na and $na2) ? "$na,$na2" : $na.$na2;
#   }
    
  } elsif($tb eq "h") {
    my ($id,@v)= @v;
    my %v= vmap(\@hh,\@v); # ortholog,paralog,uniref now ?
    my $attr= getattr($id); 
     map{ $v{$_}="" if($v{$_} eq "na") } (@HOKEY);
   
    # Dbxref: drop prefix poptr/TAIR
    $v{Dbxref} =~ s,^\w+/(\w+:),$1,;
    if($v{genegroup}){ $v{genegroup} =~ s/^(\d+%),//; }
     
    ## fixme: ortholog=0%,0,UniProt:na,Arp:na << drop if 0 bits, and/or "na"
    $$attr{ortholog} = $v{ortholog}; # ($v{ortholog} =~ /0%,0,/) ? "" :
    $$attr{paralog}  = $v{paralog}; # ($v{paralog} =~ /0%,0,/ ) ? "" :
    $$attr{uniprot}  = $v{uniprot};  # FIXME tag names
    $$attr{Dbxref}   = $v{Dbxref};   # FIXME tag names
    $$attr{genegroup} = $v{genegroup};   # FIXME tag names
    
    # Dbxref,uniprot: output only id; or id first;  input now= nn%,bb,ID : do for all 4?
    #old#map{ $$attr{$_} =~ s/(\d+%),([\d\.]+),(.+)/$3#$1/;  $$attr{$_} =~ s/,.*//; $$attr{$_} =~ s/#/,/; } 
    map{ $$attr{$_} =~ s/(\d+%),\d[^,]*,(.+)/$2#$1/;  $$attr{$_} =~ s/,.*//; $$attr{$_} =~ s/#/,/; } 
        qw(ortholog paralog Dbxref uniprot);
    # map{ $$attr{$_}  =~ s/^(\d+%),(\d+),(.+)/$3,$1,$2/; } qw(ortholog paralog Dbxref uniprot);
    
    if( $$attr{uniprot} and not $formatTABLE) { # merge Dbxref for .gff
      map{ $$attr{$_} =~ s=,=/=g; } qw(Dbxref uniprot);
      $$attr{Dbxref} .= "," if($$attr{Dbxref});
      $$attr{Dbxref} .= delete $$attr{uniprot};
    }
    if($formatTABLE) { map{ $$attr{$_}="na" unless( $$attr{$_} ) } @HOKEY; }
    
    my $pho= $v{pbest};    
    my $qho= ($pho =~ /O/)?"Ortholog":($pho =~ /P/)?"Paralog":"";
    $pho =~ s/\D+//;
    $qho .= pqual($pho);
    $$attr{qualHomology} = $qho;

      ## FIXME: have some StrongOrtho + StrongerParalog genes. MOVED to fixClass
#     if($qho =~ /Ortholog(Strong|Medium)/ and $$attr{qualProtein} =~ /poor|noncod/i) {
#       $$attr{qualProtein} =~ s/(poor|noncode)[_]?//i;
#       }

  } elsif($tb eq "eg") {
    my($id,@v)= @v;
    my %v= vmap(\@heg,\@v);
    if(my $eg= $v{estgroup}) { # may be in express.tab
      $eg=~s/B/Bean/; $eg=~s/L/Leaf/; $eg=~s/P/Pistil/; # FIXME input tab
      my $attr= getattr($id);
      $$attr{estgroup}= $eg;
    }
  } elsif($tb eq "x") {
    my($id,@v)= @v;
    my %v= vmap(\@hx,\@v);
    my $attr= getattr($id);
    
    my $trlen= $$attr{trlen}||1; # ** should compute after all tables read
    
    if( my $exp= $v{express} ) {
      $exp =~ s/^express[:=]//; 
      my($px,$pasm)= $exp =~ m/(\d+)%,eq:([\d\.]+)/;
      $$attr{express} = ($simpletab) ? "$px%,$pasm" : $exp;  
      my $q= pqual( _max($px,$pasm) ); 
      $$attr{qualExpress} = $q;
    } else {
      # express: this is null value : 1       0       1,0r    
      # gives: express=0%,0r; .. is this right or should be no value?
      my $px= pct($v{xcover} / $trlen);  # *** FIXME for missing trlen ..
      my $ec= $v{est_cover}; # cov est
      my $rc= $v{rna_cover}; # cov rseq
      my $pasm= $v{xasm_equal}; # est/rna assembly equiv
      
      $$attr{express} = ($simpletab) ? "$px%,$pasm" : "$px\%,eq:$pasm";  
      my $q= pqual( _max($px,$pasm) ); 
      $$attr{qualExpress} = $q;
      }
    
    if(my $eg= $v{estgroup}) { # may be in express.tab
      $eg=~s/B/Bean/; $eg=~s/L/Leaf/; $eg=~s/P/Pistil/; # FIXME input tab
      $$attr{estgroup}= $eg;
      }
    
  }

  } # while table input


  # foreach my $id (@ids) { fixClass( $gattr{$id}); }  #before simpletab/subcols

  #? split cols w/ x,y,z to subcolumns
  if($simpletab) {
    my @cols= (1) x scalar(@ATKEY);
    
    foreach my $id (@ids) {
      my $attr= $gattr{$id};
      fixClass( $attr);  # geneClass( $attr); ## MOVED to fixClass
      foreach my $ic (0..$#ATKEY) {
        my $k=$ATKEY[$ic];
        my $v=$$attr{$k} or next;  
        my $nc= $v =~ tr/,/,/;
        if($nc>0) { $cols[$ic]= _max($cols[$ic],1+$nc); }
        }
    }
    
    foreach my $ic (0..$#ATKEY) {
      my $k=$ATKEY[$ic];
      my $nc=$cols[$ic];
      if($nc>1) { $subcols{$k}= $nc; } # make names??
    }
  }
  
  #? as table or gff attr?
  foreach my $id (@ids) {
    putt( $id);
  }

}

sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }

sub fixClass
{
  my ($attr) = @_;
  return unless(ref($attr));  # FIXME: Name == na needs correction still, below
  
  my $qual=""; 
  unless($$attr{'fixed'}) { #  and not $$attr{'fixed'}
    ## fixme missing Homology: Class:Poor-utrpoor,Express:Weak,Intron:none,Protein:poor_complete 
    ## ensure all qual parts exist? : none default
    ## .. but special cases: ExpertChoice, FusionMaybe
    # Class:Strong,Express:Strong,Homology:OrthologStrong,Intron:Strong,Protein:complete
    # my @QUAL_CLASS=qw( Express Homology Intron Protein);
  
    map{ $$attr{$_}="None" unless($$attr{$_}) } @QUAL_CLASS; # excludes gene qualClass
  
        ## FIXME: have some StrongOrtho + StrongerParalog genes.
    if( $$attr{qualProtein} =~ /poor|noncod/i and
      ($$attr{qualHomology} =~ /Ortholog(Strong|Medium)/ or $$attr{ortholog} >= pSTRONG) ) {
        $$attr{qualProtein} =~ s/(utrpoor|poor|noncode)[_]?//i; #? maybe leave noncode flag?
      }
  
    geneClass( $attr); # do before fix name
    
    my @ak= sort keys %$attr;
    my @qk= grep /^qual/, @ak;
    $qual= join ",", map{  (my $k=$_ ) =~ s/^qual//; "$k:$$attr{$_}" } @qk; # $didk{$_}++;
    # my $qual = $$attr{'qualClass'};
    $$attr{'quality'}= $qual;
    $$attr{'fixed'}=1;
  } else {
    $qual= $$attr{'quality'} || "None";
  }
  
  if($$attr{'fixed'} != 1) {
    #? here set/reset Name = hypothetical protein depending on quality of evidence
    #    hypoth = weak evid;  uncharacterized = strong evd;  ? expressed uncharact = strong express
    my $na= $$attr{'Name'};
    if($na eq "na" or $na !~ /\w/) { $na= $NAME_NONE; } #"hypothetical protein"
    if($na =~ /$NAME_NONE/i) {
      $na =~ s/$NAME_NONE/$NAME_NOFUN/i if($qual =~ /Class:(Strong|Medium)/); ## (Express|Homology):\w*(Strong|Medium)/i);
      $$attr{'Name'}= $na;  
    }
  }
  
}


sub geneClass
{
  my ($attr) = @_;
  my ($perf,$st,$md,$wc,$poor,$upoor,$partial,$must,$te)=(0)x10;
  my @qk= grep /^qual/, sort keys %$attr;

  # FIXME: modify Intron qual with attr{intron} == nintron, e.g. intron==2 + inqual=100 not strong but medium
  foreach my $q (@qk) {
    ##next if($c =~ /qualClass/);
    my $at= $$attr{$q};
    #old# $perf++ if(($q =~ /Intron/ and $at > 89));
    
    do{ $must++; $st++; } if($q =~ /Expertchoice/);
    $st++ if($at =~ /Strong/ ); # old: or ($q =~ /Intron/ and $at > pSTRONG));
    $md++ if($at =~ /Medium|good/); # old: or ($q =~ /Intron/ and $at > pMEDIUM));
    $wc++ if($at =~ /Weak|ok/);  # ok was good; cdsok
    # Protein:[noncode_|poor_] == utrpoor
    if($at =~ /part/){ $partial++; }  # for Protein=partial[35], Class=partof, other?
    if($at =~ /utrpoor/ or ($q =~ /Protein/ and $at =~ /poor/)) { $upoor++; }
    elsif($at =~ /partof|poor|noncod/) { $poor++ ; }# cdspoor; part|
    $te++ if($q =~ /Transposon/ or $at =~ /transposon/);
  }

  # *FIX poor should not outweigh strong ho+ex, eg Thecc1EG021095t1; esp if it is aa% < min.

  ## |transposon is separate class; drop extra Transposon:1 flag ?
  ## change gene Class values?  Excellent,VeryGood,Good,Weak/Ok,Poor,Transposon ?

  ## for perf/excellent, use more than count of Strong: use ovpro/ovrna evidence equivalence >= 90?
  ## perfect == ovpro > 89 and ovrna > 89 == close match to two primary evidence models.
  
  #?Not now# $te=0 if($$attr{express} >= 33); # cancel TE class for express above weak threshol
  
  ## cancel poor if high orthology or what?
  ## add Incomplete/Partial class instead of Poor?  esp for Strong/Medium express/homol/intron
  
  $poor=0 if($must>0); # but save te?
  if($upoor) { if($st>0) { $st--; $md++; } elsif($md>0) { $md--; $wc++; } else { $wc++; } }
  $perf++ if($st>2 and $poor==0 and $te==0);
  
  my $cl= ($perf>1)? "VeryStrong" #? "Perfect" 
    : ($te>0)?"Transposon"
    : ($st>0 and $poor==0)?"Strong"
    : ($st+$md>0 and $poor==0)?"Medium"
    : ($poor>0)?"Poor"
    : ($wc>0)?"Weak"
    :"None";
 
  $cl.="Partial" if($partial>0 and $cl =~ /Strong|Medium|Weak/); #?
  
  my $ocl= $$attr{qualClass};
  if($cl =~ /Transposon/) { $ocl=""; delete $$attr{qualTransposon}; }
  $ocl=~ s/(good|ok|cdsok|poorcoding|transposon|None),?//g; 
  $ocl=~ s/(cdspoor|utrpoor),?//g if($ocl=~/part/); 
  if($$attr{qualExpertchoice}) { $ocl="expertchoice,$ocl";  delete $$attr{qualExpertchoice}; }
  $ocl=~ s/,$//; $ocl=~s/,/-/g;
  $cl .= "-$ocl" if($ocl);
  $$attr{qualClass}= $cl;
}

=item geneClass

cut -f4 cacao11pub3e.attr.tbl2 | sed 's/,.*//; s/partof:.*/partof/; s/-.*//;' | sort | uniq -c
22893 Class:Strong
 5459 Class:Medium
 7501 Class:Transposon
 8601 Class:Poor     : check, drop/keep
 1323 Class:Weak     : check: looks right, express/homol/intron are none or weak.
                      some Weak of interest: Homology:ParalogWeak/Cupredoxin superfamily protein (24%T) 
  879 Class:None     : drop



 879 Class:None == these have no evidence? scorevec=0,0,0,0,0,0,0,0,0,0,te,utr,cds: drop?
      -- 200 have weak est,rseq scores; should have been dropped.
      
 548 Class:Poor
1196 Class:Poor-cdspoor
 133 Class:Poor-cdspoor-expresste
  18 Class:Poor-expresste
  35 Class:Poor-expresste-partof
1918 Class:Poor-partof
4726 Class:Poor-utrpoor   : any good ones here?
  27 Class:Poor-utrpoor-expresste


    ## modify good class for evidence: strong Ho / Ex == Good, good+middle/weak == Ok, poor/part == Poor
    ## or Strong/Medium/Weak/None
    
Class:partof:Thecc1EG000035t1,poorcoding,Express:Strong,Homology:None,Intron:100,Protein:partial  
Class:good,Express:Strong,Homology:OrthologStrong,Intron:100,Protein:complete

=cut


sub putt {
  my($g)= @_;
  my @v; my %didk;
  my $attr= $gattr{$g};

  fixClass( $attr);  # Class and Name check; skip if attr{'fixed'}
  # geneClass( $attr); ## MOVED to fixClass
  
  my @ak= sort keys %$attr;
  my @qk= grep /^qual/, @ak;  map{ $didk{$_}++; } @qk;
  #now in attr# my $qual= join ",", map{ $didk{$_}++; (my $k=$_ ) =~ s/^qual//; "$k:$$attr{$_}" } @qk;
  
#  ## MOVED to fixClass()
#  #? here set/reset Name = hypothetical protein depending on quality of evidence
#  #    hypoth = weak evid;  uncharacterized = strong evd;  ? expressed uncharact = strong express
#   my $na= $$attr{'Name'};
#   $na= $NAME_NONE  if($na eq "na" or $na !~ /\w/); #"hypothetical protein"
#   if($na =~ /$NAME_NONE/i) {
#     $na =~ s/$NAME_NONE/$NAME_NOFUN/i
#       if($qual =~ /(Express|Homology):\w*(Strong|Medium)/i);
#     $$attr{'Name'}= $na;  
#   }
   
  foreach my $k (@ATKEY) {
    my $v=$$attr{$k}; $didk{$k}++;
    my $kl= $k;
    # if($k eq "quality") { $v=$qual; } # now in attr
    
    # $v= $MISSING if($formatTABLE and not defined $v); # fixme "defined $v" should be printable $v
    $v= $MISSING if($formatTABLE and $v !~ /\S/); # fixme "defined $v" should be printable $v
    
    if($simpletab and $subcols{$k}) { 
      my $nc=$subcols{$k}; my @nv=split",",$v; 
      push @nv, $MISSING while(@nv < $nc); push @v, @nv; 
    } else { 
      push @v, (($formatTABLE) ? $v : "$kl=$v") if($formatTABLE or $v); 
    }
    }

## FIXMEx: formatTABLE : split all? multi-value columns with ',' and / into separate cols (for excel)
## need new colhead
## maybe easiest to produce this way, then reformat that after study columns (e.g. Quality has variable subcols)

  my $tab= ($formatTABLE) ? "\t" : ";";
  if($formatTABLE and 1 > $didhead++) {  
    my @colhead = map { 
      my $k= $recode_key{$_} || $_;
      if($simpletab and $subcols{$_}) { 
        my $nc= $subcols{$_}; my @kv= ($k) x $nc; 
        foreach my $i (1..$nc) { $kv[$i-1].=$i }
        @kv;
      } else { $k; }
    } @ATKEY;
    
    print join( $tab, @colhead),"\n"; 
    }

  print join( $tab, @v),"\n";
  
}



# YEATS family protein vs YEATS family protein. sym:GAS41,TAF14B 
sub namediff {
  my($na,$nb)= @_;
  # check for substantial diff 
  map { s/[,:\(\)].*//; s/\W+/ /g; $_= lc($_); } ($na,$nb);
  return ($na eq $nb or $na =~ /$nb/ or $nb =~ /$na/) ? 0 : 1;
}


## add TE protein name classing; also quality flag Transposon

sub nameclean {
  local $_= shift;
  study;
  return $NAME_NONE if($_ eq "na" or not m/\w/);
  if($NO_NAMECLEAN) { # but replace na missing
    s=[,]=.=g if($simpletab);
    return $_;
  }
  
  my($id,$na);
  #($id,$na)= split"\t",$_,2;
  #$_=$na if($id and defined($na)); # for test only?
  
  s/^Name[=:]//; #  
  s/\s+(src|ref)[=:]\S+//ig;  s/alt-transcript[=:]\S+//;

  my $sym = (s/\s*sym:([^\n;]+)//)?$1:""; # tair, save from dels
  s/\&apos[;,]/'/g;   ## FIXME:  3&apos,-5&apos, exoribonuclease
  s/;/,/g; s/=/:/g; 
  s/PREDICTED:\s*//i;  s/\(predicted\)//i;
  s/,GN:[^\s,;]+//;  # ? drop/keep EC?
  s/,EC:[^\s,;]+//;
  s/Isoform [\w-]+ of //i; s/\s*(isoform) [\w-]+//i; #|isozyme.*
  s/\s*transcript variant \w+//i;
  s/,\s*partial\s*\w*//;  
  # s/(Acyrthosiphon pisum|Nasonia vitripennis) similar\s*//i; 
  s/similar to\s*//i;
  # s/conserved expressed\s*//i; s/conserved\s*//i;  #??
  s/whole genome shotgun sequence.*//;
  s/Putative role in.*//i;
  s/\s\(Fragment\)//i;   

  ## handle these punctuation variants w/ majority vote
  s/G.protein.coupled/G-protein-coupled/i;
  s/Signal\-/Signal /i;
  
  s/GROwth/growth/; # odd problems
  $_= lc($_) if(/PROTEIN|UNKNOWN|CYTOCHROME/);
  # COIL-COILED (MYOSIN-LIKE) PROTEIN ; DOMAIN OF UNKNOWN FUNCTION 
  
  if( s/\]\]/\]/) { s/\[//; }# [Pyruvate dehydrogenase [lipoamide]] kinase
  
  # CG9915 CG9915-PB < drop duplicate
  my($nb,$nc)= m/(\w+) (\w+)\-\w+/; if($nc and $nb eq $nc) { s/ $nc\-\w+//; }  

  if(/Protein of unknown function/i){ # Kinase-related|Arabidopsis|Plant protein of unknown function
    s/,.*//; # maybe something useful before
    s/.*Protein of unknown function.*/$NAME_NONE/i; 
  }
  s/protein family/family protein/i; # to uniprot syntax
  if(s/^protein //i) { s/\s*$/ protein/; } #?
    
  s/[,.]?\s+$//; s/^\s+//; s/^\W(\w)/$1/; 
  s/predicted protein/$NAME_NONE/i; 
  s/unknown protein.*/$NAME_NONE/i;
  s/Putative uncharacterized protein/$NAME_NOFUN/i;
  s/hypothetical protein LOC.*/$NAME_NONE/i; 
  if(s/\-like//) { s/$/, putative/; }
  if(/^(probable|putative)\s+\w+/i) { my $p=$1; s/^$p\s*//i; s/$/, putative/; }
  s/[,.]?\s+$//; s/^\s+//;
  
  unless(m/uncharacterized protein\s+\w+/i) { s/uncharacterized protein/$NAME_NONE/i; }


  $_= $NAME_NONE if( m/^(na|protein|unknown)$/i or ! m/\w/);
  $_.=" sym:$sym" if($sym);
  #$_="$id\t$_" if($id);
  s=[,]=.=g if($simpletab);
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
    s/(homolog|putative|probable|partial|fragment)\s*//g;  s/-like//g;
    # maybe drop trailing nums: family form num:  protein 1|2|3...
    s/\s+\d+$//;
    s/\s*\([^\)]+.//g;  s/\s*\[[^\]]+\]\s*$//;
    s/[,;.]*\s*$//; s/^\s+//;
    
    #? collapse ids:  CGnnnn ; GA27495; AGAP009532-PA ; ACYPI001775 ; apis_ncbi_hmm3391 ; daphnia_hxAUG25s13g261t1
    
#     s/^CG\d+\b/CG-drosmel/; # -PA sometimes
#     s/^G[A-K]\d+\b/drosophila gene/; # 
#     s/^AGAP\d+-P.\b/anopheles gene/;
#     s/^ACYPI\d+\b/aphid gene/;
#     s/^(tribolium|apis|daphnia|pediculus|drosmel)_[\w-]+/$1 gene/;
    
    my $iste= isTEname($_); # FIXME: from input table
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
sub pqual { my $p=shift; 
return(($p >= pSTRONG)? "Strong" : ($p >= pMEDIUM) ? "Medium" : ($p >= pWEAK) ? "Weak" : "None"); }


sub isTEname {
  my $na= shift;
  return 0 unless($USE_TENAME);
  foreach my $te (@TEnames) { return 1 if($na =~ /$te/i); }
  return 0;
}


=item TE vs Intron 

** TEname should use namepct to decide if valid.  EG. 9% homology to name is too little.
flamingo2.% cat cacao11_bestgenes.pub3e.tename.tab | egrep -v 'Intron:(Weak|none)' 

Thecc1EG004932t1        Class:Transposon,Intron:Strong,Protein:complete,Transposon:1    MuDR family transposase (9%T)   Unknown (44%U)
Thecc1EG006345t1        Class:Transposon,Intron:Medium,Protein:complete,Transposon:1    Mutator transposase, putative (57%U)    MuDR family transposase (13%T)
Thecc1EG007962t1        Class:Transposon,Intron:Medium,Protein:complete,Transposon:1    Micronuclear linker histone polyprotein, putative (21%U)        na
Thecc1EG012806t1        Class:Transposon,Intron:Medium,Protein:complete,Transposon:1    Transposase tnp2 (Fragment) (66%U)      na
Thecc1EG016057t1        Class:Transposon,Intron:Strong,Protein:complete,Transposon:1    hAT transposon superfamily (32%T)       Unknown (Fragment) (82%U)
Thecc1EG019572t1        Class:Transposon,Intron:Strong,Protein:complete,Transposon:1    Ac transposase protein, putative (24%U) BED zinc finger ,hAT family dimerisation domain (22%T)
Thecc1EG021500t1        Class:Transposon,Intron:Strong,Protein:complete,Transposon:1    MuDR family transposase (10%T)  Unknown (67%U)
Thecc1EG021501t1        Class:Transposon,Intron:Strong,Protein:complete,Transposon:1    MuDR family transposase (9%T)   Unknown (47%U)
Thecc1EG026266t1        Class:Transposon,Intron:Medium,Protein:complete,Transposon:1    Mutator transposase, putative (80%U)    MuDR family transposase (11%T)
Thecc1EG034299t1        Class:Transposon,Intron:Strong,Protein:poor_complete,Transposon:1       reverse transcriptase, putative (43%U)  na
Thecc1EG035107t1        Class:Transposon,Intron:Strong,Protein:complete,Transposon:1    Mutator transposase, putative (57%U)    MuDR family transposase (12%T)

=cut

sub getTEnames
{
# likely transposon gene

# *** Integrase is not always TE:  Thecc1EG000844t1 Integrase-type DNA-binding superfamily 
# .. MuDR family transposase not always TE gene.
# .. TEname class should be weaker, express>33 overrides; other evid? ..

my @TEnames= map{ s/^\d+\s*//; $_; } grep /^\d/, split "\n", <<'EOTE';
999 Retrotransposon
999 transposon
999 Ty3.gypsy
999 Ty1.copia
999 Gag.pro
999 Gag.pol
999 RNA.directed DNA polymerase
999 Polyprotein
999 Transposase
999 Transposable element
999 reverse transcriptase
999 [\b_]Copia[\b_]
999 [\b_]Gypsy[\b_]
EOTE
# not? 999 \bIntegrase\b

warn "# TEnames: @TEnames[0..9,-1] \n" if($debug);
return @TEnames;
}

# ## Transposon names from homology
# 685	Retrotransposon protein
# 341	Retrotransposon protein
# 153	Retrotransposon protein
# 139	Retrotransposon protein
# 71	Retrotransposon gag protein
# 67	Retrotransposon protein
# 271	RNA-directed DNA polymerase
# # ^?? is this always TE
# 341 Ty3-gypsy subclass
# 67  Ty3-gypsy sub-class
# 153 Ty1-copia subclass
# 39	Ty3-gypsy retrotransposon protein
# 184	Gag-pro
# 134	Gag-Pol polyprotein
# 112	Gag/pol protein
# 46	Gag-pol polymerase
# 97	Gag protease polyprotein
# 16	Gag-protease-integrase-RT-RNaseH polyprotein
# 121	Polyprotein
# 75	Copia LTR rider
# 27	Copia retrotransposable element
# 41	copia-type pol polyprotein
# 75	Integrase, catalytic region
# 70	Integrase
# 66	integrase, identical
# 69	CCHC-type integrase
# 29	Integrase core domain
# 61	Transposase
# 53	retroelement pol polyprotein
# 42	MuDRA transposase
# 35	LINE-type retrotransposon LIb DNA
# 25	Non-LTR retroelement reverse transcriptase
# 4   non-LTR retrolelement reverse transcriptase
# 22	Mutator transposase
# 17	Retrovirus-related Pol polyprotein from transposon TNT 1-94
# 16	Transposase tnp2
# 11	transposon protein
# 15	Reverse transcriptase


# ##other Transposon gene names for TEgenes identified by repbase transposons
# 17      HAT family dimerisation domain containing protein
# 10      AT hook motif-containing protein
# 6       TdcA1-ORF1-ORF2 protein
# 5       PIF transposase
# 5       RNase H family protein
# 5       polyprotein, related
# 5       transposon protein
# 4       Reverse transcriptase
# 4       retroelement polyprotein
# 3       Pol
# 3       Pol polyprotein
# 3       RNA-directed DNA polymerase, Ribonuclease H, Endonuclease/exonuclease/phosphatase
# 3       Retrotransposon Tto1 DNA
# 3       Retrotransposon like protein
# 3       Uncharacterized protein At2g14040
# 2       Ac transposase
# 2       Copia pol polyprotein
# 2       Gag protein
# 2       Genome polyprotein
# 2       LINE-type retrotransposon LIb DNA, complete sequence, Insertion at the S11 site
# 2       Polyprotein, , 77260-80472
# 2       RNA-directed DNA polymerase, .., Retrotransposon gag protein
# 2       Retrotransposon protein, , Ty3-gypsy subclass, expressed
# 2       Retrotransposon protein, , unclassified, expressed
# 2       Transcription factor
# 2       Transposon protein, , CACTA, En/Spm sub-class
# 2       Zinc knuckle family protein
# 2       gag-pol polyprotein, identical
# 2       gypsy-type retrotransposon polyprotein
# 2       polyprotein, identical
# 1       Ac transposase THELMA13
# 1       Copia-type polyprotein
# 1       Gag/pol polyprotein
# 1       HAT dimerisation
# 1       HAT family dimerisation domain, 
# 1       HAT transposase
# 1       Mutator transposase
# 1       RNA-directed DNA polymerase 
# 1       Retrotransposon Opie-2
# 1       Retrotransposon protein, , Ty1-copia subclass, expressed
# 1       Retrovirus-related like polyprotein
# 1       Reverse transcriptase-beet retrotransposon
# 1       Transposable element protein, , Retrotrans_gag
# 1       Transposon protein, , Mutator sub-class
# 1       Ty1/copia-element polyprotein
# 1       Ty3-gypsy retroelement pol polyprotein
# 1       copia protein
# 1       copia-type retrotransposon protein
# 1       non-LTR retroelement reverse transcriptase
# 1       transposase protein


# BEGIN { @TEnames= map{ s/^TE\s+\d+\s*//; } split"\n", TEnames(); }
sub getTEnames_Aphid
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


