#!/usr/bin/env perl
# evgpuban_kfish2.pl 
# cut down from bestgenes_puban_kfish.pl
# extensive revision, now just for kfish2 annot tables .. generalize hopefully
# sort -k1,1 -k2,2  kf2rae5g.main.scores.tab kf2rae5g.main.{names,homolog}.kvtab kf2rae5g.main.alttr \
#  kf2rae5g.main.ov* kf2rae5g.main.rnaxtab | env sorted=1 table=1 ../evgpuban_kfish2.pl > kf2rae5g.main.attr.tbl3 

use strict;
use Getopt::Long;
use constant VERSION => '2014.01.06'; #   

use constant MUST_NOT_EXPERT => 1; #kfish2: many must=77 from asmrna x bestgenes, not Expertchoice
use constant qualProteinDropUTRflags => 1;     # Protein:utrpoor_unavailable;Protein:utrpoor_complete-utrbad 
use constant LASTID=>1;
use constant FIRSTID=>-1;

my @HOKEY= qw(ortholog paralog genegroup Dbxref); #f  uniprot 
my @ATKEY= ( qw(ID gene isoform quality aaSize cdsSize Name oname groupname),
   @HOKEY, qw( intron express mapCover location equiv1 oid score ) ); #f equiv1 equiv2 ## drop: scorevec ; ADDED mapq/mapCover
my @naKEYS= qw(Dbxref Name oname genegroup groupname); # table output, empty = "na"

my @QUAL_CLASS=qw(qualExpress qualHomology qualIntron qualProtein); # ensure these tags exist for all
# add qualMap, 

## FIXME: process_onekvset : keys for each get_ should be constant for all process_one, for defaults
my %INPUT_KEYSET= ( 
  score   => "aalen|cxlen|locus|location",
  name    => "Name|nameln|namepct|ggname|groupname",
  homolog => "homolog|paralog|orclass",
  express => "xcover|rnax|ovrna|ovest",
  refgene => "ref1gene|ref2gene",
  class   => "class",
);
my $DEFAULT_INPUTS="score|name|homolog|express";


use constant { pSTRONG => 66, pMEDIUM => 33, pWEAK => 5, pWEAKNAME => 15, pTEWEAK => 19, pUTRPOOR => 33 };
#  change qual % for qualMap: strong => 95, medium => 90, poor => 80 or split, none/weak => NOPATH

# gequal.tab columns: use col header as REFtag ??
## kfish2: gequal == kfish1 over ids
my $REF1tag = "kfish1ref"; # "ref1"; #f "nasviOGS1";  # CGDnnnn   //"APHIDBASE";  # ref1gene REF1tag = APHIDBASE
my $REF2tag = "ref2"; #f "nasviRefSeq2";  # Tc01_t000030 .. "RefSeq";     # ref2gene REF2tag = RefSeq

    #  change col header for Microstupid excelNOT:
my %recode_key = ( ID => "transcriptID", gene => "geneID",
                  equiv1 => $REF1tag, equiv2 => $REF2tag,); # 

## add getOptions()
my $debug= $ENV{debug} || 1; # DEBUGing

my $IDSORTED= $ENV{sorted}||0; # ??** REQUIRE this
my $KEYCOUNT= $ENV{kcount} || $debug;
my $KSKIP= $ENV{kskip}||""; 

## kskip='OLD\w+,aaref,align,alignx,cdsover,nover,oaaln,ocds,omid,pct_support,splicemix,svec,score,tebp,homolog,nameln,Name,Dbxref,altc,utrx'
## FIXME: Name/nameln from score table is bad, drop here or in .score.tab?
$KSKIP= $ENV{kskip} ||
  'OLD\w+,aaref,align,alignx,cdsover,nover,oaaln,ocds,omid,pct_support,splicemix,svec,tebp,altc,utrx'
  if($debug); ## add: homolog,nameln,Name,Dbxref, score == scoresum

my $MISSING= $ENV{missing} || 0; #"."; # or na, for TABLE only
my $formatTABLE = $ENV{table} || 0;
my $simpletab= $ENV{simtable} || 0; # convert ',' in cols to subcols
my $USE_TENAME= $ENV{tename}; ##x || $namecount; # FIXME?
my $UTRPOORisNONCODING= $ENV{noncoding} || 0; # new flag for utr>>cds == noncoding; maybe drop noncode class
#x nameclean.drop# my $recase= $ENV{recase} || 0;
#x my $namecount= $ENV{count} || 0;
my $NAME_NONE  = "Unknown"; ##"hypothetical protein";
my $NAME_NOFUN = "Uncharacterized protein"; # same as Uniprot

## ugpxml fields..........
my $doxml= $ENV{xml}||0;
my $date="20140107"; # fixme
my $CLADE= $ENV{species}||"";
my $GENE_TITLE=$ENV{title} || "$CLADE gene";  
my $GENE_SOURCE="euGenes";
my $GENE_TYPE="Gene_Summary"; # a SO type



my $optok= GetOptions( 
  "date=s", \$date, "title=s", \$GENE_TITLE, 
  "species|CLADE=s", \$CLADE, 
  "MISSING=s", \$MISSING, 
  "KSKIP|keyskip=s", \$KSKIP,  
  "table!", \$formatTABLE,
  "simtable|simpletable!", \$simpletab,
  "sorted_by_id!", \$IDSORTED,
  "xml|ugpxml!", \$doxml,
  "debug!", \$debug,
);

die "EvidentialGene evgpuban_kfish2.pl VERSION ",VERSION,"
Usage: evgpuban_kfish2 input.tables > genes_ann.table | genes_ugp.xml
  opts: -sorted_by_id -table|-simpletable|-xml -keyskip='$KSKIP'  ..   
" unless($optok); 


my($didhead,$didfoot)=(0) x 9;
my ( %gattr, @ids, %subcols);
if($simpletab or $doxml) { $KEYCOUNT=$debug=0; }

$KSKIP=~s/[,\s]+/|/g;
$formatTABLE=1 if($simpletab or $doxml); # $simpletab=0 unless($formatTABLE);

@ATKEY = grep !/^(location|locus)/, @ATKEY unless($formatTABLE);
@ATKEY = grep !/^scorevec/, @ATKEY if($simpletab && $formatTABLE);
my @TEnames= ($USE_TENAME) ? getTEnames() : (); # drop this ; use annot table

# revise input tables : old are header, value rows
# replace w/ key=value tables of preceding inputs:  ID homolog=xxx; ID oid=xxx aalen=zzz locus=lll score=yyy

process_kvtables();
if($doxml and not $didfoot) { print "\n</GeneSummaries>\n"; $didfoot++; } # data BUG FIX..

#OLD processtables(); 
#x if($namecount) { namecount($namecount); } else { processtables() } ## drop here
#...................................................

sub getattr { 
  my($id)= @_;
  my $attr= $gattr{$id};  
  unless(ref $attr) { 
    $attr={}; $gattr{$id}= $attr; 
    if($formatTABLE) { map{ $$attr{$_}= "na"; } @naKEYS; }
  }
  return $attr;
}

# sub hasval { my $v=shift; return ($v and $v ne "na" and $v =~ /\w/)?1:0; }
sub hasval { return (defined($_[0]) and $_[0] ne "na" and $_[0] =~ m/\w/)?1:0; } # isname

sub process_onekvset
{
  my($id, $kvhash, $dupkvhash, $lastput, $whichtab)= @_;
  return 0 unless($id or $lastput == LASTID);
  
  vclean($kvhash);
  vclean($dupkvhash); #?? use this? for some keys??
  $kvhash->{ID}=$id;  
 
  # if($debug) { print "#dONE"; for my $k (sort keys %$kvhash) { print "\t$k=",$kvhash->{$k}; } print "\n"; }
  
  unless($whichtab) {
    $whichtab||="";
    my @k=sort keys %$kvhash;
    for my $key (sort keys %INPUT_KEYSET) {
      my $kpatt= $INPUT_KEYSET{$key};
      $whichtab .= "$key|" if(grep(/^($kpatt)/,@k)); ## $key =~ m/$DEFAULT_INPUTS/ or 
      }
  }
  
  # Fixme: bad now for whichtab always gets ID even for dropped score.tab records, from 2ndary tables.
  # require other score.tab key here? all have: locus 
  get_scoretab(2,$kvhash,$dupkvhash) if($kvhash->{locus}||$kvhash->{location}||$kvhash->{oid});
    #o if($whichtab =~ /score/); 
    #oo grep{$kvhash->{$_}} qw(aalen cxlen locus location));  # default for many keys .. need allowed key list
  # separate scoretab sub for intron,exon keys?

  # if($debug) { my $attr= $gattr{$id}; print "#dATTR1=$attr"; for my $k (sort keys %$attr) { print "\t$k=",$attr->{$k}; } print "\n"; }
  #dATTR1=HASH(0x7ff86a046be8) : ok

  get_nametab(2,$kvhash,$dupkvhash)  if($whichtab =~ /name/); # if(grep{$kvhash->{$_}} qw(Name nameln namepct));
  # ^^ bad here gattr ID=Funhe2EKm000003t1 Name=Funhe2EKm000004t1 : BAD nameclean() called $_ global
  
  get_homologtab(2,$kvhash,$dupkvhash) if($whichtab =~ /homolog/); # if(grep{$kvhash->{$_}} qw(homolog orclass));
  get_expresstab(2,$kvhash,$dupkvhash) if($whichtab =~ /express/); # if(grep{$kvhash->{$_}} qw(xcover rnax ovrna ovest));  # ovrna ovest ? are precursors to xcover express.tab
  get_classtab(2,$kvhash,$dupkvhash) if($whichtab =~ /class/); # if(grep{$kvhash->{$_}} qw(class));  # ?? need to make class.tab from others
  get_refgenetab(2,$kvhash,$dupkvhash) if($whichtab =~ /refgene/); # if(grep{$kvhash->{$_}} qw(ref1gene ref2gene));  # ?? need to make class.tab from others
 ## 
 
  ## Problem case: have $id from subsid table (rnax) but not in scores/valid table ID field.
  ## .. report such as error ..
  if($debug) { my $attr= $gattr{$id}; unless($$attr{ID}) { my @k=sort keys %$attr;
    if(@k) { print "#ERR.noID,lid=$id"; for my $k (@k) { print "\t$k=",$attr->{$k}; } print "\n"; }  }
  }
  
  #FIXME: simtabFixes() if($simpletab); ## add simpltab fixes ..  need all records for this??
  if($simpletab) {
    if($lastput == LASTID) { 
      simtabFixes(); my $nput=0; foreach my $aid (@ids) { $nput+= putattr( $aid); } return $nput;   
    } else { return 0; }
  }
  
  return ($doxml)? put_ugpxml($id,$lastput,$kvhash,$dupkvhash) : putattr($id,$lastput); ## NOT ,getattr($id));
}

=item fixme simtab cols

** paralog list: chop to one only.
** Name : problem, should be 1? or 2,3 column.
** Dbxref: 1 col
** express: 2 cols, chop rnax to max val
** oid: 1 col
transcriptID    geneID  isoform 
quality1        quality2        quality3        quality4        quality5        quality6
aaSize  cdsSize1        cdsSize2        cdsSize3        
Name1   Name2   Name3   Name4   Name5   Name6   Name7   Name8   Name9   
oname1  oname2  oname3  groupname       ortholog1       ortholog2       
paralog1        paralog2        paralog3        paralog4        paralog5        paralog6        paralog7        paralog8        paralog9        paralog10
genegroup       Dbxref1 Dbxref2 Dbxref3 Dbxref4 Dbxref5 Dbxref6 Dbxref7 
intron1 intron2 
express1        express2        express3        express4        
mapCover1       mapCover2       location        
kfish1ref       oid1    oid2    oid3    score

=cut

sub simtabFixes {
  my @cols= (1) x scalar(@ATKEY);
  foreach my $id (@ids) { #?? need all @ids hre?
    my $attr= $gattr{$id}; # my $attr= getattr($id); # safer?
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
  
sub process_kvtables
{
  my(%gval,%dupkey,%kvset,%kcount,$lid);
  %gval= %kcount= %kvset= %dupkey= (); $lid="";
  my $whichtab="";
  my $atid=0;
  
  die "ERR: need env sorted=1 for input tables, sorted by ID col1\n" unless($IDSORTED);
  ## for this form, sort input by ID col1, then all attrs in one batch
  while(<>) {
    next unless(/^\w/); chomp; 
    my($id,@val)=split"\t";    
    ## print join("\t","#d",$id,@val)."\n" if($debug);
    
    if($lid ne $id and $lid) { # $IDSORTED and 
 
      unless($whichtab) { 
        # for now, assume MOST input tables, skip class unless found.
        my @k=sort keys %kvset;
        for my $key (sort keys %INPUT_KEYSET) {
          my $kpatt = $INPUT_KEYSET{$key};
          $whichtab .= "$key|" if($key =~ m/$DEFAULT_INPUTS/ or grep(/^($kpatt)/,@k)); 
          }
      }

      process_onekvset($lid,\%kvset,\%dupkey,($atid==0)?FIRSTID:0,$whichtab);  $atid++;
      for my $k (keys %kvset) { $kcount{$k}++; }
      ## delete $gattr{$lid}; #? or not?
      %kvset=(); %dupkey=();
    }
    
    if($KSKIP){ @val= grep { not m/^($KSKIP)=/ } @val; }    
    
    foreach my $kv (grep /=/, @val) { 
      my($k,$v)= split /=/,$kv,2; 
      next unless($v=~/\S/); #? vclean may set null, keep or not?
      # next if($KSKIP and $k =~ m/^($KSKIP)$/);
      ## Dup keys now are common from Split=1/2 pairs, other? dont assume 1st is best value; preclean??
      if($kvset{$k}) { $dupkey{$k}.="$v\t"; } 
      else { $kvset{$k}=$v; } # some dupkeys ok??
    } 
   
    $lid=$id;
  }
  
  #above# use constant LASTID=>1;
  process_onekvset($lid,\%kvset,\%dupkey, LASTID, $whichtab); 
  for my $k (keys %kvset) { $kcount{$k}++; }
  
  # summary outputs .. key-count, ..
  if($KEYCOUNT) { foreach my $k (sort keys %kcount) { print "#KC\t$k\t$kcount{$k}\n"; } }
}


=item FIXME express.tab

kfish2rae5d.main.express.tab :: revert to precursor ovrna,ovest + rnax.tab ??
gene_id	xcover	est_cover	rna_cover	xasm_equal	trspan
Funhe2EKm000003t1	679	679	679	100	679
Funhe2EKm000004t1	4888	3570	4888	100	4888
Funhe2EKm000005t1	4096	3075	4096	100	4096

sort $pt.{ovest,ovrna} | perl -ne \
'($pd,$ot,$ov,$xv)=split/[\s=,]/;  if($pd ne $lpd) { puts() if($lpd);\
%xv=%ov=(); } $ov{$ot}=$ov; $xv{$ot}=$xv;  $lpd=$pd; END { puts(); } \
sub puts { ($xe,$xw)=split"/",$xv{ovest}; $rp=$ov{ovrna}||0; \
$xr=int($xw*$rp/100); $mx=($xr>$xe)?$xr:$xe; \
print join("\t",$lpd,$mx,$xe,$xr,$rp,$xw)."\n"; } \
BEGIN{ print join("\t",qw(gene_id xcover est_cover rna_cover xasm_equal trspan))."\n"; } ' \
> $pt.express.tab

=cut

=item FIXME homolog.tab mix of key=,nokeyval

names/kf2rae5g.main.homolog.tab : 
  add keys to all?  orclass=orlog1  ggroup= gref=same .. two homolog= keys bad > humlog=xxx_HUMAN
Funhe2EKm000003t1	orlog1	FISH11G14131	same	homolog=358,platyfish:ENSXMAP00000008562	homolog=130,UniProt:CDX4_HUMAN
Funhe2EKm000004t1	orlog1	FISH11G998	same	homolog=1889,platyfish:ENSXMAP00000008599	homolog=887,UniProt:PGFRB_HUMAN

=cut

=item kvtables

cat kf2rae5g.main.scores.tab kf2rae5g.main.ov* | env keys=1 kskip=$kskip perl -ne \
'next unless(/^\w/); chomp; ($gid,@val)=split"\t";  @val=  grep
/=/,@val; if($kskip){ @val= grep { not m/^($kskip)=/ } @val; } map{
$nc=tr/,/,/; if($nc>2){ @t=split",",$_; $_=join",",@t[0,1,2];} } @val;
if($lval= $gval{$gid}) { unshift(@val, $lval); }
$gval{$gid}=join"\t",@val if($gid and @val);
BEGIN{$kskip=$ENV{kskip}||""; $kskip=~s/[,\s]/|/g;} END{
$dok=$ENV{keys}||0; %k=(); for $g (sort keys %gval) { unless($dok) {
print "$g\t$gval{$g}\n"; } else { $k{ID}++; map{ ($k,$v)=split"=";
$k{$k}++;} split"\t",$gval{$g}; } } for $k (sort keys %k) { print
"$k\t$k{$k}\n";} }'

=cut


=item scoretab

  want: 
    alttr  ?? not now
    tere nintron inqual
    osrc ovpro ovrna
    
  annots from mRNA.gff, some should be replaced?
../score/kf2rae5g.mainalt.rnax3.tab   ../score/kf2rae5g.main.ovest          ../score/kf2rae5g.main.human.bltab@
../score/kf2rae5g.main.ovrna          ../score/kf2rae5g.main.ovtere         ../score/kf2rae5g.main.fish11.bltab@
../score/kf2rae5g.main.ovpro          ../score/kf2rae5g.main.ovintr         ../score/kf2rae5g.main.scores.tab
  
score/kf2rae5g.main.scores.tab
Funhe2EKm000003t1       oid=Funhe2Exx11m027882t1,Funhe2Emap3m022605t1  
  altc=noclass    cov=100.0       nexon=3 aalen=214,94%,complete 
  Dbxref=CDD:200956,TrEMBL:UniRef50_Q90423,TrEMBL:HXB1B_DANRE,   
  Name=Homeobox protein Hox-B1b   nameln=74%,227/307,214 
  homolog=229,UniRef50_I3K5P6     inqual=50       nintron=2/4    
  ovpro=83,bp7tilapia:ENSONIP00000016441/83.00   
  ovrna=100,Funhe2Exx11m027882t1/I100     cxlen=663/679,97%      
  scoresum=5464   must=77 osrc=kf2x11gmap locus=Scaffold0:5067-10366:-
Funhe2EKm000004t1       oid=Funhe2Exx11m002607t1,Fungr1EG3m001115t1    
  altc=main       cov=91%,4796/5255       nexon=23       
  aalen=1103,63%,complete
  Dbxref=CDD:173624,TrEMBL:UniRef50_P35968,TrEMBL:VGFR2_HUMAN,   
  Name=Vascular endothelial growth factor receptor 2     
  nameln=85%,1156/1356,1103       homolog=1313,UniRef50_M4A2A5   
  inqual=90       nintron=40/44  
  ovpro=98,bp7tilapia:ENSONIP00000016438/C98.00  
  ovrna=100,Funhe2Exx11m002607t1/I100,Funhe2Exx11m002607t4/89.87 
  cxlen=3283/4888,67%     scoresum=30678  must=77 osrc=kf2x11gspl
  locus=Scaffold0:25799-58393:+

=cut

  # while(<$AllInputTables>)  ..
  #   chomp; my @v=split"\t"; 
  #   if(/^\W/) { next; }
  #   elsif(/^ID\tosrc/){ $tb="a"; @ha=@v; push @hd,@v; } 
  #   case $tb in 
  #     "a": get_scoretab(\@v,\@ha); ..
  #     "h": get_hotab(\@v,\@ha); ..

=item scoretab

Required/desired/created keys
  ID 
  aalen  cxlen trlen, aaSize  cdsSize  : aa/cdsSize are output, others input
  alttr : skip now, get from IDt1,t2,.. see isoform
  chimera chim1 chim2  egover  : special gmap, replace from align.tab .. in mRNA.gff ??
  class   gene 
  location : locus alias?? or rename input locus=
  inqual intron nintron : nintron,inqual are input, intron output field name.
  must oid osrc 
  score scoresum scorevec tere  

create/update:
  isoform  : from IDt2,.. or alttr=
  intron   : from nintron, inqual
  terepeat :? need this, alias of tere,qualTransposon
  qualHomology qualIntron qualProtein qualTransposon 
    qualMap << add
  qualExpertchoice << FIXME: many/most now have must= from asmrna x bestgenes, not Expertchoice tho..
          .. all are must=77 ; recode this? parse updates for Expert vs compute"Expert" ?
  qualExpress : expresstab

=item Split gene attr

  ** need to sum some of dupkey values:
    dupkey add:
      nintron, nexon, >> recalc inqual? = nintron / (2*nexon - 2)
      locus : append both
      cov-maybe, or use align.tab summed cov?
      ovpro : not used now?
      ovrna : not if 100%, not from split overlap
      
Funhe2EKm000011t1	alttr=3
Funhe2EKm000011t1	inqual=100	nintron=12/12
Funhe2EKm000011t1	inqual=100	nintron=2/2
Funhe2EKm000011t1	oid=Funhe2Exx11m006564t1	altc=althi	Split=2/2	cov=97.6	nexon=7	osrc=kf2rae5	locus=Scaffold43:402782-424715:-
Funhe2EKm000011t1	oid=Funhe2Exx11m006564t1,Funhe2E6bm006583t2	altc=althi	Split=1/2	cov=2.1	nexon=2	aalen=683,98%,complete	
  Dbxref=CDD:214567,TrEMBL:UniRef50_O60285,TrEMBL:NUAK1_HUMAN,	
  Name=NUAK family SNF1-like kinase 1	nameln=100%,689/661,683	osrc=kf2rae5	locus=Scaffold0:133590-135018:-
Funhe2EKm000011t1	ovest=70,1461/2079
Funhe2EKm000011t1	ovpro=5,bp7platyfish:ENSXMAP00000008675/5.00
Funhe2EKm000011t1	ovpro=94,bp7zfish:ENSDARP00000029561/C94.00
Funhe2EKm000011t1	ovrna=100,Funhe2Exx11m006564t1/I100,Funhe2Exx11m006564t10/68.80,Funhe2Exx11m006564t2/56.98,Funhe2Exx11m006564t4/47.83,Funhe2Exx11m006564t7/40.97
Funhe2EKm000011t1	ovrna=100,Funhe2Exx11m006564t1/I100,Funhe2Exx11m022615t4/57.61,Funhe2Exx11m019538t11/8.08,Funhe2Exx11m019538t4/6.06,Funhe2Exx11m019538t1/5.03,Funhe2Exx11m019538t2/5.03,Funhe2Exx11m019538t3/5.03,Funhe2Exx11m019538t5/4.05
Funhe2EKm000011t1	rnax=node,4a/10e/3g
      
=cut

sub get_scoretab
{
  my($inew,$rha,$rvec)= @_;  
  
  my($id,@v,%v,$dupv); $dupv={};
  if($inew == 2) { #rha = keyval hash; if(ref($rha) =~ /HASH/) ..
    %v= %$rha; $id= $v{ID} or return -1; 
    $dupv=$rvec if(ref($rvec) =~ /HASH/);
    # $rvec == %dupkeyref here
  } elsif($inew == 1) {
    # what? # chomp($rha); my @ha=split"\t",$rha;
  } else {  
    ($id,@v)= @$rvec;
    %v= vmap($rha,\@v); 
  }
  
  #Funhe2EKm000011t1	oid=Funhe2Exx11m006564t1	altc=althi	Split=2/2	cov=97.6	nexon=7	osrc=kf2rae5	locus=Scaffold43:402782-424715:-
  #Funhe2EKm000011t1	oid=Funhe2Exx11m006564t1,Funhe2E6bm006583t2	altc=althi	Split=1/2	cov=2.1	nexon=2	aalen=683,98%,complete	
  ## add dups, only for Splits ?
  if($v{Split}) { 
    for my $k (qw(nexon nintron locus location cov oid)) {
    if(my $dv=$dupv->{$k}) {
      ($dv)=split"\t",$dv; # tab sep vals
      my $v=$v{$k};
      if($v and $dv) {
        if($k=~/location|locus/) { $v.="/$dv"; } 
        elsif($k=~/nexon/) { $v= $v + int($dv); }
        elsif($k=~/nintron/) { $v= _max($v,$dv); } ## nintron=12/12 , nintron=2/2
        elsif($k=~/cov/) { $v= int($v + $dv); }
        elsif($k=~/oid/) { $v= $dv if(length($dv)>length($v)); }
        $v{$k}=$v;
      } elsif($dv) { $v{$k}= $dv; }
    }
    }
  }
  
  map{ $v{$_}=0  if($v{$_} eq "na") } (qw(alttr tere nintron inqual));  # is alttr obsolete? no?
  map{ $v{$_}="" if($v{$_} eq "na") } (qw(osrc ovpro));   # ovrna??
  $v{location}= delete $v{locus} if($v{locus}); #? change key names or not
  
  #?was# $id= $v{ID}; 
  push(@ids, $id); #? not here? do in putattr ?
  
  my $attr= getattr($id); # $gattr{$id};   
  $$attr{ID}= $id;   # RECODE for Table: Damn Microstupid Excel barfs on ID column name
  my $gid; unless($gid=$v{gene}) { ($gid=$id) =~ s/t\d+$//; } $$attr{gene}=$gid; 
  ## $$attr{gene}= $v{gene}; # can be null; fix?  # recode for table : gene_id:

  ## fixme: aliases for aalen, cxlen:
  ## aaSize=163	cdsSize=65%,492/750 ; 
  my $aalen= $v{aalen} || $v{aaSize};  # aalen,pctcode,quality : 145,35%,curated-complete
  my $cxlen= $v{cxlen};  # cxlen=3283/4888,67% cds/tr,%pctcode  ** MISSING in some

  if( not $cxlen and $v{cdsSize} =~ /,/ ) {
    $cxlen=$v{cdsSize}; 
    #?no? $cxlen=~s/^(\d+%),(\S+)/$2,$1/; 
  }
  if( not $cxlen and $aalen =~ /,/ ) {  
    my($aal,$aap,$aaq)= split",",$aalen;
    $aap =~ s/\%//; 
    my $cdsl= 3*$aal;  
    my $trl= ($aap > 0 and $aap < 100) ? int($cdsl * (100/$aap)) : $cdsl;
    # $cxlen="$cdsl/$trl,$aap%"; # wrong way now..
    $cxlen="$aap%,$cdsl/$trl";
  }

  my $pcds=0;
  my($clen,$trlen)= $cxlen=~ m/(\d+).(\d+)/;  
  if($cxlen =~ /(\d+)%/) { $pcds=$1; }
  elsif($aalen =~ m/,(\d+)%/) { $pcds=$1; }
  unless($pcds>0 or $trlen==0) { $pcds= int(0.5 + 100*$clen/$trlen); }
  
  $aalen =~ s/,.*//; 
  unless( $cxlen=~s/(\S+),(\d+%)$/$2,$1/ or $cxlen=~/%/) { $cxlen = "$pcds%,$cxlen"; }
  # $cxlen = "$pcds%,$cxlen" ; #was# $aalen = "$aalen,$pcds%";
  
  $$attr{aaSize} = $aalen;
  $$attr{cdsSize}= $cxlen;
  
  # FIXME : missinc trlen, cxlen: get from where?  aalen*3 + utrfudge ?
  if($trlen<1 and $aalen>0) { $trlen= 3 * $aalen; $trlen += int(0.1*$trlen); }
  $$attr{trlen}= $trlen; #$trlen||=1; 
  
  my $oid=$v{oid}; $$attr{oid}= ($v{osrc}) ? "$v{osrc}:$oid" : $oid;
  if($simpletab) { $$attr{oid}=~s/,.*//; }
  my $scr=$v{scoresum}||$v{score}; $scr=~s/,.*//; $$attr{score}= $scr; #buggy data
  $$attr{location}= $v{location} if($v{location} and ($formatTABLE or $doxml)); ## locus == location alias
  
  # FIXME: isoform == 0 for no alts, == 1 if has alts? change main-attr.tab input  **
  $$attr{isoform} = 0; # fix above:  and $v{alttr} ne "na"
  if($v{alttr}) { my($t)= $v{ID} =~ m/t(\d+)$/; $$attr{isoform} = $t||$v{alttr}; } #?? this can be wrong?
  
    # note: missing nintron has no introns; should have that flagged.
    ## inqual=90	nintron=40/44
    ## inqual=-50	nintron=0/2	inerr=-1; inqual=88	nintron=24/26	inerr=-1
  my $inqual="None";  # use none, not poor for no evidence introns
  my $nin= $v{nintron} || 0; ## nintron	29204
  my $inqual1= $v{inqual} || 0; # inqual	29204 ; use this if nintron ni/nx missing?
  my $inerr= $v{inerr} || 0; # inerr 3002; add as flag?
  
  ## fix: Thecc1EG005832t3        2/12,longerr:332434,
  if($simpletab) { $nin =~ s/,(\w*err)/-$1/g;  }
  $nin=~s/,$//;
  $$attr{intron}= $nin; # was inexon=in/ex
  
  my($ni,$ns,$pin); 
  if($nin =~ m,(\d+)/(\d+),) { 
    ($ni,$ns)=($1,$2); ## split"/",$nin;
    $pin= int(0.5+100*$ni/$ns) if($ns>0); ##? same as $inqual1
    } 
  elsif($inqual1) { $ns=1; $pin=abs($inqual1); } # inqual may be neg for inerr
  if($ns>0) {
    if($ns>9) { $inqual=($pin>74)?"Strong":($pin>66)?"Medium":($pin>0)?"Weak":"None"; }
    elsif($ns>3) { $inqual=($pin>99)?"Strong":($pin>84)?"Medium":($pin>0)?"Weak":"None"; }
    elsif($ns>1) { $inqual=($pin>99)?"Medium":($pin>0)?"Weak":"None"; }          
    $$attr{intron}= "$pin%,$nin";
  }
  $inqual.="-inerr" if($inerr<0);
  $$attr{qualIntron} = $inqual;
  
  $$attr{terepeat} = $v{tere} if $v{tere}; # transposon bases; keep or not?
  
  if( $v{aalen} =~ /,/ ) { # ** FIXME: do after all tables, qualHomology not read yet?
    my($aal,$aap,$aaq)= split",",$v{aalen};
    $aap =~ s/\%//;  
    
    unless($aaq) { # BUG: missing this aalen annot for 1700 genes .. better to fix this for input scores
      $aaq= "unavailable"; #?? drop this? or not? : 80 Protein:unavailable; 162 Protein:utrpoor_unavailable
    } 
      
    ## FIXME: have some StrongOrtho + StrongerParalog genes. MOVE to fixClass
    if(not qualProteinDropUTRflags and $aap =~ /^\d/ and not $aaq =~ /utr|noncode|poor|bad/) {  #  and $$attr{qualHomology} !~ /Homology:Ortholog(Strong|Medium)/
      if($aap < $UTRPOORisNONCODING) { $aaq= "noncode_$aaq"; }  #? drop this, utrpoor_ for all?
      elsif($aap < pUTRPOOR) { $aaq= "utrpoor_$aaq"; } #?? << utrpoor instead; FIXME: use -class.tab instead
    }      
    $aaq =~ s/curated-complete/curated_complete/; # is this still valid?
    $$attr{qualProtein}= $aaq;
    ## FIXME: qualProt .. leave off utrpoor/bad/.., messes like 
    ##   Protein:utrpoor_utrpoor_partial5,Protein:utrpoor_unavailable
    ##   Protein:utrpoor_utrpoor_complete
  }
   
  unless(MUST_NOT_EXPERT) {
  $$attr{'qualExpertchoice'} = 1 if($v{must} > 0);  # must= 77 or 69 (alt)
  }
  

  ## Split should be attached to location/locus **
  ## Cov field should be somewhere in annot output: attach to locus?
  ## Separate new mapqual field maybe best: include Split, coverage, percent ident if wanted
  ## .. values for qualMap statement.  Leave Split out of location= field but keep loc1/loc2 or loc1,loc2
  ## qualMap : cov, Split, .. pctCover cuts for Strong/Med/Weak/Poor/None
  ##  change qual % for qualMap: strong => 95, medium => 90, poor => 80 or split, none/weak => NOPATH
  my $qualmap=""; # statement from mapqual...
  my $cov= defined($v{cov}) ? $v{cov} : 100; # missing == perf map?
  $cov=~s/,.*$//; $cov=~s/%//; $cov =~ s/\.\d+//; # cov=91%,4796/5255
  
  if( $v{location} =~ /^(NOPATH|NOMAP)/ or $cov < 5) { $qualmap = "None"; $cov=0; } #  ??
  elsif( $cov < 80 ) { $qualmap = "Poor"; }
  elsif( $cov < 90 ) { $qualmap = "Okay"; } # Weak ??
  elsif( $cov < 95 ) { $qualmap = "Medium"; }
  else { $qualmap = "Strong"; }
  
  my $mapqual=$cov."%";
  
  ## $$attr{Split}=0; # chimera: change this to Split; retain old usage, rename Split
  # see below: replace attr{Split} with attr{mapCover}
  if(my $chi= $v{Split}) {
    $chi =~ s/,.*//; #chi=1 or 2, place chim1/2 here
    #Above# my $dspl= $dupv->{Split};
    #Above# my $dloc= $dupv->{locus} || $dupv->{location};
    $mapqual.=",Split".$chi; ## $$attr{Split}= $chi; 
  }
  elsif(my $chi= $v{chimera}) {
    $chi  =~ s/,.*//; #chi=1 or 2, place chim1/2 here
    $chi .= ",".$v{chim1} if($v{chim1});
    $chi .= ",".$v{chim2} if($v{chim2});
    my $eg= ($v{egover}); $eg=~s/:\d+//g; # keep this or drop?
    # $chi .= ",".$eg if ($eg); 
    $mapqual.=",Split".$chi; ## $$attr{Split}= $chi; 
  }
  if( $mapqual =~ /Split/ ) {  
    $qualmap .= "-Split";  
    ## my $spl=$mapqual; ## "Split" . $$attr{Split}; # .$spl;
    #Above.now# my $dloc= $dupv->{locus} || $dupv->{location}; chomp($dloc); $spl.=",$dloc" if($dloc);
    #NOT NOW# $$attr{location} .= ":$spl" if($$attr{location}); 
  }
  $$attr{mapCover}= $mapqual;  $$attr{mapCover}=~s=,=/=g if($simpletab);
  $$attr{qualMap} = $qualmap;
  
  ## FIXME for score=0 for many genes not via bestgene scoring..
  #DROP# my $scorevec= $v{scorevec}; 
  #DROP# if($scorevec and not $simpletab) { $$attr{scorevec}= join",", map{ s,/.*,,; $_ } split ",", $scorevec; }

  my $iste= ($v{class} and $v{class} =~ /transposon/) ? 1 : 0;  # NOT here ?? or {class} transposon?
  $iste=1 if($v{tere} > pTEWEAK); #?? or  or $$attr{terepeat}
  $$attr{qualTransposon} = 1 if($iste); 

} ## scoretab


sub get_classtab { # elsif($tb eq "c") 
  my($inew,$rha,$rvec)= @_;  

  my($id,@v,%v);
  if($inew == 2) { #rha = keyval hash; if(ref($rha) =~ /HASH/) ..
    %v= %$rha; $id= $v{ID};  
  } elsif($inew == 1) {
    # what? # chomp($rha); my @ha=split"\t",$rha;
  } else {  
    ($id,@v)= @$rvec;
    %v= vmap($rha,\@v); 
  }

  my $attr= getattr($id); 
  ## vmap does this:  $w =~ s/=/:/g; $w =~ s/;/,/g;
  
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
  
#     my $eg= $v{estgroup}; # not for wasp; switch to tilegroup?
#     $$attr{estgroup}= $eg if($eg);

}



=item nametab

# kf2rae5g.main.nametab
gene_id          name    namepct nameref oname   onamepct        onameref        ggname  ggnamepct       ggnameref       orclass
Funhe2EKm000003t1       na      0       0       Homeobox protein CDX-4  130     UniProt:CDX4_HUMAN      Homeobox protein CDX-1  69      FISH11G_G14131  orlog1
Funhe2EKm000004t1       platelet-derived growth factor receptor, beta polypeptide       1889    platyfish:ENSXMAP00000008599    Platelet-derived growth factor receptor beta    887     UniProt:PGFRB_HUMAN     Platelet-derived growth factor receptor, beta polypeptide       52      FISH11G_G998.s3 orlog1
Funhe2EKm000005t1       colony stimulating factor 1 receptor    1817    platyfish:ENSXMAP00000008640    Macrophage colony-stimulating factor 1 receptor 730     UniProt:CSF1R_HUMAN     Vascular endothelial growth factor receptor kdr 48      FISH11G_G553.s4 orlog1

  ** added orclass, capture it 

  - added ggname  ggnamepct       ggnameref == orthogroup name
  - updated cols: name    namepct nameref oname   onamepct        onameref
  - preserve oname/opct/oref if found as new field?

=cut

sub get_nametab {  # elsif($tb eq "n")
  my($inew,$rha,$rvec)= @_;  
  my($id,@v,%v);
  if($inew == 2) { #rha = keyval hash; if(ref($rha) =~ /HASH/) ..
    %v= %$rha; $id= $v{ID}; 
#   } elsif($inew == 1) {
#     # what? # chomp($rha); my @ha=split"\t",$rha;
#   } elsif($inew == 0) {  
#     ($id,@v)= @$rvec;
#     %v= vmap($rha,\@v); 
  } else {
    die "ERR: inew=$inew";
  }  

# FIXME: input name,homolog tables have messy vals for this
#   1. homolog score/namepct score now often bitscore, need %align from nameln= in *.bltab > homolog.tab > name.tab
#   2. name, oname, ggname now too confused.  
#       2a. Use ggname = gene name for most cases, if nameln align high.
#       2b. keep oname = human name where available.
#       2c. add 3rd name = best fish ho name for certain cases where diff from ggname, huname ?
#       Maybe add CDD domain names, elsewhere : new field? later... want all CDD/gene, not just best.

  my $attr= getattr($id); 
  
  # which first? dup keys, Name= in scoretab, name= in name.tab; should remove scoretab names when have nametab.
  #  and namepct/nameln
  ## FIXME HERE DECIDE which of name.kvtab to use
  ## .. name=fish-ortholog, not best; oname=human-ortholog, not best but keep; ggname= best unless missing..
  
  sub isname { return ($_[0] and $_[0] ne "na" and $_[0] ne $NAME_NONE)?1:0; }
  ## see sub hasval()
  
  my($na,$pna,$naref,$nabest)=("","",0);
  
  # need more consistent bestname selector
  #? maybe change ggnamepct >= pMEDIUM == 33 to lower criteria? got mix of nameref choices.. or better nameref
  if(1) {
  my($pnn,$pon,$pgn) = map{ my $p=$v{$_}; $p=~s/%//; $p=0 if($p eq "na"); $p; } qw(namepct onamepct ggnamepct);
  my($isnn,$ison,$isgn)= map{ isname($v{$_}) } qw(name oname ggname);
  
  if($isgn and $pgn>=pWEAKNAME) { # was pMEDIUM
    $na=$v{ggname}; 
    $naref= $v{nameref} || $v{ggnameref}; 
    $pna= _max($pgn,$pnn); 
    $nabest="ggname:".$v{ggnameref};
  } 

  unless(isname($na) and $pna>=pWEAKNAME) {
    if($ison and $pon >= pWEAKNAME and ($pon > $pnn or not $isnn)) { 
      $na= $v{oname}; $naref= $v{onameref}; $pna= $v{onamepct};  $nabest="oname:".$naref;
    } else { 
      $na= $v{name}; $naref= $v{nameref}; $pna= $v{namepct}; $nabest="name:".$naref;
    }
  }
  }
  #....
  # unless(isname($na) and $pna>1) { $na= $v{name}; $naref= $v{nameref}; $pna= $v{namepct}; }
  # unless(isname($na) and $pna>1) { $na= $v{oname}; $naref= $v{onameref}; $pna= $v{onamepct}; }

  $na =~ s/LOW QUALITY PROTEIN:\s*//; # dang NCBI junk should fix inputs..
  unless(isname($na) and $pna >= pWEAKNAME) { $na= $NAME_NONE; $nabest="none"; } ## defer NAME_NONE ??
  #OFF# $na= nameclean($na); # use qualities to decide for nonames: 
  
  # FIXME: input -name.tbl has bad namepct; replace with ortholog or uniprot pct
  ## ALIAS namepct == nameln=100%,690/661,683
  ## DANG; names.tab has namepct=bitscore or nalign, not %align
  
  #above# $pna= $v{namepct} || $v{nameln}; #? need both?
  if($pna and $na ne $NAME_NONE) { # not pna >= pWEAKNAME, need to see lowqual names
    $pna =~ s/,.*//; $pna =~ s/\.\d+//; 
    $pna=99 if($pna>100); # wrong .. need proper namepct value.
    $pna =~ s/[CI]//; $pna =~ s/$/\%/ unless($pna =~ m/\%/);
    ## Fix this nameref for UniProt:xxx_HUMAN to H ?
    $pna .= uc(substr($naref,0,1)) if($naref);
    # if($naref) { $pna .= ($naref=~/_HUMAN/) ? "H": uc(substr($naref,0,1)); } # presumes naref is species:id
    $na .=  " ($pna)"; #  (55%U) or (66%T) ?
    }
  $$attr{Name} = $na; 
  $$attr{namebest} = $nabest; # add which
  
  ## skip tename when namepct < 20%
  my $iste= 0;
  $iste= isTEname($na) if($pna > pTEWEAK); # FIXME: from input table ..
 
  my $ona= $v{oname};  # Human now; keep always; "na" is no name  
  ## what if oname picked best above? switch to name here? or not
  #OFF# $ona= nameclean($ona) if($ona); # use qualities to decide for nonames: 
  ## FIXME: hasval($ona) wrong here, can have valid onamepct,onameref but no name..  (others also?)
  if( $ona ) { # NOT hasval($ona)
    my $ndiff= isname($na) ? namediff($ona,$na) : isname($ona);
    unless($ndiff) { $ona="same"; } # same/ditto/?
    elsif(isname($ona)) {
      if(my $opct= $v{onamepct}) {
        $opct =~ s/\.\d+//;  $opct=99 if($opct>100); 
        my $oref= $v{onameref};   # include naref ? abbrev? should already have ID
        $opct .= substr($oref,0,1) if($oref);
        $ona  .=  " ($opct)";  
        $iste= isTEname($ona) if($opct > pTEWEAK and !$iste);
        }
      }
    elsif( hasval($v{onameref}) and $v{onamepct} >= pMEDIUM) {
      # $ona=$v{onameref}; #? stick in Dbxref instead...
    }
    $$attr{'oname'} = $ona;
    # $$attr{'Human_name'} = $ona; ##  rename this field??
  }
  
  my $gna=$v{'ggname'}; # ?? FIXME this is now main name? drop this fld?
  if( isname($gna) ) { #   # dont add if same as Name, oname?
    my $ndiff= namediff($gna,$na);
    if($ndiff) { $ndiff= namediff($gna,$ona); }
    if($simpletab) { $gna =~ s=[,]=.=g; } #?? fixme groupname 1..5 
    if($ndiff or $formatTABLE or $doxml) { $gna="same" unless($doxml or $ndiff); $$attr{groupname} = $gna;  } # homol has genegroup source
  }

  if($simpletab) { $$attr{Name}=~s/,/+/g; $$attr{oname} =~ s/,/+/g; }
  
  $$attr{qualTransposon} = 1 if($iste); 
}


sub get_refgenetab {  # elsif($tb eq "g")
  my($inew,$rha,$rvec)= @_;  
  my($id,@v,%v);
  if($inew == 2) { #rha = keyval hash; if(ref($rha) =~ /HASH/) ..
    %v= %$rha; $id= $v{ID}; 
  } elsif($inew == 1) {
    # what? # chomp($rha); my @ha=split"\t",$rha;
  } else {  
    ($id,@v)= @$rvec;
    %v= vmap($rha,\@v); 
  }
  
  my $attr= getattr($id); 
  my $ref1 = $v{ref1gene}; # maybe na REF1tag = APHIDBASE
  my $ref2 = $v{ref2gene}; # maybe na REF2tag = RefSeq
# gene_id ref1gene        ref2gene
# Thecc1EG000001t1        CGD0000036/78.83        Tc01_t000030/10.23,Tc01_t000020/10.15

  ## fixme for simpletab
  if($simpletab) { $ref1 =~ s/,/+/g; $ref2 =~ s/,/+/g; }
  $$attr{equiv1}= $ref1 if(hasval($ref1));   # print colname = REF1tag
  $$attr{equiv2}= $ref2 if(hasval($ref2));  # print colname = REF2tag
}        


=item expresstab

* revise this, add rnax=xpressclass,3groupscores
kf2rae5g.main.attr.tab1 annot : ok? dont care about xasm_equal, or est vs rna cover now.
  rx: or dex: or xde:?
express=100%,rx:nodiff,0.02a/0.1e/0.04g 
express=100%,rx:adultenv,1000a/2e/200g

../score5d/kfish2rae5d.main.express.tab
gene_id xcover  est_cover       rna_cover       xasm_equal      trspan
Funhe2EKm000003t1       679     679     679     100     679
Funhe2EKm000004t1       4888    3570    4888    100     4888
Funhe2EKm000005t1       4096    3075    4096    100     4096

=cut

sub get_expresstab {  # elsif($tb eq "x")
  my($inew,$rha,$rvec)= @_;  
  my($id,@v,%v);
  if($inew == 2) { #rha = keyval hash; if(ref($rha) =~ /HASH/) ..
    %v= %$rha; $id= $v{ID}; 
  } elsif($inew == 1) {
    # what? # chomp($rha); my @ha=split"\t",$rha;
  } else {  
    ($id,@v)= @$rvec;
    %v= vmap($rha,\@v); 
  }
  
  my $attr= getattr($id);
  my $trlen= $$attr{trlen}||1; # ** should compute after all tables read
  my($px,$ov)= (0) x 10;
  $ov= _max( $v{xcover}, _max( $v{est_cover}, $v{rna_cover}));
  $px= pct( $ov/ $trlen) if($ov);
  if( $ov= $v{ovest} ) { $ov=~s/,.*//; $px= _max($px,int($ov)); }
  if( $ov= $v{ovrna} ) { $ov=~s/,.*//; $px= _max($px,int($ov)); }
  
  my $pasm= $v{xasm_equal}; # est/rna assembly equiv; drop?
  my $rnax= $v{rnax}; # new diffx class,scores  rnax=node,0a/0.06e/0g
  if($rnax) {
    $rnax =~ s/^node\b/nodiff/; # fix poor diffx code
    my($rvx,@rv)= sort{$b<=>$a} $rnax =~ m,([\d\.]+)\w,g;
    $pasm= ($rvx >= 5)?pSTRONG :($rvx>=1)?pMEDIUM: ($rvx>=0.1)? pWEAK: $pasm;
    # check for high rpkm vals vs low/0 px .. replace pasm score w/ rpkm pseudo%
    # rpkm >=5 pasm=pSTRONG  >=1 pasm=pMEDIUM ; >= 0.5 pasm=pWEAK 
    if($simpletab) { $rnax=~s=,=/=g; }
  }
  my $xat="$px\%";
  $xat.=",rx:$rnax" if($rnax);
  #drop# $xat.=",eq:$pasm" if($pasm); ## DROP this?
  $$attr{express} = $xat; 
  
  my $q= pqual( _max($px,$pasm) ); 
  $$attr{qualExpress} = $q;
}  

  
sub get_homologtab {  # elsif($tb eq "h")
  my($inew,$rha,$rvec)= @_;  
  my($id,@v,%v);
  if($inew == 2) { #rha = keyval hash; if(ref($rha) =~ /HASH/) ..
    %v= %$rha; $id= $v{ID}; 
  } elsif($inew == 1) {
    # what? # chomp($rha); my @ha=split"\t",$rha;
  } else {  
    ($id,@v)= @$rvec;
    %v= vmap($rha,\@v); 
  }

  my $attr= getattr($id); 
  # map{ $v{$_}="" if($v{$_} eq "na") } (@HOKEY);  
  map{ delete $v{$_} if($v{$_} eq "na") } (@HOKEY);  
  
  ## FIXED: now ortholog=# if($v{homolog} and not $v{ortholog}) { my $ho=$v{homolog}; $ho=~s/^[\d\.]+,//; $v{ortholog}= $ho; }
  ## FIXED: drop this, now uniprot=# if(my $hu=$v{humalog}) { $hu=~s/^[\d\.]+,//; $$attr{uniprot} = $hu; }
  
  ## paralog now comes from ? orclass=upar,Id1,Id2,...
  if($v{orclass} =~/inpar|upar/ and not $v{paralog}) {
    my $pd=$v{orclass}; $pd=~s/^\w+,//; $v{paralog}=$pd;  
  }
  ## replace pbest with orclass ***
  # my $pho= $v{pbest};    ## orig: pbest=80%P =66%O ..
  # my $qho= ($pho =~ /O/)?"Ortholog":($pho =~ /P/)?"Paralog":"";
  
  ## ggname=Platelet-derived growth factor receptor, beta polypeptide	ggnamepct=52	ggnameref=FISH11G_G998.s3
  $v{genegroup}= $v{ggnameref}  if($v{ggnameref} and not $v{genegroup}); # fixup
  $v{genegroup} =~ s/^(\d+%),// if($v{genegroup});
  
  ## add orthoweak flag == sig blast hit but fails orthomcl inclusion as ortholog; ie fragment of others or ..
  ## .. need to adjust bits/selfbits or other part of orthodata to say it is NOT strong enough for validation
  ## fish ortholog/paralog tables: from orthomcl, lacking bitscores now (need for bestgenes but not annot table)
  ## as ortho/paralog=pctid%,geneid,genegroup,ngrpgene  
  # Funhe5EG000138t1 91%,242,tilapia:ENSONIP00000005011,weak 53%,141,Funhe5EG014545t1,weak   91%Oweak        na
  
  ## fixme: ortholog=0%,0,UniProt:na,Arp:na << drop if 0 bits, and/or "na"
  
  # $$attr{ortholog} = $v{ortholog}; # ($v{ortholog} =~ /0%,0,/) ? "" :
  if(hasval($v{ortholog})) { $$attr{ortholog} = $v{ortholog}; }
  elsif(hasval($v{uniprot})) { $$attr{ortholog} = $v{uniprot}; }
  for my $k (qw(paralog Dbxref genegroup)) { $$attr{$k}= $v{$k} if(hasval($v{$k})); } ;
  # $$attr{paralog}  = $v{paralog} if($v{paralog}); # ($v{paralog} =~ /0%,0,/ ) ? "" :
  #f    $$attr{uniprot}  = $v{uniprot};  # FIXME tag names
  # $$attr{Dbxref}   = $v{Dbxref};   # FIXME tag names
  # $$attr{genegroup} = $v{genegroup};   # FIXME tag names
  ## add orthomcl group here? as genegroup, or ggroup or gcluster or orthomcl ..
  
  # @HOKEY= qw( Dbxref ortholog paralog uniprot genegroup);
  # Dbxref,uniprot: output only id; or id first;  input now= nn%,bb,ID : do for all 4?

  ## reformat Hopct,Hoval2,HoID > HoID,Hopct
  for my $k (qw(ortholog uniprot)){ if(my $h=$$attr{$k}) { 
   $h =~ s/(\d+%),\d[^,]*,(.+)/$2#$1/;  $h =~ s/,.*//; $h =~ s/#/,/; $$attr{$k}=$h;  
   } } 

  ## FIXUP Dbxref should have all ref IDs..
  ## more fixups needed: 
  # old CDD:197700,TrEMBL:UniRef50_B4DXW7,TrEMBL:B4DXW7_HUMAN,
  # new CDD:197700,TrEMBL:UniRef50_B4DXW7,TrEMBL:B4DXW7_HUMAN,  << TrEMBL vs UniProt same id
  #  .. platyfish:ENSXMAP00000008658,UniProt:HMGX3_HUMAN,100%,1786,platyfish:ENSXMAP00000008658,64%,875,UniProt:HMGX3_HUMAN

  use constant DBXREF_UPDATE => 0; # probably want this, later...
if(DBXREF_UPDATE) {
  my $dx= (hasval($v{Dbxref}))? $v{Dbxref} : "";
  my @dx= map{ s/TrEMBL:/UniProt:/; $_;} split",",$dx; ## should do this in input tables
  my %dx=map{$_,1}@dx;
  for my $k (qw(nameref onameref ortholog uniprot)) { ## ggnameref ??
    ## ortholog=100%,1817,platyfish:ENSXMAP00000008640	uniprot=100%,730,UniProt:CSF1R_HUMAN
    if(hasval($v{$k})) { my $v=$v{$k}; 
      my($d)= $v=~m/(\w+:[^,\s;]+)/;  # $v=~s/,.*//; 
      if($d and not $dx{$d}) { push @dx,$d; $dx{$d}=1; } 
    }
  }
  $dx=join",",@dx; if(hasval($dx)) { $dx.=","; $$attr{Dbxref}= $v{Dbxref}= $dx; } # extra ',$' to match old vers
}

  if(not $formatTABLE and hasval($$attr{uniprot})) { # merge Dbxref for .gff
    $$attr{uniprot} =~ s=,=/=g;  
    if($$attr{Dbxref}) { $$attr{Dbxref}=~s/,+$//; $$attr{Dbxref} .= "," unless($$attr{Dbxref}=~m/,$/); }
    $$attr{Dbxref} .= delete $$attr{uniprot};
  }
  # if($formatTABLE) { map{ $$attr{$_}="na" unless( $$attr{$_} ) } @HOKEY; } # defer to end..
  # if($simpletab) { $$attr{paralog}=~s/,.*//; $$attr{Dbxref} =~ s/,/+/g; }
  
  my $qho=""; my $pho=0;
  $pho= _max( $v{ggnamepct}, _max($v{namepct}, $v{onamepct})); ## fixup1

  if(my $orcl= $v{orclass}) { 
    my($orcn,$parids)= split",",$orcl,2;  #$orcn=~s/\d+$//;
    ## notor = Homology:Unique/None which?  Homology:Uniparalog/Unique ?  Homology:ParalogNone/None
    $qho=($orcn=~/^inpar/)?"Inparalog":($orcn=~/^orlog|orpar/)?"Ortholog":
         ($orcn=~/^upar/)?"Uniparalog":($orcn=~/^notor/)?"Unique":"Unique";
         
    ## fix2: use ggnamepct if larger here 
    my $pho1=0;    
    if($$attr{ortholog}) { ($pho1)= $$attr{ortholog}=~m/(\d+)%/; }    # skip: $qho =~ /Ortholog|Inparalog/ and   
    unless($pho1) { $pho1 = $v{nameln} || pWEAK; } # $v{namepct} << Bad data? FIXME here, get new homolog.tab aln% scores
    $pho= _max($pho, $pho1);
  } elsif(my $pho1= $v{pbest}) {
    $qho= ($pho1 =~ /O/)?"Ortholog":($pho1 =~ /P/)?"Paralog":"";
    $pho1=pWEAK if($pho1=~/weak/i and $pho1>pWEAK); # ortho/para flag from omcl-miss
    $pho1 =~ s/\D+//;
    $pho= _max($pho, $pho1);
  } elsif(hasval($$attr{ortholog})) {
    my($pho1)= $$attr{ortholog}=~m/(\d+)%/; 
    $pho= _max($pho, $pho1);
  }
  
  # this may be bad set pWEAK for everything? "na" bug, below now
  if( (not $qho or $qho =~ /Unique/) and ($pho>=pWEAK 
    or hasval($$attr{ortholog}) or hasval($$attr{Dbxref}) or hasval($$attr{genegroup})) ) { 
    $qho="Ortholog";
    $pho= pWEAK unless($pho>pWEAK); # some homology ..
  }
  
  ## FIXUP $pho when genegroup exists: None not valid..
  $pho= pMEDIUM if($pho < pMEDIUM and hasval($v{genegroup})); # homol validated by omcl, cant be weak
  $qho="" if($pho < pWEAK); ## drop xxxxNone classes, only None
  $qho .= pqual($pho); #?? yes or no: unless($qho=~/^Unique/);
  $$attr{qualHomology} = $qho;

  # defer to end..
  if($formatTABLE) { 
    #xx for my $k (qw(ortholog paralog)) { $$attr{$k}=0 unless( hasval($$attr{$k}) ); }
    for my $k (@HOKEY) { $$attr{$k}="na" unless( hasval($$attr{$k}) ); } 
    } 
  if($simpletab) { $$attr{paralog}=~s/,.*//; $$attr{Dbxref} =~ s/,/+/g; }
}
  

#------- ugp.xml -------------------------------

my $xin= 0;
sub xtab { my $t= ($xin>0) ? ("  ") x $xin : ""; print $t; }
sub unravel { my $t=shift; my @q=split",",shift; for my $q (@q) { ptag( $t,$q,0); } }

sub ptag {
  my($t,$v,$ln,$tpar)= @_;
  return unless(defined $v);
  xtab() unless($ln);
  ## tpar="id=xxx" or "url=yyy" or "type=zzz";
  my $tb=($tpar)?"$t $tpar":$t;
  print ($doxml? "<$tb>":"$tb: ");
  print ($doxml? "$v</$t>" : "$v;"); ## escape xml ***
  print (($ln) ? " " : "\n");
}

sub btag {
  my($t, $ln)= @_;
  xtab(); $xin++; print ($doxml?"<$t>":"$t:"); print (($ln) ? " " : "\n");
}

sub etag {
  my($t, $ln)= @_;
  $xin--; xtab() unless($ln or !$doxml);
  print ($doxml?"</$t>\n":"");  print "\n" if($ln and !$doxml); #??
}

=item ugp.xml

  inputs add, as ID key=value tables
    - proteins
    - kfish_omcl-all ortho fish scores
    - 
  FIXME:  
    ** MISSING FINAL etag "GeneSummaries"; </GeneSummaries> ** force it..
    remove _ uscores from tags for lucene default parser
      Family_name, Family_dbxref > FamilyName, FamilyDbxref
      
    unravel id1,id2,.. as sep sub fields;
      for Dbxref 
      - Dbxref not as Protein_domains but ExternalLinks
    Express/FPKM a,e,g change, remove name == ExpressGroup bad for search

  <EXPRESSION>   New ...
    <Introns>100%,2/2</Introns>  : break out or add inline val <Introns val=2>100%,2/2</Introns>
    <RNACoverage>100%</RNACoverage>
    <ExpressGroup>adultenv</ExpressGroup>
    <ExpressLevel group=adultenv>2</ExpressLevel>
    <ExpressLevel group=embryo>0.2</ExpressLevel>
    <ExpressLevel group=grandis>0.09</ExpressLevel>
    <RNAassembly>kf2x11gmap:Funhe2Exx11m009257t1</RNAassembly>
  </EXPRESSION>

      
  OPTION: add print_ugpxml() here, alt output format. need %attr not @v input.

=cut

sub put_ugpxml {
  my($ogid,$firstlast,$kvhash,$dupkvhash)= @_; # NOT:, $attr
  my $attr= $gattr{$ogid};   ##NOT my $attr= getattr($id); # makes empty attr; not here; 
  unless(ref $attr and $$attr{ID}) {
    if($firstlast == LASTID and not $didfoot) { etag "GeneSummaries"; $didfoot++; }
    return 0;
  }
  
  fixClass( $attr);     
  
  $xin= 0;
  if(1 > $didhead++) {  # $firstlast == FIRSTID 
    print "<?xml version=\"1.0\"?>\n";
    btag "GeneSummaries"; 
    # print_ugpdummy
    btag "GeneSummary id=\"header\""; 
    ptag "source","$GENE_SOURCE";
    etag "GeneSummary";  print "\n";
  } 
  
  sub isname { return ($_[0] and $_[0] ne "na" and $_[0] ne $NAME_NONE)?1:0; }
  sub namenopct { my $na=shift; $na=~s/\s*\(([\w%]+)\)$//; my $pn=$1; return ($na eq "same")?undef:$na; }  # (27%U)
  
  $xin= 0;
  btag "GeneSummary id=\"euGenes:$ogid\"";  # fixme: id= tag?
  ptag "Title", "$GENE_TITLE $ogid";
  ptag "Source","$GENE_SOURCE";
  ptag "Type","$GENE_TYPE";
  
  btag "BASIC_INFORMATION"; 
    ptag "GeneID", $ogid;
    ## name here instead/also ?
    ptag "Name",  namenopct($$attr{Name});  
    ptag "isoform",  $$attr{isoform}||0; # FIXME need isoform 0 value ; need field in all recs

    # ptag "Quality",  $$attr{quality}; # qual has subfields?
    # 	Class:Strong,Express:Strong,Homology:OrthologStrong,Intron:Strong,Map:Medium,Protein:complete
    #  split Ho val? Ho:Ortholog-Strong or Ho:Strong-Ortholog , so sep searchable
    my $qv=$$attr{quality};
    if($qv) { my($ho)= $qv=~m/Homology:([^,;:]+)/; 
      (my $hn=$ho)=~s/^(Inparalog|Ortholog|Uniparalog)(.+)/$2-$1/;
      $qv =~ s/Homology:$ho/Homology:$hn/;
      }
    btag "Quality"; unravel 'qual', $qv; etag "Quality";
 
    ptag "Species", $CLADE; #?? or not
  etag "BASIC_INFORMATION"; 
  
  btag "LOCATION"; 
    # ptag "Chromosome",  $date;
    my $loc=$$attr{location};
    # $loc=~s/:[+-\.]//g; # end strand causes problems for gbrowse url, leave or not?
    ptag "Genome_map",  $loc; # Scaffold0:25799-58393:+
    ptag "MapCoverage",  $$attr{mapCover}; # 91% ; Split ;
  etag "LOCATION"; 
  
  btag "GENE_PRODUCT"; 
  
    # FIXME here? dig in $kvhash for other flds: name,namepct,nameref; oname,onamepct,onameref
    my $didna=0; my $didana=0;
    if(1) {
      ## NOT good as attr{Name} chooser picks best from among all .. if this, need flag best
      if(hasval($$kvhash{nameref})) {
        btag "Name";
        ptag "naval", $$kvhash{name};  
        ptag "db_xref", $$kvhash{nameref};  ptag "align",$$kvhash{namepct};
        etag "Name"; $didna++;
      }
      if(hasval($$kvhash{onameref})) {
        btag "OtherName";
        ptag "naval", $$kvhash{oname};  
        ptag "db_xref", $$kvhash{onameref};  ptag "align",$$kvhash{onamepct};
        etag "OtherName"; $didana++;
      }
    }
    ptag "Name",  namenopct($$attr{Name}) unless($didna);  # Homeobox protein CDX-1 (91%P)
    ptag "OtherName",  namenopct($$attr{oname}) unless($didana);    
     
    if(1) {
    if(hasval($$kvhash{ggnameref})) {  
        btag "Family";
        ptag "naval", $$kvhash{ggname};  
        ptag "db_xref", $$kvhash{ggnameref};  ptag "align",$$kvhash{ggnamepct};
        etag "Family";  
      } 
    } else {
    if(isname($$attr{groupname})) {  
        btag "FamilyName";
        ptag "naval", $$attr{groupname};  
        ptag "db_xref",$$attr{genegroup};  #not there:ptag "align",$$attr{ggnamepct};
        etag "FamilyName";  
      }
    }
    ptag "NameSource", $$attr{namebest}; ## name,oname,ggname or none
    
    ptag "aaSize",  $$attr{aaSize};  # 214 ? here or basic?
    ptag "cdsSize",  $$attr{cdsSize}; # 94%,663/679
    
    # Transcript IDs/oids = mRNA assemblies here? or Expression ?

    ## ?? add aaseq here ...    
    if(my $aaseq=$$kvhash{protein}) {
      $aaseq =~ s/(.{1,60})/$1\n/g; chomp($aaseq); $aaseq="\n".$aaseq;
      btag "Sequence"; ptag "protein", $aaseq, 0, "id='$ogid'"; etag "Sequence"; 
    }
  
  etag "GENE_PRODUCT"; 
  
  btag "EXPRESSION"; # or FUNCTION/Expression
    ptag "Introns",  $$attr{intron}||0;  # 50%,2/4 : FIXME?? these are splice counts, /2 for nintrons
        # ^ fraction n/2 or int(n/2) or as is?
        
    my $ex=  $$attr{express};   ## 100%,rx:nodiff,0.02a/0.1e/0.04g << break out all parts
    if(hasval($ex)) {    
    my($xcov,$xgr,@fgr)= split m=[,/]=,$ex;
    $xcov||=0; $xgr=~s/rx://; ## $xgr||="na";
    ptag "RNAcover", $xcov; # 100%,rx:nodiff,0.02a/0.1e/0.04g
    ptag "ExprGroup", $xgr; # 100%,rx:nodiff,0.02a/0.1e/0.04g
    
    ## redo..
    # map { s/a$/ adultenv/; s/e$/ embryo/; s/g$/ grandis/; ptag "ExpressLevel", $_; } @fgr;
    btag "ExprLevel";
    foreach (@fgr) { my $tp="rgroup";
      if(s/a$//) { $tp="rAdultenv"; }
      elsif(s/e$//) { $tp="rEmbryo"; }
      elsif(s/g$//) { $tp="rGrandis"; }
      ptag "exprval", $_, 0, "group='$tp'";
      }
    etag "ExprLevel";
    }
    
    my $oid=  $$attr{oid}; my $xid= join",",grep(/Funhe2Exx11/,split(",",$oid));
    ptag "RNAassembly", $xid if($xid); # kf2x11gmap:Funhe2Exx11m027882t1,Funhe2Emap3m022605t1
    ## add Funhe2EKm000100t1 ovest=52,2145/4111 ?
    ##?? add all x11 ids from ovrna ? or from evg mixx11/rnaequiv tab?
    ## not all ovrna are good mappings .. should use also/instead evg equiv tables
    ## Funhe2EKm000100t1 ovrna=100,Funhe2Exx11m001640t1/I100,Funhe2Exx11m001640t2/C98.99,Funhe2Exx11m001640t5/88.95,
    # my $ovrna= $$kvhash{ovrna}; my $xov=join",",grep/Funhe2Exx11/, split",",$ovrna;
    # ptag "RNAoverlap", $xov if($xov);  # mostly alts, some others
    
    #?? ptag "ESTdb_xref",  $$attr{groupname};
  etag "EXPRESSION"; 
  
  
  btag "FUNCTION";  # enclose Expression, Gene_Ontology ?
    # btag "Protein_domains"; # CDD here ??
    my $dbx=$$attr{Dbxref};
    if(hasval($dbx)) {
      $dbx=~s/TrEMBL/UniProt/g; # dang NCBI url foulup
      btag "ExternalLinks";  unravel "db_xref", $dbx; etag "ExternalLinks"; 
    }
   
#     btag "Gene_ontology";  # maybe or not?
#       btag "Molecular_function"; 
#         ptag "goterm",  $$attr{groupname};
#       etag "Molecular_function"; 
#       btag "Biological_process"; 
#         ptag "goterm",  $$attr{groupname};
#       etag "Biological_process"; 
#       btag "Cellular_component"; 
#         ptag "goterm",  $$attr{groupname};
#       etag "Cellular_component"; 
#     etag "Gene_ontology";  

  etag "FUNCTION"; 
  
  
  btag "SIMILAR_GENES";  
     ## FIXME: add all orlogs from fish11g_kfish_orpar.tab
      # all orthos: use ocla/oref/oval from $kvhash,$dupkvhash
    if(hasval($$kvhash{ocla})) {
      my @ocla= ($$kvhash{ocla}, split"\t", $$dupkvhash{ocla});
      my @oref= ($$kvhash{oref}, split"\t", $$dupkvhash{oref});
      my @oval= ($$kvhash{oval}, split"\t", $$dupkvhash{oval});
      for my $i (0..$#ocla) {
        my $ocl=$ocla[$i]; # inpar, orlog, opar, upar, 
        (my $otype=$ocl)=~s/\d+$//;
        $ocl=~s/(\d+)$/ $1/; # drop number or not?
        $ocl=~s/^orlog/ortholog/;
        $ocl=~s/^opar/orparalog/;
        $ocl=~s/^inpar/inparalog/;
        $ocl=~s/^upar/uniparalog/;
        btag "Ortholog type='$otype'", 1;  
        ptag "class",$ocl,1;
        ptag "acc",$oref[$i],1;  
        ptag "align",$oval[$i],1;
        etag "Ortholog", 1;
      }
      
    } else { ## dont need both ..
      # ortholog: platyfish:ENSXMAP00000008562,91%
      # paralog: Funhe2EKm000100t1,91%
    my $olog= $$attr{ortholog};
    if(hasval($olog)) { # many??
      my($spp,$sid,$pal)= split/[:,]/,$olog;
      btag "Ortholog", 1; # was Similarity
      # ptag "spp",$spp,1; # Species
      ptag "acc","$spp:$sid",1; # Symbol/db_xref
      # ptag "def",$gdef,1 if($gdef); # Description
      ptag "align",$pal,1;
      etag "Ortholog", 1;
    }

    my $val= $$attr{paralog};
    for my $plog (split",",$val) {
      my($pid,$pal)= split "/",$plog;
      if(hasval($pid)) {
        btag "Paralog", 1;
        ptag "acc",$pid,1; # Symbol/db_xref
        ptag "align",$pal,1;
        etag "Paralog", 1;
        }
      }
    }
  
  etag "SIMILAR_GENES"; 
  
  btag "ADDITIONAL_INFORMATION";
    if( my $val= $$attr{equiv1} ) { # Funhe2EKm000100t1	ref1gene=Funhe5EG009378t1/91,
      ## $val=~s,/.*,,; ## not this; drops 2nd+ ids; 
      my $vid=join",", map{ s,/.*,,; $_ } split",",$val;
      ptag "PriorModels", $vid; }  
    if( my $val= $$attr{oid} ) { ptag "ObjectIds",  $val; }  
    ptag "date",  $date;
  etag "ADDITIONAL_INFORMATION"; 
  etag "GeneSummary";  print "\n";
  
  if($firstlast == LASTID) { # maybe only 1.. **FIXME missing this, sometimes .. data input effect?
    unless($didfoot) { etag "GeneSummaries"; $didfoot++ ; }
  }
}


#------- end ugp.xml -------------------------------

sub putattr {
  my($gid)= @_; # NOT:, $attr
  my $attr= $gattr{$gid};   ##NOT my $attr= getattr($id); # makes empty attr; not here; 
  return 0 unless(ref $attr);

  ## 0	0	0	Class:Strong,Express:Strong,Homology:None,Intron:None,Protein:None	0	0	Uncharacterized protein	na	na	0	0	na	na	0	0%,rx:nodiff,10a/20e/2g,eq:66	0	0	00
  ## Funhe2EKm023583t1 << single missing scores ID, in rnax.tab only
  return 0 unless($$attr{ID}); # report err above.
  
  fixClass( $attr);     
  my @v; my %didk;
  my @ak= sort keys %$attr;
  my @qk= grep /^qual/, @ak;  map{ $didk{$_}++; } @qk;
  #now in attr# my $qual= join ",", map{ $didk{$_}++; (my $k=$_ ) =~ s/^qual//; "$k:$$attr{$_}" } @qk;
     
  foreach my $k (@ATKEY) {
    my $v=$$attr{$k}; $didk{$k}++;
    my $kl= $k;
    # if($k eq "quality") { $v=$qual; } # now in attr
    
    # $v= $MISSING if($formatTABLE and not defined $v); # fixme "defined $v" should be printable $v
    $v= $MISSING if($formatTABLE and $v !~ /\S/);  
    
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
  return 1;
}


sub fixClass
{
  my ($attr) = @_;
  return unless(ref($attr) and not $$attr{'fixed'});
  
  ## fixme missing Homology: Class:Poor-utrpoor,Express:Weak,Intron:none,Protein:poor_complete 
  ## ensure all qual parts exist? : none default
  ## .. but special cases: ExpertChoice, FusionMaybe
  # Class:Strong,Express:Strong,Homology:OrthologStrong,Intron:Strong,Protein:complete
  # my @QUAL_CLASS=qw( Express Homology Intron Protein);

  map{ $$attr{$_}="None" unless($$attr{$_}) } @QUAL_CLASS; # excludes gene qualClass

      ## FIXME: have some StrongOrtho + StrongerParalog genes.
  if(qualProteinDropUTRflags) {
    if( $$attr{qualProtein} =~ /bad|poor|noncod/i) {
      $$attr{qualProtein} =~ s/[_-]?(utrpoor|utrbad|bad|poor|noncode)[_-]?//gi; #? maybe leave noncode flag? or drop utrjazz here
    }
  } else {
    if( $$attr{qualProtein} =~ /poor|noncod/i and
    ($$attr{qualHomology} =~ /Ortholog(Strong|Medium)/ or $$attr{ortholog} >= pSTRONG) 
    ) {
      $$attr{qualProtein} =~ s/(utrpoor|poor|noncode)[_-]?//gi; #? maybe leave noncode flag? or drop utrjazz here
    }
  }
  
  geneClass( $attr); # do before fix name
  
  my @ak= sort keys %$attr;
  my @qk= grep /^qual/, @ak;
  my $qual= join ",", map{  (my $k=$_ ) =~ s/^qual//; "$k:$$attr{$_}" } @qk; # $didk{$_}++;
  # my $qual = $$attr{'qualClass'};
  $$attr{'quality'}= $qual;
  
  #? here set/reset Name = hypothetical protein depending on quality of evidence
  #    hypoth = weak evid;  uncharacterized = strong evd;  ? expressed uncharact = strong express
  my $na= $$attr{'Name'};
  $na= $NAME_NONE  if($na eq "na" or $na !~ /\w/); #"hypothetical protein"
  if($na =~ /$NAME_NONE/i) {
    $na =~ s/$NAME_NONE/$NAME_NOFUN/i 
      if($qual =~ /Class:(Strong|Medium)/); ## (Express|Homology):\w*(Strong|Medium)/i);
    $$attr{'Name'}= $na;  
  }
  
  $$attr{'fixed'}=1;
}


sub geneClass
{
  my ($attr) = @_;
  my ($perf,$st,$md,$wc,$poor,$upoor,$partial,$must,$te)=(0)x10;
  my @qk= grep /^qual/, sort keys %$attr;

  # FIXME: modify Intron qual with attr{intron} == nintron, e.g. intron==2 + inqual=100 not strong but medium
  foreach my $q (@qk) {
    ##next if($q =~ /qualClass/); # /qualMap/ ??
    next if($q =~ /qualMap/); # genome quality, not relevant to gene quality
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
  
  my $cl= ($perf>10000)? "VeryStrong" #? "Perfect" never here...
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


# YEATS family protein vs YEATS family protein. sym:GAS41,TAF14B 
sub namediff {
  my($na,$nb)= @_;
  # check for substantial diff 
  # fixme: Metalloproteinase (58%M) <> Metalloproteinase
  map { s/\s*\([^\)]+\)\W*$//; s/[,:\(\)].*//; s/\W+/ /g; $_= lc($_); } ($na,$nb);
  return 0 unless(hasval($na) and hasval($nb)); # missing not diff, or ret -1 ?
  return ($na eq $nb or $na =~ /$nb/ or $nb =~ /$na/) ? 0 : 1;
}


## add TE protein name classing; also quality flag Transposon
## nameclean STUB now, use evg/prot/protein_names.pm if need be, presume input names are already clean?

sub nameclean {  
  my($na,$hasid)= @_;  my $id="";
  if($hasid and $na=~/^\w\S+\s/){ ($id,$na)= split" ",$na,2;  } 
  return (wantarray)?($na,$id):$na; 
}

sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }

sub vmap { 
  my($h,$v)=@_; my %r=(); 
  for my $i (0..$#$h){ 
    # some cleanup should be elsewhere: no [=;] in fields; drop null='.'
    my $w=$$v[$i] || ""; $w="" if($w eq ".");  $w =~ s/=/:/g; $w =~ s/;/,/g;
    $r{$$h[$i]}= $w; 
    } 
  return %r; 
}

sub vclean { # for hashed vals
  my($h)= @_; return unless(ref($h));
  for my $k (keys %$h) { 
    # some cleanup should be elsewhere: no [=;] in fields; drop null='.'
    # NOTE: dupkey hash has \t separator : change?
    my $w=$h->{$k} || ""; $w="" if($w eq ".");  $w =~ s/=/:/g; $w =~ s/;/,/g;
    $h->{$k}= $w; 
  } 
  return $h; # or new copy of %$h; ??
}

sub pct { my $p=shift; return int(0.5 + 100*$p); }
sub pqual { my $p=shift; $p=~s/\D+//; ## $p=pWEAK if($p=~/weak/i and $p>pWEAK); # ortho/para flag from omcl-miss
  return(($p >= pSTRONG)? "Strong" : ($p >= pMEDIUM) ? "Medium" : ($p >= pWEAK) ? "Weak" : "None"); 
}

sub isTEname {
  my $na= shift;
  return 0 unless($USE_TENAME);
  foreach my $te (@TEnames) { return 1 if($na =~ /$te/i); }
  return 0;
}


sub getTEnames
{
my @TEnames= map{ s/^\d+\s*//; $_; } grep /^\d/, split "\n", <<'EOTE';
999 Retrotransposon
999 transposon
999 Gag.pro
999 Gag.Pol
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


__END__

=item kfish2rae5g.main.attr.tbl, 2jan2014

24282 Class:Strong          : 300 more
7360 Class:StrongPartial    : same
1910 Class:Transposon       : 11 more
 940 Class:Medium           : 100 less
 379 Class:Weak             : 150 less
  39 Class:MediumPartial
  13 Class:WeakPartial
   2 Class:None             : 123 less, what are 2 left?
      Funhe2EKm006177t1 None, 59aa, Scaffold231:292973-296284:-,kf4bAUGpie6s231g13t1
        paralogs=Funhe2EKm011865t1,Funhe2EKm016624t1,Funhe2EKm007752t1,Funhe2EKm011692t1,Funhe2EKm015634t1,Funhe2EKm008230t1,Funhe2EKm011752t1,Funhe2EKm017991t1,Funhe2EKm006173t1,Funhe2EKm006218t1
      Funhe2EKm031266t1 None, 641aa, 0%,rx:0.02a/0e/0.02g; Scaffold9967:1214288-1253693:+  Funhe5EG009185t1/6,Funhe5EG010026t1/8,  kf4bAUGpia9cs9967g38t1
        ^^ reclass Weak for some rx: expression
        
#KC     Dbxref  31993
#KC     ID      34998
#KC     Name    32345
#KC     Split   2251
#KC     aaSize  19      << WRONG, or make all aalen > aaSize ??
#KC     aalen   34908   << missing some ?
#KC     alttr   20036   ? isoform
#KC     cdsSize 19
#KC     cov     27948   ? need more?
#KC     cxlen   28841   ? missing?
#KC     ggname  32326
#KC     ggnamepct       32326
#KC     ggnameref       32326
#KC     homolog 30559
#KC     inerr   3004
#KC     inqual  29233
#KC     locus   34925
#KC     must    19376
#KC     name    32326
#KC     nameln  29886
#KC     namepct 32326
#KC     nameref 32326
#KC     nexon   25682
#KC     nintron 29233
#KC     oid     34925
#KC     oname   32326
#KC     onamepct        32326
#KC     onameref        32326
#KC     orclass 32326
#KC     ortholog        32326
#KC     osrc    34925
#KC     ovest   34920
#KC     ovpro   23833
#KC     ovrna   29424
#KC     ref1gene        26310
#KC     rnax    34903
#KC     score   437     # drop
#KC     scoresum        28839
#KC     sense   936
#KC     srcf    22
#KC     tere    6161
#KC     uniprot 32326
#KC     upd2    567
#KC     upd4    175
#KC     upd5g   33
#KC     upr5    3

=item kf2rae5g.main.attr.tbl2

pt=kf2rae5g.main
pt=kfish2rae5g.main

sort -k1,1 -k2,2  kf2rae5g.main.scores.tab kf2rae5g.main.{names,homolog}.kvtab kf2rae5g.main.alttr \
  kf2rae5g.main.ov* kf2rae5g.main.rnaxtab | env sorted=1 table=1 ../evgpuban_kfish2.pl \
  > kf2rae5g.main.attr.tbl2 


transcriptID	geneID	isoform	quality	aaSize	cdsSize	Name	oname	groupname	ortholog	paralog	genegroup	Dbxref	intron	express	mapCover	location	oid	score
Funhe2EKm000003t1	Funhe2EKm000003	0	Class:Strong,Express:Strong,Homology:OrthologStrong,Intron:Weak,Map:Strong,Protein:complete	214	94%,663/679	na (74%)	Homeobox protein CDX-4 (99U)	Homeobox protein CDX-1	platyfish:ENSXMAP00000008562,91%	na	FISH11G_G14131	CDD:200956,TrEMBL:UniRef50_Q90423,TrEMBL:HXB1B_DANRE,	50%,2/4	100%,rx:nodiff,0.02a/0.1e/0.04g	100%	Scaffold0:5067-10366:-	kf2x11gmap:Funhe2Exx11m027882t1,Funhe2Emap3m022605t1	5464
Funhe2EKm000004t1	Funhe2EKm000004	1	Class:Strong,Express:Strong,Homology:OrthologStrong,Intron:Strong,Map:Medium,Protein:complete	1103	63%,3283/4888	platelet-derived growth factor receptor, beta polypeptide (99%p)	samPlatelet-derived growth factor receptor, beta polypeptide	platyfish:ENSXMAP00000008599,100%	na	FISH11G_G998.s3	CDD:173624,TrEMBL:UniRef50_P35968,TrEMBL:VGFR2_HUMAN,	91%,40/44	100%,rx:nodiff,3a/7e/2g	91%	Scaffold0:25799-58393:+	kf2x11gspl:Funhe2Exx11m002607t1,Fungr1EG3m001115t1	30678
Funhe2EKm000005t1	Funhe2EKm000005	1	Class:Strong,Express:Strong,Homology:OrthologStrong,Intron:Strong,Map:Medium,Protein:complete	1010	67%,2952/4096	colony stimulating factor 1 receptor (99%p)	Macrophage colony-stimulating factor 1 receptor (99U)	Vascular endothelial growth factor receptor kdr	platyfish:ENSXMAP00000008640,100%	na	FISH11G_G553.s4	CDD:219530,TrEMBL:UniRef50_P35916,TrEMBL:VGFR3_HUMAN,	87%,40/46	100%,rx:nodiff,3a/5e/0.9g	90%Scaffold0:62577-77948:+	kf2x11gspl:Funhe2Exx11m003308t5,Funhe2E6bm003260t4	24193
Funhe2EKm000006t1	Funhe2EKm000006	1	Class:Strong,Express:Strong,Homology:OrthologStrong,Intron:Strong,Map:Strong,Protein:complete	1257	85%,3773/4416	HMG box domain containing 3 (99%p)	HMG domain-containing protein 3 (99U)	HMG box domain containing	platyfish:ENSXMAP00000008658,100%	na	FISH11G_G10666	CDD:197700,TrEMBL:UniRef50_B4DXW7,TrEMBL:B4DXW7_HUMAN,	100%,40/40	100%,rx:nodiff,1a/2e/0.2g	100%	Scaffold0:85415-96983:-	kf2x11gmap:Funhe2Exx11m001896t4,Funhe2E6bm001856t1	34214
Funhe2EKm000007t1	Funhe2EKm000007	0	Class:Strong,Express:Strong,Homology:OrthologStrong,Intron:Strong,Map:Strong,Protein:complete	325	34%,941/2830	solute carrier family 35, member A4 (99%p)	UDP-sugar transporter protein SLC35A4, putative (99U)	Solute carrier family 35, member A4	platyfish:ENSXMAP00000019533,100%	na	FISH11G_G8936	CDD:217924,TrEMBL:UniRef50_A4IHW3,TrEMBL:S35A4_XENTR,	100%,4/4	100%,rx:nodiff,8a/7e/2g	96%	Scaffold0:97595-107565:+	kf2x11gspl:Funhe2Exx11m019598t1,Fungr1EG3m012092t1	7046
Funhe2EKm000008t1	Funhe2EKm000008	1	Class:Strong,Express:Strong,Homology:OrthologStrong,Intron:Strong,Map:Strong,Protein:complete	543	83%,1632/1942	kinesin-like protein KIF17-like (99%m)	Kinesin-like protein KIF16B (75U)	Kinesin protein KIF17	mayzebr:XP_004571788.1,100%	na	FISH11G_G19221	UniProt:F8WA61_HUMAN	100%,28/28	81%,rx:nodiff,1a/1e/0.2g	100%	Scaffold0:107343-127504:-	kf4best:kf4bAUGarip1p1s0g6t1	8547
Funhe2EKm000009t1	Funhe2EKm000009	1	Class:StrongPartial,Express:Strong,Homology:UniparalogStrong,Intron:Strong,Map:Strong,Protein:partial	228	99%,686/686	zgc:163143 (68%z)	THAP domain-containing protein 4 (60U)	THAP domain-containing protein	zfish:ENSDARP00000098613,34%	Funhe2EKm024213t1,Funhe2EKm004412t1,Funhe2EKm005328t1	FISH11G_G14274	CDD:214951,TrEMBL:UniRef50_G5DY87,G5DY87_9PIPI,	100%,4/4	100%,rx:nodiff,1a/2e/0.3g	100%	Scaffold0:121590-124252:-	kf2x11gmap:Funhe2Exx11m026463t1,Funhe2Eq7m053673t1	3263
Funhe2EKm000010t1	Funhe2EKm000010	1	Class:Strong,Express:Strong,Homology:OrthologStrong,Intron:Strong-inerr,Map:Strong,Protein:complete	285	62%,858/1362	CD74 molecule, major histocompatibility complex, class II invariant chain (99%p)	HLA class II histocompatibility antigen gamma chain (99U)	CD74 molecule, major histocompatibility complex, class II invariant chain	platyfish:ENSXMAP00000008675,100%	na	FISH11G_G9283	CDD:238114,TrEMBL:UniRef50_Q6J613,TrEMBL:Q6J613_CHICK,	89%,16/18	100%,rx:adultenv,1000a/2e/200g	99%	Scaffold0:130568-136869:-	kf2x11gmap:Funhe2Exx11m019538t1,Fungr1EG3m013533t7	46773
Funhe2EKm000011t1	Funhe2EKm000011	1	Class:Strong,Express:Strong,Homology:OrthologStrong,Intron:Strong,Map:Strong-Split,Protein:complete	683	98%,2049/2090	NUAK family, SNF1-like kinase, 1 (99%p)	same	NUAK family SNF1 kinase 1	platyfish:ENSXMAP00000004773,100%	na	FISH11G_G1236.s3	CDD:214567,TrEMBL:UniRef50_O60285,TrEMBL:NUAK1_HUMAN,	100%,12/12	100%,rx:nodiff,4a/10e/3g	99%,Split2/2	Scaffold43:402782-424715:-/Scaffold0:133590-135018:kf2rae5:Funhe2Exx11m006564t1,Funhe2E6bm006583t2	0

grep -v '^#' kf2rae5g.main.attr.tbl2 | cut -f4 | sed 's/,.*//; s/partof:.*/partof/; s/-.*//;' | sort | uniq -c | sort -k1,1nr

tbl3b: names fixed, maybe.
23985 Class:Strong
7349 Class:StrongPartial
1899 Class:Transposon
1075 Class:Medium
  50 Class:MediumPartial
 524 Class:Weak
  15 Class:WeakPartial
 125 Class:None

tbl2b:
23985 Class:Strong
7349 Class:StrongPartial
1899 Class:Transposon   << check expressTE
1075 Class:Medium
  50 Class:MediumPartial
 526 Class:Weak         : fixed qual scores for Express/rnax, Homol/orthologs
  15 Class:WeakPartial
 123 Class:None

tbl2a:
23535 Class:Strong
7306 Class:StrongPartial
1899 Class:Transposon   << check expressTE
 929 Class:Medium       ? Weak vs Medium?
  70 Class:MediumPartial
1112 Class:Weak         ? ok or mistakes?
  38 Class:WeakPartial
 133 Class:None

Class:Weak includes
  -- should be Medium class at least, boost Homology and/or Express qual scores for p >= 33%, add high rnax score to expr qual
  Funhe2EKm000072t1, mayzebr:XP_004575089.1,34%   kf4bAUGarip1p1s0g82t1
  Funhe2EKm000205t1, tilapia:ENSONIP00000012903,61%; exp=27%,rx:nodiff,0.04a/0.1e/0.02g; kf4best:kf4bAUGerip2p4s0g14t1
  Funhe2EKm000246t1  tilapia:ENSONIP00000008417,38% exp=0%,rx:nodiff,0.06a/0.2e/0.1g; kf4best:kf4bAUGpia9cp1s1g33t1
  Funhe2EKm000505t1  UniProt:H3AY42_LATCH exp=32%,rx:grandis,3a/6e/20g <<
  
  
#KC	Dbxref	31967
#KC	ID	35022
#KC	Name	32320
#KC	Split	2250
#KC	aalen	34902
#KC	alttr	20084
#KC	cov	27939
#KC	cxlen	28848
#KC	ggname	32299
#KC	ggnamepct	32299
#KC	ggnameref	32299
#KC	homolog	30549
#KC	inerr	3002
#KC	inqual	29204
#KC	locus	34902
#KC	must	19379
#KC	name	32299
#KC	nameln	29879
#KC	namepct	32299
#KC	nameref	32299
#KC	nexon	25673
#KC	nintron	29204
#KC	oid	34902
#KC	oname	32299
#KC	onamepct	32299
#KC	onameref	32299
#KC	orclass	32299
#KC	ortholog	32299
#KC	osrc	34902
#KC	ovest	34897
#KC	ovpro	23807
#KC	ovrna	29020
#KC	rnax	34903
#KC	score	418
#KC	scoresum	28846
#KC	sense	939
#KC	srcf	22
#KC	tere	6161
#KC	uniprot	32299
#KC	upd2	568

=cut
