#!/usr/bin/env perl
# evigene/scripts/genes/evganngff.pl 

=item usage

  evganngff.pl  -genes genes.evgmap.gff -chr chr.fasta  [ -mrna genes.mrna -names genes.names ]

=item about

  inputs: genes.gff  : evigene format gff3
          genes.mrna : evigene mrna seq
          genes.names: naming table: ID \t Name \t NameIdent \t Dbxref ..
          genome.fasta : chromosome asm seq of gff
          genes.intronok.tab : valid introns/gene, upd1804, add mrna annot: nintron=xxx or
            
  outputs:  genes.evgann.gff3, suited to some gff3 viewers
  
  from anoph/geneval/findcds.info
  
=item see also 

  evigene/scripts/genefindcds.pl : does most of work
  evigene/scripts/genes/gmap2evgff.pl, assumed genes.gff format
  evigene/scripts/genes/genefixgaps.pl
  
=item fixme 
   -- change 'findcds' to option, nofindcds default: now evganngff.pl -findcds
   -- allow io-pipe: now evganngff.pl -names xxx.names -genes stdin  < evgmap.gff > evgann.gff
   -- merge w/  bestgenes_update.pl

   $evigene/scripts/bestgenes_update.pl -act=sort -debug -vers evg5at -conf evigene_pubgff.conf -in $pt.gff -out $pt.gff3c 

=cut 

use FindBin;
use lib ("$FindBin::Bin","$FindBin::Bin/../", "$FindBin::Bin/../../lib/"); #  evigene/scripts/genes/ path

use strict;
use warnings;
use Getopt::Long;

my $VERSION='2017.12.22'; # '2017.05.01'; 
our $EVIGENES="$FindBin::Bin/..";  

my $debug= $ENV{debug}||0; 
my $pFULL= $ENV{full}||0.75;
my $AGENE= $ENV{addgene}||0; # opt
my $mrnaType= $ENV{mrnatype} || "mRNA"; # "match" alt; others possible: 'mRNA|ncRNA|transcribed_RNA...'
# my $exonType= $ENV{exontype} || "exon"; # not needed?
# our $RNATYPES='mRNA|ncRNA'; # mrnatypes exontypes genetypes
# our $EXONTYPES='exon|CDS'; # mrnatypes exontypes genetypes

my $reTarget= 1; # Target= > trg= avoid hassles from other software

my($genesgff,$mrna,$genenames,$chrfasta,$dofindcds,$oid2pubid,$intronoktab); # was $cdsdone

my $optok= GetOptions(
  "genesgff=s", \$genesgff,
  "chrfasta|chromosomes=s", \$chrfasta, # dofindcds
  "mrna|cdna=s", \$mrna, # dofindcds
  "mrnatype=s", \$mrnaType, 
  "names|annot=s", \$genenames, # OR publicset/xxx.ann.txt
  "intronoktab|inoktab=s", \$intronoktab,  
  "pFULL=s", \$pFULL,
  "debug!", \$debug, 
  "reTarget!", \$reTarget, 
  "oid2pubid!", \$oid2pubid, # using ann.txt, change ID=OrigID to ID=PubID
  "findcds!", \$dofindcds, #was# "cdsdone|nofindcds!", \$cdsdone, 
  "addgene!", \$AGENE, 
);

die "usage: evganngff.pl  -genes genes.evgmap.gff -[no]findcds -chr chr.fasta [ -mrna genes.mrna -names genes.names ]"
  unless($optok and $genesgff and ($chrfasta or not $dofindcds));

$mrnaType =~ s/[,; ]+/\|/g; #$exonType =~ s/[,; ]+/\|/g;

my($ferr, $fcgff)= (0,0);
if($dofindcds){ ($ferr, $fcgff)= findcds($genesgff,$chrfasta,$mrna); die if($ferr); }
else { $fcgff=$genesgff; }

my $geneann= {}; my $oidpubids={};
if($genenames) {
  if($genenames =~ /ann\.txt/) { ($geneann,$oidpubids)= readannot($genenames,); }
  else { ($geneann)= readnames($genenames); $oidpubids=undef; } # allow mrna as genenames source?
}

my($nintrons, $inokerr, $introns)= 
  ($intronoktab) ? readIntronOk($intronoktab,$geneann) : (0);

my($anngff) = annotgff($fcgff,$geneann,$oidpubids);

#------------

sub findcds {
  my($genesgff,$chrfa,$mrna)=@_;
  
  my $genefixgff= $genesgff; $genefixgff=~s/\.gff.*//; $genefixgff.="_fcds.gff";
  $pFULL /= 100 if($pFULL>1);
  
  my $findcdsopt="-fixcds -completecds=2 -nogoodlen -full=$pFULL " . ($debug ? " -debug " : "");
  # outopts:  -outgff x.gff -outaa x.aa -outcds x.cds -outcdna x.cdna : want all??
  my $cmd="$EVIGENES/genefindcds.pl $findcdsopt -dna $chrfasta -genes $genesgff -cdna $mrna -outgff $genefixgff -outaa  ";
    # >& $pt-gasm$spp.fcds2.log

  warn "#CMD: $cmd\n" if($debug);
  my $err= system($cmd);
  warn "#ERR: $err, $EVIGENES/genefindcds.pl .. $genesgff " if($err);
  return($err,$genefixgff);
}

=item readannot

  add sub readannot : format of evg/publicset/xxx.ann.txt
publicset/evg3whitefly.ann.txt 
PublicID	OrigID	TrLen	CDSoff	AAqual	TrGaps	Dbxref	Namealign	Product_Name	CDD_Name
Bemtab3EVm005486t9	BemtabEGm000204t1	6704	1270-2871	533,23%,complete-utrbad	125	UniProt:C4G15_DROME,75,	0	Cytochrome P450 4g15	na
Bemtab3EVm008017t2	BemtabEGm000208t1	2240	372-1586	404,54%,complete-utrpoor	2	UniProt:E0VXV7_PEDHC,75,	0	Zinc binding dehydrogenase, putative	na
  
  * ^^Bug: Dbxref includes Namealign val
   UniProt:C4G15_DROME,75,<tab>0  >> UniProt:C4G15_DROME <tab> 75
=cut

sub readannot {  
  my($natab)= @_;
  return {} unless($natab and -s $natab);
  my($ok,%ann,%oidpubids,@hdr);
  if($natab=~/.gz$/){ $ok=open(F,"gunzip -c $natab|"); } else { $ok= open(F,$natab); }
  warn "#ERR: cant read annot table: $natab\n" unless $ok;
  while(<F>) {
    next if(/^\W/); chomp; my @v=split"\t";
    my ($pd,$oid,$trlen,$cdsoff,$aaqual,$trgaps,$nadx,$nap,$nam,$cddnam);
    if(/^PublicID/){ 
      @hdr= @v; next;
    } else { # check @hdr/%hdr field names.. may change
      map{ $_="" if($_ eq 'na'); } @v;
      ($pd,$oid,$trlen,$cdsoff,$aaqual,$trgaps,$nadx,$nap,$nam,$cddnam)= @v;
      if($nadx =~ /,\d/ and $nap<1) { ($nap)= ($nadx=~s/,(\d+)//)?$1:0; } # fix bug nadx,nap:
    }
    
    ## add separate hash fields for oid, trlen, cdsoff, aaqual, cddnam ?
    my $ann="";
    $ann .= "aalen=$aaqual;" if($aaqual);
    $ann .= "cdsoff=$cdsoff;" if($cdsoff); # trlen also, or expect it exists?
    if($nadx) { $ann .= "namealn=$nap;Dbxref=$nadx;Name=$nam;"; }
    elsif($nam =~ /\w/) { $ann  .= "Name=$nam;"; }
    $ann .= "cddname=$cddnam;" if($cddnam);
    $ann{$pd}= $ann||"";
    $oidpubids{$oid}=$pd;
  } close(F);
  return (\%ann,\%oidpubids);
}


sub readnames { # allow mrna hdr input ?
  my($natab)= @_;
  return {} unless($natab and -s $natab);
  my($ok,%ann);
  if($natab=~/.gz$/){ $ok=open(F,"gunzip -c $natab|"); }
  else { $ok= open(F,$natab); }
  warn "#ERR: cant read name table: $natab\n" unless $ok;
  while(<F>) {
    next if(/^\W/);
    chomp; my @v=split"\t";
    my ($td,$nam,$nap,$nadx,@more)= @v;
    my $ann=""; next unless($nam);
    if($nadx) { $ann="namealn=$nap;Dbxref=$nadx;Name=$nam;"; }
    elsif($nam =~ /\w/) { $ann= "Name=$nam;"; }
    $ann{$td}= $ann||"";
  } close(F);
  return \%ann;
}

sub annotgff {
  my($iname,$annh,$oidpubids)= @_;
  
  # mrna/gene drop annot : change these to config opts?
  my $MDROP='match|protein|indels|cdsindel|cdsspan|chimera|ocds|aaold|xtrim|offold|cdsfix|cdnabest'; #not qlen
  my $GDROP='Target|cov|indels|nexon|pid|path|oid|cxlen|aalen|cdsoff|utrx|chimera|chim[12]';
  my $haveoid2pub= (defined($oidpubids) and ref($oidpubids))?1:0;
  my $haveann=  (defined($annh) and ref($annh))?1:0;
  
  ## fixme: allow iopipe: annot < in.gff > out.ann.gff
  my($ok,$oname,$inh,$outh,%didg)=(0,0,undef,undef);
  my $ispipe=($iname=~/^-|stdin/i)?1:0;
  if($ispipe) { 
    $oname="STDOUT";  $inh=*STDIN; $outh=*STDOUT; $ok=1;
  } else {
    $oname= $iname; $oname=~s/\.gff.*//; $oname.="ann.gff3";
    if($ok= open(GIN,$iname)) { $inh=*GIN; 
      if($ok= open(GOUT,'>',$oname)) { $outh=*GOUT; }
    }
  }
  return (undef) unless($ok);
  
  my($pubid,$inid)=(0,0);
  print $outh "##gff-version 3\n";
  while(<$inh>) {
    next if(/^\W/); #?? maybe keep gmap2gff #i..intron' rows 
    #below# s/\%/p/g; # no % allowed w/o dang URLescaping
    my @v= split"\t";  
    if(@v == 9 and $v[4]=~/^\d/) { 
      if($v[2]=~/$mrnaType/) { # was $v[2]=~/mRNA/ # FIXME: other RNA types
        my $at=$v[8]; chomp($at); 
        my($td)= $at=~m/ID=([^;\s]+)/; 
        #old# $td =~ s/_C\d$//; ## _G nnn also ?!
        $td =~ s/_[CG]\d$//; 
        $inid= $td;  $pubid= $td; 
        if($haveoid2pub and (my $d= $oidpubids->{$td})) { $pubid=$d; $td=$pubid if($oid2pubid); }

        my($oid)= ($at=~m/;oid=([^;\s]+)/)?$1:""; 
        my @oid= map{ s/_C\d//; $_ } split",",$oid;

        $at =~ s/;($MDROP)=[^;\n]+//g; 
        if($at =~ m/;clen=\d/){ $at =~ s/;qlen=[^;\n]+//; } else { $at =~ s/;qlen=/;clen=/; } 
        #ugh: leave out: s/;Name=/;Note=/; s/;/;Name=$td;/; 
        #x $at =~ s/;aalen=(\d+),/;aalen=$1;aaq=/; # this is only for one viewer?
        $at =~ s/;(cdnaorf)=[^;\n]+//g; # leave Target ??
        
        (my $gd=$td)=~s/t\d+//; 
        if($AGENE) {
          unless( $didg{$gd}++ ) { 
            my @gv= @v;    $gv[8]=$at;
            $gv[2]="gene"; $gv[8]=~s/ID=$inid/ID=$gd/; 
            $gv[8]=~s/;($GDROP)=[^;\n]+//g; 
            print $outh join("\t",@gv)."\n"; 
          }
          $at =~ s/;/;Parent=$gd/; # mRNA add Parent= or gene=
        } 
        my $ann="";
        if($haveann) { 
          $ann= $annh->{$inid} || $annh->{$pubid};  
          unless($ann) { for (@oid) { last if( $ann= $annh->{$_}); } }
        }
        
        ## FIXME: new gmap2evgff retains Name,Dbxref,.. annots .. check/replace if differ, dont dupl
        if($at and $ann) {
          ($at)= checkAddAnn($td,$at,$ann);  
          $at=~s/=,/=/g; $at=~s/;;/;/g;  
          $v[8]="$at\n";  
        }
        $v[8]=~s/\%/p/g; # no % allowed w/o dang URLescaping
        if($oid2pubid and $inid ne $pubid) { $v[8]=~s/ID=$inid/ID=$pubid/; }
        
      } else {  # check for exon|CDS or other type
        my($pard)= $v[8]=~m/Parent=([^;\s]+)/; 
        $pard =~ s/_[CG]\d$//;   
        if($pard ne $inid) { } # error
        if($oid2pubid and $inid ne $pubid) { $v[8]=~s/Parent=$inid/Parent=$pubid/; }

        $v[8]=~s/\%/p/g; # no % allowed w/o dang URLescaping
        # drop Target= or no? >> change to trg=?
      } 
      if($reTarget) { $v[8]=~s/Target=/trg=/; }
      print $outh join("\t",@v); 
    } 
  } 
  close($inh); close($outh);

  return($oname);
}

## dang regex Name clashes; FIX: $v=~s/\W/./g; for s/$k=$v//
## Quantifier follows nothing in regex; marked by <-- HERE in 
#  m/\bName=Dihydropyrimidine dehydrogenase (NADP(+ <-- HERE ))-like protein[;]?/ 
# at /home/ux455375/bio/evigene/scripts/genes/evganngff.pl line 241, <STDIN> line 6342.

sub checkAddAnn {
  my($td,$gann,$tann)=@_;
  return($gann) unless($tann and $tann=~m/\w=/);
  $gann||="";
  ## check gff.ann for same as annh vals
  ## from ann.txt:  aalen cdsoff namealn Dbxref Name; fix cdsoff => offs 
  ## from gmapgff:  aalen offs Dbxref Name
  my (@tadd,@gdrop);
  my %gann=map{ my($k,$v)=split"=",$_,2; $k=>$v; } grep(/=/, split(";",$gann));   
  my %tann=map{ my($k,$v)=split"=",$_,2; $k="offs" if($k eq "cdsoff"); push @tadd,$k; $k=>$v; } grep(/=/, split(";",$tann));
  for my $tk (sort keys %tann) {
    my $tv=$tann{$tk};
    if($gann{$tk}) { if($gann{$tk} eq $tv) { delete $tann{$tk}; } else { push @gdrop,$tk; } }
  }
  for my $k (@gdrop) { my $v=$gann{$k}; $v=~s/\W/./g; $gann=~s/\b$k=$v[;]?//; }
  for my $k (@tadd) { if(my $v=$tann{$k}) {  $gann .=";$k=$v"; } }
  return $gann;
}

sub readIntronOk {
  my($inf,$annh)= @_;
  my(%inokerr, %introns); # == inokerr, introns
  my $ANTISENSE_INTRON_IS_ERROR= 0; # option now, as many daphnia cds appear valid but reverse of rna-intron strand
  my ($nin,$nerr,$nval,$ltd)=(0) x 9;
  my $ok= open(F,$inf); # gz?
  warn "#ERR: cant read intronok table: $inf\n" unless $ok;
  while(<F>) {
    next if(/^\W|^nogene/); 
    my @v=split;
    my($td,$ind,$inw,$iv,$xtype)=@v; # intronok.tab
    $xtype ||= "other"; # bug?
    
    my($vi) = ($iv=~m/valintron=(\d+)/)?$1:0; 
    my($as) = ($iv=~m/antisense=(\d+)/)?$1:0; 
    my($obe)= ($iv=~m/offby=(\d+)/)?$1:0; # upd1707, gsplign +1/-1 splice site when antisense
    my($ler)= ($iv=~/inlongerr=([1-9]\d*)/)?$1:0; 
    
    my $iflag="";
    $iflag.="${xtype}ok," if($vi>0 and not $ler); # exclusive of errs  
    if($as>0) { $iflag.= (($ANTISENSE_INTRON_IS_ERROR) ? "antierr," : "antisense,");  }
    if($obe>0) { $iflag.= (($ANTISENSE_INTRON_IS_ERROR) ? "offbyerr," : "offby,");  }
    $iflag.="longerr," if($ler>0);
    # $iflag ||= "noevd";
    
    $inokerr{$td} .= $iflag;
    #add count val/inval for gene annot nintron=ok/err ?
    #x $introns{$td}{$ind}= join"\t",$iflag,$inw,$iv,$xtype;
    $introns{$td}{'n'}++;
    if($vi>0) { $introns{$td}{'val'}++; } else { $introns{$td}{'val'}||= 0; }  # count only valid introns here now?
    $introns{$td}{'err'}++ if( $as or $ler );
    $nin++; $nerr++ if($iflag=~/err/); $nval++ if($iflag=~/ok/);
    
    $ltd=$td;
  } close(F);
  
  # fill %ann{id}{nintron}, other instead here?
  if(ref($annh)) {
    for my $td (sort keys %introns) {
      my $vin= $introns{$td}{'val'}||0;
      my $nin= $introns{$td}{'n'}||0;
      my $iflag= $inokerr{$td}||""; # == CDSok,,CDSok,,CDSok,,CDSok,,..
      my $ian= "nintron=$vin/$nin"; # to match overbest need, nin ~= nexon, what of err in?
      #OR: my $nxn=$nin+1; my $ian= "inexon=$vin/$nxn/$nin";
      # $ian.=";inflag=$iflag" if $iflag;
      (my $tdc=$td) =~ s/_[CG]\d$//; 
      my $pubid= $tdc; 
      # if($haveoid2pub and (my $d= $oidpubids->{$tdc})) { $pubid=$d; $tdc=$pubid if($oid2pubid); }
      my($ann,$anid)=("","");  
      if( $ann= $annh->{$pubid} ) { $anid=$pubid; } 
      elsif( $ann= $annh->{$tdc} ) { $anid= $tdc; }  
      else { $ann= $annh->{$td}||""; $anid= $td; } 
      $ann.=";" if($ann); $ann.=$ian; $annh->{$anid}= $ann;
    }
  
  }
  
  warn  "# readIntronOk n=$nin,ok=$nval,errs=$nerr from $inf\n" if $debug;
  return($nin, \%inokerr, \%introns);
}
 
__END__

=item findcds.info

anoph/geneval/findcds.info

# findcds to replace poor gmap cds ; genefindcds == genefindcds2
# * gmap mappings are worse than blastn for many.
# BUT -fixcds forces match to genome, can chop large portion of aa/cds that doesnt map
# -completecds looks for partial ends on genome;

$evigene/scripts/genefindcds.pl -fixcds -debug -nogoodlen -full=0.75 -completecds=2 \
 -dna genome/gasm$spp.fa -genes $pt-gasm$spp.gmap.gff -cdna $pt.mrna.gz \
 -outgff $pt-gasm$spp.fcds2.gff -outaa >& $pt-gasm$spp.fcds2.log

# annotate with Name/Note/..

gunzip -c $pt.fcds2.gff.gz | cat  $pt.nametab - | perl -ne \
'next if(/^\W/); @v=split"\t";  
if(@v == 9 and $v[4]=~/^\d/) { if($v[2]=~/mRNA/) { $at=$v[8]; $at=~s/;($DRAT)=[^;\n]+//g; ($td)=m/ID=(\w+)/; 
$td=~s/_C\d$//;  $ann=$ann{$td}; chomp($at); $at.=";$ann" if $ann; $at=~s/;;/;/g;  $v[8]="$at\n";  } 
print join("\t",@v); } 
elsif($v[0]=~/^(Aed|Ano)/) { ($td,$oid,$cla,$aw,$dbx,$nam)=@v; chomp($ann); $ann{$td}="Dbxref=$dbx;Name=$nam;"; } 
BEGIN{ $DRAT="match|protein|offs|cdsindel|qlen|cdsspan|ocds|aaold|xtrim|offold|cdsfix|cdnabest"; } ' \
 > $pt.fann.gff

# convert to viewable gff3

env addgene=0 perl -ne 'BEGIN{ print "##gff-version 3\n"; $AGENE=$ENV{addgene}||0; } 
s/\%/p/g; @v=split"\t"; if(/\tmRNA/) { ($td)=m/ID=(\w+)/; ($gd=$td)=~s/_C\d//; s/;Name=/;Note=/; s/;/;Name=$td;/; 
$gd=~s/t\d+//; s/;aalen=(\d+),/;aalen=$1;aaq=/; s/;(Target|cdnaorf)=[^;\n]+//g; 
if($AGENE and not $did{$gd}++) { @gv=split"\t"; $gv[2]="gene"; $gv[8]=~s/ID=$td/ID=$gd/; $gv[8]=~s/Name=$td/Name=$gd/;
$gv[8]=~s/;(Target|cov|indels|nexon|pid|path|oid|cxlen|aalen|cdsoff|utrx|chimera|chim[12])=[^;\n]+//g; 
print join("\t",@gv); s/;/;Parent=$gd/; } s/=,/=/; } s/;Target=\w+ \d+ \d+//;  print; ' \
  $pt.fann.gff > $pt.fann.gff3

# anofunevg24m.hoho.fann.gff > anofunevg24m.hoho.fann.gff3
# evg12aedes-gasmaedes-sc575a.fann.gff3  > evg12aedes-sc575b.fann.gff3
