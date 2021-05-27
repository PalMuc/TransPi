#!/usr/bin/env perl
# pubseqgff.pl for evg2tribol; merge w/ pubseq.sh
# see bugs/beetlet/tcas4evg/beetle_tcas4evg1.info

=item notes

 work from keep,cull annotation tables:
 evg2tribol/publicset/evg2tribol.ann2.keep
 PublicID	OrigID	TrLen	CDSoff	AAqual	TrGaps	MapCov	MapIdn	MapInExon	MapLocus	MapPath	DbXref	OGenes	NamePct	ProductName
 Tribca2aEVm000001t1	tcas4sb2p8nmvelvk47Loc2069t1	56083	349-55173	18274,97%,complete	0	100	9988/93	tcas3NC_007424:13110520-13178971:-	0	RefID:ARP7f_G2613,dromel:FBgn0053196,CDD:249587,CDD:224280,CDD:240526,CDD:238011/2,CDD:236766	tcas3nc:XP_008197897/90	100%,19215/11023,18624	Neurogenic locus notch protein
 Tribca2aEVm000001t2	tcas4cr2p8nmvelvk47Loc1388t12	55135	76-54210	18044,98%,complete	0	100	9989/94	tcas3NC_007424:13110509-13230420:-	0	RefID:ARP7f_G2613,dromel:FBgn0053196,CDD:249587,CDD:238011/17,CDD:214542,CDD:236304,CDD:240526,CDD:251869,CDD:201524	tcas3nc:XP_008197897/77	100%,16323/10279,18003	Neurogenic locus notch protein

=cut

use strict;

my $debug=$ENV{debug}||1;
my $vclass=$ENV{class}||"main"; # alt frag|cull
my $evgname=$ENV{name}||"evg2tribol";
my $VUPD=$ENV{upd}||"2";
my $IDPRE=$ENV{idpre}||"Tribca2aEVm"; 
 
my $pubsrc="tcas4evg${VUPD}loc"; 
my $oname="$evgname.fin${VUPD}loc"; 
my $idtab="$evgname.ann2.keep";
my $gmapgff="$evgname.gmap.gff.gz";

if($vclass eq "frag" or $vclass eq "cull") {
  ## ** Dont put drops, only cull/demote class
  #$idtab="`ls $evgname.ann2.cull[1235]`"; # FIXME..
  $idtab="$evgname.ann2.allculls";  
  $oname=~s/loc/cull/; ## cull better than frag now
  $pubsrc=~s/loc/cull/;

} elsif($vclass eq "alt") {
  $idtab="$evgname.ann2.keep";
  $oname=~s/loc/alt/;
  $pubsrc=~s/loc/alt/;

} elsif($vclass eq "main") {
  #$idtab="$evgname.ann2.keep";
  #$oname="$evgname.$tabtag";
  #$pubsrc="tcas4evg$VUPD"."loc"; 
} else {
  die "I got no class for: $vclass \n"; # or main default
}


#------------------

my $VCLA=$vclass; # $VCLA=$ENV{vclass}; $IDPRE=$ENV{idpre};
my ($vd,$ok,$nok)=(0) x 9;
my (%did,@cols,%vid,%vod,%id,%puban);

#.. id class and annotations for keep/cull items ..............
warn "# vclass=$vclass read $idtab\n" if($debug);
open(F,$idtab) or die $idtab; 
# FIXME: also print input idtab == annot.txt for each vclass to 
open(O,">$oname.ann.txt");  

while(<F>) { 
  chomp; my @v= split"\t"; 
  if(/^PublicID/) {  @cols= @v;  
    splice(@v,2,0,"ClassV");  print O join("\t",@v)."\n"; 
    next;
  } elsif(/^$IDPRE/) {
  my($vd,$oids,$aaq)=@v[0,1,4,]; my $vupd=$VUPD;
  ## PublicID	OrigID	TrLen	CDSoff	AAqual	TrGaps	MapCov	MapIdn	MapInExon	MapLocus	MapPath	DbXref	OGenes	NamePct	ProductName

  my($vg,$vi)= $vd=~m/(\w+t)(\d+)$/;
  my $vact=($VCLA=~/frag|cull/)?"cull":($vi == 1)?"main":"alt";
  
  if($vact=~/^drop/) { next; }
  elsif($vact=~/^(cull|demote)/i) { next unless($VCLA=~/^frag|cull/); }
  elsif($vact=~/^(alt)/i) { next unless($VCLA=~/^alt/); }
  else { next if($VCLA=~/^(alt|frag|cull)/); }
  $vid{$vd}="$vact,$vupd"; $vod{$vd}=$oids; $id{$vd}=$vd; map{ $id{$_}=$vd; } split",",$oids; 

  my $an=""; for my $i (0..$#cols) { my $c=$cols[$i]; if($c =~ /DbXref|ProductName|OGenes/){ 
    my $v=$v[$i]; $an.="$c=$v;" if($v and $v ne "na"); } }
  $an=~s/ProductName/Name/; $an=~s/OGenes/oGenes/; 
  $puban{$vd}=$an;
  
  splice(@v,2,0,"$vact,$vupd"); print O join("\t",@v)."\n";
  }

} close(F); close(O);
$nok= scalar(keys %vid); warn "# vclass=$vclass annots nok=$nok in $idtab\n" if($debug);

exit(0) if($ENV{onlyann});

#.. pubseq ..............
## eg1:
# >Tribca2aEVm007098t1 type=protein; aalen=370,78%,complete; clen=1411; offs=165-1277;
#  oid=tcas4sb2g9sr118g9velvk65Loc5373t3; organism=Tribolium_castaneum; 
#  evgclass=althi1,okay,match:tcas4cr1c6soapk25loc30507t2,pct:99/98/./tcas4cr2p8nmsoapk45loc23381t4;; 
#  DbXref=RefID:ARP7f_G4281,daphplx:hxAUG26us142g174t1,CDD:238105,CDD:253361;oGenes=tcas3nc:XP_008198342/100;Name=Syntaxin-1A;; vers=main,1

my $facut="evgclass";
foreach my $st (qw( aa cds mrna)) {
  my $inf="$evgname.${st}_pub.fa.gz";
  warn "# vclass=$vclass seq.$st read $inf\n" if($debug);
  open(F,"gunzip -c $inf |") or die "open $inf";
  open(O,">$oname.$st");  
  $nok=0;
  while(<F>) {
   if(/^>(\S+)/) { $ok=0; my $td=$1; my($od)=m/oid=(\w+)/; $vd= $id{$td}||$id{$od}; 
     if($vd and not $did{$st}{$vd}) { 
      if($td ne $vd) { s/>$td/>$vd/; unless(s/oid=/oid=$td,/) { s/$/ oid=$td;/; } }
      ## drop this? evgclass=althi1,okay,match:tcas4cr1c6soapk25loc30507t2,pct:99/98/./tcas4cr2p8nmsoapk45loc23381t4;
      s/\s+($facut)=[^;\n]+[;]?//g; s/; *$//;
      if(my $ann=$puban{$vd}) { $ann=~s/;/; /g; $ann=~s/; $//; s/$/; $ann/; }
      if(my $vid=$vid{$vd}) { s/$/; vers=$vid/; }
      $did{$st}{$vd}++; $ok=1; $nok++;
      } 
    } 
    print O if $ok;  
  } close(F); close(O);
  # $nok= scalar(keys %{$did{$st}}); 
  warn "# vclass=$vclass seq.$st nok=$nok to $oname.$st\n" if($debug);
}

#.. pubgff .............
($vd,$ok,$nok)=(0) x 9;
my $st="gff";
my $acut="trg|Target|match|chimera|cdnabest|nover|cdsover|ocds|svec"; 
my $aspl="aalen|offs|oid"; 
my($isp,%sat);

warn "# vclass=$vclass seq.$st read $gmapgff\n" if($debug);
open(F,"gunzip -c $gmapgff |") or die "open $gmapgff";
open(O,">$oname.gff");  
while(<F>) {
  if(/^\W/){ print O unless(/^#i/); next; } # next if(/^\W/);
  my($ref,$src,$ty,$rb,$re,$rv,$ro)=split"\t";

  # splits
  $isp=0; 
  if(/\t(ID|Parent)=(\w+)_C([12])/){ 
    my($at,$d)=($1,$2); $isp=$3;
    s/${d}_C$isp/$d;Split=$isp/; 
    if(/mRNA/) { my $sat;
      if($sat=$sat{$d}) { s/aalen=\d+;//; s/$/;$sat/; } 
      else { $sat=join";", m/((?:$aspl)=[^;\s]+)/g; $sat{$d}=$sat if($sat); } 
    }
  }
  
  if($ty eq "mRNA") { $ok=0;
    my($td)=m/ID=(\w+)/; my($od)=m/oid=(\w+)/; 
    $vd= $id{$td}||$id{$od}; $ok=($vd)?1:0; 
    # FIXME td_C[12] > td;Split=[12] .. other gmap.gff fixes....
    if($ok) { 
      s/ID=/ID=$vd;sid=/ if($td ne $vd); $nok++; 
      my @al= m/;(aalen=)/g; if(@al>1) { s/aalen=\d+;//; } 
      s/;($acut)=[^;\n]+//g;
      
      if(my $ann=$puban{$vd}) { s/$/;$ann/; }
      if(my $vid=$vid{$vd}) { s/$/;vers=$vid/; }
      $did{$st}{$vd}++;
    }
  } elsif($ok and /Parent=(\w+)/) { 
     my $pd=$1; s/Parent=$pd/Parent=$vd;sid=$pd/ if($vd ne $pd); 
     s/;($acut)=[^;\n]+//g;
  }
  if($ok) { s/\t$src/\t$pubsrc/; print O $_; }
} close(F); 
close(O);
warn "# vclass=$vclass seq.$st nok=$nok to $oname.$st\n" if($debug);

## fixme: want %did for each of seqtype and gff...
# foreach my $st (qw( aa cds mrna gff)) {
# open(O,">$oname.$st.missids"); 
# for $vd (grep{not $did{$st}{$_}} sort keys %vid) { $vo=$vod{$vd}||"noid"; print O "$vd\t$vo\n"; } 
# close(O);  }
  

__END__

