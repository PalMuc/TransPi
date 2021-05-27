#!/usr/bin/env perl
# kfish2_geneclass.pl

=item run

 env note=1 mapok=1 perl kfish2_geneclass.pl ...
 env dotab=0 hoval=10 note=0 mapok=0 $kfish2/submitf/pubgenes/kfish2_geneclass.pl \
  $kfish2/submitf/pubgenes/kfish2rae5h.main.pub.cds.qual \
  $kfish2/genes/an5d/score/kfish2rae5g.orpar.kvtab \
  $kfish2/ncbifish/kfishloc2fishmapngene_ordekas.tab \
  $kfish2/submitf/pubgenes/kfish2rae5g.main.attr.tbl5 

=item outclassn

env crel=1 dotab=0 hoval=20 note=0 mapok=0 $kfish2/submitf/pubgenes/kfish2_geneclass.pl  \
   $kfish2/submitf/pubgenes/kfish2rae5h.main.pub.cds.qual  \
   $kfish2/genes/an5d/score/kfish2rae5g.orpar.kvtab  \
   $kfish2/ncbifish/kfishloc2fishmapngene_ordekas.tab  \
   $kfish2/submitf/pubgenes/kfish2rae5g.main.attr.tbl5

# nloci        :	34924
# method_source:	4904 gmodel,	5 othr,	27512 rnasm,	2503 v1mix,
# ortho_class  :	3679 inpar,	7691 nocl,	21122 orlog,	2432 upar,
# homology_val :	30231 protein,	31171 relgenomes,	25528 sigKaKs align>=20
# gene_class   :	761 CodingNew,	30894 CodingOld,	436 NoncodeNew,	697 NoncodeOld,	2136 TE,
# expression   :	34518 FPKM>=0.02,	28646 FPKM>=1
# gmapping     :	25716 ok>=80%,	2251 split-map,	4678 partial,	2267 none 
# single-exon  :	5692,	10272 introns<=1 
# cds%(complete):	15647 cd>65%,	7846 cd>=40%,	3503 cd<40% of 26996, or	5750 cd<50%

  variants of homology alignment val (note all are significant blastp/blastn evalues)
# homology_val :	28497 protein,	29975 relgenomes,	25528 sigKaKs, for align>=33
# gene_class   :	1090 CodingNew,	30448 CodingOld,	608 NoncodeNew,	642 NoncodeOld,	2136 Transposon,

# homology_val :	31602 protein,	32083 relgenomes,	25528 sigKaKs, for align>=10
# gene_class   :	529 CodingNew,	31229 CodingOld,	275 NoncodeNew,	755 NoncodeOld,	2136 Transposon,

=item inputs

 inputs: .. add relfish align/kaks tab, ...
 $kfish2/pubdir/genes/kfish2rae5/kfish2rae5g.main.pub.aa.qual  
 $kfish2/genes/an5d/score/kfish2rae5g.orpar.kvtab 
 kfish2/pubdir/genes/kfish2rae5/kfish2rae5g.main.attr.txt 
    or $kfish2/submitf/pubgenes/kfish2rae5g.main.attr.tbl5 
 add $kfish2/genes/an5d/score/kfish2rae5g.main.homolog.kvtab for all homol scores, 
 or use attr fields? ortholog, paralog?, Dbxref=CDD,UniProt if missing ortholog

# add ncbifish/kfishloc2fishmapngene_ordekas.tab : key = '(orlog|inpar|upar|uniq|noorde),'
#QueryID	camolly.pal	camolly.paloc	cguppy.pal	cguppy.paloc	cturkf.pal	cturkf.paloc	czfish.pal	czfish.paloc	amolly.cds	guppy.cds	turkf.cds	cmolly.kas	cguppy.kas	cturkf.kas	orlogdtrlen	aalen
#Funhe2EKm000003	77	77	92	0	100	100	28	28	96	95	86	0.163,0.0	0.204,0.0	0.115,0.0	orlog,funct,node,ncbi	679	214,94%,complete
#Funhe2EKm000004	100	100	100	100	90	82	48	40	97	98	88	0.185,0.0	0.214,0.0	0.135,0.0	orlog,funct,node,ncbi	5255	1103,63%,complete
#Funhe2EKm000005	98	98	98	98	91	91	41	41	99	99	98	0.129,0.0	0.137,0.0	0.072,0.0	orlog,funct,node,ncbi	4524	1010,67%,complete

# add  codingpot score in pubgenes/kfish2rae5h.{main,alt}.pub.cds.qual

# add genes/an5d/rnaxf/rxkf2mrna5f_rna3sample.tpm, recalc TPM score for all samples, 
     from genes/an5d/rnaxf/rxkf2mrna5f_rna.xcount3
  738 have zero reads, TPM<0.01  1377, <0.02 2033 loci, <0.1 6228 
 15512 TPM>=1 , 20184 TPM>= 0.5, 28675 TPM>=0.1, 32870 TPM>=0.02, 34357 TPM>0
 22054 TMP>=1 for 1 samp; 19282 TPM>=1 for 2 samples .. problem is combining samples?
 ?? is genes/an5d/rnaxf/rdkf2cds5g_rna.xgenec3 best calc, max fpk of 3 samples 
 gene_id	length	wh.ti	wh.tot	wh.fpk	wh.oid	md.ti	md.tot	md.fpk	md.oid	gr.ti	gr.tot	gr.fpk	gr.oid
   28774 FPKM>=1, 31573 FPKM>=0.5, 33947 FPKM>=0.1, 34213 FPKM>=0.02 .. near same as for attr.txt score

=cut

#BEGIN
$MAPOK=$ENV{mapok}||0;
$noTE=$ENV{note}||0; 
$DOID=$ENV{doid}||0;
$DOCLASSTAB=$ENV{dotab}||0;
$MINHOVAL=$ENV{hoval}||10;
$CRELCLASS= $ENV{crel}||0;
$NOBAD=1;

# nloci        :	34924
# method_source:	4904 gmodel,	5 othr,	27512 rnasm,	2503 v1mix,
# ortho_class  :	3679 inpar,	7691 nocl,	21122 orlog,	2432 upar,
# gene_class   :	1140 CodingNew,	30618 CodingOld,	1030 Noncoding,	2136 TE,
# expression   :	34518 FPKM>=0.02,	28646 FPKM>=1
# gmapping     :	25716 ok>=80%,	2251 split-map,	4678 partial,	2267 none 
# single-exon  :	5692,	10272 introns<=1 
# cds%(complete):	15647 cd>65%,	7846 cd>=40%,	3503 cd<40% of 26996, or	5750 cd<50%
# /bio/bio-grid/kfish2/submitf/pubgenes/kfish2_geneclass.tab1

# MINHOVAL effect:
# hoval=10 gene_class   :	1321 CodingNew,	30619 CodingOld,	665 Noncoding,	2136 TE,	184 Unknown,
# hoval=20 gene_class   :	1532 CodingNew,	30313 CodingOld,	712 Noncoding,	2136 TE,	232 Unknown,
# codep CF/CH:Noncode call, hoval=10 .. use this? are Unknowns really unclassifiable?
# gene_class   :	1140 CodingNew,	30619 CodingOld,	846 Noncoding,	2136 TE,	184 Unknown,



while(<>) { # change to per input file?
  ($td)= @v =split"\t"; 
  next unless(/^Funhe2EKm/);

  if($v[1]=~/ocla=(\w+)/) { # orpar.kvtab
    $ocl=$1; $ocl=~s/\d+//; $orpar{$td}=$ocl unless($orpar{$td});  # 1st is best
    } 
    
  elsif($v[1]=~/\tortholog=(\w+)/) { # kfish2rae5g.main.homolog.kvtab, same as attr.tab/ortholog?
    $hoval=$1; ($hop,$hov,$hoid)=split",",$hoval; 
    if($hoid and $hop>0) { $hoval{$td}="$hoid,$hop" unless($hoval{$td}); }
    } 
    
    # many cols:
  elsif($v[1]=~/^\d/ and /\t(orlog|inpar|upar|uniq|noorde),/) { # kfishloc2fishmapngene_ordekas.tab
    @crel=@v[1,3,5,7];  # palign/paloc for molly,guppy,turkf,zfish; use only palign
    @grel=@v[9,10,11];  # gene align for m,g,t
    @ckas=@v[12,13,14]; # chr kaks mol,g,t
    @orde=@v[15..17];   # orlog,func,.. \t trlen \t aalen
    $gd=$td; $td=$gd."t1";
    $camax=0; map{ $camax=$_ if($_>$camax) } @crel; 
    $crel{$td}= $camax; 
    ($kmax,$kmaxp)=(0,1); 
    map{ ($kas,$pk)=split","; if($kas>0 and $pk<$kmaxp) { $kmaxp=$pk; $kmax=$_; } } @ckas;
    $ckas{$td}= $kmax if($kmaxp<=0.05); # only signif kaks pk<=0.05 ??
    }
    
  elsif($v[1]=~/^\d/ and $v[2]=~/^\d/) { #  aa or cds.qual
    my($td,$al,$gp,$aaw,$tw,$cp)=@v; 
    $aaw{$td}=$aaw; # aa.qual ** replace w/ cds.qual including Codingpot?
    if($cp=~/Code|Noncod|Unkn/) { $codep{$td}=$cp; }
    }
    
  elsif(/^Fun/ and @v>15) { # main attr table, kfish2rae5g.main.attr.tbl5 == 
    ($td,$gd,$ti,$ql,$aw,$cdw,$nam,$onam,$gnam,$orlog,$plog,$omcl,$dbx,$inx,
      $rxp,$gcov,$loc,$kf1id,$oid,$vscore)=@v; 
    
    $hoval{$td}=$orlog unless($hoval{$td});
    
    # TE class
    $isTE=($ql=~/Transposon/ or $nam=~/TE:/)?"TE":0; 
    next if($isTE and $noTE>0); 
    next if($noTE<0 and not $isTE); 

    # prot quals
    $aaw=$aaw{$td}||"";
    $codep=$codep{$td}||"nocp"; $codep=~s,/[\d\.-]+,,g; ##$codep=~s/,.*//; 
    if($aaw) { $isbad=0; $ok{$td}=$_; } else { $isbad=1; $bad{$td}=$_;  }
    next if($NOBAD and $isbad); 
    
    # gmap quals
    ($pcov,$xcov)=split",",$gcov; $pcov=~s/%//; 
    next if($MAPOK and $pcov<80); 
    
    $gmap=0;
    if($loc=~m,/Scaf, or $gcov=~/Split/) { $gmap="split"; $gmapsp++; }
    elsif($pcov>=80) { $gmap="ok80"; $gmapok++; } elsif($pcov>1) { $gmap="part";  $gmaplo++; }
    elsif(/NOPATH/) { $gmap="none"; $gmapnone++; } 
    # $gmapv{$gmap}++;
    
    ($pix,$inn,$xnn)=split/[,\/]/,$inx; # intron/exon
    $xone= ($xnn==1)?1:0;
    $onex++ if($xnn<2); 
    $nonei++ if($inn<1 and $xnn>0); # was inn<2 .. only count 0 introns
    $onei++ if($inn>0 and $xnn>0); # was inn<2 .. only count 0 introns
    
    # xpression score, should this be sum of xa+xe+xg? NO, use max of 3 samples
    ($xco,$xde,$xa,$xe,$xg)=split/[,\/]/,$rxp; 
    $xpkm=0; map{ s/[aeg]//; $xpkm=$_ if($_>$xpkm); } ($xa,$xe,$xg); 
    $fpkm2++ if($xpkm>=0.02); $fpkm1++ if($xpkm>=1); # add fpkm01 if($xpkm>=0.1) == 33947
    
    # method: trasm, gmodel, v1mix 
    $meth="othr";
    if($oid=~/Exx11/) { $meth="rnasm"; } 
    elsif($oid=~/Funhe5EG/) { $meth="v1mix"; }  
    elsif($oid=~/AUG/){ $meth="gmodel";  } 
    else{ $meth="othr"; $gmoda{$td}=$_; } # 5 of these, are ref protein map models
    $gmod{$meth}++; 
    
    # homology quals
    $orpar= $orpar{$td}||"noor";
    $hoval= $hoval{$td}||"noho"; 
    if( ($hovalid,$hovalp)=split",",$hoval ) {
      ($hovalspp=$hovalid) =~ s/:.*//;
      $hovalp||=0;
      $hoval="$hovalp,$hovalspp";
      $hovalp=~s/%//; 
    }
    $crel=$crel{$td}||0; $ckas=$ckas{$td}||0;
    $thoval++ if($hovalp>=$MINHOVAL);
    $trel++ if($crel>=$MINHOVAL);
    $tkas++ if($ckas);
 
    # prot quals
    # $aaw=$aaw{$td}||"";
    # $codep=$codep{$td}||"nocp"; $codep=~s,/[\d\.-]+,,g; ##$codep=~s/,.*//; 
    # if($aaw) { $isbad=0; $ok{$td}=$_; } else { $isbad=1; $bad{$td}=$_;  }
    # next if($NOBAD and $isbad); 
   
    # prot quals,2
    ($acomp)= ($aaw=~/(complete|partial)/)?$1:"pother";
    ($pcds)= $cdw=~m/^(\d+)%,/; # code/utr score
    if($pcds>0 and $aaw=~/complete/) { 
      $pc5++ if($pcds<50); 
      if($pcds>65){ $pc66++; } elsif($pcds>=40){ $pc40++; } else{ $pc1++; } 
      }
    
    # Gene class rule
    # test also crel align >= MINHOVAL ? allow NoncodeOld class?
    my($tclass,$oldcode,$newcode,$nocode)=(0) x 9;
    if($isTE) { $tclass="Transposon"; } # "ExpressedTransposon";  # $isTE;  # this excludes all others? exclude orlog?
    elsif($orpar=~/orlog/) { $tclass=$oldcode="CodingOld1"; } # Old1,2,3,3 types?
    elsif($orpar=~/inpar/) { $tclass=$oldcode="CodingOld2"; } #  
    elsif($hovalp>=$MINHOVAL) {  $tclass=$oldcode="CodingOld3"; } # change this from 10%?
    elsif($ckas) { $tclass=$oldcode="CodingOld4"; }  
    # .. does $crel align to relfish affect gene class? yes
    elsif($CRELCLASS and $crel >= $MINHOVAL) {
      if($codep =~ /^Noncode/ or ($codep =~ m/CF:Noncode/ and $codep =~ m/CH:Noncode/)) { $tclass=$nocode="NoncodeOld"; }
      elsif($codep =~ /^Code/) { $tclass=$olcode="CodingOld5"; }
      elsif($codep =~ m/(CF|CH|CA):Noncode/) { $tclass=$nocode="NoncodeOld"; }
      else { $tclass=$oldcode="CodingOld5"; }
    }
    ## new loci, using codepot Code|Noncode|Unknown
    elsif($codep =~ /^(Code|Noncod)/) {
      if($codep =~ /^Noncode/ or ($codep =~ m/CF:Noncode/ and $codep =~ m/CH:Noncode/)) { 
        $tclass=$nocode="NoncodeNew";
      } 
      elsif($codep =~ /^Code/) { $tclass=$newcode="CodingNew"; }
    }
    else { 
      if($codep =~ m/(CF|CH|CA):Noncode/) { $tclass=$nocode="NoncodeNew"; }
      else { $tclass="Unknown"; } # what are unclassifiable?
    }
    
   if($tclass =~ /^Noncod/ and not ($codep =~ /^Noncod/) ) { $codep =~ s/^\w+/Noncode/; }

#.. Unknowns: use (CF|CH|CA):Noncode, get all classified
# eg Unk class .. are hard to classify
#Funhe2EKm000020t1	m:v1mix	x:0.1	o12:noor	o3:0,na	o4:79	o5:0	g1:ok80	ga:100	p12:partial	p3:40p4:Unknown,CF:Unknown,CH:Code,CA:Noncode	t:Unknown	t1:0	t1a:0	t2:0	t3:0	t4:0
#Funhe2EKm000077t1	m:gmodel	x:0.6	o12:noor	o3:0,na	o4:64	o5:0	g1:ok80	ga:100	p12:complete	p3:36	p4:Unknown,CF:Unknown,CH:Code,CA:Noncode	t:Unknown	t1:0	t1a:0	t2:0	t3:0	t4:0
#Funhe2EKm001561t1	m:rnasm	x:3	o12:noor	o3:0,na	o4:0	o5:0	g1:ok80	ga:96	p12:partial	p3:67p4:Unknown,CF:Unknown,CH:Code,CA:Noncode	t:Unknown	t1:0	t1a:1	t2:0	t3:0	t4:0
    
    $tclass =~ s/CodingOld./CodingOld/; # dont  need 1234 here.
    $tclass{$tclass}++; 
    
    $nok++; 
    if($DOCLASSTAB) {
      # table class quals: m.meth, x.xpkm, 
      #  o1,2.orpar + o3.hoval + o4.relf.align/kaks, 
      #  p.acomp/pcds/codep
      #  g.gmap,pcov
      #  t2.isTE, t1a.$xone, t3,4: codep - hoval
      $simtab= ($DOCLASSTAB>1)?1:0;
      #?? merge t2/t3/t4? same as t:Class
      #?? separate combo vals?  o3: palign,pspp  codep:Code,CF/CH/CA p4CodPot
      $codept=$codep; 
      $hovalt=$hoval; $o3hd="o3ProtHomology";
      $ckast=$ckas;

      $gmap=~s/ok80/okay/;
      if($simtab) {
	$codept=~s/,.*//; 
        $hovalt=~s/,/\t/;  $o3hd="o3ProtHomology\to3Species";
        $ckast||="NA"; $ckast=~s/,.*//;
      }

      @val=($td,"m:$meth","x:$xpkm",
        "o12:$orpar","o3:$hovalt","o4:$crel","o5:$ckast",
        "g1:$gmap", "ga:$pcov",
        "p12:$acomp","p3:$pcds","p4:$codept",
        "t:$tclass", "t1a:$xone"); 
        #old "t1:$oldcode", "t1a:$xone","t2:$isTE","t3:$nocode","t4:$newcode");
      map{ s/^\w+://; } @val; 

      unless($header++) {
        @hd=(qw(GeneID Method xFPKM o12Ortho), $o3hd, qw( o4GenomHomol o5KaKs
          g1GenoMap gGenoCov p12Complete p3CDS p4CodePot
          tClass t1aOneExon ));
          #old: tClass t1CodeOld t1aOneExon t2Transposon t3Noncode t4CodeNew );
        print join("\t",@hd)."\n";
      }
      print join("\t",@val)."\n"; 
    }
  } # attr.tab
}


#END{
print "# nloci        :","\t$nok\n"; 
print "# method_source:";  map{ print "\t$gmod{$_} $_,"; } sort keys %gmod; print "\n";
my %ocl=(); for $td (sort keys %ok) { $ocl=$orpar{$td}||"nocl"; $ocl{$ocl}++; }
print "# ortho_class  :"; map{ print "\t$ocl{$_} $_,"; } sort keys %ocl;  print "\n";
print "# homology_val :","\t$thoval protein,\t$trel relgenomes,\t$tkas sigKaKs, for align>=$MINHOVAL\n";
print "# gene_class   :";  map{ print "\t$tclass{$_} $_,"; } sort keys %tclass; print "\n";
print "# expression   :","\t$fpkm2 FPKM>=0.02,\t$fpkm1 FPKM>=1\n"; 
print "# gmapping     :","\t$gmapok ok>=80%,\t$gmapsp split-map,\t$gmaplo partial,\t$gmapnone none \n";
print "# single-exon  :","\t$onex,\t$nonei none / $onei some evidence introns \n"; 
## add protein complete/part qual counts
my $pct=$pc66+$pc40+$pc1; 
print "# cds%(complete):\t$pc66 cd>65%,\t$pc40 cd>=40%,\t$pc1 cd<40% of $pct, or\t$pc5 cd<50%\n"; 

# optional id dump
if($DOID){ @id=sort (keys %ok, keys %bad); print join("\n",@id)."\n"; } 

#} 


__END__

=item outs

nok=32789 noTE
aug	4088
othr	5
rna	26415
ver1	2281
ocl inpar	3180
ocl nocl	6815
ocl orlog	20801
ocl upar	1992
FPKM>=0.02 32496, FPKM>=1 27485
onex=4946, onei=8478
gmap ok80: 24003, split: 2145, part: 4363, none: 2266
pcds(compl): 14972 >65%, 7178 >= 40%, 2970 < 40% of 25120, 4986 < 50%
---
nok=25730  mapOK, noTE
onex=1952, onei=4645
pcds(compl): 12122 >65%, 6285 >= 40%, 2512 < 40% of 20919, 4245 < 50%
---

nok=2136  onlyTE
aug	816
rna	1098
ver1	222
ocl inpar	499
ocl nocl	876
ocl orlog	321
ocl upar	440
FPKM>=0.02 2023, FPKM>=1 1162
onex=746, onei=1794
gmap ok80: 1714, split: 106, part: 315, none: 1
pcds(compl): 612 >65%, 613 >= 40%, 346 < 40% of 1571, 553 < 50%
---

=cut
