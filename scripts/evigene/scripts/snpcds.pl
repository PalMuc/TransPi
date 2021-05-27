#!/usr/bin/env perl
# snpcds.pl

=item about snpcds

  tabulate snp changes in rnaseq reads aligned by gsnap to gene coding sequence 
  
=item output

  table of gene cds locations with snp changes, including codon and amino change, and snp ids with group tag
  in this format:
GeneID__________        bCDS    cdn     scod    qcod    dAA     SNPid   SNPgrp
Thecc1EG019914t1        558     3       GCc     GCA     A=A     s4b25709437     cc,po,uf
Thecc1EG019914t1        567     3       AGt     AGC     S=S     s4b25709446     cc,po,ts,uf
Thecc1EG029926t1        1698    3       TTg     TTA     L=L     s6b24968307     po,uf
  
=item input

  a. gsnap full output including "snps:.." ids, only these lines are processed

  b. table of snps as used for gsnap snp database: 
      >snpid  geneid:base  refbase/snpbase
    snptable is needed to distinguish which change is which snp
    where snpid by pref has 2-char group prefix, may be genomic location   
      >ccs1b1586      Thecc1EG000001t1:79      AG
      >pos1b1586      Thecc1EG000001t1:79      AG
      >ccs1b2094      Thecc1EG000001t1:359     CT

  ** update option: select only 2+ snp in same mate pair

=item gsnap output as input

 # aligned read  
  >TTTTTAGTGCCACTTTCCTTTATTAGTGCAGCAGCATCTTCTTCTGCAGTGCCTCCTATCTCACCAATTATAACG    1 concordant   qual   IRIS:4:CacaoPistilTruseq:5:1:11164:7278
 # codingseq with lower char for changed positions: mismatches + snp  
   TTTTTAGTGCCACTTTCCTTTATTAGTGCAGCAGCATCTTCTTCTGCAGTGCCcCCTATCTCACCAATTATAACG    1..75   
   -Thecc1EG007715t1:851..777      start:0..end:0,matches:75,sub:0+1=1,snps:54@tss2b8519699|54@ccs2b8519699|54@ufs2b8519699
          segs:1,align_score:0,mapq:40    pair_score:0,insert_length:110

  # sub:0+1=1 says 0 mismatch + 1 snp = 1 total changes
  # snps:54@tss2b8519699|54@ccs2b8519699|54@ufs2b8519699 are ids of matched snp positions, here 3 groups at same location
  
  <ACAAAATTCCTTGACGATCCTCAAACAGAAGGTATCGTTATAATTGGTGAGATAGGAGGCACTGCAGAAGAAGAT    1 concordant     qual    IRIS:4:CacaoPistilTruseq:5:1:11164:7278
   ACAAAATTCCTTGACGATCCTCAAACAGAAGGTATCGTTATAATTGGTGAGATAGGgGGCACTGCAGAAGAAGAT    1..75   
   +Thecc1EG007715t1:742..816      start:0..end:0,matches:75,sub:0+1=1,snps:19@tss2b8519699|19@ccs2b8519699|19@ufs2b8519699
          segs:1,align_score:0,mapq:40    pair_score:0,insert_length:110
  
  >CTTTGATAAGGCTCCTACATTCTTCCAATGTCTTCTTGGAAGTCCTTTGACAAAATTTGTTGTTATATCTCTGCT    1 concordant         IRIS:4:CacaoPistilTruseq:5:1:6334:9510
   CTTTGATAAaGCTCtTAaATTCTTCCAATGTCTTCTTGGAAGTCCTTTGACAAAATTTGTTGTaATATCTCTGCT    1..75   
   -Thecc1EG005664t1:307..233      start:0..end:0,matches:75,sub:0+4=4,snps:18@tss1b37698143|15@tss1b37698140|10@tss1b37698135|64@ccs1b37698189|18@ccs1b37698143|15@ccs1b37698140|10@ccs1b37698135|18@pos1b37698143|15@pos1b37698140|10@pos1b37698135      segs:1,align_score:0,mapq:40    pair_score:0,insert_length:137
   #  sub:0+4=4 = 4 snp changes, with ids listed (some multiple groups at same position)  
   
  <GTTGATTTCAATAGGGCCTATTCATCATGGTAACACAAATTTGGCTCGAATGGAGAGGCAAAAGCAGAGATATAA    1 concordant         IRIS:4:CacaoPistilTruseq:5:1:6334:9510
   GTTGATTTCAATAGGGCCTATTCATCATGGTgACACAAATTTGGCTCGAATGGAGAGGCAAAAGCAGAGATATtAC   1..75   
   +Thecc1EG005664t1:171..245      start:0..end:0,matches:75,sub:0+2=2,snps:2@ccs1b37698189|44@ufs1b37698231       segs:1,align_score:0,mapq:40    pair_score:0,insert_length:137
            
=cut


my $MINW=   $ENV{minw}||50; # should get from read-length, eg $rlen * 0.75
my $MAXMIS= $ENV{maxmis}||10; # max total mismatches; get some garbage cases w/ scores of mis+snp
my $OPT2SNP=$ENV{snp2}||0;
my $debug=  $ENV{debug}||0;
our %codon_table;
our %snpmap; # my $snpmap=$ENV{snpmap};

if($ENV{snpmap} and open(F,$ENV{snpmap})) { 
  my $ns=0;
  while(<F>){ 
    s/^>//; s/:/\t/; my($sid,$gid,$gb,$ac)=split;  $ns++;
    $snpmap{$gid}{$sid}="$gb\t$ac"; #? both or gb
    # sid here should match gsnap outputs: {gp}s1b1234
  } close(F);
  warn "#snpmap read n=$ns\n" if $debug;
}

while(<>){

  if(/^([<>])/) { $ro=$1; ($rs,$rn,$rtyp,$rqual,$rid)=split; $rs=~s/^.//; }
  elsif(/snps:\w/ and /matches:/ and $rs) { 
    ($cs,$cl,$cid,$ms)=split; 
    ($lb,$le)= $cl=~m/(\d+)..(\d+)/; $lw=1+$le-$lb;
    ($go,$gid,$gb,$ge)= $cid =~ m/^(.)(\w+):(\d+)..(\d+)/;
    $rev=($go eq "-")?1:0; 
    ($mm,$sm,$tm)=m/sub:(\d+).(\d+)=(\d+)/;

    ($sn)= m/snps:(\S+)/;
    @sn= map{ s/^\d+\@//; $_ } split/[\|]/,$sn;

    %snpd= %snpg= (); 
if(1) { # %snpmap
    # good only w/ snpmap
    map{ my($sg,$sp)=m/^(\w\w)(\S+)/; my($sl,$sc)= split"\t",$snpmap{$gid}{$_};
         if($sl and $sp) { $snpg{$gid}{$sl}{$sg}++; $snpd{$gid}{$sl}{$sp}++; }
       } @sn;
} else {
    #fixme: need map of snp-genomeloc to snp-cdsloc, see gsnap input table
    map{ s/^(\w\w)//; my $sg=$1; $snpg{$_}.="$sg," } @sn; # drop if snpmap
    @snpg= sort keys %snpg; # sort{$a<=>$b} keys %snpg; # s4b000 < not numeric; cut s.b
    foreach $p (@snpg) { %sg1=map{ $_,1 } split",",$snpg{$p}; $snpg{$p}=join",",sort keys %sg1; } ;
}
    
    if($lw < length($cs)) {
      $cs=substr($cs,$lb-1,$lw);
      $rs=substr($rs,$lb-1,$lw);
      # $gb += $lb-1; $ge -= $lmax - $le; # after rev?
    }
    $skip= ($tm>$MAXMIS or $lw < $MINW) ? 1: 0;
    #OFF?# $skip=1 if($mm>0); # for now, till can split mismat from snp
    if($skip) { $rs=""; next; }
    if($rev) {
      ($gb,$ge)= ($ge,$gb);  ## gb,ge  #? need $lb,$le adjust?
      map{ $_= reverse $_; $_ =~ tr/acgtACGT/tgcaTGCA/; } ($cs,$rs); 
    }
    
    @s= split"",$cs; @r=split"",$rs;
    @li=(); for my $i (0..$#s){ push @li, $i if($s[$i] =~ /[acgt]/); }

    # insert opt to pick reads w/ 2+snp from same group
    if($OPT2SNP) {
      my %sgn=(); my %sgl=();
      foreach my $j (0..$#li) {  
        my $i=$li[$j];  my $o= $i+$gb;  
        next unless($snpg{$gid}{$o});
        my @sg= sort keys %{$snpg{$gid}{$o}};
        map{ $sgn{$_}++; $sgl{$_}.="$i,"; } @sg;
      }
     my %l2=(); foreach my $sg (keys %sgn) { if($sgn{$sg}>1) { map{ $l2{$_}++ } split",", $sgl{$sg}; } }
     @li= sort{$a<=>$b} keys %l2;
    }
    
    foreach my $j (0..$#li) {  
      $i=$li[$j]; 
      $o= $i+$gb; $co=$o % 3; $co=3 if($co==0);
if(1) { # %snpmap
      $snpd= ($snpd{$gid}{$o}) ? join",", sort keys %{$snpd{$gid}{$o}} : 0; #snp-glocid
      $snpg= ($snpg{$gid}{$o}) ? join",", sort keys %{$snpg{$gid}{$o}} : 0; #snp-groupid
} else {
      $snpd=$snpg[$j]; $snpg=$snpg{$snpd}; # drop; not always, mismatches.. need snp-loc>cdsloc convert
}      
      $scod= substr($cs,$i-($co-1),3);  
      $rcod= substr($rs,$i-($co-1),3); #? should be same length, not
      #? next unless($snpd); # mismatch loc?
      next if(length($scod) < 3 or length($rcod) < 3); # endofseq: skip if scod/rcod < 3?
      my $saa= $codon_table{uc($scod)};
      my $raa= $codon_table{uc($rcod)}; $seq=($saa eq $raa)?"=":">";
      print join("\t",$gid,$o,$co,$scod,$rcod,"$saa$seq$raa",$snpd,$snpg)."\n";
BEGIN{print join("\t",qw(GeneID__________ bCDS cdn scod qcod dAA SNPid SNPgrp))."\n"; } # EG029926t1
      # "$gb-$ge:$go", "$mm+$sm", # save up ; count all reads.
    }
     
    $rs="";
   } 
}


BEGIN{
 # stops: TAG  TGA  TAA
  %codon_table = (    
  TTT => 'F', TTC => 'F', TTA => 'L', TTG => 'L',
  CTT => 'L', CTC => 'L', CTA => 'L', CTG => 'L',
  ATT => 'I', ATC => 'I', ATA => 'I', ATG => 'M',
  GTT => 'V', GTC => 'V', GTA => 'V', GTG => 'V',
  TCT => 'S', TCC => 'S', TCA => 'S', TCG => 'S',
  CCT => 'P', CCC => 'P', CCA => 'P', CCG => 'P',
  ACT => 'T', ACC => 'T', ACA => 'T', ACG => 'T',
  GCT => 'A', GCC => 'A', GCA => 'A', GCG => 'A',
  TAT => 'Y', TAC => 'Y', TAA => '*', TAG => '*',
  CAT => 'H', CAC => 'H', CAA => 'Q', CAG => 'Q',
  AAT => 'N', AAC => 'N', AAA => 'K', AAG => 'K',
  GAT => 'D', GAC => 'D', GAA => 'E', GAG => 'E',
  TGT => 'C', TGC => 'C', TGA => '*', TGG => 'W',
  CGT => 'R', CGC => 'R', CGA => 'R', CGG => 'R',
  AGT => 'S', AGC => 'S', AGA => 'R', AGG => 'R',
  GGT => 'G', GGC => 'G', GGA => 'G', GGG => 'G' 
  );
}

