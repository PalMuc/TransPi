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


my $MINW=   $ENV{minw}||40; # should get from read-length, eg $rlen * 0.75
my $MAXMIS= $ENV{maxmis}||10; # max total mismatches; get some garbage cases w/ scores of mis+snp
my $OPT2SNP=$ENV{snp2}||0;
my $debug=  $ENV{debug}||0;
my $MISUSE= $ENV{misuse}||0; # use read/genome mismatches not snpmap for tabulation

our %codon_table;
our %snpmap; # my $snpmap=$ENV{snpmap};
my %snpd=(); my %snpg= (); my %snpc=(); # GLOBAL=snpmap; NOT LOCAL for sumloc
my %snpa=(); # LOCAL per read with snps:

if($ENV{snpmap} and open(F,$ENV{snpmap})) { 
  my $ns=0;
  while(<F>){ 
    s/^>//; s/:/\t/; 
    my($sid,$gid,$gb,$ac)=split;  $ns++;
    my($rc,$sc)= split "",$ac;
    $snpmap{$gid}{$sid}="$gb\t$rc\t$sc"; #? both or gb
    # sid here should match gsnap outputs: {gp}s1b1234
    my($sg,$sp)= $sid =~ m/^(\w\w)(\S+)/; 
    $snpg{$gid}{$gb}{$sg}++; 
    $snpd{$gid}{$gb}{$sp}++; 
    $snpc{$gid}{$gb}{$sg}=$sc;
    ##$snpa{$gid}{$gb}{$sg}++ if($sngot{$sn});
  } close(F);
  warn "#snpmap read n=$ns\n" if $debug;
}

my($ro,$readseq,$rn,$rtyp,$rqual,$rid);
my %sumloc; ##      $sumloc{$gid}{$o}{$locdat} ++;
my %sumgot;

while(<>) {

  # add opt to count non-snp reads, but only at snp positions? (snpmap)
  # do as new prog; count all reads at genomic snp pos + count of snp changes.
  # FIXME: 1 read can have align to 2+ cds seqs, rn=#cds align; (alttr): Thecc1EG024602t1,Thecc1EG024602t2..
  # .. dont clear $rs == $readseq until blank line or new <> readseq
  # .. maybe limit to $rtyp =~ /concordant|../
  # .. unpaired sometimes are end-of-coding reads; keep if long enough 

  if(/^([<>])/) { $ro=$1; 
    ($readseq,$rn,$rtyp,$rqual,$rid)=split; $readseq=~s/^.//; 

  } elsif(/matches:\d/ and $readseq) { ## /snps:\w/ and 
    s/^(.)//; my $contc=$1; # FIX: ^(.) is blank or ',' for gap/intron continuation
# TaTTCAGCTAACTGCTGCagagtgttcctgc---------------------------------------------------------       1..18   -Thecc1EG016803t1:1091..1074    
#,------------------TCTCTTCTAGTGCTTCAGATTGAAATCTAGTCATAGATCGAAGACCCTCAAgTGCTG    19..75  -Thecc1EG016803t1:1060..1004    del:13..end:0,matches:56,sub:1+0=1
    my($cs,$cl,$cid,$ms)=split; # FIXME: ^(.) is blank or ',' for gap/intron continuation
    my($lb,$le)= $cl=~m/(\d+)..(\d+)/; my $lw=1+$le-$lb;
    my($go,$gid,$gb,$ge)= $cid =~ m/^(.)(\w+):(\d+)..(\d+)/;
    my $rev=($go eq "-")?1:0; 
    my($mm,$sm,$tm)=m/sub:(\d+).(\d+)=(\d+)/;
    my $rs= $readseq;
    
    my @sngot=(); my %sngot=();
    if( my($sngot)= m/snps:(\S+)/ ) { # this may be flaky, doesnt always agree with changes for gsn
      @sngot= map{ s/^\d+\@//; $_ } split/[\|]/,$sngot;
      map{ $sngot{$_}++ } @sngot;
    }
    
    # these are all snpmap GLOBAL, but for sngot = snpa hash
    my %snpa=(); # LOCAL for read w/ snps:
#    ## my %snpd=(); my %snpg= (); my %snpc=(); # NOT LOCAL for sumloc
#     my @gsn= sort keys %{$snpmap{$gid}}; # use all genomic snp now
#     # good only w/ snpmap
#     foreach my $sn (@gsn) { # was sn=sngot
#       my($sg,$sp)= $sn =~ m/^(\w\w)(\S+)/; 
#       my($sl,$rc,$sc)= split"\t",$snpmap{$gid}{$sn};
#       if($sl and $sp) { 
# #         $snpg{$gid}{$sl}{$sg}++; 
# #         $snpd{$gid}{$sl}{$sp}++; 
# #         $snpc{$gid}{$sl}{$sg}=$sc;
#         $snpa{$gid}{$sl}{$sg}++ if($sngot{$sn});
#       }
#     } 
    
    if($lw < length($cs)) {
      $cs=substr($cs,$lb-1,$lw);
      $rs=substr($rs,$lb-1,$lw);
      # $gb += $lb-1; $ge -= $lmax - $le; # after rev?  
    }
    
    my $skip= ($tm>$MAXMIS or $lw < $MINW) ? 1: 0; # change to $lw - $tm < $MINW ?
    #OFF?# $skip=1 if($mm>0); # for now, till can split mismat from snp
    if($skip) {  next; } # NOT: $rs="";
    
    if($rev) {
      ($gb,$ge)= ($ge,$gb);  ## gb,ge  #? need $lb,$le adjust?
      map{ $_= reverse $_; $_ =~ tr/acgtACGT/tgcaTGCA/; } ($cs,$rs); 
    }
    
    my @src= split"",$cs; 
    my @read=split"",$rs;
    # my @si=(); for my $i (0..$#src){ push @si, $i if($src[$i] =~ /[acgt]/); } # any mismatch to c-seq
    ## add Read-count at SNP locs without mismatch/snp change...
    my @ri=(); #  # my @rdsnpi=();
    my @misi=(); # do both ways? i.e. dont trust snpmap, check positions of all mismatches?
    for my $i (0..$#src) { # any read (no mismatch or +mis?)
      my $o= $i+$gb; 
      # my $imis=($src[$i] =~ /[acgt]/)?1:0;
      push @misi, $i if($src[$i] =~ /[acgt]/);
      if($snpg{$gid}{$o}) {
        push @ri, $i;
        # my $imis=($src[$i] =~ /[acgt]/)?1:0;
        # push @misi, $imis;  
        # my @sg= sort keys %{$snpc{$gid}{$o}};
        # my $irdsnp= ($snpc{$gid}{$o}{$sg} eq $read[$i]) ? 1 : 0; #  << test this snp base for read agreement?
        # push @rdsnpi, $irdsnp; #? need for each group; test if $read[$i] eq $snpbase{$gid}{$o}{xxx}
        }
    }  
    my @li= ($MISUSE) ? @misi : @ri;   # all positions with genomic-snp now
    #my @li= @si; # or only changed positions (snp+mis)
    
    # insert opt to pick reads w/ 2+snp from same group
    if($OPT2SNP) {
      my %sgn=(); my %sgl=();
      foreach my $j (0..$#li) {  
        my $i=$li[$j];  my $o= $i+$gb;  
        next unless($snpg{$gid}{$o});  # NOW or $snpa{$gid}{$o} == snpgot
        my @sg= sort keys %{$snpg{$gid}{$o}};
        map{ $sgn{$_}++; $sgl{$_}.="$i,"; } @sg;
      }
     my %l2=(); foreach my $sg (keys %sgn) { if($sgn{$sg}>1) { map{ $l2{$_}++ } split",", $sgl{$sg}; } }
     @li= sort{$a<=>$b} keys %l2;
    }
    
    foreach my $j (0..$#li) {
      my $i= $li[$j]; 
      # my $imis= $misi[$j]; # 1 if read changed from ref base
      my $readc= $read[$i];
      my $o= $i+$gb; 
      my $co=$o % 3; $co=3 if($co==0);
      my $scod= substr($cs,$i-($co-1),3);  
      my $rcod= substr($rs,$i-($co-1),3); #? should be same length, not
      #? next unless($snpd); # mismatch loc?
      next if(length($scod) < 3 or length($rcod) < 3); # endofseq: skip if scod/rcod < 3?

      # $snpa{$gid}{$o}{$sg}++ if($sngot{$sn});
 
      #.. need to sum:  $gid, $i/$o, $readc, $scod, $rcod
      #.. defer at sumout: co, snpd, snpg, snpa, saa, raa
      #...
if(1) {
      my $locdat= join("\t", $scod, $rcod, $readc);
      $sumloc{$gid}{$o}{$locdat} ++;

      # add snpgot list here?
      # my $sngot= join ",",sort keys %{$snpa{$gid}{$o}};
      # $sumgot{$gid}{$o}{$sngot} ++;
      
} else {      
      # %snpmap
      my $snpd= ($snpd{$gid}{$o}) ? join",", sort keys %{$snpd{$gid}{$o}} : 0; #snp-glocid
      #.. check read-base agreement w/ snp-group snp-base? # do ditto with snpa/snpgot
      my @snpg= ($snpg{$gid}{$o}) ? sort keys %{$snpg{$gid}{$o}} : (0); #snp-groupid
      my $snpg= join",", map{ ( $snpc{$gid}{$o}{$_} eq $readc ) ? uc($_) : $_; } @snpg;  
      
      #? drop this, local, others global = snpmap
      # my @snpa= ($snpa{$gid}{$o}) ? sort keys %{$snpa{$gid}{$o}} : (0); #snpgot-groupid
      my $snpa= "na"; # join",", map{ ( $snpc{$gid}{$o}{$_} eq $readc ) ? uc($_) : $_; } @snpa;  
      
      my $saa= $codon_table{uc($scod)};
      my $raa= $codon_table{uc($rcod)}; my $seq=($saa eq $raa)?"=":">";
      #! change this print to count reads, then print sum per gid,o,.. change position
      my $nread=1;
      print join("\t",$gid,$o,$nread,$co,$scod,$rcod,"$saa$seq$raa",$snpd,$snpa,$snpg)."\n";
# BEGIN{print join("\t",qw(GeneID__________ bCDS Nrd Cdn Gcdn Tcdn dAA SNPid SNPgot SNPgrp))."\n"; }
      # "$gb-$ge:$go", "$mm+$sm", # save up ; count all reads.
}

    }
     
    #NOT# $rs="";
   } 
}

sumout(); # for %sumloc

sub sumout {
  foreach my $gid (sort keys %sumloc) {
    foreach my $o (sort{$a <=> $b} keys %{$sumloc{$gid}}) {
      my $co=$o % 3; $co=3 if($co==0);
      my $snpd= ($snpd{$gid}{$o}) ? join",", sort keys %{$snpd{$gid}{$o}} : 0; #snp-glocid
      my @snpg= ($snpg{$gid}{$o}) ? sort keys %{$snpg{$gid}{$o}} : (0); #snp-groupid

      ## my $sngot= join ",",sort keys %{$snpa{$gid}{$o}};
#      my @sngot= sort keys %{ $sumgot{$gid}{$o} };
#      my @snpa= ($snpa{$gid}{$o}) ? sort keys %{$snpa{$gid}{$o}} : (0); #snpgot-groupid
      # ^^ %snpa is not global; can it append to locdat?
      
      my @locdat= sort keys %{ $sumloc{$gid}{$o} };
      foreach my $locdat (@locdat) {
        my($scod,$rcod,$readc)=split"\t",$locdat;
        my $nread= $sumloc{$gid}{$o}{$locdat};
        #.. check read-base agreement w/ snp-group snp-base? # do ditto with snpa/snpgot
        my $snpg= join",", map{ ( $snpc{$gid}{$o}{$_} eq $readc ) ? uc($_) : $_; } @snpg;  
        
        # my $snpa= "na"; #join",", map{ ( $snpc{$gid}{$o}{$_} eq $readc ) ? uc($_) : $_; } @snpa;  
        #.. maybe drop SNPgot as not informative? subset of true snp gotten?

        my $saa= $codon_table{uc($scod)};
        my $raa= $codon_table{uc($rcod)}; 
        my $seq= ($saa eq $raa)?"=":">";
        
        print join("\t",$gid,$o,$nread,$co,$scod,$rcod,"$saa$seq$raa",$snpd,$snpg)."\n"; #  $snpa,
BEGIN{  print join("\t",qw(GeneID__________ bCDS Nrd Cdn Gcdn Tcdn dAA SNPid SNPgrp))."\n"; } #  SNPgot
        
        }
    }
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

