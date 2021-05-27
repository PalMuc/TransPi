#!/usr/bin/env perl
# cdsalnpatt.pl 

=item about

  evigene/scripts/prot/cdsalnpatt.pl
  cdsalnpatt.pl : measure cds align patterns, stats to compare/discriminate
    patterns of paralogs and alt-transcripts.
    simple case: pairwise clustal/muscle -clw output align
    -- per align row pair, set base array of same|diff|gap|N
    -- base 1 = cds start ATG assumed, or check for ATG
    -- summarize patterns from total and 1/3, 2/3, 3/3 partitions of longest? cds
        count ndiff, ngap, ignore N or count?
        output per pairalign:  Tds/dn/g/n,Bds/dn/g/n,Mds/dn/g/n,Eds/dn/g/n
        for T=total, B=begin5, M=middle, E=end3, 
        ds=diff-synon(codon3), dn=diff-nonsyn(codon1,2),g=gaps,n=basetotal
           ^ should these be percents or counts?
        
=cut

use strict;
use Getopt::Long;

my $debug=$ENV{debug}||0;
my $NOPCT=$ENV{nopct}||0;
my $NOPARTS=$ENV{noparts}||0;
my ($id1,$id2,$nullop)= (0) x 10;

my $optok= GetOptions(
#   "cds=s", \$incds, "c2aa=s", \$incdsaa, # add output default: inaa.diffa
#   "stopcheck!", \$STOPCHECK,   "gapoffbyone!", \$GAPOFFBYONE,  "ok!", \$SHOWOK,
  "id1=s", \$id1,"id2=s", \$id2,
  "nopercent|nopct!", \$NOPCT, "pct", \$nullop, #Old
  "noparts!", \$NOPARTS, ## is nonoparts valid off?
  "debug!", \$debug,
  );
  
die "usage: cdsalnpatt.pl cdspairalign_genes*.clw > cdsalign.stats
  cat gene[12].cds | muscle -quiet -clw | cdsalnpatt.pl > cds12align.stats
  opts: -[no]debug;  -percent
  \n"  unless($optok);

  # base aln codes: 1=same, -1=diff base, 0=N1 or 2, -2=gap1, -3=gap2
use constant { bEQ => 1, bNE => -1, bNN => 0, bGAP1 => -2, bGAP2 => -3, };

readalign();

sub readalign
{
  my($lid,$lal,$nt,$len1,$len2)= (0) x 9;
  my(@c);
  while(<>) {
    next if(/^(MUSCLE|CLUST|\s|\W)/);
    my($id,$al)=split; next unless($al);
    unless($id1) { $id1=$id; }
    if($id eq $id1) { $lid=$id; $lal=$al; }
    elsif($lid eq $id1) { 
      unless($id2) { $id2=$id; } elsif($id ne $id2) { next; }
      if($nt==0) {
        my $mat=0; my $m1=index($lal,'ATG'); my $m2=index($al, 'ATG');
        if($m1>=0) { $mat=$m1; } if($m2>=0 and $m2<$mat) { $mat=$m2; }
        if($mat>0) { $lal=substr($lal,$mat); $al=substr($al,$mat); }  
      }
      
      my @al=split "",$lal; my @bl= split "",$al;   
      for my $i (0..$#al) { 
        my($ab,$bb)=($al[$i],$bl[$i]);
        my $c=($ab eq "-")?bGAP1: ($bb eq "-")?bGAP2: ($ab eq $bb)?bEQ: ($ab eq "N" or $bb eq "N")?bNN: bNE;
        $c[$nt]= $c; $nt++;
        $len1++ unless($c==bGAP1); $len2++ unless($c==bGAP2);
      }
      $lid=$lal=$al="";
    }
    ## $lid=$id; ## $lal=$al;
  }
  return 0 unless($nt>1);
  
  # output:  Tds/dn/g/n,Bds/dn/g/n,Mds/dn/g/n,Eds/dn/g/n
  # add max eq runs, foreach of T,B,M,E?
  my $b1=1+int($nt/3); my $m2=$b1+$b1; my $e3=$nt;
  my (%sum,%maxeq); my $eqr=0;
  for my $i (0..$nt-1) {
    my $c=$c[$i];
    my $isyn=($i % 3 == 2)?1:0;
    my $p= ($i<$b1)?'B': ($i<$m2)? 'M': 'E';
    my $s='.';
    if($c==bGAP1 or $c==bGAP2) { $s='gp'; } 
    elsif($c==bEQ) { $eqr++; if($isyn) { $s='es';} else { $s='en'; } }
    elsif($c==bNE) { if($isyn) { $s='ds'; } else { $s='dn'; } }
    # elsif($c==bNN) { next; } else { next; } 
    if($s ne '.') { for my $t ('T',$p) { $sum{$t}{$s}++; $sum{$t}{'tn'}++; } }
    if($eqr and (($c != bEQ and $s ne '.') or $i == $nt-1)) { $sum{'T'}{'er'}= $eqr if($eqr > $sum{'T'}{'er'}); $eqr=0; }
  }
  
  my @pt=($NOPCT && $NOPARTS)?qw(N):($NOPARTS)?qw(N T):qw(N T B M E); # option to drop B,M,E ?
  my @st=qw(ds dn gp er tn);
  my $out="$id1\t$id2\t";
  for my $pp (@pt) {
    $out.="$pp:"; my $p=($pp eq 'N')?'T':$pp;
    my $tn= $sum{$p}{'tn'}; my $tng= $tn - $sum{$p}{'gp'};
    for my $s (@st) {
      my $c= $sum{$p}{$s}||0; my $cp=0; my $n=$tn;
      if($s eq "ds") { $n=$c+$sum{$p}{'es'}; }  # ds/ds+es == Ks
      elsif($s eq "dn") { $n=$c+$sum{$p}{'en'}; } # dn/dn+en == Ka
      elsif($s eq "er") {  next unless($p eq "T"); 
         $n=$tng; $n=$len1 if($n>$len1); $n=$len2 if($n>$len2);  #  FIXME? %er min of len1,len2; prob not needed
         } 
      ## maybe add es,en counts to N output for later recalc Ka,Ks;
      if($pp eq 'N' or $NOPCT) { $out.= $c.$s.","; }
      else { my $cp=($n>0)? int(0.5+100*$c/$n) : 0; $out.= $cp.$s.","; }
    }
    $out=~s/\W$//; $out.="\t";
  }
  $out.="len:$len1,$len2"; #chomp($out);
  print $out,"\n"; # one line per pairalign ?
}


__END__

=item input eg

  MUSCLE (3.8) multiple sequence alignment
  
  Dapma5xEVm001968t1      ATGACAACTCCCACAACAAATCGGATTCACGTAACCTATAAAGTAGCAGACGAAATTATT
  Dapma5xEVm001953t1      ATG---------------GATCGAATC---------------------------------
                          ***                **** **                                  
  Dapma5xEVm001968t1      ACCAGGAGTGATCAACGTGCACGTCGTTCCTCATACTCACGATGATGTCGGATGGTTGAA
  Dapma5xEVm001953t1      ACCGGGCTTCATCAACGTCCACCTGGTTCCCCACACTCACGATGATGTCGGCTGGCTCAA
                          *** **  * ******** *** * ***** ** ***************** *** * **
=item output eg

3534 Jan 28 23:53 arp475.Dapma5xEVm001953t1.cds 20alts
3503 Jan 28 23:53 arp475.Dapma5xEVm001963t1.cds 26alts
3509 Jan 28 23:53 arp475.Dapma5xEVm001968t1.cds 17alts
3508 Jan 28 23:53 arp475.Dapma5xEVm001993t1.cds  6alts

dpx, genome models, no alts
1754 Jan 29 13:09 arp475.hxAUG25s63g47t1.cds
3049 Jan 29 13:09 arp475.hxAUG25s63g57t1.cds
1731 Jan 29 13:09 arp475.hxAUG26res63g33t1.cds
3236 Jan 29 13:09 arp475.hxAUG26up1s5g243t1.cds
3215 Jan 29 13:09 arp475.hxNCBI_GNO_118484.cds
3028 Jan 29 13:09 arp475.hxNCBI_GNO_504144.cds
2381 Jan 29 13:09 arp475.hxNCBI_GNO_65294.cds

cat arp475.Dapma5xEVm0019{68,53}t1.cds | $gs/muscle/muscle -quiet -clw | $evigene/scripts/prot/cdsalnpatt.pl -pct
-nopct
Dapma5xEVm001968t1	Dapma5xEVm001953t1	N:349ds,779dn,483gp,3414tn	T:349ds,779dn,483gp,3414tn	B:91ds,254dn,138gp,1139tn	M:160ds,274dn,97gp,1139tn	E:98ds,251dn,248gp,1136tn	len:3168,3177
-pct
Dapma5xEVm001968t1	Dapma5xEVm001953t1	N:349ds,779dn,483gp,3414tn	T:36ds,40dn,14gp,100tn	B:27ds,38dn,12gp,100tn	M:46ds,40dn,9gp,100tn	E:33ds,42dn,22gp,100tn	len:3168,3177

# dplx paralogs:
pt=arp475.hxAUG26up1s5g243t1.cds
pt=arp475.hxAUG25s63g57t1.cds
for ix in arp475.hx*.cds; do { if [ $ix = $pt ]; then continue; fi; 
cat $pt $ix | $gs/muscle/muscle -quiet -clw | $evigene/scripts/prot/cdsalnpatt.pl; 
} done

hxAUG26up1s5g243t1	hxAUG25s63g47t1 	N:211ds,332dn,1557gp,3111tn	T:41ds,32dn,50gp,100tn	B:38ds,35dn,5gp,100tn	M:47ds,29dn,51gp,100tn	E:26ds,8dn,95gp,100tn	
hxAUG26up1s5g243t1	hxAUG25s63g57t1  	N:348ds,770dn,292gp,3128tn	T:37ds,41dn,9gp,100tn	  B:42ds,32dn,10gp,100tn	M:32ds,46dn,12gp,100tn	E:36ds,44dn,6gp,100tn	
hxAUG26up1s5g243t1	hxAUG26res63g33t1	N:184ds,418dn,1499gp,3085tn	T:35ds,40dn,49gp,100tn	B:30ds,33dn,70gp,100tn	M:33ds,34dn,70gp,100tn	E:37ds,44dn,6gp,100tn	
hxAUG26up1s5g243t1	hxNCBI_GNO_118484	N:335ds,748dn,278gp,3187tn	T:34ds,39dn,9gp,100tn	  B:26ds,37dn,6gp,100tn	M:40ds,41dn,11gp,100tn	E:38ds,39dn,9gp,100tn	
hxAUG26up1s5g243t1	hxNCBI_GNO_504144	N:378ds,710dn,430gp,3185tn	T:41ds,39dn,14gp,100tn	B:49ds,32dn,8gp,100tn	M:38ds,40dn,16gp,100tn	E:37ds,44dn,17gp,100tn	
hxAUG26up1s5g243t1	hxNCBI_GNO_65294	N:338ds,527dn,917gp,3115tn	T:46ds,36dn,29gp,100tn	B:38ds,38dn,13gp,100tn	M:53ds,36dn,12gp,100tn	E:48ds,31dn,64gp,100tn	

'daphplx_(hxAUG26up1s5g243t1|hxAUG25s63g47t1|hxAUG25s63g57t1|hxAUG26res63g33t1|hxNCBI_GNO_118484|hxNCBI_GNO_504144|hxNCBI_GNO_65294)\t'
dplx inparalogs:
daphplx_hxNCBI_GNO_65294	daphplx_hxAUG25s63g57t1	ARP6d475	0.0	98	inpar1	80.9
daphplx_hxAUG26res63g33t1	daphplx_hxAUG25s63g57t1	ARP6d475	0.0	99	inpar1	71
daphplx_hxAUG25s63g47t1	daphplx_hxAUG25s63g57t1	ARP6d475	0.0	97	inpar1	65.9

inpar:  note low ds/dn
hxAUG25s63g57t1	hxAUG25s63g47t1    	N:5ds,8dn,1505gp,3010tn 	T:1ds,1dn,50gp,100tn	B:2ds,1dn,10gp,100tn	M:0ds,0dn,40gp,100tn	E:0ds,40dn,99gp,100tn	
hxAUG25s63g57t1	hxAUG26res63g33t1 	N:1ds,1dn,1293gp,2907tn 	T:0ds,0dn,44gp,100tn	B:0ds,0dn,99gp,100tn	M:0ds,0dn,34gp,100tn	E:0ds,0dn,0gp,100tn	
hxAUG25s63g57t1	hxNCBI_GNO_65294  	N:23ds,60dn,695gp,2929tn	T:3ds,4dn,24gp,100tn	B:7ds,10dn,8gp,100tn	M:1ds,0dn,0gp,100tn 	E:0ds,0dn,63gp,100tn	
oldpar:
hxAUG25s63g57t1	hxAUG26up1s5g243t1	N:348ds,770dn,292gp,3128tn	T:37ds,41dn,9gp,100tn	B:42ds,32dn,10gp,100tn	M:32ds,46dn,12gp,100tn	E:36ds,44dn,6gp,100tn	
hxAUG25s63g57t1	hxNCBI_GNO_118484 	N:461ds,485dn,252gp,3099tn	T:49ds,26dn,8gp,100tn	B:51ds,15dn,10gp,100tn	M:39ds,32dn,10gp,100tn	E:56ds,29dn,4gp,100tn	
hxAUG25s63g57t1	hxNCBI_GNO_504144 	N:435ds,657dn,362gp,3076tn	T:48ds,36dn,12gp,100tn	B:54ds,32dn,7gp,100tn	M:42ds,39dn,13gp,100tn	E:48ds,38dn,16gp,100tn	

# dmag paralogs: all oldpar
pt=Dapma5xEVm0019; j=68;
for i in 53 63 93; do {  cat arp475.${pt}${j}t1.cds dmcdsf/${pt}${i}t1.cds | $gs/muscle/muscle -quiet -clw | $evigene/scripts/prot/cdsalnpatt.pl; } done

Dapma5xEVm001968t1	Dapma5xEVm001953t1	N:349ds,779dn,483gp,3414tn	T:36ds,40dn,14gp,100tn	B:27ds,38dn,12gp,100tn	M:46ds,40dn,9gp,100tn	E:33ds,42dn,22gp,100tn	
Dapma5xEVm001968t1	Dapma5xEVm001963t1	N:418ds,702dn,393gp,3366tn	T:42ds,35dn,12gp,100tn	B:51ds,30dn,11gp,100tn	M:38ds,38dn,12gp,100tn	E:38ds,38dn,11gp,100tn	
Dapma5xEVm001968t1	Dapma5xEVm001993t1	N:377ds,732dn,529gp,3422tn	T:39ds,38dn,15gp,100tn	B:37ds,35dn,17gp,100tn	M:45ds,38dn,12gp,100tn	E:36ds,41dn,17gp,100tn	


alts:
mkdir dmcdsf;
cat arp475.dmagalts.cds | perl -ne'if(/^>(\S+)/) { $d=$1; close(F) if($fn); $fn="dmcdsf/$d.cds"; open(F,">$fn"); } print F $_;'

for i in 2 3 4 5 6 7 8 9; do {  cat arp475.${pt}t1.cds dmcdsf/${pt}t$i.cds | $gs/muscle/muscle -quiet -clw | $evigene/scripts/prot/cdsalnpatt.pl; } done

pt=Dapma5xEVm001953; 
Dapma5xEVm001953t1	Dapma5xEVm001953t2	N:0ds,8dn,191gp,3205tn  	T:0ds,0dn,6gp,100tn 	B:0ds,0dn,2gp,100tn   	M:0ds,0dn,0gp,100tn 	E:0ds,1dn,16gp,100tn	
Dapma5xEVm001953t1	Dapma5xEVm001953t3	N:6ds,11dn,1530gp,3177tn	T:1ds,1dn,48gp,100tn	B:32ds,22dn,94gp,100tn	M:0ds,0dn,36gp,100tn	E:0ds,0dn,14gp,100tn	
Dapma5xEVm001953t1	Dapma5xEVm001953t4	N:4ds,6dn,1128gp,3177tn 	T:1ds,0dn,36gp,100tn	B:14ds,11dn,92gp,100tn	M:0ds,0dn,0gp,100tn 	E:0ds,0dn,14gp,100tn	
Dapma5xEVm001953t1	Dapma5xEVm001953t5	N:0ds,1dn,539gp,3333tn  	T:0ds,0dn,16gp,100tn	B:0ds,0dn,35gp,100tn  	M:0ds,0dn,0gp,100tn 	E:0ds,0dn,14gp,100tn	
Dapma5xEVm001953t1	Dapma5xEVm001953t6	N:8ds,18dn,1257gp,3177tn	T:1ds,1dn,40gp,100tn	B:23ds,25dn,90gp,100tn	M:0ds,0dn,14gp,100tn	E:0ds,0dn,14gp,100tn	
Dapma5xEVm001953t1	Dapma5xEVm001953t7	N:5ds,9dn,474gp,3177tn  	T:1ds,1dn,15gp,100tn	B:0ds,0dn,11gp,100tn  	M:0ds,0dn,0gp,100tn 	E:2ds,2dn,33gp,100tn	
Dapma5xEVm001953t1	Dapma5xEVm001953t8	N:0ds,0dn,489gp,3369tn  	T:0ds,0dn,15gp,100tn	B:0ds,0dn,2gp,100tn   	M:0ds,0dn,0gp,100tn 	E:0ds,0dn,42gp,100tn	
Dapma5xEVm001953t1	Dapma5xEVm001953t9	N:0ds,0dn,663gp,3324tn  	T:0ds,0dn,20gp,100tn	B:0ds,0dn,46gp,100tn  	M:0ds,0dn,0gp,100tn 	E:0ds,0dn,14gp,100tn	

pt=Dapma5xEVm001963;
Dapma5xEVm001963t1	Dapma5xEVm001963t2	N:10ds,9dn,1315gp,3267tn	T:2ds,1dn,40gp,100tn	B:1ds,0dn,56gp,100tn	M:1ds,0dn,7gp,100tn	E:3ds,3dn,58gp,100tn	
Dapma5xEVm001963t1	Dapma5xEVm001963t3	N:1ds,2dn,258gp,3276tn	  T:0ds,0dn,8gp,100tn	  B:0ds,0dn,24gp,100tn	M:0ds,0dn,0gp,100tn	E:0ds,0dn,0gp,100tn	
Dapma5xEVm001963t1	Dapma5xEVm001963t4	N:2ds,1dn,258gp,3276tn	  T:0ds,0dn,8gp,100tn	  B:0ds,0dn,24gp,100tn	M:0ds,0dn,0gp,100tn	E:0ds,0dn,0gp,100tn	
Dapma5xEVm001963t1	Dapma5xEVm001963t5	N:4ds,6dn,1506gp,3171tn	  T:1ds,1dn,47gp,100tn	B:18ds,15dn,95gp,100tn	M:0ds,0dn,27gp,100tn	E:0ds,0dn,20gp,100tn	
Dapma5xEVm001963t1	Dapma5xEVm001963t6	N:30ds,61dn,401gp,3187tn	T:3ds,3dn,13gp,100tn	B:9ds,9dn,10gp,100tn	M:0ds,0dn,0gp,100tn	E:0ds,0dn,28gp,100tn	
Dapma5xEVm001963t1	Dapma5xEVm001963t7	N:2ds,1dn,1162gp,3295tn	  T:0ds,0dn,35gp,100tn	B:0ds,0dn,45gp,100tn	M:0ds,0dn,0gp,100tn	E:1ds,0dn,60gp,100tn	
Dapma5xEVm001963t1	Dapma5xEVm001963t8	N:15ds,13dn,537gp,3171tn	T:2ds,1dn,17gp,100tn	B:1ds,0dn,18gp,100tn	M:0ds,0dn,0gp,100tn	E:9ds,5dn,33gp,100tn	
Dapma5xEVm001963t1	Dapma5xEVm001963t9	N:16ds,25dn,1444gp,3173tn	T:3ds,2dn,46gp,100tn	B:0ds,0dn,100gp,100tn	M:2ds,0dn,30gp,100tn	E:4ds,4dn,6gp,100tn	

pt=Dapma5xEVm001993;
Dapma5xEVm001993t1	Dapma5xEVm001993t2	N:0ds,1dn,1163gp,3152tn 	T:0ds,0dn,37gp,100tn	B:0ds,0dn,100gp,100tn	M:0ds,0dn,11gp,100tn	E:0ds,0dn,0gp,100tn	
Dapma5xEVm001993t1	Dapma5xEVm001993t3	N:12ds,25dn,657gp,3147tn	T:1ds,2dn,21gp,100tn	B:9ds,10dn,63gp,100tn	M:0ds,0dn,0gp,100tn	E:0ds,0dn,0gp,100tn	
Dapma5xEVm001993t1	Dapma5xEVm001993t4	N:0ds,0dn,1560gp,3270tn 	T:0ds,0dn,48gp,100tn	B:0ds,0dn,51gp,100tn	M:0ds,0dn,0gp,100tn	E:0ds,0dn,92gp,100tn	
Dapma5xEVm001993t1	Dapma5xEVm001993t5	N:2ds,3dn,69gp,3147tn   	T:0ds,0dn,2gp,100tn 	B:1ds,1dn,7gp,100tn 	M:0ds,0dn,0gp,100tn	E:0ds,0dn,0gp,100tn	
Dapma5xEVm001993t1	Dapma5xEVm001993t6	N:3ds,11dn,305gp,3180tn 	T:0ds,1dn,10gp,100tn	B:1ds,3dn,29gp,100tn	M:0ds,0dn,0gp,100tn	E:0ds,0dn,0gp,100tn	

pt=Dapma5xEVm001968;
Dapma5xEVm001968t1	Dapma5xEVm001968t2	N:4ds,14dn,420gp,3168tn 	T:0ds,1dn,13gp,100tn	B:0ds,0dn,11gp,100tn	M:0ds,0dn,0gp,100tn	E:2ds,2dn,29gp,100tn	
Dapma5xEVm001968t1	Dapma5xEVm001968t3	N:1ds,2dn,801gp,3435tn  	T:0ds,0dn,23gp,100tn	B:0ds,0dn,9gp,100tn 	M:0ds,0dn,7gp,100tn	E:0ds,0dn,53gp,100tn	
Dapma5xEVm001968t1	Dapma5xEVm001968t4	N:11ds,17dn,1044gp,3168tn	T:2ds,1dn,33gp,100tn	B:0ds,0dn,11gp,100tn	M:0ds,0dn,0gp,100tn	E:27ds,17dn,88gp,100tn	
Dapma5xEVm001968t1	Dapma5xEVm001968t5	N:12ds,17dn,114gp,3168tn	T:1ds,1dn,4gp,100tn 	B:4ds,3dn,11gp,100tn	M:0ds,0dn,0gp,100tn	E:0ds,0dn,0gp,100tn	
Dapma5xEVm001968t1	Dapma5xEVm001968t6	N:9ds,14dn,654gp,3168tn 	T:1ds,1dn,21gp,100tn	B:0ds,0dn,11gp,100tn	M:0ds,0dn,0gp,100tn	E:5ds,3dn,51gp,100tn	
Dapma5xEVm001968t1	Dapma5xEVm001968t7	N:11ds,14dn,933gp,3168tn	T:1ds,1dn,29gp,100tn	B:0ds,0dn,0gp,100tn 	M:0ds,0dn,0gp,100tn	E:27ds,17dn,88gp,100tn	
Dapma5xEVm001968t1	Dapma5xEVm001968t8	N:11ds,17dn,594gp,3168tn	T:1ds,1dn,19gp,100tn	B:0ds,0dn,11gp,100tn	M:0ds,0dn,0gp,100tn	E:6ds,4dn,46gp,100tn	
Dapma5xEVm001968t1	Dapma5xEVm001968t9	N:2ds,11dn,1053gp,3273tn	T:0ds,1dn,32gp,100tn	B:0ds,1dn,17gp,100tn	M:0ds,0dn,0gp,100tn	E:3ds,6dn,79gp,100tn	


=item hemoglobins

aabugs5/omcl/arp6bor2daph/heglobinf/

2470  heglobin.dmagmain.Dapma5xEVm004206t1.cds
2341  heglobin.dmagmain.Dapma5xEVm004686t1.cds
1593  heglobin.dmagmain.Dapma5xEVm009245t1.cds
1330  heglobin.dmagmain.Dapma5xEVm011559t1.cds
1174  heglobin.dmagmain.Dapma5xEVm013514t1.cds
1153  heglobin.dmagmain.Dapma5xEVm013664t1.cds
1115  heglobin.dmagmain.Dapma5xEVm014094t1.cds

pt=heglobin.dmagmain.Dapma5xEVm0
j=04206;

for i in 04686 09245 11559 13514 13664 14094; do {  
cat ${pt}${j}t1.cds ${pt}${i}t1.cds | $gs/muscle/muscle -quiet -clw | $evigene/scripts/prot/cdsalnpatt.pl -pct; } done

Dapma5xEVm004206t1	Dapma5xEVm004686t1	N:31ds,41dn,126gp,1767tn  	T:6ds,4dn,7gp,100tn   	B:3ds,2dn,0gp,100tn   	M:20ds,15dn,32gp,100tn	E:3ds,1dn,0gp,100tn	len:2148,2022
Dapma5xEVm004206t1	Dapma5xEVm009245t1	N:104ds,74dn,1350gp,2391tn	T:30ds,11dn,56gp,100tn	B:0ds,17dn,98gp,100tn 	M:28ds,9dn,55gp,100tn 	E:32ds,11dn,17gp,100tn	len:2148,1284
Dapma5xEVm004206t1	Dapma5xEVm011559t1	N:122ds,96dn,1106gp,2149tn	T:35ds,14dn,51gp,100tn	B:0ds,0dn,100gp,100tn 	M:36ds,12dn,54gp,100tn	E:35ds,14dn,0gp,100tn	len:2148,1044
Dapma5xEVm004206t1	Dapma5xEVm013514t1	N:93ds,185dn,1422gp,2223tn	T:35ds,35dn,64gp,100tn	B:28ds,33dn,71gp,100tn	M:40ds,35dn,64gp,100tn	E:35ds,35dn,57gp,100tn	len:2148,876
Dapma5xEVm004206t1	Dapma5xEVm013664t1	N:94ds,213dn,1288gp,2150tn	T:33ds,37dn,60gp,100tn	B:32ds,36dn,58gp,100tn	M:35ds,38dn,54gp,100tn	E:34ds,35dn,67gp,100tn	len:2148,864
Dapma5xEVm004206t1	Dapma5xEVm014094t1	N:103ds,181dn,1356gp,2062tn	T:44ds,39dn,66gp,100tn	B:32ds,37dn,78gp,100tn	M:54ds,42dn,48gp,100tn	E:37ds,35dn,68gp,100tn	len:2148,828

1215  heglobin.plx.hxAUG25p2s4g109t1.cds
1163  heglobin.plx.hxAUG25p2s4g110t1.cds
1141  heglobin.plx.hxAUG25p2s4g111t1.cds
1144  heglobin.plx.hxAUG25p2s4g113t1.cds
1144  heglobin.plx.hxAUG25p2s4g115t1.cds
 777  heglobin.plx.hxAUG25s13g137t1.cds
 718  heglobin.plx.hxAUG26res13g1t1.cds
 692  heglobin.plx.hxJGI_V11_301607.cds
1146  heglobin.plx.hxJGI_V11_92880.cds
1140  heglobin.plx.hxJGI_V11_93831.cds
1142  heglobin.plx.hxNCBI_GNO_1448043.cds
 858  heglobin.plx.hxNCBI_GNO_264133.cds
1276  heglobin.plx.hxNCBI_GNO_920044.cds

pt=heglobin.plx.
j=hxAUG25p2s4g109t1;

for i in hxAUG25p2s4g110t1 hxAUG25p2s4g111t1 hxAUG25p2s4g113t1 hxAUG25p2s4g115t1 hxAUG25s13g137t1 hxAUG26res13g1t1 hxJGI_V11_301607 hxJGI_V11_92880 hxJGI_V11_93831 hxNCBI_GNO_1448043 hxNCBI_GNO_264133 hxNCBI_GNO_920044; 
do { cat ${pt}${j}.cds ${pt}${i}.cds | $gs/muscle/muscle -quiet -clw | $evigene/scripts/prot/cdsalnpatt.pl -pct; } done

hxAUG25p2s4g109t1	hxAUG25p2s4g110t1	N:78ds,55dn,69gp,1116tn 	T:22ds,8dn,6gp,100tn  	B:20ds,9dn,18gp,100tn 	M:22ds,6dn,0gp,100tn  	E:25ds,9dn,0gp,100tn	len:1107,1056
hxAUG25p2s4g109t1	hxAUG25p2s4g111t1	N:95ds,98dn,72gp,1107tn 	T:28ds,14dn,7gp,100tn 	B:27ds,13dn,19gp,100tn	M:28ds,14dn,0gp,100tn 	E:28ds,15dn,1gp,100tn	len:1107,1035
hxAUG25p2s4g109t1	hxAUG25p2s4g113t1	N:107ds,104dn,69gp,1107tn	T:31ds,15dn,6gp,100tn 	B:32ds,13dn,18gp,100tn	M:30ds,15dn,0gp,100tn 	E:31ds,17dn,1gp,100tn	len:1107,1038
hxAUG25p2s4g109t1	hxAUG25p2s4g115t1	N:107ds,97dn,69gp,1107tn	T:31ds,14dn,6gp,100tn 	B:32ds,13dn,18gp,100tn	M:28ds,13dn,0gp,100tn 	E:34ds,17dn,1gp,100tn	len:1107,1038
hxAUG25p2s4g109t1	hxAUG25s13g137t1	N:71ds,163dn,564gp,1176tn	T:34ds,40dn,48gp,100tn	B:39ds,48dn,20gp,100tn	M:26ds,35dn,68gp,100tn	E:33ds,29dn,56gp,100tn	len:1107,681
hxAUG25p2s4g109t1	hxAUG26res13g1t1	N:69ds,155dn,548gp,1141tn	T:36ds,39dn,48gp,100tn	B:39ds,45dn,14gp,100tn	M:36ds,32dn,59gp,100tn	E:26ds,29dn,72gp,100tn	len:1107,627
hxAUG25p2s4g109t1	hxJGI_V11_301607	N:75ds,149dn,586gp,1142tn	T:41ds,40dn,51gp,100tn	B:54ds,41dn,19gp,100tn	M:28ds,40dn,46gp,100tn	E:0ds,31dn,89gp,100tn	len:1107,591
hxAUG25p2s4g109t1	hxJGI_V11_92880 	N:123ds,103dn,69gp,1107tn	T:36ds,15dn,6gp,100tn 	B:41ds,15dn,18gp,100tn	M:29ds,11dn,0gp,100tn 	E:38ds,19dn,1gp,100tn	len:1107,1038
hxAUG25p2s4g109t1	hxJGI_V11_93831 	N:122ds,117dn,78gp,1110tn	T:35ds,17dn,7gp,100tn 	B:37ds,19dn,19gp,100tn	M:33ds,12dn,0gp,100tn 	E:36ds,20dn,2gp,100tn	len:1107,1035
hxAUG25p2s4g109t1	hxNCBI_GNO_1448043	N:97ds,123dn,86gp,1114tn	T:28ds,18dn,8gp,100tn	B:37ds,20dn,20gp,100tn	M:24ds,10dn,0gp,100tn 	E:25ds,24dn,3gp,100tn	len:1107,1035
hxAUG25p2s4g109t1	hxNCBI_GNO_264133	N:105ds,199dn,402gp,1134tn	T:43ds,41dn,35gp,100tn	B:40ds,39dn,48gp,100tn	M:36ds,40dn,42gp,100tn	E:48ds,43dn,16gp,100tn	len:1107,759
hxAUG25p2s4g109t1	hxNCBI_GNO_920044	N:110ds,105dn,84gp,1179tn	T:30ds,14dn,7gp,100tn 	B:32ds,14dn,19gp,100tn	M:26ds,12dn,0gp,100tn 	E:33ds,18dn,2gp,100tn	len:1107,1167


# highest ident inpar
daphplx_hxAUG25p2s4g113t1	daphplx_hxNCBI_GNO_920044	ARP6d4386	0.0	95	inpar1	92.7
daphplx_hxNCBI_GNO_920044	daphplx_hxAUG25p2s4g113t1	ARP6d4386	0.0	95	inpar1	92.7
daphplx_hxAUG25p2s4g115t1	daphplx_hxAUG25p2s4g113t1	ARP6d4386	0.0	87	inpar1	100
daphplx_hxJGI_V11_92880	daphplx_hxJGI_V11_93831	ARP6d4386	0.0	87	inpar1	100
daphplx_hxJGI_V11_93831	daphplx_hxJGI_V11_92880	ARP6d4386	0.0	87	inpar1	100
daphplx_hxAUG25p2s4g111t1	daphplx_hxNCBI_GNO_920044	ARP6d4386	0.0	88	inpar1	93.7
daphplx_hxAUG25p2s4g110t1	daphplx_hxAUG25p2s4g109t1	ARP6d4386	0.0	85	inpar1	97.4

j=hxAUG25p2s4g113t1
for i in hxAUG25p2s4g109t1 hxAUG25p2s4g110t1 hxAUG25p2s4g111t1 hxAUG25p2s4g115t1 hxAUG25s13g137t1 hxAUG26res13g1t1 hxJGI_V11_301607 hxJGI_V11_92880 hxJGI_V11_93831 hxNCBI_GNO_1448043 hxNCBI_GNO_264133 hxNCBI_GNO_920044; 
do { cat ${pt}${j}.cds ${pt}${i}.cds | $gs/muscle/muscle -quiet -clw | $evigene/scripts/prot/cdsalnpatt.pl -pct; } done

# highest id
hxAUG25p2s4g113t1	hxNCBI_GNO_920044	N:31ds,18dn,141gp,1173tn	T:9ds,3dn,12gp,100tn  	B:6ds,3dn,34gp,100tn  	M:13ds,3dn,0gp,100tn  	E:7ds,2dn,2gp,100tn	len:1038,1167
hxAUG25p2s4g113t1	hxAUG25p2s4g111t1	N:58ds,45dn,9gp,1041tn  	T:17ds,7dn,1gp,100tn  	B:18ds,7dn,3gp,100tn  	M:17ds,3dn,0gp,100tn  	E:16ds,9dn,0gp,100tn	len:1038,1035
hxAUG25p2s4g113t1	hxAUG25p2s4g115t1	N:44ds,50dn,2gp,1039tn  	T:13ds,7dn,0gp,100tn  	B:3ds,1dn,0gp,100tn   	M:18ds,6dn,0gp,100tn  	E:18ds,15dn,1gp,100tn	len:1038,1038
hxAUG25p2s4g113t1	hxJGI_V11_93831 	N:39ds,169dn,17gp,1045tn	T:11ds,25dn,2gp,100tn 	B:11ds,17dn,4gp,100tn 	M:15ds,30dn,0gp,100tn 	E:9ds,27dn,1gp,100tn	len:1038,1035
hxAUG25p2s4g113t1	hxJGI_V11_92880 	N:111ds,70dn,0gp,1038tn 	T:32ds,10dn,0gp,100tn 	B:23ds,7dn,0gp,100tn  	M:35ds,13dn,0gp,100tn 	E:37ds,10dn,0gp,100tn	len:1038,1038
hxAUG25p2s4g113t1	hxNCBI_GNO_1448043	N:80ds,148dn,11gp,1042tn	T:23ds,22dn,1gp,100tn	B:33ds,16dn,3gp,100tn 	M:26ds,22dn,1gp,100tn 	E:11ds,27dn,0gp,100tn	len:1038,1035

hxAUG25p2s4g113t1	hxAUG25p2s4g109t1	N:107ds,104dn,69gp,1107tn	T:31ds,15dn,6gp,100tn 	B:32ds,13dn,18gp,100tn	M:30ds,15dn,0gp,100tn 	E:31ds,17dn,1gp,100tn	len:1038,1107
hxAUG25p2s4g113t1	hxAUG25p2s4g110t1	N:118ds,110dn,22gp,1058tn	T:34ds,16dn,2gp,100tn 	B:37ds,13dn,4gp,100tn 	M:34ds,19dn,0gp,100tn 	E:32ds,16dn,2gp,100tn	len:1038,1056
hxAUG25p2s4g113t1	hxAUG25s13g137t1	N:77ds,167dn,435gp,1077tn	T:36ds,39dn,40gp,100tn	B:31ds,39dn,46gp,100tn	M:40ds,41dn,39gp,100tn	E:36ds,38dn,37gp,100tn	len:1038,681
hxAUG25p2s4g113t1	hxAUG26res13g1t1	N:72ds,140dn,521gp,1093tn	T:38ds,37dn,48gp,100tn	B:32ds,20dn,73gp,100tn	M:30ds,37dn,47gp,100tn	E:45ds,42dn,23gp,100tn	len:1038,627
hxAUG25p2s4g113t1	hxJGI_V11_301607	N:69ds,145dn,541gp,1085tn	T:38ds,40dn,50gp,100tn	B:45ds,44dn,36gp,100tn	M:30ds,34dn,54gp,100tn	E:35ds,40dn,60gp,100tn	len:1038,591
hxAUG25p2s4g113t1	hxNCBI_GNO_264133	N:100ds,204dn,361gp,1079tn	T:42ds,43dn,33gp,100tn	B:42ds,47dn,26gp,100tn	M:32ds,41dn,53gp,100tn	E:49ds,39dn,22gp,100tn	len:1038,759

=cut
