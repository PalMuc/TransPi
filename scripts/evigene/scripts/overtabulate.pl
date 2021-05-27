#!/usr/bin/env perl
# overtabulate.pl

=item usage

** assumes no overlapped genes, bad if alt tr

set refg=../genes/aphid2_mix3.gff.gz

set rnag=inch/aphid.trinity.1st.gff.gz
set rnag=velmapt7/aphid_vel7asm1st.gff.gz

gzgrep CDS $refg | sed 's/CDS/exon/' | egrep -v '^#|\.mrna' | sort -k1,1 -k4,4n -k5,5n |\
$td/overlapfilter -pass exon -nostrand -over $rnag -in stdin -act markid -mark rnaa -pct 50 | \
env mark=rnaa $td/overtabulate.pl

=item FIXME bugs

  -- assumes input reference locations (e.g. CDS, or exon) are location sorted
     but also grouped as genes, w/ Parent=geneid; 
     Gene joins measured from id transitions.
  -- doesn't handle genes inside genes

  >> should probably collect all input gff rows, 
      then sort by gene-group then location
  >> filter out any alt-tr, e.g. same/overlapped exons    
  
=cut

## should this use base counts?  overlap -act markidbase

use strict;

my $MARK= $ENV{mark} || "rnaa";

my($nx0, $nx1, $nx, $ngene, $nexon, $nglast, $ngnext, $ngcover, $covlast, $covnext, $cover,
    $d1, $d0, $d00) = (0) x 20;
my(%cx0, %cx1, %cx);

my $f= ("%7s ") x 10; $f.="\n"; 
printf $f, qw(nGene neXon GCover GJoin XCover XJoin pGene  pGJoin peXon pXJoin);

while(<>) {
  my($d)=m/Parent=([^;\s]+)/; 
  my($c)=m/$MARK=([^;\s]+)/; 
  if($d1 ne $d) {  
    putd($d0, $d1, $d00) if $d0;	 
    %cx0= %cx1; %cx1= %cx; %cx=();  
    $nx0= $nx1; $nx1= $nx; $nx=0; 
    $d00=$d0; $d0= $d1; $d1= $d; 
    }  
  $nx++; 
  ## this isnt quite right; need to use exon id,
  # map{ $cx{$_}++ } (split",",$c); 
  map{ $cx{$d}{$nx}{$_}++ } (split",",$c); 
}

#END
## putd($d1) if ($d eq $d1);	 

$ngcover||=1; $cover||=1; $ngene||=1; $nexon||=1; 

my $ngo=($nglast+$ngnext);    my $pgoc= int(1000*$ngo/$ngcover)/10; 
my $xgo=($covlast+$covnext);  my $pxoc= int(1000*$xgo/$cover)/10; 
my $pgc=int(1000*$ngcover/$ngene)/10; my $pxc=int(1000*$cover/$nexon)/10; 

my $f= ("%7d ") x 6; $f.= ("%7.1f ") x 4; $f.="\n"; 
printf $f, $ngene, $nexon, $ngcover, $ngo, $cover, $xgo,
           $pgc, $pgoc, $pxc, $pxoc; 


sub putd { 
  my($thisid, $nextid, $lastid)= @_;
  $ngene++;  $nexon += $nx1;  
  
  my @xoverid=();
  my @xi= sort keys %{$cx1{$thisid}};
  foreach my $xi (@xi) { 
    my @xo= sort keys %{$cx1{$thisid}{$xi}};
    $cover++ if(@xo); # this exon is covered
    push( @xoverid, @xo); 
  }
  $ngcover++ if(@xoverid); # this gene is covered, at least partly
  return unless(@xoverid);
  
  my %xoverid= map{ $_,1 } @xoverid;
  
  my $over=0; 
  my @xk= sort keys %{$cx0{$lastid}};
  foreach my $xi (@xk) { 
    my @xo= sort keys %{$cx0{$lastid}{$xi}};
    my $no=0; map{ $no++ if($_); } @xoverid{@xo}; $over++ if($no>0);  
  }
  $nglast++ if $over; $covlast+= $over;
  
  @xk= sort keys %{$cx{$nextid}};
  foreach my $xi (@xk) { 
    my @xo= sort keys %{$cx{$nextid}{$xi}};
    my $no=0; map{ $no++ if($_); } @xoverid{@xo}; $over++ if($no>0);   
  }
  $ngnext++ if $over; $covnext+= $over;
 
}
  

sub putd_OLD { 
  my($thisid)= @_;
  $ngene++;  $nexon += $nx1;  
  my($c1)= sort{$cx1{$b}<=>$cx1{$a}} keys %cx1; 
  if($c1) { 
    if(my $nc1= $cx1{$c1}) { 
      $ngcover++; $cover += $nc1; 
      if(my $c= $cx0{$c1}) { $nglast++;  $covlast += $c; } 
      if(my $c= $cx{$c1}) { $ngnext++;  $covnext += $c; }  
    } 
  } 
}
 
__END__

=item output table

# ** use refseq/gnomon for full genome gene set
set refg=../genes/aphid2_mix3.gff.gz

##set refg=../refseq/acyr1-ncbignomon.gmap.gff.gz << has some alt tr??
set refg=acyr1-ncbignomon.uniqcds.gff  # removed overlapped cds exons
gzcat ../refseq/acyr1-ncbignomon.gmap.gff.gz | grep '^Scaff' | grep CDS | sort -k1,1 -k4,4n -k5,5n | \
perl -ne'($r,$b,$e)=(split)[0,3,4]; ($d)=m/Parent=([^;\s]+)/; if($skip{$d}) { $ndup++; } \
elsif($lr eq $r and $b < $le and $e > $lb) { $ndup++; $skip{$d}++; } else { print;  \
($lr,$lb,$le,$ld)=($r,$b,$e,$d); } $ll=$_; ' \
> acyr1-ncbignomon.uniqcds.gff

##set rnag=aphid_rnaseq.all11cuff8.best1.gff.gz
##set rnag=aphid_rnaseq.velbest1.gff.gz
##set rnag=velmapt6/aphid_vel6asm.gff

set rnag=aphid_rnaseq.velasm5.best1.gff.gz
set rnag=inch/aphid_rnaseq.iwormg70first.gff.gz

set rnag=inch/aphid.trinity.1st.gff.gz
set rnag=inch/aphid.trinity.gff.gz
set rnag=velmapt7/aphid_vel7asm1st.gff.gz
set rnag=velmapt7/aphid_vel7asm.gff.gz 
set rnag=velmapt7/aphid_rnaseq.all17cuff8.1st.gff.gz
set rnag=../epasa2/pasa_out2/pasa2_aphid2.asmbl_bestgenes.gff.gz

gzcat $refg  | grep '^Scaff'| grep CDS | sed 's/CDS/exon/' | sort -k1,1 -k4,4n -k5,5n |\
$td/overlapfilter -pass "exon" -nostrand -over $rnag -in stdin -act markid -mark rnaa -pct 50 | \
env mark=rnaa $td/overtabulate.pl

B1.  refg=../genes/aphid2_mix3.gff.gz, CDS exons

  nGene   neXon   GOver   GJoin  XOver    XJoin   pGene  pGJoin   peXon  pXJoin 
  34806  168966   26724    7247  111917   27020    76.7+   27.1    66.2    24.1-  : cufflinks87(1st-map)

  34806  168966   22863    5797  104580   17569    65.6    25.3    61.8    16.7   : velasm5(best)
  34806  168966   24941    3955  109868    4476    71.6+   15.8    65.0+    4.0+  : vel7asm(all-map)
  34806  168966   20800    3709   99956    4177    59.7    17.8    59.1     4.1   : vel7asm(1st-map) //

  34806  168966   24486    8607  113059   36517    70.3+   35.1    66.9    32.2-  : inchworm(1st-map)  
  34806  168966   17796    3034   88780    3225    51.1    17.0    52.5     3.6+  : trinity(all-map)
  34806  168966   15315    2897   83557    3070    44.0    18.9    49.4     3.6   : trinity(1st-map) //

  34806  168966   23651    6261   99953   21326    67.9    26.4    59.1    21.3   : pasa2vel7-best


B2.  refg=../refseq/acyr1-ncbignomon.gmap.gff.gz, CDS exons   
   
  nGene   neXon   GOver   GJoin  XOver    XJoin   pGene  pGJoin   peXon  pXJoin 
  35990  131268   22061    6342   87805   24746    61.2    28.7    66.8    28.1   : cufflinks87(1st-map)
  35990  131268   18780    2294   86194    2736    52.1+   12.2    65.6     3.1+  : vel7asm(all-map)
  35990  131268   20387    7272   91993   32946    56.6    35.6    70.0    35.8   : inchworm(1st-map)  
  35990  131268   13469    1645   73511    1813    37.4    12.2    56.0     2.4+  : trinity(all-map)
  35990  131268   17459    5073   78632   18789    48.5    29.0    59.9    23.8   : pasa2vel7-best


A1. Refseq protein blastp, n=10467
      Method  Found   BitAve  pFull 
cufflinks87   10423   807     89    
vel7asm       10409   676     77    
inchworm      10431   816     91    
trinity       10378   791     89    
pasa2-velv7   10385   746     86    

B1. ref. genes= aphid2_mix3,  nGene=34806  nCDSexon=168966 

      Method  pGene  pGJoin   peXon  pXJoin 
cufflinks87    76.7+   27.1    66.2    24.1- 
vel7asm        71.6+   15.8    65.0+    4.0+ 
inchworm       70.3+   35.1    66.9    32.2- 
trinity        51.1    17.0    52.5     3.6+ 
pasa2-velv7    67.9    26.4    59.1    21.3  


#....... Summary of best aphid rna-seq assemblies ......................

Matching reference gene sequences from aphid version 1.

              A1. Refseq protein         A2.  Validated refseq transcripts
                 (blastp, n=10467)           (blastn,   n=1217)
 Rna-assembly Found   BitAve  pFull      Found   BitAve  pFull 
cufflinks87   10423   807+    89         1213    1651    73    
vel7asm       10409   676     77         1217    1357    63    
vel5asm       10414   664     78         1216    1429    65 
inchworm      10431   816+    91         1217    1699+   75     
trinity       10378   791+    89         1215    1719+   76    
pasa2-velv7   10385   746     86         1211    1749+   78    
  Found = # found at e <=1e-5, BitAve = ave. bitscore, 
  pFull = average percent of gene matched

Matching all gene locations and joining 2+ genes

               B1. Ref genes=aphid2_mix3        B2.  Ref genes=acyr1-ncbignomon.unique
               (nGene=34806, nCDS=168966)        (nGene=35990, nCDS=131268) 
 Rna-assembly  pGene  pGJoin   peXon  pXJoin     pGene  pGJoin   peXon  pXJoin 
cufflinks87    76.7+   27.1    66.2    24.1-     61.2+   28.7    66.8    28.1-
vel7asm        71.6+   15.8    65.0     4.0+     52.1+   12.2    65.6     3.1+ 
vel5asm        65.6    25.3    61.8    16.7      -- not done
inchworm       70.3+   35.1    66.9    32.2-     56.6+   35.6    70.0    35.8-
trinity        51.1    17.0    52.5     3.6+     37.4    12.2    56.0     2.4+ 
pasa2-velv7    67.9    26.4    59.1    21.3      48.5    29.0    59.9    23.8  

    pGene= % genes found, pGJoin= % found with overlap to neighbor gene CDS
    peXon= % CDS exon found, pXJoin = % found w/ overlap to neighbor gene CDS

Assembly methods:
  cufflinks87 : cufflinks v0.8 using 2011.02 aphid rna-seq data (400 Mill mapped reads)
  vel7asm     : velvet/oases assembly of 2011.02 aphid rna-seq+EST data, 
  vel5asm     : velvet using 2011.01 aphid rna-seq+EST (200 Mill mapped reads)
  inchworm    : inchworm assembly of 2011.02 aphid rna-seq data
  trinity     : trinity assembly  of 2011.02 aphid rna-seq+EST 
  pasa2-velv7 : PASA 2 assembly of vel7asm
  
#........................................

=cut