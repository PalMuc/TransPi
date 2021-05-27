#!/usr/bin/env perl
# rxpastecounts.pl

# egrep '^HS_.._X_' ../dmrna5.groups.info | head
# HS_BN_X_r1_45.bam  << group names not same as xprs folder names, need num or rename
# HS_BN_X_r2_74.bam

## FIXME: two files/1 group name, either add or choose one; now picks last
# read HS_CO_X_r4 = 77, rxdmag5xau13c2011-Dman_77_TTAGGC_L003/results.xprs
# read HS_CO_X_r4 = 77, rxdmag5xau13c2011-Dman_77_TTAGGC_L004/results.xprs

# head -5  rxdmag5xau13c2011-Dman_32_ACAGTG_L003/results.xprs
# bundle_id       target_id       length  eff_length      tot_counts      uniq_counts     est_counts      eff_counts      ambig_distr_alpha       ambig_distr_beta        fpkm    fpkm_conf_low   fpkm_conf_high  solvable
# 1       m8AUGep24bs03129g212t1  1521    4454.065741     6       6       6.000000        2.048915        0.000000e+00    0.000000e+00    2.515649e-02    2.515649e-02    2.515649e-02    T
# 2       dmag4vel4xc3k45Loc10985t1       3216    8615.101039     32      8       31.297773       11.683396       5.027594e-01    1.515387e-02    6.784353e-02    6.723559e-02    6.845148e-02    T
# 2       m8AUGepir2s02581g145t1  4899    10804.602792    27      3       3.702227        1.678656        4.747783e-02    1.575170e+00    6.398964e-03    6.133130e-03    6.664799e-03    T
# 3       dmag5nun1soapk25loc179659t1     367     72.548421       0       0       0.000000        0.000000        0.000000        0.000000        0.000000e+00    0.000000e+00    0.000000e+00    T

## optional idmap=file of newid,oldid to replace target_id

$usage="env xcol=tot_counts gpatt='^HS_.._X_' gfile=../dmrna5.groups.info rxpastecounts.pl"
  ." rxdmag5xau13c2011-Dman*/results.xprs\n";

$debug=$ENV{debug}||0;
$gfile=$ENV{gfile};
$idmap=$ENV{idmap};
$gpatt=$ENV{gpatt};
$fpatt=$ENV{fpatt}; ## || "Dman_"; ## "nodamx"
$idcol=$ENV{idcol} || 2; # fixme for rna-express, other count tables have id in col1
$xcol= $ENV{xcol};

die "ERR: missing opt\n".$usage unless($gfile and $gpatt and $xcol);
open(F,$gfile) or die "ERR: read $gfile\n$usage";
while(<F>) { chomp; if(m/$gpatt/) { s/\.bam//; ($gn,$fn)= m/^(\w+)_(\d+)$/; $fng{$fn}= $gn; 
  # leave to file rename: $fn=~s/^0+(\d)/$1/; $fng{$fn}= $gn; # fixme: fpatt/01 vs fpatt/1 
  warn "group $gn = $fn\n" if $debug; } }
close(F);

if($idmap) {
open(F,$idmap) or die "ERR: read $idmap\n$usage";
while(<F>) { chomp; next unless(/^\w/); ($id,$oid)=split; $idmap{$oid}= $id; } close(F);
}

@xpr= @ARGV;
foreach $xpr (@xpr) {
  # ($xd,$xf)= $xpr =~ m,(.+)/([^/]+)$,;
  # ($xfn)= $xpr =~ m/Dman_(\d+)/; ## FIXME: option # m/nodamx(\d+)/; ## FIXME: option
  ($xfn)= $xpr =~ m/$fpatt(\d+)/; ## FIXME: option .. fixme2: fpatt/01 vs fpatt/1
  ## rxdmag5xau13c2011-nodamx10
  # warn "at $xfn, $xpr\n" if $debug;
  if($gn=$fng{$xfn}) { 
    warn "read $gn = $xfn, $xpr\n" if $debug;
    push @gn, $gn; # group names; check for dups! HS_CO_X_r1  HS_CO_X_r2 HS_CO_X_r4  HS_CO_X_r4<<
    open(X,$xpr);  $hd=<X>; @hd=split"\t",$hd; 
    if($idcol=~/^\d/) { $dcol=$idcol-1; }
    else { for $i (0 .. @hd - 1) { if ($hd[$i] =~ m/$idcol/) { $dcol=$i; last; } } }
    if($xcol=~/^\d/) { $icol=$xcol-1; }
    else { for $i (0 .. @hd - 1) { if ($hd[$i] =~ m/$xcol/) { $icol=$i; last; } } }
    while(<X>) { chomp; @v=split"\t"; $id=$v[$dcol]||0; $x=$v[$icol]||0;  
      if($idmap) { $id=$idmap{$id}||$id; }
      $counts{$id}{$gn}=$x; $counts{'#TOTAL'}{$gn} +=$x; }
    close(X);
  }
}

die "ERR: no groups\n".$usage unless(@gn);
#TOTAL sorts 1st 
@sgn= sort @gn;
print join("\t","GeneID",@sgn)."\n";
for $id (sort keys %counts) {
  ## if($idmap) { $id=$idmap{$id}||$id; } # do before sort
  print $id;
  for $gn (@sgn) { $x=$counts{$id}{$gn}||0; print "\t",$x; }
  print "\n";
}

__END__

=item output

env idmap=dmag5xau13c2011.pubids xcol=tot_counts fpatt=nodamx gpatt='^ND_' gfile=../dmrna5.groups.info \
 ../rxpastecounts.pl rxdmag5xau13c2011-nodamx*/results.xprs > rxdmag5x_ndx.tocounts

env idmap=dmag5xau13c2011.pubids xcol=tot_counts  fpatt='Dman_'  gpatt='^HS_.._X_' gfile=../dmrna5.groups.info \
 ../rxpastecounts.pl rxdmag5xau13c2011-Dman*/results.xprs > rxdmag5x_hsx.tocounts

env idmap=dmag5xau13c2011.pubids xcol=uniq_counts  fpatt='Dman_'  gpatt='^HS_.._X_' gfile=../dmrna5.groups.info \
 ../rxpastecounts.pl rxdmag5xau13c2011-Dman*/results.xprs > rxdmag5x_hsx.uqcounts

head rxdmag5x_ndx.tocounts | cut -f1-7
GeneID  ND_CD_X_r1      ND_CD_X_r2      ND_CD_X_r3      ND_CO_X_r1      ND_CO_X_r2      ND_CO_X_r3
#TOTAL  433309775       735983579       204316354       194838010       149076930       260075392
Dapma5xEG0000001t1      29      61      21      23      24      41
Dapma5xEG0000001t2      17      48      13      13      19      23
Dapma5xEG0000001t3      11      33      10      9       14      14
Dapma5xEG0000001t4      20      50      17      11      19      28
Dapma5xEG0000001t5      20      54      16      16      22      31
Dapma5xEG0000001t6      20      50      17      11      19      27
Dapma5xEG0000002t1      21      25      4       9       9       15
Dapma5xEG0000003t1      4       0       0       1       0       0

head rxdmag5x_hsx.uqcounts | cut -f1,10-16 | head
GeneID  HS_CA_X_r3      HS_CO_X_r1      HS_CO_X_r2      HS_CO_X_r4      HS_CR_X_r1      HS_CR_X_r2      HS_CR_X_r3
#TOTAL  8903504 9350827 9652636 11112256        12555922        9678730 15761631
Dapma5xEG0000001t1      2       4       6       13      5       4       9
Dapma5xEG0000001t2      0       0       0       0       0       0       0
Dapma5xEG0000001t3      1       0       0       1       0       0       0
Dapma5xEG0000001t4      0       0       0       0       0       0       0
Dapma5xEG0000001t5      3       4       2       4       4       2       2
Dapma5xEG0000001t6      0       0       0       0       0       0       0
Dapma5xEG0000002t1      0       0       0       0       0       1       0
Dapma5xEG0000003t1      6       8       10      3       15      11      13

=cut
