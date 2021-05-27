#!/usr/bin/env perl
# tall2agp.pl
# in: $dmag/submitf/pubgenome/reads454f/asmnew/nwbdmag24g7c_allctg-dmag7finall9b.tall6
# out: dmag7gene_ctgmap.agp

while(<>) { 
  next if(/^\W/); s/dmag\w+_contig/contig/; chomp;
  my @tcval=($td,$cd,$tw,$cw,$tsp,$csp)=(split"\t")[0,1,5,6,7,8]; 
  if($td ne $ltd) { putg($ltd,@tcmap) if($ltd); @tcmap=(); }
  push @tcmap, \@tcval;
  $ltd= $td;
}
putg($ltd,@tcmap) if($ltd);

sub _smap { 
  my($ai)= $a->[4] =~ m/(\d+)/;
  my($bi)= $b->[4] =~ m/(\d+)/;
  return ($ai <=> $bi);
}

sub putg {
  my($tdr,@tmap)=@_;
  my $ix=0; $ign++;
  my($sdc,$sdt,$nd)= (0) x 10; 
  my (@ctd,%ctg);
  for my $tv (sort _smap @tmap) { 
    my($td,$cd,$tw,$cw,$tsp,$csp)=@$tv;
    push @ctd,$cd;
    my($tsw,$tsx)=split"/",$tsp; my($tsb,$tse)=split"-",$tsw; my $tsww=1+$tse-$tsb;
    my($csx,$cso)=split":",$csp;
    my @tsx=split",",$tsx; my @csx=split",",$csx;
    my($csb,$cse);
    if($cso eq "-") {
      ($cse)= $csx[0]=~m/(\d+)$/; 
      ($csb)= $csx[-1]=~m/^(\d+)/;
    } else {
      ($csb)= $csx[0]=~m/^(\d+)/; 
      ($cse)= $csx[-1]=~m/(\d+)$/;
    }
    ($csb,$cse)=($cse,$csb) if($csb>$cse); 
    my $csw=1+$cse-$csb;
    my $nx=@tsx - 1; my($ltxb,$lcxb)=(0,0);
    print "#gcm$ign\t$td:$tsw,$tsww/$tw\t$cd:$csb-$cse:$cso,$csw/$cw\n";
    for my $i (0..$nx) { 
      $ix++; my $tx=$tsx[$i]; my $cx=$csx[$i]; 
      print join("\t",$td,$cd,$ix,$tx,$cx,$cso)."\n"; 
      my($txb,$cxb)= map{ m/(\d+)/; $1; } ($tx,$cx);
      if($ltxb>0) { my $dt= 1+$txb-$ltxb; $sdt+= $dt; $dt=1+abs($cxb-$lcxb); $sdc+=$dt; $nd++; }
      ($ltxb,$lcxb)=($txb,$cxb);
      if($i==0 or $i==$nx) { push @{$ctg{$cd}}, $cxb,$cso,$cw,$txb; }
    } 
  }
  
  my($adt,$acd)=(0,100);
  if($nd>0) { $adt= int($sdt/$nd); $adc= int($sdc/$nd); }
  my @cd= @ctd; # (sort keys %ctg); 
  my $ncd=@cd - 1;
  print "#cjoin$ign\t$ncd,$adt,$adc\t@cd\n";
  for my $i (0..$ncd-1) { 
    $icd=$cd[$i]; my @ic= @{$ctg{$icd}};
    $k=(@ic<6)?0:4;
    my($icb,$ico,$ilen,$itb)=@ic[$k,$k+1,$k+2,$k+3];
    my $j=$i+1;
    $jcd= $cd[$j]; my @jc= @{$ctg{$jcd}}; 
    $k=0;
    my($jcb,$jco,$jlen,$jtb)= @jc[$k,$k+1,$k+2,$k+3];
    my $cdist= ($ico eq "-")? $icb : $ilen-$icb;
    $cdist+= $adc;
    $cdist+= ($jco eq "-")? $jlen-$jcb : $jcb;
    print join("\t",$icd,$icb,$ico,$jcd,$jcb,$jco,$cdist)."\n";
   
    #   if(0) {
    #   $itw= $ic[5] - $ic[2]; # tw[n]-tw[0]
    #   $icw= $ic[3] - $ic[0]; # cw[n] - cw[0]
    #   #for my $j ($i+1..$#cd) 
    #   {  
    #     my(@tic,@tjc,$cdist,$jcd);
    #     $jcd= $cd[$j]; my @jc= @{$ctg{$jcd}};
    #     $jtw= $jc[5] - $jc[2]; 
    #     $jcw= $jc[3] - $jc[0];
    #     #   $cdist= $icw + ($jc[2] - $ic[5]); # $icw + $jcw + $adc;
    #     #   if($cdist<0) { @tic= ($jcd,@jc); @tjc= ($icd,@ic); $cdist= -$cdist; }
    #     #   else { @tic= ($icd,@ic); @tjc= ($jcd,@jc);}
    #     my $rev= ($ic[2] > $jc[2])?1:0;
    #     if($rev) {
    #       $cdist= abs($ic[0] - 0) + abs($jcw) + $adc;
    #       @tic= ($jcd,@jc); @tjc= ($icd,@ic);
    #     } else {
    #       $cdist= abs($ic[3] - $ic[0]) + abs($jc[0] - 0) + $adc;
    #       @tic= ($icd,@ic); @tjc= ($jcd,@jc);
    #     }
    #     print join("\t",@tic[0..2],@tjc[0..2],$cdist)."\n";
    #   }
    #   }
    
  }
  print "\n";
  
}
__END__

=item in tall6 gene x contig map

Query	Source	Bits	Ident	Align	Qlen	Slen	Qspan	Sspan

Dapma7bEVm000001t1	dmag24nwb7c_contig05591	10717	5844	5861	8560	8328	2741-8560/2741-2851,2850-2999,2998-3148,3147-3448,3447-3668,3665-3910,3910-4066,4066-4195,4194-4621,4620-4853,4849-5031,5031-5198,5196-5359,5358-5517,5515-5657,5655-5823,5823-6235,6235-6303,6309-7884,7884-8136,8134-8298,8297-8560	8218-8328,7996-8145,7791-7941,7433-7734,7145-7366,6838-7083,6604-6760,6415-6544,5929-6356,5628-5861,5390-5572,5154-5321,4739-4902,4517-4676,4010-4152,3574-3742,3093-3505,2754-2822,969-2546,646-898,428-592,110-373:-
Dapma7bEVm000001t1	dmag24nwb7c_contig05592	2662	1447	1451	8560	9807	1304-2745/1304-1447,1446-1702,1703-1848,1846-2272,2272-2558,2560-2681,2678-2745	9490-9633,9140-9396,8905-9050,8420-8846,4630-4916,132-253,1-68:-
Dapma7bEVm000001t1	dmag24nwb7c_contig05593	2425	1315	1317	8560	30074	1-1304/1-280,280-414,412-626,623-763,761-922,922-1072,1072-1304	14488-14767,14183-14317,10013-10227,9467-9607,6612-6773,5907-6057,4999-5231:-
Dapma7bEVm000002t1	dmag24nwb7c_contig06689	15462	8388	8390	8697	18045	1-8697/1-796,783-1130,1130-1242,1496-1738,1738-1906,1906-2055,2054-7283,7276-8325,8407-8697	11229-12022,9919-10266,9480-9592,8054-8296,7827-7995,7216-7365,1829-7056,423-1472,75-365:-
Dapma7bEVm000003t1	dmag24nwb7c_contig04234	6563	3579	3591	6789	26705	1-3537/1-635,630-705,702-793,793-926,922-1057,1057-1211,1211-1321,1320-1444,1443-1528,1528-1683,1678-2031,2027-2344,2344-2546,2541-2944,2941-3088,3087-3232,3231-3373,3371-3537	20676-21312,21383-21458,21815-21906,22047-22180,22241-22376,22458-22612,22670-22780,22849-22973,23058-23143,23211-23366,23424-23777,23838-24155,24770-24972,25046-25449,25514-25661,25776-25921,26243-26385,26539-26705:+
Dapma7bEVm000003t1	dmag24nwb7c_contig04235	3078	1682	1691	6789	2702	3562-5228/3562-3752,3748-3985,3985-4114,4113-4255,4253-4468,4468-4606,4605-4764,4763-4864,4864-5123,5117-5228	1-191,471-708,774-903,1026-1168,1237-1452,1664-1802,1880-2039,2109-2210,2270-2529,2591-2702:+
Dapma7bEVm000003t1	dmag24nwb7c_contig04236	2736	1512	1520	6789	9023	5247-6789/5247-6352,6384-6500,6496-6789	162-1268,1463-1579,1635-1929:+
Dapma7bEVm000004t1	dmag24nwb7c_contig15334	4116	2227	2227	5895	3892	3237-5895/3237-3464,3463-3785,3782-4140,4138-4856,5298-5895	3665-3892,3259-3581,2826-3184,1741-2459,1039-1636:-
Dapma7bEVm000004t1	dmag24nwb7c_contig16939	3109	1684	1685	5895	2349	1218-3314/1218-1326,1321-1478,1471-1915,1911-2407,2841-3164,3163-3314	68-176,259-416,505-949,1162-1658,1792-2115,2198-2349:+
Dapma7bEVm000004t1	dmag24nwb7c_contig14687	1949	1062	1064	5895	5128	172-1223/172-609,605-837,831-1223	3200-3636,2108-2340,1648-2040:-

=cut
