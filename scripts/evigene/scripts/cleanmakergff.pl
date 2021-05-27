#!/usr/bin/env perl
# cleanmakergff.pl
# need also expand exon,CDS with 2+ parent IDs, convert to complete mRNA records
# oops, not good enough, need to collect all of gene record and expand mRNAs as sequential mRNA1/exons1, mRNA2/exons2..

while(<>) {
  unless(/^\w/) { print;next; }
  chomp; @r=split"\t"; ($ref,$attr)=@r[0,-1];   
  ($id)=$attr=~m/ID=([^;\s]+)/g;  
  ($pa)=$attr=~m/Parent=([^;\s]+)/g;  
  @pa= ($pa) ? split",",$pa : ();  
  map{ s/\W//g; } ($ref,$id,@pa); # remove non-alphanumerics from IDs
  $attr=~s/ID=([^;\s]+)/ID=$id/; 
  if(@pa) {
    foreach $pa (@pa) {  # print exon for each Parent record, should be only 1 parent for mRNA
      $attr=~s/Parent=([^;\s]+)/Parent=$pa/; 
      push @{$gene{$pa}}, join("\t",$ref,@r[1..7],$attr)."\n";
      }
    }
  else {
    putg() if(%gene); %gene=();
    push @{$gene{1}}, join("\t",$ref,@r[1..7],$attr)."\n";
  }
}

sub putg {
  for my $p (sort keys %gene) {
    my @r= @{$gene{$p}}; print @r;
  }
}

=item new variant to regularize maker.gff to evigene

 2015.12 for turquoise killifish

perl -ne 'BEGIN{$IDP=$ENV{idp}||"Notfur1g"; } if(/^\W/){ next; }  
@v=split"\t"; $ok=0; ($r,$s,$t,$b,$e,$v,$o,$xp,$an)=@v; 
if($t eq "mRNA") { ($not)=m/;Note=([^;\n]+)/; ($cl)= $not=~m/Tier\-(\w+)/; ($nad)=m/Name=([^;\n]+)/; 
$ti= ($nad =~ s/.mRNA\-(\w+)//)? $1 :0; ($td)= m/ID=(\w+)/; ($gd)=m/Parent=(\w+)/; $nad=~s/.\d+of\d+.*//; 
$gid=$IDP.$gd."t$ti"; $v[8]="ID=$IDP$td;gid=$gid;alt=$ti;tscore=$cl;Name=$nad\n"; $ok=1; } 
elsif($t =~/^(exon|CDS)/) { ($pd)=m/Parent=(\w+)/; $v[8]="Parent=$IDP$pd\n"; $ok=1; } 
$v[0]=~s/GapFilled//; print join("\t",@v) if($ok);' \
  NotFur1_protein_coding_gene_models_15-07-2014.gff3 > NotFur1_mrna1407.gff


  mrna.maker.gff type count
247070 CDS
20710 contig
257582 exon
26792 five_prime_UTR
28494 gene
30028 mRNA
17183 start_codon
23068 stop_codon
20752 three_prime_UTR

=cut

__DATA__

scaffold6204.1|size42157        maker   gene    8074    19935   .       -       .       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9;Name=maker-scaffold6204%252E1%257Csize42157-augustus-gene-0.9;Note
=Similar to GFPT2: Glucosamine--fructose-6-phosphate aminotransferase [isomerizing] 2 (Homo sapiens)
scaffold6204.1|size42157        maker   mRNA    8074    19935   .       -       .       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-1;Parent=maker-scaffold6204%2E1%7Csize42157-augustus-gene-0.9
;Name=maker-scaffold6204%252E1%257Csize42157-augustus-gene-0.9-mRNA-1;_AED=0.07;_eAED=0.07;_QI=30|0.90|0.
91|1|0.90|0.83|12|60|714;Note=Similar to GFPT2: Glucosamine--fructose-6-phosphate aminotransferase [isome
rizing] 2 (Homo sapiens)
scaffold6204.1|size42157        maker   mRNA    8557    19935   .       -       .       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-2;Parent=maker-scaffold6204%2E1%7Csize42157-augustus-gene-0.9
;Name=maker-scaffold6204%252E1%257Csize42157-augustus-gene-0.9-mRNA-2;_AED=0.06;_eAED=0.06;_QI=30|1|1|1|1
|1|11|168|695;Note=Similar to GFPT2: Glucosamine--fructose-6-phosphate aminotransferase [isomerizing] 2 (
Bos taurus)
scaffold6204.1|size42157        maker   exon    19650   19935   1.01    -       .       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-1:exon:60619;Parent=maker-scaffold6204%2E1%7Csize42157-august
us-gene-0.9-mRNA-1,maker-scaffold6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-2
scaffold6204.1|size42157        maker   exon    18461   18633   1.01    -       .       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-1:exon:60620;Parent=maker-scaffold6204%2E1%7Csize42157-august
us-gene-0.9-mRNA-1,maker-scaffold6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-2
scaffold6204.1|size42157        maker   exon    14823   14931   1.01    -       .       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-1:exon:60621;Parent=maker-scaffold6204%2E1%7Csize42157-august
us-gene-0.9-mRNA-1,maker-scaffold6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-2
scaffold6204.1|size42157        maker   exon    13316   13461   1.01    -       .       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-1:exon:60622;Parent=maker-scaffold6204%2E1%7Csize42157-august
us-gene-0.9-mRNA-1,maker-scaffold6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-2
scaffold6204.1|size42157        maker   exon    12000   12164   1.01    -       .       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-1:exon:60623;Parent=maker-scaffold6204%2E1%7Csize42157-august
us-gene-0.9-mRNA-1,maker-scaffold6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-2
scaffold6204.1|size42157        maker   exon    10301   10664   1.01    -       .       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-1:exon:60624;Parent=maker-scaffold6204%2E1%7Csize42157-august
us-gene-0.9-mRNA-1,maker-scaffold6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-2
scaffold6204.1|size42157        maker   exon    9956    10206   1.01    -       .       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-1:exon:60625;Parent=maker-scaffold6204%2E1%7Csize42157-august
us-gene-0.9-mRNA-1,maker-scaffold6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-2
scaffold6204.1|size42157        maker   exon    9771    9885    1.01    -       .       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-1:exon:60626;Parent=maker-scaffold6204%2E1%7Csize42157-august
us-gene-0.9-mRNA-1,maker-scaffold6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-2
scaffold6204.1|size42157        maker   exon    9396    9651    1.01    -       .       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-1:exon:60627;Parent=maker-scaffold6204%2E1%7Csize42157-august
us-gene-0.9-mRNA-1,maker-scaffold6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-2
scaffold6204.1|size42157        maker   exon    9143    9311    1.01    -       .       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-1:exon:60628;Parent=maker-scaffold6204%2E1%7Csize42157-august
us-gene-0.9-mRNA-1,maker-scaffold6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-2
scaffold6204.1|size42157        maker   exon    8729    8808    1.01    -       .       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-1:exon:60629;Parent=maker-scaffold6204%2E1%7Csize42157-august
us-gene-0.9-mRNA-1
scaffold6204.1|size42157        maker   exon    8074    8194    1.01    -       .       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-1:exon:60630;Parent=maker-scaffold6204%2E1%7Csize42157-august
us-gene-0.9-mRNA-1
scaffold6204.1|size42157        maker   exon    8557    8808    1.01    -       .       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-2:exon:60631;Parent=maker-scaffold6204%2E1%7Csize42157-august
us-gene-0.9-mRNA-2
scaffold6204.1|size42157        maker   CDS     8134    8194    .       -       1       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-1:cds:61197;Parent=maker-scaffold6204%2E1%7Csize42157-augustu
s-gene-0.9-mRNA-1
scaffold6204.1|size42157        maker   CDS     8729    8808    .       -       0       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-1:cds:61196;Parent=maker-scaffold6204%2E1%7Csize42157-augustu
s-gene-0.9-mRNA-1
scaffold6204.1|size42157        maker   CDS     9143    9311    .       -       1       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-1:cds:61195;Parent=maker-scaffold6204%2E1%7Csize42157-augustu
s-gene-0.9-mRNA-1
scaffold6204.1|size42157        maker   CDS     9396    9651    .       -       2       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-1:cds:61194;Parent=maker-scaffold6204%2E1%7Csize42157-augustu
s-gene-0.9-mRNA-1
scaffold6204.1|size42157        maker   CDS     9771    9885    .       -       0       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-1:cds:61193;Parent=maker-scaffold6204%2E1%7Csize42157-augustu
s-gene-0.9-mRNA-1
scaffold6204.1|size42157        maker   CDS     9956    10206   .       -       2       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-1:cds:61192;Parent=maker-scaffold6204%2E1%7Csize42157-augustu
s-gene-0.9-mRNA-1
scaffold6204.1|size42157        maker   CDS     10301   10664   .       -       0       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-1:cds:61191;Parent=maker-scaffold6204%2E1%7Csize42157-augustu
s-gene-0.9-mRNA-1
scaffold6204.1|size42157        maker   CDS     12000   12164   .       -       0       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-1:cds:61190;Parent=maker-scaffold6204%2E1%7Csize42157-augustu
s-gene-0.9-mRNA-1
scaffold6204.1|size42157        maker   CDS     13316   13461   .       -       2       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-1:cds:61189;Parent=maker-scaffold6204%2E1%7Csize42157-augustu
s-gene-0.9-mRNA-1
scaffold6204.1|size42157        maker   CDS     14823   14931   .       -       0       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-1:cds:61188;Parent=maker-scaffold6204%2E1%7Csize42157-augustu
s-gene-0.9-mRNA-1
scaffold6204.1|size42157        maker   CDS     18461   18633   .       -       2       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-1:cds:61187;Parent=maker-scaffold6204%2E1%7Csize42157-augustu
s-gene-0.9-mRNA-1
scaffold6204.1|size42157        maker   CDS     19650   19905   .       -       0       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-1:cds:61186;Parent=maker-scaffold6204%2E1%7Csize42157-augustu
s-gene-0.9-mRNA-1
scaffold6204.1|size42157        maker   CDS     8725    8808    .       -       0       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-2:cds:61208;Parent=maker-scaffold6204%2E1%7Csize42157-augustu
s-gene-0.9-mRNA-2
scaffold6204.1|size42157        maker   CDS     9143    9311    .       -       1       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-2:cds:61207;Parent=maker-scaffold6204%2E1%7Csize42157-augustu
s-gene-0.9-mRNA-2
scaffold6204.1|size42157        maker   CDS     9396    9651    .       -       2       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-2:cds:61206;Parent=maker-scaffold6204%2E1%7Csize42157-augustu
s-gene-0.9-mRNA-2
scaffold6204.1|size42157        maker   CDS     9771    9885    .       -       0       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-2:cds:61205;Parent=maker-scaffold6204%2E1%7Csize42157-augustu
s-gene-0.9-mRNA-2
scaffold6204.1|size42157        maker   CDS     9956    10206   .       -       2       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-2:cds:61204;Parent=maker-scaffold6204%2E1%7Csize42157-augustu
s-gene-0.9-mRNA-2
scaffold6204.1|size42157        maker   CDS     10301   10664   .       -       0       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-2:cds:61203;Parent=maker-scaffold6204%2E1%7Csize42157-augustu
s-gene-0.9-mRNA-2
scaffold6204.1|size42157        maker   CDS     12000   12164   .       -       0       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-2:cds:61202;Parent=maker-scaffold6204%2E1%7Csize42157-augustu
s-gene-0.9-mRNA-2
scaffold6204.1|size42157        maker   CDS     13316   13461   .       -       2       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-2:cds:61201;Parent=maker-scaffold6204%2E1%7Csize42157-augustu
s-gene-0.9-mRNA-2
scaffold6204.1|size42157        maker   CDS     14823   14931   .       -       0       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-2:cds:61200;Parent=maker-scaffold6204%2E1%7Csize42157-augustu
s-gene-0.9-mRNA-2
scaffold6204.1|size42157        maker   CDS     18461   18633   .       -       2       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-2:cds:61199;Parent=maker-scaffold6204%2E1%7Csize42157-augustu
s-gene-0.9-mRNA-2
scaffold6204.1|size42157        maker   CDS     19650   19905   .       -       0       ID=maker-scaffold
6204%2E1%7Csize42157-augustus-gene-0.9-mRNA-2:cds:61198;Parent=maker-scaffold6204%2E1%7Csize42157-augustu
s-gene-0.9-mRNA-2
