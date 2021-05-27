#!/usr/bin/env perl
# evaltab.pl

=item stats

#Prediction: genes/nvit2_ogs12.an2
#Prediction: gene2/cacao9_epir6-augmap.an6ho
#Evidence: est/est
#Evidence: prot/protein
#Evidence: rnas/rnaseq

=item add counts

Coding bases (Mb)                40.4    40.6      41.4     33.9
Exon bases (Mb)                  60.2    61.3      62.4     51.7
Gene count                       30307   31137     34997    28497 + 6500 TEgenes

=cut

my $KEEPCOUNTS=$ENV{count}||0;

while(<>) {
  if(/^Evidence: (\S+)/) { putr() if $ev; $ev=$1; %no=(); }
  elsif(m,^Prediction: gene[s\d]/(\S+),) { $pr=$1; 
    $pr=~s/aphid[a-z\d]*[_-]?//; $pr=~s/acyr[a-z\d]*[_-]?//; #$pr=~s/cacao[a-z\d]*[_]?//; 
    $pr=~s/view//; $pr=~s/\.an\w+//; $pr=~s/\-augmap\S*//; $pr=~s/\.gmap\S*//; $pr =~ s/\.gff.gz//;
    $pr{$pr}++; 
## what? genes/acyr1-dgilbest.gmap.an3.gff.gz > dgil9best
    
  }elsif(/^# base statistics: (.+)/) { my $v=$1; 
  @v=split /\s*,\s*/,$v; map{ my($k,$v)=split" n=",$_; $no{$k}{$pr}=$v; } @v;

  }elsif(/^# ave_baseover=([\d\.]+), ave_pctover=([\d\.]+)/) { 
  $no{poverbase}{$pr}=$1; $no{poverlap}{$pr}=$2;  

  }elsif(/^# sum_basetotal=([\d\.]+)/) {
   $no{basetotal}{$pr}=$1; 

  }elsif(/^ngene=(\d+); genehit=(\d+); gperfect=(\d+)/){
   $no{ngene}{$pr}=$1;  $no{genehit}{$pr}=$2; $no{gperfect}{$pr}=$3;

  }elsif(/^BestTr: Sn=([\d\.]+); Sp=([\d\.]+)/){
   my($sn,$sp)=($1,$2); $sn=sprintf "%.3f",$sn/100; $sp=sprintf "%.3f",$sp/100;
   $no{gsens}{$pr}=$sn; $no{gspec}{$pr}=$sp;

  }elsif(/^protein_homol (\S+) n=(\d+), bits.aa=([\d\.]+), (\S+) n=(\d+), bits.aa=([\d\.]+)/){
   my($t1,$n1,$b1,$t2,$n2,$b2)=($1,$2,$3,$4,$5,$6);
   $no{"ho1_n"}{$pr}=$n1; $no{"ho1_b"}{$pr}=$b1; $no{"ho1_t"}{$pr}=$t1;
   $no{"ho2_n"}{$pr}=$n2; $no{"ho2_b"}{$pr}=$b2; $no{"ho2_t"}{$pr}=$t2;

  }elsif(/^Coding bases: (\S+)/) { $no{cbases}{$pr}=$1; 
  }elsif(/^Exon bases: (\S+)/) { $no{ebases}{$pr}=$1; 
  }elsif(/^Gene.count: (\S+)/) { $no{gcount}{$pr}=$1; 

  }

}

=item add geneacc

Evidence: cDNA_gene_accuracy
ditto: : Prot_gene_accuracy
Prediction: genes/cacao9_mix2
#collect_overlaps=194833
#overlaps found=62131
ngene=16818; genehit=12040; gperfect=4492; nexon=69346; exonhit=62131
                                ^^^^
Exon: Sn=76.7; Sp=76.5; Sp2=25.9; Ave=76.6; hit=16989463; miss=5206378; exonbase=22126534; allover=65375994
BestTr: Sn=76; Sp=56.7; Sp2=25.7; Ave=66.35; hit=16824737; miss=12829214; allbase=22126534;
        ^^^^^^^^^^^^^
=cut

END{ final_putr();}

sub pbase { my $nb=shift;
$Kb=1024; $Mb=$Kb*$Kb; $Gb=$Kb*$Mb;
my @lv=($Gb, $Mb, $Kb); my @lb=("Gb","Mb","Kb");
foreach my $i (0..$#lv) { return sprintf "%.1f".$lb[$i], $nb/$lv[$i] if($nb>$lv[$i]); }
return $nb;
}

sub putr {
  ## sum to %tno, %tpr
  push(@evkeys, $ev) unless( grep /$ev/, @evkeys);
  foreach my $k (sort keys %no) { 
    foreach my $p (sort keys %{$no{$k}}) {
	$tno{$ev}{$k}{$p}= $no{$k}{$p};
	}
  }
}

sub final_putr {
   putr();
   foreach my $ev (@evkeys) {
   %no= %{$tno{$ev}};
   unless($hd++){ @pr= map{s/^\d//;$_} sort map{
          s/^ogs/0ogs/;
	  s/^(aphid0|peaaphid0)/0$1/;
          $_ } keys %pr; 
          print join("\t","Evid.", "Nevd", "Statistic", @pr),"\n"; }

   my $ni=$no{input}{$pr[0]};
   my $nb=pbase( $no{basetotal}{$pr[0]} );

        $ev =~ s,^.*/,,;
        $ev =~ s,_uniq,,;
        $ev =~ s,tarexons,TAR,;
        $ev =~ s,^(est|rna),\U$1,; 
        $ev =~ s,^prot,Prot,; 
        $ev =~ s,^all_evd_specif,Specif,; 
        $ev =~ s,^cDNA_gene_accuracy,ESTgene,;
        $ev =~ s,^Prot_gene_accuracy,Progene,; 
        $ev =~ s,^Gene_coverage,GCover,;
        $ev =~ s,^Protein_Homology,Homolog,; 
        $ev =~ s,^transposon,T'poson,; 
        
# add: Evidence: Protein_Homology
# protein_homol best n=25651, bits/aa=1.061, arabid n=20948, bits/aa=0.966 

  if($ev =~ /Homol/){
  foreach my $k (qw(ho1_n ho1_b ho2_n ho2_b)) {
        my $nv=1; ## $no{ngene}{$pr[0]};
        (my $kt=$k) =~ s/_\w/_t/;
        my $tv=$no{$kt}{$pr[0]};
        my $kv=$k;
        $kv =~ s,\w+_b,$tv.bits/aa,;
        $kv =~ s,\w+_n,$tv.Nfound,;
        print "$ev\t$nv\t$kv";  
        foreach my $p (@pr) { my $v=$no{$k}{$p}; print "\t$v"; }
        print "\n";
        }

  } elsif($ev =~ /GCover/){
  foreach my $k (qw( cbases ebases gcount )) {
        my %kv=(cbases=>"Coding bases", ebases=>"Exon bases", gcount=>"Gene count");
        my $kv= $kv{$k}|| $k; my $nv="--";
        print "$ev\t$nv\t$kv";
        foreach my $p (@pr) { my $v=$no{$k}{$p}; print "\t$v"; }
        print "\n";
        }

  } elsif($ev =~ /ESTgene|Progene/){
  foreach my $k (qw(gperfect genehit gsens gspec)) {
        my $nv=$no{ngene}{$pr[0]};
        my %kv=( gperfect => "Perfect ", genehit=>"Mostly  ", gsens=>"Sensitv.", gspec=>"Specifc.");
        my $kv= $kv{$k} || $k;
        print "$ev\t$nv\t$kv";  
        foreach my $p (@pr) { my $v=$no{$k}{$p}; print "\t$v"; }
        print "\n";
        }

  } else {
  ## fixme: drop poverlap overlaps output unless requested
  foreach my $k (qw(poverlap poverbase overlaps)) { # sort keys %no
        next if ($k =~/input|overset/);
        next if (!$KEEPCOUNTS && $k =~/poverlap|overlaps/);
        # my %kv=( poverlap => "pOver", poverbase=>"baseOv", overlaps=>"nOver");
        my $kv= $k; # $kv{$k} || $k;
        my $nv=($k=~/base/)?$nb:$ni;
        print "$ev\t$nv\t$kv";  # changed ni to nb here; but only for k ~ base
        foreach my $p (@pr) { my $v=$no{$k}{$p}; print "\t$v"; }
        print "\n";
        }
  }
  print "\n";
}
}
