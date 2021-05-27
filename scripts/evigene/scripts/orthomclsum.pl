#!/usr/bin/env perl
# orthomclsum.pl

=item about

Evigene script combining a collection of small scripts to process
OrthoMCL gene group primary results to summary tables.


=cut

use strict;
use Getopt::Long;
use File::Basename;

my $DEBUG=1;

my $MCLPATH="/bio/bio-grid/mb/mcl9/"; # config

# per run configs
my $orun="Jan_25";
my $phyla="arp11u11";
my $spl="antc,anth,aphid,apis2ref,bombusimp,bombusterr,daphnia,drosmel,human,trica,wasp";  
my $idprefix="ARP";


my $optok= GetOptions(
  "config=s", \$config,
  "orun=s", \$orun,
  "phyla=s", \$phyla,
  "idprefix=s", \$idprefix,
#   "changelist=s", \$changelist, # or from config
#   "version|MSRC=s", \$MySRC,
#   "cadd=s", \@configadd,
#   "notbl2asn!", \$notbl2asn, 
  "debug!", \$DEBUG, 
  );

die "usage: orthomclsum.pl -conf=evigene.conf -orun Jan_25    ...
  opts:  ... \n using orun/all_orthomcl.out\n"
  unless($optok and -d $orun and -f "$orun/all_orthomcl.out"); #  and $orun

# evigene_config($config, \@configadd); # always even if $config null

#........................ subs .........................

sub orun_cleanid
{

=item orun_cleanid

  #.. pretty run-specific fixes.. put in config
  perl -pi.old -e'if(/Nasvi2EG/) { s/(Nasvi2EG\w+)p(\d+)/${1}t$2/g; }' ${phyla}_mcl.{gg,bpo}
  perl -pi.old -e's/ \w+:/ /g; if(/Nasvi2EG/) { s/(Nasvi2EG\w+)p(\d+)/${1}t$2/g; }' $orun/all_orthomcl.out
  ## remove dupl prefix  spp:spp_ID
  # perl -pi.old -e's/ \w+:/ /g;' $orun/all_orthomcl.out

=cut

}



sub groupnames
{

=item groupnames

  cat names/*.name.tab plant9_omclgn.tab | env gtag=PLA idprefix=PLA9_G ../arp_condef.pl > plant9_omclgn.consensus_def.txt

  cat names/*.names ${phyla}_omclgn.tab | env gtag=ARP idprefix=ARP9_G ../omclw/arp_condef.pl \
  > ${phyla}_omclgn.consensus_def.txt 
  
  cat ${phyla}_omclgn.consensus_def.txt | env recase=1 count=withID debug=1 \
  $evigene/scripts/bestgenes_puban_wasp.pl | sed 's/^ARP9_G//' | sort -k1,1n | sed 's/^/ARP9_G/; s/TE:TE:/TE:/' \
    > ${phyla}_omclgn.consensus_def.rename.txt
  mv ${phyla}_omclgn.consensus_def.txt  ${phyla}_omclgn.consensus_def.txt0
  mv ${phyla}_omclgn.consensus_def.rename.txt ${phyla}_omclgn.consensus_def.txt 

=cut

}

sub genegroup_ugpdoc
{

=item genegroup_ugpdoc

cat \
plant9_omclgn2sum.tab \
plant9_omclgn.consensus_def.txt  \
names/*.name.tab \
$orun/all_orthomcl.out |\
env xml=0 date=20111112 clade=Eudicotyledons title='Plant gene group' gtag=PLA idprefix=PLA9_G perl ../genegroupbpo.pl \
> plant9_genes.ugp.txt

cat plant9_genes.ugp.txt | egrep 'GeneID|ntaxa|ngene|occur|descript|Thecc1' | perl -pe \
's/^ +//; if(m/similarity:/ and /Thecc/){ ($id)=m/iden: (\d+)/; ($d)=m/(Thecc\w+)/; $_="$d/$id,"; s/^/cacao: / if $s; $s=0;} \
else{ $s=1; print "\n" if s/\s*GeneID:\s//; s/\n/ /;} END{print"\n";} ' \
>  plant9_genes.ugp_brief.txt2

=cut

}


sub genecount
{
  open(F,"$orun/all_orthomcl.out");
  open(OUT,">${phyla}-orthomcl-count.tab");
  # cat $orun/all_orthomcl.out | env spl=$spl perl -ne  
  # $spl=$ENV{spl};
  my @spl=split",",$spl; my@drop{@spl}=(0)x20; 
  print OUT join("\t","OID","Nt","Ng",@spl),"\n";
  while(<F>) {
    if(/^ORTHOMCL/) { 
    my($om,$gn)=split /:\s+/,$_,2; my  @gn= split" ",$gn; @spc{@spl}= (0)x20; 
    my %didgn=(); map{ $sg=$_; ($g1=$sg)=~s/[pt\.]\d+$//; my $isaltNOT=($didgn{$g1}++ > 0)?1:0; 
    my($sp,$gn1)=split"_",$sg,2; $spc{$sp}++ unless($isalt); } @gn; 
    my($og,$g,$t)= $om =~ m/ORTHOMCL(\d+).(\d+) genes,\s*(\d+) taxa/; 
    my @sc= @spc{@spl}; print OUT join("\t",$og,$t,$g,@sc),"\n"; 
    } 
  }
  close(OUT); close(F);
  # > ${phyla}-orthomcl-count.tab
}


sub groupidlists
{

=item groupidlists

# do this not for 11 taxa (all) but drop human, daphnia?
cat ${phyla}-orthomcl-count.tab | env nt=11 gtag=ARP perl -ne\
'($od,$nt,$ng,@c)=split; print "$ENV{gtag}$od\t\n" if($nt==$ENV{nt});' \
	> ${phyla}_omcl11.allgr.oids

ggrep -F -f ${phyla}_omcl11.allgr.oids ${phyla}_omclgn.tab | cut -f2 | sed 's/^/gid /' > ${phyla}_omcl11.allgr.gids

=cut

}

sub groupclasses
{
  open(F,"${phyla}-orthomcl-count.tab");
  open(OUT,">${phyla}-orthomcl-gclass.tab");
# cat ${phyla}-orthomcl-count.tab | env iskip=11 perl -ne
  my $iskip=-1; ## 11; # config; $iskip=$ENV{iskip}; 
  
  ## version 2; but reorder columns:  Groups ...   Genes ...   Grminmax ...
  print OUT join("\t",qw(species nGene Uniq1 UDup Orth1 OrDup nGroup UniqGrp OrGrp OrMis1 Guniq Gmax Gmin)),"\n";
  my( %sg,%sgu,%sgh1,%sghp,%sgOrG,%sgMIS,%suni,%smax,%smin );
  while(<F>) {
    my($se,$nt,$ngbad,$vm,$smax,$vm2,$smin, $vn);
    my @v=split; $se=$#v; if(/^OID/){ @hd=@v; @spp=@v[3..$se]; $MIS1= @spp - 1;  next; } 
    elsif(/^\D/){ next; }
    $nt=$v[1]; $ngbad=$v[2]; @s= map{ $hd[$_]."=".$v[$_]; } (3..$se); 
    $vm=$smax=$vm2=$smin=0; $vn=99999; 
    map{ my($s,$v)=split"="; $sc{$s}++ if($v>0); $sg{$s} += $v; 
    $sgu{$s} += $v if($nt==1); $sgh1{$s} += $v if($nt>1 and $v==1); 
    $sghp{$s} += $v if($nt>1 and $v>1); $sgOrG{$s}++ if($nt>1 and $v>=1); 
    $sgMIS{$s}++ if($v==0 and ($nt==$MIS1 or ($iskip>=0 and $nt==$MIS1-1 and $v[$iskip]==0))); 
    if($v<=$vn){ $smin=($v==$vn)?"$smin,$s":$s; $vn=$v;} 
    if($v>=$vm){ $smax=($v==$vm)?"$smax,$s":$s; $vm2=$vm; $vm=$v;}  }@s; 
    $suni{$smax}++ if($nt==1); 
    $smax{$smax}++ unless($nt < 2 or $smax =~ /,/);  
    $smin{$smin}++ unless($nt < 4 or $smin =~ /,/);  
    }

  my @s=sort keys %sc; foreach my $s (@s) { 
    my($u,$h,$d,$c,$g,$ngu,$ngo1,$ngo2,$ngo0,$nog);
    $u=$suni{$s}||0; $h=$smax{$s}||0; $d=$smin{$s}||0; $c=$sc{$s}; $g=$sg{$s}; 
    $ngu=$sgu{$s}; $ngo1= $sgh1{$s};  $ngo2= $sghp{$s}; $ngo0=$sgMIS{$s}; $nog= $sgOrG{$s};  
    print OUT join("\t",$s, $g,"na",$ngu,$ngo1,$ngo2,  $c,$c-$nog,$nog,$ngo0,  $u,$h,$d),"\n"; 
  }

  close(OUT); close(F);
  
  #> ${phyla}-orthomcl-gclass.tab
  # nother fixup for config
  # perl -pi -e's/bombusimp/bombi/; s/bombusterr/bombt/;s/trica/tribol/; s/apis2ref/apis/; ' ${phyla}-orthomcl-gclass.tab

=item  R:gclass.tab
gtab<- read.delim("plant9-orthomcl-gclass.tab",header=T)
gtab[,"Uniq1"] <- gu <- gtab[,"nGene"] - rowSums(gtab[,c("UDup","Orth1","OrDup")])
gd <- rowSums(gtab[,c("UDup","OrDup")])
gs <- rowSums(gtab[,c("Uniq1","Orth1")])

gtab2 <- cbind(gtab,gd=gd,gs=gs, dp.gd=round(gd/12000,1), dp.gs=round(gs/15000, 1))
rownames(gtab2) <- gtab2[,1]; gtab2 <- gtab2[,-1]
gtab2[rev(order(gtab2$nGene)),]
=cut

}

sub orun_out2gntab
{
  open(F,"$orun/all_orthomcl.out");
  open(OUT,">${phyla}_omclgn.tab");
  my(%alts,$isalt);
  while(<F>) { 
    if(/^ORTHOMCL/) { 
    my($om,$gn)=split /:\s+/,$_,2;  my @gn= split" ",$gn; 
    my($og,$ng,$nt)= $om =~ m/ORTHOMCL(\d+).(\d+) genes,\s*(\d+) taxa/;  
    my %didg=(); map{ $gn=$_; ($sp)=split"_",$gn,2; (my $g1=$gn)=~s/[pt\.]\d+$//;  
    $didgn{$g1}++; my $isaltNOT=($didgn{$g1}>1)?1:0; $alts{$sp}++ if($isalt);
    print OUT "ARP$og\t$gn\n" unless($isalt); } @gn; }
  }
  if(%alts){ print OUT "#alttr-drops: "; foreach my $s (sort keys %alts){ print OUT "$s=$alts{$s}, "; } print OUT "\n"; } 
  close(OUT); close(F);

# cat ${phyla}_omclgn.tab ${phyla}*_mcl.bpo | env idprefix=ARP perl ../omclw/pastebpo.pl > ${phyla}_omclgns2.tab 
# open(F2,"${phyla}*_mcl.bpo");
# open(F1,"${phyla}_omclgn.tab");
#  open(OUT,">${phyla}_omclgns2.tab");

  open(F,">${phyla}_omclgns2.tab");
  open(OUT,">${phyla}_omclgn2sum.tab");
  my($lg,$log,$ng,$aev,$api,$mg,$mev,$mpi,$sev,$spi);
  
  sub dumpg{ if($ng>0){ $aev=sprintf"%.4g",$sev/$ng; $api=int($spi/$ng);
    print OUT join("\t",$lg,$log,$ng,$aev,$api,$mg,$mev,$mpi),"\n"; }
    $mpi=$mg=$mev=$lg=$ng=$sev=$spi=0; } 
  
  while(<F>) {
    my($g1,$g2,$og,$ev,$pi)=split; dumpg() unless($g1 eq $lg);
    $lg=$g1; $log=$og; $ng++; $sev+= $ev; $spi += $pi;
    if($pi>$mpi){ $mpi=$pi; $mg=$g2; $mev=$ev; } 
  }
  dumpg(); 
  close(OUT); close(F);

}


sub mclgrouptree
{


=item 7.	Create MCL tree joining of gene groups

  cd $orun
  $mc/bin/mclcm tmp/all_ortho.mtx -a "-I 2 --shadow=vl -te 2" > & log.mcm2 &
  #OR# 
  # $mc/bin/mclcm tmp/all_ortho.mtx -a "-I 3 --shadow=vl -te 2" > & log.mcm1 &
  $mc/bin/mcxdump -imx-tree mcl.cone --newick -o ${phyla}.newick -tab tmp/all_ortho.idx
  
  # cat ${phyla}_omclgn.tab ${phyla}.newick | perl -ne .. > ${phyla}.newick4
  set newk=${phyla}.newicki2
  
  ## improve these perls, mcltreesplit.pl
  cat ${phyla}_omclgn.tab  ${phyla}_omclgn.consensus_def.txt  ${newk} |\
  env gtag=ARP perl -ne 'BEGIN{$gtag=$ENV{gtag};} \
  $p=1; if(/^$gtag\d\D+(\d+)/){ $g=$gtag.$1; s/\(LOC\d+\)//;  s/src=\S+//; \
  $gd{$g}= (m/\s(\S.+)$/)? $1.";" : "";  $p=0; }  \
  elsif(/^($gtag\d+)\s+(\S+)$/){$ag{$2}= $1; $p=0;} \
  elsif(m/\(([^\)]+)\)/) { $gs=$1; @g=split",",$gs; %ga=(); $de="";  \
  map{ $ag=$ag{$_}; $ga{$ag}++ if($ag); } @g; \
  map{ $de .= $gd{$_} } sort keys %ga; $ga= join ",", sort keys %ga;\
  s/\([^\)]+\)/\($ga\)/; s/$/ # $de/ if $de; } print if $p; '\
   > ${newk}b
  
  cat ${newk}b | env gtag=ARP perl -ne\
  '$s=$_; s/\s*\#.*$//; s/[\,\s]+$//; $s=~s/ARP/ARP9_G/g; print $s unless($_ eq $ll); $ll=$_;' \
   > ${newk}c
   
  cat ${newk}c | env tag=ARP  perl ../omclw/mcltreesplit.pl > ${phyla}.clusters.tab
  
  #   add to genes.ugp.xml
  cat  ${phyla}.clusters.tab  ${phyla}_genes.ugp.xml |\
  perl -ne'if(/^ARC1_/){ chomp; ($c,$gc)=split"\t"; while( $gc =~ m/(ARP\w+)/g) { $gc{$1}= $gc; } } \
  elsif(/^ARDE_/){next;} else { if(m/<GeneSummary/) { $gc= (m/id=.\w+:(\w+)/) ? $gc{$1} : ""; } \
  elsif( $gc and m,</GeneSummary>,){ print "<related_gene_groups>\n$gc\n</related_gene_groups>\n"; }  \
  print; }'\
  >  ${phyla}_genes.ugp.xml2


=cut

}