#!/usr/bin/env perl
# blastab2sam.pl

=item blastab2sam.pl usage

 gunzip -c $dmag/genome/reads454f/asmnew/nwbdmag24g7d_asm-dmag7finall9b.mblastn.gz | grep -v '^#' |  \
  blastab2sam.pl -nosamechr -joins=1 -minid=99.9 -minaln=120  -minlow=0.10 \
  -onegenome 1 -no2ndary -debug \
  -outformat sam|tab -fasize nwbdmag24g7d_asm.facount - \
 > nwbdmag24g7d_asm-dmag7mrnabl_jp2i99a120nosame.sam

  from evigene/scripts/makeblastscore3.pl
  d.g.gilbert, 2015.11

=item possible bug excluding split-gene blasthit collection

  .. need to check mouse ref test set vs NCBI mouse split-gene list (~180)
  .. appear to be missing ~100 of those in this output, though small exam sez the blast hits
     show contig-split genes.
     
     
     
=item mods for gene-link lachesis scaffolding  

  -- current lachesis inputs from this:
    assembly.fasta@   assembly.fasta.names@  == input contigs.fa selecting only genelinks contig ids         
    assembly.fasta.counts_GENENEG.txt@    == putScafhits() table of ngenes/contig
    assembly_genelinks.sam  == putPairend() genelink.sam             
    
  ** NOTE: assembly.fasta.names ID order must be same as assembly_genelinks.sam.header,
    due to crude c++ string handling methods, current lachesis.c++ uses only input order int index of both,
    assumes same order.
    
  -- write full set of inputs to lachesis_evg ?

  -- cat outname.sam.header outname.sam > outname-all.sam for input to lachesis
  -- contig.fasta subset contained in links.sam
    perl -ne 'if(/^\@SQ\tSN:(\S+)/){$ok{$1}=1;} else{ if(/^>(\S+)/){$ok=$ok{$1};} print if $ok;} ' \
      $outname.sam.header $input.fasta > $outname.fasta
    grep '^>' $outname.fasta | sed 's/>//;' > $outname.fasta.names
    
  -- write lachesis.ini with proper outname values?
        
=item test outs

    /bio/bio-grid/daphmag/genome/reads454f/asmnew ; lachesis dmagf/
  # output=nwbdmag24g7d_asm_dmag7mrna.sam, nputgene=43596, nsplitgene=5644, ndupgene=11
    /bio/bio-grid/kfish2/genome/scafevg/drafta
  # output=funhe302scaf_kf2mrna_lascaf1a.sam, nputgene=41160, nsplitgene=3902, ndupgene=45

=cut

use strict;
use warnings;
use Getopt::Long ;

my $NJOIN= $ENV{joins}||1; # pair join step, ie adjacent, 2-off, 3-off
my $OVSLOP= 19; # need opt; FIXME: need overlapslop ~ < 0.01 of max part span ; or < 0.02..0.05 of min part span?
my $pctOVSLOP=$ENV{pctover} || 0.02; # was 0.06, too high; need opt
my $pMINLOW= $ENV{pmin} || 0.10; # want this low, valid exon align parts vary in size which this filters
my $MINALN=  $ENV{minalign} || 99; # align bases
my $pMINIDENT= $ENV{minident} || 95; # pctident
my $ONEGENOME=$ENV{oneref}||0; # 2015.02 ?? best spans ignoring scaffold splits .. require showSPAN ?
my $INTRONSIZE= $ENV{intron}||90; # guess, need?
my $SAMECHR=1;
# my $READSIZE=100; # size to cut from ends of long hsp aligns

my($ok, $fasta, $blastin, $outname, $no2ndary, $skiptypes, $outformat, $debug, $sizetab)= (0)x10;
my $dropdit=1; # dang id mess \.1 at end on some
$outformat="sam";

my $optok= &GetOptions (
  "blastin=s"=>\$blastin,  
  "fasize=s"=>\$sizetab, # facount ?
  "outname=s" => \$outname, # for .fasta instead of .sam, but location added for splitting.
  "outformat=s" => \$outformat, # for .fasta instead of .sam, but location added for splitting.
  "pairstep|joins=i"=>\$NJOIN,
  "MINALN|minalign=i"=>\$MINALN,  
  "MINIDENT|pMINIDENT=s"=>\$pMINIDENT,  
  "MINLOW|pMINLOW=s"=>\$pMINLOW,  
  "overlap=s"=>\$OVSLOP,  
  "pctOVSLOP=s"=>\$pctOVSLOP,  
  "oneref|ONEGENOME=s"=>\$ONEGENOME,  
  "samechr!"=>\$SAMECHR, 
  "no2ndary!"=>\$no2ndary, 
  "debug!"=>\$debug, 
  );

$blastin= shift @ARGV unless($blastin);

die "USAGE: $0 -blastin my.tabular.blastn > my.sam\n"
  unless($optok and $blastin);

## sam line; dont need as global now?
# my($qid, $flag, $chr, $cloc, $mapq, $cigar, $matechr, $mateloc, $isize, $seq, $qual, @opt)=(0) x 11;
#   $mapq=255;
#   $isize=0; #? pair insert size
#   $seq="*"; $qual="*"; @opt=();
#   $matechr="*"; $mateloc=0; # unless have mate info.. see paired below

sub MAINstub {}

##  $haveqlen now is count of size ids
my($haveqlen,$sizeh,$idorder,$trsizeh,$cdspanh)= readSizes($sizetab); 

my $inh = *STDIN;
if($blastin and not($blastin =~ /stdin|^-/)) { open($inh,$blastin) or die "reading $blastin"; }

my $outh= *STDOUT;
if($outname) { ## append .$outformat ??
  $outh=undef; open($outh, ">$outname") or die "ERR: writing $outname";
}

## sub putSamhead() .. FIXME: out.sam needs @SQ header lines == scafid\tscafsize
## @SQ	SN:scaffold_1	LN:3355761

# 1mary hit filters
$pctOVSLOP=$pctOVSLOP/100 if($pctOVSLOP >= 1);
$pMINLOW=$pMINLOW/100 if($pMINLOW >= 1);
$pMINIDENT=$pMINIDENT*100 if($pMINIDENT<1);

# 2ndary hit filters, lower than 1st/top hit per query
my $MINALN2   = _min(80,$MINALN);
my $pMINIDENT2= _min(95, 0.95*$pMINIDENT); # really depends on data, option?
if($no2ndary) { $MINALN2= $MINALN; $pMINIDENT2=$pMINIDENT; }

# NOT USED yet, opt for lower qual aligns that fill missing links?
my $MINALN3   = _min(90,$MINALN);
my $pMINIDENT3= _min(95,$pMINIDENT); # really depends on data, option?

my(%didid, %didscaf, %scafhit, %bspans);
my($lq, $bmax, $bmaxcut, $blerr, $nputout, $ndupgene, $nsplitgene)= (0) x 9; 
while(<$inh>) {
  next unless(/^\w/); chomp; my @v=split"\t"; 
  unless(@v==12){ warn"ERR: blasttab not 12 cols:'@v'\n"; $blerr++; die if ($blerr>9); next; }
  
  my($q,$t,$bits,$pid,$aln,$mis,$indl,@bspan)=  @v[0,1,11,2,3,4,5, 6,7,8,9]; # 6-9 =  q. start, q. end, s. start, s. end, 

  #all# my($q,$t,$pid,$aln,$mis,$indl,$rb,$re,$tb,$te,$ev,$bits)=@v; # blast table columns, outfmt=6/7
  ## expect q = geneid, t= scaffold/contig/genome ref id
  # @bspan = ($qb,$qe,$sb,$se);

  #NOW BELOW# next if($aln < $MINALN or $pid < $pMINIDENT); ## filter tiny aligns
  #^ ADJUST: strong filter for 1st q.t hsp, weaker for 2nds to fill in full qspan
  
  # if($swapqt) { ($q,$t)= ($t,$q); @bspan= @bspan[2,3,0,1];}
  if($q ne $lq) {
    $nputout+= putPairends($lq) if($lq);  
    %bspans=(); $bmaxcut= $bmax=0;
  }
  
  $bits= bint($bits);
  if($q ne $t) { #? assume $q ne $t
    my $aident= _max(0,$aln - $mis);  # other way to calc, maybe better? $pctident * $aln
    my $havespan=($ONEGENOME)?1:exists($bspans{$t});
    my $okfilter=0;
    my $pid2= ($aln >= 0)?(100*$aident/$aln):0; # overall maybe not an improvement..
    $okfilter= ($aln >= $MINALN and $pid2 >= $pMINIDENT and $bits >= $bmaxcut)?1:0; # filter 1st
    #x#$okfilter= ($aln >= $MINALN and $pid >= $pMINIDENT and $bits >= $bmaxcut)?1:0; # filter 1st
    #^^ problems w/ this pMINIDENT filter, knocking out top bitscore/eval, longest aligns
    # ignore pid for $pid2=100*$aident/$aln ?
    if( $havespan and not $okfilter ) { # filter 2nd
      $okfilter= ($aln >= $MINALN2 and $pid2 >= $pMINIDENT2)?1:0;
      #x#$okfilter= ($aln >= $MINALN2 and $pid >= $pMINIDENT2)?1:0;
    }
    if($okfilter) {
      sumscore( $q, $t, $bits,$aln,$aident, @bspan); ## if(($bits >= $bmaxcut) or exists($bspans{$t})); 
      if($bits > $bmax) { $bmax= $bits; $bmaxcut= $pMINLOW * $bmax; }
    }
  }
  
  $lq= $q; 
}
$nputout+= putPairends($lq) if($lq);
close($outh) if($outh);

#finiOutput()..
## ADD? 2nd rate hits of missed links? for not SAMECHR, $okfilter2=($aln >= $MINALN3 and $pid >= $pMINIDENT3);

putSamhead($outname) if($outformat=~/sam/);
putScafhits($outname) if($outformat=~/sam/); # samlachesis ??
warn "# output=$outname, nrows=$nputout, nsplitgene=$nsplitgene, ndupgene=$ndupgene\n"; # if($debug);
## output=nwbdmag24g7d_asm_dmag7mrna.sam, nputgene=43596, nsplitgene=5644, ndupgene=11
## .. nrows contain only the 5644 nsplitgene, nsp * 2pair * joins/gene
#------ subs -------------------


sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }
sub bint { my $b=shift; return ($b<0)?0 : ($b=~/e\+/)? int($b) : $b; }

## FIXME: out.sam needs @SQ header lines == scafid\tscafsize
## @SQ  SN:scaffold_1   LN:3355761
# perl -ne '($d,$c)=split; print "@SQ\tSN:$d\tLN:$c\n" if(/^\w/ and $d ne "total");' \
#  nwbdmag24g7d_asm.facount > nwbdmag24g 7d_asm.samhd
# cat asm.samhd asm.geneblast.sam > asm.genelinks.sam

sub putScafhits {
  my($outname)=@_;
  # simple table: scafid  genecount  genelinkcount for lachesis replacement RE-table
  # "GENENEG" is RE-motif, palindrome, desired by lachesis.
  my $nsc=0;
  my %idord=(); if(ref $idorder){ my $i=0; map{ $idord{$_}= ++$i; } @$idorder;}
  $outname ||="output.sam";
  $outname =~ s/\.sam/.fasta/; # lachesis now wants name: assembly.fasta.counts_RE.txt
  $outname.=".counts_GENENEG.txt" unless($outname=~/\.counts/);
  open(OUT,'>',$outname) or warn "#ERR writing $outname";
  for my $sc (sort{ $idord{$a} <=> $idord{$b} or $a cmp $b} keys %didscaf) {
    my $nlink= $didscaf{$sc}||0;
    my $nhit= $scafhit{$sc}||0;
    print OUT join("\t",$sc,$nhit,$nlink)."\n";
  }
  close(OUT); 
  return ($nsc,$outname);
}

sub putSamhead {
  my($outname)=@_;
  return 0 unless($haveqlen); # and outformat =~ /sam/
  my $nsc=0;
  ## dgg: DAMN really obscure data bug n Lachesis: SAM header contig IDS must be same order as contig.fasta IDs
  ## save input order from readSizes() ?? == @$idorder
  my %idord=(); if(ref $idorder){ my $i=0; map{ $idord{$_}= ++$i; } @$idorder;}
  # my $outf=($outname)?"$outname.header":"output.sam.header";
  $outname||="output.sam";
  $outname.=".header" unless($outname=~/\.header/);
  open(OUT,'>',$outname) or warn "#ERR writing $outname";
  for my $sc (sort{ $idord{$a} <=> $idord{$b} or $a cmp $b} keys %didscaf) {
     my $len=$sizeh->{$sc}||0; $nsc++;
     print OUT join("\t","\@SQ","SN:$sc","LN:$len\n");
  }
  #?? maybe combine both? cat output.sam.header output.sam > output-all.sam, unlink header?
  close(OUT); 
  return ($nsc,$outname);
}

#samflag  1=pair, 2=pairok, 4=mismatch 8=mismate, 16=rev, 32=revmate, 64=firstmate, 128=secondmate
use constant {
  samPAIR=> 0x0001, samPAIROK => 0x0002,
  samMISALN => 0x0004, samMISMATE => 0x0008, #  samMISALN f_mismatch == $chr eq "*"
  samREV => 0x0010, samREVMATE => 0x0020,
  samFIRST => 0x0040, samSECOND => 0x0080, # mate ids
  samSECALN=> 0x0100, samDUPRD => 0x0400,
  
  samMAPQ => 255, samSEQ => '*', samQUAL => '*',
};

sub pctScafEnd {
  my($chr,$cloc,$cend)=@_;
  # pctScafEnd: location as percent of total span, to select chr end points
  my $pTINYSCAF=0.40; # gene align span >= pTINY * chrspan
  if($haveqlen) { # add flag for picking scaf-end aligns
    my $csize= $sizeh->{$chr}||0; 
    my $qsize= 1+abs($cend-$cloc); #? or use gene span?
    my $cploc= ($csize>0)? int(0.5 + 100*$cloc/$csize):50;
    my $celoc= ($csize>0)? int(0.5 + 100*$cend/$csize):50;
    $cploc=$celoc if(abs($celoc - 50) > abs($cploc - 50));
    ## FIXME: size matters, if csize is small, cploc always *near* end.
    if($qsize > $pTINYSCAF*$csize) { $cploc= ($cploc<=50)?1:99; }
    $cploc=99 if($cploc>99); $cploc=1 if($cploc<1);
    return($cploc);
  } else {
    return(50);
  }
}

sub putsam {
  return puttab(@_) unless($outformat=~/sam/); #??
  my($qid,$flag,$chr,$cloc,$cend,$orient, $matechr,$mateloc, $insize, $aln,$mis,@xopt)=@_;  
  map{ $_||=0; } ($flag,$cloc,$cend,$mateloc,$insize,$aln,$mis);
  my $cigar= $aln."M"; # add indels, introns ?
  $matechr="=" if($matechr eq $chr);
  my @opt=();
  push(@opt,"NM:i:".$mis); # always?
  push(@opt,"XS:A:".(($orient<0)?"-":"+")); # strand tag
  push(@opt,"XE:i:".$cend); # keep end point of align : should also reconstruct cigar from exons
  push(@opt,"XP:i:".pctScafEnd($chr,$cloc,$cend)) if($haveqlen); # flag for picking scaf-end aligns
  push(@opt,@xopt) if(@xopt);
  ##global constants: $mapq=255; $seq="*"; $qual="*"; # dont need  
  print $outh join("\t",($qid, $flag, $chr, $cloc, samMAPQ, $cigar, $matechr, $mateloc, $insize, '*', '*', @opt)),"\n";
  return 1;
}
# x.8244810       16      Scaffold2       867     255     36M     *       0       0       CTAGACACCAAAAAAATAGCAGCGGCTATATGATTT      *
# x.1350890       0       Scaffold2       868     255     15M1I20M        *       0       0       TAGACACCAAAAAAAATAGCAGCGGCTATATGATTT      *

sub puttab {
  my($qid,$flag,$chr,$cloc,$cend,$orient, $matechr,$mateloc, $insize, $aln,$mis,@xopt)=@_;  
  map{ $_||=0; } ($flag,$cloc,$cend,$mateloc,$insize,$aln,$mis);
  ## flag useless here? == REV/FWD == orient?
  ## $qid.= ($flag && samREV)?"r":"f"; OR $cor.=($flag && samREV)?"r":"f"; as "+r", "+f", .. ?
  ## is $insize useful?   my $insert= abs($nxe - $lxb) + $INTRONSIZE;
  ## xopt[0] == "XX:A:$lxb-$lxe" useful here
  $matechr="=" if($matechr eq $chr);
  my($tspan)= $xopt[0] =~ m/XX:A:(\S+)/;
  my $pscafend= pctScafEnd($chr,$cloc,$cend); ## if($haveqlen)
  my $cor= ($orient<0)?"-":"+";
  my $aident= $aln - $mis;
  print $outh join("\t",$qid, $tspan, "$aident/$aln", $chr, $pscafend."%", $cloc, $cend, $cor,  $matechr, $mateloc),"\n";
  return 1;
}

sub putpair {
  my($geneid,$lhsp,$hsp)=@_;   
  # hsp == [$xb,$xe,$tb,$te,$xbit,$aln,$aident,$or,$tid]; 
  # should be lhsp.xb < hsp.xb ; or diff tid
  my $SAM_REVMATE= ($outformat=~/sam/)?1:0; #? drop this? lachesis doesn't need, esp swap locs
  
  my($lxb,$lxe,$ltb,$lte,$lxbit,$laln,$laident,$lor,$ltid)=@$lhsp;
  my($nxb,$nxe,$ntb,$nte,$nxbit,$naln,$naident,$nor,$ntid)=@$hsp;
  
  $scafhit{$ltid}++; $scafhit{$ntid}++; # add output table for lachesis pseudo-RE-count tab, include SAMECHR hits
  return 0 unless($SAMECHR or ($ltid ne $ntid));
  # if($ONEGENOME) .. ltid / ntid are true tid, ok here
  
  ## FIXME check blast locs not outabounds for chr sizes.. due to other error, split contigs..
  my $lerr=0;
  if($haveqlen) {
    # REQUIRE each scaf id to have size? maybe yes, ie filter out unwanted scafs from data
    my $len=$sizeh->{$ltid}||0; return 0 unless($len);
    map{ if($_<1) { $lerr++; $_=1; } elsif($_>$len) { $lerr++; $_=$len; } } ($ltb,$lte);
    $len=$sizeh->{$ntid}||0;  return 0 unless($len);
    map{ if($_<1) { $lerr++; $_=1; } elsif($_>$len) { $lerr++; $_=$len; } } ($ntb,$nte);
  }

  my $lmis= _max(0,$laln - $laident); # revert laident to mismatch count??
  my $nmis= _max(0,$naln - $naident); # revert laident to mismatch count??
  # also, if aln != span(xe-xb) != span(te-tb), have indels, add to cigar?

  my $insert= abs($nxe - $lxb) + $INTRONSIZE;
  my $gpid=$geneid; for(my $p=0; $p<9999; $p++) { $gpid="$geneid.$p"; last unless($didid{$gpid}); }
  $didid{$gpid}++;
  $didscaf{$ltid}++; $didscaf{$ntid}++;
 
  my($nput,$flag1,$flag2)=(0) x 9;
  $flag1= samPAIR + samPAIROK + samFIRST;  
  $flag2= samPAIR + samPAIROK + samSECOND;
  if($SAM_REVMATE) { ## is this right for rev-mate? change before 1st putsam() ?
    ($nte,$ntb)=($ntb,$nte); $nor=($nor < 0)?0:-1;
    $flag1+= samREVMATE; # flop rev/revmate using tor
    $flag2+= samREV;
  }

  $nput+= putsam($gpid,$flag1,$ltid,$ltb,$lte,$lor,$ntid,$ntb,$insert,$laln,$lmis,"XX:A:$lxb-$lxe");
  # need to rev before 1st put, not here .. screws up mate loc check
  $nput+= putsam($gpid,$flag2,$ntid,$ntb,$nte,$nor,$ltid,$ltb,-$insert,$naln,$nmis,"XX:A:$nxb-$nxe");
  return $nput;
}

sub _sortQlocAlign {
  # $a,$b == [$xb,$xe,$tb,$te,$xbit,$aln,$aident,$or,$ttrue1]; 
  # sort on gene start: $xb > xbit/aln/aident
  return ($a->[0] <=> $b->[0] or $b->[4] <=> $a->[4] or $a->[8] cmp $b->[8]);
}


  
=item getAllHsp() methods  

  -- redo sumscore() with ONEGENOME set to pick out split spans ..
     .. save %bspan, ONEG flag, new sumscore( $q, $t, $bits,$aln, $aident, $qb,$qe,$sb,$se);

  -- maybe add dup hsp filter, esp for daphnia 
      .. exclude genes or hsps that have near equal exon aligns to diff scaf
      .. need care in 'near equal' cut-offs, so that proper uniq aligns are not removed.
      .. ie., check both identical exon hsps, and total aident/tid,
         remove when have 2+ ident on both scores, or else lowish overall score for 2+ (ambiguous)
=cut
    
sub getAllHsp {
  my($qgeneid)= @_;
  my @allhsp=();
  my @tid= sort keys %bspans;
  if($ONEGENOME or @tid<2) {
    #? screen for isdupgene?
    foreach my $tid (@tid) { push @allhsp, @{$bspans{$tid}}; } 
  } else {
    my (%tspan);
    foreach my $tid (@tid) { 
      my @bsp= @{$bspans{$tid}}; 
      my($tsp,$taid,$tdups)=(0) x 9; my($tb,$te)= @{$bsp[0]};
      for my $bsp (@bsp) { 
        ## duphsp from sumscore() not right yet.. need @allhsp exons test
        my($xb,$xe,$sb,$se,$xbit,$xaln,$xaid,$xor,$xtrueid,$xduphsp)= @$bsp; 
        $tsp += 1 + $xe - $xb; $taid += $xaid; $tdups+= $xduphsp;
        $tb= _min($tb,$xb); $te= _max($te,$xe); 
        }
      $tspan{$tid}=[$tsp,$tb,$te,$taid,$tdups];
    }     
    my($topid,@xtid)= sort{ $tspan{$b}->[0]<=>$tspan{$a}->[0] } keys %tspan; #tspan sort, try taid?
    my($tsp,$tb,$te,$taid,$tdups)= @{$tspan{$topid}};
    my($issplit,$isdupgene)=(0) x 9;
    for my $xid (@xtid) {  
      my($xsp,$xb,$xe,$xaid,$xdups)= @{$tspan{$xid}}; 
      #y# if(($xdups>0 or $tdups>0) and ($xaid >= 0.985*$taid)) { $isdupgene++; } # need dupfilter option?
      #x# if($xdups>0 and $tdups>0 and ($xaid >= 0.99*$taid)) { $isdupgene++; } # need dupfilter option?
      # ^^ sumscore() may have already filtered out xdups? check only tdups? xdups or tdups ?
      if($xb>= $te-$OVSLOP or $xe <= $tb+$OVSLOP) { $issplit++; } # is this robust enough? dont care about inside splits?
    }
    
    if($isdupgene) { # skip ?
      $ndupgene++;
      # 1st try w/ daphmag genes, only 1 dupgene
      # 2nd try, xdups OR tdups, 11 dupgene .. not as much as expected.
      # output=nwbdmag24g7d_asm_dmag7mrna.sam, nputgene=43596, nsplitgene=5644, ndupgene=11

    } elsif(not $issplit) {
      push @allhsp, @{$bspans{$topid}}; # only best align, this winnows out many 2ndary aligns
    } else {
      my $SONEG=$ONEGENOME; $ONEGENOME=1; # flag sumscore() to add only gene-distinct hsps over all scafs
      my %bspansave= %bspans; %bspans=();
      foreach my $tid ($topid,@xtid) { 
        for my $bsp (@{$bspansave{$tid}}) { 
          my($xb,$xe,$sb,$se,$xbit,$xaln,$xaid,$xor,$xtrueid,$xduphsp)= @$bsp;
          if($xor<0) { ($sb,$se)=($se,$sb); } # recover xor this way
          sumscore( $qgeneid, $xtrueid, $xbit, $xaln,$xaid,$xb,$xe,$sb,$se); 
        }
      }
      foreach my $tid (sort keys %bspans) { push @allhsp, @{$bspans{$tid}}; } 
      $nsplitgene++;
      $ONEGENOME= $SONEG; %bspans= %bspansave; # bspans dont care now?
    }
  }
  
  # return sort _sortQlocAlign @allhsp; # leave to caller
  return @allhsp;
}

sub putPairends {
  my($qgeneid)= @_;
    # pairends: collect all spans, sort by qstart,align,
    #   then pick out adjacent, off-2, off-3 paired aligns of q to all t,
    #   cut READSIZE ends of each align pair for output
    
  my $nput=0;

  # FIXME2: ONEGENOME making mistakes, drop bspans tid == onegenome, but add here 2nd filter
  #  a. look for max tid covering tspan, 
  #  b. if a. fails full tspan coverage, look for split-part tids covering tspan
  # @{$bspans{$t}}, [$qb,$qe,$sb,$se,$bits,$aln,$aident,$or,$ttrue]  
  
  my @allhsp= getAllHsp($qgeneid);
  return($nput) unless(@allhsp);
  #old# foreach my $tid (sort keys %bspans) { push @allhsp, @{$bspans{$tid}}; }
  
#  ## UGH* Some of these are duplicate gene mappings, now output as cross-scaf pairs at diff gene exons..
#  eg. both scaf locs have identical aligns to same gene.. no surpise for daph.
#  ans is to add back ONEGENOME methods of makeblastscore3 ..
#   my $ONEGENOME=$ENV{oneref}||0; # 2015.02 ?? best spans ignoring scaffold splits .. require showSPAN ?
#  sub sumscore() keeps both due to identical scores.
# grep '^Dapma7bEVm029002t1 '  $dmag/genome/reads454f/asmnew/nwbdmag24g7d_asm-dmag7finall9b.tall6     
# Dapma7bEVm029002t1	dmag24nwb7d_scaffold00020	1462	791	791	798	1101401	12-797/12-158,154-797	709447-709593,709685-710328:+
# Dapma7bEVm029002t1	dmag24nwb7d_scaffold00004	1462	791	791	798	2396214	12-797/12-158,154-797	1867808-1867954,1867073-1867716:-
  
  ## change this algo, count to NJOIN only if putpair()
  use constant bNjoinPuts => 1;
  
  my @sexons= sort _sortQlocAlign @allhsp;
  my $nx= @sexons;
  
  my @lhsp=();
  my $MAXLOOK=9;  
  my $MINX= 100; # option
  use constant bNjoinPuts => 1;
  use constant bPickJunctions => 0; # maybe not as good
  
  ## better? method, from genelinkgff2sam.pl, just now for $NJOIN == 1
  if(bPickJunctions and not $SAMECHR) { # pick out join junctions, nearest exons
    # this way gets same ngenes, same nchr but fewer rows of sam links .. want more per gene x chr ?
    # note this screws scafhit count as no SAMECHR counted here
    my @junc;
    for(my $i=1; $i<$nx; $i++) {
      my $j=$i-1;
      
      ## should make these to same order records
      # my($lchr,$lgb,$lge,$lor,$lxb,$lxe,$lgid)=@{$sexons[$j]}; #gff
      # my($nchr,$ngb,$nge,$nor,$nxb,$nxe,$ngid)=@{$sexons[$i]}; #gff
      my($lxb,$lxe, $lgb,$lge,$lxbit,$laln,$laident,$lor,$lchr)=@{$sexons[$j]};  #blast hsp
      my($nxb,$nxe, $ngb,$nge,$nxbit,$naln,$naident,$nor,$nchr)=@{$sexons[$i]};
      
      ## ok4case of multi-exons per join part, ok4?case of 1-exon in join part, both lx>one and one>nx
      if($lchr ne $nchr) {
        my $lw=$lxe-$lxb; my $nw=$nxe-$nxb;
        if($lw<$MINX) { $j-- if($j>0 and $sexons[$j-1]->[0] eq $lchr); }
        if($nw<$MINX) { $i++ if($i+1<$nx and $sexons[$i+1]->[0] eq $nchr); }
        push @junc, $sexons[$j], $sexons[$i];
      }
    }
    # FIXME: $NJOIN > 1 case
    my($j,$i,$nj)=(0,1,scalar(@junc));
    while($i<$nj) {
      my $pput= putpair($qgeneid,$junc[$j],$junc[$i]);
      $nput += $pput;
      $i++; $j++; # get both ends.
    }
    return($nput);
  }
  
if(bNjoinPuts) {
  for my $hsp (@sexons) {
    my $nj= _min($MAXLOOK,scalar(@lhsp));
    my $njdid=0;
    #o# for(my $j=$nj-1; $j>=0; $j--)  # bad order, j=0..nj instead
    for(my $j=0; $j<$nj; $j++)  # want last 2 hsp as 1st pair
    {
      my $pput= putpair($qgeneid,$lhsp[$j],$hsp);
      $nput+= $pput;
      $njdid++ if($pput>0);
      last if($njdid>=$NJOIN);
    }
    pop @lhsp if(@lhsp>$MAXLOOK);
    unshift @lhsp, $hsp;
  }
} else {  
  for my $hsp (@sexons) {
    my $nj= _min($NJOIN,scalar(@lhsp));
    for(my $j=$nj-1; $j>=0; $j--) {
      $nput += putpair($qgeneid,$lhsp[$j],$hsp);
    }
    pop @lhsp if(@lhsp>$NJOIN);
    unshift @lhsp, $hsp;
  }
}
  
  ## ADD? not here? unless($nput) { $savegene{$qgeneid}=\@allhsp; } # opt to save
  return($nput);
}
  

sub sumscore {
  my( $q, $t, $bits,$aln, $aident, $qb,$qe,$sb,$se) = @_;
  my $or=0;
  my $ttrue=$t; if($ONEGENOME) { $t=$ONEGENOME.$q; } # 2015.02 ?? best spans ignoring scaffold splits ?
  if($qb > $qe) { ($qb,$qe)= ($qe,$qb); $or--; } # FIXME record $or in bspans
  if($sb > $se) { ($sb,$se)= ($se,$sb); $or--; }

  unless($bspans{$t}) { 
    $bspans{$t}=[]; 
    $bspans{$ttrue}=[] if($ONEGENOME); # need for above test bspans{torig} for sumscore()
    push( @{$bspans{$t}}, [$qb,$qe,$sb,$se,$bits,$aln,$aident,$or,$ttrue,0]); 
    return; }
    
  my $ov=0; my $duphsp=0;
  my $qlen=1+$qe-$qb; my $slen=1+$se-$sb;
  my $qslop= _max($OVSLOP, int($pctOVSLOP*$qlen));
  my $sslop= _max($OVSLOP, int($pctOVSLOP*$slen));
    
  foreach my $sp (@{$bspans{$t}}) {
    my($xb,$xe,$tb,$te)= @$sp;  my $ttrue1= $$sp[8];
    # my($qbtrim,$qetrim)=(0,0);
    
    # query gene spans
    if($qe < $xb or $qb > $xe) { }
    elsif($qe > $xe and $qb >= $xe - $qslop) { 
      if($qb <= $xe) { $qb=$xe+1; } # FIXME 201511: trim out slop so no overlap for final spans
    }
    elsif($qb < $xb and $qe <= $xb + $qslop) { 
      if($qe >= $xb) { $qe=$xb-1; } # FIXME 201511: trim out slop so no overlap for final spans
    }
    else { 
      # ADD? mark this hsp/query-exon as dup/multimap, for qual filter?
      # is this ok? may need to check all for duphsp, after over test, 
      #** not useful here, need to look at all @target exons for overlaps **
      use constant pDUPEXON => 0.75;      
      #x# if( (_min($xe,$qe) - _max($xb,$qb)) > pDUPEXON*$qlen) { $sp->[9]++; $duphsp++; } 
      $ov=1; last; 
    }
    # target chr spans
    if($ONEGENOME and $ttrue ne $ttrue1) { } # 2015.02
    elsif($se < $tb or $sb > $te) { }
    elsif($se > $te and $sb >= $te - $sslop) { 
      if($sb <= $te) { $sb=$te+1; } # FIXME 201511: trim out slop so no overlap for final spans
    }
    elsif($sb < $tb and $se <= $tb + $sslop) { 
      if($se >= $tb) { $se=$tb-1; } # FIXME 201511: trim out slop so no overlap for final spans
    }
    else { $ov=2; last; }
  }  
  
  unless($ov) { 
    push( @{$bspans{$t}}, [$qb,$qe,$sb,$se,$bits,$aln,$aident,$or,$ttrue,$duphsp]); 
  }
}

sub readSizes {
  my(@inf)= @_; 
  my $nt=0;
  my (%alen,%trlen,%cdspan,%aaqual); %alen=(); 
  my $AAGAP=0;
  my $CDSSPAN=0;
  my $hasgap=($AAGAP)?1:0;#  && $iscount .. my $iscount=1; 
  my $hasspan=$CDSSPAN;#  collect %trlen,%cdspan ?
## DAMN bug lachesis needs same contig ID order for this ctg.fasta and sam.header
  my @idin=();
  return (0) unless(@inf and $inf[0]);
  foreach my $aaf (grep/\w/, map{ split(/[,\s]+/,$_) }@inf) {
    if(open(F,$aaf)) { my $n=0;
      while(<F>){ 
        my($id,$aw,@ac)=split; 
        if($hasgap) { if(@ac>1 or @ac==0){ $hasgap=0; } else { $aw -= $ac[0]; }} 
        if($hasspan) { if(@ac>3 and $ac[3]=~/\d/){ $aaqual{$id}=$ac[1]; $trlen{$id}=$ac[2]; $cdspan{$id}=$ac[3]; } else { $hasspan=0; }} 
        $alen{$id}=$aw; $n++;  push @idin, $id;
      } close(F); 
      $nt+=$n; warn  "# read n=$n from $aaf\n" if $debug;
    } else {
      warn "# cant read sizes from $aaf\n" if($aaf);# if $debug
    }
  }
  ## $haveqlen=1 if($nt>0);
  return($nt,\%alen,\@idin,\%trlen,\%cdspan,\%aaqual);
}
