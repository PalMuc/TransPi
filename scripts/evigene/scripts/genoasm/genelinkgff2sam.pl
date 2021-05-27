#!/usr/bin/env perl
# genelinkgff2sam.pl

=item about genelinkgff2sam.pl 

* FIXME bad choice of exon pairs at split, _exonSort problem..
* pick 2 nearest exons, UNLESS too tiny, then next nearest, at split points

  make sam align table for input to scaffolder software

  Input should be genes aligned, carefully, to fragmented chromosome assembly
  e.g. with ncbi splign and/or gmap, having splice signals, 
  and full gene paths that are often split over contigs/scaffolds.
  exon.gff is primary input, by preference these have "Target=geneid start end" attribs.
  Expect mRNA, exons to have same ID for split genes, and be ordered together. 
  Desirable to  have Split= and Target= attributes (or trg= that evigene uses in some gff, should revert)
  
  output normally has the paired exons of split genes to link contigs,
  with reliability scores.
  
  from evigene/scripts/genoasm/blastab2sam.pl,
  evigene/scripts/ests/gff2aligntab.pl  and ests/fagff2sam.pl
  d.g.gilbert, 2015.11

=item usage

  genelinkgff2sam.pl -nosamechr|splitonly -joins=1 \
  -outformat sam|tab -fasize nwbdmag24g7d_asm.facount genes.gff \
     > myname.sam

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

=item needs basic test w/ mouse ref
  
  using ncbi mouse ref genome + gene set, and contig asm derived
  from that w/ chr-contig.agp, convert genes.gff to contig locs and
  test run thru this + lachesis.  Appears that lach. clustering works right
  for mouse ref, but in-cluster ordering/orient has near 50% mistakes,
  from using genelink blast2sam data.  Input of this likely wont cure those
  order/orient mistakes.
        
=item test outs

    /bio/bio-grid/daphmag/genome/reads454f/asmnew ; lachesis dmagf/
  # output=nwbdmag24g7d_asm_dmag7mrna.sam, nputgene=43596, nsplitgene=5644, ndupgene=11
    /bio/bio-grid/kfish2/genome/scafevg/drafta
  # output=funhe302scaf_kf2mrna_lascaf1a.sam, nputgene=41160, nsplitgene=3902, ndupgene=45

=cut

use strict;
use warnings;
use Getopt::Long ;

my $NJOIN=$ENV{joins}||1; # pair join step, ie adjacent, 2-off, 3-off
my $MINALN= $ENV{minalign} || 99; # align bases
my $pMINIDENT= $ENV{minident} || 95; # pctident
my $READSIZE=100; # size to cut from ends of long hsp aligns
my $INTRONSIZE= $ENV{intron}||90; # guess, need?
my $SAMECHR=0;

my($ok, $fasta, $blastin, $outname, $no2ndary, $skiptypes, $outformat, $debug, $sizetab)= (0)x10;
my $dropdit=1; # dang id mess \.1 at end on some
$outformat="sam";

my $optok= &GetOptions (
  "gffin=s"=>\$blastin,  
  "fasize=s"=>\$sizetab, # facount ?
  "outname=s" => \$outname, # for .fasta instead of .sam, but location added for splitting.
  "outformat=s" => \$outformat, # for .fasta instead of .sam, but location added for splitting.
  "samechr!"=>\$SAMECHR, 
  "pairstep|joins=i"=>\$NJOIN,
  "debug!"=>\$debug, 
  );
  #"MINALN|minalign=i"=>\$MINALN,  
  #"MINIDENT|pMINIDENT=s"=>\$pMINIDENT,  
  #"MINLOW|pMINLOW=s"=>\$pMINLOW,  

$blastin= shift @ARGV unless($blastin);

die "USAGE: $0 -gff mygenes.gff -fasize contigsizes.table > my.sam\n"
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
# $pctOVSLOP=$pctOVSLOP/100 if($pctOVSLOP >= 1);
# $pMINLOW=$pMINLOW/100 if($pMINLOW >= 1);
# $pMINIDENT=$pMINIDENT*100 if($pMINIDENT<1);

# 2ndary hit filters, lower than 1st/top hit per query
# my $MINALN2   = _min(80,$MINALN);
# my $pMINIDENT2= _min(98,$pMINIDENT); # really depends on data, option?
# if($no2ndary) { $MINALN2= $MINALN; $pMINIDENT2=$pMINIDENT; }

## gff align attribs
my @ATK   = qw(match qlen cov pid path indels nexon splice aalen offs aamap sense Split oid tag); 
my @ATKGS = qw(gescore clen gaps chimera); # my @hd=grep{ not/Split/ } @ATK;
my @ATKOUT= qw(cov pid nexon qlen path Split);

my(@exon, %geneat, %didid, %didscaf, %scafhit, %bspans, %split, );
my($ltid,$tag,$blerr, $nputout, $ndupgene, $nsplitgene, $geneatval, $geneerr, $exonbegin)= (0) x 19; 

# sub readGenesGff {}
while(<$inh>) {
  next unless(/^\w/); chomp; my @v=split"\t"; 
  unless(@v==9){ warn"ERR: gff not 9 cols:'@v'\n"; die if (++$blerr>9); next; }
  
  my($chr,$src,$typ,$gb,$ge,$gsc,$gor,$xph,$at)=@v;

  next unless($typ =~ /^(mRNA|exon)/);
  # expect mRNA, exons to have same ID for split genes, and be ordered together

  my $kid=($typ=~/mRNA/)?"ID":"Parent"; my $needtbe=0;
  my($tid,$tb,$te)= $at=~m/(?:trg|Target)=(\S+) (\d+) (\d+)/; # may not be there
  my($id)= $at=~m/$kid=([^;\s]+)/;
  ## fixme: have OLD oids in Target for some, check vs ID= idtag ?
  ## also this problem: swapmain=Funhe2EKm005226t1,Funhe2EKm005226t12 << where Target=t1, ID=t12 
  ## always use ID|Parent= but Target tb,te
  unless($tid) { $tid=$id; $needtbe=1; ## $tb=1; $te=1+ $ge - $gb; 
    unless($tid) { warn"ERR: $typ missing ID/Parent for '@v'\n";  die if (++$blerr>9);  next; } 
  } 
  $tid=$id; ## now always
    # my($tp,$ip)=map{ m/^(\D+)/; $1; } ($tid,$id); # bad
    # my($tp,$ip)=map{ my $t=$_; $t=~s/\d+t\d+$//; $t; } ($tid,$id);
    # $tid=$id if($tp ne $ip); 

  if($ltid ne $tid) {
    $nputout+= putPairends($ltid,\@exon,\%geneat) if($ltid);
    @exon=(); %geneat=(); $geneerr=$geneatval=""; $exonbegin=0;
    ## proper calcs require gene rows be all together, including split parts, either with Target= locs
    ## or with exons in ascending order of mRNA position
  }
  
  my $issplit=(/(?:Split|chimera)=([^;\s]+)/)?$1:0;  # use at{Split} # 

  if($typ eq "mRNA") {
    ## @ATK=qw(match qlen cov pid path indels nexon splice aalen offs aamap sense Split oid tag);
    ## sam.attr useful: cov,pid,nexon,Split/path,
    my %at=(); map{ my $v=($at=~m/\b$_=([^;\s]+)/)?$1:0; $at{$_} = $v; } ("ID",@ATK);
    my %ags=(); map{ if($at=~m/\b$_=([^;\s]+)/) { $ags{$_}=$1; } } @ATKGS;
    my $id= $at{ID}; # may not == tid
    
    if($issplit and not $needtbe) { # check for messy target overlap, cancel bad case like below
      if($split{$tid} and (my @sp= @{$split{$tid}})) { 
        my $si=1 + @sp;
        push @{$split{$tid}}, [$tb,$te,$si];
        my $overs= overlap($split{$tid});
        if($overs) { 
         $geneerr="ERR: $tid has overlap split parts $issplit, ov:$overs;";
         $geneat{ERROR}.=$geneerr; $ndupgene++; 
         }
      } else { 
        push @{$split{$tid}}, [$tb,$te,1];
      }
    }
    
    if($needtbe) { $tb=$exonbegin; $te=$tb + ($ge - $gb); } # NOt for mrna: $exonbegin=$te+1;
    if(/gescore=/) { ## gsplign
      $tag="gspl" unless($tag);
      if($at{splice}) { $at{splice}= 2 + int($at{splice}/2); }
      my($pcov,$match,$ql2)= $at{cov} =~ m/(\d+).,(\d+).(\d+)/;
      $at{pid}||= 99; if($pcov and $match){ $at{cov} = $pcov; $at{match}= $match; }
      $at{qlen}= $ags{clen}||$ql2; # $gaps=$ags{gaps}||0; 
    } else { # gmap
      $tag="gmap" unless($tag);
      my $aamap= $at{aamap}||$at{aalen}||0;
      my($chi)= $id =~ m/_C(\d)$/; if($ags{chimera} or $chi) { $chi||=1; $at{Split}="$chi/2" unless($at{Split}); }
      if($at =~ /aalen=(\d+,\d+[^;\s]+)/) { my $aaq=$1; $at{aalen}=$aaq; $at{aamap}=$aamap if($aamap ne $aaq); }
    } 

    if(m/cxlen=(\d+).(\d+)/) { my($cw,$xw)=($1,$2); $at{qlen}=$xw; } #$at{cdslen}=$cw;
    # $at{nexon}.=".".$at{ncds}; # print only nexon col, .ncds as fraction
    # $at{splice}= 2+$inspl if($inspl>0 and not $at{splice});
    if(not $at{cov} and $at{match}>0 and $at{qlen}>0){ $at{cov}=int(0.5+100*$at{match}/$at{qlen}); }
    $at{tag}= $tag;
    map{ $geneat{$_}= $at{$_}; } keys %at; # %geneat= %at; # last Split part is recorded, want other?
    $geneat{GENEID}= $tid; #? ensure geneat{ID} == tid ?
    $geneatval= geneatval($tid, \%geneat);
    $nsplitgene++ if($issplit or $at{Split}); # counts gene parts
    
  } elsif($typ eq "exon") {
    if($needtbe) { $tb=$exonbegin; $te=$tb + ($ge - $gb); $exonbegin=$te+1;}
    push @exon, [$chr,$gb,$ge,$gor,$tb,$te,$tid];
    # blast hsp == [$xb,$xe,$tb,$te,$xbit,$aln,$aident,$or,$tid]; 
  }
  
  $ltid= $tid; 
}

$nputout+= putPairends($ltid,\@exon,\%geneat) if($ltid);
close($outh) if($outh);

putSamhead($outname) if($outformat=~/sam/);
putScafhits($outname) if($outformat=~/sam/); # samlachesis ??
warn "# output=$outname, nrows=$nputout, nsplitgene=$nsplitgene, ndupgene=$ndupgene\n"; # if($debug);
## output=nwbdmag24g7d_asm_dmag7mrna.sam, nputgene=43596, nsplitgene=5644, ndupgene=11
## .. nrows contain only the 5644 nsplitgene, nsp * 2pair * joins/gene


#------ subs -------------------

sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }
sub bint { my $b=shift; return ($b<0)?0 : ($b=~/e\+/)? int($b) : $b; }

sub geneatval{ #my($tid,$geneat)=@_;
  my @v= map{ ($geneat{$_}) ? "$_=$geneat{$_}" : "" } ('ERROR',@ATKOUT); 
  return join",", grep /=/, @v; 
}  

sub isover{ 
  my($lb,$le,$nb,$ne)=@_; my $OVSLOP=19; # 9?
  return 0 if($nb+$OVSLOP >= $le or $ne-$OVSLOP < $lb);
  if($nb<$lb and $ne < $le) { ($lb,$le,$nb,$ne)= ($nb,$ne,$lb,$le); } # sorted?
  return ($nb < $le and $ne > $lb)?1:0;
}

sub overlap { # detect bad split parts by target overlaps
  my($partspans)=@_;
  my $np=@$partspans; my $nov=0;
  for(my $i=1; $i<$np; $i++) {
    my($lb,$le,$li)= @{$partspans->[$i-1]};
    my($nb,$ne,$ni)= @{$partspans->[$i]};
    my $ov= isover($lb,$le,$nb,$ne); $nov++ if($ov); # return 1 if($ov);
  }
  return $nov;
}

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
    my $nhit = $scafhit{$sc}||0;
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
  push(@opt,"XS:A:".$orient); # blast: (($orient<0)?"-":"+")); # strand tag
  push(@opt,"XE:i:".$cend); # keep end point of align : should also reconstruct cigar from exons
  push(@opt,"XP:i:".pctScafEnd($chr,$cloc,$cend)) if($haveqlen); # flag for picking scaf-end aligns
  push(@opt,@xopt) if(@xopt); # now holds GA:A:$geneatval
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
  my $gat= $geneatval||"noat"; # or from @xopt ??
  my $pscafend= pctScafEnd($chr,$cloc,$cend); ## if($haveqlen)
  my $cor= $orient; # blast: ($orient<0)?"-":"+";
  my $aident= $aln - $mis;
  print $outh join("\t",$qid, $tspan, "$aident/$aln", $chr, $pscafend."%", $cloc, $cend, $cor,  $matechr, $mateloc,$gat),"\n";
  return 1;
}

sub putpair {
  my($geneid,$lhsp,$hsp)=@_;   
  # exon == [$chr,$gb,$ge,$gor,$xb,$xe,$gid]; # gff exon, NOTE chr == tid, target id
  ### NOTE changed exon/hsp field order, exon target tb/te == xb,xe in hsp, gb/ge = tb/te in hsp
  # hsp  == [$xb,$xe,$tb,$te,$xbit,$aln,$aident,$or,$tid]; blast hsp
  # should be lhsp.tb < hsp.tb ; or diff tid
  my $SAM_REVMATE= 0; # drop rev ($outformat=~/sam/)?1:0; #? drop this? lachesis doesn't need, esp swap locs
  
  # my($lxb,$lxe,$ltb,$lte,$lxbit,$laln,$laident,$lor,$ltid)=@$lhsp;
  # my($nxb,$nxe,$ntb,$nte,$nxbit,$naln,$naident,$nor,$ntid)=@$hsp;
  my($ltid,$ltb,$lte,$lor,$lxb,$lxe,$lgid)=@$lhsp;
  my($ntid,$ntb,$nte,$nor,$nxb,$nxe,$ngid)=@$hsp;
  
  $scafhit{$ltid}++; $scafhit{$ntid}++; # add output table for lachesis pseudo-RE-count tab, include SAMECHR hits
  return 0 unless($SAMECHR or ($ltid ne $ntid));
  # $nsplitgene++; # exon pair, not gene..

  my $lerr=0;
  if($haveqlen) {
    # REQUIRE each scaf id to have size? maybe yes, ie filter out unwanted scafs from data
    my $len=$sizeh->{$ltid}||0; return 0 unless($len);
    map{ if($_<1) { $lerr++; $_=1; } elsif($_>$len) { $lerr++; $_=$len; } } ($ltb,$lte);
    $len=$sizeh->{$ntid}||0;  return 0 unless($len);
    map{ if($_<1) { $lerr++; $_=1; } elsif($_>$len) { $lerr++; $_=$len; } } ($ntb,$nte);
  }

  # my $lmis= _max(0,$laln - $laident); # revert laident to mismatch count??
  # my $nmis= _max(0,$naln - $naident); # revert laident to mismatch count??
  my($laln,$lmis)=(1+$lxe-$lxb,0); # fake it for now, use gmap quals
  my($naln,$nmis)=(1+$nxe-$nxb,0); # fake it for now, use gmap quals
  
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

  $nput+= putsam($gpid,$flag1,$ltid,$ltb,$lte,$lor,$ntid,$ntb,$insert,$laln,$lmis,"XX:A:$lxb-$lxe","GA:A:$geneatval");
  # need to rev before 1st put, not here .. screws up mate loc check
  $nput+= putsam($gpid,$flag2,$ntid,$ntb,$nte,$nor,$ltid,$ltb,-$insert,$naln,$nmis,"XX:A:$nxb-$nxe","GA:A:$geneatval");
  return $nput;
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

=item old blastin  
  
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
      my $SONEG=$ONEGENOME; $ONEGENOME=1;
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

=cut

sub _sortQlocAlign {
  # $a,$b == [$xb,$xe,$tb,$te,$xbit,$aln,$aident,$or,$ttrue1]; 
  # sort on gene start: $xb > xbit/aln/aident
  return ($a->[0] <=> $b->[0] or $b->[4] <=> $a->[4] or $a->[8] cmp $b->[8]);
}

sub _sortExons {
  ## @exon   [$chr,$gb,$ge,$gor,$tb,$te,$tid];
  # sort on gene start for one gene tid: $tb, independent of genome locations chr,gb  
  return ($a->[6] cmp $b->[6] or $a->[4] <=> $b->[4] or $b->[5] <=> $a->[5]);
}

sub putPairends {
  my($qgeneid,$exons,$geneat)= @_;
    # pairends: collect all spans, sort by qstart,align,
    #   then pick out adjacent, off-2, off-3 paired aligns of q to all t,
    #   cut READSIZE ends of each align pair for output
    
  my $nput=0;
  return($nput) unless(ref $exons and @$exons>1);
  if(ref $geneat and $geneat->{ERROR}) { 
    warn $geneat->{ERROR}."\n" if($debug);
    return(0); } # record problem cases?
  
# * FIXME bad choice of exon pairs at split, _exonSort problem..
# * pick 2 nearest exons, UNLESS too tiny, then next nearest, at split points
  
  my @sexons = sort _sortExons @$exons; # ascending gene parts, ignoring chr starts
  my $nx= @sexons;
  
  my @lhsp=();
  my $MAXLOOK=9;  
  my $MINX= 100; # option
  use constant bPickJunctions => 1;
  
  if(bPickJunctions and not $SAMECHR) { # pick out join junctions, nearest exons
    # this way gets same ngenes, same nchr but fewer rows of sam links .. want more per gene x chr ?
    # note this screws scafhit count as no SAMECHR counted here
    my @junc;
    for(my $i=1; $i<$nx; $i++) {
      my $j=$i-1;
      my($lchr,$lgb,$lge,$lor,$lxb,$lxe,$lgid)=@{$sexons[$j]};
      my($nchr,$ngb,$nge,$nor,$nxb,$nxe,$ngid)=@{$sexons[$i]};
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

  for my $hsp (@sexons) {
    my $nj= _min($MAXLOOK,scalar(@lhsp));
    my $njdid=0;
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
  
  return($nput);
}

=item problem split case, tiny exons outof order, bad splitgene.

cat kfish2sub5fc18wex.sam | cut -f1,3,4,6,7,8,13,15,16 | grep Funhe2EKm024075t1
Funhe2EKm024075t1.0     KN805694.1      437679  11M     KN805592.1      634080  XS:A:-  XP:i:74 XX:A:423-433
Funhe2EKm024075t1.0     KN805592.1      634080  20M     KN805694.1      437679  XS:A:+  XP:i:86 XX:A:872-891
Funhe2EKm024075t1.1     KN805592.1      630689  17M     KN805694.1      422533  XS:A:+  XP:i:85 XX:A:802-818
Funhe2EKm024075t1.1     KN805694.1      422533  11M     KN805592.1      630689  XS:A:-  XP:i:72 XX:A:915-925
Funhe2EKm024075t1.2     KN805694.1      423171  11M     KN805916.1      281447  XS:A:-  XP:i:72 XX:A:904-914
      ^^ tiny exons,
Funhe2EKm024075t1.2     KN805916.1      281447  217M    KN805694.1      423171  XS:A:+  XP:i:77 XX:A:1086-1302
Funhe2EKm024075t1.3     KN805916.1      281447  217M    KN805694.1      416962  XS:A:+  XP:i:77 XX:A:1086-1302
Funhe2EKm024075t1.3     KN805694.1      416962  206M    KN805916.1      281447  XS:A:-  XP:i:71 XX:A:1091-1296
Funhe2EKm024075t1.4     KN805694.1      416962  206M    KN805592.1      645860  XS:A:-  XP:i:71 XX:A:1091-1296
Funhe2EKm024075t1.4     KN805592.1      645860  226M    KN805694.1      416962  XS:A:+  XP:i:87 XX:A:1103-1328
Funhe2EKm024075t1.5     KN805592.1      645860  226M    KN805694.1      413710  XS:A:+  XP:i:87 XX:A:1103-1328
Funhe2EKm024075t1.5     KN805694.1      413710  13M     KN805592.1      645860  XS:A:-  XP:i:70 XX:A:1329-1341
Funhe2EKm024075t1.6     KN805694.1      413710  13M     KN805916.1      312634  XS:A:-  XP:i:70 XX:A:1329-1341
>> KN805916.1:281447 XX:A:1086-1302 here
Funhe2EKm024075t1.6     KN805916.1      312634  28M     KN805694.1      413710  XS:A:+  XP:i:85 XX:A:1519-1546

.. should be 2 junctions for Split=3 parts
Funhe2EKm024075t1.4     KN805694.1      416962  206M    KN805592.1      645860  XS:A:-  XP:i:71 XX:A:1091-1296
Funhe2EKm024075t1.4     KN805592.1      645860  226M    KN805694.1      416962  XS:A:+  XP:i:87 XX:A:1103-1328

Funhe2EKm024075t1.3b    KN805592.1      645860  226M    KN805916.1      281447  XS:A:+  XP:i:87 XX:A:1103-1328
Funhe2EKm024075t1.3b    KN805916.1      281447  217M    KN805592.1      645860  XS:A:+  XP:i:77 XX:A:1086-1302

.. gff is messy
  Sp1: trg=Funhe2EKm024075t1 380 1341; Sp2: trg=Funhe2EKm024075t1 802 1328; Sp3: Funhe2EKm024075t1 1086 1546;
  insrc=kf3n:gspl3n5h; should cancel just because Split part mRNA targets overlap much
  
KN805694.1      mRNA    413710  438195  -       ID=Funhe2EKm024075t1;locustag=Funhe2EKm024075C;gene=Funhe2EKm024075;Split=1
KN805694.1 tiny exon    413710  413722  -       Parent=Funhe2EKm024075t1;Split=1;trg=Funhe2EKm024075t1 1329 1341;splice=TATC-
KN805694.1      exon    416962  417167  -       Parent=Funhe2EKm024075t1;Split=1;trg=Funhe2EKm024075t1 1091 1296;splice=ACGG-
KN805694.1      exon    422533  422543  -       Parent=Funhe2EKm024075t1;Split=1;trg=Funhe2EKm024075t1 915 925;splice=AGGT-
KN805694.1      exon    423171  423181  -       Parent=Funhe2EKm024075t1;Split=1;trg=Funhe2EKm024075t1 904 914;splice=AGGT-
  >> tiny part 2 exons go here
KN805694.1      exon    434454  434465  -       Parent=Funhe2EKm024075t1;Split=1;trg=Funhe2EKm024075t1 500 511;splice=CAAT-
KN805694.1      exon    437679  437689  -       Parent=Funhe2EKm024075t1;Split=1;trg=Funhe2EKm024075t1 423 433;splice=GCGT-
KN805694.1 tiny exon    438178  438195  -       Parent=Funhe2EKm024075t1;Split=1;trg=Funhe2EKm024075t1 380 397;splice=CAGT-

KN805592.1      mRNA    630689  646070  +       ID=Funhe2EKm024075t1;locustag=Funhe2EKm024075A;gene=Funhe2EKm024075;Split=2
KN805592.1 tiny exon    630689  630704  +       Parent=Funhe2EKm024075t1;Split=2;trg=Funhe2EKm024075t1 802 818;splice=TATC+
KN805592.1 tiny exon    634080  634100  +       Parent=Funhe2EKm024075t1;Split=2;trg=Funhe2EKm024075t1 872 891;splice=AGGC+
    .. ^^ AND target spans out-of-order here, may indicate scaf assembly problem, but should skip?
KN805592.1      exon    645860  646070  +       Parent=Funhe2EKm024075t1;Split=2;trg=Funhe2EKm024075t1 1103 1328;splice=AGGT+
    .. ^^ partial overlap part1 exon
    
KN805916.1      mRNA    281447  312663  +       ID=Funhe2EKm024075t1;locustag=Funhe2EKm024075B;gene=Funhe2EKm024075;Split=3
KN805916.1      exon    281447  281663  +       Parent=Funhe2EKm024075t1;Split=3;trg=Funhe2EKm024075t1 1086 1302;splice=GGTN+ 
    .. dang this^^ exon is also partly in part 2, shouldnt be both places, also in part1
KN805916.1 tiny exon    312634  312663  +       Parent=Funhe2EKm024075t1;Split=3;trg=Funhe2EKm024075t1 1519 1546;splice=AGTT+

=cut
  
  
  
=item blast sumscore

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

=cut

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


