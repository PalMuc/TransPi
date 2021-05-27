#!/usr/bin/perl

# package evigenegff;
package main;

use strict;
use warnings;

=item about evigenegff

  
  ******* IN PROGRESS ******** May discard ........
  
  collected subroutines for EvidentialGene gene gff handling, attributes
  from various perls
  
=cut

#? require Exporter;
# our @ISA = qw (Exporter);
# our @EXPORT = qw (xxxx);
# use Bio::DB::Fasta; # now as get_dna() require  Bio::DB::Fasta, so can use this w/o bioperl  

use constant evigenegff_VERSION  => '20170725'; 
our $DEBUG=0;

# use cdna_proteins;
## from cdna_proteins:
use constant CDNA_TRIMNNN => 1;
our $KEEPSAMECDS=0;
our $RNATYPES='mRNA|ncRNA'; # mrnatypes exontypes genetypes
our $EXONTYPES='exon|CDS'; # mrnatypes exontypes genetypes
our $GENETYPES='UTR|polypeptide|transposon'; # any other gff gene component

sub testgene { # (\@generec, \@otherft)

  # package generic sub, chain to methods desired
}


sub splitPart {
  my($xattr)=@_;  # gff attribute col 9
  my($xid) = $xattr =~ m/(?:ID|Parent)=([^;\s]+)/;  
  ## which split key has precedence? Split=\d or ID_C\d ?
  my($spl)= ($xid =~ /_C(\d+)$/)? $1 : ($xattr =~ m/;Split=(\d+)/)? $1 : 0;
  return($spl,$xid);
}

# fixme for "." strand compare
sub _eqstrand { return($_[0] eq $_[1] || $_[0] eq '.' || $_[1] eq '.' )?1:0; }

sub _sortgene  
{
  # ($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid) == @gff record
  my($ta,$tb)= map{ (m/gene/)?1:(m/mRNA/)?2:(m/CDS|exon/)?3:4; } ($a->[2],$b->[2]);
  return ($ta <=> $tb)
      || ($a->[0] cmp $b->[0]) # ref : should all be same for gene record
      || ($a->[3] <=> $b->[3]) # begin
      || ($a->[4] <=> $b->[4]) # end
      || ($b->[2] cmp $a->[2]) # type: exon > CDS
      ;
}

sub _sortSplitgene {
  # ($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid) == @gff record
  #FIXME? for split genes, Split=1,2,3/3 order by split-part? all mRNA at top or not?
  my($spla,$aid)= splitPart($a->[8]);
  my($splb,$bid)= splitPart($b->[8]);
  
  my($ta,$tb)= map{ (m/gene/)?1:(m/mRNA/)?2:(m/CDS|exon/)?3:4; } ($a->[2],$b->[2]);
  return ($spla <=> $splb) || ($ta <=> $tb)
      || ($a->[0] cmp $b->[0]) # ref : should all be same for gene record
      || ($a->[3] <=> $b->[3]) # begin
      || ($a->[4] <=> $b->[4]) # end
      || ($b->[2] cmp $a->[2]) # type: exon > CDS
      ;
}

sub __revSplitGene {  # rev-sort -strand exons per part
  # ($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr,$gid) == @gff record
  my($spla,$aid)= splitPart($a->[8]);
  my($splb,$bid)= splitPart($b->[8]);
  
  my $ora= ($a->[6] eq '-' or $a->[6] eq -1)?-1:1;
  my $orb= ($b->[6] eq '-' or $b->[6] eq -1)?-1:1;
  return ($spla <=> $splb) # split part
    || ($a->[0] cmp $b->[0]) # ref
    || ($ora * $a->[3] <=> $orb * $b->[3]); # strand*begin
}


sub cds_span {
  my($geneid,$mrnaattr,$cdspanh)= @_;
  my $cdsor= 0; # want Target orient, mRNA is +1 always, but antisense, or cDNA targ can be -1
  my $cdsoff= $cdspanh->{$geneid} if( ref $cdspanh );
  # unless($cdsoff) { ($cdsoff)= $mrnaattr =~ m/offs=(\d+.\d+)/; }
  unless($cdsoff) { ($cdsoff)= ($mrnaattr =~ m/(?:cdsoff|offs)=([^;\s]+)/) ? $1 : ""; }
  my($cdsb,$cdse)= ($cdsoff and $cdsoff=~m/(\d+).(\d+)/) ? ($1,$2) : (0,0);
  if($cdse>0 and $cdse < $cdsb) { ($cdsb,$cdse)= ($cdse,$cdsb); $cdsor= -1; }
  elsif($cdsoff =~ m/(\d+).(\d+):(.)/) { my $cor=$3; $cdsor=($cor eq '-')?-1:+1; } # :. is common, same as :+
  else { $cdsor= +1; }
  return ($cdsb,$cdse,$cdsor);
}

sub phaseCDS {
  my($tb,$te,$cb,$ce,$cdsor, $cdslen)=@_;
  # tb,te= target_begin, target_end; cb,ce,cdsor = cds_begin,_end,_orient(-1,1), running cdslen
  
  # phase calc; FIXME partial5 may need start phase
  my($xend5,$xend3)= ($cdsor < 0)? ($te,$tb) : ($tb,$te); # if input is mRNA, tb == xend5
  my $d5= ($tb >= $cb) ? 0 : $cb - $tb; # pos
  my $c5= ($xend5 > $xend3) ? $xend5 - $d5 : $xend5 + $d5;    				  
  my $d3= ($te <= $ce) ? 0 : $ce - $te; # neg
  my $c3= ($xend5 > $xend3) ? $xend3 - $d3 : $xend3 + $d3;

  my $xwide  = 1 + abs($c3 - $c5);
  $cdslen   += $xwide;
  my $inc3   = $cdslen % 3;
  my $inc5   = ($xwide - $inc3) % 3; # only care about this one
  if ($inc5 == -1) { $inc5 = 2; }
  return( $inc5, $cdslen); # is this right?
}

sub phaseorf {
  my($xonseq,$cdnafull)=@_;
  return(0,0,"",$cdnafull,undef) unless($xonseq);
  my @lorfs;
  my @phase= ($cdnafull and length($xonseq)>3) ? (0,1,2) : (0);
  for my $j (@phase) { 
    my $cdp= ($j==0) ? $xonseq : substr($xonseq,$j);
    my $workseq= ($cdnafull) ? $cdnafull.$cdp : $cdp;
    # cdna_proteins.pm methods
    my @stops  = identify_putative_stops($workseq);
    my @starts = identify_putative_starts($workseq,\@stops);
    my @orfs = get_orfs (\@starts, \@stops, $workseq, '+');
    ## use goodlen not length, no gaps here
    my($longorf) = sort {$b->{goodlen} <=> $a->{goodlen} or $b->{complete} <=> $a->{complete}} @orfs;
    if($longorf) { push @lorfs, [$longorf->{goodlen},$j,$cdp,$workseq,$longorf]; }
    }
  my ($lorf)= sort{ $$b[0] <=> $$a[0] or $$a[1] <=> $$b[1] } @lorfs;
  return ($lorf) ? @$lorf : (0,0,"",$cdnafull,undef);
} 


sub orientTargetCDS {
  my($cb,$ce,$cdsor,$xons)= @_; 
	my($targor, $xrevorder)= (0) x 9; 
	return($cdsor,0,0,0) unless(ref($xons) and $xons->[1] and $xons->[1]->[3]);
  if($ce > 0 and $ce < $cb) { 
    ($cb,$ce)= ($ce,$cb); $cdsor= -$cdsor; #?
  }
  $xrevorder= ($xons->[0]->[3] > $xons->[1]->[3]) ? -1 : 0; # exon sort order, fwd or rev?
  my($tb1)= $xons->[0]->[8] =~ m/(?:Target|trg)=\S+\s(\d+)/; # Target MISSING for CDS exons
  my($tb2)= $xons->[1]->[8] =~ m/(?:Target|trg)=\S+\s(\d+)/;
	if($xrevorder) { $targor=($tb2 and $tb1 > $tb2)? 1 : ($tb1 and $tb1 < $tb2) ? -1 : 0; }
  else { $targor=($tb2 and $tb1 > $tb2)? -1 : ($tb1 and $tb1 < $tb2) ? 1 : 0; }
	
	my $chrstrand= 
	  ($targor > 0 and $cdsor >= 0) ? '+' :
	  ($targor > 0 and $cdsor  < 0) ? '-' :  
	  ($targor < 0 and $cdsor  < 0) ? '+' :  # CDSb == chr-start
	  ($targor < 0 and $cdsor >= 0) ? '-' :  # CDSb == chr-end
	  ''; # ($chror < 0)?'-':'+' ; # cant tell, use input
	  	 
  return($cdsor,$targor,$xrevorder,$chrstrand); #??  
}

sub makeCDSexons {
  my($mrna,$cb,$ce,$cdsor,$gstrand,$xons,$AS_STRING)= @_; 
  
	return() unless($cb and $ce and ref($xons)); #$ce > $cb
  ## UPD17: this should check, return new chr-strand matching CDS-on-genome strand
  ## .. cdsor is needed, but input chror should be ignored. Calc proper chr-or from
  ## cdsor and chr-trans align orient (exon Target order vs exon chr order). See above notes
  ## FIXME5: 16.10.11 : sense=-1; flipor CDS can be bad, wrong calc from exons, offset

	my @cdsgff=(); 
	my ($cdslen,$phase, $cdsfirst, $targor, $xrevorder, $chrstrand)=(0) x 9; 
  $cdsfirst= 1;  
  if($ce > 0 and $ce < $cb) { 
    ($cb,$ce)= ($ce,$cb); $cdsor= -$cdsor; #?
  }
  
  ($cdsor,$targor,$xrevorder,$chrstrand)= orientTargetCDS($cb,$ce,$cdsor,$xons);
  $chrstrand ||= $gstrand;
	  	 
  foreach my $x (@$xons) {
    my @c= @$x; # only this way here
    # if(ref($x) =~ /ARRAY/) { @c= @$x; } # ?? ASARRAY, return same as xons form?
    # elsif($x =~ /\texon\t/) { @c= split"\t",$x; chomp($c[-1]); }	    
    
    my($chrb,$chre,$chrorx)= @c[3,4,6];
    my($tid,$tb,$te)= $c[8] =~ m/(?:Target|trg)=(\S+)\s(\d+)\s(\d+)/;
    next unless($te); # warn/error
    
    if($te > $cb and $tb < $ce) {
      ## fixme again; -strand needs c3,c4 swap also
      # FIXME: still bad, targor < 0, cdsor > 0, +gstrand , from gmap antisense; not xrevorder
      # CDSe hits this:  if($te >= $ce) { .. NOT ($xrevorder and $cdsfirst) > $c[3] += $d; should be c[4] += d

      if($tb <= $cb) { my $d=$cb-$tb;  ## chror < 0 implies chre == tb, col[4]
        if($xrevorder) { if($cdsfirst) { $c[4] -= $d; } else { $c[3] += $d; } }
        else { if($cdsfirst) { $c[3] += $d; } else { $c[4] -= $d; } }
        $cdsfirst= 0;
        }
      if($te >= $ce) { my $d=$te-$ce; ## chror < 0 implies chrb == te, col[3]
        if($xrevorder) { if($cdsfirst) { $c[4] -= $d; } else { $c[3] += $d; } }
        else { if($cdsfirst) { $c[3] += $d; } else { $c[4] -= $d; } }
        $cdsfirst= 0;
        }
       
      ($phase,$cdslen)= phaseCDS($tb,$te,$cb,$ce,$cdsor,$cdslen);

      $c[2]= "CDS";
      # $c[5]=1; # score.. leave as exon score (alignment?)
      # $c[6]= $chrstrand; # do this here? should change exons also
      $c[7]= $phase;  
      $c[8] =~ s/;Target=.*//;
      
      my $xout= ($AS_STRING) ? join("\t",@c)."\n" : \@c;  # which? array or string?
      push @cdsgff, $xout;
    }
  }
 
	return ($chrstrand, \@cdsgff);
}


sub mrna_fixspan {
  my($mrna,$exons)= @_;
  my $changed=0;
  
  my($mpart,$geneid)= splitPart($mrna->[8]);
  my($mref,$oldstart,$oldstop,$oldor)= ($mrna->[0],$mrna->[3],$mrna->[4],$mrna->[6]);
  my($newstart,$newstop,$nrev,$nfwd,$nxok)= (0) x 9;
  
  foreach my $ex (@$exons) {
    my($xr,$xb,$xe,$xo,$xat)= @{$ex}[0,3,4,6,8];
    next if($xr ne $mref); #?? bug
    my($xpart)= splitPart($xat);
    next if($xpart ne $mpart); # Split fix
    if($xo eq '-') { $nrev++; } elsif($xo eq '+') { $nfwd++; }
    ## bug split part mrna get all parts span on same chr .. match Split=n from mrna/exon?
    $newstart= $xb if($newstart == 0 or $xb < $newstart);
    $newstop = $xe if($newstop == 0 or $xe > $newstop);
    $nxok++;
    }
  
  my $newor= ($nrev > $nfwd)?'-':'+';
  unless($newstop==0 or ($newstart == $oldstart and $newstop == $oldstop and $newor eq $oldor)) { 
    # ** Record changes in mrna attr
    $changed++; 
    $mrna->[3] = $newstart; $mrna->[4] = $newstop; 
    $mrna->[6] = $newor; # unless($newor eq $oldor);
    $mrna->[8] =~ s/$/;fixoldspan=$oldstart-$oldstop:$oldor/;
    $mrna->[8] =~ s/;nexon=\d+/;nexon=$nxok/; # add if missing?
    my $hasmix= ($mrna->[8]=~/strandmix=/)?1:0;
    if($nrev>0 and $nfwd>0) {  
      $mrna->[8] =~ s,$,;strandmix=+$nfwd/-$nrev, unless($hasmix);
    } elsif($hasmix) {
      $mrna->[8] =~ s,;strandmix=[^;\s]+,,;
    }
    }
  return $changed;
}

sub generecfix {
  my($generec,$exongff,$exonfix,)= @_;
  my $mrnachanged=0;
  
  my @oldexon= map { my @xnew= @$_; \@xnew; } @$exongff; # clone all
  foreach my $ex (@oldexon) { $ex->[2] = "oldexon"; $ex->[0] =~ s/^/#/; }

  my @generecfix= grep{ $_->[2] ne "exon" } @$generec;
  my @mrna = grep{ $_->[2] eq "mRNA" } @generecfix; # fixme issplit/Splitgene has more mrna parts
  my $issplit= (@mrna>1)? @mrna : 0;
  
  push @generecfix, @$exonfix;
  
  if($issplit) {
    @generecfix= sort _sortSplitgene  @generecfix;  
  } else {
    @generecfix= sort _sortgene  @generecfix; # sorts by genome loc, not stranded
  }
        
  # FIXME Split:  BUG, sets all part-mrna to same enclosing span 
  foreach my $mrna (@mrna) {
    $mrnachanged++ if( mrna_fixspan($mrna,$exonfix) );
  }
  
  $generec= \@generecfix;
  return( $generec, $exonfix, $mrnachanged);
}
  
our $NoStopCodon= 0; # from cdna_proteins

sub putseq {
  my($outh,$id,$seq,$attr)= @_;
  return -1 unless($id and $seq and $outh);
  if($NoStopCodon) { $seq =~ s/\*$//; }
  my $slen=length($seq); # none have \n or \s ?
  $seq =~ s/(.{60})/$1\n/g; $seq.="\n" unless($seq=~m/\n$/);
  $attr||=""; #  aalen=$aalen,$pcds%,$compl; clen=$clen; strand=$crev; offs=$prostart5-$proend3; $cdnah
  $attr.=" len=$slen;" unless($attr =~ m/len=$slen/); # qlen? clen?
  print $outh ">$id $attr\n", $seq;
}


sub putgene {
  my ($houtgff, $generec,$otherft,$flags)= @_; 
  ##  my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid]; 
  my $cc= ($flags and $flags =~ /skip=/)? "#x." : "";
  $otherft ||= [];

  my $checkGeneSpan= ($flags and $flags =~ /mrnachanged/)? 1:0;
  if($checkGeneSpan) { $flags =~ s/mrnachanged=.//;
    my($gene) = grep{ $_->[2] eq "gene" } @$generec;
    my($mrna) = grep{ $_->[2] eq "mRNA" } @$generec; # FIXME Split gene @mrna
    my $didfix= ($gene and $mrna) ? mrna_fixspan($gene,[$mrna]) : 0; 
  }
  
  foreach my $ft (@$generec, @$otherft) # resort together ?
  { 
    if(ref $ft) { 
      my @v= @$ft; 
      $v[8]=~s/$/;$flags/ if($flags and $v[2] eq "mRNA");
      print $houtgff  $cc.join("\t",@v[0..8])."\n" if(@v>4); 
      }
  }
  print $houtgff "\n"; #? spacer row
}

# sub putgene_OLD {
#   my($generec, $geneother, $nchanges) = @_;
#   $geneother ||= [];  my $nput=0;
#   # note generec contains gene, mrna, exons, cds, in orig order; updated locs ..
#   # my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid]; 
#   foreach my $rloc (@$generec, @$geneother) {
#     next unless(ref $rloc); # gene row may be undef
#     print $OUTH join("\t", @$rloc[0..8])."\n"; $nput++;
#   }
#   return $nput; # or nput?
# }


#........ FIXGENE subs / here or not .................................


=item DO_CDSFIX work

  # $DO_CDSFIX works something like SPLITGENEPROTFIX, but per cds-exon test orf with added exondna
  #  .. trim or drop cds-exons that damage orig protein, ie inside  oldoffs span, ignore utr exons?
  #  .. where "damage" is fuzzy measure, 
  #  .. replace all of below SPLIT.FIX, return ($cdna,$cdnalen, $cdnaexonseq, $xfixfix, $xtrimmed)
  # .. should use this only on trasm genomap identified as bad cds2aa but w/ full-ish cdna map

 .. see genefindcds2.pl

  # * first try past bugs worked for 4 cases of long pubaa > short cds .. got middle way to pubaa
 test1.bestof6top5
 Funhe2EKm010574t1	-2433.da	1924,41%,complete	4423,84%,complete,gaps:66  
 Funhe2EKm013880t1	-3011.da	1908,51%,complete	4983,95%,complete,gaps:64  
 Funhe2EKm027743t1	-3194.da	1469,34%,complete	4663,97%,complete	89       
 Funhe2EKm038182t1	-2268.da	1516,40%,complete	3784,97%,partial3	95

    ## this works only for some cases, eg fails if map cover is low.
    ## check and cancel if neworfaalen << oldaalen ? or need prior/nonfix aalen to see if improved enough to keep
    ## and odd cases improve beyond orig pubaa, presumably removing trasm errors or possibly got utrorf * BAD result
# Funhe2EKm036243t1	1980,28%,complete	3329,45%,complete-utrpoor	92  << nofix: cds2aa << pubaa
# Funhe2EKm036243t1 aalen=5292,58%,complete; << cdsfix improved too much??, possibly picked utrorf; this is DSCAM **
#  ID=Funhe2EKm036243t1;cov=92%,20249/22083;pid=98.6;nexon=75;splice=268;trg=Funhe2EKm036243t1 1 22083;
#  gaps=1529,MGap:12449-12561,MGap:12294-12435,MGap:6230-6395,MGap:6181-6203,MGap:12043-12282,MGap:20542-20783,MGap:20904-21410,MGap:20796-20891,;gescore=85;
#  clen=22083;offs=3714-13703;Name=Down syndrome cell adhesion molecule (54%M);
#  cxlen=15876/27197;aalen=5292,58%,complete;cdsoff=491-16369;utrx=3,41;
#  cdsfix=xtrimb25:6194-6474/35,xtrimb49:12987-14485/191,xtrimb49:12836-14294/40,xtrimb49:13102-14254/306,xtrimb50:12983-13586/186,xtrimb50:12873-13400/76,xtrimb50:12891-13324/94,;xtrim=1;
#  aaold=3329,45%,complete-utrpoor
    
#   * fixme genefindcds2.pl update -fixcds needs gap check/trim
#   gap err for Funhe2EKm017647t1
#    /protein_id="Funhe2GD:Funhe2EKm017647p1"  /translation="XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXLLEKAYAKLNG
# from kfish2n3gmap_test4bfcds5.gff
# Funhe2EKm017647t1 Target=Funhe2EKm017647t1 733 2479;cov=70.3;;offs=42-2402;
#  cxlen=1581/1583;aalen=527,99%,partial;protein=XXXXXXXXXXXXXXXXXX
#  cdsoff=3-1583;inshort=xcut12:7;cdsfix=xdrop17:2376-2479/1646,;xtrim=3;aaold=786,95%,complete
# Funhe2EKm017647t1 kfish2n3gmap_test4bfcds5.cdna
# >Funhe2EKm017647t1 aalen=527,99%,partial; cdsoff=3-1583; cxlen=1581/1583; qlen=2485; aaold=786,95%,complete; oid=Funhe2Exx11m006520t2,Funhe2E6bm006529t3;  len=1583;
# TNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
# NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCTCCTGGAGAAAGCTTACG
#   BUT getcdna() below calls this.  
#     my($didtrim,$texondna,$tft,$tstart,$tstop)= trimcdnaexon($j,$exondna,$ft,$rev,$nx1);
# * This bug because of dropped terminal exon .. dropped here? ix=19 missing, ix=18 has TNNN gap, but trim() skipped as it was 2nd, not 1st (rev)


=cut

use constant DOREVCOMP => 1;

sub fixCDS {
  my($oldoffs, $exongff)= @_;

  my ($oldcdsb,$oldcdse)= $oldoffs=~m/(\d+).(\d+)/; # fail unless have oldoffs?
  my ($cdna, $cdnalen, $cdnaexonseq, $xfixfix, $xtrimmed)= (0) x 9;
 
  my (@cdnap,@cdnaxp);
  my ($cdnafull, $cdnaflen, $cdnafend)=("",0,0); 
  my $nexon1= scalar(@$exongff) - 1;
  my $xfixed="";
  $nexon1= -1 unless($oldcdse and $oldcdse>0);
  
  for my $ix (0..$nexon1) {
    my $xft= $exongff->[$ix];
    my @xclone= @$xft; 
    my $cdnaexoni= \@xclone;  # will be changes to this here; clone?

    my $xfixi="";
    my($xtgd,$xtgb,$xtge)=(0,0,0);
    if($$xft[8] =~ /(?:Target|trg)=(\S+).(\d+).(\d+);/) { ($xtgd,$xtgb,$xtge)=($1,$2,$3); }
    
    my($cdnai,$rev,$cleni)= getexondna($cdnaexoni, DOREVCOMP); 
    
    $cdnafend= length($cdnafull);
    my $cdnafbeg= 1 + $cdnafend;
    $cdnafend += $cleni;
    ($xtgb,$xtge)= ($cdnafbeg,$cdnafend) unless($xtge);
    
    #?? should use xtgb,xtge here .. mapping target spans of exon, not local cdna string span
    my $overcds=   ($xtge > $oldcdsb and $xtgb < $oldcdse) ? 1 : 0;
    my $insidecds= ($xtgb >= $oldcdsb and $xtge <= $oldcdse) ? 1 : 0;
    my $utrbeg= ($cdnafbeg < $oldcdsb)?1:0;
    my $utrend= ($cdnafend > $oldcdse)?1:0;
    
    my $okcds=0; # need oldaalen, oldoffs cds offs / span to validate
    unless($overcds) { 
      $okcds=1;  # utr??, only append cdnai
    }
    for(my $try=0; $okcds == 0 and $try < 3; $try++) {
      $okcds=1; # if(longorf spans orig cds offset ?)

      my($workseq,$longorf,$llen,$lj,$lexon);
      ($llen,$lj,$lexon,$workseq,$longorf)= phaseorf( $cdnai, $cdnafull);
      $cdnai= $lexon; $cleni= length($lexon);
      if($lj) { $cdnafbeg += $lj;
        $xtrimmed++; #<< need to register exon change for mrna update; also? flag in exon attr?
        if($rev) { $$cdnaexoni[4] -= $lj; } else { $$cdnaexoni[3] += $lj; } 
      }
      
      my $BEGCUT= 0.33; # 0.25; # : 0.15,0.25,0.33
      my $ENDCUT= 0.66; # 0.75; # : 0.85,0.75,0.66

      $cdnafend= length($workseq);
      # orf fields: length, goodlen, complete, start,stop, orient, innerstop, protein, sequence == cds
      my $lend= ($longorf) ? $longorf->{stop} : 0; $lend||=0; # stop too soon?
      if($lend < $cdnafbeg) { # this means exon is utr, after orf, treat like not overcds ?
        $lend= $cdnafbeg + int(0.10 * $cleni); # punt?
      }  
      
      my $xokspan= $lend - $cdnafbeg; #? +2? is lend/stop 1st of codon3 ?
      $xokspan=0 if($xokspan<0);  ## ^^ NEGative xokspan ?? cdnafbeg supposed to be exon start in workseq. lend orf stop in workseq 

      if(not $overcds) { $okcds=1; } #? insidecds 
      elsif($lend <= $cdnafbeg) { #bug where?
        $okcds=1; # this is utr exon, not overcds
        
      } elsif($lend < $cdnafend-2) { #  - 3?
        if($cleni < 39) { #? bad?
          $okcds= -2; # drop, exon too short to cut
        
        } elsif($cdnafend > $oldcdse and $xokspan >= $BEGCUT*$cleni) { 
          $okcds= 1;   # last cds exon, dont cut unless start is poor
          
        } elsif($xokspan <= $BEGCUT*$cleni) { #? $lend < 9 + $cdnafbeg : 0.15,0.25,0.33 give diff effects. what?
          $okcds=0; 
          my $trimb= $xokspan + 3; #?? need to do phase checks w/ trims; simple orf check here w/ trim?
          my($llen,$lj,$lexon,$lcdna,$longorf)= phaseorf( substr($cdnai,$trimb), $cdnafull);
          $cdnai= $lexon; $trimb += $lj; # lj=0,1,2
          $xtgb+= $trimb;
          if($rev) { $$cdnaexoni[4] -= $trimb; } else { $$cdnaexoni[3] += $trimb; }
          $xfixi.="xtrimb$ix:$xtgb-$xtge/$trimb,";
          
        } elsif($xokspan >= $ENDCUT*$cleni) {  ## 0.75,0.66 what?
          $okcds=0; 
          my $trime= $cleni - ($xokspan + 3); #?  
          $trime=0 if($trime < 0);
          $cdnai= substr($cdnai,0,$xokspan-3);  
          $xtge -= $trime;
          if($rev) { $$cdnaexoni[3] += $trime;} else { $$cdnaexoni[4] -= $trime; } # rev problems! # change exon.stop
          $xfixi.="xtrime$ix:$xtgb-$xtge/$trime,";
          
        } else { # drop exon if error in middle?
          $okcds= -1; 
        }
        $cleni= length($cdnai);
      } else {
        $okcds=1;
      }
              
      if($okcds>0) { $cdnafull= $workseq; }
      elsif($okcds<0) {
       if($okcds == -1 and ($ix == 0 or $ix == $nexon1)) { $okcds=1; $cdnafull= $workseq; $xfixi.="xdropc$ix:$xtgb-$xtge/$lend,";}#??
       else { $cdnai=""; $cdnaexoni= undef; $xfixi.="xdrop$ix:$xtgb-$xtge/$lend,"; }
       }
    }
    
    #......
    # split bug: xdrop exon on other part, leaving only mRNA of split=2 .. should convert to utr-exon not drop
    
    $xfixed.=$xfixi if($xfixi);
    if($cdnai and ref $cdnaexoni) {  
      $$cdnaexoni[8].= ";cdsfix=$xfixi" if($xfixi); 
      push @cdnap, $cdnai ; # use cdnafull only?
      push @cdnaxp, $cdnaexoni; # exon feat to keep, maybe changed start/stop
    }
    
  }  # each exon
  
  if(@cdnaxp>0) {
    my $xtrim1;
    ($cdna,$cdnalen, $cdnaexonseq, $xfixfix, $xtrim1)= getcdna2( \@cdnaxp, 1, 0, 0, CDNA_TRIMNNN);
    $xtrimmed += $xtrim1; # xtrimmed used above
    $xtrimmed++ if($xfixed);  
  }

  # my($cdna, $cdnalen, $cdnaexonseq, $xfixfix, $xtrimmed)= (0) x 9;
  # $addattr .= ";cdsfix=$xfixed" if($xfixed); #? leave to caller
  return($cdna, $cdnalen, $cdnaexonseq, $xfixfix, $xtrimmed, $xfixed); # addattr xfixed xtrimmed
}

# use constant SPLITGENEPROTFIX => 1;  
sub get_splitgenecds {
  my($cdnalen, $issplit, $exongff, $addattrh)= @_;
  my ($cdna, $cdnaexonseq, $xfixfix, $xtrimmed)= (0) x 9; # returns, same as fixCDS
  my $nexon= @$exongff;
  
  if($cdnalen>0) { } # done DO_CDSFIX now
  elsif(not $issplit) { # not split, needs cdnalen
    ($cdna, $cdnalen, $cdnaexonseq, $xfixfix, $xtrimmed)= getcdna2( $exongff, 1, 0, 0, CDNA_TRIMNNN);

  } elsif($issplit) {
    my ($lp,$lid)= (0,0);
    our(@xp,@cdnap,@cdnaxp,@partexons);
    
    sub pushpart {
      our(@xp,@cdnap,@cdnaxp,@partexons);
      my($cdnai,$cleni,$cdnaexonseq,$xfixfix, $xtrimmed)= getcdna2( \@xp, 1, 0, 0, CDNA_TRIMNNN);
      my $cdnaexoni= [@xp];
      $cdnaexoni= $xfixfix if($xtrimmed); ## NEED mRNA/gene span update also ****
      push @cdnap,$cdnai; push @cdnaxp,@$cdnaexoni; @xp=(); 
      push @partexons, $cdnaexoni;
    }
    
    for( my $i=0; $i<$nexon; $i++) {
      my $x= $exongff->[$i];
      my($spla,$aid)= splitPart($x->[8]);
      pushpart() if($spla ne $lp and @xp);
      push @xp,$x;  $lp=$spla; $lid=$aid; 
      }
    pushpart() if(@xp);
      
    # assume forward mrna, find orfs adding each part w/ phase NNN tests
    my $cdnafull="";  
    for my $ip (0..$#cdnap) {
      my $cdp= $cdnap[$ip];
      my $partexons= $partexons[$ip]; # @$partexons
      my($llen,$lj,$lexon,$lcdna,$longorf)= phaseorf( $cdp, $cdnafull); 
      $cdnafull= $lcdna;   
      if($lj) {  $$addattrh .= ";sphase=$ip:$lj";
        $xtrimmed++; #<< need to register exon change for mrna update; also? flag in exon attr?
        my $rev= ($partexons->[0]->[6] eq '-')?1:0;
        if($rev) { $partexons->[0]->[4] -= $lj; } else { $partexons->[0]->[3] += $lj; } 
      }
    }
    $cdna= $cdnafull; $cdnalen= length($cdnafull);
    $cdnaexonseq= \@cdnap; # is this right? need ?
    $xfixfix= \@cdnaxp; 
  } # issplit, SPLITGENEPROTFIX
  
  return($cdna, $cdnalen, $cdnaexonseq, $xfixfix, $xtrimmed); # addattr xfixed xtrimmed
  # $cdna, $cdnalen, $cdnaexonseq, $xfixfix, $xtrimmed
}  
  
  # 201511. add completeCDS extend-ends option, ix==0 and ix==$nexon1, look for completing cdsbases of partial cds

sub completeCDSb { 
  my($geneid, $cdnafull, $exongff, $revgene, $oldProstart5, $oldProend3, $oldprot)=@_;
    # drop params: $oldoffs,$oldStart_b,$oldStart_e; ADD: oldprot
    ## drop bad ,$EXTcdna return
    # ($extbest,$EXTexongff,$extaddattr,$ext5,$ext3)= 
    #   completeCDSb($geneid, $cdna, $cdnaexons, $revgene, $oldProstart5,$oldProend3);
  
  my $addattr="";
  my $XOFF=270;  # 180 too short ??
  my $UTRSLOP=27; # dont offset where UTR exists beyond partial CDS end ?    

  ## bug? reduced complete count for -rev; check/reset?
  my @xbloc= @{$exongff->[ 0]}[0,3,4,6,8]; # ($xr,$xb,$xe,$xo,$xa)
  my @xeloc= @{$exongff->[-1]}[0,3,4,6,8];
  if($xbloc[3] eq "-") { $addattr.=",rerev" unless($revgene); $revgene=1; }
  else { $addattr.=",unrev" if($revgene); $revgene=0; }
  if( ($revgene and $xbloc[1] < $xeloc[1]) or (not $revgene and $xbloc[1] > $xeloc[1]) ) { 
    $addattr.=",resort";  my @xrev= reverse @$exongff; $exongff=\@xrev; }
  if( ($revgene and $oldProstart5<$oldProend3) or (not $revgene and $oldProstart5>$oldProend3) ) { 
    $addattr.=",revcdsbe"; ($oldProstart5,$oldProend3)= ($oldProend3,$oldProstart5); }
  
  my $nexon1= scalar(@$exongff) - 1;
  our $cdnaext= $cdnafull;  # only extend ends..
  my @extongff= @$exongff; # copy all
  my ($ext5try,$ext3try)=(0,0); 

  #NOTE exons here are rev-sorted, exon[0] = end5 of CDS
  sub exoffb{ my($ex,$ic,$xof)=@_; $ex->[$ic]+= $xof; $ex->[8]=~s/$/;cexoff=$xof/; }
  sub exoffadd{ 
    my($ex,$ic,$jc,$xof)=@_; 
    our $cdnaext;
    my($ob,$oe,$oo)= @{$ex}[3,4,6];
    my $xend = ($jc < $ic)?1:0; 
    my $xapend=($oo eq "-")? !$xend : $xend;
    $ex->[$ic]+= $xof; 
    $ex->[$jc]= ($xend)?$oe+1:$ob-1; # shift to end point
    my($addseq)= getexondna($ex, DOREVCOMP); 
    my $ret= length($addseq); $ret= -$ret if($xof<0);
    if($xapend) { $cdnaext= $cdnaext . $addseq;  } # bad for revgene.. needs opposite
    else { $cdnaext= $addseq . $cdnaext;  } 
    $ex->[$jc]= ($xend)?$ob:$oe; # revert non-end
    $ex->[8]=~s/$/;cexoff=$xof/; 
    return $ret;
    }

  # only these: if($ix==0 or $ix==$nexon1)  
  for my $ix (0,$nexon1) {
    my $xft= $exongff->[$ix];
    my($xb,$xe)= @{$xft}[3,4];  
    my $cdnaexoni= [@$xft]; # @xclone= @$xft; $cxi= \@xclone;  # will be changes to this here; clone?
    $extongff[$ix]= $cdnaexoni;
    if($revgene) {
      if($ix==0 and ($xe <= $oldProstart5 + $UTRSLOP)) {
        $ext3try= exoffadd($cdnaexoni,4,3,$XOFF);
      }
      if($ix==$nexon1 and ($xb >= $oldProend3 - $UTRSLOP)) { 
        $ext5try= exoffadd($cdnaexoni,3,4,-$XOFF);
      }
   
    } else {
      if($ix==0 and ($xb >= $oldProstart5 - $UTRSLOP)) {
        $ext5try= exoffadd($cdnaexoni,3,4,-$XOFF);
      }
      if($ix==$nexon1 and ($xe <= $oldProend3 + $UTRSLOP)) { 
        $ext3try= exoffadd($cdnaexoni,4,3,$XOFF);        
      }
    }
  }

  if($cdnaext ne $cdnafull) {
    # oldStart_b needs cdna_proteins:KEEPSAMECDS; prefer keep same CDS exons but can extend/shorten protein bounds
    $KEEPSAMECDS=1;
    my($orfprot, $prostart5, $proend3, $bestorf)= 
      getBestProt("partial", $cdnafull, $exongff); ##, $oldStart_b,$oldStart_e);
    my($orfproti, $prostart5i, $proend3i, $bestorfi)= 
      getBestProt("partial", $cdnaext, \@extongff );  

    if($oldprot) {
      ## FIX for bestaa=pubaa, oldprot .. must test agains that; BUT problem pubaa: XXXX, other?
      my $sameold= index($orfproti,substr($oldprot,2));
      if($sameold<0 and $oldprot=~/XX/) { # check gaps?
        my @oldaa= split /X+/, $oldprot; my @oai;
        for my $oaa (@oldaa) { my $i=index($orfproti,$oaa); push @oai, $i if($i>=0); }
        if(@oai == @oldaa) {
	        my($ob,$oe)= @oai[0,-1]; my $olen=$oe + length($oldaa[-1]) - $ob;
          $oldprot= substr($orfproti,$ob,$olen);
	        }
      }
      $orfprot=$oldprot;
    }

    ## test stats
    my($aawi,$aaw)= (length($orfproti),length($orfprot));
    my $aadif= $aawi - $aaw; 
    my($trwi,$trw)= (length($cdnaext) , length($cdnafull)); 
    my($aaci,$aac)= ($bestorfi->{complete}, $bestorf->{complete});
    my($relstart5i,$relend3i)= ($prostart5i + $ext5try, $proend3i + $ext5try); # relative to -XOFF shift
    my $ooff= join "-", $prostart5, $proend3;
    my $noff= join "-", $relstart5i,$relend3i;
    if($relstart5i <= 0) { $noff="EXT$noff"; }
    
    my $extbest=0;
    $extbest=1 if($aawi>$aaw); # not enough, for part[53] need other offset same
    $extbest=0 if($relend3i < $proend3-$UTRSLOP); # diff locations, happens w/ for same prot
    $extbest=0 if($relstart5i > $prostart5+$UTRSLOP);
    my $sameprot= (index($orfproti,substr($orfprot,2))<0)?0:"AASAME";  # substr() cuts possible part5 codon
    $extbest=0 unless($sameprot);
    # if($extbest and $sameprot) { $extbest=0; } # substr() cuts possible part5 codon

    warn "#cds.complete $geneid fix=$extbest aaquals=$aaci/$aac, aadif=$aadif,$aawi/$aaw,$sameprot, "
        . "coff=$noff/$ooff, trw=$trwi/$trw, at=$addattr \n"; ## if $DEBUG; , oldaaq=$oldaalen, oldoffs=$oldoffs; 

    if($extbest) { 
      ## trim off extra XOFF at ends, leave NO UTR or tiny UTR ? FIXME revor swap *
      my($ext5,$ext3)=(0,0); # use $ext5try,$ext3try above?
      my($cdsOfExt, $cdsattrOfExt)= (undef,"");
      ## try again:
      # my ($cdslong, $attrL, $pcodeL, $maxutrL)= getCDSgff( $exongff, orfParts($longorf));
      ($cdsOfExt, $cdsattrOfExt)= getCDSgff( \@extongff, orfParts($bestorfi)); # see above
      my($cdsb,$cdse)= $cdsattrOfExt=~m/cdsspan=(\d+).(\d+)/;
      if($cdse) {
        my($oxb,$oxe);
        my $crev=0; if($cdse < $cdsb) { $crev=-1; ($cdsb,$cdse)=($cdse,$cdsb); }
        if($revgene) {  
          $oxb= $exongff->[-1]->[3];  $oxe= $exongff->[0]->[4]; # rev sort
          if($cdsb<$oxb) { $ext3= -($cdsb - 3 - $oxb); } ## -() for rel..5i compat, -ext offset
          if($cdse>$oxe) { $ext5= -($cdse + 3 - $oxe); }
          }
        else { 
          $oxb= $exongff->[0]->[3];  $oxe= $exongff->[-1]->[4]; 
          if($cdsb<$oxb) { $ext5= $cdsb - 3 - $oxb; }
          if($cdse>$oxe) { $ext3= $cdse + 3 - $oxe; }
          }
      }
      
      if($ext5 == 0 and $ext3 == 0) {
        # NOTE exongff are rev order here; 5' is first
        if($relstart5i <= 0) { $ext5= $relstart5i - 3; } # 
        if($relend3i >= $trw){ $ext3= $relend3i + 3 - $trw; } # is this proper offset?
      }
      
      if($ext5 != 0) { if($revgene){ exoffb($exongff->[0],4,-$ext5); } else { exoffb($exongff->[0],3,$ext5); } }
      if($ext3 != 0) { if($revgene){ exoffb($exongff->[-1],3,-$ext3); } else { exoffb($exongff->[-1],4,$ext3); } }
 
      $addattr = "cdsfix=ext:$aaci/$aac,$noff/$ooff,$aawi/$aaw$addattr;";
      return($extbest,$exongff,$addattr,$ext5,$ext3);
      }
  }
  #?? debug return cdsnofix= attr??
  
  return(0);
}


sub cdnain_prot {
  my( $geneid, $geneidfix, $cdnaseq, $cdnalen, $mrnaat, $oldaalen, $orfprot, $exongff)= @_; # cdnagfflen  geneid 
  
  ## FIXME: cdnain is strand-less ; rev for fixstrand ?
  ## FIXME2: test cdnain before / after intronfix, intronfix can be bad.
  
  my $cdnain  = ""; 
  my $cdnainlen= 0; my $cdnainNotPrimaryId= 0; my $cdnainIsBest=0;
  my $cdnagfflen= $cdnalen;
  $cdnain= $cdnaseq; #??
  # if($cdnaseq) {  # add 2011
  #   # $cdnain= $cdnaseq{$geneid} || $cdnaseq{$geneidfix}; 
  #   $cdnain= $cdnaseq->{$geneid};
  #   if( not $cdnain and ($geneid ne $geneidfix) ) { $cdnain= $cdnaseq->{$geneidfix}; $cdnainNotPrimaryId=1 if($cdnain); }
  #  } 
   
  return(0) unless($cdnain);
  
  $cdnainlen= length($cdnain);

  use constant CIN_SAMESTRANDasMAP => 1;
  use constant CIN_BOTHSTRANDS => 0;

  my $cdnain_isrev=0; # doesnt mean samestrand-as-map, but just that cdnain was reversed.
  my $cstrandIsOpposite= 0; # 1 if antisense opposite, 0 == sense same as genome map; add -1 == dont know/cant align?
     ## FIXME6: bug in cstrandIsOpposite, reporting same Sense strand when opposite.
     ## .. probable failure of $cdnaindex not finding cdnagff inside cdnain.
     
  my($orfproti, $prostart5i, $proend3i,$bestorfi, $utrorfi);
  my $cdnaindex= -1;
  
  ORFCDNAIN: { 
  my $cdnainrev= revcomp($cdnain);
  my $validgffstrand= ( @$exongff > 1 ) ? 1 : 0;

  ## FIX: 160530.  # gmap wrong for mRNA mapping; fixme 1605
  my $orfaalen= length($orfprot);
  my ($oldaalend)= ($oldaalen=~m/(\d+)/)?$1:0;
  my $badgfforf= ($mrnaat =~ m/;sense=-1/ or ($oldaalend * 0.80 > $orfaalen))?1:0;
  if($badgfforf) { $validgffstrand=0; }

  my $cstrand_gstrand= 0;
  my $cdnagffgoodlen= $cdnagfflen;
  my $cdnagood= $cdna; # only for index cstrand_gstrand test
  if(index($cdnagood,'N') >= 0) { ## No ($cdnagffgoodlen < $cdnagfflen)
    $cdnagood =~ s/^N+//; $cdnagood =~ s/N+$//; 
    $cdnagffgoodlen=length($cdnagood);
    my $xi= index($cdnagood,"NN"); 
    if($xi>0) { 
      my $xe= rindex($cdnagood,"NN");
      if($xi>$cdnagffgoodlen/2) { $cdnagood= substr($cdnagood,$xe+2); }
      else { $cdnagood= substr($cdnagood,0,$xi); }
    }
  } 
  
  if( ($cdnaindex= index($cdnain,$cdnagood)) >=0 ) {  $cstrand_gstrand= 1; }
  elsif( ($cdnaindex= index($cdnainrev,$cdnagood)) >=0 ) { $cstrand_gstrand= -1; }
  if($cdnaindex < 0) { 
    ## FIXME6: try harder to match strands; indel problems
    my ($xfwd,$xrev)=(0,0);
    foreach my $xs (@$cdnaexonseq) {
      if(index($cdnain,$xs) >=0 ) { $xfwd += length($xs); }
      elsif(index($cdnainrev,$xs) >=0 ) { $xrev  += length($xs); }
    }
    $cstrand_gstrand= ($xfwd>$xrev)? 1 : ($xrev>$xfwd)? -1 : 0;
  }
  
  ## here use goodlength tests for cdnain,cdnagff
  ## testboth is problem now: getting many Anti of gene-mashup (fwd+rev), want only cases of gapfill tr?
  my $testboth= 0; #  2=best of both 
  if($fixstrand) { } # DONT test both in this strand loop
  elsif( $validgffstrand and $cdnainlen < 1.75 * $cdnagffgoodlen) { $testboth=0; }
  elsif( $cdnainNotPrimaryId or ($cdnainlen >= 1.75 * $cdnagffgoodlen)) { $testboth= 2; } # or not $validgffstrand
  elsif( $cstrand_gstrand == 0 ) { $testboth= 2; } #? or leave 0
  
  # gmap wrong for mRNA mapping; fixme 1605
  if($badgfforf) { $testboth=2; }

  if(CIN_SAMESTRANDasMAP) { }   
  elsif(CIN_BOTHSTRANDS) { $testboth=2; }

  my ($longorf,$longfull,$orfs);
  if( $testboth == 2 ) {    
    ($longorf,$longfull,$orfs) = getAllOrfs($cdnain,"fwd"); # ,"dropnnn"
    my ($rlongorf,$rlongfull,$rorfs) = getAllOrfs($cdnainrev,"fwd"); # ,"dropnnn"
    ($cdnain_isrev)= bestorf_test($longorf,$rlongorf);
    if($cdnain_isrev) {
      # $cstrand_gstrand= -1 if($cstrand_gstrand == 0); 
      ## report if cstrand_gstrand unknown: cstrandIsOpposite=-1 ??
      $cstrandIsOpposite=1 if($cstrand_gstrand == 1);
      $cdnain= $cdnainrev;
      ($longorf,$longfull,$orfs)= ($rlongorf,$rlongfull,$rorfs);
    } else {
      $cstrandIsOpposite=1 if($cstrand_gstrand == -1);
    }
  } elsif( $cstrand_gstrand == -1 ) {
    $cdnain_isrev=1;
    $cdnain= $cdnainrev;
    ($longorf,$longfull,$orfs) = getAllOrfs($cdnain,"fwd"); # ,"dropnnn"
  } else {  ##  $cstrand_gstrand == 1
    $cdnain_isrev=0;
    ($longorf,$longfull,$orfs) = getAllOrfs($cdnain,"fwd"); # ,"dropnnn" 
  }
  
  ($orfproti, $prostart5i, $proend3i, $bestorfi, $utrorfi)= 
    getBestProtOfOrfs($longorf,$longfull,$orfs, 
      "partial", $cdnain, $exongff, $oldStart_b,$oldStart_e);
  }
  

  if($orfproti ne $orfprot) { 
    ($cdnainIsBest,undef,undef) = bestorf_test($bestorf,$bestorfi);
     # FIXME bestorf_test: option too high -ratiocdna 1.25; at least should replace cdna XXXX gaps with cdnain perfect
     # test for near-sameprot but for XXX gaps? 
     
    my $cdnainProblems = $cdnainNotPrimaryId or ($mrnaat =~ /chimera/);
    
    if($cdnainIsBest and $cdnainProblems) { 
      if($cdnainNotPrimaryId) { ## and $mrna->[8] =~ /chim\d=(\w+[^;\s]+)/
        my $ingood= $bestorfi->{goodlen} || 1;
        my $bgood = $bestorf->{goodlen};
        $cdnainIsBest= 0 if($bgood/$ingood < 0.25); # bad test.. really need other part's stats here          
      }
    }
  }
    
  if($DEBUG) {  #  and not $cdnainIsBest; debug set note; add orig cxlen;aalen here if cdnainIsBest?
    my $trd = ( $cstrandIsOpposite ) ? "trAnti": "trSense";
    $trd .= ($cdnain_isrev) ? "Rev" : "Fwd"; # less useful than anti/sense
    my( $uprot, $ustart, $uend) = orfParts($bestorfi);
    my $ulen=  $bestorfi->{length};  
    my $ualen= length($uprot); $ualen-- if(substr($uprot,-1) eq '*');
    my $upcds  = ($cdnalen>0 && $ulen>0) ? int(100*$ulen/$cdnalen) : 0;
    my $ucompl= $bestorfi->{complete};
    $ucompl= ($ucompl==3)?"complete":($ucompl==2)?"partial5":($ucompl==1)?"partial3":"partial";
    $addattr .=";" if($addattr); 
    $addattr .= "cdnaorf=$ualen,$upcds%,$ucompl,trbest:$cdnainIsBest,$trd"; #? ;aacdnaoff=$ustart-$uend;aacdnaprot=$uprot";  
  }
  
  # UPD1707: also test mrna-cds seq align to chr-cds, not cdnainIsBest (as option?)
  # .. need to know if chr-cds is covering at least parts of same cds as mrna .. if not, diff gene
  
  if($cdnainIsBest) {
    # replace?, notice, may need to change CDS like intronfix : depends on diff
        
    my $aadif= length($orfproti) - length($orfprot); 
    my $gadif= $bestorfi->{goodlen} - $bestorf->{goodlen}; 
    # my $trdif= length($cdnain) - length($cdna); 
    my $trdif= $cdnainlen - $cdnagfflen;  
    $cdnalen= $cdnainlen; # change it
    
    my $trdif1=$trdif;
    map{ $_="+".$_ if($_>0); $_= (($_==0)?"eq:":"NE:").$_; } ($aadif,$gadif,$trdif);
    my $aac= $bestorfi->{complete} - $bestorf->{complete}; # complete == 3,2,1,0
    
    # also compare cdna vs cdnain : size, where diff (inside, ends?)
    my $nx0= @$exongff; # no way to count exons in cdnain here.
    my $trin=0; my $nxeq=0;
    
    ## change cdnain_isrev to cstrandIsOpposite here?
    my $trd= ($cstrandIsOpposite) ? "trdAnti" : "trdSense"; # "trSens" ?
    $trd .= ($cdnain_isrev) ? "Rev" : "Fwd"; # less useful than anti/sense
    
    if($cdna eq $cdnain) { $trd .= "eq"; $nxeq=$nx0; }
    elsif( $cdnaindex >=0 ) { # from index($cdnain|cdnainrev, $cdna)
      $trd .= "$trdif,in$trin";
      $nxeq= ($nx0==1)?1:$nx0-1;
      if($cstrandIsOpposite) {}
    }
    else { 
      $trd .= $trdif; # * check each \@exongff > exseq index cdnain ? report if ends or inner x diff
      my ($ix,$le,$xfwd,$xrev)=(0,1,0,0); 
      foreach my $xs (@$cdnaexonseq) {
        $ix++; my $xi= index($cdnain, $xs);  my $xr=""; 
        if($xi>=0) { $xfwd+= length($xs); }
        else { $xi= index($cdnain, revcomp($xs)); if($xi>=0) { $xr="c"; $xrev+=length($xs); }  }
        if($xi>=0) { my $xe= length($xs)+$xi; my $xb=$xi+1; $trd.=",xeq$ix:$xr$xb-$xe"; $nxeq++; 
           if( $le>1 and (my $g= $xb - $le) > 0 ) { $trd.= "g$g"; } $le=1+$xe; }
        else { $trd.=",xne$ix"; }
        my $xgap= $xs =~ tr/N/N/;  $trd .="N$xgap" if($xgap>0);
      }
      if($xrev > $xfwd and not $cstrandIsOpposite) { $trd =~ s/trdSense/trdAntix/; } # ERROR??
    }

    #check further if $nxeq > 0; xinner change?  count NNN gaps in cdna, cdnain; look for gaps at end if inner match
    if($nxeq > 0) {
      my $gaps = $cdna =~ tr/N/N/;
      my $gapi = $cdnain =~ tr/N/N/;
      if($trdif1>25) { 
        my($cdna2, $cdnalen2)= getcdna2( $exongff, 0, 100, _max(100,$trdif1));
        $gaps = $cdna2 =~ tr/N/N/;
      }
      if($gaps > $gapi) { $trd.= ",tN:$gaps"; }
    }       
    
    my $cov= ($mrnaat =~ m/(cov|Coverage)=(\d+)/) ? $2 : 0;
    $trd .= ",tcov:$cov" if($cov); # cov < 99 explains
    my $changenote= "cdnabest=aa$aadif,gd$gadif,dfull:$aac,$trd,nx0:$nx0";
    $addattr .=";" if($addattr); $addattr .= $changenote; $changed++;
    
    # FIXME HERE >> want CDS exons computed before cdnain . this changes to oldCDS prostart/end
    ($orfprot, $prostart5, $proend3, $utrorf)= ($orfproti, $prostart5i, $proend3i, $utrorfi);
    
  } #cdnabest
    
  return(xxxxxxxxxxx); # changed changenote $orfprot, $prostart5, $proend3, $utrorf
}

my $USE_CHRMAP_CDS= 1;    
sub add_cdsprot {
  my($mrna, $exongff, $orfprot, $prostart5, $proend3, $cdnalen, $oldCDS)= @_;
 
  my ($prostart5COMP, $proend3COMP)= ($prostart5, $proend3); # vals BEFORE cdnainIsBest
 
  if($orfprot) {

    my($annew,$anold)= stripOldAnnot($mrna) if($REANNOTATE); # also changes mrna
    
    #o my ($cdsgff, $cdsattr)= getCDSgff( \@exongff, $orfprot, $prostart5, $proend3, $cdnalen);    
    #^^ Change to makeCDSexons() + extras
    #** Change to cdna_protein.pm getCDSgff2() for now ..
    my ($cdsgff, $cdsattr);    
    if($USE_CHRMAP_CDS) {
      # UPD: computed CDS not cdnainIsBest CDS
      ($cdsgff, $cdsattr)= getCDSgff2( \@exongff, newOrf(protein=>$orfprot, start=>$prostart5COMP, stop=>$proend3COMP, cdnalen=>$cdnalen) );
    } else { # 
      ($cdsgff, $cdsattr)= getCDSgff2( \@exongff, newOrf(protein=>$orfprot, start=>$prostart5, stop=>$proend3, cdnalen=>$cdnalen) );
      # was getCDSgff(\@exongff, $orfprot, $prostart5, $proend3, $cdnalen);
    }

    
    # test $prostart5, $proend3 vs old
    
    my $diff=0; 
    $diff=1 if($REANNOTATE or $cdnainIsBest);
    
    my $diffcancel=0; my $oldstat="";  
    if(ref $cdsgff and @$cdsgff > 0) {
      my @newCDS; # = sort _sortgene @$cdsgff; # do here.
      if($issplit) {
      @newCDS= sort _sortSplitgene @$cdsgff; # do here.
      } else {
      @newCDS= sort _sortgene @$cdsgff; # do here.      
      }

      my @oldsave=();
      ## reversed CDS look wrong at end exons (always both ends? for nc>2)    
      if(@oldCDS > 0) {
        my($oldal, $newal, $oldprostart5, $oldproend3)= (0,0,0,0);
        $oldprostart5= $oldProstart5; $oldproend3= $oldProend3; # above, dont need both
        # compare, report in $addattr
        if($issplit) {
        @oldCDS= sort _sortSplitgene @oldCDS; # dang2    
        } else {
        @oldCDS= sort _sortgene @oldCDS; # dang
        }
        
        my $newprot= ($cdsattr =~ m/protein=([^;\s]+)/) ? $1 : "";
        if (@oldCDS == @newCDS) { $oldstat="ocds=eqn"; } else {  $oldstat="ocds=NEn"; $diff++; }
        
        if($oldprot) { 
          (my $op=$oldprot) =~ s/\*$//; (my $np=$newprot) =~ s/\*$//; 
          $oldal= length($op);
          $newal= length($np);
          my $oldap= ($cdnalen>0) ? int(0.5 + 300 * $oldal / $cdnalen) : 0;

          my $eq=($op eq $np)?1:0; # count STOP also
          if($eq and ($newprot =~ /\*$/) and not ($oldprot =~ /\*$/)) { $eq=0; }
          
          $diff++ unless($eq); 
          my $da= $newal - $oldal; $da="+$da" if($da>0);
          $oldstat .= ($eq) ? ",eqaa" : ",NEaa:$da"; 
          
          my $astat=0;
          if($oldprot =~ /^M/) { $astat |= 1; }
          if($oldprot =~ /\*$/) { $astat |= 2;} # ** AUG proteins lack '*' unless pasa-updated
          elsif($geneid =~ /AUG/) { $astat |= 2; } #  also check for inner X == augustus-fake for stop ?

          my $istop= index($oldprot,'*');
          ## find single 'X' inside prot, but allow this for NNN genome: SFXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXQP
          if($geneid =~ /AUG/ and (my $ix=index($oldprot,'X')) >0 ) { 
            if(substr($oldprot,$ix,3) eq "XXX") { } # ignore 2+; really should look at genome dna to decide
            else { $istop= $ix if($istop < 1 or $ix<$istop); }
            }
          # if($istop < 1 and $geneid =~ /AUG/) { $istop= index($oldprot,'X'); }

          # my $prostat = ($p5 and $p3) ? "partial": ($p5)? "partial5" : ($p3)? "partial3" :"complete";
          if($istop > 0 and $istop < length($oldprot)-1) { $astat= "partialinner"; } # innerstop ?
          elsif($astat == 3) { $astat="complete"; }
          elsif($astat == 1) { $astat="partial3"; }
          elsif($astat == 2) { $astat="partial5"; }
          else { $astat="partial"; }
          $astat =~ s/^(.)/o/; # dont confuse w/ new in scans
          $oldstat .= ";oaaln=$oldal,$oldap%,$astat";
          
        } else { $diff++; }

   # FIXME2: exist CDS : allow for problem cases like end-of-scaffold partials
        my $partialAtEndOk= 0;
        for(my $i=0; $i<@oldCDS; $i++) { 
          my $oc= $oldCDS[$i];
          $partialAtEndOk=1 if($oc->[3] < 450); # dont know high end to check
          unless($i<@newCDS and $oc->[3] == $newCDS[$i]->[3] and $oc->[4] == $newCDS[$i]->[4]) { 
            $diff++; $oldstat .= ",NE$i";  ## BUG: clone oc
            if($DEBUG>1) { my @doc= @$oc; $doc[2]="oldCDS"; $doc[0] =~ s/^/#/; push @oldsave, \@doc; }
            }
          }
          
       
        do{ $diff=0; $diffcancel=1; } if(!$NODIFCAN and $partialAtEndOk and $oldstat =~ /artial/ and $newal < $oldal);

   # FIXME4: cancel changes that replace all CDS with completely new CDS, eg for transposon-spans
        $oldprostart5= ($gstrand eq "-") ? $oldCDS[-1]->[4] : $oldCDS[0]->[3];  #above# 
        $oldproend3  = ($gstrand eq "-") ? $oldCDS[0]->[3] : $oldCDS[-1]->[4];  #above# 
        $diff++ if($diff==0 and ($oldproend3 != $proend3 or $oldprostart5 !=  $prostart5));   # off by -1,-2 bugs   

        #  oldc:   [-----]
        #  newc:           [-----]  ; bad = shifted too much
        #  newc:        [-----]     ; bad "
        #  newc:     [---]          ; ok  = not shifted, same stop
        
        if($diff and $oldal>0) { # and $oldstat !~ /artialinner/
          ## .. bugs here ??
          #my($tb,$te)=($prostart5,$proend3); ($tb,$te)= ($te,$tb) if($tb>$te);
          #my($lb,$le)=($oldprostart5,$oldproend3); ($lb,$le)= ($le,$lb) if($lb>$le);
          my($tb,$te)=($newCDS[0]->[3],$newCDS[-1]->[4]); ($tb,$te)= ($te,$tb) if($tb>$te);
          my($lb,$le)=($oldCDS[0]->[3],$oldCDS[-1]->[4]); ($lb,$le)= ($le,$lb) if($lb>$le);
          if($lb==0 or $le==0 or $tb==0 or $te==0) { $tb=$te=$lb=$le=0; } # bug somewhere ...
          
          # argg; location stats not best here, introns have big effect.
          # .. this cancels most span increases...
          if(($tb < $lb and $te < $le) or ($tb > $lb and $te > $le)) { # CDS shifted along exons
            my ($bb,$be)= ( _max($tb,$lb), _min($te,$le) );
            my $maxo= 1 + $be - $bb; # ** not abs
            my $leno= _min( abs($le - $lb), abs($te - $tb)) || 1;
            my $pover= $maxo/$leno; # neg desired: cancel any case of maxo < 0 
            #??do{ $diff=0; $diffcancel=1; } if($pover < 0.33); # cancel change, too little cds similar
            do{ $diff=0; $diffcancel=1; } if(!$NODIFCAN and $pover < 0); # cancel change, too little cds similar
          }
        }

        $cdsattr .=";$oldstat";
        $changed++ if($diff>0);
      } else {
        $changed++; $diff=1;
      }

      if($diff) {
        push @generec, @newCDS; 
        push @generec, @oldsave; 
      } else {
        push @generec, @oldCDS;  
      }
    }
    # else { $diff++; } ## is this error? cdnainIsBest = cdnain best?
    #? else { push @generec, @newCDS;  } # no @$cdsgff no oldCDS
  
    ## add to addattr: count utrs (and span?) and add utrx= if >2
    # change getCDSgff() to return utr-exons?
    # if($u5>2 or $u3>2) { $un=$u5+$u3; $g[0]=~s/$/;utrx=$u5,$u3/; } 

    ## fixme: diff==0, dont remove old attr, add new or not? 
    #     or make option? sometimes always replace? use ocds=eqn,eqaa as flag for no change
    # .. except for some changes above result in NEaa but are cancelled.
    
    if($utrorf) {
      my( $uprot, $ustart, $uend) = orfParts($utrorf);
      my $ulen=  $utrorf->{length};  
      my $ualen= length($uprot);  
      my $upcds  = ($cdnalen>0 && $ulen>0) ? int(100*$ulen/$cdnalen) : 0;
      my $ucompl= $utrorf->{complete};
      $ucompl= ($ucompl==3)?"complete":($ucompl==2)?"partial5":($ucompl==1)?"partial3":"partial";
      $addattr .=";" if($addattr); 
      $addattr .= "aautrlen=$ualen,$upcds%,$ucompl;utroff=$ustart-$uend;utrprot=$uprot";  
    }
    
    if($cdsattr and $diff==0) {
      # check prot, add '*' if missing and complete
      # remove prot, new aalen if  $diffcancel
      unless( $mrnaat =~ m/;protein=/ or $diffcancel) {  ## should set diff=1 for this
        $diff=1; # $mrna->[8] =~ s/$/;$cdsattr/; 
        } ## add if missing
      if($diffcancel) { 
        my ($st)= $oldstat =~ m/ocds=([^;\s]+)/; 
        my ($nal)= $cdsattr =~ m/aalen=([^;\s]+)/;
        $mrna->[8] =~ s/$/;ocds=cancel:$st,nal:$nal/; 
        
      } elsif(not $NoStopCodon and $cdsattr =~ /complete/ and $oldprot =~ /\w/ and $oldprot !~ /\*$/) {
        my $ops= $oldprot.'*';  # this WAS bug for protein=*Maaaa* 
        $mrna->[8] =~ s/;protein=$oldprot/;protein=$ops/; # no end ;
      }
      
    } # elsif()
    
    if($cdsattr and $diff>0) {  
        ## cdsattr keys: cxlen,aalen,protein,utrx
      my @keys= $cdsattr =~ m/(\w+)=/g;
      my $keyp= join('|',@keys);
      $mrna->[8] =~ s/[;]?($keyp)=([^;\n]+)//g; # OR/opt rename keyp_old ?
      $mrna->[8] =~ s/$/;$cdsattr/;
    }
    
    $addattr.=";aaold=$oldaalen" if($oldaalen); ## add always; and $mrna->[8] !~ m/$oldaalen/
    $addattr.=";offold=$oldoffs" if($oldoffs and not($cdsattr=~/cdsoff=$oldoffs/)); ## add ??
    $mrna->[8] =~ s/$/;$addattr/ if($addattr);
    $mrnaat= $mrna->[8]; # update
    
  } else { # no orfprot : can happen, but error if oldCDS << need to restore those here ??
    push @generec, @oldCDS if(@oldCDS);
    my $addattr="ocds=cancel:NEn,aatrans-failed";
    $mrna->[8] =~ s/$/;$addattr/;
    $mrnaat= $mrna->[8]; # update
  }

  return(xxxxxxxxxxx); # generec, addattr changeflags
}    


sub fixgene_cds2prot {
  my($geneid, $generecIN)= @_;
  my($addattr,$changed,$mrnachanged)=("",0,0);

  my @mrna= grep{ $_->[2] eq "mRNA" } @$generecIN;
  my $issplit= (@mrna>1) ? @mrna : 0;
  my $mrna= $mrna[0]; 

  my @oldCDS = grep{ $_->[2] eq "CDS" } @$generecIN;
  my @generec;
  if($issplit) {
    @generec= sort _sortSplitgene grep{ $_->[2] ne "CDS" } @$generecIN;  
  } else {
    @generec= sort _sortgene grep{ $_->[2] ne "CDS" } @$generecIN; # sorts by genome loc, not stranded
  }

  my @exongff= grep{ $_->[2] eq "exon" } @generec;

  #error cases: NOPATH, no @exons
  if(not $mrna or $mrna->[0] =~ /^NOPATH/) {
    $addattr="err=Missing-mrna-location";
    return ( $changed , \@generec, $addattr); 
  }
  unless(@exongff) {
    $addattr="err=Missing-mrna-exon";
    return ( $changed , \@generec, $addattr); 
  }

  my $gstrand= $mrna->[6]; # FIXME for "." ; add CDS bestpro strand (always +??)
  my $mrnaat= $mrna->[8];
  # my $geneid= ($mrnaat =~ m/ID=([^;\s]+)/)? $1 : ""; # make one up?
  my $splitpart= 0;
  ($splitpart,$geneid)= splitPart($mrnaat); # see evigene_idparts ..

  # ... fixgene methods .......
  # my($cdnaoutg,$orfcdsg,$orfprotg)=("") x 3; #output seqs outside dang fixstrand loop

  ## change here, use cdnain to recalc orf protein
  my $oldprot= ($mrnaat =~ m/protein=([^;\s]+)/) ? $1 : "";
  my $oldisbest=($oldprot =~ /\w\w/ and $mrnaat =~ m/(bestaa=pubaa)/)?$1:"";
  my $oldaalen= ## evg variants and gmap.gff stutter
    ($mrnaat =~ m/aalen=(\d+,[^;\s]+)/) ? $1 :
    ($mrnaat =~ m/aaSize=(\d+)/) ? $1 :
    ($mrnaat =~ m/aalen=([^;\s]+)/) ? $1 : ""; # FIXME gmap dup aalen and aaSize<>aalen
  # my $oldoffs= ($mrnaat =~ m/(?:cdsoff|offs)=([^;\s]+)/) ? $1 : "";
  # my($oldcdsb,$oldcdse,$oldcdsor0)= cds_span($geneid,$mrnaat);

  my $revgene=($gstrand eq "-")?1:0; # NOTE exongff are rev order here; 5' is first
  if($issplit) { # sort by parts, by exon?
    @exongff= sort __revSplitGene @exongff;
  } elsif($revgene) {
    @exongff= reverse @exongff; #?? issplit
  }

  my($cdna, $cdnalen, $cdnaexonseq, $xfixfix, $xtrimmed)= (0) x 9;

  if($DO_CDSFIX) { 
    my $xfixed="";
    ($cdna, $cdnalen, $cdnaexonseq, $xfixfix, $xtrimmed, $xfixed)= fixCDS($oldoffs, \@exongff); 
    $addattr .= ";cdsfix=$xfixattr" if($xfixed); #? leave to caller
  }
  
  if($cdnalen == 0) {
    ($cdna, $cdnalen, $cdnaexonseq, $xfixfix, $xtrimmed)= get_splitgenecds($issplit, \@exongff, \$addattr);
  }
  
  my $docomp=$DO_CDSCOMPLETE;
  $docomp= $docomp and ($mrnaat=~/cdsfix=complete/ or ($oldaalen=~/partial/ and $DO_CDSCOMPLETE>1));
  $docomp=0 if($issplit); # too many problems, cdsb,cdse on diff scafs
  if($docomp) {
    # (xxx)= get_competeCDS(xxxx);
    # DO_CDSCOMPLETE completeCDS() should follow SPLITGENEPROTFIX call?
    # 201511. add extend-ends option, ix==0 and ix==$nexon1, look for completing cdsbases of partial cds
    # can do that inside DO_CDSFIX? NO, do separate... with what flags/opts?

    my($extbest,$EXTcdna,$EXTexongff,$extaddattr,$ext5,$ext3)=(0) x 9;
    my $cdnaexons= ($xtrimmed)? $xfixfix : \@exongff;
    my $oldprotPar= ($oldisbest)?$oldprot:"";
    ($extbest,$EXTexongff,$extaddattr,$ext5,$ext3)=
      completeCDSb($geneid, $cdna, $cdnaexons, $revgene, $oldProstart5,$oldProend3,$oldprotPar);
    if($extbest) {
      ## DANG, this is bad for splits, changing exon start/end to wrong split part in EXTexongff
      ($cdna,$cdnalen)= getcdna2( $EXTexongff, 1, 0, 0, CDNA_TRIMNNN);
      @exongff= @$EXTexongff; # gosh darn too many same vars
      $xfixfix= $EXTexongff; # gosh darn too many same vars
      $xtrimmed++; # * REGISTER prior exon change for  FIXMEd, CDS extended NOT exons ..
      $addattr .= ";$extaddattr" if($extaddattr); 
    }  
  } # docomplete
  
  if($xtrimmed) { 
    my($igenerec, $iexongff, $mrnachanged1)= generecfix( \@generec, \@exongff, $xfixfix); 
    @generec= @$igenerec; 
    @exongff= @$iexongff; ## == xfixfix usually
    $mrnachanged += $mrnachanged1;    
    $changed ++; $addattr .=";" if($addattr); $addattr .= "xtrim=$xtrimmed"; 
    $oldStart_b=$oldStart_e=0; # dont use now for prot/cds
  }

  #.........
  # add samecds opt using oldStart > need start,stop of 1st CDS exon
  $KEEPSAMECDS=1; ## PACKAGE VAR; need for oldStart_b
  my($orfprot, $prostart5, $proend3, $bestorf, $utrorf)= 
      getBestProt("partial", $cdna, \@exongff, $oldStart_b,$oldStart_e);
  
  my($xxx) = cdnain_prot($geneid,$geneid,$cdna,$cdnalen,$mrnaat,$oldaalen,$orfprot,\@exongff); #( $geneid, $geneidfix, $cdnalen, $cdnaseq) .. fixme messy
  
  my($yxxx) = add_cdsprot(xxx); # if($orfprot) .. adds new CDS , oldCDS/oldprot checks
        #  ($yxx)= add_cdsprot($mrna, $exongff, $orfprot, $prostart5, $proend3, $cdnalen, $oldCDS);

  # ... end methods ...........
  
  # my($generec,$revbest)= (\@generec,0);
  # if($fixstrand) { # push to other sub? fixgene_strandcds2prot() ??
  #   ($generec,$revbest)= fixstrand( \@genefwd, \@generev);
  #   $changed= ($revbest) ? $changed1 : $changed0;
  # }
  
  $addattr .= ($mrnachanged)?";mrnachanged=$mrnachanged":"";

  # defer _sortgene to caller ?
  # if($issplit) {  @generec= sort _sortSplitgene  @$generec; }
  # } else {  @generec= sort _sortgene  @$generec; }

  return ( $changed , \@generec, $addattr); 

}


sub fixgene
{
  my($actions, $geneid, $generecIN, $geneother)= @_;
  my($changes, $didchange)=(0,0);
  my $changelog="";

  my @mrna = grep{ $_->[2] eq "mRNA" } @$generecIN; # fixme issplit/Splitgene has more mrna parts
  my $issplit= (@mrna>1)? @mrna : 0;
  my @generec;
  if($issplit) { @generec= sort _sortSplitgene  @$generecIN;  } 
  else {  @generec = sort _sortgene  @$generecIN;  }
  my $generec= \@generec;

  # .. fixgene_ actions here .............

  if($actions =~ /cds2prot/) {  
    ($didchange, $generec, $changeflags)=  fixgene_cds2prot($geneid, $generec);
    $changelog.="cds2prot:$changeflags," if($didchange); 
    $changes+= $didchange;
  }
  
  # ... actions end ..........................
  
  if($changes) {  
    my @mrnanew = grep{ $_->[2] eq "mRNA" } @$generec; # must have
    my @exons = grep{ $_->[2] eq "exon" } @$generec;
    if(@mrnanew > 1) { # FIXME.. need to separate exons/part
    
    } else {
      $changes += mrna_fixspan( $mrnanew[0], \@exons);
    }
  }
  
  $changelog ||= "none";
  push @changelog, "$geneid\t$changelog"; # if($DEBUG) 
  warn "#fixg\t$geneid\t$changelog\n" if($DEBUG>1);
  
  my $rowsput= putgene($generec, $geneother, $changes);

  ## separate this ..
  my $seqattr= seqAnnot($mrna); #? add aaold
  putseq($houtcdna,$geneid,$cdnaoutg,$seqattr) if($houtcdna);
  putseq($houtaa,$geneid,$orfprotg,$seqattr) if($houtaa);
  putseq($houtcds,$geneid,$orfcdsg,$seqattr) if($houtcds);
  
  return ($changes) ? 1 : 0;
}


# sub filter_gff_OLD
# {
#   my($inh)= @_;
#   my ($ng,$nr,$nchange)= (0,0,0);
#   my $printpass=1;
#   # $printpass=0 if($actid == ACT_DROP or $actid == ACT_KEEP);
#   my @generec=(); my @otherft=(); my @geneft;
#   
#   ## add own header ?? version?
#   my $version=VERSION;
#   print $houtgff <<"EOGFF";
# ##gff-version 3
# # appl: genefindcds
# # vers: $version
# 
# EOGFF
#   
#   my ($lasplit,$lgid)=(0,0);
#   while(<$inh>){
#     unless(/^\w/){ next if(/^(##gff-ver|#n |$)/);  print $houtgff $_ and next; }
#     my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr)= split"\t";
#     $nr++; chomp($tattr);
#     if($passtypes and "$typ.$src" !~ m/$passtypes/) { print $houtgff $_ if $printpass; next; } 
# 
#     
#     my($gid,$pid); 
#     if($tattr =~ m/\bID=([^;]+)/) {  $gid=$1; }
#     if($tattr =~ m/\bParent=([^;]+)/) {  $pid=$1; 
#       $gid=$pid unless($gid and ($typ =~ /^($RNATYPES)$/)); 
#       }
#     unless(defined $gid) { $gid = "N".$ng; }
#     my $gidfix= $gid; $gidfix =~ s/_C(\d+)$//;
#     
#     my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid]; 
# 
#     if($typ =~ /^($RNATYPES)$/) {  
#       # allow gene and mRNA types .. and/or keep other types in generec
#       
#       #ov2: Split=1/3 .. gene fix here, need to collect all parts (mRNA:1,2,.. exons..) for testgene()
#       my $issplit=($tattr =~ m/;Split=([^;\s]+)/)? $1 : ($gid =~ /_C(\d+)$/)? $1 : 0;
#       if($issplit and $gidfix eq $lgid) {
#         push @generec, $rloc; # add mRNA part2..     
#       
#       } else {
#         $nchange += testgene(\@generec, \@otherft) if(@generec);
#         @generec = ($rloc); $ng++; # maybe best store as array of [gffcols], sorted by type-loc
#         @otherft=(); # ** drops prior gene ft
#         if(@geneft) { unshift @generec, @geneft; @geneft=(); } #?? put into generec ? will cause problems
#       }
#       
#       $lgid= $gidfix; $lasplit= $issplit;
#       
#     } elsif($typ =~ /^($EXONTYPES)$/) {
#       push @generec, $rloc;         
#     } elsif($typ =~ /^($GENETYPES)$/) {
#       @geneft =( $rloc ); # keep only latest
#       
#     } elsif($tb>0 and $te>0) {
#       push @otherft, $rloc;         
#     }
#       
#   }
#   
#   # FIXME lost 1st gene row before mrna; FIXME change gene span w/ mrna span changes
#   
#   $nchange += testgene(\@generec, \@otherft) if(@generec);
#   
#   return ($nchange,$ng);
# }

sub filter_gff
{
  my($inh)= @_;
  my ($ng,$nx,$nr,$nfix,$nhit,$errord)= (0) x 10;
  my $nocomm= 1; ##($actid == ACT_NULL) ? 1 : 0;
  my @generec=();
  my @geneother=();
  my $geneid=""; my $lgeneid="";
  my $generow=undef;
  
  # while(<$inh>) 
  while(<>) { # perl special, STDIN or @ARGV
    unless(/^\w/){ print unless($nocomm or /^(#n |$)/); next; }
    my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr)= split"\t";
    $nr++; chomp($tattr);
   
    my($gid,$pid); 
    if($tattr =~ m/\bID=([^;]+)/) {  $gid=$1; }
    if($tattr =~ m/\bParent=([^;]+)/) {  $pid=$1; 
      $gid=$pid unless($gid and ($typ =~ /^($RNATYPES)$/)); 
      }
    unless(defined $gid) { $gid = "N".$ng; }

    # my $issplit= ($gid =~ m/_C(\d+)/)?$1 : (m/;Split=(\d+)/)?$1:0;
    my($issplit,$gidsp)= splitPart($tattr);

    if($issplit) {
      $gid =~ s/_C\d+//; #?? maybe not here, want split-part ID to process each exon subset
      $pid =~ s/_C\d+//; 
    }
    
    my $rloc= [$ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr,$gid,$issplit,$gidsp]; 

    if($typ =~ /^gene$/) { $generow= $rloc; }
    elsif($typ =~ /^($RNATYPES)$/) {  
      # allow gene and mRNA types .. and/or keep other types in generec
      
      if($gid eq $geneid and $issplit) {
        push @generec, $rloc; $ng++; # check Parent == $geneid ?
      
      } else {
      
        $nfix += fixgene($actions, $geneid, \@generec, \@geneother) if(@generec);
        # ($changes, $changelog, $rowsput)=  fixgene($actions, $geneid, \@generec, \@geneother)
        
        $ng++;
        $geneid= $gid;  # parse for gene vs tr id/alttr num ?
        @generec = ($rloc); # maybe best store as array of [gffcols], sorted by type-loc
        @geneother= (); 
        $generow=undef;
      }
      
    } elsif($typ =~ /^($EXONTYPES)$/) {
      if($pid ne $geneid) { warn "#ERR: Out-of-order GFF $typ:$pid in mRNA:$geneid\n"; $errord++; next; }
      push @generec, $rloc; $nx++; # check Parent == $geneid ?
      
    } elsif( ! $passtypes or "$typ.$src" =~ m/$passtypes/ ) {
      push @geneother, $rloc; # check Parent == $geneid ?
    }
    
  }
  
  $nfix += fixgene($actions, $geneid, \@generec, \@geneother) if(@generec);
  return ($ng,$nx,$nfix);
}


sub getcdna2 { 
  my($exons, $asexons, $expand, $expend, $trimNNN)= @_;  # add flag: trimNNNends. drop expand,expend ..
  my $xchange=0;
  my $nexons=  @$exons;
  my $nx1= $nexons - 1;
  # @$exons; _sortgene, always start<stop ?? maybe bug here in sort
  # NO:  my @sexon= sort _sortloc @$exons;
  
  my @j= (0..$nx1); 
  my $je= pop(@j); unshift @j, $je;
  my @exok= (0) x $nexons;
  my @exseq= ("") x $nexons;
  my $nx0= 0;
  
  for my $j (@j) {
    my $ft= $exons->[$j]; # $sexon[$j]; ??
    my($exondna,$rev,$xlen)= getexondna($ft,0); # DOES NOT revcomp(dna) if rev; DOREVCOMP does it
    $exok[$j]= $ft; $exseq[$j]=$exondna;  
    if($trimNNN) {
      my $whichj= ($j == $nx1) ? "last" : ($j == $nx0) ? "first" : "middle";
      # not bug: trimcdnaexon expects rev/fwd always w/ fwd(dna)
      my($didtrim,$texondna,$tft,$tstart,$tstop)= trimcdnaexon($j,$exondna,$ft,$rev,$nx1,$whichj);
      if($didtrim) {
        $exondna=$texondna; $ft=$tft; $xchange++;
        if($didtrim == kEXONDROP) { 
          $exok[$j]=0; $exseq[$j]="";
          if($j==$nx0) { $nx0++; }  
          elsif($j==$nx1) { }
          $nx1--;
        } else {      
          $exok[$j]= $tft; $exseq[$j]= $texondna;  
        }
      }
    }
  }

  my(@exonft,@exonseq);
  my $cdna= ""; my $cdnalen=0; 
  for my $j (0..$nexons-1) {
    my $exondna= $exseq[$j] or next;
    my $exft= $exok[$j];
    my $rev= ($exft->[6] eq "-")?1:0;
    $exondna = revcomp($exondna) if($rev);  # $exondna = reverse $exondna;  $exondna =~ tr/ACGTacgt/TGCAtgca/;
    $cdna .= $exondna;  
    push @exonseq, $exondna;
    push @exonft, $exft;
  }
  
  $cdnalen= length($cdna);
  return($cdna,$cdnalen,\@exonseq,\@exonft,$xchange);
}

sub getexondna {
  my($exonft,$dorevcomp)= @_; 
  my($ref,$start,$stop,$strand,$phase,$xgid)= @{$exonft}[0,3,4,6,7,9];
  my $isrev=($strand eq '-')?1:0;
  my $exondna  = get_dna( $chrfasta, $ref, $start, $stop);
  $exondna ||="";
  if($DEBUG and not $exondna) {
    warn "#getexondna=0 gid=$xgid loc=$ref:$start-$stop\n"; # loc=NOPATH:1-69
  }
  $exondna = uc($exondna); # always?
  my $xlen = length($exondna); # 1+$stop-$start; # add even if no dna
  if($dorevcomp and $isrev) {  $exondna=revcomp($exondna); }
  return ($exondna,$isrev,$xlen);
}

sub revcomp { 
  my($dna)=@_; 
  $dna = reverse $dna; $dna =~ tr/gatcGATC/ctagCTAG/;
  return($dna);
}

sub get_dna {
  my($fasta, $ref, $start, $stop)= @_; #, $fasta_db_ref
  unless( $ref && $stop>0) { warn "need ref-ID, start, stop location\n"; return; }
 
  require Bio::DB::Fasta;  # FIXME: not portable w/o parts of BioPerl !
  my $havedb= ref $fasta_db;
  unless($havedb) {
    my $db = eval { Bio::DB::Fasta->new($fasta); } or die "$@\nCan't open sequence file(s). "; # return ?
    $fasta_db= $db;  
    }
  
  my $seq = $fasta_db->seq($ref, $start => $stop)  or return;  ## die "cant locate seq of $ref:$start-$stop\n"; 
  $seq= $seq->seq if(ref $seq); #  weird bioperl change  
  return $seq;
}

##------------------------------------------------------
## non-gff gene info methods from evgs/prot/trclass2mainalt.pl.. to  other pm?

our $IDPREFIX= $ENV{idprefix} || 'EVGm'; 
our $SHOWDROPS=0; 
our $SIZESORT=1; 
my ($trclass,$output,$logfile);  
my ($trpath,$trname, $sradatah, %settings);
our ($pubid_format,$altid_format,$GBPROID);
our $pubidnum_start=0;
#? our $CULLXEQ=0;
our $preserveOldIds=0; # change to ok pattern, IDPREOK

#upd1708: preserveOldIds option ..  keep input gene ids as best feasible, including altnum?
#  need all alts per main, check evg gene ids for 1 x many, also need check across loci for 2+ new loci w/ old geneids
#  where cant preserve geneid for all of locus, do what? new geneid, or modify old?

sub evigene_idparts {  
  my($id)=@_;
  my($gpre,$gnum,$gd,$ti,$gdup,$isplit)=(0) x 9;
  # idparts: idprefix, genenum, geneid, altnum, mapdupnum, mapsplitnum
  # MyorgEVm000001t9_G2_C2 => MyorgEVm, 000001, t9, _G2, _C2
  # shouldnt have both _G and _C map suffices
  ## FIXME for _G2 _C2? tags, G2,n is diff locus id
  
  $gd=$id; 
  ($isplit)= ($gd=~s/_(C\d+)$//)?$1:0;
  ($gdup) = ($gd=~s/_(G\d+)$//)?$1:0;
  ($ti)   = ($gd=~s/t(\d+)$//)?$1:0; # other patts?  '\d([a-z])(\d+)$' 
  ($gnum) = ($gd=~m/(\d+)$/) ? $1:0;
  ($gpre=$gd) =~ s/$gnum$//;
  if($gdup) { $gd.=$gdup; $gnum=0; } # gnum invalid?
  #orig: return($gpre,$gnum,$gd,$ti,$isplit); #? return gdup as part? not gd?
  return($gpre,$gnum,$ti,$isplit,$gdup); #?
}

sub evigene_idOfparts {  
  my($gpre,$gnum,$ti,$isplit,$gdup)= @_; # from evigene_idparts
  $ti||=1;
  $gnum= sprintf( '%06d',$gnum); # leading 000
  my $id= $gpre . $gnum . "t$ti"; 
  $id.="_C$isplit" if($isplit); $id.="_G$gdup" if($gdup);
  return($id);
}

# sub preserveOldIds {
#   my($mainlist, $altsOfMain, $drop, $altsize)= @_;
#   my(%gids, %gnums, %gdone, %newids, $gprefix);
#   my $nids=0;
#   my $NEWIDpre='n';
#   
#   ## prefix, gnum problem: have mixed gprefices .. diff gnum set for each
#   my $GNEXTNUM=0;
#   foreach my $md (@$mainlist) {
#     my @okd = grep{ not($drop->{$_}) } ($md, keys %{$altsOfMain->{$md}} );
#     for my $id (@okd) {
#       my($gpre,$gnum,$gd,$ti,$isplit)= evigene_idparts($id);
#       next unless($gnum);
#       $gnums{$gnum}++; # $gids{$gd}{$ti}= $id;
#       $gprefix= $gpre unless($gprefix);
#       }
#   }
#   my($glast)= sort{ $b <=> $a } keys %gnums;
#   $GNEXTNUM= 9 + $glast;
#   
#   my $idformat= $NEWIDpre . $gprefix . '%06d';
#   my(%havepubg);
# 
#   foreach my $md (@$mainlist) {
#   
#     # my @ad= sort{ $altsOfMain->{$md}{$a} cmp $altsOfMain->{$md}{$b} # class sort
#     #  or $altsize->{$b} <=> $altsize->{$a} or $a cmp $b } keys %{$altsOfMain->{$md}}; 
#     my @ad= sort{ $a cmp $b } keys %{$altsOfMain->{$md}}; 
#     my @okd = grep{ not($drop->{$_}) } ($md,@ad);
#   
#     %gnums=(); %gids=();
#     my $gnumfirst=0; my $gprefat=0;
#     for my $id (@okd) {
#       my($gpre,$gnum,$gd,$ti,$isplit)= evigene_idparts($id);
#       next unless($gnum); # for _G2,n.. new loci
#       $gnums{$gnum}++; $gids{$gnum}{$ti}= $id;
#       $gnumfirst=$gnum unless($gnumfirst or $gdone{$gnum}); 
#       $gprefat= $gpre unless($gprefat); #?? problems? yes
#     }
#     
#     my @gnums= sort{ $gnums{$b}<=>$gnums{$a}  or $a <=> $b } keys %gnums;
#     unless($gnumfirst) { # pick most or first in (main) <<
#       for(my $i=0; $i<=$#gnums; $i++) {
#         unless($gdone{$gnums[$i]}) { $gnumfirst=$gnums[$i]; last; }
#       }
#     }
#     # usage assumption: -keepids=IdPreA -idpre IdPreA, so inval gprefat will go to ++GNEXTNUM 
#     unless($gprefat and $gprefat =~ /$IDPREOK/){ $gprefat=$IDPREFIX;  $gnumfirst=0; }
#     $idformat= $NEWIDpre . $gprefat . '%06d'; # maybe new for each main id??
#     my($pubgn);
#     do { 
#       unless($gnumfirst) { $gnumfirst= ++$GNEXTNUM; } #??
#       $pubgn = sprintf( $idformat, $gnumfirst);
#       if($havepubg{$pubgn}) { $gnumfirst=0; }
#     } while ($havepubg{$pubgn});
#     $havepubg{$pubgn}=1; 
# 
#     # try2
#     my %tidone=(); my $timax= @okd; #? not same using orig ti > @okd
#     for my $id (@okd) {
#       my($gpre,$gnum,$gd,$ti,$isplit)= evigene_idparts($id);
#       my $tinew=0;
#       if($gnum eq $gnumfirst and not $tidone{$ti}) { $tinew=$ti; }
#       if($tinew==0) { do { $tinew++; } while( $tidone{$tinew} ); } # can put lowqual alts at ti top
#       $tidone{$tinew}++; 
#       $gdone{$gnumfirst}++;
#       #above# my $pubgn = sprintf( $idformat, $gnumfirst); 
#       my $pubti = sprintf( $altid_format, $tinew); # same ti or new? cant use same ti w/o checks
#       my $pubid=  $pubgn . $pubti;
#       $newids{$id}= $pubid;
#       $nids++;
#     }
#   } # mainlist
#   
#   return($nids,\%newids);
# } # sub preserveOldIds
# 
# 
# sub make_pubid # add preseveOld %newids here? use altnum w/ it or preserve old altnum?
# {
#   my($oid, $pubidnum, $altnum, $preservedIds)= @_;
#   $pubidnum_start= $pubidnum; #?
# 
#   my($pubid,$pubgene);
#   if($preservedIds and ref($preservedIds)) {
#     if(my $pid= $preservedIds->{$oid}) {
#       if(my($pgene,$palt)= $pid=~m/(\w+)t(\d+)$/) {
#         $pubgene=$pgene;
#         # if($palt ne $altnum) ..
#         $pubid   = $pubgene . sprintf( $altid_format, $altnum); # or palt ?
#         return($pubid,$pubgene);
#       }
#     }
#   }
#   $pubgene = sprintf( $pubid_format, $pubidnum); 
#   $pubid   = $pubgene . sprintf( $altid_format, $altnum);
#   return($pubid,$pubgene);
# }

sub make_IDPREFIX
{
  my($existingID)= @_;
  my $digits=6; # or 7?
  my $ChangeDefaultIdPrefix= 0; # ($IDPREFIX eq $DEFAULT_SETTINGS{'IDPREFIX'}) ? 1:0;
  
  if($existingID and $existingID !~ m/^$IDPREFIX/) {
    $existingID=~ s/t\d+$//;
    my($prefix,$nums)= $existingID =~ m/^(\w+)(\d\d\d+)$/; ## prefix may have numbers, how many? 
    if($prefix) { $IDPREFIX= $prefix; $digits= length($nums) || 6; }
    $ChangeDefaultIdPrefix=0;
  }

  my $nd= ( $IDPREFIX =~ s/(0\d+)$// ) ? length($1) : $digits;
  my $pubid_format = $IDPREFIX.'%0'.$nd.'d';  
  my $altid_format = 't%d'; 
  return($pubid_format,$altid_format);  
}


1;

__END__