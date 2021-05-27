# cdna_proteins.pm

# package cdna_proteins;
package main;

use strict;
use warnings;

#? require Exporter;
# our @ISA = qw (Exporter);
# our @EXPORT = qw (translate_sequence get_protein reverse_complement);

use vars qw ($DEBUG $MINAA $MINEXON $MINGOOD $USEGOODLEN 
      $AA_cdna_GT_genome $ORF_FULLvPART $KEEPSAMECDS $NoStopCodon $InnerStopToX
      $pCDSbad $pCDSpoor $MINUTRORF  $USESelenocysteine $TRIMCDSENDGAP
      );

use constant cdna_proteins_VERSION  => '20141231'; # Selc annots added
# '20140317';  # '20131124' 20130818 0228'20120721'; 

our $DEBUG=0;
our $MINAA= 30;  # used 
our $MINEXON= 60; # for cut
our $MINGOOD= 0.75; # filter out prots w/ fewer good aminos
# our $USE_CDSEXONS = 0;
our $USEGOODLEN=1;
our $AA_cdna_GT_genome= 1.50;  # for prot(cdna) > prot(genome) test; option?
# FIXME: adjust when to take partial vs complete, eg partial5 often is a few aa longer than complete M
our $ORF_FULLvPART = 0.85;
our $KEEPSAMECDS= 0; # prefer keep same CDS exons but can extend/shorten protein bounds
our $NoStopCodon=0; 
our $InnerStopToX=0;  # 3 values?  0=leave *, 1= tr/*/X/, -1 or 2=skip orf
    ## need global params for utrbad/poor
our $pCDSbad = $ENV{pcdsbad}  ||30; # adjust DOWN *?
our $pCDSpoor= $ENV{pcdspoor} ||60;
our $MINUTR= $ENV{minutr}||300; # ~ 300b "fixed" average utr sizes, maybe too low, 
#our $BAD_GAPS= $ENV{aagapmax} || $ENV{BAD_GAPS} || 15;  # % gaps in AA

our $MINUTRORF= $ENV{minutrorf}||300; # was 300; #?
our $TRIMCDSENDGAP= 1; #201503 add; some bugs here

# our @stop_codons = qw(TAA TAG TGA); # allow changes, esp TGA => SelC
# parts from PASA/PasaLib/Nuc_translater.pm 
use vars qw ($currentCode %codon_table @stop_codons %amino2codon_table);
use constant backtrans_AmbiguousNucs => 0;  # default off
our (%amino2codon_ambig,%amino2codon_fixed); # see backtranslate_protein & BEGIN

# convert to array handling? _max(@list) 
sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }


sub bestorf_test
{
  my($bestorf,$nextorf) = @_;
  my ($bestprot)= orfParts($bestorf);
  my ($nextprot)= orfParts($nextorf);
  my $MinCDS= 3*$MINAA;

       # FIXME bestorf_test: option too high -ratiocdna 1.25; at least should replace cdna XXXX gaps with cdnain perfect
       # test for near-sameprot but for XXX gaps? 
  
  my $oki= ( $nextprot =~ /\w/ 
    && ($nextorf->{goodlen} >= $MinCDS)  # dang, goodlen,leng are cdslen not aalen
    && ($nextorf->{goodlen}/$nextorf->{length} >= $MINGOOD)) ? 1 : 0;
    
  if($oki) {
    if( $bestprot ) {
      ## problem here for rev=0,1 missing huge partial rev for short fwd .. but get huge if just rev
      my $arat = $nextorf->{goodlen} /  $bestorf->{goodlen}; # ok,> MINAA
      my $adiff= $nextorf->{goodlen} -  $bestorf->{goodlen};  

      if( $arat > 1.0 and  $arat < $AA_cdna_GT_genome and $bestprot =~ m/XX/) { ##  $bestorf->{goodlen}/$bestorf->{length} < 0.99
        (my $bp= $bestprot) =~ s/X/./g; 
        # if(length($bp) > length($nextprot) { } # chomp some?
        return(1, $nextprot, $nextorf) if($nextprot =~ m/$bp/);
      }
      
      if( $arat > $AA_cdna_GT_genome
        or ( $adiff >  0 and ($nextorf->{complete}==3 or $bestorf->{complete} <3) ) 
        or ( $adiff >= 0 and $nextorf->{complete}==3 and $bestorf->{complete} <3 ) )
       { 
        return(1, $nextprot, $nextorf);
       }
    } else {
      return(1, $nextprot, $nextorf);
    }
  }
  return( 0, $bestprot, $bestorf);
}

sub utrorf_test
{
  my($bestorf, $allorfs, $cdnasize) = @_;

  my ($utrorf,$utrosize)= getUtrOrf($bestorf, $allorfs, $cdnasize);  
  if ($utrorf) {
    my ($orfok) = bestorf_test(undef,$utrorf);
    return($utrorf,$utrosize) if($orfok);  
  }
  return(0,0);  
}

sub revorf_report
{
  my($bestorf, $allorfs, $cdnasize, $issorted) = @_;
  my ($revorf,$revosize)= getRevOrf($bestorf, $allorfs, $cdnasize, $issorted);  
  if ($revorf) {
    my($aalen,$pcds,$compl,$orflen,$fahead)  
     	= proteindoc($revorf,$cdnasize); # ,$strand,0
		 ## note: fahead= "aalen=99,90%,complete; clen=350; strand=-; offs=9-309;";
		 ## aarev=99,90%,complete,strand:-,offs:9-309;
		 ## aarev=99,90%,complete,-,9-309;  << use this?  or aarev=99,...,9-309:-; 
		 ## add pct-of-longest to report?  maybe drop strand,offset to simplify
		 ## aarev=55%,aa99,90%,complete,s-,o9-309;  << use this?  or aarev=99,...,9-309:-; 
		my $longsize= $bestorf->{goodlen} || 1; # revosize == goodlen; orflen == length
		my $plong= int(0.5 + 100*$revosize / $longsize);
		$fahead =~ s/aalen=/aa/; $fahead =~ s/clen=d+;//; 
		$fahead =~ s/strand=/s/; $fahead =~ s/offs=/o/; 
		$fahead =~ s/=/:/g; $fahead =~ s/ //g; $fahead =~ s/;\s*$//;  $fahead =~ s/;/,/g; 
		return ($aalen,"aarev=$plong%,$fahead"); #?
	}
  return(0,0);  
}

sub proteindoc
{
  my($orf, $cdnalen, $cdnarev, $forFASTA) = @_;
  my($aalen, $compl, $pcds, $istop, $Selc) = (0) x 10;
  
  my( $orfprot, $prostart, $proend, $orflen, $orient) = orfParts($orf); # [qw(protein start stop length orient cdnalen complete)]
  $cdnarev ||= $orient; # fill in blank; use always? shouldnt need to pass cdnarev as param.
    # Ooops, orf->{start,stop} are reversed for -strand; start>stop : ($lb,$le)=($le,$lb) if($lb>$le);
  if($orfprot) {
    $pcds  = ($cdnalen>0 && $orflen>0) ? int(100*$orflen/$cdnalen) : 0;
    my $urev= ($prostart>$proend)?1:0;
    my $u1len= ($urev) ? $cdnalen - $prostart : $prostart - 1; 
    my $u2len= ($urev) ? $proend - 1 : $cdnalen - $proend;
    
    $aalen= length($orfprot); # $bestorf->{length}; # this is cds-len
    if(substr($orfprot,-1,1) eq '*') { $aalen--; if($NoStopCodon) { $orfprot =~ s/\*$//; } }
    $istop= $orf->{innerstop} || 0; # add 201403:
    $compl= $orf->{complete}; #? trust? or check ^M..*$
    $compl= ($compl==3)?"complete":($compl==2)?"partial5":($compl==1)?"partial3":"partial";
    
    ##? not bad if partial? if u1len or u2len == 0
    ## need global params for utrbad/poor
    if($cdnalen - $orflen <= $MINUTR) { } # ignore pcds if utr small
    elsif($pcds < $pCDSbad or $u1len > $orflen or $u2len > $orflen) { $compl.="-utrbad"; }  
    elsif($pcds < $pCDSpoor) { $compl.="-utrpoor";  } #?? maybe change to use prostart OR protend3 > 35%? 40% ?
    ##? add istop flag to compl ??

    # 2014.12 add Selc flags 
    # Funhe2Exx11m009903t6 aalen=556,80%,complete,selcstop; Selcstop=index; .. Name=Selenoprotein N 
    if($USESelenocysteine and index($orfprot,'u')>0) { 
      $compl.=",selcstop"; 
      my $isel= $orf->{Selc}||1; $Selc="Selcstop=$isel"; # want mRNA/CDS index * 1-origin,may be list: 123,456,888
      ## NCBI needs mRNA position, as range (3b codon), .. do any have >1 u, yes
      ## CDS  /transl_except=(pos:1002..1004,aa:Sec);
      # >Funhe2Exx11m009903t6 aalen=556,80%,complete; Selcstop=1; clen=2067; strand=+; offs=128-1798; pubid=Funhe2EKm010847t1; oid=Funhe2EKm010847t1,Funhe2Emap3m010942t3; organism=Fundulus_heteroclitus; type=protein; isoform=1; Name=Selenoprotein N (100%P); genegroup=FISH11G_G7797; Dbxref=TrEMBL:UniRef50_Q9NZV5,TrEMBL:SELN_HUMAN,; tblerr=E:InternalStop,E:MisMatchAA,
      # GALDDQSC>>u<<GSGRTLRETVLESSPVLALLNQSFVSSWSLVRELENMQADEENPALSEKAR
      ## 2+ u cases: Funhe2Exx11m049881t1/Funhe2EKm029327t1 SelM, aalen=135,99%,partial; Selcstop=1; tho 2nd 'uu' may be real stop.
      ## Funhe2Exx11m129535t7/Funhe2EKm029571t1 aalen=321,92%,complete; Selcstop=1; Name=SelP; 
    }
    
    if($forFASTA) { $orfprot =~ s/(.{60})/$1\n/g; } #? only for fasta output
  } else {
    $orfprot="X"; # ? make dummy orfprot?
  }
  
  ## fixme: add goodlen or gaps count to aadoc    
  # my $aagap= $orf->{length} - $orf->{goodlen};
  # my $aagap= $orfprot =~ tr/X/X/;
  
  my $fahead= "aalen=$aalen,$pcds%,$compl; clen=$cdnalen; strand=$cdnarev; offs=$prostart-$proend;";
  $fahead .= " $Selc;" if($Selc);
  $fahead .= " innerstop=$istop;" if($istop>0);
  if(my $orflags= $orf->{flags}) { $fahead .= " orflags=$orflags;"; }
  return($aalen,$pcds,$compl,$orflen,$fahead,$orfprot);
}

sub proteinqual # short version proteindoc
{
  my($orf, $cdnalen) = @_;
  $cdnalen||=0;
  my( $orfprot, $prostart, $proend, $orflen, $orient, $ocdnalen, $complete, $istop) 
    = orfParts($orf, [qw(protein start stop length orient cdnalen complete innerstop)]);
 
  $cdnalen= $ocdnalen unless($cdnalen > $ocdnalen);

  my $pcds = ($cdnalen>0 && $orflen>0) ? int(100*$orflen/$cdnalen) : 0;
  my $urev= ($prostart>$proend)?1:0; # or $orient eq '-' ??
  my $u1len= ($urev) ? $cdnalen - $prostart : $prostart - 1; 
  my $u2len= ($urev) ? $proend - 1 : $cdnalen - $proend;
    # $istop= $orf->{innerstop} || 0; # add 201403:
  my $compl= ($complete==3)?"complete":($complete==2)?"partial5":($complete==1)?"partial3":"partial";
  if($cdnalen - $orflen <= $MINUTR) { } # ignore pcds if utr small
  elsif($pcds < $pCDSbad or $u1len > $orflen or $u2len > $orflen) { $compl.="-utrbad"; }  
  elsif($pcds < $pCDSpoor) { $compl.="-utrpoor";  } #?? maybe change to use prostart OR protend3 > 35%? 40% ?
  my $cdsoff="$prostart-$proend"; $cdsoff.=":$orient" if($orient);
  
  my $Selc=0;
  my $aalen= int($orflen/3); # length($orfprot); 
  if($orfprot) {
    $aalen= length($orfprot); # $bestorf->{length}; # this is cds-len
    if(substr($orfprot,-1,1) eq '*') { $aalen--; if($NoStopCodon) { $orfprot =~ s/\*$//; } }
    # Funhe2Exx11m009903t6 aalen=556,80%,complete,selcstop; Selcstop=index; .. Name=Selenoprotein N 
    if($USESelenocysteine and index($orfprot,'u')>0) { 
      $compl.=",selcstop"; 
      my $isel= $orf->{Selc}||1; $Selc="Selcstop=$isel"; # want mRNA/CDS index * 1-origin,may be list: 123,456,888
      ## NCBI needs mRNA position, as range (3b codon), .. do any have >1 u, yes
      ## CDS  /transl_except=(pos:1002..1004,aa:Sec);
    }
  }

  $compl.=",innerstop$istop" if($istop>0);
  if(wantarray) {
    return($aalen,$pcds,$compl,$orflen,$cdnalen,$cdsoff,$Selc);# ,$orfprot
  } else {
    my $fahead= "aalen=$aalen,$pcds%,$compl;cxlen=$orflen/$cdnalen;cdsoff=$cdsoff;";
    $fahead .= "$Selc;" if($Selc);
    return $fahead;
  }
}


sub getBestProt2
{
  my($ptype, $cdna, $exongff, $oldStart_b, $oldStart_e)= @_;
  my($orfprot,$prostart5,$proend3)=("",0,0,"");
  $oldStart_b||=0; $oldStart_e||=0;
  
  # fix this to return longest full prot, and longest partial (if longer)
  # .. test which is best.
  # FIX2: add test intron overlaps : retained = intron inside exon; err = intron rev at splice/inside
  #   my $longorf= $longest_orf_finder->get_longest_orf($cdna); # == hash
  # FIXME3: option to mark/return 2ndary orf(s) in aberrant long-utr transcripts, that dont overlap 1st orf
  
  # my ($longorf,$longfull,$orfs) = getAllOrfs($cdna); 
  
  my $strands= ($ptype =~ /(both|fwd|for|rev)/i)? $1 : "fwd";
  my ($longorf,$longfull,$orfs) = getAllOrfs($cdna, $strands, $ptype);  
  
  ## 201402 option: ($longorf)= getOneOrf($cstart,$cend,$seq,$strand); # return (lor,for,@orfs)?
  
  #  where strands=null/fwd/rev/both, and strands=rev means do revcomp on cdna
  # Ooops, orf->{start,stop} are reversed for strand=rev; start>stop
  
  if(ref($longorf)) {
    # $orfprot= $longest_orf_finder->get_peptide_sequence();
    # ($prostart5,$proend3)= $longest_orf_finder->get_end5_end3();

    my $lookmore= ($ptype =~ /long/)?0:1; 
    $lookmore=0 if($longorf->{flags} && $longorf->{flags}=~/bestspan/);# 201402 orthobest; fix2: return not bestspan..
    
    if($KEEPSAMECDS and $oldStart_b > 0) { 
      # may not be right yet; try this?
      # { $_->{start} < $oldStart_e && $_->{stop} > $oldStart_b }
      my($samestartorf);
      if($ORF_FULLvPART <= 0.8) {
        ## oldStart_e is end point of orf.start NOT orf.stop
      #BUG:($samestartorf) = grep { $_->{complete} == 3 and $_->{start} >= $oldStart_b and  $_->{stop} <= $oldStart_e } @$orfs;
      ($samestartorf) = grep { $_->{complete} == 3 and $_->{start} >= $oldStart_b and  $_->{start} <= $oldStart_e } @$orfs;
      } else {
      ($samestartorf) = grep { $_->{start} >= $oldStart_b and  $_->{start} <= $oldStart_e } @$orfs;
      #BUG:($samestartorf) = grep { $_->{start} >= $oldStart_b and  $_->{stop} <= $oldStart_e } @$orfs;
     }
      if(ref $samestartorf) { $longorf= $samestartorf; $lookmore=0; } #NOT# else { return (undef); } # not found == no change, here
    } 
    
    if($lookmore and $longorf->{complete} < 3 and ref($longfull) ) {
      my $keylen=($USEGOODLEN)?"goodlen":"length";
      my $lsize= $longorf->{$keylen}; 
      my $fsize= $longfull->{$keylen};
      
      # FIXME: adjust when to take partial vs complete, eg partial5 often is a few aa longer than complete M
      if($fsize >= $ORF_FULLvPART * $lsize) {
        if(ref($exongff)) {
        my ($cdslong, $attrL, $pcodeL, $maxutrL)= getCDSgff2( $exongff, $longorf);
        my ($cdsfull, $attrF, $pcodeF, $maxutrF)= getCDSgff2( $exongff, $longfull);
        $longorf= $longfull if( $maxutrF < 3 or ($pcodeF >= $ORF_FULLvPART * $pcodeL));
        #?? save/return cds or just cds.begin-end ?
        } else {
	      $longorf= $longfull;
	      }
      }
    }
    
    #old# ($orfprot,$prostart5,$proend3)= orfParts($longorf);
  }

  return($longorf, $orfs); # return allorfs now
  # old# return($orfprot, $prostart5, $proend3, $longorf, $orfs); ##, $utrorf 
}

sub getBestProt # old 
{
  # my($ptype, $cdna, $exongff, $oldStart_b, $oldStart_e)= @_;
  my ($longorf, $orfs)= getBestProt2( @_ );
  my ($orfprot,$prostart5,$proend3)= orfParts($longorf);
    # my $getutrorf=1; # always test for utr protein? moved out of here
    # my ($utrorf,$utrosize)= getUtrOrf($longorf, $orfs, length($cdna)); # ($getutrorf) ? xxx() : (0,0);
    # ^^ remove from here, add to utrorf_test
  return($orfprot, $prostart5, $proend3, $longorf); ##, $utrorf 
}

sub getUtrOrf
{
  my ($longorf, $orfs, $cdnalen)= @_;
  my ($utrorf,$utrosize)=(0,0);
  my $lsize= $longorf->{length};
  # my $lgood= $longorf->{goodlen};
  return($utrorf,$utrosize) unless(($cdnalen - $lsize >= $MINUTRORF));  # test even if lsize/cdna > 60% ?  
  # my $dotest= (($cdnalen - $lsize > $MINUTRORF) || ((100*$lsize/$cdnalen) < $pCDSpoor));
  #  ^^ this is bad test; some large tr with pCDS > pCDSpoor are hiding other orf, usually in 3+utr exons on 1 end
  if(1) { 
    # Ooops, orf->{start,stop} are reversed for -strand; start>stop : ($lb,$le)=($le,$lb) if($lb>$le);
    my($lb,$le)= ($longorf->{start},$longorf->{stop});  ($lb,$le)=($le,$lb) if($lb>$le);
    foreach my $orf (@$orfs) {  # orfs can be rev of longorf, start/stop are fwd tho
      my($ob,$oe,$ogood,$osize)= ($orf->{start},$orf->{stop},$orf->{goodlen},$orf->{length},);
      ($ob,$oe)=($oe,$ob) if($ob>$oe);
      if(($ob > $le or $oe < $lb) and ($ogood>=$MINUTRORF) and $ogood>$utrosize) { # size or $osize > 0.5*$lsize ??
        $utrorf= $orf; $utrosize= $ogood;        
        }
      }
    }
  
  # 201402 fixme: do splitutrorf mrnaseq here, see also introncut mrnaseq ..
  # longorf and utrorf should both have ->{mrnaseq} with proper portions cut..
  # this also changes quality utrbad, etc. and output of mrnaseq
    
  return($utrorf,$utrosize);    
}

sub getRevOrf # or revorf_test  # find/report on longest/best orf in reverse of longorf
{
  my ($longorf, $orfs, $cdnalen, $issorted)= @_;
  my $MinCDS= 3*$MINAA; $issorted||=0;
  my ($revorf,$revosize)=(0,0);
  my $lorient= $longorf->{orient};
  # my $lsize= $longorf->{length}; # goodlen ?
	##  orf->{start,stop} are reversed for -strand; start>stop : ($lb,$le)=($le,$lb) if($lb>$le);
	## my($lb,$le)= ($longorf->{start},$longorf->{stop});  ($lb,$le)=($le,$lb) if($lb>$le);
	foreach my $orf (@$orfs) {  # expect long sorted, or check all? 
		my($oorient,$ogood,$osize)= ($orf->{orient},$orf->{goodlen},$orf->{length});
		# $ob,$oe, = $orf->{start},$orf->{stop}, # ($ob,$oe,$oorient)=($oe,$ob,'-') if($ob>$oe);
		if($oorient ne $lorient and $ogood>=$MinCDS and $ogood > $revosize) {  
			$revorf= $orf; $revosize= $ogood; 
			last if $issorted;
			}
   }  
  return($revorf,$revosize);    
}

# replace old sub getCDSgff($exons,$orfprot,$prostart5,$proend3) w/ getCDSgff2($exons,$orf)
sub getCDSgff2 # cdna_proteins.pm
{
  my($exons,$orf)= @_;

  my ($orfprot,$prostart5,$proend3,$orflen0)= orfParts($orf);
  my $cdnalen= $orf->{cdnalen}||0;
#   my($exons,$orfprot,$prostart5,$proend3,$cdnalen) = @_;
#   # FIXME: need trlen= length(cdnain) for -cdna, and/or use gmap qlen= tag
#   # FIXMEs: exon Split= needs copy to CDS .. losing it, here???
    
  my ($cds5,$cds3,$cds3last,$strand)=(0,0,0,'.');
  my @cds= ();
  my @utr= ();
  ## for phase; need reverse @exons
  my ($cdna1,$inc5,$inc3,$nt_length, $nu5, $nu3)= (0) x 10;
  
  $cdna1= 0; # was 1; # offby1 at end?
  $nt_length= 0; # $prostart5 % 3; #??
  # my $KEEPAN='Split|err|gaps'; #? |Target|trg |gapfill|gapfix ?
  my $DROPAN='Target|trg|gapfill|gapfix|splice';
  
  # ** FIXME 2011Dec : stopcodon split intron >> CDS ends w/o final 1,2 bases ** WRONG
  # .. must make next exon(if exists) part of CDS stop
  
  #?? rev bug here? got bad CDS for good prot/cdna after completeCDSb, revgene
  my $addat="";
  if($exons->[0]->[6] eq "-") {
    my @xbeg= @{$exons->[0]}[3,4,6]; # b,e,o
    my @xend= @{$exons->[-1]}[3,4,6];
    if($xbeg[0] < $xend[0]) {
      $addat.=",badrev"; 
      my @xrev= reverse @$exons; $exons= \@xrev; 
      }
  }
  
  foreach my $exon (@$exons) {
    my ($ref,$src,$xtyp,$xend5, $xend3,$xv,$xo,$xph,$xattr,$gid)= @{$exon};
    
    my $cdsattr=$xattr; 
    $cdsattr=~s/;($DROPAN)=[^;\n]+//g; #** Target|trg has spaces **
    $cdsattr=~s/Parent=[^;\s]+[;]?//; #? or leave on should be same as $gid
    $cdsattr="" unless($cdsattr=~/\w+=/);

    $strand= $xo;
    ($xend5, $xend3)= ($xend3,$xend5) if($xo eq "-"); #patch rev?
    my $xd= abs($xend3 - $xend5); # ?? +1 for width
    
    $cdna1++; # add 1 here, not end loop 
    my $cdna2= $cdna1 + $xd; # ?? +1 for width
    # ** offby1 here ?? YES, need <=, >= to get full CDS stop,start split by intron; see below cdna1=cdna2+1
    #OLD.if($cdna1 < $proend3 and $cdna2 > $prostart5) 
    if($cdna1 <= $proend3 and $cdna2 >= $prostart5) 
    { # overlap
                
      my $d5= ($cdna1 >= $prostart5) ? 0 : $prostart5 - $cdna1; # pos
      my $c5= ($xend5 > $xend3) ? $xend5 - $d5 : $xend5 + $d5;    				  
      
      my $d3= ($cdna2 <= $proend3) ? 0 : $proend3 - $cdna2; # neg
      my $c3= ($xend5 > $xend3) ? $xend3 - $d3 : $xend3 + $d3;
  
      my $elength = 1 + abs($c3 - $c5);
      $nt_length  += $elength;
      $inc3        = $nt_length % 3;
      $inc5        = ($elength - $inc3) % 3; # only care about this one
      # $frame       = ($c5 + $inc5) % 3;
      if ($inc5 == -1) { $inc5 = 2; }
      my $phase= $inc5; # is this right?
      
      my($cb,$ce)= ($c5 > $c3) ? ($c3,$c5): ($c5,$c3); #? rev patch
      # exon Split=, other xattr need copy to CDS HERE ***
      my $rloc= [$ref,$src,"CDS",$cb,$ce,".",$xo,$phase,"Parent=$gid;$cdsattr",$gid]; 
      push @cds, $rloc;
     
      if($cdna1 <= $prostart5) { $cds5=$c5; }
      if($cdna2 >= $proend3) { $cds3=$c3; } else { $cds3last= $c3; } # cdna2 < proend3 here
 
    } elsif(1) { # $addutr .. not used?
      my $d5= ($cdna1 >= $proend3) ? 0 : $proend3 - $cdna1; # pos
      my $u5= ($xend5 > $xend3) ? $xend5 - $d5 : $xend5 + $d5;    				  
      my $d3= ($cdna2 <= $prostart5) ? 0 : $prostart5 - $cdna2; # neg
      my $u3= ($xend5 > $xend3) ? $xend3 - $d3 : $xend3 + $d3;

      my($ub,$ue)= ($u5 > $u3) ? ($u3,$u5): ($u5,$u3); 
      my $up= ($cdna1 < $prostart5) ? "five" : ($cdna2 > $proend3) ? "three" : "odd";
      if($cdna1 < $prostart5) { $nu5++; } elsif($cdna2 > $proend3) { $nu3++; }
      my $rloc= [$ref,$src, $up."_prime_utr",$ub,$ue,".",$xo,0,"Parent=$gid;$cdsattr",$gid]; 
      push @utr, $rloc;
    }
    
    ##$cdna1= $cdna2+1;  # is this off-by-1 now? yes, dont +1 here, do above
    $cdna1= $cdna2;  
  }       
  $cds3||=$cds3last; # partial3
  
  use constant forFASTA => 0;
  my $trlen= ($cdnalen>$cdna1) ? $cdnalen : $cdna1;  
  my($aalen,$pcds,$compl,$orflen)  
     = proteindoc($orf,$trlen,$strand,forFASTA);
  # my($aalen,$pcds,$compl,$orflen,$ocdnalen,$ocdsoffs,$oSelc) 
  #   = proteinqual($orf,$trlen);

  my $mattr="cxlen=$orflen/$trlen;aalen=$aalen,$pcds%,$compl";
  $mattr.= ";protein=$orfprot" if($orfprot);
  $mattr.= ";cdsoff=$prostart5-$proend3"; #? as per ;utroff=$ustart-$uend
  $mattr.= ";cdsspan=$cds5-$cds3$addat"; #?add for other uses? rev if needed? or not?
  $mattr.= ";utrx=$nu5,$nu3" if($nu5 > 2 or $nu3 > 2); # ;utrx=$u5,$u3
  ## mattr keys: cxlen,aalen,protein,utrx
  
  #? resort @cds by loc, not reversed. : let caller do
  # @cds = sort _sortgene @cds;

  ## return also: $clen, $trlen or $utrlen or $ap, $nu5+$nu3, 
  ## ($cdslong, $attrL, $pcodeL, $maxutrL)
  return (\@cds, $mattr, $pcds, _max($nu5,$nu3), \@utr); 
}

# sub getCDSgff2_OLD
# {
#   my($exons,$orf)= @_;
#   
#   # my($exons,$orfprot,$prostart5,$proend3,$addutr) = @_;
#   my ($orfprot,$prostart5,$proend3)= orfParts($orf);
#   # Ooops, orf->{start,stop} are reversed for strand=rev; start>stop
#     
#   my ($cds5,$cds3)=(0,0);
#   my @cds= ();
#   my @utr= ();
#   ## for phase; need reverse @exons
#   my ($cdna1,$inc5,$inc3,$nt_length, $nu5, $nu3, $strand)= (0) x 10;
#   $cdna1= 1;
#   $nt_length= 0; # $prostart5 % 3; #??
#   
#   # ** FIXME 2011Dec : stopcodon split intron >> CDS ends w/o final 1,2 bases ** WRONG
#   # .. must make next exon(if exists) part of CDS stop
#   
#   foreach my $exon (@$exons) {
#     my ($ref,$src,$xtyp,$xend5, $xend3,$xv,$xo,$xph,$xattr,$gid)= @{$exon};
#     $strand= $xo;
#     ($xend5, $xend3)= ($xend3,$xend5) if($xo eq "-"); #patch rev?
#     
#     my $xd= abs($xend3 - $xend5); # ?? +1 for width
#     
#     my $cdna2= $cdna1 + $xd; # ?? +1 for width
#     # ** offby1 here ?? YES, need <=, >= to get full CDS stop,start split by intron; see below cdna1=cdna2+1
#     #OLD.if($cdna1 < $proend3 and $cdna2 > $prostart5) 
#     if($cdna1 <= $proend3 and $cdna2 >= $prostart5) 
#     { # overlap
#                 
#       my $d5= ($cdna1 >= $prostart5) ? 0 : $prostart5 - $cdna1; # pos
#       my $c5= ($xend5 > $xend3) ? $xend5 - $d5 : $xend5 + $d5;    				  
#       
#       my $d3= ($cdna2 <= $proend3) ? 0 : $proend3 - $cdna2; # neg
#       my $c3= ($xend5 > $xend3) ? $xend3 - $d3 : $xend3 + $d3;
#   
#       my $elength = 1 + abs($c3 - $c5);
#       $nt_length  += $elength;
#       $inc3        = $nt_length % 3;
#       $inc5        = ($elength - $inc3) % 3; # only care about this one
#       # $frame       = ($c5 + $inc5) % 3;
#       if ($inc5 == -1) { $inc5 = 2; }
#       
#       my $phase= $inc5; # is this right?
#       
#       my($cb,$ce)= ($c5 > $c3) ? ($c3,$c5): ($c5,$c3); #? rev patch
#       my $rloc= [$ref,$src,"CDS",$cb,$ce,".",$xo,$phase,"Parent=$gid",$gid]; 
#       push @cds, $rloc;
#      
#       $cds5=$c5 if($cdna1 <= $prostart5);
#       $cds3=$c3 if($cdna2 >= $proend3);
#       
#     } elsif(1) { 
#       my $d5= ($cdna1 >= $proend3) ? 0 : $proend3 - $cdna1; # pos
#       my $u5= ($xend5 > $xend3) ? $xend5 - $d5 : $xend5 + $d5;    				  
#       my $d3= ($cdna2 <= $prostart5) ? 0 : $prostart5 - $cdna2; # neg
#       my $u3= ($xend5 > $xend3) ? $xend3 - $d3 : $xend3 + $d3;
# 
#       my($ub,$ue)= ($u5 > $u3) ? ($u3,$u5): ($u5,$u3); 
#       my $up= ($cdna1 < $prostart5) ? "five" : ($cdna2 > $proend3) ? "three" : "odd";
#       if($cdna1 < $prostart5) { $nu5++; } elsif($cdna2 > $proend3) { $nu3++; }
#       my $rloc= [$ref,$src, $up."_prime_utr",$ub,$ue,".",$xo,0,"Parent=$gid",$gid]; 
#       push @utr, $rloc;
#     }
#     
#     $cdna1= $cdna2+1;
#   }       
#   
#   
#   my $trlen= $cdna1; # = $nt_length
#   use constant forFASTA => 0;
#   my($aalen,$pcds,$compl,$orflen) ## ,$fahead,$orfprot2
#      = proteindoc($orf,$trlen,$strand,forFASTA);
#   
# #  #....
# #   my $aalen=length($orfprot); 
# #   $aalen-- if(substr($orfprot,-1) eq '*');
# #   my $clen= $aalen * 3; # can be off by -1,-2 here. 
# #   my $ap=int(100 * $clen/$trlen);
# #  #....
# 
#     
#   my $mattr="cxlen=$orflen/$trlen;aalen=$aalen,$pcds%,$compl";
#   $mattr.= ";protein=$orfprot" if($orfprot);
# #   if($orfprot) {
# #     my $p5= (substr($orfprot,0,1) eq 'M')?0:1;
# #     my $p3= (substr($orfprot,-1,1) eq '*')?0:1;
# #     my $prostat = ($p5 and $p3) ? "partial": ($p5)? "partial5" : ($p3)? "partial3" :"complete";
# #     $mattr.= ",$prostat;protein=$orfprot";
# #   }
#   
#   $mattr.= ";utrx=$nu5,$nu3" if($nu5 > 2 or $nu3 > 2); # ;utrx=$u5,$u3
#   ## mattr keys: cxlen,aalen,protein,utrx
#   
#   # FIXME: resort @cds by loc, not reversed. : let caller do
#   # @cds = sort _sortgene @cds;
# 
#   ## return also: $clen, $trlen or $utrlen or $ap, $nu5+$nu3, 
#   ## ($cdslong, $attrL, $pcodeL, $maxutrL)
#   ## add: $cds5,$cds3 ??
#   return (\@cds, $mattr, $pcds, _max($nu5,$nu3), \@utr); 
# }

  
sub getAllOrfs { # added strands to return both, or rev
  my ($input_sequence, $strands,$flags) = @_;
  $strands ||= "fwd";
  $flags ||="";
  # add fwd,rev orfs here ? gff-caller wants only 1 strand at a time, need option 
  
  return undef unless ($input_sequence or length ($input_sequence) >= 3) ;
  $input_sequence = uc ($input_sequence); # was lc() change all to uc()
  
#  ## fixme: this screws seq indices: start,stop; instead of chomp, restrict @starts,@stops
#   if($flags =~ /dropnnn|chomp/i) {
#     $input_sequence =~ s/^N+//;
#     $input_sequence =~ s/N+$//; 
#   }

  my $seqlen=length($input_sequence);
  
  # 201402 fix: orthobest input cds span for bestorf..
  # .. innerstop problem: retained introns looks like cause,
  # .. can add ref/tr hsp list, which gives some info on introns in un-aligned spans..
  # .. then need messy orf-calling using tr-hsp list with extensions up to find intron
  # .. or some other intron splice site seek to cut it out of trasm for mRNA/CDS/aa
  # .. try scan for splice sites fwd:AG/GT (rev:AC/CT )  ?
  
  my $CheckCutIntrons= 1; # debug
  
  my @bestspan=(); my @besthsps=(); # add @besthsps,@refhsps for innerstop,introncuts ?
  my $bestspanorf=undef; my $checkspanframe=0; my $hasHSPs=0;
  my $bestspanflag="";
  if($flags =~ /bestspan=([\w,\:\+\-]+)/) { # \d to \w to get all parts?
    $bestspanflag=$1; # my $bspan=$1;  ## now "bestspan=$tralnspan:$cla:$rid:$refaaspan";
    $checkspanframe= ($flags=~/frameck/)?1:0;
    # FIXME: spans here are 1origin, beststart/getOneOrf want 0origin ??
    my($cbe,$co)=split":", $bestspanflag; ## 2341-1364:-
    ## allow hsps here?  22-333,444-555:+ ?
    if($cbe=~/,\d/) {
      my($mcb,$mce,$lb,$le)=(-1,0,0,0); 
      for my $ch (split",",$cbe) { 
        my($cb,$ce)= split '-',$ch; next unless($ce>0);
        if($co eq '-') { ($cb,$ce)=($seqlen-$ce,$seqlen-$cb); } 
        else { $cb--; $ce--; } # 0origin fix, BUT not rev !!
        if($le>0 and $co eq '-' and $ce>=$lb) { $ce=$le; $besthsps[-2]=$cb if($cb<$lb); }
        elsif($le>0 and $cb<=$le) { $cb=$lb; $besthsps[-1]=$ce if($ce>$le); }
        else { push @besthsps, $cb,$ce; }
        ($lb,$le)=($cb,$ce);
        $mcb=$cb if($mcb==-1 or $cb<$mcb); $mce=$ce if($mce<$ce);
      }
      $hasHSPs=(@besthsps>2)? scalar(@besthsps) : 0;
      ## push @besthsps,$co;       
      @bestspan=($mcb,$mce,$co);
    } else {
    my($cb,$ce)= split /[\-]/, $cbe; 
    if($cb>$ce) { $co='-'; ($cb,$ce)=($ce,$cb); } # doublerev, or (cb,ce)=(seqlen-cb,seqlen-ce);
    if($co eq '-') { ($cb,$ce)=($seqlen-$ce,$seqlen-$cb); } 
    else { $cb--; $ce--; } # 0origin fix !!! BUT not rev !!
    @bestspan=($cb,$ce,$co); # but need to extend to @start,@stop !
    }
  }
  
  my (@starts, @stops, @orfs);
  my $workseq= $input_sequence;
  unless($strands =~ /^(r|\-)/) {  # forward_strand_only();
    @stops  = identify_putative_stops($workseq);
    @starts = identify_putative_starts($workseq,\@stops);
    push @orfs, get_orfs (\@starts, \@stops, $workseq, '+');
    if(@bestspan and $bestspan[2] eq '+') { 
      my ($bb,$be,$bo)=@bestspan; 
      
      # BUG: got shorter aa from bestspan, but same start .. found innerstop here, but not w/o span
      # .. seems to be problem of extended aln, from possibly alt-exon in same trasm: t8453-8854/r2703-2836 adds stops
      # def  =socatfishv1k39loc150483t1 aalen=2617,74%,complete; clen=10510; strand=+; offs=213-8066; 
      # bspan=socatfishv1k39loc150483t1 aalen=2593,74%,complete; clen=10510; strand=+; offs=213-7994;
      #   orflags=bestspan:213-4487,4905-6608,8453-8854:+:clmiss:ref:Funhe2EKm027765t1:1-1314,1399-1972,2703-2836,incut:20-stops:1751-bp; 
      ## aln: Funhe2EKm027765t1  socatfishv1k39loc150483t1  213-8854:+      over:213-8066:+
      
      # FIXME: need to find best frame, bb can be offby 1,2.. maybe not, blastn may be on-frame
      ## fix for stopatend: (substr($orfprot,-1) eq '*')  
      my $aa0= getOneOrf($bb,$be,$workseq,$bo); 
      my $nstop= $aa0->{innerstop}; 
      my $aa0len= $aa0->{goodlen};
      
      # # .. try scan for splice sites fwd:AG/GT (rev:AC/CT )  ?
      if($nstop>0 and $CheckCutIntrons and $hasHSPs ) {  
        my @icut=();
        for(my $i=1; $i<$hasHSPs-1; $i+=2) {
          my($b1,$e1,$b2,$e2)= @besthsps[$i-1,$i,$i+1,$i+2];
          my $ospan= getOneOrf($b1,$e2,$workseq,$bo);  # check if bad segment!
          my $instop= $ospan->{innerstop};
          next unless($instop);

          my $es= index($workseq,'AG',$e1); 
          my $bs= rindex($workseq,'GT',$b2);
          if($es>$e1 and $bs>$es) { 
            push @icut,$es,$bs+2;  
          } else {
            push @icut,$e1+1,$b2;  # yes, this doubles n unstopped cases
          }
        }
        if(@icut) {
          my $bat=0; my $mrna=""; my $ncut= $#icut; my $spancut=0;
          my $cb=$besthsps[0]; my $ce=$besthsps[-1]; # needs adjust..
          for(my $i=0; $i<$ncut; $i+=2) { 
            my($sb,$se)= @icut[$i,$i+1]; 
            $mrna .= substr($workseq,$bat,$sb-$bat);
            $bat = $se; my $cut=$se-$sb; $ce -= $cut; $spancut+=$cut;
            }
          $mrna .= substr($workseq,$bat);  
          my $aac= getOneOrf($cb,$ce,$mrna,$bo); 
          my $nsc= $aac->{innerstop};  
          my $aa1len= $aac->{goodlen};  
          # DEBUG out here..
          warn "#DBG CutIntrons: aa0w=$aa0len, aa1w=$aa1len, stop0=$nstop, stop1=$nsc, hsps=@besthsps, icut=@icut, flags=$flags\n"
            if($DEBUG);
          ## this works for ~50/500 cases in catfish1evg8.. enough to use? any other fixes?
          ## using full cut e1..b2 ups that to 120/500 fixed no innerstop, usable
          ## .. possibly test for stops in each mrnasegment add, scan over istop? ie shift se=>bat above?
          ## .. some of crap is not full intron insert, but partial .. test cut e1..b2 w/o AG/GT scan?
          
          if($nsc < 1) {
            $bestspanflag.=",incut:$nstop-stops:$spancut-bp";
            ($bb,$be)= ($cb,$ce); $nstop=$nsc; 
            $workseq= $mrna; $seqlen=length($mrna);
            @stops  = identify_putative_stops($workseq);
            @starts = identify_putative_starts($workseq,\@stops);
          }
        }
      }   
      
      if($nstop>0 and $checkspanframe) {
        my $aa1= getOneOrf($bb+1,$be,$workseq,$bo);  my $ns1= $aa1->{innerstop};  
        my $aa2= getOneOrf($bb+2,$be,$workseq,$bo);  my $ns2= $aa2->{innerstop};  
        my $f=0; 
        if($ns1 < $nstop) { $f=1; $nstop=$ns1; } 
        if($ns2 < $nstop) { $f=2; }
        $bb += $f;
      }
      
      $bb=beststart($bb,0,\@starts,\@stops);
      $be=beststop($be,$bb,$seqlen,\@stops);
      $bestspanorf = getOneOrf($bb,$be,$workseq,$bo);   @bestspan= ($bb,$be,$bo);
      $bestspanorf->{flags}= "bestspan:$bestspanflag";
      # unshift @orfs,$bestspanorf;  # check all orfs for this tho, no dupl.
      }
    }
    
  if($strands =~ /^(b|r|\-)/) { # rev|r|both
    my $revseq= revcomp($input_sequence); 
    $workseq= $revseq;
    @stops  = identify_putative_stops($workseq);
    @starts = identify_putative_starts($workseq,\@stops);
    push @orfs, get_orfs (\@starts, \@stops, $workseq, '-');
    if(@bestspan and $bestspan[2] eq '-') { 
      my ($bb,$be,$bo)=@bestspan; 
      # FIXME: need to find best frame, bb can be offby 1,2..
      my $aa0= getOneOrf($bb,$be,$workseq,$bo);
      my $nstop= $aa0->{innerstop}; # {protein} =~ tr/*/*/;
      if($nstop>0 and $checkspanframe) {
        my $aa1= getOneOrf($bb+1,$be,$workseq,$bo);  my $ns1= $aa1->{innerstop}; #{protein} =~ tr/*/*/;
        my $aa2= getOneOrf($bb+2,$be,$workseq,$bo);  my $ns2= $aa2->{innerstop}; #{protein} =~ tr/*/*/;
        my $f=0; if($ns1 < $nstop) { $f=1; $nstop=$ns1; } 
        if($ns2 < $nstop) { $f=2; }
        $bb += $f;
      }

      $bb= beststart($bb,0,\@starts,\@stops);
      $be= beststop($be,$bb,$seqlen,\@stops);
      $bestspanorf = getOneOrf($bb,$be,$workseq,$bo);  @bestspan= ($bb,$be,$bo);
      $bestspanorf->{flags}= "bestspan:$bestspanflag";
      # unshift @orfs,$bestspanorf; # unless grep{$_ eq $bestorf}@orfs; # check all orfs for this tho, no dupl.
      }
  }

  
  if (@orfs or $bestspanorf) {  # get both longest complete, longest partial
    ## set in order of decreasing length
    my $keylen=($USEGOODLEN)?"goodlen":"length";
    @orfs = sort {$b->{$keylen} <=> $a->{$keylen} or $b->{complete} <=> $a->{complete}} @orfs;
    
    # FIXME: for bestspan.. but cancel if innerstop > 0 or > 1 .. esp. scrambled bestspan's w/ ## stops
    # FIXME2: cancel if bestspan goodlen << longest goodlen, happening w/ odd hsps, maybe alt exons; need ortho vals
    if($bestspanorf) {
      my $nstop= $bestspanorf->{innerstop}; 
      if($nstop > 0 and $InnerStopToX<0) { 
         # skip but add info to orf0
        my $ncds= $bestspanorf->{length};
        if(@orfs) { $orfs[0]->{flags}= "cancelspan:innerstop.$nstop,cdsw.$ncds,$bestspanflag"; } #??
      } else { 
        if($nstop > 0 and $InnerStopToX==1) { substr( $bestspanorf->{protein}, 0,-1) =~ tr/*/X/; }
        unshift @orfs,$bestspanorf; 
        # $longest = $bestspanorf;
      }
    }
    
    # FIXME: option to mark/return 2ndary orf(s) in aberrant long-utr transcripts, that dont overlap 1st orf
    # FIXME 1505, partial 1aa > partial3 bug from utrorf.mrna cut just past partial3 cds.
    # ?? look for 2nd longest, complete > longest.complete ?

    my $longest = $orfs[0];   
    if($longest->{innerstop}) { # 201504 fix: drop innerstops unless desired
      unless($InnerStopToX > 0) { ($longest) = grep { $_->{innerstop} == 0 } @orfs; }
      #?? if($InnerStopToX <= 0) { ($longest) = grep { $_->{innerstop} == 0 } @orfs; }
      # else { } # convert to X
    }
    my($longfull) = grep { $_->{complete} == 3 && $_->{innerstop} == 0 } @orfs; # add  $_->{innerstop} == 0
    
    ## upd1712: cancel this for $flag =~ /long/ 
    if($longest->{complete} == 0 and not ($flags =~ m/long/)) { # utrorf bug 1505
      my($loc)= grep { $_->{complete} > 0 &&  $_->{innerstop} == 0 } @orfs; 
      if(not $loc) { } elsif($longfull and $longfull eq $loc) { }
      elsif($loc->{goodlen} > $ORF_FULLvPART * $longest->{goodlen}) { $longest=$loc; }
    }

    return ($longest, $longfull, \@orfs);
  } else {
    return undef;
  }
}

=item newOrf

  my $orf= newOrf( sequence => 'acgtcgat', protein => "MABCD", complete =>  3, start => 111, stop => 999, orient => '+', length => 1999 );

  see get_orfs() makes these from data
        my $orf = { sequence => $orfSeq, protein => $protein,
          start=>$start, stop=>$stop, orient=>$direction,
          length=>$orflen, goodlen=>$goodlen, cdnalen=>$seq_length,
          complete=>$isfull, innerstop=>$innerstop,
          };

   name? orfOf(xxx) or  newOrf(xxx)
   
=cut 

my @ORFPARTS= qw(sequence protein start stop orient length goodlen cdnalen complete innerstop);
# add Selc ? other sometimes keys?

sub newOrf { 
  my(%parts)= @_; 
  my %orf=(); unless(%parts){ %parts= (); }
  my $ORP=join '|', @ORFPARTS;
  for my $k (@ORFPARTS){ $orf{$k}= $parts{$k}||0; }
  for my $k (grep{ not m/$ORP/ } keys %parts) { $orf{$k}= $parts{$k}; } #? add all parts not std?

  # check data: start,stop required, stop < start  when -orient 
  # minimal: newOrf( protein=>"xxxx", start=>99,  stop=>999, cdnalen=>1999);
  #  NoStopCodon problems, *should* only chop for output
  #   if(substr($orfprot,-1,1) eq '*') { $aalen--; if($NoStopCodon) { $orfprot =~ s/\*$//; } }

  my $len= ($orf{stop} < $orf{start})? 1 + $orf{start} - $orf{stop} : 1 + $orf{stop} - $orf{start}; 
  $orf{length}= $len if($orf{length} < $len);
  $len= $orf{length};
  $orf{goodlen} ||= $len;
  #? $orf{cdnalen}= $len if($orf{cdnalen} < $len); #? or 0 for missing?
  $orf{orient} ||= '.'; $orf{protein} ||= '';  $orf{sequence} ||= '';
  
  if(my $protein= $orf{protein}) {
    my $isfull= $orf{complete}; # either comp|2 or protein* here
    if(substr($protein,0,1) eq 'M') { $isfull |= 1; }
    if(substr($protein,-1,1) eq '*') { $isfull |= 2; } elsif($NoStopCodon) { } #? need to trust {complete}
    my $innerstop= $protein =~ tr/*/*/;   
    $innerstop-- if($innerstop>0 and ($isfull & 2));  
    $orf{complete}= $isfull if($isfull > $orf{complete});
    $orf{innerstop}= $innerstop if($innerstop>0 and $orf{innerstop}<1);
    }
  
  return(\%orf);
}

sub orfParts
{
  my($orf,$parts) = @_;
  # $parts= [ qw(protein start stop length orient) ] unless(ref $parts); 
  if(ref $orf) {
    if(ref $parts) { return @$orf{ @$parts }; }
    return($orf->{protein}, $orf->{start}, $orf->{stop}, $orf->{length}, $orf->{orient});  
  } else {
    return("",0,0,0,0);
  }
}


sub beststart {
  my( $instart,$minpos,$starts_ref,$stops_ref) = @_;
  my $inframe = $instart % 3; 
  my $spat=$instart;
  ## fixme need to stop at stops before starts
  my @startrev= reverse grep { $_ <=  $instart } @{$starts_ref};
  my @stoprev = reverse grep { $_ <=  $instart } @{$stops_ref};
  my $stopat=0;
  foreach my $ss (@stoprev) { if($ss % 3 == $inframe) { $stopat= $ss; last; } }
  foreach my $sp (@startrev) { if($sp<$stopat) { last; } elsif($sp % 3 == $inframe) { $spat= $sp; last; } }
  # $spat=$stopat+3 if($spat < $stopat);
  return $spat;
#  
#  #......
#   foreach my $sp (@{$starts_ref}) {
#     last if($sp > $instart); next if($sp<$minpos);
#     $spat= $sp if($sp % 3 == $inframe);
#   }
#  return $spat;
}

sub beststop {
  my( $inend,$instart,$maxpos,$stops_ref) = @_;
  my $inframe = $instart % 3; 
  my $spat=$inend;
  foreach my $sp (@{$stops_ref}) {
    last if($sp > $maxpos); next if($sp < $inend - 3);
    if($sp % 3 == $inframe) { $spat= $sp; last; }
  }
  return $spat;
}

sub getOneOrf {
  my( $start_pos, $stop_pos, $mrnaseq, $direction, $flags) = @_;
  
  #* Need option to extend start,stop; ** Need to look for start codon, stop codon if start/stop not on these.
  #* Assume seq == revcomp(seq) if dir eq '-' ??
  
  my $seq_length = length ($mrnaseq);
  my $sp_frame = $start_pos % 3;
  
  ## FIXME: ?? bad for dir=-
#  my $isfull5= ($start_pos+3 <= $seq_length and substr($mrnaseq, $start_pos, 3) eq "ATG")?1:0;
#  
#  #?? or use: foreach my $start_pos (@{$starts_ref})
#  ## use also: foreach my $stop_pos (@{$stops_ref}) 
#   if(!$isfull5 and $start_pos > 2 and $flags =~ /full|seekstart/i ) { #??
#     for (my $b=$start_pos - 3; $b > 2 && !$isfull5; $b -= 3) {
#       if(substr($mrnaseq, $b, 3) eq "ATG") { $isfull5=1; $start_pos=$b; last; }
#       }
#   }
#  
#   my $start_pos_frame = $start_pos % 3;
#   if($isfull5) { $last_stop_full{$start_pos_frame} = $stop_pos; }
#   else { $last_stop_part{$start_pos_frame} = $stop_pos; }

  my $stopplus3= _min( $stop_pos+3, $seq_length); #dgg patch; getting overruns.
  my $orflen = ($stopplus3 - $start_pos);
  my $stopoff = $orflen % 3;
  if($stopoff>0) {  $orflen -= $stopoff; $stopplus3 -= $stopoff; } # do before pull seq..
  my ($start_pos_adj, $stop_pos_adj) = ( ($start_pos+1), $stopplus3);
  my ($start, $stop) = ($direction eq '-') 
    ? (revcomp_coord($start_pos_adj, $seq_length), revcomp_coord($stop_pos_adj, $seq_length))
    : ($start_pos_adj, $stop_pos_adj);
  
  my $orfSeq =  substr ($mrnaseq, $start_pos, $orflen); #include the stop codon too.
  my $protein= translate_sequence($orfSeq, 1); ## from Nuc_translator
  $orflen= length($orfSeq);

  my $isfull= 0;
  $isfull |= 1 if(substr($protein,0,1) eq 'M');  
  $isfull |= 2 if(substr($protein,-1,1) eq '*');  
  my $innerstop= $protein =~ tr/*/*/;   
  $innerstop-- if($innerstop>0 and ($isfull & 2));  
  my $nxxx= $innerstop + $protein =~ tr/X/X/;
  my $goodlen= $orflen - 3*$nxxx; # goodlen to avoid XXXXXXX* crap
  
  if($DEBUG>1){  
    my $aalen = length($protein); my $u=index($protein,'u');
    print STDERR "#DBG: orf1,olen=$orflen,alen=$aalen,at=$start-$stop,selc=$u,aa1=",
      substr($protein,0,9),"..",substr($protein,-3,3),"\n";
  }
  
  # orf maybe add: mRNAseq => $mrnaseq, as output help, e.g. for rev and intron-cut seqs
  my $orf = { sequence => $orfSeq, protein => $protein, 
    mrna => $mrnaseq, #? only for getOneOrf 
    start=>$start, stop=>$stop, orient=>$direction,
    length=>$orflen, goodlen=>$goodlen, cdnalen=>$seq_length,
    complete=>$isfull, innerstop=>$innerstop,
    };

  # USESelenocysteine add Selc position in orf for exception reporting, NCBI /transl_except=(pos:1002..1004,aa:Sec); in mRNA span
  if($USESelenocysteine and (my $iu=index($protein,'u')) > 0) {
    my $Selc=""; ## do any have >1 u ?? yes
    while($iu >= 0) { 
      my $uoff=$start_pos_adj + 3*$iu;
      $uoff= revcomp_coord($uoff, $seq_length) if($direction eq '-');
      $Selc.= "$uoff,"; $iu= index($protein,'u',$iu+1);  
      }
    $Selc=~s/,$//; $orf->{Selc}=$Selc;
  }
    
  return $orf;
  # return (wantarray) ? ($orf, $orf, [$orf]) : $orf;
}




sub get_orfs {
  my( $starts_ref, $stops_ref, $mrnaseq, $direction) = @_;
  
  # $direction here only affects start,stop relative to forward cdna seq .. use it.
  # .. input seq is already rev(seq) for -dir
  # Ooops, orf->{start,stop} are reversed for -strand; start>stop : ($lb,$le)=($le,$lb) if($lb>$le);
  # .. fix some places, leave that way others?
    
  # unless ($starts_ref && $stops_ref && $mrnaseq && $direction) { warn "get_orfs: params not appropriate"; }  
  # want only max orf, complete + partial
  # FIXME: option to mark/return 2ndary orf(s) in aberrant long-utr transcripts, that dont overlap 1st orf
   
	my %last_stop_part = ( 0=>-1, 1=>-1,  2=>-1); # position of last chosen stop codon in spec reading frame.
	my %last_stop_full = ( 0=>-1, 1=>-1,  2=>-1); # position of last chosen stop codon in spec reading frame.
  my @orfs;
  my $seq_length = length ($mrnaseq);
  my $norf=0;
  
  ##  $mrnaseq .= "###" if($DEBUG); # why're we getting bad prot at partial end? isnt prot but cds range needs adjust for -2,-1
  
  foreach my $start_pos (@{$starts_ref}) {
		my $start_pos_frame = $start_pos % 3;
		# includes partial5: starts at 0,1,2 unless real start
		# ** PROBLEM? -- for internal mis-assemblies, should
		#   test partial-starts also following each stop? dont require ATG start inside?
		
		my $isfull5= ($start_pos+3 <= $seq_length and substr($mrnaseq, $start_pos, 3) eq "ATG")?1:0;
		
		foreach my $stop_pos (@{$stops_ref}) {
		  if ( ($stop_pos > $start_pos) && #end3 > end5
				 ( ($stop_pos - $start_pos) % 3 == 0) && #must be in-frame
				   # ($start_pos > $last_delete_pos{$start_pos_frame}) # dgg: bad for partial5 + full
				 ($isfull5 ? $start_pos > $last_stop_full{$start_pos_frame} : $start_pos > $last_stop_part{$start_pos_frame} )
				 ) #only count each stop once.
			{
		    # includes partial3: stops at end- 0,1,2 unless real stop
				
				#dgg.no# $last_delete_pos{$start_pos_frame} = $stop_pos;
			  my $oldlast_stop_full= $last_stop_full{$start_pos_frame}; # TRIMCDS fixup
				if($isfull5) { $last_stop_full{$start_pos_frame} = $stop_pos; }
			  else { $last_stop_part{$start_pos_frame} = $stop_pos; }
			
			  my $stopplus3= _min( $stop_pos+3, $seq_length); #dgg patch; getting overruns.
				my $orflen = ($stopplus3 - $start_pos);
				my $stopoff = $orflen % 3;
				if($stopoff>0) {  $orflen -= $stopoff; $stopplus3 -= $stopoff; } # do before pull seq..
			  
				my ($start_pos_adj, $stop_pos_adj) = ( ($start_pos+1), $stopplus3);
				# my ($start, $stop) = ($direction eq '+') ? ($start_pos_adj, $stop_pos_adj) 
				#	: (&revcomp_coord($start_pos_adj, $seq_length), &revcomp_coord($stop_pos_adj, $seq_length));
        my ($start, $stop) = ($direction eq '-') 
          ? (revcomp_coord($start_pos_adj, $seq_length), revcomp_coord($stop_pos_adj, $seq_length))
          : ($start_pos_adj, $stop_pos_adj);

        ## reuse sub getOneOrf() above here ?				
	      my $orfSeq =  substr ($mrnaseq, $start_pos, $orflen); #include the stop codon too.
				
        # FIXME here? gap trim NNN/XXX at start/stop of coding seq, adjusting start/stop
        # .. this is causing bug for MxxxxxABCD cases.. > innerstops despite block; phase problem; make i%3 == 0 
	      # * should do this before/in getStarts() otherwise trim can miss M a few aa after gap.

			  if($TRIMCDSENDGAP) {
			    my $at= index($orfSeq,'NNN');
			    if($at >= 0 and $at < 4) {
			      my $i=$at+2; while($i<$orflen and substr($orfSeq,$i,1) eq 'N') { $i++; }
 			      $i++ while($i % 3 > 0); # 1505 patch solves innerstop bug
			      $orfSeq= substr($orfSeq,$i); 
			      $orflen= length($orfSeq);
			      $start += $i; $start_pos+= $i;
			## more fixups:
			      my $newfull5= (substr($orfSeq, 0, 3) eq "ATG")?1:0;
                              if($isfull5 and not $newfull5) { 
			         $last_stop_full{$start_pos_frame}= $oldlast_stop_full; ##not delete $last_stop_full{$start_pos_frame};  
                                 $last_stop_part{$start_pos_frame} = $stop_pos; $isfull5=$newfull5;
				}

			    }
			    # $at= rindex($orfSeq,'NNN'); # problems w/ stop codon change ..
			    # if($at >= $orflen - 3) { }
			  }
			  	
				my $protein= translate_sequence($orfSeq, 1); ## from Nuc_translator
			  $orflen= length($orfSeq);
  
				my $isfull= 0;
        $isfull |= 1 if(substr($protein,0,1) eq 'M'); # $protein =~ /^M/
        $isfull |= 2 if(substr($protein,-1,1) eq '*'); # $protein =~ /\*$/
 	      my $innerstop= $protein =~ tr/*/*/;   
        $innerstop-- if($innerstop>0 and ($isfull & 2));  
        ## flag problem if XXX count > non-X count
        my $nxxx= $innerstop + $protein =~ tr/X/X/;
        my $goodlen= $orflen - 3*$nxxx;
        $norf++;
        # opt to drop *stopcodon, here or where, 
        
        if($DEBUG>1){ 
  				my $aalen = length($protein); my $u=index($protein,'u');
          print STDERR "#\n" if($norf==1);  
          #no,wrong# my $dc=  ($direction eq '+') ? 'fwd' : 'rev';
          print STDERR "#DBG: orf$norf,olen=$orflen,alen=$aalen,at=$start-$stop,selc=$u,aa1=",
            substr($protein,0,9),"..",substr($protein,-3,3),"\n";
        }
        
        # my $orf= newOrf( sequence => $orfSeq, protein => $protein,
        #   start=>$start, stop=>$stop, orient=>$direction,
        #   length=>$orflen, goodlen=>$goodlen, cdnalen=>$seq_length,
        #   complete=>$isfull, innerstop=>$innerstop, );
        my $orf = { sequence => $orfSeq, protein => $protein,
          start=>$start, stop=>$stop, orient=>$direction,
          length=>$orflen, goodlen=>$goodlen, cdnalen=>$seq_length,
          complete=>$isfull, innerstop=>$innerstop,
          };
          
        # add Selc position in orf for exception reporting, NCBI /transl_except=(pos:1002..1004,aa:Sec); in mRNA span
        if($USESelenocysteine and (my $iu=index($protein,'u')) > 0) {
          my $Selc=""; ## do any have >1 u ?? yes
          while($iu >= 0) { 
            my $uoff=$start_pos_adj + 3*$iu;
            $uoff= revcomp_coord($uoff, $seq_length) if($direction eq '-');
            $Selc.= "$uoff,"; $iu= index($protein,'u',$iu+1);  
            }
          $Selc=~s/,$//; $orf->{Selc}=$Selc;
        }
          
				push (@orfs, $orf);
				last;   
			}
		}
  }
  return (@orfs);
}



sub identify_putative_starts {
  my ( $seq, $stops_aref) = @_;
  my %starts;
  my %stops;
  foreach my $stop (@$stops_aref) {
		$stops{$stop} = 1;
    }
	
  if(1) {  # ($self->{ALLOW_5PRIME_PARTIALS} || $self->{ALLOW_NON_MET_STARTS}) 
    my $i=0;
    if($seq =~ /^N/) {  # FIXME: skip leading NNN
      my $n= length($seq);
      $i++ while( $i<$n && substr($seq,$i,1) eq 'N');
    }
		$starts{$i} = 1 unless $stops{$i};
		++$i; $starts{$i} = 1 unless $stops{$i};
		++$i; $starts{$i} = 1 unless $stops{$i};
    }
    
  if (1) {  # ! $self->{ALLOW_NON_MET_STARTS}  #Look for ATG start codons.
		my $start_pos = index ($seq, "ATG");
		# cases of MxxxxxABCD gaps are problem for a few.. cancel if ATGnnnnnnn ?
		while ($start_pos != -1) {
			$starts{$start_pos} = 1;
			$start_pos = index ($seq, "ATG", ($start_pos + 1));
		  }
    } 

  my @starts = sort {$a<=>$b} keys %starts;
  return (@starts);
}


sub identify_putative_stops {
  my ($seq) = @_;
  my %stops;
  
  if(1){  # $self->{ALLOW_3PRIME_PARTIALS}
	## count terminal 3 nts as possible ORF terminators.
	my $e = length ($seq);
  if($seq =~ /N$/) { # FIXME: skip trailing NNN
    $e-- while( $e>1 && substr($seq,$e-1,1) eq 'N');
  }
	$stops{$e} = 1;
	$e--; $stops{$e} = 1;
	$e--; $stops{$e} = 1;
  }
  
  # my @stop_codons = @{$self->{stop_codons}};
  # my @stop_codons = &Nuc_translator::get_stop_codons(); # live call, depends on current genetic code.
  #global# my @stop_codons = qw(TAA TAG TGA);
## add frameshift detect option here?  only if remaining seq >> nnn, only if %cds/utr falls below pCDSpoor ?
## need indel max to test: offby -2,-1,1,2 only to shift $i
  
  foreach my $stop_codon (@stop_codons) {
    my $stop_pos = index ($seq, $stop_codon);
    while ($stop_pos != -1) {
      $stops{$stop_pos} = 1;
      $stop_pos = index ($seq, $stop_codon, ($stop_pos + 1)); #include the stop codon too.
      }
  }
  my @stops = sort {$a<=>$b} keys %stops;
  return (@stops);
}


sub revcomp {
    my ($seq) = @_;
    my $reversed_seq = reverse ($seq);
    $reversed_seq =~ tr/ACGTacgtyrkm/TGCAtgcarymk/;
    return ($reversed_seq);
}


sub revcomp_coord {
    my ($coord, $seq_length) = @_;
    return ($seq_length - $coord + 1);
}

sub backtranslate_protein {
  my ($sequence, $useAmbiguousNucs) = @_;
  $sequence = uc ($sequence); # redundant ??
  my $seq_length = length ($sequence);
  my $cds_sequence="";
  if(defined $useAmbiguousNucs) {
    if( $useAmbiguousNucs ) { %amino2codon_table= %amino2codon_ambig ;  }
    else { %amino2codon_table= %amino2codon_fixed; } 
  }
  
  for (my $i = 0; $i < $seq_length; $i++) {
      my $codon;
      my $amino = substr($sequence, $i, 1);  # deal with non-alpha
      if (exists($amino2codon_table{$amino})) {
        $codon = $amino2codon_table{$amino};
      } elsif($amino =~ /[A-Z]/) {
        $codon = 'NNN'; # fixme
      } else {
        $codon= ''; # eat it?
      }
      $cds_sequence .= $codon;
  }
  return($cds_sequence);
}
      

sub translate_sequence {
  my ($sequence, $frame) = @_;
    
  $sequence = uc ($sequence); # redundant now
	$sequence =~ tr/U/T/;
  my $seq_length = length ($sequence);
  unless ($frame >= 1 and $frame <= 6) { 
		warn "Frame $frame is not allowed. Only between 1 and 6"; # die?
		return -1; # 
	}
	
	if ($frame > 3) {
		# on reverse strand. Revcomp the sequence and reset the frame
		$sequence = revcomp($sequence);
		if ($frame == 4) {
			$frame = 1;
		}
		elsif ($frame == 5) {
			$frame = 2;
		}
		elsif ($frame == 6) {
			$frame = 3;
		}
	}
	
  # $sequence =~ tr/T/U/; # dont need this; change codon_table
  my $start_point = $frame - 1;
  my $protein_sequence="";
  for (my $i = $start_point; $i < $seq_length; $i+=3) {
      my $codon = substr($sequence, $i, 3); # problem here for i>seq_length-3 ?? or in caller getting +1,+2 > true len

## add frameshift detect option here?  if codon = stop_codons ; need it above other places
## need indel max to test: offby -2,-1,1,2 only to shift $i

      my $amino_acid;
      if (exists($codon_table{$codon})) {
          $amino_acid = $codon_table{$codon};
      } else {
          if (length($codon) == 3) {
              $amino_acid = 'X';
          } else {
              $amino_acid = "";
          }
      }
      $protein_sequence .= $amino_acid;
  }
  return($protein_sequence);
}

sub useSelenocysteine
{
  my ($turnon)= @_;
  $turnon=1 unless(defined $turnon);
  my $lastu= $USESelenocysteine; 
  unless($turnon) { # ?? dont reset global $USESelenocysteine .. BUT above uses,
    @stop_codons = qw(TAA TAG TGA); # or push @stops, 'TGA' unless(grep/TGA/,@stops);
    $codon_table{'TGA'} = '*';
    $currentCode = "universal"; $USESelenocysteine=0;
  } else {
    @stop_codons = grep !/TGA/, @stop_codons; #  qw(TAA TAG); # not TGA
    $codon_table{'TGA'} = 'u';
    $currentCode = "universalSelC"; $USESelenocysteine=1;
  }
  return $lastu;
}

our $USE_TGAw=0; 
sub useTGAw
{
  my ($turnon)= @_;
  $turnon=1 unless(defined $turnon);
  unless($turnon) {
    @stop_codons = qw(TAA TAG TGA); 
    $codon_table{'TGA'} = '*';
    $currentCode = "universal"; $USE_TGAw=0;
  } else {
    @stop_codons = qw(TAA TAG); # not TGA
    $codon_table{'TGA'} = 'w';
    $currentCode = "universal_TGAw"; $USE_TGAw=1;
  }
  return $turnon;
}


BEGIN {
  $currentCode = "universal"; 
  @stop_codons = qw(TAA TAG TGA);
  # stops: TAG  TGA  TAA; Note: TGA also == Selenocysteine valid non-stop
  ## add/allow N in silent positions; need other codon table w/ extended AA vals for ambiguous
  %codon_table = (    
  TTT => 'F', TTC => 'F', TTA => 'L', TTG => 'L',
  CTT => 'L', CTC => 'L', CTA => 'L', CTG => 'L',  CTN => 'L',
  ATT => 'I', ATC => 'I', ATA => 'I', ATG => 'M',
  GTT => 'V', GTC => 'V', GTA => 'V', GTG => 'V',  GTN => 'V',
  TCT => 'S', TCC => 'S', TCA => 'S', TCG => 'S',  TCN => 'S',
  CCT => 'P', CCC => 'P', CCA => 'P', CCG => 'P',  CCN => 'P',
  ACT => 'T', ACC => 'T', ACA => 'T', ACG => 'T',  ACN => 'T',
  GCT => 'A', GCC => 'A', GCA => 'A', GCG => 'A',  GCN => 'A',
  TAT => 'Y', TAC => 'Y', TAA => '*', TAG => '*',
  CAT => 'H', CAC => 'H', CAA => 'Q', CAG => 'Q',
  AAT => 'N', AAC => 'N', AAA => 'K', AAG => 'K',
  GAT => 'D', GAC => 'D', GAA => 'E', GAG => 'E',
  TGT => 'C', TGC => 'C', 
    TGA => '*',  # alternate UGA => 'U' Sec / SeC Se-Cys
  TGG => 'W',
  CGT => 'R', CGC => 'R', CGA => 'R', CGG => 'R',  CGN => 'R',
  AGT => 'S', AGC => 'S', AGA => 'R', AGG => 'R',
  GGT => 'G', GGC => 'G', GGA => 'G', GGG => 'G',  GGN => 'G', 
  );
  

  #DEBUG# $USESelenocysteine=1; #  # fixme, need to define it off? or expect caller to?
  $USESelenocysteine=0 unless(defined $USESelenocysteine); # fixme, need to define it off? or expect caller to?
  useSelenocysteine($USESelenocysteine) if($USESelenocysteine); # no?
  # useTGAw(0); # DEBUG
  
  # should use freq/codon or other choice of best codon per aa.
  # our %amino2codon_ambig,%amino2codon_fixed
  
  # use constant backtrans_P2C => 1;  
  # backtrans pep2codon from pal2nal.pl
  # ambiguous dna codes: R=AG, Y=TC, 
  %amino2codon_ambig = (
    "F" => "TTy", # (T|C|Y)
    "L" => "YTn", # "(CT.)|(TTR)", # A|G|R
    "I" => "ATn", # (T|C|Y|A)
    "M" => "ATG",
    "V" => "GTn",
    "S" => "TCn", # "(TC.)|(AGY)", # (T|C|Y)
    "P" => "CCn",
    "T" => "ACn",
    "A" => "GCn",
    "Y" => "TAy", # (T|C|Y)
    "*" => "Trn", # "(TAR)|(TGA)", # (A|G|R)
    #"_" => "(TA(A|G|R))|(TGA)",
    "H" => "CAy", # (T|C|Y)
    "Q" => "CAr", # (A|G|R)
    "N" => "AAy", # (T|C|Y)
    "K" => "AAr", # (A|G|R)
    "D" => "GAy", # (T|C|Y)
    "E" => "GAr", # (A|G|R)
    "C" => "TGy", # (T|C|Y)
    "W" => "TGG",
    "R" => "CGn", # "(CG.)|(AGR)", # (A|G|R)
    "G" => "GGn",
    "X" => "nnn",
    );

  foreach my $codon (sort keys %codon_table) {
    my $aa= $codon_table{$codon};
    $amino2codon_fixed{$aa}= $codon unless($amino2codon_fixed{$aa}); # multiple codons/aa, what to do?
  }
 
if(backtrans_AmbiguousNucs) { # default ; now opt to backtranslate_protein
  %amino2codon_table= %amino2codon_ambig;  
} else {
  %amino2codon_table= %amino2codon_fixed;  
}  
  
}

1;
