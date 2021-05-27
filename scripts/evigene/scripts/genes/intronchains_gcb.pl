#!/usr/bin/env perl
# intronchains_gcb.pl 

=item about

  intronchains_gcb : calculate common intron chains of genes (mrna) blastn located on chr assembly
  
  derives from parts of evigene/scripts:
    makeblastscore3.pl : blast table parsing (input loop, subs sumscore,bestscore,tallOutput)
    
    altbest2.pl : intron chains parsing
    
    rnaseq/asmrna_altreclass.pl : related to or may add to from here,
      reclassify best/lesser alts per locus (after locus classified)
      from gene qualities (aaqual, aahomology, chrmap,..)
    
    rnaseq/asmrna_dupfilter  : may revise to take inchains as major input of locus identity
      along with self-cds/rna align table; replace gene.gff > eqgene input? simpler calc?
    
  rename : what?
  
=cut

use strict;
use Getopt::Long;  

use constant VERSION => '2016.07.24'; # new, derived

use constant DUPSPANS => 0; # paralog near equal alt mappings, need?

my $ONEGENOME="g1"; # only this case
my $debug= $ENV{debug}||0;
my $OVSLOP=6; # FIXME: need overlapslop ~ < 0.01 of max part span ; or < 0.02..0.05 of min part span?
my $pctOVSLOP=$ENV{pctover} || 0.01; # need opt?

my $pMINLOW   = $ENV{pmin} || 0.70; # use to block 2ndary other-locus lower hits, was 0.3 
my $pMINIDENT = $ENV{minident}||0; # add 15.12, as percent of 100
my $pMINBEST2 = 0.75; # ignore here? for TBEST if below 0.5*topscore
my $pDUPEXON  = 0.95; # or higher?

my $MINSPLICES= $ENV{minsplice} || 6; # 4/6? option ! but adjust down for limited nsplices/chain
  ## minsplice tradeoff seems to be more accurate + more nolocus versus fewer nolocus, mostly right loc, but more wrong/new splice loci
  
my $skipho= $ENV{skipho} || ""; # pattern to skip homol match, e.g. aphid x swiss _ACYPI = same spp
my $keepho= $ENV{keepho} || ""; # pattern to skip homol match, e.g. aphid x swiss _ACYPI = same spp
my $onlyquery= $ENV{onlyquery} || ""; # species query to keep, or all
my $skipquery= $ENV{skipquery} || "";

my $swapqt= $ENV{swap} || 0;
my $geneaa= $ENV{aa1} || $ENV{aa} ||"";
my $tgeneaa= $ENV{aa2} ||"";
my $showSPAN= 0; # extended exons list
my $CDSSPAN=0; my $CDSTRIM=0;
my $NEEDLEN=0; 
my $INCHAINS=1; # not an option here; 
my $NOLOCUSID= 999999;
my $LOCUSMOD= $ENV{modlocus}||5; # 2..10 now works, make default ~5, eg slop; adjust base locations to common modulus to remove some blast fuzziness
# my $TALL = 1; # unused, only this for this app, 

my $optok= GetOptions(
  "sizes|aasize|aa1size=s", \$geneaa,  "aa2size=s", \$tgeneaa, 
  "skipho=s", \$skipho, "keepho=s", \$keepho, 
  "skipquery=s", \$skipquery, "onlyquery=s", \$onlyquery, 
  "swapqt!", \$swapqt,
  "INCHAINS!", \$INCHAINS,
  "skipnolen|NEEDLEN!", \$NEEDLEN,
  "pctover=s", \$pctOVSLOP, 
  "pMINLOW=s", \$pMINLOW,  # not used
  "pIDENTMIN=s", \$pMINIDENT, 
  "MINSPLICES=s", \$MINSPLICES, 
  "MODLOCUS=s", \$LOCUSMOD, 
  "showSPAN|spans:s", \$showSPAN,
  "CDSSPAN!", \$CDSSPAN, "CDSTRIM!", \$CDSTRIM,  
  "debug!", \$debug, 
  # "altids=s", \$altids, "ALTKEY=s", \$ALTKEY, 
  # "gapsizes!", \$AAGAP, # changed this opt name
  # "TALL!", \$TALL,  "align!", \$ADDALIGN, #?? only -notall
  # "oneref|ONEGENOME:s", \$ONEGENOME,
  # "TPRE=s", \$TPRE, "TBEST=s", \$TBEST, 
  );

my( $bother)= @ARGV; # no bself here..
# ( $bother, $bself) = ( $bself, $bother) if($bother =~ /self/);

die "usage: intronchains.pl genes-chr.blastn > genes.inchains \n"
  unless($optok and $bother); # or -f $bself));

$pMINIDENT= $pMINIDENT*100 if($pMINIDENT <= 1);   # percent
$pctOVSLOP= $pctOVSLOP/100 if($pctOVSLOP > 0.99); # proportion
$pMINLOW  = $pMINLOW/100 if($pMINLOW > 0.99);       # proportion
$showSPAN=1 if($INCHAINS);
$CDSSPAN=1 if($CDSTRIM);

#$ONEGENOME="g1" if(defined($ONEGENOME) and not($ONEGENOME=~/\w/));
map{ s/[,\s]+/\|/g; } ($skipho,$keepho,$onlyquery,$skipquery);

#------------
  
my (%bspans, @bspanlist, %tother, %tspan, %dupspans, %dupspan, $lq, $bmax); # %blen, 
my (%inchr,%ingene,%ingeneall,%geneinfo);

sub MAINstub {}

my($haveqlen,$sizeh,$trsizeh,$cdspanh)= readSizes($geneaa,$tgeneaa); #($ENV{aa}||$ENV{sizes});  
my %blen= %$sizeh;  
my $hascdspan= scalar(%$cdspanh);

# my($naltid,$altidh)= readAlts($altids);# 

# main loop, read/process blast table
my ($nbest,$blerr)=(0,0);
if(1) {
  my $QTcheck=($skipquery or $onlyquery or $skipho or $keepho)?1:0;

  my $inh= undef;
  if($bother =~ /^(stdin|\-)$/i) { $inh= *STDIN; }
  else { my $cmd= ($bother =~ /\.gz/) ? "gunzip -c $bother |" : $bother; # or STDIN
    open(GSCORE, $cmd) or die "# ERROR: $cmd\n"; # is it file or list? 
    $inh= *GSCORE; }

# FIXME: loosing some valid IDs have hits in blastn input .. where?
#  eg: Zeamay4EVm000002t8|Zeamay4EVm000015t10|Zeamay4EVm000027t1
## input evg4corn_pub.mrna-genoasm.blastn.gz
# Query: Zeamay4EVm000027t1 type=mRNA; Name=ILITYHIA; Dbxref=Sobic.001G110600.1.p,102; aalen=2692,93%,complete; clen=8678; offs=16-8094; oid=cornlo2m9slvelvk105Loc5729t1; organism=Zea_mays;
# Zeamay4EVm000027t1      NC_024462.1     100.00  414     0       0       8265    8678    238491899       238492312       0.0       765
# Zeamay4EVm000027t1      NC_024462.1     100.00  386     0       0       649     1034    238415028       238415413       0.0       713
# Zeamay4EVm000027t1      NC_024462.1     100.00  376     0       0       144     519     238414300       238414675       0.0       695
# Zeamay4EVm000027t1      NC_024462.1     100.00  301     0       0       2682    2982    238472776       238473076       4e-155    556
# Zeamay4EVm000027t1      NC_024462.1     100.00  215     0       0       8052    8266    238491587       238491801       2e-107    398
## ok for evg4cornmrna-genoasm.inv4chainct.tab, cds.blastn,
## miss for evg4cornmrna-genoasm.inv[46]chainct.tab, mrna.blastn

  %bspans=();  @bspanlist=(); ## replace %bspans w/ @bspanlist, but for split chr list
  %dupspans=(); $bmax=0; $lq="";
  while(<$inh>) { 
    unless(/^\w/) { 
     # if(/^# Query: (\S+)/) { my $q=$1; my($al)=m/len=(\d+)/; $blen{$q}=$al if($al); $haveqlen++ if($al); }
     next; } 

    chomp; my @v=split"\t"; 
    unless(@v==12){ warn"ERR: blasttab not 12 cols:'@v'\n"; $blerr++; die if ($blerr>9); next; }
    
    #all# my($q,$t,$pi,$aln,$mis,$indl,$rb,$re,$tb,$te,$ev,$bits)=@v; # blast table columns, outfmt=6/7
    my($q,$t,$bits,$pidn,$aln,$mis,@bspan)=  @v[0,1,11,2,3,4, 6,7,8,9]; # 6-9 =  q. start, q. end, s. start, s. end, 
    # @bspan = ($qb,$qe,$sb,$se);
    
    if($swapqt) { ($q,$t)= ($t,$q); @bspan= @bspan[2,3,0,1];}
    # $t=~s/^gnl\|CDD\|/CDD:/; # hack fix CDD:
    
    $ingeneall{$q}++; # bug hunt
    next if($NEEDLEN and not( $blen{$q} )); #drop: and ($ONEGENOME or $blen{$t}) ));
    
    next if($pMINIDENT and $pidn < $pMINIDENT); # percentage of 100
    # FIXME: ^ pMINIDENT should not apply here to all, or else change to pMINLOW, below
    # .. use 2nd filter to block non-ONEGENOME chr/scaff 2ndary matches, possibly highish ident,
    # .. after aquiring first/top chr/scaf locus hits, collect further same-scaf hits, possible lower qual
    
    
    if($QTcheck) {
    next if($skipquery and $q =~ m/$skipquery/);  
    next if($onlyquery and $q !~ m/$onlyquery/);
    next if($skipho and $t =~ m/$skipho/);
    next if($keepho and $t !~ m/$keepho/);
    }
    
    if($lq and $q ne $lq) {
      my($lbits,$lt,$maxa,$maxi)= bestscore($lq); $nbest++;
      # my($tbits, $tt)= split",", $bother{$lq};
      # $bother{$lq}="$lbits,$lt,$maxi,$maxa" if($lt and ($lbits > $tbits));
         # or ($TPRE and $lbits > $pMINBEST2 * $tbits));  # dont change for same score
      %bspans=(); @bspanlist=(); %dupspans=(); $bmax=0; $lq="";
    }
  
    $bits= bint($bits);
    # if($CDSTRIM and my $cspan= $cdspanh->{$q}) .. see sumscore()
    # BUT need to change tspan, .. other vals to match cdstrim
    
    if(1) { ## if($q ne $t) # only this?
    
      # my $aident= _max(0,$aln-$mis);    # other way to calc, maybe better? $pctident * $aln
      my $aident= int($aln * $pidn/100);  # other way 
      
      my $dosum= ($bmax<1)? 1:0;
      unless($dosum) {
        my $havespan= exists($bspans{$t}); # ignore ONEGENOME here.. 
        $dosum= ($havespan or $bits > $pMINLOW * $bmax)?1:0;
        #o# my $havespan=($ONEGENOME)?1:exists($bspans{$t});
        ## use to collect same-location lower qual hits
        ## use pMINLOW to block other-chr/scaf/locus lower qual hits
      }
        
      my $didsum= sumscore( $q, $t, $bits,$aln,$aident, @bspan) if($dosum);
        # if($havespan or $bits > $pMINLOW * $bmax); # dont collect from new chr/scaf unless highish qual
        
      # $bother{$q}="$bits,$t,$aident,$aln" unless($bother{$q}); 
      if($didsum) { $bmax= $bits if($bits > $bmax); }
    }
     
  
    $lq= $q; 
   } close($inh);

  my($lbits,$lt,$maxa,$maxi)= bestscore($lq); $nbest++;
  warn  "# scored n=$nbest from $bother\n" if $debug;
  #my($tbits, $tt)= split",", $bother{$lq};
  #$bother{$lq}="$lbits,$lt,$maxi,$maxa" if($lt and ($lbits > $tbits));
    # or ($TPRE and $lbits > $pMINBEST2 * $tbits));  # dont change for same score
}      

if($nbest > $NOLOCUSID) { $NOLOCUSID= int( ($nbest+1001)/1000 ) * 1000; }

## intableOut is making these tables..
# my($inidh, $incounth, $inid2chrh)= makeChrSpliceIds(\%inchr);
# my($ichainh, $xchainh, $nspliceh)= makeInchains(\%ingene, $inidh);
# my($nloci, $locusgenes, $ingenelocus, $inidlocus, $inchainlocus)
#     = lociOfInchains($ichainh, $nspliceh, $inid2chrh);

intableOutput();

#------ subs -------------------

sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }
sub bint { my $b=shift; return ($b<0)?0 : ($b=~/e\+/)? int($b) : $b; }

sub readSizes {
  my(@inf)= @_; 
  my $nt=0;
  my (%alen,%trlen,%cdspan,%aaqual); %alen=(); 
  my $hasgap=0; # ($AAGAP)?1:0;#  && $iscount .. my $iscount=1; 
  my $hasspan=$CDSSPAN;#  collect %trlen,%cdspan ?
  
  foreach my $aaf (grep/\w/, map{ split(/[,\s]+/,$_) }@inf) {
    if(open(F,$aaf)) { my $n=0;
      while(<F>){ 
        next if(/^\W/); chomp; my($id,$aw,@ac)=split"\t"; 
        if($hasgap) { if(@ac>1 or @ac==0){ $hasgap=0; } else { $aw -= $ac[0]; }} 
        ## dang new cds.qual has Code/Noncode col before cdspan col ..
        if($hasspan) { 
          if(@ac>=3 and $ac[2]=~/^\d/){ 
	          my $isutrorf=(m/utrorf/)?1:0; # key in $oid may be missing
	          my($gp,$aq,$tw,$csp,$cspx)=@ac; # csp == span OR Code/Noncode col ..
	          $csp=$cspx=""  if($isutrorf); # data bug: bad cds-offset for utrorfs, fixme
            $aaqual{$id}=$aq; $trlen{$id}=$tw;
            $cdspan{$id}= ($csp=~/^\d/) ? $csp : ($cspx=~/^\d/) ? $cspx : 0; # Code-col=3?  
            } 
          else { $hasspan=0; }
        } 
        $alen{$id}=$aw; $n++; 
      } close(F); 
      $nt+=$n; warn  "# readSizes n=$n from $aaf\n" if $debug;
    } else {
      warn "# cant read sizes from $aaf\n" ;# if $debug
    }
  }
  return($nt,\%alen,\%trlen,\%cdspan,\%aaqual);
}

sub sumscore {
  my( $q, $t, $bits, $aln, $aident, $qb,$qe,$sb,$se) = @_;
  my $or=0;
  my $ttrue=$t; $t=$ONEGENOME.$q if($ONEGENOME);  # 2015.02 ?? best spans ignoring scaffold splits ?
  if($qb > $qe) { ($qb,$qe)= ($qe,$qb); $or=-1; } # FIXME record $or in bspans
  if($sb > $se) { ($sb,$se)= ($se,$sb); $or=($or<0)?0:-1; }

  ## OPTION mRNA2CDSspan filter: if option mrna2cds=offset and have cdspan 
  ## then assume query = full mrna, filter out qb,qe that dont overlap cds start,stop
  ## CDSTRIM problem: top hit to false UTR join can/does mess up paralog locating..
  ##   eg. Zeamay4EVm000674t(3|33) = paralogs on diff chr, t33 CDS best located at NC_024466, t3 on NC_024461
  ##   but t33 has false long utr located on 2nd chr NC_024461 (next to t3),
  ## to correct this w/ mRNA aligns, need skip UTR hits, capture DUPSPANS 2ndary cds aligns

  use constant TRIM2CDS => 0;

  if($CDSTRIM and my $cspan= $cdspanh->{$q}) { #  201511; only for query ??
    #o? my($cb,$ce)= split /\D+/,$cspan; # should be 32-889:+; may be 889-32:-
    my($cb,$ce)= $cspan =~ m/(\d+)\D+(\d+)/; # should be 32-889:+; may be 889-32:-
    ## BUG HERE?, return 0 for valid data? NO, but above after return 0, set bmax=false bitscore
    if($cb>0 and $ce>$cb) { # data ok
      #?? ($cb,$ce)=($ce,$cb) if($cb>$ce);
      #bug? if($cb>$qe or $ce<$qb) { return 0; } # skip hit
      if( ($cb>$qe or $ce<$qb) and $qb > 0 and $qe > $qb) { return 0; } # skip hit
      elsif(TRIM2CDS) {
      if($cb>$qb) { my $d=$cb-$qb; $aln-=$d; $aident-=$d; $sb+=$d; $qb=$cb; }
      if($ce<$qe) { my $d=$qe-$ce; $aln-=$d; $aident-=$d; $se-=$d; $qe=$ce; }
      $sb=1 if($sb<1); $qb=1 if($qb<1); # sanity fix
      }
    }
  }

  unless(@bspanlist) { #o $bspans{$t}
    my $val=  [$qb,$qe,$sb,$se,$bits,$aln,$aident,$or,$ttrue];
    push( @bspanlist, $val); 
    #o $bspans{$t}=[];  ## need this
    #o push( @{$bspans{$t}}, $val); 
    $bspans{$ttrue}=[]; ## if($ONEGENOME); ##??  need for above test bspans{torig} for sumscore()
    return 1; }
    
  my $ov=0;
  my $qlen=1+$qe-$qb; my $slen=1+$se-$sb;
  my $qslop= _max($OVSLOP, int($pctOVSLOP*$qlen));
  # my $sslop= _max($OVSLOP, int($pctOVSLOP* _min($qlen,$slen))); # FIXMEd: 1% of large genome spans not useful; use qlen?

  my($qbnew,$qenew,$sbnew,$senew)=($qb,$qe,$sb,$se);
  foreach my $sp (@bspanlist) { #o @{$bspans{$t}}
    my($xb,$xe,$tb,$te,$xbit,$xaln,$xidn,$xor,$ttrue1)= @$sp; # my $ttrue1= $$sp[8];
    my($d,$qbtrim,$qetrim)=(0,0,0);
    if($qe < $xb or $qb > $xe) { }
    elsif($qe > $xe and $qb >= $xe - $qslop) { 
      if($qb <= $xe) { $d=($xe+1)-$qb; $qbnew=$xe+1; $sbnew= $sb+$d; } #  FIXME 201511: trim out slop so no overlap for final spans
    }
    elsif($qb < $xb and $qe <= $xb + $qslop) { 
      if($qe >= $xb) { $d= $qe-($xb+1); $qenew=$xb-1; $senew= $se-$d; } # FIXME 201511: trim out slop so no overlap for final spans
    }
    else { 
      #o# if($ONEGENOME and $bits >= $pDUPEXON*$xbit and $aln >= $pDUPEXON*$xaln and $ttrue1 ne $ttrue) 
      if(DUPSPANS and $aident >= $pDUPEXON*$xidn and $ttrue1 ne $ttrue) 
      {
        push( @{$dupspans{$t}}, [$qb,$qe,$sb,$se,$bits,$aln,$aident,$or,$ttrue, $xb,$xe]);
      }
      $ov=1; last; 
    }
    #if($ONEGENOME and $ttrue ne $ttrue1) { } # 2015.02
  }  
  unless($ov) { 
    my $val=  [$qbnew,$qenew,$sbnew,$senew,$bits,$aln,$aident,$or,$ttrue];
    push( @bspanlist, $val); 
    #o push( @{$bspans{$t}}, $val); 
    unless($bspans{$ttrue}) { $bspans{$ttrue}=[]; } # if($ONEGENOME); #?  for above test bspans{torig} for sumscore()
    return 1;
  }
  return 0;
}


sub bestscore {
  my($q)= @_;
  my($maxb,$maxt,$maxa,$maxi)= (0) x 9;
  my(%ttrues, %maxtb, %maxte);
  
  ## FIXED ONEGENOME changes:
  # $tspan{$q}{$t} >> $tspan{$q}
  # $tother{$q}{$t} >> $tother{$q}
  # dupspans{$t}, dupspan{$q}{$t} ???
  my $t=$ONEGENOME.$q; ## if($ONEGENOME); # this is FIXED

  unless(@bspanlist) { ## $bspans{$t}
    warn "# bestscore($q) = 0, err\n" if $debug;
    return ($maxb, $maxt, $maxa, $maxi) 
  }
  
  if(1) #o foreach my $t (sort keys %bspans) 
  {
    my ($tbit,$taln,$tidn,$mxb,$mxe,$mtb,$mte,$mor,$axs,$ats,$ttrue)= (0) x 19; 
    my (@axs,@ats); $axs=$ats="";

    my @spans= @bspanlist; #o @{$bspans{$t}};  
    foreach my $sp (@spans) {
      my($xb,$xe,$tb,$te,$xbit,$aln,$aident,$or,$ttrue1)= @$sp;
      $tbit += $xbit; $taln+= $aln; $tidn+= $aident;
      
      $ttrue=$ttrue1 unless($ttrue); # change this, set ttrue == commonest ttrue1
      $ttrues{$ttrue1}+=$aident; # maximal ident-align

      ## BAD HERE, see below makeChrSpliceIds
      #x if($LOCUSMOD) {  map{ $_= 1 + $LOCUSMOD*int($_/$LOCUSMOD); } ($xb,$xe,$tb,$te); }
      
      if(1) { #$ONEGENOME only this
        push @axs,"$xb-$xe"; 
        push @ats,"$ttrue1:$tb-$te"; # $showSPAN == 2 
        
        if(1) { #$INCHAINS only this
          $inchr{$ttrue1}{$tb}++; $inchr{$ttrue1}{$te}++;
          $ingene{$q}{$xb}= "$ttrue1:$tb"; #?
          $ingene{$q}{$xe}= "$ttrue1:$te"; #?
        }
          
      } else {
        push @axs,"$xb-$xe"; push @ats,"$tb-$te"; # $showSPAN == 2 
      }
      
      $mxe=$xe if($xe>$mxe); 
      $mte=$te if($te>$mte and $ttrue1 eq $ttrue);
      $mxb=$xb if($xb<$mxb or $mxb==0); 
      $mtb=$tb if(($tb<$mtb or $mtb==0) and $ttrue1 eq $ttrue);
      my $mte1=$maxte{$ttrue1}||0;  $maxte{$ttrue1}=$te if($te>$mte1);
      my $mtb1=$maxtb{$ttrue1}||0;  $maxtb{$ttrue1}=$tb if($tb<$mtb1 or $mtb1<1);
      
      $mor+=$or; #? sum of 0/-1 ?
    }
    
    ($ttrue) = sort { $ttrues{$b} <=> $ttrues{$a} } keys %ttrues; # change to max aligned chr/scaf
    $mtb=$maxtb{$ttrue}||$mtb;
    $mte=$maxte{$ttrue}||$mte;
    
    #if($TALL)  # only this
    my $nxon= @spans;
    $tother{$q}="$tbit,$tidn,$taln,$ttrue,$nxon"; #o: if(@spans);  # $tother{$q}{$t}
    
    if($showSPAN) { # always?
      my @ox= sort{ $axs[$a]<=>$axs[$b] }(0..$#axs);
      my $axs= join",",@axs[@ox]; 
      my $ats= join",",@ats[@ox];  
      my $oc=($mor<0)?'-':'+';
     
      if(1) {  # $ONEGENOME
      my $lr=""; my $allr=""; my @atsx=(); 
      for my $ts (@ats[@ox]) { 
        my($r,$v)= split /:/,$ts,2; 
        if($lr eq $r) { push @atsx,$v; } else { push @atsx,$ts; $allr.="$r,"; }
        $lr=$r; 
        } 
      $ats= join",", @atsx;  
      $allr=~s/,$//; 
      if(0 and $showSPAN<2 and $allr=~/,/) { $mtb="$allr:$mtb"; } #??? bad allr
      }

     if($showSPAN==2) { $tspan{$q}="$mxb-$mxe/$axs\t$ats:$oc"; } #?{$t} add $mxb-$mxe\t$mtb-$mte also?
     else { $tspan{$q}="$mxb-$mxe\t$mtb-$mte:$oc"; } #{$t} need flip orient ? $oc has info
    }   
   
if(DUPSPANS) {
    my @dups= ($dupspans{$t}) ? @{$dupspans{$t}} : (); ## output add  15.12; only for ONEGENOME ?
    if(@dups) { ## 201512: separate output rows? only $ONEGENOME now?
      my $lgr="";
      my($mxb,$mxe)=(0,0);
      my(@dupaxs,@dupats); 
      for my $dsp (@dups) {
        my($xb,$xe,$tb,$te,$xbit,$aln,$aident,$or,$ttrue1,$oxb,$oxe)= @$dsp; # oxb,oxe match to @sp:xb,xe
        $mxe=$xe if($xe>$mxe); $mxb=$xb if($xb<$mxb or $mxb==0); 
        push @dupaxs,"$xb-$xe"; 
        push @dupats,"$ttrue1:$tb-$te";
      }
      my $ndup=@dupaxs;
      my @ox= sort{ $dupaxs[$a]<=>$dupaxs[$b] or $a<=>$b }(0..$#dupaxs);
      my $ncut=0; if($ndup>20) { @ox=@ox[0..19]; $ncut=$ndup - @ox;}
      my $axs= join",",@dupaxs[@ox]; 
      my $ats= join",", 
        map{ my($gr)=m/^([^:]+):/; s/^$gr:// if($gr eq $lgr); $lgr=$gr; $_; } @dupats[@ox];
      if($ncut>0) { $axs.=",..n$ndup";  $ats.=",..n$ndup";}
      my $oc=($mor<0)?'-':'+';
      $dupspan{$q}="$mxb-$mxe/$axs\t$ats:$oc"; ## $dupspan{$q}{$t}
    } 
}

    if($tbit > $maxb) { $maxb=$tbit; $maxt= $t; $maxa=$taln; $maxi=$tidn;}
  }

 return($maxb, $maxt, $maxa, $maxi);
}


# OUTPUT.....................................

sub _sortGID {
  my($ag,$at)= split "t",$a,2; my($bg,$bt)=split"t",$b,2;
  return($ag cmp $bg or $at<=>$bt);
}

sub makeChrSpliceIds {
  #my($inchrh)= @_;
  my(%inid, %incount, %inid2chr);
  my $INID=0;
  for my $chr (sort keys %inchr) {
    my @cbe= sort{$a <=> $b} keys %{$inchr{$chr}};
    for(my $i=0; $i<@cbe; $i++) # for my $cb (@cbe) 
    {
      # LOCUSMOD here? look for +1,+2,.. next loc, give it same INID
      # this INID reset way may work, above mod to b,e locs fails
      # for ZeamEVM 1..99, mod=5 added 100/700 partof, reduced Ig999999 by 40
      my $cb=$cbe[$i]; 
      if($LOCUSMOD and $i>0) { 
        ++$INID unless($cb - $cbe[$i-1] <= $LOCUSMOD); 
      } else {
        ++$INID;
      }
      my $nid= 'N'.$INID;
      my $ki="$chr:$cb"; 
      $inid{$ki}= $nid; $inid2chr{$nid}= $ki; # revmap; ugh locusmod allows 2+ ki/nid
      my $ni= $inchr{$chr}{$cb}||0; 
      $incount{$ki}= $ni; $incount{$nid}+= $ni; #?
    }
  }
  warn  "# insplice ids n=$INID\n" if $debug;
  return( \%inid, \%incount, \%inid2chr);
}

sub makeInchains {
  my($ingeneh,$intronidh)= @_; #, $incounth
  my (%ichain,%xchain,%nsplice); my $nchain=0;
  # $ingeneh == \%ingene global
  foreach my $gid (sort _sortGID keys %ingene) {
    #$inchr{$ttrue1}{$tb}++; $inchr{$ttrue1}{$te}++;
    #$ingene{$q}{$xb}= "$ttrue1:$tb"; #?
    #$ingene{$q}{$xe}= "$ttrue1:$te"; #?
    my (@xchain,@inchain);
    @xchain= sort{ $a<=>$b } keys %{$ingene{$gid}};
    for my $inx (@xchain) {
      my $ki= $ingene{$gid}{$inx};
      my $inid= $intronidh->{$ki}||$ki; #*? need reverse map $inidloc{$inid}= $ki
      ##if($debug) { my $inc= $incounth->{$ki}||1; $inid.=":$inc"; }
      # push @xchain, $inx; # dont need, @inx == @xchain
      push @inchain, $inid;  
      }
    my $inchain= join",",@inchain;
    my $xchain= join",",@xchain;
    my $nsplice= @xchain;
    $ichain{$gid}= $inchain; $nchain++;
    $xchain{$gid}= $xchain;
    $nsplice{$gid}= $nsplice;
    }
  warn  "# makeInchains n=$nchain\n" if $debug;
  return (\%ichain,\%xchain,\%nsplice);
}

sub exonsof {
  my($gid)=@_;
  my @xchain= sort{ $a<=>$b } keys %{$ingene{$gid}};
  my @cchain;
  for my $inx (@xchain) {
    my $chrb= $ingene{$gid}{$inx} || "miss:0"; 
    push @cchain,$chrb;
    # my($chr,$cb)=split":",$chrb;
  }
  return(scalar(@xchain),\@xchain,\@cchain);
}

sub isover {
  my($cchain,$ochrlocs)=@_;
  my $isov=0;
  #my(%chrloc);
  #for my $chrb (@$ochrlocs) {
  #  my($chr,$cb)=split":",$chrb; $chrloc{$chr}{$cb}++;
  #}
  
  my $nc=@$cchain;
  for(my $i=0; $i<$nc; $i+=2) {# test *pairs* of xb-xe
    my $chrb= $$cchain[$i];   my($chr,$cb)=split":",$chrb;
    my $chre= $$cchain[$i+1]; my($chre,$ce)=split":",$chre;
    ## next unless($chrloc{$chr});
    ## my @ox=sort keys %{$chrloc{$chr}}; my $nox=@ox;
    my @ox=@$ochrlocs;  my $nox=@ox;
    for(my $j=0; $j<=$nox; $j+=2) {
      my($ochr,$ocb)=split":",$ox[$j];
      my($ochre,$oce)=split":",$ox[$j+1];
      next unless($ochr eq $chr);
      if($cb < $oce and $ce > $ocb) { $isov++; } # any other info? size of overlap?
    }
  }
  return($isov); # other info?
}

our $XODEB=0; 
sub exon_overlap # for weak intron chain hits, nor other
{
  my($gid,$locid,$locusgids)=@_;
  my ($isover,$overid)=(0) x 9;

  ## return 0 unless(ref($locusgids) and $ingene{$gid});
  ## $locusgids ref of  @ogids= @{$locusgenes{$locid}};
  
  if(ref($locusgids) and $ingene{$gid}) {
  my($nx,$xchain,$cchain)= exonsof($gid);
  for my $og (@$locusgids) { 
    my($nox,$oxlocs,$ochrlocs)= exonsof($og);
    my $isov= isover($cchain,$ochrlocs);
    if($isov) { $isover++; $overid=$og; last; }
  }
  }
  if(0) { # working now
  my $ogids= join",",@$locusgids;
  warn "#exon_overlap($gid,l$locid) = $isover,$overid of $ogids\n" if($debug and $XODEB++ < 9);
  }
  return($isover,$overid);
}

# ($issplit,$chr1,$inidchr1,$chrall)= is_splitlocus(\@inid, $inid2chrh, \%inidlocus); # want %inidlocus and %inchainlocus ?
sub is_splitlocus {
  my($inids, $inid2chrh, $inidlocush)= @_;
  my $issplit= 0;
  my(%locs,%chrid);
  for my $id (@$inids) {
    my $chrb= $inid2chrh->{$id} or next;
    my($chr,$cb)=split":",$chrb;
    $locs{$chr}++; # deal w/ base spans? $cb - othercb ?
    push @{$chrid{$chr}},$id;
  }
  my @chr= sort{ $locs{$b} <=> $locs{$a} or $a cmp $b } keys %locs;
  my $chr1= $chr[0];
  my $inidchr1= $chrid{$chr1};
  if(@chr>1) { $issplit=@chr; }
  return($issplit,$chr1,$inidchr1,\@chr);
}


sub lociOfInchains {
  my($ichainh, $nspliceh, $inid2chrh)= @_; 
  my(%locusgenes, %ingenelocus, %inidlocus, %inchainlocus);
  my($NEWLOCID,$noloc,$nloci,$nsplit)=(0) x 9;
  
  #above# my $MINSPLICES= $ENV{minsplice}||4; # option ! but adjust down for limited nsplices/chain

  # FIXME: split-genes are problem, should not use those to set locid of unsplit genes
  # .. affects inidlocus, need inid => chr/scaf hash map to help decide on splits
  
  # 1st? scan ichains for subchains ? hash ichains, hash @inid, ichains ?
  # problems calling locus for nsplice < ~4 ; skip those in 1st pass?

  # ** diff results: using TEMP set for 1st loop is more accurate .. not sure why, 2nd test of exonover?
  ## my (@chains); ## no %TEMPlocusgenes, %TEMPingenelocus ?? drop TEMP set
  my (@chains, %inchaingids, %TEMPlocusgenes, %TEMPingenelocus);  # keep TEMP set, more accurate
  
  # combine w/ following @gid loop
  my @gids= sort{ $nspliceh->{$b} <=> $nspliceh->{$a} or $a cmp $b }keys %$ichainh; # or %$nspliceh
  
  #o my @lgid= grep{ $nspliceh->{$_} >= $MINSPLICES } keys %$nspliceh;
  #o @lgid= sort{ $nspliceh->{$b} <=> $nspliceh->{$a} or $a cmp $b } @lgid; # sort by nsplice not length(chain)
  #dont need# my @lchains= map{ $ichainh->{$_} } @lgid;
  #o foreach my $ic (sort{ length($b)<=>length($a) or $a cmp $b } values %$ichainh) 
  #o foreach my $ic (sort{ length($b)<=>length($a) or $a cmp $b } @lchains) 
    #....... old ...
    # my @csuper= grep /$ic/, @chains; # slow for 500k chains?
    # unless(@csuper or $nsplice<5) {   # try also removing end vals of $ic
    #   my $icm=$ic; $icm=~s/^\w+,//; $icm=~s/,\w+$//; 
    #   if($icm=~/,/) { @csuper= grep /$icm/, @chains; }
    # }
    # my $csuper=(@csuper)?$csuper[0]:0;
    #....... old ...
    # foreach my $c (@chains) { my $i=index($c,$ic); if($i>=0) { if($c =~ m/$ic\b/) { $csuper=$c; last; } } }
    # unless($csuper or $nsplice<5) {
    #   # my $icm=$ic; $icm=~s/^\w+,//; $icm=~s/,\w+$//; 
    #   my $icm= join",", @inid[1..($nsplice-2)];
    #   foreach my $c (@chains) { my $i=index($c,$icm); if($i>=0) { if($c =~ m/$icm\b/) { $csuper=$c; last; }} } 
    # }
    
  
  for(my $i=0; $i < @gids; $i++)
  { 
    my $igid= $gids[$i]; 
    next if($nspliceh->{$igid} < $MINSPLICES); # filter here now
    
    ##x my $ic= $lchains[$i]; 
    my $ic= $ichainh->{$igid} or next;
    my @inid= split",", $ic; #?
    my $nsplice=@inid;
    my $minsplices=($nsplice>$MINSPLICES)?$MINSPLICES:$nsplice-2;
 
    my $locid=0;
    my $csuper=0; 
    if($nsplice>4) {
      my $icm= join",", @inid[1..($nsplice-2)];
      foreach my $c (@chains) { my $i=index($c,$icm); if($i>=0) { if($c =~ m/$icm\b/) { $csuper=$c; last; } } } 
    } else {
      foreach my $c (@chains) { my $i=index($c,$ic); if($i>=0) { if($c =~ m/$ic\b/) { $csuper=$c; last; } } }
    }
    
    if($csuper) {
      $locid= $inchainlocus{$csuper}||0; # missing locid?
      ## set locid hash vals below for all
      my ($sgid,@sgid)=sort keys %{$inchaingids{$csuper}}; $geneinfo{$igid}{partof}= $sgid if($sgid); # == altsame or altfrag class
    } 
    
    unless($locid) {
      push @chains, $ic; #? add if csuper also?
      $inchaingids{$ic}{$igid}++; #?
      
      my %inl=(); map{ my $id=$inidlocus{$_}; $inl{ $id }++ if($id); } @inid; 
      my @inloc= sort{ $inl{$b}<=>$inl{$a} or $a cmp $b } keys %inl;
      
      if(@inloc >= 1) { 
        $locid= $inloc[0]; #?? others
        if(@inloc>1) { $geneinfo{$igid}{locids}=join(",",@inloc); } 
        
        ## these cancels may include end-exons with weak 1/2 splice attachment; other test?
        if($inl{$locid} < $minsplices) { 
          my($isover,$overid)=(0,0);
          ($isover,$overid)= exon_overlap($igid, $locid, $TEMPlocusgenes{$locid}); # locusgenes not yet set
          $locid=0 unless($isover);  # cancel
          }   
       }  
      ## set locid hash vals below for all
    
      unless($locid) {
        #>> this part bad, gives locid to each distinct chain, need inid loci first
        $locid= ++$NEWLOCID;
      }
    }
    
    if($locid) {
      $inchainlocus{$ic}= $locid;  # this is final
      push @{$TEMPlocusgenes{$locid}}, $igid; # o TEMPlocusgenes
      $TEMPingenelocus{$igid}= $locid;  # old: TEMP hashes, see below gene loop for final

      my($issplit,$chr1,$inidchr1,$chrall)= is_splitlocus(\@inid, $inid2chrh, \%inidlocus); # want %inidlocus and %inchainlocus ?

      if($issplit) {
        #x $nsplit++;
        map{ $inidlocus{$_}=$locid unless($inidlocus{$_}); } @$inidchr1; 
      } else {
        map{ $inidlocus{$_}=$locid unless($inidlocus{$_}); } @inid; # what of clashes?
      }
    }
    
  }
  warn  "# lociOfInchains inchainloci=$NEWLOCID\n" if $debug;
  
  ##osort: length($ichainh->{$b}) <=> length($ichainh->{$a}) or $a cmp $b 
  ## same gid sort list as above? but for nsplice filter
   
  foreach my $gid (@gids) # nsplice-sorted @gids
  #above# sort { $nspliceh->{$b} <=> $nspliceh->{$a} or $a cmp $b } keys %$ichainh) 
  {
    next if($ingenelocus{$gid}); # all MINSPLICE done above now
    
    my $locid=0;
    my $ic= $ichainh->{$gid};
    my @inid= split",", $ic; ## N111,N113,N124,..
    my $nsplice=@inid;
    my $minsplices=($nsplice>$MINSPLICES)?$MINSPLICES:$nsplice-2;

    my %inl=(); map{ if(my $id=$inidlocus{$_}) { $inl{ $id }++; } } @inid; 
    my ($inloclid,@inloclidx)= sort{ $inl{$b}<=>$inl{$a} or $a cmp $b } keys %inl;

    if(my $lid= $inchainlocus{$ic}) {  # probably not seen, above
      $locid=$lid; 
      #? check cmp $inchainlocus{$ic} lid to inidlocus max lid?
      if($inloclid and $lid ne $inloclid) { $locid= $inloclid if($inloclid<$lid); } # no effect?
      
    } else {
      # my %inl=(); map{ my $id=$inidlocus{$_}; $inl{ $id }++ if($id); } @inid; 
      # my @inloc= sort{ $inl{$b}<=>$inl{$a} or $a cmp $b } keys %inl;
      
      if($inloclid) { # if(@inloc>0)  
        $locid= $inloclid; # $inloc[0]; #?? others
        if($inl{$locid} < $minsplices) { 
          my($isover,$overid)= exon_overlap($gid, $locid, $locusgenes{$locid});
          $locid=0 unless($isover);  # cancel
          }
      }
      # FIXME: $nsplice < 4? are problem cases, shouldn't assign newlocid unless they match something longer?
      unless($locid) { 
        if($nsplice < $MINSPLICES) { # probably all here, as above loop gets all >= MINSPLICES
          $noloc++;
        } else {
          my $newlocid= ++$NEWLOCID;
          $locid= $newlocid;
        }
      }
    }
    
    # loci with @inid .. or newloc
    unless($locid) { # bug? or problem cases?
      $ingenelocus{$gid}= $NOLOCUSID;
      push @{$locusgenes{$NOLOCUSID}}, $gid;
      
    } else {   
      push @{$locusgenes{$locid}}, $gid;
      $ingenelocus{$gid}= $locid;
      $inchainlocus{$ic}= $locid;
      
      # add ingeneflags{$gid}= "$issplit,$exonover,..";
      my($issplit,$chr1,$inidchr1,$chrall)= is_splitlocus(\@inid, $inid2chrh, \%inidlocus);  
      if($issplit) {
        $nsplit++;
        map{ $inidlocus{$_}=$locid unless($inidlocus{$_}); } @$inidchr1; 
        my $chrs=join",",@$chrall; $geneinfo{$gid}{split}= "$issplit,$chrs";
      } else {
        map{ $inidlocus{$_}=$locid unless($inidlocus{$_}); } @inid; # what of clashes?
      }
    }
  }
  
  $nloci= scalar(keys %locusgenes);  # report $noloc
  warn  "# lociOfInchains nloci=$nloci/$NEWLOCID, nsplit=$nsplit, id$NOLOCUSID/nolocus=$noloc\n" if $debug;
  return($nloci, \%locusgenes, \%ingenelocus, \%inidlocus, \%inchainlocus);
}

sub xchainpairs { 
  my($xc)=@_;
  #reformat: 1,128,268,437,600,773,.. => 1-128,268-437,600-773
  my @x=split",",$xc; my @c=();
  for(my $i=0; $i<@x; $i+=2) { my($x,$y)= @x[$i,$i+1]; push @c,"$x-$y"; }
  return (@c>0)?join(",",@c):$xc;
}

sub intableOutput { 
  # inputs: %tother, %tspan, ..
  # just ho table for all q x t
  my $hasDups= (DUPSPANS and scalar(%dupspan))?1:0; #? want dupspan for paralogs w/ 2+ near eq mappings?
  my @tcols= qw(Query Source Bits Ident Align Qlen Slen Nexon);
  if($INCHAINS) { push @tcols, qw(ILocus Qexons Sintrons); }
  elsif($showSPAN) { push @tcols, qw(Qspan Sspan); }
  push @tcols, qw(Qdups Sdups) if($hasDups);
  push @tcols, "Geneinfo"; #?
  print join("\t",@tcols)."\n";

  # if($INCHAINS) .. always true
  my($inidh, $incounth, $inid2chrh)= makeChrSpliceIds(\%inchr);
  my($ichainh, $xchainh, $nspliceh)= makeInchains(\%ingene, $inidh);
  my($nloci, $locusgenes, $ingenelocus, $inidlocus, $inchainlocus)
      = lociOfInchains($ichainh, $nspliceh, $inid2chrh);

  ## FIXED ONEGENOME changes:
  # $tspan{$q}{$t} >> $tspan{$q}
  # $tother{$q}{$t} >> $tother{$q}
  # my $t=$ONEGENOME.$q; ## if($ONEGENOME); # this is FIXED
  
  # FIXME: MISSING valid gid, input blastn ok .. lost where?
  # are both these equal? both set in bestscore()
  my @gidi= keys %ingeneall; # debug
  ##my @gidi= keys %ingene;
  my @gidt= keys %tother; 
  if(scalar(@gidi) > scalar(@gidt)) {
    warn  "# error geneids: ingene n=",scalar(@gidi),", tspans n=",scalar(@gidt),"\n";
    @gidt=@gidi;
  }
  
  foreach my $gid (sort _sortGID @gidt) {
    my $inchain= $ichainh->{$gid} || "noichain";
    my $xchain = $xchainh->{$gid} || "noxchain";
    my $glocus = $ingenelocus->{$gid} || 0; ## if($glocus eq $NOLOCUSID) ..
    #? my @locusgenes= $locusgenes->{$glocus} ? @{$locusgenes->{$glocus}} : ();

    my $glocid= "Ig".sprintf "%06d",$glocus;
    ## add $gid issplit, other info to report below? intronlocus is  unreliable for nexon 1,2,3 ..
    ## some, many? of short alts to longer genes appear as detatched/partial end exon transcripts
    ## .. add exon-overlap measure of locus identity; test mrna.blastn instead of cds.blastn

    my $geneinfo="";
    if($geneinfo{$gid}) { 
      my @gk=sort keys %{$geneinfo{$gid}};
      for my $gk (@gk) { my $v=$geneinfo{$gid}{$gk}; $geneinfo.="$gk=$v;"; }
    } 
    # $geneinfo ||="noginfo"; # let be empty
    
    # this sort on tbitscore maybe bad? replace w/ tident sort? but keep output column order
    #o $tother{$gid}{$b}->[1] <=> $tother{$gid}{$a}->[1]  # ident sort, NO good, val is "bits,idn,aln,.."

    ## only 1 tg here? rewrite to drop key2
    #o foreach my $tg (sort{  $tother{$gid}{$b} <=> $tother{$gid}{$a} } keys %{$tother{$gid}} ) 
    my $tg=$ONEGENOME.$gid;   # # this is FIXED
    if(1) {
      my $bia= $tother{$gid} || "0,0,0,0,0"; ## {$tg}; 
      my($bit,$idn,$aln,$tgtrue,$nxon)=split",",$bia; 

      my $tgval= $tgtrue; ## ($ONEGENOME)? $tgtrue : $tg; # problem $tg not true val here, change? multi-t gets what? 1st or most common?
      my $qlen= $blen{$gid}||0; 
      my $tlen= $blen{$tgval}||0; 
      ## FIXME: qlen off for CDSTRIM, use what? t align span from tspan?
      
      if(0) { #?? is this right .. blen{gid} is cds-length, if input is genes.cds.qual ?
      if($CDSTRIM and my $cspan= $cdspanh->{$gid}) { #  201511; only for query ??
        my($cb,$ce)= split /\D/,$cspan; # should be 32-889:+; may be 889-32:-
        ($cb,$ce)=($ce,$cb) if($cb>$ce);
        $qlen=1+$ce-$cb; # replace blen{$q} with cspan for length column ??
        }
      }

      my @vcols=($bit,$idn,$aln,$qlen,$tlen,$nxon);
      
      if($INCHAINS) {  ## FIXED
        my $span= $tspan{$gid} || "0\t0"; ## $tspan{$gid}{$tg}
        my($xs,$gs)=split/\t/,$span;
        
        my($xsb,$xse)=split"-",$xs; 
        my $xsw=1+$xse-$xsb; $vcols[3]=$xsw if($xsw>$qlen);  # FIX for mrna align span>cds span
        
        my $xc= xchainpairs($xchain); my $ic= xchainpairs($inchain);
        if($ENV{nochain}){ $xc=""; $ic=""; } #debug simpler out
        push @vcols, $glocid, "$xs/$xc", "$gs/$ic";  #? add tspan{$gid}{$tg} "$mxb-$mxe\t$mtb-$mte:$oc"; 
      } 
      elsif($showSPAN) { my $span= $tspan{$gid} || "0\t0"; push @vcols, $span; } ## {$tg}

      if(DUPSPANS and $hasDups){ my $span=$dupspan{$gid}||"0d\t0d"; push @vcols, $span; } #? keep? {$tg}

      push @vcols, $geneinfo if($geneinfo); #?
      print join("\t",$gid,$tgval,@vcols)."\n";
    }
  }
}

__END__
