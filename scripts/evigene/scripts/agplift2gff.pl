#!/usr/bin/env perl

=head1 NAME  
  
  agpliftgff.pl

=head1 DESCRIPTION

  convert GFF locations from old to new assembly that differs
  in minor ways identified by oldAgp, newAgp scaffold assembly tables.
  after Jim Kent's agpLift
  
  *** NOTE: This version will handle only same contigs rearranged to new scaffolds

  my $ok=&GetOptions( 
  "oldagp=s" => \$oldagp,
  "newagp=s" => \$newagp,
  "dirnew=s" => \$dirnew,
  "gff=s" => \@gff,
  "test!" => \$test,
  "verbose!" => \$verbose,
    );
  
=head1 AUTHOR
  
  Don Gilbert; gilbertd near indiana edu
  
=cut  

use strict;
use Getopt::Long;
use Digest::MD5;

use constant VERSION => '20150527'; #upd 20150425, prior 2013..
my $test = 0;
my $verbose = 1;
my $refmap= 1;
my $domd5 = 1;
my $oldagp= undef;
my $newagp= undef;
my $dirnew= undef;
my $SAME_SCAFID_ARE_IDENTICAL= 0; # maybe this should be on default, should be no-change if really same 
my @gff= ();

# no longer true: warn "This ONLY handles case of changing N spacer fragments\n";
# but will handle only same contigs rearranges to new scaffolds.


my $ok=&GetOptions( 
  "oldagp=s" => \$oldagp,
  "newagp=s" => \$newagp,
  "dirnew=s" => \$dirnew,
  "gff=s" => \@gff,
  "sameIdsAreIdentical!" => \$SAME_SCAFID_ARE_IDENTICAL, ## -new_scaffold_id_same_as_old means identical agp in old/new.agp
  "test!" => \$test,
  ## "refmap!" => \$refmap,
  "verbose!" => \$verbose, # int for vv or v>1 
  );


push(@gff, @ARGV);
die "$0 -oldagp old.agp -newagp new.agp -gff genes.gff other.gff .. 
 output to .gff.lft
 -verbose : warns
 -sameIdsAreIdentical : new scaffold ids are same as old, in old/new.agp
 version ".VERSION."; See perldoc $0 for help \n"
  unless($ok && $oldagp && $newagp && ($test || @gff));

$dirnew =~ s,/$,, if($dirnew);
if($dirnew and !(-d $dirnew)) { warn "mkdir $dirnew\n"; mkdir($dirnew); }

my $agpold= readAgp($oldagp);
my $agpnew= readAgp($newagp);

foreach my $gff (@gff) {
  (my $gffnew= $gff) =~ s/\.gz//; $gffnew.=".lft";
  $gffnew =~ s,^.*/,$dirnew/, if($dirnew);
  liftByAgp($agpold,$agpnew,$gff,$gffnew);
}
exit;

#--------------------------------------------------------------------------
# index agp row [] array
use constant SCAFID => 0;
use constant SCAF_B => 1;
use constant SCAF_E => 2;
use constant AGPI   => 3;
use constant WN_TYPE=> 4;
use constant CTGID  => 5;
use constant CTG_B  => 6;
use constant CTGLEN => 7;
use constant CTG_OR => 8;

sub readAgp {
  my($agpfile)= @_;
  my %agp; 
  my($ngap,$ncontig, $ok);
  if($agpfile =~ /\.gz/ && -f $agpfile) { $ok= open(F,"gunzip -c $agpfile|"); } 
  else { $ok= open(F,$agpfile); }
  die "cant read $agpfile" unless ($ok);

## agp2 adds GAPEVD column: paired-ends, ... for N|U WN_TYPE 
## agp2 Error unknown line: 'U' == gap, asm spacer std size, same as 'N' 
## KN806039.1	101409	101508	22	U	100	scaffold	yes	paired-ends
  my @agpzero= (0) x 1;
  
## ?? bug for these scaff == 1contig entries, same new scaf==contig id??
## all are 'cant find JXM..'  ** bug below for scafid == ctdid in %agp{ scafid }[0..n] and agp{ ctgid } = data
# new0: Scaffold1464    1       31188   1       W       JXMV01047516.1  1       31188   +
# old1: JXMV01047516.1  1       31188   1       W       JXMV01047516.1  1       31188   +

  while(<F>){
    next unless(/^\w/); chomp;
    my @v= split"\t";
    # WN_TYPE == W,N,F,... what all are these codes?  add U == N variant
    my ($scafid, $scaf_b, $scaf_e, $agpi, $WN_TYPE, 
        $ctgid, $ctgstart, $ctglen, $gaps, $FRAGMENT, $YES, $PLUS, $GAPEVD);

    # if($v[4] eq 'N' or $v[4] eq 'U' )
    if($v[4] =~ /^[UN]/ or $v[6] eq 'fragment') {
      ($scafid, $scaf_b, $scaf_e, $agpi, $WN_TYPE, $gaps, $FRAGMENT, $YES, $GAPEVD)= @v;
      $agpi--; # 0-origin
      #WAS $agpr= [ $scaf_b, $scaf_e, $WN_TYPE, $gaps, $agpi, 0, $scafid, ];
      my $agpr= [ $scafid, $scaf_b, $scaf_e, $agpi, $WN_TYPE, 0, 0, $gaps, 0, 0, ];
      $agp{$scafid}= [@agpzero] unless(ref $agp{$scafid});
      $agp{$scafid}[$agpi]= $agpr;
      $ngap++;

    } elsif($v[4] =~ /^[WF]/) {
      ($scafid, $scaf_b, $scaf_e, $agpi, $WN_TYPE, $ctgid, $ctgstart, $ctglen, $PLUS)= @v;
      $agpi--; # 0-origin
      $ctgid .="_ctg"; #1505 bugfix for scafid == ctgid, dont need orig cgtid ??
      ## $PLUS should only be '+' or '-', bug if not?
      unless($PLUS eq '+' or $PLUS eq '-') { # warn??
         warn "Bad Plus/Minus='$PLUS' : $_\n";  $PLUS=0 ;
      }
      #WAS $agpr= [  $scaf_b, $scaf_e, $WN_TYPE, $ctglen, $agpi, $ctgid, $scafid, ];
      my $agpr= [ $scafid, $scaf_b, $scaf_e, $agpi, $WN_TYPE, $ctgid, $ctgstart, $ctglen, $PLUS, 0 ];
      $agp{$scafid}= [@agpzero] unless(ref $agp{$scafid});
      $agp{$scafid}[$agpi]= $agpr;
      $agp{$ctgid}= $agpr; # for crossref by ctgid
	## ** ^^ BUG when ctgid == scafid ***
      #?# $agp{$ctgid}[0]= $agpr; #bad ??? revert, above form wipes out scafid array for cid = sid
      $ncontig++;
      
    } else {
      warn "Error unknown line: $_\n"; next;
    }
  }
  warn "# read $agpfile: n_scaf=",scalar(keys(%agp))," n_contig=",$ncontig," ngap=",$ngap,"\n" if ($verbose);
  return \%agp;
}

=item bugs here wrong ctgid/ari

  .. problem seems that ar array hits gaps, cant decide which before/after contig is right w/ changed locs
  .. maybe want thse same agpi nums, or ctgids
  
     $$ar[AGPI], $$astart[AGPI]
     $$arnew[AGPI], $$arnew0[AGPI]
     
  my ($arnew0)= anewIndex($astart, $aoldr, $anew,  1, $start,$stop);
  my ($arnew) = anewIndex($ar, $aoldr, $anew, 0, $start,$stop); #1=backup or not? problems..

  my $aiold= $$ar[AGPI] - $$astart[AGPI];
  my $ainew= $$arnew[AGPI] - $$arnew0[AGPI];
  
  my $baseold= sumbase( $$ar[AGPI], $$astart[AGPI], $ref, $aold);
  my $basenew= sumbase( $$arnew[AGPI], $$arnew0[AGPI], $$arnew[SCAFID], $anew);
  baseold ne basenew is problem for lift locs
  
=cut

sub anewIndex {
  my($ar,$aoldr,$anew, $isfirst, $start,$stop)= @_;
  ## new has diff index :((((
  # my $ari= $$ar[AGPI]; my $arnew= $$anewr[$ari]; ## BAD index often
  ## $backup == isfirst
  
  my $ctgid= $$ar[CTGID];
  
  unless($ctgid) { # ar is gap row
    my $ari= $$ar[AGPI]; # bug if 0 ?
    if($isfirst and $ari>0) { # back1
      my $ar0=  $$aoldr[$ari-1];
      $ctgid= $$ar0[CTGID];
      my $anewc= $anew->{ $ctgid };
      my $anewi= $$anewc[AGPI]; my $anews= $$anewc[SCAFID];
      my $anewr= $anew->{ $anews }[$anewi + 1]; # gap following new ctg
      return $anewr if($anewr);
      
    } else { #fwd 1
      my $ar0=  $$aoldr[$ari+1];
      $ctgid= $$ar0[CTGID];
      my $anewc= $anew->{ $ctgid };
      my $anewi= $$anewc[AGPI]; my $anews= $$anewc[SCAFID];
      my $anewr= $anew->{ $anews }[$anewi - 1]; # gap before new ctg
      return $anewr if($anewr);
    }
  }
  
  # unless(0 and $ctgid) {  
  #   # ar is gap row .. see below, want to find next anew scaf, then backup to same gap spacer
  #   my $ari= $$ar[AGPI]; # bug if 0 ?
  #   ## bugs here.. always backup check 1st? need to check input start,stop span w/ ar1
  #   my $ar0= ($ari>0) ? $$aoldr[$ari-1] : 0; 
  #   my $ar1= $$aoldr[$ari+1];  
  #   if($ar0) {
  #     $ar1= $ar0 unless($ar1);
  #     if($stop) {
  #       my $anew0= $anew->{$$ar0[CTGID]};
  #       my $ani= $$anew0[AGPI];
  #       $isfirst=1 if($ani == $ari or $ani == $ari-1); # no good?
  #       $ar1= $ar0 if($isfirst); # and $$ar0[CTGID];  not all right but fewer mistakes?
  #     }
  #   }
  #   
  #   # if($isfirst){ $ar1= $$aoldr[$ari-1]; }
  #   # unless($ar1){ $ar1= $$aoldr[$ari+1]; } # next in list always is contig ?
  #   
  #   $ctgid= $$ar1[CTGID];
  # }
    
  return unless($ctgid);
  my $anewr= $anew->{$ctgid};
  return $anewr; #?? $$anewr[0]; # @($anew->{$ctgid}); #?? [0];
}

sub sumbase {
  my($agi1, $agi0, $scafid, $agpr)= @_;
  my $nb=0;
  for (my $i=$agi0; $i<=$agi1; $i++) {
    my $ar= $agpr->{$scafid}[$i];
    next if($$ar[WN_TYPE] =~ /^[UN]/);
    $nb += $$ar[CTGLEN]; # all of it or just start - stop ?    
  }
  return $nb;
}



sub findOldNewAgp {
   my($aold,$anew,$ref,$start,$stop)= @_;
   my($astart,$astop);

   ## dgg, jun07: fixme to lift from scaffold.agp to chrom.agp : diff refs **
   ## use contigids to map from scaf to chr id?

   my $aoldr = $aold->{$ref};
   # my $anewr = $anew->{$ref};
   return (undef,undef,$start,$stop) unless(defined $aoldr); #?? && defined $anewr);
   
   ## also check for -dna ends far past old scaff end .. correct here?
   # my $stopOffBy1 = $stop;
   my $aend= $$aoldr[-1]; return (undef,undef,$start,$stop) unless(ref $aend); 
   $stop  = $$aend[SCAF_E] if($$aend[SCAF_E] < $stop);
# bug: ^Can't use string ("0") as an ARRAY ref while "strict refs" evigene/scripts/agplift2gff.pl line 249
   $start = $$aend[SCAF_E] - 1 if($$aend[SCAF_E] <= $start); # got this bad case also !
   my $arlast; my $iar=0;
   foreach my $ar ( @$aoldr ) {
      #OLD $ar == [  $scaf_b, $scaf_e, $WN_TYPE, $ctglen, $agpi, $ctgid, ];
      #NEW [ $scafid, $scaf_b, $scaf_e,  $agpi, $WN_TYPE, $ctgid, $ctgstart, $ctglen, $PLUS, 0 ]
      
      #NO# next if($$ar[WN_TYPE] =~ /^[UN]/); # skip gaps here?? upd 1504
      my $isgap= ($$ar[WN_TYPE] =~ /^[UN]/);
      
      # this avo, avn may be bad now >> astart cross-contig problems?
      
      my($arb,$are)= ($$ar[SCAF_B],$$ar[SCAF_E]);
      if($start < $are) {

        if($start >= $arb && $stop <= $are) { 
          ## contained in one contig
          my $avo= [ $arb, $arb, $are, $$ar[WN_TYPE], $$ar[CTGLEN], $$ar[SCAFID], $$ar[CTG_OR]  ];
          my ($arnew)= anewIndex($ar, $aoldr, $anew);
          return (undef,undef,$start,$stop) unless($arnew);  # error ??
          my $avn= [ $$arnew[SCAF_B], $$arnew[SCAF_B], $$arnew[SCAF_E], $$arnew[WN_TYPE], $$arnew[CTGLEN], $$arnew[SCAFID], $$arnew[CTG_OR] ];
          return($avo,$avn, $start, $stop);
        } # done contained in one contig
      
        unless($astart) { 
          $astart= ($isgap) ? $arlast : $ar; # upd1504
        }
        
        if ($stop <= $are) {
          my $alen= 1 + $are - $$astart[SCAF_B]; #? is this bad 
          my $avo= [ $$astart[SCAF_B], $arb, $are, "C", $alen, $$ar[SCAFID], $$ar[CTG_OR] ];

          my($arnew0, $arnew );
if(1) { 
          # find old,new pair by contig id in agp.new hash seems right.
          my $cidb= $$astart[CTGID];        
          $arnew0= $anew->{ $cidb };

      ## Problem is when start,stop hits gap .. next agi.contig can be mistake, 
      ## allow arnew as gap span?
          my($cide,$agie); $agie=$$ar[AGPI];
          if($isgap and $arlast) { $cide= $$arlast[CTGID]; } #  $agie=$$arlast[AGPI];
          else { $cide=$$ar[CTGID];  }
          if($isgap and ($cide eq $cidb or not $cide)) { my $anext= $aoldr->[$iar+1]; $cide= $$anext[CTGID];  }
          $arnew= $anew->{ $cide };
          if($$arnew[AGPI] > $agie) {  ## bug? this solves it.. gap aligned cases w/ agp shift
            my $ain1= $$arnew[AGPI];
            my $scafn= $$arnew[SCAFID];
            my $arnewgap= $anew->{ $scafn }[$ain1 - 1];
            $arnew= $arnewgap if(ref $arnewgap and $$arnewgap[WN_TYPE] =~ /[NU]/);
          }
          
} else { # this anewIndex() seems problematic
          ($arnew0)= anewIndex($astart, $aoldr, $anew,  1, $start,$stop);
          ($arnew) = anewIndex( $ar, $aoldr, $anew, 0, $start,$stop); #1=backup or not? problems..
          ## (($isgap and $arlast ne $astart) ? $arlast : $ar)
}          
          return (undef,undef,$start,$stop) unless($arnew);  # error ??
          $arnew0= $arnew unless($arnew0); # is this ok for missed astart index?
          my $alennew= 1 + $$arnew[SCAF_E] - $$arnew0[SCAF_B];  #? is this bad 
          my $avn= [ $$arnew0[SCAF_B], $$arnew[SCAF_B], $$arnew[SCAF_E], $$arnew[WN_TYPE], $alennew, $$arnew[SCAFID], $$arnew[CTG_OR] ];

          my $aiold= $$ar[AGPI] - $$astart[AGPI];
          my $ainew= $$arnew[AGPI] - $$arnew0[AGPI];

          # check that arnew0,arnew SCAFID is same; otherwise problem
          ## MISSING arnew0 bug here..
          if($$arnew0[SCAFID] ne $$arnew[SCAFID]) {
            my $err= "Problem:split-scaf-$ref:$start-$stop,to:$$arnew0[SCAFID]:$$arnew0[SCAF_B]..$$arnew[SCAFID]:$$arnew[SCAF_E]";
            warn "\n#$err\n" if $verbose ;
            ## need to write this to gff output
            ## should pick longest scaffold part here?
            my $diff1 = $$ar[SCAF_B] - $start; 
            my $diff0 = $stop - $$astart[SCAF_E]; 
            if($diff0 > $diff1) {
            $avn= [ $$arnew0[SCAF_B], $$arnew0[SCAF_B], $$arnew0[SCAF_E], $$arnew0[WN_TYPE], $$arnew0[CTGLEN], $$arnew0[SCAFID], $$arnew0[CTG_OR] ];            
            } else {
            $avn= [ $$arnew[SCAF_B] - $diff1, $$arnew[SCAF_B], $$arnew[SCAF_E], $$arnew[WN_TYPE], $$arnew[CTGLEN], $$arnew[SCAFID], $$arnew[CTG_OR] ];
            }
            $$avn[3]= "E:".$err; # WN_TYPE

          } elsif(abs($ainew - $aiold) > 0) { 
              ### $$arnew[AGPI] > $$arnew0[AGPI] + 2 or $$arnew[AGPI] < $$arnew0[AGPI])
              # what?? 1 or 0 ; abs or not ??
            
            my $baseold= sumbase( $$ar[AGPI], $$astart[AGPI], $ref, $aold);
            my $basenew= sumbase( $$arnew[AGPI], $$arnew0[AGPI], $$arnew[SCAFID], $anew);
            if($baseold == $basenew) {
            
            } else {
            # also need to check if arnew0 - arnew are separated by other contigs now
            my $err= "Problem:contig-change-$ref:$start-$stop,"
                   . "old-ctg$$astart[AGPI]..$$ar[AGPI],w=$baseold,"
                   . "new-ctg$$arnew0[AGPI]..$$arnew[AGPI],w=$basenew";
            $$avn[3]= "E:".$err; # WN_TYPE
  ## BUG, start > stop after lift, and -strand changed when shouldnt
  ## Problem is when start,stop hits gap .. next agi.contig can be mistake, want to use gap span?
  #lft:  KN805742.1	 exon	454139	454079	96	-	.	Parent=Funhe2EKm008580t1;Split=1;trg=Funhe2Exx11m006549t1 1582 1974;ix=2;lold=KN805742.1:453687-454079:+;
  #   error=E:Problem: contigs rearranged at scaffold KN805742.1:453687-454079  contigs 138 .. 139 (nb=1398) to contigs 138 .. 140 (nb=3352)
  #orig: KN805742.1	 exon	453687	454079	96	+	.	Parent=Funhe2EKm008580t1;Split=1;trg=Funhe2Exx11m006549t1 1582 1974;ix=2
     
            my $cido= $$astart[CTGID];
            my $cidn= $$arnew0[CTGID];
            
            $err .= ", $cido/$$astart[SCAFID]:$$astart[SCAF_B] .. $$ar[SCAF_E] to "
                  . " $cidn/$$arnew0[SCAFID]:$$arnew0[SCAF_B] .. $$arnew[SCAF_E] ";
            warn "\n#$err\n"  if $verbose ;
            }
            # need to test old,new sum contig width b/n start,stop ? to see if real error?
          }
          
#Problem: contigs rearranged on scaffold: NW_001814680.1:48812-49065  index 1 .. 2 to contig index 2 .. 2
#Problem: contigs rearranged on scaffold: NW_001820749.1:425591-425794  index 40 .. 41 to contig index 40 .. 42
# ... is this problem or added NNN gap?

          return($avo,$avn, $start, $stop);
        } else {
          # in between; save what?
        }
      }
      $arlast= $ar; $iar++;
    }
    
  return (undef,undef,$start,$stop);  # error ??
}

#Problem: split scaffold: NW_001815682.1:5829929-5836330 to  GL340867.1:5819127 .. GL341023.1:3114
  # >> SCAFFOLD1	nvogs12	mRNA	5829929	5836330	.	-	.	ID=XM_001607701.1;gene_id=LOC100123965, alt=NV10340-RA
  # .. looks like split is in gap in middle intron of this gene
#Problem: split scaffold: NW_001820638.1:935599-938800 to  GL340870.1:934808 .. GL340957.1:14284
#Problem: split scaffold: NW_001820638.1:938675-938800 to  GL340870.1:934808 .. GL340957.1:14284
#Problem: split scaffold: NW_001820638.1:938675-938797 to  GL340870.1:934808 .. GL340957.1:14284
 # >> SCAFFOLD8	nvogs12	mRNA	935599	938800; ID=XM_001605254.1;gene_id=LOC100121692, alt=NV12488-RA
 # .. last 3' CDSexon of gene is split in gap;
         
                  

# sub OLD_offsetBy {
#   my($aold,$anew,$start,$stop)= @_;
#   # return($start,$stop) unless(defined $aold && defined $anew);
#   
#   # $start= $$aold[0] if($start < $$aold[0]); # fix bad data ??
#   # $stop = $$aold[2] if($stop > $$aold[2]); # fix bad data ??
#   
#   $start = $start - $$aold[0] + $$anew[0];
#     ## tricky: needs last match not first
#   $stop = $stop - $$aold[1] + $$anew[1];
#   my $refnew= $$anew[5];
#   return ($start,$stop,$refnew);
# }

# ** FIXME: need strand changes check aold, anew
# see ~/Desktop/dspp-work/genomesoft/genoperls/agpctg2chr.pl
# ($nref,$nstart,$nstop,$nstrand)= offsetBy($aold,$anew,$ref,$start,$stop,$strand);

sub offsetBy {
  my( $aold, $anew, $ref, $start, $stop, $orient)= @_;  
  ## Uhhggg; aold, anew are NOT same agp structure array
  #  but [ $$arnew0[SCAF_B], $$arnew[SCAF_B], $$arnew[SCAF_E], $$arnew[WN_TYPE], $$arnew[CTGLEN], $$arnew[SCAFID], $$arnew[CTG_OR]  ];
  use constant { Iscafbo=>0, Iscafbn=>1, Iscafen=>2, Iwtype=>3, Ictglen=>4, Iscafid=>5, Ictgor=>6 };
  
  my $refnew= $$anew[5]; ## 5==SCAFID in offset aold,anew record
  if($ref eq $refnew) {  # maybe this should never be reason to return same range..
    return($ref, $start, $stop, $orient) if($SAME_SCAFID_ARE_IDENTICAL); # no change 
  }
  my @orig=($ref, $start, $stop, $orient);
    
  ##FIXME: is orient same for old, new?
  my $ctgor0= $$aold[6]; ##NOT CTG_OR;
  my $ctgor1= $$anew[6]; ##NOT CTG_OR;
  ## bug bad data aold, anew? or==null
  
  $ref= $refnew;
  if($ctgor0 and $ctgor1 and $ctgor0 ne $ctgor1) { # reversed orient
    my $ctg_e = $$aold[2]; ## $$aold[0] + $$aold[4] ;
    $start = $ctg_e - $stop  + $$anew[1]; #?? is this right
    $stop  = $ctg_e - $start + $$anew[1]; 
    $orient= ($orient eq "-") ? "+" : "-";

  } else {
    if(0) { # E:Problem: contigs rearranged bug, diff begin,end spans cant handle
      my $bspan= $$anew[0] - $$aold[0];
      my $espan= $$anew[1] - $$aold[1];
      if($bspan != $espan) { # this is allowable? dont change spans? but need to fix start > stop
        # warn "#Problem: lift $ref bspan:$bspan ne espan:$espan\n"; # if $verbose ;
        # if($bspan>$espan) { $bspan=$espan; } else { $espan=$bspan; }
      } 
      $start += $bspan;
      $stop  += $espan;
    } else { # old way
      $start = $start - $$aold[0] + $$anew[0];
      $stop  = $stop  - $$aold[1] + $$anew[1];
    }
    
      ## bad nums here, start == -1, only when orient=- ? stop ok
      ## tricky: needs last match not first
  }

  warn "#Problem: lift $ref stop:$stop < start:$start\n" if($stop < $start); # if $verbose ;
  # return @orig if($stop < $start);
  return($ref, $start, $stop, $orient);
}



sub liftByAgp {
  my($agpold,$agpnew,$gff,$gffnew)= @_;
 
  warn "# liftByAgp $gff\n" if ($verbose);
  my ($nchange,$nread,$ok,$nmissed,$nfound,$nproblem);
  if($gff =~ /\.gz/ && -f $gff) { $ok= open(F,"gunzip -c $gff|"); } 
  else { $ok= open(F,$gff); }
  die "cant read $gff" unless ($ok);
  open(OUT,">$gffnew") or die "cant write $gffnew";
  my $outh= *OUT;
  while(<F>){
    unless(/^\w/){ print $outh $_; next; }
    my @v= split"\t";
    my($ref,$gffsource,$type,$start,$stop,$eval,$strand,$offs,$attr)= @v;
    
    ## need to correct some bad endpoints past old scaffolds ?* -dmel-dna aligns
    my($aold,$anew,$start,$stop)= findOldNewAgp($agpold,$agpnew,$ref,$start,$stop);
    my($nref,$nstart,$nstop,$nstrand)= ($ref,$start,$stop,$strand);
    
    if(defined $aold && defined $anew) {
      ($nref,$nstart,$nstop,$nstrand)= offsetBy($aold,$anew,$ref,$start,$stop,$strand); # check for errs in offset? start>stop ..
      $nfound++;
    } else {
      warn "# cant find $ref:$start,$stop\n" ;
      $nref="old$ref";  $nmissed++;
    }

    my $changed=($nstart != $start || $nstop != $stop || $nstrand ne $strand)?1:0;
    $nchange++ if($changed);
    $nread++;
    
    $attr =~ s/$/;lold=$ref:$start-$stop:$strand/ if($changed and $verbose); # only if lold ne lnew
    ## add problem comments from $anew[nnn]
    if($$anew[3] =~ /^E/) { $attr =~ s/$/;agperror=$$anew[3]/;  $nproblem++; }
    
    print $outh join("\t", $nref,$gffsource,$type,$nstart,$nstop,$eval,$nstrand,$offs,$attr);
    if($verbose && $nread % 1000 == 0){ print STDERR "."; print STDERR "\n" if($nread % 50000 == 0); }
  }
  close(F); close(OUT);
  warn "\n# wrote $gffnew: nfound=",$nfound," nchange=",$nchange,
    " nproblem=",$nproblem," nmissed=",$nmissed," nread=",$nread,"\n" if ($verbose);
}




__END__

# AGP OLD
scaffold_1696	1	2098	1	W	contig_1918	1	2098	+
scaffold_1696	2099	2130	2	N	32	fragment	yes
scaffold_1696	2131	7778	3	W	contig_1919	1	5648	+
scaffold_1696	7779	7803	4	N	25	fragment	yes
scaffold_1696	7804	12220	5	W	contig_1920	1	4417	+

# AGP NEW
scaffold_1696	1	2098	1	W	contig_1918	1	2098	+
scaffold_1696	2099	2123	2	N	25	fragment	yes	    ** lost 7
scaffold_1696	2124	7771	3	W	contig_1919	1	5648	+
scaffold_1696	7772	7839	4	N	68	fragment	yes	    ** gained cum. 36
scaffold_1696	7840	12256	5	W	contig_1920	1	4417	+

*** 535,546 **** OLD
  scaffold_1683 DGIL_SNP        exon    1985    3136    15.865  -       .       Parent=GG_DGIL_SNP_28100197
  scaffold_1696 DGIL_SNP        gene    854     1633    .       -       .       ID=GG_DGIL_SNP_28100198
  scaffold_1696 DGIL_SNP        exon    854     1633    23.557  -       .       Parent=GG_DGIL_SNP_28100198
! scaffold_1696 DGIL_SNP        gene    4543    5259    .       +       .       ID=GG_DGIL_SNP_28100199
! scaffold_1696 DGIL_SNP        exon    4543    5259    -29.981 +       .       Parent=GG_DGIL_SNP_28100199
! scaffold_1696 DGIL_SNP        gene    9287    10111   .       +       .       ID=GG_DGIL_SNP_28100200
! scaffold_1696 DGIL_SNP        exon    9287    9369    -2.787  +       .       Parent=GG_DGIL_SNP_28100200
! scaffold_1696 DGIL_SNP        exon    9440    9725    3.004   +       .       Parent=GG_DGIL_SNP_28100200
! scaffold_1696 DGIL_SNP        exon    10079   10111   10.237  +       .       Parent=GG_DGIL_SNP_28100200
  scaffold_1704 DGIL_SNP        gene    911     1114    .       -       .       ID=GG_DGIL_SNP_28100201
  scaffold_1704 DGIL_SNP        exon    911     1114    2.787   -       .       Parent=GG_DGIL_SNP_28100201
  scaffold_1705 DGIL_SNP        gene    358     807     .       -       .       ID=GG_DGIL_SNP_28100202
--- 535,546 ---- NEW
  scaffold_1683 DGIL_SNP        exon    1985    3136    15.865  -       .       Parent=GG_DGIL_SNP_28100197
  scaffold_1696 DGIL_SNP        gene    854     1633    .       -       .       ID=GG_DGIL_SNP_28100198
  scaffold_1696 DGIL_SNP        exon    854     1633    23.557  -       .       Parent=GG_DGIL_SNP_28100198
! scaffold_1696 DGIL_SNP        gene    4536    5252    .       +       .       ID=GG_DGIL_SNP_28100199
! scaffold_1696 DGIL_SNP        exon    4536    5252    -29.981 +       .       Parent=GG_DGIL_SNP_28100199
   .. lost 7 ^^
! scaffold_1696 DGIL_SNP        gene    9323    10147   .       +       .       ID=GG_DGIL_SNP_28100200
! scaffold_1696 DGIL_SNP        exon    9323    9405    -2.787  +       .       Parent=GG_DGIL_SNP_28100200
   .. gained 36
! scaffold_1696 DGIL_SNP        exon    9476    9761    3.004   +       .       Parent=GG_DGIL_SNP_28100200
! scaffold_1696 DGIL_SNP        exon    10115   10147   10.237  +       .       Parent=GG_DGIL_SNP_28100200
  scaffold_1704 DGIL_SNP        gene    911     1114    .       -       .       ID=GG_DGIL_SNP_28100201
  scaffold_1704 DGIL_SNP        exon    911     1114    2.787   -       .       Parent=GG_DGIL_SNP_28100201
  scaffold_1705 DGIL_SNP        gene    358     807     .       -       .       ID=GG_DGIL_SNP_28100202
