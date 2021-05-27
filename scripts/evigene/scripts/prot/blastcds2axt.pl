#!/usr/bin/env perl
# blastcds2axt.pl : convert blastn -q xxx.cds -db other to axt align format

=item about blastcds2axt

  perl blastcds2axt.pl anofunevg24m_minimus.dcblast > anofun0875t1cds_amacu.dcblast.axt
  see KaKs_Calculator1.2/2.0: ./KaKs_Calculator -i example.axt -o testex.kaks

  $nbin/blastn -evalue 1e-9 -task dc-megablast -template_type coding -template_length 18 \
    -db ano_${spp} -query anofunevg24m.cds -out anofunevg24m_$spp.dcblast

=cut

my $OVSLOP = $ENV{ovslop}||39; # what? not great, should use blast out order, keep all but major overlap
my $pOVSLOP= $ENV{povslop}||0.49; # use this for now
my $SKIPALT= $ENV{skipalt}||0; # skip query alts w/ ID pat > =~ /t1$/
# FIXME: upd skipalt skips if both query/targ gids are >t1, misses some valid matches (cross species)
#   change SKIPALT to (a) skip only query.t2+, a/o (b) keep best query.t1 x targ.tany
# already have this: didq# my $ONEQUERY=1; # skip 2ndary aligns, eg. db = all-alt cds db
my $ONEONLY=$ENV{oneonly}||0; #  self-cds can have  2+ paralogs, keep all aligns >= MINDUP*maxaln
my $MINDUP=$ENV{mindup}||0.50;
my(%didq,%didqr);

# FIXME4: multi-aligns, pick only 1st?? check for query parts? doesnt puta() handle that?
while(<>) {
  if(/^Query=\s*(\S+)/){ $qd=cleanid($1); $qlen=0; } 
  elsif(/^Length=(\d+)/){ $qlen=$1 unless($qlen); } 
  elsif(/^Lambda/) { if(%al) { puta(); putb($qd,$rd,$qlen); } $rd=$ina=0; } 
   # fix mult align here? putb() only first rd? ignore splits?
  elsif(/^> (\S+)/) { if(%al) { puta(); putb($qd,$rd,$qlen); } $rd=cleanid($1); $ina=1; }
  elsif($ina) { 
    if(/^(Query|Sbjct)/) { $st=$1;
      my($sc,$bi,$al,$be)=split; 
      $bb{$st}=$bi if($bb{$st}<1); $ee{$st}=$be;
      $al{$st}.=$al; 
    } elsif(/^\s*Score =\s+(\d+)/) { 
      puta() if(%al); $bs=1; $ina++; 
    } 
  } 
}

sub cleanid { my $d=$_[0]; $d=~s/^gi\|\d+\|//; $d=~s/\|$//; $d=~s/^(\w+)\|/$1:/; $d=~s/[^\w\.]/_/g; $d; } 
# sub cleanidb { 
#   unless( $_[0] =~ s/$faprefix//) { $_[0] =~ s/^gi\|\d+\|(\S+)/$1/; }  # drop gi nums for ids
#   $_[0] =~ s/\|$//; ## some ids have '|' at end of id; chomp
#   $_[0] =~ s/\|/:/; #?? change pipe to db:id format
# }

sub puta { 
  return unless(%al);
  my($aq,$qb,$qe, $as,$rb,$re)= map{ ($al{$_},$bb{$_},$ee{$_}) } qw(Query Sbjct); 
  $ia||=0; my $or=0; 
  if($qb>$qe) { $or++; ($qb,$qe)=($qe,$qb);}
  if($rb>$re) { $or++; ($rb,$re)=($re,$rb);} 
  $or=($or==1)?"-":"+";
  $bal{$qb}= [$qb,$qe,"$rd:$rb:$re:$or",$aq,$as]
    unless($bal{$qb}); 
  if(0) {
    print "x${ia}_${qd} $qb $qe $rd $rb $re $or\n";  
    print "$aq\n$as\n\n"; } 
  $ia++;
  %al=%bb=%ee=(); 
} 

sub _min { return($_[0]<$_[1])?$_[0]:$_[1]; }
sub _max { return($_[0]>$_[1])?$_[0]:$_[1]; }

sub isover {
  my($qb,$qe,$spans)=@_;
  for my $xbe (@$spans) { 
    my($xb,$xe)=@$xbe;
    if($qb<$xe and $qe>$xb) {
      my $qw=$qe - $qb; my $ov= _min($xe,$qe) - _max($xb,$qb);
      return 1 if($ov/$qw > $pOVSLOP);
      # return 1 unless($qb+$OVSLOP > $xe or $qe-$OVSLOP < $xb);
    }
  }
  return 0;
}

## include fixup for kakscalc bug: saq,sas length must be n % 3 == 0 codon-size whole num
## also check, skip  overlap/dup aligns
## FIXME: trim to codon start, b % 3 == 0
## FIXME3: trim out '-' aligns in query (and same bases in ref), these are not CDS parts

sub putb {
  my($qd,$rd,$qlen)=@_;
  my($lb,$le,$bb,$ml,$nskip)=(0) x 9;
  my($hd,$saq,$sas)=("") x 9;
  my @alnspan;
  $qlen||=0; 
  if($qd eq $rd) { %bal=(); return 0; } # SKIPSELF always, for self-cdsblast paralog tests
  # FIXME: upd skipalt skips if both query/targ gids are >t1, misses some valid matches (cross species)
  #   change SKIPALT to (a) skip only query.t2+, a/o (b) keep best query.t1 x targ.tany
  ## SKIPALT bug2: 1st blastin not best align (need all hsp, not best/partial hsp).
    # use: $didq{$qd} = $oaln
  my($qg,$qi)= ($qd=~m/(\w+)t(\d+)$/) ? ($1,$2) : ($qd,1); 
  my($rg,$ri)= ($rd=~m/(\w+)t(\d+)$/) ? ($1,$2) : ($rd,1); 
  if(1) { #($SKIPALT) .. SKIPSELF-ALT also
    if($qg eq $rg or ($SKIPALT and ($qi>1))) { %bal=(); return 0; } 
    #ob if($qg eq $rg or ($SKIPALT and ($qi>1 or $didqr{$qg.$rg}))) { %bal=(); return 0; } 
    #o if($qg eq $rg or ($SKIPALT and ($qi>1 or $ri>1))) { %bal=(); return 0; } 
  }
  ## below now
  # if(my $oaln=$didq{$qd}) { if($ONEONLY or $aln < $MINDUP*$oaln) { %bal=(); return 0; } }

  for my $i (sort{$a<=>$b}keys %bal) { 
    my $bx=$bal{$i}; 
    my($qb,$qe,$rloc,$aq,$as)=@$bx; my $d=0; 
    if(isover($qb,$qe,\@alnspan)) { $nskip++; next; }
    if($le>0 and $qb<=$le) { $d=1+$le-$qb; } 
    elsif($le>0 and $qb>$le+1) { $d=$qb - (1+$le); $d= 3 - ($d % 3); }
    if($d > 0) { $as=substr($as,$d); $aq=substr($aq,$d); $qb+=$d; }
    $sas.=$as; $saq.=$aq; $hd.="$qb:$qe/$rloc,"; 
    $bb=$qb if($bb==0);
    ($lb,$le)=($qb,$qe); push @alnspan, [$qb,$qe]; 
  } 

  # trim cds.query '-' align crap; after size checks?
  if($saq =~ m/\-/) {
    my(@nq,@ns); my @saq=split"",$saq; my @sas=split"",$sas;
    for(my $i=0; $i<=$#saq; $i++) { if($saq[$i] ne '-') { push @nq,$saq[$i]; push @ns,$sas[$i]; } }
    $saq=join"",@nq; $sas=join"",@ns;
    # while( my $i= rindex($saq,'-') >= 0 ) { #?? bad
    #  $saq=substr($saq,0,$i).substr($saq,$i+1); 
    #  $sas=substr($sas,0,$i).substr($sas,$i+1); 
    #}
  }
  
  $ml= ($bb-1) % 3;
  if( $ml > 0) { 
    my $d= 3-$ml; map{ $_= substr($_,$d); } ($saq,$sas); $bb+=$d;  # trim to codon start
  }
  my($lq,$ls)= map{length($_)} ($saq,$sas);
  if($lq ne $ls) {
    if($ls < $lq) { $lq=$ls; $saq=substr($saq,0,$lq); } 
    elsif($ls > $lq) { $sas=substr($sas,0,$lq); } 
  } 
  unless( ($ml= $lq % 3) == 0) { $lq-=$ml; map{ $_= substr($_,0,$lq); } ($saq,$sas); } 
  $hd=~s/,$//;
  my $aln= $saq =~ tr/A-Za-z/A-Za-z/;
  my $paln= ($qlen<1)?0:int(100*$aln/$qlen);
  my $oaln=$didq{$qd}||0;
  %bal=(); 
    # FIXME: self-cds can have  2+ paralogs, didq no good... need did-qd-rd?
  if($oaln) { return 0 if($ONEONLY or $aln < $MINDUP*$oaln); }
  print "$qd\tlen=$qlen,aln=$aln,$paln%\t$hd\n";
  print "$saq\n$sas\n\n";
  $didq{$qd}= $aln if($aln > $oaln);  #old: $didq{$qd}++;
  $didqr{$qg.$rg}= $aln;
  ## DANG missing saq or sas in large subset..
}  
 
__END__

=item output

Anofunz4kEVm009125t2 1:420/mini:KB664255.1:5702819:5703230:-,
AATTTCATTCGCATAATCAACACCAATTGCTGCAATCTGTCCGTGGTGCTTTGCTGCATCGAATTGCAGTACACGTTGCGGGttt
 tttttGTTGGTAAGGCTGGTTGCACAATGAGAAAGTTGGTGGTGATACTAAAATCCGAAAGCGATAACGCGGACAACTATGGCA
 CACTGCTGGAGAAGCACGGTTTCGTACCCGTGTTCATACCGACGCTAGATTTTTGCTTCAAAAATTTGGAAGTGTTGCGGGATC
 ACTTACTGTCGCCTTACAAGTACTCAGGCCTTATCTTTACCAGTCCCCGAAGTATAACGGCCGTCCGCGATGCGGTGCAGGGGC
 AGAAGCTGAAGGACGACTGGAAAACGCTGGAAAACTACAGCGTCGGTGAGACGTCACGGGAGCTAATCCAGCGAACGCTTGAC
AATTTCGTTCGCATAATCAACACCAATTGCTGCAATCTGTCCGTGGTGCTTTGCTGCATCGAATTGCAGTACACGTTGCGGGTTT
 --------GGTAAGGCTAGTTCCACGATGAGGAAGGTAGTGGTGATACTAAAATCCGAAAGCGACAATGCGGACAACTACGGCA
 CACTGCTGGAGAAGCACGGTTTCGTGCCCGTGTTCATTCCGACGCTAGATTTTTCCTTCAAAAACTTGGAAGTGTTGCGGGATC
 ACTTACTGTCGCCTTACAAGTACTCAGGCCTTATCTTTACCAGCCCCCGAAGCATAACGGCCGTGCGCGATGCGGTACAGGGCC
 AGAAGCTGAAGGACGATTGGAAAACGCTGGAAAACTACAGCGTTGGTGAGACGTCACGGGAGCTAATCCAGCGGACGCTTGAC

=cut
