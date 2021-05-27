#!/usr/bin/env perl
# gff2aligntab.pl : tabulate alignment annots from gmap.gff and gsplign.gff
#..............
# $nam.align.tab columns
# GenomeID        gespan  geor    AQueryID        quspan  match   qlen    cov     pid     path    indels  nexon 
#	splice  aalen   offs    aamap   sense   oid     tag

# 201406 fix Split paths, for cov1 > cov2, swap; row2/split2 is now always output even for minor cov2%
## allow for some attr renames, aliases:  cov > mapCover  nexon > intron

## 201711: other annot tags from evgpub.gff: 
# mapqual=100c,98i,22869l,30x; (cov,pident,clen,nexon)
# intron=27/29; (ival/itot)

# $dosplit=1; 
my $OIDisID=$ENV{oid}||0;
my $SRCKEY=$ENV{src}||$ENV{srctag}||""; # for tag column, parse this mrna.annot
# 20151106: add .ncds to nexon col; skip (opt?) mRNA vals for qlen,nexon,offs,splice, as may be wrong
my $SKIPATTR= $ENV{noattr}||$ENV{badattr}||0; # just for exon/cds computable set

my $tMRNA=$ENV{mrnatype}||$ENV{rna}||"mRNA|ncRNA|transcribed_RNA"; # fix for ncRNA; FIXME ENV{mrnatype} better
# ^ maybe allow any type with 'RNA|rna' 
my $CDSPAN=$ENV{cdspan}||0;
my $KEEP_G=$ENV{keepg}||0;

# BEGIN
my $SKIPATK='match|nexon';
 # cant recompute cdsoffset yet: offs|cdsoff'; #? add cov?=match/qlen; qlen|cxlen << cant recompute qlen from aligned exons
my @atk=qw(match qlen cov pid path indels nexon splice aalen offs aamap sense Split oid tag); 
push @atk, "cdspan" if($CDSPAN);
## alternate oid= tags if missing, Target|trg=oid , obID= special case, maybe ncbigff: db_xref|Dbxref
my @gs= qw(gescore clen gaps chimera); @hd=grep{ not/Split/ } @atk;
## gescore is splign-only tag, but may be removed as useless.. check clen vs qlen
#x %rekey=(mapCover => cov, intron => nexon);
print join("\t","GenomeID","gespan","geor","AQueryID","quspan",@hd)."\n"; 

while(<>) {
  @v=split"\t"; 
  if($v[2] =~ /$tMRNA/) { # was /\tmRNA\t/
    #b# if($SKIPATTR) { %skipat=(); while(s/;($SKIPATK)=([^;\s]+)/;${1}old=$2/g){ $skipat{$1}=$2; } }
    
    ## 201711: other annot tags from evgpub.gff: 
    # mapqual=100c,98i,22869l,30x; (cov,pid,clen,nexon)
    # intron=27/29; (ival/itot)
    if(m/mapqual=(\d+[^;\s]+)/) { 
      my $mc=$1; my($mcov,$mpid,$mlen,$mxn)= (0);
      if($mc=~m/(\d+)[ca]/){ $mcov=$1; s/$/;cov=$mcov/; }
      if($mc=~m/(\d+)i/){ $mpid=$1; s/$/;pid=$mpid/; }
      if($mc=~m/(\d+)l/){ $mlen=$1; s/$/;qlen=$mlen/; }
      if($mc=~m/(\d+)x/){ $mxn=$1; } # s/$/;nexon=$mxn/;# no, count exon rows nxf
      if(/\bintron=(\d+),(\d+);/) { my($mival,$mint)=($1,$2);
        $mxn||=1+$mint; s,intron=\d+.\d+,intron=$mint/$mxn,; } # replace new syntax w/ old intron=nint/nexn
    }
    elsif(not m/\bcov=/ and m/mapCover=([^;\s]+)/) { 
      my $mc=$1; $mc=~s/,.*//; $mc=~s/%//; s/$/;cov=$mc/; 
    }
    
    if(/cdsoff=/ and not m/\boffs=/) { s/cdsoff=/offs=/; }
    # if(/intron=/ and not m/nexon=/) { s/intron=/nexon=/; }
    ## ^ do both?  intron=100%,42/42; nexon=22;
    ## cxlen=2952/3979 >> cdslen=2952, qlen=3979

  }
  chomp; s/^#[si]\.//; next if(/^\W/); @v=split"\t"; 
  ($gid,$src,$typ,$gb,$ge,$gv,$go,$at)=@v[0,1,2,3,4,5,6,8]; # 8 was -1
  if($typ eq "CDS") { $ncds++; my $cw=1+$ge-$gb; $cdsw+=$cw; $ocb=$gb if($ocb==0 or $gb<$ocb); $oce=$ge if($ge>$oce); }
  elsif($typ eq "exon") { $nxf++; my $xw=1+$ge-$gb; $xmatch+=$xw; }
  elsif($typ eq "intron"){ $inspl += ($gv>69)?2:($gv>44)?1:0; }
  elsif($typ =~ m/$tMRNA/) { # old: ($typ eq $tMRNA) 
    my($mid)= $at=~m/ID=([^;\s]+)/;
    $issplit=(/Split=|chimera=/)?1:0; #FIXME 18.01: evganngff drops chimera, doesnt add Split; has chim[12] and ID_C[12]
    $issplit=1 if(/;chim[12]=/ or $mid =~ /_C[12]$/);
    ($tid,$tb,$te)=m/(?:trg|Target)=(\S+) (\d+) (\d+)/; 
    unless($tid) { ($tid)=$at=~m/ID=([^;\s]+)/; $tb=1; $te=1+ $ge - $gb; }
    putg() if($ltid); 
    if($SKIPATTR) { %skipat=(); while(s/;($SKIPATK)=([^;\s]+)/;${1}old=$2/g){ $skipat{$1}=$2; } }
    $oce=$ocb=$ncds=$nxf=$inspl=$cdsw=$xmatch=0; 
    $mgb=$gb; $mge=$ge;
    
    %at= map{ my $v=($at=~m/\b$_=([^;\s]+)/)?$1:0; $_ => $v; } ("ID",@atk);
    $at{nocov}=1 unless($at=~m/cov=\d/); #^ fixme: cov=0 is valid NOPATH value, but no cov= missing
    #bad.here# ($tid,$tb,$te)=m/(?:trg|Target)=(\S+) (\d+) (\d+)/; 
    $id=$at{ID}; $tid ||= $id;  # $at{Split}=~s/^C//;
    #x %ags= map{ if($at=~m/\b$_=([^;\s]+)/) { $_=>$1; } } @gs;
    %ags=(); for my $gk (@gs){ if($at=~m/\b$gk=([^;\s]+)/) { $ags{$gk}=$1; } }
    # tag == insrc=kf2a:gmap2a5u; and/or osrc=gmap2a5u;
    $tag= ($SRCKEY and m/$SRCKEY=([^;\s]+)/) ? $1 : $SRCKEY; #was "";
    #xx my $issplign= (m/gescore=/)? 1 : (m/;clen=\d/ and not m/;qlen=\d/)?1:0;
    $at{pid}||= 99; 
    $at{oid}= $tid unless($at{oid} or $tid eq $id);
    #o# my $issplign= ($ags{gescore} or ($ags{clen} and not $at{qlen})) ? 1 : 0; # bad for dropped qlen
    my $issplign= ($ags{gescore} or $SRCKEY =~ /spl/) ? 1 : 0; 
    if($issplign) { #was /gescore=/
      $tag="gspl" unless($tag);
      if($at{splice}) { $at{splice}= 2 + int($at{splice}/2); }
      #x$at{pid}||= 99; 
      my($pcov,$match,$ql2)= $at{cov} =~ m/(\d+).,(\d+).(\d+)/;
      if($match and $ql2){ $at{cov} = $pcov; $at{match}= $match; $at{qlen}= $ql2; }
      $at{qlen} ||= $ags{clen}; $gaps=$ags{gaps}||0; 
    } else { $tag="gmap" unless($tag);
      $aamap= $at{aamap}||$at{aalen}||0;
      $at{qlen} ||= $ags{clen};
      ($chi)= $id =~ m/_C(\d)$/; if($ags{chimera} or $chi) { $chi||=1; $at{Split}="$chi/2" unless($at{Split}); }
      if($at =~ /aalen=(\d+,\d+[^;\s]+)/) { $aaq=$1; $at{aalen}=$aaq; $at{aamap}=$aamap if($aamap ne $aaq); }
    } 
    if(m/cxlen=(\d+).(\d+)/) { my($cw,$xw)=($1,$2); $at{qlen}=$xw; } #$at{cdslen}=$cw;
    if(m/intron=([^;\s]+)/) { my $inv= $1; 
      # FIXME: diff syntax for evg17: intron=ivalid/itotal .. hard to parse here, do above
      if($inv =~ m/^(\d+),(\d+)$/) { my($mival,$mint)=($1,$2); 
        $at{splice}= 2 * $mint unless($at{splice});
        # no, count exons $nxf # $at{nexon}= 1+$mint unless($at{nexon} or $SKIPATTR);
      }
      elsif( my($nsp,$nxn)= $inv=~m,(\d+)/(\d+), ) {  #?? always intron=100%,42/42; nexon=22;
        $at{splice}=$nsp unless($at{splice});
        # no, count exons $nxf # $at{nexon}= int($nxn/2) unless($at{nexon} or $SKIPATTR);
        }
    }
    
    $at{tag}=$tag; ($ltid,$ltb,$lte)=($tid,$tb,$te);
    ($lgid,$lgb,$lge,$lgv,$lgo,$lat)=($gid,$gb,$ge,$gv,$go,$at);
  } 
}

#END 
$issplit=0; putg(); 
 
sub putg { 
  # was below 3 places..
  $at{nexon}= $nxf if(not $at{nexon}); # $nxf and 
  $at{ncds} = $ncds if(not $at{ncds}); # $ncds and .. splits should have done this above.. do all above?
  $at{match}= $xmatch if($xmatch and not $at{match});
  ## ^^ handling split for all added parts should be done above: match nexon ncds

  ## fixme: nSplit>2 .. splign up to 4 now
  if($issplit and $ltid eq $tid) { 
    #a $at{nexon}= $nxf if($nxf and not $at{nexon}); # fix?? losing split CDS, not exon?
    #a $at{ncds}= $ncds if($ncds and not $at{ncds});
    %lat=%at; $cov1=$at{cov}; 
    if($gid ne $lgid) { $lat{Split}="C3:$lgid:$lgb-$lge,$lgo,$cov1"; }
    elsif($lge < $gb-1000 or $lgb > $ge+1000) { $lat{Split}="C2:$lgb-$lge,$lgo,$cov1"; } 
    else { $lat{Split}="C1:$lgb-$lge,$lgo,$cov1"; } return; }
  
  elsif($at{Split} and $lat{Split}) { 
    # 201406 fix Split paths, for cov1 > cov2, swap; here..
    my $cov2=$at{cov}||0;  # FIXME cov=0 is valid /NOPATH
    my $cov1=$lat{cov}||0;
    if($cov2 and $cov1 and $cov2 < $cov1) {
      my($og,$ob,$oe,$oo,$nsp);
      my $sp=$lat{Split}; 
      if($sp=~/^C3:/) { 
        ($og,$ob,$oe,$oo)= $sp=~m/C3:([^:]+):(\d+).(\d+),(.)/; 
        $nsp="C3:$lgid:$lgb-$lge,$lgo,$cov2";
      } elsif($sp=~/^C[12]:/) { 
        $og=$lgid; ($ob,$oe,$oo)= $sp=~m/C[12]:(\d+).(\d+),(.)/; 
        $nsp=$sp; $nsp=~s/:.*/:$lgb-$lge,$lgo,$cov2/;
      }
      if($nsp and $oe) { $lat{Split}=$nsp; $lgid=$og; $lgb=$ob; $lge=$oe; $lgo=$oo;  } #?? ltb,lte
    }
    #a $at{nexon}= $nxf if($nxf and not $at{nexon});
    #a $at{ncds}= $ncds if($ncds and not $at{ncds});
    %att=%lat; map{ $att{$_}+=$at{$_}; } qw(cov match nexon ncds aamap);  # add ncds here
    $att{cov}=100 if($att{cov}>100); ## fixme: cov>100 not allowed.. gsplign split bug
    %at=%att; %lat=(); 
  } 
  if($at{Split}) { $at{path}= $at{Split}; }
  $at{cdspan}="$ocb-$oce" if($CDSPAN);
  if($ocb>=$mgb){
    if(not $at{offs}) { 
    if($lgo eq "-") { my $ocbr=1+$mge-$oce; $oce=1+$mge-$ocb; $ocb=$ocbr; }
    else { $ocb=1+$ocb-$mgb; $oce=1+$oce-$mgb;  }
    $at{offs}="$ocb-$oce"; # this is bad, oce end has intron spans, genome span, not cds-offs on transcript
    ## need more complex cdsloc - exonutr loc to recalc cds-off of transcript
    }
  }
  #a $at{nexon}= $nxf if(not $at{nexon}); # $nxf and 
  #a $at{ncds} = $ncds if(not $at{ncds}); # $ncds and .. splits should have done this above.. do all above?
  #a $at{match}= $xmatch if($xmatch and not $at{match});
  ## ^^ handling split for all added parts should be done above: match nexon ncds
  $at{nexon}.=".".$at{ncds}; # print only nexon col, .ncds as fraction
  $at{splice}= 2+$inspl if($inspl>0 and not $at{splice});
  if($at{nocov} and $at{match}>0 and $at{qlen}>0){ $at{cov}=int(0.5+100*$at{match}/$at{qlen}); }
  $at{aamap}= int($cdsw/3) if($cdsw>0 and not $at{aamap});
  $lpid=($OIDisID)?$ltid:$at{ID}||$ltid; 
  if($KEEP_G) { $lpid=~s/_C\d+$//; } 
  else { $lpid=~s/_[CG]\d+$//; } # 1708.FIXME: keep _G sometimes, separate row per _G multimap
  @at= @at{@hd}; print join("\t",$lgid,"$lgb-$lge",$lgo,$lpid,"$ltb-$lte",@at)."\n"; 
} 

## sort not required ..
# | sort -k1,1 -k2,2n -k4,4 > $nam.align.tab


