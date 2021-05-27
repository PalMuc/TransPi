#!/usr/bin/env perl
# blast92orthomcl10.pl

=item update 10
  updated 2010: 
    clean up code, remove old,unused eugenes blast92hgsum parts
    
    make -orthomcl flag default (any -noortho uses?)
    
    add -bitscoreeval option to recalc e-values from bitscore
  
=item update 11

  bug: doesnt add partial matches of same gene; matters for big prots
  
=item usage

  gtar -Ozxf $po/prot22_modDM.blout.tgz | \
  ../blast92orthomcl.pl -deflines=$pp/prot22.h -in=stdin -out=prot22mclDM -orthomcl -org=modDM,modSC,modCE,ensAG

=cut

use strict;
use warnings; #?

use Getopt::Long;    
use POSIX;

my $VERSION = "2011.11.07";  

# FIXME: need more ID prefix > Org mappings:
#  sp|XXX_SPECIES = uniprot w/ species tail
#  hxXXXX = special org gene set, hx -> myspecies
my $hasorgprefix= 1; # default?  daphnia_ID format gene ids
my $dbpipeid= 1; # FIXME; what is blast/fasta id syntax?
my $faprefix= 'gnl\|'; # other patts 

use constant SUMSCORE => 1; # 2011.11
my $OVSLOP=6; # FIXME: need overlapslop ~ < 0.01 of max part span ; or < 0.02..0.05 of min part span?
my $pctOVSLOP= 0.05; # 2014.07, was 0.02; # use -pctover 0.05 opt
## 201408: new filter opt, this is bad fixed bitscore filter, due to split aligns vs fullalign
## $bit_score > 0.25 * $bmax or $bspans{$sid}
my $BITSMAXFILT= 0.05; # was 0.25


my %IDprefixmap=();
my $mapids="";
# %IDprefixmap=( "dpx19" => '^hx', 'human' => '^sp:' ); # need as option

my ( $infile, $outfile, $append, @blastcmd,  $orthomcl, );
my ( $allorgs, @wantorg, %wantorg,  @skiporgs, %skips); 
my ( $deflines, $mineval, $minbits, $minident,  $minalign, $simidStart,  ) = (0) x 10;
my %skiporg= ();
my $debug= 1;
my $reciprocalfake= 0;
my $checkreciprocal=0; my %reciprocal; # idhash; may get TOO big, all*all
my $ggonly=0;

## $mineval= 1e-3;
$outfile="omcldefault";
$orthomcl= 1;
my $calcevalbits=0;
my($qorg,$sorg)=('','');


my $optok= Getopt::Long::GetOptions( 
'infile=s' => \$infile, # @infile ?
'outfile=s' => \$outfile,
'append!' => \$append,
'ggonly!' => \$ggonly,
'debug!' => \$debug,
'identmin=n' => \$minident,
##'alignmin=n' => \$minalign, # NOT ADD 2013oct, minident is doing this
'evalmin=s' => \$mineval,
'bitsmin=n' => \$minbits,
'pctover=n' => \$pctOVSLOP,
'calc_eval_from_bitscore=n' => \$calcevalbits,
'deflines|fasta=s' => \$deflines,
'simidStart=s' => \$simidStart,
'organism=s' => \@wantorg,
'skiporgs=s' => \@skiporgs, # drop? == non-eugenes org
'mapids=s' => \$mapids,
'reciprocalfake!' => \$reciprocalfake,
'checkreciprocal!' => \$checkreciprocal,
);
#'orthomcl!' => \$orthomcl, # make default?


$infile= shift(@ARGV) unless($infile);

die "USAGE: $0
 -in=infile|stdin -out=outfile|stdout -fasta=fasta-deflines
Options:  
  -append                append output
  -identmin=$minident    [min %ident to keep]
  -bitsmin=$minbits      [min bitscore to keep]
  -evalmin=$mineval      [min e-value to keep]
  -calc_eval_from_bitscore=database_seqlen_size  
  -organism=orgtag of source 
  -skiporg=orgtag to skip 
  -simidstart=$simidStart [start sim. id]
  -mapids='daphnia19=^hx,human=^sp:'
  -checkreciprocal  check for reciprocal blast values
  -reciprocalfake  fake reciprocal blast values
"
unless($optok && $infile);
#  -orthomcl        [write in:protblast-type9 to orthomcl .bpo,.gg] 

$pctOVSLOP=$pctOVSLOP/100 if($pctOVSLOP > 0.9);

if($mapids) { 
# %IDprefixmap=( "daphnia19" => '^hx', 'human' => '^sp:' ); # need as option
# -mapids='daphnia19=^hx,human=^sp:'
  my @ids= split /[,\s]+/, $mapids;
  map{ my($k,$v)=split"="; $IDprefixmap{$k}=$v; } @ids;
}

# fixme: bug input quotes on wantorg ??
if (@wantorg) { @wantorg= map{ s/^\W+//; s/\W+$//; split(/[\s,]+/,$_); } @wantorg;  }
else { $allorgs=1; } 
$allorgs=1 if grep{$_ eq "all"} @wantorg;

if (@skiporgs) { map { @skiporg{ split(/[\s\,]/,$_) } = (1)x99; } @skiporgs; }

my $protsize={}; 
my($ggh, $bpoh, $outpath);
MAIN: {
  
  my $nid=0;
  if (-d $deflines) {
    opendir(D,$deflines); 
    my @df= grep(/^\w/, readdir(D)); closedir(D); # scan for .aa/.fa ?
    foreach my $in1 (@df) {
      my $nid1=0;
      ($protsize,$nid1)= readDeflines("$deflines/$in1", $protsize);
      $nid+= $nid1;
      }
  } else {
   ($protsize,$nid)= readDeflines($deflines, $protsize);
  }

  # set db size
  $calcevalbits= $nid 
    if(($calcevalbits>0 and $calcevalbits < 99) or (SUMSCORE && $calcevalbits == 0)); 


  my @didorgs;
  @wantorg= grep{ !$skiporg{$_} } @wantorg;
  %wantorg= map { $_ => 1; } @wantorg;

  $outpath= "$outfile.gg"; # blast genes groups
  if (!$append && -e $outpath) { rename ($outpath,"$outpath.old"); }
  if ($outfile && open(OUT,">$outpath")) { $ggh= *OUT; }
  else { $ggh= *STDOUT;  warn "bad outfile: $outpath"; }
  printOrthGG($ggh,$protsize); 
  close $ggh if ($outfile);
  exit(0) if($ggonly); # hack fix failed .gg why? bad .fa name 

  $outpath= "$outfile.bpo"; # blast prot? out
  if (!$append && -e $outpath) { rename ($outpath,"$outpath.old"); }
  if ($outfile && $append && open(OUT,">>$outpath")) { $bpoh= *OUT; }
  elsif ($outfile && open(OUT,">$outpath")) { $bpoh= *OUT; }
  else { $bpoh= *STDOUT;  warn "bad outfile: $outpath"; }

  if (-d $infile) {
    opendir(D,$infile); 
    my @df= grep(/^\w/, readdir(D)); closedir(D); # scan for .blastp ?
    foreach my $in1 (sort @df) {
      my @orgset= parseblastout9("$infile/$in1", $protsize, $bpoh); 
      push(@didorgs, @orgset); # get dups in @didorgs this way
      }
  } else {
    @didorgs= parseblastout9($infile,$protsize, $bpoh);  
  }
  
  close $bpoh if ($outfile);
}

if(scalar(%skips)) {  
  foreach my $t (sort keys %skips) { my $n= $skips{$t}; warn "#WARN: skipped:\t$t\t$n\n"; }
}

# this is likely time/memory sink; use only when needed? there are alot also, valid non-recip .. or blastp missed?
if($checkreciprocal) { 
  my @ids= sort keys %reciprocal; my $ni= @ids;
  for(my $i=0; $i<$ni; $i++) {
    my $qid=$ids[$i]; 
    my @sid= sort keys %{$reciprocal{$qid}};  my $nj= @sid; # small subset @sid/qid of all @ids
    for(my $j=0; $j<$nj; $j++) {
      my $sid=$sid[$j];
      next if($sid eq $qid);
      my $nqs=$reciprocal{$qid}{$sid}||0;
      my $nsq=$reciprocal{$sid}{$qid}||0;
      next if($nqs == $nsq and $nsq < 2);
      my $have= 0;
      $have |= 1 if $nqs; # ? check count?
      $have |= 2 if $nsq;
      if($have == 2){ warn "#ERR: reciprocalmiss:\t$qid\t$sid\n" ; }
      elsif($have == 1) { warn "#ERR: reciprocalmiss:\t$sid\t$qid\n" ; }
      if($nqs>1) { warn "#WARN: reciprocaldup:\t$qid\t$sid\t$nqs\n" ; }
      if($nsq>1) { warn "#WARN: reciprocaldup:\t$sid\t$qid\t$nsq\n" ; }
      }
    }
}



#------------------------

sub _min_not0 { return ($_[0] < 1 || ($_[1] > 0 && $_[1] < $_[0])) ? $_[1] : $_[0]; }
sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }
sub bint { local $_=shift; return (/e\+/) ? int($_) : $_; }

sub printOrthGG {
  my($ggh,$protsize)= @_;
  
  foreach my $db (sort keys %{$protsize}) {
    next unless($allorgs or $wantorg{$db});
    my $ng= 0;
    foreach my $id (sort keys %{$protsize->{$db}} ) {
      print $ggh "$db: " if (($ng % 50) == 0);
      print $ggh "$id ";
      $ng++;
      print $ggh "\n" if (($ng % 50) == 0);
      }
    print $ggh "\n";
  }
}  

# read either fasta deflines with size= or full fasta and count sizes
sub readDeflines {
  my($infile, $protsize)= @_;
  my ($inh,$lid,$asize, $aa);
  my $nid=0;
  
	#my %protsize=();
  $protsize={} unless(ref($protsize));
  
	warn "\nread readProts $infile for orthomcl\n" if $debug;
  
  my $ok=1;
  if (!$infile || $infile =~ m/^(stdin|\-)$/i) { $inh= *STDIN; $ok=1; }
  elsif ($infile =~ /\.(gz|Z)$/){ $ok= open(F,"gunzip -c $infile |");  $inh=*F; }
	else { $ok= open(F,$infile); $inh= *F; }
	return () unless($ok);
	
  while(<$inh>) {
    chomp;
    if (/^>(\S+)/){
      my $id=$1;  $nid++;
      if($lid && $asize) {
        if($aa =~ /\w/){ $aa=~s/\*$//; $asize=length($aa); }
        my($db,$id1,$idmore)= cleanSplitId($lid);
        $protsize->{$db}{$id1}= $asize;
      }
      if(m/(size|length)=(\d+)/) { # trust this? input may be only deflines
        $asize= $2;
        my($db,$id1,$idmore)= cleanSplitId($id);
        $protsize->{$db}{$id1}= $asize;
      }
      
      $aa=""; $asize= 0; $lid= $id;
      }
    else{ $aa.=$_; $asize += length($_); } # drop final '*' for consistency, if present
    }
  close($inh);
  if($lid && $asize) {
    if($aa =~ /\w/){ $aa=~s/\*$//; $asize=length($aa); }
    my($db,$id1,$idmore)= cleanSplitId($lid);
    $protsize->{$db}{$id1}= $asize;
  }
  return ($protsize, $nid);
}





sub cleanid { 
  unless( $_[0] =~ s/$faprefix//) { $_[0] =~ s/^gi\|\d+\|(\S+)/$1/; } 
  $_[0] =~ s/[\|,;]$//; #s/\|$//; ## ncbEC, ncbAt have '|' at end of id; chomp here?

  $_[0] =~ s/\|/:/ if($dbpipeid); #?? change to db:id ? NOT ALWAYS

  #? $_[0] =~ s/\|/ /; # change any 3rd pipe to space?
  ## problem now with modRR defline - query != subject in blout for same ID (99% paralogs!)
  ## blast query id has all |||| parts; blast subject (-oT) has only gnl|db|xxx 3 parts
  ## trunc here? only care about for ident compare: ($sid eq $qid)
}

sub cleanSplitId {
  my($id)= @_;
  cleanid($id);
  my($db,$did,$more)=('',$id,'');
  if(%IDprefixmap) {
    foreach my $dbpre (sort keys %IDprefixmap) {
      my $pre=$IDprefixmap{$dbpre};
      if($id =~ /$pre/) { 
        $db=$dbpre; $did=$id; $did=~s/:/_/;
        if($hasorgprefix and $did !~ /^$db/) { $did= $db."_".$did; } # keep it in id also :((
        last; 
        }
    }
  }
  unless($db) {
  my $dsplit= ($hasorgprefix) ? '[_:]' : '[:]';
  ($db,$did,$more)= ($id =~ /$dsplit/) ? split(/$dsplit/, $id, 2) : ('NULL',$id);
  if($hasorgprefix) { $did= $db."_".$did; } # keep it in id also :((
  }
  # if($db eq 'NULL' && $id =~ /_/ && $hasorgprefix) { ($db,$did)= split(/[_]/,$id,2); }
  ($did,$more)= ($did =~ /[\|]/) ? split(/[\|]/, $did, 2) : ($did,'');
  return ($db,$did,$more);
}

sub dbOfId {
  my($id)= @_;
  my $dsplit= ($hasorgprefix) ? '[_:]' : '[:]';
  my($db,$did,$more)= ($id =~ /$dsplit/) ? split(/$dsplit/, $id, 2) : ('NULL',$id);
  return $db;
}

# rev of this? calcbitsOfeval ?  2**bits = seqlen*dblen / eval ; bits= log(seqlen*dblen / eval)
sub calcevalbits {
  my($bits, $m_seqlen, $n_dblen)= @_;
  # bitscore to e-value :  e = m_seqlength * n_dblength / (2^S_bits)
  return 1 if($bits < 0 or $m_seqlen <= 0 or $n_dblen <= 0);
  my $eval= sprintf "%.0e", $m_seqlen * $n_dblen / (2**$bits);  # note 2**0 == 1
  # precision?  .2e or .0e ?
  return $eval;
}

## 2011.aug BUG here, need to test sb-se outside tb-te spans also
sub sumscore {
  my( $bspans, $q, $t, $pctident, $alen, $eval, $bits, $qb,$qe,$sb,$se) = @_;
  my $nident= $pctident * $alen;

  my $or=0;
  if($qb > $qe) { ($qb,$qe)= ($qe,$qb); $or--; }
  if($sb > $se) { ($sb,$se)= ($se,$sb); $or--; }
  unless($bspans->{$t}) { 
    $bspans->{$t}=[]; 
    push( @{$bspans->{$t}}, [$qb,$qe,$sb,$se,$bits, $alen, $nident, $eval]); 
    return; }
  my $ov=0;
  ## 2011oct overslop pct fix
  my $qlen=1+$qe-$qb; my $slen=1+$se-$sb;
  my $qslop= _max($OVSLOP, int($pctOVSLOP*$qlen));
  my $sslop= _max($OVSLOP, int($pctOVSLOP*$slen));
  
  foreach my $sp (@{$bspans->{$t}}) {
    my($xb,$xe,$tb,$te,$xbit)= @$sp;
    if($qe < $xb or $qb > $xe) { }
    elsif($qe > $xe and $qb >= $xe - $qslop) { }
    elsif($qb < $xb and $qe <= $xb + $qslop) { }
    else { $ov=1; last; }
  ## add 2011.aug
    if($se < $tb or $sb > $te) { }
    elsif($se > $te and $sb >= $te - $sslop) { }
    elsif($sb < $tb and $se <= $tb + $sslop) { }
    else { $ov=1; last; }
  }  
  push( @{$bspans->{$t}}, [$qb,$qe,$sb,$se,$bits, $alen, $nident, $eval]) unless($ov);
}



sub parseblastout9 {
	my($infile, $protsize, $bpoh)= @_;
	local(*R);
	my ($inh);
	my %didorgs=();
	
	warn "\nparsing $infile to orthoMCL .bpo, .gg\n" if $debug;
  if (!$infile || $infile =~ m/^(stdin|\-)$/i) { $inh= *STDIN; }
	else {
	  my $ok;
	  if ($infile =~ /\.(gz|Z)$/){ $ok= open(R,"gunzip -c $infile |"); }
	  else { $ok= open(R,$infile); }
	  unless ($ok) { warn "cannot read $infile"; return undef; }
    $inh= *R;
    }
  
  my $s_org= 'undef';
    
  my($l_qid,$l_sid,$l_qdb,$l_sdb,$l_prob);
  my($simid,$simspan,$simspanid,$sumident,$sumalign);
	my($nl,$query, $idmore, $sdb, $qdb, $lsh_id, $nqid, $bmax, $npart)=(0) x 20;;
	my @qlocs;
	my %bspans=();
	
	$simid= $simidStart;
	
  while(<$inh>) {
  
		if (/^#/) {
      #? drop all comment processing, expect only format 8 table?
			# parse some of this; this only comes w/ blout format 9, not 8 ..
		  if ($nl < 5 && $_ !~ /^# (Query|Fields)/) { chomp(); push(@blastcmd,$_); }
			if (/^# (\w*BLAST\w*)/) { 
			  $query= ''; 
			  }
			elsif (/^# (Database)/) { 
        if (m,/?([\w\.]+)\s*$,) {
	        $s_org= $1; $didorgs{$s_org}++;
	        warn "Subject: $s_org \n" if $didorgs{$s_org}==1 && $debug;
	        }
			  }
			elsif (/^# Query:\s*(\S+)(.*)/) { 
			  my $qid= $1;
			  $query= $1.$2; $nqid++; $npart= 0; @qlocs=(); 
        ($qdb,$qid,$idmore)= cleanSplitId($qid);
        $didorgs{$qdb}++;
        if($didorgs{$qdb}==1 && $debug) {
          my $keep= ($allorgs or  $wantorg{$qdb})?1:0;
          warn "Query: $qdb keep=$keep\n" ;
        }
			  }
			next;
	  }
    next unless (/^\w/); ## skip blanks, comments...
    
    my ($qid0, $sid0, $pctident, $alignment_length, $mismatches, $gap_openings, 
       $q_start, $q_end, $s_start, $s_end, $prob, $bit_score ) = split "\t"; #@v;

    my($qid,$sid);
    ($qdb,$qid,$idmore)= cleanSplitId($qid0);
    ($sdb,$sid,$idmore)= cleanSplitId($sid0);
    
    next unless($allorgs or ($wantorg{$sdb} && $wantorg{$qdb}));
    next unless($protsize->{$qdb}{$qid} and $protsize->{$sdb}{$sid}); ## ** CANT have zero-len data, skip on

    $bit_score= bint($bit_score);

# 2011.11: add loop here, collect all lines w/ same qid, assume input sorted by query

if( SUMSCORE ) {
    if($l_qid and $qid ne $l_qid) {
      processQuery( $l_qid, \%bspans, );
    
      %bspans=(); $bmax=0;
    }

    # ?? always use calcevalbits here if score is sum of parts.
    sumscore( \%bspans, $qid, $sid, $pctident, $alignment_length, $prob, $bit_score, $q_start, $q_end, $s_start, $s_end) 
      if($bit_score > $BITSMAXFILT * $bmax or $bspans{$sid}); # should drop this filter? 201408: YES, bad w/o sumbitscore, for hsp parts vs whole
      ## should use this BITSMAXFILT only after sumscore() in processQuery() ; not same as minbits..
      
    #old#$bmax= $bit_score if($bit_score > $bmax);
    $bmax= $bit_score if($bit_score > $bmax and $qid0 ne $sid0); # 12jan; dont let self match stop add
    $l_qid = $qid; $l_sid = $sid; $l_qdb= $qdb; $l_sdb= $sdb;

} else {
  # OLD here....................
  
    if($calcevalbits>0) {
	    # bitscore to e-value for orthomcl:  e = m_seqlength * n_dblength / (2^S_bits)
      my $qlen= $protsize->{$qdb}{$qid} || 0; ## or what?
      my $slen= $protsize->{$sdb}{$sid} || 0; ## or what?
      my $proborig= $prob;
      
      $prob= calcevalbits($bit_score, _min_not0($qlen,$slen), $calcevalbits) ;
    }
    
#    ## not here, may be part of gene
#     next if ($minident && $pctident < $minident);
#     next if ($minbits && $bit_score < $minbits);
#     next if ($mineval && $prob > $mineval);

    ++$npart;
    my $nident= $pctident * $alignment_length;
    
    ## fixme 2011.11; check via sumscore() for overlap hits
    if ($qid eq $l_qid && $sid eq $l_sid && $sdb eq $l_sdb) {
      # collect hsp for 1 line output
      $simspanid++;
      $sumident += $nident;
      $sumalign += $alignment_length;
		  $simspan.="$simspanid:$q_start-$q_end:$s_start-$s_end.";
      next;
      }
    
    printBPO($bpoh,$protsize,$simid,$l_qdb,$l_qid,$l_sdb,$l_sid,$sumident,$sumalign,$l_prob,$simspan);
    printBPO($bpoh,$protsize,++$simid,$l_sdb,$l_sid,$l_qdb,$l_qid,$sumident,$sumalign,$l_prob,revspan($simspan))
     if($reciprocalfake);
    
    $l_qid = $l_sid = $l_qdb= $l_sdb= 0;
    next if ($minident && $pctident < $minident); 
      # FIXME: add pctalign filter: $pctident * $alignment_length / min(trlen1,trlen2) to remove tiny aligns
    next if ($minbits && $bit_score < $minbits);
    next if ($mineval && $prob > $mineval);

	  $simid++;
    $simspanid= 1;
    $sumident = $nident;
    $sumalign = $alignment_length;
		$simspan="$simspanid:$q_start-$q_end:$s_start-$s_end.";
    $l_qid = $qid; $l_sid = $sid; $l_qdb= $qdb; $l_sdb= $sdb;
    $l_prob= $prob;
} # OLD

  }


if( SUMSCORE ) {
  processQuery( $l_qid, \%bspans, );
  $simid= $simidStart;
  
} else {  
  printBPO($bpoh,$protsize,$simid,$l_qdb,$l_qid,$l_sdb,$l_sid,$sumident,$sumalign,$l_prob,$simspan);
  printBPO($bpoh,$protsize,++$simid,$l_sdb,$l_sid,$l_qdb,$l_qid,$sumident,$sumalign,$l_prob,revspan($simspan))
     if($reciprocalfake);
}
	warn "last simid=$simid\n" if $debug;

  return sort keys %didorgs;
}


sub processQuery {
  my($qid, $bspans)= @_;

  my $simid= $simidStart;
  my $qdb= dbOfId($qid);
  my($maxb,$maxt)= (0,0);
  
  ## FIXME sid order should be sorted by max bits (or nident, or alen)
  my (%scored,%bpo);
  my $qlen= $protsize->{$qdb}{$qid} or return; ## || 0; ## or what? ** CANT have zero-len data, skip on
  my $docheckre= ($checkreciprocal and !$reciprocalfake) ? 1 : 0;
  
  foreach my $sid (sort keys %$bspans) {
    my @sp= @{$bspans->{$sid}}; 
    my $sdb= dbOfId($sid);
    my $np= @sp;
    my $simspan="";
    my ($simspanid, $sumbits, $sumident, $sumalign, $eval1) = (0) x 10;
    foreach my $sp (@sp) {
      my($qb,$qe,$sb,$se,$bits, $alen, $nident, $eval)= @$sp;
      # [$qb,$qe,$sb,$se,$bits, $alen, $nident] # add blast prob
      $eval1= $eval unless($simspanid);
      $sumbits  += $bits;
      $sumident += $nident;
      $sumalign += $alen;
      $simspanid++;
      $simspan .= "$simspanid:$qb-$qe:$sb-$se.";
    }

    my $slen= $protsize->{$sdb}{$sid} or next; ## || 0; ## CANT have zero-len data
    my $prob= ($np<=1) ? $eval1 : calcevalbits($sumbits, _min_not0($qlen,$slen), $calcevalbits) ;
    #^ fixme, use blastp prob when 1 sp part
    my $percentIdent= ($sumalign>0)? int($sumident/$sumalign) : 0;  

    # my $checkreciprocal=0; my %reciprocal; # idhash; may get TOO big, all*all
    unless( $docheckre and $reciprocal{$sid}{$qid} and !$reciprocal{$qid}{$sid}) {    
    ## dont need all these.. use minident > minalign?
    # next if ($minident && $percentIdent < $minident);
    ## count skips?
    if($minbits  && $sumbits < $minbits) { $skips{minbits}++; next; }
    if($mineval  && $prob > $mineval){ $skips{mineval}++; next; }

## 201408: YES, bad w/o sumbitscore, for hsp parts vs whole
## .. move here from above? need maxbits per query tho..
##      if($bit_score > $BITSMAXFILT * $bmax or $bspans{$sid}); # should drop this filter? 
## should use this BITSMAXFILT only after sumscore() in processQuery() ; not same as minbits..

    ## cant have qlen==0 && slen==0; also in bpo, 0-len bad; forbid nolength data? yes, above
    my $maxlen=_max($qlen,$slen); # _max(max,1);
    if($minident &&  (100*$sumalign/$maxlen) < $minident){ $skips{minident}++; next; } # filterbpo.pl WEAKMAT => 0.40;
        # ^^ this same as minalign filter?
        
    if($checkreciprocal && $reciprocal{$qid}{$sid} ) { $skips{reciprocaldup}++; $reciprocal{$qid}{$sid}++; next; } # have already; count/warn?
    }
    
    $bpo{$sid}= "$qid;$qlen;$sid;$slen;$prob;$percentIdent;$simspan";
    $reciprocal{$qid}{$sid}++ if($checkreciprocal);
    $scored{$sid}= $sumbits;
  }
  
  foreach my $sid (sort { $scored{$b} <=> $scored{$a} or $a cmp $b } keys %scored) {
    my $bpo= $bpo{$sid};
    $simid++;
    print $bpoh "$simid;$bpo\n";
  }
  
  if($reciprocalfake) {
    foreach my $sid (sort { $scored{$b} <=> $scored{$a} or $a cmp $b } keys %scored) {
      my $bpo= $bpo{$sid};
      my @bpo= split";",$bpo;
      $bpo[-1]= revspan($bpo[-1]);
      @bpo[0,1, 2,3]= @bpo[2,3, 0,1];
      $bpo= join";",@bpo;
      $simid++;
      #never here# $reciprocal{$sid}{$qid}++ if($checkreciprocal);
      print $bpoh "$simid;$bpo\n";
    }
  }
  
  $simidStart= $simid;
  #return($maxb, $maxt);
}


sub revspan {
  my($span)= @_;
  my $revspan="";
  my @hsp=split /[\.]/, $span;
  foreach my $h (@hsp) {
    my($sid,$qs,$ss)=split":",$h;
    $revspan.="$sid:$ss:$qs.";
  }
  return $revspan;
}

##? should warn/not-write if cant find protsize qid or sid ?
#1169636;1302933|CAE51908|RGD|Mcpt1l3|mast;249;1302933;0;2e-135;100;1:1-249:1-249.
#1169637;1302933|CAE51908|RGD|Mcpt1l3|mast;249;1303300;0;4e-125;94;1:1-249:1-247.
## need some db-tags on odd ids '1303300' is what? ncbi-ginum? -- no, RGD-id
## >gnl|modRR|1303300|CAE48391|RGD|Mcpt1l4|mast CRC64=AFDD2657A6573FF1; size=247;

sub printBPO {  # OLD unused
  my($bpoh,$protsize,$simid,$qdb,$qid,$sdb,$sid,$sumident,$sumalign,
    $prob,$simspan)= @_;
  return unless($simspan && $sumalign);
  # fixup to use blast align size where protsize == 0
  if($qid eq $sid and $qdb eq $sdb and $protsize->{$qdb}{$qid} == 0) {
	  $protsize->{$qdb}{$qid}= $sumalign;
	}
  my $qlen= $protsize->{$qdb}{$qid} || 0; ## or what?
  my $slen= $protsize->{$sdb}{$sid} || 0; ## or what?
  my $percentIdent=int($sumident/$sumalign); # not single pctident 

  print $bpoh "$simid;$qid;$qlen;$sid;$slen;$prob;$percentIdent;$simspan\n";
}




1;

__END__

=item old notes

# sort -t ';' -k2,2 -k7,7nr prot4mcl.bpo > prot4mcl.s
^^ sort by qid,pctident not good, need sort by -k6,6 = e.val (1e-10) = sort wont do reals

# cat prot*.bpo | sort -t ';' -k2,2 -k7,7nr |\
# perl -pe '$i++; s/^(\d+)/$i/;' > $op/prot4mcl.bpo 
# ...

cat prot*.bpo |\
perl -ne'@v=split";"; ($m,$p)=split"e-",$v[5]; $m=10-$m; $p=99999 if($v[5] == 0.0);print "$p.$m;$_";' |\
sort -t ';' -k3,3 -k1,1nr | perl -pe '$i++; s/^[\d\.]+;\d+/$i/;' > $op/prot4mcl.bpo 

./orthomcl.pl --mode 4 --usr_bpo_file=$op/prot4mcl.bpo --usr_gg_file=$op/prot4mcl.gg

FIXME: orthomcl.pl assumes .bpo/blast input sorted by query; our input is
  separated by source db; also need sort with best match 1st in query set.
  
------ all prot22 - modDP -------
set po=$bg/shared/out/prots1

foreach blout ($po/prot22_*.blout.tgz)
if ( $blout =~ *modDP* ) then
echo skip  modDP
elif ( $blout =~ *sanPF* ) then
echo skip  modDP
else
gtar -Ozxf $blout | \
../blast92orthomcl.pl -deflines=$pp/prot22.h -in=stdin -out=prot22mcl \
   -append -orthomcl -skip=modDP  
endif
end

## sanPF blastout is useless - got bogus 'Pfa3D7' gene id for all
cat prot22mcl.bpo | grep -v Pfa3D7 |\
perl -ne'@v=split";"; ($m,$p)=split"e-",$v[5]; $m=10-$m; $p=99999 if($v[5] == 0.0);print "$p.$m;$_";' |\
sort -t ';' -k3,3 -k1,1nr | perl -pe '$i++; s/^[\d\.]+;\d+/$i/;' > $op/prot22mcl.bpo 

See blast92hgsum.pl

=cut
