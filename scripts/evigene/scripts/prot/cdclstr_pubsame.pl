#!/usr/bin/env perl
# clstr_pubsame.pl

my $pMINLEN = $ENV{minlen} || 0.95;
my $pMINID  = $ENV{minid} || 0.95;
my $pSAMELEN  = $ENV{samelen} || 0; # or? 0.95;
my $pubid   = $ENV{pubid} || "Thecc";
my $maxlist = $ENV{max} || 10;
my $doclnum = $ENV{clnum} || 0;
my $noref2  = $ENV{noref2} || 0;
my $doself  = $ENV{self} ? 1 : 0; # add 1506, default noself
my $MARKBIGT = $ENV{mark} || 1; # always on

# iteratively tabulate clusters w/ pubid and matching trasm id, 
# then lower ident/len criteria for next go w/ reduced unclustered.aa

my $refid = "";
my $rsize=0;
my $no=0; 
my (@rrow,@trow);
my ($nclin, $nclboth, $nclpub, $ncltr, $ngpub, $ngtr)= (0) x 10;

sub putsum {
  $ncltr -= $nclboth; $nclpub -= $nclboth;
  print "#Summary: ncluster=$nclin, clboth=$nclboth, clpubonly=$nclpub, cltronly=$ncltr, ngenepub=$ngpub, ngenetr=$ngtr\n";
}

sub putclstr {
  ## sum stats
  $nclin++;
  if($refid){ $nclpub++; $ngpub+=@rrow; } # ngpub not both now
  if(@trow) { $ncltr++; $nclboth++ if($refid); $ngtr+= @trow; }
  return unless($refid and @trow);

  ## FIXME2: refid=1st pubid not longest .. pick from rrow..
  ##below# unshift @rrow, join("\t",$refid,$rsize,$rident);
  @rrow= sort { my($ad,$aw)=split"\t", $a; my($bd,$bw)=split"\t", $b; $bw<=>$aw or $a cmp $b; } @rrow;
  ($refid,$rsize,$rident)= split"\t",$rrow[0];

  ## sort @trow by max size: 
  @trow= sort { my($ad,$aw)=split"\t", $a; my($bd,$bw)=split"\t", $b; $bw<=>$aw or $a cmp $b; } @trow;
  shift(@rrow) unless($doself);
  @rrow=() if($noref2);
  my $i=0;
  foreach my $tr (@rrow,@trow) { 
     my($tid,$tsize,$pid)= split"\t",$tr;
     if($pid == 0) { $pid=$rident; } # 
     else { $pid=~s/100.00/100/; $pid=~s/(\.\d)\d/$1/; }
     my $pdif= $tsize/$rsize; $pdif= 1.0/$pdif if($pdif>1); $pdif= int(1000*$pdif)/10;
     my $mark="";
     my $msame= ($tsize >= $rsize || ($pSAMELEN && ($tsize>=$rsize*$pSAMELEN)));
     if($MARKBIGT and $msame) { 
       $mark= ($tsize > 1.05 * $rsize)? "t>r **" : ($tsize>$rsize) ? "t>r *" : "t=r";
         #dont need# ($tsize == $rsize || ($pSAMELEN && ($tsize>=$rsize*$pSAMELEN))) ? "t=r" : "";
         #was# ($tsize == $rsize) ? "t=r" : "";
     }
     $mark.="\t$nclin" if($doclnum);
     print join("\t",$refid,$rsize,$tid,$tsize,$pid, $pdif, $mark)."\n";
     last if(++$i >= $maxlist);
  }
}

while($ll=<>){
  if ($ll =~ /^>/) {
    putclstr() if($refid or @trow > 0); # always now# if($refid and @trow > 0); # but ignore only @rrow?
    @rrow= @trow=();
    $refid= "";
    $rident=$rsize= $no = 0;

  } else {
    chomp($ll);
    $pident= $tsize= $tid= 0;
    if ($ll =~ /(\d+)(?:aa|nt), >(.+)\.\.\./) {
       $tsize = $1; $tid = $2;
       ($pident)= $ll =~ /at ([\d\.]+)/; # 0 for \*
    } else {
      die "format error $ll";
    }

    if($tid =~ m/^$pubid/) {
      # $refid=$tid; $rsize=$tsize; $rident=$pident;
      # ^^ Oops, need to collect many possible refid/cluster
      push(@rrow, join("\t",$tid,$tsize,$pident)); # do for all refid + other pubid 
      ##if($refid) { push(@rrow, join("\t",$tid,$tsize,$pident)); }
      unless($refid) { $refid=$tid; $rsize=$tsize; $rident=$pident; } ## push @rrow also ?? yes
    } else {
      push(@trow, join("\t",$tid,$tsize,$pident));
    }
  }
}
putclstr() if($refid or @trow > 0);
putsum();
## add summary stats: nclust, nbothclust, nonly1(pub/notpub)

