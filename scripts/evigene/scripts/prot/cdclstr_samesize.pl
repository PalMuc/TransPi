#!/usr/bin/env perl
# evigene/scripts/prot/cdclstr_samesize.pl : cd-hit .clstr parser, from cd-hit scripts

my $pMINLEN = $ENV{minlen} || 0.95;  # as prop 0.90 
   $pMINLEN=$pMINLEN/100 if($pMINLEN>1);
my $pMINID  = $ENV{minid}  || 0; # as pct 90% 
   $pMINID = 100*$pMINID if($pMINID > 0 and $pMINID < 1);
my $doTINY  = $ENV{tiny} || 0;
# my $doSAME  = $ENV{same} || 1;
## add no-same, self only rows ..
my $dosingle= $ENV{nosingle} ? 0 : 1;
my $doself= $ENV{self} ? 1 : 0;
my $noalt= $ENV{noalt} ? 1 : 0;
my $nosamegr= $ENV{groups} ? 1 : 0;
$noalt=1 if $nosamegr;

use constant DOSORT=>1;

my ($no,$nclin,$ncltr,$ngtr,$nitem,$nput,$topid,$topsize) = (0) x 9;
my @trow;

sub putclstr {
  my($topid,$topsize)= @_;
  $nclin++;
  my %gid;

if(DOSORT and @trow) {
   my @srow= sort{ my($aid,$aw)=split"\t",$a; my($bid,$bw)=split"\t",$b; $bw<=>$aw or $aid cmp $bid; } @trow;
   @trow=@srow;
}

  unshift(@trow, join("\t","self",$topsize,100)) if($doself);
  if(@trow) { $ncltr++; $ngtr+= @trow; }
  elsif($dosingle) { push(@trow, join("\t","self",$topsize,100)); }

  foreach my $tr (@trow) {
	my($tid,$tsize,$pid)=split"\t",$tr;
	my $pdif= $tsize/$topsize; $pdif= 1.0/$pdif if($pdif>1 and not $doTINY);
	my $putrow= ($doTINY) ? $pdif < $pMINLEN : $pdif >= $pMINLEN;
	$putrow=0 if($pid < $pMINID);
	if($noalt and $putrow) { my $gid=($tid eq "self")?$topid:$tid; 
	   $gid=~s/utrorf//; $gid=~ s/t\d+$//; 
           if($nosamegr) { $gid=~s/k\d\d[Ll]oc.*//; $gid=~s/[Ll]oc.*//; }
           $putrow=0 if($gid{$gid}++); }
	if($putrow) 
	{
		 $pdif= int(1000*$pdif)/10;  $nput++;
		 print join("\t",$nclin,$topid,$topsize,$tid,$tsize, $pid, $pdif)."\n";
	}
    }
}

sub putsum {
  #?? ncl2=$ncltr, , nitem2nd=$ngtr
  my $type=($doTINY)?"Tiny":"Samesize";
  print "#Summary: $type nout=$nput, ncluster=$nclin, in_items=$nitem\n";
}

=item clstr

>Cluster 0
0       15315nt, >zeamay:ncbig103630487t1... at +/100.00%
1       16278nt, >Zeamay5fEVm000001t1... *
2       7605nt, >Zm00001d038475_T006... at +/100.00%

=cut

print join("\t",qw(IClstr TopID TopLen QueryID QueryLen pIden pSize))."\n";
while(my $ll=<> ) {
  if ($ll =~ /^>/) {
  	putclstr($topid,$topsize) if($topid or @trow>0); @trow=();
   	$topid = ""; $topsize= $no = 0;
  } else {
    chomp($ll);
    my ($pid,$tid,$tsize) = (0) x 9;
    if ($ll =~ /(\d+)(?:aa|nt), >(.+)\.\.\./) {
       $tsize = $1; $tid = $2;
       ($pid)= ($ll =~ m/ at ([\d\.]+)/)?$1
              :($ll =~ m, at ./([\d\.]+),)?$1 : 0; # nt: has [+-]/ orient leader now
    } else {
      die "cd-hit clstr format error: $ll";
    }

    if ($ll =~ /\*$/) {
      $topsize = $tsize; $topid= $tid;
    } else {
    	$pid=~s/100.00/100/; $pid=~s/(\.\d)\d/$1/;
    	push(@trow, join("\t",$tid,$tsize,$pid));
	}
    $nitem++;
  }
}

putclstr($topid,$topsize) if($topid or @trow>0);
putsum();

