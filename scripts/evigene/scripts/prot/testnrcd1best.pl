#!/usr/bin/env perl
# testnrcd1best.pl : evg tr2aacds failing here.

use strict;
use warnings;

## no gzip for test
my( $cdhitclstr, $cdsnrseq, $cdsinseq, $aaqual)= @ARGV;

die "input args: cdhitclstr, cdsnrseq, cdsinseq, aaqual \n"
 unless( $aaqual and  -s $cdhitclstr and -s $aaqual);

warn "# nofrag_reassignbest( $cdhitclstr, $cdsnrseq, $cdsinseq, $aaqual);\n";
my @ret= nofrag_reassignbest( $cdhitclstr, $cdsnrseq, $cdsinseq, $aaqual);
warn "# nofrag_reassignbest = @ret\n";

sub nofrag_reassignbest
{
  my( $cdhitclstr, $cdsnrseq, $cdsinseq, $aaqual)=@_;
  my (%aq, %better, %aqv ); 
  my ($nbetter,$nrec, $ntest, $nerr)=(0) x 9;
  my $flagfile= "$cdsnrseq.isbest";
  unless( -s $cdsnrseq and -s $aaqual and -s $cdhitclstr) {
    # loggit(1,"ERR: nofrag_reassignbest missing cdsnr:$cdsnrseq, $cdhitclstr or aaqual:$aaqual"); 
    return($nbetter,$nrec);
  }
  # return($nbetter,$nrec) if( -f $flagfile or $dryrun);
  
  # runcmd("touch $flagfile"); push @tmpfiles, $flagfile; #??
  open(F,$aaqual) or die; # loggit(1,"ERR: fail aaqual:$aaqual");  
  #bad# while(<F>) { my($id,$aw,$gp,$aq,$tw)=split; $aq{$id}="$aw,$aq" if($aq); } close(F);
  while(<F>) { my($id,$aw,$gp,$aq,$tw)=split; 
    my($aww,$ap,$ac)=split",",$aq;  $ap=~s/\%//; 
    my $acv=($aq=~/utr/ or $id=~/utr/)?1:2; $acv++ if($aq=~/complete/);  
    $aq{$id}="$aw,$ap,$acv,$ac"; $aqv{$id}= $aw * $ap * $acv;
  } close(F);

  use constant fPCDS_DIFF => 1; # what? should change any small diff?
  my(@fragids,$mainid,$mainlen);  # %fragids, > better
  open(F, $cdhitclstr); ## or return($nbetter,$nrec);
  while(<F>) { # cd-hit clstr format
  
    if(/^>/) { # new cluster, decide if last mainid is good, or replace by 2nd best ?
      $nrec++; 
      if($mainid and @fragids) {
        my $aqmain= $aq{$mainid}; my($awm,$apm,$acm,$aqm)=split",",$aqmain;
do{ $nerr++; warn "ERR: $mainid, aq=$aqmain\n"; } unless($aqm);
do{ warn "INFO: nr=$nrec, $mainid, aq=$aqmain\n"; } if($nrec < 4);;
        my $okmain= ($aqmain=~/utrbad|utrpoor|partial/ or $mainid =~ /utrorf/)?0:1;
        unless($okmain) { 
          my $awmMIN= 0.90 * $awm; # CONFIG
          @fragids= sort{ $aqv{$b} <=> $aqv{$a} or $a cmp $b } @fragids;  # bad aq <>
          ## 2nd longest/aagood replaces 1st longest/aapoor
          for my $fid (@fragids) {
            my $aqfrag= $aq{$fid}; my($awf,$apf,$acf,$aqf)=split",",$aqfrag;
            my $okfrag= ($aqfrag=~/utrbad|utrpoor|partial/ or $fid =~ /utrorf/)?0:1;
            my $pdif= $apf - $apm;
do{ $nerr++; warn "ERR: $fid, aq=$aqfrag\n"; } unless($aqf);
# do{ warn "INFO: nr=$nrec, $mainid/$aqmain vs $fid, aq=$aqfrag, pdif=$pdif, okfrag=$okfrag\n"; } if( $nbetter < 4);;
            if($okfrag and $awf >= $awmMIN and $pdif >= fPCDS_DIFF) { 
              delete $better{$mainid}; $better{$fid}= $mainid;
do{ warn "INFO: nr=$nrec, $mainid/$aqmain vs $fid, aq=$aqfrag, pdif=$pdif, okfrag=$okfrag\n"; } if( $nbetter < 4);;
              $nbetter++; last;
              }
            } 
          }
        } 
      @fragids=(); $mainid=0; 
die "ERRS = $nerr" if($nerr>9);
       
    } elsif(/^(\d+)/) { my $i=$1;
      m/(\d+)(nt|aa), >(.+)\.\.\. (.+)$/;  
      my($tlen,$typ,$tid,$pinfo)=($1,$2,$3,$4); 
      my $ismain=($pinfo =~ /\*/)?1:0;  # FIXME: drop /$/; new merge fnum at end of line now;
      unless($tid) {
      } elsif($ismain) {
        $mainid= $tid; $mainlen=$tlen; # $pi=100;  
        $better{$mainid}= $mainid; # always set
      } else {
        push @fragids, $tid; #save tlen in hash for sort?; need to collect all, may not have mainlen yet.
      }
    } 
  } close(F);
   
  if($nbetter>0) { # replace cds longest seq w/ better frag seq
    my($ok,$nok)= (0,0); ## check nok == nrec ??
    $ok= open(INF,$cdsinseq);
    $ok= open(OUTF,">$cdsnrseq.best") if($ok);
    if($ok) { $ok=0;
      while(<INF>) { if(/^>(\S+)/) { my $id=$1; $ok=$better{$id}; $nok++; } print OUTF $_ if($ok); } 
      close(OUTF); close(INF); 
      if($nok < $nrec) { # ERROR?? skip rename
        warn("ERR: nofrag_reassignbest fewer than input: reassign=$nok, orig=$nrec for $cdsnrseq.best"); 
        # push @tmpfiles, "$cdsnrseq.best";
      } else {
        my $cmd="mv $cdsnrseq $cdsnrseq.old; mv $cdsnrseq.best $cdsnrseq";
        warn("OK: $cmd\n");
        # system($cmd);  
        # push @tmpfiles, "$cdsnrseq.old";
      }
    }
  }
   
  return($nbetter,$nrec); # what?
}

