#!/usr/bin/env perl
# trformat.pl

=item about
 
  reformat tr assembler transcript fasta to common evigene format
  per assembler: velvet/oases, soap, trinity, idba, pacbio, others as added
  inputs: list of files transcript.fasta[.gz] 
  detect input format, else leave as is
  outputs: one transcripts.tr

=item see also

  rewrite of various scripts from evigene/scripts/rnaseq/  
  veltrmake.sh
  processtr.sh
  
=item author
  
  don gilbert, gilbertd near indiana edu, 2012
  part of EvidentialGene, evigene/scripts/

=cut

use strict;
use Getopt::Long;

use constant VERSION => '20170211'; # add idba, pacbio

my $MINTR = $ENV{MINTR} || 180; # what?
my $prefix = $ENV{prefix} || ""; # what? NEED some default
my $filenameIsPrefix=0;
my $SeqFormat= undef;
my ($debug,$format)= (0) x 9;  
my ($output,$logfile)= (undef,undef); my @input=();
my $TRSUFPATT='tr|cdna|fasta|fa|fna';
my $TRINPATH=$ENV{path}||0; # clip trinity path=... very long hdrs; make default?
my $IDCHAR_OK='\w:\.';
my $IDCLEAN=1; # for asis, default?

# ? add test-only option, output known type, unknown
my $optok= GetOptions( 
  "input=s", \@input,
  "output:s", \$output, 
  "logfile:s", \$logfile,
  "format=s", \$format, 
  "prefix|idprefix=s", \$prefix, 
  "MINTR=i", \$MINTR,  
  "debug!", \$debug,
  "cleanid!", \$IDCLEAN,
  "filenameIsPrefix!", \$filenameIsPrefix,
  "SeqFormat:n", \$SeqFormat,
  );

push @input, @ARGV;

die "usage: trformat.pl  -output all.tr  -input transcripts.fasta[.gz]  tr2.fasta  tr3.fasta
  OR gunzip -c xxx*.tr.gz | trformat.pl -input stdin > allxxx.tr 
  opts:  -prefix=name -format=velvet|idba|soapt|trinity|pacbio -MINTR=$MINTR -debug
" unless($optok and @input);

my $outh= *STDOUT;
if(not $output and defined $output) { # use input name 
  ( $output= $input[0] ) =~ s/\.($TRSUFPATT).*//; $output.=".tr";
}
if($output and $output!~/stdout|^-/) { open(OUT, ">$output") or die $output; $outh= *OUT; }

my $logh= undef;
sub loggit{ my $dowarn=shift; my $s= join(' ',@_); chomp($s);
  if($logh){ print $logh "#trf: $s\n"; } elsif($dowarn||$debug){ warn "#trf: $s\n"; }}

if(not $logfile and defined $logfile) { # use output name
  $logfile= $output||$input[0];
  $logfile=~ s/\.($TRSUFPATT).*//; $logfile.=".trformat.log";
}
if($logfile) { open(LOG, ">>$logfile") or die $logfile; $logh= *LOG; }

my $inform;
if($format =~ /vel/i) { $inform="velvet"; }
elsif($format =~ /trin/i) { $inform="trinity"; }
elsif($format =~ /soap/i) { $inform="soap"; } 
elsif($format =~ /idba/i) { $inform="idba"; } # 
elsif($format =~ /pacbio/i) { $inform="pacbio"; }
else { $inform=""; } # unknown

my $prefixin= $prefix || "evg";
my $kmer=""; # global ? for testformat...
my $ntotal=0;

my $LINELEN=($SeqFormat>9)?$SeqFormat:80; 
if(defined $SeqFormat) { $SeqFormat=$LINELEN; }

my %ids; my ($lastdupid,$ndupid)=(0,0);

=item dupids fixme

## FIXMEd3: 2013june21; 1 dupid left; fatal to tr2aacds !!!!!
## >> this is 1st id in file == dupid bug: kfish2wvelvk65Loc1t1

#t2ac: bestorf_cds= kfish2pt1all.cds nrec= 3351555
#t2ac: CMD= /home/ux455375/bio/evigene/scripts/prot/../prot/aaqual.sh kfish2pt1all.aa
#t2ac: ERR: dup id:kfish2wvelvk65Loc1t1 in kfish2pt1all.aa
#t2ac: FATAL ERR: 1 duplicate ids in kfish2pt1all.aa
#..
# [ux455375@trestles-login2 subsets]$ grep '^>' kfish2wvel_pt1.tr | sed 's/ .*//; s/>//;' | sort | uniq -c | grep -v ' 1 ' | head
#       2 kfish2wvelvk65Loc1t1
# ## dupid syntax: kfish2wvelvk65Loc71620d89568t1
# perl -pi.old -e \
# 'if(/^>kfish2wvelvk65Loc1t1/ and /len=2405/) { s/kfish2wvelvk65Loc1t1/kfish2wvelvk65Loc1d170000t1/; }' \
#   kfish2pt1all.{aa,cds,tr}

# /home/ux455375/scratchn/chrs/kfish/rnas/trsets/subsets
# [ux455375@trestles-login2 subsets]$ ls
# kfish2bveli_pt1.tr	      kfish2msop_pt1.trformat.log  kfish2wvel_pt1.tr
# kfish2bveli_pt1.trformat.log  kfish2mvel_pt1.tr		   kfish2wvel_pt1.trformat.log
# kfish2msop_pt1.tr	      kfish2mvel_pt1.trformat.log
# kfish2wvel_pt1.trformat.log
#trf: read format=velvet, in=vel2w/vel2w_65/transcripts.fa
#trf: nout=150768, nok=150769, ndupid=0
#trf: read format=velvet, in=vel3wn/vel3wn_65/transcripts.fa
#trf: WARNING changed  97266 duplicates to unique ids   << ** MISSED 1 kfish2wvelvk65Loc1t1 ****
# [ux455375@trestles-login2 subsets]$ grep ndupid= *mat.log | grep -v ndupid=0
# kfish2wvel_pt1.trformat.log:#trf: nout=138956, nok=138956, ndupid=97266
# kfish2wvel_pt1.trformat.log:#trf: nout=98856, nok=98856, ndupid=169574
#..............
## FIXMEd2: check for dupids over all @input; auto-reid as needed, with warning..
## works ok
#trf: read format=velvet, in=vel1_95/transcripts.fa   ## vel1_95 and vel3_95 giving same ids
#trf: nout=11782, nok=11782, ndupid=0
#trf: read format=velvet, in=vel3_95/transcripts.fa
#trf: WARNING changed  928 duplicates to unique ids
#trf: nout=982, nok=982, ndupid=928   
# >ztickb3velvk95Loc882t1 nt=1; cf=1.000; len=1635
# >ztickb3velvk95Loc882d925t1 nt=1; cf=1.000; len=411 ## undupid

=cut

foreach my $inf (@input) {
  my($inh,$nok,$nout,$nin)=(undef,0,0,0);
  ## loggit($debug, "read in=$inf\n") unless($inform);

  ## FIXME : velv, soap? pull kmer# from infile, add to prefix 
  ## in=velv5nun1/vel5nun1_75/transcripts.fa
  ## in=trsoap4set/trsoap4xco/sod29/sodmag4xco.k29.scafSeq.gz
  ## >dmag5Loc1t1 nt=6; cf=0.167; len=534 << missing prefix info from infile
  ## update:
  ## trformat.pl -log -out dmag5vel5nun1.tr -prefix dmag5nun1 -in velv5nun1/vel*/transcripts.fa
  ## >dmag5nun1velvk55Loc1t1 nt=6; cf=0.167; len=534

  ## 201607: k\d\d\d fix, longer reads, longer kmer
  $kmer="";
  if($inf =~ m,\b(k\d\d|k\d\d\d)\b,) { $kmer=$1; }
  elsif($inf =~ m,_(\d\d|\d\d\d)/(transcript|contig),) { $kmer='k'.$1; }
  
  $prefix= $prefixin;
  if($filenameIsPrefix) {  my $pn=$inf; $pn=~s,^.*/,,g;  $pn=~s/\W.*$//; $prefix=$pn; }
  # if($filenameIsPrefix) {  my $pn= `basename $inf`;  $pn=~s/\W.*$//; $prefix=$pn; }
  # $prefix .= $kmer if($kmer);
  
  if($inf =~ /stdin|^-/) {
    ($inform,$inh)= testformat(1,$inform,$inf);
  } elsif( -d $inf) { # handle dir of files?
    opendir(D,$inf) or loggit(1, "cant read indir=$inf\n");
    my @intr= map{ chomp; "$inf/$_"; } grep /\.($TRSUFPATT)/, readdir(D);
    closedir(D); 
    push @input, @intr;
  } elsif( -f $inf ) {
    ($inform,$inh)= testformat(0,$inform,$inf);
  } else {
    loggit(1, "cant read in=$inf\n");
  }

  if(defined $inh) {
    loggit(0, "read format=$inform, in=$inf\n");
    ($nok,$nout,$nin)= reformat($inform,$inh); $ntotal+= $nout;

    loggit(1, "WARNING changed ",$ndupid-$lastdupid,"duplicates to unique ids\n") if($ndupid>$lastdupid);
    loggit(($ndupid>0)?1:0, "nout=$nout, nok=$nok, ndupid=$ndupid\n");
    $lastdupid= $ndupid;
  }
}
loggit(0, "ntotal=$ntotal, idprefix=$prefix, output=$output\n");

#-----------------------------------------------------

sub testformat {
  my($isSTDIN, $inform, $inf)= @_;
  my ($isformat,$nrec,$inline,$outline);
  
  my($inh);
  if($isSTDIN) { $inh= *STDIN; }
  elsif($inf =~ /\.gz/) { open(IN,"gunzip -c $inf|") or warn"ERR: reading $inf"; $inh= *IN; }
  else { open(IN,$inf) or warn"ERR: reading $inf"; $inh= *IN; }
  
  my $prefixpre= $prefixin;
  if($filenameIsPrefix) {  my $pn=$inf; $pn=~s,^.*/,,g;  $pn=~s/\W.*$//; $prefixpre=$pn; }
  #if($filenameIsPrefix) {  my $pn= `basename $inf`;  $pn=~s/\W.*$//; $prefixpre=$pn; }

  if($inform) # presume known
  {
  $prefix= $prefixpre;
  $prefix .= substr($inform,0,4) if($inform); 
  $prefix .= $kmer if($kmer);
  return($inform,$inh);
  }
  
  # loggit(0, "read format=$inform, in=$inf\n");
  my $TESTPREFIX="PREFIX";  
  # while( ($inline = <$inh>) and $inline !~ /^>/) { } # test only 1st >header
  while( $inline = <$inh> ) { if($inline =~ /^>/) { last; } } # test only 1st >header
  loggit($debug, "testline1=$inline");
  
  ($isformat,$nrec,$outline)= trvelvet(1,$TESTPREFIX,$inline); # outh
  if($isformat) { $inform="velvet";} 
  else {
    ($isformat,$nrec,$outline)= trtrinity(1,$TESTPREFIX,$inline);
    if($isformat) { $inform="trinity";  } 
    else { 
      ($isformat,$nrec,$outline)= trsoap(1,$TESTPREFIX,$inline);
      if($isformat) { $inform="soap";  } 
      else {
        ($isformat,$nrec,$outline)= tridba(1,$TESTPREFIX,$inline);
        if($isformat) { $inform="idba";  } 
        else {
          ($isformat,$nrec,$outline)= trpacbio(1,$TESTPREFIX,$inline);
          if($isformat) { $inform="pacbio";  } 
          else { $outline= $inline; } ## $inform="asis"; ?
          }
        }
    }
  }
  
  $prefix= $prefixpre;
  $prefix .= substr($inform,0,4) if($inform); 
  $prefix .= $kmer if($kmer);
  if($outline) {
   	$outline =~ s/^>$TESTPREFIX/>$prefix/;
    if(/^>/) { 
    	my($id)= $outline =~ m/^>(\S+)/; use constant NOTEST=>0;
    	if(my $id1= undupid($id,$prefix,NOTEST) ) { $outline =~ s/>$id/>$id1/; } # FIX miss-1st-dup-id bug ; above testformat=1 prevents undupid
   	}	
    print $outh $outline; #?? here or defer? fix prefix? ** trsoap wants hdr defered to seq output
    $ntotal += 1;
  }
  return($inform,$inh);
}


sub reformat {
  my($inform,$inh)= @_;  # add 1st inline for troutput handler ..
  if( $inform eq "velvet" ) { return trvelvet(0, $prefix, $inh, $outh); }
  if( $inform eq "trinity" ) { return trtrinity(0, $prefix, $inh, $outh); }
  if( $inform eq "soap" ) { return trsoap(0, $prefix, $inh, $outh); }
  if( $inform eq "idba" ) { return tridba(0, $prefix, $inh, $outh); }
  if( $inform eq "pacbio" ) { return trpacbio(0, $prefix, $inh, $outh); }
  return trasis(0, $prefix, $inh, $outh);
  ## all:   return($isformat, $nput, $nrec);
}


sub formseq {
  my($s)=@_;
  #soap: $sq=~s/(.{60})/$1\n/g; print $HOUT $hd,$sq,"\n"; 
  if(length($s) > $LINELEN) { $s=~ s/(.{$LINELEN})/$1\n/g; }
  return $s;
}

sub trvelvet {
  my($testformat,$prefix,$inh,$outh)= @_;
  my ($isformat,$nrec,$nput,$ok)=(0) x 10; 
  # local $_= ($testformat) ? $inh : <$inh>;
  if($testformat) { $_= $inh; } else { $_= <$inh>; }
  
  ## >dmag5Loc1t1 nt=6; cf=0.167; len=534 << missing prefix info from infile
  ## add velvet/contigs.fa
  ## >NODE_1_length_1402_cov_195.409409   : ctg.fa
  ## >Locus_1_Transcript_1/1_Confidence_1.000_Length_1496   : tr.fa

  $ok=1; ## need after testformat ! skip no good first call
  do {
    if(/^>/) { 
      my $len=0;
      if(s,/(\d+)_Confidence_, nt=$1; cf=,) { $isformat++;
        s/_Length_(\d+)/; len=$1/; $len=$1; s/Locus_/Loc/; s/_Transcript_/t/; 
      } elsif(s/>NODE_(\d+)_length_(\d+)_cov_(\S+)/>Loc${1}ct1 /) {  $isformat++;
        $len=$2; my $cov=$3; $cov=~s/(\.\d)\d+/$1/; s/$/ cov=$cov; len=$len/;
      }
      s/>/>$prefix/;
      $ok= ($len >= $MINTR)?1:0;
      $nrec++;  $nput++ if($ok);
      my($id)= m/>(\S+)/; if($ok and my $id1= undupid($id,$prefix,$testformat) ) { s/>$id/>$id1/; }
      return($isformat,$nrec,$_) if($testformat); # return input line for STDIN ??
    } else {
      #  $_= formseq($_) if($SeqFormat);
    }
    print $outh $_ if($ok); # testformat got here .. bad
  } while(<$inh>);
  return($isformat, $nput, $nrec);
}


sub trtrinity {
  my($testformat,$prefix,$inh,$outh)= @_;
  my ($isformat,$nrec,$nput,$ok)=(0) x 10; 
  ## local $_= ($testformat) ? $inh : <$inh>;
  ## UPD 2015.06 for new variant of trin ids: TR45|c0_g1_i1 TR1069|c0_g1_i1 TR9549|c0_g2_i1 TR9553|c1_g1_i1
  ## dang.dang new id format 2014: >c0_g1_i1  ; _i1 > t1 is alt num ??
  ### double dang.dang TRINITY_DN39037_c0_g1_i1
  if($testformat) { $_= $inh; } else { $_= <$inh>; }
  do {
    if(/^>/) { $nrec++; $nput++; 
      #o# $isformat++ if( s/comp/loc/ and s/_seq/t/ ); 
      if( s/comp/loc/ and s/_seq/t/ ) { $isformat++; } 
      elsif( s/TR(\d+)\|c(\d+)_g(\d+)_i(\d+)/Loc$1c$2g$3t$4/ ) { $isformat++; } ## UPD 15 change
      elsif( s/TRINITY.(\w+).c(\d+).g(\d+).i(\d+)/Loc$1c$2g$3t$4/ ) { $isformat++; } ## UPD 15 change
      elsif( s/>c(\d+)_g(\d+)_i(\d+) />c$1g$2t$3 /) { $isformat++; } # 2014up
      unless($TRINPATH) { s/\s*path=.*$//; } ## path=[4779:0-175 11481:176-250 7363:251-384] 
      else { my($tp)=m/path=.([^\[\]\n]+)/; $tp=~s/\s/,/g; s/path=.([^\[\]\n]+)./$tp/; }
      s/_c/c/; s/>/>$prefix/; 
      my($id)= m/>(\S+)/;  if(my $id1= undupid($id,$prefix,$testformat) ) { s/>$id/>$id1/; }
      return($isformat,$nrec,$_) if($testformat); # return input line for STDIN ??
    } else {
      # $_= formseq($_) if($SeqFormat);
    }
    print $outh $_;
  }  while(<$inh>);

  return($isformat, $nput, $nrec);
}

=item idba_tran upd 1702
  
  see makeidbatr.sh
  ## each transcript-kmer.fa
  cat $subd/transcript-??.fa | env tag=${nam}rid perl -ne \
  'if(/^>/) { s/transcript-(\d+)_/$ENV{tag}k${1}Loc/;  $ok=1;  } print if($ok); '  
  
    # final contig.fa
  cat $subd/contig.fa | env tag=${nam}fid perl -ne \
  'if(/^>/) { s/contig-\d+_//; s/>_/>/; s/>/>$ENV{tag}k01Loc/; ($w,$rc)=m/(?:length|read_count)_(\d+)/g;
   s/(length|read_count)_(\d+)/$1=$2;/g; $ok=1; $maybeok=($w>=200 and $rc>0)?1:0; } print if($ok); '  

=cut

sub tridba {
  my($testformat,$prefix,$inh,$outh)= @_;
  my($isformat,$nput,$nrec,$maxl,$cldid)=(0) x 9; 
  if($testformat) { $_= $inh; } else { $_= <$inh>; }
  do {
    if(/^>/) { $nrec++; $nput++; 
      my($id);
      ($id)= m/>(\S+)/;
      if($id=~m/^transcript-(\d+)_/) { $isformat++;
        s/transcript-(\d+)_/idbtk${1}Loc/; 
      } elsif($id=~/^contig-\d+_/) { $isformat++;
        s/contig-\d+_//; s/>_/>/; s/>/>idbck01Loc/; 
        my($w,$rc)= m/(?:length|read_count)_(\d+)/g;
        s/(length|read_count)_(\d+)/$1=$2;/g;  
        # my $maybeok=($w>=200 and $rc>0)?1:0;
      }
      
      s/>/>$prefix/; 
      ($id)= m/>(\S+)/;  if(my $id1= undupid($id,$prefix,$testformat) ) { s/>$id/>$id1/; }
      return($isformat,$nrec,$_) if($testformat); # return input line for STDIN ??
    } else {
      # $_= formseq($_) if($SeqFormat);
    }
    print $outh $_;
  }  while(<$inh>);

  return($isformat, $nput, $nrec);
}

sub trasis { ## use instead trsoap() ?
  my($testformat,$prefix,$inh,$outh)= @_;
  my ($isformat,$nrec,$nput,$ok, $cldid)=(0) x 10; 
  while( <$inh> ) {
    if(/^>/) { $nrec++; $nput++; $isformat++; 
      if($IDCLEAN) { my($id)= m/>(\S+)/; 
        if($id=~/[^$IDCHAR_OK]/) { $id=~s/[^$IDCHAR_OK]/_/g; s/>\S+/>$id/; $cldid=1; } 
      }
      s/>/>$prefix/; 
      my($id)= m/>(\S+)/; if(my $id1= undupid($id,$prefix,$testformat) ) { s/>$id/>$id1/; }
      return($isformat,$nrec,$_) if($testformat); # return input line for STDIN ??
    } else {
      # $_= formseq($_) if($SeqFormat);
    }
    print $outh $_;
  }
  return($isformat, $nput, $nrec);
}


=item pacbio smrt isoseq  2017.02 upd
  
  messy ids 
  smrt _CSS and classifier ouputs have raw read ID like this, change/cut rawid
  >m150514_204533_42207_c100759272550000001823162807221566_s1_p0/15/1361_56_CCS strand=-;fiveseen=1;polyAseen=1;threeseen=1;fiveend=31;polyAend=1336;threeend=1363;primer=1;chimera=0

  isoseq after clustering: "cNNN" is ID, /nn/nn are info addons, not needed
  but "cNNN" is relative to single input file, as per other trasm sets
  >c5/3/476 isoform=c5;full_length_coverage=3;isoform_length=476

=cut

sub trpacbio {
  my($testformat,$prefix,$inh,$outh)= @_;
  my($isformat,$nput,$nrec,$maxl,$cldid)=(0) x 9; 
  if($testformat) { $_= $inh; } else { $_= <$inh>; }
  do {
    if(/^>/) { $nrec++; $nput++; 
      my($id);
      ($id)= m/>(\S+)/;
      if($id=~m,/, and $id=~m/_CCS/) { $isformat++; }
      elsif($id=~m,/, and $id=~m,^c\d+/\d+/,) { $isformat++; } # >cNNN/ccc/lll
      if($IDCLEAN and $isformat) {
        if(length($id) > 20 and $prefix) { $id=~s,_\w+/,,; } 
        if($id=~/[^$IDCHAR_OK]/) { $id=~s/[^$IDCHAR_OK]/_/g; s/>\S+/>$id/; $cldid=1; } 
      }
      s/>/>$prefix/; 
      ($id)= m/>(\S+)/;  if(my $id1= undupid($id,$prefix,$testformat) ) { s/>$id/>$id1/; }
      return($isformat,$nrec,$_) if($testformat); # return input line for STDIN ??
    } else {
      $_= formseq($_) if($SeqFormat);
    }
    print $outh $_;
  }  while(<$inh>);

  return($isformat, $nput, $nrec);
}

sub trsoap {
  my($testformat,$prefix,$inh,$outh)= @_;
  my($isformat,$nput,$nrec,$maxl)=(0) x 9; 
  my($hd,$sq)=("","");
  
  our $HOUT=$outh; ## UNDEF sym: outh 
  sub putrsoap{ my($hd,$sq)=@_; my $sl= length($sq); our $HOUT;
    if($sl>=$MINTR){ $hd=~s/$/ len=$sl;/ if($hd); # hd may be missing
      $sq=~s/(.{$LINELEN})/$1\n/g; print $HOUT $hd,$sq,"\n"; return 1; } else { return 0; }
  }
  
## BAD rec1 output missing seq..
# ==> dmag4xco1soap.tr <==
# >dmag4xcosoapk21loc0t1 i=1; val=26.7; flag=COMPLEX;
# >dmag4xcosoapk21loc2t1 i=2; val=50.0; flag=COMPLEX; len=15465;
# TTTTTTTTTTTTTTTACATGTAAAGTATTTATTTTTTAAGTTAACTTTCAGGTTGAAGAT

  ## local $_= ($testformat) ? $inh : <$inh>;
  if($testformat) { $_= $inh; } else { $_= <$inh>; }
  do {
    if(/^>/) { 
      $nput+= putrsoap($hd,$sq) if($hd or $sq);  ## FIXME here bad rec1, >hd1 printed, then calls back w/ just seq input..
      $hd=$sq="";  $nrec++; 
      my($i,$l,$t,$val,$fl,$ok)=(0) x 9;
      if(m,>scaffold(\d+)\s+Locus_(\d+)_(\d+)\s+(\S+)\s*(\S*),) { 
        ($i,$l,$t,$val,$fl)=($1,$2,$3,$4,$5);  $ok=1;
        $maxl=$l if($l>$maxl); ++$t; $fl="; flag=$fl" if($fl); 
      } elsif(m/^>C(\d+)\s*(\S*)/) { 
        ($i,$val)=($1,$2);   $ok=1;
        $l=++$maxl; $t=1; $fl=""; 
      } # elsif(not $testformat) { loggit(1, "trsoap.funny.format: $_"); } # not soap format ?
      # old1712: ${prefix}loc${l}t${t} >> change to Loc${l} for consistency
      if($ok) { $isformat++; s/>.*$/>${prefix}Loc${l}t${t} i=$i; val=$val$fl;/; }
      else { s/>/>$prefix/; } # what? flag bad format?
      my($id)= m/>(\S+)/; if(my $id1= undupid($id,$prefix,$testformat) ) { s/>$id/>$id1/; }
      $hd= $_; 
      return($isformat,$nrec,$_) if($testformat); # return input line for STDIN ??
    } else { 
      chomp; $sq .= $_;
    } 
  }  while(<$inh>);

  $nput+= putrsoap($hd,$sq) if($hd or $sq); 
  return($isformat, $nput, $nrec);
}

sub undupid {
  my($id,$prefix,$testformat)= @_;
  if($testformat) { return 0; } ## BAD miss-1st-id-dup bug; But also need for current test-many-formats; 
  	## FIX in testformat: need 2nd call here after testformat goes to output 1st line == >id-dup
  elsif($ids{$id}++) { 
    my $id1;
    do { 
      $id1= $id; $ndupid++; 
       # change where? prefix? append after tNN will cause altid confusion
      unless( $id1 =~ s/(t\d+)$/d$ndupid$1/) { $id1 =~ s/$prefix/${prefix}d$ndupid/; }
    } while( $ids{$id1} );
    $ids{$id1}++;
    #cant do in sub# s/>$id/>$id1/; 
    return $id1;
  } else { return 0; }
}

