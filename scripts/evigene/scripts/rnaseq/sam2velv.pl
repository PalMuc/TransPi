#!/usr/bin/env perl
# sam2velv.pl

=item about

  extract seq from sam into fasta (paired/unpaired) for velvet,

=item velvet call : add here? gen script?

  velbin=/bio/bio-grid/mb/rnaseq/velvet/bin
  # kmer=17 uses more reads than kmer=21 ... shorter is better here?
  # kmer=15 uses more reads than kmer=17
  
  $velbin/velveth vel 17 -fasta -shortPaired *.fa2 -short *.fa1 *.fa  -long *.longfa 
  $velbin/velvetg vel -read_trkg yes 
  $velbin/oases vel -ins_length 200 
    # alt opts:
    oases -min_pair_count 2 -cov_cutoff 2 -ins_length 150 -ins_length_sd 50
  
=cut


use strict;
use warnings;
use Getopt::Long ;
use File::Basename;

use constant PAIR1DUP => 1; # test for NO_SEQDUP for pairs; drop if 1 of pair done

my $do_sort2= 1; # debug??
my $WARN_NOPATH=0;
my $MAX_BADSEQ_NNN= 0.33;
my $MAX_MISMATCH = 0.05;  # 0.05 may be too high for assembly, 0.03?
my $MAX_ENDCLIP  = 0.33;   # separate from mismatch err
my $NO_SEQDUP= 0;
my $NO_IDCHECK=0;
my $NDUP = 1; # == number dupl seqs to keep for NO_SEQDUP, ~coverage; 
   # should this correspond to velvet cov_cutoff=3, min_pair_count=4


my( $sreadsuffix, $pairtmpsuffix, $pairsuffix, $unpairsuffix, $estsuffix)
    = (".fa",".tmpfa2",".fa2",".fa1", ".longfa");

my(  $intype, $inname, $format, $samoutput, $ispaired, $debug);

$format="fasta"; # both output? change to fa2,fa,sam == outsuffix
$intype="sam"; # or bam, only for .fa1
# $inname="reads"; #??
my $onlysplit=0;
my $hardclip=0;

my $optok= &GetOptions (
  'name=s', \$inname, # for stdin, outname
  'type=s', \$intype, # EST/longfa | paired | single ...
  'format=s', \$format, # fasta,fastq,.. add samoutput here?
  "samoutput:s"=>\$samoutput,  
  "noduplicateseq:+"=>\$NO_SEQDUP,  
  "noidcheck!"=>\$NO_IDCHECK, # bad option??
  "mismatchmax=s"=>\$MAX_MISMATCH, # is this bp or portion?
  "endclipmax=s"=>\$MAX_ENDCLIP, # is this bp or portion?
  "onlysplit!"=>\$onlysplit,
  "hardclip!"=>\$hardclip,
  "debug!"=>\$debug,  
  );

die 
"usage: sam2velv -name reads in.sam [in2.sam .. inN.sam]
 options: 
    -noduplicateseq -noidcheck 
    -samoutput:name
    -mismatchmax=$MAX_MISMATCH (%length)
    -endclipmax=$MAX_ENDCLIP (%length)  -hardclip 
    -type EST|paired|singleread -debug\n" 
unless($optok and $inname);

$NDUP=$NO_SEQDUP if($NO_SEQDUP>1);

$MAX_MISMATCH= $MAX_MISMATCH/100 if($MAX_MISMATCH>1);
$MAX_ENDCLIP= $MAX_ENDCLIP/100 if($MAX_ENDCLIP>1);
$MAX_MISMATCH=1.0 if ($MAX_MISMATCH <= 0.0);
$MAX_ENDCLIP=1.0 if ($MAX_ENDCLIP <= 0.0);

if($format =~ /fastq/) { # is this ok?
map{ s/fa/fq/ } ($sreadsuffix, $pairtmpsuffix, $pairsuffix, $unpairsuffix, $estsuffix);
}

my $singlesuf= ($intype =~ /EST|long/i) ? $estsuffix : $sreadsuffix;

if($onlysplit) {
  my $pairtmp= shift(@ARGV);
  if(-f $pairtmp) { split_fasta( $pairtmp, $inname) ; }
  else { warn "missing pairtmp input file \n"; }
} else {
  sam2fasta( $inname, $singlesuf); # # <> not right here... 
}

#............


# sam2fasta simpler than _partitioner, read .sam per part, write .fa, .tmpfa2 // .fa2, .fa1
sub sam2fasta {
  my($outnam, $singlesuf)= @_;  # $insamh,
  
  my ($nin,$npaired,$nunpaired, $ndup, $nbad)=(0) x 10; 

  $outnam=~ s/\.(sam|bam)$//;
  my $tmppairfile  = "$outnam$pairtmpsuffix";  # ".tmpfa2";
  my $unpairfile   = "$outnam$singlesuf"; #".fa";
  my $samfile      = "$outnam.sam"; 
  
  ## FIXME: need both out pair, unpair here
  open(OP,">$tmppairfile") or die "writing $tmppairfile\n";   
  open(OU,">$unpairfile") or die "writing $unpairfile\n"; 
  
  my $dosam=0;
  if(defined $samoutput) {
    if($samoutput =~ /\w/) { $samfile=$samoutput; }
    if(-f $samfile) { $samfile.=".velv"; }
    open(OSAM,">$samfile") or die "writing $samfile\n"; 
    $dosam=1;
  }
  
  my (%didid, %didseq);
  while(<>) { 
    next unless(/^\w/); # if(/^\@/ or /^\W/); # comments, sam.headers
    chomp;
    my($dm,$flag,$chr,$cloc,$qq,$cig,$mc,$mb,$pd,$seq,$qual,@opts)=split"\t";  $nin++; 
    my($softclip, $nmismat, $nclip)=(0,0,0);
    my $sl= length($seq) || 1;
    
    if(/NM:i:(\d+)/) { $nmismat=$1; } #?? got 2 bams w/o this; warn if no NM:i ?
    elsif($MAX_MISMATCH <= 0.05) { $nbad++; next; } # skip; warn?? per lib problem
    if($nmismat and $nmismat/$sl > $MAX_MISMATCH) { $nbad++; next; }
    
    ## change this: separate NM test and SH test;
    ## change: make MAX_ENDCLIP nclip/slen portional
    while($cig =~ m/(\d+)([SH])/g) { $nclip += $1; $softclip=1 if($2 eq 'S'); } 
    if( $nclip and $nclip/$sl > $MAX_ENDCLIP ) { $nbad++; next; }
   
    my $sn= $seq=~tr/Nn/Nn/; 
    if($sn and $sn/$sl > $MAX_BADSEQ_NNN or $seq=~/null/) { $nbad++; next; }

    if($softclip and $hardclip) { ($cig, $seq, $qual)= hardclip($cig,$seq,$qual); }
    my $loc= ($chr eq "*")? " loc=$chr" : " loc=$chr:$cloc";
    
    if($flag & 0x01) { # paired
      my $m=($flag & 0x80)?2:1; 
      my $dmate="$dm/$m"; 
      my $seqo= $seq;
      if($m == 1) { 
        $seqo=reverse($seqo); $seqo=~tr/ACGTacgt/TGCAtgca/; 
        $qual=reverse($qual);  
        }   
        
unless($NO_IDCHECK) { 
  ## MUST?? unless filter out dups in caller
  ## without this get messed up .fa2:  idx/1, idx/1, idx/1, idx/2, .. idy/1, idy/2, idy/2...
      next if($didid{$dmate}++); # drop multiples
}
      $npaired++;
      if($format =~ /fastq/) {
      print OP join "\t", '@'.$dmate,$seqo,'+'.$dmate,$qual."\n";
      } else { 
      print OP ">$dmate$loc\t$seqo\n"; 
      }
    } else { # single
      if($NO_SEQDUP) { do { $ndup++; next; } if($NDUP <= $didseq{$seq}++); }
      ##NO, not for singles ## else { unless($NO_IDCHECK){ next if($didid{$dm}++); }} # drop multiples

      $nunpaired++; #? write as >id\tseq and later sort, change \t ??
      if($format =~ /fastq/) {
      print OU join "\n", '@'.$dm,$seq,'+'.$dm,$qual."\n";
      } else { 
      print OU ">$dm$loc\n$seq\n";  
      }
      
    }
     
    if($dosam) {
      my $outline= join("\t", $dm,$flag,$chr,$cloc,$qq,$cig,$mc,$mb,$pd,$seq,$qual,@opts);
      print OSAM $outline,"\n"; 
    }
  } 
      
  close(OU);  close(OP);
  if($dosam) { close( OSAM); }
  warn "# sam2fasta $outnam  nin=$nin, paired=$npaired, single=$nunpaired, ndupskip=$ndup, npoor=$nbad\n" if $debug;    
  %didid=(); %didseq=();
  unlink($unpairfile) unless($nunpaired>0); #??
  
  my $nsplit=0;
  $nsplit= split_fasta( $tmppairfile, $outnam) if($npaired>0);
  unlink($tmppairfile) if($nsplit>0 or $npaired==0); # unless($debug>1); #??
  
}


sub hardclip {
  my($cigar, $seq, $qual)= @_;
  if($cigar =~ m/\dS/) { 
    my ($oldcig,$cl,$cr)= ($cigar,0,0);
    if($cigar =~ s/(\d+)S$/$1H/) { $cr=$1; }
    if($cigar =~ s/^(\d+)S/$1H/) { $cl=$1; }
    if($cr or $cl) {
      if($seq and $seq ne "*") { 
       $seq=substr($seq,0,-$cr) if $cr;  
       $seq=substr($seq,$cl) if $cl; 
       }
      if($qual and $qual ne "*") { 
       $qual=substr($qual,0,-$cr) if $cr; 
       $qual=substr($qual,$cl) if $cl; 
       }
     }
   }
  return($cigar, $seq, $qual);
}


sub split_fasta {
  my($infile,$outnam)= @_;
    
  unless( -f $infile ) { warn "Missing $infile\n"; return; }
  
  $outnam=~ s/\.(fq\w+|fa\w+)$//; 
  my $pairfile  = "$outnam$pairsuffix";  # ".pairs";
  my $unpairfile= "$outnam$unpairsuffix"; #".unpair";
  open(OP,">$pairfile") or die "writing $pairfile\n";   
  open(OU,">$unpairfile") or die "writing $unpairfile\n"; 

  my ($nin,$npair,$nunpair, $lv, $ld, $lpair, $ndup)= (0) x 10;
  my %didseq;
  my $insort="$infile.sort";
  warn "sort $infile > $insort\n" if $debug;
  my $ok=system("sort $infile > $insort");
  ##open(IN, "sort $infile |"); 
  ## ^ is this memory pig?? ; try sort to infile.sort, read infile.sort ?
  open(IN, $insort);
  while(<IN>){  
    # assumes format is all lines of: ">id<tab>seq<endline>"; warn of others?
    unless(/^[>\@]/ and /\t/) { warn "#split_fasta: input not '>id.tab.seq' \n" if $debug; next; }
    my($d)= m,[>\@]([^/\s]+),; 
    s/\t/\n/g;  $nin++;
    
    my $skip=0;
    my $pair=($d eq $ld)?1:0;
    if($pair) { 
      if($NO_SEQDUP) {
        my($lsq)= $lv =~ m/\n(\S+)/;
        my($nsq)=  $_ =~ m/\n(\S+)/;
if(PAIR1DUP) {
        $skip=1 if($NDUP <= $didseq{$lsq}++ or $NDUP <= $didseq{$nsq}++); # should this test $nsq.$lsq also?
} else {
        $skip=1 if($NDUP <= $didseq{$lsq.$nsq}++); # should this test $nsq.$lsq also?
}        
      }
      unless($skip) { print OP $lv,$_; $npair++; } else { $ndup++; } 
      
    } elsif($lv and not $lpair) { 
      if($NO_SEQDUP) { my($lsq)= $lv =~ m/\n(\S+)/; $skip=1 if($NDUP <= $didseq{$lsq}++); }
      unless($skip) { print OU $lv; $nunpair++; } else { $ndup++; }  
    }
    ($lv,$ld,$lpair)=($_,$d,$pair);
  }
  if($lv and not $lpair) { print OU $lv;  $nunpair++; }
  close(OP); close(OU); close(IN);

  unlink($insort) unless($npair == 0); ## ($debug); #??
  #? unlink pairfile or unpairfile if empty ?
  unlink($pairfile) unless($npair>0); #  unless($debug); #??
  unlink($unpairfile) unless($nunpair>0); #  unless($debug); #??
  
  warn "# split $outnam ; nin=$nin, npair=$npair, nunpair=$nunpair; ndupskip=$ndup\n" if $debug;    
  return ($npair+$nunpair); #ok
}


__END__

=item velvet call : add here? gen script?

velbin=/bio/bio-grid/mb/rnaseq/velvet/bin

# echo "#.. start rnavelv : `date`"
$velbin/velveth vel 21 -fasta -shortPaired *.fa2 -short *.fa1 *.fa  -long *.longfa 
$velbin/velvetg vel -read_trkg yes 
$velbin/oases vel -ins_length 200 
# echo "#.. end rnavelv : `date`"

# Final graph has 560 nodes and n50 of 71, max 327, total 20323, using 5375/8572 reads
# ... try to get more reads used...  velveth vel2 17 ??

@ vel  21: oases Finished heuristic approach, used 6237/8572 reads
# vel2 17: oases Finished heuristic approach, used 6454/8572 reads
# vel4 15: oases Finished heuristic approach, used 6744/8572 reads
# ... these are all pretty puky results, vs cufflinks or apparent read data
# .. not long, use shorter kmer .. 17 or 15
# vel3 29: oases Finished heuristic approach, used 4390/8572 reads
# ..
# vel5 17: Finished heuristic approach, used 6602/8572 reads
#  + oases vel5 -ins_length 150 -min_pair_count 2 -cov_cutoff 2
#   > n=72 trans.fa   vs 60? for others
#  ^^ somewhat better exon joining here...

# vel6 17: Finished heuristic approach, used 6604/8572 reads; ntr=76
#  oases vel6 -min_pair_count 1 -cov_cutoff 2 -degree_cutoff 2
# .. Not as good as vel5

# vel6b: oases vel6 -min_pair_count 1 -cov_cutoff 6,  used 6092/8572 reads
# .. many fewer tra.fa;

# vel6c: oases vel6 -min_pair_count 1 -cov_cutoff 2 -degree_cutoff 6;  69 mRNA, used 6604/8572 reads

# vel6d: oases vel6 -min_pair_count 1 -cov_cutoff 1 -degree_cutoff 6;  94 mRNA, used 6690/8572 reads
# vel6e: oases vel6 -min_pair_count 1 -cov_cutoff 1 -degree_cutoff 1, 89 mRNA, used 5071/8572 reads
# .. these last two have poor exon joining, not good.

# vel6f: oases vel6 -min_pair_count 2 -cov_cutoff 2 -paired_cutoff 0.05, 69 mRNA, used 6425/8572 reads
# vel6g: oases vel6 -min_pair_count 2 -cov_cutoff 2 -paired_cutoff 0.3, 69 mRNA, used 6425/8572 reads

# vel6h: oases vel6 -min_pair_count 2 -cov_cutoff 2 -ins_length 150 -ins_length_sd 50, 62 mRNA, used 6602/8572 reads


# save vel/ transcripts.fa stats.txt splicing_events.txt contig-ordering.txt Log
# drop vel/ Graph2 LastGraph PreGraph Roadmaps  Sequences
# Graph2     Log       Roadmaps   contig-ordering.txt  splicing_events.txt  transcripts.fa
# LastGraph  PreGraph  Sequences  contigs.fa           stats.txt


=item rRNA read selection / assembly

  for abundant rRNA,
  -- pick only perfect reads, longest paired
  -- if not too many, keep all dupls for improved velvet asm
  -- else sam2velv -nodup=4 or =5
  -- key to good velvet assembly is using only perfectly aligned reads?
  
 test cacao: got 8 transcripts, 18S, 23S rRNA fragments + intergene spacer, ...
    -- size of 300bp to 1,000 bp (18S = ~1000, 23S = ~3000)
    -- missing end reads from selection methods, should expand rrna-gff-filter to get more?

# cacao{01,09}.rrna.sam ...
cat cacao*.rrna.sam | grep 'NM:i:0' | grep 'M   =' | egrep -v '[0-9](H|I|D)' | $rnas/rgsoft/sam2velv.pl -debug -nodup=4 -name cacao_perfrrna
cat cacao*.rrna.sam | grep 'NM:i:0' | grep 'M   =' | egrep -v '[0-9](H|I|D)' | $rnas/rgsoft/sam2velv.pl -debug -name cacao_perf2rrna

 test aphid:
 ll -h {aphidrs_SRR064408,aphidrs_SRR064409,aphidpe_SRR071347,aphidpe_SRR075802,aphidpe_SRR075803}.rrna.sam.gz
  11G Jan 25 04:38 aphidpe_SRR071347.rrna.sam.gz
 634M Jan 25 05:17 aphidpe_SRR075802.rrna.sam.gz
 547M Jan 25 05:52 aphidpe_SRR075803.rrna.sam.gz
  89M Jan 23 05:13 aphidrs_SRR064408.rrna.sam.gz
 624M Jan 22 23:13 aphidrs_SRR064409.rrna.sam.gz
 
 gzcat aphidpe_SRR071347.rrna.sam.gz | grep 'NM:i:0' | grep 'M =' | egrep -v '[0-9](I|D|H)' | \
   $rnas/rgsoft/sam2velv.pl -debug -nodup=4 -name aphidpe_perfrrna1

 gzcat aphidrs_SRR06440?.rrna.sam.gz | grep 'NM:i:0' | grep 'M  *' | egrep -v '[0-9](I|D|H)' |\
   $rnas/rgsoft/sam2velv.pl -debug -nodup=4 -name aphidrs_perfrrna3 > & log.perfrrna3 &

==> log.perfrrna1 <==
# sam2fasta aphidpe_perfrrna1  nin=77986660, paired=1997897, single=0; ndupskip=0
# split aphidpe_perfrrna1 ; nin=1997897, npair=5577, nunpair=38466; ndupskip=1438111

==> log.perfrrna2 <==
# sam2fasta aphidpe_perfrrna2  nin=31973172, paired=687604, single=0; ndupskip=0
# split aphidpe_perfrrna2 ; nin=687604, npair=7546, nunpair=31047; ndupskip=447472

==> log.perfrrna3 <==
# sam2fasta aphidrs_perfrrna3  nin=18174601, paired=0, single=34837; ndupskip=328618
# split aphidrs_perfrrna3 ; nin=0, npair=0, nunpair=0; ndupskip=0


grep -c '^>' aphid??_perfrrna?.{fa,fa1,fa2,tmpfa2}
aphidrs_perfrrna3.fa:34837
aphidpe_perfrrna1.fa1:38466
aphidpe_perfrrna2.fa1:31047
aphidpe_perfrrna1.fa2:11154
aphidpe_perfrrna2.fa2:15092

$velbin/velveth velrrna 21 -fasta -shortPaired aphidpe_perfrrna?.fa2 -short aphidpe_perfrrna?.fa1 aphidrs_perfrrna?.fa > & log.velh5 &
$velbin/velvetg velrrna -read_trkg yes > & log.velg5 &
$velbin/oases velrrna -ins_length 150 -ins_length_sd 50 -min_trans_lgth 150 > & log.velo5 &

#..........
# clean set; much larger, how much mem to do all at once?

grep -c '^>' aphid??_perfclean?.*fa*
aphidpe_perfclean1.fa1:4177998
aphidpe_perfclean1.fa2:6539088
aphidpe_perfclean2.fa1:6794101
aphidpe_perfclean2.fa2:11088520
aphidrs_perfclean3.fa:13713722

$velbin/velveth velclean 21 -fasta -shortPaired aphidpe_perfclean?.fa2 -short aphidpe_perfclean?.fa1 aphidrs_perfclean?.fa > & log.velh7 &
$velbin/velvetg velclean -read_trkg yes > & log.velg7 &
$velbin/oases   velclean -ins_length 150 -ins_length_sd 50 -min_trans_lgth 40 > & log.velo7 &

# clean set too big for velvet + 236 GB mem; 
# subset by large scaffolds? iterate over scaffold groups? by csize, readcount?

# Csize > 1,999,999
Chr             CSize   Nread
Scaffold10      3073041 2259762 0
Scaffold1       2622152 1379300 0
Scaffold2       2384549 1232937 0
Scaffold7       2848108 1049056 0
Scaffold21      2482828 921282  0
Scaffold19      2144968 876799  0
Scaffold16      2289837 800832  0
Scaffold17      2239900 795334  0
Scaffold23      2106074 723920  0
Scaffold18      2560806 691223  0
Scaffold4       2478080 605805  0
Scaffold9       2216851 425411  0
...

egrep 
'(Scaffold10|Scaffold1|Scaffold2|Scaffold7|Scaffold21|Scaffold19|Scaffold16|Scaffold17|Scaffold23|Scaffold18|Scaffold4|Scaffold9):'

cat midsub.list aphid*.fa2 | env suf=fa2 perl -ne'if(/^(mid\w+)\t(\S+)/){($n,$v)=($1,$2); push @n,$n; m
ap{ $sc{$n}{$_}=1; }split",",$v; } elsif(/^>/){ ($s)=m/(Scaffold\w+)/; $p=0; foreach $n (@n) { if($sc{$n}{$s}) {
$p=$n; last; } } if($p) { $fp=$fps{$p}; unless($fp){ my $fh;  open($fh,">subset.$p.$suf"); $fps{$p}=$fp= $fh; } p
rint $fp $_; } } elsif($p) { print $fp $_; } BEGIN{$suf=$ENV{suf}; }'


## add loc=Scaffoldxxx:a-b to sam2velv.fa and then split those ?
## -nodup=4  ??

gzcat aphidpe_SRR071347.clean.sam.gz | grep 'NM:i:0' | grep 'M =' | egrep -v '[0-9](I|D|H)' | \
 $rnas/rgsoft/sam2velv.pl -debug -nodup=4 -name aphidpe_perfclean1

gzcat aphidpe_SRR07580?.clean.sam.gz | grep 'NM:i:0' | grep 'M =' | egrep -v '[0-9](I|D|H)' | \
 $rnas/rgsoft/sam2velv.pl -debug -nodup=4 -name aphidpe_perfclean2

gzcat aphidrs_SRR06440*.clean.sam.gz | grep 'NM:i:0' | grep 'M *' |  egrep -v '[0-9](I|D|H)' | \
 $rnas/rgsoft/sam2velv.pl -debug -nodup=4 -name aphidrs_perfclean3

foreach cset ($csets)

 # set clist=$clists[$cset]
 # egrep $clist |...

end

#...........

velvetg:
Final graph has 192 nodes and n50 of 150, max 557, total 10817, using 129380/130596 reads

oases:
Counted 12 mRNA loci
Finished heuristic approach, used 130444/130596 reads

Locus_1_Transcript_1/5_Confidence_0.286 w=674
Locus_1_Transcript_2/5_Confidence_0.429 w=1599
Locus_1_Transcript_3/5_Confidence_0.571 w=2252    << 18S ?
Locus_1_Transcript_4/5_Confidence_0.214 w=480
Locus_1_Transcript_5/5_Confidence_0.500 w=2231

Locus_2_Transcript_1/5_Confidence_0.375 w=3303
Locus_2_Transcript_2/5_Confidence_0.250 w=2736
Locus_2_Transcript_3/5_Confidence_0.750 w=4813     
Locus_2_Transcript_4/5_Confidence_0.562 w=1796
Locus_2_Transcript_5/5_Confidence_0.750 w=4839    << 23S ?

Locus_4_Transcript_1/1_Confidence_1.000 w=269    #.. intergenic ?
Locus_5_Transcript_1/1_Confidence_1.000 w=348
Locus_6_Transcript_1/1_Confidence_1.000 w=207
Locus_8_Transcript_1/1_Confidence_1.000 w=173


..............
peaaphid asm2 matches (lots)

# 1
peaaphid_rRNAl2t5       Scaffold922     97.57   1151    24      4       80      1228    169     1317    0.0     2056
peaaphid_rRNAl2t5       Scaffold922     99.89   3526    0       3       1318    4839    1407    4932    0.0     6954

peaaphid_rRNAl1t3       Scaffold922     99.78   921     1       1       1       920     5480    6400    0.0     1808
peaaphid_rRNAl1t3       Scaffold922     99.91   1165    0       1       1011    2174    6491    7655    0.0     2299


# 2
peaaphid_rRNAl1t3       Scaffold922     99.89   920     0       1       1       920     43581   44499   0.0     1814
peaaphid_rRNAl1t3       Scaffold922     100.00  368     0       0       1011    1378    44590   44957   0.0      729

peaaphid_rRNAl2t5       Scaffold922     98.86   615     4       3       42      653     37067   37681   0.0     1160
peaaphid_rRNAl2t5       Scaffold922     99.92   1221    0       1       1517    2736    39877   41097   0.0     2410
peaaphid_rRNAl2t5       Scaffold922     98.70   1152    1       14      2905    4044    38677   39826   0.0     2151


# peaaphid_rRNAl1t3       Scaffold922     99.54   650     0       3       1527    2174    26628   27276   0.0     1261
# peaaphid_rRNAl2t5       Scaffold922     99.38   650     4       0       3959    4608    23999   24648   0.0     1257
# peaaphid_rRNAl2t5       Scaffold922     99.23   651     3       2       42      690     11395   12045   0.0     1249

....

# 2
peaaphid_rRNAl1t3       Scaffold196     100.00  920     0       0       1       920     164308  163389  0.0     1824
peaaphid_rRNAl1t3       Scaffold196     98.54   1165    16      1       1011    2174    163298  162134  0.0     2173

peaaphid_rRNAl2t5       Scaffold196     95.63   1143    35      8       80      1220    169549  168420  0.0     1853
peaaphid_rRNAl2t5       Scaffold196     96.55   493     16      1       1318    1809    168321  167829  0.0      840
peaaphid_rRNAl2t5       Scaffold196     98.32   2970    46      3       1874    4839    167825  164856  0.0     5487

# 4  Scaffold96:43803-49515
peaaphid_rRNAl1t3       Scaffold196     99.91   1165    0       1       1011    2174    44967   43803   0.0     2299
peaaphid_rRNAl1t3       Scaffold196     100.00  920     0       0       1       920     45977   45058   0.0     1824

peaaphid_rRNAl1t3       Scaffold196     99.20   503     1       3       419     920     49515   49015   0.0      961
peaaphid_rRNAl1t3       Scaffold196     98.80   502     0       6       1011    1506    48924   48423   0.0      942
peaaphid_rRNAl2t5       Scaffold196     99.89   1843    0       2       2999    4839    48372   46530   0.0     3635


# peaaphid_rRNAl2t5       Scaffold196     99.56   1373    1       5       3425    4792    8568    7196    0.0     2668
# peaaphid_rRNAl2t5       Scaffold196     99.41   681     0       4       4163    4839    4063    3383    0.0     1314



..........
peaaphid_rRNAl1t3 NRmatch

>gb|AY216697.1|  Toxoptera citricida 18S ribosomal RNA gene, complete sequence
(brown citrus aphid)
Length=2480

 Score = 3674 bits (4074),  Expect = 0.0
 Identities = 2100/2139 (99%), Gaps = 2/2139 (0%)
 Strand=Plus/Minus
................


>peaaphid_rRNAl1t3 len=2252 name=18S rRNA pea aphid (reversed, mostly complete)
TCAGTGTAGCGCGCGTGCGGCCCAGAACATCTAAGGGCATCACAGACCTGTTATCGCTCA
GTCTCGTGCGGCTATGTTCGTCCGCCGCCTGTCCCTCTAAGAAGAGTTTAAGCTCCTGGG
AGCCGGCGGTAGCCCTAGTAACGTATCGTGATCCGCCGGCGACGGCCGCGAACACGGCGC
GCTTCACCACGGAGCCCGACGCACGACGGACGCAGGCCGACCGGAGGTTGCCCCCCGGCC
GACATACGCCCGCCGAACGCCGGGACGCGGATGGCGCGGCCGCGCCCGACGACCGCCGTG
GACGACGGGCGAACGCGTTCGGGGATACCGGGCCGAGCCGACGGGTACGCGAACACTGAC
GGCCGAAACCGCCAGCGCCGCGCACACGCCGACCCGGGTTTACCCGCCTAGTTAGCAGGA
CAGAGTCTCGTTCGTTATCGGAATTAACCAGACAGATCGCTCCACCAACTAAGAACGGCC
ATGCACCACCACCCACCGAATCAAGAAAGAGCTCTCAATCTGTCAATCTTTCCGGTGTCC
GGGCCTGGTGAGGTTTCCCGTGTTGAGTCAAATTAAGCCGCAGGCTCCACTCCTGGTGGT
GCCCTTCCGTCAATTCCTTTAAGTTTCAACTTTGCAATCATACTTCCCCCGGAACCGAAA
AGCTTCGGTTTCCCGGAAGCTGCCCGCCGGGTCGTTAATGAAACGCCGGCGGATCGCTAG
CTGGCATCGTTTACAGTTAGAACTAGGGCGGTATCTGATCGCCTTCGAACCTCTAACTTT
CGTTCTTGATCATACGAGAACGTACTTGGCAAATGCTTTCGCGTCAGTTCGTCTCGAGAC
GATCCAAGAATTTCACCTCTAACGTCTCGGTACGAATGCCCCCGCCCGTCTCTGTTGATC
ATTACCTCCGGTCCCGAAAACCGGCCCGGCGGGACGCGCGGGCGACTGGCGCCCGCGGCC
CGGCGGCCCGCGCACGGAAACGCCCCGGAGGGCGATTTCGCGCGCCCGCGAAGGGCGGAG
ATGCGCGGGACCGAGGTCTTGTTCCATTATTCCATGCGACCAGTATTCAGGGCCTTTTGA
CGAGACGGCCGTGAAGCCGCCCCGCCAGATTCGAGCCTGCTTTGAGCACTCTAATTTGTT
CAAAGTAAACGTGTCGGCCCGCCGACGGCACTCGGTGAAGAGCACCGCGCAGCAAGATTG
GAGTAGGCGGCCGCCGTCGTCGAACACCGACGGCCGCGCGACGCGTGGCCGCGCGGCGCG
CCGGAAGCACGAGACACGTGTCCGCCTGCCGACAATACGTCCGGCCGACGTGCCGGTAAC
TAACACCCCGAGACGGCTGACGAGCGCGACGGCGACCGCGCCGACGGGACACGTGGTCCC
GCGACACGGCCGGCCGCGACACGACGGACCGCCAGGGTGGTCCGGCACCGTCCAGACACA
GATCCGACTACGAGCTTTTTAACCGCAACAACTTTAATATACGCTATTGGAGCTGGAATT
ACCGCGGCTGCTGGCACCAGACTTGCCCTCCAATTGATCCTCGTTTAAAGGTTTTAAAGT
GTTCTCATTCCGATTACGGGCCTCGGATGAGTCCCGTATCGTTATTTTTCGTCACTACCT
CCCCGTTCCGGGAGTGGGTAATTTGCGCGCCTGCTGCCTTCCTTGGATGTGGTAGCCGTT
TCTCAGGCTCCCTCTCCGGAATCGAACCCTGATTCCCCGTTACCCGTTACCACCATGGTA
GGCATGGAACCTACCATCGACAGTTGATAAGGCAGACATTTGAAAGATGCGTCGCCGGTA
CGGAGACCGTGCGATCAGCTCGAAGTTATTCAGAGTCACCAGGTCTTTGCGCGGGCCGCG
ATCGGGACGCGCACGAAGCGCGCCGCGACGGGCCGGTTTTGATCTAATAAAAGCGTTCCT
CCCGCGAGCGACGCCCGAGGGCGCCGCGACGGTCGGAACTCTGTCGGCATGTATTAGCTC
TAGAATTACCACAGTTATCCAAGTAACTTGGGTACGATCTAAGGAACCACAACTGATTTA
ATGAGCCTTTCGCGGTTTCACCTTAATGCGGCTTGCACTGAGACATGCATGGCTTAATCT
TTGAGACAAGCATATGACTACTGGCAGGATCAACCAGGGATCTCGGTTTTACAACAACTG
CGGCCCGCGACGGCGCGGACACTCGTCGCGGGACGCGTGCGTCCGAGCGCCGGTCCGCCC
GTGCGGGCGAACCGACGCCGAGGAAGGTCGGA

>peaaphid_rRNAl2t5   len=4839 name=23S rRNA pea aphid (reversed, mostly complete?)
TTCGCGGTACGCGGCGCCGGAGCGCCACGCACCGATCATGCGCGCCCGACGCGTGGAACG
CGACCGGTGGGCCGATCCGTTCACGCATCGGTCACTGGGGGTCCTGTCCAACCGACAAGA
CGAACCCCCGAGGCAAAGGGCAGTCTTAACAGATCGCAGCGTGGTAACTGCTCTGCCGAG
TACAACACCCAGCCCGGTACTTAAGTCGTCTGCAGACGATTCCGAAAACCCACACCGTTT
GCCGCGGGATCACCGCGTTCACCGTTGATGCACGGCCAGCGGGCCCGAGAGCCCGCGCGG
CCGAGAATCCGGCTCGTCGTGGACCCGAAGGTCCGTGCTTTGACGCCTTCCCGATGTATA
CTGGGTTCTCCTCCGCCGTACGGACGATCGTTTTCCCGAGACCGGCCCGAAAGCCGGCCC
GGCTGTTTTCTCCAGAGTGGATACGGCCTTAGAGGCGTTCAGGCGTAATCCAACGGATGG
TAGCCTCGCACCAACGCCCGCTCGGGCGAGTGCCGAACCAAATGTCCGAACCTGCCGTTC
CTCTCGTACTGGGCAGGATTACTATCGTAACGACTGCCGCCGTAAAGGCGGTTTCGATCG
CGTGCGACCTTGCATCAGTAGGGTAAAACTAACCTGTCTCACGACGGTCTAAACCCAGCT
CACGTTCCCTTAAGCGGGTGAACAATCCGACGCTTGGCGAATTTTGCTTCGCAATGATAG
GAAGAGCCGACATCGAAGGATCAAAAAGCGACGTCGCTATGAACGCTTGGCCGCCACAAG
CCAGTTATCCCTGTGGTAACTTTTCTGACACCTCTCGCTGAAAACTCTTCAGCGGACGAG
AGGATCGAGAGGCCGATGCTTTCGCAGTCCCTACGCGTACTGAGCGTCCGGGATCAAGCC
AGCATTTGCCCTTTTGCTCTACGCGAGGTTTCCGTCCTCGCTGAGCTGGCCTTAGGACAC
CTGCGTTATTTTTTGACAGATGTACCGCCCCAGTCAAACTCCCCACCTGGCGGTGTCCTC
GGAATACGGATCGCACCAGGGACCGGGCCCGCGGATACGCCGCGAACGCGGACGGCGACT
TTTGACGGTCGGCCGTACGCGGACGGACGCGGTCCGCGGTTTGACGCCCGGATCCCTCGG
CTGGTGCTTAACGCTACCGGAATATCAACCGCGGCCCGGCACGAACGGGCACAGGCCGAC
GGCACGGCGACCACGCCGCGCGCCAGTAGCGCGCCGGCGAACCGACGGCAAACCGACGAC
GCGACGGGACGACGACCGCCGGCCGGACGCGCCCGCGGCCGGAACCGCGCGTTCCGCTCT
ACCGAGTAAGTGGGGAAACGATGCGAGTAGTGGTATTTCAAGGTCGGCCCGGAGACGAAC
GGCCGAAACCGCCCGCCTTGGTCCGGGTCTCCCACGTATGCTACACCTCGGCATGTCTCC
GAACAATGCCAGATTAGAGTCAAGCTCAACAGGGTCTTCTTTCCCCGCTGATTTTTCCAA
GCCCGTTCCCTTGGCTGTGGTTTCGCTAGATAGTAGATAGGGACAGTGAGAATCTCGTTA
ATCCATTCATGCGCGTCACTAATTAGATGACGAGGCATTTGGCTACCTTAAAAGAGTCAT
AGTTACTCCTGCCGTTTACCCGCGCTTGCTTGAATTTCTTCACGTTGACATTCAGAGCAC
TGGGCAGAAATCACATCGCGTCAACACCCGTCCCGGGCCATCGCGATGCTTTGTTTTAAT
TAGACAGTCGGATTCCCCTGGTCCGTGCCAGTTCTGAGTTGACCGTTACATGGCTGTCGA
TCCGGCTACGCGACCGACGGCGCGGCGGCCGGACCGCCCGGTGACGGCGAACCGACACCA
AGGCGGTCGGCCGGCCGCGACGCGGCCAGCGACCGAAGCTCGGTGGTTCCACGGTCGGCG
GACGGCGACGGGCCCGCGGCCGCCTCGAGGCTCCCGTTCCGAAGACCGGGCGCACGGGCG
ACGGCCGGGGTGTCGCCAAAAACGCCTACGCTTCTACGACGACCCGAACCCGGCAGCCAC
GCTCCTCAGAGCCAATCCTTATCCCGAAGTTACGGATCGGTTTTGCCGACTTCCCTTACC
TACATTATTCTATGCGGCTAGAGGCTGCTCACCTCGGAGACCTGCTGCGGATATCGGTAC
GAACCGACGCGAAGACTCCGCGTGGCCCTCTCTCGAATTTTCAAGGTCCGCGTGGGGATC
ACGGACACCGCCGCAACGAGCGGTGCTCTTCGCGCTCGCGTCCCTATCGCCCGGCTAGAG
GATTCCAGGGACGTACAACGCTCACAGAGAAAAGAGAACTCTACCCAGATCCCCCGGCGG
CTTCTTCGAGTTCATTCTGGTTACCCAGACGAGACAAAAGAGCCCCGAACACTAGGGAGC
GGTTCCGCGTCGGGTTCCGGAATACGAACCGGATTCCCTCTCGCCCCAAGGGCGATTCGA
AAACGCCTTCGCCCGTGGTACGAGGATCTCTCCCCGGGCTTAGGATCGACTGACTCTTGG
ACAACGGCTGTTCACAAGAAACCCTTCTCCACGGCAGCCCCCGAGGGCCCCTCTCGAGTA
TTTGCTACTACCACCAAGATCTGCACCGACGGCGGCTCCAGGCGGTCTCGCGACCTGCCC
TTCGACGCACACCGCCGCGCCATCCTACTCGTCGAGGCTTGCCGGGACCGGCCGCCGCGA
GCGACGACTGGCCCCACTATGCCGTCGACGGCCGAGTATAGGCACGACGCTCCAGCGCCA
TCCATTTTCAGAGCTAGTTGCTTCGGCAGGTGAGTTGTTACACACTCCTTAGCGGATTCC
GACTTCCATGGCCACCGTCCTGCTGTCATCAGCAACCAACGCCTTTCATGGTCTCTGAAT
ACGCGTCGATTTCGGCGCCTTAACTCGGCGTTCGGTTCATCCCGCAGCGCCAGTTCTGCT
TACCAAAAATGGCCCACTAAGCGCCTAGGGTTCCGTCGCCGGCTTCGCACGCGGTTCACG
CGGTGTACCAGGTAAAGCCGGCGATCTCACTCATTTATAGTTTGAGAATAGGTTGAGGTC
GTTTCGGCCCCAATGTCCTCTAATCATTCGCTTTACCGTATGAGACATCCTCCGTGTTAC
GAGCGTGTTGTTACCGCGCCGTACAGCGCGGAGTCCCACCGGCGTCCGTACGCCCGCGAA
CGGGCGGGACACATCGGCGCCAGCTATCCTGAGGGAAACTTCGGATGGAACCAGCTACTA
GATGGTTCGATTAGTCTTTCGCCCCTATACCCAGCTCAGACGATCGATTTGCACGTCAGA
ATCGCTGCGGACCTCCATCAGGGTTTCCCCTGACTTCATCCTGGCCAGGCATAGTTCACC
ATCTTTCGGGTACCAACGTGTGCGCTAAGGGTGCGCCCCCAGCCGGCCGAAGCCGACTTG
GCGGAGACGCCCCCGGACTGCGGACCAGCGCGACTTTGAACGCCGGTGGGTTCGGTCGCT
AGGCCATCGTCCGCACCGTAAACGGTTCACTTTCATTGCGCCAAACGTGGTTTTTCGTAA
GATCACCGTTGACTCGCGCACACGTTAGACTCCTTGGTCCGTGTTTCAAGACGGGTCGGA
AAGTAGCTCGAAACACGTCGCCGACCGACGGGACCCCGCCTGCGACGGGGCCGCGTCGAC
GGTTCTTGCCGTCGGGCGAGCCGACCATGAGCACCGGGGCTCGGCGCGAACGGTCACGCC
GAAACGCGACCGACGCGACGAACGTCACCCGGGCCATTCGGCCGGCGCCCAACGGGTCGC
GACGTCCGCTAACCGAGCGAAAGATCCGGCCGGCGGTTAGACCGACCGTGAATTCGCCCG
GCGGGGATTCGCGAGCTCTGTCCGTTTACAACCGAGCGGTTTCACGTTCTTATGAACTCT
CTCTTCAAAGTTCTTTTCAACTTTCCCTCACGGTACTTGTAAACTATCGGTCTCGTGGCC
GTATTTAGCCTTAGATGGAGTTTACCACCCGCTTCGGGCTGCACTCTCAAGCAACCCGAC
TCGAAGGCGCGGTCCGATCCCGGAACGCTGTGACGGCCTCTACTGGCCTGGCACCATCTG
CGGGCCATGGCCCCGTTCAAGGGAGACTTGGACCGTCACGTGCGCCCCGGCAACGAGGCC
GGCCCGTACGCCACAACTCTCATACGTGCCGCGACACAAAAGTGCGCGGCGAGATTCGGC
GATGGGCTTTTCCCGGTTCGCTCGCCGCTACTAAGGGAATCACGGTTGTTTTCTTTTCCT
CCGCTTATTAATATGCTTAAATCCGGCGGGTAGTCCCGCCTGATCTGAGGTCGGAACGTA
TCGGTTGTTTGTTTCTTAAATTACACGGCAGCCTGCGGTCCGCGCTCCGCCCGTACGGTG
ACCGCCCGCGTCCGCTTCTCGGAAGCCGGAGAGGTCTCGTGGACTCTCAGCCGGCCCGGC
CGCGCTCTCGCGCGGCGGGCGGACGGGGACGTGCCGTTTGACGACGCGAGAACCCGTGAC
TGTTTGCCGTCACAACTTGGGCGGACGACCGGAGGCGGCCGCGCGAGCGGCTCTCCCCGG
TCGAACCGCCAATCTCGCACCACCGGAGCGGCCGACGGGCGTGCCGTCGGACGCCTGCCG
GCGGCGCATCGGTCTGGCCGAGTATCCGGTATACGACCCTCAGACAGGCGTGGCCCGGGA
CCCGCGGGTCACCGAGGCCGCAATGTGCGTTCGACTGGTCGATGTTCGTATAACCTGCGG
ATTACACGACGACGCGCAGATAGCTGCGGTCTTCATCGATCCACGAGCCAAGTGATCCGC
CGCTCGGTGTCGGTTGTTTTTTGTTTTGTTCGACGAAAC

>peaaphid_rRNAl4t1  NRmatch=18S rRNA partial, 100% Ascomycota fungal match
GCTCGAATACATTAGCATGGAATAATAGAATAGGACGTGTGGTTCTATTTTGTTGGTTTC
TAGGACCGCCGTAATGATTAATAGGGATAGTCGGGGGCATCAGTATTCAAGCGTCAGAGG
TGAAATTCTTGGATTGCTTGAAGACTAACTACTGCGAAAGCATTTGCCAAGGATGTTTTC
ATTAATCAGTGAACGAAAGTTAGGGGATCGAAGACGATCAGATACCGTCGTAGTCTTAAC
CATAAACTATGCCGACTAGGGATCGGGCG
>peaaphid_rRNAl5t1  NRmatch=18S rRNA partial, 100% Ascomycota fungal match
GCCGTTCTTAGTTGGTGGAGTGATTTGTCTGCTTAATTGCGATAACGAACGAGACCTTAA
CCTGCTAAATAGCCAGGCCCGCTTTGGCGGGTCGCCGGCTTCTTAGAGGGACTATCGGCT
CAAGCCGATGGAAGTTTGAGGCAATAACAGGTCTGTGATGCCCTTAGATGTTCTGGGCCG
CACGCGCGCTACACTGACAGAGCCAACGAGTTCATTTCCTTGTCCGAAAGGTCTGGGTAA
TCTTGTTAAACTCTGTCGTGCTGGGGATAGAGCATTGCAATTATTGCTCTTCAACGAGGA
ATGCCTAGTAAGCGCATGTCATCAGCATGCGTTGATTACGTCCCTGCC
>peaaphid_rRNAl6t1  NRmatch=28S ribosomal RNA partial 3'41bp, various weak match
CGGTCTGCCCACGTGACAAGGCCGGGAATATGACGGTGACTCGGTATTAGCAAACAACGG
CATAGCCGTTGTCAGTATCCAATCTGCCGCGAGCGGGGACGCTCCGAGCCTTCGGCGAGC
GCCGCGGACTGTCAACCCCGCGACGCCCAGCGAAACCCCAGCCCGATGGCTGTGGTTTCG
CTAGATAGTAGATAGGGACAGTGAGAA
>peaaphid_rRNAl8t1  NRmatch=none
GTAGATAGGGACATTCGTAATCTTGTAACGGGAGGTATAAGGCCGAACGTCCCGGTGGGG
TGACCCCCTTTTTCCAGATAAGCGGGACTGGCTGCTTCTCGCCTTTGAGTGGCGCAGTCG
TCCCTCCCAGACGAGTACACTGGGCGCCTGGGCGCCAGACAGTGTTTCGCTGG

=cut

