#!/usr/bin/env perl
# subsetids.pl

use strict;
# my $outdir=$ENV{outdir} || ""; $outdir.="/" if($outdir);
# $outdir="" unless($outdir and -d $outdir);

my $bindir=$ENV{bindir} || ""; $bindir.="/" if($bindir);
my $samtools=`which "${bindir}samtools"`; chomp($samtools);
die "missing ENV{bindir}/samtools" unless( -x $samtools);
# set TMPDIR for sort? and or do sort elsewhere.

my $vprefix="dmg3."; # 
my ($subfile)= grep {not m/\.bam/} @ARGV;
my @bamlist= grep {/\.bam$/} @ARGV;

# subfile3: tag, bpsize, nreads, nparts, parts : last part has 4000+ items = small scafs
# sub1	7262140	58458046	2	scaffold00512 scaffold00024 
# sub2	5060408	70964689	2	scaffold01581 scaffold01361 

open(L,$subfile) or die "open subfile=$subfile"; 
while(<L>) {
  chomp; my ($part,$bsize,$nread,$nparts,$sublist)=split"\t"; 
  $sublist=~s/,/ /g;

  my $err;
  my $ib=0;
  foreach my $bam (@bamlist) {
    (my $nab=$bam) =~ s/.bam//; 
    my $tmpname="$nab.$part.ids"; $ib++;
    my $cmd1="$samtools view $bam $sublist | cut -f1 | sort -u > $tmpname"; # : sort -u  later ?
    $err= system($cmd1);
    if($err) { die "ERROR: $err on '$cmd1'"; }
  }
}

# sub subsetbams
# my $outname="$nab.$part.bam";
# my $cmd1="$samtools view -b -o $outname $bam $sublist"; # : sort -u  later ?
# .. or samtools view -1 = fast compress, newer samtools
# env row=10 ./bamsub.pl &
# env row=11 ./bamsub.pl &

sub bamsub  # to be called by shell script for subsets..
{
# export samtools=$bindir/samtools
# export bam=allbam3.bam
# export subtab=subset50m.tab
# env row=9  bamsub.pl
  my $bam=$ENV{bam};
  my $samtools=$ENV{samtools};
  my $subfile=$ENV{subtab};
  my $row=$ENV{row};
  (my $nab=$bam) =~ s/\.bam$//; 
  my $i=0; 
  open(L,$subfile) or die "open subfile=$subfile"; 
  while(<L>) {
    chomp; my ($part,$bsize,$nread,$nparts,$sublist)=split"\t"; 
    next unless($sublist); $i++; $sublist=~s/,/ /g;
    if($i == $row) { 
      my @sublist= sort split" ",$sublist; 
      my @opt=("view", "-1", "-o", "$nab.$part.bam", $bam);
      exec( $samtools, @opt, @sublist);
      exit;
    }
   }
}

sub subsetfq
{
  # my($fq,@idf)=@_;
  my $outdir= $ENV{outdir} || ""; $outdir="" unless($outdir and -d $outdir); 
  my @fqs= grep {m/\.gz$/} @ARGV;
  my @ids= grep {/\.ids$/} @ARGV;
  my %ids=(); my %subs; my $i=0;
  # add skipid set ..
  foreach my $idf (@ids) { 
    my ($sub)= $idf =~ m/(\w+)\.ids$/; ## m/(sub\d+)/;
    $i++; $sub=$i unless($sub);
    $subs{$sub}++;
    open(F,$idf); while(<F>) { chomp; $ids{$_}= $sub; } close(F);
    ## Out of memory! for %ids in 8GB system, 34,088,210 ids
    ## ** read can be in 2+ subsets.. what to do? how many? 
    ##   134 of 2.5mil each 16-CR3_pe006_1.sub[12].ids .. up to 700 for sub1 x small contigs
  }
 
  my @subs=sort keys %subs;
  push @subs, "nomap"; # should be option?
  foreach my $fq (@fqs) {
    my %subh;
    my $fs= $fq; 
    if($outdir) { ($fs)= $fq =~ m,([^/]+)$,; $fs="$outdir/$fs"; }
    $fs=~s/.gz//; 
    
    foreach my $s (@subs) {
      my $sh; 
      open( $sh,">$fs.$s") or die "$fs.$s";
      $subh{$s}=$sh;
      }
    open(F,"gunzip -c $fq |") or die "$fq"; 
    while(<F>) {
      my($ih,$is,$ij,$iq); 
      $ih=$_; $is=<F>; $ij=<F>; $iq=<F>;
      my($id)= $ih =~ m/^.(\S+)/;  $id =~ s,/[12]$,,;
      my $sid= $ids{$id} || "nomap"; #? option
      if( my $sh= $subh{$sid} ) { print $sh $ih,$is,$ij,$iq; }
    }
    foreach my $s (@subs) { if(my $sh= $subh{$s}) { close $sh; } }
    close(F);
  }  
}

# input: pair1/2.fq files; output: velpair2.fa, removing extra n duplicates
# use system sort for dup removal
# ../pairduprm.pl 16-CR3_pe006_?.sanfastq.sub8
#pairduprm: fastq pairs 16-CR3_pe006_1.sanfastq.sub8,16-CR3_pe006_2.sanfastq.sub8 to 16-CR3_pe006.sanfastq.sub8.fa2
#delete tempfiles: rm 16-CR3_pe006_p.sanfastq.sub8.unsort 16-CR3_pe006_p.sanfastq.sub8.sort
#pairduprm:  inpairs=1977407; outpairs=1846265; dropdups=131142
#..
# ../pairduprm.pl  16-C_R1*.sub8
#pairduprm: fastq pairs 16-C_R1_1.txt.sub8,16-C_R1_2.txt.sub8 to 16-C_R1.txt.sub8.fa2
#delete tempfiles: rm 16-C_R1_p.txt.sub8.unsort 16-C_R1_p.txt.sub8.sort
#pairduprm:  inpairs=716570; outpairs=596852; dropdups=119718
# ../pairduprm.pl c6.sub8.fq?.gz
#pairduprm: fastq pairs c6.sub8.fq1.gz,c6.sub8.fq2.gz to c6.sub8.fq1.fa2
#pairduprm:  inpairs=6776270; outpairs=4723926; dropdups=2052344
# ../pairduprm.pl  c3.sub8.fq?.gz
#pairduprm: fastq pairs c3.sub8.fq1.gz,c3.sub8.fq2.gz to c3.sub8.fq1.fa2
#pairduprm:  inpairs=8340390; outpairs=4918626; dropdups=3421764
#  ../pairduprm.pl  c4.sub8.fq?.gz
#pairduprm: fastq pairs c4.sub8.fq1.gz,c4.sub8.fq2.gz to c4.sub8.fq1.fa2
#pairduprm:  inpairs=8153103; outpairs=5780801; dropdups=2372302
# ../pairduprm.pl cl.sub8.fq?.gz
#pairduprm: fastq pairs cl.sub8.fq1.gz,cl.sub8.fq2.gz to cl.sub8.fq1.fa2
#pairduprm:  inpairs=6724374; outpairs=5881692; dropdups=842682
# /home/ux455375/scratchg/chrs/daphmag/rnas/tri2set
# c3.sub8.fq1.fa2:9837252
# c4.sub8.fq1.fa2:11561602
# c6.sub8.fq1.fa2:9447852
# cl.sub8.fq1.fa2:11763384



sub pairduprm
{
  # export TMPDIR=xxx for sort
  my ($nout,$ndup,$npair,$idoff)=(0)x10;
  my $MAXDUP= $ENV{max} || 3; 
  my $MAXIDERR= $ENV{errmax} || 99;
  my $outfile= $ENV{out} || "";
  my ($fq1,$fq2)= @ARGV;
  if($fq2 =~ m/_1/ and $fq1 =~ m/_2/) { ($fq1,$fq2)=($fq2,$fq1); }
  elsif($fq2 =~ m/fq1/ and $fq1 =~ m/fq2/) { ($fq1,$fq2)=($fq2,$fq1); }
  my $ftmp=$fq1; $ftmp=~s/.gz$//; $ftmp =~ s/_1/_p/; 
  my $ftmp2= "$ftmp.sort"; $ftmp .=".unsort"; 
  unless($outfile) { $outfile=$fq1; $outfile=~s/.gz$//; $outfile=~s/_1//; $outfile.=".fa2"; }
  warn "#pairduprm: fastq pairs $fq1,$fq2 to $outfile\n";
  
  if($fq1=~/\.gz$/ or $fq2=~/\.gz$/) {
  open(F,"gunzip -c $fq1|") or die "$fq1"; open(R,"gunzip -c $fq2|") or die "$fq2";
  } else {
  open(F,$fq1) or die "$fq1"; open(R,$fq2) or die "$fq2";
  }
  open(O,">$ftmp") or die $ftmp;
  while(<F>) {
    my($ih,$is,$ij,$iq); 
    my($rh,$rs,$rj,$rq); 
    $ih=$_; $is=<F>; $ij=<F>; $iq=<F>;
    $rh=<R>; $rs=<R>; $rj=<R>; $rq=<R>;
    chomp($is); chomp($rs); chomp($ih); chomp($rh); 
    unless($is and $rs) { warn "#EOFearly: lh=$ih, rh=$rh for lfq=$fq1, rfq=$fq2\n"; last; } # die? last?
    my($id)= $ih =~ m/^.(\S+)/;  $id =~ s,/[12]$,,;
    my($rd)= $rh =~ m/^.(\S+)/;  $rd =~ s,/[12]$,,;
    if($id ne $rd) {  warn "#IDmismatch: lh=$id, rh=$rd\n"; last if(++$idoff > $MAXIDERR); next; } 
    print O "$id\t$is,$rs\n";  $npair++;
  }
  close(O); close(F); close(R);
  my $err= system("sort -k2,2 -k1,1 -o $ftmp2 $ftmp");
  my $np=0; my $lpair="";
  open(F,$ftmp2); open(O,">$outfile") or die $outfile;
  while(<F>){
    my($id,$pairs)=split;
    if($pairs eq $lpair) { $np++; } else { $np=1; }
    if($np>$MAXDUP) { $ndup++; } 
    else {  
      my($ls,$rs)=split",",$pairs; 
      print O ">$id/1\n$ls\n>$id/2\n$rs\n"; $nout++; 
      }
    $lpair=$pairs;
  }
  close(O); close(F);
  system("rm $ftmp $ftmp2\n");# warn "#delete tempfiles:"; 
  warn "#pairduprm:  inpairs=$npair; outpairs=$nout; dropdups=$ndup\n";
}

__END__

=item ids x fq matchup...
    
    subsetq.pl fq6c/*.txt.gz bam3/subsets/16-CR?_1-dmag2.sub*.ids
    need ?? 16GB per idset ?

[ux455375@gordon-ln1 rnas]$ ls bam3/subsets/*.sub1.ids
bam3/subsets/16-1R1_1-dmag2.sub1.ids	    bam3/subsets/3-CR1_1-dmag2.sub1.ids
bam3/subsets/16-1R2_1-dmag2.sub1.ids	    bam3/subsets/3-CR2_1-dmag2.sub1.ids
bam3/subsets/16-1R3_1-dmag2.sub1.ids	    bam3/subsets/3-CR3_pe012_1-dmag2.sub1.ids
bam3/subsets/16-8R1_1-dmag2.sub1.ids	    bam3/subsets/4-1R1_1-dmag2.sub1.ids
bam3/subsets/16-8R2_1-dmag2.sub1.ids	    bam3/subsets/4-1R2_1-dmag2.sub1.ids
bam3/subsets/16-8R3_1-dmag2.sub1.ids	    bam3/subsets/4-1R3_pe001_1-dmag2.sub1.ids
bam3/subsets/16-CR1_1-dmag2.sub1.ids	    bam3/subsets/4-8R2_1-dmag2.sub1.ids
bam3/subsets/16-CR2_1-dmag2.sub1.ids	    bam3/subsets/4-8R3_1-dmag2.sub1.ids
bam3/subsets/16-CR3_1-dmag2.sub1.ids	    bam3/subsets/4-CR1_1-dmag2.sub1.ids
bam3/subsets/16-CR3_pe006_1-dmag2.sub1.ids  bam3/subsets/4-CR2_1-dmag2.sub1.ids
bam3/subsets/3-1R1_1-dmag2.sub1.ids	    bam3/subsets/4-CR3_1-dmag2.sub1.ids
bam3/subsets/3-1R1_pe015_1-dmag2.sub1.ids   bam3/subsets/merge-3-1R3_1-dmag2.sub1.ids
bam3/subsets/3-1R2_1-dmag2.sub1.ids	    bam3/subsets/merge-3-8R3_1-dmag2.sub1.ids
bam3/subsets/3-1R3_pe004_1-dmag2.sub1.ids   bam3/subsets/merge-3-CR3_1-dmag2.sub1.ids
bam3/subsets/3-8R1_1-dmag2.sub1.ids	    bam3/subsets/merge-4-1R3_1-dmag2.sub1.ids
bam3/subsets/3-8R2_1-dmag2.sub1.ids	    bam3/subsets/merge-4-8R1_1-dmag2.sub1.ids
bam3/subsets/3-8R3_pe020_1-dmag2.sub1.ids

fq01long:
4-1R3_pe001_1.sanfastq.gz  4-1R3_pe001_2.sanfastq.gz

fq04long:
3-1R3_pe004_1.sanfastq.gz  3-1R3_pe004_2.sanfastq.gz

fq06long:
16-CR3_pe006_1.sanfastq.gz  16-CR3_pe006_2.sanfastq.gz

fq12long:
3-CR3_pe012_1.sanfastq.gz  3-CR3_pe012_2.sanfastq.gz

fq15long:
3-1R1_pe015_1.sanfastq.gz  3-1R1_pe015_2.sanfastq.gz

fq20long:
3-8R3_pe020_1.sanfastq.gz  3-8R3_pe020_2.sanfastq.gz

fq31:
101022_3-1R3_1.txt.gz  101118_3-1R3_1.txt.gz  3-1-R1_1.txt.gz  3-1-R2_1.txt.gz	3-1R3_1.txt.gz
101022_3-1R3_2.txt.gz  101118_3-1R3_2.txt.gz  3-1-R1_2.txt.gz  3-1-R2_2.txt.gz	3-1R3_2.txt.gz

fq38:
101015_3-8R3_1.txt.gz  101022_3-8R3_2.txt.gz  3-8-R1_1.txt.gz  3-8-R2_2.txt.gz
101015_3-8R3_2.txt.gz  101118_3-8R3_1.txt.gz  3-8-R1_2.txt.gz
101022_3-8R3_1.txt.gz  101118_3-8R3_2.txt.gz  3-8-R2_1.txt.gz

fq3c:
101015_3-CR3_1.txt.gz  101022_3_CR3_2.txt.gz  3-C-R1_1.txt.gz  3-C-R2_2.txt.gz
101015_3-CR3_2.txt.gz  101118_3-CR3_1.txt.gz  3-C-R1_2.txt.gz
101022_3_CR3_1.txt.gz  101118_3-CR3_2.txt.gz  3-C-R2_1.txt.gz

fq41:
101011_4-1_R1_1.txt.gz	101015_4-1R3_2.txt.gz  101118_4-1R3_1.txt.gz  4-1-R2_2.txt.gz
101011_4-1_R1_2.txt.gz	101022_4-1R3_1.txt.gz  101118_4-1R3_2.txt.gz
101015_4-1R3_1.txt.gz	101022_4-1R3_2.txt.gz  4-1-R2_1.txt.gz

fq48:
100827_4-8-R1_1.txt.gz	101001_4-8-R1_1.txt.gz	4-8-R2_1.txt.gz  4-8-R3_1.txt.gz
100827_4-8-R1_2.txt.gz	101001_4-8-R1_2.txt.gz	4-8-R2_2.txt.gz  4-8-R3_2.txt.gz

fq4c:
101001_4-C-R1_1.txt.gz	4-C-R2_1.txt.gz  4-C-R3_1.txt.gz
101001_4-C-R1_2.txt.gz	4-C-R2_2.txt.gz  4-C-R3_2.txt.gz

fq61:
16-1_R1_1.txt.gz  16-1-R2_1.txt.gz  16-1-R3_1.txt.gz
16-1_R1_2.txt.gz  16-1-R2_2.txt.gz  16-1-R3_2.txt.gz

fq68:
16-8-R1_1.txt.gz  16-8-R2_1.txt.gz  16-8-R3_1.txt.gz
16-8-R1_2.txt.gz  16-8-R2_2.txt.gz  16-8-R3_2.txt.gz

fq6c:
16-C_R1_1.txt.gz  16-C-R2_1.txt.gz  16-C-R3_1.txt.gz
16-C_R1_2.txt.gz  16-C-R2_2.txt.gz  16-C-R3_2.txt.gz


=cut
