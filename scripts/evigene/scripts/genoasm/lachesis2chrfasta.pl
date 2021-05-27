#!/usr/bin/env perl
# lachesis2chrfasta.pl

=item lachesis2chrfasta.pl usage


  $genoasm/lachesis2chrfasta.pl -id kfish2la2s_ -contigs ../funhe302scaf.fasta  -lachesisdir laout2 \
  -outfa kfish2la2superall.fasta
  
  * fixme: all -contigs should go into outfa, with unplaced/ungrouped id/annot if not in grouping
  
=cut

use strict;
use warnings;
use Getopt::Long ;

my $GAPSIZE= 100; ## unknown_gap_size
my $FAWIDTH= 60; 
my $IDPREFIX= 'LachesisScaffold'; 
my $SCAFTAG= 'SG'; # SG{number} is always part of scafid, or see UG ungrouped set?
my $DIEonERR= 1; # debug opt
my($infasta,$lachesisdir,$outfasta,$outagp,$debug)= (0) x 10;

my $optok= &GetOptions (
  "IDPREFIX=s"=>\$IDPREFIX,  
  "GAPSIZE=i"=>\$GAPSIZE,
  "contigsfasta|infasta=s" => \$infasta, # input.contig.fasta file name
  "lachesisdir=s" => \$lachesisdir, # output.assembly.fasta file name
  "outfasta=s" => \$outfasta, # output.assembly.fasta file name
  "outagp=s" => \$outagp, # output.assembly.agp file name; make this default? not ready
  "debug!"=>\$debug, 
  );

$infasta= shift @ARGV unless($infasta);
$lachesisdir= shift @ARGV unless($lachesisdir);

die "usage: lachesis2chrfasta.pl contigs.fasta lachesisoutputdir > lascaffolds.fasta
opts: -contigs contigs.fasta -lachesisdir outdir 
  -outfasta lascaffolds.fasta -outagp lascaffolds.agp
  -gapsize=$GAPSIZE -idprefix=$IDPREFIX\n"
unless($optok and -s $infasta and -d $lachesisdir);

$DIEonERR=0 if($debug);
$outfasta = "$lachesisdir/${IDPREFIX}asm.fasta" unless($outfasta); ## add or $outagp
$outagp   = "$lachesisdir/${IDPREFIX}asm.agp" unless($outagp); ## added

my $clusters_file = "$lachesisdir/main_results/clusters.by_name.txt";
die "ERROR: Couldn't find clusters file ('clusters.by_name.txt') in $lachesisdir/main_results/. \n"
  unless -e $clusters_file;

my $N_orderings = 0;
while ( 1 ) {
  last unless -e "$lachesisdir/main_results/group$N_orderings.ordering";
  $N_orderings++;
}
die "ERROR: Couldn't find any ordering files ('group*.ordering') in $lachesisdir/main_results/. \n"
  unless ( $N_orderings );

my(%contigs_used, %facount); # %facount for simple size table and N50 summary
my($incontigids,$fasta_seqs)= readFasta($infasta);
# die "ERROR: xxx\n" unless(@$incontigids>0);

open OUT, '>', $outfasta or die "writing $outfasta";
open AGP, '>', $outagp or die "writing $outagp"; 
# add $outfasta.count table? or leave to AGP

my($nout,$nord,$nunord,$nunplaced)= (0) x 10;

foreach my $i ( 0..$N_orderings-1 ) { my($nscaf,$nctg)= scaffoldGroup($i); $nord+=$nctg; $nout+=$nscaf; }
warn "# nscaffolded=$nord\n" if($debug);

$nunord   = scaffoldUnorderedContigs($clusters_file);
$nout+= $nunord;
warn "# nunordered=$nunord\n" if($debug);

$nunplaced= unscaffoldedContigs($incontigids, $nout); ## add numbering consistant w/ above SG[1..nscaf]
$nout+= $nunplaced;
warn "# nunplaced=$nunplaced\n" if($debug);

my $n50val= basicN50($outfasta,\%facount);
warn "# $n50val\n";
print AGP "# $n50val\n";

close OUT;
close AGP;

warn "# Wrote $nout scaffolds to $outfasta, AGP to $outagp\n",
     "# $nord contigs ordered in scaffolds, $nunord unordered but grouped, $nunplaced unplaced contigs\n"
  ; #if $debug;

#------------ subs ---------------


=item basicN50 perlcmd

  -- needs facount or equiv here.. can count output bases ok.
  
  env nhalf=0 phalf=50 gaps=0 perl -ne \
'BEGIN{$PH=$ENV{phalf}||50; $NH=$ENV{nhalf}||0; $GAPS=$ENV{gaps}||0;}
if($ARGV ne $lfn) { putc($lfn); $lfn=$ARGV; }
($id,$w,$ca,$cc,$cg,$ct,$cn,$cpg)=split; $wn=($GAPS)?$w:$w-$cn;
if(/^\W/){} elsif(/^total/){} else { $wn{$id}=$wn; $swn+=$wn; $n++; }
END{ putc($lfn); } 
sub putc{ my($fn)=@_; return unless($fn and $n);
my $hwn=int($swn*$PH/100); $hwn=$NH if($NH>0);  
my @wn=sort{$wn{$b}<=>$wn{$a} or $a cmp $b}keys %wn;
my($id,$in,$wn,$wnt,@n5,@top); for $id (@wn) { $wn=$wn{$id}; $in++;
$wnt+=$wn; my @nv=($in,$id,$wn); push @top,\@nv; 
if($wnt>=$hwn and not @n5) { @n5=@nv unless(@n5); last if($in>$n5[0]+10); } }
$NN=($GAPS)?"+NNN":"-NNN";  $fn=~s/\.\w+$//; print "$fn\tn=$n,
totlen=$swn ($NN), n$PH=@n5 of sumlen${PH}=$hwn\n"; $in=$n=$swn=0; %wn=(); } ' \
    *.facount ../*.facount > kfish2la2superall.n50

=cut

# use AGP %facount table for this, not facount file
sub basicN50 {
  my($fn,$facountr)=@_;
  my $PH=$ENV{phalf}||50; my $NH=$ENV{nhalf}||0; my $GAPS=$ENV{gaps}||0;
  my($n,$swn,%wn)=(0,0);
  for my $id (sort{$facountr->{$b}[0] <=> $facountr->{$a}[0] or $a cmp $b} keys %{$facountr})
  {
    # my($id,$w,$ca,$cc,$cg,$ct,$cn,$cpg)=split; 
    my $w = $facountr->{$id}[0];
    my $cn= $facountr->{$id}[1];
    my $wn= ($GAPS)?$w:$w-$cn;
    $wn{$id}=$wn; $swn+=$wn; $n++; 
  }
  
  # putc($lfn);   
  my $hwn=int($swn*$PH/100); $hwn=$NH if($NH>0);  
  my @wn=sort{$wn{$b}<=>$wn{$a} or $a cmp $b}keys %wn;
  my($id,$in,$wn,$wnt,@n5,@top); 
  for $id (@wn) { 
    $wn=$wn{$id}; $in++;
    $wnt+=$wn; my @nv=($in,$id,$wn); push @top,\@nv; 
    if($wnt>=$hwn and not @n5) { @n5=@nv unless(@n5); last if($in>$n5[0]+10); } 
  }
  
  my $NN=($GAPS)?"+NNN":"-NNN"; $fn=~s/\.\w+$//; 
  return "N50:$fn\tn=$n, totlen=$swn ($NN), n$PH=@n5 of sumlen${PH}=$hwn"; 
}

sub MissingSeq {
  my($name,$where)=@_;
  my $err="ERROR: Contig '$name' referenced in $where, is not found in input $infasta\n";
  if($DIEonERR) { die $err; } else { warn $err; }
}

sub unscaffoldedContigs {
  my($contigids, $scafnum)=@_;
  $scafnum||=0;
  my $N_unclustered_contigs = 0;
  foreach my $ctgname ( @$contigids ) {
    next if $contigs_used{$ctgname};
    $contigs_used{$ctgname}++;
    MissingSeq($ctgname, "Clusters file $clusters_file") unless exists $fasta_seqs->{$ctgname};
    $scafnum++;
    my $scafname = $IDPREFIX.$SCAFTAG.$scafnum."ungroup_".$ctgname; 
    #which?# my $scafname = $IDPREFIX."UC${scafnum}_".$ctgname;  
    # my $scafinfo="contigs=1; length=$len;";
    my $faseq=$fasta_seqs->{$ctgname}; my $clen=length($faseq);
    putFasta($scafname, $faseq); # ,$scafinfo

    putAGP($scafname, [[$scafname,1,$clen,1,'W',$ctgname,1,$clen,'+']]);
    
    $N_unclustered_contigs++;
  } 
  return $N_unclustered_contigs;
}

sub scaffoldUnorderedContigs {
  my($clusters_file)=@_;
  # PHASE 2: For each contig that is clustered into a group but not ordered within that group, write it to file with a name indicating its group.
  open IN, '<', $clusters_file or die "Can't find file $clusters_file: $!";
  my $ID = 0;
  my $N_unordered_contigs = 0;
  while (<IN>) {
      # Each (non-commented) line simply contains a tab-delimited list of contig names in this cluster.  The cluster ID is implied by the line order.
      next if /^\#/; # skip commented lines
      my @names = split;
      
      foreach my $ctgname (@names) {
        next if $contigs_used{$ctgname};
        $contigs_used{$ctgname}++;
        MissingSeq($ctgname, "Clusters file $clusters_file") unless exists $fasta_seqs->{$ctgname};

        # *** use 1-origin IDS ***
        my $scafname = $IDPREFIX.$SCAFTAG.($ID+1)."rand_".$ctgname; # old style name for unordered = rand.
        # my $scafinfo="contigs=1; length=$len;"; #? add groupid as info
        my $faseq=$fasta_seqs->{$ctgname}; my $clen=length($faseq);
        putFasta($scafname, $faseq);
        
        putAGP($scafname, [[$scafname,1,$clen,1,'W',$ctgname,1,$clen,'+']]) ;

        $N_unordered_contigs++;
      }
      $ID++;
  }
  close IN;
  return $N_unordered_contigs;
}

sub scaffoldGroup {
  my($i)= @_;
  my $scaffold = '';
  my $N_contigs = 0;
  my($sb,$se,$si,$clen,$cor)=(0) x 5; my(@agp);
  my $ordering_file = "$lachesisdir/main_results/group$i.ordering";
  open IN, '< ', $ordering_file or die "Can't find file $ordering_file: $!";
  while (<IN>) {
    next if /^\#/; # skip commented lines
    my( $ID, $ctgname, $rc, $qual, $gap_size ) = split;
    $contigs_used{$ctgname}++;
    MissingSeq($ctgname, "Ordering file $ordering_file") unless exists $fasta_seqs->{$ctgname};

    my $faseq= $fasta_seqs->{$ctgname};
    $faseq = revcomp( $faseq ) if $rc;
    $clen= length($faseq);
    $cor= ($rc)?'-':'+';
    $scaffold .= $faseq;
    
    $gap_size = $GAPSIZE unless($gap_size =~ /\d/ and $gap_size>0); ## if $gap_size eq '.';
    my $gap = 'N' x $gap_size;
    $scaffold .= $gap;
    
    $sb=1+$se; $se=$sb+$clen; # sb=1 first
    push(@agp, [0,$sb,$se,++$si,'W',$ctgname,1,$clen,$cor]); ## W or ..
    $sb=1+$se; $se=$sb+$gap_size; 
    push(@agp, [0,$sb,$se,++$si,'N',$gap_size,'scaffold','yes','paired-ends']);  ## N or U type; evsource what?
    # push(@agp, contigrow); ## add this opt
    # push(@agp, gaprow);  # scaffold        yes     paired-ends
    
    $N_contigs++;
    }
    
  $scaffold =~ s/N+$//; # Remove the last gap.
  my $len = length $scaffold;
  return (0,$N_contigs) unless $N_contigs; # we may get an empty scaffold if we have a group in which no contigs were ordered (e.g., a singleton)
  
  #badname# my $scaffold_name = "Lachesis_group${i}__${N_contigs}_contigs__length_$len";
  my $scafname = $IDPREFIX.$SCAFTAG.($i+1); # *** use 1-origin IDS ***
  my $scafinfo="contigs=$N_contigs; length=$len;"; #? add groupid as info
  putFasta($scafname, $scaffold, $scafinfo);
  putAGP($scafname,\@agp);

  close IN;
  return (1,$N_contigs);
}

sub putFasta {
  my($scafname,$scaffold,$scafinfo)=@_;
  # my $scafname = $IDPREFIX."SG".$i;
  # my $scafinfo="contigs=$N_contigs; length=$len;";
  unless($scafinfo) { $scafinfo="length=".length($scaffold); }
  $scaffold =~ s/(.{$FAWIDTH})/$1\n/g; $scaffold=~s/\n$//;
  
  print OUT ">$scafname $scafinfo\n$scaffold\n";
}

sub putAGP {
  my($scafname,$agp)=@_;
  return 0 unless($outagp and ref($agp));
  
  # add: table scaf sizes, gaps for facount/n50 summary
  # contig cols: my($cid,$cb,$ce,$ci,$FN,$sid,$sb,$se,$so);
  # agp: chrMT   1       16596   1       O       NC_002333.2     1       16596   +
  # gap cols: my($cid,$cb,$ce,$ci,"N",gaplen,xxx
  my $nout=0;
  for my $row (@$agp) {
    my($cid,$cb,$ce,$ci,$FN,$sid,$sb,$se,$so)= @$row;
    unless($cid and $cid=~/^\w+/) { $$row[0]= $scafname; }
    print AGP join("\t",@$row)."\n"; $nout++;
    my $len=1+$ce-$cb;
    if($FN =~ /^[NU]/) { $facount{$scafname}[1]+= $len; } # gap count
    $facount{$scafname}[0] += $len;  # bases + gap count total, nongaps = total - gaps
  }
  return $nout;
}

sub readFasta {
  my($infa)= @_;
  my($nin,$cid,@contigids,%contigseq)=(0,0);
  open IN, $infa or die "reading $infa";
  while(<IN>) {
    if(/^>(\S+)/) { $cid=$1; push @contigids, $cid; $nin++; }
    else { chomp; $contigseq{$cid}.=$_; }
  } close(IN);
  warn "# read $nin contigs from $infa\n"; #if $debug;
  return(\@contigids,\%contigseq);
}

sub revcomp {
  my $seq = reverse $_[0];
  $seq =~ tr/ACGTacgtRYrySWswKMkmBDHVbdhv/TGCAtgcaYRyrWSwsMKmkVHDBvhdb/; #all IUPAC codes 
  return $seq;
}
