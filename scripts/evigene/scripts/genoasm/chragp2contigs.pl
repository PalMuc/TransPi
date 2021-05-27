#!/usr/bin/env perl
# chragp2contigs.pl

use strict;
use warnings;
use Getopt::Long ;

my $GAPSIZE= 100; ## unknown_gap_size
my $FAWIDTH= 60; 
my $IDPREFIX= 'xxx'; 
#o# my $SCAFTAG= 'SG'; # SG{number} is always part of scafid, or see UG ungrouped set?
my $DIEonERR= 1; # debug opt
my($chrdir,$agpdir,$contigs,$debug)= (0) x 10;

my $optok= &GetOptions (
  "chrdir=s" => \$chrdir, # input chr dir OR all.fasta
  "agpdir=s" => \$agpdir, # input agp dir OR all.agp
  "contigs|ctgdir=s" => \$contigs, # output contigs folder
  # "IDPREFIX=s"=>\$IDPREFIX,  # not used?
  "debug!"=>\$debug, 
  );

# for cfa in seq/mm_ref_GRCm38.p3_*.fa.gz; do { .. }
$chrdir= shift @ARGV unless($chrdir); # or list of fa files?
$agpdir= shift @ARGV unless($agpdir); # or list of agp files?

$contigs="contigs" unless($contigs);
mkdir($contigs) unless(-e $contigs);

my(@agp,@chr,%agp,%agpc,%did);
if( -d $chrdir) {
  opendir(D,$chrdir); @chr= sort grep /\.fa/, map{ chomp; $_; } readdir(D);  closedir(D); 
} elsif( -f $chrdir ) {
  @chr=($chrdir); $chrdir=".";
}
if( -d $agpdir) {
  opendir(D,$agpdir); @agp= sort grep /\.agp/, map{ chomp; $_; } readdir(D);  closedir(D); 
} elsif( -f $agpdir ) {
  @agp=($agpdir); $agpdir=".";
}
 
for my $agpin (@agp) {
  my($pt,$chrin,$ino,$ok,$nid,$id,$fa,$infa,$nout)=(0) x 19;
  ($pt= $agpin)=~s/\.agp.*//; 
  ($chrin)= grep/$pt\.fa/,@chr; # fixme for not same name?
  unless($chrin) { if(@agp==1 and @chr==1) { $chrin=$chr[0]; } }
  
  warn "# extract contigs: $agpdir/$agpin, $chrdir/$chrin to $contigs/$pt.ctg.fa\n" if($debug);
  die "ERR:missing $agpdir/$agpin" unless(-f "$agpdir/$agpin");
  die "ERR:missing $chrdir/$chrin" unless(-f "$chrdir/$chrin");
    
  $ino= ($agpin=~/\.gz/) ? "gunzip -c $agpdir/$agpin|" : "$agpdir/$agpin";
  $ok= open(AGP,$ino) or die "open $ino";
  %agp=%agpc=();
  while(<AGP>) {
    #if(/^#/){  next;} #$ag=1; $infa=0; 
    if(/^\w\S+\t\d/) { 
      my($cid,$cb,$ce,$ci,$FN,$sid,$sb,$se,$so)=split; 
      next if($FN =~ m/^[NU]/ or not($sid=~/\w/ and $se>0));
      push @{$agp{$cid}},[$cb,$ce,$ci,$sid,$sb,$se,$so];
      push @{$agpc{$sid}},[$cb,$ce,$ci,$sid,$sb,$se,$so];# for problem cases?
    }
  } close(AGP); 

## problem case:
# agp: chrMT   1       16596   1       O       NC_002333.2     1       16596   +
# seq: >gixxxx|NC_002333.2 << ctgid not chrid of agp

  $ino= ($chrin=~/\.gz/) ? "gunzip -c $chrdir/$chrin |" : "$chrdir/$chrin";
  $ok= open(CHR,$ino) or die "open $ino";
  $ok= open(OUT,">$contigs/$pt.ctg.fa");
  while(<CHR>) {
    if(/^>/) { 
      $nout+= puts($id,$nid,$fa) if($id and $fa); 
      ($id)= m/^>(\S+)/;
      $id=~s/gi.\d+.ref.//; $id=~s/\|$//; 
      my($cn)=m/chromosome (\w+)/; 
      $nid=$id; $id="chr$cn" if($cn); 
      $infa=1;  $fa=""; 
    } elsif($infa) { chomp; $fa.=$_; } 
  } close(CHR); 
  $nout+= puts($id,$nid,$fa) if($fa); 
  close(OUT);
  warn "# extract contigs: n=$nout to $contigs/$pt.ctg.fa\n" if($debug);
    
}
# END

#--------
sub revc{my $s=reverse($_[0]); $s=~tr/ACGTacgtRYrySWswKMkmBDHVbdhv/TGCAtgcaYRyrWSwsMKmkVHDBvhdb/; $s;}

sub puts { 
  my($id,$nid,$fa)=@_; my $nout=0;
  my $agp=$agp{$id}||$agp{$nid}||$agpc{$nid};  
  ##? unless($agp) { $agp=$agpc{$nid}; }
  return 0 unless($agp); 
  for my $ag (@$agp) { 
    my($cb,$ce,$ci,$sid,$sb,$se,$so)=@$ag; 
    my $s=substr($fa,$cb-1,1+$ce-$cb); $s=revc($s) if($so eq "-"); 
    my $slen=length($s);
    #x#if($did{$sid}){ my $sp="a"; for $sp (qw(a b c d e f g h i j)){ last unless($did{$sid.$sp});} $sid.=$sp; }
    if($did{$sid}){ for my $sp (qw(a b c d e f g h i j)){ unless($did{$sid.$sp}){ $sid.=$sp; last; } } }
    $s=~s/(.{60})/$1\n/g; 
    print OUT ">$sid ofs=$sb-$se:$so; chrloc=$id:$cb-$ce; chrid=$nid; ci=$ci; length=$slen\n";
    print OUT "$s\n"; $did{$sid}++; $nout++;
  } 
  return $nout;
}


__END__

#!/bin/bash
# makectg.sh
# FIXME: some contigs split in agp,chrasm .. need to join, at least give new id: AC087780.19, AC087780.19b
# any splits across chr are more of a problem here.
# FIXME2: need contig revcomp for $so eq "-"

for cfa in seq/mm_ref_GRCm38.p3_*.fa.gz; do {
  pt=`basename $cfa .fa.gz`
  echo make contigs/$pt.ctg.fa
 
  gunzip -c agp/$pt.agp.gz seq/$pt.fa.gz | perl -ne \
'if(/^#/){ $ag=1; $infa=0;  next;} elsif(/^\w\S+\t\d/) { 
($cid,$cb,$ce,$ci,$FN,$sid,$sb,$se,$so)=split; if($FN ne "N") { push @{$agp{$cid}},[$cb,$ce,$sid,$sb,$se,$so]; } }
elsif(/^>(\S+)/) { puts($id,$nid) if($id and $fa); $nid=$id=$1; $id=~s/gi.\d+.ref.//; $id=~s/\|$//; 
($cn)=m/chromosome (\w+)/; $nid=$id; $id="chr$cn" if($cn); $infa=1;  $fa=""; } elsif($infa) { chomp; $fa.=$_; } 
END{ puts($id,$nid) if($fa); } 
sub revc{my $s=reverse($_[0]); $s=~tr/ACGTacgtRYrySWswKMkmBDHVbdhv/TGCAtgcaYRyrWSwsMKmkVHDBvhdb/; $s;}
sub puts{my($id,$nid)=@_; my $agp=$agp{$id}||$agp{$nid};  return 0 unless($agp); for my $ag (@$agp) { 
my($cb,$ce,$sid,$sb,$se,$so)=@$ag; $s=substr($fa,$cb-1,1+$ce-$cb); $s=revc($s) if($so eq "-"); 
if($did{$sid}){ for $sp (qw(a b c d e f g h i j)){ last unless($did{$sid.$sp});} $sid.=$sp; }
$s=~s/(.{60})/$1\n/g; print ">$sid ofs=$sb-$se:$so chrloc=$id:$cb-$ce\n$s\n"; $did{$sid}++; }}' \
  > contigs/$pt.ctg.fa

} done

