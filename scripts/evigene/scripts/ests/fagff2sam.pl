#!/usr/bin/env perl
# fagff2sam.pl

use strict;
use warnings;
use Getopt::Long ;

my($ok, $fasta, $gff, $skiptypes, $outform, $debug)= (0)x10;
my $dropdit=1; # dang id mess \.1 at end on some
$outform="sam";

my $optok= &GetOptions (
  "gff=s"=>\$gff,  
  "fasta|sam=s"=>\$fasta,
  "outformat=s" => \$outform, # for .fasta instead of .sam, but location added for splitting.
  "skiptypes=s"=>\$skiptypes,
  "n|debug!"=>\$debug, 
  );

$fasta= shift @ARGV unless($fasta);
$gff  = shift @ARGV unless($gff);
if($outform =~ /sampe/i) {
  
}

die "USAGE: $0 -fasta my.fasta[.gz] -gff stdin|my.gff[.gz] > my.sam\n"
  unless($optok and (($fasta and $gff) or ($outform =~ /sampe/i)));

$skiptypes =~ s/[,\s]+/\|/g;

#warn: cant do both fasta gff on stdin
my $fain;
if($fasta =~ /\.gz/){ $ok= open(FA,"gunzip -c $fasta |"); $fain=*FA;}
elsif($fasta =~ /stdin|-/){ $ok= 1; $fain=*STDIN; }
else { $ok= open(FA,$fasta); $fain=*FA; }
die $fasta unless $ok;

if($outform =~ /sampe/i) {
  makepe($fain);  
  close($fain); exit(0);
}

my $gfh;
if($gff =~ /\.gz/){ $ok= open(GF,"gunzip -c $gff |"); $gfh=*GF; }
elsif($gff =~ /stdin|-/){ $ok= 1; $gfh=*STDIN; }
else { $ok= open(GF,$gff); $gfh=*GF; }
die $gff unless $ok;

# dang .fa has ID=xxx.1 ; .gff has ID=xxx
my($id, %fa, %fah);
while(<$fain>){ chomp; 
 if(/^>(\S+)/) { $id=$1; $id=~s/\.\d$// if($dropdit);
   if($outform =~ /fa/) { 
    # drop some header junk, e.g. : rank=0001767 x=2613.0 y=1202.0 length=222 trimmed=11-204
    s/^>\S+\s*//; s/\b(trimmed|length|len|size|x|y|rank)=[^\s;,]+[ ;,]?//g;
    $fah{$id}= $_; }
 } else { $fa{$id}.=$_; } }
close($fain);

# handle paired seqs? if sameid, but for  some id suffix like .fwd .rev

# sam line
my($qid, $flag, $chr, $cloc, $mapq, $cigar, $matechr, $mateloc, $isize, $seq, $qual, @opt);
    $mapq=255;
    $isize=0; #? pair insert size
    $qual="*"; @opt=();
    $matechr="*"; $mateloc=0; # unless have mate info.. see paired below
my ($cend, $orient, $len);
my(%did);

while(<$gfh>){
  next if($skiptypes and m/\t($skiptypes)\t/);
  if(/\bID=([^;\s]+)/) { 
    $qid=$1;  $qid=~s/\.\d$// if($dropdit);
    next if($did{$qid}++);
    $seq= $fa{$qid} or next;
    my @gff=split"\t";
    ($chr,$cloc,$cend,$orient)= @gff[0,3,4,6];

    ## what of chimeric aligns? skip?
    next if(/;chimera=/); # or path=\d+.\d+
# ID=inb1.FGXZ4J201ETFXN_C1;Target=inb1.FGXZ4J201ETFXN 1 137;aalen=46;cov=59.8;match=136;nexon=1;pid=99.3;qlen=229;path=1/2;chimera=exon-exon boundary at 137;chim2=scaffold00655:2151-2242:.
# ID=inb1.FGXZ4J201ETFXN_C2;Target=inb1.FGXZ4J201ETFXN 138 229;aalen=7;cov=40.2;match=88;nexon=1;pid=95.7;qlen=229;path=2/2;chimera=exon-exon boundary at 137;chim1=scaffold01348:468157-468293:.;cdsindel=26

    ## add mismatch info from gmap.gff annot: @opt NM:i:2 ?? max(qlen-match, indels) ?
# NM:i:1 for Target=xxx 1 62;aalen=21;cov=100.0;indels=0/1;match=61;nexon=1;pid=96.8;qlen=62;cdsindel=-1
# NM:i:5 for targ=1CGYXY 1 176;aalen=35;cov=98.3;indels=0/5;match=172;nexon=1;pid=95.0;qlen=179;cdsindel=-5
    my($mat,$qlen,$mis)=(0,0,0);
    if( ($mat)=m/;match=(\d+)/ ) { ($qlen)=m/;qlen=(\d+)/; $mis=($qlen && $mat)?$qlen-$mat:0; }

    $flag=0;   #1=pair, 2=pairok, 4=mismatch 8=mismate, 16=rev, 32=revmate, 64=firstmate, 128=secondmate
    # $flag |= 16 if($orient eq "-"); # is this right? rev may be only for mates
    $len= length($seq);
    if($outform =~ /fa/) {
      my $loc="$chr:$cloc-$cend:$orient";
      my $fah= $fah{$qid} || "";
      # $seq =~ s/(.{60})/$1\n/g; 
      print ">$qid loc=$loc; len=$len; $fah\n",$seq,"\n";
    } else {
    $cigar= $len."M"; # add introns ? need gff match_part/exon/...
    @opt=();
    push(@opt,"XS:A:".$orient) if($orient =~ m/[+\-]/); # strand tag
    push(@opt,"XE:i:".$cend); # keep end point of align : should also reconstruct cigar from exons
    push(@opt,"NM:i:".$mis); # always?
    # ... if($paired and $qid =~ m/$pairpatt/) { }
    print join("\t",($qid, $flag, $chr, $cloc, $mapq, $cigar, $matechr, $mateloc, $isize, $seq, $qual, @opt)),"\n";
    }
  } 
}

sub putpe {
  my($fs,$rs)= @_;
  unless($rs) {  print $fs,"\n"; return; } 
  my @f=split"\t",$fs; my @r=split"\t",$rs; 
  
  my($fc,$fb)=@f[2,3]; my($rc,$rb)=@r[2,3]; 
  my($fe)= $fs=~m/XE:i:(\d+)/; my($re)= $rs=~m/XE:i:(\d+)/; 
  my($orf)= $fs=~/\t(XS:A:.)/; my($orr)= $rs=~/\t(XS:A:.)/;
  my $in= ($fb>$rb) ? $rb - $fe : $re - $fb; # can be -neg
  $f[8]=$in; $r[8]=-$in;
  if($fc eq $rc) { $fc=$rc="="; } 
  @f[6,7]=($rc,$rb); @r[6,7]=($fc,$fb); 
     if($f[0]=~/fwd$/) { $f[1] |= 0x0063; $r[1] |= 0x0093;  } 
  elsif($f[0]=~/rev$/) { $f[1] |= 0x0093; $r[1] |= 0x0063;  } 
  $f[0]=~s/\.(fwd|rev)//; $r[0]=~s/\.(fwd|rev)//; 
  push @f,$orr if($orr and not $orf); push @r,$orf if($orf and not $orr); 
  print join("\t",@f)."\n"; print join("\t",@r)."\n"; 
}
  
sub makepe
{
  my($inh)= @_;
  my($lg,$fs,$rs)= ("") x 3;
  while(<$inh>) {
    chomp; my($d,$f,$c,$cb,$mq,$cg,$mc,$mb,$in,$sq,$ql,@opt)=split"\t";
    (my $g=$d)=~s/\.(fwd|rev)//; 
    if($g eq $lg) { $rs=$_; } else { putpe($fs,$rs) if($fs); $rs=""; $fs=$_; } 
    $lg=$g; 
    }
  putpe($fs,$rs);
}

__END__

# x.8244810       16      Scaffold2       867     255     36M     *       0       0       CTAGACACCAAAAAAATAGCAGCGGCTATATGATTT      *
# x.1350890       0       Scaffold2       868     255     15M1I20M        *       0       0       TAGACACCAAAAAAAATAGCAGCGGCTATATGATTT      *

=item  est.sampe

sort daphpulex_estjgi05.sam | perl -ne'($d,$f,$c,$cb,$mq,$cg,$mc,$mb,$in,$sq,$ql,@opt)=split"\t"; (
$g=$d)=~s/\.(fwd|rev)//; if($g eq $lg) { $rs=$_; } else { putp() if($fs); $rs=""; $fs=$_; } $lg=$g; sub putp{ u
nless($rs) {  print $fs; return; } ($orf)= $fs=~/\t(XS:A:.)/; ($orr)= $rs=~/\t(XS:A:.)/; chomp($fs); chomp($rs)
;  @f=split"\t",$fs;  @r=split"\t",$rs; ($fc,$fb)=@f[2,3]; ($rc,$rb)=@r[2,3]; ($fe)= $fs=~m/XE:i:(\d+)/; ($re)=
 $rs=~m/XE:i:(\d+)/; if($fc eq $rc) { $fc=$rc="="; } @f[6,7]=($rc,$rb); @r[6,7]=($fc,$fb); if($f[3]>$r[3]) { $in=
$r[3]-$fe; } else { $in=$re - $f[3];} if($f[0]=~/fwd$/) { $f[1] |= 0x0063; $r[1]|= 0x0093; $f[8]=$in; $r[8]=-$i
n; } elsif($f[0]=~/rev$/){ $f[1] |= 0x0093; $r[1] |= 0x0063; $f[8]=-$in; $r[8]=$in; } $r[0]=~s/\.(fwd|rev)//; $
f[0]=~s/\.(fwd|rev)//; push @f,$orr if($orr and not $orf); push @r,$orf if($orf and not $orr); print join("\t",
@f)."\n"; print join("\t",@r)."\n"; }' > daphpulex_estjgi05.sam.pe


=cut