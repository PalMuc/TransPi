#!/usr/bin/env perl
# samstats.pl

=item usage

  samtools view my.bam | $evigene/scripts/rnaseq/samstats.pl
  samtools view my.bam | env keepids=0 table=1 name=myrna nreads=123456000 TMPDIR=/var/tmp \
    $evigene/scripts/rnaseq/samstats.pl

=cut

# melon2.% head  log.bs1
# Out of memory!
# Out of memory!
# Out of memory!
#... problem from saving all those %ids
#... given sorted sam, could drop %ids each change in chr, summing  id stats over chr, but
#    that means counting dups.


use strict;
use Getopt::Long;

# calc medians for len_ ; use median of medians to avoid huge mem
use constant ID_SIZE  => 1000000;
use constant MED_SIZE =>   50000;
use constant IDFILE => 1;

my $MIN_IDENT= 0.95; # or 90?  # 2err/37bp = 95%
my $KEEPIDS = (defined $ENV{keepids}) ? $ENV{keepids} : 0;
my $TABLEFORM= $ENV{table} || 0;
my $NAME= $ENV{name} || 0;
my $CHRLEN= $ENV{chrlen} || 0;
my $ALLOW_SOFTCLIP = (defined $ENV{softclip}) ? $ENV{softclip} : 1;

# FIXME: input read total for cases where sam lacks unmapped reads...
my $IN_NREADS= $ENV{nreads} || 0;

my(%ids, %oneid, @lastids);

## replace this list with hash of stats
my($n_in, $n_read, $n_aln, $n_aln_mate1, $n_pair, $n_pairok, $n_pairnear, # $n_pairfar, 
   $n_strand, $n_map0, $n_map1, $n_mapn, $n_loident,
   $n_delete, $n_insert, $n_mapread,
   $n_mate1, $n_mate2, $n_secondary, $n_duplicate, $n_softclip, $len_softclip,
   $n_idbreaks, $n_intron, $len_intron,  $len_mate, $len_insert, $len_read, $len_chr) = (0) x 50;

# my $optok= GetOptions(
#   "sam|hits|input=s", \$input, 
#   "identity|MIN_IDENT=i", \$MIN_IDENT, 
#   "source=s", \$SOURCE, 
#   "sorted!", \$insorted, 
#   "seqsplice!", \$seqsplice, #? format=
#   "softclipok!", \$ALLOW_SOFTCLIP, # for MIN_IDENT error skips
#   "debug!", \$debug, 
#   #"test=s", \@testin, # chr:b-e ?
#   );
#
# die "opt error" unless($optok);

my(%medlen, %med2len, @MEDIANS);
BEGIN{
@MEDIANS= qw(len_intron len_mate len_insert len_read len_softclip);
push( @MEDIANS, 'len_chr') if($CHRLEN);
foreach my $name ( @MEDIANS ) {  $med2len{$name}= $medlen{$name}= []; }
}

if($ENV{grouptable}) { grouptable(); exit; }
if($ENV{mergetable}) { mergetables(); exit; }
# otherwise fall thru to MAIN below ..

sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }

sub minMedianMax
{
  my($aref)= @_;
  my @ma = sort{ $a <=> $b } @$aref;
  my $mid= int ( scalar(@ma) / 2);
  return @ma[0, $mid, $#ma];
}

# FIXME: save total count of items with medlen{name}
sub median2median
{
  my($name, $last)=@_;
  # push( @{ $medlen{$name} }, $len);
  # median2median($name) if($counter % MED_SIZE == 0);
  return unless(ref $medlen{$name});
  
  # if($last or @{$medlen{$name}} >= MED_SIZE) 
  if(1) {
    if(@{$medlen{$name}} > 1) {
    my ($min,$med,$max)= minMedianMax( $medlen{$name} );
    $medlen{$name} = [];
if(1) {
  if( @{ $med2len{$name} } < 3 ) {
    @{ $med2len{$name} }= ($min, $med, $max);
  } else {
    my $min1= @{ $med2len{$name} }[0];
    my $max1= pop(@{ $med2len{$name} });
    @{ $med2len{$name} }[0]= _min( $min, $min1);
    push( @{ $med2len{$name} }, $med, _max( $max, $max1));
  }
} else {    
    # this not right? get too many min,max; not enough med in med2len ?
    unshift( @{ $med2len{$name} }, $min);
    push( @{ $med2len{$name} }, $med);
    push( @{ $med2len{$name} }, $max);
}    
    }
 
    if( $last or @{ $med2len{$name} } >= MED_SIZE) {
      my @aful= @{ $med2len{$name} };      
      @aful= sort{ $a <=> $b } @aful ;
      my $min= shift @aful;
      my $max= pop @aful;
      my $n= int(MED_SIZE / 2);
      $n= @aful if($n > @aful);      
      my $lo=int($n/2); 
      my @acut= splice( @aful, $lo, $n);
      unshift(@acut, $min); push(@acut, $max);
      $med2len{$name} = \@acut;
    }
  }
  
  #return  minMedianMax( $med2len{$name} );
}


#?? add internal samtools call instead of this pipe for samtools view my.bam?
MAIN:
while(<>) {
  my ($qid, $flag, $chr, $cloc, $mapq, $cigar, $matechr, $mateloc, $isize, $seq, $qual, @opt)
      = split"\t";
  next unless($flag =~ /\d/ and /^\w/);    
  $n_in++;
  
  my $f_pair    = ($flag & 0x0001)?1:0;
  my $f_pairok  = ($flag & 0x0002)?1:0; #??
  my $f_mismatch= ($flag & 0x0004)?1:0; #  f_mismatch == $chr eq "*"
  my $f_mismate = ($flag & 0x0008)?1:0;
  my $f_isrev   = ($flag & 0x0010)?1:0;     # 16
  # my $f_revmate= ($flag & 0x0020); # 32
  my $f_first = ($flag & 0x0040)?1:0; # 64; mate ids
  my $f_second= ($flag & 0x0080)?1:0; # 128
  my $f_secondary= ($flag & 0x100)?1:0;
  my $f_duplicate= ($flag & 0x400)?1:0;

## revise pe id read/id counts to do both of pair: append /1 /2 to ids
#      p=0x1 (paired), P=0x2 (properly paired), u=0x4 (unmapped),
#      U=0x8 (mate unmapped), r=0x10 (reverse), R=0x20 (mate reverse)
#      1=0x40 (first), 2=0x80 (second), 's' = 0x100 (not primary), 
#      f=0x200 (failure) and d=0x400 (duplicate). 
     
  $n_pair++  if($f_pair);  # should this count if no align?
  $qid.="/2" if($f_second);
  
  if($f_mismatch) {
    $ids{$qid}= 0 unless(defined $ids{$qid});

  }  else { 
    $n_aln++;
    $n_aln_mate1++ if($f_pair and $f_first); # ? replace with $n_first $n_second ?
    $ids{$qid}++; 

    # my $lens  = length($seq); # change here, count from cigar
    my $cend = $cloc; # calc correct end from align/intron spans
    my(@aln,@intr,@itype,$instrand,$inmismat);
    my $lenc=0;
    my $alen=0; 
    my $softclip=0;
     
  # Op BAM Description : SAM Cigar ops, 2010
  # M 0 alignment match (can be a sequence match or mismatch) 
  # I 1 insertion to the reference 
  # D 2 deletion from the reference 
  # N 3 skipped region from the reference == intron
  # S 4 soft clipping (clipped sequences present in SEQ) 
  # H 5 hard clipping (clipped sequences NOT present in SEQ) 
  # P 6 padding (silent deletion from padded reference) 
  # = 7 sequence match    << M eqiv
  # X 8 sequence mismatch << M eqiv
  # ** H is treated odd cuz seq is clipped to it
  # ** S softclip also odd, $cloc is at NON-clipped seq start.. and cend = cloc + seqlen - end S
  # ** add stats for S, H other cig flags
  
  # 2012.7.2 : Insert bug for intron position  1-I shifts intron by +1

    if($cigar =~ m/^(\d+)M/ ) { $lenc= $alen=$1;  $cend += $1; } #  @aln=($1);
    ## elsif($cigar =~ /^(\d+)[HSN]/) {  } #$lenc+= $1; @aln=(0); need for intron with nSnM cigar
    while($cigar =~ m/(\d+)([DIHNSP])(\d+)M/g) { 
      my($bi,$bt,$bx)=($1,$2,$3);
      
      ## H/S/I mean read longer than genome at align; D/N/P other way
      if($bt eq 'N') { 
        $n_intron++; $len_intron += $bi;  $cend += $bi; # cend not changed for HISP; 
        push( @{ $medlen{'len_intron'} }, $bi);
      } elsif($bt eq 'H') { 
        $lenc += $bi; $bi=0; 
      } elsif($bt eq 'S') {
        $softclip += $bi; $lenc += $bi; $bi=0; 
      } elsif($bt eq 'I') {
        $n_insert++; $lenc += $bi; $bi=0; 
      } elsif($bt eq 'D') {  # P also?
        $n_delete++; $cend += $bi;
      } elsif($bt eq 'P') {  
        $cend += $bi;
      }
            
      $alen+= $bx; $lenc += $bx;  $cend += $bx;  
      # push(@intr,$bi);  push(@itype,$bt);  push(@aln,$bx); 
    }
    if($cigar =~ /(\d+)S$/) { $lenc+= $1; $softclip += $1; } # end
    elsif($cigar =~ /(\d+)H$/) { $lenc+= $1; } # end
  
    $len_read += $lenc;
    push( @{ $medlen{'len_read'} }, $lenc);
    push( @{ $medlen{'len_chr'} }, $cloc, $cend) if($CHRLEN);
    
    # my $nx= scalar(@aln);
    # my $alen=0; map{ $alen+=$_ }@aln; 
    my $amis= 0;
    foreach (@opt) { 
      if(m/NM:i:(\d+)/){ $amis += $1; }  # distinguish align-length and mismatches inside align
      elsif(m/XS:A:([\+\-])/){ $instrand=$1; $n_strand++;  } 
      # elsif(m/NS:i:(\d+)/) { $amis += $1; } # $inmismat=$1;  #new, if XS, = Mismatches within min_anchor_len of a splice junction, > 0 is poor
      # elsif(m/RG:A:(\S+)/){ $readgroup=$1; } #new
      }
   

    $alen= $alen-$amis; #? or not, report both?
    if($softclip>0) { push( @{ $medlen{'len_softclip'} }, $softclip); $len_softclip+=$softclip; $n_softclip++; }
    if($ALLOW_SOFTCLIP and $softclip>0) {
      $n_loident++ if( $lenc < 1 || $MIN_IDENT > $alen/($lenc-$softclip));  
    } else { 
      $n_loident++ if( $lenc < 1 || $MIN_IDENT > $alen/$lenc );
    }
    $n_secondary++ if($f_secondary);
    $n_duplicate++ if($f_duplicate);
 
    if($f_pair and not $f_mismate) {
      $n_pairok++ if($f_pairok);
      if($matechr eq "=") {
        $n_pairnear++;
        #Not inner?# # my $lenm= ($mateloc < $cloc) ? $cloc - ($mateloc+$lenc) : $mateloc - $cend; # inner
        my $lenm= ($mateloc <= $cloc) ? $cend - $mateloc : ($mateloc+$lenc) - $cloc; # outer
        # my $lenm= ($mateloc < $cloc) ? $cloc - $mateloc : $mateloc - $cloc; # mid
        $len_mate += $lenm;        
        push( @{ $medlen{'len_mate'} }, $lenm);
        if($f_pairok && $isize>0) {  $len_insert += $isize; push( @{ $medlen{'len_insert'} }, $isize); }  
      }
    }
  
  }

  ##  # maybe fix id breaks by cache last 10,000 ids for next count?
  # shift(@lastids) if(@lastids > 10000); push(@lastids, $qid);
  # unshift(@lastids, $qid);

  if($n_in % MED_SIZE == 0) { 
    foreach my $name ( @MEDIANS ) { median2median( $name ); }
    }

  if( $n_in % ID_SIZE == 0 ) {
    countIds(1);
    # my $nid=scalar(%ids); $nid =~ s,\d+/,,;
    # countIds(1) if($nid >= ID_SIZE);
  }
  
}


## ids hash is chewing too much memory; count & drop at intervals of ? 1Mil reads?
# counts from unique %ids

## redo using tmp files of ids, still too big for mem.
my($fid0, $fid1, $fid2, $fid01);
BEGIN{
  my $tmpdir= $ENV{TMPDIR} || "/var/tmp"; # system dependent : default /tmp ??
  my $nam = $ENV{samname} || "samst"; 
  my $fidp="$tmpdir/$nam.sst$$"; # use $ENV{TMPDIR}
  $fid0= "$fidp.id0"; # also .id0.s
  $fid1= "$fidp.id1"; # also .id1.s
  $fid2= "$fidp.id2"; 
  $fid01= "$fidp.idups";
}

END{
## FIXME: gone are $fid1.s $fid0.s ; only fid1, fid0
if($KEEPIDS) {
print "# Duplicate id file: $fid01\n";
print "# Single map id file: $fid1.s (includes dups, uniq == comm -23 $fid1.s $fid01 \n";
print "# No     map id file: $fid0.s (includes dups, uniq == comm -23 $fid0.s $fid01\n";
print "\n";
  # my $fidp="samstats.$$"; # caller name opt
  # rename("$fid0.s", "$fidp.id0");  # bad usage, across filesys; use mv ?
  # rename("$fid1.s", "$fidp.id1");
  # rename($fid01, "$fidp.idups");
} else {
  unlink("$fid0.s"); unlink("$fid1.s"); unlink($fid01);
}
  unlink($fid0); unlink($fid1); unlink($fid2); 
}

# version to keep id1 in ids is very slow. maybe above is calling this every 100k lines?
# keep in oneid hash instead;

sub countIds
{
  my($step)=@_;

if(IDFILE) {
    open(F0,">>$fid0"); open(F1,">>$fid1"); open(F2,">>$fid2");
    while ( my($id, $vid) = each %ids) {
       if($vid > 1) { print F2 $id,"\n"; }
       elsif($vid == 1) { print F1 $id,"\n"; }
       else { print F0 $id,"\n"; }
    }
    close(F0); close(F1); close(F2);
    %ids=(); undef %ids; $n_idbreaks++; 
}

  if($step) {
unless(IDFILE) {
    ## my %oneid=();
    while ( my($id, $vid) = each %ids) {
      if($vid == 1) {
        if( $oneid{$id} ) { $vid++; delete $oneid{$id}; }
        else { $oneid{$id}= 1; }
      } 
      
      if($vid != 1) { # count and drop     
        $n_read++;  
        if($n_pair > 0 and $id =~ m,/2$,) { $n_mate2++; } else { $n_mate1++; }
        if($vid == 0) { $n_map0++; }  
        elsif($vid > 1) { $n_mapn++; }
      }
    }
    undef %ids; $n_idbreaks++; # 
    ## %ids= %oneid;
}
  } else { # last

if(IDFILE) {
  ## ?? use $TMPDIR for sort, No Need; 
  # my $tmpdir= $ENV{TMPDIR} || "/var/tmp"; # system dependent : default /tmp ??
  # my $sortc=($tmpdir)? "sort -T $tmpdir" : "sort"; # No Need; ENV{TMPDIR} is what sort looks for
  
  my($map0,$map1,$map2)=(0,0,0);
  
## NOTE: add fid2 twice to make dups
  system("sort $fid0 $fid1 $fid2 $fid2 | uniq -d > $fid01");  # file of dup ids
  $map2 = `wc -l $fid01`; chomp($map2);
  if(1) { 
  $map0 = `fgrep -c -v -f $fid01 $fid0`; 
  $map1 = `fgrep -c -v -f $fid01 $fid1`;  # problems here for HUGE fid1
  } else { 
  system("sort -u $fid0 > $fid0.s");    #? not needed?
  $map0 = `fgrep -c -v -f $fid01 $fid0.s`; # $map0 = `comm -23 $fid0.s $fid01 | wc -l`; 
  system("sort -u $fid1 > $fid1.s");   #? is this sort -u needed? fgrep should be removing  dups
  $map1 = `fgrep -c -v -f $fid01 $fid1.s`; # $map1 = `comm -23 $fid1.s $fid01 | wc -l`;  # problems here for HUGE fid1
  }
  chomp($map0);
  chomp($map1); 

  # ? add hack correction if failed ; map1 == 0 ?
  my $map1c = $n_aln - ($map2 + $n_secondary + $n_duplicate); 
  if($map1 < 2 or $map1*2 < $map1c) { $map1 = $map1c; }
  
## 2012jun FIXME ***
## ** solaris largefile: comm is not able to process large files > 2GB; fgrep, others are ok
# solaris comm: /var/tmp/samst.sst16417.id1.s: Value too large for defined data type : at 5.1G
# is gnu-comm any better? ; result is map1 nreads=0; can we get by subtraction?
# .. IF counts of map0, mapn are low (fid0,fid01), could use other methods to remove from fid1
#  bamq5gs2/Dman_75_ATCACG_L006_1-dmag2.sstat2
#    fastq read counts: Dman_75_1 = 63636165; _2=63636165 == 127272330 total mates; too many more than nin
# .. this fastq set was updated with more reads since april *** need to remap
#  n=99940963 = n_aln + n_map0
#  n_aln   97440800    
#  n_map0  2500163
#  n_map1  0          0.0% of n_read
#    n_map1 = n_aln - n_mapn * num-mutli ? < dont have this count, nmulti >= 2
#      = 96585250 for multi=2; 96371362 for multi=2.5
#  n_mapn  427775
#  n_read       2927938 << also bad, from below sum
#  n_secondary   504613 << is this count of all dup alns (2ndary mappings)?  504613 2nd + 427775 = 932388 multi-aligns 
#     n_map1  = 97440800 - (504613 + 427775) = 96508412 ?
#     n_reads = 96508412 + 427775 + 2500163 = 99436350 .. about 10,000 too many, not enough multaligns above?
#  gzgrep -c '^\@HW' fastq4/Dman_75_ATCACG_L006_R1_001.fastq.gz = 49712961/1  + 49712961/2 = 99425922
#....................................

  $n_map0 = $map0;
  $n_map1 = $map1;
  $n_mapn = $map2;
  $n_read = $n_mapn + $n_map1 + $n_map0;
  ## use $IN_NREADS if given,  add diff to n_map0 unmapped
  if($IN_NREADS > $n_read) { $n_map0= $IN_NREADS - ($n_mapn + $n_map1); $n_read=$IN_NREADS;  }

} else {
    while ( my($id, $vid) = each %ids) {
      if($oneid{$id}) { $vid++; delete $oneid{$id}; }
      $n_read++;  
      if($n_pair > 0 and $id =~ m,/2$,) { $n_mate2++; } else { $n_mate1++; }
      if($vid == 1) { $n_map1++; }
      elsif($vid == 0) { $n_map0++; }
      elsif($vid > 1) { $n_mapn++; }
    }
    while ( my($id, $vid) = each %oneid) {
      next unless($id);
      $n_read++;  
      if($n_pair > 0 and $id =~ m,/2$,) { $n_mate2++; } else { $n_mate1++; }
      if($vid == 1) { $n_map1++; }
    }
# if($IN_NREADS > $n_read) { $n_map0= $IN_NREADS - ($n_mapn + $n_map1); $n_read=$IN_NREADS;  }
}
    undef %ids;  undef %oneid;
  }
}

# ** NOT Counting n_mapn right; need to keep all singles in hash til end.
# sub countIds_OLD
# {
#   my @savid= @_;
#   my @savn=();
#   if(@savid) { 
#     my %savid= map{ $_,1 } @savid;
#     @savid= sort keys %savid;
#     @savn = delete @ids{@savid}; 
#     }
# 
#   while ( my($id, $vid) = each %ids) {
#     next unless($id); #? undefs?
#     $n_read++;  
#     if($n_pair > 0 and $id =~ m,/2$,) { $n_mate2++; } else { $n_mate1++; }
#     if($vid == 1) { $n_map1++; }
#     elsif($vid == 0) { $n_map0++; }
#     elsif($vid > 1) { $n_mapn++; }
#   }
#   undef %ids; $n_idbreaks++; # 
#   if(@savid) { @ids{@savid}= @savn; }
# }

countIds(0);
foreach my $name ( @MEDIANS ) { median2median( $name, 1 ); }

$n_mapread= ($n_read - $n_map0) || 1;
$len_read= ($n_aln<1) ? 0: sprintf "%.1f", $len_read / $n_aln;
$len_intron= ($n_intron<1) ? 0: sprintf "%.0f", $len_intron / $n_intron;
$len_mate= ($n_pairnear<1) ? 0: sprintf "%.0f", $len_mate / $n_pairnear;
$len_insert= ($n_pairnear<1) ? 0: sprintf "%.0f", $len_insert / $n_pairnear;
$len_softclip= ($n_softclip<1) ? 0: sprintf "%.0f", $len_softclip / $n_softclip;

putout();

#................


sub putout
{
  my $n_map1c = $n_aln - ($n_mapn + $n_secondary + $n_duplicate); # print for check?
  
  my @colvar= qw(n_in len_read n_read n_aln n_strand n_intron len_intron
        n_map0 n_map1 n_mapn n_pair n_pairok n_pairnear len_mate len_insert n_loident 
        n_softclip len_softclip n_delete n_insert
        n_secondary n_duplicate);
  
  push(@colvar,"len_chr") if($CHRLEN); # want min,max always here

if($TABLEFORM) {
  foreach my $dorow (1..3) { # what? was 1..2, no pct
  next if($dorow==1 and $TABLEFORM > 1);
  my $na=$NAME || "count";
  print (($dorow == 1) ? "column" : ($dorow==2) ? $na : "percent");
  # print "# SAMstats ; input n=$n_in\n";
  foreach my $k (@colvar) 
  {
    my $v= eval("\$".$k);
    my $p=0;
    if($k =~ /n_map/ and $n_read>0) { $p= $v / $n_read; }
    ## elsif($k =~ /n_aln/ and $n_pair>0) { $p= $v / $n_mapread; } # was $p= $n_aln_mate1 / $n_mapread;
    elsif($k =~ /n_aln/ and $n_mapread>0) {  $p= $v / $n_mapread; }
    elsif($k =~ /n_pairok|n_pairnear/ and $n_aln_mate1>0) { $p= $v / (2*$n_aln_mate1); } # was /$n_pair;
    elsif($k =~ /n_pair$/ and $n_in>0) { $p= $v / $n_in; } # was /$n_pair;
    elsif($k =~ /^len_/) { 
      my ($min,$med,$max)= minMedianMax( $med2len{$k} );
      $p=$med; # "$min,$med,$max" if($med);
      $p= $v= 1 + $max - $min if($k =~ /len_chr/); # want span
    }
    elsif($n_aln>0) { $p = $v / $n_aln; }
    
    if( $dorow == 1 ) { print "\t$k"; }
    elsif( $dorow == 2 ) { printf "\t%7d", $v; }
    elsif( $dorow == 3 ) {
    if($k=~/^len_/) { print "\t$p"; }
    else { printf "\t%5.1f%%", 100*$p; }
    }
  }  
  print "\n";
  }

} else {

  print "# SAMstats ; NAME=$NAME; input n=$n_in\n";
  foreach my $k (@colvar) 
  {
    my $v= eval("\$".$k);
    my $p=0;
    if($k =~ /n_map/ and $n_read>0) { $p= $v / $n_read; }
    ## elsif($k =~ /n_aln/ and $n_pair>0) { $p= $v / $n_mapread; } # was $p= $n_aln_mate1 / $n_mapread;
    elsif($k =~ /n_aln/ and $n_mapread>0) {  $p= $v / $n_mapread; }
    elsif($k =~ /n_pairok|n_pairnear/ and $n_aln_mate1>0) { $p= $v / (2*$n_aln_mate1); } # was /$n_pair;
    elsif($k =~ /n_pair$/ and $n_in>0) { $p= $v / $n_in; } # was /$n_pair;
    elsif($k =~ /^len_/) { 
      my ($min,$med,$max)= minMedianMax( $med2len{$k} );
      $p="$min,$med,$max" if($med);
      $v= 1 + $max - $min if($k =~ /len_chr/); # want span
    }
    elsif($n_aln>0) { $p = $v / $n_aln; }
    
    printf "%10s\t%8d",$k, $v;
    printf "\t%5.1f%%", 100*$p unless($k =~ /n_read|^len_/);
    if($k=~/^len_/) { print "\t$p  min,median,max" if $p; }
    elsif($k=~/n_aln/) { print " of mapped_read"; }
    elsif($k=~/n_map/) { print " of n_read";
      if($k eq "n_map1") { print "\tn_map1check=$n_map1c"; }
      }
    elsif($k=~/n_pairok|n_pairnear/) { print " of n_aln_pair"; }
    elsif($k=~/n_pair/) { print " of n_input"; }
    elsif($k=~/loident|_pair|n_strand|^n_intron|n_insert|n_delete|n_second|n_duplic/) { print " of n_aln"; }
    print "\n";
  }  
  print "\n";

}

} # putout



sub mergetables
{
  # long version above: SAMstats
  while(<>) {
    if(/^# SAMstats ; NAME=([^;]+)/) { $NAME=$1 unless($NAME); }
    elsif(/^\s*(\w+)\s+(\d+)\s*(\S*)(.*)$/) { 
      my($k,$v,$p,$comm)=($1,$2,$3,$4);
      if($k =~ /^len_/) { 
        eval("\$".$k." = $v;"); 
        my($min,$med,$max)=split",",$p;
        @{ $med2len{$k} }= ($min, $med, $med, $med, $max);
        }
      else { eval("\$".$k." += $v;"); } # if($k =~ /^n_/);

    }
  }

  # $n_read = $n_mapn + $n_map1 + $n_map0; #?? should get from input ! no p,comm
  $n_mapread= ($n_read - $n_map0) || 1;
  # $n_mapread = $n_map1 + $n_mapn;
  $n_aln_mate1= int($n_pair/2); #?

  putout(); #seems ok  
}


sub grouptable
{
  # drop n_strand for n_intron
  my %drops=(n_in=>1,n_duplicate=>1,n_strand=>1,n_pair=>1,n_pairnear=>1,len_mate=>1,
             n_softclip=>1,len_softclip=>1); # n_pairnear == pairok nearly same

  my($na,@na,%kv,%pv,%lv,%mv,%lb,@k,%kf);
  while(<>) {
    if(/=== (\S+) ===/){ $na=$1; $na=~s/\.(sam|bam)//; push(@na,$na); } 
    elsif(/NAME=([^\s;]+)/) { my $nb=$1; $nb=~s/\.(sam|bam)//; unless($nb eq $na) { $na=$nb; push(@na,$na); } } 
    elsif(/^(\s*(?:n|len)_\w+)\s+\d/) { 
      my $kf=$1; my($k,$v,$p)=split; next if($drops{$k}); 
      my $lb = (/\% of (.+)/) ? $1 : /(median)/ ? $1 : "n"; 
      if($k=~/len_/) { $lv{$k}{$na}=$v; $mv{$k}{$na}= $p; } 
      else { $kv{$k}{$na}=$v;  $p=~s/%//; $p=int($p) if($p>100); $pv{$k}{$na}=$p; } 
      push(@k,$k) unless($kf{$k}); $kf{$k}=$kf;  $lb{$k}= $lb;  
    } 
  }

sub putt{ my %kv=@_; 
  print "statistic"; foreach my $na (@na) { print "\t$na"; } print "\n"; 
  foreach my $k (@k) { my @val= grep /\w/, @{$kv{$k}}{@na}; next unless(@val); 
  my $kf= $kf{$k}; print $kf; foreach my $na (@na) { my $v=$kv{$k}{$na};
  $v =~ s/^\d+,//; $v=~s/,.*//;  print "\t",$v; } print "\n"; } 
}  

  print "\nCounts\n"; putt(%kv); print "\nPercents\n"; putt(%pv); 
  print "\nSizeAverage\n"; putt(%lv); print "\nSizeMedian\n"; putt(%mv); 
  printkey();
}

sub printkey
{
  print (('-') x 70);
  print "\n\n";
  print <<'EOK';
Counts:  
  n_aln    = alignments
  n_read   = reads input,        n_map1 = reads uniquely mapped, 
  n_map0   = reads not mapped,   n_mapn = reads multiply mapped
  n_intron = spliced alignments, n_pairok = properly mapped pairs 
  n_delete = deleted read bases, n_insert = inserted read bases
  n_loident = alignments below 95% identity
  n_secondary = second+ of multiple alignments
  len_read,_intron,_insert: span in bp for reads, introns, and pair inserts
Percents: 
  p_map0, map1, mapn = % of n_read; p_aln = % of mapped_read (map1+mapn);
  p_intron, loident, delete, insert = % of n_aln;
  p_pairok, pairnear = % of n_aln_pair;
EOK

#   n_pair  = number of paired alignments
#   n_pairfar = number of pairs on separate scaffolds
#  n_duplicate = number aligns flagged as duplicate

}

__END__

#...........................................

=item key

Key:  
  n_read = number of reads input, 
  n_map0 = reads not mappable, 
  n_map1 = reads w/ single map location, 
  n_mapn = reads multiply mapped.
  n_aln  = number of alignments (> n_map1 + n_mapn due to n_mapn * x mappings)
  n_intron = number of spliced alignments
  n_pair  = number of paired alignments
  n_pairok =  " with both mates mapped
  n_pairfar = number of pairs on separate scaffolds
  n_loident = number of alignments with < 95% identity (e.g. 6+ mismatches for 109 bp read)
  n_secondary = number aligns flagged as second+ of multiple alignments
  n_duplicate = number aligns flagged as duplicate
  len_read,_intron,_mate: average span in bp (not valuable due to large range in spans)
  
 Percent of: 
 p_aln = % of mapped_read, p_strand = % of n_aln, p_intron = % of
 n_aln, p_map0 = % of n_read, p_map1 = % of n_read, p_mapn = % of
 n_read, p_pair = % of n_aln, p_pairok = % of n_pair, p_pairfar =
 % of n_pair, p_loident = % of n_aln,
 Percents can be > 100 due to choice of denominator, e.g.
 %p_aln is all alignments of n_mapn reads + single n_map1 reads / (n_mapn+n_map1).
 %p_pair is number of alignments + unmapped reads that have a mate pair / n_alignments,
      better percent would be n_pair / (n_aln + n_map0)

=item tabulate the long tables

cat 4-*.bam.sstats | perl -ne \
'if(/=== (\S+) ===/){ $na=$1; $na=~s/\.(sam|bam)//; push(@na,$na); } \
elsif(/NAME=([^\s;]+)/) { $nb=$1; $nb=~s/\.(sam|bam)//; unless($nb eq $na) { $na=$nb; push(@na,$na); } } \
elsif(/^(\s*(?:n|len)_\w+)\s+\d/) { $kf=$1; ($k,$v,$p)=split; next if($drops{$k}); \
$lb = (/\% of (.+)/) ? $1 : /(median)/ ? $1 : "n"; \
if($k=~/len_/) { $lv{$k}{$na}=$v; $mv{$k}{$na}= $p; } \
else { $kv{$k}{$na}=$v;  $p=~s/%//; $p=int($p) if($p>100); ($kp=$k)=~s/n_/p_/; $pv{$k}{$na}=$p; } \
push(@k,$k) unless($kf{$k}); $kf{$k}=$kf;  $lb{$k}= $lb;  } \
sub putt{ my %kv=@_; print "statistic"; foreach my $na (@na) { print "\t$na"; } print "\n"; \
foreach my $k (@k) { my @val= grep /\w/, @{$kv{$k}}{@na}; next unless(@val); \
my $kf= $kf{$k}; print $kf; foreach $na (@na) { $v=$kv{$k}{$na};\
$v =~ s/^\d+,//; $v=~s/,.*//;  print "\t",$v; } print "\n"; } }  \
END{ print "\nCounts\n"; putt(%kv); print "\nPercents\n"; putt(%pv); \
print "\nSizeAverage\n"; putt(%lv); print "\nSizeMedian\n"; putt(%mv); }\
BEGIN{ %drops=(n_in=>1,n_duplicate=>1,n_intron=>1,n_pairnear=>1,n_pair=>1,len_mate=>1,n_softclip=>1,len_softclip=>1); } ' \

> samstats.4.table

# drop these: n_in?; n_duplicate; n_intron (= n_strand); n_pairnear, n_pair; 
#  len_mate (for len_insert); n_softclip,len_softclip

# NAME=10_1-dmag2
cat bamq5gs3/stats/*.sstat | sed 's/_1-dmag2//; s/NAME=/NAME=nd/;' | perl -ne\

cat bamq5gs3/stats/*.sstat | sed 's/_1-dmag2//; s/NAME=/NAME=nd/;' | env grouptable=1 $evigene/scripts/rnaseq/samstats.pl

# NAME=Dman_03_CGATGT_L005_1-dmag2
cat bamq5gs/stats/*.sstat | sed 's/_1-dmag2//; s/NAME=Dman_/NAME=hi/; s/_[ACGT][ACGT][ACGT]*_.*//;' |\
 env grouptable=1 $evigene/scripts/rnaseq/samstats.pl

=item groups table

Counts
statistic       4-1R1   4-1R2   4-1R3   4-8R1   4-8R2   4-8R3   4-CR1   4-CR2   4-CR3
  len_read      51      50      51      51      51      51      51      50      50
    n_read      52865278        36691656        39965592        84308104        26683512        55037494        40369626        50430490        52514126
     n_aln      50347239        35152218        32564974        79794602        24727492        52457282        37563535        46108452        49429718
  n_strand      2139300 1429336 912311  3204466 795816  2141378 1468780 1750909 1831069
  n_intron      2139300 1429336 912311  3204466 795816  2141378 1468780 1750909 1831069
len_intron      182     203     288     182     243     180     177     180     191
    n_map0      3956716 2541447 9640520 7050661 2961902 4087683 4056344 5765101 4629694


=item samples

   head -40 daphmag3/bamq5gs/stats/Dman_03*sstat
  # SAMstats ; NAME=Dman_03_CGATGT_L005_1-dmag2; input n=127607078
        n_in      127607078       103.9%
    len_read           101        97,101,101  min,median,max
      n_read      126900715
       n_aln      122800321       100.5% of mapped_read
    n_strand      27554711         22.4% of n_aln
    n_intron      28348478         23.1% of n_aln
  len_intron           209        1,83,199048  min,median,max
      n_map0       4770769          3.8% of n_read
      n_map1      121597861        95.8% of n_read        n_map1check=121673020
      n_mapn        532085          0.4% of n_read
      n_pair      127607078       100.0% of n_input
    n_pairok      117203147        94.3% of n_aln_pair
  n_pairnear      117192906        94.3% of n_aln_pair
    len_mate           557        25,319,3703188  min,median,max
  len_insert           234        1,324,3503723  min,median,max
   n_loident       5471295          4.5% of n_aln
  n_softclip      28351072         23.1%
  len_softclip          18        1,37,67  min,median,max
    n_delete       2162674          1.8% of n_aln
    n_insert       2879080          2.3% of n_aln
  n_secondary       595216          0.5% of n_aln
  n_duplicate            0          0.0% of n_aln

  head -40 $dmag/rnas/bamm/3-1R1.bam.sstats : older gsnap runs
  #==== 3-1R1.bam =======================
  # SAMstats ; input n=14141039
    len_read            49        1,50,50  min,median,max
      n_read      13817836
       n_aln      12267151        102.7% of mapped_read
    n_strand        473822          3.9% of n_aln
    n_intron        473822          3.9% of n_aln
  len_intron           182        4,79,188030  min,median,max
      n_map0       1873888         13.6% of n_read
      n_map1      11695031         84.6% of n_read
      n_mapn        248917          1.8% of n_read
      n_pair      14141039        100.0% of n_input
    n_pairok       4784730         38.9% of n_aln_pair
  n_pairnear      11389355         92.7% of n_aln_pair
    len_mate           429        16,91,3309042  min,median,max
   n_loident       2123396         17.3% of n_aln
  n_secondary       323203          2.6% of n_aln
  n_duplicate            0          0.0% of n_aln


=cut
