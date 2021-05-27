#!/usr/bin/env perl
# sam2seqparts.pl


=item about

  extract seq from sam into (a) genome partitions and (b) formsts for velvet, other
  see also partition_GFF_genome.pl
  see also partition_GFF_qnd.pl : combine this w/ that?
  need also new parts.list
  
  examples:

	sam2seqparts.pl -type bam_single -in rnaseq_sr.bam -part rparts.list 
	samtools view -F 4 rnasr.bam | sam2seqparts.pl -name rnasr -type sam_single -in stdin -part rparts.list 
	
  fagff2sam.pl -fa aphid_est.fa.gz -gff aphid_est.gff.gz > aphid_est.sam
	sam2seqparts.pl -type sam_est -in aphid_est.sam -part rparts.list 
			# for .longfa output
			
	/not/ sam2seqparts.pl -nam=$nam -type sam_pair -in paired.sam -part rparts.list -sort
	sam2seqparts.pl -type bam_pair -in rnaseq_pe.bam -part rparts.list 
	sam2seqparts.pl -type splitpairs -in=rnaseq_pe.bam  -part rparts.list 

=item  location filter for huge rRNA read counts

  ** main problem w/ velvet assembly is some scafs are full of rRNA reads, millions
  -- these choke denovo assembly. Need option to remove/segregate these to rrna.sam
     from input rrna.gff locations

  -- also filter out/skip small scaffolds < gene size, esp w/ millions of reads .. cant assemble
  
  see 
    $workd/scripts/samfiltergff.pl  -mate -in $sam -over $workd/misc/repeat_rrna.gff \
     -filter $nam.sam.rrna -notfilter $nam.norna.sam

  gzgrep 'rRNA.repeat' $workd/misc/repeatmasker/aphid2asm.repmask.gff.gz | \
  perl -ne'($r,$b,$e)=(split)[0,3,4]; $w=1+$e-$b; $rs{$r}+=$w; $n++; \
  END{ foreach $r (sort keys %rs) { $w=$rs{$r}; print "$r/\t$w\n"; }}'\
  | ggrep -F -f velout/notdone4.big.list - | sort 

  cat velout/notdone4.big.list | perl -pe's,/,\t,' | ggrep -F -f - $workd/rnas/bams/allbam.scaf.counts | sort
  Failed or undone big-read scaffolds in velvet assembly, all w/ high rRNA counts
  
  .............   Size    Reads       rRNA_bp
  Scaffold1027    87757   31716848       8821
  Scaffold1045    27010   26262268       10143
  Scaffold1059    45617   16296905       5093
  Scaffold12      1798989 13913482       5220
  Scaffold1361    14675   18685979       6365
  Scaffold1393    12057   14427335       4222
  Scaffold1402    53150   10926465       3370
  Scaffold1534    10111   14468896       4819
  Scaffold1732    9585    11980908       3983
  Scaffold1776    9503    9083407        2634
  Scaffold2105    5120    4611128        1943
  Scaffold3395    4537    5733166        1965
  Scaffold60      1187582 16860448       5307
  Scaffold904     41307   30560967       11579
  Scaffold922     44957   29853581       10765
    
=item sam_partitioner  -in .sam

  this is for pe reads, handles pair .fa2 output including mappings to alt scaffolds;
  -- still has problems: memory overload, file handles overload, slow

=item postprocess pairs in sam_part .fa2

  maybe all via do bam_part, then postprocess this way to get paired .fa2, as
    sam2seqparts -type splitpairs -part rparts.list -name xxx  
  bam_part however will miss pairs split over scaffolds and pairs w/ one unmapped
  
      my $ppairs='($d)=m,>([^/\s]+),; s/\t/\n/; print $lv,$_ if ($d eq $ld); ($lv,$ld)=($_,$d);' ;
      my $scmd="sort $fa2 | perl -ne '$ppairs' > $fa2.sort;";
      my $ok= system($scmd);  #?? this died or failed part way thru list.. dont know why
 
      # reverse this for .fa1 unpaired
      my $punpair='($d)=m,>([^/\s]+),; s/\t/\n/; $lpair=$pair; $pair=($d eq $ld)?1:0;
          putu(); ($lv,$ld)=($_,$d);
          END{ putu(); } sub putu{ unless($pair or $lpair) { print $lv if $lv; } }' ;
      my $scmd="sort $fa2 | perl -ne '$punpair' > $fa2.fa1;";
      my $ok= system($scmd);  #?? this died or failed part way thru list.. dont know why

   ^^^ sam2seqparts.pl -nam=$nam -type splitpair -part rparts.list 
      
   .. fixparts: names, split vs old split for aphidpe_SRRxxx
      per rparts/

gfind rparts/ -name aphidpe_SRR071347.bam.fa2.pairs -execdir $aphid2/rnas/rpartfix2renam.sh \;

#!/bin/bash
# rpartfix2renam.sh
# for nam in ( aphidpe_SRR075803 aphidpe_SRR075802 )
mv aphidpe_SRR075803.fa2 aphidpe_SRR075803.fa2.tmp
mv aphidpe_SRR075803.fa2.pairs aphidpe_SRR075803.fa2
cat aphidpe_SRR075803.fa2.unpair >> aphidpe_SRR075803.fa1

mv aphidpe_SRR075802.fa2 aphidpe_SRR075802.fa2.tmp
mv aphidpe_SRR075802.fa2.pairs aphidpe_SRR075802.fa2
cat aphidpe_SRR075802.fa2.unpair >> aphidpe_SRR075802.fa1

# for nam in (aphidpe_SRR071347)
mv aphidpe_SRR071347.bam.fa2 aphidpe_SRR071347.bam.fa2.tmp
mv aphidpe_SRR071347.bam.fa2.pairs aphidpe_SRR071347.fa2
mv aphidpe_SRR071347.bam.fa2.unpair aphidpe_SRR071347.fa1
      

      
=item bam_partitioner  -in .bam

  for single reads: quicker, simpler: 'samtools view $bam $loc'
  also use this for parts of .sam format; both at once?
  
  foreach bam ( $bamset )
    set nam=`basename $bam .bam`
    sam2seqparts.pl -nam=$nam -type bam -in $bam -part rparts.list -sort -debug
  end

=item parts.list

  * not working well enough: get some very large rnaseq.fa (>1GB) combo_ sets that
    overwhelm velvet/oases assembly on e.g. 48GB memory systems.
  * need partitions based on rnaseq read counts?
  
cat ../genome/aphid2asm.fa.count | grep '^Scaff' | perl -ne'($r,$n)=split; if($n<=499999){ $sn+=$n; $nt++;} else
{ $sb+=$n; $nb++;} END{print "nsmall=$nt; bsmall=$sn; nbig=$nb; bbig=$sb\n";}' | head
nsmall=23630; bsmall=263246589; nbig=294; bbig=278428882

parts= ($accession, $base_dir, $partitioned, $partition_dir, $combolist)

cat ../genome/aphid2asm.fa.count | grep '^Scaff' | env base=rparts cut=999999 perl -ne\
'BEGIN{ $base=$ENV{base}||"parts"; $cut=$ENV{cut}||999999; $cutlo=int($cut/2); $cuthi=int($cut*2); $icom=0; } \
END{ putp($cnam,$clist) if($clist); } \
($r,$len)=split; if($len>$cutlo) { putp($r, ""); } \
else { $cw= $len + $com[$icom]; \
if($cw > $cuthi) { putp($cnam,$clist) if($clist);  $cnam=$clist=""; \
$icom++; $com[$icom]= $cw=$len; } else { $com[$icom]= $cw; } \
$cnam="combo_".$icom; $clist.="," if($clist); $clist.="$r"; } \
sub putp { my($r,$clist)=@_; \
if($r =~ /combo_/ and $clist and not $clist =~ /,/) { $r=$clist; $clist=""; } \
$rpath="$base/$r"; \
print join("\t",$r, $rpath, "N", "", $clist),"\n"; } '\
> rparts.list

set nam=aphidpe_SRR075802 ; samtools view bams/$nam.bam | \
../scripts/sam2seqparts.pl -nam=$nam -in stdin -part tparts.list -sort -debug > & log.tpart1 &

##.........................................
## replace partition by scaff size w/ partition by rnaseq count, for easier denovo assembly

foreach bam (bams/*.bam) samtools idxstats $bam > $bam.count
perl -e'print "cut -f 1,2"; $c=3; for $i (1..20) { print ",$c"; $c+=4; }print "\n"; '

paste bams/*.bam.count | \
cut -f 1,2,3,7,11,15,19,23,27,31,35,39,43,47,51,55,59,63,67,71,75,79 |\
perl -ne'($r,$w,@c)=split; $sc=0; map{$sc+=$_}@c; print join("\t",$r,$w,$sc),"\n" if($sc>0);' \
> allbam.scaf.counts

cat allbam.scaf.counts | cut -f 1,3 | env base=sparts cut=1999999 perl -ne\
...
> sparts.list

#...........................................

  #!/bin/bash
  # makeparts.sh
  
  workd=/export/udisk3/work/aphid/
  scripts=$workd/scripts/
  bamset="bams/*.bam *est.sam"
  partlist=sparts.list
  
  for bam in $bamset; {
    case "$bam" in
    *aphidpe_*)
    $scripts/sam2seqparts.pl -nodupl  -format "fasta,sam" -type bam_pair -in $bam -part $partlist -debug
    $scripts/sam2seqparts.pl -nodupl  -type splitpairs -in $bam -part $partlist -debug  
    ;;
    *est.sam)
    $scripts/sam2seqparts.pl -nodupl -format "fasta" -type sam_est -in $bam -part $partlist -debug 
    ;;
    *)
    $scripts/sam2seqparts.pl -nodupl -format "fasta,sam" -type bam_single -in $bam -part $partlist -debug 
    ;;
    esac
  }


=cut

use strict;
use warnings;
use Getopt::Long ;
use File::Basename;

my $MAX_OPENFILE = 199; # what? 256 descriptors?
my $domake= 1; # parts dirs
my $do_sort2= 1; # debug??
my $WARN_NOPATH=0;
my $NO_SEQDUP= 0;

my(%ids);
my($n_in, $n_read, $n_aln, $n_aln_mate1, $n_pair, $n_pairok, $n_pairfar, 
   $n_strand, $n_map0, $n_map1, $n_mapn, $n_loident,
   $n_intron, $len_intron, $n_matenear, $len_mate, $len_read) = (0) x 30;

my($ok, $intype, $inname, $format, $dosam, $ispaired, $sorted, $partitions_file, $accession_list, $debug);
my( $outsuffix, $pairtmpsuffix, $pairsuffix, $unpairsuffix, $estsuffix)
    = (".fa",".tmpfa2",".fa2",".fa1", ".longfa");
my @IN;
my %getmate;

# $format="fasta2,sam"; # both output?
$format="fasta"; # both output? change to fa2,fa,sam == outsuffix
# $intype="sam"; # or bam, only for .fa1

my $optok= &GetOptions (
  "in=s"=>\@IN,  
  'type=s', \$intype,  # use also for pair/single ?
  'suffix=s', \$outsuffix, # use format instead? : .fa, .fa2, .fa1, .longfa
  'name=s', \$inname, # for stdin, outname
  "format=s"=>\$format,  # expect only: "fasta" or "fasta,sam" ; use also for pair/single ?
  'partitions_file|parts=s' => \$partitions_file,
  "sorted!"=>\$sorted,   
    # these for combo_ parts? not used here
  'accession_list=s' => \$accession_list,
  "noduplicates!"=>\$NO_SEQDUP,  
  "n|debug!"=>\$debug,  
  );

push(@IN, @ARGV); # any remainder;

die "USAGE: $0 [-name bob -sorted] -part partitions.list -type sam|bam|splitpairs -in stdin|my1.sam my2.sam\n"
unless($optok and -f $partitions_file and $intype);

## should these be set for each @IN infile below, in @files_to_process ?
$ispaired= ($intype =~ m/fa2|fasta2|pair|pe/ or $format =~ m/fa2|fasta2|pair|pe/) ? 1 : 0;
$ispaired=0 if($intype =~ m/fa1|fasta1|single|sr/ or $format =~ m/fa1|fasta1|single|sr/);

$outsuffix= $estsuffix if($outsuffix eq ".fa" and not $ispaired # FIXME: long paired est
  and ($IN[0] =~ /est/i or $intype =~ m/est|long/i or $format =~ m/est|long/i)); #? or use seq length

$dosam=1 if($format =~ /sam/); # always do fasta ??


my %base_directories_to_partitions;
my %base_directories_info;
open (my $fh, $partitions_file) or die "Error, cannot open $partitions_file";
while (<$fh>) {
    chomp; next unless(/^\w/); # if(/^#/);
    my ($accession, $base_dir, $partitioned, $partition_dir, $combolist) = split (/\t/);
    $base_directories_info{$base_dir}{haspart}= ($partitioned eq 'Y');
    $base_directories_info{$base_dir}{accession}= $accession;
    $base_directories_info{$base_dir}{combolist}= $combolist;
    if ($partitioned eq 'Y') {
        $base_directories_to_partitions{$base_dir}->{$partition_dir} = 1;
    } else {
       #? $base_directories_noparts{$base_dir} = 1;
    }
}
close $fh;


my @files_to_process = ();
if($intype =~ /sam|bam|split/) { # gff ?
  foreach my $infile (@IN) {
    #? stdin ok# unless(-f $infile) { die "missing IN=$infile\n"; next; }
    my $oname= ($infile =~ /stdin|-/ and $inname) ? $inname : basename($infile);
    my $ftype= $intype; # set from filename ?
    #? should set here ispaired, outsuffix? in @files_to_process
      
    push(@files_to_process, 
      {  type => $ftype,
         basename => $oname,
         file => $infile,
     } ) ;
  }
}

my $partinfo= getPartInfo(); # base_directories_info


foreach my $file_struct (@files_to_process) {
  my $in_filename = $file_struct->{file};
  my $basename = $file_struct->{basename};
  my $type = $file_struct->{type};
  # my $in_fasta_list = $type eq 'fagff' ? $file_struct->{fasta_list} : undef;
  
  warn "# process $type input $in_filename\n" if $debug;
  
  if($type =~ /sam2fa/i) {  #?
    sam2fasta( $partinfo, $basename);
    
  } elsif($type =~ /sam/i) {
    sam_partitioner( $in_filename, $partinfo, $basename);
    
  } elsif($type =~ /bam/i) {
    bam_partitioner( $in_filename, $partinfo, $basename);

  } elsif($type =~ /split/i) {
    split_fasta( $partinfo, $basename);
  }
  
#   } elsif($type =~ /^gff/i) {
#     gff_partitioner( $in_filename, $partinfo, $basename);

}

#................

sub split_fasta {
  my($partinfo,$basename)= @_;
  
  die "Need -name basename\n" unless($basename);
  # only use for $ispaired
  die "Need -type pair\n" unless($ispaired);
  # use -name=basename.fa2 w/ suffix here?
  
  my $suffix=($ispaired)? $pairtmpsuffix : $outsuffix;
  my %didin;
  my ($nin,$npair,$nunpair, $nfile)=(0) x 10; 
  foreach my $acc (sort keys %{$partinfo}) {  
    my $partition_dir;
        
    if( $partinfo->{$acc}{combo}) {  #  or $acc =~ /^combo_/
      $partition_dir= $partinfo->{$acc}{partition_dir}{ 1 };
    } else {
      $partition_dir= $partinfo->{$acc}{partition_dir}{ 1 };
## ** fixme later
#     if($partinfo->{$contig}{cstarts})  
    }

    my $infile = "$partition_dir/$basename";
    unless( -f $infile ) { $infile = "$partition_dir/$basename$suffix"; }
    next if($didin{$infile}++);  
    unless( -f $infile ) { warn "Missing $infile\n"; next; }
    
    my $outnam    = $infile; $outnam=~ s/$suffix$//;  # pairtmpsuffix  
    my $pairfile  = "$outnam$pairsuffix";  # ".pairs";
    my $unpairfile= "$outnam$unpairsuffix"; #".unpair";
    $nfile++;
    
    open(OP,">$pairfile") or die "writing $pairfile\n";
    open(OU,">$unpairfile") or die "writing $unpairfile\n";
    open(IN, "sort $infile |"); 
    my ($lv, $ld, $lpair)= (0) x 10;
    my %didseq;
    while(<IN>){  
      # assumes format is all lines of: ">id<tab>seq<endline>"; warn of others?
      unless(/^>/ and /\t/) { warn "#split_fasta: input not '>id.tab.seq' \n" if $debug; next; }
      my($d)= m,>([^/\s]+),; 
      s/\t/\n/;  $nin++;
      my $pair=($d eq $ld)?1:0;
      if($pair) { 
        my $skip=0;
        if($NO_SEQDUP) {
          my($lsq)= $lv =~ m/\n(\S+)/;
          my($nsq)=  $_ =~ m/\n(\S+)/;
          $skip=1 if(1 < $didseq{$lsq.$nsq}++); # should this test $nsq.$lsq also?
        }
        unless($skip) { print OP $lv,$_; $npair++; }
        
      } elsif($lv and not $lpair) { 
        print OU $lv; $nunpair++; 
      }
      ($lv,$ld,$lpair)=($_,$d,$pair);
    }
    if($lv and not $lpair) { print OU $lv;  $nunpair++;}
    close(OP); close(OU); close(IN);
  }
  warn "# split $basename ; nfile=$nfile; nin=$nin, npair=$npair, nunpair=$nunpair\n" if $debug;    
}


# see gff_partitioner
# this one works well; revise to handle paired fa, including pull mates on alt scaffold and unmapped
# --- use samtools view -f 0x04 + 0x01 = unmapped + paired


sub bam_partitioner {
  my($bamfile,$partinfo,$basename)= @_;

  my $oksamtools= `which samtools`;
  unless($oksamtools) { die "ERR: cannot find samtools to read $bamfile\n"; }
  
  my $suffix=($ispaired)? $pairtmpsuffix : $outsuffix;
  ## use pairtmpsuffix for "id<tab>seq" line temp format
  
  my %dido;
  foreach my $acc (sort keys %{$partinfo}) {  
    my $partition_dir;
    my @acc;
        
    if( $partinfo->{$acc}{combo}) {  #  or $acc =~ /^combo_/

      $partition_dir= $partinfo->{$acc}{partition_dir}{ 1 };
            
      my $combolist= $base_directories_info{$partition_dir}{combolist};
      ### assume combolist exists ...  
      ### samtools view bam ref1 ref2 ... << this works
      ### it can handle combo_106 with zillions/3037 of small scaffs ok; nreads=107113
      $combolist =~ s/[,;\s]+/ /g;
      @acc=($combolist);
      
    } else {
      ## skip if acc is part of combolist, done above
      next if($partinfo->{$acc}{combo}); # is combo
      
      @acc=($acc);
      $partition_dir= $partinfo->{$acc}{partition_dir}{ 1 };

## ** fixme later
#     if($partinfo->{$contig}{cstarts}) { # {haspart}
#       my @cstarts= @{ $partinfo->{$contig}{cstarts} };
#       my @cends  = @{ $partinfo->{$contig}{cends} }; # bad
#       for(my $i=0; $i <= $#cstarts; $i++) { 
#         if( $start >= $cstarts[$i] && $end <= $cends[$i] ) {
#          $min_start= $cstarts[$i];
#          $partition_dir= $partinfo->{$contig}{partition_dir}{ $min_start }; 
#          last;
#          }
#         }
#     } else 
      
    }
  
    # add .sam and .fa1 ?
    #  push @suffi, ".sam" if ($format =~ /sam/);
    my $onam= "$partition_dir/$basename$suffix";
    my $osam= "$partition_dir/$basename.sam";
    next if($dido{$onam}++); # from combo lists
    
    # warn "# write $onam\n" if $debug;
    open(O,">$onam") or do { warn "fail write $onam\n"; next; };
    open(SAM,">$osam") if($dosam);
  
    my ($nin,$nout)=(0,0); 
    my $sep= ($ispaired) ? "\t" : "\n";
    
    foreach my $acc2 (@acc) { # now only one acc2
      my (%didid, %didseq);
      open(S,"samtools view $bamfile $acc2 |"); 
      while(<S>) { 
        my $inline=$_;
        my($d,$f,@v)=split"\t";  $nin++; 
        my $sq= $v[7]; 
        my $dm= $d;
        
        ## add optional filter by dupl sequence; not paired here; need both mate seq for that
        if($f & 0x01) { # paired
          # $sep= "\t";
          my($chr, $matec)= @v[0,4];
          if($matec ne "=" and  $matec ne $chr ) { # $matec ne "*" and : keep these
           ## $getmate{$d} = [ $chr, $matec, $onam] ; #fix:
           $getmate{$d} = $onam; ## only need onam; save mem here
           }
          my $m=($f & 0x80)?2:1; 
          $dm="$d/$m";
          if($m == 2) { $sq=reverse($sq); $sq=~tr/ACGTacgt/TGCAtgca/; }   
        } else { # single
          next if($NO_SEQDUP and 1 < $didseq{$sq}++);
        }
         
        next if($didid{$dm}++);
        print O ">$dm$sep$sq\n";  $nout++; #? write as >id\tseq and later sort, change \t ??
        print SAM $inline if ($dosam);
        } 
      } 
    close(O);  close(SAM) if($dosam);
    warn "# wrote $onam; nin=$nin, nout=$nout\n" if $debug;    
    }
    
  bam_add_mispair_parts($bamfile,$partinfo,$basename) if($ispaired);  # %getmate
}



sub bam_add_mispair_parts {
  my($bamfile,$partinfo,$basename)= @_;

  my $unpairflag=5; ## 0x01 + 0x04
  my %outhands;  
  my ($nin,$nout,$nopen)=(0,0,0); 
  my $sep= ($ispaired) ? "\t" : "\n";
  
  my %didid;
  warn( "# unmapped mates: samtools view -f $unpairflag $bamfile\n") if $debug; 
  open(S,"samtools view -f $unpairflag $bamfile |"); 
  while(<S>) { 
    my($d,$f,@v)=split"\t";  $nin++; 
    if($getmate{$d}) {
      my $inline= $_;
      # my( $mchr, $thischr, $monam)= @{$getmate{$d}}; 
      my $monam= $getmate{$d};
      (my $mosam= $monam) =~ s/\.\w+$/.sam/; # dosam
      
      my $outh= $outhands{$monam};
      my $samh= $outhands{$mosam};
      unless(ref $outh) {
        if($nopen >= $MAX_OPENFILE) {
          foreach my $file (sort keys %outhands) { 
            my $outh= delete $outhands{$file};
            close($outh) if $outh; $nopen--;
            }        
          }
        warn "# unpair $d reopen $monam\n" if $debug;
        open( $outh,">>$monam") or do { warn "ERR: cant write $monam\n"; return; };
        $outhands{$monam}= $outh;
        $nopen++;
        if($dosam) {
        open( $samh,">>$mosam");
        $outhands{$mosam}= $samh;
        $nopen++;
        }
      }
      
      my $sq= $v[7]; 
      my $dm= $d;
      if($f & 0x01) { # paired
       my $m=($f & 0x80)?2:1; 
       $dm="$d/$m";
       if($m == 2) { $sq=reverse($sq); $sq=~tr/ACGTacgt/TGCAtgca/; }   
       }
       
      next if($didid{$dm}++);
      print $outh ">$dm$sep$sq\n";  $nout++; #? write as >id\tseq and later sort, change \t ??
      print $samh $inline if ($dosam);  
      }
    } 
    
  foreach my $file (sort keys %outhands) { 
    my $outh= delete $outhands{$file};
    close($outh) if $outh; $nopen--;
    }
     
  warn "# wrote unmapped ; nin=$nin, nout=$nout\n" if $debug;    
}


# sam2fasta simpler than _partitioner, read .sam per part, write .fa, .tmpfa2 // .fa2, .fa1
sub sam2fasta {
  my($partinfo,$basename)= @_;
  
  die "Need -name basename\n" unless($basename);
  ## die "Need -type pair\n" unless($ispaired);
  
  my $insuffix= ".sam"; # $samsuffix; << may not be .sam here ?? : .sam.rrna ??
  my $suffix=($ispaired)? $pairtmpsuffix : $outsuffix;
  my $sep= ($ispaired) ? "\t" : "\n";

  my %didin;
  my ($nin,$nout)=(0) x 10; 
  foreach my $acc (sort keys %{$partinfo}) {  
    my $partition_dir;
        
    if( $partinfo->{$acc}{combo}) {  # $acc =~ /^combo_/
      $partition_dir= $partinfo->{$acc}{partition_dir}{ 1 };
    } else {
      $partition_dir= $partinfo->{$acc}{partition_dir}{ 1 };
#     if($partinfo->{$contig}{cstarts}) {}  ## ** fixme later
    }

    my $infile = "$partition_dir/$basename";
    unless( -f $infile ) { $infile = "$partition_dir/$basename$insuffix"; }
    next if($didin{$infile}++);  
    unless( -f $infile ) { warn "Missing $infile\n"; next; }

    my $onam= "$partition_dir/$basename$suffix";
    ## next if($dido{$onam}++); # from combo lists
    open(O,">$onam") or do { warn "fail write $onam\n"; next; };
    
    my (%didid, %didseq);
    open(S, $infile); 
    while(<S>) { 
      my $inline=$_;
      my($d,$f,@v)=split"\t";  $nin++; 
      my $sq= $v[7]; 
      my $dm= $d;
      
      if($f & 0x01) { # paired
        # my($chr, $matec)= @v[0,4];
        my $m=($f & 0x80)?2:1; 
        $dm="$d/$m";
        if($m == 2) { $sq=reverse($sq); $sq=~tr/ACGTacgt/TGCAtgca/; }   
      } else { # single
        next if($NO_SEQDUP and 1 < $didseq{$sq}++);
      }
       
      next if($didid{$dm}++); #? drop
      print O ">$dm$sep$sq\n";  $nout++; #? write as >id\tseq and later sort, change \t ??
      } 
        
    close(O);  close(S);
    warn "# wrote $onam; nin=$nin, nout=$nout\n" if $debug;    
    }

}




sub sam_partitioner {
  my($infile,$partinfo,$basename)= @_;
  
  my %outhands= ();
  my %outparts= (); # return file list?
  my($lcontig, $last_dir);
  my (%ocache,%didid);
  
  my $insam;
  if($infile =~ /stdin|-/i) { $insam= *STDIN; }
  else { open(SAM,$infile) or die "bad input: $infile"; $insam=*SAM; }
  
  while (<$insam>) {
    next unless (/^\w/); # in case of ^@SQ .. sam header

    my $inline=$_;    
    my ($qid, $flag, $chr, $cloc, $mapq, $cigar, $matechr, $mateloc, $isize, $seq, $qual, @opt)
      = split"\t";
    $n_in++;
    
    my $f_pair= ($flag & 0x0001);
    my $f_pairok= ($flag & 0x0002); #??
    my $f_mismatch= ($flag & 0x0004); #  f_mismatch == $chr eq "*"
    my $f_mismate= ($flag & 0x0008);
    # my $f_isrev= ($flag & 0x0010);     # 16
    # my $f_revmate= ($flag & 0x0020); # 32
    my $f_first = ($flag & 0x0040); # 64; mate ids
    my $f_second= ($flag & 0x0080); # 128

    ## NOW locate right file handle from chr, cloc  and/or matechr/mateloc
    ## for paired reads to .fa2, want to use both locations for output 

    my($nout,$suffix,$outval)=(0, "","");
    
    # if($f_mismatch) {  next; } # not here, print seq if pairok
    my $mate=($f_second)?2:1;
    my $faid= "$qid/$mate";
    
    # * FIXME: save all mapped mate ids, capture unmapped pair mates; put both in .fa2
    # ? is that being done now, to .fa1 ??
    
      #  $f_pairok and not $f_mismate : no, only need 1 mapped ok
    # if($f_pair and (!$f_mismatch or !$f_mismate)) 
    if($f_pair and $f_pairok) 
    { 
      # options 1. print as encountered, with ID /1 /2 for mate, then sort by ID for velvet
      #  2. expect input sam sorted by ID ?? no
      #  3. hash-cache all pair id, seq, then at end print to files... not good
      #   if(pairok) print fa2 >ID/1 <tab> seq ; >ID/2 <tab> seq; ... sort fa2 ... perl -pie's/\t/\n/' .fa2
      
      $matechr= $chr if($matechr eq "=");
      
      # revcomp mate for velvet format: only?
      if( $f_second ) { $seq= reverse($seq); $seq =~ tr/ACGTacgt/TGCAtgca/; }
      
      $nout=2;
      $suffix= $pairsuffix; # ".fa2";  # need option for suffix ; or rename later
      $outval= ">$faid\t$seq\n";
      
    } elsif($f_pair) {
      $nout=1;
      $suffix= $unpairsuffix; # ".fa1"; # need option for suffix ; or rename later
      $outval= ">$faid\n$seq\n";
      
    } else {
      $nout=1;
      $suffix= $outsuffix; # ".fa1"; # need option for suffix ; or rename later
      $outval= ">$faid\n$seq\n";
    }
   
    my $didout=0;
    for(my $iout=0; $iout<$nout; $iout++) 
    {
      my ($contig, $start, $end);  
      if($iout > 0) {
        ($contig, $start, $end)= ($matechr, $mateloc, $mateloc+length($seq)); 
      } else {
        ($contig, $start, $end)= ($chr, $cloc, $cloc+length($seq));  
      }
      next if($contig eq '*'); # this one not mapped but mate may be
      ## next unless($partinfo->{$contig}); # see below
      
      my $partition_dir="";
      my $min_start= 0;
      
      if($partinfo->{$contig}{cstarts}) { # {haspart}
        my @cstarts= @{ $partinfo->{$contig}{cstarts} };
        my @cends  = @{ $partinfo->{$contig}{cends} }; # bad
        for(my $i=0; $i <= $#cstarts; $i++) { 
          if( $start >= $cstarts[$i] && $end <= $cends[$i] ) {
           $min_start= $cstarts[$i];
           $partition_dir= $partinfo->{$contig}{partition_dir}{ $min_start }; 
           last;
           }
          }
      } else {
        $min_start = 1;
        $partition_dir= $partinfo->{$contig}{partition_dir}{ $min_start };       
      }
      
      unless($partition_dir) { warn "no path for $contig:$start-$end\n" if $WARN_NOPATH; next; }# these are scafs not in rparts.list, because of assembly mismatch       
      next if($didout > 0 and $partition_dir eq $last_dir); # next/last : skip, same part

      # reduce open file handles here
      my (%moreout, @moresuf);
      if($iout > 0 and $partition_dir ne $last_dir)  { # dont open new file here??
        push( @{$ocache{$partition_dir}{$suffix}}, $outval); 
        next; 
        
      } elsif( $iout == 0 and $ocache{$partition_dir}) {
        @moresuf= sort keys %{$ocache{$partition_dir}};
        foreach my $suf (@moresuf) {
          my $olist= delete $ocache{$partition_dir}{$suf};
          push  @{$moreout{$suf}}, @$olist; 
        }
        delete $ocache{$partition_dir};
      }
      
      unless( -d $partition_dir ) { 
        warn "bad path for $contig:$start-$end: $partition_dir\n";
        next; ## return; # next/return/die ?
      } 
  
      # close outh after pass contig/ref ? assume ingff sorted by ref?
      if($sorted and $iout == 0 and $lcontig and $contig ne $lcontig) {
        my $comacc= $partinfo->{$lcontig}{combo} || $lcontig;
        foreach my $file (sort keys %{$outhands{$comacc}}) { 
          next if ($file =~ m,combo,);  # this is problem; can we close some of these? reduce number?
          warn "# close $file\n" if $debug;
          my $outh= delete $outhands{$comacc}{$file};
          close($outh) if $outh; 
        }
      }
      
      if( $iout == 0 ) { $lcontig= $contig; }
      $last_dir= $partition_dir;
      
      # fixme: change basename.suffix for out type: pairs = .fa2, single = .fa1, other?
      
      my @suffi=($suffix); # if ($format =~ /fasta/);
      if(@moresuf) { push(@suffi, grep { $_ ne $suffix } @moresuf); }
      
      push @suffi, ".sam" if ($dosam and $iout == 0); # $format =~ /sam/ 
      foreach my $suf (@suffi) {
      my $part_file = $partition_dir . "/$basename" . $suf;
      next if ($suf ne ".sam" and $didid{$part_file}{$faid}); 

      my $comacc= $partinfo->{$contig}{combo} || $contig;
  
      ## FIXME: too many files open errors; reopen as >>part_file
      ## FIXME: need contig -> combo_n map here for files
      my $outh= $outhands{$comacc}{$part_file};
      unless($outh) {
        warn "# at $comacc, open $part_file\n" if $debug;
        open( $outh,">>$part_file") or do { die "ERR: at $comacc, cant write $part_file\n"; return; };
        $outhands{$comacc}{$part_file}= $outh;
        $outparts{$part_file}++;
      }
  
      if($suf =~ /sam/) { print $outh $inline if($iout==0); }
      else { 
        if($suf eq $suffix and not $didid{$part_file}{$faid}++) { print $outh $outval;  $didout++; } 
        if($moreout{$suf}) { foreach my $moreout (@{$moreout{$suf}}) { print $outh $moreout; } }
        }
      }
      
    }  # iout
      
  }
  
  close($insam);
  foreach my $contig (sort keys %outhands) { 
    foreach my $file (sort keys %{$outhands{$contig}}) { 
      warn "# close $file\n" if $debug;
      my $outh= delete $outhands{$contig}{$file};
      close($outh) if $outh; 
      }
    }
    
  if($do_sort2) {
    #   if(pairok) print fa2 >ID/1 <tab> seq ; >ID/2 <tab> seq; ... sort fa2 ... perl -pie's/\t/\n/' .fa2
    foreach my $fa2 (sort keys %outparts) {
      next unless($fa2 =~ /$pairsuffix$/ and -f $fa2);
      ## fixme: need to filter out mistake unpaired mistakes : ~2,500 of 51,000 pairs in one case
      ## my $scmd="sort $fa2 | perl -pe 's/\t/\n/;' > $fa2.sort; mv $fa2.sort $fa2";
      my $perle='($d)=m,>([^/\s]+),; s/\t/\n/; print $lv,$_ if ($d eq $ld); ($lv,$ld)=($_,$d);' ;
      my $scmd="sort $fa2 | perl -ne '$perle' > $fa2.sort;";
      ## later: mv $fa2.sort $fa2
      my $ok= system($scmd);  #?? this died or failed part way thru list.. dont know why
    }  
  }
  
}




sub getPartInfo
{
  my $partinfo= {};
  
  foreach my $base_dir (sort keys %base_directories_info) {
  
    my $accession = $base_directories_info{$base_dir}{accession};      
    unless (-d $base_dir) { my $ok=0;
      $ok= mkdir($base_dir) if($domake);
      die "Error, missing $accession folder: $base_dir" unless $ok;
    }
  
    my $haspart= $base_directories_info{$base_dir}{haspart};
  
    if( $haspart ) {
      my $partition_dirs_href = $base_directories_to_partitions{$base_dir};
      my $max=0;
      my @cstarts=();
      my @cends=();
      foreach my $partition_dir (sort keys %$partition_dirs_href) {
        $partition_dir =~ /(\w+)_(\d+)-(\d+)$/ or die "Error, cannot extract coords from $partition_dir";
        my ($segref, $min_start, $max_end) = ($1,$2,$3);
        # urk can have overlap where next min_start < last max_end
        push(@cstarts, $min_start);
        push(@cends, $max_end);
        $max= $max_end if($max_end > $max);
        $partinfo->{$accession}{partition_dir}{ $min_start } =  $partition_dir;
        } 
        ## need to have @{parts} ordered by min_start ... or sort below
      # @cstarts = sort { $a <=> $b } (@cstarts,$max);
      my @ord = sort { $cstarts[$a] <=> $cstarts[$b] } (0..$#cstarts);
      @cstarts= @cstarts[@ord];
      @cends  = @cends[@ord];
      
      $partinfo->{$accession}{cstarts} = \@cstarts;
      $partinfo->{$accession}{cends} = \@cends;
    }
  
          
    else {  # no parts
      ## my $part_file = "$base_dir/$basename";
      my %accs=();
      
      if($accession =~ /^combo_/) {
        my $lh;
        my $combolist= $base_directories_info{$base_dir}{combolist};
        
        if( $combolist ) {
          %accs= map{ $_,1 } split( /[,;\s]+/, $combolist );
          
        } elsif( $accession_list && open( $lh, "$base_dir/$accession_list") ) {  
          while(<$lh>){  s/^>//; if(m/^(\S+)/){ $accs{$1}++;} } close($lh); 

#         } elsif( $genome && open( $lh, "$base_dir/$genome") ) { 
#           while(<$lh>){ if(m/^>(\S+)/){ $accs{$1}++;} } close($lh); 

        } else {
          die "missing --accession_list or --genome for combo";
        }
        
        foreach my $acc (sort keys %accs) {
          $partinfo->{$acc}{partition_dir}{ 1 }= $base_dir;
          $partinfo->{$acc}{combo}= $accession;
        }
        
      } else {
        $partinfo->{$accession}{partition_dir}{ 1 }= $base_dir;
      }
      
    }
    
  }
  
  return $partinfo;
}


## revise pe id read/id counts to do both of pair: append /1 /2 to ids
#      p=0x1 (paired), P=0x2 (properly paired), u=0x4 (unmapped),
#      U=0x8 (mate unmapped), r=0x10 (reverse), R=0x20 (mate reverse)
#      1=0x40 (first), 2=0x80 (second), 's' = 0x100 (not primary), 
#      f=0x200 (failure) and d=0x400 (duplicate). 
#     
#   $n_pair++ if($f_pair);
# 
#   if($f_mismatch) { $ids{$qid}= 0; } # below: $n_map0++; 
#   else { 
#     $n_aln++;
#     unless($f_pair and $f_second and $ids{$qid}) { # $f_isrev 
#       $n_aln_mate1++;
#       $ids{$qid}++; # is this right for pairs w/ same id
#     }
#     
#     my $len  = length($seq);
#     $len_read += $len;
#     my $cend = $cloc; # calc correct end from align/intron spans
#     my(@aln,@intr,@itype,$instrand,$inmismat);
#     
#     if($cigar =~ m/^(\d+)M/ ) { @aln=($1); $cend += $1; }
#     elsif($cigar =~ /^\d+[HSN]/) { @aln=(0); } # need for nSnM cigar
#     while($cigar =~ m/(\d+)([DHNS])(\d+)M/g) { 
#       my($bi,$bt,$bx)=($1,$2,$3);
#       push(@intr,$bi);  push(@itype,$bt);
#       push(@aln,$bx); $cend += $bi + $bx; 
#       if($bt eq "N"){ $n_intron++; $len_intron+=$bi; }
#       }
#   
#     my $nx= scalar(@aln);
#     my $alen=0; map{ $alen+=$_ }@aln; 
#     foreach (@opt) { 
#       if(m/NM:i:(\d+)/){ $alen -= $1; } 
#       elsif(m/XS:A:([\+\-])/){ $instrand=$1; $n_strand++;  } 
#       # elsif(m/NS:i:(\d+)/) { $inmismat=$1; } #new, if XS, = Mismatches within min_anchor_len of a splice junction, > 0 is poor
#       # elsif(m/RG:A:(\S+)/){ $readgroup=$1; } #new
#       }
# 
#     # my $score= $alen/$len;
#     $n_loident++ if( $len<1 || $alen/$len < $MIN_IDENT);
# 
#     if($f_pair and not $f_mismate) {
#       $n_pairok++ if($f_pairok);
#       if($matechr ne "=" and $matechr ne "*") { $n_pairfar++; }
#       elsif($matechr eq "=") {
#         $n_matenear++;
#         $len_mate += ($mateloc < $cloc) ? $cloc - ($mateloc+$len) : $mateloc - $cend;
#       }
#     }
#       
#   
#   }

 


#...............

## convert to sam_partitioner

# sub gff_partitioner {
#   my($inh,$partinfo,$basename)= @_;
#   
#   my %outhands= ();
#   my($lcontig);
#   
#   # open(my $inh,"$ingff") or do { warn "cant read $ingff\n"; return; };  
#   
#   while (<$inh>) {
#   
#       unless (/^\w/) { 
#         next;
#       }
#     
#       my @x = split (/\t/);
#       my ($contig, $lend, $rend) = ($x[0], $x[3], $x[4]);
#             
#       my $partition_dir="";
#       my $adjust_to_1= 0;
#       my $adjust_coord = 0;
#       my $min_lend= 0;
#       
#       if($partinfo->{$contig}{cstarts}) { # {haspart}
#         my @cstarts= @{ $partinfo->{$contig}{cstarts} };
#         my @cends  = @{ $partinfo->{$contig}{cends} }; # bad
#         for(my $i=0; $i <= $#cstarts; $i++) { 
#           if( $lend >= $cstarts[$i] && $rend <= $cends[$i] ) {
#            $min_lend= $cstarts[$i];
#            $partition_dir= $partinfo->{$contig}{partition_dir}{ $min_lend }; 
#            $adjust_to_1= 1;
#            $adjust_coord = $min_lend - 1;
#            last;
#            }
#           }
#       } else {
#         $min_lend = 1;
#         $partition_dir= $partinfo->{$contig}{partition_dir}{ $min_lend };       
#       }
#       
#       unless( -d $partition_dir ) { 
#         warn "bad path for $x[2]:$x[1] $contig:$lend-$rend: $partition_dir\n";
#         next; ## return; # next/return/die ?
#       } 
#   
#       # close outh after pass contig/ref ? assume ingff sorted by ref?
#       if($lcontig and $contig ne $lcontig) {
#         my $comacc= $partinfo->{$lcontig}{combo} || $lcontig;
#         foreach my $file (sort keys %{$outhands{$comacc}}) { 
#           next if ($file =~ m,combo,);
#           ## this is bad for combo's
#           warn "# close $file\n" if $debug;
#           my $outh= delete $outhands{$comacc}{$file};
#           close($outh) if $outh; 
#         }
#       }
#       
#       $lcontig= $contig;
#       
#       my $part_file = $partition_dir . "/$basename";
#       my $comacc= $partinfo->{$contig}{combo} || $contig;
# 
#       ## FIXME: need contig -> combo_n map here for files
#       my $outh= $outhands{$comacc}{$part_file};
#       unless($outh) {
#         warn "# at $comacc, open $part_file\n" if $debug;
#         open( $outh,">$part_file") or do { die "ERR: at $comacc, cant write $part_file\n"; return; };
#         $outhands{$comacc}{$part_file}= $outh;
#       }
# 
#       if ($adjust_to_1) {
#         $x[3] -= $adjust_coord;
#         $x[4] -= $adjust_coord;
#       }
#       
#       print $outh join ("\t", @x); # gff out
#       
#   }
#   
#   close($inh);
#   foreach my $contig (sort keys %outhands) { 
#     foreach my $file (sort keys %{$outhands{$contig}}) { 
#       warn "# close $file\n" if $debug;
#       my $outh= delete $outhands{$contig}{$file};
#       close($outh) if $outh; 
#       }
#     }
# }
# 



# print "# SAMstats ; input n=$n_in\n";
# 
# # counts from unique %ids
# ($n_read, $n_map0, $n_map1, $n_mapn)= (0) x 10;  
# 
# ## dont use values() makes big array; use ?? hash iterate
# # foreach my $vid ( values( %ids ) ) 
# while ( my($id, $vid) = each %ids) {
#   $n_read++;
#   if($vid == 1) { $n_map1++; }
#   elsif($vid == 0) { $n_map0++; }
#   elsif($vid > 1) { $n_mapn++; }
# }
# 
# $n_read||=1; 
# my $n_mapread= ($n_read - $n_map0) || 1;
# $len_read= ($n_aln<1) ? 0: sprintf "%.1f", $len_read / $n_aln;
# $len_intron= ($n_intron<1) ? 0: sprintf "%.0f", $len_intron / $n_intron;
# $len_mate= ($n_matenear<1) ? 0: sprintf "%.0f", $len_mate / $n_matenear;
# 
# foreach my $k ( qw(len_read n_read n_aln n_strand n_intron len_intron
#       n_map0 n_map1 n_mapn n_pair n_pairok n_pairfar len_mate n_loident ) ) 
# {
#   my $v= eval("\$".$k);
#   my $p=0;
#   if($k =~ /n_map/) { $p= $v / $n_read; }
#   elsif($k =~ /n_aln/ and $n_pair>0) { $p= $n_aln_mate1 / $n_mapread; }
#   elsif($k =~ /n_aln/) {  $p= $v / $n_mapread; }
#   elsif($k =~ /n_pairok|n_pairfar/ and $n_pair>0) { $p= $v / $n_pair; }
#   elsif($k =~ /^len_/) { }
#   elsif($n_aln>0) { $p = $v / $n_aln; }
#   
#   printf "%10s\t%8d",$k, $v;
#   printf "\t%5.1f%%", 100*$p unless($k =~ /n_read|^len_/);
#   if($k=~/^len_/) { }
#   elsif($k=~/n_aln/) { print " of mapped_read"; }
#   elsif($k=~/n_map/) { print " of n_read"; }
#   elsif($k=~/n_pairok|n_pairfar/) { print " of n_pair"; }
#   elsif($k=~/loident|_pair|n_strand|^n_intron/) { print " of n_aln"; }
#   print "\n";
# }  
# print "\n";


#...........................................
