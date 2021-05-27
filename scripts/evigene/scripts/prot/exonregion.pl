#!/usr/bin/env perl
# exonregion.pl

=item about

  preprocess tblastn to exoner.fagff : gff/query.aa/genome.na for each region of interest

  $evigene/scripts/exonregion.pl -v -gff nasvit1-bestprot1.gff -aa nvit11prots.aa -genome $nas11/genome/nasvit1asm.fa -out nasvit1-bestarpod2.tab


=item FIXME: overlapped spans

  revise to compress overlapped regions; 2+ aa per genome span: comma-sep in aaid, aa columns?
  require loc-sorted mRNA inputs

  grep mRNA nasvit1-bestprot1.gff | sort -k1,1 -k4,4n -k5,5nr |\
  $evigene/scripts/exonregion.pl -v -sorted \
    -gff stdin -aa nvit11prots.aa -genome $nas11/genome/nasvit1asm.fa -out nasvit1-bestarpod2.tab

=cut

use strict;
use lib("/bio/argos/common/env perl/lib/");
use Bio::DB::Fasta;
use Getopt::Long;

my $MIN_CDS_LEN = 60;
my $MIN_AA_LEN = 20;
my $MIN_ESTSPAN = 60; ## == MIN_CDS_LEN ?
my $GENE_EXPAND = 1000;  
my $EXONGAP = 'nnn'; # or other? for EST pasting
use constant BIGXSPAN => 9999; # for extra-genic genes; cant handle long inter-genic introns w/o evidence      

my $verbose= 0;
my $outh = *STDOUT;
my $gffin= *STDIN;
my $genome_db= undef;
my $protein_db= undef;

my %args=();
my $ok= &GetOptions( \%args,
  "gff=s" ,
  "aa|protein=s" ,
  "genome|dnafastapath=s" ,
  "output=s" ,
  "mrnatype=s" ,"exontype=s" , 
  "xtype=s" , # extract/exonerate/... = Protein or EST/cDNA
  "o|offset=s" , 
  "stranded!" ,
  "sorted!" ,
  "MIN_CDS_LEN=i" , # not used
  "dropnnn!" ,  "v3!", "prettyfasta!", # not used
  "v|verbose!" , # == debug
  "debug!" ,
  );

my $xtype= lc($args{xtype}) || "protein";
my $gff= $args{gff} || shift(@ARGV);
my $genomefasta= $args{genome} || shift(@ARGV);
my $proteinfasta= ($xtype =~ /EST|cDNA/i) ? "" : ( $args{aa} || shift(@ARGV) );

die "Usage:  exonregion.pl [-verbose] -gff=annot.gff[.gz] -genome=genome.fasta/\n"
  ." opts: -protein=proteins.aa OR -xtype=EST|cDNA and gff=mRNA+exons, -exon=exon, -mrna=mRNA"
  unless($ok && $gff && -e $genomefasta and (-e $proteinfasta or $xtype =~ /EST|cDNA/i));

if($gff =~ /\.gz$/) { $ok= open(GFF, "gunzip -c $gff|"); $gffin = *GFF; }
elsif($gff =~ /stdin/) { $ok=1; $gffin = *STDIN; }
else { $ok= open(GFF,$gff); $gffin = *GFF; }

if($args{output}) {
  my $outf=$args{output};
  $ok= open(OUT,">$outf") or die "error opening $outf"; $outh = *OUT; 
}

my $debug= $args{debug};
my $sorted= $args{sorted};
my $stranded= $args{stranded}; # need this for cDNA seq revcomp ***

my $mrnaType= $args{mrnatype} || "mRNA,match";  $mrnaType= join"|", split(/[,;.\|]+/, $mrnaType);
my $exonType= $args{exontype} || "exon,HSP";  $exonType= join"|", split(/[,;.\|]+/, $exonType);

my @offset= ($args{o}) ? split(/[,;.]+/,$args{o}) : ();
   if(@offset == 1) { my $o=$offset[0]; @offset=( -$o, $o); }

my $prettyfa=($args{prettyfasta})?1:0;  ## not used here
my $DROPNNN= $args{dropnnn}?1:0;        ## not used here
$MIN_CDS_LEN= $args{MIN_CDS_LEN} || $MIN_CDS_LEN; ## not used here
$MIN_AA_LEN= $args{MIN_AA_LEN} || $MIN_AA_LEN;    ## not used here

if($args{v}){ $verbose= $args{v}; }
$verbose= 1 if($debug);

warn "# exonregion: gff=$gff, aa=$proteinfasta, seq=$genomefasta\n" if $verbose;

## FIXME sort out by > aa size : big ones take long time to exonerate ; do at end
## add end col of aa-size; dna-size ??

my $doaa = (-e $proteinfasta and $xtype !~ /EST|cDNA/i)? 1 : 0;
my $doest= ($xtype =~ /EST|cDNA/i)? 1: 0;

processGff($outh,$gffin,$genomefasta,$proteinfasta);

close($gffin);
close($outh);


sub processGff {
  my($outh,$gffHandle,$genomefasta,$proteinfasta)= @_;

  my($offa,$offb)= (@offset>1) ? @offset : (-$GENE_EXPAND, $GENE_EXPAND);
  my (@spanid, @spanaa, @xspanid, @xspanaa);
  my($lastref, $laststart, $laststop, $lastor, $lastxref, $lastxstart, $lastxstop, $geneid, $xgeneid, $nexon)= (0) x 20;
  
## FIXME: input sort by loc bad for EST with mrna inside mrna (long spans); exons misplaced after 2nd inner genes
## ok for prot: only mrna input, aa attached by id.
# cat $xgenes  ti.$tig.gff | sort -k1,1 -k4,4n -k5,5nr | exonregion ...
## instead collect all mrna 1st, then attach each incoming exon to one or more

  while(my $gffin= <$gffHandle>) {
  
    unless($gffin =~ m/^\w/) { next; }  
      
    my @gff= split "\t",$gffin,9;
    my ($ref,$type,$rb,$re,$ro) = @gff[0,2,3,4,6];
    #chomp($gff[8]);
    
    if($type =~ m/($mrnaType)/) {

      my $id= $gff[8];  
      if($id =~ m/Target=([^;\s]+)/) { $id=$1; } 
      elsif($id =~ m/ID=([^;\s]+)/) { $id=$1; }
      elsif($id =~ m/Parent=([^;\s]+)/) { $id=$1; }
      elsif($id =~ m/gene_id \W?([^;"\s]+)/) { $id=$1; } ## gtf/gff2
      else { $id =~ s/[;,\s].*$//; }
      $geneid= $id;
      
## problem here w/ some LONG aa that partly match to SHORT na region
## expand to fit all of aa?  should see Target begin/end to know where to expand

      my $geneaa = ($doaa) ? get_protein( $proteinfasta, $id) : "";
      my $aalen= length($geneaa);
      my $clen= 3 * $aalen;
      my $glen= 1 + $re - $rb;

      my($atb, $ate)= ($rb, $re);
      my($off1,$off2)=($offa,$offb);
      if($clen > $glen) { 
        my $d= $clen - $glen;
        $atb -= $d; $ate += $d;
        $d += 500; #?
        ($off1,$off2)= (-$d, $d);
      }
      
#      ## dont expand till putspan()
#      # ?? use only input rb,re or expanded span ? dont want to make huge span of adjacent genes
      
      if( $sorted and $lastref eq $ref and $atb < $laststop and $ate > $laststart) {
        $atb= $laststart; 
        $ate= $laststop if($laststop > $ate); 
        
      } else {
        putspan( $lastref, $laststart + $offa, $laststop + $offb, \@spanid, \@spanaa, $lastor) 
          if($lastref and @spanaa);
        putspan( $lastxref, $lastxstart + $offa, $lastxstop + $offb, \@xspanid, \@xspanaa, "") 
          if($lastxref and @xspanaa);
        @xspanid= @xspanaa=();
        @spanid= @spanaa=(); $nexon=0;
      }
      ($lastref, $laststart, $laststop, $lastor)= ($ref, $atb, $ate, $ro);
      $lastxstop= $lastxstart = 0;  $lastxref= $ref;
      
      my $regionid = "$ref:$atb-$ate";
      if($doest) { }
      elsif($geneaa =~ /\w/) { push( @spanid, $id); push( @spanaa, $geneaa); } # else error !
      else { print $outh join("\t","#errorspan",$id,($geneaa?"okaa":"missaa"),$regionid,"na"),"\n"; }
      
      
    } elsif( $doest and $type =~ m/($exonType)/) {
      # exon-gff per geneID extract to genecdna-fasta
      # for now assume exon lacks Parent=, i.e. anonymous exon span inside gene span      
      # if exon in gene span, put in @spanid, @spanaa ...
      # else put as intergeneic est, new id
      
      ## exon contained-in gene test:
      my $ingene= ( $ref eq $lastref and $rb >= $laststart + $offa and $re <= $laststop + $offb) ? 1:0;
      ## fix this, using overlap location bins for all mrna
      
      my $exondna = get_dna( $genomefasta, $ref, $rb, $re, $ro); # use ro  for revcomp ??
      
      if($ingene) { 
        # **?? add "NNN" gap between exon HSP append here? or not .. needs test
        if(@spanaa) { $spanaa[-1] .= $EXONGAP . $exondna; } 
        else { @spanaa=($exondna); @spanid=($geneid); }
        $nexon++;
        
      } else {
        ## save or not ?? should also force new $xgeneid for large extra-genic spans, eg: 
        #?? need to change laststart, laststop for xspan, doesn't belong in it
        $lastxstart= $rb if($lastxstart == 0);
        
        if( @xspanaa and ($ref ne $lastxref or $rb > $lastxstart + BIGXSPAN) ) { 
          putspan( $lastxref, $lastxstart + $offa, $lastxstop + $offb, \@xspanid, \@xspanaa, "") 
           if($lastref and @xspanaa);
          @xspanaa= @xspanid= ();
          $lastxstart= $rb;         
         } 
        if( @xspanaa ) { $xspanaa[-1] .= $EXONGAP . $exondna; }        
        else { @xspanaa=($exondna); @xspanid=( "xgene".++$xgeneid); }
        $lastxstop= $re;  $lastxref= $ref;
      }
    }

  }  # gffin
  
  putspan($lastref, $laststart + $offa, $laststop + $offb, \@spanid, \@spanaa, $lastor) if(@spanaa);  
  ## dont put trivial,small xspan
  putspan($lastxref, $lastxstart + $offa, $lastxstop + $offb, \@xspanid, \@xspanaa, "") if(@xspanaa);  
}


sub putspan {
  my($ref, $spanstart, $spanstop, $spanid, $spanaa, $spanor)= @_;

  # my($spanstart, $spanstop) = ($rb + $off1, $re + $off2); # add optional expansion
  if($spanstart<1) { $spanstart= 1; }
  
  my $regionid = "$ref:$spanstart-$spanstop"; # better
  $regionid.=":$spanor" if($spanor); #? $stranded and  OR always add

  if($spanstart > $spanstop) { ## error
    print $outh join("\t","#errorspan","na","na",$regionid,"badspan"),"\n";
    return;
  }

  ## spanor == strand should be linked to gene not region : new column or add to $id ? or aalen? or col1?
  if($doest and $spanor eq "-") { #? need to do the spanaa/cdna revcomp here?
    foreach (@$spanaa) {  $_=reverse($_); $_=~tr/ACGTacgt/TGCAtgca/; }
  }
  
  my $genedna  = get_dna( $genomefasta, $ref, $spanstart, $spanstop);
  
  my $id = join(",", @$spanid);
  my $geneaa = join(",", @$spanaa);
  my $aalen= length($geneaa) - (scalar(@$spanaa) - 1); 
  
  if($doest and $aalen < $MIN_ESTSPAN) {
    print $outh join("\t","#errorspan",$id,"tiny",$regionid,"na"),"\n";
    return; ## no note?
  }
  
  if($genedna =~ /\w/ and $geneaa =~ /\w/) {
    # need is to split out aa,dna for each to exonerate query
    my $glen= length( $genedna );
    #? add exonr type info?  est2genome vs protein2genome
    print $outh join("\t","exonerate-$xtype", $id, $geneaa, $regionid, $genedna, $aalen, $glen),"\n";
    
  } else {
    print $outh join("\t","#errorspan",$id,($geneaa?"okaa":"missaa"),$regionid,($genedna?"okna":"missna")),"\n"
    # write some error to outh
  }

}


sub get_dna {
  my($fasta, $ref, $start, $stop, $strand)= @_; 
  unless( $ref && $stop>0) {  warn "need ref-ID, start, stop location\n"; return; }
  unless( ref $genome_db ) {
    my $db = eval { Bio::DB::Fasta->new($fasta); }
      or die "$@\nCan't open sequence file(s). "; # and return;
    $genome_db= $db;  
    }
  
  my $seq = $genome_db->seq($ref, $start => $stop) 
      or return; ## die "cant locate seq of $ref:$start-$stop\n";#? and return;
  # $seq = $seq->revcomp() if($strand eq "-");
  $seq= $seq->seq if(ref $seq); # is this weird bioperl change here or not
  return $seq;
}

sub get_protein {
  my($fasta, $ref)= @_; 
  unless( $ref ) { warn "need ref-ID\n"; return; }
  unless(ref $protein_db) {
    my $db = eval { Bio::DB::Fasta->new($fasta); }
      or die "$@\nCan't open sequence file(s). "; # and return;
    $protein_db= $db;  
    }
  my $seq = $protein_db->seq($ref) or return;  
  $seq= $seq->seq if(ref $seq); # is this weird bioperl change here or not
  return $seq;
}

__END__


=item split to exonerate inputs

  nasvit1-protexonregion.tab n=21153
  
  cat nasvit1-protexonregion.tab | grep -v '^#' | head -10 | perl -ne\
  'chomp; ($exapp,$aid,$aa,$nid,$na)=split"\t"; $qt++;\
  puto("xr$qt.aa",$aid,$aa); puto("xr$qt.na",$nid,$na); \
  print "exonerate --query xr$qt.aa --target xr$qt.na \> xr$qt.exonr.gff\n";\
  sub puto{ my($a,$b,$c)=@_; open(F,">$a"); print F ">$b\n$c\n"; close(F); }'\
  > exoner.list


=item exonerate_script
  
  $exbin/exonerate \
    --model protein2genome \
    --refine region  --refineboundary 500 \
    --minintron 20 --maxintron 10000 \
    --showtargetgff --showvulgar 0 --showalignment 0 \
    --ryo '#qi %qi length=%ql alnlen=%qal\\n#ti %ti length=%tl alnlen=%tal\\n' \
    --query xr10.aa --target xr10.na > xr10.exonr.gff
  
  ##  --query $proteins  --target $genome  > $output_file_name.gff
  
  # fixme, other bindir here = $HOME/bio/augustus/scripts
  cat  $output_file_name.gff | perl -pe's,\\n,\n,g;' |\
  $evigene/scripts/process_exonerate_gff3.perl -sou exonr > $output_file_name.gff3
  
=item exonrcluster.sh : do both above for ncpu

  # /bin/bash
  ### env exontab=nasvit1-protexonregion.tab qsub -q normal exonrclust.sh
  
  ncpu=32
  
  exbin=$HOME/bio/exonerate/bin
  workd=$HOME/scratch/chrs/nasv1
  
  exonrcmd=$exbin/exonerate \
    --model protein2genome \
    --refine region  --refineboundary 500 \
    --minintron 20 --maxintron 10000 \
    --showtargetgff --showvulgar 0 --showalignment 0 \
    --ryo '#qi %qi length=%ql alnlen=%qal\\n#ti %ti length=%tl alnlen=%tal\\n' \
  
  cd $workd/prot/
  
  if ! test -f $exontab ; then  echo "missing input $exontab"; exit;  fi
  nam=`echo $exontab | sed 's/\..*//'`
  
  
  i=0; while [ $i != $ncpu ]; do 
  {
    out=$nam.exonr$i.gff ; touch $out
    $bindir/exonrclust.pl -exoner=$exonrcmd -in=$exontab -out=$out -i=$i -n=$ncpu >> $out &
    i=$(( $i + 1 ))
  }
  wait

#...
  #!/usr/bin/env perl
  # exonrclust.pl
  use strict;
  use Getopt::Long;
  
  my %args=();
  my $ok= &GetOptions( \%args,
    "exonerate=s" , "input=s" , "output=s" ,
    "i|islice=i" , "n|ncpu=i" , );
  
  my $exonrcmd="exonerate --model protein2genome "
  ." --refine region --refineboundary 500  --minintron 20 --maxintron 10000 "
  ." --showtargetgff --showvulgar 0 --showalignment 0 "
  ." --ryo '#qi %qi length=%ql alnlen=%qal\\n#ti %ti length=%tl alnlen=%tal\\n'";
  
  my $exonerate = $args{exonerate} || $exonrcmd;
  my $input = $args{input} || shift(@ARGV);
  my $output= $args{output} || "";
  my $ncpu= $args{n} || 1;
  my $icpu= $args{i} || 0;
  my $q=0;
  
  sub puto{ my($a,$b,$c)=@_; open(F,">$a"); print F ">$b\n$c\n"; close(F); } 
  
  die "bad opt or input " unless($ok and $input);
  my $inh=undef;
  if($input =~ /stdin|-/) { $inh=*STDIN; }
  else { open($inh, $input) or die "cant read $input"; }
  
  system("touch $output") if $output;
  while(<$inh>) {
    next unless(/^\w/ and ($q % $ncpu) == $icpu); 
    chomp; my($exp,$aid,$aa,$nid,$na)=split"\t"; 
    next unless($aid =~ /\w/ and $aa =~ /\w/ and $nid =~ /\w/ and $na =~ /\w/); #err?
    my @tf=("xr$q.aa","xr$q.na"); $q++;
    puto($tf[0],$aid,$aa); puto($tf[1],$nid,$na);
    my $cmd= "$exonerate --query $tf[0] --target $tf[1]"; 
    $cmd .= " >> $output" if $output;
    my $ok= system($cmd."\n"); unlink($tf[0]); unlink($tf[1]);
  } close($inh);


=cut

=item  use with expressed tile data to make exonerate est models

1. tile expression gff reformat:
# tile input needs to be per-group, then location-sorted:
# ** ?? should spans be clipped +/- n bp? n=9..18? depending on span

cat tst.hints | env minv=2 cut2=15 perl -ne \
'BEGIN{ $minv=$ENV{minv}||3; $cut1=$ENV{cut1}||6; $cut2=$ENV{cut2}||15; } \
($r,$src,$b,$e,$v)=(split)[0,1,3,4,5]; ($g)=m/grp=(\w+)/; $g||=$src; $v=int(2*$v); if($v < $minv){ next;} \
elsif($r eq $lr and $b<=$le and $e>$lb) { $b=$lb; $e=$le if($le > $e); $v=$lv if($lv>$v); } \
else { putg() if($le); } ($lr,$lb,$le,$lv,$lg)=($r,$b,$e,$v,$g); END{putg(); } \
sub putg{ my $fh; unless( $fh=$fh{$lg}){ open($fh, ">ti.$lg.gff") or die "ti.$lg.gff"; $fh{$lg}= $fh;\
warn "#write ti.$lg\n"; } my $ow= ($le - $lb > 70) ? $cut2 : $cut1; \
print $fh join("\t",$lr,$lg,"exon",$lb+$ow,$le-$ow,$lv,".",".",""),"\n"; }' \

2a.  exonregion to make genes tile-est x genome table
2b.  exonrclust to turn those into exonerate gene models

2. Notes: 
   - this requires input of mRNA spans from other source to
   delimit tile.gff exons into gene-ests, using exonr -model est2genome
   including -strand annotation file from mRNA strand.

    - test exonr --model cdna2genome 
      may or may not work well.  If input mRNA are mostly true, can use CDS bounds
      in annotation file to improve cdna2genome
      
    - test --splice3 splice3.mat --splice5 splice5.mat to improve --forcegtag 1

    
  set tigs=(adultf adultm pupaf pupam tesovao wingf wingm)
  
  foreach tig ($tigs)
    echo $tig
    
    cat nvit_epi4.sc25.mrna   ti.$tig.gff | sort -k1,1 -k4,4n -k5,5nr |\
    $evigene/scripts/exonregion.pl -v -sorted -strand -offset=300 \
        -xtype EST -exontype=exon -mrnatype=mRNA \
        -gff stdin -genome $nas11/genome/nasvit1asm.fa -out exonr.$tig.tab  
  
    cat exonr.$tig.tab | $evigene/scripts/exonrclust.pl -debug  -strand \
        -exonerate=EST -in=stdin -out=exonr.$tig.6.gff -i=0 -n=1 >& log.exonr6.$tig
  
  end

3. reprocess exonerate.gff

  here is made exon.ggb for test viewing
  # >> process_exonerate_gff3 
  #     -minalign=10  default=20 means skipping partial but good aligns?
  #     -dupalign=95 is this ok?
  
  touch xrtiles.ggb
  foreach tig ($tigs)
    echo $tig
    cat exonr.$tig.5.gff | $evigene/scripts/process_exonerate_gff3.perl \
    -refpart -t EST -minalign=10 -keep 'gene,exon,intron' -sou x$tig | \
   grep 'exon' | env s=$tig $evigene/scripts/gff2ggb.pl >> xrtiles.ggb
  end




=cut