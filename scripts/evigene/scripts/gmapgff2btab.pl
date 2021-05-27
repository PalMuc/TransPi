#!/usr/bin/env perl
# gmapgff2btab.pl

=item about gmapgff2btab.pl

This reads exons from GMAP produced GFF EST alignments, and writes 
BTAB alignment format for input to PASA.

  gmapgff2btab.pl  mydata.gff  mydata.fasta
  output: mydata.btab, mydata.btab.fa
  
Could instead use Launch_PASA_pipeline --IMPORT_CUSTOM_ALIGNMENTS_GFF3 
but that needs other changes to parse gmap.gff.

After genrating mydata.btab, and mydata.btab.fa, 
run PASA pipeline WITHOUT steps 2,3:
   
   step0. create_mysql_cdnaassembly_db
   step1. upload_transcript_data -t mydata.btab.fa
   
   skip step2. process_GMAP_alignments.pl (runs gmap est x genome)
   skip step3. gmap_to_btab.pl
   
 continue with 
   step4. populate_alignments_via_btab.dbi -b mydata.btab

fixed: 
  input gmap.gff has + genome but - target orient
  unusual, but this is correct: the target cDNA is antisense of real transcript
  and needs to be revcomp for use with PASA

todo: 
  handle duplicate alignments here? need that info on ESTs/rnaseq asm that
  match equally well to multi locations.  
  GMAP.gff marks these as ID=trID.mrna1,.mrna2,... and has percent align,ident on mRNA lines  
  For PASA, need new ID, and new copy of transcript/EST seq w/ same ID
  Should keep .mrna2..n where %align * %ident >= 0.99/0.98/.. of .mrna1

todo2: add in parts of rnagff2btab.sh

=cut


use strict;

my $debug=1;
my $RNAGFF2BTAB=  $ENV{rnagff2btab} ||0; # option
my $IN_BTAB= $ENV{btab} || 0;
my $SRCPROG="gmap";
my $ExonType='exon|HSP'; # need opt.
my $mRNAType='mRNA|match';
my $noTFLIP= $ENV{noflip} || 0; # ** NO this doesn't work right, but does work like PASA gmap_to_btab
my $revbtabfa= $ENV{dorev} || 0; # use with RNAGFF2BTAB
my $useTargSense= $ENV{sense} || 0;
my $orientonly= $ENV{oronly} || 0;
my $idtag  = $ENV{idtag} || 0;
my $scoreIsPI= $ENV{piscore} || 0;

my $gff   = shift @ARGV; #or stdin ; daphmag_velsubs512
my $fasta = shift @ARGV;
my ($inh, $ok, $nam, $lid, %idtor, %iddup, %idmap, $xi, $ci, $tw, $err, @exons, %facount);

$IN_BTAB=1 if($gff =~ /\.btab$/);

if($IN_BTAB) {  # write only btab.fa

  if($gff =~ m/^(stdin|\-)$/) { $inh= *STDIN; $nam="inputbtab"; }
  else { open(F,"$gff") or die "input $gff"; $inh=*F; ($nam = $gff) =~ s/\.btab//; }
  $fasta= "$nam.fa" unless($fasta);

  my($lid, $lb, $le)=(0,0,0);
  while(<F>) {
    my @v=split"\t"; 
    my $pid= $v[5];
    (my $tid=$pid) =~ s/\.mrna\d+$|_G\d+$|_C\d+$//;
    if($revbtabfa) {
      my $or = ($lid eq $pid and $v[6] < $lb) ? "-" : "+";
      $idtor{$tid}= $or;
    } else {
      $idtor{$tid}= "+";
    }
    $iddup{$tid}{$pid}=1 if($tid ne $pid);
    ($lid,$lb,$le)= @v[5,6,7];
  } close(F);

} elsif($RNAGFF2BTAB) {

  #above# my $scoreIsPI= $ENV{piscore} || 0;

  if($gff =~ m/^(stdin|\-)$/) { $inh= *STDIN; $nam="inputgff"; $ok=1; }
  elsif($gff =~ /\.gz$/){ open(F,"gunzip -c $gff|") or die "in $gff"; $inh=*F; ($nam= $gff) =~ s/\.gff.gz//; }
  else { open(F,"$gff") or die "input $gff"; $inh=*F; ($nam = $gff) =~ s/\.gff//; }

  open(OUT, ">$nam.btab") or die "writing $nam.btab";
  rnagff2btab($inh,$scoreIsPI);
  close(OUT);

## add also this bit for gff > est.fa
#   if( $maketr != 0 ) then
#     $workd/scripts/cdsgff2genbank.pl -MIN_CDS_LEN=1 -t=$feature -a dna \
#      -gff $trf -fasta $workd/genome/$dgenome.fa > $nam.fatmp
#     env dorev=1 $workd/scripts/gmapgff2btab.pl  $nam.btab $nam.fatmp
#   endif

  
} else {
  
  if($gff =~ m/^(stdin|\-)$/) { $inh= *STDIN; $nam="inputgff"; $ok=1; }
  elsif($gff =~ /\.gz$/){ open(F,"gunzip -c $gff|") or die "in $gff"; $inh=*F; ($nam= $gff) =~ s/\.gff.gz//; }
  else { open(F,"$gff") or die "input $gff"; $inh=*F; ($nam = $gff) =~ s/\.gff//; }
  open(OUT, ">$nam.btab") or die "writing $nam.btab";
  $fasta= "$nam.fa" unless($fasta);
  
  unless($noTFLIP) {
  my($h,$fa,$rev,$sid,%seqid);
  $ok=0;
  if($fasta =~ /\.gz$/){ $ok=open(FA, "gunzip -c $fasta |"); }
  else { $ok= open(FA, $fasta); }
  while(<FA>) {
    if(/^>(\S+)/) { my $sidt=$1; 
      $facount{$sid}= length($fa) if($fa);
      $h=$_; $fa=""; $sid= $sidt;
    } else { 
      chomp; $fa.=$_; 
    } 
  } close(FA);
  $facount{$sid}= length($fa) if($fa);
  }
  
  my($mrnaid,$tsense,$mrnaor)=(0,0,0);
  while(<$inh>) {
    next unless(/^\w/);
    my($r,$src,$ty,$b,$e,$v,$o)= split"\t";   
    
    if($ty =~ /^($mRNAType)$/) {
      ($mrnaid)= m/ID=([^;\s]+)/;
      if($idtag) { my($nid)= m/$idtag=([^;\s]+)/; $mrnaid=$nid if($nid); }
      $tsense= (m/;sense=-1/) ? "-" : "+"; # only found for sense=-1
      ###$tsense= (m/;sense=([^;\s]+)/) ? $1 : 0; # only found for sense=-1
    }
    next unless($ty =~ /^($ExonType)$/); #exon, HSP, CDS?
    next if($orientonly and $o eq ".");
    
    my($pid)= m/Parent=([^;\s]+)/;
    if($idtag) { my($nid)= m/$idtag=([^;\s]+)/;  $pid=$nid if($nid); }
    
    #fixme for Target lacking orient, but gff from gmap2gff with sense=- on match/mRNA line
    
    # my($tid,$tb,$te,$to)= m/Target=(\S+)\s(\d+)\s(\d+)\s(.)/;  
    # my($tid,$tb,$te,$to)= m/Target=(\S+)\s(\d+)\s(\d+)\s?(\S?)/;  
    my ($trg)= m/Target=([^;\n]+)/;
    my($tid,$tb,$te,$to)= split /[\s\+]/, $trg;
    if($te and !$to and $pid eq $mrnaid) { $to= $tsense; } ## ($tsense < 0) ? "-" : "+";
    ## and $useTargSense ?
    
    unless($te and $to) { 
      warn "No Target=ID start end orient for exon: $_";
      $err++; die if $err>10; next; 
      }
  
    $to="+" if $noTFLIP;
    
    my $ipath=1;
    my $id1= $tid;
    my $id= ($idtag and $pid) ? $pid : $tid;  # IdType == Target || Parent
    if($pid and $pid ne $tid) {
      $ipath=($pid =~ /(\.mrna|_G|_C)(\d+)$/) ? $2 : 1; 
      if($idtag) { ($id1=$id) =~ s/(\.mrna|_G|_C)(\d+)$//; }
    }  
    if($ipath>1) { $iddup{$id1}{$pid}=1; $id= $pid; $idtor{$id}= $to unless($idtor{$id1}); } 
    else { $idtor{$id}= $to; } # is idtor{} ok w/ dup ids?  
  
    if($lid ne $id) { putx($lid, @exons); @exons=(); $xi=0; $ci++; $tw=$te+1; } 
    $xi++; $lid=$id; 
    ($b,$e)=($e,$b) if($o eq "-"); 
    
    # if($to eq "-") { ($tb,$te)= ($tw-$te, $tw-$tb); }  ## for gmap.gff, not for gmap2gff ??
    ## this depends on sort order of exons; bad idea;
  
    push @exons, [  $r,"","",$SRCPROG,"",$id,$b,$e,$tb,$te,$v,"",$ipath,$ci,$xi, $to ];
    # print OUT join("\t",$r,"","",$SRCPROG,"",$id,$b,$e,$tb,$te,$v,"",$ipath,$ci,$xi),"\n"; 
  }
  putx($lid,@exons);
  close(OUT);
  
} # in_btab


$ok=0;
if($fasta =~ /\.gz$/){ $ok=open(FA, "gunzip -c $fasta |"); }
else { $ok= open(FA, $fasta); }

unless($ok) { 
  my $no= 0; map{ $no++ if $idtor{$_} eq "-"; } keys %idtor;
  if($no>0) {
    my $fidor="$nam.idor";
    warn "No sequence file $nam.fa to flip n=$no antisense target sequences.\nWriting id-orient table $fidor.";
    open( ID, ">$fidor") or die "writing $fidor";
    foreach my $id (sort keys %idtor) { print ID join("\t",$id, $idtor{$id}),"\n"; }
    close(ID);
    }
  exit;
  };

my($h,$fa,$rev,$sid,%seqid);
open(FOUT, ">$nam.btab.fa") or die "writing $nam.btab.fa";
while(<FA>) {
  if(/^>(\S+)/) { my $sidt=$1; 
    putfa($sid, $h,$fa,$rev) if($fa);
    $rev=($idtor{$sid} eq "-") ? 1:0;  
    $h=$_; $fa=""; $sid= $sidt;
  } else { 
    chomp; $fa.=$_; 
  } 
}
putfa($sid,$h,$fa,$rev) if($fa);
close(FOUT);

if(scalar(%seqid)) {
  my @tabno= grep{ not $seqid{$_} } sort keys %idtor; # @tabid;
  my $notab= scalar(@tabno);
  if($notab > 0) {
    my $outf="$nam.btab.nofa.ids"; 
    warn "More IDs in .btab than .btab.fa sequences\nWriting $notab extra IDs to $outf\n";
    open(OUT,">$outf"); print OUT join("\n",@tabno),"\n"; close(OUT);
  }
}

#-----------------------------------------------------------

# from epasa/rnagff2btab.sh

sub rbput { 
  my($lo, $piscore, $xa)= @_;
  my @x= @$xa;
  my $rev=($lo eq "-")?1:0; 
  my $x0=0; my $ix=0;   
  if($rev) { @x=sort{$b->[2]<=>$a->[2]} @x;} else { @x=sort{$a->[2]<=>$b->[2]} @x;}  
  foreach my $xr (@x) { my($r,$id,$b,$e,$v,$o,$w,$ci,$xi)=@$xr; 
    my $br=1 + $x0; my $er= $x0 + $w; $x0 += $w; $ix++;
    my $pi= $piscore ? $v :100; 
    print OUT join("\t",$r,"","",$SRCPROG,"",$id,$b,$e,$br,$er,$pi,"",1,$ci,$ix),"\n"; 
    } 
}
  
sub rnagff2btab {
  my($inh,$piscore)= @_;
  my($lid,$lo,$xi,$ci)= (0) x 10;
  my @x;

  # $CAT $trf | grep "    $feature" | env piscore=$SCOREISIDENT perl -ne 
  while(<$inh>) {
    next unless(/^\w/);
    my @v=split"\t";  
    my($r,$s,$t,$b,$e,$v,$o,$px,$at)= @v;
    next unless($t =~ /^($ExonType)$/); # 201110 fix
    next if($orientonly and $o eq ".");
    
    my($pid)= $at=~ m/Parent=([^;\s]+)/;
    if($idtag) { my($nid)= m/$idtag=([^;\s]+)/;  $pid=$nid if($nid); }
    (my $tid=$pid) =~ s/\.mrna\d+$|_G\d+$|_C\d+$//;
     # if($revbtabfa) .. my $or = ($lid eq $pid and $v[6] < $lb) ? "-" : "+";
    $idtor{$tid}= "+"; # unless($revbtabfa) ..
    $iddup{$tid}{$pid}=1 if($tid ne $pid);

    my $w=1+$e-$b; 
    ($b,$e)=($e,$b) if($o eq "-");  
    if($lid eq $pid) { $xi++; } else { rbput($lo,$piscore,\@x) if(@x); @x=(); $ci++; $xi=1; }  
    push(@x,[$r,$pid,$b,$e,$v,$o,$w,$ci,$xi]);  
    $lo=$o; $lid=$pid;  
    }
  rbput($lo,$piscore,\@x) if(@x);
}


sub putx {
  my ($sid, @xons)=@_ ;
  return unless(@xons);
  
  my $to= $xons[0]->[-1]; # == [15]
  my $tmax=0; 
  $tmax= $facount{$sid} || 0; $tmax++; # off by 1
  
  #?? maybe this is wrong .. or new KEEPantisense=1 gmap_to_gff/btab ??
  # ** problem here is that tmax is not length of transcript. align stops before end.
  # -- need true seq length to get this working
  if($to eq "-") { 
    # foreach my $x (@xons) { my($tb,$te)= @{$x}[8,9]; $tmax=$te if($te>$tmax); }
    # $tmax++;
    # @xons = sort{ ($tmax - $a->[9]) <=> ($tmax - $b->[9]) }  @xons;
  } else {
    # @xons = sort{  $a->[8] <=> $b->[8] }  @xons;  # no, want this for pasa compat.
  }

  #?? do we need to resort exons by tb, renumber xi? ONLY for to -
  # pasa ERROR:  Incontiguous alignment now common ...
  my $ix= 0;
  foreach my $x (@xons) {
    my( $r,$x1,$x2,$x3,$x4,$id,$b,$e,$tb,$te,$v,$x5,$ipath,$ci,$xi, $to )= @$x;
    # $ix++; $xi= $ix;
    if($to eq "-") { ($tb,$te)= ($tmax-$te, $tmax-$tb); }  ## for gmap.gff, not for gmap2gff ??
    print OUT join("\t",$r,"","",$SRCPROG,"",$id,$b,$e,$tb,$te,$v,"",$ipath,$ci,$xi),"\n"; 
  }
}

sub putfa { 
  my($sid,$h,$fa)= @_;  # ,$rev

  # return unless($fa and $h and $idtor{$sid});
  return unless($fa and $h);
  my $ok=0; 
  unless( $ok= $idtor{$sid} ) { 
    if(ref $iddup{$sid}) {  $ok=1; } # more?
    }
  return unless($ok);
  
  my $rev=($idtor{$sid} eq "-") ? 1:0;  
  if($rev){ $fa=reverse($fa); $fa=~tr/acgtACGT/tgcaTGCA/; $h=~s/$/; rev=1/; } 
  $fa =~ s/(.{60})/$1\n/g; 
  print FOUT $h,$fa,"\n"; 
  $seqid{$sid}++;
  
  if($iddup{$sid}) {
    foreach my $di (sort keys %{$iddup{$sid}} ) {
    $h =~ s/>\S+/>$di/;
    print FOUT $h,$fa,"\n"; 
    $seqid{$di}++;
    }
   }
}
