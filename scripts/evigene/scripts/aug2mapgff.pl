#!/usr/bin/env perl
# aug2mapgff.pl

=item notes

  update 2013.jul: add -cleanid, always using this post-gff.
  update 2017.mar: add fasta header info on genes

=cut

use strict;
use Getopt::Long; # qw(:config no_ignore_case bundling);

my $exp= $ENV{n} || "0";
my $idprefix= $ENV{idprefix} || ''; ## || "AUG${exp}_";
my $noconstraint= defined $ENV{con} ? $ENV{con} : 0;
my $nofasta= defined $ENV{fa} ? $ENV{fa} : 0;
my $aa_attribute=0;
my $uniqueid= 1; #?
my $cleanid= 0; #? 1 default?
my $fastaonly=0;
my $cdsonly=0;
my $offset=0;
my $dropft="gene|intron|start_codon|stop_codon|transcription_start_site|transcription_end_site";
my $keepft="exon";

my $optok= &GetOptions (
    'source|exp=s' => \$exp,  
    'idprefix=s' => \$idprefix,
    'noconstraint' => \$noconstraint, #dgg
    'nofasta' => \$nofasta, ## make default, use withaa vs aa_attribute
    'fastaonly' => \$fastaonly,
    'cdsonly' => \$cdsonly,
    'offset=s' => \$offset, # for parted genome
    'keepft=s' => \$keepft, 
    'dropft=s' => \$dropft, 
    'aa_attribute!' => \$aa_attribute,
    'uniqueid!' => \$uniqueid,
    'cleanid!' => \$cleanid,
    );
    
die "usage: aug2mapgff  < augustus.gff > map.gff 
  options: --exp|source=expname --idprefix=AUGexp --nouniqueid --nocleanid
      --aa_attribute --nofasta --fastaonly --cdsonly --noconstraint
      --offset=120000  --keepft='intron|gene'
  protein aa will be added at end as ##FASTA unless -aa_att or -nofasta; -fastaonly outputs protein.fasta
  example:  -exp=1 -aa -nocon
SCAFFOLD33      AUGUSTUS.1      mRNA    33647   34810   0.68    +       .       ID=AUG1_g2.t1;evd_pG=LOC10
0115055,FGPP_ORF_2_from_SCAFFOLD33;evd_pP=Tcas:GLEAN_15479;protein=MSCNLCKCSSDGNYAACTFMQCFDFNFEEEQRSKRSTNE
VVAKLSSDIPRISGYTQGDQCPSKSFYNDCNMCVCGPDDASAACTMMMCMPGETQQPSKIVPAKLNDIARIGPGMPRRRVLPR
SCAFFOLD33      AUGUSTUS.1      exon    33647   33656   .       +       .       Parent=AUG1_g2.t1
SCAFFOLD33      AUGUSTUS.1      exon    33734   33769   .       +       .       Parent=AUG1_g2.t1
SCAFFOLD33      AUGUSTUS.1      exon    33865   34005   .       +       .       Parent=AUG1_g2.t1
SCAFFOLD33      AUGUSTUS.1      CDS     33912   34005   0.72    +       0       Parent=AUG1_g2.t1
SCAFFOLD33      AUGUSTUS.1      CDS     34116   34175   1       +       2       Parent=AUG1_g2.t1
\n" unless($optok);

my @drops= grep { $_ !~ m/$keepft/} split(/[,\| ]/,$dropft);
my $drops= join("|", @drops) or "NADA";

my (@gff,@attr,%hints,%prots,%cds,%fattr,%sups,%seengene, %nref, %ngene, %maxgene, %warned);
my ($aa,$cds,$spt,$psup,$supexons,$supintrons,$suputr5,$suputr3,$hintfull,$hintpart);
my ($ingene, $incon, $trid, $oldid, $nwarn,
    $inpart, $inpartloc, $inpartfile, $inaugfile) = (0) x 20;

$idprefix= "AUG${exp}_" unless($idprefix);
my $SRCTAG= ($exp =~ /AUG/) ? $exp : "AUG${exp}"; #<< make default# ($cleanid)?"AUG${exp}":"AUGUSTUS.${exp}"; 

#o $nofasta=1 if($aa_attribute);
## upd1703: cdsonly should set fastaonly
if($aa_attribute) { $nofasta=1; } elsif($cdsonly) { $fastaonly=1; }
if($fastaonly) { $nofasta=0; $aa_attribute=0; $noconstraint=1; }

=item cleanid

	update 2013.jul: add -cleanid, always using this post-gff.
	ID=p1.Scaffold0.g1.t1;Parent=p1.Scaffold0.g1 << fix messy IDs, from -uniqid
	
	cat $runset-augrun.gff | \
	$evigene/scripts/aug2mapgff.pl -aa -exp=$runt | perl -pe\
	'{s/[_\.](scaffold|super)[_]?/s/i; s/[_\.]Chr/s/; s/[_\.]p/p/; s/[_\.]g/g/; s/[_\.]t/t/; s/AUGUSTUS./AUG/;}' \
	> $runset-augmap.gff

=cut

## change \. to _ or drop; . causes problems
#old# sub cleanid { my $id=shift; $id =~ s/\./_/g; return($id); } 
sub cleanid {
	local $_=shift;
	s/[_\.](scaffold|super)[_]?/s/i; s/[_\.]Chr/s/i; # scaffold/chr ..
	s/[_\.](contig)[_]?/cn/i; #  add contig? 
	s/[_\.]p/p/; # partnum
	s/[_\.]g/g/; # genenum
	s/[_\.]t/t/; # trnum
	s/\./_/g; # any other dots
	return $_;
}

sub _idparts{ my($id)=@_;
 my(@ip); #parts:($idpre,$as,$ap,$ag,$at);
 if( @ip = $id =~ m/^(\w+)(?:s|cn)(\d+)p(\d+)g(\d+)t(\d+)/) { return @ip; }
 if( @ip = $id =~ m/^(\w+)(?:s|cn)(\d+)g(\d+)t(\d+)/) { splice(@ip,2,0,0); return @ip; }
 return($id,0,0,0,0);
}

sub _idsort { 
 my($apre,$as,$ap,$ag,$at)= _idparts($a);
 my($bpre,$bs,$bp,$bg,$bt)= _idparts($b);
 return ($as <=> $bs or $ap <=> $bp or $ag <=> $bg or $at <=> $bt or $a cmp $b);
}

sub put_transcript {

  if($incon) { # not transcript
    @gff=() if($noconstraint);
    
  } else {
  
    push(@attr,"pct_support=$psup") if(defined $psup);
    foreach my $k (sort keys %sups) {
      my $v= $sups{$k}; push(@attr, "sup_${k}=$v"); #?? want this
      }
    foreach my $k (sort keys %hints) {
      my $v= $hints{$k}; push(@attr, "evd_${k}=$v");
      }
    
    if($aa)  { $prots{$trid}= $aa ;  $aa=""; }
    if($cds) { $cds{$trid}= $cds ;  $cds=""; }
    
    # addattr local
    my($trloc,$tlen,$cdlen,$nx,$ncx)=("",0,0,0);
    # my($mrna)= grep /\tmRNA/, @gff;
    # if($mrna){ my($c,$b,$e,$o)= (split"\t",$mrna)[0,3,4,6]; $trloc="$c:$b-$e:$o"; }
    foreach (@gff) { 
      my($c,$t,$b,$e,$o)=(split)[0,2,3,4,6]; 
      if($t eq 'mRNA') { $trloc="$c:$b-$e:$o"; }
      elsif($t eq 'exon') { $tlen += 1+$e-$b; $nx++; } 
      elsif($t eq 'CDS') { $cdlen += 1+$e-$b; $ncx++; } 
    } 
    my $aaseq= $prots{$trid}||""; 
    my $aalen= length($aaseq)||0;
    $tlen=1+300*$aalen if($tlen<1); # missing exons?
  
    my $pcds= ($tlen<1) ? 0 :int(300*$aalen/$tlen);
    ## dont have '*' stopcodon by choice/default? assume or check cds for stop codons?
    #x my $ac=($aaseq=~/^M/ and $aaseq=~/\*$/)?"complete":"partial";
    my $hasstop=($cdlen > 3*$aalen)?1:0;
    my $ac= ($aaseq=~/^M/)?1:0; $ac += 2 if($hasstop);
    $ac= ($ac==3) ? "complete": ($ac==2)? "partial5": ($ac==1)?"partial3":"partial";
    push(@attr,"aalen=$aalen,$pcds\%,$ac"); # always 
    push(@attr,"protein=".$aaseq) if($aaseq and $aa_attribute); 
    push(@attr,"clen=$tlen");
    push(@attr,"exons=$nx,${ncx}cds");
    
    my $attr= join(";",@attr);
    foreach (@gff) { s/$/;$attr/ if(/\tmRNA/); } 
  
    ## subset of attr for put_fasta hdrs: location, aalen=, trlen=, pct_support=??
    my @fat= grep /aalen|clen|exons/, @attr;
    unshift(@fat, "locus=$trloc");
    $fattr{$trid}= join("; ",@fat);
  }
  
  unless($fastaonly) {
  if($cleanid){ map{ while(/(ID|Parent)=([^;\s]+)/g){ my($k,$d)=($1,$2); my $c=cleanid($d); s/$k=$d/$k=$c/;} } @gff; }
  print @gff;
  }
 
  @gff=(); @attr=(); %hints=(); %sups=();
  ($aa,$cds,$spt,$psup,$supexons,$supintrons,$suputr5,$suputr3,$hintfull,$hintpart)= (0) x 20;
}

sub put_fasta {
  my($seqh)= @_;
  print "#\n##FASTA\n"  unless($fastaonly);
  foreach my $id (sort _idsort keys %$seqh) {
    my $fa= $seqh->{$id};
    my $fhdr= $fattr{$id}||""; $fhdr=" $fhdr" if($fhdr);
    $id= cleanid($id) if($cleanid);
    $fa =~ s/(.{1,60})/$1\n/g; $fa=~s/\n$//;
    print ">$id$fhdr\n$fa\n";
    }
}

sub put_aafasta { put_fasta(\%prots); }
sub put_cdsfasta { put_fasta(\%cds); }

sub put_warnings {
  return unless(%warned); #? or do even if no warnings?
  my $errh= ($fastaonly) ? *STDERR : *STDOUT;
  print $errh "# WARNINGS gene data discrepencies:\n";
  print $errh "# ",join("\t",qw(iRef nFound Models Warns Scaffolds)),"\n";
  my($sngene,$smaxgene,$swarns)=(0) x 10;
  foreach my $iref (sort keys %ngene) {  # ref now == augfilenum
    my $refs    = $nref{$iref}     || 0;
    my $warns   = $warned{$iref}   || 0;   $swarns += $warns;
    my $ngene   = $ngene{$iref}    || 0;   $sngene+= $ngene;
    my $maxgene = $maxgene{$iref}  || 0;   $smaxgene+= $maxgene;
    print $errh "# ",join("\t",$iref,$ngene,$maxgene,$warns,$refs),"\n" if($warns or $ngene != $maxgene);
  }
  print $errh "# ",join("\t","TOTAL",$sngene,$smaxgene,$swarns,""),"\n";
}


#      G:  19 (GLEAN_18911,GH_DGIL_SNO_28281162,GH_NCBI_GNO_32001285,TRdgri_22378,GH_BREN_NSC_50051886,...)
#      P:   3 (CG31663_G1,CG31663_G2,CG7295_G1)
# hint groups fully obeyed: 7
#      T:   7 (tar16108,tar16109,tar16110,tar16111,tar16112,tar16113,tar16114)
# incompatible hint groups: 1
#      T:   1 (tar16115)

sub hintids {
  my($h)= @_;
  my($t, $n)= $h =~ m/^#\s+(\w+):\s+(\d+)/;
  my($v)= $h =~ m/\s+\(([^\)]+)\)/;
  $v =~ s/,\.\.\.//;
  #my @v=split ",",$v;
  return ($t,$n,$v);
}

sub end_gene {
  put_transcript() if(@gff);
  $incon= $ingene=0;
}  


while(<>){

    ## outputs may be mangled .. check here. esp check gene count
  if(/^#/) {
        
    if(/### gene (\S+)/) { end_gene(); $ingene=$1;  $incon=0; }
    elsif(/# start gene (\S+)/) { end_gene();  $ingene=$1;  $incon=0; } # new aug210 format
    elsif(/# Constraints/) { $incon=1; } # gone
    elsif(/# end gene/) {  end_gene(); }
    elsif(/# This output was generated with/){ print unless($fastaonly); }

    elsif(/# Evidence for and against this transcript/){ 
      $prots{$trid}= $aa if($aa); $aa="";
      $cds{$trid}= $cds if($cds); $cds="";
      }

##gff-version 3
#part_location: SCAFFOLD14:1950001-3535154
#part_file: 2 genoparts/SCAFFOLD14/SCAFFOLD14_1950001-3535154/nvit_epi4a-augrun.gff
    elsif(/##gff-version/){ print unless($fastaonly);  $inaugfile++; $inpart=0; $inpartloc= $inpartfile=""; }
    elsif(/#part_location: (\S+)/) { $inpartloc=$1; } 
    elsif(/#part_file:\s+(\d+)\s+(\S+)/) { $inpart=$1; $inpartfile=$2; } 

# ----- prediction on sequence number 1 (length = 2000000, name = SCAFFOLD1) -----
## add this for combo files?

    elsif(/transcript supported by hints .any source.: (\S+)/){ $psup=$1; $spt=""; }
    ## add these to output and counts per hint type X: n
    elsif(/CDS exons: (\S+)/){ $spt="cdsx"; $sups{"$spt"}=$1; $supexons=$1; }
    elsif(/CDS introns: (\S+)/){ $spt="cdsi"; $sups{"$spt"}=$1; $supintrons=$1; }
    elsif(/5.UTR exons and introns: (\S+)/){  $spt="utr5"; $sups{"$spt"}=$1;  $suputr5=$1; }
    elsif(/3.UTR exons and introns: (\S+)/){  $spt="utr3"; $sups{"$spt"}=$1; $suputr3=$1; }
    ### need to collect type counts for each sup type?
    elsif(/hint groups fully obeyed: (\S+)/){ $hintfull=$1; $hintpart=0; $spt=""; }
    elsif(/incompatible hint groups: (\S+)/){ $hintpart=$1; $hintfull=0; $spt=""; }
    
    elsif(($hintfull or $hintpart) and m/^#\s+(\w+):\s+(\d+)/) {
      my($t,$n,$v)= hintids($_);
      my $c= $hintfull ? "f" : "p";
      $hints{$c.$t}= $v || $n;  #? should this be $n/$v ?
    }

    elsif( $spt and m/^#\s+(\w+):\s+(\d+)/) {
      $sups{"$spt.$1"}=$2;
    }

    elsif(/# protein sequence = .(\w+)/) { $aa=$1; }
    elsif($aa and /# ([A-Z]+)/) { $aa.=$1; } # last patt?
    elsif(/# coding sequence = .(\w+)/) { $cds=$1; }
    elsif($cds and /# ([a-z]+)/i) { $cds.=$1; } # last patt?
    
  } elsif(/^\w/){
  
    my $keep= 1;
    
    ## trap error messages: Error in addGene; mangled outputs ..
    # start gene SCAFFOLD1.g3
    # SCAFFOLD1       AUGUSTUS        gene    161941  227674  0.02    +       .       ID=p1.SCAFFOLD1.g3
    # .. problem w/ ref.part counts ..
#part_location: SCAFFOLD1:1-2000000
#part_file: 1 genoparts/SCAFFOLD1/SCAFFOLD1_1-2000000/nvit_epit5_xmale-augrun.gff

    # ** OOOOh, this is mostly on purpose: recombine_partial_outputs.pl is dropping OVERLAP sections,
    # .. occurs at start of .p2, .p3, etc, 2nd part files where they overlap prior.
    # .. how to tell on purpose drops and data mangling?
    
    ##need on input augrun.gff, but where to avoid foulup?: 
    #if($cleanid){ s/[_\.](scaffold|super)[_]?/s/i; s/[_\.]Chr/s/; s/[_\.]p/p/; s/[_\.]g/g/; s/[_\.]t/t/; s/AUGUSTUS./AUG/; }
    
    if(/^(\S+)\tAUGUSTUS\tgene\t.*\tID=([^;\s]+)/) {     
      my $inref=$1; my $incheck=$2;
      if($inpart and $incheck =~ /^p$inpart./ and $ingene and $ingene !~ /^p/) { $ingene="p$inpart.$ingene"; }
      if($inpart) { $inref="p$inpart.$inref"; }
      if($ingene ne $incheck) { end_gene(); $ingene= $incheck; }
      
      ## combo file numbering: aug numbers genes from 1-end of combo, not per scaffold? YES
      ## need to use aug-file num, not inref, to match gnum
      
      #my $ngene= ++$ngene{$inref}; 
      $nref{$inaugfile} .= "$inref," unless($nref{$inaugfile} and $nref{$inaugfile} =~ /$inref/);
      my $ngene= ++$ngene{$inaugfile}; #??
      my ($gnum)= $incheck =~ m/(\d+)$/;   
      $maxgene{$inaugfile}= $gnum if(not $maxgene{$inaugfile} or $maxgene{$inaugfile} < $gnum);
      
      # do this warn 1/inref, and summarize at end? so know which scaffolds have bad gene counts
      if($ngene != $gnum and not $warned{$inaugfile}++)   ## bad check for ref.parts
      { warn "# Warning: gene count off, $inref:$ngene not $incheck\n" unless($nwarn++ > 9); }
    }
    unless($ingene or $incon) { warn "# Warninng: not gene, mangled input? $_" unless( $nwarn++ > 9); } #error
    
    if($incon) {  # this is obsolete w/ new augustus v
      s/grp=/ID=/;
      s/\tep/\texon_region/;  # should use SO terms : region
      s/\tip/\tintron_region/;
      s/\tirpart/\tintergenic_region/;
      s/ "[^"]+"//;
      $keep=0 if($noconstraint);
      
    } elsif($ingene) {
      s/\tAUGUSTUS\t/\t$SRCTAG\t/; # was s/x/AUGUSTUS.$exp/; # see above cleanid shortens to AUG$exp
      my @cols= split"\t";
      my $typ= $cols[2];
      
      if($typ eq "CDS"){ s/ID=[^;\s]+;?//; } # only extra ID on CDS? exon,.. ?
      ## if(type !~ /transcript|gene/) {  s/ID=[^;\s]+;?//; } 
      
      my $prefix= $idprefix;
      if($uniqueid) { # check for/add ref/scaffold/chr num prefix?
        my($ref)= m/^(\w+)/;
        (my $sref=$ref) =~ s/^(\w)\D+/$1/; # Scaffold123 -> S123
        my($id)= m/(?:Parent|ID)=([^;\s]+)/;
        unless($id =~ m/$ref/ or $id =~ m/$sref/ or $prefix =~ m/$sref/) {
          #$ref =~ s/^(\w)(\D+)/$1/;
          $prefix= $idprefix.$sref;
          }
        }
      s/(ID|Parent)=/$1=$prefix/g if($prefix);
      
      if($typ eq "transcript") {
        put_transcript() if(@gff); ## end_gene(); $ingene=$insav; 
        s/\ttranscript\t/\tmRNA\t/;
        s/;Parent=[^;\s]+//; 
        ($trid)= m/ID=([^;\s]+)/; $oldid= $trid;
        $trid= cleanid($trid); #? here or on output
        if($uniqueid && $seengene{$trid}++) { 
          my $did=  "d" . $seengene{$oldid};
          unless($trid =~ s/$prefix/$prefix$did/) {
            $trid = $did.$trid;
          }
          #s/=$oldid\b/=$trid/g;
        }
        s/=$oldid\b/=$trid/g unless($trid eq $oldid);
        
      } else {
        if($oldid && $trid ne $oldid) {
          s/=$oldid\b/=$trid/g;
          }
      }
        
      $keep= ($typ =~ m/^($drops)/) ? 0:1; #? change to ($typ ~= m/^($keeps)/)?1:0;
    }
      
    $keep=0 if(/^Error/);
    if($offset>0) {  # for parts collating
       my @v = split"\t";
       my($b,$e)=@v[3,4]; if($e>0) { $v[3]= $offset+$b; $v[4]= $offset+$e; }
       $_= join"\t", @v;
     }

    push(@gff,$_) if $keep; 
  }

}

put_transcript() if(@gff);
put_aafasta()  unless($nofasta or $cdsonly); # this is GFF + FASTA, dont do unless requested..
put_cdsfasta() if($cdsonly);
put_warnings();  # if(%warned);

__END__

## add in hints used comments for evidence Dbxref attribs
# Evidence for and against this transcript:
# % of transcript supported by hints (any source): 87.1
# CDS exons: 14/14
#      G:  14 
#      P:   6 
# CDS introns: 13/14
#      G:  13 
# 5'UTR exons and introns: 0/2
# 3'UTR exons and introns: 0/1
# hint groups fully obeyed: 2
#      G:   1 (FGPP_ORF_1_from_SCAFFOLD15)
#      P:   1 (Amel:GB12599-PA)
# incompatible hint groups: 4
#      G:   1 (276244)
#      P:   3 (Amel:GB11078-PA,Human:NP_478144,Tcas:GLEAN_02629)


### gene g316
scaffold_4      AUGUSTUS        gene    1420609 1422424 0.21    +       .       ID=g316
scaffold_4      AUGUSTUS        transcript      1420609 1422424 0.21    +       .       ID=g316.t1;Parent=g316
scaffold_4      AUGUSTUS        exon    1420609 1422424 .       +       .       Parent=g316.t1
scaffold_4      AUGUSTUS        start_codon     1421469 1421471 .       +       0       Parent=g316.t1
scaffold_4      AUGUSTUS        CDS     1421469 1421723 0.71    +       0       ID=g316.t1.cds;Parent=g316.t1
# protein sequence = [MKELMKQLPVVKPLSTKEPAEKIALFQCRLLLPKSSPLQGEVVGDPMLTKRLAKRAAALRACERLHQLKELDDLHLLP
# VSHRKR]
# Evidence for and against this transcript:
# % of transcript supported by hints (any source): 100
# CDS exons: 1/1
#      T:   1 
# CDS introns: 0/0
# 5'UTR exons and introns: 1/1
#      T:   1 
# 3'UTR exons and introns: 1/1
#      T:   1 
# hint groups fully obeyed: 5
#      T:   5 (tar15253,tar15254,tar15255,tar15256,tar15257)
# incompatible hint groups: 0
### end gene g316
