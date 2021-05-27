#!/usr/bin/env perl
# trasmfixup.pl

## FIXed: add both strands.an.gff | overbest > bestofboth.gff
##   .. run as: cat  $name.{fwd,rev}.an.gff | trasmfixup -best -name $name -gff stdin ..

## FIXed: option for trsplitjoins, then remap trgmap, so input is only trasm.tr seq


use strict;
use Getopt::Long;

# .. CONFIG for trmap : fixme ..
my $gmap="/bio/bio-grid/mb/gmap/bin/gmap";
my $evigene="/bio/bio-grid/mb/evigene";
my $genomed="/bio/bio-grid/daphmag/genome/";
my $dgenome="dmagna20100422assembly";
my $gmapdb="$genomed/gmap";
my $introns="intronall12_good.gff"; ##  # FIXME ... caller Path Changes **
my $bestscores= "inqual:5,cov:8,pid:5,CDS:2,UTR:1,cdsindel:-2";
my $dropscores="";
my $SKIP2BEST=0;

my %args=();
my $optok= &GetOptions( \%args,
  "gff=s", "trfa|name=s", "strand=s", 
  "config=s", 
  "chimerapick!", "bestonly!", "splitjoins!", # nobest opt for skipping best output?
  # "output=s" , # use this for outputdir for all files created?
  # "logfile=s" , "pid=s", "i|part=i" , 
  "version=s", "debug:s", "nodebug" );

my $debug=$args{debug};   # input gmap.gff before strand,protein filters
my $gff=$args{gff};   # input gmap.gff before strand,protein filters
my $trfa=$args{trfa}; # or name, for -best
my $vers=$args{version} || 1;
my $strand= $args{strand} || ""; # allowed nostrand

my $CHIMERPICK= $args{chimerapick} || 0;
my $SPLITJOINS= $args{splitjoins} || 0;
my $NOBEST= 0;
$SKIP2BEST=1 if($debug =~ /best/);
if( defined $args{bestonly}) {
  $SKIP2BEST=$args{bestonly}; 
  $NOBEST= not $SKIP2BEST; ## not $args{bestonly};
}  

readconfig($args{config}) if($args{config});

sub strandt { my $s=shift; return ($s =~ /fwd|for|^f|\+/) ? 'fwd' : ($s =~ /rev|^r|\-/) ? 'rev' : 'noo'; }
sub strandc { my $s=shift; return ($s =~ /fwd|for|^f|\+/) ? '+' : ($s =~ /rev|^r|\-/) ? '-' : ''; }
my $strandc= strandc($strand);

die "usage: $0 -gff $gff -trfa $trfa -strand $strand\n" 
  unless($optok and 
    ($SKIP2BEST 
    or ($SPLITJOINS and -f $trfa)
    or (-f $gff and -f $trfa) )
  ); ##  and $strandc : skip strand filter; FIX: allow only -gff stdin for -best
  # CHIMERPICK wants both gff and trfa, precedes genefindcds
  
my ($bestgff) = trfixup( $gff, $trfa, $strandc);

# exit;

#......... subs ................
 
sub readconfig {
  my($cf)= @_;
  open(F, $cf) or die "Missing config: $cf";
  while(<F>) { 
    next unless(/^\w/); chomp;
    s/#.*$//;  s/\s*$//;   # s/=/ /;
    my($k,$v)= split /[=\s]+/,$_,2;
    CFG: for (lc($k)) {
      #/^velbin$/  && do { $velbin=$v; last CFG; };
      #/^kmers/  && do { $kset=$v; last CFG; };
      #/^optvelvet/ && do { $vopts=$v; last CFG; };
      #/^optoases/ && do { $oopts=$v; last CFG; };
      #/^samtools$/ && do { $samtools=$v; last CFG; };  
      /^gmap$/    && do { $gmap=$v; last CFG; };   # was bad for gmapdb xxx
      /^gmapdb$/  && do { $gmapdb=$v; last CFG; };
      /^evigene$/ && do { $evigene=$v; last CFG; };
      /^genomed$/ && do { $genomed=$v; last CFG; };
      /^dgenome$/ && do { $dgenome=$v; last CFG; };
      /^intron/   && do { $introns=$v; last CFG; };
      /^bestscore/  && do { $bestscores=$v; last CFG; };
      /^dropscore/  && do { $dropscores=$v; last CFG; };
      #/^log$/     && do { $LOGFILE=$v; last CFG; };
      #/^minlen$/  && do { $MINLEN=$v; last CFG; };
      ## more configs: input bam name patt >> velvet insertsize config
      ## pairinsize=1=grHS,2=grND   ... grND:380,grHS:475
      #/^pairinsize/ && do { $pairinsize=$v; last CFG; };
      
    }
  }
}

sub runcmd
{
  my (@cmd)= @_;  # options?  ## $dieOrNot
  my $cmd= join(" ",@cmd);
  warn "#CMD: ",$cmd,"\n" if($debug);

#  #.. sigh .. this works
#   open OLDOUT, ">&", \*STDOUT;
#   open OLDERR, ">&", \*STDERR;
#   open STDOUT, ">>$LOGFILE" ;
#   open STDERR, ">>$LOGFILE" ;
#   warn "\n#CMD: ",$cmd,"\n";

  my $err= system(@cmd);
  
#   open STDOUT, ">&OLDOUT";
#   open STDERR, ">&OLDERR";
  
  # FIXME here: dont die here when velvet pukes out on 1 kmer.. put this into LOGFILE ?
  if($err) { 
    warn "#CMD: ",$cmd,"\n" unless($debug);  
    warn "#FAIL: $cmd[0]: $?";
    ## die if($dieOrNot);
    }
  return $err;
}


sub gffsuf{ my($f,$s,$d)= @_; my $fn=$f; $fn=~s/\.gz//; $fn=~s/\.$d// if($d); 
  unless($fn =~ s/\.gff/.$s.gff/) { $fn.=".$s"; } return $fn; }


sub trsplitjoins   
{
  my($trfa)= @_;
  
  # FIXME: expected bug, some utrorf need further splitting, ie trasm was 3+ genes mashup, not just 2
  # .. correct cdna_bestorf.pl ; also change id, or transfer annot, of main cdna that is split
  
  my $trfasplit= $trfa.".orfsplit";
  my $cmd= "$evigene/scripts/cdna_bestorf.pl -splitutrorf -act cdnaonlyfasta -cdna $trfa > $trfasplit";
  my $err= runcmd($cmd);
  if($err) { warn "#ERROR: trsplitjoins($trfa): $err\n"; return ($trfa, $err); }
  else { return($trfasplit); }  #  push @cleanup, $trfa; 
  ## ? add log note: num in, num out=splits
}

sub trmap   
{
  my($trfa)= @_;
  my $cmd;
  my $name= $trfa; $name=~ s/\.(tr|fasta|fa)$//; 
  my $trgff= "$name.gff";
  
  $cmd="$gmap --npaths=0 --min-intronlength=30 -S -D $gmapdb -d $dgenome $trfa | ".
  "env src=$name noerrspan=1 intron=-1 best=0 $evigene/scripts/gmap_to_gff.pl > $trgff";
  my $err= runcmd($cmd); # runCmdOrDie($cmd);
  if($err) { warn "#ERROR: trmap($trfa): $err\n"; }
  return( $trgff, $err );
}


=item trchimerpick

daphmag/rnas/asmrna4/asmfull/chimbad.info
many trmap locations are splits of tr, often of genejoins, even after cdna_bestofr -splitutrorf.
.. often only 1 has useful mapping, 2nd is long utr w/o valid gene model, but maps back onto same/near locus
choose best and drop 2nd partial mapping: always? 
input is gmap.gff with chimera= annots. ?? do we need genefindcds annot for cdnabest? problematic now

  eg:
  grep ID=veln4ptr001k33Loc11t7 spltrasmn4r/trasmn4r.tr.split.1.fa.split.gff
  scaffold00512   trasmn4r.tr.split.1.fa.split    mRNA    344699  348848  74      -       .       ID=veln4ptr001k33Loc11t7cdna_C1;Target=veln4ptr001k33Loc11t7cdna 112 2625;
    aalen=265;cov=74.6;indels=3/0;match=2507;nexon=4;pid=99.7;qlen=3368;path=1/2;chimera=breakpoint at 2624..2626;chim2=scaffold00512:344559-345302:.;cdsindel=3
  scaffold00512   trasmn4r.tr.split.1.fa.split    mRNA    344559  345302  22      .       .       ID=veln4ptr001k33Loc11t7cdna_C2;Target=veln4ptr001k33Loc11t7cdna 2626 3368;
    aalen=21;cov=22.1;indels=0/1;match=743;nexon=1;pid=99.9;qlen=3368;path=2/2;chimera=breakpoint at 2624..2626;chim1=scaffold00512:344699-348848:-;cdsindel=152
  scaffold00512   trasmn4r.tr.split.1.fa.split    mRNA    342998  344074  36      .       .       ID=veln4ptr001k33Loc11t7utrorfcdna;Target=veln4ptr001k33Loc11t7utrorfcdna 552 1628;
    aalen=197;cov=66.2;match=589;nexon=1;pid=99.5;qlen=1628;path=1/2
  
=cut

sub trchimerpick
{
  my($trgff)= @_;
  
  my $trfixgff= gffsuf($trgff,"chifix","gmap"); # $trgff; $trgffan =~ s/.gff/.chifix.gff/;

  # need to read thru full gff to pick out chimer pairs; then reread/print, dropping bad cases.
  my (%chims); 
  my $err=0; my $nchi=0;
  if($trgff =~ m/\.gz/) { open(F,"gunzip -c $trgff |") or $err++; } else { open(F,"$trgff") or $err++; }
  while(<F>) {
    if(/^\w/ and /\tmRNA/ and /chimera=/) { 
      my @v=split"\t"; my $ga=$v[8]; chomp($ga);
      my($id)=m/ID=([^;\s]+)/; # expect _C[12] id suffix
      my($trid)=m/Target=([^;\s]+)/;  unless($trid) { ($trid= $id) =~ s/_C\d$//; }
      my($cn)= $id =~ m/_C(\d+)/; unless($cn) { ($cn)= $ga =~ m/path=(\d+)/; }
      my %ga= map{ my($k,$v)=split/[=,]/;  $k => $v; } split";",$ga;
      if($cn) { $chims{$trid}{$cn}= \%ga; $nchi++; }
    }
  } close(F);
  return($trgff) if($err or $nchi==0);

  ## FIXME: bug here? dropping some exons of kept mRNA : err=Missing-mrna-exon
  ## .. or below in strandc filtering.. all are chimera and '.' strand, prior were cdnabest w/ strand.
  ## ^^ yes here, chomp at mrna..
  
  my $keep=1; my $ndrop=0; 
  if($trgff =~ m/\.gz/) { open(F,"gunzip -c $trgff |") or $err++; } else { open(F,"$trgff") or $err++; }
  open(O,">$trfixgff") or $err++;
  while(<F>) {
    my $inrow=$_;
    if(/^\w/ and /\tmRNA/) { 
      $keep=1;
      if(/chimera=/) {
        # chomp; #BUG here for print
        my @v=split"\t"; my $ga=$v[8]; chomp($ga);
        my($id)=  m/ID=([^;\s]+)/; 
        my($trid)=m/Target=([^;\s]+)/; unless($trid) { ($trid= $id) =~ s/_C\d$//; }
        my($cn)= $id =~ m/_C(\d+)/; unless($cn) { ($cn)= $ga =~ m/path=(\d+)/; }
        if($cn) { 
          my $cb= ($cn==2) ? 1 : 2;
          my $ga= $chims{$trid}{$cn};
          my $gb= $chims{$trid}{$cb};
          if($ga and $gb) {
            my $gal= $ga->{aalen}; $gal -= $ga->{cdnabest} if($ga->{cdnabest});
            my $gbl= $gb->{aalen}; $gbl -= $gb->{cdnabest} if($gb->{cdnabest});
            my $bad=0; 
            if($gal < $gbl) { $bad=1; } elsif($gal > $gbl) { $bad=2; } 
            else { $bad=($ga->{cov} < $gb->{cov})?1 :($ga->{cov} > $gb->{cov})?2 :3; } 
            $keep= ($bad==1) ? 0 : 1;  $ndrop++ unless($keep);
          }
        }
        
      } # no else { }
    }
    print O $inrow if($keep);
  }
  close(O); close(F);
  if($err or $ndrop==0) { unlink($trfixgff); } # warn?
  else { 
    $trgff= $trfixgff; 
    # LOG result : 
    warn "#NOTE: trchimerpick ndrop=$ndrop of $nchi mRNA in $trgff\n";
    }
  return($trgff);
}

 
sub trfixup
{
  my($trgff, $trfa, $strandc, )= @_;
  my $cmd;
  my $name= $trfa || "trfixin"; $name =~ s/.gz//; $name=~ s/\.(tr|fasta|fa)$//;  

  unless($SKIP2BEST) {  

  if($SPLITJOINS) {
  my($trsplitfa)= trsplitjoins( $trfa, ); # test
  if($trsplitfa ne $trfa) { #? trmap() even if split fails?
    $trfa= $trsplitfa;   
    ($trgff)= trmap( $trsplitfa, ); #... 
    } 
  }
  
  if($CHIMERPICK) {
    ($trgff)= trchimerpick($trgff);
  }
  
  if($strandc) {
  ## FIXME: bug here? dropping some exons of kept mRNA : err=Missing-mrna-exon
  ## .. or below in strandc filtering.. all are chimera and '.' strand, prior were cdnabest w/ strand.
  ## .. maybe due to ';cdsindel=nnn\n' last mRNA annot, now has scaffold catenated == exon line?
  
    ## FIXME here: allow for unstranded 1-exon genes, use genefindcds prot to set strand
    my $err=0; my $keep=1; my $nok=0; 
    my $strf= strandt($strandc);
    my $trgffor= gffsuf($trgff,$strf,"gz");  
    if($trgff =~ m/\.gz/) { open(F,"gunzip -c $trgff |") or $err++; } else { open(F,"$trgff") or $err++; }
    open(O,">$trgffor") or $err++;
    while(<F>) {
      my $inrow= $_;
      if(/^\w/ and /\tmRNA/) { 
        my @v=split"\t"; my $or=$v[6]; 
        if($or eq $strandc) { $keep=1; }
        elsif($or eq "." and m/nexon=1;/) { 
          my %at=(); 
          while(m/(aalen|cov|pid|cdsindel)=(\d+)/g) { $at{$1}=$2; }
          my $pc= $at{cov} * $at{pid} / 100;
          $keep= ($pc>97 and $at{aalen}>=40 and $at{cdsindel} < 1)? 1: 0; # use global params here?
        } else { $keep=0; }
        $nok++ if($keep);
      }
      print O $inrow if($keep);
    }
    close(O); close(F);
    $err++ unless($nok);
    unless($err) { $trgff= $trgffor; }
  }
  
  my $trgffan= gffsuf($trgff,"an","gmap"); # $trgff; $trgffan =~ s/.gff/.an.gff/;
  my $dna= "$genomed/$dgenome.fa";
  unless( -f $dna ) { $dna="$genomed/$dgenome.fasta"; }
  unless( -f $dna ) { die "#FAIL: MISSING genome dna $dna"; }
  $cmd= "$evigene/scripts/genefindcds.pl -ratiocdna 1.25 -genes $trgff -cdna $trfa -dna $dna";
  
  # ** FIXME: introns use full path or use ../$introns for chdir outputdir
  if($introns) {
    if(-f $introns) { $cmd.= " -nofixintron -intron $introns "; }
    elsif( -f "../$introns") {  $cmd.= " -nofixintron -intron ../$introns "; }
    else { warn "#ERROR: MISSING introns $introns\n"; }
  }
  $cmd.= " > $trgffan";
  my $err= runcmd($cmd); # runCmdOrDie($cmd);

  $trgff= $trgffan;
  
  # LOG result : 
  my $ntrmap= `grep -c mRNA $trgff`; chomp($ntrmap);
  warn "#NOTE: trmap n=$ntrmap mRNA in $trgff\n";

  } # unless($SKIP2BEST) 
  
  unless($NOBEST) {
  # sub trfilterjunk .................
  # my $trgffbest= gffsuf($trgff,"best$vers","an"); # best1 >> best$vers   drop .fwd.an.gff
  my $trgffbest= ($trgff =~ /stdin|^-/) ? "${name}best$vers.gff" : gffsuf($trgff,"best$vers","an");
  
  my $cmd;
  ## fixme: inqual=$inok/$inerr/$inzip ; need to score -inerr > +inok
  # score: inqual:10  from genefindcds. should score aalen=size,pCDS,complete/partial
  # score: cdsindel:-9 from genefindcds << BUGgers; cdsindel neg ; need some scorewt fix: abs-cdsindel:-5 ??
  # score: chimera:-1 reduce? use path=n/m as proxy
  # score: cov=, pid= from gmap; score both?
  
  ##my $tscore= "inqual:10,cdsindel:-5,pid:2,path:-1,CDS:2,UTR:2";
  ##new# my $bestscores= "inqual:5,cov:8,pid:5,CDS:2,UTR:1,cdsindel:-2";
  my $tscore= $bestscores || "inqual:3,cov:5,pid:2,CDS:1,UTR:1";
  my $opadd="-dropscore $dropscores" if($dropscores);
  $cmd="$evigene/scripts/overbestgene1.perl -alttr -pct 10 -strand -score '$tscore' $opadd -in $trgff > $trgffbest";
  my $err= runcmd($cmd);

  my $ntrmap= `grep -c mRNA $trgffbest`; chomp($ntrmap);
  warn "#NOTE: trbest n=$ntrmap mRNA in $trgffbest\n";
  return($trgffbest);
  } else {
  return($trgff);
  }
  
}
