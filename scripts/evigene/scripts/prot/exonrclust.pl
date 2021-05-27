#!/usr/bin/env perl
# exonrclust.pl

use strict;
use Getopt::Long;

my %args=();
my $ok= &GetOptions( \%args,
  "exonerate=s" , "input=s" , "output=s" , "tmpdir=s",
  "i|islice=i" , "n|ncpu=i" , "stranded!", "debug!" );

## added tmpdir 201208
## FIXME: input table now can have 2+ aa per row/genome span, comma-sep ?
## FIXME2: add ID-strand info at request to force exonerate-stranded

## add these opts for big prots?
#   --percent 90  << bad bad bad
#   --bestn 3 is ok
#   --subopt 0  # turn off?, default on
#   --refineboundary 100 # ?dont need huge for cds refinement ??
# bad was:  --minintron 20 --maxintron 10000; default mini = 30, maxi = 200k good enough

## prot optionally try: --model protein2genome:bestfit --exhaustive 1
##   .. guaranteed to be much,much slower .. but 'the entire protein is included in the alignment'
## ** (process:18237): WARNING **: Exhaustively generating suboptimal alignments will be VERY SLOW
## =>> PBS: job killed: walltime 50142 exceeded limit 50100
## **> drop subopt for bestfit/exhaustive

## est try: --splice3 splice3.mat --splice5 splice5.mat to improve --forcegtag 1
## est try: --model cdna2genome instead, models CDS also

## ryo add: et,ei,es = Equivalenced {total,id,similarity,mismatches}
##     add: pi,ps = percent eqiv id,similarity over total

my $xshow=" --showtargetgff --showvulgar 0 --showalignment 0 ";
my $xryo= " --ryo '#qi %qi length=%ql alnlen=%qal\\n#ti %ti length=%tl alnlen=%tal\\n#pi %pi,%ps total=%et ident=%ei sim=%es\\n'";

## prot: try --forcegtag for proper splice sites?
## prot: try  --softmasktarget
#  --refineboundary 50 ?? or use default 32 ?
# my $exonrcmd_prot="exonerate --model protein2genome --bestn 3 --refine region " . $xshow . $xryo;
my $exonrcmd_prot="exonerate --model protein2genome --subopt 0 --refine region " 
  . "--forcegtag 1 --softmasktarget 1 " . $xshow . $xryo;

## ** (process:18237): WARNING **: Exhaustively generating suboptimal alignments will be VERY SLOW
## **> drop subopt for bestfit/exhaustive
my $exonrcmd_bestprot="exonerate --model protein2genome:bestfit --exhaustive 1 --subopt 0 " 
  . " --forcegtag 1 --softmasktarget 1 ". $xshow . $xryo;

## for tilex, try  --joinrangeint 24 --joinrangeext 24 (default=12) --refine region (32b)
## slower, but maybe what is needed.  does refine do same or not?
my $exonrcmd_est="exonerate --model est2genome --forcegtag 1 --refine region "
 . " --joinrangeint 24 --joinrangeext 24 ". $xshow . $xryo;

## not useful:
my $exonrcmd_cdna="exonerate --model cdna2genome --forcegtag 1 " . $xshow . $xryo;

my $exonerate = $exonrcmd_prot;
if($args{exonerate}) {
  if($args{exonerate} =~ /^(prot)/i) { $exonerate= $exonrcmd_prot; }
  elsif($args{exonerate} =~ /^(bestprot)/i) { $exonerate= $exonrcmd_bestprot; }
  elsif($args{exonerate} =~ /^(EST)/i) { $exonerate= $exonrcmd_est; }
  elsif($args{exonerate} =~ /^(cDNA)/i) { $exonerate= $exonrcmd_cdna; }
  else { $exonerate= $args{exonerate}; } 
}

my $input = $args{input} || shift(@ARGV);
my $output= $args{output} || "";
my $tmpdir= $args{tmpdir} || "";
my $ncpu= $args{n} || 1;
my $icpu= $args{i} || 0;
my $debug= $args{debug} || 0;
my $stranded= $args{stranded};
my $q=0;
my $pid=$$;  # process id for temp file names


sub puto { 
  my($a,$b,$c)=@_; my @b= split",",$b ; my @c= split",",$c; 
  open(F,">$a"); foreach $b (@b) { $c=shift @c; print F ">$b\n$c\n"; } close(F); 
} 

warn "# in=$input; out=$output; ncpu=$ncpu; icpu=$icpu;\n# command=$exonerate \n" if($debug);

if($tmpdir) {
  if(-d $tmpdir) { $tmpdir.="/" unless($tmpdir =~ m,/$,); }
  else { warn "ERR: missing tmpdir=$tmpdir \n"; $tmpdir=""; }
}

die "bad opt or input " unless($ok and $input);
my $inh=undef;
if($input =~ /stdin|-/) { $inh=*STDIN; } 
else { open(IN, $input) or die "cant read $input"; $inh=*IN; } # add aatab.gz input

## need output write locking???? among procs 
system("touch $output") if $output;
while(<$inh>) {
  next unless(/^\w/ and ($q++ % $ncpu) == $icpu); 
  chomp; my($exp,$aid,$aa,$nid,$na)=split"\t"; 
  next unless($aid =~ /\w/ and $aa =~ /\w/ and $nid =~ /\w/ and $na =~ /\w/); #err?
  
  my @tf= map{ "${tmpdir}xr$pid.$q.".$_ } qw(aa na out annot);  
  warn "# q=$q: $aid x $nid\n" if $debug;

  my $cmd= $exonerate;
  
  puto($tf[0],$aid,$aa); 
  puto($tf[1],$nid,$na);
  
  ## $nid == location-id, chr:b-e:strand
  ## per option, use $aid\tstrand table for exonerate
  if($stranded) {
    # if using cdna2genome, check for gene model CDS start,end ; add to .annot
    my ($or)= $nid =~ m/([+-])$/;
    if($or) { 
      my ($cdsb,$cdse)= $aid =~ m/:(\d+)-(\d+)/;
      my $cbe= ($cdse > 0) ? "\t$cdsb\t$cdse" : "";
      open(F,">".$tf[3]); print F "$aid\t$or$cbe\n"; close(F);   
      $cmd .= " --annotation $tf[3]"; 
    }
  }  
  
  $cmd .= " --query $tf[0] --target $tf[1] > $tf[2]"; 
  my $ok= system($cmd); 
  # write out $tf[2] here
  if($output) { system("cat $tf[2] >> $output"); }
  else { `cat $tf[2]`; } # to this stdout
  
  foreach (@tf) { unlink($_); } # ok if missing $_ ? 
} 
close($inh);

__END__

=item usage cluster call exonerate


  #! /bin/bash
  ### env exontab=nasvit1-protexonregion.tab qsub -q normal exonrclust.sh
  #PBS -N exonr1
  #PBS -A xxxxx
  #PBS -l nodes=1:ppn=32,walltime=13:55:00
  #PBS -o exonr1.$$.out
  #PBS -e exonr1.$$.err
  #PBS -V
  
  ncpu=32
  
  exbin=$HOME/bio/exonerate/bin
  export PATH=${exbin}:${PATH}
  workd=$HOME/scratch/chrs/nasv1
  dgenome=nasvit1asm
  
  cd $workd/prot/
  
  if ! test -f $exontab ; then  echo "missing input $exontab"; exit;  fi
  nam=`echo $exontab | sed 's/\..*//' `
  echo "# start exonrclust: $nam `date`"
  
  i=0 
  while [ $i -lt $ncpu ]; do {
    out=$nam.exonr$i.gff ;  touch $out
  
    echo "# perl ./exonrclust.pl -in=$exontab -out=$out -i=$i -n=$ncpu "
    # Fails # perl exonrclust.pl -debug -in=$exontab -out=$out -i=$i -n=$ncpu 
    ( cat $exontab | perl exonrclust.pl -debug -in=stdin -out=$out -i=$i -n=$ncpu ) &
  
    i=$(( $i + 1 ))
  }
  done
  
  wait
  echo "# end exonrclust: $nam `date`"

=cut
