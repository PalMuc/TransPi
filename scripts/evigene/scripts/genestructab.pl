#!/usr/bin/env perl
# introntab.pl

use strict;
use warnings;
use Getopt::Long ;
# use File::Basename;

my $addex=1;
my $debug=0;
my $format=""; #"gtf"; #?
my $fixcds=0; # only for dromel_fb ??
my $predictor="";

my (%exons, %allexons, %start, %stop, %utr, %rsize, %badcall);
%exons= %allexons= %start= %stop= %utr= %rsize= %badcall= ();

my $orgfile= $ENV{orgfile} || "organism";

my $optok= &GetOptions (
  "format=s"=>\$format, # gtf or gff.3
  "organism=s"=>\$orgfile,
  "fixcds!"=>\$fixcds,
  "predictor=s"=>\$predictor,
  "debug!"=>\$debug,
  );

die "usage: introntab [-fixcds] [-predictor=augustus|..] --format gff|gtf < input.gff > output.tab\n"
unless($optok and $format);

my $isgtf= ($format =~ /gtf/i)? 1: 0;

## read: fill %exons{trid}, %start{trid}, %stop{trid}, %utr{trid}

if($isgtf) { readGtf(); }
else { readGff(); }

putTable();

#........................................................


sub readGtf {
  # my($gtfh)= @_;  
  # what is perlish for input all @ARGV + STDIN handle ? 
  warn "readGtf\n" if $debug;
  my ($n);
  while(<>) {
  next unless(/^\w/ and /\t/);

  my($r,$s,$t,$b,$e,$v,$o,$p,$a)=split"\t"; 
  ($b,$e)= ($e,$b) if ($b>$e);
  $p=0 if($p eq ".");
  my($gid)=m/gene_id "([^"]+)"/; 
  my($tid)=m/transcript_id "([^"]+)"/;  # use this instead of gid
  # fixme got dmel gene_id == tr_id; drop -Rx; 
  # ditto? worm
  ## ^^^^ FIXME: transcript_id Gene1,Gene2,Gene3,...
  
  if($t=~/start_codon/ and $a =~ m/Parent=Gene:(\w+)/) { $gid=$1; }
  elsif($t=~/start_codon/ and $a =~ m/Parent=(\w+)/) { $gid=$1; }
  
  if($gid =~ /,/) { ($gid)= split(/,/, $gid); }
  if($gid =~ /\-R\w/) { $gid =~ s/\-R.//; } # fly alt-tr ids
  
  ## Dang Fly GFF/GTF has multiple Parent IDs per CDS/exon/UTR ...
  my @tid= ($tid =~ /,/) ? split(/,/,$tid) : ($tid);
  
  foreach my $td (@tid) {
  if($t =~/CDS/) { 
    $badcall{$td}++ if(abs($e-$b)<3);
    push( @{$exons{$td}},[$b,$e,$o,$p])  unless(abs($e-$b)<3); 
      # unless($b == $e); # some odd predictor bug? do we keep these? skip any with $e-$b < 3
    } 
  elsif($t =~/start_codon/){     $n++;
    $badcall{$td}++ if($start{$td});
    $start{$td}=[$b,$e,$o,$p,$gid,$r] unless($start{$td}); 
    #   aberrant cases w/ 2 start_codon, due to 0 or 1-bp CDS ; leave as is?
    } 
  elsif($t =~/stop_codon/){ $stop{$td}=[$b,$e,$o,$p,$gid,$r]; } 
  elsif($t =~/UTR/){ push( @{$utr{$td}},[$b,$e,$o,$p,$gid,$r]) ;}
  }  
  
  unless($rsize{$r} and $rsize{$r}>$e) { $rsize{$r}= $e; }
  }
  warn "readGtf: ng=$n\n" if $debug;
  
}

sub readGff {
  # my($gtfh)= @_;  
  # what is perlish for input all @ARGV + STDIN handle ? 

  warn "readGff\n" if $debug;
  my %didg=(); my %didtr=(); my %newtid=();
  my($gid,$tid,$n,$trnum,$hasgene)=("","",0,0,0,0,0);
  while(<>) {

  unless(/^\w/ and /\t/) {
    $predictor="augustus" if(/^# This output was generated with AUGUSTUS/i && !$predictor);  
    next;
  }
  
  my($r,$s,$t,$b,$e,$v,$o,$p,$a)=split"\t"; 
  ($b,$e)= ($e,$b) if ($b>$e);
  $p=0 if($p eq ".");

  if($t =~ /^gene$/) {
    if($a =~ m/ID=([^;\s]+)/) { $gid=$1; $gid =~ s/^\w+://; } # drop gbrowse/wormbase prefix
    if($gid =~ /,/) { ($gid)= split(/,/, $gid); }
    if($gid =~ /\-R\w/) { $gid =~ s/\-R.//; } # fly alt-tr ids
    $hasgene=1;
    
  } elsif ($t =~ /^mRNA$/) {
    $gid=0 unless($hasgene); # should reset gid to 0 unless found t=gene
    if($a =~ m/ID=([^;\s]+)/) { $tid=$1; } #? $tid =~ s/^\w+://;
    if($a =~ m/Parent=([^;\s]+)/) { $gid=$1; $gid =~ s/^\w+://; }
    # ^^^^ FIXME: Parent=Gene1,Gene2,Gene3,...
    #? $tid= $gid unless($tid);
    # die/warn unless $tid here?
    
    $predictor= $s if(!$predictor && $s =~ /^(Gnomon|AUGUSTUS|SNAP|GLEAN|FGenesh)/i);
    
    if($predictor =~ /^aug/i) { $gid=$tid; $gid =~ s/\.?t\d+$//; }
    elsif($predictor =~ m/^Gnomon/i && $a !~ m/Parent=/ && $a =~ m/gene=(\w+)/) { $gid=$1; }
    
    unless($gid) {
      # else ?? leave empty
      }
    
# Fixme for PASA update gff that throws in alt-tr, new mRNA/gene w/o changing mRNA/CDS IDs
    $trnum= ++$didg{$gid}; # should use exon overlaps to detect alt-tr
    if($trnum>1 && $didtr{$gid}{$tid}) {
      my $newtid="$tid.$trnum"; # need to fix CDS,exon also
      $newtid{$tid}= $newtid;
      $tid= $newtid;
    }
    $didtr{$gid}{$tid}++;
     
    $n++;  
    if($gid =~ /,/) { ($gid)= split(/,/, $gid); }
    if($gid =~ /\-R\w/) { $gid =~ s/\-R.//; } # fly alt-tr ids
    
    # note this is mRNA start,stop not CDS/protein start,stop .. fixme
    ## also swap start,stop for - $o
    if($o eq '-') {
    $stop{$tid} =  [$b,$b+2,$o,$p,$gid,$r]; #  unless($start{$tid}); 
    $start{$tid}=  [$e-2,$e,$o,$p,$gid,$r];
    } else {
    $start{$tid}= [$b,$b+2,$o,$p,$gid,$r]; #  unless($start{$tid}); 
    $stop{$tid}=  [$e-2,$e,$o,$p,$gid,$r];
    }
    
  } elsif ($t =~/CDS/) {
    if($a =~ m/Parent=([^;\s]+)/) { $tid=$1;  } # $tid=~s/^\w+://;
    $tid= $newtid{$tid} || $tid;
    ## Dang Fly GFF/GTF has multiple Parent IDs per CDS/exon/UTR ...
    my @tid= ($tid =~ /,/) ? split(/,/,$tid) : ($tid);
    foreach my $td (@tid) {
    $badcall{$td}++ if(abs($e-$b)<3);
    push( @{$exons{$td}},[$b,$e,$o,$p]) unless(abs($e-$b)<3); 
    }

  } elsif ($t =~/exon/) {
    if($a =~ m/Parent=([^;\s]+)/) { $tid=$1;  } # $tid=~s/^\w+://;
    $tid= $newtid{$tid} || $tid;
    my @tid= ($tid =~ /,/) ? split(/,/,$tid) : ($tid);
    foreach my $td (@tid) {
    push( @{$allexons{$td}},[$b,$e,$o,$p]) unless(abs($e-$b)<3); 
    }
  }
  
  unless($rsize{$r} and $rsize{$r}>$e) { $rsize{$r}= $e; }
  }
  warn "readGff: ng=$n\n" if $debug;
  
  # need to calc utr from allexons & cds ...
  #  elsif($t =~/UTR/){ push( @{$utr{$tid}},[$b,$e,$o,$p,$gid,$r]) ;}
  if(%allexons && %exons) {
    foreach my $tid (sort keys %exons) {
      my @cds= @ { $exons{$tid} };
      ref( $allexons{$tid}) or next;
      my @exons= @ { $allexons{$tid} };
      my($mstart,$mstart1,$o,$p,$gid,$r)= @ { $start{$tid} };
      
      my ($cmin,$cmax)= (-1,0);
      foreach my $ex (@cds) { foreach my $i (0..1) { 
        my $v= $ex->[$i]; 
        $cmax=$v if($v>$cmax); 
        $cmin=$v if($cmin<0 || $v<$cmin); 
        } }

      foreach my $ex (@exons) {
        my($b,$e,$o,$p)= @$ex;
        if($b >= $cmin and $e <= $cmax) { # skip, cds
          next;
        } elsif($b < $cmin and $e > $cmin) { # split utr
          $e= $cmin-1;
        } elsif($b < $cmax and $e > $cmax) { # split utr
          $b= $cmax+1;
        } else { # outside cds
        }
        push( @{$utr{$tid}}, [$b,$e,$o,$p,$gid,$r]);
      }
      
      if($o eq '-') {
      $stop{$tid} =  [$cmin,$cmin+2,$o,$p,$gid,$r]; #  unless($start{$tid}); 
      $start{$tid}=  [$cmax-2,$cmax,$o,$p,$gid,$r];
      } else {
      $start{$tid}= [$cmin,$cmin+2,$o,$p,$gid,$r]; #  unless($start{$tid}); 
      $stop{$tid}=  [$cmax-2,$cmax,$o,$p,$gid,$r];
      }
       
    }
  }
  
}


# for fly, worm, have true gene id in start_codon notes: worm: Parent=Gene:WBGene
# fly: Parent=CGxxx
# add UTR stats (for worm, fly have xUTR feature, some in daph; not mouse)
# fixme: set min,max intron size;
# fixme: add "mean" table option replacing size array w/ mean (offsets, sizes, phases)

# x add intergenic value; need sorted start/stop
# add 'gene density'  measure of # genes >> total cds size / genome_size


sub putTable {
  print "# Intron table for $orgfile\n";
  my $gnosize=0; map{$gnosize+= $rsize{$_}}sort keys %rsize; 
  print "# Genome_size: $gnosize\n";
  print join("\t",
    qw(TranscrID GeneID Strand N_Intron Gene_size CDS_size UTR_size Interg_size
	     Intron_offsets  Intron_sizes Intron_phases),
       ($addex ? "Exon_sizes" : ""),
       "Trnum"),"\n";

  my %didg;
  my @startids = sort{
    ($start{$a}->[5] cmp $start{$b}->[5] 
    or $start{$a}->[0] <=> $start{$b}->[0] 
    or $start{$b}->[1] <=> $start{$a}->[1])
    } keys %start;
  my @stopids = sort{
    ($stop{$a}->[5] cmp $stop{$b}->[5] 
    or $stop{$a}->[0] <=> $stop{$b}->[0] 
    or $stop{$b}->[1] <=> $stop{$a}->[1])
    } keys %stop;
  my $laststop=0; my $lastref=0;

  my $nin= @startids;
  my $nbad= scalar keys %badcall;
  warn "putTable: nin=$nin, nbad=$nbad\n" if $debug;
  my $nout=0;
  
  foreach my $igene (0..$#startids) { #sort keys %start
    my $id= $startids[$igene];
    my $flag="";
    # always have start# $start{$id} or $flag.="nostart;";
    $stop{$id} or next; # $flag.="nostop;"; # partial gene models from Gnomon, mostly; should have pre-filtered
    next if $badcall{$id};
    
    my ($lmin,$lmax)= (-1,0);
    my($bb,$be,$geneo,$pb,$gid,$ref)= @{$start{$id}};
    my($eb,$ee,$go,$pe)= @{$stop{$id}}; # NOTE GTF stop_codon is part of last CDS !
    my @ex = ref($exons{$id}) ? @{$exons{$id}} : (); 
    $gid ||= $id; # no gene id?
    
    if($eb < $bb) { $lmin= $eb; $lmax= $be; }
    else { $lmin= $bb; $lmax= $ee; }
    
    ## FIXME for genes_drosmel_fb42: must remove UTR from exons
    if($fixcds) {
      my @cds=();
      foreach my $ex (@ex) {  
        ##my($b,$e,$o,$p)= @$ex;
        my ($vb,$ve)= @$ex;
        if( $ve < $lmin || $vb > $lmax ) { next; }# UTR exons
        elsif($vb < $lmin) { $ex->[0]= $vb= $lmin; } # cut out UTR part
        elsif($ve > $lmax) { $ex->[1]= $ve= $lmax; }
        push(@cds, $ex);
        }
      @ex= @cds;
    }
    next unless(@ex);  # $nbad++; Mouse has some genes/start_codon w/ no exons; skip
   
    foreach my $ex (@ex) { foreach my $i (0..1) { 
      my $v=$ex->[$i]; 
      $lmax=$v if($v>$lmax); 
      $lmin=$v if($lmin<0 || $v<$lmin); 
      } }
    
    my $trnum= ++$didg{$gid}; # should use exon overlaps to detect alt-tr
    
    my $igenic= 0;
    if($igene>0) { 
    ## need to filter out alt tr
      my $lid= $startids[$igene-1]; 
      
      if($start{$lid} and $stop{$lid}) {
      my($lbb,$lbe,$lgn,$lpp,$lgg,$llref)= @{$start{$lid}};
      my($leb,$lee,$lgo,$lpe,$lgid,$lastref)= @{$stop{$lid}};
      my ($llmax)= $lee;
      $llmax= $leb if($llmax < $leb);
      $llmax= $lbb if($llmax < $lbb);
      $llmax= $lbe if($llmax < $lbe);
      my $thisb= $lmin;
      ## if $llb < $ee and $lle > $bb
      $igenic= ($ref eq $lastref && $llmax>0 && $llmax < $thisb) ? $thisb - $llmax : 0;
      }
      }
      
    ## geneo may be bad sometimes? check each exon strand?
    my $no=0; foreach (@ex) { $geneo= $$_[2] and $no++ if($$_[2] =~ /[+-]/ && $geneo ne $$_[2]); }
    $no == 0 or $flag.="nostrand;";
    my $rev= ($geneo eq "-" || ($geneo =~ /\d/ and $geneo < 0));
    
    # always have start# unless($start{$id}) { ($bb,$be)= ($ex[0]->[0], $ex[0]->[1]); } #? or find min
    my $geneb = $rev ? $be : $bb;
    my $genee = $rev ? $eb : $ee;
    
    my $nintron= @ex - 1;
    my $genespan= 1 + abs($genee - $geneb); # or rev
    my $cdslen= 0;
    my (@inoffset, @insize, @inphase, @exsize);
    
    @ex= sort{ # by distance from gene start
      my $be= $rev ? 1 : 0; # end, begin
      my $ad= abs($a->[$be] - $geneb); 
      my $bd= abs($b->[$be] - $geneb);  
      $ad <=> $bd } @ex;
    
    my($lb,$le,$lo,$lp);
    my($b,$e,$o,$p)= @{$ex[0]};
    $cdslen += 1 + abs($e - $b) - $p; push(@exsize, 1+abs($e-$b));
    foreach my $i (1..$nintron) {
      ($lb,$le,$lo,$lp)= ($b,$e,$o,$p);
      ($b,$e,$o,$p)= @{$ex[$i]};
      if($isgtf and $i == $nintron) { if($rev) { $b -= 3; } else { $e += 3; } } # GTF: add stop_codon length always 3
      $cdslen += 1 + abs($e-$b) - $p; # dont forget phase
      push(@exsize, 1+abs($e-$b));
      my $inoffset= $rev ? abs($geneb - $e) : abs($b - $geneb); # why the sign errs?
      my $insize  = $rev ? abs($lb - $e) : abs($b - $le); # make this abs() ?
      push(@inoffset, $inoffset);
      push(@insize, $insize);
      push(@inphase, $p);
      }

    my $utrlen=0;
    my @utr = (ref $utr{$id}) ? @{$utr{$id}} : (); 
    foreach my $i (0..$#utr) {
      ($b,$e,$o,$p)= @{$utr[$i]};
      $utrlen += 1 + abs($e-$b);
    }
    
    ## double check that $cdslen + sum($insize) = $genespan;  
    ## most are off by 1 errs here, or 3, or similar.
    my $insum=0; foreach (@insize) { $insum+=$_; }
    if($insum + $cdslen != $genespan) {
       # my $gs= $cdslen + $insum; #$flag .="badspan;"; 
       # warn "genespan $genespan != $gs ($cdslen + $insum, $geneo)\n";
       }
    $genespan= $insum + $cdslen; # force it to match
    
    foreach (@inoffset) { $_ = sprintf("%.3f", $_ / $genespan); } #?
    
    $flag = "err:$flag" if $flag; # leave out but for debug
    print join("\t",$id,$gid, $geneo,$nintron,$genespan,$cdslen,$utrlen, $igenic,
               join(",",@inoffset), join(",",@insize), join(",",@inphase),
               ($addex ? join(",",@exsize) : ""),
               $trnum,
	       ), "\n";
	  $nout++;     
    }
  
  $nbad= $nin - $nout;
  print "# Totals: n_in: $nin, n_bad: $nbad, n_out: $nout\n";
    
}

__END__


=item newer stats

  see daphwork/dpxgenestruc/daphnia-genostats2.txt

=item convert to two tabs...

  cat $genes.gff | perl genestructab.pl --format gff > $genes.introns.tab

  # this is mess, revise genestructab.pl to do this..
foreach intab (*.introns.tab)
  set otab=`echo $intab | sed -e 's/.tab/.mean/'`
  cat $intab | perl -ne\
  'chomp; @v=split"\t"; if(/#/ or not $v[3]=~/\d/){print "$_\n";} \
  else { $maxin=9999; $maxex=9999; $is=$v[9]; $es=$v[11]; \
  @is=split",",$is; $n=0; $im=0; map{ if($_ < $maxin){ $n++; $im+=$_;}}@is; $n||=1; $im=sprintf "%.1f",$im/$n; \
  @es=split",",$es; $n=0; $em=0; map{ if($_ < $maxex){ $n++; $em+=$_;}}@es; $n||=1; $em=sprintf "%.1f",$em/$n; \
  $v[9]=$im; $v[11]=$em; print join("\t",@v),"\n";}'\
  > $otab

  set otab=`echo $intab | sed -e 's/.tab/.parts/'`
  cat $intab | perl -ne\
  'chomp; @v=split"\t"; if(/#/ or not $v[3]=~/\d/){print "$_\n";} \
  else { $gi=$v[1]; $is=$v[9]; $es=$v[11]; $v[8]=$v[10]=0; \
  if($gi ne $lastgi) { @is=split",",$is; @es=split",",$es; \
  for $i (0..$#es) { $v[9]=$is[$i]||0; $v[11]=$es[$i]||0;  print join("\t",@v),"\n"; } \
  } $lastgi=$gi; }'\
  > $otab
end


=item genestruc R stats
  
  options(digits=4)

  genomes <- c( "celegans_wb167","daphnia_ncbi1", "daphnia_jgi11", "drosmel_fb42","mouse_mgi3")
  
  gfile <- c(
  culex_cpip12 = "culex_cpip12",culex_pasavalid = "culex_pasa",
  aphid_augmap5c = "aphid_augmap5c",
  ); #...
  
  
  gfname <- function (gname) {
   return( ifelse( is.na(gfile[gname]), gname,  gfile[gname]) );
  } 
  
  gsize <- c(
  aphid_acyr1= 348299353,
  aphid_augmap4a= 364191497,
  ); #....

==> poptr.introns.tab <==
# Intron table for organism
# Genome_size: 399762734

==> sorbi.introns.tab <==
# Intron table for organism
# Genome_size: 690458708

==> soybn.introns.tab <==
# Intron table for organism
# Genome_size: 955433471

==> vitvi.introns.tab <==
# Intron table for organism
# Genome_size: 485480514
> gsize
    cacao3gi     aratht10    cacaocir1 cacao3eg2cir 
   335474211    119667750    318866582    318866582 

  gsize <- c(
  cacao3gi= 335474211,
  aratht10= 119667750,
  cacaocir1= 318866582,
  cacao3eg2cir= 318866582,
  poptr= 399762734,
  vitvi= 485480514,
  soybn= 955433471,
  sorbi= 690458708
  ); #....


  #.... stats two .......

  # gas <- c()
  gas <- gas.all
  # gas <- gas.all[gas.all$species != "ixodes_scaf10",]
  
  for (gname in genomes) {
    #gfn <- paste("genes_",gfname(gname),".introns.mean.gz",sep="");
    #ga <- read.delim(gzfile(gfn),comment.char="#"); 
    gfn <- paste(gfname(gname),".introns.mean",sep="");
    ga <- read.delim(gfn,comment.char="#"); 
    cat("gname=",gname,", dim=",dim(ga),"\n")
    ga$Interg_size[ ga$Interg_size == 0] <- NA
    ga$UTR_size[ ga$UTR_size < 20 | ga$UTR_size > 9999 ] <- NA
    ga$Intron_sizes[ ga$Intron_sizes < 30 | ga$Intron_sizes > 29999 ] <- NA
    ga$Exon_sizes[ ga$Exon_sizes < 20 | ga$Exon_sizes > 9999 ] <- NA
    gas <- rbind(gas, cbind( species=rep(gname,dim(ga)[1]) ,ga) )
  }
  gas.all<- gas
  
  # gas <- c()
  gas <- gas.inex
  # gas <- gas.inex[gas.inex$species != "ixodes_scaf10",]
  
  for (gname in genomes) {
    #gfn <- paste("genes_",gfname(gname),".introns.parts.gz",sep="");
    #ga <- read.delim(gzfile(gfn),comment.char="#"); 
    gfn <- paste(gfname(gname),".introns.parts",sep="");
    ga <- read.delim(gfn,comment.char="#"); 
    cat("gname=",gname,", dim=",dim(ga),"\n")
    ga <- ga[,c("TranscrID","GeneID","Strand","N_Intron","Gene_size","Intron_sizes","Exon_sizes","Trnum")]
    ga$Intron_sizes[ ga$Intron_sizes < 30 | ga$Intron_sizes > 29999 ] <- NA
    ga$Exon_sizes[ ga$Exon_sizes < 20 | ga$Exon_sizes > 9999 ] <- NA
    gas <- rbind(gas, cbind( species=rep(gname,dim(ga)[1]) ,ga) )
  }
  gas.inex<- gas
  
  stats <-  c("N_Intron","Gene_size","CDS_size","Exon_sizes","Intron_sizes","UTR_size","Interg_size")
  oneexon<- 0   # 0,1,2
  
  gall<- c()
  for (gname in genomes) {
    gout<- gstructstat(gname, oneexon)
    gall <- rbind(gall, gout)
    print(format(gout,scientific=F))
    cat( "#------------------------------\n")
  }

  gstructab2( gall)
  # print(format( gstructab2( gall),scientific=F))

  > gstructab2( gall)[,c(1:4,6)]
                N_Intron          Gene_size           CDS_size        Exon_sizes         UTR_size
  aratht10  2 4.06 +/- 0.031  1538 1849 +/- 9.1  1032 1204 +/- 5.5  135 241 +/- 0.87  316 352 +/- 1.3
  cacao3gi  2 3.73 +/- 0.028  1993 3091 +/- 62    935 1148 +/- 5.5  138 246 +/- 0.93  562 624 +/- 1.8
  cacaocir1 2 4.02 +/- 0.029  1938 2679 +/- 16    945 1170 +/- 5.7  137 237 +/- 0.88  403 663 +/- 6.3
    
  # struct    spp1   spp2   spp3 ...
  #  cds       mean +- sem (median)
  #  exon  ...
  
  gstructab2 <- function( gall, gstats=stats)
  {
    gnames<- levels( factor(gall[,"Genome"]))
    # subst <- c("Median","Mean","SEM")
    gmat <- matrix( 0, ncol=length(gstats), nrow=length(gnames)) # need 3/stat: mean, sem, median
    gsout<- data.frame(gmat, row.names=gnames)
    colnames(gsout) <- gstats
    for (gn in gnames) {
      rg <- (gall[,"Genome"] == gn)
      ga<- gall[rg,]
      for (st in gstats) {
      rs <- (ga[,"Statistic"] == st)
      n  <- ga[rs,"N"]
      md<- ga[rs,"Median"]
      mn<- ga[rs,"Mean"]
      sem<- ga[rs,"SD"] / sqrt(n)
      gsout[gn,st] <- paste( format(md), format(mn), "+/-", format(sem), sep=" ")
      }
    }
    return(gsout)
  }

  gstructab1 <- function( gall, gstats=stats)
  {
    gnames<- levels( factor(gall[,"Genome"]))
    subst <- c("Median","Mean","SEM")
    gstatx<- paste( rep(gstats,each=3), rep(subst,length(gstats)), sep=".")
    gmat <- matrix( 0, nrow=length(gstatx), ncol=length(gnames)) # need 3/stat: mean, sem, median
    gsout<- data.frame(gmat, row.names=gstatx)
    colnames(gsout) <- gnames
    for (gn in gnames) {
      rg <- (gall[,"Genome"] == gn)
      ga<- gall[rg,]
      for (st in gstats) {
      rs <- (ga[,"Statistic"] == st)
      n  <- ga[rs,"N"]
      gsout[paste(st,"Median",sep="."), gn] <- ga[rs,"Median"]
      gsout[paste(st,"Mean",sep="."), gn] <- ga[rs,"Mean"]
      gsout[paste(st,"SEM",sep="."), gn]  <- ga[rs,"SD"] / sqrt(n)
      }
    }
    return(gsout)
  }
    
  # sem <- function( sd, n) { sem <- sd/sqrt(n); return(sem) }
  
  gstructstat <- function(gname, oneexon=0) 
  {
    gbases <- gsize[ gname]
    if(oneexon==1) { cat("One exon only\n");} else if(oneexon>1) { cat("Two+ exons only\n"); }
    out<- c()
    for (stat in stats) {
      if( stat %in% c("Intron_sizes","Exon_sizes")) {
        gas <- gas.inex ;  
      } else {
        gas <- gas.all ;  
      }
      statna <- stat
      sel <- gas$species==gname & !is.na(gas[,stat])  
      if( !(stat %in% c("Intron_sizes","Exon_sizes","UTR_size")))  sel <- sel & (gas$Trnum == 1) 
      if(oneexon==1) { sel<- sel & gas$N_Intron == 0;  } 
      else if(oneexon==2) { sel<- sel & gas$N_Intron > 0;}
    
      ga1 <- gas[ sel, stat]  
      gn  <- length( ga1)
      gm  <- ifelse(gn<2,NA,mean(ga1,na.rm=T))
      gmd <- ifelse(gn<2,NA,median(ga1,na.rm=T))
      gsd <- ifelse(gn<2,NA,sd(ga1,na.rm=T))
      
      selb <- sel & (gas$Trnum == 1)
      gt  <- sum( gas[selb,stat],na.rm=T)
      gtb <- gt / gbases
      out<- rbind(out, c( gname, statna, gn, gmd, gm, gsd, gtb, gt))
    }
    out<- data.frame(out,stringsAsFactors=F)
    for (i in 3:8) out[,i] <- as.numeric(out[,i])
    colnames(out) <- c( "Genome", "Statistic", "N", "Median", "Mean", "SD", "PerGenome", "Sum")
    return(out)  
  }


=item intronmodes R stats

see Burges Pnas 2001 paper: fit 2 lognormal distr to find intron size modes?
Lim LP, Burge CB. (2001). A computational analysis of sequence features involved in 
recognition of short introns. Proc. Natl. Acad. Sci. USA 98, 11193-11198. doi:10.1073/pnas.201407298 

  # fixme: xcut=200 option should be measured.
  
  library(MASS)
  set.seed(123)
  
  gins <- c()
  for (gname in levels(gas.inex$species)) {
    gins <- rbind(gins, c(species=gname, gintronlogn(gname)) )
  }
  gins <- data.frame(gins)
  rownames(gins) <- gins[,1]; gins<- gins[,-1]
  gins
  #                    short.mode short.est short.sd long.mode long.est long.sd pLong      n
  #   celegans_wb167           51        61     1.48       502      553       2 0.327 100639
  #   mouse_mgi3              107       112     1.37      1588     1661    3.08 0.856 159339
  #   daphnia_pasaj            70        75     1.27       464      597    2.54 0.096  49391

  gintronlogn <- function(gname, xcut=200, gas=gas.inex) 
  {
    sel<- gas$species == gname & !is.na( gas$Intron_sizes) 
    introns <- gas[sel,"Intron_sizes"]
    gintron1logn(introns, xcut)
  }
  
  gintron1logn <- function(introns, xcut=200) 
  {
    nt <- length(introns)
    
    md.lo <- median(introns[introns < xcut],na.rm=T)
    fit <- fitdistr(introns[introns < xcut], "lognormal"); 
    est.lo <- exp(fit$estimate); 
    
    nhi <- length(introns[introns >= xcut])
    md.hi <- median(introns[introns >= xcut],na.rm=T)
    fit <- fitdistr(introns[introns >= xcut], "lognormal"); 
    est.hi <- exp(fit$estimate); 
    
    return( c(short.mode=round( md.lo), short.est=round(est.lo[[1]]), short.sd=round(est.lo[[2]],2),
          long.mode=round( md.hi), long.est=round(est.hi[[1]]), long.sd=round(est.hi[[2]],2),
          pLong = round(nhi / nt,3), n=nt))
  }


=cut




######################################
# 
# =item OLDstats
# 
# # this is for GTF format (sigh...) genes data; should add GFF handling.
# # reads gene feature data, writes table with Intron values per gene
# #   GeneID Strand N_Intron CDS_size Intron_offsets  Intron_sizes Intron_phases
# 
# # dgg.090708: ** this is a bad bug: map{ next; } causes exit of top (print) loop, so we lost any genes w/ > maxin
# # caused mouse exon count to drop by 1+; others also ??  $maxin=9999; should set maxin > 20000 ? real drosmel introns > 100k **
# #bad#  @is=split",",$is; $n=0; $im=0; map{ next if($_ >$maxin); $n++; $im+=$_;}@is; $n||=1; $im=sprintf "%.1f",$im/$n; \
# #bad#  @es=split",",$es; $n=0; $em=0; map{ next if($_ >$maxex); $n++; $em+=$_;}@es; $n||=1; $em=sprintf "%.1f",$em/$n; \
# 
# foreach gtf (genes*.gtf.gz)
#   set orgfile=`echo $gtf | sed -e 's/.gtf.gz//'`
#   setenv orgfile $orgfile
#   gzcat $gtf | perl introntab.pl > $orgfile.introns.tab
# 
#   cat $orgfile.introns.tab | perl -ne\
#   'chomp; @v=split"\t"; if(/#/ or not $v[3]=~/\d/){print "$_\n";} \
#   else { $maxin=9999; $maxex=9999; $is=$v[9]; $es=$v[11]; \
#   @is=split",",$is; $n=0; $im=0; map{ if($_ <$maxin){ $n++; $im+=$_;}}@is; $n||=1; $im=sprintf "%.1f",$im/$n; \
#   @es=split",",$es; $n=0; $em=0; map{ if($_ <$maxex){ $n++; $em+=$_;}}@es; $n||=1; $em=sprintf "%.1f",$em/$n; \
#   $v[9]=$im; $v[11]=$em; print join("\t",@v),"\n";}'\
#   > $orgfile.introns.mean
# end
# 
# 
# R stats for these data:
# 
# genomes <- c("genes_apis_ncbi4","genes_celegans_wb167","genes_daphnia_ncbi1",
#   "genes_drosgri_caf1","genes_drosmel_fb42","genes_mouse_mgi3")
# 
# gsize <- c(
# genes_apis_ncbi4= 217085197,
# genes_celegans_wb167= 100255891,
# genes_daphnia_ncbi1= 169048782,
# genes_drosgri_caf1= 161589839,
# genes_drosmel_fb42= 117614329,
# genes_mouse_mgi3= 2592998464
# );
# 
# for (gname in genomes) {
# ga <- read.delim(paste(gname,".introns.mean",sep=""),comment.char="#"); 
# 
# # set Interg_size == 0 <- NA ; set UTR_size == 0 <- NA
# ga$Interg_size[ ga$Interg_size == 0] <- NA
# ga$UTR_size[ ga$UTR_size == 0] <- NA
# 
# gbases <- gsize[ gname]
# cat( "Genome", "Statistic", "N", "Median", "Mean", "SD", "PerGenome", "Sum", "\n")
# for (stat in c("N_Intron","Gene_size","CDS_size","UTR_size","Interg_size","Intron_sizes","Exon_sizes")) {
#   gn  <- sum( ! is.na(ga[,stat] ))
#   gt  <- sum(ga[,stat],na.rm=T)
#   gtb <- gt / gbases
#   gm  <- mean(ga[,stat],na.rm=T)
#   gmd <- median(ga[,stat],na.rm=T)
#   gsd <- sd(ga[,stat],na.rm=T)
#   # fixme: alt-tr mess up gtb, gt counts: need to remove same exons ...
#   cat( gname, stat, gn, gmd, gm, gsd, 
#   gtb, gt, "\n", sep="\t")
# }
# 
# cat( "#------------------------------\n")
# }
# 
# 
# =cut

