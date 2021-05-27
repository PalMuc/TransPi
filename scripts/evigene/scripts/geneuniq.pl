#!/usr/bin/env perl
# geneuniq.pl 

=item about

  geneuniq : make unique gene.gff subset from input many alt model
  assumes input.gff(s) are gene-record ordered (mRNA/exon,cds all together, can have gene row)
  output is first input of unique model, order inputs by prefered output
  
=cut

use strict;
use Getopt::Long;

my $debug=1;
my $mrnatypes='mRNA';
my $exontypes='exon|CDS';
my ( $ingff, $outgff )= ("","");
my $aabest=0;

my $optok= GetOptions(
  "exontypes=s", \$exontypes,
  "aabest!", \$aabest,
  "debug!", \$debug, 
  );

my ($ngene, $nmrna, $nexon, $nother, $ndup, $nuniq, $predid, 
    $modrows, $modscore, $cdsw, $exonw) = (0) x 10;
my (%models, @modid, %modstore);

$exontypes= join('|', split(/[,;\s|]+/,$exontypes));

sub putif {
  my($modid, $modrows, $modscore, $flag)= @_;
  
  my $lastscore= $models{$modid};
  if(!defined $lastscore or $lastscore < $modscore) { 
    $models{$modid}= $modscore;
    $modstore{$modid}= $modrows;
  }
  if(defined $lastscore) { $ndup++; } else { $nuniq++; }
  
  if($flag =~ /LAST/) {
    foreach my $id (sort keys %modstore) { print $modstore{$id}; }
  }
  
#   if($models{$modid}++) { $ndup++; }
#   else { print $modrows; $nuniq++; }  
}


while(<>) {

  unless(/^\w/) {
    print unless($nmrna); # print header up to 1st data?
    next; # dont print all comments..
  }

  my $row= $_;
  my @col= split"\t";  
  my $locid= join"_",@col[0,3,4,6,2]; # location,location,location,type
  
#  # need to collect all of gene model before know if it is uniq
#   if( $col[2] eq "gene" ) { # assume none? or save to keep w/ 1st mrna?
#   }    

  if( $col[2] =~ m/^($mrnatypes)$/) { 
    
    if($modscore == 1 and $cdsw > 0 and $exonw > 0) { 
      $modscore= int(0.5 + 100*$cdsw/$exonw); 
    }
    putif( join(",", sort @modid), $modrows, $modscore, "") if(@modid);

    @modid= (); #($locid); ## $modid= $locid; ## MAYBE drop mRNA locid and use only exon locs for id?    
    $modrows= $row;
    ($predid)= m/\bID=([^;\s]+)/;    
    $nmrna++; 
    $exonw= $cdsw= 0;
    $modscore=1;
    if($aabest) { 
      if( my($al,$ap)= m/aalen=(\d+).(\d+)/ ) { $modscore= $ap; }
      # percent cds, bigger = best in the case of false utrs
      # better = 60% cds/exon, where >90% is poor and >40% is poor
      ## better score below from sum(CDS)/sum(EXON)
      }
        
    
  } elsif( $col[2] =~ m/^($exontypes)$/ ) {
    my ($pid)= m/\bParent=([^;\s]+)/;

    if($pid ne $predid) {
      die "# out of order gene rec: $pid.$col[2] not in mRNA:$predid\n"; # and do what? any of these?
      # next; #??? better to die
    }
    
    if($col[2] eq "exon") { $exonw += 1 - $col[4] - $col[3]; }
    elsif($col[2] eq "CDS") { $cdsw += 1 - $col[4] - $col[3]; }

    ## urk, need to sort locid by type/location here so the modelid has std ordering
    push(@modid, $locid); ## $modid= $locid; 
    $modrows .= $row;
    $nexon++;
    
  } else { 
    $nother++;
    if($col[2] eq "exon") { $exonw += 1 - $col[4] - $col[3]; }
    elsif($col[2] eq "CDS") { $cdsw += 1 - $col[4] - $col[3]; }
    
    if( @col > 6 ) { # and $col[0] eq $lchr
      #NO# $modid   .= $locid; 
      $modrows .= $row; # preserve data but not location id    
    }
  }

}  # while in

putif( join(",", sort @modid), $modrows, $modscore, "LAST") if(@modid);
warn "#geneuniq: n.mrna=$nmrna, n.kept=$nuniq, n.dupl=$ndup \n";



__END__

gzcat \
genejc/tiling.adult.female.augmap.an5.gff.gz \
genejc/tiling.adult.male.augmap.an5.gff.gz \
genejc/tiling.embryo_10hr.female.augmap.an5.gff.gz \
genejc/tiling.embryo_10hr.male.augmap.an5.gff.gz \
genejc/tiling.embryo_18hr.female.augmap.an5.gff.gz \
genejc/tiling.embryo_18hr.male.augmap.an5.gff.gz \
genejc/tiling.larvae51hr.female.augmap.an5.gff.gz \
genejc/tiling.larvae51hr.male.augmap.an5.gff.gz \
genejc/tiling.tesova.ovaries.augmap.an5.gff.gz \
genejc/tiling.tesova.testes.augmap.an5.gff.gz \
genejc/tiling.wing.female.augmap.an5.gff.gz \
genejc/tiling.wing.male.augmap.an5.gff.gz \
genejc/tiling.yellow.female.augmap.an5.gff.gz \
genejc/tiling.yellow.male.augmap.an5.gff.gz \
genejc/rnaseq.xu004.augmap.an5.gff.gz \
genejc/rnaseq.s567.augmap.an5.gff.gz \
genejc/rnaseq.s3t3.augmap.an5.gff.gz \
genejc/rnaseq.016r1.augmap.an5.gff.gz \
genejc/rnaseq.009r1.augmap.an5.gff.gz \
genejc/pasa.augmap.an5.gff.gz \
genejc/homology.both.augmap.an5.gff.gz \
genejc/golden.augmap.an5.gff.gz \
| $evigene/scripts/geneuniq.pl \
> genes/alljc_uniq.augmap.an5.gff

#geneuniq: n.mrna=482360, n.kept=188314, n.dupl=294046 


gzcat \
genejc/tiling.adult.female.augmap.an5.gff.gz \
genejc/tiling.adult.male.augmap.an5.gff.gz \
genejc/tiling.embryo_10hr.female.augmap.an5.gff.gz \
genejc/tiling.embryo_10hr.male.augmap.an5.gff.gz \
genejc/tiling.embryo_18hr.female.augmap.an5.gff.gz \
genejc/tiling.embryo_18hr.male.augmap.an5.gff.gz \
genejc/tiling.larvae51hr.female.augmap.an5.gff.gz \
genejc/tiling.larvae51hr.male.augmap.an5.gff.gz \
genejc/tiling.tesova.ovaries.augmap.an5.gff.gz \
genejc/tiling.tesova.testes.augmap.an5.gff.gz \
genejc/tiling.wing.female.augmap.an5.gff.gz \
genejc/tiling.wing.male.augmap.an5.gff.gz \
genejc/tiling.yellow.female.augmap.an5.gff.gz \
genejc/tiling.yellow.male.augmap.an5.gff.gz \
genejc/pasa.augmap.an5.gff.gz \
genejc/homology.both.augmap.an5.gff.gz \
genejc/golden.augmap.an5.gff.gz \
genejc/rnaseq.xu004.augmap.an5.gff.gz \
genejc/rnaseq.s567.augmap.an5.gff.gz \
genejc/rnaseq.s3t3.augmap.an5.gff.gz \
genejc/rnaseq.016r1.augmap.an5.gff.gz \
genejc/rnaseq.009r1.augmap.an5.gff.gz \
| $evigene/scripts/geneuniq.pl -exon CDS \
> genes/alljc_cdsuniq.augmap.an5.gff

#CDS-uniq:
#geneuniq: n.mrna=482360, n.kept=58599, n.dupl=423761 

#.......
melon2.% ls {genes,genejc}/*.{gff,gff.gz}
## uniqify subset of these, not mixes

genes/nvit_epit3-augmap.gff.gz
genes/nvit_epi4-augmap.gff.gz
#? genes/bestgenes.mix7a.gff.gz@   
# genes/ogs12.gff.gz@

genes/ogs12.an2.gff.gz    
genes/nasv1_augustus0802.gff.gz
genes/nasv_pred_gnomon.gff.gz

# genes/ogs12.an1.gff.gz
# genejc/bestgenes_of7.an7a.gff.gz
# genejc/bestnvit2_mix6.an7a.gff.gz
# genejc/bestgenes_of5.an7rna.gff.gz
# genejc/bestgenes_of4.an7tifem.gff.gz
# genejc/bestgenes_of6.an7tiem.gff.gz
# genejc/bestgenes_of4.an7timal.gff.gz
# genejc/bestgenes_of7.an6best.gff.gz
# genejc/ogs12.gff.gz@

genejc/tiling.adult.female.augmap.an5.gff.gz
genejc/tiling.adult.male.augmap.an5.gff.gz
genejc/tiling.embryo_10hr.female.augmap.an5.gff.gz
genejc/tiling.embryo_10hr.male.augmap.an5.gff.gz
genejc/tiling.embryo_18hr.female.augmap.an5.gff.gz
genejc/tiling.embryo_18hr.male.augmap.an5.gff.gz
genejc/tiling.larvae51hr.female.augmap.an5.gff.gz
genejc/tiling.larvae51hr.male.augmap.an5.gff.gz
genejc/tiling.tesova.ovaries.augmap.an5.gff.gz
genejc/tiling.tesova.testes.augmap.an5.gff.gz
genejc/tiling.wing.female.augmap.an5.gff.gz
genejc/tiling.wing.male.augmap.an5.gff.gz
genejc/tiling.yellow.female.augmap.an5.gff.gz
genejc/tiling.yellow.male.augmap.an5.gff.gz
genejc/rnaseq.xu004.augmap.an5.gff.gz
genejc/rnaseq.s567.augmap.an5.gff.gz
genejc/rnaseq.s3t3.augmap.an5.gff.gz
genejc/rnaseq.016r1.augmap.an5.gff.gz
genejc/rnaseq.009r1.augmap.an5.gff.gz
genejc/pasa.augmap.an5.gff.gz
genejc/homology.both.augmap.an5.gff.gz
genejc/golden.augmap.an5.gff.gz

# genes/bestgenes_of7.an5.gff
# genes/all.tiling.adult.male.augmap.an5.gff
# genes/all.tiling.adult.female.augmap.an5.gff
# genes/all.rnaseq.xu004.augmap.an5.gff
# genes/all.rnaseq.s567.augmap.an5.gff
# genes/all.rnaseq.s3t3.augmap.an5.gff
# genes/all.rnaseq.016r1.augmap.an5.gff
# genes/all.rnaseq.009r1.augmap.an5.gff
# genes/nasvit_mix4a.gff
# genes/nvit2_mix7asm1.gff.gz
# genejc/bestnvit2_mix6.gff.gz
# genes/nvit2_mix6asm1.an3.gff.gz
# genes/rnaseq.s3t3.augustus.an3.best.gff.gz
# genes/pasa.augustus.an3.best.gff.gz
# genes/all.rnaseq.xu004.augmap.gff.gz
# genes/all.rnaseq.s567.augmap.gff.gz
# genes/all.rnaseq.s3t3.augmap.gff.gz
# genes/all.rnaseq.016r1.augmap.gff.gz
# genes/all.rnaseq.009r1.augmap.gff.gz
# genes/all.tiling.adult.male.augmap.gff.gz
# genes/all.tiling.adult.female.augmap.gff.gz
