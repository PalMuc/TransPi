#!/usr/bin/env perl
# ncbi_mapview2gff.pl

=item usage

  gunzip -c seq_gene.md.gz | grep Primary_Assembly | \
  ncbi_mapview2gff.pl -sou hsa37 > ncbiseq_gene.gff

  see also http://promotion.molgen.mpg.de/gbrowse/contrib/import_ncbi_mv_hs.pl

=cut

use strict;
use warnings;
use Getopt::Long;

my $gffsource= "ncbimap";
my $ADDEXON=defined($ENV{exon})?$ENV{exon}:1; # only CDS,UTR in data?
my $ADDEVD=$ENV{evd}||0;

my @ncbifields =  # match your data
  qw(   
    taxid 
    chr start stop strand 
    contig cstart cstop cstrand
    fname fid ftype 
    group_label
    besthit_acc  
    evd_code  
    );
    
#    besthit_acc  besthit_gi
#    evd_code  has5p has3p

## new columns?
#9606    1       11874   14409   +       NT_077402.2     1874    4409    +       LOC100287102    GeneID:100287102        GENE    Primary_Assembly        -       mRNA;identical;N
#9606    1       11874   14409   +       NT_077402.2     1874    4409    +       XM_002342010.1  GeneID:100287102        RNA     Primary_Assembly        XM_002342010.1  mRNA;identical;N
#9606    1       11874   12189   +       NT_077402.2     1874    2189    +       -       GeneID:100287102        UTR     Primary_Assembly        -       -


## update here, now have RNA => mRNA, GENE => gene, was GENE => mRNA
my %NCBI2SO_TYPE = (
  RNA => 'mRNA', # or ncRNA
  GENE => 'gene', ##was 'mRNA',
  CDS => 'CDS',
  UTR => 'UTR', # five_prime_UTR or three_prime_UTR ... fixme
  PSEUDO => 'pseudogene', # this is flag or mRNA ? duplicates mRNA, evd=pseudo ?
);
my $lastxe=0; my @lexon;
our(%trs,$tb,$te);

my $optok = GetOptions(
  "source=s"=>\$gffsource,
  );

die << "USAGE" unless($optok);
  ncbi_mapview2gff [-source $gffsource] < hmm_models.md > hmm_models.gff
USAGE

ncbi_mapview2gff(*STDIN);

sub ncbi_mapview2gff {
  my($inh)= @_;
  print "##gff-version 3\n";
  my @save=(); my $glid; my %ids;
  while(<$inh>){
    next unless(/^\w/); # is it always ^\d taxid ?
    chomp;
    
    my %row = ();
    @row{ @ncbifields } = split "\t"; # [@field_positions]
    
    my $gfftype= $NCBI2SO_TYPE{ $row{ftype} } || $row{ftype};
    my $ischild= ($gfftype =~ /CDS|exon|UTR/)?1:0;
    my $nam= $row{fname}; $nam="" if($nam eq "-"); # was $id, not working now, NM_nnn and NP_nnn same mRNA, diff ids
    $gfftype="ncRNA" if($gfftype =~ /mRNA/ and $nam =~ /XR_|NR_/);

    my $lid= $row{fid};  ##BELOW: $glid=$lid; # GeneID:nnnnnn
    if($nam and $lid) { $ids{$lid}= $nam; } elsif($lid and !$nam){ $nam=$ids{$lid}||$lid; }
    
    my $trid= $lid; # GeneID, not good for alt-tr
    my $altid= $row{"besthit_acc"}; $altid="" if($altid eq "-");
    $trid=$altid if($altid);
    
    my $attr= ($ischild ? "Parent" : "ID") . "=$trid";
    $attr .= ";Name=$nam" if($nam =~ /\w/ and $nam ne $trid);
    
    my $dbx="";
    $dbx  .= "$lid," if($lid ne $trid);
    $dbx  .= "$altid," if($altid and $altid ne $trid);
    # $dbx  .= join ",", map{ $row{$_} if ($row{$_} && $row{$_} ne "-"); } qw(besthit_acc); ## besthit_gi
    $dbx =~ s/,$//;
    $attr .= ";Dbxref=". $dbx if ($dbx =~ /\w/);
    
    $row{evd_code} =~ s/\;+/,/g;
    $row{evd_code} =~ s/\-//g;
    $row{evd_code} =~ s/ab initio/ab_initio/g;
    $attr .= ";evd=$row{evd_code}" if ($row{evd_code} and $ADDEVD);
    
    my @gff= ($row{chr}, $gffsource, $gfftype, $row{start}, $row{stop}, '.', $row{strand}, '.', $attr);

    if($gfftype =~ /gene/) {  $trs{$lid}= [ $row{start}, $row{stop} ]; } #was /RNA/
    elsif($gfftype =~ /RNA/) {  $trs{$trid}= [ $row{start}, $row{stop} ]; }  
    # NOTE: GeneID not trID, need alt tr from name ..
    
    sub same { my($ga,$gb)=@_; for my $i (0..7) { return 0 if($$ga[$i] ne $$gb[$i]); } return 1; }
 
    if( $ischild && same( \@gff, \@save )) {
      $save[8] =~ s/Parent=([^;]+)/Parent=$1,$trid/;
    } else {
      putgff($glid,\@save,\@gff) if @save;
      @save= @gff; $glid=$trid;
    }  
    # need to combine CDS, UTR, .. of same location, diff parents, as one row
  }
  putgff($glid,\@save,[]) if (@save); ##print join("\t",@save), "\n" if @save;

}

sub putgff {
   my($glid,$save,$gff)=@_; my @save=@$save; my @gff=@$gff;

   ## our($tb,$te);
   if(@save) {
      #see above if($save[2]=~/RNA/) { my @tbe= @save[3,4]; $trs{$glid}= \@tbe; }
      if($ADDEXON and $save[2]=~/UTR|CDS/) {
	      my $haslexon=@lexon;
        my $csplit=($gff[2] =~ /UTR/ and $save[2] =~ /CDS/) ? 2: ($save[2] =~ /UTR/ and $gff[2] =~ /CDS/)? 1 :0;
        my @exon=($haslexon) ? @lexon :  @save; 
        $exon[2]="exon";  my($xb,$xe)= @exon[3,4];
        if($xe == $lastxe and not $haslexon) { 
          $lastxe= -1;
        } 
        else {
          $lastxe= -1;
 	        if($xe+1 == $gff[3] and $csplit) {
          $xe= $lastxe= $exon[4]= $gff[4]; # add skip this gff on next go
          }

	      ## FIXME: many glids possible, see Parent=x,y,z
        ($xb,$xe)= @exon[3,4];
	 	    $exon[-1]=~s/;xt=[\w\.,]*//; 
        my @trids=($glid);
	      my($trids)= $exon[-1] =~ m/Parent=([^;\s]+)/;
	      if($trids =~ /,/){  @trids= split( /,/, $trids); }
	      my %xtyp;
	      foreach my $tid (@trids) {
          my $xtyp=""; ## need mRNA span to know if exon first,last,single
          my($tb,$te)= ($trs{$tid}) ? @{$trs{$tid}} : (0,0);
          if($xb == $tb and $xe == $te) { $xtyp = "single"; }
          elsif($xb == $tb) { $xtyp = ($exon[6] eq "+")?"first":"last"; }
          elsif($xe == $te) { $xtyp = ($exon[6] eq "-")?"first":"last"; }	
          elsif($xb > $tb and $xe < $te) { $xtyp = "inner"; }
          $xtyp{$xtyp}++ if($xtyp);
          }
	      my $xtyp= join",", sort keys %xtyp;
	      ## $xtyp.="$tb.$xb,$te.$xe"; #DEBUG bug
	      $exon[-1]=~s/$/;xt=$xtyp/; # got some nulls, why?

	      if($haslexon or !$csplit) { print join("\t",@exon), "\n" ; @lexon=(); }
        else { @lexon= @exon; }
 	      # print join("\t",@exon), "\n" ; # need to save if single exon gene, UTR/CDS/UTR to get 2nd utr
        }
      } elsif(@lexon) {
	      print join("\t",@lexon), "\n"; @lexon=();
      }
      print join("\t",@save), "\n";
  }
}


__END__

microbe% zmore dmel5_ncbi_hmm_models.md.gz
------> dmel5_ncbi_hmm_models.md.gz <------
#tax_id chromosome      chr_start       chr_stop        chr_orient      contig  ctg_start       ctg_st
op      ctg_orient      feature_name    feature_type    group_label     best_hit_acc    best_hit_gi   
        evidence_code   is_5p_complete  is_3p_complete
7227    2L      6774    9489    +       NT_033779.4     6774    9489    +       hmm4134014      GENE  
        reference       AAZ86793.1      GI:73853446     psCDS;; N       N
