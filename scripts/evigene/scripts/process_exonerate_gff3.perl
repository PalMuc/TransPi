#!/usr/bin/perl -w

=head1 NAME 

process_exonerate_gff3 -t EST file1.e2g.EXONERATE file2.e2g.EXONERATE > out.gff3

=head1 USAGE

Command line arguments
 -gf|gff_version  gff version string (1,2,2.5,3) default 3
 -t|type EST or Protein typically, the string in the ID= field default EST 

 -idprefix XXX
 -source   exonr.xyz
 -dropfeat  intron,splice_acceptor,splice_donor
 
=head1 DESCRIPTION

Turns EXONERATE gff output into GFF for Gbrowse use.

You need to have run exonerate with at least the following options

 --showtargetgff yes
 --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n"

Although I often use this for fungal protein mapping

 --showvulgar yes --softmaskquery yes --softmasktarget yes --minintron 20 --maxintron 3000 --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" --showalignment no

=head1 AUTHOR

Jason Stajich, jason-at-bioperl-dot-org

=item dgg notes

   revised the -ryo prefix to '#ti, #qi';
   want to drop some of these feats: splice_donor intron splice_acceptor
    .. grep out or here?
    
   change ID prefix from Protein/EST to choice.. e.g modmm
   change exonerate source field to choice

  gzcat exonerate-dmel.gff.gz | \
  $aug/scripts/process_exonerate_gff3.perl \
   -t Protein -a -id FB  -source exonerate.modDM \
   -drop 'exon,intron,splice,HSP' -rename 'gene=mRNA' \
   > exonerate-dmel.gff3

=item qucker oneliner

  see sub quickexonr2gff
  avoids all those bioperl objects that slow it down
  custom set for protein2genome

=item strand bug; now fixed see STRAND_FLIP

  values like -11129364 
  < in
  scaffold_4770   exonerate:protein2genome   gene    13128542        13129365        724     -    . gene_id 1 ; sequence CG14394-PA ; gene_orientation +
  
  > out
  scaffold_4770   exonerate.modDM mRNA    -11129364       -11128541       724     -       .       ID=FB:CG14394-PA

 check exonerate's odd use of -strand; does it mean all of target/chr was
 flipped and 2nd scan for exons done? all of these follow all +strand matches.
 
 STRAND_FLIP is  wrong for exonerate 1.4.0 has correct gff loc for - by default
 src/c4/alignment.c
 ArgumentSet_add_option(as, '\0', "forwardcoordinates", NULL,
        "Report all coordinates on the forward strand", "TRUE",
   
 * also checked output matches to known tblastn align: correct for no STRAND_FLIP

=cut

use strict;
use Getopt::Long;

## only for slow bp_ version > required
# use IO::String;
# use Bio::Tools::GFF;
# use Bio::Seq;
# use Env;
# use File::Spec;

my $state = 0; my $ct = 0;
my $buffer = '';
my (%dat,%lengths,%alnlens,%counter,$vulgar,$lastvulgar);


use constant STRAND_FLIP => 0;

my $type = 'Protein'; # EST or Protein
my $gff_ver = 3;
my $show_aln = 0; #? on by default or not

my ($source,$dropfeat,$renamefeat,$show_aln_hsp)=("") x 9;
my %renamefeat=();
my %didcomm;
my $idprefix= undef;
my ($part_ref,$part_start,$part_end)=(0) x 9; # dgg; genome partitions
my $skipaaloc=0;

my $fast = 1;
my $MINALIGN=$ENV{MINALIGN} || 20; 
my $SOURCE=  $ENV{SOURCE}; ## ||"exonerate.modDM"; # fixme
my $keepfeat=$ENV{keepfeat} || 'gene|cds';
my $refispart=$ENV{refpart} || 0;
my $MINDUPALIGN= 0.95 ; # suboptimal aligns of same query, drop below this
my $usegscore= 0;
my $usegeneor= undef;

my $optok= GetOptions(
  'fast!'      => \$fast,  # not used now
  'alignshow|showaln!' => \$show_aln, # or reverse to skipaln : aaloc=
  't|type=s'   => \$type,
  'source=s'   => \$SOURCE,
  'minalign=i' => \$MINALIGN,
  'dupalign=i' => \$MINDUPALIGN,
  'dupscore!' => \$usegscore, # fixme ; only for dupalign ?
  'orientation!' => \$usegeneor, 
  'keepfeat=s'  => \$keepfeat,
  # 'dropfeat=s'  => \$dropfeat, # not used now
  # 'renamefeat=s'   => \$renamefeat, # not used now
  'idprefix:s' => \$idprefix,  # not used now
  'refpart!'   => \$refispart,
  'skipaaloc!' => \$skipaaloc,
  'gff_version|gf|v|:s' => \$gff_ver, # not used now
  'h|help'   => sub { exec('perldoc', $0); exit }
  );

die "usage: process_exonerate_gff3 -help [...] < exonerate.gff > exonr_clean.gff3\n" 
unless($optok);

$idprefix= $type unless(defined $idprefix); # allow no prefix ""
$idprefix.=":" if($idprefix); # not used now: want?

$MINDUPALIGN= $MINDUPALIGN/100 if($MINDUPALIGN > 1);

## *** usegeneor IS BAD FOR PROT, default=1 only for EST
unless( defined $usegeneor ) { $usegeneor = ( $type =~ /prot/i ) ? 0 : 1; }

# fixme: which?
$keepfeat =~ s/mRNA/gene/;
$keepfeat =~ s/[,\s]+/\|/g if($keepfeat);
$dropfeat =~ s/[,\s]+/\|/g;  # not use now
if ($renamefeat) { %renamefeat= map{ split "=" } split /[,\s]/, $renamefeat;  } # not used now

quickexonr2gff(); # good

# if($fast) {  quickexonr2gff(); } 
# else {  bp_exonerate2gff(); }
	   


## this info is needed: align regions on protein:
# chr2L	exonerate:protein2genome	similarity	64892	65335	354	+	.	alignment_id 3 ;
#   Query dmoj_GLEAN_02256 ; Align 64892 225 24 ; Align 64916 234 420

# need add options to remove 'redundant' 2ndary matches to same location
#   -- same prot, diff prot ... need some quality cutoff to keep different ones?
#   -- need to sort by locations 1st; do after this step : process_exonr --reduce < step1.gff > step2.gff


#SCAFFOLD1:105958-120915 exonerate:protein2genome:local  similarity      1001    13961   2634    +       .       
# alignment_id 1 ; Query apis2gno_456253 ; Align 1001 294 165 ; Align 2050 349 102 ; Align 2152 384 111 ; Align 5690 421 108 ; 
# Align 5982 457 57 ; Align 6039 477 54 ; Align 6186 495 54 ; Align 6243 513 66 ; Align 6853 536 93 ; Align 9141 567 222 ;
# Align 9491 642 66 ; Align 9557 665 141 ; Align 11101 713 198 ; Align 11657 779 243 ; Align 13295 861 246 ; Align 13881 94381

# --- END OF GFF DUMP ---
#
#qi apis2gno_456253 length=1128 alnlen=676
#ti SCAFFOLD1:105958-120915 length=14958 alnlen=12961
#
# -- completed exonerate analysis
# Command line
#>>>>>>>>>>>>>>>>>>>

# SCAFFOLD1       exonr   mRNA    106958  119918  2634    +       .       
# ID=apis2gno_456253;
# >> drop aaloc= unless requested
# aaloc=880-1044:1,1045-1146:1050,1150-1260:1152,1261-1368:4690,1369-1425:4982,1429-1482:5039,1483-1536:5186,1537-1602:5243,1606-1698:5853,
# 1699-1920:8141,1924-1989:8491,1993-2133:8557,2137-2334:10101,2335-2577:10657,2581-2826:12295,2827-2907:12881;
# align=59; <<< this is PERCENT, alignlen = 676; FIXED
# aalen=1128
#<<<<<<<<<<<<<<<<<

# ** FIXME: drop suboptimal aligns at same locus; should tell exonerate to do that
#  skip() if($qid eq $lastqid and $pid < 0.9 * $lastpid);

# FIXME maybe: add scores to CDS from similarity; need to match locations: hstart,hlen : 
#   CDS score = sum(hlen) from similar hstart-hend overlaps (can be several/cds)

sub quickexonr2gff 
{
  
  my(@gvec,%gs,%didcom,%glocs);
  my($gn,$gid,$gor,$gscore,$glocid,$typ,$flip,$partlen,$partaln,
     $part_ref,$part_start,$part_end,$gene_start,$didpartcom,
     $ngene, $nlowalign, $nbador,  $lastpid, $lastgscore, $lastqid) = (0) x 30;

print"##gff-version 3\n"; 
while(<>){

  if(/^\w/ && /\t/) { 
    chomp; 
    my @v=split"\t"; 
    $typ= $v[2]; 
    my $keep=($typ=~/$keepfeat/i)?1:0;  # only feats to keep
    
    if($typ=~/^gene/) { 
      ($gid)= m/sequence (\S+)/; $gn=$gid; if( my $n= $gs{$gid}++ ){ $gn.=".".(1+$n); }
      
      #PROBLEM.NOW, but useful: gor # 
      ($gor)= m/gene_orientation (\S+)/;
      
      #not-used# $flip= ($v[6] eq "-")? 1:0; # exonerate -strand; no need to correct
      $v[8]="ID=$gn"; 
      $v[2]=~s/gene/mRNA/; 
      $gene_start= $v[3];
      $gscore    = $v[5];
      $ngene++;
      $glocid="";
      ## process refispart here at startgene ?? ; will be same to endgene
      }
      
    elsif($typ=~/^(cds|exon|intron)/i) { $v[8]="Parent=$gn";  $v[2]=~s/cds/CDS/;} 
 
    elsif($typ=~/^similarity/) {
      my $ploc="";
      while(m/Align\s+(\d+)\s+(\d+)\s+(\d+)/g){
        my ($hstart, $qstart, $hlen)=($1,$2,$3);
			  if( $type =~ /Protein/i ) { $qstart = $qstart * 3 - 2; } # dna or prot-based loc?
			  $hstart = 1 + $hstart + $part_start - $gene_start; # - genestart ?
			  my $qend= $qstart + $hlen - 1;
			  $ploc.= "$qstart-$qend:$hstart,";
        }
      $ploc=~s/,$//; 
      
      #?? instead add "Target=gid astart aend;" ??
      $gvec[0] =~ s/$/;aaloc=$ploc/ if(@gvec and $ploc and $show_aln); # put in gene record
      }
      
    $v[1]= $SOURCE if($SOURCE); 
    push @gvec, join("\t",@v) if $keep; 
    $glocid .= join("\t",@v[0,3,4,6,2])."\n" if $keep; # add type?
  }
    
  #qi apis2gno_456253 length=1128 alnlen=676
  elsif(/^#qi (\S+)\s+length=(\d+)\s+alnlen=(\d+)/) { 
    my($qid, $qln, $qal)= ($1,$2,$3);  $qln||=1;
    my $pal= int( 0.5 + 100*$qal/$qln); 
    if(@gvec) {  # was align=pct;aalen=qln
      unless( $gvec[0] =~ s+;aaloc=+;align=$pal\%,$qal/$qln;aaloc=+ ) {
      $gvec[0] =~ s+$+;align=$pal\%,$qal/$qln+;  
      }
    } 
    
    # see also gene gscore; align is not always useful: :bestfit always 100% prot align but diff gscore;
    # * now from bestfit have identicals, E2BZ94_9HYME and E2BZ94_9HYME.2 same CDS
    # * and 2ndary lower qual
# SCAFFOLD25  367232  368543  132 -  ID=UniRef50_Q9NAR7;align=100%,258/258;alignx=26%,51%,66,129/253
# SCAFFOLD25  367232  368406  294 -  ID=UniRef50_E2BHN0;align=100%,261/261;alignx=38%,61%,96,155/255 < best

    if($glocid and $glocs{$glocid}++) { $nlowalign++; @gvec=(); $gs{$qid}--;} # filter out dup identicals    
    elsif($pal<$MINALIGN) { $nlowalign++; @gvec=(); $gs{$qid}--;} # filter out crappy matches; need OPTION for this
    elsif($qid eq $lastqid and $usegscore and  $gscore < $MINDUPALIGN * $lastgscore) { $nlowalign++; @gvec=(); $gs{$qid}--;}
    elsif($qid eq $lastqid and $pal < $MINDUPALIGN * $lastpid) { $nlowalign++; @gvec=(); $gs{$qid}--;}
    else { $lastpid= $pal; $lastqid= $qid; $lastgscore= $gscore; }
    }
    
  #ti SCAFFOLD1:105958-120915 length=14958 alnlen=12961
  elsif(/^#ti (\S+)\s+length=(\d+)\s+alnlen=(\d+)/) { 
    my $tid; ($tid, $partlen, $partaln)= ($1,$2,$3);   #? not saved
    # $gvec[0] =~ s,$,;alignt=$partaln/$partlen, if(@gvec); # dont need?
    }

  ## add new ryo scores:
  #pi 49.59,66.12 total=121 ident=60 sim=80    # total == alnlen-1, using same tot as :bestfit will be lower pi,ps
  elsif(/^#pi ([\d.]+).([\d.]+)\s+total=(\d+)\s+ident=(\d+)\s+sim=(\d+)/) { 
    my($pi, $ps, $tn, $ti, $ts)= ($1,$2,$3,$4,$5); # want all these bits? only ti,ts/tn ?
    $pi=int(0.5+$pi); $ps=int(0.5+$ps); 
    # above 'qi pcti' has diff denom; this pi=ti/tn where tn=aligned query bases, not total query
    $gvec[0] =~ s+$+;alignx=$pi\%,$ps\%,$ti,$ts/$tn+ if(@gvec); # put in gene record
    }
    
  elsif(/^#part_location: (\w+):(\d+)-(\d+)/){  # dgg; genome partition header info
    unless($didpartcom) {
      ($part_ref,$part_start,$part_end)=($1,$2,$3);
      print $_; # dont do both this and below .. prefer this one?
      $part_start-- if($part_start>0);
      $didpartcom=1;
    }
    }
    
  elsif(/^Command line:/){  
    ($part_ref,$part_start,$part_end)=(0) x 3; $didpartcom=0;
    if(/\-target.(\S+)/) {
      my $targ=$1; 
      if($targ =~ m,(\w+)[:_](\d+)\-(\d+)/[^/]+$,) { # problems?
        ($part_ref,$part_start,$part_end)=($1,$2,$3);
        print "#part_location: $part_ref:$part_start-$part_end\n";
        $part_start-- if($part_start>0); 
        $didpartcom=1;
        }
      } 
    }
  elsif(/^##(source|date)/){ my $c=$1; s/#//; print unless($didcom{$c}++); }
  
  elsif(/START OF GFF|completed exonerate/) { 
    if(@gvec) {
      ##  flip needs to change for exonerate 1.4.0 has correct gff loc for - by default

      ## -- $gor flip is problem now, exonr 2.2, some cases - all?
      ## -- exonr is mixing up --annotation strand and splice strand (GT/AG = +; AC/CT = -)
      ## .. reports --annot strand even when splices are opposite ..
      ## also use $gor eq "." for no-strand
      ## *** usegeneor IS BAD FOR PROT, use only for EST
      
      if($usegeneor and $gor =~ /[+\.-]/) {
        my ($xref,$xor)= (split"\t",$gvec[0])[0,6];
        my ($refor,$bador)= (0,0);
        if($refispart and $xref =~ m/:([+-])$/) { 
          $refor=$1; $bador= ($gor ne "." and $refor ne $gor)?1:0; 
          if($bador) { @gvec=($gvec[0]); $gvec[0] =~ s/^/#bador#/; $nbador++; } #?
          }
        if($xor ne $gor and not $bador) {
          $gvec[0] =~ s+$+;reor=$xor/$gor+;
          for(my $i=0; $i<@gvec; $i++) {  
          my @v=split"\t",$gvec[$i];  $v[6]= $gor; $gvec[$i]= join("\t",@v);
          }
        }
      }
      
     if($refispart) { 
       #? do as default if $gvec[0] =~ /^\w+[:_]\d+\-\d+/
       #? need before endgene# ($part_ref,$part_start,$part_end)=(0) x 3; $didpartcom=0;
       map{ my @gr=split"\t";
         if( my($rr,$rb,$re)= $gr[0] =~ m/^(\w+)[:_](\d+)\-(\d+)/ ) {
          #? need before endgene# ($part_ref,$part_start,$part_end)=($rr,$rb,$re); $didpartcom=1;
          $rb-- if($rb>0); @gr[0,3,4]= ($rr,$rb+$gr[3],$rb+$gr[4]); $_= join("\t",@gr);  
          }
       } @gvec;
      #? need before endgene# print "#part_location: $part_ref:$part_start-$part_end\n";
      }

      print join("\n",@gvec),"\n"; 
    }
    @gvec=(); $gene_start=0;
    }
    
  }

  print '# ',join(", ", "ngene=$ngene", "skipped-lowalign($MINALIGN\%)=$nlowalign", "nbador=$nbador" ),"\n" ; 

}


### dont need, quickr is good n quick
###
# sub bp_exonerate2gff {
# 
# require IO::String;
# require Bio::Tools::GFF;
# require Bio::Seq;
# require Env;
# require File::Spec;
# 
# my $out = Bio::Tools::GFF->new(-gff_version => $gff_ver);
# 
# $out->_print("##gff-version $gff_ver\n"); # dang lack of control; want comments below
# $out->{'_first'}=0;
# 
# while(<>) {
#     if( $state == 0 && s/^vulgar:\s+// ) {
# 	chomp($vulgar = $lastvulgar = $_);
#     }
#     if( $state > 0 ) {	 
# 	if( $state == 2 ) {
# 	    if ( s/^vulgar:\s+//) {
# 		$lastvulgar = $vulgar;
# 		chomp($vulgar = $_);
# 	    }
# 	    my $in = Bio::Tools::GFF->new
# 		(-gff_version => 2,
# 		 -fh          => IO::String->new($buffer));
# 	    my ($gene_name,$gene)=("") x 9;
# 	    my $length=0;
# 	    while( my $f = $in->next_feature ) {
# 	        my $srctag  = $f->source_tag;
# 		$srctag =~ s/\:/_/g;
# 		$srctag= $SOURCE if $SOURCE;
# 		$f->source_tag($srctag);
# 		if( $f->primary_tag eq 'gene' ) {
# 		    $length = $lengths{$f->seq_id};		    
# 
# if(STRAND_FLIP) {		
# 		    if(  ! defined $length ) {
# 			die("unknown length for '",$f->seq_id,"'\n"); #why die?
# 		    }
# }
# 		    ($gene) = $f->get_tag_values('sequence');
# 		    if ( $f->has_tag('gene_orientation')) {
# 			if( $dat{$gene} ) {
# 			    $gene .= $counter{$gene};
# 			}
# 			($dat{$gene}) = $f->get_tag_values('gene_orientation');
# 			$dat{$gene} = 0 if $dat{$gene} eq '.';
# 			$f->remove_tag('gene_orientation'); 
# 		    }
# 		    for my $t ( qw(gene_id sequence) ) {
# 			$f->remove_tag($t) if $f->has_tag($t); 
# 		    }
# 		    $gene_name = $gene;
# 		    if( $counter{$gene_name}++ ) {
# 			$gene_name .= ".$counter{$gene_name}";
# 		    }
# 
# 		    $f->add_tag_value('ID',"$idprefix$gene_name");
# 		    $f->add_tag_value('vulgar', $lastvulgar) if $lastvulgar;
# 		    
# 		} elsif( $f->primary_tag eq 'similarity') {
# 		    next unless $show_aln;
# 		    # dgg ; keep stats on Gene line if no show
# 		    if( $f->has_tag('Align') ) { 
# 			my @hsps = $f->get_tag_values('Align');
# 			$f->primary_tag('match');
# 			$f->remove_tag('alignment_id');
# 			
# 			my ($min,$max);
# 			while ( @hsps ) {
# 			    my ($hstart, $qstart,$hlen) = splice(@hsps,0,3);
# 			    my $hlen_mod = $hlen;
# 			    if( $type =~ /Protein/i ) {
# 				$hlen_mod /= 3;
# 			    }
# 
#                            unless($dropfeat and "HSP" =~ m/$dropfeat/) { 
# 			    my $newf = Bio::SeqFeature::Generic->new
# 				( -seq_id => $f->seq_id,
# 				  -start => $hstart,
# 				  -end   => $hstart + $hlen,
# 				  -strand=> $dat{$gene},
# 				  -source_tag => $srctag,
# 				  -primary_tag => 'HSP',
# 				  -tag => {
# 				      'Target' => sprintf("$idprefix%s %d %d",
# 							  $gene_name,
# 							  $qstart,
# 							  $qstart+$hlen_mod),
# 				  }
# 				  );
# 			    $out->write_feature($newf);
# 			    }
# 			    $min = $qstart unless defined $min;
# 			    $max = $qstart + $hlen_mod;
# 
# 			}
# 			$f->add_tag_value('Target',
# 				    sprintf("$idprefix%s %d %d",
# 					    $gene_name, $min, $max,
# 					    ));
# 			my $genelen = $lengths{$gene};
# 			$f->add_tag_value('length', $genelen);
# 			$f->add_tag_value('coverage', int(100 * ($max-$min)/$genelen));  
# 			$f->remove_tag('Align');
# 		    }
# 		} elsif( $f->primary_tag eq 'splice5' ) {
# 		    $f->primary_tag('splice_donor');
# 		} elsif( $f->primary_tag eq 'splice3' ) {
# 		    $f->primary_tag('splice_acceptor');
# 		} elsif( $f->primary_tag eq 'frameshift' ) {
# 		    $f->primary_tag('frameshift');
# 		} elsif( $f->primary_tag eq 'cds' ) {
# 		    $f->primary_tag('CDS');
# # 		} elsif( $f->primary_tag eq 'exon' ) {
# # 		    $f->primary_tag('CDS');
# 		}
# 		for my $tag ( qw(intron_id deletions insertions) ) {
# 		    $f->remove_tag($tag) if $f->has_tag($tag); 
# 		}
# 
# if(STRAND_FLIP) {		
# 		if( $f->strand < 0 ) { 
# 		#  dgg; is this correct? NO, see above; causing -start/stop
# 		# 'gene length' from here is segment exonerate ran on not genome length!
# 		# and exonerate1.4 has own swap orient for -strand
# 		    my $topend= $part_start + $length;
# 		    my $s = $topend - $f->end + 1;
# 		    my $e = $topend - $f->start +1;
# 		    $f->start($s);
# 		    $f->end  ($e);
# 		}
# }
# 
# 		if( $gene && 
# 		    $dat{$gene} &&
# 		    $dat{$gene} eq '-' ) 
# 		{
# 		    $f->strand(-1); # dgg: shouldnt this be a flip? -1 if 1; 1 if -1 ?
# 		}
# 		if( $f->primary_tag ne 'gene' && $f->primary_tag ne 'match') {
# 		    $f->add_tag_value('Parent', "$idprefix$gene_name");
# 		}
# 		
# 		if(my $newtag= $renamefeat{$f->primary_tag}) {
# 		  $f->primary_tag($newtag); # gene => mRNA
# 		}
# 		
# 		unless($dropfeat and $f->primary_tag =~ m/$dropfeat/) { 
# 		$out->write_feature($f);
# 		# dont write here; push(@flist, $f) 
# 		# and write gene 1st, with align info attribs
# 		}
# 	    }
# 	    $state = 0;
# 	    $buffer = '';
# 	    %dat = ();
# 	    
# 	   ## next;
# 	} 
# 	
# 	# elsif(/^\#/) { next; }
# 	# orig# if( /^>(\S+)\s+length=(\d+)\s+alnlen=(\d+)/ ) 
# # dgg: 
# # --- END OF GFF DUMP ---
# #
# #qi CG12766-PA length=320 alnlen=285
# #ti scaffold_4786 length=76091 alnlen=1018
# # --- START OF GFF DUMP ---
# 	if( /^#.. (\S+)\s+length=(\d+)\s+alnlen=(\d+)/ ) {
# 	    $lengths{$1} = $2;
# 	    $alnlens{$1} = $3;
# 	    $ct++;
# 	} elsif(/^\w/ && /\t/) {	
# 	    $buffer .= $_;
# 	    #$buffer .= join("\t",split(/\s+/,$_,9));
# 	}
# 	if( $ct == 2 ) { 
# 	    $state = 2;
# 	    $ct = 0;
# 	}
#     } 
#     
#     if( /^\# --- START OF GFF DUMP ---/ ) {
# 	$state = 1; 
# 	$ct =0;
# 	$buffer = '';
# 	
#     } elsif(/^##(source-version|date)/){ 
#         my $com=$1;
# ##source-version exonerate:protein2genome 1.4.0
# ##date 2008-01-12
#         s/#//; 
#         print unless($didcomm{$com}++); # drop one # for gff3
#         # $buffer .= $_ unless($didcomm{$com}++); # no good; parser eats these
#         
#     } elsif(/^#part_location: (\w+):(\d+)-(\d+)/){  # dgg; genome partition header info
#       ($part_ref,$part_start,$part_end)=($1,$2,$3);
#       $part_start-- if($part_start>0);
# 
# ##      # can also parse exonerate Command line for partition info
# # Command line: [/N/u/gilbertd/BigRed/bio/exonerate/bin/exonerate --model protein2genome --minintron 20 --maxintr
# # on 5000 --showtargetgff --showvulgar 0 --showalignment 0 --ryo #qi %qi length=%ql alnlen=%qal\n#ti %ti length=%
# # tl alnlen=%tal\n --query /N/gpfsbr/gilbertd/chrs/aug/dere3/scaffold_4929/scaffold_4929_1-2000000/modDM.fa --tar
# # get /N/gpfsbr/gilbertd/chrs/aug/dere3/scaffold_4929/scaffold_4929_1-2000000/dere_caf060210.fa]
#     } elsif(/^Command line:/){  # dgg; genome partition info  
#       if(/\-target.(\S+)/) {
#         my $t=$1; 
#         my($p)= $t =~ m,([^/]+)/[^/]+$,; # /scaffold_4929_1-2000000/xxxx
#         my($r,$b,$e)= $p =~ m/^(\w+)_(\d+)\-(\d+)$/;
#         if($e>0) { ($part_ref,$part_start,$part_end)= ($r,$b,$e); }
#         }
#       
#     }
# }
# 
# } # bp_exonerate2gff
