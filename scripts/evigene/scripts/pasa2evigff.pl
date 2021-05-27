#!/usr/bin/env perl
# pasa2evigff.pl

=item about
 
 convert pasa-geneupdate.gff3 to evigene format
  -- move protein comment to mRNA protein= attribute
  -- rename alternate transcript IDs, "xxxx_t1.1" > "xxxx_t2"
  -- rename exon,CDS ID= to xid=, leaving Parent=mRNAid; 
      ID=parent.exon1 .. parent.cds1
  -- drop five,three_prime_utr feat, unless requested
  -- drop/comment gene feature, unless requested; change mRNA Parent=geneid to gene=geneid
  
  -- handle MERGE genes: ids joined : leave as is?
    scaffold00512   .       gene    103964  124725  .  
        ID=AUGep24b_p1s00512g32t1_AUGep24b_p1s00512g33t1;
        Name=MERGED%3A%20AUGep24b_p1s00512g32t1%3B%20AUGep24b_p1s00512g33t1

=cut

use strict;
use Getopt::Long; # qw(:config no_ignore_case bundling);

my $exp= $ENV{n} || "pasaupd";
my $idprefix= $ENV{idprefix} || ''; 
my $aa_attribute=1;
my $uniqueid= 1; #?

my $dropft="gene|intron|five_prime_utr|three_prime_utr|start_codon|stop_codon|transcription_start_site|transcription_end_site";
my $keepft="exon";

my $optok= &GetOptions (
    'source|exp=s' => \$exp,  
    'idprefix=s' => \$idprefix,
    'keepft=s' => \$keepft, 
    'dropft=s' => \$dropft, 
    'aa_attribute|protein!' => \$aa_attribute,
    #'uniqueid!' => \$uniqueid,
    );
    
die "usage: pasa2evigff  < pasa_geneupdate.gff > evigene_update.gff 
  options:  -source=$exp --noprotein
\n" unless($optok);

my @drops= grep { $_ !~ m/$keepft/} split(/[,\| ]/,$dropft);
my $drops= join("|", @drops) or "NADA";

my ($gid,$aa, $ismergedgene, $isoriginal) = (0) x 10;
my (@gff,@trids,%idmap,%prots,%status);
my ($ingene, $incon, $oldid);

sub cleanid { 
  my ($id,$altflag)= @_; # shift; 
  #? handle MERGE ? ID1_ID2
  if($id =~ m/t(\d+)\.(\d+)$/) {
    my ($oldalt,$pa)= ($1,$2);
    my $alt= $oldalt + $pa;
    $id =~ s/t$oldalt\.$pa/t$alt/;
  } 
  elsif( $altflag and $id =~ m/\.(\d+)$/) { ## handle; cufflinks ID problem: xxx.nnn.1 = 1st tr, dont bump to t2
#   orig md8aphid_cuf8r27Gsc278.10918.1 > 18t2
#   pasa md8aphid_cuf8r27Gsc278.10918.1.1 > 18.1t2 
    my $alt= 1 + $1;
    $id =~ s/\.(\d+)$/t$alt/;
  }
  # else{ $id=~s/\W/_/g; } #?? bad for orig like cuf.n.n

  $id=$idprefix.$id if($idprefix);
 
  return($id); 
} 

sub putgene {
  
  # @gff contains all mRNA transcripts for gene
  # my @trids= sort keys %status unless(@trids); # or always get here?
  
  foreach my $trid (@trids) {
    my @attr=();
    my $stat= delete $status{$trid};
    push(@attr, "upstatus=$stat") if($stat);
    
    my $tlen=0; 
    foreach my $ex (grep /\texon/, @gff) {
      next unless($ex =~ /Parent=$trid/);
      my($b,$e)=(split"\t",$ex)[3,4]; $tlen += 1+$e-$b; 
      } 
    $tlen ||= 1;
    
    my $aa= delete $prots{$trid}; 
    my $al= length($aa)||0; $al-- if($aa =~ m/\*$/);
    my $pt= int($al*300/$tlen);
    push(@attr,"aalen=$al,$pt\%"); # always 
    push(@attr,"protein=".$aa) if($aa and $aa_attribute); 
    
    if(@attr) {
      my $attr= join(";",@attr);
      #foreach my $tr (grep /\tmRNA/, @gff) { $tr =~ s/$/;$attr/ if($tr =~ /ID=$trid\b/); } 
      foreach my $tr (grep /\tmRNA/, @gff) { $tr =~ s/$/;$attr/ if($tr =~ /ID=$trid[;\s]/); } 
      }
      
    }
  
  # @gff=() if $fastaonly;
  print @gff if(@gff); 
  
  @gff= @trids= (); 
  #? %prots= %status=();  
  ($aa)= (0) x 20;
}

my $trid; #? or not global
my $altflag=0;
$ingene=1;
print "##gff-version 3\n";
while(<>){

  if(/^#/) {
    my $p=1;
    
# PASA_UPDATE: AUGep24b_p2s00512g106t1, 
#   single gene model update, valid-1, status:[pasa:asmbl_433,status:12],[pasa:asmbl_1224,status:12],[pasa:asmbl_2809,status:12], valid-1
# PASA_UPDATE: AUGep24b_p2s00512g106t1.1,
#    alt-splice addition, valid-1, status:[pasa:asmbl_432,status:25], valid-1
# PASA_UPDATE: AUGep24b_p2s00512g106t1.2, alt-splice addition, valid-1, status:[pasa:asmbl_434,status:25], valid-1

    if(/# PASA_UPDATE:\s+([^\s,]+),\s*(.+)$/) {  
      my($trid,$upd)=($1,$2);
      putgene() if(@gff);

      $altflag=(/alt-splice/)?1:0; # no good here, all # before all mRNA, mix of main/alt
      my $newtrid= cleanid($trid, $altflag);
      $idmap{$trid}= $newtrid; $trid= $newtrid;
      $isoriginal= 0;
      my ($uptype, $val, @upd)= split", *", $upd;
      my %stat;
      foreach my $up (@upd) { if($up =~ m/status:(\d+)/) { $stat{$1}++; } }
      my $stat= "codes:" . join "/", sort{ $stat{$b} <=> $stat{$a} } keys %stat;
      $status{$trid}= "$uptype/$val/$stat";
    }
    elsif(/# ORIGINAL:\s+([^\s,]+)(.+)$/) { 
      my($trid,$upd)=($1,$2);
      putgene() if(@gff);

      $altflag=0;
      my $newtrid= cleanid($trid, $altflag);
      $idmap{$trid}= $newtrid; $trid= $newtrid;
      $isoriginal= 1;
      $status{$trid}= "original";
      # change source to what if this is set?
    }
    
      # this comes AFTER gene gff; all in 1 line?
    elsif(/#PROT (\S+)\s+(\S+)\s+(\S+)/) {  
      my($trid,$gid,$aa)=($1,$2,$3);  $p=0; 
      ## $trid= cleanid($trid,$altflag);
      $trid= $idmap{$trid} || cleanid($trid); 
      $prots{$trid}= $aa if($aa); $aa="";
      } 
    # elsif($aa and /# ([A-Z]+)/) { $aa.=$1; $p=0; }  
    elsif(/##gff-version/){ $p=0; } # not seen, do self
    
    print if($p); ## ( and not $fastaonly);
    
  } elsif(/^\w/){
    my $p= 1;

    if($ingene) {
    
      s/^(\S+)\t\./$1\t$exp/; # source col2  
      
      #my $prefix= $idprefix;
      #in.cleanid# s/(ID|Parent)=/$1=$idprefix/g if($idprefix);
      
      if(m/\tgene\t/){ 
        putgene() if(@gff);
        $ismergedgene= (m/Name=MERGED/)?1:0;  # ** need merge info
        
      } elsif(m/\tmRNA\t/) { 
        ## put_transcript() if(@gff);
        
        s/;Parent=([^;\s]+)//; $gid=$1; # drop gene id or not?
        ($trid)= m/ID=([^;\s]+)/; $oldid= $trid;
        
        ## $trid= cleanid($trid, $altflag); #? here or on output
        $trid= $idmap{$trid} || cleanid($trid); 

        push @trids, $trid;
        
#         if($uniqueid && $seengene{$trid}++) { 
#           my $did=  "d" . $seengene{$oldid};
#           unless($trid =~ s/$prefix/$prefix$did/) {
#             $trid = $did.$trid;
#           }
#           #s/=$oldid\b/=$trid/g;
#         }
        
        s/=$oldid\b/=$trid/g unless($trid eq $oldid);
        
      } else {  # exon, CDS
        my($parid)= m/Parent=([^;\s]+)/; # use this or use last mRNA ?
        if($oldid && $trid ne $oldid) {  s/=$oldid\b/=$trid/g;  }
        if(s/\bID=([^;\s]+);?//) { my $xid=$1; $xid=~s/^$trid\.//; s/$/;xid=$xid/; } # move to last
        }
        
      $p= m/\t($drops)/ ? 0:1; 
      }
      
    # $p=0 if(/^Error/);

    push(@gff,$_) if $p; 
  }

}

putgene() if(@gff);


__END__

=item sample pasa.gff

# PASA_UPDATE: AUGep24b_p2s00512g106t1, single gene model update, valid-1, status:[pasa:asmbl_433,status:12],[pasa:asmbl_
1224,status:12],[pasa:asmbl_2809,status:12], valid-1
# PASA_UPDATE: AUGep24b_p2s00512g106t1.1, alt-splice addition, valid-1, status:[pasa:asmbl_432,status:25], valid-1
# PASA_UPDATE: AUGep24b_p2s00512g106t1.2, alt-splice addition, valid-1, status:[pasa:asmbl_434,status:25], valid-1
scaffold00512   .       gene    2676970 2679148 .       -       .       ID=AUGep24b_p2s00512g106t1;Name=AUGep24b_p2s00512
g106t1
scaffold00512   .       mRNA    2676970 2679148 .       -       .       ID=AUGep24b_p2s00512g106t1;Parent=AUGep24b_p2s005
12g106t1
scaffold00512   .       five_prime_utr  2678864 2679148 .       -       .       ID=AUGep24b_p2s00512g106t1.utr5p1;Parent=
AUGep24b_p2s00512g106t1
scaffold00512   .       exon    2678601 2679148 .       -       .       ID=AUGep24b_p2s00512g106t1.exon1;Parent=AUGep24b_
p2s00512g106t1
...

scaffold00512   .       mRNA    2676970 2679148 .       -       .       ID=AUGep24b_p2s00512g106t1.1;Parent=AUGep24b_p2s0
0512g106t1
scaffold00512   .       five_prime_utr  2678864 2679148 .       -       .       ID=AUGep24b_p2s00512g106t1.1.utr5p1;Paren
b_p2s00512g106t1.1
scaffold00512   .       CDS     2678601 2678863 .       -       0       ID=AUGep24b_p2s00512g106t1.1.cds1;Parent=AUGep24b
_p2s00512g106t1.1
scaffold00512   .       exon    2678340 2678538 .       -       .       ID=AUGep24b_p2s00512g106t1.1.exon2;Parent=AUGep24
b_p2s00512g106t1.1
...
scaffold00512   .       mRNA    2676970 2679148 .       -       .       ID=AUGep24b_p2s00512g106t1.2;Parent=AUGep24b_p2s0
0512g106t1
scaffold00512   .       five_prime_utr  2678601 2679148 .       -       .       ID=AUGep24b_p2s00512g106t1.2.utr5p1;Paren
t=AUGep24b_p2s00512g106t1.2
scaffold00512   .       five_prime_utr  2678529 2678542 .       -       .       ID=AUGep24b_p2s00512g106t1.2.utr5p2;Paren
t=AUGep24b_p2s00512g106t1.2
scaffold00512   .       exon    2678601 2679148 .       -       .       ID=AUGep24b_p2s00512g106t1.2.exon1;Parent=AUGep24
b_p2s00512g106t1.2
scaffold00512   .       exon    2678340 2678542 .       -       .       ID=AUGep24b_p2s00512g106t1.2.exon2;Parent=AUGep24
b_p2s00512g106t1.2
scaffold00512   .       CDS     2678340 2678528 .       -       0       ID=AUGep24b_p2s00512g106t1.2.cds2;Parent=AUGep24b

...
        mRNAid
#PROT AUGep24b_p2s00512g106t1 AUGep24b_p2s00512g106t1   MNLLSESISP..*

#PROT AUGep24b_p2s00512g106t1.1 AUGep24b_p2s00512g106t1 MNLLSES..*

#PROT AUGep24b_p2s00512g106t1.2 AUGep24b_p2s00512g106t1 MLDGG..




...
# PASA_UPDATE: AUGep24b_p2s00512g148t1, single gene model update, valid-1, status:[pasa:asmbl_459,status:12],[pasa:asmbl_
700,status:12],[pasa:asmbl_710,status:13], valid-1
scaffold00512   .       gene    2988049 2991675 .       -       .       ID=AUGep24b_p2s00512g148t1;Name=AUGep24b_p2s00512
g148t1
scaffold00512   .       mRNA    2988049 2991675 .       -       .       ID=AUGep24b_p2s00512g148t1;Parent=AUGep24b_p2s005
12g148t1
scaffold00512   .       five_prime_utr  2991521 2991675 .       -       .       ID=AUGep24b_p2s00512g148t1.utr5p1;Parent=
AUGep24b_p2s00512g148t1
scaffold00512   .       exon    2991445 2991675 .       -       .       ID=AUGep24b_p2s00512g148t1.exon1;Parent=AUGep24b_
p2s00512g148t1
scaffold00512   .       CDS     2991445 2991520 .       -       0       ID=AUGep24b_p2s00512g148t1.cds1;Parent=AUGep24b_p
2s00512g148t1
scaffold00512   .       exon    2991282 2991374 .       -       .       ID=AUGep24b_p2s00512g148t1.exon2;Parent=AUGep24b_
p2s00512g148t1


#PROT AUGep24b_p2s00512g148t1 AUGep24b_p2s00512g148t1   MASPAKS....

=cut
