#!/usr/bin/env perl
# evigene/namecleangb.pl  from evigene2genbanktbl.pl

use FindBin;
use lib ("$FindBin::Bin","$FindBin::Bin/prot/"); # == evigene/scripts/

use strict;
use Getopt::Long;
use protein_names; # Evigene nameclean() etc

use constant VERSION  => '20130504'; # '20130326';  

## moved to protein_names; all evigene.conf options
# my $NAME_NONE = "Unknown|Uncharacterized|Hypothetical";  
# my $NAME_UNK  = "Uncharacterized protein"; # uniprot  
# my $NAME_UNKADDLOCUS= 1; # policy, add/not the gene id to end of UNK names
# ## my $MIN_ESTIDENT = 10; # min align% to keep  EST/Rna evidence ; not used for naming
# my $MIN_NAMEIDENT = 35;  # min similar% for protein evid align/equivalence, for naming; JCVI uses 35%  
#   # -- note MIN_ID is used for both naming and keeping prot evidence, split this? use lower IDLIKE for evid?
#   # MIN_IDENTITY == MIN_NAMEIDENT now
# my $MIN_IDLIKE   = 15;  # low for now; 15..20 seems right;add Note with loqualname
# my $MIN_CERTAIN  = 0; ## 60;  # putative naming, what?
# my $MIN_PROTIDENT= 0; # $MIN_IDENTITY; .. only for putative, no prot id info here...

our $NAME_KEEPID= 0; # change default
our $USE_TENAME; # in protein_names
our $DO_CONSENSUS; # ?? add to protein_names
my $FILECLEAN=0;
my $DEBUG=0;
my $IDCOL1=1; # =1 default?
my($goodname,$badname);

my $optok= GetOptions(
  "KEEPID!", \$NAME_KEEPID, 
  "TENAMES!", \$USE_TENAME, 
  "FILECLEAN=s", \$FILECLEAN, 
  "consensus!", \$DO_CONSENSUS, 
  "goodname=s", \$goodname, 
  "badname=s", \$badname, 
  "IDCOL1!", \$IDCOL1, # use -noidcol1 to turn off..
  "debug!", \$DEBUG, 
  );

die "usage: app mygene.names > mygene.nameclean 
 opts: -fileclean=uniprot|cdd.names -noIDCOL1 -KEEPID -TENAMES"
 unless($optok);

our $NAME_GOODSPP= $goodname if($goodname);
our $NAME_POORSPP= $badname if($badname);

## sub cleannamefile() handles CDD also .. change to that?
## but specialized for uniprot/cdd naming table w/ id,size cols
if($FILECLEAN) {
  my $infile= $FILECLEAN; # or all of @ARGV;
  my $iscdd= ($infile =~ /cdd/i)?1:0;
  if($infile and -f $infile) {
    my($nrefname,$outnames)= cleannamefile($infile,"",$iscdd) ;
    warn "#cleannamefile(infile=$infile,,iscdd=$iscdd) => ($nrefname,$outnames)\n";
  } else {
    warn "#cleannamefile(infile=$infile,,iscdd=$iscdd) => missing infile\n";
  }
  
} else {  

  my($lid,@names);
  my $idcol=-1; # use same state for all input, resolve name="" problems
while(<>) {
  unless(/^\w/) { print; } #? or bad names?
  chomp;
  #old# my($id,$name)= split"\t",$_,2; # maybe maybenot
  my($id,$name,@other);
  if($IDCOL1) { ($id,$name,@other)= split"\t",$_; $id ||= "noid"; } 
  else { ($name,@other)= split"\t",$_; $id=""; }
  
  # unless($name and $id=~/^[\w.]+/) { 
  #  #OLD# $name=$id; $id="";  # bad for name = "" but id valid
  # } 
  
  if(not $name or $name eq "na") { $name= $NAME_UNK; }
  
  ## pi is not present in group name sets: tair, uniref, omcl .. need to dig out of genes.gff
  # my $pi= ( $name =~ s/\s+\((\d+)%.*\)// ) ? $1 : $MIN_CERTAIN;  # trailing pctident
  ## pi: look in @other for '33%,nnn/nnn,nnn' ?
  my $pi= $MIN_CERTAIN;
  if( $name =~ s/\s+\((\d+)%.*\)// ) { $pi=$1; }
  elsif(@other) {
  	my($po)= grep /^\d+%/, @other; if($po) { ($pi=$po) =~ s/%.*//; }
  }
  
  ##201414 add consensus for multiple names/id ..
  # dang consensname calls nameclean() dont do 2x
  #x my ($newna,$lowqualname,$diff)= nameclean( $name, $pi );
  if($DO_CONSENSUS and $id) {
    if($lid and $lid ne $id) { 
      my($cname,$oldna)= consensname(@names);  
      putna($lid||"",$cname,$oldna); 
      @names=();
    }
    unless($other[0]=~m/^\d+%/){ $other[0]="$pi%"; }
    push @names, join("\t",$name,@other); 
    $lid= $id;
  } else {
    my ($newna,$lowqualname,$diff)= nameclean( $name, $pi );
    putna($id,$newna,$name,\@other,$lowqualname,$diff);
  }
}

  if($DO_CONSENSUS) {
    my($cname,$oldna)= consensname(@names);  
    # my ($newna,$lowqualname,$diff)= nameclean( $cname, $pi );
    putna($lid,$cname,$oldna); 
  }
}

#old..
#   print "$id\t" if($id);
#   print $newna;
#   ##print "\tlowqual=$lowqualname" if($lowqualname); ## FIXME: this adds column; need to preserve input @other
#   print "\t",join("\t",@other) if(@other);
#   print "\t", ($lowqualname ? "lowqual=$lowqualname" : ""); # add 1 col always at end
#   print "\t", (($diff) ?  "update" : "same") if($DEBUG);  
#   print "\told:$name" if($DEBUG and $diff > 1);
#   print "\n";

sub putna {
  my($id,$newna,$oldname,$other,$lowqualname,$diff)= @_;
  print "$id\t" if($id);
  print $newna;
  print "\t",join("\t",@$other) if(ref $other and @$other);
  print "\t", ($lowqualname ? "lowqual=$lowqualname" : ""); # add 1 col always at end
  print "\t", (($diff) ?  "update" : "same") if($DEBUG);  
  print "\told:$oldname" if($DEBUG and $diff > 1);
  print "\n";
}

=item eg input for consensus

Tribca2aEVm000098t1     Brain chitinase and chia/ARP2_G4666     100%,2402/1391,2267     RefID:ARP7f_G2846       nasvit:Nasvi2EG005118t1
Tribca2aEVm000098t1     chitinase 10 isoform X1 60%     RefID:na        tcas3nc:XP_008198138
Tribca2aEVm000098t1     Probable chitinase 3-like Protein       60%     RefID:na        tribcas4a:TC012734tx1

==> names/tribcas4evg2-arpref.names <==
Tribca2aEVm000001t1     Neurogenic locus notch protein  100%,19215/11023,18624  RefID:ARP7f_G2613       dromel:FBgn0053196
==> names/tribcas4evg2-ncbi1406.names <==
Tribca2aEVm000001t1     uncharacterized protein 89%     RefID:na        tcas3nc:XP_008197897

=cut 

sub consensname {
  my(@lnames)= @_;
  my($cname,$oldna)=("","");
  if(@lnames>0) {
    ($oldna)= split"\t",$lnames[0];
    # lnames == nameidlist = "name\tpalign\tid1\tid2"
    ($cname)= consensusname(\@lnames) ; 
    }
  return($cname,$oldna);
}
 
#  # see protein_names::blast2names() and consensusname()
#   sub put3name{ our($outh,$desc1,$ld,@lname,%tdspan); my($cname); 
#     if(@lname>1) { ($cname)= consensusname(\@lname,1,CDD_NOCONSENSUS,\%tdspan); }
#     elsif(@lname==1) { $cname=$lname[0]; } else { $cname=$desc1; }
#     print $outh "$ld\t$cname\n"; }


#.........

=item   moved subs to evigene/scripts/prot/protein_names.pm
  see also nameclean.perl
  sub nameclean {}

=cut

=item add config

# evigene_config($config, \@configadd); # always even if $config null
# evigene_cacao3_gbsubmit.conf
#   nameless    = Unknown|Uncharacterized|Hypothetical  
#   nameunknown = Uncharacterized protein   # uniprot 2011
#   nameunkaddid = 0 # turn off Unc.. locus Thecc11111 addition for gbsubmit complaints
#   nameuncertain = putative # at end; 
#   pctuncertain  = 60  #  min ident% for nameuncertain
#   pctunknown    = 35  #  min similar% for protein evid align/equivalence, for naming; JCVI uses 35% MIN_PROIDENT
#   pctproevidence = 10 #  min similar% for protein evid ; not used for naming 
#   pctrnaevidence = 10 #  min align% to keep  EST/Rna evidence ; not used for naming MIN_ESTIDENT
#   nameidpatt  = (Os|At|AT)\\d{1,2}g\\d{3,} << fix, see below
#   namedrops   = Arabidopsis|thaliana|yeast|complete sequence|complete|genome|pseudogene
 
=cut

