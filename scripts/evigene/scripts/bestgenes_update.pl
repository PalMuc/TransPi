#!/usr/bin/env perl
# bestgenes_update.pl 

=item about

 evigene/bestgenes_update.pl or other name
  : from bestnotrivial.sh, addcheckprot.sh mustkeepgenes.sh


=item usage

  input is various:
    
=item configuration

  puboption:
    publicid = EG2AP0000000 # or CG00000000 # or mygene0000000
  
    keepsource = acypi|ars17trinity|ars27cuf8|pasa2_aphid3|ref_aphid2
    # ^ this really means keep from trivial reaper as all have rna-seq evidence
    keepids = bestgenes.DGILmix8d.addback1.ids

=cut

use constant VERSION => '2018.04.16'; # misc updates
# '2013.08.31'; # ... way back 

use FindBin;
use lib ("$FindBin::Bin"); # assume evigene/scripts/ layout: main, prot/, rnaseq/, ...
our $EVIGENES="$FindBin::Bin";  

use strict;
use Getopt::Long;

my $overlapfilter="$EVIGENES/overlapfilter";
my $overlapeval=""; # not here
my $overbestgenes="$EVIGENES/overbestgene2.perl";

my $debug=1;
my @annotkeys= (); # qw(est pro rseq ref intr pasa tar terepeat); # config OPTION
my @evalkeys=();
my (@moreannotkeys,@moreevalkeys);

my %evidence= ();
my %evaluate_options = ();
my %general_config = ();
my %annotate_options = ();
my %public_options = ();
my %geneset=(); # from config
my %programs= (); # add to other evigene set

my $BINSIZE   = 5000 ; 
my $mrnatypes=$ENV{'mrnatypes'} || 'mRNA|ncRNA'; # upd18 to ncRNA
my $exontypes='exon|CDS';
my $overlaplist= {};

my $config;
my @configadd;
my $action;

my $pubidnum_start=0;
my $ALTKEY="t"; # but see config altid_format
my $MID=""; # make-id = temp/prelim bestgenes id prefix
my $MSRC="evigene"; ## DGILmix8d == anv ?
my ( $ingff, $outgff )= ("","");

#my $annotate_genes = 1; # what??
#my $make_bestgenes= undef;
#my $rescore_only= 0;# = overbestgenes -rescore

my $optok= GetOptions(
  "config=s", \$config,
  "action=s", \$action,
  "input=s", \$ingff,
  "output=s", \$outgff,
  "mrnatypes=s", \$mrnatypes,
  "MID=s", \$MID,  # change/extend this mean ID prefix, evd3ap0000000, to replace config:publicid
                   # usage?  MID=mx1 == prefix for current ID; MID=mx0000000 == replacement ID patt
  "pubidnum=i", \$pubidnum_start,
  "version|MSRC=s", \$MSRC,
  "cadd=s", \@configadd,
  "debug!", \$debug, 
  );
  
$ingff= shift @ARGV unless($ingff);

die "usage: bestupdate -action=xxx -conf=evigene.conf  in.gff ...
  action=trivial_remove,mustkeep_update,pasa_update,protein_add,sortgenes,public_update
  opts: -version $MSRC -out out.gff -mid prefixID ... \n"
  unless($optok and $action and (-f $ingff or $ingff =~ /stdin|-/));

evigene_config($config, \@configadd); # always even if $config null

# new usage  -MID=mx2 is prefix for current ID; -MID=mx2000000 is replacement ID patt
if($MID and $MID =~ /^[a-zA-Z]\w*0000$/) { 
  $public_options{'publicid'} = $MID; $MID=0;
}

update_genes($ingff, $outgff);

#.......................................................


sub update_genes 
{
  my ( $ingff, $outgff )= @_;
  
  my $genedir="";
  my $grp= $ingff;
  if($grp =~ s,^(.*/)([^/]+)$,$2,) { $genedir=$1;  }
  $grp =~ s,\.gz,,; 
  $grp =~ s,\.gff,,; 
  # $grp =~ s/\W?augmap//;  #?
  $grp =~ s,\.an\w+.*,,; #? or not; need to make sure have uniq grp names
  
  $outgff="${genedir}$grp.$MSRC.gff" unless($outgff);
  
  my @tmpgff=($ingff);
  my $tmpgff= $ingff;

  # open(my $inh, $ingff) or die "reading $ingff\n";  
  # open(my $outh, ">$tmpgff") or die "writing $ingff\n";
  
  if($action =~ /trivial/) {
  $tmpgff= remove_trivial($tmpgff[-1]);
  push @tmpgff, $tmpgff if(-f $tmpgff);
  }

  # FIXME: mustkeep
  #  .. should be able to add new must to public.gff, 
  #  .. new pubids, hide/drop replaced loci
  #  .. install/copy locus old-new annotate if valid
  #  .. handle must alt-tr properly
  if($action =~ /mustkeep/) {
  $tmpgff= merge_mustkeepgenes($tmpgff[-1]); # needs mustgenes.gff
  push @tmpgff, $tmpgff if(-f $tmpgff);
  }
  
  if($action =~ /pasa/) {
  $tmpgff= merge_pasaupdate($tmpgff[-1]);  # needs pasaupd.gff
  push @tmpgff, $tmpgff if(-f $tmpgff);
  }

  if($action =~ /protein/) { # protein_add
  $tmpgff= add_missprotein($tmpgff[-1]); # FIXME: new soft, need test/report on prot qual.
  push @tmpgff, $tmpgff if(-f $tmpgff);
  }

  ## FIXME: add generec checks as per evigene/scripts/evigene2genbanktbl.pl  
  # warno "#ERROR: mrna missing gene at ",$dat,"\n"; 
  # warno "#WARN: fix missing strand: gene $gid\n";  
  # warno "#ERROR: mrna ne gene =$er, at ",$dat,"\n";
  # if($gene->[6] eq ".") { $gene->[6]= $mrna->[6]; warno "#WARN: fix missing strand: gene $gid\n"; }  
  # warno "#ERROR: $id CDSexon, exon ne mRNA: $cdserr\n" if($cdserr);
  # warno "#WARN: $tid has wrong gene=$pgene; should be $pgenet\n";
  # warno "#WARN: $tid is wrong? gene=$pgene; not $pgenet\n";
  # warno "#WARN: $tid strandfix: $onew \n";
  # warno "#WARN: $tid spanfix: $loc \n";
  # warno "#WARN: $tid rename $ntid\n";
  # unless(@generec) { warno "#NOTE: skipped exon, no parent $typ=$pid\n"; }  
  # if($rid ne $pid and not $IgnoreOutOfOrder) { warno "#ERROR: out of order gene record: mRNA=$rid, $typ=$pid\n" if($nwarn++ < 9); } # limit warns
  
  #x add gene record enclosing mRNAs; 
  if($action =~ /sort/) { # adds generec if config
  $tmpgff= genesort_update($tmpgff[-1]);
  push @tmpgff, $tmpgff if(-f $tmpgff);
  }
  
  # change mRNA to ncRNA if flagged
  #x improved alt-tr info
  if($action =~ /public/) {
  $tmpgff= public_update($tmpgff[-1]);
  push @tmpgff, $tmpgff if(-f $tmpgff);
  }

  if($action =~ /annotate/) {
#  $tmpgff= public_annotate($tmpgff[-1]); # see bestgenes_puban.pl
#  push @tmpgff, $tmpgff if(-f $tmpgff);
  }
  
  my $lastgff= pop(@tmpgff);
  if($lastgff ne $ingff and -f $lastgff) {
  print "# rename $lastgff $outgff\n" if $debug;
  # rename genes/bestgenes.DGILpub8d3.gff, genes/bestgenes.DGILpub8m3.gff-vers  
  if( -f $outgff) { rename( $outgff, "$outgff.old"); }
  rename($lastgff, $outgff);
  }
  
  foreach my $tmp (@tmpgff) {
    if($tmp and -f $tmp and $tmp ne $ingff and $tmp ne $lastgff) { # and not $debug ??
      if($debug) { print "unlink $tmp\n"; } else { unlink $tmp; }
    }  
  }
  
}



sub openio {
  my($ingff, $outgff)= @_;
  my($inh, $outh, $ok)= (undef) x 3;
  # ADD input gff.gz, stdin ..

  if($ingff =~ /.gz/) { $ok= open($inh, "gunzip -c $ingff |"); }
  elsif($ingff =~ /stdin|^-/) { $inh= *STDIN; $ok=1; }
  else { $ok= open( $inh, $ingff); }
  die "ERROR: read $ingff" unless($ok);

  if(!$outgff or $outgff =~ /^stdout|^-/) { $outh= *STDOUT; $ok=1; }
  else { $ok= open( $outh, ">$outgff"); }
  die "ERROR: write $outgff" unless($ok);

  return($inh, $outh);
}

sub setprogram {
  my($evcmd, $defname, $defprog)= @_;
  my $progna=(split " ",$evcmd)[0]; 
  my $prog= $defprog;
  if($progna and $prog= $programs{ $progna }) {  }
  elsif($defname and $defprog) { $progna=$defname; $prog= $defprog;  }
  $evcmd =~ s/$progna/$prog/;
  unless( -x $prog) { die "ERROR: missing program $progna: $prog\n"; } # or die
  return $evcmd;
}

sub merge_pasaupdate 
{
  my($ingff, $outgff)= @_;
  
  # **FIXME ** need to copy exons also as some pasa updates change exons
  #   == AUGepir16bs950g39t1 pasa code:16/EST assembly properly stitched into gene structure
  #   .. but changed exon/CDS locations
  
  # FIXME: from config or options
  my $addgff= $evidence{"pasa_updates"}; 
  unless( $addgff and -f $addgff) {
    print "# missing pasa_updates\n"; return $ingff;
  }
  
  # do this pro=  mark before or here? use evigene.conf anoption settings
  # .. this remarking here is a hassle, now have to do exon evid marks also : leave aside
  # .. do that separately to pasaupd.gff

use constant  EVMARK_PASAUPD => 0;

  my $cmd="";
if( EVMARK_PASAUPD ) {
  my $pev= "pro";
  my $evfile  = $evidence{$pev};
  my $evcmd   = $annotate_options{$pev};
  $evcmd =~ s/overlapfilter/$overlapfilter/;
  $evcmd= ($evcmd and -f $evfile) ? "| $evcmd -mark $pev -in stdin -over $evfile" : "";
  # overlapfilter -strand -pass CDS -pct 10 -act markidbase -mark pro -in stdin -over prot/protein_uniq.gff.gz

  $cmd= ($addgff =~ /\.gz/) ? "gunzip -c " : "cat ";
  ## $cmd .= " $addgff | egrep '(mRNA|CDS)' $evcmd |"; # no, pass all types, keep for alttr
  $cmd .= " $addgff $evcmd |";
} else {
  $cmd= ($addgff =~ /\.gz/) ? "gunzip -c " : "cat ";
  $cmd .= " $addgff |";
}

  use constant { PACODE_SKIP => 1, PACODE_ORIGINAL => 2, PACODE_NOCHANGE => 12, 
    PACODE_UTRCHANGE => 13, PACODE_CDSCHANGE => 14, PACODE_OTHERCHANGE => 16,
    PACODE_ALT => 25,  PACODE_MERGE => 36, PACODE_MERGE_PART => 63, };

  my( %mattr, %mstat, %mcds, %mgene, %malt, %mmerge, $paskip, $pacode, $nalt, $nupd);
  print "# merge_pasaupdate: $cmd\n" if($debug);
  open(P,$cmd) or die "ERR: $cmd";
  while(<P>) {
    
    unless(/^\w/) { 
      # nada
      
    } elsif(/\tmRNA/) { 
      my($id)=m/ID=([^;\s]+)/; 
      my($at)=(split"\t")[-1];  chomp($at); $at =~ s/ID=[^;\n]+;//;
      
      $paskip= $pacode= 0;
      unless( $at =~ /protein=/ ) { $pacode= $paskip= PACODE_SKIP; next; } # skip here
      if( m/status=original/ ) { $pacode= $paskip= PACODE_ORIGINAL; } # next if($paskip); #no, save protein=
      
      unless($pacode) {
        my ($pacodes)= m,status=[^/]+/valid-./codes:([^;\s]+),;
        ($pacode) = sort { $b<=>$a } split( /[,\/]+/, $pacodes); # pick highest
      } 
      
      $paskip=PACODE_NOCHANGE if($paskip == 0 and $pacode <= PACODE_NOCHANGE);
      
      
#      ## FIXME: merged genes need to be added, replacing un-merged ones.  e.g. this horror, but right
# # ID=md8aphid2Trinity39373p1_md8aphid_cuf8r27Gsc122.2105.1_md8AUGepi4s122g62t1_md8AUGepir2s122g46t1_md8AUGepir9s122g60t1
#  single gene model update, valid-1, status:[pasa:asmbl_22820,status:36], valid-1
#  ID=aphid_cuf8r27Gsc36.18012.1_md8aphid_cuf8r27Gsc36.18012.2_md8AUGepir3s36g36t1;
# upstatus=single gene model update/valid-1/codes:36
#  codes:36 == merge
#   count code    desc
#   26348 12      EST assembly incorporated.     << no change
#    7020 13      EST assembly extends UTRs.     << change exons
#    1922 14      EST assembly alters protein sequence, passes validation. < change CDS,exons,prot
#    1140 16      EST assembly properly stitched into gene structure.   < can change CDS,exons
#    8726 25      EST assembly stitched into Gene model requires alternative splicing isoform.  < alttr
#     424 36      EST-assembly found capable of merging multiple genes.  < replace mRNA,CDS,exons old gene recs
#      # also look for status=novel ??

      # check codes: for no-change validate
      # my $ismerge=0;
      if( $pacode == PACODE_MERGE ) {
        (my $td= $id) =~ s/(aphid|asmbl)_(\w)/${1},${2}/g; # aphid_cuf, asmbl_ .. have to parse messy ids
        my @mergeids= split(/_/, $td); # BUT asmbl_nnnn
        # $ismerge= @mergeids;  # fix MESSY ids? no, leave to pubid
        map { s/,/_/; $mstat{$_}=PACODE_MERGE_PART;  $mmerge{$_} = $id } @mergeids;
      }
      
      my $isalt= ($pacode == PACODE_ALT or m/status=alt-splice/) ? PACODE_ALT : 0;  # save this, if not isalt, reset any existing alttr=1 flag
      if($isalt) {
        (my $id0=$id) =~ s/t\d+$//;  ## or s/\.[2-9]\d*$//;
        $malt{$id0} .= "$id,";
        }

      $mattr{$id}= $at; 
      $mgene{$id}= $_;  
      $mstat{$id}= $pacode; ## ($paskip > 0) ? $paskip : ($ismerge>0)?5: ($isalt>0) ? 4 : 3;
      
    } elsif(/\t(CDS|exon)/ and $paskip == 0) { # $paskip == $pacode <= PACODE_NOCHANGE
      my($p)=m/Parent=([^;\s]+)/; s/ID=[^;]+;//; #? drop xid=exon1|cds1 ?
      $mcds{$p} .= $_ if(/\t(CDS|exon)/);  
      $mgene{$p}.= $_;
    } 
    
  } close(P);

  
  my($at,$id, $mergedrop, @addalt);
  $outgff="$ingff.pa" unless($outgff);
  my($inh, $outh)= openio($ingff, $outgff);
  
  while(<$inh>) {
  
    if(/^\w/ and /\tmRNA/) { 
      if(@addalt) {  $nalt += @addalt; print $outh @addalt; @addalt=(); }
      ($id)=m/ID=([^;\s]+)/; $at=""; 
      
      $mergedrop=0;
      if( my $idmerge= $mmerge{$id} ) {  # and $mstat{$id} == PACODE_MERGE
        my $newgene= $mgene{$idmerge};  
        if($newgene) {
          my $didmerge= ($mattr{$idmerge} == PACODE_MERGE)? 1 : 0;
          push @addalt, $newgene unless($didmerge); ## $_= ($newgene) ? $newgene : ""; 
          $mattr{$idmerge}= PACODE_MERGE;
          $_= ""; # drop this
          $mergedrop= 1;
          }
        
      } elsif( $at= delete $mattr{$id} ) { 
        my ($aold)= m/protein=([^;\n]+)/;  $aold||="";
        my ($anew)= $at =~ m/protein=([^;\n]+)/; 
        
        my $lold=length($aold)||0; $lold-- if($aold =~ m/\*$/);
        my $lnew=length($anew)||0; $lnew-- if($anew =~ m/\*$/);                      
        my $ldif= $lnew - $lold; $ldif="+$ldif" if($ldif>0); 
        
        #? dont update if ldif == 0 and $aold eq $anew -*
        # do update if aold protein= is missing but found in status=original
        if($mstat{$id} > PACODE_NOCHANGE or ($lold == 0 and $lnew > 0)) {        
          s/;?(aalen|protein|checkprot)=[^;\n]+//g; # drop cxlen or not ? only if readd here
          $at .= ";checkprot=$ldif";
          s/$/;$at/;   $nupd++;
          
          ## need to add in CDS and exon, many have exon changes;
          my $cd= $mcds{$id};
          if($cd) { 
            my($o)=(split"\t")[6];
            if($o eq "-") { chomp($cd);  my @cd=reverse split"\n",$cd; $cd=join"\n",@cd; $cd.="\n"; }
            $_.= $cd; 
            }
          }
        } 
        
      if(my $hasalt= delete $malt{$id}) {
        $hasalt =~ s/,$//; 
        s/;?alttr=1//g; # a few of these; primary == alt mixup
        foreach my $alt (split",",$hasalt) {
          push @addalt, $mgene{$alt} if $mgene{$alt};          
        }
      }
      
    } elsif(($at or $mergedrop) and /^\w/ and /\t(CDS|exon)/) { # need to drop exon also (if changed)
      my($p)=m/Parent=([^;\s]+)/; 
      #?? if($p ne $id){ warn "misorder CDS.$p ne $id\n"; } 
      $_ = "" if($mergedrop or $mcds{$p});  # or '#x' comment out ?
    }
    
    print $outh  $_; # OUT; to where ?? should be pipe-able
  } close($inh);
  
  if(@addalt) {  print $outh @addalt; @addalt=(); }

  print $outh "\n# merge_pasaupdate for $MSRC:  updated: $nupd, added alt-tr: $nalt\n";
  if($debug) {
  print       "\n# merge_pasaupdate for $MSRC:  updated: $nupd, added alt-tr: $nalt\n";
  }
  
  close($outh);  
  return $outgff;
}


# FIXME: mustkeep
#  .. add new must to public.gff, as well as pre-public bestgenes input
#  .. new pubids : expect to run this after "final" complete gene set
#  .. handle hide/drop of replaced loci
#  .. install/copy locus old-new annotate if valid
#  .. handle must alt-tr properly
#  .. ?? create complete file set mustkeep.gff,.aa,tr.fasta for merge
#  .. what of gene annots, dont expect same as full set: homolog=, express=, quality=

# **** for config
sub mustkeep_files
{
  my $hasgeneset= (scalar(keys %geneset) > 0)?1:0;
  # my $genetab  = $evidence{"mustkeepdrop"}; # change from processed ids to pickgenes.results
  my $genetab   = $evidence{"mustinput"} || "";  # public_options{mustxxx}
  ##^ replace with config/geneset:
  
  unless($genetab or $hasgeneset) {
  $genetab = <<"EOT";
AUGepi4     genes/aphid2_epi4.an8.gff 
AUGepir2    genes/aphid2_epir2.an8.gff 
AUGepir16b  genes/aphid2_epir16b.an8.gff 
augadd8     genes/aphid2_add4best.an8.gff
AUGepir3    genes/aphid2_epir3.an7.gff 
AUGepi5     genes/aphid2_epi5.an6.gff 
AUGepir9    genes/aphid2_epir9.an6.gff 
AUGepir10   genes/aphid2_epir10.an6.gff 
acyp2eg     genes/bestgenes.DGILpub8d3.gff  
aphid_cuf8    rnagene/aphid_rnaseq.all27cuff8mecdso.an8.gff  
aphid2Trinity rnagene/aphid_trinity.ident1.an2.gff
PASAgasmbl    epasa3/pasa_outr/pasa2_aphid3.asmbl_genes.gff.gz
aphid_vel7    rnagene/aphid_vel7asm.ident1.addprot.gff 
ACYPI       genes/acyr1-ACYPImRNA.an7.gff
XM_         refseq/ncbi2/acyr2_ncbirefgene.gff 
NM_         refseq/ncbi2/acyr2_ncbirefgene.gff 
EOT
  }

# PASAgasmbl epasa3/pasa_outr/pasa2_aphid3.asmbl_genes.gff.gz
# PASAgasmbl rnagene/pasa2_aphid3.asmbl_bestgenes.an8.gff  : missing many choices above
#..
#  aphid2Trinity rnagene/aphid_trinity.ident1.an2.gff
#  aphid2Trinity rnagene/aphid_trinity.ident1.an8.gff  ## missing many from .an2

  my %genefiles = %geneset;
  foreach my $in (split( /\n/, $genetab) ) {
    $in =~ s/^\s+//; $in =~ s/\s+$//;
    my($k,$v)= split /[=\s]+/, $in, 2;
    $genefiles{$k}= $v if($k and $v);  # unless( -f $v) { warn ... }
  }  
  return %genefiles;
}




sub make_mustgenes
{
  my (%gact,%gord, %act, %have, %hgene, %found, %miss, $nh, $nm, %nout, $nout);
  my %genefiles= mustkeep_files();
  my @idp= sort keys %genefiles;
  my $idp= join("\|",@idp); 
  
  # my @idp=qw(AUGepi4 AUGepir2 AUGepir16b AUGepir3 AUGepi5 AUGepir9 AUGepir10 acyp2eg); 
  # my @idr= qw(aphid_cuf8 aphid2Trini PASAgasmbl aphid_vel7  ACYPI XM_ NM_); 
  # my $idp= join("\|",@idp); 
  # my %idi; my $i=0; map{ $idi{$_}=$i++; } @idp;
  
  my $mustfile  = $evidence{"mustkeeptab"}; # was {"mustkeepdrop"}; # change from processed ids to pickgenes.results
  my $mustout   = $evidence{"mustgff"} || "mustgenes.gff";   

# gpick=  173 altmodel 426 best 224 drop 1 none 4 skiplocus

## fixme ids:
#  mustlist= ACYPI29102-RAp1  << where from p1
#  acypi.gff = ACYPI29102-RA

# mustfile == map/pickgene-aphid2x.results
# Mon/May/30/16:18:43/2011      69.136.4.212    org=aphid2x     gid=PASAgasmbl_9770     gpick=best      loc=Scaffold107:438031-882081

  my $iord=0;
  print "# make_mustgenes from $mustfile\n" if($debug);
  unless( $mustfile and  open(F,$mustfile) ) { warn "# Missing $mustfile\n"; return; }
  
  while(<F>) { 
    my($date,$ipaddr,$org,$gid,$act,$loc); # split"\t" or not?
    ($date)= (split)[0]; $date =~ s,^\w+/,,; #? save to check update w/ mustout, parse date?
    ($gid)= m/\tgid=(\S+)/;
    ($act)= m/\tgpick=(\S+)/;
    ($loc)= m/\tloc=(\S+)/;  

    $gid =~  s/^[xm]\w+($idp)/$1/;  # remove mix prefix
    $gid =~  s/\-RAp1/-RA/ if($gid =~ /ACYPI/);  # remove acypi mismatch

    # handle act = skiplocus, none = no best model ?
    # also need locus test for best2 following best1 => drop best1
    
    $gact{$gid}= $act; # "replace any old";
    $gord{$gid}= ++$iord; # record input order, last is best
    } close(F);
    
  foreach my $gid (sort keys %gact) { 
    my $act= $gact{$gid}; $act{$act}++; 
    if($gid =~ m/^($idp)/) { 
      my $p=$1;  
      $have{$p}{$act}++;  $nh++; 
      # my $genefile= $genefiles{$p};  #dont need here# $hfile{$p}= $genefile;
      push( @{$hgene{$p}}, $gid);
      } 
    else {
      my($p)= $gid =~ m/^([^\d\W]+)/; $p||=$gid; 
      $miss{$p}{$act}++; $nm++; 
    }  
  } 
  print "# make_mustgenes from $mustfile, have_infile=$nh, miss_in=$nm\n" if $debug; 


    # * FIXME: scan existing $mustout for %gact first, keep those, write only updates
    # some cannot be found in genefiles, write that to mustout ?
      # make_mustgenes to genes/mustgenes.gff, have=673, miss=112
      # make_mustgenes found=644, not found=29
   # if( date($mustfile) < date($mustout) ) ...
   
  my $have_mustout=0;
  my $need_update= 1;
  my $otmp= "$mustout.tmp"; # need to rewrite mustout or not? sometimes not
  
  if( -f $mustout and open(F,$mustout)) {
    while(<F>) {
      next unless(/^\w/);
      my @v= split"\t";
      if( $v[2] =~ m/$mrnatypes/ and /ID=([^;\s]+)/) { 
        my $gid=$1; 
        my $mustact= $v[1];
        my ($iact)= m/;must=(\w+)/;
        
        if($mustact !~ /^must/ and $iact > 0) {
          $mustact= ($iact == 77) ? "mustbest"
            : ($iact == 69) ? "mustalt"
            : ($iact == 66) ? "mustdrop"
            : ($iact == 33) ? "mustskiplocus"
            : ($iact == 22) ? "mustnone" : ""; # mustother ??          
          }
        if($mustact =~ /^must(\w+)/) { my $act=$1; 
          if($gact{$gid} eq $act){ $found{$gid}= $act; $have_mustout++; $nout++; $nout{$act}++;  }
          }
        }
      } close(F);       


    #    -M  Script start time minus file modification time, in days.
    #    -A  Same for access time.  -C  Same for inode change time
    my $M_in  = -M $mustfile; # (stat($mustfile))[9];
    my $M_out = -M $mustout;  # (stat($mustout))[9];
    
    $need_update = 0 if($have_mustout / $nh > 0.85 and $M_in > $M_out);
    $nout= $have_mustout; $M_in= int($M_in); $M_out= int($M_out);
    print "# make_mustgenes in $mustout, have_in=$nh, have_out=$have_mustout, need_update=$need_update, age_in=$M_in, age_out=$M_out\n" if $debug; 
  }
  
  
  if($need_update) {
    $nout=0; %nout=();
    my $altflag = $public_options{'altflag'} || "alttr";
    my $splitflag = $public_options{'splitflag'} || "Split";
    
    my $mustmiss = $nh - $have_mustout;
    print "# make_mustgenes update $mustout, have_infile=$nh, have_mustout=$have_mustout, need=$mustmiss\n" if $debug; 
    
    rename($mustout, "$mustout.old") if(-f $mustout);
    open(MUSTOUT,">$mustout") or die "writing $mustout";
    
    foreach my $p (sort keys %have) {
      my $genefile= $genefiles{$p};
      my @gid= @{$hgene{$p}};
      my %haveg= map{$_,1} @gid;
      unless( -f $genefile) { warn "# Missing must genes.gff=$genefile\n"; next; }
      my $cmd= ($genefile =~ /\.gz/) ? "gunzip -c " : "cat ";
      $cmd .= "$genefile |"; # grep -F -f $geneids - 
      my $nfrom=0;
      open(F, $cmd);
      while(<F>) {
        next unless(/^\w/);
        my @v= split"\t";
        next if($v[2] =~ /^(gene|(five|three)_prime_utr)/i);
        my $g;        
        if(/\tmRNA/ and /ID=([^;\s]+)/) { $g=$1; }
        elsif(/Parent=([^;\s]+)/) { $g=$1; }
        else { next; }
        if($haveg{$g}) { 
          my $act= $gact{$g}; # code as int for  best must=77 ; alt must=69 ; drop must=??          
          my $ord= $gord{$g};  #how to use this??
          my $iact=($act =~ /best/) ? 77 
            : ($act =~ /alt/) ? 69 
            : ($act =~ /drop/) ? 66 
            : ($act =~ /skip/) ? 33 
            : ($act =~ /none/) ? 22 : 11;

 # FIXME: alts need combine w/ primary gene model in genes.gff
 #      : only way seems CDS exon-align test
          my $alt=($act =~ /alt/)?";$altflag=1":"";

          s/\t\S+/\tmust$act/; # put action in gff.source ? or ;act=$act
          s/$/;must=$iact$alt;mustord=$ord/;
          print MUSTOUT $_;
          if(/\tmRNA/) { 
            $found{$g}= $act;  # save and return %gact_have{$g}= $act; for next step
            $nfrom++; $nout++; $nout{$act}++; 
            }
          }
      } close(F);
      
      my $nlook= @gid;
      print "# make_mustgenes from $genefile, want=$nlook, found=$nfrom\n" if $debug; 
      my @gmiss= map { $_.":".substr($gact{$_},0,1); } grep { ! $found{$_} } @gid;
      print "#     missed in $p: @gmiss\n" if(@gmiss); # if $nlook > $nfrom
      # ? write to mustout so we don't keep looking?
    }
    close(MUSTOUT);
  }
  
  my $havenot= $nh - $nout;
  my $nouta = join ",", map{ "$_:$nout{$_}" } sort keys %nout;
  print "# make_mustgenes found=$nout, not found=$havenot, as $nouta\n" if $debug; 
  
#       # debug table
#     if( 0 and $debug ) {    
#       my @ac=sort keys %act;  
#       print "# Must have=$nh; miss=$nm\n"; 
#       print join("\t","Type",@ac),"\n"; 
#       print "# Have ID types\n"; foreach my $p (sort keys %have) { 
#         printf "%10s", $p; foreach my $ac (@ac) { my $v=$have{$p}{$ac}||0; print "\t$v"; } 
#         print "\n"; } 
#       print "\n# Miss ID types\n"; foreach my $p (sort keys %miss) { 
#         printf "%10s", $p; foreach my $ac (@ac) { my $v=$miss{$p}{$ac}||0; print "\t$v"; } 
#         print "\n"; } 
#     }

  return ($mustout, \%found);
}


# $evigene/scripts/bestgenes_update.pl  -vers EG2APm -conf genes/evigene_aphid2.conf -in genes/mustgenes.gff 
# -out genes/mustgenes_pub.gff -act=sort,public -debug
#
# genesort_update: genes/mustgenes.gff TO genes/mustgenes.gff.so
# gene count for EG2APm: ngene=542, nmrna=644, nexon=15211, nother=0
#
# source count for EG2APm:
# mustaltmodel  99
# mustbest      340
# mustdrop      201
# mustskiplocus 4
# rename genes/mustgenes.gff.so.pb genes/mustgenes_pub.gff
# unlink genes/mustgenes.gff.so

sub read_mustgenes
{
  my($mustgenes)= @_;
  # read mustgenes 1st,save all in mem (should be small)
  
  my( $predid, $lchr, $nmrna, $nexon, $nother);
  my( %featlist, %found, %nout, $nout);
  
  open(F, $mustgenes) or die "reading $mustgenes";
  while(<F>) {
    next unless(/^\w/);
    
    my @col= split"\t";  
    my $rloc= \@col; 
    
    if( $col[2] eq "gene" ) { 
      # save it?    
      
    } elsif(  $col[2] =~ m/^($mrnatypes)/ ) {
      ($predid)= m/\bID=([^;\s]+)/;  # die if not found?
      $lchr= $col[0];
      
      my $mustact= $col[1];
      my ($iact)= m/;must=(\w+)/;
      
      if($mustact !~ /^must/ and $iact > 0) {
        $mustact= ($iact == 77) ? "mustbest"
          : ($iact == 69) ? "mustalt"   # FIXME: alts need combine w/ main parents !
          : ($iact == 66) ? "mustdrop"
          : ($iact == 33) ? "mustskiplocus"
          : ($iact == 22) ? "mustnone" : ""; # mustother ??          
        }
        
      if($mustact =~ /^must(\w+)/) { my $act=$1; 
        $found{$predid}= $act;  $nout++; $nout{$act}++; 
      } else {
        $predid="missing$iact:$predid"; next;
      }
      
      $nmrna++; 
      # push( @{$geneidmap{$gid}}, $pubid); # pubid?
      
      # these should be only keepers? act == best|alt ?
      $featlist{$predid}= [$rloc]; #new mrna
      # $commlist{$predid}.= $comsave if($comsave); $comsave=""; # really want to push preceding comms to following pubid
      
    } elsif( $col[2] =~ m/^($exontypes)/ ) {
      my ($pid)= m/\bParent=([^;\s]+)/;
      ## my $ppid= $pid; ## $pubidmap{$pid} || $pubid;
      
      $rloc->[8]=~s/\bID=([^;\s]+);?//; # clean up here or before
      
      if($pid ne $predid and not $featlist{$pid}) {
        next;  #? dont care
        # warn "# out of order gene rec: $pid.$col[2] not in mRNA:$predid\n"; # and do what? any of these?
      }

      $nexon++;
      push( @{$featlist{$pid}},$rloc);
      
    } else { # ?? save in gene rec or not?
      $nother++;
      if(@col > 6 and $col[0] eq $lchr and $featlist{$predid}) {
        push( @{$featlist{$predid}},$rloc);
      }
    }
      
    
  } close(F);
   
  return( \%found, \%featlist); # both hash by predid   
}


# $generec= mustalt_check( $actid, $generec, \%maingenes); 
sub mustalt_check
{
  my($aid, $agene, $maingenes)= @_;
  # need overlap tests 1. mRNA, 2. CDS
  my ($achr,$atb,$ate) = @{$agene->[0]}[0,3,4]; # mrna(b,e)
  my @acds= map{ ($_->[3],$_->[4]) } grep { $_->[2] eq "CDS" } @$agene;
  my($neq, $ideq, $geneeq)= (0,0,0);
  
  foreach my $mid (sort keys %{$maingenes->{$achr}}) #fixme: speed up loop using achr subset
  {
    my $mgene= $maingenes->{$achr}{$mid};
    my ($mchr,$mb,$me) = @{$mgene->[0]}[0,3,4]; # mrna(b,e)
    next unless($mchr eq $achr and $mb < $ate and $me > $atb);
    
    my $iseq=0;
    foreach my $x (@$mgene) { 
      my($xt,$xb,$xe)= @{$x}[2,3,4]; ## ($x->[2],$x->[3],$x->[4]);
      next unless($xt eq "CDS");
      for(my $i=0; $i<@acds; $i+=2) {
        if($xb == $acds[$i] and $xe == $acds[$i+1]) { $iseq++; $ideq=$mid; $geneeq=$mgene; last; }
        }
      last if($iseq > 0); #? only need 1 cds match?
    }
  
    if($iseq > 0) {
      # check if this mgene is alttr, count alttr/gene , set altnum
      $neq++;
      ## last; # wait for all alttr
     } elsif($neq>0) { 
      last; 
     }
  }
  
  if($neq > 0 and $ideq) {
    # change $agene ID to be alt for $geneeq; ** NEED geneeq to record num alts?
    my $ialt= 1 + $neq;  #? or check ideq =~ m/t(\d+)/; 
    (my $altid= $ideq) =~ s/t\d+$//; 
    $altid .= "t".$ialt;
    my $aat= $agene->[0]->[8];
    $aat =~ s/ID=([^;\s]+)/ID=$altid;oidalt=$1/; my $oid=$1; 
    unless($oid eq $aid) { } # what ?
    $agene->[0]->[8]= $aat;
    # oops, need to change Parent=id for rest of @$agene
    foreach my $x (@$agene) { $x->[8] =~ s/Parent=$oid/Parent=$altid/; }
    $aid= $altid;
  }
  
  return ($aid,$agene); # and new id?
}



sub merge_mustkeepgenes
{
  my($ingff, $outgff)= @_;

  ## mustkeepdrop  genes/mustkeepdrop.list         # for bestgenes, expert selections
  # my $mustfile  = $evidence{"mustkeepdrop"}; # change from processed ids to pickgenes.results

  my $mustin   = $evidence{"mustkeeptab"}; # was {"mustkeepdrop"}; # change from processed ids to pickgenes.results
  my $mustout  = $evidence{"mustgff"} || "mustgenes.gff";   
  my $M_in  = -M $mustin; # (stat($mustin))[9];
  my $M_out = -M $mustout;  # (stat($mustout))[9];
    
  my $need_update = ($M_in > $M_out) ? 0 : 1; # and $have_mustout / $nh > 0.85 
  my ($mustgenes, $mustidref)=($mustout, undef);
  
  if($need_update) {
    # die "CALL: ($mustgenes, $mustidref)= make_mustgenes() " if $debug; # TESTING
   ($mustgenes, $mustidref)= make_mustgenes(); # opt: no-file-update ?
   
  } else {
    $mustgenes= $mustout;
    print "# merge_mustkeepgenes using $mustout, age_in=$M_in, age_out=$M_out\n" if $debug; 
  }
  
  # now update ingff with mustgenes to outgff
  # .. need overlap tests for in/must ??  check also for already-have must
  # .. assume ingff unsorted, append must-keep to end  
  # .. assume ingff == public output;  redo that part for mustgenes updates
  # .. added best genes: write sep file, run thru -act=protein,sort,public to update; need -pubid0= option
  
  #? output here = mustgenes.update.gff, 
  #  only those not in ingff, 
  #  include drop/skiplocus as comment? not gene rec
  # THEN run must.update.gff thru this -act protein,sort,public -pubidnum=99999 ... for append
  #     but also process drops

  my($mustidref, $mustfeats)= read_mustgenes($mustgenes);
  #?? foundref == mustidref == hash{predid} = action (best/alt/drop/skiplocus)
  
  #.....
  my( %didup, %haveup, $nup, $nhave, $nnochange, $naltup) = (undef,0,0,0,0);
  my ($lastid,$lastidnum)=("",0);
  my %maingenes;
  
  $outgff="$ingff.mustup" unless($outgff);
  my($inh, $outh)= openio($ingff, $outgff);
  while(<$inh>) { 
    next unless(/^\w/);
    my $inline= $_;
    my @col= split"\t";
    my $rloc= \@col;
    my $chr= $col[0];
    
    if( $col[2] =~ m/^($mrnatypes)/ ) {    
      my($gid)= m/\bID=([^;\s]+)/;
      my($oid)= m/\boid=([^;\s]+)/; # $oid||= "nada";

      # acyp2eg0000001t1 << problem with acyp[2] and t1
      # my $pubid_format = $public_options{'publicid'} || "CG0000000";
      # my $nd= ( $pubid_format =~ s/(0\d+)// ) ? length($1) : 6;
      my ($idnum)= $gid =~ m/0+([1-9]\d+)/;
      if($idnum > $lastidnum) { $lastidnum=$idnum ; $lastid= $gid; }
      
      my $actid= $gid;
      my $act= $mustidref->{$actid};
      if( !$act and $oid) { $actid=$oid;  $act= $mustidref->{$actid}; }
      unless( $act) {
        # write current in gene? no, skip for now, only updates
        $nnochange++;
        
      } elsif($act =~ /drop|skiplocus|none/) {
        # write drop comment ??
        $inline =~ s/;protein=.*$//;
        print $outh "#update\t$act\t$actid\t$inline\n";
        $nup++; $didup{$actid}= $act;
        
      } elsif($act =~ /best|^alt/) {
        # FIXME: new alts to current mRNA : need to match gene ID w/ exon overlap 
        
        # best here means input already has update : dont rewrite
        # alt may be updone, but need check?
        
        $inline =~ s/;protein=.*$//;
        print $outh "#updone\t$act\t$actid\t$inline\n"; #?? or not
        $nhave++; $haveup{$actid}= $act;
      }
 
      $maingenes{$chr}{$gid}= [$rloc]; #new mrna
     
    } elsif( @col > 6 ) {
      ### skip all other input feats ?
      my($pid)= m/\bParent=([^;\s]+)/;
      push( @{$maingenes{$chr}{$pid}},$rloc) if($col[2] eq "CDS" and $maingenes{$chr}{$pid}); #? need only this type
    }
  }
  
  # now pull mustidref not handled
  foreach my $actid (sort keys %$mustidref) {
    next if ($didup{$actid} or $haveup{$actid} or not $actid);
    my $act= $mustidref->{$actid} or next; # bug here?
    if($act =~ /drop/) { $nhave++; next; } #? already done?    
    print $outh "#update\t$act\t$actid\n";
    unless($act =~ /skiplocus|none/) { 
      # FIXME: new alts to current mRNA : need to match gene ID w/ exon overlap 
      my $generec= $mustfeats->{$actid} or do { warn "# Missing must $actid gene\n"; next; };
      
      if($act =~ /alt/) { 
        my $altid=$actid; 
        #** this will be slow: iterates all geneid in maingenes: hash by chr/ref?
        ($altid,$generec)= mustalt_check( $actid, $generec, \%maingenes); 
        if($altid ne $actid) { 
          my $achr= $generec->[0]->[0];
          $maingenes{$achr}{$altid}= $generec;
          $naltup++;
         } # what? add to maingenes so next alt will pick it up
        }
      
      my @rows= @$generec;
      foreach my $r (@rows) { print $outh join "\t", @$r; } # has \n 
    }
    $nup++; $didup{$actid}= $act;
  }
  
  print $outh "\n#lastid\t$lastidnum\t$lastid\n"; # can this be used for public update?  

  my @gupd= map { $_.":".substr($didup{$_},0,1); } sort keys %didup;
  print $outh "\n# merge_mustkeepgenes for $MSRC: changed=$nup, haveup=$nhave, nochange=$nnochange, naltu=$naltup\n# merge_updates: @gupd\n";
  print       "\n# merge_mustkeepgenes for $MSRC: changed=$nup, haveup=$nhave, nochange=$nnochange, naltu=$naltup\n# merge_updates: @gupd\n"
    if $debug;
  close($outh);  

#.. final merge; print this, dont run ?
# $ingff is public source
# this method creates, uses mustgenes.gff from must.table
# $outgff here is ONLY mustupdate, needs -act protein,sort,public before final merge
  
  my $NOTE = <<"EON";  
# final merge: run this AFTER -act protein,sort,public  of $outgff
# * FIXME: need step to merge mustalt with primary gene models: overlap-cds test?, reassign alt ids
# * FIXME: dropids - misses gene row if no alt kept

egrep '^#update.(drop|skiplocus)' $outgff | perl -ne \\
 '(\$ac)=(split)[1]; if(m/\\tID=([^;\s]+)/) { print "\$ac\\t\$1\\n";}' > $outgff.dropids

cat $outgff.dropids $ingff | perl -ne \\
'if(/^(drop|skiplocus)\\s+(\\S+)/) { \$drop{\$2}=1; } else { if(/^\\w/){ 
if(/\\tgene/) { (\$g)=m/ID=([^;\\s]+)/; \$d=\$g."t1"; \$drop=(\$drop{\$d})?1:0; } 
elsif(/\\tmRNA/) { (\$d)=m/ID=([^;\\s]+)/; \$drop=(\$drop{\$d})?1:0; } 
elsif( (\$d)=m/Parent=([^;\\s]+)/ ) { \$drop=(\$drop{\$d})?1:0; } 
s/^/#drop./ if(\$drop); } print; }' > $ingff.upd1

# .. add in .mustchimer, other .must to end of .upd2
egrep -v '^#|^$' $outgff | cat $ingff.upd1 - | grep -v '^#drop\.' > $ingff.upd2
# end final merge

EON
  print $NOTE;

  return $outgff;
}

=item mustchimera work

## add in separate from mustgenes.gff (??)
cp an8valid/acypi-chimera.egover.top.gff bestgenes.DGILpub8d4.gff.mustchimer

perl -pi -e'BEGIN{$i=37330;} s/\tacypi1split2/\tEG2APu/; if(/\tmRNA/)\
{ $eg="acyp2eg00".($i++)."t1"; s/ID=([^;\s]+)/ID=$eg;oid=$1/; $od=$1; \
s/^/#update\tbest\t$od\tchimera\n/; s/$/;must=77;osrc=acypi1split2/; } \
else { s/Parent=/Parent=$eg;oid=/; s/;Target=[^;\n]+//; } ' \
bestgenes.DGILpub8d4.gff.mustchimer


cat bestgenes.DGILpub8d4.gff.mustup2.dropids  bestgenes.DGILpub8d3.gff | perl -ne \
'if(/^(drop|skiplocus)\s+(\S+)/) { $drop{$2}=1; } else { if(/^\w/){ \
if(/\tgene/) { ($g)=m/ID=([^;\s]+)/; $d=$g."t1"; $drop=($drop{$d})?1:0; } \
elsif(/\tmRNA/) { ($d)=m/ID=([^;\s]+)/; $drop=($drop{$d})?1:0; } \
elsif( ($d)=m/Parent=([^;\s]+)/ ) { $drop=($drop{$d})?1:0; } \
s/^/#drop./ if($drop); } print; }' > bestgenes.DGILpub8d4.gff.upd1


cat bestgenes.DGILpub8d4.gff.mustup2 bestgenes.DGILpub8d4.gff.mustchimer | \
egrep -v '^#|^$' | cat bestgenes.DGILpub8d4.gff.upd1 - | \
grep -v '^#drop.' > bestgenes.DGILpub8d4.gff.upd2

#.......

# redo w/ new mustgenes 2011.jun.05

$evigene/scripts/bestgenes_update.pl -conf genes/evigene_aphid2.conf \
-vers EG2APu -act=mustkeep -debug \
-in genes/bestgenes.DGILpub8d3.gff \
-out genes/bestgenes.DGILpub8d5.gff \
> & genes/mustgenes5.log

# make_mustgenes found=898, not found=33, as altmodel:138,best:408,drop:348,skiplocus:4
#lastid 36386   acyp2eg0036386t1
#last.upd2      acyp2eg0037265t1 ; chimera start 37330

$evigene/scripts/bestgenes_update.pl -conf genes/evigene_aphid2.conf \
-vers EG2APu -act=protein,sort,public -pubid=37000 -debug \
-in genes/bestgenes.DGILpub8d5.gff.mustup \
-out genes/bestgenes.DGILpub8d5.gff.mustup2 \
> & genes/mustgenes5p.log

# add_missprotein for EG2APu: naa=364, nadd=148, nfull=284, npart=80, aasize=371601
# gene count for EG2APu: ngene=251, nmrna=365, nexon=10392, nother=0
# source count for EG2APu: mustaltmodel  114  mustbest      251

# final merge: run this AFTER -act protein,sort,public  of genes/bestgenes.DGILpub8d3.gff.mustup

#   egrep '^#update.(drop|skiplocus)' genes/bestgenes.DGILpub8d3.gff.mustup | perl -ne \
#    '($ac)=(split)[1]; if(m/\tID=([^;s]+)/) { print "$ac\t$1\n";}' > genes/bestgenes.DGILpub8d3.gff.mustup.dropids
#   
#   cat genes/bestgenes.DGILpub8d3.gff.mustup.dropids genes/bestgenes.DGILpub8d3.gff | perl -ne \
#   'if(/^(drop|skiplocus)\s+(\S+)/) { $drop{$2}=1; } else { if(/^\w/){ 
#   if(/\tgene/) { ($g)=m/ID=([^;\s]+)/; $d=$g."t1"; $drop=($drop{$d})?1:0; } 
#   elsif(/\tmRNA/) { ($d)=m/ID=([^;\s]+)/; $drop=($drop{$d})?1:0; } 
#   elsif( ($d)=m/Parent=([^;\s]+)/ ) { $drop=($drop{$d})?1:0; } 
#   s/^/#drop./ if($drop); } print; }' > genes/bestgenes.DGILpub8d3.gff.upd1
#   
#   # .. add in .mustchimer, other .must to end of .upd2
#   egrep -v '^#|^ genes/bestgenes.DGILpub8d3.gff.mustup | cat genes/bestgenes.DGILpub8d3.gff.upd1 - | grep -v '^#drop.' 
#   > genes/bestgenes.DGILpub8d3.gff.upd2



# ?? new idnums?  bestgenes.DGILpub8d4.gff.mustchimer at acyp2eg0037330t1 ; from  an8valid/acypi-chimera.egover.top.gff

cat bestgenes.DGILpub8d5.gff.mustup2 bestgenes.DGILpub8d4.gff.mustchimer | egrep -v '^#|^$' |\
cat bestgenes.DGILpub8d5.gff.upd1 - | grep -v '^#drop.' > bestgenes.DGILpub8d5.gff.upd2

#...........

=item mustkeep work

$evigene/scripts/bestgenes_update.pl -conf genes/evigene_aphid2.conf -in genes/bestgenes.DGILpub8d3.gff
 -out genes/bestgenes.pub8d3_upd1.gff -vers EG2APm -act=mustkeep -debug

# merge_mustkeepgenes using genes/mustgenes.gff, age_in=0.775081018518519, age_out=0.773344907407407
# merge_mustkeepgenes for EG2APm: changed=372, lastup=272, nochange=40966, updates: 
# :0 0:0 ACYPI33128-RA:b ACYPI38045-RA:b ACYPI53287-RA:b
# ACYPI53893-RA:b AUGepi4p1s10g138t1:a AUGepi4p1s10g28t1:b
# AUGepi4p1s1g6t1:b AUGepi4p1s1g7t1:b AUGepi4p1s1g89t1:d
# AUGepi4p1s1g90t1:b AUGepi4p1s2g65t1:b AUGepi4p1s7g65t1:a
# AUGepi4p2s10g9t1:s AUGepi4p2
#lastid 36386   acyp2eg0036386t1


#.......

$evigene/scripts/bestgenes_update.pl -conf genes/evigene_aphid2.conf \
 -in genes/bestgenes.DGILpub8d3.gff.mustup -out genes/bestgenes_pub8d3must1.gff \
 -vers EG2APm -act=protein,sort,public -pubid=37000 -debug

# # add_missprotein: scripts/pa2dgg_gff3proteins.pl genes/bestgenes.DGILpub8d3.gff.mustup.mp29747.tmp genome/aphid2asm.fa prot |
# 
# # add_missprotein for EG2APm: naa=253, nadd=98, nfull=194, npart=59, aasize=259786
# # genesort_update: genes/bestgenes.DGILpub8d3.gff.mustup.mp TO genes/bestgenes.DGILpub8d3.gff.mustup.mp.so
# 
# # gene count for EG2APm: ngene=176, nmrna=253, nexon=7115, nother=0
# 
# # source count for EG2APm:
# # mustaltmodel  77
# # mustbest      176
# # rename genes/bestgenes.DGILpub8d3.gff.mustup.mp.so.pb genes/bestgenes_pub8d3must1.gff
# unlink genes/bestgenes.DGILpub8d3.gff.mustup.mp
# unlink genes/bestgenes.DGILpub8d3.gff.mustup.mp.so

# .. final merge of 
#  genes/bestgenes.DGILpub8d3.gff and genes/bestgenes_pub8d3must1.gff = genes/bestgenes.DGILpub8d3.gff.mustup.mp.so.pb 

egrep '^#update.(drop|skiplocus)' genes/bestgenes.DGILpub8d3.gff.mustup | perl -ne\
 '($ac)=(split)[1]; if(m/\tID=([^;\s]+)/) { print "$ac\t$1\n";}' > genes/bestgenes.DGILpub8d3.gff.mustup.dropids
 
cat genes/bestgenes.DGILpub8d3.gff.mustup.dropids genes/bestgenes.DGILpub8d3.gff | perl -ne\
'if(/^(drop|skiplocus)\t(\S+)/) { $drop{$2}=1; } else { if(/^\w/){ \
if(/\tgene/) { ($g)=m/ID=([^;\s]+)/; $d=$g."t1"; $drop=($drop{$d})?1:0; } \
elsif(/\tmRNA/) { ($d)=m/ID=([^;\s]+)/; $drop=($drop{$d})?1:0; } \
elsif( ($d)=m/Parent=([^;\s]+)/ ) { $drop=($drop{$d})?1:0; } \
s/^/#drop./ if($drop); } print; }' \
> genes/bestgenes.DGILpub8d3.gff.upd1

egrep -v '^#|^$' genes/bestgenes.DGILpub8d3.gff.mustup.mp.so.pb | \
cat genes/bestgenes.DGILpub8d3.gff.upd1 - | grep -v '^#drop\.' > genes/bestgenes.DGILpub8d3.gff.upd2

#.. final merge; print this, dont run ?

mustout=genes/bestgenes.DGILpub8d5.gff.mustup2
mustchimer=genes/bestgenes.DGILpub8d4.gff.mustchimer
ingff=genes/bestgenes.DGILpub8d3.gff 
outgff=genes/bestgenes.DGILpub8d5.gff

  egrep '^#update.(drop|skiplocus)' $mustout | perl -ne\
    '($ac)=(split)[1]; if(m/\tID=([^;\s]+)/) { print "$ac\t$1\n";}' > $mustout.dropids

  cat $mustout.dropids $ingff | perl -ne\
  'if(/^(drop|skiplocus)\t(\S+)/) { $drop{$2}=1; } else { if(/^\w/){ \
  if(/\tgene/) { ($g)=m/ID=([^;\s]+)/; $d=$g."t1"; $drop=($drop{$d})?1:0; } \
  elsif(/\tmRNA/) { ($d)=m/ID=([^;\s]+)/; $drop=($drop{$d})?1:0; } \
  elsif( ($d)=m/Parent=([^;\s]+)/ ) { $drop=($drop{$d})?1:0; } \
  s/^/#drop./ if($drop); } print; }' \
  > $outgff.upd1
  
  cat $mustout $mustchimer | egrep -v '^#|^$' | cat $outgff.upd1 - | grep -v '^#drop\.' > $outgff.upd2
    
#..

 grep -c mRNA  genes/bestgenes.DGILpub8d3.gff.upd?
genes/bestgenes.DGILpub8d3.gff:41268
genes/bestgenes.DGILpub8d3.gff.upd1:41268 # comments, same
genes/bestgenes.DGILpub8d3.gff.upd2:41406

genes/bestgenes.DGILpub8d4.gff.upd2:41587
genes/bestgenes.DGILpub8d5.gff.upd2:41559  # must combined several fragment genes

 grep -c '      gene' genes/bestgenes.DGILpub8d3.gff.upd?
genes/bestgenes.DGILpub8d3.gff:36386
genes/bestgenes.DGILpub8d3.gff.upd1:36386
genes/bestgenes.DGILpub8d3.gff.upd2:36531  # more .. alt-tr added w/ sep ids

genes/bestgenes.DGILpub8d4.gff.upd2:36666  #  bug
genes/bestgenes.DGILpub8d5.gff.upd2:36490  # fixed bug alt-tr>newgene


=cut

#........................................

#* add missing protein=  n=1655, for some source = pasa2_aphid3, acypi, ars17trinity, ..
#  1487 of 1655 have prot= in pasaupdate.gff, as upstatus=original; skipped here
#  most of rest (450) are in rnagene/*{pasa,cuf8,trinity}*.aa(.gz); some in refseq/acypiprot.fa.gz
# .. use $evigene/scripts/pa2dgg_gff3proteins.pl to add all missing prots ? are they same as source?
# .. pa2dgg_gff3proteins.pl  genes.gff genome/$dgenome.fa addprot
#
#  added: in add_missprotein; check protein and set status in aalen=aasize,aapct,aastatus
#   pasa status=complete,5prime_partial,35prime_partial,internal 

# .. use $evigene/scripts/pa2dgg_gff3proteins.pl to add all missing prots ? are they same as source?
# only update mRNA missing protein= field; 
#  1.  grep -v protein= from mRNA, but also same ID of exon + CDS
#  2.  pa2dgg_gff3proteins.pl  missprotgenes.gff genome/$dgenome.fa addprot
#  3.  grep mRNA w/ protein= from result, add to this mRNA 
# OR modify pa2dgg_gff3proteins to simplify this?
# .. all pasa gff utils seem to use index_GFF3_gene_objs(gfffile), no way to embed this here
# .. these need PASA/PerlLib/ in perl path

## wrong perl path:
##  Can't locate Gene_obj.pm in @INC (@INC contains: /export/udisk2/work/aphid/scripts/../PerlLib 

sub add_missprotein
{
  my($ingff, $outgff)= @_;
  my $missgff="$ingff.mp$$.tmp"; # temp, delete here

  # do this pro=  mark before or here? use evigene.conf anoption settings
  my $evcmd   = $public_options{'addprot'}; # this is useless, build cmd here
  my $genome  = $evidence{'genome'}; #  need this

  #my $progname='addproteins'; 
  #my $progpath= $programs{$progname} || $progname;
  # $evcmd= "$progpath $missgff $genome prot |"; 
  # $evcmd =~ s/$progname/$progpath/;  
  # $evcmd= ($evcmd) ? "$evcmd $missgff $genome prot |" : ""; 
  # my $ok= `which $progpath`; chomp($ok);

use constant NEW_ADDPROT => 1;
#  addprot addproteins  -full=1 -dna $genome  -genes $genes  #  > $newgenes  
  $evcmd= setprogram( $evcmd, "addproteins");
  
  my $ok=1;
  if($ok) { $ok= (-f $genome) ? 1 : 0; }
  unless( $ok ) { warn "# ERR in command: $evcmd\n"; return }
  
  my( $miss, $id, $aa, $nupd, $nfull, $npart, $naa, $sumalen, ) = (0) x 10;
  my (%aa);
  
if(NEW_ADDPROT) {
  $outgff="$ingff.mp" unless($outgff);

  $evcmd =~ s/\$genome/$genome/; #?? ok
  $evcmd =~ s/\$genes/$ingff/;
  $evcmd .= " > $outgff";
  # sub docommand { verbose(@_); my $err= ($norun) ? 0 : system(@_); return $err; }
  my $err= system($evcmd);

} else { # pasa gff2proteins ... old
  
  my($inh, $outh)= openio($ingff, $missgff);
  while(<$inh>) {
    next unless(/^\w/);
    if(/\tmRNA\t/) {  $miss= (m/;protein=\w+/)? 0 : 1; } 
    print $outh $_ if($miss);
  }
  close($inh); close($outh);

  ## $evcmd .= '|'; # ok??
  $evcmd= ($evcmd) ? "$evcmd $missgff $genome prot |" : ""; 
  print "# add_missprotein: $evcmd\n" if($debug);
  open(P,$evcmd) or die "ERR: $evcmd";
  ($id,$aa)= ("","");
  while(<P>) {
    # NOT NOW: record here is >fastaID AAAA ..."
    if(/^>(\S+)/){ $aa{$id}=$aa if($id and $aa); $id=$1; $aa=""; }
    elsif(/^\w/) { chomp; $aa.=$_; }
  } close(P);
  $aa{$id}=$aa if($id and $aa); 
  unlink($missgff);
  
  $outgff="$ingff.mp" unless($outgff);
  my($inh, $outh)= openio($ingff, $outgff);
  while(<$inh>) {
    if(/^\w/ and /\tmRNA\t/) { 
      my($id)=m/ID=([^;\s]+)/; 
      if( my $anew= delete $aa{$id} ) {    
        my $lnew= length($anew)||0; $lnew-- if($anew =~ m/\*$/); 
        my($clen,$xlen)= m/cxlen=(\d+).(\d+)/; # have all this?
        my $ap= ($xlen>0) ? int(0.5 + 300 * $lnew / $xlen) : 0;

        my ($aold)= m/protein=([^;\n]+)/;  # should be missing
        if($aold) { } #what?
        
        s/;?(aalen|protein)=[^;\n]+//g; 
        my $at="aalen=$lnew,$ap%;protein=$anew";
        s/$/;$at/;  $nupd++;
      }

   #* add: in add_missprotein or public_update, check protein= and set 
   #   pasa status=complete,5prime_partial,3prime_partial,internal 
      my ($aa)= m/protein=([^;\n]+)/;  # should be missing
      if($aa) { 
        my $al= length($aa)||0; $al-- if($aa =~ m/\*$/); 
        my($al0,$ap)= m/aalen=(\d+).(\d+)/; $ap ||=0; # some have this, not cxlen=
        my($clen,$xlen)= m/cxlen=(\d+).(\d+)/; # have all this?
        $ap= ($xlen>0) ? int(0.5 + 300 * $al / $xlen) : $ap;
        my $astat=0;
        if($aa =~ /^M/) { $astat |= 1; }
        if($aa =~ /\*$/) { $astat |= 2;} # ** AUG proteins lack '*' unless pasa-updated
        elsif($id =~ /AUG/) { $astat |= 2; } #  also check for inner X == augustus-fake for stop ?
      
        my $istop= index($aa,'*');
        if($istop < 1 and $id =~ /AUG/) { $istop= index($aa,'X'); }
        if($istop > 0 and $istop < $al-1) { $astat="partialinner"; }
        elsif($astat == 3) { $astat="complete"; }
        elsif($astat == 1) { $astat="partial3prime"; }
        elsif($astat == 2) { $astat="partial5prime"; }
        else { $astat="partial"; }
        if($astat =~ /complete/) { $nfull++; } else { $npart++; }
        $sumalen+= $al; $naa++;
        unless( s/;aalen=[^;\s]+/;aalen=$al,$ap%,$astat/) {
          s/$/;aalen=$al,$ap%,$astat/;
        }
        s/;status=(complete|3prime_partial|5prime_partial|internal)//; # old pasa prot status
      }  
      
    }
    print $outh  $_; # OUT; to where ?? should be pipe-able
  } close($inh);

  print $outh "\n# add_missprotein for $MSRC: naa=$naa, nadd=$nupd, nfull=$nfull, npart=$npart, aasize=$sumalen\n"; # outh ?
  print       "\n# add_missprotein for $MSRC: naa=$naa, nadd=$nupd, nfull=$nfull, npart=$npart, aasize=$sumalen\n"
    if($debug); 
  
  close($outh);
} # old pasa addprot

  return $outgff;
}


sub nontriv_scores
{
  my($scoretype, $ev, $ew)= @_;
  my @scores= split",",$scoretype;
  foreach my $i (0 .. $#scores) {
    my($k,$v)=split ":",$scores[$i];
    unless($k =~ /CDS|UTR|terepeat|paralog/) {  # FIXME
      push @$ev, $i; 
      my $w=($k=~/nintron/)?10:1; 
      push @$ew, $w; 
      }
  }
}

sub remove_trivial 
{
  my($ingff, $outgff)= @_;
  my %kid;
  my $keeps    = $public_options{"keepsource"} || "nzAtAll";
  my $keepids  = $public_options{"keepids"};
  $keeps =~ s/[,\s]+/|/g;
  if($keepids and open(F,$keepids)) { while(<F>){ chomp; $kid{$_}=1;} close(F); }

  my $MINSCORE= 5; # FIXME: opt, use instead pubopt:trivial.dropscore=XXX:n,YYY:n ??
  my @ev=(); ## (0,2,3,4,5,6,7); ## FIXME: get from in.gff or 
  my @ew=(); ## (1,1,1,1,1,1,10);
#scoretype: homolog:9,paralog:1,ref:2,est:3,pro:3,rseq:2,intr:20,nintron:40,inqual:20,maqual:5,terepeat:-3,UTR:3,CDS:1
#dropscore: homolog:30,ref:119,est:119,pro:119,rseq:119,intr:1,inqual:10,maqual:99,+CDS:180
  
  my $scoretype = $public_options{"scoretype"} || ""; # dropscore ??
  nontriv_scores($scoretype, \@ev, \@ew) if($scoretype);
  
  my ($skip, %drops);
  $outgff="$ingff.tv" unless($outgff);
  my($inh, $outh)= openio($ingff, $outgff);
  while(<$inh>) {

    if(/^#scoretype:\s*(\S+)/) { # alwasy use inh value? dropscore ??
      $scoretype=$1; @ev= @ew= ();
      nontriv_scores($scoretype, \@ev, \@ew);
    }
    
    if(/^\w/) { 
      if(/\tmRNA/) { $skip=0; 
      my($sr,$pv)=(split"\t")[1,5]; my($g)=m/ID=([^;\s]+)/; 
      unless( ($sr =~ m/$keeps/ or $kid{$g}) and @ev) {  
        my @p= map{s,/\d+,,; $_} split",", $pv; 
        my ($al,$ap)=m/aalen=(\d+).(\d+)/; 
        my($ph)=m/pHOBEST=(\d+)/; 
        
        my $ev=0; for my  $i (@ev) { $ev += $ew[$i] * $p[$i]; } 
        if( $ev < $MINSCORE ) { $skip=($p[1] > 0 and $al > 120 and $ap >= 66 and $ph >= 66)?0:2; } 
        
         if ($skip) { s/$/;skip=$skip/;  $drops{$sr}++; }
        } 
      } 
    
    s/^/#x./ if($skip); # but see below public_update that removes these
    ## option write these to 2nd file? just keep $ingff.tv output
    } 
    
    print $outh $_; # ; to where ?? should be pipe-able
  }
  
  print $outh "\n# trivial drops for $MSRC:\n";
  foreach my $s (sort keys %drops) { print $outh "# $s\t$drops{$s}\n"; }
  if($debug) {
  print "\n# trivial drops for $MSRC:\n";
  foreach my $s (sort keys %drops) { print "# $s\t$drops{$s}\n"; }
  }
  
  close($outh);
  return $outgff;
}


#?? use for mustgene updates, or run that thru pubgene_update as sep. file?
#?  need also to addgene row via genesort_update
##  easiest is to write new update.gff then process that for sort,public
#   append to pubgene.gff

# sub pubgene_rowupdate
# {
#   my($inrow)= @_;
#   local $_= $inrow;
#   
# #  # caller params  ($predgene,$pubgene,$pubid_format)
# #   my ($pubid, $predid, $predgene, $pubgene, $pubidnum)= (0) x 10;
# #   my (%pubidmap, %srccount); # include alt-tr ?
# #  
# #   my $dropannot = $public_options{'dropannot'} || "";
# #   my $dropexonannot = $public_options{'dropexonannot'} || "";
# #   $dropannot =~ s/[,\s]+/|/g;
# #   $dropexonannot =~ s/[,\s]+/|/g;
# # 
# #   my $movescorevec = $public_options{'movescorevec'} || 0;
# #   my $pubid_format = $public_options{'publicid'} || "CG0000000";
# #   my $altid_format = $public_options{'altid'} || "t00";
# #   my $altflag = $public_options{'altflag'} || "alttr";
# #   my $oidtag= "oid"; # sid?
# #   
# #   my $nd= ( $pubid_format =~ s/(0\d+)// ) ? length($1) : 6;
# #   $pubid_format .= '%0'.$nd.'d'; # CG%06d
# # 
# #   $nd= ( $altid_format =~ s/(\d+)// ) ? length($1) : 1;
# #   $ALTKEY= $altid_format; # e.g. 't'
# #   $altid_format .= '%0'.$nd.'d'; # t%02d
# 
#   
#   my @col=split"\t";  
#   s/\t(\S+)\t/\t$MSRC\t/;  # source column
#   my $src=$1;
#   
#   if( $col[2] eq "gene" ) { 
#     # save it* also handle pubid update
#     ($predgene)= m/\bID=([^;\s]+)/;
#     unless($MID) {
#       $pubgene= sprintf( $pubid_format, ++$pubidnum); # always add tNNN?
#       s/\bID=/ID=$pubgene;$oidtag=/;  
#       $pubidmap{$predgene}= $pubgene;      
#       $pubidmap{$pubgene}= $predgene;
#       }
#     if($dropannot) { s/\b($dropannot)=[^;\n]+;?//g ; s/\b($dropannot)=;//g; }
#     
#   } elsif( $col[2] =~ m/^($mrnatypes)/) { 
#     ($predid)= m/\bID=([^;\s]+)/;
#     
#     my($pgid)= m/\bgene=([^;\s]+)/;
#     if($pgid and $predgene) {
#       if( $pgid ne $predgene) { }          # error
#       elsif($pubgene) { s/\bgene=$predgene/gene=$pubgene/; }
#     }
#     
#     $srccount{$src}++;
#     #? change oid=XXX;osrc=SSS to oid=SSS:XXX ? but then do for each feature
#     unless( s/;osrc=[^;\s]+/;osrc=$src/ ) { s/$/;osrc=$src/; }
#     
#     $isalt= (m/\b$altflag=(\d+)/) ? $1 : 0; #? or use predid =~ /t(\d+)$/;
#     $isalt=1 if(m/status=alt-splice/); # upstatus=alt-splice addition # pasa
#     if($altcomm) { $altcomm=0; $isalt=1 unless($isalt); }
#     
#     if($isalt) { 
#       unless( m/;$altflag=$isalt/ or s/;$altflag=[^;\s]+/;$altflag=$isalt/ ) { s/$/;$altflag=$isalt/; }
# 
#       ## pasa upd alt-tr before prime tr; problems here w/ no prime pubid
#       my $gid = $predid;
#       my $altnum= ( $gid =~ s/t(\d+)$// ) ? $1 : $isalt+1; # "666"; #?? not working right
#       
#       my $aid= $pubidmap{$predid} || $pubidmap{$gid} || $pubidmap{$gid."t1"} || $pubidmap{$gid."t2"};
#       # got a handful of dupl alt ids .. mix up b/n pasa update and userchoice isalt > pasa primary
#       if($aid) { $aid =~ s/t\d+$//; } ## fixme below/here
#       
#       unless($aid) { 
#         $aid= ($MID) ? $MID.$predid : ($pubgene) ? $pubgene : sprintf( $pubid_format, ++$pubidnum); 
#       }
#       $aid .= sprintf( $altid_format, $altnum); #? always append tNNN to prime id?
#       
#       $pubid= $aid;
#       s/\bID=/ID=$pubid;$oidtag=/;  
#       
#     } elsif($MID) {
#       $pubid= $MID.$predid;      
#       s/\bID=/ID=$MID/;  
#       
#     } else {
# 
#       ## pasa upd alt-tr before prime tr; problems here w/ no prime pubid
#       ## but problems here w/ non-alt but altid: AUGepir2s234g37t1,AUGepir2s234g37t10 > both EG2ap014820t1
#       ##   AUGepir2s671g15t1,AUGepir2s671g15t15 > EG2ap030523t1
#       ## .. not marked as alttr=1 but from augustus-alt-predict run
#   
#       (my $gid = $predid)  =~ s/t(\d+)$//;
#       $pubid= $pubidmap{$predid} || $pubidmap{$gid} || $pubidmap{$gid."t1"} || $pubidmap{$gid."t2"} || "";
# 
#       if($pubid) {
#         if($pubid =~ /t1$/) { $pubid=""; } else { $pubid =~ s/t(\d+)$/t1/; }
#         $pubid="" if( $pubidmap{$pubid} ); 
#       }
#       
#       if($pubid) { 
#         # nada
#       } elsif($pubgene) {
#         $pubid = $pubgene . sprintf( $altid_format, 1);
#       } else {
#         $pubid= sprintf( $pubid_format.$altid_format, ++$pubidnum, 1); # always add tNNN?
#       }
#       
#       s/\bID=/ID=$pubid;$oidtag=/;  
#     }
#     
#     $pubidmap{$predid}= $pubid;      
#     $pubidmap{$pubid}= $predid;
#     
#     if($dropannot) { s/\b($dropannot)=[^;\n]+;?//g ; s/\b($dropannot)=;//g; }
#     # rename annots here?  pct_support > predsupport=
#     
#   } else {  # assume gene parts are contiguous?
#     my ($pid)= m/\bParent=([^;\s]+)/;
#     my $ppid= $pubidmap{$pid} || $pubid;
#     if($MID) {
#       s/\bParent=/Parent=$MID/;  
#     } else {
#       s/\bParent=/Parent=$ppid;$oidtag=/; 
#     }
#     if($dropexonannot) { s/\b($dropexonannot)=[^;\n]+;?//g ; s/\b($dropexonannot)=;//g; }
#   }
# 
#   # option: move scorevec to ;scorevec= and scoresum= to  col.5
#   if($movescorevec) {
#     my $sv=$col[5];
#     if($sv=~/,/) { 
#       my ($ss) = m/scoresum=(\d+)/;        
#       unless($ss) { map { $ss += $_ } split",",$sv; }
#       s/$col[3]\t$col[4]\t$col[5]/$col[3]\t$col[4]\t$ss/;
#       s/$/;scorevec=$sv/ unless($ss == 0 or $col[2] =~ /exon|CDS/);
#     }
#   }
# 
#   s/;;/;/g; # crap
#   return $_;
# }



## do more here than pubid
  # add opts here? to add gene record enclosing mRNAs; change mRNA to ncRNA if flagged
  # improve alt-tr info
  # maybe: re-sort gene records by location before assigning numid AND opt add gene record
  
sub public_update 
{
  my($ingff, $outgff)= @_;

  my ($pubid, $predid, $predgene, $pubgene, $pubidnum)= (0) x 10;
  my (%pubidmap, %srccount); # include alt-tr ?
  $pubidnum= $pubidnum_start;
  
  ## my $call_ncrna = $public_options{'call_ncrna'} || ""; #harder, maybe not here: see overlapbestgenes

  my $dropannot = $public_options{'dropannot'} || "";
  my $dropexonannot = $public_options{'dropexonannot'} || "";
  $dropannot =~ s/[,\s]+/|/g;
  $dropexonannot =~ s/[,\s]+/|/g;

  my $movescorevec = $public_options{'movescorevec'} || 0;
  
  my $pubid_format = $public_options{'publicid'} || "CG0000000";
  my $altid_format = $public_options{'altid'} || "t00";
  my $altflag = $public_options{'altflag'} || "alttr";
  my $splitflag = $public_options{'splitflag'} || "Split";
  my $oidtag= "oid"; # sid?
  
  my $nd= ( $pubid_format =~ s/(0\d+)// ) ? length($1) : 6;
  my $pubid_prefix= $pubid_format;
  $pubid_format .= '%0'.$nd.'d'; # CG%06d

  # *** FIXME: replace 't\d+' patts thruout w/ global $ALTKEY from $altid_format ***
  $nd= ( $altid_format =~ s/(\d+)// ) ? length($1) : 1;
  $ALTKEY= $altid_format; # e.g. 't'
  $altid_format .= '%0'.$nd.'d'; # t%02d

  
  # x: preserve source predictor col. as mRNA attribute;
  # fixme annots for alt-tr, pasa updates, mustkeep updates

  $outgff="$ingff.pb" unless($outgff);
  my($inh, $outh)= openio($ingff, $outgff);
  my ($altcomm,$isalt,$issplit)=(0) x 9;
  
  while(<$inh>) {
  
    if(/^#x/) { next; } # drop excluded gff here
    elsif(/^#a\./) { $altcomm=1; s/^#a\.//; } # alttr comment ??    
    
      # ?? Make this a sub{} for mustgenes to use ?
    if(/^\w/) {
      
      my @col=split"\t";  
  
      s/\t(\S+)\t/\t$MSRC\t/;  # source column
      my $src=$1;
      
      if( $col[2] eq "gene" ) { 
        # save it* also handle pubid update
        ($predgene)= m/\bID=([^;\s]+)/;
        unless($MID) {
          if($predgene =~ /^$pubid_prefix/) { $pubgene= $predgene; } # already correct?? from mustgenes?
          else {
          $pubgene= sprintf( $pubid_format, ++$pubidnum); # always add tNNN?
          s/\bID=/ID=$pubgene;$oidtag=/;  
          }
          $pubidmap{$predgene}= $pubgene;      
          $pubidmap{$pubgene}= $predgene;
          }
        if($dropannot) { s/\b($dropannot)=[^;\n]+;?//g ; s/\b($dropannot)=;//g; }
        
      } elsif( $col[2] =~ m/^($mrnatypes)/) { 
        ($predid)= m/\bID=([^;\s]+)/;
        
        my($pgid)= m/\bgene=([^;\s]+)/;
        if($pgid and $predgene) {
          if( $pgid ne $predgene) { }          # error or may be pubgene
          elsif($pubgene) { s/\bgene=$predgene/gene=$pubgene/; }
  #*** FIXME  2011sep, have two gene=ID in mRNA:  gene=evigID;gene=oidID; dont want 2nd
  # ... from orig annot? strip out 2nd in dropannot
  # ID=nasv2eg0000089t1;oid=r8nvit1v2Svelbig0Loc3247t1;gene=nasv2eg0000089;gene=nvit1v2Svelbig0Loc3247t1;
        }
        
        $srccount{$src}++;
        #? change oid=XXX;osrc=SSS to oid=SSS:XXX ? but then do for each feature
        unless( s/;osrc=[^;\s]+/;osrc=$src/ ) { s/$/;osrc=$src/; }
        
        $isalt= (m/\b$altflag=(\d+)/) ? $1 : 0; #? or use predid =~ /t(\d+)$/;
        $isalt=1 if(m/status=alt-splice/); # upstatus=alt-splice addition # pasa
        if($altcomm) { $altcomm=0; $isalt=1 unless($isalt); }

        $issplit= (m/\b$splitflag=(\d+)/) ? $1 : 0; #? or use predid =~ /t(\d+)$/;

## FIXME: mustalt ids need to be preserved when corrected above (mustalt_check)
# bestgenes.DGILpub8d5.gff.mustup2:
# Scaffold1       EG2APu  gene    2254858 2288394   ID=acyp2eg0037017;oid=acyp2eg0000188 << oid is right here
# update altmodel        PASAgasmbl_1103
# mRNA ID=acyp2eg0037017t2;oid=acyp2eg0000188t2;gene=acyp2eg0037017;oidalt=PASAgasmbl_1103; << oid is right

        if($isalt) { 
          unless( m/;$altflag=$isalt/ or s/;$altflag=[^;\s]+/;$altflag=$isalt/ ) { s/$/;$altflag=$isalt/; }

      ## pasa upd alt-tr before prime tr; problems here w/ no prime pubid
          my $gid = $predid;
          my $altnum= ( $gid =~ s/t(\d+)$// ) ? $1 : $isalt+1; # "666"; #?? not working right
          #??  check altnum not already used per gene?
          
          my $aid= $pubidmap{$predid} || $pubidmap{$gid} || $pubidmap{$gid."t1"} || $pubidmap{$gid."t2"};
          # got a handful of dupl alt ids .. mix up b/n pasa update and userchoice isalt > pasa primary
          if($aid) { $aid =~ s/t\d+$//; } ## fixme below/here
          
          unless($aid) { 
            $aid= ($MID) ? $MID.$predid : ($pubgene) ? $pubgene : sprintf( $pubid_format, ++$pubidnum); 
          }
          $aid .= sprintf( $altid_format, $altnum); #? always append tNNN to prime id?
          
          $pubid= $aid;
          s/\bID=/ID=$pubid;$oidtag=/;  
          
        } elsif($MID) {
          $pubid= $MID.$predid;      
          s/\bID=/ID=$MID/;  
          
        } else {

      ## pasa upd alt-tr before prime tr; problems here w/ no prime pubid
      ## but problems here w/ non-alt but altid: AUGepir2s234g37t1,AUGepir2s234g37t10 > both EG2ap014820t1
      ##   AUGepir2s671g15t1,AUGepir2s671g15t15 > EG2ap030523t1
      ## .. not marked as alttr=1 but from augustus-alt-predict run
      
          (my $gid = $predid)  =~ s/t(\d+)$//;
          $pubid= $pubidmap{$predid} || $pubidmap{$gid} || $pubidmap{$gid."t1"} || $pubidmap{$gid."t2"} || "";

          if($pubid) {
            if($pubid =~ /t1$/) { $pubid=""; } else { $pubid =~ s/t(\d+)$/t1/; }
            $pubid="" if( $pubidmap{$pubid} ); 
          }
          
          if($pubid) { 
            # nada
          } elsif($pubgene) {
            $pubid = $pubgene . sprintf( $altid_format, 1);
          } else {
            $pubid= sprintf( $pubid_format.$altid_format, ++$pubidnum, 1); # always add tNNN?
          }
          
          s/\bID=/ID=$pubid;$oidtag=/;  
        }
        
        $pubidmap{$predid}= $pubid;      
        $pubidmap{$pubid}= $predid;
        
        if($dropannot) { s/\b($dropannot)=[^;\n]+;?//g ; s/\b($dropannot)=;//g; }
        # rename annots here?  pct_support > predsupport=
        
      } else {  # assume gene parts are contiguous?
        my ($pid)= m/\bParent=([^;\s]+)/;
        my $ppid= $pubidmap{$pid} || $pubid;
        if($MID) {
          s/\bParent=/Parent=$MID/;  
        } else {
          s/\bParent=/Parent=$ppid;$oidtag=/; 
        }
        if($dropexonannot) { s/\b($dropexonannot)=[^;\n]+;?//g ; s/\b($dropexonannot)=;//g; }
      }
 
    # option: move scorevec to ;scorevec= and scoresum= to  col.5
    if($movescorevec) {
      my $sv=$col[5];
      if($sv=~/,/) { 
        my ($ss) = m/scoresum=(\d+)/;        
        unless($ss) { map { $ss += $_ } split",",$sv; }
        s/$col[3]\t$col[4]\t$col[5]/$col[3]\t$col[4]\t$ss/;
        s/$/;scorevec=$sv/ unless($ss == 0 or $col[2] =~ /exon|CDS/);
      }
    }

    s/;;/;/g; # crap

    }
    
    print $outh $_; # ; to where ?? should be pipe-able
  } close($inh);

  print $outh "\n# source count for $MSRC:\n";
  foreach my $s (sort keys %srccount) { print $outh "# $s\t$srccount{$s}\n"; }
  if($debug) {
  print "\n# source count for $MSRC:\n";
  foreach my $s (sort keys %srccount) { print "# $s\t$srccount{$s}\n"; }
  }
  
  close($outh);
  return $outgff;
}



sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }

#   my @generec= sort _sortgene @$generecIN;
#  do 3-level sort: mRNA rec sorted by type; gene>mRNA sorted by alt-tr num?; chr>genes sorted by location
sub _sortgene  
{
  # rec == ($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr)
  # ($b->[2] cmp $a->[2]) # typ: reverse order by mRNA, exon, CDS  
  my($ta,$tb)= map{  (/gene/)?1:(/mRNA/)?2:(/CDS|exon/)?3:4; } ($a->[2],$b->[2]);
  return ($a->[0] cmp $b->[0]) # ref : should all be same for gene record; NOT for splitgene
      || ($ta <=> $tb)       # type-level
      || ($a->[3] <=> $b->[3]) # begin
      || ($b->[4] <=> $a->[4]) # end, largest 1st
      ;
}

sub _sortexonCDS  
{
  # rec == ($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr)
  # ($b->[2] cmp $a->[2]) # typ: reverse order by mRNA, exon, CDS  
  my($ta,$tb)= map{  (/gene/)?1:(/mRNA/)?2:(/exon/)?3:(/CDS/)?4:5; } ($a->[2],$b->[2]);
  return ($a->[0] cmp $b->[0]) # ref : should all be same for gene record; NOT for splitgene
      || ($a->[3] <=> $b->[3]) # begin
      || ($b->[4] <=> $a->[4]) # end, largest 1st
      || ($ta <=> $tb)       # type-level
      ;
}

sub _sortgeneOLD
{
  # rec == ($ref,$src,$typ,$tb,$te,$tscore,$to,$tph,$tattr)
  # ($b->[2] cmp $a->[2]) # typ: reverse order by mRNA, exon, CDS  
  my($ta,$tb)= map{  (/gene/)?1:(/mRNA/)?2:(/CDS|exon/)?3:4; } ($a->[2],$b->[2]);
  return ($ta <=> $tb)
      || ($a->[0] cmp $b->[0]) # ref : should all be same for gene record; NO, for splitgene
      || ($a->[3] <=> $b->[3]) # begin
      || ($b->[4] <=> $a->[4]) # end, largest 1st
      ;
}

sub _sortaltid {
  my($ta,$tb)= map{ (m/t(\d+)/) ? $1 : 0  } ($a,$b); 
  return ($ta <=> $tb || $a cmp $b);
}

# my $CHECK_MRNASPANS=1; # debug
# checkfixmrnaspan 

sub cds_fixspan {
  my($mrna,$exons)= @_;
  # FIXME here? CDS outside of exon span
  # 13 ERROR:   SEQ_FEAT.CDSmRNAXrefLocationProblem .. looks like these are CDS span exceeds parent exon.
  ## need to trim CDS to exon span,  before this mrna span check.
  
  # 1509upd: have dangling CDS from tiny exon drops, add opt to remove trailing/leading CDS w/o exon
  # need to change $exons array
  my $FIX_DANGLING_CDS=1; # debug
  my @cdrop=();
  my $newexons= undef;
  
  my @sexon= sort _sortexonCDS @$exons;
  my $nx=@sexon; my $didtrim=0;
  my($lex,$lxr,$lxt,$lxb,$lxe,$lxo)= (0) x 9;
  for(my $i=0; $i<$nx; $i++) {
    my $ex=$sexon[$i];
    my $nex= ($i<$nx-1)? $sexon[$i+1]: undef;
    my($xr,$xt,$xb,$xe,$xo)= @{$ex}[0,2,3,4,6];
    my($nxr,$nxt,$nxb,$nxe,$nxo)=(0) x 0;
    if($xt eq 'CDS') {
      my $bad=1; my ($trimb,$trime)=(0,0);
      if($lxt eq 'exon' and $xb>=$lxb and $xe <= $lxe) { $bad=0; }
      else {
        if($lxt eq 'exon' and $xe>$lxb and $xb<$lxb) { $trimb=$lxb; }
        if($lxt eq 'exon' and $xb<$lxe and $xe>$lxe) { $trime=$lxe; }
      }
      if($bad and $nex) {
        ($nxr,$nxt,$nxb,$nxe,$nxo)= @{$nex}[0,2,3,4,6];
        if($nxt eq 'exon' and $xb>=$nxb and $xe <= $nxe) { $bad=0; }
        elsif($trimb or $trime) { }
        else { 
        ## fixme: cds ovr both ends of exon.. trim=3,6?
          if($nxt eq 'exon' and  $xe>=$nxb and $xb<$nxb) { $trimb=$nxb; }
          if($nxt eq 'exon' and  $xb<=$nxe and $xe>$nxe) { $trime=$nxe; }
          }
      }
      if($bad) {
        $bad="bad";
        my($oxb,$oxe)=($xb,$xe); 
        if($trimb or $trime) {
          $xb=$trimb if($trimb); $xe=$trime if($trime);
          # if($trim==1) { $xb=$lxb; } elsif($trim==2) { $xe=$lxe; }
          # elsif($trim==3) { $xb=$nxb; } elsif($trim==4) { $xe=$nxe; }
          $ex->[3]= $xb; $ex->[4]= $xe; $didtrim++;   $bad="trimexcessCDS";
          
        } elsif($FIX_DANGLING_CDS) {
          # maybe more than end1? or internal exon drops?  
          # if($i < 1stexon or $i > lastexon or i-1==CDS and i+1==CDS) ??
          ## change to drop any CDS outside of bounding exon spans ? ie
          my $dodrop= ($i == 0 or $i == $nx-1)?1:0;
          $dodrop= 1 if($xb > $lxe and ($nxb==0 or $xe < $nxb)); # lxe=0 for i==0; nxb=0 for no next exon
          if($dodrop) { 
            push @cdrop, $i; 
            $trimb=$xb; $trime=$xe; $didtrim++; $bad="dropdanglingCDS";
          }
        }
        ##?? None fixed?
        my $xat=$ex->[8]; chomp($xat);
        warn "# cds_fixspan: $bad, trim=$trimb,$trime/$didtrim $xt.xbe=$oxb-$oxe, $lxt.lxbe=$lxb-$lxe, $xat\n";
        #.. 17 of these, 7 have asn (new) errs after trim (and before..)
      }
    }
    ($lex,$lxr,$lxt,$lxb,$lxe,$lxo)=($ex,$xr,$xt,$xb,$xe,$xo)
      if($xt eq 'exon'); #? only exons
  }
  if($FIX_DANGLING_CDS and @cdrop>0) {
    my @xnew=(); my %cdrop=map{ $_=>1 } @cdrop;
    for(my $i=0; $i<$nx; $i++) {
      my $ex= $sexon[$i]; 
      # debug alternate: keep as "#drop.$ex" ? need other changes
      push @xnew, $ex unless($cdrop{$i});
    }
    $newexons= \@xnew;
  }
  
  return($didtrim,$newexons);
}

sub mrna_fixspan {  # see findcds.pl 
  my($mrna,$exons)= @_;
  my $changed=0;
  return -1 unless(ref $mrna);
  return -2 unless(ref($exons) =~ /ARRAY/); # some mrna missing exons << TRACK, CURE? these
    ## oddly, none of miss-exon set has asn errors?? or no id mRNA errors?
  my($mref,$oldstart,$oldstop)= ($mrna->[0],$mrna->[3],$mrna->[4]);
  my($newstart,$newstop)= (0,0);
  
  ## UPD 1509 change exon array w/ dropped CDS
  ## .. for changed exons array, need to change in exon hash store:
  ##  mrna_fixspan($m,$exonfeat{$m}); 
  ## update: $exonfeat{$mrna}= $newexons; but caller holds exonfeat hash ..
  my ($ncfix,$newmrna,$newexons)= (0,undef,undef);
  ($ncfix, $newexons)= cds_fixspan($mrna,$exons);
  if($newexons and ref($newexons) =~ /ARRAY/) { $exons=$newexons; }
  $changed++ if($ncfix);
   
  foreach my $ex (@$exons) {
    my($xr,$b,$e)= @{$ex}[0,3,4];
    next if($xr ne $mref);
    ## bug split part mrna get all parts span on same chr .. match Split=n from mrna/exon?
    ##?? bug cds span > exon span, from findcds ?? Funhe2EKm004726t1
    $newstart= $b if($newstart == 0 or $b < $newstart);
    $newstop = $e if($newstop == 0 or $e > $newstop);
    }
  unless($newstop==0 or ($newstart == $oldstart and $newstop == $oldstop)) {
    $changed++; $mrna->[3] = $newstart; $mrna->[4] = $newstop; }
  #o# return $changed;
  return($changed,$newmrna,$newexons);
}


sub sortoutchr 
{
  my($outh, $atchr, $addgene, $geneidmap, $featlist, $commlist)= @_;
  my $header= (ref $commlist) ? delete $commlist->{"0header"} : "";
  my $chrsortlast  = $public_options{'chrsortlast'} || ""; # put these at end
  $chrsortlast =~ s/[,\s]+/\|/g if($chrsortlast);
  ## FIXME: sort Scaffold0002 1st vs contig0001 2nd; add chrsortsecond flag?
  my $chrsort2nd  = $public_options{'chrsort2nd'} || ""; # put these after others, eg contig 2nd, scaf 1st
    $chrsort2nd =~ s/[,\s]+/\|/g if($chrsort2nd);

  # 2012.jan: need new sort for genbanktbl submit: splitgene sort by chr, i.e. break gene bundle
  my $chrsortgene= $public_options{'chrsortgene'} || 0; # put these at end
  my $geneannot  = $public_options{'geneannot'} || ""; $geneannot =~ s/[\s,;]+/|/g;
  
  my $addloctag= $public_options{'addlocus'} || $general_config{'addlocus'} || 0; # NCBI special
  # $locustag{$loctag}{$gr} ||= $loctag; ..
  my %locustag;
  
  # fixme sigh for gene/mrna split parts.. 
  # repackage featlist into mrna+gene records only, + hash of exon/cds/other per mrna-partid
  my(%mrnafeat,%exonfeat);
  
  my @chrset=();
  foreach my $gid (sort keys %$geneidmap) 
  {
    my($gchr,$gchr1st,$gchr2nd,$gb,$ge,$gsc)=(0) x 10;
    my $gene= undef; # fixme: array for splits? got stutter now from split alts
    my %genechr=(); # one gene/gchr split fix?
    my @generec=();
    my %mrnaid;
    
    my @mid= @{$geneidmap->{$gid}}; # splitgene+alt problems .. need uniq mid if doing this way.
    my %mid= map{ $_,1 } @mid; # compress splitgene dup ids 
    @mid = sort _sortaltid keys %mid; #? use input order of mRNA, or sort on alt-id
     ## ^^? splitC1 case needs _sortaltidSplitC1 to separate altsA, altsB of split loc parts

    ## repackage featlist into mrna+gene records only, + hash of exon/cds/other per mrna-partid
    ## add bug patch here for bad mrna spans? from findcds splitgene exonfix, put same span on all same-chr parts
    ## my(%mrnafeat,%exonfeat);
    my($nsplit, $maxpart,%pchr)=(0,0);
    foreach my $mid (@mid) { 
      my @mall =  @{$featlist->{$mid}};
      #bug#my @mrna = grep { $_->[2] eq "mRNA" } @mall;
      my @mrna = grep { $_->[2] =~ m/$mrnatypes/ } @mall;
      if(@mrna>1) { $nsplit++; $maxpart= @mrna if(@mrna > $maxpart); }
      my ($mrnaid,$ip)=(0,0); 
      for my $xm (@mall) { 
        #bug#if( $xm->[2] eq "mRNA" ) 
        if( $xm->[2] =~ m/$mrnatypes/) 
        { ## ignore gene feat??
          $ip++; 
          # $mrnaid= join",",$gid,$mid,$ip; $mrnaid{$xm}= $mrnaid;  
          $mrnaid= $xm; # just use record ref as id?
          $mrnafeat{$gid}{$mid}{$ip}= $xm; 
	        #? $exonfeat{$mrnaid}=[];
          my $c= $xm->[0]; $pchr{$c}++; 
        } elsif( $xm->[2] ne "gene" ) { 
          push( @{$exonfeat{$mrnaid}}, $xm); 
        }
      }

      ## add bug patch here for bad mrna spans? from findcds splitgene exonfix, put same span on all same-chr parts
      ## .. only splits? or check/fix all?
      my $CHECK_MRNASPANS=1; # debug, use public_option
      if($CHECK_MRNASPANS) {
        my $nfix=0; my $nerr=0;
        for my $m (@mrna) { 
          #o# my $fixe= mrna_fixspan($m,$exonfeat{$m}); 
          my($fixe,$newmrna,$newexons)= mrna_fixspan($m,$exonfeat{$m});  # UPD 1509
          if($fixe<0) { $nerr++; } else { $nfix+=$fixe; } 
          # if($newmrna and ref($newmrna)) { $m= $newmrna; } #?? what else to change?  @mrna ..
          if($newexons and ref($newexons) =~ /ARRAY/) { $exonfeat{$m}= $newexons; }
          }
	      warn "# mrna_fixspan: nfix=$nfix, nerr=$nerr for mrna=$mid, gene=$gid\n" if($nfix>0 or $nerr>0);
      }
    }

    my $pchr= scalar(keys %pchr); $maxpart=$pchr if($pchr>$maxpart);
     
    my $loctag= $gid; $loctag =~ s/t\d+$//; 
    my($lastloctag,$lastgenepart)=(0,0);
    if($addloctag) { # prepare loctags, dang splits 
    
      ## DOUBLE DANG : this is a mess
      ## for splits need very careful sorting of alts + splits for proper ordering
      ## need a. separate chr/scafs, b. separate splits/chr, c. all alts of same chr (split1 or no split but overlap locus)
      ## .. simplify, dont expect/want cases of both a. & b. allow only a or b?
      
      # my($nsplit, $maxpart,%pchr)=(0,0);
      # foreach my $mid (@mid) { 
      #   my @mrna = grep { $_->[2] eq "mRNA" } @{$featlist->{$mid}};
      #   if(@mrna>1) { $nsplit++; $maxpart= @mrna if(@mrna > $maxpart); }
      #   for my $m (@mrna) { my $c= $m->[0]; $pchr{$c}++; }
      # }
      # my $pchr= scalar(keys %pchr); $maxpart=$pchr if($pchr>$maxpart);

      # if(1) { ## ready??...
      
      if($maxpart>1) { # dont bother any more for nonsplits 
        my $nspl=0; my(%pchr, %mchr, %psame); 
        foreach my $mid (@mid) { 
          my $chr1= 0; 
          #bug#my @mrna = grep { $_->[2] eq "mRNA" } @{$featlist->{$mid}};
          my @mrna = grep { $_->[2]  =~ m/$mrnatypes/ } @{$featlist->{$mid}};
          #^^ use now  $mrnafeat{$gid}{$mid}{$ip} ??
          @mrna= sort _sortgene @mrna; # may not be  locsorted on input, need this to ensure impart ordered?
          for my $i (0..$#mrna) {
            my $chr= $mrna[$i]->[0]; 
            $pchr{$chr}++; $psame{$i}++; 
            unless($chr1) { $chr1=$chr; } 
            else { $nspl++ if($i>0 and $chr eq $chr1); } #? count only splits on same chr?
            #? $psame{$i}{$chr}++; 
            #? $mchr{$mid}{$chr}++;
            # unless($chr1) { $chr1=$chr; } elsif($i>0) { $psame{$chr1}{$i}++ if($chr eq $chr1);  }
          }
        }
        
        my @c=keys %pchr;
        if(@c > 1) {
          @c= sort{ $pchr{$b} <=> $pchr{$a} } @c; # sort by most common chrs
          # if($nspl>0) { } # ignore this problem? keep only chr splits?
          for my $i (0..$#c) { 
            my $c= $c[$i];
            my $t=substr('ABCDEFGHIJKLMNOPQRSTUVWXYZ',$i,1);
            $locustag{$loctag}{$c} = $loctag.$t;
            }
        } elsif($nspl>0) {
          # MISSING gene idB from this set, splitC1/2 samescaf .. need to force new gene rec 
          my $c= $c[0];
          # my @i= sort { $psame{$b}{$c} <=> $psame{$a}{$c} } keys %psame;
          my @i= sort { $a<=>$b } keys %psame; # sort by i-index of @mrna seen above == impart
          for my $impart (@i) { 
            my $t=substr('ABCDEFGHIJKLMNOPQRSTUVWXYZ',$impart,1);
            $locustag{$loctag}{$impart}= $loctag.$t;
          }
        } # else { } # never here = commonest not split, default loctag
      }
      
      #} else {
      #         
      #   #/// old ///
      #   foreach my $mid (@mid) { 
      #     my $feat= $featlist->{$mid}; # splitgene, feat here should be several mRNA ?
      #     my @mrna = grep { $_->[2] eq "mRNA" } @$feat;
      #     ## BUG for Split type C[12] same scaf, $c == $c for each
      #     my $impart=0; for my $mrna (@mrna) { 
      #       #oo# $impart++; $locustag{$loctag}{$impart} = $loctag; 
      #       my $chr= $mrna->[0]; $locustag{$loctag}{$chr} = $loctag;   # ignore chr?
      #     }
      #   }
      #   
      #   my @ltg= sort keys %{$locustag{$loctag}};
      #   if(@ltg>1) {
      #     for(my $i=1; $i <= @ltg; $i++) {  
      #       my $lgr= $ltg[$i-1];
      #       my $c=substr('ABCDEFGHIJKLMNOPQRSTUVWXYZ',$i-1,1);
      #       $locustag{$loctag}{$lgr}= $loctag.$c;
      #     }
      #   }
      # } # ///old///

    }
    
    # FIXME: gene for splitgene; same mRNA ID diff scaffolds: need new gene rec for each?
    foreach my $mid (@mid) {
      my $feat= $featlist->{$mid}; # splitgene, feat here should be several mRNA ?
      ##  $mrnafeat{$gid}{$mid}{$ip}
      unless(ref $feat) { warn "# bad generec $gid:$mid\n"; next; }
      my $comm= $commlist->{$mid};
      push( @generec, $comm) if($comm);
      
      my @feat= sort _sortgene @$feat; # multiple mRNA/gene for splits
      #.. ^ bad sort for SplitC1, need resort by locustag
      #below.aftergene# push (@generec, @feat); # add mrna record
      
      # my $mrna= $feat[0]; # 1st always after _sortgene ? splitgene !
      #bug# my @mrna = grep { $_->[2] eq "mRNA" } @feat;
      my @mrna = grep { $_->[2] =~ m/$mrnatypes/ } @feat;
      #..^^.. splitgene probs here .. need uniq mid if doing this way.
 
      my $impart= -1; # use 0-index not 1-? 
      foreach my $mrna (@mrna) { 
        $impart++; # use 0-index not 1-? use impart to index locustag{} ?? or is it problem
        
        my ($ss) = $mrna->[8] =~ m/scoresum=(\d+)/;        
        unless($ss) { $ss=0; map { $ss += $_ } split",",$mrna->[5]; }
        $gsc=$ss if($ss>$gsc); # what is good gene score: max mRNA score?
        
        #FIXME: splitgene, mRNA ID on multiple ref locs; need gene same
        #x unless($gchr) { $gchr1st= $gchr= $mrna->[0]; $gb= $mrna->[3]; $ge= $mrna->[4]; }
        unless($gchr) { $gchr= $mrna->[0]; $gb= $mrna->[3]; $ge= $mrna->[4]; $gchr1st=$gchr unless($gchr1st); }
        if($gchr ne $mrna->[0]) { # only NOT-SO-rare splitgene ?
          if($gene) { $gene->[3]= $gb; $gene->[4]= $ge; $gene->[5]= $gsc; }
          $gchr= $mrna->[0]; $gb= $mrna->[3]; $ge= $mrna->[4];
          #x unless($gchr1st) { $gchr1st= $gchr; } else
          if($gchr ne $gchr1st) { $gchr2nd= $gchr; }
          $gene= $genechr{$gchr}{$gid}; 
          $gene=undef unless(ref $gene); #?? || undef; 
            #?? ^^ is this bad for gene overlaps need gid ?? but all genechr should be same gid

        } else { # gchr == mrna.chr
          # SPLIT C1/2 bug here, new mrna part, same gchr .. needs new gene rec AND do not extend gb,ge spans
          my $didgenespan=0;
          if( my $lt= $locustag{$loctag}{$impart} ) {
            my($genepart)= $lt =~ m/([A-Z])$/; $lastloctag= $lt;
            if($lastgenepart and $genepart ne $lastgenepart) { ## act like new gchr for genechr{gchr}
              if($gene) { $gene->[3]= $gb; $gene->[4]= $ge; $gene->[5]= $gsc; }
              if($lastgenepart and $genepart ne $lastgenepart and not $gchr2nd) { $gchr2nd= $genepart; } #?? yes or no
              $gene= $genechr{$genepart}{$gid}; 
              unless(ref $gene){ $gene= undef; $gb= 0; $ge= 0; }  
              else { $gb=$gene->[3]; $ge=$gene->[4];} # recover last span, fall thru to not didgrene
              #NO# $gb= $mrna->[3]; $ge= $mrna->[4]; $didgenespan++; # new part gene span
            }
            $lastgenepart=$genepart; 
          } 
          
          unless($didgenespan) {  # normal gene, expand gene span to cover mrnas          
            $gb= $mrna->[3] if($gb == 0 or $gb > $mrna->[3]);
            $ge= $mrna->[4] if($ge == 0 or $ge < $mrna->[4]);
          }
        }
        
        if($addgene) { $mrna->[8] =~ s/;gene=[^;\s]+//; $mrna->[8] =~ s/;/;gene=$gid;/; } # upd15:before copy2gene; BUT remove old**
   
        if($addloctag) { 
          ## BUG for Split type C[12] same scaf, $gchr.split1 == $gchr.split2 
          #o# my $lt=$locustag{$loctag}{$gchr}||$loctag;
          #oo# my $lt=$locustag{$loctag}{$impart}||$loctag;
          my $lt= $locustag{$loctag}{$gchr}; # this works when splits are only cross-chr; both ways possible
          $lt= $locustag{$loctag}{$impart} unless($lt); # impart is problematic, need better part-id tag on each mrna part
          $lt= $loctag unless($lt);
          $lastloctag= $lt;
          
          $mrna->[8] =~ s/;locustag=[^;\s]+//;  $mrna->[8] =~ s/;/;locustag=$lt;/; 

use constant FIXLOCUSTAGSPLIT => 1; #UPD1509. correct Split= for locustag part
	  if(FIXLOCUSTAGSPLIT) {
	    my($genepart)= $lt =~ m/([A-Z])$/;
            my($spl)= $mrna->[8] =~ m/;Split=(\d[^;\s]*)/;
            if($spl and not $genepart) { 
	      $mrna->[8] =~ s/;Split=$spl/;unspl=$spl/; 
            } elsif($genepart and not $spl) { 
	      ## not always right, have genepart for unsplit mrna when other alt mrna is split
	      # my $ipart= 1 + index('ABCDEFGHIJKLMNOPQRSTUVWXYZ',$genepart);
	      # $mrna->[8] =~ s/locustag=$lt;/locustag=$lt;Split=$ipart;/;	
            }
	  }

        }
          
        unless($gene) { # dangit this needs to be outside (@mrna) ?? or not, 1 mrna == @mrna also
          my @gener= @{$mrna}; $gene= \@gener;
          $gene->[2]="gene"; 
          $gene->[5]=1; # score 
          my $ga="";
          if($geneannot) {
            my @ga= $gene->[8] =~ m/\b((?:$geneannot)=[^;\n]+);?/g;
            $ga= (@ga) ? ";".join(";",@ga) : "";
          }
          $gene->[8]= "ID=$gid$ga\n"; # FIXME: Split=1/2 splitgene needed on gene from mRNA
          
          $genechr{$gchr}{$gid}=$gene unless($genechr{$gchr}{$gid});
          if($lastgenepart) { $genechr{$lastgenepart}{$gid}= $gene; }
          #?? handle other split case C1/C2 same gchr, using locustag not gchr? or locustag instead of gid?
          #?? bug here for genepart sort order, 2 gene parts before 2 mrna parts bad. also for mid altsA > altsB
          push( @generec, $gene) if($addgene); #here, allow multiple same ID for splitgene
          }
          
        push( @generec, $mrna); # NEW; only mrna,gene in @generec, pull ordered exons from %exonfeat
        } # @mrna splitgene loop
        
      #OLD# push (@generec, @feat); # add mrna record; should this be @generec, \@feat ? keep each @mrna-part together?
    }
    
    if($gene) {
      $gene->[3]= $gb; $gene->[4]= $ge; $gene->[5]= $gsc;  #? bad for splits??
      #above.beforeft# unshift( @generec, $gene) if($addgene);
    }
    
    if(@generec) {  
      my $chr = $gchr1st; # $gene->[0]; # or use $atchr param? 
      
      ## FIXME: splitgene sort bug:  sc1:99990-99999/sc9:111-222 sorts before sc1:333-444, should be after
      ## FIXME2: BAD sort here, splits, breaks up gene records ??? eg. Funhe2EKm037931,Funhe2EKm012442 split & colocated
      ## .. or maybe generec packaging bad.. eg diff locus ids w/ same exon locs > interleaved output sort
      
      #o# if($chrsortgene and $gchr1st ne $gchr) ## THIS is bad test for splits + alts, ie chr1 > chr2 > chr1 == gchr1st
      ## BUG dupl gene records in split set, needs chrsortgene before adding gene recs ??
      ## sort bug for $genechr{$lastgenepart}{$gid} now, put gene idA, idB both before mrnaA, mrnaB
      ## .. need to treat genepart == A,B,C.. like chr for sort >> ichr, chr, ipart, part?
      ## separate sortSplitC1 method needed .. if $gchr2nd =~ /^[A-Z]$/ .. implies locustag on each mrna/gene?
    
      use constant { ICHR2ND => 5000000, ICROTHER => 8888888, ICRLAST => 9999999 };          

      if($chrsortgene and $gchr2nd)  
      {  
        my $issplitC1= ($gchr2nd =~ m/^[A-Z]$/)?1:0; # handle above same way?
        
        ## # this way is right.
        my %grparts; # need to collect all parts/chr before push to chrset !!
        foreach my $gr (@generec) {
          next unless(ref($gr) =~ /ARRAY/); # warn "# bad split-gene rec $gid = '$gr'\n"; 
          my($cr,$cb,$ce,$cat)= @{$gr}[0,3,4,8]; 
          my $part=$cr;
          if($issplitC1){ my($lt)= $cat=~/locustag=(\w+)/; ($part)= $lt=~m/([A-Z])$/; $part||=1; }
          push( @{$grparts{$part}}, $gr);
        }
        for my $cr (sort keys %grparts) {
          my $grpart= $grparts{$cr};
          my @gr= @{$grpart};
          my($chr,$gb,$ge)= @{$gr[0]}[0,3,4]; # other parts gb,ge ?? should be enclosing gene span/part
          
          my $ichr= ($chrsortlast and $chr =~ m/^($chrsortlast)$/)? ICRLAST : ($chr =~ m/(\d+)/) ? $1 : ICROTHER;  
          $ichr += ICHR2ND if($chrsort2nd and $chr =~ m/^($chrsort2nd)/);
          push(@chrset, [ $ichr,$chr,$gb,$ge,$grpart] ); # break @generec by chr
        }  
        @generec=(); # done
        
        # } else { # causes bad sort interleaved gene records
        #    my @gr=();  ($chr,$gb,$ge)=(0) x 3;
        #    foreach my $gr (@generec) {
        #      ## BUG here: @{$gr} NOT ARRAY .. what data is this
        #      unless(ref($gr) =~ /ARRAY/) { next; } # warn "# bad split-gene rec $gid = '$gr'\n"; 
        #      # # bad split-gene rec Funhe2EKm024308 = '  what is this crap?? extra newlines?
        # 
        #      my($cr,$cb,$ce)= @{$gr}[0,3,4];
        #      if($chr and $chr ne $cr) { 
        #        my $ichr= ($chrsortlast and $chr =~ m/^($chrsortlast)$/)? 9999999 : ($chr =~ m/(\d+)/) ? $1 : 8888888;  
        #        my @grclone= @gr;
        #        push(@chrset, [ $ichr,$chr,$gb,$ge,\@grclone] ); # break @generec by chr
        #        @gr=(); ($chr,$gb,$ge)=(0) x 3;
        #      }
        #      push @gr, $gr;
        #      ($chr,$gb,$ge)=($cr,$cb,$ce) unless($chr); #?? chr/cr change resets this
        #    }
        #    @generec= @gr;
        # }        
      }
      
      if(@generec) {
        my $ichr= ($chrsortlast and $chr =~ m/^($chrsortlast)$/)? ICRLAST : ($chr =~ m/(\d+)/) ? $1 : ICROTHER;  
        $ichr += ICHR2ND if($chrsort2nd and $chr =~ m/^($chrsort2nd)/);
        push(@chrset, [ $ichr,$chr,$gb,$ge,\@generec] ); # sort these for output
      }
      
    } else {
      # missing mrna data ??
    }
  }
 
  # my $chrsortlast  = $public_options{'chrsortlast'} || ""; # put these at end
  # output 3-tier sorted genes
  print $outh $header if($header);
  foreach my $gene (
    sort { $a->[0] <=> $b->[0] # ichr
      || $a->[1] cmp $b->[1]  #  same ichr, diff chr
      || $a->[2] <=> $b->[2] || $b->[3] <=> $a->[3] } # by location; but add real chr cmp
    @chrset) {
    my $feats= $gene->[4];
    
    foreach my $feat (@$feats) {
      ## insert %mrnafeat,%exonfeat here, exonfeat has ordered feats for each mrna/part
      ## cant now handle chr/part resort of mrnafeat with exonfeats ..
      ## @$feats should now become only ordered gene,mrna. pick up exons from each mrna.partid
      # my $mrnaid= $mrnaid{$feat}; 
      
      if(ref $feat) {
        print $outh join("\t", @$feat); ## ,"\n";
        if(my $exons= $exonfeat{$feat}) { #NEW
          for my $x (@$exons) { print $outh join("\t", @$x); }
        }
      } elsif($feat) {
        print $outh $feat; # comments? assume newlines? any errs?
        print $outh "\n" unless($feat =~ m/\n$/);
      }  
    }
  }

}

# sub addgene
# {
#   my($generecIN)= @_;  # , $overlaplist, $addgene ?
#   
#   # addgene: scan for all mRNA in generec, collect max span, geneid; update mRNA w/ gene= or Parent=
#   
#   my @generec= sort _sortgene @$generecIN; #? defer sort til output
#   my $generec= \@generec;
# 
#   my $mrna= $generec->[0]; #?? or check type? maybe gene ? addgene here?
#   
#   my($ref,$src,$typ,$tb,$te,$tp,$to,$tph,$tattr)= @$mrna;
#   my @bins= (int($tb/$BINSIZE) .. int($te/$BINSIZE));
#   foreach my $ib (@bins) { push( @{$overlaplist->{$ref}{$ib}}, $generec); }
# }

sub evigene_idparts {  # upd18apr evgIdParts($id)
  #     my($gpre,$gnum,$gd,$ti,$isplit)= evigene_idparts($id);
  my($id)=@_;
  my($gpre,$gnum,$gd,$ti,$gdup,$isplit)=(0) x 9;
  $gd=$id; 
  ## FIXME for _G2 _C2? tags, G2,n is diff locus id
  ($isplit)= ($gd=~s/_(C\d+)$//)?$1:0;
  ($gdup)  = ($gd=~s/_(G\d+)$//)?$1:0;
  
  # bug: Danrer6pEVm003029t41/main Danrer6pEVm003029t41t2, Danrer6pEVm003029t41t3 altids
  if($gd =~ m/^(\w+[A-Za-su-z])(\d\d+)t(\d+)$/) { 
    ($gpre,$gnum,$ti)=($1,$2,$3);
    $gd= $gpre.$gnum;
  } else {
    $gd =~ s/[a-su-z]\d+$//; # non-t suffix ??? ie. t2d33
    ($ti) = ($gd=~s/t(\d+)$//)?$1:0;  # BUG: m00001t12t3 hack format
    if( $gd =~m/^(\w+[A-Za-su-z])(\d\d+)/) { ($gpre,$gnum)=($1,$2); $gd=$gpre.$gnum; }
    else { # shouldnt be here?
      $gd =~ s/t\d+$//; # extra tNtN..
      ($gnum) = ($gd=~m/(\d+)$/) ? $1:0;  # BUG.. ? set gnum=0 meaning need new id?
      ($gpre=$gd) =~ s/$gnum$//;
      $gnum=0; # BUG fix??
    }
  }
  
  if($gdup) { $gd.="_".$gdup; $gnum=0; } # gnum invalid?
  return($gpre,$gnum,$gd,$ti,$isplit);
  #OR:   return($gpre,$gnum,$ti,$isplit,$gdup); #?
}

sub genesort_update 
{
  my($ingff, $outgff)= @_;

  my ($pubid, $predid, $inid, $gid, $pubidnum, $ngene,$nmrna,$nexon,$nother)= (0) x 10;
  my (%pubidmap, %geneidmap, %featlist, %commlist); # include alt-tr ?
 
  my $addgene= $public_options{'addgene'} || $general_config{'addgene'} || 0; # i.e. add gene record enclosing mRNA (for alttr)

  #?? my $chrindex = $evidence{'chrindex'} || ""; # sort chr if this exists
  my $chrsort  = $public_options{'chrsort'} || 0; # sort chr by numbering
  # 2011.sep: chrsort fails (for wasp genes)... why?
  # 2012.jan: need new sort for genbanktbl submit: splitgenes sort by chr, i.e. break gene bundle
  
use constant DEBUG_SORT => 0;  
  
  my $pubid_format = $public_options{'publicid'} || "CG0000000";
  my $altid_format = $public_options{'altid'} || "t00";
  my $altflag = $public_options{'altflag'} || "alttr";
  my $splitflag=  $public_options{'splitflag'} || "Split";
  ## more split flags: splitgenesort = by_gene/by_chr; idsplit_is_done=0/1; idsplit_suffix=_Cn|Sn|.n
  my $splitsort=  $public_options{'splitgenesort'} || "by_gene";
  my $splitsuf =  $public_options{'idsplit_suffix'} || "_C";
  if($splitsort =~ /chr/) { $public_options{'chrsortgene'}=1; } # see above in sortoutchr
  my $oidtag= "oid";  


  #**** Problem here ; royal id screwup from asmrna with t12345 ids, but not alts here
  #***  pps, have two gene=ID in mRNA:  gene=evigID;gene=oidID; dont want 2nd
  my $ignore_source_altid= $public_options{'ignore_source_altid'} || 0; # 2011sep problem w/ source tNNN 

  my $nochangeid = $public_options{'nochangeid'} || 0;
  $ignore_source_altid=0 if($nochangeid); # not both
  
  my $nd= ( $pubid_format =~ s/(0\d+)// ) ? length($1) : 6;
  $pubid_format .= '%0'.$nd.'d'; # CG%06d
  $nd= ( $altid_format =~ s/(\d+)// ) ? length($1) : 1;
  $altid_format .= '%0'.$nd.'d'; # t%02d
  
  $outgff="$ingff.so" unless($outgff);
  my($inh, $outh)= openio($ingff, $outgff);
  
  my ($altcomm, $isalt, $issplit, $lchr, $comsave)= (0) x 10;
  $pubid= "0header";
  
  ## this drops comments at end: recover them
  print "# genesort_update: $ingff TO $outgff\n" if($debug);
  while(<$inh>) {
  
    if(/^#x/) { next; } # drop excluded gff here
    elsif(/^#a\./) { $altcomm=1; s/^#a\.//; } # alttr comment ??    
    #?? possibly save comments, output with gene record occurence ?
    
    if(/^\w/) {
 
      my @col= split"\t";  
      
      #? collect genes/chr then sort/dump at new chr ? dont assume gene-rec sorted? but chr sorted?
      # output location sorted gene-recs per chr
      # ? add chr sorting? using conf: chrsort = genome/aphid2asm.fa.count

      if( ! $chrsort and $lchr and $lchr ne $col[0]) {
        sortoutchr( $outh, $lchr, $addgene, \%geneidmap, \%featlist, \%commlist);
        %geneidmap= %featlist= %commlist= ();
      }

      $lchr= $col[0];
      my $rloc= \@col; 
      $issplit= (m/\b$splitflag=(\w+)/) ? $1 : 0; # do all items have this now?
      
      if( $col[2] eq "gene" ) { 
        # save it?, use it
        
      } elsif( $col[2] =~ m/^($mrnatypes)/) { 
        ($predid)= m/\bID=([^;\s]+)/; $inid=$predid;

        my($genepre,$genenum,$geneid,$ti,$idsplit)= evigene_idparts($predid); #upd18apr
        if($idsplit) { $issplit=$idsplit; $predid=$geneid."t$ti"; } #??
        
        $predid =~ s/t(\d+)$/u$1/ if($ignore_source_altid); # unless m/;must=69/;
        
        $isalt= (m/\b$altflag=(\d+)/) ? $1 : 0; #? or use predid =~ /t(\d+)$/;
        $isalt=1 if($isalt==0 and m/status=alt-splice/); # upstatus=alt-splice addition # pasa
        $isalt=0 if($ignore_source_altid and not m/;must=69/);
        if($altcomm) { $altcomm=0; $isalt=1 unless($isalt); }
        
        # $issplit= (m/\b$splitflag=(\w+)/) ? $1 : 0;
        
        if($nochangeid) {
        
          $pubid= $predid;
          $gid  = $geneid; # $predid; $gid =~ s/t(\d+)$//; my $ti=$1;
          $ngene++ unless($ti > 1);
          
        } elsif($isalt) { 

          $gid = $predid;
          my $altnum= ( $gid =~ s/t(\d+)$// ) ? $1 : $isalt+1; # "666"; #?? not working right
          #??  check altnum not already used per gene?
          
          my $aid= $pubidmap{$predid} || $pubidmap{$gid} || $pubidmap{$gid."t1"} || $pubidmap{$gid."t2"};
          if($aid) { $aid =~ s/t\d+$//; } ## fixme below/here
          
          unless($aid) { $aid= ($MID) ? $MID.$predid : sprintf( $pubid_format, ++$pubidnum); }  
          $aid .= sprintf( $altid_format, $altnum); #? always append tNNN to prime id?          
          $pubid= $aid;
          #? $pubid ||= $predid;
          
        } else {
          $ngene++;
          
          $gid = $predid;
          $gid =~ s/t(\d+)$//;
          $pubid= $pubidmap{$predid} || $pubidmap{$gid} || $pubidmap{$gid."t1"} || $pubidmap{$gid."t2"} || "";

          if($pubid) {
            if($pubid =~ /t1$/) { $pubid=""; } else { $pubid =~ s/t(\d+)$/t1/; }
            $pubid="" if( $pubidmap{$pubid} ); 
          }
          
          $pubid= sprintf( $pubid_format.$altid_format, ++$pubidnum, 1) unless($pubid); # always add tNNN?
          #? $pubid ||= $predid;
        }
        
        # if($issplit and $pubid !~ /$splitsuf\d/) {  # see evigene2genotbl_kfish2.pl IDSplitIsDone
        #   $pubid .= $splitsuf.$issplit;  # FIX 201504 ncbi needs split parts w/  uniq ids
        # }
        
        $pubidmap{$predid}= $pubid; # ok for same. # splits? need predid changed?
        $pubidmap{$pubid}= $predid;
        $nmrna++; 
        #push( @{$geneidmap{$gid}}, $pubid); # pubid? splitgene dup ids here **
        
        if($featlist{$pubid}) {  # splitgene w/ alts.. problems
          # FIXME: uniq ids here for splitgene parts ?? no good for exons
          # my $pubid1="$pubid.".(++$pubidnum); # $nmrna  ## dang; exons need this also.
          # $featlist{$pubid1}= [$rloc]; 
          push( @{$featlist{$pubid}},$rloc);
          push( @{$geneidmap{$gid}}, $pubid); # pubid? splitgene dup ids here **
          warn "# genesort: mRNA duplicate for $pubid \n" unless($issplit); 
        } else { 
          $featlist{$pubid}= [$rloc]; 
          push( @{$geneidmap{$gid}}, $pubid); # pubid? splitgene dup ids here **
        } #new mrna;  
        
        $commlist{$pubid}.= $comsave if($comsave); $comsave=""; # really want to push preceding comms to following pubid
        
      } elsif( $col[2] =~ m/^($exontypes)/ ) {
        my ($pid)= m/\bParent=([^;\s]+)/;

        my($genepre,$genenum,$geneid,$ti,$idsplit)= evigene_idparts($pid); #upd18apr
        if($idsplit) { $issplit=$idsplit; $pid=$geneid."t$ti"; } #??

        $pid =~ s/t(\d+)$/u$1/ if($ignore_source_altid);
        
        ## but see evigene2genotbl_kfish2.pl
        # if($issplit and $pid !~ /$splitsuf\d/) { 
        #   $pid .= $splitsuf.$issplit;  # FIX 201504 ncbi needs split parts w/  uniq ids
        # }
        
        my $ppid= ($nochangeid) ? $pid : $pubidmap{$pid} || $pubid;
        
        if($pid ne $predid) {
          #? dont care# warn "# out of order gene rec: $pid.$col[2] not in mRNA:$predid\n"; # and do what? any of these?
        }

        $nexon++;
        push( @{$featlist{$ppid}},$rloc);
        
      } else { # ?? save in gene rec or not?
        $nother++;
        if(@col > 6 and $col[0] eq $lchr) {
          push( @{$featlist{$pubid}},$rloc);
        }
      }
 
    } else { # comment: push to generec
      if($pubid eq "0header") { $commlist{"0header"}.= $_; }
      else { $comsave.= $_;}
      # really want to push preceding comms to following pubid
    }
    
  } close($inh);
  
  sortoutchr( $outh, $lchr, $addgene, \%geneidmap, \%featlist, \%commlist);
  print $outh $comsave if($comsave);
  
  print $outh "\n# gene count for $MSRC: ngene=$ngene, nmrna=$nmrna, nexon=$nexon, nother=$nother\n";
  print       "\n# gene count for $MSRC: ngene=$ngene, nmrna=$nmrna, nexon=$nexon, nother=$nother\n"
    if $debug;
  
  close($outh);
  return $outgff;
}



sub evigene_config {  # this belongs in evigene.pm package
  my($cfile, $addoptions)= @_;
  my $ctype=0;

  use constant{ kEVFILE => 1, kEVOPT => 2, kANOPT => 3, kEVPROG => 4, kPUBOPT => 5, kEVGENES => 6, };
  return unless($cfile or $addoptions); #  and -f $cfile
  
  my @CONFIG=();
  if($cfile) {
    open(F,$cfile) or die "ERROR reading config: $cfile";
    @CONFIG= <F>; close(F);
  }
  if($addoptions) {
    # FIXME: allow addoptions only, parse many: @addopt= split /[,;|],$addopt
    # -cadd 'pubopt=1,addgene=1,chrsort=1,nochangeid=1,chrsort2nd=GRC'
    my @aconf= (ref $addoptions) ? @$addoptions : ($addoptions);
    for my $ac (@aconf) {
      my @ac= split /[,;|]/, $ac;
      push @CONFIG, @ac ;
    }
  }
  
  my ($lastkey, $lastval);
  foreach (@CONFIG) {  # key => value
    s/^\s+//;  s/\#.*$//; s/\s+$//; # end of line comments 

  ## need now to handle continuation lines, end with \
  ## or prefix with "+" << use this

    my($key,$val);
    if(/^[+\.]/ and $lastkey) { $key= $lastkey; s/^.//; $val= $lastval.$_; }
    elsif($lastval =~ m,\\$, and $lastkey) { $key= $lastkey; $lastval=~s,\\$,,; $val= $lastval.$_; }
    else { ($key,$val)= split(/[=\s]+/, $_, 2); }
    # ($key,$val)= split" ",$_,2; # space in option values : or split /[=\s]+/, $_, 2
    
    next unless($key =~ /^\w/); 
    $val =~ s/\s*$//;
    # verbose "# config k=v: $key=$val";
    # FIXME: need FindBin / sub findevigeneapp() for these now
    if($ctype == kEVPROG and $val =~ m,^scripts,) { $val =~ s,scripts,$EVIGENES,; }

    if($key =~ /^evidence/) { $ctype= kEVFILE; }
    elsif($key =~ /^evoption/) { $ctype= kEVOPT; }
    elsif($key =~ /^anoption/) { $ctype= kANOPT; }
    elsif($key =~ /^pubopt/) { $ctype= kPUBOPT; }
    elsif($key =~ /^geneset/) { $ctype= kEVGENES; }
    elsif($key =~ /^program/) { $ctype= kEVPROG; }
    elsif($key =~ /^end$/) { $ctype= 0; }
    ##elsif($ctype == kEVFILE or $val =~ /gff/) { $evidence{$key}= $val; } 
    elsif($ctype == kEVFILE ) { $evidence{$key}= $val; } # /gff/ was bad for geneset
    elsif($ctype == 0 and $val =~ /\.gff/) { $evidence{$key}= $val; } 
    elsif($ctype == kEVOPT) { $evaluate_options{$key}= $val; } 
    elsif($ctype == kANOPT ) { $annotate_options{$key}= $val; } 
    elsif($ctype == kPUBOPT ) { $public_options{$key}= $val; } 
    elsif($ctype == kEVGENES) { $geneset{$key}= $val; }
    
    elsif($key =~ /^evkey|^evalkey/) {  @evalkeys= split( /[\s,;]+/, $val); }
    elsif($key =~ /^ankey|^annotkey/) {  @annotkeys= split( /[\s,;]+/, $val); }
    elsif($key =~ /^evmore|^evalmore/) {  @moreevalkeys= split( /[\s,;]+/, $val); }
    elsif($key =~ /^anmore|^annotmore/) {  @moreannotkeys= split( /[\s,;]+/, $val); }

      ## these now in %evidence hash
    #elsif($key =~ /^genescore/) { $genescoredir= $val; }  # for annotation add_genescore()
    #elsif($key =~ /^blastself/) { $blastself= $val; }  # for annotation add_blastscores()
    #elsif($key =~ /^blastother/) { $blastother= $val; }  # for annotation add_blastscores()

    #?? change to program hash list? use eval_options to set program?
    # elsif($ctype == kEVPROG  ) 
    elsif($val =~ /overlapfilter/ or $key =~ /^overlapfilter/) { $overlapfilter=$val; }
    elsif($val =~ /overlapeval/ or $key =~ /^overlapeval/) { $overlapeval=$val; }
    elsif($val =~ /overbestgenes/ or $key =~ /^overbestgenes/) { $overbestgenes=$val; }

    # generic keys: name date genome .. other?
    elsif($key =~ /\w/ and $val =~ /\w/) { $general_config{$key}= $val; }
    
    # also for now : overlap, other progs
    if($ctype == kEVPROG) { $programs{$key}= $val; }
    
    ($lastkey, $lastval)=($key, $val);
    
  } 
  
  
#   die "ERROR: missing overlapfilter program: $overlapfilter"
#       unless( -x $overlapfilter);
#   die "ERROR: missing overlapeval program: $overlapeval"
#       if($overlapeval and not -x $overlapeval); # not used by annot

}
