# evigene_pubsets.pm

# package evigene_pubsets;
package main;

=item about package

  reusable subs from 
     evigene/scripts/evgmrna2tsa2.pl
  
=cut

use constant evigene_pubsets_VERSION => '2018.02.20'; # create from evgmrna2tsa2.pl subs

use FindBin;
use lib ("$FindBin::Bin"); # assume evigene/scripts/ layout: main, prot/, rnaseq/, ...

use strict;
use warnings;
# use FindBin;
# use File::Basename qw(basename dirname fileparse);

use cdna_evigenesub;  

#? require Exporter;
# our @ISA = qw (Exporter);
# our @EXPORT = qw (xxxx);

# cdna_evigenesub globals: 
# vars qw ( $EVIGENES $EGAPP $GDB_PREFIX $DEBUG $dryrun
#   $MIN_NAMEIDENT $MIN_IDLIKE $EVGLOGH 
#   %genenames %genedbxref %genenamepct %namedgenes %cddnames 
#   %pubids %pubidinfo $APPtraa2cds
#   $AAQUALF $AAQUALH $BAD_GAPS
#   );
# 

##  productname : use protein_names.pm
use constant FOR_NCBI => 1; # format flags
use constant FOR_NOTE => 2;
our $NAME_NONE = "Unknown|Uncharacterized conserved|Uncharacterized|Hypothetical";  
our $NAME_UNK  = 'Uncharacterized protein'; # uniprot  
our $NAME_UNKNCBI= 'hypothetical protein';  # NCBI prefers this to Uncharacterized .. special handling
our $RNATYPES='mRNA|ncRNA'; # transcribed_RNA mrnatypes exontypes 
our $EXONTYPES='exon|CDS'; # mrnatypes exontypes 

# see genes/evganngff.pl
# our $MDROP='match|protein|indels|cdsindel|cdsspan|chimera|ocds|aaold|xtrim|offold|cdsfix|cdnabest'; #not qlen
our $ANNCUTRNA='pho|lspan|scorevec|match|indels|cdsindel|cdsspan|ocds|aaold|xtrim|offold|cdsfix'; # = MDROP mrna rec
our $ANNCUTEXON='intr'; #  

#fixme chomp * on prot ends .. from where? uvcut= subset
use constant AA_NOSTOPCODON => 1;

our $IDPREFIX; # caller sets = $ENV{idprefix} || 'EVGm'; 
our $ORGANISM; # caller sets
sub okORGANISM{ return($ORGANISM and $ORGANISM =~ m/\w\w/ and $ORGANISM ne "Noname"); }
our $APPgenesupdate; # = findevigeneapp("bestgenes_update.pl")   

our $DATE; # caller sets
# $DATE=`date '+%Y%m%d'`; chomp($DATE); # default is today; use perl func?

our $CULLSEP= 1; # option = separate output cull from pub, seq + gff
our $SORT_PUBSET= 1; # option..
our $FAHDRisValid; #= 0; # for annotab2tblinfo annot.tab vs mRNA header conflicts
our $skipdropseqs=1; # make_pubseq fixup flag, default on/off? special cases of changing trclass 
our($pubid_format,$altid_format, $pubidnum_start, $GBPROID);

our (@tmpfiles,@publicset); # tidyup file sets, also @okayset,@dropset,@inputset,,@erasefiles
our (%pubids,%pubidinfo,%puboidc,%puballoids); ## cdna_evigenesub globals 

use constant { kdDROP => -999, kdCULL => -1, kdOK => 1, kdMUSTKEEP => 2, kdOther => 0 };  
our ($KEEPDROPX,%keepdropxtab); # evigene_pubsets.pm
# my($nkeepdrop,$keepdrop_idhash)= readKeepDrop($keepdropin) if($keepdropin); # fill  global %keepdrop
    #   $keepdropxtab{act}{$moid}=$act;  # use act =~ m/(keep|ok|drop|cull)/ to filter dup pubids
    #   $keepdropxtab{mapid}{$moid}=$pubid;
    #   $keepdropxtab{act}{$oid}=$act;
    #   $keepdropxtab{seqid}{$oid}=$pubid;

use constant NAME_NOEDITS => 1;
use constant CUT_IDCG => 0; # (my $soid=$oid) =~ s/_[CG]\d+$//; causes problems, xxxt1_G2 diff pubid than xxxt1
use constant splitKEEP_C => 1; # for now output gff IDs with _C split tags
use constant nogffDUPIDS => 1; # dup input.gff records? filter here?

use constant ANNO_ADD_MAPCOLS => 1;
#our @ANNO_COLS= qw(PublicID OrigID TrLen CDSoff AAqual TrGaps Dbxref Namealign Product_Name CDD_Name Class);
#our @ANNO_TBLINFO = qw(pubid oid trlen cdsoff aaqual trgaps dbxref namepct name cdd evgclass);
our @ANNO_COLS= qw(PublicID OrigID TrLen CDSoff AAqual TrGaps Dbxref Namealign Product_Name CDD_Name Class Maploc Mapqual);
our @ANNO_TBLINFO = qw(pubid oid trlen cdsoff aaqual trgaps dbxref namepct name cdd evgclass maploc mapqual);

# BEGIN { #? package begin fails here
#   if(ANNO_ADD_MAPCOLS) {
#     push @ANNO_COLS, qw(Maploc Mapqual);
#     push @ANNO_TBLINFO, qw(maploc mapqual);
#   }
# }

=item read_pubids

  read gene annots from ONE OF pubidtab or pubseq (if no pubtab)
  ($nred, $pubids_hashref)= read_pubids(pubidtab,pubseq)
  
  %pubids is global in cdna_evigenesub
  -- move read_pubids to cdna_evigenesub ?
  -- extend to read from mrna/aa publicset seq files instead, if no .pubids?
  
  realtids has extended fields, new drop/keep classing .. dmagset56i.realt_pubids
  now publicset/name.pubids looks this way, preserve cols 1..4 as prior version; 5=Class is also useful

  #Public_mRNA_ID           originalID                      PublicGeneID            AltNum  Class   AAqual  pIdAln  Notes
  # Dapma6iEVm000001t1      dmag4vel4ifik65Loc753t140       Dapma6iEVm000001        1       main    18640,93%,complete,44X  99/99/. .
  # Dapma6iEVm000001t2      hsICO2931no9gvelvk53Loc431t59   Dapma6iEVm000001        2       althi   17246,98%,complete      99/67/./dmag4vel4ifik65Loc753t141       oldid:Dapma6iEVm000001t10,
  # Dapma6iEVm000001t28     dmag4vel4ipak35Loc6068t3        Dapma6iEVm000001        28      althi   51,73%,complete 99/83/./dmag4vel4ifik65Loc753t141       oldid:Dapma6iEVm000001t8,
  # Dapma6iEVm000001t29     dmag4vel4ipak65Loc1003t146      Dapma6iEVm000001        29      dropalthi1      13944,65%,partial3      99/99/./dmag4vel4ifik65Loc753t141       oldid:Dapma6iEVm000001t18,

=cut

sub read_pubids 
{
  my($pubidtab,$pubseq,$skipclass)= @_;
  # -- extend to read from mrna/aa publicset seq files instead, if no .pubids?
  #  read_pubids( makename($cdnaseq,".pubids"), $cdnaseq); # CHECK for file, or read from annotab ??

  $skipclass ||= 'drop'; # 'drop|cull'; # drop only?
  $pubseq||="";
  
  my($ok,$hin,$isseq,$infile,$nred,$ndupskip,$ndup1,$ndropskip,$isrealt,$realtdrops,$nalltr)=(0) x 19;
  if($pubseq and not -s $pubidtab and -s $pubseq) {
    ($ok,$hin)= openRead($pubseq); $isseq=$ok;
    $infile=$pubseq;
  } else {
    ($ok,$hin)= openRead($pubidtab); 
    $infile=$pubidtab;
  }
	loggit(1,"ERR read_pubids fail from $infile") unless($ok); 
  
  # my %pubids=(); # STILL global .. or read 1st to local then decide..
  %pubids=(); %pubidinfo=();  %puboidc=(); %puballoids=();
  # my %newpubids=();
  
  while(<$hin>) { 
    unless($nred or $isseq or /^\w/) { $isseq=1 if(/^>/); }
    if($isseq) {
      if(/^>(\S+)/) { 
        my $pubid=$1; my($oid)= (m/\boid=([^;,\s]+)/)?$1:"noid";
        my($gid,$alti)= $pubid=~m/(\w+t)(\d+)$/; unless($alti) { $gid||=$pubid; $alti||=0; }
        my($aqual)=(m/aalen=([^;,\s]+)/)?$1:0;
        my($evgclass)=(m/evgclass=([^;,\s]+)/)?$1:0; 
        my($reclass)= ($evgclass =~ /alt|main|cull|noclass/) ? $evgclass
          : ($alti>1)?"alt":"main"; 
        
        #?? split this global %pubids into 2 so no oid/pubid confusion
        $pubids{$oid}= $pubid; # PROBLEMS, need only pubids for user.
        $pubids{$pubid}= $pubid; # PROBLEMS, need only pubids for user.
        $pubidinfo{$pubid}=join("\t", $oid,$gid,$alti,$reclass,$aqual); #  
        
        if(++$nred == 1 and $nalltr == 0) { # reset to existing IDPREFIX
          ($pubid_format,$altid_format,$GBPROID)= make_IDPREFIX($pubid); 
        }
      }
      next;
    }
    unless(/^\w/) { 
      if(not $nred and /^#Public_mRNA_ID/) { $isrealt=1 if(/AltNum\tClass/); }
      next;
    } 
    chomp;  my @v=split"\t"; next unless(@v > 1);
    
    ##  save extended data: reclass, aqual?
    my($pubid,$oid,$gid,$alti,$reclass,$aqual,$pialn,$notes,$alloids)=@v;  # extended for realt_pubids, add Oids col
    $alloids||="";
    map{ $_="" if(not $_ or $_ eq "na" or $_ eq "noid"); } ($oid,$alloids,$notes); # change @v ?
    
    ## need opt to save cull class
    if($reclass =~ /^($skipclass)/) { # $isrealt and .. dont need isrealt set?
      $realtdrops++; # should these go into new publicset/name.culled.files ?? keep from primary publicset files
      next;
    } 

    # my($gpre,$gnum,$ti,$isplit,$gdup)= evgIdParts($oid);
    ## this is prblem xxxt1_G2 diff pubid than xxxt1
    if(CUT_IDCG and $oid=~/_[CG]\d+$/){ (my $oidc=$oid) =~ s/_[CG]\d+$//; $pubids{$oidc}= $pubid;  } # fixme this is map oid,need s/_[CG]\d+// for seqids 
    if(not CUT_IDCG and $oid=~/_[CG]\d+$/){ (my $oidc=$oid) =~ s/_[CG]\d+$//; $puboidc{$oidc}= "$pubid,$oid"; }
    
    my $pisdup=0;
    my $pinfo= join("\t", splice(@v,1)); # all? alloids = (split)[7]

    my($kdxact, $kdxid)=(0,0);
    if($KEEPDROPX) {
      ($kdxact, $kdxid)= check_keepdrop('mapid,seqid',"$oid,$alloids"); # mapid and seqid ??
	    #xverbose: loggit(LOG_DEBUG,"read_pubids keepdropx $pubid $oid,$alloids/$kdxid $kdxact") if($kdxact); 
	    #x # FIXME, use KEEPDROPX/keepdropxtab
      #?? if(nogffDUPIDS and $kdxact =~ /drop|cull/) { $ndropskip++; next; } #?? here keep only keepers
      #OLD: if($kdxact =~ m/drop|cull/) { $ndropskip++; next; } #?? #?? NOT next; for culls.out; here keep only keepers
    }

    if($pubids{$pubid}) {
      $pisdup=1;
      if($kdxid and $pubid eq $kdxid) { $pisdup=1; 
	      # loggit(LOG_DEBUG,"read_pubids keepdropx $pubid dup $kdxact $kdxid/$oid");  # FIXME, use KEEPDROPX/keepdropxtab
        } #?? handle via KEEPDROPX ?
      elsif($pubidinfo{$pubid} eq $pinfo) { $pisdup=9; } # full dup
      elsif($puballoids{$oid}{$pubid}) { $pisdup=2; } # same oid, skip?
      else{ $pisdup=1; }
    }
    
    if($pisdup) {
	    loggit(0,"read_pubids dupid $pubid/$oid dn=$pisdup");  # FIXME, use KEEPDROPX/keepdropxtab
	    if($pisdup>1){ $ndupskip++; next; } else { $ndup1++; }
	    $pinfo= $pubidinfo{$pubid}; # orig info? both?
	    $pinfo.= ",dupid:$oid";
    }
    
    $pubids{$oid}= $pubid; 
    $pubids{$pubid}= $pubid; 
    
    ## buggers, cant put both into pubid key. only oids and pubid=pubid ..
    $pubidinfo{$pubid}= $pinfo; # join("\t", splice(@v,1)); # all? alloids = (split)[7]
    #o $pubidinfo{$pubid}=join("\t", splice(@v,1,5)); # $oid,$gid,$alti,$reclass,$aqual 
    
    ## FIXME: need alloids links?  may have many to many links
    map{
      $puballoids{$_}{$pubid}++;
      if(not CUT_IDCG and $_=~/_[CG]\d+$/){ (my $oidc=$_) =~ s/_[CG]\d+$//; $puballoids{$oidc}{$pubid}++; }
    } grep /\w\w/, split(",","$oid,$alloids"); # alloids contains oid, i think
    
    if(++$nred == 1 and $nalltr == 0) { # reset to existing IDPREFIX
      ($pubid_format,$altid_format,$GBPROID)= make_IDPREFIX($pubid); 
    }
  } close($hin);

  # %pubids=(); # clear global .. or read 1st to local then decide..
  # map{ $pubids{$_}= $newpubids{$_} } keys %newpubids;
  # #$ndupskip,$ndup1,$ndropskip

	loggit(0,"read_pubids n=$nred, dropskip=$ndropskip,dupskip=$ndupskip,dup1=$ndup1, form $pubid_format from $infile"); 
  return($nred, \%pubids);
}


# readKeepDrop globals: our ($KEEPDROPX,%keepdropxtab);  
sub readKeepDrop { # find other subs for this
  my($intable)=@_;
  my($nin,$ndup)= (0) x 9; 
  my(%lkeepdrop); # return to global %keepdrop_idhash
  my($nkdx,%lkeepdropx); # special case dup pubids x oid sets; for global package %keepdropxtab
  my($ok,$hin)= openRead($intable); 
  unless($ok) { loggit(1, "#ERR: missing readKeepDrop($intable)"); return 0; }
  while(<$hin>) {
    next if(/^\W/); chomp;
    my($id,@v)=split; 
    my $act=$v[0];  ## ID here can have _Cn split tag
    my ($oid,$moid)= ($id,$id); # allow dup pub ids, diff oid in col2, act in col3, need dup id fix for keepdrop{$id}
    
    ## upd to for dup pubid table, now use 3 id cols: pubid, seqid, mapid (_C1/2), action, ..
    ## pass this special keepdrop table to evigene_pubsets.pm, make_pubgff, make_pubseq, ..
    unless($act =~ /^(ok|keep|drop|cull)/) {
      my $no=0;
      if($v[2] =~ /^(ok|keep|drop|cull)/) { # and $v[0]=~/^\w+/ no: $v[0]=~/$IDPREFIX/
        $no=2; ($oid)= shift @v; ($moid)= shift @v; $act= $v[0];  
      } elsif($v[1] =~ /^(ok|keep|drop|cull)/) { # and $v[0]=~/^\w+/ no: $v[0]=~/$IDPREFIX/
        $no=1; ($oid)= shift @v; $moid= $oid; $act= $v[0];  
      }
      if($no>0) {
        $nkdx++;  
        $lkeepdropx{act}{$moid}=$act;
        $lkeepdropx{mapid}{$moid}=$id;
        $lkeepdropx{act}{$oid}=$act unless($lkeepdropx{act}{$oid});
        $lkeepdropx{seqid}{$oid}=$id unless($lkeepdropx{seqid}{$oid});
      }
    }
    
    my $kdok= ($act=~/^drop/i) ? kdDROP : ($act=~/^cull/i) ? kdCULL :
             ($act=~/^(ok|keep)/i) ? kdOK : kdOther; # other action verbs?
          
    # NOTE dup id entries allowed, last is active
    # FIXME: need dup pubid / diff oid keep/drop handling, need oid column (3?) defined in keepdrop table
    my $kdval= join"\t",$kdok,$oid,@v; # @v has oid now if id,oid,act,.. table
    my $dupskip=0;
    if(my $oldv= $lkeepdrop{$id}){ # warn?
      $dupskip=1 if($kdok ne kdOK and $oldv =~ m/\b(ok|keep)/);
      $ndup++; } 
    $nin++;
    $lkeepdrop{$id} = $kdval unless($dupskip); 
    $lkeepdrop{$oid}= $kdval if($oid ne $id and not $lkeepdrop{$oid}); #? skip for keepdropxtab
  } close($hin);
  
  if($nkdx>0) {
    $KEEPDROPX= $nkdx;  %keepdropxtab= %lkeepdropx; #?? bugs?
    #   for my $s (qw(act mapid seqid)) {
    #     for my $k (keys %{$lkeepdropx{$s}}) { $keepdropxtab{$s}{$k}= $lkeepdropx{$s}{$k}; }
    #   }
  } else {
    $KEEPDROPX=0; %keepdropxtab= (); 
  }
  loggit(0, "readKeepDrop($intable), n=$nin, dupentries=$ndup, keepdropx=$nkdx");
  return($nin,\%lkeepdrop);
}

sub check_keepdrop {
  my($forseqormap,$alloids)= @_;
  if($KEEPDROPX) {
    #   $keepdropxtab{act}{$moid}=$act;  # use act =~ m/(keep|ok|drop|cull)/ to filter dup pubids
    #   $keepdropxtab{mapid}{$moid}=$pubid;
    $forseqormap='mapid' unless($forseqormap and $forseqormap =~ /^(seqid|mapid)/);
    my @fsm= split",", $forseqormap;
    for my $fsm (@fsm) {
    for my $od (split",",$alloids) { 
      my $xact= $keepdropxtab{act}{$od} or next;
      my $xid = $keepdropxtab{$fsm}{$od}; # if xact but not xid, check both seqid,mapid?
      if($xid and $xact) { return( $xact, $xid); } 
    }
    }
  }
  return(0,0);
}



sub read_annotab
{
  # my($annotab,$savePubids)=@_;
  my($nan,$anh)= read_annotab2(@_);
  return %$anh; 
}

sub read_annotab2
{
  my($annotab,$savePubids)=@_;
  # uses globals: %pubids,%pubidinfo : change?
  # OPT to fill annothash from %pubids/info if no annotab
  
  my %annothash=(); ## make this global? use same for 2+ subs
  my($ok,$annoth)= ($annotab and -f $annotab) ? openRead($annotab):(0,0);
  my($nan,$addmapcols)= (0) x 9;
	loggit(0,"read_annotab($annotab): ok=$ok"); 
  unless($ok) { loggit(1, "#ERR: missing read_annotab($annotab)"); }# not return, add pubids

	if($ok) {
		while(<$annoth>) { 
		  next unless(/^\w/); chomp; 
			my @cols = split"\t";  
			my($pubid,$oid)= @cols;
			  # PublicID	OrigID	TrLen	CDSoff	AAqual	TrGaps	Dbxref	Namealign	Product_Name	CDD_Name
			  # add EvgClass:10 after CDD:9
			  ##? maybe change annotab format to kvtab, ncbiannot=value, allow special fields like note=xxxx

			##FIXME1809: ANNO_ADD_MAPCOLS for @ANNO_COLS  
			if(ANNO_ADD_MAPCOLS) {
			  if($pubid eq "PublicID") {
			    if(@cols + 2 == @ANNO_COLS) {
			      $addmapcols=1;
			      push @cols, qw(Maploc Mapqual);
			      s/$/\tMaploc\tMapqual/;
			    }
			  } elsif($addmapcols) {
          my($oid,$gid,$alti,$class,$aqual,$pia,$notes,$alloids)= split"\t", $pubidinfo{$pubid};
          # notes == chrmap:73a,99i,22691l,73x,chr24:36045828-36126968:+ 
          my($mapqual,$maploc)= ($notes) ? mapinfostr($notes) : (0,0);
			    push @cols, $maploc, $mapqual;
			    s/$/\t$maploc\t$mapqual/;
			  }
			}

			# OPT to fill in %pubids from this table
			#? is annot tab bad? should save pubids, use input mrna instead??
			if($savePubids) {
			  my($gid,$alti)= $pubid=~m/(\w+t)(\d+)$/;
        $pubids{$oid}= $pubid; 
        $pubids{$pubid}= $pubid; 
        my $aqual= $cols[4]; 
        my $class= (@cols>9)? $cols[10] : "";
        unless($class and $class=~/alt|main|noclass|cull/) { $class=($alti>1)?"alt":"main"; } # dont have here
        $pubidinfo{$pubid}=join("\t", $oid,$gid,$alti,$class,$aqual); #  
        }
			$annothash{$oid}= $annothash{$pubid}= $_;  #? both? mrnaseq in should be >oid; might be >pubid
			
			## _[CG] is problem, nnn_G2 has diff pubid than nnn
      if(CUT_IDCG and $oid=~/_[CG]\d+$/){ (my $soid=$oid) =~ s/_[CG]\d+$//; $annothash{$soid}= $_;  } # fixme this is map oid,need s/_[CG]\d+// for seqids 
			$nan++;
		}	close($annoth);
	}
	
  my $npd= scalar(keys %pubids);
  my $addpubinfo= (not $ok or $nan < $npd);
  if($addpubinfo) { # do this anyway to fill in missed annoth ?
    # pubid: my($pubid,$oid,$gid,$alti,$reclass,$aqual,$pialn,$notes)=@v;  # extended for realt_pubids
    # $pubidinfo{$pubid}=join("\t", $oid,$gid,$alti,$reclass,$aqual); #  
    # annline=		($pubid,$oid,$trlen,$cdsoff,$aaqual,$annogaps,$dbxref,$namepct,$gname,$cddname)
    my ($trlen,$coff)=(0,0); # trlen from elsewhere?
    for my $pubid (sort keys %pubidinfo) {
      next if($annothash{$pubid});
      #x my($oid,$gid,$alti,$class,$aqual)= split"\t", $pubidinfo{$pubid};
      my($oid,$gid,$alti,$class,$aqual,$pia,$notes,$alloids)= split"\t", $pubidinfo{$pubid};
      # ANNO_ADD_MAPCOLS fixme here..
      my $annline=join"\t", $pubid, $oid, $trlen, $coff, $aqual, 0, 0, 0, 0, 0;
      $annothash{$oid}= $annothash{$pubid}= $annline; $nan++;
        # fixme this is map oid,need s/_[CG]\d+// for seqids 
      if(CUT_IDCG and $oid=~/_[CG]\d+$/){ (my $soid=$oid) =~ s/_[CG]\d+$//; $annothash{$soid}= $annline;  }       
    }
  }
	
  return($nan,\%annothash);
  # return %annothash; # or \%hash ?
}

# sub annotblinfo {
#   my($oid,$hd,$fa,$pubid,$aahdr)= @_; 
#   # our(%aahdr);
#   # my $aahdr= $aahdr{$oid};
#   
#   my $lenfa= length($fa);
#   my $tblinfo= parse_evgheader($oid,$hd,$lenfa); # uses %genenames ...
#   if($pubid){ $tblinfo->{'pubid'}= $pubid; }
#   
#   if(($oid =~ /utrorf/ or $hd=~/uvcut=/) and $aahdr) {
#     if($hd=~/uvcut=([^;\s]+)/){ $tblinfo->{'uvcut'}=$1; }
#     my $aainfo= parse_evgheader($oid,$aahdr,$lenfa);
#     ## DANG, this is bad for orig ok.aa set (rev.tr), only do this for utrorf.aa  AND? trim/uvcut.aa  
#     ## .. bad case for strand=- // cdsor=; check aaqual=\d+ is same on both?
#     if($aainfo->{'cdsor'} eq '+') { 
#       map{ $tblinfo->{$_}= $aainfo->{$_}; } qw(aaqual cdsoff cdsor);
#     }
#   }
#   
#   my $trgaps= $fa =~ tr/Nn/Nn/; # for gff hd, gaps= ??
#   $tblinfo->{trgaps}= $trgaps; #? 
#   
#   return $tblinfo;
# }  

sub putPubAnnot {
  my($outh,$oid,$hd,$tblinfo,$didseq,$itr)= @_;
  return 0 unless($oid and $tblinfo);
  
  (my $soid=$oid)=~ s/_[CG]\d+$//; # use this or use aoids ?
  my $pid= $pubids{$oid} || $pubids{$soid}; # problem here is pubids.oid may have _G/C, not seq oid
  my($aoids)= $hd =~ m/oid=([^;\s]+)/?$1:"";
  
  my($kdxact, $kdxid)=(0,0);
  if($KEEPDROPX) {
    ($kdxact, $kdxid)= check_keepdrop('mapid',"$oid,$aoids"); # mapid,seqid ?
    if(nogffDUPIDS and $kdxact =~ /drop/) { return 0 ; } #??
  }

  unless($pid) { 
    ## add $puballoids{$oid}{$pubid} check?
    # my($aoids)= $hd =~ m/oid=([^;\s]+)/?$1:"";
    for my $d (split ",", "$oid,$aoids") {
      my($p,$od)= split",", ($puboidc{$d}||""); if($p){ $pid=$p; $oid=$od; last; } 
    }
  }  
  $pid=0 if($pid and $didseq->{$pid});
  if($pid) {  # ignore:  or not $skipdropseqs
    # my $tblinfo= annotblinfo($oid,$hd,$fa,$pid); #  $aahdr{$oid};
    my $oidout=$oid;
    $tblinfo->{'pubid'}= $pid; 
    if($oid =~ m/^$pid/) {  # dont stutter like this in Oid column
      my($oidlast)=  grep{ not m/$oid/ } split ",", $aoids;
      $oidlast ||= $oid; $oidout=$oidlast;
      $tblinfo->{'oid'}= $oidlast;
    }
    ## fixme pid may have diff class from tblinfo (which may hve no evgclass)
    my($xxoid,$xxgid,$xxalti,$evgclass,$xxaqual,$xxpia,$notes,$alloids)= split"\t", $pubidinfo{$pid};
    if($kdxact =~ /cull|drop/) { $evgclass=$kdxact.$evgclass unless($evgclass =~ /^(cull|drop)/); }
    $tblinfo->{'evgclass'}= $evgclass if($evgclass); 
    
    # if(ANNO_ADD_MAPCOLS) {}      
    # if(MAPQ1805)  chrmap:99a,99i,1234l,9x,chr9:123-456:+
    if( $notes and my($mapq)= $notes=~m/chrmap[:=]([^;\s]+)/ ) {
      ## $mapq =~ s/,(pflag|tscore).*//;
      my @mq= split",", $mapq;
      my ($mloc)= grep/:/, @mq; # $mq[4] maybe
      if($mloc) { $mapq=~s/,$mloc.*//; }
      $tblinfo->{'mapqual'}= $mapq;
      $tblinfo->{'maploc'}= $mloc;
    }
      
    my($no)= putAnnot($outh,$oidout,$tblinfo,$itr); 
    $didseq->{$oid}++; $didseq->{$oidout}++; $didseq->{$pid}++; # add back?
    return ($oid,$pid); # soid ?
  } else {
    return 0;
  }
}
       	    
sub make_annot_from  # @ingff or @inmrna
{
  my($fromtype, $annotab, $inputlist, $genenames, )=@_;

  my($outh, %didseq); %didseq=();
  my($inh,$ok,$notr,$nocds)= (0) x 9;

  return($annotab,$notr,$nocds) if( -s $annotab ); # or $skiprun

  my($in1)= (ref($inputlist) =~ /ARRAY/) ? $inputlist->[0] : $inputlist;
  $ok = $in1 and -f $in1;
  
  # if($ok and -f $annotab) { rename($annotab,"$annotab.old"); }
  if($ok) { $ok= open($outh,'>',$annotab); }
  unless($ok) { loggit(1,"ERR: make_annot_from($fromtype, $in1) TO $annotab"); return; }

  $genenames||="";
  my($nnamed,$namin)= parse_genenames($genenames,NAME_NOEDITS);
  loggit(0, "names $genenames n=$nnamed\n"); 

## gff Split gene stutter now: have 3 annot rows for 1 split gene;
## add all split oid_C1,2 variants to didseq?
# grep vDanrer6pEVm000009t18 chreset-map6na.ann.txt | head
# vDanrer6pEVm000009t18	Danrer6pEVm000009t16_C1	22154	170-17215	5681,76p,complete	0	zfish16nc:XP_017209375.1,	64	Obscurin	na	alt
# vDanrer6pEVm000009t18	Danrer6pEVm000009t16_C2	22154	170-17215	5681,76p,complete	0	zfish16nc:XP_017209375.1,	64	Obscurin	na	alt
# vDanrer6pEVm000009t18	Danrer6pEVm000009t16	0	0	na	0	na	0	Uncharacterized protein	na	alt
 
  if($fromtype =~ /gff/) {
    ($notr)= make_annot_from_gff($outh,$inputlist,\%didseq);
  } else {
    ($notr)= make_annot_from_seq($outh,$inputlist,\%didseq);
  }
  
  ## check  $didseq{$oid} vs pubids{} 
  ## BUG: %pubids has both pid and oid, @undone has dups
  my @undone= grep { my $p=$pubids{$_};  not( m/^$IDPREFIX/ or $didseq{$p}) } keys %pubids;
  if(@undone) {
    my $itr=$notr;
 	 	foreach my $oid (sort @undone) {  
 	 	  my $tblinfo= annotab2tblinfo($oid,"",0,0,0);
 	 	  my($no)= putAnnot($outh,$oid,$tblinfo,++$itr); 
	 	  $nocds+= $no; $didseq{$oid}++; ## NOT: $notr+= $otr; 
		}
  }
  
  close($outh);
  push @publicset, $annotab; # push @tmpfiles, $inseq; #??
  
  #add summary output?
  # my ($table_g1,@moretabs)= geneset_summary(\%pubids, $annotab);
  # if($ok) { $ok= open($outh,'>',$genesum); print $outh $table_g1; close($outh); }
  
  return($notr,$nocds); 
}

sub make_annot_from_gff  # @ingff or @inmrna
{
  my($outh, $inputlist, $didseq)=@_;
  
  my @inf;
  if($inputlist and ref($inputlist) =~ /ARRAY/){ @inf= @$inputlist; }
  else { @inf= ($inputlist); }
  return(0,"missingdata") unless( $inputlist and  -f $inf[0] );

  sub gff_tblinfo {
    my($oid,$hd,@gff)= @_; 
    #  unless(splitKEEP_C){ $oid=~ (s/_C(\d+)//; } # want this?
    (my $seqoid=$oid) =~ s/_[CG]\d+$//; # seq oid for names
    my $tblinfo= parse_evgheader($oid,$hd,0,$seqoid); # uses %genenames ... FIXME alloids check of genenames
    # if($pubid){ $tblinfo->{'pubid'}= $pubid; }
    # if(($oid =~ /utrorf/ or $hd=~/uvcut=/) and $aahdr) {}
    
    my($trgaps)= ($hd =~ m/;gaps=(\d+)/)?$1:0; #? maybe 
    $tblinfo->{trgaps}= $trgaps; #? 
    return $tblinfo;
  }  

## gff Split gene stutter now: have 3 annot rows for 1 split gene;
## add all split oid_C1,2 variants to didseq?
# grep vDanrer6pEVm000009t18 chreset-map6na.ann.txt | head
# vDanrer6pEVm000009t18	Danrer6pEVm000009t16_C1	22154	170-17215	5681,76p,complete	0	zfish16nc:XP_017209375.1,	64	Obscurin	na	alt
# vDanrer6pEVm000009t18	Danrer6pEVm000009t16_C2	22154	170-17215	5681,76p,complete	0	zfish16nc:XP_017209375.1,	64	Obscurin	na	alt
# vDanrer6pEVm000009t18	Danrer6pEVm000009t16	0	0	na	0	na	0	Uncharacterized protein	na	alt
  
  my($ngenes,$itr)= (0,0);
  for my $inf (@inf) {
    my($ok,$inh)= openRead($inf);
    unless($ok) { loggit(1,"err reading $inf"); next; }
    my($oid,$pid,$inid,$gann,$gsrc,$tblinfo,@gene);
    while (<$inh>) {
      if(/^\W/){ next; }  
      elsif(/\t($RNATYPES)/) { 
        if($inid and $gann and not $didseq->{$inid}) {
        # FIXME: is always inid eq oid ?
        
        # new bug for ingff == last pubgff, inid == pubid, dont call this new oid .. need option?
        
        $tblinfo= gff_tblinfo($inid,$gann,@gene);
        ($oid,$pid)= putPubAnnot($outh,$inid,$gann,$tblinfo,$didseq,++$itr); # ,$gsrc,@gene
        $ngenes++ if($oid);
        #x if($oid){ $didseq->{$oid}++; $didseq->{$inid}= $didseq->{$pid}= $didseq->{$oid};  $ngenes++; }
        }
        
        $inid=$gann=0;
        chomp; @gene=split"\t";
        $gsrc= $gene[1]; #?  $outsource || use input gff[1] col unless defined? not used here now?
        $gann= $gene[8];
        my($id)= ($gann =~ m/\bID=([^;\s]+)/)?$1:0;
        $inid= $id;        
        #  unless(splitKEEP_C){ $id=~ (s/_C(\d+)//; } # want this?
      }
    }
    
    if($inid and $gann and not $didseq->{$inid}) {
      $tblinfo= gff_tblinfo($inid,$gann,@gene);
      ($oid,$pid)= putPubAnnot($outh,$inid,$gann,$tblinfo,$didseq,++$itr); # ,$gsrc,@gene
      $ngenes++ if($oid);
      #x if($oid){ $didseq->{$oid}++; $didseq->{$inid}= $didseq->{$pid}= $didseq->{$oid};  $ngenes++; }
    }
    close($inh);
  } 

  return($ngenes); 
}

sub make_annot_from_seq  # @ingff or @inmrna
{
  my($outh, $inputlist, $didseq)=@_;

  my @inf;
  if($inputlist and ref($inputlist) =~ /ARRAY/){ @inf= @$inputlist; }
  else { @inf= ($inputlist); }
  return(0,"missingdata") unless( $inputlist and  -f $inf[0] );

  sub seq_tblinfo {
    my($oid,$hd,$fa,$aahdr)= @_; 
    my $lenfa= length($fa);
    #not for seq info: (my $seqoid=$oid) =~ s/_[CG]\d+$//; 
    my $tblinfo= parse_evgheader($oid,$hd,$lenfa); # uses %genenames ...
    
    if(($oid =~ /utrorf/ or $hd=~/uvcut=/) and $aahdr) {
      if($hd=~/uvcut=([^;\s]+)/){ $tblinfo->{'uvcut'}=$1; }
      my $aainfo= parse_evgheader($oid,$aahdr,$lenfa);
      ## DANG, this is bad for orig ok.aa set (rev.tr), only do this for utrorf.aa  AND? trim/uvcut.aa  
      ## .. bad case for strand=- // cdsor=; check aaqual=\d+ is same on both?
      if($aainfo->{'cdsor'} eq '+') { 
        map{ $tblinfo->{$_}= $aainfo->{$_}; } qw(aaqual cdsoff cdsor);
      }
    }
    my $trgaps= $fa =~ tr/Nn/Nn/; # for gff hd, gaps= ??
    $tblinfo->{trgaps}= $trgaps; #? 
    return $tblinfo;
  }  

  
  my($ngenes,$itr)= (0,0);
  for my $inf (@inf) {
    my($ok,$inh)= openRead($inf);
    unless($ok) { loggit(1,"err reading $inf"); next; }
    my($oid,$pid,$inid,$gann,$gsrc,$tblinfo,$fa);

    my %aahdr=(); 
    my $aaseq  =  makename($inf,".aa");  
    $aaseq="$aaseq.gz" if(! -f $aaseq and -f "$aaseq.gz");
    if(-f $aaseq) {  
      my($aok,$aah)= openRead($aaseq);
      while(<$aah>) { if(/^>(\S+)\s+(.+)$/) { $aahdr{$1}= $2; } } 
      close($aah);
    }

    while(<$inh>) { 
      if(/^>(\S+)/) { my $d=$1; 
        if($inid and $gann and not $didseq->{$inid}) {
          $tblinfo= seq_tblinfo($inid,$gann,$fa,$aahdr{$inid});
          ($oid,$pid)= putPubAnnot($outh,$inid,$gann,$tblinfo,$didseq,++$itr); # ,$gsrc,@gene
          $ngenes++ if($oid); 
          #x if($oid){ $didseq->{$oid}++; $didseq->{$inid}= $didseq->{$pid}= $didseq->{$oid};  $ngenes++; }
          }
        $inid=$d; $gann=$_; chomp($gann); 
        $fa=""; # $itr++;
        }
      elsif(/^\w/) { chomp; $fa.=$_; }
    } 
    close($inh); 
  
    if($inid and $gann and not $didseq->{$inid}) {
      $tblinfo= seq_tblinfo($inid,$gann,$fa,$aahdr{$inid});
      ($oid,$pid)= putPubAnnot($outh,$inid,$gann,$tblinfo,$didseq,++$itr); # ,$gsrc,@gene
      $ngenes++ if($oid); 
      #x if($oid){ $didseq->{$oid}++; $didseq->{$inid}= $didseq->{$pid}= $didseq->{$oid};  $ngenes++; }
      }
  
  }
  
  return($ngenes); 
}
  

sub make_annotab # or make_publicset
{
  my($mrnaseq,$genenames,$skiprun, $outname)=@_;
  $outname ||= $mrnaseq;
  # globals: %pubids, %genenames,.. 
  
  #FIXME* input genes.gff instead of mrnaseq for annot info; parse_evgheader works on both?
	#FIXME: also re-write publicset/mrna,aa,cds w/ >pubid  old=oldid; name=xxx, other annots in header ...
  
  my($notr,$nocds)=(0,0);
	my $annotab =  makename($outname,".ann.txt"); ## FIXME: change suffix: was .annotab
	my $pubdir="publicset";
  if(not -f $annotab and $pubdir and -d $pubdir) {
  	 my($pubd,$ft)= getFileset($pubdir,'.ann.txt');  $annotab=$ft if($ft);
  }
  return($annotab,$notr,$nocds) if( -s $annotab or $skiprun); # or dryrun ..
  return("",0,0) unless( $mrnaseq and -f $mrnaseq); # alternate annot data, gene.gff ?

  my($nnamed,$namin)= parse_genenames($genenames,NAME_NOEDITS);
  loggit(0, "names $genenames n=$nnamed\n"); 

  our($outh, %didseq); %didseq=();
  my ($inh,$ok,$hd,$oid,$fa)= (0) x 9;
  my($itr,$otr,$ocds,$oerr)= (0) x 10;

	($ok,$inh)= openRead($mrnaseq);
  $ok= open($outh,'>',$annotab) if($ok);
  unless($ok) { loggit(1,"ERR: make_annotab $mrnaseq TO $annotab"); return; }

  ## BUG 201405  here missed utrorf aaqual, from missing input aaseq next to mrnaseq .. makename() suffix bug?
  #FOXME: annotab had WRONG aalen/qual, from mrna but not updated uvcut.aa or utrorf.aa
  # .. need to parse pubset/*.aa not or as well as .mrna

  our %aahdr=(); 
  my $aaseq  =  makename($mrnaseq,".aa"); # aa.gz ?
  $aaseq="$aaseq.gz" if(! -f $aaseq and -f "$aaseq.gz");
  if(-f $aaseq) { # warn if missing ??
	  my($aok,$aah)= openRead($aaseq);
    while(<$aah>) { if(/^>(\S+)\s+(.+)$/) { $aahdr{$1}= $2; } } close($aah);
  } else {
    loggit(1,"make_annotab missing $aaseq");
  }

  sub maketblinfo {
    my($oid,$hd,$fa,$pubid)= @_; our(%aahdr);
    my $lenfa= length($fa);
    my $tblinfo= parse_evgheader($oid,$hd,$lenfa); # uses %genenames ...
    if($pubid){ $tblinfo->{'pubid'}= $pubid; }
    if(($oid =~ /utrorf/ or $hd=~/uvcut=/) and (my $aahdr= $aahdr{$oid})) {
      if($hd=~/uvcut=([^;\s]+)/){ $tblinfo->{'uvcut'}=$1; }
      my $aainfo= parse_evgheader($oid,$aahdr,$lenfa);
      ## DANG, this is bad for orig ok.aa set (rev.tr), only do this for utrorf.aa  AND? trim/uvcut.aa  
      ## .. bad case for strand=- // cdsor=; check aaqual=\d+ is same on both?
      if($aainfo->{'cdsor'} eq '+') { 
        map{ $tblinfo->{$_}= $aainfo->{$_}; } qw(aaqual cdsoff cdsor);
      }
    }
    my $trgaps= $fa =~ tr/Nn/Nn/; $tblinfo->{trgaps}= $trgaps; #? 
    return $tblinfo;
  }  
  
  sub putpubann {
    my($oid,$hd,$fa,$itr)= @_;
	  our(%didseq,$outh); # no dup oids
	  return 0 unless($oid and $hd);
    my $pid= $pubids{$oid}; # problem here is pubids.oid may have _G/C, not seq oid
    unless($pid){ 
      ## add $puballoids{$oid}{$pubid} check?
      my($aoids)= $hd =~ m/oid=([^;\s]+)/?$1:"";
      for my $d (split ",", "$oid,$aoids") {
        my($p,$od)= split",", ($puboidc{$d}||""); if($p){ $pid=$p; $oid=$od; last; } 
      }
    }  
    if($pid or not $skipdropseqs) {
      ## return 0 if($didseq{$oid}); #?? any dups
      my $tblinfo= maketblinfo($oid,$hd,$fa,$pid);
      my($no)= putAnnot($outh,$oid,$tblinfo,$itr); 
      $didseq{$oid}++;
      return $no;
    } else {
      return 0;
    }
  }
       	    
  
  while(<$inh>) { 
    if(/^>(\S+)/) { my $d=$1; 
    	#x if($oid and ($pubids{$oid} or not $skipdropseqs))
      if($oid) { $notr+= putpubann($oid,$hd,$fa,$itr); }
      $oid=$d; $hd=$_; chomp($hd); 
      $fa=""; $itr++;
      }
 		elsif(/^\w/) { chomp; $fa.=$_; }
  } 
  
	# if($oid and ($pubids{$oid} or not $skipdropseqs))
  if($oid) { $notr+= putpubann($oid,$hd,$fa,$itr); }
  close($inh); 
  
	## FIXME for utrorf.aa also in pubids{}, but not in mrnaseq ...
  ## check  $didseq{$oid} vs pubids{} 
  my @undone= grep { not( m/^$IDPREFIX/ or $didseq{$_}) } keys %pubids;
  if(@undone) {
		# my $pubaa =  makename($mrnaseq,".aa"); 
  	# if( ! -f $pubaa and -f "$pubaa.gz" ) { $pubaa.=".gz"; }
 		# ($ok,$inh)= openRead($pubaa); .. parse @undone aa and putAnnot ??
 	 	foreach $oid (sort @undone) {  
 	 	  my $tblinfo= annotab2tblinfo($oid,"",0,0,0);
 	 	  my($no)= putAnnot($outh,$oid,$tblinfo,++$itr); 
	 	  $nocds+= $no; $didseq{$oid}++; ## NOT: $notr+= $otr; 
		}
  }
  
  close($outh);
  push @publicset, $annotab; # push @tmpfiles, $inseq; #??
  return($annotab,$notr,$nocds); 
}

sub putAnnot {
  my ($annoth, $oid, $tblinfo, $itr)=@_;  ## $hdr, 

	# FIXME: annorec becomes tblinfo, not from faseq.hdr parsing here
  # my $tblinfo= parse_evgheader($oid,$hdr,length($faseq));
  # upd1803: evgclass ; add Class col, esp cull vs ok
  # my($cdsoff,$gname,$dbxref,$aaqual,$trlen,$trgaps,$protid,$lotag,$namepct, $cddname)= 
  #    @{$tblinfo}{qw(cdsoff name dbxref aaqual trlen trgaps protid locustag namepct cdd)};
  
  my $pubid= $tblinfo->{'pubid'} || $oid; ## or what? err?
  my($cdsoff,$gname,$dbxref,$aaqual,$trlen,$trgaps,$protid,$lotag,$namepct, $cddname, $evgclass)= 
      map{ $tblinfo->{$_} || 0 } 
      qw(cdsoff name dbxref aaqual trlen trgaps protid locustag namepct cdd evgclass);
  # my @tblval= map{ $tblinfo->{$_} || 0 } @ANNO_TBLINFO;
      
  # @ANNO_COLS= qw(PublicID OrigID TrLen CDSoff AAqual TrGaps Dbxref Namealign Product_Name CDD_Name Class Maploc Mapqual);
  # @ANNO_TBLINFO = qw(pubid oid trlen cdsoff aaqual trgaps dbxref namepct name cdd evgclass maploc mapqual);

  ($gname,$namepct)= productname($gname,$aaqual,$namepct);
  
  my @endvals=($evgclass); my @endc=qw(Class);
  # UPD add mapqual/locus cols to ann.txt, from tblinfo keys
  if(ANNO_ADD_MAPCOLS) {
    push @endc, qw(Maploc Mapqual);
    push @endvals, map{ $tblinfo->{$_} || 0 } qw(maploc mapqual);
  }
    
  # pubid cols: Public_mRNA_ID	originalID	PublicGeneID	AltNum	Class	AAqual	pIdAln	Notes	Oids
  # same: PublicID=Public_mRNA_ID; OrigID=originalID; Class=Class; AAqual=AAqual; 
	#x.print $annoth join("\t",qw(PublicID OrigID TrLen CDSoff AAqual TrGaps Dbxref Namealign Product_Name CDD_Name Class) )."\n" 

  if($itr==1) {
    my @cols= qw( PublicID OrigID TrLen CDSoff AAqual TrGaps Dbxref Namealign Product_Name CDD_Name );
    push @cols, @endc; @ANNO_COLS= @cols; # FIXME !**!*!
	  print $annoth join("\t", @cols)."\n";
	}
	print $annoth join("\t",$pubid,$oid,$trlen,$cdsoff,$aaqual,$trgaps,$dbxref,$namepct,$gname,$cddname, @endvals)."\n";
 	return(1);
}


sub putPubSeq # for mrna,aa,cds.public w/ rewritten headers????
{
  my ($outh, $stype, $oid, $hdr, $faseq, $tblinfo, $pubidin, $idupout)=@_; 
	## stype = mRNA|cDNA|ncRNA, CDS, amino|protein|polypeptide
  my($ntrout,$ncdsout,$nerr)=(0) x 10;
  my $pubid= $tblinfo->{'pubid'} || $pubidin || $oid; ## or what? err?
  # OPTION : return/drop any seq lacking pubid .. ie source okayset has extras..
  $stype= $tblinfo->{'seqtype'} || $stype;
  
  ## add   $tblinfo{'evgclass'}= $evgclass if($evgclass);
  ## replace tr2aacds evgclass w/ pubids evgclass
  ## <evgclass=althi,okay,match:tridba1a_sNn12l1SRR1524240ridk81Loc0,pct:100/68/.;
  ## >evgclass=cullalt
  
  my($cdsoff,$gname,$dbxref,$aaqual,$trlen,$protid,$lotag,$namepct, $cddname, $evgclass)= 
      @{$tblinfo}{qw(cdsoff name dbxref aaqual trlen protid locustag namepct cdd evgclass)};
  map{ $_='' if($_ eq "na"); } ($protid,$lotag,$gname,$cddname,$dbxref,$aaqual);  
	$gname='' if($gname =~ /^($NAME_NONE)/i);

  ## preserve input oid, add new if needed: |oid
  my $DROPxtrakeys= 'i|val|flag|path|strand|species';
	my $DROPevgkeys = 'type|Name|Dbxref|organism|db_xref|aalen|aaqual|aaSize|clen|offs|cdsoff';  # evigene hdr keys, will replace
  # evigene set: do below# aalen|aaqual|aaSize|clen|offs|strand';
	
 	$hdr =~ s/>\S+//; 
 	#?? change to keep keys when FAHDRisValid  and no replacement ?
 	my $keepkeys= ($FAHDRisValid and not ($trlen and $cdsoff))?1:0;
 	$hdr =~ s/\s*\b($DROPxtrakeys)=[^;\n]+;?//g;
 	
 	## should step thru keys and replace 1b1 if FAHDRisValid..
	our $hdr2="";
 	if(1) {
	  our $hdr1= $hdr;
 	  sub updkey{ my($k,$val)=@_; our($hdr1,$hdr2); if($val) { $hdr2.="$k=$val; "; 
 	    if($hdr1=~m/\b$k=/) { $hdr1 =~ s/\s*\b$k=[^;\n]+;?//; } } }
 	  updkey('type',$stype);
 	  updkey('Name',$gname);
 	  updkey('Dbxref',$dbxref); # add namepct ? cdd?
 	  updkey('aalen',$aaqual);
 	  updkey('clen',$trlen);
 	  updkey('offs',$cdsoff);
 	  updkey('organism',$ORGANISM) if(okORGANISM());
    updkey('evgclass',$evgclass);
	  updkey('copy',1+$idupout) if($idupout>0); # dupl=2 seqdup=2 copy=2 ?
 	  $hdr= $hdr1;
 	}
  # else {
  #   unless($keepkeys) { $hdr =~ s/\s*\b($DROPevgkeys)=[^;\n]+;?//g; }
  #   $hdr2 .="type=$stype; ";
  #   $hdr2 .="Name=$gname; " if($gname);
  #   $hdr2 .="Dbxref=$dbxref; " if($dbxref);
  #   $hdr2 .="aalen=$aaqual; " if($aaqual);
  #   $hdr2 .="clen=$trlen; " if($trlen);
  #   $hdr2 .="offs=$cdsoff; " if($cdsoff);
  #   # $hdr2 .="oid=$oid; " if($oid);
  #   $hdr2 .="organism=$ORGANISM; " if(okORGANISM()); #?? and $organism ne $DEFAULT_SETTINGS{'organism'});
  # }
 	
 	if(1) { # clean/uniq oids
 	  my($oida)= $hdr =~ m/\boid=([^;\s]+)/?$1:"";
 	  my @oh= split",",$oida;
 	  unshift(@oh,$oid) if($oid and $oid ne $pubid);
 	  my %oh= map{ $_ => 1 } @oh;  # keep order?
 	  my $oh= join",", grep { 1 == $oh{$_}++ } @oh; # keep in order
 	  $hdr2.="oid=$oh; "; $hdr =~ s/\boid=([^;\s]+)//;
 	} else {
	  if($oid and $oid ne $pubid) {
	    if($hdr =~ m/\boid=/){ $hdr=~s/oid=/oid=$oid,/; } else { $hdr2.="oid=$oid; " }
	  }
	}
	
	if(AA_NOSTOPCODON and $stype =~ /protein/ and substr($faseq,-1) eq '*') {
	  chop($faseq); ## or $faseq= substr($faseq,0,-1);
	}
	
	$hdr ="" unless($hdr =~ /\w/); # or drop all?
	$hdr2 =~ s/;\s*$//; $hdr=~s/^\W+/ /;
  $faseq =~ s/(.{60})/$1\n/g; 
  print $outh ">$pubid $hdr2;$hdr\n$faseq\n";  $ntrout++;
  return($pubid, $ntrout ); # ,$ncdsout,$nerr,
}

=item make_pubgff
  sub make_pubgff(\@inputgenegff,$annothash) { }
  read ingff, filter by pubids(pd,oid), annotate, write outgff
  see evigene/scripts/genes/evganngff.pl for annot handling
  see evigene/scripts/genes/evigenegff.pm : in progress.. use that here? OR better move this sub to there.

=cut


sub make_pubgff {
  my($inputgenegff,$outsource,$annothash,$outname)=@_;
  
  ## FIXME: make_pubgff should change some attr as for make_pubseq, eg. Name updates
  ## .. rely on annot.tab ? which probably was created from input.gff, but may be updated
  
  my @ingff=();  
  if($inputgenegff and ref($inputgenegff) =~ /ARRAY/){ @ingff= @$inputgenegff; }
  elsif($inputgenegff) { @ingff= ($inputgenegff); }
  return("",0,"missingdata") unless( @ingff > 0 and -f $ingff[0] );

	my $pubdir="publicset"; my $pubd=undef; my $ft;
	my $outgff= $outname || $ingff[0]; # array or 1 file?
	$outgff =~ s/.gz$//; $outgff=~s/\.gff.*//; 
  ## add _cull.gff as per pubseq
  my $cullgff= ($CULLSEP) ? $outgff."_cull.gff" : "";
  $outgff.="_pub.gff"; 
  my($pubsuf)= "_pub.gff"; # $outgff=~m/(\.\w+.gff)$/;
  unless(-s $outgff or not -d $pubdir) { ($pubd,$ft)= getFileset($pubdir,$pubsuf,$pubd);  $outgff=$ft if($ft); }
  return($outgff,0,"exists") if( -s $outgff); # or dryrun ..

  our($outh,$cullh,$annoth,%didseq);  %didseq=();
  my($ok,$ngenes,$nexons,$nsplit)=(0) x 9;
  # our links to sub 
  $annoth = $annothash;

  sub putpubgff {
    my($oid,$mattr,$source,@gff)= @_;
	  our(%didseq,$annoth,$outh,$cullh); # no dup oids
	  return 0 unless($oid and $mattr);
    my $pid= $pubids{$oid}; # problem here is pubids.oid may have _G/C, not seq oid

    my($kdxact, $kdxid)=(0,0);
    if($KEEPDROPX) {
      my($xod)= ($mattr =~ m/;oid=([^;\s]+)/)?$1:""; 
      ($kdxact, $kdxid)= check_keepdrop('mapid',"$oid,$xod");
      if(nogffDUPIDS and $kdxact =~ /drop/) { return 0 ; }
    }
        
    # FIXME splits _C, cut other way here, as gff2pubidtab cut them out: KEEP_C => 0;
    #  unless(KEEP_C){ $id=~s/_C\d// if($isplit); $mainid=~s/_C\d$//; } # want this?
    unless($pid) { 
      ## add $puballoids{$oid}{$pubid} check?
      (my $oidc= $oid) =~ s/_C\d//; $pid= $pubids{$oidc};  $oid=$oidc if($pid); 
      #OR: if(not CUT_IDCG and $oid=~/_[CG]\d+$/){ (my $oidc=$oid) =~ s/_[CG]\d+$//; $pid= $pubids{$oidc}; }
      }
      
    # unless($pid) { 
    #   my($aoids)= $mattr =~ m/oid=([^;\s]+)/?$1:"";
    #   for my $d (split ",", "$oid,$aoids") {
    #     my($p,$od)= split",", ($puboidc{$d}||""); if($p){ $pid=$p; $oid=$od; last; } 
    #   }
    # }  
    
    if($pid or not $skipdropseqs) {
      # NOT IF _C split;  return 0 if($didseq{$pid});  #? yes? or have _C split parts same pid?
      # my $annorec= annotab2tblinfo( $oid, $annoth->{$pid}, $mattr, $pid); #??
      # my($pubid,$otr)= putPubGff($outh,$intypeo,$oid,$hd,$fa,$annorec); #??
      
      my($poid,$gid,$alti,$class,$aqual)= split"\t", $pubidinfo{$pid};
      if($kdxact =~ /cull|drop/) { $class=$kdxact.$class unless($class =~ /^(cull|drop)/); }
      
      if($class =~ /^alt/) { $source .= "alt"; }
      elsif($class =~ /^(main|noclass)/) {}
      else { my $st= $class; $st=substr($st,0,4) if(length($st)>4); $source .= $st; }
      
      ## FIXME: ncRNA class=noclass ? maybe add  $source .= "nc" as per alts
      #x else { my $st= $class; $st=substr($st,0,5) if(length($st)>5); $source .= $st; }
          #^uck: zf17evgm6ncull{alt}, zf17evgm6ncullm{main}, zf17evgm6nculln{noclass} ; keep this or no?

      ## UPD for gene names
      # (my $seqoid=$oid) =~ s/_[CG]\d+$//; # seq oid for names
      my $tblinfo= {};
      if( scalar(%genenames) ) { $tblinfo= parse_evgheader($oid,$mattr,0,$poid); } # uses %genenames ... FIXME alloids check of genenames

      #FIXME.culls:
      my $out1h= ($class =~ /^(cull|drop)/ and ref($cullh)) ? $cullh : $outh;

      my ($pidc,$isplit,$oidispid)=(0) x 9;
      for my $ex (@gff) { 
        #x $ex =~ s/(ID|Parent)=([^;\s]+)/$1=$pid/; 
        my($lid)= $ex =~ m/\b(?:ID|Parent)=([^;\s]+)/;
        unless($pidc) {
          $pidc= $pid || $lid;
          if(splitKEEP_C) {
            if($lid=~m/_C(\d+)/) { $isplit=$1; $pidc .= "_C$isplit" unless($pidc=~m/_C/); }
          } elsif($lid=~m/_C(\d+)$/){ $isplit=$1; }
          if(nogffDUPIDS) { return 0 if($didseq{$pidc}); }
          }
        $ex =~ s/(ID|Parent)=$lid/$1=$pidc/; 
        
        # strip excess annots, add some from annorec/tblinfo ?
        if($ex=~/\t($RNATYPES)/) { 
          chomp($ex); # NO LF
          $ex =~ s/;($ANNCUTRNA)=[^;\n]+//g; 
          # $ex =~ s/$/;$anadd/ if($anadd);
          ## FIXME: insert attr update from annorec= annotab2tblinfo(); which? options?
          
          ## FIXME2: updated genenanes now in tblinfo
          if(my $tdbx= $tblinfo->{'nameref'}) {
            my $tnap=$tblinfo->{'namepct'};
            my $tnam=$tblinfo->{'name'};
            my ($gdbx)= ($ex =~ m/\b[Dd]bxref=([^\s;]+)/)?$1:0; # db_xref ?
            my ($gnap)= ($ex =~ m/\b(?:namealn|namepct)=([^\s;]+)/)?$1:0;
            my ($gnam)= ($ex =~ m/\b[Nn]ame=([^=\n;]+)/)?$1:0;
            if($tdbx ne $gdbx or $tnap ne $gnap) {
              $ex.=";Dbxref=na" unless($gdbx); $ex =~ s/\b([Dd]bxref)=([^\s;]+)/$1=$tdbx/;
              $ex.=";namealn=na" unless($gnap); $ex =~ s/\b(namealn|namepct)=([^\s;]+)/$1=$tnap/;
              $ex.=";Name=na" unless($gnam); $ex =~ s/\b([Nn]ame)=([^\s;]+)/$1=$tnam/;
            }
          }
          
          my($xod)= ($ex =~ m/;oid=([^;\s]+)/)?$1:""; unless($xod){ $ex =~ s/$/;oid=/; }
          $oidispid= ($poid eq $pidc)?1:0;
          unless($oidispid or $xod=~m/\b$poid\b/) { my $nd="$poid,$xod"; $ex=~s/oid=$xod/oid=$nd/; $xod=$nd; } 
          $oidispid= ($lid eq $pidc)?1:0;
          unless($oidispid or $xod=~m/\b$lid\b/) { my $nd="$lid,$xod"; $ex=~s/oid=$xod/oid=$nd/; $xod=$nd; } 
          if($isplit) { # splitKEEP_C # formats: Split=2/3; | Split=2; | Split=C2; ..
            $ex =~ s/;Split=[^;\s]+//; $ex =~ s/;/;Split=$isplit;/; 
          }
          $ex.="\n"; # ADD BACK LF
        } else {
          $ex =~ s/;($ANNCUTEXON)=[^;\n]+//g; 
        }
        $ex =~ s/\t\S+\t/\t$source\t/;
        print $out1h $ex;
      }
      $didseq{$pidc}++; # if pidc ne pid ?  record splits in didseq? want in checktab?
      #? $didseq{$pid}++; # which? both?
      return 1; # return 1 + $isplit ?
    } else {
      return 0;
    }
  }
  
  my ($inh,$tblh,$inid,$gann,$gsrc);
  my @gene;
  open($outh,'>',$outgff) or return("",0,"error.$outgff");
  $cullh= undef; if($cullgff){ my $ok= open($cullh,'>',$cullgff); }

  print $outh "##gff-version 3\n"; # add hdr info
  print $cullh "##gff-version 3\n" if($cullh); # add hdr info
  
  for my $ingff (@ingff) {
    ($ok,$inh)= openRead($ingff);
    unless($ok) { loggit(1,"err reading $ingff"); next; }
    while (<$inh>) {
      if(/^\W/){ next; }  
      if(/\t($RNATYPES)/) { 
        $ngenes += putpubgff($inid,$gann,$gsrc,@gene) if($inid and $gann);
        @gene=(); $inid=$gann=0;
        my @v=split"\t";
        $gann= $v[8];
        $gsrc= $outsource || $v[1]; #? use input gff[1] col unless defined? not used here now?
        ## get from pubidinfo
        # if($isalt) { $gsrc.="alt"; } elsif($xxx) { $gsrc.="xxx"; }
        
        my($id)= m/\bID=([^;\s]+)/;  $inid= $id;
        # my($gpre,$gnum,$gd,$ti,$isplit)= evigene_idparts($id);
        push @gene, $_;
        
      } else {
        my($id)= m/\bParent=([^;\s]+)/; if($id ne $inid) { } # what?
        push @gene, $_;
      }
    
    }
    $ngenes += putpubgff($inid,$gann,$gsrc,@gene) if($inid and $gann);
    close($inh);
  } 
  # print $outh "# stats..\n";
  close($outh);
  close($cullh) if($cullh);

  #upd: -sort option; 
  #** THIS bestgenes_update sort is memory pig (reads all gff..); clear tmp hashes?
  # or defer gffsort to end of main caller? altbest2pubset.pl makePublicset()
  if($SORT_PUBSET) { 
    my $sortgff= $outgff.".sorted";
    my ($ns)= 0; # gffsort($outgff,$sortgff);  #via bestgenes_update.pl
		$APPgenesupdate= findevigeneapp("bestgenes_update.pl"); # unless($APPgenesupdate); #  
    
    ## drops longish logfile output: 
    my $dbf=($DEBUG)?"-debug":"";
    my $gffsortcmd="$APPgenesupdate $dbf -act sortgenes -cadd 'pubopt=1,addgene=1,chrsort=1,nochangeid=1'"
      ." -in $outgff -out $sortgff >& $outgff.sort.log "; 
      ## ,chrsort2nd=GRC ...    -mrnatype 'mRNA|ncRNA'
    my $err= runcmd($gffsortcmd);
    $err=666 if($gffsortcmd=~/MISSING/); # silly findapp=echo MISSING_xxx
    if($err == 0 and -s $sortgff) { 
      rename($outgff,"$outgff.unsort");  
      rename($sortgff,$outgff);
      push @tmpfiles, "$outgff.unsort", "$outgff.sort.log";
    } else {
      loggit(1,"sortgenes failed for $outgff")
    }
  } 

  my($nok,$nmiss,$misslist)= checkputids(\%didseq,$outgff,splitKEEP_C); # nsplit?
  
  push @publicset, $outgff; # push @tmpfiles, $inseq; #??
  return($outgff,$ngenes,"ok:$nok,miss:$nmiss"); 
}  


=item make_pubseq

=cut

sub make_pubseq
{
  my($inseq,$intype,$annothash,$outname)=@_;  # call for each input mrna,aa,cds
  my($notr,$nocds,$pubid)=(0) x 9;
  #above# our $CULLSEP= 1; # option
	
	# as for pubgff, allow inseq to be ref(array)
  my @inseq=();
  if($inseq and ref($inseq) =~ /ARRAY/){ @inseq= @$inseq; } 
  elsif($inseq =~ /\w/) { @inseq= ($inseq); }
  return("",0,"missing$intype") unless( @inseq > 0 and -f $inseq[0]); # or dryrun ..
	
	my $pubdir="publicset"; my $pubd=undef; my $ft;
#   if(not -f $annotab and $pubdir and -d $pubdir) {
#   	($pubd,$ft)= getFileset($pubdir,'.ann.txt',$pubd);  $annotab=$ft if($ft);
#   }

  # fixme: look in publicset/ for annot, submitset/ for outfa,tblout
  # FIXME/opt: separate output of _pub = ok seqs and _cull = cull seqs
  # FIXME: ncRNA now also, intype = mRNA only? separate seq files or not? pubids should have that in class info or aaqual
  
	my $outfa= $outname || $inseq[0]; $outfa =~ s/.gz$//; 
	my $cullfa= ($CULLSEP) ? $outfa."_cull.fa" : "";
	$outfa.="_pub.fa"; 
	my($pubsuf)= $outfa=~m/(\.\w+.fa)$/;
  unless(-s $outfa or not -d $pubdir) {	($pubd,$ft)= getFileset($pubdir,$pubsuf,$pubd);  $outfa=$ft if($ft); }
  return($outfa,$notr,"exists") if( -s $outfa); # or dryrun ..

	our($outh,$cullh,%didseq,$annoth,$intypeo);  %didseq=(); # no dup oids
  my ($ok,$inh,$oid,$fa,$hd,$itr,$otr,$ocds,$oerr)= (0) x 10;
  # our links to sub putpubseq
  $intypeo= $intype;
  $annoth = $annothash;
 
	# ($ok,$inh)= openRead($inseq); # loop over @inseq here
  $ok= open($outh,'>',$outfa); #  if($ok);
  ## unless($ok) { loggit(1,"ERR: make_pubseq $inseq TO $outfa"); return; }
  return("",0,"missing$intype") unless($ok); # or dryrun ..
  $cullh= undef; if($cullfa){ $ok= open($cullh,'>',$cullfa); }

#   sub putpubseq_OLD {
#     my($oid,$hd,$fa)= @_;
# 	  our(%didseq,$annoth,$intypeo,$outh); # no dup oids
# 	  return 0 unless($oid and $fa);
#     my $pid= $pubids{$oid}; # problem here is pubids.oid may have _G/C, not seq oid
#     
#     unless($pid) { 
#       my($aoids)= $hd =~ m/oid=([^;\s]+)/?$1:"";
#       my @aoids= split ",", $aoids; push @aoids, $oid;
#       
#       # FIXME: @p = sort keys %{$puballoids{$d}};  need to print dup seq for each @p
#       for my $d (grep { $puboidc{$_} } @aoids) {
#         my($p,$od)= split",", $puboidc{$d}; if($p){ $pid=$p; $oid=$od; last; }         
#       }
#       ## add $puballoids{$oid}{$pubid} check?  got some pid, wrong one?
#       if($pid and $didseq{$pid}) { $pid=0; } # fix dup map miss?
#       unless($pid) {
#       for my $d (grep { $puballoids{$_} } @aoids) {
#         my @p = sort keys %{$puballoids{$d}}; 
#         ## multimap seqs end up skipped here, 2+ @pid for same seq oid
#         @p= grep{ not $didseq{$_} } @p; #?? better, still missing some dup map _G seqs (get 1 copy)
#         if(@p){ $pid=$p[0]; $oid=$d; last; } 
#       }
#       }
#     }
#     
#     if($pid or not $skipdropseqs) {
#       return 0 if($didseq{$pid});  # this is skipping _G2 multimap seqs; want that? or dup seq output
#       my $annorec= annotab2tblinfo( $oid, $annoth->{$pid}, $hd, $pid); # have missing annorec
#       my($pubid,$otr)= putPubSeq($outh,$intypeo,$oid,$hd,$fa,$annorec,$pid); 
#       $didseq{$pubid}++; # any pubid != pid ??
#       return $otr;
#     } else {
#       return 0;
#     }
#   }
    
  sub putpubseq {
    my($oid,$hd,$fa)= @_;
	  our(%didseq,$annoth,$intypeo,$outh,$cullh); # no dup oids
	  return 0 unless($oid and $fa);
    my $pid= $pubids{$oid}; # problem here is pubids.oid may have _G/C, not seq oid
      # FIXME: @p = sort keys %{$puballoids{$d}};  need to print dup seq for each @p
    my @pid=();
    my($aoids)= $hd =~ m/oid=([^;\s]+)/?$1:"";
    
    my($kdxact, $kdxid)=(0,0);
    if($KEEPDROPX) {
      ($kdxact, $kdxid)= check_keepdrop('seqid',"$oid,$aoids");
      if(nogffDUPIDS and $kdxact =~ /drop/) { return 0 ; } #??
    }
    
    if($pid) { @pid=($pid); }    
    else { 
      my @aoids= split ",", $aoids; push @aoids, $oid;
        ## multimap seqs end up skipped here, 2+ @pid for same seq oid
      for my $d (grep { $puballoids{$_} } @aoids) {
        my @p = sort keys %{$puballoids{$d}}; 
        @p= grep{ not $didseq{$_} } @p; 
        if(@p){ @pid=@p; $oid=$d; last; } 
      }
      
      # if(0) { unless(@pid) {  # oid_Gn problem, use puboidc
      #   # losing many id_G2,3,4 here:
      #   #    if(not CUT_IDCG and $oid=~/_[CG]\d+$/){ (my $oidc=$oid) =~ s/_[CG]\d+$//; $puboidc{$oidc}= "$pubid,$oid"; }
      #   for my $d (grep { $puboidc{$_} } @aoids) {
      #     my($p,$od)= split",", $puboidc{$d}; if($p){ @pid=($p); $oid=$od; last; }         
      #   }
      # } }
    }
    
    my $nout=0;
    for my $pid (@pid) {
      next if($didseq{$pid});  # this is skipping _G2 multimap seqs; want that? or dup seq output
      my $annorec= annotab2tblinfo( $oid, $annoth->{$pid}, $hd, $pid); # have missing annorec
      #FIXME.culls: 
      my $evgclass= ""; # $annorec->{'evgclass'}||""; # FIXME this is missing for some, use pubidinfo{$pid}
      unless($evgclass) {
	     if(my $pin= $pubidinfo{$pid}){ my @pi=split"\t",$pin; $annorec->{'evgclass'}=$evgclass= $pi[3]; }
       else { $evgclass= $annorec->{'evgclass'}||""; }# 2nd choice now; mistakes from oldpub.aa ..
     }

      if($kdxact =~ /cull|drop/) { $evgclass=$kdxact.$evgclass unless($evgclass =~ /^(cull|drop)/); }
      my $out1h= ($evgclass =~ /^(cull|drop)/ and ref($cullh)) ? $cullh : $outh;
      my($pubid,$otr)= putPubSeq($out1h,$intypeo,$oid,$hd,$fa,$annorec,$pid,$nout); 
      $nout += $otr; 
      $didseq{$pubid}= $nout; 
    } 
    return $nout;
  }
  
  for my $inf (@inseq) {
  	($ok,$inh)= openRead($inf); # loop over @inseq here
    unless($ok) { loggit(1,"ERR reading $inf"); next; }
    while(<$inh>) { 
      if(/^>(\S+)/) { my $d=$1; # is this oid or pubid ???
        #x if($fa and $oid and ($pubids{$oid} or not $skipdropseqs)) #  and not $didseq{$oid}??
        if($fa and $oid) { $notr+= putpubseq($oid,$hd,$fa); }
        $oid=$d; $hd=$_; chomp($hd); 
        $fa=""; $itr++;
        }
      elsif(/^\w/) { chomp; $fa.=$_; }
    } 
    
    #x if($fa and $oid and ($pubids{$oid} or not $skipdropseqs)) #  and not $didseq{$oid} ??
    if($fa and $oid) { $notr+= putpubseq($oid,$hd,$fa); }
    close($inh); 
  } # @ inseq
  close($outh);
  close($cullh) if($cullh);
  
  #upd: -sort option
  if($SORT_PUBSET) { 
    my $sortfa= $outfa.".sorted";
    my($ns)= fasort($outfa,$sortfa);  # not for cull.fa
    if($ns>0 and -s $sortfa){ 
      rename($outfa,"$outfa.unsort");
      rename($sortfa,$outfa);
      push @tmpfiles, "$outfa.unsort"; 
    } else { # warn fail
      loggit(1,"sortfasta failed for $outfa")
    }
  } 
  
  ## count pubidinfo vs didseq
  my($nok,$nmiss,$misslist)= checkputids(\%didseq, $outfa);
  # write misslist to $outfa.misslist, or both ok/miss table
  
  push @publicset, $outfa; 
  push @publicset, $cullfa if($cullfa);  #? same public folder
  push @tmpfiles, $inseq; #??
  return($outfa,$notr,"ok:$nok,miss:$nmiss"); 
}

# 			if($fa and $oid) 
#       {
#         my $pid= $pubids{$oid}; # problem here is pubids.oid may have _G/C, not seq oid
#         unless($pid){ 
#           # my($p,$od)= split",", ($puboidc{$oid}||""); if($p){ $pid=$p; $oid=$od; } 
#           my($aoids)= $hd =~ m/oid=([^;\s]+)/?$1:"";
#           for my $d (split ",", "$oid,$aoids") {
#             my($p,$od)= split",", ($puboidc{$d}||""); if($p){ $pid=$p; $oid=$od; last; } 
#           }
#         } # also need oid= oidofpubid{pid} ?
#         # if(not $pid and $oid=~m/_[CG]/) { 
#         #  (my $soid=$oid) =~ s/_[CG]\d+//; $pid= $pubids{$soid}; $oid=$soid if($pid); }
#         if($pid or not $skipdropseqs) {
#      		my $annorec= annotab2tblinfo( $oid, $annothash->{$pid}, $hd, $pid); 
#      	  ($pubid,$otr)= putPubSeq($outh,$intype,$oid,$hd,$fa,$annorec); 
# 	      $notr+= $otr; $didseq{$pubid}++;
# 	      }
#      	}


sub checkputids {
  my($dididh,$listout,$checksplits)= @_;
  my @pubids= (sort keys %pubidinfo); # oids?
  # my @pok = grep{ $$dididh{$_} } @pubids;
  # my @pmiss= grep{ not $$dididh{$_} } @pubids;
  
  # FIXME need to correct for splitKEEP_C gff ids, not in pubidinfo, but in dididh
  #    if(splitKEEP_C and $lid=~m/(_C\d+)/) { my $sn=$1; $pidc .= $sn unless($pidc=~m/_C/); }
  my(@pok,@pmiss,%issplit);
  if($checksplits) {
    my @splid= map{ s/_C\d+//; $_; } (grep /_C/, sort keys %{$dididh});
    for (grep{ $pubidinfo{$_} and not $$dididh{$_} } @splid){ $issplit{$_}++; }  
    for (@pubids) { if($$dididh{$_} or $issplit{$_}) { push @pok, $_; } else { push @pmiss, $_; } }
  } else {
    for (@pubids) { if($$dididh{$_}) { push @pok, $_; } else { push @pmiss, $_; } }
  }
  
  my $nok= @pok; my $nmiss= @pmiss;
  if($listout) {
    loggit(LOG_DEBUG, "checkputids $listout.checktab ok=$nok,miss=$nmiss\n"); 
    open(LO,'>',"$listout.checktab"); # listout.check ?
    my $checkerr= $nmiss;
    for my $pd (@pubids) {
      #? use $$dididh{$pd}>1 as dup count? "okdup"
      #x my $nd= $$dididh{$pd} || 0;
      my $isplit= $issplit{$pd}||0;
      my $nd= $$dididh{$pd} || $isplit; # call split okdup?
      my $ok= ($nd > 1) ? "okdup$nd" : ($nd) ? "ok" : "miss";
      $checkerr++ unless($ok eq "ok");
      my $info= $pubidinfo{$pd} || "na";
      print LO join("\t",$pd,$ok,$info)."\n";
    } close(LO);
      #? same public folder ; if no err, put in @tmpfiles
    if($checkerr == 0) { push @tmpfiles,"$listout.checktab" ;}
    else { push @publicset,"$listout.checktab" ; }  #? same public folder ; if no err, put in @tmpfiles
  }
  return($nok,$nmiss,\@pmiss);
}

=item readAlignTab

  ($nmap,\%maps,\%alntab)= readAlignTab($alnfile);
  
  update to read gmap.attr (brief) and/or gmap.align (wide) tables
	evg align.tab
    GenomeID        gespan  geor    AQueryID        quspan  match   qlen    cov     pid     path    indels  nexon   splice  
    aalen   offs    aamap   sense   oid     tag
  brief map.attr
    AQueryID	cov	pid	splice/nexon	GenomeID:gespan:geor	path	oid

  ------> evg4corn_tgok2pub.align.tab.gz <------
  GenomeID	gespan	geor	AQueryID	quspan	match	qlen	cov	pid	path	indels	nexon	splice	aalenoffs	aamap	sense	oid	tag
  NC_001666.2	18-1190	.	Zeamay4gEVm016889t1	1-1173	1169	1173	100.0	99.7	0	0	1.1	0	353,90%,complete	41-1102	353	0	Zeamay4EVm014354t1,cornlo2m9slvelvk125Loc2249t1	gmap
  NC_001666.2	316-2184	.	Zeamay4gEVm016889t3	1-2007	1731	3298	60.9	91.2	0	144/6	1.1	230,21%,complete-utrbad	1587-2279	278	0	Zeamay4EVm014354t5,cornhi8ms9sgvelvk37Loc3472t1	gmap

=item mapqual_brief

  2018.05 update, new short mapqual statement as per
    evigene_pubsets.pm: sub gene_annot_brief()
    
  annotbrief is brief string of gene attributes,
  used in various evigene tables (eg  trclass2mainalt.pl):
    chrmap:alignpct[a],identpct[i],mrnalength[l],nexons[x],chrlocation
  Reproduces this annotation of trclass2mainalt.pl
  aaref:5767,dapsim:Dapsim1EVm000004t1,chrmap:100a,98i,25555l,33x,scaffold29:324762-356329:+ 
  aaref:496,zfish16nc:NP_938183.2,chrmap:99a,99i,1758l,9x,chr23:42516712-42528502:+
    
=cut


sub readAlignTab {
	my($maptab)= @_;
	
	my($ok,$hin)= openRead($maptab);  
	unless($ok) { loggit(1,"ERR readAlignTab( $maptab)"); return; }

	my $nmap=0;
	my (@hd,%hset,%alntab,%maps); 
	
	# my @alset=qw(QueryID qlen cov pid path nexon splice aalen sense); # align.tab
	# my @maset=qw(QueryID cov pid splice path); # map.attr SKIP FOR NOW ?
	my @geloc=qw( GenomeID gespan geor);
  my @alset=qw( QueryID qlen cov pid path nexon aalen sense); # align.tab
	my @hset= (@geloc,@alset);
	
	while(<$hin>) {
		unless(/^\w/) {
		  if(not @hd and /^#\w*Query\t/) { s/^#//; } else { next; }
		}
		chomp; my @v=split"\t";		
		my($td,$ql,$cov,$npath,$nx,$nspl,$aw,$sens,$oid,
		  $chr,$cbe,$cor,$pid,$alnrow) = (0) x 19;
		
		if(/^[A-Z]\w+ID\t/) { # ^GenomeID|AQueryID header
      @hd=@v; 
      map{ s/AQuery/Query/; s,splice/nexon,nexon,;  s/GenomeID:gespan:geor/GenomeID/; } @hd; #??
      %hset= map{ $_=>1 } @hset; 
      for my $h (@hd) { $hset{$h}=2 if($hset{$h}); } # check fields
      ##my @miss= grep{ $hset{$_} == 1 } @hset;
      ##if(@miss) { die "ERR: missing Mapinfo table fields: @miss\n"; }
      next;
		} elsif(@hd) {
			my %v=(); my(@kv);
			for my $i (0..$#v) { my $h=$hd[$i]; $v{$h}=$v[$i]; push @kv,"$h=".$v[$i]; }
		  $alnrow= join"\t",@kv;
		  ($chr,$cbe,$cor)= map{ $_||0 } @v{ @geloc }; 
		  if($chr=~/:/ and not $cbe) { ($chr,$cbe,$cor)= split":",$chr; }
		  ($td,$ql,$cov,$pid,$npath,$nx,$aw,$sens)= map{ $_||0 } @v{ @alset };  
			$oid=$v{oid}||$td;
		} else {
		  # guess or error?
		  next;
		}
		
		## split flag:   add NOPATH flag? or use cov == 0
		my $msplit= ($npath=~m/(C[123]):/)?$1:0;
		my $mpath = ($msplit or $npath eq "0/0")?0:($npath=~m,^(\d+/\d+),)?$1:($cov==0)?-1:0; # NOPATH=nocov=-1 ?
    # my $nsx= ($nx<2)?0:$nspl/$nx; 
    # my $as1=($sens<0 and $nx>3 and $nsx> 1.5 and $cov>90)?1:0; # certain?
    # my $as2=($sens<0 and not $as1 and (($nsx >= 1.5 and $cov>95) or ($nsx >= 1.8 and $nx > 5 and $cov>85)))?1:0; # likely
    # my $asense=(($as1 or $as2) and not $msplit)?"antisense":0;
		
		my $quals="";
    # if(MAPQ1805)  chrmap:99a,99i,1234l,9x,..
    $quals= sprintf "%da,%di,%dl,%dx,%s:%s:%s", ($cov,$pid,$ql,$nx,$chr,$cbe,$cor);
		#$quals.=",sense=$asense" if($asense);
		#$quals.=",split=$msplit" if($msplit); # or Spl:$npath if($msplit);
		#$quals.=",paths=$mpath" if($mpath);
		
		$maps{$td}= $quals; $alntab{$td}=  $alnrow; $nmap++;
		if($oid and $td ne $oid){ $maps{$oid}= $quals; $alntab{$oid}= $alnrow; }
  } close($hin);
		
  loggit(0,"readAlignTab($maptab)= $nmap");  
	return ($nmap>0)? ($nmap,\%maps,\%alntab) : (0,{},{}); 
}


=item gene_annot_brief

  ($annotbrief,$tblinfo)= gene_annot_brief($id,@gff_mrna_row);
  
  @gff_mrna_row is evigene-annotated mRNA row from  genes.gff
  with attributes (cov,pid,clen,nexon,cdsoff,aaqual,Dbxref,namealn,Name,.. )
  
  tblinfo is from parse_evgheader()
  
  annotbrief is brief string of gene attributes,
  used in various evigene tables (eg  trclass2mainalt.pl):
    aaref:alignlength,Dbxref,
    chrmap:alignpct[a],identpct[i],mrnalength[l],nexons[x],chrlocation
  with possibly more annots appended, eg. moltype when not mRNA
    mol:ncRNA when input type is ncRNA  

  Reproduces this annotation of trclass2mainalt.pl
  aaref:5767,dapsim:Dapsim1EVm000004t1,chrmap:100a,98i,25555l,33x,scaffold29:324762-356329:+ 
  aaref:496,zfish16nc:NP_938183.2,chrmap:99a,99i,1758l,9x,chr23:42516712-42528502:+

=cut

sub gene_annot_brief {
  my($id, @gff)=@_;
  my($stype,$at,$locus)=  ("","","");
  chomp(@gff);
  if(@gff>7) { ($stype,$at)=@gff[2,8]; $locus= sprintf "%s:%d-%d:%s",@gff[0,3,4,6]; } 
  elsif(@gff>1) { ($stype,$at)= @gff; ($locus)= $at=~m/\b(?:location|locus|maploc|loc)=([^;\s]+)/; }
  # my $stype= $gff[2]; # stype == mRNA|ncRNA.. add to tblinfo? to an
  
  # tblinfo: $pubid,$oid,$trlen,$cdsoff,$aaqual,$annogaps,$dbxref,$namepct,$gname,$cddname
  my $tblinfo= parse_evgheader($id,$at,0); # this works right; parses evg.gff or mrna attr
  
  my @mapq= map{ my($v)= $at=~m/\b$_=([^;\s]+)/?$1:0; int($v); } qw(cov pid clen nexon);
  my $mok= scalar( grep{ $_ ne "0" } @mapq );
  
  my($ndx,$nal)=(0,0);
  if($tblinfo) {
    $ndx= $$tblinfo{'dbxref'}||0; $ndx=~s/,$//;
    $nal= $$tblinfo{'namealn'}||0;
  } else {
    ($ndx)= $at =~ m/Dbxref=([^;\s]+)/?$1:0; $ndx=~s/,$//;
    if( my($nap)= $at =~ m/namealn=([^;\s]+)/) {  # namealn=58p,17578/30278,17622;
      ($nal)= $nap=~m/^\d+..(\d+)/; unless($nal) { ($nal)= $nap=~m/(\d+)/; }
    }
  }
  
  my $an= ($nal>0)?"aaref:$nal,$ndx,":"0,0,";
  my $cm= "";
  #x.if(@mapq and @gff>7){ $cm= sprintf "chrmap:%da,%di,%dl,%dx,%s:%d-%d:%s", @mapq, @gff[0,3,4,6]; }
  if($mok>1){ $cm= sprintf "chrmap:%da,%di,%dl,%dx", @mapq;  $cm.=",$locus" if($locus); }
  elsif($at=~m/(chrmap\W[^;\s]+)/){ $cm=$1; }
  #  $cm.=",$locus" if($locus); ??
  
  unless($stype=~/mRNA/){ 
    $cm .= ",mol:$stype"; #? stype= seqt= rnat= rtype= mol= ?
    $tblinfo->{'seqtype'}= $stype;
    } 
  $an .= $cm;
  return ($an, $tblinfo);
}

sub mapquals{ my($mq)=@_; my @m= $mq=~m/(\d+)a,(\d+)i,(\d+)l,(\d+)x/; return (@m>3)?@m:(0,0,0,0); }
sub mapinfostr{ my($mq)=@_;
  my($mapqual,$maploc)= (0,0);
  if( $mq=~m/(\d+a,\d+i,\d+l,\d+x).([\w:\.\+\-]+)/) { ($mapqual,$maploc)=($1,$2); }
  elsif( $mq=~m/(\d+a,\d+i,\d+l,\d+x)/ ) { ($mapqual)=($1); }
  return ($mapqual,$maploc);
}

=item annotab2tblinfo 

  parse putAnnot table row; make same rec as parse_evgheader()
  FIXME: allow variations in cols, use header, ID=PublicID
	#print $annoth join("\t",qw(PublicID OrigID TrLen CDSoff AAqual TrGaps Dbxref Namealign Product_Name CDD_Name))."\n" if($itr==1);
	#print $annoth join("\t",$pubid,$oid,$trlen,$cdsoff,$aaqual,$annogaps,$dbxref,$namepct,$gname,$cddname)."\n";
 ?? ERR if no annorec? == utrorf in .aa,cds not in mrna
 
 upd1803: add EvgClass:10 after CDD:9; pubinfo[class] supercedes annotab[class]
 
=cut

sub annotab2tblinfo 
{
  my($oidin,$tabrow,$fahdr,$pubidin)= @_;
	my($pubid,$oid,$trlen,$cdsoff,$aaqual,$annogaps,$dbxref,$namepct,$gname,$cddname,$evgclass)
			= (0) x 19; 
  $fahdr||="";

  # my $lotagpre= $settings{'LOCUSTAG'} || $ENV{'LOCUSTAG'} || ""; # global settings .. change for pack
  my $lotagpre= $ENV{'LOCUSTAG'} || ""; # global settings .. change for pack
  
  my $noann=0;
  $oid= $oidin; 
  $pubid= $pubids{$oid} || $pubidin || $oid; # ERR if missing? FIXME, oidin may be pubid, global %pubids has both pub/oid
  $cddname= $aaqual= $dbxref= 'na'; 
  my $fatinfo= ($fahdr) ? parse_evgheader($oidin,$fahdr) : 0; 
  my $genenameinfo= 0; # (ref($fatinfo) and $fatinfo->{'nameref'}) ? $fatinfo : 0;
  
  ## but see upd parse_evgheader CHECK_NAMEOIDS
  if( $genenames{$oid} ) { # give these precedence ??
    $gname  =  $genenames{$oid};
    $namepct=  $genenamepct{$oid} || 0;
    $dbxref =  $genedbxref{$oid}||"na";
    $cddname=  $cddnames{$oid}||"na";
    my %tblinfo= (pubid => $pubid, oid => $oid, name => $gname, 
      namepct => $namepct, nameref => $dbxref, cdd =>  $cddname);
    $genenameinfo= \%tblinfo;
  } elsif(scalar(%genenames)) {
    if(ref($fatinfo) and $fatinfo->{'nameref'}) { $genenameinfo= $fatinfo; }
  }

  use constant USE_ANNOCOLS => 1;
  # @ANNO_COLS= qw(PublicID OrigID TrLen CDSoff AAqual TrGaps Dbxref Namealign Product_Name CDD_Name Class Maploc Mapqual);
  # @ANNO_TBLINFO = qw(pubid oid trlen cdsoff aaqual trgaps dbxref namepct name cdd evgclass maploc mapqual);
  # ^^ use this instead of fixed set of col vars
  # 	my %tblinfo= (pubid => $pubid, oid => $oid,  protid => 0, locustag => 0,
  #       aaqual => $aaqual, trlen => $trlen, trgaps => $annogaps, cdsoff => $cdsoff, cdsor => $cdsor, 
  #       name => $gname, namepct => $namepct, dbxref => $dbxref, cdd =>  $cddname,
  #       evgclass => $evgclass, );  # added trgaps =>  $annogaps, evgclass

  # note: oid.in often == pubid, want instead original id from pubids/ann.txt tables
  my $pubidt= $pubid; # UPD: is tabrow:pubid always same as pubids val? which is certain?
  my %tblinfo= ( pubid => $pubid, oid => $oid, protid => 0, locustag => 0);
  $tabrow ||=""; # undef for some?
  
if(USE_ANNOCOLS) {
  my @tabr= (ref($tabrow) =~ /ARRAY/) ? @$tabrow : ($tabrow =~ /\t/) ? split("\t",$tabrow) : ();
  my $ncol= @tabr;
  if($ncol == 0) { $noann=1; } 
  elsif($ncol == @ANNO_TBLINFO) { # ncol <= @tblinfo ? 
    for my $i (0..($ncol-1)) { my $k=$ANNO_TBLINFO[$i]; $tblinfo{ $k } = $tabr[$i]; } # unless($k=~m/^(pubid|oid)$/)
  } else {
    my $nc= @ANNO_TBLINFO; $nc=$ncol if($ncol<$nc); # _min()
    for my $i (0..($nc-1)) { my $k=$ANNO_TBLINFO[$i]; $tblinfo{ $k } = $tabr[$i]; } #?  unless($k=~m/^(pubid|oid)$/)
  }
  # keep old way alive
  unless($noann) {
  ($pubidt,$oid,$trlen,$cdsoff,$aaqual,$annogaps,$dbxref,$namepct,$gname,$cddname,$evgclass)= 
    map{ $tblinfo{$_} || 0 } qw(pubid oid trlen cdsoff aaqual trgaps dbxref namepct name cdd evgclass);
  }
    
} else {  
	if(ref($tabrow) =~ /ARRAY/) {
		($pubidt,$oid,$trlen,$cdsoff,$aaqual,$annogaps,$dbxref,$namepct,$gname,$cddname,$evgclass)
 			= @$tabrow;
	} elsif($tabrow =~ /\t/) {
		chomp($tabrow);
		($pubidt,$oid,$trlen,$cdsoff,$aaqual,$annogaps,$dbxref,$namepct,$gname,$cddname,$evgclass)
 			= split"\t",$tabrow;
	} else {
	  $noann=1; 
	}		
}

  if($pubidt ne $pubid) {  # what?
    if($pubidt and $pubid eq $oidin){ $pubid= $pubidt; }
    else { } # should check errs
  } 
  
  # UPD: maybe add pubidinfo{pubid} : alloids? aaref? evgclass?
  my($pubclass,$alloids,$pubnotes,$seqtype,$maploc,$mapqual)=(0) x 9;
  if(my $pinfo= $pubidinfo{$pubid}) {
    # pubidinfo: ($oid,$gid,$alti,$class,$aqual,$pialn,$notes,$alloids); extended for realt_pubids, add Oids col
    # notes: aaref:146,sacavefish16nc:XP_016297840.1,chrmap:...
    my @pi=split"\t",$pinfo; 
    ($pubclass,$pubnotes,$alloids)= map{ $_ || "" } @pi[3,6,7];
    if($pubclass) {
      unless($evgclass) { $evgclass= $pubclass; }
      elsif($evgclass ne $pubclass) { $evgclass= $pubclass; } # always replace? check what?
    }
    if($pubnotes =~ m/:(\d+a,\d+i,\d+l,\d+x),([\w\.]+:\d+-\d+:.)/) { $mapqual=$1; $maploc=$2; } 
    if($pubnotes =~ m/,mol:(\w+)/) { $seqtype=$1; } # ncRNA fix.. 
    #no: my ($aaref)= $pi[6] =~ m/aaref:(\S+)/; $aaref=~s/,chrmap.*//; ..
  }
  
  my $cdsor=1; if($cdsoff =~ s/:([+.-])$//) { $cdsor=$1; } # drop strand if there. set cdsor  ?
  $dbxref =~ s/,$//;

if(USE_ANNOCOLS) {
   # check  $gname genenameinfo below
   # FIXME: noann needs to update tblinfo fields from pubidinfo
   if($noann) {
    # @ANNO_TBLINFO = qw(pubid oid trlen cdsoff aaqual trgaps dbxref namepct name cdd evgclass maploc mapqual);
    my @av= ($pubidt,$oid,$trlen,$cdsoff,$aaqual,$annogaps,$dbxref,$namepct,$gname,$cddname,$evgclass, $maploc, $mapqual); # vals we have 
    for my $i (2..$#ANNO_TBLINFO){
      my $k=$ANNO_TBLINFO[$i];
      unless(exists $tblinfo{$k}) { $tblinfo{$k}= $av[$i]; }
      }
    # $tblinfo{'pubid'}= $pubidt; $tblinfo{'oid'}= $oid;   
    $tblinfo{'evgclass'}= $evgclass if($evgclass);   
    # $tblinfo{'aaqual'}= $aaqual if($aaqual);   
    # $tblinfo{'dbxref'}= $dbxref if($dbxref);   
    # $tblinfo{'cdd'}= $cddname if($cddname);   
   }
} else {  
  %tblinfo= (pubid => $pubid, oid => $oid,  protid => 0, locustag => 0,
      aaqual => $aaqual, trlen => $trlen, trgaps => $annogaps, cdsoff => $cdsoff, cdsor => $cdsor, 
      name => $gname, namepct => $namepct, dbxref => $dbxref, cdd =>  $cddname,
      evgclass => $evgclass, );  # added trgaps =>  $annogaps, evgclass
}

  $tblinfo{'alloids'}= $alloids if($alloids);
  my $anninfo= \%tblinfo;
  
    # FIXME2: ** check annotrec is uptodate, match oid at least w/ evgheader
    # .. add this check to annotab2tblinfo() .. option to return faheader version?
    ## FIXME: Selcstop=  only in fahdr now ?
    # >Funhe2Exx11m129535t7 aalen=321,92%,complete,selcstop; Selcstop=655,676,694,742,748,754,802,832,838,883,889,910,934,940,961,967,; clen=1039; strand=+; offs=16-981;  pubid=Funhe2EKm029571t1; Name=Selenoprotein P,  
  if($fahdr) {
    #above# my $fatinfo= parse_evgheader($oidin,$fahdr); 
    if($fatinfo->{Selcstop}) { 
      $anninfo->{Selcstop}=$fatinfo->{Selcstop}; 
      $anninfo->{aaqual}=  $fatinfo->{aaqual}; # qual,selcstop
      }
    if(my $ft=$fatinfo->{'seqtype'}) { $seqtype= $ft; } # type=mRNA; trust this?
    
    my $hdif="";
    my @KS=qw(oid pubid cdsoff aaqual trlen); # fckn.mess cdsoff/off has :+/- sometimes or not
    @KS= grep{ not m/cdsoff/ }@KS if($fahdr=~/strand=/); # if($fatinfo->{strand} eq "-");
    ## mrna: aalen=1884,83%,complete; clen=6736; strand=+; offs=99-5753; 
    ## ann.txt: 6736    99-5753:+       1884,83%,complete 
    ## ** cdsoff ignore :+ end or fix fahdr reader
      # which are required to match? cdsoff, aaqual?, trlen? 
      # * add aaqual, need proper complete/partial flag
    for my $ks (@KS) {
      my($avs,$fvs)= ($anninfo->{$ks},$fatinfo->{$ks});
      if($fvs and $fvs ne $avs) {
        my $dok=($ks =~ /^(oid|cdsoff)/ and $avs and ($fvs=~m/$avs/ or $avs=~m/$fvs/))?1:0;
        $hdif .="$ks=$fvs/$avs," unless($dok);
        }
    }

  ## annot mismatch problems: non-mrna strand- cdsoff diff mrna cdsoff
  ## utrorf trlen differs from input.. should be fixed somewhere, but not in seqfile.fa hdr
  ## ?? okayset/xxxx.utrorf.mrna
      
    if($hdif) {
      my $fatok=($FAHDRisValid or ($noann and $hdif=~/cdsoff=/))?1:0;
      if($fatok) { 
        my $annold=$anninfo; # or use %tblinfo
        ## $FAHDRisValid >1 or == 2 means merge only missing vals !!!
        
        $anninfo=$fatinfo; #?? merge best of both ?? where hdif is missing/0
        $anninfo->{pubid}= $pubid; #ensure same 
        ## fix2: copy bits from wrong anninfo: dbxref,  name???? missing genenameinfo ?
        if($FAHDRisValid > 1) {
          for my $ks (sort keys %$annold){ 
            my($anv,$fav)= ($annold->{$ks},$fatinfo->{$ks});
            if( $anv and (not $fav or length($anv)>length($fav)) ) { $anninfo->{$ks}=$anv; }
          }
        }
        if( length($annold->{dbxref}) > length($anninfo->{dbxref}) and $annold->{dbxref}=~m/:/) {
          $anninfo->{dbxref}= $annold->{dbxref}; # ok?
          }
      }
      ## both info are missing cdsoff for some .. recalc?
      #** cancel this err for mrna uvcut= annots
      #Off: my $nodif=($noann or $fahdr=~/uvcut=/)?1:0; #?? $fatinfo->{'uvcut'}
      #Off: loggit(1,"ERR: annot mismatch $pubid/$oidin:",$hdif) unless($nodif); 
        #**TOO MUCH annot mismatch n=686532 : due to tr<>mrna revcomp ? cdsoff changed 
    }
  }
  
  $anninfo->{'seqtype'}= $seqtype if($seqtype); # after both pubid + fahdr check

  ## fix3: names precedence; copy cdd, name? namepct? nameref=dbxref
  if(ref $genenameinfo) {
    for my $k (qw(name namepct nameref cdd)) {
      if($genenameinfo->{$k}=~/\w/) { $anninfo->{$k}= $genenameinfo->{$k}; }
    }
    
    ## this needs to be in parse_evgheader, also
    my $adx= $anninfo->{dbxref} || ""; my @adx=split",",$adx;
    my $ndx= $genenameinfo->{nameref} || ""; my @ndx=split",",$ndx;
    ## fixme.. drop adx when gnameref is similar:  gnameref=SwissProt:TDRD7_HUMAN; adxref=TrEMBL:TDRD3_HUMAN,TrEMBL:UniRef50_Q9H7E2
    my @addx=(); my @dropx=();
    for my $nd (reverse @ndx) { 
      unless($adx =~ /$nd/) { unshift(@addx,$nd); 
        if($nd=~/SwissProt|TrEMBL|UniProt/) { my($sp)= $nd=~m/(_\w+)$/; 
          if($sp and $adx =~ /$sp/){ push(@dropx, grep( /$sp/, @adx)); @adx= grep{ not m/$sp/ } @adx; }
        }
      } 
      # unless($adx =~ /$nd/) { unshift(@adx,$nd); } 
    }
    unshift(@adx, @addx);
    $anninfo->{dbxref}= join",",@adx;
  }

  ## add pubid to locustagid for ncbi tsa
  #m2t.LOCUSTAG =     Locus tag prefix:       RC70 (SAMN02116571)    <<dgg is this ok? conflict w/ other kfish EST project?
  unless($anninfo->{locustag}) {
    if($lotagpre) {  # my $lotagpre= $settings{'LOCUSTAG'}
      my($pubidnum)= $pubid =~ m/$IDPREFIX(\d+)/;
      $anninfo->{locustag}= $lotagpre.'_'.$pubidnum if($pubidnum);
    }
  }
  
	return $anninfo;
}



sub evgIdParts { # or evigene_idparts 
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
      # $gnum=0; # BUG fix??
    }
  }
  
  return ($gpre,$gnum,$ti,$isplit,$gdup); #?
  #OR:   if($gdup) { $gd.=$gdup; $gnum=0; } # gnum invalid?
  #OR:   return($gpre,$gnum,$gd,$ti,$isplit);
}

sub evgIdOfParts {  # was evigene_idOfparts(@parts)
  my($gpre,$gnum,$ti,$isplit,$gdup)= @_; # new evgIdParts
  #o my($gpre,$gnum,$gd,$ti,$isplit)= @_; # old evgIdParts
  $gpre ||= $IDPREFIX;
  $ti ||=1;
  $gnum= sprintf( '%06d', int($gnum)); # leading 000
  my $id= $gpre . $gnum . "t$ti"; 
  $id.="_C$isplit" if($isplit); 
  # $id.="_G$gdup" if($gdup);
  return($id);
}

sub make_pubid # add preseveOld %newids here? use altnum w/ it or preserve old altnum?
{
  my($oid, $pubidnum, $altnum, $preservedIds)= @_;
  $pubidnum_start= $pubidnum; #?
  my $alti=$altnum;

  my($pubid,$pubgene);
  if($preservedIds and ref($preservedIds)) {
    if(my $pid= $preservedIds->{$oid}) {
      ($pubid,$pubgene)=($pid,$pid); #? yes? no?
      if(my($pgene,$palt)= $pid=~m/(\w+)t(\d+)$/) {
        $pubgene= $pgene;
        $alti= $altnum; # $palt; #? preserve altnum? or caller's altnum 
        # if($palt ne $altnum) ..
        $pubid   = $pubgene . sprintf( $altid_format, $alti);  
        return($pubid,$pubgene,$altnum);
      }
    }
  }
  $pubgene = sprintf( $pubid_format, $pubidnum); 
  $pubid   = $pubgene . sprintf( $altid_format, $alti);
  return($pubid,$pubgene,$alti);
}

# sub init_PUBSET # make_IDPREFIX handles these?
# {
#   # global var defaults, if unset by caller
#   # our($pubid_format,$altid_format, $pubidnum_start, $GBPROID, );
#   $IDPREFIX ||= "EVGm";
#   $pubid_format = $IDPREFIX.'%06d' unless($pubid_format);  # digits == 6
#   $altid_format = 't%d' unless($altid_format);  
#   $pubidnum_start ||= 0;
#   unless($DATE) { $DATE=`date '+%Y%m%d'`; chomp($DATE); } # default is today; use perl func?
# }

# my($pubid_format,$altid_format,$GBPROID)= make_IDPREFIX($sampleid||undef);
sub make_IDPREFIX
{
  my($existingID)= @_;
  my $digits=6; # or 7?
  # my $nd= ( $IDPREFIX =~ s/(0\d+)$// ) ? length($1) : $digits;
  # my $nd= $digits; 
  
  my $ChangeDefaultIdPrefix= 1; #? ($IDPREFIX eq $DEFAULT_SETTINGS{'IDPREFIX'}) ? 1:0; # def: EVGm
  $IDPREFIX ||= "EVGm"; #upd? EVm?
  
  if($existingID and $existingID !~ m/^$IDPREFIX/) {
    ## ($gpre,$gnum,$ti,$isplit,$gdup)= evgIdParts($existingID);
    $existingID=~ s/_[CG]\d+$//; 
    $existingID=~ s/t\d+$//;
    my($prefix,$nums)= $existingID =~ m/^(\w+[a-zA-Z])(\d\d\d+)$/; ## prefix may have numbers, how many? 
    if($prefix) { $IDPREFIX= $prefix; my $nd= length($nums); $digits= $nd if($nd > $digits); }
    $ChangeDefaultIdPrefix=0;
  }
  
  #  get species/organism from that, make IDPREFIX from Spp.org. abbev. 
  if($ChangeDefaultIdPrefix and okORGANISM()) {
    my $prefix="";
    my($gen,$spp)=split /[\s\_]/, $ORGANISM, 2;  
    if($spp and length($gen)>1) { $prefix= ucfirst(substr($gen,0,3)) . lc(substr($spp,0,3)) . 'EVm'; }
    else { $prefix= ucfirst(substr($ORGANISM,0,6)) . 'EVm'; } # was "EGm"
    $IDPREFIX= $prefix;
  }
  
  my $pubid_format = $IDPREFIX.'%0'.$digits.'d'; 
  my $altid_format = 't%d';  
  #^ set globals of same name?
  
  unless($DATE) { $DATE=`date '+%Y%m%d'`; chomp($DATE); } # default is today; use perl func?
  my $GBPROID= $IDPREFIX."_".$DATE; # "cacao11evigene_20120827";
  
  return($pubid_format,$altid_format,$GBPROID);
}



=item CDD2name

  convert CD name to ncbi/uniprot acceptable protein name
  NCBI eukgenosub_annotation "<domain|repeat>-containing protein". e.g. "PAS domain-containing protein 5".

   PH_IRS, Insulin receptor substrate (IRS) pleckstrin (PH) domain protein, putative       lcl|Funhe2Exx11m084853t1:24-236 
   ^^ need new CDD to NCBI name routine: cut out leading symbol XXX_yyy, for note?
      .. remove fragment,partial,-like,putative and such terms
      .. add ' domain-containing protein' unless has domain
   FOG: bug, whatisit? CDD: FOG: Zn-finger   CDD:227381 
   bad: "ZP, Zona pellucida  domain-like fragment domain-containing protein"

=cut

sub CDD2name
{
  local $_= shift;  my $fornote= shift or 0; 
  my $csym="";
  s/CDD:\s*//; s/FOG:\s*//;
  if(s/^(\S+),\s*(\w)/$2/) { $csym=$1; }
  s/(fragment|putative|\-like)\s*/ /g; s/_/ /g;
  ## Changed 'Domain of unknown function DUF4218' to 'protein of unknown function DUF4218' ..
  ## better productname: Unknown function DUFnnn domain-containing protein
  if(/^Domain of unknown function/i and not $fornote) {
    s/Domain of unknown function (\w.*)/Unknown function $1 domain-containing protein/i ; # use key $UNK_CDD ??
    s/Domain of unknown function/Protein of unknown function/i ; # use key $UNK_CDD ??
  }
  
  if(m/\bdomain/i) { s/\bdomain protein// if($fornote); }
  else { s/\bprotein//; $_.=" domain-containing protein" unless($fornote); }
  return (wantarray) ? ($_,$csym) : $_;
}



=item  productname ncbi tbl2asn bugs

  * Uhggg, fixed problem in evigene/gbrnaseq parser of names,
    was getting product names from other than CDS subrecords, with
    other-than-ncbi-proper gene names. ie. all these 'fatal' things
    
  since almost all of names for pig evg set are NCBI RefSeq,
  the name warns/errors are case of NCBI disagreeing with NCBI
  e.g.
  FATAL: 1 features contains 'pseudogene'(FATAL)
  human18nc:NP_001342187.1   histone cluster 2, H3c pseudogene       gid:391769      t1      NM_001355258    sym:H3  ncbirna:NM_001355258
  human18nc:NP_001342212.1   ribosomal protein SA pseudogene 58      gid:388524      t1      NM_001355283    sym:RPSAP58     HGNC:HGNC:36809,ncbirna:NM_001355283
  human18nc:NP_001345618.1        tubulin beta 8 pseudogene 12    gid:260334      t1      NM_001358689    sym:TUBB8P12    HGNC:HGNC:24983,ncbirna:NM_001358689
    ^^ pseudogene is really in many protein names from NCBI RefSEq Human 2018 set
 
  FATAL: 1 features contains 'open reading frame'(FATAL)
  human18nc:NP_001340422.1        ASNSD1 upstream open reading frame protein      gid:110599588   t1      NM_001353493    sym:ASDURF      HGNC:HGNC:53619,ncbirna:NM_001353493
   
  FATAL: Remove organism from product name
            1 features contains 'staphylococcal'
  human18nc:NP_055205.2   staphylococcal nuclease domain-containing       gid:27044       t1      NM_014390       sym:SND1        CCDS:CCDS34747.1,HGNC:HGNC:30646,MIM:602181,ncbirna:NM_014390
        ^^ NCBI RefSeq human name
            
  FATAL: Possible parsing error or incorrect formatting; remove inappropriate symbols
            17 features Contains unbalanced brackets or parentheses
  CUTLEADCRAP bug here?  ^^
    DiscRep_SUB:SUSPECT_PRODUCT_NAMES::19 features Contains unbalanced brackets or parentheses
    mrna_pub.split9:CDS     nesprin-1 isoform 1 (nesprin-1 giant or Susscr4EVm000002t135:7-13170    
    mrna_pub.split9:CDS     F-actin]-monooxygenase MICAL2 isoform e Susscr4EVm000595t41:291->2984   
    mrna_pub.split9:CDS     F-actin]-monooxygenase MICAL1 isoform 3 Susscr4EVm001749t35:316-2871    
    mrna_pub.split9:CDS     Pyruvate dehydrogenase (acetyl-transferring)]   Susscr4EVm008796t4:296-1444     
 
  9 features contains '. ' << check these
    where is '{ECO:0000250}' from? NCBI Human refseq
    ^^ this is mis-parse by evg script from ncbi.rna.gbff, i.e. {ECO:..} is NOT CDS name, but from other misc feature product name
  human18nc:NP_001078851.1	Saposin D-like. {ECO:0000250}	gid:768239	t1	NM_001085382	sym:PSAPLCCDS:CCDS47009.1,HGNC:HGNC:33131,ncbirna:NM_001085382
  human18nc:NP_000154.1	Growth hormone-binding protein. {ECO:0000250}	gid:2690	t1	NM_000163	sym:GHR	CCDS:CCDS3940.1,HGNC:HGNC:4263,MIM:600946,ncbirna:NM_000163
  human/ncbigenes/human18ncbi_grch38_rna.gbff.gz 
  LOCUS       NM_000163               4563 bp    mRNA    linear   PRI 26-FEB-2018
  DEFINITION  Homo sapiens growth hormone receptor (GHR), transcript variant 1,
            mRNA.
     CDS             193..2109
                     /gene="GHR"
                     /gene_synonym="GHBP; GHIP"
                     /note="isoform 1 precursor is encoded by transcript
                     variant 1; growth hormone binding protein; somatotropin
                     receptor; serum binding protein; GH receptor"
                     /codon_start=1
                     /product="growth hormone receptor isoform 1 precursor" << NOT THIS NAME
                     /protein_id="NP_000154.1"
                     /db_xref="CCDS:CCDS3940.1"
                     /db_xref="GeneID:2690"
                     /db_xref="HGNC:HGNC:4263"
                     /db_xref="MIM:600946"
                     .....
     mat_peptide     247..960
                     /gene="GHR"
                     /gene_synonym="GHBP; GHIP"
                     /product="Growth hormone-binding protein. {ECO:0000250}" << THIS NAME collected
                     /experiment="experimental evidence, no additional details
                     recorded"
                     /note="propagated from UniProtKB/Swiss-Prot (P10912.1)"
    
    DiscRep_SUB:SUSPECT_PRODUCT_NAMES::18 features contains '. '
    mrna_pub.split9:CDS     Agrin N-terminal 110 kDa subunit. {ECO:0000250} Susscr4EVm000363t18:105-6248    
    mrna_pub.split9:CDS     Notch 4 intracellular domain. {ECO:0000250}     Susscr4EVm000354t56:132->5819   
    mrna_pub.split9:CDS     Presenilin-2 NTF subunit. {ECO:0000250} Susscr4EVm005619t11:548-1891    
    mrna_pub.split9:CDS     Transmembrane protein. {ECO:0000250}    Susscr4EVm006911t169:24-440     
  
=cut

sub productname 
{
  my($gname,$aaqual,$napct, $FOR_NCBI)=@_;  #? add FOR_NCBI flag?

  my $cdsym="";
  $napct ||= 0; $FOR_NCBI ||=0; 
  # FIXME: input napct = 83%,228/276,232 , not numeric; return orig but need num
  my($dnapct)= ($napct =~ m/^(\d+)/)?$1:0;
  local $_= $gname || ""; 
  
  ## this sub belongs with evigenes/prot/protein_names.pm
  ## see evigene/scripts/prot/protein_names.pm : nameclean() .. BUT these are supposed to have been run thru that.
  ## add option to call nameclean() here .. again?
  
  ## FIXME1806: names w/ leading special chars, from NCBI Refseq so must be allowed..
  # (E3-independent) E2 ubiquitin-conjugating enzyme
  # >> thissub changed to 'E3-independent) E2 ubiquitin-conjugating enzyme'
  #...
  # [Pyruvate dehydrogenase  << problem unbalanced -[ no-]
  # [F-actin]-monooxygenase MICAL1 isoform 1

  use constant NAMEFIX1_FRAG => 1; # FIXME: NCBI hates fragment now, wants ', partial' ??
  ## fix3? -like should never end name, or follow punctuation .. always "blah blah-like protein" ??

  #noname? should this be done already?# 
  $_='' if($dnapct>0 and $dnapct < $MIN_IDLIKE);
  
  my $BUSCOID= 'EOG'; # hack fix input names, EOG090B07VD EOG090B08S3 .. some are 'Uncharacterized protein'
  # FIXME: need something better than 'hypothetical protein' for naming BUSCO conserved genes.
  # "vertebrate|CLADE conserved protein" ?
  
  if($_ eq 'na' or /^($NAME_NONE)$/i or /^($NAME_NONE) protein\W*$/i) { 
    $_='';
  } elsif(/^($NAME_NONE)\s*\w+/) { 
    if(/^$NAME_UNK\s*\w+/) { s/^$NAME_UNK\s*//; }
    else { s/^($NAME_NONE)\s*//; s/^protein\s*//i; }
  } 
  
  if(/\w\w/) {
    if(s/\s*\((\d+)\%\w*\)\W*$//) { my $p=$1; $napct=$p if($p>$dnapct); } # dang remove namepct: (99%P)
    s/LOW QUALITY PROTEIN:\s*//; # not part of true name
    if($FOR_NCBI and s/\s(LOC\d\d+)\b//) { my $cutid=$1; } ## uncharacterized protein LOC106505107  << should LOC ids be removed?
    
    if(/^($BUSCOID\d\d\w+)$/) { my $cutid=$1; $_= "vertebrate conserved protein"; } # drop id here, is dbxref
    
    ## TE: CDD: RP_RTVL_H_like, Retropepsin .. names bug, should leave off TE: 
    s/^TE:\s*//; # names bug, leave out TE: here
    my $fornote = ($FOR_NCBI & FOR_NOTE)?1:0;
    ($_,$cdsym)= CDD2name($_,$fornote) if(/^CDD:/);
    # if(/^CDD:/) { s/CDD:\s*//; s/\-like\s*/ /; unless(m/\bdomain/i) { s/\bprotein//; $_.=" domain-containing protein"; } }  
      
    # 45 WARNING: SEQ_FEAT.BadTrailingCharacter  Funhe2Exx11m001464t1
      # WARNING: valid [SEQ_FEAT.BadTrailingCharacter] Protein name ends with undesired character FEATURE:
      # "Roundabout, axon guidance receptor," n=32 are Roundabout
      #  "LATS, large tumor suppressor," << ','  n=13 rest are LATS,
      
    # 8 WARNING: SEQ_FEAT.ProteinNameEndsInBracket .. Funhe2Exx11m008798t7/Funhe2EKm026612t + other alts of same
      # Protein name ends with bracket and may contain organism name FEATURE: 
      # Dimethylaniline monooxygenase [N-oxide-forming]  << this appears to be Human HUGO approved name FMO3_HUMAN/FMO5_HUMAN
      # Funhe2EKm026612t6/Funhe2Exx11m008798t7 == TrEMBL:UniRef50_G7MFB6,TrEMBL:G7MFB6_MACMU,mayzebr:XP_004557286.1,Omcl:FISH11G_G24186

    if(NAMEFIX1_FRAG) {  
      s/[,]?\s*(partial|fragment)\b//i; # always drop these? add , partial below     
    } else {
      # 4 WARNING: SEQ_DESCR.InconsistentProteinTitle : "partial" name on complete aa
      if($aaqual =~ /complete/ and m/\b(partial|fragment)\b/i) { 
        s/[,]?\s*(partial|fragment)\b//i; 
      } elsif($aaqual =~ /partial/ and not m/(partial|fragment)/i) {
        $_ .= " fragment"; # FIXME: NCBI wants ', partial' ??  fragment is Uniprot's preferred syntax?
      }
    }    

    if($FOR_NCBI and /isoform/) { s/[ ,;\.]*\s+isoform \w+//; } # isoform cut: always, or only for NCBI or other usage?
    
    # see protein_names:nameclean()
    use constant CUTLEADCRAP => 0; # should be option?
    use constant CUTLEADCRAPb => 0; # should be option?
    ## this is still buggy, leave in all leading symbols for now
    if(CUTLEADCRAPb and m/^\W/){
      s/^[^A-Za-z0-9_\(\[]+//;
      
    } elsif(CUTLEADCRAP and s/^(\W+)//) {  # no leading crap 
      my $cutb=$1;
      # upd1806: allow some leading crap (ie refseq like '[bobo] was here'
      # cant cut leading '([' without cutting following closing '])'
      if($cutb eq '[') { s/\]//; } elsif($cutb eq '(') { s/\)//; } elsif($cutb eq '{') { s/\}//; }
      ##?? not working, why not?
      ## 'E3-independent) E2 ubiquitin-conjugating enzyme' Susscr4EVm116640t1
    }
    
    s/\s*[\/\.,;_-]+\s*$//; # trailing punc; BadTrailingCharacter
    if(NAMEFIX1_FRAG) { # add after trim; BUT NOT FOR_NCBI .. they add it ?? 
     unless($FOR_NCBI) { $_ .= ", partial" if(/\w\w/ and $aaqual =~ /partial/); } # FIXME: NCBI hates fragment now, wants ', partial' ??
    }    
  }
	## change b/n NAME_UNKUNIP and NAME_UNKNCBI depending on outputs?  prefer NAME_UNKUNIP but for TSA submit
	unless(m/\w\w/) { $_= ($FOR_NCBI)? $NAME_UNKNCBI: $NAME_UNK; }  # SEQ_FEAT.MissingCDSproduct; need some name
	return($_,$napct,$cdsym);
}


1;

__END__
