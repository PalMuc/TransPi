#!/usr/bin/env perl -w
## AUTHORS: Li Li, Feng Chen <fengchen@sas.upenn.edu>
## ORTHOMCL [2007-04-04] Version 1.4
## mod for evigene uses, d.g.gilbert, 2013 : evigene/scripts/omcl/orthomcl_evg.pl

## Copyright (C) 2004~2006 by University of Pennsylvania, Philadelphia, PA USA.
## All rights reserved.

## Before orthomcl.pl can be used, some variables (including directory variables
## or parameter variables) in orthomcl_modevg.pm need to be set, as described in
## README.

my $starttime = `date`;
use strict;
use Getopt::Long;
use File::Basename;
use orthomcl_modevg;  # dgg

my ($mode,$fa_files,$pv_cutoff,$pi_cutoff,$pmatch_cutoff,%blast_flag,$inflation,$maximum_weight);
my ($usr_blast_file,$usr_bpo_file,$usr_gg_file,$usr_taxa_file,$former_run_dir);         # For Mode 2, 3 or 4

my $DEBUGUNDEF=0; #dgg
my $DEBUGMAT=0;

my $command=basename($0)." ".join(' ',@ARGV)."\n";

&GetOptions(
			"mode=s"              => \$mode, # dgg: was =i, use named modes
			"fa_files=s"          => \$fa_files,
			"pv_cutoff=s"         => \$pv_cutoff,
			"pi_cutoff=f"         => \$pi_cutoff,
			"pmatch_cutoff=f"     => \$pmatch_cutoff,
			"inflation=f"         => \$inflation,
			"maximum_weight=i"    => \$maximum_weight,
			"blast_file=s"        => \$usr_blast_file,
			"bpo_file=s"          => \$usr_bpo_file,
			"gg_file=s"           => \$usr_gg_file,
			"taxa_file=s"         => \$usr_taxa_file,
			"former_run_dir=s"    => \$former_run_dir
);

if (!defined $mode) {printHelp();}

#set the default
$pv_cutoff      = $pv_cutoff      ? $pv_cutoff      : $BLAST_PVALUE_CUTOFF_DEFAULT;
$pi_cutoff      = $pi_cutoff      ? $pi_cutoff      : $PERCENT_IDENTITY_CUTOFF_DEFAULT;
$pmatch_cutoff  = $pmatch_cutoff  ? $pmatch_cutoff  : $PERCENT_MATCH_CUTOFF_DEFAULT;
$inflation      = $inflation      ? $inflation      : $MCL_INFLATION_DEFAULT;
$maximum_weight = $maximum_weight ? $maximum_weight : $MAX_WEIGHT_DEFAULT;

if ($BLAST_FORMAT eq 'full') {
	%blast_flag=( 'm'      =>0,
				  'format' =>'blast',
				  'hsp'=>1
		);
} elsif ($BLAST_FORMAT eq 'compact') {
	%blast_flag=( 'm'      =>8,
				  'format' =>'blasttable',
				  'hsp'=>1
		);
} else {
	die "\$BLAST_FORMAT can only be 'full' or 'compact'!\n";
}

my (%connect, %ortho); # dgg: these two globals must be filled for COMPUTE_MCL_MATRIX;  %weight also **

## dgg : skip for now, require .bpo,.gg inputs
# if ($mode == 1) { ## or mode =~ /blastall/
# 	if (defined $fa_files) {
# 		&constructDirectory($starttime);
# 		my %seq_len=%{constructAllFasta($fa_files,$all_fa_file)};                  #construct all.fa file
# 		&executeFORMATDB($all_fa_file);                                            # and run blast
# 		&executeBLASTALL($all_fa_file,$blast_file,$all_fa_file,$pv_cutoff,\%blast_flag);
# 		&blast_parse($blast_file,$bpo_file,$pv_cutoff,\%seq_len,\%blast_flag) unless (-e $bpo_file);
# 	} else {dieWithUnexpectedError("In Mode 1, NAMES OF FASTA FILES need to be given!");}
# }
# elsif ($mode == 2) {  ## or mode =~ /ggfile/
# 	if (defined $former_run_dir) {
# 		&constructDirectory($starttime,$former_run_dir);
# 		&read_ggfile($genome_gene_file);
# 	} else {dieWithUnexpectedError("In Mode 2, FORMER RUN DIRECTORY needs to be given!");}
# }
# elsif ($mode == 3) {  ## or mode =~ /blastparse/
# 	if ((defined $usr_blast_file) && (defined $usr_gg_file)) {
# 		&constructDirectory($starttime);
# 		$all_fa_file      = 'N/A';
# 		$genome_gene_file = $usr_gg_file;
# 		read_ggfile($genome_gene_file);
# 		$blast_file       = $usr_blast_file;
# 		&blast_parse($blast_file,$bpo_file,$pv_cutoff,{},\%blast_flag) unless (-e $bpo_file);
# 	} else {dieWithUnexpectedError("In Mode 3, BLAST OUT FILE and GENOME-GENE FILE are required!");}
# }
# elsif(...)

if ($mode =~ /^4/ or $mode =~ /^par/) {  ## or mode =~ /inputbpo/
  $mode=~s/par/parallel/ unless($mode=~/parallel/);
	if ((defined $usr_bpo_file) && (defined $usr_gg_file)) {
		
		if($former_run_dir) { # FIXME: need construct($former_run_dir) at par_start ..
		  reuseDirectory($starttime,$former_run_dir); #   need param to set ORTHOMCL_WORKING_DIR name
		} else {  
		  &constructDirectory($starttime,$former_run_dir); # par_start only, need param to set ORTHOMCL_WORKING_DIR name
		}
		
		$all_fa_file      = 'N/A';
		$genome_gene_file = $usr_gg_file;
		$blast_file       = 'N/A';
		read_ggfile($genome_gene_file);  #dgg.parallel: need for all steps
		$bpo_file         = $usr_bpo_file;
		if ($usr_bpo_file =~ m/(\S+)\.(\S+)/) {
			$bpo_idx_file     = $1.'_bpo.idx';
			$bpo_se_file      = $1.'_bpo.se';
		} else {
			$bpo_idx_file     = $usr_bpo_file.'_bpo.idx';
			$bpo_se_file      = $usr_bpo_file.'_bpo.se';
		}
	} else {dieWithUnexpectedError("In Mode 4, BPO (BLAST PARSE OUT) FILE and GG (GENOME-GENE RELATION) FILE are required!");}
} 
else {dieWithUnexpectedError("Mode 4 or Mode parallel needs to be given!");}

# elsif ($mode == 5) {  ## or mode =~ /mclrun/
# 	if ((defined $former_run_dir) && (defined $usr_taxa_file)) {
# 		mode5($starttime,$command,$former_run_dir,$usr_taxa_file,$inflation);
# 		my $endtime = `date`;
# 		&write_endtime_in_parameter_log($endtime);
# 		write_log("\nStart Time: $starttime\nEnd Time:   $endtime\n");
# 		warn "\nStart Time: $starttime\nEnd Time:   $endtime\n"; exit(0); # dgg: this should be exit(0) success, not die
# 		##die "\nStart Time: $starttime\nEnd Time:   $endtime\n"; # dgg: this should be exit(0) success, not die
# 
# 	} else {dieWithUnexpectedError("In Mode 5, FORMER RUN DIR and TAXA LIST FILE are required!");}
# } 
# else {dieWithUnexpectedError("Mode 1,2,3,4 or 5 needs to be given!");}

sub MAIN_START{} # dgg

unless($mode =~ /parallel/) {
&write_parameter_log($starttime,$command,$mode,$pv_cutoff,$pi_cutoff,$pmatch_cutoff,$inflation,$maximum_weight);
}

#dgg: step1. 1cpu as before
&constructIDX_for_bpofile($bpo_file,$bpo_idx_file) unless (-e $bpo_idx_file);
&constructSE_for_bpofile($bpo_file,$bpo_se_file) unless (-e $bpo_se_file);

## dgg.parallel: need open_bpo/retrieve for all steps
&open_bpofile($bpo_file);
&retrieve_from_file($bpo_idx_file,$bpo_se_file);

#........... redo main steps .................
# parallel driver script:
#   orthomcl_evg.pl -mode par_start; 
#   i=0; while [ $i -lt $ntaxa ]; do { orthomcl_evg.pl -mode par_part$i &; i=$(( $i + 1 )); } done; wait; ## parallel part
#   orthomcl_evg.pl -mode par_end;

if($mode =~ /parallel/) {

  if($mode =~ /start|first/) { 
    # exit? first step = constructIDX_for_bpofile
    &write_parameter_log($starttime,$command,$mode,$pv_cutoff,$pi_cutoff,$pmatch_cutoff,$inflation,$maximum_weight);
    my $endtime = `date`;
   # write_log("\nMODE:$mode; Start: $starttime, End: $endtime\n");
    exit;
    
  } elsif($mode =~ /part(\d+)/) {  # -mode=parallel_part0 .. part10 for 11 species
    my $itaxa=$1; 
    makepart_Inparalog($itaxa);
    makepart_Ortholog($itaxa);
    
    # NO: makepart_Connect($itaxa); ## move here? cant do, need all ortholog.xxx-*.bbh for connect

    my $endtime = `date`;
    # write_log("\nMODE:$mode; Start: $starttime, End: $endtime\n");
    exit;
    

  } elsif($mode =~ /connect(\d+)/) {  # -mode=par_conn0 .. conn10 for 11 species
    my $itaxa=$1; 
    readall_Inparalog();  # these two fill up %connect{sppA_sppB} hash
    readall_Ortholog();   # for makepart_conn dont need readall_ here ? but simpler
    makepart_Connect($itaxa);
    my $endtime = `date`;
    exit;

  } elsif($mode =~ /finish|last|end/) {  # -mode=parallel_combine / _finish / _connect / _end
    ## ** THIS TAKES LOONNNGGGG TIME now, for kfish2.  Probable step where time zoomed from 11hr to 48hr, 1cpu
    ## Why? any par fix?  Mem use is small, 8 GB on gordon.sdsc.  Steps start, part0..10 took only 1hr (16cpu).
    write_log("MODE $mode at: ".`date`."\n");
    readall_Inparalog();  # these two fill up %connect{sppA_sppB} hash
    readall_Ortholog();   # 
    
    write_log("connect_OrthoParalogs at: ".`date`."\n");
    ## ADD/Replace here makepart_Connect > readall_Connect
    my @spconn= readall_Connect();
    if(@spconn==0) {  connect_OrthoParalogs(); } # error if nconn < ntaxa ?
    elsif(scalar(@spconn) < scalar(@taxa) - 1) { warn "ERR: readall_connect: got only @spconn of @taxa\n"; } # what? die?
    
    write_log("COMPUTE_MCL_MATRIX at: ".`date`."\n");
    COMPUTE_MCL_MATRIX();
  }
  
} else { # original steps
  makeall_Inparalog();
  makeall_Ortholog();
  COMPUTE_MCL_MATRIX();
}
## dgg: MAIN ends here now.

#.......... parallel parts ........................

our %bbhpart_fh=();
our $CUR_PARTBBH="all_blast.bbh"; # default $bbh_file in orthomcl_modevg.pm

sub open_partbbh {
	my($bbhpart,$dowrite) = @_;
	return 1 if(ref $bbhpart_fh{$bbhpart});
  my $bbh_file=$ORTHOMCL_WORKING_DIR.$bbhpart;
  my ($fh,$ok)=(undef,0);
  if($dowrite and $dowrite =~ /^w/) { $ok=open ($fh,">",$bbh_file); }
  else { $ok=open ($fh, $bbh_file); }
  $bbhpart_fh{$bbhpart}= $fh;
  # $CUR_PARTBBH=$bbhpart; # always set this? or not?
  return $ok;
}

sub close_partbbh {
	my($bbhpart) = @_;
  if(ref $bbhpart_fh{$bbhpart}) {
    my $fh= delete $bbhpart_fh{$bbhpart}; close($fh);
  }
}

sub write_partbbh {
	my($bbhpart,$data) = @_;
	# open_partbbh($bbhpart,"write") unless(ref $bbhpart_fh{$bbhpart});
	my $fh= $bbhpart_fh{$bbhpart};
	print $fh $data;
}

sub makepart_Inparalog {
	my($itaxa) = @_; # my($taxon) = @_;
  my $itaxname= $taxa[$itaxa];
  $CUR_PARTBBH="inparalog.$itaxname.bbh";
  write_log("makepart_Inparalog: $CUR_PARTBBH\n");  
  open_partbbh($CUR_PARTBBH,"write");
  my @matrix= makeInparalog($itaxname); ## return(\%edge, \%weight, $avgw, $sumw, $count);
}

sub readall_Inparalog {
  foreach my $taxon (@taxa) {
    my $partbbh="inparalog.$taxon.bbh";
    write_log("readall_Inparalog: $partbbh\n");  
	  @{$connect{$taxon.' '.$taxon}}  = read_bbhmatrix($partbbh);  
  }
}

# matrix: write_partbbh($CUR_PARTBBH,"$query\t$subject\t".$pvalue{$query.' '.$subject}."	".$pvalue{$subject.' '.$query}."\n");
# matrix:	return (\%edge, \%weight, $avgw, $sumw, $count);

sub makepart_Ortholog {
	my($itaxa) = @_;
  my $itaxname= $taxa[$itaxa];
  # for(my $i=0;$i<scalar(@taxa)-1;$i++) { makepart_Ortholog($i); }
	for(my $j=$itaxa+1;$j<scalar(@taxa);$j++) {
	  my $jtaxname= $taxa[$j];
    my $partbbh="ortholog.$itaxname-$jtaxname.bbh"; # put all jtax in itax.bbh
	  write_log("makepart_Ortholog: $partbbh\n");  
    $CUR_PARTBBH=$partbbh; # put all jtax in itax.bbh ?? no, not for readall_
    open_partbbh($partbbh,"write");
    my @matrix= makeOrtholog($itaxname,$jtaxname); ## (\%edge, \%weight, $avgw, $sumw, $count);
    ##  @{$connect{$taxa[$i].' '.$taxa[$j]}}  = @matrux
    }
}

sub readall_Ortholog {
  my $ntaxa=scalar(@taxa);
  for(my $i=0; $i<$ntaxa-1; $i++) {
    my $itaxname= $taxa[$i];
    for(my $j=$i+1; $j<$ntaxa; $j++) {
      my $jtaxname= $taxa[$j];
      my $partbbh="ortholog.$itaxname-$jtaxname.bbh"; # put all jtax in itax.bbh
      write_log("readall_Ortholog: $partbbh\n");  
      @{$connect{$taxa[$i].' '.$taxa[$j]}}  = read_bbhmatrix($partbbh);  
      }
    }
}


sub readall_Connect {
  my $ntaxa=scalar(@taxa); my $nok=0; my @sok=();
  for(my $i=0; $i<$ntaxa-1; $i++) {
    my $itaxname= $taxa[$i];
    my $nread= read_connpart($i);
    do{ push @sok,$itaxname; $nok++; } if($nread>0);
     ## does: @{$connect{$taxai.' '.$taxaj}}  = (\%edge1, \%weight1, $avgw, $sumw, $count);
  }
  return (wantarray)?@sok:$nok;
}

=item makepart_Connect connect.slow.for.smallest/last  

  THIS step really SLOW for LAST itaxa 1 pair: tilapia-zfish, should be fastest there ????
  what bug is this? .. is it slow due to spp sort order? tilapia-zfish are both last
  blastqueryab() calls to disk .bpo only likely problem.  
  .. why last 2 tilapia-zfish spp in gene lists?
  ** ~12hr runtime for par_connect ..
  
First gene.bpo indices for each spp x zfish:
-- bpo needs to have all 2nd IDs together for 1st id, but 1st IDs can be random?

  9;catfish_IctpunEGm021478t1;129;zfish_ENSDARP00000056560;99;6e-58;89;1:33-129:1-97.
  404140;human_UniRef50_P62258;255;zfish_ENSDARP00000032575;255;2e-180;96;1:1-255:1-255.
  2065590;kfish2_Funhe2EKm032398t1;1049;zfish_ENSDARP00000073993;1207;0.0;51;1:19-1031:10-1074.
  4000953;mayzebr_XP_004538613.1;478;zfish_ENSDARP00000083712;881;0.0;65;1:1-465:1-474.
  5363243;medaka_ENSORLP00000000001;441;zfish_ENSDARP00000109103;300;1e-09;25;1:12-202:10-226.
  6630518;platyfish_ENSXMAP00000000003;214;zfish_ENSDARP00000032941;213;7e-59;56;1:1-214:1-213.
  7701479;spotgar_ENSGACP00000014531_1;293;zfish_ENSDARP00000006025;298;3e-178;86;1:1-290:1-290.
  9038032;stickleback_ENSGACP00000000002;525;zfish_ENSDARP00000080637;611;0.0;67;1:1-474:44-522.
  10362994;tetraodon_ENSTNIP00000006514;226;zfish_ENSDARP00000059027;235;4e-88;57;1:1-226:1-235.
  11436607;tilapia_ENSONIP00000000001;893;zfish_ENSDARP00000090596;1013;0.0;59;1:1-872:1-1000.
    ^^ lineno in bpo
  1;catfish_IctpunEGm021478t1;129;catfish_IctpunEGm021478t1;129;2e-91;100;1:1-129:1-129.
    ^^ 1st line
  15351476;catfish_IctpunEGm021477t1;564;human_UniRef50_H0YL91;368;3e-138;60;1:3-336:10-347.
    ^^ last line, nlines=15351476 fish11c_omcl.bpo
    why is catfish 1st and last? probably due to blastp input order, from split-cat, somewhat random.

** ortholog.*-zfish.bbh is largest for slowest; 
    1973674 refish11c/ortholog.catfish-zfish.bbh
    1428514 refish11c/ortholog.human-zfish.bbh
    2171649 refish11c/ortholog.kfish2-zfish.bbh       * 5th slowest
    2452609 refish11c/ortholog.mayzebr-zfish.bbh      * 4th slowest
    1641055 refish11c/ortholog.medaka-zfish.bbh       + 3rd fast
    1774085 refish11c/ortholog.platyfish-zfish.bbh    + 2nd fastest
    2278025 refish11c/ortholog.spotgar-zfish.bbh      * 2nd slowest
    1926218 refish11c/ortholog.stickleback-zfish.bbh  + 1st fastest
    1697235 refish11c/ortholog.tetraodon-zfish.bbh    + 4th fast
    3053540 refish11c/ortholog.tilapia-zfish.bbh      * slowest && largest bbh, 1.5x 2x others
    
-- ortholog*-stickback.bbh has 1.3M-1.9M size
   ortholog*-tetraodon.bbh has 1.2M-1.7M size
-- however gene count for tilapia,zfish (main.aa) is in same range as others, 20k-25k

-- ortholog size is problem, can reduce w/ these filters:
      pmatch_cutoff: makebpo already has -identmin=33 align-identity cut, equal to pmatch_cutoff here
      -- boost blast2omcl -identmin=40  likely would speed up much.. maybe used that for kfish1, 12hr runs
      
			next if($bla[1].'e'.$bla[2] > $pv_cutoff || $bla[3]< $pi_cutoff);
		  next if(&simspan($bla[4]) < $pmatch_cutoff); 
      * pmatch_cutoff will reduce ortholog counts, removing tiny aligns, as some added compute cost
        -- uses simspan() that reads from disk bpo.
        
        
  ## rev-order output .. not directly due to size of j-taxa list
  connect.tilapia.con : largest ortholog.til-zfish.bbh by 1.5/2x causes this slow
  -- 4hr step here --  
  connect.spotgar.con   stickleback tetraodon tilapia zfish
  -- 3hr step here --
  connect.human.con	    kfish2 mayzebr medaka platyfish spotgar stickleback tetraodon tilapia zfish
  connect.mayzebr.con	  medaka platyfish spotgar stickleback tetraodon tilapia zfish
  connect.kfish2.con	  mayzebr medaka platyfish spotgar stickleback tetraodon tilapia zfish
  connect.catfish.con	  human kfish2 mayzebr platyfish spotgar stickleback tetraodon tilapia zfish
  -- 4hr step here -- vv finish in ~1hr
  connect.tetraodon.con	tilapia zfish
  connect.medaka.con	  platyfish spotgar stickleback tetraodon tilapia zfish
  connect.platyfish.con	spotgar stickleback tetraodon tilapia zfish
  connect.stickleback.con	tetraodon tilapia zfish
  connect.zfish.con	== zero

  ## still no .con after 7hr run
  tilapia-zfish : 10hr still running .. what the f*ck?
  spotgar-stickleback tetraodon tilapia  zfish : done after 8hr

=cut
  
sub makepart_Connect {
	my($itaxa) = @_;
  my $ntaxa=scalar(@taxa);
  my $itaxname= $taxa[$itaxa];

  my $i=$itaxa; 
	for(my $j=$itaxa+1;$j<$ntaxa;$j++) {
		write_log("makepart_Connect co-ortholog pairs between $taxa[$i] and $taxa[$j]: ");  
		my $c_coortholog=0;
     ## connect{i-j}== (\%edge, \%weight, $avgw, $sumw, $count);
     
		my %e;
		my $edge_ref=$connect{$taxa[$i].' '.$taxa[$j]}->[0];
		foreach my $pi (keys %$edge_ref) { @{$e{$pi}}=@{$edge_ref->{$pi}}; }  #make a copy of current edge data structure into %e

		my %w =  %{$connect{$taxa[$i].' '.$taxa[$j]}->[1]};
		my $sumw =  $connect{$taxa[$i].' '.$taxa[$j]}->[3];
		my $c_ortholog = $connect{$taxa[$i].' '.$taxa[$j]}->[4];
		my %p1 = %{$connect{$taxa[$i].' '.$taxa[$i]}->[0]};
		my %p2 = %{$connect{$taxa[$j].' '.$taxa[$j]}->[0]};
		my %para;
		foreach my $p (keys %{$connect{$taxa[$i].' '.$taxa[$i]}->[0]}) {
			push (@{$para{$p}}, @{$p1{$p}}); }
		foreach my $p (keys %{$connect{$taxa[$j].' '.$taxa[$j]}->[0]}) {
			push (@{$para{$p}}, @{$p2{$p}}); }
		
		foreach my $n (keys %e) {
			$ortho{$n} = 1;
			my (@nodes1, @nodes2);

			if (exists($para{$n})) {push (@nodes1, $n, @{$para{$n}});}
			else {push (@nodes1, $n);}

			foreach (@{$e{$n}}) {
				if (exists($para{$_})) {push (@nodes2, $_, @{$para{$_}});}
				else {push (@nodes2, $_);}
			}

			@nodes1=@{nonredundant_list(\@nodes1)}; #can be commented
			@nodes2=@{nonredundant_list(\@nodes2)};

			for(my $k=0;$k<scalar(@nodes1);$k++) {
				for(my $l=0;$l<scalar(@nodes2);$l++) {
					
					next if(exists($w{$nodes1[$k].' '.$nodes2[$l]}));
					my ($pv1, $pv2);

          my @blab= blastqueryab($nodes1[$k],$nodes2[$l]);
					if (@blab > 4) {
						my ($s,$pm,$pe,$pi)=@blab[0,3,4,5];
		        # blab == (&getline_from_bpofile($i))[0,1,3,5,6,7]; == simid;qid;//qlen;tid;//tlen;5:eval-pm;eval-pe;pid;//offs
						next if($pm.'e'.$pe > $pv_cutoff || $pi< $pi_cutoff);
						if($pmatch_cutoff) {
							next if(&simspan($s) < $pmatch_cutoff);
						}
						if($pm==0) { $pv1 = $maximum_weight;} else { $pv1 = -log($pm.'e'.$pe)/log(10); }
					} else {next;}

          @blab= blastqueryab($nodes2[$l],$nodes1[$k]);
					if (@blab > 4) {
						my ($s,$pm,$pe,$pi)=@blab[0,3,4,5];
						next if($pm.'e'.$pe > $pv_cutoff || $pi < $pi_cutoff);
						if($pmatch_cutoff) {
							next if(&simspan($s) < $pmatch_cutoff);
						}
						if($pm==0) { $pv2 = $maximum_weight;} else { $pv2 = -log($pm.'e'.$pe)/log(10); }
						
						push (@{$edge_ref->{$nodes1[$k]}}, $nodes2[$l]);
						push (@{$edge_ref->{$nodes2[$l]}}, $nodes1[$k]);
						my $wt = ($pv1+$pv2)/2;
						# use averaged score as edge weight
						$sumw += $wt;
						$c_coortholog++;
						my $wtf= sprintf("%.3f", $wt);
						$w{$nodes1[$k].' '.$nodes2[$l]} = $w{$nodes2[$l].' '.$nodes1[$k]} = $wtf;
					}
				}
			}
		}
	  write_log("$c_coortholog pairs\n");
		my $avgw = 'N/A';
		if ($c_ortholog+$c_coortholog) {
			$avgw = $sumw/($c_ortholog+$c_coortholog);
		  foreach my $p (keys %w) { $w{$p} = sprintf("%.3f", $w{$p}/$avgw); }
		} else { # never here?
		  my $AVW= $sumw || 1;
		  foreach my $p (keys %w) { $w{$p} = sprintf("%.3f", $w{$p}/$AVW); }
		}
		write_log("$taxa[$i] and $taxa[$j] average weight: $avgw\n");
		$connect{$taxa[$i].' '.$taxa[$j]}->[1] = \%w;
  }
  write_connpart($itaxa);
}



sub connect_OrthoParalogs {
  my $ntaxa=scalar(@taxa);

## THIS step taking hours vs prior, can it be parallelzd ?? over i-taxa?
## maybe, step $i doesnt refer to i - 1; need to write slices of matrix connect[ itax - jtax ] then reread
## %connect == { i x j } (\%edge, \%weight, $avgw, $sumw, $count) PLUS add below
## this is only change to connect? $connect{$taxa[$i].' '.$taxa[$j]}->[1] = \%w;

  for(my $i=0; $i<$ntaxa-1; $i++) {
	 for(my $j=$i+1; $j<$ntaxa; $j++) {

    ## done already..
		# write_log("\nIdentifying ortholog pairs between $taxa[$i] and $taxa[$j]\n");  
		# @{$connect{$taxa[$i].' '.$taxa[$j]}} = &makeOrtholog($taxa[$i],$taxa[$j]); # identification of orthologs

		write_log("Appending co-ortholog pairs between $taxa[$i] and $taxa[$j]: ");  
		my $c_coortholog=0;

		my %e;
		my $edge_ref=$connect{$taxa[$i].' '.$taxa[$j]}->[0];
		foreach my $pi (keys %$edge_ref) { @{$e{$pi}}=@{$edge_ref->{$pi}}; }  #make a copy of current edge data structure into %e

		my %w =  %{$connect{$taxa[$i].' '.$taxa[$j]}->[1]};
		my $sumw =  $connect{$taxa[$i].' '.$taxa[$j]}->[3];
		my $c_ortholog = $connect{$taxa[$i].' '.$taxa[$j]}->[4];
		my %p1 = %{$connect{$taxa[$i].' '.$taxa[$i]}->[0]};
		my %p2 = %{$connect{$taxa[$j].' '.$taxa[$j]}->[0]};
		my %para;
		foreach my $p (keys %{$connect{$taxa[$i].' '.$taxa[$i]}->[0]}) {
			push (@{$para{$p}}, @{$p1{$p}}); }
		foreach my $p (keys %{$connect{$taxa[$j].' '.$taxa[$j]}->[0]}) {
			push (@{$para{$p}}, @{$p2{$p}}); }
		
		foreach my $n (keys %e) {
			$ortho{$n} = 1;
			my (@nodes1, @nodes2);

			if (exists($para{$n})) {push (@nodes1, $n, @{$para{$n}});}
			else {push (@nodes1, $n);}

			foreach (@{$e{$n}}) {
				if (exists($para{$_})) {push (@nodes2, $_, @{$para{$_}});}
				else {push (@nodes2, $_);}
			}

			@nodes1=@{nonredundant_list(\@nodes1)}; #can be commented
			@nodes2=@{nonredundant_list(\@nodes2)};

			for(my $k=0;$k<scalar(@nodes1);$k++) {
				for(my $l=0;$l<scalar(@nodes2);$l++) {
					
					next if(exists($w{$nodes1[$k].' '.$nodes2[$l]}));
					my ($pv1, $pv2);

          my @blab= blastqueryab($nodes1[$k],$nodes2[$l]);
					if (@blab > 4) {
						my ($s,$pm,$pe,$pi)=@blab[0,3,4,5];
		        # blab == (&getline_from_bpofile($i))[0,1,3,5,6,7]; == simid;qid;//qlen;tid;//tlen;5:eval-pm;eval-pe;pid;//offs
						next if($pm.'e'.$pe > $pv_cutoff || $pi< $pi_cutoff);
						if($pmatch_cutoff) {
							next if(&simspan($s) < $pmatch_cutoff);
						}
						if($pm==0) { $pv1 = $maximum_weight;} else { $pv1 = -log($pm.'e'.$pe)/log(10); }
					} else {next;}

          @blab= blastqueryab($nodes2[$l],$nodes1[$k]);
					if (@blab > 4) {
						my ($s,$pm,$pe,$pi)=@blab[0,3,4,5];
						next if($pm.'e'.$pe > $pv_cutoff || $pi < $pi_cutoff);
						if($pmatch_cutoff) {
							next if(&simspan($s) < $pmatch_cutoff);
						}
						if($pm==0) { $pv2 = $maximum_weight;} else { $pv2 = -log($pm.'e'.$pe)/log(10); }
						
						push (@{$edge_ref->{$nodes1[$k]}}, $nodes2[$l]);
						push (@{$edge_ref->{$nodes2[$l]}}, $nodes1[$k]);
						my $wt = ($pv1+$pv2)/2;
						# use averaged score as edge weight
						$sumw += $wt;
						$c_coortholog++;
						my $wtf= sprintf("%.3f", $wt);
						$w{$nodes1[$k].' '.$nodes2[$l]} = $w{$nodes2[$l].' '.$nodes1[$k]} = $wtf;
					}
				}
			}
		}
	  write_log("$c_coortholog pairs\n");
		my $avgw = 'N/A';
		if ($c_ortholog+$c_coortholog) {
			$avgw = $sumw/($c_ortholog+$c_coortholog);
		  foreach my $p (keys %w) { $w{$p} = sprintf("%.3f", $w{$p}/$avgw); }
		} else { # never here?
		  my $AVW= $sumw || 1;
		  foreach my $p (keys %w) { $w{$p} = sprintf("%.3f", $w{$p}/$AVW); }
		}
		write_log("$taxa[$i] and $taxa[$j] average weight: $avgw\n");
		#bad# foreach my $p (keys %w) { $w{$p} = sprintf("%.3f", $w{$p}/$avgw); }
		$connect{$taxa[$i].' '.$taxa[$j]}->[1] = \%w;
	}
}

} # dgg: connect_Ortholog


#..................................................
#.......... original parts ........................

#dgg: step2. n-taxa cpu change, write inparalog.$taxon.bbh for each taxon
# -- table is bbh_file? result of sub matrix()

sub makeall_Inparalog {
open_partbbh($CUR_PARTBBH,"write"); # all_blast.bbh here
foreach my $taxon (@taxa) {
	write_log("\nIdentifying inparalogs from $taxon\n");
	@{$connect{$taxon.' '.$taxon}}  = &makeInparalog($taxon);                      # identification of inparalogs
}

}

#dgg: step3. n-taxa cpu change, write ortholog.$taxon.bbh for each taxon
# -- table is bbh_file? result of sub matrix()
# dgg: this could be parallel loop; each taxa x taxa independent test ??
# dgg: died silently here during 250k gene, 13 taxa run (1.6GB .bpo file)

sub makeall_Ortholog {

for(my $i=0;$i<scalar(@taxa)-1;$i++) {
	for(my $j=$i+1;$j<scalar(@taxa);$j++) {
		write_log("\nIdentifying ortholog pairs between $taxa[$i] and $taxa[$j]\n");  

		@{$connect{$taxa[$i].' '.$taxa[$j]}} = &makeOrtholog($taxa[$i],$taxa[$j]); # identification of orthologs

		write_log("Appending co-ortholog pairs between $taxa[$i] and $taxa[$j]: ");  
		my $c_coortholog=0;

		my %e;
		my $edge_ref=$connect{$taxa[$i].' '.$taxa[$j]}->[0];
		foreach my $pi (keys %$edge_ref) {@{$e{$pi}}=@{$edge_ref->{$pi}};}  #make a copy of current edge data structure into %e

		my %w =  %{$connect{$taxa[$i].' '.$taxa[$j]}->[1]};
		my $sumw =  $connect{$taxa[$i].' '.$taxa[$j]}->[3];
		my $c_ortholog = $connect{$taxa[$i].' '.$taxa[$j]}->[4];
		my %p1 = %{$connect{$taxa[$i].' '.$taxa[$i]}->[0]};
		my %p2 = %{$connect{$taxa[$j].' '.$taxa[$j]}->[0]};
		my %para;
		foreach my $p (keys %{$connect{$taxa[$i].' '.$taxa[$i]}->[0]}) {
			push (@{$para{$p}}, @{$p1{$p}}); }
		foreach my $p (keys %{$connect{$taxa[$j].' '.$taxa[$j]}->[0]}) {
			push (@{$para{$p}}, @{$p2{$p}}); }

		
		foreach my $n (keys %e) {
			$ortho{$n} = 1;
			my (@nodes1, @nodes2);

			if (exists($para{$n})) {push (@nodes1, $n, @{$para{$n}});}
			else {push (@nodes1, $n);}

			foreach (@{$e{$n}}) {
				if (exists($para{$_})) {push (@nodes2, $_, @{$para{$_}});}
				else {push (@nodes2, $_);}
			}

			@nodes1=@{nonredundant_list(\@nodes1)}; #can be commented
			@nodes2=@{nonredundant_list(\@nodes2)};

			for(my $k=0;$k<scalar(@nodes1);$k++) {
				for(my $l=0;$l<scalar(@nodes2);$l++) {
					
					next if(exists($w{$nodes1[$k].' '.$nodes2[$l]}));
					my ($pv1, $pv2);

					if (blastqueryab($nodes1[$k],$nodes2[$l])) {
						my ($s,$pm,$pe,$pi)=(blastqueryab($nodes1[$k],$nodes2[$l]))[0,3,4,5];
						next if($pm.'e'.$pe > $pv_cutoff || $pi< $pi_cutoff);
						if($pmatch_cutoff) {
							next if(&simspan($s) < $pmatch_cutoff);
						}
						if($pm==0) { $pv1 = $maximum_weight;} else { $pv1 = -log($pm.'e'.$pe)/log(10); }
					} else {next;}

					if (blastqueryab($nodes2[$l],$nodes1[$k])) {
						my ($s,$pm,$pe,$pi)=(blastqueryab($nodes2[$l],$nodes1[$k]))[0,3,4,5];
						next if($pm.'e'.$pe > $pv_cutoff || $pi < $pi_cutoff);
						if($pmatch_cutoff) {
							next if(&simspan($s) < $pmatch_cutoff);
						}
						if($pm==0) { $pv2 = $maximum_weight;} else { $pv2 = -log($pm.'e'.$pe)/log(10); }
						push (@{$edge_ref->{$nodes1[$k]}}, $nodes2[$l]);
						push (@{$edge_ref->{$nodes2[$l]}}, $nodes1[$k]);
						my $wt = ($pv1+$pv2)/2;
						# use averaged score as edge weight
						$w{$nodes1[$k].' '.$nodes2[$l]} = sprintf("%.3f", $wt);
						$w{$nodes2[$l].' '.$nodes1[$k]} = sprintf("%.3f", $wt);
						$sumw += $wt;
						$c_coortholog++;
					}
				}
			}
		}
		write_log("$c_coortholog pairs\n");
		my $avgw = 'N/A';
		if ($c_ortholog+$c_coortholog) {
			$avgw = $sumw/($c_ortholog+$c_coortholog);
		}
		write_log("$taxa[$i] and $taxa[$j] average weight: $avgw\n");
		foreach my $p (keys %w) {
			$w{$p} = sprintf("%.3f", $w{$p}/$avgw);
		}
		$connect{$taxa[$i].' '.$taxa[$j]}->[1] = \%w;
	}
}

} # dgg: makeall_Ortholog


#........... main continues here ............................

sub COMPUTE_MCL_MATRIX
{
# evg par problem: all_ortho.mtx : zeros are missing from weight{ q x s }

%blastquery=();
%gindex=();

foreach my $taxon (@taxa) {
	write_log("\ncalculate average weight from $taxon\n");
	my %e = %{$connect{$taxon.' '.$taxon}->[0]};
	my %w = %{$connect{$taxon.' '.$taxon}->[1]};

	my $count=0; my $sum=0;
	my $count_all=0; my $sum_all = 0;
	foreach my $pair (keys %w) {
		my ($n,$p) = split(' ',$pair);
		$count_all++; $sum_all += $w{$n.' '.$p};
		if ($ortho{$n} || $ortho{$p}) {
			$count++;
			$sum += $w{$n.' '.$p};
		}
	}
	my $avgw;
	# normalize the in-paralog weights by the average weight of inparalogs which have orthologs in other species
	# common case, for eukaryotes and most prokaryotes
	if ($count) {
		$avgw = $sum/$count;
	}
	# OR normalize the in-paralog weights by the average weight of all inparalogs
	# not common, useful for prokaryotes or pathogens
	elsif ($count_all) {
		$avgw = $sum_all/$count_all;
		write_log("taxon average weight is calculated based on all inparalog pairs\n");
	}
	# OR no normalization since $count_all=0 and there is nothing stored in %weight
	# not common, useful for prokaryotes or pathogens 
	else {
		$avgw = 'N/A';
		write_log("taxon average weight is not calculated because there's no inparalog pairs\n");
	}
	write_log("$taxon average weight: $avgw\n");
	if($avgw>0) {
	foreach my $p (keys %w) {
		$w{$p} = sprintf("%.3f", $w{$p}/$avgw); ## dgg: fails here for avew==0
	} }
	$connect{$taxon.' '.$taxon}->[1] = \%w; 
}
%ortho=();


foreach my $p (keys %connect) {
	my %e = %{$connect{$p}->[0]};
	my %w =  %{$connect{$p}->[1]};
	
	foreach my $n (keys %e) {
		push(@{$graph{$n}}, @{$e{$n}});
		delete $e{$n};
	}
	%e=();
	foreach my $n (keys %w) {
		$weight{$n} = $w{$n};
		delete $w{$n};
	}
	%w=();
	delete $connect{$p};
}
%connect=();

write_matrix_index($matrix_file,$index_file);

%graph=();
%weight=();

executeMCL($matrix_file,$mcl_file,$inflation);
mcl_backindex($mcl_file,$mcl_bi_file);
%gindex2=();

my $endtime = `date`;
write_log("\nStart Time: $starttime\nEnd Time:   $endtime\n");
&write_endtime_in_parameter_log($endtime);

} # end COMPUTE_MCL_MATRIX / main


#######################################SUBROUTINES###########################################

# This subroutine is an important part of OrthoMCL, used to
# look for inparalog (recent paralog) which is defined as 
# reciprocal better hits here.
# Please refer to the OrthoMCL paper for more details.
# One Arguments:
# 1. String Variable: Taxon name
# Last modified: 10/02/06

sub makeInparalog {
	my $taxon = $_[0];
	my (%seqs, %inbest, %pvalue,%sim);
	foreach (@{$gindex{$taxon}}) {$seqs{$_} = 1;}
	foreach my $qid (keys %seqs) {
		my ($sStart,$sEnd);
		if (defined $blastquery{$qid}) {
			($sStart,$sEnd)=split (";",$blastquery{$qid});
		} else {next;}
		
		## FIXME dgg : inefficient here, pvtie_sort() reads same getline_bpo as below..
		## combine both here, read once; pvtie_sort used only here
		## FIXME2: inparalogs computed here but not reported.  comp is sort order by evalue (pm-pe)
		##   but not as good as bitscore + align-size/pct-align .. comp both here, use bits-align as better para/orth score
		
		my @sorted_simid=pvtie_sort($sStart,$sEnd,$taxon);
		LINE:foreach (0..$#sorted_simid) {
			my ($s,$sid,$pm,$pe,$pi)=(&getline_from_bpofile($sorted_simid[$_]))[0,3,5,6,7];
			if ($sid ne $qid) {
				last LINE unless ($seqs{$sid});                                  ## better hit not from the same species
				last LINE if($pm.'e'.$pe > $pv_cutoff || $pi < $pi_cutoff);      ## better hit not meet the cutoff
				if($pmatch_cutoff) {
					next LINE if(&simspan($s) < $pmatch_cutoff);
				}
				push(@{$inbest{$qid}}, $sid);
				$pvalue{$qid.' '.$sid} = $pm.'e'.$pe; 
			}
		}
	}
	my @b = keys %inbest;
	write_log(scalar(@b)." sequences have better hits within species\n");
	return &matrix(\%inbest, \%pvalue);
} ##makeInparalog


# This subroutine is used by the subroutine makeInparalog,
# to solve the pv_tie problem. It rearranges the pv-tied blast hits
# so that the hits from a specific taxon are moved higher than hits from  other species.
# Three Arguments:
# 1. Number Variable: starting line id (or similarity id) of bpo file 
# 2. Number Variable: ending line id (or similarity id) of bpo file 
# 3. String Variable: taxon
# Last modified: 07/20/04, then by dgg...

sub pvtie_sort {
	my ($sStart,$sEnd,$taxon)=@_;
	my (@sorted_simid,@tmp);
	my ($lastpm,$lastpe)=(getline_from_bpofile($sStart))[5,6];
	foreach ($sStart..$sEnd) {
	  local $^W = 0; #dgg no warnings missing $gindex2{sid}
		my ($s,$sid,$pm,$pe,$pi)=(&getline_from_bpofile($_))[0,3,5,6,7];
		
		if (($lastpm==$pm) && ($lastpe==$pe)) {
			if ($gindex2{$sid} eq $taxon) { push (@sorted_simid,$s); } 
			else { push (@tmp,$s); }
		} else {
			if (scalar(@tmp)>0) { push (@sorted_simid,@tmp); @tmp=(); } 
			if ($gindex2{$sid} eq $taxon) { push (@sorted_simid,$s); } 
			else { push (@tmp,$s); }
		}
		$lastpm=$pm;$lastpe=$pe;
	}
	if (scalar(@tmp)>0) { push (@sorted_simid,@tmp);}
	return @sorted_simid;
} ## pvtie_sort



# This subroutine is an important part of OrthoMCL, used to
# look for ortholog which is defined as the reciprocal best
# hit between two species.
# Please refer to the OrthoMCL paper for more details.
# Two Arguments:
# 1. String Variable: Taxon name
# 2. String Variable: Taxon name
# Last modified: 10/02/06
sub makeOrtholog {
	my ($ta,$tb) = @_;
	my (@seqs,%best,%sim,%pvalue);
	foreach my $qid (@{$gindex{$ta}}) {
		my ($sStart,$sEnd);
		if (defined $blastquery{$qid}) {
			($sStart,$sEnd)=split (";",$blastquery{$qid});
		} else {next;}
		my ($lastpm,$lastpe);
		my $hit_id=0;
		LINE:foreach ($sStart..$sEnd) {
			my ($s,$sid,$pm,$pe,$pi)=(&getline_from_bpofile($_))[0,3,5,6,7];
			if (defined $gindex2{$sid}) {
				if ($gindex2{$sid} eq $tb) {
					$hit_id++;
					if ($hit_id==1) {
						push(@{$sim{$qid}},"$sid,$pm,$pe,$pi,$s");
						$lastpm=$pm;$lastpe=$pe;
					}
					else {
						if (($lastpm==$pm) && ($lastpe==$pe)) {
							push(@{$sim{$qid}},"$sid,$pm,$pe,$pi,$s");
						}
						else {last LINE;}
					}
				}
			} else {write_log("$sid gindex2 not defined; lineid: $_\n");}
		}
	}
	foreach my $qid (@{$gindex{$tb}}) {
		my ($sStart,$sEnd);
		if (defined $blastquery{$qid}) {
			($sStart,$sEnd)=split (";",$blastquery{$qid});
		} else {next;}
		my ($lastpm,$lastpe);
		my $hit_id=0;
		LINE:foreach ($sStart..$sEnd) {
			my ($s,$sid,$pm,$pe,$pi)=(&getline_from_bpofile($_))[0,3,5,6,7];
			if (defined $gindex2{$sid}) {
				if ($gindex2{$sid} eq $ta) {
					$hit_id++;
					if ($hit_id==1) {
						push(@{$sim{$qid}},"$sid,$pm,$pe,$pi,$s");
						$lastpm=$pm;$lastpe=$pe;
					}
					else {
						if (($lastpm==$pm) && ($lastpe==$pe)) {
							push(@{$sim{$qid}},"$sid,$pm,$pe,$pi,$s");
						}
						else {last LINE;}
					}
				}
			} else {write_log("$sid gindex2 not defined; lineid: $_\n");}
		}
	}

	foreach my $q (keys %sim) {
		foreach (@{$sim{$q}}) {
			my @bla=split (',',$_);
			next if($bla[1].'e'.$bla[2] > $pv_cutoff || $bla[3]< $pi_cutoff);
			if($pmatch_cutoff) {
				next if(&simspan($bla[4]) < $pmatch_cutoff);
			}
			push(@{$best{$q}}, $bla[0]);
			$pvalue{$q.' '.$bla[0]} = $bla[1].'e'.$bla[2];

		}
	}
	my @b = keys %best;
	write_log(scalar(@b)." sequences have best hits from the other species\n");
	return &matrix(\%best, \%pvalue);
} ## makeOrtholog




# This subroutine is used to choose two-way hits among one-way hits (best
# hits between two species or better hits within one species), 
# calculate the weight between two nodes (minus logrithm of the p-value, 
# or $MAX_WEIGHT_DEFAULT for p-value 0 ), and calculate average
# weight among all inparalogs within one species or all orthologs between
# two species. (Weighting process takes place in the main script)
# Two Arguments:
# 1. Reference Variable: reference to a hash which stores all the possible
#    gene pairs (one-way best hit, or better hit).
# 2. Reference Variable: reference to a hash which stores the pvalue for
#    the gene pairs.
# Last modified: 10/02/06
sub matrix {
	my %best      = %{$_[0]};
	my %pvalue    = %{$_[1]};
	my (%edge, %weight);
	my $count=0;
	my $sumw=0;

	foreach my $query (sort keys %best) {
		foreach my $subject (@{$best{$query}}) {
			next if($weight{$query.' '.$subject});
			my $flag = 0;
			foreach my $q (@{$best{$subject}}) {
				if($q eq $query) { $flag = 1; }
			}
			if($flag == 1) {
				push (@{$edge{$query}}, $subject);
				push (@{$edge{$subject}}, $query);
				#use -logP as weights and treat P=0 as -logP=$maximum_weight (DEFAULT=300)
				my ($pv1, $pv2);
				if($pvalue{$query.' '.$subject} == 0) {
					$pv1 = $maximum_weight;
				}else { 
					$pv1 = -log($pvalue{$query.' '.$subject})/log(10);
				}	    
				if($pvalue{$subject.' '.$query} == 0) {
					$pv2 = $maximum_weight;
				}else {
					$pv2 = -log($pvalue{$subject.' '.$query})/log(10);
				}

				write_partbbh($CUR_PARTBBH,"$query\t$subject\t".$pvalue{$query.' '.$subject}."	".$pvalue{$subject.' '.$query}."\n");

				my $w = ($pv1+$pv2)/2;
				$sumw += $w;
				$count++;
				# use averaged score as edge weight
				$weight{$query.' '.$subject} = sprintf("%.3f", $w);
				$weight{$subject.' '.$query} = sprintf("%.3f", $w);
			}
		}
	}
	my $avgw = 'N/A';
	if ($count) {
		$avgw = $sumw/$count;
	}
	my $no_tmp = scalar(keys %weight)/2;
	write_log("$no_tmp sequence pairs were identified as Reciprocal Better/Best Hit\n");
	return (\%edge, \%weight, $avgw, $sumw, $count);
} ## matrix


sub read_bbhmatrix { ## dgg add
	my ($bbhpart)= @_;
  my (%best,%pvalue);
	my (%edge, %weight);
	my $count=0;
	my $sumw=0;

  ## $CUR_PARTBBH= $bbhpart; #??
  open_partbbh($bbhpart) or die "ERR: cant read $bbhpart bbh";
  my $bbhin= $bbhpart_fh{$bbhpart};
  while (<$bbhin>) {
    next unless(/^\w/);
    my($query,$subject,$pv1,$pv2)=split;
    
    ## WRONG >>> NOTE: pv already is  : -log($pvalue{$query.' '.$subject})/log(10); NONONO !! THIS IS BUG
    ## see above: write_partbbh($CUR_PARTBBH,"$query\t$subject\t".$pvalue{$query.' '.$subject}."	".$pvalue{$subject.' '.$query}."\n");
    
    $pvalue{$query.' '.$subject}= $pv1;
    $pvalue{$subject.' '.$query}= $pv2;
    push @{$best{$query}}, $subject;
    push @{$best{$subject}}, $query;
  }
  close_partbbh($bbhpart);
  
  #..... and redo matrix() calc here
	foreach my $query (sort keys %best) {
		foreach my $subject (@{$best{$query}}) {
			next if($weight{$query.' '.$subject});
			my $flag = 0;
			foreach my $q (@{$best{$subject}}) { if($q eq $query) { $flag = 1; last; } }
			if($flag == 1) {
				push (@{$edge{$query}}, $subject);
				push (@{$edge{$subject}}, $query);
				#use -logP as weights and treat P=0 as -logP=$maximum_weight (DEFAULT=300)
				my ($pv1, $pv2);
				$pv1= $pvalue{$query.' '.$subject};
				$pv2= $pvalue{$subject.' '.$query};
    ## WRONG >>> NOTE: pv already is  : -log($pvalue{$query.' '.$subject})/log(10); NONONO !! THIS IS BUG
				if($pvalue{$query.' '.$subject} == 0) {
					$pv1 = $maximum_weight;
				}else { 
					$pv1 = -log($pvalue{$query.' '.$subject})/log(10);
				}	    
				if($pvalue{$subject.' '.$query} == 0) {
					$pv2 = $maximum_weight;
				}else {
					$pv2 = -log($pvalue{$subject.' '.$query})/log(10);
				}
				#NOT# write_partbbh($CUR_PARTBBH,"$query\t$subject\t".$pvalue{$query.' '.$subject}."	".$pvalue{$subject.' '.$query}."\n");
				my $w = ($pv1+$pv2)/2;
				$sumw += $w;
				$count++;
				my $wtf= sprintf("%.3f", $w); # use averaged score as edge weight
				$weight{$query.' '.$subject} = $wtf; $weight{$subject.' '.$query} = $wtf;
			}
		}
	}
	my $avgw = 'N/A';
	if ($count) {
		$avgw = $sumw/$count;
	}
	my $no_tmp = scalar(keys %weight)/2;
	write_log("read_bbhmatrix $bbhpart: $no_tmp sequence pairs are Reciprocal Better/Best Hit\n");
	return (\%edge, \%weight, $avgw, $sumw, $count);
} ## read_bbhmatrix


sub read_connpart { ## dgg add
	my ($itaxa)= @_;
  my ($avgw, $sumw, $count,$key,$taxai,$taxaj);
  my $itaxname= $taxa[$itaxa];
  my $partconn="connect.$itaxname.con"; # put all jtax in itax.con
  my $nread=0;
  
  ## maybe par_conn calc bug is from %edge, %weight reading .. need local/my full nested arrays?
  ## 		foreach my $pi (keys %$edge_ref) { @{$e{$pi}}=@{$edge_ref->{$pi}}; }  #make a copy of current edge data structure into %e
  ## replace edge w/ edger={} foreach TAX ?
  my( $edger, $weightr); # old: my(%edge, %weight);
  
  write_log("read_connpart: $partconn\n");  
  open_partbbh($partconn) or return 0; # die "ERR: cant read $partconn conn"; # or return -1 ?
  
  my $bbhin= $bbhpart_fh{$partconn};
  while (<$bbhin>) {
    if(/^TAX:/) {
      my %edge_local=(); $edger=\%edge_local; 
      my %weight_local=(); $weightr=\%weight_local; # BUG was here global package %weight zeroed .. bad news
      ($avgw, $sumw, $count,$key,$taxai,$taxaj)= (0) x 10;
      ($key, $taxai, $taxaj)=split;    
    } elsif(/^AV:/) {
      ($key, $avgw, $sumw, $count)=split;    
    } elsif(/^EG:/) {
      my($key1,$query,@subjects)=split;
      push (@{$edger->{$query}}, @subjects);
      $ortho{$query}= 1; # see write_
    } elsif(/^WT:/) {
      my($key1,$query,$subject,$wtf)=split;
      $weightr->{$query.' '.$subject} = $wtf;
      #?? does it need recip: wt{sb x qu}= wtf ?? shoudnt, this is all written from weight_ref hash
    } elsif(/^END/) {
      @{$connect{$taxai.' '.$taxaj}}  = ( $edger, $weightr, $avgw, $sumw, $count);
      ## 		foreach my $n (keys %e) { $ortho{$n} = 1; .. } ## ** %ortho is global also
      $nread++; write_log("read_connpart$nread $taxai-$taxaj: $avgw, $sumw, $count\n");  
      ($avgw, $sumw, $count,$key,$taxai,$taxaj)= (0) x 10;
    }
  }
  close_partbbh($partconn);
  return $nread;
}

sub write_connpart { ## dgg add
	my ($itaxa)= @_;
  my $itaxname= $taxa[$itaxa];
  my $partconn="connect.$itaxname.con"; # put all jtax in itax.con

  write_log("write_connpart: $partconn\n");  
  open_partbbh($partconn,"write") or die "ERR: cant write $partconn";
  my $bbhin= $bbhpart_fh{$partconn};

	for(my $j=$itaxa+1;$j<scalar(@taxa);$j++) { #?? problem here for j>last taxa? lastspp.con is empty
	  my $jtaxname= $taxa[$j];
	  my ($edger,$weightr,$avew,$sumw,$count)= @{$connect{$itaxname.' '.$jtaxname}};

    print $bbhin join("\t","TAX:",$itaxname,$jtaxname)."\n";
    print $bbhin join("\t","AV:",$avew,$sumw,$count)."\n";
    foreach my $qs (sort keys %$edger) { ## @{$edge{$query}}, @subjects
      my @ss= @{$edger->{$qs}};
      print $bbhin join("\t","EG:",$qs,@ss)."\n";
    }
    foreach my $qs (sort keys %$weightr) { ## ($query,$subject,$wtf)
      my $wtf= $weightr->{$qs};
      print $bbhin join("\t","WT:",$qs,$wtf)."\n";
    }
    print $bbhin "ENDTAX:\n";
  }
  close_partbbh($partconn);
}





#d # This subroutine is used by the subroutine makeInparalog,
#d # to solve the pv_tie problem. It rearranges the pv-tied blast hits
#d # so that the hits from a specific taxon are moved higher than hits from
#d # other species.
#d # Three Arguments:
#d # 1. Number Variable: starting line id (or similarity id) of bpo file (blast
#d #    parse out file)
#d # 2. Number Variable: ending line id (or similarity id) of bpo file (blast
#d #    parse out file)
#d # 3. String Variable: taxon
#d # Last modified: 07/20/04
#d sub pvtie_sort {
#d 	my ($sStart,$sEnd,$taxon)=@_;
#d 	my (@sorted_simid,@tmp);
#d 	my ($lastpm,$lastpe)=(getline_from_bpofile($sStart))[5,6];
#d 	foreach ($sStart..$sEnd) {
#d 	  local $^W = 0; #dgg no warnings missing $gindex2{sid}
#d 		my ($s,$sid,$pm,$pe,$pi)=(&getline_from_bpofile($_))[0,3,5,6,7];
#d 		if (($lastpm==$pm) && ($lastpe==$pe)) {
#d 			if ($gindex2{$sid} eq $taxon) {
#d 				push (@sorted_simid,$s);
#d 			}
#d 			else {
#d 				push (@tmp,$s);
#d 			}
#d 		}
#d 		else {
#d 			if (scalar(@tmp)>0) {
#d 				push (@sorted_simid,@tmp);
#d 				@tmp=();
#d 			} 
#d 			if ($gindex2{$sid} eq $taxon) {
#d 				push (@sorted_simid,$s);
#d 			}
#d 			else {
#d 				push (@tmp,$s);
#d 			}
#d 		}
#d 		$lastpm=$pm;$lastpe=$pe;
#d 	}
#d 	if (scalar(@tmp)>0) {push (@sorted_simid,@tmp);}
#d 	return @sorted_simid;
#d } ## pvtie_sort





# This subroutine, together with matchlen, are used to calculate
# how much of the query sequences match each other.
# One Argument:
# 1. Number Variable: line id (or similarity id) of bpo file (blast
# parse out file)
# Last modified: 07/21/04
sub simspan {
	my $s = $_[0];
	my (%sub_start, %sub_length, %query_start, %query_length);
	my @bporow= &getline_from_bpofile($s); #dgg; improve efficiency here
	my @hsp=split ('\.',(@bporow)[8]);
	foreach (@hsp) {
		if (/(\d+)\:(\d+)\-(\d+)\:(\d+)\-(\d+)/) {
			$sub_start{$1}=$4; 
			$sub_length{$1}=$5-$4+1;
			$query_start{$1}=$2;
			$query_length{$1}=$3-$2+1;
		}
	}
	my $match_lengths = &matchlen(\%sub_start,\%sub_length);
	my $match_lengthq = &matchlen(\%query_start,\%query_length);			
	my ($lengthq,$lengths)=(@bporow)[2,4];   # June 3
	if($lengths >= $lengthq) {
		return 100*$match_lengthq/$lengthq;
	}else{
		return 100*$match_lengths/$lengths;
	}
} ##simspan





# This subroutine, together with simspan, are used to calculate
# how much of the query sequences match each other.
# Two Arguments:
# 1. Reference Variable: reference to an hash which stores the starting
#    position of each HSP.
# 2. Reference Variable: reference to an hash which stores the length
#    of each HSP.
# Last modified: 07/19/04
sub matchlen {

	my %start        = %{$_[0]}; 
	my %length       = %{$_[1]};

	my @starts = sort{$start{$a}<=>$start{$b}} (keys %start);
	return $length{$starts[0]} if(scalar(@starts)==1);
	my $i=1; 
	my  $match_length = $length{$starts[0]}; 
	my $pos = $length{$starts[0]} + $start{$starts[0]} ;
	while($i<scalar(@starts)) {

	if($length{$starts[$i]} + $start{$starts[$i]} <= $pos) {
		$i++;
		next;
	}
	if($start{$starts[$i]}> $pos) {
		$match_length += $length{$starts[$i]};
		$pos = $start{$starts[$i]} + $length{$starts[$i]};
	}else {
		$match_length += $length{$starts[$i]} - ($pos - $start{$starts[$i]});
		$pos = $start{$starts[$i]} + $length{$starts[$i]};
	}
	$i++;
	}

	return $match_length;
} ## matchlen



# Last modified: 07/22/04
sub printHelp {
	my (@foo) = <DATA>;
	print STDERR "OrthoMCL V$VERSION\n";
	print STDERR @foo;
	exit 1;
}



######################################USAGE OF ORTHOMCL.PL###################################
__DATA__

Copyright (C) 2004-2006 by University of Pennsylvania,
Philadelphia, PA USA. All rights reserved.

Before orthomcl.pl can be used, some variables 
(including directory variables or parameter variables)
in orthomcl_module.pm need to be set, as described in
README.

Usage: orthomcl.pl --mode 1,2,3,4 or 5 <tagged arguments>

Modes:
~~~~~~

 1: OrthoMCL analysis from FASTA files
% orthomcl.pl --mode 1 --fa_files Ath.fa,Hsa.fa,Sce.fa

 2: OrthoMCL analysis based on former OrthoMCL run. No BLAST or BLAST
 parsing performed.
% orthomcl.pl --mode 2 --former_run_dir Sep_8 --inflation 1.4

 3: OrthoMCL analysis from user-provided BLAST result. No BLAST 
 performed.
% orthomcl.pl --mode 3 --blast_file AtCeHs_blast.out --gg_file 
AtCeHs.gg

 4: OrthoMCL analysis based on user-provided BPO (BLAST PARSE OUT) 
 file and GG (Genome-Gene Index) file
% orthomcl.pl --mode 4 --bpo_file AtCeHs.bpo --gg_file AtCeHs.gg

 5: OrthoMCL analysis based on matrix of former OrthoMCL run, but with
 LESS genomes
% orthomcl.pl --mode 5 --former_run_dir Sep_8 --taxa_file AtCeHs.gg
--inflation=1.1

Arguments:
~~~~~~~~~~

 fa_files=<String>       Protein FASTA file names, with each file 
                         containing protein sequences from one species,
                         separated by comma(e.g. "Eco.fa,Sce.fa,Afu.fa")
 pv_cutoff=<Float>       P-Value or E-Value Cutoff in BLAST search and/or
                         ortholog clustering, 1e-5 (DEFAULT).
 pi_cutoff=<Int>         Percent Identity Cutoff <0-100> in ortholog 
                         clustering, 0 (DEFAULT).
 pmatch_cutoff=<Int>     Percent Match Cutoff <0-100> in ortholog
                         clustering, 0 (DEFAULT).
 inflation=<Float>       Markov Inflation Index, used in MCL algorithm,
                         1.5 (DEFAULT). Increasing this index increases
                         cluster tightness, and the number of clusters.
 former_run_dir=<String> Former run directory, required in Mode 2, e.g. 
                         "July_21". Then the blast result file and bpo 
                         file in former run directory will be used 
                         instead of running from the very beginning.
 blast_file=<String>     Blast out file provided by user, required in
                         Mode 3. It will be parsed into BPO file, and
                         further used for ortholog clustering.
 bpo_file=<String>       BPO (Blast Parse Out) file provided by user,
                         required in Mode 4. Please refer to README 
						 about its format.
 gg_file=<String>        GG (Genome Gene mapping) file provided by user,
                         required in Mode 3 & 4. Please refer to 
                         README about its format.
 taxa_file=<String>      TAXA file provided by user, required in Mode 5. 
                         Please refer to README about its
                         format.
 
