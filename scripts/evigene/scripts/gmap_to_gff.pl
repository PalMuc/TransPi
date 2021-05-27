#!/usr/bin/env perl

use strict;
use warnings;

=item about

=item notes

 dgg mods for multiple-path EST gmap
 'Alignment for path 1' preceding align line
 add parse of match span...
 dgg add dobest to btab output ... 2010sep
 .. updated for gmap.out 2011.08 ; chimera default; GFF: add intron output; add CDS output

=item CLEAN ME UP
	
	this script has gotten messy .. needs cleanup.
	add GetOptions instead of env params=xyz gmap2gff ..

=cut

my $VERSION='2014.06.07'; # Apis wacky 'Group9.123' scaf ref names '.' threw out ref.begin.end parser
# '2013.07.18'; # .. many tiny changes since "2011.10.27"; # "2011.08.24";

my $asGFF=$ENV{gff}||0; #dgg
$asGFF=1 if($0 =~ /gff/); # gmap_to_gff or gmap_to_btab ..
my $dobest=1; my $bestoff=1; # was 5 = 5% lower match, too low;
if(defined $ENV{best}) { $bestoff=$ENV{best}; $dobest=1; } # allow bestoff=0 == same quality
elsif(defined $ENV{nobest}) { $dobest=0; }
my $NOPATHGFF= $ENV{nopath}||0;

## FIXME: need options, turn-off bad skip-low-pctalign in asgff()
my $LOWIDENT=$ENV{low} || 90;
my $SKIPLOW=  (defined $ENV{skiplow}) ? $ENV{skiplow} : 1;

my $MAX_SPAN_ERR=100; #? error if exon span genomic - transcript > this, too big an error
my $MAX_SPAN_ERR_HIDE= 9999; # inner span error too big and cannot fix, hide this as #err

my $source  = $ENV{src}  || "gmap"; # PASA.btab needs this
my $mrnaType= $ENV{mrna} || "mRNA"; # "match";
my $exonType= $ENV{exon} || "exon"; # "HSP";
my $dointron= (defined $ENV{intron}) ? $ENV{intron} : 2; # gff only
my $docds   = (defined $ENV{cds}) ? $ENV{cds} : 1; # gff only
my $noerrspan= $ENV{noerrspan}||0;  # cancel badspan changes
my $GSTRANDED= $ENV{strand}||0;  # when "genome" is stranded, like mRNA ref

# special: print ONLY chimeras: 2+ paths, high ident, partial cover
# rename: dochim => ONLYchimera
my $dochim=$ENV{chim} || $ENV{chimera} ||0;
if($dochim) { $dobest=0; $asGFF=1; } #??

my $pathnum= 0;
my $transcript_acc = "";
my $chain_counter = 0;
my $segment_counter = 0;
my ($ipath,$ischim,$npath,%path); #dgg
my %flipor=( "+" => "-", "-" => "+", "."=>"." );
my ($lref,$lgend5,$lgend3,$ltend5,$ltend3);
my ($badspan, $err_setorient); # see below note, buggy spans 2011.08
my $cdstrlen = 0; # see below note, buggy spans 2011.08
my @aligno; # hold align output till check for error spans, correct mRNA span.

my @keepattr=(); # from >header attributes
my $KEEPATTR=$ENV{keepat}||"";
if($KEEPATTR) {
 my @at= split /[\s\|,;]+/, $KEEPATTR; my $ap="";
 if(@at>1) { $ap= join '|', @at; $ap= '(?:'. $ap .')=[^;=\s]+';  } 
 elsif(@at==1) { $ap= $at[0] . '=[^;=\s]+'; }
 $KEEPATTR=$ap; # warn "#DBG: keepat='",$KEEPATTR,"'\n";
}

#.....................
sub MAIN {}

print "##gff-version 3\n# gmap_to_gff $VERSION\n# authors: Brian Haas (PASA) + Don Gilbert (updates)\n" 
  if($asGFF);


while (<STDIN>) {
	## KEEPATTR OPTION: parse more of >faheader for attribs to keep in mRNA attrrib
	## here want: aalen, offs, oid
	## FunhetEGm031420t1 type=mRNA; aalen=414,100%,partial; clen=1242; offs=1-1242; oid=kfish2qf7soapk21loc1090t1; 

    if (/^>(\S+)/) {
        putalign(@aligno) if @aligno;  @aligno=(); 
        $transcript_acc = $1;
				if($KEEPATTR) { @keepattr= m/($KEEPATTR)/g; }
	 				##warn "#DBG: $transcript_acc atr: ",join(", ",@keepattr),"\n";
        $chain_counter++;
        $segment_counter = 0; #reinit
        %path=(); 
        $err_setorient= $badspan= $ischim=$npath=$ipath=0;
    }
    elsif(/^Paths .(\d+)/){ $npath=$1;
	      if(/chimera with (.*)$/i) { $ischim=$1; # \(
	        $ischim =~ s/\s*\([^\)]+\)//g;  # new # 
	        $ischim =~ s/ .donor_prob.*$//; # old
	        }
        #OPTION: npath == 0, print comment or fake mRNA line with NOPATH ref..
        if($npath == 0 and $asGFF) {
	  			print "#NOPATH: $transcript_acc\n";
          if($NOPATHGFF) {  # add fake mRNA for further processing
						# ID=Nasvi2EG000074t1;Target=Nasvi2EG000074t1 1 7841;aalen=1627;cov=100.0;match=7841;nexon=20;pid=100.0;qlen=7841
            my @xout;
	    			my $ar="Target=$transcript_acc 1 69;cov=0;match=0;pid=1;qlen=99";
            $ar.=";path=0/0"; # if ($showpath);
            push @xout, join("\t","NOPATH",$source,$mrnaType,1,69,0,".",".","ID=$transcript_acc;$ar")."\n";
    	    	$ar="Target=$transcript_acc 1 69";
            push @xout, join("\t","NOPATH",$source,$exonType,1,69,0,".",".","Parent=$transcript_acc;$ar")."\n";
	    			putalign(@xout);
	  			}
	 			}
    	}
    elsif(/^\s*Alignment for path (\d+)/) { 
      putalign(@aligno) if @aligno;  @aligno=();
      $pathnum=$1; $segment_counter=0; 
      $err_setorient= $badspan=0; 
      ($lref,$lgend5,$lgend3,$ltend5,$ltend3)=(0) x 5;
      } # dgg
    
    elsif (/\s*([\+\-])([^\:]+):(\d+)-(\d+)\s+\((\d+)-(\d+)\)\s+(\d[^\%]+\%)/) {
			# FIXME^ complex regex may break, need alternate reporting for  Alignment line not matching
			
## parsing this Alignment line..
#     -scaffold_3:29457826-29457762  (766-830)   100% ->   ...471...  0.995, 0.994
#     -scaffold_3:29457290-29457050  (831-1071)   99%
#     +scaffold_2:38213061-38214311  (1347-2597)   99% <-   ...1808...  0.995, 0.993
#     +scaffold_2:38216120-38216516  (2598-2994)   100% (-   ...328...  0.838, 0.833


        my $orient = $1;
        my $genome_acc = $2;
        my $genome_end5 = $3;
        my $genome_end3 = $4;
        my $transcript_end5 = $5;
        my $transcript_end3 = $6;
        my $percent_identity = $7;
        
        my @x; # = ("") x 15;
        $#x = 14;  foreach my $ele (@x) {  $ele = ""; }  #prealloc # init list

# parse out introns *
#  variations in -> <-  are  -) -] ==
#     -scaffold_9:35497474-35497115  (1-353)   95% (-   ...96...  0.993, 0.721
#     +scaffold_9:35902128-35902279  (1999-2150)   100% [-   ...542...  0.003, 0.005
#     -scaffold_8:5039941-5039925  (502-518)   94% ==   ...371...  0.877, 0.018
#     +scaffold_2:41253101-41253133  (404-437)   97% ==   ...30...   ***query_skip:5***  0.000, 0.188

        my @intron=(0,0,0);
        if($dointron and m/$percent_identity\s+(\S+)\s+\.+(\d+)\.+/) {  # ...\s+([\d\.]+), ([\d\.]+)
          my($iqual, $ispan)= ($1,$2);
          if($iqual eq "==") { }
          else {
            my $ip=0;
            if( my($ip1,$ip2)= m/([\d\.]+), ([\d\.]+)$/) {  $ip= int(50*($ip1+$ip2)); }
            if($genome_end5 > $genome_end3) { ## NOT $orient eq "-"  # end3 < end5 here
              @intron= ($genome_end3 - $ispan, $genome_end3 - 1, $ip); # abs start,stop
            } else {
              @intron= ($genome_end3 + 1, $genome_end3 + $ispan, $ip);
            }
          }
        }
        
        $percent_identity =~ s/\%//;
        $segment_counter++;

        # need to collect full transcript and comment out if badspan?
        # or correct exon genome span to tr span?
        # ** Fixme: allow *shorter* genome span? = partial match
        my $gspan= abs($genome_end5 - $genome_end3);
        my $tspan= abs($transcript_end5 - $transcript_end3);
        if(abs($gspan - $tspan) > $MAX_SPAN_ERR) { 
	        my $badspan1= ($gspan - $tspan);
          $badspan +=  abs($badspan1);
          #** these all are for chimera? some are, not all
          $x[11] = "ERROR.span:genome_span:$gspan,tr_span:$tspan,$genome_acc:$genome_end5-$genome_end3";
          # fix if know last/next span to adjust proper end
          # ** also need to fix mRNA span (for gff) but too late here, print that at first align line
	        if($noerrspan) {  # cancel changes, but report in attr
	          $x[11].=",notfixed=$badspan1";
            $badspan=0;	
          } elsif($lgend5 > 0 and $lgend3 > 0) {  # OR segment_counter > 1
            my $d3= abs($genome_end3 - $lgend3) ;
            my $d5= abs($genome_end5 - $lgend5) ;
            if($d5 < $d3) { # adjust end3
              if($genome_end3 > $genome_end5) { $genome_end3= $genome_end5 + $tspan; }
              elsif($genome_end3 < $genome_end5) { $genome_end3= $genome_end5 - $tspan; }
            } else { # adjust end5
              if($genome_end3 > $genome_end5) { $genome_end5= $genome_end3 - $tspan; }
              elsif($genome_end3 < $genome_end5) { $genome_end5= $genome_end3 + $tspan; }
            }
            $x[11].=",spanfix=$tspan";
            
          } else { # 1st of align spans ; wait for next to adjust? skip?
            if($genome_end3 > $genome_end5) { $genome_end5= $genome_end3 - $tspan; }
            elsif($genome_end3 < $genome_end5) { $genome_end5= $genome_end3 + $tspan; }
            $x[11].=",spanguess=$tspan";
          }
        } #... skip?
        
        $x[0] = $genome_acc;
        $x[3] = $source;
        $x[5] = $transcript_acc;
        $x[6] = $genome_end5;
        $x[7] = $genome_end3;
        $x[8] = $transcript_end5;
        $x[9] = $transcript_end3;
        $x[10] = $percent_identity;
#dgg:   $x[11] = what? $x[12] == segmentscore (blat); or change chain_counter?
        $x[12] = $pathnum;
        $x[13] = $chain_counter;
        $x[14] = $segment_counter;
        
        ($lref,$lgend5,$lgend3,$ltend5,$ltend3)= 
          ($genome_acc, $genome_end5, $genome_end3,$transcript_end5, $transcript_end3);
        
        # my $btab_line = join ("\t", @x);
        # print "$btab_line\n" if($btab_line);
        my @exons = ($asGFF)? as_gff($orient, @x, @intron) : as_btab($orient,@x);  
        push @aligno, @exons if(@exons);
        
    } elsif($chain_counter and $segment_counter == 0) {
      if(/^\s+Path (\d+):/) { $ipath=$1; 
        my($targ)= m/ query (\S+)/; $targ =~ s/[\.\-]+/ /; $path{$ipath}{Target}=$transcript_acc." ".$targ;
	my($ref)= m/ =\> \S+ (\S+)/;  my $or="+";
        #o#my($r,$b,$e)= split/[:\.\-]+/, $ref; # Problem here Dang Apis: Group3.45 << dingbat '.' in scaf names
        my($r,$be)= split/:/,$ref;  my($b,$e)= split/\.\./,$be;
        # Apis:   Path 1: query 1..196 (196 bp) => genome Group6.16:957,359..957,164 (-196 bp)
        if($e) { $b=~s/,//g; $e=~s/,//g; ($b,$e,$or)=($e,$b,"-") if($b>$e); 
          $path{$ipath}{ref}=$r; $path{$ipath}{b}=$b; $path{$ipath}{e}=$e; $path{$ipath}{orient}=$or;
	}
      }
      
      elsif($ipath and /direction: indeterminate/ and not $GSTRANDED) { $path{$ipath}{orient}="."; }
        ##^^ 2011.08 bug this set when have exons, introns w/ strand !
      
      	## ** BUG in orient, cDNA direction: antisense  << this means flip orient, cDNA is comprev
      elsif($ipath and /direction: antisense/) { 
        $path{$ipath}{orient}=$flipor{ $path{$ipath}{orient} }; $path{$ipath}{sense}=-1; }
      elsif($ipath and /Number of exons: (\d+)/) { $path{$ipath}{nexon}=$1; }
      elsif($ipath and /Translation: (\d+)..(\d+) .(\d+) aa/) { 
        $path{$ipath}{cdsb}=$1; $path{$ipath}{cdse}=$2; $path{$ipath}{aalen}=$3; }
      elsif($ipath and /Coverage: (\S+)/) { $path{$ipath}{cov}=$1;
	      if( m/query length: (\d+)/){ $path{$ipath}{qlen}=$1;} }
      elsif($ipath and /Percent identity: (\S+)/) {  $path{$ipath}{pid}=$1; 
        if(m/(\d+) matches/){ $path{$ipath}{match}=$1;} }

      ## 2012jun: add count of 'Non-intron gaps:' to mRNA annot. also flag when cds-tr span != cds-genome span
      # Non-intron gaps: 0 openings, 0 bases in cdna; 1 openings, 1 bases in genome
      # Non-intron gaps: 3 openings, 12 bases in cdna; 2 openings, 9 bases in genome
      elsif($ipath and /Non-intron gaps:/) { 
        my($gapt)=m/(\d+) bases in cdna/; my($gapg)=m/(\d+) bases in genome/;
        $path{$ipath}{indels}="$gapt/$gapg" if($gapt>0 or $gapg>0);  
        }
        
    }
}

putalign(@aligno) if @aligno; @aligno=();

#-------------------------------------------------------------------------

sub putalign {
  my(@exons)= @_;
  if($asGFF) {
  
    ## calc if cds-transcript size <> cds-genomic size
    ## need to calc cdsglen from all cds exons in @aligno
    if($docds) {
      # my $cdstrlen=0; from below...
      my $cdsglen= 0;
      foreach my $xn (@exons) { 
        next unless($xn =~ /\tCDS\t/);
        my @xn= split"\t", $xn; 
        my $w=1+ $xn[4] - $xn[3];
        $cdsglen+= $w;
        }
      unless($cdsglen == $cdstrlen) { 

## FIXME: cdsindel can be way wrong; on chimera2 only? use $path{$ipath}{indels} as max value?
## my $pathindel= $path{$ipath}{indels}; my ($in1,$del1)= ($pathindel)? split "/", $pathindel : (0,0);
## my $indels= ($in1>$del1)? $in1 : $del1;
## $d=$indels if($d>$indels);

        my $d= $cdstrlen - $cdsglen; $exons[0]=~s,$,;cdsindel=$d,;  # :$cdstrlen/$cdsglen
        }
    }
    
    if(@keepattr) {
      my($m,$mi,$i)=(0) x 3;
      foreach my $xn (@exons) { if($xn =~ /\t$mrnaType\t/) { $m=$xn; $mi=$i; last; } $i++; }
      if($m) { my $kat=join(";",@keepattr); $m =~ s/$/;$kat/; $exons[$mi]= $m;  @keepattr=(); }
    } # 2013jul add

    if($badspan or $err_setorient) {
      my ($i, $mi, $m, $mb, $me, $nexon)= (0) x 10;
      my @ierr=();
      $mi= -1;
      
# ** some bad spans at inner exon; cannot be corrected this way;
# .. split to chimera ?? or drop? or what

      foreach my $xn (@exons) { 
      
        my @xn= split"\t", $xn; 
        if($err_setorient) { 
          $xn[6]= $err_setorient; 
          $xn= join "\t", @xn; 
          $exons[$i]= $xn;
          }
          
        if($xn =~ /\t$mrnaType\t/) {  
	  			$m=$xn; $mi=$i;  ## @mr= @xn;
        } elsif($xn =~ /\t$exonType\t/) { 
          my($r,$b,$e)= @xn[0,3,4]; ## (split"\t",$xn)[0,3,4];
          $mb= $b if( $mb == 0 or $b < $mb); 
          $me= $e if( $e > $me );
          $nexon= 1 + $i; push @ierr, 1+$i if( $xn =~ m/ERROR.span/); # check for inner exon err
          }
        $i++;
      }
      
    if($badspan and $m) { ## add mrna attr from exon ERROR:...
      my $innerbad=0;
      if(@ierr) { 
        my @inner= grep { $_ > 1 and $_ < $nexon } @ierr; 
        $innerbad=join",",@inner if(@inner); # cant handle by mrna span adjust
      }
      
      my @m= split "\t", $m; 
      my($mr,$mbo,$meo)= @m[0,3,4]; 

      if($innerbad or $mb != $mbo or $me != $meo) {
        $m[3]= $mb; $m[4]= $me; 
        $m[-1] =~ s/$/;errspaninner=$innerbad/ if($innerbad);
        $m[-1] =~ s/$/;errspan=$mr:$mbo-$meo/;
        $exons[$mi]= join "\t", @m;
        }
        
      if($innerbad and $badspan >= $MAX_SPAN_ERR_HIDE) { map{ s/^/#err./ } @exons; } #?? or not; check spanerr size; only if huge err
      }      
      
    }
    print @exons;
  } else {
    print @exons;
  }
}


sub pathscore {
 my($ipath)= @_;
 my $sc = 100;
 if($path{$ipath}) {
   # my %maink= map{ $_,1 } qw(ref b e orient);
   my($r,$b,$e,$or)= @{$path{$ipath}}{qw(ref b e orient)};
   my($ma,$ql)= ($path{$ipath}{match}||0, $path{$ipath}{qlen}||0);
   $sc= ($ql > 0 and $ma > 0) ? int(100 * $ma / $ql) : 1;
  }
 return $sc;
}


sub as_gff {
  return as_what("gff",@_);
}
sub as_btab {
  return as_what("btab",@_);
}

sub as_what {
 my ($what,$orientaln,@x)=@_;
 
 my $orient= $orientaln; # sense vs align orient
 my @xout=();
 my @g= @x[0,3,3,6,7,10,3,3,5];
 my($b,$e)=@g[3,4]; @g[3,4]=($e,$b) if($b>$e);
 my $ipath= $x[12];
 my $error= $x[11]; # flag in unused slot
 my $id=$x[5]; 
 if($ischim) { $id .= "_C".$ipath; }
 elsif($ipath>1) { $id .= "_G".$ipath; } # if ($ipath > 1 or $ischim);
 my ($isbad,$nextbad,$showpath, $nput)=(0) x 4;

 return @xout if($dochim and $npath < 2);

 if($path{$ipath}){
   # my %maink= map{ $_,1 } qw(ref b e orient);
   my %skipk= map{ $_,1 } qw(ref b e orient cdsb cdse);   ## 2012: added indels/Non-intron gaps as qual flag
   my($r,$b,$e,$or)= @{$path{$ipath}}{qw(ref b e orient)};
   my($ma,$ql)= ($path{$ipath}{match}||0, $path{$ipath}{qlen}||0);
   my $sc= ($ql > 0 and $ma > 0) ? int(100 * $ma / $ql) : 1;
   my $achim="";
   if($dochim or $ischim) {
     my($cov,$pid)= ($path{$ipath}{cov}||0, $path{$ipath}{pid}||0);     
       ## BAD here drops nearly full match cov > 95 as bad; do other end;  pid option??
       #BAD# $isbad = ($pid > 95 and $cov < 95) ? 0 : 1; # ?? change pct here...
		 ## fixme: $LOWIDENT and $SKIPLOW;  $SKIPLOW maybe NOT default ..
     $isbad = ($pid > 94 and $cov > 5) ? 0 : 1; # ?? change pct here...
     if($isbad and not $dochim) { $isbad=0 if($pid>=$LOWIDENT and $cov>20); }
     return @xout if($isbad and $SKIPLOW);

     my $i2= ($ipath==1) ? 2 : 1;
     my($xr,$xb,$xe,$xor)= @{$path{$i2}}{qw(ref b e orient)};
     $achim=";chim$i2=$xr:$xb-$xe:$xor";
   }
   elsif($dobest and $npath>1) {
     my $sc0 = ($ipath > 1) ? pathscore(1) : 0; 
     $isbad= $sc < $sc0 - $bestoff; #??
     return @xout if($isbad);
     my $sc1 = ($ipath < $npath) ? pathscore($ipath + 1) : $sc; 
     $nextbad = $sc1 < $sc - $bestoff;
   }
   my @k= grep { !$skipk{$_} } sort keys %{$path{$ipath}};
   my $ar= join ";", map{ $_."=".$path{$ipath}{$_} } @k;
   $showpath= ($npath>1 and ($ischim or $ipath>1 or !$nextbad));
   $orient= $or; #?? always reset, need path{}{sense}

   if($what =~ /gff/) {
     $ar.=";path=$ipath/$npath" if ($showpath);
     $ar.=";chimera=$ischim" if($ischim);
     $ar.= $achim if($achim);
     do{ 
      # print
      push @xout, join("\t",$r,$source,$mrnaType,$b,$e,$sc,$or,".","ID=$id;$ar")."\n";  
      $nput++; 
      } unless($path{$ipath}{did});
   }
   $path{$ipath}{did}++; # delete $path{$ipath};
 }
 
  if($what =~ /gff/) {
    $g[2]=$exonType; $g[6]= $orient; $g[7]=".";
    $g[8]="Parent=$id;Target=$x[5] $x[8] $x[9];ix=$x[14]";
    $g[8] .= ";error=$error" if $error;
    # print 
    push @xout, join("\t",@g)."\n"; $nput++; # return join ("\t", @g);

    if($docds) {
      # should add phase to cds, calc over exons..
      my($cdsb,$cdse)= ($path{$ipath}{cdsb}||0, $path{$ipath}{cdse}||0);
      ($cdsb,$cdse)= ($cdse,$cdsb) if($cdsb > $cdse);

 ## 201206: offby1,2 ? at CDS end for +/- strand, +/-sense ; cdsb,cdse from gmap.out are ok; tb,te should be ok
 ## Problem is not this CDS calc, but indels in genome vs transcript:
 ## ** count cds-trspan and cds-genomicspan and flag difference
 ## .. Non-intron gaps: 0 openings, 0 bases in cdna; 3 openings, 3 bases in genome >> inner-stop
 ## .. Non-intron gaps: 1 openings, 1 bases in cdna; 3 openings, 4 bases in genome >> inner-stop
 
 ## wrong cds offs when -sense rev ... exon order is backass
 # (- strand) + antisense == forward, but exons, cds rev on input
 ## fix: use orientaln instead of orient (-sense)
 ## fix111111: offby -1 CDS-end base; for +or, 1 CDS 
 ## Translation: 156..857 in exon 1-856 << bug here 857>856 .. this may be gmap bug (fixed?)
 ##   got 860069  860769 = 701 bp, should be 702 (1+857-156)
 
      $cdstrlen= 1+$cdse-$cdsb; # global for putalign!
      
      my($tb,$te)= @x[8,9]; # transcript start,end ; always forward?
      ($tb,$te)= ($te,$tb) if($tb > $te);
      if($tb < $cdse and $te > $cdsb) { 
        my($xb,$xe)= @g[3,4]; # genomic start,end
        my @cds= @g;  $cds[2]= "CDS"; $cds[8]= "Parent=$id"; 
        unless($tb >= $cdsb and $te <= $cdse) {
          # need rev 
          if($tb < $cdsb) { my $d= $cdsb - $tb; 
            if($orientaln eq "-") { $xe -= $d; } else { $xb += $d;} }
          if($te > $cdse) { my $d= $te - $cdse; 
            if($orientaln eq "-") { $xb += $d; } else { $xe -= $d;} }
          @cds[3,4]= ($xb,$xe);
        }
       # my $cdsglen= 1+$xe-$xb; # not right, need to remove introns.
       ## need to calc cdsglen from all cds exons in @aligno
       # unless($cdsglen == $cdstrlen) { my $d= $cdstrlen - $cdsglen; $xout[0]=~s,$,;cdsindel=$d,; } # :$cdstrlen/$cdsglen
       push @xout, join("\t",@cds)."\n"; $nput++;
      }
    }
    
    if($dointron and $x[16]>0) { # @intron addon to  @x
      # my @intron= @x[15,16]; # only gff
      my @in= @g;  $in[2]= "intron"; $in[8]= "Parent=$id";
      @in[3,4,5]= @x[15,16,17];
      if($orient eq ".") { # BAD strand error, fix
        my($gb,$ge)= @x[6,7];
        my $or= ($gb > $ge) ? "-" : "+"; # ($genome_end5 > $genome_end3)
        $in[6]= $or;  $err_setorient = $or;
      }
      # print 
      # * want intron but not always as gff row: add #i. comment if $dointron == 2
      my $ci= ($dointron == 2 or $dointron == -1)? "#i." : "";
      push @xout, $ci . join("\t",@in)."\n"; $nput++;
    }
    
  } else {  # btab format
    $x[5]= $id; # changed 
    $x[11]= ""; # drop err if there
    #print 
    push @xout, join("\t",@x)."\n"; $nput++; # return join("\t",@x);
  }
  return @xout; # return $nput;
}


__END__

=item 2014, gmap 2014 ; apis strand wrong, got + need -

# ** Fixme.revmrna+Group.num = 2 bugs in apis gmap out: 
#   1. many mrna wrong strand, -Group.scaf in gmap.out, +Group in gff
#   2. wacky bee scaffolds: Group 12.34 < .dot parsed bad, convert to _ on gmap.out input
#1 is likely result of #2 problem, see above .. orient came from begin <> end, begin wrong == orient wrong
  
gzgrep '=Apimel3aEVm000293t1;' evg3hbee_pub-amel4s.gmap.gff.gz 
Group9  << chopped '.12' bad, fix above
                evg3hbee_pub    mRNA    12      1373114 99      +       .       ID=Apimel3aEVm000293t1;Target=Apimel3aEVm000293t1 1 6252;aalen=1292;cov=100.0;indels=4/0;match=6236;nexon=29;pid=99.7;qlen=6252;cdsindel=3;aalen=2061,98%,complete;offs=57-6242;oid=hbee1uinorm4trinloc60214c0t398
Group9.12       evg3hbee_pub    exon    1372625 1373114 99      +       .       Parent=Apimel3aEVm000293t1;Target=Apimel3aEVm000293t1 1 492;ix=1
Group9.12       evg3hbee_pub    exon    1372415 1372533 99      +       .       Parent=Apimel3aEVm000293t1;Target=Apimel3aEVm000293t1 493 611;ix=2
Group9.12       evg3hbee_pub    exon    1371946 1372153 100     +       .       Parent=Apimel3aEVm000293t1;Target=Apimel3aEVm000293t1 612 819;ix=3
Group9.12       evg3hbee_pub    exon    1371607 1371851 99      +       .       Parent=Apimel3aEVm000293t1;Target=Apimel3aEVm000293t1 820 1064;ix=4
Group9.12       evg3hbee_pub    exon    1370905 1371463 100     +       .       Parent=Apimel3aEVm000293t1;Target=Apimel3aEVm000293t1 1065 1623;ix=5
Group9.12       evg3hbee_pub    exon    1369482 1369835 99      +       .       Parent=Apimel3aEVm000293t1;Target=Apimel3aEVm000293t1 1624 1977;ix=6
 .. etc.
 
>Apimel3aEVm000293t1 type=mRNA; Name=LOW QUALITY PROTEIN: protein still life, isoform; Dbxref=na,apimel:XP_006567523; aalen=2061,98%,complete; clen=6252; offs=57-6242; oid=hbee1uinorm4trinloc60214c0t398; organism=Daphnia_magna;
Paths (1):
  Path 1: query 1..6252 (6252 bp) => genome Group9.12:1,373,114..1,350,006 (-23109 bp)
    cDNA direction: sense
    Genomic pos: Amel_4.5_scaffolds:199,250,653..199,227,545 (- strand)
    Accessions: Group9.12:1,350,006..1,373,114 (out of 1871057 bp)
    Number of exons: 29
    Coverage: 100.0 (query length: 6252 bp)
    Trimmed coverage: 100.0 (trimmed length: 6252 bp, trimmed region: 1..6252)
    Percent identity: 99.7 (6236 matches, 12 mismatches, 4 indels, 0 unknowns)
    Non-intron gaps: 3 openings, 4 bases in cdna; 0 openings, 0 bases in genome
    Translation: 255..4134 (1292 aa)
    Amino acid changes: 

Alignments:
  Alignment for path 1:

    -Group9.12:1373114-1372625  (1-492)   99% ->   ...91...  0.997, 0.996
    -Group9.12:1372533-1372415  (493-611)   99% ->   ...261...  0.999, 0.926
    -Group9.12:1372153-1371946  (612-819)   100% ->   ...94...  0.782, 0.996
    -Group9.12:1371851-1371607  (820-1064)   99% ->   ...143...  0.999, 0.647
    -Group9.12:1371463-1370905  (1065-1623)   100% ->   ...1069...  0.981, 0.922
    -Group9.12:1369835-1369482  (1624-1977)   99% ->   ...91...  0.998, 0.997
    -Group9.12:1369390-1369334  (1978-2034)   100% ->   ...218...  0.999, 0.995
    -Group9.12:1369115-1368950  (2035-2200)   100% -)   ...139...  0.909, 0.673
    -Group9.12:1368810-1368583  (2201-2428)   100% ->   ...215...  0.979, 0.999
    -Group9.12:1368367-1368238  (2429-2558)   100% ->   ...5785...  0.999, 0.861
    -Group9.12:1362452-1362307  (2559-2704)   100% ->   ...137...  0.999, 0.961
    -Group9.12:1362169-1361902  (2705-2972)   99% ->   ...717...  0.999, 0.994
    -Group9.12:1361184-1361041  (2973-3116)   100% ->   ...137...  0.999, 1.000
    -Group9.12:1360903-1360721  (3117-3299)   100% ->   ...1255...  0.999, 0.999
    -Group9.12:1359465-1359211  (3300-3554)   100% ->   ...791...  0.999, 0.995
    -Group9.12:1358419-1358249  (3555-3725)   99% ->   ...2147...  0.995, 0.966
    -Group9.12:1356101-1355915  (3726-3912)   100% ->   ...602...  0.989, 0.999
    -Group9.12:1355312-1355143  (3913-4082)   100% -]   ...100...  0.114, 0.000
    -Group9.12:1355042-1354947  (4083-4179)   98% ->   ...326...  0.971, 0.986
    -Group9.12:1354620-1354340  (4180-4460)   99% ->   ...78...  0.985, 0.960
    -Group9.12:1354261-1354067  (4461-4655)   99% ->   ...92...  0.995, 0.989
    -Group9.12:1353974-1353867  (4656-4764)   99% ->   ...90...  0.999, 0.995
    -Group9.12:1353776-1353597  (4765-4944)   100% ->   ...102...  0.991, 0.975
    -Group9.12:1353494-1353376  (4945-5063)   99% ->   ...1207...  0.996, 0.994
    -Group9.12:1352168-1351992  (5064-5240)   99% ->   ...66...  0.988, 0.949
    -Group9.12:1351925-1351738  (5241-5428)   99% ->   ...596...  0.998, 0.961
    -Group9.12:1351141-1350729  (5429-5841)   99% ->   ...109...  0.999, 0.900
    -Group9.12:1350619-1350398  (5842-6063)   100% ->   ...203...  0.998, 0.955
    -Group9.12:1350194-1350006  (6064-6252)   100%

=cut


=item chimera update 2011

Paths (2): *** Possible chimera with breakpoint at 526
  Path 1: query 1..526 (526 bp) => genome Scaffold71:1,005,091..1,027,620 (22530 bp)
  Path 2: query 527..1131 (605 bp) => genome Scaffold71:1,027,402..1,028,300 (899 bp)
>ACYPI000251-RA RefSeq transcript XM_001947776 PREDICTED: Acyrthosiphon pisum similar to prolyl 4-hydroxylase alpha sub
unit 1, putative (LOC100158825), partial mRNA; len=1483
Paths (3): *** Possible chimera with exon-exon boundary (sense) at 234 (donor_prob = 0.990, acceptor_prob = 0.959)

=cut

=item parse for match span

gmap 2007:
Paths (1):
  Path 1: query 51--590 (540 bp) => chr scaffold00107:42,869--42,327 (-543 bp)
    cDNA direction: indeterminate
    Genomic pos: dmagna_asm2_090119.fa:1,311,605--1,311,063 (- strand)
    Accessions: scaffold00107:42,327--42,869 (out of 61814 bp)
    Number of exons: 1
    Coverage: 38.1 (query length: 1416 bp, trimmed query: 1..1416)
    Percent identity: 85.5 (464 matches, 76 mismatches, 3 indels, 0 unknowns)
    Non-intron gaps: 0 openings, 0 bases in cdna; 2 openings, 3 bases in genome
    Translation: 52..591 (181 aa)
    Amino acid changes: P5H [60], E14D [87], R21K [108], M33V [144], M40L [165], R88H [309], E89G [312], C158V [51
9], L166I [543], Q175T [570], R178Q [579], S179V [582], delA181 [588]
==============

gmap 2009
  Path 1: query 64..398 (335 bp) => genome scaffold_32:417,582..417,913 (332 bp)
    cDNA direction: indeterminate
    Genomic pos: dpulex_jgi060905.repeatmasked.fa:111,581,071..111,581,402 (+ strand)
    Accessions: scaffold_32:417,582..417,913 (out of 1144766 bp)
    Number of exons: 1
    Coverage: 84.2 (query length: 398 bp)
    Trimmed coverage: 84.2 (trimmed length: 398 bp, trimmed region: 1..398)
    Percent identity: 70.3 (242 matches, 81 mismatches, 21 indels, 0 unknowns)
    Non-intron gaps: 3 openings, 12 bases in cdna; 2 openings, 9 bases in genome
    Translation: 198..395 (63 aa)

Alignments:
  Alignment for path 1:
    +scaffold_3:9192-9422  (1-231)   98%

  Alignment for path 2:
    -scaffold_226:114532-114323  (1-219)   84%
===========

gmap 2011.aug

Paths (2): *** Possible chimera with exon-exon boundary (antisense) at 1071 (dinucl = TT-AG, donor_prob = 0.000, acceptor_prob = 0.907)
  Path 1: query 1..1071 (1071 bp) => genome scaffold_3:29,460,409..29,457,050 (-3360 bp)
    cDNA direction: sense
    Genomic pos: cacao11allasm:145,909,592..145,906,233 (- strand)
    Accessions: scaffold_3:29,457,050..29,460,409 (out of 34397752 bp)
    Number of exons: 7
    Coverage: 29.4 (query length: 3641 bp)
    Trimmed coverage: 29.4 (trimmed length: 3641 bp, trimmed region: 1..3641)
    Percent identity: 99.9 (1070 matches, 1 mismatches, 0 indels, 0 unknowns)
    Translation: 253..924 (223 aa)
    Amino acid changes: 

  Path 2: query 1072..3641 (2570 bp) => genome scaffold_2:38,211,380..38,217,491 (6112 bp)
    cDNA direction: antisense
    Genomic pos: cacao11allasm:109,771,070..109,777,181 (+ strand)
    Accessions: scaffold_2:38,211,380..38,217,491 (out of 42436413 bp)
    Number of exons: 4
    Coverage: 70.6 (query length: 3641 bp)
    Trimmed coverage: 70.6 (trimmed length: 3641 bp, trimmed region: 1..3641)
    Percent identity: 99.8 (2574 matches, 1 mismatches, 4 indels, 0 unknowns)
    Non-intron gaps: 0 openings, 0 bases in cdna; 2 openings, 4 bases in genome
    Amino acid changes: 

Alignments:
  Alignment for path 1:

    -scaffold_3:29460409-29460038  (1-372)   100% ->   ...695...  0.990, 0.994
    -scaffold_3:29459342-29459210  (373-505)   100% ->   ...280...  0.707, 0.998
    -scaffold_3:29458929-29458859  (506-576)   100% ->   ...101...  0.990, 0.989
    -scaffold_3:29458757-29458679  (577-655)   100% ->   ...85...  0.995, 1.000
    -scaffold_3:29458593-29458484  (656-765)   100% ->   ...657...  0.994, 1.000
    -scaffold_3:29457826-29457762  (766-830)   100% ->   ...471...  0.995, 0.994
    -scaffold_3:29457290-29457050  (831-1071)   99%

  Alignment for path 2:

    +scaffold_2:38211380-38211658  (1072-1346)   98% <-   ...1402...  0.997, 0.999
    +scaffold_2:38213061-38214311  (1347-2597)   99% <-   ...1808...  0.995, 0.993
    +scaffold_2:38216120-38216516  (2598-2994)   100% (-   ...328...  0.838, 0.833
    +scaffold_2:38216845-38217491  (2995-3641)   100%

=item more changes

gmap v2011.08 : another bug, no strand when have 10 exons !!

    cDNA direction: indeterminate << ??? is this bad call?

>cacao3v1Svsc3Loc682t6 nt=10; cf=0.714; len=2064
Paths (1):
  Path 1: query 1..2064 (2064 bp) => genome scaffold_3:24,875,681..28,250,050 (3374370 bp)
    cDNA direction: indeterminate
    Genomic pos: cacao11allasm:141,324,864..144,699,233 (+ strand)
    Accessions: scaffold_3:24,875,681..28,250,050 (out of 34397752 bp)
    Number of exons: 10
    Coverage: 100.0 (query length: 2064 bp)
    Trimmed coverage: 100.0 (trimmed length: 2064 bp, trimmed region: 1..2064)
    Percent identity: 96.3 (2031 matches, 18 mismatches, 60 indels, 0 unknowns)

  Alignment for path 1:

    +scaffold_3:24875681-28246054  (1-453)   96% <-   ...103...  0.991, 0.997  # BAD SPAN
    +scaffold_3:28246158-28246219  (454-513)   82% <-   ...1549...  0.999, 0.935
    +scaffold_3:28247769-28247837  (514-577)   92% <-   ...97...  0.997, 0.988
    +scaffold_3:28247935-28248131  (578-772)   98% <-   ...127...  0.946, 0.974
    +scaffold_3:28248259-28248371  (773-884)   95% <-   ...81...  0.984, 0.992
    +scaffold_3:28248453-28248525  (885-956)   95% <-   ...108...  0.962, 0.981
    +scaffold_3:28248634-28248678  (957-998)   93% <-   ...90...  0.999, 0.978
    +scaffold_3:28248769-28248876  (999-1104)   98% <-   ...87...  0.562, 0.000
    +scaffold_3:28248964-28249556  (1105-1695)   95% <-   ...119...  0.989, 0.999
    +scaffold_3:28249676-28250050  (1696-2064)   97%


#** bad span inner exon ; cant be corrected to reasonable mrna span; keep as problem?

>cacao3v1Svsc9Loc544t31 nt=35; cf=0.618; len=5323
Paths (1):
  Path 1: query 1..5327 (5327 bp) => genome scaffold_9:141,007..37,832,529 (37691521 bp)
    cDNA direction: antisense
    Genomic pos: cacao11allasm:304,473,609..342,165,131 (+ strand)
    Accessions: scaffold_9:141,007..37,832,529 (out of 42035188 bp)
    Number of exons: 9
    Coverage: 100.0 (query length: 5327 bp)
    Trimmed coverage: 100.0 (trimmed length: 5327 bp, trimmed region: 1..5327)
    Percent identity: 97.0 (4956 matches, 19 mismatches, 132 indels, 92 unknowns)
    Non-intron gaps: 15 openings, 72 bases in cdna; 40 openings, 60 bases in genome
    Translation: 4931..4006 (306 aa)
    Amino acid changes: R52G [4777]

Alignments:
  Alignment for path 1:

    -scaffold_9:141007-140312  (1-729)   90% <-   ...258...  0.992, 0.953
    -scaffold_9:140053-139839  (730-936)   82% <-   ...143...  0.988, 0.992
    -scaffold_9:139695-37836232  (937-3680)   99% ==   ...753...   ***query_skip:191***  0.001, 0.000
        ^^^ bad span; 1st 3 exons are 37 MB away from last 6 a problem.
    -scaffold_9:37835478-37835081  (3872-4274)   98% <-   ...116...  0.976, 0.996
    -scaffold_9:37834964-37834876  (4275-4363)   100% <-   ...89...  0.333, 0.986
    -scaffold_9:37834786-37834722  (4364-4428)   100% <-   ...698...  0.995, 0.895
    -scaffold_9:37834023-37833907  (4429-4545)   100% <-   ...200...  0.962, 0.998
    -scaffold_9:37833706-37833622  (4546-4630)   100% <-   ...389...  0.988, 0.070
    -scaffold_9:37833232-37832529  (4631-5327)   98%

#.........
gmap v2011.08 : have some bad exon spans (e.g. 500 bp exon matches > 20 Mb of genome w/o intron)
  -- solve by check that exon genome span =~ transcript span

>P1_g00863t00001  org=isotig02360  gene=isogroup00863  length=1285  numContigs=2
Paths (1):
  Path 1: query 9..1285 (1277 bp) => genome scaffold_3:31,187,299..1,743,136 (-29444164 bp)
    cDNA direction: sense
    Genomic pos: cacao11allasm:147,636,482..118,192,319 (- strand)
    Accessions: scaffold_3:1,743,136..31,187,299 (out of 34397752 bp)
    Number of exons: 5
    Coverage: 99.4 (query length: 1285 bp)
    Trimmed coverage: 99.4 (trimmed length: 1285 bp, trimmed region: 1..1285)
    Percent identity: 99.1 (1271 matches, 3 mismatches, 8 indels, 0 unknowns)
    Non-intron gaps: 1 openings, 8 bases in cdna; 0 openings, 0 bases in genome
    Translation: 10..432 (140 aa)
    Amino acid changes: 

Alignments:
  Alignment for path 1:

    -scaffold_3:31187299-31187263  (9-45)   100% ->   ...190...  0.928, 0.971
    -scaffold_3:31187072-31187014  (46-104)   100% -)   ...406...  0.924, 0.278
    -scaffold_3:31186607-31186526  (105-186)   100% ->   ...88...  0.998, 0.936
    -scaffold_3:31186437-31186366  (187-258)   100% ->   ...263...  0.994, 0.994
    -scaffold_3:31186102-1743136  (259-1285)   98%   << BAD SPAN

>P1_g02609t00001  org=isotig05191  gene=isogroup02609  length=1516  numContigs=1
Paths (1):
  Path 1: query 1..1500 (1500 bp) => genome scaffold_1:15,795,847..37,239,716 (21443868 bp)
    cDNA direction: sense
    Genomic pos: cacao11allasm:16,511,377..37,955,246 (+ strand)
    Accessions: scaffold_1:15,795,847..37,239,716 (out of 38988864 bp)
    Number of exons: 3
    Coverage: 98.9 (query length: 1516 bp)
    Trimmed coverage: 98.9 (trimmed length: 1516 bp, trimmed region: 1..1516)
    Percent identity: 99.5 (1503 matches, 3 mismatches, 4 indels, 0 unknowns)
    Non-intron gaps: 0 openings, 0 bases in cdna; 1 openings, 4 bases in genome
    Translation: 3..1205 (400 aa)
    Amino acid changes: 

Alignments:
  Alignment for path 1:

    -scaffold_1:15795847-15795210  (1-638)   100% ->   ...210...  0.973, 0.961
    -scaffold_1:15794999-15794745  (639-893)   100% ->   ...323...  0.995, 0.987
    -scaffold_1:15794421-37239716  (894-1500)   98%   << BAD SPAN

--------------------

gmap_to_gff.pl is off for gmap.v2011oct:
  sense=-1 bug: rev-rev for CDS
  
 grep ccng06240t00001 ../estgff_mar3/cgb*CC*ff
scaffold_4      cgba100.CCN51   mRNA    17119988        17124134        99      +       .       
  ID=ccng06240t00001;Target=ccng06240t00001 1 2104;aalen=529;cov=100.0;match=2094;nexon=4;pid=99.5;qlen=2104;sense=-1

scaffold_4      cgba100.CCN51   exon    17119988        17120438        99      +       .       Parent=ccng06240t00001;Target=ccng06240t00001 1654 2104;ix=4
scaffold_4      cgba100.CCN51   CDS     17119988        17120437<       99      +       .       Parent=ccng06240t00001
#i.scaffold_4   cgba100.CCN51   intron  17120439        17121245        99      +       .       Parent=ccng06240t00001

scaffold_4      cgba100.CCN51   exon    17121246        17121946        99      +       .       Parent=ccng06240t00001;Target=ccng06240t00001 953 1653;ix=3
scaffold_4      cgba100.CCN51   CDS     17121246        17121946        99      +       .       Parent=ccng06240t00001
#i.scaffold_4   cgba100.CCN51   intron  17121947        17122035        95      +       .       Parent=ccng06240t00001

scaffold_4      cgba100.CCN51   exon    17122036        17122191        100     +       .       Parent=ccng06240t00001;Target=ccng06240t00001 797 952;ix=2
scaffold_4      cgba100.CCN51   CDS     17122036        17122191        100     +       .       Parent=ccng06240t00001
#i.scaffold_4   cgba100.CCN51   intron  17122192        17123338        99      +       .       Parent=ccng06240t00001

scaffold_4      cgba100.CCN51   exon    17123339        17124134        99      +       .       Parent=ccng06240t00001;Target=ccng06240t00001 1 796;ix=1
scaffold_4      cgba100.CCN51   CDS     17123852<       17124134        99      +       .       Parent=ccng06240t00001

>g06240t00001  org=isotig14911  gene=isogroup06240  length=2104  numContigs=1
Paths (1):
  Path 1: query 1..2104 (2104 bp) => genome scaffold_4:17,124,134..17,119,988 (-4147 bp)
    cDNA direction: antisense
    Genomic pos: cacao11allasm:169,606,432..169,602,286 (- strand)
    Accessions: scaffold_4:17,119,988..17,124,134 (out of 33492547 bp)
    Number of exons: 4
    Coverage: 100.0 (query length: 2104 bp)
    Trimmed coverage: 100.0 (trimmed length: 2104 bp, trimmed region: 1..2104)
    Percent identity: 99.5 (2094 matches, 10 mismatches, 0 indels, 0 unknowns)
    Translation: 2103..514 (529 aa)
    Amino acid changes: R346G [1067], S274P [1283]

Alignments:
  Alignment for path 1:

    -scaffold_4:17124134-17123339  (1-796)   99% <-   ...1147...  0.996, 0.993
    -scaffold_4:17122191-17122036  (797-952)   100% (-   ...89...  0.997, 0.905
    -scaffold_4:17121946-17121246  (953-1653)   99% <-   ...807...  0.993, 0.997
    -scaffold_4:17120438-17119988  (1654-2104)   99%

=cut
