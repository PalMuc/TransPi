#!/usr/bin/env perl
# samintrons.pl

use strict;
use Getopt::Long;

use constant DOT => '.';
use constant HALF_INTRONS => 0; # gsnap 20M10S cigar types, stranded
my $MIN_IDENT = 95;
my $ALLOW_SOFTCLIP=1; # default yes/no ?
my $SOURCE= "rs";
my $insorted=1;
my $debug=1;
my $idebug=0;
my ($ok,$input,$seqsplice);
my %ins;
my @testin;

my $optok= GetOptions(
  "sam|hits|input=s", \$input, 
  "identity|MIN_IDENT=i", \$MIN_IDENT, 
  "source=s", \$SOURCE, 
  "sorted!", \$insorted, 
  "seqsplice!", \$seqsplice, #? format=
  "softclipok!", \$ALLOW_SOFTCLIP, # for MIN_IDENT error skips
  "debug!", \$debug, 
  "idebug=s", \$idebug, 
  #"test=s", \@testin, # chr:b-e ?
  );

die "opt error" unless($optok);


my %testin=();
# @testin = map{ split /[,;]/ } @testin;
# foreach my $t (@testin) {
#   my($c,$b,$e)=split /[,:\-\s]/,$t;
#   $testin{c}{$c}++;
#   $testin{b}{$b}++;
#   $testin{e}{$e}++;
# }

my $inh= *STDIN;
if($input) {
$ok = ($input =~ /.gz$/) ? open($inh,"gunzip -c $input |") 
      : ($input =~ /^(stdin|-)/) ? $inh= *STDIN
      : open($inh,$input);
die "bad -input=$input" unless($ok);
}

my( $nin, $nerr, $nomap, $nstrand ) = 
  readSam($inh);

warn "#samintrons.$SOURCE: nin=$nin, nerr=$nerr, nstrand=$nstrand, nomap=$nomap\n";
  # if($debug);
  
#...................................................

sub readSam 
{
	my($inh)= @_;
  
  my ( $nin, $nerr, $nomap,  $nstranded ) = (0) x 10;
  my $lchr="";

  while(<$inh>) {
  next unless(/^\w/); # comments?
  chomp;
  my $linein=$_;  $nin++;
  my($qid, $flag, $chr, $cloc, $mapq, $cigar, $matechr, $mateloc, $isize, $seq, $qual, @opt)
    = split"\t";

#   my $f_pair= ($flag & 0x0001);
#   my $f_pairok= ($flag & 0x0002);
  my $f_mismatch= ($flag & 0x0004);  # also set when chr eq '*'
#   my $f_mismate= ($flag & 0x0008);
  my $f_isrev= ($flag & 0x0010);
#   my $f_revmate= ($flag & 0x0020);
  
  putIntrons() if($insorted and $lchr and $chr ne $lchr);
  
  $lchr= $chr; 
  do{ $nomap++; next;} if($chr eq '*'); # no map

  my $len  = length($seq);
  my $cend = $cloc; # calc correct end from align/intron spans
  my $intype="intron";
  my(@aln,@intr,@itype);
  my $softclip=0;
  # fixme for: ^nSnM ... was m/^(\d+)M/
#   if( $cigar =~ m/^(\d+)M/ ) { @aln=($1); $cend += $1; }
#   while($cigar =~ m/(\d+)N(\d+)M/g) { my($bi,$bx)=($1,$2); push(@intr,$bi); push(@aln,$bx); $cend += $bi + $bx; }

  # Op BAM Description : SAM Cigar ops, 2010
  # M 0 alignment match (can be a sequence match or mismatch) 
  # I 1 insertion to the reference 
  # D 2 deletion from the reference 
  # N 3 skipped region from the reference == intron
  # S 4 soft clipping (clipped sequences present in SEQ) 
  # H 5 hard clipping (clipped sequences NOT present in SEQ) 
  # P 6 padding (silent deletion from padded reference) 
  # = 7 sequence match    << M eqiv
  # X 8 sequence mismatch << M eqiv
  # ** H is treated odd cuz seq is clipped to it
  # ** S softclip also odd, $cloc is at NON-clipped seq start.. and cend = cloc + seqlen - end S
  # ** add stats for S, H other cig flags
    # ** H is treated odd cuz seq is clipped to it
    # ** S softclip also odd, $cloc is at NON-clipped seq start.. and cend = cloc + seqlen - end S
    
  if($cigar =~ m/^(\d+)M/ ) { @aln=($1); $cend += $1; }
  elsif($cigar =~ /^(\d+)[HSN]/) { @aln=(0); } # need for nSnM cigar
  while($cigar =~ m/(\d+)([DIHNSP])(\d+)M/g) { 
    my($bi,$bt,$bx)=($1,$2,$3);
    $softclip+=$bi if($bt eq "S");
    $bi=0 if($bt eq "H" or $bt eq "S" or $bt eq "I");  # 2012.7 I fix
    push(@intr,$bi);  push(@itype,$bt);
    push(@aln,$bx); 
    $cend += $bi + $bx; 
    }
  # if($cigar =~ /(\d+)[HS]$/) { $lenc+= $1; } # end

# 2012.7.2 : Insert bug for intron position  1-I shifts intron by +1

#  # 2010.7: add gsnap half-intron (S) from cigar?? 11S25M; 27M9S; doesnt look right...
#   if(HALF_INTRONS and $cigar =~ m/S/){
#   if( $cigar =~ m/^(\d+)M(\d+)S/ ) {  my($al,$bi,$bx,$off)=($1,2,1,$2); 
#       @aln=($al); $cend += $al; push(@intr,$bi); push(@aln,$bx); $intype="inhalf5";}
#   elsif( $cigar =~ m/^(\d+)S(\d+)M/ ) { my($al,$bi,$bx,$off)=($2,2,1,$1); 
#       @aln=(1); $cend += $al + $off; push(@intr,$bi); push(@aln,$al); $intype="inhalf3";}
#   }

  my $instrand= DOT; my $inmismat=0;
  my $score= DOT; 

  my $nx= scalar(@aln);
  my $alen=0; map{ $alen+=$_ }@aln; 
  foreach (@opt) { 
    if(m/NM:i:(\d+)/){ $alen -= $1; } 
    elsif(m/XS:A:([\+\-])/){ $instrand=$1;  $nstranded++; } 
    elsif(m/NS:i:(\d+)/) { $inmismat=$1; } #new, if XS, = Mismatches within min_anchor_len of a splice junction, > 0 is poor
    # elsif(m/RG:A:(\S+)/){ $readgroup=$1; } #new
    }
    
  if($ALLOW_SOFTCLIP and $softclip>0) {
  $score=int(0.5 + 100*$alen/($len-$softclip));  
  } else {
  $score=int(0.5 + 100*$alen/$len);
  }
  
  do{ $nerr++;  next;} if($score < $MIN_IDENT); # 2err/37bp = 35/37 = 95%; 36/37 = 97%

  my ($sid, $matesid, $matestrand)=(0,0,DOT);
  my $xb= $cloc; my $lin=0; my $qb=0;
  for (my $i=0; $i < $nx; $i++) {
    my $aln= $aln[$i]; 
    my $intr=$intr[$i]||0;
    my $isintron= ($itype[$i] eq "N")?1:0;

    if($aln > 1 and $isintron and $intr>0) {
      my $xe= $xb + $aln - 1; # this is right
      my ($sb,$se);
      $sb= $xe + 1;  $se= $xe + $intr; 
      $qb += $aln;

my $SPLICELEN = 20;      
      if($seqsplice) {
        my $ab= $qb - $SPLICELEN; $ab=0 if $ab<0;
        my $aw= $qb - $ab;
        my $bb= $qb;
        my $bw= ($i+1 < $nx) ? $aln[ $i+1] : 0; $bw=$SPLICELEN if($bw>$SPLICELEN);
        my $ssb= substr($seq,$ab,$aw); $aw=length($ssb); 
          $ssb= substr("nnnnnnnnnnnn",0,$SPLICELEN-$aw).$ssb if($aw<$SPLICELEN);
        my $sse= substr($seq,$bb,$bw); $bw=length($sse); 
          $sse.=substr("nnnnnnnnnnnn",0,$SPLICELEN-$bw) if($bw<$SPLICELEN);
        my $ins= $ssb .".". $sse;
        if($instrand eq "-") { $ins=reverse($ins); $ins=~tr/acgtACGT/tgcaTGCA/; }
        print $ins,"\n";

      } else {   
        # addIntron( $chr, $sb, $se, $instrand, $intype);  
my $inatt= ($idebug==1) ? "$intype.$qid" : $intype;
        my $in= join("\t",$chr,$sb,$se,$instrand,$inatt);
        $ins{$in}++;
      }
#       print "# $chr:$sb-$se:$instrand/$f_isrev\n# $linein\n" 
#         if($debug and $testin{c}{$chr} and $testin{b}{$sb} and $testin{e}{$se});
    
    }
    
    $xb += $aln + $intr ; # this is right, not offby1
    $lin= $intr;
 	  }
 	  
  }
  
  putIntrons();
  return( $nin, $nerr, $nomap,  $nstranded );
}

# sub addIntron
# {
#   my( $chr, $tb, $te, $tor,$typ)= @_;  
#   my $in= join("\t",$chr,$tb,$te,$tor,$typ);
#   $ins{$in}++;
# }

sub putIntrons {
  foreach my $in (sort keys %ins) {
    my $n=$ins{$in};
    my @v=split"\t",$in;
    my $typ=$v[4]||"intron";
    my $vend=".";
    if($debug) { my $w=1+$v[2] - $v[1]; $vend.="\tw=$w"; }
    print join("\t",$v[0],$SOURCE,$typ,$v[1],$v[2],$n,$v[3],$vend),"\n"; # gff? or bed-like?
  }
  %ins=();
}

# sub _strandlet  { my $s=$_[0]; return ($s eq "-") ? "r" : ($s eq "+") ? "f" : "u"; } # want \w char
# 
# sub _locid
# {
# 	my ($tag,$ref,$tb,$te,$to)= @_; 
# 	return join("_",$tag,$ref,$tb,$te,_strandlet($to));
# }
# 
