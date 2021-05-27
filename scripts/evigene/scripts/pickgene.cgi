#!/usr/bin/perl

use strict;

## need dsn = aphid2x.conf
my $GBROOT="/bio/argos/work/ggb169/databases/aphid2";
my $JAVACP=" /bio/argos/work/ggb169/lib/java/:/bio/argos/work/ggb169/lib/java/lucene.jar";
my $luindex="$GBROOT/aphid2prot,aphid2pasa,aphid2gene,aphid2rna3,aphid2rna3b,aphid2gene6,aphid2gene7";

## ** this one has to be web writable **
my $updatesindex="$GBROOT/webupdates";

## here result in url params: QUERY_STRING="org=aphid2x&gid=m7AUGepir16bs1017g40t1&gpick=best" 
my $qs= $ENV{QUERY_STRING};
my $ref= $ENV{HTTP_REFERER}; # redirect to "http://server7.eugenes.org:8091/gbrowse/cgi-bin/gbrowse/aphid2x/" 
my $fromad=$ENV{REMOTE_ADDR}; ## eg ="69.136.4.212" 

exit(-1) unless($qs);

my($org,$id,$act)=("","","");
my @qs= split /[\&\;]/, $qs;
map{ $_="mangled" if(length($_)>99 or /[^\w=%,.+-]/); } @qs;
foreach my $q (@qs) { my($k,$v)= split"=", $q; 
  $org=$v if($k eq "org"); $id=$v if($k =~ /id/); $act=$v if($k =~ /pick/);
}
$org||="org";
$ref =~ s/\?.*//;

my $date=localtime();  $date=~s,\s+,/,g;

#my $comment="?comment=PickGene:$id:$act";
#print "Location: $ref$comment\n"; # this works
#print "Content-type: text/html\n\n";
#print "comment=PickGene:$id:$act\n";

open(SAVETO,">>/var/tmp/pickgene-$org.results");
print SAVETO join("\t","#u",$date,$fromad,@qs),"\n";
close(SAVETO);

my @lu= lucegene_get($id, $act); # can recompute index from above table of ids
lucegene_index( \@lu) if @lu;

my $comment="?comment=PickGene:$id:$act";
print "Location: $ref$comment\n"; # this works
print "Content-type: text/html\n\n";
print "comment=PickGene:$id:$act\n";

sub lucegene_get {
  my( $id, $action)= @_;
  my @ret;
  return unless($id and $action);
  my $cmd="echo 'find id:$id parent:$id' | java -cp $JAVACP LuceneSearcher2 -index $luindex |";
  # warn "# CMD: $cmd\n"; #DEBUG
  open(LU, $cmd ) or return; ## warn "# ERROR: pickgene : bad search on $id";
  while(<LU>) { 
    next unless(/^\w/);
    s/$/\tupaction=$action/;
    push @ret, $_;
  }
  return @ret;
}

sub lucegene_index {
  my($aref)= @_;

  my @gcol = qw( ref start stop source method score strand phase 
  gclass gname tstart tstop feature_id group_id bin);

  my $UPSOURCE="USERMOD";

  my $ludata=$ENV{TMPDIR} || "/tmp";
  my $uperr="/dev/null"; ## "$ludata/luperr.$$"; ## or /dev/null
  $ludata.="/lupdata.$$";
  open(LU,">$ludata") or return;
  foreach my $row (@$aref) {
    chomp($row);
    my @val= split"\t",$row;
    $val[3]= $UPSOURCE;
    my $type="$val[4]:$UPSOURCE"; # method:source 
    print LU "i type:$type\n"; 
    for my $i (0..$#gcol) { 
       my ($g,$v)=($gcol[$i],$val[$i]);
       print LU  "i $g:$v\n" if $v;
    }
    $row = join"\t",@val;
    print LU "t contents:$row\n"; 
    print LU "add\n"; 
  } close(LU);

      ## create ONE TIME ONLY
  my $create= (-f "$updatesindex/segments.gen") ? "" : "-create";
  my $ixcmd = "java -cp $JAVACP LuceneIndexer -index $updatesindex $create -data $ludata 1> $uperr 2>&1";
  ## warn "#CMD: $ixcmd\n"; ##DEBUG
  my $res=`$ixcmd`; # does STDERR outpt ##  system($ixcmd);
  unlink($ludata);
  ##return $err;
}

__END__

## get full gff this way and add map view of user curated genes
# ( echo 'find id:mu7AUGepir10p1s2g88t1' ; ) | java -cp /bio/argos/work/ggb169/lib/java/:/bio/argos/work/ggb169/
#  lib/java//lucene.jar LuceneSearcher2 -index /export/udisk2/work/aphid/map/aphid2gene7
## contents are gff db fields:  @gff_vals{@gff_fields}= (
#      $ref, $start, $stop, $source, $method, $score, $strand, $phase, 
#      $gclass, $gname, $target_start, $target_stop, $fid, $gid, $bin);
# + attribute fields: id, name, ...

