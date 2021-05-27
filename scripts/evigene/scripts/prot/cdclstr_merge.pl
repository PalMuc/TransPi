#!/usr/bin/env perl
# cdclstr_merge.pl
# usage:  cat $pt.{ok_cd,poor_cd}.aa.clstr | cdclstr_merge.pl > $pt-merged.aa.clstr
# ?? add keepids=file.ids to leave out useless rows of merge

use strict;

my $nf=1; # infile
my($mainid,$nc,$nf,$lc,%id2main,%clmain,%cls,%keepid,$keepset);
if($ENV{keepid} and -f $ENV{keepid}) {
  open(F,$ENV{keepid}); while(<F>) { my($id)=split; $keepid{$id}=1; $keepset++; } close(F);
}

my(@id,@cl);  
while(<>) {
  if(/^>Cluster\s*(\w+)/) { 
    my $cl=$1; 
    puto($nc,$mainid,\@id,\@cl) if(@cl); # add nf file param?
    $nc++; $nf++ if($cl<$lc); $lc=$cl;  
    @id=@cl=();
  } elsif(/^(\d+)/) {
    my $i=$1;
    if(m/(\d+)(aa|nt), >(.+)\.\.\. (.*)$/) {
      my($aw,$typ,$id,$pi)=($1,$2,$3,$4); 
      # $pi=~s/at //; 
      # s/^(\d+)/${1}f$nf/; # bad location; put file info at end of line, new columns; add matchid?
      s/$/\tf$nf/; # makes problems with parse of m/\*$/ in pinfo ..
      ## change " at xxx" to " atf2 xxx" ?
      my $ok= ($keepset) ? $keepid{$id} : 1;  
      if($pi =~ m/\*/) { $mainid=$id; } # need main if not ok but others ok
      elsif($ok) { push @id, $id; } 
      push @cl, $_ if($ok);  
      }
      # mainid not 1st/unsorted: $id2main{$id}=$mainid;
    } 
}
puto($nc,$mainid,\@id,\@cl) if(@cl); 
  
foreach my $clid (sort{$a<=>$b} keys %cls) { 
  my $oc= $cls{$clid}; 
  print ">Cluster $clid\n"; 
  my $i=0; foreach my $ll (@$oc) { $ll=~s/^\d+/$i/; $i++; print $ll; } 
} 

sub puto { 
  my($clid,$mainid,$ids,$clin)= @_;
  ## my $clid=$nc; 
  if(ref $ids) { map{ $id2main{$_}=$mainid } @$ids; }
  my $onc=$clmain{$mainid} || $clmain{ $id2main{$mainid} } || 0; 
  my $oc; my @cl;
  ## ?? append infile:mainid to each clin line for merge reference?? this.mainid becomes prior altid..
  if($onc and $oc=$cls{$onc}) { $clid=$onc; my %dido=(); 
    @cl= grep/\w/,map{ my($d)=m/ >(.+)\.\./; ($dido{$d}++)?"":$_; } (@$oc,@$clin); 
  } else {
    @cl= @$clin;
  }
  $cls{$clid}= \@cl; $clmain{$mainid}= $clid; 
} 

