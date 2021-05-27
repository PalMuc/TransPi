#!/usr/bin/env perl
# pindel2gff.pl

while(<>) {
  if(/^###/) { $hd=1; }
  elsif($hd and (/^\d/ or /^ChrID/)) { 
    $hd=0; chomp; $in++;
    @v=split"\t";
    
    # damn, not fixed num cols, need parse NT|ChrID.. col prefixes !!
    if($v[0] =~ /^ChrID/) { # fixBP >> BRK
      my($chr,$bp,$up,$dn)=@v[0..2]; $dn||=""; if($up=~/^-/) { ($dn,$up)=($up,$dn); }
      my @samp=splice(@v,3); my $ns=0; # @samp;
      map{ my($sn,$sc)=split; $ns++ if($sc>0); } @samp;  ## ccn51 3 pound7 0  tsh11 1  xxx 2
      ## BP here is chomped below for other BP field !!
      @v= ($in-1,"BK 1","NT 1 ",$chr,"BP ".$bp,$bp,"BP_range ".$bp,$bp,
        "Supports 1",1,$up,1,$dn,1,"S1 1","SUM_MS 100",1,"NumSupSamples ".$ns,$ns,@samp);
        
    } elsif($v[1] eq "LI") { # fixLI
      my($i,$typ,$chr,$bp,$up,$be,$dn)=@v[0..6]; 
      my @samli=splice(@v,7); my @samp; my $ns=0; # @samp; # wrong
      map{ my($sn,$u,$nu,$d,$nd)=split; my $nx=($nu<$nd)?$nd:$nu; 
          push @samp, "$sn $nx";  $ns++ if($nx>0);} @samli; ## ccn51 + 3 - 2  pound7 + 0 - 0 
      # LI length is unknown > read length, use 1 or readlen?
      @v= ($i,$typ." 1","NT 1 ",$chr,"BP ".$bp,$be,"BP_range ".$bp,$be,
        "Supports 1",1,$up,1,$dn,1,"S1 1","SUM_MS 100",1,"NumSupSamples ".$ns,$ns,@samp);
    }
    
    @samp=splice(@v,19);
    map{ s/^(NT|ChrID|BP_range|BP|Supports|S1|SUM_MS|NumSupSamples)\s+//g; } @v;
    my($i,$TYP,$NT,$CHR,$BP,$BPe,$BPR,$BPRe,$SUP,$SUPU,
       $UP,$UPu,$DN,$DNu,$SCO,$MAPQ,$NSA,$NSS,$NSSu)= @v; #19 col before samp
    $samp= join",",map{ ($sn,$ru,$uu,$rd,$ud)=split; $rt=$ru+$rd; "$sn:$rt"; } @samp;
    $NT=~s/"//g; $NT=~s/ +/,/g;
    ($TYP,$TYPl)=split" ",$TYP;
    $or="+"; if($BPR>$BPRe) { ($BPR,$BPRe,$or)= ($BPRe,$BPR,"-"); }
    my $id="$TYP$CHR"."i$i"; $id=~s/([Ss]caffold|[Ss]uper)[_]?/s/; $id=~s/([Cc]ontig)[_]?/c/;
    $attr="ID=$id;ns=$NSS;samples=$samp;len=$TYPl;nt=$NT";
    print join("\t",$CHR,"pindel".$TYP,"strucvar",$BPR,$BPRe,$SCO,$or,".",$attr)."\n";
  }
}

