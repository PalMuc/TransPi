
gunzip -c evg7pubdir/Genes/earlyaccess/dmagset7finloc9b.puban.mrna.gz | \
env end=all perl mrnaparts.pl

# mrnaparts.pl
BEGIN{ $PART=$ENV{end}||"all"; } 
while(<>) {
  if(/^>(\S+)/){ putpart($id,$offb,$offe,$fa) if($id and $fa);
    $id=$1; $fa=""; ($offb,$offe)=m/offs=(\d+)-(\d+)/; ($len)=m/clen=(\d+)/; } 
  else { chomp; $fa.=$_; } 
}
END{ putpart($id,$offb,$offe,$fa) if($id and $fa); } 

sub putpart{ 
my($id,$offb,$offe,$fa)=@_; 
$flen=length($fa); $cdsw=1+$offe-$offb; 
$u5w=($offb>1)? $offb-1 : 0; $u3w=($offe<$flen)?$flen-$offe:0;
($u5)= substr($fa,0,$u5w); 
($cds)=substr($fa,$offb-1,$cdsw);  
($u3)=substr($fa,$offe,$u3w);  
$hdr="cdsoff=$offb-$offe,$flen";
putone($id,"utr5", $hdr, $u5) if($PART=~/utr5|all/);  
putone($id,"cds", $hdr, $cds) if($PART=~/cds|all/); 
putone($id,"utr3", $hdr, $u3) if($PART=~/utr3|all/); 
} 

sub putone{ my($id,$typ,$hdr,$fa)=@_; my $len=length($fa); 
print ">$id.$typ type=$typ len=$len $hdr\n$fa\n" if($len>0); } 