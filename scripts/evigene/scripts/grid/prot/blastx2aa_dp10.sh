# blastx2aa_dp10.sh: bash shell script for gg_job processing
#protbase=daph10arp2 #protbase=daph10arp2c_bx  ## FIXME protbase=daph10arp2c_bx2
protbase=daph10arp2d_bx

cat $part_dir/$protbase.blastx | \
$bin_dir/blast92gff3.pl -swap -sorted -minbits 50 -lowscore 0.5 -target -match |\
$bin_dir/overbestgene2.perl -in stdin -pctover 75 |\
perl -ne 'if(/^\w/){ ($pid)=m/Target=([^;\s]+)/;  print "protid:$pid\n" if $pid;}' |\
sort | uniq | cat - $proteins | perl -ne\
'if(/^>(\S+)/){$d=$1; $p=$d{$d}?1:0;}elsif(/^protid:(\S+)/){$d{$1}++;$p=0;} print if($p);' \
> $part_dir/$protbase.blastx.aa

