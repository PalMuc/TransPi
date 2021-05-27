# a bash shell script for gg_job processing
# one-off for acyr-tcas3 blastx 2 exonerate

protbase=acyr-tcas3_gleanaa
# bin_dir=$HOME/bio/augustus/scripts
# proteins=tcas3_gleangenes.aa

cat $part_dir/$protbase.blastx | sort -k2,2 |\
$bin_dir/blast92gff3.pl -swap -minbits 50 -lowscore 0.5 -target -match |\
$bin_dir/overbestgene2.perl -in stdin -pctover 75 |\
perl -ne 'if(/^\w/){ ($pid)=m/Target=([^;\s]+)/; $pid=~s/tcas3_//; print "protid:$pid\n" if $pid;}' |\
sort | uniq | cat - $proteins | perl -ne\
'if(/^>(\S+)/){$d=$1; $p=$d{$d}?1:0;}elsif(/^protid:(\S+)/){$d{$1}++;$p=0;} print if($p);' \
> $part_dir/$protbase.blastx.aa

