# blastx2gff.sh: bash shell script for gg_job processing
protbase=arp5
cat $part_dir/$protbase.aa.blastx  | \
$bin_dir/blast92gff3.pl -swap -low 0.9 -sou bx$protbase > $part_dir/$protbase-blastx.gff

