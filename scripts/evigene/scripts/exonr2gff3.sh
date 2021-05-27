# exonr2gff3.sh: bash shell script for gg_job processing
protbase=plants
cat $part_dir/$protbase-exonr.gff | perl -pe's,\\n,\n,g;' |\
$bin_dir/process_exonerate_gff3.perl -sou exonr > $part_dir/$protbase-exonr.gff3

