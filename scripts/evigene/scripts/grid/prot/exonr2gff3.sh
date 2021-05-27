# exonr2gff3.sh:
cat $part_dir/$output_file_name.gff | perl -pe's,\\n,\n,g;' |\
$bin_dir/process_exonerate_gff3.perl -sou exonr > $part_dir/$output_file_name.gff3

