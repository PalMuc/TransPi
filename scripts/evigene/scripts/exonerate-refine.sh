# exonerate-refine.sh : exonerate with region refinement , script for gg_job processing
# change refineb, maxin to speed up .. refineb=500 : 10x faster than =5000;
# preprocess tblastn to exoner.fagff : gff/query.aa/genome.na for each region of interest

$bin_dir/exonerate \
  --model protein2genome \
  --refine region  --refineboundary 500 \
  --minintron 20 --maxintron 10000 \
  --showtargetgff --showvulgar 0 --showalignment 0 \
  --ryo '#qi %qi length=%ql alnlen=%qal\\n#ti %ti length=%tl alnlen=%tal\\n' \
  --query  $part_dir/$proteins \
  --target $part_dir/$genome \
  > $part_dir/$output_file_name.gff

# fixme, other bindir here = $HOME/bio/augustus/scripts
cat $part_dir/$output_file_name.gff | perl -pe's,\\n,\n,g;' |\
$HOME/bio/augustus/scripts/process_exonerate_gff3.perl -sou exonr > $part_dir/$output_file_name.gff3

