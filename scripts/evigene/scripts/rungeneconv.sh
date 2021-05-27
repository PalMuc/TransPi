# rungeneconv.sh; bash shell

cd $part_dir/
for fa in *.fa; {
  $bin_dir/muscle -quiet -clwstrict -in $fa -out $fa.mus
  $bin_dir/geneconv $fa.mus /g2 /w123 /lp >& /dev/null
}

