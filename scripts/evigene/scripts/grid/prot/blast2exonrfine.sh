# blast2exonrfine.sh : blastx locate proteins, exonerate with region refinement 
# a shell script for gg_job processing

# fixme: need ncbi_bin variable, exonerate bin_dir variable
ncbi_bin=$HOME/bio/ncbi/bin
protbase=`basename $proteins`

$ncbi_bin/blastall -p blastx  \
  -m 9 -e 0.00001 -v 999999 -b 999999 \
  -d $proteins \
  -i $part_dir/$genome \
  -o $part_dir/$protbase.blastx

## fixme: filter lower scored hits at same location; 
cat $part_dir/$protbase.blastx |\
perl -ne 'if(/^\w/){ ($gref,$pid)=split; print "prot:$pid\n";}' |\
sort | uniq | cat - $proteins | perl -ne\
'if(/^>(\S+)/){$d=$1; $p=$d{$d}?1:0;}elsif(/^prot:(\S+)/){$d{$1}++;$p=0;} print if($p);' \
> $part_dir/$protbase.blastx.aa

$bin_dir/exonerate \
  --model protein2genome \
  --refine region  --refineboundary 5000 \
  --minintron 20 --maxintron 25000 \
  --showtargetgff --showvulgar 0 --showalignment 0 \
  --ryo '#qi %qi length=%ql alnlen=%qal\\n#ti %ti length=%tl alnlen=%tal\\n' \
  --query  $part_dir/$protbase.blastx.aa \
  --target $part_dir/$genome \
  > $part_dir/$output_file_name.gff

# fixme: the ryo \\n were not converted to newlines
# perl -pi -e's,\\n,\n,g;' output.gff

