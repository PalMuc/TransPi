# blastx.sh : blastx locate proteins
#   -U  Use lower case filtering of FASTA sequence [T/F]  Optional

ncbi_bin=/opt/bio/ncbi/bin
protbase=`basename $proteins`

$ncbi_bin/blastall -p blastx  \
  -U T -m 9 -e 0.00001 -v 999999 -b 999999 \
  -d $proteins \
  -i $part_dir/$genome \
  -o $part_dir/$protbase.blastx

