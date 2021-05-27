# blast2exonrfine.sh : blastx locate proteins, exonerate with region refinement 

## add exonr --softmasktarget  for genome_mask.fa

# cacao9asm_mask.fa  plant2.aa.blastx  plants.aa.blastx
# group 2: all but arabidopsis
scripts=$HOME/bio/augustus/scripts
protbase=plant3
# in plants.aa + plant2.aa

cat $part_dir/plant?.aa.blastx | \
grep -v 'arabidopsis:' |\
$scripts/blast92gff3.pl -swap -sorted -minbits 50 -lowscore 0.5 -target -match |\
$scripts/overbestgene2.perl -in stdin -pctover 75 |\
perl -ne 'if(/^\w/){ ($pid)=m/Target=([^;\s]+)/;  print "protid:$pid\n" if $pid;}' |\
sort | uniq | cat - $proteins | perl -ne\
'if(/^>(\S+)/){$d=$1; $p=$d{$d}?1:0;}elsif(/^protid:(\S+)/){$d{$1}++;$p=0;} print if($p);' \
> $part_dir/$protbase.blastx.aa

$bin_dir/exonerate \
  --model protein2genome \
  --refine region  --refineboundary 5000 \
  --minintron 20 --maxintron 25000 \
  --showtargetgff --showvulgar 0 --showalignment 0 \
  --ryo '#qi %qi length=%ql alnlen=%qal\\n#ti %ti length=%tl alnlen=%tal\\n' \
  --query  $part_dir/$protbase.blastx.aa \
  --softmasktarget 1 --target $part_dir/$genome \
  > $part_dir/$output_file_name.gff

