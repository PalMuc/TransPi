# this is shell script for ll_snapjob processing

proteins=/N/u/gilbertd/BigRed/scratch/chrs/aug/dmoj_glean1_good.aa
partpath=$part_dir
outname=$output_file_name   # dmoj_gleangenes
myfasta=$genome            # dyak_caf051213.fa

/N/soft/linux-sles9-ppc64/ncbi-2.2.16/bin/blastall -p blastx  \
  -m 9 -e 0.001 -v 999999 -b 999999 \
  -d $proteins \
  -i $partpath/$myfasta \
  -o $partpath/$outname.blastx

cat $partpath/$outname.blastx |\
perl -ne 'if(/^\w/){ ($gref,$pid)=split; print "prot:$pid\n";}' |\
sort | uniq | cat - $proteins | perl -ne\
'if(/^>(\S+)/){$d=$1; $p=$d{$d}?1:0;}elsif(/^prot:(\S+)/){$d{$1}++;$p=0;} print if($p);' \
> $partpath/$outname.blastx.aa

/N/u/gilbertd/BigRed/bio/exonerate/bin/exonerate --model protein2genome \
  --minintron 20 --maxintron 5000 \
  --showtargetgff --showvulgar 0 --showalignment 0 \
  --ryo '#qi %qi length=%ql alnlen=%qal\n#ti %ti length=%tl alnlen=%tal\n' \
  --query  $partpath/$outname.blastx.aa \
  --target $partpath/$myfasta \
  > $partpath/$outname.exonr.gff

#.......... end of job ...........
# echo $aug/scripts/ll_snapjob.pl  --job script --bin $aug/scripts --script blast2exonr.bash \
#  --out blast2exr --part partitions.list --genome *.fa  --ACCOUNT 'TG-xxxxx' --debug


