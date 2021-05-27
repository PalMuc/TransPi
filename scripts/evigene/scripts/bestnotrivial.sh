#!/bin/tcsh
# bestnotrivial.sh
## filter out no-evid but paralogs (~5,000)
## beware @ev matches bestgenes scoretype

# * rewrite this to do more, merging old/new genes
# *  as bestgenes_update.pl; best_predictions.pl ; finalize_predictions ; update_predictions; ...

#  -- merge pasa updates (CDS only) to genes, marking
#     .. see genes/addcheckprot.sh
#  -- add pasa alt-tr, better marking of alt-tr ??
#  -- add/update mustkeepgenes from pickgene.results;
#     .. see genes/mustkeepgenes.sh
#     .. make simple enough to do this often
#  -- add/option PublicID numbering
#  -- maybe add merge annotations, e.g. Name=, Note=, desc=, Dbxref=, GOterms=, ...

#  -- get scoretype= from input bestgenes, not fixed here,
#  -- use keepids, keepsource to avoid kicking out those as trivial
#  -- need to know sources.gff to pull new musteepgenes, replacing dropped/overlapped old genes


## FIXME: no rna-evd score on rnagenes: cant drop those however, all have +rna scores
## fix here? or should add "perfect" rna evd score to those
# keep sources:  acypi ars17trinity ars27cuf8 pasa2_aphid3 ref_aphid2

# set bestg=bestgenes_of5.an7i.gff
#scoretype: many.homolog:30,paralog:4,ref:5,est:5,pro:5,rseq:5,intr:20,nintron:100,insplit:-10,terepeat:-1,UTR:1,CDS:2
# ev=(0,2,3,4,5,6,7) # ev6,7 same here

# grep -v '^#x' bestgenes_of5.an7k.gff | perl -pe'if(/^\w/){ s/\t\w+/\tDGILmix7k/; s/(Parent|ID)=/$1=mk7/;}' \
#  > bestgenes_cdsrna.mix7k.gff
# set bestg=bestgenes_of5.an7k.gff
# scoretype: many.homolog:30,paralog:4,ref:5,est:5,pro:5,rseq:5,intr:20,nintron:100,insplit:-10,terepeat:-1,UTR:1,CDS:2

#set bestg=bestgenes_of6.an7m.gff
#set MSRC=DGILmix7m
#set MID=xm7
# scoretype=homolog:7,paralog:1,ref:3,est:3,pro:3,rseq:3,intr:7,nintron:10,inqual:10,terepeat:-1,UTR:1,CDS:2, 

#set bestg=bestgenes_of6.an7n.gff
# scoretype=homolog:7,paralog:1,ref:3,est:2,pro:2,rseq:2,intr:19,nintron:50,inqual:50,terepeat:-1,UTR:1,CDS:1, 

#set bestg=bestgenes_of6.an7p.gff
#set bestg=bestgenes_of6.an7q.gff
#set bestg=bestgenes_of6.an7r.gff
#set MSRC=DGILmix7r
#set MID=mr7

set bestg=bestgenes_of11.an8d.gff
set keeps="acypi|ars17trinity|ars27cuf8|pasa2_aphid3|ref_aphid2"
set keepids=bestgenes.DGILmix8d.addback1.ids
#score=homolog:7,paralog:1,ref:3,est:3,pro:3,rseq:2,intr:10,nintron:40,inqual:15,maqual:5,terepeat:-3,UTR:1,CDS:1, 
set MSRC=DGILmix8d
set MID=md8

cat  $bestg | env keepf=$keepids keeps=$keeps perl -pe \
'BEGIN{ $keeps=$ENV{keeps}||"nzAtAll"; $kf=$ENV{keepf}; if($kf and -f $kf) { open(F,$kf); while(<F>){ chomp; $kid{$_}=1;}}} \
if(/^\w/) { if(/\tmRNA/) { $skip=0; \
($sr,$pv)=(split"\t")[1,5]; ($g)=m/ID=([^;\s]+)/; unless($sr =~ m/$keeps/ or $kid{$g}) {  \
@p=map{s,/\d+,,; $_} split",", $pv; \
($al,$ap)=m/aalen=(\d+).(\d+)/; ($ph)=m/pHOBEST=(\d+)/; \
@ev=(0,2,3,4,5,6,7); @ew=(1,1,1,1,1,1,10); \
$ev=0; for $i (@ev) {$ev+=$ew[$i]*$p[$i];} \
if( $ev < 5 ) { $skip=($p[1] > 0 and $al > 120 and $ap >= 66 and $ph >= 66)?0:2; } \
s/$/;skip=$skip/ if ($skip); } } s/^/#x./ if($skip); } ' \
> $bestg.notriv

echo "# Predictor totals after trival removal ......." >> $bestg.notriv
grep -v '^#x' $bestg.notriv | grep mRNA | cut -f2 | sort | uniq -c | sed 's/^/# /' >> $bestg.notriv

# FIXME: preserve source predictor col. as mRNA attribute;
# fixme for alt-tr, pasa updates, mustkeep updates

grep -v '^#x' $bestg.notriv | perl -pe\
 'if(/^\w/){ s/\t\w+/\t'$MSRC'/; s/(Parent|ID)=/$1='$MID'/;}' > bestgenes.$MSRC.gff


#..........
#!/bin/tcsh
# addcheckprot.sh

# env genes=xxxx.gff add=xxx.checkprot.gff  addcheckprot.sh
# ** need to add pro= scores to checkprot.gff CDS
# overlapfilter.perl -strand -pass CDS -pct 10 -act markidbase -mark pro -in stdin -over prot/protein_uniq.gff.gz

if( ! ( -f $genes & -f $add ) ) then
  echo "MISSING genes=$genes or add=$add"; exit -1
endif

set fout=`echo $genes | sed 's/.gff/.ckup.gff/'`

egrep '(mRNA|CDS)' $add |\
scripts/overlapfilter.perl -strand -pass CDS -pct 10 -act markidbase -mark pro \
-in stdin -over prot/protein_uniq.gff.gz |\
sed 's/^/#ck/' | cat - $genes | perl -ne\
'if(s/^#ck//){ if(/\tmRNA/){ ($d)=m/ID=([^;\s]+)/; ($at)=m/(cxlen=.*)$/; $ma{$d}=$at; } \
elsif(/\tCDS/){ ($p)=m/Parent=([^;\s]+)/; s/ID=[^;]+;//; $mc{$p}.=$_; } } \
else { if(/\tmRNA/) { ($d)=m/ID=([^;\s]+)/; $at=""; if( $at= $ma{$d}){ \
s/;?(cxlen|aalen|protein)=[^;\s]+//g; s/$/;$at/; $cd=$mc{$d}; ($o)=(split)[6]; \
if($o eq "-") { chomp($cd);  @cd=reverse split"\n",$cd; $cd=join"\n",@cd; $cd.="\n"; } $_.=$cd if($cd); } } \
elsif($at and /\tCDS/) { ($p)=m/Parent=([^;\s]+)/; if($p ne $d){warn "misorder CDS.$p ne $d\n"; } \
$_ = "" if($mc{$p});  } print; }' \
> $fout

