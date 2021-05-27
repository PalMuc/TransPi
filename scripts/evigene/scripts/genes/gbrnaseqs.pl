#!/usr/bin/env perl
# gbrnaseqs.pl : ncbi genomes/SPECIES/RNA/rna.gbk parser for gene transcripts, proteins, annots

=item about evigene/scripts/genes/gbrnaseqs.pl

 Evigene parse script to extract cdna, prots from NCBI genome data
 for Evigene genes x NCBI genome genes comparisons
 input source data: ftp.ncbi.nih.gov/genomes/SPECIES/RNA/rna.gbk.gz

=item usage

 gunzip -c $ngeno/M_musculus/RNA/rna.gbk.gz | env rna=1 prefix=mouse ./gbrnaseqs.pl > mouse_ncbi1503.cdna  
 gunzip -c $ngeno/M_musculus/RNA/rna.gbk.gz | env aa=1 prefix=mouse ./gbrnaseqs.pl > mouse_ncbi1503.aa  

 gunzip -c $ngeno/D_rerio/RNA/rna.gbk.gz | env aa=1 prefix=zfish  gbrnaseqs.pl > zfish_ncbi1409.aa  
 gunzip -c $ngeno/D_rerio/RNA/rna.gbk.gz | env rna=1 prefix=zfish gbrnaseqs.pl > zfish_ncbi1409.rna 

=cut 

$dorna=$ENV{rna}||$ENV{cdna}||0;
$doaa= $ENV{prot}||$ENV{aa}|| not $dorna;
$idpre=$ENV{prefix}||""; $idpre.=":" if($idpre);
%alts=();
## add /^LOCUS/ or not? YES, LOCUS always exists, ACC may be empty
## added source moltype for mRNA, ncRNA 
## ncRNA ncRNA_class="lncRNA" mol_type="transcribed RNA" 

while(<>) {
if(/^LOCUS\s+(\S+)/) { $accid=$1; @dbx=();
  $gid=$pid=$locid=$mtype=$infa=$inaa=$incds=$trb=$tre=$cdb=$cde=$tlex=$inam=$nam=$sp=$flags=$gname=0; } 
elsif(/^ACCESSION\s+(\S+)/) { $accid=$1; }
elsif(/^\s+ORGANISM\s+(.+)$/) { $sp=$1; $sp=~s/ /_/g; } 
elsif(m,^\s+(gene|mRNA)\s,) { putprot() if($aa and $doaa); $inam=$inaa=$incds=0; 
  ($trb,$tre)= m/\s([\d\<\>]+)..([\d\<\>]+)/; } 
elsif(m,^\s+CDS\s,) { $incds=1; ($cdb,$cde)= m/\s([\d\<\>]+)..([\d\<\>]+)/;  } 
elsif(m,^\s+/product=.([^"\n]+),) { my $pn=$1;
  my $nok=($nam)?0 :($incds or $dorna)?1 : 0; # rna doesnt have product, except via CDS subrec
  if($nok){ $nam=$pn; $inam=(m/\"$/)?0:1;} } # product name only from CDS record?
#b.elsif(m,^\s+/product=.([^"\n]+),) { $nam=$1; } # BUG, get product name only from CDS record ..
elsif(m,^\s+/mol_type=.([^"\n]+),) { $mtype=$1; } 
elsif(m,^\s+/locus_tag=.([\w\.]+),) { $locid=$1; } 
elsif(m,^\s+/protein_id=.([\w\.]+),) { $pid=$1; } 
elsif(m,^\s+/gene=.(\w+),) { $gname=$1; } 
elsif(m,^\s+/db_xref=.GeneID:(\w+),) { $gid=$1; } 
elsif(m,^\s+/db_xref=(\S+),) { $dx=$1; $dx=~s/"//g; 
  push @dbx,$dx unless($dx=~/^(GI|taxon):/ or ($doaa and not $incds)); } 
# ugh, add /ribosomal_slippage and parse CDS disjunct spans: CDS join(90..182,184..651)  
# eg zfish NM_001080090 gene=oaz2b, antizyme retina, /note="protein translation is dependent on +1 polyamine-induced ribosomal frameshift
elsif(m,^\s+/(ribosomal_slippage|partial|exception),) { my $f=$1; $flags.="$f,"; } 
# add: /transl_except=(pos:3787..3789,aa:OTHER)
elsif(m,^\s+/transl_except=(\S+),) { $tlex=$1; $tlex=~s/pos://g; $tlex=~s/[\(\)]//g; }
elsif(m,^\s+/translation=.([^"\n]+),) { $inaa=1; $aa=$1; $inaa=0 if(m/\"$/); }
# FIXME1804.one line: /translation="MERSTQELFINFTVVLITVLLMWLLVRSYQY" ; inaa=0
elsif(/^ORIGIN/) { putprot() if($aa and $doaa); $inam=$inaa=0; $infa=1; $fa=""; } 
elsif(m,^/,) { putrna() if($fa and $dorna); }
elsif($infa and /^\s*\d+\s+(.+)/) { $seq=$1; $seq=~s/ //g; $fa.=$seq; }
elsif($inaa) { $inaa=0 if(/"/);  s/["\s\d]+//g; $aa.=$_; } # multiline prot in attr
elsif($inam) { $inam=0 if(/"/);  s/^\s*/ /; s/["\n]*$//; $nam.=$_; } # multiline in attr
elsif(m,^     \w,){ $incds=0; } # '     mat_peptide  ' any ^5spaceKEY cancels last 5key. . has name
}

sub putrna { push @dbx,"ncbiprot:$pid" if($pid); putseq($accid, $fa, $mtype||"rna"); $fa=""; }
sub putprot { if($accid){ $pid||=$accid; push @dbx,"ncbirna:$accid"; } putseq($pid, $aa,"protein"); $aa=""; }

sub putseq { 
  my($sid,$aa,$stype)=@_; my($aw,$genev,$klen,$dbxv); 
  return 0 unless($aa); $aw=length($aa); $aa=~s/(.{60})/$1\n/g; $nam=~s/;/,/g;
  $gid= $gid || $locid || $sid || $gname; #..solves cases of no gene/geneid.
  $ialt=++$alts{$gid};
  $genev= ($gname)?" gene=$gname;" : "";
  my $offs= "";
  if($cde) { 
    #below# $offs=" offs=$cdb-$cde;"; #!?? drop <> symbols in offs, use full/partial flag
    $trb=~s/[<>]//; $tre=~s/[<>]//; 
    my $aac=3; # aacomplete, dont need has partial flag
    if($cdb=~s/[<>]//) { $aac -= 1; } if($cde=~s/[<>]//) { $aac -= 2; } 
    my $full=($aac==3)?"complete":($aac==1)?"partial3":($aac==2)?"partial5":"partial";
    my $pcds= ($tre>$trb) ? int(0.5 + 100*($cde-$cdb)/($tre-$trb)) : 99;
    $aw.=",$pcds%,$full"; 
    $offs=" offs=$cdb-$cde;"; #! drop <> symbols in offs, use full/partial flag
  }  
  $offs .= " orfexcept=$tlex;" if($tlex);
  my $flagv=($flags)?" flags=$flags;":"";
  my $dbx="";
  if(@dbx) { my %dbx=map{ $_,1 } grep /\w/, @dbx;  $dbx=join",",sort keys %dbx; }
  $dbxv = ($dbx=~/\w/)?" db_xref=$dbx;":"";
  $klen = ($stype=~/^aa|prot/)?"aalen":"len"; $stype=~s/ /_/g; # ncRNA == "transcribed RNA"
  print
">$idpre$sid type=$stype; $klen=$aw;$genev geneid=$gid;$dbxv isoform=$ialt;$offs$flagv Name=$nam; species=$sp;\n$aa\n";
  return 1;
}


__END__

=item input ncbi.rna.gbk 

 Mouse15EVm000072t2 = XM_006544295 but much shorter prot.
 ** Note this CDS/Prot has inner-stop ignored by NCBI/Gnomon model **
 ** ?? add to Evigene option to auto detect/correct such also?
    /transl_except=(pos:3787..3789,aa:OTHER)
 *? add /note info 
 ** add CDS span == evg standard key offs=499-12021
 
LOCUS       XM_006544295           12067 bp    mRNA    linear   ROD 09-FEB-2015
DEFINITION  PREDICTED: Mus musculus dynein, axonemal, heavy chain 3 (Dnah3),
            transcript variant X1, mRNA.
ACCESSION   XM_006544295
VERSION     XM_006544295.2  GI:755479576
DBLINK      BioProject: PRJNA16113
KEYWORDS    RefSeq; corrected model; includes ab initio.
SOURCE      Mus musculus (house mouse)
  ORGANISM  Mus musculus
            Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;
            Mammalia; Eutheria; Euarchontoglires; Glires; Rodentia;
            Sciurognathi; Muroidea; Muridae; Murinae; Mus.
COMMENT     MODEL REFSEQ:  This record is predicted by automated computational
            analysis. This record is derived from a genomic sequence
            (NW_001030877.1) annotated using gene prediction method: Gnomon,
            supported by mRNA and EST evidence.
            Also see:
                Documentation of NCBI's Annotation Process
            
            On Feb 9, 2015 this sequence version replaced gi:568894745.
            
            ##RefSeq-Attributes-START##
            ab initio            :: 1% of CDS bases
            internal stop codons :: corrected 1 genomic stop codon
            ##RefSeq-Attributes-END##
            
            ##Genome-Annotation-Data-START##
            Annotation Provider         :: NCBI
            Annotation Status           :: Full annotation
            Annotation Version          :: Mus musculus Annotation Release 105
            Annotation Pipeline         :: NCBI eukaryotic genome annotation
                                           pipeline
            Annotation Software Version :: 6.2
            Annotation Method           :: Best-placed RefSeq; Gnomon
            Features Annotated          :: Gene; mRNA; CDS; ncRNA
            ##Genome-Annotation-Data-END##
FEATURES             Location/Qualifiers
     source          1..12067
                     /organism="Mus musculus"
                     /mol_type="mRNA"
                     /strain="mixed"
                     /db_xref="taxon:10090"
                     /chromosome="7"
     gene            1..12067
                     /gene="Dnah3"
                     /note="Derived by automated computational analysis using
                     gene prediction method: Gnomon. Supporting evidence
                     includes similarity to: 1 mRNA, 3 ESTs, 13 Proteins, and
                     99% coverage of the annotated genomic feature by RNAseq
                     alignments"
                     /db_xref="GeneID:381917"
                     /db_xref="MGI:MGI:2683040"
     CDS             499..12021
                     /gene="Dnah3"
                     /note="The sequence of the model RefSeq protein was
                     modified relative to its source genomic sequence to
                     represent the inferred CDS: substituted 1 base at 1
                     genomic stop codon"
                     /codon_start=1
                     /transl_except=(pos:3787..3789,aa:OTHER)
                     /product="LOW QUALITY PROTEIN: dynein heavy chain 3,
                     axonemal isoform X1"
                     /protein_id="XP_006544358.1"
                     /db_xref="GI:568894746"
                     /db_xref="GeneID:381917"
                     /db_xref="MGI:MGI:2683040"
                     /translation="MSDTNCSAQKLDKSDSVHHMSHSQARPELPPLPVSANEEPSELY
                     KTVMSHSFYPPLMQRTSWTLAVPFKEQDHHRGPSDSIGNNYSLTARDMKLKDLLKVYQ
                     PVTISIPRDKTSQGLPLAFTASKTSTEPSKKKMKFNLKAKDDVTGMPFVCKFSSSLSI
                     KNTTDSSVTHPESRPMSPEQQMDVMLQQEMEIESKEQKPSELDLERYYYYLTNGIRKD
                     MIAPENEEVMMRIYKLIPKTLLTTPALEPLQVSLRSEKESDYYYSLMKSIVDYILMDP
                     MEKKRLFIKSIPRLFPHRVIRAPVPWHNIYQSTKKWNEEHLHTVNPMMYKLKELWFAE
                     FQNLRFVRTADLLAGKLPLLPHEYKEVVQKHCREARHILLTKWIPTCAQLFVTQKEHW
                     VHFAPKNDYDSSRNIEEYFASVASFMSLQLRDLVIKSLRDLVSFFMIHKDGNDFKEPY
                     QEMDFFIPQLIMIKLEVRDPIIVFNPTFDDCWQLIKNSFLEIIKNSDGIPKVESILFP
                     DLKGYNMILGTVNPEESLVSDFLDQTLEVFKKNQVGPYKYLNVYKKYDDLLDNMAEKG
                     ISEFLKEKHEIEDFVTSINSIKKRKNEIASMHITVPLAMFCLDAVFLNYDLCERAQNL
                     KDNLILYQVDVNRETNTSICNQYSTIADKVSEIPANTAELVALIEYLKKSSDVTVFKL
                     RRQLRDASERLEFLMDYADLPNEDIKLNSNLFLWPDQIEDVFESSRNLLMSKRDQAEM
                     DLIKRCSEFESRLEGYSKELEMFRKREVMTTEEMKNNVEKLHDLSKNLDLALTEFELI
                     NKEEELLEKEKSSFPLLQTLMINKIPYEQLWVTAYEFSTKSEEWMNGPLYLMNAEQIA
                     EEIGNMWRTTYKLTKTLIDVPAPKRLAENVKLKIEKFKQHIPILNIACNPGMKDRHWQ
                     QISDIVGYEIKPTETTCLANMLEYGFGKFVDKLEPIGAAASKEYSLEKNLEKMKADWV
                     NMCFSFVKYRDTDTSILCAVDDIQLILDDHVIKTQTMCGSVFIKPIEAECRKWEEKLV
                     RVQENLDAWLKCQVTWLYLEPIFSSEDIIAQMPEEGKKFTTVDTYWKSLMAQASEEEC
                     PVTICSLLXNRFFFLSNDELLEILSETKDPLRVQPHLKKCFEGIAKLEFTDNLEIKGM
                     ISSEKETVPFIQTIYPVKAKGMVEKWLQQVEQVMLASMRQVIENGIEAYVQVPRNAWV
                     LEWPGQVVICVSSIFWTREVSEALVEDTLTDFLKKSNDQIAQIVELVRGKLSSGARLT
                     LGALTVIDVHARDVVAKLRHDHINSLNDFQWISQLRYYWTEKNVHVQMITTEALYGYE
                     YLGNSPRLVITPLTDRCYRTLMGALKLNLGGAPEGPAGTGKTETTKDLAKALAKQCVV
                     FNCSDGLDYKAMGKFFKGLAQAGAWACFDEFNRIEVEVLSVVAQQILSIQQAIIRKLK
                     RFIFEGTELSLNPTCAVFITMNPGYAGRAELPDNLKALFRTVAMMVPDYALIGEISLY
                     SMGFLDSRSLAQKIVATYRLCSEQLSSQHHYDYGMRAVKSVLTAAGNLKLKYPEENES
                     VLLLRALLDVNLAKFLAQDVPLFQGIISDLFPGVVLPKPDYEVFLEALNNNIRKMKLQ
                     PVPWFIGKIIQIYEMMLVRHGYMIVGDPMGGKTSAYKVLAAALGDLHAANQMEEFAVE
                     FKIINPKAITMGQLYGCFDAVSHEWTDGVLANAFREQASSITDDRKWIIFDGPVDAVW
                     IENMNTVLDDNKKLCLMSGEIIQMSSKMSLIFEPADLEQASPATVSRCGMIYMEAHQL
                     GWKPLKDSYMDTLPRCLTKEHTELVEDMFTWLVQPCLDFSRLHCKFVVQTSPIHLAFS
                     MMRLYSSLLGKTQCVWLWLPSTDDLNMPAKEVYGAQPPIELLRQWIDHGYWFDKKDTN
                     RLDIVDVLLVTAMGPPGGGRNDITGRFTRHLNIISINAFEDEILTKIFSSIADWHFGK
                     GFDVMFLRYGKMLVQATQTIYRAAVENFLPTPSKSHYVFNLRDFSRVIQGVLLCPHTH
                     LQDLEKFIRLWIHEVYRVFYDRLIDNDDRQTFFNLVKETTSNCFKQTMEKVLIHLSPT
                     GKITDDNIRSLFFGDYLKPESDQKIYDEIIDLRGLTVVMEYYLDEFNSVSKAPMSLVM
                     FKFAIEHISRICRVLKQKKGHLLLVGIGGSGRQSATKLSTFMNSYELYQIEITKNYTN
                     SDWREDLKKIMLQSGVATKSTVFLFSDNQIKHESFVEDINMLLNTGDVPNIFPADEKA
                     DLVEKMQTAARTEGEKVEATPLSMYNFFIERVRKNLHIVLAMSPIGDAFRTRLRMFPS
                     LINCCTIDWFQSWPTDALELVANKFLEDVELDDNIRAEVVSMCKYFQESVKKLSVDYY
                     NTLLRHNYVTPTSYLELILTFKTLLNSKRQEVDTIRNRYLAGLQKLEFASSQVAVMQV
                     ELTALQPQLIQTSEDTAMMMVKIELETKEADAKKLLVQADEKEANAAAAISQAIKNEC
                     EGDLAEAMPALEAALAALDTLNPSDITLVKSMQNPPGPVKLVMESICVMKGLKPERKP
                     DPSGSGKMIEDYWGVSRKILGDLKFLESLKTYDKDNIPSVIMKRIRERFIDHPDFQPA
                     VIKNVSSACEGLCKWVRAMEVYDRVAKVVAPKRERLREAEGKLEIQMQKLNQKRAELK
                     LVEDRLQDLNDEFELMNRKKNSLEKNIEICSQKLVRAEKLISGLGGEKDRWTEAARQL
                     GIRYDNLTGDVLLASGTVAYLGAFTVDYRAQCQNEWLVSCKDKVIPGSVDFSLSNTLG
                     DPIKIRAWQIAGLPVDSFSVDNGIIVSNSRRWPLMIDPQGQANKWVKNMEKTNKLSVI
                     KFSDTNYVRTLENALQFGTPVLLENVGEELDAFIEPILLKATFKQQGVEYMRLGENII
                     EYSREFKFYITTRLRNPHYLPEVAVKVCLLNFMITPLGLQDQLLGIVAAKEKPELEEK
                     KNKLILESAQNKKQLKEIEDKILEVLSLCEGNILEDETAIKILSSSKVLSEEISEKQE
                     IASVTETQIDETRMGYKPVAVHSAAIFFCISDLAHIEPMYQYSLTWFINLYVQSLANS
                     NKSDELDLRIEYIIEHFTLSIYNNVCRSLFEKDKLLFSLLLTVGLLKERKAIDEEVWY
                     FLLTGGVALDNPFPNPAPEWLSEKSWGEIVRASSLQKLKGLMEDVMQNIKVWKDIYDS
                     AWPHEESLPSPWFFLQTLEKIAILRCLRPDKIVPAIQNFICETMGKIFIEAPTFDLQG
                     SYGDSSCCVPLIFILSPGADPMAGLLKFADDVSMGGTKTQTISLGQGQGPIAANMINK
                     AIHEGTWVVLQNCHLATSWMPALEKICEEVIVPENTNSEFRLWLTSYPSEKFPVSILQ
                     NGIKMTNEPPKGLRANLLRSYLNDPISDPLFFQSCTKPVIWQKLLFGLCFFHAIVQER
                     RNYGALGWNIPYEFNESDLRISMRQIQMFLNDYEEVPFEALTYLTGECNYGGRVTDDK
                     DRRLLLSLLSMFYCKEIETDNYHIAPGDAYVIPPYGSYQSYIEYLRTLPITAHPEVFG
                     LHENADITKDNQETNQLFQGVLLTLPRQSGGSGKSPQEVVEELAQDILSKLPNDFNLE
                     EVMKKYPVVYKESMNTVLRQELIRFNRLTKVVRRSLIDLGRAIKGQVLMSSELEEVFN
                     SMLVGKVPAMWAAKSYPSLKPLGGYVADLLARLTFFQEWIDKGPPVVFWISGFYFTQS
                     FLTGVSQNYARKYTIPIDHIGFEFEVTPKETTMENIPEDGAYIKGLFLEGARWDRSTS
                     QIGESLPKILYDPLPIIWLKPGESASFLHQDIYVCPVYKTSARRGILSTTGHSTNYVL
                     SIELPTDMPQKHWINRGVASLCQLDN"
ORIGIN      
        1 ctgggtgcaa gctccattac ccttagaact gatgagggaa gagggagagg aatctagact
       61 cccaggggtg ttatttggga cagaaaatcc tcttttcata ttattatgta actgtgcttc
  ..
    12001 ctgtgccagc tggataattg atggcagggg tcacaaatgc gagagtaaat ggtgtcattc
    12061 tccatga
//
  
=item evigene processing of NCBI genome-genes.cdna

 step1: mouse/evg15mouse/run_tr2aacds.sh
   step2a: parse gene names from genome-genes.cdna
 step2:	mouse/evg15mouse/run_evgmrna2tsa.sh

=item  step1: mouse/evg15mouse/run_tr2aacds.sh <==

  #! /bin/bash
  ## RUN-LOCAL: env trset=arath_TAIR10_20101214up.cdna.gz datad=`pwd` ./tr2aacds_qsub.sh
  ## RUN-CLUSTER: env trset=arath_TAIR10_20101214up.cdna.gz datad=`pwd` qsub -q normal tr2aacds_qsub.sh
  
  ## update 2015, use tr2aacds2.pl improvement, replace orig tr2aacds.pl?
  ncpu=4
  maxmem=4000
  bapps=/bio/bio-grid/mb
  evigene=$bapps/evigene/scripts
  
  #t2ac: app=cd-hit-est, path= echo MISSING_cd-hit-est
  export PATH=$bapps/bin/:$PATH
  #t2ac: app=fastanrdb, path= echo MISSING_fastanrdb
  export fastanrdb=$bapps/bin/fastanrdb
  #t2ac: app=blastn, path= echo MISSING_blastn
  export PATH=$bapps/ncbi229/bin:$PATH
  
  if [ "X" = "X$trset" ]; then echo "Missing env trset=xxxx.tr"; usage; fi
  if [ "X" = "X$datad" ]; then echo "Missing env datad=/path/to/data"; usage; fi
  
  cd $datad/
  echo $evigene/prot/tr2aacds2.pl -tidy -NCPU $ncpu -MAXMEM $maxmem -log -cdna $trset
  $evigene/prot/tr2aacds2.pl -tidy -NCPU $ncpu -MAXMEM $maxmem -log -cdna $trset

=item  step2a: pull gene names from .cdna 

  #! /bin/bash
  ## pullnames.sh arath_TAIR10_20101214up.cdna.gz

  bapps=/bio/bio-grid/mb
  evigene=$bapps/evigene/scripts

  cdna=$1
  nogz=`echo $cdna | sed 's/.gz//;'`
  name=`echo $cdna | sed 's/\.gz//; s/\.cdna//; s/\.rna//; s/\.fasta//; s/\.fa$//;'`
  CAT="gunzip -c"; if [ $cdna = $nogz ]; CAT=cat; fi
  $CAT $cdna | grep '^>' | perl -ne \
  'if(/^>(\S+)/){ $id=$1; chomp; s/>//; @v=split /\s*\|\s*/; ($id2,$sym,$na,$loc)=@v; 
  $na=~s/;\s.*//; $na="noname" unless($na=~/\w\w/); print join("\t",$id,$na)."\n"; }  ' \
    > $name.names.orig;  
  ##need nameclean for ncbi checker
  $evigene/scripts/namecleangb.pl $name.names.orig > $name.names

=item  step2: mouse/evg15mouse/run_evgmrna2tsa.sh <==

  #! /bin/bash
  ### env idprefix=MysppEGm trclass=myspp.trclass datad=`pwd` ./run_evgmrna2tsa.sh
  # $evigene/scripts/evgmrna2tsa2.pl -species Daphnia_magna -idprefix Dapma5xEVm -novectrim -class $pt.trclass -log
  # add  -names $trset.names;   arath_tair10.cdna test data has names in hdr;
  
  if [ "X" = "X$species" ]; then species=Arabidopsis_thaliana; fi
  if [ "X" = "X$idprefix" ]; then idprefix=Arath10EVm; fi
  
  ncpu=2; # most 100K trsets run on 8cpu in 10m-20min; 300K set can take 1hr.
  ## common option: skip vecscreen and tbl2asn == -novectrim (not def)  -noruntbl2asn (default)
  
  ## local path to bio-apps
  biopath=/bio/bio-grid/mb
  export EVIGENES=$biopath/evigene/scripts
  
  ## NCBIc++ version: /bio/bio-grid/mb/ncbi229/bin/vecscreen
  ## NCBIcc version: /bio/bio-grid/mb/ncbic/bin/vecscreen # both now work w/ evigene ..
  export vecscreen=$biopath/ncbi229/bin/vecscreen
  export tbl2asn=$biopath/ncbic/bin/tbl2asn
  
  opts="-debug  -dropshow -skipdropseq -NCPU $ncpu"
  
  if [ "X" = "X$datad" ]; then echo "missing env datad=path/to/data"; exit -1; fi
  if [ "X" = "X$trclass" ]; then "echo env trclass=path/to/name.trclass"; exit -1; fi
  ## .. these are now read via sra_result.csv, species => idprefix
  if [ "X" != "X$idprefix" ]; then opts="$opts -idprefix $idprefix"; fi
  if [ "X" = "X$vectrim" ]; then opts="$opts -novectrim"; fi
  if [ "X" != "X$tbl2asn" ]; then opts="$opts -runtbl2asn"; fi
  if [ "X" != "X$species" ]; then spp=`echo $species | sed 's/ /_/g;'`; opts="$opts -species=$spp"; fi
  #Fixme# if [ "X" != "X$sra" ]; then opts="$opts -sraids=$sra"; fi
  
  if [ "X" = "X$names" ]; then names=`echo $trclass | sed 's/.gz//; s/\.trclass/.names/;'`; fi
  if [ -s $names ]; then opts="$opts -names $names"; fi
  
  cd $datad/
  ##FIXME: template files; need these or generate defaults? # if [ ! -f evgr_tsasubmit.sbt ]; then
  
  echo $EVIGENES/evgmrna2tsa.pl  $opts -log -class $trclass
  if [ ! -f $trclass ]; then echo "ERR: missing -class $trclass"; exit -1; fi
  $EVIGENES/evgmrna2tsa.pl  $opts -log -class $trclass


=item follow on scripts for Evigene tr2aacds, evgmrna2tbl outputs vs NCBI rna inputs
  
  pt=mouse_ncbi1503
  
  gzgrep '^>' tmpfiles/${pt}nr.cds.gz | perl -ne \
  's/>//; @d=grep{ not/utrorf/} split; print join("\t",@d)."\n" if(@d>1);' \
   > tmpfiles/${pt}nr.cds.dupids
    
  gzgrep '^>' $pt.cdna.gz | \
  env spp=mouse MID=MGI perl -ne 'BEGIN{ $MID=$ENV{MID}||"nada"; $SPP=$ENV{spp}||"nospp"; } \
  ($d)=m/>(\S+)/; ($t,$w,$loc,$alt)=m/(?:type|len|geneid|isoform)=([^;\n]+)/g; print join("\t",$d,$w,$t,"LOC$loc",$alt); s/MGI:MGI/MGI/;  \ 
  for $k (qw(gene ncbiprot $MID Name)) { ($v)=m/$k[:=]([^=:;,\n]+)/; $v||="na"; print "\t$k:$v"; } print "\n"; ' \
   > $pt.attr.tab
  
  env spp=mouse \
  perl -ne 'BEGIN{ $MID=$ENV{MID}||"nada"; $SPP=$ENV{spp}||"nospp"; } \
  next if(/^\W/); chomp; @v=split"\t"; if($v[1]=~/^(main|noclass|NOMAIN)/) { ($md,$mc,$ad)=split; @ad=grep/\w\w/,split",",$ad; \
  ($mda)=grep { m,$md/,} @ad; unshift @ad,"$md/$mc" unless($mda);  \
  for $dc (@ad) { ($d,$c)=split"/",$dc;  next if($d=~/utrorf/ and $c=~/drop/); \
  $dup=$dupid{$d}||"nod"; $at=$tda{$d}||"na";  print join("\t",$d,$c,$md,$at,$dup)."\n" unless($did{$d}++); } } \
  elsif($v[2] =~ m/RNA/) { ($td,$tw,$ty,$gid,$alt,$gna,$pid)=@v; $tda{$td}="$ty,$tw,$gid.t$alt,$pid"; } \
  elsif($v[1] =~ m/^$SPP/) { $td=shift @v; $dupid{$td}=join",",@v; } \
  BEGIN{ print join("\t","#TranscriptID","TrClass","MainID","Attr","DupCDSid")."\n"; } ' \
   tmpfiles/${pt}nr.cds.dupids $pt.attr.tab publicset/$pt.mainalt.tab > publicset/$pt.allin.tab
  
  #................... locus eq .............
  
  env spp=mouse PUBID=Mouse perl -ne 'BEGIN{ $PUBID=$ENV{PUBID}||"nada"; $SPP=$ENV{spp}||"nospp"; } \
  next if(/^\W/); chomp; @v=split"\t"; \
  if(/^$PUBID/) {($pd,$od,$gd,$ti,$cla,$aw,$pia)=@v; $pod{$od}="$gd\t$ti"; $aw{$od}=$aw; } \
  elsif(/^\w/) { ($od,$cla,$md,$at,$dupd)=@v; $pd=$pod{$od}||"nopd\t0";\
  $aw=$aw{$od}||"na";  @at=split",",$at; ($lc)=grep/^LOC/,@at;\
  ($mt,$mw)=@at; $lc=~s/\.t/\t/; $od=~s/^$SPP://; \
  print join("\t",$od,$cla,$pd,$lc,$mt,"$mw,$aw")."\n"; } ' \
    $pt.pubids $pt.allin.tab > $pt.loctab
  
  env spp=mouse PUBID=Mouse perl -ne 'BEGIN{ $PUBID=$ENV{PUBID}||"nada"; $SPP=$ENV{spp}||"nospp"; } \
  chomp; ($od,$cl,$pd,$ti,$lc,$li,$ty,$maa)=@v=split"\t"; \
  $lc||="nolc"; $lcn{$lc}++; $pdn{$pd}++;  \
  if($cl=~/drop/) { $ldrop{$lc}++; }  \
  else { $pdl{$pd}{$lc}++; $lpd{$lc}{$pd}++; } $odv{$od}=[@v];  \
  END{ @pd=sort keys %pdl; @lc=sort keys %lpd; for $lc (@lc) { @p=sort keys %{$lpd{$lc}};  \
  $c=$lcn{$lc}; $d=$ldrop{$lc}||0; print "$lc:$c";  \
  for $p (@p) { $c=$lpd{$lc}{$p}; print "\t$p:$c";  \
  @np=sort keys %{$pdl{$p}}; print ",".scalar(@np) if(@np>1); }  \
  print "\tdrop:$d"; if(@np>1) { $np=join",",@np; print "\tdup:$np"; } \
  print "\n"; } } ' \
    $pt.loctab > $pt.loceq.tab

=item mouse_ncbi1503 cdna evigene class table

  * most mRNA type are either retained (locus reps = main, noclass),
   or dropped as high-identity, redundant alts of retained transcripts

  * $pt.loceq.tab tabulates agreement b/n NCBI LOC ids and Evigene locus IDs
    for model organisms mouse and zebrafish.  Disagreement is approx 
      5% NCBI paralogs (separate loci) are mis-classed as Evigene alternates (one locus)
      1% NCBI alternates (one locus) are mis-classed as Evigene paralogs (separate loci)
    This agrees with earlier evigene tr2aacds development testing and parameter adjustment
    with other species gene sets.

        mouse: 1265/22388, 5.65% LOC par>alt misclass as 1 Evigene locus
                215/22388, 0.96% LOC alt>par misclass as 2+ Evigene loci
            
        zfish: 1324/24820, 5.33% LOC par>alt misclass as 1 Evigene locus
                137/24820, 0.55% LOC alt>par misclass as 2+ Evigene loci
            

  pt=mouse_ncbi1503
  grep mRNA  publicset/$pt.allin.tab | cut -f2 | sort | uniq -c | sort -k1,1nr | head -20
  16350 dropalthi1
  11050 dropperfectdup
  10555 main
  9190 noclass
  5833 dropperfectfrag
  5184 althi1
  4115 dropalthi
  2877 althi
  1013 dropnoclass
   756 dropparthi
   399 dropparthi1
   374 altmid
   356 dropmain
    39 altmidfrag
    23 dropaltmid
     9 dropaltmidfrag
     1 NOMAIN
  
    * most noncoding RNA type are dropped; tiny or missing ORF proteins lead to
      many "perfect duplicates" of coding sequences that are dropped.
  
  grep -v mRNA  publicset/$pt.allin.tab | cut -f2 | sort | uniq -c | sort -k1,1nr | head -20
  17531 dropperfectdup
  13935 dropnoclass
  3258 dropalthi1
  2511 dropalthi
  2353 dropmain
  1429 noclass
   817 main
   770 dropparthi
   670 dropperfectfrag
   515 althi
   502 althi1
   155 dropaltmid
   150 dropparthi1
    83 altmid
    56 dropaltmidfrag
    26 altmidfrag
     4 NOMAIN

=item original: pull proteins from cds of gene records in genbank format

     gene            <62188..>75177
                     /gene="EYA"
                     /locus_tag="DAPPUDRAFT_204955"
     mRNA            join(62188..62346,62647..62736,67111..67169,67365..67525,
                     67667..67781,69784..69905,70068..70168,73875..75177)
                     /gene="EYA"
                     /locus_tag="DAPPUDRAFT_204955"
                     /product="eye-absent 1-like protein"
                     /note="Has EST support. PASA Supported: pasa:asmbl_1"
     CDS             join(62314..62346,62647..62736,67111..67169,67365..67525,
                     67667..67781,69784..69905,70068..70168,73875..73955)
                     /gene="EYA"
                     /locus_tag="DAPPUDRAFT_204955"
                     /note="member of eyes absent gene family, confirmed by
                     maximum likelihood tree analysis implemented in PHYML;
                     GO_process: GO:7275 - development"
                     /codon_start=1
                     /product="eye-absent 1-like protein"
                     /protein_id="EFX89861.1"
                     /db_xref="GI:321478905"
                     /db_xref="JGIDB:Dappu1_204955"
                     /translation="MTGSYASRFGKDVPGSIQLGVGMEELIFNLSDTHFFFNDLEECD
                     QVHIDDVSSDDNGQDLSNYNFRTDGFHAAATNASLCLATGVRGGVDWMRKLAFRYRAI
                     KEIYNRYRNSVGGLLAPAKREQWIQLRMEIESLTDNWLTLVTKCLELINSRPNAVNVL
                     VTTTQLVPALAKVALYGLGGVFAIENIYSATKIGKESCFERIVSRFGRKCTYVVVGDG
                     RDEESAAKQLNFPFWRIASHQDSAALYNALDLGYM"


=cut
