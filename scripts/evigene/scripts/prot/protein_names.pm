# protein_names.pm

# package protein_names;
package main; # include symbols in main w/o Exporter/ISA.. hassles

=item about

  my $pi= ( $name =~ m/\((\d+)%.\)/ ) ? $1 : $MIN_PROTIDENT; 
  my $locusadd= ($NAME_UNKADDLOCUS) ? "locus $gid" : "";

  my($newna,$lowqualname)= nameclean( $name, $pi, $locusadd );

=cut

use strict;
use warnings;

use vars qw (
 $NAME_NONE $NAME_UNK $NAME_NOLIKE $NAME_NOPUTATIVE $NAME_NODIGITS $NAME_UNKADDLOCUS  $NAME_IDPATT $NAME_KEEPID $NAME_UCFIRST
 $MIN_NAMEIDENT $MIN_IDLIKE  $MIN_CERTAIN $NAMED_MINALIGN $NAMED_ISCLEAN $SPPinNAME
 @NAMETAG_CUT @NAMETAG_KEEP
 $NAME_GOODSPP $NAME_POORSPP
 @TEnames  $USE_TENAME
);

use constant { NAMEDIFF_MINOR => 1, NAMEDIFF_MAJOR => 2, };
use constant MAXNAMELEN => 70;
use constant CDD_MAXNAMELEN =>  MAXNAMELEN - 10; # was 39;

our $NAME_NONE = "Unknown|Uncharacterized conserved|Uncharacterized|Hypothetical|noname";  
our $NAME_UNKUNIP= 'Uncharacterized protein'; # Uniprot std
our $NAME_UNKNCBI= 'hypothetical protein';  # NCBI prefers this to Uncharacterized .. special handling
our $NAME_UNK  = $NAME_UNKUNIP;
our $NAME_UNKCDD = 'Domain of unknown function'; # DUF nnn stuf

our $NAME_UNKADDLOCUS= 1; # policy, add/not the gene id to end of UNK names
our $MIN_NAMEIDENT = 35;  # min similar% for protein evid align/equivalence, for naming; JCVI uses 35%  
  # -- note MIN_ID is used for both naming and keeping prot evidence, split this? use lower IDLIKE for evid?
  # MIN_IDENTITY == MIN_NAMEIDENT now
our $MIN_IDLIKE   = 15;  # low for now; 15..20 seems right;add Note with loqualname
our $MIN_CERTAIN  = 66;  # putative naming, what?

our $NAME_AddBackProtprefix=1; # FIXME.1510 need opt  
our $NAME_NOUSCORE= 1; # FIXME.1510 ncbi name opt

our $NAMED_MINALIGN = 0.60; # or 60% or $MIN_CERTAIN ??

## ^CG\d/ || AGAP .. allow for -PA,-PB, other isoform suffix?   ACYPI009998
## zfish:  Zgc:0000 Zgc:112221 ENSP00000382042
## human: cDNA FLJ30174 fis, clone BRACE2000975, terms-here...
## human: C6orf10 ? becomes Uncharacterized C6orf10 protein; bad;  DKFZp434B061 ?  DKFZp\d+\w+
## Macaca fascicularis brain cDNA clone: QtrA-18802, similar to .. UGGGGGHHHHHHHH
## Testis cDNA clone: QtsA-13078, similar to ...
## arabid: Domain of unknown function (DUF1767) << keep DUFid? change to Unchar? Unchar conserved? NCBI sez: "protein of unknown function"

our $NAME_IDPATT= 
      '(?:At|at|AT|Os)\d{1,2}[Gg]\d{3,}|DUF\d\d+|OSJNB\w\d+[\w.]+|' # plants 
      .'(:?clone|[cC]DNA) \w+\d\d+[,]?|C\d+orf\d+\w?|DKFZp\d+\w+|(?:ENS\w+\d|Zgc:|DDB_G|LOC|FLJ|KIAA|LINC)\d\d+|' # human/vert, may have several..
      .'(?:AGAP|AAEL|ACYPI|Phum_PHUM|PF11_|BcDNA.[CG][A-Z]|CG|G[A-Z]|GLEAN_|NM_|NP_|XP_)\d\d\d+'; # bugs;
  ## add: DDB_G0283697 ;                
  ## special case mess: 'HERV\-\w\S+|' #HERV-K_5q33.3 provirus ancestral Pol protein  
  ## weird NCBI-ID in name: KIF5B-RET(NM_020975)_K16R12 fusion protein
  ## .. cut all of name/id: '\w[\w-]+\WNM_\d+\W\w+' or just ncbi name/id ? (NM_020975)  '\(NM_\d+\)'

  # 1510: more yuck to chop:  at5g43400/mwf20_9 < case change from At5g..; Phum_PHUM454910; PF11_0207; Peptide XP_001950702 similar hypothetical protein;
    
  ## cut SPPinName, need species dictionary now
  ## * put dual names 1st for cut precedence : Escherichia coli|Escherichia
our $SPPinNAME='H.\s*sapiens|Homo sapiens|Mus musculus|human and mouse|human|mouse|rat|Drosophila|melanogaster'
            . '|C. elegans|yeast|S. cerevisiae|cerevisiae|S. pombe|Schizosaccharomyces pombe'
            . '|Escherichia coli|E. coli|Escherichia|ECOLI|mycoplasma|pylori'
            . '|pseudomonas|Staphylococcal|streptococcal|Streptococcus|streptomyces'
            . '|Aspergillus|Bacteroides|Helicobacter|Mycobacterium|Plasmodium|Tuberculosis'
            . '|Chlamydomonas|Arabidopsis|thaliana|xenopus|avian';
  
  # s/\b($SPPinNAME)\b\s*//ig;  # case insense ok?
  
our $NAME_KEEPID= 1; # yes or no? caller choice FIXME: 2 meanings 
# FIXME: NAME_KEEPID has 2 meanings: a. Name == ID instead of Unknown; b. ID of protein
our $NAMED_ISCLEAN= 0; #? call nameclean() or not for blast,.. need opt.
our $NAME_NOLIKE=0;
our $NAME_NOPUTATIVE=0;
our $NAME_NODIGITS=0;
our $NAME_UCFIRST= 1; # UniProt std
our $USE_TENAME= 0; # $NAME_USETELIST
our @TEnames; # @NAME_TELIST;

our  @NAMETAG_CUT= qw(sym n Tax len Gene); # is Gene= from swissprot keep or cut?
our  @NAMETAG_KEEP= qw(RepID); # qw(sym N Tax RepID DbXref); # ignorecase?

our($NAME_GOODSPP,$NAME_POORSPP)= ('HUMAN|MOUSE|ARATH|CAEEL|DROME','BADSPECIES'); # consensusname options ..


=item nameclean

  #NOTE rename Thecc1EG026905t1: 'Uncharacterized protein locus Thecc1EG026905' < 'Prefoldin chaperone subunit family protein (19%T)'
  #NOTE rename Thecc1EG026905t2: 'Uncharacterized protein locus Thecc1EG026905 isoform 2' < 'Prefoldin chaperone subunit family protein (24%T)'
    # -- should these be putative instead of Unk? or -like? or family?
    # Proteins of unknown function which exhibit significant sequence similarity
    # to a defined protein family have been named in accordance with other members
    # of that family... e.g. "Holliday junction resolvase family endonuclease".
    # It is also possible to use "-like" in the name, for outliers

See update in evigene/scripts/namecleangb.pl
   -- use for genes.gff naming before here;
   -- fix below to use above as package sub?
   
See also longer method in evigene/scripts/bestgenes_puban_cacao.pl : nameclean()

=item more tbl2asn name fixes

  dependant  initation  glcosylase  monooxigenase phopholipase  cholinephosphotranferase # bad spelling
  dependent  initiation glycosylase monooxygenase phospholipase cholinephosphotransferase # corrected
  fibre sulphate sulphoglucosamine sulphohydrolase   # britglish
  fiber sulfate  sulf..     sulf..  # amerglish
  
cat $namefix | grep '^Changed' | sed 's/ for CDS.*//;' | sort -u
 'CAMP-dependant protein kinase' to 'CAMP-dependent protein kinase'
 'CDD: SUI1_rel, Translation initation factor SUI1-like' to 'CDD: SUI1_rel, Translation initiation factor SUI1-like'
 'CDD: UDG_F2_MUG, G:T/U mismatch specific DNA glcosylase fragment' to 'CDD: UDG_F2_MUG, G:T/U mismatch specific DNA glycosylase fragment'
 'CDD: sulP, High affinity sulphate transporter 1-like fragment' to 'CDD: sulP, High affinity sulfate transporter 1-like fragment'
 'Cytochrome P450 monooxigenase CYP4Q3, putative fragment' to 'Cytochrome P450 monooxygenase CYP4Q3, putative fragment'
 'Cytochrome P450 monooxigenase, putative' to 'Cytochrome P450 monooxygenase, putative'
 'Giant fibre A fragment' to 'Giant fiber A fragment'
 'Giant fibre A' to 'Giant fiber A'
 'N-sulphoglucosamine sulphohydrolase fragment' to 'N-sulfoglucosamine sulfohydrolase fragment'
 'N-sulphoglucosamine sulphohydrolase' to 'N-sulfoglucosamine sulfohydrolase'
 'Phopholipase D-like fragment' to 'phospholipase D-like fragment'
 'sn-1,2-diacylglycerol ethanolamine-and cholinephosphotranferase, putative fragment' to 'sn-1,2-diacylglycerol ethanolamine-and cholinephosphotransferase, putative fragment'
 'sn-1,2-diacylglycerol ethanolamine-and cholinephosphotranferase-like fragment' to 'sn-1,2-diacylglycerol ethanolamine-and cholinephosphotransferase-like fragment'
   
=cut

sub nameclean
{
  my ($namin, $pi, $locusname)= @_;
  my $lowqualname=undef;
  my ($isunk,$isput,$islike,$iste)= (0) x 9;
  # $locusname="" unless(defined $locusname);
  
  #*** ?? IS THIS FAILING .. missing haem, Tumour
  local $_= $namin; 
  # my $_ = $namin;   # perl 5.9.1

  ## s/\s+\(\d+%.*\)//; # trailing pctident
  my $namepi= (s/\s+\((\d+)%.*\)//)?$1:0;
  unless(defined $pi and $pi =~ /\d/) { $pi= $namepi || $MIN_CERTAIN; }

  ## FIXME: add NAME_IgnoreParts = '(sym|Tax|RepID)[:=]'; 
  ## option to remove/reattach parts to name/description, preserve from cleaning...
  ## these are all extra fields of format 'key[=:]value[.;]'
  ## fixme: tag=words with spaces;delim

  my @nametags=();
  foreach my $tag (@NAMETAG_KEEP) {
    if(/\b$tag[:=]/) { my $val='';
    if(s/\s*$tag[:=]\s*([^;,=]+)[;,]//) { $val=$1; }
    elsif(s/\s*$tag[:=]\s*(\S+)//) { $val=$1; }
    push @nametags, "$tag:$val"; 
    }
  }
  foreach my $tag (@NAMETAG_CUT) {
    if(/\b$tag[:=]/) { my $val='';
    if(s/\s*$tag[:=]\s*([^;,=]+)[;,]//) { $val=$1; }
    elsif(s/\s*$tag[:=]\s*(\S+)//) { $val=$1; }
    }
  }

  my $nameprechange= $_;  # but after cut tags; before CDD cut?
  my $nameidpatt= $NAME_IDPATT; ## $public_options{nameidpatt} || $NAME_IDPATT;
  my $isid=0;
  if(s/\-like\b//) { $islike++; } # remove like early, messes other parsing

  ## FIXME: ^CDD: needs to be cut/added back to front of name as CDD key .. old stuff
  my $iscdd=(s/^CDD:\s*//)?1:0;
  my $cddsym='';
  if($iscdd) { # CDD: symbol, description.. temp cut symbol
    if(s/^(\S+),\s+(\w)/$2/) { $cddsym=$1; $cddsym='' if($cddsym=~/domain/i); } # like $isid, not same; stutter here? Domain of .. Domain ; not here,
    }
   
  #  horribles die quickly...
  if(/The protein encoded by this gene was identified/i) { $_=$NAME_UNK; $isunk=1; }
  elsif(/^(Protein|Domain|Family) of unknown function\W+((?:DUF|UPF)\d+)/i) { $isid=$1; $isunk=1; $_=($iscdd)?$NAME_UNKCDD:$NAME_UNK; } 
  elsif(/^(Protein|Domain|Family) of unknown function/i and $cddsym) { $isid=$cddsym; $isunk=1; $_=($iscdd)?$NAME_UNKCDD:$NAME_UNK; }  
  #o#elsif(/([A-Za-z]\w+\d+), (Protein|Domain|Family) of unknown function/i) { $isid=$1; $isunk=1; $_=($iscdd)?$NAME_UNKCDD:$NAME_UNK; }  
  #o#elsif($iscdd and /([A-Za-z]\w+\d+), .* of unknown function/i) { $isid=$1; $isunk=1; $_=$NAME_UNKCDD; }  
  elsif($iscdd and /unknown function\W+((?:DUF|UPF)\d+)/i) { $isid=$1; $isunk=1; $_=$NAME_UNKCDD; }  
  ##elsif($iscdd and /of unknown function\W+((?:DUF|UPF)\d+)/i) { $isid=$1; $isunk=1; $_=$NAME_UNKCDD; }  
  ## CDD: Unknown function DUF3886 << isid/stutter bug, ingore 'of '
  elsif($iscdd and /of unknown function/i and $cddsym) { $isid=$cddsym; $isunk=1; $_=$NAME_UNKCDD; }  
  # ?? preserve DUFnnn if $iscdd ?? 
  #  Change to 'DUFnnn, domain of unknown function' for $iscdd ?
  #  Change to 'Unknown function DUFnnn domain-containing protein' for prot product name ?
  # also: Family of unknown function ; Plant protein of unknown function; 
  
  ## not this 'CBS domain-containing protein with a domain of unknown function (DUF21)'
  elsif(/Encodes a close of the Cauliflower OR/i) { $_='Orange protein'; }
  ## Dang this monster got thru again: TAIR:AT5G61670.2  Encodes a close of the Cauliflower OR (Orange) protein,..
  ## okay in updated arath_TAIR10.aa.nameclean; problem was uncleaned names
  
  ## FIXME: turn into BEGIN hash list of translations ; need spelling checker :(
  # typos and britglish > amglish; need config list
  ## tetramerisation utilisation polarisation oligomerisation
  # s/Uncharacterise/uncharacterize/ig; s/\brecognise/recognize/ig;
  s/(Uncharacter|Recogn)ise/${1}ize/ig; # both isa,ise ??  
  s/(Organ|Util|Local|Polar|Dimer|Oligomer|Tetramer)isa/${1}iza/ig; # both isa,ise ?? ..zation; suffixes assumed
  s/fibre\b/fiber/ig; s/monoxygenase/monooxygenase/ig; # [Mm]onox..
  s/signalling/signaling/ig; # Two Ls in British English, one in American English. 
  s/tumour/tumor/ig;   
  if(/haem/i) { s/haemo/hemo/ig; s/\bhaem/heme/ig; s/haem\b/heme/ig;} # Dihaem>Diheme? haemoprotein Quinohaemoprotein>Quinohemoprotein
  if(/sulph/i) {
    #  sulphur sulphide sulphate sulphoglucosamine sulphohydrolase   # britglish
    s/\bsulph(\w+)/sulf$1/ig; 
    s/\b(tri|di|bi)sulph(\w+)/$1sulf$2/ig; # FIXME:Disulphide ; other xxxsulph di- bi- mono- tri- ?
    }
    
  # Typos/bad spelling (Uniref/Uniprot names)
  ## 'Ptotein (Zgc:112221)' << ugh, new zfish typo
  s/characteris\b/characters/ig;  s/\bcontaing/containing/ig;  s/\bcomponenet/\bcomponent/ig;
  s/dependant/dependent/ig;  s/initation/initiation/ig;
  s/monooxigenase/monooxygenase/ig; s/phopholipase/phospholipase/ig;  
  s/glcosylase/glycosylase/ig; s/aluminium/aluminum/ig;  
  s/cholinephosphotranferase/cholinephosphotransferase/ig;
  s/(proteine|ptotein|\bprotei\b)/protein/ig;  
  
#  # UPDATE: add to species.config some of these nameclean patterns
#   my $nameclean = ($config{nameclean} and ref($config{nameclean})) ? $config{nameclean} : {};
#   foreach my $nk (sort keys %$nameclean) {
#     my $npatt= $nameclean->{$nk};
#     # if($npatt =~ m/^s(.)/) { my $nc=$1; my($ncut,$nto)= split/[$nc]/,$npatt; s/$ncut/$nto/; }
#     my $eok= eval $npatt; # is this ok? prefer s/$nacut/$nato/ ; 
#     # any use for m/$napatt/ ; need specific keys for this, m/$ISGENEID/ =  m/^(Os|At|AT)\d{1,2}g\d{3,}/
#   }
  
  # FIXME here: XXXX-garbage similar to (NAME HERE) needs special cut..
  # FIXME2 : move above all m/^patt/ cases ?
  # .. mRNA at end is common indicating name from mRNA transcript..
  # fixme3: leading garbage maybe problem: cDNA FLJ55953, highly similar to ADP-dependent glucokinase (EC 2.7.1.147) > ADP.., highly
  #o# if(s/^(.*)(?:(weakly|moderately|highly) similar|Similarity|related) to\s*//i) 
  if(s/^(.*)(weak|moderate|high)(ly | )(similar|similarity|related) to\s*//i) 
  { 
    my $garbagemaybe=$1; 
    my $isgarbage= ($garbagemaybe =~ /\b($nameidpatt)\b/)?1:0;
    if($isgarbage) { $isid=$garbagemaybe; } else { s/^/$garbagemaybe /; }
    $islike++; ## set islike=1 ?
  } elsif(s/\b(similar|Similarity|related) to\s*//i) {
    $islike++; 
    s/^(cDNA|mRNA) clone with//i; # mRNA clone with similarity to Real Thing Here ..
  }
  
  if(s/^TE://i) {  $iste=1; }  #?? set this here or not
  ## FIXME: add TE check option: $iste= $USE_TENAME && isTEname($_); 
  s,/sw$,,; s,/ARP[\w]+$,,; # my garbage: /sw = swissprot, /ARPnnn arpod group id; bug:'.. subunit/sw-like protein'
  
  s/\bLOW QUALITY PROTEIN:\s*//i; 
  s/^[Ii]nvolved in\b\s*//; # [Cc]ontains handled below..
  ## not all# s/\s*\b[Ii]nvolved in\b\s*/ /; # CDD: Contains xxx; Involved in ...
  ## bad prior CDD clean of 'AidB, Proteins involved in DNA damage response, similar to the AidB gene product'  
  ##  >> CDD: AidB, s involved in DNA damage response

  s/\bPredicted\b\s*/ /i;
  if(s/^(Conserved|Expressed|Novel|Alternative protein)\s+//i) { }  # maybe set $isunk? # s/PREDICTED[:]?\s*//i;
  #no# if(s/\b(Conserved|Novel)\b\s+/ /i) { }    
  ## bad? Conserved in .. valid Conserved: Evolutionarily conserved signaling..
  ## CDD: conserved common.. some with DUFnnn; conserved region; conserved HD motif; ..
  #no# if(s/^(Expressed|Alternative protein)\s+//i) { }  # maybe set $isunk? # s/PREDICTED[:]?\s*//i;
  if(s/^(Possible|Candidate|Putative|Potential|Probable)\s+//i) { $isput=1; }  
  # if(s/hypothetical\s+//) { $isput=1; } # defer after NAME_NONE
  
  ## domain syntax changes: ncbi wants XXX domain-containing, change all '^CDD:' to that here?
  #others: Containing a NACHT domain; Containing a small cytokines (Intecrine/chemokine), interleukin-8-like domain
  #  : Containing a membrane-associating domain
  if(/^Contain/i and /\bdomain/i){ s/^Contain\w+\s*//i; s/^(an|a|the)\s+//; s/\bdomain[s]?/domain-containing/i; }
  else { s/^[Cc]ontai(ns|ning)\bs*//; }
  
  #  # various crap, parts below: Full.length cDNA|Partial cDNA sequence
  #  ## Full-length cDNA;  Full-length cDNA clone CS0DI058YO07 of
  if(/^(Full.length cDNA clone|Full.length cDNA 5.PRIME end of clone)/i) {  $_= $NAME_UNK; $isunk=1;  }
  ## all 'cDNA clone' names starting w/ Macaca or Testis ..  ?? Testis cDNA, clone: QtsA-16864
  elsif(/^(Macaca|Testis).*cDNA clone/i) {  $_= $NAME_UNK; $isunk=1;  }
  ## uniprot vert junk names: "Chromosome undetermined SCAF10187", Chromosome xxx
  ## Chromosome 18 open reading frame; Chromosome 21 SCAF25229
  elsif(/^Chromosome \S+ open reading frame/i ) {  $_= $NAME_UNK; $isunk=1;  }
  elsif(/^[Cc]hromosome .* (reading frame|undetermined|SCAF\d+)/ ) {  $_= $NAME_UNK; $isunk=1;  }
  elsif(/^si:(ch|d|rp|zfos)|^si:\S+ novel|^Zgc:\d/i){ $_= $NAME_UNK; $isunk=1;  } ## zebrafish crap names: Si:ch211-229l10.1; si:dkey-15h8.5
  elsif(/^DNA\s+for\s*/i){ $_= $NAME_UNK; $isunk=1;  } ##  #?? DNA for LINE-1 transposable element ORFI

  s/^HERV\-\w\S+\s+//; ## HERV-T_19q13.11 provirus ancestral Env polyprotein
  s/^([Cc]DNA|[Mm]RNA)[,;:\.]\s+//; 
  s/,\s*(cDNA|mRNA)\s*$//;  # human garbage; but not 'mRNA-binding' ...
  ## what is 'BcDNA' ; ?? Alternative protein C22orf29
  
  #o# s/\(\s*includes:[^\)]+\)//i;  # problem missing ) at end
  s/\(\s*includes:[^\)]+[\)]?//i;   
  
  s/\s*[^\w\s]?(fragment|partial)[^\w\s]?//i; #ok? #was# s/\s\(Fragment[s]?\)//; 
  s/^Isoform \S+ of\s*//; 
  s/\s+[Ii]soform\s*.*$//; ## ^Isoform 2 of; Isoform VEGF-B167 of ..
  s/[,]?\s*transcript variant \w+//i;
  
  # * FIXME: some of UniP valid names are "Protein XXX", eg. Protein NEDD1 /NEDD1_HUMAN; Protein SERAC1 / SRAC1_HUMAN; 
  # .. any words following Protein are ok?? 'Protein tyrosine phosphatase' ..
  # .. need to distinguish valid UniP from other junk names?
  my $isnamedprot=0; ## treat like $isid ??  no, keep as-is
  # if(/^Protein ([A-Za-z][A-Za-z]\S*)/) { $isnamedprot=$1; }
  if(/^[Pp]rotein ([A-Za-z]\S*)/) { $isnamedprot=$1; } # should be [Pp]rotein allow lc
  elsif(not m/protein\-protein/) { s/^(ORF|protein)\s*//i; } # some cuts restore leading Protein
    ## highly similar to Homo sapiens protein xxxx >> becomes Protein xxx, should be Xxx
    ## .. keep leading Protein on these also ??? ie almost anything following /^Protein / ???
    # C10_HUMAN       C10     Q99622  C12orf57                update  old:Protein C10
    # ETS1_HUMAN      C-ets-1 P14921  ETS1            update  old:Protein C-ets-1
   
  #old# s/protein family/family protein/i; # to uniprot syntax
  s/\bprotein family/family protein/i; # to uniprot syntax
   ##FIXME: bad for 'glycoprotein family protein' and similar
  s/\bgene protein\b/protein/; # seems redundant..
 
  # no-no species names: Arabidopsis yeast  human
  # >> staphylococcal? = Staphylococcal nuclease ue, TAIR name
  # and: genome complete  pseudogene? = TAIR name
  # my $namedrops= $public_options{namedrops} || 'Arabidopsis|thaliana|yeast|complete sequence|complete|genome|pseudogene'; #plants;
  # s/\b($namedrops)[,\.\s]*//ig;
  
  s/Whole genome shotgun.*//i;  # W g s (sequence|assembly)..
  ## eg Whole genome shotgun assembly, reference scaffold set, scaffold scaffold_1
  s/\s*\Wgene\Wpseudogene\W//i; ## ugh: prosaposin 1 (gene/pseudogene) .. namer cant make up mind.
  if(s/\b($SPPinNAME)\b\s*//ig) { s/^[Pp]rotein\s+//; }
  s/\b(complete sequence|Partial \w+ sequence|complete|genome|pseudogene|Gene supported by)[,\.\s]*//ig;
  s/\s*\((?:InterPro|TAIR):[\w\.-]+\)//ig; #  (TAIR:AT1G22000.1); (InterPro:IPR010678) # << move to IDPATT? or DbxrefPATT ?
  if(s/paralog of //i) { $islike++; }
  #above# if(s/\-like//) { $islike++; } # remove like above, messes other parsing
  
  ## FIXME: NCBI tbl2asn doesnt like EC nums now: EC 1.14.14.1 .. save it?  
  ## (EC 6.3.2.-) ; (EC 3.1.-.-) ;  also got some 'EC1.2.3.'
  ## Warn SEQ_FEAT.EcNumberProblem  Apparent EC number in protein title ..
  ## 1505: EC=nnnn also, EC=1.-.-.-/sw also; Ugh.. one as TWO EC=1,EC=2 "Dual specificity..", add s///g ??
  #o if(/\bEC\s*(\d+\.(?:[\d\.\-]*))/) { my $ec=$1; s/\W*EC\s*$ec\W*/ /; 
  if(/\bEC[=\s]*(\d+\.(?:[\d\.\-]*))/) { my $ec=$1; s,/sw,,; s/\W*EC[=\s]$ec\W*/ /;  s/\W*EC[=\s]\d[\d\.\-]+//;
    if(grep{ $_ eq 'EC'} @NAMETAG_KEEP) { push @nametags, "EC:$ec"; }
  }
  
  # OPTION: s/\s*\([^\)]+\) *$//;  # for (Species) trailers in names; can be (any thing); but may be real name part
  
  # horrible names:
  # hydrolases, acting on acid anhydrides, in phosphorus-containing anhydrides,ATP-dependent helicases,nucleic acid binding,ATP bi...
  # tRNA (guanine-N(1)-)-methyltransferase, metazoa , tRNA (guanine-N1-)-methyltransferase, eukaryotic == duplicated
  # mannose-1-phosphate guanylyltransferase (GDP)s,GDP-galactose:mannose-1-phosphate guanylyltransferases,GDP-galactose:glucose-1-
  # ATP-dependent peptidases,nucleotide binding,serine-type endopeptidases,DNA helicases,ATP binding,damaged DNA binding,nucleosid>
  # serine/threonine kinases,protein kinases,ATP binding,sugar binding,kinases,carbohydrate binding 
  # "The protein encoded by this gene was identified as a part of pollen proteome by mass spec analysis, It has weak LEA proteins, Encodes protein phosphatase 2A B'gamma subunit, Targeted to nucleus and cytosol"
  
  if(length($_) > MAXNAMELEN) {
     my $nc= tr/[.,]/[.,]/; 
     if($nc > 1) { # keep only 1st two phrases.
     my $i= _min(index($_,',',0),index($_,'.',0));  
     while($i>0 and $i<30) { $i= _min(index($_,',',$i+1),index($_,'.',$i+1)); }
     $_= substr($_,0,$i) if($i>20);
     } 
  }       

  ## should drop these geneid == name cases: At3g18210 for arabid genes, Os12g0628400, ...
  ## replace w/ special UNK name? Unchar prot gene id
  
## ^CG\d/ || AGAP .. allow for -PA,-PB, other isoform suffix?   ACYPI009998
## zfish:  Zcg?0000 ; DDB_G0283697 ??
## human: cDNA FLJ30174 fis, clone BRACE2000975, terms-here...
## Macaca fascicularis brain cDNA clone: QtrA-18802, similar to .. UGGGGGHHHHHHHH
## Testis cDNA clone: QtsA-13078, similar to ...
## Brain cDNA, clone: QflA-14336, gene product-like protein
## mRNA, similar to triabin-like lipocalin 4a, clone: M2A11 

## .. how to deal w/ above crap?  NAME_IDPATT or other?
##   patt: m/^... clone: [\w-]+\d+, similar to (NAME HERE)

  if(not $isnamedprot and m/\b($nameidpatt)\b/) { $isid=$1; # may be 2+ patts
    if(m/^($nameidpatt)/ or  m/^($NAME_NONE)/i or m/hypothetical protein/i) { $_= $NAME_UNK; $isunk=1; } 
    elsif(m/^[cC]DNA\W* FLJ\d+ fis, clone \w+\W*$/) { $_= $NAME_UNK; $isunk=1; } # human garbage doesnt fit idpatt
    ## ^cDNA cut above ..
    ## crap: "Peptide similar hypothetical protein" from "Peptide XP_001950702 similar hypothetical protein"
    else { s/\s*\b($nameidpatt)\W*/ /g; } # FIXME:  CG123-PA, AGA123-XX, ..
    ## else { s/\s*[\(\[]?($nameidpatt)[\)\]]?\W*/ /g; } # FIXME:  CG123-PA, AGA123-XX, ..
  }  
    
  $iste= $USE_TENAME && isTEname($_);  ## add TE check option: 
  ##?? cancel isunk, but keep isput, for iste .. ignore $pi ??
  ## or maybe cancel iste if pi < minident..
  ## TE: Uncharacterized protein     update  old:Pol polyprotein     lowqual=Pol polyprotein 13%,229/1738,232
  ## iscdd: turn off $MIN_NAMEIDENT
  if($iscdd or $iste) {
  $isunk= (m/^($NAME_NONE)/i )?1:$isunk; # Unknown|Uncharacterized|Hypothetical
  } else {
  $isunk= ($pi < $MIN_NAMEIDENT or m/^($NAME_NONE)/i )?1:$isunk;  
  $isput= (!$isunk and ($pi < $MIN_CERTAIN or $isput))?1:0;
  }
  
  ## ugh: TRIGALACTOSYLDIACYLGLYCEROL 1,
  my $nuc= tr/[A-Z]/[A-Z]/;  #  uc($_) eq $_
  if(m/[A-Z]\s+[A-Z0-9]/ and ($nuc > 19 or uc($_) eq $_ ) ) { $_= ucfirst(lc($_)); } # SHOUTING phrase ...  TRICHOME BIREFRINGENCE-LIKE 19
  
  s/\s*[Hh]omolo(gy|gue|g)\s+\d+//g; s/\b[Hh]omolo(gy|gue|g)\s*//g; # ? set $isput ? set -like? # add 'ortholog'
  s/\s*[Oo]rtholo(gy|gue|g)\s+\d+//g; s/\b[Oo]rtholo(gy|gue|g)\s*//g; 
  s/^(of|with)\s+//ig; # bad leading words
  s/ and related protein[s]?//;
  ## not all# s/\s*\b[Ii]nvolved in\b\s*/ /; # CDD:
  ## bad prior CDD clean of 'AidB, Proteins involved in DNA damage response, similar to the AidB gene product'  
  ##  >> CDD: AidB, s involved in DNA damage response

  # Add protein to the end when ending in:  'binding|domain|like|related'
  unless($iscdd) { s/\b(binding|domain|like|related)\s(\W*)$/$1 protein $2/; }
  if( s/[,\s]*putative//g ) { $isput=1; } #s/putative, putative/putative/;
  
  # punctuation; move after bracket balance ??
  s/\|/:/g; s/#/n/g; s/\@/ /g; # dingbats
  s/([,;])\s*[,;]+/$1/g; ## duplicate punc left from cuts
  # s/_/ /g; # or leave uscores ?? not for ncbi ! see below NAME_NOUSCORE
  s/^\s$//; s/  +/ /g; # lead/end, extra spaces
  s/[\s\/\.,;_-]+$//; # trailing punc,spaces
  s/[.] /, /g; # no sentences?
  s/^\W+//; # no leading crap
    
  # SEQ_FEAT.ProteinNameEndsInBracket:  Phosphoenolpyruvate carboxykinase [ATP] << change [] for ncbi
  if(/\]$/) {  s/\[/\(/g; s/\]/\)/g;}
  
  # unbalanced brackets; # add {} ? not used; <> ? not brackets
  # Bug here: 'blah blah (XXX) bob bob (ZZZ' >> 'blah blah XXX) bob bob (ZZZ' unbalancd..
  # bug2: single bracket from cuts .. 'Zinc finger protein 3 (A8-5'
  # bug3:  blah ) xxx ( blah << remove )(
  if(/[\(\)\[\]]/) {
    my($nb,$ne,$d,$i,$j);
    $nb= tr/\[/\[/; $ne= tr/\]/\]/; $d=$nb - $ne;
    # while($d>0) { s/\[//; $d--; } while($d<0) { s/\]//; $d++; }
    while($d>0){ $i= rindex($_,'['); substr($_,$i,1,' '); $d--; }
    while($d<0){ $i= index($_,']'); substr($_,$i,1,' '); $d++; }
    $nb= tr/\(/\(/; $ne= tr/\)/\)/; 
    if($nb and $ne) { $i= index($_,')'); $j= index($_,'('); if($i<$j) { substr($_,$i,1,' '); $ne--; } }
    $d=$nb - $ne;
    # while($d>0) { s/\(//; $d--; } while($d<0) { s/\)//; $d++; } # BUG here? trim LAST first, 
    while($d>0){ $i= rindex($_,'('); substr($_,$i,1,' '); $d--; }
    while($d<0){ $i= index($_,')'); substr($_,$i,1,' '); $d++; }
    s/\[\s*\]//g; s/\(\s*\)//g; s/  / /g; # empties; and trail punc
    s/[\s\/\.,;_-]+$//; # trailing punc,spaces
  }   

  ## somewhere add this option:
  s/\s+\d+$// if($NAME_NODIGITS); # or s/[-\s]+\d+$//; or this? 'Subfamily 2 member A25'
  
  unless(/\w\w/) { $_= $NAME_UNK; $isunk=1; } # 'Conserved protein' becomes blank
  ## what of 'Uncharacterized N-acetyltransferase C16C4.12' >> 'Uncharacterized protein N-acetyltransferase C16C4.12'?
  if($isunk) { # regularize, but check/keep some additions
    $islike=$isput=0; 
    
    # iscdd bug: not a protein but domain, use Domain of unknown function ??
    # new: CDD: Uncharacterized protein ; old: CDD: Domain of Unknown Function DUF1041
    # new: CDD: Domain of Unknown Function; old: CDD: Topoisomerase II-associated PAT1-like protein
    # ... $islike and $iscdd bug
    
    if($iscdd or /$NAME_UNKCDD/) {  # bad now, replacing good CDD
      $_= $NAME_UNKCDD; 
      if($cddsym){ $_.=' '.$cddsym unless($cddsym=~/domain/i); $cddsym=''; } 
      elsif($isid){ $_.=' '.$isid  unless($isid=~/domain/i); $isid=0; } # append DUFnnn ?? .. append only, not prepend ??
        
    } elsif(/^($NAME_NONE)$/i or /^($NAME_NONE) protein\W*$/i) { 
      $_= $NAME_UNK;  
      
    } elsif(/^($NAME_NONE)\s*\w+/) { 
      
      if(/^$NAME_UNK\s*\w+/) {  s/^$NAME_UNK\s*//; }
      else { s/^($NAME_NONE)\s*//; s/^protein\s*//i; }
      # my $hasterm= (m/[a-z]{6,}/)?1:0;   # NCBI tsa2asn doesnt agree here.. any string is hasterm?
      my $hasterm= (m/[A-Za-z]{2,}/)?1:0;   # NCBI tsa2asn :(
      if($hasterm) { $islike++; $isunk=0; } 
      else { (my $unknoprot= $NAME_UNK) =~ s/\s*protein//; s/^/$unknoprot /; s/$/ protein/ unless(/\bprotein/); }
      # s/^($NAME_NONE)\s*/$unknoprot /;  
      # s/$/ protein/ unless(/\bprotein/); # bad?

    ## 'Uncharacterized N-acetyltransferase' becomes 'Uncharacterized protein N-acetyltransferase'
    ## .. protein should be at end
    ## UniRef50_P32084 Uncharacterized protein update  
    ##   old:Uncharacterized protein HIT-like protein Synpcc7942_1390 << nameclean.bad
    ##  orig:UniRef50_P32084 Uncharacterized HIT-like protein Synpcc7942_1390 
    ##  >> should clean to: 'HIT-like protein Synpcc7942_1390' or 'HIT-like protein' drop Unchar
    
    } else { 
      # ?? save old,clean-name as Note in some/all cases.
      if($pi >= $MIN_IDLIKE) { $islike++; $isunk=0;  } ##  for > 10-15% ident? or any?
      # Name=Mitogen-activated protein kinase kinase kinase 7-interacting protein 2
      # >> Mitogen-activated kinase kinase kinase 7-interacting protein 2-like protein
      else { $lowqualname=$_; $_= $NAME_UNK; }  # replace entirely? or not; NOT yet
    }  
  }
  #?? if(s/hypothetical\s+//i) { $isput=1; } # defer after NAME_NONE

  ## uck: old:Pol polyprotein-like new:Pol poly-like protein
  ## yuck: Uncharacterized encoded by-like protein; RepID:CU100_HUMAN; id-cut problem
  ## Cysteine rich..(chordin-like) instead of Cysteine rich..(chordin)-like ???
  ## orig: 'cysteine rich..(chordin like)' from zfish
  $islike++ if(/\blike\b/);
  $islike=$isput=0 if($iscdd or $isunk);
  if($islike and not $NAME_NOLIKE) { 
    unless(/family/ or /\blike/ or /protein kinase/ or /\Sprotein/) { s/\s+protein//; $_ .= '-like protein'; }  
    #s/\s+\-like/-like/; 
    s/\s*\bprotein.like protein/-like protein/;  
    s/\s*encoded by-like protein$/ protein/;
    ## protein like protein; above nameidpatt cut changes '-like' to ' like'
    s/[\s\-]+like\b/-like/;  
    $_=$NAME_UNK if(m/^-like/);
  }
  $isput=0 if($NAME_NOPUTATIVE or /\blike\b/);
  s/\bprotein protein/protein/ig; # other stutters? 'ribonucleoprotein protein' ?
  # .. UniProt uses: MPP10_HUMAN; U3 small nucleolar ribonucleoprotein protein MPP10 ; O00566  MPHOSPH10
  s/\s*\bprotein$// if(/^[Pp]rotein/); # and? $isnamedprot ; dont need both  
  s/^\W+//; # no leading crap, again

  $_.=' '.$isid if($isid and $NAME_KEEPID);
  # FIXME: NAME_KEEPID has 2 meanings: a. Name == ID instead of Unknown; b. ID of protein
  # FIXME: putative drop/add messes diff

    ##.. FIXME.1510: $isnamedprot should keep leading 'protein XXX' to defeat NCBI name checka whine about no dbids..
    ##.. now is $isid? problem prior-cleaned "Protein XXX" are now missing Protein prefix, add back if $name~=m/^[A-Z][A-Z\d.-]+/ ie no lc(name)
    ##.. should use dbxref check of Uniprot valid names
    ## Ugh: "CDD: Protein Domain of unknown function Domain Domain" here, from 
  if($NAME_AddBackProtprefix and not($iscdd or $iste or $isunk or m/protein|domain/i) and ($isid or m/^[A-Z][A-Z_\d\s\.\-]+$/)) {
    # above# $_.=' '.$isid if($isid);  $isid=0;
    $_= 'Protein '.$_;
  }
  
  if($NAME_NOUSCORE and m/_/) { # ncbi hates these :( mostly dbids
    # if($whatcheckherefordbid) { $isunk=1; $_= $NAME_UNK; }
    s/_/ /g;
  }
  

  s/\s*$//; s/^\s$//; s/  +/ /g; # lead/end, extra spaces, again..
  s/[\s\/\.,;_-]+protein$/ protein/; # " punc,spaces protein" 
  if($iscdd){ $_= $cddsym.', '.$_ if($cddsym); $_= 'CDD: '.$_; }   # $_ = "CDD: $_" if($iscdd);  
 
  my $diff=($nameprechange ne $_) ? NAMEDIFF_MINOR :0;
  if($diff and $isput and ($nameprechange eq "$_, putative")) { $diff=0; }
  if($diff) {
    my($nat,$net)=($nameprechange,$_); 
    map{ $_=lc($_); s/\W+$//; tr/[a-z0-9]/ /c; } ($nat,$net);    
    $diff= NAMEDIFF_MAJOR if($nat ne $net); # $diff++ 
  }  
    
  s/^([a-z])([a-z][a-z])/\u$1$2/ if($NAME_UCFIRST); # upcase 1st let
  $_ = 'TE: '.$_ if($USE_TENAME and $iste); # want this?? or return $iste ? add before diff test
  $_ .= ', putative'  if($isput and not $isunk); ## add $NAME_NOPUTATIVE option
  $_ .= ' '.$locusname if($isunk and $locusname);
  foreach my $kv (@nametags) {  $_.='; '.$kv; }

  ## ncbi complains about this ^^ locusname addition; leave out at least w/ config.
  return wantarray ? ($_, $lowqualname,$diff) : $_;
}


sub idnameclean { 
  my ($namin, $pi, $locusname)= @_;
  my $id=''; if($namin =~ s/^(\w\S+)\t//) { $id=$1; }
  my($naclean,$loqualname,$diff)= nameclean($namin,$pi,$locusname);
  return($id,$naclean,$loqualname,$diff);
}
  


=item blast2names
  was unipblast2names

  input blasttab == blast table FROM blastp -outfmt 7 via evigene/scripts/makeblastscore2.pl
  ?? change to allow blastp fmt 6,7 input, but will miss query,ref sizes.

  QueryID  RefID  Bitscore Ident Align QueryLen RefLen
  sorted by QueryID >> +Bitscore by default; optional sort QueryID >> +Align
  
  cat $refna-$pt.tall4 | grep -v '^Query' | sort -k1,1 -k5,5nr -k3,3nr -k2,2 | env na=$refnames perl -ne\
  'BEGIN{open(F,$ENV{na}); while(<F>){ chomp; ($rd, $rdesc)= split" ",$_,2; $rdesc{$rd}=$rdesc; } } \
  ($td,$rd,$bits,$iden,$aln, $tlen, $rlen)=split; \
  if($ld and $td ne $ld) { put1(); $desc1=$lcd=""; } \
  $rdesc=$rdesc{$rd}; $named=(not $rdesc or $rdesc =~ /uncharacterized protein/i)?0:1; \
  if($td ne $ld) { $lcd = "$rdesc; aln=$aln" if($named); $aln1=$aln; $desc1="$rdesc; aln=$aln"; } \
  elsif( not $lcd and $named and $aln >= 0.6*$aln1) { $lcd= "$rdesc; aln=$aln ;; $desc1"; } \
  $ld=$td; END{put1()} sub put1{ $lcd=$desc1 if(not $lcd and $desc1); print "$ld\t$lcd\n"; }' > $refna-$pt.namea

=cut

use constant { BLSCORE_TABLE => 4, BLAST_TABLE => 7,
        NAMEOUT_CONSENSUS => 3, NAMEOUT_TABLE => 2, NAMEOUT_LIST => 1, 
        CDD_NOCONSENSUS => 1,
        };

sub blast2names # was unipblast2names
{
  my( $blasttab, $refnames, $sortblast, $outnames, $cddinfo, $blformat, $oformat )= @_;  
  our($ok,$inh,$outh,$oform,$desc1,$aln1,$ld,$lcd,$nin,$nout,$lastref);
  our(@lname,@lcdd,%tdspan);
  
  #o my $NOCLEAN= not $NAMED_ISCLEAN;  #2017:wrong way? 
  my $NOCLEAN= $NAMED_ISCLEAN;   
  $NOCLEAN=1 if($refnames =~ /clean/); # fixme hack
  my $ONENAME=1; # need opt, samename, put single name (+ CDD if there), first align, despite align/other best stats for 2nd+
  
  $blformat  ||= 0; # BLSCORE_TABLE; need input option
  $sortblast ||= ""; # FIXME: sort now only for BLSCORE_TABLE 
  $oform = NAMEOUT_TABLE; # new default 2017
  if($oformat) { 
    $oform= ($oformat =~ /con|3/) ? NAMEOUT_CONSENSUS 
          : ($oformat =~ /tab|2/) ? NAMEOUT_TABLE 
          : ($oformat =~ /list|1/) ? NAMEOUT_LIST 
          : NAMEOUT_TABLE; #? NAMEOUT_TABLE default instead, was NAMEOUT_LIST
  }
  $nin=$nout=$aln1=0; $lastref=$desc1=$lcd=$ld="";  @lname=(); %tdspan=(); @lcdd=();
  
  #FIXME: add CDDnames for deltablast results w/ -rpsdb $cddb -show_domain_hits
  ##  cddnames = CDD:idnum  messyname...
  
  my($ncdd,$cddesc,$cdlen,$cdrepid)=(0,{},{},{});
  if($cddinfo) {
    ($ncdd,$cddesc,$cdlen,$cdrepid)= readnames($cddinfo, 1, $NOCLEAN); # = cddnames($cddinfo);
  }
  
  my($nrefname,$refdesc,$reflen,$rrefid)=(0,{},{},{});
  ($nrefname,$refdesc,$reflen,$rrefid)= readnames($refnames, 0, $NOCLEAN);
    # = readnames($refnames); # need opt to call nameclean or expect cleaned names
    
  unless($ncdd > 0 or $nrefname > 0) {
    warn "ERR: No names for blast2names, in $refnames or $cddinfo\n"; 
    return -1; # or cancel sort? $sortblast=0
  }  

  # blasttab == blast table FROM blastp -outfmt 7
  # want sort options here? default blastp is sorted right:  queryid > +bitscore 
  #  option for queryid > alignscore
  # sort ONLY for $blformat == BLSCORE_TABLE : open and test format?
  # FIXME: dont need to sort all of file, only each query set.  also that can be done to blastp.out
  # .. align sort == best naming looks likely.
  
  # if($blformat == 0 and $blasttab =~ /\.blastp/) { $blformat= BLAST_TABLE; } # ?? safe from name?
  # if($blformat == 0 and $blasttab =~ /\.tall|\.bltab/) { $blformat= BLSCORE_TABLE; } # ?? safe from name?
  
  my $isgzinput= ($blasttab =~ /\.gz/)?1:0;
  my $CAT = ($isgzinput) ? "gunzip -c " : "cat ";
  if($sortblast and $blformat != BLSCORE_TABLE) {
    warn "BUG: need input blasttable format to sort, for now"; 
    return -1; # or cancel sort? $sortblast=0
    #? die "BUG: now need input blasttable format to sort";
  }
  if($sortblast =~ /align/i) {
    $ok= open($inh, "$CAT $blasttab | sort -k1,1 -k5,5nr -k3,3nr -k2,2 |"); # QueryID > +Align > +Bitscore > RefID
  } elsif($sortblast =~ /iden/i) {
    $ok= open($inh, "$CAT $blasttab | sort -k1,1 -k4,4nr -k3,3nr -k2,2 |"); # QueryID > +Ident > +Bitscore > RefID
  } elsif($sortblast) {
    $ok= open($inh, "$CAT $blasttab | sort -k1,1 -k3,3nr -k2,2 |"); # QueryID > +Bitscore > RefID
  } else {
    if($isgzinput) { $ok= open($inh, "$CAT $blasttab |"); }
    elsif($blasttab =~ /stdin|^-/) { $inh= *STDIN; $ok=1; }
    else { $ok= open( $inh, $blasttab); }
  }
  unless($ok) { warn "ERR: blast2names open $blasttab"; return -1; }

  #? $outnames="stdout" unless(defined($outnames)); 
  unless($outnames) { $outnames=$blasttab; $outnames =~ s/\.*//; $outnames.=".named"; }
  if($outnames =~ /stdin|^-/) { $outh= *STDOUT; $ok=1; }
  else { 
    rename($outnames,$outnames.".old") if(-f $outnames); 
    $ok= open($outh, '>', $outnames); 
  }
  unless($ok) { warn "ERR: blast2names write $outnames"; return -1; }

  ## FIXME: output 3+ columns, add pctalign/ident AND RefID;use syntax as before for namepi at end of name
  ##   my $namepi= (s/\s+\((\d+)%.*\)//)?$1:0; == Name (66%U) .. 
  ## maybe use multi rows per gene for alternate names rather than ;; concat, as with CDD:
  ##  TrID1 \t somedefhere (55%CD) \t  120align/240rlen,300qlen \t  CDD:xxx
  ##  TrID1 \t somedefhere (25%CD) \t   25align/100rlen,300qlen \t  CDD:yyy
  ##  TrID1 \t somedefhere (50%UP) \t  200align/400rlen,300qlen \t  UniRef50:zzz
  ##  TrID1 \t Uncharprot  (83%UP) \t  250align/300rlen,300qlen \t  UniRef90:qqq
  ## should use diff method for CDD names, not consensus.
  ##  .. want all non-overlap domains, or longest align if "consensus"?
  # 2015..17 expected names table format, changed default outformat
  #  my($td,$name,$alnscore,$rd,@more)=split"\t"; 
  #  # alnscore format expected: 72%,3270/4555,3282 ;  may be '72' or '72%' only
  
  ## FIXME/Option?: add @lcdd list of all CDD ids/prot; no, more work w/ %tdspan needed
  sub put1name{ our($outh,$desc1,$ld,$lcd); $lcd=$desc1 unless($lcd); print $outh "$ld\t$lcd\n"; }
  sub put2name{ our($outh,$desc1,$ld,@lname); push @lname,$desc1 unless(@lname); map{ print $outh "$ld\t$_\n"; } @lname; } 
  sub put3name{ our($outh,$desc1,$ld,@lname,%tdspan); my($cname); 
    if(@lname>1) { ($cname)= consensusname(\@lname,1,CDD_NOCONSENSUS,\%tdspan); }
    elsif(@lname==1) { $cname=$lname[0]; } else { $cname=$desc1; }
    print $outh "$ld\t$cname\n"; }
    
  sub putname { our($oform,$nout,@lname,%tdspan);  
  
  	    ## add non-overlap collect ids, for CDD only?
  	    ## this isn't working or no nonover cdd in tests ??
  	my($cdd1,$cid1,$icd); 
  	for($icd=0;$icd<@lname;$icd++) { if($lname[$icd]=~/CDD:/) { 
  		$cdd1=$lname[$icd]; ($cid1)=$cdd1=~/(CDD:[\d.]+)/; last; } }
  	if($cid1 and $tdspan{$cid1}) {  my @addid=();
  		my($spanb,$spane)=split",",$tdspan{$cid1};# @{$tdspan{$cid1}};
  		foreach my $rd (keys %tdspan) {
  			my($tb,$te)= split",",$tdspan{$rd}; # @{$tdspan{$rd}}; 
  			if($tb >= $spane-5 or $te <= $spanb+5) {
  			push @addid, $rd; $spanb= _min($spanb,$tb); $spane= _max($spane,$te);
  			}
  		if(@addid) {
  		  my $cid2= join ",", @addid; $cdd1=~s/$cid1/$cid1,$cid2/;
  		  $lname[$icd]= $cdd1;
  		  }
  		} 
  	}
  	# if(@lcdd>0 and @lname) { my $cdd2=join(",",@lcdd);  # cdd2 add 13aug12 : NOT here, @lcdd empty, @lname has all
  	#  foreach (@lname) { s/(\tCDD:[\d.]+)/$1,$cdd2/ if(m/\tCDD:\d/); } 
  	# }

    ($oform == NAMEOUT_CONSENSUS) ? put3name() :
        ($oform == NAMEOUT_TABLE) ? put2name() : put1name(); $nout++; 
	}
  
  my($td,$rd,$bits,$iden, $aln, $tlen, $rlen, $tlenid)=(0) x 10;
  while(<$inh>) {
    ## if missing tlen, pick tlen from blast # Query
    ## # Query: dmag4vel4xcak45Loc1t17560 aalen=51,67%,complete; clen=231; strand=-; offs=201-46;
    
    if(/^\W|^Query/) {
      ## $tlen= 0;
      if($blformat != BLSCORE_TABLE and /^# Query:/) {
        $tlenid= $tlen=0; my($qid)= m/# Query: (\S+)/;
        if(m/aalen=(\d+)/) { $tlen=$1; $tlenid=$qid; } elsif(m/\blen=(\d+)/) { $tlen=$1; $tlenid=$qid; }
      }
      next;
    }
    
    ($td,$rd,$bits,$iden, $aln, $rlen)=(0) x 10; # zero all but tlen/row
    my @row= split;
    if($blformat == 0) { my $ncol= @row;
      if(($ncol == 7 or $ncol == 5) and $row[2] =~ /\d/) { $blformat= BLSCORE_TABLE; }
      elsif($ncol == 12 and $row[11] =~ /\d/) { $blformat= BLAST_TABLE; }
      warn "# blast2names detected blast format $blformat\n"; # if $DEBUG
    }
    if($blformat == BLAST_TABLE) {
      my($mis,$pctid,@bspan);  
      ##FIXME? add bspan overlap checker to allow non-over hits (CDD at least); BLSCORE_TABLE (tall4) has that.
      ($td,$rd,$bits,$aln,$mis,$pctid,@bspan)= @row[0,1,-1,3,4,2, 6,7,8,9]; # 6-9 =  q. start, q. end, s. start, s. end, 
      $iden= _max(0,$aln-$mis);  ## FIXME: check, is (aln-mis)/aln same as pctid or not, may depend on blast
      $tlen=0 unless($td eq $tlenid);

      # FIXME for CDD : collect all non-overlap domains; ditto nonCDD ??
      ## not working, via consensus() .. add here, putname ? not working there or no such cases?
      my($tb,$te,$rb,$re)= @bspan;  ($tb,$te)=($te,$tb) if($tb>$te);
      ## pick only 1st span??
      unless($tdspan{$rd}) { $tdspan{$rd}="$tb,$te"; } ## [$tb,$te]; old
      # else { $tdspan{$rd}->[0]= $tb if($tb < $tdspan{$rd}->[0]);
      #	$tdspan{$rd}->[1]= $te if($te > $tdspan{$rd}->[1]); }
      
    } elsif($blformat == BLSCORE_TABLE) {
      ($td,$rd,$bits,$iden,$aln, $tlen, $rlen)=@row;
    } else {
      ($td,$rd,$bits,$iden,$aln, $tlen, $rlen)=@row;
      if($rd and $bits =~ /\d/) { $blformat= BLSCORE_TABLE; }
      else { warn "ERR: blast2names input blast format unknown $blasttab"; return -1; }
    }
 
    $rd =~ s/^gnl\|CDD\|/CDD:/;
    my $samename= ($td eq $ld) ? 1: 0;
    ## FIXME: CDD is blocking rdesc name as td eq ld..
    if($lastref =~ /^CDD:/ or $rd =~ /^CDD:/) {
      $samename=0 if(($lastref =~ /^CDD:/ and $rd !~ /^CDD:/) or ($lastref !~ /^CDD:/ and $rd =~ /^CDD:/));
    }
    
    unless($samename) { putname() if($ld); $aln1=$desc1=$lcd=""; @lname=(); %tdspan=(); @lcdd=(); $nin++; }

    ## FIXME: add aln/tlen aln/rlen tests for min score..
    ## FIX2: annotate with "aln=$aln/$rlen,$tlen" ?? aln/max(rlen,tlen) ??
    ## FIXME: handle sorting here; if($samename) .. collect all input rows, keep $bits,$iden,$aln,
    
    my ($rdesc,$repid); 
    if($rd =~ /^CDD:/) { ## and $ncdd > 0
      $rdesc= $cddesc->{$rd} || $rd;  $repid= $cdrepid->{$rd}||"";
      $rlen= $cdlen->{$rd} || $rlen;
      # need to treat these diff; separate output set.. 2 names per td gene?
      # .. collect all CDD names per gene for output? or just top score?
     } else {
      $rdesc= $refdesc->{$rd} || "";  $repid= $rrefid->{$rd}||"";
      $rlen= $reflen->{$rd} || $rlen;
    }
    # if($rlen==0 and $rdesc =~ /\blen=(\d+)/) { $rlen=$1; } # drop this
    
    my $palign=0;
    if(0) { # this is not good; for larger query prot, paln is same hiding best ref
      my $mlen= _max($rlen,$tlen); 
      $palign= ($mlen > 0) ? int(0.5 + 100*$aln/$mlen) : 1; # $NAMED_MINALIGN; ?? no default
    } else {  # ref-palign, best percent but >aln size for <palign hides true best ref sometimes.
      $palign= ($rlen > 0) ? int(0.5 + 100*$aln/$rlen) : 1; # $NAMED_MINALIGN; # use only rlen?
    }
    $palign=100 if($palign>100);
    # FIXME: NAME_KEEPID has 2 meanings: a. Name == ID instead of Unknown; b. ID of protein
     
    my $named=(not $rdesc or $rdesc =~ /$NAME_UNK/i)?0:1; 
    my $refid= ($rd =~ /:/)? $rd : "RefID:$rd";
    unless($rdesc) { $rdesc = $NAME_UNK." $refid" ; } #? or not
    elsif($NAME_KEEPID) { $rdesc .= " $refid"; } # FIXME: what option
    
    if($oform > NAMEOUT_LIST) { ## NAMEOUT_TABLE
      ##  TrID1 \t somedefhere (50%UP) \t  200align/400rlen,300qlen \t  UniRef50:zzz
      my $desc= join("\t",$rdesc,"$palign%,$aln/$rlen,$tlen",$refid,$repid);
      unless($samename) { $aln1=$aln; $desc1=$desc; }
      if($named) {
        $ok= ((not $samename) or ($palign >= $MIN_CERTAIN) or (! $lcd and $aln >= $NAMED_MINALIGN * $aln1)) ? 1:0;
        $ok=0 if($ONENAME and $samename and @lname);
        if($ok) { push @lname, $desc; $lcd=$desc unless($lcd); } 
        ## elsif($samename and $refid=~/CDD:/ and ($palign >= $MIN_CERTAIN)) { push @lcdd, $refid; } # these add to 1st cdd above
        ## ^^ not so good here, need non-overlapped domains, most are overlap same prot frag.
      } 
    } else {
      unless($samename) { $lcd = "$rdesc; aln=$aln" if($named); $aln1=$aln; $desc1="$rdesc; aln=$aln"; } 
      elsif( not $lcd and $named and $aln >= $NAMED_MINALIGN * $aln1) { $lcd= "$rdesc; aln=$aln ;; $desc1"; } 
    }
    
    $ld=$td;  $lastref=$rd;  
  } 
  putname() if($ld);
  close($inh);close($outh);   
  return($nout, $outnames);
}


sub cleannamefile  
{
  my( $namefile,  $outnames, $iscdd)= @_;  
  unless($outnames) { $outnames= $namefile; $outnames .= ".clean"; }
  use constant NOCLEAN => 0;
  my($nrefname,$refdesc,$reflen,$rrepid) = readnames($namefile, $iscdd, NOCLEAN);
  open( my $outh, '>', $outnames) or do { warn "ERR:writing $outnames"; return; };
  foreach my $id (sort keys %$refdesc) {
    my $de= $refdesc->{$id} || $NAME_UNK; # change "" to Unchar
    my $repid= $rrepid->{$id}||"";
    my $rlen= $reflen->{$id}||0;
    
    ## cdd input format is this:  preserve it ? or change
    ##  ($rid,$rdesc) = "id\s+desc" ;  my($roid,$rshort,$rlong)= split /,\s*/,$rdesc,3;  
    if($iscdd) { # or m/^CDD: /
    $de =~ s/CDD: //; # is this right?
    $de .="; len=$rlen" if($rlen); # should this have tab not ; ?
    print $outh "$id\t$repid, $de\n"; #? tab not comma, ?
    } else {
    $de .="\t$repid" if($repid); # tab not space
    $de .="; len=$rlen" if($rlen);
    print $outh "$id\t$de\n";
    }
  } close($outh);
  return($nrefname,$outnames);
}

=item cleannamefile bugs

< DRERI:ENSDARG00000017562      Uncharacterized protein; len=508
> DRERI:ENSDARG00000017562      UniProt:IPR018379; len=508 << missing Unchar, why?
  ...
  
=cut

# sub cddnames: reuse for uniprotnames
sub readnames {
  my($namefile,$iscdd, $NOCLEAN)= @_;
  my($inh,%rdesc,%rlen,%repid,$nin,$ok); $nin=0; $iscdd||=0;
  %rdesc= %rlen= ();
  unless($namefile and -f $namefile) { $ok=0; }
  elsif($namefile =~ /.gz/) { $ok= open($inh, "gunzip -c $namefile |"); }
  else { $ok= open( $inh, $namefile); }
  unless($ok) { warn "ERR: readnames open $namefile"; return($nin, \%rdesc, \%rlen); }
  
  while(<$inh>) { 
    next if(/^#/); chomp; 
    #o# my($rd, $rdesc)= split" ",$_,2; 
    my($rd,$rdesc,@xtra)= split /\t/,$_; #look for tabs first# 
    ($rd, $rdesc)= split(" ",$_,2) unless($rdesc);
    $rd=~s/^>//; # may be fasta hdr
    my $rlen=0; 
    for ($rdesc,@xtra) { if(m/\blen=(\d+)/) { $rlen=$1; s/\W*len=$rlen\W*/ /; last; } }
    #o if($rdesc =~ m/\blen=(\d+)/) { $rlen=$1; $rdesc =~ s/\W*len=$rlen\W*/ /;  }
    $rlen{$rd}= $rlen; # even if missing?
    $nin++;
   
    if($iscdd) {
      my($roid,$rshort,$rlong)= split /,\s*/,$rdesc,3;  
      $rshort||=""; $rlong||=""; 
      $rlong =~ s/\s*[\.\;\[].*$//;

      if($rlong) {
        (my $rsex=$rshort) =~ s/[\(\)\[\]]/./g;
        if($rlong eq $rshort) { $rlong=""; }
        elsif($rlong) { $rlong =~ s/[\(]$rsex[\)]//; }
      # DUF3552, Domain of unknown function (DUF3552) ; but not SEA, SEA domain
      ## fixme CDD special cases per source?
      ## cd0nnn have some wordy preambles to drop: n=38
      ##  This family is most closely related to the (good stuff here)
      ##  This group of proteins belong to a xxx family of;
      ## A group of hypothetical ..

      $rlong =~ s/^This family is most closely related to the //;
      $rlong =~ s/^This group of proteins belong to a \w+ family of //;
      $rlong = nameclean($rlong) unless($NOCLEAN);
       
      if(length($rlong) > CDD_MAXNAMELEN) { $rlong =~ s/\s*[,].*$//; }
      my $toolong= (length($rlong) > CDD_MAXNAMELEN)?1:0;
      while(length($rlong) > CDD_MAXNAMELEN) {
        my $sp= rindex($rlong,' '); last if($sp<1);
        $rlong= substr($rlong,0,$sp);
        }
      if($toolong) { # chomp trailing crap words: on|is a|a|or|found in|in| of the|  
        $rlong =~ s/\s+(are a|is a|a|found in|in|of the|of|on|includes all|all|are)$//;
      } 
        if($rlong =~ m/\b$rsex\b/) { $rdesc=$rlong; }
        elsif($roid eq $rshort and $rlong) { $rdesc=$rlong; }
        else { $rdesc=$rshort; $rdesc.=", $rlong" if($rlong); } ## maybe dont want rshort ..
      } else { # not rlong
        $rdesc= $rshort; 
      }
      # $rdesc{$rd}="CDD: $rdesc, $roid";  # put roid in %repid instead
      $rdesc{$rd}="CDD: $rdesc";  #?? is this CDD: prefix desired?
      $repid{$rd}=$roid;
      
    } else {
      # chop this: RepID:DPOG1_HUMAN or keep? want _HUMAN part ; all?
      ## FIXME: here, where? change useless RepID:, RefID: dbx to standard tags: UniProt: .. 
      ## need option or check ID for known patts: 'UniRefnnxxxx' 'xxxxx_ZZZZZ' ?
      ## FIXME: alt ids follow \tab for some refnames ..
      my $repid="";
      #? my $REFTAG='RepID|DbXref';
      for ($rdesc,@xtra) { if(m/\b(RepID:[\w\.-]+)/) { $repid=$1; s/\W*(RepID:[\w\.-]+)\W*//; last; } }
      unless($repid) { for (@xtra) { if(m/^(\w+:[\w\.-]+)/) { $repid=$1; last; } } }
      $repid{$rd}= $repid if($repid);
      #o if($rdesc =~ s/\W*(RepID:\w+)\W*//) { $repid{$rd}=$1; }
      #o elsif($rdesc =~ s/\t(\w+:[\w\.-]+).*//) { $repid{$rd}=$1; }
      
      $rdesc = nameclean($rdesc) unless($NOCLEAN);
      $rdesc{$rd}= $rdesc;
    }
  } close($inh);
  return($nin, \%rdesc, \%rlen, \%repid);
}

sub cddnames {
  my($cddinfo)= @_;
  return readnames($cddinfo,1);
}  

=item cddblast2names : merge w/ above

  add CDD naming from same deltablast results?
  deltablast -rpsdb $cddb -show_domain_hits -evalue 1e-5 -outfmt 7 -db $refdb -query $qfile -out $onam.deblastp
  namegenes -cddnames=info.cdd.txt -refnames uniref.names -blast xxx-uniref.deblastp

## CDD
    ## CDD:176950 CHL00005, rps16, ribosomal protein S16.; len=82
    cat $refna-$pt.tall4 | env na=../aaset/info.cdd.txt perl -ne\
'BEGIN{open(F,$ENV{na}); while(<F>){ chomp; ($rd, $rdesc)= split" ",$_,2;
($roid,$rshort,$rlong)= split /,\s*/,$rdesc,3;  $rlong =~ s/\s*[\.\;\[].*$//;
if(length($rlong)>29) { $rlong =~ s/\s*[,].*$//; $rlong=substr($rlong,0,30).".." if(length($rlong)>29); }
if($roid eq $rshort){ $rdesc=$rlong; } else { $rdesc=$rshort; $rdesc.=", $rlong" if($rlong); }
$rdesc{$rd}="$rdesc, $roid"; } }
($td,$cid,@val)=split; $rdesc=$rdesc{$cid}; if($ld and $td ne $ld) { print "$ld\t$lcd\n"; $lcd=""; }
$lcd .= "$rdesc; " if($rdesc); $ld=$td; END{  print "$ld\t$lcd\n"; }' > $refna-$pt.named

=cut

=item various names == gene IDS

  my $nameidpatt= $public_options{nameidpatt} || '(Os|At|AT)\d{1,2}g\d{3,}|(GLEAN_\d+|G[A-Z]\d\d+)'; #plants;
  my $nameidpatt= '(?:At|AT|Os)\d{1,2}[Gg]\d{3,}|OSJNB\w\d+[\w.]+|GLEAN_\d+|G[A-Z]\d\d+'; # plants + bugs;  
  my $nameidpatt= '(GLEAN_\d+|G[A-Z]\d\d+)'; #bugs;  

=cut

#.........
sub _min { return ($_[1] < $_[0]) ? $_[1] : $_[0]; }
sub _max { return ($_[1] > $_[0]) ? $_[1] : $_[0]; }


# consensusname from evigene/scripts/omcl/arp_condef2.pl
sub consensusname
{
  my($nameidlist,$cleanednames,$cddnoconsensus,$tdspan)= @_;
  ## flag for CDD:? dont try this consensus for those?
  $tdspan ||= {};    # not used now
  $cleanednames||=0; # not used now
  $cddnoconsensus||=0;
  
  our $NAME_NONE_PLUS = $NAME_NONE.'|putative|like|protein';
  our($goodannots,$poorannots)= ($NAME_GOODSPP,$NAME_POORSPP); # options ..

  # should have option/config for goodannot weights:  HUMAN=>1.5, MOUSE=>1.4, DROME=>1.1, ... JunkSpp => 0.5  
  # 2017upd: fix so default sppwt not used (HUMAN) when caller supplies good/bad set; got too many human..
  my %sppwt=(); # 2017upd
  my %OLD_sppwt= (HUMAN => 2.0, MOUSE => 1.5, DANRE => 1.2,
              DROME => 1.3, CAEEL => 1.1, ARATH => 1.1,
              #old# IXOSC => 0.5, # bad is what? <1 or <0
              );
  while( $poorannots =~ m/(\w+)[:=]([\d\.]+)/g) { my($k,$v)=($1,$2); $sppwt{$k}=$v; $poorannots =~ s/$k[:=]$v/$k/; }
  while( $goodannots =~ m/(\w+)[:=]([\d\.]+)/g) { my($k,$v)=($1,$2); $sppwt{$k}=$v; $goodannots =~ s/$k[:=]$v/$k/; }
  my $sppweights= join('|',sort keys %sppwt);
  if($poorannots) { $poorannots.="|CDD:"; } else { $poorannots="CDD:"; } # not consensus unless only annot
  
  my(%de,%palign,%balign,%inputname,@nonoverids); 
  my $firstid=""; my($spanb,$spane)=(0,0);
  foreach my $natab (@$nameidlist) {
      ##  natab == somedefhere (50%UP) \t  200align/400rlen,300qlen \t  UniRef50:zzz \t RepID:xxxx
      ## Lipase  80%,375/375,470 RefID:UniRef50_Q658L8   RepID:Q658L8_HUMAN
    my($de,$palign,$id1,$id2)=split"\t",$natab;
    map{ $_||="" } ($palign,$id1,$id2);
    my $id= ($id1=~/CDD:/) ? $id1 : ($id2)? $id2 : $id1; # for goodannots: uniprot ID expected in id2, not for CDD:id
    my $iscd=($id=~/CDD:/)?1:0;
    next if($inputname{$id}); #CDD many same
    $inputname{$id}= $natab; # return best from this
    $firstid=$id unless($firstid); # id or id1 here?
    
## not working.    ## add non-overlap collect ids, for CDD only?
#   	if($tdspan->{$id1}) { 
#   		my($tb,$te)= @{$tdspan->{$id1}}; 
#   		if($spane==0) {
#   			push @nonoverids, $id; ($spanb,$spane)=($tb,$te);
#   		} elsif($tb >= $spane-5 or $te <= $spanb+5) {
#   			push @nonoverids, $id; $spanb= _min($spanb,$tb); $spane= _max($spane,$te);
#   		}
#   	}
    
    my $piname= $palign;
    # CDD has '99%,180/181,115-294' namepct, trspan replaces trlen
    if($iscd) {
      # skip balign ? CDD superceeding good other annots .. why?
      $balign{$id}=1; $piname= $palign{$id}= ($palign=~/^(\d+)/)?$1:$MIN_NAMEIDENT
      
    } elsif($palign=~/^(\d+)%,(\d+).(\d+),(\d+)/) {
       my($pi,$aln,$rl,$tl)=($1,$2,$3,$4); 
       my $wt= $aln;
       #?or# my $wt= int(100*$aln/$rl);
       #? $wt = $wt * 1.2 if($id =~ /$goodannots/); # probably want this
       $balign{$id}= $aln;
        ## need rl,tl to avoid high palign, low size agreement.. or use _max(rl,tl) here, w/ balign?
       $tl=$rl if($iscd); # not here tho..
       $piname= $palign{$id}= _min(100, int(100*$aln/_max($tl,$rl)));
       # $palign{$id}= int(100*$aln/$rl);
       
    } elsif($palign=~/^(\d+)/) { 
      $piname= $balign{$id}= $palign{$id}= $1;  
    } else {
      $piname= $balign{$id}= $palign{$id}= $MIN_NAMEIDENT; #? missing data
    }
    
    $de =~ s/^[A-Z]+:\s+//; # CDD: other prefix?
    
    ## FIXME: call nameclean() always, with palign, to add -like, .. other flags from align
    # ($de)= nameclean($de) unless($cleanednames); 
    my $loqualname;
    ($de,$loqualname)= nameclean($de,$piname);
    
    $de{$id}= $de if($de =~ /\w+/);  
  }
  
	## CDD fixme: collect all CDD:ids of non-overlapped aligns.. need more work on blast bspan parsing
  if($firstid =~ /^CDD:/ and $cddnoconsensus == CDD_NOCONSENSUS) {
  	my $name= $inputname{$firstid};
#   	if(@nonoverids>1) { 
#   		my $cdd2= join",", grep{ $_ ne $firstid } @nonoverids; $name=~s/$firstid/$firstid,$cdd2/; 
#   	}
  	return $name;
	}
	
## revise this: best name should be some compound of goodspecies * max-alignbases * common-terms 
##  -- dont pick name by just one criteria, i.e. max-alignbases can be close spp w/ poor names 
##  -- common-terms wt should not prefer long names to short, but should prefer terms that are common to several
##  -- where goodspecies align is high, should prefer that name; for weaker align, prefer common-terms name?
##..........
## FIXME want 2 pids to break ties, or just nbase aligned to break ties?
# 80%,375/375,470 RefID:UniRef50_Q658L8   RepID:Q658L8_HUMAN
# 80%,375/398,470 RefID:UniRef50_P07098   RepID:LIPG_HUMAN
# 80%,376/423,470 RefID:UniRef50_Q5VYY2   RepID:LIPM_HUMAN << longest align
#.. which is best? max-alignbases or max-percentalign? probably max bases
# < dmag4vel4xbnk65Loc10201t7     Coagulation factor XI   71%,175/210,245 RefID:UniRef50_H0Y596   RepID:H0Y596_HUMAN
# > dmag4vel4xbnk65Loc10201t7     Coagulation factor XI   34%,193/573,245 RefID:UniRef50_E9PGP2   RepID:E9PGP2_HUMAN
#.. smad2 vs smad4
# < dmag4vel4xbnk65Loc1901t2      Mothers against decapentaplegic 100%,489/467,484  RefID:UniRef50_Q15796   RepID:SMAD2_HUMAN
# > dmag4vel4xbnk65Loc1901t2      Mothers against decapentaplegic 93%,512/552,484 RefID:UniRef50_Q13485   RepID:SMAD4_HUMAN    

## BUT top align picks close relatives that often have poor qual names .. just best align != best name
# < dmag4vel4xbnk65Loc3665t3      Aminopeptidase N   94%,965/967,1028  RefID:UniRef50_P15144   RepID:AMPN_HUMAN
# > dmag4vel4xbnk65Loc3665t3      Aminopeptidase-like protein 98%,1009/948,1028   RefID:UniRef50_D6WCL5   RepID:D6WCL5_TRICA

  ## pull one readable/valid name from these; not all agree
  # Alpha 1,3-fucosyltransferase, putative
  # Alpha-(1,3)-fucosyltransferase 10
  my (%kw, %kwd, %kwcomm, %wtd, %didspp);
  foreach my $id (sort keys %de) {
    my @w= split /\W+/, lc($de{$id}); # words
    # my @w= grep /^[a-z]\w/, split /\W+/, lc($de{$id}); # words, skip numbers or not?
    my $wt= $palign{$id} || 1; # palign or balign ?
    my %wseen=(); foreach (@w) { 
      unless($wseen{$_}++ or m/^($NAME_NONE_PLUS)/i) { $kwcomm{$_}++; $kw{$_}+=$wt; $kwd{$id}{$_}+=$wt; }
    }
  }

## new algo problem using algn-bases w/o reflen : longer aln wins when minor part of reflen
## add: bok=0 if($align/$reflen < mingoodaln) OR bok=1 if(aln < bestalign but aln/reflen >> bestaln/bestref and baln/bref << xxx)

## v15a
#old: dmag4vel4xbnk65Loc12036t1       Alpha-carbonic anhydrase        81%,286/354,352 RefID:UniRef50_E9FX53   RepID:E9FX53_DAPPU
#new: dmag4vel4xbnk65Loc12036t1       Receptor-type tyrosine-protein phosphatase zeta 13%,290/2315,352        RefID:UniRef50_P23471   RepID:PTPRZ_HUMAN
## v15b : this is not improvement (balign/palign effect?) 100% palgn to tiny ref not as good as >> balign to right-sized ref
## also should weight IXOSC => 0.5
#old < dmag4vel4xbnk75Loc11295t1     Transcription factor Sox-2, putative    70%,320/455,498 RefID:UniRef50_E0VH23   RepID:E0VH23_PEDHC
#new > dmag4vel4xbnk75Loc11295t1     SOX transcription factor, putative      100%,111/109,498        RefID:UniRef50_B7PTB7       RepID:B7PTB7_IXOSC

  my($bestid2,$bestaln,$bestpaln,$bestgood,$bestcomm,$bdone)=(0) x 10;
  foreach my $id (sort { $balign{$b} <=> $balign{$a} } keys %kwd) # balign or palign ?
  {
    my $balign= $balign{$id} || 1;
    my $palign= $palign{$id} || $MIN_NAMEIDENT;
    my $isgood= ($id =~ /$goodannots/)?2:1; # should have good/bad-weight
    $isgood=0 if($id =~ /$poorannots/);
    # isgood change to >1 .. 0 .. < -1; or >1 for good, 0.01..0.99 for bad, 1 for no opinion
    if($id =~ m/($sppweights)/) { $isgood= $sppwt{$1} || $isgood; }
    my $iscd=($id=~/CDD:/)?1:0;
    
    my @kw= sort keys %{$kwd{$id}};
    my ($tuniq,$tcomm)=(0,0);
    foreach my $kw (@kw) { if($kwcomm{$kw} > 1) { $tcomm++; } else { $tuniq++; } }
    
    my $bok=0;
    if($tcomm > 0 and $tuniq == 0) { $bok=1; }
    elsif($isgood > 1 and $tcomm > $tuniq) { $bok=1; }
    elsif($tcomm - $tuniq > $bestcomm) { $bok=1; }
    elsif($bestaln < 1 and $palign > $MIN_NAMEIDENT) { $bok=1; }
    if($bok and $bestid2) { # have 1, check more things
      if($palign < $MIN_NAMEIDENT) { $bok=0; }
      elsif($iscd and not($bestid2=~/CDD:/)) { $bok=0; } #? problems w/ CDD now replacing better
      elsif($palign > 2*$bestpaln) { } # keep bok ?
      elsif($bestgood >= $isgood) { $bok=0; }
      elsif($bestgood < 1) { } # keep bok ?
      elsif($isgood > 1 and ($bestcomm - 1 > $tcomm - $tuniq)) { $bok=0; }
      elsif($bestaln > $balign and $palign > 1.8*$bestpaln ) { $bok=0; }
      elsif($bestaln > $balign and $bestpaln > 0.8*$palign and ($bestcomm >= $tcomm - $tuniq)) { $bok=0; }
      $bdone++ unless($bok);
    }
    if($bok) { 
      $bestid2=$id; $bestaln=$balign; $bestpaln= $palign; $bestgood=$isgood; $bestcomm=$tcomm-$tuniq; 
    } else {
      last if($bdone > 1);
    }
  }
  
  $bestid2 ||= $firstid;
  return $inputname{$bestid2}; # always valid? NO..
  #.................
}
  
#   my @kws= sort{ $kw{$b} <=> $kw{$a} } keys %kw;
#  
#   foreach my $id (keys %kwd) {
#     my $wt= $palign{$id} || 1;
#     foreach my $w (@kws) { my $kw=$kw{$w}; if($kwd{$id}{$w}) { $wt += $kw; } else { $wt -= $kw; } }
#     $wt = $wt * 1.2 if($id =~ /$goodannots/); # down/up weight others? : add more good annots
#     $wtd{$id}= $wt;
#   }
#   
#   #DEBUG test pick by top align .. but w/ keywords
#   my($bestid1)= sort{ $wtd{$b} <=> $wtd{$a} or $palign{$b} <=> $palign{$a} } keys %wtd;
#   # my($bestid1)= sort{ $palign{$b} <=> $palign{$a} } keys %kwd; 
#   $bestid1 ||= $firstid;
#   return $inputname{$bestid1}; # always valid? NO..
#   #.................
 
#  ## ? DROP this more complex term weighting?  not doing well, though okay used in omcl/arpconsensus
#   foreach my $id (sort{$palign{$b} <=> $palign{$a}} keys %kwd) { #not sort by id .. useless
#     #skip# my $spp= ($id =~ m/^([a-zA-Z]+)/) ? $1 : $id; # skip this for now
#     #or#   my $spp= ($id =~ m/_([A-Z]+)$/) ? $1 : $id;
#     my $palign= $palign{$id} || 1;
#     my $wt= $palign; # pa == 100% 
#     
#   if(1) {
#     my $wi= @kws;
#     my $nw= int( 0.7 * $wi);
#     foreach my $w (@kws[0..$nw]) { 
#       if($kwd{$id}{$w}) { $wt += $kw{$w}; }
#     }
#   } else {
#     my $wi= @kws;
#     my $nw= int( 0.7 * $wi);
#     ## this now is wting many words, low align >> few words, high align .. bad choice
#     ## try giving each word sum wt of palign{id}?
#     $wi *= 3; # wt by frequency rank order, but not by word frequency?
#     foreach my $w (@kws[0..$nw]) { 
#       # my $kwn= $kw{$w}; # use word count wt, or unique words better?
#       $wt += $wi if( $kwd{$id}{$w} ); 
#       $wi -= 3;
#       }
#     $wt = $wt * $palign/100; # is this right? or too much weight?
#   }
#        
#     #skip# $wt = $wt * 0.3 if($didspp{$spp}++);  
#     # $wt = $wt * 0.4 if($id =~ /$poorannots/); # down/up weight others?
#     $wt = $wt * 1.2 if($id =~ /$goodannots/); # down/up weight others? : add more good annots
#     $wtd{$id}= $wt;
#   }
#   
#   # should require at least 2 cases of agreement, otherwise decline consensus desc
#   my($bestid)= sort{ $wtd{$b} <=> $wtd{$a} or $palign{$b} <=> $palign{$a} } keys %wtd;
#   return $inputname{$bestid}; # always valid?
#   # return($bestid, $de{$bestid});





my (%casewords, %casefirst);

sub recase {
  my($na,$namehash)= @_;  # wont always work, case is context dependent
  if(ref $namehash) {
    #my $D=($DEBUG)?'#':"";
    foreach my $n (sort keys %$namehash) {
      my @w= grep /[A-Z]/, split(/\W+/, $n); 
      if($n =~ /^[A-Z][a-z]+/) { my $w= shift @w; $casefirst{ lc($w) }= $w; }
      map{ $casewords{ lc($_) }= $_ unless(/^Protein$/i); } @w;  # $D.
      
      #? COIL-COILED PROTEIN problem
    }
  }
  
  my @uw; my $i=0;
  my @w= split( /(\W+)/, $na);
  foreach my $w (@w) { 
    my $uw= $casewords{ lc($w) };
    unless($uw or $i>0) { $uw= $casefirst{ lc($w) }; }   
    push(@uw, $uw || $w);  $i++;
    }
  return join("",@uw);
}



# add sort/count names, removing extras: -like putative homolog
sub namecount {
  my($naflags,$recase)=@_;
  my $withID= ($naflags =~ /withID/i)?1:0;
  
  my (%nac, %nav, %naid, %ncase, $n, $id, $src);
  while(<>) {  
    $src= (s/\s+(src[=:]\S+)//) ? $1 : "";
    s/\s*\(\d+[CI%]*\)$//; # namepct
    s/\s*\[[^\]]+\]\s*$//; # ncbi [species]
    
    #old# ($_,$id) = nameclean($_,1); ## old hasid flag
    ($id,$_) = idnameclean($_);
    
    s/(homolog|putative|probable|partial|fragment)\s*//g;  
    s/-like//g;
    # maybe drop trailing nums: family form num:  protein 1|2|3...
    s/\s+\d+$//;
    s/\s*\([^\)]+.//g;  s/\s*\[[^\]]+\]\s*$//;
    s/[,;.]*\s*$//; s/^\s+//;
    
    my $iste= isTEname($_); # FIXME: from input table
    my $lna= lc($_);  
    if($withID) {
      $lna= $id || $lna;
      $lna = $lna."_zzte" if($iste);
      s/$/; $src/ if($src); # preserve for withID
      $nav{$lna}= $_;
    } else {
      $lna = "zzte_".$lna if($iste);
      $nav{$lna}= $_; $nac{$lna} ++;  
    }
    $n++;
    $ncase{$_}++;
  }
  
  if($recase) { 
    recase("xxx", \%ncase);
    while( my($k,$v) = each(%nav) ) { $nav{$k}= recase($v); }
  }
  
  if($withID) {
    foreach my $id (sort keys %nav) { 
      my $nav=$nav{$id};
      my $tep= ($id =~ s/_zzte$//) ? "TE:" : "";
      print  $id, "\t", $tep.$nav, "\n";
    }
  } else {
    foreach my $na (sort keys %nac) { 
      my $tep= ($na =~ /^zzte_/) ? "TE:" : "";
      print  $nac{$na}, "\t", $tep.$nav{$na}, "\n";
    }
  }
  
}


sub isTEname {
  my $na= shift;
  return 0 unless($USE_TENAME);
  ## FIXME: this is notTE: Telomerase reverse transcriptase
  return 0 if($na =~ m/Telomerase/i);
  # return 1 if(/^TE:/); ## add check ^TE: ??
  foreach my $te (@TEnames) { return 1 if($na =~ m/$te/i); }
  return 0;
}

sub getTEnames
{
## likely transposon gene
## see also getTEnames_Aphid() ..
# *** Integrase is not always TE:  Thecc1EG000844t1 Integrase-type DNA-binding superfamily 
# .. MuDR family transposase not always TE gene.
# .. TEname class should be weaker, express>33 overrides; other evid? ..

my @TEnames= map{ s/^\d+\s*//; $_; } grep /^\d/, split "\n", <<'EOTE';
999 Retrotransposon
999 transposon
999 Gag.pro
999 Gag.Pol
999 RNA.directed DNA polymerase
999 Polyprotein
999 Transposase
999 Transposable element
999 reverse transcriptase
998 retropepsin
997 (\b_)Copia(\b_)
997 (\b_)Gypsy(\b_)
EOTE

## !! This bad regex: [\b_]Gypsy[\b_]; should be (\b|_)Gypsy(\b|_)
# not? 999 \bIntegrase\b
## 2013jul additions: esp CDD TE names
## long terminal repeat \bLTR\b RT_nLTR PiggyBac
## CDD: retropepsin_like, Retropepsins 
## CDD: pfam00665, rve, Integrase core domain << is this TE domain? with 'Pol polyprotein'

# warn "# TEnames: @TEnames[0..9,-1] \n" if($DEBUG);
return @TEnames;
}

BEGIN {
  @TEnames=  getTEnames(); # drop this? ; use annot table
}


1;

__END__


=item ncbi tbl2asn whines ; take 2

grep -h '^Change'  okayset/litova1all3.mrna_tsasubmit/*.fixedproducts | sed 's/ for CDS.*//;' | sort | uniq -c
1 Changed 'CDD: SLC5sbd_vSGLT, Vibrio parahaemolyticus Na(+)/galactose-like fragment' to 'CDD: SLC5sbd_vSGLT, Vibrio parahemolyticus Na(+)/galactose-like fragment'
1 Changed 'CDD: SUI1_rel, Translation initation factor SUI1-like' to 'CDD: SUI1_rel, Translation initiation factor SUI1-like'
1 Changed 'CDD: sulP, High affinity sulphate transporter 1-like fragment' to 'CDD: sulP, High affinity sulfate transporter 1-like fragment'
1 Changed 'CDD: UDG_F2_MUG, G:T/U mismatch specific DNA glcosylase fragment' to 'CDD: UDG_F2_MUG, G:T/U mismatch specific DNA glycosylase fragment'
# ^^ not updated nameclean

3 Changed 'Disulphide isomerase fragment' to 'disulfide isomerase fragment'
1 Changed 'Disulphide isomerase' to 'disulfide isomerase'
   #^ britglish bug?
1 Changed 'Potential phospholipid-transporting ATPase IB-like protein fragment' to 'putative phospholipid-transporting ATPase IB-like protein fragment'
1 Changed 'Potential phospholipid-transporting ATPase IB-like protein' to 'putative phospholipid-transporting ATPase IB-like protein'
2 Changed 'Probable ATP-dependent RNA helicase DHX36-like protein' to 'putative ATP-dependent RNA helicase DHX36-like protein'
    # ^^ ? bug here should have cut Probable, Potential ..    
4 Changed 'Uncharacterized HI_0882 protein fragment' to 'putative HI_0882 protein fragment'
2 Changed 'Uncharacterized pp10122-like protein fragment' to 'putative pp10122-like protein fragment'
3 Changed 'Uncharacterized RFS2-like protein fragment' to 'putative RFS2-like protein fragment'
3 Changed 'Uncharacterized RFS2-like protein' to 'putative RFS2-like protein'
    # ^^ ? disagree on Unchar vs puta..  that 'HI_0882/pp10122/RFS' is likely problem for Uncharacterized
    
=cut
