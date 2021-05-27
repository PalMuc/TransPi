#!/usr/bin/env perl
# ncbigff2evg.pl

# kfish ncbi.gff to usable genes gff
# gunzip -c GFF/ref_Fundulus_heteroclitus-3.0.2_scaffolds.gff3.gz | perl ncbigff2evb.pl > funhe302ncbi.gff

$addchr="chr";

use constant { cRNA => 1, cGENE => 2, cEXON => 3, cCDS => 4, cOTHER => 5 };
while(<>) {
if(/^#Assembly/) { 
  next;  
} elsif(/^#/) { 
  print;  
} elsif(/^\w/) { 
  @v=split"\t"; 
  ($nws,$nsr,$nt,$nb,$ne,$nv,$no,$nx,$nat)=@v; 
  if($nt eq "region") {
    if(my($chr)=m/chromosome=(\w+)/) { 
      if($chr eq "Unknown" or not m/\bID=/) { }
      else { if($chr=~/^\d+/) { $chr=$addchr.$chr; } 
      $nsc{$nws}=$chr; s/ID=\w+/ID=$nws/;  }
    }
    # next;
  }
  $isrna=($nt =~ m/RNA|transcript/)?1:0;
  $ft= ($isrna)? cRNA : ($nt =~ /^gene|pseudogene/)? cGENE 
       :($nt=~/^exon/)? cEXON :($nt=~/^CDS/)?cCDS : cOTHER; 
  $nsc=$nsc{$nws}||$nws; s/^$nws/$nsc/; 
  s/;(gbkey)=[^;\n]+//g; 
  ($lid)=m/ID=(\w+)/; ($lpd)=m/Parent=(\w+)/; 

     if($ft == cRNA) { ($td)=m/Genbank:([\w\.]+)/; $tld{$lid}=$td; s/;(model_evidence)=[^;\n]+//; } 
  elsif($ft == cEXON) { ($td)=m/Genbank:([\w\.]+)/; s/;(product|gene|exception|partial|inference|Dbxref|Note)=[^;\n]+//g; } 
  elsif($ft == cCDS) { $td=$tld{$lpd}||$lpd;  s/;(product|gene|Name|exception|partial|inference|Dbxref|Note)=[^;\n]+//g; } 
  elsif($ft == cGENE) { $td=""; } 
  if($td) { 
     if($ft == cRNA) { $nd="ID=$td" } else { $nd="Parent=$td"; }  
     s/Parent=/lpar=/; s/ID=/$nd;lid=/; } 
  print; 
}

}

=item ncbi.gff parts

NC_010444.4     RefSeq  region  1       151935994       .       +       .       ID=id91877;Dbxref=taxon:9823;Name=2;breed=Duroc;chromosome=2;gbkey=Src;genome=chromosome;isolate=TJ Tabasco;mol_type=genomic DNA;sex=female

drop attr
 exon exception partial transcript_id 
  cds Name exception inference partia

types:
4219 CDS
5133 exon
 176 gene
  83 lnc_RNA
 337 mRNA
   5 miRNA
   5 primary_transcript
   7 pseudogene
   1 region
  33 transcript

=cut

