OrthoMCL: Ortholog Clustering Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
AUTHOR: Feng Chen <fengchen@sas.upenn.edu>, Li Li
ORTHOMCL [2006-10-02] Version 1.4 + MCL-02-063

Copyright (C) 2004~2006 by University of Pennsylvania, Philadelphia, PA 
USA. All rights reserved.

Before orthomcl.pl can be used, some variables (including directory 
variables or parameter variables) in orthomcl_module.pm need to be set, 
as described in README.

Reference: 
1. Li Li, Christian J. Stoeckert, Jr. and David S. Roos. OrthoMCL: 
Identification of Ortholog Groups for Eukaryotic Genomes. Genome Research
13:2178-2189, 2003
http://www.genome.org/cgi/content/full/13/9/2178
2. Feng Chen, Aaron J. Mackey, Christian J. Stoeckert, Jr., and David S.
Roos OrthoMCL-DB: querying a comprehensive multi-species collection of 
ortholog groups. Nucleic Acids Res. 2006 34: D363-8.
http://nar.oxfordjournals.org/cgi/content/full/34/suppl_1/D363

Preconditions: A Unix System, Perl and some experience with Unix & Perl

1. Installation of required softwares and perl modules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
OrthoMCL is a perl script which doesn't need compilation. However, it 
requires some softwares and perl modules to run, as listed below:

Softwares:      1. BLAST (NCBI-BLAST, WU-BLAST, etc.)
               *2. MCL (Markov Clustering algorithm), available 
                   at http://micans.org/mcl/. NOTE: MCL changed the output 
                   format recently which is not compatible with OrthoMCL. 
                   Please use the MCL-02-063 version enclosed with this 
                   package, which has been the default for all test analysis.
Perl Modules:   1. Bio::SearchIO (part of BioPerl, http://bioperl.org)
                2. Storable


2. Setting the variables in "OrthoMCL/orthomcl_module.pm"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Most global variables used in "OrthoMCL/orthomcl.pl" need to be set
in the perl module "OrthoMCL/orthomcl_module.pm".

--- $PATH_TO_ORTHOMCL: the orthomcl directory itself
    (example: $PATH_TO_ORTHOMCL = "/disk3/fengchen/orthomcl/";)
--- $BLASTALL: your BLAST software
    (example: $BLASTALL  = "/genomics/share/bin/blastall";)
--- $BLAST_FORMAT: how BLAST result is stored
    Options: a) "compact" corresponds to NCBI-BLAST's -m 8
             b) "full" corresponds to NCBI-BLAST's -m 0
    # for WU-BLAST, make changes on subroutine executeBLASTALL
--- $BLAST_NOCPU: the number of CPUs
    For multi-processor machine, setting it higher than 1
    will significantly save time in BLAST step
--- $FORMATDB: your FORMATDB software
    (example: $FORMATDB = "/genomics/share/bin/formatdb";)
--- $MCL: your MCL software
    (example: $MCL = "/disk2/fengchen/mcl-02-063/shmcl/mcl";)
--- $MAX_WEIGHT_DEFAULT: Weight used for protein pairs whose BLAST p-value is
    zero (0). This depends on the algorithm you use: if the second smallest 
    p-value is in the order of -99, maximum_weight should be 100; if -299, 
    maximum_weight should be 300 <DEFAULT>.

Now you can run orthomcl.pl on a three-species test set, since
the variables $PATH_TO_ORTHOMCL, $BLASTALL, $FORMATDB and $MCL
are set.
 % orthomcl.pl --mode 1 --fa_files "Ath.fa,Hsa.fa,Sce.fa"

 Note: Here the test set Ath.fa, Hsa.fa and Sce.fa only contain
 15, 16 and 11 sequences, respectively. Such a test set is selected
 to make sure you have everything set and OrthoMCL can run on your
 machine. Since it takes OrthoMCL long time to finish clustering
 a big data set, from BLAST to MCL, it's wise to try a very small
 set first.
 
To use OrthoMCL on your data, you need to collect protein fasta
files ".fa" (with each ".fa" file representing one species only, and
having a simple name, e.g. "Eco.fa") and put them in the directory 
"data" or reset the following variable:

--- $ORTHOMCL_DATA_DIR: the data directory to store the fasta files
    $ORTHOMCL_DATA_DIR = $PATH_TO_ORTHOMCL."/data/"; (DEFAULT)


3. Running OrthoMCL
~~~~~~~~~~~~~~~~~~~

The COMPLETE proteome data for each species should be chosen, theoretically.
And you should have enough memory (>=800MB) if you have around 100,000
sequences to cluster, because this stand-alone version tries to read 
BLAST information into memory.

There are five modes to run OrthoMCL, with each mode having a different
process. We strongly suggest you to use MODE 4 for very big set, since
BLAST was not programmed to run parallelly. You can simply prepare two
files for mode 4, BPO file and GG file. And it's very fast, for our test
set of 200,000 sequences on a Mac G5 computer, it took 8 hours to finish.

The five modes of OrthoMCL are:

  Mode 1: OrthoMCL analysis from FASTA files. OrthoMCL starts from
          the beginning BLAST to final MCL.

     Example:  % orthomcl.pl --mode 1 --fa_files Ath.fa,Hsa.fa,Sce.fa

  Mode 2: OrthoMCL analysis based on former OrthoMCL run (former run
          directory needs to be given), if you want to change the 
          inflation parameter, p-value cutoff (can only be lower than
          your former run BLAST p-value cutoff), percent identity cutoff
          or percent match cutoff. No BLAST or BLAST parsing performed.

     Example:  % orthomcl.pl --mode 2 --former_run_dir Sep_8 --inflation 1.4

  Mode 3: OrthoMCL analysis from user-provided BLAST result BLAST out file
          and genome gene relation file telling which genome has which gene
          (Please refer to 5. File Formats). No BLAST performed.

     Example:  % orthomcl.pl --mode 3 --blast_file AtCeHs_blast.out --gg_file 
                 AtCeHs.gg


  Mode 4: OrthoMCL analysis from user-provided BPO (BLAST PARSING OUT) file
          and GG (genome gene relation) file telling which genome has which gene
          (Please refer to 5. File Formats). No BLAST or BLAST parsing performed.

     Example:  % orthomcl.pl --mode 4 --bpo_file AtCeHs.bpo --gg_file AtCeHs.gg

  Mode 5: OrthoMCL analysis based on previous run, but with less taxa included
          or with only inflation value changed (FASTER than mode 2, no selection
          on reciprocal best/better hits performed).

     Example:  % orthomcl.pl --mode 5 --former_run_dir Sep_8 --taxa_file AtCeHs.gg
                 --inflation=1.1

4. OrthoMCL Arguments
~~~~~~~~~~~~~~~~~~~~~

 fa_files=<String>       Protein FASTA file names, with each file 
                         containing protein sequences from one species,
                         separated by comma(e.g. "Eco.fa,Sce.fa,Afu.fa")
 pv_cutoff=<Float>       P-Value or E-Value Cutoff in BLAST search and/or
                         ortholog clustering, 1e-5 (DEFAULT).
 pi_cutoff=<Int>         Percent Identity Cutoff <0-100> in ortholog 
                         clustering, 0 (DEFAULT).
 pmatch_cutoff=<Int>     Percent Match Cutoff <0-100> in ortholog
                         clustering, 0 (DEFAULT).
 inflation=<Float>       Markov Inflation Index, used in MCL algorithm,
                         1.5 (DEFAULT). Increasing this index increases
                         cluster tightness, and the number of clusters.
 former_run_dir=<String> Former run directory, required in Mode 2, e.g. 
                         "July_21". Then the blast result file and bpo 
                         file in former run directory will be used 
                         instead of running from the very beginning.
 blast_file=<String>     Blast out file provided by user, required in
                         Mode 3. It will be parsed into BPO file, and
                         further used for ortholog clustering.
 bpo_file=<String>       BPO (Blast Parse Out) file provided by user,
                         required in Mode 4. Please refer to 
                         ORTHOMCL_INSTALL about its format.
 gg_file=<String>        GG (Genome Gene mapping) file provided by user,
                         required in Mode 3 & 4. Please refer to 
                         ORTHOMCL_INSTALL about its format.
 taxa_file=<String>      TAXA file provided by user, required in Mode 5. 
                         Please refer to ORTHOMCL_INSTALL about its
                         format.

5. File Format in OrthoMCL
~~~~~~~~~~~~~~~~~~~~~~~~~~

There are some files generated during OrthoMCL running:
   all_orthomcl.out          Final result showing how genes are clustered
   all_blast.bbh             List of reciprocal best/better hit, for further reference
   parameter.log             Parameter Log file recording all parameters
   orthomcl.log              OrthoMCL Log file
   tmp/all.fa                Fasta file containing protein sequences from all species
   tmp/all.gg                Genome gene relation file, telling which genome contains
                             which gene, format is shown afterwards
   tmp/all_blast.out         All-against-all BLAST out file
   tmp/all.bpo               BLAST PARSING OUT file from all_blast.out
   tmp/all_bpo.idx           Storing a data array, used in blast query
   tmp/all_bpo.se            Storing a hash, used in blast query
   tmp/all_ortho.mtx         Matrix file storing all weights
   tmp/all_ortho.idx         Mapping information between Gene ID's and Index'es used in Matrix
   tmp/all_ortho.mcl         MCL clustering result

FORMAT of all.gg or "usr_gg_file"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ath: At1g01190 At1g01280 At1g04160 ...
Hsa: Hs10834998 Hs10835119 Hs10835271 ...
Sce: YAL029c YAR009c YAR010c YHR023w ...

Each line stands for each genome. Each line starts with genome name, followed by a 
colon ":", and then followed by all the gene id's separated by space key " ".

FORMAT of all.bpo or "usr_bpo_file"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1;At1g01190;535;At1g01190;535;0.0;97;1:1-535:1-535.
2;At1g01190;535;At1g01280;510;2e-56;29;1:69-499:28-474.
3;At1g01190;535;At1g11600;510;1e-45;27;1:59-531:21-509.

Each line represents each query-subject similarity relation. And all the info is
separated by ";", which are, in order, similarity id, query id, query length, 
subject id, subject length, BLAST E-value, percent identity, HSP info (each HSP
is in the format of HSP_id:query_start-query_end:subject_start-subject_end. 
different HSP info are seperated by "." )
IMPORTANT: 1. Similarity ID represents BPO file line id, so it should start 
              from 1 for the first line, and be consecutive for the whole file.
           2. BPO file is a parsing result from BLAST, so for each query gene
              id, its hits can't be scattered in the file, but should be listed 
              in ajacent lines.
           3. For BLAST m8 format (i.e. $BLAST_FORMAT="compact"), sequence length
              information is not stored. So when running OrthoMCL in mode 3,
              the corresponding columns (i.e. query length, and subject length)
              will be 0. Please do not use Percent match cutoff in this case.

FORMAT of user taxa list file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This file is simple: you can use similar file like GG file, or use a file with
each taxa id listed in different lines. (Each taxa id should be defined in GG file)


I attached sample GG and BPO files in directory "sample_data".

Please report any bug to us. Thanks a lot.

**********************BLAST result truncation problem***********************
I got several emails reporting "division by zero" error during BLAST parsing
process. If you happen to meet this as well, that means some of your BLAST 
result is truncated, which causes some hits are not listed in the alignment
part. In this case, set $blast_flag{'hsp'}=0 in orthomcl.pl lines 48-60,
which will not result in HSP parsing. Besides, you can not use any percent
identity cutoff or percent coverage cutoff, since your data lack such 
information.
****************************************************************************

Feng Chen fengchen@sas.upenn.edu

Roos Lab
Depart of Biology
University of Pennsylvania