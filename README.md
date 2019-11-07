# TransPi

### Requirenments  
- Directories

  reads = with paired-end reads (e.g. IndA_R1.fastq.gz, IndA_R2.fastq.gz).  
          Make sure reads end with _R2.fastq.gz  
          
  uniprot_db = custom database from Uniprot (e.g metazoans, corals, etc.). Do not include here swiss prot. Swiss prot will be         install in anpther directory.

### STEPS (Palmuc)
1- Run the precheck  
2- Modify the `nextflow.config` for the desire kmers, read length, and PATH for signalP, tmhmm and rnammer.   
3- Run the pipeline


### STEPS (VM - in progress)
1- Run the precheck  
2- Modify the `nextflow.config` for the desire kmers, read length, and PATH for signalP, tmhmm and rnammer.   
3- Run the pipeline

