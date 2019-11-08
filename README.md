# TransPi

## REQUIREMENTS   
- Directories *before* the precheck is run  

  reads = with paired-end reads (e.g. IndA_R1.fastq.gz, IndA_R2.fastq.gz).    
  		**Make sure reads end with _R1.fastq.gz and _R1.fastq.gz**  
          
  uniprot_db = custom database from Uniprot (e.g metazoans, corals, etc.). Do not include here swiss prot. Swiss prot will be         install in another directory.  
		**Example: uniprot-taxonomy_metazoaA33208.fasta**   


## STEPS (VM testing)   
1- After the requirentments are satisfied you can now run the precheck. Precehck run needs a PATH as narguments for running installing (locally) all the databse and programs it needs.  

```bash precheck_TransPi.sh /home/ubuntu/TransPi```  

2- If the precheck run was succesful, a file called `nextflow.config` will be created. You just need to edit now the `PATH` of rnammer, tmhmm and siglnalp. This file is created using the `template.nextflow.config`. It is convenient then to modify the `PATH` of rnammer, tmhmm and siglnalp in the `template.nextflow.config` so the only changes need to be in the kmers and reads length.   

3- To run the pipeline, first activate the conda environment before running TransPi.

```conda activate TransPi```  

4- Run the pipeline.   

```./nextflow TransPi_part1_test.nf```

5- If an error occur and you need to resume the run just include the `-resume` when calling the pipeline.  

```./nextflow TransPi_part1_test.nf -resume```  

Possible errors:  

1- The pipeline is written to use the installation of conda stored at `~/anaconda3`. If your installation differs from this `PATH` you need to modify the pipeline like this: `sed -i "s|/OLD/PATH/HERE|/NEW/PATH/HERE|g" TransPi_part1_test.nf`    

Example: 
```sed -i "s|\~/anaconda3|\~/tools/python3/anaconda3|g" TransPi_part1_test.nf```  
  

## NOTES
The precheck run is designed to create a new `nextflow.config` file after each run with the respectives `PATH`.   


## OPTIONAL   
Make an account with your email at [Nextflow Tower](https://tower.nf/login) for monitoring the pipeline on any internet browser without connecting to the server.  
