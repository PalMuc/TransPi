# TransPi

## REQUIREMENTS   
- Directories **before** the precheck is run  

  reads = with paired-end reads (e.g. IndA_R1.fastq.gz, IndA_R2.fastq.gz).    
  		**Make sure reads end with _R1.fastq.gz and _R1.fastq.gz**  
          
  uniprot_db = custom database from Uniprot (e.g metazoans, corals, etc.). Do not include here swiss prot. Swiss prot will be         install in another directory.  
		**Example: uniprot-taxonomy_metazoaA33208.fasta**   


## STEPS before running the pipeline    
1- After the requirentments are satisfied you can now run the precheck. Precheck run needs a PATH as an argument for running and installing (locally) all the databases and programs the pipieline needs.   

```
bash precheck_TransPi.sh /home/ubuntu/TransPi
```  

2- If the precheck run was successful, a file called `nextflow.config` will be created. Modify in the `nextflow.config` the kmer list for your reads and read lengths. Example:    

```  
    // kmers list (depends on read length!)
    k="25,37,45,53"

    //#maximal read length
    max_rd_len="76"

``` 

3-  In the same file, `nextflow.config`, edit now the `PATH` of rnammer, tmhmm and siglnalp. Example:     

```
    //rnammer
    rnam="/home/ubuntu/pipe/rnammer/rnammer"
    //tmhmm
    tmhmm="/home/ubuntu/pipe/tmhmm-2.0c/bin/tmhmm"
    //signalP
    signalp="/home/ubuntu/pipe/signalp-4.1/signalp"

```   

The file `nextflow.config` is created using the `template.nextflow.config`. It is convenient then to modify the `PATH` of rnammer, tmhmm and siglnalp in the `template.nextflow.config` so the every time you run the precheck just the PATH for the databases (*e.g.* BUSCO) are changed.     


4- If your installation of conda is not located on `~/anaconda3` modify the `nextflow.config` for the `PATH` of your conda installation. Example:  

```

    //Examples: ~/anaconda3...    ~/tools/anaconda3...   ~/tools/py3/anaconda3...
    condash="~/anaconda3"

```  

## RUNNING the pipeline (VM testing)  

1- To run the pipeline, first activate the conda environment before of TransPi.  

```
conda activate TransPi
```  

2- Run the pipeline.   

```
./nextflow TransPi_part1_test.nf
```

3- If an error occur and you need to resume the run just include the `-resume` when calling the pipeline.  

```
./nextflow TransPi_part1_test.nf -resume

```    
 

## NOTES
1- The precheck run is designed to create a new `nextflow.config` every time is run with with the respectives `PATH` to the databases. You can modify the values that do not need editing for your analysis on the `template.nextflow.config` to avoid doing the changes after the precheck run.  

2- The `template.nextflow.config` file has different configuration for the each program of the pipeline (*e.g.* some with a lot of CPUs, others with a small amount of CPUs). You can modify this depending on the resources you have in your system. Example:

```
process {
    withLabel: big_cpus {
        cpus='30'
        memory='15 GB'
        clusterOptions='-p lemmium --qos=normal'
        executor='slurm'
    }

```

In this case, the processes using the label `big_cpus` will use 30 CPUs. If your system only has 20 please modify this values accordingly to avoid errors. Also, you will notice that we are using `SLURM` as our job manager in our server. If you do not need this specification just simply get rid of the followiing lines that are custom for our system.    

```
	//erase this lines
        clusterOptions='-p lemmium --qos=normal'
        executor='slurm'
```  


3- To avoid calling the pipeline using `./nextflow` you can modify the nextflow command like this `chmod 777 nextflow`. For running the pipeline you just need to use:  

```
nextflow TransPi_part1_test.nf

```    


## OPTIONAL   
Make an account with your email at [Nextflow Tower](https://tower.nf/login) for monitoring the pipeline on any internet browser without connecting to the server.  

The run the pipeline adding the `-with-tower` option.

```
nextflow TransPi_part1_test.nf -with-tower 
```  
