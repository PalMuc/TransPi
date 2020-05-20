# TransPi - TRanscriptome ANalysiS PIpeline (in development)

```
 _______                                 _____   _ 
|__   __|                               |  __ \ (_)
   | |     _ __    __ _   _ __    ___   | |__) | _ 
   | |    |  __|  / _  | |  _ \  / __|  |  ___/ | |
   | |    | |    | (_| | | | | | \__ \  | |     | |
   |_|    |_|     \__,_| |_| |_| |___/  |_|     |_|
 ```

# Table of contents
* [General info](#General-info)
* [Detailed description](#Detailed-description)  
    * [Precheck run](#Precheck-run)  
    * [Main Pipeline](#Main-pipeline)  
* [Why use TransPi?](#Why-use-TransPi)  
* [What TransPi uses for the analyses?](#What-TransPi-uses-for-the-analyses)  
* [How to use TransPi](#How-to-use-TransPi)
    * [Requirements](#Requirements)
    * [Downloading TransPi](#Downloading-TransPi)
    * [Installation and configuration](#Installation-and-configuration)
    * [Running TransPi](#Running-TransPi)
* [Extra notes](#Extra-notes)
* [Future work](#Future-work)

# General info
TransPi is TRanscriptome ANalysiS PIpeline based on the scientific workflow manager [Nextflow](https://www.nextflow.io). It is designed to help researchers get the best reference transcriptome assembly for their organisms of interest. It performs multiple assemblies with different parameters to then get a non-redundant consensus assembly. It also performs other valuable analyses such as quality assessment of the assembly, BUSCO scores, Transdecoder (ORFs), and gene ontologies (Trinotate). All these with minimum input from the user but without losing the potential of a comprehensive analysis.

# Detailed description

### Precheck run (`precheck_TransPi.sh`)

TransPi depends on several programs and databases (e.g. SwissProt, PFAM, etc.) to generate the reference transcriptome and its annotation. For that reason, we designed a simple shell script (i.e. precheck_TransPi.sh) to help ease the setup of all the dependencies currently needed by TransPi to run. The script will essentially check for the necessary directories for running the tool (e.g. reads directory), install the Conda package management system (if needed), install all dependencies, configure databases, automatically download and configure BUSCO databases (Simão et al., 2015, v3 and v4), download Nextflow, among others. The script is designed to recognize when a previous run of the script was done, thus skipping steps that do not need to be repeated. Another advantage of the precheck script is that it will automatically create the configuration file needed by Nextflow to execute the pipeline with all the necessary information . As a result, the user will only have to make some minor changes to the file (e.g. select kmers) before running the pipeline. Essentially, the precheck has to be run entirely only one time for the dependencies and databases installation. Subsequent pipeline runs can be done with the same configuration file.  


### Main pipeline (`TransPi.nf`)
A detailed diagram of the complete TransPi pipeline is presented in **Figure 1**. First, reads are checked for adapter presence and/or errors with Fastqc (Andrews, 2010). Reads are then filtered (Q>25) and adapters are removed (if existent) with fastp (Chen et al., 2018). Filtered reads are then normalized before being assembled using a combination of five different assemblers and kmers. The assemblers used by TransPi are rnaSPADES (Bushmanova et al., 2019), Trans-ABySS (Robertson et al., 2010), SOAP (Luo et al., 2012), Trinity (Grabherr et al., 2011) and Velvet/Oases (Zerbino & Birney, 2008, Schulz et al., 2012). Together they generate an over-assembled transcriptome, which is then reduced with EvidentialGene (Gilbert, 2013, 2019). EvidentialGene aims to keep the most valid biological transcript, discard the less valid and, reduce the redundancy of the multiple assemblers to arrive at the best non-redundant consensus transcriptome assembly. This is done by merging perfect duplicates, clustering of proteins, and local similarities searches with BLAST (Altschul et al., 1997) (for more details see Gilbert, 2019).   

![TransPi flowchart](https://sync.palmuc.org/index.php/s/FNYmpbgtziFkZte/preview)

**Figure 1. TransPi flowchart**   

Next, TransPi uses this non-redundant reference transcriptome to run several downstream analyses commonly applied to de novo transcriptomes projects: (1) rnaQUAST for quality assessment of the transcriptome (Bushmanova et al., 2016). (2) BUSCO (Simao et al., 2015; v3 and v4) to quantitatively assess the completeness in terms of expected gene content of the consensus transcriptome assembly. (3) TransDecoder (https://transdecoder.github.io) to identify open reading frames (ORFs) within the consensus assembly. It scores them according to their sequence composition and retains the ones that are consistent with coding transcripts. (4) Homology searches of all TransDecoder predicted open reading frames (ORFs) to known proteins via BLAST in order to retain ORFs that may have functional significance but don’t pass the coding likelihood scores. (5) Trinotate (Bryant et al., 2017) to provide automatic functional annotation of the reference transcriptome assembly. We make use of Diamond (Buchfink et al., 2015) to accelerate the search similarities to the SwissProt and custom UniProt databases (chosen by the user).  RNAmmer (Lagesen et al., 2007), TMhmm (Krogh et al., 2001), SignalP (Petersen et al., 2011) are used to search for ribosomal RNA, signal peptide proteins and transmembrane domain prediction, respectively. Protein domain searches are done with HMMER (Finn et al., 2011) against the last version of the PFAM database. All this information is combined by Trinotate for the creation of an annotation report. These reports will also include information on Gene Ontology (GO), eggNOG, and KEGG together with the similarity search done to SwissProt and the custom UniProt database.  


# Why use TransPi?
Even though most transcriptome analyses are only based on Trinity, this does not mean is the "*one fits all*" solution for transcriptome assembly. Gene sizes, heterozygosity, read length, kmers, quality, and many others affect the performance of the assemblers. Our method is designed to avoid such problems by using multiple assemblers and different parameters for the assemblies. All combined transcripts are reduced to a reference transcriptome without the redundancy created by over assembly of the reads (i.e. multiple programs and parameters).

Example: Long genes are not always assembled correctly by Trinity. However, Velvet (with long kmers) is more suitable for such genes, thus compensating the missing genes in the transcriptome generated by Trinity.  

We have tested the pipeline with various non-model organisms that often present difficulty in these analyses. **More info here**    

**INFO here about the test we made and the improvement**  


# What TransPi uses for the analyses?  

List of programs use by TransPi:
- rnaSPADES  
- SOAP  
- Trinity  
- Velvet  
- EvidentialGene  
- CD-Hit  
- BUSCO  
- Trandecoder  
- Trinotate  
- Diamond  
- SQLite  
- Hmmer  
- Exonerate  
- Blast  
- Bowtie2  
- rnammer
- tmhmm
- signalP
- R  
- Python  
- Java  

Databases used by TransPi:
- Swissprot
- Uniprot custom database (e.g. all metazoan proteins)
- Pfam

Installing all these programs can be confusing and time consuming. Thus, apart from the automatization for the transcriptome analysis by TransPi, we also developed a script to install all dependencies and databases needed for the pipeline. Also, it will create all the necessary config files needed to run the analysis. More info on the "How to use TransPi" section below.  

**NOTE**  
Rnammer, tmhmm and signalP are for academic users; other users are requested to contact CBS Software Package Manager at software@cbs.dtu.dk.  

For internal use of **PALMUC** users this is provided in the repository. No need to get a new license.  

# How to use TransPi

## Requirements   

- System: Linux OS  

- Data type: Paired-end reads   

- Directory called `reads` with the compressed FASTQ files to analyze  

    Create this directory before all other steps.   

    Example: IndA_R1.fastq.gz, IndA_R2.fastq.gz

    NOTES:  
    - Make sure reads end with `_R1.fastq.gz` and `_R2.fastq.gz`    

    - Multiple individuals can be added to the same directory.  


## Downloading TransPi   

1- Clone the repository   

```

git clone https://github.com/rivera10/TransPi.git  

```  

2- Move to the TransPi directory  

```

cd TransPi

```  


## Installation and configuration

1- After the requirements are satisfied you can run the precheck script. The precheck run needs a PATH as an argument for running and installing (locally) all the databases and programs the pipeline needs.   

```

bash precheck_TransPi.sh /home/ubuntu/TransPi

```  

**NOTE:**   
This process may take a while depending on the options you select. Step that takes longer is if you choose to download the entire metazoan proteins from UniProt (6Gb). Other processes and databases are relatively fast depending on internet connection.   

2- If the precheck run was successful all the programs and databases for running TransPi were installed (locally) in the specified PATH you ran the script. Proceed to modify the file `nextflow.config` for the kmer list to be used for your reads and your reads length.  

Example for read lengths of 100bp:    

```  
    // kmers list (depends on read length!)
    k="25,37,43,55,63"

    //#maximal read length
    max_rd_len="100"
    //[LIB]
    //#maximal read length in this lib
    rd_len_cutof="100"

```

**NOTES**  
1- If you combined multiple libraries of the same individual to create a reference transcriptome, which will be later use in downstream analyses (e.g. Differential Expression), make sure the kmer list is based on the length for the shortest read library.  

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Example: Combining reads of 100bp with 125bp  

2- There are other values you can modify such as insert size and minimum aligned length. If not sure leave the values as it is.   

3- Lines that start with `//` are comments that can provide information of the parameter that follows the comment.  

## Running TransPi   

#### Complete pipeline
We recommend to run TransPi with the option `--all` where it will do the complete analysis, from raw reads to annotation.


To run the complete pipeline.    
```

./nextflow run TransPi.nf --all -with-conda ~/anaconda3/envs/TransPi


```

In this example the whole pipeline will be run (`--all`) and with are using the conda installation from the precheck (`-with-conda CONDA_FROM_PRECHECK`)  

#### Other options  
Currently, the pipeline has other features depending on the user's needs.    

Mandatory arguments (--all or --onlyAsm or --onlyEvi or --onlyAnn):
```
        --all           Run the entire pipeline (Assemblies, EvidentialGene, Annotation, etc.)

        --onlyAsm       Run only the Assemblies and EvidentialGene analysis

        --onlyEvi       Run only the Evidential Gene analysis

        --onlyAnn       Run only the Annotation analysis (starting from a final assembly)
```
Other options:
```
        -with-conda     To run with a local conda installation (generated with the precheck) and not installed by nextflow
                        This is the preferred method of running TransPi

                        Example:
                            nextflow run transpi.nf --all -with-conda /home/ubuntu/anaconda3/envs/TransPi

        -profile        Configuration profile to use. Can use multiple (comma separated)

                        Available:
                            test (Run TransPi with test dataset)
                            conda (Use a conda env created by nextflow. Not recommended, not all programs are installed by conda. Use the precheck)
                            docker (in development - Run TransPi with a docker container with all the neccesary tools)
                            singularity (in development - Run TransPi with a singularity container with all the neccesary tools)

```
</br>  

## Extra notes

1- If an error occurs and you need to resume the pipeline just include the `-resume` option when calling the pipeline.  

```

./nextflow run TransPi.nf --all -with-conda ~/anaconda3/envs/TransPi -resume


```   

</br>  

2- The `template.nextflow.config` file has different configurations for the each program of the pipeline (*e.g.* some with a lot of CPUs, others with a small amount of CPUs). You can modify this depending on the resources you have in your system.   

Example:

```

process {
    withLabel: big_cpus {
        cpus='30'
        memory='15 GB'
        clusterOptions='-p lemmium --qos=normal'
        executor='slurm'
    }


```  

In this case, all the processes using the label `big_cpus` will use 30 CPUs. If your system only has 20 please modify this values accordingly to avoid errors. Also, you will notice that we are using [SLURM](https://slurm.schedmd.com/documentation.html) as our job manager in our server. If you do not need this specification just simply erase of the following lines that are custom for our system.    

```
	//erase this lines
        clusterOptions='-p lemmium --qos=normal'
        executor='slurm'
```  
In the contrary, if your system is using `SLURM`, `SGE`, `PBS/Torque`, etc., change the `executor` part for the one use by your system.

Example:
```

process {
    withLabel: big_cpus {
        cpus='30'
        memory='15 GB'
        clusterOptions='-p lemmium --qos=normal'
        executor='pbs'
    }


```


The line `clusterOptions` can be used to add any other option that you will usually use for your job submission.  

Example of `SLURM`:

```

#SBATCH --get-user-env                                                                                                           
#SBATCH --clusters=inter                                                                                                  
#SBATCH --partition=teramem_inter                                                                                              
#SBATCH --mail-user=user@mail.com                                                                                            
#SBATCH --mail-type=ALL

```
Can be used like this:

```

clusterOptions='--get-user-env --clusters=inter --partition=teramem_inter --mail-user=user@mail.com --mail-type=ALL'


```

If you run this on a local or virtual machine (VM) with no need for a job scheduler you can erase the lines `clusterOptions` and `executor`.  

Example:
```

process {
    withLabel: big_cpus {
        cpus='30'
        memory='15 GB'
    }


```

This way nextflow will know the resources for each process with that specific label. This is really important, because nextflow will start as many jobs as possible if the resources are available. If you are in a VM with 120 CPUs, nextfow will be able to start four processes with this configuration.   

</br>

3- The precheck run is designed to create a new `nextflow.config` every time is run with with the all the `PATH` to the databases. You can modify the values that do not need editing for your analysis on the `template.nextflow.config` to avoid doing the changes after the precheck run.

Example: Modify the `template.nextflow.config` with your cluster info to avoid issues in the future.  

</br>

4- To avoid calling the pipeline using `./nextflow` you can modify the nextflow command like this `chmod 777 nextflow`. For running the pipeline you just need to use:  

```

nextflow run TransPi.nf

```    

5- To monitor your pipeline remotely without connecting to the server via ssh use [Nextflow Tower](https://tower.nf/login). Make an account with your email and follow their instructions. After this, you can now run the pipeline adding the `-with-tower` option. **Not to use in PALMUC server.**   

```

./nextflow run TransPi.nf --all -with-tower -with-conda ~/anaconda3/envs/TransPi  


```


## Future work
- Docker and Singularity for running the pipeline (in progress)
- Create summary file with all results
- Cloud deployment of the tool
