#!/usr/bin/env nextflow
/*
========================================================================================
                                    TransPi
========================================================================================
                       Transcriptome Analysis Pipeline
                       Author: Ramón E. Rivera-Vicéns
                       GitHub: rivera10
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info """
    ==================================================
      TransPi - Transcriptome Analysis Pipeline v${workflow.manifest.version}
    ==================================================

    Steps:
        1- Run the `precheck_TransPi.sh` to set up the databases and tools (if neccesary) used by TransPi
        2- Run TransPi

        Usage:
            nextflow run TransPi.nf TransPi_analysis_option other_options

            Example usage:
                nextflow run TransPi.nf --all --reads "/PATH/TO/READS/*_R[1,2].fastq.gz" --k 25,41,53 --maxReadLen 75

        Manadatory arguments:
            --all                   Run the entire pipeline (Assemblies, EvidentialGene, Annotation, etc.)
                                    This option also requires arguments --reads, --k, and --maxReadLen
                                    Example:
                                        --reads "/PATH/TO/READS/*_R[1,2].fastq.gz" --k 25,35,55,75,85 --maxReadLen 150
                                        NOTE: Use of quotes is needed for the reads PATH. Kmer list depends on read length.

            --onlyAsm               Run only the Assemblies and EvidentialGene analysis
                                    This option also requires arguments --reads, --k, --maxReadLen

            --onlyEvi               Run only the Evidential Gene analysis
                                    Transcriptome expected to be in a directory called "onlyEvi"

            --onlyAnn               Run only the Annotation analysis (starting from a final assembly)
                                    Transcriptome expected to be in a directory called "onlyAnn"

        Other options:
            -profile                Configuration profile to use. Can use multiple (comma separated)
                    test                Run TransPi with a test dataset
                    conda               Run TransPi with conda.
                    docker              Run TransPi with docker container
                    singularity         Run TransPi with singularity container with all the neccesary tools
                    TransPiContainer    Run TransPi with a single container with all tools

            --help                  Display this message
            --fullHelp              Display this message and examples for running TransPi

        Output options:
            --outdir 	            Name of output directory. Default "results"
            -w, -work 	            Name of working directory. Default "work". Only one dash is needed for -work since it is a nextflow function.
            --tracedir              Name for directory to save pipeline trace files. Default "pipeline_info"

        Additional analyses:
            --rRNAfilter            Remove rRNA from sequences. Requires option --rRNAdb
                --rRNAdb                PATH to database of rRNA sequences to use for filtering of rRNA. Default ""

            --filterSpecies         Perform psytrans filtering of transcriptome. Default "false" Requires options --host and --symbiont
                --host 	                PATH to host (or similar) protein file. Default ""
                --symbiont     	        PATH to symbionts (or similar) protein files. Default ""

            --psyval       	        Psytrans value to train model. Default "160"
            --allBuscos       	    Run BUSCO analysis in all assemblies. Default "false"
            --rescueBusco 	        Generate BUSCO distribution analysis. Default "false"
            --minPerc        	    Mininmum percentage of assemblers require for the BUSCO distribution. Default ".70"
            --shortTransdecoder     Run Transdecoder without the homology searches. Default "false"
            --withSignalP 	        Include SignalP for the annotation. Needs manual installation of CBS-DTU tools. Default "false". Requires --signalp
                --signalp               PATH to SignalP software. Default ""

            --withTMHMM 	        Include TMHMM for the annotation. Needs manual installation of CBS-DTU tools. Default "false". Requires --tmhmm
                --tmhmm                 PATH to TMHMM software. Default ""

            --withRnammer 	        Include Rnammer for the annotation. Needs manual installation of CBS-DTU tools. Default "false". Requires --rnam
                --rnam                  PATH to Rnammer software. Default ""

        Skip options:
            --skipEvi 	            Skip EvidentialGene run in --onlyAsm option. Default "false"
            --skipQC 	            Skip FastQC step. Default "false"
            --skipFilter 	        Skip fastp filtering step. Default "false"
            --skipKegg              Skip kegg analysis. Default "false"
            --skipReport 	        Skip generation of final TransPi report. Default "false"

        Others:
            --minQual 	            Minimum quality score for fastp filtering. Default "25"
            --pipeInstall           PATH to TransPi directory. Default "". If precheck is used this will be added to the nextflow.config automatically.
            --myCondaInstall        PATH to local conda environment of TransPi. Default "". Requires use of --myConda.
                --myConda           Make TransPi use a local conda environemt with all the tools (generated with precheck)

            --envCacheDir           PATH for environment cache directory (either conda or containers). Default "Launch directory of pipeline"
            --getVersions           Get software versions. Default "false"

    """.stripIndent()
}
def fullHelpMessage() {
    log.info """
    ==================================================
      TransPi - Transcriptome Analysis Pipeline v${workflow.manifest.version}
    ==================================================

    Steps:
        1- Run the `precheck_TransPi.sh` to set up the databases and tools (if neccesary) used by TransPi
        2- Run TransPi

        Usage:
            nextflow run TransPi.nf TransPi_analysis_option other_options

            Example usage:
                nextflow run TransPi.nf --all --reads "/PATH/TO/READS/*_R[1,2].fastq.gz" --k 25,41,53 --maxReadLen 75

        Manadatory arguments:
            --all                   Run the entire pipeline (Assemblies, EvidentialGene, Annotation, etc.)
                                    This option also requires arguments --reads, --k, and --maxReadLen
                                    Example:
                                        --reads "/PATH/TO/READS/*_R[1,2].fastq.gz" --k 25,35,55,75,85 --maxReadLen 150
                                        NOTE: Use of quotes is needed for the reads PATH. Kmer list depends on read length.

            --onlyAsm               Run only the Assemblies and EvidentialGene analysis
                                    This option also requires arguments --reads, --k, --maxReadLen

            --onlyEvi               Run only the Evidential Gene analysis
                                    Transcriptome expected to be in a directory called "onlyEvi"

            --onlyAnn               Run only the Annotation analysis (starting from a final assembly)
                                    Transcriptome expected to be in a directory called "onlyAnn"

        Other options:
            -profile                Configuration profile to use. Can use multiple (comma separated)
                    test                Run TransPi with a test dataset
                    conda               Run TransPi with conda.
                    docker              Run TransPi with docker container
                    singularity         Run TransPi with singularity container with all the neccesary tools
                    TransPiContainer    Run TransPi with a single container with all tools

            --help                  Display this message
            --fullHelp              Display this message and examples for running TransPi

        Output options:
            --outdir 	            Name of output directory. Default "results"
            -w, -work 	            Name of working directory. Default "work". Only one dash is needed for -work since it is a nextflow function.
            --tracedir              Name for directory to save pipeline trace files. Default "pipeline_info"

        Additional analyses:
            --rRNAfilter            Remove rRNA from sequences. Requires option --rRNAdb
                --rRNAdb                PATH to database of rRNA sequences to use for filtering of rRNA. Default ""

            --filterSpecies         Perform psytrans filtering of transcriptome. Default "false" Requires options --host and --symbiont
                --host 	                PATH to host (or similar) protein file. Default ""
                --symbiont     	        PATH to symbionts (or similar) protein files. Default ""

            --psyval       	        Psytrans value to train model. Default "160"
            --allBuscos       	    Run BUSCO analysis in all assemblies. Default "false"
            --rescueBusco 	        Generate BUSCO distribution analysis. Default "false"
            --minPerc        	    Mininmum percentage of assemblers require for the BUSCO distribution. Default ".70"
            --shortTransdecoder     Run Transdecoder without the homology searches. Default "false"
            --withSignalP 	        Include SignalP for the annotation. Needs manual installation of CBS-DTU tools. Default "false". Requires --signalp
                --signalp               PATH to SignalP software. Default ""

            --withTMHMM 	        Include TMHMM for the annotation. Needs manual installation of CBS-DTU tools. Default "false". Requires --tmhmm
                --tmhmm                 PATH to TMHMM software. Default ""

            --withRnammer 	        Include Rnammer for the annotation. Needs manual installation of CBS-DTU tools. Default "false". Requires --rnam
                --rnam                  PATH to Rnammer software. Default ""

        Skip options:
            --skipEvi 	            Skip EvidentialGene run in --onlyAsm option. Default "false"
            --skipQC 	            Skip FastQC step. Default "false"
            --skipFilter 	        Skip fastp filtering step. Default "false"
            --skipKegg              Skip kegg analysis. Default "false"
            --skipReport 	        Skip generation of final TransPi report. Default "false"

        Others:
            --pipeInstall           PATH to TransPi directory. Default "". If precheck is used this will be added to the nextflow.config automatically.
            --minQual 	            Minimum quality score for fastp filtering. Default "25"
            --myCondaInstall        PATH to local conda environment of TransPi. Default "". Requires use of --myConda.
                --myConda           Make TransPi use a local conda environemt with all the tools (generated with precheck)

            --envCacheDir           PATH for environment cache directory (either conda or containers). Default "Launch directory of pipeline"
            --getVersions           Get software versions. Default "false"

        #################################################################################################

                                    Various examples on how to deploy TransPi

        #################################################################################################

        I. Steps for running on a local cluster and local conda installation of TransPi

            1- Run the `precheck_TransPi.sh` to install tools with conda, set up the databases and directories for TransPi
            2- Run TransPi

            Usage:

                nextflow run TransPi.nf --all -profile conda --myConda OTHER_PARAMETERS_HERE

        #################################################################################################

        II. Steps for running on a local cluster and conda installation by Nextflow

            1- Run the `precheck_TransPi.sh` to set up the databases for TransPi
            2- Run TransPi

            Usage:

                nextflow run TransPi.nf --all -profile conda OTHER_PARAMETERS_HERE

            NOTE:
                A conda environment will be created for each process.

        #################################################################################################

        III. Steps for running with docker

            1- Run the `container_precheck_TransPi.sh` to set up the databases for TransPi
            2- Run TransPi

            Usage:

                nextflow run TransPi.nf --all -profile docker OTHER_PARAMETERS_HERE

            NOTE:
                All necesary tools for running TransPi are pre-installed in the container

        #################################################################################################

        IV. Steps for running with singualarity

            1- Run the `container_precheck_TransPi.sh` to set up the databases for TransPi
            2- Run TransPi

            Usage:

                nextflow run TransPi.nf --all -profile singularity OTHER_PARAMETERS_HERE

            NOTE:
                All necesary tools for running TransPi are pre-installed in the container

        #################################################################################################

        V. Steps for running with a container engine and TransPi container.

            1- Run the `precheck_TransPi.sh` to install tools, set up the databases and directories for TransPi
            2- Run TransPi

            Usage:

                nextflow run TransPi.nf --all -profile singularity,TransPiContainer

            NOTE:
                This will run TransPi using a single container only.

        #################################################################################################

        VI. Steps for running with multiple profiles.

            1- Run the `precheck_TransPi.sh` to set up the databases for TransPi
            2- Run TransPi

            Usage:

                nextflow run TransPi.nf --all -profile singularity,TransPiContainer,test

            NOTE:
                This will run TransPi using a test dataset, with singularity and the TransPi container.

        #################################################################################################
    """.stripIndent()
}

def workdir = System.getProperty("user.dir");

// Show help message
if (params.help) {
    helpMessage()
    exit 0
} else if (params.fullHelp) {
    fullHelpMessage()
    exit 0
}

def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

def checkArgs() {
    if (!params.k) {
        println("\n\t\033[0;31mKmer list not specified.\n\tFor more info use `nextflow run TransPi.nf --help`\n\033[0m")
        exit 0
    }
    if (!params.reads && !params.readsTest) {
        println("\n\t\033[0;31mReads mandatory argument not specified.\n\tFor more info use `nextflow run TransPi.nf --help`\n\033[0m")
        exit 0
    }
    if (!params.maxReadLen) {
        println("\n\t\033[0;31mMax read length argument not specified.\n\tFor more info use `nextflow run TransPi.nf --help`\n\033[0m")
        exit 0
    }
}

if (params.condaActivate && params.myConda && params.myCondaInstall == "") {
    println("\n\t\033[0;31mNeed to specify the local conda installation in parameter \"myCondaInstall\" in the config file.\n\tFor more info use `nextflow run TransPi.nf --help`\n\033[0m")
    exit 0
}

if (params.rRNAfilter) {
    if (params.rRNAdb == "" || params.rRNAdb == true) {
        println("\n\t\033[0;31mNeed to provide the PATH to the rRNA database.\n\tFor more info use `nextflow run TransPi.nf --help`\n\033[0m")
        exit 0
    }
}

if (params.all) {
    log.info """\
            ==================================================
              TransPi - Transcriptome Analysis Pipeline v${workflow.manifest.version}
            ==================================================
            TransPi.nf Directory:   ${projectDir}
            Launch Directory:       ${launchDir}
            Results Directory:      ${launchDir}/${params.outdir}
            Work Directory:         ${workDir}
            TransPi DBs:            ${params.pipeInstall}
            Uniprot DB:             ${params.uniprot}
            Busco DB:               ${params.busco4db}
            Reads Directory:        ${params.reads}
            Read Length:            ${params.maxReadLen}
            Kmers:                  ${params.k}
            """.stripIndent()
    checkArgs()
} else if (params.onlyAnn) {
    log.info """\
            ==================================================
              TransPi - Transcriptome Analysis Pipeline v${workflow.manifest.version}
            ==================================================
            TransPi.nf Directory:   ${projectDir}
            Launch Directory:       ${launchDir}
            Results Directory:      ${launchDir}/${params.outdir}
            Work Directory:         ${workDir}
            TransPi DBs:            ${params.pipeInstall}
            Uniprot DB:             ${params.uniprot}
            """.stripIndent()
} else if (params.onlyAsm) {
    log.info """\
            ==================================================
              TransPi - Transcriptome Analysis Pipeline v${workflow.manifest.version}
            ==================================================
            TransPi.nf Directory:   ${projectDir}
            Launch Directory:       ${launchDir}
            Results Directory:      ${launchDir}/${params.outdir}
            Work Directory:         ${workDir}
            TransPi DBs:            ${params.pipeInstall}
            Busco DB:               ${params.busco4db}
            Reads Directory:        ${params.reads}
            Read Length:            ${params.maxReadLen}
            Kmers:                  ${params.k}
            """.stripIndent()
    checkArgs()
} else if (params.onlyEvi){
    log.info """\
            ==================================================
              TransPi - Transcriptome Analysis Pipeline v${workflow.manifest.version}
            ==================================================
            TransPi.nf Directory:   ${projectDir}
            Launch Directory:       ${launchDir}
            Results Directory:      ${launchDir}/${params.outdir}
            Work Directory:         ${workDir}
            TransPi DBs:            ${params.pipeInstall}
            Busco DB:               ${params.busco4db}
            """.stripIndent()
}

if (params.readsTest) {
    println("\n\tRunning TransPi with TEST dataset\n")
    Channel
        .from(params.readsTest)
        .map{ row -> [ row[0], [ file(row[1][0],checkIfExists: true),file(row[2][0],checkIfExists: true) ] ] }
        .ifEmpty{ exit 1, "params.readsTest was empty - no input files supplied" }
        .into{ reads_ch; reads_qc_ch; report_reads }
} else {
    println("\n\tRunning TransPi with your dataset\n")
    if ( ( params.all || params.onlyAsm ) && !params.readsTest) {
        Channel
            .fromFilePairs("${params.reads}", checkIfExists: true)
            .into{ reads_ch; reads_qc_ch; report_reads }
    }
}

if (params.onlyAsm) {

    println("\n\tRunning assemblies and Evidential Gene analysis only \n")

    if (!params.skipQC) {

        process fasqc_OAS {

            label 'low_cpus'

            tag "${sample_id}"

            publishDir "${launchDir}/${params.outdir}/fastqc", mode: "copy", overwrite: true

            conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::fastqc=0.11.9=0" : null)
            if (params.oneContainer){ container "${params.TPcontainer}" } else {
            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0" : "quay.io/biocontainers/fastqc:0.11.9--0")
            }

            input:
                tuple sample_id, file(reads) from reads_qc_ch

            output:
                tuple sample_id, file("*_fastqc.{zip,html}") into fastqc_results_OAS

            script:
                """
                fastqc --quiet --threads $task.cpus $reads
                """
        }
    }

    if (!params.skipFilter) {

        process fastp_OAS {

            label 'med_cpus'

            tag "${sample_id}"

            publishDir "${launchDir}/${params.outdir}/filter", mode: "copy", overwrite: true, pattern: "*.fastp.{json,html}"

            conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge conda-forge::pigz=2.3.4=hed695b0_1 conda-forge::jq=1.6=h14c3975_1000 bioconda::fastp=0.20.1=h8b12597_0" : null)
            if (params.oneContainer){ container "${params.TPcontainer}" } else {
            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/fastp:0.20.1--h8b12597_0" : "quay.io/biocontainers/fastp:0.20.1--h8b12597_0")
            }

            input:
                tuple sample_id, file(reads) from reads_ch

            output:
                tuple sample_id, file("*.fastp.{json,html}") into fastp_results_OAS
                tuple sample_id, file("${sample_id}.fastp.json") into fastp_stats_OAS_ch
                tuple sample_id, file("*${sample_id}.filter.fq") into reads_ass_ch_OAS
                tuple sample_id, file("${sample_id}_filter.R1.fq.gz"), file("${sample_id}_filter.R2.fq.gz") into save_filter_reads

            script:
                """
                echo ${sample_id}

                fastp -i ${reads[0]} -I ${reads[1]} -o left-${sample_id}.filter.fq -O right-${sample_id}.filter.fq --detect_adapter_for_pe \
                --average_qual ${params.minQual} --overrepresentation_analysis --html ${sample_id}.fastp.html --json ${sample_id}.fastp.json \
                --thread ${task.cpus} --report_title ${sample_id}

                cp left-${sample_id}.filter.fq ${sample_id}_filter.R1.fq
                cp right-${sample_id}.filter.fq ${sample_id}_filter.R2.fq

                pigz --best --force -p ${task.cpus} -r ${sample_id}_filter.R1.fq
                pigz --best --force -p ${task.cpus} -r ${sample_id}_filter.R2.fq
                """
        }

        process fastp_stats_OAS {

            label 'low_cpus'

            tag "${sample_id}"

            publishDir "${launchDir}/${params.outdir}/filter", mode: "copy", overwrite: true, pattern: "*.fastp.{json,html}"

            conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge conda-forge::pigz=2.3.4=hed695b0_1 conda-forge::jq=1.6=h14c3975_1000 bioconda::fastp=0.20.1=h8b12597_0" : null)
            if (params.oneContainer){ container "${params.TPcontainer}" } else {
            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/fastp:0.20.1--h8b12597_0" : "quay.io/biocontainers/fastp:0.20.1--h8b12597_0")
            }

            input:
                tuple sample_id, file(json) from fastp_stats_OAS_ch

            output:
                tuple sample_id, file("*.csv") into fastp_csv_OAS

            script:
                """
                echo ${sample_id}
                bash get_readstats.sh ${json}
                bash get_readqual.sh ${json}
                """

        }

        if (params.saveReads) {

            process save_filter_reads_OAS {

                tag "${sample_id}"

                publishDir "${launchDir}/${params.outdir}/saveReads/filtering", mode: "copy", overwrite: true, pattern: "*_R{1,2}.filter.fq.gz"

                input:
                    tuple sample_id, file(r1), file(r2) from save_filter_reads

                output:
                    tuple sample_id, file("*.filter.fq.gz") into save_filter_reads_out

                script:
                    """
                    cat $r1 >${sample_id}_filter.R1.fq.gz
                    cat $r2 >${sample_id}_filter.R2.fq.gz
                    """
            }
        }

    } else {
        reads_ch
            .set{ reads_ass_ch_OAS }
        fastp_results_OAS = Channel.empty()
        reads_ass_ch_OAS = Channel.empty()
        fastp_csv_OAS = Channel.empty()
    }

    if (!params.skipNormalization) {

        process normalize_reads_OAS {

            label 'med_mem'

            tag "${sample_id}"

            conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge -c bioconda trinity=2.9.1 pigz=2.3.4" : null)
            if (params.oneContainer){ container "${params.TPcontainer}" } else {
            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/mulled-v2-13733e73a40b40f427f5bc7edcfbd4f0dbd27ae0:fcba8b2e03b4c752f3076b94d7c315f77345e143-0" : "quay.io/biocontainers/mulled-v2-13733e73a40b40f427f5bc7edcfbd4f0dbd27ae0:fcba8b2e03b4c752f3076b94d7c315f77345e143-0")
            }

            input:
                tuple sample_id, file(reads) from reads_ass_ch_OAS

            output:
                tuple sample_id, file("left-${sample_id}.norm.fq"), file("right-${sample_id}.norm.fq") into ( norm_reads_soap_OAS, norm_reads_velvet_OAS, norm_reads_trinity_OAS, norm_reads_spades_OAS, norm_reads_transabyss_OAS, reads_rna_quast_OAS )
                tuple sample_id, file("left-${sample_id}.norm.fq"), file("right-${sample_id}.norm.fq") into ( mapping_reads_trinity_OAS, mapping_reads_evi_OAS )
                tuple sample_id, file("${sample_id}_norm.R1.fq.gz"), file("${sample_id}_norm.R2.fq.gz") into ( save_norm_reads )
                tuple sample_id, file("${sample_id}_normStats.txt") into norm_report_OAS

            script:
                //def mem=(task.memory)
                //def mem_MB=(task.memory.toMega())
            if (!params.skipFilter) {
                """
                echo ${sample_id}

                echo -e "\\n-- Starting Normalization --\\n"

                mem=\$( echo ${task.memory} | cut -f 1 -d " " )

                insilico_read_normalization.pl --seqType fq -JM \${mem}G --max_cov 100 --min_cov 1 --left ${reads[0]} --right ${reads[1]} --pairs_together --PARALLEL_STATS --CPU ${task.cpus}

                echo -e "\\n-- DONE with Normalization --\\n"

                cat .command.out | grep "stats_file" -A 3 | tail -n 3 >${sample_id}_normStats.txt

                cp left.norm.fq left-"${sample_id}".norm.fq
                cp right.norm.fq right-"${sample_id}".norm.fq

                mv left.norm.fq ${sample_id}_norm.R1.fq
                mv right.norm.fq ${sample_id}_norm.R2.fq

                pigz --best --force -p ${task.cpus} -r ${sample_id}_norm.R1.fq
                pigz --best --force -p ${task.cpus} -r ${sample_id}_norm.R2.fq
                """
            } else {
                """
                echo ${sample_id}
                zcat ${reads[0]} >left-${sample_id}.fq &
                zcat ${reads[1]} >right-${sample_id}.fq

                echo -e "\\n-- Starting Normalization --\\n"

                mem=\$( echo ${task.memory} | cut -f 1 -d " " )

                insilico_read_normalization.pl --seqType fq -JM \${mem}G --max_cov 100 --min_cov 1 --left left-${sample_id}.fq --right right-${sample_id}.fq --pairs_together --PARALLEL_STATS --CPU ${task.cpus}

                echo -e "\\n-- DONE with Normalization --\\n"

                cat .command.out | grep "stats_file" -A 3 | tail -n 3 >${sample_id}_normStats.txt

                cp left.norm.fq left-"${sample_id}".norm.fq
                cp right.norm.fq right-"${sample_id}".norm.fq

                mv left.norm.fq ${sample_id}_norm.R1.fq
                mv right.norm.fq ${sample_id}_norm.R2.fq

                pigz --best --force -p ${task.cpus} -r ${sample_id}_norm.R1.fq
                pigz --best --force -p ${task.cpus} -r ${sample_id}_norm.R2.fq
                """
            }
        }

        if (params.saveReads) {

            process save_norm_reads_OAS {

                tag "${sample_id}"

                publishDir "${launchDir}/${params.outdir}/saveReads/normalization", mode: "copy", overwrite: true, pattern: "*_R{1,2}.norm.fq.gz"

                input:
                    tuple sample_id, file(r1), file(r2) from save_norm_reads

                output:
                    tuple sample_id, file("*.norm.fq.gz") into save_norm_reads_out

                script:
                    """
                    cat $r1 >${sample_id}_norm.R1.fq.gz
                    cat $r2 >${sample_id}_norm.R2.fq.gz
                    """
            }
        }
    }

    if (params.skipFilter && params.skipNormalization) {

        process prepare_reads_OAS {

            label 'low_cpus'

            tag "${sample_id}"

            input:
                tuple sample_id, file(reads) from reads_ass_ch_OAS

            output:
                tuple sample_id, file("left-${sample_id}.fq"), file("right-${sample_id}.fq") into ( norm_reads_soap_OAS, norm_reads_velvet_OAS, norm_reads_trinity_OAS, norm_reads_spades_OAS, norm_reads_transabyss_OAS, reads_rna_quast_OAS )
                tuple sample_id, file("left-${sample_id}.fq"), file("right-${sample_id}.fq") into ( mapping_reads_trinity_OAS, mapping_reads_evi_OAS )

            script:
            if (hasExtension(params.reads, 'gz')) {
                """
                echo ${sample_id}
                zcat ${reads[0]} >left-${sample_id}.fq &
                zcat ${reads[1]} >right-${sample_id}.fq
                """
            } else {
                """
                echo ${sample_id}
                cat ${reads[0]} >left-${sample_id}.fq &
                cat ${reads[1]} >right-${sample_id}.fq
                """
            }
        }
    }

    process trinity_assembly_OAS {

        label 'med_mem'

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/assemblies", mode: "copy", overwrite: true

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::trinity=2.9.1=h8b12597_1" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/trinity:2.9.1--h8b12597_1" : "quay.io/biocontainers/trinity:2.9.1--h8b12597_1")
        }

        input:
            tuple sample_id, file(left), file(right) from norm_reads_trinity_OAS

        output:
            tuple sample_id, file("${sample_id}.Trinity.fa") into ( assemblies_ch_trinity_OAS, busco3_ch_trinity_OAS, busco4_ch_trinity_OAS, mapping_trinity_OAS )

        script:
            """
            mem=\$( echo ${task.memory} | cut -f 1 -d " " )

            Trinity --max_memory \${mem}G --seqType fq --left ${left} --right ${right} --CPU ${task.cpus} --no_normalize_reads --full_cleanup --output trinity_out_dir

            mv trinity_out_dir.Trinity.fasta ${sample_id}.Trinity.fa
            """
    }

    process soap_assembly_OAS {

        label 'med_mem'

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/assemblies", mode: "copy", overwrite: true

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::soapdenovo-trans=1.04=ha92aebf_2" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/soapdenovo-trans:1.04--ha92aebf_2" : "quay.io/biocontainers/soapdenovo-trans:1.04--ha92aebf_2")
        }

        input:
            val k from "${params.k}"
            tuple sample_id, file(left), file(right) from norm_reads_soap_OAS

        output:
            tuple sample_id, file("${sample_id}.SOAP.fa") into assemblies_ch_soap_OAS

        script:
            """
            echo -e "\\n-- Generating SOAP config file --\\n"
            echo "maxReadLen="${params.maxReadLen} >>config.txt
            echo "[LIB]" >>config.txt
            echo "rd_len_cutof="${params.rd_len_cutof} >>config.txt
            #echo "avg_ins="${params.avg_ins} >>config.txt
            echo "reverse_seq="${params.reverse_seq} >>config.txt
            echo "asm_flags="${params.asm_flags} >>config.txt
            echo "map_len="${params.map_len} >>config.txt
            echo "q1="${left} >>config.txt
            echo "q2="${right} >>config.txt

            echo -e "\\n-- Starting SOAP assemblies --\\n"

            for x in `echo $k | tr "," " "`;do
                echo -e "\\n-- SOAP k\${x} --\\n"
                SOAPdenovo-Trans-127mer all -s config.txt -K \${x} -o output\${x} -p ${task.cpus}
                sed -i "s/>/>SOAP.k\${x}./g" output\${x}.scafSeq
            done

            echo -e "\\n-- Finished with the assemblies --\\n"

            cat output*.scafSeq >${sample_id}.SOAP.fa

            rm -rf output*
            """
    }

    process velvet_oases_assembly_OAS {

        label 'med_mem'

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/assemblies", mode: "copy", overwrite: true

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge -c bioconda velvet=1.2.10 oases=0.2.09" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/mulled-v2-8ce10492777ba3fb1db6e6e13fa9b78ac116db2f:f54a9246f1216443f2e0f6de9ec5908ca882f710-0" : "quay.io/biocontainers/mulled-v2-8ce10492777ba3fb1db6e6e13fa9b78ac116db2f:f54a9246f1216443f2e0f6de9ec5908ca882f710-0")
        }

        input:
            val k from "${params.k}"
            tuple sample_id, file(left), file(right) from norm_reads_velvet_OAS

        output:
            tuple sample_id, file("${sample_id}.Velvet.fa") into assemblies_ch_velvet_OAS

        script:
            """
            echo -e "\\n-- Starting with Velveth --\\n"
            for x in `echo $k | tr "," " "`;do
                echo -e "\\n-- k\${x} --\\n"
                velveth oases.\${x} \${x} -shortPaired -fastq -separate ${left} ${right}
            done

            echo -e "\\n-- Starting with Velvetg --\\n"
            for x in `echo $k | tr "," " "`;do
                echo -e "\\n-- vg \${x} --\\n"
                velvetg oases.\${x} -read_trkg yes
            done

            echo -e "\\n-- Starting with Oases --\\n"
            for x in `echo $k | tr "," " "`;do
                echo -e "\\n-- oases \${x} --\\n"
                oases oases.\${x}
            done

            echo -e "\\n-- Finished with Velvet/Oases assemblies --\\n"

            for x in `echo $k | tr "," " "`;do
                sed -i "s/>/>Velvet.k\${x}./g" oases.\${x}/contigs.fa
            done

            cat oases.*/contigs.fa >${sample_id}.Velvet.fa

            rm -rf oases.*
            """
    }

    process rna_spades_assembly_OAS {

        label 'med_mem'

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/assemblies", mode: "copy", overwrite: true

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::spades=3.14.0=h2d02072_0" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/spades:3.14.0--h2d02072_0" : "quay.io/biocontainers/spades:3.14.0--h2d02072_0")
        }

        input:
            val k from "${params.k}"
            tuple sample_id, file(left), file(right) from norm_reads_spades_OAS

        output:
            tuple sample_id, file("${sample_id}.SPADES.fa") into assemblies_ch_spades_OAS

        script:
            """
            echo -e "\\n-- Starting rnaSPADES assemblies --\\n"

            mem=\$( echo ${task.memory} | cut -f 1 -d " " )

            for x in `echo $k | tr "," " "`;do
                echo -e "\\n-- rnaSPADES k\${x} --\\n"
                rnaspades.py -1 ${left} -2 ${right} -o ${sample_id}_spades_\${x} -t ${task.cpus} -k \${x} -m \${mem}
            done

            echo -e "\\n-- Finished with the assemblies --\\n"

            for x in `echo $k | tr "," " "`;do
                sed -i "s/>/>SPADES.k\${x}./g" ${sample_id}_spades_\${x}/transcripts.fasta
            done

            cat ${sample_id}_spades_*/transcripts.fasta >${sample_id}.SPADES.fa

            rm -rf ${sample_id}_spades_*
            """
    }

    process transabyss_assembly_OAS {

        label 'med_mem'

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/assemblies", mode: "copy", overwrite: true

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::transabyss=2.0.1=py_6" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/transabyss:2.0.1--py_6" : "quay.io/biocontainers/transabyss:2.0.1--py_6")
        }

        input:
            val k from "${params.k}"
            tuple sample_id, file(left), file(right) from norm_reads_transabyss_OAS

        output:
            tuple sample_id, file("${sample_id}.TransABySS.fa") into assemblies_ch_transabyss_OAS

        script:
            """
            echo -e "\\n-- Starting Trans-ABySS assemblies --\\n"

            for x in `echo $k | tr "," " "`;do
                echo -e "\\n-- Trans-ABySS k\${x} --\\n"
                transabyss -k \${x} --pe ${left} ${right} --outdir ${sample_id}_transabyss_\${x} --name k\${x}.transabyss.fa --threads ${task.cpus} -c 12 --length 200
            done

            echo -e "\\n-- Finished with the assemblies --\\n"

            for x in `echo $k | tr "," " "`;do
                sed -i "s/>/>TransABySS.k\${x}./g" ${sample_id}_transabyss_\${x}/k\${x}.transabyss.fa-final.fa
            done

            cat ${sample_id}_transabyss_*/k*.transabyss.fa-final.fa >${sample_id}.TransABySS.fa

            rm -rf ${sample_id}_transabyss_*
            """
    }

    all_assemblies = Channel.create()
    assemblies_ch_trinity_OAS.mix( assemblies_ch_transabyss_OAS, assemblies_ch_spades_OAS, assemblies_ch_velvet_OAS, assemblies_ch_soap_OAS ).groupTuple(by:0,size:5).into(all_assemblies)

    process evigene_OAS {

        label 'med_mem'

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/evigene", mode: "copy", overwrite: true

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::cd-hit=4.8.1 bioconda::exonerate=2.4 bioconda::blast=2.2.31" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/mulled-v2-962eae98c9ff8d5b31e1df7e41a355a99e1152c4:5aac6d1d2253d47aee81f01cc070a17664c86f07-0" : "quay.io/biocontainers/mulled-v2-962eae98c9ff8d5b31e1df7e41a355a99e1152c4:5aac6d1d2253d47aee81f01cc070a17664c86f07-0")
        }

        input:
            tuple sample_id, file(assemblies) from all_assemblies

        output:
            tuple sample_id, file("*.combined.okay.fa") into ( evigene_ch_busco3_OAS, evigene_ch_busco4_OAS, evigene_ch_rna_quast_OAS, mapping_evigene_OAS )
            tuple sample_id, file("${sample_id}.combined.fa"), file("${sample_id}.combined.okay.fa") into evigene_summary_OAS

        script:
            def mem_MB=(task.memory.toMega())

            """
            echo -e "\\n-- Starting EviGene --\\n"

            cat ${assemblies} >${sample_id}.combined.fa

            $evi/scripts/prot/tr2aacds.pl -tidy -NCPU ${task.cpus} -MAXMEM ${mem_MB} -log -cdna ${sample_id}.combined.fa

            echo -e "\\n-- DONE with EviGene --\\n"

            cp okayset/*combined.okay*.fa ${sample_id}.combined.okay.fa

            if [ -d tmpfiles/ ];then
                rm -rf tmpfiles/
            fi
            """
    }

    // check groupTuple
    rna_quast_OAS = Channel.create()
    reads_rna_quast_OAS.join( evigene_ch_rna_quast_OAS ).into( rna_quast_OAS )

    process rna_quast_OAS {

        label 'med_mem'

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/rnaQuast", mode: "copy", overwrite: true

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::rnaquast=2.0.1=0" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/rnaquast:2.0.1--0" : "quay.io/biocontainers/rnaquast:2.0.1--0")
        }

        input:
            tuple sample_id, file(r1), file(r2), file(assembly) from rna_quast_OAS

        output:
            tuple sample_id, file("${sample_id}.rna_quast") into rna_quast_sum_OAS

        script:
            """
            rnaQUAST.py --transcripts ${assembly} -1 ${r1} -2 ${r2} -o ${sample_id}.rna_quast -t ${task.cpus} --blat
            """
    }

    mapping_evigene_in_OAS=Channel.create()
    mapping_evigene_OAS.mix( mapping_reads_evi_OAS ).groupTuple(by:0,size:2).into(mapping_evigene_in_OAS)

    process mapping_evigene_OAS {

        label 'big_cpus'

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/mapping", mode: "copy", overwrite: true

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::bowtie2=2.3.5.1=py36he513fc3_0" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/bowtie2:2.3.5.1--py36he513fc3_0" : "quay.io/biocontainers/bowtie2:2.3.5.1--py36he513fc3_0")
        }

        input:
            tuple sample_id, file(files), file(files2) from mapping_evigene_in_OAS

        output:
            tuple sample_id, file("log*") into mapping_evi_results

        script:
            """
            a=\$( echo $files $files2 )
            ass=\$( echo \$a | tr " " "\\n" | grep ".combined.okay.fa" )
            r1=\$( echo \$a | tr " " "\\n" | grep "left-${sample_id}.norm.fq" )
            r2=\$( echo \$a | tr " " "\\n" | grep "right-${sample_id}.norm.fq" )
            bowtie2-build \${ass} \${ass} --threads ${task.cpus}
            bowtie2 -x \${ass} -1 \${r1} -2 \${r2} -p ${task.cpus} -S \${ass}.sam 2>&1 | tee -a log_\${ass}.txt
            rm *.sam
            """

    }

    mapping_trinity_in_OAS=Channel.create()
    mapping_trinity_OAS.mix( mapping_reads_trinity_OAS ).groupTuple(by:0,size:2).into(mapping_trinity_in_OAS)

    process mapping_trinity_OAS {

        label 'big_cpus'

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/mapping", mode: "copy", overwrite: true

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::bowtie2=2.3.5.1=py36he513fc3_0" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/bowtie2:2.3.5.1--py36he513fc3_0" : "quay.io/biocontainers/bowtie2:2.3.5.1--py36he513fc3_0")
        }

        input:
            tuple sample_id, file(files), file(files2) from mapping_trinity_in_OAS

        output:
            tuple sample_id, file("log*") into mapping_trinity_results

        script:
            """
            a=\$( echo $files $files2 )
            ass=\$( echo \$a | tr " " "\\n" | grep ".Trinity.fa" )
            r1=\$( echo \$a | tr " " "\\n" | grep "left-${sample_id}.norm.fq" )
            r2=\$( echo \$a | tr " " "\\n" | grep "right-${sample_id}.norm.fq" )
            bowtie2-build \${ass} \${ass} --threads ${task.cpus}
            bowtie2 -x \${ass} -1 \${r1} -2 \${r2} -p ${task.cpus} -S \${ass}.sam 2>&1 | tee -a log_\${ass}.txt
            rm *.sam
            """

    }

    process summary_evigene_individual_OAS {

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/stats", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.combined.fa"), file("${sample_id}.combined.okay.fa") from evigene_summary_OAS

        output:
            tuple sample_id, file("${sample_id}_sum_preEG.txt"), file("${sample_id}_sum_EG.txt") into final_sum_1_OAS
            tuple sample_id, file("*.csv") into summary_evi_csv_OAS

        script:
            """
            #Summary of total number of transcripts
            echo -e "- Number of transcripts before Evidential Genes\\n" >>${sample_id}_sum_preEG.txt
            echo -e "- Individual ${sample_id} \\n" >>${sample_id}_sum_preEG.txt
            echo -e "\\t Total transcripts:" >>${sample_id}_sum_preEG.txt
            num=\$( cat ${sample_id}.combined.fa | grep -c ">" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}_sum_preEG.txt
            echo -e "\\t Trinity" >>${sample_id}_sum_preEG.txt
            num=\$( cat ${sample_id}.combined.fa | grep -c ">TRINITY" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}_sum_preEG.txt
            echo -e "\\t SOAP" >>${sample_id}_sum_preEG.txt
            num=\$( cat ${sample_id}.combined.fa | grep -c ">SOAP" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}_sum_preEG.txt
            echo -e "\\t Velvet/Oases" >>${sample_id}_sum_preEG.txt
            num=\$( cat ${sample_id}.combined.fa | grep -c ">Velvet" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}_sum_preEG.txt
            echo -e "\\t rna-SPADES" >>${sample_id}_sum_preEG.txt
            num=\$( cat ${sample_id}.combined.fa | grep -c ">SPADES" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}_sum_preEG.txt
            echo -e "\\t Trans-ABySS" >>${sample_id}_sum_preEG.txt
            num=\$( cat ${sample_id}.combined.fa | grep -c ">TransABySS" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}_sum_preEG.txt

            # csv report
            echo "Total,Trinity,SOAP,Velvet,SPADES,TransBySS" >${sample_id}_sum_preEG.csv
            total=\$( cat ${sample_id}.combined.fa | grep -c ">" )
            trinity=\$( cat ${sample_id}.combined.fa | grep -c ">TRINITY" )
            soap=\$( cat ${sample_id}.combined.fa | grep -c ">SOAP" )
            velvet=\$( cat ${sample_id}.combined.fa | grep -c ">Velvet" )
            spades=\$( cat ${sample_id}.combined.fa | grep -c ">SPADES" )
            transabyss=\$( cat ${sample_id}.combined.fa | grep -c ">TransABySS" )
            echo "\${total},\${trinity},\${soap},\${velvet},\${spades},\${transabyss}" >>${sample_id}_sum_preEG.csv

            #Summary of transcripts after EvidentialGenes
            echo -e "- Number of transcripts by individual after EvidentialGenes\\n" >>${sample_id}_sum_EG.txt
            echo -e "- Individual ${sample_id} \\n" >>${sample_id}_sum_EG.txt
            echo -e "\\t Total transcripts:" >>${sample_id}_sum_EG.txt
            num=\$( cat ${sample_id}.combined.okay.fa | grep -c ">" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}_sum_EG.txt
            echo -e "\\t Trinity" >>${sample_id}_sum_EG.txt
            num=\$( cat ${sample_id}.combined.okay.fa | grep -c ">TRINITY" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}_sum_EG.txt
            echo -e "\\t SOAP" >>${sample_id}_sum_EG.txt
            num=\$( cat ${sample_id}.combined.okay.fa | grep -c ">SOAP" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}_sum_EG.txt
            echo -e "\\t Velvet/Oases" >>${sample_id}_sum_EG.txt
            num=\$( cat ${sample_id}.combined.okay.fa | grep -c ">Velvet" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}_sum_EG.txt
            echo -e "\\t rna-SPADES" >>${sample_id}_sum_EG.txt
            num=\$( cat ${sample_id}.combined.okay.fa | grep -c ">SPADES" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}_sum_EG.txt
            echo -e "\\t Trans-ABySS" >>${sample_id}_sum_EG.txt
            num=\$( cat ${sample_id}.combined.fa | grep -c ">TransABySS" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}_sum_EG.txt

            # csv report after evigene
            echo "Total,Trinity,SOAP,Velvet,SPADES,TransBySS" >${sample_id}_sum_EG.csv
            total=\$( cat ${sample_id}.combined.okay.fa | grep -c ">" )
            trinity=\$( cat ${sample_id}.combined.okay.fa | grep -c ">TRINITY" )
            soap=\$( cat ${sample_id}.combined.okay.fa | grep -c ">SOAP" )
            velvet=\$( cat ${sample_id}.combined.okay.fa | grep -c ">Velvet" )
            spades=\$( cat ${sample_id}.combined.okay.fa | grep -c ">SPADES" )
            transabyss=\$( cat ${sample_id}.combined.okay.fa | grep -c ">TransABySS" )
            echo "\${total},\${trinity},\${soap},\${velvet},\${spades},\${transabyss}" >>${sample_id}_sum_EG.csv
            """

    }

    process busco3_OAS {

        label 'med_cpus'

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/busco3", mode: "copy", overwrite: true

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::busco=3.0.2=py_13" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/busco:3.0.2--py_13" : "quay.io/biocontainers/busco:3.0.2--py_13")
        }

        input:
            tuple sample_id, file("${sample_id}.combined.okay.fa") from evigene_ch_busco3_OAS

        output:
            tuple sample_id, file("run_${sample_id}.TransPi.bus3") into busco3_ch
            tuple sample_id, file("*${sample_id}.TransPi.bus3.txt") into ( busco3_summary_OAS, busco3_comp_1_OAS )

        script:
            """
            echo -e "\\n-- Starting BUSCO --\\n"

            run_BUSCO.py -i ${sample_id}.combined.okay.fa -o ${sample_id}.TransPi.bus3 -l ${params.busco3db} -m tran -c ${task.cpus}

            echo -e "\\n-- DONE with BUSCO --\\n"

            cp run_${sample_id}.TransPi.bus3/short_summary_${sample_id}.TransPi.bus3.txt .
            """
    }

    process busco3_tri_OAS {

        label 'med_cpus'

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/busco3", mode: "copy", overwrite: true

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::busco=3.0.2=py_13" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/busco:3.0.2--py_13" : "quay.io/biocontainers/busco:3.0.2--py_13")
        }

        input:
            tuple sample_id, file("${sample_id}.Trinity.fa") from busco3_ch_trinity_OAS

        output:
            tuple sample_id, file("*${sample_id}.Trinity.bus3.txt") into ( busco3_ch_trinity_sum_OAS, busco3_comp_2_OAS )

        script:
            """
            echo -e "\\n-- Starting BUSCO --\\n"

            run_BUSCO.py -i ${sample_id}.Trinity.fa -o ${sample_id}.Trinity.bus3 -l ${params.busco3db} -m tran -c ${task.cpus}

            echo -e "\\n-- DONE with BUSCO --\\n"

            cp run_${sample_id}.Trinity.bus3/short_summary_${sample_id}.Trinity.bus3.txt .
            """
    }

    process busco4_OAS {

        label 'med_cpus'

        tag "${sample_id}"

        // change container in oneContainer option
        conda (params.condaActivate && params.myConda ? params.cenv : params.condaActivate ? "-c conda-forge bioconda::busco=4.1.4=py_0" : null)
        if (params.oneContainer){ container "${params.v4container}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/busco:4.0.5--pyr36_0" : "ezlabgva/busco:v4.0.5_cv1")
        }

        publishDir "${launchDir}/${params.outdir}/busco4", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.combined.okay.fa") from evigene_ch_busco4_OAS

        output:
            tuple sample_id, file("${sample_id}.TransPi.bus4") into busco4_ch
            tuple sample_id, file("*${sample_id}.TransPi.bus4.txt") into ( busco4_summary_OAS, busco4_comp_1_OAS )

        script:
            """
            echo -e "\\n-- Starting BUSCO --\\n"

            busco -i ${sample_id}.combined.okay.fa -o ${sample_id}.TransPi.bus4 -l ${params.busco4db} -m tran -c ${task.cpus} --offline

            echo -e "\\n-- DONE with BUSCO --\\n"

            cp ${sample_id}.TransPi.bus4/short_summary.*.${sample_id}.TransPi.bus4.txt .
            """
    }

    process busco4_tri_OAS {

        label 'med_cpus'

        tag "${sample_id}"

        // change container in oneContainer option
        conda (params.condaActivate && params.myConda ? params.cenv : params.condaActivate ? "-c conda-forge bioconda::busco=4.1.4=py_0" : null)
        if (params.oneContainer){ container "${params.v4container}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/busco:4.0.5--pyr36_0" : "ezlabgva/busco:v4.0.5_cv1")
        }

        publishDir "${launchDir}/${params.outdir}/busco4", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.Trinity.fa") from busco4_ch_trinity_OAS

        output:
            tuple sample_id, file("*${sample_id}.Trinity.bus4.txt") into ( busco4_ch_trinity_sum_OAS, busco4_comp_2_OAS )

        script:
            """
            echo -e "\\n-- Starting BUSCO --\\n"

            busco -i ${sample_id}.Trinity.fa -o ${sample_id}.Trinity.bus4 -l ${params.busco4db} -m tran -c ${task.cpus} --offline

            echo -e "\\n-- DONE with BUSCO --\\n"

            cp ${sample_id}.Trinity.bus4/short_summary.*.${sample_id}.Trinity.bus4.txt .
            """
    }

    busco3_sum_OAS = Channel.create()
    busco3_summary_OAS.mix(busco3_ch_trinity_sum_OAS).groupTuple(by:0,size:2).into(busco3_sum_OAS)

    process summary_busco3_individual_OAS {

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/stats", mode: "copy", overwrite: true

        input:
            tuple sample_id, file(files) from busco3_sum_OAS

        output:
            tuple sample_id, file("${sample_id}.sum_busco3.txt") into final_sum_2v3_OAS

        script:
            """
            tri=\$( echo $files | tr " " "\\n" | grep ".Trinity.bus3.txt" )
            trans=\$( echo $files | tr " " "\\n" | grep ".TransPi.bus3.txt" )
            #Summary of BUSCO scores for the final_assemblies
            echo -e "Summary of BUSCO V3 \n" >>${sample_id}.sum_busco3.txt
            echo "-- TransPi BUSCO V3 scores -- " >>${sample_id}.sum_busco3.txt
            cat \${trans} >>${sample_id}.sum_busco3.txt
            echo -e "\\n-- Trinity BUSCO V3 scores --" >>${sample_id}.sum_busco3.txt
            cat \${tri} >>${sample_id}.sum_busco3.txt
            """
    }

    busco4_sum_OAS = Channel.create()
    busco4_summary_OAS.mix(busco4_ch_trinity_sum_OAS).groupTuple(by:0,size:2).into(busco4_sum_OAS)

    process summary_busco4_individual_OAS {

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/stats", mode: "copy", overwrite: true

        input:
            tuple sample_id, file(files) from busco4_sum_OAS

        output:
            tuple sample_id, file("${sample_id}.sum_busco4.txt") into final_sum_2v4_OAS

        script:
            """
            tri=\$( echo $files | tr " " "\\n" | grep ".Trinity.bus4.txt" )
            trans=\$( echo $files | tr " " "\\n" | grep ".TransPi.bus4.txt" )
            #Summary of BUSCO scores for the final_assemblies
            echo -e "Summary of BUSCO V4 \n" >>${sample_id}.sum_busco4.txt
            echo "-- TransPi BUSCO V4 scores -- " >>${sample_id}.sum_busco4.txt
            cat \${trans} >>${sample_id}.sum_busco4.txt
            echo -e "\\n-- Trinity BUSCO V4 scores --" >>${sample_id}.sum_busco4.txt
            cat \${tri} >>${sample_id}.sum_busco4.txt
            """
    }

    busco3_comp_OAS = Channel.create()
    busco3_comp_1_OAS.mix(busco3_comp_2_OAS).groupTuple(by:0,size:2).into(busco3_comp_OAS)

    process get_busco3_comparison_OAS {

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/figures/BUSCO3", mode: "copy", overwrite: true

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge -c bioconda r-reshape2=1.4.4 r-plotly=4.9.2.1 plotly-orca=3.4.2 r-ggplot2=3.3.0 r-svglite=1.2.3 r-ggthemes=4.2.0 r-knitr=1.29 r-rmarkdown=2.3 r-kableextra=1.1.0" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/mulled-v2-3f431f5f8e54df68ea0029c209fce3b154f6e186:94cad00b5306639ceab6aaf211f45740560abb90-0" : "quay.io/biocontainers/mulled-v2-3f431f5f8e54df68ea0029c209fce3b154f6e186:94cad00b5306639ceab6aaf211f45740560abb90-0")
        }

        input:
            tuple sample_id, file(files) from busco3_comp_OAS

        output:
            tuple sample_id, file("${sample_id}_BUSCO3_comparison.pdf"), file("${sample_id}_BUSCO3_comparison.svg") into busco3_fig_OAS
            tuple sample_id, file("*.csv") into busco3_OAS_csv

        script:
            """
            set +e
            tri=\$( echo $files | tr " " "\\n" | grep ".Trinity.bus3.txt" )
            trans=\$( echo $files | tr " " "\\n" | grep ".TransPi.bus3.txt" )
            bash get_busco_val.sh \${tri} \${trans} v3 ${sample_id}
            cp ${params.pipeInstall}/bin/busco_comparison.R .
            a=\$( cat final_spec )
            sed -i "s/MYSPEC/\${a}/" busco_comparison.R
            b=\$( cat final_perc )
            sed -i "s/MYPERC/\${b}/" busco_comparison.R
            c=\$( cat final_num )
            sed -i "s/MYVAL/\${c}/" busco_comparison.R
            Rscript busco_comparison.R ${sample_id}
            mv ${sample_id}_BUSCO_comparison.pdf ${sample_id}_BUSCO3_comparison.pdf
            mv ${sample_id}_BUSCO_comparison.svg ${sample_id}_BUSCO3_comparison.svg
            # csv
            sed -i 's/\$/\\n/g' final_*
            cat final_spec final_perc final_num | tr -d "'" >${sample_id}_busco3.csv
            """
    }

    busco4_comp_OAS = Channel.create()
    busco4_comp_1_OAS.mix(busco4_comp_2_OAS).groupTuple(by:0,size:2).into(busco4_comp_OAS)

    process get_busco4_comparison_OAS {

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/figures/BUSCO4", mode: "copy", overwrite: true

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge -c bioconda r-reshape2=1.4.4 r-plotly=4.9.2.1 plotly-orca=3.4.2 r-ggplot2=3.3.0 r-svglite=1.2.3 r-ggthemes=4.2.0 r-knitr=1.29 r-rmarkdown=2.3 r-kableextra=1.1.0" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/mulled-v2-3f431f5f8e54df68ea0029c209fce3b154f6e186:94cad00b5306639ceab6aaf211f45740560abb90-0" : "quay.io/biocontainers/mulled-v2-3f431f5f8e54df68ea0029c209fce3b154f6e186:94cad00b5306639ceab6aaf211f45740560abb90-0")
        }

        input:
            tuple sample_id, file(files) from busco4_comp_OAS

        output:
            tuple sample_id, file("${sample_id}_BUSCO4_comparison.pdf"), file("${sample_id}_BUSCO4_comparison.svg") into busco4_fig_OAS
            tuple sample_id, file("*.csv") into busco4_OAS_csv

        script:
            """
            set +e
            tri=\$( echo $files | tr " " "\\n" | grep ".Trinity.bus4.txt" )
            trans=\$( echo $files | tr " " "\\n" | grep ".TransPi.bus4.txt" )
            bash get_busco_val.sh \${tri} \${trans} v4 ${sample_id}
            cp ${params.pipeInstall}/bin/busco_comparison.R .
            a=\$( cat final_spec )
            sed -i "s/MYSPEC/\${a}/" busco_comparison.R
            b=\$( cat final_perc )
            sed -i "s/MYPERC/\${b}/" busco_comparison.R
            c=\$( cat final_num )
            sed -i "s/MYVAL/\${c}/" busco_comparison.R
            Rscript busco_comparison.R ${sample_id}
            mv ${sample_id}_BUSCO_comparison.pdf ${sample_id}_BUSCO4_comparison.pdf
            mv ${sample_id}_BUSCO_comparison.svg ${sample_id}_BUSCO4_comparison.svg
            # csv
            sed -i 's/\$/\\n/g' final_*
            cat final_spec final_perc final_num | tr -d "'" >${sample_id}_busco4.csv
            """
    }

} else if (params.onlyAnn) {

    println("\n\tRunning only annotation analysis\n")

    Channel
        .fromFilePairs("${launchDir}/onlyAnn/*.{fa,fasta}", size: -1, checkIfExists: true)
        .into{ annotation_ch_transdecoder_OA; annotation_ch_transdecoderB_OA; assembly_ch_rnammer_OA }

    if (params.shortTransdecoder) {

        process transdecoder_short_OA {

            label 'med_cpus'

            tag "${sample_id}"

            publishDir "${launchDir}/${params.outdir}/transdecoder", mode: "copy", overwrite: true

            conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::transdecoder=5.5.0=pl526_2" : null)
            if (params.oneContainer){ container "${params.TPcontainer}" } else {
            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/transdecoder:5.5.0--pl526_2" : "quay.io/biocontainers/transdecoder:5.5.0--pl526_2")
            }

            input:
                tuple sample_id, file(assembly) from annotation_ch_transdecoder_OA

            output:
                tuple sample_id, file("${assembly}"), file("${sample_id}.*.transdecoder.pep") into ( transdecoder_ch_diamond_OA, transdecoder_ch_diamond_custom_OA )
                tuple sample_id, file("${sample_id}.*.transdecoder.pep") into ( transdecoder_ch_trinotate_OA, transdecoder_ch_hmmer_OA, transdecoder_ch_signalp_OA, transdecoder_ch_tmhmm_OA )
                tuple sample_id, file("${assembly}") into transdecoder_assembly_ch_trinotate_OA
                tuple sample_id, file("${sample_id}_transdecoder.stats") into transdecoder_summary_OA
                tuple sample_id, file("${sample_id}_transdecoder.csv") into transdecoder_csv_OA
                tuple sample_id, file("${sample_id}.*.transdecoder.{cds,gff,bed}") into transdecoder_files_OA

            script:
                """
                echo -e "\\n-- TransDecoder.LongOrfs... --\\n"

                TransDecoder.LongOrfs -t ${assembly} --output_dir ${sample_id}.transdecoder_dir

                echo -e "\\n-- Done with TransDecoder.LongOrfs --\\n"

                echo -e "\\n-- TransDecoder.Predict... --\\n"

                TransDecoder.Predict -t ${assembly} --output_dir ${sample_id}.transdecoder_dir

                echo -e "\\n-- Done with TransDecoder.Predict --\\n"

                echo -e "\\n-- Calculating statistics... --\\n"
                #Calculate statistics of Transdecoder
                echo "- Transdecoder (short,no homolgy) stats for ${sample_id}" >${sample_id}_transdecoder.stats
                orfnum=\$( cat ${sample_id}.*.transdecoder.pep | grep -c ">" )
                echo -e "Total number of ORFs: \$orfnum \\n" >>${sample_id}_transdecoder.stats
                orfnum=\$( cat ${sample_id}.*.transdecoder.pep | grep -c "ORF type:complete" )
                echo -e "\\t ORFs type=complete: \$orfnum \\n" >>${sample_id}_transdecoder.stats
                orfnum=\$( cat ${sample_id}.*.transdecoder.pep | grep -c "ORF type:5prime_partial" )
                echo -e "\\t ORFs type=5prime_partial: \$orfnum \\n" >>${sample_id}_transdecoder.stats
                orfnum=\$( cat ${sample_id}.*.transdecoder.pep | grep -c "ORF type:3prime_partial" )
                echo -e "\\t ORFs type=3prime_partial: \$orfnum \\n" >>${sample_id}_transdecoder.stats
                orfnum=\$( cat ${sample_id}.*.transdecoder.pep | grep -c "ORF type:internal" )
                echo -e "\\t ORFs type=internal: \$orfnum \\n">>${sample_id}_transdecoder.stats
                # csv for report
                echo "Sample,Total_orf,orf_complete,orf_5prime_partial,orf_3prime_partial,orf_internal" >${sample_id}_transdecoder.csv
                total=\$( cat ${sample_id}.*.transdecoder.pep  | grep -c ">" )
                complete=\$( cat ${sample_id}.*.transdecoder.pep  | grep -c "ORF type:complete" )
                n5prime=\$( cat ${sample_id}.*.transdecoder.pep  | grep -c "ORF type:5prime_partial" )
                n3prime=\$( cat ${sample_id}.*.transdecoder.pep  | grep -c "ORF type:3prime_partial" )
                internal=\$( cat ${sample_id}.*.transdecoder.pep  | grep -c "ORF type:internal" )
                echo "${sample_id},\${total},\${complete},\${n5prime},\${n3prime},\${internal}" >>${sample_id}_transdecoder.csv
                echo -e "\\n-- Done with statistics --\\n"

                cp ${assembly} ${assembly}.tmp
                rm ${assembly}
                mv ${assembly}.tmp ${assembly}

                echo -e "\\n-- DONE with TransDecoder --\\n"
                """
            }
        } else {

            process transdecoder_longorf_OA {

                label 'med_cpus'

                tag "${sample_id}"

                publishDir "${launchDir}/${params.outdir}/transdecoder", mode: "copy", overwrite: true

                conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::transdecoder=5.5.0=pl526_2" : null)
                if (params.oneContainer){ container "${params.TPcontainer}" } else {
                container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/transdecoder:5.5.0--pl526_2" : "quay.io/biocontainers/transdecoder:5.5.0--pl526_2")
                }

                input:
                    tuple sample_id, file(assembly) from annotation_ch_transdecoder_OA

                output:
                    tuple sample_id, file("${sample_id}.longest_orfs.pep") into transdecoder_diamond_OA, transdecoder_hmmer_OA

                script:
                    """
                    cp ${assembly} ${sample_id}_asssembly.fasta

                    echo -e "\\n-- TransDecoder.LongOrfs... --\\n"

                    TransDecoder.LongOrfs -t ${assembly} --output_dir ${sample_id}.transdecoder_dir

                    cp ${sample_id}.transdecoder_dir/longest_orfs.pep ${sample_id}.longest_orfs.pep

                    echo -e "\\n-- Done with TransDecoder.LongOrfs --\\n"
                    """
            }

            process transdecoder_diamond_OA {

                label 'med_cpus'

                tag "${sample_id}"

                publishDir "${launchDir}/${params.outdir}/transdecoder", mode: "copy", overwrite: true

                conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::diamond=0.9.30=h56fc30b_0" : null)
                if (params.oneContainer){ container "${params.TPcontainer}" } else {
                container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/diamond:0.9.30--h56fc30b_0" : "quay.io/biocontainers/diamond:0.9.30--h56fc30b_0")
                }

                input:
                    tuple sample_id, file(pep) from transdecoder_diamond_OA

                output:
                    tuple sample_id, file("${sample_id}.diamond_blastp.outfmt6") into transdecoder_predict_diamond_OA

                script:
                    """
                    dbPATH=${params.pipeInstall}/DBs/diamonddb_custom/

                    echo -e "\\n-- Starting Diamond (blastp) --\\n"
                    if [ ! -d \${dbPATH} ];then
                        echo "Directory \${dbPATH} not found. Run the precheck to fix this issue"
                        exit 0
                    elif [ -d \${dbPATH} ];then
                        if [ ! -e \${dbPATH}/${params.uniname} ];then
                            echo "File \${dbPATH}/${params.uniname} not found. Run the precheck to fix this issue"
                            exit 0
                        elif [ -e \${dbPATH}/${params.uniname} ];then
                            if [ ! -e \${dbPATH}/${params.uniname}.dmnd ];then
                                cp \${dbPATH}/${params.uniname} .
                                diamond makedb --in ${params.uniname} -d ${params.uniname} -p ${task.cpus}
                                diamond blastp -d ${params.uniname}.dmnd -q ${pep} -p ${task.cpus} -f 6 -k 1 -e 0.00001 >${sample_id}.diamond_blastp.outfmt6
                            elif [ -e \${dbPATH}/${params.uniname}.dmnd ];then
                                cp \${dbPATH}/${params.uniname}.dmnd .
                                diamond blastp -d ${params.uniname}.dmnd -q ${pep} -p ${task.cpus} -f 6 -k 1 -e 0.00001 >${sample_id}.diamond_blastp.outfmt6
                            fi
                        fi
                    fi
                    echo -e "\\n-- Done with Diamond (blastp) --\\n"
                    """
            }

            process transdecoder_hmmer_OA {

                label 'med_cpus'

                tag "${sample_id}"

                publishDir "${launchDir}/${params.outdir}/transdecoder", mode: "copy", overwrite: true

                conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::hmmer=3.3=he1b5a44_0" : null)
                if (params.oneContainer){ container "${params.TPcontainer}" } else {
                container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/hmmer:3.3--he1b5a44_0" : "quay.io/biocontainers/hmmer:3.3--he1b5a44_0")
                }

                input:
                    tuple sample_id, file(pep) from transdecoder_hmmer_OA

                output:
                    tuple sample_id, file("${sample_id}.pfam.domtblout") into transdecoder_predict_hmmer_OA

                script:
                    """
                    dbPATH=${params.pipeInstall}/DBs/hmmerdb/
                    echo -e "\\n-- Starting HMMER --\\n"
                    if [ ! -d \${dbPATH} ];then
                        echo "Directory \${dbPATH} not found. Run the precheck to fix this issue"
                        exit 0
                    elif [ -d \${dbPATH} ];then
                        if [ ! -e \${dbPATH}/${params.pfname} ];then
                            echo "File \${dbPATH}/${params.pfname} not found. Run the precheck to fix this issue"
                            exit 0
                        elif [ -e \${dbPATH}/${params.pfname} ];then
                            if [ ! -e \${dbPATH}/${params.pfname}.h3f ] && [ ! -e \${dbPATH}/${params.pfname}.h3i ] && [ ! -e \${dbPATH}/${params.pfname}.h3m ] && [ ! -e \${dbPATH}/${params.pfname}.h3p ];then
                                cp \${dbPATH}/${params.pfname} .
                                hmmpress ${params.pfname}
                                hmmscan --cpu ${task.cpus} --domtblout ${sample_id}.pfam.domtblout ${params.pfname} ${pep}
                            elif [ -s \${dbPATH}/${params.pfname}.h3f ] && [ -s \${dbPATH}/${params.pfname}.h3i ] && [ -s \${dbPATH}/${params.pfname}.h3m ] && [ -s \${dbPATH}/${params.pfname}.h3p ];then
                                cp \${dbPATH}/${params.pfname}.* .
                                hmmscan --cpu ${task.cpus} --domtblout ${sample_id}.pfam.domtblout ${params.pfname} ${pep}
                            else
                                cp \${dbPATH}/${params.pfname} .
                                hmmpress ${params.pfname}
                                hmmscan --cpu ${task.cpus} --domtblout ${sample_id}.pfam.domtblout ${params.pfname} ${pep}
                            fi
                        fi
                    fi
                    echo -e "\\n-- Done with HMMER --\\n"
                    """
            }

            transdecoder_predict_OA_ch=Channel.create()
            annotation_ch_transdecoderB_OA.mix(transdecoder_predict_diamond_OA,transdecoder_predict_hmmer_OA).groupTuple(by:0,size:3).into(transdecoder_predict_OA_ch)

            process transdecoder_predict_OA {

                label 'med_cpus'

                tag "${sample_id}"

                publishDir "${launchDir}/${params.outdir}/transdecoder", mode: "copy", overwrite: true

                conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::transdecoder=5.5.0=pl526_2" : null)
                if (params.oneContainer){ container "${params.TPcontainer}" } else {
                container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/transdecoder:5.5.0--pl526_2" : "quay.io/biocontainers/transdecoder:5.5.0--pl526_2")
                }

                input:
                    tuple sample_id, file(files) from transdecoder_predict_OA_ch

                output:
                    tuple sample_id, file("${assembly}"), file("${sample_id}*.transdecoder.pep") into ( transdecoder_ch_diamond_OA, transdecoder_ch_diamond_custom_OA )
                    tuple sample_id, file("${sample_id}*.transdecoder.pep") into ( transdecoder_ch_trinotate_OA, transdecoder_ch_hmmer_OA, transdecoder_ch_signalp_OA, transdecoder_ch_tmhmm_OA )
                    tuple sample_id, file("${sample_id}_asssembly.fasta") into transdecoder_assembly_ch_trinotate_OA
                    tuple sample_id, file("${sample_id}_transdecoder.stats") into transdecoder_summary_OA
                    tuple sample_id, file("${sample_id}_transdecoder.csv") into transdecoder_csv_OA
                    tuple sample_id, file("${sample_id}*.transdecoder.{cds,gff,bed}") into transdecoder_files_OA

                script:
                    """
                    ass=\$( echo $files | tr " " "\\n" | grep -v "pfam.domtblout" | grep ".fa" )
                    dia=\$( echo $files | tr " " "\\n" | grep ".diamond_blastp.outfmt6" )
                    pfa=\$( echo $files | tr " " "\\n" | grep ".pfam.domtblout" )

                    echo -e "\\n-- TransDecoder.LongOrfs... --\\n"

                    TransDecoder.LongOrfs -t \${ass} --output_dir ${sample_id}.transdecoder_dir

                    echo -e "\\n-- Done with TransDecoder.LongOrfs --\\n"

                    echo -e "\\n-- TransDecoder.Predict... --\\n"

                    TransDecoder.Predict -t \${ass} --retain_pfam_hits \${pfa} --retain_blastp_hits \${dia} --output_dir ${sample_id}.transdecoder_dir

                    echo -e "\\n-- Done with TransDecoder.Predict --\\n"

                    echo -e "\\n-- Calculating statistics... --\\n"
                    #Calculate statistics of Transdecoder
                    echo "- Transdecoder (long, with homology) stats for ${sample_id}" >${sample_id}_transdecoder.stats
                    orfnum=\$( cat ${sample_id}*.transdecoder.pep | grep -c ">" )
                    echo -e "Total number of ORFs: \$orfnum \\n" >>${sample_id}_transdecoder.stats
                    echo -e "\\t Of these ORFs" >>${sample_id}_transdecoder.stats
                    orfnum=\$( cat ${sample_id}*.transdecoder.pep | grep ">" | grep -c "|" )
                    echo -e "\\t\\t with annotations: \$orfnum \\n" >>${sample_id}_transdecoder.stats
                    orfnum=\$( cat ${sample_id}*.transdecoder.pep | grep ">" | grep -v "|" | grep -c ">" )
                    echo -e "\\t\\t no annotation: \$orfnum \\n" >>${sample_id}_transdecoder.stats
                    orfnum=\$( cat ${sample_id}*.transdecoder.pep | grep -c "ORF type:complete" )
                    echo -e "\\t ORFs type=complete: \$orfnum \\n" >>${sample_id}_transdecoder.stats
                    orfnum=\$( cat ${sample_id}*.transdecoder.pep | grep "ORF type:complete" | grep -c "|" )
                    echo -e "\\t\\t with annotations: \$orfnum \\n" >>${sample_id}_transdecoder.stats
                    orfnum=\$( cat ${sample_id}*.transdecoder.pep | grep -c "ORF type:5prime_partial" )
                    echo -e "\\t ORFs type=5prime_partial: \$orfnum \\n" >>${sample_id}_transdecoder.stats
                    orfnum=\$( cat ${sample_id}*.transdecoder.pep | grep "ORF type:5prime_partial" | grep -c "|" )
                    echo -e "\\t\\t with annotations: \$orfnum \\n" >>${sample_id}_transdecoder.stats
                    orfnum=\$( cat ${sample_id}*.transdecoder.pep | grep -c "ORF type:3prime_partial" )
                    echo -e "\\t ORFs type=3prime_partial: \$orfnum \\n" >>${sample_id}_transdecoder.stats
                    orfnum=\$( cat ${sample_id}*.transdecoder.pep | grep "ORF type:3prime_partial" | grep -c "|" )
                    echo -e "\\t\\t with annotations: \$orfnum \\n" >>${sample_id}_transdecoder.stats
                    orfnum=\$( cat ${sample_id}*.transdecoder.pep | grep -c "ORF type:internal" )
                    echo -e "\\t ORFs type=internal: \$orfnum \\n" >>${sample_id}_transdecoder.stats
                    orfnum=\$( cat ${sample_id}*.transdecoder.pep | grep "ORF type:internal" | grep -c "|" )
                    echo -e "\\t\\t with annotations: \$orfnum \\n" >>${sample_id}_transdecoder.stats
                    # csv for report
                    echo "Sample,Total_orf,orf_complete,orf_5prime_partial,orf_3prime_partial,orf_internal" >${sample_id}_transdecoder.csv
                    total=\$( cat ${sample_id}*.transdecoder.pep  | grep -c ">" )
                    complete=\$( cat ${sample_id}*.transdecoder.pep  | grep -c "ORF type:complete" )
                    n5prime=\$( cat ${sample_id}*.transdecoder.pep  | grep -c "ORF type:5prime_partial" )
                    n3prime=\$( cat ${sample_id}*.transdecoder.pep  | grep -c "ORF type:3prime_partial" )
                    internal=\$( cat ${sample_id}*.transdecoder.pep  | grep -c "ORF type:internal" )
                    echo "${sample_id},\${total},\${complete},\${n5prime},\${n3prime},\${internal}" >>${sample_id}_transdecoder.csv
                    echo -e "\\n-- Done with statistics --\\n"

                    cp ${assembly} ${assembly}.tmp
                    rm ${assembly}
                    mv ${assembly}.tmp ${assembly}

                    echo -e "\\n-- DONE with TransDecoder --\\n"
                    """
            }
        }

    process swiss_diamond_trinotate_OA {

        label 'big_cpus'

        tag "${sample_id}"

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::diamond=0.9.30=h56fc30b_0" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/diamond:0.9.30--h56fc30b_0" : "quay.io/biocontainers/diamond:0.9.30--h56fc30b_0")
        }

        input:
            tuple sample_id, file(assembly), file(transdecoder) from transdecoder_ch_diamond_OA

        output:
            tuple sample_id, file("${sample_id}.diamond_blastx.outfmt6") into trinotate_ch_diamondX_OA
            tuple sample_id, file("${sample_id}.diamond_blastp.outfmt6") into trinotate_ch_diamondP_OA

        script:
            """
            dbPATH=${params.pipeInstall}/DBs/diamonddb_swiss/

            echo -e "\\n-- Starting Diamond --\\n"
            if [ ! -d \${dbPATH} ];then
                echo "Directory \${dbPATH} not found. Run the precheck to fix this issue"
                exit 0
            elif [ -d \${dbPATH} ];then
                if [ ! -e \${dbPATH}/uniprot_sprot.pep ];then
                    echo "File \${dbPATH}/uniprot_sprot.pep not found. Run the precheck to fix this issue"
                    exit 0
                elif [ -e \${dbPATH}/uniprot_sprot.pep ];then
                    if [ ! -e \${dbPATH}/uniprot_sprot.pep.dmnd ];then
                        cp \${dbPATH}/uniprot_sprot.pep .
                        diamond makedb --in uniprot_sprot.pep -d uniprot_sprot.pep -p ${task.cpus}
                        echo -e "\\n-- Starting with Diamond (blastx) --\\n"
                        diamond blastx -d uniprot_sprot.pep.dmnd -q ${assembly} -p ${task.cpus} -f 6 -k 1 -e 0.001 >${sample_id}.diamond_blastx.outfmt6
                        echo -e "\\n-- Done with Diamond (blastx) --\\n"
                        echo -e "\\n-- Starting with Diamond (blastp) --\\n"
                        diamond blastp -d uniprot_sprot.pep.dmnd -q ${transdecoder} -p ${task.cpus} -f 6 -k 1 -e 0.001 >${sample_id}.diamond_blastp.outfmt6
                        echo -e "\\n-- Done with Diamond (blastp)  --\\n"
                    elif [ -e \${dbPATH}/uniprot_sprot.pep.dmnd ];then
                        cp \${dbPATH}/uniprot_sprot.pep.dmnd .
                        echo -e "\\n-- Starting with Diamond (blastx) --\\n"
                        diamond blastx -d uniprot_sprot.pep.dmnd -q ${assembly} -p ${task.cpus} -f 6 -k 1 -e 0.001 >${sample_id}.diamond_blastx.outfmt6
                        echo -e "\\n-- Done with Diamond (blastx) --\\n"
                        echo -e "\\n-- Starting with Diamond (blastp) --\\n"
                        diamond blastp -d uniprot_sprot.pep.dmnd -q ${transdecoder} -p ${task.cpus} -f 6 -k 1 -e 0.001 >${sample_id}.diamond_blastp.outfmt6
                        echo -e "\\n-- Done with Diamond (blastp)  --\\n"
                    fi
                fi
            fi
            echo -e "\\n-- Done with Diamond --\\n"
            """
    }

    process custom_diamond_trinotate_OA {

        label 'big_cpus'

        tag "${sample_id}"

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::diamond=0.9.30=h56fc30b_0" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/diamond:0.9.30--h56fc30b_0" : "quay.io/biocontainers/diamond:0.9.30--h56fc30b_0")
        }

        input:
            tuple sample_id, file(assembly), file(transdecoder) from transdecoder_ch_diamond_custom_OA

        output:
            tuple sample_id, file("${sample_id}.custom.diamond_blastx.outfmt6") into trinotate_ch_diamondX_custom_OA
            tuple sample_id, file("${sample_id}.custom.diamond_blastp.outfmt6") into trinotate_ch_diamondP_custom_OA

        script:
            """
            dbPATH=${params.pipeInstall}/DBs/diamonddb_custom/

            echo -e "\\n-- Starting Diamond --\\n"
            if [ ! -d \${dbPATH} ];then
                echo "Directory \${dbPATH} not found. Run the precheck to fix this issue"
                exit 0
            elif [ -d \${dbPATH} ];then
                if [ ! -e \${dbPATH}/${params.uniname} ];then
                    echo "File \${dbPATH}/${params.uniname} not found. Run the precheck to fix this issue"
                    exit 0
                elif [ -e \${dbPATH}/${params.uniname} ];then
                    if [ ! -e \${dbPATH}/${params.uniname}.dmnd ];then
                        cp \${dbPATH}/${params.uniname} .
                        diamond makedb --in ${params.uniname} -d ${params.uniname} -p ${task.cpus}
                        echo -e "\\n-- Starting with Diamond (blastx) --\\n"
                        diamond blastx -d ${params.uniname}.dmnd -q ${assembly} -p ${task.cpus} -f 6 -k 1 -e 0.001 >${sample_id}.custom.diamond_blastx.outfmt6
                        echo -e "\\n-- Done with Diamond (blastx) --\\n"
                        echo -e "\\n-- Starting with Diamond (blastp) --\\n"
                        diamond blastp -d ${params.uniname}.dmnd -q ${transdecoder} -p ${task.cpus} -f 6 -k 1 -e 0.001 >${sample_id}.custom.diamond_blastp.outfmt6
                        echo -e "\\n-- Done with Diamond (blastp)  --\\n"
                    elif [ -e \${dbPATH}/${params.uniname}.dmnd ];then
                        cp \${dbPATH}/${params.uniname}.dmnd .
                        echo -e "\\n-- Starting with Diamond (blastx) --\\n"
                        diamond blastx -d ${params.uniname}.dmnd -q ${assembly} -p ${task.cpus} -f 6 -k 1 -e 0.001 >${sample_id}.custom.diamond_blastx.outfmt6
                        echo -e "\\n-- Done with Diamond (blastx) --\\n"
                        echo -e "\\n-- Starting with Diamond (blastp) --\\n"
                        diamond blastp -d ${params.uniname}.dmnd -q ${transdecoder} -p ${task.cpus} -f 6 -k 1 -e 0.001 >${sample_id}.custom.diamond_blastp.outfmt6
                        echo -e "\\n-- Done with Diamond (blastp)  --\\n"
                    fi
                fi
            fi
            echo -e "\\n-- Done with Diamond --\\n"
            """
    }

    process hmmer_trinotate_OA {

        label 'low_cpus'

        tag "${sample_id}"

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::hmmer=3.3=he1b5a44_0" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/hmmer:3.3--he1b5a44_0" : "quay.io/biocontainers/hmmer:3.3--he1b5a44_0")
        }

        input:
            tuple sample_id, file(transdecoder_pep) from transdecoder_ch_hmmer_OA

        output:
            tuple sample_id, file("${sample_id}.TrinotatePFAM.out") into trinotate_ch_hmmer_OA

        script:
            """
            dbPATH=${params.pipeInstall}/DBs/hmmerdb/

            echo -e "\\n-- Starting HMMER --\\n"
            if [ ! -d \${dbPATH} ];then
                echo "Directory \${dbPATH} not found. Run the precheck to fix this issue"
                exit 0
            elif [ -d \${dbPATH} ];then
                if [ ! -e \${dbPATH}/${params.pfname} ];then
                    echo "File \${dbPATH}/${params.pfname} not found. Run the precheck to fix this issue"
                    exit 0
                elif [ -e \${dbPATH}/${params.pfname} ];then
                    if [ ! -e \${dbPATH}/${params.pfname}.h3f ] && [ ! -e \${dbPATH}/${params.pfname}.h3i ] && [ ! -e \${dbPATH}/${params.pfname}.h3m ] && [ ! -e \${dbPATH}/${params.pfname}.h3p ];then
                        cp \${dbPATH}/${params.pfname} .
                        hmmpress ${params.pfname}
                        hmmscan --cpu ${task.cpus} --domtblout ${sample_id}.TrinotatePFAM.out ${params.pfname} ${transdecoder_pep} >pfam.log
                    elif [ -s \${dbPATH}/${params.pfname}.h3f ] && [ -s \${dbPATH}/${params.pfname}.h3i ] && [ -s \${dbPATH}/${params.pfname}.h3m ] && [ -s \${dbPATH}/${params.pfname}.h3p ];then
                        cp \${dbPATH}/${params.pfname}.* .
                        hmmscan --cpu ${task.cpus} --domtblout ${sample_id}.TrinotatePFAM.out ${params.pfname} ${transdecoder_pep} >pfam.log
                    else
                        cp \${dbPATH}/${params.pfname} .
                        hmmpress ${params.pfname}
                        hmmscan --cpu ${task.cpus} --domtblout ${sample_id}.TrinotatePFAM.out ${params.pfname} ${transdecoder_pep} >pfam.log
                    fi
                fi
            fi
            echo -e "\\n-- Done with HMMER --\\n"
            """
    }

    if (params.withSignalP) {

        process signalP_trinotate_OA {

            label 'low_cpus'

            tag "${sample_id}"

            input:
                tuple sample_id, file(transdecoder_pep) from transdecoder_ch_signalp_OA

            output:
                tuple sample_id, file("${sample_id}.signalp.out") into trinotate_ch_signalp_OA

            script:
                """
                #signalP to predict signal peptides

                echo -e "\\n-- Starting with SignalP --\\n"

                ${params.signalp} -f short -n ${sample_id}.signalp.out ${transdecoder_pep}

                echo -e "\\n-- Done with SignalP --\\n"
                """
        }
    } else {

        process skip_signalP_OA {
            tag "${sample_id}"

            input:
                tuple sample_id, file(transdecoder_pep) from transdecoder_ch_signalp_OA

            output:
                tuple sample_id, file("${sample_id}.signalp.out") into trinotate_ch_signalp_OA

            script:
                """
                touch ${sample_id}.signalp.out
                """
        }
    }

    if (params.withTMHMM) {

        process tmhmm_trinotate_OA {

            label 'low_cpus'

            tag "${sample_id}"

            input:
                tuple sample_id, file(transdecoder_pep) from transdecoder_ch_tmhmm_OA

            output:
                tuple sample_id, file("${sample_id}.tmhmm.out") into trinotate_ch_tmhmm_OA

            script:
                """
                #tmHMM to predict transmembrane regions

                echo -e "\\n-- Starting with tmHMM --\\n"

                ${params.tmhmm} --short < ${transdecoder_pep} >${sample_id}.tmhmm.out

                echo -e "\\n-- Done with tmHMM --\\n"
                """
        }
    } else {

        process skip_tmhmm_OA {
            tag "${sample_id}"

            input:
                tuple sample_id, file(transdecoder_pep) from transdecoder_ch_tmhmm_OA

            output:
                tuple sample_id, file("${sample_id}.tmhmm.out") into trinotate_ch_tmhmm_OA

            script:
                """
                touch ${sample_id}.tmhmm.out
                """
        }
    }

    if (params.withRnammer) {

        process rnammer_trinotate_OA {

            label 'low_cpus'

            tag "${sample_id}"

            input:
                tuple sample_id, file(assembly) from assembly_ch_rnammer_OA

            output:
                tuple sample_id, file("${sample_id}.rnammer.gff") into trinotate_ch_rnammer_OA

            script:
                """
                set +e
                #RNAMMER to identify rRNA transcripts

                echo -e "\\n-- Starting with RNAMMER --\\n"

                RnammerTranscriptome.pl --transcriptome ${assembly} --path_to_rnammer ${params.rnam}

                mv ${assembly}*rnammer.gff ${sample_id}.rnammer.gff

                echo -e "\\n-- Done with RNAMMER --\\n"
                """
        }
    } else {

        process skip_rnammer_OA {
            tag "${sample_id}"

            input:
                tuple sample_id, file(assembly) from assembly_ch_rnammer_OA

            output:
                tuple sample_id, file("${sample_id}.rnammer.gff") into trinotate_ch_rnammer_OA

            script:
                """
                touch ${sample_id}.rnammer.gff
                """
        }
    }

    trinotate_ch = Channel.create()
    transdecoder_ch_trinotate_OA.mix( transdecoder_assembly_ch_trinotate_OA,trinotate_ch_diamondX_OA,trinotate_ch_diamondP_OA,trinotate_ch_diamondX_custom_OA,trinotate_ch_diamondP_custom_OA,trinotate_ch_hmmer_OA,trinotate_ch_signalp_OA,trinotate_ch_tmhmm_OA,trinotate_ch_rnammer_OA ).groupTuple(by:0,size:10).into(trinotate_ch)

    process trinotate_OA {

        label 'low_cpus'

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/trinotate", mode: "copy", overwrite: true

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::trinotate=3.2.1=pl526_0" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/trinotate:3.2.1--pl526_0" : "quay.io/biocontainers/trinotate:3.2.1--pl526_0")
        }

        input:
            tuple sample_id, file(files) from trinotate_ch

        output:
            tuple sample_id, file("${sample_id}.GO.terms.txt") into trinotate_summary_OA
            tuple sample_id, file("${sample_id}.trinotate_annotation_report.xls") into ( trinotate_ch_OA, custom_uniprot_ch_OA )
            tuple sample_id, file("*.terms.txt") into other_files_OA

        script:
            """
            for x in ${files};do
                echo \${x} >>.vars.txt
            done

            assembly=\$( cat .vars.txt | grep "${sample_id}_asssembly.fasta" | grep -v "transdecoder" )
            transdecoder=\$( cat .vars.txt | grep -E "${sample_id}.*.transdecoder.pep" )
            diamond_blastx=\$( cat .vars.txt | grep "${sample_id}.diamond_blastx.outfmt6" )
            diamond_blastp=\$( cat .vars.txt | grep "${sample_id}.diamond_blastp.outfmt6" )
            custom_blastx=\$( cat .vars.txt | grep "${sample_id}.custom.diamond_blastx.outfmt6" )
            custom_blastp=\$( cat .vars.txt | grep "${sample_id}.custom.diamond_blastp.outfmt6" )
            pfam=\$( cat .vars.txt | grep "${sample_id}.TrinotatePFAM.out" )
            signalp=\$( cat .vars.txt | grep "${sample_id}.signalp.out" )
            tmhmm=\$( cat .vars.txt | grep "${sample_id}.tmhmm.out" )
            rnammer=\$( cat .vars.txt | grep "${sample_id}.rnammer.gff" )

            #Generate gene_trans_map
            #Not using get_Trinity_gene_to_trans_map.pl since all the names are uniq
            cat \${assembly} | awk '{print \$1}' | grep ">" | cut -c 2- >a.txt

            paste a.txt a.txt >\${assembly}.gene_trans_map

            #Get Trinotate.sqlite from folder (original)
            cp ${params.Tsql} .
            sqlname=`echo ${params.Tsql} | tr "\\/" "\\n" | grep "\\.sqlite"`

            echo -e "\\n-- Running Trinotate --\\n"

            Trinotate \$sqlname init --gene_trans_map \${assembly}.gene_trans_map --transcript_fasta \${assembly} --transdecoder_pep \${transdecoder}

            echo -e "\\n-- Ending run of Trinotate --\\n"

            echo -e "\\n-- Loading hits and predictions to sqlite database... --\\n"

            #Load protein hits
            Trinotate \$sqlname LOAD_swissprot_blastp \${diamond_blastp}

            #Load transcript hits
            Trinotate \$sqlname LOAD_swissprot_blastx \${diamond_blastx}

            #Load custom protein hits
            Trinotate \$sqlname LOAD_custom_blast --outfmt6 \${custom_blastp} --prog blastp --dbtype ${params.uniname}

            #Load custom transcript hits
            Trinotate \$sqlname LOAD_custom_blast --outfmt6 \${custom_blastx} --prog blastx --dbtype ${params.uniname}

            #Load Pfam domain entries
            Trinotate \$sqlname LOAD_pfam \${pfam}

            #Load transmembrane domains
            if [ -s \${tmhmm} ];then
                Trinotate \$sqlname LOAD_tmhmm \${tmhmm}
            else
                echo "No transmembrane domains (tmhmm)"
            fi

            #Load signal peptide predictions
            if [ -s \${signalp} ];then
                Trinotate \$sqlname LOAD_signalp \${signalp}
            else
                echo "No Signal-P"
            fi

            #Load rnammer results
            if [ -s \${rnammer} ];then
                Trinotate \$sqlname LOAD_rnammer \${rnammer}
            else
                echo "No rnammer results"
            fi

            echo -e "\\n-- Loading finished --\\n"

            #Report

            echo -e "\\n-- Generating report... --\\n"

            Trinotate \$sqlname report >${sample_id}.trinotate_annotation_report.xls

            echo -e "\\n-- Report generated --\\n"

            #Extract info from XLS file

            echo -e "\\n-- Creating GO file from XLS... --\\n"

            extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls ${sample_id}.trinotate_annotation_report.xls --trans >${sample_id}.GO.terms.txt

            echo -e "\\n-- Done with the GO --\\n"

            echo -e "\\n-- Creating KEGG file from XLS... --\\n"

            cat ${sample_id}.trinotate_annotation_report.xls | cut -f 1,14 | grep "KEGG" | tr "\\`" ";" | grep "KO:K" | sed 's/\\tKEGG/\\t#KEGG/g' | sed 's/KO:/KO:#/g' | cut -f 1,3 -d "#" | tr -d "#" >${sample_id}.KEGG.terms.txt

            echo -e "\\n-- Done with the KEGG --\\n"

            echo -e "\\n-- Creating eggNOG file from XLS... --\\n"

            cat ${sample_id}.trinotate_annotation_report.xls | cut -f 1,13 | grep "OG" | tr "\\`" ";" | sed 's/^/#/g' | sed 's/;/\\n;/g' | cut -f 1 -d "^" | tr -d "\\n" | tr "#" "\\n" | grep "OG" >${sample_id}.eggNOG_COG.terms.txt

            echo -e "\\n-- Done with the eggNOG --\\n"

            echo -e "\\n-- Creating PFAM file from XLS... --\\n"

            cat ${sample_id}.trinotate_annotation_report.xls | cut -f 1,10 | grep "PF" | tr "\\`" ";" | sed 's/^/#/g' | sed 's/;PF/\\n;PF/g' | cut -f 1 -d "^" | tr -d "\\n" | tr "#" "\\n" | grep "PF" | tr ";" "," >${sample_id}.PFAM.terms.txt

            echo -e "\\n-- Done with the PFAM --\\n"

            echo -e "\\n-- DONE with Trinotate --\\n"
            """
    }

    process summary_transdecoder_individual_OA {

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/stats", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}_transdecoder.stats") from transdecoder_summary_OA

        output:
            tuple sample_id, file("${sample_id}.sum_transdecoder.txt") into final_sum_3_OA

        script:
            """
            #Summary of Transdecoder stats
            echo -e "Summary of Transdecoder \\n" >>${sample_id}.sum_transdecoder.txt
            cat ${sample_id}_transdecoder.stats >>${sample_id}.sum_transdecoder.txt
            echo -e "##### \\n" >>${sample_id}.sum_transdecoder.txt
            """
    }

    process summary_trinotate_individual_OA {

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/stats", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.GO.terms.txt") from trinotate_summary_OA

        output:
            tuple sample_id, file("${sample_id}.sum_GO.txt") into final_sum_4_OA

        script:
            """
            #Summary of Trinotate (Gene Ontologies)
            echo -e "Summary of Trinotate/Gene Ontologies \\n" >>${sample_id}.sum_GO.txt
            echo "- Individual "${sample_id} >>${sample_id}.sum_GO.txt
            echo -e "\\t Total transcripts with GO:" >>${sample_id}.sum_GO.txt
            num=\$( cat ${sample_id}.GO.terms.txt | wc -l )
            echo -e "\\t\\t \$num" >>${sample_id}.sum_GO.txt
            tnum=\$num

            echo -e "\\t Total transcripts with only one GO:" >>${sample_id}.sum_GO.txt
            a=0
            while read lines;do
                if [[ `echo \$lines | awk '{print \$2}' | tr "," "\\n" | wc -l` -eq 1 ]];then
                    a=\$((a+1))
                fi
            done <${sample_id}.GO.terms.txt
            num=\$a
            echo -e "\\t\\t \$num" >>${sample_id}.sum_GO.txt
            onum=\$num

            echo -e "\\t Total transcripts with multiple GO:" >>${sample_id}.sum_GO.txt
            num=\$( echo \$tnum-\$onum | bc )
            echo -e "\\t\\t \$num" >>${sample_id}.sum_GO.txt

            echo -e "\\t Total GO in the file:" >>${sample_id}.sum_GO.txt
            num=\$( cat ${sample_id}.GO.terms.txt | awk '{print \$2}' | tr "," "\\n" | wc -l )
            echo -e "\\t\\t \$num" >>${sample_id}.sum_GO.txt

            echo -e "\\t Total uniq GO in the file:" >>${sample_id}.sum_GO.txt
            num=\$( cat ${sample_id}.GO.terms.txt | awk '{print \$2}' | tr "," "\\n" | sort --parallel=10 | uniq | wc -l )
            echo -e "\\t\\t \$num" >>${sample_id}.sum_GO.txt
            """
    }

    process get_GO_comparison_OA {

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/figures/GO", mode: "copy", overwrite: true

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge -c bioconda r-reshape2=1.4.4 r-plotly=4.9.2.1 plotly-orca=3.4.2 r-ggplot2=3.3.0 r-svglite=1.2.3 r-ggthemes=4.2.0 r-knitr=1.29 r-rmarkdown=2.3 r-kableextra=1.1.0" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/mulled-v2-3f431f5f8e54df68ea0029c209fce3b154f6e186:94cad00b5306639ceab6aaf211f45740560abb90-0" : "quay.io/biocontainers/mulled-v2-3f431f5f8e54df68ea0029c209fce3b154f6e186:94cad00b5306639ceab6aaf211f45740560abb90-0")
        }

        input:
            tuple sample_id, file("${sample_id}.trinotate_annotation_report.xls") from trinotate_ch_OA

        output:
            tuple sample_id, file("*.svg"), file("*.pdf"), file("*.txt") into go_fig_OA
            tuple sample_id, file("*.csv") into go_csv

        script:
            """
            set +e
            cp ${params.pipeInstall}/bin/GO_plots.R .

            cat ${sample_id}.trinotate_annotation_report.xls | awk 'FS="\\t",OFS="#" {print \$1,\$15,\$16,\$17}' | grep -v "gene_id" >all_GOs.txt

            while read line;do
                echo \${line} | cut -f 2,3,4 -d "#" | grep "GO:" | tr "#" "\\n" | tr "\\`" "\\n" | sed 's/\\. /,/g' | tr "," "\\n" | grep "GO:" | sort -u >>final_GOs.txt
            done<all_GOs.txt

            cat final_GOs.txt | tr [a-z] [A-Z] | grep "CELLULAR_COMPONENT" | cut -f 3 -d "^" | sort | uniq -c | sort -nr | head -n 15 | sed 's/^ *//g' \
            | sed 's/\\([0-9] \\)/\\1#/g' | tr "#" "\\t" >GO_cellular.txt

            cat final_GOs.txt | tr [a-z] [A-Z] | grep "BIOLOGICAL_PROCESS" | cut -f 3 -d "^" | sort | uniq -c | sort -nr | head -n 15 | sed 's/^ *//g' \
            | sed 's/\\([0-9] \\)/\\1#/g' | tr "#" "\\t" >GO_biological.txt

            cat final_GOs.txt | tr [a-z] [A-Z] | grep "MOLECULAR_FUNCTION" | cut -f 3 -d "^" | sort | uniq -c | sort -nr | head -n 15 | sed 's/^ *//g' \
            | sed 's/\\([0-9] \\)/\\1#/g' | tr "#" "\\t" >GO_molecular.txt

            Rscript GO_plots.R ${sample_id}

            rm all_GOs.txt final_GOs.txt

            mv GO_cellular.txt ${sample_id}_GO_cellular.txt
            mv GO_biological.txt ${sample_id}_GO_biological.txt
            mv GO_molecular.txt ${sample_id}_GO_molecular.txt

            cat ${sample_id}_GO_cellular.txt | sed -r 's/^[^0-9]*([0-9]+)/\\1,/g' >${sample_id}_GO_cellular.csv
            cat ${sample_id}_GO_biological.txt | sed -r 's/^[^0-9]*([0-9]+)/\\1,/g' >${sample_id}_GO_biological.csv
            cat ${sample_id}_GO_molecular.txt | sed -r 's/^[^0-9]*([0-9]+)/\\1,/g' >${sample_id}_GO_molecular.csv
            """
    }

    process summary_custom_uniprot_OA {

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/figures/CustomUniProt", mode: "copy", overwrite: true

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge -c bioconda r-reshape2=1.4.4 r-plotly=4.9.2.1 plotly-orca=3.4.2 r-ggplot2=3.3.0 r-svglite=1.2.3 r-ggthemes=4.2.0 r-knitr=1.29 r-rmarkdown=2.3 r-kableextra=1.1.0" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/mulled-v2-3f431f5f8e54df68ea0029c209fce3b154f6e186:94cad00b5306639ceab6aaf211f45740560abb90-0" : "quay.io/biocontainers/mulled-v2-3f431f5f8e54df68ea0029c209fce3b154f6e186:94cad00b5306639ceab6aaf211f45740560abb90-0")
        }

        input:
            tuple sample_id, file("${sample_id}.trinotate_annotation_report.xls") from custom_uniprot_ch_OA

        output:
            tuple sample_id, file("${sample_id}_custom_uniprot_hits.txt") into custom_uniprot_sum_OA
            tuple sample_id, file("${sample_id}_custom_uniprot_hits.svg"), file("${sample_id}_custom_uniprot_hits.pdf") into custom_uniprot_fig_OA
            tuple sample_id, file("*.csv") into uniprot_OA_csv

        script:
            """
            #get custom blast hits
            cat ${sample_id}.trinotate_annotation_report.xls | cut -f 8 | grep [A-Z] | grep "|" | tr "\\`" "\n" | \
                cut -f 1 -d "^" | cut -f 3 -d "|" | cut -f 2 -d "_" >a.txt

            cat ${sample_id}.trinotate_annotation_report.xls | cut -f 9 | grep [A-Z] | grep "|" | tr "\\`" "\n" | \
                cut -f 1 -d "^" | cut -f 3 -d "|" | cut -f 2 -d "_" >b.txt

            cat a.txt b.txt | sort | uniq -c | sort -nr | head -n 20 | awk 'OFS="," {print \$1,\$2}' >${sample_id}_custom_uniprot_hits.txt

            rm a.txt b.txt

            cp ${params.pipeInstall}/conf/uni_tax.txt .
            cp ${sample_id}_custom_uniprot_hits.txt ${sample_id}_custom_uniprot_hits

            while read line;do
                a=\$( echo \${line} | cut -f 2 -d "," )
                b=\$( cat uni_tax.txt | grep "\${a}" | cut -f 2 -d "," | wc -l )
                if [ "\${b}" == "1" ];then
                    c=\$( cat uni_tax.txt | grep "\${a}" | cut -f 2 -d "," )
                    sed -i "s/\${a}/\${c}/" ${sample_id}_custom_uniprot_hits
                fi
            done <${sample_id}_custom_uniprot_hits.txt

            rm ${sample_id}_custom_uniprot_hits.txt uni_tax.txt
            mv ${sample_id}_custom_uniprot_hits ${sample_id}_custom_uniprot_hits.txt

            cp ${params.pipeInstall}/bin/custom_uniprot_hits.R .
            Rscript custom_uniprot_hits.R ${sample_id}

            cp ${sample_id}_custom_uniprot_hits.txt ${sample_id}_custom_uniprot_hits.csv
            """
    }

} else if (params.all) {

    println("\n\tRunning the full TransPi analysis\n")

    if (!params.skipQC) {

        process fasqc {

            label 'low_cpus'

            tag "${sample_id}"

            publishDir "${launchDir}/${params.outdir}/fastqc", mode: "copy", overwrite: true

            conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::fastqc=0.11.9=0" : null)
            if (params.oneContainer){ container "${params.TPcontainer}" } else {
            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0" : "quay.io/biocontainers/fastqc:0.11.9--0")
            }

            input:
                tuple sample_id, file(reads) from reads_qc_ch

            output:
                tuple sample_id, file("*_fastqc.{zip,html}") into fastqc_results

            script:
                """
                fastqc --quiet --threads $task.cpus $reads
                """
        }
    }

    // create channels for skip processes for the report
    skip_filter_ch = Channel.create()
    skip_filter_only_ch = Channel.create()
    skip_norm_ch = Channel.create()
    skip_busco_dist = Channel.create()
    report_reads.into( skip_filter_ch, skip_filter_only_ch, skip_norm_ch, skip_busco_dist )

    if (!params.skipFilter) {

        process fastp {

            label 'med_cpus'

            tag "${sample_id}"

            publishDir "${launchDir}/${params.outdir}/filter", mode: "copy", overwrite: true, pattern: "*.fastp.{json,html}"

            conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge conda-forge::pigz=2.3.4=hed695b0_1 conda-forge::jq=1.6=h14c3975_1000 bioconda::fastp=0.20.1=h8b12597_0" : null)
            if (params.oneContainer){ container "${params.TPcontainer}" } else {
            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/fastp:0.20.1--h8b12597_0" : "quay.io/biocontainers/fastp:0.20.1--h8b12597_0")
            }

            input:
                tuple sample_id, file(reads) from reads_ch

            output:
                tuple sample_id, file("*.fastp.{json,html}") into fastp_results
                tuple sample_id, file("${sample_id}.fastp.json") into fastp_stats_ch
                tuple sample_id, file("*${sample_id}.filter.fq") into reads_rna_ch
                tuple sample_id, file("left-${sample_id}.filter.fq"), file("right-${sample_id}.filter.fq") into save_filter_reads

            script:
                """
                echo ${sample_id}

                fastp -i ${reads[0]} -I ${reads[1]} -o left-${sample_id}.filter.fq -O right-${sample_id}.filter.fq --detect_adapter_for_pe \
                --average_qual ${params.minQual} --overrepresentation_analysis --html ${sample_id}.fastp.html --json ${sample_id}.fastp.json \
                --thread ${task.cpus} --report_title ${sample_id}
                """
        }

        process fastp_stats {

            label 'low_cpus'

            tag "${sample_id}"

            publishDir "${launchDir}/${params.outdir}/filter", mode: "copy", overwrite: true, pattern: "*.fastp.{json,html}"

            conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "conda-forge::jq=1.6=h14c3975_1000 " : null)
            if (params.oneContainer){ container "${params.TPcontainer}" } else {
            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/jq:1.6" : "quay.io/biocontainers/jq:1.6")
            }

            input:
                tuple sample_id, file(json) from fastp_stats_ch

            output:
                tuple sample_id, file("*.csv") into fastp_csv

            script:
                """
                echo ${sample_id}
                bash get_readstats.sh ${json}
                bash get_readqual.sh ${json}
                """

        }

        if (params.saveReads) {

            process save_filter_reads {

                label 'med_cpus'

                tag "${sample_id}"

                publishDir "${launchDir}/${params.outdir}/saveReads/filter", mode: "copy", overwrite: true, pattern: "*_R{1,2}.filter.fq.gz"

                conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "conda-forge::pigz=2.3.4=hed695b0_1" : null)
                if (params.oneContainer){ container "${params.TPcontainer}" } else {
                container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/pigz:2.3.4" : "quay.io/biocontainers/pigz:2.3.4")
                }

                input:
                    tuple sample_id, file(r1), file(r2) from save_filter_reads

                output:
                    tuple sample_id, file("*.filter.fq.gz") into save_filter_reads_out

                script:
                    """
                    cat $r1 >${sample_id}_filter.R1.fq
                    cat $r2 >${sample_id}_filter.R2.fq
                    pigz --best --force -p ${task.cpus} -r ${sample_id}_filter.R1.fq
                    pigz --best --force -p ${task.cpus} -r ${sample_id}_filter.R2.fq
                    """
            }
        }

    } else {
        reads_ch
            .set{ reads_rna_ch }
        fastp_results = Channel.empty()
    }

    if (params.rRNAfilter) {

        process remove_rrna {

            label 'med_mem'

            tag "${sample_id}"

            publishDir "${launchDir}/${params.outdir}/rRNA_reads", mode: "copy", overwrite: true, pattern: "*.log"

            conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "bioconda::sortmerna=4.2.0" : null)
            if (params.oneContainer){ container "${params.TPcontainer}" } else {
            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/sortmerna:4.2.0--0" : "quay.io/biocontainers/sortmerna:4.2.0--0")
            }

            input:
                tuple sample_id, file(reads) from reads_rna_ch

            output:
                tuple sample_id, file("*rRNA.R*.fq") into reads_rna_ass
                tuple sample_id, file("${sample_id}_rRNA_reads.R1.fq"), file("${sample_id}_rRNA_reads.R2.fq"), file("${sample_id}_no_rRNA.R1.fq"), file("${sample_id}_no_rRNA.R2.fq") into save_rrna_reads
                tuple sample_id, file("${sample_id}_remove_rRNA.log") into remove_rrna_sum

            script:
            if (!params.skipFilter) {
                """
                sortmerna --ref ${params.rRNAdb} --reads ${reads[0]} --reads ${reads[1]} --threads ${task.cpus} --aligned rRNAreads --other nonrRNAreads --paired_in --out2 --fastx --workdir .
                mv rRNAreads.log ${sample_id}_remove_rRNA.log
                mv rRNAreads_fwd* ${sample_id}_rRNA_reads.R1.fq
                mv rRNAreads_rev* ${sample_id}_rRNA_reads.R2.fq
                mv nonrRNAreads_fwd* ${sample_id}_no_rRNA.R1.fq
                mv nonrRNAreads_rev* ${sample_id}_no_rRNA.R2.fq
                """
            } else {
                """
                sortmerna --ref ${params.rRNAdb} --reads ${reads[0]} --reads ${reads[1]} --threads ${task.cpus} --aligned rRNAreads --other nonrRNAreads --paired_in --out2 --fastx --workdir .
                mv rRNAreads.log ${sample_id}_remove_rRNA.log
                mv rRNAreads_fwd* ${sample_id}_rRNA_reads.R1.fq
                mv rRNAreads_rev* ${sample_id}_rRNA_reads.R2.fq
                mv nonrRNAreads_fwd* ${sample_id}_no_rRNA.R1.fq
                mv nonrRNAreads_rev* ${sample_id}_no_rRNA.R2.fq
                """
            }
        }

        if (params.saveReads) {

            process save_rrna_reads {

                label 'med_cpus'

                tag "${sample_id}"

                publishDir "${launchDir}/${params.outdir}/saveReads/rRNA", mode: "copy", overwrite: true, pattern: "*.gz"

                conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "conda-forge::pigz=2.3.4=hed695b0_1" : null)
                if (params.oneContainer){ container "${params.TPcontainer}" } else {
                container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/pigz:2.3.4" : "quay.io/biocontainers/pigz:2.3.4")
                }

                input:
                    tuple sample_id, file(r1), file(r2), file(r3), file (r4) from save_rrna_reads

                output:
                    tuple sample_id, file("*.gz") into save_rrna_reads_out

                script:
                    """
                    cat $r1 >${sample_id}_rRNA_reads.R1.fq
                    pigz --best --force -p ${task.cpus} -r ${sample_id}_rRNA_reads.R1.fq
                    cat $r2 >${sample_id}_rRNA_reads.R2.fq
                    pigz --best --force -p ${task.cpus} -r ${sample_id}_rRNA_reads.R2.fq
                    cat $r3 >${sample_id}_no_rRNA.R1.fq
                    pigz --best --force -p ${task.cpus} -r ${sample_id}_no_rRNA.R1.fq
                    cat $r4 >${sample_id}_no_rRNA.R2.fq
                    pigz --best --force -p ${task.cpus} -r ${sample_id}_no_rRNA.R2.fq
                    """
            }
        }

    } else {
        reads_rna_ch
            .into{ reads_rna_ass; skip_rna_ch }

        process skip_rrna_removal {
            tag "${sample_id}"

            input:
                tuple sample_id, file(files) from skip_rna_ch

            output:
                tuple sample_id, file("rrna_removal.txt") into remove_rrna_sum

            script:
                """
                echo "rRNA removal step was skipped" >rrna_removal.txt
                """
        }
    }

    if (!params.skipNormalization) {

        process normalize_reads {

            label 'med_mem'

            tag "${sample_id}"

            conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge -c bioconda trinity=2.9.1=h8b12597_1" : null)
            if (params.oneContainer){ container "${params.TPcontainer}" } else {
            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/trinity:2.9.1--h8b12597_1" : "quay.io/biocontainers/trinity:2.9.1--h8b12597_1")
            }

            input:
                tuple sample_id, file(reads) from reads_rna_ass

            output:
                tuple sample_id, file("left-${sample_id}.norm.fq"), file("right-${sample_id}.norm.fq") into ( norm_reads_soap, norm_reads_velvet, norm_reads_trinity, norm_reads_spades, norm_reads_transabyss, reads_rna_quast )
                tuple sample_id, file("left-${sample_id}.norm.fq"), file("right-${sample_id}.norm.fq") into ( mapping_reads_trinity, mapping_reads_evi, mapping_symbiont )
                tuple sample_id, file("left-"${sample_id}".norm.fq"), file("right-"${sample_id}".norm.fq") into save_norm_reads
                tuple sample_id, file("${sample_id}_normStats.txt") into norm_report

            script:
                //def mem=(task.memory)
                //def mem_MB=(task.memory.toMega())
            if (!params.skipFilter) {
                """
                echo ${sample_id}

                echo -e "\\n-- Starting Normalization --\\n"

                mem=\$( echo ${task.memory} | cut -f 1 -d " " )

                insilico_read_normalization.pl --seqType fq -JM \${mem}G --max_cov ${params.normMaxCov} --min_cov ${params.normMinCov} --left ${reads[0]} --right ${reads[1]} --pairs_together --PARALLEL_STATS --CPU ${task.cpus}

                echo -e "\\n-- DONE with Normalization --\\n"

                cat .command.out | grep "stats_file" -A 3 | tail -n 3 >${sample_id}_normStats.txt

                cp left.norm.fq left-"${sample_id}".norm.fq
                cp right.norm.fq right-"${sample_id}".norm.fq

                rm \$(readlink -e left.norm.fq) \$(readlink -e right.norm.fq ) left.norm.fq right.norm.fq
                """
            } else if (params.skipFilter && !params.rRNAfilter) {
                """
                echo ${sample_id}
                zcat ${reads[0]} >left-${sample_id}.fq &
                zcat ${reads[1]} >right-${sample_id}.fq

                echo -e "\\n-- Starting Normalization --\\n"

                mem=\$( echo ${task.memory} | cut -f 1 -d " " )

                insilico_read_normalization.pl --seqType fq -JM \${mem}G --max_cov ${params.normMaxCov} --min_cov ${params.normMinCov} --left left-${sample_id}.fq --right right-${sample_id}.fq --pairs_together --PARALLEL_STATS --CPU ${task.cpus}

                echo -e "\\n-- DONE with Normalization --\\n"

                cat .command.out | grep "stats_file" -A 3 | tail -n 3 >${sample_id}_normStats.txt

                cp left.norm.fq left-"${sample_id}".norm.fq
                cp right.norm.fq right-"${sample_id}".norm.fq

                rm \$(readlink -e left.norm.fq) \$(readlink -e right.norm.fq ) left.norm.fq right.norm.fq
                """
            }
        }

        if (params.saveReads) {

            process save_norm_reads {

                label 'med_cpus'

                tag "${sample_id}"

                publishDir "${launchDir}/${params.outdir}/saveReads/normalization", mode: "copy", overwrite: true, pattern: "*.gz"

                conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "conda-forge::pigz=2.3.4=hed695b0_1" : null)
                if (params.oneContainer){ container "${params.TPcontainer}" } else {
                container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/pigz:2.3.4" : "quay.io/biocontainers/pigz:2.3.4")
                }

                input:
                    tuple sample_id, file(r1), file(r2) from save_norm_reads

                output:
                    tuple sample_id, file("*.gz") into save_norm_reads_out

                script:
                    """
                    cat $r1 >${sample_id}_norm.R1.fq
                    pigz --best --force -p ${task.cpus} -r ${sample_id}_norm.R1.fq
                    cat $r2 >${sample_id}_norm.R2.fq
                    pigz --best --force -p ${task.cpus} -r ${sample_id}_norm.R2.fq
                    """
            }
        }
    }

    if (params.skipFilter && params.skipNormalization) {

        process prepare_reads {

            label 'low_cpus'

            tag "${sample_id}"

            input:
                tuple sample_id, file(reads) from reads_rna_ass

            output:
                tuple sample_id, file("left-${sample_id}.fq"), file("right-${sample_id}.fq") into ( norm_reads_soap, norm_reads_velvet, norm_reads_trinity, norm_reads_spades, norm_reads_transabyss, reads_rna_quast )
                tuple sample_id, file("left-${sample_id}.fq"), file("right-${sample_id}.fq") into ( mapping_reads_trinity, mapping_reads_evi, mapping_symbiont )

            script:
            if (hasExtension(params.reads, 'gz')) {
                """
                echo ${sample_id}
                zcat ${reads[0]} >left-${sample_id}.fq &
                zcat ${reads[1]} >right-${sample_id}.fq
                """
            } else {
                """
                echo ${sample_id}
                cat ${reads[0]} >left-${sample_id}.fq &
                cat ${reads[1]} >right-${sample_id}.fq
                """
            }
        }

        process skip_filter {
            tag "${sample_id}"

            input:
                tuple sample_id, file(files) from skip_filter_ch

            output:
                tuple sample_id, file("filter_reads.txt") into fastp_csv

            script:
                """
                echo "Filter step was skipped" >filter_reads.txt
                """
        }

        process skip_normalization {
            tag "${sample_id}"

            input:
                tuple sample_id, file(files) from skip_norm_ch

            output:
                tuple sample_id, file("norm_reads.txt") into norm_report

            script:
                """
                echo "Normalization step was skipped" >norm_reads.txt
                """
        }

    } else if (!params.skipFilter && params.skipNormalization) {
        //reads_rna_ass
        //    .into{ norm_reads_soap; norm_reads_velvet; norm_reads_trinity; norm_reads_spades; norm_reads_transabyss; reads_rna_quast; mapping_reads_trinity; mapping_reads_evi; mapping_symbiont }

        process skip_normalization_only {
            tag "${sample_id}"

            input:
                tuple sample_id, file(reads) from reads_rna_ass

            output:
                tuple sample_id, file("norm_reads.txt") into norm_report
                tuple sample_id, file("left-${sample_id}.fq"), file("right-${sample_id}.fq") into norm_reads_soap, norm_reads_velvet, norm_reads_trinity, norm_reads_spades, norm_reads_transabyss, reads_rna_quast, mapping_reads_trinity, mapping_reads_evi, mapping_symbiont

            script:
                """
                echo "Normalization step was skipped" >norm_reads.txt
                echo ${sample_id}
                cat ${reads[0]} >left-${sample_id}.fq &
                cat ${reads[1]} >right-${sample_id}.fq
                """
        }

    } else if (params.skipFilter && !params.skipNormalization) {
        process skip_filter_only {
            tag "${sample_id}"

            input:
                tuple sample_id, file(files) from skip_filter_only_ch

            output:
                tuple sample_id, file("filter_reads.txt") into fastp_csv

            script:
                """
                echo "Filter step was skipped" >filter_reads.txt
                """
        }
    }

    process trinity_assembly {

        label 'med_mem'

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/assemblies", mode: "copy", overwrite: true

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::trinity=2.9.1=h8b12597_1" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/trinity:2.9.1--h8b12597_1" : "quay.io/biocontainers/trinity:2.9.1--h8b12597_1")
        }

        input:
            tuple sample_id, file(left), file(right) from norm_reads_trinity

        output:
            tuple sample_id, file("${sample_id}.Trinity.fa") into ( assemblies_ch_trinity, busco3_ch_trinity, busco4_ch_trinity, assemblies_ch_trinity_busco3, assemblies_ch_trinity_busco4, mapping_trinity )

        script:
            """
            mem=\$( echo ${task.memory} | cut -f 1 -d " " )

            Trinity --max_memory \${mem}G --seqType fq --left ${left} --right ${right} --CPU ${task.cpus} --no_normalize_reads --full_cleanup --output trinity_out_dir

            mv trinity_out_dir.Trinity.fasta ${sample_id}.Trinity.fa
            """
    }

    process soap_assembly {

        label 'med_mem'

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/assemblies", mode: "copy", overwrite: true

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::soapdenovo-trans=1.04=ha92aebf_2" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/soapdenovo-trans:1.04--ha92aebf_2" : "quay.io/biocontainers/soapdenovo-trans:1.04--ha92aebf_2")
        }

        input:
            val k from "${params.k}"
            tuple sample_id, file(left), file(right) from norm_reads_soap

        output:
            tuple sample_id, file("${sample_id}.SOAP.fa") into assemblies_ch_soap
            tuple sample_id, file("${sample_id}.SOAP.k*") into assemblies_ch_soap_busco3, assemblies_ch_soap_busco4

        script:
            """
            echo -e "\\n-- Generating SOAP config file --\\n"
            echo "maxReadLen="${params.maxReadLen} >>config.txt
            echo "[LIB]" >>config.txt
            echo "rd_len_cutof="${params.rd_len_cutof} >>config.txt
            #echo "avg_ins="${params.avg_ins} >>config.txt
            echo "reverse_seq="${params.reverse_seq} >>config.txt
            echo "asm_flags="${params.asm_flags} >>config.txt
            echo "map_len="${params.map_len} >>config.txt
            echo "q1="${left} >>config.txt
            echo "q2="${right} >>config.txt

            echo -e "\\n-- Starting SOAP assemblies --\\n"

            for x in `echo $k | tr "," " "`;do
                echo -e "\\n-- SOAP k\${x} --\\n"
                SOAPdenovo-Trans-127mer all -s config.txt -K \${x} -o output\${x} -p ${task.cpus}
                sed -i "s/>/>SOAP.k\${x}./g" output\${x}.scafSeq
            done

            echo -e "\\n-- Finished with the assemblies --\\n"

            cat output*.scafSeq >${sample_id}.SOAP.fa

            for x in `echo $k | tr "," " "`;do
                cp output\${x}.scafSeq ${sample_id}.SOAP.k\${x}.fa
            done

            rm -rf output*
            """
    }

    process velvet_oases_assembly {

        label 'med_mem'

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/assemblies", mode: "copy", overwrite: true

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge -c bioconda velvet=1.2.10 oases=0.2.09" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/mulled-v2-8ce10492777ba3fb1db6e6e13fa9b78ac116db2f:f54a9246f1216443f2e0f6de9ec5908ca882f710-0" : "quay.io/biocontainers/mulled-v2-8ce10492777ba3fb1db6e6e13fa9b78ac116db2f:f54a9246f1216443f2e0f6de9ec5908ca882f710-0")
        }

        input:
            val k from "${params.k}"
            tuple sample_id, file(left), file(right) from norm_reads_velvet

        output:
            tuple sample_id, file("${sample_id}.Velvet.fa") into assemblies_ch_velvet
            tuple sample_id, file("${sample_id}.Velvet.k*") into assemblies_ch_velvet_busco3, assemblies_ch_velvet_busco4

        script:
            """
    	    echo -e "\\n-- Starting with Velveth --\\n"
            for x in `echo $k | tr "," " "`;do
                echo -e "\\n-- k\${x} --\\n"
                velveth oases.\${x} \${x} -shortPaired -fastq -separate ${left} ${right}
            done

            echo -e "\\n-- Starting with Velvetg --\\n"
            for x in `echo $k | tr "," " "`;do
                echo -e "\\n-- vg \${x} --\\n"
                velvetg oases.\${x} -read_trkg yes
            done

            echo -e "\\n-- Starting with Oases --\\n"
            for x in `echo $k | tr "," " "`;do
                echo -e "\\n-- oases \${x} --\\n"
                oases oases.\${x}
            done

            echo -e "\\n-- Finished with Velvet/Oases assemblies --\\n"

            for x in `echo $k | tr "," " "`;do
                sed -i "s/>/>Velvet.k\${x}./g" oases.\${x}/contigs.fa
            done

            cat oases.*/contigs.fa >${sample_id}.Velvet.fa

            for x in `echo $k | tr "," " "`;do
                cp oases.\${x}/contigs.fa ${sample_id}.Velvet.k\${x}.fa
            done

            rm -rf oases.*
            """
    }

    process rna_spades_assembly {

        label 'med_mem'

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/assemblies", mode: "copy", overwrite: true

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::spades=3.14.0=h2d02072_0" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/spades:3.14.0--h2d02072_0" : "quay.io/biocontainers/spades:3.14.0--h2d02072_0")
        }

        input:
            val k from "${params.k}"
            tuple sample_id, file(left), file(right) from norm_reads_spades

        output:
            tuple sample_id, file("${sample_id}.SPADES.fa") into assemblies_ch_spades
            tuple sample_id, file("${sample_id}.SPADES.k*") into assemblies_ch_spades_busco3, assemblies_ch_spades_busco4

        script:
            """
            echo -e "\\n-- Starting rnaSPADES assemblies --\\n"

            mem=\$( echo ${task.memory} | cut -f 1 -d " " )

            for x in `echo $k | tr "," " "`;do
                echo -e "\\n-- rnaSPADES k\${x} --\\n"
                rnaspades.py -1 ${left} -2 ${right} -o ${sample_id}_spades_\${x} -t ${task.cpus} -k \${x} -m \${mem}
            done

            echo -e "\\n-- Finished with the assemblies --\\n"

            for x in `echo $k | tr "," " "`;do
                sed -i "s/>/>SPADES.k\${x}./g" ${sample_id}_spades_\${x}/transcripts.fasta
            done

            cat ${sample_id}_spades_*/transcripts.fasta >${sample_id}.SPADES.fa

            for x in `echo $k | tr "," " "`;do
                cp ${sample_id}_spades_\${x}/transcripts.fasta ${sample_id}.SPADES.k\${x}.fa
            done

            rm -rf ${sample_id}_spades_*
            """
    }

    process transabyss_assembly {

        label 'med_mem'

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/assemblies", mode: "copy", overwrite: true

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::transabyss=2.0.1=py_6" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/transabyss:2.0.1--py_6" : "quay.io/biocontainers/transabyss:2.0.1--py_6")
        }

        input:
            val k from "${params.k}"
            tuple sample_id, file(left), file(right) from norm_reads_transabyss

        output:
            tuple sample_id, file("${sample_id}.TransABySS.fa") into assemblies_ch_transabyss
            tuple sample_id, file("${sample_id}.TransABySS.k*") into assemblies_ch_transabyss_busco3, assemblies_ch_transabyss_busco4

        script:
            if ((workflow.containerEngine == 'singularity' || workflow.containerEngine == 'docker') && !params.oneContainer) {
            """
            mkdir -p /usr/local/lib/python3.8/site-packages/bin
            cp /usr/local/bin/skip_psl_self* /usr/local/lib/python3.8/site-packages/bin/
            echo -e "\\n-- Starting Trans-ABySS assemblies --\\n"

            for x in `echo $k | tr "," " "`;do
                echo -e "\\n-- Trans-ABySS k\${x} --\\n"
                transabyss -k \${x} --pe ${left} ${right} --outdir ${sample_id}_transabyss_\${x} --name k\${x}.transabyss.fa --threads ${task.cpus} -c 12 --length 200
            done

            echo -e "\\n-- Finished with the assemblies --\\n"

            for x in `echo $k | tr "," " "`;do
                sed -i "s/>/>TransABySS.k\${x}./g" ${sample_id}_transabyss_\${x}/k\${x}.transabyss.fa-final.fa
            done

            cat ${sample_id}_transabyss_*/k*.transabyss.fa-final.fa >${sample_id}.TransABySS.fa

            for x in `echo $k | tr "," " "`;do
                cp ${sample_id}_transabyss_\${x}/k\${x}.transabyss.fa-final.fa ${sample_id}.TransABySS.k\${x}.fa
            done

            rm -rf ${sample_id}_transabyss_*
            """
            } else {
            """
            echo -e "\\n-- Starting Trans-ABySS assemblies --\\n"

            for x in `echo $k | tr "," " "`;do
                echo -e "\\n-- Trans-ABySS k\${x} --\\n"
                transabyss -k \${x} --pe ${left} ${right} --outdir ${sample_id}_transabyss_\${x} --name k\${x}.transabyss.fa --threads ${task.cpus} -c 12 --length 200
            done

            echo -e "\\n-- Finished with the assemblies --\\n"

            for x in `echo $k | tr "," " "`;do
                sed -i "s/>/>TransABySS.k\${x}./g" ${sample_id}_transabyss_\${x}/k\${x}.transabyss.fa-final.fa
            done

            cat ${sample_id}_transabyss_*/k*.transabyss.fa-final.fa >${sample_id}.TransABySS.fa

            for x in `echo $k | tr "," " "`;do
                cp ${sample_id}_transabyss_\${x}/k\${x}.transabyss.fa-final.fa ${sample_id}.TransABySS.k\${x}.fa
            done

            rm -rf ${sample_id}_transabyss_*
            """
            }
    }

    all_assemblies = Channel.create()
    assemblies_ch_trinity.mix( assemblies_ch_transabyss, assemblies_ch_spades, assemblies_ch_velvet, assemblies_ch_soap ).groupTuple(by:0,size:5).into(all_assemblies)

    process evigene {

        label 'med_mem'

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/evigene", mode: "copy", overwrite: true

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::cd-hit=4.8.1 bioconda::exonerate=2.4 bioconda::blast=2.2.31" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/mulled-v2-962eae98c9ff8d5b31e1df7e41a355a99e1152c4:5aac6d1d2253d47aee81f01cc070a17664c86f07-0" : "quay.io/biocontainers/mulled-v2-962eae98c9ff8d5b31e1df7e41a355a99e1152c4:5aac6d1d2253d47aee81f01cc070a17664c86f07-0")
        }

        input:
            tuple sample_id, file(assemblies) from all_assemblies

        output:
            tuple sample_id, file("${sample_id}.combined.okay.fa") into ( evigene_ch_busco3, evigene_ch_busco4, evigene_ch_transdecoder, evigene_ch_transdecoderB, evigene_ch_rnammer, evigene_ch_trinotate, evi_dist, evi_filt, evigene_ch_rna_quast, mapping_evi )
            tuple sample_id, file("${sample_id}.combined.fa"), file("${sample_id}.combined.okay.fa") into evigene_summary

        script:
            def mem_MB=(task.memory.toMega())
            if ((workflow.containerEngine == 'singularity' || workflow.containerEngine == 'docker') && params.oneContainer) {
                """
                echo -e "\\n-- Starting EviGene --\\n"

                cat ${assemblies} >${sample_id}.combined.fa

                tr2aacds.pl -tidy -NCPU ${task.cpus} -MAXMEM ${mem_MB} -log -cdna ${sample_id}.combined.fa

                echo -e "\\n-- DONE with EviGene --\\n"

                cp okayset/*combined.okay*.fa ${sample_id}.combined.okay.fa
                cp okayset/*combined.okay*.cds ${sample_id}.combined.okay.cds

                if [ -d tmpfiles/ ];then
                    rm -rf tmpfiles/
                fi
                """
            } else {
                """
                echo -e "\\n-- Starting EviGene --\\n"

                cat ${assemblies} >${sample_id}.combined.fa

                $evi/scripts/prot/tr2aacds.pl -tidy -NCPU ${task.cpus} -MAXMEM ${mem_MB} -log -cdna ${sample_id}.combined.fa

                echo -e "\\n-- DONE with EviGene --\\n"

                cp okayset/*combined.okay*.fa ${sample_id}.combined.okay.fa
                cp okayset/*combined.okay*.cds ${sample_id}.combined.okay.cds

                if [ -d tmpfiles/ ];then
                    rm -rf tmpfiles/
                fi
            	"""
            }
    }

    // check groupTuple
    rna_quast = Channel.create()
    reads_rna_quast.join( evigene_ch_rna_quast ).into( rna_quast )

    process rna_quast {

        label 'med_mem'

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/rnaQuast", mode: "copy", overwrite: true

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::rnaquast=2.0.1=0" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/rnaquast:2.0.1--0" : "quay.io/biocontainers/rnaquast:2.0.1--0")
        }

        input:
            tuple sample_id, file(r1), file(r2), file(assembly) from rna_quast

        output:
            tuple sample_id, file("${sample_id}.rna_quast") into rna_quast_sum
            tuple sample_id, file("${sample_id}_rnaQUAST.csv") into rna_quast_report

        script:
            """
            rnaQUAST.py --transcripts ${assembly} -1 ${r1} -2 ${r2} -o ${sample_id}.rna_quast -t ${task.cpus} --blat
            echo "Category,Value" >${sample_id}_rnaQUAST.csv
            cat ${sample_id}.rna_quast/*_output/basic_metrics.txt | grep -v "METRICS" |  sed 's/\\(\\ \\)* \\([0-9]\\)/,\\2/g' | sed 's/>,/>/g' | grep [0-9] >>${sample_id}_rnaQUAST.csv
            cat Sponge_sample.rna_quast/*_output/sensitivity.txt | grep "Genes" | sed 's/\\(\\ \\)* \\([0-9]\\)/,\\2/g' >>${sample_id}_rnaQUAST.csv
            """
    }

    mapping_evi_in=Channel.create()
    mapping_evi.mix( mapping_reads_evi ).groupTuple(by:0,size:2).into(mapping_evi_in)

    process mapping_evigene {

        label 'big_cpus'

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/mapping", mode: "copy", overwrite: true

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::bowtie2=2.3.5.1=py36he513fc3_0" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/bowtie2:2.3.5.1--py36he513fc3_0" : "quay.io/biocontainers/bowtie2:2.3.5.1--py36he513fc3_0")
        }

        input:
            tuple sample_id, file(files), file(files2) from mapping_evi_in

        output:
            tuple sample_id, file("log*") into mapping_evi_results
            tuple sample_id, file("*") into mapping_evi_results_bam

        script:
        if (params.saveBam) {
            """
            a=\$( echo $files $files2 )
            ass=\$( echo \$a | tr " " "\\n" | grep ".combined.okay.fa" )
            r1=\$( echo \$a | tr " " "\\n" | grep "left-" )
            r2=\$( echo \$a | tr " " "\\n" | grep "right-" )
            bowtie2-build \${ass} \${ass} --threads ${task.cpus}
            bowtie2 -x \${ass} -1 \${r1} -2 \${r2} -p ${task.cpus} 2>log_\${ass}.txt | samtools view -@ ${task.cpus} -bS - >\${ass}.bam
            rm *.bt2
            """
        } else {
            """
            a=\$( echo $files $files2 )
            ass=\$( echo \$a | tr " " "\\n" | grep ".combined.okay.fa" )
            r1=\$( echo \$a | tr " " "\\n" | grep "left-" )
            r2=\$( echo \$a | tr " " "\\n" | grep "right-" )
            bowtie2-build \${ass} \${ass} --threads ${task.cpus}
            bowtie2 -x \${ass} -1 \${r1} -2 \${r2} -p ${task.cpus} 2>log_\${ass}.txt | samtools view -@ ${task.cpus} -bS - >\${ass}.bam
            rm *.bam *.bt2
            """
        }
    }

    mapping_trinity_in=Channel.create()
    mapping_trinity.mix( mapping_reads_trinity ).groupTuple(by:0,size:2).into(mapping_trinity_in)

    process mapping_trinity {

        label 'big_cpus'

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/mapping", mode: "copy", overwrite: true

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::bowtie2=2.3.5.1=py36he513fc3_0" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/bowtie2:2.3.5.1--py36he513fc3_0" : "quay.io/biocontainers/bowtie2:2.3.5.1--py36he513fc3_0")
        }

        input:
            tuple sample_id, file(files), file(files2) from mapping_trinity_in

        output:
            tuple sample_id, file("log*") into mapping_trinity_results

        script:
            """
            a=\$( echo $files $files2 )
            ass=\$( echo \$a | tr " " "\\n" | grep ".Trinity.fa" )
            r1=\$( echo \$a | tr " " "\\n" | grep "left-" )
            r2=\$( echo \$a | tr " " "\\n" | grep "right-" )
            bowtie2-build \${ass} \${ass} --threads ${task.cpus}
            bowtie2 -x \${ass} -1 \${r1} -2 \${r2} -p ${task.cpus} 2>log_\${ass}.txt | samtools view -@ ${task.cpus} -bS - >\${ass}.bam
            rm *.bam
            """

    }

    if (params.filterSpecies) {

        if (params.host && params.symbiont) {

            File file1 = new File("${params.host}")
            String hostPath = file1.getCanonicalPath()

            File file2 = new File("${params.symbiont}")
            String symbiontPath = file2.getCanonicalPath()

            host_sequences = file(hostPath)
            symbiont_sequences = file(symbiontPath)
            transcriptome_sequences1 = Channel.create()
            transcriptome_sequences2 = Channel.create()
            evi_filt.into( transcriptome_sequences1, transcriptome_sequences2 )

            process update_host_ids {

                label 'med_cpus'

            	input:
            	   file(host) from host_sequences

            	output:
            	   file("updated_host.fasta") into updated_host

                script:
                if (hasExtension(host_sequences, 'gz')) {
                    """
                	zcat ${host} | awk '{if (\$1 ~ /^ *>/ ) {print substr(\$1,0,1)"species1_" substr(\$1,2,length(\$1))} else {print \$1} }' >updated_host.fasta
                	"""
                } else {
                	"""
                	awk '{if (\$1 ~ /^ *>/ ) {print substr(\$1,0,1)"species1_" substr(\$1,2,length(\$1))} else {print \$1} }' < ${host} >updated_host.fasta
                	"""
                }
            }

            process update_symbio_ids {

                label 'med_cpus'

                input:
                    file(symbiont) from symbiont_sequences

                output:
                    file("updated_symbio.fasta") into updated_symbiont

                script:
                if (hasExtension(symbiont_sequences, 'gz')) {
                    """
                    zcat ${symbiont} | awk '{if (\$1 ~ /^ *>/ ) {print substr(\$1,0,1)"species2_" substr(\$1,2,length(\$1))} else {print \$1} }' >updated_symbio.fasta
                    """
                } else {
                    """
                    awk '{if (\$1 ~ /^ *>/ ) {print substr(\$1,0,1)"species2_" substr(\$1,2,length(\$1))} else {print \$1} }' < ${symbiont} >updated_symbio.fasta
                    """
                }
            }

            process concatenate_sets {

                label 'exlow_cpus'

            	input:
                	file(host) from updated_host
                	file(symbiont) from updated_symbiont

            	output:
            	   file("concatenated_file.fasta") into file_for_database

                script:
                	"""
                	cat ${host} ${symbiont} >concatenated_file.fasta
                	"""
            }

            process create_diamond_db {

                label 'low_mem'

                conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::diamond=0.9.30=h56fc30b_0" : null)
                if (params.oneContainer){ container "${params.TPcontainer}" } else {
                container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/diamond:0.9.30--h56fc30b_0" : "quay.io/biocontainers/diamond:0.9.30--h56fc30b_0")
                }

            	input:
            	   file(seqs) from file_for_database

            	output:
            	   file("*.dmnd") into diamond_db

                script:
                	"""
                	diamond makedb --in ${seqs} -d diamond_database -p ${task.cpus}
                	"""
            }

            process diamond_run {

                label 'med_cpus'

                tag "${sample_id}"

                conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::diamond=0.9.30=h56fc30b_0" : null)
                if (params.oneContainer){ container "${params.TPcontainer}" } else {
                container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/diamond:0.9.30--h56fc30b_0" : "quay.io/biocontainers/diamond:0.9.30--h56fc30b_0")
                }

            	input:
                	tuple sample_id, file(transcriptome) from transcriptome_sequences1
                	file(database) from diamond_db

            	output:
            	   file("*.tabular") into diamond_output

                script:
                	"""
                	diamond blastx -d ${database} -q ${transcriptome} -p ${task.cpus} -f 6 -o diamond_output.tabular
                	"""
            }

            process psytrans_run {

                label 'med_cpus'

                tag "${sample_id}"

            	publishDir "${launchDir}/${params.outdir}/psytrans_output/", mode: "copy", overwrite: true

                conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge -c bioconda biopython=1.78 pandas=1.1.2 numpy=1.18.1" : null)
                if (params.oneContainer){ container "${params.TPcontainer}" } else {
                container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/mulled-v2-1e9d4f78feac0eb2c8d8246367973b3f6358defc:41ffac721ff9b03ca1121742e969d0e7d78e589f-0" : "quay.io/biocontainers/mulled-v2-1e9d4f78feac0eb2c8d8246367973b3f6358defc:41ffac721ff9b03ca1121742e969d0e7d78e589f-0")
                }

            	input:
                	tuple sample_id, file(transcriptome) from transcriptome_sequences2
                	file(host) from updated_host
                	file(symbiont) from updated_symbiont
                	file(database) from diamond_output

            	output:
            	   file("*.fasta") into psytrans

                script:
                	"""
                	mkdir temp
                	psytrans.py ${transcriptome} -A ${host} -B ${symbiont} -b ${database} -t temp -n ${params.psyval}
                    mv species1_${transcriptome} ${sample_id}_only.fasta
                    mv species2_${transcriptome} ${sample_id}_removed.fasta
                	"""
            }
        } else {
            println("\n\t\033[0;31mNeed to provide a host and symbiont sequence.\n\tFor more info use `nextflow run TransPi.nf --help`\n\033[0m")
            exit 0
        }

        // mapping
        if (params.mapping) {

            if (params.symbiont) {

                File file2 = new File("${params.symbiont}")
                String symbiontPath = file2.getCanonicalPath()
                symbiont_sequences = file(symbiontPath)

                process mapping_reads {

                    label 'med_cpus'

                    tag "${sample_id}"

                	publishDir "${launchDir}/${params.outdir}/mapping_output/", mode: "copy", overwrite: true

                    conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::bowtie2=2.3.5.1=py36he513fc3_0" : null)
                    if (params.oneContainer){ container "${params.TPcontainer}" } else {
                    container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/bowtie2:2.3.5.1--py36he513fc3_0" : "quay.io/biocontainers/bowtie2:2.3.5.1--py36he513fc3_0")
                    }

                	input:
                        file(symbiont) from symbiont_sequences
                        tuple sample_id, file(r1), file(r2) from mapping_symbiont

                	output:
                	   tuple sample_id, file("${sample_id}.symbiont.sam") into clean_reads_ch

                    script:
                    if (hasExtension(symbiont_sequences, 'gz')) {
                        """
                        zcat ${symbiont} >symbiont.seqs
                    	bowtie2-build symbiont.seqs symbiont.db --threads ${task.cpus}
                        bowtie2 -x symbiont.db -1 ${r1} -2 ${r2} --no-unal -p ${task.cpus} -S ${sample_id}.symbiont.sam
                    	"""
                    } else {
                        """
                    	bowtie2-build ${symbiont} symbiont.db --threads ${task.cpus}
                        bowtie2 -x symbiont.db -1 ${r1} -2 ${r2} --no-unal -p ${task.cpus} -S ${sample_id}.symbiont.sam
                    	"""
                    }

                }

                process clean_reads {

                    label 'med_cpus'

                    tag "${sample_id}"

                    publishDir "${launchDir}/${params.outdir}/mapping_output/", mode: "copy", overwrite: true

                    conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge -c bioconda biopython=1.78 pandas=1.1.2 numpy=1.18.1" : null)
                    if (params.oneContainer){ container "${params.TPcontainer}" } else {
                    container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/mulled-v2-1e9d4f78feac0eb2c8d8246367973b3f6358defc:41ffac721ff9b03ca1121742e969d0e7d78e589f-0" : "quay.io/biocontainers/mulled-v2-1e9d4f78feac0eb2c8d8246367973b3f6358defc:41ffac721ff9b03ca1121742e969d0e7d78e589f-0")
                    }

                    input:
                        tuple sample_id, file("${sample_id}.symbiont.sam") from clean_reads_ch

                    output:
                       file("*.fasta") into ch_FINISH

                    script:
                        """
                        clean_seqs.py ${sample_id}.symbiont.sam
                        """

                }
            }
        }
    //
    }

    def slist="${params.k}".split(',').collect{it as int}
    def m=slist.size()
    def n=4*m+1

    if (params.allBuscos) {

        // estimate the number of files, either 5 or n
        busco3_all = Channel.create()
        assemblies_ch_soap_busco3.mix( assemblies_ch_velvet_busco3, assemblies_ch_spades_busco3, assemblies_ch_transabyss_busco3, assemblies_ch_trinity_busco3 ).groupTuple(by:0,size:5).into(busco3_all)

        process busco3_all {

            label 'med_cpus'

            tag "${sample_id}"

            publishDir "${launchDir}/${params.outdir}/busco3_all", mode: "copy", overwrite: true, pattern: "*.{tsv,txt,bus3}"

            conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::busco=3.0.2=py_13" : null)
            if (params.oneContainer){ container "${params.TPcontainer}" } else {
            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/busco:3.0.2--py_13" : "quay.io/biocontainers/busco:3.0.2--py_13")
            }

            input:
                tuple sample_id, file(files) from busco3_all

            output:
                tuple sample_id, file("*.bus3") into busco3_all_ch
                tuple sample_id, file("*.Trinity.bus3.txt") into ( busco3_ch_trinity_sum, busco3_comp_2 )
                tuple sample_id, file("*.txt"), file("*.tsv") into busco3_all_sum_ch
                tuple sample_id, file("${sample_id}_all_busco3.tsv"), file("${sample_id}_all_assemblers.fa") into busco3_all_tsv

            script:
                """
                cat input.1 | sed 's/, /\\n/g' | tr -d "[" | tr "]" "\\n" >list.txt
                cat input.2 | sed 's/, /\\n/g' | tr -d "[" | tr "]" "\\n" >>list.txt
                cat input.3 | sed 's/, /\\n/g' | tr -d "[" | tr "]" "\\n" >>list.txt
                cat input.4 | sed 's/, /\\n/g' | tr -d "[" | tr "]" "\\n" >>list.txt

                for x in `cat list.txt`;do
                    ln -s \$x \$( basename \$x )
                done

                find . -maxdepth 1 -type l -ls | grep "Trinity" | awk -F "-> " '{print \$2}' >>list.txt

                for x in `cat list.txt`;do

                    name=\$( basename \$x .fa )

                    echo -e "\\n-- Starting BUSCO --\\n"

                    run_BUSCO.py -i \${name}.fa -o \${name}.bus3 -l ${params.busco3db} -m tran -c ${task.cpus}

                    cp run_\${name}.bus3/short* .
                    cp run_\${name}.bus3/full_table_* .

                    echo -e "\\n-- DONE with BUSCO --\\n"

                done

                echo "Busco_id,Status,Sequence,Score,Length" >.header.txt
                cat full_table_*.tsv | grep -v "#" | tr "\t" "," >.busco_names.txt
                cat .header.txt .busco_names.txt >${sample_id}_all_busco3.tsv
                rm .header.txt .busco_names.txt

                cat *.fa >${sample_id}_all_assemblers.fa
                """
        }

        busco4_all = Channel.create()
        assemblies_ch_soap_busco4.mix( assemblies_ch_velvet_busco4, assemblies_ch_spades_busco4, assemblies_ch_transabyss_busco4, assemblies_ch_trinity_busco4 ).groupTuple(by:0,size:5).into(busco4_all)

        process busco4_all {

            label 'med_cpus'

            tag "${sample_id}"

            publishDir "${launchDir}/${params.outdir}/busco4_all", mode: "copy", overwrite: true, pattern: "*.{tsv,txt,bus4}"

            // change container in oneContainer option
            conda (params.condaActivate && params.myConda ? params.cenv : params.condaActivate ? "-c conda-forge bioconda::busco=4.1.4=py_0" : null)
            if (params.oneContainer){ container "${params.v4container}" } else {
            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/busco:4.0.5--pyr36_0" : "ezlabgva/busco:v4.0.5_cv1")
            }

            input:
                tuple sample_id, file(files) from busco4_all

            output:
                tuple sample_id, file("*.bus4") into busco4_all_ch
                tuple sample_id, file("*.Trinity.bus4.txt") into ( busco4_ch_trinity_sum, busco4_comp_2 )
                tuple sample_id, file("*.txt"), file("*.tsv") into busco4_all_sum_ch
                tuple sample_id, file("${sample_id}_all_busco4.tsv"), file("${sample_id}_all_assemblers.fa") into busco4_all_tsv

            script:
                """
                cat input.1 | sed 's/, /\\n/g' | tr -d "[" | tr "]" "\\n" >list.txt
                cat input.2 | sed 's/, /\\n/g' | tr -d "[" | tr "]" "\\n" >>list.txt
                cat input.3 | sed 's/, /\\n/g' | tr -d "[" | tr "]" "\\n" >>list.txt
                cat input.4 | sed 's/, /\\n/g' | tr -d "[" | tr "]" "\\n" >>list.txt

                for x in `cat list.txt`;do
                    ln -s \$x \$( basename \$x )
                done

                find . -maxdepth 1 -type l -ls | grep "Trinity" | awk -F "-> " '{print \$2}' >>list.txt

                for x in `cat list.txt`;do

                    name=\$( basename \$x .fa )

                    echo -e "\\n-- Starting BUSCO --\\n"

                    busco -i \${name}.fa -o \${name}.bus4 -l ${params.busco4db} -m tran -c ${task.cpus} --offline

                    cp \${name}.bus4/short_summary.* .
                    cp \${name}.bus4/run_*/full_table.tsv full_table_\${name}.tsv

                    echo -e "\\n-- DONE with BUSCO --\\n"

                done

                echo "Busco_id,Status,Sequence,Score,Length" >.header.txt
                cat full_table_*.tsv | grep -v "#" | tr "\t" "," >.busco_names.txt
                cat .header.txt .busco_names.txt >${sample_id}_all_busco4.tsv
                rm .header.txt .busco_names.txt

                cat *.fa >${sample_id}_all_assemblers.fa
                """
        }

    } else {

        process busco3_tri {

            label 'med_cpus'

            tag "${sample_id}"

            publishDir "${launchDir}/${params.outdir}/busco3", mode: "copy", overwrite: true

            conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::busco=3.0.2=py_13" : null)
            if (params.oneContainer){ container "${params.TPcontainer}" } else {
            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/busco:3.0.2--py_13" : "quay.io/biocontainers/busco:3.0.2--py_13")
            }

            input:
                tuple sample_id, file("${sample_id}.Trinity.fa") from busco3_ch_trinity

            output:
                tuple sample_id, file("*${sample_id}.Trinity.bus3.txt") into ( busco3_ch_trinity_sum, busco3_comp_2 )
                file("run_${sample_id}.Trinity.bus3")
                tuple sample_id, file("*tsv") into busco3_trinity_rescue

            script:
                """
                echo -e "\\n-- Starting BUSCO --\\n"

                run_BUSCO.py -i ${sample_id}.Trinity.fa -o ${sample_id}.Trinity.bus3 -l ${params.busco3db} -m tran -c ${task.cpus}

                echo -e "\\n-- DONE with BUSCO --\\n"

                cp run_${sample_id}.Trinity.bus3/short_summary_${sample_id}.Trinity.bus3.txt .
                cp run_${sample_id}.Trinity.bus3/full_table_${sample_id}.Trinity.bus3.tsv .
                """
        }

        process busco4_tri {

            label 'med_cpus'

            tag "${sample_id}"

            publishDir "${launchDir}/${params.outdir}/busco4", mode: "copy", overwrite: true

            // change container in oneContainer option
            conda (params.condaActivate && params.myConda ? params.cenv : params.condaActivate ? "-c conda-forge bioconda::busco=4.1.4=py_0" : null)
            if (params.oneContainer){ container "${params.v4container}" } else {
            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/busco:4.0.5--pyr36_0" : "ezlabgva/busco:v4.0.5_cv1")
            }

            input:
                tuple sample_id, file("${sample_id}.Trinity.fa") from busco4_ch_trinity

            output:
                tuple sample_id, file("*${sample_id}.Trinity.bus4.txt") into ( busco4_ch_trinity_sum, busco4_comp_2 )
                file("${sample_id}.Trinity.bus4")
                tuple sample_id, file("*tsv") into busco4_trinity_rescue

            script:
                """
                echo -e "\\n-- Starting BUSCO --\\n"

                busco -i ${sample_id}.Trinity.fa -o ${sample_id}.Trinity.bus4 -l ${params.busco4db} -m tran -c ${task.cpus} --offline

                echo -e "\\n-- DONE with BUSCO --\\n"

                cp ${sample_id}.Trinity.bus4/short_summary.*.${sample_id}.Trinity.bus4.txt .
                cp ${sample_id}.Trinity.bus4/run_*/full_table.tsv full_table_${sample_id}.Trinity.bus4.tsv
                """
        }

    }

    process busco3 {

        label 'med_cpus'

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/busco3", mode: "copy", overwrite: true

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::busco=3.0.2=py_13" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/busco:3.0.2--py_13" : "quay.io/biocontainers/busco:3.0.2--py_13")
        }

        input:
            tuple sample_id, file("${sample_id}.combined.okay.fa") from evigene_ch_busco3

        output:
            tuple sample_id, file("run_${sample_id}.TransPi.bus3") into busco3_ch
            tuple sample_id, file("*${sample_id}.TransPi.bus3.txt") into ( busco3_summary, busco3_comp_1 )
            tuple sample_id, file("*tsv") into busco3_transpi_tsv

        script:
            """
            echo -e "\\n-- Starting BUSCO --\\n"

            run_BUSCO.py -i ${sample_id}.combined.okay.fa -o ${sample_id}.TransPi.bus3 -l ${params.busco3db} -m tran -c ${task.cpus}

            echo -e "\\n-- DONE with BUSCO --\\n"

            cp run_${sample_id}.TransPi.bus3/short_summary_${sample_id}.TransPi.bus3.txt .
            cp run_${sample_id}.TransPi.bus3/full_table_* .
            """
    }

    process busco4 {

        label 'med_cpus'

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/busco4", mode: "copy", overwrite: true

        // change container in oneContainer option
        conda (params.condaActivate && params.myConda ? params.cenv : params.condaActivate ? "-c conda-forge bioconda::busco=4.1.4=py_0" : null)
        if (params.oneContainer){ container "${params.v4container}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/busco:4.0.5--pyr36_0" : "ezlabgva/busco:v4.0.5_cv1")
        }

        input:
            tuple sample_id, file("${sample_id}.combined.okay.fa") from evigene_ch_busco4

        output:
            tuple sample_id, file("${sample_id}.TransPi.bus4") into busco4_ch
            tuple sample_id, file("*${sample_id}.TransPi.bus4.txt") into ( busco4_summary, busco4_comp_1 )
            tuple sample_id, file("*tsv") into busco4_transpi_tsv

        script:
            """
            echo -e "\\n-- Starting BUSCO --\\n"

            busco -i ${sample_id}.combined.okay.fa -o ${sample_id}.TransPi.bus4 -l ${params.busco4db} -m tran -c ${task.cpus} --offline

            echo -e "\\n-- DONE with BUSCO --\\n"

            cp ${sample_id}.TransPi.bus4/short_summary.*.${sample_id}.TransPi.bus4.txt .
            cp ${sample_id}.TransPi.bus4/run_*/full_table.tsv full_table_${sample_id}.TransPi.bus4.tsv
            """
    }

    // default perhaps
    if (params.buscoDist && params.allBuscos) {

        busco3_dist_ch=Channel.create()
        busco3_transpi_tsv.join( busco3_all_tsv ).into( busco3_dist_ch )

        process busco3_dist {

            label 'low_cpus'

            tag "${sample_id}"

            publishDir "${launchDir}/${params.outdir}/busco3_dist", mode: "copy", overwrite: true, pattern: "*.{tsv,fasta}"

            conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge -c bioconda biopython=1.78 pandas=1.1.2 numpy=1.18.1" : null)
            if (params.oneContainer){ container "${params.TPcontainer}" } else {
            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/mulled-v2-1e9d4f78feac0eb2c8d8246367973b3f6358defc:41ffac721ff9b03ca1121742e969d0e7d78e589f-0" : "quay.io/biocontainers/mulled-v2-1e9d4f78feac0eb2c8d8246367973b3f6358defc:41ffac721ff9b03ca1121742e969d0e7d78e589f-0")
            }

            input:
                tuple sample_id, file(transpi_tsv), file(all_busco), file(assembly) from busco3_dist_ch

            output:
                tuple sample_id, file("*.fasta"), file("*_table.tsv") into busco3_dist_sum
                tuple sample_id, file("*_table.tsv") into busco3_heatmap

            script:
                """
                cat $transpi_tsv | grep -v "#" | tr "\\t" "," >>$all_busco
                SOS_busco.py -input_file_busco $all_busco -input_file_fasta $assembly -min ${params.minPerc} -kmers ${params.k}
                mv Complete_comparison_table ${sample_id}_complete_BUSCO3_table.tsv
                mv TransPi_comparison_table ${sample_id}_TransPi_missing_BUSCO3_table.tsv
                if [ -e sequences_to_add.fasta ];then
                    mv sequences_to_add.fasta ${sample_id}_rescued_BUSCO4.fasta
                else
                    touch ${sample_id}_rescued_BUSCO4.fasta
                fi
                """

        }

        busco4_dist_ch=Channel.create()
        busco4_transpi_tsv.join( busco4_all_tsv ).into( busco4_dist_ch )

        process busco4_dist {

            label 'low_cpus'

            tag "${sample_id}"

            publishDir "${launchDir}/${params.outdir}/busco4_dist", mode: "copy", overwrite: true, pattern: "*.{tsv,fasta}"

            conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge -c bioconda biopython=1.78 pandas=1.1.2 numpy=1.18.1" : null)
            if (params.oneContainer){ container "${params.TPcontainer}" } else {
            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/mulled-v2-1e9d4f78feac0eb2c8d8246367973b3f6358defc:41ffac721ff9b03ca1121742e969d0e7d78e589f-0" : "quay.io/biocontainers/mulled-v2-1e9d4f78feac0eb2c8d8246367973b3f6358defc:41ffac721ff9b03ca1121742e969d0e7d78e589f-0")
            }

            input:
                tuple sample_id, file(transpi_tsv), file(all_busco), file(assembly) from busco4_dist_ch

            output:
                tuple sample_id, file("*.fasta"), file("*_table.tsv") into busco4_dist_sum
                tuple sample_id, file("*_table.tsv") into busco4_heatmap

            script:
                """
                cat $transpi_tsv | grep -v "#" | tr "\\t" "," >>$all_busco
                SOS_busco.py -input_file_busco $all_busco -input_file_fasta $assembly -min ${params.minPerc} -kmers ${params.k}
                mv Complete_comparison_table ${sample_id}_complete_BUSCO4_table.tsv
                mv TransPi_comparison_table ${sample_id}_TransPi_missing_BUSCO4_table.tsv
                if [ -e sequences_to_add.fasta ];then
                    mv sequences_to_add.fasta ${sample_id}_rescued_BUSCO4.fasta
                else
                    touch ${sample_id}_rescued_BUSCO4.fasta
                fi
                """

        }

    } else {
        process skip_busco_dist {
            tag "${sample_id}"

            input:
                tuple sample_id, file(files) from skip_busco_dist

            output:
                tuple sample_id, file("busco3_dist.txt") into busco3_heatmap
                tuple sample_id, file("busco4_dist.txt") into busco4_heatmap

            script:
                """
                echo "BUSCO3 distribution analysis was skipped" >busco3_dist.txt
                echo "BUSCO4 distribution analysis was skipped" >busco4_dist.txt
                """
        }
    }

    if (params.shortTransdecoder) {

        process transdecoder_short {

            label 'med_cpus'

            tag "${sample_id}"

            publishDir "${launchDir}/${params.outdir}/transdecoder", mode: "copy", overwrite: true

            conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::transdecoder=5.5.0=pl526_2" : null)
            if (params.oneContainer){ container "${params.TPcontainer}" } else {
            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/transdecoder:5.5.0--pl526_2" : "quay.io/biocontainers/transdecoder:5.5.0--pl526_2")
            }

            input:
                tuple sample_id, file(assembly) from evigene_ch_transdecoder

            output:
                tuple sample_id, file("${sample_id}.combined.okay.fa.transdecoder.pep") into ( transdecoder_ch_hmmer, transdecoder_ch_signalp, transdecoder_ch_tmhmm, transdecoder_ch_trinotate )
                tuple sample_id, file("${assembly}"), file("${sample_id}.combined.okay.fa.transdecoder.pep") into ( transdecoder_ch_diamond, transdecoder_ch_diamond_custom )
                tuple sample_id, file("${sample_id}_transdecoder.stats") into transdecoder_summary
                tuple sample_id, file("${sample_id}_transdecoder.csv") into transdecoder_csv
                tuple sample_id, file("${sample_id}*.transdecoder.{cds,gff,bed}") into transdecoder_files

            script:
                """
                echo -e "\\n-- TransDecoder.LongOrfs... --\\n"

                TransDecoder.LongOrfs -t ${assembly} --output_dir ${sample_id}.transdecoder_dir

                echo -e "\\n-- Done with TransDecoder.LongOrfs --\\n"

                echo -e "\\n-- TransDecoder.Predict... --\\n"

                TransDecoder.Predict -t ${assembly} --output_dir ${sample_id}.transdecoder_dir

                echo -e "\\n-- Done with TransDecoder.Predict --\\n"

                echo -e "\\n-- Calculating statistics... --\\n"
                #Calculate statistics of Transdecoder
                echo "- Transdecoder (short,no homolgy) stats for ${sample_id}" >${sample_id}_transdecoder.stats
                orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep -c ">" )
                echo -e "Total number of ORFs: \$orfnum \\n" >>${sample_id}_transdecoder.stats
                orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep -c "ORF type:complete" )
                echo -e "\\t ORFs type=complete: \$orfnum \\n" >>${sample_id}_transdecoder.stats
                orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep -c "ORF type:5prime_partial" )
                echo -e "\\t ORFs type=5prime_partial: \$orfnum \\n" >>${sample_id}_transdecoder.stats
                orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep -c "ORF type:3prime_partial" )
                echo -e "\\t ORFs type=3prime_partial: \$orfnum \\n" >>${sample_id}_transdecoder.stats
                orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep -c "ORF type:internal" )
                echo -e "\\t ORFs type=internal: \$orfnum \\n">>${sample_id}_transdecoder.stats
                # csv for report
                echo "Sample,Total_orf,orf_complete,orf_5prime_partial,orf_3prime_partial,orf_internal" >${sample_id}_transdecoder.csv
                total=\$( cat ${sample_id}*.transdecoder.pep  | grep -c ">" )
                complete=\$( cat ${sample_id}*.transdecoder.pep  | grep -c "ORF type:complete" )
                n5prime=\$( cat ${sample_id}*.transdecoder.pep  | grep -c "ORF type:5prime_partial" )
                n3prime=\$( cat ${sample_id}*.transdecoder.pep  | grep -c "ORF type:3prime_partial" )
                internal=\$( cat ${sample_id}*.transdecoder.pep  | grep -c "ORF type:internal" )
                echo "${sample_id},\${total},\${complete},\${n5prime},\${n3prime},\${internal}" >>${sample_id}_transdecoder.csv
                echo -e "\\n-- Done with statistics --\\n"

                cp ${assembly} ${assembly}.tmp
                rm ${assembly}
                mv ${assembly}.tmp ${assembly}

                echo -e "\\n-- DONE with TransDecoder --\\n"
                """
            }

    } else {

        process transdecoder_longorf {

            label 'med_cpus'

            tag "${sample_id}"

            publishDir "${launchDir}/${params.outdir}/transdecoder", mode: "copy", overwrite: true

            conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::transdecoder=5.5.0=pl526_2" : null)
            if (params.oneContainer){ container "${params.TPcontainer}" } else {
            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/transdecoder:5.5.0--pl526_2" : "quay.io/biocontainers/transdecoder:5.5.0--pl526_2")
            }

            input:
                tuple sample_id, file(assembly) from evigene_ch_transdecoder

            output:
                tuple sample_id, file("${sample_id}.longest_orfs.pep") into transdecoder_diamond, transdecoder_hmmer

            script:
                """
                cp ${assembly} ${sample_id}_asssembly.fasta

                echo -e "\\n-- TransDecoder.LongOrfs... --\\n"

                TransDecoder.LongOrfs -t ${assembly} --output_dir ${sample_id}.transdecoder_dir

                cp ${sample_id}.transdecoder_dir/longest_orfs.pep ${sample_id}.longest_orfs.pep

                echo -e "\\n-- Done with TransDecoder.LongOrfs --\\n"
                """
        }

        process transdecoder_diamond {

            label 'med_cpus'

            tag "${sample_id}"

            publishDir "${launchDir}/${params.outdir}/transdecoder", mode: "copy", overwrite: true

            conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::diamond=0.9.30=h56fc30b_0" : null)
            if (params.oneContainer){ container "${params.TPcontainer}" } else {
            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/diamond:0.9.30--h56fc30b_0" : "quay.io/biocontainers/diamond:0.9.30--h56fc30b_0")
            }

            input:
                tuple sample_id, file(pep) from transdecoder_diamond

            output:
                tuple sample_id, file("${sample_id}.diamond_blastp.outfmt6") into transdecoder_predict_diamond

            script:
                """
                dbPATH=${params.pipeInstall}/DBs/diamonddb_custom/

                echo -e "\\n-- Starting Diamond (blastp) --\\n"
                if [ ! -d \${dbPATH} ];then
                    echo "Directory \${dbPATH} not found. Run the precheck to fix this issue"
                    exit 0
                elif [ -d \${dbPATH} ];then
                    if [ ! -e \${dbPATH}/${params.uniname} ];then
                        echo "File \${dbPATH}/${params.uniname} not found. Run the precheck to fix this issue"
                        exit 0
                    elif [ -e \${dbPATH}/${params.uniname} ];then
                        if [ ! -e \${dbPATH}/${params.uniname}.dmnd ];then
                            cp \${dbPATH}/${params.uniname} .
                            diamond makedb --in ${params.uniname} -d ${params.uniname} -p ${task.cpus}
                            diamond blastp -d ${params.uniname}.dmnd -q ${pep} -p ${task.cpus} -f 6 -k 1 -e 0.00001 >${sample_id}.diamond_blastp.outfmt6
                        elif [ -e \${dbPATH}/${params.uniname}.dmnd ];then
                            cp \${dbPATH}/${params.uniname}.dmnd .
                            diamond blastp -d ${params.uniname}.dmnd -q ${pep} -p ${task.cpus} -f 6 -k 1 -e 0.00001 >${sample_id}.diamond_blastp.outfmt6
                        fi
                    fi
                fi
                echo -e "\\n-- Done with Diamond (blastp) --\\n"
                """
        }

        process transdecoder_hmmer {

            label 'med_cpus'

            tag "${sample_id}"

            publishDir "${launchDir}/${params.outdir}/transdecoder", mode: "copy", overwrite: true

            conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::hmmer=3.3=he1b5a44_0" : null)
            if (params.oneContainer){ container "${params.TPcontainer}" } else {
            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/hmmer:3.3--he1b5a44_0" : "quay.io/biocontainers/hmmer:3.3--he1b5a44_0")
            }

            input:
                tuple sample_id, file(pep) from transdecoder_hmmer

            output:
                tuple sample_id, file("${sample_id}.pfam.domtblout") into transdecoder_predict_hmmer

            script:
                """
                dbPATH=${params.pipeInstall}/DBs/hmmerdb/
                echo -e "\\n-- Starting HMMER --\\n"
                if [ ! -d \${dbPATH} ];then
                    echo "Directory \${dbPATH} not found. Run the precheck to fix this issue"
                    exit 0
                elif [ -d \${dbPATH} ];then
                    if [ ! -e \${dbPATH}/${params.pfname} ];then
                        echo "File \${dbPATH}/${params.pfname} not found. Run the precheck to fix this issue"
                        exit 0
                    elif [ -e \${dbPATH}/${params.pfname} ];then
                        if [ ! -e \${dbPATH}/${params.pfname}.h3f ] && [ ! -e \${dbPATH}/${params.pfname}.h3i ] && [ ! -e \${dbPATH}/${params.pfname}.h3m ] && [ ! -e \${dbPATH}/${params.pfname}.h3p ];then
                            cp \${dbPATH}/${params.pfname} .
                            hmmpress ${params.pfname}
                            hmmscan --cpu ${task.cpus} --domtblout ${sample_id}.pfam.domtblout ${params.pfname} ${pep}
                        elif [ -s \${dbPATH}/${params.pfname}.h3f ] && [ -s \${dbPATH}/${params.pfname}.h3i ] && [ -s \${dbPATH}/${params.pfname}.h3m ] && [ -s \${dbPATH}/${params.pfname}.h3p ];then
                            cp \${dbPATH}/${params.pfname}.* .
                            hmmscan --cpu ${task.cpus} --domtblout ${sample_id}.pfam.domtblout ${params.pfname} ${pep}
                        else
                            cp \${dbPATH}/${params.pfname} .
                            hmmpress ${params.pfname}
                            hmmscan --cpu ${task.cpus} --domtblout ${sample_id}.pfam.domtblout ${params.pfname} ${pep}
                        fi
                    fi
                fi
                echo -e "\\n-- Done with HMMER --\\n"
                """
        }

        transdecoder_predict_ch=Channel.create()
        transdecoder_predict_diamond.mix( transdecoder_predict_hmmer, evigene_ch_transdecoderB ).groupTuple(by:0,size:3).into(transdecoder_predict_ch)

        process transdecoder_predict {

            label 'med_cpus'

            tag "${sample_id}"

            publishDir "${launchDir}/${params.outdir}/transdecoder", mode: "copy", overwrite: true

            conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::transdecoder=5.5.0=pl526_2" : null)
            if (params.oneContainer){ container "${params.TPcontainer}" } else {
            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/transdecoder:5.5.0--pl526_2" : "quay.io/biocontainers/transdecoder:5.5.0--pl526_2")
            }

            input:
                tuple sample_id, file(files) from transdecoder_predict_ch

            output:
                tuple sample_id, file("${sample_id}.combined.okay.fa.transdecoder.pep") into ( transdecoder_ch_hmmer, transdecoder_ch_signalp, transdecoder_ch_tmhmm, transdecoder_ch_trinotate )
                tuple sample_id, file("${assembly}"), file("${sample_id}.combined.okay.fa.transdecoder.pep") into ( transdecoder_ch_diamond, transdecoder_ch_diamond_custom )
                tuple sample_id, file("${sample_id}_transdecoder.stats") into transdecoder_summary
                tuple sample_id, file("${sample_id}_transdecoder.csv") into transdecoder_csv
                tuple sample_id, file("${sample_id}*.transdecoder.{cds,gff,bed}") into transdecoder_files

            script:
                """
                ass=\$( echo $files | tr " " "\\n" | grep ".combined.okay.fa" )
                dia=\$( echo $files | tr " " "\\n" | grep ".diamond_blastp.outfmt6" )
                pfa=\$( echo $files | tr " " "\\n" | grep ".pfam.domtblout" )

                echo -e "\\n-- TransDecoder.LongOrfs... --\\n"

                TransDecoder.LongOrfs -t \${ass} --output_dir ${sample_id}.transdecoder_dir

                echo -e "\\n-- Done with TransDecoder.LongOrfs --\\n"

                echo -e "\\n-- TransDecoder.Predict... --\\n"

                TransDecoder.Predict -t \${ass} --retain_pfam_hits \${pfa} --retain_blastp_hits \${dia} --output_dir ${sample_id}.transdecoder_dir

                echo -e "\\n-- Done with TransDecoder.Predict --\\n"

                echo -e "\\n-- Calculating statistics... --\\n"
                #Calculate statistics of Transdecoder
                echo "- Transdecoder (long, with homology) stats for ${sample_id}" >${sample_id}_transdecoder.stats
                orfnum=\$( cat ${sample_id}*.transdecoder.pep | grep -c ">" )
                echo -e "Total number of ORFs: \$orfnum \\n" >>${sample_id}_transdecoder.stats
                echo -e "\\t Of these ORFs" >>${sample_id}_transdecoder.stats
                orfnum=\$( cat ${sample_id}*.transdecoder.pep | grep ">" | grep -c "|" )
                echo -e "\\t\\t with annotations: \$orfnum \\n" >>${sample_id}_transdecoder.stats
                orfnum=\$( cat ${sample_id}*.transdecoder.pep | grep ">" | grep -v "|" | grep -c ">" )
                echo -e "\\t\\t no annotation: \$orfnum \\n" >>${sample_id}_transdecoder.stats
                orfnum=\$( cat ${sample_id}*.transdecoder.pep | grep -c "ORF type:complete" )
                echo -e "\\t ORFs type=complete: \$orfnum \\n" >>${sample_id}_transdecoder.stats
                orfnum=\$( cat ${sample_id}*.transdecoder.pep | grep "ORF type:complete" | grep -c "|" )
                echo -e "\\t\\t with annotations: \$orfnum \\n" >>${sample_id}_transdecoder.stats
                orfnum=\$( cat ${sample_id}*.transdecoder.pep | grep -c "ORF type:5prime_partial" )
                echo -e "\\t ORFs type=5prime_partial: \$orfnum \\n" >>${sample_id}_transdecoder.stats
                orfnum=\$( cat ${sample_id}*.transdecoder.pep | grep "ORF type:5prime_partial" | grep -c "|" )
                echo -e "\\t\\t with annotations: \$orfnum \\n" >>${sample_id}_transdecoder.stats
                orfnum=\$( cat ${sample_id}*.transdecoder.pep | grep -c "ORF type:3prime_partial" )
                echo -e "\\t ORFs type=3prime_partial: \$orfnum \\n" >>${sample_id}_transdecoder.stats
                orfnum=\$( cat ${sample_id}*.transdecoder.pep | grep "ORF type:3prime_partial" | grep -c "|" )
                echo -e "\\t\\t with annotations: \$orfnum \\n" >>${sample_id}_transdecoder.stats
                orfnum=\$( cat ${sample_id}*.transdecoder.pep | grep -c "ORF type:internal" )
                echo -e "\\t ORFs type=internal: \$orfnum \\n" >>${sample_id}_transdecoder.stats
                orfnum=\$( cat ${sample_id}*.transdecoder.pep | grep "ORF type:internal" | grep -c "|" )
                echo -e "\\t\\t with annotations: \$orfnum \\n" >>${sample_id}_transdecoder.stats
                # csv for report
                echo "Sample,Total_orf,orf_complete,orf_5prime_partial,orf_3prime_partial,orf_internal" >${sample_id}_transdecoder.csv
                total=\$( cat ${sample_id}*.transdecoder.pep  | grep -c ">" )
                complete=\$( cat ${sample_id}*.transdecoder.pep  | grep -c "ORF type:complete" )
                n5prime=\$( cat ${sample_id}*.transdecoder.pep  | grep -c "ORF type:5prime_partial" )
                n3prime=\$( cat ${sample_id}*.transdecoder.pep  | grep -c "ORF type:3prime_partial" )
                internal=\$( cat ${sample_id}*.transdecoder.pep  | grep -c "ORF type:internal" )
                echo "${sample_id},\${total},\${complete},\${n5prime},\${n3prime},\${internal}" >>${sample_id}_transdecoder.csv
                echo -e "\\n-- Done with statistics --\\n"

                cp ${assembly} ${assembly}.tmp
                rm ${assembly}
                mv ${assembly}.tmp ${assembly}

                echo -e "\\n-- DONE with TransDecoder --\\n"
                """
        }
    }

    process swiss_diamond_trinotate {

        label 'big_cpus'

        tag "${sample_id}"

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::diamond=0.9.30=h56fc30b_0" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/diamond:0.9.30--h56fc30b_0" : "quay.io/biocontainers/diamond:0.9.30--h56fc30b_0")
        }

        input:
            tuple sample_id, file(assembly), file(transdecoder) from transdecoder_ch_diamond

        output:
            tuple sample_id, file("${sample_id}.diamond_blastx.outfmt6") into trinotate_ch_diamondX
            tuple sample_id, file("${sample_id}.diamond_blastp.outfmt6") into trinotate_ch_diamondP

        script:
            """
            dbPATH=${params.pipeInstall}/DBs/diamonddb_swiss/

            echo -e "\\n-- Starting Diamond --\\n"
            if [ ! -d \${dbPATH} ];then
                echo "Directory \${dbPATH} not found. Run the precheck to fix this issue"
                exit 0
            elif [ -d \${dbPATH} ];then
                if [ ! -e \${dbPATH}/uniprot_sprot.pep ];then
                    echo "File \${dbPATH}/uniprot_sprot.pep not found. Run the precheck to fix this issue"
                    exit 0
                elif [ -e \${dbPATH}/uniprot_sprot.pep ];then
                    if [ ! -e \${dbPATH}/uniprot_sprot.pep.dmnd ];then
                        cp \${dbPATH}/uniprot_sprot.pep .
                        diamond makedb --in uniprot_sprot.pep -d uniprot_sprot.pep -p ${task.cpus}
                        echo -e "\\n-- Starting with Diamond (blastx) --\\n"
                        diamond blastx -d uniprot_sprot.pep.dmnd -q ${assembly} -p ${task.cpus} -f 6 -k 1 -e 0.001 >${sample_id}.diamond_blastx.outfmt6
                        echo -e "\\n-- Done with Diamond (blastx) --\\n"
                        echo -e "\\n-- Starting with Diamond (blastp) --\\n"
                        diamond blastp -d uniprot_sprot.pep.dmnd -q ${transdecoder} -p ${task.cpus} -f 6 -k 1 -e 0.001 >${sample_id}.diamond_blastp.outfmt6
                        echo -e "\\n-- Done with Diamond (blastp)  --\\n"
                    elif [ -e \${dbPATH}/uniprot_sprot.pep.dmnd ];then
                        cp \${dbPATH}/uniprot_sprot.pep.dmnd .
                        echo -e "\\n-- Starting with Diamond (blastx) --\\n"
                        diamond blastx -d uniprot_sprot.pep.dmnd -q ${assembly} -p ${task.cpus} -f 6 -k 1 -e 0.001 >${sample_id}.diamond_blastx.outfmt6
                        echo -e "\\n-- Done with Diamond (blastx) --\\n"
                        echo -e "\\n-- Starting with Diamond (blastp) --\\n"
                        diamond blastp -d uniprot_sprot.pep.dmnd -q ${transdecoder} -p ${task.cpus} -f 6 -k 1 -e 0.001 >${sample_id}.diamond_blastp.outfmt6
                        echo -e "\\n-- Done with Diamond (blastp)  --\\n"
                    fi
                fi
            fi
            echo -e "\\n-- Done with Diamond --\\n"
            """
    }

    process custom_diamond_trinotate {

        label 'big_cpus'

        tag "${sample_id}"

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::diamond=0.9.30=h56fc30b_0" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/diamond:0.9.30--h56fc30b_0" : "quay.io/biocontainers/diamond:0.9.30--h56fc30b_0")
        }

        input:
            tuple sample_id, file(assembly), file(transdecoder) from transdecoder_ch_diamond_custom

        output:
            tuple sample_id, file("${sample_id}.custom.diamond_blastx.outfmt6") into trinotate_ch_diamondX_custom
            tuple sample_id, file("${sample_id}.custom.diamond_blastp.outfmt6") into trinotate_ch_diamondP_custom

        script:
            """
            dbPATH=${params.pipeInstall}/DBs/diamonddb_custom/

            echo -e "\\n-- Starting Diamond --\\n"
            if [ ! -d \${dbPATH} ];then
                echo "Directory \${dbPATH} not found. Run the precheck to fix this issue"
                exit 0
            elif [ -d \${dbPATH} ];then
                if [ ! -e \${dbPATH}/${params.uniname} ];then
                    echo "File \${dbPATH}/${params.uniname} not found. Run the precheck to fix this issue"
                    exit 0
                elif [ -e \${dbPATH}/${params.uniname} ];then
                    if [ ! -e \${dbPATH}/${params.uniname}.dmnd ];then
                        cp \${dbPATH}/${params.uniname} .
                        diamond makedb --in ${params.uniname} -d ${params.uniname} -p ${task.cpus}
                        echo -e "\\n-- Starting with Diamond (blastx) --\\n"
                        diamond blastx -d ${params.uniname}.dmnd -q ${assembly} -p ${task.cpus} -f 6 -k 1 -e 0.001 >${sample_id}.custom.diamond_blastx.outfmt6
                        echo -e "\\n-- Done with Diamond (blastx) --\\n"
                        echo -e "\\n-- Starting with Diamond (blastp) --\\n"
                        diamond blastp -d ${params.uniname}.dmnd -q ${transdecoder} -p ${task.cpus} -f 6 -k 1 -e 0.001 >${sample_id}.custom.diamond_blastp.outfmt6
                        echo -e "\\n-- Done with Diamond (blastp)  --\\n"
                    elif [ -e \${dbPATH}/${params.uniname}.dmnd ];then
                        cp \${dbPATH}/${params.uniname}.dmnd .
                        echo -e "\\n-- Starting with Diamond (blastx) --\\n"
                        diamond blastx -d ${params.uniname}.dmnd -q ${assembly} -p ${task.cpus} -f 6 -k 1 -e 0.001 >${sample_id}.custom.diamond_blastx.outfmt6
                        echo -e "\\n-- Done with Diamond (blastx) --\\n"
                        echo -e "\\n-- Starting with Diamond (blastp) --\\n"
                        diamond blastp -d ${params.uniname}.dmnd -q ${transdecoder} -p ${task.cpus} -f 6 -k 1 -e 0.001 >${sample_id}.custom.diamond_blastp.outfmt6
                        echo -e "\\n-- Done with Diamond (blastp)  --\\n"
                    fi
                fi
            fi
            echo -e "\\n-- Done with Diamond --\\n"
            """
    }

    process hmmer_trinotate {

        label 'low_cpus'

        tag "${sample_id}"

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::hmmer=3.3=he1b5a44_0" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/hmmer:3.3--he1b5a44_0" : "quay.io/biocontainers/hmmer:3.3--he1b5a44_0")
        }

        input:
            tuple sample_id, file(transdecoder_pep) from transdecoder_ch_hmmer

        output:
            tuple sample_id, file("${sample_id}.TrinotatePFAM.out") into trinotate_ch_hmmer

        script:
            """
            dbPATH=${params.pipeInstall}/DBs/hmmerdb/

            echo -e "\\n-- Starting HMMER --\\n"
            if [ ! -d \${dbPATH} ];then
                echo "Directory \${dbPATH} not found. Run the precheck to fix this issue"
                exit 0
            elif [ -d \${dbPATH} ];then
                if [ ! -e \${dbPATH}/${params.pfname} ];then
                    echo "File \${dbPATH}/${params.pfname} not found. Run the precheck to fix this issue"
                    exit 0
                elif [ -e \${dbPATH}/${params.pfname} ];then
                    if [ ! -e \${dbPATH}/${params.pfname}.h3f ] && [ ! -e \${dbPATH}/${params.pfname}.h3i ] && [ ! -e \${dbPATH}/${params.pfname}.h3m ] && [ ! -e \${dbPATH}/${params.pfname}.h3p ];then
                        cp \${dbPATH}/${params.pfname} .
                        hmmpress ${params.pfname}
                        hmmscan --cpu ${task.cpus} --domtblout ${sample_id}.TrinotatePFAM.out ${params.pfname} ${transdecoder_pep} >pfam.log
                    elif [ -s \${dbPATH}/${params.pfname}.h3f ] && [ -s \${dbPATH}/${params.pfname}.h3i ] && [ -s \${dbPATH}/${params.pfname}.h3m ] && [ -s \${dbPATH}/${params.pfname}.h3p ];then
                        cp \${dbPATH}/${params.pfname}.* .
                        hmmscan --cpu ${task.cpus} --domtblout ${sample_id}.TrinotatePFAM.out ${params.pfname} ${transdecoder_pep} >pfam.log
                    else
                        cp \${dbPATH}/${params.pfname} .
                        hmmpress ${params.pfname}
                        hmmscan --cpu ${task.cpus} --domtblout ${sample_id}.TrinotatePFAM.out ${params.pfname} ${transdecoder_pep} >pfam.log
                    fi
                fi
            fi
            echo -e "\\n-- Done with HMMER --\\n"
            """
    }

    if (params.withSignalP) {

        process signalP_trinotate {

            label 'low_cpus'

            tag "${sample_id}"

            input:
                tuple sample_id, file(transdecoder_pep) from transdecoder_ch_signalp

            output:
                tuple sample_id, file("${sample_id}.signalp.out") into trinotate_ch_signalp

            script:
                """
                #signalP to predict signal peptides

                echo -e "\\n-- Starting with SignalP --\\n"

                ${params.signalp} -f short -n ${sample_id}.signalp.out ${transdecoder_pep}

                echo -e "\\n-- Done with SignalP --\\n"
                """
        }
    } else {

        process skip_signalP {
            tag "${sample_id}"

            input:
                tuple sample_id, file(transdecoder_pep) from transdecoder_ch_signalp

            output:
                tuple sample_id, file("${sample_id}.signalp.out") into trinotate_ch_signalp

            script:
                """
                touch ${sample_id}.signalp.out
                """
        }
    }

    if (params.withTMHMM) {

        process tmhmm_trinotate {

            label 'low_cpus'

            tag "${sample_id}"

            input:
                tuple sample_id, file(transdecoder_pep) from transdecoder_ch_tmhmm

            output:
                tuple sample_id, file("${sample_id}.tmhmm.out") into trinotate_ch_tmhmm

            script:
                """
                #tmHMM to predict transmembrane regions

                echo -e "\\n-- Starting with tmHMM --\\n"

                ${params.tmhmm} --short < ${transdecoder_pep} >${sample_id}.tmhmm.out

                echo -e "\\n-- Done with tmHMM --\\n"
                """
        }
    } else {

        process skip_tmhmm {
            tag "${sample_id}"

            input:
                tuple sample_id, file(transdecoder_pep) from transdecoder_ch_tmhmm

            output:
                tuple sample_id, file("${sample_id}.tmhmm.out") into trinotate_ch_tmhmm

            script:
                """
                touch ${sample_id}.tmhmm.out
                """
        }
    }

    if (params.withRnammer) {

        process rnammer_trinotate {

            label 'low_cpus'

            tag "${sample_id}"

            input:
                tuple sample_id, file(transcriptome) from evigene_ch_rnammer

            output:
                tuple sample_id, file("${sample_id}.rnammer.gff") into trinotate_ch_rnammer

            script:
                """
                set +e
                #RNAMMER to identify rRNA transcripts

                echo -e "\\n-- Starting with RNAMMER --\\n"

                RnammerTranscriptome.pl --transcriptome ${transcriptome} --path_to_rnammer ${params.rnam}

                mv ${sample_id}.combined.okay.fa.rnammer.gff ${sample_id}.rnammer.gff

                echo -e "\\n-- Done with RNAMMER --\\n"
                """
        }
    } else {

        process skip_rnammer {
            tag "${sample_id}"

            input:
                tuple sample_id, file(transcriptome) from evigene_ch_rnammer

            output:
                tuple sample_id, file("${sample_id}.rnammer.gff") into trinotate_ch_rnammer

            script:
                """
                touch ${sample_id}.rnammer.gff
                """
        }
    }

    trinotate_ch = Channel.create()
    evigene_ch_trinotate.mix( transdecoder_ch_trinotate,trinotate_ch_diamondX,trinotate_ch_diamondP,trinotate_ch_diamondX_custom,trinotate_ch_diamondP_custom,trinotate_ch_hmmer,trinotate_ch_signalp,trinotate_ch_tmhmm,trinotate_ch_rnammer ).groupTuple(by:0,size:10).into(trinotate_ch)

    process trinotate {

        label 'low_cpus'

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/trinotate", mode: "copy", overwrite: true

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::trinotate=3.2.1=pl526_0" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/trinotate:3.2.1--pl526_0" : "quay.io/biocontainers/trinotate:3.2.1--pl526_0")
        }

        input:
            tuple sample_id, file(files) from trinotate_ch

        output:
            tuple sample_id, file("${sample_id}.GO.terms.txt") into trinotate_summary
            tuple sample_id, file("${sample_id}.trinotate_annotation_report.xls") into ( trinotate_out_ch, custom_uniprot_ch )
            tuple sample_id, file("*.terms.txt") into other_files
            tuple sample_id, file("${sample_id}.KEGG.terms.txt") into kegg_paths

        script:
            """
            for x in `echo ${files}`;do
                echo \${x} >>.vars.txt
            done

            assembly=\$( cat .vars.txt | grep "${sample_id}.combined.okay.fa" | grep -v "${sample_id}.combined.okay.fa.transdecoder.pep" )
            transdecoder=\$( cat .vars.txt | grep -E "${sample_id}.*.transdecoder.pep" )
            diamond_blastx=\$( cat .vars.txt | grep "${sample_id}.diamond_blastx.outfmt6" )
            diamond_blastp=\$( cat .vars.txt | grep "${sample_id}.diamond_blastp.outfmt6" )
            custom_blastx=\$( cat .vars.txt | grep "${sample_id}.custom.diamond_blastx.outfmt6" )
            custom_blastp=\$( cat .vars.txt | grep "${sample_id}.custom.diamond_blastp.outfmt6" )
            pfam=\$( cat .vars.txt | grep "${sample_id}.TrinotatePFAM.out" )
            signalp=\$( cat .vars.txt | grep "${sample_id}.signalp.out" )
            tmhmm=\$( cat .vars.txt | grep "${sample_id}.tmhmm.out" )
            rnammer=\$( cat .vars.txt | grep "${sample_id}.rnammer.gff" )

            #Generate gene_trans_map
            #Not using get_Trinity_gene_to_trans_map.pl since all the names are uniq
            cat \${assembly} | awk '{print \$1}' | grep ">" | cut -c 2- >a.txt

            paste a.txt a.txt >\${assembly}.gene_trans_map

            #Get Trinotate.sqlite from folder (original)
            cp ${params.Tsql} .
            sqlname=`echo ${params.Tsql} | tr "\\/" "\\n" | grep "\\.sqlite"`

            echo -e "\\n-- Running Trinotate --\\n"

            Trinotate \$sqlname init --gene_trans_map \${assembly}.gene_trans_map --transcript_fasta \${assembly} --transdecoder_pep \${transdecoder}

            echo -e "\\n-- Ending run of Trinotate --\\n"

            echo -e "\\n-- Loading hits and predictions to sqlite database... --\\n"

            #Load protein hits
            Trinotate \$sqlname LOAD_swissprot_blastp \${diamond_blastp}

            #Load transcript hits
            Trinotate \$sqlname LOAD_swissprot_blastx \${diamond_blastx}

            #Load custom protein hits
            Trinotate \$sqlname LOAD_custom_blast --outfmt6 \${custom_blastp} --prog blastp --dbtype ${params.uniname}

            #Load custom transcript hits
            Trinotate \$sqlname LOAD_custom_blast --outfmt6 \${custom_blastx} --prog blastx --dbtype ${params.uniname}

            #Load Pfam domain entries
            Trinotate \$sqlname LOAD_pfam \${pfam}

            #Load transmembrane domains
            if [ -s \${tmhmm} ];then
                Trinotate \$sqlname LOAD_tmhmm \${tmhmm}
            else
                echo "No transmembrane domains (tmhmm)"
            fi

            #Load signal peptide predictions
            if [ -s \${signalp} ];then
                Trinotate \$sqlname LOAD_signalp \${signalp}
            else
                echo "No Signal-P"
            fi

            #Load rnammer results
            if [ -s \${rnammer} ];then
                Trinotate \$sqlname LOAD_rnammer \${rnammer}
            else
                echo "No rnammer results"
            fi

            echo -e "\\n-- Loading finished --\\n"

            #Report

            echo -e "\\n-- Generating report... --\\n"

            Trinotate \$sqlname report >${sample_id}.trinotate_annotation_report.xls

            echo -e "\\n-- Report generated --\\n"

            #Extract info from XLS file

            echo -e "\\n-- Creating GO file from XLS... --\\n"

            extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls ${sample_id}.trinotate_annotation_report.xls --trans >${sample_id}.GO.terms.txt

            echo -e "\\n-- Done with the GO --\\n"

            echo -e "\\n-- Creating KEGG file from XLS... --\\n"

            cat ${sample_id}.trinotate_annotation_report.xls | cut -f 1,14 | grep "KEGG" | tr "\\`" ";" | grep "KO:K" | sed 's/\\tKEGG/\\t#KEGG/g' | sed 's/KO:/KO:#/g' | cut -f 1,3 -d "#" | tr -d "#" >${sample_id}.KEGG.terms.txt

            echo -e "\\n-- Done with the KEGG --\\n"

            echo -e "\\n-- Creating eggNOG file from XLS... --\\n"

            cat ${sample_id}.trinotate_annotation_report.xls | cut -f 1,13 | grep "OG" | tr "\\`" ";" | sed 's/^/#/g' | sed 's/;/\\n;/g' | cut -f 1 -d "^" | tr -d "\\n" | tr "#" "\\n" | grep "OG" >${sample_id}.eggNOG_COG.terms.txt

            echo -e "\\n-- Done with the eggNOG --\\n"

            echo -e "\\n-- Creating PFAM file from XLS... --\\n"

            cat ${sample_id}.trinotate_annotation_report.xls | cut -f 1,10 | grep "PF" | tr "\\`" ";" | sed 's/^/#/g' | sed 's/;PF/\\n;PF/g' | cut -f 1 -d "^" | tr -d "\\n" | tr "#" "\\n" | grep "PF" | tr ";" "," >${sample_id}.PFAM.terms.txt

            echo -e "\\n-- Done with the PFAM --\\n"

            echo -e "\\n-- DONE with Trinotate --\\n"
            """
    }

    process summary_evigene_individual {

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/stats", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.combined.fa"), file("${sample_id}.combined.okay.fa") from evigene_summary

        output:
            tuple sample_id, file("${sample_id}_sum_preEG.txt"), file("${sample_id}_sum_EG.txt") into final_sum_1
            tuple sample_id, file("*.csv") into summary_evi_csv

        script:
            """
            #Summary of total number of transcripts
            echo -e "- Number of transcripts before Evidential Genes\\n" >>${sample_id}_sum_preEG.txt
            echo -e "- Individual ${sample_id} \\n" >>${sample_id}_sum_preEG.txt
            echo -e "\\t Total transcripts:" >>${sample_id}_sum_preEG.txt
            num=\$( cat ${sample_id}.combined.fa | grep -c ">" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}_sum_preEG.txt
            echo -e "\\t Trinity" >>${sample_id}_sum_preEG.txt
            num=\$( cat ${sample_id}.combined.fa | grep -c ">TRINITY" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}_sum_preEG.txt
            echo -e "\\t SOAP" >>${sample_id}_sum_preEG.txt
            num=\$( cat ${sample_id}.combined.fa | grep -c ">SOAP" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}_sum_preEG.txt
            echo -e "\\t Velvet/Oases" >>${sample_id}_sum_preEG.txt
            num=\$( cat ${sample_id}.combined.fa | grep -c ">Velvet" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}_sum_preEG.txt
            echo -e "\\t rna-SPADES" >>${sample_id}_sum_preEG.txt
            num=\$( cat ${sample_id}.combined.fa | grep -c ">SPADES" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}_sum_preEG.txt
            echo -e "\\t Trans-ABySS" >>${sample_id}_sum_preEG.txt
            num=\$( cat ${sample_id}.combined.fa | grep -c ">TransABySS" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}_sum_preEG.txt

            # csv report
            echo "Total,Trinity,SOAP,Velvet,SPADES,TransBySS" >${sample_id}_sum_preEG.csv
            total=\$( cat ${sample_id}.combined.fa | grep -c ">" )
            trinity=\$( cat ${sample_id}.combined.fa | grep -c ">TRINITY" )
            soap=\$( cat ${sample_id}.combined.fa | grep -c ">SOAP" )
            velvet=\$( cat ${sample_id}.combined.fa | grep -c ">Velvet" )
            spades=\$( cat ${sample_id}.combined.fa | grep -c ">SPADES" )
            transabyss=\$( cat ${sample_id}.combined.fa | grep -c ">TransABySS" )
            echo "\${total},\${trinity},\${soap},\${velvet},\${spades},\${transabyss}" >>${sample_id}_sum_preEG.csv

            #Summary of transcripts after EvidentialGenes
            echo -e "- Number of transcripts by individual after EvidentialGenes\\n" >>${sample_id}_sum_EG.txt
            echo -e "- Individual ${sample_id} \\n" >>${sample_id}_sum_EG.txt
            echo -e "\\t Total transcripts:" >>${sample_id}_sum_EG.txt
            num=\$( cat ${sample_id}.combined.okay.fa | grep -c ">" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}_sum_EG.txt
            echo -e "\\t Trinity" >>${sample_id}_sum_EG.txt
            num=\$( cat ${sample_id}.combined.okay.fa | grep -c ">TRINITY" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}_sum_EG.txt
            echo -e "\\t SOAP" >>${sample_id}_sum_EG.txt
            num=\$( cat ${sample_id}.combined.okay.fa | grep -c ">SOAP" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}_sum_EG.txt
            echo -e "\\t Velvet/Oases" >>${sample_id}_sum_EG.txt
            num=\$( cat ${sample_id}.combined.okay.fa | grep -c ">Velvet" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}_sum_EG.txt
            echo -e "\\t rna-SPADES" >>${sample_id}_sum_EG.txt
            num=\$( cat ${sample_id}.combined.okay.fa | grep -c ">SPADES" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}_sum_EG.txt
            echo -e "\\t Trans-ABySS" >>${sample_id}_sum_EG.txt
            num=\$( cat ${sample_id}.combined.okay.fa | grep -c ">TransABySS" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}_sum_EG.txt

            # csv report after evigene
            echo "Total,Trinity,SOAP,Velvet,SPADES,TransBySS" >${sample_id}_sum_EG.csv
            total=\$( cat ${sample_id}.combined.okay.fa | grep -c ">" )
            trinity=\$( cat ${sample_id}.combined.okay.fa | grep -c ">TRINITY" )
            soap=\$( cat ${sample_id}.combined.okay.fa | grep -c ">SOAP" )
            velvet=\$( cat ${sample_id}.combined.okay.fa | grep -c ">Velvet" )
            spades=\$( cat ${sample_id}.combined.okay.fa | grep -c ">SPADES" )
            transabyss=\$( cat ${sample_id}.combined.okay.fa | grep -c ">TransABySS" )
            echo "\${total},\${trinity},\${soap},\${velvet},\${spades},\${transabyss}" >>${sample_id}_sum_EG.csv
            """
    }

    busco3_sum = Channel.create()
    busco3_summary.mix(busco3_ch_trinity_sum).groupTuple(by:0,size:2).into(busco3_sum)

    process summary_busco3_individual {

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/stats", mode: "copy", overwrite: true

        input:
            tuple sample_id, file(files) from busco3_sum

        output:
            tuple sample_id, file("${sample_id}.sum_busco3.txt") into final_sum_2v3

        script:
            """
            tri=\$( echo $files | tr " " "\\n" | grep ".Trinity.bus3.txt" )
            trans=\$( echo $files | tr " " "\\n" | grep ".TransPi.bus3.txt" )
            #Summary of BUSCO scores for the final_assemblies
            echo -e "Summary of BUSCO V3 \\n" >>${sample_id}.sum_busco3.txt
            echo "-- TransPi BUSCO V3 scores -- " >>${sample_id}.sum_busco3.txt
            cat \${trans} >>${sample_id}.sum_busco3.txt
            echo -e "\\n-- Trinity BUSCO V3 scores --" >>${sample_id}.sum_busco3.txt
            cat \${tri} >>${sample_id}.sum_busco3.txt
            """
    }

    busco4_sum = Channel.create()
    busco4_summary.mix(busco4_ch_trinity_sum).groupTuple(by:0,size:2).into(busco4_sum)

    process summary_busco4_individual {

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/stats", mode: "copy", overwrite: true

        input:
            tuple sample_id, file(files) from busco4_sum

        output:
            tuple sample_id, file("${sample_id}.sum_busco4.txt") into final_sum_2v4

        script:
            """
            tri=\$( echo $files | tr " " "\\n" | grep ".Trinity.bus4.txt" )
            trans=\$( echo $files | tr " " "\\n" | grep ".TransPi.bus4.txt" )
            #Summary of BUSCO scores for the final_assemblies
            echo -e "Summary of BUSCO V4 \\n" >>${sample_id}.sum_busco4.txt
            echo "-- TransPi BUSCO V4 scores -- " >>${sample_id}.sum_busco4.txt
            cat \${trans} >>${sample_id}.sum_busco4.txt
            echo -e "\\n-- Trinity BUSCO V4 scores --" >>${sample_id}.sum_busco4.txt
            cat \${tri} >>${sample_id}.sum_busco4.txt
            """
    }

    process summary_transdecoder_individual {

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/stats", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}_transdecoder.stats") from transdecoder_summary

        output:
            tuple sample_id, file("${sample_id}.sum_transdecoder.txt") into final_sum_3

        script:
            """
            #Summary of Transdecoder stats
            echo -e "Summary of Transdecoder \\n" >>${sample_id}.sum_transdecoder.txt
            cat ${sample_id}_transdecoder.stats >>${sample_id}.sum_transdecoder.txt
            echo -e "##### \\n" >>${sample_id}.sum_transdecoder.txt
            """
    }

    process summary_trinotate_individual {

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/stats", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.GO.terms.txt") from trinotate_summary

        output:
            tuple sample_id, file("${sample_id}.sum_GO.txt") into final_sum_4

        script:
            """
            #Summary of Trinotate (Gene Ontologies)
            echo -e "Summary of Trinotate/Gene Ontologies \\n" >>${sample_id}.sum_GO.txt
            echo "- Individual "${sample_id} >>${sample_id}.sum_GO.txt
            echo -e "\\t Total transcripts with GO:" >>${sample_id}.sum_GO.txt
            num=\$( cat ${sample_id}.GO.terms.txt | wc -l )
            echo -e "\\t\\t \$num" >>${sample_id}.sum_GO.txt
            tnum=\$num

            echo -e "\\t Total transcripts with only one GO:" >>${sample_id}.sum_GO.txt
            a=0
            while read lines;do
                if [[ `echo \$lines | awk '{print \$2}' | tr "," "\\n" | wc -l` -eq 1 ]];then
                    a=\$((a+1))
                fi
            done <${sample_id}.GO.terms.txt
            num=\$a
            echo -e "\\t\\t \$num" >>${sample_id}.sum_GO.txt
            onum=\$num

            echo -e "\\t Total transcripts with multiple GO:" >>${sample_id}.sum_GO.txt
            num=\$( echo \$tnum-\$onum | bc )
            echo -e "\\t\\t \$num" >>${sample_id}.sum_GO.txt

            echo -e "\\t Total GO in the file:" >>${sample_id}.sum_GO.txt
            num=\$( cat ${sample_id}.GO.terms.txt | awk '{print \$2}' | tr "," "\\n" | wc -l )
            echo -e "\\t\\t \$num" >>${sample_id}.sum_GO.txt

            echo -e "\\t Total uniq GO in the file:" >>${sample_id}.sum_GO.txt
            num=\$( cat ${sample_id}.GO.terms.txt | awk '{print \$2}' | tr "," "\\n" | sort --parallel=10 | uniq | wc -l )
            echo -e "\\t\\t \$num" >>${sample_id}.sum_GO.txt
            """
    }

    busco3_comp = Channel.create()
    busco3_comp_1.mix(busco3_comp_2).groupTuple(by:0,size:2).into(busco3_comp)

    process get_busco3_comparison {

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/figures/BUSCO3", mode: "copy", overwrite: true

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge -c bioconda r-reshape2=1.4.4 r-plotly=4.9.2.1 plotly-orca=3.4.2 r-ggplot2=3.3.0 r-svglite=1.2.3 r-ggthemes=4.2.0 r-knitr=1.29 r-rmarkdown=2.3 r-kableextra=1.1.0" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/mulled-v2-3f431f5f8e54df68ea0029c209fce3b154f6e186:94cad00b5306639ceab6aaf211f45740560abb90-0" : "quay.io/biocontainers/mulled-v2-3f431f5f8e54df68ea0029c209fce3b154f6e186:94cad00b5306639ceab6aaf211f45740560abb90-0")
        }

        input:
            tuple sample_id, file(files) from busco3_comp

        output:
            tuple sample_id, file("${sample_id}_BUSCO3_comparison.pdf"), file("${sample_id}_BUSCO3_comparison.svg") into busco3_fig
            tuple sample_id, file("*.csv") into busco3_csv

        script:
            """
            set +e
            tri=\$( echo $files | tr " " "\\n" | grep ".Trinity.bus3.txt" )
            trans=\$( echo $files | tr " " "\\n" | grep ".TransPi.bus3.txt" )
            bash get_busco_val.sh \${tri} \${trans} v3 ${sample_id}
            cp ${params.pipeInstall}/bin/busco_comparison.R .
            a=\$( cat final_spec )
            sed -i "s/MYSPEC/\${a}/" busco_comparison.R
            b=\$( cat final_perc )
            sed -i "s/MYPERC/\${b}/" busco_comparison.R
            c=\$( cat final_num )
            sed -i "s/MYVAL/\${c}/" busco_comparison.R
            Rscript busco_comparison.R ${sample_id}
            mv ${sample_id}_BUSCO_comparison.pdf ${sample_id}_BUSCO3_comparison.pdf
            mv ${sample_id}_BUSCO_comparison.svg ${sample_id}_BUSCO3_comparison.svg
            # csv
            sed -i 's/\$/\\n/g' final_*
            cat final_spec final_perc final_num | tr -d "'" >${sample_id}_busco3.csv
            """
    }

    busco4_comp = Channel.create()
    busco4_comp_1.mix(busco4_comp_2).groupTuple(by:0,size:2).into(busco4_comp)

    process get_busco4_comparison {

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/figures/BUSCO4", mode: "copy", overwrite: true

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge -c bioconda r-reshape2=1.4.4 r-plotly=4.9.2.1 plotly-orca=3.4.2 r-ggplot2=3.3.0 r-svglite=1.2.3 r-ggthemes=4.2.0 r-knitr=1.29 r-rmarkdown=2.3 r-kableextra=1.1.0" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/mulled-v2-3f431f5f8e54df68ea0029c209fce3b154f6e186:94cad00b5306639ceab6aaf211f45740560abb90-0" : "quay.io/biocontainers/mulled-v2-3f431f5f8e54df68ea0029c209fce3b154f6e186:94cad00b5306639ceab6aaf211f45740560abb90-0")
        }

        input:
            tuple sample_id, file(files) from busco4_comp

        output:
            tuple sample_id, file("${sample_id}_BUSCO4_comparison.pdf"), file("${sample_id}_BUSCO4_comparison.svg") into busco4_fig
            tuple sample_id, file("*.csv") into busco4_csv

        script:
            """
            set +e
            tri=\$( echo $files | tr " " "\\n" | grep ".Trinity.bus4.txt" )
            trans=\$( echo $files | tr " " "\\n" | grep ".TransPi.bus4.txt" )
            bash get_busco_val.sh \${tri} \${trans} v4 ${sample_id}
            cp ${params.pipeInstall}/bin/busco_comparison.R .
            a=\$( cat final_spec )
            sed -i "s/MYSPEC/\${a}/" busco_comparison.R
            b=\$( cat final_perc )
            sed -i "s/MYPERC/\${b}/" busco_comparison.R
            c=\$( cat final_num )
            sed -i "s/MYVAL/\${c}/" busco_comparison.R
            Rscript busco_comparison.R ${sample_id}
            mv ${sample_id}_BUSCO_comparison.pdf ${sample_id}_BUSCO4_comparison.pdf
            mv ${sample_id}_BUSCO_comparison.svg ${sample_id}_BUSCO4_comparison.svg
            # csv
            sed -i 's/\$/\\n/g' final_*
            cat final_spec final_perc final_num | tr -d "'" >${sample_id}_busco4.csv
            """
    }

    process get_GO_comparison {

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/figures/GO", mode: "copy", overwrite: true

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge -c bioconda r-reshape2=1.4.4 r-plotly=4.9.2.1 plotly-orca=3.4.2 r-ggplot2=3.3.0 r-svglite=1.2.3 r-ggthemes=4.2.0 r-knitr=1.29 r-rmarkdown=2.3 r-kableextra=1.1.0" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/mulled-v2-3f431f5f8e54df68ea0029c209fce3b154f6e186:94cad00b5306639ceab6aaf211f45740560abb90-0" : "quay.io/biocontainers/mulled-v2-3f431f5f8e54df68ea0029c209fce3b154f6e186:94cad00b5306639ceab6aaf211f45740560abb90-0")
        }

        input:
            tuple sample_id, file("${sample_id}.trinotate_annotation_report.xls") from trinotate_out_ch

        output:
            tuple sample_id, file("*.svg"), file("*.pdf") , file("*.txt") into go_fig
            tuple sample_id, file("*.csv") into go_csv

        script:
            """
            set +e
	        cp ${params.pipeInstall}/bin/GO_plots.R .

            cat ${sample_id}.trinotate_annotation_report.xls | awk 'FS="\\t",OFS="#" {print \$1,\$15,\$16,\$17}' | grep -v "gene_id" >all_GOs.txt

            touch final_GOs.txt

            while read line;do
                echo \${line} | cut -f 2,3,4 -d "#" | grep "GO:" | tr "#" "\\n" | tr "\\`" "\\n" | sed 's/\\. /,/g' | tr "," "\\n" | grep "GO:" | sort -u >>final_GOs.txt
            done<all_GOs.txt

            cat final_GOs.txt | tr [a-z] [A-Z] | grep "CELLULAR_COMPONENT" | cut -f 3 -d "^" | sort | uniq -c | sort -nr | head -n 15 | sed 's/^ *//g' \
            | sed 's/\\([0-9] \\)/\\1#/g' | tr "#" "\\t" >GO_cellular.txt

            cat final_GOs.txt | tr [a-z] [A-Z] | grep "BIOLOGICAL_PROCESS" | cut -f 3 -d "^" | sort | uniq -c | sort -nr | head -n 15 | sed 's/^ *//g' \
            | sed 's/\\([0-9] \\)/\\1#/g' | tr "#" "\\t" >GO_biological.txt

            cat final_GOs.txt | tr [a-z] [A-Z] | grep "MOLECULAR_FUNCTION" | cut -f 3 -d "^" | sort | uniq -c | sort -nr | head -n 15 | sed 's/^ *//g' \
            | sed 's/\\([0-9] \\)/\\1#/g' | tr "#" "\\t" >GO_molecular.txt

            Rscript GO_plots.R ${sample_id}

            rm all_GOs.txt final_GOs.txt

            mv GO_cellular.txt ${sample_id}_GO_cellular.txt
            mv GO_biological.txt ${sample_id}_GO_biological.txt
            mv GO_molecular.txt ${sample_id}_GO_molecular.txt

            cat ${sample_id}_GO_cellular.txt | sed -r 's/^[^0-9]*([0-9]+)/\\1,/g' >${sample_id}_GO_cellular.csv
            cat ${sample_id}_GO_biological.txt | sed -r 's/^[^0-9]*([0-9]+)/\\1,/g' >${sample_id}_GO_biological.csv
            cat ${sample_id}_GO_molecular.txt | sed -r 's/^[^0-9]*([0-9]+)/\\1,/g' >${sample_id}_GO_molecular.csv
            """
    }

    process summary_custom_uniprot {

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/figures/CustomUniProt", mode: "copy", overwrite: true

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge -c bioconda r-reshape2=1.4.4 r-plotly=4.9.2.1 plotly-orca=3.4.2 r-ggplot2=3.3.0 r-svglite=1.2.3 r-ggthemes=4.2.0 r-knitr=1.29 r-rmarkdown=2.3 r-kableextra=1.1.0" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/mulled-v2-3f431f5f8e54df68ea0029c209fce3b154f6e186:94cad00b5306639ceab6aaf211f45740560abb90-0" : "quay.io/biocontainers/mulled-v2-3f431f5f8e54df68ea0029c209fce3b154f6e186:94cad00b5306639ceab6aaf211f45740560abb90-0")
        }

        input:
            tuple sample_id, file("${sample_id}.trinotate_annotation_report.xls") from custom_uniprot_ch

        output:
            tuple sample_id, file("${sample_id}_custom_uniprot_hits.txt") into custom_uniprot_sum
            tuple sample_id, file("${sample_id}_custom_uniprot_hits.svg"), file("${sample_id}_custom_uniprot_hits.pdf") into custom_uniprot_fig
            tuple sample_id, file("*.csv") into uniprot_csv

        script:
            """
            #get custom blast hits
            cat ${sample_id}.trinotate_annotation_report.xls | cut -f 8 | grep [A-Z] | grep "|" | tr "\\`" "\\n" | \
                cut -f 1 -d "^" | cut -f 3 -d "|" | cut -f 2 -d "_" >a.txt

            cat ${sample_id}.trinotate_annotation_report.xls | cut -f 9 | grep [A-Z] | grep "|" | tr "\\`" "\\n" | \
                cut -f 1 -d "^" | cut -f 3 -d "|" | cut -f 2 -d "_" >b.txt

            cat a.txt b.txt | sort | uniq -c | sort -nr | head -n 20 | awk 'OFS="," {print \$1,\$2}' >${sample_id}_custom_uniprot_hits.txt

            rm a.txt b.txt

            cp ${params.pipeInstall}/conf/uni_tax.txt .
            cp ${sample_id}_custom_uniprot_hits.txt ${sample_id}_custom_uniprot_hits

            while read line;do
                a=\$( echo \${line} | cut -f 2 -d "," )
                b=\$( cat uni_tax.txt | grep "\${a}" | cut -f 2 -d "," | wc -l )
                if [ "\${b}" == "1" ];then
                    c=\$( cat uni_tax.txt | grep "\${a}" | cut -f 2 -d "," )
                    sed -i "s/\${a}/\${c}/" ${sample_id}_custom_uniprot_hits
                fi
            done <${sample_id}_custom_uniprot_hits.txt

            rm ${sample_id}_custom_uniprot_hits.txt uni_tax.txt
            mv ${sample_id}_custom_uniprot_hits ${sample_id}_custom_uniprot_hits.txt

            cp ${params.pipeInstall}/bin/custom_uniprot_hits.R .
            Rscript custom_uniprot_hits.R ${sample_id}

            cp ${sample_id}_custom_uniprot_hits.txt ${sample_id}_custom_uniprot_hits.csv
            """
    }

    if (!params.skipKegg) {
        process get_kegg {

            tag "${sample_id}"

            publishDir "${launchDir}/${params.outdir}/figures/kegg", mode: "copy", overwrite: true

            input:
                tuple sample_id, file(kegg) from kegg_paths

            output:
                tuple sample_id, file("${sample_id}_kegg.svg") into kegg_report

            script:
                """
                awk '{print \$2}' ${kegg} >kegg_terms
                curl -X POST --data-urlencode "selection@kegg_terms" -d "export_type=svg" -d "default_opacity=.5" -d "default_width=2" -d "default_radius=5" https://pathways.embl.de/mapping.cgi >${sample_id}_kegg.svg
                """
        }
    } else {
        process skip_kegg {
            tag "${sample_id}"

            input:
                tuple sample_id, file(kegg) from kegg_paths

            output:
                tuple sample_id, file("${sample_id}_kegg.svg") into kegg_report

            script:
                """
                touch ${sample_id}_kegg.svg
                """
        }
    }

    process get_transcript_dist {

        tag "${sample_id}"

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge -c bioconda biopython=1.78 pandas=1.1.2 numpy=1.18.1" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/mulled-v2-1e9d4f78feac0eb2c8d8246367973b3f6358defc:41ffac721ff9b03ca1121742e969d0e7d78e589f-0" : "quay.io/biocontainers/mulled-v2-1e9d4f78feac0eb2c8d8246367973b3f6358defc:41ffac721ff9b03ca1121742e969d0e7d78e589f-0")
        }

        input:
            tuple sample_id, file(dist) from evi_dist

        output:
            tuple sample_id, file("${sample_id}_sizes.txt") into size_dist

        script:
            """
            len.py ${dist} >all_transcript_sizes.txt
            bash get_sizes.sh all_transcript_sizes.txt
            mv final_sizes.txt ${sample_id}_sizes.txt
            """
    }

    if (!params.skipReport) {

        report_ch = Channel.create()
        fastp_csv.mix( norm_report, remove_rrna_sum, mapping_evi_results, mapping_trinity_results, rna_quast_report, size_dist, summary_evi_csv, busco3_csv, busco4_csv, busco3_heatmap, busco4_heatmap, transdecoder_csv, go_csv, uniprot_csv, kegg_report ).groupTuple(by:0,size:16).flatten().toList().into(report_ch)

        process get_report {

            label 'low_cpus'

            publishDir "${launchDir}/${params.outdir}/report", mode: "copy", overwrite: true, pattern: "*.{html,pdf}"

            conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge -c bioconda r-reshape2=1.4.4 r-plotly=4.9.2.1 plotly-orca=3.4.2 r-ggplot2=3.3.0 r-svglite=1.2.3 r-ggthemes=4.2.0 r-knitr=1.29 r-rmarkdown=2.3 r-kableextra=1.1.0" : null)
            if (params.oneContainer){ container "${params.TPcontainer}" } else {
            container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/mulled-v2-3f431f5f8e54df68ea0029c209fce3b154f6e186:94cad00b5306639ceab6aaf211f45740560abb90-0" : "quay.io/biocontainers/mulled-v2-3f431f5f8e54df68ea0029c209fce3b154f6e186:94cad00b5306639ceab6aaf211f45740560abb90-0")
            }

            input:
                file(files) from report_ch
                    .collect()

            output:
                file("*html") into final_report

            script:
                """
                sample_id=\$( cat input.1 )
                cp ${params.pipeInstall}/bin/TransPi_Report_Ind.Rmd .
                Rscript -e "rmarkdown::render('TransPi_Report_Ind.Rmd',output_file='TransPi_Report_\${sample_id}.html')" \${sample_id} ${params.skipFilter} ${params.skipNormalization} ${params.rRNAfilter} ${params.buscoDist} ${params.allBuscos} ${params.skipKegg}
                """
        }
    }

} else if (params.onlyEvi) {

    println("\n\tRunning Evidential Gene analysis only \n")

    Channel
        .fromFilePairs("${launchDir}/onlyEvi/*.{fa,fasta}", size: -1, checkIfExists: true)
        .set{ evigene_ch_OE }

    process evigene_OE {

        label 'med_mem'

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/evigene", mode: "copy", overwrite: true, pattern: "*.combined.okay.fa"

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::cd-hit=4.8.1 bioconda::exonerate=2.4 bioconda::blast=2.2.31" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/mulled-v2-962eae98c9ff8d5b31e1df7e41a355a99e1152c4:5aac6d1d2253d47aee81f01cc070a17664c86f07-0" : "quay.io/biocontainers/mulled-v2-962eae98c9ff8d5b31e1df7e41a355a99e1152c4:5aac6d1d2253d47aee81f01cc070a17664c86f07-0")
        }

        input:
            tuple sample_id, file(assembly) from evigene_ch_OE

        output:
            tuple sample_id, file("${sample_id}.combined.okay.fa") into ( evigene_ch_busco3_OE, evigene_ch_busco4_OE )

        script:
            def mem_MB=(task.memory.toMega())
            """
            echo -e "\\n-- Starting EviGene --\\n"

            cat ${assembly} >${sample_id}.combined.fa

            $evi/scripts/prot/tr2aacds.pl -tidy -NCPU ${task.cpus} -MAXMEM ${mem_MB} -log -cdna ${sample_id}.combined.fa

            echo -e "\\n-- DONE with EviGene --\\n"

            cp okayset/*combined.okay*.fa ${sample_id}.combined.okay.fa

            if [ -d tmpfiles/ ];then
                rm -rf tmpfiles/
            fi
            """
    }

    process busco3_OE {

        label 'med_cpus'

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/busco3", mode: "copy", overwrite: true

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::busco=3.0.2=py_13" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/busco:3.0.2--py_13" : "quay.io/biocontainers/busco:3.0.2--py_13")
        }

        input:
            tuple sample_id, file("${sample_id}.combined.okay.fa") from evigene_ch_busco3_OE

        output:
            tuple sample_id, file("*.TransPi.bus3") into busco3_ch_OE
            tuple sample_id, file("*.bus3.txt") into busco3_summary_OE

        script:
            """
            echo -e "\\n-- Starting BUSCO --\\n"

            run_BUSCO.py -i ${sample_id}.combined.okay.fa -o ${sample_id}.TransPi.bus3 -l ${params.busco3db} -m tran -c ${task.cpus}

            echo -e "\\n-- DONE with BUSCO --\\n"

            cp run_${sample_id}.TransPi.bus3/short_summary_${sample_id}.TransPi.bus3.txt .
            """
    }

    process busco4_OE {

        label 'med_cpus'

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/busco4", mode: "copy", overwrite: true

        // change container in oneContainer option
        conda (params.condaActivate && params.myConda ? params.cenv : params.condaActivate ? "-c conda-forge bioconda::busco=4.1.4=py_0" : null)
        if (params.oneContainer){ container "${params.v4container}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/busco:4.0.5--pyr36_0" : "ezlabgva/busco:v4.0.5_cv1")
        }

        input:
            tuple sample_id, file("${sample_id}.combined.okay.fa") from evigene_ch_busco4_OE

        output:
            tuple sample_id, file("*.TransPi.bus4") into busco4_ch_OE
            tuple sample_id, file("*.bus4.txt") into busco4_summary_OE

        script:
            """
            echo -e "\\n-- Starting BUSCO --\\n"

            busco -i ${sample_id}.combined.okay.fa -o ${sample_id}.TransPi.bus4 -l ${params.busco4db} -m tran -c ${task.cpus} --offline

            echo -e "\\n-- DONE with BUSCO --\\n"

            cp ${sample_id}.TransPi.bus4/short_summary.*.${sample_id}.TransPi.bus4.txt .
            """
    }

    busco4_comp_OE = Channel.create()
    busco4_summary_OE.mix(busco3_summary_OE).groupTuple(by:0,size:2).into(busco4_comp_OE)

    process summary_busco_OE {

        tag "${sample_id}"

        publishDir "${launchDir}/${params.outdir}/stats", mode: "copy", overwrite: true

        input:
            tuple sample_id, file(files) from busco4_comp_OE

        output:
            file("*.sum_busco.txt") into final_sum_2v3_OE

        script:
            """
            trans4=\$( echo $files | tr " " "\\n" | grep "short_summary\\." )
            trans3=\$( echo $files | tr " " "\\n" | grep "short_summary_" )
            echo -e "Summary of BUSCO V4 \\n" >>${sample_id}.sum_busco.txt
            echo "-- TransPi BUSCO V4 scores -- " >>${sample_id}.sum_busco.txt
            cat \${trans4} >>${sample_id}.sum_busco.txt
            echo -e "\\nSummary of BUSCO V3 \\n" >>${sample_id}.sum_busco.txt
            echo "-- TransPi BUSCO V3 scores -- " >>${sample_id}.sum_busco.txt
            cat \${trans3} >>${sample_id}.sum_busco.txt
            """
    }

} else {
    println("\n\t\033[0;31mMandatory argument not specified.\n\tFor more info use `nextflow run TransPi.nf --help`\n\033[0m")
    exit 0
}

if (params.getVersions) {
    process get_run_info {

        publishDir "${launchDir}/${params.outdir}/", mode: "copy", overwrite: true

        conda (params.condaActivate && params.myConda ? params.localConda : params.condaActivate ? "-c conda-forge bioconda::NAME-HERE" : null)
        if (params.oneContainer){ container "${params.TPcontainer}" } else {
        container (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container ? "https://depot.galaxyproject.org/singularity/NAME-HERE" : "quay.io/biocontainers/NAME-HERE")
        }

        output:
           file("versions.txt") into run_info

        script:
            """
            echo "==================================================" >>versions.txt
            echo "  TransPi - Transcriptome Analysis Pipeline v${workflow.manifest.version}" >>versions.txt
            echo -e "==================================================\\n" >>versions.txt
            echo -e "\\t\\tRUN INFO\\n" >>versions.txt
            echo "-- Kmers used --" >>versions.txt
            echo ${params.k} >>versions.txt

            echo -e "\\n-- Databases name and last update --" >>versions.txt

            v=\$( echo ${params.uniname} )
            echo "Uniprot_DB: \$v" >>versions.txt

            if [ -f ${params.pipeInstall}/DBs/uniprot_db/.lastrun.txt ];then
                v=\$( cat ${params.pipeInstall}/DBs/uniprot_db/.lastrun.txt )
            else
                v="No info available. Check Instructions on README."
            fi
            echo -e "Uniprot_DB last update: \$v \\n" >>versions.txt

            if [ -f ${params.pipeInstall}/DBs/hmmerdb/.lastrun.txt ];then
                v=\$( cat ${params.pipeInstall}/DBs/hmmerdb/.lastrun.txt )
            else
                v="No info available. Check Instructions on README."
            fi
            echo -e "PfamA last update: \$v \\n" >>versions.txt

            v=\$( echo ${params.busco3db} | tr "/" "\\n" | tail -n 1 )
            echo "BUSCO_v3_DB: \$v" >>versions.txt

            v=\$( echo ${params.busco4db} | tr "/" "\\n" | tail -n 1 )
            echo "BUSCO_v4_DB: \$v" >>versions.txt

            echo -e "\\n-- Program versions --" >>versions.txt

            v=\$( SOAPdenovo-Trans-127mer --version | grep "version" | awk '{print \$2,\$3}' | cut -f 1 -d ":" | cut -f 2 -d " " )
            echo "SOAP: \$v" >>versions.txt

            v=\$( velveth | grep "Version" | cut -f 2 -d " " )
            echo "Velveth: \$v" >>versions.txt

            v=\$( velvetg | grep "Version" | cut -f 2 -d " " )
            echo "Velvetg: \$v" >>versions.txt

            v=\$( oases | grep "Version" | cut -f 2 -d " " )
            echo "Oases: \$v" >>versions.txt

            v=\$( rnaspades.py -v )
            echo "rna-SPADES: \$v" >>versions.txt

            v=\$( transabyss --version )
            echo "Trans-ABySS: \$v" >>versions.txt

            v=\$( Trinity --version | grep "version" | head -n 1 | cut -f 2 -d "-" )
            echo "Trinity: \$v" >>versions.txt

            v=\$( diamond --version 2>&1 | tail -n 1 | cut -f 3 -d " " )
            echo "Diamond: \$v" >>versions.txt

            v=\$( hmmsearch -h | head -n 2 | cut -f 3 -d " " | grep [0-9] )
            echo "HMMER: \$v" >>versions.txt

            v=\$( echo "2019.05.14" )
            echo "EvidentialGene: \$v" >>versions.txt

            v=\$( TransDecoder.LongOrfs --version | cut -f 2 -d " " )
            echo "Transdecoder: \$v" >>versions.txt

            v=\$( run_BUSCO.py -v | cut -f 2 -d " " )
            echo "BUSCO3: \$v" >>versions.txt

            v=\$( echo "4.0.5" )
            echo "BUSCO4: \$v" >>versions.txt

            v=\$( cat ${params.pipeInstall}/transpi_env.yml | grep trinotate  | cut -f 2 -d "=" )
            echo "Trinotate: \$v" >>versions.txt

            v=\$( cd-hit -h | head -n1 | cut -f 1 -d "(" | cut -f 2 -d "n" )
            echo "CD-HIT:\$v" >>versions.txt

            v=\$( exonerate -v | head -n1 | cut -f 5 -d " " )
            echo "Exonerate: \$v" >>versions.txt

            echo -e "\\n-- Programming Languages --" >>versions.txt

            v=\$( R --version | grep "R version" | awk '{print \$3}' )
            echo "R: \$v" >>versions.txt

            v=\$( python --version | cut -f 2 -d " " )
            echo "Python: \$v" >>versions.txt

            v=\$( perl -v | head -n2 | grep version | cut -f 1 -d ")" | cut -f 2 -d "(" | tr -d "v" )
            echo "Perl: \$v" >>versions.txt
            """
    }
}
workflow.onComplete {
    log.info ( workflow.success ? \
        "---------------------------------------------------------------------------------" \
        + "\n\033[0;32mDone! Open the following reports in your browser\033[0m" \
        + "\n\033[0;32mPipeline performance report: ${launchDir}/${params.outdir}/${params.tracedir}/transpi_report.html\033[0m" \
        + "\n\033[0;32mTransPi (--all) interactive report: ${launchDir}/${params.outdir}/report/TransPi_Report_*.hmtl\033[0m" \
        : \
        "---------------------------------------------------------------------------------" \
        + "\n\033[0;31mSomething went wrong. Check error message below and/or log files.\033[0m" )
}
