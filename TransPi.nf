#!/usr/bin/env nextflow

/*
========================================================================================
                                    TransPi
========================================================================================
                       Transcriptomes Analysis Pipeline
                       Author: Ramon E. Rivera-Vicens
                       Version: 1.0 (dev)
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info """
    ==========================================
    TransPi - Transcriptome Analysis Pipeline
    ==========================================

        Steps:
            1- Run the `precheck_TransPi.sh` to install tools, set up the databases and directories for TransPi
            2- Run TransPi

            Usage:

                nextflow run TransPi.nf --all (other_options_here)

        Mandatory arguments (--all or --onlyAsm or --onlyEvi or --onlyAnn):

                --all           Run the entire pipeline (Assemblies, EvidentialGene, Annotation, etc.)

                --onlyAsm       Run only the Assemblies and EvidentialGene analysis (testing)

                --onlyEvi       Run only the Evidential Gene analysis (testing)

                --onlyAnn       Run only the Annotation analysis (starting from a final assembly)

        Other options:

                -with-conda     To run with a local conda installation (generated with the precheck) and not installed by nextflow
                                This is the preferred method of running TransPi

                                Example:
                                    nextflow run transpi.nf --all -with-conda /home/ubuntu/anaconda3/envs/TransPi

                -profile        Configuration profile to use. Can use multiple (comma separated)

                                Available:
                                    test (ready - Run TransPi with test dataset)
                                    conda (ready - Not recommended, not all programs are installed by conda. Use the precheck)
                                    docker (in development - Run TransPi with a docker container with all the neccesary tools)
                                    singularity (in development - Run TransPi with a singularity container with all the neccesary tools)

                --help          Display this message

                --fullHelp      Display more options and examples for running TransPi

    """.stripIndent()
}
def fullHelpMessage() {
    log.info """
    ==========================================
    TransPi - Transcriptome Analysis Pipeline
    ==========================================

        Steps:
            1- Run the `precheck_TransPi.sh` to install tools, set up the databases and directories for TransPi
            2- Run TransPi

            Usage:

                nextflow run TransPi.nf --all (other_options_here)

        Mandatory arguments (--all or --onlyAsm or --onlyEvi or --onlyAnn):

                --all           Run the entire pipeline (Assemblies, EvidentialGene, Annotation, etc.)

                --onlyAsm       Run only the Assemblies and EvidentialGene analysis (testing)

                --onlyEvi       Run only the Evidential Gene analysis (testing)

                --onlyAnn       Run only the Annotation analysis (starting from a final assembly)

        Other options:

                -with-conda     To run with a local conda installation (generated with the precheck) and not installed by nextflow
                                This is the preferred method of running TransPi

                                Example:
                                    nextflow run transpi.nf --all -with-conda /home/ubuntu/anaconda3/envs/TransPi

                -profile        Configuration profile to use. Can use multiple (comma separated)

                                Available:
                                    test (ready - Run TransPi with test dataset)
                                    conda (ready - Not recommended, not all programs are installed by conda. Use the precheck)
                                    docker (in development - Run TransPi with a docker container with all the neccesary tools)
                                    singularity (in development - Run TransPi with a singularity container with all the neccesary tools)

        #################################################################################################

                                    Various examples on how to deploy TransPi

        #################################################################################################

        I. Steps for running on a local cluster and local installation of TransPi

            1- Run the `precheck_TransPi.sh` to install tools, set up the databases and directories for TransPi
            2- Run TransPi

            Usage:

                nextflow run TransPi.nf --all -with-conda ~/anaconda3/envs/TransPi

        #################################################################################################

        II. Steps for running on a local cluster and conda installation by nextflow

            1- Run the `precheck_TransPi.sh` to install tools, set up the databases and directories for TransPi
            2- Run TransPi

            Usage:

                nextflow run TransPi.nf --all -profile conda

            NOTE:
                Not recommended, not all programs are installed by conda. Use if other dependencies are manually installed

        #################################################################################################

        III. Steps for running on docker (in development)

            1- Run the `container_precheck_TransPi.sh` to install tools, set up the databases and directories for TransPi
            2- Run TransPi

            Usage:

                nextflow run TransPi.nf --all -profile docker

            NOTE:
                All necesary tools for running TransPi are pre-installed in the container
                `container_precheck_TransPi.sh` will install only the database used by TransPi

        #################################################################################################

        IV. Steps for running on singualarity (in development)

            1- Run the `container_precheck_TransPi.sh` to install tools, set up the databases and directories for TransPi
            2- Run TransPi

            Usage:

                nextflow run TransPi.nf --all -profile singularity

            NOTE:
                All necesary tools for running TransPi are pre-installed in the container
                `container_precheck_TransPi.sh` will install only the database used by TransPi

        #################################################################################################

        V. Steps for running with multiple profiles (in development - depending on profile selected)

            1- Run the `precheck_TransPi.sh` to install tools, set up the databases and directories for TransPi
            2- Run TransPi

            Usage:

                nextflow run TransPi.nf --all -profile conda,test

            NOTE:
                This will run TransPi using a test dataset and a conda environment created by Nextflow
                To run with containers first run the `container_precheck_TransPi.sh` and use -profile docker,test

        #################################################################################################

    """.stripIndent()
}

def workDir = System.getProperty("user.dir");

// Show help message
if (params.help) {
    helpMessage()
    exit 0
} else if (params.fullHelp) {
    fullHelpMessage()
    exit 0
}

log.info """\
        =========================================
        TransPi - Transcriptome Analysis Pipeline
        =========================================
        TransPi Installation:       ${params.mypwd}
        Reads Length:               ${params.max_rd_len}
        Kmers:                      ${params.k}
        Working Directory:          ${workDir}
        Uniprot DB:                 ${params.uniprot}
        Busco DB:                   ${params.busco4db}
        """.stripIndent()

if (params.readsTest) {
    println("\n\tRunning TransPi with TEST dataset\n")
    Channel
        .from(params.readsTest)
        .map{ row -> [ row[0], [ file(row[1][0],checkIfExists: true),file(row[2][0],checkIfExists: true) ] ] }
        .ifEmpty{ exit 1, "params.readsTest was empty - no input files supplied" }
        .into{ reads_ch; reads_qc_ch }
} else {
    println("\n\tRunning TransPi with your dataset\n")
    Channel
        .fromFilePairs("${params.reads}", checkIfExists: true)
        .into{ reads_ch; reads_qc_ch }
}

if (params.onlyAsm) {

    //testing
    println("\n\tRunning assemblies and Evidential Gene analysis only \n")

    if (!params.skipQC) {

        process fasqc_OAS {

            label 'low_cpus'

            tag "${sample_id}"

            publishDir "${workDir}/${params.outdir}/fastqc", mode: "copy", overwrite: true

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

            label 'med_cpus' // check this later

            tag "${sample_id}"

            publishDir "${workDir}/${params.outdir}/filter", mode: "copy", overwrite: true, pattern: "*.fastp.{json,html}"

            input:
                tuple sample_id, file(reads) from reads_ch

            output:
                tuple sample_id, file("*.fastp.{json,html}") into fastp_results_OAS
                tuple sample_id, file("*.filter.fq") into reads_ass_ch_OAS

            script:
                """
                echo ${sample_id}

                fastp -i ${reads[0]} -I ${reads[1]} -o left-${sample_id}.filter.fq -O right-${sample_id}.filter.fq --detect_adapter_for_pe \
                --average_qual ${params.minQual} --overrepresentation_analysis --html ${sample_id}.fastp.html --json ${sample_id}.fastp.json \
                --thread ${task.cpus} --report_title ${sample_id}
                """
        }
    } else {
        reads_ch
            .set{ reads_ass_ch_OAS }
        fastp_results_OAS = Channel.empty()
    }

    process normalize_reads_OAS {

        label 'med_mem'

        tag "${sample_id}"

        input:
            tuple sample_id, file(reads) from reads_ass_ch_OAS

        output:
            tuple sample_id, file("left-${sample_id}.norm.fq"), file("right-${sample_id}.norm.fq") into ( norm_reads_soap_OAS, norm_reads_velvet_OAS, norm_reads_trinity_OAS, norm_reads_spades_OAS, norm_reads_transabyss_OAS )

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

            mv left.norm.fq left-"${sample_id}".norm.fq
            mv right.norm.fq right-"${sample_id}".norm.fq
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

            mv left.norm.fq left-"${sample_id}".norm.fq
            mv right.norm.fq right-"${sample_id}".norm.fq

            rm left-${sample_id}.fq right-${sample_id}.fq
            """
        }
    }

    process trinity_assembly_OAS {

        label 'med_mem'

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/assemblies", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("left-${sample_id}.norm.fq"), file("right-${sample_id}.norm.fq") from norm_reads_trinity_OAS

        output:
            tuple sample_id, file("${sample_id}.Trinity.fa") into ( assemblies_ch_trinity_OAS, busco3_ch_trinity_OAS, busco4_ch_trinity_OAS )

        script:
            """
            mem=\$( echo ${task.memory} | cut -f 1 -d " " )

            Trinity --max_memory \${mem}G --seqType fq --left left-${sample_id}.norm.fq --right right-${sample_id}.norm.fq --CPU ${task.cpus} --no_normalize_reads --full_cleanup --output trinity_out_dir

            mv trinity_out_dir.Trinity.fasta ${sample_id}.Trinity.fa
            """
    }

    process soap_assembly_OAS {

        label 'med_mem'

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/assemblies", mode: "copy", overwrite: true

        input:
            val k from "${params.k}"
            tuple sample_id, file("left-${sample_id}.norm.fq"), file("right-${sample_id}.norm.fq") from norm_reads_soap_OAS

        output:
            tuple sample_id, file("${sample_id}.SOAP.fa") into assemblies_ch_soap_OAS

        script:
            """
            echo -e "\\n-- Generating SOAP config file --\\n"
            echo "max_rd_len="${params.max_rd_len} >>config.txt
            echo "[LIB]" >>config.txt
            echo "rd_len_cutof="${params.rd_len_cutof} >>config.txt
            #echo "avg_ins="${params.avg_ins} >>config.txt
            echo "reverse_seq="${params.reverse_seq} >>config.txt
            echo "asm_flags="${params.asm_flags} >>config.txt
            echo "map_len="${params.map_len} >>config.txt
            echo "q1="left-${sample_id}.norm.fq >>config.txt
            echo "q2="right-${sample_id}.norm.fq >>config.txt

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

        publishDir "${workDir}/${params.outdir}/assemblies", mode: "copy", overwrite: true

        input:
            val k from "${params.k}"
            tuple sample_id, file("left-${sample_id}.norm.fq"), file("right-${sample_id}.norm.fq") from norm_reads_velvet_OAS

        output:
            tuple sample_id, file("${sample_id}.Velvet.fa") into assemblies_ch_velvet_OAS

        script:
            """
            echo -e "\\n-- Starting with Velveth --\\n"
            for x in `echo $k | tr "," " "`;do
                echo -e "\\n-- k\${x} --\\n"
                velveth oases.\${x} \${x} -shortPaired -fastq -separate left-${sample_id}.norm.fq right-${sample_id}.norm.fq
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

        publishDir "${workDir}/${params.outdir}/assemblies", mode: "copy", overwrite: true

        input:
            val k from "${params.k}"
            tuple sample_id, file("left-${sample_id}.norm.fq"), file("right-${sample_id}.norm.fq") from norm_reads_spades_OAS

        output:
            tuple sample_id, file("${sample_id}.SPADES.fa") into assemblies_ch_spades_OAS

        script:
            """
            echo -e "\\n-- Starting rnaSPADES assemblies --\\n"

            mem=\$( echo ${task.memory} | cut -f 1 -d " " )

            for x in `echo $k | tr "," " "`;do
                echo -e "\\n-- rnaSPADES k\${x} --\\n"

                rnaspades.py -1 left-${sample_id}.norm.fq -2 right-${sample_id}.norm.fq -o ${sample_id}_spades_\${x} -t ${task.cpus} -k \${x} -m \${mem}

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

        publishDir "${workDir}/${params.outdir}/assemblies", mode: "copy", overwrite: true

        input:
            val k from "${params.k}"
            tuple sample_id, file("left-${sample_id}.norm.fq"), file("right-${sample_id}.norm.fq") from norm_reads_transabyss_OAS

        output:
            tuple sample_id, file("${sample_id}.TransABySS.fa") into assemblies_ch_transabyss_OAS

        script:
            """
            echo -e "\\n-- Starting Trans-ABySS assemblies --\\n"

            for x in `echo $k | tr "," " "`;do
                echo -e "\\n-- Trans-ABySS k\${x} --\\n"

                transabyss -k \${x} --pe left-${sample_id}.norm.fq right-${sample_id}.norm.fq --outdir ${sample_id}_transabyss_\${x} --name k\${x}.transabyss.fa --threads ${task.cpus} -c 12 --length 200

            done

            echo -e "\\n-- Finished with the assemblies --\\n"

            for x in `echo $k | tr "," " "`;do
                sed -i "s/>/>TransABySS.k\${x}./g" ${sample_id}_transabyss_\${x}/k\${x}.transabyss.fa-final.fa
            done

            cat ${sample_id}_transabyss_*/k*.transabyss.fa-final.fa >${sample_id}.TransABySS.fa

            rm -rf ${sample_id}_transabyss_*
            """
    }

    process evigene_OAS {

        label 'med_mem'

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/evigene", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.Trinity.fa") from assemblies_ch_trinity_OAS
            tuple sample_id, file("${sample_id}.SOAP.fa") from assemblies_ch_soap_OAS
            tuple sample_id, file("${sample_id}.Velvet.fa") from assemblies_ch_velvet_OAS
            tuple sample_id, file("${sample_id}.SPADES.fa") from assemblies_ch_spades_OAS
            tuple sample_id, file("${sample_id}.TransABySS.fa") from assemblies_ch_transabyss_OAS

        output:
            tuple sample_id, file("*.combined.okay.fa") into ( evigene_ch_busco3_OAS, evigene_ch_busco4_OAS )
            tuple sample_id, file("${sample_id}.combined.fa"), file("${sample_id}.combined.okay.fa") into evigene_summary_OAS

        script:
            def mem_MB=(task.memory.toMega())

            """
            echo -e "\\n-- Starting EviGene --\\n"

            cat *.fa >${sample_id}.combined.fa

            $evi/scripts/prot/tr2aacds.pl -tidy -NCPU ${task.cpus} -MAXMEM ${mem_MB} -log -cdna ${sample_id}.combined.fa

            echo -e "\\n-- DONE with EviGene --\\n"

            cp okayset/${sample_id}.combined.okay.combined.fa ${sample_id}.combined.okay.fa

            if [ -d tmpfiles/ ];then
                rm -rf tmpfiles/
            fi
            """
    }

    process summary_evigene_individual_OAS {

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/stats", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.combined.fa"), file("${sample_id}.combined.okay.fa") from evigene_summary_OAS

        output:
            tuple sample_id, file("${sample_id}.sum_preEG.txt"), file("${sample_id}.sum_EG.txt") into final_sum_1_OAS

        script:
            """
            #Summary of total number of transcripts
            echo -e "- Number of transcripts before Evidential Genes\\n" >>${sample_id}.sum_preEG.txt
            echo -e "- Individual ${sample_id} \\n" >>${sample_id}.sum_preEG.txt
            echo -e "\\t Total transcripts:" >>${sample_id}.sum_preEG.txt
            num=\$( cat ${sample_id}.combined.fa | grep -c ">" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}.sum_preEG.txt
            echo -e "\\t Trinity" >>${sample_id}.sum_preEG.txt
            num=\$( cat ${sample_id}.combined.fa | grep -c ">TRINITY" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}.sum_preEG.txt
            echo -e "\\t SOAP" >>${sample_id}.sum_preEG.txt
            num=\$( cat ${sample_id}.combined.fa | grep -c ">SOAP" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}.sum_preEG.txt
            echo -e "\\t Velvet/Oases" >>${sample_id}.sum_preEG.txt
            num=\$( cat ${sample_id}.combined.fa | grep -c ">Velvet" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}.sum_preEG.txt
            echo -e "\\t rna-SPADES" >>${sample_id}.sum_preEG.txt
            num=\$( cat ${sample_id}.combined.fa | grep -c ">SPADES" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}.sum_preEG.txt
            echo -e "\\t Trans-ABySS" >>${sample_id}.sum_preEG.txt
            num=\$( cat ${sample_id}.combined.fa | grep -c ">TransABySS" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}.sum_preEG.txt

            #Summary of transcripts after EvidentialGenes
            echo -e "- Number of transcripts by individual after EvidentialGenes\\n" >>${sample_id}.sum_EG.txt
            echo -e "- Individual ${sample_id} \\n" >>${sample_id}.sum_EG.txt
            echo -e "\\t Total transcripts:" >>${sample_id}.sum_EG.txt
            num=\$( cat ${sample_id}.combined.okay.fa | grep -c ">" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}.sum_EG.txt
            echo -e "\\t Trinity" >>${sample_id}.sum_EG.txt
            num=\$( cat ${sample_id}.combined.okay.fa | grep -c ">TRINITY" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}.sum_EG.txt
            echo -e "\\t SOAP" >>${sample_id}.sum_EG.txt
            num=\$( cat ${sample_id}.combined.okay.fa | grep -c ">SOAP" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}.sum_EG.txt
            echo -e "\\t Velvet/Oases" >>${sample_id}.sum_EG.txt
            num=\$( cat ${sample_id}.combined.okay.fa | grep -c ">Velvet" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}.sum_EG.txt
            echo -e "\\t rna-SPADES" >>${sample_id}.sum_EG.txt
            num=\$( cat ${sample_id}.combined.okay.fa | grep -c ">SPADES" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}.sum_EG.txt
            echo -e "\\t Trans-ABySS" >>${sample_id}.sum_EG.txt
            num=\$( cat ${sample_id}.combined.fa | grep -c ">TransABySS" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}.sum_EG.txt
            """
    }

    process busco3_OAS {

        label 'med_cpus'

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/busco3", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.combined.okay.fa") from evigene_ch_busco3_OAS

        output:
            tuple sample_id, file("run_${sample_id}.fa.bus") into busco3_ch
            tuple sample_id, file("*${sample_id}.fa.bus.txt") into ( busco3_summary_OAS, busco3_comp_1_OAS )

        script:
            """
            echo -e "\\n-- Starting BUSCO --\\n"

            run_BUSCO.py -i ${sample_id}.combined.okay.fa -o ${sample_id}.fa.bus -l ${params.busco3db} -m tran -c ${task.cpus}

            echo -e "\\n-- DONE with BUSCO --\\n"

            cp run_${sample_id}.fa.bus/short_summary_${sample_id}.fa.bus.txt .
            """
    }

    process busco3_tri_OAS {

        label 'med_cpus'

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/busco3", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.Trinity.fa") from busco3_ch_trinity_OAS

        output:
            tuple sample_id, file("*${sample_id}.Trinity.fa.bus.txt") into ( busco3_ch_trinity_sum_OAS, busco3_comp_2_OAS )

        script:
            """
            echo -e "\\n-- Starting BUSCO --\\n"

            run_BUSCO.py -i ${sample_id}.Trinity.fa -o ${sample_id}.Trinity.fa.bus -l ${params.busco3db} -m tran -c ${task.cpus}

            echo -e "\\n-- DONE with BUSCO --\\n"

            cp run_${sample_id}.Trinity.fa.bus/short_summary_${sample_id}.Trinity.fa.bus.txt .
            """
    }

    process busco4_OAS {

        conda "${params.cenv}"

        label 'med_cpus'

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/busco4", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.combined.okay.fa") from evigene_ch_busco4_OAS

        output:
            tuple sample_id, file("${sample_id}.fa.bus") into busco4_ch
            tuple sample_id, file("*${sample_id}.fa.bus.txt") into ( busco4_summary_OAS, busco4_comp_1_OAS )

        script:
            """
            echo -e "\\n-- Starting BUSCO --\\n"

            busco -i ${sample_id}.combined.okay.fa -o ${sample_id}.fa.bus -l ${params.busco4db} -m tran -c ${task.cpus} --offline

            echo -e "\\n-- DONE with BUSCO --\\n"

            cp ${sample_id}.fa.bus/short_summary.*.${sample_id}.fa.bus.txt .
            """
    }

    process busco4_tri_OAS {

        conda "${params.cenv}"

        label 'med_cpus'

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/busco4", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.Trinity.fa") from busco4_ch_trinity_OAS

        output:
            tuple sample_id, file("*${sample_id}.Trinity.fa.bus.txt") into ( busco4_ch_trinity_sum_OAS, busco4_comp_2_OAS )

        script:
            """
            echo -e "\\n-- Starting BUSCO --\\n"

            busco -i ${sample_id}.Trinity.fa -o ${sample_id}.Trinity.fa.bus -l ${params.busco4db} -m tran -c ${task.cpus} --offline

            echo -e "\\n-- DONE with BUSCO --\\n"

            cp ${sample_id}.Trinity.fa.bus/short_summary.*.${sample_id}.Trinity.fa.bus.txt .
            """
    }

    process summary_busco3_individual_OAS {

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/stats", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.fa.bus.txt") from busco3_summary_OAS
            tuple sample_id, file("${sample_id}.Trinity.fa.bus.txt") from busco3_ch_trinity_sum_OAS

        output:
            tuple sample_id, file("${sample_id}.sum_busco3.txt") into final_sum_2v3_OAS

        script:
            """
            #Summary of BUSCO scores for the final_assemblies
            echo -e "Summary of BUSCO V3 \n" >>${sample_id}.sum_busco3.txt
            echo "-- TransPi BUSCO V3 scores -- " >>${sample_id}.sum_busco3.txt
            cat ${sample_id}.fa.bus.txt >>${sample_id}.sum_busco3.txt
            echo -e "\\n-- Trinity BUSCO V3 scores --" >>${sample_id}.sum_busco3.txt
            cat ${sample_id}.Trinity.fa.bus.txt >>${sample_id}.sum_busco3.txt
            """
    }

    process summary_busco4_individual_OAS {

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/stats", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.fa.bus.txt") from busco4_summary_OAS
            tuple sample_id, file("${sample_id}.Trinity.fa.bus.txt") from busco4_ch_trinity_sum_OAS

        output:
            tuple sample_id, file("${sample_id}.sum_busco4.txt") into final_sum_2v4_OAS

        script:
            """
            #Summary of BUSCO scores for the final_assemblies
            echo -e "Summary of BUSCO V4 \n" >>${sample_id}.sum_busco4.txt
            echo "-- TransPi BUSCO V4 scores -- " >>${sample_id}.sum_busco4.txt
            cat ${sample_id}.fa.bus.txt >>${sample_id}.sum_busco4.txt
            echo -e "\\n-- Trinity BUSCO V4 scores --" >>${sample_id}.sum_busco4.txt
            cat ${sample_id}.Trinity.fa.bus.txt >>${sample_id}.sum_busco4.txt
            """
    }

    process get_busco3_comparison_OAS {

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/figures/BUSCO3", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.fa.bus.txt") from busco3_comp_1_OAS
            tuple sample_id, file("${sample_id}.Trinity.fa.bus.txt") from busco3_comp_2_OAS

        output:
            tuple sample_id, file("${sample_id}_BUSCO3_comparison.pdf"), file("${sample_id}_BUSCO3_comparison.svg") into busco3_fig_OAS

        script:
            """
            set +e
            bash get_busco_val.sh ${sample_id}.Trinity.fa.bus.txt ${sample_id}.fa.bus.txt v3
            cp ${params.mypwd}/bin/busco_comparison.R .
            a=\$( cat final_spec )
            sed -i "s/MYSPEC/\${a}/" busco_comparison.R
            b=\$( cat final_perc )
            sed -i "s/MYPERC/\${b}/" busco_comparison.R
            c=\$( cat final_num )
            sed -i "s/MYVAL/\${c}/" busco_comparison.R
            Rscript busco_comparison.R ${sample_id}
            mv ${sample_id}_BUSCO_comparison.pdf ${sample_id}_BUSCO3_comparison.pdf
            mv ${sample_id}_BUSCO_comparison.svg ${sample_id}_BUSCO3_comparison.svg
            """
    }

    process get_busco4_comparison_OAS {

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/figures/BUSCO4", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.fa.bus.txt") from busco4_comp_1_OAS
            tuple sample_id, file("${sample_id}.Trinity.fa.bus.txt") from busco4_comp_2_OAS

        output:
            tuple sample_id, file("${sample_id}_BUSCO4_comparison.pdf"), file("${sample_id}_BUSCO4_comparison.svg") into busco4_fig_OAS

        script:
            """
            set +e
            bash get_busco_val.sh ${sample_id}.Trinity.fa.bus.txt ${sample_id}.fa.bus.txt v4
            cp ${params.mypwd}/bin/busco_comparison.R .
            a=\$( cat final_spec )
            sed -i "s/MYSPEC/\${a}/" busco_comparison.R
            b=\$( cat final_perc )
            sed -i "s/MYPERC/\${b}/" busco_comparison.R
            c=\$( cat final_num )
            sed -i "s/MYVAL/\${c}/" busco_comparison.R
            Rscript busco_comparison.R ${sample_id}
            mv ${sample_id}_BUSCO_comparison.pdf ${sample_id}_BUSCO4_comparison.pdf
            mv ${sample_id}_BUSCO_comparison.svg ${sample_id}_BUSCO4_comparison.svg
            """
    }

} else if (params.onlyAnn) {

    println("\n\tRunning only annotation analysis\n")

    Channel
        .fromFilePairs("${params.workDir}/onlyAnn/*.{fa,fasta}", size: -1, checkIfExists: true)
        .into{ annotation_ch_transdecoder_OA; assembly_ch_diamond_OA; assembly_ch_diamond_custom_OA; assembly_ch_rnammer_OA; assembly_ch_trinotate_OA }

    process custom_diamond_db_OA {
        script:
            """
            cd ${params.mypwd}
            echo -e "-- Checking if Diamond database folder is present --\\n"
            if [ ! -d DBs/diamonddb_custom/ ];then
                echo -e "-- Folder is not present, creating one and the Diamond database --\\n"
                mkdir -p DBs/diamonddb_custom/
                cd DBs/diamonddb_custom
                cp ${params.uniprot} .
                diamond makedb --in ${params.uniname} -d ${params.uniname}
                export unidb=`pwd`/${params.uniname}
                cd ../
            elif [ -d DBs/diamonddb_custom/ ];then
                echo -e "-- Folder is present. Checking if Diamond database is built --\\n"
                cd DBs/diamonddb_custom
                if [ ! -e ${params.uniname}.dmnd ];then
                    echo -e "-- Diamond database not present, creating one --\\n"
                    cp ${params.uniprot} .
                    diamond makedb --in ${params.uniname} -d ${params.uniname}
                    export unidb=`pwd`/${params.uniname}
                elif [ -e ${params.uniname}.dmnd  ];then
                    echo -e "-- Diamond database already created --\\n"
                    export unidb=`pwd`/${params.uniname}
                fi
                cd ../
            fi
            """
    }

    process hmmer_db_OA {
        script:
            """
            cd ${params.mypwd}
            echo -e "-- Checking if HMMER database folder is present --\\n"
            if [ -d DBs/hmmerdb/ ];then
                echo -e "-- Folder is present. Checking if HMMER database is built --\\n"
                cd DBs/hmmerdb
                if [ ! -e ${params.pfname}.h3f ] && [ ! -e ${params.pfname}.h3i ] && [ ! -e ${params.pfname}.h3m ] && [ ! -e ${params.pfname}.h3p ];then
                    echo -e "-- HMMER database not present, creating one --\\n"
                    hmmpress ${params.pfname}
                    export pf=`pwd`/${params.pfname}
                elif [ -s ${params.pfname}.h3f ] && [ -s ${params.pfname}.h3i ] && [ -s ${params.pfname}.h3m ] && [ -s ${params.pfname}.h3p ];then
                    echo -e "-- HMMER database already created --\\n"
                    export pf=`pwd`/${params.pfname}
                fi
                cd ../
            fi
            """
    }

    process swiss_diamond_db_OA {
        script:
            """
            cd ${params.mypwd}/DBs/sqlite_db
            if [ -e uniprot_sprot.pep ];then
                cd ${params.mypwd}
                if [ ! -d DBs/diamonddb_swiss/ ];then
                    echo -e "-- Folder is not present, creating one and the Diamond database --\\n"
                    mkdir -p DBs/diamonddb_swiss
                    cd DBs/diamonddb_swiss
                    cp ${params.mypwd}/DBs/sqlite_db/uniprot_sprot.pep .
                    diamond makedb --in uniprot_sprot.pep -d uniprot_sprot.pep
                    export swissdb=`pwd`/uniprot_sprot.pep
                elif [ -d DBs/diamonddb_swiss/ ];then
                    cd DBs/diamonddb_swiss
                    if [ ! -e uniprot_sprot.pep.dmnd ];then
                        echo -e "-- Diamond database not present, creating one --\\n"
                        cp ${params.mypwd}/DBs/sqlite_db/uniprot_sprot.pep .
                        diamond makedb --in uniprot_sprot.pep -d uniprot_sprot.pep
                        export swissdb=`pwd`/uniprot_sprot.pep
                    elif [ -e uniprot_sprot.pep.dmnd ];then
                        echo -e "-- Diamond database already created --\\n"
                        export swissdb=`pwd`/uniprot_sprot.pep
                    fi
                fi
            elif [ ! -e uniprot_sprot.pep ];then
                cd ${params.mypwd}
                if [ ! -d DBs/diamonddb_swiss/ ];then
                    echo -e "-- Folder is not present, creating one and the Diamond database --\\n"
                    mkdir -p DBs/diamonddb_swiss
                    cd DBs/diamonddb_swiss
                    wget http://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
                    EMBL_swissprot_parser.pl uniprot_sprot.dat.gz ind
                    rm ind.*
                    mv uniprot_sprot.dat.gz.pep uniprot_sprot.pep
                    diamond makedb --in uniprot_sprot.pep -d uniprot_sprot.pep
                    export swissdb=`pwd`/uniprot_sprot.pep
                elif [ -d DBs/diamonddb_swiss/ ];then
                    echo -e "-- Folder is present. Checking if Diamond database is built --\\n"
                    cd DBs/diamonddb_swiss
                    if [ ! -e uniprot_sprot.pep.dmnd ];then
                        if [ ! -e uniprot_sprot.pep ];then
                            echo -e "-- Diamond database not present, creating one --\\n"
                            wget http://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
                            EMBL_swissprot_parser.pl uniprot_sprot.dat.gz ind
                            rm ind.*
                            mv uniprot_sprot.dat.gz.pep uniprot_sprot.pep
                            diamond makedb --in uniprot_sprot.pep -d uniprot_sprot.pep
                            export swissdb=`pwd`/uniprot_sprot.pep
                        elif [ -e uniprot_sprot.pep ];then
                            diamond makedb --in uniprot_sprot.pep -d uniprot_sprot.pep
                            export swissdb=`pwd`/uniprot_sprot.pep
                        fi
                    elif [ -e uniprot_sprot.pep.dmnd ];then
                        echo -e "-- Diamond database already created --\\n"
                        export swissdb=`pwd`/uniprot_sprot.pep
                    fi
                fi
            fi
            """
    }

    process transdecoder_OA {

        label 'med_cpus'

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/transdecoder", mode: "copy", overwrite: true

        input:
            tuple sample_id, file(assembly) from annotation_ch_transdecoder_OA

        output:
            tuple sample_id, file("${sample_id}.*.transdecoder.pep") into ( transdecoder_ch_diamond_OA, transdecoder_ch_hmmer_OA, transdecoder_ch_signalp_OA, transdecoder_ch_tmhmm_OA, transdecoder_ch_trinotate_OA, transdecoder_ch_diamond_custom_OA )
            tuple sample_id, file("${sample_id}.transdecoder.stats") into transdecoder_summary_OA
            tuple sample_id, file("${sample_id}.*.transdecoder.{cds,gff,bed}") into transdecoder_files_OA

        script:
        if (params.shortTransdecoder) {
            """
            echo -e "\\n-- TransDecoder.LongOrfs... --\\n"

            TransDecoder.LongOrfs -t ${assembly}

            echo -e "\\n-- Done with TransDecoder.LongOrfs --\\n"

            echo -e "\\n-- TransDecoder.Predict... --\\n"

            TransDecoder.Predict -t ${assembly}

            echo -e "\\n-- Done with TransDecoder.Predict --\\n"

            echo -e "\\n-- Calculating statistics... --\\n"

            #Calculate statistics of Transdecoder
            echo "- Transdecoder (short,no homolgy) stats for ${sample_id}" >>${sample_id}.transdecoder.stats
            orfnum=\$( cat ${sample_id}.*.transdecoder.pep | grep -c ">" )
            echo -e "Total number of ORFs: \$orfnum \\n" >>${sample_id}.transdecoder.stats
            orfnum=\$( cat ${sample_id}.*.transdecoder.pep | grep -c "ORF type:complete" )
            echo -e "\\t ORFs type=complete: \$orfnum \\n" >>${sample_id}.transdecoder.stats
            orfnum=\$( cat ${sample_id}.*.transdecoder.pep | grep -c "ORF type:5prime_partial" )
            echo -e "\\t ORFs type=5prime_partial: \$orfnum \\n" >>${sample_id}.transdecoder.stats
            orfnum=\$( cat ${sample_id}.*.transdecoder.pep | grep -c "ORF type:3prime_partial" )
            echo -e "\\t ORFs type=3prime_partial: \$orfnum \\n" >>${sample_id}.transdecoder.stats
            orfnum=\$( cat ${sample_id}.*.transdecoder.pep | grep -c "ORF type:internal" )
            echo -e "\\t ORFs type=internal: \$orfnum \\n">>${sample_id}.transdecoder.stats

            echo -e "\\n-- Done with statistics --\\n"

            echo -e "\\n-- DONE with TransDecoder --\\n"
            """
        } else {
            """
            unidb=${params.mypwd}/DBs/diamonddb_custom/${params.uniname}

            echo -e "\\n-- TransDecoder.LongOrfs... --\\n"

            TransDecoder.LongOrfs -t ${assembly}

            echo -e "\\n-- Done with TransDecoder.LongOrfs --\\n"

            fname=${assembly}

            echo -e "\\n-- Starting Diamond (blastp) --\\n"

            diamond blastp -d \$unidb -q \$fname.transdecoder_dir/longest_orfs.pep -p ${task.cpus} -f 6 -k 1 -e 0.00001 >diamond_blastp.outfmt6

            echo -e "\\n-- Done with Diamond (blastp) --\\n"

            echo -e "\\n-- Starting HMMER --\\n"

            hmmscan --cpu ${task.cpus} --domtblout pfam.domtblout ${params.pfloc} \$fname.transdecoder_dir/longest_orfs.pep

            echo -e "\\n-- Done with HMMER --\\n"

            echo -e "\\n-- TransDecoder.Predict... --\\n"

            TransDecoder.Predict -t ${assembly} --retain_pfam_hits pfam.domtblout --retain_blastp_hits diamond_blastp.outfmt6

            echo -e "\\n-- Done with TransDecoder.Predict --\\n"

            echo -e "\\n-- Calculating statistics... --\\n"

            #Calculate statistics of Transdecoder
            echo "- Transdecoder (long, with homology) stats for ${sample_id}" >>${sample_id}.transdecoder.stats
            orfnum=\$( cat ${sample_id}.*.transdecoder.pep | grep -c ">" )
            echo -e "Total number of ORFs: \$orfnum \\n" >>${sample_id}.transdecoder.stats
            echo -e "\\t Of these ORFs" >>${sample_id}.transdecoder.stats
            orfnum=\$( cat ${sample_id}.*.transdecoder.pep | grep ">" | grep -c "|" )
            echo -e "\\t\\t with annotations: \$orfnum \\n" >>${sample_id}.transdecoder.stats
            orfnum=\$( cat ${sample_id}.*.transdecoder.pep | grep ">" | grep -v "|" | grep -c ">" )
            echo -e "\\t\\t no annotation: \$orfnum \\n" >>${sample_id}.transdecoder.stats
            orfnum=\$( cat ${sample_id}.*.transdecoder.pep | grep -c "ORF type:complete" )
            echo -e "\\t ORFs type=complete: \$orfnum \\n" >>${sample_id}.transdecoder.stats
            orfnum=\$( cat ${sample_id}.*.transdecoder.pep | grep "ORF type:complete" | grep -c "|" )
            echo -e "\\t\\t with annotations: \$orfnum \\n" >>${sample_id}.transdecoder.stats
            orfnum=\$( cat ${sample_id}.*.transdecoder.pep | grep -c "ORF type:5prime_partial" )
            echo -e "\\t ORFs type=5prime_partial: \$orfnum \\n" >>${sample_id}.transdecoder.stats
            orfnum=\$( cat ${sample_id}.*.transdecoder.pep | grep "ORF type:5prime_partial" | grep -c "|" )
            echo -e "\\t\\t with annotations: \$orfnum \\n" >>${sample_id}.transdecoder.stats
            orfnum=\$( cat ${sample_id}.*.transdecoder.pep | grep -c "ORF type:3prime_partial" )
            echo -e "\\t ORFs type=3prime_partial: \$orfnum \\n" >>${sample_id}.transdecoder.stats
            orfnum=\$( cat ${sample_id}.*.transdecoder.pep | grep "ORF type:3prime_partial" | grep -c "|" )
            echo -e "\\t\\t with annotations: \$orfnum \\n" >>${sample_id}.transdecoder.stats
            orfnum=\$( cat ${sample_id}.*.transdecoder.pep | grep -c "ORF type:internal" )
            echo -e "\\t ORFs type=internal: \$orfnum \\n" >>${sample_id}.transdecoder.stats
            orfnum=\$( cat ${sample_id}.*.transdecoder.pep | grep "ORF type:internal" | grep -c "|" )
            echo -e "\\t\\t with annotations: \$orfnum \\n" >>${sample_id}.transdecoder.stats

            echo -e "\\n-- Done with statistics --\\n"

            echo -e "\\n-- DONE with TransDecoder --\\n"
            """
        }
    }

    process swiss_diamond_trinotate_OA {

        label 'big_cpus'

        tag "${sample_id}"

        input:
            tuple sample_id, file(assembly) from assembly_ch_diamond_OA
            tuple sample_id, file("${sample_id}.*.transdecoder.pep") from transdecoder_ch_diamond_OA

        output:
            tuple sample_id, file("${sample_id}.diamond_blastx.outfmt6"), file("${sample_id}.diamond_blastp.outfmt6") into trinotate_ch_diamond_OA

        script:
            """
            swissdb=${params.mypwd}/DBs/diamonddb_swiss/uniprot_sprot.pep

            #Diamond (BLAST) Homologies

            echo -e "\\n-- Starting with Diamond (blastx) --\\n"

            diamond blastx -d \$swissdb -q ${assembly} -p ${task.cpus} -f 6 -k 1 -e 0.001 >${sample_id}.diamond_blastx.outfmt6

            echo -e "\\n-- Done with Diamond (blastx) --\\n"

            echo -e "\\n-- Starting with Diamond (blastp) --\\n"

            diamond blastp -d \$swissdb -q ${sample_id}.*.transdecoder.pep -p ${task.cpus} -f 6 -k 1 -e 0.001 >${sample_id}.diamond_blastp.outfmt6

            echo -e "\\n-- Done with Diamond (blastp)  --\\n"
            """
    }

    process custom_diamond_trinotate_OA {

        label 'big_cpus'

        tag "${sample_id}"

        input:
            tuple sample_id, file(assembly) from assembly_ch_diamond_custom_OA
            tuple sample_id, file("${sample_id}.*.transdecoder.pep") from transdecoder_ch_diamond_custom_OA

        output:
            tuple sample_id, file("${sample_id}.custom.diamond_blastx.outfmt6"), file("${sample_id}.custom.diamond_blastp.outfmt6") into trinotate_ch_diamond_custom_OA

        script:
            """
            unidb=${params.mypwd}/DBs/diamonddb_custom/${params.uniname}

            #Diamond (BLAST) Homologies

            echo -e "\\n-- Starting with Diamond (blastx) --\\n"

            diamond blastx -d \$unidb -q ${assembly} -p ${task.cpus} -f 6 -k 1 -e 0.001 >${sample_id}.custom.diamond_blastx.outfmt6

            echo -e "\\n-- Done with Diamond (blastx) --\\n"

            echo -e "\\n-- Starting with Diamond (blastp) --\\n"

            diamond blastp -d \$unidb -q ${sample_id}.*.transdecoder.pep -p ${task.cpus} -f 6 -k 1 -e 0.001 >${sample_id}.custom.diamond_blastp.outfmt6

            echo -e "\\n-- Done with Diamond (blastp)  --\\n"
            """
    }

    process hmmer_trinotate_OA {

        label 'low_cpus'

        tag "${sample_id}"

        input:
            tuple sample_id, file("${sample_id}.*.transdecoder.pep") from transdecoder_ch_hmmer_OA

        output:
            tuple sample_id, file("${sample_id}.TrinotatePFAM.out") into trinotate_ch_hmmer_OA

        script:
            """

            echo -e "\\n-- Starting with HMMER --\\n"

            hmmscan --cpu ${task.cpus} --domtblout ${sample_id}.TrinotatePFAM.out ${params.pfloc} ${sample_id}.*.transdecoder.pep >pfam.log

            echo -e "\\n-- Done with HMMER --\\n"
            """
    }

    process signalP_trinotate_OA {

        label 'low_cpus'

        tag "${sample_id}"

        input:
            tuple sample_id, file("${sample_id}.*.transdecoder.pep") from transdecoder_ch_signalp_OA

        output:
            tuple sample_id, file("${sample_id}.signalp.out") into trinotate_ch_signalp_OA

        script:
            """
            #signalP to predict signal peptides

            echo -e "\\n-- Starting with SignalP --\\n"

            ${params.signalp} -f short -n ${sample_id}.signalp.out ${sample_id}.*.transdecoder.pep

            echo -e "\\n-- Done with SignalP --\\n"
            """
    }

    process tmhmm_trinotate_OA {

        label 'low_cpus'

        tag "${sample_id}"

        input:
            tuple sample_id, file("${sample_id}.*.transdecoder.pep") from transdecoder_ch_tmhmm_OA

        output:
            tuple sample_id, file("${sample_id}.tmhmm.out") into trinotate_ch_tmhmm_OA

        script:
            """
            #tmHMM to predict transmembrane regions

            echo -e "\\n-- Starting with tmHMM --\\n"

            ${params.tmhmm} --short < ${sample_id}.*.transdecoder.pep >${sample_id}.tmhmm.out

            echo -e "\\n-- Done with tmHMM --\\n"
            """
    }

    process rnammer_trinotate_OA {

        label 'low_cpus'

        tag "${sample_id}"

        input:
            tuple sample_id, file(assembly) from assembly_ch_rnammer_OA

        output:
            tuple sample_id, file("${assembly}.rnammer.gff") into trinotate_ch_rnammer_OA

        script:
            """
            set +e
            #RNAMMER to identify rRNA transcripts

            echo -e "\\n-- Starting with RNAMMER --\\n"

            RnammerTranscriptome.pl --transcriptome ${assembly} --path_to_rnammer ${params.rnam}

            echo -e "\\n-- Done with RNAMMER --\\n"
            """
    }

    process trinotate_OA {

        label 'low_cpus'

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/trinotate", mode: "copy", overwrite: true

        input:
            tuple sample_id, file(assembly) from assembly_ch_trinotate_OA
            tuple sample_id, file("${sample_id}.*.transdecoder.pep") from transdecoder_ch_trinotate_OA
            tuple sample_id, file("${sample_id}.diamond_blastx.outfmt6"), file("${sample_id}.diamond_blastp.outfmt6") from trinotate_ch_diamond_OA
            tuple sample_id, file("${sample_id}.custom.diamond_blastx.outfmt6"), file("${sample_id}.custom.diamond_blastp.outfmt6") from trinotate_ch_diamond_custom_OA
            tuple sample_id, file("${sample_id}.TrinotatePFAM.out") from trinotate_ch_hmmer_OA
            tuple sample_id, file("${sample_id}.signalp.out") from trinotate_ch_signalp_OA
            tuple sample_id, file("${sample_id}.tmhmm.out") from trinotate_ch_tmhmm_OA
            tuple sample_id, file("${assembly}.rnammer.gff") from trinotate_ch_rnammer_OA

        output:
            tuple sample_id, file("${sample_id}.GO.terms.txt") into trinotate_summary_OA
            tuple sample_id, file("${sample_id}.trinotate_annotation_report.xls") into ( trinotate_ch_OA, custom_uniprot_ch_OA )
            tuple sample_id, file("*.terms.txt") into other_files_OA

        script:
            """
            #Generate gene_trans_map
            #Not using get_Trinity_gene_to_trans_map.pl since all the names are uniq
            cat ${assembly} | awk '{print \$1}' | grep ">" | cut -c 2- >a.txt

            paste a.txt a.txt >${assembly}.gene_trans_map

            #Get Trinotate.sqlite from folder (original)
            cp ${params.Tsql} .
            sqlname=`echo ${params.Tsql} | tr "\\/" "\\n" | grep "\\.sqlite"`

            echo -e "\\n-- Running Trinotate --\\n"

            Trinotate \$sqlname init --gene_trans_map ${assembly}.gene_trans_map --transcript_fasta ${assembly} --transdecoder_pep ${sample_id}.*.transdecoder.pep

            echo -e "\\n-- Ending run of Trinotate --\\n"

            echo -e "\\n-- Loading hits and predictions to sqlite database... --\\n"

            #Load protein hits
            Trinotate \$sqlname LOAD_swissprot_blastp ${sample_id}.diamond_blastp.outfmt6

            #Load transcript hits
            Trinotate \$sqlname LOAD_swissprot_blastx ${sample_id}.diamond_blastx.outfmt6

            #Load custom protein hits
            Trinotate \$sqlname LOAD_custom_blast --outfmt6 ${sample_id}.custom.diamond_blastp.outfmt6 --prog blastp --dbtype ${params.uniname}

            #Load custom transcript hits
            Trinotate \$sqlname LOAD_custom_blast --outfmt6 ${sample_id}.custom.diamond_blastx.outfmt6 --prog blastx --dbtype ${params.uniname}

            #Load Pfam domain entries
            Trinotate \$sqlname LOAD_pfam ${sample_id}.TrinotatePFAM.out

            #Load transmembrane domains
            if [ -s ${sample_id}.tmhmm.out ];then
                Trinotate \$sqlname LOAD_tmhmm ${sample_id}.tmhmm.out
            else
                echo "No transmembrane domains (tmhmm)"
            fi

            #Load signal peptide predictions
            if [ -s ${sample_id}.signalp.out ];then
                Trinotate \$sqlname LOAD_signalp ${sample_id}.signalp.out
            else
                echo "No Signal-P"
            fi

            #Load rnammer results
            if [ -s ${sample_id}.combined.okay.fa.rnammer.gff ];then
                Trinotate \$sqlname LOAD_rnammer ${assembly}.rnammer.gff
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

            cat ${sample_id}.trinotate_annotation_report.xls | cut -f 1,14 | grep "KEGG" | tr "\\`" ";" | grep "KO:K" | sed 's/\\tKEGG/\\t#KEGG/g' | sed 's/KO:/KO:#/g' | cut -f 1,3 -d "#" | tr -d "#" | sed 's/;KO://g' >${sample_id}.KEGG.terms.txt

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

        publishDir "${workDir}/${params.outdir}/stats", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.transdecoder.stats") from transdecoder_summary_OA

        output:
            tuple sample_id, file("${sample_id}.sum_transdecoder.txt") into final_sum_3_OA

        script:
            """
            #Summary of Transdecoder stats
            echo -e "Summary of Transdecoder \\n" >>${sample_id}.sum_transdecoder.txt
            cat ${sample_id}.transdecoder.stats >>${sample_id}.sum_transdecoder.txt
            echo -e "##### \\n" >>${sample_id}.sum_transdecoder.txt
            """
    }

    process summary_trinotate_individual_OA {

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/stats", mode: "copy", overwrite: true

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

        publishDir "${workDir}/${params.outdir}/figures/GO", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.trinotate_annotation_report.xls") from trinotate_ch_OA

        output:
            tuple sample_id, file("*.svg"), file("*.pdf"), file("*.txt") into go_fig_OA

        script:
            """
            cp ${params.mypwd}/bin/GO_plots.R .

            cat ${sample_id}.trinotate_annotation_report.xls | awk 'FS="\\t",OFS="#" {print \$1,\$15,\$16,\$17}' | grep -v "gene_id" >all_GOs.txt

            while read line;do
                echo ${line} | cut -f 2,3,4 -d "#" | grep "GO:" | tr "#" "\\n" | tr "\\`" "\\n" | sed 's/\\. /,/g' | tr "," "\\n" | grep "GO:" | sort -u >>final_GOs.txt
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
            """
    }

    process summary_custom_uniprot_OA {

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/figures/CustomUniProt", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.trinotate_annotation_report.xls") from custom_uniprot_ch_OA

        output:
            tuple sample_id, file("${sample_id}_custom_uniprot_hits.txt") into custom_uniprot_sum_OA
            tuple sample_id, file("${sample_id}_custom_uniprot_hits.svg"), file("${sample_id}_custom_uniprot_hits.pdf") into custom_uniprot_fig_OA

        script:
            """
            #get custom blast hits
            cat ${sample_id}.trinotate_annotation_report.xls | cut -f 8 | grep [A-Z] | grep "|" | tr "\\`" "\n" | \
                cut -f 1 -d "^" | cut -f 3 -d "|" | cut -f 2 -d "_" >a.txt

            cat ${sample_id}.trinotate_annotation_report.xls | cut -f 9 | grep [A-Z] | grep "|" | tr "\\`" "\n" | \
                cut -f 1 -d "^" | cut -f 3 -d "|" | cut -f 2 -d "_" >b.txt

            cat a.txt b.txt | sort | uniq -c | sort -nr | head -n 20 | awk 'OFS="," {print \$1,\$2}' >${sample_id}_custom_uniprot_hits.txt

            rm a.txt b.txt

            cp ${params.mypwd}/conf/uni_tax.txt .
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

            cp ${params.mypwd}/bin/custom_uniprot_hits.R .
            Rscript custom_uniprot_hits.R ${sample_id}
            """
    }

} else if (params.all) {

    println("\n\tRunning the full TransPi analysis\n")

    process custom_diamond_db {
        script:
            """
            cd ${params.mypwd}
            echo -e "-- Checking if Diamond database folder is present --\\n"
            if [ ! -d DBs/diamonddb_custom/ ];then
                echo -e "-- Folder is not present, creating one and the Diamond database --\\n"
                mkdir -p DBs/diamonddb_custom/
                cd DBs/diamonddb_custom
                cp ${params.uniprot} .
                diamond makedb --in ${params.uniname} -d ${params.uniname}
                export unidb=`pwd`/${params.uniname}
                cd ../
            elif [ -d DBs/diamonddb_custom/ ];then
                echo -e "-- Folder is present. Checking if Diamond database is built --\\n"
                cd DBs/diamonddb_custom
                if [ ! -e ${params.uniname}.dmnd ];then
                    echo -e "-- Diamond database not present, creating one --\\n"
                    cp ${params.uniprot} .
                    diamond makedb --in ${params.uniname} -d ${params.uniname}
                    export unidb=`pwd`/${params.uniname}
                elif [ -e ${params.uniname}.dmnd  ];then
                    echo -e "-- Diamond database already created --\\n"
                    export unidb=`pwd`/${params.uniname}
                fi
                cd ../
            fi
            """
    }

    process hmmer_db {
        script:
            """
            cd ${params.mypwd}
            echo -e "-- Checking if HMMER database folder is present --\\n"
            if [ -d DBs/hmmerdb/ ];then
                echo -e "-- Folder is present. Checking if HMMER database is built --\\n"
                cd DBs/hmmerdb
                if [ ! -e ${params.pfname}.h3f ] && [ ! -e ${params.pfname}.h3i ] && [ ! -e ${params.pfname}.h3m ] && [ ! -e ${params.pfname}.h3p ];then
                    echo -e "-- HMMER database not present, creating one --\\n"
                    hmmpress ${params.pfname}
                    export pf=`pwd`/${params.pfname}
                elif [ -s ${params.pfname}.h3f ] && [ -s ${params.pfname}.h3i ] && [ -s ${params.pfname}.h3m ] && [ -s ${params.pfname}.h3p ];then
                    echo -e "-- HMMER database already created --\\n"
                    export pf=`pwd`/${params.pfname}
                fi
                cd ../
            fi
            """
    }

    process swiss_diamond_db {
        script:
            """
            cd ${params.mypwd}/DBs/sqlite_db
            if [ -e uniprot_sprot.pep ];then
                cd ${params.mypwd}
                if [ ! -d DBs/diamonddb_swiss/ ];then
                    echo -e "-- Folder is not present, creating one and the Diamond database --\\n"
                    mkdir -p DBs/diamonddb_swiss
                    cd DBs/diamonddb_swiss
                    cp ${params.mypwd}/DBs/sqlite_db/uniprot_sprot.pep .
                    diamond makedb --in uniprot_sprot.pep -d uniprot_sprot.pep
                    export swissdb=`pwd`/uniprot_sprot.pep
                elif [ -d DBs/diamonddb_swiss/ ];then
                    cd DBs/diamonddb_swiss
                    if [ ! -e uniprot_sprot.pep.dmnd ];then
                        echo -e "-- Diamond database not present, creating one --\\n"
                        cp ${params.mypwd}/DBs/sqlite_db/uniprot_sprot.pep .
                        diamond makedb --in uniprot_sprot.pep -d uniprot_sprot.pep
                        export swissdb=`pwd`/uniprot_sprot.pep
                    elif [ -e uniprot_sprot.pep.dmnd ];then
                        echo -e "-- Diamond database already created --\\n"
                        export swissdb=`pwd`/uniprot_sprot.pep
                    fi
                fi
            elif [ ! -e uniprot_sprot.pep ];then
                cd ${params.mypwd}
                if [ ! -d DBs/diamonddb_swiss/ ];then
                    echo -e "-- Folder is not present, creating one and the Diamond database --\\n"
                    mkdir -p DBs/diamonddb_swiss
                    cd DBs/diamonddb_swiss
                    wget http://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
                    EMBL_swissprot_parser.pl uniprot_sprot.dat.gz ind
                    rm ind.*
                    mv uniprot_sprot.dat.gz.pep uniprot_sprot.pep
                    diamond makedb --in uniprot_sprot.pep -d uniprot_sprot.pep
                    export swissdb=`pwd`/uniprot_sprot.pep
                elif [ -d DBs/diamonddb_swiss/ ];then
                    echo -e "-- Folder is present. Checking if Diamond database is built --\\n"
                    cd DBs/diamonddb_swiss
                    if [ ! -e uniprot_sprot.pep.dmnd ];then
                        if [ ! -e uniprot_sprot.pep ];then
                            echo -e "-- Diamond database not present, creating one --\\n"
                            wget http://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
                            EMBL_swissprot_parser.pl uniprot_sprot.dat.gz ind
                            rm ind.*
                            mv uniprot_sprot.dat.gz.pep uniprot_sprot.pep
                            diamond makedb --in uniprot_sprot.pep -d uniprot_sprot.pep
                            export swissdb=`pwd`/uniprot_sprot.pep
                        elif [ -e uniprot_sprot.pep ];then
                            diamond makedb --in uniprot_sprot.pep -d uniprot_sprot.pep
                            export swissdb=`pwd`/uniprot_sprot.pep
                        fi
                    elif [ -e uniprot_sprot.pep.dmnd ];then
                        echo -e "-- Diamond database already created --\\n"
                        export swissdb=`pwd`/uniprot_sprot.pep
                    fi
                fi
            fi
            """
    }

    if (!params.skipQC) {

        process fasqc {

            label 'low_cpus'

            tag "${sample_id}"

            publishDir "${workDir}/${params.outdir}/fastqc", mode: "copy", overwrite: true

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

    if (!params.skipFilter) {

        process fastp {

            label 'med_cpus' // check this later

            tag "${sample_id}"

            publishDir "${workDir}/${params.outdir}/filter", mode: "copy", overwrite: true, pattern: "*.fastp.{json,html}"

            input:
                tuple sample_id, file(reads) from reads_ch

            output:
                tuple sample_id, file("*.fastp.{json,html}") into fastp_results
                tuple sample_id, file("*.filter.fq") into reads_ass_ch

            script:
                """
                echo ${sample_id}

                fastp -i ${reads[0]} -I ${reads[1]} -o left-${sample_id}.filter.fq -O right-${sample_id}.filter.fq --detect_adapter_for_pe \
                --average_qual ${params.minQual} --overrepresentation_analysis --html ${sample_id}.fastp.html --json ${sample_id}.fastp.json \
                --thread ${task.cpus} --report_title ${sample_id}
                """
        }
    } else {
        reads_ch
            .set{ reads_ass_ch }
        fastp_results = Channel.empty()
    }

    process normalize_reads {

        label 'med_mem'

        tag "${sample_id}"

        input:
            tuple sample_id, file(reads) from reads_ass_ch

        output:
            tuple sample_id, file("left-${sample_id}.norm.fq"), file("right-${sample_id}.norm.fq") into ( norm_reads_soap, norm_reads_velvet, norm_reads_trinity, norm_reads_spades, norm_reads_transabyss )

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

            mv left.norm.fq left-"${sample_id}".norm.fq
            mv right.norm.fq right-"${sample_id}".norm.fq
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

            mv left.norm.fq left-"${sample_id}".norm.fq
            mv right.norm.fq right-"${sample_id}".norm.fq

            rm left-${sample_id}.fq right-${sample_id}.fq
            """
        }
    }

    process trinity_assembly {

        label 'med_mem'

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/assemblies", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("left-${sample_id}.norm.fq"), file("right-${sample_id}.norm.fq") from norm_reads_trinity

        output:
            tuple sample_id, file("${sample_id}.Trinity.fa") into ( assemblies_ch_trinity, busco3_ch_trinity, busco4_ch_trinity )

        script:
            """
            mem=\$( echo ${task.memory} | cut -f 1 -d " " )

            Trinity --max_memory \${mem}G --seqType fq --left left-${sample_id}.norm.fq --right right-${sample_id}.norm.fq --CPU ${task.cpus} --no_normalize_reads --full_cleanup --output trinity_out_dir

            mv trinity_out_dir.Trinity.fasta ${sample_id}.Trinity.fa
            """
    }

    process soap_assembly {

        label 'med_mem'

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/assemblies", mode: "copy", overwrite: true

        input:
            val k from "${params.k}"
            tuple sample_id, file("left-${sample_id}.norm.fq"), file("right-${sample_id}.norm.fq") from norm_reads_soap

        output:
            tuple sample_id, file("${sample_id}.SOAP.fa") into assemblies_ch_soap

        script:
            """
            echo -e "\\n-- Generating SOAP config file --\\n"
            echo "max_rd_len="${params.max_rd_len} >>config.txt
            echo "[LIB]" >>config.txt
            echo "rd_len_cutof="${params.rd_len_cutof} >>config.txt
            #echo "avg_ins="${params.avg_ins} >>config.txt
            echo "reverse_seq="${params.reverse_seq} >>config.txt
            echo "asm_flags="${params.asm_flags} >>config.txt
            echo "map_len="${params.map_len} >>config.txt
            echo "q1="left-${sample_id}.norm.fq >>config.txt
            echo "q2="right-${sample_id}.norm.fq >>config.txt

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

    process velvet_oases_assembly {

        label 'med_mem'

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/assemblies", mode: "copy", overwrite: true

        input:
            val k from "${params.k}"
            tuple sample_id, file("left-${sample_id}.norm.fq"), file("right-${sample_id}.norm.fq") from norm_reads_velvet

        output:
            tuple sample_id, file("${sample_id}.Velvet.fa") into assemblies_ch_velvet

        script:
            """
    	    echo -e "\\n-- Starting with Velveth --\\n"
            for x in `echo $k | tr "," " "`;do
                echo -e "\\n-- k\${x} --\\n"
                velveth oases.\${x} \${x} -shortPaired -fastq -separate left-${sample_id}.norm.fq right-${sample_id}.norm.fq
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

    process rna_spades_assembly {

        label 'med_mem'

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/assemblies", mode: "copy", overwrite: true

        input:
            val k from "${params.k}"
            tuple sample_id, file("left-${sample_id}.norm.fq"), file("right-${sample_id}.norm.fq") from norm_reads_spades

        output:
            tuple sample_id, file("${sample_id}.SPADES.fa") into assemblies_ch_spades

        script:
            """
            echo -e "\\n-- Starting rnaSPADES assemblies --\\n"

            mem=\$( echo ${task.memory} | cut -f 1 -d " " )

            for x in `echo $k | tr "," " "`;do
                echo -e "\\n-- rnaSPADES k\${x} --\\n"

                rnaspades.py -1 left-${sample_id}.norm.fq -2 right-${sample_id}.norm.fq -o ${sample_id}_spades_\${x} -t ${task.cpus} -k \${x} -m \${mem}

            done

            echo -e "\\n-- Finished with the assemblies --\\n"

            for x in `echo $k | tr "," " "`;do
                sed -i "s/>/>SPADES.k\${x}./g" ${sample_id}_spades_\${x}/transcripts.fasta
            done

            cat ${sample_id}_spades_*/transcripts.fasta >${sample_id}.SPADES.fa

            rm -rf ${sample_id}_spades_*
            """
    }

    process transabyss_assembly {

        label 'med_mem'

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/assemblies", mode: "copy", overwrite: true

        input:
            val k from "${params.k}"
            tuple sample_id, file("left-${sample_id}.norm.fq"), file("right-${sample_id}.norm.fq") from norm_reads_transabyss

        output:
            tuple sample_id, file("${sample_id}.TransABySS.fa") into assemblies_ch_transabyss

        script:
            """
            echo -e "\\n-- Starting Trans-ABySS assemblies --\\n"

            for x in `echo $k | tr "," " "`;do
                echo -e "\\n-- Trans-ABySS k\${x} --\\n"

                transabyss -k \${x} --pe left-${sample_id}.norm.fq right-${sample_id}.norm.fq --outdir ${sample_id}_transabyss_\${x} --name k\${x}.transabyss.fa --threads ${task.cpus} -c 12 --length 200

            done

            echo -e "\\n-- Finished with the assemblies --\\n"

            for x in `echo $k | tr "," " "`;do
                sed -i "s/>/>TransABySS.k\${x}./g" ${sample_id}_transabyss_\${x}/k\${x}.transabyss.fa-final.fa
            done

            cat ${sample_id}_transabyss_*/k*.transabyss.fa-final.fa >${sample_id}.TransABySS.fa

            rm -rf ${sample_id}_transabyss_*
            """
    }

    process evigene {

        label 'med_mem'

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/evigene", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.Trinity.fa") from assemblies_ch_trinity
            tuple sample_id, file("${sample_id}.SOAP.fa") from assemblies_ch_soap
            tuple sample_id, file("${sample_id}.Velvet.fa") from assemblies_ch_velvet
            tuple sample_id, file("${sample_id}.SPADES.fa") from assemblies_ch_spades
            tuple sample_id, file("${sample_id}.TransABySS.fa") from assemblies_ch_transabyss

        output:
            tuple sample_id, file("${sample_id}.combined.okay.fa") into ( evigene_ch_busco3, evigene_ch_busco4, evigene_ch_transdecoder, evigene_ch_diamond, evigene_ch_rnammer, evigene_ch_trinotate, evigene_ch_trinotate_custom )
            tuple sample_id, file("${sample_id}.combined.fa"), file("${sample_id}.combined.okay.fa") into evigene_summary

        script:
            def mem_MB=(task.memory.toMega())

            """
            echo -e "\\n-- Starting EviGene --\\n"

            cat *.fa >${sample_id}.combined.fa

            $evi/scripts/prot/tr2aacds.pl -tidy -NCPU ${task.cpus} -MAXMEM ${mem_MB} -log -cdna ${sample_id}.combined.fa

            echo -e "\\n-- DONE with EviGene --\\n"

            cp okayset/${sample_id}.combined.okay.combined.fa ${sample_id}.combined.okay.fa
            cp okayset/${sample_id}.combined.okay.combined.cds ${sample_id}.combined.okay.cds


            if [ -d tmpfiles/ ];then
                rm -rf tmpfiles/
            fi
        	"""
    }

    process busco3 {

        label 'med_cpus'

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/busco3", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.combined.okay.fa") from evigene_ch_busco3

        output:
            tuple sample_id, file("run_${sample_id}.fa.bus") into busco3_ch
            tuple sample_id, file("*${sample_id}.fa.bus.txt") into ( busco3_summary, busco3_comp_1 )

        script:
            """
            echo -e "\\n-- Starting BUSCO --\\n"

            run_BUSCO.py -i ${sample_id}.combined.okay.fa -o ${sample_id}.fa.bus -l ${params.busco3db} -m tran -c ${task.cpus}

            echo -e "\\n-- DONE with BUSCO --\\n"

            cp run_${sample_id}.fa.bus/short_summary_${sample_id}.fa.bus.txt .
            """
    }

    process busco3_tri {

        label 'med_cpus'

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/busco3", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.Trinity.fa") from busco3_ch_trinity

        output:
            tuple sample_id, file("*${sample_id}.Trinity.fa.bus.txt") into ( busco3_ch_trinity_sum, busco3_comp_2 )

        script:
            """
            echo -e "\\n-- Starting BUSCO --\\n"

            run_BUSCO.py -i ${sample_id}.Trinity.fa -o ${sample_id}.Trinity.fa.bus -l ${params.busco3db} -m tran -c ${task.cpus}

            echo -e "\\n-- DONE with BUSCO --\\n"

            cp run_${sample_id}.Trinity.fa.bus/short_summary_${sample_id}.Trinity.fa.bus.txt .
            """
    }

    process busco4 {

        conda "${params.cenv}"

        label 'med_cpus'

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/busco4", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.combined.okay.fa") from evigene_ch_busco4

        output:
            tuple sample_id, file("${sample_id}.fa.bus") into busco4_ch
            tuple sample_id, file("*${sample_id}.fa.bus.txt") into ( busco4_summary, busco4_comp_1 )

        script:
            """
            echo -e "\\n-- Starting BUSCO --\\n"

            busco -i ${sample_id}.combined.okay.fa -o ${sample_id}.fa.bus -l ${params.busco4db} -m tran -c ${task.cpus} --offline

            echo -e "\\n-- DONE with BUSCO --\\n"

            cp ${sample_id}.fa.bus/short_summary.*.${sample_id}.fa.bus.txt .
            """
    }

    process busco4_tri {

        conda "${params.cenv}"

        label 'med_cpus'

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/busco4", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.Trinity.fa") from busco4_ch_trinity

        output:
            tuple sample_id, file("*${sample_id}.Trinity.fa.bus.txt") into ( busco4_ch_trinity_sum, busco4_comp_2 )

        script:
            """
            echo -e "\\n-- Starting BUSCO --\\n"

            busco -i ${sample_id}.Trinity.fa -o ${sample_id}.Trinity.fa.bus -l ${params.busco4db} -m tran -c ${task.cpus} --offline

            echo -e "\\n-- DONE with BUSCO --\\n"

            cp ${sample_id}.Trinity.fa.bus/short_summary.*.${sample_id}.Trinity.fa.bus.txt .
            """
    }

    process transdecoder {

        label 'med_cpus'

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/transdecoder", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.combined.okay.fa") from evigene_ch_transdecoder

        output:
            tuple sample_id, file("${sample_id}.combined.okay.fa.transdecoder.pep") into ( transdecoder_ch_diamond, transdecoder_ch_hmmer, transdecoder_ch_signalp, transdecoder_ch_tmhmm, transdecoder_ch_trinotate, transdecoder_ch_diamond_custom )
            tuple sample_id, file("${sample_id}.transdecoder.stats") into transdecoder_summary
            tuple sample_id, file("${sample_id}.*.transdecoder.{cds,gff,bed}") into transdecoder_files

        script:
        if (params.shortTransdecoder) {
            """
            echo -e "\\n-- TransDecoder.LongOrfs... --\\n"

            TransDecoder.LongOrfs -t ${sample_id}.combined.okay.fa"

            echo -e "\\n-- Done with TransDecoder.LongOrfs --\\n"

            echo -e "\\n-- TransDecoder.Predict... --\\n"

            TransDecoder.Predict -t ${sample_id}.combined.okay.fa"

            echo -e "\\n-- Done with TransDecoder.Predict --\\n"

            echo -e "\\n-- Calculating statistics... --\\n"

            #Calculate statistics of Transdecoder
            echo "- Transdecoder (short,no homolgy) stats for ${sample_id}" >>${sample_id}.transdecoder.stats
            orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep -c ">" )
            echo -e "Total number of ORFs: \$orfnum \\n" >>${sample_id}.transdecoder.stats
            orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep -c "ORF type:complete" )
            echo -e "\\t ORFs type=complete: \$orfnum \\n" >>${sample_id}.transdecoder.stats
            orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep -c "ORF type:5prime_partial" )
            echo -e "\\t ORFs type=5prime_partial: \$orfnum \\n" >>${sample_id}.transdecoder.stats
            orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep -c "ORF type:3prime_partial" )
            echo -e "\\t ORFs type=3prime_partial: \$orfnum \\n" >>${sample_id}.transdecoder.stats
            orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep -c "ORF type:internal" )
            echo -e "\\t ORFs type=internal: \$orfnum \\n">>${sample_id}.transdecoder.stats

            echo -e "\\n-- Done with statistics --\\n"

            echo -e "\\n-- DONE with TransDecoder --\\n"
            """

        } else {
            """
            unidb=${params.mypwd}/DBs/diamonddb_custom/${params.uniname}

            echo -e "\\n-- TransDecoder.LongOrfs... --\\n"

            TransDecoder.LongOrfs -t ${sample_id}.combined.okay.fa

            echo -e "\\n-- Done with TransDecoder.LongOrfs --\\n"

            fname=${sample_id}.combined.okay.fa

            echo -e "\\n-- Starting Diamond (blastp) --\\n"

            diamond blastp -d \$unidb -q \$fname.transdecoder_dir/longest_orfs.pep -p ${task.cpus} -f 6 -k 1 -e 0.00001 >diamond_blastp.outfmt6

            echo -e "\\n-- Done with Diamond (blastp) --\\n"

            echo -e "\\n-- Starting HMMER --\\n"

            hmmscan --cpu ${task.cpus} --domtblout pfam.domtblout ${params.pfloc} \$fname.transdecoder_dir/longest_orfs.pep

            echo -e "\\n-- Done with HMMER --\\n"

            echo -e "\\n-- TransDecoder.Predict... --\\n"

            TransDecoder.Predict -t ${sample_id}.combined.okay.fa --retain_pfam_hits pfam.domtblout --retain_blastp_hits diamond_blastp.outfmt6

            echo -e "\\n-- Done with TransDecoder.Predict --\\n"

            echo -e "\\n-- Calculating statistics... --\\n"

            #Calculate statistics of Transdecoder
            echo "- Transdecoder (long, with homology) stats for ${sample_id}" >>${sample_id}.transdecoder.stats
            orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep -c ">" )
            echo -e "Total number of ORFs: \$orfnum \\n" >>${sample_id}.transdecoder.stats
            echo -e "\\t Of these ORFs" >>${sample_id}.transdecoder.stats
            orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep ">" | grep -c "|" )
            echo -e "\\t\\t with annotations: \$orfnum \\n" >>${sample_id}.transdecoder.stats
            orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep ">" | grep -v "|" | grep -c ">" )
            echo -e "\\t\\t no annotation: \$orfnum \\n" >>${sample_id}.transdecoder.stats
            orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep -c "ORF type:complete" )
            echo -e "\\t ORFs type=complete: \$orfnum \\n" >>${sample_id}.transdecoder.stats
            orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep "ORF type:complete" | grep -c "|" )
            echo -e "\\t\\t with annotations: \$orfnum \\n" >>${sample_id}.transdecoder.stats
            orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep -c "ORF type:5prime_partial" )
            echo -e "\\t ORFs type=5prime_partial: \$orfnum \\n" >>${sample_id}.transdecoder.stats
            orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep "ORF type:5prime_partial" | grep -c "|" )
            echo -e "\\t\\t with annotations: \$orfnum \\n" >>${sample_id}.transdecoder.stats
            orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep -c "ORF type:3prime_partial" )
            echo -e "\\t ORFs type=3prime_partial: \$orfnum \\n" >>${sample_id}.transdecoder.stats
            orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep "ORF type:3prime_partial" | grep -c "|" )
            echo -e "\\t\\t with annotations: \$orfnum \\n" >>${sample_id}.transdecoder.stats
            orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep -c "ORF type:internal" )
            echo -e "\\t ORFs type=internal: \$orfnum \\n" >>${sample_id}.transdecoder.stats
            orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep "ORF type:internal" | grep -c "|" )
            echo -e "\\t\\t with annotations: \$orfnum \\n" >>${sample_id}.transdecoder.stats

            echo -e "\\n-- Done with statistics --\\n"

            echo -e "\\n-- DONE with TransDecoder --\\n"
            """
        }
    }

    process swiss_diamond_trinotate {

        label 'big_cpus'

        tag "${sample_id}"

        input:
            tuple sample_id, file("${sample_id}.combined.okay.fa") from evigene_ch_diamond
            tuple sample_id, file("${sample_id}.combined.okay.fa.transdecoder.pep") from transdecoder_ch_diamond

        output:
            tuple sample_id, file("${sample_id}.diamond_blastx.outfmt6"), file("${sample_id}.diamond_blastp.outfmt6") into trinotate_ch_diamond

        script:
            """
            swissdb=${params.mypwd}/DBs/diamonddb_swiss/uniprot_sprot.pep

            #Diamond (BLAST) Homologies

            echo -e "\\n-- Starting with Diamond (blastx) --\\n"

            diamond blastx -d \$swissdb -q ${sample_id}.combined.okay.fa -p ${task.cpus} -f 6 -k 1 -e 0.001 >${sample_id}.diamond_blastx.outfmt6

            echo -e "\\n-- Done with Diamond (blastx) --\\n"

            echo -e "\\n-- Starting with Diamond (blastp) --\\n"

            diamond blastp -d \$swissdb -q ${sample_id}.combined.okay.fa.transdecoder.pep -p ${task.cpus} -f 6 -k 1 -e 0.001 >${sample_id}.diamond_blastp.outfmt6

            echo -e "\\n-- Done with Diamond (blastp)  --\\n"
            """
    }

    process custom_diamond_trinotate {

        label 'big_cpus'

        tag "${sample_id}"

        input:
            tuple sample_id, file("${sample_id}.combined.okay.fa") from evigene_ch_trinotate_custom
            tuple sample_id, file("${sample_id}.combined.okay.fa.transdecoder.pep") from transdecoder_ch_diamond_custom

        output:
            tuple sample_id, file("${sample_id}.custom.diamond_blastx.outfmt6"), file("${sample_id}.custom.diamond_blastp.outfmt6") into trinotate_ch_diamond_custom

        script:
            """
            unidb=${params.mypwd}/DBs/diamonddb_custom/${params.uniname}

            #Diamond (BLAST) Homologies

            echo -e "\\n-- Starting with Diamond (blastx) --\\n"

            diamond blastx -d \$unidb -q ${sample_id}.combined.okay.fa -p ${task.cpus} -f 6 -k 1 -e 0.001 >${sample_id}.custom.diamond_blastx.outfmt6

            echo -e "\\n-- Done with Diamond (blastx) --\\n"

            echo -e "\\n-- Starting with Diamond (blastp) --\\n"

            diamond blastp -d \$unidb -q ${sample_id}.combined.okay.fa.transdecoder.pep -p ${task.cpus} -f 6 -k 1 -e 0.001 >${sample_id}.custom.diamond_blastp.outfmt6

            echo -e "\\n-- Done with Diamond (blastp)  --\\n"
            """
    }

    process hmmer_trinotate {

        label 'low_cpus'

        tag "${sample_id}"

        input:
            tuple sample_id, file("${sample_id}.combined.okay.fa.transdecoder.pep") from transdecoder_ch_hmmer

        output:
            tuple sample_id, file("${sample_id}.TrinotatePFAM.out") into trinotate_ch_hmmer

        script:
            """

            echo -e "\\n-- Starting with HMMER --\\n"

            hmmscan --cpu ${task.cpus} --domtblout ${sample_id}.TrinotatePFAM.out ${params.pfloc} ${sample_id}.combined.okay.fa.transdecoder.pep >pfam.log

            echo -e "\\n-- Done with HMMER --\\n"
            """
    }

    process signalP_trinotate {

        label 'low_cpus'

        tag "${sample_id}"

        input:
            tuple sample_id, file("${sample_id}.combined.okay.fa.transdecoder.pep") from transdecoder_ch_signalp

        output:
            tuple sample_id, file("${sample_id}.signalp.out") into trinotate_ch_signalp

        script:
            """
            #signalP to predict signal peptides

            echo -e "\\n-- Starting with SignalP --\\n"

            ${params.signalp} -f short -n ${sample_id}.signalp.out ${sample_id}.combined.okay.fa.transdecoder.pep

            echo -e "\\n-- Done with SignalP --\\n"
            """
    }

    process tmhmm_trinotate {

        label 'low_cpus'

        tag "${sample_id}"

        input:
            tuple sample_id, file("${sample_id}.combined.okay.fa.transdecoder.pep") from transdecoder_ch_tmhmm

        output:
            tuple sample_id, file("${sample_id}.tmhmm.out") into trinotate_ch_tmhmm

        script:
            """
            #tmHMM to predict transmembrane regions

            echo -e "\\n-- Starting with tmHMM --\\n"

            ${params.tmhmm} --short < ${sample_id}.combined.okay.fa.transdecoder.pep >${sample_id}.tmhmm.out

            echo -e "\\n-- Done with tmHMM --\\n"
            """
    }

    process rnammer_trinotate {

        label 'low_cpus'

        tag "${sample_id}"

        input:
            tuple sample_id, file("${sample_id}.combined.okay.fa") from evigene_ch_rnammer

        output:
            tuple sample_id, file("${sample_id}.combined.okay.fa.rnammer.gff") into trinotate_ch_rnammer

        script:
            """
            set +e
            #RNAMMER to identify rRNA transcripts

            echo -e "\\n-- Starting with RNAMMER --\\n"

            RnammerTranscriptome.pl --transcriptome ${sample_id}.combined.okay.fa --path_to_rnammer ${params.rnam}

            echo -e "\\n-- Done with RNAMMER --\\n"
            """
    }

    process trinotate {

        label 'low_cpus'

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/trinotate", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.combined.okay.fa") from evigene_ch_trinotate
            tuple sample_id, file("${sample_id}.combined.okay.fa.transdecoder.pep") from transdecoder_ch_trinotate
            tuple sample_id, file("${sample_id}.diamond_blastx.outfmt6"), file("${sample_id}.diamond_blastp.outfmt6") from trinotate_ch_diamond
            tuple sample_id, file("${sample_id}.custom.diamond_blastx.outfmt6"), file("${sample_id}.custom.diamond_blastp.outfmt6") from trinotate_ch_diamond_custom
            tuple sample_id, file("${sample_id}.TrinotatePFAM.out") from trinotate_ch_hmmer
            tuple sample_id, file("${sample_id}.signalp.out") from trinotate_ch_signalp
            tuple sample_id, file("${sample_id}.tmhmm.out") from trinotate_ch_tmhmm
            tuple sample_id, file("${sample_id}.combined.okay.fa.rnammer.gff") from trinotate_ch_rnammer

        output:
            tuple sample_id, file("${sample_id}.GO.terms.txt") into trinotate_summary
            tuple sample_id, file("${sample_id}.trinotate_annotation_report.xls") into ( trinotate_ch, custom_uniprot_ch )
            tuple sample_id, file("*.terms.txt") into other_files

        script:
            """
            #Generate gene_trans_map
            #Not using get_Trinity_gene_to_trans_map.pl since all the names are uniq
            cat ${sample_id}.combined.okay.fa | awk '{print \$1}' | grep ">" | cut -c 2- >a.txt

            paste a.txt a.txt >${sample_id}.combined.okay.fa.gene_trans_map

            #Get Trinotate.sqlite from folder (original)
            cp ${params.Tsql} .
            sqlname=`echo ${params.Tsql} | tr "\\/" "\\n" | grep "\\.sqlite"`

            echo -e "\\n-- Running Trinotate --\\n"

            Trinotate \$sqlname init --gene_trans_map ${sample_id}.combined.okay.fa.gene_trans_map --transcript_fasta ${sample_id}.combined.okay.fa --transdecoder_pep ${sample_id}.combined.okay.fa.transdecoder.pep

            echo -e "\\n-- Ending run of Trinotate --\\n"

            echo -e "\\n-- Loading hits and predictions to sqlite database... --\\n"

            #Load protein hits
            Trinotate \$sqlname LOAD_swissprot_blastp ${sample_id}.diamond_blastp.outfmt6

            #Load transcript hits
            Trinotate \$sqlname LOAD_swissprot_blastx ${sample_id}.diamond_blastx.outfmt6

            #Load custom protein hits
            Trinotate \$sqlname LOAD_custom_blast --outfmt6 ${sample_id}.custom.diamond_blastp.outfmt6 --prog blastp --dbtype ${params.uniname}

            #Load custom transcript hits
            Trinotate \$sqlname LOAD_custom_blast --outfmt6 ${sample_id}.custom.diamond_blastx.outfmt6 --prog blastx --dbtype ${params.uniname}

            #Load Pfam domain entries
            Trinotate \$sqlname LOAD_pfam ${sample_id}.TrinotatePFAM.out

            #Load transmembrane domains
            if [ -s ${sample_id}.tmhmm.out ];then
                Trinotate \$sqlname LOAD_tmhmm ${sample_id}.tmhmm.out
            else
                echo "No transmembrane domains (tmhmm)"
            fi

            #Load signal peptide predictions
            if [ -s ${sample_id}.signalp.out ];then
                Trinotate \$sqlname LOAD_signalp ${sample_id}.signalp.out
            else
                echo "No Signal-P"
            fi

            #Load rnammer results
            if [ -s ${sample_id}.combined.okay.fa.rnammer.gff ];then
                Trinotate \$sqlname LOAD_rnammer ${sample_id}.combined.okay.fa.rnammer.gff
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

        publishDir "${workDir}/${params.outdir}/stats", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.combined.fa"), file("${sample_id}.combined.okay.fa") from evigene_summary

        output:
            tuple sample_id, file("${sample_id}.sum_preEG.txt"), file("${sample_id}.sum_EG.txt") into final_sum_1

        script:
            """
            #Summary of total number of transcripts
            echo -e "- Number of transcripts before Evidential Genes\\n" >>${sample_id}.sum_preEG.txt
            echo -e "- Individual ${sample_id} \\n" >>${sample_id}.sum_preEG.txt
            echo -e "\\t Total transcripts:" >>${sample_id}.sum_preEG.txt
            num=\$( cat ${sample_id}.combined.fa | grep -c ">" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}.sum_preEG.txt
            echo -e "\\t Trinity" >>${sample_id}.sum_preEG.txt
            num=\$( cat ${sample_id}.combined.fa | grep -c ">TRINITY" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}.sum_preEG.txt
            echo -e "\\t SOAP" >>${sample_id}.sum_preEG.txt
            num=\$( cat ${sample_id}.combined.fa | grep -c ">SOAP" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}.sum_preEG.txt
            echo -e "\\t Velvet/Oases" >>${sample_id}.sum_preEG.txt
            num=\$( cat ${sample_id}.combined.fa | grep -c ">Velvet" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}.sum_preEG.txt
            echo -e "\\t rna-SPADES" >>${sample_id}.sum_preEG.txt
            num=\$( cat ${sample_id}.combined.fa | grep -c ">SPADES" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}.sum_preEG.txt
            echo -e "\\t Trans-ABySS" >>${sample_id}.sum_preEG.txt
            num=\$( cat ${sample_id}.combined.fa | grep -c ">TransABySS" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}.sum_preEG.txt

            #Summary of transcripts after EvidentialGenes
            echo -e "- Number of transcripts by individual after EvidentialGenes\\n" >>${sample_id}.sum_EG.txt
            echo -e "- Individual ${sample_id} \\n" >>${sample_id}.sum_EG.txt
            echo -e "\\t Total transcripts:" >>${sample_id}.sum_EG.txt
            num=\$( cat ${sample_id}.combined.okay.fa | grep -c ">" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}.sum_EG.txt
            echo -e "\\t Trinity" >>${sample_id}.sum_EG.txt
            num=\$( cat ${sample_id}.combined.okay.fa | grep -c ">TRINITY" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}.sum_EG.txt
            echo -e "\\t SOAP" >>${sample_id}.sum_EG.txt
            num=\$( cat ${sample_id}.combined.okay.fa | grep -c ">SOAP" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}.sum_EG.txt
            echo -e "\\t Velvet/Oases" >>${sample_id}.sum_EG.txt
            num=\$( cat ${sample_id}.combined.okay.fa | grep -c ">Velvet" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}.sum_EG.txt
            echo -e "\\t rna-SPADES" >>${sample_id}.sum_EG.txt
            num=\$( cat ${sample_id}.combined.okay.fa | grep -c ">SPADES" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}.sum_EG.txt
            echo -e "\\t Trans-ABySS" >>${sample_id}.sum_EG.txt
            num=\$( cat ${sample_id}.combined.okay.fa | grep -c ">TransABySS" )
            echo -e "\\t\\t \$num \\n" >>${sample_id}.sum_EG.txt
            """
    }

    process summary_busco3_individual {

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/stats", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.fa.bus.txt") from busco3_summary
            tuple sample_id, file("${sample_id}.Trinity.fa.bus.txt") from busco3_ch_trinity_sum

        output:
            tuple sample_id, file("${sample_id}.sum_busco3.txt") into final_sum_2v3

        script:
            """
            #Summary of BUSCO scores for the final_assemblies
            echo -e "Summary of BUSCO V3 \\n" >>${sample_id}.sum_busco3.txt
            echo "-- TransPi BUSCO V3 scores -- " >>${sample_id}.sum_busco3.txt
            cat ${sample_id}.fa.bus.txt >>${sample_id}.sum_busco3.txt
            echo -e "\\n-- Trinity BUSCO V3 scores --" >>${sample_id}.sum_busco3.txt
            cat ${sample_id}.Trinity.fa.bus.txt >>${sample_id}.sum_busco3.txt
            """
    }

    process summary_busco4_individual {

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/stats", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.fa.bus.txt") from busco4_summary
            tuple sample_id, file("${sample_id}.Trinity.fa.bus.txt") from busco4_ch_trinity_sum

        output:
            tuple sample_id, file("${sample_id}.sum_busco4.txt") into final_sum_2v4

        script:
            """
            #Summary of BUSCO scores for the final_assemblies
            echo -e "Summary of BUSCO V4 \\n" >>${sample_id}.sum_busco4.txt
            echo "-- TransPi BUSCO V4 scores -- " >>${sample_id}.sum_busco4.txt
            cat ${sample_id}.fa.bus.txt >>${sample_id}.sum_busco4.txt
            echo -e "\\n-- Trinity BUSCO V4 scores --" >>${sample_id}.sum_busco4.txt
            cat ${sample_id}.Trinity.fa.bus.txt >>${sample_id}.sum_busco4.txt
            """
    }

    process summary_transdecoder_individual {

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/stats", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.transdecoder.stats") from transdecoder_summary

        output:
            tuple sample_id, file("${sample_id}.sum_transdecoder.txt") into final_sum_3

        script:
            """
            #Summary of Transdecoder stats
            echo -e "Summary of Transdecoder \\n" >>${sample_id}.sum_transdecoder.txt
            cat ${sample_id}.transdecoder.stats >>${sample_id}.sum_transdecoder.txt
            echo -e "##### \\n" >>${sample_id}.sum_transdecoder.txt
            """
    }

    process summary_trinotate_individual {

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/stats", mode: "copy", overwrite: true

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

    process get_busco3_comparison {

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/figures/BUSCO3", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.fa.bus.txt") from busco3_comp_1
            tuple sample_id, file("${sample_id}.Trinity.fa.bus.txt") from busco3_comp_2

        output:
            tuple sample_id, file("${sample_id}_BUSCO3_comparison.pdf"), file("${sample_id}_BUSCO3_comparison.svg") into busco3_fig

        script:
            """
            set +e
            bash get_busco_val.sh ${sample_id}.Trinity.fa.bus.txt ${sample_id}.fa.bus.txt v3
            cp ${params.mypwd}/bin/busco_comparison.R .
            a=\$( cat final_spec )
            sed -i "s/MYSPEC/\${a}/" busco_comparison.R
            b=\$( cat final_perc )
            sed -i "s/MYPERC/\${b}/" busco_comparison.R
            c=\$( cat final_num )
            sed -i "s/MYVAL/\${c}/" busco_comparison.R
            Rscript busco_comparison.R ${sample_id}
            mv ${sample_id}_BUSCO_comparison.pdf ${sample_id}_BUSCO3_comparison.pdf
            mv ${sample_id}_BUSCO_comparison.svg ${sample_id}_BUSCO3_comparison.svg
            """
    }

    process get_busco4_comparison {

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/figures/BUSCO4", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("{sample_id}.fa.bus.txt") from busco4_comp_1
            tuple sample_id, file("{sample_id}.Trinity.fa.bus.txt") from busco4_comp_2

        output:
            tuple sample_id, file("${sample_id}_BUSCO4_comparison.pdf"), file("${sample_id}_BUSCO4_comparison.svg") into busco4_fig

        script:
            """
            set +e
            bash get_busco_val.sh {sample_id}.Trinity.fa.bus.txt {sample_id}.fa.bus.txt v4
            cp ${params.mypwd}/bin/busco_comparison.R .
            a=\$( cat final_spec )
            sed -i "s/MYSPEC/\${a}/" busco_comparison.R
            b=\$( cat final_perc )
            sed -i "s/MYPERC/\${b}/" busco_comparison.R
            c=\$( cat final_num )
            sed -i "s/MYVAL/\${c}/" busco_comparison.R
            Rscript busco_comparison.R ${sample_id}
            mv ${sample_id}_BUSCO_comparison.pdf ${sample_id}_BUSCO4_comparison.pdf
            mv ${sample_id}_BUSCO_comparison.svg ${sample_id}_BUSCO4_comparison.svg
            """
    }

    process get_GO_comparison {

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/figures/GO", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.trinotate_annotation_report.xls") from trinotate_ch

        output:
            tuple sample_id, file("*.svg"), file("*.pdf"), file("*.txt") into go_fig

        script:
            """
	        cp ${params.mypwd}/bin/GO_plots.R .

            cat ${sample_id}.trinotate_annotation_report.xls | awk 'FS="\\t",OFS="#" {print \$1,\$15,\$16,\$17}' | grep -v "gene_id" >all_GOs.txt

            while read line;do
                echo ${line} | cut -f 2,3,4 -d "#" | grep "GO:" | tr "#" "\\n" | tr "\\`" "\\n" | sed 's/\\. /,/g' | tr "," "\\n" | grep "GO:" | sort -u >>final_GOs.txt
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
            """
    }

    process summary_custom_uniprot {

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/figures/CustomUniProt", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.trinotate_annotation_report.xls") from custom_uniprot_ch

        output:
            tuple sample_id, file("${sample_id}_custom_uniprot_hits.txt") into custom_uniprot_sum
            tuple sample_id, file("${sample_id}_custom_uniprot_hits.svg"), file("${sample_id}_custom_uniprot_hits.pdf") into custom_uniprot_fig

        script:
            """
            #get custom blast hits
            cat ${sample_id}.trinotate_annotation_report.xls | cut -f 8 | grep [A-Z] | grep "|" | tr "\\`" "\n" | \
                cut -f 1 -d "^" | cut -f 3 -d "|" | cut -f 2 -d "_" >a.txt

            cat ${sample_id}.trinotate_annotation_report.xls | cut -f 9 | grep [A-Z] | grep "|" | tr "\\`" "\n" | \
                cut -f 1 -d "^" | cut -f 3 -d "|" | cut -f 2 -d "_" >b.txt

            cat a.txt b.txt | sort | uniq -c | sort -nr | head -n 20 | awk 'OFS="," {print \$1,\$2}' >${sample_id}_custom_uniprot_hits.txt

            rm a.txt b.txt

            cp ${params.mypwd}/conf/uni_tax.txt .
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

            cp ${params.mypwd}/bin/custom_uniprot_hits.R .
            Rscript custom_uniprot_hits.R ${sample_id}
            """
    }

} else if (params.onlyEvi) {

        if (params.nameEvi == "") {
            println("\n\t\033[0;31mMandatory argument \"nameEvi\" not specified. Argument needed for the \"onlyEvi\" option of TransPi. \n\033[0m")
            exit 0
        }
        //testing
        println("\n\tRunning Evidential Gene analysis only \n")

        Channel
            .fromPath("${params.workDir}/onlyEvi/*.{fa,fasta}", checkIfExists: true)
            .collect()
            .set{ evigene_ch_OE }

        process evigene_OE {

            label 'med_mem'

            tag "${params.nameEvi}"

            publishDir "${workDir}/${params.outdir}/evigene", mode: "copy", overwrite: true, pattern: "*.combined.okay.fa"

            input:
                file(assembly) from evigene_ch_OE

            output:
                file("${params.nameEvi}.combined.okay.fa") into ( evigene_ch_busco3_OE, evigene_ch_busco4_OE )

            script:
                def mem_MB=(task.memory.toMega())
                """
                echo -e "\\n-- Starting EviGene --\\n"

                cat ${assembly} >${params.nameEvi}.combined.fa

                $evi/scripts/prot/tr2aacds.pl -tidy -NCPU ${task.cpus} -MAXMEM ${mem_MB} -log -cdna ${params.nameEvi}.combined.fa

                echo -e "\\n-- DONE with EviGene --\\n"

                cp okayset/${params.nameEvi}.combined.okay.combined.fa ${params.nameEvi}.combined.okay.fa

                if [ -d tmpfiles/ ];then
                    rm -rf tmpfiles/
                fi
                """
        }

        process busco3_OE {

            label 'med_cpus'

            tag "${params.nameEvi}"

            publishDir "${workDir}/${params.outdir}/busco3", mode: "copy", overwrite: true

            input:
                file("${params.nameEvi}.combined.okay.fa") from evigene_ch_busco3_OE

            output:
                file("run_${params.nameEvi}.fa.bus") into busco3_ch_OE
                file("short_summary_${params.nameEvi}.fa.bus.txt") into busco3_summary_OE

            script:
                """
                echo -e "\\n-- Starting BUSCO --\\n"

                run_BUSCO.py -i ${params.nameEvi}.combined.okay.fa -o ${params.nameEvi}.fa.bus -l ${params.busco3db} -m tran -c ${task.cpus}

                echo -e "\\n-- DONE with BUSCO --\\n"

                cp run_${params.nameEvi}.fa.bus/short_summary_${params.nameEvi}.fa.bus.txt .
                """
        }

        process busco4_OE {

            conda "${params.cenv}"

            label 'med_cpus'

            tag "${params.nameEvi}"

            publishDir "${workDir}/${params.outdir}/busco4", mode: "copy", overwrite: true

            input:
                file("${params.nameEvi}.combined.okay.fa") from evigene_ch_busco4_OE

            output:
                file("${params.nameEvi}.fa.bus") into busco4_ch_OE
                file("*.${params.nameEvi}.fa.bus.txt") into ( busco4_summary_OE )

            script:
                """
                echo -e "\\n-- Starting BUSCO --\\n"

                busco -i ${params.nameEvi}.combined.okay.fa -o ${params.nameEvi}.fa.bus -l ${params.busco4db} -m tran -c ${task.cpus} --offline

                echo -e "\\n-- DONE with BUSCO --\\n"

                cp ${params.nameEvi}.fa.bus/short_summary.*.${params.nameEvi}.fa.bus.txt .
                """
        }

        process summary_busco_OE {

            tag "${params.nameEvi}"

            publishDir "${workDir}/${params.outdir}/stats", mode: "copy", overwrite: true

            input:
                file("*.${params.nameEvi}.fa.bus.txt") from busco4_summary_OE
                file("short_summary_${params.nameEvi}.fa.bus.txt") from busco3_summary_OE

            output:
                file("${params.nameEvi}.sum_busco.txt") into final_sum_2v3_OE

            script:
                """
                echo -e "Summary of BUSCO V4 \\n" >>${params.nameEvi}.sum_busco.txt
                echo "-- TransPi BUSCO V4 scores -- " >>${params.nameEvi}.sum_busco.txt
                cat *.${params.nameEvi}.fa.bus.txt >>${params.nameEvi}.sum_busco.txt
                echo -e "\\nSummary of BUSCO V3 \\n" >>${params.nameEvi}.sum_busco.txt
                echo "-- TransPi BUSCO V3 scores -- " >>${params.nameEvi}.sum_busco.txt
                cat short_summary_${params.nameEvi}.fa.bus.txt >>${params.nameEvi}.sum_busco.txt
                """
        }

} else {
    println("\n\t\033[0;31mMandatory argument not specified. For more info use `nextflow run TransPi.nf --help`\n\033[0m")
    exit 0
}

process get_run_info {

    publishDir "${workDir}/${params.outdir}/", mode: "copy", overwrite: true

    output:
       file("versions.txt") into run_info

    script:
        """
        echo "==========================================" >>versions.txt
        echo "TransPi - Transcriptome Analysis Pipeline" >>versions.txt
        echo -e "==========================================\\n" >>versions.txt
        echo -e "\\t\\tRUN INFO\\n" >>versions.txt
        echo "-- Kmers used --" >>versions.txt
        echo ${params.k} >>versions.txt

        echo -e "\\n-- Databases name and last update --" >>versions.txt

        v=\$( echo ${params.uniname} )
        echo "Uniprot_DB: \$v" >>versions.txt

        if [ -f ${params.mypwd}/DBs/uniprot_db/.lastrun.txt ];then
            v=\$( cat ${params.mypwd}/DBs/uniprot_db/.lastrun.txt )
        else
            v="No info available. Check Instructions on README."
        fi
        echo -e "Uniprot_DB last update: \$v \\n" >>versions.txt

        if [ -f ${params.mypwd}/DBs/hmmerdb/.lastrun.txt ];then
            v=\$( cat ${params.mypwd}/DBs/hmmerdb/.lastrun.txt )
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

        v=\$( echo "1.2" )
        echo "RNAmmer: \$v" >>versions.txt

        v=\$( TransDecoder.LongOrfs --version | cut -f 2 -d " " )
        echo "Transdecoder: \$v" >>versions.txt

        v=\$( run_BUSCO.py -v | cut -f 2 -d " " )
        echo "BUSCO3: \$v" >>versions.txt

        v=\$( echo "4.0.5" )
        echo "BUSCO4: \$v" >>versions.txt

        v=\$( cat ${params.mypwd}/transpi_env.yml | grep trinotate  | cut -f 2 -d "=" )
        echo "Trinotate: \$v" >>versions.txt

        v=\$( echo "4.1" )
        echo "SignalP: \$v" >>versions.txt

        v=\$( echo "2.0" )
        echo "tmhmm: \$v" >>versions.txt

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

workflow.onComplete {
    log.info ( workflow.success ? \
        "---------------------------------------------------------------------------------" \
        + "\n\033[0;32mDone! Open the following report in your browser --> ${params.workDir}/${params.outdir}/${params.tracedir}/transpi_report.html\033[0m" : \
        "---------------------------------------------------------------------------------" \
        + "\n\033[0;31mSomething went wrong. Check error message below and/or log files.\033[0m" )
}
