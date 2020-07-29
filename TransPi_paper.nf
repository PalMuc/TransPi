#!/usr/bin/env nextflow
/*
========================================================================================
                                    TransPi PAPER TEST
========================================================================================
                       Transcriptomes Analysis Pipeline
                       Author: Ramon E. Rivera-Vicens
                       Version: 1.0 (dev)
----------------------------------------------------------------------------------------
*/

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
} else if (params.onlyEvi) {
    println("\n\tRunning TransPi with your dataset\n")
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
                tuple sample_id, file("*.csv") into fastp_csv_OAS

            script:
                """
                echo ${sample_id}

                fastp -i ${reads[0]} -I ${reads[1]} -o left-${sample_id}.filter.fq -O right-${sample_id}.filter.fq --detect_adapter_for_pe \
                --average_qual ${params.minQual} --overrepresentation_analysis --html ${sample_id}.fastp.html --json ${sample_id}.fastp.json \
                --thread ${task.cpus} --report_title ${sample_id}

                bash get_readstats.sh ${sample_id}.fastp.json
                bash get_readqual.sh ${sample_id}.fastp.json
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
            tuple sample_id, file("left-${sample_id}.norm.fq"), file("right-${sample_id}.norm.fq") into ( mapping_reads_trinity_OAS, mapping_reads_evi_OAS )

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
            tuple sample_id, file("${sample_id}.Trinity.fa") into ( assemblies_ch_trinity_OAS, busco3_ch_trinity_OAS, busco4_ch_trinity_OAS, mapping_trinity_OAS )

        script:
            """
            mem=\$( echo ${task.memory} | cut -f 1 -d " " )

            Trinity --max_memory \${mem}G --seqType fq --left left-${sample_id}.norm.fq --right right-${sample_id}.norm.fq --CPU ${task.cpus} --no_normalize_reads --full_cleanup --output trinity_out_dir

            mv trinity_out_dir.Trinity.fasta ${sample_id}.Trinity.fa
            """
    }

    mapping_trinity=Channel.create()
    mapping_trinity_OAS.mix( mapping_reads_trinity_OAS ).groupTuple(by:0,size:2).into(mapping_trinity)

    process mapping_trinity_assembly_OAS {

        label 'big_cpus'

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/mapping", mode: "copy", overwrite: true

        input:
            tuple sample_id, file(files), file(files2) from mapping_trinity

        output:
            tuple sample_id, file("log*") into mapping_results1

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

    process soap_assembly_OAS {

        label 'med_mem'

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/assemblies", mode: "copy", overwrite: true

        input:
            val k from "${params.k}"
            tuple sample_id, file("left-${sample_id}.norm.fq"), file("right-${sample_id}.norm.fq") from norm_reads_soap_OAS

        output:
            tuple sample_id, file("${sample_id}.SOAP.fa") into assemblies_ch_soap_OAS
            tuple sample_id, file("${sample_id}.SOAP.k*") into assemblies_ch_soap_busco3_OAS, assemblies_ch_soap_busco4_OAS

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

            for x in `echo $k | tr "," " "`;do
                cp output\${x}.scafSeq ${sample_id}.SOAP.k\${x}.fa
            done

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
            tuple sample_id, file("${sample_id}.Velvet.k*") into assemblies_ch_velvet_busco3_OAS, assemblies_ch_velvet_busco4_OAS

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

            for x in `echo $k | tr "," " "`;do
                cp oases.\${x}/contigs.fa ${sample_id}.Velvet.k\${x}.fa
            done

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
            tuple sample_id, file("${sample_id}.SPADES.k*") into assemblies_ch_spades_busco3_OAS, assemblies_ch_spades_busco4_OAS

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

            for x in `echo $k | tr "," " "`;do
                cp ${sample_id}_spades_\${x}/transcripts.fasta ${sample_id}.SPADES.k\${x}.fa
            done

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
            tuple sample_id, file("${sample_id}.TransABySS.k*") into assemblies_ch_transabyss_busco3_OAS, assemblies_ch_transabyss_busco4_OAS

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

            for x in `echo $k | tr "," " "`;do
                cp ${sample_id}_transabyss_\${x}/k\${x}.transabyss.fa-final.fa ${sample_id}.TransABySS.k\${x}.fa
            done

            rm -rf ${sample_id}_transabyss_*
            """
    }

    all_assemblies = Channel.create()
    assemblies_ch_trinity_OAS.mix( assemblies_ch_transabyss_OAS, assemblies_ch_spades_OAS, assemblies_ch_velvet_OAS, assemblies_ch_soap_OAS ).groupTuple(by:0,size:5).into(all_assemblies)

    process evigene_OAS {

        label 'med_mem'

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/evigene", mode: "copy", overwrite: true

        input:
            tuple sample_id, file(assemblies) from all_assemblies

        output:
            tuple sample_id, file("*.combined.okay.fa") into ( evigene_ch_busco3_OAS, evigene_ch_busco4_OAS, mapping_evi_OAS )
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

    mapping_evi=Channel.create()
    mapping_evi_OAS.mix( mapping_reads_evi_OAS ).groupTuple(by:0,size:2).view().into(mapping_evi)

    process mapping_evi_OAS {

        label 'big_cpus'

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/mapping", mode: "copy", overwrite: true

        input:
            tuple sample_id, file(files), file(files2) from mapping_evi

        output:
            tuple sample_id, file("log*") into mapping_results2

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

    process summary_evigene_individual_OAS {

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/stats", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.combined.fa"), file("${sample_id}.combined.okay.fa") from evigene_summary_OAS

        output:
            tuple sample_id, file("${sample_id}.sum_preEG.txt"), file("${sample_id}.sum_EG.txt") into final_sum_1_OAS
            tuple sample_id, file("${sample_id}.sum_preEG.csv"), file("${sample_id}.sum_EG.csv") into summary_evi_csv_OAS

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

            # csv report
            echo "Total,Trinity,SOAP,Velvet,SPADES,TransBySS" >${sample_id}.sum_preEG.csv
            total=\$( cat ${sample_id}.combined.fa | grep -c ">" )
            trinity=\$( cat ${sample_id}.combined.fa | grep -c ">TRINITY" )
            soap=\$( cat ${sample_id}.combined.fa | grep -c ">SOAP" )
            velvet=\$( cat ${sample_id}.combined.fa | grep -c ">Velvet" )
            spades=\$( cat ${sample_id}.combined.fa | grep -c ">SPADES" )
            transabyss=\$( cat ${sample_id}.combined.fa | grep -c ">TransABySS" )
            echo "\${total},\${trinity},\${soap},\${velvet},\${spades},\${transabyss}" >>${sample_id}.sum_preEG.csv

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

            # csv report after evigene
            echo "Total,Trinity,SOAP,Velvet,SPADES,TransBySS" >${sample_id}.sum_EG.csv
            total=\$( cat ${sample_id}.combined.okay.fa | grep -c ">" )
            trinity=\$( cat ${sample_id}.combined.okay.fa | grep -c ">TRINITY" )
            soap=\$( cat ${sample_id}.combined.okay.fa | grep -c ">SOAP" )
            velvet=\$( cat ${sample_id}.combined.okay.fa | grep -c ">Velvet" )
            spades=\$( cat ${sample_id}.combined.okay.fa | grep -c ">SPADES" )
            transabyss=\$( cat ${sample_id}.combined.okay.fa | grep -c ">TransABySS" )
            echo "\${total},\${trinity},\${soap},\${velvet},\${spades},\${transabyss}" >>${sample_id}.sum_EG.csv
            """

    }

    busco3_all_OAS = Channel.create()
    assemblies_ch_soap_busco3_OAS.mix( assemblies_ch_velvet_busco3_OAS, assemblies_ch_spades_busco3_OAS, assemblies_ch_transabyss_busco3_OAS ).into(busco3_all_OAS)

    process busco3_all_OAS {

        label 'med_cpus'

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/busco3_all", mode: "copy", overwrite: true

        input:
            tuple sample_id, file(files) from busco3_all_OAS

        output:
            tuple sample_id, file("*.bus3") into busco3_all_ch

        script:
            """
            echo "$files" | tr " " "\\n" | awk -F '\\.fa' '{print \$1}' >list.txt

            for x in `cat list.txt`;do

                echo -e "\\n-- Starting BUSCO --\\n"

                run_BUSCO.py -i \${x}.fa -o \${x}.bus3 -l ${params.busco3db} -m tran -c ${task.cpus}

                echo -e "\\n-- DONE with BUSCO --\\n"

            done
            """
    }

    process busco3_OAS {

        label 'med_cpus'

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/busco3", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.combined.okay.fa") from evigene_ch_busco3_OAS

        output:
            tuple sample_id, file("run_${sample_id}.TransPi.bus") into busco3_ch
            tuple sample_id, file("*${sample_id}.TransPi.bus.txt") into ( busco3_summary_OAS, busco3_comp_1_OAS )

        script:
            """
            echo -e "\\n-- Starting BUSCO --\\n"

            run_BUSCO.py -i ${sample_id}.combined.okay.fa -o ${sample_id}.TransPi.bus -l ${params.busco3db} -m tran -c ${task.cpus}

            echo -e "\\n-- DONE with BUSCO --\\n"

            cp run_${sample_id}.TransPi.bus/short_summary_${sample_id}.TransPi.bus.txt .
            """
    }

    process busco3_tri_OAS {

        label 'med_cpus'

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/busco3", mode: "copy", overwrite: true

        input:
            tuple sample_id, file("${sample_id}.Trinity.fa") from busco3_ch_trinity_OAS

        output:
            tuple sample_id, file("*${sample_id}.Trinity.bus.txt") into ( busco3_ch_trinity_sum_OAS, busco3_comp_2_OAS )

        script:
            """
            echo -e "\\n-- Starting BUSCO --\\n"

            run_BUSCO.py -i ${sample_id}.Trinity.fa -o ${sample_id}.Trinity.bus -l ${params.busco3db} -m tran -c ${task.cpus}

            echo -e "\\n-- DONE with BUSCO --\\n"

            cp run_${sample_id}.Trinity.bus/short_summary_${sample_id}.Trinity.bus.txt .
            """
    }

    busco4_all_OAS = Channel.create()
    assemblies_ch_soap_busco4_OAS.mix( assemblies_ch_velvet_busco4_OAS, assemblies_ch_spades_busco4_OAS, assemblies_ch_transabyss_busco4_OAS ).into(busco4_all_OAS)

    process busco4_all_OAS {

        conda "${params.cenv}"

        label 'med_cpus'

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/busco4_all", mode: "copy", overwrite: true

        input:
            tuple sample_id, file(files) from busco4_all_OAS

        output:
            tuple sample_id, file("*.bus4") into busco4_all_ch

        script:
            """
            echo "$files" | tr " " "\\n" | awk -F '\\.fa' '{print \$1}' >list.txt

            for x in `cat list.txt`;do

                echo -e "\\n-- Starting BUSCO --\\n"

                busco -i \${x}.fa -o \${x}.bus4 -l ${params.busco4db} -m tran -c ${task.cpus} --offline

                echo -e "\\n-- DONE with BUSCO --\\n"

            done
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
            tuple sample_id, file("${sample_id}.TransPi.bus") into busco4_ch
            tuple sample_id, file("*${sample_id}.TransPi.bus.txt") into ( busco4_summary_OAS, busco4_comp_1_OAS )

        script:
            """
            echo -e "\\n-- Starting BUSCO --\\n"

            busco -i ${sample_id}.combined.okay.fa -o ${sample_id}.TransPi.bus -l ${params.busco4db} -m tran -c ${task.cpus} --offline

            echo -e "\\n-- DONE with BUSCO --\\n"

            cp ${sample_id}.TransPi.bus/short_summary.*.${sample_id}.TransPi.bus.txt .
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
            tuple sample_id, file("*${sample_id}.Trinity.bus.txt") into ( busco4_ch_trinity_sum_OAS, busco4_comp_2_OAS )

        script:
            """
            echo -e "\\n-- Starting BUSCO --\\n"

            busco -i ${sample_id}.Trinity.fa -o ${sample_id}.Trinity.bus -l ${params.busco4db} -m tran -c ${task.cpus} --offline

            echo -e "\\n-- DONE with BUSCO --\\n"

            cp ${sample_id}.Trinity.bus/short_summary.*.${sample_id}.Trinity.bus.txt .
            """
    }

    busco3_sum_OAS = Channel.create()
    busco3_summary_OAS.mix(busco3_ch_trinity_sum_OAS).groupTuple(by:0,size:2).into(busco3_sum_OAS)

    process summary_busco3_individual_OAS {

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/stats", mode: "copy", overwrite: true

        input:
            tuple sample_id, file(files) from busco3_sum_OAS

        output:
            tuple sample_id, file("${sample_id}.sum_busco3.txt") into final_sum_2v3_OAS

        script:
            """
            tri=\$( echo $files | tr " " "\\n" | grep ".Trinity.bus.txt" )
            trans=\$( echo $files | tr " " "\\n" | grep ".TransPi.bus.txt" )
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

        publishDir "${workDir}/${params.outdir}/stats", mode: "copy", overwrite: true

        input:
            tuple sample_id, file(files) from busco4_sum_OAS

        output:
            tuple sample_id, file("${sample_id}.sum_busco4.txt") into final_sum_2v4_OAS

        script:
            """
            tri=\$( echo $files | tr " " "\\n" | grep ".Trinity.bus.txt" )
            trans=\$( echo $files | tr " " "\\n" | grep ".TransPi.bus.txt" )
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

        publishDir "${workDir}/${params.outdir}/figures/BUSCO3", mode: "copy", overwrite: true

        input:
            tuple sample_id, file(files) from busco3_comp_OAS

        output:
            tuple sample_id, file("${sample_id}_BUSCO3_comparison.pdf"), file("${sample_id}_BUSCO3_comparison.svg") into busco3_fig_OAS
            tuple sample_id, file("*.csv") into busco3_OAS_csv

        script:
            """
            set +e
            tri=\$( echo $files | tr " " "\\n" | grep ".Trinity.bus.txt" )
            trans=\$( echo $files | tr " " "\\n" | grep ".TransPi.bus.txt" )
            bash get_busco_val.sh \${tri} \${trans} v3 ${sample_id}
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
            # csv
            sed -i 's/\$/\n/g' final_*
            cat final_spec final_perc final_num | tr -d "'" >${sample_id}_busco3.csv
            """
    }

    busco4_comp_OAS = Channel.create()
    busco4_comp_1_OAS.mix(busco4_comp_2_OAS).groupTuple(by:0,size:2).into(busco4_comp_OAS)

    process get_busco4_comparison_OAS {

        tag "${sample_id}"

        publishDir "${workDir}/${params.outdir}/figures/BUSCO4", mode: "copy", overwrite: true

        input:
            tuple sample_id, file(files) from busco4_comp_OAS

        output:
            tuple sample_id, file("${sample_id}_BUSCO4_comparison.pdf"), file("${sample_id}_BUSCO4_comparison.svg") into busco4_fig_OAS
            tuple sample_id, file("*.csv") into busco4_OAS_csv

        script:
            """
            set +e
            tri=\$( echo $files | tr " " "\\n" | grep ".Trinity.bus.txt" )
            trans=\$( echo $files | tr " " "\\n" | grep ".TransPi.bus.txt" )
            bash get_busco_val.sh \${tri} \${trans} v4 ${sample_id}
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
            # csv
            sed -i 's/\$/\n/g' final_*
            cat final_spec final_perc final_num | tr -d "'" >${sample_id}_busco4.csv
            """
    }

}
