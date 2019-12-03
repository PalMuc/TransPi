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

process custom_diamond_db {
    script:
        """
        set +u
        source ${params.condash}/etc/profile.d/conda.sh
        conda activate TransPi

        cd ${params.mypwd}
        echo -e "-- Checking if Diamond database folder is present --\n"
        if [ ! -d diamonddb_custom/ ];then
            echo -e "-- Folder is not present, creating one and the Diamond database --\n"
            mkdir diamonddb_custom/
            cd diamonddb_custom
            cp ${params.uniprot} .
            diamond makedb --in ${params.uniname} -d ${params.uniname}
            export unidb=`pwd`/${params.uniname}
            cd ../
        elif [ -d diamonddb_custom/ ];then
            echo -e "-- Folder is present. Checking if Diamond database is built --\n"
            cd diamonddb_custom
            if [ ! -e ${params.uniname}.dmnd ];then
                echo -e "-- Diamond database not present, creating one --\n"
                cp ${params.uniprot} .
                diamond makedb --in ${params.uniname} -d ${params.uniname}
                export unidb=`pwd`/${params.uniname}
            elif [ -e ${params.uniname}.dmnd  ];then
                echo -e "-- Diamond database already created --\n"
                export unidb=`pwd`/${params.uniname}
            fi
            cd ../
        fi
        """
}

process hmmer_db {
    script:
        """
        set +u
        source ${params.condash}/etc/profile.d/conda.sh
        conda activate TransPi

        cd ${params.mypwd}
        echo -e "-- Checking if HMMER database folder is present --\n"
        if [ -d hmmerdb/ ];then
    	    echo -e "-- Folder is present. Checking if HMMER database is built --\n"
    		cd hmmerdb
    		if [ ! -e ${params.pfname}.h3f ] && [ ! -e ${params.pfname}.h3i ] && [ ! -e ${params.pfname}.h3m ] && [ ! -e ${params.pfname}.h3p ];then
                echo -e "-- HMMER database not present, creating one --\n"
        	    hmmpress ${params.pfname}
        	    export pf=`pwd`/${params.pfname}
    	    elif [ -s ${params.pfname}.h3f ] && [ -s ${params.pfname}.h3i ] && [ -s ${params.pfname}.h3m ] && [ -s ${params.pfname}.h3p ];then
                echo -e "-- HMMER database already created --\n"
        	    export pf=`pwd`/${params.pfname}
    	    fi
            cd ../
        fi
        """
}

process swiss_diamond_db {
    script:
        """
        set +u
        source ${params.condash}/etc/profile.d/conda.sh
        conda activate TransPi

        cd ${params.mypwd}/sqlite_db
        if [ -e uniprot_sprot.pep ];then
            cd ${params.mypwd}
            if [ ! -d diamonddb_swiss/ ];then
                echo -e "-- Folder is not present, creating one and the Diamond database --\n"
                mkdir diamonddb_swiss
                cd diamonddb_swiss
                cp ${params.mypwd}/sqlite_db/uniprot_sprot.pep .
                diamond makedb --in uniprot_sprot.pep -d uniprot_sprot.pep
                export swissdb=`pwd`/uniprot_sprot.pep
            elif [ -d diamonddb_swiss/ ];then
                cd diamonddb_swiss
                if [ ! -e uniprot_sprot.pep.dmnd ];then
                    echo -e "-- Diamond database not present, creating one --\n"
                    cp ${params.mypwd}/sqlite_db/uniprot_sprot.pep .
                    diamond makedb --in uniprot_sprot.pep -d uniprot_sprot.pep
                    export swissdb=`pwd`/uniprot_sprot.pep
                elif [ -e uniprot_sprot.pep.dmnd ];then
                    echo -e "-- Diamond database already created --\n"
                    export swissdb=`pwd`/uniprot_sprot.pep
                fi
            fi
        elif [ ! -e uniprot_sprot.pep ];then
            cd ${params.mypwd}
            if [ ! -d diamonddb_swiss/ ];then
                echo -e "-- Folder is not present, creating one and the Diamond database --\n"
                mkdir diamonddb_swiss
                cd diamonddb_swiss
                wget http://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
                EMBL_swissprot_parser.pl uniprot_sprot.dat.gz ind
                rm ind.*
                mv uniprot_sprot.dat.gz.pep uniprot_sprot.pep
                diamond makedb --in uniprot_sprot.pep -d uniprot_sprot.pep
                export swissdb=`pwd`/uniprot_sprot.pep
            elif [ -d diamonddb_swiss/ ];then
                echo -e "-- Folder is present. Checking if Diamond database is built --\n"
                cd diamonddb_swiss
                if [ ! -e uniprot_sprot.pep.dmnd ];then
                    if [ ! -e uniprot_sprot.pep ];then
                        echo -e "-- Diamond database not present, creating one --\n"
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
                    echo -e "-- Diamond database already created --\n"
                    export swissdb=`pwd`/uniprot_sprot.pep
                fi
            fi
        fi
        """
}

process results_dir {
    script:
        """
        cd ${params.mypwd}
        if [ ! -d results/ ];then
            mkdir results
        fi
        """
}

reads_ch=Channel.fromFilePairs("${params.mypwd}/reads/*_R{1,2}.fastq.gz")

process normalize_reads {

    label 'med_mem'

    tag "${sample_id}"

    input:
        set sample_id, input_read from reads_ch

    output:
        set sample_id, file("left-${sample_id}.norm.fq"), file("right-${sample_id}.norm.fq") into ( norm_reads_soap, norm_reads_velvet, norm_reads_trinity, norm_reads_idba )

    script:
        //def mem=(task.memory)
        //def mem_MB=(task.memory.toMega())
        """
        set +u
        source ${params.condash}/etc/profile.d/conda.sh
        conda activate TransPi

        echo ${sample_id}
        zcat ${input_read[0]} >left-${sample_id}.fq &
        zcat ${input_read[1]} >right-${sample_id}.fq

        echo -e "\n-- Starting Normalization --\n"

        mem=\$( echo ${task.memory} | cut -f 1 -d " " )

        ${params.in_norm} --seqType fq -JM \${mem}G --max_cov 100 --min_cov 1 --left left-${sample_id}.fq --right right-${sample_id}.fq --pairs_together --PARALLEL_STATS --CPU ${task.cpus}

        echo -e "\n-- DONE with Normalization --\n"

        mv left.norm.fq left-"${sample_id}".norm.fq
        mv right.norm.fq right-"${sample_id}".norm.fq

        rm left-${sample_id}.fq right-${sample_id}.fq
        """
}

process trinity_assembly {

    label 'med_mem'

    tag "${sample_id}"

    input:
        set sample_id, file("left-${sample_id}.norm.fq"), file("right-${sample_id}.norm.fq") from norm_reads_trinity

    output:
        set sample_id, file("${sample_id}.Trinity.fa") into assemblies_ch_trinity

    script:
        """
        set +u
        source ${params.condash}/etc/profile.d/conda.sh
        conda activate TransPi

        mem=\$( echo ${task.memory} | cut -f 1 -d " " )

        ${params.tr} --max_memory \${mem}G --seqType fq --left left-${sample_id}.norm.fq --right right-${sample_id}.norm.fq --CPU ${task.cpus} --no_normalize_reads --full_cleanup --output trinity_out_dir

        mv trinity_out_dir.Trinity.fasta ${sample_id}.Trinity.fa

        cp ${sample_id}.Trinity.fa ${params.mypwd}/results/${sample_id}.Trinity.fa
        """
}

process soap_assembly {

    label 'med_mem'

    tag "${sample_id}"

    input:
        val k from "${params.k}"
        set sample_id, file("left-${sample_id}.norm.fq"), file("right-${sample_id}.norm.fq") from norm_reads_soap

    output:
        set sample_id, file("${sample_id}.SOAP.fa") into assemblies_ch_soap

    script:
        """
        set +u
        source ${params.condash}/etc/profile.d/conda.sh
        conda activate TransPi

        echo -e "\n-- Generating SOAP config file --\n"
        echo "max_rd_len="${params.max_rd_len} >>config.txt
        echo "[LIB]" >>config.txt
        echo "rd_len_cutof="${params.rd_len_cutof} >>config.txt
        echo "avg_ins="${params.avg_ins} >>config.txt
        echo "reverse_seq="${params.reverse_seq} >>config.txt
        echo "asm_flags="${params.asm_flags} >>config.txt
        echo "map_len="${params.map_len} >>config.txt
        echo "q1="left-${sample_id}.norm.fq >>config.txt
        echo "q2="right-${sample_id}.norm.fq >>config.txt

        echo -e "\n-- Starting SOAP assemblies --\n"

        for x in `echo $k | tr "," " "`;do
            echo -e "\n-- SOAP k\${x} --\n"
            ${params.soap} all -s config.txt -K \${x} -o output\${x} -p ${task.cpus}
            sed -i "s/>/>SOAP.k\${x}./g" output\${x}.scafSeq
        done

        echo -e "\n-- Finished with the assemblies --\n"

        cat output*.scafSeq >${sample_id}.SOAP.fa

        cp ${sample_id}.SOAP.fa ${params.mypwd}/results/${sample_id}.SOAP.fa

        rm -rf output*
        """
}

process velvet_oases_assembly {

    label 'med_mem'

    tag "${sample_id}"

    input:
        val k from "${params.k}"
        set sample_id, file("left-${sample_id}.norm.fq"), file("right-${sample_id}.norm.fq") from norm_reads_velvet

    output:
        set sample_id, file("${sample_id}.Velvet.fa") into assemblies_ch_velvet

    script:
        """
	    set +u
        source ${params.condash}/etc/profile.d/conda.sh
        conda activate TransPi

	    echo -e "\n-- Starting with Velveth --\n"
        for x in `echo $k | tr "," " "`;do
            echo -e "\n-- k\${x} --\n"
            ${params.vh} oases.\${x} \${x} -shortPaired -fastq -separate left-${sample_id}.norm.fq right-${sample_id}.norm.fq
        done

        echo -e "\n-- Starting with Velvetg --\n"
        for x in `echo $k | tr "," " "`;do
            echo -e "\n-- vg \${x} --\n"
            ${params.vg} oases.\${x} -read_trkg yes
        done

        echo -e "\n-- Starting with Oases --\n"
        for x in `echo $k | tr "," " "`;do
            echo -e "\n-- oases \${x} --\n"
            ${params.oa} oases.\${x}
        done

        echo -e "\n-- Finished with Velvet/Oases assemblies --\n"

        for x in `echo $k | tr "," " "`;do
            sed -i "s/>/>Velvet.k\${x}./g" oases.\${x}/contigs.fa
        done

        cat oases.*/contigs.fa >${sample_id}.Velvet.fa

        cp ${sample_id}.Velvet.fa ${params.mypwd}/results/${sample_id}.Velvet.fa

        rm -rf oases.*
        """
}

process idba_assembly {

    label 'med_mem'

    tag "${sample_id}"

    input:
        val k from "${params.k}"
        set sample_id, file("left-${sample_id}.norm.fq"), file("right-${sample_id}.norm.fq") from norm_reads_idba

    output:
        set sample_id, file("${sample_id}.IDBA.fa") into assemblies_ch_idba

    script:
        """
        set +u
        source ${params.condash}/etc/profile.d/conda.sh
        conda activate TransPi

        echo -e "\n-- Starting IDBA assemblies --\n"

        echo -e "\n-- Converting reads for IDBA --\n"
        fq2fa --merge left-${sample_id}.norm.fq right-${sample_id}.norm.fq ${sample_id}_reads.fa

        for x in `echo $k | tr "," " "`;do
            echo -e "\n-- IDBA k\${x} --\n"
            idba_tran --out ${sample_id}_idba_\${x} --read ${sample_id}_reads.fa --num_threads ${task.cpus} --mink \${x} --maxk \${x}
        done

        echo -e "\n-- Finished with the assemblies --\n"

        for x in `echo $k | tr "," " "`;do
            sed -i "s/>/>IDBA.k/g" ${sample_id}_idba_\${x}/contig.fa && sed -i "s/contig-//g" ${sample_id}_idba_\${x}/contig.fa
        done

        cat ${sample_id}_idba_*/contig.fa >${sample_id}.IDBA.fa

        cp ${sample_id}.IDBA.fa ${params.mypwd}/results/${sample_id}.IDBA.fa

        rm ${sample_id}_reads.fa

        rm -rf ${sample_id}_idba_*
        """
}

process evigene {

    label 'med_mem'

    tag "${sample_id}"

    publishDir "${params.mypwd}/evigene", mode: "copy", overwrite: true

    input:
        set sample_id, file("${sample_id}.Trinity.fa") from assemblies_ch_trinity
        set sample_id, file("${sample_id}.SOAP.fa") from assemblies_ch_soap
        set sample_id, file("${sample_id}.Velvet.fa") from assemblies_ch_velvet
        set sample_id, file("${sample_id}.IDBA.fa") from assemblies_ch_idba

    output:
        set sample_id, file("${sample_id}.combined.okay.fa"), file("${sample_id}.combined.okay.cds"), file("${sample_id}.combined.okay.aa") into ( evigene_ch_busco, evigene_ch_transdecoder, evigene_ch_transdecoder_short, evigene_ch_diamond, evigene_ch_rnammer, evigene_ch_trinotate, evigene_ch_trinotate_custom )
        // for test
        set sample_id, file("${sample_id}.combined.okay.fa"), file("${sample_id}.combined.okay.aa") into ( evigene_ch_diamond_evigene, evigene_ch_diamond_evigene_custom, evigene_ch_hmmer_trinotate_evigene, evigene_ch_signalP_trinotate_evigene, evigene_ch_tmhmm_trinotate_evigene, evigene_ch_trinotate_evigene_aa )
        set sample_id, file("${sample_id}.combined.okay.fa") into ( evigene_ch_rnammer_trinotate_evigene, evigene_ch_trinotate_evigene )
        // end test
        set sample_id, file("${sample_id}.combined.fa"), file("${sample_id}.combined.okay.cds"), file("${sample_id}.combined.okay.aa") into evigene_summary_ind
        set sample_id, file("${sample_id}.combined.fa"), file("${sample_id}.combined.okay.cds"), file("${sample_id}.combined.okay.aa") into evigene_summary_cds_aa

    script:
        def mem_MB=(task.memory.toMega())

        """
        echo -e "\n-- Starting EviGene --\n"

        cat *.fa >${sample_id}.combined.fa

        tr2aacds.pl -tidy -NCPU ${task.cpus} -MAXMEM ${mem_MB} -log -cdna ${sample_id}.combined.fa

        echo -e "\n-- DONE with EviGene --\n"

        cp okayset/${sample_id}.combined.okay.combined.fa ${sample_id}.combined.okay.fa
        cp okayset/${sample_id}.combined.okay.combined.cds ${sample_id}.combined.okay.cds
        cp okayset/${sample_id}.combined.okay.combined.aa ${sample_id}.combined.okay.aa

	    if [ -d tmpfiles/ ];then
	        rm -rf tmpfiles/
	    fi
	    """
}

process busco_evigene {

    label 'big_cpus'

    tag "${sample_id}"

    publishDir "${params.mypwd}/busco_evigene", mode: "copy", overwrite: true

    input:
        set sample_id, file("${sample_id}.combined.okay.fa"), file("${sample_id}.combined.okay.cds"), file("${sample_id}.combined.okay.aa") from evigene_ch_busco

    output:
        set sample_id, file("run_${sample_id}.fa.bus"), file("run_${sample_id}.cds.bus"), file("run_${sample_id}.aa.bus") into busco_ch
        set sample_id, file("short_summary_${sample_id}.fa.bus.txt"), file("short_summary_${sample_id}.cds.bus.txt"), file("short_summary_${sample_id}.aa.bus.txt") into busco_summary

    script:
        """
        set +u
        source ${params.condash}/etc/profile.d/conda.sh
        conda activate TransPi

        echo -e "\n-- Starting with BUSCO --\n"

        run_BUSCO.py -i ${sample_id}.combined.okay.fa -o ${sample_id}.fa.bus -l ${params.buscodb} -m tran -c ${task.cpus}

        run_BUSCO.py -i ${sample_id}.combined.okay.cds -o ${sample_id}.cds.bus -l ${params.buscodb} -m tran -c ${task.cpus}

        ## for test
        run_BUSCO.py -i ${sample_id}.combined.okay.aa -o ${sample_id}.aa.bus -l ${params.buscodb} -m prot -c ${task.cpus}

        echo -e "\n-- DONE with BUSCO --\n"

        cp run_${sample_id}.fa.bus/short_summary_${sample_id}.fa.bus.txt ${params.mypwd}/results/short_summary_${sample_id}.fa.bus.txt
        cp run_${sample_id}.cds.bus/short_summary_${sample_id}.cds.bus.txt ${params.mypwd}/results/short_summary_${sample_id}.cds.bus.txt

        ## for test
        cp run_${sample_id}.pep.bus/short_summary_${sample_id}.aa.bus.txt ${params.mypwd}/results/short_summary_${sample_id}.aa.bus.txt

        cp run_${sample_id}.fa.bus/short_summary_${sample_id}.fa.bus.txt .
        cp run_${sample_id}.cds.bus/short_summary_${sample_id}.cds.bus.txt .

        ## for test
        cp run_${sample_id}.pep.bus/short_summary_${sample_id}.aa.bus.txt .
        """
}

process transdecoder_long {

    label 'low_cpus'

    tag "${sample_id}"

    input:
        set sample_id, file("${sample_id}.combined.okay.fa") from evigene_ch_transdecoder

    output:
        set sample_id, file("${sample_id}.combined.okay.fa.transdecoder.pep") into transdecoder_ch_diamond
        set sample_id, file("${sample_id}.combined.okay.fa.transdecoder.pep") into transdecoder_ch_hmmer
        set sample_id, file("${sample_id}.combined.okay.fa.transdecoder.pep") into transdecoder_ch_signalp
        set sample_id, file("${sample_id}.combined.okay.fa.transdecoder.pep") into transdecoder_ch_tmhmm
        set sample_id, file("${sample_id}.combined.okay.fa.transdecoder.pep") into transdecoder_ch_trinotate
        set sample_id, file("${sample_id}.combined.okay.fa.transdecoder.pep") into transdecoder_ch_diamond_custom
        set sample_id, file("${sample_id}.transdecoder.stats") into transdecoder_summary
        // for test
        set sample_id, file("${sample_id}.combined.okay.fa.transdecoder.pep") into transdecoder_ch_busco

    script:
        """
        set +u
        source ${params.condash}/etc/profile.d/conda.sh
        conda activate TransPi

        unidb=${params.mypwd}/diamonddb_custom/${params.uniname}
        pf=${params.mypwd}/hmmerdb/${params.pfname}

        echo -e "\n-- TransDecoder.LongOrfs... --\n"

        ${params.torf} -t ${sample_id}.combined.okay.fa

        echo -e "\n-- Done with TransDecoder.LongOrfs --\n"

        fname=${sample_id}.combined.okay.fa

        echo -e "\n-- Starting Diamond (blastp) --\n"

        ${params.diam} blastp -d \$unidb -q \$fname.transdecoder_dir/longest_orfs.pep -p ${task.cpus} -f 6 -k 1 -e 0.00001 >diamond_blastp.outfmt6

        echo -e "\n-- Done with Diamond (blastp) --\n"

        echo -e "\n-- Starting HMMER --\n"

        ${params.hmsc} --cpu ${task.cpus} --domtblout pfam.domtblout \$pf \$fname.transdecoder_dir/longest_orfs.pep

        echo -e "\n-- Done with HMMER --\n"

        echo -e "\n-- TransDecoder.Predict... --\n"

        ${params.tpred} -t ${sample_id}.combined.okay.fa --retain_pfam_hits pfam.domtblout --retain_blastp_hits diamond_blastp.outfmt6

        echo -e "\n-- Done with TransDecoder.Predict --\n"

        echo -e "\n-- Calculating statistics... --\n"

        #Calculate statistics of Transdecoder
        echo "- Transdecoder stats for "${sample_id} >>${sample_id}.transdecoder.stats
        orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep -c ">" )
        echo "Total number of ORFs: "\$orfnum >>${sample_id}.transdecoder.stats
        echo -e "\t Of these ORFs" >>${sample_id}.transdecoder.stats
        orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep ">" | grep -c "|" )
        echo -e "\t\t with annotations: "\$orfnum >>${sample_id}.transdecoder.stats
        orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep ">" | grep -v "|" | grep -c ">" )
        echo -e "\t\t no annotation: "\$orfnum >>${sample_id}.transdecoder.stats
        echo >>${sample_id}.transdecoder.stats
        orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep -c "ORF type:complete" )
        echo -e "\t ORFs type=complete: "\$orfnum >>${sample_id}.transdecoder.stats
        orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep "ORF type:complete" | grep -c "|" )
        echo -e "\t\t with annotations: "\$orfnum >>${sample_id}.transdecoder.stats
        echo >>${sample_id}.transdecoder.stats
        orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep -c "ORF type:5prime_partial" )
        echo -e "\t ORFs type=5prime_partial: "\$orfnum >>${sample_id}.transdecoder.stats
        orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep "ORF type:5prime_partial" | grep -c "|" )
        echo -e "\t\t with annotations: "\$orfnum >>${sample_id}.transdecoder.stats
        echo >>${sample_id}.transdecoder.stats
        orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep -c "ORF type:3prime_partial" )
        echo -e "\t ORFs type=3prime_partial: "\$orfnum >>${sample_id}.transdecoder.stats
        orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep "ORF type:3prime_partial" | grep -c "|" )
        echo -e "\t\t with annotations: "\$orfnum >>${sample_id}.transdecoder.stats
        echo >>${sample_id}.transdecoder.stats
        orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep -c "ORF type:internal" )
        echo -e "\t ORFs type=internal: "\$orfnum >>${sample_id}.transdecoder.stats
        orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep "ORF type:internal" | grep -c "|" )
        echo -e "\t\t with annotations: "\$orfnum >>${sample_id}.transdecoder.stats
        orfnum=0

        echo -e "\n-- Done with statistics --\n"

        echo -e "\n-- DONE with TransDecoder --\n"

        ## added _long for test
        cp ${sample_id}.transdecoder.stats ${params.mypwd}/results/${sample_id}.transdecoder.stats_long
        """

}

process transdecoder_short {

    label 'low_cpus'

    tag "${sample_id}"

    input:
        set sample_id, file("${sample_id}.combined.okay.fa") from evigene_ch_transdecoder_short

    output:
        set sample_id, file("${sample_id}.combined.okay.fa.transdecoder.pep") into transdecoder_ch_diamond_short
        set sample_id, file("${sample_id}.combined.okay.fa.transdecoder.pep") into transdecoder_ch_hmmer_short
        set sample_id, file("${sample_id}.combined.okay.fa.transdecoder.pep") into transdecoder_ch_signalp_short
        set sample_id, file("${sample_id}.combined.okay.fa.transdecoder.pep") into transdecoder_ch_tmhmm_short
        set sample_id, file("${sample_id}.combined.okay.fa.transdecoder.pep") into transdecoder_ch_trinotate_short
        set sample_id, file("${sample_id}.combined.okay.fa.transdecoder.pep") into transdecoder_ch_diamond_custom_short
        // added _short
        set sample_id, file("${sample_id}.transdecoder.stats_short") into transdecoder_summary_short
        // for test
        set sample_id, file("${sample_id}.combined.okay.fa.transdecoder.pep_short") into transdecoder_ch_busco_short

    script:
        """
        set +u
        source ${params.condash}/etc/profile.d/conda.sh
        conda activate TransPi

        unidb=${params.mypwd}/diamonddb_custom/${params.uniname}
        pf=${params.mypwd}/hmmerdb/${params.pfname}

        echo -e "\n-- TransDecoder.LongOrfs... --\n"

        ${params.torf} -t ${sample_id}.combined.okay.fa

        echo -e "\n-- Done with TransDecoder.LongOrfs --\n"

        echo -e "\n-- TransDecoder.Predict... --\n"

        ${params.tpred} -t ${sample_id}.combined.okay.fa

        echo -e "\n-- Done with TransDecoder.Predict --\n"

        echo -e "\n-- Calculating statistics... --\n"

        #Calculate statistics of Transdecoder
        echo "- Transdecoder stats for "${sample_id} >>${sample_id}.transdecoder.stats
        orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep -c ">" )
        echo "Total number of ORFs: "\$orfnum >>${sample_id}.transdecoder.stats
        echo -e "\t Of these ORFs" >>${sample_id}.transdecoder.stats
        orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep ">" | grep -c "|" )
        echo -e "\t\t with annotations: "\$orfnum >>${sample_id}.transdecoder.stats
        orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep ">" | grep -v "|" | grep -c ">" )
        echo -e "\t\t no annotation: "\$orfnum >>${sample_id}.transdecoder.stats
        echo >>${sample_id}.transdecoder.stats
        orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep -c "ORF type:complete" )
        echo -e "\t ORFs type=complete: "\$orfnum >>${sample_id}.transdecoder.stats
        orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep "ORF type:complete" | grep -c "|" )
        echo -e "\t\t with annotations: "\$orfnum >>${sample_id}.transdecoder.stats
        echo >>${sample_id}.transdecoder.stats
        orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep -c "ORF type:5prime_partial" )
        echo -e "\t ORFs type=5prime_partial: "\$orfnum >>${sample_id}.transdecoder.stats
        orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep "ORF type:5prime_partial" | grep -c "|" )
        echo -e "\t\t with annotations: "\$orfnum >>${sample_id}.transdecoder.stats
        echo >>${sample_id}.transdecoder.stats
        orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep -c "ORF type:3prime_partial" )
        echo -e "\t ORFs type=3prime_partial: "\$orfnum >>${sample_id}.transdecoder.stats
        orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep "ORF type:3prime_partial" | grep -c "|" )
        echo -e "\t\t with annotations: "\$orfnum >>${sample_id}.transdecoder.stats
        echo >>${sample_id}.transdecoder.stats
        orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep -c "ORF type:internal" )
        echo -e "\t ORFs type=internal: "\$orfnum >>${sample_id}.transdecoder.stats
        orfnum=\$( cat ${sample_id}.combined.okay.fa.transdecoder.pep | grep "ORF type:internal" | grep -c "|" )
        echo -e "\t\t with annotations: "\$orfnum >>${sample_id}.transdecoder.stats
        orfnum=0

        echo -e "\n-- Done with statistics --\n"

        echo -e "\n-- DONE with TransDecoder --\n"

        ## added _short for test
        cp ${sample_id}.transdecoder.stats ${params.mypwd}/results/${sample_id}.transdecoder.stats_short

        ## added _short for test
        cp ${sample_id}.combined.okay.fa.transdecoder.pep ${sample_id}.combined.okay.fa.transdecoder.pep_short
        """

}

// busco for test
process busco_transdecoder {

    label 'big_cpus'

    tag "${sample_id}"

    publishDir "${params.mypwd}/busco_transdecoder", mode: "copy", overwrite: true

    input:
        set sample_id, file("${sample_id}.combined.okay.fa.transdecoder.pep") from transdecoder_ch_busco
        set sample_id, file("${sample_id}.combined.okay.fa.transdecoder.pep_short") from transdecoder_ch_busco_short

    output:
        set sample_id, file("run_${sample_id}.transdecoder.bus"), file("run_${sample_id}.transdecoder_short.bus") into busco_ch_transdecoder
        set sample_id, file("short_summary_${sample_id}.transdecoder.bus.txt"), file("short_summary_${sample_id}.transdecoder_short.bus.txt") into busco_summary_transdecoder

    script:
        """
        set +u
        source ${params.condash}/etc/profile.d/conda.sh
        conda activate TransPi

        echo -e "\n-- Starting with BUSCO --\n"

        run_BUSCO.py -i ${sample_id}.combined.okay.fa.transdecoder.pep -o ${sample_id}.transdecoder.bus -l ${params.buscodb} -m prot -c ${task.cpus}

        run_BUSCO.py -i ${sample_id}.combined.okay.fa.transdecoder.pep_short -o ${sample_id}.transdecoder_short.bus -l ${params.buscodb} -m prot -c ${task.cpus}

        echo -e "\n-- DONE with BUSCO --\n"

        cp run_${sample_id}.transdecoder.bus/short_summary_${sample_id}.transdecoder.bus.txt ${params.mypwd}/results/short_summary_${sample_id}.transdecoder.bus.txt
        cp run_${sample_id}.transdecoder_short.bus/short_summary_${sample_id}.transdecoder_short.bus.txt ${params.mypwd}/results/short_summary_${sample_id}.transdecoder_short.bus.txt


        cp run_${sample_id}.transdecoder.bus/short_summary_${sample_id}.transdecoder.bus.txt .
        cp run_${sample_id}.transdecoder_short.bus/short_summary_${sample_id}.transdecoder_short.bus.txt .

        cp run_${sample_id}.transdecoder.bus/short_summary_${sample_id}.transdecoder.bus.txt ${params.mypwd}/results/short_summary_${sample_id}.transdecoder.bus.txt
        cp run_${sample_id}.transdecoder_short.bus/short_summary_${sample_id}.transdecoder_short.bus.txt ${params.mypwd}/results/short_summary_${sample_id}.transdecoder_short.bus.txt
        """
}

process swiss_diamond_trinotate {

    label 'big_cpus'

    tag "${sample_id}"

    input:
        set sample_id, file("${sample_id}.combined.okay.fa") from evigene_ch_diamond
        set sample_id, file("${sample_id}.combined.okay.fa.transdecoder.pep") from transdecoder_ch_diamond

    output:
        set sample_id, file("${sample_id}.diamond_blastx.outfmt6"), file("${sample_id}.diamond_blastp.outfmt6") into trinotate_ch_diamond

    script:
        """
        set +u
        source ${params.condash}/etc/profile.d/conda.sh
        conda activate TransPi

        swissdb=${params.mypwd}/diamonddb_swiss/uniprot_sprot.pep

        #Diamond (BLAST) Homologies

        echo -e "\n-- Starting with Diamond (blastx) --\n"

        ${params.diam} blastx -d \$swissdb -q ${sample_id}.combined.okay.fa -p ${task.cpus} -f 6 -k 1 -e 0.001 >${sample_id}.diamond_blastx.outfmt6

        echo -e "\n-- Done with Diamond (blastx) --\n"

        echo -e "\n-- Starting with Diamond (blastp) --\n"

        ${params.diam} blastp -d \$swissdb -q ${sample_id}.combined.okay.fa.transdecoder.pep -p ${task.cpus} -f 6 -k 1 -e 0.001 >${sample_id}.diamond_blastp.outfmt6

        echo -e "\n-- Done with Diamond (blastp)  --\n"
        """
}

process custom_diamond_trinotate {

    label 'big_cpus'

    tag "${sample_id}"

    input:
        set sample_id, file("${sample_id}.combined.okay.fa") from evigene_ch_trinotate_custom
        set sample_id, file("${sample_id}.combined.okay.fa.transdecoder.pep") from transdecoder_ch_diamond_custom

    output:
        set sample_id, file("${sample_id}.custom.diamond_blastx.outfmt6"), file("${sample_id}.custom.diamond_blastp.outfmt6") into trinotate_ch_diamond_custom

    script:
        """
        set +u
        source ${params.condash}/etc/profile.d/conda.sh
        conda activate TransPi

        unidb=${params.mypwd}/diamonddb_custom/${params.uniname}

        #Diamond (BLAST) Homologies

        echo -e "\n-- Starting with Diamond (blastx) --\n"

        ${params.diam} blastx -d \$unidb -q ${sample_id}.combined.okay.fa -p ${task.cpus} -f 6 -k 1 -e 0.001 >${sample_id}.custom.diamond_blastx.outfmt6

        echo -e "\n-- Done with Diamond (blastx) --\n"

        echo -e "\n-- Starting with Diamond (blastp) --\n"

        ${params.diam} blastp -d \$unidb -q ${sample_id}.combined.okay.fa.transdecoder.pep -p ${task.cpus} -f 6 -k 1 -e 0.001 >${sample_id}.custom.diamond_blastp.outfmt6

        echo -e "\n-- Done with Diamond (blastp)  --\n"
        """
}

process hmmer_trinotate {

    label 'low_cpus'

    tag "${sample_id}"

    input:
        set sample_id, file("${sample_id}.combined.okay.fa.transdecoder.pep") from transdecoder_ch_hmmer

    output:
        set sample_id, file("${sample_id}.TrinotatePFAM.out") into trinotate_ch_hmmer

    script:
        """
	    set +u
        source ${params.condash}/etc/profile.d/conda.sh
        conda activate TransPi

        pf=${params.mypwd}/hmmerdb/${params.pfname}

        echo -e "\n-- Starting with HMMER --\n"

        ${params.hmsc} --cpu ${task.cpus} --domtblout ${sample_id}.TrinotatePFAM.out \$pf ${sample_id}.combined.okay.fa.transdecoder.pep >pfam.log

        echo -e "\n-- Done with HMMER --\n"
        """
}

process signalP_trinotate {

    label 'low_cpus'

    tag "${sample_id}"

    input:
        set sample_id, file("${sample_id}.combined.okay.fa.transdecoder.pep") from transdecoder_ch_signalp

    output:
        set sample_id, file("${sample_id}.signalp.out") into trinotate_ch_signalp

    script:
        """
        #signalP to predict signal peptides

        echo -e "\n-- Starting with SignalP --\n"

        ${params.signalp} -f short -n ${sample_id}.signalp.out ${sample_id}.combined.okay.fa.transdecoder.pep

        echo -e "\n-- Done with SignalP --\n"
        """
}

process tmhmm_trinotate {

    label 'low_cpus'

    tag "${sample_id}"

    input:
        set sample_id, file("${sample_id}.combined.okay.fa.transdecoder.pep") from transdecoder_ch_tmhmm

    output:
        set sample_id, file("${sample_id}.tmhmm.out") into trinotate_ch_tmhmm

    script:
        """
        #tmHMM to predict transmembrane regions

        echo -e "\n-- Starting with tmHMM --\n"

        ${params.tmhmm} --short < ${sample_id}.combined.okay.fa.transdecoder.pep >${sample_id}.tmhmm.out

        echo -e "\n-- Done with tmHMM --\n"
        """
}

process rnammer_trinotate {

    label 'low_cpus'

    tag "${sample_id}"

    input:
        set sample_id, file("${sample_id}.combined.okay.fa") from evigene_ch_rnammer

    output:
        set sample_id, file("${sample_id}.combined.okay.fa.rnammer.gff") into trinotate_ch_rnammer

    script:
        """
        set +ue
        source ${params.condash}/etc/profile.d/conda.sh
        conda activate TransPi

        #RNAMMER to identify rRNA transcripts

        echo -e "\n-- Starting with RNAMMER --\n"

        ${params.rnaTri} --transcriptome ${sample_id}.combined.okay.fa --path_to_rnammer ${params.rnam}

        echo -e "\n-- Done with RNAMMER --\n"
        """
}

process trinotate {

    label 'low_cpus'

    tag "${sample_id}"

    input:
        set sample_id, file("${sample_id}.combined.okay.fa") from evigene_ch_trinotate
        set sample_id, file("${sample_id}.combined.okay.fa.transdecoder.pep") from transdecoder_ch_trinotate
        set sample_id, file("${sample_id}.diamond_blastx.outfmt6"), file("${sample_id}.diamond_blastp.outfmt6") from trinotate_ch_diamond
        set sample_id, file("${sample_id}.custom.diamond_blastx.outfmt6"), file("${sample_id}.custom.diamond_blastp.outfmt6") from trinotate_ch_diamond_custom
        set sample_id, file("${sample_id}.TrinotatePFAM.out") from trinotate_ch_hmmer
        set sample_id, file("${sample_id}.signalp.out") from trinotate_ch_signalp
        set sample_id, file("${sample_id}.tmhmm.out") from trinotate_ch_tmhmm
        set sample_id, file("${sample_id}.combined.okay.fa.rnammer.gff") from trinotate_ch_rnammer

    output:
        set sample_id, file("${sample_id}.GO.terms.txt"), file("${sample_id}.trinotate_annotation_report.xls") into trinotate_ch
        set sample_id, file("${sample_id}.GO.terms.txt") into trinotate_summary

    script:
        """
	    set +u
        source ${params.condash}/etc/profile.d/conda.sh
        conda activate TransPi

        #Generate gene_trans_map
        #Not using get_Trinity_gene_to_trans_map.pl since all the names are uniq
        cat ${sample_id}.combined.okay.fa | awk '{print \$1}' | grep ">" | cut -c 2- >a.txt

        paste a.txt a.txt >${sample_id}.combined.okay.fa.gene_trans_map

        #Get Trinotate.sqlite from folder (original)
        cp ${params.Tsql} .
        sqlname=`echo ${params.Tsql} | tr "\\/" "\n" | grep "\\.sqlite"`

        echo -e "\n-- Running Trinotate --\n"

        Trinotate \$sqlname init --gene_trans_map ${sample_id}.combined.okay.fa.gene_trans_map --transcript_fasta ${sample_id}.combined.okay.fa --transdecoder_pep ${sample_id}.combined.okay.fa.transdecoder.pep

        echo -e "\n-- Ending run of Trinotate --\n"

        echo -e "\n-- Loading hits and predictions to sqlite database... --\n"

        #Load protein hits
        Trinotate \$sqlname LOAD_swissprot_blastp ${sample_id}.diamond_blastp.outfmt6

        #Load transcript hits
        Trinotate \$sqlname LOAD_swissprot_blastx ${sample_id}.diamond_blastx.outfmt6

        #Load custom protein hits
        Trinotate \$sqlname LOAD_custom_blast --outfmt6 ${sample_id}.custom.diamond_blastp.outfmt6 --prog blastp --dbtype ${sample_id}_custom_uniprot

        #Load custom transcript hits
        Trinotate \$sqlname LOAD_custom_blast --outfmt6 ${sample_id}.custom.diamond_blastx.outfmt6 --prog blastx --dbtype ${sample_id}_custom_uniprot

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

        echo -e "\n-- Loading finished --\n"

        #Report

        echo -e "\n-- Generating report... --\n"

        Trinotate \$sqlname report >${sample_id}.trinotate_annotation_report.xls

        echo -e "\n-- Report generated --\n"

        #Extract info from XML file

        echo -e "\n-- Creating GO file from XML... --\n"

        ${params.XMLtoGO} --Trinotate_xls ${sample_id}.trinotate_annotation_report.xls --trans >${sample_id}.GO.terms.txt

        echo -e "\n-- Done with the GO --\n"

        cp ${sample_id}.trinotate_annotation_report.xls ${params.mypwd}/results/${sample_id}.trinotate_annotation_report.xls
        cp ${sample_id}.GO.terms.txt ${params.mypwd}/results/${sample_id}.GO.terms.txt

        echo -e "\n-- DONE with Trinotate --\n"
        """
}

process summary_evigene_individual {
    tag "${sample_id}"

    input:
        set sample_id, file("${sample_id}.combined.fa"), file("${sample_id}.combined.okay.fa") from evigene_summary_ind

    output:
        set sample_id, file("${sample_id}.sum_preEG.txt"), file("${sample_id}.sum_EG.txt") into final_sum_1

    script:
        """
        #Summary of total number of transcripts
        echo -e "- Number of transcripts before Evidential Genes\\n" >>${sample_id}.sum_preEG.txt
        echo "- Individual "${sample_id} >>${sample_id}.sum_preEG.txt
        echo -e "\\t Total transcripts:" >>${sample_id}.sum_preEG.txt
        num=\$( cat ${sample_id}.combined.fa | grep -c ">" )
        echo -e "\\t\\t \$num" >>${sample_id}.sum_preEG.txt
        echo -e "\\t Trinity" >>${sample_id}.sum_preEG.txt
        num=\$( cat ${sample_id}.combined.fa | grep -c ">TRINITY" )
        echo -e "\\t\\t \$num" >>${sample_id}.sum_preEG.txt
        echo -e "\\t SOAP" >>${sample_id}.sum_preEG.txt
        num=\$( cat ${sample_id}.combined.fa | grep -c ">SOAP" )
        echo -e "\\t\\t \$num" >>${sample_id}.sum_preEG.txt
        echo -e "\\t Velvet/Oases" >>${sample_id}.sum_preEG.txt
        num=\$( cat ${sample_id}.combined.fa | grep -c ">Velvet" )
        echo -e "\\t\\t \$num" >>${sample_id}.sum_preEG.txt
        echo -e "\\t IDBA_tran" >>${sample_id}.sum_preEG.txt
        num=\$( cat ${sample_id}.combined.fa | grep -c ">IDBA" )
        echo -e "\\t\\t \$num \\n" >>${sample_id}.sum_preEG.txt

        #Summary of transcripts after EvidentialGenes
        echo -e "- Number of transcripts by individual after EvidentialGenes\\n" >>${sample_id}.sum_EG.txt
        echo -e "- Individual "${sample_id} >>${sample_id}.sum_EG.txt
        echo -e "\\t Total transcripts:" >>${sample_id}.sum_EG.txt
        num=\$( cat ${sample_id}.combined.okay.fa | grep -c ">" )
        echo -e "\\t\\t \$num" >>${sample_id}.sum_EG.txt
        echo -e "\\t Trinity" >>${sample_id}.sum_EG.txt
        num=\$( cat ${sample_id}.combined.okay.fa | grep -c ">TRINITY" )
        echo -e "\\t\\t \$num" >>${sample_id}.sum_EG.txt
        echo -e "\\t SOAP" >>${sample_id}.sum_EG.txt
        num=\$( cat ${sample_id}.combined.okay.fa | grep -c ">SOAP" )
        echo -e "\\t\\t \$num" >>${sample_id}.sum_EG.txt
        echo -e "\\t Velvet/Oases" >>${sample_id}.sum_EG.txt
        num=\$( cat ${sample_id}.combined.okay.fa | grep -c ">Velvet" )
        echo -e "\\t\\t \$num" >>${sample_id}.sum_EG.txt
        echo -e "\\t IDBA_tran" >>${sample_id}.sum_EG.txt
        num=\$( cat ${sample_id}.combined.okay.fa | grep -c ">IDBA" )
        echo -e "\\t\\t \$num \\n" >>${sample_id}.sum_EG.txt

        cp ${sample_id}.sum_preEG.txt ${params.mypwd}/results/${sample_id}.sum_preEG.txt
        cp ${sample_id}.sum_EG.txt ${params.mypwd}/results/${sample_id}.sum_EG.txt
        """
}

// for test
process summary_evigene_cds_aa {
    tag "${sample_id}"

    input:
        set sample_id, file("${sample_id}.combined.cds"), file("${sample_id}.combined.okay.aa") from evigene_summary_cds_aa

    output:
        set sample_id, file("${sample_id}.sum_EG_cds.txt"), file("${sample_id}.sum_EG_aa.txt") into final_sum_aa

    script:
        """
        #Summary of total number of transcripts
        echo -e "- Number of transcripts before Evidential Genes\\n" >>${sample_id}.sum_EG_cds.txt
        echo "- Individual "${sample_id} >>${sample_id}.sum_EG_cds.txt
        echo -e "\\t Total transcripts:" >>${sample_id}.sum_EG_cds.txt
        num=\$( cat ${sample_id}.combined.cds | grep -c ">" )
        echo -e "\\t\\t \$num" >>${sample_id}.sum_EG_cds.txt
        echo -e "\\t Trinity" >>${sample_id}.sum_EG_cds.txt
        num=\$( cat ${sample_id}.combined.cds | grep -c ">TRINITY" )
        echo -e "\\t\\t \$num" >>${sample_id}.sum_EG_cds.txt
        echo -e "\\t SOAP" >>${sample_id}.sum_EG_cds.txt
        num=\$( cat ${sample_id}.combined.cds | grep -c ">SOAP" )
        echo -e "\\t\\t \$num" >>${sample_id}.sum_EG_cds.txt
        echo -e "\\t Velvet/Oases" >>${sample_id}.sum_EG_cds.txt
        num=\$( cat ${sample_id}.combined.cds | grep -c ">Velvet" )
        echo -e "\\t\\t \$num" >>${sample_id}.sum_EG_cds.txt
        echo -e "\\t IDBA_tran" >>${sample_id}.sum_EG_cds.txt
        num=\$( cat ${sample_id}.combined.cds | grep -c ">IDBA" )
        echo -e "\\t\\t \$num \\n" >>${sample_id}.sum_EG_cds.txt

        #Summary of transcripts after EvidentialGenes
        echo -e "- Number of transcripts by individual after EvidentialGenes\\n" >>${sample_id}.sum_EG_aa.txt
        echo -e "- Individual "${sample_id} >>${sample_id}.sum_EG_aa.txt
        echo -e "\\t Total transcripts:" >>${sample_id}.sum_EG_aa.txt
        num=\$( cat ${sample_id}.combined.okay.aa | grep -c ">" )
        echo -e "\\t\\t \$num" >>${sample_id}.sum_EG_aa.txt
        echo -e "\\t Trinity" >>${sample_id}.sum_EG_aa.txt
        num=\$( cat ${sample_id}.combined.okay.aa | grep -c ">TRINITY" )
        echo -e "\\t\\t \$num" >>${sample_id}.sum_EG_aa.txt
        echo -e "\\t SOAP" >>${sample_id}.sum_EG_aa.txt
        num=\$( cat ${sample_id}.combined.okay.aa | grep -c ">SOAP" )
        echo -e "\\t\\t \$num" >>${sample_id}.sum_EG_aa.txt
        echo -e "\\t Velvet/Oases" >>${sample_id}.sum_EG_aa.txt
        num=\$( cat ${sample_id}.combined.okay.aa | grep -c ">Velvet" )
        echo -e "\\t\\t \$num" >>${sample_id}.sum_EG_aa.txt
        echo -e "\\t IDBA_tran" >>${sample_id}.sum_EG_aa.txt
        num=\$( cat ${sample_id}.combined.okay.aa | grep -c ">IDBA" )
        echo -e "\\t\\t \$num \\n" >>${sample_id}.sum_EG_aa.txt

        cp ${sample_id}.sum_EG_cds.txt ${params.mypwd}/results/${sample_id}.sum_EG_cds.txt
        cp ${sample_id}.sum_EG_aa.txt ${params.mypwd}/results/${sample_id}.sum_EG_aa.txt
        """
}
// end test

process summary_busco_individual {
    tag "${sample_id}"

    input:
        set sample_id, file("short_summary_${sample_id}.fa.bus.txt"), file("short_summary_${sample_id}.cds.bus.txt"), file("short_summary_${sample_id}.aa.bus.txt")from busco_summary

    output:
        set sample_id, file("${sample_id}.sum_busco.txt") into final_sum_2

    script:
        """
        #Summary of BUSCO scores for all the final_assemblies
        echo -e "Summary of BUSCO \n" >>${sample_id}.sum_busco.txt
        echo "Using transcript" >>${sample_id}.sum_busco.txt
        cat short_summary_${sample_id}.fa.bus.txt >>${sample_id}.sum_busco.txt
        echo -e "\n Using CDS" >>${sample_id}.sum_busco.txt
        cat short_summary_${sample_id}.cds.bus.txt >>${sample_id}.sum_busco.txt
        echo -e "\n Using AA" >>${sample_id}.sum_busco.txt
        cat short_summary_${sample_id}.aa.bus.txt >>${sample_id}.sum_busco.txt

        cp ${sample_id}.sum_busco.txt ${params.mypwd}/results/${sample_id}.sum_busco.txt
        """
}

process summary_transdecoder_individual {
    tag "${sample_id}"

    input:
        set sample_id, file("${sample_id}.transdecoder.stats") from transdecoder_summary
        // for test
        set sample_id, file("${sample_id}.transdecoder.stats_short") from transdecoder_summary_short

    output:
        set sample_id, file("${sample_id}.sum_transdecoder.txt") into final_sum_3

    script:
        """
        #Summary of Transdecoder stats
        echo -e "Summary of Transdecoder \n" >>${sample_id}.sum_transdecoder.txt
        echo -e "- Long run (diamond and hmmer) \n" >>${sample_id}.sum_transdecoder.txt
        cat ${sample_id}.transdecoder.stats >>${sample_id}.sum_transdecoder.txt
        echo -e "- Short run \n" >>${sample_id}.sum_transdecoder.txt
        cat ${sample_id}.transdecoder.stats_short >>${sample_id}.sum_transdecoder.txt
        echo -e "##### \n" >>${sample_id}.sum_transdecoder.txt
        """
}

process summary_trinotate_individual {
    tag "${sample_id}"

    input:
        set sample_id, file("${sample_id}.GO.terms.txt") from trinotate_summary

    output:
        set sample_id, file("${sample_id}.sum_GO.txt") into final_sum_4

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

process get_combined_sum {
    publishDir "${params.mypwd}/results/summaries/", mode: "copy", overwrite: true

    input:
        set sample_id, file("${sample_id}.sum_preEG.txt"), file("${sample_id}.sum_EG.txt") from final_sum_1
        set sample_id, file("${sample_id}.sum_busco.txt") from final_sum_2
        set sample_id, file("${sample_id}.sum_transdecoder.txt") from final_sum_3
        set sample_id, file("${sample_id}.sum_GO.txt") from final_sum_4

    output:
        set sample_id, file("all_sum_preEG.txt"), file("all_sum_EG.txt"), file("all_sum_busco.txt"), file("all_sum_transdecoder.txt"), file("all_sum_GO.txt") into all_sum

    script:
        """
        cat *.sum_preEG.txt >all_sum_preEG.txt
        cat *.sum_EG.txt >all_sum_EG.txt
        cat *.sum_busco.txt >all_sum_busco.txt
        cat *.sum_transdecoder.txt >all_sum_transdecoder.txt
        cat *.sum_GO.txt >all_sum_GO.txt
        """
}

process get_run_info {
    script:
        """
        echo -e "-- Kmers used --" >>run_info.txt
        echo ${params.k} >>run_info.txt
        echo -e "\n-- Program versions --" >>run_info.txt

        v=\$( ${params.soap} --version | grep "version" | awk '{print \$2,\$3}' | cut -f 1 -d ":" | cut -f 2 -d " " )
        echo "SOAP:"\$v >>run_info.txt

        v=\$( ${params.vh} | grep "Version" | cut -f 2 -d " " )
        echo "Velveth:"\$v >>run_info.txt

        v=\$( ${params.vg} | grep "Version" | cut -f 2 -d " " )
        echo "Velvetg:"\$v >>run_info.txt

        v=\$( ${params.oa} | grep "Version" | cut -f 2 -d " " )
        echo "Oases:"\$v >>run_info.txt

        v=\$( echo "1.1.3" )
        echo "IDBA:"\$v >>run_info.txt

        v=\$( ${params.tr} --version | grep "version" | head -n 1 | cut -f 2 -d "-" )
        echo "Trinity:"\$v >>run_info.txt

        v=\$( diamond --version | cut -f 3 -d " " )
        echo "Diamond:"\$v >>run_info.txt

        v=\$( hmmsearch -h | head -n 2 | cut -f 3 -d " " | grep [0-9] )
        echo "HMMER:"\$v >>run_info.txt

        v=\$( echo "2019.01.01" )
        echo "EvidentialGene:"\$v >>run_info.txt

        v=\$( echo "1.2" )
        echo "RNAmmer:"\$v >>run_info.txt

        v=\$( TransDecoder.LongOrfs --version | cut -f 2 -d " " )
        echo "Transdecoder:"\$v >>run_info.txt

        v=\$( echo ${params.uniname} )
        echo "Uniprot_DB:"\$v >>run_info.txt

        v=\$( head -n 1 ${params.pfloc} | cut -f 2 -d "[" | cut -f 1 -d "]" | tr "|" "-" | tr -d " " )
        echo "Pfam_A:"\$v >>run_info.txt

        v=\$( run_BUSCO.py -v | cut -f 2 -d " " )
        echo "BUSCO:"\$v >>run_info.txt

        v=\$( echo ${params.buscodb} | tr "/" "\n" | tail -n 1 )
        echo "BUSCO_DB:"\$v >>run_info.txt

        v=\$( echo "TRI" )
        echo "Trinotate:"\$v >>run_info.txt

        v=\$( echo "4.1" )
        echo "SignalP:"\$v >>run_info.txt

        v=\$( echo "2.0" )
        echo "tmhmm:"\$v >>run_info.txt

        cp run_info.txt ${params.mypwd}/results/run_info.txt
        """
}




// Trinotate test

process swiss_diamond_trinotate_evigene {

    label 'big_cpus'

    tag "${sample_id}"

    input:
        set sample_id, file("${sample_id}.combined.okay.fa"), file("${sample_id}.combined.okay.aa") from evigene_ch_diamond_evigene

    output:
        set sample_id, file("${sample_id}.evigene.diamond_blastx.outfmt6"), file("${sample_id}.evigene.diamond_blastp.outfmt6") into trinotate_ch_diamond_evigene

    script:
        """
        set +u
        source ${params.condash}/etc/profile.d/conda.sh
        conda activate TransPi

        swissdb=${params.mypwd}/diamonddb_swiss/uniprot_sprot.pep

        #Diamond (BLAST) Homologies

        echo -e "\n-- Starting with Diamond (blastx) --\n"

        ${params.diam} blastx -d \$swissdb -q ${sample_id}.combined.okay.fa -p ${task.cpus} -f 6 -k 1 -e 0.001 >${sample_id}.evigene.diamond_blastx.outfmt6

        echo -e "\n-- Done with Diamond (blastx) --\n"

        echo -e "\n-- Starting with Diamond (blastp) --\n"

        ${params.diam} blastp -d \$swissdb -q ${sample_id}.combined.okay.aa -p ${task.cpus} -f 6 -k 1 -e 0.001 >${sample_id}.evigene.diamond_blastp.outfmt6

        echo -e "\n-- Done with Diamond (blastp)  --\n"
        """
}

process custom_diamond_trinotate_evigene {

    label 'big_cpus'

    tag "${sample_id}"

    input:
        set sample_id, file("${sample_id}.combined.okay.fa"), file("${sample_id}.combined.okay.aa") from evigene_ch_diamond_evigene_custom

    output:
        set sample_id, file("${sample_id}.evigene.custom.diamond_blastx.outfmt6"), file("${sample_id}.evigene.custom.diamond_blastp.outfmt6") into trinotate_ch_diamond_custom_evigene

    script:
        """
        set +u
        source ${params.condash}/etc/profile.d/conda.sh
        conda activate TransPi

        unidb=${params.mypwd}/diamonddb_custom/${params.uniname}

        #Diamond (BLAST) Homologies

        echo -e "\n-- Starting with Diamond (blastx) --\n"

        ${params.diam} blastx -d \$unidb -q ${sample_id}.combined.okay.fa -p ${task.cpus} -f 6 -k 1 -e 0.001 >${sample_id}.evigene.custom.diamond_blastx.outfmt6

        echo -e "\n-- Done with Diamond (blastx) --\n"

        echo -e "\n-- Starting with Diamond (blastp) --\n"

        ${params.diam} blastp -d \$unidb -q ${sample_id}.combined.okay.aa -p ${task.cpus} -f 6 -k 1 -e 0.001 >${sample_id}.evigene.custom.diamond_blastp.outfmt6

        echo -e "\n-- Done with Diamond (blastp)  --\n"
        """
}

process hmmer_trinotate_evigene {

    label 'low_cpus'

    tag "${sample_id}"

    input:
        set sample_id, file("${sample_id}.combined.okay.aa") from evigene_ch_hmmer_trinotate_evigene

    output:
        set sample_id, file("${sample_id}.evigene.TrinotatePFAM.out") into trinotate_ch_hmmer_evigene

    script:
        """
	    set +u
        source ${params.condash}/etc/profile.d/conda.sh
        conda activate TransPi

        pf=${params.mypwd}/hmmerdb/${params.pfname}

        echo -e "\n-- Starting with HMMER --\n"

        ${params.hmsc} --cpu ${task.cpus} --domtblout ${sample_id}.evigene.TrinotatePFAM.out \$pf ${sample_id}.combined.okay.aa >pfam.log

        echo -e "\n-- Done with HMMER --\n"
        """
}

process signalP_trinotate_evigene {

    label 'low_cpus'

    tag "${sample_id}"

    input:
        set sample_id, file("${sample_id}.combined.okay.aa") from evigene_ch_signalP_trinotate_evigene

    output:
        set sample_id, file("${sample_id}.evigene.signalp.out") into trinotate_ch_signalp_evigene

    script:
        """
        #signalP to predict signal peptides

        echo -e "\n-- Starting with SignalP --\n"

        ${params.signalp} -f short -n ${sample_id}.evigene.signalp.out ${sample_id}.combined.okay.aa

        echo -e "\n-- Done with SignalP --\n"
        """
}

process tmhmm_trinotate_evigene {

    label 'low_cpus'

    tag "${sample_id}"

    input:
        set sample_id, file("${sample_id}.combined.okay.aa") from evigene_ch_tmhmm_trinotate_evigene

    output:
        set sample_id, file("${sample_id}.evigene.tmhmm.out") into trinotate_ch_tmhmm_evigene

    script:
        """
        #tmHMM to predict transmembrane regions

        echo -e "\n-- Starting with tmHMM --\n"

        ${params.tmhmm} --short < ${sample_id}.combined.okay.aa >${sample_id}.evigene.tmhmm.out

        echo -e "\n-- Done with tmHMM --\n"
        """
}

process rnammer_trinotate_evigene {

    label 'low_cpus'

    tag "${sample_id}"

    input:
        set sample_id, file("${sample_id}.combined.okay.fa") from evigene_ch_rnammer_trinotate_evigene

    output:
        set sample_id, file("${sample_id}.combined.okay.fa.rnammer.gff") into trinotate_ch_rnammer_evigene

    script:
        """
        set +ue
        source ${params.condash}/etc/profile.d/conda.sh
        conda activate TransPi

        #RNAMMER to identify rRNA transcripts

        echo -e "\n-- Starting with RNAMMER --\n"

        ${params.rnaTri} --transcriptome ${sample_id}.combined.okay.fa --path_to_rnammer ${params.rnam}

        echo -e "\n-- Done with RNAMMER --\n"
        """
}

process trinotate_evigene {

    label 'low_cpus'

    tag "${sample_id}"

    input:
        set sample_id, file("${sample_id}.combined.okay.fa") from evigene_ch_trinotate_evigene
        set sample_id, file("${sample_id}.combined.okay.aa") from evigene_ch_trinotate_evigene_aa
        set sample_id, file("${sample_id}.evigene.diamond_blastx.outfmt6"), file("${sample_id}.evigene.diamond_blastp.outfmt6") from trinotate_ch_diamond_evigene
        set sample_id, file("${sample_id}.evigene.custom.diamond_blastx.outfmt6"), file("${sample_id}.evigene.custom.diamond_blastp.outfmt6") from trinotate_ch_diamond_custom_evigene
        set sample_id, file("${sample_id}.evigene.TrinotatePFAM.out") from trinotate_ch_hmmer_evigene
        set sample_id, file("${sample_id}.evigene.tmhmm.out") from trinotate_ch_tmhmm_evigene
        set sample_id, file("${sample_id}.evigene.signalp.out") from trinotate_ch_signalp_evigene
        set sample_id, file("${sample_id}.combined.okay.fa.rnammer.gff") from trinotate_ch_rnammer_evigene

    output:
        set sample_id, file("${sample_id}.evigene.GO.terms.txt"), file("${sample_id}.evigene.trinotate_annotation_report.xls") into trinotate_ch_evigene
        set sample_id, file("${sample_id}.evigene.GO.terms.txt") into trinotate_summary_evigene

    script:
        """
	    set +u
        source ${params.condash}/etc/profile.d/conda.sh
        conda activate TransPi

        #Generate gene_trans_map
        #Not using get_Trinity_gene_to_trans_map.pl since all the names are uniq
        cat ${sample_id}.combined.okay.fa | awk '{print \$1}' | grep ">" | cut -c 2- >a.txt

        paste a.txt a.txt >${sample_id}.combined.okay.fa.gene_trans_map

        #Get Trinotate.sqlite from folder (original)
        cp ${params.Tsql} .
        sqlname=`echo ${params.Tsql} | tr "\\/" "\n" | grep "\\.sqlite"`

        echo -e "\n-- Running Trinotate --\n"

        Trinotate \$sqlname init --gene_trans_map ${sample_id}.combined.okay.fa.gene_trans_map --transcript_fasta ${sample_id}.combined.okay.fa --transdecoder_pep ${sample_id}.combined.okay.aa

        echo -e "\n-- Ending run of Trinotate --\n"

        echo -e "\n-- Loading hits and predictions to sqlite database... --\n"

        #Load protein hits
        Trinotate \$sqlname LOAD_swissprot_blastp ${sample_id}.evigene.diamond_blastp.outfmt6

        #Load transcript hits
        Trinotate \$sqlname LOAD_swissprot_blastx ${sample_id}.evigene.diamond_blastx.outfmt6

        #Load custom protein hits
        Trinotate \$sqlname LOAD_custom_blast --outfmt6 ${sample_id}.evigene.custom.diamond_blastp.outfmt6 --prog blastp --dbtype ${sample_id}_custom_uniprot

        #Load custom transcript hits
        Trinotate \$sqlname LOAD_custom_blast --outfmt6 ${sample_id}.evigene.custom.diamond_blastx.outfmt6 --prog blastx --dbtype ${sample_id}_custom_uniprot

        #Load Pfam domain entries
        Trinotate \$sqlname LOAD_pfam ${sample_id}.evigene.TrinotatePFAM.out

        #Load transmembrane domains
        if [ -s ${sample_id}.tmhmm.out ];then
            Trinotate \$sqlname LOAD_tmhmm ${sample_id}.evigene.tmhmm.out
        else
            echo "No transmembrane domains (tmhmm)"
        fi

        #Load signal peptide predictions
        if [ -s ${sample_id}.signalp.out ];then
            Trinotate \$sqlname LOAD_signalp ${sample_id}.evigene.signalp.out
        else
            echo "No Signal-P"
        fi

        echo -e "\n-- Loading finished --\n"

        #Report

        echo -e "\n-- Generating report... --\n"

        Trinotate \$sqlname report >${sample_id}.evigene.trinotate_annotation_report.xls

        echo -e "\n-- Report generated --\n"

        #Extract info from XML file

        echo -e "\n-- Creating GO file from XML... --\n"

        ${params.XMLtoGO} --Trinotate_xls ${sample_id}.evigene.trinotate_annotation_report.xls --trans >${sample_id}.evigene.GO.terms.txt

        echo -e "\n-- Done with the GO --\n"

        cp ${sample_id}.evigene.trinotate_annotation_report.xls ${params.mypwd}/results/${sample_id}.evigene.trinotate_annotation_report.xls
        cp ${sample_id}.evigene.GO.terms.txt ${params.mypwd}/results/${sample_id}.evigene.GO.terms.txt

        echo -e "\n-- DONE with Trinotate --\n"
        """
}

// End trinotate test
