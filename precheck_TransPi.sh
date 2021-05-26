#!/usr/bin/env bash -e
export mypwd="$1"
os_c() {
    OS="$(uname)"
    if [ "$OS" == "Linux" ]; then
        echo -e "\n\t -- Downloading Linux Anaconda3 installation -- \n"
        curl -o Anaconda3-2020.11-Linux-x86_64.sh https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh
    else
        echo -e "\n\t\e[31m -- ERROR: Are you in a Linux system? Please check requirements and rerun the pre-check --\e[39m\n"
        exit 0
    fi
}
source_c() {
    if [ -f ~/.bashrc ];then
        source ~/.bashrc
    fi
}
conda_only() {
    source_c
    #Check conda and environment
    check_conda=$( command -v conda )
    if [ "$check_conda" != "" ];then #&& [ "$ver" -gt "45" ];then
        echo -e "\n\t -- Conda seems to be installed in your system --\n"
        ver=$( conda -V | awk '{print $2}' | cut -f 1,2 -d "." )
        vern=4.8
        if [ $( echo "$ver >= $vern" | bc -l ) -eq 1 ];then
            echo -e "\n\t -- Conda is installed (v4.8 or higher). Checking environment... --\n"
        fi
    else
        echo -e "\n\t -- Conda is not intalled --\n"
        os_c
        echo -e "\n\t -- Starting Anaconda installation -- \n"
        bash Anaconda3-20*.sh
        echo -e "\n\t -- Installation done -- \n"
        rm Anaconda3-20*.sh
        source_c
    fi
}
conda_c() {
    source_c
    #Check conda and environment
    check_conda=$( command -v conda )
    if [ "$check_conda" != "" ];then #&& [ "$ver" -gt "45" ];then
        echo -e "\n\t -- Conda seems to be installed in your system --\n"
        ver=$( conda -V | awk '{print $2}' | cut -f 1,2 -d "." )
        vern=4.8
        if [ $( echo "$ver >= $vern" | bc -l ) -eq 1 ];then
            echo -e "\n\t -- Conda is installed (v4.8 or higher). Checking environment... --\n"
            #Check environment
            check_env=$( conda info -e | awk '$1 == "TransPi" {print $2}' | wc -l )
	        if [ "$check_env" -eq 0 ];then
                echo -e "\n\t -- TransPi environment has not been created. Checking environment file... --\n"
                if [ -f ${confDir}/transpi_env.yml ];then
                    echo -e "\n\t -- TransPi environment file found. Creating environment... --\n"
                    conda env create -f ${confDir}/transpi_env.yml
                else
                    echo -e "\n\t\e[31m -- ERROR: TransPi environment file not found (transpi_env.yml). Please run the precheck in the TransPi directory. See manual for more info --\e[39m\n"
                    exit 0
                fi
            elif [ "$check_env" -eq 1 ];then
                echo -e "\n\t -- TransPi environment is installed and ready to be used --\n"
            fi
        fi
    else
        echo -e "\n\t -- Conda is not intalled --\n"
        os_c
        echo -e "\n\t -- Starting Anaconda installation -- \n"
        bash Anaconda3-20*.sh
        echo -e "\n\t -- Installation done -- \n"
        rm Anaconda3-20*.sh
        source_c
        if [ -f ${confDir}/transpi_env.yml ];then
            echo -e "\n\t -- TransPi environment file found. Creating environment... --\n"
            conda env create -f ${confDir}/transpi_env.yml
        else
            echo -e "\n\t\e[31m -- ERROR: TransPi environment file not found (transpi_env.yml). Please run the precheck in the TransPi directory. See manual for more info --\e[39m\n"
            exit 0
        fi
    fi
}
dir_c () {
    cd $mypwd
    if [ ! -d scripts/ ];then
        mkdir scripts
    fi
    if [ ! -d DBs ];then
        mkdir DBs
    fi
}
bus_dow () {
    name=$1
    cd $mypwd
    if [ ! -d DBs/busco_db/ ];then
        echo -e "\n\t -- Creating directory for the BUSCO V4 database --\n"
        mkdir -p DBs/busco_db
        cd DBs/busco_db
        bname=$( echo $name | cut -f 1 -d "_" )
        if [ `cat ${confDir}/conf/busV4list.txt | grep "${bname};" | wc -l` -eq 1 ];then
            echo -e "\n\t -- Downloading BUSCO V4 \"$name\" database --\n";wait
            wname=$( cat ${confDir}/conf/busV4list.txt | grep "${bname};" | cut -f 2 -d ";" )
            wget $wname
            echo -e "\n\t -- Preparing files ... --\n";wait
            tname=$( cat ${confDir}/conf/busV4list.txt | grep "${bname};" | cut -f 1 -d ";" | tr [A-Z] [a-z] )
            tar -xf ${tname}*.tar.gz
            rm ${tname}*.tar.gz
            echo -e "\n\t -- DONE with BUSCO V4 database --\n";wait
        fi
        dname=$( cat ${confDir}/conf/busV4list.txt | grep "${bname};" | cut -f 1 -d ";" | tr [A-Z] [a-z] )
        if [ -d ${dname}_odb10 ];then
            export busna=${dname}_odb10
        fi
    elif [ -d DBs/busco_db/ ];then
        cd DBs/busco_db
        bname=$( echo $name | cut -f 1 -d "_" )
        dname=$( cat ${confDir}/conf/busV4list.txt | grep "${bname};" | cut -f 1 -d ";" | tr [A-Z] [a-z] )
        if [ -d ${dname}_odb10 ];then
            echo -e "\n\t -- BUSCO V4 \"$name\" database found -- \n"
            export busna=${dname}_odb10
        else
            bname=$( echo $name | cut -f 1 -d "_" )
            if [ `cat ${confDir}/conf/busV4list.txt | grep "${bname};" | wc -l` -eq 1 ];then
                echo -e "\n\t -- Downloading BUSCO V4 \"$name\" database --\n";wait
                wname=$( cat ${confDir}/conf/busV4list.txt | grep "${bname};" | cut -f 2 -d ";" )
                wget $wname
                echo -e "\n\t -- Preparing files ... --\n";wait
                tname=$( cat ${confDir}/conf/busV4list.txt | grep "${bname};" | cut -f 1 -d ";" | tr [A-Z] [a-z] )
                tar -xvf ${tname}*.tar.gz
                rm ${tname}*.tar.gz
                echo -e "\n\t -- DONE with BUSCO V4 database --\n";wait
            fi
            dname=$( cat ${confDir}/conf/busV4list.txt | grep "${bname};" | cut -f 1 -d ";" | tr [A-Z] [a-z] )
            if [ -d ${dname}_odb10 ];then
                export busna=${dname}_odb10
            fi
        fi
    fi
}
bus_c () {
    cd $mypwd
    echo -e "\n\t -- Selecting BUSCO V4 database -- \n"
    PS3="
    Please select one (1-5): "
    if [ -f ${confDir}/conf/busV4list.txt ];then
    select var in `cat ${confDir}/conf/busV4list.txt | grep "###" | tr -d "#"`;do
    case $var in
        BACTERIA)
            echo -e "\n\t You selected BACTERIA. Which specific database? \n"
            PS3="
	    Please select database: "
            select var1 in `cat ${confDir}/conf/busV4list.txt | sed -n "/##BACTERIA/,/#MAIN/p" | grep -v "##" | tr -d "#"`;do
    	    case $var1 in
    	        MAIN_MENU)
                    bus_c
                ;;
                *)
                if [ "$var1" != "" ];then
                    if [ `cat ${confDir}/conf/busV4list.txt | grep -c "$var1"` -ge 1 ];then
                        bus_dow $var1
                    fi
                else
                    echo -e "\n\t Wrong option. Try again \n"
                    bus_c
                fi
            ;;
            esac
            break
            done
	   ;;
       EUKARYOTA)
            echo -e "\n\tYou selected EUKARYOTA. Which specific database? \n"
            PS3="
	    Please select database: "
            select var1 in `cat ${confDir}/conf/busV4list.txt | sed -n "/##EUKARYOTA/,/#MAIN/p" | grep -v "##" | tr -d "#"`;do
        	case $var1 in
        	    MAIN_MENU)
                    bus_c
                ;;
                Arthropoda_\(Phylum\))
                    select var2 in `cat ${confDir}/conf/busV4list.txt | sed -n "/##ARTHROPODA/,/#MAIN/p" | grep -v "##" | tr -d "#"`;do
                    case $var2 in
                    MAIN_MENU)
                        bus_c
                    ;;
                    *)
                    if [ "$var2" != "" ];then
                        if [ `cat ${confDir}/conf/busV4list.txt | grep -c "$var2"` -ge 1 ];then
                            bus_dow $var2
                        fi
                    else
                        echo -e "\n\t Wrong option. Try again \n"
                        bus_c
                    fi
                    esac
                    break
                    done
                ;;
                Fungi_\(Kingdom\))
                    select var2 in `cat ${confDir}/conf/busV4list.txt | sed -n "/##FUNGI/,/#MAIN/p" | grep -v "##" | tr -d "#"`;do
                    case $var2 in
                    MAIN_MENU)
                        bus_c
                    ;;
                    *)
                    if [ "$var2" != "" ];then
                        if [ `cat ${confDir}/conf/busV4list.txt | grep -c "$var2"` -ge 1 ];then
                            bus_dow $var2
                        fi
                    else
                        echo -e "\n\t Wrong option. Try again \n"
                        bus_c
                    fi
                    esac
                    break
                    done
                ;;
                Plants_\(Kingdom\))
                    select var2 in `cat ${confDir}/conf/busV4list.txt | sed -n "/##PLANTS/,/#MAIN/p" | grep -v "##" | tr -d "#"`;do
                    case $var2 in
                    MAIN_MENU)
                        bus_c
                    ;;
                    *)
                    if [ "$var2" != "" ];then
                        if [ `cat ${confDir}/conf/busV4list.txt | grep -c "$var2"` -ge 1 ];then
                            bus_dow $var2
                        fi
                    else
                        echo -e "\n\t Wrong option. Try again \n"
                        bus_c
                    fi
                    esac
                    break
                    done
                ;;
                Protists_\(Clade\))
                    select var2 in `cat ${confDir}/conf/busV4list.txt | sed -n "/##PROTIST/,/#MAIN/p" | grep -v "##" | tr -d "#"`;do
                    case $var2 in
                    MAIN_MENU)
                        bus_c
                    ;;
                    *)
                    if [ "$var2" != "" ];then
                        if [ `cat ${confDir}/conf/busV4list.txt | grep -c "$var2"` -ge 1 ];then
                            bus_dow $var2
                        fi
                    else
                        echo -e "\n\t Wrong option. Try again \n"
                        bus_c
                    fi
                    esac
                    break
                    done
                ;;
                Vertebrata_\(Sub_phylum\))
                    select var2 in `cat ${confDir}/conf/busV4list.txt | sed -n "/##VERTEBRATA/,/#MAIN/p" | grep -v "##" | tr -d "#"`;do
                    case $var2 in
                    MAIN_MENU)
                        bus_c
                    ;;
                    *)
                    if [ "$var2" != "" ];then
                        if [ `cat ${confDir}/conf/busV4list.txt | grep -c "$var2"` -ge 1 ];then
                            bus_dow $var2
                        fi
                    else
                        echo -e "\n\t Wrong option. Try again \n"
                        bus_c
                    fi
                    esac
                    break
                    done
                ;;
                *)
                if [ "$var1" != "" ];then
                    if [ `cat ${confDir}/conf/busV4list.txt | grep -c "$var1"` -ge 1 ];then
                        bus_dow $var1
                    fi
                else
                    echo -e "\n\t Wrong option. Try again \n"
                    bus_c
                fi
                ;;
            esac
            break
            done
        ;;
        ARCHAEA)
            echo -e "\n\tYou selected ARCHAEA. Which specific database? \n"
            PS3="
	    Please select database: "
            select var1 in `cat ${confDir}/conf/busV4list.txt | sed -n "/##ARCHAEA/,/#MAIN/p" | grep -v "##" | tr -d "#"`;do
            case $var1 in
            	MAIN_MENU)
                    bus_c
                ;;
                *)
                if [ "$var1" != "" ];then
                    if [ `cat ${confDir}/conf/busV4list.txt | grep -c "$var1"` -ge 1 ];then
                        bus_dow $var1
                    fi
                else
                    echo -e "\n\t Wrong option. Try again \n"
                    bus_c
                fi
                ;;
            esac
            break
            done
        ;;
        EXIT)
            echo -e "\n\t Exiting \n"
            exit 0
        ;;
        *)
            echo -e "\n\t Wrong option. Try again \n"
            bus_c
            ;;
    esac
    break
    done
    else
        echo -e "\n\t\e[31m -- ERROR: Please make sure that file \"busV4list.txt\" is available. Please run the precheck in the TransPi directory. See manual for more info --\e[39m\n\n"
	    exit 0
    fi
}
uni_c () {
    PS3="
    Please select UNIPROT database to use: "
    select var in `ls *`;do
        if [ "$var" != "" ];then
            if [ `echo $var | grep ".gz" | wc -l` -eq 1 ];then
                echo -e "\n\n\t -- File is compressed -- \n"
                echo -e "\n\n\t -- Uncompressing file ... -- \n"
                gunzip $var
                echo -e "\n\t -- UNIPROT database selected: \"${var%.gz}\" --\n"
                export unina=${var%.gz}
            else
                echo -e "\n\t -- UNIPROT database selected: \"$var\" --\n"
                export unina=${var}
            fi
        else
            echo -e "\n\t Wrong option. Try again \n"
            uni_c
        fi
    break
    done
}
unicomp_c () {
    echo -e -n "\n\t    Do you want to uncompress the file(s)? (y,n,exit): "
    read ans
    case $ans in
        [yY] | [yY][eE][sS])
            echo -e "\n\n\t -- Uncompressing file(s) ... -- \n"
            gunzip *.gz
        ;;
        [nN] | [nN][oO])
            echo -e "\n\n\t\e[31m -- ERROR: Please uncompress the file(s) and rerun the pre-check  --\e[39m\n"
            exit 0
        ;;
        exit)
            echo -e "\n\t -- Exiting -- \n"
            exit 0
        ;;
        *)
            echo -e "\n\n\t\e[31m -- Yes or No answer not specified. Try again --\e[39m\n"
            unicomp_c
        ;;
    esac
}
uniprot_user_DB(){
    echo -e -n "\n\t -- Provide the PATH where to locate your proteins file: "
    read -e ans
    if [ -d ${ans} ];then
        echo -e "\n\t -- Directory ${ans} found -- \n"
        cd ${ans}
        uni_c
    elif [ -d $( dirname ${ans} ) ];then
        echo -e "\n\t -- Directory "$( dirname ${ans} )" found -- \n"
        cd $( dirname ${ans} )
        uni_c
    else
        echo -e "\n\t\e[31m -- Directory ${ans} not found --\e[39m\n"
        uniprot_meta
    fi
}
uniprot_taxon_DB(){
    echo -e "\n\t -- Input the Taxon ID (Taxonomy ID, NCBI txid) of your interest. TransPi will download the proteins from UNIPROT --"
    echo -e "\t    Example: metazoan TaxID = 33208 -- \n"
    echo -e -n "\n\t -- Your Taxon ID (only the numbers): "
    read ans
    echo -e "\n\t -- Downloading UNIPROT proteins from Taxon ID: $ans -- \n"
    curl -o uniprot_${ans}.fasta.gz "https://www.uniprot.org/uniprot/?query=taxonomy:${ans}&format=fasta&compress=yes&include=no"
    gunzip uniprot_${ans}.fasta.gz
    date -u >.lastrun.txt
    uni_c
}
uniprot_meta () {
    myuni=$( pwd )
    echo -e "\n\t -- TransPi uses a custom protein database (one of many) from UNIPROT for the annotation -- \n"
    echo "
        Options available:

            1- Download metazoan proteins from UNIPROT

            2- Provide the PATH of my DB

            3- Provide UNIPROT Taxon ID

            4- Skip for now

    "
    echo -e -n "\t Which option you want? "
    read ans
    case $ans in
        1)
            echo -e "\n\n\t -- Downloading current metazoan protein dataset from UNIPROT -- \n"
            echo -e "\n\t -- This could take a couple of minutes depending on connection. Please wait -- \n"
            curl -o uniprot_metazoa_33208.fasta.gz "https://www.uniprot.org/uniprot/?query=taxonomy:33208&format=fasta&compress=yes&include=no"
            echo -e "\n\t -- Uncompressing uniprot_metazoa_33208.fasta.gz ... -- \n"
            gunzip uniprot_metazoa_33208.fasta.gz
            date -u >.lastrun.txt
            uni_c
        ;;
        2)
            uniprot_user_DB
        ;;
        3)
            uniprot_taxon_DB
        ;;
        4)
            echo -e "\n\t -- Skipping UNIPROT DB -- \n"
        ;;
        *)
            echo -e "\n\t\e[31m -- Wrong option. Try again --\e[39m\n"
            uniprot_meta
        ;;
    esac
}
uniprot_c () {
    #Check UNIPROT
    cd $mypwd
    if [ ! -d DBs/uniprot_db/ ];then
        echo -e "\n\t -- Creating directory for the UNIPROT database --\n"
        mkdir -p DBs/uniprot_db/
        cd DBs/uniprot_db/
        uniprot_meta
    elif [ -d DBs/uniprot_db/ ];then
        cd DBs/uniprot_db/
        myuni=$( pwd )
        echo -e "\n\t -- UNIPROT database directory found at: $myuni -- \n"
        myfasta=$( ls -1 | grep -v ".gz" | egrep ".fasta|.fa" | wc -l )
        myfastagz=$( ls -1 | egrep ".fasta.gz|.fa.gz" | wc -l )
        if [ $myfasta -eq 0 ] && [ $myfastagz -eq 0 ];then
            echo -e "\n\t -- Directory \"$myuni\" is empty --\n"
            uniprot_meta
        else
            echo -e "\n\t -- Here is the list of UNIPROT files found at: $myuni -- \n"
            uni_c
        fi
    fi
}
java_c () {
	export NXF_VER=20.10.0 && curl -s https://get.nextflow.io | bash 2>.error_nextflow
	check_err=$( head -n 1 .error_nextflow | grep -c "java: command not found" )
	if [ $check_err -eq 1 ];then
		echo -e "\n\t\e[31m -- ERROR: Please install Java 1.8 (or later). Requirement for Nextflow --\e[39m\n"
		exit 0
	fi
	rm .error_nextflow
}
nextflow_c () {
    #Check Nextflow
    cd $mypwd
    check_next=$( command -v nextflow | wc -l )
    if [ $check_next -eq 1 ];then
        echo -e "\n\t -- Nextflow is installed -- \n"
    elif [ $check_next -eq 0 ];then
	check_next=$( ls -1 | grep -v "nextflow.config" | grep -c "nextflow" )
        if [ $check_next -eq 1 ];then
            echo -e "\n\t -- Nextflow is installed -- \n"
	    else
            echo -e -n "\n\t    Do you want to install Nextflow? (y or n): "
            read ans
            case $ans in
                [yY] | [yY][eE][sS])
                    echo -e "\n\t -- Downloading Nextflow ... -- \n"
                    java_c
		    		echo -e "\n\t -- Nextflow is now installed on $mypwd (local installation) -- \n"
                ;;
                [nN] | [nN][oO])
                    echo -e "\n\n\t\e[31m -- ERROR: Download and Install Nextflow. Then rerun the pre-check  --\e[39m\n"
                    exit 0
                ;;
                *)
                    echo -e "\n\n\t\e[31m -- Yes or No answer not specified. Try again --\e[39m\n"
                    nextflow_c
                ;;
            esac
	    fi
    fi
}
evi_c () {
	cd ${confDir}
    check_evi=$( command -v tr2aacds.pl | wc -l )
    if [ $check_evi -eq 0 ];then
        if [ ! -d ${confDir}/scripts/evigene/ ];then
            echo -e "\n\t -- EvidentialGene is not installed -- \n"
            mkdir -p ${confDir}/scripts && cd ${confDir}/scripts
            echo -e "\n\t -- Downloading EvidentialGene -- \n"
            wget http://arthropods.eugenes.org/EvidentialGene/other/evigene_older/evigene19may14.tar
            tar -xf evigene19may14.tar
            rm evigene19may14.tar
            echo -e "\n\t -- Done with EvidentialGene -- \n"
        else
            echo -e "\n\t -- EvidentialGene directory was found at ${confDir}/scripts (local installation) -- \n"
        fi
    elif [ $check_evi -eq 1 ];then
        echo -e "\n\t -- EvidentialGene is already installed and in the PATH  -- \n"
    fi
}
trisql_container () {
    if [ ! -e *.sqlite ];then
        echo -e "\n\n\t -- Custom sqlite database for Trinotate is not installed -- \n"
        echo -e "\n\t -- This could take a couple of minutes depending on connection. Please wait -- \n"
        rm -rf *
        wget https://github.com/Trinotate/Trinotate/archive/Trinotate-v3.2.1.tar.gz
        tar -xf Trinotate-v3.2.1.tar.gz
        mv Trinotate-Trinotate-v3.2.1/ Trinotate_build_scripts/
        ./Trinotate_build_scripts/admin/Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate
        rm uniprot_sprot.dat.gz Pfam-A.hmm.gz
        date -u >.lastrun.txt
    elif [ -e *.sqlite ];then
        echo -e "\n\t -- Custom sqlite database for Trinotate found at "${mypwd}/DBs/sqlite_db" -- \n"
        DB=$( if [ -f ${mypwd}/DBs/sqlite_db/.lastrun.txt ];then cat .lastrun.txt;else echo "N/A";fi )
        echo -e "\n\t -- Databases (PFAM,SwissProt,EggNOG,GO) last update: ${DB} --\n "
    fi
}
trisql_c () {
    source ~/.bashrc
    check_conda=$( command -v conda )
    if [ "$check_conda" == "" ];then
        echo -e "\n\t\e[31m -- Looks like conda is not installed--\e[39m\n"
        exit 0
    fi
    if [ ! -e *.sqlite ];then
        echo -e "\n\t -- Custom sqlite database for Trinotate is not installed -- \n"
        echo -e "\n\t -- This could take a couple of minutes depending on connection. Please wait -- \n"
        condaRoot=$( conda info --json | grep "CONDA_ROOT" | cut -f 2 -d ":" | tr -d "," | tr -d " " | tr -d "\"" )
        if [ -f ${condaRoot}/etc/profile.d/conda.sh ];then
            source ${condaRoot}/etc/profile.d/conda.sh
            conda activate TransPi
            check_sql=$( command -v Build_Trinotate_Boilerplate_SQLite_db.pl | wc -l )
            if [ $check_sql -eq 0 ];then
                echo -e "\n\t -- Script \"Build_Trinotate_Boilerplate_SQLite_db.pl\" from Trinotate cannot be found -- \n"
                echo -e "\n\t\e[31m -- Verify your conda installation --\e[39m\n"
                exit 0
            elif [ $check_sql -eq 1 ];then
                Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate
                rm uniprot_sprot.dat.gz Pfam-A.hmm.gz
                date -u >.lastrun.txt
            fi
        fi
    elif [ -e *.sqlite ];then
        echo -e "\n\t -- Custom sqlite database for Trinotate found at "${mypwd}/DBs/sqlite_db" -- \n"
        DB=$( if [ -f ${mypwd}/DBs/sqlite_db/.lastrun.txt ];then cat .lastrun.txt;else echo "N/A";fi )
        echo -e "\n\t -- Databases (PFAM,SwissProt,EggNOG,GO) last update: ${DB} --\n "
    fi
}
buildsql_c () {
    cd ${mypwd}
    if [ -d DBs/sqlite_db/ ];then
        cd DBs/sqlite_db/
    else
        mkdir -p DBs/sqlite_db/
        cd DBs/sqlite_db/
    fi
}
pfam_c() {
    #Check PFAM files
    cd $mypwd
    if [ ! -d DBs/hmmerdb/ ];then
        echo -e "\n\t -- Creating directory for the HMMER database --\n"
        mkdir -p DBs/hmmerdb/
        cd DBs/hmmerdb/
        echo -e "-- Downloading Pfam-A files ... --\n"
        wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
        echo -e "-- Preparing Pfam-A files ... --\n"
        gunzip Pfam-A.hmm.gz
        date -u >.lastrun.txt
    elif [ -d DBs/hmmerdb/ ];then
        echo -e "\n\t -- Directory for the HMMER database is present --\n"
        cd DBs/hmmerdb/
        if [ -f Pfam-A.hmm ];then
            echo -e "\n\t -- Pfam file is present and ready to be used --\n"
            DB=$( if [ -f ${mypwd}/DBs/hmmerdb/.lastrun.txt ];then cat .lastrun.txt;else echo "N/A";fi )
            echo -e "\n\t -- Pfam last update: ${DB} --\n"
        else
            echo -e "-- Downloading Pfam-A files ... --\n"
            wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
            echo -e "-- Preparing Pfam-A files ... --\n"
            gunzip Pfam-A.hmm.gz
            date -u >.lastrun.txt
        fi
    fi
}
#temporary for buscoV4
bus_env4 () {
    echo -e "\n\t -- Creating BUSCO V4 environment --\n"
    conda create -n busco4 -c conda-forge -c bioconda busco=4.1.4=py_0 -y
    conda clean -a -y
}
bus4 () {
    cd $mypwd
    cpath=$( conda info --json | sed -n '/\"envs\":/,/\],/p' | grep "busco4" | tr -d "," | tr -d " " | tr -d "\"" )
    if [ "$cpath" == "" ];then
        echo -e "\n\t -- Cannot find the BUSCO V4 environment -- \n"
        echo -e "\n\t -- Cleaning conda before installing BUSCO V4 environment -- \n"
        conda clean -a -y
        bus_env4
    else
        echo -e "\n\t -- BUSCO V4 environment found --\n"
    fi
}
pfam_u() {
    cd $installDir
    if [ ! -d DBs/hmmerdb/ ];then
        echo -e "\n\t -- Creating directory for the HMMER database --\n"
        mkdir -p DBs/hmmerdb/
        cd DBs/hmmerdb/
        echo -e "-- Downloading Pfam-A files ... --\n"
        wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
        echo -e "-- Preparing Pfam-A files ... --\n"
        gunzip Pfam-A.hmm.gz
        date -u >.lastrun.txt
    elif [ -d DBs/hmmerdb/ ];then
        echo -e "\n\t -- Directory for the HMMER database is present --\n"
        cd DBs/hmmerdb/
        rm -rf *
        echo -e "-- Downloading Pfam-A files ... --\n"
        wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
        echo -e "-- Preparing Pfam-A files ... --\n"
        gunzip Pfam-A.hmm.gz
        date -u >.lastrun.txt
    fi
}
sqld(){
    rm -rf *
    source ~/.bashrc
    check_conda=$( command -v conda )
    if [ "$check_conda" == "" ];then
        echo -e "\n\t\e[31m -- Looks like conda is not installed --\e[39m\n"
        echo -e "\n\t\e[31m -- Install conda and rerun this script --\e[39m\n"
        exit 0
    fi
    if [ ! -e *.sqlite ];then
        echo -e "\n\t -- Custom sqlite database for Trinotate is not installed -- \n"
        echo -e -n "\n\t    Do you want to install the custom sqlite database? (y or n): "
        read ans
        case $ans in
            [yY] | [yY][eE][sS])
                condaRoot=$( conda info --json | grep "CONDA_ROOT" | cut -f 2 -d ":" | tr -d "," | tr -d " " | tr -d "\"" )
                if [ -f ${condaRoot}/etc/profile.d/conda.sh ];then
                    source ${condaRoot}/etc/profile.d/conda.sh
                    conda activate TransPi
                    check_sql=$( command -v Build_Trinotate_Boilerplate_SQLite_db.pl | wc -l )
                    if [ $check_sql -eq 0 ];then
                        echo -e "\n\t -- Script \"Build_Trinotate_Boilerplate_SQLite_db.pl\" from Trinotate cannot be found -- \n"
                        echo -e "\n\t\e[31m -- Verify your conda installation --\e[39m\n"
                        exit 0
                    elif [ $check_sql -eq 1 ];then
                        echo -e "\n\t -- This could take a couple of minutes depending on connection. Please wait -- \n"
                        Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate
                        rm uniprot_sprot.dat.gz Pfam-A.hmm.gz
                        date -u >.lastrun.txt
                    fi
                fi
            ;;
            [nN] | [nN][oO])
                echo -e "\n\n\t\e[31m -- ERROR: Generate the custom trinotate sqlite database at "${mypwd}/DBs/sqlite_db". Then rerun the pre-check  --\e[39m\n"
                exit 0
            ;;
            *)
                echo -e "\n\n\t\e[31m -- Yes or No answer not specified. Try again --\e[39m\n"
                sqld
            ;;
        esac
    elif [ -e *.sqlite ];then
        echo -e "\n\t -- Custom sqlite database for Trinotate found at "${mypwd}/DBs/sqlite_db" -- \n"
        DB=$( if [ -f ${mypwd}/DBs/sqlite_db/.lastrun.txt ];then cat .lastrun.txt;else echo "N/A";fi )
        echo -e "\n\t -- Databases (PFAM,SwissProt,EggNOG,GO) last update: ${DB} --\n "
    fi
    pfam_u
}
ddate() {
    if [ ! -e .lastrun.txt ];then
        echo -e "\n\t -- No info about when the databse was created -- \n"
        echo -e -n "\n\t -- Do you want rerun script and update the databases? (y or n): "
        read ans
        case $ans in
            [yY] | [yY][eE][sS])
                sqld
            ;;
            [nN] | [nN][oO])
                echo -e "\n\n\t -- Exiting program -- \n"
                exit 0
            ;;
            *)
                echo -e "\n\n\t\e[31m -- Yes or No answer not specified. Try again --\e[39m\n"
                ddate
            ;;
        esac
    elif [ -e .lastrun.txt ];then
        a=$( cat .lastrun.txt )
        echo -e "\n\t -- Database was created on \e[32m${a}\e[39m -- \n"
        echo -e -n "\n\t -- Do you want rerun script and update the databases? (y or n): "
        read ans
        case $ans in
            [yY] | [yY][eE][sS])
                sqld
            ;;
            [nN] | [nN][oO])
                echo -e "\n\n\t -- Exiting program -- \n"
                exit
            ;;
            *)
                echo -e "\n\n\t\e[31m -- Yes or No answer not specified. Try again --\e[39m\n"
                ddate
            ;;
        esac
    fi
}
downd() {
    cd $mypwd
    if [ ! -d DBs/sqlite_db/ ];then
        echo -e "\n\t -- SQLite directory not found at ${mypwd}/DBs -- \n"
        echo -e "\n\t -- Creating ${mypwd}/DBs -- \n"
        mkdir -p $mypwd/DBs/sqlite_db/
        downd
    elif [ -d DBs/sqlite_db/ ];then
        echo -e "\n\t -- SQLite directory found at ${mypwd}/DBs -- \n"
        cd DBs/sqlite_db/
        if [ ! -e *.sqlite ];then
            sqld
        elif [ -e *.sqlite ];then
            echo -e "\n\t -- Custom sqlite database for Trinotate is installed -- \n"
            echo -e "\n\t -- Verifying when scripts was last run -- \n"
            ddate
        fi
    fi
}
get_var_container () {
    cd $mypwd
    echo "busco4db=$mypwd/DBs/busco_db/$busna" >>${mypwd}/.varfile.sh
    echo "uniname=$unina" >>${mypwd}/.varfile.sh
    echo "uniprot=$mypwd/DBs/uniprot_db/$unina" >>${mypwd}/.varfile.sh
    echo "pfloc=$mypwd/DBs/hmmerdb/Pfam-A.hmm" >>${mypwd}/.varfile.sh
    echo "pfname=Pfam-A.hmm" >>${mypwd}/.varfile.sh
    echo "nextflow=$mypwd/nextflow" >>${mypwd}/.varfile.sh
    echo "Tsql=$mypwd/DBs/sqlite_db/*.sqlite" >>${mypwd}/.varfile.sh
    echo "unpdate=\"$( if [ -f ${mypwd}/DBs/uniprot_db/.lastrun.txt ];then cat ${mypwd}/DBs/uniprot_db/.lastrun.txt;else echo "N/A";fi )\"" >>${mypwd}/.varfile.sh
    echo "pfdate=\"$( if [ -f ${mypwd}/DBs/hmmerdb/.lastrun.txt ];then cat ${mypwd}/DBs/hmmerdb/.lastrun.txt;else echo "N/A";fi )\"" >>${mypwd}/.varfile.sh
    echo "dbdate=\"$( if [ -f ${mypwd}/DBs/sqlite_db/.lastrun.txt ];then cat ${mypwd}/DBs/sqlite_db/.lastrun.txt;else echo "N/A";fi )\"" >>${mypwd}/.varfile.sh
    vpwd=$mypwd
    echo "mypwd=$mypwd" >>${vpwd}/.varfile.sh
    source .varfile.sh
    echo -e "\n\t -- INFO to use in TransPi --\n"
    echo -e "\t Installation PATH:\t $mypwd"
    echo -e "\t BUSCO V4 database:\t $busco4db"
    echo -e "\t UNIPROT database:\t $uniprot"
    echo -e "\t UNIPROT last update:\t $unpdate"
    echo -e "\t PFAM files:\t\t $pfloc"
    echo -e "\t PFAM last update:\t $pfdate"
    echo -e "\t SQL DB last update: \t $dbdate"
    echo -e "\t NEXTFLOW:\t\t $nextflow \n\n"
    cat ${confDir}/template.nextflow.config | sed -e "s|pipeInstall|pipeInstall=\"${mypwd}\"|" -e "s|busco4db|busco4db=\"${busco4db}\"|" -e "s|uniprot|uniprot=\"${uniprot}\"|" \
        -e "s|uniname|uniname=\"${uniname}\"|" -e "s|pfloc|pfloc=\"${pfloc}\"|" -e "s|pfname|pfname=\"${pfname}\"|" -e "s|Tsql|Tsql=\"${Tsql}\"|" >nextflow.config
    rm .varfile.sh
}
get_var () {
    cd $mypwd
    echo "busco4db=$mypwd/DBs/busco_db/$busna" >>${mypwd}/.varfile.sh
    echo "uniname=$unina" >>${mypwd}/.varfile.sh
    echo "uniprot=$mypwd/DBs/uniprot_db/$unina" >>${mypwd}/.varfile.sh
    echo "pfloc=$mypwd/DBs/hmmerdb/Pfam-A.hmm" >>${mypwd}/.varfile.sh
    echo "pfname=Pfam-A.hmm" >>${mypwd}/.varfile.sh
    echo "nextflow=$mypwd/nextflow" >>${mypwd}/.varfile.sh
    echo "Tsql=$mypwd/DBs/sqlite_db/*.sqlite" >>${mypwd}/.varfile.sh
    echo "unpdate=\"$( if [ -f ${mypwd}/DBs/uniprot_db/.lastrun.txt ];then cat ${mypwd}/DBs/uniprot_db/.lastrun.txt;else echo "N/A";fi )\"" >>${mypwd}/.varfile.sh
    echo "pfdate=\"$( if [ -f ${mypwd}/DBs/hmmerdb/.lastrun.txt ];then cat ${mypwd}/DBs/hmmerdb/.lastrun.txt;else echo "N/A";fi )\"" >>${mypwd}/.varfile.sh
    echo "dbdate=\"$( if [ -f ${mypwd}/DBs/sqlite_db/.lastrun.txt ];then cat ${mypwd}/DBs/sqlite_db/.lastrun.txt;else echo "N/A";fi )\"" >>${mypwd}/.varfile.sh
    echo "tenv=$( conda info --json | sed -n '/\"envs\":/,/\],/p' | grep -w "TransPi\"" | tr -d "," | tr -d " " )" >>${mypwd}/.varfile.sh
    echo "cenv=$( conda info --json | sed -n '/\"envs\":/,/\],/p' | grep "busco4" | tr -d "," | tr -d " " )" >>${mypwd}/.varfile.sh
    vpwd=$mypwd
    echo "mypwd=$mypwd" >>${vpwd}/.varfile.sh
    source .varfile.sh
    echo -e "\n\t -- INFO to use in TransPi --\n"
    echo -e "\t Installation PATH:\t $mypwd"
    echo -e "\t BUSCO V4 database:\t $busco4db"
    echo -e "\t UNIPROT database:\t $uniprot"
    echo -e "\t UNIPROT last update:\t $unpdate"
    echo -e "\t PFAM files:\t\t $pfloc"
    echo -e "\t PFAM last update:\t $pfdate"
    echo -e "\t SQL DB last update: \t $dbdate"
    echo -e "\t NEXTFLOW:\t\t $nextflow \n\n"
    cat ${confDir}/template.nextflow.config | sed -e "s|pipeInstall|pipeInstall=\"${mypwd}\"|" -e "s|busco4db|busco4db=\"${busco4db}\"|" -e "s|uniprot|uniprot=\"${uniprot}\"|" \
        -e "s|uniname|uniname=\"${uniname}\"|" -e "s|pfloc|pfloc=\"${pfloc}\"|" -e "s|pfname|pfname=\"${pfname}\"|" -e "s|Tsql|Tsql=\"${Tsql}\"|" \
        -e "s|myCondaInstall=\"\"|myCondaInstall=\"${tenv}\"|" -e "s|cenv=\"\"|cenv=\"${cenv}\"|" >nextflow.config
    rm .varfile.sh
}
get_var_user() {
    cd $mypwd
    echo "busco4db=${busco4db}" >>${mypwd}/.varfile.sh
    echo "uniname=${uniname}" >>${mypwd}/.varfile.sh
    echo "uniprot=${uniprot}" >>${mypwd}/.varfile.sh
    echo "pfloc=${pfloc}" >>${mypwd}/.varfile.sh
    echo "pfname=${pfname}" >>${mypwd}/.varfile.sh
    echo "nextflow=$mypwd/nextflow" >>${mypwd}/.varfile.sh
    echo "Tsql=${Tsql}" >>${mypwd}/.varfile.sh
    vpwd=$mypwd
    echo "mypwd=$mypwd" >>${vpwd}/.varfile.sh
    source .varfile.sh
    echo -e "\n\t -- INFO to use in TransPi --\n"
    echo -e "\t Installation PATH:\t $mypwd"
    echo -e "\t Using your DBs\t\t"
    echo -e "\t BUSCO V4 database:\t $busco4db"
    echo -e "\t UNIPROT database:\t $uniprot"
    echo -e "\t PFAM files:\t\t $pfloc"
    echo -e "\t NEXTFLOW:\t\t $nextflow \n\n"
    cat ${confDir}/template.nextflow.config | sed -e "s|pipeInstall|pipeInstall=\"${mypwd}\"|" -e "s|busco4db|busco4db=\"${busco4db}\"|" -e "s|uniprot|uniprot=\"${uniprot}\"|" \
        -e "s|uniname|uniname=\"${uniname}\"|" -e "s|pfloc|pfloc=\"${pfloc}\"|" -e "s|pfname|pfname=\"${pfname}\"|" -e "s|Tsql|Tsql=\"${Tsql}\"|" \
        -e "s|myCondaInstall=\"\"|myCondaInstall=\"${tenv}\"|" -e "s|cenv=\"\"|cenv=\"${cenv}\"|" >nextflow.config
    rm .varfile.sh
}
container_pipeline_setup() {
    if [ "${userVar}" == "y" ];then
        nextflow_c
        evi_c
        echo -e "\n\t -- If no \"ERROR\" was found and all the neccesary databases are installed proceed to run TransPi -- \n"
        get_var_user
    else
        echo -e "\n\t -- Installing databases only -- \n"
        dir_c
        bus_c
        uniprot_c
        nextflow_c
        evi_c
        buildsql_c
        trisql_container
        pfam_c
        echo -e "\n\t -- If no \"ERROR\" was found and all the neccesary databases are installed proceed to run TransPi -- \n"
        get_var_container
    fi
}
conda_pipeline_setup() {
    if [ "${userVar}" == "y" ];then
        echo -e "\n\t -- Installing conda --\n"
        conda_only
        nextflow_c
        evi_c
        echo -e "\n\t -- If no \"ERROR\" was found and all the neccesary databases are installed proceed to run TransPi -- \n"
        get_var_user
    else
        echo -e "\n\t -- Installing conda and the databases -- \n"
        conda_only
        dir_c
        bus_c
        uniprot_c
        nextflow_c
        evi_c
        buildsql_c
        trisql_c
        pfam_c
        bus4
        echo -e "\n\t -- If no \"ERROR\" was found and all the neccesary databases are installed proceed to run TransPi -- \n"
        get_var
    fi
}
user_buscoDBv4(){
    echo -e "\n\t -- PATH where to locate your BUSCO v4 file -- "
    echo -e "\n\t -- Example: /home/ubuntu/myDB/metazoa_odb10 -- "
    echo -e -n "\n\t -- Provide the PATH where to locate your BUSCO v4 file: "
    read -e ans
    if [ -d ${ans} ];then
        echo -e "\n\t -- File ${ans} found -- \n"
        export busco4db=${ans}
    elif [ ! -d ${ans} ];then
        echo -e "\n\t\e[31m -- File ${ans} not found -- \e[39m\n"
        user_buscoDBv4
    fi
}
user_uniDB(){
    echo -e "\n\t -- PATH where to locate your UNIPROT file -- "
    echo -e "\n\t -- Example: /home/ubuntu/myDB/uniprot_proteins.fasta -- "
    echo -e -n "\n\t -- Provide the PATH where to locate your UNIPROT file: "
    read -e ans
    if [ -f ${ans} ];then
        echo -e "\n\t -- File ${ans} found -- \n"
        export uniprot=${ans}
        export uniname=$( basename ${ans} )
    elif [ ! -f ${ans} ];then
        echo -e "\n\t\e[31m -- File ${ans} not found -- \e[39m\n"
        user_uniDB
    fi
}
user_pfDB(){
    echo -e "\n\t -- PATH where to locate your PFAM file -- "
    echo -e "\n\t -- Example: /home/ubuntu/myDB/Pfam-A.hmm -- "
    echo -e -n "\n\t -- Provide the PATH where to locate your PFAM file: "
    read -e ans
    if [ -f ${ans} ];then
        echo -e "\n\t -- File ${ans} found -- \n"
        export pfloc=${ans}
        export pfname=$( basename ${ans} )
    elif [ ! -f ${ans} ];then
        echo -e "\n\t\e[31m -- File ${ans} not found -- \e[39m\n"
        user_pfDB
    fi
}
user_sqlDB(){
    echo -e "\n\t -- PATH where to locate your SQL file -- "
    echo -e "\n\t -- Example: /home/ubuntu/myDB/Trinotate.sqlite -- "
    echo -e -n "\n\t -- Provide the PATH where to locate your SQL file: "
    read -e ans
    if [ -f ${ans} ];then
        echo -e "\n\t -- File ${ans} found -- \n"
        export Tsql=${ans}
    elif [ ! -f ${ans} ];then
        echo -e "\n\t\e[31m -- File ${ans} not found -- \e[39m\n"
        user_sqlDB
    fi
}
userDBs(){
    user_buscoDBv4
    user_uniDB
    user_pfDB
    user_sqlDB
    userVar=y
}
dbs(){
    echo -e "\n\t -- Either TransPi install the DBs for you or you provide the PATH of the DBs -- \n"
    echo -e -n "\t Do you want TransPi to handle the DBs installation? (y,n,exit): "
    read ans
    case $ans in
        [yY] | [yY][eE][sS])
            echo -e "\n\n\t -- TransPi will handle the installation -- \n"
        ;;
        [nN] | [nN][oO])
            echo -e "\n\n\t -- Using your DBs -- \n"
            userDBs
        ;;
        exit)
            echo -e "\n\n\t -- Exiting --\n"
        ;;
        *)
            echo -e "\n\n\t\e[31m -- Yes or No answer not specified. Try again --\e[39m\n"
            dbs
        ;;
    esac
    if [ "$1" == "1" ];then
        conda_pipeline_setup
    elif [ "$1" == "2" ];then
        container_pipeline_setup
    fi
}
message(){
    echo "

    #########################################################################################
    #                                                                                       #
    #                             TransPi precheck script                                   #
    #                                                                                       #
    #   Options available:                                                                  #
    #                                                                                       #
    #        1- Install conda (if neccesary) and DBs                                        #
    #                                                                                       #
    #               Runs of TransPi using individual conda enviroments per process          #
    #                                                                                       #
    #        2- Install DBs for containers (docker or singularity)                          #
    #                                                                                       #
    #               Runs of TransPi with individual containers per process                  #
    #                                                                                       #
    #        3- Install DBs for conda enviroments                                           #
    #                                                                                       #
    #               Runs of TransPi with individual enviroments per process                 #
    #                                                                                       #
    #        4- Update DBs                                                                  #
    #                                                                                       #
    #               SwissProt, PFAM, SQL DB used for annotation (requires conda)            #
    #                                                                                       #
    #        5- Exit                                                                        #
    #                                                                                       #
    #########################################################################################

    "
}
moption(){
    echo -e -n "\t Which option you want? "
    read ans
    case $ans in
        1 | 2 | 3)
            dbs $ans
        ;;
        4)
            echo -e "\n\t -- Updating DBs -- \n"
            downd
        ;;
        5)
            echo -e "\n\t -- Exit -- \n"
            exit 0
        ;;
        *)
            echo -e "\n\t\e[31m -- Wrong option. Try again --\e[39m\n"
            moption
        ;;
    esac
}
main(){
    if [ "$mypwd" == "" ] || [ "$mypwd" == "-h" ] || [ "$mypwd" == "-help" ] || [ "$mypwd" == "--help" ];then
        echo -e "\n\t Script for checking the requirements of TransPi \n"
        echo -e "\t Usage:\n\n\t\t bash precheck_TransPi.sh WORK_PATH \n"
        echo -e "\t\t\t WORK_PATH = PATH to download requirements and databases used by TransPi \n\n\t\t\t Example: /home/bioinf/run/ \n"
        exit 0
    elif [ ! -d "$mypwd" ];then
        echo -e "\n\t -- Directory "${mypwd}" is not found -- \n"
        echo -e "\n\t -- Creating "${mypwd}" -- \n"
        mkdir -p ${mypwd}
        if [ -d "$mypwd" ];then
            echo -e "\n\t -- Directory created succesfully -- \n"
            main
        else
            echo -e "\n\t -- Please provide a valid PATH to run TransPi -- \n"
            exit 0
        fi
    elif [ -d "$mypwd" ];then
        if [ ${mypwd} == "." ];then
            mypwd=$(pwd)
            confDir=$(pwd)
        elif [ ${mypwd} == $(pwd) ]; then
            confDir=$(pwd)
        else
            cd ${mypwd} && mypwd=$(pwd) && cd -
            confDir=$( dirname ${BASH_SOURCE} )
            if [ ${confDir} == "." ];then
                 confDir=$(pwd)
            else
                cd ${confDir} && confDir=$(pwd)
            fi
        fi
        message
        moption
    fi
}
main
