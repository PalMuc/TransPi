#!/usr/bin/env bash -e
export mypwd="$1"
read_c() {
    #Check reads
    cd $mypwd
    if [ ! -d reads/ ];then
        echo -e "\n\t\e[31m -- ERROR: Directory \"reads\" not present. Please create a directory \"reads\" and put the reads files there --\e[39m\n"
        exit 0
    elif [ -d reads/ ];then
        ls -1 reads/*.gz 2>&1 | head -n 1 >.readlist.txt
        if [ `cat .readlist.txt | grep -c "ls:\ cannot"` -eq 1 ];then
            echo -e "\n\t\e[31m -- ERROR: Directory \"reads\" is present but empty. Please copy your reads files --\e[39m\n"
            echo -e "\n\t -- Example: IndA_R1.fastq.gz IndA_R2.fastq.gz -- \n"
            exit 0
        else
            ls -1 reads/*.gz >.readlist.txt
            if [ $(( `cat .readlist.txt | wc -l` % 2 )) -eq 0 ];then
                echo -e "\n\t -- Reads found in $( pwd )/reads/ and in pairs -- \n"
                nind=$( cat .readlist.txt | cut -f 2 -d "/" | sed 's/R\{1,2\}.*//g' | sort -u | wc -l )
                echo -e "\n\t -- Number of samples: $nind -- \n"
            elif [ $(( `cat .readlist.txt | wc -l` % 2 )) -eq 1 ];then
                echo -e "\n\t\e[31m -- ERROR: Reads found in $( pwd )/reads/ but not in pairs. Make sure you have an R1 and R2 for each sample --\e[39m\n"
                echo -e "\n\t -- Example: IndA_R1.fastq.gz IndA_R2.fastq.gz -- \n"
                nind=$( cat .readlist.txt | cut -f 2 -d "/" | sed 's/R\{1,2\}.*//g' | sort -u | wc -l )
                echo -e "\n\t\e[31m -- Number of samples: $nind --\e[39m\n"
                exit 0
            fi
        fi
        rm .readlist.txt
    fi
}
os_c() {
    OS="$(uname)"
    if [ "$OS" == "Linux" ]; then
        echo -e "\n\t -- Downloading Linux Anaconda3 installation -- \n"
        curl -o Anaconda3-2020.02-Linux-x86_64.sh https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh
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
conda_c() {
    source_c
    cd $mypwd
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
                if [ -f transpi_env.yml ];then
                    echo -e "\n\t -- TransPi environment file found. Creating environment... --\n"
                    conda env create -f transpi_env.yml
                else
                    echo -e "\n\t\e[31m -- ERROR: TransPi environment file not found (transpi_env.yml). Please check requirements and rerun the pre-check --\e[39m\n"
                    exit 0
                fi
            elif [ "$check_env" -eq 1 ];then
                echo -e "\n\t -- TransPi environment is installed and ready to be used --\n"
            fi
        fi
    else
        echo -e "\n\t -- Conda is not intalled. Please install Anaconda (https://www.anaconda.com) and rerun this script --\n"
        echo -e -n "\n\t    Do you want to install Anaconda? (y,n,exit): "
        read ans
        case $ans in
            [yY] | [yY][eE][sS])
                os_c
                echo -e "\n\t -- Starting Anaconda installation -- \n"
                bash Anaconda3-20*.sh
                echo -e "\n\t -- Installation done -- \n"
                rm Anaconda3-20*.sh
                source_c
                if [ -f transpi_env.yml ];then
                    echo -e "\n\t -- TransPi environment file found. Creating environment... --\n"
                    conda env create -f transpi_env.yml
                else
                    echo -e "\n\t\e[31m -- ERROR: TransPi environment file not found (transpi_env.yml). Please check requirements and rerun the pre-check --\e[39m\n"
                    exit 0
                fi
            ;;
            [nN] | [nN][oO])
                echo -e "\n\t\e[31m -- ERROR: Download and Install Anaconda. Then rerun the pre-check  --\e[39m\n"
                exit 0
            ;;
            exit)
	           echo -e "\n\t -- Exiting -- \n"
               exit 0
            ;;
            *)
                echo -e "\n\n\t\e[31m -- Yes or No answer not specified. Try again --\e[39m\n"
	            conda_c
            ;;
        esac
    fi
}
dir_c () {
    if [ ! -d scripts/ ];then
        mkdir scripts
    fi
    if [ ! -d DBs ];then
        mkdir DBs
    fi
}
#temporary function for busco V3
busv3_get () {
    v3name=$1
    if [ `cat ${confDir}/conf/busV3list.txt | grep "${v3name}" | wc -l` -eq 1 ];then
        if [ -d ${v3name}_odb9 ];then
            export busnaV3=${v3name}_odb9
        else
            tname=$( cat ${confDir}/conf/busV3list.txt | grep "${v3name}" )
            wget $tname
            tar -xf ${v3name}_odb9.tar.gz
            export busnaV3=${v3name}_odb9
            rm ${v3name}_odb9.tar.gz
        fi
    else
        echo -e "\n\t -- No BUSCO V3 available for ${v3name} --\n"
        exit 0
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
        #get buscov3
        busv3_get $dname
        if [ -d ${dname}_odb10 ];then
            export busna=${dname}_odb10
        fi
    elif [ -d DBs/busco_db/ ];then
        cd DBs/busco_db
        bname=$( echo $name | cut -f 1 -d "_" )
        dname=$( cat ${confDir}/conf/busV4list.txt | grep "${bname};" | cut -f 1 -d ";" | tr [A-Z] [a-z] )
        if [ -d ${dname}_odb10 ];then
            #get buscov3
            busv3_get $dname
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
            #get buscov3
            busv3_get $dname
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
        echo -e "\n\t\e[31m -- ERROR: Please make sure that file \"busV4list.txt\" is available. Please check requirements and rerun the pre-check --\e[39m\n\n"
	    exit 0
    fi
}
uni_c () {
    PS3="
    Please select UNIPROT database to use: "
    select var in `ls *`;do
        if [ "$var" != "" ];then
            echo -e "\n\t -- UNIPROT database selected: \"$var\" --\n"
            export unina=${var}
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
            echo -e "\n\t\e[31m -- ERROR: Please uncompress the file(s) and rerun the pre-check  --\e[39m\n"
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
uniprot_meta () {
    myuni=$( pwd )
    echo -e "\n\t -- TransPi uses a custom protein database (one of many) from UNIPROT for the annotation -- \n"
    echo -e -n "\n\t    Do you want to download the current metazoan proteins from UNIPROT? (y,n,skip,exit): "
    read ans
    case $ans in
        [yY] | [yY][eE][sS])
            echo -e "\n\n\t -- Downloading current metazoan protein dataset from UNIPROT -- \n"
            echo -e "\n\t -- This could take a couple of minutes depending on connection. Please wait -- \n"
            curl -o uniprot_metazoa_33208.fasta.gz "https://www.uniprot.org/uniprot/?query=taxonomy:33208&format=fasta&compress=yes&include=no"
            echo -e "\n\t -- Uncompressing uniprot_metazoa_33208.fasta.gz ... -- \n"
            gunzip uniprot_metazoa_33208.fasta.gz
            date -u >.lastrun.txt
            uni_c
        ;;
        [nN] | [nN][oO])
            echo -e "\n\t\e[31m -- ERROR: Please download your desire UNIPROT database and save it at \"$myuni\". Rerun the pre-check  --\e[39m\n"
            exit 0
        ;;
        skip)
            echo -e "\n\t -- Skipping download of UniProt metazoan proteins. Remember to add the proteins before running TransPi  -- \n"
            exit 0
        ;;
        exit)
            echo -e "\n\t -- Exiting -- \n"
            exit 0
        ;;
        *)
            echo -e "\n\n\t\e[31m -- Yes or No answer not specified. Try again --\e[39m\n"
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
        myfasta=$( ls -1 * | egrep '.fasta|.fa' | wc -l )
        if [ $myfasta -eq 0 ];then
            myfastagz=$( ls -1 * | egrep '.fasta.gz|.fa.gz' | wc -l )
            if [ $myfastagz -eq 0 ];then
                echo -e "\n\t -- Directory \"$myuni\" is empty --\n"
                uniprot_meta
            else
                echo -e "\n\t\e[31m -- Directory \"$myuni\" is available but UNIPROT database is compressed --\e[39m\n"
                unicomp_c
                uni_c
            fi
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
                    echo -e "\n\t\e[31m -- ERROR: Download and Install Nextflow. Then rerun the pre-check  --\e[39m\n"
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
	cd $mypwd
    check_evi=$( command -v tr2aacds.pl | wc -l )
    if [ $check_evi -eq 0 ];then
        if [ ! -d scripts/evigene/ ];then
        echo -e "\n\t -- EvidentialGene is not installed -- \n"
        echo -e "\n\t -- If you will use TransPi container (-profile TransPiContainer) you do not need it."
        echo -e "\t -- Otherwise, please answer yes to the following question -- \n"
        echo -e -n "\n\t    Do you want to install EvidentialGene? (y or n): "
        read ans
        case $ans in
            [yY] | [yY][eE][sS])
                cd scripts
                echo -e "\n\t -- Downloading EvidentialGene -- \n"
                wget http://arthropods.eugenes.org/EvidentialGene/other/evigene_older/evigene19may14.tar
                tar -xf evigene19may14.tar
                rm evigene19may14.tar
                echo -e "\n\t -- Done with EvidentialGene -- \n"
            ;;
            [nN] | [nN][oO])
                echo -e "\n\t\e[31m -- ERROR: Download and Install EvidentialGene. Then rerun the pre-check  --\e[39m\n"
                exit 0
            ;;
            *)
                echo -e "\n\n\t\e[31m -- Yes or No answer not specified. Try again --\e[39m\n"
                evi_c
            ;;
        esac
        else
            echo -e "\n\t -- EvidentialGene directory was found at ${mypwd}/scripts (local installation) -- \n"
        fi
    elif [ $check_evi -eq 1 ];then
        echo -e "\n\t -- EvidentialGene is already installed -- \n"
    fi
}
trisql_container () {
    if [ ! -e *.sqlite ];then
        echo -e "\n\t -- Custom sqlite database for Trinotate is not installed -- \n"
        echo -e -n "\n\t    Do you want to install the custom sqlite database? (y or n): "
        read ans
        case $ans in
            [yY] | [yY][eE][sS])
                echo -e "\n\t -- This could take a couple of minutes depending on connection. Please wait -- \n"
                rm -rf *
                wget https://github.com/Trinotate/Trinotate/archive/Trinotate-v3.2.1.tar.gz
                tar -xf Trinotate-v3.2.1.tar.gz
                mv Trinotate-Trinotate-v3.2.1/ Trinotate_build_scripts/
                ./Trinotate_build_scripts/admin/Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate
                rm uniprot_sprot.dat.gz Pfam-A.hmm.gz
                date -u >.lastrun.txt
            ;;
            [nN] | [nN][oO])
                echo -e "\n\t\e[31m -- ERROR: Generate the custom trinotate sqlite database at "${mypwd}/DBs/sqlite_db". Then rerun the pre-check  --\e[39m\n"
                exit 0
            ;;
            *)
                echo -e "\n\n\t\e[31m -- Yes or No answer not specified. Try again --\e[39m\n"
                trisql_container
            ;;
        esac
    elif [ -e *.sqlite ];then
        echo -e "\n\t -- Custom sqlite database for Trinotate found at "${mypwd}/DBs/sqlite_db" -- \n"
        DB=$( cat .lastrun.txt )
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
                echo -e "\n\t\e[31m -- ERROR: Generate the custom trinotate sqlite database at "${mypwd}/DBs/sqlite_db". Then rerun the pre-check  --\e[39m\n"
                exit 0
            ;;
            *)
                echo -e "\n\n\t\e[31m -- Yes or No answer not specified. Try again --\e[39m\n"
                trisql_c
            ;;
        esac
    elif [ -e *.sqlite ];then
        echo -e "\n\t -- Custom sqlite database for Trinotate found at "${mypwd}/DBs/sqlite_db" -- \n"
        DB=$( cat .lastrun.txt )
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
            DB=$( cat .lastrun.txt )
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
get_var_container () {
    cd $mypwd
    echo "busco3db=$mypwd/DBs/busco_db/$busnaV3" >${mypwd}/.varfile.sh
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
    echo -e "\t Pipeline PATH:\t\t $mypwd"
    echo -e "\t BUSCO V4 database:\t $busco4db"
    echo -e "\t UNIPROT database:\t $uniprot"
    echo -e "\t UNIPROT last update:\t $unpdate"
    echo -e "\t PFAM files:\t\t $pfloc"
    echo -e "\t PFAM last update:\t $pfdate"
    echo -e "\t SQL DB last update: \t $dbdate"
    echo -e "\t NEXTFLOW:\t\t $nextflow \n\n"
    cat template.nextflow.config | sed -e "s|pipeInstall|pipeInstall=\"${mypwd}\"|" -e "s|busco4db|busco4db=\"${busco4db}\"|" -e "s|uniprot|uniprot=\"${uniprot}\"|" \
        -e "s|uniname|uniname=\"${uniname}\"|" -e "s|pfloc|pfloc=\"${pfloc}\"|" -e "s|pfname|pfname=\"${pfname}\"|" -e "s|Tsql|Tsql=\"${Tsql}\"|" \
        -e "s|busco3db|busco3db=\"${busco3db}\"|" >nextflow.config
    rm .varfile.sh
}
get_var () {
    cd $mypwd
    echo "busco3db=$mypwd/DBs/busco_db/$busnaV3" >${mypwd}/.varfile.sh
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
    echo "tenv=$( conda info --json | sed -n '/\"envs\":/,/\],/p' | grep "TransPi" | tr -d "," | tr -d " " )" >>${mypwd}/.varfile.sh
    echo "cenv=$( conda info --json | sed -n '/\"envs\":/,/\],/p' | grep "busco4" | tr -d "," | tr -d " " )" >>${mypwd}/.varfile.sh
    vpwd=$mypwd
    echo "mypwd=$mypwd" >>${vpwd}/.varfile.sh
    source .varfile.sh
    echo -e "\n\t -- INFO to use in TransPi --\n"
    echo -e "\t Pipeline PATH:\t\t $mypwd"
    echo -e "\t BUSCO V4 database:\t $busco4db"
    echo -e "\t UNIPROT database:\t $uniprot"
    echo -e "\t UNIPROT last update:\t $unpdate"
    echo -e "\t PFAM files:\t\t $pfloc"
    echo -e "\t PFAM last update:\t $pfdate"
    echo -e "\t SQL DB last update: \t $dbdate"
    echo -e "\t NEXTFLOW:\t\t $nextflow \n\n"
    cat template.nextflow.config | sed -e "s|pipeInstall|pipeInstall=\"${mypwd}\"|" -e "s|busco4db|busco4db=\"${busco4db}\"|" -e "s|uniprot|uniprot=\"${uniprot}\"|" \
        -e "s|uniname|uniname=\"${uniname}\"|" -e "s|pfloc|pfloc=\"${pfloc}\"|" -e "s|pfname|pfname=\"${pfname}\"|" -e "s|Tsql|Tsql=\"${Tsql}\"|" \
        -e "s|busco3db|busco3db=\"${busco3db}\"|" -e "s|myCondaInstall=\"\"|myCondaInstall=\"${tenv}\"|" -e "s|cenv=\"\"|cenv=\"${cenv}\"|" >nextflow.config
    rm .varfile.sh
}
pipeline_setup() {
    echo -e "\n\t -- Conda or containers (e.g. docker, singularity) -- \n"
    echo -e -n "\n\t    Do you plan to run TransPi with containers (i.e. Singularity or Docker)? (y or n): "
    read ans
    case $ans in
        [yY] | [yY][eE][sS])
            echo -e "\n\n\t -- Installing databases only -- \n"
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
        ;;
        [nN] | [nN][oO])
            echo -e "\n\n\t -- Installing conda and the databases -- \n"
            conda_c
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

        ;;
        *)
            echo -e "\n\n\t\e[31m -- Yes or No answer not specified. Try again --\e[39m\n"
            pipeline_setup
        ;;
    esac
}
#Main
if [ "$mypwd" == "" ] || [ "$mypwd" == "-h" ] || [ "$mypwd" == "-help" ] || [ "$mypwd" == "--help" ];then
    echo -e "\n\t Script for checking the requirements of TransPi \n"
    echo -e "\t Usage:\n\n\t\t bash precheck_TransPi.sh WORK_PATH \n"
    echo -e "\t\t\t WORK_PATH = PATH to download requirements and databases used by TransPi \n\n\t\t\t Example: /home/bioinf/run/ \n"
    exit 0
elif [ ! -d "$mypwd" ];then
    echo -e "\n\t -- Please provide a valid PATH to run TransPi -- \n"
    exit 0
elif [ -d "$mypwd" ];then
    if [ ${mypwd} == "." ];then
        mypwd=$(pwd)
        confDir=$(pwd)
        cd $mypwd
    elif [ ${mypwd} == $(pwd) ]; then
        confDir=$(pwd)
        cd $mypwd
    else
        confDir=$( dirname ${BASH_SOURCE} )
        if [ ${confDir} == "." ];then
             confDir=$(pwd)
        fi
    fi
    pipeline_setup
fi
