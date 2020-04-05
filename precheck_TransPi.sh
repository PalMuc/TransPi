#!/usr/bin/env bash
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
    if [ -f /etc/os-release ];then
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
    source ~/.bashrc
    cd $mypwd
    #Check conda and environment
    check_conda=$( command -v conda )
    if [ "$check_conda" != "" ];then #&& [ "$ver" -gt "45" ];then
        echo -e "\n\t -- Conda seems to be installed in your system environment --\n"
        ver=$( conda -V | cut -f 2 -d " " | cut -f 1,2 -d "." | tr -d "." )
        if [ "$ver" -gt 45 ];then
            echo -e "\n\t -- Conda is installed (v4.5 or higher). Checking environment... --\n"
            #Check environment
            check_env=$( conda info -e | grep -c "TransPi" )
	    if [ "$check_env" -eq 0 ];then
                echo -e "\n\t -- TransPi environment has not been created. Checking environment file... --\n"
                if [ -f transpi_env.yml ];then
                    echo -e "\n\t -- TransPi environment file found. Creating environment... --\n"
                    conda env create -f transpi_env.yml
                else
                    echo -e "\n\t\e[31m -- ERROR: TransPi environment file not found \(transpi_env.yml\). Please check requirements and rerun the pre-check --\e[39m\n"
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
    if [ `cat ${mypwd}/conf/busV3list.txt | grep "${v3name}" | wc -l` -eq 1 ];then
        if [ -d ${v3name}_odb9 ];then
            export busnaV3=${v3name}_odb9
        else
            tname=$( cat ${mypwd}/conf/busV3list.txt | grep "${v3name}" )
            wget $tname
            tar -xvf ${v3name}_odb9.tar.gz
            export busnaV3=${v3name}_odb9
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
        if [ `cat ${mypwd}/conf/busV4list.txt | grep "${bname};" | wc -l` -eq 1 ];then
            echo -e "\n\t -- Downloading BUSCO V4 \"$name\" database --\n";wait
            wname=$( cat ${mypwd}/conf/busV4list.txt | grep "${bname};" | cut -f 2 -d ";" )
            wget $wname
            echo -e "\n\t -- Preparing files ... --\n";wait
            tname=$( cat ${mypwd}/conf/busV4list.txt | grep "${bname};" | cut -f 1 -d ";" | tr [A-Z] [a-z] )
            tar -xvf ${tname}*.tar.gz
            rm ${tname}*.tar.gz
            echo -e "\n\t -- DONE with BUSCO V4 database --\n";wait
        fi
        dname=$( cat ${mypwd}/conf/busV4list.txt | grep "${bname};" | cut -f 1 -d ";" | tr [A-Z] [a-z] )
        #get buscov3
        busv3_get $dname
        if [ -d ${dname}_odb10 ];then
            export busna=${dname}_odb10
        fi
    elif [ -d DBs/busco_db/ ];then
        cd DBs/busco_db
        bname=$( echo $name | cut -f 1 -d "_" )
        dname=$( cat ${mypwd}/conf/busV4list.txt | grep "${bname};" | cut -f 1 -d ";" | tr [A-Z] [a-z] )
        if [ -d ${dname}_odb10 ];then
            #get buscov3
            busv3_get $dname
            echo -e "\n\t -- BUSCO V4 \"$name\" database found -- \n"
            export busna=${dname}_odb10
        else
            bname=$( echo $name | cut -f 1 -d "_" )
            if [ `cat ${mypwd}/conf/busV4list.txt | grep "${bname};" | wc -l` -eq 1 ];then
                echo -e "\n\t -- Downloading BUSCO V4 \"$name\" database --\n";wait
                wname=$( cat ${mypwd}/conf/busV4list.txt | grep "${bname};" | cut -f 2 -d ";" )
                wget $wname
                echo -e "\n\t -- Preparing files ... --\n";wait
                tname=$( cat ${mypwd}/conf/busV4list.txt | grep "${bname};" | cut -f 1 -d ";" | tr [A-Z] [a-z] )
                tar -xvf ${tname}*.tar.gz
                rm ${tname}*.tar.gz
                echo -e "\n\t -- DONE with BUSCO V4 database --\n";wait
            fi
            dname=$( cat ${mypwd}/conf/busV4list.txt | grep "${bname};" | cut -f 1 -d ";" | tr [A-Z] [a-z] )
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
    if [ -f ${mypwd}/conf/busV4list.txt ];then
    select var in `cat ${mypwd}/conf/busV4list.txt | grep "###" | tr -d "#"`;do
    case $var in
        BACTERIA)
            echo -e "\n\t You selected BACTERIA. Which specific database? \n"
            PS3="
	    Please select database: "
            select var1 in `cat ${mypwd}/conf/busV4list.txt | sed -n "/##BACTERIA/,/#MAIN/p" | grep -v "##" | tr -d "#"`;do
    	    case $var1 in
    	        MAIN_MENU)
                    bus_c
                ;;
                *)
                if [ "$var1" != "" ];then
                    if [ `cat ${mypwd}/conf/busV4list.txt | grep -c "$var1"` -ge 1 ];then
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
            select var1 in `cat ${mypwd}/conf/busV4list.txt | sed -n "/##EUKARYOTA/,/#MAIN/p" | grep -v "##" | tr -d "#"`;do
        	case $var1 in
        	    MAIN_MENU)
                    bus_c
                ;;
                Arthropoda_\(Phylum\))
                    select var2 in `cat ${mypwd}/conf/busV4list.txt | sed -n "/##ARTHROPODA/,/#MAIN/p" | grep -v "##" | tr -d "#"`;do
                    case $var2 in
                    MAIN_MENU)
                        bus_c
                    ;;
                    *)
                    if [ "$var2" != "" ];then
                        if [ `cat ${mypwd}/conf/busV4list.txt | grep -c "$var2"` -ge 1 ];then
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
                    select var2 in `cat ${mypwd}/conf/busV4list.txt | sed -n "/##FUNGI/,/#MAIN/p" | grep -v "##" | tr -d "#"`;do
                    case $var2 in
                    MAIN_MENU)
                        bus_c
                    ;;
                    *)
                    if [ "$var2" != "" ];then
                        if [ `cat ${mypwd}/conf/busV4list.txt | grep -c "$var2"` -ge 1 ];then
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
                    select var2 in `cat ${mypwd}/conf/busV4list.txt | sed -n "/##PLANTS/,/#MAIN/p" | grep -v "##" | tr -d "#"`;do
                    case $var2 in
                    MAIN_MENU)
                        bus_c
                    ;;
                    *)
                    if [ "$var2" != "" ];then
                        if [ `cat ${mypwd}/conf/busV4list.txt | grep -c "$var2"` -ge 1 ];then
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
                    select var2 in `cat ${mypwd}/conf/busV4list.txt | sed -n "/##PROTIST/,/#MAIN/p" | grep -v "##" | tr -d "#"`;do
                    case $var2 in
                    MAIN_MENU)
                        bus_c
                    ;;
                    *)
                    if [ "$var2" != "" ];then
                        if [ `cat ${mypwd}/conf/busV4list.txt | grep -c "$var2"` -ge 1 ];then
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
                    select var2 in `cat ${mypwd}/conf/busV4list.txt | sed -n "/##VERTEBRATA/,/#MAIN/p" | grep -v "##" | tr -d "#"`;do
                    case $var2 in
                    MAIN_MENU)
                        bus_c
                    ;;
                    *)
                    if [ "$var2" != "" ];then
                        if [ `cat ${mypwd}/conf/busV4list.txt | grep -c "$var2"` -ge 1 ];then
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
                    if [ `cat ${mypwd}/conf/busV4list.txt | grep -c "$var1"` -ge 1 ];then
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
            select var1 in `cat ${mypwd}/conf/busV4list.txt | sed -n "/##ARCHAEA/,/#MAIN/p" | grep -v "##" | tr -d "#"`;do
            case $var1 in
            	MAIN_MENU)
                    bus_c
                ;;
                *)
                if [ "$var1" != "" ];then
                    if [ `cat ${mypwd}/conf/busV4list.txt | grep -c "$var1"` -ge 1 ];then
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
    select var in `ls *fasta`;do
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
            gunzip *fasta.gz
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
uniprot_c () {
    #Check UNIPROT
    cd $mypwd
    if [ ! -d DBs/uniprot_db/ ];then
        echo -e "\n\t -- Creating directory for the UNIPROT database --\n"
        mkdir -p DBs/uniprot_db/
        cd DBs/uniprot_db/
        myuni=$( pwd )
        echo -e "\n\t -- TransPi uses customs protein databases from UNIPROT for the annotation -- \n"
        echo -e -n "\n\t    Do you want to download the current metazoan proteins from UNIPROT? (y,n,exit): "
        read ans
        case $ans in
            [yY] | [yY][eE][sS])
                echo -e "\n\n\t -- Downloading metazoa protein dataset from UNIPROT -- \n"
                echo -e "\n\t -- This could take a couple of minutes depending on connection. Please wait -- \n"
                curl -o uniprot_metazoa_33208.fasta.gz "https://www.uniprot.org/uniprot/?query=taxonomy:33208&format=fasta&compress=yes&include=no"
                gunzip uniprot_metazoa_33208.fasta.gz
                date -u >.lastrun.txt
                uni_c
            ;;
            [nN] | [nN][oO])
                echo -e "\n\t\e[31m -- ERROR: Please download your desire UNIPROT database and save it at \"$myuni\". rerun the pre-check  --\e[39m\n"
                exit 0
            ;;
            exit)
                echo -e "\n\t -- Exiting -- \n"
                exit 0
            ;;
            *)
                echo -e "\n\n\t\e[31m -- Yes or No answer not specified. Try again --\e[39m\n"
                uniprot_c
            ;;
        esac
    elif [ -d DBs/uniprot_db/];then
        cd DBs/uniprot_db/
        myuni=$( pwd )
        echo -e "\n\t -- UNIPROT database directory found at: $myuni -- \n"
        ls -1 *.fasta 2>&1 | head -n 1 >.unilist.txt
        if [ `cat .unilist.txt | grep -c "ls:\ cannot"` -eq 1 ];then
            ls -1 *.fasta.gz 2>&1 | head -n 1 >.unilist.txt
            if [ `cat .unilist.txt | grep -c "ls:\ cannot"` -eq 1 ];then
                echo -e "\n\t\e[31m -- ERROR: Directory \"$myuni\" is empty. Please download a UNIPROT database and rerun the pre-check --\e[39m\n"
                rm .unilist.txt
                uniprot_c
            else
                echo -e "\n\t\e[31m -- Directory \"$myuni\" is available but UNIPROT database is compressed --\e[39m\n"
                unicomp_c
                uni_c
                rm .unilist.txt
            fi
        else
            echo -e "\n\t -- Here is the list of UNIPROT files found at: $myuni -- \n"
            uni_c
            rm .unilist.txt
        fi
    fi
}
java_c () {
	export NXF_VER=20.01.0-edge && curl -s https://get.nextflow.io | bash 2>.error_nextflow
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
	check_next=$( ./nextflow info | head -n 1 | wc -l )
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
        echo -e -n "\n\t    Do you want to install EvidentialGene? (y or n): "
        read ans
        case $ans in
            [yY] | [yY][eE][sS])
                cd scripts
                echo -e "\n\t -- Downloading EvidentialGene ... -- \n"
                wget http://arthropods.eugenes.org/EvidentialGene/other/evigene_old/evigene_older/evigene19may14.tar
                tar -xf evigene19may14.tar
                rm evigene19may14.tar
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
#change BuildSQL PATH (temporary)
sql_path () {
    a=$( which Build_Trinotate_Boilerplate_SQLite_db.pl )
    sed -i 's/$SPROT_DAT_URL = \"http:/$SPROT_DAT_URL = \"ftp:/' $a
    sed -i 's/$EGGNOG_DAT_URL = \".*/$EGGNOG_DAT_URL = \"http://eggnog5.embl.de/download/eggnog_5.0/e5.og_annotations.tsv\";/' $a
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
                echo -e "\n\t -- This could take a couple of minutes depending on connection. Please wait -- \n"
                if [ ! -f ~/anaconda3/etc/profile.d/conda.sh ];then
                    echo -e -n "\n\t    Provide the full PATH of your Anaconda main installation, not an environment (Examples: /home/bioinf/anaconda3 ,  ~/tools/anaconda3 ,  ~/tools/py3/anaconda3): "
                    read ans
                    source ${ans}/etc/profile.d/conda.sh
                    conda activate TransPi
                    check_sql=$( command -v Build_Trinotate_Boilerplate_SQLite_db.pl | wc -l )
                    if [ $check_sql -eq 0 ];then
                        echo -e "\n\t -- Script "Build_Trinotate_Boilerplate_SQLite_db.pl" from Trinotate cannot be found -- \n"
                        echo -e "\n\t\e[31m -- Verify your conda installation --\e[39m\n"
                        exit 0
                    elif [ $check_sql -eq 1 ];then
                        Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate
                        rm uniprot_sprot.dat.gz
                        date -u >.lastrun.txt
                        sql_path
                    fi
                elif [ -f ~/anaconda3/etc/profile.d/conda.sh ];then
                    source ~/anaconda3/etc/profile.d/conda.sh
                    conda activate TransPi
                    check_sql=$( command -v Build_Trinotate_Boilerplate_SQLite_db.pl | wc -l )
                    if [ $check_sql -eq 0 ];then
                        echo -e "\n\t -- Script "Build_Trinotate_Boilerplate_SQLite_db.pl" from Trinotate cannot be found -- \n"
                        echo -e "\n\t\e[31m -- Verify your conda installation --\e[39m\n"
                        exit 0
                    elif [ $check_sql -eq 1 ];then
                        Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate
                        rm uniprot_sprot.dat.gz
                        date -u >.lastrun.txt
                        sql_path
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
    else
        echo -e "\n\t -- Custom sqlite database for Trinotate found at "${mypwd}/DBs/sqlite_db" -- \n"
        DB=$( cat .lastrun.txt )
        echo -e "\n\t -- Databases (PFAM,SwissProt,EggNOG,GO) last update: ${DB} --\n "
    fi
}
buildsql_c () {
    cd ${mypwd}
    if [ -d DBs/sqlite_db/ ];then
        cd DBs/sqlite_db/
        trisql_c
    else
        mkdir -p DBs/sqlite_db/
        cd DBs/sqlite_db/
        trisql_c
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
cbs_c () {
    #rnammer
    cd cbs-dtu-tools/rnammer/
    name=$( pwd )
    sed -i "s|/home/ubuntu/pipe/rnammer|$name|g" rnammer
    cd ..
    #signalP
    cd signalp-4.1/
    name=$( pwd )
    sed -i "s|/home/ubuntu/pipe/signalp-4.1|$name|g" signalp
    cd ..
}
#later to be removed or to scripts directory
cbs_dtu_c () {
    cd $mypwd
    if [ -f cbs-dtu-tools.tar.gz ] && [ ! -d cbs-dtu-tools/ ];then
        echo -e "\n\t -- Preparing scripts of CBS-DTU -- \n"
        echo -e "\n\t -- Uncompressing files -- \n"
        tar -xvf cbs-dtu-tools.tar.gz
        cbs_c
    elif [ -f cbs-dtu-tools.tar.gz ] && [ -d cbs-dtu-tools/ ];then
        cbs_c
    elif [ ! -f cbs-dtu-tools.tar.gz ] && [ ! -d cbs-dtu-tools/ ];then
        echo -e "\n\t\e[31m -- ERROR: Please make sure the cbs-dtu-tools.tar.gz is available. Then rerun the pre-check  --\e[39m\n"
        exit 0
    fi
}
util_c () {
    source ~/.bashrc
    cpath=$( conda info -e | grep "TransPi" | awk '{print $2}' )
    if [ -f ${cpath}/bin/RnammerTranscriptome.pl ];then
        sed -i "s|RealBin/util|RealBin|g" ${cpath}/bin/RnammerTranscriptome.pl
    else
        echo -e "\n\t -- Cannot find the TransPi environment -- \n"
        echo -e -n "\n\t    Provide the PATH of TransPi environment (Examples: /home/bioinf/anaconda3/envs/TransPi ,  ~/tools/anaconda3/.conda/envs/TransPi): "
        read ans
        sed -i "s|RealBin/util|RealBin|g" ${ans}/bin/RnammerTranscriptome.pl
    fi
}
#temporary for buscoV4
bus_conf () {
    head -n 56 config.ini >.56.txt
    rm config.ini
    #get the .57.txt
    cpath=$( conda info -e | grep "TransPi" | awk '{print $2}' )
    tail -n +45 ${cpath}/config/config.ini >.57.txt
    cat .56.txt .57.txt >config.ini
    rm .56.txt .57.txt
}
bus_dow4 () {
    wget https://gitlab.com/ezlab/busco/-/archive/4.0.5/busco-4.0.5.tar.gz
    tar -xf busco-4.0.5.tar.gz
    cd busco-4.0.5
    python3 setup.py install --user
    cd ..
    cp busco-4.0.5/bin/busco .
    cp busco-4.0.5/config/config.ini .
    rm -rf busco-4.0.5/
}
bus4 () {
    cd $mypwd
    if [ ! -d scripts/ ];then
        mkdir -p scripts/busco4
        cd scripts/busco4
        bus_dow4
        # here modify the config.ini
        bus_conf
    elif [ -d scripts/ ];then
        cd scripts
        if [ ! -d busco4 ];then
            mkdir busco4
            cd busco4
            cpath=$( conda info -e | grep "TransPi" | awk '{print $2}' )
            if [ "$cpath" == "" ];then
                echo -e "\n\t -- Cannot find the TransPi environment -- \n"
                echo -e -n "\n\t    Provide the PATH of TransPi environment ( Examples: /home/bioinf/anaconda3/envs/TransPi ,  ~/tools/anaconda3/.conda/envs/TransPi ): "
                read ans
                conda activate ${ans}
                bus_dow4
                # here modify the config.ini
                bus_conf
            else
                conda activate TransPi
                bus_dow4
                # here modify the config.ini
                bus_conf
            fi
        elif [ -d busco4 ];then
            cd busco4
            if [ -f busco ] && [ -f config.ini ];then
                echo "BUSCO V4 is ready to use"
            else
                rm -rf *
                bus_dow4
                # here modify the config.ini
                bus_conf
            fi
        fi
    fi
}
get_var () {
    cd $mypwd
    #echo "=$mypwd/" >${mypwd}/.varfile.sh
    echo "busco3db=$mypwd/DBs/busco_db/$busnaV3" >${mypwd}/.varfile.sh
    echo "busco4db=$mypwd/DBs/busco_db/$busna" >${mypwd}/.varfile.sh
    echo "uniname=$unina" >>${mypwd}/.varfile.sh
    echo "uniprot=$mypwd/DBs/uniprot_db/$unina" >>${mypwd}/.varfile.sh
    echo "pfloc=$mypwd/DBs/hmmerdb/Pfam-A.hmm" >>${mypwd}/.varfile.sh
    echo "pfname=Pfam-A.hmm" >>${mypwd}/.varfile.sh
    echo "nextflow=$mypwd/nextflow" >>${mypwd}/.varfile.sh
    echo "Tsql=$mypwd/DBs/sqlite_db/*.sqlite" >>${mypwd}/.varfile.sh
    echo "rnam=$mypwd/cbs-dtu-tools/rnammer/rnammer" >>${mypwd}/.varfile.sh
    echo "tmhmm=$mypwd/cbs-dtu-tools/tmhmm-2.0c/bin/tmhmm" >>${mypwd}/.varfile.sh
    echo "signalp=$mypwd/cbs-dtu-tools/signalp-4.1/signalp" >>${mypwd}/.varfile.sh
    #echo "unpdate=$( cat ${mypwd}/uniprot_db/.lastrun.txt )" >>${mypwd}/.varfile.sh
    echo "pfdate=$( cat ${mypwd}/DBs/hmmerdb/.lastrun.txt )" >>${mypwd}/.varfile.sh
    echo "dbdate=$( cat ${mypwd}/DBs/sqlite_db/.lastrun.txt )" >>${mypwd}/.varfile.sh
    vpwd=$mypwd
    echo "mypwd=$mypwd" >>${vpwd}/.varfile.sh
    source .varfile.sh
    echo -e "\n\n\t -- INFO to use in TransPi --\n"
    echo -e "\t Pipeline PATH:\t\t $mypwd"
    echo -e "\t BUSCO V4 database:\t $busco4db"
    echo -e "\t UNIPROT database:\t $uniprot"
    #echo -e "\t UNIPROT last update:\t $unpdate"
    echo -e "\t PFAM files:\t\t $pfloc"
    echo -e "\t PFAM last update:\t $pfdate"
    echo -e "\t SQL database (SwissProt,EggNOG,GO,PFAM) last update: \t $dbdate"
    echo -e "\t NEXTFLOW:\t\t $nextflow \n\n"
    cat template.nextflow.config | sed -e "s|mypwd|mypwd=\"${mypwd}\"|" -e "s|busco4db|busco4db=\"${busco4db}\"|" -e "s|uniprot|uniprot=\"${uniprot}\"|" \
        -e "s|uniname|uniname=\"${uniname}\"|" -e "s|pfloc|pfloc=\"${pfloc}\"|" -e "s|pfname|pfname=\"${pfname}\"|" -e "s|Tsql|Tsql=\"${Tsql}\"|" \
        -e "s|reads=|reads=\"${mypwd}|" -e "s|rnam|rnam=\"${rnam}\"|" -e "s|tmhmm|tmhmm=\"${tmhmm}\"|" -e "s|signalp|signalp=\"${signalp}\"|" \
        -e "s|busco3db|busco3db=\"${busco3db}\"|" >nextflow.config
    rm .varfile.sh
}
#Main
if [ "$mypwd" == "" ] || [ "$mypwd" == "-h" ] || [ "$mypwd" == "-help" ] || [ "$mypwd" == "--help" ];then
    echo -e "\n\t Script for checking the requirenments of TransPi \n"
    echo -e "\t Usage:\n\n\t\t bash precheck_TransPi.sh WORK_PATH \n"
    echo -e "\n\t\t WORK_PATH = PATH to run TransPi and download the requirenments \n\n\t\t Example: /home/bioinf/run/ \n"
    exit 0
elif [ ! -d "$mypwd" ];then
    echo -e "\n\t -- Please provide a valid PATH to run TransPi -- \n"
    exit 0
elif [ -d "$mypwd" ];then
    cd $mypwd
    read_c
    conda_c
    dir_c
    bus_c
    uniprot_c
    nextflow_c
    evi_c
    buildsql_c
    pfam_c
    cbs_dtu_c
    util_c
    bus4
    echo -e "\n\t -- If no \"ERROR\" was found and all the neccesary databases are installed proceed to run TransPi -- \n"
    get_var
fi
