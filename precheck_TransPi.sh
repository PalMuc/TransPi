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
conda_c() {
    #Check conda and environment
    check_conda=$( command -v conda )
    ver=$( conda -V | cut -f 2 -d " " | cut -f 1,2 -d "." | tr -d "." )
    if [ -f "$check_conda" ] && [ "$ver" -gt "45" ];then
        echo -e "\n\t -- Conda is intalled (v4.5 or higher). Checking environment... --\n"
        #Check environment
        check_env=$( conda env list | grep "TransPi" )
        if [ `echo check_env | wc -l` -eq 0 ];then
            echo -e "\n\t -- TransPi environment has not been created. Checking environment file... --\n"
            if [ -f transpi_env.yml ];then
                echo -e "\n\t -- TransPi environment file found. Creating environment... --\n"
                conda env create -f transpi_env.yml
            else
                echo -e "\n\t\e[31m -- ERROR: TransPi environment file not found \(transpi_env.yml\). Please check requirements and rerun the pre-check --\e[39m\n"
                exit 0
            fi
        elif [ `echo check_env | wc -l` -eq 1 ];then
            echo -e "\n\t -- TransPi environment is installed and ready to be used --\n"
        fi
    else
        echo -e "\n\t -- Conda is not intalled. Please install Anaconda or Miniconda (https://www.anaconda.com) and rerun this script --\n"
        echo -e -n "\n\t    Do you want me to try to install Anaconda for you? (y,n,exit): "
        read ans
        case $ans in
            [yY] | [yY][eE][sS])
                echo -e "\n\t -- Downloading Anaconda ... -- \n"
                wget https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh
                echo -e "\n\t -- Starting anaconda installation -- \n"
                bash Anaconda3-2019.10-Linux-x86_64.sh
                echo -e "\n\t -- Installation done -- \n"
                rm Anaconda3-2019.10-Linux-x86_64.sh
                source ~/.bashrc
                if [ -f transpi_env.yml ];then
                    echo -e "\n\t -- TransPi environment file found. Creating environment... --\n"
                    conda env create -f transpi_env.yml
                else
                    echo -e "\n\t\e[31m -- ERROR: TransPi environment file not found (transpi_env.yml). Please check requirements and rerun the pre-check --\e[39m\n"
                    exit 0
                fi
            ;;
            [nN] | [nN][oO])
                echo -e "\n\t\e[31m -- ERROR: Download and Install Anaconda. Then rerun the pre-check again --\e[39m\n"
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
pfam_c() {
    #Check PFAM files
    cd $mypwd
    if [ ! -d hmmerdb/ ];then
        echo -e "\n\t -- Creating directory for the HMMER database --\n"
        mkdir hmmerdb
        cd hmmerdb
        echo -e "\n\t -- Downloading current release of PFAM for the HMMER database --\n"
        wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
        echo -e "-- Preparing files ... --\n"
        gunzip Pfam-A.hmm.gz
    elif [ -d hmmerdb/ ];then
        echo -e "\n\t -- Directory for the HMMER database is present --\n"
        cd hmmerdb
        if [ -f Pfam-A.hmm ];then
            echo -e "\n\t -- Pfam file is present iand ready to be used --\n"
        else
            echo -e "\n\t -- Downloading current release of PFAM for the HMMER database --\n"
            wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
            echo -e "-- Preparing files ... --\n"
            gunzip Pfam-A.hmm.gz
        fi
    fi
}
bus_dow () {
    name=$1
    cd $mypwd
    if [ ! -d busco_db/ ];then
        echo -e "\n\t -- Creating directory for the BUSCO database --\n"
        mkdir busco_db
        cd busco_db
        if [ `echo $name | cut -f 1 -d "/"` != "prerelease" ];then
            echo -e "\n\t -- Downloading BUSCO \"$name\" database --\n";wait
            wname=$( echo "https://busco.ezlab.org/datasets/${name}.tar.gz" )
            wget $wname
            echo -e "\n\t -- Preparing files ... --\n";wait
            tar -xvf ${name}.tar.gz
            rm ${name}.tar.gz
            echo -e "\n\t -- DONE with BUSCO database --\n";wait
            if [ -d ${name} ];then
                export busna=${name}
                echo $busna
            fi
        elif [ `echo $name | cut -f 1 -d "/"` == "prerelease" ];then
            echo -e "\n\t -- BUSCO \"$name\" database not found -- \n"
            echo -e "\n\t -- Downloading BUSCO \"$name\" database --\n";wait
            wname=$( echo "https://busco.ezlab.org/datasets/${name}.tar.gz" )
            wget $wname
            echo -e "\n\t -- Preparing files ... --\n";wait
            name2=$( echo $name | cut -f 2 -d "/" )
            tar -xvf ${name2}.tar.gz
            rm ${name2}.tar.gz
            echo -e "\n\t -- DONE with BUSCO database --\n";wait
            if [ -d ${name} ];then
                export busna=${name2}
            fi
        fi
    elif [ -d busco_db/ ];then
        cd busco_db
        if [ `echo $name | cut -f 1 -d "/"` != "prerelease" ];then
            if [ -d ${name} ];then
                echo -e "\n\t -- BUSCO \"$name\" database found -- \n"
                export busna=${name}
            else
                echo -e "\n\t -- BUSCO \"$name\" database not found -- \n"
                echo -e "\n\t -- Downloading BUSCO \"$name\" database --\n";wait
                wname=$( echo "https://busco.ezlab.org/datasets/${name}.tar.gz" )
                wget $wname
                echo -e "\n\t -- Preparing files ... --\n";wait
                tar -xvf ${name}.tar.gz
                rm ${name}.tar.gz
                echo -e "\n\t -- DONE with BUSCO database --\n";wait
                if [ -d ${name} ];then
                    export busna=${name}
                fi
            fi
        elif [ `echo $name | cut -f 1 -d "/"` == "prerelease" ];then
            name2=$( echo $name | cut -f 2 -d "/" )
            if [ -d ${name2} ];then
                echo -e "\n\t -- BUSCO \"$name\" database found -- \n"
                export busna=${name2}
            else
                echo -e "\n\t -- BUSCO \"$name\" database not found -- \n"
                echo -e "\n\t -- Downloading BUSCO \"$name\" database --\n";wait
                wname=$( echo "https://busco.ezlab.org/datasets/${name}.tar.gz" )
                wget $wname
                echo -e "\n\t -- Preparing files ... --\n";wait
                tar -xvf ${name2}.tar.gz
                rm ${name2}.tar.gz
                echo -e "\n\t -- DONE with BUSCO database --\n";wait
                if [ -d ${name2} ];then
                    export busna=${name2}
                fi
            fi
        fi
    fi
}
bus_c () {
    cd $mypwd
    echo -e "\n\t -- Selecting BUSCO database -- \n"
    PS3="
    Please select one (1-5): "
    if [ -f buslist.txt ];then
    select var in `cat buslist.txt | grep "##" | tr -d "#"`;do
    case $var in
        BACTERIA)
            echo -e "\n\t You selected BACTERIA. Which specific database? \n"
            PS3="
	    Please select database: "
            select var1 in `cat buslist.txt | sed -n "/##BACTERIA/,/#/p" | grep -v "##" | tr -d "#" | cut -f 5- -d "/" | cut -f 1 -d "."`;do
    	    case $var1 in
    	        MAIN_MENU)
                    bus_c
                ;;
                *)
                if [ "$var1" != "" ];then
                    if [ `cat buslist.txt | grep -c "$var1"` -ge 1 ];then
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
            select var1 in `cat buslist.txt | sed -n "/##EUKARYOTA/,/#/p" | grep -v "##" | tr -d "#" | cut -f 5- -d "/" | cut -f 1 -d "."`;do
        	case $var1 in
        	    MAIN_MENU)
                    bus_c
                ;;
                *)
                if [ "$var1" != "" ];then
                    if [ `cat buslist.txt | grep -c "$var1"` -ge 1 ];then
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
        FUNGI)
            echo -e "\n\tYou selected FUNGI. Which specific database? \n"
            PS3="
	    Please select database: "
            select var1 in `cat buslist.txt | sed -n "/##FUNGI/,/#/p" | grep -v "##" | tr -d "#" | cut -f 5- -d "/" | cut -f 1 -d "."`;do
            case $var1 in
            	MAIN_MENU)
                    bus_c
                ;;
                *)
                if [ "$var1" != "" ];then
                    if [ `cat buslist.txt | grep -c "$var1"` -ge 1 ];then
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
        PLANTS)
            echo -e "\n\tYou selected PLANTS. Which specific database? \n"
            PS3="
	    Please select database: "
            select var1 in `cat buslist.txt | sed -n "/##PLANTS/,/#/p" | grep -v "##" | tr -d "#" | cut -f 5- -d "/" | cut -f 1 -d "."`;do
            case $var1 in
                MAIN_MENU)
                    bus_c
                ;;
                *)
                if [ "$var1" != "" ];then
                    if [ `cat buslist.txt | grep -c "$var1"` -ge 1 ];then
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
        echo -e "\n\t\e[31m -- ERROR: Please make sure that file \"buslist.txt\" is available. Please check requirements and rerun the pre-check --\e[39m\n\n"
        #
	#
	#
	# Try to get it automatically from GitHub  ############################################################################################################
	#
	#
	#
	exit 0
    fi
}
uni_c () {
    cd $mypwd
    cd uniprot_db
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
    echo -e -n "\n\t    Do you want me to uncompress the file(s)? (y,n,exit): "
    read ans
    case $ans in
        [yY] | [yY][eE][sS])
            echo -e "\n\n\t -- Uncompressing file(s) ... -- \n"
            gunzip *fasta.gz
        ;;
        [nN] | [nN][oO])
            echo -e "\n\t\e[31m -- ERROR: Please uncompress the file(s) and rerun the pre-check again --\e[39m\n"
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
    if [ ! -d uniprot_db/ ];then
        echo -e "\n\t -- Creating directory for the UNIPROT database --\n"
        mkdir uniprot_db
        cd uniprot_db/
        myuni=$( pwd )
        echo -e "\n\t -- Before running TransPi, please download the desire UNIPROT database to compare and annotate your transcriptome -- \n"
        echo -e "\n\t -- PATH for UNIPROT database at: $myuni -- \n"
        echo -e "\n\t -- Example: \"uniprot-taxonomy_metazoaA33208.fasta\" -- \n"
        echo -e "\n\t\e[31m -- ERROR: Download UNIPROT database at \"$myuni\" and rerun the pre-check --\e[39m\n\n"
        exit 0
    elif [ -d uniprot_db/ ];then
        cd uniprot_db/
        myuni=$( pwd )
        echo -e "\n\t -- UNIPROT database directory found at: $myuni -- \n"
        ls -1 *.fasta 2>&1 | head -n 1 >.unilist.txt
        if [ `cat .unilist.txt | grep -c "ls:\ cannot"` -eq 1 ];then
            ls -1 *.fasta.gz 2>&1 | head -n 1 >.unilist.txt
            if [ `cat .unilist.txt | grep -c "ls:\ cannot"` -eq 1 ];then
                echo -e "\n\t\e[31m -- ERROR: Directory \"$myuni\" is empty. Please download a UNIPROT database and rerun the pre-check --\e[39m\n"
                rm .unilist.txt
                exit 0
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
    jav=$( java -version 2>&1 | awk '{print $3}' | grep [0-9] | cut -f 1,2 -d "." | tr -d "." | tr -d "\"" )
    if [ $jav -eq 18 ] || [ $jav -eq 110 ];then
        echo -e "\n\t -- Java 1.8 (or later) is installed -- \n"
	rep=yes
    elif [ $jav -eq 17 ];then
        echo -e "\n\t\e[31m -- ERROR: Please install Java 1.8 (or later). Requirement for Nextflow --\e[39m\n"
        exit 0
    fi
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
            echo -e -n "\n\t    Do you want me to try to install Nextflow for you? (y or n): "
            read ans
            case $ans in
                [yY] | [yY][eE][sS])
                    java_c
                    if [ "$rep" == "yes" ];then
                        echo -e "\n\t -- Downloading Nextflow ... -- \n"
                        curl -s https://get.nextflow.io | bash
		                echo -e "\n\t -- Nextflow is now installed on $mypwd (local installation) -- \n"
                    fi
                ;;
                [nN] | [nN][oO])
                    echo -e "\n\t\e[31m -- ERROR: Download and Install Nextflow. Then rerun the pre-check again --\e[39m\n"
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
evi_bash () {
    if [ `cat ~/.bashrc | grep -c "evigene"` -eq 0 ];then
        echo -e "\n\t -- Adding info and sourcing .bashrc file -- \n"
        echo -e "# EvidentialGene\nexport PATH=\"\$PATH:${mypwd}/evigene/scripts/prot/\"\n# EvidentialGene(other scripts)\nexport PATH=\"\$PATH:${mypwd}/evigene/scripts/\"\n \
        export PATH=\"\$PATH:${mypwd}/evigene/scripts/ests/\"\nexport PATH=\"\$PATH:${mypwd}/evigene/scripts/genes/\"\nexport PATH=\"\$PATH:${mypwd}/evigene/scripts/genoasm/\"\n \
        export PATH=\"\$PATH:${mypwd}/evigene/scripts/omcl/\"\nexport PATH=\"\$PATH:${mypwd}/evigene/scripts/rnaseq/\"\n" >> ~/.bashrc
        source ~/.bashrc
    else
        source ~/.bashrc
    fi
}
evi_c () {
    echo -e "\n\t -- EvidentialGene installation -- \n"
    check_evi=$( command -v tr2aacds.pl | wc -l )
    if [ $check_evi -eq 0 ];then
        echo -e "\n\t -- EvidentialGene is not installed -- \n"
        echo -e -n "\n\t    Do you want me to try to install EvidentialGene for you? (y or n): "
        read ans
        case $ans in
            [yY] | [yY][eE][sS])
                echo -e "\n\t -- Downloading EvidentialGene ... -- \n"
                wget http://arthropods.eugenes.org/EvidentialGene/other/evigene_old/evigene_older/evigene19jan01.tar
                tar -xf evigene19jan01.tar
                mv evigene19jan01/ evigene/
                rm evigene19jan01.tar
            ;;
            [nN] | [nN][oO])
                echo -e "\n\t\e[31m -- ERROR: Download and Install EvidentialGene. Then rerun the pre-check again --\e[39m\n"
                exit 0
            ;;
            *)
                echo -e "\n\n\t\e[31m -- Yes or No answer not specified. Try again --\e[39m\n"
                evi_c
            ;;
        esac
    elif [ $check_evi -eq 1 ];then
        echo -e "\n\t -- EvidentialGene is already installed -- \n"
    fi
}
get_var () {
    cd $mypwd
    #echo "=$mypwd/" >${mypwd}/.varfile.sh
    echo "buscodb=$mypwd/busco_db/$busna" >${mypwd}/.varfile.sh
    echo "uniname=$unina" >>${mypwd}/.varfile.sh
    echo "uniprot=$mypwd/uniprot_db/$unina" >>${mypwd}/.varfile.sh
    echo "pfloc=$mypwd/hmmerdb/Pfam-A.hmm" >>${mypwd}/.varfile.sh
    echo "pfname=Pfam-A.hmm" >>${mypwd}/.varfile.sh
    echo "nextflow=$mypwd/nextflow" >>${mypwd}/.varfile.sh
    vpwd=$mypwd
    echo "mypwd=$mypwd" >>${vpwd}/.varfile.sh
    source .varfile.sh
    echo -e "\n\t -- INFO to use in the pipeline --\n"
    echo -e "\t Pipeline PATH:\t\t $mypwd"
    echo -e "\t BUSCO database:\t $buscodb"
    echo -e "\t UNIPROT database:\t $uniprot"
    echo -e "\t PFAM files:\t\t $pfloc"
    echo -e "\t NEXTFLOW:\t\t $nextflow \n\n"
    cat sample.nextflow.config | sed -e "s|mypwd|mypwd=\"${mypwd}\"|" -e "s|buscodb|buscodb=\"${buscodb}\"|" -e "s|uniprot|uniprot=\"${uniprot}\"|" \
        -e "s|uniname|uniname=\"${uniname}\"|" -e "s|pfloc|pfloc=\"${pfloc}\"|" -e "s|pfname|pfname=\"${pfname}\"|" >nextflow.config
    evi_bash
}
#Main
if [ "$mypwd" == "" ] || [ "$mypwd" == "-h" ] || [ "$mypwd" == "-help" ] || [ "$mypwd" == "--help" ];then
    echo -e "\n\t Script for checking the requirenments of TransPi \n"
    echo -e "\t Usage:\n\n\t\t pre-check_TransPi.sh WORK_PATH \n"
    echo -e "\n\t\t\t WORK_PATH = PATH to run TransPi and download the requirenments \n\n\t\t\t\t Example: /home/bioinf/run/ \n"
    exit 0
elif [ ! -d "$mypwd" ];then
    echo -e "\n\t -- Please provide a valid PATH to run TransPi -- \n"
    exit 0
elif [ -d "$mypwd" ];then
    cd $mypwd
    read_c
    conda_c
    pfam_c
    bus_c
    uniprot_c
    nextflow_c
    evi_c
    echo -e "\n\t -- If no \"ERROR\" was found and all the neccesary databases are installed proceed to run TransPi -- \n"
    get_var
fi
