#!/usr/bin/env bash

export mypwd="$1"
export busdb="$2"
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
            exit 0
        else
            ls -1 reads/*.gz >.readlist.txt
            if [ $(( `cat .readlist.txt | wc -l` % 2 )) -eq 0 ];then
                echo -e "\n\t -- Reads found in $( pwd )/reads/ and in pairs -- \n"
                nind=$( cat .readlist.txt | cut -f 2 -d "/" | sed 's/R\{1,2\}.*//g' | sort -u | wc -l )
                echo -e "\n\t -- Number of samples: $nind -- \n"
            elif [ $(( `cat .readlist.txt | wc -l` % 2 )) -eq 1 ];then
                echo -e "\n\t\e[31m -- ERROR: Reads found in $( pwd )/reads/ but not in pairs. Make sure you have an R1 and R2 for each sample --\e[39m\n"
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
    if [ "$check_conda" == "conda" ] && [ `conda -V | cut -f 2 -d " "` > 4.5 ];then
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
        echo -e -n "\n\t    Do you want me to try to install Anaconda for you? (y or n): "
        read ans
        if [ "$ans" == "y" ] || [ "$ans" == "yes" ] || [ "$ans" == "Y" ] || [ "$ans" == "YES" ] || [ "$ans" == "Yes" ];then
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
        elif [ "$ans" == "n" ] || [ "$ans" == "no" ] || [ "$ans" == "N" ] || [ "$ans" == "NO" ] || [ "$ans" == "No" ];then
            echo -e "\n\t\e[31m -- ERROR: Download and Install Anaconda. Then rerun the pre-check again --\e[39m\n"
            exit 0
        else
            echo -e "\n\t\e[31m -- ERROR: Yes or No answer not specified. Rerun the pre-check again --\e[39m\n"
            exit 0
        fi
    fi
}
pfam_c() {
    #Check PFAM files
    cd $mypwd
    if [ ! -d hmmerdb/ ];then
        echo -e "\n\t -- Creating directory for the HMMER database --\n"
        mkdir hmmerdb
        cd hmmerdb
        echo -e "\t -- Downloading current release of PFAM for the HMMER database --\n"
        wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
        gunzip Pfam-A.hmm.gz
    elif [ -d hmmerdb/ ];then
        cd hmmerdb
        if [ -f Pfam-A.hmm ];then
            echo -e "\n\t -- Pfam file is present --\n"
        else
            echo -e "-- Downloading current release of PFAM for the HMMER database --\n"
            wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
            echo -e "-- Preparing files ... --\n"
            gunzip Pfam-A.hmm.gz
        fi
    fi
}
bus_down () {
    name=$1
    cd $mypwd
    cd busco_db
    echo -e "\n\t -- Downloading BUSCO \"$name\" database --\n";wait
    wname=$( echo "https://busco.ezlab.org/datasets/${name}.tar.gz" )
    wget $wname
    echo -e "\n\t -- Preparing files ... --\n";wait
    tar -xvf ${name}.tar.gz
    echo -e "\n\t -- DONE with BUSCO database --\n";wait
    name2=$( echo $name | cut -f 1 -d "_" )
    if [ -d ${name2}_odb9 ];then
        export busna=${name2}_odb9
    elif [ -d ${name2}_odb10 ];then
        export busna=${name2}_odb10
    fi
}
bus_c () {
    PS3="Please select one (1-5): "                             
    select var in `cat buslist.txt | grep "##" | tr -d "#"`;do
    case $var in
        BACTERIA)
            echo -e "\n\t You selected BACTERIA. Which specific database? \n"
            PS3="Please select database: "
            select var1 in `cat buslist.txt | sed -n "/##BACTERIA/,/#/p" | grep -v "##" | tr -d "#" | cut -f 5- -d "/" | cut -f 1 -d "."`;do
            if [ $var1 == "MAIN_MENU" ];then
                bus_c
            else
                bus_down $var1
                exit 0
            fi
            done
        ;;
        EUKARYOTA)
            echo -e "\n\tYou selected EUKARYOTA. Which specific database? \n"
            PS3="Please select database: "
            select var1 in `cat buslist.txt | sed -n "/##EUKARYOTA/,/#/p" | grep -v "##" | tr -d "#" | cut -f 5- -d "/" | cut -f 1 -d "."`;do
            if [ $var1 == "MAIN_MENU" ];then
                bus_c
            else
                bus_down $var1
                exit 0
            fi
            done
        ;;
        FUNGI)
            echo -e "\n\tYou selected FUNGI. Which specific database? \n"
            PS3="Please select database: "
            select var1 in `cat buslist.txt | sed -n "/##FUNGI/,/#/p" | grep -v "##" | tr -d "#" | cut -f 5- -d "/" | cut -f 1 -d "."`;do
            if [ $var1 == "MAIN_MENU" ];then
                bus_c
            else
                bus_down $var1
                exit 0
            fi
            done
        ;;
        PLANTS)
            echo -e "\n\tYou selected PLANTS. Which specific database? \n"
            PS3="Please select database: "
            select var1 in `cat buslist.txt | sed -n "/##PLANTS/,/#/p" | grep -v "##" | tr -d "#" | cut -f 5- -d "/" | cut -f 1 -d "."`;do
            if [ $var1 == "MAIN_MENU" ];then
                bus_c
            else
                bus_down $var1
                exit 0
            fi
            done
        ;;
        EXIT)
            echo -e "\n\t Exiting \n"
            exit 0
        ;;
        *)
            echo -e "\n\t Wrong option. Start again \n"
            exit 1
            ;;
    esac
    done
}
buscodb () {
    #Check BUSCO database
    name=$1
    cd $mypwd
    if [ -f buslist.txt ];then
        wname=$( cat buslist.txt | grep "/$name" )
        if [ `echo "$wname" | wc -l` -eq 1 ];then
            cd busco_db
            echo -e "-- Downloading BUSCO \"$name\" database --\n"
            wget $wname
            echo -e "-- Preparing files ... --\n"
            tar -xvf $wname
            name2=$( echo $name | cut -f 2 -d "/" )
            if [ -d ${name}_odb9 ];then
                export busna=${name}_odb9
            elif [ -d ${name2}_odb10 ];then
                export busna=${name2}_odb10
            fi
        else
            echo -e "\n\t\e[31m -- ERROR: Did you provide a valid name for the BUSCO database? -- \n\t -- For more info: \"pre-check_TransPi.sh -busco\" --\e[39m\n\n"
            exit 0
        fi
    else
        echo -e "\n\t\e[31m -- ERROR: Please make sure that file \"buslist.txt\" is available. Please check requirements and rerun the pre-check --\e[39m\n\n"
        exit 0
    fi
}
buscodb2 () {
    #Check BUSCO database2
    name=$1
    cd $mypwd
    if [ `echo $name | cut -f 1 -d "/"` != "prerelease" ];then
        cd busco_db
        if [ -d ${name}_odb9 ];then
            echo -e "\n\t -- BUSCO \"$name\" database found -- \n"
            export busna=${name}_odb9
        else
            echo -e "\n\t -- BUSCO \"$name\" database not found -- \n"
            rep="yes"
        fi
    elif [ `echo $name | cut -f 1 -d "/"` == "prerelease" ];then
        cd busco_db
        name2=$( echo $name | cut -f 2 -d "/" )
        if [ -d ${name2}_odb10 ];then
            echo -e "\n\t -- BUSCO \"$name\" database found -- \n"
            export busna=${name2}_odb10
        else
            echo -e "\n\t -- BUSCO \"$name\" database not found -- \n"
            rep="yes"
        fi
    fi
}
busco_c () {
    #BUSCO database
    cd $mypwd
    if [ ! -d busco_db/ ];then
        echo -e "\n\t -- Creating directory for the BUSCO database --\n"
        mkdir busco_db
        cd busco_db
        buscodb $busdb
    elif [ -d busco_db/ ];then
        cd busco_db
        buscodb2 $busdb
        if [ $rep == "yes" ];then
            buscodb $busdb
        fi
    fi
}
unicomp_c () {
    echo -e -n "\n\t    Do you want me to uncompress the file(s)? (y or n): \n"
    read ans
    if [ "$ans" == "y" ] || [ "$ans" == "yes" ] || [ "$ans" == "Y" ] || [ "$ans" == "YES" ] || [ "$ans" == "Yes" ];then
        echo -e "\n\t -- Uncompressing file(s) ... -- \n"
        gunzip *fasta.gz
        uniprot_c
    elif [ "$ans" == "n" ] || [ "$ans" == "no" ] || [ "$ans" == "N" ] || [ "$ans" == "NO" ] || [ "$ans" == "No" ];then
        echo -e "\n\t\e[31m -- ERROR: Please uncompress the file(s) and rerun the pre-check again --\e[39m\n"
        exit 0
    else
        echo -e "\n\t\e[31m -- ERROR: Yes or No answer not specified. Rerun the pre-check again --\e[39m\n"
        exit 0
    fi
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
        echo -e "\n\t -- Example: \"uniprot-taxonomy_metazoaA33208.fasta\" (file needs to be uncompressed) -- \n"
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
                echo -e "\n\t\e[31m -- ERROR: Directory \"$myuni\" is empty. Please download a UNIPROT database and rerun the pre-check files --\e[39m\n"
                rm .unilist.txt
                exit 0
            else
                echo -e "\n\t\e[31m -- Directory \"$myuni\" is available but UNIPROT database is compressed --\e[39m\n"
                unicomp_c
            fi
        else
            echo -e "\n\t -- Here is the list of UNIPROT files found at: $myuni -- \n"
            ls *fasta
            echo -e -n "\n\t    Which UNIPROT database you want to use? (paste answer here): "
            read ans
            echo -e -n "\n\t -- UNIPROT database selected: \"$ans\" --\n"
            export unina=${ans}
            rm .unilist.txt
        fi
    fi
}
java_c () {
    jav=$( java -version 2>&1 | awk '{print $3}' | grep "1.8" | wc -l )
    if [ $jav -eq 1 ];then
        echo -e "\n\t -- Java 1.8 (or later) is installed -- \n"
        rep=yes
    elif [ $jav -eq 0 ];then
        echo -e "\n\t\e[31m -- ERROR: Please install Java 1.8 (or later). Requirement of Nextflow --\e[39m\n"
        exit 0
    fi
}
nextflow_c () {
    #Check Nextflow
    cd $mypwd
    check_next=$( command -v nextflow )
    if [ `echo $check_next | wc -l` -eq 1 ];then
        echo -e "\n\t -- Nextflow is intalled -- \n"
    elif [ `echo $check_next | wc -l` -eq 1 ];then
        echo -e "\n\t -- Nextflow is not intalled -- \n"
        echo -e -n "\n\t    Do you want me to try to install Nextflow for you? (y or n): "
        read ans
        if [ "$ans" == "y" ] || [ "$ans" == "yes" ] || [ "$ans" == "Y" ] || [ "$ans" == "YES" ] || [ "$ans" == "Yes" ];then
            java_c
            if [ $rep == "yes" ];then
                echo -e "\n\t -- Downloading Nextflow ... -- \n"
                curl -s https://get.nextflow.io | bash
                echo -e "\n\t -- Nextflow is now installed on $mypwd -- \n"
            fi
        elif [ "$ans" == "n" ] || [ "$ans" == "no" ] || [ "$ans" == "N" ] || [ "$ans" == "NO" ] || [ "$ans" == "No" ];then
            echo -e "\n\t\e[31m -- ERROR: Download and Install Nextflow. Then rerun the pre-check again --\e[39m\n"
            exit 0
        else
            echo -e "\n\t\e[31m -- ERROR: Yes or No answer not specified. Rerun the pre-check again --\e[39m\n"
            exit 0
        fi
    fi
}
get_var () {
    cd $mypwd
    #echo "=$mypwd/" >${mypwd}/.varfile.sh
    busname=$busn
    echo "buscodb=$mypwd/busco_db/$busname" >${mypwd}/.varfile.sh
    echo "uniname=$unina" >>${mypwd}/.varfile.sh
    echo "uniprot=$mypwd/uniprot_db/$uniname" >>${mypwd}/.varfile.sh
    echo "pfloc=$mypwd/hmmerdb/Pfam-A.hmm" >>${mypwd}/.varfile.sh
    echo "pfname=Pfam-A.hmm" >>${mypwd}/.varfile.sh
    echo "nextflow=$mypwd/nextflow" >>${mypwd}/.varfile.sh
    vpwd=$mypwd
    echo "mypwd=$mypwd" >>${vpwd}/.varfile.sh
}
#Main
if [ "$mypwd" == "" ] || [ "$mypwd" == "-h" ] || [ "$mypwd" == "-help" ] || [ "$mypwd" == "--help" ];then
    echo -e "\n\t Script for checking the requirenments of TransPi \n"
    echo -e "\t Usage:\n\n\t\t pre-check_TransPi.sh WORK_PATH BUSCO_DB \n"
    echo -e "\n\t\t\t WORK_PATH = PATH to run TransPi and download the requirenments \n\n\t\t\t\t Example: /home/bioinf/run/ \n"
    echo -e "\n\t\t\t BUSCO_DB = BUSCO database to run the analysis \n\n\t\t\t\t Example: metazoa,bacteria,fungi,etc. \n\n"
    echo -e "\n\t\t\t To view the entire list of BUSCO database available run \"pre-check_TransPi.sh -busco\" \n\n"
    exit 0
elif [ "$mypwd" == "-busco" ] || [ "$mypwd" == "--busco" ];then
    echo -e "\n\t Script for checking the requirenments of TransPi \n"
    echo -e "\t Usage:\n\n\t\t pre-check_TransPi.sh WORK_PATH BUSCO_DB \n"
    echo -e "\n\t\t\t WORK_PATH = PATH to run TransPi and download the requirenments \n\n\t\t\t\t Example: /home/bioinf/run/ \n"
    echo -e "\n\t\t\t BUSCO_DB = BUSCO database to run the analysis \n\n\t\t\t\t Example: metazoa, bacteria, fungi, prerelease/solanaceae, etc. \n"
    echo -e "\n\t
#BACTERIA                       #EUKARYOTA              #FUNGI                      #PLANTS prerelease = odb10
    bacteria                        eukaryota               fungi                       embryophyta
    proteobacteria                  protists                microsporidia               prerelease/viridiplantae
    rhizobiales                     alveolata               dikarya                     prerelease/chlorophyta
    betaproteobacteria              metazoa                 ascomycota                  prerelease/embryophyta
    gammaproteobacteria             nematoda                pezizomycotina              prerelease/liliopsida
    enterobacteriales               arthropoda              eurotiomycetes              prerelease/eudicotyledons
    deltaepsilonsub                 insecta                 sordariomyceta              prerelease/solanaceae
    actinobacteria                  endopterygota           saccharomyceta
    cyanobacteria                   hymenoptera             saccharomycetales
    firmicutes                      diptera                 basidiomycota
    clostridia                      vertebrata
    lactobacillales                 actinopterygii
    bacillales                      tetrapoda
    bacteroidetes                   aves
    spirochaetes                    mammalia
    tenericutes                     euarchontoglires
    laurasiatheria
    \n"
    exit 0
elif [ ! -d "$mypwd" ] && [ "$busdb" != "" ];then
    echo -e "\n\t -- Please provide a valid PATH to run TransPi -- \n"
    exit 0
elif [ ! -d "$mypwd" ] && [ "$busdb" == "" ];then
    echo -e "\n\t -- Please provide a valid PATH to run TransPi and a BUSCO database name -- \n"
    exit 0
elif [ -d "$mypwd" ] && [ "$busdb" == "" ];then
    echo -e "\n\t -- Please provide the name of the BUSCO dataset to use (e.g. metazoa) -- \n"
    echo -e "\n\t -- To view the entire list of BUSCO database available run \"pre-check_TransPi.sh -busco\" -- \n"
    exit 0
elif [ -d "$mypwd" ] && [ "$busdb" != "" ];then
    cd $mypwd
    read_c
    conda_c
    pfam_c
    busco_c
    uniprot_c
    nextflow_c
    echo -e "\n\t -- If no \"ERROR\" was found and all the neccesary databases are installed \(incluiding UNIPROT\) proceed to run TransPi -- \n"
    get_var
    source .varfile.sh
fi
