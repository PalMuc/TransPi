#!/usr/bin env bash
export mypwd="$1"
pfam_c() {
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
    pfam_c
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
    cd $installDir
    if [ ! -d DBs/sqlite_db/ ];then
        echo -e "\n\t -- SQLite directory not found at ${installDir}/DBs -- \n"
        echo -e -n "\n\t -- Do you want to create the directory and istall the last version of the databses? (y or n): "
        read ans
        case $ans in
            [yY] | [yY][eE][sS])
                mkdir -p DBs/sqlite_db/
                cd DBs/sqlite_db/
                sqld
            ;;
            [nN] | [nN][oO])
                echo -e "\n\n\t\e[31m -- Exiting program --\e[39m\n"
                exit 0
            ;;
            *)
                echo -e "\n\n\t\e[31m -- Yes or No answer not specified. Try again --\e[39m\n"
                downd
            ;;
        esac
    elif [ -d DBs/sqlite_db/ ];then
        echo -e "\n\t -- SQLite direcotry found at ${installDir}/DBs -- \n"
        cd DBs/sqlite_db/
        if [ ! -e *.sqlite ];then
            echo -e "\n\t -- Custom sqlite database for Trinotate is not installed -- \n"
            echo -e -n "\n\t    Do you want to install the custom sqlite database? (y or n): "
            read ans
            case $ans in
                [yY] | [yY][eE][sS])
                    echo -e "\n\t -- This could take a couple of minutes depending on connection. Please wait -- \n"
                    sqld
                ;;
                [nN] | [nN][oO])
                    echo -e "\n\t\e[31m -- Custom trinotate sqlite database not created --\e[39m\n"
                    exit 0
                ;;
                *)
                    echo -e "\n\n\t\e[31m -- Yes or No answer not specified. Try again --\e[39m\n"
                    downd
                ;;
            esac
        elif [ -e *.sqlite ];then
            echo -e "\n\t -- Custom sqlite database for Trinotate is installed -- \n"
            echo -e "\n\t -- Verifying when scripts was last run -- \n"
            ddate
        fi
    fi
}
message(){
    echo -e "\n###################################################################\n"
    echo -e "\n  Script for updating the databases used by TransPi \n"
    echo -e "\n  - SwissProt, PFAM, eggNOG, GO, and Trinotate SQL database - \n"
    echo -e "\n  You need the TransPi conda environment installed and the PATH to install/update the DBs \n"
    echo -e "\n  Example: /home/ubuntu/TransPi/ \n"
    echo -e "\n###################################################################\n"
}
if [ "$mypwd" == "" ] || [ "$mypwd" == "-h" ] || [ "$mypwd" == "-help" ] || [ "$mypwd" == "--help" ];then
    message
elif [ ! -d "$mypwd" ];then
    echo -e "\n\t -- Please provide a valid PATH to install/update the DBs -- \n"
    exit 0
elif [ -d "$mypwd" ];then
    if [ ${mypwd} == "." ];then
        mypwd=$(pwd)
        installDir=$(pwd)
    elif [ ${mypwd} == $(pwd) ]; then
        installDir=$(pwd)
    else
        installDir=${mypwd}
    fi
    downd
fi
