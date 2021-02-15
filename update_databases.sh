#!/usr/bin env bash
mypwd=$( pwd )
db_c () {
    if [ ! -d DBs ];then
        mkdir DBs
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
    echo -e "\n\t -- Custom sqlite database for Trinotate will be installed -- \n"
    echo -e "\n\t -- This could take a couple of minutes depending on connection. Please wait -- \n"
    wget https://github.com/Trinotate/Trinotate/archive/Trinotate-v3.2.1.tar.gz
    tar -xf Trinotate-v3.2.1.tar.gz
    mv Trinotate-Trinotate-v3.2.1/ Trinotate_build_scripts/
    ./Trinotate_build_scripts/admin/Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate
    rm uniprot_sprot.dat.gz
    date -u >.lastrun.txt
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
                exit 0
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
        echo -e "\n\t -- SQLite direcotry found at ${mypwd}/DBs -- \n"
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
    echo -e "\n###################################################################\n"
}
message
db_c
downd
