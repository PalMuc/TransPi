#!/usr/bin env bash 
mypwd=$( pwd )
sqld(){
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
            rm Trinotate*
            Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate
            rm Pfam-A.hmm.gz uniprot_sprot.dat.gz
            date -u >.lastrun.txt 
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
            rm Trinotate*
            Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate
            rm Pfam-A.hmm.gz uniprot_sprot.dat.gz
            date -u >.lastrun.txt 
        fi 
    fi     
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
    if [ ! -d sqlite_db/ ];then
        echo -e "\n\t -- SQLite direcotry not found at $mypwd -- \n"
        echo -e -n "\n\t -- Do you want to create the directory and istall the last version of the databses? (y or n): "
        read ans
        case $ans in
            [yY] | [yY][eE][sS])
                mkdir sqlite_db/
                cd sqlite_db/
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
    elif [ -d sqlite_db/ ];then
        echo -e "\n\t -- SQLite direcotry found at $mypwd -- \n"
        cd sqlite_db/ 
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
    echo -e "\n  Script for checking when the Trninotate SQL databse was created \n"
    echo -e "\n  Make sure the conda env TransPi is installed and working \n"
    echo -e "\n###################################################################\n"
}
message
downd
