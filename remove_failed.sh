#!/bin/bash
#Script to remove directories of FAILED and ABORTED processes in a nextflow pipeline
#INPUT = filename_trace.txt
file=$1
if [ "$file" == "" ];then
    echo -e "\n\t Provide a trace file as input (e.g. filename_trace.txt)"
    echo -e "\n\t Usage: bash remove_failed.sh filename_trace.txt\n"
    exit 0
else
    cat $file | grep "ABORTED" >.erase.txt
    cat $file | grep "FAILED" >>.erase.txt
    while read line;do
        a=$( echo $line | awk '{print $2}' )
        echo $a
        if [ -d work/${a}* ];then
            rm -rf work/${a}*
        fi
    done <.erase.txt
    rm .erase.txt
fi
