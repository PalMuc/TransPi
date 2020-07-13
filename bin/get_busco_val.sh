name_tri=$1
name_transpi=$2
version=$3
a=$4
if [ "$version" == "v3" ];then
    #trinity
    for x in $name_tri;do
        echo "'${a}','${a}','${a}','${a}'," >>tspec.txt
        b=$( cat $x | grep "(C)" -A5 | awk '{print $1}' | awk -v RS= -v OFS=, '{$1 = $1} 1' | cut -f 2,3,4,5 -d "," )
        echo "${b}," >>tnum.txt
        c=$( cat $x | grep "C:" | cut -f 2 -d "[" | cut -f 1,2,3,4 -d "," | tr -d "%" | tr -d "]" | tr -d [A-Z] | tr -d ":" )
        echo "${c}," >>tperc.txt
    done
    #transpi
    for x in $name_transpi;do
        echo "'${a}_TP','${a}_TP','${a}_TP','${a}_TP'" >>pspec.txt
        b=$( cat $x | grep "(C)" -A5 | awk '{print $1}' | awk -v RS= -v OFS=, '{$1 = $1} 1' | cut -f 2,3,4,5 -d "," )
        echo "${b}" >>pnum.txt
        c=$( cat $x | grep "C:" | cut -f 2 -d "[" | cut -f 1,2,3,4 -d "," | tr -d "%" | tr -d "]" | tr -d [A-Z] | tr -d ":" )
        echo "${c}" >>pperc.txt
    done
    paste tspec.txt pspec.txt | tr "\t" "\n" | tr -d "\n" >final_spec
    paste tnum.txt pnum.txt | tr "\t" "\n" | tr -d "\n" >final_num
    paste tperc.txt pperc.txt | tr "\t" "\n" | tr -d "\n" >final_perc
    rm tnum.txt tperc.txt tspec.txt
    rm pnum.txt pperc.txt pspec.txt
elif [ "$version" == "v4" ];then
    #trinity
    for x in $name_tri;do
        echo "'${a}','${a}','${a}','${a}'," >>tspec.txt
        b=$( cat $x | grep "(C)" -A5 | awk '{print $1}' | awk -v RS= -v OFS=, '{$1 = $1} 1' | cut -f 2,3,4,5 -d "," )
        echo "${b}," >>tnum.txt
        c=$( cat $x | grep "C:" | cut -f 2 -d "[" | cut -f 1,2,3,4 -d "," | tr -d "%" | tr -d "]" | tr -d [A-Z] | tr -d ":" )
        echo "${c}," >>tperc.txt
    done
    #transpi
    for x in $name_transpi;do
        echo "'${a}_TP','${a}_TP','${a}_TP','${a}_TP'" >>pspec.txt
        b=$( cat $x | grep "(C)" -A5 | awk '{print $1}' | awk -v RS= -v OFS=, '{$1 = $1} 1' | cut -f 2,3,4,5 -d "," )
        echo "${b}" >>pnum.txt
        c=$( cat $x | grep "C:" | cut -f 2 -d "[" | cut -f 1,2,3,4 -d "," | tr -d "%" | tr -d "]" | tr -d [A-Z] | tr -d ":" )
        echo "${c}" >>pperc.txt
    done
    paste tspec.txt pspec.txt | tr "\t" "\n" | tr -d "\n" >final_spec
    paste tnum.txt pnum.txt | tr "\t" "\n" | tr -d "\n" >final_num
    paste tperc.txt pperc.txt | tr "\t" "\n" | tr -d "\n" >final_perc
    rm tnum.txt tperc.txt tspec.txt
    rm pnum.txt pperc.txt pspec.txt
fi
