jfile="$1"
sampleid=$( echo $jfile | cut -f 1 -d "." )
r1bn=$( jq '.read1_before_filtering.quality_curves.mean' ${jfile} | grep -c [0-9] )
r1bq=$( jq '.read1_before_filtering.quality_curves.mean' ${jfile} | grep [0-9] | tr -d "\n" | tr -d " " )
r2bn=$( jq '.read2_before_filtering.quality_curves.mean' ${jfile} | grep -c [0-9] )
r2bq=$( jq '.read2_before_filtering.quality_curves.mean' ${jfile} | grep [0-9] | tr -d "\n" | tr -d " " )
r1an=$( jq '.read1_after_filtering.quality_curves.mean' ${jfile} | grep -c [0-9] )
r1aq=$( jq '.read1_after_filtering.quality_curves.mean' ${jfile} | grep [0-9] | tr -d "\n" | tr -d " " )
r2an=$( jq '.read2_after_filtering.quality_curves.mean' ${jfile} | grep -c [0-9] )
r2aq=$( jq '.read2_after_filtering.quality_curves.mean' ${jfile} | grep [0-9] | tr -d "\n" | tr -d " " )
echo -e "${r1bn}\n${r1bq}\n${r2bn}\n${r2bq}\n${r1an}\n${r1aq}\n${r2an}\n${r2aq}" >${sampleid}_reads_qual.csv
