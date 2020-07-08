jfile="$1"
sampleid=$( echo $jfile | cut -f 1 -d "." )
tb=$( jq '.summary.before_filtering.total_reads' $jfile )
r1b=$( jq '.read1_before_filtering.total_reads' $jfile )
r1bl=$( jq '.summary.before_filtering.read1_mean_length' $jfile )
r2b=$( jq '.read2_before_filtering.total_reads' $jfile )
r2bl=$( jq '.summary.before_filtering.read2_mean_length' $jfile )
ta=$( jq '.summary.after_filtering.total_reads' $jfile )
r1a=$( jq '.read1_after_filtering.total_reads' $jfile )
r1al=$( jq '.summary.after_filtering.read1_mean_length' $jfile )
r2a=$( jq '.read2_after_filtering.total_reads' $jfile )
r2al=$( jq '.summary.after_filtering.read2_mean_length' $jfile )
loss=$( echo "${tb}-${ta}" | bc )
gcb=$( jq '.summary.before_filtering.gc_content' $jfile )
gca=$( jq '.summary.after_filtering.gc_content' $jfile )
echo "Sample_name,Total_before,R1_before,R1_before_length,R2_before,R2_before_length,GC_before,Total_after,R1_after,R1_after_length,R2_after,R2_after_length,GC_after,Reads_discarded" >${sampleid}_reads_stats.csv
echo "${sampleid},${tb},${r1b},${r1bl},${r2b},${r2bl},${gcb},${ta},${r1a},${r1al},${r2a},${r2al},${gca},${loss}" >>${sampleid}_reads_stats.csv
