filename="$1"
cat $filename | cut -f 1 -d "," | awk '{print $3}' >nam
cat $filename | cut -f 2 -d "," | awk '{print $2}' >len
paste nam len | awk '$2<2500 {print $0}' >lt_2500
he=$( paste nam len | awk '$2>=2500 {print $0}' | wc -l )
for x in `seq 1 $he`;do
    echo data >>temp_1
    echo $x >>temp_2
done
paste temp_1 temp_2 >he_2500
cat lt_2500 he_2500 >final_sizes.txt
rm nam len temp_1 temp_2 lt_2500 he_2500 $filename
