#!/bin/bash
for j in {1..5}
do
    for i in {1..20}
    do
        min=`cat Positives_Blosum50_${i}_${j}.txt | sort -k 3,3 -n -r | head -35 | tail -1 | cut -f 3`
        echo -ne $i $j $min ' '
        cat Negatives_Blosum50_${i}_${j}.txt | sort -k 3,3 -n -r | awk -v min=${min} 'BEGIN {x=0}; {if($3<min){print $3,x; exit 1}; x++}'
        
    done
done