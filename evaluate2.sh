#!/bin/bash
for i in ../matrices/*
do
    name=`basename ${i}`
    min=`cat Positives_${name}.txt | sort -k 3,3 -n -r | head -35 | tail -1 | cut -f 3`
    echo -ne $i $j $min ' '
    cat Negatives_${name}.txt | sort -k 3,3 -n -r | awk -v min=${min} 'BEGIN {x=0}; {if($3<min){print $3,x; exit 1}; x++}'      
done
