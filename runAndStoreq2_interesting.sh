#!/bin/bash
mkdir q2output_interesting
for i in matrices/*
do
    name=`basename ${i}`
    python smith-waterman.py -C ${i} Pospairs.txt q2output_interesting/Positives_${name}.txt 8 3
    python smith-waterman.py -C ${i} Negpairs.txt q2output_interesting/Negatives_${name}.txt 8 3
done
cd q2output_interesting
bash ../evaluate2.sh > evaluation.tsv
for sm in ../matrices/*
do
    name=`basename ${sm}`
    cat Positives_${name}.txt | cut -f 3 | sort -n -r > pos
    cat Negatives_${name}.txt | cut -f 3 | sort -n -r > neg
    ../ucsf-roc/a.out pos neg ${name}
done
cd ..
