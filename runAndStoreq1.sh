#!/bin/bash
mkdir q1output
for j in {1..5}
do
    for i in {1..20}
    do
        python smith-waterman.py  -C matrices/BLOSUM50 Pospairs.txt q1output/Positives_Blosum50_${i}_${j}.txt ${i} ${j}
        python smith-waterman.py  -C matrices/BLOSUM50 Negpairs.txt q1output/Negatives_Blosum50_${i}_${j}.txt ${i} ${j}
    done
done
cd q1output
bash ../evaluate1.sh > evaluation.tsv
cd ..
