#!/bin/bash
mkdir q3output
python -C smith-waterman.py matrices/BLOSUM50 Pospairs.txt q3output/Positives.txt 7 3
python -C smith-waterman.py matrices/BLOSUM50 Negpairs.txt q3output/Negatives.txt 7 3

