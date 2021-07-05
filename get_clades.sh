#!/bin/bash

source /Users/brie/anaconda3/etc/profile.d/conda.sh
conda activate nextstrain

cd pango-next; pwd

python scripts/assign_clades.py --sequences ../$1_consensus2.fasta

mv clade_assignment.tsv ../

rm *tmp*

echo "Clade assignment complete"
