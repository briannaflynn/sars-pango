#!/usr/bin/bash

source /Users/brie/anaconda3/etc/profile.d/conda.sh
conda activate bioinformatics

python tree-maker.py $1-results/tree/tree_raw.nwk $1-results/tree/clades.json

mv tree-plot.pdf $1-results/tree/