#!/usr/bin/bash

source $2/etc/profile.d/conda.sh
conda activate bioinformatics

python tree-maker.py $1-results/tree/tree_raw.nwk $1-results/tree/clades.json
echo "If successful tree pdf now available: $1-results/tree/tree-plot.png"
mv tree-plot.png $1-results/tree/
#echo "If successful tree pdf now available: $1-results/tree/tree-plot.pdf"
#mv tree-plot.pdf $1-results/tree/
