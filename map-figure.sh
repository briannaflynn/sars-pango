#!/usr/bin/bash

source $2/etc/profile.d/conda.sh
conda activate bioinformatics

python maps.py $1

mv map.pdf $1-results/
