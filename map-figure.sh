#!/usr/bin/bash

source $2/etc/profile.d/conda.sh
conda activate bioinformatics

python maps.py $1
rm $1_map.pdf
mv $1_result.pdf $1-results/$1_worldmap.pdf
