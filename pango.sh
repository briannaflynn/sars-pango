#!/usr/bin/bash

source $2/etc/profile.d/conda.sh
conda activate pangolin
# enter the virus fasta as the first positional arg 
pangolin $1_consensus2.fasta
echo "Pangolin complete, results in lineage_report.csv"
