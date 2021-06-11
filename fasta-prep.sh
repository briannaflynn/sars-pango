#!/usr/bin/bash

# replace name with SAMPLE
var=">SAMPLE"
sed -i "1s/.*/$var/" $1_consensus.fasta

# prepend it to the fasta with other sars-cov-2 viral genomes from genbank
cat $1_consensus.fasta ./pango-next/data/pruned_gb_sequences.fasta > $1_combined.fasta 

