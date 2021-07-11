#!/usr/bin/bash

source $2/etc/profile.d/conda.sh
conda activate nextstrain

echo "Prepping fasta sequence"
bash fasta-prep.sh $1

echo "Creating fasta sequence dataset for Nextstrain"
# prepend it to the fasta with other sars-cov-2 viral genomes from genbank
cat $1_consensus2.fasta ./pango-next/data/GISAID_3000.fasta ./pango-next/data/pruned_gb_sequences.fasta > $1_combined.fasta
cp $1_combined.fasta combined_sequences.fasta
mv combined_sequences.fasta pango-next/data/

echo "combined_sequences.fasta is now available in pango-next/data"

cd pango-next; pwd

snakemake --cores 2 --profile my_profiles/example -p
