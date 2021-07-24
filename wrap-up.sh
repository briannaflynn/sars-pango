#!/usr/bin/bash

cp pango-next/data/combined_sequences.fasta pango-next/results
cp -r pango-next/results ./nextstrain_results
rm -r pango-next/results

mkdir $1-results; mkdir $1-results/genetic_data
mv *.pdf $1-results
#mv *$1_* $1-results/genetic_data/
#mv *$1.* $1-results/genetic_data/
mv $1_consensus2.fasta $1-results/genetic_data/
mv $1_calls.vcf.gz.tmp.csi $1-results/genetic_data/
mv $1_calls.vcf.gz.tmp $1-results/genetic_data/
mv $1.vcf.tmp $1-results/genetic_data/
mv $1_sorted.bam.bai $1-results/genetic_data/
mv $1_sorted.bam $1-results/genetic_data/
mv $1.bam.tmp $1-results/genetic_data/
mv $1.sam.tmp $1-results/genetic_data/
mv $1_read2.trim.fastq.tmp $1-results/genetic_data/
mv $1_read1.trim.fastq.tmp $1-results/genetic_data/
mv $1_interleaved.tmp $1-results/genetic_data/
mv $1_combined.fasta $1-results/genetic_data/

mv nextstrain_results $1-results/nextstrain_results
mv lineage_report.csv $1-results
mv clade_assignment.tsv $1-results

mv $1_variants.vcf $1-results/genetic_data/
mv _combined.fasta $1-results/genetic_data/

mkdir $1-results/tree
cp $1-results/nextstrain_results/global/tree_raw.nwk $1-results/tree
cp $1-results/nextstrain_results/global/clades.json $1-results/tree



