#!/usr/bin/bash

cp pango-next/data/combined_sequences.fasta pango-next/results
cp -r pango-next/results ./nextstrain_results
rm -r pango-next/results

mkdir $1-results; mkdir $1-results/genetic_data
mv *.pdf $1-results
mv *$1_* $1-results/genetic_data
mv *$1.* $1-results/genetic_data

mv nextstrain_results $1-results/nextstrain_results
mv lineage_report.csv $1-results
mv clade_assignment.tsv $1-results

mkdir $1-results/tree
cp $1-results/nextstrain_results/global/tree_raw.nwk $1-results/tree
cp $1-results/nextstrain_results/global/clades.json $1-results/tree

#rm $1_consensus2.fasta
rm $1_calls.vcf.gz.tmp.csi
rm $1_calls.vcf.gz.tmp
rm $1.vcf.tmp
rm $1_sorted.bam.bai
rm $1_sorted.bam
rm $1.bam.tmp
rm $1.sam.tmp
rm $1_read2.trim.fastq.tmp
rm $1_read1.trim.fastq.tmp
rm $1_interleaved.tmp

