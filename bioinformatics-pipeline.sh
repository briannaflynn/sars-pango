#!/bin/bash

source $4/etc/profile.d/conda.sh
conda activate bioinformatics

# trim raw fastq data
cutadapt -j 12 -g GTTTCCCAGTCACGATA -G GTTTCCCAGTCACGATA -a TATCGTGACTGGGAAAC -A TATCGTGACTGGGAAAC -g ACACTCTTTCCCTACACGACGCTCTTCCGATCT -G ACACTCTTTCCCTACACGACGCTCTTCCGATCT -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -n 3 -m 75 -q 25 --interleaved *"$2" *"$3" > "$1_interleaved.tmp"
cutadapt -j 12 -m 75 -u 25 -u -25 -U 25 -U -25 --interleaved "$1_interleaved.tmp" -o "$1_read1.trim.fastq.tmp" -p "$1_read2.trim.fastq.tmp"

# map to reference, get vcf
bowtie2 -p 12 -x "./index/Sars2" -1 "$1_read1.trim.fastq.tmp" -2 "$1_read2.trim.fastq.tmp" -S "$1.sam.tmp"
echo "Completed bowtie2 alignment"
samtools view -S "$1.sam.tmp" -b > "$1.bam.tmp"
echo "Completed samtools view"
samtools sort --threads 12 "$1.bam.tmp" -o "$1_sorted.bam"
echo "Completed samtools sort"
samtools index "$1_sorted.bam"
echo "Samtools index"
bcftools mpileup --threads 12 -f "./index/Sars2.fna" "$1_sorted.bam" > "$1.vcf.tmp"
echo "bcftools to vcf"

# vcf to fasta
bcftools call -Mc --ploidy 1 --threads 12 -Oz "$1.vcf.tmp" -o "$1_calls.vcf.gz.tmp"
echo "bcftools call completed"
bcftools index "$1_calls.vcf.gz.tmp" 
echo "bcftools index completed"
cat "./index/Sars2.fna" | bcftools consensus "$1_calls.vcf.gz.tmp" > "$1_consensus2.fasta"
echo "vcf to fasta completed"

# get calls only
bcftools call -Mc --ploidy 1 -v --threads 12 -Ov $1.vcf.tmp -o $1_variants.vcf 
