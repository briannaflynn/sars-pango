#!/usr/bin/bash

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
