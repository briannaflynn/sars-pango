#!/usr/bin/bash

bowtie2 -p 12 -x "./index/Sars2" -1 "O25799359_read1.trim.fastq.tmp" -2 "O25799359_read2.trim.fastq.tmp" -S "O25799359.sam.tmp"
echo "Completed bowtie2 alignment"
samtools view -S "O25799359.sam.tmp" -b > "O25799359.bam.tmp"
echo "Completed samtools view"
samtools sort --threads 12 "O25799359.bam.tmp" -o "O25799359_sorted.bam"
echo "Completed samtools sort"
samtools index "O25799359_sorted.bam"
echo "Samtools index"
bcftools mpileup --threads 12 -f "./index/Sars2.fna" "O25799359_sorted.bam" > "O25799359.vcf.tmp"
echo "bcftools to vcf"
