#!/usr/bin/bash

bcftools call -Mc --ploidy 1 --threads 12 -Oz "$1.vcf.tmp" -o "$1_calls.vcf.gz.tmp"
echo "bcftools call completed"
bcftools index "$1_calls.vcf.gz.tmp" 
echo "bcftools index completed"
cat "./index/Sars2.fna" | bcftools consensus "$1_calls.vcf.gz.tmp" > "$1_consensus2.fasta"
echo "vcf to fasta completed"
