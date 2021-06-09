#!/usr/bin/bash

bcftools call -Mc --ploidy 1 --threads 12 -Oz "O25799359.vcf.tmp" -o "O25799359_calls.vcf.gz.tmp"
echo "bcftools call completed"
bcftools index "O25799359_calls.vcf.gz.tmp" 
echo "bcftools index completed"
cat "./index/Sars2.fna" | bcftools consensus "O25799359_calls.vcf.gz.tmp" > "O25799359_consensus2.fasta"
echo "vcf to fasta completed"
