#!/bin/bash

bcftools call -Mc --ploidy 1 -v --threads 12 -Ov $1.vcf.tmp -o $1_variants.vcf 
