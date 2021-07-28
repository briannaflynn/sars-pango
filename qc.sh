#!/bin/bash

PASSQ=30
EXCLUDE_FLAG=2316

#mean depth
samtools depth -a $1 | awk '{c++;s+=$3}END{print s/c}' > $1.meandepth

#min q30, read or mate unmapped, not primary alignment, and supplementary alignmetns excluded and counted https://broadinstitute.github.io/picard/explain-flags.html
samtools view -q ${PASSQ} -F ${EXCLUDE_FLAG} $1 | wc -l > $1.passingcount

#convert passing reads into bed file, merged
samtools view -bh -q ${PASSQ} -F ${EXCLUDE_FLAG} $1 | bedtools bamtobed -i stdin | bedtools merge -i stdin > $1.bed

#coverage histogram with bedtools
bedtools coverage -hist -b $1 -a $1.bed | grep "all" > $1.coverage

#coverage levels
#bedtools genomecov -ibam $1 -g /Users/gareth/covid/sars-pango/index/Sars2.fna > $1.coverage

#build json
python qc.py $1.meandepth $1.passingcount $1.coverage $1.qc.json