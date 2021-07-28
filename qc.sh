#!/bin/bash

#mean depth
samtools depth -a $1 | awk '{c++;s+=$3}END{print s/c}' > $1.meandepth

#min q30, read or mate unmapped, not primary alignment, and supplementary alignmetns excluded and counted https://broadinstitute.github.io/picard/explain-flags.html
samtools view -q 30 -F 2316 $1 | wc -l > $1.passingcount

#coverage levels
bedtools genomecov -ibam $1 -g /Users/gareth/covid/sars-pango/index/Sars2.fna > $1.coverage

#build json
python qc.py $1.meandepth $1.passingcount $1.coverage