#!/usr/bin/bash

# replace name with SAMPLE
# -e is for mac osx version of sed
# for linux, remove the -e flag
var=">SAMPLE"
sed -i '' -e "1s/.*/$var/" $1_consensus2.fasta
