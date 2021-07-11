#!/usr/bin/bash

# replace name with SAMPLE
# -e is for mac osx version of sed
# for linux, remove the -e flag
var=">SAMPLE"
sed -i '' -e "1s/.*/$var/" $1_consensus2.fasta

# sed is used here to replace the first line of the $1_consensus2.fasta file with ">SAMPLE"

# for Linux, use this command:
# sed -i "1s/.*/$var/" $1_consensus2.fasta
