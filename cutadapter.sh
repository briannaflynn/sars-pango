#!/usr/bin/bash

cutadapt -j 12 -g GTTTCCCAGTCACGATA -G GTTTCCCAGTCACGATA -a TATCGTGACTGGGAAAC -A TATCGTGACTGGGAAAC -g ACACTCTTTCCCTACACGACGCTCTTCCGATCT -G ACACTCTTTCCCTACACGACGCTCTTCCGATCT -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -n 3 -m 75 -q 25 --interleaved *"$1_S1_L001_R1_001.fastq.gz" *"$1_S1_L001_R2_001.fastq.gz" > "$1_interleaved.tmp"
cutadapt -j 12 -m 75 -u 25 -u -25 -U 25 -U -25 --interleaved "$1_interleaved.tmp" -o "$1_read1.trim.fastq.tmp" -p "$1_read2.trim.fastq.tmp"
