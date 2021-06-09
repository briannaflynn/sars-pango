#!/usr/bin/bash

conda activate nextstrain

snakemake --cores 1 --profile pango-next/my_profiles/example -p
