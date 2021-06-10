#!/usr/bin/bash

conda activate bioinformatics
#uncommment if you do not have toytree installed
#pip install toytree

# specificy the tree.nwk and clades.json file paths as the first and second positional arguments
tree = $1
clades = $2

python toytree.py tree clades
