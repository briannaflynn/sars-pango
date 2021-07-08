# sars-pango

usage: process-runner.py [-h] [--f1 F1] [--f2 F2] --id ID --conda_path
                         CONDA_PATH

SARS-CoV-2 strain identification and phylogenetic analysis pipeline

optional arguments:
  -h, --help            show this help message and exit
  --f1 F1               Path to fastq.gz 1/2. Required for bioinformatics
                        pipeline, not required for phylogenetic tree
                        visualization and mapping (if the data is already
                        available in id-results directory)
  --f2 F2               Path to fastq.gz 2/2. Required for bioinformatics
                        pipeline, not required for phylogenetic tree
                        visualization and mapping (if the data is already
                        available in id-results directory)

required named arguments:
  --id ID               Unique name used to label all files
  --conda_path CONDA_PATH
                        Path to anaconda directory required. Eg.
                        /Users/name/anaconda3. To locate your conda path, type
                        `which conda` on your terminal, and copy the path up
                        until `anaconda` or `anaconda3`.
