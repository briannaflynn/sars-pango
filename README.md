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

```
# example
python process-runner --id 1234 --f1 1234_R1.fastq.gz --f2 1234_R2.fastq.gz --conda_path /Users/name/anaconda
```
# Results

Results should appear in a directory with the id name in it: id-results.
Ie. if ID is 1234, results directory will be formatted: 1234-results.

Within id-results directory there is:

- lineage_report.csv: from pangolin

- clade_assignment.tsv: for the sample produced from Nextstrain

- id_map.pdf: world map produced from plotly

- id_result.pdf: cropped world map

- genetic_data: all genetic sequence data

- nextstrain_results: raw nextstrain result data

- tree: tree_raw.nwk is the main tree file used to plot phylogeny, clades.json is used to assign clades to branches and nodes, and tree-plot.pdf is a pdf of the phylogenetic tree

#### If using Linux, please see note within the commented portion of fasta-prep.sh

To update the pango-next submodule, follow these steps:
``` 
cd pango-next

git pull origin master

cd ../; git status

if submodule updated, will see:

# Not currently on any branch.
# Changed but not updated:
#   (use "git add ..." to update what will be committed)
#   (use "git checkout -- ..." to discard changes in working directory)
#
#       modified:   pango-next (new commits)

Then run:

git add pango-next

git commit -m "update submodule commit message"
```

# Anaconda Installation

if using Mac OSX: 

``` 
wget https://repo.anaconda.com/archive/Anaconda3-2021.05-MacOSX-x86_64.sh

bash Anaconda3-2021.05-MacOSX-x86_64.sh

# specify the path you want at the end or use the default path provided
# default example:

/Users/garth/anaconda3

# (base) should appear next to your command line prompt indicating you're in the conda base environment
# type the following command to check the install went through and anaconda is installed in the directory you expect

which conda
# this should print out the path to your conda binary file
# example:

/Users/garth/anaconda3/bin/conda
```
If using Linux:

```
# Almost all the same instructions as for Mac OSX, except the path to the Anaconda.sh file
# Copy and paste the command below and run with the Mac OSX instructions

wget https://repo.anaconda.com/archive/Anaconda3-2021.05-Linux-x86_64.sh
```

You can find more information about the Anaconda package installer here: https://www.anaconda.com/products/individual

# Environment Setup

After installing Anaconda, begin setting up all necessary environments. 

Commands needed to set up all virtual environments can be found in the ./Installation/README/ directory in the form of three text files, each corresponding to the three virtual environments that will be used to run the pipeline. 
