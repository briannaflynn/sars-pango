#!/usr/bin/python

import os
import re
import subprocess
import sys
import argparse
import time
import logging
import logging.config
logging.config.fileConfig('logging.conf')
logger = logging.getLogger('main')


def run(a: list):
    args = None
    for j in a:
        args = ' '.join(a)

    return subprocess.run(args, shell=True, stdout=subprocess.PIPE)


def nextrun(a: list):
    args = None
    for j in a:
        args = ' '.join(a)

    args = args + f' > {id}_nextstrain.log'

    return subprocess.run(args, shell=True, stderr=subprocess.PIPE)


def bioinformatix(*args):
    b = os.path.abspath('bioinformatics-pipeline.sh')
    trun = 'bash ' + b
    cmd = [trun] + [*args]

    try:
        process = run(cmd)
    except subprocess.CalledProcessError:
        logger.error(f"Bio pipeline crashed: {process.stderr}")
        raise
    return process.stdout


def pangolin(*args):
    b = os.path.abspath('pango.sh')
    trun = 'bash ' + b
    cmd = [trun] + [*args]

    try:
        process = run(cmd)
    except subprocess.CalledProcessError:
        logger.error(f"Pangolin pipeline crashed: {process.stderr}")
        raise
    return process.stdout


def clade(*args):
    b = os.path.abspath('get_clades.sh')
    trun = 'bash ' + b
    cmd = [trun] + [*args]

    try:
        process = run(cmd)
    except subprocess.CalledProcessError:
        logger.error(f"Get clades crashed: {process.stderr}")
        raise
    return process.stdout


def nextstrain(*args):
    b = os.path.abspath('nextstrain.sh')
    trun = 'bash ' + b
    cmd = [trun] + [*args]

    try:
        process = run(cmd)
    except subprocess.CalledProcessError:
        
        #log.error(f"Nextstrain runner crashed: {process.stderr}")
        #raise
        pass
    # 	print(process.stderr, file=open(f"{id}_nextstrainlog.txt", 'w'))
    return process.stdout


def cleanup(args):
    b = os.path.abspath('wrap-up.sh')
    trun = 'bash ' + b
    cmd = [trun] + [args]

    try:
        process = run(cmd)
    except subprocess.CalledProcessError:
        logger.error(f"Cleanup crashed: {process.stderr}")
        raise
    return process.stdout


def toytree(*args):
    b = os.path.abspath('toytree-run.sh')
    trun = 'bash ' + b
    cmd = [trun] + [*args]

    try:
        process = run(cmd)
    except subprocess.CalledProcessError:
        logger.error(f"Toytree crashed: {process.stderr}")
        raise
    return process.stdout


def worldmap(*args):
    b = os.path.abspath('map-figure.sh')
    trun = 'bash ' + b
    cmd = [trun] + [*args]

    try:
        process = run(cmd)
    except subprocess.CalledProcessError:
        logger.error(f"World map figure crashed: {process.stderr}")
        raise
    return process.stdout


def full_run(id, f1, f2, conda):
    logger.debug("\nStarting full pipeline run ...")
    if f1 and f2 != None:
        logger.debug("\nBeginning bioinformatics pipeline\n")
        bioinformatix(id, f1, f2, conda)
    else:
        raise ValueError("Missing fastq file input! Please provide two (paired) fastq files.")

    logger.debug("\nStarting pangolin analysis ...\n")
    pangolin(id, conda)
    logger.debug("\nExtracting virus clade ...")
    clade(id, conda)
    logger.debug("\nBeginning nextstrain pipeline\nHeads up! This may take awhile ...")
    nextstrain(id, conda)
    logger.debug("\nNextstrain pipeline completed!")
    logger.debug("\nCleaning up ...\n")
    cleanup(id)
    logger.debug("Done!\n")


def visualizer(id, conda):
    logger.debug("\nStarting visualization: Phylogenetics ...")
    toytree(id, conda)
    loc = str(id)
    logger.debug("\nStarting visualization: World Map ...")
    worldmap(id, conda)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="SARS-CoV-2 strain identification and phylogenetic analysis pipeline")
    parser.add_argument('--f1', action='store', type=str,
                        help='Path to fastq.gz 1/2. Required for bioinformatics pipeline, not required for phylogenetic tree visualization and mapping (if the data is already available in id-results directory)')
    parser.add_argument('--f2', action='store', type=str,
                        help='Path to fastq.gz 2/2. Required for bioinformatics pipeline, not required for phylogenetic tree visualization and mapping (if the data is already available in id-results directory)')

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('--id', action='store', type=str, required=True,
                               help='Unique name used to label all files')
    requiredNamed.add_argument('--conda_path', action='store', type=str, required=True,
                               help='Path to anaconda directory required. Eg. /Users/name/anaconda3. To locate your conda path, type `which conda` on your terminal, and copy the path up until `anaconda` or `anaconda3`.')
    args = parser.parse_args()
    fastq_1 = args.f1
    fastq_2 = args.f2
    id = args.id
    conda = args.conda_path

    full_run(id, fastq_1, fastq_2, conda)

    try:
        visualizer(id, conda)
    except FileNotFoundError:
        logger.error(
            f"tree_raw.nwk and clades.json files not found in the expected directory ({id}-results/tree). Check Nextstrain error log and search within directory for these files.")
