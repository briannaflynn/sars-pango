#!/usr/bin/python

import logging
import os
import re
import subprocess
import sys
import argparse
import time

# log = logging.getLogger("bioinformatix")

def run(a:list):
    
    args = None
    for j in a:
    	args = ' '.join(a)
   
    return subprocess.run(args, shell=True, stdout=subprocess.PIPE)
    
def bioinformatix(*args):
	b = os.path.abspath('bioinformatics-pipeline.sh')
	trun = 'bash ' + b
	cmd = [trun] + [*args]
	
	try:
		process = run(cmd)
	except subprocess.CalledProcessError:
		log.error(f"Bio pipeline crashed: {proc.stderr}")
		raise
	return process.stdout
        
def pangolin(args):

	b = os.path.abspath('pango.sh')
	trun = 'bash ' + b
	cmd = [trun] + [args]
		
	try:
		process = run(cmd)
	except subprocess.CalledProcessError:
		log.error(f"Pangolin pipeline crashed: {proc.stderr}")
		raise
	return process.stdout
    
def clade(args):

	b = os.path.abspath('get_clades.sh')
	trun = 'bash ' + b
	cmd = [trun] + [args]
		
	try:
		process = run(cmd)
	except subprocess.CalledProcessError:
		log.error(f"Get clades crashed: {proc.stderr}")
		raise
	return process.stdout        

def nextstrain(args):

	b = os.path.abspath('nextstrain.sh')
	trun = 'bash ' + b
	cmd = [trun] + [args]
		
	try:
		process = run(cmd)
	except subprocess.CalledProcessError:
		log.error(f"Nextstrain runner crashed: {proc.stderr}")
		raise
	return process.stdout
	
def cleanup(args):

	b = os.path.abspath('wrap-up.sh')
	trun = 'bash ' + b
	cmd = [trun] + [args]
		
	try:
		process = run(cmd)
	except subprocess.CalledProcessError:
		log.error(f"Cleanup crashed: {proc.stderr}")
		raise
	return process.stdout
	
def toytree(args):

	b = os.path.abspath('toytree-run.sh')
	trun = 'bash ' + b
	cmd = [trun] + [args]
		
	try:
		process = run(cmd)
	except subprocess.CalledProcessError:
		log.error(f"Toytree crashed: {proc.stderr}")
		raise
	return process.stdout	
	
def worldmap(args):

	b = os.path.abspath('map-figure.sh')
	trun = 'bash ' + b
	cmd = [trun] + [args]
		
	try:
		process = run(cmd)
	except subprocess.CalledProcessError:
		log.error(f"World map figure crashed: {proc.stderr}")
		raise
	return process.stdout	
		
	
	
def full_run(id, f1, f2):

	print("\nStarting full pipeline run ...")
	time.sleep(5)
	
	if f1 and f2 != None:
	
		print("\nBeginning bioinformatics pipeline\n")
		bioinformatix(id, f1, f2)
		
	else:
		
		raise ValueError("Missing fastq file input! Please provide two (paired) fastq files.")
	
	time.sleep(3)
	print("\nStarting pangolin analysis ...\n")
	pangolin(id)
	
	time.sleep(4)
	print("\nExtracting virus clade ...")
	clade(id)
	
	print("\nBeginning nextstrain pipeline\nHeads up! This may take awhile ...")
	
	time.sleep(4)
	nextstrain(id)
	print("\nNextstrain pipeline completed!")
	
	print("\nCleaning up ...\n")
	cleanup(id)
	time.sleep(1)
	
	print("Done!\n")
	
def visualizer(id):
	
# 	print("\nStarting visualization: Phylogeneitcs ...")
# 	time.sleep(3)
# 	
# 	toytree(id)
	loc = str(id)
# 	
# 	print("\nTree pdf now available: " + loc + "-results/tree/tree-plot.pdf")
	
	print("\nStarting visualization: World Map ...")
	time.sleep(3)
	
	worldmap(id)
	print("\nWorld map pdf now available: " + loc + "-results/map.pdf")
	
	 
if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description="SARS-CoV-2 strain identification and phylogenetic analysis pipeline")
	parser.add_argument('--f1', action='store', type=str, help ='Path to fastq.gz 1/2. Required for bioinformatics pipeline, not required for phylogenetic tree visualization and mapping (if the data is already available in id-results directory)')
	parser.add_argument('--f2', action='store', type=str, help ='Path to fastq.gz 2/2. Required for bioinformatics pipeline, not required for phylogenetic tree visualization and mapping (if the data is already available in id-results directory)')
	
	requiredNamed = parser.add_argument_group('required named arguments')
	requiredNamed.add_argument('--id', action='store', type=str, required = True, help ='Unique name used to label all files')
	
	args = parser.parse_args()
	fastq_1 = args.f1
	fastq_2 = args.f2
	id = args.id
	
# 	full_run(id, fastq_1, fastq_2)	
	
	visualizer(id)
	