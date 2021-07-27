import click
import re
import csv
import importlib
import processrunner
from manifest import covid_json
import logging
import logging.config
logging.config.fileConfig('logging.conf')
logger = logging.getLogger('main')

"""
First some functions to load demographics and lab data for report 
"""
def demos_load(demos):
    """
    Load the demographics for the input batch, return ordered dict
    Input format expected:
    accession,last,first,dob,provider_id,ordering_physician,medications,icd_codes,collection_date,received_date,gender,contracted_party_name,contracted_party_address,contracted_party_city,contracted_party_state,contracted_party_zip,contracted_party_clia
    Required: Accession, last name, first name, dob
    """
    with open(demos, newline='') as csvfile:
        out_dict = {}
        reader = csv.DictReader(csvfile)
        for row in reader:
            #if row['accession'] == '' or row['last name'] == '' or row['dob'] == '' or row['first name'] == '':
            if row['accession'] == '':
                raise RuntimeError("Missing one of the required demographics in line {}".format(row))
            if row['accession'] in out_dict.keys():
                raise RuntimeError("Identical accession IDs detected at {}".format(row['accession']))

            out_dict[row['accession']] = {
                'last': row['last name'],
                'first': row['first name'],
                'dob': row['dob'],
                'provider_id': row['provider id'],
                'ordering_physician': row['ordering physician'],
                'medications': row['medications'],
                'icd_codes': row['icd codes'],
                'collection_date': row['collection date'],
                'received_date': row['received date'],
                'gender': row['gender'],
                'contracted_party_name': row['account name'],
                'contracted_party_address': row['account address'],
                'contracted_party_city': row['account city'],
                'contracted_party_state': row['account state'],
                'contracted_party_zip': row['account zip'],
                'contracted_party_clia': row['account clia'],
                'order_code': None
            }
            if "order code" in row.keys():
                out_dict[row['accession']]["order_code"] = row["order code"]
        return out_dict


def lab_load(lab_info):
    """
    Load the information for lab white labeling
    See test_data for lab_info.py
    the string should look like a pythong import (test_data.lab_info not test_data/lab_info.py)
    """
    lab_out = importlib.import_module(lab_info.replace('.py', ''))
    global lab_methods
    lab_methods = lab_out
    return lab_out


@click.group()
def cli():
    pass

@cli.command('full_report')
@click.argument('demos', type=click.Path())
@click.argument('lab_info')
@click.argument('fastq', type=click.Path())
@click.argument('fastq2', type=click.Path())
@click.argument('output')
@click.argument('conda')
def full_report(demos, lab_info, fastq, fastq2, output, conda):
    """
    FASTQ and FASTQ2 - Paths to the paired end FASTQ samples for your covid sample
    OUTPUT - the name you would like to designate sample by in output files. Output will be in folder <output>-results. It is expected that the sample name is accession ID,
    and if it is a path, then the string following last / will be this ID.
    CONDA - the path to your anaconda directory (e.g. /Users/myname/anaconda3)
    """
    loaded_demos = demos_load(demos)
    loaded_lab = lab_load(lab_info)
    processrunner.full_run(output, fastq, fastq2, conda)
    processrunner.visualizer(output, conda)
    #processrunner.cleanup(output)
    logger.info("Done!")


@cli.command('visuals')
@click.argument('output')
@click.argument('conda')
def visuals(output, conda):
    """
    Visualizer only. Requires output from full_run prior
    """
    processrunner.visualizer(output, conda)

@cli.command('pdf')
@click.argument('demos', type=click.Path())
@click.argument('lab_info')
@click.argument('output')
def pdf(demos, lab_info, output):
    """
    Run only the pdf for a given sample
    Requires that the output be given as if it were the sample name from the full_report function, NOT the full folder name Example: scratch/SAMPLEID

    Example: strain-reporter.py scratch/example_demos.csv clients.sample scratch/SAMPLE
    """
    loaded_demos = demos_load(demos)
    loaded_lab = lab_load(lab_info)
    sample_name = re.sub(".+\/", "", output)
    if sample_name not in loaded_demos.keys():
        raise KeyError("The sample ID {} was expected to be found in the demographics, and was not".format(sample_name))
    data = {
        "output": output,
        "first_name": loaded_demos[sample_name]["first"],
        "last_name": loaded_demos[sample_name]["last"],
        "collection_date": loaded_demos[sample_name]["collection_date"],
        "received_date": loaded_demos[sample_name]["received_date"],
        "id": sample_name,
        "dob": loaded_demos[sample_name]["dob"],
        "gender": loaded_demos[sample_name]["gender"]
    }
    covid_json(data, output+"-results/{}".format(sample_name), lab_info)


if __name__ == '__main__':
    cli()