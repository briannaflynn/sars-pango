import click
import re
import processrunner
import logging
import logging.config
logging.config.fileConfig('logging.conf')
logger = logging.getLogger('main')


@click.group()
def cli():
    pass


@cli.command('full_report')
@click.argument('fastq', type=click.Path())
@click.argument('fastq2', type=click.Path())
@click.argument('output')
@click.argument('conda')
def full_report(fastq, fastq2, output, conda):
    """
    FASTQ and FASTQ2 - Paths to the paired end FASTQ samples for your covid sample
    OUTPUT - the name you would like to designate sample by in output files. Output will be in folder <output>-results
    CONDA - the path to your anaconda directory (e.g. /Users/myname/anaconda3)
    """
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

if __name__ == '__main__':
    cli()