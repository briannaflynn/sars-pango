from weasyprint import HTML, CSS
from weasyprint.pdf import PDFFile, pdf_format
from io import BytesIO
import json
import jinja2
import os
import re
import csv
from datetime import date
from operator import itemgetter
from copy import deepcopy
import logging
import logging.config

logging.config.fileConfig('logging.conf')
logger = logging.getLogger('main')

taxon = "NC_045512.2"


def sortAlleles(genotype):
    """

    Sort the genotype alleles and return the sorted list

    """
    #res = sorted(genotype, key=lambda ele: (0, int(re.sub("\*", "", ele))) if re.sub("\*", "", ele).isdigit() else (1, ele))
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    res = sorted(genotype, key=alphanum_key)
    #print("INPUT {} SORTER {}".format(genotype, res))
    return res

def process_vcf(vcf):
    """
    Parse VCF file and run anything necessary to produce the variants table:
    Dict with position as key and value as dict {'reference': str, 'sample': str}
    """
    output = {}
    with open(vcf, 'r') as f:
        for line in f.readlines():
            if line.startswith("#"):
                continue
            line = line.strip()
            line = line.split()
            position = line[1]
            reference = line[3]
            alleles = [reference]
            alleles.extend(line[4].split(","))
            gt = line[-1].split(":")[0]
            sample = alleles[int(gt)]
            if position in output.keys():
                raise KeyError("Position {} already existed when trying to create variants table. You may need to look at the VCF for an edge case.")

            output[position] = {
                "reference": reference,
                "sample": sample
            }
    return output

def process_results(pangolin, clades):
    """
    Parse the file called "lineage_report.csv" in results folder
    Return dict
    """
    output = {}
    with open(clades) as c:
        for line in c.readlines():
            line = line.strip()
            line = line.split("\t")
            output["clade"] = line[1]
            output["parent clade"] = line[2]

    with open(pangolin, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row["status"] != "passed_qc":
                raise RuntimeError("This sample {} did not pass internal qc for N content 50% and min length (10kbp), as required by the pangolin algorithm.")
            output = {
                "lineage": row["lineage"],
                "taxon": taxon,
                "confidence": row["ambiguity_score"],
                "qc": "",
                "other": re.sub("scorpio call:", "Functional alleles:", row["note"])
            }
    return output


def covid_json(data, pdf_name, lab_info):
    """

    Input JSON-loaded dict, output PDF made from HTML template to file
    JSON file is given as string
    """
    pdfcss_abspath = os.path.abspath("manifest/static/covid_pdf.css")
    css_abspath = os.path.abspath("manifest/static/covid_full.css")
    logo_abspath = os.path.abspath(lab_info.logo_path)
    logger.info("CSS {} LOGO {}".format(css_abspath, logo_abspath))

    #paths to sample files
    sample_name = data["output"].split("/")[-1]
    results_dir = data["output"] + "-results"
    vcf_file = os.path.abspath(results_dir + "/genetic_data/" + sample_name + "_variants.vcf")
    map_file = os.path.abspath(results_dir + "/{}_result.png".format(sample_name))
    tree_file = os.path.abspath(results_dir + "/tree/tree-plot.png")
    qc_file = os.path.abspath(results_dir + "/qc.json")
    legend_file = os.path.abspath("manifest/static/vertical_legend.png")
    results_file = os.path.abspath(results_dir + "/lineage_report.csv")
    clades_file = os.path.abspath(results_dir + "/clade_assignment.tsv")

    #load qc json
    qc = json.load(open(qc_file, 'r'))

    #results summaries
    print(clades_file)
    summary = process_results(results_file, clades_file)
    summary["qc"] = qc["qc"]

    #load variants dict
    vardict = process_vcf(vcf_file)

    templateLoader = jinja2.FileSystemLoader(searchpath="./manifest/static/")
    templateEnv = jinja2.Environment(loader=templateLoader)
    TEMPLATE_FILE = "covid_full.html"
    template = templateEnv.get_template(TEMPLATE_FILE)
    title = lab_info.report_title
    output_html = template.render(css_file=css_abspath,
                                  lablogo=logo_abspath,
                                  first_name=data['first_name'],
                                  last_name=data['last_name'],
                                  speciment_type=lab_info.specimen_type,
                                  collection_date=data["collection_date"],
                                  received_date=data["received_date"],
                                  report_date=str(date.today().strftime("%m/%d/%Y")),
                                  sample_id=data["id"],
                                  dob=data["dob"],
                                  sex=data["gender"],
                                  lab_name=lab_info.address["name"],
                                  lab_street=lab_info.address["street"],
                                  lab_city=lab_info.address["city"],
                                  lab_state=lab_info.address["state"],
                                  lab_zip=lab_info.address["zip"],
                                  report_title=title,
                                  lab_dir=lab_info.clia["lab_dir"],
                                  disclaimer=lab_info.legal["disclaimer"],
                                  method=lab_info.legal["methodology"],
                                  limitations=lab_info.legal["limitations"],
                                  cert=lab_info.legal["labcert"],
                                  confident=lab_info.legal["confident"],
                                  recipient=lab_info.legal["recipient"],
                                  approved=lab_info.legal["approval"],
                                  qc=qc,
                                  summary=summary,
                                  gts=vardict,
                                  tree_path=tree_file,
                                  legend_path=legend_file,
                                  map_path=map_file
                                  )

    #write html to /tmp
    logger.debug("writing HTML")
    h = open(pdf_name + ".html", 'w')
    h.write(output_html)
    html_abspath = os.path.abspath(pdf_name+".html")
    #html_abspath = pdf_name + ".html"
    logger.debug("html {}".format(html_abspath))


    #weasyprint begins here
    #start loading HTML
    html = HTML(string=output_html)
    logger.debug("HTML loaded")
    css = CSS(pdfcss_abspath)
    css_main = CSS(css_abspath)
    logger.debug("CSS loaded")


    # write the pdf from the html
    content = BytesIO(html.write_pdf(stylesheets=[css, css_main]))
    pdf_file = PDFFile(content)
    params = pdf_format('/OpenAction [0 /FitV null]')
    pdf_file.extend_dict(pdf_file.catalog, params)
    pdf_file.finish()
    pdf = pdf_file.fileobj.getvalue()
    open(pdf_name + '.pdf', 'wb').write(pdf)


def weasyprint(file_html):
    """
    input html, output pdf with weasyprint engine
    """
    pdfcss_abspath = os.path.abspath("manifest/static/med_pdf.css")
    css_abspath = os.path.abspath("manifest/static/med_pgx.css")
    pdf_name = re.sub("\.html", "", file_html)

    # weasyprint begins here
    # start loading HTML
    html = HTML(filename=file_html)
    logger.debug("HTML loaded")
    css = CSS(pdfcss_abspath)
    css_main = CSS(css_abspath)
    logger.debug("CSS loaded")
    # write the pdf from the html
    content = BytesIO(html.write_pdf(stylesheets=[css, css_main]))
    pdf_file = PDFFile(content)
    params = pdf_format('/OpenAction [0 /FitV null]')
    pdf_file.extend_dict(pdf_file.catalog, params)
    pdf_file.finish()
    pdf = pdf_file.fileobj.getvalue()
    open(pdf_name + '.pdf', 'wb').write(pdf)