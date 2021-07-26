from weasyprint import HTML, CSS
from weasyprint.pdf import PDFFile, pdf_format
from io import BytesIO
import json
import jinja2
import os
import re
from datetime import date
from operator import itemgetter
from copy import deepcopy
import logging
import logging.config

logging.config.fileConfig('logging.conf')
logger = logging.getLogger('main')

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



def med_json(data, pdf_name, lab_info):
    """

    Input JSON-loaded dict, output PDF made from HTML template to file
    JSON file is given as string
    """
    pdfcss_abspath = os.path.abspath("manifest/static/covid_pdf.css")
    css_abspath = os.path.abspath("manifest/static/covid_full.css")
    logo_abspath = os.path.abspath(lab_info.logo_path)
    logger.info("CSS {} LOGO {}".format(css_abspath, logo_abspath))

    #paths to sample files
    sample_name = data["output"].split("\/")[-1]
    results_dir = data["output"]
    vcf_file = os.path.abspath(results_dir + "/genetic_data/" + sample_name + "_variants.vcf")
    map_file = os.path.abspath(results_dir + "/{}_result.png".format(sample_name))
    tree_file = os.path.abspath(results_dir )
    #load drug content so we can add brand names and other info in
    drug_dir = content_dir + "/drugs"
    drug_content = {}
    for i in os.scandir(drug_dir):
        if i.path.endswith('json') and i.is_file():
            content = json.load(open(i, "r"))
            if content['chemName'] in drug_content.keys():
                logger.error("For some reason drug {} is present as chem name more than once.".format(content['chemName']))
                raise RuntimeError
            else:
                drug_content[content['chemName']] = content

    #prepare medication recommendations
    drugs = []
    genetics = {}
    alerts_count = 0
    info_count = 0
    normal_count = 0
    for rec in data['combined']['drugs']:
        #check to see if the gene is present in the scope of the order code
        if data["order_code"] is not None and lab_info.order_codes[data["order_code"]]["genes"] != "full":
            track = True
            for gene in rec['genes']:
                if gene not in lab_info.order_codes[data["order_code"]]["genes"]:
                    track = False
            if track is False:
                continue
        if rec['guidanceLevel'] == "1":
            rec['guidanceLevel'] = normal_icon_abspath
            rec['guidanceInt'] = 1
            normal_count += 1
        elif rec['guidanceLevel'] == "2":
            rec['guidanceLevel'] = info_icon_abspath
            rec['guidanceInt'] = 2
            info_count += 1
        elif rec['guidanceLevel'] == "3":
            rec['guidanceLevel'] = alert_icon_abspath
            rec['guidanceInt'] = 3
            alerts_count += 1
        else:
            logger.error("Not a valid guidance level {}".format(rec['guidanceLevel']))
            raise RuntimeError

        # if brand names become too long, cut them off for the record
        if len(drug_content[rec['name']]['brandNames']) > 3:
            branders = drug_content[rec['name']]['brandNames']
            rec['brand_names'] = ', '.join(branders[0:4])
        else:
            rec['brand_names'] = ', '.join(drug_content[rec['name']]['brandNames'])

        rec['brand_names_full'] = ', '.join(drug_content[rec['name']]['brandNames'])
        # pull out genotypes for the table
        gts = []
        for gene in rec['genes']:
            if "singleton" in data["report_genos"][gene]:
                for key, val in data["report_genos"][gene].items():
                    if key != "singleton":
                        gts.append({key: val})
                gts.append(data["report_genos"][gene])
                continue
            gt = '/'.join(data["report_genos"][gene])
            gts.append(gt)
        rec['gts'] = gts
        rec['genelist'] = ", ".join(rec['genes'])
        #use the drug content dict to pull in the category
        rec['category'] = drug_content[rec['name']]['category']

        rec['evidenceLevel'] = re.sub('evidence', '', rec['evidenceLevel'] )
        drugs.append(rec)

    #prepare disease management recommendations
    diseases = []
    for rec in data['combined']['diseases']:
        # check to see if the gene is present in the scope of the order code
        if data["order_code"] is not None and lab_info.order_codes[data["order_code"]]["genes"] != "full":
            track = True
            for gene in rec['genes']:
                if gene not in lab_info.order_codes[data["order_code"]]["genes"]:
                    track = False
            if track is False:
                continue
        if rec['guidanceLevel'] == "1":
            rec['guidanceLevel'] = normal_icon_abspath
            rec['guidanceInt'] = 1
            normal_count += 1
        elif rec['guidanceLevel'] == "2":
            rec['guidanceLevel'] = info_icon_abspath
            rec['guidanceInt'] = 2
            info_count += 1
        elif rec['guidanceLevel'] == "3":
            rec['guidanceLevel'] = alert_icon_abspath
            rec['guidanceInt'] = 3
            alerts_count += 1
        else:
            logger.error("Not a valid guidance level {}".format(rec['guidanceLevel']))
            raise RuntimeError

        # pull out genotypes for the table
        gts = []
        for gene in rec['genes']:
            if "singleton" in data["report_genos"][gene]:
                for key, val in data["report_genos"][gene].items():
                    if key != "singleton":
                        gts.append({key: val})
                gts.append(data["report_genos"][gene])
                continue
            gt = '/'.join(data["report_genos"][gene])
            gts.append(gt)
        rec['gts'] = gts
        rec['genelist'] = ", ".join(rec['genes'])
        rec['evidenceLevel'] = re.sub('evidence', '', rec['evidenceLevel'])
        diseases.append(rec)

    #prepare the genetics data
    genetics = prepare_genetics(data, drugs, diseases, lab_info)

    #sort drugs by alphabetical order
    drugs = sorted(drugs, key=itemgetter('name'))
    diseases = sorted(diseases, key=itemgetter('name'))

    # load the current medications
    currents = current_build(data["current_meds"], drugs)
    current_drugs = currents["current"]
    medtable = medication_guide(drugs, drug_content)
    templateLoader = jinja2.FileSystemLoader(searchpath="./manifest/static/")
    templateEnv = jinja2.Environment(loader=templateLoader)
    TEMPLATE_FILE = "med_full.html"
    template = templateEnv.get_template(TEMPLATE_FILE)
    if len(lab_info.order_codes.keys()) > 0:
        title = lab_info.order_codes[data["order_code"]]["title"]
    else:
        title = lab_info.report_title
    output_html = template.render(css_file=css_abspath, lablogo=logo_abspath, drugs=drugs, first_name=data['first_name'], last_name=data['last_name'], speciment_type=lab_info.specimen_type, collection_date=data["collection_date"], received_date=data["received_date"], report_date=str(date.today().strftime("%m/%d/%Y")), sample_id=data["id"], dob=data["dob"], sex=data["gender"], lab_name=lab_info.address["name"], lab_street=lab_info.address["street"], lab_city=lab_info.address["city"], lab_state=lab_info.address["state"], lab_zip=lab_info.address["zip"], report_title=title, lab_dir=lab_info.clia["lab_dir"], gts=genetics, user_svg=user_abspath, vial_svg=vial_abspath, chartbar_svg=chart_abspath, alert_icon_svg=alert_icon_abspath, info_icon_svg=info_icon_abspath, normal_icon_svg=normal_icon_abspath, current_drugs=current_drugs, medication_guide=medtable, sigpresent=lab_info.signature, diseases=diseases, disclaimer=lab_info.legal["disclaimer"], method=lab_info.legal["methodology"], limitations=lab_info.legal["limitations"], cert=lab_info.legal["labcert"], confident=lab_info.legal["confident"], recipient=lab_info.legal["recipient"], approved=lab_info.legal["approval"])

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