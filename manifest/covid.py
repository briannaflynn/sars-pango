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


def current_build(current, drugs):
    """
    take the current medications of the patient and search for them in the drug annotations

    returns { "current": [], "not_found": [<list of chem names or the input ] }
    """
    output = {}
    lcurr = [x.lower() for x in current]
    not_found = deepcopy(lcurr)

    """
    # search for every drug result possible and create a list with it
    for entry in drugs:
        if entry["name"].lower() in lcurr:
            output.append(entry)
            try:
                not_found.remove(entry["name"].lower())
            except:
                logger.debug("Tried to cut {} but it was already gone.".format(entry["name"]))
    logger.debug("The current meds for this sample are {}".format(output))
    """
    # put only the highest guidance levels for each drug into the current meds section
    for entry in drugs:
        if entry["name"].lower() in lcurr or bool(set([j.lower() for j in entry["brand_names_full"]]) & set(lcurr)):
            if entry["name"].lower() not in output.keys():
                output[entry["name"].lower()] = entry
            else:
                if int(output[entry["name"].lower()]["guidanceInt"]) < int(entry["guidanceInt"]):
                    output[entry["name"].lower()] = entry

            try:
                not_found.remove(entry["name"].lower())
            except:
                logger.debug("Tried to cut {} but it was already gone.".format(entry["name"]))

    #turn the output dict into a list
    final_output = output.values()
    logger.debug("The current meds for this sample are {}".format(final_output))
    return {
        "current": final_output,
        "not_found": not_found
    }


def medication_guide(drugs, drug_content):
    """
    Create a dict that is ready to be turned into html

    return {
        category name: {
            drug class: {
                '1': [], '2': [], '3': [],
                'biggest': 0
                }
            }
        }

    """
    output = {}
    for drug in drugs:
        category = drug_content[drug['name']]['category']
        dclass = drug_content[drug['name']]['class']
        logger.debug("Dclass {} category {}".format(dclass, category))
        name = drug['name']
        level = str(drug['guidanceInt'])
        if category not in output.keys():
            output[category] = {}

        # just making sure the dict is ready prior to data entry
        if dclass not in output[category].keys():
            output[category][dclass] = {
                    '1': [],
                    '2': [],
                    '3': [],
                    'biggest': 0
                }
        # if there are repeated chemnames, only keep the one on highest alert
        done = False
        for lev, lis in output[category][dclass].items():
            if lev == "biggest":
                continue

            if int(lev) < int(level):
                if name in lis:
                    output[category][dclass][lev].remove(name)
            else:
                if name in lis:
                    done = True # avoid repeating names in same category

        if done is False:
            output[category][dclass][level].append(name)

    #add a field that contains the length of the longest chem list
    for cat,subs in output.items():
        for fdclass in output[cat].keys():
            maxer = [len(output[cat][fdclass]['1']), len(output[cat][fdclass]['2']), len(output[cat][fdclass]['3'])]
            maxd = max(maxer)
            output[cat][fdclass]['biggest'] = maxd

    logger.debug("The medication_guide dict is {}".format(output))
    return output


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


def prepare_genetics(sample, drugs, diseases, lab_info):
    """

    Go through the drugs and diseases and return the ready genetics dict
    genetics[str_genes] = {
                "info": string,
                "alleles": string
            }

    to-do: more intelligently add CNVs to genotypes

    Rules:
    - if multiple genes
        - genes are not already covered in dict on their own
            - the haplotypes are not numbered
                    -- prepare the info statement as "GENE1 rsX N/N, GENE2 rsX N/N + <phenotype>"
                    -- prepare the alleles tested as list of RSIDs
            - the haplotypes are numbered
                    -- prepare the info statement as "GENE1 N/N, GENE2 N/N + <phenotype>"
                    -- prepare the alleles tested as list of RSIDs
        - at least one of the genes are already covered in keys individually
            -- prepare the info as "<phenotype>"
            -- prepare the alleles tested as "See tests for GENE1, GENEN"
        - the exact combo of genes has already been reported
            - hasn't switched between drugs/diseases
                -- raise error that the same gene combo was found in same annotation group
            - has switched between drugs/diseases
                -- add phenotype to the existing patient info "GENE1 N/N, GENE2 N/N <phenotype>, <phenotype>"
    - if single gene
        - gene is not already covered in dict
            - the haplotypes are not numbered
                -- prepare the info statement as "rs1 N/N, rsX N/N <phenotype>"
                -- prepare the alleles tested as list of RSIDs
            - the haplotypes are numbered
                -- prepare the info statement as "N/N <phenotype>"
                -- prepare the alleles tested as list of RSIDs
        - gene is already covered in dict
            - hasn't switched between drugs/diseases
                -- raise error that the same gene combo was found in same annotation group
            - has switched
                -- add phenotype to the existing patient info "N/N <phenotype>, <phenotype>"
    """
    labgen = {}
    for id, d in lab_info.assays.items():
        if re.search("\*", d[4]):
            added = d[4]
        elif d[2] == "APOE":
            added = d[4]
        elif d[4] == "XN":
            added = "CNV XN"
        else:
            added = d[3]

        if d[2] not in labgen.keys():
            labgen[d[2]] = {
                "alleles": [added]
            }
        else:
            if added not in labgen[d[2]]["alleles"]:
                labgen[d[2]]["alleles"].append(added)

    genetics = {}
    solos = []
    #preconstruct a list of genes that have solo rows, for deciding when to omit gene alleles in multi-gene rows
    for anno in [drugs, diseases]:
        for drug in anno:
            str_genes = ','.join(drug['genes'])
            if str_genes not in solos:
                solos.append(str_genes)

    for anno in [drugs, diseases]:
        for drug in anno:
            str_genes = ', '.join(drug['genes'])
            if str_genes in genetics.keys():
                if anno == diseases:
                    genetics[str_genes]["info"] += drug['phenotype'] + ". "
                continue
            else:
                genetics[str_genes] = {
                    "info": [],
                    "alleles": []
                }

            if len(drug['genes']) > 1:
                adder = ""
                alleles = []
                if "CYP2D6" in drug['genes']:
                    adder = ""
                    alleles = ["See CYP2D6 and CYP2C19"]
                else:

                    for gene in drug['genes']:
                        logger.debug("Sample ")
                        adder += "{} ".format(gene)
                        sample['report_genos'][gene] = sortAlleles(sample['report_genos'][gene])
                        if "singleton" in sample['report_genos'][gene]:
                            for rsid, bases in sample['report_genos'][gene].items():
                                if rsid == "singleton":
                                    continue
                                adder += "{} {} ".format(rsid, "/".join(bases))

                        elif sample['report_genos'][gene] == ["Uncertain genotype."]:
                                adder += "Uncertain genotype, "

                        elif re.search("[0-9]", sample['report_genos'][gene][0]):
                            if re.match("[A-Z]", sample['report_genos'][gene][0]):
                                gt = sample['report_genos'][gene][0] + "/" + sample['report_genos'][gene][1]
                            else:
                                gt = "*" + sample['report_genos'][gene][0] + "/" + "*" + sample['report_genos'][gene][1]

                            adder += "{}, ".format(gt)
                        else: # should only be one gt in gts
                            #if not star alleles, they should have RSIDs in labgen dict
                            rsids = labgen[gene]["alleles"]
                            for rsid in rsids:
                                logger.debug("Sample {} gene {} rsid {}".format(sample['id'], gene, rsid))
                                #check to make sure the RSID didn't drop out along the way
                                if rsid not in sample['rsids'].keys():
                                    continue
                                forward = sample['rsids'][rsid]['gt']
                                final = []
                                if STRANDS[gene] == 'minus':
                                    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
                                    for g in forward:
                                        r = "".join(complement.get(base, base) for base in g)
                                        final.append(r)
                                else:
                                    final = forward

                                adder += "{} {}, ".format(rsid, "/".join(final))

                        if gene in solos:
                            alleles.append("See {}".format(gene))
                            #genetics[str_genes]["alleles"].append("{}".format(gene))
                        else:
                            strung = ', '.join(labgen[gene]["alleles"])
                            alleles.append(strung)
                            #genetics[str_genes]["alleles"].append(strung)

                phenotype = drug['phenotype']
                adder += phenotype + ". "
                genetics[str_genes]["info"] = adder
                genetics[str_genes]["alleles"] = ", ".join(sortAlleles(alleles))

            if len(drug['genes']) == 1:
                adder = ""
                gene = drug['genes'][0]
                #we go through the same complex process as multigene to handle multiple RSID/no stars/single gene cases like MTFHR
                #the phenotype for the solo row on this gene will be whatever comes first, since genes present in the dict already are skipped
                if "singleton" in sample['report_genos'][gene]:
                    for rsid, bases in sample['report_genos'][gene].items():
                        if rsid == "singleton":
                            continue
                        adder += "{} {} ".format(rsid, "/".join(bases))
                elif sample['report_genos'][gene] == ["Uncertain genotype."]:
                    continue
                elif re.search("[0-9]", sample['report_genos'][gene][0]):
                    sample['report_genos'][gene] = sortAlleles(sample['report_genos'][gene])
                    if re.match("[A-Z]", sample['report_genos'][gene][0]):
                        gt = sample['report_genos'][gene][0] + "/" + sample['report_genos'][gene][1]
                    else:
                        gt = "*" + sample['report_genos'][gene][0] + "/" + "*" + sample['report_genos'][gene][1]

                    if gene == "CYP2D6" and int(sample['cnv']) > 2:
                        gt = gt + " XN" + sample['cnv']

                    adder += "{}, ".format(gt)
                else:
                    rsids = labgen[gene]["alleles"]
                    for rsid in rsids:
                        # check to make sure the RSID didn't drop out along the way
                        if rsid not in sample['rsids'].keys():
                            continue
                        forward = sample['rsids'][rsid]['gt']
                        final = []
                        if STRANDS[gene] == 'minus':
                            complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
                            for g in forward:
                                r = "".join(complement.get(base, base) for base in g)
                                final.append(r)
                        else:
                            final = forward

                        adder += "{} {}, ".format(rsid, "/".join(final))

                phenotype = drug['phenotype']
                adder += phenotype + ". "
                genetics[str_genes]["info"] = adder
                genetics[str_genes]["alleles"] = ', '.join(sortAlleles(labgen[gene]["alleles"]))


    #Handling genes that could not be given a genotype despite having markers tested, genes that were N/A
    final_list = []
    for name in genetics.keys():
        if re.search(',', name):
            extender = name.split(', ')
            final_list += extender
        else:
            final_list.append(name)

    for gene, values in sample["report_genos"].items():
        if values == ["Uncertain genotype."]:
            if gene in genetics.keys():
                raise RuntimeError("For sample {} gene {} there was an uncertain genotype but this gene has been used in genetics table already\ngenetics {}".format(sample["id"], gene, genetics))
            genetics[gene] = {
                "info": [],
                "alleles": []
            }
            genetics[gene]["info"] = "Uncertain genotype."
            genetics[gene]["alleles"] = ', '.join(sortAlleles(labgen[gene]["alleles"]))
            final_list.append(gene)
        #now determine if there were tested genes left out of the table because their partner genes were uncertain genotypes, and they didn't have solo rows
        if gene not in final_list:
            adder = ""
            sample['report_genos'][gene] = sortAlleles(sample['report_genos'][gene])
            if "singleton" in sample['report_genos'][gene]:
                for rsid, bases in sample['report_genos'][gene].items():
                    if rsid == "singleton":
                        continue
                    adder += " {} {}".format(rsid, "/".join(bases))
            elif re.search("[0-9]", sample['report_genos'][gene][0]):
                if re.match("[A-Z]", sample['report_genos'][gene][0]):
                    gt = sample['report_genos'][gene][0] + "/" + sample['report_genos'][gene][1]
                else:
                    gt = "*" + sample['report_genos'][gene][0] + "/" + "*" + sample['report_genos'][gene][1]

                if gene == "CYP2D6" and int(sample['cnv']) > 2:
                    gt = gt + " XN" + sample['cnv']

                adder += "{}".format(gt)
            else:
                rsids = labgen[gene]["alleles"]
                for rsid in rsids:
                    # check to make sure the RSID didn't drop out along the way
                    if rsid not in sample['rsids'].keys():
                        continue
                    forward = sample['rsids'][rsid]['gt']
                    final = []
                    if STRANDS[gene] == 'minus':
                        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
                        for g in forward:
                            r = "".join(complement.get(base, base) for base in g)
                            final.append(r)
                    else:
                        final = forward
                    adder += " {} {}".format(rsid, "/".join(final))
            genetics[gene] = {
                "info": [],
                "alleles": []
            }
            adder += ". "
            genetics[gene]["info"] = adder
            genetics[gene]["alleles"] = ', '.join(sortAlleles(labgen[gene]["alleles"]))
            final_list.append(gene)


    return genetics



def med_json(data, pdf_name, lab_info):
    """

    Input JSON-loaded dict, output PDF made from HTML template to file
    JSON file is given as string
    """
    pdfcss_abspath = os.path.abspath("manifest/static/covid_pdf.css")
    css_abspath = os.path.abspath("manifest/static/covid_full.css")
    logo_abspath = os.path.abspath(lab_info.logo_path)
    user_abspath = os.path.abspath("manifest/static/user.svg")
    vial_abspath = os.path.abspath("manifest/static/vial.svg")
    chart_abspath = os.path.abspath("manifest/static/chart-bar.svg")
    alert_icon_abspath = os.path.abspath("manifest/static/exclamation-triangle.svg")
    info_icon_abspath = os.path.abspath("manifest/static/info-circle.svg")
    normal_icon_abspath = os.path.abspath("manifest/static/check-square.svg")


    logger.info("CSS {} LOGO {}".format(css_abspath, logo_abspath))

    if data['kit_type'] == "twentythree":
        data['kit_type'] = "23AndMe"
    elif data['kit_type'] == "ancestry":
        data['kit_type'] = "Ancestry"

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