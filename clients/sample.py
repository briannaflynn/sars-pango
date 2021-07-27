address = {
    "name": "Example Laboratory",
    "street": "1 Example St, Ste 333",
    "city": "Houston",
    "state": "TX",
    "zip": "77054",
    "fax": ""
}

clia = {
    "number": "45D2135859",
    "medical_dir": "Example Doctor, M.D., FCAP",
    "lab_dir": "Example Director",
    "lab_notes": ""
}

legal = {
    "disclaimer": "This is an example disclaimer, and requires more information about client",
    "limitations": "This is an example set of limitations, and requires more information about client",
    "methodology": "This is an example set of methods, and requires more information about client",
    "labcert": "{} is regulated under the Clinical Laboratory Improvement Amendments of 1988 and is qualified to perform high-complexity testing. The laboratory is certified by the CLIA program. CLIA # {}".format(address["name"], clia["number"]),
    "confident": 'The information contained in this document is confidential and intended solely for the use of the specific recipient(s) addressed above. To the extent the information in this document contains protected health information, as defined by the Health Insurance Portability and Accountability Act of 1996 ("HIPAA"), 45 C.F.R. Parts 160 & 164, such information is subject to specific privacy and security requirements under HIPAA, and may also be subject to protection under state privacy and privilege laws. Any use, disclosure, dissemination, distribution, or copying of protected health information contained in this document by anyone other than the intended recipient(s) is only permitted in accordance with federal and state law.',
    "recipient": "If you are not the intended recipient, or an authorized representative of the intended recipient, and have received this document in error, please notify the sender immediately by fax or e-mail and delete this document. Do not deliver, distribute, or copy this document, and do not use or disclose its contents or take any action in reliance on the information it contains.",
    "approval": "Approved by {}".format(clia["lab_dir"])
}

logo_path = "clients/Sample_Logo.png"

forward_stranded = True

specimen_type = "Buccal swab"

report_title = "SARS-Cov-2 Analysis"

signature = False

demo_converter = "sample"