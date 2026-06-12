import os
import csv
import zipfile

from bigecyhmm import HMM_FOLDER, HMM_TEMPLATE_FILE, PATHWAY_TEMPLATE_FILE, CUSTOM_HYDROGEN_TABLE

def test_template_file_db():
    """ Checks that HMM in hmm folder correspond to HMM in template file.
    """
    hmms_in_db = []
    for hmm_filebasename in os.listdir(HMM_FOLDER):
        hmm_filename = os.path.join(HMM_FOLDER, hmm_filebasename)
        if hmm_filename.endswith('.hmm') and 'check' not in hmm_filename:
            hmm_filebasename = os.path.basename(hmm_filename)
            hmms_in_db.append(hmm_filebasename)

    hmms_in_template = []
    with open(HMM_TEMPLATE_FILE, 'r') as open_hmm_template:
        csvreader = csv.DictReader(open_hmm_template, delimiter='\t')
        for line in csvreader:
            for hmm_file in line['Hmm file'].split(', '):
                hmms_in_template.append(hmm_file)

    assert set(hmms_in_template).issubset(set(hmms_in_db))


def test_pathway_file_db():
    """ Checks that HMM in hmm folder correspond to HMM in template file.
    """
    hmms_in_db = []
    for hmm_filebasename in os.listdir(HMM_FOLDER):
        hmm_filename = os.path.join(HMM_FOLDER, hmm_filebasename)
        if hmm_filename.endswith('.hmm') and 'check' not in hmm_filename:
            hmm_filebasename = os.path.basename(hmm_filename)
            hmms_in_db.append(hmm_filebasename)

    hmms_in_template = []
    with open(HMM_TEMPLATE_FILE, 'r') as open_hmm_template:
        csvreader = csv.DictReader(open_hmm_template, delimiter='\t')
        for line in csvreader:
            for hmm_file in line['Hmm file'].split(', '):
                hmms_in_template.append(hmm_file)

    hmms_in_template_pathway = []
    with open(PATHWAY_TEMPLATE_FILE, 'r') as open_hmm_template:
        csvreader = csv.DictReader(open_hmm_template, delimiter='\t')
        for line in csvreader:
            tmp_pathway_hmms = [hmm.replace('(', '').replace(')', '') for hmm in line['HMMs'].split(' ')]
            hmms_in_template_pathway.extend([hmm for hmm in tmp_pathway_hmms if hmm not in ["and", "or", "not"]])

    assert set(hmms_in_template_pathway).issubset(set(hmms_in_template))

    assert set(hmms_in_template_pathway).issubset(set(hmms_in_db))
