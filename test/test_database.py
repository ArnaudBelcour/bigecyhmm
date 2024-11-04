import os
import csv
import zipfile

from bigecyhmm.hmm_search import HMM_COMPRESS_FILE, HMM_TEMPLATE_FILE, PATHWAY_TEMPLATE_FILE

def test_template_file_db():
    """ Checks that HMM in compressed file scorrespond to HMM in template file.
    """
    hmms_in_compress_db = []
    with zipfile.ZipFile(HMM_COMPRESS_FILE, 'r') as zip_object:
        for hmm_filename in zip_object.namelist():
            if hmm_filename.endswith('.hmm') and 'check' not in hmm_filename:
                hmm_filebasename = os.path.basename(hmm_filename)
                hmms_in_compress_db.append(hmm_filebasename)

    hmms_in_template = []
    with open(HMM_TEMPLATE_FILE, 'r') as open_hmm_template:
        csvreader = csv.DictReader(open_hmm_template, delimiter='\t')
        for line in csvreader:
            for hmm_file in line['Hmm file'].split(', '):
                hmms_in_template.append(hmm_file)

    assert set(hmms_in_compress_db) == set(hmms_in_template)


def test_pathway_file_db():
    """ Checks that HMM in compressed file scorrespond to HMM in template file.
    """
    hmms_in_compress_db = []
    with zipfile.ZipFile(HMM_COMPRESS_FILE, 'r') as zip_object:
        for hmm_filename in zip_object.namelist():
            if hmm_filename.endswith('.hmm') and 'check' not in hmm_filename:
                hmm_filebasename = os.path.basename(hmm_filename)
                hmms_in_compress_db.append(hmm_filebasename)

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
            for and_hmms in line['HMMs'].split('; '):
                for hmm_file in and_hmms.split(', '):
                    hmms_in_template_pathway.append(hmm_file.replace('NO|', ''))

    assert set(hmms_in_template_pathway).issubset(set(hmms_in_template))

    assert set(hmms_in_template_pathway).issubset(set(hmms_in_compress_db))