import os
import csv
import json
import subprocess
import shutil
import pyhmmer
import zipfile

from bigecyhmm.hmm_search import search_hmm, check_motif_regex, check_motif_pair
from bigecyhmm import HMM_TEMPLATE_FILE, PATHWAY_TEMPLATE_FILE, MOTIF, MOTIF_PAIR, HMM_COMPRESSED_FILE

EXPECTED_RESULTS = {'Q08582': ('Thermophilic specific', None, 'TIGR01054.hmm'), 'P50457': ('4-aminobutyrate aminotransferase and related aminotransferases', 'C-S-01:Organic carbon oxidation', 'K00823.hmm'),
                   'P06292': ('CBB cycle - Rubisco', 'C-S-02:Carbon fixation', 'rubisco_form_I.hmm'), 'P11766': ('Alcohol utilization', 'C-S-03:Ethanol oxidation', 'K00001.hmm'),
                   'Q9NR19': ('Acetate to acetyl-CoA', 'C-S-04:Acetate oxidation', 'TIGR02188.hmm'), 'P29166': ('FeFe hydrogenase', 'C-S-05:Hydrogen generation', 'K00532.hmm'),
                   'P80900': ('Pyruvate oxidation', 'C-S-06:Fermentation', 'K00169.hmm'), 'P11562': ('Methane production', 'C-S-07:Methanogenesis', 'TIGR03259.hmm'),
                   'Q607G3': ('Methane oxidation - Partculate methane monooxygenase', 'C-S-08:Methanotrophy', 'pmoA.hmm'), 'P0ACD8': ('Ni-Fe Hydrogenase', 'C-S-09:Hydrogen oxidation', 'nife-group-1.hmm'),
                   'P13419': ('Wood Ljungdahl pathway (methyl branch)', 'C-S-10:Acetogenesis WL', 'K01938.hmm')}

EXPECTED_FUNCTIONS = {'org_1': ['C-S-01:Organic carbon oxidation', 'C-S-02:Carbon fixation', 'C-S-03:Ethanol oxidation', 'C-S-04:Acetate oxidation', 'C-S-05:Hydrogen generation', 'C-S-06:Fermentation', 'C-S-07:Methanogenesis',
                                'C-S-08:Methanotrophy', 'C-S-09:Hydrogen oxidation', 'C-S-10:Acetogenesis WL', 'P-S-02:Mineralisation'],
                        'org_2': ['C-S-02:Carbon fixation', 'P-S-01:Immobilisation (P-poor)', 'P-S-01:Immobilisation (P-rich)'],
                        'org_3': ['C-S-02:Carbon fixation', 'P-S-01:Immobilisation (P-poor)', 'P-S-01:Immobilisation (P-rich)']}

def extract_hmm_to_function():
    with open(HMM_TEMPLATE_FILE, 'r') as open_hmm_template:
        csvreader = csv.DictReader(open_hmm_template, delimiter='\t')
        hmm_to_function = {}
        for line in csvreader:
            for hmm_file in line['Hmm file'].split(', '):
                function_name = line['Function']
                hmm_to_function[hmm_file] = function_name
    return hmm_to_function


def extract_hmm_to_pathway():
    hmm_to_pathways = {}
    with open(PATHWAY_TEMPLATE_FILE, 'r') as open_r_pathways:
        csvreader = csv.DictReader(open_r_pathways, delimiter = '\t')
        for line in csvreader:
            if '; ' in line['HMMs']:
                for combination in line['HMMs'].split('; '):
                    for hmm_id in combination.split(', '):
                        if hmm_id not in hmm_to_pathways:
                            hmm_to_pathways[hmm_id] = [line['Pathways']]
                        else:
                            hmm_to_pathways[hmm_id].append(line['Pathways'])
            else:
                for hmm_id in line['HMMs'].split(', '):
                    if hmm_id not in hmm_to_pathways:
                        hmm_to_pathways[hmm_id] = [line['Pathways']]
                    else:
                        hmm_to_pathways[hmm_id].append(line['Pathways'])

    return hmm_to_pathways


def test_check_motif():
    gene_name = 'dsrA'
    sequence = 'MSETPLLDELEKGPWPSFVKEIKKTAELMEKAAAEGKDVKMPKGARGLLKQLEISYKDKKTHWKHGGIVSVVGYGGGVIGRYSDLGEQIPEVEHFHTMRINQPSGWFYSTKALRGLCDVWEKWGSGLTNFHGSTGDIIFLGTRSEYLQPCFEDLGNLEIPFDIGGSGSDLRTPSACMGPALCEFACYDTLELCYDLTMTYQDELHRPMWPYKFKIKCAGCPNDCVASKARSDFAIIGTWKDDIKVDQEAVKEYASWMDIENEVVKLCPTGAIKWDGKELTIDNRECVRCMHCINKMPKALKPGDERGATILIGGKAPFVEGAVIGWVAVPFVEVEKPYDEIKEILEAIWDWWDEEGKFRERIGELIWRKGMREFLKVIGREADVRMVKAPRNNPFMFFEKDELKPSAYTEELKKRGMW'
    check_bool = check_motif_regex(gene_name, sequence, MOTIF)
    assert check_bool == True


def test_check_motif_json_file():
    gene_name = 'dsrA'
    sequence = 'MSETPLLDELEKGPWPSFVKEIKKTAELMEKAAAEGKDVKMPKGARGLLKQLEISYKDKKTHWKHGGIVSVVGYGGGVIGRYSDLGEQIPEVEHFHTMRINQPSGWFYSTKALRGLCDVWEKWGSGLTNFHGSTGDIIFLGTRSEYLQPCFEDLGNLEIPFDIGGSGSDLRTPSACMGPALCEFACYDTLELCYDLTMTYQDELHRPMWPYKFKIKCAGCPNDCVASKARSDFAIIGTWKDDIKVDQEAVKEYASWMDIENEVVKLCPTGAIKWDGKELTIDNRECVRCMHCINKMPKALKPGDERGATILIGGKAPFVEGAVIGWVAVPFVEVEKPYDEIKEILEAIWDWWDEEGKFRERIGELIWRKGMREFLKVIGREADVRMVKAPRNNPFMFFEKDELKPSAYTEELKKRGMW'
    custom_motif_file = os.path.join('input_data', 'motif.json')
    with open(custom_motif_file, 'r') as open_motif_json:
        motif_data = json.load(open_motif_json)
    check_bool = check_motif_regex(gene_name, sequence, motif_data)
    assert check_bool == True


def test_check_motif_pair_pmoA():
    gene_name = 'pmoA'
    check_hmms = {'pmoA': 'hmm_files/pmoA.check.hmm', 'amoA': 'hmm_files/amoA.check.hmm'}

    input_protein_fasta = os.path.join('input_data', 'motif_test_data', 'pmoA.fasta')
    with pyhmmer.easel.SequenceFile(input_protein_fasta, digital=True) as seq_file:
        sequences = list(seq_file)
        gene_sequence = [sequence for sequence in sequences]

    first_check_hmm = check_hmms[gene_name]
    second_check_hmm = check_hmms[MOTIF_PAIR[gene_name]]

    with zipfile.ZipFile(HMM_COMPRESSED_FILE, 'r') as zip_object:
        check_bool = check_motif_pair(gene_sequence, first_check_hmm, second_check_hmm, zip_object)
    assert check_bool == True


def test_check_motif_pair_pmoA_json():
    gene_name = 'pmoA'
    check_hmms = {'pmoA': 'hmm_files/pmoA.check.hmm', 'amoA': 'hmm_files/amoA.check.hmm'}

    input_protein_fasta = os.path.join('input_data', 'motif_test_data', 'pmoA.fasta')
    with pyhmmer.easel.SequenceFile(input_protein_fasta, digital=True) as seq_file:
        sequences = list(seq_file)
        gene_sequence = [sequence for sequence in sequences]

    motif_pair_json = os.path.join('input_data', 'motif_pair.json')
    with open(motif_pair_json, 'r') as open_motif_pair_json:
        motif_pair_data = json.load(open_motif_pair_json)

    first_check_hmm = check_hmms[gene_name]
    second_check_hmm = check_hmms[motif_pair_data[gene_name]]

    with zipfile.ZipFile(HMM_COMPRESSED_FILE, 'r') as zip_object:
        check_bool = check_motif_pair(gene_sequence, first_check_hmm, second_check_hmm, zip_object)
    assert check_bool == True


def test_check_motif_pair_pmoA_negative_amoA():
    gene_name = 'amoA'
    check_hmms = {'pmoA': 'hmm_files/pmoA.check.hmm', 'amoA': 'hmm_files/amoA.check.hmm'}

    input_protein_fasta = os.path.join('input_data', 'motif_test_data', 'pmoA.fasta')
    with pyhmmer.easel.SequenceFile(input_protein_fasta, digital=True) as seq_file:
        sequences = list(seq_file)
        gene_sequence = [sequence for sequence in sequences]

    motif_pair_json = os.path.join('input_data', 'motif_pair.json')
    with open(motif_pair_json, 'r') as open_motif_pair_json:
        motif_pair_data = json.load(open_motif_pair_json)

    first_check_hmm = check_hmms[gene_name]
    second_check_hmm = check_hmms[motif_pair_data[gene_name]]

    with zipfile.ZipFile(HMM_COMPRESSED_FILE, 'r') as zip_object:
        check_bool = check_motif_pair(gene_sequence, first_check_hmm, second_check_hmm, zip_object)
    assert check_bool == False


def test_check_motif_pair_pmoA_negative_amoA_json():
    gene_name = 'amoA'
    check_hmms = {'pmoA': 'hmm_files/pmoA.check.hmm', 'amoA': 'hmm_files/amoA.check.hmm'}

    input_protein_fasta = os.path.join('input_data', 'motif_test_data', 'pmoA.fasta')
    with pyhmmer.easel.SequenceFile(input_protein_fasta, digital=True) as seq_file:
        sequences = list(seq_file)
        gene_sequence = [sequence for sequence in sequences]

    first_check_hmm = check_hmms[gene_name]
    second_check_hmm = check_hmms[MOTIF_PAIR[gene_name]]

    with zipfile.ZipFile(HMM_COMPRESSED_FILE, 'r') as zip_object:
        check_bool = check_motif_pair(gene_sequence, first_check_hmm, second_check_hmm, zip_object)
    assert check_bool == False


def test_check_motif_pair_amoA():
    gene_name = 'amoA'
    check_hmms = {'pmoA': 'hmm_files/pmoA.check.hmm', 'amoA': 'hmm_files/amoA.check.hmm'}

    input_protein_fasta = os.path.join('input_data', 'motif_test_data', 'amoA.fasta')
    with pyhmmer.easel.SequenceFile(input_protein_fasta, digital=True) as seq_file:
        sequences = list(seq_file)
        gene_sequence = [sequence for sequence in sequences]

    first_check_hmm = check_hmms[gene_name]
    second_check_hmm = check_hmms[MOTIF_PAIR[gene_name]]

    with zipfile.ZipFile(HMM_COMPRESSED_FILE, 'r') as zip_object:
        check_bool = check_motif_pair(gene_sequence, first_check_hmm, second_check_hmm, zip_object)
    assert check_bool == True


def test_check_motif_pair_amoA_negative_pmoA():
    gene_name = 'pmoA'
    check_hmms = {'pmoA': 'hmm_files/pmoA.check.hmm', 'amoA': 'hmm_files/amoA.check.hmm'}

    input_protein_fasta = os.path.join('input_data', 'motif_test_data', 'amoA.fasta')
    with pyhmmer.easel.SequenceFile(input_protein_fasta, digital=True) as seq_file:
        sequences = list(seq_file)
        gene_sequence = [sequence for sequence in sequences]

    first_check_hmm = check_hmms[gene_name]
    second_check_hmm = check_hmms[MOTIF_PAIR[gene_name]]

    with zipfile.ZipFile(HMM_COMPRESSED_FILE, 'r') as zip_object:
        check_bool = check_motif_pair(gene_sequence, first_check_hmm, second_check_hmm, zip_object)
    assert check_bool == False


def test_search_hmm():
    input_file = os.path.join('input_data', 'meta_organism_test.faa')
    output_folder = 'output_folder'

    search_hmm(input_file, output_folder)

    hmm_to_function = extract_hmm_to_function()
    hmm_to_pathways = extract_hmm_to_pathway()
    predicted_hmm_file = os.path.join(output_folder, 'hmm_results', 'meta_organism_test.tsv')
    predicted_hmms = {}
    with open(predicted_hmm_file, 'r') as open_predicted_hmm_file:
        csvreader = csv.DictReader(open_predicted_hmm_file, delimiter='\t')
        for line in csvreader:
            protein_id = line['protein'].split('|')[1]
            hmm_id = line['HMM']
            if hmm_id in hmm_to_pathways:
                pathways = hmm_to_pathways[hmm_id]
            else:
                pathways = None
            if pathways is not None:
                for pathway in pathways:
                    if protein_id not in predicted_hmms:
                        predicted_hmms[protein_id] = [(hmm_to_function[hmm_id], pathway, hmm_id)]
                    else:
                        predicted_hmms[protein_id].append((hmm_to_function[hmm_id], pathway, hmm_id))
            else:
                if protein_id not in predicted_hmms:
                    predicted_hmms[protein_id] = [(hmm_to_function[hmm_id], pathways, hmm_id)]
                else:
                    predicted_hmms[protein_id].append((hmm_to_function[hmm_id], pathways, hmm_id))

    # Check that expected proteins are matching HMMs.
    for protein_id in EXPECTED_RESULTS:
        assert EXPECTED_RESULTS[protein_id] in predicted_hmms[protein_id]

    shutil.rmtree(output_folder)


def test_search_hmm_cli():
    input_file = os.path.join('input_data', 'meta_organism_test.faa')
    output_folder = 'output_folder'

    cmd = ['bigecyhmm', '-i', input_file, '-o', output_folder]
    subprocess.call(cmd)

    hmm_to_function = extract_hmm_to_function()
    hmm_to_pathways = extract_hmm_to_pathway()
    predicted_hmm_file = os.path.join(output_folder, 'hmm_results', 'meta_organism_test.tsv')
    predicted_hmms = {}
    with open(predicted_hmm_file, 'r') as open_predicted_hmm_file:
        csvreader = csv.DictReader(open_predicted_hmm_file, delimiter='\t')
        for line in csvreader:
            protein_id = line['protein'].split('|')[1]
            hmm_id = line['HMM']
            if hmm_id in hmm_to_pathways:
                pathways = hmm_to_pathways[hmm_id]
            else:
                pathways = None
            if pathways is not None:
                for pathway in pathways:
                    if protein_id not in predicted_hmms:
                        predicted_hmms[protein_id] = [(hmm_to_function[hmm_id], pathway, hmm_id)]
                    else:
                        predicted_hmms[protein_id].append((hmm_to_function[hmm_id], pathway, hmm_id))
            else:
                if protein_id not in predicted_hmms:
                    predicted_hmms[protein_id] = [(hmm_to_function[hmm_id], pathways, hmm_id)]
                else:
                    predicted_hmms[protein_id].append((hmm_to_function[hmm_id], pathways, hmm_id))

    # Check that expected proteins are matching HMMs.
    for protein_id in EXPECTED_RESULTS:
        assert EXPECTED_RESULTS[protein_id] in predicted_hmms[protein_id]

    shutil.rmtree(output_folder)


def test_search_hmm_folder_cli():
    input_file = os.path.join('input_data', 'org_prot')
    output_folder = 'output_folder'

    subprocess.call(['bigecyhmm', '-i', input_file, '-o', output_folder])

    pathway_presence_file = os.path.join(output_folder, 'pathway_presence.tsv')
    pathway_presence_predicted = {}
    pathway_presence_predicted['org_1'] = []
    pathway_presence_predicted['org_2'] = []
    pathway_presence_predicted['org_3'] = []

    with open(pathway_presence_file, 'r') as open_predicted_hmm_file:
        csvreader = csv.DictReader(open_predicted_hmm_file, delimiter='\t')
        for line in csvreader:
            for col_header in line:
                if 'org' in col_header:
                    if line[col_header] == '1':
                        pathway_presence_predicted[col_header].append(line['function'])

    for organism in EXPECTED_FUNCTIONS:
        assert set(EXPECTED_FUNCTIONS[organism]) == set(pathway_presence_predicted[organism])

    shutil.rmtree(output_folder)