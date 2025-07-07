import os
import csv
import subprocess
import networkx as nx
import shutil

from bigecyhmm.custom_db import identify_run_custom_db_search, search_hmm_custom_db
from bigecyhmm import HMM_TEMPLATE_FILE, PATHWAY_TEMPLATE_FILE


EXPECTED_RESULTS = {'Q08582': ('Thermophilic specific', None, 'TIGR01054.hmm'), 'P50457': ('4-aminobutyrate aminotransferase and related aminotransferases', 'C-S-01:Organic carbon oxidation', 'K00823.hmm'),
                   'P06292': ('CBB cycle - Rubisco', 'C-S-02:Carbon fixation', 'rubisco_form_I.hmm'), 'P11766': ('Alcohol utilization', 'C-S-03:Ethanol oxidation', 'K00001.hmm'),
                   'Q9NR19': ('Acetate to acetyl-CoA', 'C-S-04:Acetate oxidation', 'TIGR02188.hmm'), 'P29166': ('FeFe hydrogenase', 'C-S-05:Hydrogen generation', 'K00532.hmm'),
                   'P80900': ('Pyruvate oxidation', 'C-S-06:Fermentation', 'K00169.hmm'), 'P11562': ('Methane production', 'C-S-07:Methanogenesis', 'TIGR03259.hmm'),
                   'Q607G3': ('Methane oxidation - Partculate methane monooxygenase', 'C-S-08:Methanotrophy', 'pmoA.hmm'), 'P0ACD8': ('Ni-Fe Hydrogenase', 'C-S-09:Hydrogen oxidation', 'nife-group-1.hmm'),
                   'P13419': ('Wood Ljungdahl pathway (methyl branch)', 'C-S-10:Acetogenesis WL', 'K01938.hmm')}

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

def test_search_hmm_custom_db():
    input_file = os.path.join('input_data', 'meta_organism_test.faa')
    output_folder = 'output_folder'
    custom_db_json = os.path.join('input_data', 'custom_db', 'carbon_cycle.json')
    custom_motif = os.path.join('input_data', 'motif.json')
    custom_motif_pair = os.path.join('input_data', 'motif_pair.json')

    search_hmm_custom_db(input_file, custom_db_json, output_folder, motif_json=custom_motif, motif_pair_json=custom_motif_pair)

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


def test_search_hmm_custom_db_cli():
    input_file = os.path.join('input_data', 'meta_organism_test.faa')
    output_folder = 'output_folder'
    custom_db_json = os.path.join('input_data', 'custom_db', 'carbon_cycle.json')
    custom_motif = os.path.join('input_data', 'motif.json')
    custom_motif_pair = os.path.join('input_data', 'motif_pair.json')

    subprocess.call(['bigecyhmm_custom', '-i', input_file, '-d', custom_db_json, '-o', output_folder, '-m', custom_motif, '-p', custom_motif_pair])

    hmm_to_function = extract_hmm_to_function()
    hmm_to_pathways = extract_hmm_to_pathway()
    predicted_hmm_file = os.path.join(output_folder, 'carbon_cycle', 'hmm_results', 'meta_organism_test.tsv')
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


def test_search_hmm_custom_db_abundance_cli():
    input_file = os.path.join('input_data', 'org_prot')
    output_folder = 'output_folder'
    custom_db = os.path.join('input_data', 'custom_db')
    custom_motif = os.path.join('input_data', 'motif.json')
    custom_motif_pair = os.path.join('input_data', 'motif_pair.json')
    abundance_file = os.path.join('input_data', 'proteome_tax_id_abundance.tsv')

    subprocess.call(['bigecyhmm_custom', '-i', input_file, '-d', custom_db, '-o', output_folder, '-m', custom_motif, '-p', custom_motif_pair, '--abundance-file', abundance_file])
    expected_abundance = {'Acetate oxidation': {'sample_1': 100.0, 'sample_2': 0.0, 'sample_3': 0.0}, 'Acetogenesis (WL)': {'sample_1': 100.0, 'sample_2': 0.0, 'sample_3': 0.0},
     'Carbon fixation': {'sample_1': 200.0, 'sample_2': 800.0, 'sample_3': 520.0}, 'Ethanol oxidation': {'sample_1': 100.0, 'sample_2': 0.0, 'sample_3': 0.0},
     'Fermentation': {'sample_1': 100.0, 'sample_2': 0.0, 'sample_3': 0.0}, 'Hydrogen generation': {'sample_1': 100.0, 'sample_2': 0.0, 'sample_3': 0.0},
     'Hydrogen oxidation': {'sample_1': 100.0, 'sample_2': 0.0, 'sample_3': 0.0}, 'Methanogenesis': {'sample_1': 100.0, 'sample_2': 0.0, 'sample_3': 0.0},
     'Methanotrophy': {'sample_1': 100.0, 'sample_2': 0.0, 'sample_3': 0.0}, 'Organic carbon oxidation': {'sample_1': 100.0, 'sample_2': 0.0, 'sample_3': 0.0}}

    output_abundance_network_file = os.path.join(output_folder, 'carbon_cycle', 'cycle_diagram_bipartite_abundance.graphml')
    abundance_network = nx.read_graphml(output_abundance_network_file)

    predicted_abundance = {node: abundance_network.nodes[node] for node in abundance_network.nodes}

    for function in expected_abundance:
        for sample in expected_abundance[function]:
            assert expected_abundance[function][sample] == predicted_abundance[function][sample]

    shutil.rmtree(output_folder)


def test_search_hmm_custom_db_measure_cli():
    input_file = os.path.join('input_data', 'org_prot')
    output_folder = 'output_folder'
    custom_db = os.path.join('input_data', 'mini_custom_db')
    custom_motif = os.path.join('input_data', 'motif.json')
    custom_motif_pair = os.path.join('input_data', 'motif_pair.json')
    measure_file = os.path.join('input_data', 'test_measure.tsv')

    subprocess.call(['bigecyhmm_custom', '-i', input_file, '-d', custom_db, '-o', output_folder, '-m', custom_motif, '-p', custom_motif_pair, '--measure-file', measure_file])

    expected_measure = {'Acetate': {'sample_1': 100.0, 'sample_3': 300.0}, 'H2': {'sample_1': 0.0, 'sample_2': 500.0, 'sample_3': 50.0}}

    output_abundance_network_file = os.path.join(output_folder, 'carbon_cycle', 'cycle_diagram_bipartite_abundance.graphml')
    abundance_network = nx.read_graphml(output_abundance_network_file)

    predicted_abundance = {node: abundance_network.nodes[node] for node in abundance_network.nodes}

    for metabolite in expected_measure:
        for sample in expected_measure[metabolite]:
            assert expected_measure[metabolite][sample] == predicted_abundance[metabolite][sample]

    shutil.rmtree(output_folder)


def test_search_hmm_custom_db_abundance_measure_cli():
    input_file = os.path.join('input_data', 'org_prot')
    output_folder = 'output_folder'
    custom_db = os.path.join('input_data', 'mini_custom_db')
    custom_motif = os.path.join('input_data', 'motif.json')
    custom_motif_pair = os.path.join('input_data', 'motif_pair.json')
    abundance_file = os.path.join('input_data', 'proteome_tax_id_abundance.tsv')
    measure_file = os.path.join('input_data', 'test_measure.tsv')

    subprocess.call(['bigecyhmm_custom', '-i', input_file, '-d', custom_db, '-o', output_folder, '-m', custom_motif, '-p', custom_motif_pair, '--abundance-file', abundance_file,
                     '--measure-file', measure_file])
    expected_abundance = {'Acetate oxidation': {'sample_1': 100.0, 'sample_2': 0.0, 'sample_3': 0.0}, 'Acetogenesis (WL)': {'sample_1': 100.0, 'sample_2': 0.0, 'sample_3': 0.0},
     'Carbon fixation': {'sample_1': 200.0, 'sample_2': 800.0, 'sample_3': 520.0}, 'Ethanol oxidation': {'sample_1': 100.0, 'sample_2': 0.0, 'sample_3': 0.0},
     'Fermentation': {'sample_1': 100.0, 'sample_2': 0.0, 'sample_3': 0.0}, 'Hydrogen generation': {'sample_1': 0.0, 'sample_2': 0.0, 'sample_3': 0.0},
     'Hydrogen oxidation': {'sample_1': 0.0, 'sample_2': 0.0, 'sample_3': 0.0}, 'Methanogenesis': {'sample_1': 0.0, 'sample_2': 0.0, 'sample_3': 0.0},
     'Methanotrophy': {'sample_1': 100.0, 'sample_2': 0.0, 'sample_3': 0.0}, 'Organic carbon oxidation': {'sample_1': 100.0, 'sample_2': 0.0, 'sample_3': 0.0}}

    expected_measure = {'Acetate': {'sample_1': 100.0, 'sample_3': 300.0}, 'H2': {'sample_1': 0.0, 'sample_2': 500.0, 'sample_3': 50.0}}

    output_abundance_network_file = os.path.join(output_folder, 'carbon_cycle', 'cycle_diagram_bipartite_abundance.graphml')
    abundance_network = nx.read_graphml(output_abundance_network_file)

    predicted_abundance = {node: abundance_network.nodes[node] for node in abundance_network.nodes}

    for function in expected_abundance:
        for sample in expected_abundance[function]:
            assert expected_abundance[function][sample] == predicted_abundance[function][sample]

    for metabolite in expected_measure:
        for sample in expected_measure[metabolite]:
            assert expected_measure[metabolite][sample] == predicted_abundance[metabolite][sample]

    output_abundance_network_file = os.path.join(output_folder, 'phosphorus_cycle', 'cycle_diagram_bipartite_abundance.graphml')
    abundance_network = nx.read_graphml(output_abundance_network_file)

    expected_abundance = {'Immobilisation (P-rich)': {'sample_1': 100.0, 'sample_2': 800.0, 'sample_3': 520.0}, 'Immobilisation (P-poor)': {'sample_1': 100.0, 'sample_2': 800.0, 'sample_3': 520.0},
     'Mineralisation': {'sample_1': 100.0, 'sample_2': 0.0, 'sample_3': 0.0}}

    predicted_abundance = {node: abundance_network.nodes[node] for node in abundance_network.nodes}

    for function in expected_abundance:
        for sample in expected_abundance[function]:
            assert expected_abundance[function][sample] == predicted_abundance[function][sample]

    shutil.rmtree(output_folder)


def test_search_hmm_custom_db_abundance_measure_esmecata_cli():
    input_file = os.path.join('input_data', 'esmecata_output_folder', '1_clustering', 'reference_proteins_consensus_fasta')
    output_folder = 'output_folder'
    custom_db = os.path.join('input_data', 'mini_custom_db')
    custom_motif = os.path.join('input_data', 'motif.json')
    custom_motif_pair = os.path.join('input_data', 'motif_pair.json')
    abundance_file = os.path.join('input_data', 'proteome_tax_id_abundance.tsv')
    measure_file = os.path.join('input_data', 'test_measure.tsv')
    esmecata_folder = os.path.join('input_data', 'esmecata_output_folder')

    subprocess.call(['bigecyhmm_custom', '-i', input_file, '-d', custom_db, '-o', output_folder, '-m', custom_motif, '-p', custom_motif_pair, '--abundance-file', abundance_file,
                     '--measure-file', measure_file, '--esmecata', esmecata_folder])

    expected_abundance = {'Acetate oxidation': {'sample_1': 100.0, 'sample_2': 0.0, 'sample_3': 0.0}, 'Acetogenesis (WL)': {'sample_1': 100.0, 'sample_2': 0.0, 'sample_3': 0.0},
     'Carbon fixation': {'sample_1': 200.0, 'sample_2': 800.0, 'sample_3': 520.0}, 'Ethanol oxidation': {'sample_1': 100.0, 'sample_2': 0.0, 'sample_3': 0.0},
     'Fermentation': {'sample_1': 100.0, 'sample_2': 0.0, 'sample_3': 0.0}, 'Hydrogen generation': {'sample_1': 0.0, 'sample_2': 0.0, 'sample_3': 0.0},
     'Hydrogen oxidation': {'sample_1': 0.0, 'sample_2': 0.0, 'sample_3': 0.0}, 'Methanogenesis': {'sample_1': 0.0, 'sample_2': 0.0, 'sample_3': 0.0},
     'Methanotrophy': {'sample_1': 100.0, 'sample_2': 0.0, 'sample_3': 0.0}, 'Organic carbon oxidation': {'sample_1': 100.0, 'sample_2': 0.0, 'sample_3': 0.0}}

    expected_measure = {'Acetate': {'sample_1': 100.0, 'sample_3': 300.0}, 'H2': {'sample_1': 0.0, 'sample_2': 500.0, 'sample_3': 50.0}}

    output_abundance_network_file = os.path.join(output_folder, 'carbon_cycle', 'cycle_diagram_bipartite_abundance.graphml')
    abundance_network = nx.read_graphml(output_abundance_network_file)

    predicted_abundance = {node: abundance_network.nodes[node] for node in abundance_network.nodes}

    for function in expected_abundance:
        for sample in expected_abundance[function]:
            assert expected_abundance[function][sample] == predicted_abundance[function][sample]

    for metabolite in expected_measure:
        for sample in expected_measure[metabolite]:
            assert expected_measure[metabolite][sample] == predicted_abundance[metabolite][sample]

    output_abundance_network_file = os.path.join(output_folder, 'phosphorus_cycle', 'cycle_diagram_bipartite_abundance.graphml')
    abundance_network = nx.read_graphml(output_abundance_network_file)

    expected_abundance = {'Immobilisation (P-rich)': {'sample_1': 100.0, 'sample_2': 800.0, 'sample_3': 520.0}, 'Immobilisation (P-poor)': {'sample_1': 100.0, 'sample_2': 800.0, 'sample_3': 520.0},
     'Mineralisation': {'sample_1': 100.0, 'sample_2': 0.0, 'sample_3': 0.0}}

    predicted_abundance = {node: abundance_network.nodes[node] for node in abundance_network.nodes}

    for function in expected_abundance:
        for sample in expected_abundance[function]:
            assert expected_abundance[function][sample] == predicted_abundance[function][sample]

    shutil.rmtree(output_folder)