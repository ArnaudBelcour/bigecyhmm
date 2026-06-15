import os
import pandas as pd
import subprocess
import shutil
import pytest
import networkx as nx

from bigecyhmm.visualisation import create_visualisation


def test_create_visualisation():
    bigecyhmm_output_folder = os.path.join('input_data', 'bigecyhmm_output_folder')

    output_folder = 'output_folder'
    create_visualisation(bigecyhmm_output_folder, output_folder)

    heatmap_occurrence = os.path.join(output_folder, 'function_occurrence', 'heatmap_occurrence.png')
    assert os.path.exists(heatmap_occurrence)

    shutil.rmtree(output_folder)


def test_create_visualisation_cli():
    bigecyhmm_output_folder = os.path.join('input_data', 'bigecyhmm_output_folder')

    output_folder = 'output_folder'
    subprocess.call(['bigecyhmm_visualisation', 'genomes', '--bigecyhmm', bigecyhmm_output_folder, '-o', output_folder])

    heatmap_occurrence = os.path.join(output_folder, 'function_occurrence', 'heatmap_occurrence.png')
    assert os.path.exists(heatmap_occurrence)

    shutil.rmtree(output_folder)


def test_create_visualisation_ko_cli():
    ko_abundance_file = os.path.join('input_data', 'ko_file.tsv')
    output_folder = 'output_folder'

    subprocess.call(['bigecyhmm_visualisation', 'ko', '--ko', ko_abundance_file, '-o', output_folder])

    polar_plot = os.path.join(output_folder, 'diagram_visualisation', 'sample_1_carbon_cycle.png')
    assert os.path.exists(polar_plot)
    heatmap_occurrence = os.path.join(output_folder, 'diagram_visualisation', 'sample_2_nitrogen_cycle.png')
    assert os.path.exists(heatmap_occurrence)

    shutil.rmtree(output_folder)


def test_create_visualisation_abundance():
    abundance_file_path = os.path.join('input_data', 'abundance_file_from_genomes.tsv')
    bigecyhmm_output_folder = os.path.join('input_data', 'bigecyhmm_output_folder')

    output_folder = 'output_folder'
    create_visualisation(bigecyhmm_output_folder, output_folder, abundance_file_path=abundance_file_path)

    carbon_cycle_diagram = os.path.join(output_folder, 'function_abundance', 'cycle_diagrams_abundance', 'sample_1_carbon_cycle.png')
    assert os.path.exists(carbon_cycle_diagram)
    polar_plot = os.path.join(output_folder, 'function_abundance', 'polar_plot_abundance', 'polar_plot_abundance_sample_sample_1.png')
    assert os.path.exists(polar_plot)
    polar_plot = os.path.join(output_folder, 'function_abundance', 'polar_plot_abundance', 'polar_plot_abundance_sample_sample_2.png')
    assert os.path.exists(polar_plot)
    polar_plot = os.path.join(output_folder, 'function_abundance', 'polar_plot_abundance', 'polar_plot_abundance_sample_sample_3.png')
    assert os.path.exists(polar_plot)
    heatmap_occurrence = os.path.join(output_folder, 'function_occurrence', 'heatmap_occurrence.png')
    assert os.path.exists(heatmap_occurrence)

    shutil.rmtree(output_folder)


def test_create_visualisation_abundance_cli():
    abundance_file_path = os.path.join('input_data', 'abundance_file_from_genomes.tsv')
    bigecyhmm_output_folder = os.path.join('input_data', 'bigecyhmm_output_folder')

    output_folder = 'output_folder'
    subprocess.call(['bigecyhmm_visualisation', 'genomes', '--bigecyhmm', bigecyhmm_output_folder, '--abundance-file', abundance_file_path, '-o', output_folder])

    carbon_cycle_diagram = os.path.join(output_folder, 'function_abundance', 'cycle_diagrams_abundance', 'sample_1_carbon_cycle.png')
    assert os.path.exists(carbon_cycle_diagram)
    polar_plot = os.path.join(output_folder, 'function_abundance', 'polar_plot_abundance', 'polar_plot_abundance_sample_sample_1.png')
    assert os.path.exists(polar_plot)
    polar_plot = os.path.join(output_folder, 'function_abundance', 'polar_plot_abundance', 'polar_plot_abundance_sample_sample_2.png')
    assert os.path.exists(polar_plot)
    polar_plot = os.path.join(output_folder, 'function_abundance', 'polar_plot_abundance', 'polar_plot_abundance_sample_sample_3.png')
    assert os.path.exists(polar_plot)
    heatmap_occurrence = os.path.join(output_folder, 'function_occurrence', 'heatmap_occurrence.png')
    assert os.path.exists(heatmap_occurrence)

    shutil.rmtree(output_folder)


def test_create_visualisation_abundance_from_esmecata_cli():
    esmecata_output_folder = os.path.join('input_data', 'esmecata_output_folder')
    abundance_file_path = os.path.join('input_data', 'proteome_tax_id_abundance.tsv')
    bigecyhmm_output_folder = os.path.join('input_data', 'bigecyhmm_output_folder')

    output_folder = 'output_folder'
    subprocess.call(['bigecyhmm_visualisation', 'esmecata', '--esmecata', esmecata_output_folder, '--bigecyhmm', bigecyhmm_output_folder, '--abundance-file', abundance_file_path, '-o', output_folder])

    carbon_cycle_diagram = os.path.join(output_folder, 'function_abundance', 'cycle_diagrams_abundance', 'sample_1_carbon_cycle.png')
    assert os.path.exists(carbon_cycle_diagram)
    polar_plot = os.path.join(output_folder, 'function_abundance', 'polar_plot_abundance', 'polar_plot_abundance_sample_sample_1.png')
    assert os.path.exists(polar_plot)
    polar_plot = os.path.join(output_folder, 'function_abundance', 'polar_plot_abundance', 'polar_plot_abundance_sample_sample_2.png')
    assert os.path.exists(polar_plot)
    polar_plot = os.path.join(output_folder, 'function_abundance', 'polar_plot_abundance', 'polar_plot_abundance_sample_sample_3.png')
    assert os.path.exists(polar_plot)
    heatmap_occurrence = os.path.join(output_folder, 'function_occurrence', 'heatmap_occurrence.png')
    assert os.path.exists(heatmap_occurrence)

    shutil.rmtree(output_folder)


def test_search_hmm_custom_db_abundance_cli():
    input_file = os.path.join('input_data', 'org_prot')
    output_folder = 'output_folder'
    custom_db = os.path.join('input_data', 'custom_db')
    custom_motif = os.path.join('input_data', 'motif.json')
    custom_motif_pair = os.path.join('input_data', 'motif_pair.json')
    abundance_file = os.path.join('input_data', 'proteome_tax_id_abundance.tsv')

    subprocess.call(['bigecyhmm_custom', '-i', input_file, '-d', custom_db, '-o', output_folder, '-m', custom_motif, '-p', custom_motif_pair])
    subprocess.call(['bigecyhmm_visualisation', 'genomes', '--bigecyhmm', output_folder, '-o', output_folder, '--abundance-file', abundance_file])

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

    subprocess.call(['bigecyhmm_custom', '-i', input_file, '-d', custom_db, '-o', output_folder, '-m', custom_motif, '-p', custom_motif_pair])
    subprocess.call(['bigecyhmm_visualisation', 'genomes', '--bigecyhmm', output_folder, '-o', output_folder, '--measure-file', measure_file])

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

    subprocess.call(['bigecyhmm_custom', '-i', input_file, '-d', custom_db, '-o', output_folder, '-m', custom_motif, '-p', custom_motif_pair])
    subprocess.call(['bigecyhmm_visualisation', 'genomes', '--bigecyhmm', output_folder, '-o', output_folder, '--abundance-file', abundance_file, '--measure-file', measure_file])

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

    subprocess.call(['bigecyhmm_custom', '-i', input_file, '-d', custom_db, '-o', output_folder, '-m', custom_motif, '-p', custom_motif_pair])
    subprocess.call(['bigecyhmm_visualisation', 'esmecata', '--esmecata', esmecata_folder, '--bigecyhmm', output_folder, '-o', output_folder, '--abundance-file', abundance_file, '--measure-file', measure_file])

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


def test_search_hmm_custom_db_abundance_measure_esmecata_cli_modified_motif():
    input_file = os.path.join('input_data', 'esmecata_output_folder', '1_clustering', 'reference_proteins_consensus_fasta')
    output_folder = 'output_folder'
    custom_db = os.path.join('input_data', 'mini_custom_db')
    custom_motif = os.path.join('input_data', 'motif.json')
    custom_motif_pair = os.path.join('input_data', 'motif_pair_mod.json')
    abundance_file = os.path.join('input_data', 'proteome_tax_id_abundance.tsv')
    measure_file = os.path.join('input_data', 'test_measure.tsv')
    esmecata_folder = os.path.join('input_data', 'esmecata_output_folder')

    subprocess.call(['bigecyhmm_custom', '-i', input_file, '-d', custom_db, '-o', output_folder, '-m', custom_motif, '-p', custom_motif_pair])
    subprocess.call(['bigecyhmm_visualisation', 'esmecata', '--esmecata', esmecata_folder, '--bigecyhmm', output_folder, '-o', output_folder, '--abundance-file', abundance_file, '--measure-file', measure_file])

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


def test_bigecyhmm_custom_one_file_abundance():
    input_file = os.path.join('input_data', 'org_prot')
    output_folder = 'output_folder'
    custom_db = os.path.join('input_data', 'custom_db_one_file')
    custom_motif = os.path.join('input_data', 'motif.json')
    custom_motif_pair = os.path.join('input_data', 'motif_pair.json')
    abundance_file = os.path.join('input_data', 'proteome_tax_id_abundance.tsv')

    subprocess.call(['bigecyhmm_custom', '-i', input_file, '-d', custom_db, '-o', output_folder, '-m', custom_motif, '-p', custom_motif_pair])
    subprocess.call(['bigecyhmm_visualisation', 'genomes', '--bigecyhmm', output_folder, '-o', output_folder, '--abundance-file', abundance_file])

    expected_abundance = {'Acetate metabo.dissim.': {'sample_1': 100.0, 'sample_2': 0.0, 'sample_3': 0.0}, 'Acetogen. WL': {'sample_1': 100.0, 'sample_2': 0.0, 'sample_3': 0.0},
     'Carbon fixation': {'sample_1': 200.0, 'sample_2': 800.0, 'sample_3': 520.0}, 'Ethanol oxidation': {'sample_1': 100.0, 'sample_2': 0.0, 'sample_3': 0.0},
     'Fermentation': {'sample_1': 100.0, 'sample_2': 0.0, 'sample_3': 0.0}, 'Hydrogen generation': {'sample_1': 0.0, 'sample_2': 0.0, 'sample_3': 0.0},
     'Hydrogen oxidation': {'sample_1': 0.0, 'sample_2': 0.0, 'sample_3': 0.0}, 'Methanogen.': {'sample_1': 0.0, 'sample_2': 0.0, 'sample_3': 0.0},
     'Methanotrophy': {'sample_1': 100.0, 'sample_2': 0.0, 'sample_3': 0.0}, 'Organic carbon oxidation': {'sample_1': 100.0, 'sample_2': 0.0, 'sample_3': 0.0}}

    output_abundance_network_file = os.path.join(output_folder, 'carbon_cycle_od', 'cycle_diagram_bipartite_abundance.graphml')
    abundance_network = nx.read_graphml(output_abundance_network_file)

    predicted_abundance = {node: abundance_network.nodes[node] for node in abundance_network.nodes}

    for function in expected_abundance:
        for sample in expected_abundance[function]:
            assert expected_abundance[function][sample] == predicted_abundance[function][sample]

    shutil.rmtree(output_folder)


if __name__ == "__main__":
    test_create_visualisation_abundance_from_esmecata_cli()