# Copyright (C) 2024-2025 Arnaud Belcour - Univ. Grenoble Alpes, Inria, Grenoble, France Microcosme
# Univ. Grenoble Alpes, Inria, Microcosme
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>

import argparse
import logging
import json
import os
import csv
import sys
import time
import networkx as nx
import matplotlib.pyplot as plt
import pyhmmer

from bigecyhmm.utils import is_valid_dir, file_or_folder, parse_result_files
from bigecyhmm.diagram_cycles import create_input_diagram, create_diagram_figures
from bigecyhmm.hmm_search import get_hmm_thresholds, hmm_search_write_results, create_major_functions, create_phenotypes
from bigecyhmm.utils import read_measures_file, read_esmecata_proteome_file
from bigecyhmm import __version__ as bigecyhmm_version
from bigecyhmm import HMM_COMPRESS_FILE, HMM_TEMPLATE_FILE, PHENOTYPE_TEMPLATE_FILE, MOTIF, MOTIF_PAIR

from multiprocessing import Pool
from networkx.readwrite import json_graph
from matplotlib import __version__ as matplotlib_version

MESSAGE = '''
Run bigecyhmm using a custom database (custom biogeochemical cycles with HMMs).
'''
REQUIRES = '''
Requires pyhmmer, networkx, matplotlib.
'''

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)


def search_hmm_custom_db(input_variable, custom_database_folder, output_folder, core_number=1, motif_json=None, motif_pair_json=None,
                         abundance_file=None, metabolite_measure=None, esmecata_output_folder=None):
    """Main function to use HMM search on protein sequences and write results with a custom database.

    Args:
        input_variable (str): path to input file or folder
        custom_database_folder (str): path to folder containing custom database
        output_folder (str): path to output folder
        core_number (int): number of core to use for the multiprocessing
        motif_json (str): JSON file containing gene associated with protein motifs to check for predictions
        motif_pair_json (str): JSON file containing association between two genes to check for predictions
        abundance_file_path (str): path to abundance file indicating the abundance of organisms in samples
        metabolite_measure (str): path to metaboltie measure file indicating the abundance of metabolites in samples
        esmecata_output_folder (str): path to esmecata output folder.
    """
    start_time = time.time()
    input_dicts = file_or_folder(input_variable)

    hmm_output_folder = os.path.join(output_folder, 'hmm_results')
    is_valid_dir(hmm_output_folder)

    # Get HMM database, either from internal database or custom database.
    hmm_compress_databases = [database_file for database_file in os.listdir(custom_database_folder) if database_file.endswith('.zip')]
    if len(hmm_compress_databases) > 0:
        hmm_compress_database = hmm_compress_databases[0]
        hmm_compress_database = os.path.join(custom_database_folder, hmm_compress_database)
        logger.info("Found custom HMM database {0}".format(hmm_compress_database))
    else:
        hmm_compress_database = HMM_COMPRESS_FILE
        logger.info("No custom HMM database, will use default: {0}".format(hmm_compress_database))

    # Get HMM threhsold either from internal database or custom database.
    hmm_template_files = [database_file for database_file in os.listdir(custom_database_folder) if database_file.endswith('.tsv')]
    if len(hmm_template_files) > 0:
        hmm_template_file = hmm_template_files[0]
        hmm_template_file = os.path.join(custom_database_folder, hmm_template_file)
        logger.info("Found custom HMM threshold file {0}".format(hmm_template_file))
    else:
        hmm_template_file = HMM_TEMPLATE_FILE
        logger.info("No custom HMM threshold file, will use default: {0}".format(hmm_template_file))

    hmm_thresholds = get_hmm_thresholds(hmm_template_file)

    # Get pathway cycle data from custom database.
    json_database_file = [database_file for database_file in os.listdir(custom_database_folder) if database_file.endswith('.json')][0]
    json_database_file_path = os.path.join(custom_database_folder, json_database_file)
    logger.info("Parsing cycle json file {0}".format(json_database_file_path))
    with open(json_database_file_path) as open_json_database_file_path:
        json_cycle_database = json.load(open_json_database_file_path)

    pathway_template_file = os.path.join(output_folder, 'pathway_template_file.tsv')
    already_search_function = []
    with open(pathway_template_file, 'w') as open_pathway_template_file:
        csvwriter = csv.writer(open_pathway_template_file, delimiter='\t')
        csvwriter.writerow(['Pathways', 'HMMs'])
        for cycle in json_cycle_database['edges']:
            function_name = cycle['id']
            if function_name not in already_search_function:
                function_hmms = cycle['hmm']
                csvwriter.writerow([function_name, function_hmms])
                already_search_function.append(function_name)

    # Get motif and motif_pair dictionaries.
    if motif_json is not None:
        with open(motif_json, 'r') as open_motif_json:
            motif_data = json.load(open_motif_json)
    else:
        motif_data = MOTIF

    if motif_pair_json is not None:
        with open(motif_pair_json, 'r') as open_motif_json:
            motif_pair_data = json.load(open_motif_json)
    else:
        motif_pair_data = MOTIF_PAIR

    # Get observation name and taxon names from esmecata.
    if esmecata_output_folder is not None:
        logger.info("Read EsMeCaTa proteome_tax_id file.")
        proteome_tax_id_file = os.path.join(esmecata_output_folder, '0_proteomes', 'proteome_tax_id.tsv')
        observation_names_tax_id_names, observation_names_tax_ranks = read_esmecata_proteome_file(proteome_tax_id_file)
        tax_id_names_observation_names = {}
        for observation_name in observation_names_tax_id_names:
            tax_id_name = observation_names_tax_id_names[observation_name]
            if tax_id_name not in tax_id_names_observation_names:
                tax_id_names_observation_names[tax_id_name] = [observation_name]
            else:
                tax_id_names_observation_names[tax_id_name].append(observation_name)

    hmm_search_pool = Pool(processes=core_number)

    multiprocess_input_hmm_searches = []
    for input_file in input_dicts:
        input_filename = os.path.splitext(os.path.basename(input_file))[0]
        output_file = os.path.join(hmm_output_folder, input_filename + '.tsv')

        input_file_path = input_dicts[input_file]
        multiprocess_input_hmm_searches.append([input_file_path, output_file, hmm_thresholds, hmm_compress_database, motif_data, motif_pair_data])

    hmm_search_pool.starmap(hmm_search_write_results, multiprocess_input_hmm_searches)

    hmm_search_pool.close()
    hmm_search_pool.join()

    logger.info("Create output files.")
    function_matrix_file = os.path.join(output_folder, 'function_presence.tsv')
    create_major_functions(hmm_output_folder, function_matrix_file)
    function_matrix_file = os.path.join(output_folder, 'phenotypes_presence.tsv')
    create_phenotypes(hmm_output_folder, function_matrix_file)

    input_diagram_folder = os.path.join(output_folder, 'diagram_input')
    create_input_diagram(hmm_output_folder, input_diagram_folder, output_folder, pathway_template_file)

    cycle_network = json_graph.node_link_graph(json_cycle_database, edges='edges')

    pathway_data = {}
    total_file = os.path.join(output_folder, 'Total.R_input.txt')
    with open(total_file, 'r') as open_total_file:
        csvreader = csv.reader(open_total_file, delimiter='\t')
        for line in csvreader:
            pathway = line[0]
            nb_pathway = line[1]
            percentage_pathway = line[2]
            pathway_data[pathway] = (nb_pathway, percentage_pathway)

    # Read abundance data and create a network with nodes associated with the abundance.
    if abundance_file is not None or metabolite_measure is not None:
        abundance_cycle_network = nx.DiGraph()

    if abundance_file is not None:
        sample_abundance, sample_tot_abundance = read_measures_file(abundance_file)
        bipartite_edge_abundances = []
        all_pathway_abundances = []

        pathway_abundance = {}
        cycle_pathway_presence_file = os.path.join(output_folder, 'pathway_presence.tsv')
        with open(cycle_pathway_presence_file, 'r') as open_cycle_pathway_presence_file:
            csvreader = csv.DictReader(open_cycle_pathway_presence_file, delimiter='\t')
            headers = csvreader.fieldnames
            organisms = headers[1:]
            for line in csvreader:
                function_name = line['function']
                if function_name not in pathway_abundance:
                    pathway_abundance[function_name] = {}
                if function_name not in abundance_cycle_network.nodes:
                    abundance_cycle_network.add_node(function_name, type='function')
                for sample in sample_abundance:
                    # If esmecata output are given, translate tax_id into organism names.
                    if esmecata_output_folder is not None:
                        pathway_abundance[function_name][sample] = 0
                        for tax_id in organisms:
                            pathway_abundance[function_name][sample] += sum([sample_abundance[sample][organism] for organism in tax_id_names_observation_names[tax_id] if float(line[tax_id]) > 0 and organism in sample_abundance[sample]])
                    else:
                        pathway_abundance[function_name][sample] = sum([sample_abundance[sample][organism] for organism in organisms if float(line[organism]) > 0 and organism in sample_abundance[sample]])
                    abundance_cycle_network.nodes[function_name][sample] = pathway_abundance[function_name][sample]

    pathway_names = {}
    bipartite_edges = []
    all_pathways = []
    for pathway in json_cycle_database['edges']:
        function_name = pathway['id']
        source = pathway['source']
        target = pathway['target']
        pathway_names[(source, target)] = function_name
        function_name_weighted = function_name + '\nCoverage: ' + pathway_data[function_name][1]
        bipartite_edges.append((source, function_name_weighted))
        bipartite_edges.append((function_name_weighted, target))
        all_pathways.append(function_name_weighted)
        if abundance_file is not None or metabolite_measure is not None:
            for sample in sample_abundance:
                bipartite_edge_abundances.append((source, function_name))
                bipartite_edge_abundances.append((function_name, target))
                all_pathway_abundances.append(function_name)
        cycle_network[source][target]['weight'] = pathway_data[function_name][1]

    logger.info("Generate network files.")
    # (1) Represent network as a graph.
    fig, axes = plt.subplots(figsize=(40,20))
    pos = nx.circular_layout(cycle_network)

    # Get edge weight and name:
    edge_labels=dict([((u,v,),pathway_names[u,v] + '\nCoverage:' + d['weight']) for u,v,d in cycle_network.edges(data=True)])
    nx.draw_networkx_edge_labels(cycle_network, pos, edge_labels=edge_labels, font_size=20)
    # Get node names:
    labels = {node: node for node in cycle_network.nodes}
    nx.draw_networkx_labels(cycle_network, pos, labels, font_size=20)
    nx.draw(cycle_network, pos, node_size=1000, arrowsize=20)
    network_output_file = os.path.join(output_folder, 'cycle_diagram.png')
    plt.savefig(network_output_file)

    network_json_output_file = os.path.join(output_folder, 'cycle_diagram.json')
    with open(network_json_output_file, 'w') as open_network_json_output_file:
        json.dump(json_graph.node_link_data(cycle_network, edges='edges'), open_network_json_output_file, indent=4)

    network_graphml_output_file = os.path.join(output_folder, 'cycle_diagram.graphml')
    nx.write_graphml(cycle_network, network_graphml_output_file)

    # (2) Represent network as a bipartie graph with occurrence if functions.
    bipartite_graph = nx.DiGraph()
    bipartite_graph.add_nodes_from(cycle_network.nodes, type='metabolite')
    bipartite_graph.add_nodes_from(all_pathways, type='function')
    bipartite_graph.add_edges_from(bipartite_edges)
    network_graphml_output_file = os.path.join(output_folder, 'cycle_diagram_bipartite_occurrence.graphml')
    nx.write_graphml(bipartite_graph, network_graphml_output_file)

    fig, axes = plt.subplots(figsize=(40,20))
    pos = nx.spring_layout(bipartite_graph)
    # Get node names:
    labels = {node: node for node in bipartite_graph.nodes}
    nx.draw_networkx_labels(bipartite_graph, pos, labels, font_size=20)
    nx.draw(bipartite_graph, pos, node_size=1000, arrowsize=20)
    network_output_file = os.path.join(output_folder, 'cycle_diagram_bipartite.png')
    plt.savefig(network_output_file)

    if abundance_file is not None or metabolite_measure is not None:
        abundance_cycle_network.add_nodes_from(cycle_network.nodes, type='metabolite')
        # Add sample metabolite measurements to node.
        if metabolite_measure is not None:
            sample_metabolite_measure, sample_tot_abundance = read_measures_file(metabolite_measure)
            for metabolite in abundance_cycle_network.nodes:
                for sample in sample_metabolite_measure:
                    if metabolite in sample_metabolite_measure[sample]:
                        metabolite_measure = sample_metabolite_measure[sample][metabolite]
                        if isinstance(metabolite_measure, float):
                            abundance_cycle_network.nodes[metabolite][sample] = metabolite_measure

        abundance_cycle_network.add_edges_from(bipartite_edge_abundances)
        network_graphml_output_file = os.path.join(output_folder, 'cycle_diagram_bipartite_abundance.graphml')
        nx.write_graphml(abundance_cycle_network, network_graphml_output_file)

    duration = time.time() - start_time
    metadata_json = {}
    metadata_json['tool_dependencies'] = {}
    metadata_json['tool_dependencies']['python_package'] = {}
    metadata_json['tool_dependencies']['python_package']['Python_version'] = sys.version
    metadata_json['tool_dependencies']['python_package']['bigecyhmm'] = bigecyhmm_version
    metadata_json['tool_dependencies']['python_package']['pyhmmer'] = pyhmmer.__version__
    metadata_json['tool_dependencies']['python_package']['networkx'] = nx.__version__
    metadata_json['tool_dependencies']['python_package']['matplotlib'] = matplotlib_version

    metadata_json['input_parameters'] = {'input_variable': input_variable, 'output_folder': output_folder, 'core_number': core_number}
    metadata_json['input_parameters']['custom_db'] = {'hmm_compress_database': hmm_compress_database, 'hmm_template_file': hmm_template_file,
                                                      'json_database_file_path': json_database_file_path}

    metadata_json['duration'] = duration

    metadata_file = os.path.join(output_folder, 'bigecyhmm_metadata.json')
    with open(metadata_file, 'w') as ouput_file:
        json.dump(metadata_json, ouput_file, indent=4)


def main():
    start_time = time.time()

    parser = argparse.ArgumentParser(
        'bigecyhmm_custom',
        description=MESSAGE + ' For specific help on each subcommand use: esmecata {cmd} --help',
        epilog=REQUIRES
    )
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s ' + bigecyhmm_version + '\n')

    parser.add_argument(
        '-i',
        '--input',
        dest='input',
        required=True,
        help='Input data, either a protein fasta file or a folder containing protein fasta files.',
        metavar='INPUT_FILE_OR_FOLDER')

    parser.add_argument(
        '-d',
        '--database',
        dest='custom_database',
        required=True,
        help='Custom database: a folder containing a json file for the cycle representation, a tabulated file for HMM threshold and a zip file containing HMM profiles',
        metavar='CUSTOM_DATABASE_FOLDER')


    parser.add_argument(
        '-o',
        '--output',
        dest='output',
        required=True,
        help='Output directory path.',
        metavar='OUPUT_DIR')

    parser.add_argument(
        "-c",
        "--core",
        help="Number of cores for multiprocessing",
        required=False,
        type=int,
        default=1)

    parser.add_argument(
        "-m",
        "--motif",
        dest='motif_file',
        help="JSON file containing gene associated with protein motifs to check for predictions.",
        required=False,
        default=None)

    parser.add_argument(
        "-p",
        "--motif-pair",
        dest='motif_pair_file',
        help="JSON file containing association between two genes to check for predictions.",
        required=False,
        default=None)

    parser.add_argument(
        '--abundance-file',
        dest='abundance_file',
        required=False,
        help='Abundance file indicating the abundance for each organisms in different samples.',
        metavar='INPUT_FILE',
        default=None)

    parser.add_argument(
        '--measure-file',
        dest='measure_file',
        required=False,
        help='Measure file indicating the abundance for each metabolties of the graph in different samples.',
        metavar='INPUT_FILE',
        default=None)

    parser.add_argument(
        '--esmecata',
        dest='esmecata_folder',
        required=False,
        help='EsMeCaTa output folder for the input file.',
        metavar='INPUT_FILE',
        default=None)

    args = parser.parse_args()

    # If no argument print the help.
    if len(sys.argv) == 1 or len(sys.argv) == 0:
        parser.print_help()
        sys.exit(1)

    is_valid_dir(args.output)

    # add logger in file
    formatter = logging.Formatter('%(message)s')
    log_file_path = os.path.join(args.output, f'bigecyhmm_custom.log')
    file_handler = logging.FileHandler(log_file_path, 'w+')
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    # set up the default console logger
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    logger.info("--- Launch HMM search on custom database ---")
    search_hmm_custom_db(args.input, args.custom_database, args.output, args.core, args.motif_file, args.motif_pair_file,
                         args.abundance_file, args.measure_file, args.esmecata_folder)

    duration = time.time() - start_time
    logger.info("--- Total runtime %.2f seconds ---" % (duration))
    logger.warning(f'--- Logs written in {log_file_path} ---')
