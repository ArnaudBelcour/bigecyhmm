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
import pyhmmer
import zipfile

from bigecyhmm.utils import is_valid_dir, file_or_folder
from bigecyhmm.diagram_cycles import create_input_diagram
from bigecyhmm.hmm_search import get_hmm_thresholds, hmm_search_write_results, create_major_functions
from bigecyhmm.utils import read_measures_file, read_esmecata_proteome_file
from bigecyhmm import __version__ as bigecyhmm_version
from bigecyhmm import HMM_COMPRESSED_FILE, HMM_TEMPLATE_FILE, MOTIF, MOTIF_PAIR, CUSTOM_CARBON_CYCLE_NETWORK, \
                    CUSTOM_SULFUR_CYCLE_NETWORK, CUSTOM_NITROGEN_CYCLE_NETWORK, CUSTOM_PHOSPHORUS_CYCLE_NETWORK, \
                    CUSTOM_HYDROGENOTROPHIC_CYCLE_NETWORK, CUSTOM_OTHER_CYCLE_NETWORK

from multiprocessing import Pool

MESSAGE = '''
Run bigecyhmm using a custom database (custom biogeochemical cycles with HMMs).
'''
REQUIRES = '''
Requires pyhmmer, networkx, matplotlib.
'''

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

try:
    import networkx as nx
    from networkx.readwrite import json_graph
except:
    logger.critical('networkx not installed, bigecyhmm_custom requires networkx installed: pip install networkx matplotlib')
    sys.exit(1)

try:
    import matplotlib.pyplot as plt
    from matplotlib import __version__ as matplotlib_version
except:
    logger.critical('matplotlib not installed, bigecyhmm_custom requires matplotlib installed: pip install networkx matplotlib')
    sys.exit(1)


def search_hmm_custom_db(input_variable, custom_database_json, output_folder, hmm_compressed_database=HMM_COMPRESSED_FILE,
                         hmm_template_file=HMM_TEMPLATE_FILE, core_number=1, motif_json=None, motif_pair_json=None,
                         abundance_file=None, metabolite_measure=None, esmecata_output_folder=None):
    """Main function to use HMM search on protein sequences and write results with a custom database.

    Args:
        input_variable (str): path to input file or folder
        custom_database_json (str): path to json file containing the custom database
        output_folder (str): path to output folder
        hmm_compressed_database (str): path to HMM compress database
        hmm_template_file (str): path to HMM template file
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

    # Extract thresholds from HMM template file.
    hmm_thresholds = get_hmm_thresholds(hmm_template_file)

    # Get pathway cycle data from custom database.
    logger.info("  -> Parsing cycle json file {0}".format(custom_database_json))
    with open(custom_database_json) as open_custom_database_json:
        json_cycle_database = json.load(open_custom_database_json)

    pathway_template_file = os.path.join(output_folder, 'pathway_template_file.tsv')
    already_search_function = []
    hmms_in_pathway_template = []
    with open(pathway_template_file, 'w') as open_pathway_template_file:
        csvwriter = csv.writer(open_pathway_template_file, delimiter='\t')
        csvwriter.writerow(['Pathways', 'HMMs'])
        for node in json_cycle_database['nodes']:
            if node['type'] == 'Function':
                function_name = node['id']
                if function_name not in already_search_function:
                    function_hmms = node['hmm']
                    csvwriter.writerow([function_name, function_hmms])
                    already_search_function.append(function_name)
                    if function_hmms != '':
                        hmms_in_pathway_template.extend([hmm.replace('NO|', '') for hmms in function_hmms.split(', ') for hmm in hmms.split('; ')])

    # Check that the same HMM profiles are present in compressed zip
    with zipfile.ZipFile(hmm_compressed_database, 'r') as zip_object:
        list_of_hmms = set([os.path.basename(hmm_filename) for hmm_filename in zip_object.namelist() if hmm_filename.endswith('.hmm') and 'check' not in hmm_filename])
    # Check that the HMMs in the pathway template file are present in the compressed database.
    hmms_in_pathway_template = set(hmms_in_pathway_template)
    if not hmms_in_pathway_template.issubset(list_of_hmms):
        not_found_hmms = hmms_in_pathway_template - list_of_hmms
        logger.critical("  Some HMMs present in {0} are not present in the HMM compressed database {1}: {2}".format(custom_database_json, hmm_compressed_database, not_found_hmms))
        sys.exit(1)
    # Check that the HMMs in the threshold file are present in the compressed database.
    hmm_in_threshold_file = set(hmm_thresholds.keys())
    if not hmm_in_threshold_file.issubset(list_of_hmms):
        not_found_hmms = hmm_in_threshold_file - list_of_hmms
        logger.critical("  Some HMMs present in {0} are not present in the HMM compressed database {1}: {2}".format(hmm_template_file, hmm_compressed_database, not_found_hmms))
        sys.exit(1)
    # Check that the HMMs in the pathway template file are present in the threshold file.
    if not hmms_in_pathway_template.issubset(hmm_in_threshold_file):
        not_found_hmms = hmms_in_pathway_template - hmm_in_threshold_file
        logger.critical("  Some HMMs present in {0} are not present in the HMM template file {1}: {2}".format(custom_database_json, hmm_template_file, not_found_hmms))
        sys.exit(1)

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
        logger.info("  -> Read EsMeCaTa proteome_tax_id file.")
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
    for input_filename in input_dicts:
        output_file = os.path.join(hmm_output_folder, input_filename + '.tsv')
        input_file_path = input_dicts[input_filename]
        multiprocess_input_hmm_searches.append([input_file_path, output_file, hmm_thresholds, hmm_compressed_database, motif_data, motif_pair_data])

    hmm_search_pool.starmap(hmm_search_write_results, multiprocess_input_hmm_searches)

    hmm_search_pool.close()
    hmm_search_pool.join()

    logger.info("  -> Create output files.")
    function_matrix_file = os.path.join(output_folder, 'function_presence.tsv')
    create_major_functions(hmm_output_folder, function_matrix_file, hmm_template_file)

    input_diagram_folder = os.path.join(output_folder, 'diagram_input')
    create_input_diagram(hmm_output_folder, input_diagram_folder, output_folder, pathway_template_file)

    bipartite_cycle_network = json_graph.node_link_graph(json_cycle_database, edges='edges')

    pathway_data = {}
    total_file = os.path.join(output_folder, 'Total.R_input.txt')
    with open(total_file, 'r') as open_total_file:
        csvreader = csv.reader(open_total_file, delimiter='\t')
        for line in csvreader:
            pathway = line[0]
            nb_pathway = line[1]
            percentage_pathway = line[2]
            pathway_data[pathway] = (nb_pathway, percentage_pathway)

    # Initiate a bipartite network when giving abundance or measure files.
    if abundance_file is not None or metabolite_measure is not None:
        abundance_cycle_network = nx.DiGraph()
        nodes_data = [(node, bipartite_cycle_network.nodes[node]) for node in bipartite_cycle_network.nodes]
        abundance_cycle_network.add_nodes_from(nodes_data)
        for edge in bipartite_cycle_network.edges:
            abundance_cycle_network.add_edges_from([edge], **bipartite_cycle_network.edges[edge])

    # Read abundance data and create a network with nodes associated with the abundance.
    if abundance_file is not None:
        sample_abundance, sample_tot_abundance = read_measures_file(abundance_file)

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
                    abundance_cycle_network.add_node(function_name, type='Function')
                for sample in sample_abundance:
                    # If esmecata output are given, translate tax_id into organism names.
                    if esmecata_output_folder is not None:
                        pathway_abundance[function_name][sample] = 0
                        for tax_id in organisms:
                            pathway_abundance[function_name][sample] += sum([sample_abundance[sample][organism] for organism in tax_id_names_observation_names[tax_id] if float(line[tax_id]) > 0 and organism in sample_abundance[sample]])
                    else:
                        pathway_abundance[function_name][sample] = sum([sample_abundance[sample][organism] for organism in organisms if float(line[organism]) > 0 and organism in sample_abundance[sample]])
                    abundance_cycle_network.nodes[function_name][sample] = pathway_abundance[function_name][sample]

        # Create a tsv file containing abundance information for cycle function.
        pathway_abundance_data = []
        all_cycles = list(pathway_abundance.keys())
        all_samples = list(sample_abundance.keys())
        for cycle in all_cycles:
            pathway_abundance_data.append([cycle, *[pathway_abundance[cycle][sample] for sample in all_samples]])

        cycle_pathway_abundance_file = os.path.join(output_folder, 'pathway_abundance.tsv')
        with open(cycle_pathway_abundance_file, 'w') as open_cycle_pathway_abundance_file:
            csvwriter = csv.writer(open_cycle_pathway_abundance_file, delimiter='\t')
            csvwriter.writerow(['function', *all_samples])
            for data in pathway_abundance_data:
                csvwriter.writerow(data)

    # Represent network as a bipartie graph with occurrence if functions.
    bipartite_graph = nx.DiGraph()

    mapped_new_node_ids = {}
    all_pathways = []
    for node in json_cycle_database['nodes']:
        if node['type'] == 'Function':
            function_name = node['id']
            function_name_weighted = function_name + '\nCoverage: ' + pathway_data[function_name][1]
            node_data = {node_d: node[node_d] for node_d in node if node_d != 'id'}
            bipartite_graph.add_node(function_name_weighted, **node_data)
            mapped_new_node_ids[function_name] = function_name_weighted
            all_pathways.append(function_name_weighted)
        else:
            bipartite_graph.add_node(node['id'], **node)

    for edge in json_cycle_database['edges']:
        source = edge['source']
        target = edge['target']
        if source in mapped_new_node_ids:
            function_name = source
            source = mapped_new_node_ids[source]
        if target in mapped_new_node_ids:
            function_name = target
            target = mapped_new_node_ids[target]
        bipartite_graph.add_edge(source, target, weight=pathway_data[function_name][1])

    logger.info("  -> Generate network files.")

    network_graphml_output_file = os.path.join(output_folder, 'cycle_diagram_bipartite_occurrence.graphml')
    nx.write_graphml(bipartite_graph, network_graphml_output_file)

    fig, axes = plt.subplots(figsize=(40,20))
    pos = nx.spring_layout(bipartite_graph)
    # Get node names:
    labels = {node: node for node in bipartite_graph.nodes}
    nx.draw_networkx_labels(bipartite_graph, pos, labels, font_size=20)
    nx.draw(bipartite_graph, pos, node_size=1000, arrowsize=20)
    network_output_file = os.path.join(output_folder, 'cycle_diagram_bipartite_occurrence.png')
    plt.savefig(network_output_file)

    if metabolite_measure is not None:
        # Add sample metabolite measurements to node.
        if metabolite_measure is not None:
            sample_metabolite_measure, sample_tot_abundance = read_measures_file(metabolite_measure)
            for metabolite in abundance_cycle_network.nodes:
                for sample in sample_metabolite_measure:
                    if metabolite in sample_metabolite_measure[sample]:
                        metabolite_measure = sample_metabolite_measure[sample][metabolite]
                        if isinstance(metabolite_measure, float):
                            abundance_cycle_network.nodes[metabolite][sample] = metabolite_measure

    if abundance_file is not None or metabolite_measure is not None:
        abundance_cycle_network.add_edges_from(bipartite_cycle_network.edges)
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

    metadata_json['input_parameters'] = {'input_variable': input_variable, 'output_folder': output_folder, 'core_number': core_number,
                                         'motif_json': motif_json, 'motif_pair_json': motif_pair_json, 'abundance_file': abundance_file,
                                         'metabolite_measure': metabolite_measure, 'esmecata_output_folder': esmecata_output_folder}
    metadata_json['input_parameters']['custom_db'] = {'hmm_compressed_database': hmm_compressed_database, 'hmm_template_file': hmm_template_file,
                                                      'json_database_file_path': custom_database_json}

    metadata_json['duration'] = duration

    metadata_file = os.path.join(output_folder, 'bigecyhmm_custom_metadata.json')
    with open(metadata_file, 'w') as ouput_file:
        json.dump(metadata_json, ouput_file, indent=4)


def identify_run_custom_db_search(input_variable, custom_database_folder, output_folder, core_number=1, motif_json=None, motif_pair_json=None,
                         abundance_file=None, metabolite_measure=None, esmecata_output_folder=None):
    """Main function to use HMM search on protein sequences and write results with a custom database.

    Args:
        input_variable (str): path to input file or folder
        custom_database_folder (str): path to file/folder containing custom database
        output_folder (str): path to output folder
        core_number (int): number of core to use for the multiprocessing
        motif_json (str): JSON file containing gene associated with protein motifs to check for predictions
        motif_pair_json (str): JSON file containing association between two genes to check for predictions
        abundance_file_path (str): path to abundance file indicating the abundance of organisms in samples
        metabolite_measure (str): path to metabolite measure file indicating the abundance of metabolites in samples
        esmecata_output_folder (str): path to esmecata output folder
    """
    start_time = time.time()

    # Search for json files in input custom database.
    if custom_database_folder == 'internal_carbon':
        custom_database_folder = CUSTOM_CARBON_CYCLE_NETWORK
    if custom_database_folder == 'internal_sulfur':
        custom_database_folder = CUSTOM_SULFUR_CYCLE_NETWORK
    if custom_database_folder == 'internal_nitrogen':
        custom_database_folder = CUSTOM_NITROGEN_CYCLE_NETWORK
    if custom_database_folder == 'internal_phosphorus':
        custom_database_folder = CUSTOM_PHOSPHORUS_CYCLE_NETWORK
    if custom_database_folder == 'internal_other':
        custom_database_folder = CUSTOM_OTHER_CYCLE_NETWORK
    if custom_database_folder == 'internal_hydrogenotrophic':
        custom_database_folder = CUSTOM_HYDROGENOTROPHIC_CYCLE_NETWORK
    if custom_database_folder == 'internal_all':
        all_json_network = [CUSTOM_CARBON_CYCLE_NETWORK, CUSTOM_SULFUR_CYCLE_NETWORK, CUSTOM_NITROGEN_CYCLE_NETWORK, CUSTOM_PHOSPHORUS_CYCLE_NETWORK, CUSTOM_OTHER_CYCLE_NETWORK]
        all_json_dict = {'directed': True, 'multigraph': False, 'graph': {'node_default': {}, 'edge_default': {}}, 'nodes': [], 'edges': []}
        internal_all_json_file = os.path.join(output_folder, 'internal_all.json')
        for json_network_file in all_json_network:
            with open(json_network_file) as open_json_network_file:
                json_cycle_database = json.load(open_json_network_file)
                all_json_dict['nodes'].extend(json_cycle_database['nodes'])
                all_json_dict['edges'].extend(json_cycle_database['edges'])
        with open(internal_all_json_file, 'w') as ouput_file:
            json.dump(all_json_dict, ouput_file, indent=4)
        custom_database_folder = internal_all_json_file

    json_extensions = ['.json']
    input_dicts = file_or_folder(custom_database_folder, json_extensions)
    for input_filename in input_dicts:
        logger.info("Launch HMM search on custom database {0}.".format(input_filename))
        custom_database_json = input_dicts[input_filename]
        # Check the presence of custom HMM compressed database.
        custom_hmm_compressed_database = custom_database_json.replace('.json', '.zip')
        if os.path.exists(custom_hmm_compressed_database):
            logger.info("  -> Custom HMM compressed database {0} exists, it will be used.".format(custom_hmm_compressed_database))
        else:
            logger.info("  -> Custom HMM compressed database {0} not found, bigecyhmm will use the default one.".format(custom_hmm_compressed_database))
            custom_hmm_compressed_database = HMM_COMPRESSED_FILE
        # Check the presence of custom HMM template file.
        custom_hmm_template_file = custom_database_json.replace('.json', '.tsv')
        if os.path.exists(custom_hmm_template_file):
            logger.info("  -> Custom HMM template file {0} exists, it will be used.".format(custom_hmm_template_file))
        else:
            logger.info("  -> Custom HMM template file {0} not found, bigecyhmm will use the default one.".format(custom_hmm_template_file))
            custom_hmm_template_file = HMM_TEMPLATE_FILE

        # Create specific output folder for custom database.
        output_folder_custom_db = os.path.join(output_folder, input_filename)
        is_valid_dir(output_folder_custom_db)
        search_hmm_custom_db(input_variable, custom_database_json, output_folder_custom_db, custom_hmm_compressed_database,
                            custom_hmm_template_file, core_number=core_number, motif_json=motif_json, motif_pair_json=motif_pair_json,
                            abundance_file=abundance_file, metabolite_measure=metabolite_measure, esmecata_output_folder=esmecata_output_folder)

    duration = time.time() - start_time
    metadata_json = {}
    metadata_json['tool_dependencies'] = {}
    metadata_json['tool_dependencies']['python_package'] = {}
    metadata_json['tool_dependencies']['python_package']['Python_version'] = sys.version
    metadata_json['tool_dependencies']['python_package']['bigecyhmm'] = bigecyhmm_version
    metadata_json['tool_dependencies']['python_package']['pyhmmer'] = pyhmmer.__version__
    metadata_json['tool_dependencies']['python_package']['networkx'] = nx.__version__
    metadata_json['tool_dependencies']['python_package']['matplotlib'] = matplotlib_version

    metadata_json['input_parameters'] = {'input_variable': input_variable, 'custom_database_folder': custom_database_folder, 'output_folder': output_folder,
                                         'core_number': core_number, 'motif_json': motif_json, 'motif_pair_json': motif_pair_json, 'abundance_file': abundance_file,
                                         'metabolite_measure': metabolite_measure, 'esmecata_output_folder': esmecata_output_folder}

    metadata_json['duration'] = duration

    metadata_file = os.path.join(output_folder, 'bigecyhmm_custom_metadata.json')
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
        help='Path to a json file or folder containing a representation of the custom cycle. It will also search for associated tsv and zip file. If it is a folder, it will do the same but for each json in the folder.',
        metavar='CUSTOM_DATABASE_FILE_OR_FOLDER')

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
    identify_run_custom_db_search(args.input, args.custom_database, args.output, args.core, args.motif_file, args.motif_pair_file,
                         args.abundance_file, args.measure_file, args.esmecata_folder)

    duration = time.time() - start_time
    logger.info("--- Total runtime %.2f seconds ---" % (duration))
    logger.warning(f'--- Logs written in {log_file_path} ---')
