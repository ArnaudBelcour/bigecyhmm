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

import csv
import os
import zipfile
import logging
import pyhmmer
import time
import sys
import json

from multiprocessing import Pool
from PIL import __version__ as pillow_version

from bigecyhmm.utils import is_valid_dir, file_or_folder, parse_result_files
from bigecyhmm.diagram_cycles import create_input_diagram, create_diagram_figures
from bigecyhmm import __version__ as bigecyhmm_version

ROOT = os.path.dirname(__file__)
HMM_COMPRESS_FILE = os.path.join(ROOT, 'hmm_databases', 'hmm_files.zip')
HMM_TEMPLATE_FILE = os.path.join(ROOT, 'hmm_databases', 'hmm_table_template.tsv')
PATHWAY_TEMPLATE_FILE = os.path.join(ROOT, 'hmm_databases', 'cycle_pathways.tsv')
PHENOTYPE_TEMPLATE_FILE = os.path.join(ROOT, 'hmm_databases', 'phenotypes.tsv')

logger = logging.getLogger(__name__)


def get_hmm_thresholds(hmm_template_file):
    """Extract threhsolds from HMM template file.

    Args:
        hmm_template_file (str): path of HMM template file

    Returns:
        hmm_thresholds (dict): threshold string for each HMM
    """
    with open(hmm_template_file, 'r') as open_hmm_template:
        csvreader = csv.DictReader(open_hmm_template, delimiter='\t')

        hmm_thresholds = {}
        for line in csvreader:
            for hmm_file in line['Hmm file'].split(', '):
                hmm_thresholds[hmm_file] = line['Hmm detecting threshold']

    return hmm_thresholds


def query_fasta_file(input_protein_fasta, hmm_thresholds):
    """Run HMM search with pyhmmer on protein fasta file

    Args:
        input_protein_fasta (str): path of protein fasta file
        hmm_thresholds (dict): threshold for each HMM

    Returns:
        results (list): list of result for HMM search, which are sublist containing: evalue, score and length
    """
    input_filename = os.path.splitext(os.path.basename(input_protein_fasta))[0]

    # Extract the sequence from the protein fasta files.
    with pyhmmer.easel.SequenceFile(input_protein_fasta, digital=True) as seq_file:
        sequences = list(seq_file)

    # Iterate on the HMM to query them. 
    results = []
    with zipfile.ZipFile(HMM_COMPRESS_FILE, 'r') as zip_object:
        for hmm_filename in zip_object.namelist():
            if hmm_filename.endswith('.hmm') and 'check' not in hmm_filename:
                hmm_filebasename = os.path.basename(hmm_filename)
                with zip_object.open(hmm_filename) as open_hmm_zipfile:
                    with pyhmmer.plan7.HMMFile(open_hmm_zipfile) as hmm_file:
                        for threshold_data in hmm_thresholds[hmm_filebasename].split(', '):
                            threshold, threshold_type = threshold_data.split('|')
                            threshold = float(threshold)
                            if threshold_type == 'full':
                                for hits in pyhmmer.hmmsearch(hmm_file, sequences, T=threshold, cpus=1):
                                    for hit in hits:
                                            results.append([input_filename, hit.name.decode(), hmm_filebasename, hit.evalue, hit.score, hit.length])
                            if threshold_type == 'domain':
                                for hits in pyhmmer.hmmsearch(hmm_file, sequences, domT=threshold, cpus=1):
                                    for hit in hits:
                                            results.append([input_filename, hit.name.decode(), hmm_filebasename, hit.evalue, hit.score, hit.length])

    return results


def write_results(hmm_results, output_file):
    """Write HMM results in a tsv file 

    Args:
        hmm_results (list): list of result for HMM search, which are sublist containing: evalue, score and length
        output_file (str): path to ouput tsv file
    """
    with open(output_file, 'w') as open_output_file:
        csvwriter = csv.writer(open_output_file, delimiter='\t')
        csvwriter.writerow(['organism', 'protein', 'HMM', 'evalue', 'score', 'length'])
        for result in hmm_results:
            csvwriter.writerow(result)


def create_major_functions(hmm_output_folder, output_file):
    """Map hit HMMs with list of major functions to create a tsv file showing these results.

    Args:
        hmm_output_folder (str): path to HMM search results folder (one tsv file per organism)
        output_file (str): path to the output tsv file
    """
    with open(HMM_TEMPLATE_FILE, 'r') as open_hmm_template:
        csvreader = csv.DictReader(open_hmm_template, delimiter='\t')

        hmm_functions = {}
        for line in csvreader:
            for hmm_file in line['Hmm file'].split(', '):
                function_name = line['Function'] + ' ' + line['Gene abbreviation']
                if function_name not in hmm_functions:
                    hmm_functions[function_name] = [hmm_file]
                else:
                    hmm_functions[function_name].append(hmm_file)

    hmm_list_functions = [function for function in hmm_functions]
    hmm_hits = parse_result_files(hmm_output_folder)
    org_list = [org for org in hmm_hits]
    with open(output_file, 'w') as open_output_file:
        csvwriter = csv.writer(open_output_file, delimiter='\t')
        csvwriter.writerow(['function', *org_list])
        for function in hmm_list_functions:
            present_functions = [len(set(hmm_functions[function]).intersection(set(hmm_hits[org])))/len(set(hmm_functions[function])) if len(set(hmm_functions[function]).intersection(set(hmm_hits[org]))) > 0 else 'NA' for org in org_list]
            csvwriter.writerow([function, *present_functions])


def create_phenotypes(hmm_output_folder, output_file):
    """Map hit HMMs with list of phenotypes to create a tsv file showing these results.

    Args:
        hmm_output_folder (str): path to HMM search results folder (one tsv file per organism)
        output_file (str): path to the output tsv file
    """
    with open(PHENOTYPE_TEMPLATE_FILE, 'r') as open_hmm_template:
        csvreader = csv.DictReader(open_hmm_template, delimiter='\t')

        hmm_functions = {}
        for line in csvreader:
            for hmm_file in line['HMMs'].split(', '):
                function_name = line['Phenotypes']
                if function_name not in hmm_functions:
                    hmm_functions[function_name] = [hmm_file]
                else:
                    hmm_functions[function_name].append(hmm_file)

    hmm_list_functions = [function for function in hmm_functions]
    hmm_hits = parse_result_files(hmm_output_folder)
    org_list = [org for org in hmm_hits]
    with open(output_file, 'w') as open_output_file:
        csvwriter = csv.writer(open_output_file, delimiter='\t')
        csvwriter.writerow(['function', *org_list])
        for function in hmm_list_functions:
            present_functions = [len(set(hmm_functions[function]).intersection(set(hmm_hits[org])))/len(set(hmm_functions[function])) if len(set(hmm_functions[function]).intersection(set(hmm_hits[org]))) > 0 else 'NA' for org in org_list]
            csvwriter.writerow([function, *present_functions])


def hmm_search_write_results(input_file_path, output_file, hmm_thresholds):
    """Little functions for the starmap multiprocessing to launch HMM search and result writing

    Args:
        input_file_path (str): path of protein fasta file
        output_file (str): output tsv file containing HMM search hits
        hmm_thresholds (dict): threshold for each HMM
    """
    logger.info('Search for HMMs on ' + input_file_path)
    hmm_results = query_fasta_file(input_file_path, hmm_thresholds)
    write_results(hmm_results, output_file)


def search_hmm(input_variable, output_folder, core_number=1):
    """Main function to use HMM search on protein sequences and write results

    Args:
        input_variable (str): path to input file or folder
        output_folder (str): path to output folder
        core_number (int): number of core to use for the multiprocessing
    """
    start_time = time.time()
    input_dicts = file_or_folder(input_variable)

    hmm_output_folder = os.path.join(output_folder, 'hmm_results')
    is_valid_dir(hmm_output_folder)

    hmm_thresholds = get_hmm_thresholds(HMM_TEMPLATE_FILE)

    hmm_search_pool = Pool(processes=core_number)

    multiprocess_input_hmm_searches = []
    for input_file in input_dicts:
        input_filename = os.path.splitext(os.path.basename(input_file))[0]
        output_file = os.path.join(hmm_output_folder, input_filename + '.tsv')

        input_file_path = input_dicts[input_file]
        multiprocess_input_hmm_searches.append([input_file_path, output_file, hmm_thresholds])

    hmm_search_pool.starmap(hmm_search_write_results, multiprocess_input_hmm_searches)

    hmm_search_pool.close()
    hmm_search_pool.join()

    function_matrix_file = os.path.join(output_folder, 'function_presence.tsv')
    create_major_functions(hmm_output_folder, function_matrix_file)
    function_matrix_file = os.path.join(output_folder, 'phenotypes_presence.tsv')
    create_phenotypes(hmm_output_folder, function_matrix_file)

    input_diagram_folder = os.path.join(output_folder, 'diagram_input')
    create_input_diagram(hmm_output_folder, input_diagram_folder, output_folder)

    input_diagram_file = os.path.join(output_folder, 'Total.R_input.txt')
    create_diagram_figures(input_diagram_file, output_folder)

    duration = time.time() - start_time
    metadata_json = {}
    metadata_json['tool_dependencies'] = {}
    metadata_json['tool_dependencies']['python_package'] = {}
    metadata_json['tool_dependencies']['python_package']['Python_version'] = sys.version
    metadata_json['tool_dependencies']['python_package']['bigecyhmm'] = bigecyhmm_version
    metadata_json['tool_dependencies']['python_package']['pyhmmer'] = pyhmmer.__version__
    metadata_json['tool_dependencies']['python_package']['pillow'] = pillow_version

    metadata_json['input_parameters'] = {'input_variable': input_variable, 'output_folder': output_folder, 'core_number': core_number}
    metadata_json['duration'] = duration

    metadata_file = os.path.join(output_folder, 'bigecyhmm_metadata.json')
    with open(metadata_file, 'w') as ouput_file:
        json.dump(metadata_json, ouput_file, indent=4)
