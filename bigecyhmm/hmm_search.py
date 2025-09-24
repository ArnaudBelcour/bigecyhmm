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
import re
import sys
import json

from multiprocessing import Pool
from PIL import __version__ as pillow_version

from bigecyhmm.utils import is_valid_dir, file_or_folder, parse_result_files
from bigecyhmm.diagram_cycles import create_input_diagram, create_diagram_figures
from bigecyhmm import __version__ as bigecyhmm_version
from bigecyhmm import HMM_COMPRESSED_FILE, HMM_TEMPLATE_FILE, MOTIF, MOTIF_PAIR

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


def check_motif_regex(gene_name, sequence, motif_db=MOTIF):
    """ Check the presence of a motif in a protein sequence using regex.

    Args:
        gene_name (str): gene name associated with the protein sequence
        sequence (str): string of the protein sequence
        motif_db (dict): dictionary containing gene name as key and motif to search as values

    Returns:
        boolean: True if motif found, False if not
    """
    motif_regex = motif_db[gene_name]
    # Replace X by any amino-acid.
    motif_regex_gene = re.sub(r'X', r'[ARNDCQEGHILKMFPSTWYV]', motif_regex)
    motif_found = re.findall(motif_regex_gene, sequence)
    if len(motif_found) > 0:
        return True
    else:
        False


def check_motif_pair(input_sequence, hmm_filename, pair_hmm_filename, zip_object):
    """ Check for a protein sequence and a HMM if it is not better associated with another HMM.

    Args:
        input_sequence (list): list of input sequences to check
        hmm_filename (str): path to the first HMM file
        pair_hmm_filename (str): path to the second HMM file
        zip_object (zipfile object): zip object associated with the compress HMM database

    Returns:
        boolean: True if first HMM has a better association with the sequence than the second HMM, False if not
    """
    with zip_object.open(hmm_filename) as open_hmm_zipfile:
        with pyhmmer.plan7.HMMFile(open_hmm_zipfile) as hmm_file:
            check_scores = [hit.score
                             for hits in pyhmmer.hmmsearch(hmm_file, input_sequence, cpus=1)
                             for hit in hits]
            if len(check_scores) > 0:
                motif_check_score = max(check_scores)
            else:
                motif_check_score = 0

    with zip_object.open(pair_hmm_filename) as open_hmm_zipfile:
        with pyhmmer.plan7.HMMFile(open_hmm_zipfile) as pair_hmm_file:
            anti_check_scores = [second_hit.score
                                    for second_hits in pyhmmer.hmmsearch(pair_hmm_file, input_sequence, cpus=1)
                                    for second_hit in second_hits]
            if len(anti_check_scores) > 0:
                motif_anti_check_score = max(anti_check_scores)
            else:
                motif_anti_check_score = 0

    if motif_check_score >= motif_anti_check_score and motif_check_score != 0:
        return True
    else:
        return False


def query_fasta_file(input_protein_fasta, hmm_thresholds, hmm_compressed_database=HMM_COMPRESSED_FILE, motif_db=MOTIF, motif_pair_db=MOTIF_PAIR, pyhmmer_core=1):
    """Run HMM search with pyhmmer on protein fasta file using HMM files from database.
    Use associated threshold either for full sequence or domain.

    Args:
        input_protein_fasta (str): path of protein fasta file
        hmm_thresholds (dict): threshold for each HMM
        hmm_compressed_database (str): path to HMM compress database
        motif_db (dict): dictionary containing gene name as key and motif to search as values
        motif_pair_db (dict): dictionary containing gene name as key and a second gene name as values
        pyhmmer_core (int): number of core used by pyhmmer

    Returns:
        results (list): list of result for HMM search, which are sublist containing: evalue, score and length
    """
    input_filename = os.path.splitext(os.path.basename(input_protein_fasta))[0]

    # Extract the sequence from the protein fasta files.
    with pyhmmer.easel.SequenceFile(input_protein_fasta, digital=True) as seq_file:
        sequences = pyhmmer.easel.DigitalSequenceBlock(pyhmmer.easel.Alphabet.amino(), seq_file)

    # Iterate on the HMM to query them. 
    results = []
    with zipfile.ZipFile(hmm_compressed_database, 'r') as zip_object:
        list_of_hmms = [hmm_filename for hmm_filename in zip_object.namelist() if hmm_filename.endswith('.hmm') and 'check' not in hmm_filename]
        check_hmms = {os.path.basename(hmm_filename).replace('.check.hmm', ''): hmm_filename
                      for hmm_filename in zip_object.namelist() if hmm_filename.endswith('.hmm') and 'check' in hmm_filename}
        for hmm_filename in list_of_hmms:
            hmm_filebasename = os.path.basename(hmm_filename)
            hmm_name = hmm_filebasename.replace('.hmm', '')
            with zip_object.open(hmm_filename) as open_hmm_zipfile:
                with pyhmmer.plan7.HMMFile(open_hmm_zipfile) as hmm_file:
                    for threshold_data in hmm_thresholds[hmm_filebasename].split(', '):
                        threshold, threshold_type = threshold_data.split('|')
                        threshold = float(threshold)
                        if threshold_type == 'full':
                            for hits in pyhmmer.hmmsearch(hmm_file, sequences, cpus=pyhmmer_core, Z=len(list_of_hmms), parallel="targets"):
                                for hit in hits.included:
                                    if hit.score >= threshold:
                                        gene_match = hit.name
                                        # Check the presence of specific motif in gene sequence.
                                        if hmm_name in motif_db:
                                            gene_sequence_str = [sequence for sequence in sequences if sequence.name == gene_match][0].textize().sequence
                                            if check_motif_regex(hmm_name, gene_sequence_str, motif_db):
                                                results.append([input_filename, gene_match.decode(), hmm_filebasename, hit.evalue, hit.score, hit.length])
                                        # Motif validation by checking that it is not better associated with another HMM.
                                        elif hmm_name in motif_pair_db:
                                            gene_sequence = [sequence for sequence in sequences if sequence.name == gene_match]
                                            first_check_hmm = check_hmms[hmm_name]
                                            second_check_hmm = check_hmms[motif_pair_db[hmm_name]]
                                            if check_motif_pair(gene_sequence, first_check_hmm, second_check_hmm, zip_object):
                                                results.append([input_filename, gene_match.decode(), hmm_filebasename, hit.evalue, hit.score, hit.length])
                                        else:
                                            results.append([input_filename, gene_match.decode(), hmm_filebasename, hit.evalue, hit.score, hit.length])
                        if threshold_type == 'domain':
                            for hits in pyhmmer.hmmsearch(hmm_file, sequences, cpus=pyhmmer_core, Z=len(list_of_hmms), parallel="targets"):
                                for hit in hits.included:
                                    for domain in hit.domains.included:
                                        if domain.score >= threshold:
                                            gene_match = hit.name
                                            if hmm_name in motif_db:
                                                gene_sequence_str = [sequence for sequence in sequences if sequence.name == gene_match][0]
                                                if check_motif_regex(hmm_name, gene_sequence_str):
                                                    results.append([input_filename, hit.name.decode(), hmm_filebasename, hit.evalue, domain.score, hit.length])
                                            elif hmm_name in motif_pair_db:
                                                gene_sequence = [sequence for sequence in sequences if sequence.name == gene_match]
                                                first_check_hmm = check_hmms[hmm_name]
                                                second_check_hmm = check_hmms[motif_pair_db[hmm_name]]
                                                if check_motif_pair(gene_sequence, first_check_hmm, second_check_hmm, zip_object):
                                                    results.append([input_filename, gene_match.decode(), hmm_filebasename, hit.evalue, domain.score, hit.length])
                                            else:
                                                results.append([input_filename, hit.name.decode(), hmm_filebasename, hit.evalue, domain.score, hit.length])

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


def create_major_functions(hmm_output_folder, output_file, hmm_template_file=HMM_TEMPLATE_FILE):
    """Map hit HMMs with list of major functions to create a tsv file showing these results.

    Args:
        hmm_output_folder (str): path to HMM search results folder (one tsv file per organism)
        output_file (str): path to the output tsv file
        hmm_template_file (str): path of HMM template file
    """
    with open(hmm_template_file, 'r') as open_hmm_template:
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


def hmm_search_write_results(input_file_path, output_file, hmm_thresholds, hmm_compressed_database=HMM_COMPRESSED_FILE, motif_db=MOTIF, motif_pair_db=MOTIF_PAIR, pyhmmer_core=1):
    """Little functions for the starmap multiprocessing to launch HMM search and result writing

    Args:
        input_file_path (str): path of protein fasta file
        output_file (str): output tsv file containing HMM search hits
        hmm_thresholds (dict): threshold for each HMM
        hmm_compressed_database (str): path to HMM compress database
        motif_db (dict): dictionary containing gene name as key and motif to search as values
        motif_pair_db (dict): dictionary containing gene name as key and a second gene name as values
        pyhmmer_core (int): number of core used by pyhmmer
    """
    logger.info('Search for HMMs on ' + input_file_path)
    hmm_results = query_fasta_file(input_file_path, hmm_thresholds, hmm_compressed_database, motif_db, motif_pair_db, pyhmmer_core)
    write_results(hmm_results, output_file)


def search_hmm(input_variable, output_folder, hmm_compressed_database=HMM_COMPRESSED_FILE, hmm_template_file=HMM_TEMPLATE_FILE, motif_db=MOTIF, motif_pair_db=MOTIF_PAIR, core_number=1):
    """Main function to use HMM search on protein sequences and write results

    Args:
        input_variable (str): path to input file or folder
        output_folder (str): path to output folder
        hmm_compressed_database (str): path to HMM compress database
        hmm_template_file (str): path of HMM template file
        motif_db (dict): dictionary containing gene name as key and motif to search as values
        motif_pair_db (dict): dictionary containing gene name as key and a second gene name as values
        core_number (int): number of core to use for the multiprocessing
    """
    start_time = time.time()
    input_dicts = file_or_folder(input_variable)

    logger.info('HMM compressed file: ' + hmm_compressed_database)
    logger.info('HMM template file : ' + hmm_template_file)

    # If there is only one input fasta file, use pyhmmer multiprocessing.
    if len(input_dicts) == 1:
        core_number = 1
        pyhmmer_core = 0
    # Otherwise allocate one CPU per fasta file running HMM search.
    else:
        pyhmmer_core = 1

    hmm_output_folder = os.path.join(output_folder, 'hmm_results')
    is_valid_dir(hmm_output_folder)

    hmm_thresholds = get_hmm_thresholds(hmm_template_file)

    hmm_search_pool = Pool(processes=core_number)

    multiprocess_input_hmm_searches = []
    for input_filename in input_dicts:
        output_file = os.path.join(hmm_output_folder, input_filename + '.tsv')
        input_file_path = input_dicts[input_filename]
        multiprocess_input_hmm_searches.append([input_file_path, output_file, hmm_thresholds, hmm_compressed_database, motif_db, motif_pair_db, pyhmmer_core])

    hmm_search_pool.starmap(hmm_search_write_results, multiprocess_input_hmm_searches)

    hmm_search_pool.close()
    hmm_search_pool.join()

    function_matrix_file = os.path.join(output_folder, 'function_presence.tsv')
    create_major_functions(hmm_output_folder, function_matrix_file)

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
