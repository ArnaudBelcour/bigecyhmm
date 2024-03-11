# Copyright (C) 2024 Arnaud Belcour - Univ. Grenoble Alpes, Inria, Grenoble, France Microcosme
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
import pyhmmer

from multiprocessing import Pool

from bigecyhmm.utils import is_valid_dir, file_or_folder
from bigecyhmm.diagram_cycles import create_input_diagram

ROOT = os.path.dirname(__file__)
HMM_COMPRESS_FILE = os.path.join(ROOT, 'hmm_databases', 'hmm_files.zip')
HMM_TEMPLATE_FILE = os.path.join(ROOT, 'hmm_databases', 'hmm_table_template.tsv')
DIAGRAM_TEMPLATE_FILE = os.path.join(ROOT, 'hmm_databases', 'R_diagram_pathways.tsv')


def get_hmm_thresholds(hmm_template_file):
    with open(hmm_template_file, 'r') as open_hmm_template:
        csvreader = csv.DictReader(open_hmm_template, delimiter='\t')

        hmm_thresholds = {}
        for line in csvreader:
            for hmm_file in line['Hmm file'].split(', '):
                hmm_thresholds[hmm_file] = line['Hmm detecting threshold'].split('|')[0]

    return hmm_thresholds

def query_fasta_file(input_protein_fasta):
    input_filename = os.path.splitext(os.path.basename(input_protein_fasta))[0]

    # Extract the sequence from the protein fasta files.
    with pyhmmer.easel.SequenceFile(input_protein_fasta, digital=True) as seq_file:
        sequences = list(seq_file)

    # Iterate on the HMM to query them. 
    results = []
    with zipfile.ZipFile(HMM_COMPRESS_FILE, 'r') as zip_object:
        for hmm_filename in zip_object.namelist():
            if hmm_filename.endswith('.hmm'):
                hmm_filebasename = os.path.basename(hmm_filename)
                with zip_object.open(hmm_filename) as open_hmm_zipfile:
                    with pyhmmer.plan7.HMMFile(open_hmm_zipfile) as hmm_file:
                        for hits in pyhmmer.hmmsearch(hmm_file, sequences, cpus=1):
                            for hit in hits:
                                results.append([input_filename, hit.name.decode(), hmm_filebasename, hit.evalue, hit.score, hit.length])

    return results


def write_results(hmm_results, output_file):
    with open(output_file, 'w') as open_output_file:
        csvwriter = csv.writer(open_output_file, delimiter='\t')
        csvwriter.writerow(['organism', 'protein', 'HMM', 'evalue', 'score', 'length'])
        for result in hmm_results:
            csvwriter.writerow(result)


def parse_result_files(hmm_output_folder):
    hmm_hits = {}
    for hmm_tsv_file in os.listdir(hmm_output_folder):
        hmm_output_filepath = os.path.join(hmm_output_folder, hmm_tsv_file)
        hmm_tsv_filename = hmm_tsv_file.replace('.tsv', '')
        hmm_hits[hmm_tsv_filename] = []

        with open(hmm_output_filepath, 'r') as open_result_file:
            csvreader = csv.DictReader(open_result_file, delimiter='\t')
            for line in csvreader:
                if float(line['evalue']) < 0.05:
                    hmm_hits[hmm_tsv_filename].append(line['HMM'])

    return hmm_hits


def create_major_functions(hmm_output_folder, output_file):
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
    with open(output_file, 'w') as open_output_file:
        csvwriter = csv.writer(open_output_file, delimiter='\t')
        csvwriter.writerow(['organism', *hmm_list_functions])
        for org in hmm_hits:
            present_functions = [len(set(hmm_functions[function]).intersection(set(hmm_hits[org])))/len(set(hmm_functions[function])) if len(set(hmm_functions[function]).intersection(set(hmm_hits[org]))) > 0 else 'NA' for function in hmm_list_functions]
            csvwriter.writerow([org, *present_functions])


def hmm_search_write_resutls(input_file_path, output_file):
    print('Search for HMMs on ' + input_file_path)
    hmm_results = query_fasta_file(input_file_path)
    write_results(hmm_results, output_file)

def search_hmm(input_variable, output_folder, cpu_number=1):
    input_dicts = file_or_folder(input_variable)

    hmm_output_folder = os.path.join(output_folder, 'hmm_results')
    is_valid_dir(hmm_output_folder)

    hmm_search_pool = Pool(processes=cpu_number)

    multiprocess_input_hmm_searches = []
    for input_file in input_dicts:
        input_filename = os.path.splitext(os.path.basename(input_file))[0]
        output_file = os.path.join(hmm_output_folder, input_filename + '.tsv')

        input_file_path = input_dicts[input_file]
        multiprocess_input_hmm_searches.append([input_file_path, output_file])

    hmm_search_pool.starmap(hmm_search_write_resutls, multiprocess_input_hmm_searches)

    hmm_search_pool.close()
    hmm_search_pool.join()

    function_matrix_file = os.path.join(output_folder, 'function_presence.tsv')
    create_major_functions(hmm_output_folder, function_matrix_file)

    hmm_diagram_folder = os.path.join(output_folder, 'diagram_input_folder')
    create_input_diagram(hmm_output_folder, hmm_diagram_folder)
