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

import logging
import os
import csv
import sys

logger = logging.getLogger(__name__)


def is_valid_dir(dirpath):
    """Return True if directory exists or can be created (then create it)
    
    Args:
        dirpath (str): path of directory

    Returns:
        bool: True if dir exists, False otherwise
    """
    if not os.path.isdir(dirpath):
        try:
            os.makedirs(dirpath)
            return True
        except OSError:
            return False
    else:
        return True


def file_or_folder(variable_folder_file, extension_checks=['.faa']):
    """Check if the variable is file or a folder

    Args:
        variable_folder_file (str): path to a file or a folder
        extension_checks (list): list of extension to keep

    Returns:
        dict: {name of input file: path to input file}
    """
    file_folder_paths = {}

    check_file = False
    if os.path.isfile(variable_folder_file):
        filename, file_extension = os.path.splitext(os.path.basename(variable_folder_file))
        if file_extension in extension_checks:
            file_folder_paths[filename] = variable_folder_file
            check_file = True

    check_folder = False
    # For folder, iterate through all files inside the folder.
    if os.path.isdir(variable_folder_file):
        for file_from_folder in os.listdir(variable_folder_file):
            filename, file_extension = os.path.splitext(os.path.basename(file_from_folder))
            if file_extension in extension_checks:
                file_folder_paths[filename] = os.path.join(variable_folder_file, file_from_folder)
                check_folder = True

    if check_file is False and check_folder is False:
        logger.critical('ERROR: Wrong input, either {0} does not exist, or no files with {1} extension are found.'.format(variable_folder_file, ','.join(extension_checks)))
        sys.exit(1)
    return file_folder_paths




def parse_result_files(hmm_output_folder):
    """Parse HMM search results and extract filtered hits.

    Args:
        hmm_output_folder (str): path to HMM search results folder (one tsv file per organism)

    Returns:
        hmm_hits (dict): dictionary with organism as key and list of hit HMMs as value
    """
    hmm_hits = {}
    for hmm_tsv_file in os.listdir(hmm_output_folder):
        hmm_output_filepath = os.path.join(hmm_output_folder, hmm_tsv_file)
        hmm_tsv_filename = hmm_tsv_file.replace('.tsv', '')
        hmm_hits[hmm_tsv_filename] = []

        with open(hmm_output_filepath, 'r') as open_result_file:
            csvreader = csv.DictReader(open_result_file, delimiter='\t')
            for line in csvreader:
                hmm_hits[hmm_tsv_filename].append(line['HMM'])

    return hmm_hits


def read_measures_file(measures_file_path):
    """Read measurement file (such as abundance file for samples). Expect a tsv or csv files with organisms as rows, samples as columns and abundance as values.

    Args:
        measurement_file_path (str): path to measure file

    Returns:
        column_measure (dict): for each column, subdict with the measure of the different rows.
        total_measure_per_column (dict): for each column, the total measure of all rows.
    """
    if measures_file_path.endswith('.tsv'):
        delimiter = '\t'
    elif measures_file_path.endswith('.csv'):
        delimiter = ','

    column_measure = {}
    with open(measures_file_path, 'r') as open_measures_file:
        csvreader = csv.DictReader(open_measures_file, delimiter=delimiter)
        headers = csvreader.fieldnames
        columns = headers[1:]
        first_row = headers[0]
        for row in csvreader:
            for column in columns:
                if column not in column_measure:
                    column_measure[column] = {}
                try:
                    column_measure[column][row[first_row]] = float(row[column])
                except:
                    column_measure[column][row[first_row]] = 0

    total_measure_per_column = {}
    for col in column_measure:
        total_measure_per_column[col] = sum([column_measure[col][row] for row in column_measure[col]])

    return column_measure, total_measure_per_column


def read_esmecata_proteome_file(proteome_tax_id_file):
    """Read esmecata proteome file to extract associated betwenn organism name and tax_id_name.

    Args:
        proteome_tax_id_file (str): path to proteome tax id file of esmecata

    Returns:
        observation_names_tax_id_names (dict): dictionary associating organism name with tax_id_name
        observation_names_tax_ranks (dict): dictionary associated organism name with tax_rank
    """
    observation_names_tax_id_names = {}
    observation_names_tax_ranks = {}

    with open(proteome_tax_id_file, 'r') as open_proteome_tax_id_file:
        csvreader = csv.DictReader(open_proteome_tax_id_file, delimiter='\t')
        for row in csvreader:
            observation_names_tax_id_names[row['observation_name']] = row['tax_id_name']
            observation_names_tax_ranks[row['observation_name']] = row['tax_rank']

    return observation_names_tax_id_names, observation_names_tax_ranks