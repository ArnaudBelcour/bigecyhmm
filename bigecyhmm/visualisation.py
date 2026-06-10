# Copyright (C) 2024-2026 Arnaud Belcour - Univ. Grenoble Alpes, Inria, Grenoble, France Microcosme
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

import pandas as pd
import numpy as np
from pandas import __version__ as pandas_version
import seaborn as sns
from seaborn import __version__ as seaborn_version
import matplotlib.pyplot as plt
from matplotlib import __version__ as matplotlib_version

import argparse
import logging
import csv
import json
import math
import os
import sys
import time

from bigecyhmm import __version__ as bigecyhmm_version
from bigecyhmm import PATHWAY_TEMPLATE_FILE, HMM_TEMPLATE_FILE, CUSTOM_HYDROGEN_TABLE, TEMPLATE_CUSTOM_CENTRAL_HYDROGEN, TEMPLATE_BACKGROUND_BIGECYHMM
from bigecyhmm.utils import is_valid_dir, read_measures_file, read_esmecata_proteome_file
from bigecyhmm.diagram_cycles import create_carbon_cycle, create_nitrogen_cycle, create_sulfur_cycle, create_other_cycle, create_phosphorus_cycle, get_diagram_pathways_hmms
from bigecyhmm.group_analysis import statNut_run

from esmecata.utils import get_domain_or_superkingdom_from_ncbi_tax_database

RANK_SORTED = ['isolate', 'strain', 'serotype', 'serogroup', 'forma', 'subvariety', 'varietas',
               'subspecies', 'forma specialis', 'species', 'species subgroup', 'species group',
               'subseries', 'series',
               'subsection', 'section',
               'subgenus', 'genus',
               'subtribe', 'tribe',
               'subfamily', 'family', 'superfamily',
               'parvorder', 'infraorder', 'suborder', 'order', 'superorder',
               'subcohort', 'cohort',
               'infraclass', 'subclass', 'class', 'superclass',
               'infraphylum', 'subphylum', 'phylum', 'superphylum',
               'subkingdom', 'kingdom', get_domain_or_superkingdom_from_ncbi_tax_database(),
               'clade', 'environmental samples', 'incertae sedis', 'unclassified', 'no rank', 'Not found']

FUNCTION_GROUP_TEMPLATE = {'Carbon cycle': ['C-S-01:Organic carbon oxidation', 'C-S-02:Carbon fixation', 'C-S-03:Ethanol oxidation', 'C-S-04:Acetate oxidation',
                                   'C-S-05:Hydrogen generation', 'C-S-06:Fermentation', 'C-S-07:Methanogenesis', 'C-S-08:Methanotrophy', 'C-S-09:Hydrogen oxidation',
                                   'C-S-10:Acetogenesis WL', 'Hydrogen generation', 'Hydrogen oxidation', 'Hydrogenotrophic methanogen.', 'Methanogenesis', 'Hydrogenotrophic acetogen.',
                                   'Acetate oxidation', 'Acetoclastic methanogenesis', 'Acetogen. WL', 'Carbon fixation', 'Methylotrophic methanogen.', 'Organic carbon oxidation',
                                   'CO2-dep. Methanogen.', '(Aerobic) methanotrophy', 'Fermentation', 'Ethanol oxidation', 'Fermentative H2 generation'],
                    'Sulfur cycle': ['S-S-01:Sulfide oxidation', 'S-S-02:Sulfur reduction', 'S-S-03:Sulfur oxidation', 'S-S-04:Sulfite oxidation', 'S-S-05:Sulfate reduction', 'S-S-06:Sulfite reduction',
                                    'S-S-07:Thiosulfate oxidation', 'S-S-08:Thiosulfate disproportionation 1', 'S-S-09:Thiosulfate disproportionation 2', 'S-S-10:Polysulfide reduction',
                                    'Sulfur respiration', 'Sulfate respiration', 'Sulfite oxidation', 'Sulfite reduction', 'Sulfur reduction', 'Polysulfide reduction', 'Sulfide oxidation',
                                    'Sulfur oxidation', 'Thiosulfate disproportionation 2', 'Thiosulfate oxidation', 'Thiosulfate disproportionation 1', 'Sulfate reduction (ASR or DSR)'],
                    'Nitrogen cycle': ['N-S-01:Nitrogen fixation', 'N-S-02:Ammonia oxidation', 'N-S-03:Nitrite oxidation', 'N-S-04:Nitrate reduction', 'N-S-05:Nitrite reduction', 'N-S-06:Nitric oxide reduction',
                                    'N-S-07:Nitrous oxide reduction', 'N-S-08:Nitrite ammonification', 'N-S-09:Anammox', 'N-S-10:Nitric oxide dismutase', 'Nitrate respiration / denitrification',
                                    'Assimilatory nitrate red.', 'Dissimilatory nitrate red.', 'Ammonia oxidation', 'Nitrite oxidation', 'Nitrite reduction', 'Nitric oxide dismutase',
                                    'Nitric oxide reduction', 'Nitrous oxide reduction', 'Anammox', 'Nitrite ammonification', 'Nitrogen fixation', 'Nitrate respiration or denitrification'],
                    'Other cycle': ['O-S-01:Iron reduction', 'O-S-02:Iron oxidation', 'O-S-03:Arsenate reduction', 'O-S-04:Arsenite oxidation', 'O-S-05:Selenate reduction', 'O-S-06:Aerobic respiration',
                                    'Arsenate reduction', 'Arsenite oxidation', 'e-input Metal respiration', 'Metal Respiration', 'Mn oxidation', 'Iron oxidation', 'Selenate reduction', 'Organohal. Respira. 1',
                                    'Organohal. Respira. 2', 'Hydrogen sensing'],
                    'Phosphorus cycle' : ['P-S-01:Immobilisation (P-rich)', 'P-S-01:Immobilisation (P-poor)', 'P-S-02:Mineralisation', 'P-S-03:Dissolution',
                                          'P dissolution', 'Immobilisation (P-rich)', 'Immobilisation (P-poor)', 'P mineralisation']}

MESSAGE = '''
Create figures from bigecyhmm and EsMeCaTa outputs (and optionally with an abundance file).
'''
REQUIRES = '''
Requires seaborn, pandas and matplotlib.
'''

logger = logging.getLogger()
logger.setLevel(logging.INFO)

def get_function_categories():
    """Extract function categories from HMM template file.

    Returns:
        function_categories (dict): category associated with function name
    """
    function_categories = {}
    with open(HMM_TEMPLATE_FILE, 'r') as open_hmm_template:
        csvreader = csv.DictReader(open_hmm_template, delimiter='\t')
        for line in csvreader:
            function_name = line['Function'] + ' ' + line['Gene abbreviation']
            category_name = line['Category']
            if category_name not in function_categories:
                function_categories[category_name] = [function_name]
            else:
                function_categories[category_name].append(function_name)

    return function_categories


def compute_relative_abundance_per_tax_id(sample_abundance, sample_tot_abundance, observation_names_tax_id_names):
    """For each tax_id_name selected by esmecata (from observation_names_tax_id_names) compute the relative abundace of this taxon.
    It is done by summing the abundance of all organisms in this tax_id_name and then dividing it by the total abundance in the sample.

    Args:
        sample_abundance (dict): for each sample, subdict with the abundance of the different organisms.
        sample_tot_abundance (dict): for each sample, the total abundance of all organisms in the sample.
        observation_names_tax_id_names (dict): dictionary associating organism name with tax_id_name

    Returns:
        abundance_data (dict): for each sample, contains a subdict with the relative abundance of tax_id_name in these samples.
    """
    abundance_data = {}

    for sample_name in sample_abundance:
        for observation_name in sample_abundance[sample_name]:
            if observation_name in observation_names_tax_id_names:
                tax_id_name = observation_names_tax_id_names[observation_name]
                if sample_name not in abundance_data:
                    abundance_data[sample_name] = {}
                if tax_id_name not in abundance_data[sample_name]:
                    abundance_data[sample_name][tax_id_name] = float(sample_abundance[sample_name][observation_name])
                else:
                    abundance_data[sample_name][tax_id_name] = float(sample_abundance[sample_name][observation_name]) + float(abundance_data[sample_name][tax_id_name])

        for tax_id_name in abundance_data[sample_name]:
            abundance_data[sample_name][tax_id_name] = abundance_data[sample_name][tax_id_name] / sample_tot_abundance[sample_name]

    return abundance_data


def compute_abundance_per_tax_rank(sample_abundance, observation_names_tax_ranks, sample_tot_abundance):
    """For each tax_rank selected by esmecata (from observation_names_tax_id_names) compute the relative abundace of this taxon.
    It is done by summing the abundance of all organisms in this tax_id_name and then dividing it by the total abundance in the sample.

    Args:
        sample_abundance (dict): for each sample, subdict with the abundance of the different organisms.
        observation_names_tax_id_names (dict): dictionary associating organism name with tax_id_name
        sample_tot_abundance (dict): for each sample, the total abundance of all organisms in the sample.

    Returns:
        data_abundance_taxon_sample (dict): for each sample, contains a subdict with the relative abundance of tax_rank in these samples.
        sample_abundance_tax_rank (list): list of list containing the abundance of each tax_rank in the different samples.
        data_abundance_organism_sample (list): list of list containing the abundance of each organism in the different samples.
    """
    sample_abundance_tax_rank = {}
    organism_abundance_tax_rank = {}
    for sample in sample_abundance:
        if sample not in sample_abundance_tax_rank:
            sample_abundance_tax_rank[sample] = {}
        if sample not in organism_abundance_tax_rank:
            organism_abundance_tax_rank[sample] = {}
        for organism in sample_abundance[sample]:
            abundance_organism = sample_abundance[sample][organism]
            if organism in observation_names_tax_ranks:
                tax_rank = observation_names_tax_ranks[organism]
            else:
                tax_rank = 'Not found'
            if tax_rank not in sample_abundance_tax_rank[sample]:
                sample_abundance_tax_rank[sample][tax_rank] = abundance_organism
            else:
                sample_abundance_tax_rank[sample][tax_rank] = sample_abundance_tax_rank[sample][tax_rank] + abundance_organism
            if organism not in organism_abundance_tax_rank[sample]:
                organism_abundance_tax_rank[sample][organism] = (tax_rank, abundance_organism)
            else:
                organism_abundance_tax_rank[sample][organism][1] = organism_abundance_tax_rank[sample][organism][1] + abundance_organism

    data_abundance_taxon_sample = []
    for sample in sample_abundance_tax_rank:
        for tax_rank in sample_abundance_tax_rank[sample]:
            data_abundance_taxon_sample.append([sample, tax_rank, sample_abundance_tax_rank[sample][tax_rank]/sample_tot_abundance[sample]])

    data_abundance_organism_sample = []
    for sample in organism_abundance_tax_rank:
        for organism in organism_abundance_tax_rank[sample]:
            tax_rank = organism_abundance_tax_rank[sample][organism][0]
            data_abundance_organism_sample.append([sample, organism, tax_rank, organism_abundance_tax_rank[sample][organism][1]/sample_tot_abundance[sample]])

    return data_abundance_taxon_sample, sample_abundance_tax_rank, data_abundance_organism_sample


def compute_bigecyhmm_functions_occurrence(bigecyhmm_output_file, tax_id_names_observation_names=None):
    """Read pathway_presence.tsv or function_presence.tsv created by bigecyhmm to compute the occurrence of each functions/pathways.

    Args:
        bigecyhmm_output_file (str): path to the output file of bigecyhmm (either pathway_presence.tsv or function_presence.tsv).
        tax_id_names_observation_names (dict): dictionary associating tax_id_name with organism name.

    Returns:
        function_occurrence_organisms (dict): dictionary containing function as key and subdict with organism as key and value of function in organism.
        all_studied_organisms (list): list of all organisms in community.
    """
    bigecyhmm_function_df = pd.read_csv(bigecyhmm_output_file, sep='\t')
    bigecyhmm_function_df.set_index('function', inplace=True)

    all_studied_organisms = bigecyhmm_function_df.columns
    # Get all tax_id_names in the community then the observation_names associated with them (in the case of run with esmecata results).
    if tax_id_names_observation_names is not None:
        all_studied_organisms = list(set([observation_name for tax_id in all_studied_organisms for observation_name in tax_id_names_observation_names[tax_id]]))

    # For each function, count the number of organisms predicted to have it.
    all_functions = []
    function_occurrence_organisms = {}
    for function_name, row in bigecyhmm_function_df.iterrows():
        all_functions.append(function_name)
        for organism in bigecyhmm_function_df.columns:
            if math.isnan(row[organism]):
                row[organism] = 0
            else:
                row[organism] = int(row[organism])
            # If results come from esmecata, convert tax_id_names into observation_names.
            if tax_id_names_observation_names is not None:
                observation_names = tax_id_names_observation_names[organism]
            else:
                observation_names = [organism]

            if row[organism] > 0:
                if function_name not in function_occurrence_organisms:
                    function_occurrence_organisms[function_name] = {}
                for observation_name in observation_names:
                    if observation_name not in function_occurrence_organisms[function_name]:
                        function_occurrence_organisms[function_name][observation_name] = row[organism]

    all_functions = set(all_functions)
    for function in all_functions:
        if function not in function_occurrence_organisms:
            function_occurrence_organisms[function] = {}

    return function_occurrence_organisms, all_studied_organisms


def compute_bigecyhmm_functions_abundance(bigecyhmm_output_file, sample_abundance, sample_tot_abundance, tax_id_names_observation_names=None):
    """Read pathway_presence.tsv or function_presence.tsv created by bigecyhmm to compute the occurrence of each functions/pathways.

    Args:
        bigecyhmm_output_file (str): path to the output file of bigecyhmm (either pathway_presence.tsv or function_presence.tsv).
        sample_abundance (dict): for each sample, subdict with the abundance of the different organisms.
        sample_tot_abundance (dict): for each sample, the total abundance of all organisms in the sample.
        tax_id_names_observation_names (dict): dictionary associating tax_id_name with organism name.

    Returns:
        function_abundance_samples (dict): dictionary containing sample as dict and a subdict containing function associated with abundance
        function_relative_abundance_samples (dict): dictionary containing sample as dict and a subdict containing function associated with relative abundance
        function_participation_samples (dict): dictionary containing sample as dict and a subdict containing organism associated with function and their abundance
    """
    bigecyhmm_function_df = pd.read_csv(bigecyhmm_output_file, sep='\t')
    bigecyhmm_function_df.set_index('function', inplace=True)

    sample_abundance_dataframe = pd.DataFrame(sample_abundance)

    # Compute the occurrence of functions in organism from bigecyhmm file.
    function_organisms, all_studied_organisms = compute_bigecyhmm_functions_occurrence(bigecyhmm_output_file, tax_id_names_observation_names)
    # Transpose the dataframe to get function as rows and organisms as columns.
    function_dataframe = pd.DataFrame(function_organisms).T
    # Replace nan by 0.
    function_dataframe = function_dataframe.fillna(0)

    missing_organism_abundance = list(set(function_dataframe.columns) - set(sample_abundance_dataframe.index))
    if len(missing_organism_abundance) > 0:
        logger.critical('ERROR: Several organisms ({0}) having predicted functions are not present in abundance file.'.format(','.join(missing_organism_abundance)))
        sys.exit(1)
    # Ensure that indexes in sample_abundance_dataframe are the same than the columns in function_dataframe.
    # There can be more indexes in sample_abundance_dataframe if some organisms do not have functional predictions.
    sample_abundance_dataframe = sample_abundance_dataframe.reindex(function_dataframe.columns)
    # Replace nan by 0.
    sample_abundance_dataframe = sample_abundance_dataframe.fillna(0)

    # Matrix multiplication between function matrix and abundance matrix.
    abundance_function_df = function_dataframe.dot(sample_abundance_dataframe)

    # Matrix division by the total abundance in each sample to get relative abundace
    relative_abundance_df = abundance_function_df.div(sample_tot_abundance)

    function_abundance_samples = abundance_function_df.to_dict()
    function_relative_abundance_samples = relative_abundance_df.to_dict()

    # For each sample compute the abundance of function according to the organisms.
    function_participation_samples = {}
    for sample in sample_abundance:
        # Multiply each function by the abundance of the organism.
        function_participation_df = function_dataframe.mul(sample_abundance_dataframe[sample])
        # Remove row (organism) containing only zero.
        function_participation_df = function_participation_df.loc[~(function_participation_df == 0).all(axis=1)]
        # Remove column (function) containing only zero.
        function_participation_df = function_participation_df[function_participation_df.columns[~(function_participation_df == 0).all(axis=0)]]
        function_participation_samples[sample] = function_participation_df.to_dict()

    return function_abundance_samples, function_relative_abundance_samples, function_participation_samples


def get_hmm_per_organism(bigecyhmm_output, tax_id_names_observation_names):
    """ Get for each organism, the number of matching HMMs to their proteomes.

    Args:
        bigecyhmm_output_file (str): path to the output file of bigecyhmm (either pathway_presence.tsv or function_presence.tsv).
        tax_id_names_observation_names (dict): dictionary associating tax_id_name with organism name.

    Returns:
        hmm_occurrences (dict): dictionary containing organism as value and a subdict containing HMM occurrence
    """
    hmm_found_folder = os.path.join(bigecyhmm_output, 'hmm_results')
    hmm_occurrences = {}

    for organism_result_file in os.listdir(hmm_found_folder):
        organism_name = organism_result_file.replace('.tsv', '')
        # If results come from esmecata, convert tax_id_names into observation_names.
        if tax_id_names_observation_names is not None:
            observation_names = tax_id_names_observation_names[organism_name]
        else:
            observation_names = [organism_name]

        organism_result_file_path = os.path.join(hmm_found_folder, organism_result_file)
        hmm_found_df = pd.read_csv(organism_result_file_path, sep='\t')
        hmm_found_df['HMM'] = hmm_found_df['HMM'].str.replace('.hmm', '')
        # Group the dataframe by the hmm and merged all proteins found for an HMM with a ';'.
        hmm_founds = hmm_found_df.groupby('HMM').apply(lambda x: '; '.join(x.protein), include_groups=False).to_dict()
        # Count the number of protein found for each HMM.
        hmm_occurrence_org = {hmm: len(set(hmm_founds[hmm])) for hmm in hmm_founds}

        for observation_name in observation_names:
            hmm_occurrences[observation_name] = hmm_occurrence_org

    return hmm_occurrences


def create_ko_functional_profile(bigecyhmm_output, sample_abundance, output_folder_abundance, tax_id_names_observation_names):
    """ Compute HMM abundance in each sample and write it in a tsv file..

    Args:
        bigecyhmm_output_file (str): path to the output file of bigecyhmm (either pathway_presence.tsv or function_presence.tsv).
        sample_abundance (dict): for each sample, subdict with the abundance of the different organisms.
        output_folder_abundance (str): path to write the output file.
        tax_id_names_observation_names (dict): dictionary associating tax_id_name with organism name.
    """
    hmm_occurrences = get_hmm_per_organism(bigecyhmm_output, tax_id_names_observation_names)

    hmm_abundance = {}
    for sample in sample_abundance:
        if sample not in hmm_abundance:
            hmm_abundance[sample] = {}
        for organism in sample_abundance[sample]:
            if organism in hmm_occurrences:
                for hmm_found in hmm_occurrences[organism]:
                    # Compute the abundance of function in the sample.
                    if hmm_found not in hmm_abundance[sample]:
                        if hmm_occurrences[organism][hmm_found] > 0:
                            hmm_abundance[sample][hmm_found] = 1*sample_abundance[sample][organism]
                    else:
                        if hmm_occurrences[organism][hmm_found] > 0:
                            hmm_abundance[sample][hmm_found] = 1*sample_abundance[sample][organism] + hmm_abundance[sample][hmm_found]


    hmm_abundance_df = pd.DataFrame(hmm_abundance)
    hmm_abundance_df.index.name = 'function'
    hmm_abundance_df.sort_index(inplace=True)
    hmm_abundance_df.to_csv(os.path.join(output_folder_abundance, 'hmm_functional_profile.tsv'), sep='\t')


def create_swarmplot_community(df_seaborn, output_file_name):
    """Create swarmplot from pandas dataframe with function as 'name' column' and associated ratio as a second column

    Args:
        df_seaborn (pd.DataFrame): dataframe pandas containing a column with the name of function, a second column with the ratio of organisms having it in the community and a third column for the sample
        output_file_name (path): path to the output file.
    """
    fig, axes = plt.subplots(figsize=(40,20))
    plt.rc('font', size=30)
    ax = sns.swarmplot(data=df_seaborn, x='name', y='ratio', s=10)
    [ax.axvline(x+.5,color='k') for x in ax.get_xticks()]
    plt.xticks(rotation=90)
    plt.savefig(output_file_name, bbox_inches="tight")
    plt.clf()


def create_swarmplot_sample(df_seaborn, output_file_name):
    """Create swarmplot from pandas dataframe with function as 'name' column' and associated ratio as a second column

    Args:
        df_seaborn (pd.DataFrame): dataframe pandas containing a column with the name of function, a second column with the ratio of organisms having it in the community and a third column for the sample
        output_file_name (path): path to the output file.
    """
    fig, axes = plt.subplots(figsize=(40,20))
    plt.rc('font', size=30)
    ax = sns.swarmplot(data=df_seaborn, x='name', y='ratio', hue='sample', s=10)
    [ax.axvline(x+.5,color='k') for x in ax.get_xticks()]
    plt.xticks(rotation=90)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(output_file_name, bbox_inches="tight")
    plt.clf()


def create_boxplot_sample(df_seaborn, output_file_name):
    """Create boxplot from pandas dataframe with function as 'name' column' and associated ratio as a second column

    Args:
        df_seaborn (pd.DataFrame): dataframe pandas containing a column with the name of function, a second column with the ratio of organisms having it in the community and a third column for the sample
        output_file_name (path): path to the output file.
    """
    fig, axes = plt.subplots(figsize=(40,20))
    plt.rc('font', size=30)
    ax = sns.boxplot(data=df_seaborn, x='name', y='ratio', hue='sample')
    [ax.axvline(x+.5,color='k') for x in ax.get_xticks()]
    plt.xticks(rotation=90)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(output_file_name, bbox_inches="tight")
    plt.clf()


def create_polar_plot(df_seaborn_abundance, output_file):
    fig, ax = plt.subplots(1, 1, figsize=(10, 18), subplot_kw={'projection': 'polar'})

    removed_functions = ['N-S-10:Nitric oxide dismutase', 'S-S-10:Polysulfide reduction']

    kept_functions = [name for name in df_seaborn_abundance['name']
                        if df_seaborn_abundance[df_seaborn_abundance['name']==name]['ratio'].max()>0.1]
    kept_functions = [kept_func for kept_func in kept_functions if not kept_func.startswith('P-S')]
    kept_functions = [kept_func for kept_func in kept_functions if not kept_func.startswith('P1-S')]
    kept_functions = [kept_func for kept_func in kept_functions if not 'Aerobic' in kept_func]
    kept_functions = [kept_func for kept_func in kept_functions if not 'P-S' in kept_func]

    df_seaborn_abundance = df_seaborn_abundance.sort_values(['name'], ascending=False)
    # Remove function
    df_seaborn_abundance = df_seaborn_abundance[~df_seaborn_abundance['name'].isin(removed_functions)]
    df_seaborn_abundance = df_seaborn_abundance[df_seaborn_abundance['name'].isin(kept_functions)]

    #tmp_df_seaborn_abundance = tmp_df_seaborn_abundance.groupby(['sample_nb', 'name'], as_index=False)['ratio'].mean()
    mean_data = df_seaborn_abundance.groupby('name')['ratio'].mean()
    max_data = df_seaborn_abundance.groupby('name')['ratio'].max()
    min_data = df_seaborn_abundance.groupby('name')['ratio'].min()
    color = 'blue'
    if not df_seaborn_abundance.empty and 'name' in df_seaborn_abundance.columns:
        theta = [function_name.split(':')[1] if ':' in function_name else function_name for function_name in mean_data.index]
        # Add several angles to the polar plot according to the number of functions.
        funciton_angles = np.linspace(0.05, 2 * np.pi-0.05, len(theta), endpoint=False)
        ax.set_xticks(funciton_angles)
        # Add the function labels.
        ax.set_xticklabels(theta, size=6)
        # Plot mean value.
        r = mean_data.tolist()
        ax.plot(funciton_angles, r, color='grey')
        # Plot max value.
        r_max = max_data.tolist()
        ax.plot(funciton_angles, r_max, color=color, linewidth=0.5)
        # Plot min value.
        r_min = min_data.tolist()
        ax.plot(funciton_angles, r_min, color=color, linewidth=0.5)
        # Add color fill in between r_max-r and r-r_min.
        ax.fill_between(funciton_angles, r, r_min, alpha=0.2, facecolor=color)
        ax.fill_between(funciton_angles, r_max, r, alpha=0.2, facecolor=color)
        # Set min and max for r values.
        ax.set_rmax(1)
        ax.set_rmin(0)
    plt.tight_layout()
    plt.savefig(output_file, bbox_inches="tight")
    plt.close(fig)


def generate_bubble_plot(melted_cycle_relative_abundance_samples_df, output_file, group_file=None):
    """ Generate a bubble plot showing main functions abundance in pathways.

    Args:
        melted_cycle_relative_abundance_samples_df (pd.DataFrame): abundance of pathways according to samples.
        output_file_name (str): path to the output file.
        group_file_path (str): path to group file.
    """
    function_in_cycle = {function_name: cycle_name for cycle_name in FUNCTION_GROUP_TEMPLATE for function_name in FUNCTION_GROUP_TEMPLATE[cycle_name]}
    function_in_cycle_set = set(function_in_cycle.keys())

    if set(melted_cycle_relative_abundance_samples_df['name'].unique()).issubset(function_in_cycle_set):
        melted_cycle_relative_abundance_samples_df['group'] = [function_in_cycle[function_name] if function_name in function_in_cycle else '' for function_name in melted_cycle_relative_abundance_samples_df['name']]
    else:
        melted_cycle_relative_abundance_samples_df['group'] = 'Function'
    tmp_melted_cycle_relative_abundance_samples_df = melted_cycle_relative_abundance_samples_df[melted_cycle_relative_abundance_samples_df['group']!='']

    function_groups = tmp_melted_cycle_relative_abundance_samples_df['group'].unique()
    tmp_melted_cycle_relative_abundance_samples_df = tmp_melted_cycle_relative_abundance_samples_df[tmp_melted_cycle_relative_abundance_samples_df['ratio']>0.05]
    # Keep only function groups with predictions.
    function_groups = [function_group for function_group in function_groups
                       if tmp_melted_cycle_relative_abundance_samples_df[tmp_melted_cycle_relative_abundance_samples_df['group']==function_group].empty is False]

    # If group file is not None, sort dataframe with sample according to group.
    # And generate position of group labels and line splitting samples according to group.
    if group_file is not None:
        mapping_df = pd.read_csv(group_file, sep='\t', dtype=str)

        mapping_df = mapping_df.sort_values('group')
        sample_groups = mapping_df['group'].unique()
        sorted_samples = mapping_df['sample'].tolist()
        length_sample_name = int(max([len(sample_name) for sample_name in sorted_samples])/1.5)

        x_group_position = []
        x_end_group_position = []
        for group in sample_groups:
            tmp_group_mapping_df = mapping_df[mapping_df['group']==group]
            group_samples = tmp_group_mapping_df['sample'].unique().tolist()
            # Get sample of the group at the middle of the sample position.
            if len(group_samples) > 1:
                x_group_position.append(group_samples[math.ceil(len(group_samples)/2)])
            else:
                x_group_position.append(group_samples[0])
            x_end_group_position.append(group_samples[-1])
        x_end_group_position.append(sorted_samples[-1])
        tmp_melted_cycle_relative_abundance_samples_df.sort_values(by="sample", key=lambda x: x.map(sorted_samples.index), inplace=True)

    nb_samples = len(tmp_melted_cycle_relative_abundance_samples_df['sample'].unique())
    nb_unique_functions = len(tmp_melted_cycle_relative_abundance_samples_df['name'].unique())

    # Compute height ratios according to number of function groups.
    ratios = [len(tmp_melted_cycle_relative_abundance_samples_df[tmp_melted_cycle_relative_abundance_samples_df['group']==function_group]['name'].unique()) for function_group in function_groups]

    fig_width = 15
    if nb_samples > 30:
        fig_width = nb_samples / 2
    fig_height = 15
    if nb_unique_functions > 30:
        fig_height = nb_unique_functions / 2
    fig, axes = plt.subplots(nrows=len(function_groups), ncols=1, figsize=(fig_width, fig_height), gridspec_kw={'height_ratios': ratios})

    if len(function_groups) > 1:
        for index, function_group in enumerate(function_groups):
            tmp_data = tmp_melted_cycle_relative_abundance_samples_df[tmp_melted_cycle_relative_abundance_samples_df['group']==function_group]
            # Create bubble scatter plot.
            axes[index].scatter(tmp_data['sample'], tmp_data['name'], s=tmp_data['ratio']*500, c=tmp_data['ratio'], cmap='viridis_r', alpha=0.8)
            # Show grid.
            axes[index].grid(True, color='lightgrey', linewidth=0.5)
            # Remove tick labels except for the last one.
            if index < len(function_groups)-1:
                axes[index].set_xticklabels([])
            last_axe = axes[index]
    else:
        tmp_data = tmp_melted_cycle_relative_abundance_samples_df
        # Create bubble scatter plot.
        axes.scatter(tmp_data['sample'], tmp_data['name'], s=tmp_data['ratio']*500, c=tmp_data['ratio'], cmap='viridis_r', alpha=0.8)
        # Show grid.
        axes.grid(True, color='lightgrey', linewidth=0.5)
        last_axe = axes

    if group_file is not None:
        # Add group label below sample names.
        x_group_position = [sorted_samples.index(sample_name) for sample_name in x_group_position]
        sec = last_axe.secondary_xaxis(location=0)
        group_labels = ['\n'*length_sample_name+sample_group for sample_group in sample_groups]
        sec.set_xticks(x_group_position, labels=group_labels)
        sec.tick_params('x', length=0)
        # Add lines between groups.
        x_end_group_position = [sorted_samples.index(sample_name)+0.5 if sample_name!=sorted_samples[-1] else len(sorted_samples)-0.5 for sample_name in x_end_group_position]
        sec2 = last_axe.secondary_xaxis(location=0)
        sec2.set_xticks([-0.5, *x_end_group_position], labels=[])
        sec2.tick_params('x', length=40, width=1.5)

    fig.supylabel('Functions above 5% of relative abundance')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)


def visualise_barplot_category(category, function_categories, df_seaborn_abundance, output_file_name):
    """Create bar plot for functions of the associated categories from pandas dataframe with function as 'name' column' and associated ratio as a second column

    Args:
        cateogry (str or list): name of the function category to plot or list of categories
        function_categories (dict): a dicitonary mapping function category to their respective function inferred by bigecyhmm
        df_seaborn_abundance (pd.DataFrame): dataframe pandas containing a column with the name of function, a second column with the relative abundance of organisms having it in the community and a third column for the sample
        output_file_name (str): path to the output file.
    """
    fig, axes = plt.subplots(figsize=(40,20))
    plt.rc('font', size=30)
    if isinstance(category, str):
        kept_functions = function_categories[category]
    elif isinstance(category, list):
        kept_functions = []
        for cat in category:
            kept_functions.extend(function_categories[cat])

    df_seaborn_abundance = df_seaborn_abundance[df_seaborn_abundance['name'].isin(kept_functions)]

    g = sns.barplot(data=df_seaborn_abundance, x='name', y='ratio', hue='sample')
    plt.xticks(rotation=90)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(output_file_name, bbox_inches="tight")
    plt.clf()


def create_heatmap_functions(df, output_heatmap_filepath):
    """Create heatmap of function abundances in samples

    Args:
        df (pd.DataFrame): dataframe pandas containing a column with the name of function, one column by sample and the abundance of function in sample as value.
        output_heatmap_filepath (str): path to the output file.
    """
    sns.set_theme(font_scale=0.5)
    fig, axes = plt.subplots(figsize=(35,70))
    g = sns.heatmap(data=df, xticklabels=1, cmap='viridis_r', linewidths=1, linecolor='black', square=True)
    plt.tight_layout()
    plt.savefig(output_heatmap_filepath)
    plt.clf()


def add_abundance_and_measure_to_graph(graph_file, output_folder, metabolite_measure=None, cycle_relative_abundance_samples_df=None):
    import networkx as nx
    pathway_graph = nx.read_graphml(graph_file)

    # Initiate a bipartite network when giving abundance or measure files.
    if cycle_relative_abundance_samples_df is not None or metabolite_measure is not None:
        abundance_cycle_network = nx.DiGraph()
        nodes_data = [(node, pathway_graph.nodes[node]) for node in pathway_graph.nodes]
        abundance_cycle_network.add_nodes_from(nodes_data)
        for edge in pathway_graph.edges:
            abundance_cycle_network.add_edges_from([edge], **pathway_graph.edges[edge])

    if cycle_relative_abundance_samples_df is not None:
        for function in cycle_relative_abundance_samples_df.index.tolist():
            if function not in abundance_cycle_network.nodes:
                abundance_cycle_network.add_node(function, type='Function')
            for sample in cycle_relative_abundance_samples_df.columns:
                pathway_abundance = cycle_relative_abundance_samples_df[sample].loc[function]
                abundance_cycle_network.nodes[function][sample] = pathway_abundance

    if metabolite_measure is not None:
        # Add sample metabolite measurements to node.
        if metabolite_measure is not None:
            sample_metabolite_measure, sample_tot_abundance = read_measures_file(metabolite_measure)
            for metabolite in abundance_cycle_network.nodes:
                for sample in sample_metabolite_measure:
                    if metabolite in sample_metabolite_measure[sample]:
                        metabolite_measure = sample_metabolite_measure[sample][metabolite]
                        if isinstance(metabolite_measure, float) or isinstance(metabolite_measure, int):
                            abundance_cycle_network.nodes[metabolite][sample] = metabolite_measure

    if cycle_relative_abundance_samples_df is not None or metabolite_measure is not None:
        abundance_cycle_network.add_edges_from(pathway_graph.edges)
        network_graphml_output_file = os.path.join(output_folder, 'cycle_diagram_bipartite_abundance.graphml')
        nx.write_graphml(abundance_cycle_network, network_graphml_output_file)


def generate_barplot_esmecata_taxon_abundance(sample_abundance, observation_names_tax_ranks, sample_tot_abundance,
                                              output_folder_abundance):
    # Compute and create a visualisation of the abundance of selected taxonomic rank by esmecata in the different samples
    # This allows to identify the coverage of taxon found by esmecata compared to all the organisms in the sample.
    data_abundance_taxon_sample, sample_abundance_tax_rank, data_abundance_organism_sample = compute_abundance_per_tax_rank(sample_abundance, 
                                                                                                                            observation_names_tax_ranks, sample_tot_abundance)

    df_abundance_organism_sample = pd.DataFrame(data_abundance_organism_sample, columns=['Sample', 'Organism_name', 'Taxonomic rank selected by EsMeCaTa', 'Relative abundance'])
    output_organism_abundance_file = os.path.join(output_folder_abundance, 'barplot_esmecata_found_organism_sample.tsv')
    df_abundance_organism_sample.to_csv(output_organism_abundance_file, sep='\t', index=0)
    output_missing_organism_abundance_file = os.path.join(output_folder_abundance, 'barplot_esmecata_missing_organism_sample.tsv')
    df_abundance_organism_sample = df_abundance_organism_sample[df_abundance_organism_sample['Relative abundance']>0]
    df_abundance_organism_sample[df_abundance_organism_sample['Taxonomic rank selected by EsMeCaTa'] == 'Not found'].to_csv(output_missing_organism_abundance_file, sep='\t', index=0)

    df_abundance_taxon_sample = pd.DataFrame(data_abundance_taxon_sample, columns=['Sample', 'Taxonomic rank selected by EsMeCaTa', 'Relative abundance'])
    # Sort the dataframe using taxonomic rank, from lowest (species, genus) to highest (kingdom).
    df_abundance_taxon_sample.sort_values(by="Taxonomic rank selected by EsMeCaTa", key=lambda column: column.map(lambda e: RANK_SORTED.index(e)), inplace=True)
    output_taxon_rank_abundance_plot_file = os.path.join(output_folder_abundance, 'barplot_esmecata_found_taxon_sample.png')

    df_abundance_taxon_sample = df_abundance_taxon_sample.set_index(['Sample', 'Taxonomic rank selected by EsMeCaTa'])['Relative abundance'].unstack().reset_index()
    df_abundance_taxon_sample.set_index('Sample', inplace=True)

    color_ranks = {'species': '#EF553B', 'genus': '#00CC96', 'family': '#AB63FA', 'order': '#FFA15A', 'class': '#19D3F3',
                'phylum': '#FF6692', 'Not found': 'grey'}
    unused_colors = ['#B6E880', '#FF97FF', '#FECB52' '#636EFA']
    added_color = 0
    for rank in df_abundance_taxon_sample.columns:
        if rank not in color_ranks:
            color_ranks[rank] = unused_colors[added_color]
            added_color += 1

    fig_width = 16
    nb_samples = len(df_abundance_taxon_sample.index.unique())
    if nb_samples > 150:
        fig_width = nb_samples / 10
    # Sort the dataframe using taxonomic rank, from lowest (species, genus) to highest (kingdom).
    sorted_columns = [rank for rank in RANK_SORTED if rank in df_abundance_taxon_sample.columns]
    df_abundance_taxon_sample = df_abundance_taxon_sample[sorted_columns]
    df_abundance_taxon_sample.plot(kind='bar', stacked=True, color=color_ranks, figsize=(fig_width, 14))
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    plt.style.use('default')
    plt.savefig(output_taxon_rank_abundance_plot_file, bbox_inches="tight", transparent=False)
    plt.clf()


def taxon_function_heatmap(df_cycle_occurrence_organisms, proteome_tax_id_file, sample_abundance, output_folder):
    df_abundance = pd.DataFrame(sample_abundance)
    for col in df_abundance.columns:
        df_abundance[col] = df_abundance[col] / df_abundance[col].sum()
    nb_samples = len(df_abundance.columns)
    for cycle_name in df_cycle_occurrence_organisms.columns:
        tmp_df_cycle_occurrence_organisms = df_cycle_occurrence_organisms[df_cycle_occurrence_organisms[cycle_name]>0].copy()

        cycle_organisms = tmp_df_cycle_occurrence_organisms.index.tolist()

        cycle_organisms = set(cycle_organisms)

        proteome_tax_id_df = pd.read_csv(proteome_tax_id_file, sep='\t')
        proteome_tax_id_df.set_index('observation_name', inplace=True)
        asv_to_names = proteome_tax_id_df['name'].to_dict()

        tmp_df_abundance = df_abundance[df_abundance.index.isin(cycle_organisms)].copy()
        if tmp_df_abundance.empty is False:
            tmp_df_abundance['taxon'] = [asv_to_names[index] for index in tmp_df_abundance.index]
            nb_taxa = len(tmp_df_abundance['taxon'].unique())
            tmp_df_abundance = tmp_df_abundance.groupby('taxon').sum()
            tmp_df_abundance.reset_index(inplace=True)
            tmp_df_abundance.set_index(['taxon'], inplace=True)
            tmp_df_abundance = tmp_df_abundance.loc[~(tmp_df_abundance == 0).all(axis=1)]
            tmp_df_abundance.loc['Sum ' + cycle_name] = tmp_df_abundance.sum(axis=0)
            sns.set_theme(font_scale=1.1, style='white')
            fig_height = 12
            if nb_taxa > 36:
                fig_height = nb_taxa / 3
            fig_width = 19
            if nb_samples > 133:
                fig_width = nb_samples / 7
            fig, axes = plt.subplots(figsize=(fig_width, fig_height))
            g = sns.heatmap(data=tmp_df_abundance, center=1, xticklabels=True, yticklabels=True, cmap='viridis_r', linewidths=1, linecolor='black', square=False, mask=(tmp_df_abundance==0), vmin=0, vmax=1)

            plt.tight_layout()
            cycle_name_heatmap_file = os.path.join(output_folder, cycle_name+'.png')
            plt.savefig(cycle_name_heatmap_file)
            cycle_name_dataframe_file = os.path.join(output_folder, cycle_name+'.tsv')
            tmp_df_abundance.to_csv(cycle_name_dataframe_file, sep='\t')
            plt.show()
            plt.clf()
            plt.close(fig)


def generate_graph_figure(bigecyhmm_database_folder, graph_output_file):
    """ Generate graph background figure for donut plot.

    Args:
        bigecyhmm_database_folder (str): path to bigecyhmm database output folder (containing reference graphml file).
        graph_output_file (str): path to output background image.
    """
    from networkx.drawing.nx_agraph import graphviz_layout
    import networkx as nx

    graph_file = os.path.join(bigecyhmm_database_folder, 'input_graph.graphml')

    bipartite_networkx_graph = nx.read_graphml(graph_file)

    pathway_template_file = os.path.join(bigecyhmm_database_folder, 'pathway_template_file.tsv')
    pathway_template_df = pd.read_csv(pathway_template_file, sep='\t')

    pathway_template_df['index'] = pathway_template_df.index+1
    pathway_index_labels = pathway_template_df.set_index('Pathways')['index'].to_dict()

    metabolite_nodes = []

    metabolite_nodes = [network_node for network_node in bipartite_networkx_graph.nodes if bipartite_networkx_graph.nodes[network_node]['type'] == 'Metabolite']
    function_nodes = [network_node for network_node in bipartite_networkx_graph.nodes if bipartite_networkx_graph.nodes[network_node]['type'] == 'Function']
    for metabolite_node in metabolite_nodes:
        pathway_index_labels[metabolite_node] = metabolite_node

    nitrate_nodes = ['NO2-', 'NO', 'N2O', 'N2', 'NO3-', 'NH3']
    sulfur_nodes = ['S', 'H2S', 'SO3 2-', 'S2O3 2-', 'SO4 2-']
    organic_carbon_nodes = ['CH4', 'Methyl-CoM', 'Acetate',	'CO2', 'Organic C', 'Ethanol']

    node_colors = {}
    for node in metabolite_nodes:
        if node in nitrate_nodes:
            node_colors[node] = 'green'
        if node in sulfur_nodes:
            node_colors[node] = 'orange'
        if node in organic_carbon_nodes:
            node_colors[node] = 'gray'

    metabolite_node_colors = [node_colors[node] if node in node_colors else 'lightgray' for node in metabolite_nodes]

    fig, axes = plt.subplots(figsize=(12,8))

    layout_prog = 'neato'
    pos = graphviz_layout(bipartite_networkx_graph, prog=layout_prog)
    nx.draw_networkx_nodes(bipartite_networkx_graph, pos, nodelist=metabolite_nodes,
                           node_color=metabolite_node_colors, node_shape='o', alpha=0.6, node_size=2000)
    nx.draw_networkx_nodes(bipartite_networkx_graph, pos, nodelist=function_nodes, node_shape='d', alpha=0.6, node_size=1000)
    nx.draw_networkx_labels(bipartite_networkx_graph, pos, labels=pathway_index_labels, font_size=14)
    nx.draw_networkx_edges(bipartite_networkx_graph, pos, arrows=True, width=3, node_size=1000, connectionstyle="arc3,rad=0.1")

    # Relabel nodes.
    bipartite_networkx_graph = nx.relabel_nodes(bipartite_networkx_graph, pathway_index_labels)
    graphml_output_file = graph_output_file.replace('.png', '.graphml')
    nx.write_graphml(bipartite_networkx_graph, graphml_output_file)

    plt.axis('off')
    plt.savefig(graph_output_file, bbox_inches='tight')
    plt.clf()


def create_visualisation(bigecyhmm_output, output_folder, esmecata_output_folder=None, abundance_file_path=None, group_file=None, metabolite_measure=None,
                         bigecyhmm_run_database=None, background_path_donut_plot=None):
    """Create visualisation plots from esmecata, bigecyhmm output folders

    Args:
        bigecyhmm_output (str): path to bigecyhmm output folder.
        output_folder (str): path to the output folder where files will be created.
        esmecata_output_folder (str): path to esmecata output folder.
        abundance_file_path (str): path to abundance file.
        group_file_path (str): path to group file.
        metabolite_measure (str): path to metaboltie measure file indicating the abundance of metabolites in samples.
        bigecyhmm_run_database (sttr): path to bigecyhmm run internal database (only when used with bigecyhmm_custom).
        background_path_donut_plot (str): path to background figure for donut plot.
    """
    start_time = time.time()

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    output_folder_occurrence = os.path.join(output_folder, 'function_occurrence')
    if not os.path.exists(output_folder_occurrence):
        os.mkdir(output_folder_occurrence)

    if abundance_file_path is not None:
        sample_abundance, sample_tot_abundance = read_measures_file(abundance_file_path)
        output_folder_abundance = os.path.join(output_folder, 'function_abundance')
        if not os.path.exists(output_folder_abundance):
            os.mkdir(output_folder_abundance)
    else:
        sample_abundance = None
        sample_tot_abundance = None
        output_folder_abundance = None

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
    else:
        tax_id_names_observation_names = None
        observation_names_tax_ranks = None

    # Group file must be used with abundance file argument.
    if group_file is not None and abundance_file_path is None:
        logger.info("Option --group-file/group_file must be used with option --abundance-file/abundance_file_path.")
        sys.exit(1)
    # Check that all samples from abundance file are in group_file.
    if group_file is not None:
        mapping_df = pd.read_csv(group_file, sep='\t', dtype=str)
        samples_from_group_file = set(mapping_df['sample'].tolist())
        samples_from_abundance_file = set(sample_abundance.keys())
        if samples_from_group_file != set(sample_abundance.keys()):
            missing_samples = samples_from_abundance_file - samples_from_group_file
            logger.info("Missing samples from abundance file in group file: {0}.".format(' '.join(missing_samples)))
            sys.exit(1)

    logger.info("## Compute function occurrences and create visualisation.")
    logger.info("  -> Read bigecyhmm cycle output files.")
    bigecyhmm_pathway_presence_file = os.path.join(bigecyhmm_output, 'pathway_presence.tsv')
    cycle_occurrence_organisms, all_studied_organisms = compute_bigecyhmm_functions_occurrence(bigecyhmm_pathway_presence_file, tax_id_names_observation_names)

    if abundance_file_path is not None:
        if abundance_file_path.endswith('.tsv'):
            delimiter = '\t'
        elif abundance_file_path.endswith('.csv'):
            delimiter = ','
        abundance_data_df = pd.read_csv(abundance_file_path, sep=delimiter)
        all_studied_organisms = abundance_data_df.index.tolist()

    df_cycle_occurrence_organisms = pd.DataFrame(cycle_occurrence_organisms)
    df_cycle_occurrence_organisms.index.name = 'function'
    df_cycle_occurrence_organisms.fillna(0, inplace=True)
    organism_pathway_presence_file = os.path.join(output_folder_occurrence, 'pathway_presence_in_organism.tsv')
    df_cycle_occurrence_organisms.to_csv(organism_pathway_presence_file, sep='\t')

    # Compute the relative occurrence of functions by dividing the number of organisms having it by the total number of organisms in the community.
    cycle_occurrence_community = []
    for function in cycle_occurrence_organisms:
        cycle_occurrence_community.append([function, len(cycle_occurrence_organisms[function])/len(all_studied_organisms)])

    cycle_occurrence_community_df = pd.DataFrame(cycle_occurrence_community, columns=['name', 'ratio'])
    cycle_occurrence_community_df.set_index('name', inplace=True)
    cycle_occurrence_community_df.to_csv(os.path.join(output_folder_occurrence, 'cycle_occurrence.tsv'), sep='\t')
    cycle_occurrences = cycle_occurrence_community_df['ratio'].to_dict()

    all_cycles = pd.read_csv(bigecyhmm_pathway_presence_file, sep='\t')['function'].tolist()
    all_bigecyhmm_template_cycles = pd.read_csv(PATHWAY_TEMPLATE_FILE, sep='\t')['Pathways'].tolist()

    custom_central_hydrogen_template_df = pd.read_csv(CUSTOM_HYDROGEN_TABLE, sep='\t')
    custom_central_hydrogen_template_df = custom_central_hydrogen_template_df[custom_central_hydrogen_template_df['Type']=='FUNCTION']
    all_custom_central_hydrogen_template_cycles = custom_central_hydrogen_template_df['ID'].tolist()

    if set(all_bigecyhmm_template_cycles).issubset(set(all_cycles)):
        logger.info("  -> Create diagrams.")
        diagram_data = {}
        for cycle_name in all_cycles:
            if cycle_name in cycle_occurrences:
                diagram_data[cycle_name] = (round(sum(cycle_occurrence_organisms[cycle_name].values()), 1), round(cycle_occurrences[cycle_name]*100, 1))
            else:
                diagram_data[cycle_name] = (0, 0)

        carbon_cycle_file = os.path.join(output_folder_occurrence, 'diagram_carbon_cycle.png')
        create_carbon_cycle(diagram_data, carbon_cycle_file, 'Occurrence', 'Percentage')

        nitrogen_cycle_file = os.path.join(output_folder_occurrence, 'diagram_nitrogen_cycle.png')
        create_nitrogen_cycle(diagram_data, nitrogen_cycle_file, 'Occurrence', 'Percentage')

        sulfur_cycle_file = os.path.join(output_folder_occurrence, 'diagram_sulfur_cycle.png')
        create_sulfur_cycle(diagram_data, sulfur_cycle_file, 'Occurrence', 'Percentage')

        other_cycle_file = os.path.join(output_folder_occurrence, 'diagram_other_cycle.png')
        create_other_cycle(diagram_data, other_cycle_file, 'Occurrence', 'Percentage')

    logger.info("  -> Read bigecyhmm functions output files.")
    bigecyhmm_function_presence_file = os.path.join(bigecyhmm_output, 'function_presence.tsv')
    function_occurrence_organisms, all_studied_organisms = compute_bigecyhmm_functions_occurrence(bigecyhmm_function_presence_file, tax_id_names_observation_names)

    df_function_occurrence_organisms = pd.DataFrame(function_occurrence_organisms)
    df_function_occurrence_organisms.index.name = 'function'
    df_function_occurrence_organisms.fillna(0, inplace=True)
    organism_function_presence_file = os.path.join(output_folder_occurrence, 'function_occurrence_in_organism.tsv')
    df_function_occurrence_organisms.to_csv(organism_function_presence_file, sep='\t')

    # Compute the relative abundance of organisms by dividing for a function the number of organisms having it by the total number of organisms in the community.
    function_occurrence_community = []
    for index in function_occurrence_organisms:
        function_occurrence_community.append([index, len(function_occurrence_organisms[index])/len(all_studied_organisms)])

    function_occurrence_community_df = pd.DataFrame(function_occurrence_community, columns=['name', 'ratio'])
    function_occurrence_community_df.set_index('name', inplace=True)
    function_occurrence_community_df.to_csv(os.path.join(output_folder_occurrence, 'function_occurrence.tsv'), sep='\t')

    logger.info("  -> Create heatmap.")
    output_heatmap_filepath = os.path.join(output_folder_occurrence, 'heatmap_occurrence.png')
    create_heatmap_functions(function_occurrence_community_df, output_heatmap_filepath)

    if abundance_file_path is not None:
        logger.info("## Compute function abundances and create visualisation.")
        logger.info("  -> Read abundance file.")
        sample_abundance, sample_tot_abundance = read_measures_file(abundance_file_path)
        output_folder_abundance = os.path.join(output_folder, 'function_abundance')
        if not os.path.exists(output_folder_abundance):
            os.mkdir(output_folder_abundance)

        if observation_names_tax_ranks is not None:
            generate_barplot_esmecata_taxon_abundance(sample_abundance, observation_names_tax_ranks, sample_tot_abundance, output_folder_abundance)
            specific_function_folder = os.path.join(output_folder_abundance, 'cycle_taxa_abundance')
            if not os.path.exists(specific_function_folder):
                os.mkdir(specific_function_folder)
            taxon_function_heatmap(df_cycle_occurrence_organisms, proteome_tax_id_file, sample_abundance, specific_function_folder)

        logger.info("  -> Read bigecyhmm cycle output files.")
        bigecyhmm_pathway_presence_file = os.path.join(bigecyhmm_output, 'pathway_presence.tsv')
        cycle_abundance_samples, cycle_relative_abundance_samples, cycle_participation_samples = compute_bigecyhmm_functions_abundance(bigecyhmm_pathway_presence_file, sample_abundance, sample_tot_abundance, tax_id_names_observation_names)

        cycle_raw_abundance_samples_df = pd.DataFrame(cycle_abundance_samples)
        cycle_raw_abundance_samples_df.index.name = 'name'
        cycle_raw_abundance_samples_df.sort_index(inplace=True)
        cycle_raw_abundance_samples_df.to_csv(os.path.join(output_folder_abundance, 'cycle_abundance_sample_raw.tsv'), sep='\t')

        cycle_relative_abundance_samples_df = pd.DataFrame(cycle_relative_abundance_samples)
        cycle_relative_abundance_samples_df.index.name = 'name'
        cycle_relative_abundance_samples_df.sort_index(inplace=True)
        cycle_abundance_sample_filepath = os.path.join(output_folder_abundance, 'cycle_abundance_sample.tsv')
        if set(all_custom_central_hydrogen_template_cycles).issubset(set(all_cycles)):
            cycle_relative_abundance_samples_df = cycle_relative_abundance_samples_df.reindex(all_custom_central_hydrogen_template_cycles)
        cycle_relative_abundance_samples_df.to_csv(cycle_abundance_sample_filepath, sep='\t')

        logger.info("  -> Compute function abundance participation in each sample.")
        output_folder_cycle_participation = os.path.join(output_folder_abundance, 'cycle_participation')
        if not os.path.exists(output_folder_cycle_participation):
            os.mkdir(output_folder_cycle_participation)

        for sample in cycle_participation_samples:
            data_cycle_participation = []
            index_organism_names = []
            for organism in cycle_participation_samples[sample]:
                data_cycle_participation.append([*[cycle_participation_samples[sample][organism][function_name] if function_name in cycle_participation_samples[sample][organism] else 0 for function_name in all_cycles]])
                index_organism_names.append(organism)
            data_cycle_participation_df = pd.DataFrame(data_cycle_participation, index=index_organism_names, columns=all_cycles)
            data_cycle_participation_df.index.name = 'organism'
            data_cycle_participation_df.to_csv(os.path.join(output_folder_cycle_participation, sample+'.tsv'), sep='\t')

        logger.info("  -> Create polarplot.")
        cycle_relative_abundance_samples_df.reset_index(inplace=True)
        melted_cycle_relative_abundance_samples_df = pd.melt(cycle_relative_abundance_samples_df, id_vars='name', value_vars=cycle_relative_abundance_samples_df.columns.tolist())
        melted_cycle_relative_abundance_samples_df.columns = ['name', 'sample', 'ratio']
        melted_cycle_relative_abundance_samples_df.to_csv(os.path.join(output_folder_abundance, 'cycle_abundance_sample_melted.tsv'), sep='\t')

        output_folder_polar_plot = os.path.join(output_folder_abundance, 'polar_plot_abundance')
        if not os.path.exists(output_folder_polar_plot):
            os.mkdir(output_folder_polar_plot)
        for sample in melted_cycle_relative_abundance_samples_df['sample'].unique():
            output_polar_plot = os.path.join(output_folder_polar_plot, 'polar_plot_abundance_sample_'+sample+'.png')
            sample_melted_cycle_relative_abundance_samples_df = melted_cycle_relative_abundance_samples_df[melted_cycle_relative_abundance_samples_df['sample']==sample]
            create_polar_plot(sample_melted_cycle_relative_abundance_samples_df, output_polar_plot)

        logger.info("  -> Create diagrams.")
        output_folder_cycle_diagram = os.path.join(output_folder_abundance, 'cycle_diagrams_abundance')
        if not os.path.exists(output_folder_cycle_diagram):
            os.mkdir(output_folder_cycle_diagram)

        if set(all_bigecyhmm_template_cycles).issubset(set(all_cycles)):
            for sample in cycle_relative_abundance_samples:
                diagram_data = {}
                for cycle_name in all_cycles:
                    if cycle_name in cycle_relative_abundance_samples[sample]:
                        diagram_data[cycle_name] = (round(cycle_abundance_samples[sample][cycle_name], 1), round(cycle_relative_abundance_samples[sample][cycle_name]*100, 1))
                    else:
                        diagram_data[cycle_name] = (0, 0)

                carbon_cycle_file = os.path.join(output_folder_cycle_diagram, sample + '_carbon_cycle.png')
                create_carbon_cycle(diagram_data, carbon_cycle_file, 'Abundance', 'Percentage')

                nitrogen_cycle_file = os.path.join(output_folder_cycle_diagram, sample + '_nitrogen_cycle.png')
                create_nitrogen_cycle(diagram_data, nitrogen_cycle_file, 'Abundance', 'Percentage')

                sulfur_cycle_file = os.path.join(output_folder_cycle_diagram, sample + '_sulfur_cycle.png')
                create_sulfur_cycle(diagram_data, sulfur_cycle_file, 'Abundance', 'Percentage')

                other_cycle_file = os.path.join(output_folder_cycle_diagram, sample + '_other_cycle.png')
                create_other_cycle(diagram_data, other_cycle_file, 'Abundance', 'Percentage')

                phosphorus_cycle_file = os.path.join(output_folder_cycle_diagram, sample + '_phosphorus_cycle.png')
                create_phosphorus_cycle(diagram_data, phosphorus_cycle_file, 'Abundance', 'Percentage')

        if background_path_donut_plot is None and set(all_bigecyhmm_template_cycles).issubset(set(all_cycles)):
            background_path_donut_plot = TEMPLATE_BACKGROUND_BIGECYHMM
        if background_path_donut_plot is None and set(all_custom_central_hydrogen_template_cycles).issubset(set(all_cycles)):
            background_path_donut_plot = TEMPLATE_CUSTOM_CENTRAL_HYDROGEN
        if background_path_donut_plot is not None:
            output_folder_plots = os.path.join(output_folder_abundance, 'plots')
            if not os.path.exists(output_folder_plots):
                os.mkdir(output_folder_plots)
            group_stats_file = os.path.join(output_folder_plots, 'group_stats.tsv')
            cleaned_data_file = os.path.join(output_folder_plots, 'cleaned_data.tsv')
            group_medians_donut_file = os.path.join(output_folder_plots, 'group_medians_donut.png')
            group_stats_table_file = os.path.join(output_folder_plots, 'group_stats_table.png')
            statNut_run(input_tsv=cycle_abundance_sample_filepath, sample_groups_tsv=group_file, background_path=background_path_donut_plot,
                        output_stats_csv=group_stats_file, output_cleaned_csv=cleaned_data_file, output_donut_png=group_medians_donut_file, output_table_png=group_stats_table_file)
            """Function to: 
            - resolve sample groups based on an input-tsv 
            - calculate stats on groups (kruskal-wallis with Benjamini-hochberg correction) 
            - make a table for the functions/ stats 
            - make a donut-plot showing a "drawn" network, group medians, max+min of each group for each function
            - output a combined figure with donut & table"""

        # Generate a donut graph visualisation from input graph file from custom_db. 
        if bigecyhmm_run_database is not None:
            output_folder_graph_plots = os.path.join(output_folder_abundance, 'graph_plots')
            if not os.path.exists(output_folder_graph_plots):
                os.mkdir(output_folder_graph_plots)
            # Generate background figure from graphml in database folder.
            background_graph_output_file = os.path.join(output_folder_graph_plots, 'graph_reference.png')
            generate_graph_figure(bigecyhmm_run_database, background_graph_output_file)

            # Generate donut plot using generated background image from graph.
            group_stats_file = os.path.join(output_folder_graph_plots, 'group_stats.tsv')
            cleaned_data_file = os.path.join(output_folder_graph_plots, 'cleaned_data.tsv')
            group_medians_donut_file = os.path.join(output_folder_graph_plots, 'group_medians_donut.png')
            group_stats_table_file = os.path.join(output_folder_graph_plots, 'group_stats_table.png')
            statNut_run(input_tsv=cycle_abundance_sample_filepath, sample_groups_tsv=group_file, background_path=background_graph_output_file,
                        output_stats_csv=group_stats_file, output_cleaned_csv=cleaned_data_file, output_donut_png=group_medians_donut_file, output_table_png=group_stats_table_file)

        bubble_plot_output_file = os.path.join(output_folder_abundance, 'cycle_pathways_bubble_plot.png')
        generate_bubble_plot(melted_cycle_relative_abundance_samples_df, bubble_plot_output_file, group_file)

        logger.info("  -> Read bigecyhmm function output files.")
        bigecyhmm_function_presence_file = os.path.join(bigecyhmm_output, 'function_presence.tsv')
        function_abundance_samples, function_relative_abundance_samples, function_participation_samples = compute_bigecyhmm_functions_abundance(bigecyhmm_function_presence_file, sample_abundance, sample_tot_abundance, tax_id_names_observation_names)

        function_relative_abundance_samples_df = pd.DataFrame(function_relative_abundance_samples)
        function_relative_abundance_samples_df.index.name = 'name'
        function_relative_abundance_samples_df.sort_index(inplace=True)
        function_relative_abundance_samples_df.to_csv(os.path.join(output_folder_abundance, 'function_abundance_sample.tsv'), sep='\t')

        function_absolute_abundance_samples_df = pd.DataFrame(function_abundance_samples)
        function_absolute_abundance_samples_df.index.name = 'name'
        function_absolute_abundance_samples_df.sort_index(inplace=True)
        function_absolute_abundance_samples_df.to_csv(os.path.join(output_folder_abundance, 'function_abundance_sample_absolute.tsv'), sep='\t')

        logger.info("  -> Compute function abundance participation in each sample.")
        output_folder_function_participation = os.path.join(output_folder_abundance, 'function_participation')
        if not os.path.exists(output_folder_function_participation):
            os.mkdir(output_folder_function_participation)

        all_functions = function_relative_abundance_samples_df.index.tolist()

        for sample in function_participation_samples:
            data_function_participation = []
            index_organism_names = []
            for organism in function_participation_samples[sample]:
                data_function_participation.append([*[function_participation_samples[sample][organism][function_name] if function_name in function_participation_samples[sample][organism] else 0 for function_name in all_functions]])
                index_organism_names.append(organism)
            data_function_participation_df = pd.DataFrame(data_function_participation, index=index_organism_names, columns=all_functions)
            data_function_participation_df.index.name = 'organism'
            data_function_participation_df.to_csv(os.path.join(output_folder_function_participation, sample+'.tsv'), sep='\t')

        logger.info("  -> Create heatmap.")
        output_heatmap_filepath = os.path.join(output_folder_abundance, 'heatmap_abundance_samples.png')
        function_relative_abundance_samples_df = function_relative_abundance_samples_df.replace(0, np.nan)
        create_heatmap_functions(function_relative_abundance_samples_df, output_heatmap_filepath)
        output_heatmap_filepath = os.path.join(output_folder_abundance, 'heatmap_abundance_samples.svg')
        create_heatmap_functions(function_relative_abundance_samples_df, output_heatmap_filepath)

        function_relative_abundance_samples_df.reset_index(inplace=True)
        melted_function_relative_abundance_samples_df = pd.melt(function_relative_abundance_samples_df, id_vars='name', value_vars=function_relative_abundance_samples_df.columns.tolist())
        melted_function_relative_abundance_samples_df.columns = ['name', 'sample', 'ratio']

        # Create plot for hydrogenases.
        #logger.info("  -> Create barplot of hydrogenases.")
        #output_barplot_hydrogenase_filepath = os.path.join(output_folder_abundance, 'barplot_abundance_hydrogenase.png')
        #function_categories = get_function_categories()
        #visualise_barplot_category('Hydrogenases', function_categories, melted_function_relative_abundance_samples_df, output_barplot_hydrogenase_filepath)

        # Create HMM functional profiles.
        logger.info("  -> Create HMM functional profiles.")
        create_ko_functional_profile(bigecyhmm_output, sample_abundance, output_folder_abundance, tax_id_names_observation_names)
    """
    if abundance_file_path is not None:
        output_file_name = os.path.join(output_folder_abundance, 'barplot_gene_function.png')
        visualise_barplot_category('Fermentation', gene_categories, df_seaborn_abundance, output_file_name)
        output_file_name = os.path.join(output_folder_abundance, 'barplot_gene_function_2.png')
        kept_names = [name for name in df_seaborn_abundance['name'] if 'Wood' in name]
        df_seaborn_abundance = df_seaborn_abundance[df_seaborn_abundance['name'].isin(kept_names)]
        visualise_barplot_category('Carbon fixation', gene_categories, df_seaborn_abundance, output_file_name)
        output_heatmap_filepath = os.path.join(output_folder_abundance, 'heatmap_abundance_samples.png')
        create_heatmap_functions(df_heatmap_abundance_samples, output_heatmap_filepath)
    """
    # If there is a graph in bigecyhmm folder, run graph analyses.
    graph_file = os.path.join(bigecyhmm_output, 'input_graph.graphml')
    if os.path.exists(graph_file):
        if abundance_file_path is None:
            cycle_raw_abundance_samples_df = None
        add_abundance_and_measure_to_graph(graph_file, output_folder, metabolite_measure, cycle_raw_abundance_samples_df)

    duration = time.time() - start_time
    metadata_json = {}
    metadata_json['tool_dependencies'] = {}
    metadata_json['tool_dependencies']['python_package'] = {}
    metadata_json['tool_dependencies']['python_package']['Python_version'] = sys.version
    metadata_json['tool_dependencies']['python_package']['bigecyhmm'] = bigecyhmm_version
    metadata_json['tool_dependencies']['python_package']['pandas'] = pandas_version
    metadata_json['tool_dependencies']['python_package']['matplotlib'] = matplotlib_version
    metadata_json['tool_dependencies']['python_package']['seaborn'] = seaborn_version

    metadata_json['input_parameters'] = {'esmecata_output_folder': esmecata_output_folder, 'bigecyhmm_output': bigecyhmm_output, 'output_folder': output_folder,
                                         'abundance_file_path': abundance_file_path}
    metadata_json['duration'] = duration

    metadata_file = os.path.join(output_folder, 'bigecyhmm_visualisation_metadata.json')
    with open(metadata_file, 'w') as ouput_file:
        json.dump(metadata_json, ouput_file, indent=4)


def create_visualisation_from_ko_file(ko_abundance_file, output_folder, group_file=None):
    """Create visualisation plots from abundance file with KEGG Orthologs.

    Args:
        ko_abundance_file (str): path to ko abundance file.
        output_folder (str): path to the output folder where files will be created.
    """
    start_time = time.time()

    pathway_hmms, pathway_expression, sorted_pathways = get_diagram_pathways_hmms(PATHWAY_TEMPLATE_FILE)
    df_ko_abundance = pd.read_csv(ko_abundance_file, sep='\t', index_col=0)
    kegg_ortholog_abundance_samples = df_ko_abundance.to_dict()

    all_kos = df_ko_abundance.index.tolist()

    sample_data_pathway = {}
    for sample in kegg_ortholog_abundance_samples:
        ko_sample_list = [ko+'.hmm' for ko in kegg_ortholog_abundance_samples[sample] if  kegg_ortholog_abundance_samples[sample][ko] > 0]
        for pathway in sorted_pathways:
            pathway_checks = []
            hmms_in_org = []
            # For AND group of HMMs in pathway, check if their corresponding HMMs are there.
            pathway_functions = []
            for hmm_combination in pathway_hmms[pathway]:
                pathway_check = False
                negative_hmms = [hmm.replace('NO|', '') for hmm in hmm_combination if 'NO' in hmm]
                # Check if there are negative HMMs in pathway string.
                if len(negative_hmms) > 0:
                    hmm_combination = [hmm.replace('NO|', '') for hmm in hmm_combination]
                intersection_hmms = set(hmm_combination).intersection(ko_sample_list)
                intersection_negative_hmms = set(negative_hmms).intersection(ko_sample_list)

                # First check if all combination corresponds to negative HMMs.
                if len(hmm_combination) == len(negative_hmms):
                    # If all HMMs of the combination are negative ones and are not present in organism, then this combination is checked.
                    if len(intersection_hmms) == 0:
                        pathway_check = True
                else:
                    # If all HMMs of the combination present in the organism are not negative HMMs, then this combination of HMMs is checked.
                    if len(intersection_hmms) > 0:
                        if len(intersection_negative_hmms) == 0:
                            pathway_check = True
                        # But if there are at least one negative HMMs, it is not checked.
                        elif len(intersection_negative_hmms) > 0:
                            pathway_check = False
                pathway_checks.append(pathway_check)
                positive_hmms = list(set(intersection_hmms) - set(intersection_negative_hmms))
                found_hmms = list(positive_hmms + ['NO|' + hmm for hmm in intersection_negative_hmms])
                if len(found_hmms) > 0:
                    hmms_in_org.append(', '.join(found_hmms))
                    pathway_functions.extend([kegg_ortholog_abundance_samples[sample][ko.replace('.hmm', '')] for ko in found_hmms if kegg_ortholog_abundance_samples[sample][ko.replace('.hmm', '')] > 0])

            # If all combination HMMs have been checked, keep the pathway.
            if all(pathway_checks) is True:
                if sample not in sample_data_pathway:
                    sample_data_pathway[sample] = {}
                func_abundance = max(pathway_functions)
                sample_data_pathway[sample][pathway] = func_abundance

            else:
                if sample not in sample_data_pathway:
                    sample_data_pathway[sample] = {}
                sample_data_pathway[sample][pathway] = 0

    sample_pathway_df = pd.DataFrame(sample_data_pathway)
    output_folder_cycle_tsv = os.path.join(output_folder, 'pathway_presence_abundance.tsv')
    sample_pathway_df.index.name = 'pathway'
    sample_pathway_df.to_csv(output_folder_cycle_tsv, sep='\t')

    output_folder_cycle_diagram = os.path.join(output_folder, 'diagram_visualisation')
    if not os.path.exists(output_folder_cycle_diagram):
        os.mkdir(output_folder_cycle_diagram)

    for sample in sample_data_pathway:
        diagram_data = {}
        for cycle_name in sample_data_pathway[sample]:
            if cycle_name in sample_data_pathway[sample]:
                diagram_data[cycle_name] = (round(sample_data_pathway[sample][cycle_name], 1), 0)
            else:
                diagram_data[cycle_name] = (0, 0)

        carbon_cycle_file = os.path.join(output_folder_cycle_diagram, sample + '_carbon_cycle.png')
        create_carbon_cycle(diagram_data, carbon_cycle_file)

        nitrogen_cycle_file = os.path.join(output_folder_cycle_diagram, sample + '_nitrogen_cycle.png')
        create_nitrogen_cycle(diagram_data, nitrogen_cycle_file)

        sulfur_cycle_file = os.path.join(output_folder_cycle_diagram, sample + '_sulfur_cycle.png')
        create_sulfur_cycle(diagram_data, sulfur_cycle_file)

        other_cycle_file = os.path.join(output_folder_cycle_diagram, sample + '_other_cycle.png')
        create_other_cycle(diagram_data, other_cycle_file)

        phosphorus_cycle_file = os.path.join(output_folder_cycle_diagram, sample + '_phosphorus_cycle.png')
        create_phosphorus_cycle(diagram_data, phosphorus_cycle_file)

    duration = time.time() - start_time
    metadata_json = {}
    metadata_json['tool_dependencies'] = {}
    metadata_json['tool_dependencies']['python_package'] = {}
    metadata_json['tool_dependencies']['python_package']['Python_version'] = sys.version
    metadata_json['tool_dependencies']['python_package']['bigecyhmm'] = bigecyhmm_version
    metadata_json['tool_dependencies']['python_package']['pandas'] = pandas_version
    metadata_json['tool_dependencies']['python_package']['matplotlib'] = matplotlib_version
    metadata_json['tool_dependencies']['python_package']['seaborn'] = seaborn_version

    metadata_json['input_parameters'] = {'ko_abundance_file': ko_abundance_file, 'output_folder': output_folder}
    metadata_json['duration'] = duration

    metadata_file = os.path.join(output_folder, 'bigecyhmm_visualisation_metadata.json')
    with open(metadata_file, 'w') as ouput_file:
        json.dump(metadata_json, ouput_file, indent=4)

    return sample_data_pathway


def visualisation_input_handler(bigecyhmm_output, output_folder, esmecata_output_folder=None, abundance_file_path=None,
                                group_file=None, metabolite_measure=None, background_path_donut_plot=None):
    """Create visualisation plots from esmecata, bigecyhmm output folders

    Args:
        bigecyhmm_output (str): path to bigecyhmm output folder.
        output_folder (str): path to the output folder where files will be created.
        esmecata_output_folder (str): path to esmecata output folder.
        abundance_file_path (str): path to abundance file.
        group_file_path (str): path to group file.
        metabolite_measure (str): path to metaboltie measure file indicating the abundance of metabolites in samples.
        background_path_donut_plot (str): path to background figure for donut plot.
    """
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    bigecyhmm_pathway_presence_file = os.path.join(bigecyhmm_output, 'pathway_presence.tsv')
    # Output files are directly in input folder, run bigecyhmm visualisation on it.
    if os.path.exists(bigecyhmm_pathway_presence_file):
        logger.info("|bigecyhmm|visualisation| Launch analyses on {0}.".format(bigecyhmm_output))
        create_visualisation(bigecyhmm_output, output_folder, esmecata_output_folder, abundance_file_path, group_file, metabolite_measure)
    # If there is no such files but if there are subfolders, check if bigecyhmm output files are not in subfolders. 
    else:
        bigecyhmm_run_internal_database_path = os.path.join(bigecyhmm_output, 'database')
        for bigecyhmm_input_folder in os.listdir(bigecyhmm_output):
            bigecyhmm_input_folder_path = os.path.join(bigecyhmm_output, bigecyhmm_input_folder)
            bigecyhmm_pathway_presence_file = os.path.join(bigecyhmm_input_folder_path, 'pathway_presence.tsv')
            if os.path.isdir(bigecyhmm_input_folder_path) and os.path.exists(bigecyhmm_pathway_presence_file):
                # If yes, then run each time bigecyhmm visualisation on these different subfolders.
                logger.info("|bigecyhmm|visualisation| Found one subfolder {0} from {1}, launch analysis on it.".format(bigecyhmm_input_folder_path, bigecyhmm_output))
                subfolder_output_folder = os.path.join(output_folder, bigecyhmm_input_folder)
                if os.path.exists(bigecyhmm_run_internal_database_path) is False:
                    bigecyhmm_run_internal_database_path = None
                create_visualisation(bigecyhmm_input_folder_path, subfolder_output_folder, esmecata_output_folder, abundance_file_path, group_file, metabolite_measure,
                                     bigecyhmm_run_internal_database_path, background_path_donut_plot)


def main():
    start_time = time.time()

    parser = argparse.ArgumentParser(
        'bigecyhmm_visualisation',
        description=MESSAGE + ' For specific help on each subcommand use: bigecyhmm_visualisation {cmd} --help',
        epilog=REQUIRES
    )
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s ' + bigecyhmm_version + '\n')

    parent_parser_esmecata = argparse.ArgumentParser(add_help=False)
    parent_parser_esmecata.add_argument(
        '--esmecata',
        dest='esmecata',
        required=True,
        help='EsMeCaTa output folder for the input file.',
        metavar='INPUT_FOLDER')

    parent_parser_bigecyhmm = argparse.ArgumentParser(add_help=False)
    parent_parser_bigecyhmm.add_argument(
        '--bigecyhmm',
        dest='bigecyhmm',
        required=True,
        help='Bigecyhmm output folder for the input file.',
        metavar='INPUT_FOLDER')

    parent_parser_abundance_file = argparse.ArgumentParser(add_help=False)
    parent_parser_abundance_file.add_argument(
        '--abundance-file',
        dest='abundance_file',
        required=False,
        help='Abundance file indicating the abundance for each organisms.',
        metavar='INPUT_FILE')

    parent_parser_ko_file = argparse.ArgumentParser(add_help=False)
    parent_parser_ko_file.add_argument(
        '--ko',
        dest='ko_file',
        required=True,
        help='Abundance file indicating the abundance of KO in different sampels (coming form tabigecy, picrust, ...).',
        metavar='INPUT_FILE')

    parent_parser_output_folder = argparse.ArgumentParser(add_help=False)
    parent_parser_output_folder.add_argument(
        '-o',
        '--output',
        dest='output',
        required=True,
        help='Output directory path.',
        metavar='OUPUT_DIR')

    parent_parser_group_file = argparse.ArgumentParser(add_help=False)
    parent_parser_group_file.add_argument(
        '--group-file',
        dest='group_file',
        required=False,
        help='Group file associating samples with group.',
        metavar='INPUT_FILE')

    parent_parser_measure_file = argparse.ArgumentParser(add_help=False)
    parent_parser_measure_file.add_argument(
        '--measure-file',
        dest='measure_file',
        required=False,
        help='Measure file indicating the abundance for each metabolties of the graph in different samples.',
        metavar='INPUT_FILE',
        default=None)

    parent_parser_background_figure_file = argparse.ArgumentParser(add_help=False)
    parent_parser_background_figure_file.add_argument(
        '--background-file',
        dest='background_file',
        required=False,
        help='Background figure file for donut plot with group.',
        metavar='INPUT_FILE')

    # subparsers
    subparsers = parser.add_subparsers(
        title='subcommands',
        description='valid subcommands:',
        dest='cmd')

    esmecata_parser = subparsers.add_parser(
        'esmecata',
        help='Create visualisation from runs of EsMeCaTa and bigecyhmm.',
        parents=[
            parent_parser_esmecata, parent_parser_bigecyhmm, parent_parser_abundance_file,
            parent_parser_output_folder, parent_parser_group_file, parent_parser_measure_file,
            parent_parser_background_figure_file
            ],
        allow_abbrev=False)
    genomes_parser = subparsers.add_parser(
        'genomes',
        help='Creates visualisation from runs of bigecyhmm on genomes.',
        parents=[
            parent_parser_bigecyhmm, parent_parser_abundance_file,
            parent_parser_output_folder, parent_parser_group_file,
            parent_parser_measure_file
            ],
        allow_abbrev=False)
    ko_parser = subparsers.add_parser(
        'ko',
        help='Creates visualisation from a table containing the abundances of HMM (especially KO) for different samples.',
        parents=[
            parent_parser_ko_file, parent_parser_output_folder, parent_parser_group_file
            ],
        allow_abbrev=False)

    args = parser.parse_args()

    # If no argument print the help.
    if len(sys.argv) == 1 or len(sys.argv) == 0:
        parser.print_help()
        sys.exit(1)

    is_valid_dir(args.output)

    # add logger in file
    formatter = logging.Formatter('%(message)s')
    log_file_path = os.path.join(args.output, f'bigecyhmm_visualisation.log')
    file_handler = logging.FileHandler(log_file_path, 'w+')
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    # set up the default console logger
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    logger.info("--- Create visualisation ---")

    if args.cmd in ['esmecata', 'genomes']:
        if args.abundance_file == 'false':
            abundance_file = None
        else:
            abundance_file = args.abundance_file

        if args.group_file == 'false':
            group_file = None
        else:
            group_file = args.group_file

    if args.cmd in ['esmecata']:
        visualisation_input_handler(args.bigecyhmm, args.output, esmecata_output_folder=args.esmecata, abundance_file_path=abundance_file, group_file=group_file, metabolite_measure=args.measure_file,
                                    background_path_donut_plot=args.background_file)
    elif args.cmd in ['genomes']:
        visualisation_input_handler(args.bigecyhmm, args.output, abundance_file_path=abundance_file, group_file=group_file, metabolite_measure=args.measure_file)
    elif args.cmd in ['ko']:
        create_visualisation_from_ko_file(args.ko_file, args.output)

    duration = time.time() - start_time
    logger.info("--- Total runtime %.2f seconds ---" % (duration))
    logger.warning(f'--- Logs written in {log_file_path} ---')
