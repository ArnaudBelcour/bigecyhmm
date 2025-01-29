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

import pandas as pd
from pandas import __version__ as pandas_version
import seaborn as sns
from seaborn import __version__ as seaborn_version
import matplotlib.pyplot as plt
from matplotlib import __version__ as matplotlib_version
from kaleido import __version__ as kaleido_version

from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.express as px
from plotly import __version__ as plotly_version

import argparse
import logging
import json
import math
import os
import sys
import time

from bigecyhmm import __version__ as bigecyhmm_version
from bigecyhmm.utils import is_valid_dir
from bigecyhmm.diagram_cycles import create_carbon_cycle, create_nitrogen_cycle, create_sulfur_cycle, create_other_cycle

MESSAGE = '''
Create figures from bigecyhmm and esmecata outputs.
'''
REQUIRES = '''
Requires seaborn, pandas, plotly and kaleido.
'''

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

ROOT = os.path.dirname(__file__)
HMM_TEMPLATE_FILE = os.path.join(ROOT, 'hmm_databases', 'hmm_table_template.tsv')


def read_abundance_file(abundance_file_path):
    """Read abundance file for samples. Expect a tsv or csv files with organisms as rows, samples as columns and abundance as values.

    Args:
        abundance_file_path (str): path to abundance file

    Returns:
        sample_abundance (dict): for each sample, subdict with the abundance of the different organisms.
        sample_tot_abundance (dict): for each sample, the total abundance of all organisms in the sample.
    """
    if abundance_file_path.endswith('.tsv'):
        input_data_df = pd.read_csv(abundance_file_path, sep='\t')
    elif abundance_file_path.endswith('.csv'):
        input_data_df = pd.read_csv(abundance_file_path)
    input_data_df.set_index('observation_name', inplace=True)

    sample_abundance = {}
    sample_tot_abundance = {}
    for sample_name in input_data_df.columns:
        sample_abundance[sample_name] = input_data_df[sample_name].to_dict()
        tot_abundance = input_data_df[sample_name].sum()
        sample_tot_abundance[sample_name] = tot_abundance

    return sample_abundance, sample_tot_abundance


def read_esmecata_proteome_file(proteome_tax_id_file):
    """Read esmecata proteome file to extract associated betwenn organism name and tax_id_name.

    Args:
        proteome_tax_id_file (str): path to proteome tax id file of esmecata

    Returns:
        observation_names_tax_id_names (dict): dictionary associating organism name with tax_id_name
    """
    observation_names_tax_id_names = {}

    df_proteome_tax_id = pd.read_csv(proteome_tax_id_file, sep='\t')
    for index, row in df_proteome_tax_id.iterrows():
        observation_names_tax_id_names[row['observation_name']] = row['tax_id_name']

    return observation_names_tax_id_names


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


def compute_bigecyhmm_functions_occurrence(bigecyhmm_output_file, tax_id_names_observation_names=None):
    """Read pathway_presence.tsv or function_presence.tsv created by bigecyhmm to compute the occurrence of each functions/pathways.

    Args:
        bigecyhmm_output_file (str): path to the output file of bigecyhmm (either pathway_presence.tsv or function_presence.tsv).
        tax_id_names_observation_names (dict): dictionary associating tax_id_name with organism name.

    Returns:
        function_occurrence_organisms (dict): dictionary containing functio nas key and subdict with organism as key and value of function in organism.
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

    # Compute the occurrence of functions in organism from bigecyhmm file.
    function_organisms, all_studied_organisms = compute_bigecyhmm_functions_occurrence(bigecyhmm_output_file, tax_id_names_observation_names)

    # For each sample compute the abundance of function according to the organisms.
    function_abundance_samples = {}
    function_relative_abundance_samples = {}
    function_participation_samples = {}
    for sample in sample_abundance:
        function_abundance = {}
        function_participation = {}
        for organism in sample_abundance[sample]:
            if organism not in function_participation:
                function_participation[organism] = {}
            for function_name in function_organisms:
                if organism in function_organisms[function_name]:
                    # Compute the abundance of function in the sample.
                    if function_name not in function_abundance:
                        function_abundance[function_name] = function_organisms[function_name][organism]*sample_abundance[sample][organism]
                    else:
                        function_abundance[function_name] = function_organisms[function_name][organism]*sample_abundance[sample][organism] + function_abundance[function_name]
                    # Compute the participation by each organism for the function.
                    if function_name not in function_participation[organism]:
                        function_participation[organism][function_name] = function_organisms[function_name][organism]*sample_abundance[sample][organism]
                    else:
                        function_participation[organism][function_name] = function_organisms[function_name][organism]*sample_abundance[sample][organism] + function_participation[organism][function_name]

        if sample not in function_abundance_samples:
            function_abundance_samples[sample] = {}
        if sample not in function_relative_abundance_samples:
            function_relative_abundance_samples[sample] = {}
        for function_name in function_abundance:
            if function_name not in function_abundance_samples[sample]:
                function_abundance_samples[sample][function_name] = function_abundance[function_name]
            if function_name not in function_relative_abundance_samples[sample]:
                function_relative_abundance_samples[sample][function_name] = function_abundance[function_name] / sample_tot_abundance[sample]
        function_participation_samples[sample] = function_participation

    return function_abundance_samples, function_relative_abundance_samples, function_participation_samples


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


def create_polar_plot_occurrence(function_occurrence_community_df, output_polar_plot):
    """Create polar plot from pandas dataframe with function as 'name' column' and associated ratio as a second column

    Args:
        function_occurrence_community_df (pd.DataFrame): dataframe pandas containing a column with the name of function, a second column with the relative abundance of organisms having it in the community
        output_polar_plot (path): path to the output file.
    """
    function_occurrence_community_df.reset_index(inplace=True)
    function_occurrence_community_df = function_occurrence_community_df.sort_values(['name'], ascending=True)
    function_occurrence_community_df['name'] = function_occurrence_community_df['name'].apply(lambda x: x.split(':')[1])

    fig = px.line_polar(function_occurrence_community_df, r="ratio", theta="name", line_close=True)
    fig.write_image(output_polar_plot, scale=1, width=1400, height=1200)


def create_polar_plot(df_seaborn_abundance, output_polar_plot):
    """Create polar plot from pandas dataframe with function as 'name' column' and associated ratio as a second column

    Args:
        df_seaborn_abundance (pd.DataFrame): dataframe pandas containing a column with the name of function,  a second column with the relative abundance of organisms having it in the community and a third column for the sample
        output_polar_plot (path): path to the output file.
    """

    """
    Script to make one polar plot per samples. TODO: make it more general.
    specs = [[{'type': 'polar'}]*2]*2
    fig = make_subplots(rows=2, cols=2, specs=specs)

    removed_functions = ['N-S-10:Nitric oxide dismutase', 'O-S-04:Arsenite oxidation', 'S-S-10:Polysulfide reduction']

    kept_functions = [name for name in df_seaborn_abundance['name']
                        if df_seaborn_abundance[df_seaborn_abundance['name']==name]['ratio'].max()>0]
    row = 1
    col = 1
    color = ['red', 'blue', 'green', 'purple', 'black']
    for sample in sorted(df_seaborn_abundance['sample'].unique()):
        tmp_df_seaborn_abundance = df_seaborn_abundance[df_seaborn_abundance['sample']==sample]
        tmp_df_seaborn_abundance = tmp_df_seaborn_abundance.sort_values(['name'], ascending=False)
        # Remove function
        tmp_df_seaborn_abundance = tmp_df_seaborn_abundance[~tmp_df_seaborn_abundance['name'].isin(removed_functions)]
        tmp_df_seaborn_abundance = tmp_df_seaborn_abundance[tmp_df_seaborn_abundance['name'].isin(kept_functions)]

        # Keep only name of function
        tmp_df_seaborn_abundance['name'] = tmp_df_seaborn_abundance['name'].apply(lambda x: x.split(':')[1])
        #tmp_df_seaborn_abundance = tmp_df_seaborn_abundance[tmp_df_seaborn_abundance['ratio']>0.05]

        fig.add_trace(go.Scatterpolar(
            name = sample,
            r = tmp_df_seaborn_abundance["ratio"],
            theta = tmp_df_seaborn_abundance["name"],
            ), row, col)
        if col < 2:
            col = col + 1
        else:
            col = 1
            row = row + 1

    fig.update_traces(fill='toself')
    fig.update_polars(radialaxis=dict(range=[0,1]))
    fig.write_image(output_polar_plot_1, scale=1, width=1600, height=1200)
    """
    df_seaborn_abundance = df_seaborn_abundance.sort_values(['sample', 'name'], ascending=True)
    df_seaborn_abundance['name'] = df_seaborn_abundance['name'].apply(lambda x: x.split(':')[1])

    fig = px.line_polar(df_seaborn_abundance, r="ratio", theta="name", color="sample", line_close=True)
    fig.write_image(output_polar_plot, scale=1, width=1400, height=1200)


def visualise_barplot_category(category, gene_categories, df_seaborn_abundance, output_file_name):
    """Create bar plot for functions of the associated categories from pandas dataframe with function as 'name' column' and associated ratio as a second column

    Args:
        cateogry (str): name of the function category to plot
        gene_categories (dict): adicitonary mapping function category to their respective function inferred by bigecyhmm
        df_seaborn_abundance (pd.DataFrame): dataframe pandas containing a column with the name of function,  a second column with the relative abundance of organisms having it in the community and a third column for the sample
        output_file_name (str): path to the output file.
    """
    fig, axes = plt.subplots(figsize=(40,20))
    plt.rc('font', size=30)
    kept_functions = gene_categories[category]
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
    sns.set_theme(font_scale=1)
    fig, axes = plt.subplots(figsize=(35,60))
    plt.rc('font', size=10)
    g = sns.heatmap(data=df, xticklabels=1, cmap='viridis_r', linewidths=1, linecolor='black')
    plt.tight_layout()
    plt.savefig(output_heatmap_filepath)
    plt.clf()


def create_visualisation(bigecyhmm_output, output_folder, esmecata_output_folder=None, abundance_file_path=None):
    """Create visualisation plots from esmecata, bigecyhmm output folders

    Args:
        bigecyhmm_output (str): path to bigecyhmm output folder.
        output_folder (str): path to the output folder where files will be created.
        esmecata_output_folder (str): path to esmecata output folder.
        abundance_file_path (str): path to abundance file.
    """
    start_time = time.time()

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    output_folder_occurrence = os.path.join(output_folder, 'function_occurrence')
    if not os.path.exists(output_folder_occurrence):
        os.mkdir(output_folder_occurrence)

    if abundance_file_path is not None:
        sample_abundance, sample_tot_abundance = read_abundance_file(abundance_file_path)
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
        observation_names_tax_id_names = read_esmecata_proteome_file(proteome_tax_id_file)
        tax_id_names_observation_names = {}
        for observation_name in observation_names_tax_id_names:
            tax_id_name = observation_names_tax_id_names[observation_name]
            if tax_id_name not in tax_id_names_observation_names:
                tax_id_names_observation_names[tax_id_name] = [observation_name]
            else:
                tax_id_names_observation_names[tax_id_name].append(observation_name)
    else:
        tax_id_names_observation_names = None

    logger.info("## Compute function occurrences and create visualisation.")
    logger.info("  -> Read bigecyhmm cycle output files.")
    bigecyhmm_pathway_presence_file = os.path.join(bigecyhmm_output, 'pathway_presence.tsv')
    cycle_occurrence_organisms, all_studied_organisms = compute_bigecyhmm_functions_occurrence(bigecyhmm_pathway_presence_file, tax_id_names_observation_names)
    # Compute the relative abundance of organisms by dividing for a function the number of organisms having it by the total number of organisms in the community.
    cycle_occurrence_community = []
    for index in cycle_occurrence_organisms:
        cycle_occurrence_community.append([index, len(cycle_occurrence_organisms[index])/len(all_studied_organisms)])

    cycle_occurrence_community_df = pd.DataFrame(cycle_occurrence_community, columns=['name', 'ratio'])
    cycle_occurrence_community_df.set_index('name', inplace=True)
    cycle_occurrence_community_df.to_csv(os.path.join(output_folder_occurrence, 'cycle_occurrence.tsv'), sep='\t')
    cycle_occurrences = cycle_occurrence_community_df['ratio'].to_dict()

    output_file_name = os.path.join(output_folder_occurrence, 'swarmplot_function_ratio_community.png')
    create_swarmplot_community(cycle_occurrence_community_df, output_file_name)

    logger.info("  -> Create polarplot.")
    output_polar_plot = os.path.join(output_folder_occurrence, 'polar_plot_occurrence.png')
    create_polar_plot_occurrence(cycle_occurrence_community_df, output_polar_plot)

    logger.info("  -> Create diagrams.")
    all_cycles = pd.read_csv(bigecyhmm_pathway_presence_file, sep='\t')['function'].tolist()
    diagram_data = {}
    for cycle_name in all_cycles:
        if cycle_name in cycle_occurrences:
            diagram_data[cycle_name] = (round(sum(cycle_occurrence_organisms[cycle_name].values()), 1), round(cycle_occurrences[cycle_name]*100, 1))
        else:
            diagram_data[cycle_name] = (0, 0)

    carbon_cycle_file = os.path.join(output_folder_occurrence, 'diagram_carbon_cycle.png')
    create_carbon_cycle(diagram_data, carbon_cycle_file)

    nitrogen_cycle_file = os.path.join(output_folder_occurrence, 'diagram_nitrogen_cycle.png')
    create_nitrogen_cycle(diagram_data, nitrogen_cycle_file)

    sulfur_cycle_file = os.path.join(output_folder_occurrence, 'diagram_sulfur_cycle.png')
    create_sulfur_cycle(diagram_data, sulfur_cycle_file)

    other_cycle_file = os.path.join(output_folder_occurrence, 'diagram_other_cycle.png')
    create_other_cycle(diagram_data, other_cycle_file)

    logger.info("  -> Read bigecyhmm functions output files.")
    bigecyhmm_function_presence_file = os.path.join(bigecyhmm_output, 'function_presence.tsv')
    function_occurrence_organisms, all_studied_organisms = compute_bigecyhmm_functions_occurrence(bigecyhmm_function_presence_file, tax_id_names_observation_names)
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
        sample_abundance, sample_tot_abundance = read_abundance_file(abundance_file_path)
        output_folder_abundance = os.path.join(output_folder, 'function_abundance')
        if not os.path.exists(output_folder_abundance):
            os.mkdir(output_folder_abundance)

        logger.info("  -> Read bigecyhmm cycle output files.")
        bigecyhmm_pathway_presence_file = os.path.join(bigecyhmm_output, 'pathway_presence.tsv')
        cycle_abundance_samples, cycle_relative_abundance_samples, cycle_participation_samples = compute_bigecyhmm_functions_abundance(bigecyhmm_pathway_presence_file, sample_abundance, sample_tot_abundance, tax_id_names_observation_names)

        cycle_relative_abundance_samples_df = pd.DataFrame(cycle_relative_abundance_samples)
        cycle_relative_abundance_samples_df.index.name = 'name'
        cycle_relative_abundance_samples_df.to_csv(os.path.join(output_folder_abundance, 'cycle_abundance_sample.tsv'), sep='\t')

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
        output_polar_plot = os.path.join(output_folder_abundance, 'polar_plot_abundance_samples.png')
        create_polar_plot(melted_cycle_relative_abundance_samples_df, output_polar_plot)

        logger.info("  -> Create diagrams.")
        output_folder_cycle_diagram= os.path.join(output_folder_abundance, 'cycle_diagrams_abundance')
        if not os.path.exists(output_folder_cycle_diagram):
            os.mkdir(output_folder_cycle_diagram)

        for sample in cycle_relative_abundance_samples:
            diagram_data = {}
            for cycle_name in all_cycles:
                if cycle_name in cycle_relative_abundance_samples[sample]:
                    diagram_data[cycle_name] = (round(cycle_abundance_samples[sample][cycle_name], 1), round(cycle_relative_abundance_samples[sample][cycle_name]*100, 1))
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

        logger.info("  -> Read bigecyhmm function output files.")
        bigecyhmm_function_presence_file = os.path.join(bigecyhmm_output, 'function_presence.tsv')
        function_abundance_samples, function_relative_abundance_samples, function_participation_samples = compute_bigecyhmm_functions_abundance(bigecyhmm_function_presence_file, sample_abundance, sample_tot_abundance, tax_id_names_observation_names)

        function_relative_abundance_samples_df = pd.DataFrame(function_relative_abundance_samples)
        function_relative_abundance_samples_df.index.name = 'name'
        function_relative_abundance_samples_df.to_csv(os.path.join(output_folder_abundance, 'function_abundance_sample.tsv'), sep='\t')

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
        create_heatmap_functions(function_relative_abundance_samples_df, output_heatmap_filepath)
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
    duration = time.time() - start_time
    metadata_json = {}
    metadata_json['tool_dependencies'] = {}
    metadata_json['tool_dependencies']['python_package'] = {}
    metadata_json['tool_dependencies']['python_package']['Python_version'] = sys.version
    metadata_json['tool_dependencies']['python_package']['bigecyhmm'] = bigecyhmm_version
    metadata_json['tool_dependencies']['python_package']['pandas'] = pandas_version
    metadata_json['tool_dependencies']['python_package']['plotly'] = plotly_version
    metadata_json['tool_dependencies']['python_package']['matplotlib'] = matplotlib_version
    metadata_json['tool_dependencies']['python_package']['seaborn'] = seaborn_version
    metadata_json['tool_dependencies']['python_package']['kaleido'] = kaleido_version

    metadata_json['input_parameters'] = {'esmecata_output_folder': esmecata_output_folder, 'bigecyhmm_output': bigecyhmm_output, 'output_folder': output_folder,
                                         'abundance_file_path': abundance_file_path}
    metadata_json['duration'] = duration

    metadata_file = os.path.join(output_folder, 'bigecyhmm_visualisation_metadata.json')
    with open(metadata_file, 'w') as ouput_file:
        json.dump(metadata_json, ouput_file, indent=4)


def main():
    start_time = time.time()

    parser = argparse.ArgumentParser(
        'bigecyhmm_visualisation',
        description=MESSAGE + ' For specific help on each subcommand use: esmecata {cmd} --help',
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

    parent_parser_output_folder = argparse.ArgumentParser(add_help=False)
    parent_parser_output_folder.add_argument(
        '-o',
        '--output',
        dest='output',
        required=True,
        help='Output directory path.',
        metavar='OUPUT_DIR')

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
            parent_parser_output_folder
            ],
        allow_abbrev=False)
    genomes_parser = subparsers.add_parser(
        'genomes',
        help='Creates visualisation from runs of bigecyhmm on genomes.',
        parents=[
            parent_parser_bigecyhmm, parent_parser_abundance_file,
            parent_parser_output_folder
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

    if args.abundance_file == 'false':
        abundance_file = None
    else:
        abundance_file = args.abundance_file

    if args.cmd in ['esmecata']:
        create_visualisation(args.bigecyhmm, args.output, esmecata_output_folder=args.esmecata, abundance_file_path=abundance_file)
    elif args.cmd in ['genomes']:
        create_visualisation(args.bigecyhmm, args.output, abundance_file_path=abundance_file)

    duration = time.time() - start_time
    logger.info("--- Total runtime %.2f seconds ---" % (duration))
    logger.warning(f'--- Logs written in {log_file_path} ---')
