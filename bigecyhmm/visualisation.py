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
import seaborn as sns
import matplotlib.pyplot as plt
import os
import sys

from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.express as px

import argparse
import logging
import os
import sys
import time

from bigecyhmm import __version__ as VERSION
from bigecyhmm.utils import is_valid_dir
from bigecyhmm.hmm_search import search_hmm

MESSAGE = '''
Create figures from bigecyhmm and esmecata outputs.
'''
REQUIRES = '''
Requires seaborn, pandas, plotly and kaleido.
'''

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

plt.rcParams['figure.figsize'] = [40, 20]
plt.rc('font', size=30)

def read_abundance_file(abundance_file_path):
    if abundance_file_path.endswith('.tsv'):
        input_data_df = pd.read_csv(abundance_file_path, sep='\t')
    elif abundance_file_path.endswith('.csv'):
        input_data_df = pd.read_csv(abundance_file_path)
    input_data_df.set_index('observation_name', inplace=True)
    sample_abundance = {}
    for col in input_data_df.columns:
        sample_abundance[col] = input_data_df[col].to_dict()

    return sample_abundance

def read_esmecata_proteome_file(proteome_tax_id_file):
    observation_names_tax_id_names = {}
    observation_names_tax_ids = {}
    df_proteome_tax_id = pd.read_csv(proteome_tax_id_file, sep='\t')

    for index, row in df_proteome_tax_id.iterrows():
        observation_names_tax_id_names[row['observation_name']] = row['tax_id_name']
        if row['tax_id'] not in observation_names_tax_ids:
            observation_names_tax_ids[row['tax_id']] = [row['observation_name']]
        else:
            observation_names_tax_ids[row['tax_id']].append(row['observation_name'])
    return observation_names_tax_id_names

def compute_abundance_per_tax_id(sample_abundance, observation_names_tax_id_names):
    abundance_data = {}
    for col in sample_abundance:
        for observation_name in sample_abundance[col]:
            if observation_name in observation_names_tax_id_names:
                tax_id_name = observation_names_tax_id_names[observation_name]
                if col not in abundance_data:
                    abundance_data[col] = {}
                if tax_id_name not in abundance_data[col]:
                    abundance_data[col][tax_id_name] = float(sample_abundance[col][observation_name])/100
                else:
                    abundance_data[col][tax_id_name] = float(sample_abundance[col][observation_name])/100 + float(abundance_data[col][tax_id_name])/100

    return abundance_data

def read_bigecyhmm(bigecyhmm_output, output_folder, abundance_data=None):
    data_seaborn = []
    data_seaborn_abundance = []
    data_stat = {}
    annot_folder = 'results'
    annot_table_path = os.path.join(bigecyhmm_output, 'function_presence.tsv')
    df = pd.read_csv(annot_table_path, sep='\t')
    df.set_index('function', inplace=True)
    df[annot_folder] = df.sum(axis=1)
    df = df[[annot_folder]]
    cycle_path = os.path.join(bigecyhmm_output, 'Total.R_input.txt')
    cycle_df = pd.read_csv(cycle_path, sep='\t', index_col=0, header=None)
    cycle_df.columns = ['genome', annot_folder]
    for index, row in cycle_df.iterrows():
        if index not in data_stat:
            data_stat[index] = {}
            if annot_folder not in data_stat[index]:
                data_stat[index][annot_folder] = [row[annot_folder]]
            else:
                data_stat[index][annot_folder].append(row[annot_folder])
        else:
            if annot_folder not in data_stat[index]:
                data_stat[index][annot_folder] = [row[annot_folder]]
            else:
                data_stat[index][annot_folder].append(row[annot_folder])
    cycle_df = cycle_df[annot_folder]


    function_abundance = {}
    tax_id_function = {}
    
    cycle_path = os.path.join(bigecyhmm_output, 'pathway_presence.tsv')
    pathway_presence_df = pd.read_csv(cycle_path, sep='\t', index_col=0).T

    all_tax_ids = pathway_presence_df.columns

    for index, row in pathway_presence_df.iterrows():
        for tax_id_name in pathway_presence_df.columns:
            row[tax_id_name] = int(row[tax_id_name])
            if index not in function_abundance:
                function_abundance[index] = row[tax_id_name]
            else:
                function_abundance[index] = row[tax_id_name] + function_abundance[index]
            if row[tax_id_name] > 0:
                if index not in tax_id_function:
                    tax_id_function[index] = [tax_id_name]
                else:
                    if tax_id_name not in tax_id_function[index]:
                        tax_id_function[index].append(tax_id_name)
    for index in function_abundance:
        if index in tax_id_function:
            data_seaborn.append([index, len(tax_id_function[index])/len(all_tax_ids)])
        else:
            data_seaborn.append([index, 0])
    df_seaborn_community = pd.DataFrame(data_seaborn, columns=['name', 'ratio'])
    df_seaborn_community.to_csv(os.path.join(output_folder, 'hmm_cycleboxplot_community.tsv'), sep='\t')

    if abundance_data is None:
        return df_seaborn_community, None, None
    else:
        for sample in abundance_data:
            function_abundance = {}
            all_tax_ids = []
            tax_id_function = {}
            for tax_id_name in abundance_data[sample]:
                tax_id_name_cycle_path = os.path.join(bigecyhmm_output, 'diagram_input', tax_id_name+'.R_input.txt')
                cycle_df = pd.read_csv(tax_id_name_cycle_path, sep='\t', index_col=0, header=None)
                cycle_df.columns = ['genome']
                if abundance_data[sample][tax_id_name] > 0:
                    all_tax_ids.append(tax_id_name)
                for index, row in cycle_df.iterrows():
                    if index not in function_abundance:
                        function_abundance[index] = row['genome']*abundance_data[sample][tax_id_name]
                    else:
                        function_abundance[index] = row['genome']*abundance_data[sample][tax_id_name] + function_abundance[index]
                    if row['genome']*abundance_data[sample][tax_id_name] > 0:
                        if index not in tax_id_function:
                            tax_id_function[index] = [tax_id_name]
                        else:
                            if tax_id_name not in tax_id_function[index]:
                                tax_id_function[index].append(tax_id_name)
            for index in function_abundance:
                data_seaborn_abundance.append([index, function_abundance[index], sample])
                if index in tax_id_function:
                    data_seaborn.append([index, len(tax_id_function[index])/len(all_tax_ids), sample])
                else:
                    data_seaborn.append([index, 0, sample])

        df_seaborn_sample = pd.DataFrame(data_seaborn, columns=['name', 'ratio', 'sample'])
        df_seaborn_sample_abundance = pd.DataFrame(data_seaborn_abundance, columns=['name', 'ratio',  'sample'])
        df_seaborn_sample_abundance.to_csv(os.path.join(output_folder, 'hmm_cycleboxplot_sample_abundance.tsv'), sep='\t')
        df_seaborn_sample.to_csv(os.path.join(output_folder, 'hmm_cycleboxplot_sample.tsv'), sep='\t')

        return df_seaborn_community, df_seaborn_sample, df_seaborn_sample_abundance

def create_swarmplot_community(df_seaborn, output_file_name):
    ax = sns.swarmplot(data=df_seaborn, x='name', y='ratio', s=10)
    [ax.axvline(x+.5,color='k') for x in ax.get_xticks()]
    plt.xticks(rotation=90)
    plt.savefig(output_file_name, bbox_inches="tight")
    plt.clf()

def create_swarmplot_sample(df_seaborn, output_file_name):
    ax = sns.swarmplot(data=df_seaborn, x='name', y='ratio', hue='sample', s=10)
    [ax.axvline(x+.5,color='k') for x in ax.get_xticks()]
    plt.xticks(rotation=90)
    plt.savefig(output_file_name, bbox_inches="tight")
    plt.clf()

def create_polar_plot(df_seaborn_abundance, output_polar_plot_2):
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
    fig.write_image(output_polar_plot_2, scale=1, width=1400, height=1200)


def visualisation(esmecata_output_folder, bigecyhmm_output, output_folder, abundance_file_path=None):
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    if abundance_file_path is not None: 
        sample_abundance = read_abundance_file(abundance_file_path)
    else:
        sample_abundance = None

    proteome_tax_id_file = os.path.join(esmecata_output_folder, '0_proteomes', 'proteome_tax_id.tsv')
    observation_names_tax_id_names = read_esmecata_proteome_file(proteome_tax_id_file)
    if abundance_file_path is not None: 
        abundance_data = compute_abundance_per_tax_id(sample_abundance, observation_names_tax_id_names)
    else:
        abundance_data = None

    df_seaborn_community, df_seaborn, df_seaborn_abundance = read_bigecyhmm(bigecyhmm_output, output_folder, abundance_data)

    output_file_name = os.path.join(output_folder, 'barplot_function_ratio_community.png')
    create_swarmplot_community(df_seaborn_community, output_file_name)

    if abundance_file_path is not None:
        output_file_name = os.path.join(output_folder, 'barplot_function_ratio_sample.png')
        output_file_name_abund = os.path.join(output_folder, 'barplot_function_abundance_ratio_sample.png')
        create_swarmplot_sample(df_seaborn, output_file_name)
        create_swarmplot_sample(df_seaborn_abundance, output_file_name_abund)

    if abundance_file_path is not None:
        output_polar_plot = os.path.join(output_folder, 'polar_plot_merged.png')
        create_polar_plot(df_seaborn_abundance, output_polar_plot)


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
        version='%(prog)s ' + VERSION + '\n')

    parser.add_argument(
        '--esmecata',
        dest='esmecata',
        required=True,
        help='EsMeCaTa output folder for the input file.',
        metavar='INPUT_FOLDER')

    parser.add_argument(
        '--bigecyhmm',
        dest='bigecyhmm',
        required=True,
        help='Bigecyhmm output folder for the input file.',
        metavar='INPUT_FOLDER')

    parser.add_argument(
        '--abundance-file',
        dest='abundance_file',
        required=True,
        help='Abundance file indicating the abundance for each organisms selected by EsMeCaTa.',
        metavar='INPUT_FILE')

    parser.add_argument(
        '-o',
        '--output',
        dest='output',
        required=True,
        help='Output directory path.',
        metavar='OUPUT_DIR')

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

    visualisation(args.esmecata, args.bigecyhmm, args.output, abundance_file)

    duration = time.time() - start_time
    logger.info("--- Total runtime %.2f seconds ---" % (duration))
    logger.warning(f'--- Logs written in {log_file_path} ---')
