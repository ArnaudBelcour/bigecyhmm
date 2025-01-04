import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import sys

from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.express as px

# Requires seaborn, pandas, plotly and kaleido

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

def read_bigecyhmm(bigecyhmm_output, abundance_data):
    tmp_df = []
    cycle_tmp_df = []
    data_seaborn = []
    data_seaborn_abundance = []
    data_stat = {}
    annot_folder = 'results'
    annot_table_path = os.path.join(bigecyhmm_output, 'function_presence.tsv')
    df = pd.read_csv(annot_table_path, sep='\t')
    df.set_index('function', inplace=True)
    df[annot_folder] = df.sum(axis=1)
    df = df[[annot_folder]]
    tmp_df.append(df)
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
    cycle_tmp_df.append(cycle_df)

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

    df_seaborn = pd.DataFrame(data_seaborn, columns=['name', 'ratio', 'sample'])
    df_seaborn_abundance = pd.DataFrame(data_seaborn_abundance, columns=['name', 'ratio',  'sample'])
    df_seaborn_abundance.to_csv('hmm_cycleboxplot_abundance.tsv', sep='\t')
    df_seaborn.to_csv('hmm_cycleboxplot.tsv', sep='\t')

    return df_seaborn, df_seaborn_abundance

def create_swarmplot(df_seaborn, df_seaborn_abundance, output_file_name, output_file_name_abund):
    ax = sns.swarmplot(data=df_seaborn, x='name', y='ratio', hue='sample', s=10)
    [ax.axvline(x+.5,color='k') for x in ax.get_xticks()]
    plt.xticks(rotation=90)
    plt.savefig(output_file_name, bbox_inches="tight")
    plt.clf()

    ax = sns.swarmplot(data=df_seaborn_abundance, x='name', y='ratio', hue='sample', s=10)
    [ax.axvline(x+.5,color='k') for x in ax.get_xticks()]
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xticks(rotation=90)
    plt.savefig(output_file_name_abund, bbox_inches="tight")
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


def visualisation(abundance_file_path, esmecata_output_folder, bigecyhmm_output, output_folder):
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    sample_abundance = read_abundance_file(abundance_file_path)
    proteome_tax_id_file = os.path.join(esmecata_output_folder, '0_proteomes', 'proteome_tax_id.tsv')
    observation_names_tax_id_names = read_esmecata_proteome_file(proteome_tax_id_file)
    abundance_data = compute_abundance_per_tax_id(sample_abundance, observation_names_tax_id_names)
    df_seaborn, df_seaborn_abundance = read_bigecyhmm(bigecyhmm_output, abundance_data)

    output_file_name = os.path.join(output_folder, 'barplot_function_ratio.png')
    output_file_name_abund = os.path.join(output_folder, 'barplot_function_abundance_ratio.png')

    create_swarmplot(df_seaborn, df_seaborn_abundance, output_file_name, output_file_name_abund)

    output_polar_plot = os.path.join(output_folder, 'polar_plot_merged.png')
    create_polar_plot(df_seaborn_abundance, output_polar_plot)

def main():
    abundance_file_path = sys.argv[1]
    esmecata_output_folder = sys.argv[2]
    bigecyhmm_output = sys.argv[3]
    output_folder = sys.argv[4]
    visualisation(abundance_file_path, esmecata_output_folder, bigecyhmm_output, output_folder)