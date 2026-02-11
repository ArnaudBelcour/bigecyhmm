# Copyright (C) 2026 Michael Baumgartner, Arnaud Belcour - Univ. Grenoble Alpes, Inria, Grenoble, France Microcosme
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

import os
import pandas as pd
import re
import networkx as nx

from bigecyhmm.utils import is_valid_dir


def generate_hmm_template(hmm_custom_db_df, metabolic_functions):
    """ From DataFrame containing HMMs information reconstruct format of dataframe for hmm_template.

    Args:
        hmm_custom_db_df (pandas DataFrame): dataframe containing description of metabolic pathways
        metabolic_functions (dict): dictionary associating metabolic pathways number to their names

    Returns:
        hmm_template_df (pandas DataFrame): dataframe at the format of hmm_template
    """
    hmm_custom_db_df['#Entry'] = [graph_nb.split('.')[0] for graph_nb in hmm_custom_db_df['Graph_No']]
    hmm_custom_db_df['Category'] = [metabolic_functions[graph_nb] for graph_nb in hmm_custom_db_df['#Entry']]
    hmm_custom_db_df['Function'] = ['' for graph_nb in hmm_custom_db_df['#Entry']]

    hmm_custom_db_df['Gene abbreviation'] = hmm_custom_db_df['ID']
    hmm_custom_db_df['Gene name'] = hmm_custom_db_df['enzyme_long']
    hmm_custom_db_df['Hmm file'] = hmm_custom_db_df['Implementation']
    hmm_custom_db_df['Corresponding KO'] = hmm_custom_db_df['kegg_ortholog']

    hmm_custom_db_df['Reaction'] = ['' for graph_nb in hmm_custom_db_df['#Entry']]
    hmm_custom_db_df['Substrate'] = ['' for graph_nb in hmm_custom_db_df['#Entry']]
    hmm_custom_db_df['Product'] = ['' for graph_nb in hmm_custom_db_df['#Entry']]
    hmm_custom_db_df['Hmm detecting threshold'] = hmm_custom_db_df['HMM_threshold']
    hmm_custom_db_df['source'] = hmm_custom_db_df['HMM_source']

    hmm_template_df = hmm_custom_db_df[['#Entry', 'Category','Function' , 'Gene abbreviation', 'Gene name',
                                        'Hmm file', 'Corresponding KO', 'Reaction', 'Substrate', 'Product',
                                        'Hmm detecting threshold', 'source']]

    return hmm_template_df


def translate_expression(expr: str, func_major: str, id_to_hmms: dict) -> str:
    """ Convert string of boolean expression of genes into boolean expression of HMMs as required by bigecyhmm

    Args:
        expr (str): boolean expression with AND, OR, NOT showing required genes to infer metabolic pathway
        func_major (str): name of the associated metabolic pathway
        id_to_hmms (dict): dictionary associating metabolic pathways with subdictionary linking genes to HMMs

    Returns:
        converted_boolean_expression (str): converted boolean expression into a readable one by bigecyhmm
    """
    #translate FUNCTION expressions: normalize textual operators to lowercase
    #and replace ID tokens with HMM filenames within the same func_major.
    if not expr:
        return expr

    #transform text operators to lowercase
    expr = re.sub(r"\bAND\b", "and", expr, flags=re.IGNORECASE)
    expr = re.sub(r"\bOR\b", "or", expr, flags=re.IGNORECASE)
    expr = re.sub(r"\bNOT\b", "not", expr, flags=re.IGNORECASE)

    #preserve parentheses and whitespace; process other tokens replacing IDs
    parts = re.split(r'(\(|\))', expr)
    out_parts = []
    hmms_for_func = id_to_hmms.get(func_major, {})

    for p in parts:
        if p in {'(', ')'}:
            out_parts.append(p)
            continue

        #split by whitespace but keep separators
        subparts = re.split(r'(\s+)', p)
        newsub = []
        for s in subparts:
            if s.isspace() or s == '':
                newsub.append(s)
                continue
            key = s.strip()
            #leave operator words, make sure they are lowercase (should be the case through regex earlier), and replace ID by HMM if available. 
            if key.lower() in {'and', 'or', 'not'}:
                newsub.append(key.lower())
            else:
                newsub.append(hmms_for_func.get(key, s))
        out_parts.append(''.join(newsub))

    converted_boolean_expression = ''.join(out_parts).strip()
    return converted_boolean_expression


def generate_pathway_template(function_custom_db_df, gene_name_to_hmms):
    """ From DataFrame containing metabolic pathways reconstruct format of dataframe for pathway_template.

    Args:
        function_custom_db_df (pandas DataFrame): dataframe containing description of metabolic pathways
        gene_name_to_hmms (dict): dictionary associating metabolic pathways with subdictionary linking genes to HMMs

    Returns:
        pathway_template_df (pandas DataFrame): dataframe at the format of pathway template
    """
    function_custom_db_df['Pathways'] = function_custom_db_df['ID']
    function_custom_db_df['HMMs'] = [translate_expression(row['Implementation'], row['ID'], gene_name_to_hmms) for index, row in function_custom_db_df.iterrows()]
    function_custom_db_df['reactants'] = function_custom_db_df['function_input']
    function_custom_db_df['products'] = function_custom_db_df['function_output']
    function_custom_db_df['formula'] = function_custom_db_df['function_input'] + '->' + function_custom_db_df['function_output']

    pathway_template_df = function_custom_db_df[['Pathways', 'HMMs', 'reactants', 'products', 'formula']]

    return pathway_template_df


def build_bipartite_graph(function_custom_db_df):
    """ From DataFrame containing metabolic pathways reconstruct bipartite network graph

    Args:
        function_custom_db_df (pandas DataFrame): dataframe containing description of metabolic pathways

    Returns:
        bipartie_graph (networkx DirectedGraph): biparite directed graph of metabolic pathways
    """
    bipartie_graph = nx.DiGraph()

    for index, row in function_custom_db_df.iterrows():
        function_id = row['ID']
        bipartie_graph.add_node(function_id, type='Function')
        bipartie_graph.nodes[function_id]['label'] = function_id

        for substrate in row['function_input'].split(', '):
            if not bipartie_graph.has_node(substrate):
                bipartie_graph.add_node(substrate, type='Metabolite')
                bipartie_graph.nodes[substrate]['label'] = substrate
            bipartie_graph.add_edge(substrate, function_id)

        for product in row['function_output'].split(', '):
            if not bipartie_graph.has_node(product):
                bipartie_graph.add_node(product, type='Metabolite')
                bipartie_graph.nodes[product]['label'] = product
            bipartie_graph.add_edge(function_id, product)

    return bipartie_graph


def generate_custom_db_from_tsv_one_file(custom_database_input, output_folder):
    """ From a tsv file containing all information, reconstruct the different files required by bigecyhmm

    Args:
        custom_database_input (str): filepath to tsv file containing all information
        output_folder (str): path to output folder

    Returns:
        hmm_template_path (str): path to custom hmm_template file
        pathway_template_path (str): path to custom pathway_template file
        bipartie_graph (str): path to custom bipartite graph file
    """
    database_folder = os.path.join(output_folder, 'database')
    is_valid_dir(database_folder)

    custom_db_df = pd.read_csv(custom_database_input, sep='\t')
    custom_db_df['Graph_No'] = custom_db_df['Graph_No'].astype(str)

    function_custom_db_df = custom_db_df[custom_db_df['Type']=='FUNCTION']

    function_custom_db_df['Graph_No'] = [graph_nb.split('.')[0] for graph_nb in function_custom_db_df['Graph_No']]
    metabolic_functions = function_custom_db_df.set_index('Graph_No')['ID'].to_dict()

    hmm_custom_db_df = custom_db_df[custom_db_df['Type']=='HMM']

    hmm_custom_db_df['Function_Nb'] = [graph_nb.split('.')[0] for graph_nb in hmm_custom_db_df['Graph_No']]
    metabolic_pathway_to_hmms = {}
    for metabolic_nb in metabolic_functions:
        metabolic_pathway = metabolic_functions[metabolic_nb]
        tmp_hmm_custom_db_df = hmm_custom_db_df[hmm_custom_db_df['Function_Nb']==metabolic_nb]
        metabolic_pathway_to_hmms[metabolic_pathway] = tmp_hmm_custom_db_df.set_index('ID')['Implementation'].to_dict()

    hmm_template_df = generate_hmm_template(hmm_custom_db_df, metabolic_functions)
    hmm_template_path = os.path.join(database_folder, 'hmm_template_file.tsv')
    hmm_template_df.to_csv(hmm_template_path, sep='\t', index=False)

    pathway_template_df = generate_pathway_template(function_custom_db_df, metabolic_pathway_to_hmms)
    pathway_template_path = os.path.join(database_folder, 'pathway_template_file.tsv')
    pathway_template_df.to_csv(pathway_template_path, sep='\t', index=False)

    #build bipartite graph and write as GraphML into database folder
    bipartie_graph = build_bipartite_graph(function_custom_db_df)

    graphml_path = os.path.join(database_folder, 'input_graph.graphml')
    nx.write_graphml(bipartie_graph, str(graphml_path), encoding='utf-8', named_key_ids=True)

    return hmm_template_path, pathway_template_path, bipartie_graph
