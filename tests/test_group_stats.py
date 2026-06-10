import os
import pandas as pd
import shutil


def test_generate_hmm_template():
    
    from bigecyhmm.custom_tsv_parser import generate_hmm_template

    input_file = os.path.join('input_data', 'group_stats', 'custom_hydrogen_central_cycles.tsv')
    output_folder = 'output_folder' 
    os.makedirs(output_folder, exist_ok=True)

    #set up hmm_custom_db_df and metabolic_functions, which are the inputs of generate_hmm_template. 
    #Unfortunately need to do it here to keep functions independent.
    custom_db_df = pd.read_csv(input_file, sep='\t')
    function_custom_db_df = custom_db_df[custom_db_df['Type']=='FUNCTION'].copy()
    function_custom_db_df['Graph_No'] = [graph_nb.split('.')[0] for graph_nb in function_custom_db_df['Graph_No']]
    metabolic_functions = function_custom_db_df.set_index('Graph_No')['ID'].to_dict()
    hmm_custom_db_df = custom_db_df[custom_db_df['Type']=='HMM'].copy()

    #run the function to test
    hmm_template_df = generate_hmm_template(hmm_custom_db_df, metabolic_functions)

    #convert df to and from csv to have same dtype and format as saved expected_hmm_template_df, to be able to compare them with assert_frame_equal
    hmm_template_df.to_csv(os.path.join(output_folder, 'hmm_template_df.tsv'), sep='\t')
    hmm_template_df = pd.read_csv(os.path.join(output_folder, 'hmm_template_df.tsv'), sep='\t')

    expected_hmm_template_df = pd.read_csv(os.path.join('input_data', 'group_stats', 'expected_hmm_template_df.tsv'), sep='\t').reset_index(drop=True)

    pd.testing.assert_frame_equal(hmm_template_df, expected_hmm_template_df)  

    shutil.rmtree(output_folder)



def test_translate_expression():
    
    from bigecyhmm.custom_tsv_parser import translate_expression

    expression = '((A AND B) OR (C AND D)) AND NOT E'
    func_major = 'Arbitrary metabolic function'
    id_to_hmms = {'Arbitrary metabolic function': {'A': 'a.hmm', 'B': 'b.hmm', 'C': 'c.hmm', 'D': 'd.hmm', 'E': 'e.hmm'}}

    expected_expression = '((a.hmm and b.hmm) or (c.hmm and d.hmm)) and not e.hmm'

    translated_expression = translate_expression(expression, func_major, id_to_hmms)

    assert translated_expression == expected_expression


def test_generate_pathway_template():

    from bigecyhmm.custom_tsv_parser import generate_pathway_template   

    input_file = os.path.join('input_data', 'group_stats', 'custom_hydrogen_central_cycles.tsv')
    output_folder = 'output_folder'
    os.makedirs(output_folder, exist_ok=True)

    #Needed to create variables: function_custom_db_df and metabolic_pathways_to_hmms
    custom_db_df = pd.read_csv(input_file, sep='\t')
    function_custom_db_df = custom_db_df[custom_db_df['Type']=='FUNCTION'].copy()
    function_custom_db_df['Graph_No'] = [graph_nb.split('.')[0] for graph_nb in function_custom_db_df['Graph_No']]

    metabolic_functions = function_custom_db_df.set_index('Graph_No')['ID'].to_dict()
    hmm_custom_db_df = custom_db_df[custom_db_df['Type']=='HMM'].copy()
    hmm_custom_db_df['Function_Nb'] = [graph_nb.split('.')[0] for graph_nb in hmm_custom_db_df['Graph_No']]

    metabolic_pathway_to_hmms = {}
    for metabolic_nb in metabolic_functions:
        metabolic_pathway = metabolic_functions[metabolic_nb]
        tmp_hmm_custom_db_df = hmm_custom_db_df[hmm_custom_db_df['Function_Nb']==metabolic_nb]
        metabolic_pathway_to_hmms[metabolic_pathway] = tmp_hmm_custom_db_df.set_index('ID')['Implementation'].to_dict()
        
    #run the function to test
    pathway_template_df = generate_pathway_template(function_custom_db_df, metabolic_pathway_to_hmms)

    #convert df to and from csv to have same dtype and format as saved expected_pathway_template_df, to be able to compare them with assert_frame_equal
    pathway_template_df.to_csv(os.path.join(output_folder, 'pathway_template_df.tsv'), sep="\t")
    pathway_template_df = pd.read_csv(os.path.join(output_folder, 'pathway_template_df.tsv'), sep="\t")

    expected_pathway_template_df = pd.read_csv(os.path.join('input_data', 'group_stats', 'expected_pathway_template_df.tsv'), sep='\t')

    pd.testing.assert_frame_equal(pathway_template_df, expected_pathway_template_df)

    shutil.rmtree(output_folder)


def test_build_bipartite_graph():

    from bigecyhmm.custom_tsv_parser import build_bipartite_graph

    input_file = os.path.join('input_data', 'group_stats', 'custom_hydrogen_central_cycles.tsv')

    custom_db_df = pd.read_csv(input_file, sep='\t')
    function_custom_db_df = custom_db_df[custom_db_df['Type']=='FUNCTION']

    bipartite_graph = build_bipartite_graph(function_custom_db_df)

    assert bipartite_graph.number_of_nodes() == 100
    assert bipartite_graph.number_of_edges() == 106


def test_generate_custom_db_from_tsv_one_file():

    from bigecyhmm.custom_tsv_parser import generate_custom_db_from_tsv_one_file

    input_file = os.path.join('input_data', 'group_stats', 'custom_hydrogen_central_cycles.tsv')
    output_folder = 'output_folder'    

    generate_custom_db_from_tsv_one_file(input_file, output_folder)

    assert os.path.exists(os.path.join(output_folder, 'database', 'hmm_template_file.tsv'))
    assert os.path.exists(os.path.join(output_folder, 'database', 'pathway_template_file.tsv'))
    assert os.path.exists(os.path.join(output_folder, 'database', 'input_graph.graphml'))

    shutil.rmtree(output_folder)


def test_get_group_col_names():

    from bigecyhmm.group_analysis import get_group_col_names

    input_folder = os.path.join('input_data', 'group_stats')

    expected_group_names = ['first_group', 'second_group', 'third_group']
    expected_resolved_groups = [['sample_1', 'sample_2', 'sample_3'], ['sample_4', 'sample_5', 'sample_6'], ['sample_7', 'sample_8', 'sample_9']]
    expected_groups_dict = {
        'first_group': ['sample_1', 'sample_2', 'sample_3'],   
        'second_group': ['sample_4', 'sample_5', 'sample_6'],   
        'third_group': ['sample_7', 'sample_8', 'sample_9'],    
    }

    #load cycle_abundance_sample.tsv
    df = pd.read_csv(os.path.join(input_folder, 'cycle_abundance_sample.tsv'), sep='\t', index_col=0)
    mapping_df = pd.read_csv(os.path.join(input_folder, 'assigned_groups.tsv'), sep='\t', dtype=str)

    #test get_group_col_names
    group_names, resolved_groups, groups_dict = get_group_col_names(df, mapping_df)
   
    assert group_names == expected_group_names
    assert resolved_groups == expected_resolved_groups
    assert groups_dict == expected_groups_dict

    

def test_run_statistics():

    #here, both functions are tested together, because run_statistics returns a dataclass of results which is easier to assess by transforming it into a df. 
    from bigecyhmm.group_analysis import run_statistics
    from bigecyhmm.group_analysis import results_to_dataframe
    from pandas.testing import assert_frame_equal
    

    input_folder = os.path.join('input_data', 'group_stats')
    output_folder = 'output_folder'
    os.makedirs(output_folder, exist_ok=True)
    
    df = pd.read_csv(os.path.join(input_folder, 'cycle_abundance_sample.tsv'), sep='\t', index_col=0)
    resolved_groups = [['sample_1', 'sample_2', 'sample_3'], ['sample_4', 'sample_5', 'sample_6'], ['sample_7', 'sample_8', 'sample_9']]

    stat_results = run_statistics(df, resolved_groups)
    groups_dict = {
        'first_group': ['sample_1', 'sample_2', 'sample_3'],   
        'second_group': ['sample_4', 'sample_5', 'sample_6'],   
        'third_group': ['sample_7', 'sample_8', 'sample_9'],    
    }


    display_df, stat_results_df = results_to_dataframe(stat_results, groups_dict)

    #need to convert to csv, and read csv again to have same format and dtype as expected_stat_results_df (saved as csv), to be able to compare them with assert_frame_equal
    stat_results_df.to_csv(os.path.join(output_folder, 'stat_results_df.tsv'), sep="\t")
    stat_results_df = pd.read_csv(os.path.join(output_folder, 'stat_results_df.tsv'), sep="\t")

    expected_stat_results_df = pd.read_csv(os.path.join(input_folder, 'expected_stat_results_df.tsv'), sep="\t")

    assert_frame_equal(stat_results_df, expected_stat_results_df)

    shutil.rmtree(output_folder)
    

def test_donut_plot():

    from bigecyhmm.group_plot import plot_donut

    input_folder = os.path.join('input_data', 'group_stats')
    output_folder = 'output_folder'

    groups_dict = {
        'first_group': ['sample_1', 'sample_2', 'sample_3'],
        'second_group': ['sample_4', 'sample_5', 'sample_6'],
        'third_group': ['sample_7', 'sample_8', 'sample_9'],
    }

    df = pd.read_csv(os.path.join(input_folder, 'cycle_abundance_sample.tsv'), sep='\t', index_col=0)

    metabolic_labels = df.index.astype(str).tolist()

    plot_donut(
        df,
        groups_dict,
        metabolic_labels,
        output_path=os.path.join(output_folder, 'donut_plot.png'),
        background_path=None,
        background_scale=1.0,
        background_offset=(0.0, 0.0),
        group_col_names=None
    )

    assert os.path.exists(os.path.join(output_folder, 'donut_plot.png'))

    shutil.rmtree(output_folder)


def test_donut_plot_background_image():

    from bigecyhmm.group_plot import plot_donut
    from bigecyhmm import TEMPLATE_BACKGROUND_BIGECYHMM

    input_folder = os.path.join('input_data', 'group_stats')
    output_folder = 'output_folder'

    groups_dict = {
        'first_group': ['sample_1', 'sample_2', 'sample_3'],
        'second_group': ['sample_4', 'sample_5', 'sample_6'],
        'third_group': ['sample_7', 'sample_8', 'sample_9'],
    }

    df = pd.read_csv(os.path.join(input_folder, 'cycle_abundance_sample.tsv'), sep='\t', index_col=0)

    metabolic_labels = df.index.astype(str).tolist()

    plot_donut(
        df,
        groups_dict,
        metabolic_labels,
        output_path=os.path.join(output_folder, 'donut_plot.png'),
        background_path=TEMPLATE_BACKGROUND_BIGECYHMM,
        background_scale=1.0,
        background_offset=(0.0, 0.0),
        group_col_names=None
    )

    assert os.path.exists(os.path.join(output_folder, 'donut_plot.png'))

    shutil.rmtree(output_folder)



def test_plot_table():

    from bigecyhmm.group_plot import plot_table

    input_folder = os.path.join('input_data', 'group_stats')
    output_folder = 'output_folder'

    df = pd.read_csv(os.path.join(input_folder, 'expected_stat_results_df.tsv'), sep='\t', index_col=0)

    plot_table(
        df,
        output_path=os.path.join(output_folder, 'table_plot.png')
    )

    assert os.path.exists(os.path.join(output_folder, 'table_plot.png'))

    shutil.rmtree(output_folder)


def test_combine_images_side_by_side():
    
    from bigecyhmm.group_plot import combine_images_side_by_side

    input_folder = os.path.join('input_data', 'group_stats')
    output_folder = 'output_folder'
    os.makedirs(output_folder, exist_ok=True)

    image1_path = os.path.join(input_folder, 'donut_plot.png')
    image2_path = os.path.join(input_folder, 'table_plot.png')

    combine_images_side_by_side(
        image1_path,
        image2_path,
        output_path=os.path.join(output_folder, 'combined_plot.png')
    )

    assert os.path.exists(os.path.join(output_folder, 'combined_plot.png'))

    shutil.rmtree(output_folder)

