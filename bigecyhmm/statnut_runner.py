import os
import data_stats
import donut_plot
import table_plot
import combine_plots
import tempfile

#Default file paths and groups. 
DEFAULT_INPUT_TSV = os.path.join('function_abundance', 'cycle_abundance_sample.tsv')
DEFAULT_STATS_TSV = os.path.join('plots', 'group_stats.tsv')
DEFAULT_CLEANED_TSV = os.path.join('plots', 'cleaned_data.tsv')
DEFAULT_DONUT_PNG = os.path.join('plots', 'group_medians_donut.png')
DEFAULT_TABLE_PNG = os.path.join('plots', 'group_stats_table.png')
DEFAULT_BACKGROUND_PATH = os.path.join(os.path.dirname(__file__), 'network_background_v4.png') 
#Default sample groups: expect sample_groups.tsv to be located in the pipeline input folder (one level up from visualisation outputs)
DEFAULT_SAMPLE_GROUPS = os.path.join('..', 'sample_groups.tsv')
DEFAULT_GROUPS = {'All Samples': ['']}


def statNut_run(input_tsv: str = DEFAULT_INPUT_TSV,
         sample_groups_tsv: str = DEFAULT_SAMPLE_GROUPS,
         groups: dict | None = None,
         stats_csv: str = DEFAULT_STATS_TSV,
         cleaned_csv: str = DEFAULT_CLEANED_TSV,
         donut_png: str = DEFAULT_DONUT_PNG,
         table_png: str = DEFAULT_TABLE_PNG,
         background_path: str = DEFAULT_BACKGROUND_PATH,
         background_offset: tuple = (0.016, -0.006),      #adjust background image position (right, up)
         background_scale: float = 0.60,                 #adjust background image scale relative to donut
         ):  

    # If the provided sample_groups_tsv does not exist, create a temporary one which assigns all samples to a single group called "all samples". 
    temp_groups_path = None
    created_temp_groups = False
    if not sample_groups_tsv or not os.path.exists(sample_groups_tsv):
        # create a temporary TSV that instructs the pipeline to select all samples
        tf = tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.tsv')
        tf.write('sample\tgroup\n')
        tf.write('\tAll Samples\n')
        tf.close()
        temp_groups_path = tf.name
        sample_groups_tsv = temp_groups_path
        created_temp_groups = True

    try:
        #load input data
        df = data_stats.load_data(input_tsv)

        #if no groups provided by caller, default groups are used. 
        if groups is None:
            groups = DEFAULT_GROUPS

        #compute statistics and save cleaned + rounded CSVs for later plotting
        cleaned_df, df_res = data_stats.compute_and_save(input_tsv=input_tsv,
                                                          groups_tsv=sample_groups_tsv,
                                                          output_csv=stats_csv,
                                                          cleaned_output_csv=cleaned_csv)

        # obtain group resolution and formatted display dataframe for plotting
        group_names, group_col_names, groups_dict = data_stats.get_group_col_names(df, sample_groups_tsv)
        results = data_stats.run_statistics(df, group_col_names)
        _, display_df = data_stats.results_to_dataframe(results, groups_dict)

        # use index labels for plotting
        metabolic_labels = df.index.tolist()

        # avoid changing original display_df
        display_df = display_df.copy()

        os.makedirs(os.path.dirname(donut_png), exist_ok=True)
        donut_plot.plot_donut(
            cleaned_df,
            groups_dict,
            metabolic_labels=metabolic_labels,
            group_col_names=group_col_names,
            output_path=donut_png,
            background_path=background_path,
            background_scale=background_scale,
            background_offset=background_offset,
        )

        os.makedirs(os.path.dirname(table_png), exist_ok=True)
        table_plot.plot_table(display_df, output_path=table_png)

        combined_png = os.path.join(os.path.dirname(donut_png), 'combined_donut_table.png')
        combine_plots.combine_images_side_by_side(donut_png, table_png, combined_png, padding=24, match_height='left')

    finally:
        #remove temporary groups file if it was created
        if created_temp_groups and temp_groups_path and os.path.exists(temp_groups_path):
            os.remove(temp_groups_path)
            

if __name__ == '__main__':
    statNut_run()
