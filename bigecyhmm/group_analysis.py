import os
import tempfile
import pandas as pd
import numpy as np

from dataclasses import dataclass
from typing import List, Tuple, Optional

from bigecyhmm.group_plot import plot_donut, plot_table, combine_images_side_by_side, combine_images_side_by_side, get_group_col_names

from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests

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


#DEFAULTS
ROUND_TOLERANCE = 1e-12  #small tolerance for near-zero checks 
ALPHA = 0.05             #significance level


@dataclass
class FunctionStatResults:
    name: str
    medians: List[Optional[float]]
    counts: List[int]
    iq_range: List[Optional[float]]
    percent_zero: List[Optional[float]]
    H: Optional[float]
    eps2: Optional[float]
    p: Optional[float]
    p_bh: Optional[float]
    significant: bool


#load input tsv, first column (metabo. funct. name) is index. 
def load_data(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep='\t', index_col=0)
    return df


#organize the stats pipeline: load, calculate, save, etc
def compute_and_save(input_tsv: str, groups_tsv: str, output_csv: str, cleaned_output_csv: str) -> Tuple[pd.DataFrame, pd.DataFrame]:

    df = load_data(input_tsv)
    group_names, group_col_names, groups_dict = get_group_col_names(df, groups_tsv)

    cleaned_df = df.copy()
    os.makedirs(os.path.dirname(cleaned_output_csv), exist_ok=True)
    cleaned_df.to_csv(cleaned_output_csv)
    print(f"Saved cleaned data to {cleaned_output_csv}")

    results = run_statistics(df, group_col_names)
    df_res, _ = results_to_dataframe(results, groups_dict)

    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    df_res.to_csv(output_csv, index=False)
    print(f"Saved numeric stats to {output_csv}")

    return cleaned_df, df_res


#Main part for stat calculation. 
def run_statistics(df: pd.DataFrame, group_col_names: List[List[str]]) -> List[FunctionStatResults]:

    index_names = df.index.tolist()
    results: List[FunctionStatResults] = []

    #precompute per-group Series (one Series per group) using pandas aggregations
    med_series_list = [df[cols].median(axis=1).astype(float) if cols else pd.Series(np.nan, index=df.index) for cols in group_col_names]
    count_series_list = [df[cols].count(axis=1).astype(int) if cols else pd.Series(0, index=df.index) for cols in group_col_names]
    iqr_series_list = [ (df[cols].quantile(0.75, axis=1) - df[cols].quantile(0.25, axis=1)).astype(float) if cols else pd.Series(np.nan, index=df.index) for cols in group_col_names]

    percent_zero_series_list: List[pd.Series] = []
    for cols in group_col_names:
        if cols:
            is_zero = df[cols].map(lambda v: np.isclose(v, 0.0, atol=ROUND_TOLERANCE)).sum(axis=1)
            percent_zero_series_list.append((is_zero / len(cols) * 100).astype(float))
        else:
            percent_zero_series_list.append(pd.Series(np.nan, index=df.index))

    for name in index_names:
        med_row = [float(s.loc[name]) if pd.notna(s.loc[name]) else np.nan for s in med_series_list]
        counts = [int(s.loc[name]) if pd.notna(s.loc[name]) else 0 for s in count_series_list]
        iq_range = [float(s.loc[name]) if pd.notna(s.loc[name]) else np.nan for s in iqr_series_list]
        percent_zero = [float(s.loc[name]) if pd.notna(s.loc[name]) else np.nan for s in percent_zero_series_list]

        #numeric arrays per group for Kruskal
        data_per_group = [df.loc[name, cols].dropna().values.astype(float) if len(cols) > 0 else np.array([], dtype=float) for cols in group_col_names]

        #decide whether Kruskal–Wallis can be computed
        valid_groups = [arr for arr in data_per_group if len(arr) > 0]
        do_kw = len(valid_groups) >= 2
        if do_kw:
            all_values = np.concatenate(valid_groups)
            if all_values.size == 0 or np.all(np.isclose(all_values, 0.0, atol=ROUND_TOLERANCE)) or len(np.unique(all_values)) <= 1:
                do_kw = False

        if do_kw:
            H_stat, p_value = kruskal(*valid_groups)
            n_valid = len(valid_groups)
            tot_n = sum(len(g) for g in valid_groups)
            if tot_n - n_valid <= 0:
                effect_size_eps2 = np.nan
            else:
                effect_size_eps2 = (H_stat - n_valid + 1) / (tot_n - n_valid)
                effect_size_eps2 = np.maximum(effect_size_eps2, 0)
        else:
            H_stat, p_value, effect_size_eps2 = np.nan, np.nan, np.nan

        fr = FunctionStatResults(
            name=name,
            medians=med_row,
            counts=counts,
            iq_range=iq_range,
            percent_zero=percent_zero,
            H=H_stat,
            eps2=effect_size_eps2,
            p=p_value,
            p_bh=np.nan,
            significant=False,
        )
        results.append(fr)

    #Benjamini–Hochberg correction
    p_series = pd.Series([r.p for r in results], index=[r.name for r in results], dtype='float64')
    valid_mask = p_series.notna()
    if valid_mask.any():
        valid_index = p_series[valid_mask].index.tolist()
        valid_p = p_series[valid_mask].values
        reject, pvals_corr, _, _ = multipletests(valid_p, alpha=ALPHA, method='fdr_bh')
        pvals_dict = dict(zip(valid_index, pvals_corr))
        reject_dict = dict(zip(valid_index, reject))
        for r in results:
            if r.name in pvals_dict:
                r.p_bh = float(pvals_dict[r.name]) if not np.isnan(pvals_dict[r.name]) else np.nan
                r.significant = bool(reject_dict.get(r.name, False))
            else:
                r.p_bh = np.nan
                r.significant = False
    else:
        for r in results:
            r.p_bh = np.nan
            r.significant = False

    return results


def results_to_dataframe(results: List[FunctionStatResults], groups: dict) -> Tuple[pd.DataFrame, pd.DataFrame]:

    group_names = list(groups.keys())
    rows = []
    for r in results:
        row = {
            'Metabolic Function': r.name,
            'H': r.H,
            'eps2': r.eps2,
            'p': r.p,
            'p_bh': r.p_bh,
            'significant': r.significant
        }
        for i, g in enumerate(group_names):
            row[f'Median {g}'] = r.medians[i] if i < len(r.medians) else np.nan
            row[f'n {g}'] = r.counts[i] if i < len(r.counts) else 0
            row[f'iq_range {g}'] = r.iq_range[i] if i < len(r.iq_range) else np.nan
            row[f'%0 {g}'] = r.percent_zero[i] if i < len(r.percent_zero) else np.nan
        rows.append(row)
    df_res = pd.DataFrame(rows)

    #lightweight formatting for display DataFrame
    display_df = df_res.copy()
    med_cols = [c for c in display_df.columns if c.startswith('Median ')]
    for c in med_cols:
        display_df[c] = display_df[c].apply(lambda x: f"{x:.4f}" if pd.notna(x) else "")
    display_df['H'] = display_df['H'].apply(lambda x: f"{x:.3f}" if pd.notna(x) else "")
    display_df['eps2'] = display_df['eps2'].apply(lambda x: f"{x:.4f}" if pd.notna(x) else "")
    display_df['p'] = display_df['p'].apply(lambda x: f"{x:.4g}" if pd.notna(x) else "")
    display_df['p_bh'] = display_df['p_bh'].apply(lambda x: f"{x:.4g}" if pd.notna(x) else "n/a")
    display_df['significant'] = display_df['significant'].apply(lambda x: 'Yes' if x else 'No')

    return df_res, display_df


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
        df = load_data(input_tsv)

        #if no groups provided by caller, default groups are used. 
        if groups is None:
            groups = DEFAULT_GROUPS

        #compute statistics and save cleaned + rounded CSVs for later plotting
        cleaned_df, df_res = compute_and_save(input_tsv=input_tsv,
                                                          groups_tsv=sample_groups_tsv,
                                                          output_csv=stats_csv,
                                                          cleaned_output_csv=cleaned_csv)

        # obtain group resolution and formatted display dataframe for plotting
        group_names, group_col_names, groups_dict = get_group_col_names(df, sample_groups_tsv)
        results = run_statistics(df, group_col_names)
        _, display_df = results_to_dataframe(results, groups_dict)

        # use index labels for plotting
        metabolic_labels = df.index.tolist()

        # avoid changing original display_df
        display_df = display_df.copy()

        os.makedirs(os.path.dirname(donut_png), exist_ok=True)
        plot_donut(
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
        plot_table(display_df, output_path=table_png)

        combined_png = os.path.join(os.path.dirname(donut_png), 'combined_donut_table.png')
        combine_images_side_by_side(donut_png, table_png, combined_png, padding=24, match_height='left')

    finally:
        #remove temporary groups file if it was created
        if created_temp_groups and temp_groups_path and os.path.exists(temp_groups_path):
            os.remove(temp_groups_path)
            

if __name__ == '__main__':
    statNut_run()
