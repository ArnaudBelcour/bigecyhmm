import os
import csv
import tempfile
import pandas as pd
import numpy as np
import logging

from dataclasses import dataclass
from typing import List, Tuple, Optional

from bigecyhmm.group_plot import plot_donut, plot_table, combine_images_side_by_side, combine_images_side_by_side
from bigecyhmm import PATHWAY_TEMPLATE_FILE

from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests

#Default file paths and groups. 
DEFAULT_INPUT_TSV = os.path.join('function_abundance', 'cycle_abundance_sample.tsv')
DEFAULT_STATS_TSV = os.path.join('plots', 'group_stats.tsv')
DEFAULT_CLEANED_TSV = os.path.join('plots', 'cleaned_data.tsv')
DEFAULT_DONUT_PNG = os.path.join('plots', 'group_medians_donut.png')
DEFAULT_TABLE_PNG = os.path.join('plots', 'group_stats_table.png')
DEFAULT_BACKGROUND_PATH = os.path.join(os.path.dirname(__file__), 'template_central_hydrogen.png') 
#Default sample groups: expect sample_groups.tsv to be located in the pipeline input folder (one level up from visualisation outputs)
DEFAULT_SAMPLE_GROUPS = os.path.join('..', 'sample_groups.tsv')
DEFAULT_GROUPS = {'All Samples': ['']}


#DEFAULTS
ROUND_TOLERANCE = 1e-12  #small tolerance for near-zero checks 
ALPHA = 0.05             #significance level

logger = logging.getLogger()
logger.setLevel(logging.INFO)

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


def get_group_col_names(df: pd.DataFrame, mapping_df: pd.DataFrame) -> Tuple[List[str], List[List[str]], dict]:
    """ Create list of column names per "analysis-group" as defined in the groups TSV. Also return group names and dict for later use.

    Args:
        df (pd.DataFrame): dataframe containing samples as columns and metabolic pathway as row, indicating abundance of metabolic pathway.
        mapping_df (pd.DataFrame): dataframe linking sample to group.

    Returns:
        group_names (list): list of groups
        resolved_groups (list): list of lists, each list is linked to a group and contains the samples associated with the group
        groups_dict (dict): dictionary linking group (keys) to sample (values)
    """
    mapping_df['sample'] = mapping_df['sample'].astype(str).str.strip()
    mapping_df['group'] = mapping_df['group'].astype(str).str.strip()

    groups_dict = mapping_df.groupby('group', sort=False)['sample'].apply(list).to_dict()
    cols = list(df.columns)

    group_names = list(groups_dict.keys())
    resolved_groups: List[List[str]] = []
    for samples in groups_dict.values():
        samples_clean = [str(s).strip() for s in samples if isinstance(s, str)]
        if '' in samples_clean:
            matched = cols.copy()
        else:
            matched = [s for s in samples_clean if s in df.columns]
        # preserve order and remove duplicates
        uniq = list(dict.fromkeys(matched))
        resolved_groups.append(uniq)

    return group_names, resolved_groups, groups_dict


def run_statistics(df: pd.DataFrame, group_col_names: List[List[str]]) -> List[FunctionStatResults]:
    """ Compute Kruskal-Wallis test on abundance of metabolic pathways between sample of different groups, then ran Benjamini–Hochberg correction.

    Args:
        df (pd.DataFrame): dataframe containing samples as columns and metabolic pathway as row, indicating abundance of metabolic pathway
        group_col_names (list): list of lists, each list is linked to a group and contains the samples associated with the group

    Returns:
        results (list): list of FunctionStatResults objects, indicating for each metabolic function the associated stats linked to group comparison
    """
    # Get metabolic function names.
    function_names = df.index.tolist()
    results: List[FunctionStatResults] = []

    # Compute median, count and quantile per-group Series (one Series per group) using pandas aggregations.
    med_series_list = [df[cols].median(axis=1).astype(float) if cols else pd.Series(np.nan, index=df.index) for cols in group_col_names]
    count_series_list = [df[cols].count(axis=1).astype(int) if cols else pd.Series(0, index=df.index) for cols in group_col_names]
    iqr_series_list = [(df[cols].quantile(0.75, axis=1) - df[cols].quantile(0.25, axis=1)).astype(float) if cols else pd.Series(np.nan, index=df.index) for cols in group_col_names]

    percent_zero_series_list: List[pd.Series] = []
    for cols in group_col_names:
        if cols:
            is_zero = df[cols].map(lambda v: np.isclose(v, 0.0, atol=ROUND_TOLERANCE)).sum(axis=1)
            percent_zero_series_list.append((is_zero / len(cols) * 100).astype(float))
        else:
            percent_zero_series_list.append(pd.Series(np.nan, index=df.index))

    for function_name in function_names:
        med_row = [float(s.loc[function_name]) if pd.notna(s.loc[function_name]) else np.nan for s in med_series_list]
        counts = [int(s.loc[function_name]) if pd.notna(s.loc[function_name]) else 0 for s in count_series_list]
        iq_range = [float(s.loc[function_name]) if pd.notna(s.loc[function_name]) else np.nan for s in iqr_series_list]
        percent_zero = [float(s.loc[function_name]) if pd.notna(s.loc[function_name]) else np.nan for s in percent_zero_series_list]

        # Numeric arrays per group for Kruskal
        data_per_group = [df.loc[function_name, cols].dropna().values.astype(float) if len(cols) > 0 else np.array([], dtype=float) for cols in group_col_names]

        # Decide whether Kruskal–Wallis can be computed
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
            name=function_name,
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

    # Benjamini–Hochberg correction
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
    """ Compute Kruskal-Wallis test on abundance of metabolic pathways between sample of different groups, then ran Benjamini–Hochberg correction.

    Args:
        results (list): list of FunctionStatResults objects, indicating for each metabolic function the associated stats linked to group comparison
        groups (dict): dictionary linking group (keys) to sample (values)

    Returns:
        df_res (pd.DataFrame): DataFrame containing statistical results as columns and metabolic pathway as rows
        display_df (pd.DataFrame): DataFrame containing statistical results as columns and metabolic pathway as rows, modified to display
    """
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
    display_df['p'] = display_df['p'].apply(lambda x: f"{x:.4g}" if pd.notna(x) else "n/a")
    display_df['p_bh'] = display_df['p_bh'].apply(lambda x: f"{x:.4g}" if pd.notna(x) else "n/a")
    display_df['significant'] = display_df['significant'].apply(lambda x: 'Yes' if x else 'No')

    return df_res, display_df


def compute_and_save(input_df: pd.DataFrame, mapping_df: pd.DataFrame, output_csv: str, cleaned_output_csv: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """ Retrieve gorups, compute statistics, generate stat output file.

    Args:
        input_df (pd.DataFrame): dataframe containing samples as columns and metabolic pathway as row, indicating abundance of metabolic pathway
        mapping_df (pd.DataFrame): dataframe linking sample to group
        output_stats_csv (str): path to output file containing statistical results
        cleaned_output_csv (str): path to cleaned file containg abundance of metabolic pathways in samples

    Returns:
        cleaned_df (pd.DataFrame): cleaned DataFrame containg abundance of metabolic pathways in samples
        group_col_names (list): list of lists, each list is linked to a group and contains the samples associated with the group
        groups_dict (dict): dictionary linking group (keys) to sample (values)
        display_df (pd.DataFrame): DataFrame containing statistical results as columns and metabolic pathway as rows, modified to display
    """
    group_names, group_col_names, groups_dict = get_group_col_names(input_df, mapping_df)

    cleaned_df = input_df.copy()
    os.makedirs(os.path.dirname(cleaned_output_csv), exist_ok=True)
    cleaned_df.to_csv(cleaned_output_csv)
    logger.info(f"Saved cleaned data to {cleaned_output_csv}")

    results = run_statistics(input_df, group_col_names)
    df_res, display_df = results_to_dataframe(results, groups_dict)

    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    df_res.to_csv(output_csv, index=False)
    logger.info(f"Saved numeric stats to {output_csv}")

    return cleaned_df, group_col_names, groups_dict, display_df


def statNut_run(input_tsv: str = DEFAULT_INPUT_TSV,
         sample_groups_tsv: str = DEFAULT_SAMPLE_GROUPS,
         output_stats_csv: str = DEFAULT_STATS_TSV,
         output_cleaned_csv: str = DEFAULT_CLEANED_TSV,
         output_donut_png: str = DEFAULT_DONUT_PNG,
         output_table_png: str = DEFAULT_TABLE_PNG,
         background_path: str = DEFAULT_BACKGROUND_PATH,
         background_offset: tuple = (0.016, -0.006),
         background_scale: float = 0.60):
    """ From a tsv file showing abundance of metabolic pathways in samples, another file linking sampels to group and a background file, generate a donut plot.

    Args:
        input_tsv (str): path to tsv file containing sampels as column and metabolic pathways as row indicating the abundance of pathway in samples
        sample_groups_tsv (str): path to group file indicating for each sample their associated group
        output_stats_csv (str): path to output file containing statistical results
        output_cleaned_csv (str): path to cleaned file containg abundance of metabolic pathways in samples
        output_donut_png (str): path to output donut plot png
        output_donut_png (str): path to output table png
        background_path (str): path to background image for donut plot
        background_offset (tuple): adjust background image position (right, up)
        background_scale (float): adjust background image scale relative to donut
    """
    # Load input tsv, first column (metabo. funct. name) is index.
    df = pd.read_csv(input_tsv, sep='\t', index_col=0)

    # If the provided sample_groups_tsv does not exist, create a temporary one which assigns all samples to a single group called "all samples". 
    temp_groups_path = None
    created_temp_groups = False
    if not sample_groups_tsv or not os.path.exists(sample_groups_tsv):
        # create a temporary TSV that instructs the pipeline to select all samples
        with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.tsv') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerow(['sample', 'group'])
            for sample in df.columns:
                writer.writerow([sample, 'All Samples'])
            temp_groups_path = csvfile.name
        sample_groups_tsv = temp_groups_path
        created_temp_groups = True

    # Load group tsv file.
    mapping_df = pd.read_csv(sample_groups_tsv, sep='\t', dtype=str)
    # Trim whitespace to avoid accidental mismatches.
    if 'sample' not in mapping_df.columns or 'group' not in mapping_df.columns:
        raise ValueError("Group mapping TSV must contain columns 'sample' and 'group'")

    # Compute statistics and save cleaned + rounded CSVs for later plotting.
    cleaned_df, group_col_names, groups_dict, display_df = compute_and_save(input_df=df,
                                            mapping_df=mapping_df,
                                            output_csv=output_stats_csv,
                                            cleaned_output_csv=output_cleaned_csv)

    # Use index labels for plotting.
    metabolic_labels = df.index.tolist()
    all_bigecyhmm_template_cycles = pd.read_csv(PATHWAY_TEMPLATE_FILE, sep='\t')['Pathways'].tolist()
    if set(all_bigecyhmm_template_cycles).issubset(set(metabolic_labels)):
        metabolic_labels = [label.split(':')[1] for label in metabolic_labels]

    # Avoid changing original display_df.
    display_df = display_df.copy()

    os.makedirs(os.path.dirname(output_donut_png), exist_ok=True)
    plot_donut(cleaned_df, groups_dict,
        metabolic_labels=metabolic_labels,
        group_col_names=group_col_names,
        output_path=output_donut_png,
        background_path=background_path,
        background_offset=background_offset,
        background_scale=background_scale)

    os.makedirs(os.path.dirname(output_table_png), exist_ok=True)
    plot_table(display_df, output_path=output_table_png)

    combined_png = os.path.join(os.path.dirname(output_donut_png), 'combined_donut_table.png')
    combine_images_side_by_side(output_donut_png, output_table_png, combined_png, padding=24, match_height='left')

    #remove temporary groups file if it was created
    if created_temp_groups and temp_groups_path and os.path.exists(temp_groups_path):
        os.remove(temp_groups_path)
            

if __name__ == '__main__':
    statNut_run()
