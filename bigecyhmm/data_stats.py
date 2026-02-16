import os
import numpy as np
import pandas as pd
from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests
from dataclasses import dataclass
from typing import List, Tuple, Optional


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


#create list of column names per "analysis-group" as defined in the groups TSV. Also return group names and dict for later use.
def get_group_col_names(df: pd.DataFrame, groups_tsv: str) -> Tuple[List[str], List[List[str]], dict]:

    mapping_df = pd.read_csv(groups_tsv, sep='\t', dtype=str)
    # trim whitespace to avoid accidental mismatches
    if 'sample' not in mapping_df.columns or 'group' not in mapping_df.columns:
        raise ValueError("Group mapping TSV must contain columns 'sample' and 'group'")

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





if __name__ == '__main__':
    print('This module provides data-processing utilities. Run the pipeline from main.py')
