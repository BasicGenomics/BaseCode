import polars as pl
import pandas as pd
import argparse
import yaml
import tabulate
import json
from collections import Counter
import numpy as np


def main():
    parser = argparse.ArgumentParser(description='Calculate summary statistics', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-j','--json',metavar='json', type=str, help='JSON read count file')
    parser.add_argument('-l','--long-form', metavar='long', type=str, help='Long form file')
    parser.add_argument('-s','--sample-map', metavar='sample map', type=str, help='Sample map file')
    parser.add_argument('-o','--output', help='Output summary file (csv)')

    args = parser.parse_args()

    json_file = args.json
    long_form_file = args.long_form
    sample_map_file = args.sample_map
    output_file = args.output

    with open(json_file, 'r') as f:
        read_flow_dict = json.load(f)
    
    with open(sample_map_file, 'r') as f:
        sample_map_dict = yaml.safe_load(f)
    
    df_long_form = pl.read_csv(long_form_file)

    tags_to_column_names = {'TP_read': '3\' Counts', 'internal': 'Internal Counts', 'FP_read': '5\' Counts'}
    
    sample_counts = {}
    for bc, xx_stats in read_flow_dict['XX_stats'].items():
        if bc in sample_map_dict:
            sample = sample_map_dict[bc]
            if sample not in sample_counts:
                sample_counts[sample] = Counter()
            sample_counts[sample].update(xx_stats)
        else:
            continue
    
    total_counts_dict =  {sample: sum(c.values()) for sample, c in sample_counts.items()}



    df_sample_counts_pandas = pd.DataFrame(sample_counts).T
    df_sample_counts_pandas.columns = [tags_to_column_names[c] for c in df_sample_counts_pandas.columns]

    series_total_counts = pd.Series(total_counts_dict)
    series_total_counts.name = 'Total Counts'

    df_sample_counts_pandas = df_sample_counts_pandas.join(series_total_counts)

    df_sample_counts_polars = pl.from_pandas(df_sample_counts_pandas.reset_index())

    sample_list = []
    genes_detected_list = []
    molecules_detected_TP_list = []
    molecules_detected_INT_list = []
    molecules_detected_FP_list = []
    percentage_reconstructed_list = []
    median_reconstructed_completed_list = []

    for (unique_sample_id, df_long_form_sample) in df_long_form.group_by('SM'):
        df_long_form_completed = df_long_form_sample.filter(((pl.col('TC') > 0) & (pl.col('IC') > 0) & (pl.col('FC') > 0)))
        df_long_form_sample_threep = df_long_form_sample.filter(pl.col('TC') > 0)
        df_long_form_sample_int = df_long_form_sample.filter(pl.col('IC') > 0)
        df_long_form_sample_fivep = df_long_form_sample.filter(pl.col('FC') > 0)

        sample_list.extend(unique_sample_id)

        genes_detected_sample = df_long_form_sample.unique(subset=['XT']).shape[0]
        genes_detected_list.append(genes_detected_sample)

        molecules_detected_TP_sample = df_long_form_sample_threep.shape[0]
        molecules_detected_TP_list.append(molecules_detected_TP_sample)

        molecules_detected_INT_sample = df_long_form_sample_int.shape[0]
        molecules_detected_INT_list.append(molecules_detected_INT_sample)

        molecules_detected_FP_sample = df_long_form_sample_fivep.shape[0]
        molecules_detected_FP_list.append(molecules_detected_FP_sample)

        percentage_reconstructed_sample = np.round(100*(df_long_form_completed.shape[0]/df_long_form_sample_threep.shape[0]),2)
        percentage_reconstructed_list.append(percentage_reconstructed_sample)
        
        median_reconstructed_completed_sample = df_long_form_completed.median()['QL'][0]
        median_reconstructed_completed_list.append(median_reconstructed_completed_sample)

    data = {'index': sample_list, 'Genes Detected': genes_detected_list, 'Molecules Detected (5\')': molecules_detected_FP_list, 'Molecules Detected (INT)': molecules_detected_INT_list, 'Molecules Detected (3\')': molecules_detected_TP_list, 'Percentage Completed (%)': percentage_reconstructed_list, 'Median Length Completed Molecules (bp)': median_reconstructed_completed_list }

    df_reconstruction_polars = pl.DataFrame(data)

    df_full = df_sample_counts_polars.join(df_reconstruction_polars, on = 'index')

    df_full = df_full.sort('Total Counts', descending=True)
    
    df_full_head_10 = df_full.head(10)

    values = [list(s) for s in df_full_head_10.rows()]

    header = df_full_head_10.columns
    
    print(tabulate.tabulate(values, headers=header))

    df_full.write_csv(output_file) 
    print('Full summary report for all samples found at {}'.format(output_file))

if __name__ == "__main__":
    main()
