import argparse
import polars as pl
import pandas as pd
import yaml
import pyfastx
from collections import Counter
import os

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}

def reverse_complement(seq):
    return ''.join([complement[s] for s in seq][::-1])

def make_sample_barcode_PCR_A(row, seq_order, forward_list):
    FW1_PCR_A_seqs = row['FW1_PCR_A_SEQ'].split(';')
    RV1_seqs = row['RV1_SEQ'].split(';')
    sample_barcode_list_PCR_A = []
    for index1_seq in FW1_PCR_A_seqs:
        for index2_seq in RV1_seqs:
            if seq_order == '1_2':
                seq = ''
                if forward_list[0]:
                    seq += index1_seq
                else:
                    seq += reverse_complement(index1_seq)
                
                if forward_list[1]:
                    seq += index2_seq
                else:
                    seq += reverse_complement(index2_seq)
            elif seq_order == '2_1':
                seq = ''
                if forward_list[0]:
                    seq += index2_seq
                else:
                    seq += reverse_complement(index2_seq)
                
                if forward_list[1]:
                    seq += index1_seq
                else:
                    seq += reverse_complement(index1_seq)
            else:
                raise ValueError(f"Unexpected seq_order: {seq_order!r}. Expected '1_2' or '2_1'.")
            sample_barcode_list_PCR_A.append(seq)
    return ';'.join(sample_barcode_list_PCR_A)

def make_sample_barcode_PCR_B(row, seq_order, forward_list):
    FW1_PCR_B_seqs = row['FW1_PCR_B_SEQ'].split(';')
    RV2_seqs = row['RV2_SEQ'].split(';')
    sample_barcode_list_PCR_B = []
    for index1_seq in FW1_PCR_B_seqs:
        for index2_seq in RV2_seqs:
            if seq_order == '1_2':
                seq = ''
                if forward_list[0]:
                    seq += index1_seq
                else:
                    seq += reverse_complement(index1_seq)
                
                if forward_list[1]:
                    seq += index2_seq
                else:
                    seq += reverse_complement(index2_seq)
            elif seq_order == '2_1':
                seq = ''
                if forward_list[0]:
                    seq += index2_seq
                else:
                    seq += reverse_complement(index2_seq)
                
                if forward_list[1]:
                    seq += index1_seq
                else:
                    seq += reverse_complement(index1_seq)
            else:
                raise ValueError(f"Unexpected seq_order: {seq_order!r}. Expected '1_2' or '2_1'.")
            sample_barcode_list_PCR_B.append(seq)
    return ';'.join(sample_barcode_list_PCR_B)

def make_sample_barcode_to_sample_id_map(df):
    map_dict = {}
    for idx, row in df.iterrows():
        for barcode in row['SAMPLE_BARCODES_PCR_A'].split(';'):
            if barcode == "":
                continue
            if barcode in map_dict:
                sample_name = map_dict[barcode]
                if sample_name != row['SAMPLE_ID']:
                    raise ValueError(
                        f"Barcode collision: '{barcode}' is assigned to both "
                        f"'{sample_name}' and '{row['SAMPLE_ID']}'. "
                        f"Barcodes must be unique across samples."
                    )
            else:
                map_dict[barcode] = row['SAMPLE_ID']
        for barcode in row['SAMPLE_BARCODES_PCR_B'].split(';'):
            if barcode == "":
                continue
            if barcode in map_dict:
                sample_name = map_dict[barcode]
                if sample_name != row['SAMPLE_ID']:
                    raise ValueError(
                        f"Barcode collision: '{barcode}' is assigned to both "
                        f"'{sample_name}' and '{row['SAMPLE_ID']}'. "
                        f"Barcodes must be unique across samples."
                    )
            else:
                map_dict[barcode] = row['SAMPLE_ID']
    return map_dict

def make_sample_barcode_list(barcode_set):
    l = []
    for barcodes in barcode_set:
        for barcode in barcodes.split(';'):
            if barcode == "":
                continue
            l.append(barcode)
    return list(set(l))

def fix_SAMPLE_ID(row):
    s = row['SAMPLE_ID']
    s = s.replace(" ", "_")
    s = s.replace("-", "_")
    return s

def make_sample_barcode_to_readType_map(df):
    map_dict = {}
    for idx, row in df.iterrows():
        for barcodes in row['SAMPLE_BARCODES_PCR_A'].split(';'):
            for barcode in barcodes.split(';'):
                if barcode == "":
                    continue
                if barcode in map_dict and map_dict[barcode] != '5prime_int':
                    raise ValueError(
                        f"Barcode {barcode} assigned conflicting read types "
                        f"({map_dict[barcode]} vs 5prime_int)"
                    )
                map_dict[barcode] = '5prime_int'
        for barcodes in row['SAMPLE_BARCODES_PCR_B'].split(';'):
            for barcode in barcodes.split(';'):
                if barcode == "":
                    continue
                if barcode in map_dict and map_dict[barcode] != '3prime':
                    raise ValueError(
                        f"Barcode {barcode} assigned conflicting read types "
                        f"({map_dict[barcode]} vs 3prime)"
                    )
                map_dict[barcode] = '3prime'
    return map_dict

def get_forward_and_reverse_sequences(index_sequence_map):
    index1_sequences_forward = []
    index2_sequences_forward = []
    index1_sequences_reverse = []
    index2_sequences_reverse = []
    for name, barcodes in index_sequence_map.items():
        if barcodes is None or name is None:
            continue
        if 'Fw' in name:
            for barcode in barcodes.split(';'):
                index1_sequences_forward.append(barcode)
                index1_sequences_reverse.append(reverse_complement(barcode))
        elif 'Rv' in name:
            for barcode in barcodes.split(';'):
                index2_sequences_forward.append(barcode)
                index2_sequences_reverse.append(reverse_complement(barcode))
    return index1_sequences_forward, index1_sequences_reverse, index2_sequences_forward, index2_sequences_reverse

def fastq_iteration(fastq_files: list):
    if len(fastq_files) == 1:
        fastq_file = pyfastx.Fastq(fastq_files[0], build_index=False)
        for t in fastq_file:
            yield t
    elif len(fastq_files) == 2:
        fastq_file_i1 = pyfastx.Fastq(fastq_files[0], build_index=False)
        fastq_file_i2 = pyfastx.Fastq(fastq_files[1], build_index=False)
        for i1, i2 in zip(fastq_file_i1, fastq_file_i2):
            yield (i1[0], i1[1]+i2[1], i1[2]+i2[2])
    else:
        raise ValueError(f"Expected 1 or 2 FASTQ files, got {len(fastq_files)}.")

def scan_fastq_for_order_and_orientation(fastq_files: list, index_sequence_map: dict):
    index1_sequences_forward, index1_sequences_reverse, index2_sequences_forward, index2_sequences_reverse = get_forward_and_reverse_sequences(index_sequence_map)
    first = Counter({'index1_forward': 0, 'index1_reverse': 0, 'index2_forward': 0, 'index2_reverse': 0})
    second = Counter({'index1_forward': 0, 'index1_reverse': 0, 'index2_forward': 0, 'index2_reverse': 0})
    i = 0
    for (name, seq, qual) in fastq_iteration(fastq_files):
        if len(fastq_files) == 1:
            index_seq = seq[-16:]
        elif len(fastq_files) == 2:
            index_seq = seq
        else:
            raise ValueError(f"Expected 1 or 2 FASTQ files, got {len(fastq_files)}.")
        i += 1
        if i > 10000:
            break
        if index_seq[:8] in index1_sequences_forward:
            first['index1_forward'] += 1
        if index_seq[:8] in index1_sequences_reverse:
            first['index1_reverse'] += 1
        if index_seq[:8] in index2_sequences_forward:
            first['index2_forward'] += 1
        if index_seq[:8] in index2_sequences_reverse:
            first['index2_reverse'] += 1
        if index_seq[-8:] in index1_sequences_forward:
            second['index1_forward'] += 1
        if index_seq[-8:] in index1_sequences_reverse:
            second['index1_reverse'] += 1
        if index_seq[-8:] in index2_sequences_forward:
            second['index2_forward'] += 1
        if index_seq[-8:] in index2_sequences_reverse:
            second['index2_reverse'] += 1
    first_most_common = first.most_common()[0][0]
    second_most_common = second.most_common()[0][0]
    seq_order = ''
    if 'index1' == first_most_common.split('_')[0]:
        seq_order += '1_'
    else:
        seq_order += '2_'
    if 'index1' == second_most_common.split('_')[0]:
        seq_order += '1'
    else:
        seq_order += '2'
    if seq_order == '1_1' or seq_order == '2_2':
        raise RuntimeError(
            f"Ambiguous index order: both read positions matched the same index type "
            f"(seq_order: {seq_order!r}). Check that index sequences are correct and "
            f"that the right FASTQ files were provided."
        )
    forward_list = ['forward' in first_most_common, 'forward' in second_most_common]
    return seq_order, forward_list


def main():
    parser = argparse.ArgumentParser(
        description=(
            'Parse a samplesheet (CSV or Excel) and index sequence YAML to produce '
            'sample barcode files, sample/readtype maps, and a processed samplesheet. '
            'Index order and orientation are auto-detected from the provided FASTQ index file(s).'
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '-s', '--samplesheet', metavar='Samplesheet', type=str, required=True,
        help=(
            'Samplesheet file. Accepted formats: '
            '(1) Excel (.xlsx) with a two-row header where row 1 has group labels. '
            '"PCR A" and "PCR B" above the FW1 columns, and row 2 has: '
            'SAMPLE_ID, FW1, RV1, FW1, RV2, DESC. '
            '(2) CSV with columns: SAMPLE_ID, FW1_PCR_A, RV1, FW1_PCR_B, RV2, DESC.'
        )
    )
    parser.add_argument(
        '--index-sequences', metavar='Index Sequences', type=str, required=True,
        help=(
            'YAML file mapping index names to semicolon-separated sequences.'
        )
    )
    parser.add_argument(
        '--fastq', metavar='FASTQ files', type=str, nargs='+', required=True,
        help=(
            'One or two FASTQ files containing index sequences.'
        )
    )
    parser.add_argument(
        '--sample-barcodes', metavar='Sample Barcodes', type=str,
        help='Output text file with all assembled sample barcodes, one per line.'
    )
    parser.add_argument(
        '--cell-barcodes', metavar='Cell Barcodes', type=str,
        help='Output text file with cell barcodes, one per line.'
    )
    parser.add_argument(
        '--sample-map', metavar='Sample Map', type=str,
        help='Output YAML file mapping each sample barcode to its SAMPLE_ID.'
    )
    parser.add_argument(
        '--readtype-map', metavar='Readtype Map', type=str,
        help='Output YAML file mapping each sample barcode to its read type (5prime_int or 3prime).'
    )
    parser.add_argument(
        '--samplesheet-out', metavar='Samplesheet with Barcodes', type=str,
        help='Output CSV file with the processed samplesheet including all computed barcode columns.'
    )
    parser.add_argument(
        '--ignore-none', action='store_true',
        help='Treat empty/None cells in the samplesheet as empty barcode sequences instead of raising an error.'
    )
    args = parser.parse_args()

    samplesheet_file = args.samplesheet
    index_sequences = args.index_sequences
    fastq_files = list(args.fastq)
    sample_barcodes_file = args.sample_barcodes
    cell_barcodes_file = args.cell_barcodes
    sample_map_file = args.sample_map
    readtype_map_file = args.readtype_map
    samplesheet_out = args.samplesheet_out
    ignore_none = args.ignore_none

    with open(index_sequences, 'r') as f_is:
        index_sequence_map = yaml.safe_load(f_is)

    if ignore_none:
        index_sequence_map[None] = ""

    seq_order, forward_list = scan_fastq_for_order_and_orientation(fastq_files, index_sequence_map)

    _filename, file_extension = os.path.splitext(samplesheet_file)

    if file_extension == '.csv':
        samplesheet_df = pd.read_csv(samplesheet_file, index_col=0)
    else:
        df_excel = pl.read_excel(samplesheet_file)
        df_excel[0, 'PCR A'] = 'FW1_PCR_A'
        df_excel[0, 'PCR B'] = 'FW1_PCR_B'
        df_excel.columns = df_excel.iter_rows().__next__()
        df_excel = df_excel.remove(pl.col("RV1") == "RV1")
        samplesheet_df = df_excel.to_pandas()

    samplesheet_df["FW1_PCR_A_SEQ"] = samplesheet_df.apply(lambda row: index_sequence_map[row['FW1_PCR_A']], axis=1)
    samplesheet_df["RV1_SEQ"] = samplesheet_df.apply(lambda row: index_sequence_map[row['RV1']], axis=1)
    samplesheet_df["FW1_PCR_B_SEQ"] = samplesheet_df.apply(lambda row: index_sequence_map[row['FW1_PCR_B']], axis=1)
    samplesheet_df["RV2_SEQ"] = samplesheet_df.apply(lambda row: index_sequence_map[row['RV2']], axis=1)
    samplesheet_df["SAMPLE_BARCODES_PCR_A"] = samplesheet_df.apply(lambda row: make_sample_barcode_PCR_A(row, seq_order, forward_list), axis=1)
    samplesheet_df["SAMPLE_BARCODES_PCR_B"] = samplesheet_df.apply(lambda row: make_sample_barcode_PCR_B(row, seq_order, forward_list), axis=1)
    
    sample_barcode_to_readType_map = make_sample_barcode_to_readType_map(samplesheet_df)

    with open(sample_barcodes_file, 'w') as f:
        f.writelines("\n".join(make_sample_barcode_list(set(samplesheet_df["SAMPLE_BARCODES_PCR_A"]))))
        f.write('\n')
        f.writelines("\n".join(make_sample_barcode_list(set(samplesheet_df["SAMPLE_BARCODES_PCR_B"]))))
    
    samplesheet_df['BARCODE'] = 'nan'

    with open(cell_barcodes_file, 'w') as f:
        f.writelines("\n".join(make_sample_barcode_list(set(samplesheet_df['BARCODE']))))

    samplesheet_df['SAMPLE_ID'] = samplesheet_df.apply(fix_SAMPLE_ID, axis=1)
    samplesheet_df.index = samplesheet_df.apply(lambda row: row['SAMPLE_ID']+"_"+row['BARCODE'], axis=1) 

    map_dict = make_sample_barcode_to_sample_id_map(samplesheet_df)

    with open(sample_map_file, 'w') as f:
        yaml.dump(map_dict, f)
    
    with open(readtype_map_file, 'w') as f:
        yaml.dump(sample_barcode_to_readType_map, f)
    
    samplesheet_df.to_csv(samplesheet_out)
if __name__ == "__main__":
    main()
