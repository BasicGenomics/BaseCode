import argparse
import polars as pl
import pandas as pd
import yaml
import pyfastx
from collections import Counter

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
                raise NotImplementedError
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
                raise NotImplementedError
            sample_barcode_list_PCR_B.append(seq)
    return ';'.join(sample_barcode_list_PCR_B)


def make_sample_barcode_to_sample_id_map(df):
    map_dict = {}
    for idx, row in df.iterrows():
        for barcode in row['SAMPLE_BARCODES_PCR_A'].split(';'):
            if barcode in map_dict:
                sample_name = map_dict[barcode]
                if sample_name != row['SAMPLE_ID']:
                    raise NotImplementedError
            else:
                map_dict[barcode] = row['SAMPLE_ID']
        for barcode in row['SAMPLE_BARCODES_PCR_B'].split(';'):
            if barcode in map_dict:
                sample_name = map_dict[barcode]
                if sample_name != row['SAMPLE_ID']:
                    raise NotImplementedError
            else:
                map_dict[barcode] = row['SAMPLE_ID']
    return map_dict

def make_sample_barcode_list(barcode_set):
    l = []
    for barcodes in barcode_set:
        for barcode in barcodes.split(';'):
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
                map_dict[barcode] = '5prime_int'
        for barcodes in row['SAMPLE_BARCODES_PCR_B'].split(';'):
            for barcode in barcodes.split(';'):
                map_dict[barcode] = '3prime'
    return map_dict

def get_forward_and_reverse_sequences(index_sequence_map):
    index1_sequences_forward = []
    index2_sequences_forward = []
    index1_sequences_reverse = []
    index2_sequences_reverse = []
    for name, barcodes in index_sequence_map.items():
        if 'Fw' in name:
            for barcode in barcodes.split(';'):
                index1_sequences_forward.append(barcode)
                index1_sequences_reverse.append(reverse_complement(barcode))
        elif 'Rv' in name:
            for barcode in barcodes.split(';'):
                index2_sequences_forward.append(barcode)
                index2_sequences_reverse.append(reverse_complement(barcode))
    return index1_sequences_forward, index1_sequences_reverse, index2_sequences_forward, index2_sequences_reverse

def scan_fastq_for_order_and_orientation(fastq_file, index_sequence_map):
    index1_sequences_forward, index1_sequences_reverse, index2_sequences_forward, index2_sequences_reverse = get_forward_and_reverse_sequences(index_sequence_map)

    first = Counter({'index1_forward': 0, 'index1_reverse': 0, 'index2_forward': 0, 'index2_reverse': 0})
    second = Counter({'index1_forward': 0, 'index1_reverse': 0, 'index2_forward': 0, 'index2_reverse': 0})
    i = 0
    for (name, seq, qual) in pyfastx.Fastq(fastq_file, build_index=False):
        index_seq = seq[-16:]
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

    print(first, second)
    
    if 'index1' == first_most_common.split('_')[0]:
        seq_order += '1_'
    else:
        seq_order += '2_'
    
    if 'index1' == second_most_common.split('_')[0]:
        seq_order += '1'
    else:
        seq_order += '2'

    if seq_order == '1_1' or seq_order == '2_2':
        raise RuntimeError
    
    forward_list = ['forward' in first_most_common, 'forward' in second_most_common]
    
    return seq_order, forward_list


def main():
    parser = argparse.ArgumentParser(description='', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-s','--samplesheet',metavar='samplesheet', type=str, help='Samplesheet')
    parser.add_argument('--index-sequences',metavar='indexes', type=str, help='Index Sequences')
    parser.add_argument('--fastq', metavar='fastq', type=str, help='fastq file with index sequences')
    parser.add_argument('--sample-barcodes', metavar='barcodes', type=str, help='Sample barcode output')
    parser.add_argument('--cell-barcodes', metavar='barcodes', type=str, help='Cell barcode output')
    parser.add_argument('--sample-map', metavar='samples', type=str, help='Sample map output')
    parser.add_argument('--readtype-map', metavar='readtypes', type=str, help='readtype map output')
    parser.add_argument('--samplesheet-out', metavar='samplesheet_out', type=str, help='Samplesheet output')
    args = parser.parse_args()

    samplesheet_file = args.samplesheet
    index_sequences = args.index_sequences
    fastq_file = args.fastq
    sample_barcodes_file = args.sample_barcodes
    cell_barcodes_file = args.cell_barcodes
    sample_map_file = args.sample_map
    readtype_map_file = args.readtype_map
    samplesheet_out = args.samplesheet_out

    with open(index_sequences, 'r') as f_is:
        index_sequence_map = yaml.safe_load(f_is)

    seq_order, forward_list = scan_fastq_for_order_and_orientation(fastq_file, index_sequence_map)

    df_excel = pl.read_excel(samplesheet_file)
    df_excel[0, 'PCR A'] = 'FW1_PCR_A'
    df_excel[0, 'PCR B'] = 'FW1_PCR_B'
    df_excel.columns = df_excel.iter_rows().__next__()
    df_excel = df_excel.filter(pl.col("RV1") != "RV1")

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
        f.writelines("\n".join(make_sample_barcode_list(set(samplesheet_df["SAMPLE_BARCODES_PCR_B"]))))

    samplesheet_df['BARCODE'] = 'nan'
    with open(cell_barcodes_file, 'w') as f:
        f.writelines("\n".join(make_sample_barcode_list(set(samplesheet_df['BARCODE']))))

    samplesheet_df['SAMPLE_ID'] = samplesheet_df.apply(fix_SAMPLE_ID, axis=1)
    samplesheet_df.index = samplesheet_df.apply(lambda row: row['SAMPLE_ID']+row['BARCODE'], axis=1) 

    map_dict = make_sample_barcode_to_sample_id_map(samplesheet_df)

    with open(sample_map_file, 'w') as f:
        yaml.dump(map_dict, f)
    
    with open(readtype_map_file, 'w') as f:
        yaml.dump(sample_barcode_to_readType_map, f)
    
    samplesheet_df.to_csv(samplesheet_out)
if __name__ == "__main__":
    main()
