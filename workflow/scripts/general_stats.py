import pandas as pd
import pysam
import seaborn as sns
import matplotlib.pyplot as plt
from collections import Counter
import argparse

# Exon, Intron, Ambiguous, Intergenic, Unmapped
def get_mapping_group(read):
    if read.is_unmapped:
        return 'Unmapped'
    ES_tag = read.get_tag('ES')
    IS_tag = read.get_tag('IS')
    if 'Assigned' in ES_tag:
        return 'Exon'
    if 'Assigned' in IS_tag:
        return 'Intron'
    if 'Unassigned' in ES_tag and 'Unassigned' in IS_tag:
        return 'Intergenic'
    return 'Ambiguous'

def get_read_type_stats(read_type_dicts):
    read_type_stats = {}
    for read_type, read_type_dict in read_type_dicts.items():
        if read_type == 'TP_read':
            n = 0
            for read_pair, read_pair_dict in read_type_dict.items():
                for cell, cell_counter in read_pair_dict.items():
                    n += sum(cell_counter.values())
                read_type_stats[read_type] = n
        else:
            read_type_stats[read_type] = 0
            for read_pair, read_pair_dict in read_type_dict.items():
                read_type_stats[read_type] += sum(read_pair_dict.values())
    return read_type_stats

def modify_read_type_stats(read_type_dicts):
    read_type_stats = {}
    for read_type, read_type_dict in read_type_dicts.items():
        if read_type == 'TP_read':
            threep_counter = Counter()
            for cell, cell_counter in read_type_dict.items():
                threep_counter.update(cell_counter)
            read_type_stats[read_type] = threep_counter
        else:
            read_type_stats[read_type] = read_type_dict
    return read_type_stats

def get_mapping_group_dict(bamfile, cell_set, sample_set):
    mapping_group_dict = {}

    bam = pysam.AlignmentFile(bamfile)

    for read in bam.fetch(until_eof=True):
        if read.has_tag('SB'):
            sample = read.get_tag('SB')
        else:
            sample = read.get_tag('SM')
        if sample not in sample_set:
            sample = 'UNASSIGNED'
            
        if sample not in mapping_group_dict:
            mapping_group_dict[sample] = {}
        
        mapping_group = get_mapping_group(read)
        read_type = read.get_tag('XX')
        read_pair = 'read1' if read.is_read1 else 'read2'
        
        if read_type not in mapping_group_dict[sample]:
            mapping_group_dict[sample][read_type] = {}
            
        if read_type == 'TP_read':
            if read_pair not in  mapping_group_dict[sample][read_type]:
                 mapping_group_dict[sample][read_type][read_pair] = {}
            cell = read.get_tag('BC')
            if cell not in cell_set:
                cell = 'UNASSIGNED'
            if cell not in mapping_group_dict[sample][read_type][read_pair]:
                mapping_group_dict[sample][read_type][read_pair][cell] = Counter()
            mapping_group_dict[sample][read_type][read_pair][cell].update({mapping_group:1})
        else:
            if read_pair not in  mapping_group_dict[sample][read_type]:
                 mapping_group_dict[sample][read_type][read_pair] = Counter()
            mapping_group_dict[sample][read_type][read_pair].update({mapping_group:1})
    bam.close()
    return mapping_group_dict


def main():
    parser = argparse.ArgumentParser(description='', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--input',metavar='input', type=str, help='Input .bam file')
    parser.add_argument('-c','--cells',metavar='cells', type=str, help='Cell barcodes .txt file' )
    parser.add_argument('-s','--samplesheet', help='Samplesheet file')
    parser.add_argument('-o','--out_dir', help='Output directory')
    parser.add_argument('-p', '--prefix')

    args = parser.parse_args()

    cell_set = set([line.rstrip() for line in open(args.cells)])
    sample_set = set([line.rstrip() for line in open(args.samplesheet)])

    mapping_group_dict = get_mapping_group_dict(args.input, cell_set, sample_set)

    

    

    sample_dfs = {}
    for sample, sample_dict in mapping_group_dict.items():
        if 'TP_read' in sample_dict:
            for read_pair, read_pair_dict in sample_dict['TP_read'].items():
                df = pd.DataFrame(read_pair_dict).T
                df['sample_barcode'] = sample
                df['read_pair'] = read_pair
                sample_dfs['{}_{}'.format(sample, read_pair)] = df
    
    df_cell_level = pd.concat(sample_dfs,axis=0)
    df_cell_level.index = ['_'.join([s1,s2]) for (s1,s2) in df_cell_level.index]

    df_cell_level.to_csv('{}/{}_mapping_categories_per_sample.csv'.format(args.out_dir, args.prefix))

    sample_nonthreep_dfs = {}
    for sample, sample_dict in mapping_group_dict.items():
        res_dict = {}
        for read_type, read_type_dict in sample_dict.items():
            if read_type == 'TP_read':
                continue
            for read_pair, counter in read_type_dict.items():
                if read_pair not in res_dict:
                    res_dict[read_pair] = {}
                
                res_dict[read_pair][read_type] = counter
        for read_pair, read_type_dict in res_dict.items():
            df = pd.DataFrame(read_type_dict).T
            df['sample_barcode'] = sample
            df['read_pair'] = read_pair
            sample_nonthreep_dfs['{}_{}'.format(sample, read_pair)] = df

    df_nonbarcoded = pd.concat(sample_nonthreep_dfs.values(), axis=0)

    df_nonbarcoded.to_csv('{}/{}_nonbarcoded_mapping_categories_per_sample.csv'.format(args.out_dir, args.prefix))

    read_type_per_sample_df = pd.DataFrame({sample: get_read_type_stats(sample_dict) for sample, sample_dict in mapping_group_dict.items()}).T

    read_type_per_sample_df.to_csv('{}/{}_read_type_per_sample.csv'.format(args.out_dir, args.prefix))

if __name__ == '__main__':
    main()
