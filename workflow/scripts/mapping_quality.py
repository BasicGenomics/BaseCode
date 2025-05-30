import pandas as pd
import pysam
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
            for cell, cell_counter in read_type_dict.items():
                n += sum(cell_counter.values())
            read_type_stats[read_type] = n
        else:
            read_type_stats[read_type] = sum(read_type_dict.values())
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

def get_quality_dict(bamfile, cell_set, sample_set):
    quality_dict = {}

    bam = pysam.AlignmentFile(bamfile)

    for read in bam.fetch(until_eof=True):
        if read.has_tag('SB'):
            sample = read.get_tag('SB')
        else:
            sample = read.get_tag('SM')
        if sample not in sample_set:
            sample = 'UNASSIGNED'
            
        if sample not in quality_dict:
            quality_dict[sample] = {}
        
        quality_group = read.mapping_quality
        read_type = read.get_tag('XX')
        
        if read_type not in quality_dict[sample]:
            if read_type == 'TP_read':
                quality_dict[sample][read_type] = {}
            else:
                quality_dict[sample][read_type] = Counter()
            
        if read_type == 'TP_read':
            cell = read.get_tag('BC')
            if cell not in cell_set:
                cell = 'UNASSIGNED'
            if cell not in quality_dict[sample][read_type]:
                quality_dict[sample][read_type][cell] = Counter()
            quality_dict[sample][read_type][cell].update({quality_group:1})
        else:
            quality_dict[sample][read_type].update({quality_group:1})
    bam.close()
    return quality_dict


def main():
    parser = argparse.ArgumentParser(description='', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--input',metavar='input', type=str, help='Input .bam file')
    parser.add_argument('-c','--cells',metavar='cells', type=str, help='Cell barcodes .txt file' )
    parser.add_argument('-s','--samplesheet', help='Samplesheet file')
    parser.add_argument('-o','--out_dir', help='Output directory')
    parser.add_argument('-p', '--prefix')

    args = parser.parse_args()

    cell_set = set([line.rstrip() for line in open(args.cells)])
    samplesheet_df = pd.read_csv(args.samplesheet, index_col=0)

    sample_set = set(samplesheet_df['sample_barcode'])

    mapping_group_dict = get_quality_dict(args.input, cell_set, sample_set)

    sample_dfs = {sample: pd.DataFrame(sample_dict['TP_read']).T for sample, sample_dict in mapping_group_dict.items()}

    for sample, sample_df in sample_dfs.items():
        sample_df = sample_df.div(sample_df.sum(axis=1), axis=0)
        sample_df['sample_barcode'] = sample
        sample_dfs[sample] = sample_df
    
    df_cell_level = pd.concat(sample_dfs,axis=0)
    df_cell_level.index = [s1+s2 for (s1,s2) in df_cell_level.index]

    df_cell_level.to_csv('{}/{}_quality_per_sample.csv'.format(args.out_dir, args.prefix))

    
    
    read_type_per_sample_df = pd.DataFrame({sample: get_read_type_stats(sample_dict) for sample, sample_dict in mapping_group_dict.items()}).T

    read_type_per_sample_df.to_csv('{}/{}_quality_per_sample.csv'.format(args.out_dir, args.prefix))

if __name__ == '__main__':
    main()
