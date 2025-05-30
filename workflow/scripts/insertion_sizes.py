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
def get_template_length_dict(bamfile, cell_set, sample_set):
    mapping_group_dict = {}
    read_length_dict = {}
    double_read_dict = {}

    bam = pysam.AlignmentFile(bamfile)

    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.mate_is_unmapped:
            continue
        if read.has_tag('SB'):
            sample = read.get_tag('SB')
        else:
            sample = read.get_tag('SM')
        if sample not in sample_set:
            sample = 'UNASSIGNED'
            
        if sample not in mapping_group_dict:
            mapping_group_dict[sample] = {}
            read_length_dict[sample] = {}
            double_read_dict[sample] = Counter()
        
        mapping_group = get_mapping_group(read)
        if mapping_group != 'Exon' and read.is_proper_pair and not read.is_unmapped:
            continue
        read_type = read.get_tag('XX')
        
        if read_type not in mapping_group_dict[sample]:
            mapping_group_dict[sample][read_type] = {}
            read_length_dict[sample][read_type] = {}
            
        
        
        read_start = int(read.get_tag('RS')) > 0
        if read_type == 'TP_read':
            if read.is_read2:
                if 'TCTTCTCTCCTCCTCC' in read.query_alignment_sequence:
                    double_read_dict[sample].update({'Double': 1})
                else:
                    double_read_dict[sample].update({'Single': 1})
                    
            cell = read.get_tag('BC')
            if cell not in cell_set:
                cell = 'UNASSIGNED'
            if cell not in mapping_group_dict[sample][read_type]:
                mapping_group_dict[sample][read_type][cell] = {}
                read_length_dict[sample][read_type][cell] = {}
            if read_start not in mapping_group_dict[sample][read_type][cell]:
                mapping_group_dict[sample][read_type][cell][read_start] = Counter()
                read_length_dict[sample][read_type][cell][read_start] = Counter()
            mapping_group_dict[sample][read_type][cell][read_start].update({read.template_length:1})
            if read.template_length == 0:
                read_length_dict[sample][read_type][cell][read_start].update({read.infer_query_length():1})
        else:
            if read_start not in mapping_group_dict[sample][read_type]:
                mapping_group_dict[sample][read_type][read_start] = Counter()
                read_length_dict[sample][read_type][read_start] = Counter()
            mapping_group_dict[sample][read_type][read_start].update({read.template_length:1})
            if read.template_length == 0:
                read_length_dict[sample][read_type][read_start].update({read.infer_query_length():1})
    bam.close()
    return mapping_group_dict, read_length_dict, double_read_dict


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

    inp = args.input

    template_length_dict, read_length_dict, double_read_dict = get_template_length_dict(inp, cell_set, sample_set)

    lol = []
    for sample_barcode, mapping_group_dict in template_length_dict.items():
        for mapping_group, read_start_group_dict in mapping_group_dict.items():
            if mapping_group == 'TP_read':
                for cell_barcode, actual_read_start_group_dict in read_start_group_dict.items():
                    for read_start, counter in actual_read_start_group_dict.items():
                        for k,v in counter.items():
                            if k < 0:
                                continue
                            lol.append([k,v, sample_barcode, mapping_group, cell_barcode, read_start])
            else:
                cell_barcode = 'NA'
                for read_start, counter in read_start_group_dict.items():
                    for k,v in counter.items():
                        if k < 0:
                            continue
                        lol.append([k,v, sample_barcode, mapping_group, cell_barcode, read_start])
    df = pd.DataFrame(lol, columns=['insert_size', 'count', 'sample_barcode', 'mapping_group', 'cell_barcode', 'read_start'])

    df.to_csv('{}/{}_insert_sizes_per_sample_barcode.csv'.format(args.out_dir, args.prefix))

    double_read_df = pd.DataFrame(double_read_dict).fillna(0).astype(int).T

    double_read_df.to_csv('{}/{}_double_reads_per_sample_barcode.csv'.format(args.out_dir, args.prefix))

if __name__ == '__main__':
    main()
