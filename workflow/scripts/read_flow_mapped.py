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
        return ES_tag
    if 'Assigned' in IS_tag:
        return IS_tag
    if 'Unassigned' in ES_tag and 'Unassigned' in IS_tag:
        return ES_tag
    return ES_tag

def keep_read(read):
    if read.is_mapped and read.mate_is_mapped:
        if read.is_read1:
            return True
        else:
            return False
    if read.is_mapped and read.mate_is_unmapped:
        return True
    if read.is_unmapped and read.mate_is_unmapped:
        if read.is_read1:
            return True
        else:
            return False
    if read.is_unmapped and read.mate_is_mapped:
        return False

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

def get_mapping_group_dict(bamfile, cell_set, sample_set):
    mapping_group_dict = {}
    bam = pysam.AlignmentFile(bamfile)

    for read in bam.fetch(until_eof=True):
        if not keep_read(read):
            continue
        sample = read.get_tag('SB')
        if sample not in sample_set:
            sample = 'UNASSIGNED'
        
        if sample not in mapping_group_dict:
            mapping_group_dict[sample] = {}
        
        mapping_group = get_mapping_group(read)
        quality_group = read.mapping_quality
        read_type = read.get_tag('XX')
        
        if read_type not in mapping_group_dict[sample]:
            if read_type == 'TP_read':
                mapping_group_dict[sample][read_type] = {}
            else:
                mapping_group_dict[sample][read_type] = Counter()
        
        if read_type == 'TP_read':
            cell = read.get_tag('BC')
            if cell not in cell_set:
                cell = 'UNASSIGNED'
            if cell not in mapping_group_dict[sample][read_type]:
                mapping_group_dict[sample][read_type][cell] = Counter()
            mapping_group_dict[sample][read_type][cell].update({(mapping_group, quality_group):1})
        else:
            mapping_group_dict[sample][read_type].update({(mapping_group, quality_group):1})

    bam.close()
    return mapping_group_dict


def main():
    parser = argparse.ArgumentParser(description='', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--input',metavar='input', type=str, help='Input .bam file')
    parser.add_argument('-c','--cells',metavar='cells', type=str, help='Cell barcodes .txt file' )
    parser.add_argument('-s','--sample-barcodes', help='Sample barcode .txt file')
    parser.add_argument('--mapping-group-out', help='Output mapping group')

    args = parser.parse_args()

    mapping_group_out = args.mapping_group_out
    reassignment_out = args.reassignment_out

    cell_set = set([line.rstrip() for line in open(args.cells)])
    sample_set = set([line.rstrip() for line in open(args.sample_barcodes)])

    mapping_group_dict = get_mapping_group_dict(args.input, cell_set, sample_set)

    mapping_group_list = []
    for sample, sample_dict in mapping_group_dict.items():
        for read_type, read_type_dict in sample_dict.items():
            if read_type == 'TP_read':
                for cell, cell_dict in read_type_dict.items():
                    for (mapping_group, quality_group), count in cell_dict.items():
                        mapping_group_list.append([sample, read_type, cell, mapping_group, quality_group, count])
            else:
                for (mapping_group, quality_group), count in read_type_dict.items():
                    mapping_group_list.append([sample, read_type, 'NA', mapping_group, quality_group, count])
    mapping_group_df = pd.DataFrame(mapping_group_list, columns= ['sample', 'read_type', 'cell', 'mapping_group','quality_group','count'])


    mapping_group_df.to_csv(mapping_group_out)


if __name__ == '__main__':
    main()
