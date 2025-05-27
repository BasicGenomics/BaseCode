import pandas as pd
import pysam
from collections import Counter
import argparse

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

def get_status_group_dict(bamfile, cell_set, sample_set):
    status_group_dict = {}
    bam = pysam.AlignmentFile(bamfile)

    for read in bam.fetch(until_eof=True):
        if not keep_read(read):
            continue
        sample = read.get_tag('SB')
        if sample not in sample_set:
            sample = 'UNASSIGNED'
        
        if sample not in status_group_dict:
            status_group_dict[sample] = {}
        
        status_group = read.get_tag('ST')
        quality_group = read.mapping_quality
        read_type = read.get_tag('XX')
        
        if read_type not in status_group_dict[sample]:
            if read_type == 'TP_read':
                status_group_dict[sample][read_type] = {}
            else:
                status_group_dict[sample][read_type] = Counter()
        
        if read_type == 'TP_read':
            cell = read.get_tag('BC')
            if cell not in cell_set:
                cell = 'UNASSIGNED'
            if cell not in status_group_dict[sample][read_type]:
                status_group_dict[sample][read_type][cell] = Counter()
            status_group_dict[sample][read_type][cell].update({(status_group, quality_group):1})
        else:
            status_group_dict[sample][read_type].update({(status_group, quality_group):1})
    bam.close()
    return status_group_dict

def main():
    parser = argparse.ArgumentParser(description='', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--input',metavar='input', type=str, help='Input .bam file')
    parser.add_argument('-c','--cells',metavar='cells', type=str, help='Cell barcodes .txt file' )
    parser.add_argument('-s','--sample-barcodes', help='Sample barcode .txt file')
    parser.add_argument('--status-group-out', help='Output status group')

    args = parser.parse_args()

    status_group_out = args.status_group_out

    cell_set = set([line.rstrip() for line in open(args.cells)])
    sample_set = set([line.rstrip() for line in open(args.sample_barcodes)])

    status_group_dict = get_status_group_dict(args.input, cell_set, sample_set)

    status_group_list = []
    for sample, sample_dict in status_group_dict.items():
        for read_type, read_type_dict in sample_dict.items():
            if read_type == 'TP_read':
                for cell, cell_dict in read_type_dict.items():
                    for (status_group, quality_group), count in cell_dict.items():
                        status_group_list.append([sample, read_type, cell, status_group, quality_group, count])
            else:
                for (status_group, quality_group), count in read_type_dict.items():
                    status_group_list.append([sample, read_type, 'NA', status_group, quality_group, count])
    status_group_df = pd.DataFrame(status_group_list, columns= ['sample', 'read_type', 'cell', 'status_group','quality_group','count'])

    status_group_df.to_csv(status_group_out)

if __name__ == '__main__':
    main()
