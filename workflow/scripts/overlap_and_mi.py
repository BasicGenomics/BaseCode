import argparse
import pysam
from collections import Counter
import numpy as np
import pandas as pd
from joblib import delayed,Parallel
from tools import parse_gtf, filterGeneDict, determine_gene_tag




def func(bamfile, g_dict):
    OV_counter = {}
    MI_counter = {}
    MI_OV_counter = {}
    bam = pysam.AlignmentFile(bamfile)
    for read in bam.fetch(g_dict['seqid'], g_dict['start'], g_dict['end']):
        gene_name = determine_gene_tag(read)
        if gene_name != g_dict['gene_id']:
            continue
        if not read.has_tag('OV') or not read.has_tag('MI') or not read.has_tag('ST') or not read.has_tag('BC'):
            continue
        sample = read.get_tag('SM')
        if sample not in OV_counter:
            OV_counter[sample] = {}
            MI_counter[sample] = {}
            MI_OV_counter[sample] = {}
        cell = read.get_tag('BC')
        if cell not in OV_counter[sample]:
            OV_counter[sample][cell] = {}
            MI_counter[sample][cell] = {}
            MI_OV_counter[sample][cell] = {}
        read_type = read.get_tag('XX')
        if read_type not in OV_counter[sample][cell]:
            OV_counter[sample][cell][read_type] = {}
            MI_counter[sample][cell][read_type] = {}
            MI_OV_counter[sample][cell][read_type] = {}
        status = read.get_tag('ST')
        if status not in OV_counter[sample][cell][read_type]:
            OV_counter[sample][cell][read_type][status] = Counter()
            MI_counter[sample][cell][read_type][status] = Counter()
            MI_OV_counter[sample][cell][read_type][status] = Counter()
        
        overlap = read.get_tag('OV')
        mutual_information = np.round(read.get_tag('MI'),2)
        
        OV_counter[sample][cell][read_type][status].update({overlap:1})
        MI_counter[sample][cell][read_type][status].update({mutual_information:1})
        MI_OV_counter[sample][cell][read_type][status].update({(overlap,mutual_information):1})
        
    return OV_counter, MI_counter, MI_OV_counter

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Collect overlap and mutual information (MI) for reconstruction', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--input',metavar='input', type=str, help='Input .bam file')
    parser.add_argument('-g','--gtf', metavar='gtf', type=str, help='gtf file with gene information')
    parser.add_argument('--genes-out', help='Output cell file')
    parser.add_argument('--cells-out', help='Output gene file')
    parser.add_argument('-t', '--threads', metavar='threads', type=int, default=1, help='Number of threads')

    args = parser.parse_args()

    gtffile = args.gtf
    bamfile = args.input 
    threads = args.threads

    gene_out = args.genes_out
    cell_out = args.cells_out

    gene_dicts = filterGeneDict(parse_gtf(gtffile, None), bamfile)

    res = Parallel(n_jobs=threads, verbose = 3, backend='multiprocessing')(delayed(func)(bamfile, g_dict) for g,g_dict in gene_dicts.items())

    OV_per_cell = {}
    OV_per_gene = {}

    MI_per_gene = {}
    MI_per_cell = {}

    MI_OV_per_cell = {}
    MI_OV_per_gene = {}

    for gene, (OV_counter, MI_counter, MI_OV_counter) in zip(gene_dicts.keys(), res): 
        for sample, sample_dict in OV_counter.items():
            
            MI_counter_sample = MI_counter[sample]
            MI_OV_counter_sample = MI_OV_counter[sample]
            
            if sample not in OV_per_cell:
                OV_per_cell[sample] = {}
                OV_per_gene[sample] = {}
                
                MI_per_gene[sample] = {}
                MI_per_cell[sample] = {}
                
                MI_OV_per_cell[sample] = {}
                MI_OV_per_gene[sample] = {}
            
            OV_per_gene[sample][gene] = {}
            MI_per_gene[sample][gene] = {}
            MI_OV_per_gene[sample][gene] = {}
            
            for cell, cell_dict in sample_dict.items():
                MI_counter_sample_cell = MI_counter_sample[cell]
                MI_OV_counter_sample_cell = MI_OV_counter_sample[cell]
                
                if cell not in OV_per_cell[sample]:
                    OV_per_cell[sample][cell] = {}
                    MI_per_cell[sample][cell] = {}
                    MI_OV_per_cell[sample][cell] = {}
                
                for read_type, read_type_dict in cell_dict.items():
                    MI_counter_sample_cell_read_type = MI_counter_sample_cell[read_type]
                    MI_OV_counter_sample_cell_read_type = MI_OV_counter_sample_cell[read_type]
                    
                    if read_type not in OV_per_cell[sample][cell]:
                        OV_per_cell[sample][cell][read_type] = {}
                        MI_per_cell[sample][cell][read_type] = {}
                        MI_OV_per_cell[sample][cell][read_type] = {}
                    
                    OV_per_gene[sample][gene][read_type] = {}
                    MI_per_gene[sample][gene][read_type] = {}
                    MI_OV_per_gene[sample][gene][read_type] = {}
                    
                    for status, status_dict in read_type_dict.items():
                        MI_counter_sample_cell_read_type_status = MI_counter_sample_cell_read_type[status]
                        MI_OV_counter_sample_cell_read_type_status = MI_OV_counter_sample_cell_read_type[status]
                        
                        if status not in OV_per_cell[sample][cell][read_type]:
                            OV_per_cell[sample][cell][read_type][status] = Counter()
                            MI_per_cell[sample][cell][read_type][status] = Counter()
                            MI_OV_per_cell[sample][cell][read_type][status] = Counter()
                        
                        OV_per_gene[sample][gene][read_type][status] = status_dict
                        MI_per_gene[sample][gene][read_type][status] = MI_counter_sample_cell_read_type_status
                        MI_OV_per_gene[sample][gene][read_type][status] = MI_OV_counter_sample_cell_read_type_status
                        
                        OV_per_cell[sample][cell][read_type][status].update(status_dict)
                        MI_per_cell[sample][cell][read_type][status].update(MI_counter_sample_cell_read_type_status)
                        MI_OV_per_cell[sample][cell][read_type][status].update(MI_OV_counter_sample_cell_read_type_status)
    MI_OV_per_cell_list = []
    for sample, sample_dict in MI_OV_per_cell.items():
        for cell, cell_dict in sample_dict.items():
            if cell == '':
                cell = 'NA'
            for read_type, read_type_dict in cell_dict.items():
                for status, mi_ov_counter in read_type_dict.items():
                    for (mi, ov), count in mi_ov_counter.items():
                        MI_OV_per_cell_list.append([sample, cell, read_type, status, mi, ov, count])
    MI_OV_per_gene_list = []
    for sample, sample_dict in MI_OV_per_gene.items():
        for gene, gene_dict in sample_dict.items():
            for read_type, read_type_dict in gene_dict.items():
                for status, mi_ov_counter in read_type_dict.items():
                    for (mi, ov), count in mi_ov_counter.items():
                        MI_OV_per_gene_list.append([sample, gene, read_type, status, mi, ov, count])
    
    mi_ov_per_cell_df = pd.DataFrame(MI_OV_per_cell_list, columns = ['sample', 'cell', 'read_type', 'status', 'overlap', 'adjusted_mutual_information', 'count'])
    mi_ov_per_cell_df.to_csv(cell_out)

    mi_ov_per_gene_df = pd.DataFrame(MI_OV_per_gene_list, columns = ['sample', 'gene', 'read_type', 'status', 'overlap', 'adjusted_mutual_information', 'count'])
    mi_ov_per_gene_df.to_csv(gene_out)

