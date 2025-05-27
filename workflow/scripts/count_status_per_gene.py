import argparse
import pysam
from collections import Counter
import numpy as np
import pandas as pd
from joblib import delayed,Parallel

def determine_gene_tag(read):
    if read.has_tag('GE'):
        gene_exon = read.get_tag('GE')
    else:
        gene_exon = 'Unassigned'
    if read.has_tag('GI'):
        gene_intron = read.get_tag('GI')
    else:
        gene_intron = 'Unassigned'
    # if it maps to the intron or exon of a gene
    if gene_intron != 'Unassigned' or gene_exon != 'Unassigned':
        # if it is a junction read
        if gene_intron == gene_exon:
            gene = gene_intron
            # if it's an only intronic read
        elif gene_intron != 'Unassigned' and gene_exon == 'Unassigned':
            gene = gene_intron
            # if it's an only exonic read
        elif gene_exon != 'Unassigned' and gene_intron == 'Unassigned':
            gene = gene_exon
            # if the exon and intron gene tag contradict each other
        else:
            gene = ''
    else:
        gene = ''
    return gene
def parse_gtf_gene_name(gtffile, contig, ban_set=set()):
    gene_list = []
    with open(gtffile, 'r') as f:
        for line in f:
            l = line.split('\t')
            if len(l) < 8:
                continue
            if l[2] == 'gene':
                if l[2] in ban_set:
                    continue
                if contig is not None:
                    if l[0] == contig:
                        try:
                            gene_list.append({'gene_id': l[8].split(' ')[5].replace('"', '').strip(';'), 'seqid':l[0], 'start':int(l[3]), 'end':int(l[4]), 'strand': l[6]})
                        except:
                            gene_list.append({'gene_id': l[8].split(' ')[1].replace('"', '').strip(';'), 'seqid':l[0], 'start':int(l[3]), 'end':int(l[4]), 'strand': l[6]})
                    else:
                        continue
                else:
                    try:
                        gene_list.append({'gene_id': l[8].split(' ')[5].replace('"', '').strip(';'), 'seqid':l[0], 'start':int(l[3]), 'end':int(l[4]), 'strand': l[6]})
                    except:
                        gene_list.append({'gene_id': l[8].split(' ')[1].replace('"', '').strip(';'), 'seqid':l[0], 'start':int(l[3]), 'end':int(l[4]), 'strand': l[6]})
    gene_dict = {g['gene_id']: g for g in gene_list}
    return gene_dict
def filterGeneDict(gene_dict, bam_infile):
        bam = pysam.AlignmentFile(bam_infile,'rb')
        contigs = {d['SN'] for d in bam.header['SQ']}
        new_gene_dict = {}
        i = 0
        for gene, g_dict in gene_dict.items():
            if g_dict['seqid'] not in contigs:
                i += 1
                continue
            new_gene_dict[gene] = g_dict
        bam.close()
        print("Filtered {} genes which are not found in bam file".format(i))
        return new_gene_dict
def func(bamfile, g_dict):
    bam = pysam.AlignmentFile(bamfile)
    status_counter = {}
    for read in bam.fetch(g_dict['seqid'], g_dict['start'], g_dict['end']):
        gene_name = determine_gene_tag(read)
        if gene_name != g_dict['gene_id']:
            continue
        
        sample = read.get_tag('SM')
        if sample not in status_counter:
            status_counter[sample] = {}
        
        if sample == 'Unassigned':
            sample_barcode = 'Unassigned'
        else:
            sample_barcode = read.get_tag('SB')
        
        if sample_barcode not in status_counter[sample]:
            status_counter[sample][sample_barcode] = {}
        
        if read.has_tag('ST'):
            status = read.get_tag('ST')
        else:
            status = 'NoTag'
            
        XX_tag = read.get_tag('XX')
        if XX_tag not in status_counter[sample][sample_barcode]:
            status_counter[sample][sample_barcode][XX_tag] = Counter()
        status_counter[sample][sample_barcode][XX_tag].update({status:1})
    return status_counter

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Count reconstruction status for each read per gene', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--input',metavar='input', type=str, help='Input .bam file')
    parser.add_argument('-g','--gtf', metavar='gtf', type=str, help='gtf file with gene information')
    parser.add_argument('--counts', help='Output counts per gene file')
    parser.add_argument('--sum', help='Output sum per sample barcode file')
    parser.add_argument('--fraction', help='Output status fraction per gene file')
    parser.add_argument('-t', '--threads', metavar='threads', type=int, default=1, help='Number of threads')

    args = parser.parse_args()

    gtffile = args.gtf
    bamfile = args.input
    threads = args.threads

    counts_out = args.counts
    sum_out = args.sum
    fraction_out = args.fraction

    gene_dicts = filterGeneDict(parse_gtf_gene_name(gtffile, None), bamfile)
    res = Parallel(n_jobs=threads, verbose = 3, backend='multiprocessing')(delayed(func)(bamfile, g_dict) for g,g_dict in gene_dicts.items())

    total_counts_per_gene = {}
    total_sum_over_genes = {}
    fraction_per_gene = {}
    for gene, res_dict in zip(gene_dicts.keys(), res):
        for sample, status_dicts in res_dict.items():
            if sample not in total_sum_over_genes:
                total_sum_over_genes[sample] = Counter()
                total_counts_per_gene[sample] = {}
                fraction_per_gene[sample] = {}
            
            for sample_barcode, sample_barcode_dict in status_dicts.items():
                if sample_barcode not in total_sum_over_genes[sample]:
                    total_sum_over_genes[sample][sample_barcode] = {}
                    total_counts_per_gene[sample][sample_barcode] = {}
                    fraction_per_gene[sample][sample_barcode] = {}
                total_counts_per_gene[sample][sample_barcode][gene] = 0
                for read_type, status_dict in sample_barcode_dict.items():
                    if read_type not in fraction_per_gene[sample][sample_barcode]:
                        fraction_per_gene[sample][sample_barcode][read_type] = {}
                        total_sum_over_genes[sample][sample_barcode][read_type] = Counter()
                    status_sum = sum(status_dict.values())
                    total_counts_per_gene[sample][sample_barcode][gene] += status_sum

                    total_sum_over_genes[sample][sample_barcode][read_type].update(status_dict)

                    status_fraction = {status: n/status_sum for status, n in status_dict.items()}
                    fraction_per_gene[sample][sample_barcode][read_type][gene] = status_fraction
    
    fraction_list = []
    for sample, sample_dict in fraction_per_gene.items():
        for sample_barcode, sample_barcode_dict in sample_dict.items():
            for read_type, read_type_dict in sample_barcode_dict.items():
                for gene, gene_dict in read_type_dict.items():
                    for status, fraction in gene_dict.items():
                        fraction_list.append([sample, sample_barcode, read_type, gene, status, fraction])

    fraction_df = pd.DataFrame(fraction_list, columns=['sample', 'sample_barcode', 'read_type', 'gene', 'status', 'fraction'])

    counts_per_gene_list = []
    for sample, sample_dict in total_counts_per_gene.items():
        for sample_barcode, sample_barcode_dict in sample_dict.items():
            for gene, count in sample_barcode_dict.items():
                counts_per_gene_list.append([sample, sample_barcode, gene, count])

    counts_per_gene_df = pd.DataFrame(counts_per_gene_list, columns = ['sample', 'sample_barcode', 'gene', 'count'])

    sum_over_genes_list = []
    for sample, sample_dict in total_sum_over_genes.items():
        for sample_barcode, sample_barcode_dict in sample_dict.items():
            for read_type, read_type_dict in sample_barcode_dict.items():
                for status, count in read_type_dict.items():
                    sum_over_genes_list.append([sample, sample_barcode, read_type, status, count])

    sum_over_genes_df = pd.DataFrame(sum_over_genes_list, columns = ['sample', 'sample_barcode', 'read_type', 'status', 'count'])

    fraction_df.to_csv(fraction_out)
    counts_per_gene_df.to_csv(counts_out)
    sum_over_genes_df.to_csv(sum_out)
    
