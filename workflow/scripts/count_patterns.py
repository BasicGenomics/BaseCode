import argparse
import pysam
from collections import Counter
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix
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
    patterns_per_cell = {}
    bam = pysam.AlignmentFile(bamfile)
    for read in bam.fetch(g_dict['seqid'], g_dict['start'], g_dict['end']):
        gene_name = determine_gene_tag(read)
        if gene_name != g_dict['gene_id']:
            continue
        if not read.has_tag('RM'):
            continue
        if not read.has_tag('BC'):
            continue
        if not read.has_tag('ST'):
            continue
        sample = read.get_tag('SM')
        if sample not in patterns_per_cell:
            patterns_per_cell[sample] = {}
        cell = read.get_tag('BC')
        if cell not in patterns_per_cell[sample]:
            patterns_per_cell[sample][cell] = {}
        status = read.get_tag('ST')
        if status not in patterns_per_cell[sample][cell]:
            patterns_per_cell[sample][cell][status] = {}
        RM_tag = read.get_tag('RM')
        if RM_tag not in patterns_per_cell[sample][cell][status]:
            patterns_per_cell[sample][cell][status][RM_tag] = Counter()
        XX_tag = read.get_tag('XX')
        patterns_per_cell[sample][cell][status][RM_tag].update({XX_tag: 1})
    patterns_per_cell_return = {}
    for sample, sample_dict in patterns_per_cell.items():
        patterns_per_cell_return[sample] = {}
        for cell, cell_dict in sample_dict.items():
            patterns_per_cell_return[sample][cell] = {}
            for status, status_dict in cell_dict.items():
                patterns_per_cell_return[sample][cell][status] = len(status_dict.keys())
    return patterns_per_cell_return, patterns_per_cell

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Count detected patterns and save to files', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--input',metavar='input', type=str, help='Input .bam file')
    parser.add_argument('-g','--gtf', metavar='gtf', type=str, help='gtf file with gene information')
    parser.add_argument('-m','--matrix-file', help='Matrix Output file')
    parser.add_argument('-t', '--threads', metavar='threads', type=int, default=1, help='Number of threads')

    args = parser.parse_args()

    gtffile = args.gtf
    bamfile = args.input
    matrix_file = args.matrix_file
    threads = args.threads

    gene_dicts = filterGeneDict(parse_gtf_gene_name(gtffile, None), bamfile)
    res = Parallel(n_jobs=threads, verbose = 3, backend='multiprocessing')(delayed(func)(bamfile,g_dict) for g,g_dict in gene_dicts.items())

    pattern_counting_dict = {}
    other_counting_dict = {}
    for gene, res_tuple in zip(gene_dicts.keys(), res):
        counting_res, detailed_res = res_tuple
        for sample, sample_dict in counting_res.items():
            for cell, cell_dict in sample_dict.items():
                if cell == '' or cell == 'Unassigned':
                    continue
                cellbarcode = sample+cell
                if cellbarcode not in pattern_counting_dict:
                    pattern_counting_dict[cellbarcode] = {}
                    other_counting_dict[cellbarcode] = {}
                pattern_counting_dict[cellbarcode][gene] = 0
                other_counting_dict[cellbarcode][gene] = 0
                for status, status_n in cell_dict.items():
                    if status in ['Molecule', 'ReadPair','ReadSingleton']:
                        pattern_counting_dict[cellbarcode][gene] += status_n
                    else:
                        other_counting_dict[cellbarcode][gene] += status_n
    pattern_counting_df = pd.DataFrame(pattern_counting_dict).fillna(0).astype(int)
    adata = sc.AnnData(pattern_counting_df.T)
    adata.X = csr_matrix(adata.X)
    adata.write(matrix_file)


