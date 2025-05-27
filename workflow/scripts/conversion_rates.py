import pysam
from collections import Counter
import pandas as pd
import sys
import os
import numpy as np
import argparse
import time
from collections import Counter
from joblib import delayed,Parallel
from multiprocessing import Process, Manager
from pyfaidx import Fasta
import itertools


def base_conversion_complement(base_conversions):
    complement_dict_lower = {'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
    complement_dict_upper = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complement = [''.join([complement_dict_lower[bc[0]], complement_dict_upper[bc[1]]]) for bc in base_conversions]
    return complement

def find_eligible_positions(fasta_file, base_conversions, g_dict):
    # Find position in the gene region which contain the reference bases which may be converted

    fasta_ref = Fasta(fasta_file)
    whole_sequence = fasta_ref[g_dict['seqid']][g_dict['start']:g_dict['end']].seq.upper()

    reference_bases = set([bc[0].upper() for bc in base_conversions])

    eligible_position_set = set([pos+g_dict['start']+1 for pos, base in enumerate(whole_sequence) if base in reference_bases])

    return whole_sequence, eligible_position_set

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

def get_insertions_locs(cigtuples):
    insertion_locs = []
    l = 0
    for c in cigtuples:
        if c[0] == 0:
            l += c[1]
        elif c[0] == 1:
            for i in range(c[1]):
                insertion_locs.append(l)
                l += 1
    return insertion_locs

def compare_to_reference(read_seq, read_qual, ref_seq):
    total_content = Counter({'A': 0, 'T': 0, 'C': 0, 'G': 0})
    specific_content = Counter()
    for i, (read_base, base_qual, ref_base) in enumerate(zip(read_seq, read_qual, ref_seq)):
        if base_qual < 15:
            continue
        if read_base == ref_base:
            total_content.update({ref_base: 1})
        else:
            k = ''.join((ref_base.lower(),read_base))
            specific_content.update({k: 1})
    return specific_content, total_content

def find_base_conversions(read, whole_sequence,base_conversions, start_offset):

    intervals = read.get_blocks()
    offset = 0
    read_seq = read.query_alignment_sequence
    read_qual = read.query_alignment_qualities

    cigtuples = read.cigartuples
    insertion_locs = get_insertions_locs(cigtuples)

    if len(insertion_locs) > 0:
        read_seq = "".join([char for idx, char in enumerate(read_seq) if idx not in insertion_locs])
        read_qual = [qual for idx, qual in enumerate(read_qual) if idx not in insertion_locs]
    total_content = Counter()
    specific_content = Counter()
    for interval in intervals:
        interval_start = interval[0]
        interval_end = interval[1]-1
        interval_length = interval_end-interval_start
        start_idx = interval_start-start_offset
        end_idx = interval_end-start_offset
        if start_idx < 0:
            interval_start -= start_idx
            offset -= start_idx
            interval_length += start_idx
            start_idx = 0
        if end_idx > len(whole_sequence):
            end_idx = len(whole_sequence)
            interval_length = end_idx - offset
        ref_seq = whole_sequence[start_idx:end_idx]
        read_seq_interval = read_seq[offset:offset + interval_length]
        read_qual_interval = read_qual[offset:offset + interval_length]
        offset += interval_length + 1
        specific_content_interval, total_content_interval = compare_to_reference(read_seq_interval, read_qual_interval, ref_seq)
        specific_content.update(specific_content_interval)
        total_content.update(total_content_interval)
    return specific_content, total_content

def get_conversion_stats_per_cell(bam_infile,fasta_file, base_conversions, g_dict, cell_set, only_samples):
    whole_sequence, eligible_position_set = find_eligible_positions(fasta_file, base_conversions, g_dict)
    bam = pysam.AlignmentFile(bam_infile, 'r')
    
    sc_per_cell = {}
    tc_per_cell = {}
    for read in bam.fetch(g_dict['seqid'], g_dict['start'], g_dict['end']):
        gene_name = determine_gene_tag(read)
        if gene_name != g_dict['gene_id'] or read.mapping_quality < 1 or not read.has_tag('BC'):
            continue
        cell = read.get_tag('SM')
        if not only_samples:
            cell += read.get_tag('BC')
            if cell not in cell_set:
                continue
            if cell not in sc_per_cell:
                sc_per_cell[cell] = Counter()
                tc_per_cell[cell] = Counter()
        else:
            if cell not in sc_per_cell:
                sc_per_cell[cell] = Counter()
                tc_per_cell[cell] = Counter()
        specific_conversions, total_content = find_base_conversions(read, whole_sequence, base_conversions, g_dict['start'])
        sc_per_cell[cell].update(specific_conversions)
        tc_per_cell[cell].update(total_content)
    return sc_per_cell, tc_per_cell
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
def main():
    parser = argparse.ArgumentParser(description='Calculate conversion rates', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--input',metavar='input', type=str, help='Input .bam file')
    parser.add_argument('-f','--fasta', metavar='fasta', type=str, help='Reference fasta file')
    parser.add_argument('-g','--gtf', metavar='gtf', type=str, help='gtf file with gene information')
    parser.add_argument('-p','--prefix', help='Run prefix')
    parser.add_argument('-s','--samplesheet', help='Samplesheet file')
    parser.add_argument('-o','--out_dir', help='Output directory')
    parser.add_argument('-t', '--threads', metavar='threads', type=int, default=1, help='Number of threads')
    parser.add_argument('-bc', '--base-conversions', nargs='+', type=str, default=['aG','gA'], help='Base conversions used to reconstruct')
    parser.add_argument('--only-samples', action='store_true')

    args = parser.parse_args()

    bam_infile = args.input
    fastafile = args.fasta
    gtffile = args.gtf
    prefix = args.prefix
    samplesheet_file = args.samplesheet
    threads = args.threads
    out_dir = args.out_dir
    isExist = os.path.exists(out_dir)
    only_samples = args.only_samples
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(out_dir)

   
    samplesheet_df = pd.read_csv(samplesheet_file, index_col=0)

    gene_dict = filterGeneDict(parse_gtf_gene_name(gtffile, None), bam_infile)
    
    base_conversions = args.base_conversions

    cell_set = set(samplesheet_df.index)


    res = Parallel(n_jobs=threads, verbose = 0, backend='multiprocessing')(delayed(get_conversion_stats_per_cell)(bam_infile,fastafile,base_conversions, g_dict, cell_set, only_samples) for g,g_dict in gene_dict.items())
    tc_per_cell_stranded = {s: {} for s in ['+', '-']}
    sc_per_cell_stranded = {s: {} for s in ['+', '-']}

    sample_list = []
    for g, (sc_per_cell, tc_per_cell) in zip(gene_dict.keys(), res):
        g_dict = gene_dict[g]
        for cell, sc_counter in sc_per_cell.items():
            tc_counter = tc_per_cell[cell]
            if cell not in sc_per_cell_stranded[g_dict['strand']]:
                sample_list.append(cell)
                sc_per_cell_stranded[g_dict['strand']][cell] = Counter()
                tc_per_cell_stranded[g_dict['strand']][cell] = Counter()
            sc_per_cell_stranded[g_dict['strand']][cell].update(sc_counter)
            tc_per_cell_stranded[g_dict['strand']][cell].update(tc_counter)
    sample_list = list(set(sample_list))
    
    if only_samples:
        cell_set = sample_list
    sc_per_cell_total = {cell: Counter() for cell in cell_set}
    tc_per_cell_total = {cell: Counter() for cell in cell_set}

    for strand, cell_dict in sc_per_cell_stranded.items():
        for cell, cell_counter in cell_dict.items():
            if cell in sc_per_cell_total:
                sc_per_cell_total[cell].update(cell_counter)
    for strand, cell_dict in tc_per_cell_stranded.items():
        for cell, cell_counter in cell_dict.items():
            if cell in tc_per_cell_total:
                tc_per_cell_total[cell].update(cell_counter)
    
    sc_key_to_total_key = {}
    sc_list = []
    for cell, cell_dict in sc_per_cell_total.items():
        for sc_key in cell_dict.keys():
            sc_list.append(sc_key)
    sc_list = list(set(sc_list))
    for sc_key in sc_list:
        sc_key_to_total_key[sc_key] = sc_key[0].upper()
    
    conv_rate_per_cell_total = {}
    for cell, sc_counter in sc_per_cell_total.items():
        conv_rate_per_cell_total[cell] = {}
        tc_dict = tc_per_cell_total[cell]
        for sc_key, sc_n in sc_counter.items():
            if sc_key not in sc_key_to_total_key:
                continue
            tc_key = sc_key_to_total_key[sc_key]
            if tc_key == 'N':
                continue
            if tc_dict[tc_key] == 0:
                conv_rate_per_cell_total[cell][sc_key] = 0
            else:
                conv_rate_per_cell_total[cell][sc_key] = sc_n/tc_dict[tc_key]
   
    conv_rate_per_cell_stranded = {s: {cell: {} for cell in cell_set} for s in ['+', '-']}
    
    for strand, strand_sc_per_cell in sc_per_cell_stranded.items():
        for cell, sc_counter in strand_sc_per_cell.items():
            if cell not in cell_set:
                continue
            tc_dict = tc_per_cell_stranded[strand][cell]
            for sc_key, sc_n in sc_counter.items():
                if sc_key not in sc_key_to_total_key:
                    continue
                tc_key = sc_key_to_total_key[sc_key]
                if tc_key == 'N':
                    continue
                if tc_dict[tc_key] == 0:
                    conv_rate_per_cell_total[cell][sc_key] = 0
                else:
                    conv_rate_per_cell_stranded[strand][cell][sc_key] = sc_n/tc_dict[tc_key]
    

    df_conversion_rate_total = pd.DataFrame.from_dict(conv_rate_per_cell_total, orient='index')
    df_conversion_rate_total.to_csv('{}/{}_conversion_rate_total.csv'.format(out_dir, prefix))
    df_conversion_rate_pos = pd.DataFrame.from_dict(conv_rate_per_cell_stranded['+'], orient='index')
    df_conversion_rate_pos.to_csv('{}/{}_conversion_rate_pos.csv'.format(out_dir, prefix))
    df_conversion_rate_neg = pd.DataFrame.from_dict(conv_rate_per_cell_stranded['-'], orient='index')
    df_conversion_rate_neg.to_csv('{}/{}_conversion_rate_neg.csv'.format(out_dir, prefix))

if __name__ == "__main__":
    main()
