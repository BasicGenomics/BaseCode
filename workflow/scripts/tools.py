import pysam

def parse_col_nine(col9):
    info_list = col9.split(';')
    info_dict = {t.split('=')[0]: t.split('=')[1] for t in info_list}
    return info_dict

def parse_gtf(gtffile, contig, ban_set=set()):
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
                        info_dict = parse_col_nine(l[8])
                        if 'gene_name' not in info_dict:
                            continue
                        gene_list.append({'gene_id': info_dict['gene_name'], 'seqid': l[0], 'start': int(l[3]), 'end': int(l[4]), 'strand': l[6]})
                else:
                    info_dict = parse_col_nine(l[8])
                    if 'gene_name' not in info_dict:
                        continue
                    gene_list.append({'gene_id': info_dict['gene_name'], 'seqid': l[0], 'start': int(l[3]), 'end': int(l[4]), 'strand': l[6]})
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
        return new_gene_dict

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