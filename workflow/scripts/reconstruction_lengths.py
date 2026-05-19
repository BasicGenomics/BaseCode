import pysam
import argparse
import pandas as pd
from joblib import delayed,Parallel
from multiprocessing import Process, Manager
import warnings

# CIGAR op codes (pysam): 0=M, 1=I, 2=D, 3=N, 4=S, 5=H, 6=P, 7==, 8=X
CIGAR_M_LIKE = {0, 7, 8}
CIGAR_D = 2
CIGAR_N = 3
CIGAR_S = 4


def process_cigar(read):

    aligned_bases = 0  # M/=/X — read bases on reference
    deletion_bases = 0  # D — reference bases not covered, non-intronic
    intron_bases = 0  # N — skipped (intron)
    soft_clip = 0
    del_lengths = []
    aligned_blocks = []  # absolute genomic (start, end) per M block; D gaps are implicit
    ref_pos = read.reference_start

    for op, length in read.cigartuples:
        if op in CIGAR_M_LIKE:
            aligned_blocks.append((ref_pos, ref_pos + length))
            aligned_bases += length
            ref_pos += length
        elif op == CIGAR_D:
            deletion_bases += length
            del_lengths.append(length)
            ref_pos += length
        elif op == CIGAR_N:
            intron_bases += length
            ref_pos += length
        elif op == CIGAR_S:
            soft_clip += length
            # soft clips do not advance reference position

    read_length = read.infer_read_length()  # includes soft clips; None if unavailable
    if read_length is None:
        read_length = read.query_length
    ref_span = aligned_bases + deletion_bases + intron_bases
    denom = aligned_bases + deletion_bases
    missing_fraction = deletion_bases / denom if denom > 0 else 0.0

    return [
        read_length,
        aligned_bases,
        ';'.join('{}:{}'.format(s, e) for s, e in aligned_blocks),
        deletion_bases,
        len(del_lengths),
        ';'.join(map(str, del_lengths)),
        intron_bases,
        soft_clip,
        ref_span,
        missing_fraction,
    ]

def process_read(bamfile, seqid):
    res = []
    with pysam.AlignmentFile(bamfile, 'rb') as bam:
        for mol in bam.fetch(seqid):
            qname = mol.query_name
            sample_name = '_'.join(qname.split(':')[0].split('_')[:-1])
            cigar_stats = process_cigar(mol)
            res_i = [mol.query_name,
                        sample_name, mol.get_tag('BC'),
                        mol.get_tag('NR'), mol.get_tag('ER'), mol.get_tag('IR'),
                        mol.get_tag('TC'), mol.get_tag('IC'), mol.get_tag('FC'),
                        mol.get_tag('XT'),
                        mol.reference_name, mol.reference_start, mol.reference_end,
                        mol.get_tag('CC'), mol.get_tag('SC'),
                        mol.query_length, cigar_stats[4],
                        mol.get_tag('F1'), mol.get_tag('T1')]
            res_i.extend(cigar_stats[:4])   # skip index 4 (len(del_lengths) == GAPS, already added)
            res_i.extend(cigar_stats[5:])
            res_i.append(mol.get_tag('CV'))
            res.append(res_i)
    return res


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Reconstruction length statistics', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--input',metavar='input', type=str, help='Input .bam file')
    parser.add_argument('-o','--outfile', help='Output file')
    parser.add_argument('-t', '--threads', metavar='threads', type=int, default=1, help='Number of threads')
    parser.add_argument('-m','--mode',default='basic')


    args = parser.parse_args()


    bamfile = args.input
    outfile = args.outfile
    threads = args.threads

    with pysam.AlignmentFile(bamfile, 'rb') as bam:
        refs = bam.references

    res_contig = Parallel(n_jobs=threads, verbose=0, backend='multiprocessing')(delayed(process_read)(bamfile, ref) for ref in refs)
    res = [row for res_contig in res_contig for row in res_contig]

    df = pd.DataFrame(res, columns = ['MOL_NAME',
                                    'SM','BC',
                                    'NR', 'ER', 'IR',
                                    'TC', 'IC', 'FC',
                                    'XT','SEQID','ref_start','ref_end',
                                    'CC','SC',
                                    'QL', 'GAPS', 'F1', 'T1',
                                    'read_length', 'aligned_bases', 'aligned_blocks',
                                    'deletion_bases', 'del_lengths',
                                    'intron_bases', 'soft_clip', 'ref_span', 'missing_fraction','read_depth_profile'])
    
    df['ANY_TC'] = df['TC'] > 0
    df['ANY_IC'] = df['IC'] > 0
    df['ANY_FC'] = df['FC'] > 0
    df['ANY_GAPS'] = df['GAPS'] > 0

    
    df_basic = df[['MOL_NAME','SM','BC', 'TC', 'IC', 'FC', 'XT','CC','SC','QL', 'GAPS', 'F1', 'T1','ANY_TC','ANY_IC','ANY_FC','ANY_GAPS']]
   
    df_basic.to_csv(outfile)

    if args.mode =='comprehensive':
        stem = outfile.rstrip('.csv')
        df.to_csv(f'{stem}_detailed.csv')
    
        
