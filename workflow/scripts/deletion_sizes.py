import argparse
import os

import numpy as np
import pandas as pd
import pysam

# CIGAR op codes (pysam): 0=M, 1=I, 2=D, 3=N, 4=S, 5=H, 6=P, 7==, 8=X
CIGAR_M_LIKE = {0, 7, 8}
CIGAR_D = 2
CIGAR_N = 3
CIGAR_S = 4


def _tag_or_na(read, tag):
    return read.get_tag(tag) if read.has_tag(tag) else 'NA'


def per_read_stats(bamfile):
    rows = []
    bam = pysam.AlignmentFile(bamfile)
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.cigartuples is None:
            continue

        aligned_bases = 0  # M/=/X — read bases on reference
        deletion_bases = 0  # D — reference bases not covered, non-intronic
        intron_bases = 0  # N — skipped (intron)
        soft_clip = 0
        del_lengths = []
        for op, length in read.cigartuples:
            if op in CIGAR_M_LIKE:
                aligned_bases += length
            elif op == CIGAR_D:
                deletion_bases += length
                del_lengths.append(length)
            elif op == CIGAR_N:
                intron_bases += length
            elif op == CIGAR_S:
                soft_clip += length

        read_length = read.infer_read_length()  # includes soft clips; None if unavailable
        if read_length is None:
            read_length = read.query_length
        ref_span = aligned_bases + deletion_bases + intron_bases
        denom = aligned_bases + deletion_bases
        missing_fraction = deletion_bases / denom if denom > 0 else 0.0

        rows.append([
            read.query_name,
            _tag_or_na(read, 'SM'),
            _tag_or_na(read, 'XT'),  # gene
            read.reference_name,
            read_length,
            aligned_bases,
            deletion_bases,
            len(del_lengths),
            ';'.join(map(str, del_lengths)),
            intron_bases,
            soft_clip,
            ref_span,
            missing_fraction,
            _tag_or_na(read, 'NR'),  # n reads stitched
            _tag_or_na(read, 'ER'),  # exonic reads
            _tag_or_na(read, 'IR'),  # intronic reads
            _tag_or_na(read, 'TC'),  # 3' read count
            _tag_or_na(read, 'IC'),  # internal read count
            _tag_or_na(read, 'FC'),  # 5' read count
        ])
    bam.close()
    return pd.DataFrame(rows, columns=[
        'read_name', 'sample', 'gene', 'reference',
        'read_length', 'aligned_bases', 'deletion_bases',
        'del_count', 'del_lengths',
        'intron_bases', 'soft_clip', 'ref_span', 'missing_fraction',
        'NR', 'ER', 'IR', 'TC', 'IC', 'FC',
    ])


def summary_stats(df):
    agg = df.groupby('sample').agg(
        n_reads=('read_name', 'count'),
        read_length_mean=('read_length', 'mean'),
        read_length_median=('read_length', 'median'),
        deletion_bases_mean=('deletion_bases', 'mean'),
        deletion_bases_median=('deletion_bases', 'median'),
        deletion_bases_sum=('deletion_bases', 'sum'),
        del_count_mean=('del_count', 'mean'),
        del_count_sum=('del_count', 'sum'),
        missing_fraction_mean=('missing_fraction', 'mean'),
        missing_fraction_median=('missing_fraction', 'median'),
        missing_fraction_p90=('missing_fraction', lambda s: np.percentile(s, 90)),
        reads_with_any_deletion=('deletion_bases', lambda s: int((s > 0).sum())),
    ).reset_index()
    return agg


def main():
    parser = argparse.ArgumentParser(
        description='Per-read read length and non-intronic missing bases (CIGAR D) from a BAM.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('-i', '--input', required=True, help='Input .bam file')
    parser.add_argument('-o', '--out_dir', required=True, help='Output directory')
    parser.add_argument('-p', '--prefix', required=True, help='Output filename prefix')
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    df = per_read_stats(args.input)
    df.to_csv(os.path.join(args.out_dir, f'{args.prefix}_per_read.tsv'),
              sep='\t', index=False)

    summary_stats(df).to_csv(os.path.join(args.out_dir, f'{args.prefix}_summary.tsv'),
                             sep='\t', index=False)


if __name__ == '__main__':
    main()
