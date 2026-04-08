#!/usr/bin/env python3

import argparse
import pysam
import subprocess
import tempfile
import os
import csv
from collections import defaultdict

TAIL_LEN = {
    "polya": 24,
}

def passes_filter(read, mode, site_type):
    try:
        TC = read.get_tag('TC')
        FC = read.get_tag('FC') if read.has_tag('FC') else 0
        IC = read.get_tag('IC') if read.has_tag('IC') else 0
    except KeyError:
        return False
    
    if site_type == "polya" and TC == 0:
        return False
    if site_type == "tss" and FC == 0:
        return False

    if mode == "default":
        return True 
    elif mode == "full_length":
        return TC > 0 and FC > 0
    elif mode == "long_full_length":
        return TC > 0 and IC > 0 and FC > 0
    elif mode == "short_full_length":
        return TC > 0 and IC == 0 and FC > 0

    return False

def get_site_coords(read, site_type):
    """Return (start, end, strand) for the requested site type."""
    if site_type == "polya":
        # 3' end — window downstream of cleavage site
        tail = TAIL_LEN["polya"]
        if read.is_reverse:
            start = read.reference_start - tail
            end = read.reference_start
            strand = "-"
        else:
            start = read.reference_end
            end = read.reference_end + tail
            strand = "+"
    else:
        # TSS — single base at 5' end of molecule
        if read.is_reverse:
            start = read.reference_end - 1
            end = read.reference_end
            strand = "-"
        else:
            start = read.reference_start
            end = read.reference_start + 1
            strand = "+"
    start = max(0, start)
    return start, end, strand


def main():
    parser = argparse.ArgumentParser(
        description='Generate sorted BED file with polyA or TSS positions from BAM'
    )
    parser.add_argument('-i', '--input', required=True, 
                        help='Input BAM file')
    parser.add_argument('-o', '--output', required=True, 
                        help='Output BED file prefix')
    parser.add_argument('--site-type', choices=['polya', 'tss'], default='polya',
                        help='Site type to extract: polya (3\' end) or tss (5\' end) (default: polya)')
    parser.add_argument('--per-sample', action='store_true',
                        help='Generate one BED file per SM tag')
    parser.add_argument('--metadata',
                        help='Metadata file (tab-separated)')
    parser.add_argument('--group-by',
                        help='Metadata column to group samples by')
    parser.add_argument('--sample-id-column', default='ID',
                        help='Column in metadata matching BAM SM tag (default: ID)')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--full-length', action='store_true',
                       help='Full-length molecules only (TC>0 & FC>0)')
    group.add_argument('--long-full-length', action='store_true',
                       help='Long full-length molecules only (TC>0 & IC>0 & FC>0)')
    group.add_argument('--short-full-length', action='store_true',
                       help='Short full-length molecules only (TC>0 & IC==0 & FC>0)')

    args = parser.parse_args()

    mode = "default"
    if args.full_length:
        mode = "full_length"
    elif args.long_full_length:
        mode = "long_full_length"
    elif args.short_full_length:
        mode = "short_full_length"

    site_label = "polyA" if args.site_type == "polya" else "TSS"

    sample_to_group = {}
    group_to_samples = defaultdict(list)

    if args.group_by:
        if not args.metadata:
            raise ValueError("--group-by requires --metadata")

        delimiter = "," if args.metadata.endswith(".csv") else "\t"
        with open(args.metadata) as meta:
            reader = csv.DictReader(meta, delimiter=delimiter)
            if args.group_by not in reader.fieldnames:
                raise ValueError(f"{args.group_by} not found in metadata columns")
            if args.sample_id_column not in reader.fieldnames:
                raise ValueError(f"{args.sample_id_column} not found in metadata columns")
            for row in reader:
                sample_id = row[args.sample_id_column]
                group_value = row[args.group_by]
                sample_to_group[sample_id] = group_value
                group_to_samples[group_value].append(sample_id)

    bam_in = pysam.AlignmentFile(args.input, 'rb')

    sample_files = {}
    sample_handles = {}

    try:
        for read in bam_in.fetch(until_eof=True):
            if not passes_filter(read, mode, args.site_type):
                continue
            sample = read.get_tag('SM') if read.has_tag('SM') else "NA"
            gene = read.get_tag('XT') if read.has_tag('XT') else "NA"
            if args.group_by:
                if sample not in sample_to_group:
                    continue
                key = sample_to_group[sample]
            elif args.per_sample:
                key = sample
            else:
                key = "ALL"
            if key not in sample_files:
                tmp = tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=f"_{key}.bed")
                sample_files[key] = tmp.name
                sample_handles[key] = tmp
            bed = sample_handles[key]
            chrom = bam_in.get_reference_name(read.reference_id)
            start, end, strand = get_site_coords(read, args.site_type)
            name = f"{read.query_name}|{gene}|{sample}"
            bed.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n")
        for handle in sample_handles.values():
            handle.close()
        for key, tmp_bed in sample_files.items():
            if args.group_by:
                samples = sorted(group_to_samples[key])
                if len(samples) == 1:
                    final_output = f"{samples[0]}_{key}_{site_label}.bed"
                else:
                    final_output = f"{samples[0]}_to_{samples[-1]}_{key}_{site_label}.bed"
            elif args.per_sample:
                final_output = f"{args.output}_{key}_{site_label}.bed"
            else:
                final_output = f"{args.output}_{site_label}.bed"
            with open(final_output, "w") as out_fh:
                subprocess.run(
                    ["sort", "-k1,1V", "-k2,2n", tmp_bed],
                    stdout=out_fh,
                    check=True
                )
            os.remove(tmp_bed)
    finally:
        bam_in.close()

if __name__ == '__main__':
    main()