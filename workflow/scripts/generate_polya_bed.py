#!/usr/bin/env python3

import argparse
import pysam
import subprocess
import tempfile
import os

TAIL_LEN = 24


def passes_filter(read, mode):
    """Evaluate read against selected filtering mode"""

    try:
        TC = read.get_tag('TC')
        FC = read.get_tag('FC') if read.has_tag('FC') else 0
        IC = read.get_tag('IC') if read.has_tag('IC') else 0
    except KeyError:
        return False

    if mode == "default":
        return TC > 0

    elif mode == "full_length":
        return TC > 0 and FC > 0

    elif mode == "long_full_length":
        return TC > 0 and IC > 0 and FC > 0

    elif mode == "short_full_length":
        return TC > 0 and IC == 0 and FC > 0

    return False


def main():
    parser = argparse.ArgumentParser(
        description='Generate sorted BED file with polyA tail positions from BAM'
    )

    parser.add_argument('-i', '--input', required=True, help='Input BAM file')
    parser.add_argument('-o', '--output', required=True, help='Output BED file')

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

    bam_in = pysam.AlignmentFile(args.input, 'rb')

    with tempfile.NamedTemporaryFile(delete=False, suffix=".bed") as tmp:
        tmp_bed = tmp.name

    try:
        with open(tmp_bed, "w") as bed:

            for read in bam_in.fetch(until_eof=True):

                if not passes_filter(read, mode):
                    continue

                chrom = bam_in.get_reference_name(read.reference_id)

                gene   = read.get_tag('XT') if read.has_tag('XT') else "NA"
                sample = read.get_tag('SM') if read.has_tag('SM') else "NA"

                if read.is_reverse:
                    start = read.reference_start - TAIL_LEN
                    end   = read.reference_start
                    strand = "-"
                else:
                    start = read.reference_end
                    end   = read.reference_end + TAIL_LEN
                    strand = "+"

                if start < 0:
                    start = 0

                name = f"{read.query_name}|{gene}|{sample}"

                bed.write(
                    f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n"
                )

        subprocess.run(
            f"sort -k1,1V -k2,2n {tmp_bed} > {args.output}",
            shell=True,
            check=True
        )

    finally:
        if os.path.exists(tmp_bed):
            os.remove(tmp_bed)


if __name__ == '__main__':
    main()