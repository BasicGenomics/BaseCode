#!/usr/bin/env python3

import argparse
import pysam
import subprocess
import tempfile
import os

TAIL_LEN = 24

def main():
    parser = argparse.ArgumentParser(
        description='Generate sorted BED file with polyA tail positions from BAM'
    )
    parser.add_argument('-i', '--input', required=True, help='Input BAM file')
    parser.add_argument('-o', '--output', required=True, help='Output BED file')

    args = parser.parse_args()

    bam_in = pysam.AlignmentFile(args.input, 'rb')

    with tempfile.NamedTemporaryFile(delete=False, suffix=".bed") as tmp:
        tmp_bed = tmp.name

    try:
        with open(tmp_bed, "w") as bed:
            for read in bam_in.fetch(until_eof=True):

                if read.get_tag('TC') > 0:

                    chrom = bam_in.get_reference_name(read.reference_id)

                    gene   = read.get_tag('XT')
                    sample = read.get_tag('SM')

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