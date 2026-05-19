import argparse
from tools import parse_gtf
import polars as pl

bed12_cols = [
    "chrom", "chromStart", "chromEnd", "transcript_id",
    "score", "strand", "thickStart", "thickEnd",
    "itemRgb", "blockCount", "blockSizes", "blockStarts"
]

def _add_canonical_len(bed, gff, longform,gene_identifier):
    # Read only transcript_id (col 3) and blockSizes (col 10) — skip the other 10 columns
    bed_df = (
        pl.read_csv(bed, separator="\t", has_header=False, null_values=".",
                    columns=[3, 10], new_columns=["transcript_id", "blockSizes"])
        .with_columns(
            pl.col('blockSizes')
              .str.split(',')
              .list.eval(pl.element().cast(pl.Int64, strict=False))
              .list.sum()
              .alias('canonical_transcripts_len')
        )
        .select(['transcript_id', 'canonical_transcripts_len'])
    )

    gene_dict = parse_gtf(gff, None, gene_identifier=gene_identifier, feature_type='transcript')
    mapping = {k: v['transcript_id'] for k, v in gene_dict.items() if v.get('transcript_id')}
    del gene_dict

    mapping_df = pl.DataFrame({
        'XT': list(mapping.keys()),
        'canonical_transcripts_id': list(mapping.values())
    })
    del mapping

    longform_df = pl.read_csv(longform, schema_overrides={"SEQID": pl.Utf8}).join(mapping_df, on='XT', how='left')

    return longform_df.join(bed_df, left_on='canonical_transcripts_id', right_on='transcript_id', how='inner')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Write coverage by canonical length into longform summary file', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--bed', help='bed file of canonical transcripts.')
    parser.add_argument('--gff', help='GFF file with canonical transcripts.')
    parser.add_argument('--longform', help='longform summary file')
    parser.add_argument('--gene_identifier', help='Gene identifier for parsing gff file. Either gene_name or gene_id.')

    args = parser.parse_args()

    out_df = _add_canonical_len(args.bed,args.gff,args.longform,args.gene_identifier)
    out_df.write_csv(args.longform)
