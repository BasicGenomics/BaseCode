import argparse
import math
import polars as pl
import numpy as np

def _parse_blocks(s):
    return np.array([int(x) for x in s.rstrip(",").split(",") if x.strip()], dtype=np.int64)


def _parse_cv(cv_str):
    """Parse CV tag 'start-end:count;...' → list of (start, end, count).
    end is inclusive in the tag; convert to exclusive here."""
    result = []
    for seg in cv_str.split(";"):
        if not seg:
            continue
        region, count = seg.split(":")
        start, end = region.split("-")
        result.append((int(start), int(end) + 1, int(count)))
    return result


def _m_block_to_transcript_intervals(blk_start, blk_end, abs_starts, abs_ends, cum, strand, canonical_len):
    """
    Map one genomic M block [blk_start, blk_end) to (t_lo, t_hi) transcript-coordinate intervals.
    One interval per BED12 exon the block overlaps. D gaps are excluded by not being M blocks.
    """
    intervals = []
    for i in range(len(abs_starts)):
        ol_s = max(blk_start, abs_starts[i])
        ol_e = min(blk_end,   abs_ends[i])
        if ol_s >= ol_e:
            continue
        t_lo = int(cum[i]) + int(ol_s) - int(abs_starts[i])
        t_hi = int(cum[i]) + int(ol_e) - int(abs_starts[i])   # exclusive
        if strand == "-":
            t_lo, t_hi = canonical_len - t_hi, canonical_len - t_lo
        intervals.append((t_lo, t_hi))
    return intervals


def _coverage_from_intervals(lo_arr, hi_arr, n_bins):
    """O(n) interval coverage via the diff trick."""
    diff = np.zeros(n_bins + 1, dtype=np.int64)
    np.add.at(diff, np.clip(lo_arr, 0, n_bins - 1), 1)
    np.add.at(diff, np.clip(hi_arr + 1, 0, n_bins), -1)
    return np.cumsum(diff[:n_bins])


def _weighted_coverage_from_intervals(lo_arr, hi_arr, weights_arr, n_bins):
    """Same diff trick but weighted — each interval contributes its weight instead of 1."""
    diff = np.zeros(n_bins + 1, dtype=np.float64)
    np.add.at(diff, np.clip(lo_arr, 0, n_bins - 1), weights_arr)
    np.add.at(diff, np.clip(hi_arr + 1, 0, n_bins), -weights_arr)
    return np.cumsum(diff[:n_bins])


def _to_bins(t_lo_expr, t_hi_expr, canonical_len_expr, n_bins):
    """Polars expressions: transcript coords → (b_lo, b_hi) clipped to [0, n_bins-1]."""
    raw_lo = (t_lo_expr * n_bins // canonical_len_expr).clip(0, n_bins - 1)
    raw_hi = ((t_hi_expr - 1) * n_bins // canonical_len_expr).clip(0, n_bins - 1)
    return pl.min_horizontal(raw_lo, raw_hi), pl.max_horizontal(raw_lo, raw_hi)


def compute_coverage(parquet, group_cols=None, n_bins=100, weighted=False):
    """
    Aggregate pre-computed bin intervals from parquet → per-group gene body coverage.
    Input parquet must be produced by precompute_bins().
    group_cols=None computes over the whole dataset (equivalent to RSeQC whole-BAM mode).
    Returns long-form DataFrame: <group_cols> | bin | coverage | coverage_norm
    and, if weighted=True: coverage_weighted | coverage_weighted_norm | mean_weight.
    """

    if isinstance(parquet,str):
        df = pl.read_parquet(parquet)
    else:
        df = parquet

    b_lo_expr, b_hi_expr = _to_bins(
        pl.col("mol_t_lo"), pl.col("mol_t_hi"), pl.col("canonical_len"), n_bins
    )
    sel_cols = (group_cols or []) + ["mol_t_lo", "mol_t_hi", "canonical_len"]
    mol_pdf = (
        df.select(sel_cols)
          .explode(["mol_t_lo", "mol_t_hi"])
          .filter(pl.col("canonical_len") > 0)
          .with_columns([b_lo_expr.alias("b_lo"), b_hi_expr.alias("b_hi")])
          .to_pandas()
    )

    if weighted:
        cv_b_lo, cv_b_hi = _to_bins(
            pl.col("cv_t_lo"), pl.col("cv_t_hi"), pl.col("canonical_len"), n_bins
        )
        cv_sel = (group_cols or []) + ["cv_t_lo", "cv_t_hi", "cv_count", "canonical_len"]
        cv_pdf = (
            df.select(cv_sel)
              .explode(["cv_t_lo", "cv_t_hi", "cv_count"])
              .filter(pl.col("canonical_len") > 0)
              .with_columns([cv_b_lo.alias("b_lo"), cv_b_hi.alias("b_hi")])
              .to_pandas()
        )

    if group_cols is None:
        total = len(df)
        groups = [((), mol_pdf)]
    else:
        totals = df.select(group_cols).to_pandas().groupby(group_cols).size()
        groups = [(key, grp) for key, grp in mol_pdf.groupby(group_cols)]

    rows = []
    for key, mol_grp in groups:
        key = (key,) if not isinstance(key, tuple) else key
        group_dict = dict(zip(group_cols or [], key))

        lo = mol_grp["b_lo"].to_numpy(dtype=np.int32)
        hi = mol_grp["b_hi"].to_numpy(dtype=np.int32)
        cov  = _coverage_from_intervals(lo, hi, n_bins)
        if group_cols is not None:
            total = int(totals.get(key[0] if len(key) == 1 else key, 1))
        norm = cov / total if total > 0 else cov.astype(float)

        if weighted:
            cv_grp = cv_pdf
            for col, val in group_dict.items():
                cv_grp = cv_grp[cv_grp[col] == val]
            if len(cv_grp):
                cov_w = _weighted_coverage_from_intervals(
                    cv_grp["b_lo"].to_numpy(dtype=np.int32),
                    cv_grp["b_hi"].to_numpy(dtype=np.int32),
                    cv_grp["cv_count"].to_numpy(dtype=np.float64), n_bins,
                )
            else:
                cov_w = np.zeros(n_bins, dtype=np.float64)
            total_w = cov_w.sum()
            norm_w  = cov_w / total_w if total_w > 0 else cov_w.copy()

        for b in range(n_bins):
            row = {**group_dict, "bin": int(b), "coverage": int(cov[b]), "coverage_norm": float(norm[b])}
            if weighted:
                row["coverage_weighted"]      = float(cov_w[b])
                row["coverage_weighted_norm"] = float(norm_w[b])
                row["mean_weight"]            = float(cov_w[b] / cov[b]) if cov[b] > 0 else float("nan")
            rows.append(row)

    return pl.DataFrame(rows)


def precompute_bins(longform, bed, weighted=False,CV_COL='read_depth_profile'):
    """
    Pre-compute per-molecule transcript bin intervals and write to parquet.
    Output has one row per molecule with list columns:
      mol_t_lo, mol_t_hi          — transcript coords of each M block
      cv_t_lo, cv_t_hi, cv_count  — transcript coords + counts from read_depth_profile (if weighted)
    Plus canonical_len and all metadata columns from the longform.
    """
    bed_df = pl.read_csv(
        bed, separator="\t", has_header=False, null_values=".",
        columns=[1, 3, 5, 10, 11],
        new_columns=["chromStart", "transcript_id", "strand", "blockSizes", "blockStarts"],
    )

    if isinstance(longform,str):
        df = pl.read_csv(longform, schema_overrides={"SEQID": pl.Utf8}).join(
        bed_df, left_on="canonical_transcripts_id", right_on="transcript_id", how="inner"
    )
    else:
        df = longform.join(
        bed_df, left_on="canonical_transcripts_id", right_on="transcript_id", how="inner"
    )
    
    pdf = df.to_pandas()

    def add_intervals(grp):
        r0 = grp.iloc[0]
        block_starts  = _parse_blocks(r0["blockStarts"])
        block_sizes   = _parse_blocks(r0["blockSizes"])
        abs_starts    = int(r0["chromStart"]) + block_starts
        abs_ends      = abs_starts + block_sizes
        cum           = np.concatenate([[0], np.cumsum(block_sizes)])
        strand        = r0["strand"]
        canonical_len = int(block_sizes.sum())

        def map_iv(start, end):
            return _m_block_to_transcript_intervals(
                start, end, abs_starts, abs_ends, cum, strand, canonical_len
            )

        grp = grp.copy()
        grp["canonical_len"] = canonical_len

        ivs = grp["aligned_blocks"].apply(lambda s: [
            iv for blk in s.split(";") if blk
            for iv in map_iv(*map(int, blk.split(":")))
        ])
        grp["mol_t_lo"] = ivs.apply(lambda x: [i[0] for i in x])
        grp["mol_t_hi"] = ivs.apply(lambda x: [i[1] for i in x])

        if weighted:
            def cv_ivs(row):
                offset = int(row["ref_start"])
                lo, hi, cnt = [], [], []
                for start, end, count in _parse_cv(str(row[CV_COL])):
                    for t_lo, t_hi in map_iv(start + offset, end + offset):
                        lo.append(t_lo); hi.append(t_hi); cnt.append(count)
                return lo, hi, cnt
            cv = grp[[CV_COL, "ref_start"]].apply(cv_ivs, axis=1)
            grp["cv_t_lo"]  = cv.apply(lambda x: x[0])
            grp["cv_t_hi"]  = cv.apply(lambda x: x[1])
            grp["cv_count"] = cv.apply(lambda x: x[2])

        return grp

    pdf = pdf.groupby("canonical_transcripts_id", group_keys=False).apply(add_intervals)

    # drop_after = ["aligned_blocks", "chromStart", "blockStarts", "blockSizes",
    #               "strand", "canonical_transcripts_id"] + (([CV_COL]) if weighted else [])
    # pdf = pdf.drop(columns=[c for c in drop_after if c in pdf.columns])

    return pl.from_pandas(pdf)


def _bed_genomic_percentile_positions(bed, n_bins=100, min_len=100):
    """
    For each transcript in BED12, pick n_bins evenly-spaced exonic genomic positions
    using the same sampling logic as RSeQC's percentile_list.
    Returns (chrom_pos_arr, chrom_bin_arr): per-chromosome sorted arrays of
    0-based genomic positions and their corresponding 0-indexed bin numbers.
    """
    from collections import defaultdict
    pos_by_chrom = defaultdict(list)

    with open(bed) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip().split("\t")
            if len(f) < 12:
                continue
            chrom = f[0]
            tx_start  = int(f[1])
            strand    = f[5]
            block_sizes  = _parse_blocks(f[10])
            block_starts = _parse_blocks(f[11])
            abs_starts = tx_start + block_starts
            abs_ends   = abs_starts + block_sizes

            all_bases = []
            for s, e in zip(abs_starts.tolist(), abs_ends.tolist()):
                all_bases.extend(range(int(s) + 1, int(e) + 1))  # 1-based, matches RSeQC

            if len(all_bases) < min_len:
                continue

            n = len(all_bases)
            # Match RSeQC percentile_list exactly: k = (n-1)*i/100 for i in 1..100,
            # with linear interpolation when k is non-integer.
            sampled = []
            for i in range(1, n_bins + 1):
                k = (n - 1) * i / n_bins
                f = math.floor(k)
                c = math.ceil(k)
                if f == c:
                    sampled.append(all_bases[int(k)] - 1)  # convert back to 0-based
                else:
                    gpos = int(round(all_bases[f] * (c - k) + all_bases[c] * (k - f)))
                    sampled.append(gpos - 1)  # convert back to 0-based
            if strand == "-":
                sampled = sampled[::-1]

            for b, gpos in enumerate(sampled):
                pos_by_chrom[chrom].append((gpos, b))

    chrom_pos_arr = {}
    chrom_bin_arr = {}
    for chrom, entries in pos_by_chrom.items():
        entries.sort(key=lambda x: x[0])
        chrom_pos_arr[chrom] = np.array([e[0] for e in entries], dtype=np.int64)
        chrom_bin_arr[chrom] = np.array([e[1] for e in entries], dtype=np.int32)

    return chrom_pos_arr, chrom_bin_arr


def precompute_bins_rseqc(longform, bed, n_bins=100, min_len=100):
    """
    RSeQC-equivalent precompute: for each molecule, find all transcript percentile bin
    indices it overlaps across ALL transcripts in the BED (not just its XT-assigned one).
    Duplicates are kept — a molecule covering the same percentile in N transcripts
    contributes N counts, matching RSeQC pileup semantics.
    Output parquet has a list column 'rseqc_bins' plus all metadata columns.
    """
    chrom_pos_arr, chrom_bin_arr = _bed_genomic_percentile_positions(bed, n_bins, min_len)

    if isinstance(longform, str):
        pdf = pl.read_csv(longform, schema_overrides={"SEQID": pl.Utf8}).to_pandas()
    else:
        pdf = longform.to_pandas()

    def mol_bins(row):
        seqid = str(row["SEQID"])
        if seqid not in chrom_pos_arr:
            return []
        pos_arr = chrom_pos_arr[seqid]
        bin_arr = chrom_bin_arr[seqid]
        hits = []
        for blk in str(row["aligned_blocks"]).split(";"):
            if not blk:
                continue
            s, e = map(int, blk.split(":"))
            lo = int(np.searchsorted(pos_arr, s, side="left"))
            hi = int(np.searchsorted(pos_arr, e, side="left"))
            if lo < hi:
                hits.extend(bin_arr[lo:hi].tolist())
        return hits

    drop = {"aligned_blocks", "del_lengths", "read_depth_profile"}
    keep = [c for c in pdf.columns if c not in drop]
    result = pdf[keep].copy()
    result["rseqc_bins"] = pdf.apply(mol_bins, axis=1)
    result = result[result["rseqc_bins"].apply(len) > 0].reset_index(drop=True)
    return pl.from_pandas(result)


def compute_coverage_rseqc(parquet, group_cols=None, n_bins=100):
    """
    Aggregate RSeQC-style precomputed bin hits → per-group gene body coverage.
    Input parquet must be from precompute_bins_rseqc().
    bin is 1-indexed (1–n_bins) to match RSeQC output convention.
    """

    if isinstance(parquet,str):
        df = pl.read_parquet(parquet)
    else:
        df = parquet

    sel_cols = (group_cols or []) + ["rseqc_bins"]
    flat = (
        df.select(sel_cols)
          .explode("rseqc_bins")
          .to_pandas()
          .rename(columns={"rseqc_bins": "bin"})
          .dropna(subset=["bin"])
    )

    if group_cols is None:
        total = len(df)
        groups = [((), flat)]
    else:
        totals = df.select(group_cols).to_pandas().groupby(group_cols).size()
        groups = [(key, grp) for key, grp in flat.groupby(group_cols)]

    rows = []
    for key, grp in groups:
        key = (key,) if not isinstance(key, tuple) else key
        group_dict = dict(zip(group_cols or [], key))

        cov = np.zeros(n_bins, dtype=np.int64)
        np.add.at(cov, np.clip(grp["bin"].to_numpy(dtype=np.int32), 0, n_bins - 1), 1)

        if group_cols is not None:
            total = int(totals.get(key[0] if len(key) == 1 else key, 1))
        norm = cov / total if total > 0 else cov.astype(float)

        for b in range(n_bins):
            rows.append({**group_dict, "bin": b + 1, "coverage": int(cov[b]), "coverage_norm": float(norm[b])})

    return pl.DataFrame(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Gene body coverage tools.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_pre = sub.add_parser("precompute", help="Pre-compute per-molecule transcript bin intervals → parquet.")
    p_pre.add_argument("--longform", required=True)
    p_pre.add_argument("--bed", required=True)
    p_pre.add_argument("--output", required=True, help="Output .parquet file.")
    p_pre.add_argument("--weighted", action="store_true")

    p_cov = sub.add_parser("coverage", help="Aggregate parquet → per-group gene body coverage CSV.")
    p_cov.add_argument("--parquet", required=True, help="Parquet from precompute step.")
    p_cov.add_argument("--output", required=True, help="Output CSV.")
    p_cov.add_argument("--groupby", default=None,
                       help="Comma-separated column names to group by (e.g. SM,XT or SM). "
                            "Omit for whole-dataset mode (equivalent to RSeQC whole-BAM).")
    p_cov.add_argument("--n_bins", type=int, default=100)
    p_cov.add_argument("--weighted", action="store_true")

    p_pre_r = sub.add_parser("precompute_rseqc", help="RSeQC-equivalent: precompute per-molecule bin hits across all BED transcripts → parquet.")
    p_pre_r.add_argument("--longform", required=True)
    p_pre_r.add_argument("--bed", required=True)
    p_pre_r.add_argument("--output", required=True, help="Output .parquet file.")
    p_pre_r.add_argument("--n_bins", type=int, default=100)
    p_pre_r.add_argument("--min_len", type=int, default=100, help="Skip transcripts shorter than this (bp).")

    p_cov_r = sub.add_parser("coverage_rseqc", help="RSeQC-equivalent: aggregate parquet → per-group coverage CSV (1-indexed bins).")
    p_cov_r.add_argument("--parquet", required=True, help="Parquet from precompute_rseqc step.")
    p_cov_r.add_argument("--output", required=True, help="Output CSV.")
    p_cov_r.add_argument("--groupby", default=None,
                         help="Comma-separated column names to group by. Omit for whole-dataset mode.")
    p_cov_r.add_argument("--n_bins", type=int, default=100)

    args = parser.parse_args()

    if args.cmd == "precompute":
        precompute_bins(args.longform, args.bed, args.weighted).write_parquet(args.output)
    elif args.cmd == "coverage":
        group_cols = [c.strip() for c in args.groupby.split(",")] if args.groupby else None
        compute_coverage(args.parquet, group_cols, args.n_bins, args.weighted).write_csv(args.output)
    elif args.cmd == "precompute_rseqc":
        precompute_bins_rseqc(args.longform, args.bed, args.n_bins, args.min_len).write_parquet(args.output)
    else:
        group_cols = [c.strip() for c in args.groupby.split(",")] if args.groupby else None
        compute_coverage_rseqc(args.parquet, group_cols, args.n_bins).write_csv(args.output)
