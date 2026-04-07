import re
import io
from pathlib import Path
import argparse
import math
import pandas as pd
import numpy as np
import polars as pl
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import matplotlib.colors as mcolors
import seaborn as sns
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, PageBreak, Image, Table, TableStyle
from reportlab.lib.enums import TA_LEFT, TA_CENTER
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import cm
from reportlab.lib import colors

## Options
sns.set_theme(style="whitegrid", rc={
    "text.color": "black",
    "axes.labelcolor": "black",
    "axes.edgecolor": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "axes.titlecolor": "black",
    "legend.edgecolor": "black",
})

def main(args):
    samplesheet_path = Path(args.samplesheet)
    summary_path = Path(args.summary)
    long_form_summary_path = Path(args.long_form_summary)
    outdir = Path(args.outdir)
    for p in [samplesheet_path, summary_path, long_form_summary_path]:
        if not p.exists():
            raise FileNotFoundError(f"File not found: {p}")
    background_pdf_path = Path(args.background_pdf) if args.background_pdf else None
    if background_pdf_path and not background_pdf_path.exists():
        print(f"Warning: Background image not found: {background_pdf_path}, skipping.")
        background_pdf_path = None
    readtype_schema_path = Path(args.readtype_schema) if args.readtype_schema else None
    if readtype_schema_path and not readtype_schema_path.exists():
        print(f"Warning: Read type schema not found: {readtype_schema_path}, skipping.")
        readtype_schema_path = None
    # moleculetype_schema_path = Path(args.moleculetype_schema)
    # if not moleculetype_schema_path.exists():
    #     raise FileNotFoundError(f"Molecule type schema not found: {moleculetype_schema_path}")
    outdir.mkdir(parents=True, exist_ok=True)
    run_id_samplesheet = samplesheet_path.stem.split("_")[0] 
    run_id = args.run_id if args.run_id else run_id_samplesheet
    out_pdf = outdir / f"{run_id}_run_report.pdf"
    if args.verbose:
        print("Generating BaseCode Processing Pipeline Run Report")

    ## Metadata
    meta = pd.read_csv(samplesheet_path)
    tab_samples = meta[['SAMPLE_ID', 'DESC']]
    tab_samples = tab_samples.rename(columns={
        "SAMPLE_ID": "Sample",
        "DESC": "Description"
    })
    samples = tab_samples["Sample"].tolist()
    n = len(samples)

    def natural_key(x):
        return [int(text) if text.isdigit() else text for text in re.split(r'(\d+)', x)]

    ## Summary stats
    summary = pd.read_csv(summary_path)
    summary["index"] = summary["index"].astype(str)
    summary = summary.sort_values(by="index", key=lambda col: col.map(natural_key))
    df = summary.copy()

    ### Plot #1
    def millions(x, pos):
        return f"{x/1e6:.0f}M"
    colors_counts = {
        "3' Counts": "#583092",
        "Internal Counts": "#97B3D6",
        "5' Counts": "#EC008C"
    }
    samples = df["index"]
    c5 = df["5' Counts"]
    cI = df["Internal Counts"]
    c3 = df["3' Counts"]
    fig, ax = plt.subplots(figsize=(7, max(2, min(14, len(samples)*0.5))))
    fig.patch.set_alpha(0)
    ax.barh(samples, c5, label="5' Reads", color=colors_counts["5' Counts"])
    ax.barh(samples, cI, left=c5, label="Internal Reads", color=colors_counts["Internal Counts"])
    ax.barh(samples, c3, left=c5 + cI, label="3' Reads", color=colors_counts["3' Counts"])
    ax.set_ylabel("")
    ax.set_xlabel("Read Counts", fontsize=11)
    ax.invert_yaxis()
    ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False, fontsize=10)
    ax.xaxis.set_major_formatter(FuncFormatter(millions))
    plt.grid(axis="y")
    for i, (v5, vI, v3) in enumerate(zip(c5, cI, c3)):
        if v5 > 0:
            ax.text(v5/2, i, f"{v5/1e6:.1f}M",
                    va="center", ha="center", color="white", fontsize=9)
        if vI > 0:
            ax.text(v5 + vI/2, i, f"{vI/1e6:.1f}M",
                    va="center", ha="center", color="black", fontsize=9)
        if v3 > 0:
            ax.text(v5 + vI + v3/2, i, f"{v3/1e6:.1f}M",
                    va="center", ha="center", color="white", fontsize=9)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.5)
    ax.tick_params(axis="x", labelsize=10)
    ax.tick_params(axis="y", labelsize=11)
    plt.tight_layout()
    fig_counts_sns = fig

    ### Plot #1 stats
    counts_min = df['Total Counts'].min()
    counts_max = df['Total Counts'].max()
    counts_min_M = math.ceil(counts_min / 1e5) / 10  
    counts_max_M = math.ceil(counts_max / 1e5) / 10

    ### Plot #2
    def millions(x, pos):
        v = x / 1e6
        if abs(v - round(v)) < 1e-9:
            return f"{int(round(v))}M"
        return f"{v:.1f}M"
    compartments = [
        "5' Transcript Region",
        "Internal Transcript Region",
        "3' Transcript Region"
    ]
    colors_counts = {
        "3' Transcript Region": "#583092",
        "Internal Transcript Region": "#97B3D6",
        "5' Transcript Region": "#EC008C"
    }
    plot_df = (
        df[["index",
            "Molecules Detected (5')",
            "Molecules Detected (INT)",
            "Molecules Detected (3')"]]
        .rename(columns={
            "Molecules Detected (5')": "5' Transcript Region",
            "Molecules Detected (INT)": "Internal Transcript Region",
            "Molecules Detected (3')": "3' Transcript Region"
        })
        .melt(id_vars="index", var_name="compartment", value_name="molecules")
    )
    plot_df["compartment"] = pd.Categorical(
        plot_df["compartment"],
        categories=compartments,
        ordered=True
    )
    unique_samples = plot_df["index"].unique()
    n = len(unique_samples)
    if n > 8:
        split = int(np.ceil(n / 2))
        chunks = [unique_samples[:split], unique_samples[split:]]
        ncols = 2
    else:
        chunks = [unique_samples]
        ncols = 1
    fig, axes = plt.subplots(
        ncols=ncols,
        figsize=(4 * ncols, max(2, min(14, n * 0.75))),
        sharex=True
    )
    if ncols == 1:
        axes = [axes]
    fig.patch.set_alpha(0)
    xmax = plot_df["molecules"].max() * 1.1
    for i, (ax, chunk) in enumerate(zip(axes, chunks)):
        sub_df = plot_df[plot_df["index"].isin(chunk)]
        sns.barplot(
            data=sub_df,
            y="index",
            x="molecules",
            hue="compartment",
            palette=colors_counts,
            hue_order=compartments,
            saturation=1,
            ax=ax,
            width=0.75
        )
        ax.set_ylabel("")
        ax.set_xlabel("")
        ax.set_xlim(0, xmax)
        ax.xaxis.set_major_formatter(FuncFormatter(millions))
        ax.legend_.remove()
        for container in ax.containers:
            for bar in container:
                w = bar.get_width()
                if w <= 0:
                    continue
                ax.text(
                    w + xmax * 0.01,
                    bar.get_y() + bar.get_height()/2,
                    f"{w/1e6:.2f}M",
                    va="center",
                    ha="left",
                    fontsize=9
                )
        for j in range(len(chunk)):
            if j % 2 == 0:
                ax.axhspan(j - 0.5, j + 0.5, color="grey", alpha=0.05)
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color("black")
            spine.set_linewidth(1.5)
        ax.tick_params(axis="x", labelsize=10)
        ax.tick_params(axis="y", labelsize=11)
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.55, -0.05), ncol=3, frameon=False, fontsize=10)
    fig.supxlabel("Molecules Detected", x=0.55, fontsize=11)
    plt.subplots_adjust(wspace=0.3)
    plt.tight_layout()
    fig_molecules_sns = fig

    ### Plot #2 stats
    molecules_5_mean = df["Molecules Detected (5')"].mean()
    molecules_int_mean = df["Molecules Detected (INT)"].mean()
    molecules_3_mean = df["Molecules Detected (3')"].mean()
    molecules_5_mean_M = round(molecules_5_mean / 1e6, 2)
    molecules_int_mean_M = round(molecules_int_mean / 1e6, 2)
    molecules_3_mean_M = round(molecules_3_mean / 1e6, 2)

    ### Plot #3
    def kilos(x, pos):
        return f"{x/1e3:.0f}k"
    fig, ax = plt.subplots(figsize=(max(4,  min(14, len(df) * 0.5)), 5))
    fig.patch.set_alpha(0)
    sns.barplot(
        x=df["index"],
        y=df["Genes Detected"],
        color="#583092",
        ax=ax,
        saturation=1
    )
    ax.set_xlabel("")
    ax.set_ylabel("Genes Detected", fontsize=11)
    plt.xticks(rotation=90, ha="center")
    ax.yaxis.set_major_formatter(FuncFormatter(kilos))
    for container in ax.containers:
        ax.bar_label(
            container,
            labels=[f"{v/1e3:.1f}k" for v in container.datavalues],
            padding=1,
            fontsize=9,
            color="black"
        )
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.5)
    ax.tick_params(axis='x', length=0)
    ax.tick_params(axis="x", labelsize=11)
    ax.tick_params(axis="y", labelsize=10)
    plt.tight_layout()
    fig_genes_sns = fig

    ### Plot #3 stats
    genes_detected_mean = df["Genes Detected"].mean() 
    genes_detected_mean_k = round(genes_detected_mean / 1e3, 1)

    ## Long summary stats
    df_lf = pl.read_csv(long_form_summary_path)
    sample_list = []
    genes_detected_list = []
    molecules_detected_total_list = []
    molecules_detected_TP_list = []
    molecules_detected_INT_list = []
    molecules_detected_FP_list = []
    molecules_detected_TP_FP_list = []
    molecules_detected_TP_INT_FP_list = []
    molecules_detected_TP_FP_noINT_list = []
    molecules_detected_TP_only_list = []
    molecules_detected_INT_only_list = []
    molecules_detected_FP_only_list = []
    molecules_detected_TP_INT_only_list = []
    molecules_detected_FP_INT_only_list = []
    percentage_long_full_length_list = []
    percentage_short_full_length_list = []
    percentage_reconstructed_end_complete_list = []
    percentage_reconstructed_long_full_length_list = []
    percentage_reconstructed_short_full_length_list = []
    median_reconstructed_end_complete_list = []
    median_reconstructed_long_full_length_list = []
    median_reconstructed_short_full_length_list = []
    for (unique_sample_id, df_long_form_sample) in df_lf.group_by('SM'):
        sample_list.extend(unique_sample_id)
        genes_detected_sample = df_long_form_sample.unique(subset=['XT']).shape[0]
        genes_detected_list.append(genes_detected_sample)
        molecules_detected_total_sample = df_long_form_sample.shape[0]
        molecules_detected_total_list.append(molecules_detected_total_sample)
        df_long_form_sample_threep = df_long_form_sample.filter(pl.col('TC') > 0)
        df_long_form_sample_int = df_long_form_sample.filter(pl.col('IC') > 0)
        df_long_form_sample_fivep = df_long_form_sample.filter(pl.col('FC') > 0)
        molecules_detected_TP_sample = df_long_form_sample_threep.shape[0]
        molecules_detected_TP_list.append(molecules_detected_TP_sample)
        molecules_detected_INT_sample = df_long_form_sample_int.shape[0]
        molecules_detected_INT_list.append(molecules_detected_INT_sample)
        molecules_detected_FP_sample = df_long_form_sample_fivep.shape[0]
        molecules_detected_FP_list.append(molecules_detected_FP_sample)
        df_long_form_end_complete_sample = df_long_form_sample.filter(((pl.col('TC') > 0) & (pl.col('FC') > 0)))
        df_long_form_long_full_length_sample = df_long_form_sample.filter(((pl.col('TC') > 0) & (pl.col('IC') > 0) & (pl.col('FC') > 0)))
        df_long_form_short_full_length_sample = df_long_form_sample.filter(((pl.col('TC') > 0) & (pl.col('IC') == 0) & (pl.col('FC') > 0)))
        molecules_detected_TP_FP_sample = df_long_form_end_complete_sample.shape[0]
        molecules_detected_TP_FP_list.append(molecules_detected_TP_FP_sample)
        molecules_detected_TP_INT_FP_sample = df_long_form_long_full_length_sample.shape[0]
        molecules_detected_TP_INT_FP_list.append(molecules_detected_TP_INT_FP_sample)
        molecules_detected_TP_FP_noINT_sample = df_long_form_short_full_length_sample.shape[0]
        molecules_detected_TP_FP_noINT_list.append(molecules_detected_TP_FP_noINT_sample)
        df_long_form_TP_only_sample = df_long_form_sample.filter(((pl.col('TC') > 0) & (pl.col('IC') == 0) & (pl.col('FC') == 0)))
        df_long_form_INT_only_sample = df_long_form_sample.filter(((pl.col('TC') == 0) & (pl.col('IC') > 0) & (pl.col('FC') == 0)))
        df_long_form_FP_only_sample = df_long_form_sample.filter(((pl.col('TC') == 0) & (pl.col('IC') == 0) & (pl.col('FC') > 0)))
        df_long_form_TP_INT_only_sample = df_long_form_sample.filter(((pl.col('TC') > 0) & (pl.col('IC') > 0) & (pl.col('FC') == 0)))
        df_long_form_FP_INT_only_sample = df_long_form_sample.filter(((pl.col('TC') == 0) & (pl.col('IC') > 0) & (pl.col('FC') > 0)))
        molecules_detected_TP_only_sample = df_long_form_TP_only_sample.shape[0]
        molecules_detected_TP_only_list.append(molecules_detected_TP_only_sample)
        molecules_detected_INT_only_sample = df_long_form_INT_only_sample.shape[0]
        molecules_detected_INT_only_list.append(molecules_detected_INT_only_sample)
        molecules_detected_FP_only_sample = df_long_form_FP_only_sample.shape[0]
        molecules_detected_FP_only_list.append(molecules_detected_FP_only_sample)
        molecules_detected_TP_INT_only_sample = df_long_form_TP_INT_only_sample.shape[0]
        molecules_detected_TP_INT_only_list.append(molecules_detected_TP_INT_only_sample)
        molecules_detected_FP_INT_only_sample = df_long_form_FP_INT_only_sample.shape[0]
        molecules_detected_FP_INT_only_list.append(molecules_detected_FP_INT_only_sample)
        percentage_long_full_length_sample = 100 * df_long_form_long_full_length_sample.shape[0] / df_long_form_end_complete_sample.shape[0]
        percentage_long_full_length_list.append(round(percentage_long_full_length_sample, 2))
        percentage_short_full_length_sample = 100 * df_long_form_short_full_length_sample.shape[0] / df_long_form_end_complete_sample.shape[0]
        percentage_short_full_length_list.append(round(percentage_short_full_length_sample, 2))
        percentage_reconstructed_end_complete_sample = np.round(100*(df_long_form_end_complete_sample.shape[0]/df_long_form_sample_threep.shape[0]),2)
        percentage_reconstructed_end_complete_list.append(percentage_reconstructed_end_complete_sample)
        percentage_reconstructed_long_full_length_sample = np.round(100*(df_long_form_long_full_length_sample.shape[0]/df_long_form_sample_threep.shape[0]),2)
        percentage_reconstructed_long_full_length_list.append(percentage_reconstructed_long_full_length_sample)
        percentage_reconstructed_short_full_length_sample = np.round(100*(df_long_form_short_full_length_sample.shape[0]/df_long_form_sample_threep.shape[0]),2)
        percentage_reconstructed_short_full_length_list.append(percentage_reconstructed_short_full_length_sample)
        median_reconstructed_end_complete_sample = int(df_long_form_end_complete_sample.median()['QL'][0])
        median_reconstructed_end_complete_list.append(median_reconstructed_end_complete_sample)
        median_reconstructed_long_full_length_sample = int(df_long_form_long_full_length_sample.median()['QL'][0])
        median_reconstructed_long_full_length_list.append(median_reconstructed_long_full_length_sample)
        median_reconstructed_short_full_length_sample = int(df_long_form_short_full_length_sample.median()['QL'][0])
        median_reconstructed_short_full_length_list.append(median_reconstructed_short_full_length_sample)
    data = {'index': sample_list, 'Genes Detected': genes_detected_list, 'Total Molecules Detected': molecules_detected_total_list, 
            'Molecules Detected (5\')': molecules_detected_FP_list, 'Molecules Detected (INT)': molecules_detected_INT_list, 'Molecules Detected (3\')': molecules_detected_TP_list, 
            'Molecules Detected (5\' only)': molecules_detected_FP_only_list, 'Molecules Detected (INT only)': molecules_detected_INT_only_list, 'Molecules Detected (3\' only)': molecules_detected_TP_only_list,
            'Molecules Detected (5\' and INT only)': molecules_detected_FP_INT_only_list, 'Molecules Detected (3\' and INT only)': molecules_detected_TP_INT_only_list,
            'Molecules Detected (5\' and 3\' with INT)': molecules_detected_TP_INT_FP_list, 'Molecules Detected (5\' and 3\' without INT)': molecules_detected_TP_FP_noINT_list, 'Molecules Detected (5\' and 3\' with/without INT)': molecules_detected_TP_FP_list, 
            'Percentage of Molecules with 5\' and 3\' containing INT (%)': percentage_long_full_length_list, 'Percentage of Molecules with 5\' and 3\' not containing INT (%)': percentage_short_full_length_list,
            'Percentage Completed (5\' and 3\') Relative to All 3\' Molecules (%)': percentage_reconstructed_end_complete_list, 'Percentage Completed (5\' and 3\' with INT) Relative to All 3\' Molecules (%)': percentage_reconstructed_long_full_length_list, 'Percentage Completed (5\' and 3\' without INT) Relative to All 3\' Molecules (%)': percentage_reconstructed_short_full_length_list, 
            'Median Length Completed Molecules (5\' and 3\') (bp)': median_reconstructed_end_complete_list,  'Median Length Completed Molecules (5\' and 3\' with INT) (bp)': median_reconstructed_long_full_length_list,  'Median Length Completed Molecules (5\' and 3\' without INT) (bp)': median_reconstructed_short_full_length_list
            }
    df_out = pd.DataFrame(data)
    df_out = df_out.sort_values(by="index")
    molecules_5_int_3_mean = round(df_out["Percentage Completed (5' and 3' with INT) Relative to All 3' Molecules (%)"].mean(), 1)
    molecules_5_noint_3_mean = round(df_out["Percentage Completed (5' and 3' without INT) Relative to All 3' Molecules (%)"].mean(), 1)
    molecules_5_int_3_length = round(math.ceil(df_out["Median Length Completed Molecules (5' and 3' with INT) (bp)"].mean()) / 10) * 10
    molecules_5_noint_3_length = round(math.ceil(df_out["Median Length Completed Molecules (5' and 3' without INT) (bp)"].mean()) / 10) * 10

    ### Plot #4
    def millions(x, pos):
        v = x / 1e6
        if abs(v - round(v)) < 1e-9:
            return f"{int(round(v))}M"
        return f"{v:.1f}M"
    def lighten(color, amount=0.5):
        c = np.array(mcolors.to_rgb(color))
        return tuple(c + (1 - c) * amount)
    categories = [
        "Molecules Detected (5' only)",
        "Molecules Detected (INT only)",
        "Molecules Detected (3' only)",
        "Molecules Detected (5' and INT only)",
        "Molecules Detected (3' and INT only)",
        "Molecules Detected (5' and 3' without INT)",
        "Molecules Detected (5' and 3' with INT)",
    ]
    label_map = {
        "Molecules Detected (5' only)": "5' —",
        "Molecules Detected (INT only)": "— Internal —",
        "Molecules Detected (3' only)": "— 3'",
        "Molecules Detected (5' and INT only)": "5' — Internal",
        "Molecules Detected (3' and INT only)": "Internal — 3'",
        "Molecules Detected (5' and 3' without INT)": "5' — 3'",
        "Molecules Detected (5' and 3' with INT)": "5' — Internal — 3'",
    }
    order = [
        "5' —",
        "5' — Internal",
        "5' — 3'",
        "— Internal —",
        "5' — Internal — 3'",
        "Internal — 3'",
        "— 3'"
    ]
    color_map = {
        "5' —": lighten("#EC008C", 0.8),
        "5' — Internal": lighten("#EC008C", 0.6),
        "5' — 3'": "#EC008C",
        "— Internal —": "#97B3D6",
        "5' — Internal — 3'": "#583092",
        "Internal — 3'": lighten("#583092", 0.4),
        "— 3'": lighten("#583092", 0.6),
    }
    plot_df = (
        df_out[["index"] + categories]
        .melt(id_vars="index", var_name="compartment", value_name="molecules")
    )
    plot_df["compartment"] = plot_df["compartment"].map(label_map)
    plot_df["molecules"] = pd.to_numeric(plot_df["molecules"], errors="coerce").fillna(0)
    plot_df["compartment"] = pd.Categorical(
        plot_df["compartment"], categories=order, ordered=True
    )
    samples = sorted(plot_df["index"].unique(), key=natural_key)[::-1]
    plot_df["index"] = pd.Categorical(
        plot_df["index"], categories=samples, ordered=True
    )
    plot_df["pct"] = 100 * plot_df["molecules"] / plot_df.groupby("index")["molecules"].transform("sum")
    stack_df = (
        plot_df.pivot(index="index", columns="compartment", values="molecules")
        .reindex(samples)
        .fillna(0)
    )
    fig = plt.figure(figsize=(7, max(4, min(14, len(samples)*0.75))))
    gs = fig.add_gridspec(1, 2, width_ratios=[4, 1], wspace=0.05)
    ax = fig.add_subplot(gs[0])
    ax_bar = fig.add_subplot(gs[1])
    sns.scatterplot(
        data=plot_df,
        y="index",
        x="compartment",
        hue="compartment",
        palette=color_map,
        size="pct",
        sizes=(10, 250),
        linewidth=0,
        ax=ax,
        legend=False
    )
    for _, r in plot_df.iterrows():
        if r["molecules"] <= 0:
            continue
        ax.annotate(
            f"{r['molecules']/1e6:.2f}M\n({r['pct']:.1f}%)",
            (r["compartment"], r["index"]),
            xytext=(0, 10),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=7
        )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=90, labelsize=11)
    ax.tick_params(axis="y", labelsize=11)
    ax.set_yticks(np.arange(len(samples)))
    ax.set_yticklabels(samples)
    ax.invert_yaxis()
    ax.grid(False)
    ax.set_xlim(-0.5, len(order)-0.5)
    left = np.zeros(len(samples))
    for cat in order:
        vals = stack_df[cat].values
        ax_bar.barh(
            samples,
            vals,
            left=left,
            color=color_map[cat],
            height=0.8
        )
        left += vals
    totals = stack_df.sum(axis=1).values
    max_val = stack_df.sum(axis=1).max()
    for i, total in enumerate(totals):
        if total <= 0:
            continue
        ax_bar.text(
            total + max_val * 0.01,
            i,
            f"{total/1e6:.2f}M",
            va="center",
            ha="left",
            fontsize=9,
            color="black"
        )
    ax_bar.set_xlim(0, max_val * 1.05)
    ax_bar.xaxis.set_major_formatter(FuncFormatter(millions))
    ax_bar.set_xlabel("")
    ax_bar.tick_params(axis="x", labelsize=10)
    ax_bar.set_yticks([])
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)
    ax_bar.set_ylim(ax.get_ylim())
    xticks = ax.get_xticks()
    xticklabels = [t.get_text() for t in ax.get_xticklabels()]
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.5)
    plt.subplots_adjust(left=0.1, right=1)
    fig_molecules_ext_sns = fig

    ### Plot #5
    def lighten(color, amount=0.5):
        c = np.array(mcolors.to_rgb(color))
        return tuple(c + (1 - c) * amount)
    categories = [
    "Molecules Detected (5' and 3' without INT)",
    "Molecules Detected (5' and 3' with INT)", 
    "Molecules Detected (3' and INT only)",
    "Molecules Detected (3' only)"
    ]
    label_map = {
        "Molecules Detected (5' and 3' without INT)": "5' — 3'",
        "Molecules Detected (5' and 3' with INT)": "5' — Internal — 3'",
        "Molecules Detected (3' and INT only)": "Internal — 3'",
        "Molecules Detected (3' only)": "— 3'",
    }
    order = list(label_map.values())
    color_map = {
        "5' — 3'": "#EC008C",
        "5' — Internal — 3'": "#583092",
        "Internal — 3'": lighten("#583092", 0.4),
        "— 3'": lighten("#583092", 0.6),
    }
    plot_df = (
        df_out[["index"] + categories]
        .melt(id_vars="index", var_name="compartment", value_name="molecules")
    )
    plot_df["compartment"] = plot_df["compartment"].map(label_map)
    plot_df["molecules"] = pd.to_numeric(plot_df["molecules"], errors="coerce").fillna(0)
    plot_df["compartment"] = pd.Categorical(
        plot_df["compartment"], categories=order, ordered=True
    )
    samples = sorted(plot_df["index"].unique(), key=natural_key)[::-1]
    plot_df["index"] = pd.Categorical(
        plot_df["index"], categories=samples, ordered=True
    )
    plot_df["pct"] = 100 * plot_df["molecules"] / plot_df.groupby("index")["molecules"].transform("sum")
    stack_df = (
        plot_df.pivot(index="index", columns="compartment", values="pct")
        .reindex(samples)
        .fillna(0)
    )
    fig = plt.figure(figsize=(7, max(4, min(14, len(samples)*0.75))))
    gs = fig.add_gridspec(1, 2, width_ratios=[3, 1], wspace=0.05)
    ax = fig.add_subplot(gs[0])
    ax_bar = fig.add_subplot(gs[1])
    sns.scatterplot(
        data=plot_df,
        y="index",
        x="compartment",
        hue="compartment",
        palette=color_map,
        size="pct",
        sizes=(10, 250),
        linewidth=0,
        ax=ax,
        legend=False
    )
    for _, r in plot_df.iterrows():
        ax.annotate(
            f"{r.pct:.1f}%",
            (r["compartment"], r["index"]),
            xytext=(0, 10),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=9
        )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.grid(False)
    ax.tick_params(axis="x", rotation=90, labelsize=11)
    ax.tick_params(axis="y", labelsize=11)
    ax.set_yticks(np.arange(len(samples)))
    ax.set_yticklabels(samples)
    ax.invert_yaxis()
    ax.set_xlim(-0.5, len(order)-0.5)
    left = np.zeros(len(samples))
    for cat in order:
        vals = stack_df[cat].values
        ax_bar.barh(
            samples,
            vals,
            left=left,
            color=color_map[cat],
            height=0.8
        )
        left += vals
    ax_bar.set_xlim(0, 100)
    ax_bar.set_xlabel("")
    ax_bar.set_xticks([0, 50, 100])
    ax_bar.set_xticklabels([f"{x}%" for x in [0, 50, 100]])
    ax_bar.tick_params(axis="x", labelsize=10)
    ax_bar.set_yticks([])
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)
    ax_bar.set_ylim(ax.get_ylim())
    xticks = ax.get_xticks()
    xticklabels = [t.get_text() for t in ax.get_xticklabels()]
    def get_xpos(label):
        return xticks[xticklabels.index(label)]
    ax.text(
        get_xpos("5' — Internal — 3'"),
        len(samples) + 0.15,
        "Long\nFull-Length",
        ha="center",
        va="bottom",
        fontsize=11,
        color="#583092"
    )
    ax.text(
        get_xpos("5' — 3'"),
        len(samples) + 0.15,
        "Short\nFull-Length",
        ha="center",
        va="bottom",
        fontsize=11,
        color="#EC008C"
    )
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.5)
    plt.subplots_adjust(left=0.1, right=1)
    fig_molecules_ext_full_length_sns = fig

    ### Plot #6
    def rain_panel(ax, df_lf, samples, title, completed_filter, color):
        records = []
        medians = {}
        for sm in samples:
            df_sample = df_lf.filter(pl.col("SM") == sm)
            df_completed = df_sample.filter(completed_filter)
            ql_df = df_completed.select("QL").to_pandas()
            ql_df["SM"] = sm
            records.append(ql_df)
            medians[sm] = ql_df["QL"].median() if len(ql_df) > 0 else np.nan
        df_rain = pd.concat(records, ignore_index=True)
        order = sorted(df_rain["SM"].unique(), key=natural_key)
        df_rain["SM"] = pd.Categorical(df_rain["SM"], categories=order, ordered=True)
        sns.violinplot(
            data=df_rain, x="QL", y="SM",
            order=order,
            color=color,
            inner=None, cut=0, linewidth=1, width=0.75, ax=ax,
            saturation=1
        )
        yticks = ax.get_yticks()
        for i, poly in enumerate(ax.collections[:len(order)]):
            verts = poly.get_paths()[0].vertices
            x_min, x_max = verts[:, 0].min(), verts[:, 0].max()
            y_min, y_max = verts[:, 1].min(), verts[:, 1].max()
            y_center = yticks[i]
            clip_rect = plt.Rectangle(
                (x_min, y_min),
                x_max - x_min,
                y_center - y_min,
                transform=ax.transData
            )
            poly.set_clip_path(clip_rect)
        sns.boxplot(
            data=df_rain, x="QL", y="SM",
            order=order,
            width=0.2, fliersize=0,
            boxprops=dict(facecolor="white", edgecolor="black", linewidth=1),
            medianprops=dict(color="black", linewidth=1),
            whiskerprops=dict(color="black", linewidth=1),
            capprops=dict(color="black", linewidth=1),
            ax=ax
        )
        before = len(ax.collections)
        sns.stripplot(
            data=df_rain, x="QL", y="SM",
            order=order,
            color=color,
            jitter=0.05,
            size=1.5,
            alpha=0.3,
            ax=ax
        )
        after = ax.collections[before:]
        for coll in after:
            offsets = coll.get_offsets()
            offsets[:, 1] += 0.25
            coll.set_offsets(offsets)
        xmax = df_rain["QL"].max()
        for i, sm in enumerate(order):
            med = medians[sm]
            if np.isnan(med):
                continue
            ax.text(
                xmax * 0.98,
                i + 0.25,
                f"{int(med)} bp",
                va="center",
                ha="right",
                fontsize=9,
                color=color,
                bbox=dict(
                    facecolor="white",
                    edgecolor="none",
                    alpha=0.90,
                    boxstyle="round,pad=0.2"
                )
            )
        ax.set_xlim(0, xmax * 1.05)
        ax.set_yticks(range(len(order)))
        ax.set_yticklabels(order)
        ax.tick_params(axis="x", labelsize=10)
        ax.tick_params(axis="y", labelsize=11)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_title(title, color=color, fontsize=12)
        ax.grid(axis="x")
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color("black")
            spine.set_linewidth(1.5)
    sns.set_theme(style="white")
    fig, (ax1, ax2) = plt.subplots(
        1, 2,
        figsize=(7, max(5, min(14, len(samples)*0.75))),
        sharey=True
    )
    fig.patch.set_alpha(0)
    rain_panel(
        ax1, df_lf, samples,
        "Short Full-Length",
        (pl.col("TC") > 0) & (pl.col("IC") == 0) & (pl.col("FC") > 0),
        color="#EC008C"
    )
    rain_panel(
        ax2, df_lf, samples,
        "Long Full-Length",
        (pl.col("TC") > 0) & (pl.col("IC") > 0) & (pl.col("FC") > 0),
        color="#583092"
    )
    fig.supxlabel("Mean Molecule Length (bp)",  x=0.55, fontsize=11)
    plt.tight_layout()
    fig_molecule_length_sns = fig

    ## PDF creation
    figs = [fig_counts_sns,
        fig_molecules_sns,
        fig_genes_sns,
        fig_molecules_ext_sns,
        fig_molecules_ext_full_length_sns,
        fig_molecule_length_sns
    ]      
    def get_fig_size_cm(fig, dpi=300):
        w_in, h_in = fig.get_size_inches()
        return w_in * 2.54, h_in * 2.54
    frame_width_cm = 469.88976377952764 / cm
    frame_height_cm = 688.1574803149607 / cm
    scales = []
    for fig in figs:
        w_cm, h_cm = get_fig_size_cm(fig)
        scale_w = frame_width_cm / w_cm
        scale_h = frame_height_cm / h_cm
        scales.append(min(scale_w, scale_h))
    global_scale = min(scales) * 0.95
    def fig_to_rl_image_nodistort(fig, scale=global_scale, dpi=300, max_height_cm=17):
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight", pad_inches=0.2)
        buf.seek(0)
        img = Image(buf)
        px_to_cm = 2.54 / dpi
        draw_w = img.imageWidth * px_to_cm * cm * scale
        draw_h = img.imageHeight * px_to_cm * cm * scale
        max_h = max_height_cm * cm
        if draw_h > max_h:
            factor = max_h / draw_h
            draw_h *= factor
            draw_w *= factor
        img.drawWidth = draw_w
        img.drawHeight = draw_h
        img.hAlign = "CENTER"
        return img

    BACKGROUND_PDF   = background_pdf_path
    READTYPE_SCHEMA  = readtype_schema_path
    MONA_SANS_REGULAR = args.font_regular
    MONA_SANS_BOLD    = args.font_bold

    if Path(MONA_SANS_REGULAR).exists() and Path(MONA_SANS_BOLD).exists():
        pdfmetrics.registerFont(TTFont("MonaSans",      MONA_SANS_REGULAR))
        pdfmetrics.registerFont(TTFont("MonaSans-Bold", MONA_SANS_BOLD))
        font_regular = "MonaSans"
        font_bold    = "MonaSans-Bold"
    else:
        font_regular = "Helvetica"
        font_bold    = "Helvetica-Bold"

    if READTYPE_SCHEMA:
        readtype_schema = Image(str(READTYPE_SCHEMA))
        readtype_schema.drawWidth  = 12 * cm
        readtype_schema.drawHeight = 6  * cm
        readtype_schema.hAlign     = "CENTER"
        appendix_items = [(
            "RNA BaseCode read types are defined as follows:",
            [Spacer(1, 12), readtype_schema, Spacer(1, 40)]
        )]
    else:
        appendix_items = []

    doc = SimpleDocTemplate(str(out_pdf), pagesize=A4, rightMargin=2*cm, leftMargin=2*cm, topMargin=2.5*cm, bottomMargin=2.5*cm)

    def draw_background(canvas, doc):
        canvas.saveState()
        if BACKGROUND_PDF:
            canvas.drawImage(BACKGROUND_PDF, 0, 0, width=A4[0], height=A4[1], mask="auto")
        canvas.setFont(font_regular, 10)
        canvas.setFillColor(colors.black)
        canvas.drawString(doc.leftMargin + 7.5, A4[1] - doc.topMargin + 16, f"Project: {run_id}")
        canvas.setFillColor(colors.white)
        canvas.drawRightString(A4[0] - doc.rightMargin, doc.bottomMargin - 40, f"Page {canvas.getPageNumber()}") 
        canvas.restoreState()

    styles = getSampleStyleSheet()
    story  = []
    styles.add(ParagraphStyle("SectionTitle",  fontSize=14, leading=18, spaceBefore=12, spaceAfter=12, alignment=TA_CENTER, fontName=font_bold,    textColor="#583092"))
    styles.add(ParagraphStyle("SubTitle",      fontSize=12, leading=16, spaceAfter=8,   alignment=TA_LEFT,   fontName=font_bold,    textColor="#EC008C"))
    styles.add(ParagraphStyle("SubSubTitle",   fontSize=11, leading=14, spaceAfter=6,   alignment=TA_LEFT,   fontName=font_bold))
    styles.add(ParagraphStyle("Body",          fontSize=10, leading=12, spaceAfter=6,   alignment=TA_LEFT,   fontName=font_regular))
    styles.add(ParagraphStyle("Link",          fontSize=10, leading=12, spaceAfter=6,   alignment=TA_LEFT,   fontName=font_regular, textColor="#583092"))

    def section_page(title=None, subtitle=None, subsubtitle=None, items=None, page_break=True):
        story.extend(
            ([Paragraph(title,       styles["SectionTitle"])] if title       else []) +
            ([Paragraph(subtitle,    styles["SubTitle"])]     if subtitle    else []) +
            ([Paragraph(subsubtitle, styles["SubSubTitle"])]  if subsubtitle else [])
        )
        if items:
            for body, figs in items:
                if body:
                    story.append(Paragraph(body, styles["Body"]))
                if figs:
                    story.extend(figs)
        if page_break:
            story.append(PageBreak())

    table_data = [tab_samples.columns.tolist()] + tab_samples.values.tolist()
    tbl = Table(table_data, repeatRows=1)
    tbl.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1,  0), colors.HexColor("#583092")),
        ("BACKGROUND", (0, 1), (-1, -1), colors.white),
        ("TEXTCOLOR",  (0, 0), (-1,  0), colors.white),
        ("GRID",       (0, 0), (-1, -1), 0.25, colors.black),
        ("FONT",       (0, 0), (-1,  0), font_bold),
        ("ALIGN",      (0, 0), (-1,  0), "LEFT"),
        ("ALIGN",      (1, 1), (-1, -1), "LEFT"),
    ]))

    section_page(
        title="BaseCode Processing Pipeline Run Report",
        subtitle="Experiment overview",
        items=[
            (
                "Library preparation was performed according to the BaseCode protocol "
                "and libraries were sequenced on an MGI DNBSEQ-G400 platform.",
                [Spacer(1, 12), tbl, Spacer(1, 24)],
            ),
            (
                'For full documentation and analysis guides visit: '
                '<a href="https://basecode.readthedocs.io/en/latest/" color="#583092">'
                '<u>basecode.readthedocs.io</u></a>',
                [Spacer(1, 12)]
            )
        ],
        page_break=True
    )
    section_page(
        subtitle="Data QC",
        subsubtitle="Read Counts",
        items=[
            (   f"Samples were sequenced at {counts_min_M}–{counts_max_M}M total reads.",
                [fig_to_rl_image_nodistort(fig_counts_sns.figure), Spacer(1, 12)]
            ),
            (
                "<font size='8'><i>Refer to Appendix I for BaseCode read types.</i></font>",
                [Spacer(1, 12)]
            )
        ],
        page_break=True
    )
    section_page(
        subsubtitle="Transcript Coverage",
        items=[
            (
                f"Short reads were reconstructed into synthetic long reads, i.e., reconstructed molecules. On average, approximately {molecules_5_mean_M}M reconstructed molecules covered the 5' "
                f"transcript region, {molecules_int_mean_M}M the internal transcript region, and {molecules_3_mean_M}M the 3' transcript region per sample.",
                [fig_to_rl_image_nodistort(fig_molecules_sns.figure), Spacer(1, 12)],
            ),
            (
                "<font size='8'><i>Note: Reconstructed molecules covering only internal region may originate from the same molecule.</i></font>",
                [Spacer(1, 12)]
            )
        ],
        page_break=True
    )

    section_page(
        subsubtitle="Detected Genes",
        items=[
            (
                f"On average, reconstructed molecules were detected across ~{genes_detected_mean_k}k genes per sample.",
                [fig_to_rl_image_nodistort(fig_genes_sns.figure), Spacer(1, 12)],
            )
        ],
        page_break=True
    )
    section_page(
        subtitle="Reconstruction",
        subsubtitle="Detailed Transcript Coverage",
        items=[
            (
                f"Reconstructed molecules covered different transcript regions, with molecules covering 5' — 3' and 5' — Internal — 3' regions representing full-length molecules.",
                [fig_to_rl_image_nodistort(fig_molecules_ext_sns.figure), Spacer(1, 12)],
            ),
            (
                "<font size='8'><i>Note: Reconstructed molecules covering only internal region may originate from the same molecule.</i></font>",
                [Spacer(1, 12)]
            )
        ],
        page_break=True
    )
    section_page(
        subsubtitle="Full-Length Molecules",
        items=[
            (
                f"Of all reconstructed molecules containing a 3' region, on average, approximately {molecules_5_noint_3_mean}% were classified as <font color='#EC008C'><b>short full-length molecules</b></font> "
                f"(spanning 5' — 3' regions without internal part), whereas approximately {molecules_5_int_3_mean}% were classified as <font color='#583092'><b>long full-length molecules</b></font> "
                f"(spanning 5' — Internal — 3' regions).",
                [fig_to_rl_image_nodistort(fig_molecules_ext_full_length_sns.figure), Spacer(1, 12)],
            )
        ],
        page_break=True
    )
    section_page(
        subsubtitle="Full-Length Molecule Length",
        items=[
            (
                f"Among full-length reconstructed molecules, <font color='#EC008C'><b>short full-length molecules</b></font> had on average a median length of approximately {molecules_5_noint_3_length} bp per sample, "
                f"whereas <font color='#583092'><b>long full-length molecules</b></font> had on average a median length of approximately {molecules_5_int_3_length} bp per sample.",
                [fig_to_rl_image_nodistort(fig_molecule_length_sns.figure), Spacer(1, 12)],
            )
        ],
        page_break=True
    )
    section_page(
        subtitle="Accompanying Data",
        subsubtitle="Processed data",
        items=[
            (
                "Each entry in a BAM file corresponds to a reconstructed molecule.<br/>"
                "<br></br><font face='Courier'>…stitched.molecules.sorted.bam</font><br/>"
                "<font face='Courier'>…stitched.molecules.sorted.bam.bai</font>",
                [Spacer(1, 12)]
            )
        ],
        page_break=False
    )
    section_page(
        subsubtitle="Statistics",
        items=[
            (
                "Additional files to reproduce the above plots, the long form statistics file contain information extracted from the stitched bam file. Number of reads used and the length of each molecule.<br/>"
                "<br></br><font face='Courier'>…long_form_reconstruction_stats.csv</font>",
                [Spacer(1, 12)]
            )
        ],
        page_break=False
    )
    section_page(
        subsubtitle="Other",
        items=[
            (
                "<font face='Courier'>…md5sums.txt</font>",
                [Spacer(1, 12)]
            )
        ],
        page_break=False
    )
    section_page(
        subtitle="Custom tag description",
        items=[
            (
                "The processed BAM file contains the following tags:<br/>"
                "SM: Sample name.<br/>"
                "BC: Cell barcode.<br/>"
                "XT: Gene name or ID.<br/>"
                "RM: UMI based on pattern.<br/>"
                "NR: Number of reads used to stitch.<br/>"
                "ER: Number of reads covering an exon.<br/>"
                "IR: Number of reads covering an intron.<br/>"
                "FC: Number of 5' reads.<br/>"
                "IC: Number of internal reads.<br/>"
                "TC: Number of 3' reads.<br/>",
                [Spacer(1, 12)]
            )
        ],
        page_break=True
    )
    if appendix_items:
        section_page(
            subtitle="Appendices",
            subsubtitle="Appendix I",
            items=appendix_items,
            page_break=False
        )

    doc.build(
        story,
        onFirstPage=draw_background,
        onLaterPages=draw_background,
    )

## Argparse
def parse_args():
    parser = argparse.ArgumentParser(description="Generate BaseCode Processing Pipeline run report")
    parser.add_argument("--run-id",
                        help="Run ID")
    parser.add_argument("--samplesheet", required=True,
                        help="Path to samplesheet CSV")
    parser.add_argument("--summary", required=True, 
                        help="Path to summary stats CSV")
    parser.add_argument("--long-form-summary", required=True, 
                        help="Path to long-format summary stats CSV")
    parser.add_argument("--outdir", default="./results", help="Output directory")
    parser.add_argument("--verbose", action="store_true", 
                        help="Verbose logging")
    parser.add_argument("--background-pdf", default=None,
                        help="Path to background image used in the report. "
                             "If not provided, no background will be used.")
    parser.add_argument("--readtype-schema", default=None,
                        help="Path to read type schema image. "
                             "If not provided, Appendix I will be skipped.")
    parser.add_argument("--font-regular", default=None,
                        help="Path to regular TTF font."
                             "If not provided, falls back to Helvetica.")
    parser.add_argument("--font-bold", default=None,
                        help="Path to bold TTF font. "
                             "If not provided, falls back to Helvetica-Bold.")
    return parser.parse_args()
if __name__ == "__main__":
    args = parse_args()
    main(args)
    