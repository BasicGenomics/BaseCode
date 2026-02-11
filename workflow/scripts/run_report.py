import pandas as pd
import numpy as np
import polars as pl
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.colors import LinearSegmentedColormap, to_hex
import seaborn as sns
import io
import argparse
import math
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
    background_pdf_path = Path(args.background_pdf)
    if not background_pdf_path.exists():
        raise FileNotFoundError(f"Background image not found: {background_pdf_path}")
    outdir.mkdir(parents=True, exist_ok=True)
    run_id = samplesheet_path.stem.split("_")[0]
    out_pdf = outdir / f"{run_id}_run_report.pdf"
    if args.verbose:
        print(f"Run ID: {run_id}")
        print(f"Output dir: {outdir}")


    ## Metadata
    meta = pd.read_csv(samplesheet_path)
    tab_samples = meta[['SAMPLE_ID', 'DESC']]
    tab_samples = tab_samples.rename(columns={
        "SAMPLE_ID": "Sample",
        "DESC": "Description"
    })
    samples = tab_samples["Sample"].tolist()
    n = len(samples)


    ## Palette
    cmap = LinearSegmentedColormap.from_list(
        "sample_grad",
        ["#EC008C", "#97B3D6","#583092"]
    )
    colors_samples = [to_hex(cmap(i/(n-1))) for i in range(n)]
    sample_color = dict(zip(samples, colors_samples))


    ## Summary stats
    summary = pd.read_csv(summary_path)
    summary["index"] = summary["index"].astype(str)
    summary = summary.sort_values(by="index")
    df = summary.copy()
    wid_1 = len(df) * 80
    wid_2 = len(df) * 100


    ### Plot #1
    def millions(x, pos):
        return f"{x/1e6:.0f}M"
    colors_counts = {
        "3' Counts": "#583092",
        "Internal Counts": "#97B3D6",
        "5' Counts": "#EC008C"
    }
    samples = df["index"]
    c3 = df["3' Counts"]
    cI = df["Internal Counts"]
    c5 = df["5' Counts"]
    fig, ax = plt.subplots(figsize=(wid_2/165, 5))
    fig.patch.set_alpha(0)
    #fig.subplots_adjust(right=0.75)
    ax.bar(
        samples,
        c3,
        label="3' Counts",
        color=colors_counts["3' Counts"]
    )
    ax.bar(
        samples,
        cI,
        bottom=c3,
        label="Internal Counts",
        color=colors_counts["Internal Counts"]
    )
    ax.bar(
        samples,
        c5,
        bottom=c3 + cI,
        label="5' Counts",
        color=colors_counts["5' Counts"]
    )
    ax.set_xlabel("")
    ax.set_ylabel("Read Counts")
    ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False)
    plt.xticks(rotation=90, ha="center")
    #ax.ticklabel_format(style="plain", axis="y", scilimits=(0,0))
    ax.yaxis.set_major_formatter(FuncFormatter(millions))
    plt.grid(axis='x')
    for i, (v3, vI, v5) in enumerate(zip(c3, cI, c5)):
        if v3 > 0:
            ax.text(
                i, v3/2,
                f"{v3/1e6:.1f}M",
                ha="center", va="center",
                color="white",
                fontsize=9
            )
        if vI > 0:
            ax.text(
                i, v3 + vI/2,
                f"{vI/1e6:.1f}M",
                ha="center", va="center",
                color="black",
                fontsize=9
            )
        if v5 > 0:
            ax.text(
                i, v3 + vI + v5/2,
                f"{v5/1e6:.1f}M",
                ha="center", va="center",
                color="white",
                fontsize=9
            )
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.5)
    plt.show()
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
    categories = [
        "Molecules Detected (5')",
        "Molecules Detected (INT)",
        "Molecules Detected (3')"
    ]
    label_map = {
        "Molecules Detected (5')": "5'",
        "Molecules Detected (INT)": "Internal",
        "Molecules Detected (3')": "3'"
    }
    plot_df = (
        df[["index"] + categories]
        .melt(
            id_vars="index",
            value_vars=categories,
            var_name="compartment",
            value_name="molecules"
        )
    )
    plot_df["compartment"] = plot_df["compartment"].map(label_map)
    palette = sample_color 
    fig, ax = plt.subplots(figsize=(wid_2/55, 5))
    fig.patch.set_alpha(0)
    #fig.subplots_adjust(right=0.75)
    sns.barplot(
        data=plot_df,
        x="compartment",
        y="molecules",
        hue="index",
        palette=palette,
        ax=ax
    )
    ax.set_xlabel("Transcript Region")
    ax.set_ylabel("Molecules Detected")
    ax.legend(title="Sample", bbox_to_anchor=(1, 1), loc="upper left", frameon=False)
    #ax.ticklabel_format(style="plain", axis="y", scilimits=(0,0))
    ax.yaxis.set_major_formatter(FuncFormatter(millions))
    for container in ax.containers:
        for bar in container:
            h = bar.get_height()
            if h <= 0:
                continue
            ax.text(
                bar.get_x() + bar.get_width()/2,
                h+50000,
                f"{h/1e6:.1f}M",
                ha="center",
                va="bottom",
                rotation=90,
                fontsize=9,
                color=bar.get_facecolor()
            )
    ax.margins(y=0.15)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.5)
    plt.show()

    fig_molecules_sns = fig


    ### Plot #2 stats
    molecules_5_mean = df["Molecules Detected (5')"].mean()
    molecules_int_mean = df["Molecules Detected (INT)"].mean()
    molecules_3_mean = df["Molecules Detected (3')"].mean()
    molecules_5_mean_M = round(molecules_5_mean / 1e6, 1)
    molecules_int_mean_M = round(molecules_int_mean / 1e6, 1)
    molecules_3_mean_M = round(molecules_3_mean / 1e6, 1)


    ### Plot #3
    def kilos(x, pos):
        return f"{x/1e3:.0f}k"
    palette = sample_color 
    fig, ax = plt.subplots(figsize=(wid_2/90, 5))
    fig.patch.set_alpha(0)
    fig.subplots_adjust(right=0.75)
    sns.barplot(
        x=df["index"],
        y=df["Genes Detected"],
        hue=df["index"],
        palette=palette,
        ax=ax
    )
    ax.set_xlabel("")
    ax.set_ylabel("Genes Detected")
    plt.xticks(rotation=90, ha="center")
    #ax.ticklabel_format(style="plain", axis="y", scilimits=(0,0))
    ax.yaxis.set_major_formatter(FuncFormatter(kilos))
    for container in ax.containers:
        ax.bar_label(
            container,
            labels=[f"{v/1e3:.1f}k" for v in container.datavalues],
            padding=0,
            fontsize=10,
            color="black"
        )
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.5)
    plt.axhline(df["Genes Detected"].mean(), linestyle="--", color="grey", zorder=-1)
    plt.show()
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
    molecules_5_int_3_length = round(math.ceil(df_out["Median Length Completed Molecules (5' and 3' with INT) (bp)"].mean()) / 50) * 50
    molecules_5_noint_3_length = round(math.ceil(df_out["Median Length Completed Molecules (5' and 3' without INT) (bp)"].mean()) / 50) * 50


    ### Plot #4
    def millions(x, pos):
        v = x / 1e6
        if abs(v - round(v)) < 1e-9:
            return f"{int(round(v))}M"
        return f"{v:.1f}M"
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
    plot_df = (df_out[["index"] + categories].melt(id_vars="index", value_vars=categories, var_name="compartment", value_name="molecules"))
    plot_df["compartment"] = plot_df["compartment"].map(label_map)
    plot_df["molecules"] = pd.to_numeric(plot_df["molecules"], errors="coerce").fillna(0)
    order = [label_map[c] for c in categories]
    plot_df["compartment"] = pd.Categorical(plot_df["compartment"], categories=order, ordered=True)
    plot_df["pct"] = 100 * plot_df["molecules"] / plot_df.groupby("index")["molecules"].transform("sum")
    samples = list(df_out["index"]) if "index" in df_out.columns else list(plot_df["index"].unique())
    fig, axes = plt.subplots(
        nrows=1, ncols=7,
        figsize=(10, 0.6 * len(samples)),
        sharey=True, sharex=True
    )
    fig.patch.set_alpha(0)
    fig.subplots_adjust(wspace=0.45) 
    for ax, cat in zip(axes, order):
        d = (plot_df[plot_df["compartment"] == cat]
            .set_index("index")
            .reindex(samples)
            .fillna({"molecules": 0, "pct": 0}))
        bar_colors = [sample_color.get(s) for s in samples] 
        bars = ax.barh(samples, d["molecules"].to_numpy(), color=bar_colors)
        ax.set_title(cat, fontsize=9)
        ax.xaxis.set_major_formatter(FuncFormatter(millions))
        ax.tick_params(axis="x", labelsize=9)
        ax.tick_params(axis="y", labelsize=9)
        for b, m, p in zip(bars, d["molecules"].to_numpy(), d["pct"].to_numpy()):
            if m <= 0:
                continue
            ax.text(
                b.get_width() + 0.02,
                b.get_y() + b.get_height() / 2,
                f"{m/1e6:.2f}M\n({p:.1f}%)",
                ha="left", va="center",
                fontsize=7
            )
        for spine in ax.spines.values():
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(True)
            ax.spines["bottom"].set_visible(True)
            ax.spines["left"].set_color("black");   ax.spines["left"].set_linewidth(1.5)
            ax.spines["bottom"].set_color("black"); ax.spines["bottom"].set_linewidth(1.5)
    for ax in axes:
        ax.xaxis.grid(True)
        ax.yaxis.grid(False)
    axes[0].invert_yaxis()
    for ax in axes:
        ax.set_xlabel("")
    mid = len(axes) // 2
    axes[mid].set_xlabel("Molecules Detected", fontsize=9)
    plt.show()
    fig_molecules_ext_sns = fig


    ### Plot #5
    categories = [
    "Molecules Detected (3' only)",
    "Molecules Detected (3' and INT only)",
    "Molecules Detected (5' and 3' without INT)",
    "Molecules Detected (5' and 3' with INT)"    
    ]

    label_map = {
        "Molecules Detected (3' only)": "— 3'",
        "Molecules Detected (3' and INT only)": "Internal — 3'",
        "Molecules Detected (5' and 3' without INT)": "5' — 3'",
        "Molecules Detected (5' and 3' with INT)": "5' — Internal — 3'",  
    }
    order = list(label_map.values())
    plot_df = (
        df_out[["index"] + categories]
        .melt(
            id_vars="index",
            value_vars=categories,
            var_name="compartment",
            value_name="molecules"
        )
    )
    plot_df["compartment"] = plot_df["compartment"].map(label_map)
    dfp = plot_df.copy()
    dfp["compartment"] = pd.Categorical(dfp["compartment"], categories=order, ordered=True)
    sample_tot = dfp.groupby("index")["molecules"].transform("sum")
    dfp["pct"] = 100 * dfp["molecules"] / sample_tot
    fig, ax = plt.subplots(figsize=(0.9*dfp["index"].nunique(), 4))
    fig.patch.set_alpha(0)
    sns.scatterplot(
        data=dfp,
        x="index", y="compartment",
        hue="index", palette=sample_color,
        size="pct", sizes=(10, 400),  
        linewidth=0, ax=ax, legend=False
    )
    for _, r in dfp.iterrows():
        ax.annotate(f"{r.pct:.1f}%", (r["index"], r["compartment"]),
                    xytext=(0, 12), textcoords="offset points",
                    ha="center", va="bottom", fontsize=9, color="black")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=90)
    ax.invert_yaxis()
    def y_frac_for_label(label):
        labels = [str(x) for x in order]
        i = labels.index(label)
        return (i + 0.5) / len(labels)
    ax.text(
        1.02, y_frac_for_label("5' — Internal — 3'"),
        "Long Full-Length",
        transform=ax.transAxes,
        ha="left", va="center",
        fontsize=10, color="#583092",
        clip_on=False
    )
    ax.text(
        1.02, y_frac_for_label("5' — 3'"),
        "Short Full-Length",
        transform=ax.transAxes,
        ha="left", va="center",
        fontsize=10, color="#583092",
        clip_on=False
    )
    fig.subplots_adjust(right=0.82)
    ax.margins(y=0.2)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.5)
        ax.tick_params(axis="x", labelsize=9)
        ax.tick_params(axis="y", labelsize=9)
    plt.grid(False)
    plt.show()
    fig_molecules_ext_full_length_sns = fig


    ### Plot #6
    def rain_panel(ax, df_lf, samples, sample_color, title, completed_filter):
        records = []
        for sm in samples:
            df_sample = df_lf.filter(pl.col("SM") == sm)
            df_completed = df_sample.filter(completed_filter)
            gene_means = (
                df_completed.group_by("XT")
                .mean()
                .select(["QL", "BC"])
                .to_pandas()
            )
            gene_means["SM"] = sm
            records.append(gene_means)
        df_rain = pd.concat(records, ignore_index=True)
        sns.violinplot(
            data=df_rain, x="QL", y="SM",
            hue="SM", order=samples, palette=sample_color,
            inner=None, cut=0, linewidth=1, width=0.75, ax=ax
        )
        yticks = ax.get_yticks()
        for i, poly in enumerate(ax.collections[:len(samples)]):
            verts = poly.get_paths()[0].vertices
            x_min, x_max = verts[:, 0].min(), verts[:, 0].max()
            y_min, y_max = verts[:, 1].min(), verts[:, 1].max()
            y_center = yticks[i]
            clip_rect = plt.Rectangle((x_min, y_min), x_max - x_min, y_center - y_min, transform=ax.transData)
            poly.set_clip_path(clip_rect)
        sns.boxplot(
            data=df_rain, x="QL", y="SM",
            hue="SM", order=samples, width=0.2, fliersize=0,
            boxprops=dict(facecolor="white", edgecolor="black", linewidth=1),
            medianprops=dict(color="black", linewidth=1),
            whiskerprops=dict(color="black", linewidth=1),
            capprops=dict(color="black", linewidth=1),
            ax=ax
        )
        before = len(ax.collections)
        sns.stripplot(
            data=df_rain, x="QL", y="SM",
            hue="SM", order=samples, palette=sample_color,
            jitter=0.05, size=1, alpha=0.25, ax=ax
        )
        after = ax.collections[before:]
        for coll in after:
            offsets = coll.get_offsets()
            offsets[:, 1] -= -0.3
            coll.set_offsets(offsets)
        if ax.legend_:
            ax.legend_.remove()
        ax.set_xlabel("Mean Molecule Length Per Gene (bp)")
        ax.set_ylabel("")
        ax.set_title(title)
        ax.grid(axis="x")
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color("black")
            spine.set_linewidth(1.5)
    sns.set_theme(style="white")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, wid_1/100), sharey=True)
    fig.patch.set_alpha(0)
    ax1.set_facecolor("white")
    ax2.set_facecolor("white")
    rain_panel(
        ax1, df_lf, samples, sample_color,
        "Long Full-Length Molecules",
        (pl.col("TC") > 0) & (pl.col("IC") > 0) & (pl.col("FC") > 0)
    )
    rain_panel(
        ax2, df_lf, samples, sample_color,
        "Short Full-Length Molecules",
        (pl.col("TC") > 0) & (pl.col("IC") == 0) & (pl.col("FC") > 0)
    )
    plt.show()
    fig_molecule_length_sns = fig




    ## PDF creation
    def fig_to_rl_image_nodistort(fig, max_width_cm=16, max_height_cm=9, dpi=300):
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight", pad_inches=0.2)
        buf.seek(0)
        img = Image(buf)
        aspect = img.imageWidth / img.imageHeight
        if aspect < max_width_cm / max_height_cm:
            img.drawHeight = max_height_cm * cm
            img.drawWidth = img.drawHeight * aspect
        else:
            img.drawWidth = max_width_cm * cm
            img.drawHeight = img.drawWidth / aspect
        img.hAlign = "CENTER"
        return img
    BACKGROUND_PDF = background_pdf_path 
    doc = SimpleDocTemplate(str(out_pdf), pagesize=A4, rightMargin=2*cm, leftMargin=2*cm, topMargin=2.5*cm, bottomMargin=2.5*cm)
    def draw_background(canvas, doc):
        canvas.saveState()
        canvas.drawImage(BACKGROUND_PDF, 0, 0, width=A4[0], height=A4[1], mask="auto")
        canvas.setFont("Helvetica", 10)
        canvas.setFillColor(colors.black)
        canvas.drawString(doc.leftMargin + 7.5, A4[1] - doc.topMargin + 16, f"Project: {run_id}")
        canvas.setFillColor(colors.white)
        canvas.drawRightString(A4[0] - doc.rightMargin, doc.bottomMargin - 40, f"Page {canvas.getPageNumber()}") 
        canvas.restoreState()
    styles = getSampleStyleSheet()
    story = []
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle("SectionTitle",  fontSize=16, leading=18, spaceBefore=12, spaceAfter=12, alignment=TA_CENTER, fontName="Helvetica-Bold", textColor="#583092"))
    styles.add(ParagraphStyle("SubTitle",      fontSize=14, leading=16, spaceAfter=8,  alignment=TA_LEFT, fontName="Helvetica-Bold", textColor="#EC008C"))
    styles.add(ParagraphStyle("SubSubTitle",   fontSize=12, leading=14, spaceAfter=6,  alignment=TA_LEFT, fontName="Helvetica-Bold"))
    styles.add(ParagraphStyle("Body",          fontSize=11, leading=12, spaceAfter=6,  alignment=TA_LEFT))
    def section_page(title=None, subtitle=None, subsubtitle=None, items=None, page_break=True):
        story.extend(
            ([Paragraph(title, styles["SectionTitle"])] if title else []) +
            ([Paragraph(subtitle, styles["SubTitle"])] if subtitle else []) +
            ([Paragraph(subsubtitle, styles["SubSubTitle"])] if subsubtitle else [])
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
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#583092")),
        ("BACKGROUND", (0, 1), (-1, -1), colors.white),  
        ("TEXTCOLOR",  (0, 0), (-1, 0), colors.white),
        ("GRID",       (0, 0), (-1, -1), 0.25, colors.black),
        ("FONT",       (0, 0), (-1, 0), "Helvetica-Bold"),
        ("ALIGN",      (0, 0), (-1, 0), "LEFT"),
        ("ALIGN",      (1, 1), (-1, -1), "LEFT"),
    ]))
    section_page(
        title="BaseCode Data Report",
        subtitle="Experiment overview",
        items=[(
            "Library preparation was performed according to the BaseCode protocol "
            "and libraries were sequenced "
            "on an MGI DNBSEQ-G400 platform.",
            [Spacer(1, 12), tbl]
        )]
    )
    section_page(
        subtitle="Data QC",
        subsubtitle="Read Counts",
        items=[(f"Samples were sequenced at {counts_min_M}–{counts_max_M}M total reads.",
        [fig_to_rl_image_nodistort(fig_counts_sns.figure)]
    )],
        page_break=True
    )
    section_page(
        subsubtitle="Transcript Coverage",
        items=[
            (
                f"Short reads are reconstructed into synthetic long reads, i.e reconstructed molecules. On average, approximately {molecules_5_mean_M}M reconstructed molecules covered the 5' "
                f"transcript region, {molecules_int_mean_M}M the internal transcript region, "
                f"and {molecules_3_mean_M}M the 3' transcript region per sample.",
                [fig_to_rl_image_nodistort(fig_molecules_sns.figure), Spacer(1, 12)],
            ),
            (
                "<font size='9'><i>Note: Reconstructed molecules covering only internal region may originate from the same molecule.</i></font>",
                []
            ),
        ],
        page_break=True
    )
    section_page(
        subsubtitle="Detected Genes",
        items=[
            (
                f"On average, reconstructed molecules were detected across ~{genes_detected_mean_k}k genes per sample.",
                [fig_to_rl_image_nodistort(fig_genes_sns.figure, max_height_cm=7)],
            ),
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
                "<font size='9'><i>Note: Reconstructed molecules covering only internal region may originate from the same molecule.</i></font>",
                [Spacer(1, 12)]
            ),
        ],
        page_break=True
    )
    section_page(
        subsubtitle="Full-Length Molecules",
        items=[
            (
                f"Of all reconstructed molecules containing a 3' region, on average, approximately {molecules_5_int_3_mean}% were classified as <font color='#583092'><b>long full-length molecules</b></font> "
                f"(spanning 5' — Internal — 3' regions), whereas approximately {molecules_5_noint_3_mean}% were classified as <font color='#583092'><b>short full-length molecules</b></font> "
                f"(spanning 5' — 3' regions without internal part).",
                [fig_to_rl_image_nodistort(fig_molecules_ext_full_length_sns.figure)],
            )
        ],
        page_break=True
    )
    section_page(
        subsubtitle="Full-Length Molecule Length",
        items=[
            (
                f"Among full-length reconstructed molecules, <font color='#583092'><b>long full-length molecules</b></font> had a median length of approximately {molecules_5_int_3_length} bp, "
                f"whereas <font color='#583092'><b>short full-length molecules</b></font> had a median length of approximately {molecules_5_noint_3_length} bp.",
                [fig_to_rl_image_nodistort(fig_molecule_length_sns.figure)],
            )
        ],
        page_break=True
    )
    section_page(
        subtitle="Accompanying Data",
        subsubtitle="Processed data",
        items=[(
            "Each entry corresponds to a reconstructed molecule.<br/>"
            "<br></br><font face='Courier'>…stitched.molecules.sorted.bam</font><br/>"
            "<font face='Courier'>…stitched.molecules.sorted.bam.bai</font>",
            [Spacer(1, 12)]
        )],
        page_break=False
    )
    section_page(
        subsubtitle="Statistics",
        items=[(
            "Additional files to reproduce the above plots, the long form statistics file contain information extracted from the stitched bam file. Number of reads used and the length of each molecule.<br/>"
            "<br></br><font face='Courier'>…long_form_reconstruction_stats.csv</font>",
            [Spacer(1, 12)]
        )],
        page_break=False
    )
    section_page(
        subsubtitle="Other",
        items=[(
            "<font face='Courier'>…md5sums.txt</font>",
            [Spacer(1, 12)]
        )],
        page_break=False
    )
    section_page(
        subtitle="Custom tag description",
        items=[(
            "The processed data BAM file contains the following tags:<br/>"
            "NR: Number of reads used to stitch.<br/>"
            "ER: Number of reads covering an exon.<br/>"
            "IR: Number of reads covering an intron.<br/>"
            "TC: Number of 3' reads.<br/>"
            "IC: Number of internal reads.<br/>"
            "FC: Number of 5' reads.<br/>"
            "BC: Cell barcode.<br/>"
            "XT: Gene name or ID.<br/>"
            "RM: UMI based on pattern.<br/>"
            "SM: Sample name.",
            []
        )]
    )
    doc.build(
        story,
        onFirstPage=draw_background,
        onLaterPages=draw_background,
    )

## Argparse
def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate BaseCode run report"
    )
    parser.add_argument("--samplesheet", required=True, help="Path to samplesheet CSV")
    parser.add_argument("--summary", required=True, help="Path to summary stats CSV")
    parser.add_argument("--long-form-summary", required=True, help="Path to long-format summary stats CSV")
    parser.add_argument("--outdir", default="./results", help="Output directory")
    parser.add_argument("--verbose", action="store_true", help="Verbose logging")
    parser.add_argument("--background-pdf", default="workflow/resources/Background_v1.1.png", help="Path to background image used in the report")
    return parser.parse_args()
if __name__ == "__main__":
    args = parse_args()
    main(args)