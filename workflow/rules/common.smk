import os
import glob

EXCLUDE = {"Unassigned"}

if INCLUDE_unassigned:
    EXCLUDE = {}

def get_final_output():
    all_samples = get_samples()
    only = config.get("sample_i", None) 

    if only:
        sample_in = [s for s in all_samples if s == only and s not in EXCLUDE]
    else:
        sample_in = all_samples
  
    final_output = [
    *expand("results/done/{sample}.processing.done",sample=sample_in),
    *expand("results/QC_files/{sample}_read_type_per_sample.csv",sample=sample_in),
    *expand("results/QC_files/{sample}_conversion_rate_total.csv",sample=sample_in),
    *expand("results/read_flow_files/{sample}_fastq_processed_stats.json",sample=sample_in),
    *expand("results/read_flow_files/{sample}_fastq_trimmed_stats.json",sample=sample_in),
    *expand("results/read_flow_files/{sample}_fastq_tooshort_stats.json",sample=sample_in),
    *expand("results/read_flow_files/{sample}_mapping_group.csv",sample=sample_in),
    *expand("results/read_flow_files/{sample}_status_group.csv",sample=sample_in),
    *expand("results/QC_files/{sample}_long_form_reconstruction_stats.csv",sample=sample_in),
    *expand("results/QC_files/{sample}_overlap_and_mi_genes.csv",sample=sample_in),
    *expand("results/QC_files/{sample}_overlap_and_mi_cells.csv",sample=sample_in),
    *expand("results/QC_files/{sample}_counts_per_gene.csv",sample=sample_in),
    *expand("results/QC_files/{sample}_status_sum.csv",sample=sample_in),
    *expand("results/QC_files/{sample}_status_fraction_per_gene.csv",sample=sample_in),
    *expand("results/QC_files/{sample}_insert_sizes_per_sample_barcode.csv",sample=sample_in),
    *expand("results/QC_files/{sample}_double_reads_per_sample_barcode.csv",sample=sample_in),
    *expand("results/QC_files/{sample}_summary_stats.csv",sample=sample_in)
    ]

    return final_output

def get_preprocessing_output(wildcards):
    output = [
    "results/logs/python_version.log",
     *expand("{FASTQ_DIR}/{sample}_1.fq.gz",sample=[s for s in get_samples()]),
    ]
    return output
