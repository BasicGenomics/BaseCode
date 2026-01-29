def get_final_output(wildcards):
    final_output = [
    "results/metadata/{}_cell_barcodes.txt".format(config['name']),
    *expand("results/QC_files/{sample}_read_type_per_sample.csv",sample=get_samples()),
    *expand("results/QC_files/{sample}_conversion_rate_total.csv",sample=get_samples()),
    *expand("results/read_flow_files/{sample}_fastq_processed_stats.json",sample=get_samples()),
    *expand("results/read_flow_files/{sample}_fastq_trimmed_stats.json",sample=get_samples()),
    *expand("results/read_flow_files/{sample}_fastq_tooshort_stats.json",sample=get_samples()),
    *expand("results/read_flow_files/{sample}_mapping_group.csv",sample=get_samples()),
    *expand("results/read_flow_files/{sample}_status_group.csv",sample=get_samples()),
    *expand("results/QC_files/{sample}_long_form_reconstruction_stats.csv",sample=get_samples()),
    *expand("results/QC_files/{sample}_overlap_and_mi_genes.csv",sample=get_samples()),
    *expand("results/QC_files/{sample}_overlap_and_mi_cells.csv",sample=get_samples()),
    *expand("results/QC_files/{sample}_counts_per_gene.csv",sample=get_samples()),
    *expand("results/QC_files/{sample}_status_sum.csv",sample=get_samples()),
    *expand("results/QC_files/{sample}_status_fraction_per_gene.csv",sample=get_samples()),
    *expand("results/QC_files/{sample}_insert_sizes_per_sample_barcode.csv",sample=get_samples()),
    *expand("results/QC_files/{sample}_double_reads_per_sample_barcode.csv",sample=get_samples()),
    *expand("results/logs/python_version.log",sample=get_samples()),
    *expand("results/QC_files/{sample}_summary_stats.csv",sample=get_samples())
    ]

    return final_output


def get_final_output_sc(wildcards):
    final_output = [
    "results/metadata/{}_cell_barcodes.txt".format(config['name'])]

    return final_output
