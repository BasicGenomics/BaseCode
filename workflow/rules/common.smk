def get_final_output(wildcards):
    final_output = ["results/{name}.stitched.molecules.sorted.bam".format(name=config["name"]),
    "results/QC_files/{name}_read_type_per_sample.csv".format(name=config["name"]),
    "results/QC_files/{name}_conversion_rate_total.csv".format(name=config["name"]),
    "results/read_flow_files/{name}_fastq_processed_stats.json".format(name=config["name"]),
    "results/read_flow_files/{name}_fastq_trimmed_stats.json".format(name=config["name"]),
    "results/read_flow_files/{name}_fastq_tooshort_stats.json".format(name=config["name"]),
    "results/read_flow_files/{name}_mapping_group.csv".format(name=config["name"]),
    "results/read_flow_files/{name}_status_group.csv".format(name=config["name"]),
    "results/QC_files/{name}_long_form_reconstruction_stats.csv".format(name=config["name"]),
    "results/QC_files/{name}_overlap_and_mi_genes.csv".format(name=config["name"]),
    "results/QC_files/{name}_overlap_and_mi_cells.csv".format(name=config["name"]),
    "results/QC_files/{name}_counts_per_gene.csv".format(name=config["name"]),
    "results/QC_files/{name}_status_sum.csv".format(name=config["name"]),
    "results/QC_files/{name}_status_fraction_per_gene.csv".format(name=config["name"]),
    "results/QC_files/{name}_insert_sizes_per_sample_barcode.csv".format(name=config["name"]),
    "results/QC_files/{name}_double_reads_per_sample_barcode.csv".format(name=config["name"]),
    "results/QC_files/{name}_conversion_binomial_mixture.csv".format(name=config["name"])]

    return final_output