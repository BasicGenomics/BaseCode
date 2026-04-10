rule count_lengths:
    input: bam = "results/{name}.stitched.molecules.sorted.bam".format(name=config["name"]),
           bai = "results/{name}.stitched.molecules.sorted.bam.bai".format(name=config["name"])
    output: long_form_reconstruction_stats = "results/QC_files/{name}_long_form_reconstruction_stats.csv".format(name=config["name"]),
            done = "results/dones/{name}_count_lengths.done".format(name=config["name"])
    log: "results/logs/count_lengths.log"
    shell: """
    python3 workflow/scripts/reconstruction_lengths.py -i {input.bam} -o {output.long_form_reconstruction_stats} > {log} 2>&1
    touch {output.done}
    """

rule overlap_and_mi:
    input: bam = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.reconstructed.sorted.bam".format(name=config["name"]), 
           bai = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.reconstructed.sorted.bam.bai".format(name=config["name"])
    output: overlap_and_mi_genes = "results/QC_files/{name}_overlap_and_mi_genes.csv".format(name=config["name"]), 
            overlap_and_mi_cells = "results/QC_files/{name}_overlap_and_mi_cells.csv".format(name=config["name"]),
            done = "results/dones/{name}_overlap_and_mi.done".format(name=config["name"])
    log: "results/logs/overlap_and_mi.log"
    params: gff = "{}.gff3".format(GFF)
    threads: config["threads"]
    shell: """
    python3 workflow/scripts/overlap_and_mi.py -i {input.bam} -g {params.gff} --genes-out {output.overlap_and_mi_genes} --cells-out {output.overlap_and_mi_cells} -t {threads} --gene-identifier {config[gff_gene_identifier]} > {log} 2>&1
    touch {output.done}
    """

rule status_stats:
    input: bam = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.reconstructed.sorted.bam".format(name=config["name"]), 
           bai = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.reconstructed.sorted.bam.bai".format(name=config["name"])
    output: counts_per_gene = "results/QC_files/{name}_counts_per_gene.csv".format(name=config["name"]),
            status_sum = "results/QC_files/{name}_status_sum.csv".format(name=config["name"]),
            status_fraction_per_gene = "results/QC_files/{name}_status_fraction_per_gene.csv".format(name=config["name"]),
            done = "results/dones/{name}_status_stats.done".format(name=config["name"])
    log: "results/logs/status_stats.log"
    params: gff = "{}.gff3".format(GFF)
    threads: config["threads"]
    shell: """
    python3 workflow/scripts/count_status_per_gene.py -i {input.bam} -g {params.gff} --counts {output.counts_per_gene} --sum {output.status_sum} --fraction {output.status_fraction_per_gene} -t {threads} --gene-identifier {config[gff_gene_identifier]} > {log} 2>&1
    touch {output.done}
    """
