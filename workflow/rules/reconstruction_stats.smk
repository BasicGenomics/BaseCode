rule count_lengths:
    input: bam = "results/{name}.stitched.molecules.sorted.bam".format(name=config["name"])
    output: "results/QC_files/{name}_long_form_reconstruction_stats.csv".format(name=config["name"])
    threads: config["threads"]
    params: gtf =  "{}.gff3".format(GTFFILE)
    conda: "../envs/full.yaml"
    log: "results/logs/count_lengths.log"
    shell: "python3 workflow/scripts/reconstruction_lengths.py -i {input.bam} -o {output} > {log} 2>&1"

rule overlap_and_mi:
    input: bam = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.reconstructed.sorted.bam".format(name=config["name"]), bai = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.reconstructed.sorted.bam.bai".format(name=config["name"])
    output: genes_out = "results/QC_files/{name}_overlap_and_mi_genes.csv".format(name=config["name"]), cells_out = "results/QC_files/{name}_overlap_and_mi_cells.csv".format(name=config["name"])
    threads: config["threads"]
    params: gtf =  "{}.gff3".format(GTFFILE)
    conda: "../envs/full.yaml"
    log: "results/logs/overlap_and_mi.log"
    shell: "python3 workflow/scripts/overlap_and_mi.py -i {input.bam} -g {params.gtf} --genes-out {output.genes_out} --cells-out {output.cells_out} -t {threads} > {log} 2>&1"

rule status_stats:
    input: bam = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.reconstructed.sorted.bam".format(name=config["name"]), bai = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.reconstructed.sorted.bam.bai".format(name=config["name"])
    output: counts_per_gene = "results/QC_files/{name}_counts_per_gene.csv".format(name=config["name"]),
            sum_over_genes = "results/QC_files/{name}_status_sum.csv".format(name=config["name"]),
            fraction_per_gene = "results/QC_files/{name}_status_fraction_per_gene.csv".format(name=config["name"])
    threads: config["threads"]
    params: gtf =  "{}.gff3".format(GTFFILE)
    conda: "../envs/full.yaml"
    log: "results/logs/status_stats.log"
    shell: "python3 workflow/scripts/count_status_per_gene.py -i {input.bam} -g {params.gtf} --counts {output.counts_per_gene} --sum {output.sum_over_genes} --fraction {output.fraction_per_gene} -t {threads} > {log} 2>&1"

