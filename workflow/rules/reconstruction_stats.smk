rule count_lengths:
    input: bam = "results/{sample}.stitched.molecules.sorted.bam"
    output: "results/QC_files/{sample}_long_form_reconstruction_stats.csv"
    threads: config["threads"]
    params: gtf =  "{}.gff3".format(GTFFILE)
    conda: "../envs/full.yaml"
    log: "results/QC_files/logs/{sample}.count_lengths.log"
    shell: "python3 workflow/scripts/reconstruction_lengths.py -i {input.bam} -o {output} > {log} 2>&1"


def input_bam(wc,ext):
    
        if config['params_key']=='bulk':
                suf = ".reads.aligned_trimmed_genetagged_sorted.reconstructed.sorted"
        else:
                suf = ".reads.aligned_trimmed_genetagged_sorted_umicorrected.reconstructed.sorted"

        if ext =='bai':
                ext_ = "bam.bai"
        else:
                ext_ = "bam"

        input_bam_list = os.path.join("results/intermediate/", f"{wc.sample}{suf}.{ext_}")
        
        return input_bam_list

rule overlap_and_mi:
    input:  bam = lambda wc: input_bam(wc, "bam"),
            bai = lambda wc: input_bam(wc, "bai")
    output: genes_out = "results/QC_files/{sample}_overlap_and_mi_genes.csv",
            cells_out = "results/QC_files/{sample}_overlap_and_mi_cells.csv"
    threads: config["threads"]
    params: gtf =  "{}.gff3".format(GTFFILE)
    conda: "../envs/full.yaml"
    log: "results/QC_files/logs/{sample}.overlap_and_mi.log"
    shell: "python3 workflow/scripts/overlap_and_mi.py -i {input.bam} -g {params.gtf} --genes-out {output.genes_out} --cells-out {output.cells_out} -t {threads} --gene-identifier {config[gff_gene_identifier]} > {log} 2>&1"

rule status_stats:
    input:  bam = lambda wc: input_bam(wc, "bam"),
            bai = lambda wc: input_bam(wc, "bai")
    output: counts_per_gene = "results/QC_files/{sample}_counts_per_gene.csv",
            sum_over_genes = "results/QC_files/{sample}_status_sum.csv",
            fraction_per_gene = "results/QC_files/{sample}_status_fraction_per_gene.csv"
    threads: config["threads"]
    params: gtf =  "{}.gff3".format(GTFFILE)
    conda: "../envs/full.yaml"
    log: "results/QC_files/logs/{sample}.status_stats.log"
    shell: "python3 workflow/scripts/count_status_per_gene.py -i {input.bam} -g {params.gtf} --counts {output.counts_per_gene} --sum {output.sum_over_genes} --fraction {output.fraction_per_gene} -t {threads} --gene-identifier {config[gff_gene_identifier]} > {log} 2>&1"

