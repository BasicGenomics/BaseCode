rule general:
    input: "results/{name}.reads.aligned_trimmed_genetagged_sorted.bam".format(name=config["name"])
    output: "results/QC_files/{name}_read_type_per_sample.csv".format(name=config["name"]), "results/QC_files/{name}_mapping_categories_per_sample.csv".format(name=config["name"]), "results/QC_files/{name}_nonbarcoded_mapping_categories_per_sample.csv".format(name=config["name"])
    params: cbcpath = "results/{name}_cell_barcodes.txt".format(name=config["name"]), smppath = "results/{name}_sample_barcodes.txt".format(name=config["name"])
    conda: "../envs/full.yaml"
    shell: "python3 workflow/scripts/general_stats.py -i {input} -o results/QC_files -s {params.smppath} -p {config[name]} -c {params.cbcpath}"

rule insertion_sizes:
    input: "results/{name}.reads.aligned_trimmed_genetagged_sorted.bam".format(name=config["name"])
    output: "results/QC_files/{name}_insert_sizes_per_sample_barcode.csv".format(name=config["name"]), "results/QC_files/{name}_double_reads_per_sample_barcode.csv".format(name=config["name"])
    params: cbcpath = "results/{name}_cell_barcodes.txt".format(name=config["name"]), smppath = "results/{name}_sample_barcodes.txt".format(name=config["name"])
    conda: "../envs/full.yaml"
    shell: "python3 workflow/scripts/insertion_sizes.py -i {input} -o results/QC_files -s {params.smppath} -p {config[name]} -c {params.cbcpath}"

rule conversion:
    input: bam = "results/{name}.reads.aligned_trimmed_genetagged_sorted.bam".format(name=config["name"]), bai = "results/{name}.reads.aligned_trimmed_genetagged_sorted.bam.bai".format(name=config["name"]), samplesheet = "results/{name}_samplesheet.csv".format(name=config["name"])
    output: "results/QC_files/{name}_conversion_rate_pos.csv".format(name=config["name"]),"results/QC_files/{name}_conversion_rate_total.csv".format(name=config["name"])
    params: fasta = REFFILE, gtf =  "{}.gtf".format(GTFFILE)
    threads: int(config["threads"]/2)
    conda: "../envs/full.yaml"
    shell: "python3 workflow/scripts/conversion_rates.py -i {input.bam} -f {params.fasta} -g {params.gtf} -o results/QC_files -s {input.samplesheet} -p {config[name]} -t {threads}"

rule stats:
    input: bam = "results/{name}.reads.aligned_trimmed_genetagged_sorted.bam".format(name=config["name"]), bai = "results/{name}.reads.aligned_trimmed_genetagged_sorted.bam.bai".format(name=config["name"]), samplesheet = "results/{name}_samplesheet.csv".format(name=config["name"])
    output: "results/QC_files/{name}_human_stats.csv".format(name=config["name"]), "results/QC_files/{name}_read_matrix_human.csv".format(name=config["name"])
    params: gtf = "{}.gtf".format(GTFFILE)
    threads: int(config["threads"]/2)
    conda: "../envs/full.yaml"
    shell: "python3 workflow/scripts/human_stats.py -i {input.bam}  -g {params.gtf} -o results/QC_files -s {input.samplesheet} -p {config[name]} -t {threads}"

rule conversion_binomial_mixture:
    input: bam = "results/{name}.reads.aligned_trimmed_genetagged_sorted.bam".format(name=config["name"]), bai = "results/{name}.reads.aligned_trimmed_genetagged_sorted.bam.bai".format(name=config["name"])
    output: "results/QC_files/{name}_conversion_binomial_mixture.csv".format(name=config["name"])
    params: fasta = REFFILE, gtf =  "{}.gtf".format(GTFFILE)
    threads: int(config["threads"]/2)
    conda: "../envs/full.yaml"
    shell: "python3 workflow/scripts/conversion_binomial_mixture.py -i {input.bam} -f {params.fasta} -g {params.gtf} -o results/QC_files -p {config[name]} -t {threads}"