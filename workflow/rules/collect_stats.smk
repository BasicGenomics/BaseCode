rule general:
    input: bam = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.bam".format(name=config["name"]), cbcpath = "results/metadata/{name}_cell_barcodes.txt".format(name=config["name"]), smppath = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"])
    output: "results/QC_files/{name}_read_type_per_sample.csv".format(name=config["name"]), "results/QC_files/{name}_mapping_categories_per_sample.csv".format(name=config["name"]), "results/QC_files/{name}_nonbarcoded_mapping_categories_per_sample.csv".format(name=config["name"])
    conda: "../envs/full.yaml"
    shell: "python3 workflow/scripts/general_stats.py -i {input.bam} -o results/QC_files -s {input.smppath} -p {config[name]} -c {input.cbcpath}"

rule insertion_sizes:
    input: bam = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.bam".format(name=config["name"]), cbcpath = "results/metadata/{name}_cell_barcodes.txt".format(name=config["name"]), smppath = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"])
    output: "results/QC_files/{name}_insert_sizes_per_sample_barcode.csv".format(name=config["name"]), "results/QC_files/{name}_double_reads_per_sample_barcode.csv".format(name=config["name"])
    conda: "../envs/full.yaml"
    shell: "python3 workflow/scripts/insertion_sizes.py -i {input.bam} -o results/QC_files -s {input.smppath} -p {config[name]} -c {input.cbcpath}"

rule conversion:
    input: bam = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.bam".format(name=config["name"]), bai = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.bam.bai".format(name=config["name"]), samplesheet = "results/metadata/{name}_samplesheet.csv".format(name=config["name"])
    output: "results/QC_files/{name}_conversion_rate_pos.csv".format(name=config["name"]),"results/QC_files/{name}_conversion_rate_total.csv".format(name=config["name"])
    params: fasta = REFFILE, gtf =  "{}.gtf".format(GTFFILE)
    threads: int(config["threads"]/2)
    conda: "../envs/full.yaml"
    shell: "python3 workflow/scripts/conversion_rates.py -i {input.bam} -f {params.fasta} -g {params.gtf} -o results/QC_files -s {input.samplesheet} -p {config[name]} -t {threads}"

