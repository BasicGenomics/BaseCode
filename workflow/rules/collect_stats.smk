rule general_stats:
    input:  bam = "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted.bam", 
            cbcpath = "results/metadata/{name}_cell_barcodes.txt".format(name=config["name"]), 
            smppath = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"])
    output: "results/QC_files/{sample}_read_type_per_sample.csv", 
            "results/QC_files/{sample}_mapping_categories_per_sample.csv",
            "results/QC_files/{sample}_nonbarcoded_mapping_categories_per_sample.csv"
    conda: "../envs/full.yaml"
    log: "results/QC_files/logs/{sample}.general_stats.log"
    shell: "python3 workflow/scripts/general_stats.py -i {input.bam} -o results/QC_files -s {input.smppath} -p {wildcards.sample} -c {input.cbcpath} > {log} 2>&1"

rule insertion_sizes:
    input:  bam = "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted.bam", 
            cbcpath = "results/metadata/{name}_cell_barcodes.txt".format(name=config["name"]), 
            smppath = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"])
    output: "results/QC_files/{sample}_insert_sizes_per_sample_barcode.csv", 
            "results/QC_files/{sample}_double_reads_per_sample_barcode.csv"
    conda: "../envs/full.yaml"
    log: "results/QC_files/logs/{sample}.insertion_sizes.log"
    shell: "python3 workflow/scripts/insertion_sizes.py -i {input.bam} -o results/QC_files -s {input.smppath} -p {wildcards.sample} -c {input.cbcpath} > {log} 2>&1"

rule conversion:
    input:  bam = "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted.bam",
            bai = "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted.bam.bai", 
            samplesheet = "results/metadata/{name}_samplesheet.csv".format(name=config["name"])
    output: "results/QC_files/{sample}_conversion_rate_pos.csv",
            "results/QC_files/{sample}_conversion_rate_total.csv"
    params: 
        fasta = REFFILE, 
        gtf =  "{}.gff3".format(GTFFILE),
        only_sample_flag = lambda wc: bulk_flag(config)
    threads: int(config["threads"]/2)
    conda: "../envs/full.yaml"
    log: "results/QC_files/logs/{sample}.conversion.log"
    shell: "python3 workflow/scripts/conversion_rates.py -i {input.bam} -f {params.fasta} -g {params.gtf} -o results/QC_files -s {input.samplesheet} -p {wildcards.sample} -t {threads} --gene-identifier {config[gff_gene_identifier]} {params.only_sample_flag}> {log} 2>&1"

rule summary_stats:
    input:  long_form = "results/QC_files/{sample}_long_form_reconstruction_stats.csv", 
            json = "results/read_flow_files/{sample}_fastq_processed_stats.json", 
            sample_map = "results/metadata/{name}_sample_map.yaml".format(name=config["name"]), 
    output: summary_file = "results/QC_files/{sample}_summary_stats.csv"
    log: "results/QC_files/logs/{sample}.summary_stats.log"
    threads: 1
    shell: "echo BaseCode Processing Pipeline Finished &&  python3 workflow/scripts/summary_stats.py --long-form {input.long_form} --json {input.json} --sample-map {input.sample_map} --prefix {wildcards.sample} --output {output.summary_file}"


