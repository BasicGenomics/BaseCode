OUT_DIR = "results/read_flow_files"

rule fastq_processed:
    input: r1 = "results/{name}.read1.fastq.gz".format(name=config["name"]),
           r2 = "results/{name}.read2.fastq.gz".format(name=config["name"]),
           cbcpath = "results/{name}_cell_barcodes.txt".format(name=config["name"]),
           pbcpath = "results/{name}_sample_barcodes.txt".format(name=config["name"])
    output: "results/read_flow_files/{name}_fastq_processed_stats.json".format(name=config["name"])
    threads: int(config["threads"])-1
    shell: "{config[resource_dir]}/binaries/analyze_fastq --read1 {input.r1} --read2 {input.r2} --cbcpath {input.cbcpath} --pbcpath {input.pbcpath} --threads {config[threads]} --output {output}"

rule fastq_trimmed:
    input: r1 = "results/{name}.trimmed.read1.fastq.gz".format(name=config["name"]),
           r2 = "results/{name}.trimmed.read2.fastq.gz".format(name=config["name"]),
           cbcpath = "results/{name}_cell_barcodes.txt".format(name=config["name"]),
           pbcpath = "results/{name}_sample_barcodes.txt".format(name=config["name"])
    output: "results/read_flow_files/{name}_fastq_trimmed_stats.json".format(name=config["name"])
    threads: int(config["threads"])-1
    shell: "{config[resource_dir]}/binaries/analyze_fastq --read1 {input.r1} --read2 {input.r2} --cbcpath {input.cbcpath} --pbcpath {input.pbcpath} --threads {config[threads]} --output {output}"

rule fastq_too_short:
    input: r1 = "results/{name}.tooshort.read1.fastq.gz".format(name=config["name"]),
           r2 = "results/{name}.tooshort.read2.fastq.gz".format(name=config["name"]),
           cbcpath = "results/{name}_cell_barcodes.txt".format(name=config["name"]),
           pbcpath = "results/{name}_sample_barcodes.txt".format(name=config["name"])
    output: "results/read_flow_files/{name}_fastq_tooshort_stats.json".format(name=config["name"])
    threads: int(config["threads"])-1
    shell: "{config[resource_dir]}/binaries/analyze_fastq --read1 {input.r1} --read2 {input.r2} --cbcpath {input.cbcpath} --pbcpath {input.pbcpath} --threads {config[threads]} --output {output}"

rule bam_geneassigned:
    input: bam = "results/{name}.reads.aligned_trimmed_genetagged_sorted.bam".format(name=config["name"]),
           cbcpath = "results/{name}_cell_barcodes.txt".format(name=config["name"]),
           pbcpath = "results/{name}_sample_barcodes.txt".format(name=config["name"])
    output: mapping_group = "results/read_flow_files/{name}_mapping_group.csv".format(name=config["name"]),
            reassignment = "results/read_flow_files/{name}_reassignment.csv".format(name=config["name"])
    conda: "../envs/full.yaml"
    shell: "python3 workflow/scripts/read_flow_mapped.py -i {input.bam} -c {input.cbcpath} -s {input.pbcpath} --mapping-group-out {output.mapping_group}"

rule bam_reconstructed:
    input: bam = "results/{name}.reads.aligned_trimmed_genetagged_sorted.reconstructed.sorted.bam".format(name=config["name"]),
            cbcpath = "results/{name}_cell_barcodes.txt".format(name=config["name"]), 
            pbcpath = "results/{name}_sample_barcodes.txt".format(name=config["name"])
    output: status_group = "results/read_flow_files/{name}_status_group.csv".format(name=config["name"])
    conda: "../envs/full.yaml"
    shell: "python3 workflow/scripts/read_flow_reconstructed.py -i {input.bam} -c {input.cbcpath} -s {input.pbcpath} --status-group-out {output.status_group}"