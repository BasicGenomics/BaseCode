OUT_DIR = "results/read_flow_files"

rule fastq_processed:
    input: r1 = "results/intermediate/{name}.read1.fastq.gz".format(name=config["name"]),
           r2 = "results/intermediate/{name}.read2.fastq.gz".format(name=config["name"]),
           cbcpath = "results/metadata/{name}_cell_barcodes.txt".format(name=config["name"]),
           pbcpath = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"])
    output: "results/read_flow_files/{name}_fastq_processed_stats.json".format(name=config["name"])
    threads: int(config["threads"])-1
    log: "results/logs/fastq_processed.log"
    shell: "{config[resource_dir]}/binaries/analyze_fastq --read1 {input.r1} --read2 {input.r2} --cbcpath {input.cbcpath} --pbcpath {input.pbcpath} --threads {config[threads]} --output {output}  > {log} 2>&1"

rule fastq_trimmed:
    input: r1 = "results/intermediate/{name}.trimmed.read1.fastq.gz".format(name=config["name"]),
           r2 = "results/intermediate/{name}.trimmed.read2.fastq.gz".format(name=config["name"]),
           cbcpath = "results/metadata/{name}_cell_barcodes.txt".format(name=config["name"]),
           pbcpath = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"])
    output: "results/read_flow_files/{name}_fastq_trimmed_stats.json".format(name=config["name"])
    threads: int(config["threads"])-1
    log: "results/logs/fastq_trimmed.log"
    shell: "{config[resource_dir]}/binaries/analyze_fastq --read1 {input.r1} --read2 {input.r2} --cbcpath {input.cbcpath} --pbcpath {input.pbcpath} --threads {config[threads]} --output {output}  > {log} 2>&1"

rule fastq_too_short:
    input: r1 = "results/intermediate/{name}.tooshort.read1.fastq.gz".format(name=config["name"]),
           r2 = "results/intermediate/{name}.tooshort.read2.fastq.gz".format(name=config["name"]),
           cbcpath = "results/metadata/{name}_cell_barcodes.txt".format(name=config["name"]),
           pbcpath = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"])
    output: "results/read_flow_files/{name}_fastq_tooshort_stats.json".format(name=config["name"])
    threads: int(config["threads"])-1
    log: "results/logs/fastq_too_short.log"
    shell: "{config[resource_dir]}/binaries/analyze_fastq --read1 {input.r1} --read2 {input.r2} --cbcpath {input.cbcpath} --pbcpath {input.pbcpath} --threads {config[threads]} --output {output}  > {log} 2>&1"

rule bam_geneassigned:
    input: bam = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.bam".format(name=config["name"]),
           cbcpath = "results/metadata/{name}_cell_barcodes.txt".format(name=config["name"]),
           pbcpath = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"])
    output: mapping_group = "results/read_flow_files/{name}_mapping_group.csv".format(name=config["name"])
    conda: "../envs/full.yaml"
    log: "results/logs/bam_geneassigned.log"
    shell: "python3 workflow/scripts/read_flow_mapped.py -i {input.bam} -c {input.cbcpath} -s {input.pbcpath} --mapping-group-out {output.mapping_group}  > {log} 2>&1"

rule bam_reconstructed:
    input: bam = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.reconstructed.sorted.bam".format(name=config["name"]),
            cbcpath = "results/metadata/{name}_cell_barcodes.txt".format(name=config["name"]), 
            pbcpath = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"])
    output: status_group = "results/read_flow_files/{name}_status_group.csv".format(name=config["name"])
    conda: "../envs/full.yaml"
    log: "results/logs/bam_reconstructed.log"
    shell: "python3 workflow/scripts/read_flow_reconstructed.py -i {input.bam} -c {input.cbcpath} -s {input.pbcpath} --status-group-out {output.status_group}  > {log} 2>&1"