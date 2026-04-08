rule fastq_processed:
    input: r1 = "results/intermediate/{name}.read1.fastq.gz".format(name=config["name"]),
           r2 = "results/intermediate/{name}.read2.fastq.gz".format(name=config["name"]),
           cbcpath = "results/metadata/{name}_cell_barcodes.txt".format(name=config["name"]),
           pbcpath = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"])
    output: "results/read_flow_files/{name}_fastq_processed_stats.json".format(name=config["name"])
    log: "results/logs/fastq_processed.log"
    threads: config["threads"]
    shell: "binaries/analyze_fastq --read1 {input.r1} --read2 {input.r2} --cbcpath {input.cbcpath} --pbcpath {input.pbcpath} --threads {threads} --output {output} > {log} 2>&1"

rule fastq_trimmed:
    input: r1 = "results/intermediate/{name}.trimmed.read1.fastq.gz".format(name=config["name"]),
           r2 = "results/intermediate/{name}.trimmed.read2.fastq.gz".format(name=config["name"]),
           cbcpath = "results/metadata/{name}_cell_barcodes.txt".format(name=config["name"]),
           pbcpath = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"])
    output: "results/read_flow_files/{name}_fastq_trimmed_stats.json".format(name=config["name"])
    log: "results/logs/fastq_trimmed.log"
    threads: config["threads"]
    shell: "binaries/analyze_fastq --read1 {input.r1} --read2 {input.r2} --cbcpath {input.cbcpath} --pbcpath {input.pbcpath} --threads {threads} --output {output} > {log} 2>&1"

rule fastq_too_short:
    input: r1 = "results/intermediate/{name}.tooshort.read1.fastq.gz".format(name=config["name"]),
           r2 = "results/intermediate/{name}.tooshort.read2.fastq.gz".format(name=config["name"]),
           cbcpath = "results/metadata/{name}_cell_barcodes.txt".format(name=config["name"]),
           pbcpath = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"])
    output: "results/read_flow_files/{name}_fastq_tooshort_stats.json".format(name=config["name"])
    log: "results/logs/fastq_too_short.log"
    threads: config["threads"]
    shell: "binaries/analyze_fastq --read1 {input.r1} --read2 {input.r2} --cbcpath {input.cbcpath} --pbcpath {input.pbcpath} --threads {threads} --output {output} > {log} 2>&1"

rule bam_geneassigned:
    input: bam = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.bam".format(name=config["name"]),
           bai = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.bam.bai".format(name=config["name"]),
           cbcpath = "results/metadata/{name}_cell_barcodes.txt".format(name=config["name"]),
           pbcpath = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"])
    output: "results/read_flow_files/{name}_mapping_group.csv".format(name=config["name"])
    log: "results/logs/bam_geneassigned.log"
    shell: "python3 workflow/scripts/read_flow_mapped.py -i {input.bam} -c {input.cbcpath} -s {input.pbcpath} --mapping-group-out {output} > {log} 2>&1"

rule bam_reconstructed:
    input: bam = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.reconstructed.sorted.bam".format(name=config["name"]),
           bai = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.reconstructed.sorted.bam.bai".format(name=config["name"]),
           cbcpath = "results/metadata/{name}_cell_barcodes.txt".format(name=config["name"]), 
           pbcpath = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"])
    output: "results/read_flow_files/{name}_status_group.csv".format(name=config["name"])
    log: "results/logs/bam_reconstructed.log"
    shell: "python3 workflow/scripts/read_flow_reconstructed.py -i {input.bam} -c {input.cbcpath} -s {input.pbcpath} --status-group-out {output} > {log} 2>&1"