OUT_DIR = "results/read_flow_files"

rule fastq_processed:
    input: r1 = lambda wc: os.path.join(parse_fq_dir(wc), f"{wc.sample}_1.fq.gz"),
           r2 = lambda wc: os.path.join(parse_fq_dir(wc), f"{wc.sample}_2.fq.gz"),
           cbcpath = "results/metadata/{name}_cell_barcodes.txt".format(name=config["name"]),
           pbcpath = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"])
    output: "results/read_flow_files/{sample}/{sample}_fastq_processed_stats.json"
    threads: int(config["threads"])-1
    log: "results/read_flow_files/logs/{sample}/{sample}.fastq_processed.log"
    shell: "binaries/analyze_fastq --read1 {input.r1} --read2 {input.r2} --cbcpath {input.cbcpath} --pbcpath {input.pbcpath} --threads {config[threads]} --output {output}  > {log} 2>&1"

rule fastq_trimmed:
    input: r1 = "results/intermediate/{sample}.trimmed.read1.fastq.gz",
           r2 = "results/intermediate/{sample}.trimmed.read2.fastq.gz",
           cbcpath = "results/metadata/{name}_cell_barcodes.txt".format(name=config["name"]),
           pbcpath = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"])
    output: "results/read_flow_files/{sample}/{sample}_fastq_trimmed_stats.json"
    threads: int(config["threads"])-1
    log: "results/read_flow_files/logs/{sample}/{sample}.fastq_trimmed.log"
    shell: "binaries/analyze_fastq --read1 {input.r1} --read2 {input.r2} --cbcpath {input.cbcpath} --pbcpath {input.pbcpath} --threads {config[threads]} --output {output}  > {log} 2>&1"

rule fastq_too_short:
    input: r1 = "results/intermediate/{sample}.tooshort.read1.fastq.gz",
           r2 = "results/intermediate/{sample}.tooshort.read2.fastq.gz",
           cbcpath = "results/metadata/{name}_cell_barcodes.txt".format(name=config["name"]),
           pbcpath = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"])
    output: "results/read_flow_files/{sample}/{sample}_fastq_tooshort_stats.json"
    threads: int(config["threads"])-1
    log: "results/read_flow_files/logs/{sample}/{sample}.fastq_too_short.log"
    shell: "binaries/analyze_fastq --read1 {input.r1} --read2 {input.r2} --cbcpath {input.cbcpath} --pbcpath {input.pbcpath} --threads {config[threads]} --output {output}  > {log} 2>&1"

rule bam_geneassigned:
    input: bam = "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted.bam",
           cbcpath = "results/metadata/{name}_cell_barcodes.txt".format(name=config["name"]),
           pbcpath = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"])
    output: mapping_group = "results/read_flow_files/{sample}/{sample}_mapping_group.csv"
    conda: "../envs/full.yaml"
    log: "results/read_flow_files/logs/{sample}/{sample}.bam_geneassigned.log"
    shell: "python3 workflow/scripts/read_flow_mapped.py -i {input.bam} -c {input.cbcpath} -s {input.pbcpath} --mapping-group-out {output.mapping_group}  > {log} 2>&1"


def input_bam(wc):
    
    if config['params_key']=='bulk':
        suf = ".reads.aligned_trimmed_genetagged_sorted.reconstructed.sorted.bam"
    else:
        suf = ".reads.aligned_trimmed_genetagged_sorted_umicorrected.reconstructed.sorted.bam"
    
    input_bam_list = os.path.join("results/intermediate/", f"{wc.sample}{suf}")
    
    return input_bam_list
    
rule bam_reconstructed:
    input:  bam = input_bam,
            cbcpath = "results/metadata/{name}_cell_barcodes.txt".format(name=config["name"]), 
            pbcpath = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"])
    output: status_group = "results/read_flow_files/{sample}/{sample}_status_group.csv"
    conda: "../envs/full.yaml"
    log: "results/read_flow_files/logs/{sample}{sample}.bam_reconstructed.log"
    shell: "python3 workflow/scripts/read_flow_reconstructed.py -i {input.bam} -c {input.cbcpath} -s {input.pbcpath} --status-group-out {output.status_group}  > {log} 2>&1"
