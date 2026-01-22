rule log_python_version:
    conda: "../envs/full.yaml"
    output: "results/logs/python_version.log"
    shell: "python3 --version > {output}"

rule get_cell_barcodes:
    input: fastq_r1 = config["r1"], fastq_r2 = config["r2"]
    output: r1_out = "results/barcodes/{}_R1.fq.gz".format(config['name']),
        r2_out = "results/barcodes/{}_R2.fq.gz".format(config['name']),
        log_yaml = "results/barcodes/{}_log.yaml".format(config['name']),
        cell_barcodes =  "results/barcodes/{}_whitelist.txt".format(config['name'])
    params: barcode_cfg = BARCODE_CFG,
            prefix = config["name"]
    log: "results/logs/get_cell_barcodes.log"
    benchmark: "results/barcodes/benchmarks/get_cell_barcodes.benchmark.txt"
    shell: "echo Get barcode whitelist && binaries/pipspeak_basecode -c {params.barcode_cfg} -i {input.fastq_r1} -I {input.fastq_r2} -l -u 0 -p results/barcodes/{params.prefix} > {log} 2>&1"

rule make_barcode_files:
    input: samplesheet = config["samplesheet"], fastq = config["r2"], index_sequences = "workflow/resources/index_sequences.yaml", cell_barcodes = "results/barcodes/{name}_whitelist.txt".format(name=config['name'])
    output: barcodes = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"]),
        sample_map = "results/metadata/{name}_sample_map.yaml".format(name=config["name"]),
        readtype_map = "results/metadata/{name}_readtype_map.yaml".format(name=config["name"]),
        dt_structure = "results/metadata/{name}_dt_structure.yaml".format(name=config["name"]),
        samplesheet_out = "results/metadata/{name}_samplesheet.csv".format(name=config["name"])
    params: barcode_cfg = BARCODE_CFG
    log: "results/logs/make_barcode_files.log"
    conda: "../envs/full.yaml"
    shell: "echo Process Samplesheet && python3 workflow/scripts/make_sample_files.py -s {input.samplesheet} --fastq {input.fastq} --index-sequences {input.index_sequences} --sample-barcodes {output.barcodes} --cell-barcodes {input.cell_barcodes} --sample-map {output.sample_map} --readtype-map {output.readtype_map}  --dt-structure {output.dt_structure} --samplesheet-out {output.samplesheet_out} > {log} 2>&1"

rule parse_fastq:
        input: r1_in = config["r1"], r2_in = config["r2"], pbcpath = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"]), cell_barcodes = "results/barcodes/{name}_whitelist.txt".format(name=config["name"]), sample_map = "results/metadata/{name}_sample_map.yaml".format(name=config["name"]), readtype_map = "results/metadata/{name}_readtype_map.yaml".format(name=config["name"]), dt_structure = "results/metadata/{}_dt_structure.yaml".format(config["name"])
        output: r1_out = temp("results/intermediate/{name}.read1.fastq.gz".format(name=config["name"])),
                r2_out = temp("results/intermediate/{name}.read2.fastq.gz".format(name=config["name"]))
        log: "results/logs/parse_fastq.log"
        benchmark: "results/benchmarks/parse_fastq.benchmark.txt"
        params: comp_threads = int(config["threads"]*0.2),
                proc_threads = config["threads"]-int(config["threads"]*0.2)
        shell: "echo Parse FASTQ && binaries/parse_fastq --read1 {input.r1_in} --read2 {input.r2_in} --r1-out {output.r1_out} --r2-out {output.r2_out} --cbcpath {input.cell_barcodes} --pbcpath {input.pbcpath} --readtype-structure {input.readtype_map} --dt-structure {input.dt_structure} --index-layout {config[index_layout]} --sample-structure {input.sample_map} --processing-threads {params.proc_threads} --compression-threads {params.comp_threads} --ts-sequence {config[ts_sequence]} --ts-pad {config[ts_pad]} --ts-cutoff {config[ts_cutoff]} > {log} 2>&1"