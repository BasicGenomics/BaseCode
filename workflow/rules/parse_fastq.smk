
rule parse_fastq:
        input: r1_in = config["r1"], r2_in = config["r2"], pbcpath = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"]), cell_barcodes = "results/barcodes/{name}_whitelist.txt".format(name=config["name"]), sample_map = "results/metadata/{name}_sample_map.yaml".format(name=config["name"]), readtype_map = "results/metadata/{name}_readtype_map.yaml".format(name=config["name"]), dt_structure = "results/metadata/{}_dt_structure.yaml".format(config["name"])
        output: r1_out = temp("results/intermediate/{name}.read1.fastq.gz".format(name=config["name"])),
                r2_out = temp("results/intermediate/{name}.read2.fastq.gz".format(name=config["name"]))
        log: "results/logs/parse_fastq.log"
        benchmark: "results/benchmarks/parse_fastq.benchmark.txt"
        params: comp_threads = int(config["threads"]*0.2),
                proc_threads = config["threads"]-int(config["threads"]*0.2),
                cbc_offset = config['params']['cbc_offset']
        shell: "echo Parse FASTQ && binaries/parse_fastq --read1 {input.r1_in} --read2 {input.r2_in} --r1-out {output.r1_out} --r2-out {output.r2_out} --cbcpath {input.cell_barcodes} --pbcpath {input.pbcpath} --readtype-structure {input.readtype_map} --dt-structure {input.dt_structure} --index-layout {config[index_layout]} --sample-structure {input.sample_map} --processing-threads {params.proc_threads} --compression-threads {params.comp_threads} --ts-sequence {config[ts_sequence]} --ts-pad {config[ts_pad]} --ts-cutoff {config[ts_cutoff]} --cbc-offset {params.cbc_offset}> {log} 2>&1"