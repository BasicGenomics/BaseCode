rule log_python_version:
    output:
        log = "results/logs/python_version.log",
    shell: """
    echo Log Python Version
    python3 --version > {output.log}
    """

rule log_version:
    output:
        log = "results/logs/basecode_version.log",
    shell: """
    echo Log BaseCode Version
    cat VERSION > {output.log}
    """


rule get_cell_barcodes:
    input: fastq_r1 = config["r1"], fastq_r2 = config["r2"]
    output: r1_out = "results/barcodes/{}_R1.fq.gz".format(config['name']),
        r2_out = "results/barcodes/{}_R2.fq.gz".format(config['name']),
        log_yaml = "results/barcodes/{}_log.yaml".format(config['name']),
        cell_barcodes =  "results/barcodes/{}_whitelist.txt".format(config['name']),
        done = "results/dones/{}_get_cell_barcodes.done".format(config['name'])
    params: barcode_cfg = BARCODE_CFG,
            prefix = config["name"]
    log: "results/logs/get_cell_barcodes.log"
    benchmark: "results/barcodes/benchmarks/get_cell_barcodes.benchmark.txt"
    shell: """
    echo Step 1/10 Get Cell Barcodes
    binaries/pipspeak_basecode -c {params.barcode_cfg} \
    -i {input.fastq_r1} -I {input.fastq_r2} -l -u 0 -p results/barcodes/{params.prefix} > {log} 2>&1
    mkdir -p results/metadata/
    cp {output.cell_barcodes} results/metadata/{params.prefix}_cell_barcodes.txt
    touch {output.done}
    """

rule make_barcode_files:
    input:  samplesheet = config["samplesheet"],
            fastq = config["r2"],
            index_sequences = "workflow/resources/index_sequences.yaml",
            cell_barcodes = "results/barcodes/{name}_whitelist.txt".format(name=config['name'])
    output: barcodes = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"]),
            sample_map = "results/metadata/{name}_sample_map.yaml".format(name=config["name"]),
            readtype_map = "results/metadata/{name}_readtype_map.yaml".format(name=config["name"]),
            dt_structure = "results/metadata/{name}_dt_structure.yaml".format(name=config["name"]),
            samplesheet_out = "results/metadata/{name}_samplesheet.csv".format(name=config["name"]),
            done = "results/dones/{name}_make_barcode_files.done".format(name=config["name"])
    params: ignore_none_flag =  '--ignore_none' if config["ignore_none"] else ''
    log: "results/logs/make_barcode_files.log"
    benchmark: "results/benchmarks/make_barcode_files.benchmark.txt"
    conda: "../envs/full.yaml"
    shell: """
    echo Step 2/10 Process Samplesheet
    python3 workflow/scripts/make_sample_files.py \
    -s {input.samplesheet} --fastq {input.fastq} \
    --index-sequences {input.index_sequences} --sample-barcodes {output.barcodes} \
    --cell-barcodes {input.cell_barcodes} --sample-map {output.sample_map} \
    --readtype-map {output.readtype_map}  --dt-structure {output.dt_structure} \
    --samplesheet-out {output.samplesheet_out} --cells {params.ignore_none_flag} > {log} 2>&1
    touch {output.done}
    """

if config["i1"] != "" and config["i2"] != "":
    checkpoint parse_fastq:
        input:  r1_in = config["r1"],
                r2_in = config["r2"],
                i1_in = config["i1"],
                i2_in = config["i2"],
                pbcpath = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"]),
                cell_barcodes = "results/barcodes/{name}_whitelist.txt".format(name=config["name"]),
                sample_map = "results/metadata/{name}_sample_map.yaml".format(name=config["name"]),
                readtype_map = "results/metadata/{name}_readtype_map.yaml".format(name=config["name"]),
                dt_structure = "results/metadata/{}_dt_structure.yaml".format(config["name"])
        output: fastq_dir = directory(FASTQ_DIR),
                done = "results/dones/{name}_parse_fastq.done".format(name=config["name"])
        log: "results/logs/parse_fastq.log"
        benchmark: "results/benchmarks/parse_fastq.benchmark.txt"
        params: comp_threads = int(config["threads"]*0.2),
                proc_threads = config["threads"]-int(config["threads"]*0.2),
                cbc_offset = config['params']['parse_fq']['cbc_offset'],
                prefix = FASTQ_DIR,
                demultiplex_flag = FLAG_demultiplex,
                run_name = config["name"]
        shell: """
        echo Step 3/10 Parse FASTQ
        binaries/parse_fastq \
        --read1 {input.r1_in} --read2 {input.r2_in} \
        --index1 {input.i1_in} --index2 {input.i2_in} \
        --cbcpath {input.cell_barcodes} \
        --pbcpath {input.pbcpath} \
        --readtype-structure {input.readtype_map} \
        --dt-structure {input.dt_structure} \
        --index-layout {config[index_layout]} \
        --sample-structure {input.sample_map} \
        --processing-threads {params.proc_threads} \
        --compression-threads {params.comp_threads} \
        --ts-sequence {config[ts_sequence]} \
        --ts-pad {config[ts_pad]} --ts-cutoff {config[ts_cutoff]} \
        --cbc-offset {params.cbc_offset} {params.demultiplex_flag} \
        --prefix {params.prefix} --run-name {params.run_name} > {log} 2>&1
        touch {output.done}
        """
else:
    checkpoint parse_fastq:
        input:  r1_in = config["r1"],
                r2_in = config["r2"],
                pbcpath = "results/metadata/{}_sample_barcodes.txt".format(config["name"]),
                cell_barcodes = "results/barcodes/{}_whitelist.txt".format(config["name"]),
                sample_map = "results/metadata/{}_sample_map.yaml".format(config["name"]),
                readtype_map = "results/metadata/{}_readtype_map.yaml".format(config["name"]),
                dt_structure = "results/metadata/{}_dt_structure.yaml".format(config["name"])
        output: fastq_dir = directory(FASTQ_DIR),
                done = "results/dones/{name}_parse_fastq.done".format(name=config["name"])
        log: "results/logs/parse_fastq.log"
        benchmark: "results/benchmarks/parse_fastq.benchmark.txt"
        params: comp_threads = int(config["threads"]*0.2),
                proc_threads = config["threads"]-int(config["threads"]*0.2),
                cbc_offset = config['params']['parse_fq']['cbc_offset'],
                prefix = FASTQ_DIR,
                demultiplex_flag = FLAG_demultiplex,
                run_name = config["name"]
        shell: """
        echo Step 3/10 Parse FASTQ
        binaries/parse_fastq \
        --read1 {input.r1_in} --read2 {input.r2_in} \
        --cbcpath {input.cell_barcodes} \
        --pbcpath {input.pbcpath} \
        --readtype-structure {input.readtype_map} \
        --dt-structure {input.dt_structure} \
        --index-layout {config[index_layout]} \
        --sample-structure {input.sample_map} \
        --processing-threads {params.proc_threads} \
        --compression-threads {params.comp_threads} \
        --ts-sequence {config[ts_sequence]} \
        --ts-pad {config[ts_pad]} \
        --ts-cutoff {config[ts_cutoff]} \
        --cbc-offset {params.cbc_offset} {params.demultiplex_flag} \
        --prefix {params.prefix} --run-name {params.run_name} > {log} 2>&1
        touch {output.done}
        """

def get_samples():
    checkpoint_output = checkpoints.parse_fastq.get().output.fastq_dir
    return [f.replace("_1.fq.gz", "") for f in os.listdir(checkpoint_output) if f.endswith("_1.fq.gz")]

def parse_fq_dir(wc):
    return checkpoints.parse_fastq.get().output.fastq_dir

rule trim_fastq:
    input:  r1 = lambda wc: os.path.join(parse_fq_dir(wc), f"{wc.sample}_1.fq.gz"),
            r2 = lambda wc: os.path.join(parse_fq_dir(wc), f"{wc.sample}_2.fq.gz")
    output:
            r1 = temp("results/intermediate/{sample}.trimmed.read1.fastq.gz"),
            r2 = temp("results/intermediate/{sample}.trimmed.read2.fastq.gz"),
            r1_short = temp("results/intermediate/{sample}.tooshort.read1.fastq.gz"),
            r2_short = temp("results/intermediate/{sample}.tooshort.read2.fastq.gz"),
            done = "results/dones/{sample}.trim_fastq.done"
    log: stdout = "results/logs/{sample}/{sample}.trim_fastq.log",
         summary = "results/summaries/{sample}/{sample}.cutadapt.json"
    benchmark: "results/benchmarks/{sample}/{sample}.trim_fastq.benchmark.txt"
    conda: "../envs/full.yaml"
    threads: config["threads"]
    shell: """
    echo Step 4/10 Trim FASTQ
    cutadapt -j {threads} --json {log.summary} {config[params][cutadapt]} \
    --too-short-output {output.r1_short} --too-short-paired-output {output.r2_short} \
    -o {output.r1} -p {output.r2} {input.r1} {input.r2} > {log.stdout} 2>&1
    touch {output.done}
    """


rule map_reads:
    input:  r1 = "results/intermediate/{sample}.trimmed.read1.fastq.gz",
            r2 = "results/intermediate/{sample}.trimmed.read2.fastq.gz"
    output: bam = temp("results/intermediate/{sample}.trimmed.aligned.bam"),
            done = "results/dones/{sample}.map_reads.done"
    log: stdout = "results/logs/{sample}/{sample}.map_reads.log",
        summary = "results/summaries/{sample}/{sample}.hisat2.summary.txt"
    benchmark: "results/benchmarks/{sample}.map_reads.benchmark.txt"
    params: splicesites = SPLICESITES, genomeref = GENOMEREF, cores_hisat = cores_hisat, cores_samtools = cores_samtools
    threads: min(config["threads"], 64)
    conda: "../envs/full.yaml"
    shell: """
    echo Step 5/10 Map Reads
    binaries/hisat-3n --new-summary --summary-file {log.summary} \
    {config[params][hisat3n]} -p {params.cores_hisat} \
    --known-splicesite-infile {params.splicesites} \
    -x {params.genomeref} -1 {input.r1} -2 {input.r2} | samtools view -F 256 -b -@ {params.cores_samtools} -o {output.bam} > {log.stdout} 2>&1
    touch {output.done}
    """

rule split_bam_by_strand:
    input: "results/intermediate/{sample}.trimmed.aligned.bam"
    output: pstrand = temp("results/intermediate/{sample}.trimmed.aligned.pstrand.bam"),
            mstrand = temp("results/intermediate/{sample}.trimmed.aligned.mstrand.bam"),
            nostrand = temp("results/intermediate/{sample}.trimmed.aligned.nostrand.bam"),
            done = "results/dones/{sample}.split_bam_by_strand.done"
    log: "results/logs/{sample}/{sample}.split_bam_by_strand.log"
    benchmark: "results/benchmarks/{sample}.split_bam_by_strand.benchmark.txt"
    conda: "../envs/full.yaml"
    shell: """
    echo Step 6/10 Assign Genes
    binaries/move_tags --input {input} > {log} 2>&1
    touch {output.done}
    """

rule assign_genes_exon:
    input: pstrand = "results/intermediate/{sample}.trimmed.aligned.pstrand.bam",
            mstrand = "results/intermediate/{sample}.trimmed.aligned.mstrand.bam",
            nostrand = "results/intermediate/{sample}.trimmed.aligned.nostrand.bam"
    output: pstrand = temp("results/intermediate/{sample}.trimmed.aligned.pstrand.bam.featureCounts.bam"),
            mstrand = temp("results/intermediate/{sample}.trimmed.aligned.mstrand.bam.featureCounts.bam"),
            nostrand = temp("results/intermediate/{sample}.trimmed.aligned.nostrand.bam.featureCounts.bam"),
            done = "results/dones/{sample}.assign_genes_exon.done"
    log: "results/logs/{sample}/{sample}.assign_genes.log"
    benchmark: "results/benchmarks/{sample}/{sample}.assign_genes.benchmark.txt"
    params: gtffile = "{}.gff3".format(GTFFILE),
            gtffile_positive = "{}.positive.gff3".format(GTFFILE),
            gtffile_negative = "{}.negative.gff3".format(GTFFILE)
    conda: "../envs/full.yaml"
    threads: min(config["threads"], 64)
    shell:"""
    binaries/featureCounts -t exon --primary -g {config[gff_gene_identifier]} -T {threads} -R BAM -p --countReadPairs -O -M --largestOverlap --fracOverlap 0.1 -a {params.gtffile_positive} -o results/intermediate/pos.tmp {input.pstrand} >> {log} 2>&1
    binaries/featureCounts -t exon --primary -g {config[gff_gene_identifier]} -T {threads} -R BAM -p --countReadPairs -O -M --largestOverlap --fracOverlap 0.1 -a {params.gtffile_negative} -o results/intermediate/neg.tmp {input.mstrand} >> {log} 2>&1
    binaries/featureCounts -t exon --primary -g {config[gff_gene_identifier]} -T {threads} -R BAM -p --countReadPairs -O -M --largestOverlap --fracOverlap 0.1 -a {params.gtffile} -o results/intermediate/no.tmp {input.nostrand} >> {log} 2>&1
    rm results/intermediate/pos.tmp results/intermediate/neg.tmp results/intermediate/no.tmp
    mkdir -p results/.tmp_bgab/
    touch {output.done}
    """

rule rename_tags_exon:
    input:  pstrand = "results/intermediate/{sample}.trimmed.aligned.pstrand.bam.featureCounts.bam",
            mstrand = "results/intermediate/{sample}.trimmed.aligned.mstrand.bam.featureCounts.bam",
            nostrand = "results/intermediate/{sample}.trimmed.aligned.nostrand.bam.featureCounts.bam"
    output: pstrand = temp("results/intermediate/{sample}.trimmed.aligned.pstrand.Exon.bam"),
            mstrand = temp("results/intermediate/{sample}.trimmed.aligned.mstrand.Exon.bam"),
            nostrand = temp("results/intermediate/{sample}.trimmed.aligned.nostrand.Exon.bam"),
            done = "results/dones/{sample}.rename_tags_exon.done"
    log: "results/logs/{sample}/{sample}.rename_tags_exon.log"
    benchmark: "results/benchmarks/{sample}/{sample}.rename_tags_exon.benchmark.txt"
    conda: "../envs/full.yaml"
    shell:"""
    binaries/rename_tags --input {input.pstrand} --output {output.pstrand} >> {log} 2>&1
    binaries/rename_tags --input {input.mstrand} --output {output.mstrand} >> {log} 2>&1
    binaries/rename_tags --input {input.nostrand} --output {output.nostrand} >> {log} 2>&1
    touch {output.done}
    """

rule assign_genes_intron:
    input: pstrand = "results/intermediate/{sample}.trimmed.aligned.pstrand.Exon.bam",
            mstrand = "results/intermediate/{sample}.trimmed.aligned.mstrand.Exon.bam",
            nostrand = "results/intermediate/{sample}.trimmed.aligned.nostrand.Exon.bam"
    output: pstrand = temp("results/intermediate/{sample}.trimmed.aligned.pstrand.Exon.bam.featureCounts.bam"),
            mstrand = temp("results/intermediate/{sample}.trimmed.aligned.mstrand.Exon.bam.featureCounts.bam"),
            nostrand = temp("results/intermediate/{sample}.trimmed.aligned.nostrand.Exon.bam.featureCounts.bam"),
            done = "results/dones/{sample}.assign_genes_intron.done"
    log: "results/logs/{sample}/{sample}.assign_genes.log"
    benchmark: "results/benchmarks/{sample}/{sample}.assign_genes.benchmark.txt"
    params: gtffile = "{}.gff3".format(GTFFILE),
            gtffile_positive = "{}.positive.gff3".format(GTFFILE),
            gtffile_negative = "{}.negative.gff3".format(GTFFILE)
    conda: "../envs/full.yaml"
    threads: min(config["threads"], 64)
    shell:"""
    binaries/featureCounts -t intron --primary -g {config[gff_gene_identifier]} -T {threads} -R BAM -p --countReadPairs -O -M --largestOverlap --fracOverlap 0.1 -a {params.gtffile_positive} -o results/intermediate/pos.tmp {input.pstrand} >> {log} 2>&1
    binaries/featureCounts -t intron --primary -g {config[gff_gene_identifier]} -T {threads} -R BAM -p --countReadPairs -O -M --largestOverlap --fracOverlap 0.1 -a {params.gtffile_negative} -o results/intermediate/neg.tmp {input.mstrand} >> {log} 2>&1
    binaries/featureCounts -t intron --primary -g {config[gff_gene_identifier]} -T {threads} -R BAM -p --countReadPairs -O -M --largestOverlap --fracOverlap 0.1 -a {params.gtffile} -o results/intermediate/no.tmp {input.nostrand} >> {log} 2>&1
    rm results/intermediate/pos.tmp results/intermediate/neg.tmp results/intermediate/no.tmp
    mkdir -p results/.tmp_bgab/
    touch {output.done}
    """

rule rename_tags_intron:
    input:  pstrand = "results/intermediate/{sample}.trimmed.aligned.pstrand.Exon.bam.featureCounts.bam",
            mstrand = "results/intermediate/{sample}.trimmed.aligned.mstrand.Exon.bam.featureCounts.bam",
            nostrand = "results/intermediate/{sample}.trimmed.aligned.nostrand.Exon.bam.featureCounts.bam"
    output: pstrand = temp("results/intermediate/{sample}.trimmed.aligned.pstrand.GeneTagged.bam"),
            mstrand = temp("results/intermediate/{sample}.trimmed.aligned.mstrand.GeneTagged.bam"),
            nostrand = temp("results/intermediate/{sample}.trimmed.aligned.nostrand.GeneTagged.bam"),
            done = "results/dones/{sample}.rename_tags_intron.done"
    log: "results/logs/{sample}/{sample}.rename_tags_intron.log"
    benchmark: "results/benchmarks/{sample}/{sample}.rename_tags_intron.benchmark.txt"
    conda: "../envs/full.yaml"
    shell:"""
    binaries/rename_tags --input {input.pstrand} --output {output.pstrand} --intron-mode >> {log} 2>&1
    binaries/rename_tags --input {input.mstrand} --output {output.mstrand} --intron-mode >> {log} 2>&1
    binaries/rename_tags --input {input.nostrand} --output {output.nostrand} --intron-mode >> {log} 2>&1
    touch {output.done}
    """

rule concatenate_and_sort:
    input: pstrand = "results/intermediate/{sample}.trimmed.aligned.pstrand.GeneTagged.bam",
            mstrand = "results/intermediate/{sample}.trimmed.aligned.mstrand.GeneTagged.bam",
            nostrand = "results/intermediate/{sample}.trimmed.aligned.nostrand.GeneTagged.bam"
    log: "results/logs/{sample}/{sample}.concatenate_and_sort.log"
    benchmark: "results/benchmarks/{sample}/{sample}.concatenate_and_sort.benchmark.txt"
    output: bam = temp("results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted.bam"),
            done = "results/dones/{sample}.concatenate_and_sort.done"
    conda: "../envs/full.yaml"
    threads: config["threads"]
    shell: """
    samtools cat {input.nostrand} {input.pstrand} {input.mstrand} | samtools sort -m 1000M -@ {threads} -T results/.tmp_bgab/sorttmp. -o {output.bam} >> {log} 2>&1
    touch {output.done}
    """

rule first_index:
    input: "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted.bam"
    output: bai = temp("results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted.bam.bai"),
            done = "results/dones/{sample}.first_index.done"
    threads: config["threads"]
    log: "results/logs/{sample}/{sample}.first_index.log"
    benchmark: "results/benchmarks/{sample}/{sample}.first_index.benchmark.txt"
    shell: """
    samtools index -@ {config[threads]} {input}
    touch {output.done}
    """

rule error_correct_umis:
    input: bam = "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted.bam",
        bai = "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted.bam.bai"
    output: bam = "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted_umicorrected.bam",
            done = "results/dones/{sample}.error_correct_umis.done"
    log: "results/logs/{sample}/{sample}.error_correct_umis.log"
    benchmark: "results/benchmarks/{sample}/{sample}.error_correct_umis.benchmark.txt"
    params: gtffile = GTFFILE
    conda: "../envs/full.yaml"
    shell: """
    echo Step 7/10 Error Correct UMIs
    python3 workflow/scripts/error_correct_umis.py {input.bam} {params.gtffile}.gff3 {config[threads]} &> {log}
    touch {output.done}
    """

rule second_index:
    input: "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted_umicorrected.bam"
    output: bai = "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted_umicorrected.bam.bai",
            done = "results/dones/{sample}.second_index.done"
    params:
        extra="",
    threads: config["threads"]
    benchmark: "results/benchmarks/{sample}/{sample}.second_index.benchmark.txt"
    shell: """
    samtools index -@ {threads} {params.extra} {input}
    touch {output.done}
    """

if config["reverse"]:
    rule reconstruct:
        input: bam = "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted_umicorrected.bam",
            bai = "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted_umicorrected.bam.bai",
            sample_map = "results/metadata/{name}_sample_map.yaml".format(name=config["name"])
        output: bam = temp("results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted_umicorrected.reconstructed.bam"),
                done = "results/dones/{sample}.reconstruct.done"
        params: gtffile = GTFFILE, cores_hisat = cores_hisat, cores_samtools = cores_samtools
        log: "results/logs/{sample}/{sample}.reconstruct.log"
        threads: min(config["threads"], 64)
        benchmark: "results/benchmarks/{sample}/{sample}.reconstruction.benchmark.txt"
        shell: """
        echo Step 8/10 Reconstruct Molecules
        mkdir -p results/tmp/
        binaries/basic_reconstruction --input {input.bam} --output {output.bam} --gtf {params.gtffile}.gff3 --sample-map {input.sample_map} --threads {threads} --gene-identifier {config[gff_gene_identifier]} --reverse | samtools view -F 256 -b -@ {params.cores_samtools} -o {output.bam} > {log} 2>&1
        touch {output.done}
        """
else:
    rule reconstruct:
        input: bam = "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted_umicorrected.bam",
            bai = "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted_umicorrected.bam.bai",
            sample_map = "results/metadata/{name}_sample_map.yaml".format(name=config["name"])
        output: bam = temp("results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted_umicorrected.reconstructed.bam"),
                done = "results/dones/{sample}.reconstruct.done"
        params: gtffile = GTFFILE, cores_hisat = cores_hisat, cores_samtools = cores_samtools
        log: "results/logs/{sample}/{sample}.reconstruct.log"
        threads: min(config["threads"], 64)
        benchmark: "results/benchmarks/{sample}/{sample}.reconstruction.benchmark.txt"
        shell: """
        echo Step 8/10 Reconstruct Molecules
        mkdir -p results/tmp/
        binaries/basic_reconstruction --input {input.bam} --output {output.bam} --gtf {params.gtffile}.gff3 --sample-map {input.sample_map} --threads {threads} --gene-identifier {config[gff_gene_identifier]}  | samtools view -F 256 -b -@ {params.cores_samtools} -o {output.bam} > {log} 2>&1
        touch {output.done}
        """

rule sort_reconstructed:
    input: "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted_umicorrected.reconstructed.bam"
    output: bam = "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted_umicorrected.reconstructed.sorted.bam",
            done = "results/dones/{sample}.sort_reconstructed.done"
    threads: config["threads"]
    log: "results/logs/{sample}/{sample}.sort_reconstructed.log"
    benchmark: "results/benchmarks/{sample}/{sample}.sort_reconstructed.benchmark.txt"
    params:
        extra="-m 1000M",
    conda: "../envs/full.yaml"
    shell: """
    samtools sort -@ {threads} -o {output.bam} {params.extra} -T . {input} &> {log}
    touch {output.done}
    """

rule index_reconstructed:
    input: "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted_umicorrected.reconstructed.sorted.bam"
    output: bai = "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted_umicorrected.reconstructed.sorted.bam.bai",
            done = "results/dones/{sample}.index_reconstructed.done"
    log: "results/logs/{sample}/{sample}.index_reconstructed.log"
    benchmark: "results/benchmarks/{sample}/{sample}.index_reconstructed.benchmark.txt"
    params:
        extra="",
    threads: config["threads"]
    shell: """
    samtools index -@ {threads} {params.extra} {input}
    touch {output.done}
    """

rule stitch_reconstruction:
        input: bam =  "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted_umicorrected.reconstructed.sorted.bam",
            bai = "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted_umicorrected.reconstructed.sorted.bam.bai"
        output: bam = temp("results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted_umicorrected.stitched.bam"),
                done = "results/dones/{sample}.stitch_reconstruction.done"
        threads: config["threads"]
        params: gtffile = GTFFILE
        log: "results/logs/{sample}/{sample}.stitch_reconstruction.log"
        benchmark: "results/benchmarks/{sample}/{sample}.stitch_reconstruction.benchmark.txt"
        conda: "../envs/full.yaml"
        shell: """
        echo Step 9/10 Stitch Molecules
        python3 workflow/scripts/stitcher.py --input {input.bam} --output {output.bam} --gtf {params.gtffile}.gff3 --threads {threads} --cell-tag CB --UMI-tag RM --gene-identifier {config[gff_gene_identifier]} --only-molecules >> {log} 2>&1
        touch {output.done}
        """

rule sorted_stitched:
    input: "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted_umicorrected.stitched.bam"
    output: bam = "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted_umicorrected.stitched.sorted.bam",
            done = "results/dones/{sample}.sorted_stitched.done"
    threads: config["threads"]
    log: "results/logs/{sample}/{sample}.sort_stitched.log"
    benchmark: "results/benchmarks/{sample}/{sample}.sorted_stitched.benchmark.txt"
    params:
        extra="-m 1000M",
    conda: "../envs/full.yaml"
    shell: """
    samtools sort -@ {threads} -o {output.bam} {params.extra} -T . {input} &> {log}
    touch {output.done}
    """

rule index_stitched:
    input: "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted_umicorrected.stitched.sorted.bam"
    output: bai = "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted_umicorrected.stitched.sorted.bam.bai",
            done = "results/dones/{sample}.index_stitched.done"
    params:
        extra="",
    threads: config["threads"]
    benchmark: "results/benchmarks/{sample}/{sample}.index_stitched.benchmark.txt"
    shell: """
    samtools index -@ {threads} {params.extra} {input}
    touch {output.done}
    """

rule make_molecule_bams:
    input: bam = "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted_umicorrected.stitched.sorted.bam", bai = "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted_umicorrected.stitched.sorted.bam.bai"
    output: molecules_out = temp("results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted_umicorrected.stitched.molecules.bam"),
            done = "results/dones/{sample}.make_molecule_bams.done"
    benchmark: "results/benchmarks/{sample}/{sample}.make_molecule_bams.benchmark.txt"
    conda: "../envs/full.yaml"
    shell: """
    python3 workflow/scripts/filter_stitched_bam.py --input {input.bam} --molecules-out {output.molecules_out}
    touch {output.done}
    """

rule sort_molecule_bams:
    input: molecules_out = "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted_umicorrected.stitched.molecules.bam"
    output: molecules_out = "results/{sample}.stitched.molecules.sorted.bam",
            done = "results/done/{sample}.processing.done"
    log: "results/logs/{sample}/{sample}.sort_molecules.log"
    benchmark: "results/benchmarks/{sample}/{sample}.sort_molecule_bams.benchmark.txt"
    conda: "../envs/full.yaml"
    shell:"""
    echo Step 10/10 Creating Final Output File
    samtools sort -m 1000M -@ {config[threads]} -T results/.tmp_bgab/sorttmp. -o {output.molecules_out} {input.molecules_out} &> {log}
    samtools index -@ {config[threads]} {output.molecules_out}
    touch {output.done}
    """
