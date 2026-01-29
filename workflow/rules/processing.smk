rule log_python_version:
    conda: "../envs/full.yaml"
    output: "results/logs/python_version.log"
    shell: "python3 --version > {output}"

if config["ignore_none"]:
    if config["i1"] != "" and config["i2"] != "":
        rule make_barcode_files:
            input: samplesheet = config["samplesheet"], fastq_i1 = config["i1"], fastq_i2 = config["i2"], index_sequences = "workflow/resources/index_sequences.yaml"
            output: barcodes = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"]),
                cell_barcodes = "results/metadata/{name}_cell_barcodes.txt".format(name=config["name"]),
                sample_map = "results/metadata/{name}_sample_map.yaml".format(name=config["name"]),
                readtype_map = "results/metadata/{name}_readtype_map.yaml".format(name=config["name"]),
                samplesheet_out = "results/metadata/{name}_samplesheet.csv".format(name=config["name"])
            log: "results/logs/make_barcode_files.log"
            conda: "../envs/full.yaml"
            shell: "echo Process Samplesheet && python3 workflow/scripts/make_sample_files.py -s {input.samplesheet} --fastq {input.fastq_i1} {input.fastq_i2} --index-sequences {input.index_sequences} --sample-barcodes {output.barcodes} --cell-barcodes {output.cell_barcodes} --sample-map {output.sample_map} --readtype-map {output.readtype_map} --samplesheet-out {output.samplesheet_out} --ignore-none  > {log} 2>&1"
    else:
        rule make_barcode_files:
            input: samplesheet = config["samplesheet"], fastq = config["r2"], index_sequences = "workflow/resources/index_sequences.yaml"
            output: barcodes = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"]),
                cell_barcodes = "results/metadata/{name}_cell_barcodes.txt".format(name=config["name"]),
                sample_map = "results/metadata/{name}_sample_map.yaml".format(name=config["name"]),
                readtype_map = "results/metadata/{name}_readtype_map.yaml".format(name=config["name"]),
                samplesheet_out = "results/metadata/{name}_samplesheet.csv".format(name=config["name"])
            log: "results/logs/make_barcode_files.log"
            conda: "../envs/full.yaml"
            shell: "echo Process Samplesheet && python3 workflow/scripts/make_sample_files.py -s {input.samplesheet} --fastq {input.fastq} --index-sequences {input.index_sequences} --sample-barcodes {output.barcodes} --cell-barcodes {output.cell_barcodes} --sample-map {output.sample_map} --readtype-map {output.readtype_map} --samplesheet-out {output.samplesheet_out} --ignore-none  > {log} 2>&1"
else:
    if config["i1"] != "" and config["i2"] != "":
        rule make_barcode_files:
            input: samplesheet = config["samplesheet"], fastq_i1 = config["i1"], fastq_i2 = config["i2"], index_sequences = "workflow/resources/index_sequences.yaml"
            output: barcodes = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"]),
                cell_barcodes = "results/metadata/{name}_cell_barcodes.txt".format(name=config["name"]),
                sample_map = "results/metadata/{name}_sample_map.yaml".format(name=config["name"]),
                readtype_map = "results/metadata/{name}_readtype_map.yaml".format(name=config["name"]),
                samplesheet_out = "results/metadata/{name}_samplesheet.csv".format(name=config["name"])
            log: "results/logs/make_barcode_files.log"
            conda: "../envs/full.yaml"
            shell: "echo Process Samplesheet && python3 workflow/scripts/make_sample_files.py -s {input.samplesheet} --fastq {input.fastq_i1} {input.fastq_i2} --index-sequences {input.index_sequences} --sample-barcodes {output.barcodes} --cell-barcodes {output.cell_barcodes} --sample-map {output.sample_map} --readtype-map {output.readtype_map} --samplesheet-out {output.samplesheet_out}  > {log} 2>&1"
    else:
        rule make_barcode_files:
            input: samplesheet = config["samplesheet"], fastq = config["r2"], index_sequences = "workflow/resources/index_sequences.yaml"
            output: barcodes = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"]),
                cell_barcodes = "results/metadata/{name}_cell_barcodes.txt".format(name=config["name"]),
                sample_map = "results/metadata/{name}_sample_map.yaml".format(name=config["name"]),
                readtype_map = "results/metadata/{name}_readtype_map.yaml".format(name=config["name"]),
                samplesheet_out = "results/metadata/{name}_samplesheet.csv".format(name=config["name"])
            log: "results/logs/make_barcode_files.log"
            conda: "../envs/full.yaml"
            shell: "echo Process Samplesheet && python3 workflow/scripts/make_sample_files.py -s {input.samplesheet} --fastq {input.fastq} --index-sequences {input.index_sequences} --sample-barcodes {output.barcodes} --cell-barcodes {output.cell_barcodes} --sample-map {output.sample_map} --readtype-map {output.readtype_map} --samplesheet-out {output.samplesheet_out}   > {log} 2>&1"

if config["i1"] != "" and config["i2"] != "":
    checkpoint parse_fastq:
        input: r1_in = config["r1"], r2_in = config["r2"], i1_in = config["i1"], i2_in = config["i2"], pbcpath = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"]), cell_barcodes = "results/metadata/{name}_cell_barcodes.txt".format(name=config["name"]), sample_map = "results/metadata/{name}_sample_map.yaml".format(name=config["name"]), readtype_map = "results/metadata/{name}_readtype_map.yaml".format(name=config["name"]), dt_structure = "results/metadata/{}_dt_structure.yaml".format(config["name"])
        output: directory(FASTQ_DIR)
        log: "results/logs/parse_fastq.log"
        benchmark: "results/benchmarks/parse_fastq.benchmark.txt"
        params: comp_threads = int(config["threads"]*0.2),
                proc_threads = config["threads"]-int(config["threads"]*0.2),
                cbc_offset = config['params']['parse_fq']['cbc_offset'],
                prefix = FASTQ_DIR
        shell: "echo Parse FASTQ && binaries/parse_fastq --read1 {input.r1_in} --read2 {input.r2_in} --index1 {input.i1_in} --index2 {input.i2_in} --cbcpath {input.cell_barcodes} --pbcpath {input.pbcpath} --readtype-structure {input.readtype_map} --dt-structure {input.dt_structure} --index-layout {config[index_layout]} --sample-structure {input.sample_map} --processing-threads {params.proc_threads} --compression-threads {params.comp_threads} --ts-sequence {config[ts_sequence]} --ts-pad {config[ts_pad]} --ts-cutoff {config[ts_cutoff]} --cbc-offset {params.cbc_offset} --demultiplex --prefix {params.prefix}/ > {log} 2>&1"
else:
    checkpoint parse_fastq:
        input: r1_in = config["r1"], r2_in = config["r2"], pbcpath = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"]), cell_barcodes = "results/metadata/{name}_cell_barcodes.txt".format(name=config["name"]), sample_map = "results/metadata/{name}_sample_map.yaml".format(name=config["name"]), readtype_map = "results/metadata/{name}_readtype_map.yaml".format(name=config["name"]), dt_structure = "results/metadata/{}_dt_structure.yaml".format(config["name"])
        output: directory(FASTQ_DIR)
        log: "results/logs/parse_fastq.log"
        benchmark: "results/benchmarks/parse_fastq.benchmark.txt"
        params: comp_threads = int(config["threads"]*0.2),
                proc_threads = config["threads"]-int(config["threads"]*0.2),
                cbc_offset = config['params']['parse_fq']['cbc_offset'],
                prefix = FASTQ_DIR
        shell: "echo Parse FASTQ && binaries/parse_fastq --read1 {input.r1_in} --read2 {input.r2_in} --cbcpath {input.cell_barcodes} --pbcpath {input.pbcpath} --readtype-structure {input.readtype_map} --dt-structure {input.dt_structure} --index-layout {config[index_layout]} --sample-structure {input.sample_map} --processing-threads {params.proc_threads} --compression-threads {params.comp_threads} --ts-sequence {config[ts_sequence]} --ts-pad {config[ts_pad]} --ts-cutoff {config[ts_cutoff]} --cbc-offset {params.cbc_offset} --demultiplex --prefix {params.prefix}/ > {log} 2>&1"

def get_samples(wildcards):
    checkpoint_output = checkpoints.parse_fastq.get(**wildcards).output[0]
    print("Samples to be processed:")
    for f in os.listdir(checkpoint_output):
        if f.endswith("_1.fq.gz")
            print(f.replace("_1.fq.gz", ""))
    return [f.replace("_1.fq.gz", "") for f in os.listdir(checkpoint_output) if f.endswith("_1.fq.gz")]

def parse_fq_dir(_wc):
    return checkpoints.parse_fastq.get().output[0]

rule trim_fastq:
    input:  r1=lambda wc: os.path.join(parse_fq_dir(wc), f"{wc.sample}_1.fq.gz"),
            r2=lambda wc: os.path.join(parse_fq_dir(wc), f"{wc.sample}_2.fq.gz")
    output: r1 = temp(os.path.join(FASTQ_DIR,"{sample}.trimmed.read1.fastq.gz")),
           r2 = temp(os.path.join(FASTQ_DIR,"{sample}.trimmed.read2.fastq.gz")),
           r1_short = temp(os.path.join(FASTQ_DIR,"{sample}.tooshort.read1.fastq.gz")),
           r2_short = temp(os.path.join(FASTQ_DIR,"{sample}.tooshort.read2.fastq.gz"))
    log: stdout = "results/logs/{sample}.trim_fastq.log",
         summary = "results/summaries/{sample}.cutadapt.json"
    benchmark: "results/benchmarks/{sample}.trim_fastq.benchmark.txt"
    conda: "../envs/full.yaml"
    threads: config["threads"]
    shell: "echo Trim FASTQ && cutadapt -j {threads} --json {log.summary} {config[params][cutadapt]} --too-short-output {output.r1_short} --too-short-paired-output {output.r2_short} -o {output.r1} -p {output.r2} {input.r1} {input.r2} > {log.stdout} 2>&1"

rule map_reads:
    input: r1 = os.path.join(FASTQ_DIR,"{sample}.trimmed.read1.fastq.gz"),
        r2 = os.path.join(FASTQ_DIR,"{sample}.trimmed.read2.fastq.gz")
    output: temp("results/intermediate/{sample}.trimmed.aligned.bam")
    log: stdout = "results/logs/map_reads.log",
        summary = "results/summaries/{sample}.hisat2.summary.txt"
    benchmark: "results/benchmarks/map_reads.benchmark.txt"
    params: splicesites = SPLICESITES, genomeref = GENOMEREF, cores_hisat = cores_hisat, cores_samtools = cores_samtools
    conda: "../envs/full.yaml"
    shell: "echo Map Reads && binaries/hisat-3n --new-summary --summary-file {log.summary} {config[params][hisat3n]} -p {params.cores_hisat} --known-splicesite-infile {params.splicesites} -x {params.genomeref} -1 {input.r1} -2 {input.r2} | samtools view -F 256 -b -@ {params.cores_samtools} -o {output} > {log.stdout} 2>&1"

rule split_bam_by_strand:
    input: "results/intermediate/{sample}.trimmed.aligned.bam"
    output: pstrand = temp("results/intermediate/{sample}.trimmed.aligned.pstrand.bam"),
            mstrand = temp("results/intermediate/{sample}.trimmed.aligned.mstrand.bam"),
            nostrand = temp("results/intermediate/{sample}.trimmed.aligned.nostrand.bam")
    log: "results/logs/split_bam_by_strand.log"
    benchmark: "results/benchmarks/split_bam_by_strand.benchmark.txt"
    conda: "../envs/full.yaml"
    shell: "echo Assign Genes && binaries/move_tags --input {input} > {log} 2>&1"

rule assign_genes_exon:
    input: pstrand = "results/intermediate/{sample}.trimmed.aligned.pstrand.bam",
            mstrand = "results/intermediate/{sample}.trimmed.aligned.mstrand.bam",
            nostrand = "results/intermediate/{sample}.trimmed.aligned.nostrand.bam"
    output: pstrand = temp("results/intermediate/{sample}.trimmed.aligned.pstrand.bam.featureCounts.bam"),
            mstrand = temp("results/intermediate/{sample}.trimmed.aligned.mstrand.bam.featureCounts.bam"),
            nostrand = temp("results/intermediate/{sample}.trimmed.aligned.nostrand.bam.featureCounts.bam")
    log: "results/logs/assign_genes.log"
    benchmark: "results/benchmarks/assign_genes.benchmark.txt"
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
    """

rule rename_tags_exon:
    input:  pstrand = "results/intermediate/{sample}.trimmed.aligned.pstrand.bam.featureCounts.bam",
            mstrand = "results/intermediate/{sample}.trimmed.aligned.mstrand.bam.featureCounts.bam",
            nostrand = "results/intermediate/{sample}.trimmed.aligned.nostrand.bam.featureCounts.bam"
    output: pstrand = temp("results/intermediate/{sample}.trimmed.aligned.pstrand.Exon.bam"),
            mstrand = temp("results/intermediate/{sample}.trimmed.aligned.mstrand.Exon.bam"),
            nostrand = temp("results/intermediate/{sample}.trimmed.aligned.nostrand.Exon.bam")
    log: "results/logs/rename_tags_exon.log"
    conda: "../envs/full.yaml"
    shell:"""
    binaries/rename_tags --input {input.pstrand} --output {output.pstrand} >> {log} 2>&1
    binaries/rename_tags --input {input.mstrand} --output {output.mstrand} >> {log} 2>&1
    binaries/rename_tags --input {input.nostrand} --output {output.nostrand} >> {log} 2>&1
    """

rule assign_genes_intron:
    input: pstrand = "results/intermediate/{sample}.trimmed.aligned.pstrand.Exon.bam",
            mstrand = "results/intermediate/{sample}.trimmed.aligned.mstrand.Exon.bam",
            nostrand = "results/intermediate/{sample}.trimmed.aligned.nostrand.Exon.bam"
    output: pstrand = temp("results/intermediate/{sample}.trimmed.aligned.pstrand.Exon.bam.featureCounts.bam"),
            mstrand = temp("results/intermediate/{sample}.trimmed.aligned.mstrand.Exon.bam.featureCounts.bam"),
            nostrand = temp("results/intermediate/{sample}.trimmed.aligned.nostrand.Exon.bam.featureCounts.bam")
    log: "results/logs/assign_genes.log"
    benchmark: "results/benchmarks/assign_genes.benchmark.txt"
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
    """
rule rename_tags_intron:
    input:  pstrand = "results/intermediate/{sample}.trimmed.aligned.pstrand.Exon.bam.featureCounts.bam",
            mstrand = "results/intermediate/{sample}.trimmed.aligned.mstrand.Exon.bam.featureCounts.bam",
            nostrand = "results/intermediate/{sample}.trimmed.aligned.nostrand.Exon.bam.featureCounts.bam"
    output: pstrand = temp("results/intermediate/{sample}.trimmed.aligned.pstrand.GeneTagged.bam"),
            mstrand = temp("results/intermediate/{sample}.trimmed.aligned.mstrand.GeneTagged.bam"),
            nostrand = temp("results/intermediate/{sample}.trimmed.aligned.nostrand.GeneTagged.bam")
    log: "results/logs/rename_tags_intron.log"
    conda: "../envs/full.yaml"
    shell:"""
    binaries/rename_tags --input {input.pstrand} --output {output.pstrand} --intron-mode >> {log} 2>&1
    binaries/rename_tags --input {input.mstrand} --output {output.mstrand} --intron-mode >> {log} 2>&1
    binaries/rename_tags --input {input.nostrand} --output {output.nostrand} --intron-mode >> {log} 2>&1
    """
rule concatenate_and_sort:
    input: pstrand = "results/intermediate/{sample}.trimmed.aligned.pstrand.GeneTagged.bam",
            mstrand = "results/intermediate/{sample}.trimmed.aligned.mstrand.GeneTagged.bam",
            nostrand = "results/intermediate/{sample}.trimmed.aligned.nostrand.GeneTagged.bam"
    log: "results/logs/concatenate_and_sort.log"
    benchmark: "results/benchmarks/concatenate_and_sort.benchmark.txt"
    output: temp("results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted.bam")
    conda: "../envs/full.yaml"
    threads: config["threads"]
    shell: "samtools cat {input.nostrand} {input.pstrand} {input.mstrand} | samtools sort -m 1000M -@ {threads} -T results/.tmp_bgab/sorttmp. -o {output} >> {log} 2>&1"

rule first_index:
    input: "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted.bam"
    output: temp("results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted.bam.bai")
    threads: config["threads"]
    log: "results/logs/first_index.log"
    shell: "samtools index -@ {config[threads]} {input}"

if config["reverse"]:
    rule reconstruct:
        input: bam = "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted.bam",
            bai = "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted.bam.bai",
            sample_map = "results/metadata/{sample}_sample_map.yaml"
        output: temp("results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted.reconstructed.bam")
        params: gtffile = GTFFILE, cores_hisat = cores_hisat, cores_samtools = cores_samtools
        log: "results/logs/reconstruct.log"
        threads: min(config["threads"], 64)
        benchmark: "results/benchmarks/reconstruction.benchmark.txt"
        shell:"""
        echo Reconstruct Molecules
        mkdir -p results/tmp/
        binaries/basic_reconstruction --input {input.bam} --output {output} --gtf {params.gtffile}.gff3 --sample-map {input.sample_map} --threads {threads} --gene-identifier {config[gff_gene_identifier]} --reverse --bulk | samtools view -F 256 -b -@ {params.cores_samtools} -o {output} > {log} 2>&1
        """
else:
    rule reconstruct:
        input: bam = "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted.bam",
            bai = "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted.bam.bai",
            sample_map = "results/metadata/{sample}_sample_map.yaml"
        output: temp("results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted.reconstructed.bam")
        params: gtffile = GTFFILE, cores_hisat = cores_hisat, cores_samtools = cores_samtools
        log: "results/logs/reconstruct.log"
        threads: min(config["threads"], 64)
        benchmark: "results/benchmarks/reconstruction.benchmark.txt"
        shell:"""
        echo Reconstruct Molecules
        mkdir -p results/tmp/
        binaries/basic_reconstruction --input {input.bam} --output {output} --gtf {params.gtffile}.gff3 --sample-map {input.sample_map} --threads {threads} --gene-identifier {config[gff_gene_identifier]} --bulk  | samtools view -F 256 -b -@ {params.cores_samtools} -o {output} > {log} 2>&1
        """

rule sort_reconstructed:
    input: "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted.reconstructed.bam"
    output: "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted.reconstructed.sorted.bam"
    threads: config["threads"]
    log: "results/logs/sort_reconstructed.log"
    params:
        extra="-m 1000M",
    conda: "../envs/full.yaml"
    shell: "samtools sort -@ {threads} -o {output} {params.extra} -T results/.tmp_bgab/sorttmp. {input} &> {log}"

rule index_reconstructed:
    input: "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted.reconstructed.sorted.bam"
    output: "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted.reconstructed.sorted.bam.bai"
    threads: config["threads"]
    log: "results/logs/index_reconstructed.log"
    shell: "samtools index -@ {config[threads]} {input}"

rule stitch_reconstruction:
    input: bam =  "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted.reconstructed.sorted.bam",
        bai = "results/intermediate/{sample}.reads.aligned_trimmed_genetagged_sorted.reconstructed.sorted.bam.bai"
    output: temp("results/intermediate/{sample}.stitched.bam")
    threads: config["threads"]
    params: gtffile = GTFFILE
    log: "results/logs/stitch_reconstruction.log"
    conda: "../envs/full.yaml"
    shell: "echo Stitch Molecules && python3 workflow/scripts/stitcher.py --input {input.bam} --output {output} --gtf {params.gtffile}.gff3 --threads {threads} --cell-tag CB --UMI-tag RM --gene-identifier {config[gff_gene_identifier]} >> {log} 2>&1"

rule sorted_stitched:
    input: "results/intermediate/{sample}.stitched.bam"
    output: temp("results/intermediate/{sample}.stitched.sorted.bam")
    threads: config["threads"]
    log: "results/logs/sort_stitched.log"
    params:
        extra="-m 1000M",
    conda: "../envs/full.yaml"
    shell: "samtools sort -@ {threads} -o {output} {params.extra} -T results/.tmp_bgab/sorttmp. {input} &> {log}"

rule index_stitched:
    input: "results/intermediate/{sample}.stitched.sorted.bam"
    output: temp("results/intermediate/{sample}.stitched.sorted.bam.bai")
    threads: config["threads"]
    log: "results/logs/index_stitched.log"
    shell: "samtools index -@ {config[threads]} {input}"

rule make_molecule_bams:
    input: bam = "results/intermediate/{sample}.stitched.sorted.bam", bai = "results/intermediate/{sample}.stitched.sorted.bam.bai"
    output: molecules_out = temp("results/intermediate/{sample}.stitched.molecules.bam")
    conda: "../envs/full.yaml"
    shell:"python3 workflow/scripts/filter_stitched_bam.py --input {input.bam} --molecules-out {output.molecules_out}"

rule sort_molecule_bams:
    input: molecules_out = "results/intermediate/{sample}.stitched.molecules.bam"
    output: molecules_out = "results/{sample}.stitched.molecules.sorted.bam"
    log: "results/logs/sort_molecules.log"
    conda: "../envs/full.yaml"
    shell:"""
    echo Creating final output file
    samtools sort -m 1000M -@ {config[threads]} -T results/.tmp_bgab/sorttmp. -o {output.molecules_out} {input.molecules_out} > {log} 2>&1
    samtools index -@ {config[threads]} {output.molecules_out}
    """