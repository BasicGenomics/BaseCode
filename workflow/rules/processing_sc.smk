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
    shell: "echo Get cell barcodes && binaries/pipspeak_basecode -c {params.barcode_cfg} -i {input.fastq_r1} -I {input.fastq_r2} -l -u 0 -p results/barcodes/{params.prefix} > {log} 2>&1"

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

if config["i1"] != "" and config["i2"] != "":
    checkpoint parse_fastq:
        input: r1_in = config["r1"], r2_in = config["r2"], i1_in = config["i1"], i2_in = config["i2"], pbcpath = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"]), cell_barcodes = "results/metadata/{name}_whitelist.txt".format(name=config["name"]), sample_map = "results/metadata/{name}_sample_map.yaml".format(name=config["name"]), readtype_map = "results/metadata/{name}_readtype_map.yaml".format(name=config["name"]), dt_structure = "results/metadata/{}_dt_structure.yaml".format(config["name"])
        output: directory("results/intermediate/demultiplexed_fastq")
        log: "results/logs/parse_fastq.log"
        benchmark: "results/benchmarks/parse_fastq.benchmark.txt"
        params: comp_threads = int(config["threads"]*0.2),
                proc_threads = config["threads"]-int(config["threads"]*0.2),
                cbc_offset = config['params']['cbc_offset']
        shell: "echo Parse FASTQ && binaries/parse_fastq --read1 {input.r1_in} --read2 {input.r2_in} --index1 {input.i1_in} --index2 {input.i2_in} --cbcpath {input.cell_barcodes} --pbcpath {input.pbcpath} --readtype-structure {input.readtype_map} --dt-structure {input.dt_structure} --index-layout {config[index_layout]} --sample-structure {input.sample_map} --processing-threads {params.proc_threads} --compression-threads {params.comp_threads} --ts-sequence {config[ts_sequence]} --ts-pad {config[ts_pad]} --ts-cutoff {config[ts_cutoff]} --cbc-offset {params.cbc_offset} --demultiplex --prefix results/intermediate/demultiplexed_fastq/ > {log} 2>&1"
else:
    checkpoint parse_fastq:
        input: r1_in = config["r1"], r2_in = config["r2"], pbcpath = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"]), cell_barcodes = "results/metadata/{name}_whitelist.txt".format(name=config["name"]), sample_map = "results/metadata/{name}_sample_map.yaml".format(name=config["name"]), readtype_map = "results/metadata/{name}_readtype_map.yaml".format(name=config["name"]), dt_structure = config["dt_structure"]
        output: directory("results/intermediate/demultiplexed_fastq")
        log: "results/logs/parse_fastq.log"
        benchmark: "results/benchmarks/parse_fastq.benchmark.txt"
        params: comp_threads = int(config["threads"]*0.2),
                proc_threads = config["threads"]-int(config["threads"]*0.2),
                cbc_offset = config['params']['cbc_offset']
        shell: "echo Parse FASTQ && binaries/parse_fastq --read1 {input.r1_in} --read2 {input.r2_in} --cbcpath {input.cell_barcodes} --pbcpath {input.pbcpath} --readtype-structure {input.readtype_map} --dt-structure {input.dt_structure} --index-layout {config[index_layout]} --sample-structure {input.sample_map} --processing-threads {params.proc_threads} --compression-threads {params.comp_threads} --ts-sequence {config[ts_sequence]} --ts-pad {config[ts_pad]} --ts-cutoff {config[ts_cutoff]} --cbc-offset {params.cbc_offset} --demultiplex --prefix results/intermediate/demultiplexed_fastq/ > {log} 2>&1"

def get_samples(wildcards):
    checkpoint_output = checkpoints.parse_fastq.get(**wildcards).output[0]
    print(checkpoint_output)
    return [f.replace("_1.fq.gz", "") for f in os.listdir(checkpoint_output) if f.endswith("_1.fq.gz")]

rule trim_fastq:
    input: r1 = "results/intermediate/{name}.read1.fastq.gz".format(name=config["name"]),
           r2 = "results/intermediate/{name}.read2.fastq.gz".format(name=config["name"])
    output: r1 = "results/intermediate/{name}.trimmed.read1.fastq.gz".format(name=config["name"]),
           r2 = "results/intermediate/{name}.trimmed.read2.fastq.gz".format(name=config["name"]),
           r1_short = "results/intermediate/{name}.tooshort.read1.fastq.gz".format(name=config["name"]),
           r2_short = "results/intermediate/{name}.tooshort.read2.fastq.gz".format(name=config["name"])
    log: stdout = "results/logs/trim_fastq.log",
         summary = "results/summaries/{name}.cutadapt.json".format(name=config["name"])
    benchmark: "results/benchmarks/trim_fastq.benchmark.txt"
    conda: "../envs/full.yaml"
    threads: config["threads"]
    shell: "echo Trim FASTQ && cutadapt -j {threads} --json {log.summary} {config[params][cutadapt]} --too-short-output {output.r1_short} --too-short-paired-output {output.r2_short} -o {output.r1} -p {output.r2} {input.r1} {input.r2} > {log.stdout} 2>&1"

rule map_reads:
    input: r1 = "results/intermediate/{name}.trimmed.read1.fastq.gz".format(name=config["name"]),
        r2 = "results/intermediate/{name}.trimmed.read2.fastq.gz".format(name=config["name"])
    output: "results/intermediate/{name}.trimmed.aligned.bam".format(name=config["name"])
    log: stdout = "results/logs/map_reads.log",
        summary = "results/summaries/{name}.hisat2.summary.txt".format(name=config["name"])
    benchmark: "results/benchmarks/map_reads.benchmark.txt"
    params: splicesites = SPLICESITES, genomeref = GENOMEREF, cores_hisat = cores_hisat, cores_samtools = cores_samtools
    conda: "../envs/full.yaml"
    shell: "echo Map Reads && binaries/hisat-3n --new-summary --summary-file {log.summary} {config[params][hisat3n]} -p {params.cores_hisat} --known-splicesite-infile {params.splicesites} -x {params.genomeref} -1 {input.r1} -2 {input.r2} | samtools view -F 256 -b -@ {params.cores_samtools} -o {output} > {log.stdout} 2>&1"

rule split_bam_by_strand:
    input: "results/intermediate/{name}.trimmed.aligned.bam".format(name=config["name"])
    output: pstrand = "results/intermediate/{name}.trimmed.aligned.pstrand.bam".format(name=config["name"]),
            mstrand = "results/intermediate/{name}.trimmed.aligned.mstrand.bam".format(name=config["name"]),
            nostrand = "results/intermediate/{name}.trimmed.aligned.nostrand.bam".format(name=config["name"])
    log: "results/logs/split_bam_by_strand.log"
    benchmark: "results/benchmarks/split_bam_by_strand.benchmark.txt"
    conda: "../envs/full.yaml"
    shell: "echo Assign Genes && binaries/move_tags --input {input} > {log} 2>&1"

rule assign_genes_exon:
    input: pstrand = "results/intermediate/{name}.trimmed.aligned.pstrand.bam".format(name=config["name"]),
            mstrand = "results/intermediate/{name}.trimmed.aligned.mstrand.bam".format(name=config["name"]),
            nostrand = "results/intermediate/{name}.trimmed.aligned.nostrand.bam".format(name=config["name"])
    output: pstrand = temp("results/intermediate/{name}.trimmed.aligned.pstrand.bam.featureCounts.bam".format(name=config["name"])),
            mstrand = temp("results/intermediate/{name}.trimmed.aligned.mstrand.bam.featureCounts.bam".format(name=config["name"])),
            nostrand = temp("results/intermediate/{name}.trimmed.aligned.nostrand.bam.featureCounts.bam".format(name=config["name"]))
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
    input:  pstrand = "results/intermediate/{name}.trimmed.aligned.pstrand.bam.featureCounts.bam".format(name=config["name"]),
            mstrand = "results/intermediate/{name}.trimmed.aligned.mstrand.bam.featureCounts.bam".format(name=config["name"]),
            nostrand = "results/intermediate/{name}.trimmed.aligned.nostrand.bam.featureCounts.bam".format(name=config["name"])
    output: pstrand = temp("results/intermediate/{name}.trimmed.aligned.pstrand.Exon.bam".format(name=config["name"])),
            mstrand = temp("results/intermediate/{name}.trimmed.aligned.mstrand.Exon.bam".format(name=config["name"])),
            nostrand = temp("results/intermediate/{name}.trimmed.aligned.nostrand.Exon.bam".format(name=config["name"]))
    log: "results/logs/rename_tags_exon.log"
    conda: "../envs/full.yaml"
    shell:"""
    binaries/rename_tags --input {input.pstrand} --output {output.pstrand} >> {log} 2>&1
    binaries/rename_tags --input {input.mstrand} --output {output.mstrand} >> {log} 2>&1
    binaries/rename_tags --input {input.nostrand} --output {output.nostrand} >> {log} 2>&1
    """

rule assign_genes_intron:
    input: pstrand = "results/intermediate/{name}.trimmed.aligned.pstrand.Exon.bam".format(name=config["name"]),
            mstrand = "results/intermediate/{name}.trimmed.aligned.mstrand.Exon.bam".format(name=config["name"]),
            nostrand = "results/intermediate/{name}.trimmed.aligned.nostrand.Exon.bam".format(name=config["name"])
    output: pstrand = temp("results/intermediate/{name}.trimmed.aligned.pstrand.Exon.bam.featureCounts.bam".format(name=config["name"])),
            mstrand = temp("results/intermediate/{name}.trimmed.aligned.mstrand.Exon.bam.featureCounts.bam".format(name=config["name"])),
            nostrand = temp("results/intermediate/{name}.trimmed.aligned.nostrand.Exon.bam.featureCounts.bam".format(name=config["name"]))
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
    input:  pstrand = "results/intermediate/{name}.trimmed.aligned.pstrand.Exon.bam.featureCounts.bam".format(name=config["name"]),
            mstrand = "results/intermediate/{name}.trimmed.aligned.mstrand.Exon.bam.featureCounts.bam".format(name=config["name"]),
            nostrand = "results/intermediate/{name}.trimmed.aligned.nostrand.Exon.bam.featureCounts.bam".format(name=config["name"])
    output: pstrand = temp("results/intermediate/{name}.trimmed.aligned.pstrand.GeneTagged.bam".format(name=config["name"])),
            mstrand = temp("results/intermediate/{name}.trimmed.aligned.mstrand.GeneTagged.bam".format(name=config["name"])),
            nostrand = temp("results/intermediate/{name}.trimmed.aligned.nostrand.GeneTagged.bam".format(name=config["name"]))
    log: "results/logs/rename_tags_intron.log"
    conda: "../envs/full.yaml"
    shell:"""
    binaries/rename_tags --input {input.pstrand} --output {output.pstrand} --intron-mode >> {log} 2>&1
    binaries/rename_tags --input {input.mstrand} --output {output.mstrand} --intron-mode >> {log} 2>&1
    binaries/rename_tags --input {input.nostrand} --output {output.nostrand} --intron-mode >> {log} 2>&1
    """

rule concatenate_and_sort:
    input: pstrand = "results/intermediate/{name}.trimmed.aligned.pstrand.GeneTagged.bam".format(name=config["name"]),
            mstrand = "results/intermediate/{name}.trimmed.aligned.mstrand.GeneTagged.bam".format(name=config["name"]),
            nostrand = "results/intermediate/{name}.trimmed.aligned.nostrand.GeneTagged.bam".format(name=config["name"])
    log: "results/logs/concatenate_and_sort.log"
    benchmark: "results/benchmarks/concatenate_and_sort.benchmark.txt"
    output: temp("results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.bam".format(name=config["name"]))
    conda: "../envs/full.yaml"
    threads: config["threads"]
    shell: "samtools cat {input.nostrand} {input.pstrand} {input.mstrand} | samtools sort -m 1000M -@ {threads} -T results/.tmp_bgab/sorttmp. -o {output} >> {log} 2>&1"

rule first_index:
    input: "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.bam".format(name=config["name"])
    output: temp("results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.bam.bai".format(name=config["name"]))
    threads: config["threads"]
    log: "results/logs/first_index.log"
    shell: "samtools index -@ {config[threads]} {input}"

rule error_correct_umis:
    input: bam = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.bam".format(name=config["name"]), 
        bai = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.bam.bai".format(name=config["name"])
    output: "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted_umicorrected.bam".format(name=config["name"])
    log: "results/logs/error_correct_umis.log"
    benchmark: "results/benchmarks/error_correct_umis.benchmark.txt"
    params: gtffile = GTFFILE
    conda: "../envs/full.yaml"
    shell: "python3 workflow/scripts/error_correct_umis.py {input.bam} {params.gtffile}.gff3 {config[threads]} &> {log}"

rule second_index:
    input: "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted_umicorrected.bam".format(name=config["name"])
    output: "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted_umicorrected.bam.bai".format(name=config["name"])
    params:
        extra="",
    threads: config["threads"]
    wrapper:
        "v3.3.3/bio/samtools/index"

if config["reverse"]:
    rule reconstruct:
        input: bam = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted_umicorrected.bam".format(name=config["name"]),
            bai = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted_umicorrected.bam.bai".format(name=config["name"]),
            sample_map = "results/metadata/{name}_sample_map.yaml".format(name=config["name"])
        output: temp("results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted_umicorrected.reconstructed.bam".format(name=config["name"]))
        params: gtffile = GTFFILE, cores_hisat = cores_hisat, cores_samtools = cores_samtools
        log: "results/logs/reconstruct.log"
        threads: min(config["threads"], 64)
        benchmark: "results/benchmarks/reconstruction.benchmark.txt"
        shell: """
        echo Reconstruct Molecules
        mkdir -p results/tmp/
        binaries/basic_reconstruction --input {input.bam} --output {output} --gtf {params.gtffile}.gff3 --sample-map {input.sample_map} --threads {threads} --gene-identifier {config[gff_gene_identifier]} --reverse | samtools view -F 256 -b -@ {params.cores_samtools} -o {output} > {log} 2>&1
        """
else:
    rule reconstruct:
        input: bam = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted_umicorrected.bam".format(name=config["name"]),
            bai = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted_umicorrected.bam.bai".format(name=config["name"]),
            sample_map = "results/metadata/{name}_sample_map.yaml".format(name=config["name"])
        output: temp("results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted_umicorrected.reconstructed.bam".format(name=config["name"]))
        params: gtffile = GTFFILE, cores_hisat = cores_hisat, cores_samtools = cores_samtools
        log: "results/logs/reconstruct.log"
        threads: min(config["threads"], 64)
        benchmark: "results/benchmarks/reconstruction.benchmark.txt"
        shell: """
        echo Reconstruct Molecules
        mkdir -p results/tmp/
        binaries/basic_reconstruction --input {input.bam} --output {output} --gtf {params.gtffile}.gff3 --sample-map {input.sample_map} --threads {threads} --gene-identifier {config[gff_gene_identifier]}  | samtools view -F 256 -b -@ {params.cores_samtools} -o {output} > {log} 2>&1
        """

rule sort_reconstructed:
    input: "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted_umicorrected.reconstructed.bam".format(name=config["name"])
    output: "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted_umicorrected.reconstructed.sorted.bam".format(name=config["name"])
    threads: config["threads"]
    log: "results/logs/sort_reconstructed.log"
    params:
        extra="-m 1000M",
    conda: "../envs/full.yaml"
    shell: "samtools sort -@ {threads} -o {output} {params.extra} -T . {input} &> {log}"

rule index_reconstructed:
    input: "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted_umicorrected.reconstructed.sorted.bam".format(name=config["name"])
    output: "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted_umicorrected.reconstructed.sorted.bam.bai".format(name=config["name"])
    log: "results/logs/index_reconstructed.log"
    params:
        extra="",
    threads: config["threads"]
    wrapper:
        "v3.3.3/bio/samtools/index"

rule stitch_reconstruction:
        input: bam =  "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted_umicorrected.reconstructed.sorted.bam".format(name=config["name"]),
            bai = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted_umicorrected.reconstructed.sorted.bam.bai".format(name=config["name"])
        output: temp("results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted_umicorrected.stitched.bam".format(name=config["name"]))
        threads: config["threads"]
        params: gtffile = GTFFILE
        log: "results/logs/stitch_reconstruction.log"
        conda: "../envs/full.yaml"
        shell: """
        echo Stitch Molecules 
        python3 workflow/scripts/stitcher.py --input {input.bam} --output {output} --gtf {params.gtffile}.gff3 --threads {threads} --cell-tag CB --UMI-tag RM --gene-identifier {config[gff_gene_identifier]} --only-molecules >> {log} 2>&1
        """

rule sorted_stitched:
    input: "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted_umicorrected.stitched.bam".format(name=config["name"])
    output: "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted_umicorrected.stitched.sorted.bam".format(name=config["name"])
    threads: config["threads"]
    log: "results/logs/sort_stitched.log"
    params:
        extra="-m 1000M",
    conda: "../envs/full.yaml"
    shell: "samtools sort -@ {threads} -o {output} {params.extra} -T . {input} &> {log}"

rule index_stitched:
    input: "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted_umicorrected.stitched.sorted.bam".format(name=config["name"])
    output: "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted_umicorrected.stitched.sorted.bam.bai".format(name=config["name"])
    params:
        extra="",
    threads: config["threads"]
    wrapper:
        "v3.3.3/bio/samtools/index"

rule make_molecule_bams:
    input: bam = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted_umicorrected.stitched.sorted.bam".format(name=config["name"]), bai = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted_umicorrected.stitched.sorted.bam.bai".format(name=config["name"])
    output: molecules_out = temp("results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted_umicorrected.stitched.molecules.bam".format(name=config["name"])),
    conda: "../envs/full.yaml"
    shell:"python3 workflow/scripts/filter_stitched_bam.py --input {input.bam} --molecules-out {output.molecules_out}"

rule sort_molecule_bams:
    input: molecules_out = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted_umicorrected.stitched.molecules.bam".format(name=config["name"])
    output: molecules_out = "results/{name}.reads.stitched.molecules.sorted.bam".format(name=config["name"])
    log: "results/logs/sort_molecules.log"
    conda: "../envs/full.yaml"
    shell:"""
    echo Creating final output file
    samtools sort -m 1000M -@ {config[threads]} -T results/.tmp_bgab/sorttmp. -o {output.molecules_out} {input.molecules_out} &> {log}
    samtools index -@ {config[threads]} {output.molecules_out}
    """