rule log_python_version:
    output: "results/logs/python_version.log"
    shell: "python3 --version > {output}"

rule log_version:
    output: "results/logs/basecode_version.log"
    shell: "cat VERSION > {output}"

if config["ignore_none"]:
    if config["i1"] != "" and config["i2"] != "":
        rule make_barcode_files:
            input: samplesheet = config["samplesheet"], 
                   i1 = config["i1"],
                   i2 = config["i2"],
                   index_sequences = "workflow/resources/index_sequences.yaml"
            output: barcodes = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"]),
                    cell_barcodes = "results/metadata/{name}_cell_barcodes.txt".format(name=config["name"]),
                    sample_map = "results/metadata/{name}_sample_map.yaml".format(name=config["name"]),
                    readtype_map = "results/metadata/{name}_readtype_map.yaml".format(name=config["name"]),
                    samplesheet = "results/metadata/{name}_samplesheet.csv".format(name=config["name"]),
                    done = "results/dones/{name}_make_barcode_files.done".format(name=config["name"])
            log: "results/logs/make_barcode_files.log"
            shell: """
            echo Step 1/7 Process Sample Sheet
            python3 workflow/scripts/make_sample_files.py -s {input.samplesheet} --fastq {input.i1} {input.i2} --index-sequences {input.index_sequences} --sample-barcodes {output.barcodes} --cell-barcodes {output.cell_barcodes} --sample-map {output.sample_map} --readtype-map {output.readtype_map} --samplesheet-out {output.samplesheet} --ignore-none > {log} 2>&1
            touch {output.done}
            """
    else:
        rule make_barcode_files:
            input: samplesheet = config["samplesheet"],
                   r2 = config["r2"],
                   index_sequences = "workflow/resources/index_sequences.yaml"
            output: barcodes = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"]),
                    cell_barcodes = "results/metadata/{name}_cell_barcodes.txt".format(name=config["name"]),
                    sample_map = "results/metadata/{name}_sample_map.yaml".format(name=config["name"]),
                    readtype_map = "results/metadata/{name}_readtype_map.yaml".format(name=config["name"]),
                    samplesheet = "results/metadata/{name}_samplesheet.csv".format(name=config["name"]),
                    done = "results/dones/{name}_make_barcode_files.done".format(name=config["name"])
            log: "results/logs/make_barcode_files.log"
            shell: """
            echo Step 1/7 Process Sample Sheet
            python3 workflow/scripts/make_sample_files.py -s {input.samplesheet} --fastq {input.r2} --index-sequences {input.index_sequences} --sample-barcodes {output.barcodes} --cell-barcodes {output.cell_barcodes} --sample-map {output.sample_map} --readtype-map {output.readtype_map} --samplesheet-out {output.samplesheet} --ignore-none > {log} 2>&1
            touch {output.done}
            """
else:
    if config["i1"] != "" and config["i2"] != "":
        rule make_barcode_files:
            input: samplesheet = config["samplesheet"],
                   i1 = config["i1"],
                   i2 = config["i2"],
                   index_sequences = "workflow/resources/index_sequences.yaml"
            output: barcodes = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"]),
                    cell_barcodes = "results/metadata/{name}_cell_barcodes.txt".format(name=config["name"]),
                    sample_map = "results/metadata/{name}_sample_map.yaml".format(name=config["name"]),
                    readtype_map = "results/metadata/{name}_readtype_map.yaml".format(name=config["name"]),
                    samplesheet = "results/metadata/{name}_samplesheet.csv".format(name=config["name"]),
                    done = "results/dones/{name}_make_barcode_files.done".format(name=config["name"])
            log: "results/logs/make_barcode_files.log"
            shell: """
            echo Step 1/7 Process Sample Sheet
            python3 workflow/scripts/make_sample_files.py -s {input.samplesheet} --fastq {input.i1} {input.i2} --index-sequences {input.index_sequences} --sample-barcodes {output.barcodes} --cell-barcodes {output.cell_barcodes} --sample-map {output.sample_map} --readtype-map {output.readtype_map} --samplesheet-out {output.samplesheet} > {log} 2>&1
            touch {output.done}
            """
    else:
        rule make_barcode_files:
            input: samplesheet = config["samplesheet"],
                   r2 = config["r2"],
                   index_sequences = "workflow/resources/index_sequences.yaml"
            output: barcodes = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"]),
                    cell_barcodes = "results/metadata/{name}_cell_barcodes.txt".format(name=config["name"]),
                    sample_map = "results/metadata/{name}_sample_map.yaml".format(name=config["name"]),
                    readtype_map = "results/metadata/{name}_readtype_map.yaml".format(name=config["name"]),
                    samplesheet = "results/metadata/{name}_samplesheet.csv".format(name=config["name"]),
                    done = "results/dones/{name}_make_barcode_files.done".format(name=config["name"])
            log: "results/logs/make_barcode_files.log"
            shell: """
            echo Step 1/7 Process Sample Sheet
            python3 workflow/scripts/make_sample_files.py -s {input.samplesheet} --fastq {input.r2} --index-sequences {input.index_sequences} --sample-barcodes {output.barcodes} --cell-barcodes {output.cell_barcodes} --sample-map {output.sample_map} --readtype-map {output.readtype_map} --samplesheet-out {output.samplesheet} > {log} 2>&1
            touch {output.done}
            """

if config["i1"] != "" and config["i2"] != "":
    rule parse_fastq:
        input: r1 = config["r1"],
               r2 = config["r2"],
               i1 = config["i1"],
               i2 = config["i2"], 
               pbcpath = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"]),
               cell_barcodes = "results/metadata/{name}_cell_barcodes.txt".format(name=config["name"]),
               sample_map = "results/metadata/{name}_sample_map.yaml".format(name=config["name"]),
               readtype_map = "results/metadata/{name}_readtype_map.yaml".format(name=config["name"]),
               dt_structure = config["dt_structure"]
        output: r1 = temp("results/intermediate/{name}.read1.fastq.gz".format(name=config["name"])),
                r2 = temp("results/intermediate/{name}.read2.fastq.gz".format(name=config["name"])),
                done = "results/dones/{name}_parse_fastq.done".format(name=config["name"])
        log: "results/logs/parse_fastq.log"
        benchmark: "results/benchmarks/parse_fastq.benchmark.txt"
        params: comp_threads = int(config["threads"]*0.2),
                proc_threads = config["threads"]-int(config["threads"]*0.2)
        shell: """
        echo Step 2/7 Parse FASTQ
        binaries/parse_fastq --read1 {input.r1} --read2 {input.r2} --index1 {input.i1} --index2 {input.i2} --r1-out {output.r1} --r2-out {output.r2} --cbcpath {input.cell_barcodes} --pbcpath {input.pbcpath} --readtype-structure {input.readtype_map} --dt-structure {input.dt_structure} --index-layout {config[index_layout]} --sample-structure {input.sample_map} --processing-threads {params.proc_threads} --compression-threads {params.comp_threads} --ts-sequence {config[ts_sequence]} --ts-pad {config[ts_pad]} --ts-cutoff {config[ts_cutoff]} > {log} 2>&1
        touch {output.done}
        """
else:
    rule parse_fastq:
        input: r1 = config["r1"],
               r2 = config["r2"], 
               pbcpath = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"]),
               cell_barcodes = "results/metadata/{name}_cell_barcodes.txt".format(name=config["name"]),
               sample_map = "results/metadata/{name}_sample_map.yaml".format(name=config["name"]),
               readtype_map = "results/metadata/{name}_readtype_map.yaml".format(name=config["name"]),
               dt_structure = config["dt_structure"]
        output: r1 = temp("results/intermediate/{name}.read1.fastq.gz".format(name=config["name"])),
                r2 = temp("results/intermediate/{name}.read2.fastq.gz".format(name=config["name"])),
                done = "results/dones/{name}_parse_fastq.done".format(name=config["name"])
        log: "results/logs/parse_fastq.log"
        benchmark: "results/benchmarks/parse_fastq.benchmark.txt"
        params: comp_threads = int(config["threads"]*0.2),
                proc_threads = config["threads"]-int(config["threads"]*0.2)
        shell: """
        echo Step 2/7 Parse FASTQ
        binaries/parse_fastq --read1 {input.r1} --read2 {input.r2} --r1-out {output.r1} --r2-out {output.r2} --cbcpath {input.cell_barcodes} --pbcpath {input.pbcpath} --readtype-structure {input.readtype_map} --dt-structure {input.dt_structure} --index-layout {config[index_layout]} --sample-structure {input.sample_map} --processing-threads {params.proc_threads} --compression-threads {params.comp_threads} --ts-sequence {config[ts_sequence]} --ts-pad {config[ts_pad]} --ts-cutoff {config[ts_cutoff]} > {log} 2>&1
        touch {output.done}
        """

rule trim_fastq:
    input: r1 = "results/intermediate/{name}.read1.fastq.gz".format(name=config["name"]),
           r2 = "results/intermediate/{name}.read2.fastq.gz".format(name=config["name"])
    output: r1 = temp("results/intermediate/{name}.trimmed.read1.fastq.gz".format(name=config["name"])),
            r2 = temp("results/intermediate/{name}.trimmed.read2.fastq.gz".format(name=config["name"])),
            r1_short = temp("results/intermediate/{name}.tooshort.read1.fastq.gz".format(name=config["name"])),
            r2_short = temp("results/intermediate/{name}.tooshort.read2.fastq.gz".format(name=config["name"])),
            done = "results/dones/{name}_trim_fastq.done".format(name=config["name"])
    log: log = "results/logs/trim_fastq.log",
         summary = "results/summaries/{name}.cutadapt.json".format(name=config["name"])
    benchmark: "results/benchmarks/trim_fastq.benchmark.txt"
    threads: config["threads"]
    shell: """
    echo Step 3/7 Trim FASTQ
    cutadapt -j {threads} --json {log.summary} {config[params][cutadapt]} --too-short-output {output.r1_short} --too-short-paired-output {output.r2_short} -o {output.r1} -p {output.r2} {input.r1} {input.r2} > {log.log} 2>&1
    touch {output.done}
    """

rule map_reads:
    input: r1 = "results/intermediate/{name}.trimmed.read1.fastq.gz".format(name=config["name"]),
           r2 = "results/intermediate/{name}.trimmed.read2.fastq.gz".format(name=config["name"])
    output: bam = temp("results/intermediate/{name}.trimmed.aligned.bam".format(name=config["name"])),
            done = "results/dones/{name}_map_reads.done".format(name=config["name"])
    log: log = "results/logs/map_reads.log",
         summary = "results/summaries/{name}.hisat2.summary.txt".format(name=config["name"])
    benchmark: "results/benchmarks/map_reads.benchmark.txt"
    params: splicesites = SPLICESITES,
            genomeref = GENOMEREF,
            cores_hisat = cores_hisat,
            cores_samtools = cores_samtools
    shell: """
    echo Step 4/7 Map Reads
    binaries/hisat-3n --new-summary --summary-file {log.summary} {config[params][hisat3n]} -p {params.cores_hisat} --known-splicesite-infile {params.splicesites} -x {params.genomeref} -1 {input.r1} -2 {input.r2} | samtools view -F 256 -b -@ {params.cores_samtools} -o {output.bam} > {log.log} 2>&1
    touch {output.done}
    """

rule split_bam_by_strand:
    input: "results/intermediate/{name}.trimmed.aligned.bam".format(name=config["name"])
    output: pstrand = temp("results/intermediate/{name}.trimmed.aligned.pstrand.bam".format(name=config["name"])),
            mstrand = temp("results/intermediate/{name}.trimmed.aligned.mstrand.bam".format(name=config["name"])),
            nostrand = temp("results/intermediate/{name}.trimmed.aligned.nostrand.bam".format(name=config["name"])),
            done = "results/dones/{name}_split_bam_by_strand.done".format(name=config["name"])
    log: "results/logs/split_bam_by_strand.log"
    benchmark: "results/benchmarks/split_bam_by_strand.benchmark.txt"
    threads: config["threads"]
    shell: """
    echo Step 5/7 Assign Genes
    binaries/move_tags --input {input} --threads {threads} > {log} 2>&1
    touch {output.done}
    """

rule assign_genes_exon:
    input: pstrand = "results/intermediate/{name}.trimmed.aligned.pstrand.bam".format(name=config["name"]),
           mstrand = "results/intermediate/{name}.trimmed.aligned.mstrand.bam".format(name=config["name"]),
           nostrand = "results/intermediate/{name}.trimmed.aligned.nostrand.bam".format(name=config["name"])
    output: pstrand = temp("results/intermediate/{name}.trimmed.aligned.pstrand.bam.featureCounts.bam".format(name=config["name"])),
            mstrand = temp("results/intermediate/{name}.trimmed.aligned.mstrand.bam.featureCounts.bam".format(name=config["name"])),
            nostrand = temp("results/intermediate/{name}.trimmed.aligned.nostrand.bam.featureCounts.bam".format(name=config["name"])),
            done = "results/dones/{name}_assign_genes_exon.done".format(name=config["name"])
    log: "results/logs/assign_genes_exon.log"
    benchmark: "results/benchmarks/assign_genes_exon.benchmark.txt"
    params: gff = "{}.gff3".format(GFF),
            gff_positive = "{}.positive.gff3".format(GFF),
            gff_negative = "{}.negative.gff3".format(GFF)
    threads: min(config["threads"], 64)
    shell:"""
    binaries/featureCounts -t exon --primary -g {config[gff_gene_identifier]} -T {threads} -R BAM -p --countReadPairs -O -M --largestOverlap --fracOverlap 0.1 -a {params.gff_positive} -o results/intermediate/positive.exon.FeatureCounts.txt {input.pstrand} >> {log} 2>&1
    binaries/featureCounts -t exon --primary -g {config[gff_gene_identifier]} -T {threads} -R BAM -p --countReadPairs -O -M --largestOverlap --fracOverlap 0.1 -a {params.gff_negative} -o results/intermediate/negative.exon.FeatureCounts.txt {input.mstrand} >> {log} 2>&1
    binaries/featureCounts -t exon --primary -g {config[gff_gene_identifier]} -T {threads} -R BAM -p --countReadPairs -O -M --largestOverlap --fracOverlap 0.1 -a {params.gff} -o results/intermediate/no.exon.FeatureCounts.txt {input.nostrand} >> {log} 2>&1
    touch {output.done}
    """

rule rename_tags_exon:
    input: pstrand = "results/intermediate/{name}.trimmed.aligned.pstrand.bam.featureCounts.bam".format(name=config["name"]),
           mstrand = "results/intermediate/{name}.trimmed.aligned.mstrand.bam.featureCounts.bam".format(name=config["name"]),
           nostrand = "results/intermediate/{name}.trimmed.aligned.nostrand.bam.featureCounts.bam".format(name=config["name"])
    output: pstrand = temp("results/intermediate/{name}.trimmed.aligned.pstrand.Exon.bam".format(name=config["name"])),
            mstrand = temp("results/intermediate/{name}.trimmed.aligned.mstrand.Exon.bam".format(name=config["name"])),
            nostrand = temp("results/intermediate/{name}.trimmed.aligned.nostrand.Exon.bam".format(name=config["name"])),
            done = "results/dones/{name}_rename_tags_exon.done".format(name=config["name"])
    log: "results/logs/rename_tags_exon.log"
    benchmark: "results/benchmarks/rename_tags_exon.benchmark.txt"
    threads: config["threads"]
    shell:"""
    binaries/rename_tags --input {input.pstrand} --output {output.pstrand} --threads {threads} >> {log} 2>&1
    binaries/rename_tags --input {input.mstrand} --output {output.mstrand} --threads {threads} >> {log} 2>&1
    binaries/rename_tags --input {input.nostrand} --output {output.nostrand} --threads {threads} >> {log} 2>&1
    touch {output.done}
    """

rule assign_genes_intron:
    input: pstrand = "results/intermediate/{name}.trimmed.aligned.pstrand.Exon.bam".format(name=config["name"]),
           mstrand = "results/intermediate/{name}.trimmed.aligned.mstrand.Exon.bam".format(name=config["name"]),
           nostrand = "results/intermediate/{name}.trimmed.aligned.nostrand.Exon.bam".format(name=config["name"])
    output: pstrand = temp("results/intermediate/{name}.trimmed.aligned.pstrand.Exon.bam.featureCounts.bam".format(name=config["name"])),
            mstrand = temp("results/intermediate/{name}.trimmed.aligned.mstrand.Exon.bam.featureCounts.bam".format(name=config["name"])),
            nostrand = temp("results/intermediate/{name}.trimmed.aligned.nostrand.Exon.bam.featureCounts.bam".format(name=config["name"])),
            done = "results/dones/{name}_assign_genes_intron.done".format(name=config["name"])
    log: "results/logs/assign_genes_intron.log"
    benchmark: "results/benchmarks/assign_genes_intron.benchmark.txt"
    params: gff = "{}.gff3".format(GFF),
            gff_positive = "{}.positive.gff3".format(GFF),
            gff_negative = "{}.negative.gff3".format(GFF)
    threads: min(config["threads"], 64)
    shell:"""
    binaries/featureCounts -t intron --primary -g {config[gff_gene_identifier]} -T {threads} -R BAM -p --countReadPairs -O -M --largestOverlap --fracOverlap 0.1 -a {params.gff_positive} -o results/intermediate/positive.intron.FeatureCounts.txt {input.pstrand} >> {log} 2>&1
    binaries/featureCounts -t intron --primary -g {config[gff_gene_identifier]} -T {threads} -R BAM -p --countReadPairs -O -M --largestOverlap --fracOverlap 0.1 -a {params.gff_negative} -o results/intermediate/negative.intron.FeatureCounts.txt {input.mstrand} >> {log} 2>&1
    binaries/featureCounts -t intron --primary -g {config[gff_gene_identifier]} -T {threads} -R BAM -p --countReadPairs -O -M --largestOverlap --fracOverlap 0.1 -a {params.gff} -o results/intermediate/no.intron.FeatureCounts.txt {input.nostrand} >> {log} 2>&1
    touch {output.done}
    """

rule rename_tags_intron:
    input: pstrand = "results/intermediate/{name}.trimmed.aligned.pstrand.Exon.bam.featureCounts.bam".format(name=config["name"]),
           mstrand = "results/intermediate/{name}.trimmed.aligned.mstrand.Exon.bam.featureCounts.bam".format(name=config["name"]),
           nostrand = "results/intermediate/{name}.trimmed.aligned.nostrand.Exon.bam.featureCounts.bam".format(name=config["name"])
    output: pstrand = temp("results/intermediate/{name}.trimmed.aligned.pstrand.GeneTagged.bam".format(name=config["name"])),
            mstrand = temp("results/intermediate/{name}.trimmed.aligned.mstrand.GeneTagged.bam".format(name=config["name"])),
            nostrand = temp("results/intermediate/{name}.trimmed.aligned.nostrand.GeneTagged.bam".format(name=config["name"])),
            done = "results/dones/{name}_rename_tags_intron.done".format(name=config["name"])
    log: "results/logs/rename_tags_intron.log"
    benchmark: "results/benchmarks/rename_tags_intron.benchmark.txt"
    threads: config["threads"]
    shell:"""
    binaries/rename_tags --input {input.pstrand} --output {output.pstrand} --intron-mode --threads {threads} >> {log} 2>&1
    binaries/rename_tags --input {input.mstrand} --output {output.mstrand} --intron-mode --threads {threads} >> {log} 2>&1
    binaries/rename_tags --input {input.nostrand} --output {output.nostrand} --intron-mode --threads {threads} >> {log} 2>&1
    touch {output.done}
    """

rule concatenate_and_sort:
    input: pstrand = "results/intermediate/{name}.trimmed.aligned.pstrand.GeneTagged.bam".format(name=config["name"]),
           mstrand = "results/intermediate/{name}.trimmed.aligned.mstrand.GeneTagged.bam".format(name=config["name"]),
           nostrand = "results/intermediate/{name}.trimmed.aligned.nostrand.GeneTagged.bam".format(name=config["name"])
    output: bam = temp("results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.bam".format(name=config["name"])),
            done = "results/dones/{name}_concatenate_and_sort.done".format(name=config["name"])
    log: "results/logs/concatenate_and_sort.log"
    benchmark: "results/benchmarks/concatenate_and_sort.benchmark.txt"
    threads: config["threads"]
    shell:"""
    mkdir -p results/.tmp_bgab/
    samtools cat {input.nostrand} {input.pstrand} {input.mstrand} | samtools sort -m 1000M -@ {threads} -T results/.tmp_bgab/sorttmp. -o {output.bam} > {log} 2>&1
    touch {output.done}
    """

rule first_index:
    input: "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.bam".format(name=config["name"])
    output: bai = temp("results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.bam.bai".format(name=config["name"])),
            done = "results/dones/{name}_first_index.done".format(name=config["name"])
    log: "results/logs/first_index.log"
    threads: config["threads"]
    shell: """
    samtools index -@ {threads} {input} > {log} 2>&1
    touch {output.done}
    """

if config["reverse"]:
    rule reconstruct:
        input: bam = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.bam".format(name=config["name"]),
               bai = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.bam.bai".format(name=config["name"]),
               sample_map = "results/metadata/{name}_sample_map.yaml".format(name=config["name"])
        output: bam = temp("results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.reconstructed.bam".format(name=config["name"])),
                merged_genes = "results/intermediate/{name}.merged_genes.csv".format(name=config["name"]),
                done = "results/dones/{name}_reconstruct.done".format(name=config["name"])
        log: "results/logs/reconstruct.log"
        benchmark: "results/benchmarks/reconstruct.benchmark.txt"
        threads: min(config["threads"], 64)
        params: comp_threads = int(min(config["threads"], 64)*0.2), proc_threads = min(config["threads"], 64)-int(min(config["threads"], 64)*0.2), gff = "{}.gff3".format(GFF)
        shell:"""
        echo Step 6/7 Reconstruct Molecules
        binaries/basic_reconstruction --input {input.bam} --output {output.bam} --gtf {params.gff} --sample-map {input.sample_map}  --threads {params.proc_threads} --bam-threads {params.comp_threads} --gene-identifier {config[gff_gene_identifier]} --merged-genes {output.merged_genes} --reverse --bulk > {log} 2>&1
        touch {output.done}
        """
else:
    rule reconstruct:
        input: bam = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.bam".format(name=config["name"]),
               bai = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.bam.bai".format(name=config["name"]),
               sample_map = "results/metadata/{name}_sample_map.yaml".format(name=config["name"])
        output: bam = temp("results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.reconstructed.bam".format(name=config["name"])),
                merged_genes = "results/intermediate/{name}.merged_genes.csv".format(name=config["name"]),
                done = "results/dones/{name}_reconstruct.done".format(name=config["name"])
        log: "results/logs/reconstruct.log"
        benchmark: "results/benchmarks/reconstruct.benchmark.txt"
        threads: min(config["threads"], 64)
        params: comp_threads = int(min(config["threads"], 64)*0.2), proc_threads = min(config["threads"], 64)-int(min(config["threads"], 64)*0.2),gff = "{}.gff3".format(GFF)
        shell:"""
        echo Step 6/7 Reconstruct Molecules
        binaries/basic_reconstruction --input {input.bam} --output {output.bam} --gtf {params.gff} --sample-map {input.sample_map} --threads {params.proc_threads} --bam-threads {params.comp_threads} --gene-identifier {config[gff_gene_identifier]} --merged-genes {output.merged_genes} --bulk > {log} 2>&1
        touch {output.done}
        """

rule sort_reconstructed:
    input: "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.reconstructed.bam".format(name=config["name"])
    output: bam = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.reconstructed.sorted.bam".format(name=config["name"]),
            done = "results/dones/{name}_sort_reconstructed.done".format(name=config["name"])
    log: "results/logs/sort_reconstructed.log"
    threads: config["threads"]
    shell: """
    mkdir -p results/.tmp_bgab/
    samtools sort -m 1000M -@ {threads} -o {output.bam} -T results/.tmp_bgab/sorttmp. {input} > {log} 2>&1
    touch {output.done}
    """

rule index_reconstructed:
    input: "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.reconstructed.sorted.bam".format(name=config["name"])
    output: bai = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.reconstructed.sorted.bam.bai".format(name=config["name"]),
            done = "results/dones/{name}_index_reconstructed.done".format(name=config["name"])
    log: "results/logs/index_reconstructed.log"
    threads: config["threads"]
    shell: """
    samtools index -@ {threads} {input} > {log} 2>&1
    touch {output.done}
    """

rule stitch_reconstruction:
    input: bam =  "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.reconstructed.sorted.bam".format(name=config["name"]),
           bai = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.reconstructed.sorted.bam.bai".format(name=config["name"]),
           merged_genes = "results/intermediate/{name}.merged_genes.csv".format(name=config["name"])
    output: bam = temp("results/intermediate/{name}.stitched.bam".format(name=config["name"])),
            done = "results/dones/{name}_stitch_reconstruction.done".format(name=config["name"])
    params: gff = "{}.gff3".format(GFF)
    log: "results/logs/stitch_reconstruction.log"
    benchmark: "results/benchmarks/stitch_reconstruction.benchmark.txt"
    threads: config["threads"]
    shell: """
    echo Step 7/7 Stitch Molecules
    binaries/stitcher_rs --input {input.bam} --output {output.bam} --gtf {params.gff} --threads {threads} --cell-tag CB --umi-tag RM --gene-identifier {config[gff_gene_identifier]} --merged-genes {input.merged_genes} > {log} 2>&1
    touch {output.done}
    """

rule sort_stitched:
    input: "results/intermediate/{name}.stitched.bam".format(name=config["name"])
    output: bam = temp("results/intermediate/{name}.stitched.sorted.bam".format(name=config["name"])),
            done = "results/dones/{name}_sort_stitched.done".format(name=config["name"])
    log: "results/logs/sort_stitched.log"
    threads: config["threads"]
    shell: """
    mkdir -p results/.tmp_bgab/
    samtools sort -m 1000M -@ {threads} -o {output.bam} -T results/.tmp_bgab/sorttmp. {input} > {log} 2>&1
    touch {output.done}
    """

rule index_stitched:
    input: "results/intermediate/{name}.stitched.sorted.bam".format(name=config["name"])
    output: bai = "results/intermediate/{name}.stitched.sorted.bam.bai".format(name=config["name"]),
            done = "results/dones/{name}_index_stitched.done".format(name=config["name"])
    log: "results/logs/index_stitched.log"
    threads: config["threads"]
    shell: """
    samtools index -@ {threads} {input} > {log} 2>&1
    touch {output.done}
    """

rule make_molecule_bam:
    input: bam = "results/intermediate/{name}.stitched.sorted.bam".format(name=config["name"]), 
           bai = "results/intermediate/{name}.stitched.sorted.bam.bai".format(name=config["name"])
    output: bam = temp("results/intermediate/{name}.stitched.molecules.bam".format(name=config["name"])),
            done = "results/dones/{name}_make_molecule_bam.done".format(name=config["name"])
    log: "results/logs/make_molecule_bam.log"
    shell:"""
    python3 workflow/scripts/filter_stitched_bam.py --input {input.bam} --molecules-out {output.bam} > {log} 2>&1
    touch {output.done}
    """

rule sort_molecule_bam:
    input: "results/intermediate/{name}.stitched.molecules.bam".format(name=config["name"])
    output: bam = "results/{name}.stitched.molecules.sorted.bam".format(name=config["name"]),
            done = "results/dones/{name}_sort_molecule_bam.done".format(name=config["name"])
    log: "results/logs/sort_molecule_bam.log"
    threads: config["threads"]
    shell:"""
    echo Creating final output file
    mkdir -p results/.tmp_bgab/
    samtools sort -m 1000M -@ {threads} -o {output.bam} -T results/.tmp_bgab/sorttmp. {input} > {log} 2>&1
    touch {output.done}
    """

rule index_molecule_bam:
    input: "results/{name}.stitched.molecules.sorted.bam".format(name=config["name"])
    output: bai = "results/{name}.stitched.molecules.sorted.bam.bai".format(name=config["name"]),
            done = "results/dones/{name}_index_molecule_bam.done".format(name=config["name"])
    log: "results/logs/index_molecule_bam.log"
    threads: config["threads"]
    shell:"""
    samtools index -@ {threads} {input} > {log} 2>&1
    touch {output.done}
    """