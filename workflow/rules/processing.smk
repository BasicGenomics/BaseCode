rule log_python_version:
    conda: "../envs/full.yaml"
    output: "results/logs/python_version.log"
    shell: "python3 --version > {output}"

rule make_barcode_files:
    input: samplesheet = config["samplesheet"], fastq = config["r2"]
    output: barcodes = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"]),
        cell_barcodes = "results/metadata/{name}_cell_barcodes.txt".format(name=config["name"]),
        sample_map = "results/metadata/{name}_sample_map.yaml".format(name=config["name"]),
        readtype_map = "results/metadata/{name}_readtype_map.yaml".format(name=config["name"]),
        samplesheet_out = "results/metadata/{name}_samplesheet.csv".format(name=config["name"])
    params: index_sequences = "config/index_sequences.yaml"
    log: "results/logs/make_barcode_files.log"
    conda: "../envs/full.yaml"
    shell: "echo Creating barcode files && python3 workflow/scripts/make_sample_files.py -s {input.samplesheet} --fastq {input.fastq} --index-sequences {params.index_sequences} --sample-barcodes {output.barcodes} --cell-barcodes {output.cell_barcodes} --sample-map {output.sample_map} --readtype-map {output.readtype_map} --samplesheet-out {output.samplesheet_out} > {log} 2>&1"

rule parse_fastq:
    input: r1_in = config["r1"], r2_in = config["r2"], pbcpath = "results/metadata/{name}_sample_barcodes.txt".format(name=config["name"]), cell_barcodes = "results/metadata/{name}_cell_barcodes.txt".format(name=config["name"]), sample_map = "results/metadata/{name}_sample_map.yaml".format(name=config["name"]), readtype_map = "results/metadata/{name}_readtype_map.yaml".format(name=config["name"])
    output: r1_out = temp("results/intermediate/{name}.read1.fastq.gz".format(name=config["name"])),
            r2_out = temp("results/intermediate/{name}.read2.fastq.gz".format(name=config["name"]))
    log: "results/logs/parse_fastq.log"
    benchmark: "results/benchmarks/parse_fastq.benchmark.txt"
    params: comp_threads = int(config["threads"]*0.2),
            proc_threads = config["threads"]-int(config["threads"]*0.2)
    shell: "echo Parsing reads && {config[resource_dir]}/binaries/parse_fastq --read1 {input.r1_in} --read2 {input.r2_in} --r1-out {output.r1_out} --r2-out {output.r2_out} --cbcpath {input.cell_barcodes} --pbcpath {input.pbcpath} --readtype-structure {input.readtype_map} --index-layout {config[index_layout]} --sample-structure {input.sample_map} --processing-threads {params.proc_threads} --compression-threads {params.comp_threads}  --umilen {config[umilen]} --dtlen {config[dtlen]} --dt-cutoff {config[dt_cutoff]} --ts-sequence {config[ts_sequence]} --ts-pad {config[ts_pad]} --ts-cutoff {config[ts_cutoff]} > {log} 2>&1"

rule trim_fastq:
    input: r1 = "results/intermediate/{name}.read1.fastq.gz".format(name=config["name"]),
           r2 = "results/intermediate/{name}.read2.fastq.gz".format(name=config["name"])
    output: r1 = temp("results/intermediate/{name}.trimmed.read1.fastq.gz".format(name=config["name"])),
           r2 = temp("results/intermediate/{name}.trimmed.read2.fastq.gz".format(name=config["name"])),
           r1_short = temp("results/intermediate/{name}.tooshort.read1.fastq.gz".format(name=config["name"])),
           r2_short = temp("results/intermediate/{name}.tooshort.read2.fastq.gz".format(name=config["name"]))
    log: stdout = "results/logs/trim_fastq.log",
         summary = "results/summaries/{name}.cutadapt.json".format(name=config["name"])
    benchmark: "results/benchmarks/trim_fastq.benchmark.txt"
    conda: "../envs/full.yaml"
    threads: config["threads"]
    shell: "echo Trimming reads && cutadapt -j {threads} --json {log.summary} {config[params][cutadapt]} --too-short-output {output.r1_short} --too-short-paired-output {output.r2_short} -o {output.r1} -p {output.r2} {input.r1} {input.r2} > {log.stdout} 2>&1"


rule map_reads:
    input: r1 = "results/intermediate/{name}.trimmed.read1.fastq.gz".format(name=config["name"]),
        r2 = "results/intermediate/{name}.trimmed.read2.fastq.gz".format(name=config["name"])
    output: temp("results/intermediate/{name}.trimmed.aligned.bam".format(name=config["name"]))
    log: stdout = "results/logs/map_reads.log",
        summary = "results/summaries/{name}.hisat2.summary.txt".format(name=config["name"])
    benchmark: "results/benchmarks/map_reads.benchmark.txt"
    params: splicesites = SPLICESITES, genomeref = GENOMEREF, cores_hisat = cores_hisat, cores_samtools = cores_samtools
    conda: "../envs/full.yaml"
    shell: "echo Mapping reads && {config[resource_dir]}/binaries/hisat-3n --new-summary --summary-file {log.summary} {config[params][hisat3n]} -p {params.cores_hisat} --known-splicesite-infile {params.splicesites} -x {params.genomeref} -1 {input.r1} -2 {input.r2} | samtools view -F 256 -b -@ {params.cores_samtools} -o {output} > {log.stdout} 2>&1"

rule split_bam_by_strand:
    input: "results/intermediate/{name}.trimmed.aligned.bam".format(name=config["name"])
    output: pstrand = temp("results/intermediate/{name}.trimmed.aligned.pstrand.bam".format(name=config["name"])),
            mstrand = temp("results/intermediate/{name}.trimmed.aligned.mstrand.bam".format(name=config["name"])),
            nostrand = temp("results/intermediate/{name}.trimmed.aligned.nostrand.bam".format(name=config["name"]))
    log: "results/logs/split_bam_by_strand.log"
    benchmark: "results/benchmarks/split_bam_by_strand.benchmark.txt"
    conda: "../envs/full.yaml"
    shell: "echo Assigning reads to genes && python3 workflow/scripts/split_bam_by_strand.py {input} > {log} 2>&1"

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
    threads: config["threads"]
    shell:"""
    featureCounts -t exon --primary -g gene_name -T {threads} -R BAM -p --countReadPairs --largestOverlap --fracOverlap 0.1 -a {params.gtffile_positive} -o results/intermediate/pos.tmp {input.pstrand} >> {log} 2>&1
    featureCounts -t exon --primary -g gene_name -T {threads} -R BAM -p --countReadPairs --largestOverlap --fracOverlap 0.1 -a {params.gtffile_negative} -o results/intermediate/neg.tmp {input.mstrand} >> {log} 2>&1
    featureCounts -t exon --primary -g gene_name -T {threads} -R BAM -p --countReadPairs --largestOverlap --fracOverlap 0.1 -a {params.gtffile} -o results/intermediate/no.tmp {input.nostrand} >> {log} 2>&1
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
    shell: "python3 workflow/scripts/rename_tags_exon.py {input.pstrand} {input.mstrand} {input.nostrand} {output.pstrand} {output.mstrand} {output.nostrand} > {log} 2>&1"

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
    threads: config["threads"]
    shell:"""
    featureCounts -t intron --primary -g gene_name -T {threads} -R BAM -p --countReadPairs --largestOverlap --fracOverlap 0.1 -a {params.gtffile_positive} -o results/intermediate/pos.tmp {input.pstrand} >> {log} 2>&1
    featureCounts -t intron --primary -g gene_name -T {threads} -R BAM -p --countReadPairs --largestOverlap --fracOverlap 0.1 -a {params.gtffile_negative} -o results/intermediate/neg.tmp {input.mstrand} >> {log} 2>&1
    featureCounts -t intron --primary -g gene_name -T {threads} -R BAM -p --countReadPairs --largestOverlap --fracOverlap 0.1 -a {params.gtffile} -o results/intermediate/no.tmp {input.nostrand} >> {log} 2>&1
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
    shell: "python3 workflow/scripts/rename_tags_intron.py {input.pstrand} {input.mstrand} {input.nostrand} {output.pstrand} {output.mstrand} {output.nostrand} > {log} 2>&1"

rule concatenate_and_sort:
    input: pstrand = "results/intermediate/{name}.trimmed.aligned.pstrand.GeneTagged.bam".format(name=config["name"]),
            mstrand = "results/intermediate/{name}.trimmed.aligned.mstrand.GeneTagged.bam".format(name=config["name"]),
            nostrand = "results/intermediate/{name}.trimmed.aligned.nostrand.GeneTagged.bam".format(name=config["name"])
    log: "results/logs/concatenate_and_sort.log"
    benchmark: "results/benchmarks/concatenate_and_sort.benchmark.txt"
    output: "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.bam".format(name=config["name"])
    conda: "../envs/full.yaml"
    threads: config["threads"]
    shell: "samtools cat {input.nostrand} {input.pstrand} {input.mstrand} | samtools sort -m 1000M -@ {threads} -T results/.tmp_bgab/sorttmp. -o {output} >> {log} 2>&1"

rule first_index:
    input: "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.bam".format(name=config["name"])
    output: "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.bam.bai".format(name=config["name"])
    params:
        extra="",
    threads: config["threads"]
    log: "results/logs/first_index.log"
    wrapper:
        "v3.3.3/bio/samtools/index" 

rule reconstruct:
    input: bam = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.bam".format(name=config["name"]),
        bai = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.bam.bai".format(name=config["name"]),
        sample_map = "results/metadata/{name}_sample_map.yaml".format(name=config["name"])
    output: temp("results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.reconstructed.bam".format(name=config["name"]))
    params: gtffile = GTFFILE
    log: "results/logs/reconstruct.log"
    threads: config["threads"]
    benchmark: "results/benchmarks/reconstruction.benchmark.txt"
    shell: "echo Reconstructing reads && {config[resource_dir]}/binaries/basic_reconstruction --input {input.bam} --output {output} --gtf {params.gtffile}.gff3 --sample-map {input.sample_map} --threads {threads} >> {log} 2>&1"

rule sort_reconstructed:
    input: "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.reconstructed.bam".format(name=config["name"])
    output: "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.reconstructed.sorted.bam".format(name=config["name"])
    threads: config["threads"]
    log: "results/logs/sort_reconstructed.log"
    params:
        extra="-m 1000M",
    conda: "../envs/full.yaml"
    shell: "samtools sort -@ {threads} -o {output} {params.extra} -T . {input} &> {log}"

rule index_reconstructed:
    input: "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.reconstructed.sorted.bam".format(name=config["name"])
    output: "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.reconstructed.sorted.bam.bai".format(name=config["name"])
    log: "results/logs/index_reconstructed.log"
    params:
        extra="",
    threads: config["threads"]
    log: "results/logs/index_reconstructed.log"
    wrapper:
        "v3.3.3/bio/samtools/index"

rule stitch_reconstruction:
    input: bam =  "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.reconstructed.sorted.bam".format(name=config["name"]),
        bai = "results/intermediate/{name}.reads.aligned_trimmed_genetagged_sorted.reconstructed.sorted.bam.bai".format(name=config["name"])
    output: temp("results/intermediate/{name}.stitched.bam".format(name=config["name"]))
    threads: config["threads"]
    params: gtffile = GTFFILE
    log: "results/logs/stitch_reconstruction.log"
    conda: "../envs/full.yaml"
    shell: "echo Combining reads into synthetic long reads && python3 workflow/scripts/stitcher.py --input {input.bam} --output {output} --gtf {params.gtffile}.gff3 --threads {threads} --cell-tag CB --UMI-tag RM --gene-identifier gene_id >> {log} 2>&1"

rule sorted_stitched:
    input: "results/intermediate/{name}.stitched.bam".format(name=config["name"])
    output: "results/intermediate/{name}.stitched.sorted.bam".format(name=config["name"])
    threads: config["threads"]
    log: "results/logs/sort_stitched.log"
    params:
        extra="-m 1000M",
    conda: "../envs/full.yaml"
    shell: "samtools sort -@ {threads} -o {output} {params.extra} -T . {input} &> {log}"

rule index_stitched:
    input: "results/intermediate/{name}.stitched.sorted.bam".format(name=config["name"])
    output: "results/intermediate/{name}.stitched.sorted.bam.bai".format(name=config["name"])
    params:
        extra="",
    threads: config["threads"]
    log: "results/logs/index_stitched.log"
    wrapper:
        "v3.3.3/bio/samtools/index"

rule make_molecule_bams:
    input: bam = "results/intermediate/{name}.stitched.sorted.bam".format(name=config["name"]), bai = "results/intermediate/{name}.stitched.sorted.bam.bai".format(name=config["name"])
    output: molecules_out = temp("results/intermediate/{name}.stitched.molecules.bam".format(name=config["name"]))
    conda: "../envs/full.yaml"
    shell:"python3 workflow/scripts/filter_stitched_bam.py --input {input.bam} --molecules-out {output.molecules_out}"

rule sort_molecule_bams:
    input: molecules_out = "results/intermediate/{name}.stitched.molecules.bam".format(name=config["name"])
    output: molecules_out = "results/{name}.stitched.molecules.sorted.bam".format(name=config["name"])
    log: "results/logs/sort_molecules.log"
    conda: "../envs/full.yaml"
    shell:"""
    echo Creating final output file
    samtools sort -m 1000M -@ {config[threads]} -T results/.tmp_bgab/sorttmp. -o {output.molecules_out} {input.molecules_out} > {log} 2>&1
    samtools index -@ {config[threads]} {output.molecules_out}
    """