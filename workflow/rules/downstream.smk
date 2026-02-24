rule extract_polya:
    input: bam = "results/{name}.stitched.molecules.sorted.bam".format(name=config["name"]), fai = GENOMEREFINDEX
    output: bed_file = "results/downstream/{name}_polyA.bed".format(name=config["name"]), merged_bed_file = "results/downstream/{name}_polyA_merged.bed".format(name=config["name"]), bedgraph_file = "results/downstream/{name}_polyA.bedgraph".format(name=config["name"]), genome_file = temp("results/downstream/{name}.genome.txt".format(name=config["name"]))
    shell:
        """
        python generate_polya_bed.py --input {input.bam} --output {output.bed_file}
        bedtools merge -i {output.bed_file} -c 4 -o count > {output.merged_bed_file} 
        cut -f1,2 {input.fai} | sort -k1,1V > {output.genome_file}
        bedtools genomecov -i {output.bed_file} -g {output.genome_file} -bg > {output.bedgraph_file} 
        """